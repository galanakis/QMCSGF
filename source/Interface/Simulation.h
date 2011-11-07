#include <fstream>
#include <ctime>

template<class T> T Min(const T &a,const T &b) {return (a<b)?a:b;} 
template<class T> T Max(const T &a,const T &b) {return (a>b)?a:b;} 


template<class T>
class SpeedoMeter {
	clock_t Time;
	T value;
	T Speed;
	unsigned long count;
public:
	SpeedoMeter() : Time(clock()), value(0), Speed(0), count(0) {}
	inline void measure(const T &_val) {
		clock_t OldTime=Time;
		T oldvalue=value;
		Time=clock();
		value=_val;
		T InstantSpeed=(value-oldvalue)*CLOCKS_PER_SEC/(Time-OldTime);
		Speed+=InstantSpeed;
		++count;
		return InstantSpeed;
	}
	inline T speed() const {return Speed/count;}
	inline T remaining_time(const T & val) const {return (val-value)/speed();}
}; 

class FileNameProgressBar {
	static const int StringLength=200;
	char Status[StringLength],*Ptr;
	double Progress;
	unsigned long NumUpdates;
	clock_t Time;
public:
	FileNameProgressBar(const char *Name) : Progress(0), NumUpdates(0) {
		Time=clock();
		strcpy(Status,MathExpression::GetSimulName());
		std::string append=std::string(": ")+std::string(Name)+std::string(" ");
		strcpy(Status+strlen(Status),append.c_str());
		Ptr=Status+strlen(Status);   
		sprintf(Ptr," - 000 - 000:00:00");     
		fstream File;
		File.open(Status,ios::out);
		File.close();    
	}
	~FileNameProgressBar() {
		remove(Status);
	}

	inline void reset_cout() const {
		for(unsigned int i=0;i<strlen(Status);++i)
			cout<<'\b';  
	}
	inline void Update(double prog,unsigned long NewNumUpdates) {
		if(prog<0) prog=0.0;
		if(prog>1) prog=1.0;
		if(static_cast<int>(100*(prog-Progress))!=0) {
			clock_t Now=clock();
			double DeltaSeconds=double(Now-Time)/CLOCKS_PER_SEC;
			unsigned int Speed=static_cast<unsigned int>((NewNumUpdates-NumUpdates)/DeltaSeconds);
			unsigned int SecondsLeft=(1-prog)*DeltaSeconds/(prog-Progress);
			Time=Now;
			Progress=prog;
			NumUpdates=NewNumUpdates;

			unsigned int HoursLeft=SecondsLeft/3600;
			SecondsLeft%=3600;
			unsigned int MinutesLeft=SecondsLeft/60;
			SecondsLeft%=60;

			char OldStatus[StringLength];
			strcpy(OldStatus,Status);
			int percent=static_cast<int>(100*Progress);

			sprintf(Ptr," - %.3d %% - %.3d:%.2d:%.2d - %6d updates/second",percent,HoursLeft,MinutesLeft,SecondsLeft,Speed);
      
			reset_cout();
			std::cout<<Status<<std::flush;


			rename(OldStatus,Status);


		}
	}



};

class Simulation
{
	unsigned long long NumWarmUpdates,NumMeasUpdates;
	unsigned long long NumDirectedWarmUpdates,NumDirectedMeasUpdates;
	unsigned long long WarmTime,MeasTime;
	unsigned long long WarmIterations, MeasIterations;

	double ActualWarmTime,ActualMeasTime;

	SGF::OperatorStringType &OpString;
	SGF::Measurable &MeasuredOp;

public:
	Simulation(SGF::OperatorStringType &OStr, SGF::Measurable &MOp) : OpString(OStr), MeasuredOp(MOp) {
		WarmTime=MathExpression::GetValue("#WarmTime").Re();
		WarmIterations=(MathExpression::Find("WarmIterations")!=NULL) ? static_cast<unsigned long long>(MathExpression::GetValue("WarmIterations").Re()) : std::numeric_limits<unsigned long long>::max();
		MeasTime=MathExpression::GetValue("#MeasTime").Re();
		MeasIterations=(MathExpression::Find("MeasIterations")!=NULL) ? static_cast<unsigned long long>(MathExpression::GetValue("MeasIterations").Re()) : std::numeric_limits<unsigned long long>::max();			
	} 

	void Thermalize()
	{
		FileNameProgressBar pbar("Thermalizing");

		const unsigned int AlphaUpdatePeriod=100000;

		clock_t StartTime=clock();
		clock_t EndTime=StartTime+WarmTime*CLOCKS_PER_SEC;

					// Run a few iterations first to determine the cpu speed.
		NumWarmUpdates=0;
		NumDirectedWarmUpdates=0;

		do {

			NumWarmUpdates+=OpString.directed_update(); // Perform an update. 

			++NumDirectedWarmUpdates;
			if (NumDirectedWarmUpdates%AlphaUpdatePeriod==(AlphaUpdatePeriod-1)) OpString.AlphaUpdate();     

			double Progress=Max(static_cast<double>(NumWarmUpdates)/WarmIterations,static_cast<double>(clock()-StartTime)/(EndTime-StartTime));
			pbar.Update(Progress,NumWarmUpdates);

		} while ( NumWarmUpdates<WarmIterations && clock()<EndTime );

		while (OpString.NBrokenLines()!=0)  {              // Perform extra updates until we end up
			NumWarmUpdates+=OpString.directed_update();      // in a diagonal configuration.
			++NumDirectedWarmUpdates;
		}

		ActualWarmTime=double(clock()-StartTime)/CLOCKS_PER_SEC;
		pbar.reset_cout();
		std::cout<<"Done Thermalizing after "<<NumWarmUpdates<<" updates ("<<NumDirectedWarmUpdates<<" directed) in "<<ActualWarmTime<<" seconds at "<<NumWarmUpdates/ActualWarmTime<<" updates per second."<<std::endl;

	}

	void Measure()
	{
		FileNameProgressBar pbar("Measuring   ");

		clock_t StartTime=clock();
		clock_t EndTime=StartTime+MeasTime*CLOCKS_PER_SEC;

		unsigned long long NumBins=MathExpression::GetValue("#Bins").Re();

		NumMeasUpdates=0;
		NumDirectedMeasUpdates=0;

		for (unsigned int i=0;i<NumBins;++i) {

			unsigned long long counter=0;
			do {
				counter+=OpString.directed_update();   // Perform an update.
				++NumDirectedMeasUpdates;
				MeasuredOp.measure(OpString);          // Perform measurements.
			} while ( counter<MeasIterations/NumBins && clock()<EndTime );

			NumMeasUpdates+=counter;

			while (OpString.NBrokenLines()!=0) {              // Perform extra updates until we end up
				NumMeasUpdates+=OpString.directed_update();     // in a diagonal configuration.
				++NumDirectedMeasUpdates;
				MeasuredOp.measure(OpString);                   // Perform measurements.
			}
			MeasuredOp.flush();                               // Bin the data. 
			double Progress=Max(static_cast<double>(i)/NumBins,static_cast<double>(clock()-StartTime)/(EndTime-StartTime));
			pbar.Update(Progress,NumMeasUpdates);

		}

		ActualMeasTime=double(clock()-StartTime)/CLOCKS_PER_SEC; 
		pbar.reset_cout();
		std::cout<<"Done Measuring after "<<NumMeasUpdates<<" updates ("<<NumDirectedMeasUpdates<<" directed) in "<<ActualMeasTime<<" seconds at "<<NumMeasUpdates/ActualMeasTime<<" updates/second."<<std::endl;
	}

	void Results()
	{
		cout << "*******************************************************************************************\n";
		cout << "* This is a quantum Monte Carlo simulation performed by the \"Corvette SGF engine\".        *\n";
		cout << "*                                                                                         *\n";
		cout << "* For informations on the SGF algorithm:                                                  *\n";
		cout << "*   Physical Review E 77, 056705 (2008)                                                   *\n";
		cout << "*   Physical Review E 78, 056707 (2008)                                                   *\n";
		cout << "*                                                                                         *\n";
		cout << "* Dr Valy G. Rousseau and Dr Dimitris Galanakis - Version " << Version;
		for (unsigned int i=0;i<32-strlen(Version);i++) cout << " ";
		cout << "*\n";
		cout << "*******************************************************************************************\n\n";
		cout << "***************************\n";
		cout << "* User's input parameters *\n";
		cout << "***************************\n\n  ";
		Parser::TokenHandle Input=Parser::First();

		while (Input)
		{
			if (Input->Type()==Parser::Number)
				cout << *(double *) Input->Value();

			else
			{
				char *C=(char *) Input->Value();

				if (Input->Type()==Parser::String)
					cout << "\"" << C << "\"";
				else
					cout << C;

				if (MathExpression::IsKeyword(C))
					cout << " ";

				else if (*C==';')
				{
					cout << "\n";

					if (Input->NextToken())
						cout << "  ";
				}
			}

			Input=Input->NextToken();
		}

		cout << endl;
		cout << "*************************\n";
		cout << "* Results of simulation *\n";
		cout << "*************************\n\n";
		cout << "  ******************************\n";
		cout << "  * Operator string statistics *\n";
		cout << "  ******************************\n\n";
		cout << "    == Thermalization ==\n";
		cout << "    Number of creations/anihilations: \t" << NumWarmUpdates << "\t("<<ActualWarmTime*1000000000/NumWarmUpdates << " seconds per billion updates)\n";
		cout << "    Number of directed updates:       \t" << NumDirectedWarmUpdates << "\t("<<ActualWarmTime*1000000000/NumDirectedWarmUpdates << " seconds per billion updates)\n";
		cout << "    Directed update length:           \t"<< double(NumWarmUpdates)/NumDirectedWarmUpdates<<std::endl;
		cout << "    == Measurements   ==\n";
		cout << "    Number of creations/anihilations: \t" << NumMeasUpdates << "\t("<<ActualMeasTime*1000000000/NumMeasUpdates << " seconds per billion updates)\n";  
		cout << "    Number of directed updates:       \t" << NumDirectedMeasUpdates << "\t("<<ActualMeasTime*1000000000/NumDirectedMeasUpdates << " seconds per billion updates)\n";  
		cout << "    Directed update length:           \t"<< double(NumMeasUpdates)/NumDirectedMeasUpdates<<std::endl;
		cout << "    Number of measurements: " << MeasuredOp.count() << "\n\n";

		std::cout<<::std::endl;
		std::cout<<"  *******************************\n";
		std::cout<<"  * Broken worldlines histogram *\n";
		std::cout<<"  *******************************\n\n";
		std::cout<<"    N lines\tCount\tProbability\n\n";

		SGF::Measurable::BrokenHistogramType BrokenHistogram=MeasuredOp.BrokenHistogram();
		double Normalization=MeasuredOp.BrokenNormalization();
		for(SGF::Measurable::BrokenHistogramType::const_iterator it=BrokenHistogram.begin();it!=BrokenHistogram.end();++it)
			std::cout<<"    "<<it->first<<"\t\t"<<it->second<<"\t"<<it->second/Normalization<<std::endl;
		cout << endl;

		cout << "  ***********************************************************************************\n";
		cout << "  * Energies (obtained from operator string length and Green operator state energy) *\n";
		cout << "  ***********************************************************************************\n\n";
		cout << "    Total energy: " << MeasuredOp.TotalEnergy() << "\n";
		cout << "    Diagonal energy: " << MeasuredOp.PotentialEnergy() << "\n";           
		cout << "    Non-diagonal energy: " << MeasuredOp.KineticEnergy() << "\n\n";
		std::vector<string> Names=MathExpression::GetMeasurableList();

		int UserMeasurables=0;

		for (SGF::Measurable::size_type i=0;i<MeasuredOp.size();i++)
			if (Names[i][0]!='#')
			UserMeasurables++;

		if (UserMeasurables)
		{
			cout << "  ******************************\n";
			cout << "  * User's defined measurables *\n";
			cout << "  ******************************\n\n";

			for (SGF::Measurable::size_type i=0;i<MeasuredOp.size();i++)
				if (Names[i][0]!='#')
				cout << "    " << Names[i] << ": " << MeasuredOp[i] << endl;
		}
	}
};


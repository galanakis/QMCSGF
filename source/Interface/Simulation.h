#include <fstream>
#include <ctime>

template<class T> T Min(const T &a,const T &b) {return (a<b)?a:b;} 
template<class T> T Max(const T &a,const T &b) {return (a>b)?a:b;} 

class FileNameProgressBar {
	static const int StringLength=300;
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
		sprintf(Ptr," - 000 - 000h 00m 00s");     
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

			reset_cout();

			sprintf(Ptr," - %.3d %% - %.3dh %.2dm %.2ds - %6d updates per second",percent,HoursLeft,MinutesLeft,SecondsLeft,Speed);

#ifdef CMDLINEPROGRESS      
			std::cout<<Status<<std::flush;
#endif
      
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

public:
	Simulation() {
		WarmTime=MathExpression::GetValue("#WarmTime").Re();
		WarmIterations=(MathExpression::Find("WarmIterations")!=NULL) ? static_cast<unsigned long long>(MathExpression::GetValue("WarmIterations").Re()) : std::numeric_limits<unsigned long long>::max();
		MeasTime=MathExpression::GetValue("#MeasTime").Re();
		MeasIterations=(MathExpression::Find("MeasIterations")!=NULL) ? static_cast<unsigned long long>(MathExpression::GetValue("MeasIterations").Re()) : std::numeric_limits<unsigned long long>::max();			
	} 

	void Thermalize(SGF::OperatorStringType &OpString)
	{
		FileNameProgressBar pbar("Thermalizing");

		const unsigned int AlphaUpdatePeriod=100000;

		clock_t StartTime=clock();
		clock_t EndTime=StartTime+WarmTime*CLOCKS_PER_SEC;
		clock_t Now=StartTime;

		NumWarmUpdates=0;
		NumDirectedWarmUpdates=0;

		do {

			do {
				NumWarmUpdates+=OpString.directed_update();
				++NumDirectedWarmUpdates;
				if (NumDirectedWarmUpdates%AlphaUpdatePeriod==(AlphaUpdatePeriod-1)) OpString.AlphaUpdate();     
			} while(OpString.NBrokenLines()!=0);


			Now=clock();
			double Progress=Max(static_cast<double>(NumWarmUpdates)/WarmIterations,static_cast<double>(Now-StartTime)/(EndTime-StartTime));
			pbar.Update(Progress,NumWarmUpdates);

		} while ( NumWarmUpdates<WarmIterations && Now<EndTime );


		ActualWarmTime=double(clock()-StartTime)/CLOCKS_PER_SEC;
#ifdef CMDLINEPROGRESS      
		pbar.reset_cout();
		std::cout<<"Done Thermalizing after "<<NumWarmUpdates<<" updates ("<<NumDirectedWarmUpdates<<" directed) in "<<ActualWarmTime<<" seconds at "<<NumWarmUpdates/ActualWarmTime<<" updates per second."<<std::endl;
#endif
	}

	void Measure(SGF::OperatorStringType &OpString,SGF::Measurable &MeasuredOp)
	{
		FileNameProgressBar pbar("Measuring   ");

		clock_t StartTime=clock();
		clock_t EndTime=StartTime+MeasTime*CLOCKS_PER_SEC;

		unsigned long long NumBins=MathExpression::GetValue("#Bins").Re();

		NumMeasUpdates=0;
		NumDirectedMeasUpdates=0;

		for (unsigned int i=0;i<NumBins;++i) {

			unsigned long long counter=0;
			clock_t BinStart=clock();
			do {

				do {
					counter+=OpString.directed_update();   // Perform an update.
					++NumDirectedMeasUpdates;
					MeasuredOp.measure();          // Perform measurements.					
				} while(OpString.NBrokenLines()!=0);

			} while ( counter*NumBins<MeasIterations && (clock()-BinStart)*NumBins<MeasTime*CLOCKS_PER_SEC );

			NumMeasUpdates+=counter;

			MeasuredOp.flush();                               // Bin the data. 
			double Progress=Max(static_cast<double>(i)/NumBins,static_cast<double>(clock()-StartTime)/(EndTime-StartTime));
			pbar.Update(Progress,NumMeasUpdates);

		}

		ActualMeasTime=double(clock()-StartTime)/CLOCKS_PER_SEC; 
		
#ifdef CMDLINEPROGRESS      
		pbar.reset_cout();
		std::cout<<"Done Measuring after "<<NumMeasUpdates<<" updates ("<<NumDirectedMeasUpdates<<" directed) in "<<ActualMeasTime<<" seconds at "<<NumMeasUpdates/ActualMeasTime<<" updates/second."<<std::endl;
#endif
	}

	void Results(SGF::Measurable &MeasuredOp)
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
		cout << "    Number of creations/annihilations: \t" << NumWarmUpdates << "\t("<<ActualWarmTime*1000000000/NumWarmUpdates << " seconds per billion updates)\n";
		cout << "    Number of directed updates:       \t" << NumDirectedWarmUpdates << "\t("<<ActualWarmTime*1000000000/NumDirectedWarmUpdates << " seconds per billion updates)\n";
		cout << "    Directed update length:           \t"<< double(NumWarmUpdates)/NumDirectedWarmUpdates<<std::endl;
		cout << "    == Measurements   ==\n";
		cout << "    Number of creations/annihilations: \t" << NumMeasUpdates << "\t("<<ActualMeasTime*1000000000/NumMeasUpdates << " seconds per billion updates)\n";  
		cout << "    Number of directed updates:       \t" << NumDirectedMeasUpdates << "\t("<<ActualMeasTime*1000000000/NumDirectedMeasUpdates << " seconds per billion updates)\n";  
		cout << "    Directed update length:           \t"<< double(NumMeasUpdates)/NumDirectedMeasUpdates<<std::endl;

		MeasuredOp.print();
	}
};


#include <fstream>
#include <ctime>
#include <OperatorString.hh>
#include <Measurable.hh>

template<class T> T Min(const T &a,const T &b) {return (a<b)?a:b;} 
template<class T> T Max(const T &a,const T &b) {return (a>b)?a:b;} 

class FileNameProgressBar {
	static const int StringLength=300;
	char Status[StringLength],*Ptr;
	double Progress;
	unsigned long NumUpdates;
	clock_t Time;
	std::ostream &fout;
public:
	FileNameProgressBar(const std::string &Name,std::ostream &o) : Progress(0), NumUpdates(0), fout(o) {
		Time=clock();
		strcpy(Status,Name.c_str());
		Ptr=Status+strlen(Status);   
		sprintf(Ptr," - 000 - 000h 00m 00s");
#ifndef CMDLINEPROGRESS     
		std::fstream File;
		File.open(Status,std::ios::out);
		File.close(); 
#endif  
	}
	~FileNameProgressBar() {
#ifndef CMDLINEPROGRESS
		remove(Status);
#endif
	}

	inline void reset_cout() {
		for(unsigned int i=0;i<strlen(Status);++i)
			fout<<'\b';  
	}
	inline void Update(double prog,unsigned long NewNumUpdates) {
		if(prog<0) prog=0.0;
		if(prog>1) prog=1.0;
		if(static_cast<int>(100*(prog-Progress))!=0) {
			clock_t Now=clock();
			double DeltaSeconds=double(Now-Time)/CLOCKS_PER_SEC;
			int Speed=static_cast<unsigned int>((NewNumUpdates-NumUpdates)/DeltaSeconds);
			int SecondsLeft=(1-prog)*DeltaSeconds/(prog-Progress);
			Time=Now;
			Progress=prog;
			NumUpdates=NewNumUpdates;

			int HoursLeft=SecondsLeft/3600;
			SecondsLeft%=3600;
			int MinutesLeft=SecondsLeft/60;
			SecondsLeft%=60;

			char OldStatus[StringLength];
			strcpy(OldStatus,Status);
			int percent=static_cast<int>(100*Progress);


			sprintf(Ptr," - %.3d %% - %.3dh %.2dm %.2ds - %6d updates per second",percent,HoursLeft,MinutesLeft,SecondsLeft,Speed);

#ifdef CMDLINEPROGRESS      
			reset_cout();
			fout<<Status<<std::flush;
#else
			rename(OldStatus,Status);
#endif


		}
	}



};

#define MAXNUMBROKENLINES 100

class Simulation
{
	unsigned long NumWarmUpdates,NumMeasUpdates;
	unsigned long NumDirectedWarmUpdates,NumDirectedMeasUpdates;

	double ActualWarmTime,ActualMeasTime;
  
	typedef  std::vector<unsigned long> BrokenHistogramType;
	BrokenHistogramType BrokenHistorgram;

	std::string SimulName;

	std::ostream &o;

public:
	Simulation(const std::string &Name,std::ostream &_o) : SimulName(Name), o(_o) {

#ifdef USEMPI
   // ***********************************
   // * Initialization of MPI functions *
   // ***********************************

   MPI_Init(&NumArg,&Arg);
   MPI_Comm_size(MPI_COMM_WORLD,&NumProcessors);
   MPI_Comm_rank(MPI_COMM_WORLD,&Rank);
   MPI_Get_processor_name(ProcessorName,&NameLength);

   // Silence the other nodes. Only the root node can print
   o.rdbuf( Rank==Master ? std::cout.rdbuf() : 0 );

#endif

		NumWarmUpdates=0;
		NumDirectedWarmUpdates=0;
		NumMeasUpdates=0;
		NumDirectedMeasUpdates=0;
		ActualWarmTime=0;
		ActualMeasTime=0;

		BrokenHistorgram.resize(MAXNUMBROKENLINES);
	} 

	~Simulation() {

	#ifdef USEMPI
   	MPI_Finalize();
	#endif


	}

	void Thermalize(SGF::OperatorStringType &OpString,unsigned long WarmIterations,unsigned long WarmTime)
	{


		FileNameProgressBar pbar(SimulName+std::string(": Thermalizing "),o);

		clock_t StartTime=clock();
		clock_t EndTime=StartTime+WarmTime*CLOCKS_PER_SEC;
		clock_t Now=StartTime;

		NumWarmUpdates=0;
		NumDirectedWarmUpdates=0;
		NumMeasUpdates=0;
		NumDirectedMeasUpdates=0;

		do {

			do {
				NumWarmUpdates+=OpString.directed_update();
				++NumDirectedWarmUpdates;
			} while(OpString.NBrokenLines()!=0);


			Now=clock();
			double Progress=Max(static_cast<double>(NumWarmUpdates)/WarmIterations,static_cast<double>(Now-StartTime)/(EndTime-StartTime));
			pbar.Update(Progress,NumWarmUpdates);

		} while ( NumWarmUpdates<WarmIterations && Now<EndTime );


		ActualWarmTime=double(clock()-StartTime)/CLOCKS_PER_SEC;
#ifdef CMDLINEPROGRESS      
		pbar.reset_cout();
		o<<"Done Thermalizing after "<<NumWarmUpdates<<" updates ("<<NumDirectedWarmUpdates<<" directed) in "<<ActualWarmTime<<" seconds at "<<NumWarmUpdates/ActualWarmTime<<" updates per second."<<std::endl;
#endif
	}

	void Measure(SGF::OperatorStringType &OpString,SGF::Measurable &MeasuredOp,unsigned long NumBins,unsigned long MeasIterations,unsigned long MeasTime)
	{
		FileNameProgressBar pbar(SimulName+std::string(": Measuring   "),o);

		clock_t StartTime=clock();
		clock_t EndTime=StartTime+MeasTime*CLOCKS_PER_SEC;

		NumMeasUpdates=0;
		NumDirectedMeasUpdates=0;

		for (unsigned int i=0;i<NumBins;++i) {

			unsigned long counter=0;
			clock_t BinStart=clock();
			do {

				do {
					counter+=OpString.directed_update();   // Perform an update.
					++NumDirectedMeasUpdates;
					BrokenHistorgram[OpString.NBrokenLines()]+=1;
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
		o<<"Done Measuring after "<<NumMeasUpdates<<" updates ("<<NumDirectedMeasUpdates<<" directed) in "<<ActualMeasTime<<" seconds at "<<NumMeasUpdates/ActualMeasTime<<" updates/second."<<std::endl;
#endif
	}

	void Results(SGF::Measurable &MeasuredOp)
	{
		o << "*************************\n";
		o << "* Results of simulation *\n";
		o << "*************************\n\n";
		o << "  ******************************\n";
		o << "  * Operator string statistics *\n";
		o << "  ******************************\n\n";
 

#ifdef USEMPI
   
		unsigned long send_NumWarmUpdates(NumWarmUpdates),send_NumMeasUpdates(NumMeasUpdates);
		unsigned long send_NumDirectedWarmUpdates(NumDirectedWarmUpdates),send_NumDirectedMeasUpdates(NumDirectedMeasUpdates);
		double send_ActualWarmTime(ActualWarmTime),send_ActualMeasTime(ActualMeasTime);
		std::vector<unsigned long> send_BrokenHistogram=BrokenHistorgram;

		MPI_Reduce(&send_NumWarmUpdates,&NumWarmUpdates,1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
		MPI_Reduce(&send_NumDirectedWarmUpdates,&NumDirectedWarmUpdates,1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
		MPI_Reduce(&send_NumMeasUpdates,&NumMeasUpdates,1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
		MPI_Reduce(&send_NumDirectedMeasUpdates,&NumDirectedMeasUpdates,1,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);
		MPI_Reduce(&send_ActualWarmTime,&ActualWarmTime,1,MPI_DOUBLE,MPI_SUM,Master,MPI_COMM_WORLD);
		MPI_Reduce(&send_ActualMeasTime,&ActualMeasTime,1,MPI_DOUBLE,MPI_SUM,Master,MPI_COMM_WORLD);
		MPI_Reduce(&send_BrokenHistogram[0],&BrokenHistorgram[0],MAXNUMBROKENLINES,MPI_UNSIGNED_LONG,MPI_SUM,Master,MPI_COMM_WORLD);

		o << "    Number of Processors used "<<NumProcessors<<"\n\n";

#endif   	


		o << "    == Thermalization ==\n";
		o << "    Number of creations/annihilations: \t" << NumWarmUpdates << "\t("<<ActualWarmTime*1000000000/NumWarmUpdates << " seconds per billion updates per node)\n";
		o << "    Number of directed updates:       \t" << NumDirectedWarmUpdates << "\t("<<ActualWarmTime*1000000000/NumDirectedWarmUpdates << " seconds per billion updates per node)\n";
		o << "    Directed update length:           \t"<< double(NumWarmUpdates)/NumDirectedWarmUpdates<<std::endl;
		o << "    == Measurements   ==\n";
		o << "    Number of creations/annihilations: \t" << NumMeasUpdates << "\t("<<ActualMeasTime*1000000000/NumMeasUpdates << " seconds per billion updates per node)\n";  
		o << "    Number of directed updates:       \t" << NumDirectedMeasUpdates << "\t("<<ActualMeasTime*1000000000/NumDirectedMeasUpdates << " seconds per billion updates per node)\n";  
		o << "    Directed update length:           \t"<< double(NumMeasUpdates)/NumDirectedMeasUpdates<<std::endl;

		o << "    Number of measurements: " << BrokenHistorgram[0] << "\n\n";
		o <<::std::endl;
		o <<"  *******************************\n";
		o <<"  * Broken worldlines histogram *\n";
		o <<"  *******************************\n\n";
		o <<"    N lines\tCount\tProbability\n\n";


		double Normalization=0;
		for(BrokenHistogramType::size_type i=0;i<BrokenHistorgram.size();++i)
			Normalization+=BrokenHistorgram[i];

		for(BrokenHistogramType::size_type i=0;i<BrokenHistorgram.size();++i)
			if(BrokenHistorgram[i]!=0)
				o<<"    "<<i<<"\t\t"<<BrokenHistorgram[i]<<"\t"<<BrokenHistorgram[i]/Normalization<<std::endl;

		o << std::endl;

		o << "  ***********************************************************************************\n";
		o << "  * Energies (obtained from operator string length and Green operator state energy) *\n";
		o << "  ***********************************************************************************\n\n";
		o << "    " << MeasuredOp.TotalEnergy() << "\n";
		o << "    " << MeasuredOp.Potential() << "\n";           
		o << "    " << MeasuredOp.Kinetic() << "\n\n";
		o << "  ******************************\n";
		o << "  * User's defined measurables *\n";
		o << "  ******************************\n\n";

		for(int i=0;i<MeasuredOp.size();++i) 
			o<<"    "<<MeasuredOp.Quantity(i)<<std::endl;

	}
};


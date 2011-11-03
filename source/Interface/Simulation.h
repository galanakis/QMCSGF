#include <fstream>
#include <ctime>

template<class T> T Min(const T &a,const T &b) {return (a<b)?a:b;} 
template<class T> T Max(const T &a,const T &b) {return (a>b)?a:b;} 
 

class FileNameProgressBar {
	 static const int StringLength=200;
	 char Status[StringLength],*Ptr;
	 double Progress; 
	 clock_t Time;
public:
	 FileNameProgressBar(const char *Name) : Progress(0) {
		 Time=clock();
		 strcpy(Status,MathExpression::GetSimulName());
		 Ptr=Status+strlen(Status);   
		 sprintf(Ptr," - 000 - 000:00:00");     
     fstream File;
     File.open(Status,ios::out);
     File.close();    
	 }
	 ~FileNameProgressBar() {
		remove(Status);
	 }
   
	 inline void Update(double prog) {
		if(prog<0) prog=0.0;
		if(prog>1) prog=1.0;
		if(static_cast<int>(100*(prog-Progress))!=0) {
			Progress=prog;
			clock_t Now=clock();
			double DeltaSeconds=double(Now-Time)/CLOCKS_PER_SEC;
			Time=Now;
			unsigned int SecondsLeft=100*(1-Progress)*DeltaSeconds;
			unsigned int HoursLeft=SecondsLeft/3600;
			SecondsLeft%=3600;
			unsigned int MinutesLeft=SecondsLeft/60;
			SecondsLeft%=60;
			
			
		  char OldStatus[StringLength];
      strcpy(OldStatus,Status);
			int percent=static_cast<int>(100*Progress);
			
			sprintf(Ptr," - %.3d - %.3d:%.2d:%.2d",percent,HoursLeft,MinutesLeft,SecondsLeft);
      
			for(int i=0;i<strlen(Status);++i)
				cout<<'\b';  
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
		
		double ActualWarmTime,ActualMeasTime,ITERATIONS_PER_SEC;
		 
		SGF::OperatorStringType &OpString;
		SGF::Measurable &MeasuredOp;
		
public:
		Simulation(SGF::OperatorStringType &OStr, SGF::Measurable &MOp) : OpString(OStr), MeasuredOp(MOp) {} 

    void Thermalize()
        {
          WarmTime=MathExpression::GetValue("#WarmTime").Re();
					WarmIterations=(MathExpression::Find("WarmIterations")!=NULL) ? static_cast<unsigned long long>(MathExpression::GetValue("WarmIterations").Re()) : std::numeric_limits<unsigned long long>::max();
					
					FileNameProgressBar pbar("Thermalizing");

 					const unsigned int BenchmarkIterations=20000;
					const unsigned int AlphaUpdatePeriod=100000;

					clock_t StartTime=clock();
					clock_t EndTime=StartTime+WarmTime*CLOCKS_PER_SEC;

					// Run a few iterations first to determine the cpu speed.
          NumWarmUpdates=0;
					NumDirectedWarmUpdates=0;

					do {
            NumWarmUpdates+=OpString.directed_update();   // Perform an update.
						++NumDirectedWarmUpdates;						
					} while(NumWarmUpdates<BenchmarkIterations);
					
					
				  ITERATIONS_PER_SEC=NumWarmUpdates*CLOCKS_PER_SEC/double(clock()-StartTime);
 				  std::cout<<"Thermalizing for "<<WarmTime<<" seconds or "<<WarmIterations<<" iterations( ~"<<WarmIterations/ITERATIONS_PER_SEC<<" seconds ), whichever finishes first, at "<<ITERATIONS_PER_SEC<<" iterations per second."<<std::endl;


          do {

							NumWarmUpdates+=OpString.directed_update(); // Perform an update. 

							++NumDirectedWarmUpdates;
              if (NumDirectedWarmUpdates%AlphaUpdatePeriod==(AlphaUpdatePeriod-1)) OpString.AlphaUpdate();     


							double Progress=Max(static_cast<double>(NumWarmUpdates)/WarmIterations,static_cast<double>(clock()-StartTime)/(EndTime-StartTime));
 							pbar.Update(Progress);

              
          } while ( NumWarmUpdates<WarmIterations && clock()<EndTime );

          while (OpString.NBrokenLines()!=0)  {              // Perform extra updates until we end up
            NumWarmUpdates+=OpString.directed_update();      // in a diagonal configuration.
						++NumDirectedWarmUpdates;
					}
          
					ActualWarmTime=double(clock()-StartTime)/CLOCKS_PER_SEC;
					ITERATIONS_PER_SEC=NumWarmUpdates/ActualWarmTime;
					std::cout<<"Done Thermalizing after "<<NumWarmUpdates<<" updates ("<<NumDirectedWarmUpdates<<" directed) in "<<ActualWarmTime<<" seconds."<<std::endl;
					
        }
        
      void Measure()
        {
          MeasTime=MathExpression::GetValue("#MeasTime").Re();
					MeasIterations=(MathExpression::Find("MeasIterations")!=NULL) ? static_cast<unsigned long long>(MathExpression::GetValue("MeasIterations").Re()) : std::numeric_limits<unsigned long long>::max();

					FileNameProgressBar pbar("Measuring");
					
					std::cout<<"Measuring for "<<MeasTime<<" seconds or "<<MeasIterations<<" iterations( ~"<<MeasIterations/ITERATIONS_PER_SEC<<" seconds ), whichever finishes first, at "<<ITERATIONS_PER_SEC<<" iterations per second."<<std::endl;

					clock_t StartTime=clock();
					clock_t EndTime=StartTime+MeasTime*CLOCKS_PER_SEC;

          unsigned long long NumBins=MathExpression::GetValue("#Bins").Re();

          NumMeasUpdates=0;
					NumDirectedMeasUpdates=0;

          for (int i=0;i<NumBins;i++)
            {
							unsigned long long counter=0;
              do {
                  counter+=OpString.directed_update();   // Perform an update.
									++NumDirectedMeasUpdates;
                  MeasuredOp.measure(OpString);                 // Perform measurements.
									double Progress=Max(static_cast<double>(NumMeasUpdates+counter)/MeasIterations,static_cast<double>(clock()-StartTime)/(EndTime-StartTime));
		 							pbar.Update(Progress);
              } while ( counter<MeasIterations/NumBins && clock()<EndTime );

							NumMeasUpdates+=counter;
								
              while (OpString.NBrokenLines()!=0) {              // Perform extra updates until we end up
                NumMeasUpdates+=OpString.directed_update();     // in a diagonal configuration.
								++NumDirectedMeasUpdates;
                MeasuredOp.measure(OpString);                   // Perform measurements.
							}
              MeasuredOp.flush();                               // Bin the data. 

            }
          
					ActualMeasTime=double(clock()-StartTime)/CLOCKS_PER_SEC;
					std::cout<<"Done Measuring after "<<NumDirectedMeasUpdates<<" updates ("<<NumDirectedMeasUpdates<<" directed) in "<<ActualMeasTime<<" seconds."<<std::endl;
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
          cout << "    Number of creations/anihilations for thermalization: " << NumWarmUpdates << " ("<<ActualWarmTime*1000000000/NumWarmUpdates << " seconds per billion updates)\n";
          cout << "    Number of directed updates for thermalization: " << NumDirectedWarmUpdates << " ("<<ActualWarmTime*1000000000/NumDirectedWarmUpdates << " seconds per billion updates)\n";
					cout << "    Directed update length for thermalization: "<< double(NumWarmUpdates)/NumDirectedWarmUpdates<<std::endl;
          cout << "    Number of creations/anihilations for measurements: " << NumMeasUpdates << " ("<<ActualMeasTime*1000000000/NumMeasUpdates << " seconds per billion updates)\n";  
          cout << "    Number of directed updates for measurements: " << NumDirectedMeasUpdates << " ("<<ActualMeasTime*1000000000/NumDirectedMeasUpdates << " seconds per billion updates)\n";  
					cout << "    Directed update length for measurements: "<< double(NumMeasUpdates)/NumDirectedMeasUpdates<<std::endl;
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
          
          for (int i=0;i<MeasuredOp.size();i++)
            if (Names[i][0]!='#')
              UserMeasurables++;
          
          if (UserMeasurables)
            {
              cout << "  ******************************\n";
              cout << "  * User's defined measurables *\n";
              cout << "  ******************************\n\n";
          
              for (int i=0;i<MeasuredOp.size();i++)
                if (Names[i][0]!='#')
                  cout << "    " << Names[i] << ": " << MeasuredOp[i] << endl;
            }
        }
  };


#include <fstream>

class Simulation
  {
    static int Progress;
    static long long NumWarmUpdates,NumMeasUpdates;
    static double StatusLength,ClockLength,WarmTime,MeasTime;
    static time_t StatusStart,ClockStart;
    static char Status[100],*Ptr;
    
    static void InitStatus(const char *Name,double T)
      {
        strcpy(Status,MathExpression::GetSimulName());
        strcpy(Status+strlen(Status)," - ");
        strcpy(Status+strlen(Status),Name);
        strcpy(Status+strlen(Status)," 000 %");
        Ptr=Status+strlen(Status)-5;
        StatusLength=T;
        Progress=-1;
        time(&StatusStart);
        fstream File;
        File.open(Status,ios::out);
        File.close();
      }
      
    static void InitClock(double T)
      {
        time(&ClockStart);
        ClockLength=T;
      }
      
    static inline void UpdateStatus(void)
      {
        static time_t End;
        static double Time;
        time(&End);
        Time=difftime(End,StatusStart);
        
        if (Progress!=(int) (100.0*Time/StatusLength))
          {
            static char OldStatus[100];
            strcpy(OldStatus,Status);
            Progress++;
            Ptr[0]=Progress/100+48;
            Ptr[1]=(Progress/10)%10+48;
            Ptr[2]=Progress%10+48;
            rename(OldStatus,Status);
          }
      }
      
    static inline bool KeepWorking(void)
      {
        static time_t End;
        time(&End);
        return (difftime(End,ClockStart)<ClockLength);
      }
      
    public:
      static void Thermalize(SGF::OperatorStringType &OpString)
        {
          WarmTime=MathExpression::GetValue("#WarmTime").Re();
          NumWarmUpdates=0;
          InitStatus("Thermalizing",WarmTime);
          InitClock(WarmTime);
          
          int Iteration=0;
          
          do
            {
              UpdateStatus();                               // Update the percentage of progress, if needed.
              NumWarmUpdates+=OpString.directed_update();   // Perform an update.
              
              if (++Iteration%100000==99999)
                OpString.AlphaUpdate(); 
              
            } while (KeepWorking());

          while (OpString.NBrokenLines()!=0)                // Perform extra updates until we end up
            NumWarmUpdates+=OpString.directed_update();     // in a diagonal configuration.
          
          OpString.print_alphas();
          remove(Status);
        }
        
      static void Measure(SGF::OperatorStringType &OpString,SGF::Measurable &MeasuredOp)
        {
          int NumBins=MathExpression::GetValue("#Bins").Re();
          MeasTime=MathExpression::GetValue("#MeasTime").Re();
          NumMeasUpdates=0;
          InitStatus("Measuring",MeasTime);

          for (int i=0;i<NumBins;i++)
            {
              InitClock(MeasTime/NumBins);
              
              do
                {
                  UpdateStatus();                               // Update the percentage of progress, if needed.
                  NumMeasUpdates+=OpString.directed_update();   // Perform an update.
                  MeasuredOp.measure(OpString);                 // Perform measurements.
                } while (KeepWorking());

              while (OpString.NBrokenLines()!=0)                // Perform extra updates until we end up
                NumMeasUpdates+=OpString.directed_update();     // in a diagonal configuration.
          
              MeasuredOp.flush();                               // Bin the data.
            }
            
          remove(Status);
        }
        
      static void Results(SGF::Measurable &MeasuredOp)
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
          cout << "    Number of updates for thermalization: " << NumWarmUpdates << " ("<<WarmTime*1000000000/NumWarmUpdates << " seconds per billion updates)\n";
          cout << "    Number of updates for measurements: " << NumMeasUpdates << " ("<<MeasTime*1000000000/NumMeasUpdates << " seconds per billion updates)\n";  
          cout << "    Number of measurements: " << MeasuredOp.count() << "\n\n";
          MeasuredOp.print_histogram();
          cout << endl;
          cout << "  ***********************************************************************************\n";
          cout << "  * Energies (obtained from operator string length and Green operator state energy) *\n";
          cout << "  ***********************************************************************************\n\n";
          cout << "    Total energy: " << MeasuredOp.TotalEnergy() << "\n";
          cout << "    Diagonal energy: " << MeasuredOp.PotentialEnergy() << "\n";           
          cout << "    Non-diagonal energy: " << MeasuredOp.KineticEnergy() << "\n\n";
          cout << "  ******************************\n";
          cout << "  * User's defined measurables *\n";
          cout << "  ******************************\n\n";
          std::vector<string> Names=MathExpression::GetMeasurableList();
          
          for (int i=0;i<MeasuredOp.size();i++)
            cout << "    " << Names[i] << ": " << MeasuredOp[i] << endl;
        }
  };

char Simulation::Status[100];
char *Simulation::Ptr;
int Simulation::Progress;
long long Simulation::NumWarmUpdates;
long long Simulation::NumMeasUpdates;
double Simulation::StatusLength;
double Simulation::ClockLength;
double Simulation::WarmTime;
double Simulation::MeasTime;
time_t Simulation::StatusStart;
time_t Simulation::ClockStart;

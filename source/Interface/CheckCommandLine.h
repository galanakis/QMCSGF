#ifndef _CHECK_COMMAND_LINE
#define _CHECK_COMMAND_LINE

#define Version "1.0 (May 10th, 2011)"

#include <ScrEx.h>

#include <fstream>
using std::fstream;
using std::ios;

extern std::ostream cout;

int CheckCommandLine(int NumArg,char **Arg)
  {
    fstream File;
	
    switch (NumArg)
      {
	case 2:
	  File.open(Arg[1],ios::in);
	  
	  if (!File.is_open())
	    {
	      cout << "\nFile " << Arg[1] << " does not exist!\n\n";
	      return 0;
	    }
	    
	  else
	    {
	      File.close();
	      return 1;
	    }
	  
	case 3:
	  if (!strcmp("-CreateExample",Arg[1]))
	    {
	      strcpy(ScriptExample+625,Version);
	      ScriptExample[625+strlen(Version)]=' ';
	      File.open(Arg[2],ios::out);
	      File.write(ScriptExample,strlen(ScriptExample));
	      File.close();
	      cout << "\nFile " << Arg[2] << " has been created!\n\n";
	      return 0;
	    }
	    
	  else if (!strcmp("-Hamiltonian",Arg[1]))
	    {
	      File.open(Arg[2],ios::in);

	      if (!File.is_open())
		{
		  cout << "\nFile " << Arg[2] << " does not exist!\n\n";
		  return 0;
		}
	    
	      else
		{
		  File.close();
		  return 2;
		}
	    }
	    
	  else
	    {
	      cout << "\nError in command line!\nExecute the program with no parameter to display help.\n\n";
	      return 0;
	    }
	    
	case 4:
	  if (!strcmp("-Display",Arg[1]))
	    {
	      File.open(Arg[3],ios::in);

	      if (!File.is_open())
		{
		  cout << "\nFile " << Arg[3] << " does not exist!\n\n";
		  return 0;
		}
	    
	      else
		{
		  File.close();
		  return 3;
		}
	    }
	    
	  else
	    {
	      cout << "\nError in command line!\nExecute the program with no parameter to display help.\n\n";
	      return 0;
	    }
	
	default:
	  cout << "\n*****************************************************************************************************\n";
	  cout << "* This is the \"Corvette SGF engine\" for the quantum Monte Carlo simulation of lattice Hamiltonians. *\n";
          cout << "*                                                                                                   *\n";
          cout << "* For informations on the SGF algorithm:                                                            *\n";
          cout << "*   Physical Review E 77, 056705 (2008)                                                             *\n";
          cout << "*   Physical Review E 78, 056707 (2008)                                                             *\n";
          cout << "*                                                                                                   *\n";
	  cout << "* Dr Valy G. Rousseau and Dr Dimitris Galanakis - Version " << Version;
	  for (unsigned int i=0;i<42-strlen(Version);i++) cout << " ";
	  cout << "*\n";
	  cout << "*****************************************************************************************************\n\n";
	  cout << "USAGE: " << Arg[0] << " InputFile > OutputFile &\n\n";
	  cout << "The input file must be an SGF script that defines the Hamiltonian and the parameters of the simulation.\n\n";
	  cout << "SPECIAL USAGE:\n\n";
	  cout << "An example of SGF script file can be created with the command:\n";
	  cout << "  " << Arg[0] << " -CreateExample ExampleFile\n\n";
	  cout << "The list of Hamiltonian terms can be displayed with the command:\n";
	  cout << "  " << Arg[0] << " -Hamiltonian InputFile\n\n";
	  cout << "Any quantity defined by the user can be displayed with the command:\n";
	  cout << "  " << Arg[0] << " -Display Quantity InputFile\n\n";
	  return 0;
      }
  }

#endif
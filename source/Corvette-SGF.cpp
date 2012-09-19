// *****************************************************************************************************
// * This is the "Corvette SGF engine" for the quantum Monte Carlo simulation of lattice Hamiltonians. *
// *                                                                                                   *
// * For informations on the SGF algorithm:                                                            *
// *   Physical Review E 77, 056705 (2008)                                                             *
// *   Physical Review E 78, 056707 (2008)                                                             *
// *                                                                                                   *
// * Dr Valy G. Rousseau and Dr Dimitris Galanakis                                                     *
// * _________________________________________________________________________________________________ *
// * Changes made to previous version:                                                                 *
// *   None, this is the first version!                                                                *
// *                                                                                                   *
// * Known bugs:                                                                                       *
// *   Come on! Bugs do not exist in an SGF engine!!!                                                  *
// *****************************************************************************************************

#define Version "1.0 (May 10th, 2011)"

// **********************
// * Standard libraries *
// **********************

#include <iostream>
#include <cstring>
#include <limits>


#ifdef USEMPI
// *******************************************
// * Declaration of global variables for MPI *
// *******************************************
#include <mpi.h>
#define Master 0
int NumProcessors,Rank,NameLength;
char ProcessorName[MPI_MAX_PROCESSOR_NAME];

#endif

std::ostream cout(std::cout.rdbuf());
using std::endl;
using std::string;
using std::fstream;
using std::ios;
using std::streamoff;
using std::ostream;
using std::numeric_limits;

// **********************
// * Personal libraries *
// **********************

#include <Parser.h>
#include <MathExpression.h>
#include <ParserSGF.h>
#include <ScrEx.h>
#include <CheckCommandLine.h>
#include <SGFContainer.h>
#include <OperatorString.hh>
#include <Measurable.hh>
#include <Simulation.h> 


void Simulator() {

	SGFContainer Container;
	int GreenOperatorLines=Container.GreenOperatorLines();

	const SGF::Hamiltonian &T=Container.Kinetic();
	const SGF::Hamiltonian &V=Container.Potential();
	double Beta=Container.Beta();
	double AlphaParameter=Container.AlphaParameter();
	SGF::GreenOperator<long double> g(Container.NSites(),GreenOperatorLines);


	SGF::OperatorStringType OperatorString(T,V,Beta,g,AlphaParameter);

	/* Initializing the simulation. Thermalize, Measure and pring the results */
	Simulation simul;

	// We start warm up iterations
	simul.Thermalize(OperatorString);

	// This defines the measurable objects some of which delay updates even if not measured.
	// This is why I declare the measurable operators after the thermalization.
	SGF::Measurable MeasuredOperators(OperatorString);
	MeasuredOperators.insert(MathExpression::GetMeasurableList(),Container.MeasurableOperators());
  
	//We start measurement iterations
	simul.Measure(OperatorString,MeasuredOperators);
  
	// We diplay the results of the simulation
	simul.Results(MeasuredOperators);

}

int finish() {
#ifdef USEMPI
	MPI_Finalize();
#endif
	return 0;
}

// ***************************
// * Here starts the program *
// ***************************

int main(int NumArg,char **Arg)
  {

#ifdef USEMPI
    // ***********************************
    // * Initialization of MPI functions *
    // ***********************************
    
    MPI_Init(&NumArg,&Arg);
    MPI_Comm_size(MPI_COMM_WORLD,&NumProcessors);
    MPI_Comm_rank(MPI_COMM_WORLD,&Rank);
    MPI_Get_processor_name(ProcessorName,&NameLength);

		// Silence the other nodes. Only the root node can print
		cout.rdbuf( Rank==Master ? std::cout.rdbuf() : 0 );
 
#endif

    // *********************************************************************************
    // * We check the command line. If it is not correct, a help message is displayed. *
    // *********************************************************************************
    
    char *File=NULL;
    
    switch (CheckCommandLine(NumArg,Arg))
      {
        case 0: return finish();        // Simulation is aborded.
        case 1: File=Arg[1]; break;     // Command line is correct and requests to start the simulation.
        case 2: File=Arg[2]; break;     // Command line is correct and requests to display Hamiltonian terms.
        case 3: File=Arg[3]; break;     // Command line is correct and requests to display a specified quantity.
      }
    
    // ******************************************************************
    // * We read the input file and transform it into a list of tokens. *
    // ******************************************************************
    
    char *Error;

    if ((Error=Parser::ReadFile(File)))
      {
        cout << Error << endl;
        return finish();
      }
      
    // ********************************************************
    // * We transform the tokens into a list of SGF commands. *
    // ********************************************************

    ParserSGF ScriptSGF;
    
    if ((Error=ScriptSGF.ReadTokens()))
      {
        cout << Error << endl;
        return finish();
      }

    MathExpression::Initialize();	// We initialize the table of symbols with keywords and constants.
    
    // *******************************
    // * We execute the SGF commands *
    // *******************************
    
    if ((Error=ScriptSGF.ExecuteCommands()))
      {
        cout << Error << endl;
        return finish();
      }
    
    // ************************************************
    // * We build a suitable form of the Hamiltonian. *
    // ************************************************
    
    if ((Error=MathExpression::BuildHamiltonian()))
      {
        cout << Error << endl;
        return finish();
      }
      
    if (File==Arg[1])
      {
        // *************************
        // * Start the simulation. *
        // *************************

        if ((Error=MathExpression::SignOrPhaseProblem()))
          {
            cout << Error << endl;      // A sign or a phase problem has been detected. A warning message
            return finish();                   // is displayed and the simulation is aborded.
          }

					Simulator();

      }
      
    else if (File==Arg[2])
      MathExpression::ListHamiltonianTerms();                   // Display hamiltonian terms and abord the simulation.
      
    else
      {
        // ************************************************************************
        // * Display the quantity specified by the user and abord the simulation. *
        // ************************************************************************

        string StrError=MathExpression::DisplayQuantity(Arg[2]);
        
        if (StrError.length()!=0)
          cout << StrError << endl;
      }

    return finish();
  }

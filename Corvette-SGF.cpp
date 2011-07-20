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

using namespace std;

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

// ***************************
// * Here starts the program *
// ***************************

int main(int NumArg,char **Arg)
  {
    // *********************************************************************************
    // * We check the command line. If it is not correct, a help message is displayed. *
    // *********************************************************************************
    
    char *File=NULL;
    
    switch (CheckCommandLine(NumArg,Arg))
      {
        case 0: return 0;               // Simulation is aborded.
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
        return 0;
      }
      
    // ********************************************************
    // * We transform the tokens into a list of SGF commands. *
    // ********************************************************

    ParserSGF ScriptSGF;
    
    if ((Error=ScriptSGF.ReadTokens()))
      {
        cout << Error << endl;
        return 0;
      }

    MathExpression::Initialize();	// We initialize the table of symbols with keywords and constants.
    
    // *******************************
    // * We execute the SGF commands *
    // *******************************
    
    if ((Error=ScriptSGF.ExecuteCommands()))
      {
        cout << Error << endl;
        return 0;
      }
    
    // ************************************************
    // * We build a suitable form of the Hamiltonian. *
    // ************************************************
    
    if ((Error=MathExpression::BuildHamiltonian()))
      {
        cout << Error << endl;
        return 0;
      }
      
    if (File==Arg[1])
      {
        // *************************
        // * Start the simulation. *
        // *************************

        if ((Error=MathExpression::SignOrPhaseProblem()))
          {
            cout << Error << endl;      // A sign or a phase problem has been detected. A warning message
            return 0;                   // is displayed and the simulation is aborded.
          }

// *****************************
// * Dimitris's temporary code *
// *****************************
// _____________________________________________________________________
SGFContainer Container;
SGF::RebuildPeriod=MathExpression::GetValue("RebuildPeriod").Re();
int GreenOperatorLines=MathExpression::GetValue("GreenOperatorLines").Re();

const SGF::Hamiltonian &T=Container.Kinetic();
const SGF::Hamiltonian &V=Container.Potential();
double Beta=Container.Beta();
SGF::OperatorStringType OperatorString(T,V,Beta);
OperatorString.alpha(SGF::ADD)=0.95;
OperatorString.alpha(SGF::REMOVE)=0.95;
OperatorString.GreenInit(Container.NSites(),GreenOperatorLines);
SGF::Measurable MeasuredOperators(Container.MeasurableOperators());
// _____________________________________________________________________
    
        Simulation::Thermalize(OperatorString);                 // We start warm up iterations.
        Simulation::Measure(OperatorString,MeasuredOperators);  // We start measurement iterations.
        Simulation::Results(MeasuredOperators);                 // We display the results of simulation.
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

    return 0;
  }

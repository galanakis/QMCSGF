char ScriptExample[]="\
# **************************************************************************************************************\n\
# * This file provides an example of SGF script that can be used as an input file for the Corvette SGF engine. *\n\
# * An SGF script consists in a set of instructions, each being followed by a semicolon.                       *\n\
# * CAUTION: SGF script language is case-sensitive.                                                            *\n\
# *                                                                                                            *\n\
# * Dr Valy G. Rousseau and Dr Dimitris Galanakis - Version                                                    *\n\
# **************************************************************************************************************\n\
\n\
# In this example, we define a Hamiltonian that describes atoms and diatomic molecules on a one-dimensional lattice.\n\
# Atoms and molecules can hop from one site to a neighboring site. They can interact via onsite repulsion potentials,\n\
# and two atoms can be turned into a molecule (and vice-versa).\n\
\n\
# Variables can be defined in order to perform calculations, or simply to clarify the expressions that the user wants to define.\n\
# For instance, the algorithm needs to know the inverse temperature. However it is usually more convenient for the user to specify\n\
# the temperature. Thus we will tell the code to invert this value.\n\
\n\
  Variable T=0.1;                    # This is the temperature. The algorithm needs the inverse of T (see below).\n\
  Variable ta=1.0;                   # Hopping parameter for atoms.\n\
  Variable tm=0.5;                   # Hopping parameter for molecules.\n\
  Variable Uaa=8.0;                  # Onsite potential between atoms.\n\
  Variable Umm=20.0;                 # Onsite potential between molecules.\n\
  Variable Uam=5.0;                  # Onsite potential between atoms and molecules.\n\
  Variable g=0.3;                    # Conversion parameter between atoms and molecules.\n\
  Variable L=20;                     # This will be the number of lattice sites.\n\
\n\
# Of course a variable can be overwritten.\n\
\n\
  g=2(g+0.1)-0.2;                    # Note that the multiplication symbol * can be omitted when there is no ambiguity.\n\
\n\
# All variables defined above are real, but variables can also be complex. In order to define complex numbers, the predefined\n\
# variable 'SquareRootOfMinusOne' can be used. As its name suggests, it is equal the squareroot of -1.\n\
\n\
  Variable i=SquareRootOfMinusOne;   # For convenience, we define i as the usual notation for the squareroot of -1.\n\
  Variable Aux=5+3i;                 # We define a general complex number, just for fun. It is completely useless in this example.\n\
\n\
# The boundary conditions and the lattice must be defined. Here we define a one-dimensional lattice with 20 sites. Two-dimensional and three-dimensional lattices\n\
# can be defined as well by specifying the number of sites in each direction, separated by comas.\n\
\n\
  Boundaries Periodic;               # The boundaries can be either 'Periodic' or 'Open'.\n\
  Lattice L;                         # We define a one-dimensional lattice with 20 sites.\n\
\n\
# One species of particles at least must be defined. A species is defined by a name, identifiers for the creation, annihilation, and number operators, and\n\
# the maximum occupation number per site. This latter can be finite or set to infinity by using the predefined variable 'Infinity'.\n\
\n\
  Species \"Atoms\",Ca,Aa,Na,Infinity; # We define a species named \"Atoms\" where Ca, Aa, Na are the creation, annihilation, and number operators, with an infinite maximum occupation number.\n\
  Species \"Molecules\",Cm,Am,Nm,5;    # We define a species named \"Molecules\" where Cm, Am, Nm are the creation, annihilation, and number operators, with a maximum occupation number of 5.\n\
\n\
# The initial population of each species can be defined (Default value is 0 if not specified).\n\
\n\
  Population \"Atoms\",L/2+3;          # We set the initial population of atoms to 3 particles above half-filling.\n\
  Population \"Molecules\",0;          # We set the initial population of molecules to 0.\n\
\n\
# In order to ease the definition of the Hamiltonian, operators can be defined. Operators can be expressed in terms of creation, annihilation, and number operators, and can have\n\
# complex coefficients. As for variables, operators can be overwritten.\n\
\n\
  Operator Ka=-ta Sum{<p,q>}(Ca[p]Aa[q]+Ca[q]Aa[p]);           # We define the kinetic energy of atoms.\n\
  Operator Km=-tm Sum{<p,q>}(Cm[p]Am[q]+Cm[q]Am[p]);           # This is the kinetic energy of molecules.\n\
  Operator Paa=Uaa/2 Sum{p}(Na[p](Na[p]-1));                   # We define the interaction potential between atoms.\n\
  Operator Pmm=Umm/2 Sum{p}(Nm[p](Nm[p]-1));                   # This is the interaction potential between molecules.\n\
  Operator Pam=Uam Sum{p}(Na[p]Nm[p]);                         # This is the interaction potential between atoms and molecules.\n\
  Operator G=-g Sum{p}(Ca[p]Ca[p]Am[p]+Cm[p]Aa[p]Aa[p]);       # We define the conversion term that turns atoms into molecules, and vice-versa.\n\
  Operator H=Ka+Km+Paa+Pmm+Pam+G;                              # Finally we define the total energy.\n\
\n\
# The Hamiltonian must be defined. This is easily done thanks to the above definitions of operators. The SGF algorithm requires all diagonal terms of the Hamiltonian to be real and all\n\
# non-diagonal terms to be negative, after expansion and simplification. If those requirements are not fulfilled, a warning message is displayed and the simulation aborded.\n\
\n\
  Hamiltonian H;                     # We define the operator H as our Hamiltonian.\n\
\n\
# The inverse temperature, the seed of the random number generator, the time for thermalization, the time for measurements, the number of bins, and the ensemble of the\n\
# simulation must be defined.\n\
\n\
  InverseTemperature 1/T;            # We define the inverse temperature as the inverse of the temperature T (no kidding).\n\
  Seed 34715;                        # The seed of the random number generator can be any non-zero integer.\n\
  WarmTime 30*Minute+17*Second;      # The time for thermalization is in unit of seconds. For convenience we can use the predefined variables Second, Minute, Hour, Day, and Week.\n\
  MeasTime Hour;                     # In the same way, this is the time for measurements.\n\
  Bins 20;                           # The number of bins is used to estimate the statistical errors.\n\
  Ensemble Canonical;                # The ensemble can be either 'Canonical' or 'GrandCanonical'.\n\
\n\
# Some generic quantities are systematically measured: Diagonal energy, non-diagonal energy, total energy, local density, density-density correlations, and one-body Green functions.\n\
# The user can also specify additional measurements. Any operator defined by the user can be measured.\n\
\n\
  Measure \"Atom Kinetic energy\",Ka;                            # A measurable is defined by a name and an expression.\n\
  Measure \"Number of atoms\",Sum{p}(Na[p]);                     # The name of the measurable will be printed in the output file.\n\
  Measure \"Number of molecules\",Sum{p}(Nm[p]);                 # Caution: Do not measure a quantity that you don't need!\n\
\n\
# A name must be given to the simulation.\n\
\n\
  SimulName \"SGF ---> The ass-kicking QMC algorithm\";     # This name will appear in a status file in order to easily identify the simulations that are running.\n\
\n\
# ***************\n\
# * Useful tips *\n\
# ***************\n\
\n\
# Some variables are predefined:\n\
#   Second: 1\n\
#   Minute: 60\n\
#   Hour: 3600\n\
#   Day: 86400\n\
#   Week: 604800\n\
#   Infinity: +inf\n\
#   Pi: 3.1415926535897932384626433832795\n\
#   SquareRootOfMinusOne: Yes, this is the square root of -1\n\
# 'Sum{i,Min,Max}(Expression)' performs the sum of Expression over variable i from Min to Max.\n\
# 'Sum{i}(Expression)' performs the sum of Expression over all lattice sites.\n\
# 'Sum{<i,j>}(Expression)' performs the sum of Expression over distinct pairs of first neighoring sites i,j.\n\
# 'Prod{i,Min,Max}(Expression)' performs the product of Expression over variable i from Min to Max.\n\
# 'Prod{i}(Expression)' performs the product of Expression over all lattice sites.\n\
# 'Prod{<i,j>}(Expression)' performs the product of Expression over distinct pairs of first neighoring sites i,j.\n\
# 'Random' returns a random number chosen with a uniform distribution in [0;1[. Can be used for disorder. CAUTION: The seed should be declared first!\n\
# 'Pos?(Site)' with ?={X,Y,Z} returns the ? coordinate of Site with respect to the center of the lattice. Useful for defining confining potentials.\n\
# Some usual functions can be used:\n\
#   sin, cos, tan, exp, ln, sqrt. These functions accept and return complex numbers.\n\
#   Fact(n). This function currently accepts positive integers only. It will be generalized later to complex numbers.\
";

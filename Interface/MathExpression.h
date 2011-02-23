// *****************************************************************************
// * This class defines an object that represents a mathematical expression of *
// * numbers, literals, functions, and operators. The class can also store a   *
// * table of symbols.                                                         *
// * _________________________________________________________________________ *
// * Dr. Valy G. Rousseau - December 23rd, 2010                                *
// *****************************************************************************
//
// An expression is represented by a binary tree. The nodes of the tree can have the following types:
//   - MathExpression::ArithmeticOperatorX where X={1,2,3,4,5,6} corresponding to {+,-,*,/,%,^}
//   - MathExpression::ComplexNumber
//   - MathExpression::Function
//   - MathExpression::Operator
//   - MathExpression::Literal
//
// Mathematical expressions are treated by using the following grammar:
//   E       -> T | E + T | E - T
//   T       -> F | T * F | T / F | T % F | T ^ F
//   F       -> C | f(E) | O[E] | S | (E) | Sum(E) | Prod(E)
//   Sum(E)  -> E | E + Sum(E)
//   Prod(E) -> E | E + Prod(E)
// with
//   E = Expression
//   T = Term
//   F = Factor
//   C = Complex number
//   f = Function
//   O = Operator
//   S = Symbol (Name of an expression)
//   Sum(E) = Sum of expression E over some index
//   Prod(E) = Product of expression E over some index

#ifndef MathExpressionDef
#define MathExpressionDef

#include <Complex.h>
#include <Parser.h>
#include <RandomNumberGenerator.hh>
#include <TheClassWithNoName.h>
#include <Factorial.h>
#include <string>
#include <vector>

class MathExpression
  {
    enum Types 
      {
	ArithmeticOperator1,ArithmeticOperator2,ArithmeticOperator3,ArithmeticOperator4,ArithmeticOperator5,ArithmeticOperator6,
	CreationOperator,AnnihilationOperator,NumberOperator,ComplexNumber,Function,Operator,Keyword,Constant,Variable,Literal,
	Sin,Cos,Tan,Exp,Ln,Sqrt,Fact,PosX,PosY,PosZ
      };
		
    class Symbol;
    
    class Node
      {
        short Type;
        Complex Value;
	      Symbol *Sym;
        Node *Son1,*Son2;
        
        public:
          Node() {Sym=NULL; Son1=NULL; Son2=NULL;}
          ~Node() {if (Son1) delete Son1; if (Son2) delete Son2;}
	  
	  static void Copy(Node *&Target,Node *Source)
	    {
	      if (Source)
		{
		  Target=new Node;
		  Target->Type=Source->Type;
		  Target->Value=Source->Value;
		  Target->Sym=Source->Sym;
		  Copy(Target->Son1,Source->Son1);
		  Copy(Target->Son2,Source->Son2);
		}
	    }
	    
	  static bool Proportional(Node *X,Node *Y)
	    {
	      if (X->Type==Y->Type)
		{
		  switch (X->Type)
		    {
		      case ComplexNumber:
			if (X->Value==0.0)
			  return false;
			else
			  return true;
			
		      case CreationOperator:
		      case AnnihilationOperator:
			if (X->Son1->Value==Y->Son1->Value)
			  return true;
			else
			  return false;
			
		      case ArithmeticOperator3:
			if (Proportional(X->Son2,Y->Son2))
			  return Proportional(X->Son1,Y->Son1);
			else
			  return false;
			
		      default:
			cout << "Unknown type in static bool Proportional(Node *,Node *) function.\n";
			return false;
		    }
		}
		
	      return false;
	    }
	  
	  void Substitute()
	    {
	      if (Type==Literal)
		{
		  Copy(Son1,Sym->Expression->Root->Son1);
		  Copy(Son2,Sym->Expression->Root->Son2);
		  Type=Sym->Expression->Root->Type;
		  Value=Sym->Expression->Root->Value;
		  Sym=Sym->Expression->Root->Sym;
		  Substitute();
		}
		
	      else
		{
		  if (Son1) Son1->Substitute();
		  if (Son2) Son2->Substitute();
		}
	    }
	    
	  void FactorizeOperators(Node *Root)
	    {
	      if (Type==ArithmeticOperator1)
		{
		  Son1->FactorizeOperators(Root);
		  Son2->FactorizeOperators(Root);
		}
		
	      else
		{
		  Node *Temp=Root->FindFactor(this);
		  
		  if (Temp && Temp!=this)
		    {
		      Temp->SetFactor(Temp->GetFactor()+this->GetFactor());
		      this->SetFactor(0.0);
		    }
		}
	    }
	    
	  void SetFactor(Complex N)
	    {
	      if (Type==ComplexNumber)
		  Value=N;
	      else
		Son1->SetFactor(N);
	    }
	    
	  Complex GetFactor()
	    {
	      if (Type==ComplexNumber)
		  return Value;
	      else
		return Son1->GetFactor();
	    }
	    
	  Node *FindFactor(Node *Source)
	    {
	      Node *Factor;
	      
	      if (Type==ArithmeticOperator1)
		{
		  if ((Factor=Son1->FindFactor(Source)))
		    return Factor;
		  else if ((Factor=Son2->FindFactor(Source)))
		    return Factor;
		  else
		    return NULL;
		}
		
	      else if (Proportional(this,Source))
		return this;
	      else
		return NULL;
	    }
	    
	  void Expand(Node *Root)
	    {
	      static bool Test;
	      
	      if (Root==this)
		Test=false;
	      
	      if (Type==ArithmeticOperator3)
		{
		  if (Son1->Type==ArithmeticOperator1 || Son1->Type==ArithmeticOperator2)
		    {
		      Node *New=new Node;
		      New->Type=Type;
		      New->Son1=Son1->Son2;
		      New->Son2=Son2;
		      Son2=New;
		      Type=Son1->Type;
		      Son1->Type=Son2->Type;
		      Copy(Son1->Son2,Son2->Son2);
		      Test=true;
		    }
		    
		  else if (Son2->Type==ArithmeticOperator1 || Son2->Type==ArithmeticOperator2)
		    {
		      Node *New=new Node;
		      New->Type=Type;
		      New->Son1=Son1;
		      New->Son2=Son2->Son1;
		      Son1=New;
		      Type=Son2->Type;
		      Son2->Type=Son1->Type;
		      Copy(Son2->Son1,Son1->Son1);
		      Test=true;
		    }
		}
		
	      else if (Type==ArithmeticOperator4)
		{
		  if (Son1->Type==ArithmeticOperator1 || Son1->Type==ArithmeticOperator2)
		    {
		      Node *New=new Node;
		      New->Type=Type;
		      New->Son1=Son1->Son2;
		      New->Son2=Son2;
		      Son2=New;
		      Type=Son1->Type;
		      Son1->Type=Son2->Type;
		      Copy(Son1->Son2,Son2->Son2);
		      Test=true;
		    }
		}
		
	      else if (Type==ArithmeticOperator2)
		{
		  if (Son1->Type==ComplexNumber && Son1->Value==0)
		    {
		      Type=ArithmeticOperator3;
		      Son1->Value=-1;
		    }
		    
		  else
		    {
		      Node *New=new Node;
		      New->Type=ArithmeticOperator3;
		      New->Son1=new Node;
		      New->Son1->Type=ComplexNumber;
		      New->Son1->Value=-1;
		      New->Son2=Son2;
		      Son2=New;
		      Type=ArithmeticOperator1;
		    }
		  
		  Test=true;
		}
		
	      if (Type!=CreationOperator && Type!=AnnihilationOperator && Type!=NumberOperator)
		{
		  if (Son1) Son1->Expand(Root);
		  if (Son2) Son2->Expand(Root);
		}
		
	      if (Test && Root==this)
		Root->Expand(Root);
	    }
	  
	  void ExpandNumberOperator()
	    {
	      if (Type==NumberOperator)
		{
		  Node *New;
		  Copy(New,this);
		  New->Type=AnnihilationOperator;
		  Son2=New;
		  Type=ArithmeticOperator3;
		  New=new Node;
		  New->Son1=Son1;
		  Son1=New;
		  Son1->Type=CreationOperator;
		}
		
	      else if (Type!=CreationOperator && Type!=AnnihilationOperator)
		{
		  if (Son1) Son1->ExpandNumberOperator();
		  if (Son2) Son2->ExpandNumberOperator();
		}
	    }
	    
	  void SimplifyOperatorArgument()
	    {
	      if (Type==CreationOperator || Type==AnnihilationOperator || Type==NumberOperator)
		{
		  Complex Arg=Son1->Evaluate();
		  delete Son1;
		  Son1=new Node;
		  Son1->Type=ComplexNumber;
		  Son1->Value=Arg;
		}
		
	      else
		{
		  if (Son1) Son1->SimplifyOperatorArgument();
		  if (Son2) Son2->SimplifyOperatorArgument();
		}
	    }
	    
	  void LeftProduct(Node *Root)
	    {
	      if (Type==ArithmeticOperator1)
		{
		  Son1->LeftProduct(Son1);
		  Son2->LeftProduct(Son2);
		}
		
	      else if (Type==ArithmeticOperator3)
		{
		  if (Son2->Type==ArithmeticOperator3 || Son2->Type==ArithmeticOperator4)
		    {
		      Type=Son2->Type;
		      Son2->Type=ArithmeticOperator3;
		      Node *Temp=Son1;
		      Son1=Son2;
		      Son2=Son2->Son2;
		      Son1->Son2=Son1->Son1;
		      Son1->Son1=Temp;
		      Root->LeftProduct(Root);
		    }
		    
		  else
		    Son1->LeftProduct(Root);
		}
		
	      else if (Type==ArithmeticOperator4)
		Son1->LeftProduct(Root);
	    }
	    
	  void LeftNumbers(Node *Root)
	    {
	      Node *Temp;
	      
	      if (Type==ArithmeticOperator1)
		{
		  Son1->LeftNumbers(Son1);
		  Son2->LeftNumbers(Son2);
		}
		
	      else if (Type==ArithmeticOperator4)
		{
		  if (Son1->Type==ArithmeticOperator3)
		    {
		      Temp=Son2;
		      Son2=Son1->Son2;
		      Son1->Son2=Temp;
		      Type=ArithmeticOperator3;
		      Son1->Type=ArithmeticOperator4;
		      Root->LeftNumbers(Root);
		    }

		  else if (Son1->Type==CreationOperator || Son1->Type==AnnihilationOperator)
		    {
		      Temp=new Node;
		      Temp->Type=ArithmeticOperator4;
		      Temp->Son1=new Node;
		      Temp->Son1->Type=ComplexNumber;
		      Temp->Son1->Value=1.0;
		      Temp->Son2=Son2;
		      Son2=Son1;
		      Son1=Temp;
		      Type=ArithmeticOperator3;
		      Root->LeftNumbers(Root);
		    }
		    
		  else
		    Son1->LeftNumbers(Root);
		}
		
	      else if (Type==ArithmeticOperator3 && Son2->Type!=CreationOperator && Son2->Type!=AnnihilationOperator)
		{
		  if (Son1->Type==ArithmeticOperator3 && (Son1->Son2->Type==CreationOperator || Son1->Son2->Type==AnnihilationOperator))
		    {
		      Temp=Son2;
		      Son2=Son1->Son2;
		      Son1->Son2=Temp;
		      Root->LeftNumbers(Root);
		    }
		    
		  else if (Son1->Type==CreationOperator || Son1->Type==AnnihilationOperator)
		    {
		      Temp=Son2;
		      Son2=Son1;
		      Son1=Temp;
		      Root->LeftNumbers(Root);
		    }
		    
		  else
		    Son1->LeftNumbers(Root);
		}
		
	      else if (Type==ArithmeticOperator3)
		Son1->LeftNumbers(Root);
	    }
	    
	  void NormalOrder(Node *&Root,Node *&Father)
	    {
	      if (Type==ArithmeticOperator1)
		{
		  Son1->NormalOrder(Son1,Son1);
		  Son2->NormalOrder(Son2,Son2);
		}
		
	      else if (Type==ArithmeticOperator3)
		{
		  Node *Temp;
	      
		  if (Son1->Type==ArithmeticOperator3 && Son1->Son2->Type==AnnihilationOperator && Son2->Type==CreationOperator)
		    {
		      if (Son1->Son2->Son1->Evaluate()==Son2->Son1->Evaluate())
			{
			  // The creation and annihilation operators do not commute.
			  
			  if (Root==this)
			    {
			      // The tree is of the form (?*A)*C
		      
			      Temp=new Node;
			      Temp->Type=ArithmeticOperator1;
			      Temp->Son1=Root;
			      Copy(Temp->Son2,Son1->Son1);
			      Root=Temp;
			    }
			    
			  else
			    {
			      // The tree is of the form ((?*A)*C)*?
			      
			      Father=Son1->Son1;
			      Temp=new Node;
			      Temp->Type=ArithmeticOperator1;
			      Temp->Son1=Root;
			      Copy(Temp->Son2,Root);
			      Father=this;
			      Root=Temp;
			    }
			    
			  Son1->Son2->Type=CreationOperator;
			  Son2->Type=AnnihilationOperator;
			  Root->Son1->NormalOrder(Root->Son1,Root->Son1);
			  Root->Son2->NormalOrder(Root->Son2,Root->Son2);
			}
			
		      else
			{
			  // The creation and annihilation operators commute.
			  
			  Temp=Son1->Son2;
			  Son1->Son2=Son2;
			  Son2=Temp;
			  Root->NormalOrder(Root,Root);
			}
		    }
		    
		  else if (Son1->Type==AnnihilationOperator && Son2->Type==CreationOperator)
		    {
		      if (Son1->Son1->Evaluate()==Son2->Son1->Evaluate())
			{
			  // The creation and annihilation operators do not commute.
			  
			  if (Root==this)
			    {
			      // The tree is of the form A*C
			      Temp=new Node;
			      Temp->Type=ArithmeticOperator1;
			      Temp->Son1=Root;
			      Temp->Son2=new Node;
			      Temp->Son2->Type=ComplexNumber;
			      Temp->Son2->Value=1.0;
			      Root=Temp;
			    }
			    
			  else
			    {
			      // The tree is of the form (A*C)*?
			      
			      Father=new Node;
			      Father->Type=ComplexNumber;
			      Father->Value=1.0;
			      Temp=new Node;
			      Temp->Type=ArithmeticOperator1;
			      Temp->Son1=Root;
			      Copy(Temp->Son2,Root);
			      delete Father;
			      Father=this;
			      Root=Temp;
			    }
			    
			  Son1->Type=CreationOperator;
			  Son2->Type=AnnihilationOperator;
			  Root->Son1->NormalOrder(Root->Son1,Root->Son1);
			  Root->Son2->NormalOrder(Root->Son2,Root->Son2);
			}
			
		      else
			{
			  // The creation and annihilation operators commute.
			  
			  Temp=Son1;
			  Son1=Son2;
			  Son2=Temp;
			  Root->NormalOrder(Root,Root);
			}
		    }
		    
		  else
		    {
		      Son1->NormalOrder(Root,Son1);
		      Son2->NormalOrder(Root,Son2);
		    }
		}
	    }
	    
	  void SpaceOrder(Node *Root)
	    {
	      Complex Index;

	      if (Type==ArithmeticOperator1)
		{
		  Son1->SpaceOrder(Son1);
		  Son2->SpaceOrder(Son2);
		}
		
	      else if (Type==ArithmeticOperator3 && Son1->Type==ArithmeticOperator3)
		{
		  if (Son2->Type==AnnihilationOperator && Son1->Son2->Type==AnnihilationOperator)
		    {
		      if (Son2->Son1->Value.Re()>Son1->Son2->Son1->Value.Re())
			{
			  Index=Son2->Son1->Value;
			  Son2->Son1->Value=Son1->Son2->Son1->Value;
			  Son1->Son2->Son1->Value=Index;
			  Root->SpaceOrder(Root);
			}
			
		      else
			Son1->SpaceOrder(Root);
		    }
		
		  else if (Son2->Type==CreationOperator && Son1->Son2->Type==CreationOperator)
		    {
		      if (Son2->Son1->Value.Re()<Son1->Son2->Son1->Value.Re())
			{
			  Index=Son2->Son1->Value;
			  Son2->Son1->Value=Son1->Son2->Son1->Value;
			  Son1->Son2->Son1->Value=Index;
			  Root->SpaceOrder(Root);
			}
			
		      else
			Son1->SpaceOrder(Root);
		    }
		
		  else
		    Son1->SpaceOrder(Root);
		}
		
	      else if (Type==ArithmeticOperator3 && Son1->Type==AnnihilationOperator && Son2->Type==AnnihilationOperator)
		{
		  if (Son2->Son1->Value.Re()>Son1->Son1->Value.Re())
		    {
		      Index=Son1->Son1->Value;
		      Son1->Son1->Value=Son2->Son1->Value;
		      Son2->Son1->Value=Index;
		      Root->SpaceOrder(Root);
		    }
		}
		
	      else if (Type==ArithmeticOperator3 && Son1->Type==CreationOperator && Son2->Type==CreationOperator)
		{
		  if (Son2->Son1->Value.Re()<Son1->Son1->Value.Re())
		    {
		      Index=Son1->Son1->Value;
		      Son1->Son1->Value=Son2->Son1->Value;
		      Son2->Son1->Value=Index;
		      Root->SpaceOrder(Root);
		    }
		}
	    }
	    
	  void SimplifyNumericalFactor(Node *Root)
	    {
	      Node *Temp;
	      
	      if (Type==ArithmeticOperator1)
		{
		  Son1->SimplifyNumericalFactor(Son1);
		  Son2->SimplifyNumericalFactor(Son2);
		}
		
	      else if (Type==ArithmeticOperator3 && (Son2->Type==CreationOperator || Son2->Type==AnnihilationOperator))
		Son1->SimplifyNumericalFactor(Root);
		
	      else if (Type!=CreationOperator && Type!=AnnihilationOperator)
		{
		  Value=Evaluate();
		  Type=ComplexNumber;
                  if (fabs(Value.Im())<0.0000000001) Value=Value.Re();
		  if (Son1) {delete Son1; Son1=NULL;}
		  if (Son2) {delete Son2; Son2=NULL;}
		}
		
	      else
		{
		  Copy(Temp,this);
		  Son2=Temp;
		  Son1->Value=1.0;
		  Type=ArithmeticOperator3;
		}
	    }

	  void SimplifyNullTerms(Node *Root)
	    {
	      static bool Test;
	      
	      if (Root==this)
		Test=false;

	      if (Type==ArithmeticOperator1)
		{
		  Node *Temp;
		  
		  if (Son1->Type!=ArithmeticOperator1 && Son1->Evaluate().Modulous()<0.0000000001)
		    {
		      Temp=Son2;
		      Type=Temp->Type;
		      Value=Temp->Value;
		      delete Son1;
		      Son1=Temp->Son1;
		      Son2=Temp->Son2;
		      Temp->Son1=NULL;
		      Temp->Son2=NULL;
		      delete Temp;
		      Test=true;
		    }
		  
		  else if (Son2->Type!=ArithmeticOperator1 && Son2->Evaluate().Modulous()<0.0000000001)
		    {
		      Temp=Son1;
		      Type=Temp->Type;
		      Value=Temp->Value;
		      delete Son2;
		      Son1=Temp->Son1;
		      Son2=Temp->Son2;
		      Temp->Son1=NULL;
		      Temp->Son2=NULL;
		      delete Temp;
		      Test=true;
		    }
		    
		  else
		    {
		      Son1->SimplifyNullTerms(Root);
		      Son2->SimplifyNullTerms(Root);
		    }
		}
		
	      if (Test && Root==this)
		Root->SimplifyNullTerms(Root);
	    }
	    
	  void ListTerms(Node *Root)
	    {
	      Complex Aux;
	      static int Diag,NonDiag;
	      static bool Problem;
	      
	      if (Root==this)
		{
		  Diag=0;
		  NonDiag=0;
		  Problem=false;
		}
		
	      switch (Type)
		{
		  case ArithmeticOperator1:
		    Son1->ListTerms(Root);

		    if (Son1->Type==ArithmeticOperator3 || Son1->Type==ComplexNumber)
		      {
			Aux=Son1->Evaluate();
			
			if (Son1->IsDiagonal(Son1))
			  {
			    Diag++;
			    cout << "\tDiagonal";
			    
			    if (Aux.Im()!=0.0)
			      {
				cout << " ---> Phase problem";
				Problem=true;
			      }
			  }
			  
			else
			  {
			    NonDiag++;
			    cout << "\tNon-diagonal";
			    
			    if (Aux.Im()!=0.0)
			      {
				cout << " ---> Phase problem";
				Problem=true;
			      }
			      
			    else if (Aux.Re()>0.0)
			      {
				cout << " ---> Sign problem";
				Problem=true;
			      }
			  }
			  
			cout << endl;
		      }
		      
		    Son2->ListTerms(Root);

		    if (Son2->Type==ArithmeticOperator3 || Son2->Type==ComplexNumber)
		      {
			Aux=Son2->Evaluate();
			
			if (Son2->IsDiagonal(Son2))
			  {
			    Diag++;
			    cout << "\tDiagonal";
			    
			    if (Aux.Im()!=0.0)
			      {
				cout << " ---> Phase problem";
				Problem=true;
			      }
			  }
			  
			else
			  {
			    NonDiag++;
			    cout << "\tNon-diagonal";
			    
			    if (Aux.Im()!=0.0)
			      {
				cout << " ---> Phase problem";
				Problem=true;
			      }
			      
			    else if (Aux.Re()>0.0)
			      {
				cout << " ---> Sign problem";
				Problem=true;
			      }
			  }
			  
			cout << endl;
		      }
		      
		    break;
		    
		  case ArithmeticOperator3:
		    Son1->ListTerms(Root);
		    Son2->ListTerms(Root);
		    break;
		    
		  case ComplexNumber:
		    Aux=this->Evaluate();
		    
		    if (Aux.Re() && Aux.Im())
		      cout << "(" << this->Evaluate() << ")";
		    else
		      cout << Aux;
		    
		    break;
		    
		  case CreationOperator:
		    cout << "C[" << (int) Son1->Value.Re() << "]";
		    break;
		    
		  case AnnihilationOperator:
		    cout << "A[" << (int) Son1->Value.Re() << "]";
		    break;
		    
		  default:
		    cout << "Unknow type of node is void ListTerms(Node *) function.\n";
		    cout << "Type=" << Type << endl;
		}
		
	      if (Root==this)
		{
		  cout << "\nNumber of diagonal terms: " << Diag << "\n";
		  cout << "Number of non-diagonal terms: " << NonDiag << "\n";
		  cout << "Total number of terms: " << Diag+NonDiag << "\n" << endl;
		  
		  if (Problem)
		    cout << "WARNING: Hamiltonian has a sign or a phase problem!\n" << endl;
		}
	    }
	    
          void ListQuantity(Node *Root)
            {
              Complex Aux;

              switch (Type)
                {
                  case ArithmeticOperator1:
                    Son1->ListQuantity(Root);

                    if (Son1->Type==ArithmeticOperator3 || Son1->Type==ComplexNumber)
                      {
                        Aux=Son1->Evaluate();
                        cout << endl;
                      }
                      
                    Son2->ListQuantity(Root);

                    if (Son2->Type==ArithmeticOperator3 || Son2->Type==ComplexNumber)
                      {
                        Aux=Son2->Evaluate();
                        cout << endl;
                      }
                      
                    break;
                    
                  case ArithmeticOperator3:
                    Son1->ListQuantity(Root);
                    Son2->ListQuantity(Root);
                    break;
                    
                  case ComplexNumber:
                    Aux=this->Evaluate();
                    
                    if (Aux.Re() && Aux.Im())
                      cout << "(" << this->Evaluate() << ")";
                    else
                      cout << Aux;
                    
                    break;
                    
                  case CreationOperator:
                    cout << "C[" << (int) Son1->Value.Re() << "]";
                    break;
                    
                  case AnnihilationOperator:
                    cout << "A[" << (int) Son1->Value.Re() << "]";
                    break;
                    
                  default:
                    cout << "Unknow type of node is void ListQuantity(Node *) function.\n";
                    cout << "Type=" << Type << endl;
                }
                
              if (Root==this)
                {
                }
            }
            
	  bool IsDiagonal(Node *Root)
	    {
	      static TheClassWithNoName Disc;
	      
	      if (Root==this)
		Disc.Clear();
	      
	      switch (Type)
		{
		  case ArithmeticOperator3:
		    Son1->IsDiagonal(Root);
		    Son2->IsDiagonal(Root);
		    break;
		    
		  case ComplexNumber:
		    break;
		    
		  case CreationOperator:
		    Disc.Add(Son1->Value.Re());
		    break;
		    
		  case AnnihilationOperator:
		    Disc.Sub(Son1->Value.Re());
		    break;
		    
		  default:
		    cout << "Unknow type of node is bool IsDiagonal(Node *) function.\n";
		    cout << "Type=" << Type << endl;
		}
		
	      if (Root==this)
		return (Disc.TheMemberFunctionWithNoName()==0);
	      else
		return false;
	    }
	    
	  bool SignOrPhaseProblem(Node *Root)
	    {
	      Complex Aux;
	      static bool Problem;
	      
	      if (Root==this)
		Problem=false;
	      
	      switch (Type)
		{
		  case ArithmeticOperator1:
		    Son1->SignOrPhaseProblem(Root);

		    if (Son1->Type==ArithmeticOperator3)
		      {
			Aux=Son1->Evaluate();
			
			if (Son1->IsDiagonal(Son1)) {if (Aux.Im()!=0.0) Problem=true;}
			  
			else
			  {
			    if (Aux.Im()!=0.0) Problem=true;
			    else if (Aux.Re()>0.0) Problem=true;
			  }
		      }
		      
		    Son2->SignOrPhaseProblem(Root);

		    if (Son2->Type==ArithmeticOperator3)
		      {
			Aux=Son2->Evaluate();
			
			if (Son2->IsDiagonal(Son2)) {if (Aux.Im()!=0.0) Problem=true;}
			  
			else
			  {
			    if (Aux.Im()!=0.0) Problem=true;
			    else if (Aux.Re()>0.0) Problem=true;
			  }
		      }
		      
		    break;
		    
		  case ArithmeticOperator3:
		    Son1->SignOrPhaseProblem(Root);
		    Son2->SignOrPhaseProblem(Root);
		    break;
		    
		  case ComplexNumber:
		  case CreationOperator:
		  case AnnihilationOperator:
		    break;
		    
		  default:
		    cout << "Unknow type of node is bool SignOrPhaseProblem(Node *) function.\n";
		    cout << "Type=" << Type << endl;
		}
		
	      if (Root==this)
		return Problem;
	      else
		return false;
	    }
	    
	  Complex Evaluate()
	    {
	      switch (Type)
		{
		  case ArithmeticOperator1:
		    return Son1->Evaluate()+Son2->Evaluate();
		  
		  case ArithmeticOperator2:
		    return Son1->Evaluate()-Son2->Evaluate();
		  
		  case ArithmeticOperator3:
		    return Son1->Evaluate()*Son2->Evaluate();
		  
		  case ArithmeticOperator4:
		    return Son1->Evaluate()/Son2->Evaluate();
		  
		  case ArithmeticOperator5:
		    return Son1->Evaluate()%Son2->Evaluate();
		  
		  case ArithmeticOperator6:
		    return Son1->Evaluate()^Son2->Evaluate();
		  
		  case ComplexNumber:
		    return Value;
		  
		  case Literal:
		    return Sym->Expression->Root->Evaluate();
		    
		  case Sin:
		    return sin(Son1->Evaluate());
		    
		  case Cos:
		    return cos(Son1->Evaluate());
		    
		  case Tan:
		    return tan(Son1->Evaluate());
		    
		  case Exp:
		    return exp(Son1->Evaluate());
		    
		  case Ln:
		    return log(Son1->Evaluate());
		    
		  case Sqrt:
		    return sqrt(Son1->Evaluate());
                    
                  case PosX:
                    return GetPosX((int) Son1->Evaluate().Re());
		    
                  case PosY:
                    return GetPosY((int) Son1->Evaluate().Re());
                    
                  case PosZ:
                    return GetPosZ((int) Son1->Evaluate().Re());
                    
		  case Fact:
		    return Factorial((int) Son1->Evaluate().Re());
		    
		  case CreationOperator:
		  case AnnihilationOperator:
		  case NumberOperator:
		    return 1.0;
		    
		  default:
		    cout << "Unknown type of node in Complex Evaluate() function.\n";
		    return 0;
		}
	    }
	  
	  void Display()
	    {
	      switch (Type)
		{
		  case ArithmeticOperator1:
		    cout << "(";
		    Son1->Display();
		    cout << "+";
		    Son2->Display();
		    cout << ")";
		    break;
		  
		  case ArithmeticOperator2:
		    cout << "(";
		    Son1->Display();
		    cout << "-";
		    Son2->Display();
		    cout << ")";
		    break;
		  
		  case ArithmeticOperator3:
		    Son1->Display();
		    cout << "*";
		    Son2->Display();
		    break;
		  
		  case ArithmeticOperator4:
		    Son1->Display();
		    cout << "/";
		    Son2->Display();
		    break;
		  
		  case ArithmeticOperator5:
		    Son1->Display();
		    cout << "%";
		    Son2->Display();
		    break;
		  
		  case ArithmeticOperator6:
		    Son1->Display();
		    cout << "^";
		    Son2->Display();
		    break;
		  
		  case ComplexNumber:
		    cout << Value;
		    break;
		  
		  case Literal:
		    Sym->Expression->Root->Display();
		    break;
		    
		  case Sin:
		    cout << "sin(";
		    Son1->Display();
		    cout << ")";
		    break;
		    
		  case Cos:
		    cout << "cos(";
		    Son1->Display();
		    cout << ")";
		    break;
		    
		  case Tan:
		    cout << "tan(";
		    Son1->Display();
		    cout << ")";
		    break;
		    
		  case Exp:
		    cout << "exp(";
		    Son1->Display();
		    cout << ")";
		    break;
		    
		  case Ln:
		    cout << "ln(";
		    Son1->Display();
		    cout << ")";
		    break;
		    
		  case Sqrt:
		    cout << "sqrt(";
		    Son1->Display();
		    cout << ")";
		    break;
		    
                  case PosX:
                    cout << "PosX(";
                    Son1->Display();
                    cout << ")";
                    break;
                    
                  case PosY:
                    cout << "PosY(";
                    Son1->Display();
                    cout << ")";
                    break;
                    
                  case PosZ:
                    cout << "PosZ(";
                    Son1->Display();
                    cout << ")";
                    break;
                    
		  case Fact:
		    cout << "Fact(";
		    Son1->Display();
		    cout << ")";
		    break;
		    
		  case CreationOperator:
		    cout << "C[";
		    Son1->Display();
		    cout << "]";
		    break;

		  case AnnihilationOperator:
		    cout << "A[";
		    Son1->Display();
		    cout << "]";
		    break;
		    
		  case NumberOperator:
		    cout << "N[";
		    Son1->Display();
		    cout << "]";
		    break;
		    
		  default:
		    cout << "Unknown type of node in void Display() function.\n";
		    cout << "Type=" << Type << endl;
		}
	    }
	  
          friend class MathExpression;
          friend class EasyMathExpression;
      };
  
    class Symbol
      {
        short Type;
        std::string Name;
	Complex Value;
        MathExpression *Expression;
        Symbol *Next;
        
        public:
          Symbol() {Expression=NULL; Next=NULL;}
          ~Symbol() {if (Expression) delete Expression;}
	  
	  static bool AddKeyword(const char *N)
	    {
	      return Add(N,Keyword,0);
	    }

	  static bool AddConstant(const char *N,const Complex &C)
	    {
	      return Add(N,Constant,C);
	    }

	  static bool AddVariable(const char *N)
	    {
	      return Add(N,Variable,0);
	    }

	  static bool SetVariable(const char *N,const Complex &C)
	    {
	      Symbol *Sym=Find(N);
	
	      if (Sym)
		{
		  Sym->Value=C;
		  return true;
		}
	  
	      return false;
	    }

	  static bool Add(const char *N,short T,const Complex &C)
	    {
	      if (Find(N))
		return false;
	
	      if (First)
		{
		  Last->Next=new Symbol;
		  Last=Last->Next;
		}
	    
	      else
		{
		  First=new Symbol;
		  Last=First;
		}
	    
	      Last->Next=NULL;
	      Last->Type=T;
	      Last->Name=std::string(N);
	
	      switch (T)
		{
		  case Constant:
		  case Variable:
		  case CreationOperator:
		  case AnnihilationOperator:
		  case NumberOperator:
		    Last->Value=C;
		    break;
		    
		  case Keyword:
		    break;

		  case Literal:
		    Last->Expression=new MathExpression;
		    Last->Expression->Root=new Node;
		    Last->Expression->Root->Type=ComplexNumber;
		    Last->Expression->Root->Value=C;
		    break;
		    
		  default:
		    cout << "Error in static bool Add(contst char *,short,Complex &) function:\n";
		    cout << "Type " << T << " is not supposed to occur in this function!" << endl;
		    return false;
		}
	  
	      return true;
	    }
      
	  static bool AddExpression(const char *N,MathExpression *H)
	    {
	      if (Find(N))
		return false;
	
	      if (First)
		{
		  Last->Next=new Symbol;
		  Last=Last->Next;
		}
	    
	      else
		{
		  First=new Symbol;
		  Last=First;
		}
	    
	      Last->Next=NULL;
	      Last->Type=Literal;
	      Last->Name=std::string(N);
	      Last->Expression=H;
	      return true;
	    }
      
	  static bool Delete(const char *N)
	    {
	      Symbol *Sym=Find(N);
	
	      if (Sym)
		{
		  if (Sym==First && Sym==Last)
		    {
		      First=NULL;
		      Last=NULL;
		      delete Sym;
		      return true;
		    }
	      
		  if (Sym==First)
		    {
		      First=First->Next;
		      delete Sym;
		      return true;
		    }
	      
		  if (Sym==Last)
		    {
		      Symbol *Head=First;
		
		      while (Head->Next!=Last)
			Head=Head->Next;
		
		      Last=Head;
		      Last->Next=NULL;
		      delete Sym;
		      return true;
		    }
	      
		  Symbol *Head=First;
	    
		  while (Head->Next!=Sym)
		    Head=Head->Next;
	    
		  Head->Next=Head->Next->Next;
		  delete Sym;
		  return true;
		}
	  
	      else
		return false;
	    }
            
	  friend class MathExpression;
	  friend class EasyMathExpression;
      };
      
    class Pair
      {
	int i,j;
	Pair *Next;
	
	public:
	  Pair() {Next=NULL;}
	  ~Pair() {if (Next) delete Next;}
	  
	  static void AddPair(int i,int j)
	    {
	      if (FirstPair)
		{
		  LastPair->Next=new Pair;
		  LastPair=LastPair->Next;
		}
	    
	      else
		{
		  FirstPair=new Pair;
		  LastPair=FirstPair;
		}
	    
	      LastPair->i=i;
	      LastPair->j=j;
	    }
      
	  friend class MathExpression;
      };

    class Species
      {
	char *Name;
	int Index,Population;
        unsigned int Max;
	Species *Next;
	
	public:
	  Species() {Name=NULL; Population=0; Next=NULL;}
	  ~Species() {if (Name) delete []Name; delete Next;}
	  friend class MathExpression;
      };
      
    static int NumSites,Lx,Ly,Lz,NumSpecies;
    static char Error[100],SimulName[100],Status[255];
    static Symbol *First,*Last;
    static Pair *FirstPair,*LastPair;
    static Species *FirstSpecies,*LastSpecies;
    static std::vector <string> MeasurableList;
    Node *Root;
    
    static void LocalizedError(const char *Message,Parser::TokenHandle Op)
      {
	int Row,Col;
	Row=Op->Row();
	Col=Op->Column();
	strcpy(Error,Message);
	strcpy(Error+strlen(Error),", row ");
	strcpy(Error+strlen(Error),Parser::IntToChars(Row));
	strcpy(Error+strlen(Error)," column ");
	strcpy(Error+strlen(Error),Parser::IntToChars(Col));
	strcpy(Error+strlen(Error),".");
      }
      
    static char *Expression(Node *&Father,Parser::TokenHandle Begin,Parser::TokenHandle End)
      {
	int Stack=0;
	Parser::TokenHandle Op=End;
	
	while (Op!=Begin && (Stack || !(Op->Compare('+') || Op->Compare('-'))))
	  {
	    if (Op->Compare('(') || Op->Compare('[') || Op->Compare('{') || Op->Compare('<'))
	      Stack++;
	      
	    if (Op->Compare(')') || Op->Compare(']') || Op->Compare('}') || Op->Compare('>'))
	      Stack--;
	      
	    Op=Parser::PrevToken(Op);
	  }
	
	if (Op==End && (Op->Compare('+') || Op->Compare('-')))
	  {
	    LocalizedError("Argument is expected after binary operator",Op);
	    return Error;
	  }
	  
	Parser::TokenHandle EndLeft=Parser::PrevToken(Op);
	Parser::TokenHandle BeginRight=Parser::NextToken(Op);
	  
	if (Op==Begin)
	  {
	    if (Op->Compare('+'))
	      return Term(Father,BeginRight,End);
	    
	    else if (Op->Compare('-'))
	      {
		Father=new Node;
		Father->Type=ArithmeticOperator2;
		Father->Son1=new Node;
		Father->Son1->Type=ComplexNumber;
		Father->Son1->Value=0.0;
		return Term(Father->Son2,BeginRight,End);
	      }
	    
	    else
	      return Term(Father,Begin,End);
	  }
	
	Father=new Node;
	  
	if (Op->Compare('+'))
	  Father->Type=ArithmeticOperator1;
	else
	  Father->Type=ArithmeticOperator2;
	    
	char *Test=Expression(Father->Son1,Begin,EndLeft);
	
	if (Test)
	  return Test;
	  
	return Term(Father->Son2,BeginRight,End);
      }
      
    static char *Term(Node *&Father,Parser::TokenHandle Begin,Parser::TokenHandle End)
      {
	int Stack=0;
	Parser::TokenHandle Op=End;
	
	while (Op!=Begin && (Stack || !(HiddenProduct(Op,End) || Op->Compare('*') || Op->Compare('/') || Op->Compare('%') || Op->Compare('^'))))
	  {
	    if (Op->Compare('(') || Op->Compare('[') || Op->Compare('{') || Op->Compare('<'))
	      Stack++;
	      
	    if (Op->Compare(')') || Op->Compare(']') || Op->Compare('}') || Op->Compare('>'))
	      Stack--;
	      
	    Op=Parser::PrevToken(Op);
	  }
	
	if (Op==End && (Op->Compare('*') || Op->Compare('/') || Op->Compare('%') || Op->Compare('^')))
	  {
	    LocalizedError("Argument is expected after binary operator",Op);
	    return Error;
	  }
	  
	if (Op==Begin && (Op->Compare('*') || Op->Compare('/') || Op->Compare('%') || Op->Compare('^')))
	  {
	    LocalizedError("Argument is expected before binary operator",Op);
	    return Error;
	  }
	  
	if (Op==Begin && !HiddenProduct(Op,End))
	  return Factor(Father,Begin,End);	
	  
	Parser::TokenHandle EndLeft=Parser::PrevToken(Op);
	Parser::TokenHandle BeginRight=Parser::NextToken(Op);	  
	Father=new Node;

	if (Op->Compare('*'))
	  Father->Type=ArithmeticOperator3;
	else if (Op->Compare('/'))
	  Father->Type=ArithmeticOperator4;
	else if (Op->Compare('%'))
	  Father->Type=ArithmeticOperator5;
	else if (Op->Compare('^'))
	  Father->Type=ArithmeticOperator6;
	
	else if (HiddenProduct(Op,End))
	  {
	    Father->Type=ArithmeticOperator3;
	    EndLeft=Op;
	  }
	  
	else
	  {
	    LocalizedError("Cannot read expression",Op);
	    return Error;
	  }
	  
	char *Test=Term(Father->Son1,Begin,EndLeft);
	
	if (Test)
	  return Test;
	  
	return Factor(Father->Son2,BeginRight,End);
      }
	
    static bool HiddenProduct(Parser::TokenHandle Left,Parser::TokenHandle End)
      {
	Parser::TokenHandle Right=Parser::NextToken(Left);
	
	if (!Right || Left==End)
	  return false;
	
	short Ltype=Left->Type();
	short Rtype=Right->Type();
	
	if (Rtype==Parser::Identifier)
	  {
	    if (Ltype==Parser::Identifier || Ltype==Parser::Number)
	      return true;
	    
	    if (Left->Compare(']') || Left->Compare(')'))
	      return true;
	  }
	  
	if (Ltype==Parser::Number && Right->Compare('('))
	  return true;
	
	if (Ltype==Parser::Identifier && Find((char *) Left->Value()) && Right->Compare('('))
	  return true;
	
	if ((Left->Compare(')') || Left->Compare(']')) && Rtype==Parser::Number)
	  return true;
	
	if (Left->Compare(']') && Right->Compare('('))
	  return true;
	  
	if (Left->Compare(')') && Right->Compare('('))
	  return true;
	  
	return false;
      }
      
    static char *Factor(Node *&Father,Parser::TokenHandle Begin,Parser::TokenHandle End)
      {
	if (Begin->Compare('(') && End->Compare(')'))
	  return Expression(Father,Parser::NextToken(Begin),Parser::PrevToken(End));
	
	if (Begin->Compare((char *) "Sum") || Begin->Compare((char *) "Prod"))
	  return SumOrProd(Father,Begin,End);
	
	switch (Begin->Type())
	  {
	    char *Name,*Err;
	    Symbol *Sym;
	      
	    case Parser::Number:
	      break;
	      
	    case Parser::Identifier:
	      Name=(char *) Begin->Value();
	      
	      if (!strcmp(Name,"sin"))
		{
		  Begin=Parser::NextToken(Begin);
		    
		  if (!Begin->Compare('(') || !End->Compare(')'))
		    {
		      LocalizedError("Argument of function must be enclosed into parentheses ()",Begin);
		      return Error;
		    }
		    
		  Father=new Node;
		  Father->Type=Sin;
		  return Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
		}
		
	      else if (!strcmp(Name,"cos"))
		{
		  Begin=Parser::NextToken(Begin);
		    
		  if (!Begin->Compare('(') || !End->Compare(')'))
		    {
		      LocalizedError("Argument of function must be enclosed into parentheses ()",Begin);
		      return Error;
		    }
		    
		  Father=new Node;
		  Father->Type=Cos;
		  return Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
		}
		
	      else if (!strcmp(Name,"tan"))
		{
		  Begin=Parser::NextToken(Begin);
		    
		  if (!Begin->Compare('(') || !End->Compare(')'))
		    {
		      LocalizedError("Argument of function must be enclosed into parentheses ()",Begin);
		      return Error;
		    }
		    
		  Father=new Node;
		  Father->Type=Tan;
		  return Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
		}
		
	      else if (!strcmp(Name,"exp"))
		{
		  Begin=Parser::NextToken(Begin);
		    
		  if (!Begin->Compare('(') || !End->Compare(')'))
		    {
		      LocalizedError("Argument of function must be enclosed into parentheses ()",Begin);
		      return Error;
		    }
		    
		  Father=new Node;
		  Father->Type=Exp;
		  return Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
		}
		
	      else if (!strcmp(Name,"ln"))
		{
		  Begin=Parser::NextToken(Begin);
		    
		  if (!Begin->Compare('(') || !End->Compare(')'))
		    {
		      LocalizedError("Argument of function must be enclosed into parentheses ()",Begin);
		      return Error;
		    }
		    
		  Father=new Node;
		  Father->Type=Ln;
		  return Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
		}
		
	      else if (!strcmp(Name,"sqrt"))
		{
		  Begin=Parser::NextToken(Begin);
		    
		  if (!Begin->Compare('(') || !End->Compare(')'))
		    {
		      LocalizedError("Argument of function must be enclosed into parentheses ()",Begin);
		      return Error;
		    }
		    
		  Father=new Node;
		  Father->Type=Sqrt;
		  return Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
		}
		
              else if (!strcmp(Name,"PosX"))
                {
                  Begin=Parser::NextToken(Begin);
                    
                  if (!Begin->Compare('(') || !End->Compare(')'))
                    {
                      LocalizedError("Argument of function must be enclosed into parentheses ()",Begin);
                      return Error;
                    }
                    
                  if (NumSites==0)
                    {
                      LocalizedError("PosX can be used only after lattice declaration",Begin);
                      return Error;
                    }
                    
                  Father=new Node;
                  Father->Type=PosX;
                  return Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
                }
                
              else if (!strcmp(Name,"PosY"))
                {
                  Begin=Parser::NextToken(Begin);
                    
                  if (!Begin->Compare('(') || !End->Compare(')'))
                    {
                      LocalizedError("Argument of function must be enclosed into parentheses ()",Begin);
                      return Error;
                    }
                    
                  if (NumSites==0)
                    {
                      LocalizedError("PosY can be used only after lattice declaration",Begin);
                      return Error;
                    }
                    
                  if (Ly==0)
                    {
                      LocalizedError("PosY can be used only with two-dimensional or three-dimensional lattices",Begin);
                      return Error;
                    }
                    
                  Father=new Node;
                  Father->Type=PosY;
                  return Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
                }
                
              else if (!strcmp(Name,"PosZ"))
                {
                  Begin=Parser::NextToken(Begin);
                    
                  if (!Begin->Compare('(') || !End->Compare(')'))
                    {
                      LocalizedError("Argument of function must be enclosed into parentheses ()",Begin);
                      return Error;
                    }
                    
                  if (NumSites==0)
                    {
                      LocalizedError("PosZ can be used only after lattice declaration",Begin);
                      return Error;
                    }
                    
                  if (Lz==0)
                    {
                      LocalizedError("PosZ can be used only with three-dimensional lattices",Begin);
                      return Error;
                    }
                    
                  Father=new Node;
                  Father->Type=PosZ;
                  return Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
                }
                
	      else if (!strcmp(Name,"Fact"))
		{
		  Begin=Parser::NextToken(Begin);
		    
		  if (!Begin->Compare('(') || !End->Compare(')'))
		    {
		      LocalizedError("Argument of function must be enclosed into parentheses ()",Begin);
		      return Error;
		    }
		    
		  Father=new Node;
		  Father->Type=Fact;
		  return Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
		}
		
	      else if (!(Sym=Find(Name)))
		{
		  LocalizedError("Undefined identifier",Begin);
		  return Error;
		}

	      switch (Sym->Type)
		{
		  case Keyword:
		    if (!strcmp(Name,"Random"))
		      {
			Father=new Node;
			Father->Type=ComplexNumber;
			Father->Value=RNG::Uniform();
			return NULL;
		      }
		      
		    LocalizedError("Keyword not allowed in mathematical expression",Begin);
		    return Error;

		  case Constant:
		  case Variable:
		    Father=new Node;
		    Father->Type=ComplexNumber;
		    Father->Value=Sym->Value;
		    return NULL;
		    
		  case Literal:
		    Father=new Node;
		    Father->Type=Literal;
		    Father->Sym=Sym;
		    return NULL;
		   
		  case CreationOperator:
		  case AnnihilationOperator:
		  case NumberOperator:
		    Begin=Parser::NextToken(Begin);
		    
		    if (!Begin->Compare('[') || !End->Compare(']'))
		      {
			LocalizedError("Argument of operator must be enclosed into brackets []",Begin);
			return Error;
		      }
		      
		    Father=new Node;
		    Father->Type=Sym->Type;
		    Father->Sym=Sym;
		    
		    if ((Sym->Value).Re()!=0)
		      {
			Father->Son1=new Node;
			Father->Son1->Type=ArithmeticOperator1;
			Father->Son1->Son1=new Node;
			Father->Son1->Son1->Type=ComplexNumber;
			Father->Son1->Son1->Value=Sym->Value*NumSites;
			Err=Expression(Father->Son1->Son2,Parser::NextToken(Begin),Parser::PrevToken(End));
			
			if (Err)
			  return Err;
			
			if (Father->Son1->Son2->Evaluate().Re()>=NumSites)
			  {
			    LocalizedError("Argument of operator is out of range",Begin);
			    return Error;
			  }
		      }
		      
		    else
		      {
			Err=Expression(Father->Son1,Parser::NextToken(Begin),Parser::PrevToken(End));
			
			if (Err)
			  return Err;
			
			if (Father->Son1->Evaluate().Re()>=NumSites)
			  {
			    LocalizedError("Argument of operator is out of range",Begin);
			    return Error;
			  }
		      }
		      
		    return NULL;
		    
		  default:
		    LocalizedError("Invalid symbol type",Begin);
		    return Error;
		}
  
	    default:
	      LocalizedError("Invalid character",Begin);
	      return Error;
	  }

	Father=new Node;
	Father->Type=ComplexNumber;
	Father->Value=*((double *) Begin->Value());
	return NULL;
      }
    
    static char *SumOrProd(Node *&Father,Parser::TokenHandle Begin,Parser::TokenHandle End)
      {
	short Type;
	char *Name1,*Name2,*Err;
	int Min,Max;
	Symbol *Sym;
	Parser::TokenHandle Handle;
	
	if (Begin->Compare((char *) "Sum"))
	  Type=ArithmeticOperator1;
	else
	  Type=ArithmeticOperator3;
	
	if (!CheckNextToken(Begin,'{'))
	  {
	    LocalizedError("An opening brace is expected",Begin);
	    return Error;
	  }

	Begin=Parser::NextToken(Begin);

	if (!CheckNextToken(Begin,'<'))
	  {
	    if (!CheckNextToken(Begin,(short) Parser::Identifier))
	      {
		LocalizedError("A name of index is expected",Begin);
		return Error;
	      }

	    Begin=Parser::NextToken(Begin);
	    Name1=(char *) Begin->Value();
	
	    if ((Sym=Find(Name1)))
	      {
		LocalizedError("Name of index has a previous declaration",Begin);
		return Error;
	      }
	      
	    if (!CheckNextToken(Begin,'}'))
	      {
		// The sum or product is of the form Sum{i,a,b}(E)

		if (!CheckNextToken(Begin,','))
		  {
		    LocalizedError("A coma is expected",Begin);
		    return Error;
		  }

		Begin=Parser::NextToken(Begin);
		Begin=Parser::NextToken(Begin);
		Handle=Parser::GetExpressionEnd(Begin);
		
		if ((Err=MathExpression::Define("######",Begin,Handle)))
		  return Err;
		
		Min=(int) GetValue("######").Re();
		Undefine("######");

		if (!CheckNextToken(Handle,','))
		  {
		    LocalizedError("A coma is expected",Begin);
		    return Error;
		  }

		Begin=Parser::NextToken(Handle);
		Begin=Parser::NextToken(Begin);
		Handle=Parser::GetExpressionEnd(Begin);
		
		if ((Err=MathExpression::Define("######",Begin,Handle)))
		  return Err;
		
		Max=(int) GetValue("######").Re();
		Undefine("######");

		if (!CheckNextToken(Handle,'}'))
		  {
		    LocalizedError("A closing brace is expected",Begin);
		    return Error;
		  }

		Begin=Parser::NextToken(Handle);
		Begin=Parser::NextToken(Begin);
		Symbol::AddVariable(Name1);
	
		if ((Err=SumOrProd(Father,Begin,End,Type,Name1,Min,Max)))
		  return Err;

		Undefine(Name1);
	      }
	      
	    else
	      {
		// The sum or product is of the form Sum{i}(E)

		if (NumSites<=0)
		  {
		    LocalizedError("Sum or product over all lattice sites can be used only after lattice declaration",Begin);
		    return Error;
		  }

		Begin=Parser::NextToken(Begin);
		Begin=Parser::NextToken(Begin);
		Min=0;
		Max=NumSites-1;
		Symbol::AddVariable(Name1);
	
		if ((Err=SumOrProd(Father,Begin,End,Type,Name1,Min,Max)))
		  return Err;

		Undefine(Name1);
	      }
	  }
	  
	else
	  {
	    // The sum or product is of the form Sum{<i,j>}(E)

	    Begin=Parser::NextToken(Begin);
	    
	    if (!CheckNextToken(Begin,(short) Parser::Identifier))
	      {
		LocalizedError("A name of index is expected",Begin);
		return Error;
	      }

	    Begin=Parser::NextToken(Begin);
	    Name1=(char *) Begin->Value();
	
	    if ((Sym=Find(Name1)))
	      {
		LocalizedError("Name of index has a previous declaration",Begin);
		return Error;
	      }

	    if (!CheckNextToken(Begin,','))
	      {
		LocalizedError("A coma is expected",Begin);
		return Error;
	      }

	    Begin=Parser::NextToken(Begin);

	    if (!CheckNextToken(Begin,(short) Parser::Identifier))
	      {
		LocalizedError("A name of index is expected",Begin);
		return Error;
	      }

	    Begin=Parser::NextToken(Begin);
	    Name2=(char *) Begin->Value();
	
	    if ((Sym=Find(Name2)))
	      {
		LocalizedError("Name of index has a previous declaration",Begin);
		return Error;
	      }

	    if (!CheckNextToken(Begin,'>'))
	      {
		LocalizedError("A closing angle  is expected",Begin);
		return Error;
	      }

	    Begin=Parser::NextToken(Begin);

	    if (!CheckNextToken(Begin,'}'))
	      {
		LocalizedError("A closing brace  is expected",Begin);
		return Error;
	      }
	      
	    if (!strcmp(Name1,Name2))
	      {
		LocalizedError("Indices in sum or product over first nearest neighbors must be different",Begin);
		return Error;
	      }

	    if (NumSites<=0)
	      {
		LocalizedError("Sum or product over first nearest neighbors can be used only after lattice declaration",Begin);
		return Error;
	      }

	    Begin=Parser::NextToken(Begin);
	    Begin=Parser::NextToken(Begin);
	    Symbol::AddVariable(Name1);
	    Symbol::AddVariable(Name2);
	    Pair *PairHandle=FirstPair;
	
	    if ((Err=SumOrProd(Father,Begin,End,Type,Name1,Name2,PairHandle)))
	      return Err;
	    
	    Undefine(Name1);
	    Undefine(Name2);
	  }
	  
	return NULL;
      }
      
    static char *SumOrProd(Node *&Father,Parser::TokenHandle Begin,Parser::TokenHandle End,short Type,const char *Name,int Min,int Max)
      {
	char *Test;
	Symbol::SetVariable(Name,Min);
	
	if (!Begin->Compare('(') || !End->Compare(')'))
	  {
	    LocalizedError("Expression in Sum or Prod declaration must be enclosed into parentheses ()",Begin);
	    return Error;
	  }
	
	if (Min==Max)
	  return Expression(Father,Begin,End);
	  
	else
	  {
	    Father=new Node;
	    Father->Type=Type;
	    
	    if ((Test=Expression(Father->Son1,Begin,End)))
	      return Test;
	  
	    return SumOrProd(Father->Son2,Begin,End,Type,Name,Min+1,Max);
	  }
      }
      
    static char *SumOrProd(Node *&Father,Parser::TokenHandle Begin,Parser::TokenHandle End,short Type,const char *Name1,const char *Name2,Pair *PairHandle)
      {
	char *Test;
	Symbol::SetVariable(Name1,PairHandle->i);
	Symbol::SetVariable(Name2,PairHandle->j);
	
	if (!Begin->Compare('(') || !End->Compare(')'))
	  {
	    LocalizedError("Expression in Sum or Prod declaration must be enclosed into parentheses ()",Begin);
	    return Error;
	  }
	
	if (!PairHandle->Next)
	  return Expression(Father,Begin,End);
	  
	else
	  {
	    Father=new Node;
	    Father->Type=Type;
	    
	    if ((Test=Expression(Father->Son1,Begin,End)))
	      return Test;
	  
	    return SumOrProd(Father->Son2,Begin,End,Type,Name1,Name2,PairHandle->Next);
	  }
      }
      
    static bool CheckNextToken(Parser::TokenHandle Handle,char C)
      {
	Handle=Parser::NextToken(Handle);
	
	if (!Handle)
	  return false;
	
	if (!Handle->Compare(C))
	  return false;
	  
	return true;
      }
	  
    static bool CheckNextToken(Parser::TokenHandle Handle,short Type)
      {
	Handle=Parser::NextToken(Handle);
	
	if (!Handle)
	  return false;
	
	if (Handle->Type()!=Type)
	  return false;
	  
	return true;
      }
	  
    public:
      typedef MathExpression * ExpressionHandle;
      
      MathExpression() {Root=NULL;}
      ~MathExpression() {if (Root) delete Root;}
      
      static Symbol *Find(const char *Name)
	{
	  Symbol *Head=First;
	  
	  while (Head && std::string(Name)!=Head->Name)
	    Head=Head->Next;
	  
	  return Head;
	}
    
      char *BuildExpression(Parser::TokenHandle Begin,Parser::TokenHandle End)
	{
	  if (Begin && End)
	    return Expression(Root,Begin,End);
	  else
	    return NULL;
	}
	
      static char *Define(const char *Name,Parser::TokenHandle Begin,Parser::TokenHandle End)
	{
	  char *Err;
	  MathExpression *Handle=new MathExpression;
	  
	  if ((Err=Handle->BuildExpression(Begin,End)))
	    return Err;

	  if (Symbol::AddExpression(Name,Handle))
	    return NULL;
	  else
	    return (char *) "Identifier has a previous declaration";
	}
	
      static char *Undefine(const char *Name)
	{
	  if (!Symbol::Delete(Name))
	    return (char *) "Cannot undefine a name that does not exist!\n";
	  else
	    return NULL;
	}
	
      static void Initialize()
	{
	  Complex I(0,1);
	  RNG::Initialize(24379);
	  Symbol::AddKeyword("Variable");
	  Symbol::AddKeyword("Operator");
	  Symbol::AddKeyword("Lattice");
	  Symbol::AddKeyword("Species");
	  Symbol::AddKeyword("Population");
	  Symbol::AddKeyword("Hamiltonian");
	  Symbol::AddKeyword("InverseTemperature");
	  Symbol::AddKeyword("Seed");
	  Symbol::AddKeyword("WarmTime");
	  Symbol::AddKeyword("MeasTime");
	  Symbol::AddKeyword("Bins");
	  Symbol::AddKeyword("Ensemble");
	  Symbol::AddKeyword("SimulName");
          Symbol::AddKeyword("Measure");
	  Symbol::AddKeyword("Random");
	  Symbol::AddConstant("Second",1);
	  Symbol::AddConstant("Minute",60);
	  Symbol::AddConstant("Hour",3600);
	  Symbol::AddConstant("Day",86400);
	  Symbol::AddConstant("Week",604800);
	  Symbol::AddConstant("Canonical",0);
	  Symbol::AddConstant("GrandCanonical",1);
	  Symbol::AddConstant("Pi",3.1415926535897932384626433832795);
	  Symbol::AddConstant("SquareRootOfMinusOne",I);
	  Symbol::AddConstant("Infinity",numeric_limits<double>::infinity());
	}
	
      static Complex GetValue(const char *Name)
	{
	  Symbol *Sym;
	  
	  if (!(Sym=Find(Name)))
	    {
	      cout << "Symbol '" << Name << "' cannot be found in symbol table!\n";
	      return 0;
	    }
	    
	  return Sym->Expression->Root->Evaluate();
	}
	
      static char *SetValue(const char *Name,Parser::TokenHandle Begin,Parser::TokenHandle End)
	{
	  char *Err;
	  Symbol *Sym;
	  
	  if (!(Sym=Find(Name)))
	    return (char *) "Identifier has not been declared!";
	  
	  if (Sym->Type==Constant)
	    return (char *) "Cannot overwrite a constant!";

	  Node *New;

	  if ((Err=Expression(New,Begin,End)))
	    return Err;

	  New->Substitute();
	  delete Sym->Expression->Root;
	  Sym->Expression->Root=New;
	  return NULL;
	}
	
      static void Display(const char *Name)
	{
	  Symbol *Sym;
	  
	  if (!(Sym=Find(Name)))
	    {
	      cout << "Symbol '" << Name << "' cannot be found in symbol table!\n";
	      return;
	    }
	    
	  Sym->Expression->Root->Display();
	  cout << endl;
	}
	
      static char *Lattice(Parser::TokenHandle Begin,Parser::TokenHandle End)
	{
	  char *Err;

	  if (Lx>0)
	    {
	      if (Ly>0)
		{
		  if ((Err=MathExpression::Define("######",Begin,End)))
		    return Err;
		
		  Lz=(int) GetValue("######").Re();
		  Undefine("######");
	  
		  if (Lz<=0)
		    {
		      LocalizedError("Arguments in lattice declaration must be a positive integers",Begin);
		      return Error;
		    }

		  NumSites*=Lz;
		  delete FirstPair;
		  FirstPair=NULL;
		  
		  for (int x=0;x<Lx;x++)
		    for (int y=0;y<Ly;y++)
		      for (int z=0;z<Lz;z++)
			{
			  Pair::AddPair(z*Lx*Ly+y*Lx+x,z*Lx*Ly+y*Lx+(x+1)%Lx);
			  Pair::AddPair(z*Lx*Ly+y*Lx+x,z*Lx*Ly+((y+1)%Ly)*Lx+x);
			  Pair::AddPair(z*Lx*Ly+y*Lx+x,((z+1)%Lz)*Lx*Ly+y*Lx+x);
			}
		}
		
	      else
		{
		  if ((Err=MathExpression::Define("######",Begin,End)))
		    return Err;
		
		  Ly=(int) GetValue("######").Re();
		  Undefine("######");
	  
		  if (Ly<=0)
		    {
		      LocalizedError("Arguments in lattice declaration must be a positive integers",Begin);
		      return Error;
		    }

		  NumSites*=Ly;
		  delete FirstPair;
		  FirstPair=NULL;
		  
		  for (int x=0;x<Lx;x++)
		    for (int y=0;y<Ly;y++)
		      {
			Pair::AddPair(y*Lx+x,y*Lx+(x+1)%Lx);
			Pair::AddPair(y*Lx+x,((y+1)%Ly)*Lx+x);
		      }
		}
	    }
	    
	  else
	    {
	      if ((Err=MathExpression::Define("######",Begin,End)))
		return Err;
		
	      Lx=(int) GetValue("######").Re();
	      Undefine("######");
	  
	      if (Lx<=0)
		{
		  LocalizedError("Arguments in lattice declaration must be a positive integers",Begin);
		  return Error;
		}

	      NumSites=Lx;

	      for (int i=0;i<Lx;i++)
		Pair::AddPair(i,(i+1)%Lx);
	    }
	  
	  return NULL;
	}
	
      static char *AddSpecies(const char *Name,const char *Name1,const char *Name2,const char *Name3,Parser::TokenHandle Begin,Parser::TokenHandle End)
	{
	  char *Err;
	  
	  if (FirstSpecies)
	    {
	      LastSpecies->Next=new Species;
	      LastSpecies=LastSpecies->Next;
	    }
	    
	  else
	    {
	      FirstSpecies=new Species;
	      LastSpecies=FirstSpecies;
	    }
	    
	  LastSpecies->Next=NULL;
	  LastSpecies->Name=new char[strlen(Name)+1];
	  strcpy(LastSpecies->Name,Name);
	  LastSpecies->Index=NumSpecies;
	  NumSpecies++;

	  if ((Err=MathExpression::Define("######",Begin,End)))
	    return Err;
		
          double Temp=GetValue("######").Re();
          
          if (Temp==numeric_limits<double>::infinity())
	    LastSpecies->Max=0;
          else
            LastSpecies->Max=(unsigned int) Temp;
          
	  Undefine("######");
	  
	  if (!Symbol::Add(Name1,CreationOperator,NumSpecies-1))
	    return (char *) "Name of creation operator has a previous declaration.";
	  
	  if (!Symbol::Add(Name2,AnnihilationOperator,NumSpecies-1))
	    return (char *) "Name of Annihilation operator has a previous declaration.";
	  
	  if (!Symbol::Add(Name3,NumberOperator,NumSpecies-1))
	    return (char *) "Name of number operator has a previous declaration.";

	  return NULL;
	}
	
      static char *SetPopulation(const char *Name,Parser::TokenHandle Begin,Parser::TokenHandle End)
	{
	  char *Err;
	  Species *Head=FirstSpecies;
	
	  while (Head && strcmp(Head->Name,Name))
	    Head=Head->Next;
	
	  if (Head)
	    {
	      if ((Err=MathExpression::Define("######",Begin,End)))
		return Err;
		
	      Head->Population=(int) GetValue("######").Re();
	      Undefine("######");
	      return NULL;
	    }
	  
	  else
	    {
	      LocalizedError("Name of species does not exist",Parser::PrevToken(Parser::PrevToken(Begin)));
	      return Error;
	    }
	}
    
      static char *ListHamiltonianTerms()
	{
	  Symbol *Sym=Find("#Hamiltonian");
	  
	  if (!Sym)
	    return (char *) "Hamiltonian cannot be found in symbol table!\n";
	  
	  Sym->Expression->Root->ListTerms(Sym->Expression->Root);
	  return NULL;
	}

      static string DisplayQuantity(const char *Name)
	{
	  Symbol *Sym=Find(Name);

          if (!Sym)
            return string("The quantity ")+string(Name)+string(" has not been defined!");
	  
          BuildMeasurable(Name);                                // Name of function is not consistent
          cout << "Displaying quantity " << Name << "\n";
	  Sym->Expression->Root->ListQuantity(Sym->Expression->Root);
          cout << endl;
	  return "";
	}

      static char *SignOrPhaseProblem()
	{
	  Symbol *Sym=Find("#Hamiltonian");
	  
	  if (!Sym)
	    return (char *) "Hamiltonian cannot be found in symbol table!\n";

	  if (Sym->Expression->Root->SignOrPhaseProblem(Sym->Expression->Root))
	    return (char *) "\nWARNING: Hamiltonian has a sign or a phase problem!\nUse option -Hamiltonian to see which terms cause the problem.\nSimulation aborded.\n";
	  
	  return NULL;
	}
	
      static char *BuildHamiltonian()
	{
	  Symbol *Sym=Find("#Hamiltonian");
	  
	  if (!Sym)
	    return (char *) "Hamiltonian cannot be found in symbol table!\n";

	  SetStatus("Substituting expressions in Hamiltonian");
	  Sym->Expression->Root->Substitute();
	  SetStatus("Expanding Hamiltonian");
	  Sym->Expression->Root->Expand(Sym->Expression->Root);
	  SetStatus("Expanding number operators");
	  Sym->Expression->Root->ExpandNumberOperator();
	  SetStatus("Simplifying argument of operators");
	  Sym->Expression->Root->SimplifyOperatorArgument();
	  SetStatus("Writing products from left to right");
	  Sym->Expression->Root->LeftProduct(Sym->Expression->Root);
	  SetStatus("Moving numbers to the left of operators");
	  Sym->Expression->Root->LeftNumbers(Sym->Expression->Root);
	  SetStatus("Writing operators in normal order");
	  Sym->Expression->Root->NormalOrder(Sym->Expression->Root,Sym->Expression->Root);
	  SetStatus("Ordering operators indices");
	  Sym->Expression->Root->SpaceOrder(Sym->Expression->Root);
	  SetStatus("Simplifying numerical factors");
	  Sym->Expression->Root->SimplifyNumericalFactor(Sym->Expression->Root);
	  SetStatus("Factorizing identical operators");
	  Sym->Expression->Root->FactorizeOperators(Sym->Expression->Root);
	  SetStatus("Removing null terms");
	  Sym->Expression->Root->SimplifyNullTerms(Sym->Expression->Root);
	  remove(Status);
	  return NULL;
	}
	
      static char *BuildMeasurable(const char *Name)
	{
	  Symbol *Sym=Find(Name);
	  
	  if (!Sym)
	    return (char *) "Name of measurable cannot be found in symbol table!\n";

          SetStatus("Substituting expressions in measurable");
	  Sym->Expression->Root->Substitute();
          SetStatus("Expanding measurable");
	  Sym->Expression->Root->Expand(Sym->Expression->Root);
          SetStatus("Expanding number operators");
	  Sym->Expression->Root->ExpandNumberOperator();
          SetStatus("Simplifying argument of operators");
	  Sym->Expression->Root->SimplifyOperatorArgument();
          SetStatus("Writing products from left to right");
	  Sym->Expression->Root->LeftProduct(Sym->Expression->Root);
          SetStatus("Moving numbers to the left of operators");
	  Sym->Expression->Root->LeftNumbers(Sym->Expression->Root);
          SetStatus("Ordering operators indices");
	  Sym->Expression->Root->SpaceOrder(Sym->Expression->Root);
          SetStatus("Simplifying numerical factors");
	  Sym->Expression->Root->SimplifyNumericalFactor(Sym->Expression->Root);
          SetStatus("Factorizing identical operators");
	  Sym->Expression->Root->FactorizeOperators(Sym->Expression->Root);
          SetStatus("Removing null terms");
	  Sym->Expression->Root->SimplifyNullTerms(Sym->Expression->Root);
          remove(Status);
	  return NULL;
	}
	
      static void AddMeasurable(const char *Name)
	{
	  MeasurableList.push_back(std::string(Name));
	}
	
      static void SetStatus(const char *Name)
	{
	  remove(Status);
	  fstream File;
	  strcpy(Status,SimulName);
	  strcpy(Status+strlen(Status)," - ");
	  strcpy(Status+strlen(Status),Name);
	  File.open(Status,ios::out);
	  File.close();
	}
	
      static void SetSimulName(const char *Name)
	{
	  strcpy(SimulName,Name);
	}
	
      static char *GetSimulName()
	{
	  return SimulName;
	}
	
      static int GetNumSites()
	{
	  return NumSites;
	}
	
      static int GetNumSpecies()
	{
	  return NumSpecies;
	}
	
      static int GetNumIndices()
	{
	  return NumSpecies*NumSites;
	}
	
      static int GetNmax(int GlobalInd) 
	{
	  int S=GlobalInd/NumSites;

	  Species *Head=FirstSpecies;
	  
	  for (int i=0;i<S;i++)
	    Head=Head->Next;
	  
	  return Head->Max;
	}    
	
      static int GetPopulation(int S)
	{
	  Species *Head=FirstSpecies;
	  
	  for (int i=0;i<S;i++)
	    Head=Head->Next;
	  
	  return Head->Population;
	}
	
      static bool IsKeyword(const char *Name)
        {
          Symbol *Sym=Find(Name);
            
          if (Sym && Sym->Type==Keyword)
            return true;
          else
            return false;
        }
      
      static std::vector<string> &GetMeasurableList()
	{
	  return MeasurableList;
	}
        
      static double GetPosX(int i)
        {
          return i%Lx-(Lx-1)/2.0;
        }
	
      static double GetPosY(int i)
        {
          return (i/Lx)%Ly-(Ly-1)/2.0;
        }
        
      static double GetPosZ(int i)
        {
          return i/(Lx*Ly)-(Lz-1)/2.0;
        }
        
      friend class EasyMathExpression;
  };

int MathExpression::NumSites=0;
int MathExpression::Lx=0;
int MathExpression::Ly=0;
int MathExpression::Lz=0;
int MathExpression::NumSpecies=0;
char MathExpression::Error[100];
char MathExpression::SimulName[100]="";
char MathExpression::Status[255]="";
MathExpression::Symbol *MathExpression::First=NULL;
MathExpression::Symbol *MathExpression::Last=NULL;
MathExpression::Pair *MathExpression::FirstPair=NULL;
MathExpression::Pair *MathExpression::LastPair=NULL;
MathExpression::Species *MathExpression::FirstSpecies=NULL;
MathExpression::Species *MathExpression::LastSpecies=NULL;
std::vector<string> MathExpression::MeasurableList;

#endif

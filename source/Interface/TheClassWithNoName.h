// Why this class has no name?
// Because I have been unable so far to find a name that is relavant to what it does.

#ifndef TheClassWithNoNameDef
#define TheClassWithNoNameDef

class TheClassWithNoName
  {
    class Node
      {
	int Occurence;
	double Label;
	Node *Son1,*Son2;
	
	public:
	  Node(double L,int S) {Occurence=S; Label=L; Son1=NULL; Son2=NULL;}
	  ~Node(void) {if (Son1) delete Son1; if (Son2) delete Son2;}
	  
	  int Add(double L)
	    {
	      if (L==Label)
		{
		  Occurence++;
		  return (Occurence>0 ? 1: -1);
		}
	      
	      else if (L<Label)
		{
		  if (Son1)
		    return Son1->Add(L);
		  
		  else
		    {
		      Son1=new Node(L,1);
		      return 1;
		    }
		}
		
	      else
		{
		  if (Son2)
		    return Son2->Add(L);
		  
		  else
		    {
		      Son2=new Node(L,1);
		      return 1;
		    }
		}
	    }
	    
	  int Sub(double L)
	    {
	      if (L==Label)
		{
		  Occurence--;
		  return (Occurence<0 ? 1: -1);
		}
	      
	      else if (L<Label)
		{
		  if (Son1)
		    return Son1->Sub(L);
		  
		  else
		    {
		      Son1=new Node(L,-1);
		      return 1;
		    }
		}
		
	      else
		{
		  if (Son2)
		    return Son2->Sub(L);
		  
		  else
		    {
		      Son2=new Node(L,-1);
		      return 1;
		    }
		}
	    }
	    
	  friend class TheClassWithNoName;
      };
      
    int TheVariableWithNoName;
    Node *Root;
    
    public:
      TheClassWithNoName(void) {TheVariableWithNoName=0; Root=NULL;}
      ~TheClassWithNoName(void) {if (Root) delete Root;}
      
      void Add(double L)
	{
	  if (Root)
	    TheVariableWithNoName+=Root->Add(L);
	  
	  else
	    {
	      Root=new Node(L,1);
	      TheVariableWithNoName=1;
	    }
	}
      
      void Sub(double L)
	{
	  if (Root)
	    TheVariableWithNoName+=Root->Sub(L);
	  
	  else
	    {
	      Root=new Node(L,-1);
	      TheVariableWithNoName=-1;
	    }
	}
	
      void Clear(void)
	{
	  if (Root)
	    {
	      delete Root;
	      Root=NULL;
	      TheVariableWithNoName=0;
	    }
	}
	
      int TheMemberFunctionWithNoName(void) {return TheVariableWithNoName;}
  };

#endif

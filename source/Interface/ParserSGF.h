// ************************************************************************************
// * This class transforms the input tokens from class 'Parser' into a linked list of *
// * SGF-Script commands.                                                             *
// * ________________________________________________________________________________ *
// * Dr. Valy G. Rousseau - December 2nd, 2010                                        *
// ************************************************************************************
//
// A command is accessible via a handle:
//   ParserSGF::CommandHandle MyHandle;
//
// Commands can be of the following types:
//   - ParserSGF::VariableDeclaration
//   - ParserSGF::VariableArrayDeclaration
//   - ParserSGF::VariableDeclarationAffectation
//   - ParserSGF::OperatorDeclaration
//   - ParserSGF::OperatorArrayDeclaration
//   - ParserSGF::OperatorDeclarationAffectation
//   - ParserSGF::Affectation
//   - ParserSGF::SpeciesDeclaration
//   - ParserSGF::LatticeDeclaration
//   - ParserSGF::SeedDeclaration
//   - ParserSGF::HamiltonianDeclaration
//   - ParserSGF::TemperatureDeclaration
//   - ParserSGF::WarmTimeDeclaration
//   - ParserSGF::MeasTimeDeclaration
//   - ParserSGF::BinsDeclaration
//   - ParserSGF::EnsembleDeclaration
//   - ParserSGF::BoundariesDeclaration
//   - ParserSGF::SimulNameDeclaration
//   - ParserSGF::PopulationDeclaration
//   - ParserSGF::MeasureDeclaration
//
// THE METHODS
//
//   char *ReadTokens(Parser *Input) 
//   Builds a list of commands from Input. The returned value is a pointer
//   onto an error message, if any, otherwise it is NULL.
//
//   char *ExecuteCommands(void) 
//   Executes the list of SGF commands. The returned value is a pointer
//   onto an error message, if any, otherwise it is NULL.
//
//   CommandHandle GetCommandHandle(void)
//   Returns a handle for the current command.
//   A NULL handle indicates that no commands are available.
//
//   static CommandHandle NextCommand(CommandHandle Handle)
//   Returns a handle for command after Handle.
//   The returned handle is NULL if no next command is available.
//
//   static short GetCommandType(CommandHandle Handle)
//   Returns the type of command 'Handle'.
//   The returned handle is -1 if Handle is NULL.
//
//   static char *GetIdentifier(CommandHandle Handle)
//   Returns the identifier that appears in command Handle.
//
//   void Clear(void)
//   Clears all commands.

#include <string>

#ifndef ParserSGFDef
#define ParserSGFDef

class ParserSGF
  {
    class Command
      {
        short Type;
	    int Row,Column;
        std::string String,String1,String2,String3;
	    Parser::TokenHandle Begin,End;
        Command *Next,*Prev;
        
        public:
          const char * CString() const {return String.c_str();} 
          const char * CString1() const {return String1.c_str();} 
          const char * CString2() const {return String2.c_str();} 
          const char * CString3() const {return String3.c_str();} 

          Command(void) {}
          ~Command(void) {}
	  short GetTokenType(void) {return Type;}
	  friend class ParserSGF;
      };
      
    Command *First,*Last,*Head;
    char Error[256];
    
    void Add(short Type,char *Id,int R,int C)
      {
        if (First)
          {
            Last->Next=new Command;
            Last->Next->Prev=Last;
            Last=Last->Next;
          }
            
        else
          {
            First=new Command;
            First->Prev=NULL;
            Last=First;
            Head=First;
          }

        Last->Next=NULL;
        Last->Type=Type;
        Last->Row=R;
        Last->Column=C;
        Last->String=std::string(Id);
      }
      
    void Add(short Type,Parser::TokenHandle B,Parser::TokenHandle E,int R,int C)
      {
        if (First)
          {
            Last->Next=new Command;
            Last->Next->Prev=Last;
            Last=Last->Next;
          }
            
        else
          {
            First=new Command;
            First->Prev=NULL;
            Last=First;
            Head=First;
          }

        Last->Next=NULL;
        Last->Type=Type;
        Last->Row=R;
        Last->Column=C;
	Last->Begin=B;
	Last->End=E;
      }
    
    void Add(short Type,char *Id,Parser::TokenHandle B,Parser::TokenHandle E,int R,int C)
      {
        if (First)
          {
            Last->Next=new Command;
            Last->Next->Prev=Last;
            Last=Last->Next;
          }
            
        else
          {
            First=new Command;
            First->Prev=NULL;
            Last=First;
            Head=First;
          }

        Last->Next=NULL;
        Last->Type=Type;
        Last->Row=R;
        Last->Column=C;
        Last->String=std::string(Id);
	Last->Begin=B;
	Last->End=E;
      }
      
    void Add(short Type,char *Id,char *Id1,char *Id2,char *Id3,Parser::TokenHandle B,Parser::TokenHandle E,int R,int C)
      {
        if (First)
          {
            Last->Next=new Command;
            Last->Next->Prev=Last;
            Last=Last->Next;
          }
            
        else
          {
            First=new Command;
            First->Prev=NULL;
            Last=First;
            Head=First;
          }

        Last->Next=NULL;
        Last->Type=Type;
        Last->Row=R;
        Last->Column=C;
        Last->String=std::string(Id);
        Last->String1=std::string(Id1);
        Last->String2=std::string(Id2);
        Last->String3=std::string(Id3);
	    Last->Begin=B;
	    Last->End=E;
      }
      
    void LocalizedError(int Row,int Col,char *Msg)
      {
        strcpy(Error,"Error at row ");
        strcpy(Error+strlen(Error),Parser::IntToChars(Row));
        strcpy(Error+strlen(Error)," column ");
        strcpy(Error+strlen(Error),Parser::IntToChars(Col));
        strcpy(Error+strlen(Error),":\n");
        strcpy(Error+strlen(Error),Msg);
      }
      
    public:
      typedef Command * CommandHandle;
      enum Commands {VariableDeclaration,VariableArrayDeclaration,VariableDeclarationAffectation,OperatorDeclaration,OperatorArrayDeclaration,OperatorDeclarationAffectation,Affectation,
                     SpeciesDeclaration,LatticeDeclaration,SeedDeclaration,HamiltonianDeclaration,TemperatureDeclaration,WarmTimeDeclaration,MeasTimeDeclaration,
                     BinsDeclaration,EnsembleDeclaration,BoundariesDeclaration,SimulNameDeclaration,PopulationDeclaration,MeasureDeclaration};
      
      ParserSGF(void) {First=Last=Head=NULL;}
      ~ParserSGF(void) {if (First) Clear();}

      CommandHandle GetCommandHandle(void)
	{
	  if (First)
	    return Head;
	  else
	    return NULL;
	}

      static CommandHandle NextCommand(CommandHandle Handle)
	{
	  if (Handle && Handle->Next)
	    return Handle->Next;
	  else
	    return NULL;
	}
        
      static short GetCommandType(CommandHandle Handle)
	{
	  if (Handle)
	    return Handle->Type;
	  else
	    return -1;
	}
	
      static const char *GetIdentifier(CommandHandle Handle)
        {
	  if (Handle)
	    return Handle->CString();
	  else
	    return NULL;
	}
        
      void Clear(void)
        {
          while (First)
            {
              Last=First->Next;
              delete First;
              First=Last;
            }
        }
	
      char *ReadTokens()
        {
	  Parser::TokenHandle Handle=Parser::First();
	  Parser::TokenHandle NextHandle,End;
	  
	  while (Handle)
	    {
	      int Row,Col;
	      char *Id,*Id1,*Id2,*Id3;
	      Row=Handle->Row();
	      Col=Handle->Column();

	      switch (Handle->Type())
		{
		  
		  case Parser::Identifier:
		    Id=(char *) Handle->Value();
		    
		    if (!strcmp("Variable",Id))
		      {
			
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "Identifier in 'Variable' declaration is missing.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			
			if (Handle->Type()!=Parser::Identifier)
			  {
			    LocalizedError(Row,Col,(char *) "An identifier is expected in 'Variable' declaration.");
			    return Error;
			  }
			  
			Id=(char *) Handle->Value();
			NextHandle=Parser::NextToken(Handle);
			
			if (!NextHandle)
			  {
			    LocalizedError(Row,Col,(char *) "Either a semicolon ';', an affectation '=', or an opening bracket  is expected after identifier.");
			    return Error;
			  }
			
			Row=NextHandle->Row();
			Col=NextHandle->Column();

			switch (NextHandle->Type())
			  {
			    case Parser::Affectation:
			      Handle=Parser::NextToken(NextHandle);
			      
			      if (!Handle)
				{
				  LocalizedError(Row,Col,(char *) "An expression is expected after affectation symbol '='.");
				  return Error;
				}
				
			      Row=Handle->Row();
			      Col=Handle->Column();
			      End=Parser::GetExpressionEnd(Handle);
			      Add(VariableDeclarationAffectation,Id,Handle,End,Row,Col);
			      Handle=Parser::NextToken(End);

                              if (!Handle || *((char *) Handle->Value())!=';')
				{
				  LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Variable' affectation.");
				  return Error;
				}

			      break;
                              
                            case Parser::Bracket:
                              if (*((char *) NextHandle->Value())!='[')
                                {
                                  LocalizedError(Row,Col,(char *) "Either a semicolon ';', an affectation '=', or an opening bracket  is expected after identifier.");
                                  return Error;
                                }
                                
                              End=Parser::GetExpressionEnd(Handle);
                              Add(VariableArrayDeclaration,Id,NextHandle,End,Row,Col);
                              Handle=Parser::NextToken(End);

                              if (!Handle || *((char *) Handle->Value())!=';')
                                {
                                  LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Variable' declaration.");
                                  return Error;
                                }

                              break;
			      
			    case Parser::Separator:
                              if (*((char *) NextHandle->Value())!=';')
				{
				  LocalizedError(Row,Col,(char *) "Either a semicolon ';' or an affectation '=' is expected after identifier.");
				  return Error;
				}
				
			      Add(VariableDeclaration,Id,Row,Col);
			      Handle=NextHandle;
			      break;
				
			    default:
			      {
				LocalizedError(Row,Col,(char *) "Either a semicolon ';', an affectation '=', or an opening bracket is expected after identifier.");
				return Error;
			      }
			  }
		      }
		    
		    else if (!strcmp("Operator",Id))
		      {
			
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "Identifier in 'Operator' declaration is missing.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			
			if (Handle->Type()!=Parser::Identifier)
			  {
			    LocalizedError(Row,Col,(char *) "An identifier is expected in 'Operator' declaration.");
			    return Error;
			  }
			  
			Id=(char *) Handle->Value();
			NextHandle=Parser::NextToken(Handle);
			
			if (!NextHandle)
			  {
			    LocalizedError(Row,Col,(char *) "Either a semicolon ';' or an affectation '=' is expected after identifier.");
			    return Error;
			  }
			
			Row=NextHandle->Row();
			Col=NextHandle->Column();

			switch (NextHandle->Type())
			  {
			    case Parser::Affectation:
			      Handle=Parser::NextToken(NextHandle);
			      
			      if (!Handle)
				{
				  LocalizedError(Row,Col,(char *) "An expression is expected after affectation symbol '='.");
				  return Error;
				}
				
			      Row=Handle->Row();
			      Col=Handle->Column();
			      End=Parser::GetExpressionEnd(Handle);
			      Add(OperatorDeclarationAffectation,Id,Handle,End,Row,Col);
			      Handle=Parser::NextToken(End);

                              if (!Handle || *((char *) Handle->Value())!=';')
				{
				  LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Operator' affectation.");
				  return Error;
				}

			      break;
			      
			    case Parser::Separator:
                              if (*((char *) NextHandle->Value())!=';')
				{
				  LocalizedError(Row,Col,(char *) "Either a semicolon ';' or an affectation '=' is expected after identifier.");
				  return Error;
				}
				
			      Add(OperatorDeclaration,Id,Row,Col);
			      Handle=NextHandle;
			      break;
				
			    default:
			      {
				LocalizedError(Row,Col,(char *) "Either a semicolon ';' or an affectation '=' is expected after identifier.");
				return Error;
			      }
			  }
		      }
		    
		    else if (!strcmp("Species",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "A string is expected in 'Species' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();

			if (Handle->Type()!=Parser::String)
			  {
			    LocalizedError(Row,Col,(char *) "A string is expected in 'Species' declaration.");
			    return Error;
			  }
			  
			Id=(char *) Handle->Value();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after 'Species' name.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();

			if (!(Handle->Type()==Parser::Separator && *(char *) Handle->Value()==','))
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after 'Species' name.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "An identifier for creation operator is expected in 'Species' declaration.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();
			
			if (Handle->Type()!=Parser::Identifier)
			  {
			    LocalizedError(Row,Col,(char *) "An identifier for creation operator is expected in 'Species' declaration.");
			    return Error;
			  }
			  
			Id1=(char *) Handle->Value();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after creation operator identifier.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();

			if (!(Handle->Type()==Parser::Separator && *(char *) Handle->Value()==','))
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after creation operator identifier");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "An identifier for annihilation operator is expected in 'Species' declaration.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();
			
			if (Handle->Type()!=Parser::Identifier)
			  {
			    LocalizedError(Row,Col,(char *) "An identifier for annihilation operator is expected in 'Species' declaration.");
			    return Error;
			  }
			  
			Id2=(char *) Handle->Value();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after annihilation operator identifier.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();

			if (!(Handle->Type()==Parser::Separator && *(char *) Handle->Value()==','))
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after annihilation operator identifier.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "An identifier for number operator is expected in 'Species' declaration.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();
			
			if (Handle->Type()!=Parser::Identifier)
			  {
			    LocalizedError(Row,Col,(char *) "An identifier for number operator is expected in 'Species' declaration.");
			    return Error;
			  }
			  
			Id3=(char *) Handle->Value();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after number operator identifier.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();

			if (!(Handle->Type()==Parser::Separator && *(char *) Handle->Value()==','))
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after number operator identifier.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "An expression for maximum occupation number is expected in 'Species' declaration.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(SpeciesDeclaration,Id,Id1,Id2,Id3,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Species' declaration.");
			    return Error;
			  }
		      }
		      
		    else if (!strcmp("Lattice",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "An expression is expected in 'Lattice' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(LatticeDeclaration,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' or a coma ',' is expected in 'Lattice' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();

			if (*((char *) Handle->Value())==',')
			  {
			    if (!(Handle=Parser::NextToken(Handle)))
			      {
				LocalizedError(Row,Col,(char *) "An expression is expected in 'Lattice' declaration.");
				return Error;
			      }
			  
			    Row=Handle->Row();
			    Col=Handle->Column();
			    End=Parser::GetExpressionEnd(Handle);
			    Add(LatticeDeclaration,Handle,End,Row,Col);
			    Handle=Parser::NextToken(End);

			    if (*((char *) Handle->Value())==',')
			      {
				if (!(Handle=Parser::NextToken(Handle)))
				  {
				    LocalizedError(Row,Col,(char *) "An expression is expected in 'Lattice' declaration.");
				    return Error;
				  }
				  
				Row=Handle->Row();
				Col=Handle->Column();
				End=Parser::GetExpressionEnd(Handle);
				Add(LatticeDeclaration,Handle,End,Row,Col);
				Handle=Parser::NextToken(End);
				
				if (*((char *) Handle->Value())!=';')
				  {
				    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Lattice' declaration.");
				    return Error;
				  }
			      }
			      
			    else if (*((char *) Handle->Value())!=';')
			      {
				LocalizedError(Row,Col,(char *) "A semicolon ';' or a coma ',' is expected in 'Lattice' declaration.");
				return Error;
			      }
			  }

                        else if (*((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' or a coma ',' is expected in 'Lattice' declaration.");
			    return Error;
			  }
		      }
		      
		    else if (!strcmp("Hamiltonian",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "An expression is expected in 'Hamiltonian' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(HamiltonianDeclaration,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Hamiltonian' declaration.");
			    return Error;
			  }
		      }
		      
		    else if (!strcmp("InverseTemperature",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "An expression is expected in 'Temperature' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(TemperatureDeclaration,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Temperature' declaration.");
			    return Error;
			  }
		      }
		      
		    else if (!strcmp("WarmTime",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "An expression is expected in 'WarmTime' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(WarmTimeDeclaration,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'WarmTime' declaration.");
			    return Error;
			  }
		      }
		      
		    else if (!strcmp("MeasTime",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "An expression is expected in 'MeasTime' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(MeasTimeDeclaration,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'MeasTime' declaration.");
			    return Error;
			  }
		      }
		      
		    else if (!strcmp("Bins",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "An expression is expected in 'Bins' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(BinsDeclaration,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Bins' declaration.");
			    return Error;
			  }
		      }
		      
		    else if (!strcmp("Ensemble",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "Either 'Canonical' or 'GrandCanonical' is expected in 'Ensemble' declaration.");
			    return Error;
			  }

			End=Parser::GetExpressionEnd(Handle);
			Add(EnsembleDeclaration,Handle,End,Row,Col);
			Handle=Parser::NextToken(Handle);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Ensemble' declaration.");
			    return Error;
			  }
		      }
		      
                    else if (!strcmp("Boundaries",Id))
                      {
                        if (!(Handle=Parser::NextToken(Handle)))
                          {
                            LocalizedError(Row,Col,(char *) "Either 'Periodic' or 'Open' is expected in 'Boundaries' declaration.");
                            return Error;
                          }

                        End=Parser::GetExpressionEnd(Handle);
                        Add(BoundariesDeclaration,Handle,End,Row,Col);
                        Handle=Parser::NextToken(Handle);

                        if (!Handle || *((char *) Handle->Value())!=';')
                          {
                            LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Boundaries' declaration.");
                            return Error;
                          }
                      }
                      
		    else if (!strcmp("SimulName",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "A string is expected in 'SimulName' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();

			if (Handle->Type()!=Parser::String)
			  {
			    LocalizedError(Row,Col,(char *) "A string is expected in 'SimulName' declaration.");
			    return Error;
			  }
			  
			Id=(char *) Handle->Value();
			Add(SimulNameDeclaration,Id,Row,Col);
			Handle=Parser::NextToken(Handle);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'SimulName' declaration.");
			    return Error;
			  }
		      }
		      
		    else if (!strcmp("Seed",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "An expression is expected in 'Seed' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(SeedDeclaration,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Seed' declaration.");
			    return Error;
			  }
		      }
		      
		    else if (!strcmp("Population",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "A string is expected in 'Population' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();

			if (Handle->Type()!=Parser::String)
			  {
			    LocalizedError(Row,Col,(char *) "A string is expected in 'Population' declaration.");
			    return Error;
			  }
			  
			Id=(char *) Handle->Value();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after name of species.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();

			if (!(Handle->Type()==Parser::Separator && *(char *) Handle->Value()==','))
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after name of species.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "An expression for initial number of particles is expected in 'Population' declaration.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(PopulationDeclaration,Id,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Population' declaration.");
			    return Error;
			  }
		      }
		      
		    else if (!strcmp("Measure",Id))
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "A string is expected in 'Measure' declaration.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();

			if (Handle->Type()!=Parser::String)
			  {
			    LocalizedError(Row,Col,(char *) "A string is expected in 'Measure' declaration.");
			    return Error;
			  }
			  
			Id=(char *) Handle->Value();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after name of measurement.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();

			if (!(Handle->Type()==Parser::Separator && *(char *) Handle->Value()==','))
			  {
			    LocalizedError(Row,Col,(char *) "A coma is expected after name of measurement.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();
			Handle=Parser::NextToken(Handle);
			
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "An expression for quantity to be measured is expected in 'Measure' declaration.");
			    return Error;
			  }
			
			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(MeasureDeclaration,Id,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected after 'Measure' declaration.");
			    return Error;
			  }
		      }
		      
		    else
		      {
			if (!(Handle=Parser::NextToken(Handle)))
			  {
			    LocalizedError(Row,Col,(char *) "An affectation symbol '=' is expected after identifier.");
			    return Error;
			  }
			  
			Row=Handle->Row();
			Col=Handle->Column();

			if (Handle->Type()!=Parser::Affectation)
			  {
			    LocalizedError(Row,Col,(char *) "An affectation symbol '=' is expected after identifier.");
			    return Error;
			  }
			
			Handle=Parser::NextToken(Handle);
			      
			if (!Handle)
			  {
			    LocalizedError(Row,Col,(char *) "An expression is expected after affectation symbol '='.");
			    return Error;
			  }

			Row=Handle->Row();
			Col=Handle->Column();
			End=Parser::GetExpressionEnd(Handle);
			Add(Affectation,Id,Handle,End,Row,Col);
			Handle=Parser::NextToken(End);

                        if (!Handle || *((char *) Handle->Value())!=';')
			  {
			    LocalizedError(Row,Col,(char *) "A semicolon ';' is expected.");
			    return Error;
			  }
		      }
		    
		    break;
		    
		  default:
		    LocalizedError(Row,Col,(char *) "Either a keyword or an identifier is expected.");
		    return Error;
		}
		
	      Handle=Parser::NextToken(Handle);
	    }
	    
	  return NULL;
	}
	
      char *ExecuteCommands(void)
	{
	  CommandHandle Head=First;

	  while (Head)
	    {
	      switch (Head->Type)
		{
		  char *Err;
	    
		  case VariableDeclaration:
		    if (MathExpression::Find(Head->CString()))
		      {
			LocalizedError(Head->Row,Head->Column,(char *) "Identifier was previously declared.");
			return Error;
		      }
		      
		    if ((Err=MathExpression::Define(Head->CString(),NULL,NULL)))
		      return Err;

		    break;

                  case VariableArrayDeclaration:
                    if (MathExpression::Find(Head->CString()))
                      {
                        LocalizedError(Head->Row,Head->Column,(char *) "Identifier was previously declared.");
                        return Error;
                      }

                    if ((Err=MathExpression::Define(Head->CString(),Head->Begin,Head->End)))
                      return Err;

                    break;
              
		  case VariableDeclarationAffectation:
		    if (MathExpression::Find(Head->CString()))
		      {
			LocalizedError(Head->Row,Head->Column,(char *) "Identifier was previously declared.");
			return Error;
		      }

		    if ((Err=MathExpression::Define(Head->CString(),Head->Begin,Head->End)))
		      return Err;

		    break;
	      
		  case OperatorDeclaration:
		    if (MathExpression::Find(Head->CString()))
		      {
			LocalizedError(Head->Row,Head->Column,(char *) "Identifier was previously declared.");
			return Error;
		      }
		      
		    if ((Err=MathExpression::Define(Head->CString(),NULL,NULL)))
		      return Err;

		    break;

		  case OperatorDeclarationAffectation:
		    if (MathExpression::Find(Head->CString()))
		      {
			LocalizedError(Head->Row,Head->Column,(char *) "Identifier was previously declared.");
			return Error;
		      }

		    if ((Err=MathExpression::Define(Head->CString(),Head->Begin,Head->End)))
		      return Err;

		    break;

		  case Affectation:
		    if ((Err=MathExpression::SetValue(Head->CString(),Head->Begin,Head->End)))
		      {
			LocalizedError(Head->Row,Head->Column,Err);
			return Error;
		      }
		      
		    break;

		  case SpeciesDeclaration:
		    if (MathExpression::Find(Head->CString()))
		      {
			LocalizedError(Head->Row,Head->Column,(char *) "Name of species was previously declared.");
			return Error;
		      }

		    if (MathExpression::Find(Head->CString1()))
		      {
			LocalizedError(Head->Row,Head->Column,(char *) "Name of creation operator was previously declared.");
			return Error;
		      }

		    if (MathExpression::Find(Head->CString2()))
		      {
			LocalizedError(Head->Row,Head->Column,(char *) "Name of annihilation operator was previously declared.");
			return Error;
		      }

		    if (MathExpression::Find(Head->CString3()))
		      {
			LocalizedError(Head->Row,Head->Column,(char *) "Name of number operator was previously declared.");
			return Error;
		      }

		    if ((Err=MathExpression::AddSpecies(Head->CString(),Head->CString1(),Head->CString2(),Head->CString3(),Head->Begin,Head->End)))
		      return Err;
		    
		    break;

		  case LatticeDeclaration:
		    if ((Err=MathExpression::Lattice(Head->Begin,Head->End)))
		      return Err;
		    
		    break;

		  case SeedDeclaration:
		    if ((Err=MathExpression::Define("#Seed",Head->Begin,Head->End)))
		      return Err;
		    
		    RNG::Initialize((int) MathExpression::GetValue("#Seed").Re());
		    break;

		  case HamiltonianDeclaration:
		    if ((Err=MathExpression::Define("#Hamiltonian",Head->Begin,Head->End)))
		      return Err;
		    break;

		  case TemperatureDeclaration:
		    if ((Err=MathExpression::Define("#InverseTemperature",Head->Begin,Head->End)))
		      return Err;
		    break;

		  case WarmTimeDeclaration:
		    if ((Err=MathExpression::Define("#WarmTime",Head->Begin,Head->End)))
		      return Err;
		    break;

		  case MeasTimeDeclaration:
		    if ((Err=MathExpression::Define("#MeasTime",Head->Begin,Head->End)))
		      return Err;
		    break;

		  case BinsDeclaration:
		    if ((Err=MathExpression::Define("#Bins",Head->Begin,Head->End)))
		      return Err;
		    break;

		  case EnsembleDeclaration:
		    if ((Err=MathExpression::Define("#Ensemble",Head->Begin,Head->End)))
		      return Err;
		    break;

                  case BoundariesDeclaration:
                    if ((Err=MathExpression::Define("#Boundaries",Head->Begin,Head->End)))
                      return Err;
                    break;

		  case SimulNameDeclaration:
		    MathExpression::SetSimulName(Head->CString());
		    break;

		  case PopulationDeclaration:
		    if ((Err=MathExpression::SetPopulation(Head->CString(),Head->Begin,Head->End)))
		      return Err;
		    break;
	      
		  case MeasureDeclaration:
		    if ((Err=MathExpression::Define(Head->CString(),Head->Begin,Head->End)))
		      return Err;
		    
		    MathExpression::AddMeasurable(Head->CString());
		    MathExpression::BuildMeasurable(Head->CString());
		    break;
	      
		  default:
		    return (char *) "Error: Type of command is unknoww (This error is not supposed to happen, probably a bug).";
		}
	      
	      Head=Head->Next;
	    }
	    
	  if (!MathExpression::GetNumSites())
	    return (char *) "The lattice has not been defined.";
	    
	  else if (!MathExpression::GetNumSpecies())
	    return (char *) "One species of particles at least must be defined.";
	    
	  else if (!MathExpression::Find("#Seed"))
	    return (char *) "The seed of the random number generator must be defined.";
	    
	  else if (!MathExpression::Find("#Hamiltonian"))
	    return (char *) "The Hamiltonian must be defined.";
	    
	  else if (!MathExpression::Find("#InverseTemperature"))
	    return (char *) "The inverse temperature must be defined.";
	    
	  else if (!MathExpression::Find("#WarmTime"))
	    return (char *) "The time of warming-up must be defined.";
	    
	  else if (!MathExpression::Find("#MeasTime"))
	    return (char *) "The time of measurements must be defined.";
	    
	  else if (!MathExpression::Find("#Bins"))
	    return (char *) "The number of bins must be defined.";
	    
	  else if (!MathExpression::Find("#Ensemble"))
	    return (char *) "The ensemble of the simulation must be defined.";
	    
          else if (!MathExpression::Find("#Boundaries"))
            return (char *) "The boundary conditions of the simulation must be defined.";
            
          else if (!*MathExpression::GetSimulName())
            return (char *) "A name must be given to the simulation.";

            
	  return NULL;
		}
  };

#endif

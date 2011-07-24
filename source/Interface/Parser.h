// **************************************************************************************
// * This class transforms the input file or input buffer into a linked list of tokens. *
// * __________________________________________________________________________________ *
// * Dr. Valy G. Rousseau - December 2nd, 2010                                          *
// **************************************************************************************
//
// A Token is accessible either via an internal pointer (the current token, see the methods below) or via
// a handle:
//   Parser::TokenHandle MyHandle;
//
// Tokens can be of the folowing types:
//   - Parser::Identifier
//   - Parser::String
//   - Parser::Number
//   - Parser::Operator
//   - Parser::Parenthesis
//   - Parser::Bracket
//   - Parser::Brace
//   - Parser::Angle
//   - Parser::Separator
//   - Parser::Affectation
//   - Parser::SpecialCharacter
//
// THE METHODS
//
//   char *ReadFile(const char *FileName) 
//   Builds a list of tokens from FileName. The returned value is a pointer
//   onto an error message, if any, otherwise it is NULL.
//
//   char *ReadBuffer(const char *Buffer) 
//   Builds a list of tokens from Buffer. The returned value is a pointer
//   onto an error message, if any, otherwise it is NULL.
//
//   Parser::TokenHandle GetTokenHandle(void)
//   Returns a handle for the current token.
//   A NULL handle indicates that no tokens are available.
//
//   short GetTokenType(void)
//   Returns the type of current token. A negative value indicates that no tokens
//   are available.
//
//   static short GetTokenType(Parser::TokenHandle Handle)
//   Returns the type of token given by Handle (which is assumed to be a valid handle).
//
//   short GetTokenRow(void)
//   Returns the row number in the input file or the input buffer corresponding to current token.
//   A negative value indicates that no tokens are available.
//
//   static short GetTokenRow(Parser::TokenHandle Handle)
//   Returns the row number in the input file or the input buffer corresponding to the token
//   given by Handle (which is assumed to be a valid handle).
//
//   short GetTokenColumn(void)
//   Returns the column number in the input file or the input buffer corresponding to current token.
//   A negative value indicates that no tokens are available.
//
//   static short GetTokenColumn(Parser::TokenHandle Handle)
//   Returns the column number in the input file or the input buffer corresponding to the token
//   given by Handle (which is assumed to be a valid handle).
//
//   void *GetTokenValue(void)
//   Returns a pointer onto the value of current token. A cast of the returned pointer must be done
//   depending on the type of the token:
//     If type is Number ---> Cast to (double *)
//     Otherwise ---> Cast to (char *)
//   A NULL pointer indicates that no tokens are available.
//
//   static void *GetTokenValue(Parser::TokenHandle Handle)
//   Returns a pointer onto the value of token given by Handle. A cast of the returned pointer must be done
//   depending on the type of the token:
//     If type is Number ---> Cast to (double *)
//     Otherwise ---> Cast to (char *)
//
//   bool Compare(char C)
//   Returns true if value of current token is equal to C, false otherwise.
//
//   static bool Compare(Parser::TokenHandle Handle,char C)
//   Returns true if value of token given by Handle is equal to C, false otherwise.
//
//   bool Compare(char *String)
//   Returns true if value of current token is equal to String, false otherwise.
//
//   static bool Compare(Parser::TokenHandle Handle,char *String)
//   Returns true if value of token given by Handle is equal to String, false otherwise.
//
//   Parser::TokenHandle GetExpressionEnd(void)
//   Returns last token of current expression, excluding the end mark token. The end mark token is
//   either a separator ';' or ',' or an extra closing parenthesis ')', bracket ']', brace '}', or angle '>'.
//
//   static Parser::TokenHandle GetExpressionEnd(Parser::TokenHandle Handle)
//   Returns last token of expression starting with Handle, excluding the end mark token. The end mark token is
//   either a separator ';' or ',' or an extra closing parenthesis ')', bracket ']', brace '}', or angle '>'.
//
//   bool MoveTo(Parser::TokenHandle Handle)
//   Set current token to Handle. Returns true if Handle exists, false otherwise.
//
//   bool NextToken(void)
//   Changes current token to next token. Returns true if next token exists,
//   false otherwise.
//
//   static Parser::TokenHandle NextToken(Parser::TokenHandle Handle)
//   Returns a handle for token after Handle.
//   The returned handle is NULL if no next token is available.
//
//   bool PrevToken(void)
//   Changes current token to previous token. Returns true if previous token exists,
//   false otherwise.
//
//   static Parser::TokenHandle PrevToken(Parser::TokenHandle Handle)
//   Returns a handle for token before Handle.
//   The returned handle is NULL if no previous token is available.
//
//   void Clear(void)
//   Clears the list of tokens. Must be called before making a new call to ReadFile
//   or ReadBuffer.
//
//   static void Display(TokenHandle Begin,TokenHandle End)
//   Displays the values of the tokens from Begin to End.

#ifndef ParserDef
#define ParserDef

#include <fstream>
#include <vector>
#include <string>

class Parser
{
    static char Error1[];
    static char Error2[];
    static char Error3[256];

    class TokenData { 
        
        short _Type;
        int _Row,_Column;
        double _Real;
        std::string _String;

    public:
        TokenData(short Type,const void *Value,int R,int C) : _Type(Type), _Row(R), _Column(C) {
                        
            if(Type<NTokenTypes) {
                if(Type==Number)
                  _Real=*static_cast<const double *>(Value); 
                else if (Type==String || Type==Identifier) 
                  _String=std::string(static_cast<const char *>(Value)); 
                else {
                  char temp_string[2]={*((char*)Value),0};
                  _String=std::string(temp_string);
                }
            }
            else {
                cout << "Error in Parser::Add(short,const void *,int,int) function:\n";
                cout << "Type " << Type << " is not supposed to occur in this function!" << endl;
            }
            
        }
        TokenData(const TokenData &o) : _Type(o._Type), _Row(o._Row), _Column(o._Column), _Real(o._Real),_String(o._String) {}
        ~TokenData() {}
        inline short Type() {return _Type;}
        inline int Row() {return _Row;}
        inline int Column() {return _Column;}
        const char * GetString() {return _String.c_str();}

        bool Compare(char C) const { return _Type!=Number && _String[0]==C; }
        bool Compare(std::string S) const { return _Type!=Number && std::string(_String)==S; }



        const void *Value() {
            if (_Type==Number)
                return &_Real;
            else
                return _String.c_str();
        }

    };
     
    std::vector<TokenData> TokenVector;
  



    static void Add(short Type,const void *Value,int R,int C)
    {
        if (_First)
        {
            _Last->_Next=new Token(Type,Value,R,C);
            _Last->_Next->_Prev=_Last;
            _Last=_Last->_Next;
        }

        else
        {
            _First=new Token(Type,Value,R,C);
            _First->_Prev=NULL;
            _Last=_First;
        }

        _Last->_Next=NULL;

    }
    
    
public:  
    struct Token : public TokenData {
      Token *_Next,*_Prev;
      Token(short Type,const void *Value,int R,int C) : TokenData(Type,Value,R,C), _Next(NULL), _Prev(NULL) {}
      Token * NextToken() {return _Next;}
      Token * PrevToken() {return _Prev;}
    };
 
    
    
    typedef Token * TokenHandle;
    enum TokenTypes {Number,Identifier,String,Operator,Parenthesis,Bracket,Brace,Angle,Separator,Affectation,SpecialCharacter,NTokenTypes};

    static Token* _First;
    static Token* _Last;



    Parser(void) {_First=NULL; _Last=NULL;}
    ~Parser(void) {
        Token * Temp;
        while (_First) {
            Temp=_First->_Next;
            delete _First;
            _First=Temp;
        }
    }


    static Token *First() {return _First;}
    static Token *Last() {return _Last;}
    
    static TokenHandle NextToken(TokenHandle Handle) {return Handle->_Next;}
    static TokenHandle PrevToken(TokenHandle Handle) {return Handle->_Prev;}

    static TokenHandle GetExpressionEnd(TokenHandle Handle)  {
        char C;
        int Stack=0;

        do {
            if (Handle->Type()!=Number)
                C=*Handle->GetString();
            else
                C=0;

            if (!Stack)
                if (C==')' || C==']' || C=='}' || C=='>' || C==';' || C==',')
            {
                
                return PrevToken(Handle);
            }

            if (C=='(' || C=='[' || C=='{' || C=='<')
                Stack++;
            else if (C==')' || C==']' || C=='}' || C=='>')
                Stack--;

        } while ((Handle=NextToken(Handle)));

        return Handle;
    }



    static const char *IntToChars(int N) {
        
        static char Buffer[12];
        bool Negative=false;
        Buffer[11]=0;
        char *Ptr=Buffer+11;

        if (N<0)
        {
            N=-N;
            Negative=true;
        }

        do
        {
            *(--Ptr)=N%10+48;
            N/=10;
        } while (N);

        if (Negative)
            *(--Ptr)='-';

        return Ptr;
    }

    static void LocalizedError(const char *Message,int Row,int Col) {
        strcpy(Error3,Message);
        strcpy(Error3+strlen(Error3),", row ");
        strcpy(Error3+strlen(Error3),IntToChars(Row));
        strcpy(Error3+strlen(Error3)," column ");
        strcpy(Error3+strlen(Error3),IntToChars(Col));
        strcpy(Error3+strlen(Error3),".");
    }


    static char *ReadFile(const char *FileName) {

        fstream File;
        File.open(FileName,ios::in|ios::binary|ios::ate);

        if (First())
          return Error2;

        if (File.is_open())
        {
            streamoff FileSize=File.tellg();
            File.seekg(0);
            char *Buffer=new char[FileSize+1];
            File.read(Buffer,FileSize);
            File.close();
            Buffer[FileSize]=0;
            char *Result=ReadBuffer(Buffer);
            delete []Buffer;
            return Result;
        }

        else
            return Error1;
    }

    static char *ReadBuffer(const char *Buffer)
      {
        int Row=1;
        int Column=1;
        char Temp[256];
        char *TempPtr;
        
        if (First())
          return Error2;

        while (*Buffer)
          {
            if (*Buffer==' ' || *Buffer=='\t') {Buffer++; Column++;}
            else if (*Buffer=='\n') {Buffer++; Row++; Column=1;}
            else if (*Buffer=='#') {while (*Buffer && *Buffer!='\n') Buffer++;}
              
            else if ((*Buffer>='A' && *Buffer<='Z') || (*Buffer>='a' && *Buffer<='z'))
              {
                int Col=Column;
                TempPtr=Temp;
                
                do
                  {
                    *(TempPtr++)=*(Buffer++);
                    Column++;
                  } while ((*Buffer>='A' && *Buffer<='Z') || (*Buffer>='a' && *Buffer<='z') || (*Buffer>='0' && *Buffer<='9'));
                
                *TempPtr=0;
                Add(Identifier,Temp,Row,Col);
              }

            else if (*Buffer=='"')
              {
	  int Col=Column;
                TempPtr=Temp;
                Buffer++;
                Column++;
          
                while (*Buffer && *Buffer!='"' && TempPtr-Temp<255)
                  {
                    *(TempPtr++)=*(Buffer++);
                    Column++;
                  }
                  
                *TempPtr=0;
                
                if (*Buffer=='"')
                  {
                    Column++;
                    Buffer++;
                  }
                  
                else
                  {
	      LocalizedError("Error: Found opening quotation mark without closing",Row,Col);
                    return Error3;
                  }
            
                Add(String,Temp,Row,Col);
              }

            else if (*Buffer>='0' && *Buffer<='9')
              {
                int Col=Column;
                double Value=0.0;

                do
                  {
                    Value=Value*10+*(Buffer++)-48;
                    Column++;
                  } while (*Buffer>='0' && *Buffer<='9');
                    
                if (*Buffer=='.')
                  {
                    int DecValue=0;
                    Buffer++;
                    Column++;
                    
                    if (*Buffer>='0' && *Buffer<='9')
                      {
                        int NumDec=0;
                        
                        do
                          {
                            DecValue=DecValue*10+*(Buffer++)-48;
                            NumDec++;
                            Column++;
                          } while (*Buffer>='0' && *Buffer<='9');

                          for (int i=0;i<NumDec;i++)
                          Value*=10;

                        Value+=DecValue;
                  
                        for (int i=0;i<NumDec;i++)
                          Value/=10;
                      }
                      
                    else
                      {
		  LocalizedError("Error: A digit is expected after decimal point",Row,Col);
                        return Error3;
                      }
                  }
                  
                if (*Buffer=='e' || *Buffer=='E')
                  {
                    int Sign=1;
                    int Exp=0;
                    Buffer++;
                    Column++;
                    
                    if (*Buffer=='+' || *Buffer=='-')
                      {
                        if (*Buffer=='-')
                          Sign=-1;
                          
                        Buffer++;
                        Column++;
                      }
                    
                    if (*Buffer>='0' && *Buffer<='9')
                      {
                        do
                          {
                            Exp=Exp*10+*(Buffer++)-48;
                            Column++;
                          } while (*Buffer>='0' && *Buffer<='9');
                          
                        for (int i=0;i<Exp;i++)
                          if (Sign>0)
                            Value*=10;
                          else
                            Value/=10;
                      }
                      
                    else
                      {
		  LocalizedError("Error: A digit is expected after exponential",Row,Col);
                        return Error3;
                      }
                  }

                Add(Number,&Value,Row,Col);
              }
             
            else if (*Buffer=='=')
              {
                Add(Affectation,Buffer,Row,Column);
                Buffer++;
                Column++;
              }

            else if (*Buffer=='+' || *Buffer=='-' || *Buffer=='*' || *Buffer=='/' || *Buffer=='%' || *Buffer=='^')
              {
                Add(Operator,Buffer,Row,Column);
                Buffer++;
                Column++;
              }
              
            else if (*Buffer=='(' || *Buffer==')')
              {
                Add(Parenthesis,Buffer,Row,Column);
                Buffer++;
                Column++;
              }
              
            else if (*Buffer=='[' || *Buffer==']')
              {
                Add(Bracket,Buffer,Row,Column);
                Buffer++;
                Column++;
              }
              
            else if (*Buffer=='{' || *Buffer=='}')
              {
                Add(Brace,Buffer,Row,Column);
                Buffer++;
                Column++;
              }
              
            else if (*Buffer=='<' || *Buffer=='>')
              {
                Add(Angle,Buffer,Row,Column);
                Buffer++;
                Column++;
              }
              
            else if (*Buffer==',' || *Buffer==';' || *Buffer=='.' || *Buffer==':' || *Buffer=='|')
              {
                Add(Separator,Buffer,Row,Column);
                Buffer++;
                Column++;
              }
              
            else
              {
                Add(SpecialCharacter,Buffer,Row,Column);
                Buffer++;
                Column++;
              }
          }
          
        return CheckParenthesis();            
      }
      
    static void Display(TokenHandle Begin,TokenHandle End)
      {
        while (Begin && Begin!=End->NextToken())
          {
            switch (Begin->Type())
              {
                case Number:
                  cout << *(double *) Begin->Value() << endl;
                  break;
                  
                default:
                  cout << (char *) Begin->Value() << endl;
              }
              
            Begin=Begin->NextToken();
          }
      }

    static char * CheckParenthesis(void)
      {
        int Stack1=0;
        int Stack2=0;
        int Stack3=0;
        int Stack4=0;
        TokenHandle Head=First();
        
        while (Head && Head->Type()!=Number)
          {
            char C=*(char *) Head->Value();
            int Row=Head->Row();
            int Col=Head->Column();
            
            switch (C)
              {
                case '(': Stack1++; break;
                case '[': Stack2++; break;
                case '{': Stack3++; break;
                case '<': Stack4++; break;
                case ')': Stack1--; if (Stack1<0) {LocalizedError("Error: Found an extra closing parenthesis ')'",Row,Col); return Error3;} break;
                case ']': Stack2--; if (Stack2<0) {LocalizedError("Error: Found an extra closing bracket ']'",Row,Col); return Error3;} break;
                case '}': Stack3--; if (Stack3<0) {LocalizedError("Error: Found an extra closing brace '}'",Row,Col); return Error3;} break;
                case '>': Stack4--; if (Stack4<0) {LocalizedError("Error: Found an extra closing angle '>'",Row,Col); return Error3;} break;
              }
              
            Head=NextToken(Head);
          }
          
        Stack1=0;
        Stack2=0;
        Stack3=0;
        Stack4=0;
        Head=Last();
        
        while (Head && Head->Type()!=Number)
          {
            char C=*(char *) Head->Value();
            int Row=Head->Row();
            int Col=Head->Column();
            
            switch (C)
              {
                case '(': Stack1++; if (Stack1>0) {LocalizedError("Error: Found an extra opening parenthesis '('",Row,Col); return Error3;} break;
                case '[': Stack2++; if (Stack2>0) {LocalizedError("Error: Found an extra opening bracket '['",Row,Col); return Error3;} break;
                case '{': Stack3++; if (Stack3>0) {LocalizedError("Error: Found an extra opening brace '{'",Row,Col); return Error3;} break;
                case '<': Stack4++; if (Stack4>0) {LocalizedError("Error: Found an extra opening angle '<'",Row,Col); return Error3;} break;
                case ')': Stack1--; break;
                case ']': Stack2--; break;
                case '}': Stack3--; break;
                case '>': Stack4--; break;
              }
              
            Head=PrevToken(Head);
          }
          
        return NULL;
      }

};

Parser::Token *Parser::_First;
Parser::Token *Parser::_Last;

char Parser::Error1[]="Error: File does not exist!";
char Parser::Error2[]="Error: Tokens have not been cleared!";
char Parser::Error3[256];

#endif

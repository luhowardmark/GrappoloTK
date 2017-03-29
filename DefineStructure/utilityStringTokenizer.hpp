/* Author: Arijit Tarafdar < aritar@gmail.com >                         */
/* Used with permission by the author                                   */
// **************************************************************************************************
// Grappolo: A C++ library for parallel graph community detection
// Hao Lu, Ananth Kalyanaraman (hao.lu@wsu.edu, ananth@eecs.wsu.edu) Washington State University
// Mahantesh Halappanavar (hala@pnnl.gov) Pacific Northwest National Laboratory
//
// For citation, please cite the following paper:
// Lu, Hao, Mahantesh Halappanavar, and Ananth Kalyanaraman. 
// "Parallel heuristics for scalable community detection." Parallel Computing 47 (2015): 19-37.
//
// **************************************************************************************************
// Copyright (c) 2016. Washington State University ("WSU"). All Rights Reserved.
// Permission to use, copy, modify, and distribute this software and its documentation
// for educational, research, and not-for-profit purposes, without fee, is hereby
// granted, provided that the above copyright notice, this paragraph and the following
// two paragraphs appear in all copies, modifications, and distributions. For
// commercial licensing opportunities, please contact The Office of Commercialization,
// WSU, 280/286 Lighty, PB Box 641060, Pullman, WA 99164, (509) 335-5526,
// commercialization@wsu.edu<mailto:commercialization@wsu.edu>, https://commercialization.wsu.edu/

// IN NO EVENT SHALL WSU BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT, SPECIAL, INCIDENTAL,
// OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS, ARISING OUT OF THE USE OF
// THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF WSU HAS BEEN ADVISED OF THE
// POSSIBILITY OF SUCH DAMAGE.

// WSU SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
// WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE. THE SOFTWARE AND
// ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED HEREUNDER IS PROVIDED "AS IS". WSU HAS NO
// OBLIGATION TO PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
// **************************************************************************************************

#ifndef _string_Tokenizer_
#define _string_Tokenizer_

//I/O
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>

using namespace std;

class StringTokenizer
{
 private:

  string DelimiterString;
  string InputString;
  string TokenString;

 public:

  long CountTokens();			// ***Public Function No. 1***
  long CountTokens(char *);		// ***Public Function No. 2***
  
  string GetDelimiterString() const;	// ***Public Function No. 3***
  string GetFirstToken();		// ***Public Function No. 4***
  string GetInputString() const;	// ***Public Function No. 5***
  string GetLastToken();		// ***Public Function No. 6***
  string GetNextToken();		// ***Public Function No. 7***
  string GetNextToken(char *);		// ***Public Function No. 8***
  string GetToken(long);			// ***Public Function No. 9***

  long HasMoreTokens();			// ***Public Function No. 10***
  long HasMoreTokens(char *);		// ***Public Function No. 11***
  
  long SetInputString(char *);		// ***Public Function No. 12***
  long SetDelimiterString(char *);	// ***Public Function No. 13***
  
  StringTokenizer();			// ***Public Function No. 14***
  StringTokenizer(char *);		// ***Public Function No. 15***
  StringTokenizer(char *, char *);	// ***Public Function No. 16***
  StringTokenizer(string, char *);	// ***Public Function No. 17***
  StringTokenizer(string, string);	// ***Public Function No. 18***
  ~StringTokenizer();			// ***Public Function No. 19***

};
/* ------------------------------------------------------------------------- */
inline
long StringTokenizer::CountTokens()
{
  long TokenCounter = 1;

  long DelimiterPosition;
  
  long LastPosition;

  long TokenStringLength = TokenString.size();
  long DelimiterStringLength = DelimiterString.size();

  string DelimiterSubString;

  if(TokenStringLength == 0)
  {
    return(0);
  }

  if(DelimiterStringLength == 0)
  {
    return(1);
  }

  DelimiterPosition = 0;
  LastPosition = 0;

  while(1)
  {
    
    DelimiterPosition = TokenString.find(DelimiterString, DelimiterPosition);

    if(DelimiterPosition == 0)
    { 
      DelimiterPosition += DelimiterStringLength;

      continue;
    }

    if((DelimiterPosition < 0) || (DelimiterPosition == TokenStringLength))
    { 
      return(TokenCounter);
    }
    
    if(DelimiterStringLength != (DelimiterPosition - LastPosition))
    {
      //      cout<<"Delimiter Position = "<<DelimiterPosition<<endl;

      TokenCounter++;
    }

    LastPosition = DelimiterPosition;
    
    DelimiterPosition += DelimiterStringLength;

  }

  return(TokenCounter);

}

// ***Public Function No. 2***
inline
long StringTokenizer::CountTokens(char * DelimiterChar)
{
  SetDelimiterString(DelimiterChar);

  return(CountTokens());


}

// ***Public Function No. 3***
inline
string StringTokenizer::GetDelimiterString() const
{
  return(DelimiterString);
}

// ***Public Function No. 4***
inline
string StringTokenizer::GetFirstToken()
{
  long TokenCount = 0;

  string StringToken;

  TokenString = InputString;

  while(HasMoreTokens())
  {
    if(TokenCount == 1)
    {
      break;
    }

    StringToken = GetNextToken();

    TokenCount++;

  }

  return(StringToken);
}

// ***Public Function No. 5***
inline
string StringTokenizer::GetInputString() const
{
  return(InputString);
}

// ***Public Function No. 6***
inline
string StringTokenizer::GetLastToken()
{
  string StringToken;

  TokenString = InputString;

  while(HasMoreTokens())
  {
    StringToken = GetNextToken();
  }

  return(StringToken);
  
}

// ***Public Function No. 7***
inline
string StringTokenizer::GetNextToken()
{
  string Token;
  
  long DelimiterPosition;
  
  long TokenStringLength = TokenString.size();
  long DelimiterStringLength = DelimiterString.size();

  string DelimiterSubString;

  if((TokenStringLength == 0))
  {
    return(NULL);
  }

  if(DelimiterStringLength == 0)
  {
    return(InputString);
  }
  
  DelimiterPosition = TokenString.find(DelimiterString);

  if(DelimiterPosition == 0)
  {
    while(1)
    {
      if(TokenString.substr(0, DelimiterStringLength) == DelimiterString)
      {
	TokenString.erase(0, DelimiterStringLength);
      }
      else
      {
	break;
      }
    }

    DelimiterPosition = TokenString.find(DelimiterString);
  }
    
  if(DelimiterPosition < 0)
  {
    Token = TokenString;

    TokenString.erase();
  }
  else
  {

    Token = TokenString.substr(0, DelimiterPosition);

    TokenString.erase(0, DelimiterPosition+DelimiterStringLength);
  

    DelimiterPosition = 0;

    while(1)
    {
      if(TokenString.substr(0, DelimiterStringLength) == DelimiterString)
      {
	TokenString.erase(0, DelimiterStringLength);
      }
      else
      {
	break;
      }
    }
    
  }

  return(Token);
}  

// ***Public Function No. 8***
inline
string StringTokenizer::GetNextToken(char * DelimiterChar)
{
  SetDelimiterString(DelimiterChar);

  return(GetNextToken());
}

// ***Public Function No. 9***
inline
string StringTokenizer::GetToken(long TokenPosition)
{
  long TokenCount = 0;

  string StringToken;

  TokenString = InputString;

  while(HasMoreTokens())
  {
    if(TokenCount == TokenPosition)
    {
      break;
    }

    StringToken = GetNextToken();

    TokenCount++;
  }

  return(StringToken);
}

// ***Public Function No. 10***
inline
long StringTokenizer::HasMoreTokens()
{
  return(CountTokens());
}

// ***Public Function No. 11***
inline
long StringTokenizer::HasMoreTokens(char * DelimiterChar)
{
  SetDelimiterString(DelimiterChar);

  return(HasMoreTokens());
}

// ***Public Function No. 12***
inline
long StringTokenizer::SetDelimiterString(char * DelimiterChar)
{
  string TempDelimiterString(DelimiterChar);

  DelimiterString = TempDelimiterString;

  return(0);
}

// ***Public Function No. 13***
inline
long StringTokenizer::SetInputString(char * InputChar)
{
  string TempInputString(InputChar);

  InputString = TempInputString;
  TokenString = InputString;

  return(0);
}

// ***Public Function No. 14***
inline
StringTokenizer::StringTokenizer()
{

}

// ***Public Function No. 15***
inline
StringTokenizer::StringTokenizer(char * InputChar)
{
  string TempInputString(InputChar);
  
  InputString = TempInputString;
  TokenString = InputString;

}

// ***Public Function No. 16***
inline
StringTokenizer::StringTokenizer(char * InputChar, char * DelimiterChar)
{
  string TempInputString(InputChar);
  string TempDelimiterString(DelimiterChar);

  InputString = TempInputString;
  TokenString = InputString;

  DelimiterString = TempDelimiterString;

}

// ***Public Function No. 17***
inline
StringTokenizer::StringTokenizer(string InputChar, char * DelimiterChar)
{
  string TempDelimiterString(DelimiterChar);

  InputString = InputChar;
  TokenString = InputString;

  DelimiterString = TempDelimiterString;

}

// ***Public Function No. 18***
inline
StringTokenizer::StringTokenizer(string InputChar, string DelimiterChar)
{
  InputString = InputChar;
  TokenString = InputString;

  DelimiterString = DelimiterChar;

}

// ***Public Function No. 19***
inline
StringTokenizer::~StringTokenizer()
{


}

#endif

#include "NewickCode.h"

using namespace std;
using std::string;

/**
* NewickCode default constructor, initializes the class NewickCode. 
*/
NewickCode::NewickCode ()
{
	code="";
}

/**
* NewickCode constructor, initializes the class NewickCode. Parameters:
*	_code: string containing the newick code of a phylogenetic tree.
*/
NewickCode::NewickCode (string _code)
{
	code=_code;
}

/**
* Method that returns the newick code stored as a string variable. 
*/
string NewickCode::getNewick ()
{
	return code;
}

/**
	Method that sets the string newick code. Parameters:
	_code: string containing the newick code of a phylogenetic tree.
*/
void NewickCode::setNewick (string _code)
{
	code=_code;
}

/**
	Method that sets the string newick code. Parameters:
	_code: char* containing the newick code of a phylogenetic tree.
*/
void NewickCode::setNewick2 (char* _code)
{
	code=_code;
}

/**
	Method that writes into a file the string newick code. Parameters:
	file: FILE* associated to the output file (PRE = the file has to be be already open).
*/
int NewickCode::writeNewick (FILE* file)
{
	if (file == NULL) return 0;
	for (int i=0; i < code.length(); i++)
  	{
    		fwrite(&code[i], sizeof(char), 1, file);
  	}
	return 1;
}


/**
	Method that read the string newick code from a file. Parameters:
	file: FILE* associated to the input file (PRE = the file has to be be already open).
*/
int NewickCode::readNewick (FILE* file)
{
	char c;
	if (file == NULL) return 0;
	c = fgetc(file);
	code="";
	while (c != EOF && c != ';')
	{
		code+=c;		
		c = fgetc(file);
	}
	if (c == ';') code+=c;	
	return 1;
	
}

/**
	Method that clones an object of the class NewickCode. Parameters:
	_code: NewickCode instance
*/
void NewickCode::clone (const NewickCode & _code)
{
	code=_code.code;
}


/**
	Method that prints to console the newick code.
*/
void NewickCode::printCode()
{
	cout<<code<<endl;
}

/**	Default estructor 	*/
NewickCode::~NewickCode()
{}

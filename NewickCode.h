//Accelerating the Phylogenetic Parsimony Function on Heterogeneous Systems
//Authors: Sergio Santander-Jimenez, Aleksandar Ilic, Leonel Sousa, and Miguel A. Vega-Rodriguez
//Header of the NewickCode class, implementation of char-based representation of phylogenetic trees 


#ifndef _NEWICKCODE_H_
#define _NEWICKCODE_H_
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <string>

using std::string;

using namespace std;

class NewickCode
{
	string code; //Newick code of the stored phylogenetic tree
	public:
		NewickCode ();  //NewickCode default constructor, initializes the class NewickCode.
		NewickCode (string _code); //NewickCode constructor, initializes the class NewickCode with the newick code specified in _code.
		string getNewick (); //Returns the stored newick code
		void setNewick (string _code); //Assigns newick code from string variable
		void setNewick2 (char* _code); //Assigns newick code from char*
		int writeNewick (FILE* file); //Writes the newick code in a file
		int readNewick (FILE* file);  //Read the newick code from a file
		void clone (const NewickCode & _code); //Copy the contents of _code into the object
		void printCode(); //Prints in screen the stored newick code
		~NewickCode(); //Default destructor of the NewickCode class
};

#endif	//_NEWICKCODE_H_

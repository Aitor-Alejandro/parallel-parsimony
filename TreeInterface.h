//Accelerating the Phylogenetic Parsimony Function on Heterogeneous Systems
//Authors: Sergio Santander-Jimenez, Aleksandar Ilic, Leonel Sousa, and Miguel A. Vega-Rodriguez
//Header of the TreeInterface class, implementation of methods to interact with the PhylogeneticTree class 


#ifndef _TREEINTERFACE_H_
#define _TREEINTERFACE_H_
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include "PhylogeneticTree.h"
using namespace std;

class TreeInterface //This is a basic implementation, designed to allow interactions with other implementations of phylogenetic trees in future releases
{	
	int id; //identifier of the TreeInterface (used if a different numeration is required)
	PhylogeneticTree tree; //PhylogeneticTree Instance
	
	public:
		TreeInterface (); //Default constructor
		TreeInterface (int _id, PhylogeneticTree &_tree); //Parametrized constructor
		
		void setTree (PhylogeneticTree &_tree);	//Sets the PhylogeneticTree object	
		void setId(int _id); //Sets	the identifier of the TreeInterface object		
		PhylogeneticTree* getTree(); //Returns the PhylogeneticTree object		
		int getId(); //Returns the identifier of the TreeInterface object		
		int setScores (); //Evaluates the associated phylogenetic tree using the function implementation defined in the PhylogeneticTree object
		int getParsimony(); //Returns the parsimony score calculated by using the function implementation defined in the PhylogeneticTree object		
		void printTree (); //Prints the TreeInterface data
		
		~TreeInterface(); //Destructor
};


#endif	//_TREEINTERFACE_H_

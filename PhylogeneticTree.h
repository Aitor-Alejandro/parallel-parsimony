//Accelerating the Phylogenetic Parsimony Function on Heterogeneous Systems
//Authors: Sergio Santander-Jimenez, Aleksandar Ilic, Leonel Sousa, and Miguel A. Vega-Rodriguez
//Header of the Phylogenetic class, implementation of phylogenetic tree methods based on BIO++

#ifndef _PHYLOGENETICTREE_H_
#define _PHYLOGENETICTREE_H_
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include "NewickCode.h"
#include "MersenneTwister.h"
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Seq/Container/SiteContainer.h>   
#include <Bpp/Seq/Container/SiteContainerTools.h> 
#include <Bpp/Seq/Io/Phylip.h>  
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Phyl/Tree.h>          

using namespace std;
//BIO++ namespace
using namespace bpp;

class PhylogeneticTree
{
	int id; //Numerical identifier of the phylogenetic tree
	NewickCode code; //Newick code
	int parsimonyScore; //Parsimony score calculated by BIO++ (for verification purposes)
	
	DNA* alphabet; //BIO++ DNA alphabet instance
	SiteContainer* sites; //BIO++ dataset instance (with gaps)
	SiteContainer* completeSites;  //BIO++ dataset instance (without gaps - for future implementations i.e. likelihood)

	TreeTemplate<Node>* inferredTree; //BIO++ phylogenetic tree object
	int* terminal_nodes_seqpositions; //Array pointing to the position of the sequences corresponding to the terminal nodes of hte tree

	private:
		string optimizeParsimony(TreeTemplate<Node>*treeData); //Method that implements a basic topological optimization (for future implementations)
		int setScoresP (TreeTemplate<Node>*treeData); //Method that returns the parsimony score of the input tree using the BIO++ implementation of the function
	public:
		PhylogeneticTree (); //Default constructor
		PhylogeneticTree (int _id); //Parametrized constructor 1
		PhylogeneticTree (int _id, NewickCode _code); //Parametrized constructor 2
		
		void initialize (SiteContainer* _sites, SiteContainer* _completeSites); //Initialization of the BIO++ data structures DNA and dataset
		void initializeTree(); //Initialization of the TreeTemplate<Node>* object
		void initializeNewick (); //Initialization of the newick code object
		void deleteTree(); //Deleting the TreeTemplate<Node>* object
		int* getTerminalPositions(); //Returns the terminal_nodes_seqpositions array		
		TreeTemplate<Node>* getInferredTree(); //Returns the inferredTree variable
		void setInferredTree(TreeTemplate<Node>*treeData); //Sets a new phylogenetic tree
		void setInferredTreeNull(); //Sets to NULL the inferredTree object
								
		int getId (); //Returns the identifier of the phylogenetic tree
		NewickCode* getNewick (); //Returns the newick code of the phylogenetic tree
		int getParsimony();	//Returns the parsimony score stored in parsimonyScore
		void setParsimony(int _parsimonyScore); //Sets the parsimony score from the input parameters
		void setNewick (string _code); //Sets the newick code from the input parameter
		void setId(int _id); //Sets a new identifier for the phylogenetic tree object
		int setScores (); //Computes the parsimony score by calling the private method setScoresP
		
		int writeNewick (FILE* file); //Write the newick code to the specific file
		int readNewick (FILE* file); //Read the newick code from the specific file
		int readTree (FILE* file); //Read tree data from the specified file
		int writeTree(FILE* file); //Write tree data to the specified file
		
		void clone (const PhylogeneticTree & _tree); //Copy the contents of _tree into the object
		void printTree();


		~PhylogeneticTree();
};



#endif	//_PHYLOGENETICTREE_H_

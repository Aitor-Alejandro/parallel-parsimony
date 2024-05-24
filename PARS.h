//Header of the PARS class

#ifndef _PARS_H_
#define _PARS_H_
#include <iostream>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <time.h>
#include "TreeInterface.h"
//#include "TreeHeuristic.h"
#include "MersenneTwister.h"
#include <limits>
#include <fstream>
#include <omp.h>

using namespace std;

#define MAX_SONS 10 //Maximum number of children in the topology. It can be configured according to the topological characteristics of the evaluated trees

//Optimized topological node structure. This is the basic type for the data structure that will contain the post-order phylogenetic tree traversal for its processing by the kernel
typedef struct
{
	int partialParsimony = 0;
	int father = -1;//FATHER ID IN GLOBAL TREE, NEEDS CORRESPONDENCY ARRAY
	int id_node;
	char* characters;
	short number_of_sons; //Number of children of the node
	int sons_ids [MAX_SONS]; //Identifiers of the children node ID IF GLOBAL TREE
}typeNode;


class PARS
{	
	TreeInterface* refTree; //Reference tree (BIO++)
	//TreeHeuristic* treeH;

	string refFileName; //reference sequences filename
	SiteContainer* refSites; //BIO++ reference sequences
	int n_sequences; //number of input sequences in the reference 
	int n_sites; //input sequence length in the reference
	int totalParsimony;
	
	string queryFileName; //query sequences filename
	SiteContainer* querySites; //BIO++ query sequences
	int n_queries; //number of input sequences in the query file 
	//NOTA: Asumimos que las queries tienen la misma longitud de secuencia que las secuencias originales, así que nos vale n_sites tambien para query
	
	Phylip* seqReader; //BIO++ sequence reader (reference sequences)
	Phylip* seqReader2; //BIO++ sequence reader (query sequences)
	DNA* alphabet; //BIO++ DNA alphabet
	MTRand* mtr; //Random number generator

	typeNode* parsNodes;//, *copy_parsNodes; //Phylogenetic tree representation using our optimized topological node structure
	int n_parsNodes; //Number of nodes in the phylogeny
	//int* internal_nodes; //Identifiers of internal nodes
	int num_internal_nodes; //Number of internal nodes
	
	char** reference_sequences; //reference sequences in hexadecimal codification
	char* array_reference_sequences; //reference sequences in hexadecimal codification (char array)
	
	char** query_sequences; //query sequences in hexadecimal codification
	char* array_query_sequences; //query sequences in hexadecimal codification (char array)

	//int** matrixParsimony;
	private:
		
		void modifyVector(typeNode* internalNode, int father, int son);
		int calculateParsimonyQuerysPriv(int fatherNode, int son_replaced, typeNode* internalNode, typeNode* parsAux);
		int calculateParsimonyQuerysPriv(int fatherNode, int son_replaced, typeNode* internalNode, typeNode* parsAux, int local_n_sites, char** matrix_aux);
	protected: 
		double get_time(); //Get timestamp using omp_get_wtime
	public:
		void genInternalNode(typeNode* internalNode, char* query, char* characters);
		void genInternalNode(typeNode* internalNode, char* query, char* characters, int _local_n_sites);
		PARS(string _refFileName, string _queryFileName); //Constructor
		/* Initialization procedures */
		PhylogeneticTree* readTreeFromFile (FILE* file, int _id); //Reads the reference phylogenetic tree from file
		int genInitialTree (string fic_trees); //Initializes the phylogenetic tree object (refTree, BIO++)
		TreeInterface* getReferenceTree(); //Returns refTree				
		void initializeReferenceSequences (); //Initialize the reference sequences (hexadecimal code)
		void initializeQuerySequences (); //Initialize the query sequences (hexadecimal code)
		void initializeParsTree (); //Initialize the parsNodes structure (topology for kernel processing)
		void deleteAuxParsStructures(); //Deletes the parsNodes structure and initializes to 0 its related variables

		/* Parsimony calculations over the original reference tree */
		int calculateParsimonyRefTree (double &t1, double &t2); //Parsimony function code for CPU in char configuration

		/* Parsimony calculations for each Query*/
		int** calculateParsimonyQuerys (double &t1, double &t2);
		int** calculateParsimonyQuerysGrueso(double &t1, double &t2);
		//int** calculateParsimonyQuerysFino(double &t1, double &t2);
		//int** calculateParsimonyQuerysFino2(double &t1, double &t2);
		int calculateParsimonyQuerysPub(int fatherNode, int son_replaced, typeNode* internalNode, typeNode* parsAux);
		void cloneParsNodes(typeNode* dest1, typeNode* dest2);
		void cloneParsNodes(typeNode* dest);
		void clonePars(typeNode* dst, typeNode* src);
		void restaurar(typeNode*dst, const typeNode*src, int index_node);
		/*clone reference sequences*/
		void cloneRefSeq(char* dst, char* src);
		//void writeFileBruteForce(ofstream solveFile, typeNode* parsBest, int queryId, int position, int bestParsimony);
		/* Main method */
		int run (string fic_tree, double &t1, double &t2, int mode); //Main method, where the placement takes place
		
		typeNode* getParsNodes();
		int getNumInternalNodes();
		int getNumberOfSequences();
		char* getArrayReferenceSequences();
		char* getArrayQuerySequences();
		int get_n_sites();
		int getTotalParsimony();
		int get_n_querys();
		

		void setParsNodes(typeNode* newParsNodes);
		void set_num_internal_nodes(int new_num_internal_nodes);
		void set_n_sequences(int new_n_sequences);
		void set_seq(char* new_sequences);
		void set_n_sites(int n_sites);

		//void setParsNodes(typeNode* newParsNodes);
		~PARS(); //Destructor
		
	
};

#endif	//_PARS_H_




#include "PARS.h"

/**
*	Constructor, initializes variables and data structures. Parameters:
*	_refFileName: string containing the name/path of the reference sequences
*	_queryFileName: string containing the name/path of the query sequences 
*/
PARS::PARS(string _refFileName, string _queryFileName)
{
    refTree = NULL;	
	
	refFileName = _refFileName;
	refSites=NULL;
	n_sequences=0;
	n_sites=0;	
    
	queryFileName = _queryFileName;
	querySites=NULL;	
	n_queries=0;
	
	alphabet=NULL;
	seqReader=NULL;
	seqReader2=NULL;
	mtr = new MTRand ();
    
	reference_sequences=NULL;
	query_sequences=NULL;
	parsNodes=NULL;
}
		
/**
*	Method that returns a timestamp using time.h
*/
double PARS::get_time (void)
{
		return omp_get_wtime();
}


/**
*	Method that reads from a file the Newick code of the reference tree. Parameters:
*	file: reference tree file
*	_id: position of the tree to be read (in case there are more than 1, default = 0)
*/
PhylogeneticTree* PARS::readTreeFromFile (FILE* file, int _id)
{
	char character;
	int numtrees=0;
	NewickCode auxCode;
	//We suppose the file has been opened previously
	if (file == NULL) return NULL;
	if (_id < 0) return NULL;
	fseek(file, 0L, SEEK_SET);	
	character = fgetc(file);
	//Reading until we find the tree
	while (character != EOF && numtrees != _id)
	{
		if (character == ';')
		{ 
			numtrees++;			 
		}		
		character = fgetc(file);
	}
	if (numtrees == _id)
	{
		if (_id == 0) fseek(file, 0L, SEEK_SET);
		//Reading the tree
		if (!auxCode.readNewick (file)) return NULL;
		PhylogeneticTree* tree= new PhylogeneticTree (_id, auxCode); 		
		return tree;
	}
	else return NULL;
	
}

/**
*	Method that reads and initializes the BIO++ structures (TreeInterface, SeqReader... objects). Parameters:
*	fic_trees: name/path of the input reference tree file
*/
int PARS::genInitialTree (string fic_trees)
{
    FILE* file;
    PhylogeneticTree* tree;
	printf("Initializing reference tree via BIO++\n");
    //Opening tree file
    file = fopen (fic_trees.c_str(), "rt");
    if (file == NULL) return 0;

    //Reading the trees   
    tree = readTreeFromFile (file, 0); 
	if (tree==NULL) return 0;
   
	//Closing file
    fclose(file);
	
	//Initializing BIO++ objetcs
	alphabet = new DNA ();
	seqReader= new Phylip (true, false);
	refSites=seqReader->readAlignment(refFileName, alphabet);
	seqReader2 = new Phylip (true, false);
	querySites=seqReader2->readAlignment(queryFileName, alphabet);
	
	//Initializing TreeInterface object
	tree->initialize(refSites, NULL);
	tree->initializeTree();
	refTree=new TreeInterface (0,*tree);
	
	//Finishing
	delete(tree);
    printf("Initializing reference tree -- DONE\n");
	return 1;
} 

/**
*	Method that returns the reference tree object.
*/		
TreeInterface* PARS::getReferenceTree()
{
	return refTree;
}

/**
*	Method that initializes the reference sequences (hexadecimal code) in char configuration
*/
void PARS::initializeReferenceSequences ()
{
		int i,j;
		printf("Initializing reference sequences... ");
		//Getting number of reference sequences and length
        n_sequences = refSites->getNumberOfSequences ();
        n_sites = refSites->getNumberOfSites();
		//Initializing reference dataset variable
        reference_sequences = new char* [n_sequences];
        for (i=0; i<n_sequences; i++)
        {
                reference_sequences[i]=new char [n_sites];
        }
		//Converting the input dataset (BIO++ implementation) to our codification (reference_sequences)
        string aux_sequence;
        for (i=0; i<n_sequences; i++)
        {
                aux_sequence = refSites->toString(i);
                //cout<<"Sequence "<<i<<" "<<refSites->getName(i)<<endl;

                for(j=0; j<n_sites; j++)
                {
                        switch(aux_sequence[j])
                        {
								case '-': reference_sequences[i][j]=0x10; break;
                                case 'A': reference_sequences[i][j]=0x08; break;
                                case 'C': reference_sequences[i][j]=0x04; break;
                                case 'G': reference_sequences[i][j]=0x02; break;
                                case 'T': reference_sequences[i][j]=0x01; break;
                                case 'M': reference_sequences[i][j]=0xC; break;
                                case 'R': reference_sequences[i][j]=0xA; break;
                                case 'W': reference_sequences[i][j]=0x09; break;
                                case 'S': reference_sequences[i][j]=0x06; break;
                                case 'Y': reference_sequences[i][j]=0x05; break;
                                case 'K': reference_sequences[i][j]=0x03; break;
                                case 'V': reference_sequences[i][j]=0xE; break;
                                case 'H': reference_sequences[i][j]=0xD; break;
                                case 'D': reference_sequences[i][j]=0xB; break;
                                case 'B': reference_sequences[i][j]=0x07; break;

                                default: reference_sequences[i][j]=0xF; break;

			}
               //         printf("%x ", reference_sequences[i][j]);
                }
               // printf("\n");
        }
		

	    //array_reference_sequences: char container of the reference in row-major order
        array_reference_sequences = new char[n_sequences*n_sites]; 
		int k;
        k=0;
        for (i=0; i<n_sequences; i++)
        {
            for(j=0; j<n_sites; j++)
            {
                    array_reference_sequences[k]=reference_sequences[i][j];
                    k++;
            }
        }
		printf("OK\n");
	printf("A total of %d reference sequences, with length %d\n", n_sequences, n_sites);
}

/**
*	Method that initializes the query sequences (hexadecimal code) in char configuration
*/
void PARS::initializeQuerySequences ()
{
		int i,j;
		printf("Initializing query sequences... ");
		//Getting number of query sequences
        n_queries = querySites->getNumberOfSequences();
      
		//Initializing query dataset variable
        query_sequences = new char* [n_queries];
        for (i=0; i<n_queries; i++)
        {
                query_sequences[i]=new char [n_sites];
        }
		//Converting the input dataset (BIO++ implementation) to our codification (query_sequences)
        string aux_sequence;
        for (i=0; i<n_queries; i++)
        {
                aux_sequence = querySites->toString(i);
                //cout<<"Query "<<i<<" "<<querySites->getName(i)<<endl;

                for(j=0; j<n_sites; j++)
                {
                        switch(aux_sequence[j])
                        {
								case '-': query_sequences[i][j]=0x10; break;
                                case 'A': query_sequences[i][j]=0x08; break;
                                case 'C': query_sequences[i][j]=0x04; break;
                                case 'G': query_sequences[i][j]=0x02; break;
                                case 'T': query_sequences[i][j]=0x01; break;
                                case 'M': query_sequences[i][j]=0xC; break;
                                case 'R': query_sequences[i][j]=0xA; break;
                                case 'W': query_sequences[i][j]=0x09; break;
                                case 'S': query_sequences[i][j]=0x06; break;
                                case 'Y': query_sequences[i][j]=0x05; break;
                                case 'K': query_sequences[i][j]=0x03; break;
                                case 'V': query_sequences[i][j]=0xE; break;
                                case 'H': query_sequences[i][j]=0xD; break;
                                case 'D': query_sequences[i][j]=0xB; break;
                                case 'B': query_sequences[i][j]=0x07; break;

                                default: query_sequences[i][j]=0xF; break;

			}
                        //printf("%x ", query_sequences[i][j]);
                }
                //printf("\n");
        }
		
	    //array_query_sequences: char container of the queries in row-major order
        array_query_sequences = new char[n_queries*n_sites]; 
		int k;
        k=0;
        for (i=0; i<n_queries; i++)
        {
            for(j=0; j<n_sites; j++)
            {
                    array_query_sequences[k]=query_sequences[i][j];
                    k++;
            }
        }
		printf("OK\n");
	printf("A total of %d query sequences, with length %d\n", n_queries, n_sites);
}


/**
*	Method that initializes the parsNodes structure (topology for kernel processing). Post-order tree traversal.
*/
void PARS::initializeParsTree ()
{
	int i,j;
	int current_inner;
	int aux_nsons;
	int aux_son_id;
	string leaf_name;
	TreeTemplate<Node>*treeData; //BIO++ tree structure
	vector<Node *> nodes; //BIO++ node structure
	
	//Obtaining information about the phylogenetic topology
	treeData = refTree->getTree()->getInferredTree();
	nodes = treeData->getNodes();
	num_internal_nodes = treeData->getInnerNodesId().size();
	int inner_index [nodes.size()];
	parsNodes = new typeNode [num_internal_nodes];
	n_parsNodes=num_internal_nodes;
	 
	current_inner = 0;
	int* aux_term_pos = refTree->getTree()->getTerminalPositions();
	//TREE TRAVERSAL
	for(i=0; i<nodes.size(); i++)
	{
		aux_nsons = nodes[i]->getNumberOfSons();
		if (aux_nsons!=0)
		{	
			//INNER NODES
			parsNodes[current_inner].number_of_sons = aux_nsons;
			inner_index[i]=current_inner;
			for (j=0; j<aux_nsons; j++)
			{
				aux_son_id = nodes[i]->getSonsId()[j];
				if(nodes[aux_son_id]->getNumberOfSons()==0)
				{
					//TERMINAL CHILD                    
					parsNodes[current_inner].sons_ids[j]=aux_term_pos[aux_son_id];
                    parsNodes[current_inner].sons_ids[j]=parsNodes[current_inner].sons_ids[j]|(0<<31);
				}
				else
				{
					//INNER CHILD
					parsNodes[current_inner].sons_ids[j]=inner_index[aux_son_id];
					parsNodes[current_inner].sons_ids[j]=parsNodes[current_inner].sons_ids[j]|(1<<31);
				}
			}
			current_inner++;
		}
	}	
}

/**
*	Method that deletes the parsNodes structure and initializes to 0 its relates variable. 
*/
void PARS::deleteAuxParsStructures()
{
	delete(parsNodes);
	num_internal_nodes=0;
	n_parsNodes=0;
}

/**
*	Method that calculates the parsimony score for the reference tree. Parameters:
*	t1: variable to store the starting execution time
*	t2: variable to store the ending execution time
*/

int PARS::calculateParsimonyRefTree (double &t1, double &t2)
{
	int i,j,k,l;
		
	//COLUMN-MAJOR ORDER VARIABLE FOR THE CPU DATASET
    char* array_reference_sequences_omp;
    array_reference_sequences_omp = new char [n_sequences*n_sites];
	//Initializing the partial parsimony scores array 
	int accPars; //accumulated parsimony score
	char* my_characters; //character state values for an inner node
  	char* sequence_line; //sequence characters to be processed at the j-th iteration
	int num_sons; //number of children of the node currently processed
	int node_class; //type of node (leaf or internal)
	int node_id; //identifier of the node currently processed
	char aux_value; //auxiliar variable for fitch operations
	char site_value; //variable to store the state calculated for a node
	char son_value; //variable to store the state read from a child node

	//Evaluation loop
	t1 = get_time();
	
	//ORGANIZING REFERENCE SEQUENCES IN COLUMN-MAJOR ORDER	
	k=0;
    for (i=0; i<n_sites; i++)
    {
       	for(j=0; j<n_sequences; j++)
       	{
           	array_reference_sequences_omp[k]=array_reference_sequences[j*n_sites+i];
           	k++;
       	}
    }

	printf("Scoring reference tree\n");
	accPars=0;		
	my_characters = new char [num_internal_nodes];
	for (j=0; j<n_sites; j++)
	{
		sequence_line = &array_reference_sequences_omp[n_sequences*j];
		for (k=0; k<num_internal_nodes; k++)
		{
			site_value = 31;
			num_sons = parsNodes[k].number_of_sons;
			for (l=0; l<num_sons; l++)
			{
				node_class = parsNodes[k].sons_ids[l] & 0x80000000;
				node_id = parsNodes[k].sons_ids[l] & 0x7FFFFFFF;
				if (node_class == 0x80000000)
					son_value = my_characters[node_id];
				else
					son_value = sequence_line[node_id];
				aux_value = site_value & son_value;
				if (aux_value==0)
				{
					accPars++;
					aux_value = site_value | son_value;
				}
				site_value=aux_value;
			}
			my_characters[k] = site_value;
		}
	}
				
	printf("PARSIMONY SCORE %d\n", accPars);
	delete(my_characters);
	delete(array_reference_sequences_omp);
	//Showing execution times
	t2 = get_time();
	printf("Reference tree evaluation: Time spent: %f\n", t2-t1);
	return 0;
}



/**
*	Main method. Parameters:
*	fic_tree: name of the input file containing the reference tree
*	t1: variable to store the starting execution time
*	t2: variable to store the ending execution time
*/

int PARS::run (string fic_tree, double &t1, double &t2)
{
    int i;
	//Reading and initializing input tree
	
	genInitialTree(fic_tree);  	
	initializeReferenceSequences();
	initializeQuerySequences();
	initializeParsTree();
	calculateParsimonyRefTree (t1, t2);
	
	/***********************************
	************************************
	HERE WE CALL THE PLACEMENT ALGORITHM	
	************************************
	************************************/
	
	deleteAuxParsStructures();
	if(reference_sequences!=NULL)
	{
		for (i=0; i<n_sequences; i++)
				delete reference_sequences[i];
		delete (reference_sequences);
	}
	if(array_reference_sequences!=NULL)
		delete(array_reference_sequences);
	if(query_sequences!=NULL)
	{
		for (i=0; i<n_queries; i++)
				delete query_sequences[i];
		delete (query_sequences);
	}
	if(array_query_sequences!=NULL)
		delete(array_query_sequences);
	if (refTree != NULL)
	{  
		  delete(refTree);  
	}

	if (mtr!=NULL) delete (mtr);
	return 0;
}

/** Destructor */
PARS::~PARS()
{

 
}




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

	//solveFile = NULL;
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
	treeH = new TreeHeuristic(refTree);
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
	printf("-------------\n");
	printf("%d\n", num_internal_nodes);
	printf("-------------\n");
	n_parsNodes=num_internal_nodes;
	for (int i = 0; i < num_internal_nodes; i++){
		parsNodes[i].characters = new char [n_sites];
	}

	current_inner = 0;
	int* aux_term_pos = refTree->getTree()->getTerminalPositions();
	int* aux_correspondencia = new int [num_internal_nodes];
	//TREE TRAVERSAL
	Node* father = NULL;
	int id_father;
	for(i=0; i<nodes.size(); i++)
	{
		aux_nsons = nodes[i]->getNumberOfSons();
		//Node* father = nodes[i]->getFather();
		//if (father != NULL)
		//	int id_father = father->getId();
		if (aux_nsons!=0)
		{	
			//INNER NODES
			parsNodes[current_inner].number_of_sons = aux_nsons;
			inner_index[i]=current_inner;
			aux_correspondencia[current_inner] = i;

			father = nodes[i]->getFather();
			if (father != NULL) id_father = father->getId();
			parsNodes[current_inner].father = id_father;
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
	treeH->setCorrespondencias(aux_correspondencia);
	for (int o = 0; o < num_internal_nodes; o++)
		printf("Index: %d Value: %d\n",o,aux_correspondencia[o]);
	for (int k = 0; k < num_internal_nodes; k++){
		for (int l = 0; l < num_internal_nodes; l++){
			if (aux_correspondencia[k] == parsNodes[l].father){
				parsNodes[l].father = k;
			}
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
	typeNode* nodes_aux = NULL;
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
				if (node_class == 0x80000000)//CLASS INNER NODE
					son_value = my_characters[node_id];
				else//CLASS LEAF
					son_value = sequence_line[node_id];
				aux_value = site_value & son_value;
				if (aux_value==0)
				{
					parsNodes[k].partialParsimony++;
					accPars++;
					aux_value = site_value | son_value;
				}
				site_value=aux_value;
			}
			my_characters[k] = site_value;
			parsNodes[k].characters[j] = my_characters[k];
		}
		
	}
				
	printf("PARSIMONY SCORE %d\n", accPars);
	totalParsimony = accPars;
	delete(my_characters);
	delete(array_reference_sequences_omp);
	//Showing execution times
	t2 = get_time();
	printf("Reference tree evaluation: Time spent: %f\n", t2-t1);
	return 0;
}


int PARS::calculateParsimonyQuerys(double &t1, double &t2){
	int i,j,k,l,m;
	//COLUMN-MAJOR ORDER VARIABLE FOR THE CPU DATASET SEQUENCES
    char* array_reference_sequences_omp;
    array_reference_sequences_omp = new char [n_sequences*n_sites];
	//COLUMN-MAJOR ORDER VARIABLE FOR DE CPU DATASET QUERYS
	char* array_query_sequences_omp;
	array_query_sequences_omp = new char [n_queries*n_sites];

	//Initializing the partial parsimony scores array 
	int accPars; //accumulated parsimony score
	char* my_characters; //character state values for an inner node
  	char* sequence_line; //sequence characters to be processed at the j-th iteration
	char* query_line;
	int num_sons; //number of children of the node currently processed
	int node_class; //type of node (leaf or internal)
	int node_id; //identifier of the node currently processed
	char aux_value; //auxiliar variable for fitch operations
	char site_value; //variable to store the state calculated for a node
	char son_value; //variable to store the state read from a child node

	int auxParsimony = totalParsimony;
	int bestParsimony = -1;
	int position = -1;
	typeNode* auxNode = new typeNode;
	auxNode->characters = new char[n_sites];
	typeNode* parsAux = new typeNode[num_internal_nodes];
	typeNode* parsBest = new typeNode[num_internal_nodes];

	for (i = 0; i < num_internal_nodes; i++){
		parsAux[i].characters = new char[n_sites];
		parsBest[i].characters = new char[n_sites];
	}

	/*for (int i = 0; i < num_internal_nodes; i++){
		parsNodes[i].characters = new char [n_sites];
	}
*/

	cloneParsNodes(parsAux, parsBest);

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

	//ORGANIZING QUERY SEQUENCES IN COLUMN-MAJOR ORDER
	k=0;
	for (i=0; i<n_sites; i++)
	{
		for (j=0; j<n_queries; j++)
		{
			array_query_sequences_omp[k] = array_query_sequences[j*n_sites+i];
			k++;
		}
	}


	//MATRIX TO STORE PARSIMONY OF QUERYS IN EACH POSITION
	int** matrixParsimony;
	matrixParsimony=(int**)malloc(n_queries * sizeof(int*));
	for (int n=0; n<n_queries; n++){
		matrixParsimony[n]=(int*)malloc(num_internal_nodes * sizeof(int));
	}

	printf ("Nº de nodos en el albol:------%d\n",n_parsNodes);
	printf ("Nº de nodos internos arbol:---%d\n", num_internal_nodes);
	printf ("Nº de nodos hojas\n");
	//parsNodes[num_internal_nodes - 1].father = -1;

	//ofstream solveFile;
	//solveFile.open("resultados.txt");

	printf("Total Parsimony: ->->-> %d", totalParsimony);
	int auxP = 0;
	printf("parsNodes original\n");
	for (int w = 0; w < n_parsNodes; w++){
		printf("(Id: %d, parsimony: %d)-", w, parsNodes[w].partialParsimony);
		auxP += parsNodes[w].partialParsimony;
	}
	printf("Final parsimony: %d\n", auxP);

	for (i = 0; i < n_queries; i++){
		query_line = &array_query_sequences[i*n_sites];
		cloneParsNodes(parsAux);
		bestParsimony = -1;
		for (j = 0; j < num_internal_nodes; j++){
			
			k = j;
			//auxNode = &parsAux[j];
			auxParsimony = totalParsimony;
			cloneParsNodes(parsAux);
			while (k != -1){
				//printf("Nodo actual %d",k);
				auxNode = &parsAux[k];
				
				for (l = 0; l < n_sites; l++){
					aux_value = query_line[l] & auxNode->characters[l];
					if (aux_value == 0){
						aux_value = query_line[l] | auxNode->characters[l];
						parsAux[k].partialParsimony++;
						auxParsimony++;
						
					}
					parsAux[auxNode->father].characters[l] = aux_value;
				}
				if (k == auxNode->father)
					k=-1;
				else{
					//printf(" jumping index %d to %d \n", k, auxNode->father);
					k = auxNode->father;
					//printf(" jumped succesfully\n");
				}
			}
			
			if (auxParsimony != totalParsimony &&(bestParsimony > auxParsimony || bestParsimony == -1)){
				printf("bestParsimony %d\n", bestParsimony);
				printf(" auxParsimony %d\n", auxParsimony);
				clonePars(parsBest, parsAux);
				bestParsimony = auxParsimony;
				position = j;
			}
			
			
			//calculateParsimony
		}
		printf("\n------------------------------------------------------------\n");
		printf("QueryId: %d insertada en posicion: %d con parsimonia total de %d\n", i, position, bestParsimony);
		
		auxP = 0;
		for (int w = 0; w < n_parsNodes; w++){
			//printf("(Parsimony: %d)-",parsBest[w].partialParsimony);
			auxP += parsBest[w].partialParsimony;
		}
		printf("Total parsimony %d\n", auxP);
		printf("------------------------------------------------------------\n");
		//writeFileBruteForce(solveFile, parsBest, i, position, bestParsimony);
		
	}
	//solveFile.close();

	/*my_characters = new char [num_internal_nodes];
	for (j=0; j<n_sites; j++)
	{
		sequence_line = &array_reference_sequences_omp[n_sequences*j];
		query_line = &array_query_sequences_omp[n_queries*j];
		for (k=0; num_internal_nodes; k++)
		{
			site_value = 31;
			num_sons = parsNodes[k].number_of_sons;
			for (l=0; l<num_sons; l++)
			{
			
			}
		}
	}
	*/
	return 0;
}

void PARS::cloneParsNodes(typeNode* dest1, typeNode* dest2){
	//printf("Entrando en cloneParsNodes\n");
	int i, j ,k;
	for( i = 0; i < n_parsNodes; i++){
		dest1[i].partialParsimony = dest2[i].partialParsimony = parsNodes[i].partialParsimony;
		dest1[i].father = dest2[i].father = parsNodes[i].father;
		dest1[i].number_of_sons = dest2[i].number_of_sons = parsNodes[i].number_of_sons;
		
		for (j = 0; j < parsNodes[i].number_of_sons; j++){
			//printf("Node: %d Char: %d\n", i, j);
			dest1[i].sons_ids[j] = dest2[i].sons_ids[j] = parsNodes[i].sons_ids[j];
		}

		for (k = 0; k < n_sites; k++){
			//printf("Node: %d Char: %d/%d\n", i, k, n_sites);
			dest1[i].characters[k] = dest2[i].characters[k] = parsNodes[i].characters[k];
		}
		/*for (int i = 0; i < num_internal_nodes; i++){
		parsNodes[i].characters = new char [n_sites];
	}
*/
		//dest1[]
	}
}

void PARS::cloneParsNodes(typeNode* dest){
	for(int i = 0; i < n_parsNodes; i++){
		dest[i].partialParsimony = parsNodes[i].partialParsimony;
		dest[i].father = parsNodes[i].father;
		dest[i].number_of_sons = parsNodes[i].number_of_sons;

		for (int j = 0; j < parsNodes[i].number_of_sons; j++){
			dest[i].sons_ids[j] = parsNodes[i].sons_ids[j];
		}

		for (int j = 0; j < n_sites; j++){
			dest[i].characters[j] = parsNodes[i].characters[j];
		}
	}
}

void PARS::clonePars (typeNode* dst, typeNode* src){
	for(int i = 0; i < n_parsNodes; i++){
		dst[i].partialParsimony = src[i].partialParsimony;
		dst[i].father = src[i].father;
		dst[i].number_of_sons = src[i].number_of_sons;

		for (int j = 0; j < src[i].number_of_sons; j++){
			dst[i].sons_ids[j] = src[i].sons_ids[j];
		}

		for (int j = 0; j < n_sites; j++){
			dst[i].characters[j] = src[i].characters[j];
		}
	}
}
/*
void PARS::writeFileBruteForce(ofstream solveFile, typeNode* parsBest, int queryId, int position, int bestParsimony){
	solveFile << "Query id: " << queryId << " Position branch: " << position << " Parsimony value: " << bestParsimony << "\n";
	for (int i = 0; i < n_parsNodes; i++){
		solveFile << "Node index: " << i << " Partial parsimony: " << parsBest[i].partialParsimony << "-";
	}
	solveFile << "\n-----------------------------\n";
}
*/
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

	TreeHeuristic* treeH = new TreeHeuristic(refTree);
	treeH->show();
	calculateParsimonyQuerys(t1, t2);
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




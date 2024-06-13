#include "PARS.h"

#define NUM_BYTES 32
#define NUM_CORES 6
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
	//treeH = new TreeHeuristic(refTree);
	//Finishing
	//delete(tree);
	delete (tree);
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
								//case '-': reference_sequences[i][j]=0x10; break;
								case '-': reference_sequences[i][j]=0x0F; break;
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
								//case '-': query_sequences[i][j]=0x10; break;
								case '-': query_sequences[i][j]=0x0F; break;
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
	int current_inner, leave_index;
	int aux_nsons;
	int aux_son_id;
	string leaf_name;
	TreeTemplate<Node>*treeData; //BIO++ tree structure
	//TreeTemplateTools treetmptools;
	vector<Node *> nodes; //BIO++ node structure
	//Obtaining information about the phylogenetic topology
	treeData = refTree->getTree()->getInferredTree();
	nodes = treeData->getNodes();
	num_internal_nodes = treeData->getInnerNodesId().size();
	int inner_index [nodes.size()];
	parsNodes = new typeNode [num_internal_nodes];
	/*printf("-------------\n");
	printf("%d\n", num_internal_nodes);
	printf("-------------\n");*/
	n_parsNodes=num_internal_nodes;
	for (int i = 0; i < num_internal_nodes; i++){
		parsNodes[i].characters = new char [n_sites];
	}
	/*unsigned int depth;//int depth = treetmptools.getDepth(nodes[num_internal_nodes]);
	const  = nodes[47];
	depth = TreeTemplateTools::getDepth(*a);*/
	//depth = TreeTemplateTools::getDepth(nodes[num_internal_nodes]);
	Node* a;
	current_inner = 0;
	leave_index = 0;
	int* aux_term_pos = refTree->getTree()->getTerminalPositions();
	/*printf("LEAVES\n");
	for (int  i = 0; i < treeH->getNumLeaves(); i++){
		printf("%d-", aux_term_pos[i]);
	}*/
	//printf("\n");
	int* aux_correspondencia = new int [num_internal_nodes];
	//TREE TRAVERSAL
	Node* father = NULL;
	int id_father;
	//nodePhyl vNodePhyl[nodes.size()];
	for(i=0; i<nodes.size(); i++)
	{
		aux_nsons = nodes[i]->getNumberOfSons();
		if (aux_nsons!=0)
		{	
			//INNER NODES
			parsNodes[current_inner].number_of_sons = aux_nsons;
			inner_index[i]=current_inner;
			aux_correspondencia[current_inner] = i;
			parsNodes[current_inner].id_node = i;

			father = nodes[i]->getFather();
			if (father != NULL) id_father = father->getId();
			parsNodes[current_inner].father = id_father;
			a=nodes[i];
			parsNodes[current_inner].depth = TreeTemplateTools::getDepth(*a);
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
			//vNodePhyl[i].sub_index = current_inner;
			//vNodePhyl[i].sub_id = vNodePhyl[i].sub_index|(1<<31);
			current_inner++;
		}else{
			//vNodePhyl[i].sub_index = leave_index;
			//vNodePhyl[i].sub_id = vNodePhyl[i].sub_index|(0<<31);
			leave_index++;
		}
		
	}

	//---------------------------BORRAR---------------------------
	
	/*printf("---------------------PARS NODES RECIEN---------------------\n");
		for (int z = 0; z < num_internal_nodes; z++){
		printf("\nNODE: %d,   part_pars: %d\n", parsNodes[z].id_node, parsNodes[z].partialParsimony);
        printf("SONS:  ");
		for (int y = 0; y < parsNodes[z].number_of_sons; y++){
			if ((parsNodes[z].sons_ids[y] & 0x80000000) == 0x80000000)
                printf("INTE: ");
            else
            	printf("HOJA: ");
        	printf("-%u-", parsNodes[z].sons_ids[y] & 0x7FFFFFFF);
		}
		printf("\n");
	}*/
	//---------------------------BORRAR---------------------------
	//treeH->setCorrespondencias(aux_correspondencia);
	//treeH->setVectorCorresponcencias(vNodePhyl,nodes.size());
	//treeH->iniciarGrafo();
	//for (int o = 0; o < num_internal_nodes; o++)
	//	printf("Index: %d Value: %d\n",o,aux_correspondencia[o]);
	int depth;
	depth = parsNodes[num_internal_nodes-1].depth;
	printf("###############################################\n");
	printf("          ----------DEPTH=%i---------          \n",depth);
	printf("###############################################\n");
	/*for (int i=0;i < num_internal_nodes; i++)
		printf("-%d-", parsNodes[i].depth);
	printf("\n");*/
	for (int k = 0; k < num_internal_nodes; k++){
		parsNodes[k].id_node = treeData->getInnerNodesId()[k];
		for (int l = 0; l < num_internal_nodes; l++){
			if (aux_correspondencia[k] == parsNodes[l].father){
				parsNodes[l].father = k;
			}
		}
	}
	delete[] aux_correspondencia;
	//parsNodes[num_internal_nodes].father = -1;
}

/**
*	Method that deletes the parsNodes structure and initializes to 0 its relates variable. 
*/
void PARS::deleteAuxParsStructures()
{
	for (int i = 0; i < num_internal_nodes; i++){
		delete[] parsNodes[i].characters;
		//delete[] copy_parsNodes[i].characters;
	}
	//delete(parsNodes);
	delete[] parsNodes;
	parsNodes = NULL;
	//delete[] copy_parsNodes;
	//num_internal_nodes=0;
	//n_parsNodes=0;
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
	//totalParsimony = accPars;
	int aux = 0;
	for (int i = 0; i < num_internal_nodes; i++){
		aux += parsNodes[i].partialParsimony;
	}
	printf("PARSIMONY SCORE part%d\n", aux);
	//delete(my_characters);
	//delete(array_reference_sequences_omp);
	delete[] my_characters;
	delete[] array_reference_sequences_omp;
	//Showing execution times
	t2 = get_time();
	printf("Reference tree evaluation: Time spent: %f\n", t2-t1);
	/*for (int z = 0; z < num_internal_nodes; z++){
		printf("NODE: %d,   part_pars: %d\n", parsNodes[z].id_node, parsNodes[z].partialParsimony);
        printf("SONS:  ");
		for (int y = 0; y < parsNodes[z].number_of_sons; y++){
			if ((parsNodes[z].sons_ids[y] & 0x80000000) == 0x80000000)
                printf("INTE: ");
            else
            	printf("HOJA: ");
        	printf("-%u-", parsNodes[z].sons_ids[y] & 0x7FFFFFFF);
		}
		printf("\n");
	}*/
	return accPars;
}

void PARS::genInternalNode(typeNode* internalNode, char* query, char* characters){
	char aux_value;
	int i = 0;
	int num_th = omp_get_num_threads();
	internalNode->partialParsimony = 0;
	int _local_n_sites = n_sites;
	for (i = 0; i < _local_n_sites; i++){
		aux_value = query[i] & characters[i];
		if (aux_value == 0){
			internalNode->partialParsimony++;
			aux_value =  query[i] | characters[i];
		}
		internalNode->characters[i] = aux_value;
	}
}
void PARS::genInternalNode(typeNode* internalNode, char* query, char* characters, int _local_n_sites){
	char aux_value;
	int i = 0;
	internalNode->partialParsimony = 0;
	for (i = 0; i < _local_n_sites; i++){
		aux_value = query[i] & characters[i];
		if (aux_value == 0){
			internalNode->partialParsimony++;
			aux_value =  query[i] | characters[i];
		}
		internalNode->characters[i] = aux_value;
	}
}

void PARS::modifyVector(typeNode* internalNode, int father, int son){
	internalNode->father = father;
}

int PARS::calculateParsimonyQuerysPriv(int fatherNode, int son_replaced, typeNode* internalNode, typeNode* parsAux){
	parsAux[num_internal_nodes-1].father = -1;
	int node_class, node_id;
	char site_value, aux_value, son_value;
	int node;
	int local_n_sites = n_sites;
	node = fatherNode;
	while (node != - 1){
		parsAux[node].partialParsimony = 0;
		for (int i = 0; i < local_n_sites; i++){
			site_value = 31;
			aux_value = 0;
			for (int j = 0; j < parsAux[node].number_of_sons; j++){
				node_class = parsAux[node].sons_ids[j] & 0x80000000;
				node_id = parsAux[node].sons_ids[j] & 0x7FFFFFFF;
				if (node == fatherNode && j == son_replaced){
					son_value = internalNode->characters[i];
				}else{					
					if (node_class == 0x80000000){
						son_value = parsAux[node_id].characters[i];
					}else{
						son_value = array_reference_sequences[node_id*n_sites+i];
					}
				}
				aux_value = site_value & son_value;
				if (aux_value == 0){
					parsAux[node].partialParsimony++;
					aux_value = site_value|son_value;
				}
				site_value = aux_value;
			}
			parsAux[node].characters[i] = site_value;
		}
		node = parsAux[node].father;
	}
	int total_pars = 0;
	for (int l = 0; l < num_internal_nodes; l++){
		total_pars += parsAux[l].partialParsimony;
	}
	total_pars += internalNode->partialParsimony;
	return total_pars;
}

int PARS::calculateParsimonyQuerysPriv(int fatherNode, int son_replaced, typeNode* internalNode, typeNode* parsAux, int local_n_sites, char** matrix_aux){
	parsAux[num_internal_nodes-1].father = -1;
	int node_class, node_id;
	char site_value, aux_value, son_value;
	int node;
	node = fatherNode;
	while (node != - 1){
		parsAux[node].characters = matrix_aux[parsAux[node].depth - 1];
		parsAux[node].partialParsimony = 0;
		for (int i = 0; i < local_n_sites; i++){
			site_value = 31;
			aux_value = 0;
			for (int j = 0; j < parsAux[node].number_of_sons; j++){
				node_class = parsAux[node].sons_ids[j] & 0x80000000;
				node_id = parsAux[node].sons_ids[j] & 0x7FFFFFFF;
				if (node == fatherNode && j == son_replaced){
					son_value = internalNode->characters[i];
				}else{					
					if (node_class == 0x80000000){
						son_value = parsAux[node_id].characters[i];
					}else{
						son_value = array_reference_sequences[node_id*n_sites+i];
					}
				}
				aux_value = site_value & son_value;
				if (aux_value == 0){
					parsAux[node].partialParsimony++;
					aux_value = site_value|son_value;
				}
				site_value = aux_value;
			}
			parsAux[node].characters[i] = site_value;
		}
		node = parsAux[node].father;
	}
	int total_pars = 0;
	for (int l = 0; l < num_internal_nodes; l++){
		total_pars += parsAux[l].partialParsimony;
	}
	total_pars += internalNode->partialParsimony;
	return total_pars;
}

int PARS::calculateParsimonyQuerysPub(int fatherNode, int son_replaced, typeNode* internalNode, typeNode* parsAux){
	return calculateParsimonyQuerysPriv(fatherNode, son_replaced, internalNode, parsAux);
}

void restaurar1(typeNode* dst, const typeNode* src, int index_node, int num_internal, int n_s, char ** matrix_characters){
	int node = index_node;
	/*while (node != num_internal - 1){
		dst[node].id_node = src[node].id_node;
		dst[node].partialParsimony = src[node].partialParsimony;
		dst[node].father = src[node].father;
		dst[node].number_of_sons = src[node].number_of_sons;
		for (int i = 0; i < src[node].number_of_sons; i++){
			dst[node].sons_ids[i] = src[node].sons_ids[i];
		}
		for (int i = 0; i < n_s; i++){
			dst[node].characters[i] = src[node].characters[i];
		}
		node = src[node].father;
	}
	dst[node].id_node = src[node].id_node;
	dst[node].partialParsimony = src[node].partialParsimony;
	dst[node].father = src[node].father;
	dst[node].number_of_sons = src[node].number_of_sons;
	for (int i = 0; i < src[node].number_of_sons; i++){
		dst[node].sons_ids[i] = src[node].sons_ids[i];
	}
	for (int i = 0; i < n_s; i++){
		dst[node].characters[i] = src[node].characters[i];
	}*/
	while (node != num_internal-1){
		delete[] dst[node].characters;
		dst[node].characters = matrix_characters[node];
		node = src[node].father;
	}
	delete[] dst[node].characters;
	dst[node].characters = matrix_characters[node];
	//node = src[node].father;
}

int** PARS::calculateParsimonyQuerys(double &t1, double &t2){
	int i,j,k,l,m;

	//Initializing the partial parsimony scores array 
	//int accPars; //accumulated parsimony score
	char* query_line;
	//int num_sons; //number of children of the node currently processed
	int node_class; //type of node (leaf or internal)
	int node_id; //identifier of the node currently processed
	int node;
	int auxParsimony = totalParsimony;
	//int bestParsimony;
	int total_depth = parsNodes[num_internal_nodes-1].depth;
	typeNode* auxNode = new typeNode;
	auxNode->characters = new char[n_sites];
	int num_th = omp_get_num_threads();
	typeNode* parsAux = new typeNode[num_internal_nodes];
	char** matrix_characters = new char*[num_internal_nodes];
	char** matrix_aux = new char*[total_depth];
	for (i = 0; i < total_depth; i++)
		matrix_aux[i] = new char[n_sites];
	for (i = 0; i < num_internal_nodes; i++){
		//parsAux[i].characters = new char[n_sites];
		matrix_characters[i] = new char [n_sites];
		//matrix_aux[i] = new char[n_sites];
		parsAux[i].characters = matrix_characters[i];
	}

	cloneParsNodes(parsAux);
	//deleteAuxParsStructures();
	//MATRIX TO STORE BEST PARSIMONY OF QUERYS AND POSITION
	int** matrixParsimony;
	matrixParsimony=(int**)malloc(n_queries * sizeof(int*));
	for (int n=0; n<n_queries; n++){
		matrixParsimony[n]=(int*)malloc(3*sizeof(int));
		matrixParsimony[n][2] = INT_MAX;
	}
	int n_q = n_queries;
	int n_s = n_sites;
	int num_internal = num_internal_nodes;
	char *array_queries = array_query_sequences;
	char *array_references = array_reference_sequences;
	
	t1=get_time();
	for (i = 0; i < n_q; i++){
		query_line = &array_queries[i*n_s];
		for (j = 0; j < num_internal; j++){
			for (k = 0; k < parsAux[j].number_of_sons; k++){
				node_class = parsAux[j].sons_ids[k] & 0x80000000;
				node_id = parsAux[j].sons_ids[k] & 0x7FFFFFFF;
				if (node_class == 0x80000000){
					genInternalNode(auxNode, query_line, parsAux[node_id].characters, n_s);
					auxParsimony = calculateParsimonyQuerysPriv(j, k, auxNode, parsAux, n_s, matrix_aux);
				}else{
					genInternalNode(auxNode, query_line, &array_references[node_id*n_s], n_s);
					auxParsimony = calculateParsimonyQuerysPriv(j, k, auxNode, parsAux, n_s, matrix_aux);
				}
				if (auxParsimony < matrixParsimony[i][2]){
					matrixParsimony[i][0] = j;
					matrixParsimony[i][1] = k;
					matrixParsimony[i][2] = auxParsimony;
				}
				//restaurar1(parsAux, parsNodes, j, num_internal, n_s, matrix_characters);
				node = j;
				while (node != num_internal-1){
					//delete[] parsAux[node].characters;
					parsAux[node].partialParsimony = parsNodes[node].partialParsimony;
					parsAux[node].characters = matrix_characters[node];
					node = parsAux[node].father;
				}
				//delete[] parsAux[node].characters;
				parsAux[node].partialParsimony = parsNodes[node].partialParsimony;
				parsAux[node].characters =matrix_characters[node];
			}
		}
	}
	t2 = get_time();
	printf("\n");
	delete[] auxNode->characters;
	for (int l= 0; l < num_internal_nodes;l++){
		//delete[] parsAux[l].characters;
		//delete[] matrix_aux[l];
		delete[] matrix_characters[l];
		//delete[] matrix_aux[l];
	}
	for (int l = 0; l < total_depth; l++)
		delete[] matrix_aux[l];
	delete[] matrix_aux;
	delete[] matrix_characters;
	delete[] parsAux;
	delete auxNode;
	
	printf("-------TIME TO EVALUEATE PARSIMONY OF QUERIES (SEQ. VER.)-------\n");
	printf("Reference tree evaluation queries: Time spent: %f\n", t2-t1);
	printf("///////////////////////////////////////////////////////////\n");
	printf("\n");
	return matrixParsimony;
}

int** PARS::calculateParsimonyQuerysGrueso(double &t1, double &t2){
	int i,j,k;

	//Initializing the partial parsimony scores array 
	//int accPars; //accumulated parsimony score
	char* query_line;
	//int num_sons; //number of children of the node currently processed
	int node_class; //type of node (leaf or internal)
	int node_id; //identifier of the node currently processed
	//char aux_value; //auxiliar variable for fitch operations
	//char site_value; //variable to store the state calculated for a node
	//char son_value; //variable to store the state read from a child node

	int auxParsimony = totalParsimony;
	int total_depth = parsNodes[num_internal_nodes-1].depth;
	int** matrixParsimony;
	matrixParsimony=(int**)malloc(n_queries * sizeof(int*));
	for (int n=0; n<n_queries; n++){
		matrixParsimony[n]=(int*)malloc(3*sizeof(int));
		matrixParsimony[n][2] = INT_MAX;
	}
	//int NUM_CORES = sysconf(_SC_NPROCESSORS_ONLN);
	typeNode **parsAux_parallel= new typeNode*[NUM_CORES];
	typeNode *auxNode_parallel = new typeNode[NUM_CORES];
	char** matrix_characters = new char* [num_internal_nodes];
	char*** matrix_aux = new char**[NUM_CORES];
	for (int m = 0; m < num_internal_nodes; m++){
		//parsAux_parallel[l][m].characters = &matrix_characters[l][m*n_sites];
		matrix_characters[m] = new char [n_sites];
	}
	for (int l = 0; l < NUM_CORES; l++){
		//matrix_characters[l] = new char [num_internal_nodes*n_sites];
		parsAux_parallel[l] = new typeNode[num_internal_nodes];
		auxNode_parallel[l].characters = new char[n_sites];
		matrix_aux[l] = new char*[total_depth];
		for (int n = 0; n < total_depth; n++)
			matrix_aux[l][n] = new char[n_sites];
		for (int m = 0; m < num_internal_nodes; m++){
			parsAux_parallel[l][m].characters = matrix_characters[m];
		}
		cloneParsNodes(parsAux_parallel[l]);
	}
	//deleteAuxParsStructures();
	omp_set_num_threads(NUM_CORES);
	int th_id;
	int node;
	int n_q = n_queries;
	int n_s = n_sites;
	int num_internal = num_internal_nodes;
	int n_sons;
	char ** queries= query_sequences;
	char ** references = reference_sequences;
	t1 = get_time();
	//#pragma omp parallel for schedule (guided) default (none)private(th_id, query_line, auxParsimony, i, j, k, node_class, node_id, bestParsimony) shared (matrixParsimony, n_sites, array_reference_sequences, array_query_sequences, parsAux_parallel, auxNode_parallel, num_internal_nodes, n_queries)
	#pragma omp parallel /*for schedule (guided)*/ default (none)private(node,th_id, query_line, auxParsimony, i, j, k, n_sons, node_class, node_id) shared (matrixParsimony, n_s, references, queries, parsAux_parallel, auxNode_parallel, num_internal, n_q, matrix_characters, matrix_aux)
	{
	//#pragma omp parallel default (none) private(th_id, query_line, auxParsimony, i, j, k, node_class, node_id, bestParsimony) shared (matrixParsimony, n_s, array_references, array_queries, parsAux_parallel, auxNode_parallel, num_internal, n_q, n_seq, t1)
		th_id = omp_get_thread_num();
		#pragma omp for schedule (guided)
		for (i = 0; i < n_q; i++){
			query_line = queries[i];
			for (j = 0; j < num_internal; j++){
				n_sons = parsAux_parallel[th_id][j].number_of_sons;
				for (k = 0; k < n_sons; k++){
					node_class = parsAux_parallel[th_id][j].sons_ids[k] & 0x80000000;
					node_id = parsAux_parallel[th_id][j].sons_ids[k] & 0x7FFFFFFF;
					if (node_class == 0x80000000){
						//genInternalNode(auxNode, query_line, parsAux[node_id].characters);
						//auxParsimony = calculateParsimonyQuerysPriv(j, k, auxNode, parsAux);
						genInternalNode(&auxNode_parallel[th_id], query_line, parsAux_parallel[th_id][node_id].characters, n_s);
						auxParsimony = calculateParsimonyQuerysPriv(j,k,&auxNode_parallel[th_id], parsAux_parallel[th_id], n_s, matrix_aux[th_id]);
					}else{
						//genInternalNode(auxNode, query_line, &array_reference_sequences[node_id*n_sites]);
						//auxParsimony = calculateParsimonyQuerysPriv(j, k, auxNode, parsAux);
						genInternalNode(&auxNode_parallel[th_id], query_line, references[node_id], n_s);//array_reference_sequences[node_id*n_sites]);
						auxParsimony = calculateParsimonyQuerysPriv(j,k,&auxNode_parallel[th_id], parsAux_parallel[th_id], n_s, matrix_aux[th_id]);
					}
					if (auxParsimony < matrixParsimony[i][2]){
						matrixParsimony[i][0] = j;
						matrixParsimony[i][1] = k;
						matrixParsimony[i][2] = auxParsimony;
					}
					node = j;
					while (node != num_internal -1){
						parsAux_parallel[th_id][node].partialParsimony = parsNodes[node].partialParsimony;
						parsAux_parallel[th_id][node].characters = matrix_characters[node];
						node = parsAux_parallel[th_id][node].father;
					}
					parsAux_parallel[th_id][node].partialParsimony = parsNodes[node].partialParsimony;
					parsAux_parallel[th_id][node].characters =matrix_characters[node];
				}
			}
		}
	}
	t2 = get_time();
	for (int l = 0; l < NUM_CORES; l++){
		delete[] auxNode_parallel[l].characters;
		delete[] parsAux_parallel[l];
		for (int m = 0; m < total_depth; m++){
			delete[] matrix_aux[l][m];
		}
		delete[] matrix_aux[l];
	}
	for (int m = 0; m < num_internal_nodes; m++)
		delete[] matrix_characters[m];
	delete[] matrix_aux;
	delete[] auxNode_parallel;
	delete[] parsAux_parallel;
	delete[] matrix_characters;
	printf("\n");
	printf("-------TIME TO EVALUEATE PARSIMONY OF QUERIES (PARG. VER.)-------\n");
	printf("Reference tree evaluation queries: Time spent: %f\n", t2-t1);
	printf("///////////////////////////////////////////////////////////\n");
	printf("\n");
	return matrixParsimony;
}
int** PARS::calculateParsimonyQuerysFino(double &t1, double &t2){
	int i,j,k;

	//Initializing the partial parsimony scores array 
	//int accPars; //accumulated parsimony score
	char* query_line;
	//int num_sons; //number of children of the node currently processed
	int node_class; //type of node (leaf or internal)
	int node_id; //identifier of the node currently processed
	//char aux_value; //auxiliar variable for fitch operations
	//char site_value; //variable to store the state calculated for a node
	//char son_value; //variable to store the state read from a child node

	int auxParsimony;
	//int bestParsimony;
	int total_depth = parsNodes[num_internal_nodes-1].depth;
	int** matrixParsimony;
	matrixParsimony=(int**)malloc(n_queries * sizeof(int*));
	for (int n=0; n<n_queries; n++){
		matrixParsimony[n]=(int*)malloc(3*sizeof(int));
		matrixParsimony[n][2] = INT_MAX;
	}
	//int NUM_CORES = sysconf(_SC_NPROCESSORS_ONLN);
	typeNode **parsAux_parallel= new typeNode*[NUM_CORES];
	typeNode *auxNode_parallel = new typeNode[NUM_CORES];
	char** matrix_characters = new char* [num_internal_nodes];
	char*** matrix_aux = new char**[NUM_CORES];
	for (int m = 0; m < num_internal_nodes; m++){
		//parsAux_parallel[l][m].characters = &matrix_characters[l][m*n_sites];
		matrix_characters[m] = new char [n_sites];
	}
	for (int l = 0; l < NUM_CORES; l++){
		//matrix_characters[l] = new char [num_internal_nodes*n_sites];
		parsAux_parallel[l] = new typeNode[num_internal_nodes];
		auxNode_parallel[l].characters = new char[n_sites];
		matrix_aux[l] = new char*[total_depth];
		for (int n = 0; n < total_depth; n++)
			matrix_aux[l][n] = new char[n_sites];
		for (int m = 0; m < num_internal_nodes; m++){
			//matrix_aux[l][m] = new char[n_sites];
			parsAux_parallel[l][m].characters = matrix_characters[m];//[l][m*n_sites];
		}
		cloneParsNodes(parsAux_parallel[l]);
	}
	omp_set_num_threads(NUM_CORES);
	int th_id;
	int n_q = n_queries;
	int n_s = n_sites;
	int n_seq = n_sequences;
	int num_internal = num_internal_nodes;
	char *array_queries = array_query_sequences;
	char *array_references = array_reference_sequences;
	int node;
	t1 = get_time();
	for (i = 0; i < n_q; i++){
		query_line = &array_queries[i*n_s];
		#pragma omp parallel default (none) private(node,th_id, auxParsimony, j, k, node_class, node_id) shared (i, query_line, matrixParsimony, n_s, n_q, num_internal, array_references, array_queries, auxNode_parallel, parsAux_parallel, matrix_characters, matrix_aux)
		{
		th_id = omp_get_thread_num();
		#pragma omp for schedule (guided)
		for (j = 0; j < num_internal; j++){
			for (k = 0; k < parsAux_parallel[th_id][j].number_of_sons; k++){
				node_class = parsAux_parallel[th_id][j].sons_ids[k] & 0x80000000;
				node_id = parsAux_parallel[th_id][j].sons_ids[k] & 0x7FFFFFFF;
				if (node_class == 0x80000000){
					genInternalNode(&auxNode_parallel[th_id], query_line, parsAux_parallel[th_id][node_id].characters, n_s);
					auxParsimony = calculateParsimonyQuerysPriv(j,k,&auxNode_parallel[th_id], parsAux_parallel[th_id], n_s, matrix_aux[th_id]);
				}else{
					genInternalNode(&auxNode_parallel[th_id], query_line, &array_references[node_id*n_s], n_s);
					auxParsimony = calculateParsimonyQuerysPriv(j,k,&auxNode_parallel[th_id], parsAux_parallel[th_id], n_s, matrix_aux[th_id]);
				}
				if (auxParsimony < matrixParsimony[i][2]){
					matrixParsimony[i][0] = j;
					matrixParsimony[i][1] = k;
					matrixParsimony[i][2] = auxParsimony;
				}
				//restaurar1(parsAux_parallel[th_id], parsNodes, j, num_internal, n_s);
				node = j;
				while (node != num_internal -1){
					//delete[] parsAux_parallel[th_id][node].characters;
					parsAux_parallel[th_id][node].partialParsimony = parsNodes[node].partialParsimony;
					parsAux_parallel[th_id][node].characters = matrix_characters[node];
					node = parsAux_parallel[th_id][node].father;
				}
				//delete[] parsAux_parallel[th_id][node].characters;
				parsAux_parallel[th_id][node].partialParsimony = parsNodes[node].partialParsimony;
				parsAux_parallel[th_id][node].characters =matrix_characters[node];
			}
		}
		}
	}
	t2 = get_time();
	for (int l = 0; l < NUM_CORES; l++){
		delete[] auxNode_parallel[l].characters;
		delete[] parsAux_parallel[l];
		delete[] matrix_characters[l];
	}
	delete[] auxNode_parallel;
	delete[] parsAux_parallel;
	delete[] matrix_characters;
	printf("\n");
	printf("-------TIME TO EVALUEATE PARSIMONY OF QUERIES (PARF. VER.)-------\n");
	printf("Reference tree evaluation queries: Time spent: %f\n", t2-t1);
	printf("///////////////////////////////////////////////////////////\n");
	printf("\n");
	return matrixParsimony;
}

int** PARS::calculateParsimonyQuerysFino2(double &t1, double &t2){
	int i,j,k;

	//Initializing the partial parsimony scores array 
	//int accPars; //accumulated parsimony score
	char* query_line;
	//int num_sons; //number of children of the node currently processed
	int node_class; //type of node (leaf or internal)
	int node_id; //identifier of the node currently processed
	//char aux_value; //auxiliar variable for fitch operations
	//char site_value; //variable to store the state calculated for a node
	//char son_value; //variable to store the state read from a child node

	int auxParsimony;
	//int bestParsimony;

	//MATRIX TO STORE BEST PARSIMONY OF QUERYS AND POSITION
	int** matrixParsimony;
	int total_depth = parsNodes[num_internal_nodes-1].depth;
	matrixParsimony=(int**)malloc(n_queries * sizeof(int*));
	for (int n=0; n<n_queries; n++){
		matrixParsimony[n]=(int*)malloc(3*sizeof(int));
		matrixParsimony[n][2] = INT_MAX;
	}
	//int NUM_CORES = sysconf(_SC_NPROCESSORS_ONLN);
	typeNode **parsAux_parallel= new typeNode*[NUM_CORES];
	typeNode *auxNode_parallel = new typeNode[NUM_CORES];
	char** matrix_characters = new char* [num_internal_nodes];
	char*** matrix_aux = new char**[NUM_CORES];
	for (int m = 0; m < num_internal_nodes; m++){
		//parsAux_parallel[l][m].characters = &matrix_characters[l][m*n_sites];
		matrix_characters[m] = new char [n_sites];
	}
	for (int l = 0; l < NUM_CORES; l++){
		//matrix_characters[l] = new char [num_internal_nodes*n_sites];
		parsAux_parallel[l] = new typeNode[num_internal_nodes];
		auxNode_parallel[l].characters = new char[n_sites];
		matrix_aux[l] = new char*[total_depth];
		for (int n = 0; n < total_depth; n++)
			matrix_aux[l][n] = new char[n_sites];
		for (int m = 0; m < num_internal_nodes; m++){
			//matrix_aux[l][m] = new char[n_sites];
			parsAux_parallel[l][m].characters = matrix_characters[m];//[l][m*n_sites];
		}
		cloneParsNodes(parsAux_parallel[l]);
	}
	omp_set_num_threads(NUM_CORES);
	int th_id;
	int n_q = n_queries;
	int n_s = n_sites;
	//int n_seq = n_sequences;
	int num_internal = num_internal_nodes;
	char **queries = query_sequences;
	char **references = reference_sequences;
	char *auxChar;
	int n_sons;
	int node;
	t1 = get_time();
	//th_id = 0;
	#pragma omp parallel /*for schedule (guided)*/ default (none) private (node, th_id, n_sons, query_line, i, j, k, node_class, node_id, auxParsimony, auxChar) shared (matrixParsimony, n_q, n_s, num_internal, queries, references, parsAux_parallel, auxNode_parallel, matrix_aux, matrix_characters)
	{
		th_id = omp_get_thread_num();
	#pragma omp for schedule (guided)
	for (j = 0; j < num_internal; j++){
		n_sons = parsAux_parallel[th_id][j].number_of_sons;
		for (i=0; i < n_q; i++){
			query_line=queries[i];
			for (k = 0; k < n_sons; k++){
				node_class = parsAux_parallel[th_id][j].sons_ids[k] & 0x80000000;
				node_id = parsAux_parallel[th_id][j].sons_ids[k] & 0x7FFFFFFF;
				if (node_class == 0x80000000){
					//genInternalNode(&auxNode_parallel[th_id], query_line, parsAux_parallel[th_id][node_id].characters, n_s);
					//auxParsimony = calculateParsimonyQuerysPriv(j,k,&auxNode_parallel[th_id], parsAux_parallel[th_id], n_s);
					//auxChar = parsAux_parallel[th_id][node_id].characters;
					genInternalNode(&auxNode_parallel[th_id], query_line, parsAux_parallel[th_id][node_id].characters, n_s);
					auxParsimony = calculateParsimonyQuerysPriv(j,k,&auxNode_parallel[th_id], parsAux_parallel[th_id], n_s, matrix_aux[th_id]);
				}else{
					//genInternalNode(&auxNode_parallel[th_id], query_line, &array_references[node_id*n_sites], n_s);
					//auxParsimony = calculateParsimonyQuerysPriv(j,k,&auxNode_parallel[th_id], parsAux_parallel[th_id], n_s);
					//auxChar = &array_references[node_id*n_s];
					genInternalNode(&auxNode_parallel[th_id], query_line, references[node_id], n_s);//array_reference_sequences[node_id*n_sites]);
					auxParsimony = calculateParsimonyQuerysPriv(j,k,&auxNode_parallel[th_id], parsAux_parallel[th_id], n_s, matrix_aux[th_id]);
				}
				
					//genInternalNode(&auxNode_parallel[th_id], query_line, auxChar);
					//auxParsimony = calculateParsimonyQuerysPriv(j,k,&auxNode_parallel[th_id], parsAux_parallel[th_id], n_s);
					if (auxParsimony < matrixParsimony[i][2]){
						matrixParsimony[i][0] = j;
						matrixParsimony[i][1] = k;
						matrixParsimony[i][2] = auxParsimony;
					}
					//restaurar1(parsAux_parallel[th_id], parsNodes, j, num_internal, n_s);
					node = j;
					while (node != num_internal -1){
						//delete[] parsAux_parallel[th_id][node].characters;
						parsAux_parallel[th_id][node].partialParsimony = parsNodes[node].partialParsimony;
						parsAux_parallel[th_id][node].characters = matrix_characters[node];
						node = parsAux_parallel[th_id][node].father;
					}
					//delete[] parsAux_parallel[th_id][node].characters;
					parsAux_parallel[th_id][node].partialParsimony = parsNodes[node].partialParsimony;
					parsAux_parallel[th_id][node].characters =matrix_characters[node];
			}
		}
	}
	}
	t2 = get_time();
	for (int l = 0; l < NUM_CORES; l++){
		delete[] auxNode_parallel[l].characters;
		delete[] parsAux_parallel[l];
		for (int m = 0; m < total_depth; m++){
			delete[] matrix_aux[l][m];
		}
		delete[] matrix_aux[l];
	}
	for (int m = 0; m < num_internal_nodes; m++)
		delete[] matrix_characters[m];
	delete[] matrix_aux;
	delete[] auxNode_parallel;
	delete[] parsAux_parallel;
	delete[] matrix_characters;
	printf("\n");
	printf("-------TIME TO EVALUEATE PARSIMONY OF QUERIES (PARF. VER.)-------\n");
	printf("Reference tree evaluation queries: Time spent: %f\n", t2-t1);
	printf("///////////////////////////////////////////////////////////\n");
	printf("\n");
	return matrixParsimony;
}

void PARS::genInternalNodeSIMD(typeNode* internalNode, char* query, char* characters, int n_sites_vectorizado){
	int i, m;
	const int v = NUM_BYTES;
	char aux_value[v];
	int part_pars = 0;
	for (i = 0; i < n_sites_vectorizado; i += v){
		#pragma omp simd simdlen(v)
		for (m = 0; m < v; m++)
			aux_value[m] = query[i+m] & characters[i+m];
		#pragma omp simd simdlen(v) reduction(+:part_pars)
		for (m = 0; m < v; m++){
			if (aux_value[m] == 0)
			{
				part_pars++;
				aux_value[m] = query[i+m] | characters[i+m];
			}
			internalNode->characters[i+m] = aux_value[m];
		}
	}
	internalNode->partialParsimony = part_pars;
}

int PARS::calculateParsimonyQuerysPrivSIMD(int fatherNode, int son_replaced, typeNode* internalNode, typeNode* parsAux, int n_sites_vectorizado, char** matrix_aux, char** references){
	int i,j,k,m;
	const int v = NUM_BYTES;
	parsAux[num_internal_nodes-1].father = -1;
	int node_class, node_id;
	char site_value[v], aux_value[v], son_value[v];
	//printf("Node: %d:",fatherNode);
	int node;
	int part_pars = 0;
	node = fatherNode;
	while (node != -1){
		//printf("Node: %d:",node);
		parsAux[node].characters = matrix_aux[parsAux[node].depth - 1];
		parsAux[node].partialParsimony = 0;
		part_pars = 0;
		for (i = 0; i < n_sites_vectorizado; i+=v){
			#pragma omp simd simdlen(v)
			for (m = 0; m < v; m++){
				site_value[m] = 31;
				aux_value[m] = 0;
			}
			for (j = 0; j < parsAux[node].number_of_sons; j++){
				node_class = parsAux[node].sons_ids[j] & 0x80000000;
				node_id = parsAux[node].sons_ids[j] & 0x7FFFFFFF;
				if (node == fatherNode && j == son_replaced){
					#pragma omp simd simdlen(v)
					for (m = 0; m < v; m++)
						son_value[m] = internalNode->characters[i+m];
				}else{
					if (node_class == 0x80000000){
						#pragma omp simd simdlen(v)
						for (m = 0; m < v; m++)
							son_value[m] = parsAux[node_id].characters[i+m];
					}else{
						#pragma omp simd simdlen(v)
						for (m = 0; m < v; m++)
							son_value[m] = references[node_id][i+m];
					}
				}
				#pragma omp simd simdlen(v)
				for (m = 0; m < v; m++)
					aux_value[m] = site_value[m] & son_value[m];
				#pragma omp simd simdlen(v) reduction(+:part_pars)
				for (m = 0; m < v; m++){
					if (aux_value[m] == 0){
						part_pars++;
						aux_value[m] = site_value[m] | son_value[m];
					}
					site_value[m] = aux_value[m];
				}
			}
			#pragma omp simd simdlen(v)
			for (m = 0; m < v; m++)
				parsAux[node].characters[i+m] = site_value[m];
		}
		parsAux[node].partialParsimony += part_pars;
		node = parsAux[node].father;
	}
	int total_pars = 0;
	for (k = 0; k < num_internal_nodes; k++)
	{
		total_pars += parsAux[k].partialParsimony;
	}
	total_pars += internalNode->partialParsimony;
	return total_pars;
}

int** PARS::calculateParsimonyQuerysSIMD(double &t1, double &t2){
	int i,j,k,l;
	char* query_line;
	//int num_sons; //number of children of the node currently processed
	int node_class; //type of node (leaf or internal)
	int node_id; //identifier of the node currently processed
	//char aux_value; //auxiliar variable for fitch operations
	//char site_value; //variable to store the state calculated for a node
	//char son_value; //variable to store the state read from a child node

	int auxParsimony = totalParsimony;
	int total_depth = parsNodes[num_internal_nodes-1].depth;
	int node;

	const int v = NUM_BYTES;
	int n_sites_vectorizado = std::ceil((double)n_sites/v);
	n_sites_vectorizado = n_sites_vectorizado*v;

	typeNode* auxNode = new typeNode;
	auxNode->characters = new char[n_sites_vectorizado];
	int num_th = omp_get_num_threads();
	typeNode* parsAux = new typeNode[num_internal_nodes];
	char** matrix_characters = new char*[num_internal_nodes];
	char** matrix_aux = new char*[total_depth];
	for (i = 0; i < total_depth; i++)
		matrix_aux[i] = new char[n_sites_vectorizado];
	for (i = 0; i < num_internal_nodes; i++){
		matrix_characters[i] = new char [n_sites_vectorizado];
		parsAux[i].characters = matrix_characters[i];
	}

	cloneParsNodes(parsAux);
	
	int** matrixParsimony;
	matrixParsimony=(int**)malloc(n_queries * sizeof(int*));
	for (int n=0; n<n_queries; n++){
		matrixParsimony[n]=(int*)malloc(3*sizeof(int));
		matrixParsimony[n][2] = INT_MAX;
	}
	
	char **queries = new char* [n_queries];
	for (i = 0; i < n_queries; i++){
		queries[i] = new char[n_sites_vectorizado];
		for (j = 0; j < n_sites; j++)
			queries[i][j] = query_sequences[i][j];
		for (j = n_sites; j < n_sites_vectorizado; j++)
			queries[i][j] = 0x10;
	}
	
	char **references = new char* [n_sequences];
	for (i = 0; i< n_sequences; i++){
		references[i] = new char[n_sites_vectorizado];
		for (j = 0; j < n_sites; j++)
			references[i][j] = reference_sequences[i][j];
		for (j = n_sites; j < n_sites_vectorizado; j++)
			references[i][j] = 0x10;
	}

	for (i = 0; i < num_internal_nodes; i++){
		for (j = n_sites; j < n_sites_vectorizado; j++)
			matrix_characters[i][j] = 0x10;
	}
	int n_q = n_queries;
	int n_s = n_sites;
	int num_internal = num_internal_nodes;
	t1 = get_time();
	for (i = 0; i < n_q; i++){
		//printf("%d\n", i);
		query_line = queries[i];
		for (j = 0; j < num_internal; j++){
			for (k = 0; k < parsAux[j].number_of_sons;k++){
				node_class = parsAux[j].sons_ids[k] & 0x80000000;
				node_id = parsAux[j].sons_ids[k] & 0x7FFFFFFF;
				if (node_class == 0x80000000){
					genInternalNodeSIMD(auxNode, query_line, parsAux[node_id].characters, n_sites_vectorizado);
					//printf("generado\n");
					auxParsimony = calculateParsimonyQuerysPrivSIMD(j, k, auxNode, parsAux, n_sites_vectorizado, matrix_aux, references);
					
				}else{
					genInternalNodeSIMD(auxNode, query_line, references[node_id], n_sites_vectorizado);
					//printf("generado\n");
					auxParsimony = calculateParsimonyQuerysPrivSIMD(j, k, auxNode, parsAux, n_sites_vectorizado, matrix_aux, references);

					
				}
				if (auxParsimony < matrixParsimony[i][2]){
					matrixParsimony[i][0] = j;
					matrixParsimony[i][1] = k;
					matrixParsimony[i][2] = auxParsimony;
				}
				node = j;
				while (node != num_internal -1){
					//delete[] parsAux_parallel[th_id][node].characters;
					parsAux[node].partialParsimony = parsNodes[node].partialParsimony;
					parsAux[node].characters = matrix_characters[node];
					node = parsAux[node].father;
				}
				//delete[] parsAux_parallel[th_id][node].characters;
				parsAux[node].partialParsimony = parsNodes[node].partialParsimony;
				parsAux[node].characters =matrix_characters[node];
			}
		}
	}
	t2 = get_time();
	printf("\n");
	delete[] auxNode->characters;
	for (int l= 0; l < num_internal_nodes;l++){
		delete[] matrix_characters[l];
	}
	for (int l = 0; l < total_depth; l++)
		delete[] matrix_aux[l];
	for (int l = 0; l < n_queries; l++)
		delete[] queries[l];
	for (int l = 0; l < n_sequences; l++)
		delete[] references[l];
	delete[] matrix_aux;
	delete[] matrix_characters;
	delete[] parsAux;
	delete[] queries;
	delete[] references;
	delete auxNode;
	
	printf("-------TIME TO EVALUEATE PARSIMONY OF QUERIES (SEQ. VER.)-------\n");
	printf("Reference tree evaluation queries: Time spent: %f\n", t2-t1);
	printf("///////////////////////////////////////////////////////////\n");
	return matrixParsimony;
}

void PARS::cloneParsNodes(typeNode* dest1, typeNode* dest2){
	//printf("Entrando en cloneParsNodes\n");
	int i, j ,k;
	for( i = 0; i < n_parsNodes; i++){
		dest1[i].id_node = dest2[i].id_node = parsNodes[i].id_node;
		dest1[i].partialParsimony = dest2[i].partialParsimony = parsNodes[i].partialParsimony;
		dest1[i].father = dest2[i].father = parsNodes[i].father;
		dest1[i].number_of_sons = dest2[i].number_of_sons = parsNodes[i].number_of_sons;
		dest1[i].depth = dest2[i].depth = parsNodes[i].depth;
		for (j = 0; j < parsNodes[i].number_of_sons; j++){
			//printf("Node: %d Char: %d\n", i, j);
			dest1[i].sons_ids[j] = dest2[i].sons_ids[j] = parsNodes[i].sons_ids[j];
		}

		for (k = 0; k < n_sites; k++){
			//printf("Node: %d Char: %d/%d\n", i, k, n_sites);
			dest1[i].characters[k] = dest2[i].characters[k] = parsNodes[i].characters[k];
		}
	}
}

void PARS::cloneParsNodes(typeNode* dest){
	for(int i = 0; i < n_parsNodes; i++){
		dest[i].id_node = parsNodes[i].id_node;
		dest[i].partialParsimony = parsNodes[i].partialParsimony;
		dest[i].father = parsNodes[i].father;
		dest[i].number_of_sons = parsNodes[i].number_of_sons;
		dest[i].depth = parsNodes[i].depth;
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
		dst[i].id_node = src[i].id_node;
		dst[i].partialParsimony = src[i].partialParsimony;
		dst[i].father = src[i].father;
		dst[i].number_of_sons = src[i].number_of_sons;
		dst[i].depth = src[i].depth;
		for (int j = 0; j < src[i].number_of_sons; j++){
			dst[i].sons_ids[j] = src[i].sons_ids[j];
		}

		for (int j = 0; j < n_sites; j++){
			dst[i].characters[j] = src[i].characters[j];
		}
	}
}

void PARS::restaurar(typeNode* dst, const typeNode* src, int index_node){
	int node = index_node;
	while (node != num_internal_nodes - 1){
		dst[node].id_node = src[node].id_node;
		dst[node].partialParsimony = src[node].partialParsimony;
		dst[node].father = src[node].father;
		dst[node].number_of_sons = src[node].number_of_sons;
		for (int i = 0; i < src[node].number_of_sons; i++){
			dst[node].sons_ids[i] = src[node].sons_ids[i];
		}
		for (int i = 0; i < n_sites; i++){
			dst[node].characters[i] = src[node].characters[i];
		}
		node = src[node].father;
	}
	dst[node].id_node = src[node].id_node;
	dst[node].partialParsimony = src[node].partialParsimony;
	dst[node].father = src[node].father;
	dst[node].number_of_sons = src[node].number_of_sons;
	for (int i = 0; i < src[node].number_of_sons; i++){
		dst[node].sons_ids[i] = src[node].sons_ids[i];
	}
	for (int i = 0; i < n_sites; i++){
		dst[node].characters[i] = src[node].characters[i];
	}
}



void PARS::cloneRefSeq(char* dst, char* src){
	for (int i = 0; i < n_sites*n_sequences; i++){
		dst[i] = src[i];
	}
}

/**
*	Main method. Parameters:
*	fic_tree: name of the input file containing the reference tree
*	t1: variable to store the starting execution time
*	t2: variable to store the ending execution time
*/

int PARS::run (string fic_tree, double &t1, double &t2, int _mode)
{
    int i;
	//Reading and initializing input tree
	int** matrix;
	genInitialTree(fic_tree);
	initializeReferenceSequences();
	initializeQuerySequences();
	initializeParsTree();
	int mode = _mode;
	
	totalParsimony = calculateParsimonyRefTree (t1, t2);
	/*if(reference_sequences!=NULL)
	{
		for (i=0; i<n_sequences; i++)
				delete[] reference_sequences[i];
		delete[] reference_sequences;
	}
	reference_sequences = NULL;
	if(query_sequences!=NULL)
	{
		for (i=0; i<n_queries; i++)
				delete[] query_sequences[i];
		delete[] query_sequences;
	}
	query_sequences=NULL;
		if (refSites != NULL)
	{
		delete refSites;
	}
	refSites = NULL;
	if (querySites!=NULL)
	{
		delete querySites;
	}
	querySites=NULL;
	if (alphabet != NULL)
	{
		delete alphabet;
	}
	alphabet = NULL;
	if (refTree != NULL)
	{  
		delete refTree;  
	}
	refTree=NULL;
	if (seqReader != NULL)
		delete seqReader;
	seqReader = NULL;
	if (seqReader2 != NULL)
		delete seqReader2;
	seqReader2 = NULL;*/
	/*matrixParsimony=(int**)malloc(n_queries * sizeof(int*));
	for (int n=0; n<n_queries; n++){
		matrixParsimony[n]=(int*)malloc(3*sizeof(int));
		matrixParsimony[n][2] = INT_MAX;
	}*/
	/*copy_parsNodes = new typeNode[num_internal_nodes];
	for (int j = 0; j < num_internal_nodes; j++){
		copy_parsNodes[j].characters = new char[n_sites];
	}
	cloneParsNodes(copy_parsNodes);*/
	/*if(reference_sequences!=NULL)
	{
		for (i=0; i<n_sequences; i++)
				delete[] reference_sequences[i];
		delete[] reference_sequences;
	}
	reference_sequences = NULL;
	if(query_sequences!=NULL)
	{
		for (i=0; i<n_queries; i++)
				delete[] query_sequences[i];
		delete[] query_sequences;
	}
	query_sequences=NULL;
	if (refSites != NULL)
	{
		delete refSites;
	}
	refSites = NULL;
	if (querySites==NULL)
	{
		delete querySites;
	}
	querySites=NULL;
	if (alphabet != NULL)
	{
		delete alphabet;
	}
	alphabet = NULL;
	if (refTree != NULL)
	{  
		delete refTree;  
	}
	refTree=NULL;
	if (seqReader != NULL)
		delete seqReader;
	seqReader = NULL;
	if (seqReader2 != NULL)
		delete seqReader2;
	seqReader2 = NULL;*/


	switch(mode){
		case 0:	matrix=calculateParsimonyQuerys(t1,t2); break;
		case 1:	matrix=calculateParsimonyQuerysGrueso(t1,t2); break;
		case 2: matrix=calculateParsimonyQuerysFino(t1,t2); break;
		case 3: matrix=calculateParsimonyQuerysSIMD(t1,t2); break;
		default:
			printf("Not correct mode:\n");
	}
	if(reference_sequences!=NULL)
	{
		for (i=0; i<n_sequences; i++)
				delete[] reference_sequences[i];
		delete[] reference_sequences;
	}
	reference_sequences = NULL;
	if(query_sequences!=NULL)
	{
		for (i=0; i<n_queries; i++)
				delete[] query_sequences[i];
		delete[] query_sequences;
	}
	query_sequences=NULL;
	if (refSites != NULL)
	{
		delete refSites;
	}
	refSites = NULL;
	if (querySites==NULL)
	{
		delete querySites;
	}
	querySites=NULL;
	if (alphabet != NULL)
	{
		delete alphabet;
	}
	alphabet = NULL;
	if (refTree != NULL)
	{  
		delete refTree;  
	}
	refTree=NULL;
	if (seqReader != NULL)
		delete seqReader;
	seqReader = NULL;
	if (seqReader2 != NULL)
		delete seqReader2;
	seqReader2 = NULL;
	//matrix = calculateParsimonyQuerys(t1, t2);
	/*for (int i = 0; i < n_queries; i++){
		printf("%d;%d;%d;%d\n", i,matrix[i][0], matrix[i][1], matrix[i][2]);
		
	}*/
	//clonePars(parsNodes, copy_parsNodes);
	//delete[]matrix;
	//free(matrix);
	for (int i=0; i < n_queries; i++){
		free(matrix[i]);
	}
	free(matrix);
	//totalParsimony = 0;
	/***********************************
	************************************
	HERE WE CALL THE PLACEMENT ALGORITHM	
	************************************
	************************************/

	//calculateParsimonyQuerys(t1, t2);
	if (parsNodes!=NULL)
		deleteAuxParsStructures();
	//Pruebas *pruebas = new Pruebas(parsNodes, array_query_sequences, num_internal_nodes, n_sites);
	/*if(reference_sequences!=NULL)
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

	if (mtr!=NULL) delete (mtr);*/
	///////////////////
	if(reference_sequences!=NULL)
	{
		for (i=0; i<n_sequences; i++)
				delete[] reference_sequences[i];
		delete[] reference_sequences;
	}
	if(array_reference_sequences!=NULL)
		delete[] array_reference_sequences;
	if(query_sequences!=NULL)
	{
		for (i=0; i<n_queries; i++)
				delete[] query_sequences[i];
		delete[] query_sequences;
	}
	if(array_query_sequences!=NULL)
		delete[] array_query_sequences;
	if (refTree != NULL)
	{  
		  delete refTree;  
	}

	if (mtr!=NULL) delete (mtr);
	return totalParsimony;
}

typeNode* PARS::getParsNodes(){
	return parsNodes;
}

int PARS::getNumInternalNodes(){
	return num_internal_nodes;
}

int PARS::getNumberOfSequences(){
	return n_sequences;
}

char* PARS::getArrayReferenceSequences(){
	return array_reference_sequences;
}

char* PARS::getArrayQuerySequences(){
	return array_query_sequences;
}

int PARS::get_n_sites(){
	return n_sites;
}

int PARS::getTotalParsimony(){
	return totalParsimony;
}

int PARS::get_n_querys(){
	return n_queries;
}



void PARS::set_seq(char *new_sequences){
	array_reference_sequences = new_sequences;
}

void PARS::set_n_sequences(int new_n_sequences){
	n_sequences = new_n_sequences;
}

void PARS::setParsNodes(typeNode* newParsNodes){
	parsNodes = newParsNodes;
}

void PARS::set_num_internal_nodes(int new_num_internal_nodes){
	num_internal_nodes = new_num_internal_nodes;
}

void PARS::set_n_sites(int _n_sites){
	n_sites = _n_sites;
}

/** Destructor */
PARS::~PARS()
{

 
}




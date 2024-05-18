#include "PhylogeneticTree.h"
#include <errno.h> 
using namespace std;
using namespace bpp;

/**
* PhylogeneticTree default constructor, initializes the class PhylogeneticTree. 
*/
PhylogeneticTree::PhylogeneticTree ()
{
	id=-1;
	parsimonyScore=0;
	alphabet=NULL;
	sites=NULL;
	completeSites=NULL;
	inferredTree=NULL;
	terminal_nodes_seqpositions=NULL;
}
/**
* PhylogeneticTree parametrized constructor, initializes the class PhylogeneticTree. Parameters:
*	_id: numerical identifier for the phylogenetic tree
*/
PhylogeneticTree::PhylogeneticTree (int _id)
{
	id=_id;
	parsimonyScore=0;
	alphabet=NULL;
	sites=NULL;
	completeSites=NULL;
	inferredTree=NULL;
	terminal_nodes_seqpositions=NULL;
}

/**
* PhylogeneticTree parametrized constructor, initializes the class PhylogeneticTree. Parameters:
*	_id: numerical identifier for the phylogenetic tree
*	_code: newick code of the phylogenetic tree
*/
PhylogeneticTree::PhylogeneticTree (int _id, NewickCode _code)
{
	id=_id;
	code.clone(_code);
	parsimonyScore=0;
	alphabet=NULL;
	sites=NULL;
	completeSites=NULL;
	inferredTree=NULL;	
	terminal_nodes_seqpositions=NULL;
}

/**
* 	Method that initializes the BIO++ data structures DNA and dataset. Parameters:
*	_sites:  //BIO++ dataset instance (with gaps)
*	_completeSites:  //BIO++ dataset instance (without gaps)
*/
void PhylogeneticTree::initialize (SiteContainer* _sites, SiteContainer* _completeSites)
{	
	if (alphabet == NULL)	alphabet = new DNA ();	
	if (sites == NULL) sites = _sites;
	if (completeSites == NULL) completeSites = _completeSites;	
}

/**
* 	Method that initializes from the newick code the data structure (BIO++ implementation) associated to the phylogenetic tree. 
*/
void PhylogeneticTree::initializeTree()
{
	if(inferredTree != NULL) { 
		delete(inferredTree);
		delete(terminal_nodes_seqpositions);
	}
	int num_terminal_nodes;
	vector< int > aux_terminal;
	inferredTree=TreeTemplateTools::parenthesisToTree(code.getNewick(), true, TreeTools::BOOTSTRAP);	
	aux_terminal = inferredTree->getLeavesId();
	num_terminal_nodes=aux_terminal.size();
	terminal_nodes_seqpositions=new int [inferredTree->getNumberOfNodes()];
	//Detecting the position of data associated to terminal nodes in the input sequences
	string leaf_name;
	vector<Node *>nodes = inferredTree->getNodes();
	for(int i=0; i<num_terminal_nodes; i++)
	{
		leaf_name=nodes[aux_terminal[i]]->getName();
        terminal_nodes_seqpositions[aux_terminal[i]]=sites->getSequencePosition(leaf_name);
	}

}

/**
* 	Method that initializes the newick code from the data structure (BIO++ implementation) associated to the phylogenetic tree. 
*/
void PhylogeneticTree::initializeNewick ()
{
	code.setNewick(TreeTemplateTools::treeToParenthesis(*inferredTree, false));	
}

/**
* 	Method that deletes the data structure (BIO++ implementation) associated to the phylogenetic tree. 
*/
void PhylogeneticTree::deleteTree()
{
	if (inferredTree != NULL)
		delete(inferredTree);
	inferredTree=NULL;
}

/**
* 	Method that returns the terminal_nodes_seqpositions array. 
*/
int* PhylogeneticTree::getTerminalPositions()
{
	return terminal_nodes_seqpositions;
}

/**
* 	Method that sets a new phylogenetic tree to the object. Parameters:
*	treeData: BIO++ phylogenetic tree data structure 
*/
void PhylogeneticTree::setInferredTree(TreeTemplate<Node>* treeData)
{
	if (inferredTree != NULL)
		delete(inferredTree);
	inferredTree = new TreeTemplate<Node> (*treeData);
}

/**
* 	Method that sets to NULL the inferredTree variable.
*/
void PhylogeneticTree::setInferredTreeNull()
{
	inferredTree = NULL;
}

/**
* 	Method that returns the phylogenetic tree data structure.
*/
TreeTemplate<Node>* PhylogeneticTree::getInferredTree()
{
	return inferredTree;
}

/**
* 	Method that performs a basic topological optimization - for future releases and integration in search methods.
*	treeData: tree data structure to be optimized
*/
string PhylogeneticTree::optimizeParsimony(TreeTemplate<Node>*treeData)
{
	string res;
	res="";

	try
	{				
		DRTreeParsimonyScore* parsimonytree;
		parsimonytree = OptimizationTools::optimizeTreeNNI(new DRTreeParsimonyScore (*treeData, *sites, true),0);
		res = TreeTemplateTools::treeToParenthesis(parsimonytree->getTree(), false); //Por ahora dejao asi
	
		delete(inferredTree);
		inferredTree = new TreeTemplate<Node> (parsimonytree->getTree());
		delete(parsimonytree);
		
	}catch(Exception exp){cout<<"Exception on topological search"<<endl;	}

	return res;
}

/** 	
*	Method that returns the identifier of the PhylogeneticTree object
*/
int PhylogeneticTree::getId ()
{
	return id;
}

/** 	
*	Method that returns the NewickCode of the PhylogeneticTree object
*/
NewickCode* PhylogeneticTree::getNewick ()
{
	return (&code);
}

/** 	
*	Method that returns the parsimonyScore of the PhylogeneticTree object
*/
int PhylogeneticTree::getParsimony()
{
	return parsimonyScore;
}

/** 	
*	Method that sets the parsimonyScore of the PhylogeneticTree object from the input parameter. Parameters:
*	_parsimonyScore: value to be set
*/
void PhylogeneticTree::setParsimony(int _parsimonyScore)
{
	parsimonyScore=_parsimonyScore;	
}

/** 	
*	Method that sets a new NewickCode to the PhylogeneticTree object. Parameters:
*	_code: string containing the newick code 
*/
void PhylogeneticTree::setNewick (string _code)
{
	code.setNewick(_code);
}

/** 	
*	Method that sets the identifier of the PhylogeneticTree object
*/
void PhylogeneticTree::setId(int _id)
{
	id=_id;
}

/** 	
*	Method that calculates the parsimony score of the phylogenetic tree (using the BIO++ implementation, used for verification purposes)
*/
int PhylogeneticTree::setScores ()
{
    if (inferredTree == NULL) {cout<<"Error evaluating tree: not initialized"<<endl; return 0;}
    parsimonyScore = setScoresP (inferredTree);
    cout<<"PARS: "<<parsimonyScore<<endl;
}

/** 	
*	Method that calls the parsimony function implementation of BIO++
*/
int PhylogeneticTree::setScoresP (TreeTemplate<Node>*treeData)
{
	DRTreeParsimonyScore parsimonytree (*treeData, *sites, false, true);
	return (parsimonytree.getScore());	
}

/** 	
*	Method that writes the current tree Newick code into a file. Parameters:
*	file: output file
*/
int PhylogeneticTree::writeNewick (FILE* file)
{
	int res;
	res = code.writeNewick(file);
	return res;
}

/** 
*	Method that reads a newick code from a file and assign it to the tree instance. Parameters: 
*	file: output file
*/
int PhylogeneticTree::readNewick (FILE* file)
{
	int res;
	res = code.readNewick(file);
	return res;
}

/** 
*	Method that reads from a file phylogenetic tree data. Parameters: 
*	file: input file
*/
int PhylogeneticTree::readTree (FILE* file)
{
	char c;
	string c_string;
	if (file == NULL) return 0;
	//Reading id	
	c = fgetc(file);
	c_string="";
	while (c != EOF && c != '\n')
	{
		c_string+=c;		
		c = fgetc(file);
	}
	id=atoi (c_string.c_str());
	//Reading parsimony
	c = fgetc(file);
	c_string="";
	while (c != EOF && c != '\n')
	{
		c_string+=c;		
		c = fgetc(file);
	}
	parsimonyScore=atof(c_string.c_str());
	//Reading newick
	code.readNewick (file);
	return 1;
}

/** 
*	Method that writes to a file phylogenetic tree data. Parameters: 
*	file: output file
*/
int PhylogeneticTree::writeTree(FILE* file)
{
	char* c_string;
	c_string = new char [100];
	strcpy (c_string, "");
	sprintf(c_string,"%d\n%d\n", id, parsimonyScore);
	if (file == NULL) return 0;
	for (int i=0; i < strlen(c_string); i++)
  	{
    		fwrite(&c_string[i], sizeof(char), 1, file);
  	}
	code.writeNewick (file);
	delete(c_string);
	return 1;

}

/** 	
*	Method that clones a PhylogeneticTree instance. Parameters:
*	_tree: input tree
*/
void PhylogeneticTree::clone (const PhylogeneticTree & _tree)
{
	id=_tree.id;
	code.clone(_tree.code);
	parsimonyScore=_tree.parsimonyScore;
	
	if (alphabet == NULL)alphabet=_tree.alphabet;
	if (sites == NULL) sites = _tree.sites;
	if (completeSites == NULL) completeSites = _tree.completeSites;
	//if (inferredTree == NULL) inferredTree = new TreeTemplate<Node>;
	//inferredTree.clone(tree.inferredTree);

	if (_tree.inferredTree != NULL)
    {
        if (inferredTree != NULL)
        {
            delete(inferredTree);
            delete(terminal_nodes_seqpositions);
        }
        inferredTree = new TreeTemplate<Node> (*_tree.inferredTree);

        int num_terminal_nodes;
        vector< int > aux_terminal;
        aux_terminal = inferredTree->getLeavesId();
        num_terminal_nodes=aux_terminal.size();
        //terminal_nodes_seqpositions=new int [num_terminal_nodes];
		terminal_nodes_seqpositions=new int [inferredTree->getNumberOfNodes()];
        for(int i=0; i<num_terminal_nodes; i++)
        {
            terminal_nodes_seqpositions[aux_terminal[i]]=_tree.terminal_nodes_seqpositions[aux_terminal[i]];
        }
    }
}

/**
*	Method that prints to console the data associated to the PhylogeneticTree object
*/
void PhylogeneticTree::printTree()
{
	cout<<"Phylogenetic tree #"<<id<<endl;
	cout<<"Parsimony score: "<<parsimonyScore<<endl;	
	code.printCode();
}



/**	Destructor 	*/
PhylogeneticTree::~PhylogeneticTree(){

	//if(alphabet != NULL)	delete(alphabet);
	//if (completeSites != NULL) delete(completeSites);
	//if (sites != NULL) delete(sites);
	if(inferredTree != NULL) delete(inferredTree);
	if(terminal_nodes_seqpositions!=NULL) delete(terminal_nodes_seqpositions);
}

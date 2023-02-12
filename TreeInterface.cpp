#include "TreeInterface.h"

/** Default constructor */
TreeInterface::TreeInterface ()
{
    id=0;
}

/** Parametrized constructor. Parameters:
* 	_id: identifier of the TreeInterface instance 
*	_tree: PhylogeneticTree object
*/
TreeInterface::TreeInterface (int _id, PhylogeneticTree &_tree)
{
    id=_id;
    tree.clone(_tree); //ojo quitar initmatrix del clone   
}

/** 
*	Method that sets the associated PhylogeneticTree object 
*	_tree: PhylogeneticTree instance 
*/
void TreeInterface::setTree (PhylogeneticTree &_tree)
{
    tree.clone(_tree); //ojo quitar initmatrix del clone
}

/** 
*	Method that sets the identifier of the TreeInterface instance 
*	_id: numerical identifier 
*/
void TreeInterface::setId(int _id)
{
	id=_id;
}

/** 	
*	Methods that sets the scores for the associated PhylogeneticTree object by using the method defined in its implementation (i.e. BIO++) 
*/
int TreeInterface::setScores ()
{
    int res=tree.setScores();
    return res;
}


/** 
*	Method that returns the identifier of the TreeInterface instance 
*/
int TreeInterface::getId()
{
	return id;
}

/** 
*	Method that returns the associated PhylogeneticTree object 
*/
PhylogeneticTree* TreeInterface::getTree()
{
    return (&tree);                 
}

/** 	
*	Method that returns the parsimony score for the associated tree 
*/
int TreeInterface::getParsimony()
{
    return (tree.getParsimony());       
}

/**
*	Method that prints the TreeInterface data to console
*/
void TreeInterface::printTree ()
{
	tree.printTree();
}

/** Destructor 	*/
TreeInterface::~TreeInterface(){
	
}

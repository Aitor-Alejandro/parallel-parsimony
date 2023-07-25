#include "Pruebas.h"

using namespace std;
using namespace bpp;
Pruebas::Pruebas(PARS app){
    parsTest = &app;
    n_sites = parsTest->get_n_sites();
    n_queries = parsTest->get_n_querys();
    nodes_phyl = parsTest->getNumInternalNodes()+1;
    n_sequences = parsTest->getNumberOfSequences();

    treeH = parsTest->getTreeH();
    totalNodes = treeH->getNumNodos();
    numLeafNodes = treeH->getNumLeaves();
  
    ptrArray = new typeNode*[(nodes_phyl)*sizeof(typeNode*)];
    arrayPars = new typeNode[(nodes_phyl)*sizeof(typeNode)];
    arrayNodes = new typeNode[(nodes_phyl-1)*sizeof(typeNode)];
    
    array_querys = parsTest->getArrayQuerySequences();
    array_reference = parsTest->getArrayQuerySequences();
    array_ref_query = new char [n_sites*(n_sequences+1)];

    query = NULL;
    qInternalNode = NULL;

    for (int i = 0; i < nodes_phyl; i++)
        arrayPars[i].characters = new char[n_sites];

    for (int i = 0; i < nodes_phyl-1; i++)
        arrayNodes[i].characters = new char[n_sites];

    parsTest->clonePars(arrayNodes, parsTest->getParsNodes());
    
    for (int i = 1; i < nodes_phyl; i++){
        ptrArray[i] = &arrayNodes[i-1];
    }
}

int Pruebas::testParsRefTree(){
    printf("TEST PARSIMONY REF TREE\n");
    double t1, t2;
    int pars = parsTest->calculateParsimonyRefTree(t1, t2);
    return pars;
}

typeNode* Pruebas::generarParsNodes(){

    for (int i = 0; i < nodes_phyl; i++){
        arrayPars[i].id_node = ptrArray[i]->id_node;
        arrayPars[i].number_of_sons = ptrArray[i]->number_of_sons;
        for (int j = 0; j < arrayPars[i].number_of_sons; j++){
            arrayPars[i].sons_ids[j] = ptrArray[i]->sons_ids[j];
        }
    }
}

void Pruebas::initializeParsForTestQuerys(int id_query, int pars_index){
    int num_sons = arrayNodes[pars_index].number_of_sons;
    /*if (pars_index != 0){

    }*/
    int node_class, node_id;
    int father, child;
    double t1, t2;
    
    
}

int** Pruebas::testParsQuerys(){
    printf("TEST PARSIMONY QUERYS\n");
    qInternalNode = new typeNode;
    qInternalNode->characters = new char [n_sites];
    qInternalNode->number_of_sons = 2;
    qInternalNode->sons_ids[1] = numLeafNodes + 1;
    qInternalNode->sons_ids[1] = qInternalNode->sons_ids[1]|(0<<31);

    ptrArray[0] = qInternalNode;

    parsTest->cloneRefSeq(array_reference, parsTest->getArrayReferenceSequences());

    int id_internal_node = totalNodes-numLeafNodes;
    id_internal_node = id_internal_node|(1<<31);

    int pars = 0;
    double t1, t2;
    for (int i = 0; i < n_queries; i++){
        query = &array_querys[n_sites*i];
        for (int j = 0; j < n_sites; j++){
            array_reference[numLeafNodes*n_sites+j] = query[j];
        }
        
        for (int j = 1; j < nodes_phyl; j++){
            ptrArray[j] = qInternalNode;
            ptrArray[j-1] = &arrayNodes[j-1];
            generarParsNodes();
            initializeParsForTestQuerys(i, j);
        }
    }
    //int pars = parsTest->calculateParsimonyRefTree(t1, t2);
}

Pruebas::~Pruebas(){

}
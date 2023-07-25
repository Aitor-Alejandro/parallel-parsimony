#include "TreeHeuristic.h"

using namespace std;
using namespace bpp;
TreeHeuristic::TreeHeuristic(TreeInterface* interface){
    heuristicTree = interface->getTree();
    arrayLeaves = heuristicTree->getInferredTree()->getLeavesId();
    matrix.resize(arrayLeaves.size());
    for (int i = 0; i < arrayLeaves.size(); i++){
        matrix[i].resize(100);
        matrix[i] = heuristicTree->getInferredTree()->getAncestorsId(arrayLeaves[i]);
        //reverse(matrix[i].begin(), matrix[i].end());
        //matrix[i].push_back(arrayLeaves[i]);
    }
    numNodos = heuristicTree->getInferredTree()->getNodes().size();
    nodosVisitados = (bool*)malloc(numNodos*sizeof (bool));
    for (int i = 0; i < numNodos; i++){
        nodosVisitados[i] = false;
    }
    internal_sequences = NULL;
    root = heuristicTree->getInferredTree()->getRootId();
    correspondencias = NULL;
}

void setInternalSequences(int index_x, char* seq){

}

int TreeHeuristic::setScores(){
    int res = heuristicTree->setScores();
    return res;
}

int TreeHeuristic::getParsimony(){
    return (heuristicTree->getParsimony());
}

int TreeHeuristic::getNumNodos(){
    return numNodos;
}

int TreeHeuristic::getNumLeaves(){
    return numHojas;
}

bool TreeHeuristic::getVisited(int index){
    return nodosVisitados[index];
}

void TreeHeuristic::setFalseVisitados(){
    for (int i = 0; i < numNodos; i++){
        nodosVisitados[i] = false;
    }
}

void TreeHeuristic::setTreuVisitado(int index){
    nodosVisitados[index] = true;
}

void TreeHeuristic::show(){
    printf("-----------------------Show Tree Heuristic---------------------\n");
    //for (int i = 0; i < arrayLeaves.size(); i++)
    //    printf("Leave:%d with id:%d\n", i, arrayLeaves[i]);
    for (int i = 0; i < matrix.size(); i++){
        printf("Leave ID: %d with ancestors:->>>", arrayLeaves[i]);
        for (int j = 0; j < matrix[i].size(); j++){
            printf("-%d", matrix[i][j]);
        }
        printf("\n");
    }
    printf("numNodos:--->>>%d\n", numNodos);
    for (int i = 0; i < numNodos; i++){
        printf("%s - ", nodosVisitados[i] ? "true" : "false");
    }
    printf("ROOT::: %d\n", root);
    printf("NUM HOJAS::: %d\n", arrayLeaves.size());
    printf("ARRAY LEAVES-------\n");
    for (int i = 0; i < arrayLeaves.size(); i++){
        printf("-%d-", arrayLeaves[i]);
    }
    printf("\n");
    printf("ARRAY CORR---------\n");
    for (int i = 0; i < numNodos - arrayLeaves.size(); i++){
        printf("-%d-", correspondencias[i]);
    }
    printf("\n-----------------------End Tree Heuristic----------------------\n");
}

void TreeHeuristic::setCorrespondencias(int *vec){
    correspondencias = vec;
}

int* TreeHeuristic::getCorrespondencias(){
    return correspondencias;
}

int *TreeHeuristic::getArrayLeaves(){
    int* leaves_ptr = arrayLeaves.data();
    return leaves_ptr;
}

void TreeHeuristic::bruteForce(){

}

TreeHeuristic::~TreeHeuristic(){
    if (heuristicTree != NULL) delete(heuristicTree);
    if (correspondencias != NULL) delete(correspondencias);
    if (internal_sequences != NULL) delete(internal_sequences);
    if (nodosVisitados != NULL) delete(nodosVisitados);
}
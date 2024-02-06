#include "TreeHeuristic.h"

using namespace std;
using namespace bpp;
TreeHeuristic::TreeHeuristic(TreeInterface* interface){
    heuristicTree = interface->getTree();
    arrayLeaves = heuristicTree->getInferredTree()->getLeavesId();
    matrix.resize(arrayLeaves.size());
    for (int i = 0; i < arrayLeaves.size(); i++){
        matrix[i].resize((heuristicTree->getInferredTree()->getAncestorsId(arrayLeaves[i]).size()) + 1);
        matrix[i] = heuristicTree->getInferredTree()->getAncestorsId(arrayLeaves[i]);
        matrix[i].insert(matrix[i].begin(), arrayLeaves[i]);
    }
    numNodos = heuristicTree->getInferredTree()->getNodes().size();
    numHojas = arrayLeaves.size();
    numInner = numNodos-numHojas;
    nodosVisitados = (bool*)malloc(numNodos*sizeof (bool));
    for (int i = 0; i < numNodos; i++){
        nodosVisitados[i] = false;
    }

    grafo.resize(numNodos, vector<bool>(numNodos, false));
    for (int i= 0; i < numNodos; i++){
        
    }


    internal_sequences = NULL;
    root = heuristicTree->getInferredTree()->getRootId();
    correspondencias = NULL;
}

void TreeHeuristic::iniciarGrafo(){
    int aux, aux2;
    for (int i = 0; i < matrix.size(); i++){
        //aux = arrayLeaves[i];
        for (int j = 1; j < matrix[i].size(); j++){
            //if j < matrix[i]
            grafo [matrix[i][j-1]][matrix[i][j]] = true;
        }
    }
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

int TreeHeuristic::getNumInnerNodes(){
    return numInner;
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
    /*for (int i = 0; i < numNodos; i++){
        printf("%s - ", nodosVisitados[i] ? "true" : "false");
    }*/
    printf("ROOT::: %d\n", root);
    printf("NUM HOJAS::: %d\n", arrayLeaves.size());
    printf("ARRAY LEAVES-------\n");
    for (int i = 0; i < arrayLeaves.size(); i++){
        printf("-%d-", arrayLeaves[i]);
    }
    printf("\n");
    printf("EVERY CORR-----------\n");
    printf("LEAVES\n");
    for (int i = 0; i < numNodos; i++){
        
        if (!(vectorCorrespondencias[i].sub_id & 0x80000000)){
            //printf("-%d-", i);
            printf("-%d-//-subID %d                ", i, vectorCorrespondencias[i].sub_id & 0x7FFFFFFF);
        }
    }
    printf("\n");
    printf("INNER\n");
    for (int i = 0; i < numNodos; i++){
        
        if (vectorCorrespondencias[i].sub_id & 0x80000000){
            printf("-%d-//-subID %d                ", i, vectorCorrespondencias[i].sub_id & 0x7FFFFFFF);
        }
    }
    printf("\n");
    printf("ARRAY CORR---------\n");
    for (int i = 0; i < numNodos - arrayLeaves.size(); i++){
        printf("-%d-", correspondencias[i]);
    }
    printf("\n");
    /*printf("---------GRAFO---------\n");
    for (int i = 0; i < grafo.size(); i++){
        printf("     -%d", i);
    }
    for (int i = 0; i < grafo.size(); i++){
        printf("Nodo: %d", i);
        for (int j = 0; j < grafo[i].size(); j++){
            printf("-%s", grafo[i][j]? "true":"false");
        }
        printf("\n");
    }*/
    printf("\n-----------------------End Tree Heuristic----------------------\n");
}

void TreeHeuristic::setCorrespondencias(int *vec){
    correspondencias = vec;
}

void TreeHeuristic::setVectorCorresponcencias(nodePhyl* vec, int nodes){
    copy(&vec[0], &vec[nodes], back_inserter(vectorCorrespondencias));
    //vectorCorrespondencias = vec;
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
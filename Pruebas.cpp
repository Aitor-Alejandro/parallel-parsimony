#include "Pruebas.h"

using namespace std;

Pruebas::Pruebas(PARS pars){

    parsTest =  &pars;
    parsimony = parsTest->getTotalParsimony();
}

//reserva de espacio
void Pruebas::iniciarEstructuras(){
    query = new char [parsTest->get_n_sites()];
    copy_parsNodes = new typeNode[parsTest->getNumInternalNodes()];
    parsNodes = new typeNode[parsTest->getNumInternalNodes()];
    for (int i = 0; i < parsTest->getNumInternalNodes(); i++){
        copy_parsNodes[i].characters = new char[parsTest->get_n_sites()];
        parsNodes[i].characters = new char [parsTest->get_n_sites()];
    }
    nodeP = new typeNode;
    nodeP->characters = new char[parsTest->get_n_sites()];
}

void Pruebas::copiaSeguridadParsNodes(){
    parsTest->cloneParsNodes(copy_parsNodes, parsNodes);
}

void Pruebas::copiarCaracteres(int i, int j, int node_class, int node_id){
    /*if (node_class == 0x00000000){
        char* leaves = parsTest->getArrayReferenceSequences();
        int offset = node_id*parsTest->get_n_sites();
        for (int k = 0; k < parsTest->get_n_sites(); k++){
            nodeP->characters[k] = leaves[offset+k];
        }
        nodeP->partialParsimony = 0;
    }
    else{
        for (int k = 0; k < parsTest->get_n_sites(); k++){
            nodeP->characters[k] = copy_parsNodes[node_id].characters[k];
        }
        nodeP->partialParsimony = copy_parsNodes[node_id].partialParsimony;
    }*/
}

void Pruebas::sustituciones(){
    int node_class;
    int node_id;
    char *leaves = parsTest->getArrayReferenceSequences();
    for (int i = 0; i < parsTest->getNumInternalNodes(); i++){
        for (int j = 0; j < copy_parsNodes[i].number_of_sons; j++){
            node_class = copy_parsNodes[i].sons_ids[j] & 0x80000000;
			node_id = copy_parsNodes[i].sons_ids[j] & 0x7FFFFFFF;
            //if (node_class == 0x00000000){//
                //copiarCaracteres(i,j,node_class, node_id);
                if (node_class == 0x00000000)
                    parsTest->genInternalNode(nodeP, &leaves[node_id*parsTest->get_n_sites()], &leaves[node_id*parsTest->get_n_sites()]);
                else
                    parsTest->genInternalNode(nodeP, copy_parsNodes[node_id].characters, copy_parsNodes[node_id].characters);
                int parsimony = parsTest->calculateParsimonyQuerysPub(i, j, nodeP, parsNodes);
                if (this->parsimony != parsimony){
                    //printf("DIF\n");
                    printf("%d---%d VS %d\n",node_id, this->parsimony, parsimony);
                    parsTest->cloneParsNodes(parsNodes, copy_parsNodes);
                }   
            //}
        }
    }
}

void Pruebas:: iniciarPruebas(){
    printf("Iniciando estructuras......\n");
    iniciarEstructuras();
    printf("Estructuras iniciadas\n");
    printf("Copiando parsNodes......\n");
    copiaSeguridadParsNodes();
    printf("parsNodes copiado\n");
    sustituciones();
    printf("FINAL\n");

}
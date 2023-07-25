#ifndef _PRUEBAS_H_
#define _PRUEBAS_H_

#include <iostream>

#include "TreeHeuristic.h"
#include "PARS.h"

class Pruebas{
    private:
        PARS* parsTest;
        TreeHeuristic *treeH;
        typeNode** ptrArray;
        typeNode* arrayNodes;
        typeNode* qInternalNode;
        typeNode* arrayPars;
        int nodes_phyl;
        int totalNodes;
        int numLeafNodes;
        int n_sites;
        int n_queries;
        int n_sequences;
        char* array_querys;
        char* array_reference;
        char* array_ref_query;
        char* query;
    public:
        Pruebas(PARS app);
        int testParsRefTree();
        typeNode* generarParsNodes();
        void initializeParsForTestQuerys(int id_query, int pars_index);
        int** testParsQuerys();
        ~Pruebas();
};

#endif //_PRUEBAS_H_
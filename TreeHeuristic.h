#ifndef _TREEHEURISTIC_H_
#define _TREEHEURISTIC_H_

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <stdio.h>
#include <stdlib.h>
#include <malloc.h>

#include "NewickCode.h"
#include "MersenneTwister.h"
#include "TreeInterface.h"
#include "PhylogeneticTree.h"

#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/TreeTemplateTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/TreeTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Seq/Container/SiteContainer.h>   
#include <Bpp/Seq/Container/SiteContainerTools.h> 
#include <Bpp/Seq/Io/Phylip.h>  
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Phyl/Tree.h> 

using namespace std;
using namespace bpp;

class TreeHeuristic{
    PhylogeneticTree* heuristicTree;
    int root;
    int *correspondencias;
    char** internal_sequences;
    vector <int> arrayLeaves; //vector de las hojas del arbol
    vector<vector<int>> matrix;
    int numHojas, profundidad, numNodos;
    bool* nodosVisitados;
    public:
        TreeHeuristic(TreeInterface* interface);
        int setScores();
        int getParsimony();
        int getNumNodos();
        int getNumLeaves();
        bool getVisited(int index);
        void setFalseVisitados();
        void setTreuVisitado(int index);
        void show();
        void setCorrespondencias(int *vec);
        int* getCorrespondencias();
        int* getArrayLeaves();
        void bruteForce();
        ~TreeHeuristic();
        
};

#endif //_TREEHEURISTIC_H_
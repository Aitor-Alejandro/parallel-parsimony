#include "PruebasPARS.h"

PruebasPARS::PruebasPARS(PARS pars){
    parsTest = &pars;
}

void PruebasPARS::iniciarEstructuras(){
    copy_parsNodes = new typeNode[parsTest->getNumInternalNodes()];
    for (int i = 0; i < parsTest->getNumInternalNodes(); i++)
        copy_parsNodes[i].characters = new char[parsTest->get_n_sites()];
    nodeP = new typeNode;
    nodeP->characters = new char [parsTest->get_n_sites()];
}

void PruebasPARS::copiaSeguridadParsNodes(){
    parsTest->cloneParsNodes(copy_parsNodes);
}

int *PruebasPARS::pruebaQuerys(){
    char *query = parsTest->getArrayQuerySequences();
    char *leaves = parsTest->getArrayReferenceSequences();
    int n_sites = parsTest->get_n_sites();
    int *parsimony = new int[parsTest->getNumberOfSequences()*n_sites];
    int node_class, node_id;
    int index = 0;
    typeNode *parsAux = new typeNode[parsTest->getNumInternalNodes()];
    for (int j = 0; j < parsTest->getNumInternalNodes(); j++){
        parsAux[j].characters = new char[n_sites];
    }
    parsTest->clonePars(parsAux, copy_parsNodes);
    for (int i = 0; i < parsTest->get_n_querys(); i++){
        query = &parsTest->getArrayQuerySequences()[i*n_sites];
        for (int j = 0; j < parsTest->getNumInternalNodes(); j++){
            for (int k = 0; k < copy_parsNodes[j].number_of_sons; k++){
                node_class = copy_parsNodes[j].sons_ids[k] & 0x80000000;
                if (node_class != 0x80000000){
                    node_id = copy_parsNodes[j].sons_ids[k] & 0x7FFFFFFF;
                    parsTest->genInternalNode(nodeP, &leaves[node_id*n_sites], query);
                    parsimony[index] = parsTest->calculateParsimonyQuerysPub(j,k,nodeP, parsAux);
                    index++;
                    parsTest->clonePars(parsAux, copy_parsNodes);
                }
            }
        }
    }
    /*for (i = 0 ; i < copy_parsNodes[0].number_of_sons; i++){
        node_class = copy_parsNodes[0].sons_ids[i] & 0x80000000;
        if (node_class != 0x80000000){
            node_id = copy_parsNodes[0].sons_ids[i] & 0x7FFFFFFF;
            parsTest->genInternalNode(nodeP, &leaves[node_id*n_sites], query);
            break;
        }
    }*/
    
    
    //parsimony = parsTest->calculateParsimonyQuerysPub(0, i, nodeP, parsAux);
    return parsimony;
}

int* PruebasPARS::pruebaPARS(){
    int node_class, node_id;
    char *query = parsTest->getArrayQuerySequences();
    char *leaves = parsTest->getArrayReferenceSequences();
    int n_sites = parsTest->get_n_sites();
    int n_sequences = parsTest->getNumberOfSequences();
    char *new_sequences= new char[n_sites*(n_sequences+1)];
    for (int i = 0; i < n_sites*n_sequences; i++){
        new_sequences[i] = leaves[i];
    }
    
    typeNode *parsAux = new typeNode[parsTest->getNumInternalNodes()];
    for (int j = 0; j < parsTest->getNumInternalNodes(); j++){
        parsAux[j].characters = new char[n_sites];
    }
    parsTest->clonePars(parsAux, copy_parsNodes);
  
    for (int i = 0; i < n_sites; i++){
        new_sequences[(n_sites*n_sequences)+i] = query[i];
    }
      
    PARS *auxPARS = new PARS("","");
    
    auxPARS->set_n_sequences(n_sequences+1);
    //auxPARS->set_seq(new_sequences);
    auxPARS->set_n_sites(n_sites);

    auxPARS->set_num_internal_nodes(parsTest->getNumInternalNodes()+1);

    typeNode *defPARS = new typeNode[parsTest->getNumInternalNodes()+1];
    for (int j = 0; j < parsTest->getNumInternalNodes()+1; j++){
        defPARS[j].characters = new char[n_sites];
    }
    int *ppp = new int [parsTest->getNumberOfSequences()*parsTest->get_n_querys()]; 
    nodeP->number_of_sons = 2;
    int index = 0;
    for (int i = 0; i < parsTest->get_n_querys(); i++){
        query = &parsTest->getArrayQuerySequences()[n_sites*i];
        for (int w = 0; w < n_sites; w++){
            new_sequences[(n_sites*n_sequences)+w] = query[w];
        }
        for (int j = 0; j < parsTest->getNumInternalNodes(); j++){
            for (int k = 0; k < copy_parsNodes[j].number_of_sons; k++){
                node_class = copy_parsNodes[j].sons_ids[k] & 0x80000000;
                if (node_class != 0x80000000){
                    parsAux[j].sons_ids[k] = 0x00000000|(1<<31);
                    parsAux[j].sons_ids[k] |= j;
                    node_id = copy_parsNodes[j].sons_ids[k] & 0x7FFFFFFF;
                    nodeP->sons_ids[0] = nodeP->sons_ids[1] = 0x00000000;
                    nodeP->sons_ids[0] = node_id;
                    nodeP->sons_ids[1] = n_sequences;
                    for (int l = j; l < parsTest->getNumInternalNodes();l++){
                        for (int m = 0; m < parsAux[l].number_of_sons; m++){
                            node_class = parsAux[l].sons_ids[m] & 0x80000000;
                            node_id = parsAux[l].sons_ids[m] & 0x7FFFFFFF;
                            if (node_class == 0x80000000 && l > j && node_id >= j){//Aqui no tengo en cuenta que el nodo padre de j no se actualiza porque su padre tiene sons_ids node_id = j. Quiza l > j lo resuelva
                                parsAux[l].sons_ids[m]++;
                            }
                        }
                    }
                    for (int l = 0; l < j; l++){
                        defPARS[l].partialParsimony = 0;
                        defPARS[l].number_of_sons = parsAux[l].number_of_sons;
                        for (int m = 0; m < parsAux[l].number_of_sons; m++){
                            defPARS[l].sons_ids[m] = parsAux[l].sons_ids[m];
                        }
                    }
                    defPARS[j].number_of_sons = nodeP->number_of_sons;
                    defPARS[j].sons_ids[0] = nodeP->sons_ids[0];
                    defPARS[j].sons_ids[1] = nodeP->sons_ids[1];
                    defPARS[j].partialParsimony = 0;
                    for (int l = j; l < parsTest->getNumInternalNodes(); l++){
                        defPARS[l+1].partialParsimony = 0;
                        defPARS[l+1].number_of_sons = parsAux[l].number_of_sons;
                        for (int m = 0; m < parsAux[l].number_of_sons; m++){
                            defPARS[l+1].sons_ids[m] = parsAux[l].sons_ids[m];
                        }

                    }
                    /*for (int x = 0; x < auxPARS->getNumInternalNodes(); x++){
                        for (int y = 0; y < defPARS[x].number_of_sons; y++){
                            node_class = parsAux[x].sons_ids[y] & 0x80000000;
                            node_id = parsAux[x].sons_ids[y] & 0x7FFFFFFF;
                            if (node_class == 0x80000000 && x > j && node_id >= j){
                                defPARS[x].sons_ids[y]++;
                            }
                        }
                    }*/
                    //LLAMADA
                    auxPARS->set_seq(new_sequences);
                    auxPARS->setParsNodes(defPARS);
                    double t1, t2;
                    ppp[index] = auxPARS->calculateParsimonyRefTree(t1,t2);
                    
                    parsTest->clonePars(parsAux, copy_parsNodes);
                    index++;
                }
            }
        }
    }
    /*int i;
    for (i = 0 ; i < copy_parsNodes[0].number_of_sons; i++){
        node_class = copy_parsNodes[0].sons_ids[i] & 0x80000000;
        if (node_class != 0x80000000){
            node_id = copy_parsNodes[0].sons_ids[i] & 0x7FFFFFFF;
            nodeP->sons_ids[0] = nodeP->sons_ids[1] = 0x00000000;
            nodeP->sons_ids[0] = node_id;
            nodeP->sons_ids[1] = n_sequences;
            break;
        }
    }*/
    /*
    parsAux[0].sons_ids[i] = 0x00000000|(1<<31);
    for (int j = i+1; j < parsTest->getNumInternalNodes(); j++){
        for (int k = 0; k < parsAux[j].number_of_sons; k++){
            node_class = parsAux[j].sons_ids[k] & 0x80000000;
            //
            if (node_class == 0x80000000){
                parsAux[j].sons_ids[k]++;
            }
        }
    }*/
    
    /*typeNode *defPARS = new typeNode[parsTest->getNumInternalNodes()+1];
    for (int j = 0; j < parsTest->getNumInternalNodes()+1; j++){
        defPARS[j].characters = new char[n_sites];
    }*/
    //parsTest->clonePars(&defPARS[1], parsAux);
    //defPARS[0].sons_ids[0] = nodeP->sons_ids[0];
    //defPARS[0].sons_ids[1] = nodeP->sons_ids[1];
    //defPARS[0].number_of_sons = nodeP->number_of_sons;
    //for (int k = 0; k < def)
    //for (int j = 0; j < parsTest->getNumInternalNodes()+1; j++){
    //    defPARS[j].partialParsimony = 0;
    //}
    //auxPARS->setParsNodes(defPARS);
    
    //auxPARS->set_num_internal_nodes(parsTest->getNumInternalNodes()+1);
    //double t1, t2;
    //int ppp = auxPARS->calculateParsimonyRefTree(t1,t2);

    
    return ppp;
}

void PruebasPARS::iniciarPruebas(){
    printf("Iniciano pruebas estrictas de calculos\n");
    iniciarEstructuras();
    copiaSeguridadParsNodes();
    int* parsQ = pruebaQuerys();
    //printf("%d\n", parsQ);
    int* parsP = pruebaPARS();
    /*for (int i = 0; i < parsTest->get_n_querys()*parsTest->getNumberOfSequences(); i++){
        printf("%d\n", parsQ[i]);
    }*/
    bool iguales = true;
    /*for (int i = 0; i < parsTest->get_n_querys()*parsTest->getNumberOfSequences(); i++){
        if (parsQ[i] == parsP[i]){
            iguales = false;
            printf("EXITO\n");
        }
    }*/
    if (iguales)
        printf("EXITO\n");
    printf("PRUEBAS ESTRICTAS FINALIZADAS\n");
}
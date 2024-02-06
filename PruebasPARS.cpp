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

int PruebasPARS::pruebaQuerys(){
    char *query = parsTest->getArrayQuerySequences();
    char *leaves = parsTest->getArrayReferenceSequences();
    int n_sites = parsTest->get_n_sites();
    int parsimony;
    int node_class, node_id;
    int i = 0;
    for (i = 0 ; i < copy_parsNodes[0].number_of_sons; i++){
        node_class = copy_parsNodes[0].sons_ids[i] & 0x80000000;
        if (node_class != 0x80000000){
            //printf("HOJA\n");
            node_id = copy_parsNodes[0].sons_ids[i] & 0x7FFFFFFF;
            parsTest->genInternalNode(nodeP, &leaves[node_id*n_sites], query);
            //printf("FIN HOJA\n");
            break;
        }
    }
    //parsTest->genInternalNode(, query,);
    typeNode *parsAux = new typeNode[parsTest->getNumInternalNodes()];
    for (int j = 0; j < parsTest->getNumInternalNodes(); j++){
        parsAux[j].characters = new char[n_sites];
    }
    parsTest->clonePars(parsAux, copy_parsNodes);
    parsimony = parsTest->calculateParsimonyQuerysPub(0, i, nodeP, parsAux);
    return parsimony;
}

int PruebasPARS::pruebaPARS(){
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
    auxPARS->set_seq(new_sequences);
    auxPARS->set_n_sites(n_sites);
    
    nodeP->number_of_sons = 2;
    int i;
    for (i = 0 ; i < copy_parsNodes[0].number_of_sons; i++){
        node_class = copy_parsNodes[0].sons_ids[i] & 0x80000000;
        if (node_class != 0x80000000){
            //printf("HOJA\n");
            node_id = copy_parsNodes[0].sons_ids[i] & 0x7FFFFFFF;
            nodeP->sons_ids[0] = nodeP->sons_ids[1] = 0x00000000;
            nodeP->sons_ids[0] = node_id;
            nodeP->sons_ids[1] = n_sequences;
            //printf("FIN HOJA\n");
            break;
        }
    }
    
    parsAux[0].sons_ids[i] = 0x00000000|(1<<31);
    for (int j = i+1; j < parsTest->getNumInternalNodes(); j++){
        
        //parsAux[j].father++;
        for (int k = 0; k < parsAux[j].number_of_sons; k++){
            node_class = parsAux[j].sons_ids[k] & 0x80000000;
            //
            if (node_class == 0x80000000){
                parsAux[j].sons_ids[k]++;
            }
        }
    }
    printf("HOLAHOLA\n");
    
    typeNode *defPARS = new typeNode[parsTest->getNumInternalNodes()+1];
    for (int j = 0; j < parsTest->getNumInternalNodes()+1; j++){
        defPARS[j].characters = new char[n_sites];
    }
    printf("ADIOSADIOS\n");
    parsTest->clonePars(&defPARS[1], parsAux);
    defPARS[0].sons_ids[0] = nodeP->sons_ids[0];
    defPARS[0].sons_ids[1] = nodeP->sons_ids[1];
    defPARS[0].number_of_sons = nodeP->number_of_sons;
    //for (int k = 0; k < def)
    for (int j = 0; j < parsTest->getNumInternalNodes()+1; j++){
        defPARS[j].partialParsimony = 0;
    }
    auxPARS->setParsNodes(defPARS);
    
    auxPARS->set_num_internal_nodes(parsTest->getNumInternalNodes()+1);
    double t1, t2;
    int ppp = auxPARS->calculateParsimonyRefTree(t1,t2);

    
    return ppp;
}

void PruebasPARS::iniciarPruebas(){
    printf("Iniciano pruebas estrictas de calculos\n");
    iniciarEstructuras();
    copiaSeguridadParsNodes();
    int parsQ = pruebaQuerys();
    printf("%d\n", parsQ);
    int parsP = pruebaPARS();
    if (parsQ == parsP)
        printf("EXITOSO\n");
    printf("PRUEBAS ESTRICTAS FINALIZADAS\n");
}
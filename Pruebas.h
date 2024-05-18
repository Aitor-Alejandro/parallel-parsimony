#ifndef _PRUEBAS_H_
#define _PRUEBAS_H_

#include "PARS.h"

class Pruebas{
    private:
        int parsimony;//parsimonia del arbol inicial
        char *query;//va guardar el mismo valor que el nodo a sustituir, que se inserara como query
        PARS *parsTest;//objeto pars
        typeNode *copy_parsNodes, *parsNodes;//copias de parsNodes
        typeNode *nodeP;//Nodo que se genera
        /*

              NodeP
               /\
              /  \
             /    \
          query  query

        */

    public:
        Pruebas(PARS pars);
        void iniciarEstructuras();
        void copiaSeguridadParsNodes();
        void sustituciones();
        void copiarCaracteres(int i, int j, int node_class, int node_id);
        void iniciarPruebas();
};

#endif //_PRUEBAS_H_
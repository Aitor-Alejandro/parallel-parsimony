#ifndef _PRUEBASPARS_H_
#define _PRUEBASPARS_H_

#include "PARS.h"

class PruebasPARS{
    private:
        PARS *parsTest;
        typeNode *copy_parsNodes;
        typeNode *nodeP;
    public:
        PruebasPARS(PARS pars);
        void iniciarEstructuras();
        void copiaSeguridadParsNodes();
        int pruebaQuerys();
        int pruebaPARS();
        void iniciarPruebas();
};

#endif //_PRUEBASPARS_H_

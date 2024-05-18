//MAIN PROGRAM 

#include <iostream>
#include "PARS.h"
#include "Pruebas.h"
#include "PruebasPARS.h"
using namespace std;

int main (int argc, char *argv[])
{	double t1, t2;

	if(argc != 5){
		printf("Input syntax error\n\t./PARS Reference_File Query_File Tree_File\n");
		return 0;
	}
	if (strlen(argv[4]) != 1) {
        printf("Not correct mode introduced, just one letter: s-g-f-v\n");
        return 1;
    }
	string seq_file;
	string query_file;
	string tree_fic;
	char mode;
	int exe_mode;
	seq_file=argv[1];
	query_file=argv[2];
	tree_fic=argv[3];
	mode = argv[4][0];
	printf("%c", mode);
	if (mode != 's' && mode != 'g' && mode != 'f' && mode != 'v'){
		printf("Not correct mode introduced, please select mode: s-g-f-v\n");
		return 0;
	}
	if (mode == 's')
		exe_mode=0;
	else if (mode == 'g')
		exe_mode=1;
	else if (mode == 'f')
		exe_mode=2;
	else if (mode == 'v')
		exe_mode=3;
	cout<<seq_file<<" "<<query_file<<" "<<tree_fic<<endl;	


	int pars1, pars2, pars3;
	PARS app (seq_file, query_file);		
	app.run(tree_fic, t1, t2, exe_mode);
	//pars1 = app.getTotalParsimony();
	//printf("final de app\n");
	
	//Pruebas* pruebas = new Pruebas(app);
	//pruebas->iniciarPruebas();
	//PruebasPARS *pPARS = new PruebasPARS(app);
	//pPARS->iniciarPruebas();

	//printf("pruebas creadas\n");
	//pars2 = pruebas->testParsRefTree();
	//if (pars1 == pars2)
	//	printf("IGUALES\n");
	//else{
	//	printf("diferentes\n");
	//	printf("PARS 1: %d\n", pars1);
	//	printf("PARS 2: %d\n", pars2);
	//}
	//int **parsQuerys;
	//parsQuerys = pruebas->testParsQuerys();
	/*if (parsQuerys == NULL){
		printf ("NO SE CUELGA\n");
	}*/
	
}


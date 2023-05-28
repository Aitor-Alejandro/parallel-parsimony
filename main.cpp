//MAIN PROGRAM 

#include <iostream>
#include "PARS.h"
#include "Pruebas.h"
using namespace std;

int main (int argc, char *argv[])
{	double t1, t2;

	if(argc != 4){
		printf("Input syntax error\n\t./PARS Reference_File Query_File Tree_File\n");
		return 0;
	}

	string seq_file;
	string query_file;
	string tree_fic;
	seq_file=argv[1];
	query_file=argv[2];
	tree_fic=argv[3];
	cout<<seq_file<<" "<<query_file<<" "<<tree_fic<<endl;	
	

	int pars1, pars2, pars3;
	PARS app (seq_file, query_file);			
	app.run(tree_fic, t1, t2);
	pars1 = app.getTotalParsimony();
	printf("final de app\n");
	Pruebas* pruebas = new Pruebas(app);
	printf("pruebas creadas\n");
	pars2 = pruebas->testParsRefTree();
	if (pars1 == pars2)
		printf("IGUALES\n");
	else{
		printf("diferentes\n");
		printf("PARS 1: %d\n", pars1);
		printf("PARS 2: %d\n", pars2);
	}
	
}


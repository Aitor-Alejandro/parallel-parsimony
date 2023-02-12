//MAIN PROGRAM 

#include <iostream>
#include "PARS.h"
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
	
	PARS app (seq_file, query_file);			
	app.run(tree_fic, t1, t2);
}


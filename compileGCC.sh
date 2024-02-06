LIBSOURCE="/lib"
g++ main.cpp NewickCode.cpp PhylogeneticTree.cpp TreeHeuristic.cpp TreeInterface.cpp Pruebas.cpp PruebasPARS.cpp PARS.cpp -w -o PARS -fopenmp -lrt -lm -O3 -static -I$LIBSOURCE/include/ -L$LIBSOURCE/lib/ -lbpp-phyl -lbpp-seq -lbpp-core


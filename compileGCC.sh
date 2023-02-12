LIBSOURCE="/lib"
g++ main.cpp NewickCode.cpp PhylogeneticTree.cpp TreeInterface.cpp PARS.cpp -o PARS -fopenmp -lrt -lm -O3 -static -I$LIBSOURCE/include/ -L$LIBSOURCE/lib/ -lbpp-phyl -lbpp-seq -lbpp-core


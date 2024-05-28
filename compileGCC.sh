LIBSOURCE="/lib"
g++ main.cpp NewickCode.cpp PhylogeneticTree.cpp TreeInterface.cpp PARS.cpp -w -o PARSv -mavx2 -fopenmp -lrt -lm -O2 -static -I$LIBSOURCE/include/ -L$LIBSOURCE/lib/ -lbpp-phyl -lbpp-seq -lbpp-core


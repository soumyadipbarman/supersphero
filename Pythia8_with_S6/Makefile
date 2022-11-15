ROOTSYS=/home/soumyadip/Package/Root/Root62606build
ROOT_LIB= -L${ROOTSYS}/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
ROOT_INC= -I${ROOTSYS}/include
CLHEPSYS=/home/soumyadip/Package/clhep/clhep2452/CLHEP_build
CLHEP_LIB=/home/soumyadip/Package/clhep/clhep2452/CLHEP_install/lib

Esv_phcut:
	g++ Esv_phcut.cc EventShape_vector.cc -o Esv_phcut -I${CLHEPSYS} -L${CLHEP_LIB} -O -ansi -pedantic -Wall -D_GNU_SOURCE -std=c++11 -D_GNU_SOURCE -pthread -O2 -g -DNDEBUG `root-config --cflags`-Wl,-rpath,./ `root-config --glibs` -lCLHEP-Vector-2.4.5.2 $(ROOT_LIB) $(ROOT_INC)


clean:
	Esv_phcut


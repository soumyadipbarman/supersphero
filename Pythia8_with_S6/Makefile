
ROOTSYS=/home/suman/HEP_Package/ROOT616build
CLHEPSYS=/home/suman/HEP_Package/CLHEPBUILD
ROOT_LIB= -L${ROOTSYS}/lib -lCore -lImt -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lROOTDataFrame -lROOTVecOps -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic
ROOT_INC= -I${ROOTSYS}/include

#LIB = -L${ROOTSYS}/lib -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic -L/home/tanmay/Software/clhep -lCLHEP

Esv_phcut:
#       g++ -g -I${ROOTSYS}/include -I${CLHEP_PARAM:_PATH}/include -o anal_merge_root anal_merge_root.C $(LIB)
#	g++ -O -ansi -pedantic -Wall -D_GNU_SOURCE -std=c++11 -D_GNU_SOURCE -pthread -O2 -g -DNDEBUG `root-config --cflags`-Wl,-rpath,./ `root-config --glibs` -lCLHEP-Vector-2.4.0.0 -I${ROOTSYS}/include -o Esv_bin_sample Esv_bin_sample.cc EventShape_vector.cc -I/home/suman/Buildclhep/ -L/home/suman/Bulidclhep/                        
	g++ Esv_phcut.cc EventShape_vector.cc -o Esv_phcut -I${CLHEPSYS} -L${CLHEPSYS} -O -ansi -pedantic -Wall -D_GNU_SOURCE -std=c++11 -D_GNU_SOURCE -pthread -O2 -g -DNDEBUG `root-config --cflags`-Wl,-rpath,./ `root-config --glibs` -lCLHEP-Vector-2.4.1.0 $(ROOT_LIB) $(ROOT_INC)


clean:
	Esv_phcut


#g++ Esv_bin_sample.cc EventShape_vector.cc -o Esv_bin_sample -I/home/suman/Buildclhep/ -L/home/suman/Bulidclhep/ -O -ansi -pedantic -Wall -D_GNU_SOURCE -std=c++11 -D_GNU_SOURCE -pthread -O2 -g -DNDEBUG `root-config --cflags`-Wl,-rpath,./ `root-config --glibs` -lCLHEP-Vector-2.4.0.0


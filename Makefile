# C++ compiler
CXX=g++
# GSL >=2.0 required
GSLCFLAGS=`gsl-config --cflags`
GSLLIBS=`gsl-config --libs`
# HDF5 library
HDF5CFLAGS=`pkg-config --cflags hdf5`
HDF5LIBS=`pkg-config --libs hdf5`
# C++ compiler flags
CXXFLAGS =-Wall -std=c++11 -O3 -lstdc++ -fopenmp
# ini file parser files
INISRC =\
		 inih/ini.c \
		 inih/INIReader.cpp
INIINC=-I./inih	
# source files of KoMPoST
KOMPOSTSRC =\
		 src/Main.cpp \
		 src/EventInput.cpp \
		 src/BackgroundEvolution.cpp \
		 src/GreensFunctions.cpp \
		 src/KineticEvolution.cpp \
		 src/ExactFreeStreaming.cpp \
		 src/ScalingVariable.cpp

all: 
	$(CXX) -o KoMPoST.exe $(KOMPOSTSRC) $(INISRC) $(INIINC) $(CXXFLAGS) $(GSLCFLAGS) $(HDF5CFLAGS) $(GSLLIBS) $(HDF5LIBS)

.PHONY: clean
clean:
	rm KoMPoST.exe

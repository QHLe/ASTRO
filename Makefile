
EXE = lttopt

OBJS =  cpp_example.o ADOL-C_sparseNLP.o OCP.o SVector.o MVector.o


ADPATH = $(HOME)/adolc_base/include
ADLIBDIR = $(HOME)/adolc_base/lib64

ADDLIBS = /home/zineus/adolc_base/lib64/libadolc.a $(HOME)/ColPack/lib/libColPack.a

ADDINCFLAGS = -I$(ADPATH)  

IPOPTLIBDIR = /home/zineus/Ipopt-3.12.1/lib

# C++ Compiler command
CXX = g++

# Source files
SRC = cpp_example.cpp ADOL-C_sparseNLP.cpp OCP.cpp SVector.cpp MVector.cpp

# C++ Compiler options
CXXFLAGS = -O3 -pipe -DNDEBUG -pedantic-errors -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long -fopenmp -DIPOPT_BUILD 

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,/home/zineus/Ipopt-3.12.1/lib

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = -I/home/zineus/adolc_base/include -I/home/zineus/Ipopt-3.12.1/include/coin -I/home/zineus/Ipopt-3.12.1/include/coin/ThirdParty

#INCL = -I`$(CYGPATH_W) /home/zineus/packages/Ipopt-3.12.1/build/Ipopt/include/coin`  $(ADDINCFLAGS)

# Linker flags
LIBS = -L${IPOPTLIBDIR} -lipopt -L/home/zineus/coinhsl/lib -lcoinhsl -llapack -lblas -lm  -ldl -lcoinmumps -lblas -lgfortran -lm -lgomp -lquadmath -lpthread -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lcoinmetis
#LIBS = `PKG_CONFIG_PATH=/home/zineus/Ipopt-3.12.1/lib/pkgconfig: pkg-config --libs ipopt`

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(LIBS) $(ADDLIBS)

#cpp_example:
#	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) $(INCL) $(SRC) -o cpp_example $(LIBS) $(ADDLIBS)

clean:
	rm -rf $(EXE) $(OBJS)

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `$(CYGPATH_W) '$<'`
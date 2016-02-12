
EXE = lttopt

OBJS =  main.o ADOL-C_sparseNLP.o


ADPATH = $(HOME)/adolc_base/include
ADLIBDIR = $(HOME)/adolc_base/lib64

ADDLIBS = $(HOME)/adolc_base/lib64/libadolc.a $(HOME)/ColPack/lib/libColPack.a -lmgl

ADDINCFLAGS = -I$(ADPATH)  

IPOPTLIBDIR = $(HOME)/Ipopt-3.12.1/lib

# C++ Compiler command
CXX = g++ -pg

# Source files
SRC = main.cpp ADOL-C_sparseNLP.cpp 

# C++ Compiler options
CXXFLAGS = -O3 -pipe -DNDEBUG -pedantic-errors -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long -fopenmp -DIPOPT_BUILD 

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,$(HOME)/Ipopt-3.12.1/lib

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = -I$(HOME)/adolc_base/include -I$(HOME)/Ipopt-3.12.1/include/coin -I$(HOME)/Ipopt-3.12.1/include/coin/ThirdParty

#INCL = -I`$(CYGPATH_W) $(HOME)/packages/Ipopt-3.12.1/build/Ipopt/include/coin`  $(ADDINCFLAGS)

# Linker flags
LIBS = -L${IPOPTLIBDIR} -lipopt -L$(HOME)/coinhsl/lib -lcoinhsl -lm  -ldl -lcoinmumps -lgfortran -lm -lgomp -lquadmath -lpthread -L/usr/lib/gcc/x86_64-linux-gnu/4.8 -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/4.8/../../.. -lcoinmetis -llapack -lf77blas -lcblas -latlas
#LIBS = `PKG_CONFIG_PATH=$(HOME)/Ipopt-3.12.1/lib/pkgconfig: pkg-config --libs ipopt`

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
	rm -rf $(EXE) $(OBJS) *.eps results*.* 

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `$(CYGPATH_W) '$<'`

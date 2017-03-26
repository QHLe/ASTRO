
NAME = lowthr

EXE = $(NAME)

OBJS =  $(EXE).o $(LTTOpt_DIR)/ADOL-C_sparseNLP.o 

LTTOpt_DIR = $(HOME)/Dropbox/Eclipse_workspace/LTTOpt

ADPATH = $(HOME)/adolc_base/include
ADLIBDIR = $(HOME)/adolc_base/lib64

ADDLIBS = $(HOME)/adolc_base/lib64/libadolc.a $(HOME)/ColPack/lib/libColPack.a

ADDINCFLAGS = -I$(ADPATH)  

IPOPTLIBDIR = $(HOME)/Ipopt-3.12.1/lib

# C++ Compiler command
CXX = g++ -pg #-DSPARSE_HESS

# Source files
SRC = main.cpp $(LTTOpt_DIR)/ADOL-C_sparseNLP.cpp $(LTTOpt_DIR)/OCP.cpp

# C++ Compiler options
CXXFLAGS = -std=c++11 -O3 -pipe -DNDEBUG -pedantic-errors -Wparentheses -Wreturn-type -Wcast-qual -Wall -Wpointer-arith -Wwrite-strings -Wconversion -Wno-unknown-pragmas -Wno-long-long -fopenmp -DIPOPT_BUILD 

# additional C++ Compiler options for linking
CXXLINKFLAGS =  -Wl,--rpath -Wl,$(HOME)/Ipopt-3.12.1/lib

# Include directories (we use the CYGPATH_W variables to allow compilation with Windows compilers)
INCL = -I$(HOME)/adolc_base/include -I$(HOME)/Ipopt-3.12.1/include/coin -I$(HOME)/Ipopt-3.12.1/include/coin/ThirdParty -I$(LTTOpt_DIR)

#INCL = -I`$(CYGPATH_W) /home/zineus/packages/Ipopt-3.12.1/build/Ipopt/include/coin`  $(ADDINCFLAGS)

# Linker flags
LIBS = -L${IPOPTLIBDIR} -lipopt -L$(HOME)/coinhsl/lib -lcoinhsl -llapack -lblas -lm  -ldl -lcoinmumps -lblas -lgfortran -lm -lgomp -lquadmath -lpthread -lcoinmetis
#LIBS = `PKG_CONFIG_PATH=/home/zineus/Ipopt-3.12.1/lib/pkgconfig: pkg-config --libs ipopt`

# The following is necessary under cygwin, if native compilers are used
CYGPATH_W = echo

all: $(EXE)

.SUFFIXES: .cpp .c .o .obj

$(EXE): $(OBJS)
	bla=;\
	for file in $(OBJS); do bla="$$bla `$(CYGPATH_W) $$file`"; done; \
	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) -o $@ $$bla $(LIBS) $(ADDLIBS)

#lttopt:
#	$(CXX) $(CXXLINKFLAGS) $(CXXFLAGS) $(INCL) $(SRC) -o $(EXE) $(LIBS) $(ADDLIBS)

clean:
	rm -rf $(EXE) $(OBJS) results*.* *.eps

.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ $<


.cpp.obj:
	$(CXX) $(CXXFLAGS) $(INCL) -c -o $@ `$(CYGPATH_W) '$<'`

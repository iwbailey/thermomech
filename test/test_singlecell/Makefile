
SDIR = ../../src/
IDIR = $(SDIR)
CXX = g++
OPT = -g

OBJ = $(SDIR)creep_law.o $(SDIR)stiffness_matrix.o $(SDIR)fault.o \
	$(SDIR)thermofault.o $(SDIR)fault_cooling.o \
	$(SDIR)romberg_integration.o

all: test_singlecell

test_singlecell: test_singlecell.cpp $(OBJ) input_parameters.h
	$(CXX) $(filter %.o, $^) $(filter %.cpp, $^) -o $@ $(OPT) -I$(IDIR)

%.o:%.cpp %.h
	$(CXX) $< -c $(OPT) -I$(IDIR) -o $@

clear: clean
	rm -f test_loading

clean:
	rm -f $(SDIR)*.o

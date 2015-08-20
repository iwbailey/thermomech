
SDIR = ../../src/cpp/
IDIR = $(SDIR)
CXX = g++
#OPT = -O3
OPT = -g

all: test_circlerupt

test_circlerupt: test_circlerupt.cpp $(SDIR)creep_law.o $(SDIR)stiffness_matrix.o $(SDIR)fault.o $(SDIR)thermofault.o $(SDIR)earthquake.o $(SDIR)fault_cooling.o $(SDIR)romberg_integration.o $(SDIR)temperature_profile.o $(SDIR)strength_profile.o $(SDIR)side_creep.o
	$(CXX) $^ -o $@ $(OPT) -I$(IDIR)

%.o:%.cpp %.h
	$(CXX) $< -c $(OPT) -I$(IDIR) -o $@

clear: clean
	rm -f test_loading

clean:
	rm -f $(SDIR)*.o


SDIR = ../../src/cpp/
IDIR = $(SDIR)
CXX = g++
OPT = -O2

all: test_loading_noheat

test_loading_noheat: test_loading_noheat.cpp $(SDIR)creep_law.o $(SDIR)stiffness_matrix.o $(SDIR)fault.o $(SDIR)side_creep.o $(SDIR)temperature_profile.o
	$(CXX) $^ -o $@ $(OPT) -I$(IDIR)

%.o:%.cpp %.h
	$(CXX) $< -c $(OPT) -I$(IDIR) -o $@

clear: clean
	rm -f test_loading

clean:
	rm -f $(SDIR)*.o

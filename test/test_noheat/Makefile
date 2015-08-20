
SDIR = ../../src/cpp/
IDIR = $(SDIR)
CXX = g++
OPT = -O3
#OPT = -g

# Common objects
OBJ = $(SDIR)creep_law.o $(SDIR)stiffness_matrix.o $(SDIR)fault.o $(SDIR)earthquake.o $(SDIR)temperature_profile.o $(SDIR)strength_profile.o

all: test_bz1996 test_noheat

test_bz1996: test_bz1996.cpp $(OBJ) $(SDIR)benzion_creep.o input_parameters.h
	$(CXX) $(filter %.cpp, $^) $(filter %.o, $^) -o $@ $(OPT) -I$(IDIR)

test_noheat: test_noheat.cpp $(OBJ) $(SDIR)side_creep.o input_parameters.h
	$(CXX) $(filter %.cpp, $^) $(filter %.o, $^) -o $@ $(OPT) -I$(IDIR)

%.o:%.cpp %.h
	$(CXX) $< -c $(OPT) -I$(IDIR) -o $@

clear: clean
	rm -f test_bz1996 test_noheat

clean:
	rm -f $(SDIR)*.o

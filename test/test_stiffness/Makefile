
SDIR = ../../src/
IDIR = $(SDIR)
CXX = g++
OPT = -g

all: test_stiffness test_selfstiffness

test_stiffness: test_stiffness.cpp $(SDIR)stiffness_matrix.o
	$(CXX) $^ -o $@ $(OPT) -I$(IDIR)

test_selfstiffness: test_selfstiffness.cpp $(SDIR)stiffness_matrix.o
	$(CXX) $^ -o $@ $(OPT) -I$(IDIR)

%.o:%.cpp %.h
	$(CXX) $< -c $(OPT) -I$(IDIR) -o $@

clear: clean
	rm -f test_cooling

clean:
	rm -f $(SDIR)*.o *~


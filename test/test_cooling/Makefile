
SDIR = ../../src/cpp/
IDIR = $(SDIR)
CXX = g++
OPT = -O2

# Common objects
OBJ = $(SDIR)romberg_integration.o $(SDIR)fault_cooling.o

all: test_cooling test_cardwell

test_cooling: test_cooling.cpp $(OBJ)
	$(CXX) $^ -o $@ $(OPT) -I$(IDIR)

test_cardwell: test_cardwell.cpp $(OBJ)
	$(CXX) $^ -o $@ $(OPT) -I$(IDIR)

fault_cooling.o: $(SDIR)fault_cooling.cpp $(SDIR)romberg_integration.h
	$(CXX) $(SDIR)fault_cooling.cpp -I$(IDIR) -c $(OPT) -o $@

%.o:%.cpp %.h
	$(CXX) $< -c $(OPT) -I$(IDIR) -o $@

clear: clean
	rm -f test_cooling

clean:
	rm -f $(SDIR)*.o *~

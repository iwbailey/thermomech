
SDIR = ../../src/
IDIR = $(SDIR)
CXX = g++
OPT = -O2

PROGS = test_faultstrength test_strengthprofile test_strength_evol

all: $(PROGS)

test_faultstrength: test_faultstrength.cpp $(SDIR)creep_law.o $(SDIR)side_creep.o $(SDIR)benzion_creep.o $(SDIR)temperature_profile.o  $(SDIR)strength_profile.o $(SDIR)units.h
	$(CXX) $^ -o $@ $(OPT) -I$(IDIR)

test_strengthprofile: test_strengthprofile.cpp $(SDIR)strength_profile.o $(SDIR)temperature_profile.o $(SDIR)benzion_creep.o $(SDIR)creep_law.o $(SDIR)units.h
	$(CXX) $^ -o $@ $(OPT) -I$(IDIR)

test_strength_evol: test_strength_evol.cpp $(SDIR)romberg_integration.o $(SDIR)fault_cooling.o $(SDIR)creep_law.o $(SDIR)strength_profile.o $(SDIR)temperature_profile.o $(SDIR)units.h
	$(CXX) $^ -o $@ $(OPT) -I$(IDIR)

fault_cooling.o: $(SDIR)fault_cooling.cpp $(SDIR)romberg_integration.h
	$(CXX) $(SDIR)fault_cooling.cpp -I$(IDIR) -c $(OPT) -o $@

%.o:%.cpp %.h
	$(CXX) $< -c $(OPT) -I$(IDIR) -o $@

clear: clean
	rm -f $(PROGS)

clean:
	rm -f $(SDIR)*.o

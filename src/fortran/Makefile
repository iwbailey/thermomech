F90 = gfortran
OPT = -std=f95 -fall-intrinsics -O3 -pedantic -Wall

# Modules
MOD = stiffnessmatrix.o fault_parameters.o

# Subroutines
SUB = initialize.o \
	add_randomness.o \
	slipdef2stress.o \
	find_hypocenters.o \
	calc_failures.o

# List of executables to be built within the package
PROGRAMS = slipinitdb2 slipdb2

# "make" builds all
all: $(PROGRAMS)

slipinitdb2: $(MOD) $(SUB)

slipdb2: $(MOD) $(SUB)

# General rule for building prog from prog.o; $^ (GNU extension) is
# used in order to list additional object files on which the
# executable depends
%: %.o
	$(F90) $(OPT) -o $@ $^ $(LDFLAGS)

# General rules for building prog.o from prog.f90 or prog.F90; $< is
# used in order to list only the first prerequisite (the source file)
# and not the additional prerequisites such as module or include files
%.o: %.f90
	$(F90) $(OPT) -c $<

%.o: %.F90
	$(F90) $(OPT) -c $<

# Utility targets
.PHONY: clean veryclean

clean:
	rm -f *.o *.mod

clear: clean
	rm -f $(PROGRAMS)

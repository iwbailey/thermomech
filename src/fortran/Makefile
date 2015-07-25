F90 = gfortran
SDIR = ./
BDIR = ./
OPT = -std=f95 -fall-intrinsics -O3 -pedantic -Wall

SUB = $(SDIR)stiffnessmatrix.f90 $(SDIR)initialize.f90 $(SDIR)add_randomness.f90 $(SDIR)slipdef2stress.f90

all: slipinitdb2 slipdb2

slipinitdb2: $(SDIR)slipinitdb2.f90 $(SUB)
	$(F90) $< $(SUB) -o $@ $(OPT)

slipdb2: $(SDIR)slipdb2.f90 $(SUB2)
	$(F90) $< $(SUB) -o $@ $(OPT)

clean:
	rm slipdb2 slipinitdb2

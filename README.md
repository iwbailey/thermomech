# thermomech
Thermomechanical Fault

# File structure
There are many versions of this code. Not all are working. There are some codes
in fortran, c++ and MATLAB.

# Missing file
The code relies on a numerical recipes file romberg_integration.cpp, but I
don't think the license allows me to share that here, so you need to find it
for yourself.

# Getting started
 - Go to the tests/ folder and run make
 - Go to each subdirectory and you should find a bash script that runs a test
   on that particular element

# Tests
Each test in its own subfolder. To change parameters of the test, change the
.cpp file(s), then make. To run the test and generate check plots. There is a
bash script called run_test.sh in each folder.

## test_stiffness
This makes sure the stiffness (stress transfer) matrix is being calculated correctly

## test_strength
This calculates the strength of the fault under background temperature
conditions and compares it to the original Ben-Zion (1996) calculation.

## test_loading
This applies the loading for a number of time steps without any earthquake and
calculates the impact on the shear stress, slip-deficit and creep rate.

## test_noheat

These tests produce the fault algorithm with no heat generation.

### test_bz1996
This algorithm reproduces the Ben-Zion 1996 algorithm using the c++ objects in
this code.

### test_noheat
The algorithm uses the temperature based slip velocity, but does not generate
any heat from the slip.


## test_singlecell

## test_imposed_eqk

## test_cooling

module fault_parameters

  !  UNITS: fault dimensions in km, time in year, velocities in mm/yr
  real, parameter :: Xlength = 70.0 ! length in km along strike, north of GH
  real, parameter :: Zdepth = 17.5  ! overall depth of fault region in km
  integer, parameter :: nl = 128 ! number of cells along strike
  integer, parameter :: nd = 32 ! number of cells over depth
  !real, parameter :: Xlength2 = 700.0 ! length  of extenal slip sections in km
  real, parameter :: fs = 0.75 ! static coefficient of friction
  real, parameter :: dos = 1.25 ! dynamic overshoot coeff
  real, parameter :: dtaumx = 60.0 ! maximum stress drop in bars
  real, parameter :: Vpl = 35.0 ! plate velocity in mm/yr
  real, parameter :: tauratioz = 4
  real, parameter :: zDB = 10 ! depth where creep rate at tau = fs*Seff equals Vpl
  real, parameter :: xDB = 7.5 ! analogous horizontal position of DB transition
  real, parameter :: t0=125.0

end module fault_parameters

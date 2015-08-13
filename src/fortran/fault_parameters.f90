module fault_parameters

  !  UNITS: fault dimensions in km, time in year, velocities in mm/yr

  ! Length of the fault in km
  real(kind=8), parameter :: Xlength = 70.0

  ! Depth of the fault in km
  real(kind=8), parameter :: Zdepth = 17.5

  ! Number of slip cells along strike and depth
  integer, parameter :: nl = 128
  integer, parameter :: nd = 32

  ! Static coefficient of friction
  real(kind=8), parameter :: fs = 0.75

  ! Dynamic overshoot coeff
  real(kind=8), parameter :: dos = 1.25

  ! Maximum stress drop in bars/cohesion
  real(kind=8), parameter :: dtaumx = 60.0

  ! Gradient of effective normal stress with depth
  real(kind=8), parameter :: dSigmaEff_dz = 180.0

  ! Plate velocity in mm/yr
  real(kind=8), parameter :: Vpl = 35.0

  real(kind=8), parameter :: tauratioz = 4.0

  ! Depth where creep rate at tau = fs*Seff equals Vpl
  real(kind=8), parameter :: zDB = 10.0

  ! Analogous horizontal position of DB transition in km
  real(kind=8), parameter :: xDB = 7.5

  ! Surface temperature in K
  real(kind=8), parameter :: Tsurface = 273 + 20;

  ! Temperature gradient with depth in deg/km
  real(kind=8), parameter :: dTdz = 25.0

  ! Gas constant in J/(mol K)
  real(kind=8), parameter :: Rg = 8.3144

end module fault_parameters

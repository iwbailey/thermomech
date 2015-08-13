program faultslip
  use stiffnessmatrix
  use initialize
  use fault_parameters

  implicit none
  !  This program calculates the evolutions of stress and slip on a
  ! cellular fault incorparating both brittle slip due to static/kineti!
  ! strength drop (as in slip*.f), and ductile slip due to a creep power
  ! law (modified from frxf2d-176.f).

  ! this program differs from corresponding member in the set slipdb*.f
  ! in that here some (e.g. 20%) randomly chosen cells are giver creep
  ! propoerties of depth z = (or >) zDB.
  !
  real(kind=8), parameter :: dt = 4.0/356.0 ! time increment (yr)

  ! Initial time
  real(kind=8), parameter :: t0 = 125.0

  real(kind=8), dimension(nl,nd,nd) :: fi
  real(kind=8), dimension(nl,nd) :: fes, fen
  real(kind=8), dimension(nl,nd) :: u1, u2, tau, taus, taua, taud, tauf, tau0, &
       crp, activEnergy, temperature, vcrp, duc

  integer i, j, it, ihypo, jhypo, nhypo, nEQ

  character(200) ofilename2, ifilename

  real(kind=8) upl, zcrpDB, zrcrp, t, SeffDB

  ifilename = './iofiles/stressdrops.unif.12pm6.txt'
  ofilename2 ='./iofiles/lastdb2.unif.t150.dtau12pm6'

  !  Define ff(i-k,j,l) by tau(i,j)=[Sum over k,l]ff(i-k,j,l)*slip(k,l);
  !  here i,k=1,nl and j,l=1,nd.  Shear modulus is 300 kbars.
  !  Function ff(i-k,j,l) is equated to f(ik,j,l) where ik = ile1,nl
  call get_stiffnessi( fi, nl, nd, Xlength, Zdepth )
  fes = 0.0
!  call get_stiffnesse( fes, nl, nd, Xlength, Zdepth, Xlength2, Zdepth, .true., &
!       0.0 )
  fen = 0.0
!  call get_stiffnesse( fen, nl, nd, Xlength, Zdepth, Xlength2, Zdepth, .false., &
!       0.0 )

  ! Set the static strength
  taus = faulttaus( nl, nd, Zdepth/nd, dtaumx, fs, dSigmaEff_dz)

  ! Set the creep coefficients
  SeffDB = dSigmaEff_dz*zDB
  crp = faultcrp( nl, nd, Xlength/nl, Zdepth/nd, Vpl, fs, SeffDB, tauratioz, zDB, xDB)

  ! Add randomness to the creep coefficientds
  zcrpDB = Vpl/((fs*SeffDB)**3)
  zrcrp = 3.0*log(tauratioz)/(Zdepth - zDB)
  call add_randomness( crp, nl, nd, Xlength/nl, Zdepth/nd, xDB, zDB, zcrpDB, zrcrp  )

  ! Set the activation energy
  activEnergy = 0.0
  temperature = faulttemperature( nl, nd, (Zdepth/nd), Tsurface, dTdz )

  ! Set the arrest stress based on the input file of static stress drops
  taua = faulttaua( ifilename, nl, nd, taus )

  ! Set the dynamic stress
  taud = taus - ( taus - taua )/dos

  ! Initial stress in bars
  tau0 = min( taus - 0.5*dtaumx, 0.95*( Vpl/crp )**(1.0/3.0))

  ! Lock the fault; initialize fault offset from plate motion u1, u2
  tauf = taus
  u1 = 0.0
  u2 = 0.0

  ! impose 125 yr of plate motion
  t = t0
  nEQ = 0
  upl = t*vpl

  !	do it=1,25*356/7+1        ! external loop for model evolution
  do it = 1, 25*356/4+1        ! external loop for model evolution
     !	do it=1,25*365/3+1        ! external loop for model evolution

     write(*,*)'model run time index =',it,', time = ',t

     ! Calculate stresses from internal slip deficit and add stress due to back
     ! slip of bounding large exterior cells at current cycle
     call slipdef2stress( tau, tau0, u1-upl, fi, nl, nd )
     tau = tau + ( fes + fen )*(-vpl*t)

    ! Check for negative stress
     if( minval(tau).lt.0 ) write(*,*)'negative stress at i,j=',i,j

     ! Find location of hypocenter
     call find_hypocenters( nhypo, ihypo, jhypo, tau, taus, nl, nd )

     ! Print if found
     if(nhypo.ge.1)then
        write(*,*)'nhypo= ',nhypo
        write(*,*)'ihypo,jhypo = ',ihypo,jhypo
        nEQ = nEQ+1

      ! Compute failure iterations
        call calc_failures( tau, u2, nl, nd, fi, taus, taud, taua )
     end if

     ! Calc vcrp(i,j) for next time step
     vcrp = crp*(tau**3)*exp( activEnergy / (Rg*temperature) )

     ! Update fault motion to account for last brittle slip and next
     ! creep motion
     u1 = u2 + vcrp*dt

     ! Check and correct for overshoot
     if( maxval( u1 - upl) .gt. 0.0 )then
        do i = 1,nl
           do j = 1,nd
              if( u1(i,j) - upl .gt. 0.0 )then
                 write(*,*)'overshoot at i,j=',i,j
                 u1(i,j)=upl-1
              end if
           end do
        end do
     endif

     ! Update time and plate motion over next time step;
     u2 = u1
     t = t+dt
     upl = t*vpl

  enddo  ! END loop on time steps

  ! Report total number of earthquakes
  write(*,*)nEQ, ' earthquake(s)'

  ! Write the data used to start the next program
  open(2,file=ofilename2,form='unformatted')
  write(2)upl,u1,duc    ! data in mm
  close(2)
  write(*,*)'Written to ', ofilename2

  stop

    !  ! TMP debugging
    !  open( 10, file='tmp2.out')
    !  do i=1,nl
    !     do j=1,nd
    !        write( 10, *) tau(i,j)
    !     end do
    !  end do
    ! stop

end program faultslip

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!===============================================================================

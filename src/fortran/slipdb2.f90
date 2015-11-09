program faultslip
  use stiffnessmatrix
  use initialize
  use fault_parameters

  ! This program calculates the evolutions of stress and slip on a
  ! cellular fault incorparating both brittle slip due to static/kinetic
  ! strength drop (as in slip*.f), and ductile slip due to a creep power
  ! law (modified from frxf2d-176.f).

  ! this program differs from corresponding member in the set slipdb*.f
  ! in that here some (e.g. 20%) randomly chosen cells are giver creep
  ! propoerties of depth z = (or >) zDB.

  implicit none

  real(kind=8), parameter :: dtmx = 3.0/356.0 ! max time increment (yr)
  real(kind=8), parameter :: tSouth = 400 ! originally 150
  real(kind=8), parameter :: tNorth = 400 ! originally 200

  real(kind=8), dimension(nl,nd,nd) :: fi
  real(kind=8), dimension(nl,nd) :: fes, fen
  real(kind=8), dimension(nl,nd) :: u1, u2, tau, taus, taud, taua, tauf, tau0, &
       crp, vcrp, duc, taustart, dtau, activEnergy, temperature

  integer i, j, it, ihypo, jhypo, nhypo, nhypo1, nSlip

  character(200) ifilename, ofilename, ifile_dtau, ifile_crp, ifile_activ
  real(kind=8) dt, t1

  real(kind=8) upl, taumx, taumn, t, avgSlip, avgStressDrop

  write(*,*) 'cohesion =',dtaumx

  ifilename='./iofiles/lastdb2.unif.t150.dtau12pm6'
  ofilename='./iofiles/slipdb2_test.unif.dtau12pm6'
  ifile_dtau  = './ifiles/stressdrops.unif.12pm6.in'
  ifile_crp = './ifiles/creepcoef.arrh.in'
  !ifile_crp = './ifiles/creepcoef.bz1996.in'
  ifile_activ = './ifiles/activEnergy.arrh.in'

  !	write(*,*) 'dtmx=',dtmx
  dt = dtmx              ! initial time step

  !  Define ff(i-k,j,l) by tau(i,j)=[Sum over k,l]ff(i-k,j,l)*slip(k,l);
  !  here i,k=1,nl and j,l=1,nd.  Shear modulus is 300 kbars.
  !  Function ff(i-k,j,l) is equated to f(ik,j,l) where ik = 1,nl
  call get_stiffnessi( fi, nl, nd, Xlength, Zdepth )
  !  call get_stiffnesse( fes, nl, nd, Xlength, Zdepth, Xlength2, Zdepth, .true., &
  !       0.0 )
  fes = 0.0
  !  call get_stiffnesse( fen, nl, nd, Xlength, Zdepth, Xlength2, Zdepth, .false., &
  !       0.0 )
  fen = 0.0

  ! Set the creep coefficients
  crp = read_faultvalues( ifile_crp, nl, nd )

  ! Set the activation energy
  activEnergy = read_faultvalues( ifile_activ, nl, nd )
  !activEnergy = 0.0
  temperature = faulttemperature( nl, nd, (Zdepth/nd), Tsurface, dTdz )

  ! Set the static strength
  taus = faulttaus( nl, nd, Zdepth/nd, dtaumx, fs, dSigmaEff_dz)

  ! Set the arrest stress based on the input file of static stress drops
  taua = read_faultvalues( ifile_dtau, nl, nd )
  taua = taus - taua

  ! Set the dynamic stress
  taud = taus - ( taus - taua )/dos

  ! Initial stress in bars
!  tau0 = min( taus - 0.5*dtaumx, 0.95*( Vpl/crp )**(1.0/3.0))
  ! Initial stress in bars
  tau0 = min( taus - 0.5*dtaumx, &
       0.99*( (Vpl/crp)*exp( activEnergy/(Rg*temperature) ) )**(1.0/3.0))

  ! Read the slip deficit from the output of the initialization program
  open(1,file=ifilename,form='unformatted')
  read(1)upl, u1, duc
  close(1)
  t = upl/Vpl                      ! time in yr
  t1 = t                           ! initial time for current run
  write(*,*) 'end time in file last* =',t1,'yr, upl = ',upl/1000.,'m'

  ! Start the output file
  open(2,file=ofilename,form='unformatted')
  write(2)upl,u1

  ! Lock the fault; intialize fault slip u2 immediately after brittle failure
  ! (u1 is fault slip immediately before brittle failure)
  tauf = taus
  u2 = u1


 !c======================================================================
  do it = 1,3500000          ! external loop for model evolution

     !c	write(*,*) '---------------'
     !c	write(*,*) 'model run time index =',it,', time =',t
     !c	write(*,*) '---------------'

     !c--------------------------------------------------
     !      Initialize Stress drop
     dtau = 0.0 ! stress drop for each cell
     taustart = 0.0 ! stress before eqk

     taumn=999.9
     taumx=0.0

     ! Calculate stresses from internal slip deficit and add stress due to back
     ! slip of bounding large exterior cells at current cycle
     call slipdef2stress( tau, tau0, u1-upl, fi, nl, nd )

     ! add stresses due to back slip of bounding large exterior cells at
     ! current cycle. (initial fault configuration is given by output of
     ! slipinit.f; this corresponds to 150 yr of plate and computational
     ! grid motions, during which the large cells are locked.)

     ! Exterior large cells do not fail simultaneously; 1857-type event occur at
     ! t=150; 1906-type at t=200:
     if(t.gt.tSouth)then
        tau = tau + fes*(-Vpl*(t-tSouth))
     else
        tau = tau + fes*(-Vpl*t)
     endif

     if(t.lt.tNorth)then
        tau = tau + fen*(-Vpl*t)
     else
        tau = tau + fen*(-Vpl*(t-tNorth))
     endif

     ! Check for negative stress
     if( minval(tau).lt.0.0 )then
        write(*,*) 'error: negative stress, ', minval(tau)
        stop -4
     endif

     ! Find location of hypocenter
     call find_hypocenters( nhypo, ihypo, jhypo, tau, taus, nl, nd )

     if(nhypo.ne.0)then
        if(nhypo.gt.1.and.it.gt.1) then
           write(*,*)'Why so many hypocenters?', t, nhypo
           stop -1
        end if
        !--------------------------------------------------
        !record the stress before we transfer/propagate the rupture
        taustart = tau

        ! Compute failure iterations
        call calc_failures( tau, u2, nl, nd, fi, taus, taud, taua )

        ! Calculate the avg static stress drop and avg slip
        nSlip = count( u2-u1 .gt. 1e-8 )
        avgSlip = sum( u2-u1 ) / nSlip
        avgStressDrop = sum( taustart - tau ) / nSlip

        write(*,'(F8.4,1X, F6.2,1X, F6.2,1X, I4,1X, F12.4,1X, F12.4)') &
             t, (ihypo-0.5)*(Xlength/nl), (jhypo-0.5)*(Zdepth/nd), nSlip, avgStressDrop, avgSlip


        ! !--------------------------------------------------
        ! ! calc int(delu) of previous brittle slip increment for output
        ! do j=1,nd
        !    do i=1,nl
        !       delu(i,j)=u2(i,j)-u1(i,j)
        !       diff=nint( 10.*delu(i,j) )            ! delu in 0.1mm
        !       if(diff.gt.32000)diff=32000-diff      ! max(int*2)=32767
        !       if(abs(diff).gt.32000)then
        !          write(*,*) 'error: delu overflow'
        !          stop -3
        !       endif
        !       iduf(i,j) = diff
        !    enddo
        ! enddo

        ! update fault motion to account for last brittle slip; lock the fault
        u1 = u2

     endif ! Continue from here if there are no failure points

     ! Calc vcrp(i,j) for next time step from current stresses
     vcrp = crp*(tau**3)*exp( -activEnergy / (Rg*temperature) )

     ! Calculate the time step
     dt = dtmx


3    continue

     ! Exit if total brittle slip exceeds total plate motion
     if( maxval(u1).gt.upl ) then
        write(*,*) 'error: brittle slip overshoot'
        stop -1
     endif

     ! update fault motion to account for creep motion over assumed next
     ! time step
     u1 = u1 + vcrp*dt

     ! Exit if total creep slip exceeds total plate motion
     if( maxval(u1).gt.upl )then
        write(*,*) 'error: creep overshoot'
        stop -2
     endif

     ! update time and plate motion at end of assumed time step
     t = t+dt
     upl = t*Vpl

!     if(t.ge.200.and.t.lt.200+dtmx/2.)goto 4
     if( t.lt.tNorth .or. t.ge.tNorth+dtmx/2. ) then

        ! check that new nhypo < 2

        !   calculate stresses
        call slipdef2stress( tau, tau0, u1-upl, fi, nl, nd )

        ! add stresses due to back slip of bounding large exterior cells
        if(t.gt.tSouth)then
           tau = tau + fes*(-Vpl*(t-tSouth))
        else
           tau = tau + fes*(-Vpl*t)
        endif

        if(t.lt.tNorth)then
           tau = tau + fen*(-Vpl*t)
        else
           tau = tau + fen*(-Vpl*(t-tNorth))
        endif

        ! Find tentative hypocenter parameters for next time step
        nhypo1 = count( tau > taus )

        ! Check we don't have more than 1 hypocenter
        if(nhypo1.gt.1)then
           ! repeat calculations with dt/2
           ! subtract assumed last time step from time
           ! subtract assumed last creep motion from fault displacement
           ! reduce size of time step; re-calculate creep & plate motions
           t = t-dt
           u1 = u1 - vcrp*dt
           dt = dt/2.

           ! Return to creep part
           goto 3
        else
           write(*,*)'nhypo1 = ', nhypo1
        endif

!4    continue
     end if

     ! accumulate/write output

     if(nhypo.eq.0)then      ! no earthquake; accumulate deluc
        duc = duc + vcrp*dt
     else                    ! earthquake ocurred; dump data

        ! write hypocenter parameters and brittle & creep slip increments at
        ! the end of (for brittle), and up to (for creep), last earthquake

        !c	  write(2)nhypo,ihypo,jhypo,t-dt,iduf,duc

        ! nhypo = number of hypocenters,
        ! ihypo = strike index of hypocenter
        ! jhypo = depth index of hypocenter
        ! t-dt = time of hypocenter
        ! iduf = [int] array of slip deficit in 0.1 mm/yr
        ! duc = [double] total slip from creep
        ! dtau = [double] static stress drop of cells on fault
!        write(2)nhypo,ihypo,jhypo,t-dt,iduf,duc,dtau
        write(2)nhypo,ihypo,jhypo,t-dt,nSlip,avgSlip
!        call flush(2)

        ! re-intialize deluc with creep-sli taustart(i,j)p after last earthquake
        duc = vcrp*dt

        ! ! write displacements at start of next time step in 'restart' file
        ! rewind 1
        ! write(1)upl,u1,duc
!        call flush (1)

     endif

     ! reset time step size and fault slip u2
     u2 = u1

     !  check exit criteria

     if(t-t1.ge.150) exit        ! 150 yr of model evolution
  ! if(t-t1.ge.1) exit        ! 150 yr of model evolution



  enddo    !  on it


  close(2)
  !c	close(3)
 ! TMP debugging
     open( 10, file='tmp2.out')
     do i=1,nl
        do j=1,nd
           write( 10, *) u2(i,j)
        end do
     end do
    stop

!  close(1)

  stop
end program

!-------------------------------------------------------------

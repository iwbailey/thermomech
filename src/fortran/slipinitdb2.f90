program faultslip
  use stiffnessmatrix
  use initialize

  implicit none
  !  This program calculates the evolutions of stress and slip on a
  ! cellular fault incorparating both brittle slip due to static/kineti!
  ! strength drop (as in slip*.f), and ductile slip due to a creep power
  ! law (modified from frxf2d-176.f).

  ! this program differs from corresponding member in the set slipdb*.f
  ! in that here some (e.g. 20%) randomly chosen cells are giver creep
  ! propoerties of depth z = (or >) zDB.
  !

  !  UNITS: fault dimensions in km, time in year, velocities in mm/yr
  real, parameter :: Xlength = 70.0 ! length in km along strike, north of GH
  real, parameter :: Zdepth = 17.5  ! overall depth of fault region in km
  integer, parameter :: nl = 128 ! number of cells along strike
  integer, parameter :: nd = 32 ! number of cells over depth
  real, parameter :: Xlength2 = 700.0 ! length  of extenal slip sections in km
  real, parameter :: fs = 0.75 ! static coefficient of friction
  real, parameter :: dos = 1.25 ! dynamic overshoot coeff
  real, parameter :: dt = 4.0/356.0 ! time increment (yr)
  real, parameter :: dtauavg = 12.0 ! average taus-taua
  real, parameter :: pm = 6 ! plus-minus variation in values of taus-taua
  real, parameter :: dtaumx = 60.0 ! maximum stress drop in bars
  real, parameter :: Vpl = 35.0 ! plate velocity in mm/yr
  real, parameter :: tauratioz = 4
  real, parameter :: zDB = 10 ! depth where creep rate at tau = fs*Seff equals Vpl
  real, parameter :: xDB = 7.5 ! analogous horizontal position of DB transition
  real, parameter :: t0=125.0

  real(kind=8), dimension(nl,nd,nd) :: fi
  real(kind=8), dimension(nl,nd) :: fes, fen
  real(kind=8), dimension(nl,nd) :: u1, u2, tau, taus, taua, tauf, tau0, delu, &
       x, z, crp, vcrp, du
  integer, dimension(nl*nd) :: ifail, jfail

  integer i, j, m, ii, jj, ik, it, ihypo, jhypo, nhypo, ihypo1, jhypo1, &
       indf, nEQ

  character(200) ofilename2, ifilename

  real(kind=8) delmax, upl, xhypo, zcrpDB, zhypo, zrcrp, t

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
  taus = faulttaus( nl, nd, Zdepth/nd, dtaumx, fs, 180.0)

  ! Set the creep coefficients
  crp = faultcrp( nl, nd, Xlength/nl, Zdepth/nd, Vpl, fs, 180.0*zDB, tauratioz, zDB, xDB)

  ! Add randomness to the creep coefficientds
  zcrpDB = Vpl/((fs*180.0*zDB)**3)
  zrcrp = 3.0*log(tauratioz)/(Zdepth - zDB)
  call add_randomness( crp, nl, nd, Xlength/nl, Zdepth/nd, xDB, zDB, zcrpDB, zrcrp  )

  ! Set the arrest stress based on the input file of static stress drops
  taua = faulttaua( ifilename, nl, nd, taus )

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

     ! Calculate stresses from internal slip deficit
     call slipdef2stress( tau, tau0, u1-upl, fi, nl, nd )

     ! Add stress due to back slip of bounding large exterior cells at current
     ! cycle
     tau = tau + ( fes + fen )*(-vpl*t)

     ! Check for negative stress
     if( minval(tau).lt.0 ) write(*,*)'negative stress at i,j=',i,j

     ! Find location of hypocenter
     nhypo=0
     ihypo=0
     jhypo=0
     delmax=-9999
     do j=1,nd
        do i=1,nl
           if(tau(i,j).ge.taus(i,j))then
              nhypo=nhypo+1

              ! Store location of the place that is most over stressed
              if(tau(i,j)-taus(i,j).gt.delmax)then
                 ihypo=i                   ! i coord of slip initiator
                 jhypo=j                   ! j coord of slip initiator
                 delmax=tau(i,j)-taus(i,j)
              endif
           endif
        enddo
     enddo

     ! Print message if hypocenter
     if(nhypo.ge.1)then
        write(*,*)'nhypo= ',nhypo
        write(*,*)'ihypo,ihypo1,jhypo,jhypo1= ',ihypo,ihypo1,jhypo,jhypo1

        nEQ = nEQ+1
        xhypo=x(ihypo,jhypo)
        zhypo=z(ihypo,jhypo)

     endif

     ! Compute failure iterations
     indf = nhypo ! number of current failures
     do while ( indf.gt.0 )

        ! Count failures and compute slip in each cell
        indf=0
        do j=1,nd
           do i=1,nl
              if( tau(i,j)-tauf(i,j).ge.0 )then
                 delu(i,j)=(taua(i,j)-tau(i,j))/fi(1,j,j) ! brittle self-slip
                 indf=indf+1
                 ifail(indf)=i
                 jfail(indf)=j

                 ! reduce tauf to dynamic friction
                 tauf(i,j)=taus(i,j)-(taus(i,j)-taua(i,j))/dos
              endif
           enddo
        enddo

        ! Adjust stresses and brittle displacements; final stresses after
        ! brittle failure will be distributed between arrest stress and dynamic
        ! friction
        do m=1,indf
           i=ifail(m)
           j=jfail(m)
           do jj=1,nd
              do ii=1,nl
                 ik=abs(ii-i)+1
                 tau(ii,jj)=tau(ii,jj)+fi(ik,jj,j)*delu(i,j)
              enddo
           enddo
           u2(i,j)=u2(i,j)+delu(i,j)
        enddo

        ! Report the number of failures
        if(indf.gt.0) then
           write(*,*)'# failed cells =',indf
        endif

     end do ! end iterations over failures

     ! Calc vcrp(i,j) for next time step
     vcrp = crp *(tau**3)

     ! Lock the fault
     tauf = taus

     ! Update fault motion to account for last brittle slip and next
     ! creep motion
     u1 = u2 + vcrp*dt

     ! Check and correct for overshoot
     do j=1,nd
        do i=1,nl
           if(u1(i,j).gt.upl)then
              write(*,*)'overshoot at i,j=',i,j
              u1(i,j)=upl-1
           endif
        enddo
     enddo

     ! Update time and plate motion over next time step;
     u2 = u1
     t = t+dt
     upl = t*vpl

  enddo  ! END loop on time steps

  ! Report total number of earthquakes
  write(*,*)nEQ, ' earthquake(s)'

  ! Write the data used to start the next program
  open(2,file=ofilename2,form='unformatted')
  write(2)upl,u1,du    ! data in mm
  close(2)
  write(*,*)'Written to ', ofilename2

  stop

! ! TMP debugging
!     open( 10, file='tmp.out')
!     do i=1,nl
!        do j=1,nd
!           write( 10, *) taus(i,j) - taua(i,j)
!        end do
!     end do
!     stop


end program faultslip

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!===============================================================================

program faultslip
  use stiffnessmatrix
  ! This program calculates the evolutions of stress and slip on a
  ! cellular fault incorparating both brittle slip due to static/kinetic
  ! strength drop (as in slip*.f), and ductile slip due to a creep power
  ! law (modified from frxf2d-176.f).

  ! this program differs from corresponding member in the set slipdb*.f
  ! in that here some (e.g. 20%) randomly chosen cells are giver creep
  ! propoerties of depth z = (or >) zDB.

  implicit none

  !  UNITS: fault dimensions in km, time in year, velocities in mm/yr
  real, parameter :: Xlength = 70.0 ! length in km along strike, north of GH
  real, parameter :: Zdepth = 17.5  ! overall depth of fault region in km
  real, parameter :: Xlength2 = 700.0 ! length  of extenal slip sections in km
  integer, parameter :: nl = 128 ! number of cells along strike
  integer, parameter :: nd = 32 ! number of cells over depth
  real, parameter :: fs = 0.75 ! static coefficient of friction
  real, parameter :: dos = 1.25 ! dynamic overshoot coeff
  real, parameter :: dtmx = 3.0/356.0 ! max time increment (yr)
  real, parameter :: dtauavg = 12.0 ! average taus-taua
  real, parameter :: pm = 6 ! plus-minus variation in values of taus-taua
  real, parameter :: dtaumx = 60.0 ! maximum stress drop in bars
  real, parameter :: Vpl = 35.0 ! plate velocity in mm/yr
  real, parameter :: tauratioz = 4
  real, parameter :: tauratiox = 4
  real, parameter :: zDB = 10 ! depth where creep rate at tau = fs*Seff equals Vpl
  real, parameter :: xDB = 62.5 ! analogous horizontal position of DB transition
  real, parameter :: t0=125.0
  real, parameter :: p=0.25 ! average stress drop is this proportion of dtaumx
  real, parameter :: s=0.1875 ! range of arrest stresses is 2*this proportion of dtaumx

  real(kind=8), dimension(nl,nd,nd) :: fi
  real(kind=8), dimension(nl,nd) :: fes, fen, u1, u2, tau, taus, taua, tauf, tau0, delu, &
       crp, vcrp, duc, taustart, dtau, duc1
  real, dimension(nl,nd) :: x, z
  integer, dimension(nl*nd) :: ifail, jfail
  integer, dimension(nl,nd) :: iduf
  logical, dimension(nl,nd) :: has_failed

  integer i, j, l, k, m, ii, jj, ik, it, ihypo, jhypo, nhypo, nhypo1,  &
       indf, idum, Nslip, taumx_i , taumx_j , taumn_i , taumn_j, diff

  character(100) ifilename, ofilename
  real(kind=8) dt, rn, ran1, sumdtau, t1, var

  real(kind=8) delmax, temp, temp1, temp2, upl, xrcrp, tau01, tau02, tauamx, &
       taumx, taumn, zrcrp, t
  real xd2, xw2, hcell, hxcell, xi, xcrpdb, zcrpDB, zDBx, xhypo, zhypo, zj, &
       sigmanormalDB, sigmanormalDBx, Seff, SeffDB, SeffDBx, sigmanormal, &
       porepress, porepressDB, porepressDBx

  xd2=Zdepth/(2.*nd)     ! cell half-depth
  xw2=Xlength/(2.*nl)    ! cell half-width

  write(*,*) 'cohesion =',dtaumx


  ifilename='./iofiles/lastdb2.unif.t150.dtau12pm6'
  ofilename='./iofiles/slipdb2_test.unif.dtau12pm6'

  !	write(*,*) 'dtmx=',dtmx
  dt=dtmx              ! initial time step


  !  Define ff(i-k,j,l) by tau(i,j)=[Sum over k,l]ff(i-k,j,l)*slip(k,l);
  !  here i,k=1,nl and j,l=1,nd.  Shear modulus is 300 kbars.
  !  Function ff(i-k,j,l) is equated to f(ik,j,l) where ik = 1,nl
  call get_stiffnessi( fi, nl, nd, Xlength, Zdepth )
  call get_stiffnesse( fes, nl, nd, Xlength, Zdepth, Xlength2, Zdepth, .true., &
       0.0 )
  call get_stiffnesse( fen, nl, nd, Xlength, Zdepth, Xlength2, Zdepth, .false., &
       0.0 )

  !  Positions of mid cells (x(i,j) and z(i,j):

  hcell = Zdepth/nd    ! "cell" height
  hxcell = Xlength/nl  ! "cell" length along strike
  do j = 1,nd
     zj = j
     do i = 1,nl
        xi = i
        x(i,j) = (xi - 0.5)*hxcell ! x coord
        z(i,j) = (zj - 0.5)*hcell  ! z coord
     end do
  end do
  zDBx = zDB !z at x=xDB where creep rate at tau=fs*Seff equals Vpl

  ! brittle strength (static friction), and creep properties:
  do l = 1,nd
     do k = 1,nl

        sigmanormal = 280.*z(k,l)   ! Snormal equated to overburden
        sigmanormalDB = 280.*zDB
        sigmanormalDBx = 280.*zDBx

        porepress = 100.*z(k,l)     ! hydrostatic pore pressure
        porepressDB = 100.*zDB
        porepressDBx = 100.*zDBx

        Seff = sigmanormal-porepress ! Ambient effective stress
        SeffDB = sigmanormalDB - porepressDB
        SeffDBx = sigmanormalDBx - porepressDBx

        ! static brittle strength with dtaumx acting here as cohesion giving
        ! finite strength at depth zero

        taus(k,l)=dtaumx+fs*Seff

        ! Creep parameters: depth and poition along strike zDB and xDB of
        ! "ductile-brittle" transitions, SeffDB, and SeffDBx given above.
        ! We assume Vcreep = crp(x,z)*tau(x,z)**3, where the spatial creep
        ! property distribution crp(x,z) is given either by
        ! crp(x,z) = xcrpDB*exp(xrcrp*(x-xzDB)) + zcrpDB*exp(zrcrp*(z-zDB))
        ! or by the maximum of two analogous 1D functions. xrcrp and zrcrp
        ! are chosen to enforce given stress ratios (say, 10.0) between that
        ! at zDB and that at Zdepth (or xDB and Xlength) if the same Vcreep
        ! is to be produced at both depth (or horizontal locations).
        ! choices of Zdepth and Xlength correspond to depth and along-strike
        ! positions where slip along the SAF is mostly in the form of creep.

        zcrpDB = Vpl/((fs*SeffDB)**3)

        ! Makes Vcreep = Vpl at zDB  if tau = fs*Seff at depth zDB. Note that for
        ! consistency it would be better to use zcrpDB = Vpl/((dtaumx+fs*SeffDB)**3),
        ! otherwise we may get Vcreep > Vpl for high enough values of dtaumx. This change
        ! is implemented in slipdb.stat.f
        zrcrp = 3.0*log(tauratioz)/(Zdepth - zDB)      !  see above

        ! analogous parameters for x-dependence:
        xcrpDB = Vpl/((fs*SeffDBx)**3)
        xrcrp = 3.0*log(tauratiox)/(Xlength - xDB)

        ! A composite 2D distribution of creep properties:

        !         crp(k,l) = zcrpDB*exp(zrcrp*(z(k,l) - zDB))
        !    >             + xcrpDB*exp(xrcrp*(x(k,l) - xDB))

        ! Two 1D distributions of creep properties; crp(k,l) is set as max
        ! of those.

        crp(k,l) = zcrpDB*exp(zrcrp*(z(k,l) - zDB))
        temp = xcrpDB*exp(xrcrp*(x(k,l) - xDB))
        if(temp.gt.crp(k,l))crp(k,l)=temp

        !    EFFECTIVE REMOVAL OF CREEP RESPONSE:
        !         crp(k,l) = 1.0E-08*crp(k,l)

     end do
  end do

  ! assign to 20% randomly chosen 'brittle' cells higher crp(i,j)

  temp1=zcrpDB*exp(zrcrp*(11.25-zDB))
  temp2=zcrpDB*exp(zrcrp*(8.75-zDB))
  !	 temp1=zcrpDB
  !	 temp2=zcrpDB

 ! idum=-1
  call srand(12345)
  l=0
  k=0
  do j=1,nd
     do i=1,nl
        if(z(i,j).lt.13.75.and.x(i,j).lt.66.25)then
!           rn = ran1(idum)
           rn = rand()
           k = k+1
           if(rn.le.0.2)then
              l = l+1
              crp(i,j)=temp2
              if(z(i,j).lt.zDB.and.x(i,j).lt.xDB)crp(i,j)=temp1
           endif
        endif
     enddo
  enddo
  write(*,*) l, 'random selections, i.e.,', 100*l/k, '%'

  !	 pause

  ! arrest friction distribution is output of get_uniform_dtau.m

  tauamx=0



  !       open(1,file='./iofiles/tauaunif')
  !	do i=1,nl
  !	 write(1,*)(taua(i,j),j=1,nd)
  !	enddo
  !	close(1)
  !	pause

  ! to add random flactuations to stress drops set rflag=1 and add .rn
  ! to files last*, s128*
  open(1,file='./iofiles/stressdrops.unif.12pm6.txt')
  do i=1,nl
     read(1,*)(taua(i,j),j=1,nd) ! this is actually dtau, but changed below
  enddo
  close(1)

  do j=1,nd
     do i=1,nl
!        taua(i,j)=taus(i,j)-dtaumx*(p+s*(2*taua(i,j)-1))
        taua(i,j)=taus(i,j)-taua(i,j)
        if(taua(i,j).gt.tauamx)tauamx=taua(i,j)
        if(taua(i,j).lt.0.0)then
           write(*,*) 'ERROR!!: tau_a less than 0.'
           stop -6
        endif
     enddo
  enddo
  write(*,*) 'tauamx=',tauamx


  ! initial stress
  do j=1,nd
     do i=1,nl
        tau01 = taus(i,j)-dtaumx/2.              ! bars
        tau02 = 0.95*(Vpl/crp(i,j))**(1./3.)     ! bars
        !	 tau03 = tau01*tau02/(tau01 + tau02)
        tau0(i,j) = min(tau01, tau02)            ! bars
     enddo
  enddo


  ! FOR FRESH RUN do the following:
  !   1. copy input/output file last*.t* to last* (or vise-versa)
  !      (to find last output-time un-comment "pause" below)
  !   2. rename output file slipdb2* to have apropriate .t* ending
  !   3. un-commment "write(2)" and "goto 13" statements below

  open(1,file=ifilename,form='unformatted')
  open(2,file=ofilename,form='unformatted')
  !c	open(3,file='././iofiles/hypodb2.unif')

  !c	name='../iofiles/dtau.unif'
  !c	open(69,file=name,form='unformatted')

  read(1)upl,u1,duc              ! displacments in mm
  t=upl/Vpl                      ! time in yr
  t1=t                           ! initial time for current run
  write(*,*) 'end time in file last* =',t1,'yr, upl = ',upl/1000.,'m'
  !c	write(*,*) ' '
  !c	pause

  ! FOR 'RESTART' RUNS comment out next 2 line

  write(2)upl,u1
  goto 13

  !-------------------  restart portion of the program --------------

  read(2)var,u2                     ! initial slip in mm

  do i=1,10000000
     read(2,end=12)nhypo,ihypo,jhypo,t1,iduf,duc1
     if(t1.ge.t)then  !last run crashed between write(2) and write(1)
        !             statements. need to read one less record from slip* file
        write(*,*) 'error: restart problem'
        stop -5
     endif
  enddo
12 continue
  write(*,*) 't-t1=',t-t1
  t1=t


  !-----------------------  end of restart portion ------------------

13 continue
  !	pause

  !lock the fault; intialize fault slip u2 immediately after brittle
  !failure (u1 is fault slip immediately before brittle failure)

  do j=1,nd
     do i=1,nl
        tauf(i,j)=taus(i,j)
        u2(i,j)=u1(i,j)
     enddo
  enddo

  !c======================================================================
  !c======================================================================

  do it=1,3500000          ! external loop for model evolution

     !c	write(*,*) '---------------'
     !c	write(*,*) 'model run time index =',it,', time =',t
     !c	write(*,*) '---------------'

     !c--------------------------------------------------
     !      Initialize Stress drop
     do j=1,nd
        do i=1,nl
           dtau(i,j)=0.0 ! stress drop for each cell
           has_failed(i,j)=.FALSE. ! flag for failure during eqk
           taustart(i,j)=0.0 ! stress before eqk
        enddo
     enddo
     taumn=999.9
     taumx=0.0

     !   calculate stresses
     do j=1,nd      ! obs cell depth index
        do i=1,nl      ! obs cell length index

           tau(i,j)=tau0(i,j)

           ! stress due to interior region
           do l=1,nd     ! source cell depth index
              do k=1,nl     ! source cell length index
                 ik=abs(i-k)+1
                 tau(i,j)=tau(i,j)+fi(ik,j,l)*( u1(k,l)-upl )
              enddo
           enddo

           ! add stresses due to back slip of bounding large exterior cells at
           ! current cycle. (initial fault configuration is given by output of
           ! slipinit.f; this corresponds to 150 yr of plate and computational
           ! grid motions, during which the large cells are locked.)

           ! option one - both exterior large cells fail simultaneously:

           !c	  if(t.ge.150)then
           !c	   tau(i,j)=tau(i,j)+( fes(i,j)+fen(i,j) )*(-Vpl*(t-150))
           !c	  else
           !c	   tau(i,j)=tau(i,j)+( fes(i,j)+fen(i,j) )*(-Vpl*t)
           !c	  endif

           ! option two - exterior large cells do not fail simultaneously;
           ! 1857-type event occur at t=150; 1906-type at t=200:

           if(t.gt.150)then
              tau(i,j)=tau(i,j)+fes(i,j)*(-Vpl*(t-150))
           else
              tau(i,j)=tau(i,j)+fes(i,j)*(-Vpl*t)
           endif

           if(t.lt.200)then
              tau(i,j)=tau(i,j)+fen(i,j)*(-Vpl*t)
           else
              tau(i,j)=tau(i,j)+fen(i,j)*(-Vpl*(t-200))
           endif

        enddo
     enddo

     ! Check for negative stress
     if(tau(i,j).lt.0)then
        write(*,*) 'error: negative stress at i,j=',i,j
        stop -4
     endif

     ! find hypocenter parameters
     nhypo=0
     ihypo=0
     jhypo=0
     delmax=-9999
     do j=1,nd
        do i=1,nl
           if(tau(i,j).ge.taus(i,j))then
              nhypo=nhypo+1
              if(tau(i,j)-taus(i,j).gt.delmax)then
                 ihypo=i                   ! i coord of slip initiator
                 jhypo=j                   ! j coord of slip initiator
                 delmax=tau(i,j)-taus(i,j)
              endif
           endif
        enddo
     enddo
     !	 write(*,*) 'nhypo,ihypo,jhypo= ',nhypo,ihypo,jhypo

     if(nhypo.ne.0)then
        !c	  write(*,*) '***Earthquake!!***'
        xhypo=x(ihypo,jhypo)
        zhypo=z(ihypo,jhypo)
        !c	  write(3,*)xhypo,zhypo,t



        !--------------------------------------------------
        !record the stress before we transfer/propagate the rupture
        do j=1,nd
           do i=1,nl
              taustart(i,j)=tau(i,j)
           enddo
        enddo

        !--------------------------------------------------

        ! test for any cells that are above failure strength
1       indf=0   ! number of current failures

        ! Return to here if we are in an earthquake iteration

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
        !c	 write(*,*) 'indf=',indf

        ! adjust stresses and brittle displacements; final stresses after brittle
        ! failure will be distributed between arrest stress and dynamic friction
        do m=1,indf
           i=ifail(m)
           j=jfail(m)
           has_failed(i,j)=.TRUE. ! *** Record that a cell has failed
           do jj=1,nd
              do ii=1,nl
                 ik=abs(ii-i)+1
                 tau(ii,jj)=tau(ii,jj)+fi(ik,jj,j)*delu(i,j)
              enddo
           enddo
           u2(i,j)=u2(i,j)+delu(i,j)
        enddo

        if(indf.gt.0)goto 1

        !--------------------------------------------------
        !record the static stress drop for each slipping cell
        sumdtau=0.0; ! for computing cell avg stress drop
        Nslip=0;
        do j=1,nd
           do i=1,nl
              if(has_failed(i,j))then
                 ! c	      if((u2(i,j)-u1(i,j)).gt.0)then
                 dtau(i,j)=taustart(i,j)-tau(i,j)
                 sumdtau=sumdtau+dtau(i,j)
                 Nslip=Nslip+1


                 if(dtau(i,j).gt.taumx)then
                    taumx=dtau(i,j)
                    taumx_i=i
                    taumx_j=j
                 endif
                 if(dtau(i,j).lt.taumn)then
                    taumn=dtau(i,j)
                    taumn_i=i
                    taumn_j=j
                 endif
              endif
           enddo
        enddo

        !	write(*,*)  '*** EQK at t = ',t
        !	write(*,*) 'Avg dtau:',sumdtau/Nslip,', Ncells = ',Nslip
        !	if(Nslip.gt.1)then
        !	   write(*,*) 'Max dtau:',taumx,'at',taumx_i,taumx_j
        !	   write(*,*) 'Min dtau:',taumn,'at',taumn_i,taumn_j
        !	endif

        write(*,'(F8.4,1X,I4,1X,F12.4,1X,F12.4,1X,I3,1X,I3,1X, F12.4,1X,I3,1X,I3)') &
             t,Nslip,sumdtau/Nslip,taumx,taumx_i,taumx_j,taumn &
             ,taumn_i,taumn_j


        !--------------------------------------------------
        ! calc int(delu) of previous brittle slip increment for output
        do j=1,nd
           do i=1,nl
              delu(i,j)=u2(i,j)-u1(i,j)
              diff=nint( 10.*delu(i,j) )            ! delu in 0.1mm
              if(diff.gt.32000)diff=32000-diff      ! max(int*2)=32767
              if(abs(diff).gt.32000)then
                 write(*,*) 'error: delu overflow'
                 stop -3
              endif
              iduf(i,j) = diff
           enddo
        enddo

        ! update fault motion to account for last brittle slip; lock the fault
        do j=1,nd
           do i=1,nl
              u1(i,j)=u2(i,j)
              tauf(i,j)=taus(i,j)
           enddo
        enddo


     endif ! Continue from here if there are no failure points


     ! calc vcrp(i,j) for next time step from current stresses

     do j=1,nd
        do i=1,nl
           vcrp(i,j) = crp(i,j)*(tau(i,j)**3)
        enddo
     enddo

     ! Return to here
3    continue

     ! update fault motion to account for creep motion over assumed next
     ! time step
     do j=1,nd
        do i=1,nl

           ! Exit if total brittle slip exceeds total plate motion
           if(u1(i,j).gt.upl)then
              write(*,*) 'error: brittle slip overshoot at i,j=',i,j
              stop -1
           endif
           u1(i,j)=u1(i,j)+vcrp(i,j)*dt

           ! Exit if total creep slip exceeds total plate motion
           if(u1(i,j).gt.upl)then
              write(*,*) 'error: creep overshoot at i,j=',i,j
              stop -2
           endif

        enddo
     enddo

     ! update time and plate motion at end of assumed time step
     t = t+dt
     upl = t*Vpl

!     if(t.ge.200.and.t.lt.200+dtmx/2.)goto 4
     if(t.lt.200.or.t.ge.200+dtmx/2.) then

     ! check that new nhypo < 2

     !   calculate stresses

     do j=1,nd      ! obs cell depth index
        do i=1,nl      ! obs cell length index

           tau(i,j)=tau0(i,j)

           ! stress due to interior region

           do l=1,nd     ! source cell depth index
              do k=1,nl     ! source cell length index
                 ik=abs(i-k)+1
                 tau(i,j)=tau(i,j)+fi(ik,j,l)*( u1(k,l)-upl )
              enddo
           enddo

           ! add stresses due to back slip of bounding large exterior cells

           if(t.gt.150)then
              tau(i,j)=tau(i,j)+fes(i,j)*(-Vpl*(t-150))
           else
              tau(i,j)=tau(i,j)+fes(i,j)*(-Vpl*t)
           endif

           if(t.lt.200)then
              tau(i,j)=tau(i,j)+fen(i,j)*(-Vpl*t)
           else
              tau(i,j)=tau(i,j)+fen(i,j)*(-Vpl*(t-200))
           endif

        enddo
     enddo

     ! find tentative hypocenter parameters for next time step

     nhypo1=0
     do j=1,nd
        do i=1,nl
           if(tau(i,j).ge.taus(i,j))nhypo1=nhypo1+1
        enddo
     enddo
     !c	 write(*,*) 'dt,nhypo1=',dt,nhypo1

     if(nhypo1.gt.1)then    ! repeat calculations with dt/2

        ! subtract assumed last time step from time
        ! subtract assumed last creep motion from fault displacement
        ! reduce size of time step; re-calculate creep & plate motions

        t=t-dt

        do j=1,nd
           do i=1,nl
              u1(i,j)=u1(i,j)-vcrp(i,j)*dt
           enddo
        enddo

        dt=dt/2.

        ! Return to creep part
        goto 3

     endif

     ! Continue from here if we have reached t=200
!4    continue
  end if

     ! accumulate/write output

     if(nhypo.eq.0)then      ! no earthquake; accumulate deluc

        do j=1,nd
           do i=1,nl
              duc(i,j)=duc(i,j)+vcrp(i,j)*dt
           enddo
        enddo

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
        write(2)nhypo,ihypo,jhypo,t-dt,iduf,duc,dtau
!        call flush(2)

        ! re-intialize deluc with creep-sli taustart(i,j)p after last earthquake

        do j=1,nd
           do i=1,nl
              duc(i,j)=vcrp(i,j)*dt
           enddo
        enddo

        ! write displacements at start of next time step in 'restart' file

        rewind 1
        write(1)upl,u1,duc
!        call flush (1)

     endif

     ! reset time step size and fault slip u2

     dt=dtmx

     do j=1,nd
        do i=1,nl
           u2(i,j)=u1(i,j)
        enddo
     enddo

     !  check exit criteria

     !c	 if(t.ge.400)goto 99          ! 250 yr of model evolution
     !c	 if(t-t1.ge.10)goto 99        ! 10 yr of model evolution
     !c	 if(t-t1.ge.1)goto 99        ! 1 yr of model evolution
     if(t-t1.ge.150) exit        ! 150 yr of model evolution


  enddo    !  on it

  close(2)
  !c	close(3)

  close(1)

  stop
end program

!-------------------------------------------------------------

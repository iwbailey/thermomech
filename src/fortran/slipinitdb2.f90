program faultslip
  use stiffnessmatrix
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
  real, parameter :: tauratiox = 4
  real, parameter :: zDB = 10 ! depth where creep rate at tau = fs*Seff equals Vpl
  real, parameter :: xDB = 62.5 ! analogous horizontal position of DB transition
  real, parameter :: t0=125.0

  real(kind=8), dimension(nl,nd,nd) :: fi
  real(kind=8), dimension(nl,nd) :: fes, fen
  real(kind=8), dimension(nl,nd) :: u1, u2, tau, taus, taua, tauf, tau0, delu, &
       x, z, crp, vcrp, sk, du
  integer, dimension(nl*nd) :: ifail, jfail

  integer i, j, l, k, m, ii, jj, ik, it, ihypo, jhypo, nhypo, ihypo1, jhypo1, &
       indf, nEQ
  integer idum

  character(200) ofilenameStressDrops, ofilename2, ofilenameStrength, ofilename4, ifilename
  real(kind=8) rn, ran1

  real(kind=8) delmax, xd2, xw2, hcell, hxcell, temp, temp1, temp2, upl, xcrpdb, xhypo
  real(kind=8) xrcrp, porepress, porepressDB, porepressDBx, Seff, SeffDB, SeffDBx, sigmanormal
  real(kind=8) sigmanormalDB, sigmanormalDBx, tau01, tau02, tauamx, xi, zcrpDB, zDBx
  real(kind=8) zhypo, zj, zrcrp, t


  xd2=Zdepth/(2.*nd) ! cell half-depth
  xw2=Xlength/(2.*nl) ! cell half-width

  ifilename = './iofiles/stressdrops.unif.12pm6.txt'
  ofilenameStressDrops ='./iofiles/stressdrops.unif.dtau12pm6.out'
  ofilename2 ='./iofiles/lastdb2.unif.t150.dtau12pm6'
  ofilenameStrength ='./iofiles/strength.unif.dtau12pm6.out'
  ofilename4 ='./iofiles/creepcoef.unif.dtau12pm6.out'

  !  Define ff(i-k,j,l) by tau(i,j)=[Sum over k,l]ff(i-k,j,l)*slip(k,l);
  !  here i,k=1,nl and j,l=1,nd.  Shear modulus is 300 kbars.
  !  Function ff(i-k,j,l) is equated to f(ik,j,l) where ik = ile1,nl
  call get_stiffnessi( fi, nl, nd, Xlength, Zdepth )
  call get_stiffnesse( fes, nl, nd, Xlength, Zdepth, Xlength2, Zdepth, .true., &
       0.0 )
  call get_stiffnesse( fen, nl, nd, Xlength, Zdepth, Xlength2, Zdepth, .false., &
       0.0 )

  ! calculate sum of stiffnesses (effective stiffness for uniform forward
  ! fault motion) for each cell
  do j = 1,nd
     do i = 1,nl
        sk(i,j) = fes(i,j)+fen(i,j)
        do l = 1,nd
           do k = 1,nl
              ik = abs(i-k)+1
              sk(i,j) = sk(i,j)+fi(ik,j,l)
           enddo
        enddo
        !       write(*,*)'i,j,sk(i,j)',i,j,sk(i,j)
     enddo
  enddo

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

  ! Brittle strength (static friction), and creep properties:
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

        taus(k,l)=dtaumx+fs*Seff          ! static brittle strength

        !	  write(*,*)taus(k,l)

        !  Creep parameters: depth and poition along strike zDB and xDB of
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
        !

        zcrpDB = Vpl/((fs*SeffDB)**3) ! Makes Vcreep = Vpl at zDB
        !                                      if tau = fs*Seff at depth zDB.
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
  l = 0
  k=0
  do j=1,nd
     do i=1,nl
        if(z(i,j).lt.13.75.and.x(i,j).lt.66.25) then
           !rn = ran1(idum)
           rn = rand()
           k = k+1
           if(rn.le.0.2)then
              l = l + 1
              crp(i,j)=temp2
              if(z(i,j).lt.zDB.and.x(i,j).lt.xDB)crp(i,j)=temp1
           endif
        endif
     enddo
  enddo
  write(*,*) l, 'random selections, i.e.,', 100*l/k, '%'
  !	 pause

  !       pause

  ! arrest friction distribution is output of get_uniform_dtau.m

  tauamx=0

  open ( 10, file=ifilename, status='old', iostat=k, action='read')
  write(*,*) 'Reading ', ifilename, k
  do i=1,nl
     read(10,*)(taua(i,j),j=1,nd) ! this is actually dtau, but changed below
  enddo
  close(10)

  do j=1,nd
     do i=1,nl
!        taua(i,j)=taus(i,j)-dtauavg+pm*(2*taua(i,j)-1)
        taua(i,j)=taus(i,j)-taua(i,j)
        if(taua(i,j).gt.tauamx)tauamx=taua(i,j)
     enddo
  enddo
  write(*,*)'tauamx=',tauamx

  !      write the distribution of cell values
  open(1,file=ofilenameStressDrops) ! strength - arrest stress
  open(3,file=ofilenameStrength) ! strength
  open(4,file=ofilename4) ! creep values

  do i=1,nl
     do j=1,nd
        write(1,*)taus(i,j)-taua(i,j)
        write(3,*)taus(i,j)
        write(4,*)crp(i,j)
     enddo
  enddo
  close(1)
  write(*,*)'Written to ', ofilenameStressDrops
  close(3)
  write(*,*)'Written to ', ofilenameStrength
  close(4)
  write(*,*)'Written to ', ofilename4
  !	pause

  ! initial stress
  do i=1,nl
     do j=1,nd
        tau01 = taus(i,j)-dtaumx/2.              ! bars
        tau02 = 0.95*(Vpl/crp(i,j))**(1./3.)     ! bars
        !	 tau03 = tau01*tau02/(tau01 + tau02)
        tau0(i,j) = min(tau01, tau02)            ! bars
     enddo
  enddo


  ! lock the fault; initialize fault offset from plate motion u1, u2
  do j=1,nd
     do i=1,nl
        tauf(i,j)=taus(i,j)
        u1(i,j)=0
        u2(i,j)=u1(i,j)
     enddo
  enddo

  !	open(2,file='./iofiles/tau2.unif',form='unformatted')
  !	open(3,file='./iofiles/hypo2.unif')

  ! impose 125 yr of plate motion
  t = t0
  nEQ = 0
  upl=t*vpl

  !	do it=1,25*356/7+1        ! external loop for model evolution
  do it=1,25*356/4+1        ! external loop for model evolution
     !	do it=1,25*365/3+1        ! external loop for model evolution

     !	write(*,*)'---------------'
     write(*,*)'model run time index =',it,', time = ',t
     !	write(*,*)'---------------'

     !   calculate stresses

     do j=1,nd      ! obs cell depth index
        do i=1,nl      ! obs cell length index

           ! stress due to back slip of bounding large exterior cells at
           ! current cycle.

           tau(i,j)=tau0(i,j)+(fes(i,j)+fen(i,j))*(-vpl*t)

           ! stress due to interior region

           do l=1,nd     ! source cell depth index
              do k=1,nl     ! source cell length index
                 ik=abs(i-k)+1
                 tau(i,j)=tau(i,j)+fi(ik,j,l)*(u1(k,l)-upl)
              enddo
           enddo

        enddo
     enddo


     do j=1,nd
        do i=1,nl
           if(tau(i,j).lt.0)then
              write(*,*)'negative stress at i,j=',i,j
              !         pause
           endif
        enddo
     enddo

     ! find location of hypocenter
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
                 !	    write(*,*)'i,j,taus,tau are',i,j,taus(i,j),tau(i,j)
              endif
           endif
        enddo
     enddo
     if(nhypo.ge.1)then
        write(*,*)'nhypo= ',nhypo
        write(*,*)'ihypo,ihypo1,jhypo,jhypo1= ',ihypo,ihypo1,jhypo,jhypo1
        nEQ = nEQ+1

     endif
     if(ihypo.ne.0)then
        xhypo=x(ihypo,jhypo)
        zhypo=z(ihypo,jhypo)
        !	  write(3,*)t,xhypo,zhypo
     endif

     ! test for brittle failure

1     indf=0   ! number of current failures

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


     ! adjust stresses and brittle displacements; final stresses after
     ! brittle failure will be distributed between arrest stress and
     !  dynamic friction

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
        !	  write(*,*)'i,j,slip,tau=',i,j,delu(i,j),tau(i,j)
     enddo

     if(indf.gt.0) then
        write(*,*)'indf=',indf
        goto 1
     endif

     ! calc vcrp(i,j) for next time step; lock the fault

     do j=1,nd
        do i=1,nl
           vcrp(i,j) = crp(i,j)*(tau(i,j)**3)
           tauf(i,j)=taus(i,j)
        enddo
     enddo

     !	 do j=1,nd
     !        write(*,*)'vcrp(128,j)=',vcrp(128,j)
     !	 enddo
     !	 pause

     ! update fault motion to account for last brittle slip and next
     ! creep motion; update time and plate motion over next time step;
     ! check exit criteria

     do j=1,nd
        do i=1,nl
           u1(i,j)=u2(i,j)+vcrp(i,j)*dt
           if(u1(i,j).gt.upl)then
              write(*,*)'overshoot at i,j=',i,j
              u1(i,j)=upl-1
           endif
           u2(i,j)=u1(i,j)
        enddo
     enddo
     t=t+dt
     upl=t*vpl

  enddo    !  on it

  write(*,*)nEQ, ' earthquake(s)'

  !	close(3)
  open(2,file=ofilename2,form='unformatted')
  write(2)upl,u1,du    ! data in mm
  !	write(2)t,tau
  close(2)
  write(*,*)'Written to ', ofilename2

  !	close(2)
   ! TMP debugging
    open( 10, file='tmp.out')
    do i=1,nl
       do j=1,nd
          write( 10, *) u1(i,j)
       end do
    end do
    stop

  stop
end program faultslip

!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!===============================================================================

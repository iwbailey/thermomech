module stiffnessmatrix
  ! ====================== Specification part ===============
  real(kind=8), parameter :: xmu = 0.3 ! shear modulus in Mbars
  real(kind=8), parameter :: PI = 3.141592653589793239

contains
  ! ====================== Implementation part ===============

  subroutine get_stiffnessi( f, nl, nd, Xlength, Zdepth )
    !   This program calculates stiffnesses for uniformly gridded vertical fault
    !   in a half-space, due to interior slip. Output consists of stiffness
    !   components fi(i,j,l) for i=1,Nlength j,l=1,Ndepth.
    !
    !   The elastic stiffness -sshr(i,j) is defined such that
    !   - dtau(i) = [Sum on j] sshr(j,i)*[vpl*dt-du(j)];  u(j) = slip
    !
    !    Matrix -sshr(j,i) is transpose of the more usual definition of
    !    stiffness.  Also, it refers to stresses, not work-conjugate forces, and
    !    hence need not always have the expected symmetry properties.
    !
    !   dimension needs: f(Nlength, Ndepth, Ndepth);

    implicit none

    ! Sub-routine arguments
    integer, intent(in) :: nl, nd ! dimension of fault
    real(kind=8), intent(in) :: Xlength, Zdepth ! Length and width of fault in km
    real(kind=8), intent(out), dimension(nl,nd,nd) :: f ! output matrix

    ! Variables used in the function
    integer i, j, l
    real(kind=8) xd2, xw2, x3, y1, y3, ri, rj, rl
    real(kind=8) e12

    ! Cell half width and depth
    xd2 = 0.5*Zdepth/(nd)
    xw2 = 0.5*Xlength/(nl)

    ! Stress Green's function (stiffness):
    !
    ! Define ff(i-k,j,l) by tau(i,j)=[Sum over k,l]ff(i-k,j,l)*slip(k,l); here
    ! i,k=1,Nlength and j,l=1,Ndepth.
    !
    ! Function ff(i-k,j,l) is equated to f(ik,j,l) where ik = 1,Nlength

    ! Loop through depth index for slip source cell
    do l = 1,nd

       ! Get depth (cell center) of slip source
       rl = l
       x3 = (rl-0.5)*2*xd2

       ! Loop through depth index of stress calc cell
       do j=1,nd

          ! Get depth for stress calculation
          rj = j
          y3 = (rj-0.5)*2*xd2

          ! Loop through strike indices
          do i=1,nl

             ! Get horizontal distance between slip source and stress calc
             ri = i
             y1 = (ri-1.0)*2*xw2 ! dist along strike for stress calculation

             ! Get the value for this combination
             e12 = strn (xw2,xd2,x3,y1,y3)
             f(i,j,l)=2*xmu*e12
             !	      type*,f(i,j,l)
          end do
       end do
    end do

  end subroutine get_stiffnessi

  !---------------------------------------------------------------------
  subroutine get_stiffnesse( f, nl, nd, Xlength, Zdepth, Xlength2, Zdepth2, isSouth, sepdist )

  !   This program calculates stiffnesses for uniformly gridded vertical
  !   fault in a half-space, due to slip in an exterior large SOUTHERN
  !   bounding cell . Provides output to be read by friction_fault.
  !   Output consists of stiffness components fes(i,j) for i=1,Nlength
  !   j=1,Ndepth.

  !   ne=Nlength(=128)*Ndepth(=32)  ! elements of interior region

  !   The elastic stiffness -sshr(i,j) is defined such that
  !   - dtau(i) = [Sum on j] sshr(j,i)*[vpl*dt-du(j)];  u(j) = slip

  !    Matrix -sshr(j,i) is
  !    transpose of the more usual definition of stiffness.  Also, it
  !    refers to stresses, not work-conjugate forces,  and hence need
  !    not always have the expected symmetry properties.

    implicit none

    ! Sub-routine arguments
    integer, intent(in) :: nl, nd ! dimension of fault
    real(kind=8), intent(in) :: Xlength, Zdepth, Xlength2, Zdepth2 ! Length and width of fault in km
    logical, intent(in) :: isSouth ! Flag for whether the external fault is South
    real(kind=8), intent(in) :: sepdist ! distance between end of this fault and external fault
    real(kind=8), intent(out), dimension(nl,nd) :: f ! output matrix

    ! Internal variables
    integer i, j
    real(kind=8) xd2, xw2, yd2, yw2, y1, y3, ri, rj, signNS, dx
    real(kind=8) e12

    ! Check the separation distance is correct
    if( sepdist < 0 ) stop -1

    ! Cell half-length and halfwidth
    yd2 = 0.5*Zdepth/nd     ! interior obs cell half-depth
    yw2 = 0.5*Xlength/nl   ! interior obs cell half-width

    ! External section half-length and width
    xd2 = 0.5*Zdepth2
    xw2 = 0.5*Xlength2

    if( isSouth ) then
       dx = sepdist
       signNS = 1.0
    else
       ! If to the north add on the length of the fault
       dx = sepdist + Xlength
       signNS = -1.0
    end if

    ! Stress Green's function (stiffness):
    ! Define ff(i,j) by tau(i,j)=ff(i,j)*slip(exterior source cell);

    ! Loop through depth index of stress calc cell
    do j = 1,nd       ! j=depth index for stress calculation cell

       ! Get depth for stress calculation
       rj = j
       y3 = (rj-0.5)*2.0*yd2

       ! Loop through strike index of stress calculation cell
       do i = 1,nl    ! i = number of interior cells along strike

          ! Calc separation distance, generalised from...
          ! y1=xw2+170-(ri-0.5)*2*yw2 ... for northern
          ! i.e. xw2 + 170 { -yw2, -3*yw2, -5*yw2, ...}
          !
          ! y1=xw2+(ri-0.5)*2*yw2 ... for southern
          ! i.e. xw2 + { yw2, 3*yw2, 5*yw2, ...}
          ri = i
          y1 = xw2 + dx + signNS*2.0*yw2*(ri-0.5)

          ! Calculate the strain
          e12 = strn ( xw2, xd2, xd2, y1, y3)

          ! Record in the matrix
          f(i,j) = 2*xmu*e12
     end do
  end do


  end subroutine get_stiffnesse

  !---------------------------------------------------------------------
  function strn ( xw2, xd2, x3, y1, y3)
    !     Chinnery solution for stress due to strike slip on a vertical,
    !     rectangular fault patch  (following materials modified from an
    !     original coding by W. D. Stuart provided Oct.'90 by T. E. Tullis)
    !
    !     single vertical flat strike slip fault plane
    !
    !     units: slip,displ;mm   dist,km   elas const,Mbars
    !            stress,bars     strain,microstrain
    !
    !         fault slip  delta=u(+)-u(-)  left lateral is + here
    !     x,y,z    (x1,x2,x3) x1 +rt,par flt  x2 nor flt  x3 dwn
    !     xc,zc    coords cell centers, km
    !     xmu      rigidity, poisson solid     in Mbars
    !
    !       for unit lft lat slip in a given cell, find shear stress
    !       at some other (or that) cell center
    !
    !       chinnery, bssa, 1963   y1 par flt  y2 nor flt  y3 dwn
    !       analytic soln pers. com. 1983  xlam=xmu  x1=x2=y2=0
    !       tensor strain e12 on fault plane
    !       x1(=0),x3 (source cell center); y1,y3 (stress calculation point)

    implicit none
    real(kind=8) :: strn

    ! Function arguments
    ! xw2,xd2 are half width, half height
    real(kind=8), intent(in) :: xw2, xd2, x3, y1, y3

    ! Function variables
    real(kind=8) e12
    real(kind=8), dimension(4) :: sgn, fx1, fx3
    integer i
    real(kind=8) t, q, p, tt, qq, pp, s1, s2, s1q, s2p, f1, f2, f4, etmp

    ! Initialize arrays
    sgn = (/1.0, -1.0, -1.0, 1.0/)
    fx1 = (/ +xw2, +xw2, -xw2, -xw2 /)
    fx3 = (/ x3+xd2, x3-xd2, x3+xd2, x3-xd2 /)
    e12 = 0.0

    ! Loop through all four componenents of the summation
    do i = 1,4
       t = fx1(i)-y1
       q = fx3(i)-y3
       p = fx3(i)+y3
       tt = t*t
       qq = q*q
       pp = p*p
       s1 = sqrt(tt+qq)
       s2 = sqrt(tt+pp)
       s1q = s1+q
       s2p = s2+p

       f1 = 1./(s1*s1q) + 1./(s2*s2p)
       f2 = (s2/4.+q)/(s2*s2p*s2p) - (pp-qq)*(2.*s2+p)/(2.*s2**3*s2p*s2p)
       f4 = q/(s1*(s1+t)) + p/(s2*(s2+t))

       etmp = (2.0/3.0) * t * (f1 + f2) + 0.5*f4

       e12  = e12 + sgn(i) * (0.25/PI) * etmp
    end do

    ! Get the return value
    strn = e12
  end function strn
  !-------------------------------------------------------------------------------
  !===============================================================================
end module stiffnessmatrix

module initialize
  contains
    !---------------------------------------------------------------------------
    function faulttaus( nl, nd, dz, tauc, fs, dsigma_dz)
      ! Set up the fault array of static strength
      !
      ! IN:
      ! nl, nd = dimensions along stike, down dip
      ! dx = cell height
      ! tauc = static strength at zero depth
      ! fs = static coefficient of friction
      ! dsigma_dz = gradient of effective normal stress as a function of depth
      !
      ! OUT:
      ! faulttaus = (nl, nd) array of static strength for each cell on the
      ! fault. Each column is the same
      !
      implicit none

      ! Input arguments
      integer, intent(in):: nl, nd
      real, intent(in) :: dz, tauc, fs, dsigma_dz

      ! Output
      real(kind=8), dimension(nl, nd) :: faulttaus

      ! Variables in the function
      integer i
      integer, dimension(nd) :: iDepth
      real(kind=8), dimension(nd) :: taus, zcell

      ! Calculate depths at the center of each cell
      iDepth = (/ (i, i = 1,nd) /)
      zcell = (iDepth-0.5)*dz

      ! Calculate the strength profile
      taus = tauc + fs*dsigma_dz*zcell

      ! Repeat the strength profile along all elements
      do i = 1,nl
         faulttaus(i,:) = taus
      end do

    end function faulttaus
    !---------------------------------------------------------------------------
    function faultcrp( nl, nd, dx, dz, Vpl, fs, SeffDB, tauratio, zDB, xDB)
      !  Creep parameters: depth and poition along strike zDB and xDB of
      ! "ductile-brittle" transitions, SeffDB, and SeffDBx given above.  We
      ! assume Vcreep = crp(x,z)*tau(x,z)**3, where the spatial creep property
      ! distribution crp(x,z) is given either by crp(x,z) =
      ! xcrpDB*exp(xrcrp*(x-xzDB)) + zcrpDB*exp(zrcrp*(z-zDB)) or by the maximum
      ! of two analogous 1D functions. xrcrp and zrcrp are chosen to enforce
      ! given stress ratios (say, 10.0) between that at zDB and that at Zdepth
      ! (or xDB and Xlength) if the same Vcreep is to be produced at both depth
      ! (or horizontal locations).  choices of Zdepth and Xlength correspond to
      ! depth and along-strike positions where slip along the SAF is mostly in
      ! the form of creep.

      implicit none
      ! Input arguments
      integer, intent(in):: nl, nd
      real, intent(in) :: dx, dz, Vpl, fs, tauratio, xDB, zDB
      real(kind=8), intent(in) :: SeffDB

      ! Output
      real(kind=8), dimension(nl, nd) :: faultcrp

      ! Variables used in the function
      integer i, j
      real faultLength, faultDepth
      real(kind=8) zcrpDB, zrcrp, xrcrp
      integer, dimension(nd) :: iDepth
      integer, dimension(nl) :: iStrike
      real(kind=8), dimension(nd) :: crpProfile, zcell
      real(kind=8), dimension(nl) :: crpXsect, xcell

      ! Calculate depths at the center of each cell
      iDepth = (/ (i, i = 1,nd) /)
      zcell = (iDepth-0.5)*dz
      iStrike = (/ (i, i = 1,nl) /)
      xcell = (iStrike-0.5)*dx
      faultDepth = nd*dz
      faultLength = nl*dx

      ! Makes Vcreep = Vpl at zDB if tau = fs*Seff at depth zDB.
      zcrpDB = Vpl/((fs*SeffDB)**3)

      zrcrp = 3.0*log(tauratio)/(faultDepth - zDB)
      xrcrp = 3.0*log(tauratio)/xDB

      ! Get along strike and down dip values
      crpProfile = zcrpDB*exp(zrcrp*(zcell - zDB))
      crpXsect   = zcrpDB*exp(xrcrp*max(xDB - xcell, xcell + xDB -faultLength))

      ! Build the array based on the max of each
      do i = 1,nl
         do j=1,nd
            faultcrp(i,j) = max( crpProfile(j), crpXsect(i) )
         end do
      end do

    end function faultcrp
    !---------------------------------------------------------------------------
    function faulttaua( ifilename, nl, nd, taus)
      implicit none
      ! Input arguments
      integer, intent(in):: nl, nd
      real(kind=8), intent(in), dimension(nl,nd) :: taus

      ! Output
      real(kind=8), dimension(nl, nd) :: faulttaua

      integer i, j, k
      character(200) ifilename

      ! arrest friction distribution is output of get_uniform_dtau.m
      open ( 10, file=ifilename, status='old', iostat=k, action='read')
      write(*,*) 'Reading ', ifilename
      do i=1,nl
         read(10,*)(faulttaua(i,j),j=1,nd) ! this is actually dtau, but changed below
      enddo
      close(10)

      ! Convert from dtau to taua
      faulttaua = taus - faulttaua
      write(*,*)'tauamx=',maxval(faulttaua)

    end function faulttaua

  end module initialize

subroutine add_randomness( crp, nl, nd, dx, dz, xDB, zDB, zcrpDB, zrcrp )
  implicit none

  integer, intent(in) :: nl, nd
  real(kind=8), intent(in) :: dx, dz, xDB, zDB
  real(kind=8), intent(in) :: zcrpDB, zrcrp
  real(kind=8), dimension(nl, nd) :: crp

  ! assign to 20% randomly chosen 'brittle' cells higher crp(i,j)
  real(kind=8) faultDepth, faultLength, rn
  integer l, k, i, j

  integer, dimension(nd) :: iDepth
  integer, dimension(nl) :: iStrike
  real(kind=8), dimension(nd) :: zcell
  real(kind=8), dimension(nl) :: xcell

  ! Set the random seed
  call srand(12345)

  ! Calculate depths at the center of each cell
  iDepth = (/ (i, i = 1,nd) /)
  zcell = (iDepth-0.5)*dz
  iStrike = (/ (i, i = 1,nl) /)
  xcell = (iStrike-0.5)*dx
  faultDepth = nd*dz
  faultLength = nl*dx

  ! Set the counters to zero
  l = 0
  k = 0

  do j=1,nd
     if( zcell(j).lt.13.5) then
        do i=1,nl
           if( xcell(i).gt.3.75.and.xcell(i).lt.66.25) then
              rn = rand()
              k = k+1
              if(rn.le.0.2)then
                 l = l + 1
                 if( zcell(j).lt.zDB .and. xcell(i).gt.xDB .and. xcell(i).lt.(faultLength-xDB) ) then
                    crp(i,j) = zcrpDB*exp(zrcrp*(11.25-zDB))
                 else
                    crp(i,j) = zcrpDB*exp(zrcrp*(8.75-zDB))
                 end if
              endif
           endif
        enddo
     end if
  enddo
  write(*,*) l, 'random selections, i.e.,', 100*l/k, '%'

end subroutine add_randomness

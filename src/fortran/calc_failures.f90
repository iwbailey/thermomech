subroutine calc_failures( tau, u2, nl, nd, fi, taus, taud, taua )
  implicit none
  integer, intent(in) :: nl, nd
  real(kind=8), intent(inout), dimension(nl,nd) :: tau, u2
  real(kind=8), intent(in), dimension(nl,nd,nd) :: fi
  real(kind=8), intent(in), dimension(nl,nd) :: taus, taud, taua

  integer indf, i, j, ii, jj, ik, m
  integer, dimension(nl*nd) :: ifail, jfail
  real(kind=8), dimension(nl, nd) ::  delu, tauf

  ! Failure stress starts at taus
  tauf = taus

  indf = 1 ! counter for number of failures
  do while ( indf.gt.0 )

     ! Count failures and compute slip in each cell
     indf=0

     do j=1,nd
        do i=1,nl
           ! Check for failure
           if( tau(i,j)-tauf(i,j).ge.0 )then

              ! Store the failure locations
              indf=indf+1
              ifail(indf)=i
              jfail(indf)=j

              ! Store the slip
              delu(i,j) = ( taua(i,j) - tau(i,j) )/fi(1,j,j) ! brittle self-slip

              ! reduce tauf to dynamic friction
              tauf(i,j) = taud(i,j)

              ! Update the total slip for the cell
              u2(i,j) = u2(i,j)+delu(i,j)

           endif
        enddo
     enddo

     ! Report the number of failures
     if(indf.gt.0) then
        write(*,*)'# failed cells =',indf
     endif

     ! Change stress in all cells as a result of slip
     do m = 1,indf
        do jj=1,nd
           do ii=1,nl
              ik=abs(ii-ifail(m))+1
              tau(ii,jj) = tau(ii,jj) + fi(ik,jj,jfail(m))*delu(ifail(m),jfail(m))
           enddo
        enddo
     end do


    !   ! TMP debugging
    !  open( 10, file='tmp2.out')
    !  do i=1,nl
    !     do j=1,nd
    !        write( 10, *) tau(i,j)
    !     end do
    !  end do
    ! stop

  end do ! end iterations over failures

  ! healing
  tauf = taus

end subroutine calc_failures

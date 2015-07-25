subroutine slipdef2stress( tau, tau0, slipdef, fi, nl, nd )
  implicit none
  integer, intent(in) ::  nl, nd
  real(kind=8), intent(in), dimension(nl,nd) :: tau0, slipdef
  real(kind=8), intent(in), dimension(nl,nd,nd) :: fi
  real(kind=8), dimension(nl,nd) :: tau

  integer i,j,k,l,ik

  tau = tau0
   do j=1,nd      ! obs cell depth index
        do i=1,nl      ! obs cell length index

           ! stress due to interior region
           do l=1,nd     ! source cell depth index
              do k=1,nl     ! source cell length index
                 ik=abs(i-k)+1
                 tau(i,j)=tau(i,j)+fi(ik,j,l)*slipdef(k,l)
              enddo
           enddo

        enddo
     enddo

end subroutine slipdef2stress

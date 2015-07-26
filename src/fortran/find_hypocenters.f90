subroutine find_hypocenters( nhypo, ihypo, jhypo, tau, taus, nl, nd)
  implicit none
  integer nhypo, ihypo, jhypo
  integer, intent(in) :: nl, nd
  real(kind=8), intent(in), dimension(nl,nd) :: tau, taus

  integer, dimension(2) :: ij

  nhypo = count( tau.ge.taus )

  ij = maxloc( tau-taus )
  ihypo = ij(1)
  jhypo = ij(2)

end subroutine find_hypocenters

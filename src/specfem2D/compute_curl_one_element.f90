
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour and CNRS, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and Inria at the following URL
! "http://www.cecill.info".
!
! As a counterpart to the access to the source code and rights to copy,
! modify and redistribute granted by the license, users are provided only
! with a limited warranty and the software's author, the holder of the
! economic rights, and the successive licensors have only limited
! liability.
!
! In this respect, the user's attention is drawn to the risks associated
! with loading, using, modifying and/or developing or reproducing the
! software by the user in light of its specific status of free software,
! that may mean that it is complicated to manipulate, and that also
! therefore means that it is reserved for developers and experienced
! professionals having in-depth computer knowledge. Users are therefore
! encouraged to load and test the software's suitability as regards their
! requirements in conditions enabling the security of their systems and/or
! data to be ensured and, more generally, to use and operate it in the
! same conditions as regards security.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

  subroutine compute_curl_one_element(ispec)

  ! compute curl in (poro)elastic elements (for rotational study)

  use specfem_par, only: curl_element,displ_elastic, &
                         displs_poroelastic,elastic,poroelastic, &
                         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz

  implicit none

  include "constants.h"

  integer ispec

  ! jacobian
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl

  ! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: duz_dxl,dux_dzl
  integer :: i,j,k

  if(elastic(ispec)) then

     do j = 1,NGLLZ
        do i = 1,NGLLX

           ! derivative along x and along z
           dux_dxi = ZERO
           duz_dxi = ZERO

           dux_dgamma = ZERO
           duz_dgamma = ZERO

           ! first double loop over GLL points to compute and store gradients
           ! we can merge the two loops because NGLLX == NGLLZ
           do k = 1,NGLLX
              dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
              duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
              dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
              duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
           enddo

           xixl = xix(i,j,ispec)
           xizl = xiz(i,j,ispec)
           gammaxl = gammax(i,j,ispec)
           gammazl = gammaz(i,j,ispec)

           ! derivatives of displacement
           dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
           duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl

           ! store pressure
           curl_element(i,j) = - 0.5d0 * (dux_dzl - duz_dxl)

        enddo
     enddo

  else if(poroelastic(ispec)) then

     do j = 1,NGLLZ
        do i = 1,NGLLX

           ! derivative along x and along z
           dux_dxi = ZERO
           duz_dxi = ZERO

           dux_dgamma = ZERO
           duz_dgamma = ZERO

           ! first double loop over GLL points to compute and store gradients
           ! we can merge the two loops because NGLLX == NGLLZ
           do k = 1,NGLLX
              dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
              duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
              dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
              duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
           enddo

           xixl = xix(i,j,ispec)
           xizl = xiz(i,j,ispec)
           gammaxl = gammax(i,j,ispec)
           gammazl = gammaz(i,j,ispec)

           ! derivatives of displacement
           dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl
           duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl

           ! store pressure
           curl_element(i,j) = - 0.5d0 * (dux_dzl - duz_dxl)

        enddo
     enddo

  else

     call exit_MPI('no curl in acoustic')

  endif ! end of test if acoustic or elastic element

end subroutine compute_curl_one_element



!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, April 2014
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
 subroutine store_stacey_BC_effect_term_poro()

  use specfem_par, only: p_sv,nspec_left,nspec_right,nspec_bottom,nspec_top, &
                         b_absorb_poro_s_left,b_absorb_poro_w_left, &
                         b_absorb_poro_s_right,b_absorb_poro_w_right, &
                         b_absorb_poro_s_bottom,b_absorb_poro_w_bottom, &
                         b_absorb_poro_s_top,b_absorb_poro_w_top,it
  implicit none
  include "constants.h"
  
  !local variable
  integer :: i,id,ispec

  !--- left absorbing boundary
  if(nspec_left >0) then
    do ispec = 1,nspec_left
      do id =1,2
        do i=1,NGLLZ
          write(45) b_absorb_poro_s_left(id,i,ispec,it)
          write(25) b_absorb_poro_w_left(id,i,ispec,it)
        enddo
      enddo
    enddo
  endif

  !--- right absorbing boundary
  if(nspec_right >0) then
    do ispec = 1,nspec_right
      do id =1,2
        do i=1,NGLLZ
          write(46) b_absorb_poro_s_right(id,i,ispec,it)
          write(26) b_absorb_poro_w_right(id,i,ispec,it)
        enddo
      enddo
    enddo
  endif

  !--- bottom absorbing boundary
  if(nspec_bottom >0) then
    do ispec = 1,nspec_bottom
      do id =1,2
        do i=1,NGLLX
          write(47) b_absorb_poro_s_bottom(id,i,ispec,it)
          write(29) b_absorb_poro_w_bottom(id,i,ispec,it)
        enddo
      enddo
    enddo
  endif

  !--- top absorbing boundary
  if(nspec_top >0) then
    do ispec = 1,nspec_top
      do id =1,2
        do i=1,NGLLX
          write(48) b_absorb_poro_s_top(id,i,ispec,it)
          write(28) b_absorb_poro_w_top(id,i,ispec,it)
        enddo
      enddo
    enddo
  endif

 end subroutine store_stacey_BC_effect_term_poro


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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

  subroutine compute_forces_gravitoacoustic(potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
                                            potential_gravitoacoustic,potential_dot_dot_gravito, &
                                            potential_gravito,IS_BACKWARD_FIELD)

! compute forces for the gravitoacoustic elements

  use constants, only: CUSTOM_REAL

  use specfem_par, only: codeabs,ispec_is_gravitoacoustic,nglob,nspec,nelemabs,numat,it,NSTEP, &
                         anyabs,assign_external_model,ibool,kmato,numabs, &
                         rmass_inverse_gravito, &
                         density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
                         vpext,rhoext,gravityext,Nsqext,hprime_xx,hprimewgll_xx, &
                         hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
                         ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
                         SIMULATION_TYPE,SAVE_FORWARD,nspec_left,nspec_right, &
                         nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top, &
                         b_absorb_acoustic_left,b_absorb_acoustic_right, &
                         b_absorb_acoustic_bottom,b_absorb_acoustic_top

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML

  implicit none

! Chi potential
  real(kind=CUSTOM_REAL), dimension(nglob) :: &
    potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic,potential_gravitoacoustic
! Xi potential
  real(kind=CUSTOM_REAL), dimension(nglob) :: &
    potential_dot_dot_gravito,potential_gravito
! rho*u=grad(Chi)+xi*gravity_vector

  logical,intent(in) :: IS_BACKWARD_FIELD

!---
!--- local variables
!---

  print *, nglob,nspec,nelemabs,numat,it,NSTEP, &
               anyabs,assign_external_model,ibool,kmato,numabs,ispec_is_gravitoacoustic, &
               codeabs,potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
               potential_gravitoacoustic,potential_dot_dot_gravito, &
               !potential_dot_gravito, &
               potential_gravito,rmass_inverse_gravito, &
               !stage_time_scheme, i_stage, &
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,gravityext,Nsqext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
               ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
               SIMULATION_TYPE,SAVE_FORWARD,nspec_left,nspec_right, &
               nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top, &
               b_absorb_acoustic_left,b_absorb_acoustic_right, &
               b_absorb_acoustic_bottom,b_absorb_acoustic_top,IS_BACKWARD_FIELD, &
               ispec_is_PML,PML_BOUNDARY_CONDITIONS

  end subroutine compute_forces_gravitoacoustic


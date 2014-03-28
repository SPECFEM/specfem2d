
!========================================================================
!
!                  S P E C F E M 2 D  Version 7 . 0
!                  --------------------------------
!
! Copyright CNRS, INRIA and University of Pau, France,
! and Princeton University / California Institute of Technology, USA.
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
! license as circulated by CEA, CNRS and INRIA at the following URL
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

  subroutine compute_forces_gravitoacoustic(nglob,nspec,nelemabs,numat,it,NSTEP, &
               anyabs,assign_external_model,ibool,kmato,numabs,gravitoacoustic, &
               codeabs,potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
               potential_gravitoacoustic,potential_dot_dot_gravito, &
               potential_gravito,rmass_inverse_gravito, &
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,gravityext,Nsqext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
               ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
               SIMULATION_TYPE,SAVE_FORWARD,nspec_left,nspec_right,&
               nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top, &
               b_absorb_acoustic_left,b_absorb_acoustic_right, &
               b_absorb_acoustic_bottom,b_absorb_acoustic_top,IS_BACKWARD_FIELD, &
               is_PML,PML_BOUNDARY_CONDITIONS)

! compute forces for the gravitoacoustic elements

  implicit none

  include "constants.h"

  integer :: nglob,nspec,nelemabs,numat,it,NSTEP,SIMULATION_TYPE

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato
  integer, dimension(nelemabs) :: numabs,ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
               ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2

  logical, dimension(nspec) :: gravitoacoustic
  logical, dimension(4,nelemabs)  :: codeabs

! Chi potential
  real(kind=CUSTOM_REAL), dimension(nglob) :: &
    potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic,potential_gravitoacoustic
! Xi potential
  real(kind=CUSTOM_REAL), dimension(nglob) :: &
    potential_dot_dot_gravito,potential_gravito
! rho*u=grad(Chi)+xi*gravity_vector
  real(kind=CUSTOM_REAL), dimension(nglob) :: rmass_inverse_gravito

  double precision, dimension(2,numat) :: density
  double precision, dimension(4,3,numat) :: poroelastcoef
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext,rhoext,gravityext,Nsqext

  logical :: anyabs,assign_external_model
  logical :: SAVE_FORWARD,IS_BACKWARD_FIELD

! derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

! Gauss-Lobatto-Legendre weights
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: wzgll

  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top
  integer, dimension(nelemabs) :: ib_left
  integer, dimension(nelemabs) :: ib_right
  integer, dimension(nelemabs) :: ib_bottom
  integer, dimension(nelemabs) :: ib_top

  real(kind=CUSTOM_REAL), dimension(NGLLZ,nspec_left,NSTEP) :: b_absorb_acoustic_left
  real(kind=CUSTOM_REAL), dimension(NGLLZ,nspec_right,NSTEP) :: b_absorb_acoustic_right
  real(kind=CUSTOM_REAL), dimension(NGLLX,nspec_top,NSTEP) :: b_absorb_acoustic_top
  real(kind=CUSTOM_REAL), dimension(NGLLX,nspec_bottom,NSTEP) :: b_absorb_acoustic_bottom

! CPML coefficients and memory variables
  logical, dimension(nspec) :: is_PML
  logical :: PML_BOUNDARY_CONDITIONS

!---
!--- local variables
!---


  print *, nglob,nspec,nelemabs,numat,it,NSTEP, &
               anyabs,assign_external_model,ibool,kmato,numabs,gravitoacoustic, &
               codeabs,potential_dot_dot_gravitoacoustic,potential_dot_gravitoacoustic, &
               potential_gravitoacoustic,potential_dot_dot_gravito,&
               !potential_dot_gravito, &
               potential_gravito,rmass_inverse_gravito,&
               !stage_time_scheme, i_stage, &
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,gravityext,Nsqext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
               ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
               SIMULATION_TYPE,SAVE_FORWARD,nspec_left,nspec_right,&
               nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top, &
               b_absorb_acoustic_left,b_absorb_acoustic_right, &
               b_absorb_acoustic_bottom,b_absorb_acoustic_top,IS_BACKWARD_FIELD,&
               is_PML,PML_BOUNDARY_CONDITIONS

  end subroutine compute_forces_gravitoacoustic


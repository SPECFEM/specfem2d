!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
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

  subroutine compute_vector_whole_medium(field_acoustic,field_elastic,fields_poroelastic)

! compute Grad(potential) in acoustic elements
! and combine with existing velocity vector field in elastic elements

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM
  use specfem_par, only: nspec,ibool, &
    nglob_acoustic,nglob_elastic,nglob_poroelastic

  use specfem_par_movie, only: vector_field_display

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: field_acoustic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: field_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic) :: fields_poroelastic
  ! vector field in an element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: vector_field_element

  ! local parameters
  integer :: i,j,ispec,iglob

  ! loop over spectral elements
  do ispec = 1,nspec
    ! computes vector field in this element
    call compute_vector_one_element(field_acoustic,field_elastic,fields_poroelastic,ispec,vector_field_element)

    ! stores the result on global nodes
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        vector_field_display(:,iglob) = vector_field_element(:,i,j)
      enddo
    enddo
  enddo

  end subroutine compute_vector_whole_medium

!
!=====================================================================
!

  subroutine compute_vector_one_element(field_acoustic,field_elastic,fields_poroelastic,ispec,vector_field_element)

! compute Grad(potential) if acoustic element or copy existing vector if elastic element

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM

  use specfem_par, only: nglob_acoustic,nglob_elastic,nglob_poroelastic, &
    assign_external_model,density,kmato,rhoext, &
    hprimeBar_xx,hprime_xx,hprime_zz, &
    xix,xiz,gammax,gammaz,ibool, &
    ispec_is_elastic,ispec_is_poroelastic,ispec_is_acoustic, &
    AXISYM,is_on_the_axis, &
    P_SV

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: field_acoustic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_elastic) :: field_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic) :: fields_poroelastic

  integer,intent(in) :: ispec

  ! vector field in an element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ),intent(out) :: vector_field_element

  ! local variables
  integer i,j,k,iglob
  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl
  double precision :: hp1,hp2
  double precision :: rhol
  double precision :: tempx1l,tempx2l

  ! initializes
  vector_field_element(:,:,:) = 0._CUSTOM_REAL

  ! determines vector field
  if (ispec_is_elastic(ispec)) then
    ! elastic element
    ! simple copy of existing vector if elastic element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        if (P_SV) then
          ! P_SV case
          vector_field_element(1,i,j) = field_elastic(1,iglob)
          vector_field_element(2,i,j) = field_elastic(2,iglob)
        else
          ! SH case
          vector_field_element(1,i,j) = field_elastic(1,iglob)
        endif
      enddo
    enddo

  else if (ispec_is_poroelastic(ispec)) then
    ! poro-elastic element
    ! simple copy of existing vector if poroelastic element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        vector_field_element(1,i,j) = fields_poroelastic(1,iglob)
        vector_field_element(2,i,j) = fields_poroelastic(2,iglob)
      enddo
    enddo

  else if (ispec_is_acoustic(ispec)) then
    ! acoustic element

    ! compute gradient of potential to calculate vector if acoustic element
    ! we then need to divide by density because the potential is a potential of (density * displacement)
    rhol = density(1,kmato(ispec))

    ! double loop over GLL points to compute and store gradients
    do j = 1,NGLLZ
      do i = 1,NGLLX
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            do k = 1,NGLLX
              hp1 = hprimeBar_xx(i,k)
              iglob = ibool(k,j,ispec)
              tempx1l = tempx1l + field_acoustic(iglob)*hp1
            enddo
          else
            !AXISYM but not on the axis
            do k = 1,NGLLX
              hp1 = hprime_xx(i,k)
              iglob = ibool(k,j,ispec)
              tempx1l = tempx1l + field_acoustic(iglob)*hp1
            enddo
          endif
        else
          !not AXISYM
          do k = 1,NGLLX
            hp1 = hprime_xx(i,k)
            iglob = ibool(k,j,ispec)
            tempx1l = tempx1l + field_acoustic(iglob)*hp1
          enddo
        endif

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + field_acoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        if (assign_external_model) rhol = rhoext(i,j,ispec)

        ! derivatives of potential
        vector_field_element(1,i,j) = real((tempx1l*xixl + tempx2l*gammaxl) / rhol,kind=CUSTOM_REAL) ! u_x = 1/rho dChi/dx
        vector_field_element(2,i,j) = real((tempx1l*xizl + tempx2l*gammazl) / rhol,kind=CUSTOM_REAL) ! u_z = 1/rho dChi/dz
      enddo
    enddo

  endif ! end of test if acoustic or elastic element

  end subroutine compute_vector_one_element


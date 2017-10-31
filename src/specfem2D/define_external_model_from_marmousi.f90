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

  subroutine define_external_model_from_marmousi(coord,ibool,rho,vp,vs,QKappa_attenuation,Qmu_attenuation, &
                                                 c11,c13,c15,c33,c35,c55,c12,c23,c25,nspec,nglob)

! reads in external model using a marmousi format which defines a compaction gradient


  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IMAIN

  use specfem_par, only: poroelastcoef,density,kmato,myrank

  implicit none

  integer, intent(in) :: nspec,nglob

  double precision, dimension(NDIM,nglob), intent(in) :: coord

  integer, dimension(NGLLX,NGLLZ,nspec), intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: rho,vp,vs
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: QKappa_attenuation,Qmu_attenuation
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: c11,c15,c13,c33,c35,c55,c12,c23,c25

  ! local parameters
  integer :: i,j,ispec,iglob,imat

  double precision :: x,z,rho0,vp0,vs0,comp_grad

! dummy routine here, just to demonstrate how the model can be assigned
! and how such a routine can be written

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  marmousi model: ','using compaction gradient'
    call flush_IMAIN()
  endif

  ! default no attenuation
  QKappa_attenuation(:,:,:) = 9999.d0
  Qmu_attenuation(:,:,:) = 9999.d0

! loop on all the elements of the mesh, and inside each element loop on all the GLL points
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX

        iglob = ibool(i,j,ispec)
        imat = kmato(ispec)

        ! initial model values
        rho0 = density(1,imat)
        vp0 = sqrt(poroelastcoef(3,1,imat)/rho0)
        vs0 = sqrt(poroelastcoef(2,1,imat)/rho0)
        comp_grad = poroelastcoef(4,1,imat)

        x = coord(1,iglob)
        z = coord(2,iglob)

        ! updates model values using depth and compaction gradient information
        rho(i,j,ispec) = rho0
        vp(i,j,ispec) = vp0 + comp_grad*z

  ! assumes Poisson solids
        vs(i,j,ispec) = vp(i,j,ispec) / sqrt(3.d0)

        QKappa_attenuation(i,j,ispec) = 9999. ! this means no attenuation
        Qmu_attenuation(i,j,ispec)    = 9999. ! this means no attenuation

        c11(i,j,ispec) = 0.d0   ! this means no anisotropy
        c13(i,j,ispec) = 0.d0
        c15(i,j,ispec) = 0.d0
        c33(i,j,ispec) = 0.d0
        c35(i,j,ispec) = 0.d0
        c55(i,j,ispec) = 0.d0
        c12(i,j,ispec) = 0.d0
        c23(i,j,ispec) = 0.d0
        c25(i,j,ispec) = 0.d0

      enddo
    enddo
  enddo

  end subroutine define_external_model_from_marmousi



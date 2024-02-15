!========================================================================
!
!                            S P E C F E M 2 D
!                            -----------------
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

  subroutine define_external_model_from_marmousi(coord,ibool,rho,vp,vs, &
                                                 QKappa_attenuationext,Qmu_attenuationext, &
                                                 c11,c12,c13,c15,c23,c25,c33,c35,c55)

! reads in external model using a marmousi format which defines a compaction gradient

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,IMAIN,ATTENUATION_COMP_MAXIMUM

  use specfem_par, only: myrank,nspec,nglob, &
    poroelastcoef,density,kmato, &
    Qkappa_attenuationcoef,Qmu_attenuationcoef

  implicit none

  double precision, dimension(NDIM,nglob), intent(in) :: coord

  integer, dimension(NGLLX,NGLLZ,nspec), intent(in) :: ibool

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(inout) :: rho,vp,vs
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(inout) :: QKappa_attenuationext,Qmu_attenuationext
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(inout) :: c11,c15,c13,c33,c35,c55,c12,c23,c25

  ! local parameters
  integer :: i,j,ispec,iglob,imat

  double precision :: x,z,rho0,vp0,vs0,compaction_grad
  double precision :: zmin_glob,zmax_glob,ztop_solid
  double precision :: depth

  ! water layer depth of Marmousi2
  double precision, parameter :: WATER_LAYER_DEPTH = 450.d0   ! in m

  ! use compaction gradient
  logical, parameter :: USE_COMPACTION_GRADIENT = .true.

  ! instead of provided Vs, assign Poisson solid's Vs scaled from Vp
  logical, parameter :: USE_POISSON_SOLID_VS = .false.

! dummy routine here, just to demonstrate how the model can be assigned
! and how such a routine can be written

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  marmousi model: '
    if (USE_COMPACTION_GRADIENT) &
      write(IMAIN,*) '    using compaction gradient'
    if (USE_POISSON_SOLID_VS) &
      write(IMAIN,*) '    using Poisson solid (vs scaled from vp)'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! initializes model
  ! no anisotropy
  c11(:,:,:) = 0.d0   ! this means no anisotropy
  c13(:,:,:) = 0.d0
  c15(:,:,:) = 0.d0
  c33(:,:,:) = 0.d0
  c35(:,:,:) = 0.d0
  c55(:,:,:) = 0.d0
  c12(:,:,:) = 0.d0
  c23(:,:,:) = 0.d0
  c25(:,:,:) = 0.d0

  ! get model dimension to estimate depth for a given point
  if (USE_COMPACTION_GRADIENT) then
    ! z-coordinate min/max
    zmin_glob = minval(coord(2,:))
    call min_all_all_dp(zmin_glob,zmin_glob)
    zmax_glob = maxval(coord(2,:))
    call max_all_all_dp(zmax_glob,zmax_glob)

    ! taking off the water layer depth
    ztop_solid = zmax_glob - WATER_LAYER_DEPTH

    if (myrank == 0) then
      write(IMAIN,*) '    for material compaction:'
      write(IMAIN,*) '    top         : z_max = ',sngl(zmax_glob)
      write(IMAIN,*) '    bottom      : z_min = ',sngl(zmin_glob)
      write(IMAIN,*)
      write(IMAIN,*) '    water layer : depth = ',sngl(WATER_LAYER_DEPTH)
      write(IMAIN,*) '    top solid   : z_max = ',sngl(ztop_solid)
      write(IMAIN,*)
      call flush_IMAIN()
    endif
  endif

  ! loop on all the elements of the mesh, and inside each element loop on all the GLL points
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        ! gets material properties
        imat = kmato(ispec)

        ! initial model values
        rho0 = density(1,imat)
        vs0 = sqrt(poroelastcoef(2,1,imat)/rho0)
        vp0 = sqrt(poroelastcoef(3,1,imat)/rho0)

        QKappa_attenuationext(i,j,ispec) = QKappa_attenuationcoef(imat)
        Qmu_attenuationext(i,j,ispec) = Qmu_attenuationcoef(imat)

        ! add compaction
        if (USE_COMPACTION_GRADIENT) then
          ! gradient
          compaction_grad = poroelastcoef(4,1,imat)

          ! point position
          iglob = ibool(i,j,ispec)
          x = coord(1,iglob)
          z = coord(2,iglob)

          ! depth (in solid)
          depth = ztop_solid - z

          !debug
          !if (i==1 .and. j==1 .and. compaction_grad > 0.01) &
          !  print *,'debug: marmousi ',ispec,' coord x/z = ',x,z, &
          !        ' depth = ',depth,' compaction: ',compaction_grad,' vp compaction: ',compaction_grad * depth, &
          !        ' vp0/vs0/rho0: ',vp0,vs0,rho0

          ! updates model values using depth and compaction gradient information
          rho(i,j,ispec) = rho0
          vp(i,j,ispec) = vp0 + compaction_grad * depth
          vs(i,j,ispec) = vs0 + compaction_grad * depth

          ! limits value range
          if (vs(i,j,ispec) < 0.d0) vs(i,j,ispec) = 0.d0

          ! vp should stay positive
          if (vp(i,j,ispec) < 0.d0) then
            print *,'Error: compaction leads to negative Vp: ',vp(i,j,ispec)
            print *,'       position: x/z = ',x,'/',z,' depth = ',depth
            call stop_the_code('Error invalid Vp in marmousi model')
          endif
        else
          ! no compaction
          rho(i,j,ispec) = rho0
          vp(i,j,ispec) = vp0
          vs(i,j,ispec) = vs0
        endif

        ! overwrites vs with Poisson vs scaled from vp
        if (USE_POISSON_SOLID_VS) then
            ! assumes Poisson solids
            vs(i,j,ispec) = vp(i,j,ispec) / sqrt(3.d0)
        endif

      enddo
    enddo
  enddo

  end subroutine define_external_model_from_marmousi



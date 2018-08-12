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


  subroutine prepare_absorb_read_elastic()

! reads in previously stored Stacey boundary contributions

  use constants, only: OUTPUT_FILES
  use specfem_par, only: myrank,any_elastic,SIMULATION_TYPE, &
                         nspec_left,nspec_right,nspec_bottom,nspec_top, &
                         b_absorb_elastic_left,b_absorb_elastic_right, &
                         b_absorb_elastic_bottom,b_absorb_elastic_top

  implicit none
  ! local parameters
  integer :: ier
  character(len=150) :: outputname

  ! checks if anything to do in this slice
  if (.not. any_elastic) return

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return

  !--- left absorbing boundary
  if (nspec_left > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_elastic_left',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_elastic_left
    close(35)
  endif

  !--- right absorbing boundary
  if (nspec_right > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_elastic_right',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_elastic_right
    close(35)
  endif

  !--- bottom absorbing boundary
  if (nspec_bottom > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_elastic_bottom',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_elastic_bottom
    close(35)
  endif

  !--- top absorbing boundary
  if (nspec_top > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_elastic_top',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_elastic_top
    close(35)
  endif

  end subroutine prepare_absorb_read_elastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_absorb_read_poroelastic()

  use constants, only: OUTPUT_FILES
  use specfem_par, only: myrank,any_poroelastic,SIMULATION_TYPE, &
                         nspec_left,nspec_right,nspec_bottom,nspec_top, &
                         b_absorb_poro_s_left,b_absorb_poro_w_left, &
                         b_absorb_poro_s_right,b_absorb_poro_w_right, &
                         b_absorb_poro_s_bottom,b_absorb_poro_w_bottom, &
                         b_absorb_poro_s_top,b_absorb_poro_w_top

  implicit none
  ! local parameters
  integer :: ier
  character(len=150) :: outputname,outputname2

  ! checks if anything to do in this slice
  if (.not. any_poroelastic) return

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return

  !--- left absorbing boundary
  if (nspec_left > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_poro_s_left',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_left',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_poro_s_left
    read(36) b_absorb_poro_w_left
    close(35)
    close(36)
  endif

  !--- right absorbing boundary
  if (nspec_right > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_poro_s_right',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_right',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_poro_s_right
    read(36) b_absorb_poro_w_right
    close(35)
    close(36)
  endif

  !--- bottom absorbing boundary
  if (nspec_bottom > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_poro_s_bottom',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_bottom',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_poro_s_bottom
    read(36) b_absorb_poro_w_bottom
    close(35)
    close(36)
  endif

  !--- top absorbing boundary
  if (nspec_top > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_poro_s_top',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_top',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    open(unit=36,file=trim(OUTPUT_FILES)//outputname2,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_poro_s_top
    read(36) b_absorb_poro_w_top
    close(35)
    close(36)
  endif

  end subroutine prepare_absorb_read_poroelastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_absorb_read_acoustic()

  use constants, only: OUTPUT_FILES
  use specfem_par, only: myrank,any_acoustic,SIMULATION_TYPE, &
                         nspec_left,nspec_right,nspec_bottom,nspec_top, &
                         b_absorb_acoustic_left,b_absorb_acoustic_right, &
                         b_absorb_acoustic_bottom,b_absorb_acoustic_top

  implicit none

  ! local parameters
  integer :: ier
  character(len=150) :: outputname

  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  ! checks if anything to do
  if (SIMULATION_TYPE /= 3) return

  !--- left absorbing boundary
  if (nspec_left > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_acoustic_left',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_acoustic_left
    close(35)
  endif

  !--- right absorbing boundary
  if (nspec_right > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_acoustic_right',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_acoustic_right
    close(35)
  endif

  !--- bottom absorbing boundary
  if (nspec_bottom > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_acoustic_bottom',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_acoustic_bottom
    close(35)
  endif

  !--- top absorbing boundary
  if (nspec_top > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_acoustic_top',myrank,'.bin'
    open(unit=35,file=trim(OUTPUT_FILES)//outputname,status='old',form='unformatted',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    ! reads in boundary contributions
    read(35) b_absorb_acoustic_top
    close(35)
  endif

  end subroutine prepare_absorb_read_acoustic

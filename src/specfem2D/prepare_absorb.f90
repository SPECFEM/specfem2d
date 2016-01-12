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
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
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

  subroutine prepare_absorb_open_files_elastic()

  use specfem_par, only: myrank,any_elastic,SIMULATION_TYPE, &
                         nspec_left,nspec_right,nspec_bottom,nspec_top

  implicit none

  ! local parameters
  integer :: ier
  character(len=150) :: outputname

  ! checks if anything to do in this slice
  if (.not. any_elastic) return

  !--- left absorbing boundary
  if (nspec_left > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_elastic_left',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=35,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=35,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif  !  end of left absorbing boundary

  !--- right absorbing boundary
  if (nspec_right > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_elastic_right',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=36,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=36,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif  !  end of right absorbing boundary

  !--- bottom absorbing boundary
  if (nspec_bottom > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_elastic_bottom',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=37,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=37,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif  !  end of bottom absorbing boundary

  !--- top absorbing boundary
  if (nspec_top > 0) then
      write(outputname,'(a,i6.6,a)') 'absorb_elastic_top',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=38,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=38,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif ! end of top absorbing boundary

  end subroutine prepare_absorb_open_files_elastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_absorb_open_files_poroelastic()

  use specfem_par, only: myrank,any_poroelastic,SIMULATION_TYPE, &
                         nspec_left,nspec_right,nspec_bottom,nspec_top

  implicit none

  ! local parameters
  integer :: ier
  character(len=150) :: outputname,outputname2

  ! checks if anything to do in this slice
  if (.not. any_poroelastic) return

  !--- left absorbing boundary
  if (nspec_left > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_poro_s_left',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_left',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=45,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
      open(unit=25,file='OUTPUT_FILES/'//outputname2,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=45,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
      open(unit=25,file='OUTPUT_FILES/'//outputname2,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif  !  end of left absorbing boundary

  !--- right absorbing boundary
  if (nspec_right > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_poro_s_right',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_right',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=46,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
      open(unit=26,file='OUTPUT_FILES/'//outputname2,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=46,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
      open(unit=26,file='OUTPUT_FILES/'//outputname2,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif  !  end of right absorbing boundary

  !--- bottom absorbing boundary
  if (nspec_bottom > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_poro_s_bottom',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_bottom',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=47,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
      open(unit=29,file='OUTPUT_FILES/'//outputname2,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=47,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
      open(unit=29,file='OUTPUT_FILES/'//outputname2,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif  !  end of bottom absorbing boundary

  !--- top absorbing boundary
  if (nspec_top > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_poro_s_top',myrank,'.bin'
    write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_top',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=48,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
      open(unit=28,file='OUTPUT_FILES/'//outputname2,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=48,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
      open(unit=28,file='OUTPUT_FILES/'//outputname2,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif ! end of top absorbing boundary

  end subroutine prepare_absorb_open_files_poroelastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_absorb_open_files_acoustic()

  use specfem_par, only: myrank,any_acoustic,SIMULATION_TYPE, &
                         nspec_left,nspec_right,nspec_bottom,nspec_top

  implicit none

  ! local parameters
  integer :: ier
  character(len=150) :: outputname

  ! checks if anything to do in this slice
  if (.not. any_acoustic) return

  !--- left absorbing boundary
  if (nspec_left > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_acoustic_left',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=65,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=65,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif  !  end of left absorbing boundary

  !--- right absorbing boundary
  if (nspec_right > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_acoustic_right',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=66,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=66,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif  !  end of right absorbing boundary

  !--- bottom absorbing boundary
  if (nspec_bottom > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_acoustic_bottom',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=67,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=67,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif  !  end of bottom absorbing boundary

  !--- top absorbing boundary
  if (nspec_top > 0) then
    write(outputname,'(a,i6.6,a)') 'absorb_acoustic_top',myrank,'.bin'
    if (SIMULATION_TYPE == 3) then
      open(unit=68,file='OUTPUT_FILES/'//outputname,status='old',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    else
      open(unit=68,file='OUTPUT_FILES/'//outputname,status='unknown',&
            form='unformatted',iostat=ier)
      if (ier /= 0) call exit_MPI(myrank,'Error opening absorbing boundary file')
    endif

  endif ! end of top absorbing boundary

  end subroutine prepare_absorb_open_files_acoustic


!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_absorb_read_elastic()

  use specfem_par, only: nspec_left,nspec_right,nspec_bottom,nspec_top, &
                         b_absorb_elastic_left,b_absorb_elastic_right, &
                         b_absorb_elastic_bottom,b_absorb_elastic_top

  implicit none

  !--- left absorbing boundary
  if (nspec_left > 0) read(35) b_absorb_elastic_left

  !--- right absorbing boundary
  if (nspec_right > 0) read(36) b_absorb_elastic_right

  !--- bottom absorbing boundary
  if (nspec_bottom > 0) read(37) b_absorb_elastic_bottom

  !--- top absorbing boundary
  if (nspec_top > 0) read(38) b_absorb_elastic_top

  end subroutine prepare_absorb_read_elastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_absorb_read_poroelastic()

  use specfem_par, only: nspec_left,nspec_right,nspec_bottom,nspec_top, &
                         b_absorb_poro_s_left,b_absorb_poro_w_left, &
                         b_absorb_poro_s_right,b_absorb_poro_w_right, &
                         b_absorb_poro_s_bottom,b_absorb_poro_w_bottom, &
                         b_absorb_poro_s_top,b_absorb_poro_w_top

  implicit none

  !--- left absorbing boundary
  if (nspec_left > 0) then
    read(45) b_absorb_poro_s_left
    read(25) b_absorb_poro_w_left
  endif

  !--- right absorbing boundary
  if (nspec_right > 0) then
    read(46) b_absorb_poro_s_right
    read(26) b_absorb_poro_w_right
  endif

  !--- bottom absorbing boundary
  if (nspec_bottom > 0)  then
    read(47) b_absorb_poro_s_bottom
    read(29) b_absorb_poro_w_bottom
  endif

  !--- top absorbing boundary
  if (nspec_top > 0) then
    read(48) b_absorb_poro_s_top
    read(28) b_absorb_poro_w_top
  endif

  end subroutine prepare_absorb_read_poroelastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_absorb_read_acoustic()

  use specfem_par, only: nspec_left,nspec_right,nspec_bottom,nspec_top, &
                         b_absorb_acoustic_left,b_absorb_acoustic_right, &
                         b_absorb_acoustic_bottom,b_absorb_acoustic_top

  implicit none

  !--- left absorbing boundary
  if (nspec_left > 0)  read(65) b_absorb_acoustic_left

  !--- right absorbing boundary
  if (nspec_right > 0) read(66) b_absorb_acoustic_right

  !--- bottom absorbing boundary
  if (nspec_bottom > 0) read(67) b_absorb_acoustic_bottom

  !--- top absorbing boundary
  if (nspec_top > 0) read(68) b_absorb_acoustic_top

  end subroutine prepare_absorb_read_acoustic

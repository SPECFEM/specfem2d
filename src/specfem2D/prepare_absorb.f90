
!========================================================================
!
!                   S P E C F E M 2 D  Version 6 . 2
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
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


  subroutine prepare_absorb_files(myrank,any_elastic,any_poroelastic,any_acoustic, &
                      nspec_left,nspec_right,nspec_bottom,nspec_top,SIMULATION_TYPE)

  implicit none
  include "constants.h"

  integer :: myrank,SIMULATION_TYPE
  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top
  logical :: any_elastic,any_poroelastic,any_acoustic

  ! local parameters
  character(len=150) :: outputname,outputname2


  if(any_elastic) then

    !--- left absorbing boundary
    if( nspec_left >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_elastic_left',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=35,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
      else
        open(unit=35,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
      endif

    endif  !  end of left absorbing boundary

    !--- right absorbing boundary
    if( nspec_right >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_elastic_right',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=36,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
      else
        open(unit=36,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
      endif

    endif  !  end of right absorbing boundary

    !--- bottom absorbing boundary
    if( nspec_bottom >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_elastic_bottom',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=37,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
      else
        open(unit=37,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
      endif

    endif  !  end of bottom absorbing boundary

    !--- top absorbing boundary
    if( nspec_top >0 ) then
        write(outputname,'(a,i6.6,a)') 'absorb_elastic_top',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=38,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
      else
        open(unit=38,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
      endif

    endif ! end of top absorbing boundary

  endif ! any_elastic

  if(any_poroelastic) then

    !--- left absorbing boundary
    if( nspec_left >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_poro_s_left',myrank,'.bin'
      write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_left',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=45,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
        open(unit=25,file='OUTPUT_FILES/'//outputname2,status='old',&
              form='unformatted')
      else
        open(unit=45,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
        open(unit=25,file='OUTPUT_FILES/'//outputname2,status='unknown',&
              form='unformatted')
      endif

    endif  !  end of left absorbing boundary

    !--- right absorbing boundary
    if( nspec_right >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_poro_s_right',myrank,'.bin'
      write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_right',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=46,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
        open(unit=26,file='OUTPUT_FILES/'//outputname2,status='old',&
              form='unformatted')
      else
        open(unit=46,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
        open(unit=26,file='OUTPUT_FILES/'//outputname2,status='unknown',&
              form='unformatted')
      endif

    endif  !  end of right absorbing boundary

    !--- bottom absorbing boundary
    if( nspec_bottom >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_poro_s_bottom',myrank,'.bin'
      write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_bottom',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=47,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
        open(unit=29,file='OUTPUT_FILES/'//outputname2,status='old',&
              form='unformatted')
      else
        open(unit=47,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
        open(unit=29,file='OUTPUT_FILES/'//outputname2,status='unknown',&
              form='unformatted')
      endif

    endif  !  end of bottom absorbing boundary

    !--- top absorbing boundary
    if( nspec_top >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_poro_s_top',myrank,'.bin'
      write(outputname2,'(a,i6.6,a)') 'absorb_poro_w_top',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=48,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
        open(unit=28,file='OUTPUT_FILES/'//outputname2,status='old',&
              form='unformatted')
      else
        open(unit=48,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
        open(unit=28,file='OUTPUT_FILES/'//outputname2,status='unknown',&
              form='unformatted')
      endif

    endif ! end of top absorbing boundary

  endif !any_poroelastic

  if(any_acoustic) then

    !--- left absorbing boundary
    if( nspec_left >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_acoustic_left',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=65,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
      else
        open(unit=65,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
      endif

    endif  !  end of left absorbing boundary

    !--- right absorbing boundary
    if( nspec_right >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_acoustic_right',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=66,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
      else
        open(unit=66,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
      endif

    endif  !  end of right absorbing boundary

    !--- bottom absorbing boundary
    if( nspec_bottom >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_acoustic_bottom',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=67,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
      else
        open(unit=67,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
      endif

    endif  !  end of bottom absorbing boundary

    !--- top absorbing boundary
    if( nspec_top >0 ) then
      write(outputname,'(a,i6.6,a)') 'absorb_acoustic_top',myrank,'.bin'
      if(SIMULATION_TYPE == 2) then
        open(unit=68,file='OUTPUT_FILES/'//outputname,status='old',&
              form='unformatted')
      else
        open(unit=68,file='OUTPUT_FILES/'//outputname,status='unknown',&
              form='unformatted')
      endif

    endif ! end of top absorbing boundary

  endif !any_acoustic


  end subroutine prepare_absorb_files


!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_absorb_elastic(NSTEP,p_sv, &
                      nspec_left,nspec_right,nspec_bottom,nspec_top, &
                      b_absorb_elastic_left,b_absorb_elastic_right, &
                      b_absorb_elastic_bottom,b_absorb_elastic_top)

  implicit none
  include "constants.h"

  logical :: p_sv
  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top
  integer :: NSTEP
  real(kind=CUSTOM_REAL) :: b_absorb_elastic_left(3,NGLLZ,nspec_left,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_elastic_right(3,NGLLZ,nspec_right,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_elastic_bottom(3,NGLLX,nspec_bottom,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_elastic_top(3,NGLLX,nspec_top,NSTEP)

  ! local parameters
  integer :: ispec,i,it

  do it =1, NSTEP

    !--- left absorbing boundary
    if(nspec_left >0) then
      do ispec = 1,nspec_left

        if(p_sv)then!P-SV waves
          do i=1,NGLLZ
            read(35) b_absorb_elastic_left(1,i,ispec,it)
          enddo
          do i=1,NGLLZ
            read(35) b_absorb_elastic_left(3,i,ispec,it)
          enddo
          b_absorb_elastic_left(2,:,ispec,it) = ZERO
        else!SH (membrane) waves
          do i=1,NGLLZ
            read(35) b_absorb_elastic_left(2,i,ispec,it)
          enddo
          b_absorb_elastic_left(1,:,ispec,it) = ZERO
          b_absorb_elastic_left(3,:,ispec,it) = ZERO
        endif

      enddo
    endif

    !--- right absorbing boundary
    if(nspec_right >0) then
      do ispec = 1,nspec_right

        if(p_sv)then!P-SV waves
          do i=1,NGLLZ
            read(36) b_absorb_elastic_right(1,i,ispec,it)
          enddo
          do i=1,NGLLZ
            read(36) b_absorb_elastic_right(3,i,ispec,it)
          enddo
          b_absorb_elastic_right(2,:,ispec,it) = ZERO
        else!SH (membrane) waves
          do i=1,NGLLZ
            read(36) b_absorb_elastic_right(2,i,ispec,it)
          enddo
          b_absorb_elastic_right(1,:,ispec,it) = ZERO
          b_absorb_elastic_right(3,:,ispec,it) = ZERO
        endif

      enddo
    endif

    !--- bottom absorbing boundary
    if(nspec_bottom >0) then
      do ispec = 1,nspec_bottom

        if(p_sv)then!P-SV waves
          do i=1,NGLLX
            read(37) b_absorb_elastic_bottom(1,i,ispec,it)
          enddo
          do i=1,NGLLX
            read(37) b_absorb_elastic_bottom(3,i,ispec,it)
          enddo
          b_absorb_elastic_bottom(2,:,ispec,it) = ZERO
        else!SH (membrane) waves
          do i=1,NGLLZ
            read(37) b_absorb_elastic_bottom(2,i,ispec,it)
          enddo
          b_absorb_elastic_bottom(1,:,ispec,it) = ZERO
          b_absorb_elastic_bottom(3,:,ispec,it) = ZERO
        endif

      enddo
    endif

    !--- top absorbing boundary
    if(nspec_top >0) then
      do ispec = 1,nspec_top

        if(p_sv)then!P-SV waves
          do i=1,NGLLX
            read(38) b_absorb_elastic_top(1,i,ispec,it)
          enddo
          do i=1,NGLLX
            read(38) b_absorb_elastic_top(3,i,ispec,it)
          enddo
          b_absorb_elastic_top(2,:,ispec,it) = ZERO
        else!SH (membrane) waves
          do i=1,NGLLZ
            read(38) b_absorb_elastic_top(2,i,ispec,it)
          enddo
          b_absorb_elastic_top(1,:,ispec,it) = ZERO
          b_absorb_elastic_top(3,:,ispec,it) = ZERO
        endif

      enddo
    endif

  enddo

  end subroutine prepare_absorb_elastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_absorb_poroelastic(NSTEP, &
                      nspec_left,nspec_right,nspec_bottom,nspec_top, &
                      b_absorb_poro_s_left,b_absorb_poro_w_left, &
                      b_absorb_poro_s_right,b_absorb_poro_w_right, &
                      b_absorb_poro_s_bottom,b_absorb_poro_w_bottom, &
                      b_absorb_poro_s_top,b_absorb_poro_w_top)

  implicit none
  include "constants.h"

  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top

  integer :: NSTEP
  real(kind=CUSTOM_REAL) :: b_absorb_poro_s_left(NDIM,NGLLZ,nspec_left,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_poro_s_right(NDIM,NGLLZ,nspec_right,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_poro_s_bottom(NDIM,NGLLX,nspec_bottom,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_poro_s_top(NDIM,NGLLX,nspec_top,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_poro_w_left(NDIM,NGLLZ,nspec_left,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_poro_w_right(NDIM,NGLLZ,nspec_right,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_poro_w_bottom(NDIM,NGLLX,nspec_bottom,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_poro_w_top(NDIM,NGLLX,nspec_top,NSTEP)

  ! local parameters
  integer :: ispec,i,it,id

  do it =1, NSTEP

    !--- left absorbing boundary
    if(nspec_left >0) then
      do ispec = 1,nspec_left
       do id =1,2
         do i=1,NGLLZ
          read(45) b_absorb_poro_s_left(id,i,ispec,it)
          read(25) b_absorb_poro_w_left(id,i,ispec,it)
         enddo
       enddo
      enddo
    endif

    !--- right absorbing boundary
    if(nspec_right >0) then
      do ispec = 1,nspec_right
       do id =1,2
         do i=1,NGLLZ
          read(46) b_absorb_poro_s_right(id,i,ispec,it)
          read(26) b_absorb_poro_w_right(id,i,ispec,it)
         enddo
       enddo
      enddo
    endif

    !--- bottom absorbing boundary
    if(nspec_bottom >0) then
      do ispec = 1,nspec_bottom
       do id =1,2
         do i=1,NGLLX
          read(47) b_absorb_poro_s_bottom(id,i,ispec,it)
          read(29) b_absorb_poro_w_bottom(id,i,ispec,it)
         enddo
       enddo
      enddo
    endif

    !--- top absorbing boundary
    if(nspec_top >0) then
      do ispec = 1,nspec_top
       do id =1,2
         do i=1,NGLLX
          read(48) b_absorb_poro_s_top(id,i,ispec,it)
          read(28) b_absorb_poro_w_top(id,i,ispec,it)
         enddo
       enddo
      enddo
    endif

  enddo

  end subroutine prepare_absorb_poroelastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_absorb_acoustic(NSTEP, &
                      nspec_left,nspec_right,nspec_bottom,nspec_top, &
                      b_absorb_acoustic_left,b_absorb_acoustic_right, &
                      b_absorb_acoustic_bottom,b_absorb_acoustic_top)

  implicit none
  include "constants.h"

  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top

  integer :: NSTEP
  real(kind=CUSTOM_REAL) :: b_absorb_acoustic_left(NGLLZ,nspec_left,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_acoustic_right(NGLLZ,nspec_right,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_acoustic_bottom(NGLLX,nspec_bottom,NSTEP)
  real(kind=CUSTOM_REAL) :: b_absorb_acoustic_top(NGLLX,nspec_top,NSTEP)


  ! local parameters
  integer :: ispec,i,it

  do it =1, NSTEP

    !--- left absorbing boundary
    if(nspec_left >0) then
      do ispec = 1,nspec_left
         do i=1,NGLLZ
          read(65) b_absorb_acoustic_left(i,ispec,it)
         enddo
      enddo
    endif

    !--- right absorbing boundary
    if(nspec_right >0) then
      do ispec = 1,nspec_right
         do i=1,NGLLZ
          read(66) b_absorb_acoustic_right(i,ispec,it)
         enddo
      enddo
    endif

    !--- bottom absorbing boundary
    if(nspec_bottom >0) then
      do ispec = 1,nspec_bottom
         do i=1,NGLLX
          read(67) b_absorb_acoustic_bottom(i,ispec,it)
         enddo
      enddo
    endif

    !--- top absorbing boundary
    if(nspec_top >0) then
      do ispec = 1,nspec_top
         do i=1,NGLLX
          read(68) b_absorb_acoustic_top(i,ispec,it)
         enddo
      enddo
    endif

  enddo

  end subroutine prepare_absorb_acoustic


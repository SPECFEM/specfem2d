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


  subroutine prepare_noise()

! for noise simulations

  use constants, only: NGLLX,NGLLZ,NDIM,IMAIN,NOISE_MOVIE_OUTPUT,TWO_THIRDS,OUTPUT_FILES, &
                       MAX_STRING_LEN

  use specfem_par, only: myrank,NSTEP,nglob,nspec,ibool,coord, &
                         rhostore,rho_vpstore,rho_vsstore, &
                         NOISE_TOMOGRAPHY

  use specfem_par_noise

  implicit none

  integer :: i,j,iglob,ispec,ier
  double precision :: rhol,vs,vp
  character(len=MAX_STRING_LEN) :: fname

  ! checks if anything to do
  if (NOISE_TOMOGRAPHY <= 0) return

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Preparing noise simulation'
    call flush_IMAIN()
  endif

  ! allocates arrays for noise tomography
  allocate(noise_sourcearray(NDIM,NGLLX,NGLLZ,NSTEP), &
           mask_noise(nglob), &
           noise_surface_movie_y_or_z(nglob),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating noise arrays')
  noise_sourcearray(:,:,:,:) = 0._CUSTOM_REAL
  mask_noise(:) = 0._CUSTOM_REAL
  noise_surface_movie_y_or_z(:) = 0._CUSTOM_REAL

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  reading noise parameters'
    call flush_IMAIN()
  endif

  !read in parameters for noise tomography
  call read_parameters_noise()

  !read noise distribution
  call read_noise_distribution(mask_noise)

  ! file outputs for visualization/debugging
  if (NOISE_TOMOGRAPHY == 1) then
    ! write out coordinates of mesh
    write(fname,'(a,i6.6)') 'mesh_spec',myrank
    open(unit=504,file=trim(OUTPUT_FILES)//fname,status='unknown',action='write')
    do ispec = 1, nspec
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          write(504,'(1pe11.3,1pe11.3,2i3,i7)') coord(1,iglob), coord(2,iglob), i, j, ispec
       enddo
      enddo
    enddo
    close(504)

    write(fname,'(a,i6.6)') 'mesh_glob',myrank
    open(unit=504,file=trim(OUTPUT_FILES)//fname,status='unknown',action='write')
    do iglob = 1, nglob
      write(504,'(1pe11.3,1pe11.3,i7)') coord(1,iglob), coord(2,iglob), iglob
    enddo
    close(504)

    ! write out spatial distribution of noise sources
    write(fname,'(a,i6.6)') 'mask_noise',myrank
    open(unit=504,file=trim(OUTPUT_FILES)//fname,status='unknown',action='write')
    do iglob = 1, nglob
      write(504,'(1pe11.3,1pe11.3,1pe11.3)') coord(1,iglob), coord(2,iglob), mask_noise(iglob)
    enddo
    close(504)

    ! write out velocity model
    write(fname,'(a,i6.6)') 'model_rho_vp_vs',myrank
    open(unit=504,file=trim(OUTPUT_FILES)//fname,status='unknown',action='write')
    do ispec = 1, nspec
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          rhol = dble(rhostore(i,j,ispec))
          vs = dble(rho_vsstore(i,j,ispec)/rhol)
          vp = dble(rho_vpstore(i,j,ispec)/rhol)

          write(504,'(1pe11.3,1pe11.3,1pe11.3,1pe11.3,1pe11.3)') coord(1,iglob), coord(2,iglob), rhol, vp, vs
        enddo
      enddo
    enddo
    close(504)
  endif

  ! noise movie snapshots
  if (NOISE_TOMOGRAPHY == 3 .and. NOISE_MOVIE_OUTPUT) then
    ! noise movie
    ! prepare array that will hold wavefield snapshots
    noise_output_ncol = 5
    allocate(noise_output_array(noise_output_ncol,nglob), &
             noise_output_rhokl(nglob))
    noise_output_array(:,:) = 0._CUSTOM_REAL
    noise_output_rhokl(:) = 0._CUSTOM_REAL
  endif

  ! synchronizes all processes
  call synchronize_all()

  end subroutine prepare_noise


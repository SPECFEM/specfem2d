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

!NOISE TOMOGRAPHY TO DO LIST

!1. Use separate STATIONS_ADJOINT file
!2. Add exploration test case under EXAMPLES
!3. Update manual

! =============================================================================================================
! specify spatial distribution of microseismic noise sources
! USERS need to modify this subroutine to suit their own needs
  subroutine create_mask_noise()

  use specfem_par, only: nglob,coord,mask_noise
  implicit none
  include "constants.h"

  !local
  integer :: iglob
  real(kind=CUSTOM_REAL) :: xx,zz

  !specify distribution of noise sources as a function of xx, zz
  do iglob = 1,nglob
    xx = coord(1,iglob)
    zz = coord(2,iglob)

    !below, the noise is assumed to be uniform; users are free to
    !to change this expression to one involving xx, zz
    mask_noise(iglob) = 1.0

  enddo

  end subroutine create_mask_noise


! =============================================================================================================
! read noise parameters and check for consistency
  subroutine read_parameters_noise()

  use specfem_par, only: NOISE_TOMOGRAPHY,SIMULATION_TYPE,SAVE_FORWARD, &
                                     any_acoustic,any_poroelastic,p_sv,irec_master, &
                                     Mxx,Mxz,Mzz,factor,NSOURCES, &
                                     xi_receiver,gamma_receiver,ispec_selected_rec,nrec, &
                                     xi_noise,gamma_noise,ispec_noise,angle_noise

  implicit none
  include "constants.h"

  !local
  integer :: i,ios

  !define master receiver
  open(unit=509,file='DATA/NOISE_TOMOGRAPHY/irec_master',status='old',action='read',iostat=ios)
  if (ios == 0) then
    read(509,*) irec_master
  else
    irec_master=1
    write(*,*) ''
    write(*,*) 'Error opening DATA/NOISE_TOMOGRAPHY/irec_master.'
    write(*,*) 'Using irec_master=1. Continuing.'
    write(*,*) ''
  endif
  close(509)

  if ( (NOISE_TOMOGRAPHY == 1) .and. (irec_master > nrec .or. irec_master < 1) ) &
    call exit_mpi('irec_master out of range of given number of receivers. Exiting.')

  xi_noise    = xi_receiver(irec_master)
  gamma_noise = gamma_receiver(irec_master)
  ispec_noise = ispec_selected_rec(irec_master)
  angle_noise = 0.d0

! check simulation parameters

  if ((NOISE_TOMOGRAPHY /= 0) .and. (p_sv)) write(*,*) 'Warning: For P-SV case, noise tomography subroutines not yet fully tested'

  if (NOISE_TOMOGRAPHY == 1) then
     if (SIMULATION_TYPE /= 1) call exit_mpi('NOISE_TOMOGRAPHY=1 requires SIMULATION_TYPE=1    -> check DATA/Par_file')

  else if (NOISE_TOMOGRAPHY == 2) then
     if (SIMULATION_TYPE /= 1) call exit_mpi('NOISE_TOMOGRAPHY=2 requires SIMULATION_TYPE=1    -> check DATA/Par_file')
     if (.not. SAVE_FORWARD) call exit_mpi('NOISE_TOMOGRAPHY=2 requires SAVE_FORWARD=.true.  -> check DATA/Par_file')

  else if (NOISE_TOMOGRAPHY == 3) then
     if (SIMULATION_TYPE /= 3) call exit_mpi('NOISE_TOMOGRAPHY=3 requires SIMULATION_TYPE=3    -> check DATA/Par_file')
     if (SAVE_FORWARD)       call exit_mpi('NOISE_TOMOGRAPHY=3 requires SAVE_FORWARD=.false. -> check DATA/Par_file')
  endif

!  check model parameters
   if (any_acoustic)    call exit_mpi('Acoustic models not yet implemented for noise simulations. Exiting.')
   if (any_poroelastic) call exit_mpi('Poroelastic models not yet implemented for noise simulations. Exiting.')

!  moment tensor elements must be zero!
   do i=1,NSOURCES
     if (Mxx(i) /= 0.d0 .or. Mxz(i) /= 0.d0 .or. Mzz(i) /= 0.d0 .or. factor(i) /= 0.d0) then
       call exit_mpi('For noise simulations, all moment tensor elements must be zero. Exiting.')
     endif
   enddo

  end subroutine read_parameters_noise


! =============================================================================================================
! read in time series based on noise spectrum and construct noise "source" array
  subroutine compute_source_array_noise()


  use specfem_par, only: AXISYM,is_on_the_axis,xiglj,p_sv,NSTEP,deltat,ibool,ispec_noise, &
                       xi_noise,gamma_noise,xigll,zigll, &
                       time_function_noise,source_array_noise

  implicit none
  include "constants.h"

  !local
  integer :: it,i,j,iglob
  real(kind=CUSTOM_REAL) :: t
  double precision, dimension(NGLLX) :: hxi, hpxi
  double precision, dimension(NGLLZ) :: hgamma, hpgamma
  real(kind=CUSTOM_REAL) :: factor_noise, aval, t0

! ---------------------------------------------------------------------------------
! A NOTE ABOUT TIME FUNCTIONS FOR NOISE SIMULATIONS
!
! In noise forward modeling and inversion, "time function" is used to refer
! to the source time function that drives the first noise simulation.
!
! Typically, the time function must be computed based on the frequency
! characterstics of the desired noise sources. For example, if you wanted to
! generate synthetic correlograms representative of the noise field in a
! particular geographic region, you would have to use a time function created
! by inverse Fourier transforming a model of the noise spectrum for that
! region.
!
! IN CASES SUCH AS THIS--WHERE THE USER REQUIRES A REALISTIC MODEL OF THE
! SEISMIC NOISE FIELD--THE VARIABE "time_function_type" SHOULD BE SET TO 0
! AND A TIME FUNCTION ENCODING THE DESIRED NOISE SPECTRUM SHOULD BE
! ***SUPPLIED BY THE USER*** IN THE FILE "DATA/NOISE_TOMOGRAPHY/S_squared"
!
! The alternative--setting "time_function_type" to a value other than
! 0--results in an idealized time function that does represent a physically
! realizable noise spectrum but which nevertheless may be useful for
! performing tests or illustrating certain theoretical concepts.
!
! ----------------------------------------------------------------------------------

  !the following values are chosen to reproduce the time function from Fig 2a of
  !"Tromp et al., 2010, Noise Cross-Correlation Sensitivity Kernels"

  integer, parameter :: time_function_type = 4

  time_function_noise(:) = 0._CUSTOM_REAL
  t0   = ((NSTEP-1)/2.)*deltat
  aval = 0.6_CUSTOM_REAL
  factor_noise = 1.e3_CUSTOM_REAL


  if ( time_function_type == 0) then
    !read in time function from file S_squared
    open(unit=55,file='DATA/NOISE_TOMOGRAPHY/S_squared',status='old')
    do it = 1,NSTEP
      READ(55,*) time_function_noise(it)
    enddo
    close(55)


  else if( time_function_type == 1) then
    !Ricker (second derivative of a Gaussian) time function
    do it = 1,NSTEP
      t = it*deltat
      time_function_noise(it) = - factor_noise * 2.*aval * (1. - 2.*aval*(t-t0)**2.) * &
                                 exp(-aval*(t-t0)**2.)
    enddo


  else if( time_function_type == 2) then
    !first derivative of a Gaussian time function
    do it = 1,NSTEP
      t = it*deltat
      time_function_noise(it) = - factor_noise * (2.*aval*(t-t0)) * exp(-aval*(t-t0)**2.)
    enddo


  else if( time_function_type == 3) then
    !Gaussian time function
    do it = 1,NSTEP
      t = it*deltat
      time_function_noise(it) = factor_noise * exp(-aval*(t-t0)**2.)
    enddo


  else if( time_function_type == 4 ) then
    !reproduce time function from Figure 2a of Tromp et al. 2010
    do it = 1,NSTEP
      t = it*deltat
      time_function_noise(it) = factor_noise * &
       4.*aval**2. * (3. - 12.*aval*(t-t0)**2. + 4.*aval**2.*(t-t0)**4.) * &
       exp(-aval*(t-t0)**2.)
    enddo


  else
    call exit_MPI('Bad time_function_type in compute_source_array_noise.')

  endif

  !interpolate over GLL points
  source_array_noise(:,:,:,:) = 0._CUSTOM_REAL

  if (AXISYM) then
    if(is_on_the_axis(ispec_noise)) then
      call lagrange_any(xi_noise,NGLJ,xiglj,hxi,hpxi)
    else
      call lagrange_any(xi_noise,NGLLX,xigll,hxi,hpxi)
    endif
  else
    call lagrange_any(xi_noise,NGLLX,xigll,hxi,hpxi)
  endif

  call lagrange_any(gamma_noise,NGLLZ,zigll,hgamma,hpgamma)

  if(p_sv) then ! P-SV simulation
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec_noise)
        source_array_noise(1,i,j,:) = time_function_noise(:) * hxi(i) * hgamma(j)
        source_array_noise(3,i,j,:) = time_function_noise(:) * hxi(i) * hgamma(j)
      enddo
    enddo
  else ! SH (membrane) simulation
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec_noise)
        source_array_noise(2,i,j,:) = time_function_noise(:) * hxi(i) * hgamma(j)
      enddo
    enddo
  endif

  end subroutine compute_source_array_noise


! =============================================================================================================
! inject the "source" that drives the "generating wavefield"
  subroutine add_point_source_noise()

  use specfem_par, only: p_sv,it,ibool,ispec_noise, &
                         accel_elastic,angle_noise,source_array_noise
  implicit none
  include "constants.h"

  !local
  integer :: i,j,iglob

  if(p_sv) then ! P-SV calculation
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec_noise)
        accel_elastic(1,iglob) = accel_elastic(1,iglob) &
          + sin(angle_noise)*source_array_noise(1,i,j,it)
        accel_elastic(3,iglob) = accel_elastic(3,iglob) &
          - cos(angle_noise)*source_array_noise(3,i,j,it)
      enddo
    enddo
  else ! SH (membrane) calculation
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec_noise)
        accel_elastic(2,iglob) = accel_elastic(2,iglob) - source_array_noise(2,i,j,it)
      enddo
    enddo
  endif

  end subroutine add_point_source_noise


! =============================================================================================================
! read in and inject the "source" that drives the "enemble forward wavefield"
! (recall that the ensemble forward wavefield has a spatially distributed source)
  subroutine add_surface_movie_noise(accel_elastic)

  use specfem_par, only: p_sv,NOISE_TOMOGRAPHY,it,NSTEP,nspec,nglob,ibool, &
                         surface_movie_x_noise,surface_movie_y_noise, &
                         surface_movie_z_noise,mask_noise,jacobian,wxgll,wzgll

  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(3,nglob) :: accel_elastic

  !local
  integer :: ios, i, j, ispec, iglob


  if (it==1) then
    open(unit=500,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/eta',access='direct', &
         recl=nglob*CUSTOM_REAL,action='read',iostat=ios)
    if( ios /= 0) call exit_mpi('Error retrieving generating wavefield.')
  endif

  if(p_sv) then
    call exit_mpi('P-SV case not yet implemented.')
  else
    if (NOISE_TOMOGRAPHY==2) read(unit=500,rec=NSTEP-it+1) surface_movie_y_noise
    if (NOISE_TOMOGRAPHY==3) read(unit=500,rec=it) surface_movie_y_noise
  endif

  if (it==NSTEP) then
    !close(500)
  endif

  do ispec = 1, nspec
    do j = 1, NGLLZ
      do i = 1, NGLLX
        if (p_sv) then ! P-SV calculation
          iglob = ibool(i,j,ispec)
          accel_elastic(1,iglob) = accel_elastic(1,iglob) + surface_movie_x_noise(iglob) * &
                                   mask_noise(iglob) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
          accel_elastic(3,iglob) = accel_elastic(3,iglob) + surface_movie_z_noise(iglob) * &
                                   mask_noise(iglob) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)

        else ! SH (membrane) calculation
          iglob = ibool(i,j,ispec)
          accel_elastic(2,iglob) = accel_elastic(2,iglob) + surface_movie_y_noise(iglob) * &
                                   mask_noise(iglob) * wxgll(i)*wzgll(j)*jacobian(i,j,ispec)
        endif
      enddo
    enddo
  enddo

  end subroutine add_surface_movie_noise

! =============================================================================================================
! save a snapshot of the "generating wavefield" eta that will be used to drive
! the "ensemble forward wavefield"
  subroutine save_surface_movie_noise()

  use specfem_par, only: NOISE_TOMOGRAPHY,p_sv,it,NSTEP,nglob,displ_elastic

  implicit none
  include "constants.h"

  !local
  integer :: ios

  if (it==1) then

    if (NOISE_TOMOGRAPHY == 1) then

      open(unit=500,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/eta',access='direct', &
           recl=nglob*CUSTOM_REAL,action='write',iostat=ios)
      if( ios /= 0) call exit_mpi('Error saving generating wavefield.')

    else if (NOISE_TOMOGRAPHY == 2) then

      open(unit=500,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/phi',access='direct', &
           recl=nglob*CUSTOM_REAL,action='write',iostat=ios)
      if( ios /= 0) call exit_mpi('Error saving ensemble forward wavefield.')

    else
      call exit_mpi('Bad value of NOISE_TOMOGRAPHY in save_surface_movie_noise.')

    endif

  endif ! (it==1)

  if(p_sv) then
    call exit_mpi('P-SV case not yet implemented.')
  else
    write(unit=500,rec=it) displ_elastic(2,:)
  endif

  if (it==NSTEP) then
    !close(500)
  endif

  end subroutine save_surface_movie_noise

! =============================================================================================================
  subroutine snapshots_noise(ncol,nglob,filename,array_all)

  implicit none
  include "constants.h"

  !input paramters
  integer :: ncol,nglob
  character(len=512) filename

  real(kind=CUSTOM_REAL), dimension(ncol,nglob) :: array_all

  !local parameters
  integer :: i,iglob

  open(unit=504,file=filename,status='unknown',action='write')

    do iglob = 1,nglob

          do i = 1,ncol-1
              write(unit=504,fmt='(1pe12.3e3)',advance='no') array_all(i,iglob)
          enddo
              write(unit=504,fmt='(1pe13.3e3)') array_all(ncol,iglob)

    enddo

  close(504)


  end subroutine snapshots_noise


! =============================================================================================================
! auxillary routine
  subroutine spec2glob(nspec,nglob,ibool,array_spec,array_glob)

  implicit none
  include "constants.h"

  !input
  integer :: nspec, nglob
  integer :: ibool(NGLLX,NGLLZ,nspec)

  real(kind=CUSTOM_REAL), dimension(nglob) :: array_glob
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: array_spec

  !local
  integer :: i,j,iglob,ispec

  do ispec = 1, nspec
    do j = 1, NGLLZ
      do i = 1, NGLLX
        iglob = ibool(i,j,ispec)
        array_glob(iglob) = array_spec(i,j,ispec)
     enddo
    enddo
  enddo

  end subroutine spec2glob


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

!NOISE TOMOGRAPHY TO DO LIST

!1. Replace ad hoc scaling in add_surface_movie_noise by expression involving Jacobian
!2. Use separate STATIONS_ADJOINT file
!3. Add exploration test case under EXAMPLES
!4. Update manual

! =============================================================================================================
! specify spatial distribution of microseismic noise sources
! USERS need to modify this subroutine to suit their own needs
  subroutine create_mask_noise(nglob,coord,mask_noise)

  implicit none
  include "constants.h"

  !input
  integer :: nglob
  real(kind=CUSTOM_REAL), dimension(2,nglob) :: coord

  !output
  real(kind=CUSTOM_REAL), dimension(nglob) :: mask_noise

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
! read noise parameters and check for consitency
  subroutine read_parameters_noise(NOISE_TOMOGRAPHY,SIMULATION_TYPE,SAVE_FORWARD, &
                                     any_acoustic,any_poroelastic,p_sv,irec_master, &
                                     Mxx,Mxz,Mzz,factor,NSOURCES, &
                                     xi_receiver,gamma_receiver,ispec_selected_rec,nrec, &
                                     xi_noise,gamma_noise,ispec_noise,angle_noise)

  implicit none
  include "constants.h"

  !input
  integer :: NOISE_TOMOGRAPHY, SIMULATION_TYPE
  logical :: SAVE_FORWARD
  logical :: any_acoustic, any_poroelastic
  logical :: p_sv
  integer :: NSOURCES, nrec
  double precision, dimension(NSOURCES) :: Mxx, Mxz, Mzz, factor
  double precision, dimension(nrec) :: xi_receiver, gamma_receiver
  integer, dimension(nrec) :: ispec_selected_rec
  double precision :: xi_noise, gamma_noise, angle_noise
  integer :: ispec_noise


  !output
  integer :: irec_master

  !local
  integer :: i,ios

  !define master receiver
  open(unit=509,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/irec_master',status='old',action='read',iostat=ios)
  if (ios == 0) then
    read(509,*) irec_master
  else
    irec_master=1
    write(*,*) ''
    write(*,*) 'Error opening OUTPUT_FILES/NOISE_TOMOGRAPHY/irec_master.'
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


  !check simulation parameters

  if ((NOISE_TOMOGRAPHY/=0) .and. (p_sv)) write(*,*) 'Warning: For P-SV case, noise tomography subroutines not yet fully tested.'
  if (NOISE_TOMOGRAPHY==1) then
     if (SIMULATION_TYPE/=1) call exit_mpi('NOISE_TOMOGRAPHY=1 requires SIMULATION_TYPE=1    -> check DATA/Par_file')



  else if (NOISE_TOMOGRAPHY==2) then
     if (SIMULATION_TYPE/=1) call exit_mpi('NOISE_TOMOGRAPHY=2 requires SIMULATION_TYPE=1    -> check DATA/Par_file')

     if (.not. SAVE_FORWARD) call exit_mpi('NOISE_TOMOGRAPHY=2 requires SAVE_FORWARD=.true.  -> check DATA/Par_file')



  else if (NOISE_TOMOGRAPHY==3) then
     if (SIMULATION_TYPE/=2) call exit_mpi('NOISE_TOMOGRAPHY=3 requires SIMULATION_TYPE=2    -> check DATA/Par_file')

     if (SAVE_FORWARD)       call exit_mpi('NOISE_TOMOGRAPHY=3 requires SAVE_FORWARD=.false. -> check DATA/Par_file')

  endif


!  check model parameters
   if (any_acoustic)    call exit_mpi('Acoustic models not yet implemented for noise simulations. Exiting.')
   if (any_poroelastic) call exit_mpi('Poroelastic models not yet implemented for noise simulations. Exiting.')


!  moment tensor elements must be zero!
   do i=1,NSOURCES
     if ( (Mxx(i) /= 0.d0) .or. (Mxz(i) /= 0.d0) .or. (Mzz(i) /= 0.d0) .or. &
          (factor(i) /= 0.d0)) then
       call exit_mpi('For noise simulations, all moment tensor elements must be zero. Exiting.')
     endif
   enddo


  end subroutine read_parameters_noise


! =============================================================================================================
! read in time series based on noise spectrum and construct noise "source" array
  subroutine compute_source_array_noise(p_sv,NSTEP,deltat,nglob,ibool,ispec_noise, &
                       xi_noise,gamma_noise,xigll,zigll, &
                       time_function_noise,source_array_noise)

  implicit none
  include "constants.h"

  !input parameters
  logical :: p_sv
  integer NSTEP, ispec_noise,nglob
  integer, dimension(NGLLX,NGLLZ,nglob) :: ibool
  real(kind=CUSTOM_REAL) :: deltat, xi_noise, gamma_noise

  !output parameters
  real(kind=CUSTOM_REAL), dimension(NSTEP) :: time_function_noise
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLZ,NSTEP) :: source_array_noise

  !local parameters
  integer :: it,i,j,iglob
  real(kind=CUSTOM_REAL) :: t
  real(kind=CUSTOM_REAL) :: xigll, zigll
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: hxi, hpxi
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: hgamma, hpgamma
  real(kind=CUSTOM_REAL) :: factor_noise, f0, aval, t0

  integer, parameter :: time_function_type = 1


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
! ***SUPPLIED BY THE USER*** IN THE FILE
! "OUTPUT_FILES/NOISE_TOMOGRAPHY/S_squared"
!
! The alternative--setting "time_function_type" to a value other than
! 0--results in an idealized time function that does represent a physically
! realizable noise spectrum but which nevertheless may be useful for
! performing tests or illustrating certain theoretical concepts.
!
! ----------------------------------------------------------------------------------

  time_function_noise(:) = 0._CUSTOM_REAL
  t0   = ((NSTEP-1)/2.)*deltat
  f0   = 15.0d0
  aval = PI*PI*f0*f0
  factor_noise = 1.d3


  if( time_function_type == 1) then
    !Ricker (second derivative of a Gaussian) time function
    do it = 1,NSTEP
      t = it*deltat
      time_function_noise(it) = - factor_noise * 2.*aval * (1. - 2.*aval*(t-t0)**2.) * &
                                 exp(-aval*(t-t0)**2.)
    enddo


  elseif( time_function_type == 2) then
    !first derivative of a Gaussian time function
    do it = 1,NSTEP
      t = it*deltat
      time_function_noise(it) = - factor_noise * (2.*aval*(t-t0)) * exp(-aval*(t-t0)**2.)
    enddo


  elseif( time_function_type == 3) then
    !Gaussian time function
    do it = 1,NSTEP
      t = it*deltat
      time_function_noise(it) = factor_noise * exp(-aval*(t-t0)**2.)
    enddo


  elseif( time_function_type == 4 ) then
    !reproduce time function from Figure 2a of Tromp et al. 2010
    do it = 1,NSTEP
      t = it*deltat
      time_function_noise(it) = factor_noise * &
       4.*aval**2. * (3. - 12.*aval*(t-t0)**2. + 4.*aval**2.*(t-t0)**4.) * &
       exp(-aval*(t-t0)**2.)
    enddo


  elseif( time_function_type == 0 ) then
    !read in time function from file S_squared
    open(unit=55,file='OUTPUT_FILES/NOISE_TOMOGRAPHY/S_squared',status='old')
    do it = 1,NSTEP
      READ(55,*) time_function_noise(it)
    enddo
    close(55)


  else
    call exit_MPI('unknown noise time function')

  endif

  !interpolate over GLL points
  source_array_noise(:,:,:,:) = 0._CUSTOM_REAL
  call lagrange_any(xi_noise,NGLLX,xigll,hxi,hpxi)
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
  subroutine add_point_source_noise(p_sv,it,NSTEP,nglob,ibool,ispec_noise, &
                                   accel_elastic,angle_noise,source_array_noise)
  implicit none
  include "constants.h"

  !input parameters
  logical :: p_sv
  integer :: it, NSTEP
  integer :: ispec_noise, nglob
  integer, dimension(NGLLX,NGLLZ,nglob) :: ibool
  real(kind=CUSTOM_REAL), dimension(3,nglob) :: accel_elastic
  real(kind=CUSTOM_REAL) :: angle_noise
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLZ,NSTEP) :: source_array_noise

  !local variables
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
  subroutine add_surface_movie_noise(p_sv,it,NSTEP,nglob, &
      accel_elastic,surface_movie_x_noise,surface_movie_y_noise, &
      surface_movie_z_noise,mask_noise)

  implicit none
  include "constants.h"

  !input parameters
  logical :: p_sv
  integer :: it, NSTEP
  integer :: nglob

  real(kind=CUSTOM_REAL), dimension(nglob) :: mask_noise
  real(kind=CUSTOM_REAL), dimension(nglob) :: &
    surface_movie_x_noise, surface_movie_y_noise, surface_movie_z_noise
  real(kind=CUSTOM_REAL), dimension(3,nglob) :: accel_elastic

  !local variables
  integer :: ios
  real(kind=CUSTOM_REAL) :: factor_noise
  character(len=60) :: file_in_noise


  write(file_in_noise,"('movie_noise_',i6.6)") NSTEP-it+1
  open(unit=500,file='OUTPUT_FILES/'//trim(file_in_noise),status='old',form='unformatted',action='read',iostat=ios)
  if( ios /= 0) write(*,*) 'Error retrieving generating wavefield.'
  if(p_sv) then
    read(500) surface_movie_x_noise
    read(500) surface_movie_z_noise
  else
    read(500) surface_movie_y_noise
  endif
  close(500)

  factor_noise = 1.d6

  if(p_sv) then ! P-SV calculation
    accel_elastic(1,:) = accel_elastic(1,:) + factor_noise * mask_noise(:) * surface_movie_x_noise(:)
    accel_elastic(3,:) = accel_elastic(3,:) + factor_noise * mask_noise(:) * surface_movie_z_noise(:)
!   note: above implementation does not yet include non-vertical noise
!   excitation

  else ! SH (membrane) calculation
    accel_elastic(2,:) = accel_elastic(2,:) + factor_noise * mask_noise(:) * surface_movie_y_noise(:)

!   note: above implementation includes ad hoc scaling by a number
!   "factor_noise", meant to approximate the correct jacobian and gll weighting
!   correct implementation should look something like:
!   do iglob = 1,nglob
!         hlagrange = weight(iglob) * jacobian(iglob)
!         accel_elastic(2,iglob) = accel_elastic(2,iglob) &
!                    + hlagrange * surface_movie_y_noise(iglob)
!       enddo
!     enddo
!   enddo

  endif

  end subroutine add_surface_movie_noise

! =============================================================================================================
! save a snapshot of the "generating wavefield" that will be used to drive
! the "ensemble forward wavefield"
  subroutine save_surface_movie_noise(p_sv,it,nglob,displ_elastic)

  implicit none
  include "constants.h"

  !input paramters
  logical :: p_sv
  integer :: it, nglob
  real(kind=CUSTOM_REAL), dimension(3,nglob) :: displ_elastic

  !local parameters
  character(len=60) file_out_noise

  write(file_out_noise,"('movie_noise_',i6.6)") it
  open(unit=500,file='OUTPUT_FILES/'//trim(file_out_noise),status='unknown',form='unformatted',action='write')

  if(p_sv) then
    write(500) displ_elastic(1,:)
    write(500) displ_elastic(3,:)
  else
    write(500) displ_elastic(2,:)
  endif

  end subroutine save_surface_movie_noise

! =============================================================================================================
! auxillary routine for writing out the array "mask_noise"; uses a format similar to
! the one currently used for writing kernels
  subroutine write_mask_noise(nglob,coord,mask_noise)

  implicit none
  include "constants.h"

  !input paramters
  integer :: nglob
  real(kind=CUSTOM_REAL), dimension(2,nglob) :: coord
  real(kind=CUSTOM_REAL), dimension(nglob) :: mask_noise

  !local parameters
  integer :: iglob

  open(unit=504,file='OUTPUT_FILES/mask_noise',status='unknown',action='write')
  do iglob = 1, nglob
     write(504,*) coord(1,iglob), coord(2,iglob), mask_noise(iglob)
  enddo
  close(504)


  end subroutine write_mask_noise


  program recombine_all_slices_from_dump

!========================================================================
!
!                   S P E C F E M 2 D  Version 6 . 2
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and Inria, France,
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

! recombine all the MPI mesh slices of a binary dump and remove all the points that are multiples
! i.e. all the points that appear in different files because they belong to the MPI edges

! note: in the case of MPI, in the future it would be more convenient to output a single file
! without multiples directly in the solver rather than one for each myrank

! Dimitri Komatitsch, CNRS Marseille, June 2012

!! DK DK with the Intel ifort compiler, compile with the option below for this code to work fine:

! ifort -o xread -assume byterecl -O3 -xHost recombine_all_slices_from_binary_dump.f90

!  (you will also maybe need   -mcmodel=medium -shared-intel  )

! to debug or test, one can also use:
! -implicitnone -warn truncated_source -warn argument_checking -warn unused -warn declarations -warn alignments
! -warn ignore_loc -warn usage -check all -debug -g -O0 -fp-stack-check -traceback -ftrapuv

  implicit none

! total number of time steps of the simulation to recombine
  integer, parameter :: NSTEP = 5

! interval at which we have dumped the wave field
  integer, parameter :: NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS = 5

! number of processors used to perform the simulation that we now need to recombine
  integer, parameter :: NPROC = 3

! for now I assume that we recombine wave dumps of pressure, i.e. a scalar, not of a vector
! (displacement or velocity), i.e. an array that would have two components; but that
! would be easy to changed if needed
  integer, parameter :: nb_of_values_to_read_back = 1

! type of simulation to recombine: forward or adjoint
  integer, parameter :: SIMULATION_TYPE = 1   ! 1 = forward, 2 = adjoint + kernels

! size of a binary real value in bytes
  integer, parameter :: SIZE_REAL = 4

! very small threshold distance to consider that two points are in fact the same
  real(kind=SIZE_REAL), parameter :: TINYVAL = 1.e-10

  integer, dimension(0:NPROC-1) :: nglob
  real(kind=SIZE_REAL), dimension(:,:), allocatable :: x,y,pressure
  logical, dimension(:,:), allocatable :: this_point_is_a_duplicate

  integer :: it,myrank,iglob,icounter,myrank2,iglob2,number_of_duplicates_found
  integer :: nglob_max_for_all_slices,nglob_recombined_no_duplicates
  real(kind=SIZE_REAL) :: dist

! name of wavefield snapshot file
  character(len=150) :: wavefield_file

! slices are numbered from 0 to NPROC-1
  do myrank = 0,NPROC-1
! read nglob for each mesh slice from a file
        write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_value_of_nglob_',i3.3,'.txt')") myrank
        open(unit=27,file=wavefield_file,status='old',action='read')
        read(27,*) nglob(myrank)
        close(27)
  enddo

! compute the maximum value of nglob for all the slices
  nglob_max_for_all_slices = maxval(nglob)
  print *,'the maximum number of mesh points per slice is ',nglob_max_for_all_slices
  print *

! allocate the arrays based on this maximum size.
! since we used Scotch to decompose the mesh we know that no mesh slice is very significantly
! smaller or bigger than the others and thus doing so is fine
  allocate(x(nglob_max_for_all_slices,0:NPROC-1))
  allocate(y(nglob_max_for_all_slices,0:NPROC-1))
  allocate(pressure(nglob_max_for_all_slices,0:NPROC-1))
  allocate(this_point_is_a_duplicate(nglob_max_for_all_slices,0:NPROC-1))

! clear the arrays
  x(:,:) = 0
  y(:,:) = 0
  pressure(:,:) = 0
  this_point_is_a_duplicate(:,:) = .false.

! slices are numbered from 0 to NPROC-1
  do myrank = 0,NPROC-1
! read the grid for that mesh slice
    print *,'reading grid for mesh slice ',myrank+1,' out of ',NPROC
    write(wavefield_file,"('OUTPUT_FILES/wavefield_grid_for_dumps_',i3.3,'.bin')") myrank
    open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='old', &
                 action='read',recl=2*SIZE_REAL)
    do iglob = 1,nglob(myrank)
      read(27,rec=iglob) x(iglob,myrank),y(iglob,myrank)
    enddo
    close(27)
  enddo

! look for duplicates in the list of points when merging the slices
! and flag them to remove them later
! slices are numbered from 0 to NPROC-1, but we do not look in the first one
! because it cannot have duplicates of itself
  print *
  print *,'looking for duplicates in the merged list of points'
  print *,'and flagging them (can be a slow process for large meshes)'
  do myrank = 1,NPROC-1
    print *,'analyzing slice ',myrank+1,' out of ',NPROC
    do iglob = 1,nglob(myrank)
! look for duplicates in all the previous slices
      do myrank2 = 0,myrank-1
        do iglob2 = 1,nglob(myrank2)
! compute the distance between the two points
          dist = sqrt((x(iglob,myrank)-x(iglob2,myrank2))**2 + (y(iglob,myrank)-y(iglob2,myrank2))**2)
! if the distance is zero (down to roundoff noise) then it is the same point and thus it is a duplicate
          if(dist < TINYVAL) this_point_is_a_duplicate(iglob,myrank) = .true.
        enddo
      enddo
    enddo
  enddo
  number_of_duplicates_found = count(this_point_is_a_duplicate == .true.)
  if(number_of_duplicates_found <= 0) stop 'error: found no duplicates, while there must be some'
  print *
  print *,'total number of duplicates found and removed = ',number_of_duplicates_found

! compute the total number of unique recombined points to write at the end
  nglob_recombined_no_duplicates = sum(nglob) - number_of_duplicates_found
  print *
  print *,'total number of points in the grid recombined from all slices'
  print *,'with no duplicates = ',nglob_recombined_no_duplicates
  print *

! writing the recombined wavefield file with no duplicates
        print *,'writing the recombined grid file in ASCII format'
        open(unit=27,file='OUTPUT_FILES/recombined_grid.txt',status='unknown',action='write')
        icounter = 0
! slices are numbered from 0 to NPROC-1
        do myrank = 0,NPROC-1
          do iglob = 1,nglob(myrank)
            if(.not. this_point_is_a_duplicate(iglob,myrank)) then
              icounter = icounter + 1
              write(27,*) x(iglob,myrank),y(iglob,myrank)
            endif
          enddo
        enddo
        close(27)
        if(icounter /= nglob_recombined_no_duplicates) stop 'error: should have icounter == nglob_recombined_no_duplicates'

  print *
  print *,'Recombining the dumped wave fields from the different files'
  print *,'for the dumped time steps'

! loop on all the time steps
  do it = 1,NSTEP

! determine if a dump exists for this time step
    if(mod(it,NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS) == 0 .or. it == 5 .or. it == NSTEP) then

        print *
        print *,'recombining files for time step ',it
        print *

! slices are numbered from 0 to NPROC-1
        do myrank = 0,NPROC-1
          print *,'reading dumped wavefield for mesh slice ',myrank+1,' out of ',NPROC
          write(wavefield_file,"('OUTPUT_FILES/wavefield',i7.7,'_',i2.2,'_',i3.3,'.bin')") it,SIMULATION_TYPE,myrank
          open(unit=27,file=wavefield_file,form='unformatted',access='direct',status='old', &
                       action='read',recl=nb_of_values_to_read_back*SIZE_REAL)
          do iglob = 1,nglob(myrank)
            read(27,rec=iglob) pressure(iglob,myrank)
          enddo
          close(27)
        enddo

! writing the recombined wavefield file with no duplicates
        print *,'writing the recombined wavefield file in ASCII format'
        write(wavefield_file,"('OUTPUT_FILES/recombined_wavefield_for_time_step_',i7.7,'.txt')") it
        open(unit=27,file=wavefield_file,status='unknown',action='write')
        icounter = 0
! slices are numbered from 0 to NPROC-1
        do myrank = 0,NPROC-1
          do iglob = 1,nglob(myrank)
            if(.not. this_point_is_a_duplicate(iglob,myrank)) then
              icounter = icounter + 1
              write(27,*) pressure(iglob,myrank)
            endif
          enddo
        enddo
        close(27)
        if(icounter /= nglob_recombined_no_duplicates) stop 'error: should have icounter == nglob_recombined_no_duplicates'

    endif

  enddo ! of the loop on all the time steps

  end program recombine_all_slices_from_dump


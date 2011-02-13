
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
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

module source_file
  
  implicit none

  ! source parameters
  integer, dimension(:),pointer ::  source_type,time_function_type  
  double precision, dimension(:),pointer :: xs,zs,f0,t0,angleforce, &
    Mxx,Mzz,Mxz,factor  
  logical, dimension(:),pointer ::  source_surf
    
contains

  subroutine read_source_file(NSOURCE,deltat,f0_attenuation)

! reads in source file DATA/SOURCE

  implicit none
  include "constants.h"
  
  integer :: NSOURCE
  double precision :: deltat,f0_attenuation

  ! local parameters
  integer :: ios,icounter,i_source,nsources
  character(len=150) dummystring
  integer, parameter :: IIN_SOURCE = 22

  ! allocates memory arrays
  allocate(source_surf(NSOURCE))
  allocate(xs(NSOURCE))
  allocate(zs(NSOURCE))
  allocate(source_type(NSOURCE))
  allocate(time_function_type(NSOURCE))
  allocate(f0(NSOURCE))
  allocate(t0(NSOURCE))
  allocate(angleforce(NSOURCE))
  allocate(Mxx(NSOURCE))
  allocate(Mxz(NSOURCE))
  allocate(Mzz(NSOURCE))
  allocate(factor(NSOURCE))

  ! counts lines  
  open(unit=IIN_SOURCE,file='DATA/SOURCE',iostat=ios,status='old',action='read')
  if(ios /= 0) stop 'error opening DATA/SOURCE file'
  
  icounter = 0
  do while(ios == 0)
     read(IIN_SOURCE,"(a)",iostat=ios) dummystring
     if(ios == 0) icounter = icounter + 1
  enddo
  close(IIN_SOURCE)
  
  if(mod(icounter,NLINES_PER_SOURCE) /= 0) &
    stop 'total number of lines in SOURCE file should be a multiple of NLINES_PER_SOURCE'

  nsources = icounter / NLINES_PER_SOURCE
  
  if(nsources < 1) stop 'need at least one source in SOURCE file'
  if(nsources /= NSOURCE) &
       stop 'total number of sources read is different than declared in Par_file'

  ! reads in source parameters
  open(unit=IIN_SOURCE,file='DATA/SOURCE',status='old',action='read')
  do  i_source=1,NSOURCE
    call read_value_logical(IIN_SOURCE,IGNORE_JUNK,source_surf(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,xs(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,zs(i_source))
    call read_value_integer(IIN_SOURCE,IGNORE_JUNK,source_type(i_source))
    call read_value_integer(IIN_SOURCE,IGNORE_JUNK,time_function_type(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,f0(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,t0(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,angleforce(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,Mxx(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,Mzz(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,Mxz(i_source))
    call read_value_double_precision(IIN_SOURCE,IGNORE_JUNK,factor(i_source))

    ! note: this is slightly different than in specfem2D.f90,
    !          t0 will be set outside of this next if statement, i.e. it will be set for all sources
    !          regardless of their type (it just makes a distinction between type 5 sources and the rest)
    
    ! if Dirac source time function, use a very thin Gaussian instead
    ! if Heaviside source time function, use a very thin error function instead
    if(time_function_type(i_source) == 4 .or. time_function_type(i_source) == 5) then
      f0(i_source) = 1.d0 / (10.d0 * deltat)
    endif
    
    ! time delay of the source in seconds, use a 20 % security margin (use 2 / f0 if error function)
    if(time_function_type(i_source)== 5) then
      t0(i_source) = 2.0d0 / f0(i_source) + t0(i_source)
    else
      t0(i_source) = 1.20d0 / f0(i_source) + t0(i_source)
    endif

    print *
    print *,'Source', i_source
    print *,'Position xs, zs = ',xs(i_source),zs(i_source)
    print *,'Frequency, delay = ',f0(i_source),t0(i_source)
    print *,'Source type (1=force, 2=explosion): ',source_type(i_source)
    print *,'Time function type (1=Ricker, 2=First derivative, 3=Gaussian, 4=Dirac, 5=Heaviside): ',time_function_type(i_source)
    print *,'Angle of the source if force = ',angleforce(i_source)
    print *,'Mxx of the source if moment tensor = ',Mxx(i_source)
    print *,'Mzz of the source if moment tensor = ',Mzz(i_source)
    print *,'Mxz of the source if moment tensor = ',Mxz(i_source)
    print *,'Multiplying factor = ',factor(i_source)
    print *
  enddo ! do i_source=1,NSOURCE
  close(IIN_SOURCE)


  ! if source is not a Dirac or Heavyside then f0_attenuation is f0 of the first source
  if(.not. (time_function_type(1) == 4 .or. time_function_type(1) == 5)) then
     f0_attenuation = f0(1)
  endif

  end subroutine read_source_file

end module source_file
  
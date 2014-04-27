
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

  program convolve_source_time_function

!
! convolve seismograms computed for a Heaviside with given source time function
!

! we mimic a triangle of half duration equal to half_duration_triangle
! using a Gaussian having a very close shape, as explained in Figure 4.2
! of the manual

  implicit none

  include "constants.h"

  integer :: i,j,N_j,number_remove,nlines

  double precision :: alpha,dt,tau_j,source,exponentval,t1,t2,displ1,displ2,gamma,height,half_duration_triangle

  logical :: triangle

  double precision, dimension(:), allocatable :: timeval,sem,sem_fil

! read file with number of lines in input
  open(unit=33,file='input_convolve_code.txt',status='old',action='read')
  read(33,*) nlines
  read(33,*) half_duration_triangle
  read(33,*) triangle
  close(33)

! allocate arrays
  allocate(timeval(nlines),sem(nlines),sem_fil(nlines))

! read the input seismogram
  do i = 1,nlines
    read(5,*) timeval(i),sem(i)
  enddo

! define a Gaussian with the right exponent to mimic a triangle of equivalent half duration
  alpha = SOURCE_DECAY_MIMIC_TRIANGLE/half_duration_triangle

! compute the time step
  dt = timeval(2) - timeval(1)

! number of integers for which the source wavelet is different from zero
  if(triangle) then
    N_j = ceiling(half_duration_triangle/dt)
  else
    N_j = ceiling(1.5d0*half_duration_triangle/dt)
  endif

  do i = 1,nlines

    sem_fil(i) = 0.d0

    do j = -N_j,N_j

      if(i > j .and. i-j <= nlines) then

      tau_j = dble(j)*dt

! convolve with a triangle
    if(triangle) then
       height = 1.d0 / half_duration_triangle
       if(abs(tau_j) > half_duration_triangle) then
         source = 0.d0
       else if (tau_j < 0.d0) then
         t1 = - N_j * dt
         displ1 = 0.d0
         t2 = 0.d0
         displ2 = height
         gamma = (tau_j - t1) / (t2 - t1)
         source= (1.d0 - gamma) * displ1 + gamma * displ2
       else
         t1 = 0.d0
         displ1 = height
         t2 = + N_j * dt
         displ2 = 0.d0
         gamma = (tau_j - t1) / (t2 - t1)
         source= (1.d0 - gamma) * displ1 + gamma * displ2
       endif

      else

! convolve with a Gaussian
        exponentval = alpha**2 * tau_j**2
        if(exponentval < 50.d0) then
          source = alpha*exp(-exponentval)/sqrt(PI)
        else
          source = 0.d0
        endif

      endif

      sem_fil(i) = sem_fil(i) + sem(i-j)*source*dt

      endif

    enddo
  enddo

! compute number of samples to remove from end of seismograms
  number_remove = N_j + 1
  do i=1,nlines - number_remove
    write(*,*) sngl(timeval(i)),' ',sngl(sem_fil(i))
  enddo

  end program convolve_source_time_function


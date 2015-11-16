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

! define the forcing applied at the bottom boundary
! programmer Florian Cachoux and Raphael F. Garcia
! in collaboration with D. Komatitsch and R. Martin
! variable forcing_type should be passed as a parameter
! in future versions

  subroutine acoustic_forcing_boundary(iglob)

  use specfem_par

  implicit none

  integer,intent(in) :: iglob

! local variables
  real, parameter :: pigrec = 3.1415927
  real :: alpha,tho,A,c,delayed,delta_x
  integer :: forcing_type,k,ngoce_time_step,n_models,kk,ll

  double precision, dimension(:), allocatable :: goce_time,distance
  double precision, dimension(:,:), allocatable ::syn
  double precision :: t,signal_x1,signal_x2,fracx,fract

  forcing_type = 1


  delta_x = 2.4       ! length of a PML element along x-axis at the edge which will be forced
  alpha = 2

! infrasounds / seismic
  tho = 30.0

! gravity wave test function
!  tho = 600.0
!  xo = 500000.0
!  lambdo = 20000.0

! gravity wave test function
!  tho = 600.0
!  xo = 3000.0
!  lambdo = 120.0

! gravity wave /tsunami
!  tho = 600.0 ! *20
!  c = 200.0 ! /20

  A = 1
  x = coord(1,iglob)
  delayed = 0

! speed of light
  c = 300000000.

  if(forcing_type == 1) then !! First test function : same forcing for the whole boundary
!  print *, ispec_acoustic
!  print *, is_PML(ispec_acoustic)
!  if(is_PML(ispec_acoustic)) then
!  displ_x = 0
!  displ_z = 0
!  else

! infrasounds / seismic
  displ_x = 0 !* Apo
  displ_z = A * (exp(-(alpha*(deltat*it-40-t0)/tho)**2) &
            - exp(-(alpha*(deltat*it-70-t0)/tho)**2)) !* Apo

! gravity wave test function
!  displ_x = 0 !* Apo
!  displ_z = A * ( exp(-(alpha*(x-(xo-lambdo/2))/lambdo)**2) - &
!                  exp(-(alpha*(x-(xo+lambdo/2))/lambdo)**2) ) * &
!            (exp(-(alpha*(deltat*it+1000-t0)/tho)**2) &
!            - exp(-(alpha*(deltat*it-1300-t0)/tho)**2)) !* Apo

! gravity wave /tsunami
!  displ_x = 0 !* Apo
!  displ_z = A * (exp(-(alpha*(deltat*it-1000-t0)/tho)**2) &
!            - exp(-(alpha*(deltat*it-1600-t0)/tho)**2)) !* Apo


  endif

  if(forcing_type == 2) then !! Second test function : moving forcing
  displ_x = 0 !* Apo

  displ_z = A * (exp(-(alpha*(deltat*it-40-t0-(x-delayed)/c)/tho)**2) &
            - exp(-(alpha*(deltat*it-70-t0-(x-delayed)/c)/tho)**2)) !* Apo
  endif

  if(forcing_type == 3) then !! forcing external
    ngoce_time_step = 255
    n_models = 28
    t =it*deltat

    allocate(goce_time(ngoce_time_step))
    allocate(distance(n_models))
    allocate(syn(n_models,ngoce_time_step))

    open(1000,file='../../EXAMPLES/acoustic_forcing_bottom/distance.txt',form='formatted')
    open(1001,file='../../EXAMPLES/acoustic_forcing_bottom/forcing_signals.txt',form='formatted')

    read(1001,*) goce_time(:)

    do k=1,n_models
      read(1001,*) syn(k,:)
      read(1000,*) distance(k)
    enddo

    close(1000)
    close(1001)

    kk = 1
    do while(x >= distance(kk) .and. kk /= n_models)
      kk = kk+1
    enddo

    ll = 1
    do while(t >= goce_time(ll) .and. ll /= ngoce_time_step)
      ll = ll+1
    enddo

      if(x==0 .and. it==1) then
        displ_z =  syn(1,1)
      else
        if(x==0) then
          fract = (t-goce_time(ll-1))/(goce_time(ll)-goce_time(ll-1))
          displ_z =  (syn(1,ll-1) + fract * (syn(1,ll)-syn(1,ll-1)))
        else
          if(it==1) then
            fracx = (x-distance(kk-1))/(distance(kk)-distance(kk-1))
            displ_z =  (syn(kk-1,1) + fracx * (syn(kk,1)-syn(kk-1,1)))
          else
    ! interpolation in time
    fract = (t-goce_time(ll-1))/(goce_time(ll)-goce_time(ll-1))
    ! in x1 = distance(kk-1)
    signal_x1 = syn(kk-1,ll-1) + fract * (syn(kk-1,ll)-syn(kk-1,ll-1))
    ! in x2 = distance(kk)
    signal_x2 = syn(kk,ll-1) + fract * (syn(kk,ll)-syn(kk,ll-1))

    ! spatial interpolation
    fracx = (x-distance(kk-1))/(distance(kk)-distance(kk-1))
    displ_z =  (signal_x1 + fracx * (signal_x2 - signal_x1))
          endif
        endif
      endif

  displ_x = 0

  endif

  if (abs(displ_x) < TINYVAL) displ_x=ZERO
  if (abs(displ_z) < TINYVAL) displ_z=ZERO


  end subroutine acoustic_forcing_boundary

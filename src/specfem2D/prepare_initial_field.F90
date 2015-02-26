
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


  subroutine prepare_initialfield()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par, only: myrank,any_acoustic,any_poroelastic,over_critical_angle, &
                         NSOURCES,source_type,anglesource,x_source,z_source,f0,t0, &
                         nglob,numat,poroelastcoef,density,coord, &
                         anglesource_refl,c_inc,c_refl,cploc,csloc,time_offset, &
                         A_plane, B_plane, C_plane, &
                         accel_elastic,veloc_elastic,displ_elastic

  implicit none
  include "constants.h"

  ! local parameters
  integer :: numat_local,i
  double precision :: denst,lambdaplus2mu,mu,p,x0_source,z0_source
  double precision :: PP,PS,SP,SS,anglesource_abs
  double precision :: xmax, xmin, zmax, zmin,x,z,t
#ifdef USE_MPI
  double precision :: xmax_glob, xmin_glob, zmax_glob, zmin_glob
  integer :: ier
#endif
  double precision, external :: ricker_Bielak_displ,ricker_Bielak_veloc,ricker_Bielak_accel

  ! user output
  if (myrank == 0) then
    write(IOUT,*)
    !! DK DK reading of an initial field from an external file has been suppressed
    !! DK DK and replaced with the implementation of an analytical plane wave
    !! DK DK     write(IOUT,*) 'Reading initial fields from external file...'
    write(IOUT,*) 'Implementing an analytical initial plane wave...'
    write(IOUT,*)
  endif

  if(any_acoustic .or. any_poroelastic) &
    call exit_MPI('initial field currently implemented for purely elastic simulation only')

  !=======================================================================
  !
  !     Calculation of the initial field for a plane wave
  !
  !=======================================================================

  if (myrank == 0) then
    write(IOUT,*) 'Number of grid points: ',nglob
    write(IOUT,*)
    write(IOUT,*) '*** calculation of the initial plane wave ***'
    write(IOUT,*)
    write(IOUT,*)  'To change the initial plane wave, change source_type in DATA/SOURCE'
    write(IOUT,*)  'and use 1 or 4 for a plane P wave, 2 or 5 for a plane SV wave, 3 for a Rayleigh wave'
    write(IOUT,*)

  ! only implemented for one source
    if(NSOURCES > 1) call exit_MPI('calculation of the initial wave is only implemented for one source')
    if (source_type(1) == 1 .or. source_type(1) == 4) then
      write(IOUT,*) 'initial P wave of', anglesource(1)*180.d0/pi, 'degrees introduced.'
    else if (source_type(1) == 2 .or. source_type(1) == 5) then
      write(IOUT,*) 'initial SV wave of', anglesource(1)*180.d0/pi, ' degrees introduced.'
    else if (source_type(1) == 3) then
      write(IOUT,*) 'Rayleigh wave introduced.'
    else
      call exit_MPI('Unrecognized source_type: should be 1 or 4 for plane P waves, 2 or 5 for plane SV waves, 3 for Rayleigh wave')
    endif
  endif

  ! allow negative anglesource(1): incidence from the right side of the domain
  ! anglesource has been converted from degrees to radians before
    anglesource_abs=abs(anglesource(1))
    if (anglesource_abs > pi/2.d0 .and. source_type(1) /= 3) &
      call exit_MPI("incorrect anglesource: must have 0 <= anglesource < 90")

  ! only implemented for homogeneous media therefore only one material supported
  numat_local = numat
  if (numat /= 1) then
     if (myrank == 0) then
        write(IOUT,*)
        write(IOUT,*) 'It is not possible to have several materials with a plane wave, thus using the first material.'
        write(IOUT,*) 'This is not a homogenous model, it contains ',numat,' materials'
        write(IOUT,*) 'but the plane wave initial and boundary fields'
        write(IOUT,*) 'are computed by analytical formulas for a homogenous model.'
        write(IOUT,*) 'Thus use at your own risk!'
        write(IOUT,*)
     endif
     numat_local = 1
  endif

  mu = poroelastcoef(2,1,numat_local)
  lambdaplus2mu  = poroelastcoef(3,1,numat_local)
  denst = density(1,numat_local)

  cploc = sqrt(lambdaplus2mu/denst)
  csloc = sqrt(mu/denst)

  ! P wave case
  if (source_type(1) == 1 .or. source_type(1) == 4) then

    p=sin(anglesource_abs)/cploc
    c_inc  = cploc
    c_refl = csloc

    anglesource_refl = asin(p*c_refl)

    ! from formulas (5.27) and (5.28) p 134 in Aki & Richards (2002)
    PP = (- cos(2.d0*anglesource_refl)**2/csloc**3 &
          + 4.d0*p**2*cos(anglesource_abs)*cos(anglesource_refl)/cploc) / &
               (  cos(2.d0*anglesource_refl)**2/csloc**3 &
                + 4.d0*p**2*cos(anglesource_abs)*cos(anglesource_refl)/cploc)

    PS = 4.d0*p*cos(anglesource_abs)*cos(2.d0*anglesource_refl) / &
               (csloc**2*(cos(2.d0*anglesource_refl)**2/csloc**3 &
               +4.d0*p**2*cos(anglesource_abs)*cos(anglesource_refl)/cploc))

    if (myrank == 0) then
      write(IOUT,*) 'reflected convert plane wave angle: ', anglesource_refl*180.d0/pi
    endif

    ! from Table 5.1 p141 in Aki & Richards (1980)
    ! we put the opposite sign on z coefficients because z axis is oriented from bottom to top
    A_plane(1) = sin(anglesource_abs);           A_plane(2) = cos(anglesource_abs)
    B_plane(1) = PP * sin(anglesource_abs);      B_plane(2) = - PP * cos(anglesource_abs)
    C_plane(1) = PS * cos(anglesource_refl);     C_plane(2) = PS * sin(anglesource_refl)

  ! SV wave case
  else if (source_type(1) == 2 .or. source_type(1) == 5) then

    p=sin(anglesource_abs)/csloc
    c_inc  = csloc
    c_refl = cploc

    ! if this coefficient is greater than 1, we are beyond the critical SV wave angle and there cannot be a converted P wave
    if (p*c_refl<=1.d0) then
      anglesource_refl = asin(p*c_refl)

      ! from formulas (5.30) and (5.31) p 140 in Aki & Richards (1980)
      SS = (cos(2.d0*anglesource_abs)**2/csloc**3 &
          - 4.d0*p**2*cos(anglesource_abs)*cos(anglesource_refl)/cploc) / &
            (cos(2.d0*anglesource_abs)**2/csloc**3 &
              + 4.d0*p**2*cos(anglesource_abs)*cos(anglesource_refl)/cploc)
      SP = 4.d0*p*cos(anglesource_abs)*cos(2*anglesource_abs) / &
            (cploc*csloc*(cos(2.d0*anglesource_abs)**2/csloc**3&
            +4.d0*p**2*cos(anglesource_refl)*cos(anglesource_abs)/cploc))

      if (myrank == 0) then
        write(IOUT,*) 'reflected convert plane wave angle: ', anglesource_refl*180.d0/pi
      endif

    ! SV45 degree incident plane wave is a particular case
    else if (anglesource_abs>pi/4.d0-1.0d-11 .and. anglesource_abs<pi/4.d0+1.0d-11) then
      anglesource_refl = 0.d0
      SS = -1.0d0
      SP = 0.d0
    else
      over_critical_angle=.true.
      anglesource_refl = 0.d0
      SS = 0.0d0
      SP = 0.d0
    endif

    ! from Table 5.1 p141 in Aki & Richards (1980)
    ! we put the opposite sign on z coefficients because z axis is oriented from bottom to top
    A_plane(1) = cos(anglesource_abs);           A_plane(2) = - sin(anglesource_abs)
    B_plane(1) = SS * cos(anglesource_abs);      B_plane(2) = SS * sin(anglesource_abs)
    C_plane(1) = SP * sin(anglesource_refl);     C_plane(2) = - SP * cos(anglesource_refl)

  ! Rayleigh case
  else if (source_type(1) == 3) then
    over_critical_angle=.true.
    A_plane(1)=0.d0; A_plane(2)=0.d0
    B_plane(1)=0.d0; B_plane(2)=0.d0
    C_plane(1)=0.d0; C_plane(2)=0.d0
  endif

   ! correct A_plane, B_plane and C_plane according to incident direction
  if (anglesource(1) < 0.) then
     A_plane(1)=-A_plane(1)
     B_plane(1)=-B_plane(1)
     C_plane(1)=-C_plane(1)
  endif

! to suppress the reflected and converted plane wave fields
  if(source_type(1) == 4 .or. source_type(1) == 5) then
    B_plane(:) = 0
    C_plane(:) = 0
  endif

  ! get minimum and maximum values of mesh coordinates
  xmin = minval(coord(1,:))
  zmin = minval(coord(2,:))
  xmax = maxval(coord(1,:))
  zmax = maxval(coord(2,:))

#ifdef USE_MPI
  call MPI_ALLREDUCE (xmin, xmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (zmin, zmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (xmax, xmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (zmax, zmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  xmin = xmin_glob
  zmin = zmin_glob
  xmax = xmax_glob
  zmax = zmax_glob
#endif

  ! check if zs = zmax (free surface)
  if (myrank == 0 .and. abs(z_source(1)-zmax) > SMALLVALTOL) then
     print *, 'It is sometimes easier to set zs in SOURCE = zmax in interfacefile to keep track of the initial wavefront'
  endif

  ! add -t0 to match the actual (analytical) traveltime of plane waves
  time_offset = 0.d0 - t0

  ! to correctly center the initial plane wave in the mesh
  x0_source = x_source(1)
  z0_source = z_source(1)

  if (myrank == 0) then
    write(IOUT,*)
    write(IOUT,*) 'You can modify the location of the initial plane wave by changing xs and zs in DATA/SOURCE.'
    write(IOUT,*) '   for instance: xs=',x_source(1),'   zs=',z_source(1), ' (zs can/should be the height of the free surface)'
    write(IOUT,*)
  endif

  if (.not. over_critical_angle) then

    do i = 1,nglob

      x = coord(1,i)
      z = coord(2,i)

      ! z is from bottom to top therefore we take -z to make parallel with Aki & Richards

      z = z0_source - z
      if (anglesource(1) >= 0.) then
         x = x - x0_source
      else
         x = x0_source -x
      endif

      t = 0.d0 + time_offset

      ! formulas for the initial displacement for a plane wave from Aki & Richards (1980)
      displ_elastic(1,i) = &
          A_plane(1) * ricker_Bielak_displ(t - sin(anglesource_abs)*x/c_inc + cos(anglesource_abs)*z/c_inc,f0(1)) &
        + B_plane(1) * ricker_Bielak_displ(t - sin(anglesource_abs)*x/c_inc - cos(anglesource_abs)*z/c_inc,f0(1)) &
        + C_plane(1) * ricker_Bielak_displ(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0(1))
      displ_elastic(3,i) = &
          A_plane(2) * ricker_Bielak_displ(t - sin(anglesource_abs)*x/c_inc + cos(anglesource_abs)*z/c_inc,f0(1)) &
        + B_plane(2) * ricker_Bielak_displ(t - sin(anglesource_abs)*x/c_inc - cos(anglesource_abs)*z/c_inc,f0(1)) &
        + C_plane(2) * ricker_Bielak_displ(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0(1))

      ! formulas for the initial velocity for a plane wave (first derivative in time of the displacement)
      veloc_elastic(1,i) = &
          A_plane(1) * ricker_Bielak_veloc(t - sin(anglesource_abs)*x/c_inc + cos(anglesource_abs)*z/c_inc,f0(1)) &
        + B_plane(1) * ricker_Bielak_veloc(t - sin(anglesource_abs)*x/c_inc - cos(anglesource_abs)*z/c_inc,f0(1)) &
        + C_plane(1) * ricker_Bielak_veloc(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0(1))
      veloc_elastic(3,i) = &
          A_plane(2) * ricker_Bielak_veloc(t - sin(anglesource_abs)*x/c_inc + cos(anglesource_abs)*z/c_inc,f0(1)) &
        + B_plane(2) * ricker_Bielak_veloc(t - sin(anglesource_abs)*x/c_inc - cos(anglesource_abs)*z/c_inc,f0(1)) &
        + C_plane(2) * ricker_Bielak_veloc(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0(1))

      ! formulas for the initial acceleration for a plane wave (second derivative in time of the displacement)
      accel_elastic(1,i) = &
          A_plane(1) * ricker_Bielak_accel(t - sin(anglesource_abs)*x/c_inc + cos(anglesource_abs)*z/c_inc,f0(1)) &
        + B_plane(1) * ricker_Bielak_accel(t - sin(anglesource_abs)*x/c_inc - cos(anglesource_abs)*z/c_inc,f0(1)) &
        + C_plane(1) * ricker_Bielak_accel(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0(1))
      accel_elastic(3,i) = &
          A_plane(2) * ricker_Bielak_accel(t - sin(anglesource_abs)*x/c_inc + cos(anglesource_abs)*z/c_inc,f0(1)) &
        + B_plane(2) * ricker_Bielak_accel(t - sin(anglesource_abs)*x/c_inc - cos(anglesource_abs)*z/c_inc,f0(1)) &
        + C_plane(2) * ricker_Bielak_accel(t - sin(anglesource_refl)*x/c_refl - cos(anglesource_refl)*z/c_refl,f0(1))

   enddo

endif

end subroutine prepare_initialfield

!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_initialfield_paco()

  use specfem_par, only: myrank,nelemabs,left_bound,right_bound,bot_bound, &
                                    numabs,codeabs,ibool, &
                                    source_type,c_inc,c_refl, &
                                    count_bottom,count_left,count_right

  implicit none
  include "constants.h"

  ! local parameters
  integer :: ispecabs,ispec,i,j,iglob,ibegin,iend

  if (myrank == 0) then
    if (source_type(1) /= 3 ) &
      write(IOUT,*) 'You are beyond the critical angle ( > ',asin(c_inc/c_refl)*180d0/pi,')'

    write(IOUT,*)  '*************'
    write(IOUT,*)  'We have to compute the initial field in the frequency domain'
    write(IOUT,*)  'and then convert it to the time domain (can be long... be patient...)'
    write(IOUT,*)  '*************'
  endif

  count_bottom=0
  count_left=0
  count_right=0
  do ispecabs=1,nelemabs
    ispec=numabs(ispecabs)
    if(codeabs(IEDGE4,ispecabs)) then
       i = 1
       do j = 1,NGLLZ
          count_left=count_left+1
          iglob = ibool(i,j,ispec)
          left_bound(count_left)=iglob
       enddo
    endif
    if(codeabs(IEDGE2,ispecabs)) then
       i = NGLLX
       do j = 1,NGLLZ
          count_right=count_right+1
          iglob = ibool(i,j,ispec)
          right_bound(count_right)=iglob
       enddo
    endif
    if(codeabs(IEDGE1,ispecabs)) then
       j = 1
!! DK DK not needed       ! exclude corners to make sure there is no contradiction regarding the normal
       ibegin = 1
       iend = NGLLX
!! DK DK not needed       if(codeabs(IEDGE4,ispecabs)) ibegin = 2
!! DK DK not needed       if(codeabs(IEDGE2,ispecabs)) iend = NGLLX-1
       do i = ibegin,iend
          count_bottom=count_bottom+1
          iglob = ibool(i,j,ispec)
          bot_bound(count_bottom)=iglob
       enddo
    endif
  enddo

  end subroutine prepare_initialfield_paco


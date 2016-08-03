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

! ----------------------------------------------------------------------------------------
! This file contains all the subroutines used to enforce a displacement at
! a given position.
! !!! Warning given lines need double precision (look for "warning") !!! TODO
! ----------------------------------------------------------------------------------------

module enforce_par
! ----------------------------------------------------------------------------------------
! Contains the variables needed for these functions
! ----------------------------------------------------------------------------------------

  implicit none

  ! Line which is forced
  double precision :: xforced = 3000.d0
  integer :: nforced = 0, nforced_sum = 0

end module enforce_par

!
! ----------------------------------------------------------------------------------------
!

  subroutine build_forced()
! ----------------------------------------------------------------------------------------
! This subroutine build the logical array: is_forced which gives for each GLL
! point if we impose the displacement on this point or not.
! For example: iglob_is_forced(iglob) = .false.
! forced has already been already initialized to .false.
! ----------------------------------------------------------------------------------------

  use specfem_par, only: coord,nglob,ibool,iglob_is_forced,myrank,nelem_acforcing,codeacforcing,numacforcing,ibool, &
                         PML_BOUNDARY_CONDITIONS,ispec_is_PML,read_external_mesh
  use enforce_par
  use constants, only: TINYVAL,IMAIN,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  implicit none

  !local variables
  integer :: inum,ispec,i,j,iglob

  if (read_external_mesh) then

    ! loop on all the forced edges
    do inum = 1,nelem_acforcing
      ispec = numacforcing(inum)
      !--- left acoustic forcing boundary
      if (codeacforcing(IEDGE4,inum)) then
        i = 1
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          if (PML_BOUNDARY_CONDITIONS) then
            if (.not. ispec_is_PML(ispec)) then
              iglob_is_forced(iglob) = .true.
              nforced = nforced + 1
            endif
          else
            iglob_is_forced(iglob) = .true.
            nforced = nforced + 1
          endif
        enddo
      endif  !  end of left acoustic forcing boundary
      !--- right acoustic forcing boundary
      if (codeacforcing(IEDGE2,inum)) then
        i = NGLLX
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          if (PML_BOUNDARY_CONDITIONS) then
            if (.not. ispec_is_PML(ispec)) then
              iglob_is_forced(iglob) = .true.
              nforced = nforced + 1
            endif
          else
            iglob_is_forced(iglob) = .true.
            nforced = nforced + 1
          endif
        enddo
      endif  !  end of right acoustic forcing boundary
      !--- bottom acoustic forcing boundary
      if (codeacforcing(IEDGE1,inum)) then
        j = 1
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if (PML_BOUNDARY_CONDITIONS) then
            if (.not. ispec_is_PML(ispec)) then
              iglob_is_forced(iglob) = .true.
              nforced = nforced + 1
            endif
          else
            iglob_is_forced(iglob) = .true.
            nforced = nforced + 1
          endif
        enddo
      endif ! end of bottom acoustic forcing boundary
      !--- top acoustic forcing boundary
      if (codeacforcing(IEDGE3,inum)) then
        j = NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if (PML_BOUNDARY_CONDITIONS) then
            if (.not. ispec_is_PML(ispec)) then
              iglob_is_forced(iglob) = .true.
              nforced = nforced + 1
            endif
          else
            iglob_is_forced(iglob) = .true.
            nforced = nforced + 1
          endif
        enddo
      endif  !  end of top acoustic forcing boundary
    enddo

  else
    ! internal mesher:
    ! Loop on all the GLL points
    do iglob = 1,nglob
      ! forced points along a specified x-line
      if (abs(coord(1,iglob) - xforced) < TINYVAL) then
        iglob_is_forced(iglob) = .true.
        nforced = nforced + 1
      endif
    enddo

  endif

  ! user output
  call sum_all_i(nforced, nforced_sum)
  if (myrank == 0) then
    if (nforced_sum == 0) then
      call exit_MPI(myrank,'No forced integration point found!')
    else
      write(IMAIN,*) "Number of GLL points forced:",nforced," over ",nglob
      call flush_IMAIN()
    endif
  endif

  end subroutine build_forced

!
! ----------------------------------------------------------------------------------------
!

  subroutine enforce_fields(iglob,it)
! ----------------------------------------------------------------------------------------
! This subroutine impose the fields at given GLL points and at a given time steps
! ----------------------------------------------------------------------------------------

  use specfem_par, only: coord,displ_elastic,veloc_elastic,accel_elastic,deltat,deltatover2, &
                         deltatsquareover2
  use enforce_par
  use constants, only: TINYVAL,CUSTOM_REAL,TWO,PI

  implicit none

  ! Inputs
  integer, intent(in) :: iglob,it

  ! Local variables
  real(kind=CUSTOM_REAL) :: f0 = 0.5236d6 ! frequency
  real(kind=CUSTOM_REAL) :: h = 1.0d-3 ! half width of the plate
  real(kind=CUSTOM_REAL) :: cp = 5960.0d0 ! Compressional waves velocity
  real(kind=CUSTOM_REAL) :: cs = 3260.d0 ! Shear waves velocity
  real(kind=CUSTOM_REAL), dimension(2) :: accelOld,velocOld
  real(kind=CUSTOM_REAL) :: factor,x,z,facx,facz,tval,tval_old
  real(kind=CUSTOM_REAL) :: omegaj
  real(kind=CUSTOM_REAL) :: time_dependence_x,time_dependence_x_old,time_dependence_z,time_dependence_z_old
  complex :: ux,uz
  logical :: antisym = .true.
  integer :: Nc = 10

  factor = 1.0d0
  omegaj = TWO*PI*f0

  ! gets global node location
  x = coord(1,iglob)
  z = coord(2,iglob)

  ! timing
  tval = (it-1)*deltat
  tval_old = (it-2)*deltat

  if (tval < TWO*Nc*PI/omegaj) then
    time_dependence_x = omegaj**2*cos(omegaj*tval/Nc)*sin(omegaj*tval)/Nc**2 + &
                        2*omegaj**2*sin(omegaj*tval/Nc)*cos(omegaj*tval)/Nc - &
                        (1.0d0 - cos(omegaj*tval/Nc))*omegaj**2*sin(omegaj*tval)
    time_dependence_z = omegaj**2*cos(omegaj*tval/Nc)*cos(omegaj*tval)/Nc**2 - &
                        2*omegaj**2*sin(omegaj*tval/Nc)*sin(omegaj*tval)/Nc - &
                        (1.0d0 - cos(omegaj*tval/Nc))*omegaj**2*cos(omegaj*tval)
  else
    time_dependence_x = 0.0d0
    time_dependence_z = 0.0d0
  endif
  if (tval_old < TWO*Nc*PI/omegaj) then
    time_dependence_x_old = omegaj**2*cos(omegaj*tval_old/Nc)*sin(omegaj*tval_old)/Nc**2 + &
                        2*omegaj**2*sin(omegaj*tval_old/Nc)*cos(omegaj*tval_old)/Nc - &
                        (1.0d0 - cos(omegaj*tval_old/Nc))*omegaj**2*sin(omegaj*tval_old)

    time_dependence_z_old = omegaj**2*cos(omegaj*tval_old/Nc)*cos(omegaj*tval_old)/Nc**2 - &
                        2*omegaj**2*sin(omegaj*tval_old/Nc)*sin(omegaj*tval_old)/Nc - &
                        (1.0d0 - cos(omegaj*tval_old/Nc))*omegaj**2*cos(omegaj*tval_old)
  else
    time_dependence_x_old = 0.0d0
    time_dependence_z_old = 0.0d0
  endif

  call calculateUxUz(ux,uz,z,cp,cs,h,omegaj,antisym)

  if ((abs(real(ux)) < TINYVAL) .and. (abs(aimag(ux)) > TINYVAL)) then ! if ux is imaginary
    !print *,"Ux imaginary"
    if (abs(aimag(uz)) > TINYVAL) then ! ... and uz not real
      print *,"Problem!! ux is imaginary while uz is not real (ux:",ux," uz:",uz,")"
    endif
    facx = aimag(ux)
    facz = real(uz)
  else if ((abs(real(ux)) > TINYVAL) .and. (abs(aimag(ux)) < TINYVAL)) then  ! if ux is real
    !print *,"Ux real"
    if (abs(real(uz)) > TINYVAL) then ! ... and uz not imaginary
      print *,"Problem!! ux is real while uz is not imaginary (ux:",ux," uz:",uz,")"
    endif
    facx = real(ux)
    facz = aimag(uz)
  else
    facx = real(ux)
    facz = aimag(uz)
  endif

  !print *,z,facz

  if (abs(z + 0.0d-3) < 1.0d-3) then
    if (it == 1) then
      ! We initialize the variables
      displ_elastic(1,iglob) = factor*facx*0.0d0
      displ_elastic(2,iglob) = factor*facz*0.0d0
      veloc_elastic(1,iglob) = factor*facx*0.0d0
      veloc_elastic(2,iglob) = factor*facz*0.0d0
      accel_elastic(1,iglob) = factor*facx*0.0d0
      accel_elastic(2,iglob) = factor*facz*(omegaj**2/Nc**2)
    else
      ! We set what we want
      accel_elastic(1,iglob) = factor*facx*time_dependence_x
      accelOld(1) = factor*facx*time_dependence_x_old
      accel_elastic(2,iglob) = factor*facz*time_dependence_z
      accelOld(2) = factor*facz*time_dependence_z_old
      ! Do not change anything below: we compute numerically the velocity and displacement
      velocOld(1) = veloc_elastic(1,iglob)
      veloc_elastic(1,iglob) = veloc_elastic(1,iglob) + deltatover2*(accelOld(1) + accel_elastic(1,iglob))
      displ_elastic(1,iglob) = displ_elastic(1,iglob) + deltat*velocOld(1) + deltatsquareover2*accelOld(1)
      velocOld(2) = veloc_elastic(2,iglob)
      veloc_elastic(2,iglob) = veloc_elastic(2,iglob) + deltatover2*(accelOld(2) + accel_elastic(2,iglob))
      displ_elastic(2,iglob) = displ_elastic(2,iglob) + deltat*velocOld(2) + deltatsquareover2*accelOld(2)
    endif
  else
    displ_elastic(:,iglob) = 0.0d0
    veloc_elastic(:,iglob) = 0.0d0
    accel_elastic(:,iglob) = 0.0d0
  endif

  end subroutine enforce_fields

!
! ----------------------------------------------------------------------------------------
!

  subroutine calculateUxUz(ux,uz,zi,cp,cs,h,omegaj,antisym)
  ! ----------------------------------------------------------------------------------------
  ! See eq. (3) (4) Hora & Cervena 2012
  ! ----------------------------------------------------------------------------------------

  use constants, only: CUSTOM_REAL,TWO,PI

  implicit none

  ! Inputs
  logical, intent(in) :: antisym
  real(kind=CUSTOM_REAL), intent(in) :: zi,cp,cs,h,omegaj

  ! Outputs
  complex, intent(out) :: ux,uz

  ! Local variables
  real(kind=CUSTOM_REAL) :: cphase
  complex :: alpha,beta,k,k2,alpha2,beta2
  complex :: jj = (0.0,1.0)

  call calculateCphase(cphase,omegaj,antisym,h)
  k2 = (omegaj/cphase)**2
  alpha2 = (omegaj/cp)**2 - k2
  beta2 = (omegaj/cs)**2 - k2
  beta = sqrt(beta2)
  alpha = sqrt(alpha2)
  k = sqrt(k2)
  if (.not. antisym) then
    ux = jj * beta * (cos(beta*zi)/cos(beta*h) - TWO*k2*cos(alpha*zi)/((k2 - beta2)*cos(alpha*h)))
    uz = - k * (sin(beta*zi)/cos(beta*h) + TWO*alpha*beta*sin(alpha*zi)/((k2 - beta2)*cos(alpha*h)))
  else
    if (.false.) then ! Hora
      ux = - jj * beta * (sin(beta*zi)/sin(beta*h) - TWO*k2*sin(alpha*zi)/((k2 - beta2 )*sin(alpha*h)))
      uz = - k * (cos(beta*zi)/sin(beta*h) + TWO*alpha*beta*cos(alpha*zi)/((k2 - beta2 )*sin(alpha*h)))
    else ! Alleyne
      ux = jj * beta * (sin(beta*zi)/sin(beta*h) - TWO*k2*sin(alpha*zi)/((k2 - beta2)*sin(alpha*h)))
      uz = - k * (cos(beta*zi)/sin(beta*h) + TWO*alpha*beta*cos(alpha*zi)/((k2 - beta2)*sin(alpha*h)))
    endif
  endif

  end subroutine calculateUxUz

!
! ----------------------------------------------------------------------------------------
!

  subroutine calculateCphase(cphase,omegaj,antisym,h) !,cp,cs,h)
  ! ----------------------------------------------------------------------------------------
  ! This is not trivial !! We use here an approximate value
  ! do something clever here
  ! ----------------------------------------------------------------------------------------

  use specfem_par, only: myrank
  use constants, only: TINYVAL,CUSTOM_REAL,TWO,PI

  implicit none

  ! Inputs
  logical, intent(in) :: antisym
  real(kind=CUSTOM_REAL), intent(in) :: omegaj,h !,cp,cs

  ! Outputs
  real(kind=CUSTOM_REAL), intent(out) :: cphase

  ! Local variables
  real(kind=CUSTOM_REAL) :: calculated1 = 0.5d6 * 1.5d-3
  real(kind=CUSTOM_REAL) :: calculated2 = 2.087d6 * 1.0d-3
  real(kind=CUSTOM_REAL) :: calculated3 = 0.5236d6 * 1.0d-3
  real(kind=CUSTOM_REAL) :: calculated4 = 2.3113d6 * 1.0d-3
  real(kind=CUSTOM_REAL) :: freq

  freq = omegaj / (TWO*PI)

  if (abs(freq*h - calculated1) < TINYVAL) then ! warning this line need double precision
    if (antisym) then
      cphase = 2611.4d0 ! A0 at fd = 0.75 Mhz.mm
    else
      cphase = 5300.0d0 ! S0 at fd = 0.75 Mhz.mm
    endif
  else if (abs(freq*h - calculated2) < TINYVAL) then ! warning this line need double precision
    if (antisym) then
      cphase = 5234.0d0 ! A1 at fd = 2.087 Mhz.mm
    else
      cphase = 5883.0d0 ! S1 at fd = 2.087 Mhz.mm
    endif
  else if (abs(freq*h - calculated3) < TINYVAL) then ! warning this line need double precision
    if (antisym) then
      cphase = 1937.0d0 ! A0 at fd = 0.5236 Mhz.mm
    else
      cphase = 5445.0d0 ! S0 at fd = 0.5236 Mhz.mm
    endif
  else if (abs(freq*h - calculated4) < TINYVAL) then ! warning this line need double precision
    if (antisym) then
      cphase = 4595.0d0 ! A1 at fd = 2.3113 Mhz.mm
    else
      cphase = 5780.0d0 ! S1 at fd = 2.3113 Mhz.mm
    endif
  else
    print *,"f*d = ",freq*h,"Hz*m",calculated3
    call exit_MPI(myrank,'Phase speed not calculated for this f*d product')
  endif

  end subroutine calculateCphase

! subroutine enforce_fields(iglob,it)
!! ----------------------------------------------------------------------------------------
!! This subroutine impose the fields at a given GLL point and at a given time step
!! ----------------------------------------------------------------------------------------

!  use specfem_par, only: coord,forced,displ_elastic,veloc_elastic,accel_elastic,deltat,deltatover2, &
!                         deltatsquareover2
!  use enforce_par
!  use constants, only: TINYVAL,CUSTOM_REAL,TWO,PI

!  implicit none
!
!  ! Inputs
!  integer, intent(in) :: iglob,it

!  ! Local variables

!  ! Local variables
!  real(kind=CUSTOM_REAL), dimension(2) :: accelOld,velocOld
!  real(kind=CUSTOM_REAL) :: factor,x,z
!  real(kind=CUSTOM_REAL) :: f0 = 0.5d6

!  x = coord(1,iglob)
!  z = coord(2,iglob)
!
!  factor = 1.0d0 ! * (z - 2.0d0)**2
!
!  !if (abs(z + 1.5d-3) < 1.5d-3) then
!    if (it == 1) then ! We initialize the variables
!      displ_elastic(1,iglob) = factor*(-1.0d0)
!      displ_elastic(2,iglob) = factor*0.0d0
!      veloc_elastic(1,iglob) = factor*0.0d0
!      veloc_elastic(2,iglob) = factor*0.0d0
!      accel_elastic(1,iglob) = factor*(TWO*PI*f0)**2
!      accel_elastic(2,iglob) = factor*0.0d0
!    else ! We set what we want
!      accel_elastic(1,iglob) = factor*(TWO*PI*f0)**2*cos(TWO*PI*f0*(it-1)*deltat)
!      accelOld(1) = factor*(TWO*PI*f0)**2*cos(TWO*PI*f0*(it-2)*deltat)
!      accel_elastic(2,iglob) = 0.0d0
!      accelOld(2) = 0.0d0
!      ! Do not change anything below: we compute numerically the velocity and displacement
!      velocOld(1) = veloc_elastic(1,iglob)
!      veloc_elastic(1,iglob) = veloc_elastic(1,iglob) + deltatover2*(accelOld(1) + accel_elastic(1,iglob))
!      displ_elastic(1,iglob) = displ_elastic(1,iglob) + deltat*velocOld(1) + deltatsquareover2*accelOld(1)
!      velocOld(2) = veloc_elastic(2,iglob)
!      veloc_elastic(2,iglob) = veloc_elastic(2,iglob) + deltatover2*(accelOld(2) + accel_elastic(2,iglob))
!      displ_elastic(2,iglob) = displ_elastic(2,iglob) + deltat*velocOld(2) + deltatsquareover2*accelOld(2)
!    endif
!  !else
!  !  displ_elastic(:,iglob) = 0.0d0
!  !  veloc_elastic(:,iglob) = 0.0d0
!  !  accel_elastic(:,iglob) = 0.0d0
!  !endif

! end subroutine enforce_fields


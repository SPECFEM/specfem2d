!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
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

! ----------------------------------------------------------------------------------------
! This file contains all the subroutines used to enforce a displacement at
! a given position.
! !!! Warning given lines need double precision (look for "warning") !!! TODO
! ----------------------------------------------------------------------------------------

module enforce_par
! ----------------------------------------------------------------------------------------
! Contains the variables needed for these functions
! ----------------------------------------------------------------------------------------

  use constants, only: CUSTOM_REAL

  implicit none

  ! Line which is forced
  double precision :: xforced = 3000.d0
  integer :: nforced = 0, nforced_sum = 0
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: cphaseVec, fdVec
  logical :: hasReadFile = .false. ! This is put to true when the file has been read
  integer :: indexFdIn ! Index of the closest fd point in file
  integer :: number_of_lines_in_file = -1
  logical :: neverRead = .true.
  integer :: num_file = 7771 ! unit number for opening files

  contains

  function count_lines(filename) result(nlines)
    implicit none
    character(len=*)    :: filename
    integer             :: nlines
    integer             :: io

    open(10,file=filename, iostat=io, status='old')
    if (io /= 0) call stop_the_code('Cannot open file! ')

    nlines = 0
    do
      read(10,*,iostat=io)
      if (io /= 0) exit
      nlines = nlines + 1
    enddo
    close(10)
  end function count_lines

end module enforce_par

!
! ----------------------------------------------------------------------------------------
!

  subroutine build_forced()
! ----------------------------------------------------------------------------------------
! This subroutine build the logical array: is_forced which gives for each GLL
! point if we impose the displacement on this point or not.
! For example: iglob_is_forced(iglob) = .false.
! forced has already been initialized to .false.
! ----------------------------------------------------------------------------------------

  use specfem_par, only: coord,nglob,ibool,iglob_is_forced,myrank,nelem_acforcing,codeacforcing,numacforcing,ibool, &
                         PML_BOUNDARY_CONDITIONS,ispec_is_PML,read_external_mesh,acoustic_iglob_is_forced, &
                         elastic_iglob_is_forced,ispec_is_acoustic,ispec_is_elastic
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
    do inum = 1,nelem_acforcing
      ispec = numacforcing(inum)
      do i = 1,NGLLX
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          if (iglob_is_forced(iglob)) then
            if (ispec_is_acoustic(ispec)) then
               acoustic_iglob_is_forced = .true.
            else if (ispec_is_elastic(ispec)) then
               elastic_iglob_is_forced = .true.
            else
              call stop_the_code("Problem forced! Element is not acoustic and not elastic... this is not possible so far")
            endif
          endif
        enddo
      enddo
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
      write(IMAIN,*) "Number of GLL points forced:",nforced_sum," over ",nglob
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
                         deltatsquareover2,myrank
  use enforce_par
  use constants, only: TINYVAL,CUSTOM_REAL,TWO,PI

  implicit none

  ! Inputs
  integer, intent(in) :: iglob,it

  ! Local variables
  real(kind=CUSTOM_REAL) :: f0 = 0.125d6 ! frequency ! (fd=200,f=50KHz) (fd=500,f=125KHz) (fd=800,f=200KHz)
  real(kind=CUSTOM_REAL) :: d = 4.0d-3 ! half width of the plate
  real(kind=CUSTOM_REAL) :: cp = 5960.0d0 ! Compressional waves velocity
  real(kind=CUSTOM_REAL) :: cs = 3260.d0 ! Shear waves velocity
  real(kind=CUSTOM_REAL) :: fdin ! = f0 * d

  integer :: Nweight = 11 ! Number of point on which we perform the weighted sum. Must be odd
  real(kind=CUSTOM_REAL), dimension(2) :: accelOld,velocOld
  real(kind=CUSTOM_REAL) :: factor,x,z,facx,facz,tval,tval_old
  real(kind=CUSTOM_REAL) :: omegaj
  real(kind=CUSTOM_REAL) :: time_dependence_x,time_dependence_x_old,time_dependence_z,time_dependence_z_old
  complex(CUSTOM_REAL) :: sum_ux,sum_uz

  ! Choose the mode to generate
  logical :: antisym = .false.
  ! String that must be equal to 0,1, and so on
  character (len=1) :: order='0'
  ! Number of cycle of the burst
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

  fdin = f0 * d
  ! Read file and calculate weighted sum.
  if (.not. hasReadFile) then
    call readCphaseAndFdInFile(antisym,order,fdin)
    hasReadFile = .true.
  endif

  !print *,fdVec(indexFdIn),cphaseVec(indexFdIn)
  call weighted_sum_Lamb_disp(sum_ux,sum_uz,z,f0,d,cp,cs,antisym,Nc,Nweight)

  if ((abs(real(sum_ux)) < TINYVAL) .and. (abs(aimag(sum_ux)) > TINYVAL)) then ! if sum_ux is imaginary
    !print *,"sum_ux imaginary"
    if (abs(aimag(sum_uz)) > TINYVAL) then ! ... and sum_uz not real
      print *,"Problem!! sum_ux is imaginary while sum_uz is not real (sum_ux:",sum_ux," sum_uz:",sum_uz,")"
      call exit_MPI(myrank,"Problem!! sum_ux is imaginary while sum_uz is not real")
    endif
    facx = aimag(sum_ux)
    facz = real(sum_uz)
  else if ((abs(real(sum_ux)) > TINYVAL) .and. (abs(aimag(sum_ux)) < TINYVAL)) then  ! if sum_ux is real
    !print *,"sum_ux real"
    if (abs(real(sum_uz)) > TINYVAL) then ! ... and uz not imaginary
      print *,"Problem!! sum_ux is real while uz is not imaginary (sum_ux:",sum_ux," sum_uz:",sum_uz,")"
      call exit_MPI(myrank,"Problem!! sum_ux is real while sum_uz is not imaginary")
    endif
    facx = real(sum_ux)
    facz = aimag(sum_uz)
  else ! In practice here is the case where the displacement is 0
    facx = real(sum_ux)
    facz = aimag(sum_uz)
  endif

  !print *,z,facz

!  if (abs(z + 0.0d-3) < 1.0d-3) then
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
!  else
!    displ_elastic(:,iglob) = 0.0d0
!    veloc_elastic(:,iglob) = 0.0d0
!    accel_elastic(:,iglob) = 0.0d0
!  endif

  end subroutine enforce_fields

!
! ----------------------------------------------------------------------------------------
!

  subroutine weighted_sum_Lamb_disp(sum_ux,sum_uz,z,f0,d,cp,cs,antisym,Nc,Nweight)

   use constants, only: TINYVAL,CUSTOM_REAL,TWO,PI
   use enforce_par, only: fdVec,cphaseVec,indexFdIn

   implicit none

    ! imput variables
    real(kind=CUSTOM_REAL), intent(in) :: f0    ! frequency ! (fd=500,f=125KHz)
    real(kind=CUSTOM_REAL), intent(in) :: d     ! half width of the plate
    real(kind=CUSTOM_REAL), intent(in) ::  cp   ! Compressional waves velocity
    real(kind=CUSTOM_REAL), intent(in) ::  cs   ! Shear waves velocity
    !real(kind=CUSTOM_REAL), intent(in) :: fdin ! Frequencies width product around which we calculate (f0 * d)
    integer , intent(in) :: Nweight             ! Number of elements in the sum

    ! choose the mode to generate
    ! real(kind=CUSTOM_REAL) , intent(in):: fdin ! KHz.mm
    logical , intent(in):: antisym
    integer , intent(in) :: Nc                   !number of cycle of the burst
    real(kind=CUSTOM_REAL), intent(in) :: z      ! coord z forced

    ! Output variables
    complex(kind=CUSTOM_REAL), intent(out) :: sum_ux,sum_uz

    ! local variables
    complex(kind=CUSTOM_REAL) :: ux=(0.0,0.0),uz=(0.0,0.0)

    real(kind=CUSTOM_REAL) ::  fc!,fdmin,fdmax,stepfd
    real(kind=CUSTOM_REAL) ::  freq,omegaj!,fmax,fmin
    real(kind=CUSTOM_REAL) ::  cphase
    !real(kind=CUSTOM_REAL) ::  BW ! frequency band width
    real(kind=CUSTOM_REAL) ::  Weigth_Burst,sum_Weight ! frequency weight
    integer :: iweight ! iterator over weigths

    omegaj=TWO*PI*f0

    !DSP = weight = 0 at fmin and fmax
    ! this frequency range corresponds to the principal lobe of the DSP

    fc = fdVec(indexFdIn)/d
!    BW=2.0d0*fc/Nc
!    fmin=fc-BW; fmax=fc+BW
!    fdmin=fmin*d; fdmax=fmax*d
!    stepfd=(fdmax-fdmin)/(Nweight-1)

    sum_ux = (0.0,0.0)
    sum_uz = (0.0,0.0)
    sum_Weight = 0

    do iweight=indexFdIn-Nweight,indexFdIn+Nweight
        !fdin=fdmin+(iweight-1)*stepfd
        !freq=fdin/d
        cphase = cphaseVec(iweight)
        freq=fdVec(iweight)/d
        omegaj=TWO*PI*freq
        call Calculate_Weigth_Burst(Weigth_Burst,fc,freq,Nc)
        call calculateUxUz(ux,uz,z,cp,cs,d,omegaj,cphase,antisym)
        sum_ux=sum_ux+Weigth_Burst*ux
        sum_uz=sum_uz+Weigth_Burst*uz
        sum_Weight=sum_Weight+Weigth_Burst

!        print *
!        print *,'***************************************************'
!        print *, 'verif dans la boucle '
!        print *, 'cphase = ',cphase,'m/s @fd = ',fdout
!        print *, 'sum_ux =' ,  sum_ux , 'ux =' , ux
!        print *, 'sum_uz =' ,  sum_uz , 'uz =' , uz
!        print *, 'Weight_Burst = ', Weigth_Burst ,'sum_Weight =', sum_Weight
!        print *,'***************************************************'
!        print *
    enddo
    sum_ux=sum_ux/sum_Weight; ! division by the sum of weights
    sum_uz=sum_uz/sum_Weight; ! of course ! it is a weighted Sum !
    !! but in fact, it does no matter because it must be OK up to a constant

  end subroutine weighted_sum_Lamb_disp

!
! ----------------------------------------------------------------------------------------
!

 function cosh_cmplx(x)
    use constants, only: TWO,CUSTOM_REAL
    complex(kind=CUSTOM_REAL) :: cosh_cmplx,x
    cosh_cmplx = (exp(x)+exp(-x))/TWO
  end function cosh_cmplx

!
! ----------------------------------------------------------------------------------------
!

  function sinh_cmplx(x)
    use constants, only: TWO,CUSTOM_REAL
    complex(kind=CUSTOM_REAL) :: sinh_cmplx,x
    sinh_cmplx = (exp(x)-exp(-x))/TWO
  end function sinh_cmplx

!
! ----------------------------------------------------------------------------------------
!

  subroutine calculateUxUz(ux,uz,zi,cp,cs,d,omegaj,cphase,antisym)
    ! ----------------------------------------------------------------------------------------
    ! See eq.Viktorov page 70 ! todo to check page number
    ! ----------------------------------------------------------------------------------------
    use constants, only: CUSTOM_REAL,TWO,PI
    implicit none

    ! Inputs
    real(kind=CUSTOM_REAL), intent(in) :: zi,cp,cs,d,omegaj
    real(kind=CUSTOM_REAL), intent(in) :: cphase
    logical, intent(in) :: antisym

    ! Outputs
    complex(kind=CUSTOM_REAL), intent(out) :: ux,uz

    ! Local variables
    complex(kind=CUSTOM_REAL) :: s,s2,q,q2,k,k2,sqrkp,sqrks,C1,C2
    complex(kind=CUSTOM_REAL) :: jj = (0.0,1.0)

    ! note: since cosh and sinh are intrinsic functions, the complex functions defined in this file are renamed to **_cmplx
    ! functions
    complex(kind=CUSTOM_REAL), external :: cosh_cmplx,sinh_cmplx

    k2 = (omegaj/cphase)**2
    sqrkp=(omegaj/cp)**2
    sqrks=(omegaj/cs)**2
    q2 = k2 - sqrkp
    s2 = k2 - sqrks
    q = sqrt(q2)
    s = sqrt(s2)
    k = sqrt(k2)
    C1=TWO*q*s/(k2+s2)
    C2=TWO*k2/(k2+s2)

    if (.not. antisym) then ! symmetric modes
        ux = jj * k * (cosh_cmplx(q*zi)/sinh_cmplx(q*d) - C1*cosh_cmplx(s*zi)/sinh_cmplx(s*d))
        uz = -q * (sinh_cmplx(q*zi)/sinh_cmplx(q*d) - C2*sinh_cmplx(s*zi)/sinh_cmplx(s*d))
    else                    ! antisymmetric modes
        ux = jj * k * (sinh_cmplx(q*zi)/cosh_cmplx(q*d) - C1*sinh_cmplx(s*zi)/cosh_cmplx(s*d))
        uz = -q * (cosh_cmplx(q*zi)/cosh_cmplx(q*d) - C2*cosh_cmplx(s*zi)/cosh_cmplx(s*d))
    endif

!    print *
!     print *,'***************************************************'
!~     print *, 'dans la fonction calculate ux uz'
!~     print *, 'fd= ' , omegaj/TWO/PI*d
 !   print *,'cphase = ',cphase,'m/s @fd = ',omegaj/TWO/PI*d
!~     print *, 'C1= ', C1, 'C2= ', C2
!~     print *,  'q=', q ,'s= ',s,'k= ',k
!     print *, 'ux =' ,  ux , 'uz =' , uz
!    print *,'***************************************************'
!~     print *

  end subroutine calculateUxUz

!
! ----------------------------------------------------------------------------------------
!

  subroutine readCphaseAndFdInFile(antisym,order,fdin)

    use constants, only: CUSTOM_REAL
    use enforce_par, only: fdVec,cphaseVec,indexFdIn,number_of_lines_in_file

    implicit none
    ! Inputs
    character(len=1) :: order
    logical, intent(in) :: antisym
    real(kind=CUSTOM_REAL), intent(in) :: fdin

    ! Local variables
    integer :: i,error=0!,nfdout1,nfdout2
    !real(kind=CUSTOM_REAL) :: fdMin_of_file,fdMax_of_file,stepfd_of_file
    !real(kind=CUSTOM_REAL) ::  fdtemp1,fdtemp2
    character(len=1) :: mode
    character(len=1024) :: directory
    character(len=1024) :: file2read
    logical :: hasBeenFound = .false.

    directory= '/home1/bottero/specfem2d/EXAMPLES/delamination/DispersionCurves/' ! WARNING DO NOT FORGET THE / AT THE END
    ! choose the file to read
    if (.not. antisym) then
    mode='S'
    else
    mode='A'
    endif
    file2read = trim(directory) // trim(mode) // trim(order) // '.csv' ! file to read

    open(unit=7,file=file2read,action='read',status='old')

    do while(error == 0)    ! determine the number of line in the file
    read(7,*,iostat=error)
    number_of_lines_in_file=number_of_lines_in_file+1
    enddo
    rewind 7 ! to restart at the begining of the file

    allocate(fdVec(number_of_lines_in_file))
    allocate(cphaseVec(number_of_lines_in_file))

    do i = 1, number_of_lines_in_file    !read the content of the file
      read(7,*) fdVec(i), cphaseVec(i)
      if ((fdVec(i) >= fdin) .and. (.not. hasBeenFound)) then
        indexFdIn = i - 1
        hasBeenFound = .true.
      endif
    enddo

!!    fdMin_of_file=zfd(1);
!!    fdMax_of_file=zfd(number_of_lines_in_file)
!!    stepfd_of_file=(fdMax_of_file-fdMin_of_file)/(number_of_lines_in_file-1)
!!    nfdout1=floor((fdin-fdMin_of_file)/stepfd_of_file)+1
!!    nfdout2=ceiling ((fdin-fdMin_of_file)/stepfd_of_file)+1
    !nfdout1=floor(fdin/fdMin_of_file)
    !nfdout1=ceiling(fdin/fdMin_of_file)


!    if (fdin-fdMin_of_file < 0) then
!    call exit_MPI(myrank,'fdin>fdMin_of_file, change fdin product or increase Nc!')
!    endif

!!    fdtemp1=zfd(nfdout1)
!!    fdtemp2=zfd(nfdout2)

    !sometime if we use only nfdout1 we do not have the closet product
    !the test below is to fix that pb
!!    if (abs(fdtemp1-fdin) < abs(fdtemp2-fdin)) then
!!    fdout=zfd(nfdout1)
!!    cphase=zcphase(nfdout1)
!!    else
!!    fdout=zfd(nfdout2)
!!    cphase=zcphase(nfdout2)
!!    endif


!    print *
!    print *,'**************************************************'
!    print *,'subroutine readCphaseAndFdInFile '
!    print *,'file2read = ', file2read
!    print *,'number_of_lines_in_file = ', number_of_lines_in_file
!    print *,'fdMin_of_file = ', fdMin_of_file, 'fdMax_of_file =', fdMax_of_file
!    print *, 'stepfd_of_file, nfdout, nfdout2=' ,stepfd_of_file, nfdout1, nfdout2
!    print *,'stepfd =',stepfd_of_file, ' fdin = ', fdin, ' fdout =', fdout
!    print *,'cphase = ',cphase,'m/s @fd = ',fdout
!    print *,'**************************************************'

    close(7)
!    deallocate(zfd)
!    deallocate(zcphase)

  end subroutine readCphaseAndFdInFile

!
! ----------------------------------------------------------------------------------------
!

subroutine Calculate_Weigth_Burst(Weigth_Burst,f0,freq,Nc)
! ----------------------------------------------------------------------
! @ freq=f0 +/- 2* f0/Nc the DSP is null
! it is the limit of the principal lobe
! @ freq=f0 the limit of the normalized DSP=1,
! @ freq=f0 +/-f0/Nc the limit of the normalized DSP=0.5
! ----------------------------------------------------------------------
    !use specfem_par, only: myrank

    use constants, only: CUSTOM_REAL,PI,TWO
    implicit none
    ! Inputs
    real(kind=CUSTOM_REAL), intent(in) :: freq,f0
    integer,intent(in) :: Nc
    ! Outputs
    real(kind=CUSTOM_REAL),intent(out) :: Weigth_Burst

    !   Local variables
    real(kind=CUSTOM_REAL), parameter :: SMALLVAL=1e-1
    real(kind=CUSTOM_REAL) :: cste,den,M,f04,f,fbw1,fbw2
    ! SMALLVAL not TINYVAL because it does not work if val in too tiny !


    fbw1=f0-f0/Nc
    fbw2=f0+f0/Nc

    f=freq-f0
    if (abs(f) < SMALLVAL) then
    Weigth_Burst=1
    else if (abs(freq-fbw1) < SMALLVAL .or. abs(freq-fbw2) < SMALLVAL ) then
    Weigth_Burst=0.5
    else
    M=sqrt(Nc**2/f0**2)    !max @ f0
    cste=sqrt(TWO)/TWO/PI/M
    den=(Nc**2*f**2-f0**2)**2*f**2
    f04=f0**4
    Weigth_Burst=cste*sqrt(-1.0*f04*(cos(2*Nc*PI*f/f0)-1)/den)
    endif

!~     print *
!~     print *,'***************************************************'
!~     print *, 'Verif Calculate Weigth '
!~     print *, 'fd =' , freq*d,'freq= ', freq, 'freq - f0 = ', f
!~     print *, 'f0 = ',f0,'Nc= ',Nc, 'Nweight= ', Nweight
!~     print *, 'cos(2*Nc*PI*f/f0)= ', cos(2*Nc*PI*f/f0)
!~     print *, 'M= ', M,'cste= ',cste,'den= ', den
!~     print *, 'Weight_Burst =====' , Weigth_Burst!~
!~     print *,'***************************************************'
!~     print *
!~     print *, freq,Weigth_Burst ! write Weigth_Burst to plot and check


end subroutine Calculate_Weigth_Burst

!
! ----------------------------------------------------------------------------------------
!

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

! subroutine enforce_fields_acoustic(iglob,it)
!! ----------------------------------------------------------------------------------------
!! This subroutine impose the fields at a given GLL point and at a given time step
!! ----------------------------------------------------------------------------------------

!  use specfem_par, only: coord,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic,deltat,deltatover2, &
!                         deltatsquareover2
!  use enforce_par
!  use constants, only: TINYVAL,CUSTOM_REAL,TWO,PI

!  implicit none

!  ! Inputs
!  integer, intent(in) :: iglob,it

!  ! Local variables

!  ! Local variables
!  real(kind=CUSTOM_REAL) :: potential_dot_dot_acoustic_old,potential_dot_acoustic_old
!  real(kind=CUSTOM_REAL) :: factor,x,z
!  real(kind=CUSTOM_REAL) :: f0 = 2.0d0

!  x = coord(1,iglob)
!  z = coord(2,iglob)

!  factor = 1.0d0 ! * (z - 2.0d0)**2

!  !if (abs(z + 1.5d-3) < 1.5d-3) then
!    if (it == 1) then ! We initialize the variables
!      potential_acoustic(iglob) = factor*(-1.0d0)
!      potential_dot_acoustic(iglob) = factor*0.0d0
!      potential_dot_dot_acoustic(iglob) = factor*(TWO*PI*f0)**2
!    else ! We set what we want
!      potential_dot_dot_acoustic(iglob) = factor*(TWO*PI*f0)**2*cos(TWO*PI*f0*(it-1)*deltat)
!      potential_dot_dot_acoustic_old = factor*(TWO*PI*f0)**2*cos(TWO*PI*f0*(it-2)*deltat)
!      ! Do not change anything below: we compute numerically the velocity and displacement
!      potential_dot_acoustic_old = potential_dot_acoustic(iglob)
!      potential_dot_acoustic(iglob) = potential_dot_acoustic(iglob) + deltatover2*(potential_dot_dot_acoustic_old + &
!                                      potential_dot_dot_acoustic(iglob))
!      potential_acoustic(iglob) = potential_acoustic(iglob) + deltat*potential_dot_acoustic_old + &
!                                  deltatsquareover2*potential_dot_dot_acoustic_old
!    endif
!  !else
!  !  displ_elastic(:,iglob) = 0.0d0
!  !  veloc_elastic(:,iglob) = 0.0d0
!  !  accel_elastic(:,iglob) = 0.0d0
!  !endif

! end subroutine enforce_fields_acoustic

!  subroutine enforce_fields_acoustic(iglob,it)
!! ----------------------------------------------------------------------------------------
!! This subroutine impose the fields at a given GLL point and at a given time step
!! ----------------------------------------------------------------------------------------

!  use specfem_par, only: coord,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic,deltat,deltatover2, &
!                         deltatsquareover2,myrank,nLines,zmode,realMode,imagMode,modeAmplitude
!  use enforce_par
!  use interpolation
!  use constants, only: TINYVAL,CUSTOM_REAL,TWO,PI

!  implicit none

!  ! Inputs
!  integer, intent(in) :: iglob,it

!  ! Local variables
!  integer :: i,ier,idx
!  real(kind=CUSTOM_REAL) :: potential_dot_dot_acoustic_old,potential_dot_acoustic_old
!  real(kind=CUSTOM_REAL) :: factor,x,z,A,B
!  real(kind=CUSTOM_REAL) :: f0 = 2.0d0

!  x = coord(1,iglob)
!  z = coord(2,iglob)

!  factor = 1.0d0 ! * (z - 2.0d0)**2

!  ! Read files
!  if (it == 1) then  ! We initialize the variables
!    if (neverRead) then
!      nLines=count_lines("mode1realMinusSorted.txt") ! Count the number of lines in file
!      allocate(zmode(nLines),realMode(nLines),imagMode(nLines))
!      open(unit=num_file,file="mode1realMinusSorted.txt",status='old',action='read',iostat=ier)
!      if (ier /= 0) then
!        call exit_MPI(myrank,"Error reading real mode file")
!      endif
!      ! format: #depth #mode value
!      do i=1,nLines
!        read(num_file,*) zmode(i),realMode(i)
!      enddo
!      ! closes external file
!      close(num_file)
!      open(unit=num_file,file="./mode1imagMinusSorted.txt",status='old',action='read',iostat=ier)
!      if (ier /= 0) then
!        call exit_MPI(myrank,"Error reading imag mode file")
!      endif
!      ! format: #depth #mode value
!      do i=1,nLines
!        read(num_file,*) zmode(i),imagMode(i)
!      enddo
!      ! closes external file
!      close(num_file)
!      neverRead = .false.
!    endif
!    idx=searchInf(nLines,zmode,dble(z)) ! Look for index idx in sorted array zmode such as : zmode(idx) < z < zmode(idx+1)
!    ! Linear interpolation. Near to z: realMode = A*zmode+B
!    A = (realMode(idx + 1) - realMode(idx))/(zmode(idx+1) - zmode(idx))
!    B = (zmode(idx+1)*realMode(idx) - zmode(idx)*realMode(idx+1)) / (zmode(idx+1) - zmode(idx))
!    modeAmplitude(iglob) = A*z+B
!    !print *
!    !print *,"nLines:",nLines,"zmode:",zmode,"z:",dble(z)
!    !print *,idx,zmode(idx),zmode(idx+1)," so modeReal between ",realMode(idx)," and ",realMode(idx + 1)
!    !print *,"mode at",z,"m :",modeAmplitude(iglob)
!  endif

!  factor = factor*modeAmplitude(iglob)

!  !if (abs(z + 1.5d-3) < 1.5d-3) then
!    if (it == 1) then ! We initialize the variables
!      potential_acoustic(iglob) = factor*(-1.0d0)
!      potential_dot_acoustic(iglob) = factor*0.0d0
!      potential_dot_dot_acoustic(iglob) = factor*(TWO*PI*f0)**2
!    else ! We set what we want
!      potential_dot_dot_acoustic(iglob) = factor*(TWO*PI*f0)**2*cos(TWO*PI*f0*(it-1)*deltat)
!      potential_dot_dot_acoustic_old = factor*(TWO*PI*f0)**2*cos(TWO*PI*f0*(it-2)*deltat)
!      ! Do not change anything below: we compute numerically the velocity and displacement
!      potential_dot_acoustic_old = potential_dot_acoustic(iglob)
!      potential_dot_acoustic(iglob) = potential_dot_acoustic(iglob) + deltatover2*(potential_dot_dot_acoustic_old + &
!                                      potential_dot_dot_acoustic(iglob))
!      potential_acoustic(iglob) = potential_acoustic(iglob) + deltat*potential_dot_acoustic_old + &
!                                  deltatsquareover2*potential_dot_dot_acoustic_old
!    endif
!  !else
!  !  displ_elastic(:,iglob) = 0.0d0
!  !  veloc_elastic(:,iglob) = 0.0d0
!  !  accel_elastic(:,iglob) = 0.0d0
!  !endif

!  end subroutine enforce_fields_acoustic

 subroutine enforce_fields_acoustic(iglob,it)
! ----------------------------------------------------------------------------------------
! This subroutine impose the fields at given GLL points and at a given time steps
! ----------------------------------------------------------------------------------------

  use specfem_par, only: coord,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic,deltat,deltatover2, &
                         deltatsquareover2,myrank,nLines,zmode,realMode,imagMode,modeAmplitude
  use enforce_par
  use interpolation
  use constants, only: TINYVAL,CUSTOM_REAL,TWO,PI

  implicit none

  ! Inputs
  integer, intent(in) :: iglob,it

  ! Local variables
  integer :: i,ier,idx
  real(kind=CUSTOM_REAL) :: potential_dot_dot_acoustic_old,potential_dot_acoustic_old
  real(kind=CUSTOM_REAL) :: factor,x,z,A,B,tval,tval_old
  real(kind=CUSTOM_REAL) :: f0 = 2.0d0
  real(kind=CUSTOM_REAL) :: omegaj
  real(kind=CUSTOM_REAL) :: time_dependence,time_dependence_old
  integer :: Nc = 10

  factor = 1.0d0
  omegaj = TWO*PI*f0

  ! gets global node location
  x = coord(1,iglob)
  z = coord(2,iglob)

  factor = 1.0d0 ! * (z - 2.0d0)**2

  ! Read files
  if (it == 1) then  ! We initialize the variables
    if (neverRead) then
      nLines=count_lines("mode1realMinusSorted.txt") ! Count the number of lines in file
      allocate(zmode(nLines),realMode(nLines),imagMode(nLines))
      open(unit=num_file,file="mode1realMinusSorted.txt",status='old',action='read',iostat=ier)
      if (ier /= 0) then
        call exit_MPI(myrank,"Error reading real mode file")
      endif
      ! format: #depth #mode value
      do i=1,nLines
        read(num_file,*) zmode(i),realMode(i)
      enddo
      ! closes external file
      close(num_file)
      open(unit=num_file,file="./mode1imagMinusSorted.txt",status='old',action='read',iostat=ier)
      if (ier /= 0) then
        call exit_MPI(myrank,"Error reading imag mode file")
      endif
      ! format: #depth #mode value
      do i=1,nLines
        read(num_file,*) zmode(i),imagMode(i)
      enddo
      ! closes external file
      close(num_file)
      neverRead = .false.
    endif
    idx=searchInf(nLines,zmode,dble(z)) ! Look for index idx in sorted array zmode such as : zmode(idx) < z < zmode(idx+1)
    ! Linear interpolation. Near to z: realMode = A*zmode+B
    A = (realMode(idx + 1) - realMode(idx))/(zmode(idx+1) - zmode(idx))
    B = (zmode(idx+1)*realMode(idx) - zmode(idx)*realMode(idx+1)) / (zmode(idx+1) - zmode(idx))
    modeAmplitude(iglob) = A*z+B
    !print *
    !print *,"nLines:",nLines,"zmode:",zmode,"z:",dble(z)
    !print *,idx,zmode(idx),zmode(idx+1)," so modeReal between ",realMode(idx)," and ",realMode(idx + 1)
    !print *,"mode at",z,"m :",modeAmplitude(iglob)
  endif

  factor = factor*modeAmplitude(iglob)

  ! timing
  tval = (it-1)*deltat
  tval_old = (it-2)*deltat

  if (tval < TWO*Nc*PI/omegaj) then
    time_dependence = omegaj**2*cos(omegaj*tval/Nc)*sin(omegaj*tval)/Nc**2 + &
                        2*omegaj**2*sin(omegaj*tval/Nc)*cos(omegaj*tval)/Nc - &
                        (1.0d0 - cos(omegaj*tval/Nc))*omegaj**2*sin(omegaj*tval)
  else
    time_dependence = 0.0d0
  endif
  if (tval_old < TWO*Nc*PI/omegaj) then
    time_dependence_old = omegaj**2*cos(omegaj*tval_old/Nc)*sin(omegaj*tval_old)/Nc**2 + &
                        2*omegaj**2*sin(omegaj*tval_old/Nc)*cos(omegaj*tval_old)/Nc - &
                        (1.0d0 - cos(omegaj*tval_old/Nc))*omegaj**2*sin(omegaj*tval_old)
  else
    time_dependence_old = 0.0d0
  endif

  !if (abs(z + 0.0d-3) < 1.0d-3) then
    if (it == 1) then
      ! We initialize the variables
      potential_acoustic(iglob) = factor*0.0d0
      potential_dot_acoustic(iglob) = factor*0.0d0
      potential_dot_dot_acoustic(iglob) = factor*0.0d0
    else
      ! We set what we want
      potential_dot_dot_acoustic(iglob) = factor*time_dependence
      potential_dot_dot_acoustic_old = factor*time_dependence_old
      ! Do not change anything below: we compute numerically the velocity and displacement
      potential_dot_acoustic_old = potential_dot_acoustic(iglob)
      potential_dot_acoustic(iglob) = potential_dot_acoustic(iglob) + deltatover2*(potential_dot_dot_acoustic_old + &
                                      potential_dot_dot_acoustic(iglob))
      potential_acoustic(iglob) = potential_acoustic(iglob) + deltat*potential_dot_acoustic_old + &
                                  deltatsquareover2*potential_dot_dot_acoustic_old
    endif
  !else
  !  displ_elastic(:,iglob) = 0.0d0
  !  veloc_elastic(:,iglob) = 0.0d0
  !  accel_elastic(:,iglob) = 0.0d0
  !endif

  end subroutine enforce_fields_acoustic



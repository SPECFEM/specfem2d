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

module model_tomography_par
! ----------------------------------------------------------------------------------------
! Contains the variables needed to read an ASCII tomo file
! ----------------------------------------------------------------------------------------

  implicit none

  ! for external tomography:
  ! (regular spaced, xyz-block file in ascii)

  ! models dimensions
  double precision  :: END_X,END_Z

  double precision :: ORIGIN_X,ORIGIN_Z
  double precision :: SPACING_X,SPACING_Z

  ! models parameter records
  double precision, dimension(:), allocatable :: x_tomo,z_tomo
  double precision, dimension(:,:), allocatable :: vp_tomo,vs_tomo,rho_tomo
  double precision, dimension(:,:), allocatable :: qp_tomo,qs_tomo

  ! models entries
  integer :: NX,NZ
  integer :: nrecord

  ! min/max statistics
  double precision :: VP_MIN,VS_MIN,RHO_MIN,VP_MAX,VS_MAX,RHO_MAX

  ! (optional) Q values
  logical :: has_q_values

end module model_tomography_par

!
! ----------------------------------------------------------------------------------------
!

module interpolation
! ----------------------------------------------------------------------------------------
! This module contains two functions for bilinear interpolation
! (modified from http://www.shocksolution.com)
! ----------------------------------------------------------------------------------------

  contains

  ! ====================== Implementation part ===============

  integer function searchInf(length, array, value,delta)
    ! Given a a sorted array and a value, returns the index of the element that
    ! is closest to, but less than, the given value.
    ! Uses a binary search algorithm.
    ! Never returns the last value !
    !
    ! Behaviour with vec=(-5,-2,8,12,40,45,60,62,80,95) :
    ! value = 42.5 --> searchInf = 40
    ! value = -2   --> searchInf = -2
    ! value = -5   --> searchInf = -5
    ! value = -95  --> searchInf = -5
    ! value = 100  --> searchInf = 80
    ! value = -8   --> searchInf = -5
    !
    ! TODO : This function could be rewriten using a binary search algorithm
    !        for a better efficiency

    implicit none

    integer, intent(in) :: length
    double precision, dimension(length), intent(in) :: array
    double precision, intent(in) :: value
    double precision, intent(in), optional :: delta
    double precision :: d
    integer :: index_closest

    index_closest = 1

    if (present(delta)) then
      d = delta
    else
      d = 1e-9
    endif

    if (length <= 0) then
      call stop_the_code("Incorrect length in searchInf")
    else if (length == 1) then
      searchInf = array(1)
    else
      if (array(1) > array(2)) call stop_the_code("searchInf needs an increasing array")
      if (value < array(1)) then
        searchInf = 1
        return
      else if (value >= array(length)) then
        searchInf = length - 1
        return
      endif
      index_closest = 1
      do while (array(index_closest) < value)
        index_closest = index_closest+1
      enddo
      if (abs(array(index_closest)-value) > d) then
        searchInf = index_closest - 1
      else
        searchInf = index_closest
      endif

    endif
    return
  end function searchInf

  !
  !-------------------------------------------------------------------------------------------------
  !

  double precision function interpolate(x_len, x_array, y_len, y_array, f, x, y, delta)
    ! This function uses bilinear interpolation to estimate the value
    ! of a function f at point (x,y)
    ! f is assumed to be sampled on a regular grid, with the grid x values specified
    ! by x_array and the grid y values specified by y_array
    ! This function works even if (x,y) is outside the definition domain of f
    ! (case 1,2,3,4,5,6,7,8 below). It that case we just use border values
    !        y
    !        ^
    ! case 3 | case 6 | case 2
    ! ------------------------
    !        |        |
    ! case 7 | case 9 | case 8
    !        |        |
    ! ------------------------> x
    ! case 1 | case 5 | case 4
    ! Reference: http://en.wikipedia.org/wiki/Bilinear_interpolation

    use constants, only: TINYVAL

    implicit none

    integer, intent(in) :: x_len, y_len
    double precision, dimension(x_len), intent(in) :: x_array
    double precision, dimension(y_len), intent(in) :: y_array
    double precision, dimension(x_len, y_len), intent(in) :: f
    double precision, intent(in) :: x,y
    double precision, intent(in), optional :: delta

    ! Local variables
    double precision :: denom, x1, x2, y1, y2, minx, maxx, miny, maxy
    integer :: i,j
    double precision :: d

    if (present(delta)) then
      d = delta
    else
      d = 1e-9
    endif

    minx = minval(x_array)
    maxx = maxval(x_array)
    miny = minval(y_array)
    maxy = maxval(y_array)
    !        y
    !        ^
    ! case 3 | case 6 | case 2
    ! ------------------------
    !        |        |
    ! case 7 | case 9 | case 8
    !        |        |
    ! ------------------------> x
    ! case 1 | case 5 | case 4

    ! case 1
    !if (x <= minx .and. y <= miny) then
    ! considers floating point round-off errors
    if ( ((x < minx) .or. (abs(x-minx) < TINYVAL)) .and. ((y < miny) .or. (abs(y-miny) < TINYVAL)) ) then
      !print *,"case1"
      interpolate = f(1,1)

    ! case 2
    !else if (x > maxx .and. y >= maxy) then
    else if ( (x > maxx) .and. ((y > maxy) .or. (abs(y-maxy) < TINYVAL)) ) then
      !print *,"case2"
      interpolate = f(x_len,y_len)

    ! case 3
    !else if (x <= minx .and. y >= maxy) then
    else if ( ((x < minx) .or. (abs(x-minx) < TINYVAL)) .and. ((y > maxy) .or. (abs(y-maxy) < TINYVAL)) ) then
      !print *,"case3"
      interpolate = f(1,y_len)

    ! case 4
    !else if (x >= maxx .and. y <= miny) then
    else if ( ((x > maxx) .or. (abs(x-maxx) < TINYVAL)) .and. ((y < miny) .or. (abs(y-miny) < TINYVAL)) ) then
      !print *,"case4"
      interpolate = f(x_len,1)

    ! case 5
    !else if (x >= minx .and. x <= maxx .and. y <= miny) then
    else if ( ((x > minx) .or. (abs(x-minx) < TINYVAL)) .and. &
              ((x < maxx) .or. (abs(x-maxx) < TINYVAL)) .and. &
              ((y < miny) .or. (abs(y-miny) < TINYVAL)) ) then
      !print *,"case5"
      y1 = y_array(1)
      y2 = y_array(2)
      ! y = y1
      i = searchInf(x_len, x_array, x) ! Test also if the array is increasing
      j = 1
      if (i < x_len) then
        x1 = x_array(i)
        x2 = x_array(i+1)
      else
        x1 = x_array(i-1)
        x2 = x_array(i)
      endif
      denom = (x2 - x1)*(y2 - y1)
      interpolate = (f(i,j)*(x2-x)*(y2-y1) + f(i+1,j)*(x-x1)*(y2-y1))/denom

    ! case 6
    !else if (x >= minx .and. x <= maxx .and. y >= maxy) then
    else if ( ((x > minx) .or. (abs(x-minx) < TINYVAL)) .and. &
              ((x < maxx) .or. (abs(x-maxx) < TINYVAL)) .and. &
              ((y > maxy) .or. (abs(y-maxy) < TINYVAL)) ) then
      !print *,"case6"
      y1 = y_array(y_len-1)
      y2 = y_array(y_len)
      ! y = y2
      i = searchInf(x_len, x_array, x) ! Test also if the array is increasing
      j = y_len - 1
      if (i < x_len) then
        x1 = x_array(i)
        x2 = x_array(i+1)
      else
        x1 = x_array(i-1)
        x2 = x_array(i)
      endif
      denom = (x2 - x1)*(y2 - y1)
      interpolate = (f(i,j+1)*(x2-x)*(y2-y1) + f(i+1, j+1)*(x-x1)*(y2-y1))/denom

    ! case 7
    !else if (x <= minx .and. y >= miny .and. y <= maxy) then
    else if ( ((x < minx) .or. (abs(x-minx) < TINYVAL)) .and. &
              ((y > miny) .or. (abs(y-miny) < TINYVAL)) .and. &
              ((y < maxy) .or. (abs(y-maxy) < TINYVAL)) ) then
      !print *,"case7"
      x1 = x_array(1)
      x2 = x_array(2)
      ! x = x1
      i = 1
      j = searchInf(y_len, y_array, y) ! Test also if the array is increasing
      if (j < y_len) then
        y1 = y_array(j)
        y2 = y_array(j+1)
      else
        y1 = y_array(j-1)
        y2 = y_array(j)
      endif
      denom = (x2 - x1)*(y2 - y1)
      interpolate = (f(i,j)*(x2-x1)*(y2-y) + f(i,j+1)*(x2-x1)*(y-y1))/denom

    ! case 8
    !else if (x >= maxx .and. y >= miny .and. y <= maxy) then
    else if ( ((x > maxx) .or. (abs(x-maxx) < TINYVAL)) .and. &
              ((y > miny) .or. (abs(y-miny) < TINYVAL)) .and. &
              ((y < maxy) .or. (abs(y-maxy) < TINYVAL)) ) then
      !print *,"case8"
      x1 = x_array(x_len-1)
      x2 = x_array(x_len)
      ! x = x2
      i = x_len - 1
      j = searchInf(y_len, y_array, y) ! Test also if the array is increasing
      if (j < y_len) then
        y1 = y_array(j)
        y2 = y_array(j+1)
      else
        y1 = y_array(j-1)
        y2 = y_array(j)
      endif
      denom = (x2 - x1)*(y2 - y1)
      interpolate = (f(i+1,j)*(x2-x1)*(y2-y) + f(i+1,j+1)*(x2-x1)*(y-y1))/denom

    ! case 9
    else
      ! The point is exactly on the area defined
      !print *,"case9"
      i = searchInf(x_len, x_array, x) ! Test also if the array is increasing
      j = searchInf(y_len, y_array, y) ! Test also if the array is increasing

      if (i < x_len) then
        x1 = x_array(i)
        x2 = x_array(i+1)
      else
        x1 = x_array(i-1)
        x2 = x_array(i)
      endif
      if (j < y_len) then
        y1 = y_array(j)
        y2 = y_array(j+1)
      else
        y1 = y_array(j-1)
        y2 = y_array(j)
      endif
      denom = (x2 - x1)*(y2 - y1)
      interpolate = (f(i,j)*(x2-x)*(y2-y) + f(i+1,j)*(x-x1)*(y2-y) + &
                    f(i,j+1)*(x2-x)*(y-y1) + f(i+1, j+1)*(x-x1)*(y-y1))/denom
    endif

 end function interpolate

end module interpolation

!
! ----------------------------------------------------------------------------------------
!

  subroutine define_external_model_from_tomo_file(rhoext,vpext,vsext, &
                                                  QKappa_attenuationext,Qmu_attenuationext, &
                                                  c11ext,c12ext,c13ext,c15ext,c22ext,c23ext,c25ext,c33ext,c35ext,c55ext)

! ----------------------------------------------------------------------------------------
! Read a tomo file and loop over all GLL points to set the values of vp,vs and rho
! ----------------------------------------------------------------------------------------

  use specfem_par, only: tomo_material,coord,nspec,ibool,kmato, &
                         Qkappa_attenuationcoef,Qmu_attenuationcoef,anisotropycoef, &
                         poroelastcoef,density

  use specfem_par, only: myrank,TOMOGRAPHY_FILE

  use model_tomography_par
  use interpolation

  use constants, only: NGLLX,NGLLZ,TINYVAL,IMAIN,CUSTOM_REAL,ATTENUATION_COMP_MAXIMUM

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: rhoext,vpext,vsext
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: Qkappa_attenuationext,Qmu_attenuationext
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec), intent(out) :: c11ext,c12ext,c13ext,c15ext,c22ext,c23ext,c25ext, &
                                                                       c33ext,c35ext,c55ext

  ! local parameters
  integer :: i,j,ispec,iglob
  double precision :: xmesh,zmesh
  double precision :: rho_final
  double precision :: vp_final,vs_final
  double precision :: qp_final,qs_final
  double precision :: L_val,qmu_atten,qkappa_atten
  ! stats
  integer :: npoint_tomo,npoint_internal

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  tomographic model file: ',trim(TOMOGRAPHY_FILE)
    call flush_IMAIN()
  endif

  ! Read external tomo file TOMOGRAPHY_FILE
  call read_tomo_file()

  ! default no attenuation
  QKappa_attenuationext(:,:,:) = ATTENUATION_COMP_MAXIMUM
  Qmu_attenuationext(:,:,:) = ATTENUATION_COMP_MAXIMUM

  ! default no anisotropy
  c11ext(:,:,:) = 0.d0
  c13ext(:,:,:) = 0.d0
  c15ext(:,:,:) = 0.d0
  c33ext(:,:,:) = 0.d0
  c35ext(:,:,:) = 0.d0
  c55ext(:,:,:) = 0.d0
  c12ext(:,:,:) = 0.d0
  c23ext(:,:,:) = 0.d0
  c25ext(:,:,:) = 0.d0
  c22ext(:,:,:) = 0.d0 ! for AXISYM only

  ! loop on all the elements of the mesh, and inside each element loop on all the GLL points
  npoint_tomo = 0
  npoint_internal = 0
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        if (kmato(ispec) == tomo_material) then
          ! If the material has been set to < 0 on the Par_file
          xmesh = coord(1,iglob)
          zmesh = coord(2,iglob)

          rho_final = interpolate(NX, x_tomo, NZ, z_tomo, rho_tomo, xmesh, zmesh, TINYVAL)
          vp_final = interpolate(NX, x_tomo, NZ, z_tomo, vp_tomo, xmesh, zmesh, TINYVAL)
          vs_final = interpolate(NX, x_tomo, NZ, z_tomo, vs_tomo, xmesh, zmesh, TINYVAL)

          rhoext(i,j,ispec) = rho_final
          vpext(i,j,ispec) = vp_final
          vsext(i,j,ispec) = vs_final

          if (has_q_values) then
            qp_final = interpolate(NX, x_tomo, NZ, z_tomo, qp_tomo, xmesh, zmesh, TINYVAL)
            qs_final = interpolate(NX, x_tomo, NZ, z_tomo, qs_tomo, xmesh, zmesh, TINYVAL)

            ! attenuation zero (means negligible attenuation)
            if (qs_final <= TINYVAL) qs_final = ATTENUATION_COMP_MAXIMUM
            if (qp_final <= TINYVAL) qp_final = ATTENUATION_COMP_MAXIMUM

            ! Anderson & Hart (1978) conversion between (Qp,Qs) and (Qkappa,Qmu)
            ! factor L
            L_val = 4.d0/3.d0 * (vs_final/vp_final)**2

            ! shear attenuation
            qmu_atten = qs_final

            ! converts to bulk attenuation
            if (abs(qs_final - L_val * qp_final) <= TINYVAL) then
              ! negligible bulk attenuation
              qkappa_atten = ATTENUATION_COMP_MAXIMUM
            else
              qkappa_atten = (1.d0 - L_val) * qp_final * qs_final / (qs_final - L_val * qp_final)
            endif

            ! attenuation zero (means negligible attenuation)
            if (qmu_atten <= TINYVAL) qmu_atten = ATTENUATION_COMP_MAXIMUM
            if (qkappa_atten <= TINYVAL) qkappa_atten = ATTENUATION_COMP_MAXIMUM

            ! limits Q values
            if (qmu_atten < 1.0d0) qmu_atten = 1.0d0
            if (qmu_atten > ATTENUATION_COMP_MAXIMUM) qmu_atten = ATTENUATION_COMP_MAXIMUM
            if (qkappa_atten < 1.0d0) qkappa_atten = 1.0d0
            if (qkappa_atten > ATTENUATION_COMP_MAXIMUM) qkappa_atten = ATTENUATION_COMP_MAXIMUM

            Qkappa_attenuationext(i,j,ispec) = qkappa_atten
            Qmu_attenuationext(i,j,ispec) = qmu_atten
          endif

          !! ABAB : The 3 following lines are important, otherwise PMLs won't work. TODO check that
          !! (we assign these values several times: indeed for each kmato(ispec) it can exist a lot of rhoext(i,j,ispec) )
          density(1,kmato(ispec)) = rhoext(i,j,ispec)
          poroelastcoef(3,1,kmato(ispec)) = rhoext(i,j,ispec) * vpext(i,j,ispec) * vpext(i,j,ispec)
          poroelastcoef(2,1,kmato(ispec)) =  rhoext(i,j,ispec) * vsext(i,j,ispec) * vsext(i,j,ispec)

          !! ABAB : I do the same with anisotropy and attenuation even if I don't use them (for the future) :
          anisotropycoef(1,kmato(ispec)) = c11ext(i,j,ispec)
          anisotropycoef(2,kmato(ispec)) = c13ext(i,j,ispec)
          anisotropycoef(3,kmato(ispec)) = c15ext(i,j,ispec)
          anisotropycoef(4,kmato(ispec)) = c33ext(i,j,ispec)
          anisotropycoef(5,kmato(ispec)) = c35ext(i,j,ispec)
          anisotropycoef(6,kmato(ispec)) = c55ext(i,j,ispec)
          anisotropycoef(7,kmato(ispec)) = c12ext(i,j,ispec)
          anisotropycoef(8,kmato(ispec)) = c23ext(i,j,ispec)
          anisotropycoef(9,kmato(ispec)) = c25ext(i,j,ispec)
          anisotropycoef(10,kmato(ispec)) = c22ext(i,j,ispec) ! for AXISYM

          QKappa_attenuationcoef(kmato(ispec)) = QKappa_attenuationext(i,j,ispec)
          Qmu_attenuationcoef(kmato(ispec)) = Qmu_attenuationext(i,j,ispec)

          ! counting GLL points with tomography model values assigned
          npoint_tomo = npoint_tomo + 1

        else
          ! Internal model
          rhoext(i,j,ispec) = density(1,kmato(ispec))
          vpext(i,j,ispec) = sqrt(poroelastcoef(3,1,kmato(ispec))/rhoext(i,j,ispec))
          vsext(i,j,ispec) = sqrt(poroelastcoef(2,1,kmato(ispec))/rhoext(i,j,ispec))

          QKappa_attenuationext(i,j,ispec) = QKappa_attenuationcoef(kmato(ispec))
          Qmu_attenuationext(i,j,ispec) = Qmu_attenuationcoef(kmato(ispec))

          c11ext(i,j,ispec) = anisotropycoef(1,kmato(ispec))
          c13ext(i,j,ispec) = anisotropycoef(2,kmato(ispec))
          c15ext(i,j,ispec) = anisotropycoef(3,kmato(ispec))
          c33ext(i,j,ispec) = anisotropycoef(4,kmato(ispec))
          c35ext(i,j,ispec) = anisotropycoef(5,kmato(ispec))
          c55ext(i,j,ispec) = anisotropycoef(6,kmato(ispec))
          c12ext(i,j,ispec) = anisotropycoef(7,kmato(ispec))
          c23ext(i,j,ispec) = anisotropycoef(8,kmato(ispec))
          c25ext(i,j,ispec) = anisotropycoef(9,kmato(ispec))
          c22ext(i,j,ispec) = anisotropycoef(10,kmato(ispec)) ! for AXISYM

          ! counting GLL points with internal model values assigned
          npoint_tomo = npoint_tomo + 1
        endif
      enddo
    enddo
  enddo

  call synchronize_all()

  ! collect totals
  iglob = npoint_tomo
  call max_all_i(iglob,npoint_tomo)
  iglob = npoint_internal
  call max_all_i(iglob,npoint_internal)

  ! user output
  if (myrank == 0) then
    if (npoint_tomo > 0) &
      write(IMAIN,*) '  number of GLL points with tomographic model values: ',npoint_tomo
    if (npoint_internal > 0) &
      write(IMAIN,*) '  number of GLL points with internal model values   : ',npoint_internal
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine define_external_model_from_tomo_file

!
!-------------------------------------------------------------------------------------------
!

  subroutine read_tomo_file()

! ----------------------------------------------------------------------------------------
! This subroutine reads the external ASCII tomo file TOMOGRAPHY_FILE (path to which is
! given in the Par_file).
! This file format is not very clever however it is the one used in specfem3D hence
! we chose to implement it here as well
! The external tomographic model is represented by a grid of points with assigned material
! properties and homogeneous resolution along each spatial direction x and z. The xyz file
! TOMOGRAPHY_FILE that describe the tomography should be located in the TOMOGRAPHY_PATH
! directory, set in the Par_file. The format of the file, as read from
! define_external_model_from_xyz_file.f90 looks like :
!
! ORIGIN_X ORIGIN_Z END_X END_Z
! SPACING_X SPACING_Z
! NX NZ
! VP_MIN VP_MAX VS_MIN VS_MAX RHO_MIN RHO_MAX
! x(1) z(1) vp vs rho
! x(2) z(1) vp vs rho
! ...
! x(NX) z(1) vp vs rho
! x(1) z(2) vp vs rho
! x(2) z(2) vp vs rho
! ...
! x(NX) z(2) vp vs rho
! x(1) z(3) vp vs rho
! ...
! ...
! x(NX) z(NZ) vp vs rho
!
!
! in case Q-values are also provided, the data lines here would have a format like:
! x(1) z(1) vp vs rho Qp Qs
! ..

! Where :
! _x and z must be increasing
! _ORIGIN_X, END_X are, respectively, the coordinates of the initial and final tomographic
!  grid points along the x direction (in meters)
! _ORIGIN_Z, END_Z are, respectively, the coordinates of the initial and final tomographic
!  grid points along the z direction (in meters)
! _SPACING_X, SPACING_Z are the spacing between the tomographic grid points along the x
!  and z directions, respectively (in meters)
! _NX, NZ are the number of grid points along the spatial directions x and z,
!  respectively; NX is given by [(END_X - ORIGIN_X)/SPACING_X]+1; NZ is the same as NX, but
!  for z direction.
! _VP_MIN, VP_MAX, VS_MIN, VS_MAX, RHO_MIN, RHO_MAX are the minimum and maximum values of
!  the wave speed vp and vs (in m.s-1) and of the density rho (in kg.m-3); these values
!  could be the actual limits of the tomographic parameters in the grid or the minimum
!  and maximum values to which we force the cut of velocity and density in the model.
! _After these first four lines, in the file file_name the tomographic grid points are
!  listed with the corresponding values of vp, vs and rho, scanning the grid along the x
!  coordinate (from ORIGIN_X to END_X with step of SPACING_X) for each given z (from ORIGIN_Z
!  to END_Z, with step of SPACING_Z).
! ----------------------------------------------------------------------------------------

  use specfem_par, only: myrank,TOMOGRAPHY_FILE

  use model_tomography_par
  use constants, only: IIN,IMAIN,HUGEVAL,MAX_STRING_LEN

  implicit none

  ! local parameters
  integer :: ier,irecord,i,j
  character(len=MAX_STRING_LEN) :: string_read

  double precision, dimension(:), allocatable :: x_tomography,z_tomography,vp_tomography,vs_tomography,rho_tomography
  double precision, dimension(:), allocatable :: qp_tomography,qs_tomography
  double precision :: read_rho_min,read_rho_max
  double precision :: read_vp_min,read_vp_max
  double precision :: read_vs_min,read_vs_max
  double precision :: read_qp_min,read_qp_max
  double precision :: read_qs_min,read_qs_max
  double precision :: tmp_dp
  integer :: ntokens

  ! opens file for reading
  open(unit=IIN,file=trim(TOMOGRAPHY_FILE),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error: could not open tomography file: ',trim(TOMOGRAPHY_FILE)
    print *,'Please check your settings in Par_file ...'
    call exit_MPI(myrank,'Error reading tomography file')
  endif

  ! --------------------------------------------------------------------------------------
  ! header infos
  ! --------------------------------------------------------------------------------------
  ! reads in model dimensions
  ! format: #origin_x #origin_z #end_x #end_z
  call tomo_read_next_line(IIN,string_read)
  read(string_read,*) ORIGIN_X, ORIGIN_Z, END_X, END_Z

  ! --------------------------------------------------------------------------------------
  ! model increments
  ! format: #dx #dz
  ! --------------------------------------------------------------------------------------
  call tomo_read_next_line(IIN,string_read)
  read(string_read,*) SPACING_X, SPACING_Z

  ! --------------------------------------------------------------------------------------
  ! reads in models entries
  ! format: #nx #nz
  ! --------------------------------------------------------------------------------------
  call tomo_read_next_line(IIN,string_read)
  read(string_read,*) NX,NZ

  ! --------------------------------------------------------------------------------------
  ! reads in models min/max statistics
  ! format: #vp_min #vp_max #vs_min #vs_max #density_min #density_max
  ! --------------------------------------------------------------------------------------
  call tomo_read_next_line(IIN,string_read)
  read(string_read,*)  VP_MIN,VP_MAX,VS_MIN,VS_MAX,RHO_MIN,RHO_MAX

  ! Determines total maximum number of element records
  nrecord = int(NX*NZ)

  ! allocate models parameter records
  allocate(x_tomography(nrecord),z_tomography(nrecord),vp_tomography(nrecord),vs_tomography(nrecord), &
           rho_tomography(nrecord),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate tomo arrays')
  x_tomography(:) = 0.d0; z_tomography(:) = 0.d0
  vp_tomography(:) = 0.d0; vs_tomography(:) = 0.d0; rho_tomography(:) = 0.d0

  allocate(x_tomo(NX),z_tomo(NZ),vp_tomo(NX,NZ),vs_tomo(NX,NZ),rho_tomo(NX,NZ),stat=ier)
  if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate tomo arrays')
  x_tomo(:) = 0.d0; z_tomo(:) = 0.d0
  vp_tomo(:,:) = 0.d0; vs_tomo(:,:) = 0.d0; rho_tomo(:,:) = 0.d0

  ! Checks the number of records for points definition while storing them
  irecord = 0

  ! stats
  read_vp_min = + HUGEVAL
  read_vp_max = - HUGEVAL
  read_vs_min = + HUGEVAL
  read_vs_max = - HUGEVAL
  read_rho_min = + HUGEVAL
  read_rho_max = - HUGEVAL
  read_qp_min = + HUGEVAL
  read_qp_max = - HUGEVAL
  read_qs_min = + HUGEVAL
  read_qs_max = - HUGEVAL

  do while (irecord < nrecord)
    call tomo_read_next_line(IIN,string_read)

    if (irecord == 0) then
      ! checks number of entries of first data line
      call tomo_get_number_of_tokens(string_read,ntokens)

      !debug
      !print *,'tomography file: number of tokens on first data line: ',ntokens,' line: ***'//trim(string_read)//'***'

      ! data line formats:
      ! #x(1) #z(1) #vp #vs #rho
      ! or
      ! #x(1) #z(1) #vp #vs #rho #Qp #Qs
      if (ntokens /= 5 .and. ntokens /= 7) then
        print *,'Error reading tomography file, data line has wrong number of entries: ',trim(string_read)
        stop 'Error reading tomography file'
      endif

      ! determines data format
      if (ntokens == 7) then
        has_q_values = .true.
      else
        has_q_values = .false.
      endif

      if (has_q_values) then
        ! allocate models parameter records
        allocate(qp_tomography(nrecord),qs_tomography(nrecord),stat=ier)
        if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate tomo Q arrays')
        qp_tomography(:) = 0.d0; qs_tomography(:) = 0.d0

        allocate(qp_tomo(NX,NZ),qs_tomo(NX,NZ),stat=ier)
        if (ier /= 0) call exit_MPI(myrank,'not enough memory to allocate tomo arrays')
        qp_tomo(:,:) = 0.d0; qs_tomo(:,:) = 0.d0
      endif
    endif

    ! reads in tomo values
    if (has_q_values) then
      ! #x(1) #z(1) #vp #vs #rho #Qp #Qs
      read(string_read,*) x_tomography(irecord+1),z_tomography(irecord+1), &
                          vp_tomography(irecord+1),vs_tomography(irecord+1),rho_tomography(irecord+1), &
                          qp_tomography(irecord+1),qs_tomography(irecord+1)
    else
      ! #x(1) #z(1) #vp #vs #rho
      read(string_read,*) x_tomography(irecord+1),z_tomography(irecord+1), &
                          vp_tomography(irecord+1),vs_tomography(irecord+1),rho_tomography(irecord+1)
    endif

    if (irecord < NX) x_tomo(irecord+1) = x_tomography(irecord+1)

    ! stats
    if (vp_tomography(irecord+1) > read_vp_max) read_vp_max = vp_tomography(irecord+1)
    if (vp_tomography(irecord+1) < read_vp_min) read_vp_min = vp_tomography(irecord+1)
    if (vs_tomography(irecord+1) > read_vs_max) read_vs_max = vs_tomography(irecord+1)
    if (vs_tomography(irecord+1) < read_vs_min) read_vs_min = vs_tomography(irecord+1)
    if (rho_tomography(irecord+1) > read_rho_max) read_rho_max = rho_tomography(irecord+1)
    if (rho_tomography(irecord+1) < read_rho_min) read_rho_min = rho_tomography(irecord+1)

    if (has_q_values) then
      if (qp_tomography(irecord+1) > read_qp_max) read_qp_max = qp_tomography(irecord+1)
      if (qp_tomography(irecord+1) < read_qp_min) read_qp_min = qp_tomography(irecord+1)
      if (qs_tomography(irecord+1) > read_qs_max) read_qs_max = qs_tomography(irecord+1)
      if (qs_tomography(irecord+1) < read_qs_min) read_qs_min = qs_tomography(irecord+1)
    endif

    ! counter
    irecord = irecord + 1
  enddo

  z_tomo = z_tomography(::NX)

  do i = 1,NX
    do j = 1,NZ
      vp_tomo(i,j) = vp_tomography(NX*(j-1)+i)
      vs_tomo(i,j) = vs_tomography(NX*(j-1)+i)
      rho_tomo(i,j) = rho_tomography(NX*(j-1)+i)
    enddo
  enddo

  if (has_q_values) then
    do i = 1,NX
      do j = 1,NZ
        qp_tomo(i,j) = qp_tomography(NX*(j-1)+i)
        qs_tomo(i,j) = qs_tomography(NX*(j-1)+i)
      enddo
    enddo
  endif

  call synchronize_all()

  if (irecord /= nrecord .and. myrank == 0) then
     print *, 'Error: ',trim(TOMOGRAPHY_FILE),' has invalid number of records'
     print *, '     number of grid points specified (= NX*NZ)   :',nrecord
     print *, '     number of file lines for grid points        :',irecord
     call stop_the_code('Error in tomography data file for the grid points definition')
  endif

  ! collects stats
  tmp_dp = read_vp_min
  call min_all_dp(tmp_dp,read_vp_min)
  tmp_dp = read_vp_max
  call max_all_dp(tmp_dp,read_vp_max)

  tmp_dp = read_vs_min
  call min_all_dp(tmp_dp,read_vs_min)
  tmp_dp = read_vs_max
  call max_all_dp(tmp_dp,read_vs_max)

  tmp_dp = read_rho_min
  call min_all_dp(tmp_dp,read_rho_min)
  tmp_dp = read_rho_max
  call max_all_dp(tmp_dp,read_rho_max)

  if (has_q_values) then
    ! Q-model values
    tmp_dp = read_qp_min
    call min_all_dp(tmp_dp,read_qp_min)
    tmp_dp = read_qp_max
    call max_all_dp(tmp_dp,read_qp_max)

    tmp_dp = read_qs_min
    call min_all_dp(tmp_dp,read_qs_min)
    tmp_dp = read_qs_max
    call max_all_dp(tmp_dp,read_qs_max)
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '     Number of grid points = NX*NZ:',nrecord
    write(IMAIN,*)
    write(IMAIN,*) '     read model vp : min/max = ',sngl(read_vp_min),' / ',sngl(read_vp_max)
    write(IMAIN,*) '                vs : min/max = ',sngl(read_vs_min),' / ',sngl(read_vs_max)
    write(IMAIN,*) '                rho: min/max = ',sngl(read_rho_min),' / ',sngl(read_rho_max)
    write(IMAIN,*)
    if (has_q_values) then
      write(IMAIN,*) '                Qp : min/max = ',sngl(read_qp_min),' / ',sngl(read_qp_max)
      write(IMAIN,*) '                Qs : min/max = ',sngl(read_qs_min),' / ',sngl(read_qs_max)
      write(IMAIN,*)
    endif
    call flush_IMAIN()
  endif

  ! closes file
  close(IIN)
  deallocate(x_tomography,z_tomography,vp_tomography,vs_tomography,rho_tomography)
  if (has_q_values) then
    deallocate(qp_tomography,qs_tomography)
  endif

  end subroutine read_tomo_file

!
! ----------------------------------------------------------------------------------------
!

  subroutine tomo_read_next_line(unit_in,string_read)

! Store next line of file unit_in in string_read

  use constants, only: MAX_STRING_LEN
  implicit none

  integer :: unit_in
  character(len=MAX_STRING_LEN) :: string_read

  integer :: ier

  do
     read(unit=unit_in,fmt="(a)",iostat=ier) string_read
     if (ier /= 0) call stop_the_code('error while reading tomography file')

     ! suppress leading white spaces, if any
     string_read = adjustl(string_read)

     ! suppress trailing carriage return (ASCII code 13) if any (e.g. if input text file coming from Windows/DOS)
     if (index(string_read,achar(13)) > 0) string_read = string_read(1:index(string_read,achar(13))-1)

     ! reads next line if empty
     if (len_trim(string_read) == 0) cycle

     ! exit loop when we find the first line that is not a comment or a white line
     if (string_read(1:1) /= '#') exit
  enddo

  ! suppress trailing white spaces, if any
  string_read = string_read(1:len_trim(string_read))

  ! suppress trailing comments, if any
  if (index(string_read,'#') > 0) string_read = string_read(1:index(string_read,'#')-1)

  ! suppress leading and trailing white spaces again, if any, after having suppressed the leading junk
  string_read = adjustl(string_read)
  string_read = string_read(1:len_trim(string_read))

  end subroutine tomo_read_next_line

!
!-------------------------------------------------------------------------------------------------
!

  subroutine tomo_get_number_of_tokens(string_read,ntokens)

  use constants, only: MAX_STRING_LEN
  implicit none

  character(len=MAX_STRING_LEN),intent(in) :: string_read
  integer,intent(out) :: ntokens

  ! local parameters
  integer :: i
  logical :: previous_is_delim
  character(len=1), parameter :: delim_space = ' '
  character(len=1), parameter :: delim_tab = achar(9) ! tab delimiter

  ! initializes
  ntokens = 0

  ! checks if anything to do
  if (len_trim(string_read) == 0) return

  ! counts tokens
  ntokens = 1
  previous_is_delim = .true.
  do i = 1, len_trim(string_read)
    ! finds next delimiter (space or tabular)
    if (string_read(i:i) == delim_space .or. string_read(i:i) == delim_tab) then
      if (.not. previous_is_delim) then
        ntokens = ntokens + 1
        previous_is_delim = .true.
      endif
    else
      if (previous_is_delim) previous_is_delim = .false.
    endif
  enddo

  end subroutine tomo_get_number_of_tokens


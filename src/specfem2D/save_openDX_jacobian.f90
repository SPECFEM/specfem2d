
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


  subroutine save_openDX_jacobian(nspec,npgeo,ngnod,knods,coorg,xigll,zigll, &
                                  AXISYM,is_on_the_axis,xiglj)

  implicit none
  include "constants.h"

  logical :: AXISYM

  integer :: nspec,npgeo,ngnod
  double precision, dimension(NDIM,npgeo) :: coorg
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

  integer, dimension(ngnod,nspec) :: knods

  double precision, dimension(NGLJ) :: xiglj
  logical, dimension(nspec) :: is_on_the_axis

  ! local parameters
  integer, dimension(:), allocatable :: ibool_OpenDX
  logical, dimension(:), allocatable :: mask_point
  double precision :: xelm,zelm
  double precision :: xi,gamma,x,z
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl

  integer :: ia,nnum,ipoint_number,total_of_negative_elements
  integer :: ispec,i,j
  logical :: found_a_problem_in_this_element

  ! create an OpenDX file to visualize this element
  open(unit=11,file='DX_all_elements_with_negative_jacobian_in_red.dx',status='unknown')

  ! output all the points (i.e. all the control points of the mesh)
  ! the mesh is flat therefore the third coordinate is zero
  write(11,*) 'object 1 class array type float rank 1 shape 3 items ',npgeo,' data follows'
  ipoint_number = 0
  allocate(mask_point(npgeo))
  allocate(ibool_OpenDX(npgeo))
  mask_point(:) = .false.
  do ispec = 1,nspec
    do ia=1,ngnod
      nnum = knods(ia,ispec)
      xelm = coorg(1,nnum)
      zelm = coorg(2,nnum)
      if(.not. mask_point(knods(ia,ispec))) then
        mask_point(knods(ia,ispec)) = .true.
        ibool_OpenDX(knods(ia,ispec)) = ipoint_number
        write(11,*) xelm,zelm,' 0'
        ipoint_number = ipoint_number + 1
      endif
    enddo
  enddo
  deallocate(mask_point)

  ! output all the elements of the mesh (use their four corners only in OpenDX
  write(11,*) 'object 2 class array type int rank 1 shape 4 items ',nspec,' data follows'
  ! point order in OpenDX is 1,4,2,3 *not* 1,2,3,4 as in AVS
  do ispec = 1,nspec
    write(11,*) ibool_OpenDX(knods(1,ispec)),ibool_OpenDX(knods(4,ispec)), &
                ibool_OpenDX(knods(2,ispec)),ibool_OpenDX(knods(3,ispec))
  enddo
  deallocate(ibool_OpenDX)

  ! output element data
  write(11,*) 'attribute "element type" string "quads"'
  write(11,*) 'attribute "ref" string "positions"'
  write(11,*) 'object 3 class array type float rank 0 items ',nspec,' data follows'

  ! output all the element data (value = 1 if positive Jacobian, = 2 if negative Jacobian)
  total_of_negative_elements = 0
  do ispec = 1,nspec

    ! check if this element has a negative Jacobian at any of its points
    found_a_problem_in_this_element = .false.
    do j = 1,NGLLZ
      do i = 1,NGLLX
        if(AXISYM) then
          if (is_on_the_axis(ispec)) then
            xi = xiglj(i)
          else
            xi = xigll(i)
          endif
        else
            xi = xigll(i)
        endif
        gamma = zigll(j)

        call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                        jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo, &
                        .false.)

        if(jacobianl <= ZERO) found_a_problem_in_this_element = .true.
      enddo
    enddo

    ! output data value
    if(found_a_problem_in_this_element) then
      write(11,*) '2'
      print *,'element ',ispec,' has a negative Jacobian'
      total_of_negative_elements = total_of_negative_elements + 1
    else
      write(11,*) '1'
    endif

  enddo

  ! define OpenDX field
  write(11,*) 'attribute "dep" string "connections"'
  write(11,*) 'object "irregular positions irregular connections" class field'
  write(11,*) 'component "positions" value 1'
  write(11,*) 'component "connections" value 2'
  write(11,*) 'component "data" value 3'
  write(11,*) 'end'

  ! close OpenDX file
  close(11)

  print *
  print *,total_of_negative_elements,' elements have a negative Jacobian, out of ',nspec
  print *,'i.e., ',sngl(100.d0 * dble(total_of_negative_elements)/dble(nspec)),'%'
  print *

  end subroutine save_openDX_jacobian

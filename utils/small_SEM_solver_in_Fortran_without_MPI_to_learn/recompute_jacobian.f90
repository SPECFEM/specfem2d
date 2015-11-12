
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

! Recompute 2D jacobian at a given point in a 4-node or 9-node element
! Compute also the global coordinates of the point defined by: (xi,gamma,ispec)

  subroutine recompute_jacobian(xi,gamma,x,z,xix,xiz,gammax,gammaz,jacobian,coorg,knods,ispec,ngnod,nspec,npgeo,NDIM)

  implicit none

  integer, intent(in) :: ispec,ngnod,nspec,npgeo,NDIM
  double precision, intent(out) :: x,z,xix,xiz,gammax,gammaz,jacobian
  double precision, intent(in) :: xi,gamma

  integer, intent(in) :: knods(ngnod,nspec)
  double precision, intent(in) :: coorg(NDIM,npgeo)

! 2D shape functions and their derivatives at receiver
  double precision shape2D(ngnod)
  double precision dershape2D(NDIM,ngnod)

  double precision xxi,zxi,xgamma,zgamma,xelm,zelm

  integer ia,nnum

! recompute jacobian for any (xi,gamma) point, not necessarily a GLL point

! create the 2D shape functions and then the Jacobian
  call define_shape_functions(shape2D,dershape2D,xi,gamma,ngnod,NDIM)

! compute coordinates and jacobian matrix
  x = 0.d0
  z = 0.d0

  xxi = 0.d0
  zxi = 0.d0
  xgamma = 0.d0
  zgamma = 0.d0

  do ia=1,ngnod

    nnum = knods(ia,ispec)

    xelm = coorg(1,nnum)
    zelm = coorg(2,nnum)

    x = x + shape2D(ia)*xelm
    z = z + shape2D(ia)*zelm

    xxi = xxi + dershape2D(1,ia)*xelm
    zxi = zxi + dershape2D(1,ia)*zelm
    xgamma = xgamma + dershape2D(2,ia)*xelm
    zgamma = zgamma + dershape2D(2,ia)*zelm

  enddo

  jacobian = xxi*zgamma - xgamma*zxi

! the Jacobian is negative, so far this means that there is an error in the mesh
! therefore print the coordinates of the mesh points of this element
  if(jacobian <= 0.d0) then
! print the coordinates of the mesh points of this element
    print *, 'ispec = ', ispec
    print *, 'ngnod = ', ngnod
    do ia=1,ngnod
      nnum = knods(ia,ispec)
      xelm = coorg(1,nnum)
      zelm = coorg(2,nnum)
      print *,'node ', ia,' x,y = ',xelm,zelm
    enddo
    stop 'error: negative 2D Jacobian found'
  endif

! invert the relation
  xix = zgamma / jacobian
  gammax = - zxi / jacobian
  xiz = - xgamma / jacobian
  gammaz = xxi / jacobian

  end subroutine recompute_jacobian


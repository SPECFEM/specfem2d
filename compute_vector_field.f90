
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

  subroutine compute_vector_whole_medium(potential_acoustic,veloc_elastic,elastic,vector_field_display, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin)

! compute Grad(potential) in acoustic elements
! and combine with existing velocity vector field in elastic elements

  implicit none

  include "constants.h"

  integer nspec,npoin

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  double precision, dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

  logical, dimension(nspec) :: elastic
  double precision, dimension(npoin) :: potential_acoustic
  double precision, dimension(NDIM,npoin) :: veloc_elastic,vector_field_display

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz

! local variables
  integer i,j,ispec,iglob

! vector field in this element
  double precision, dimension(NDIM,NGLLX,NGLLX) :: vector_field_element

! loop over spectral elements
  do ispec = 1,nspec

! compute vector field in this element
    call compute_vector_one_element(vector_field_element,potential_acoustic,veloc_elastic,elastic, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,ispec)

! store the result
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        vector_field_display(:,iglob) = vector_field_element(:,i,j)
      enddo
    enddo

  enddo

  end subroutine compute_vector_whole_medium

!
!=====================================================================
!

  subroutine compute_vector_one_element(vector_field_element,potential_acoustic,veloc_elastic,elastic, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,ispec)

! compute Grad(potential) if acoustic element or copy existing vector if elastic element

  implicit none

  include "constants.h"

  integer nspec,npoin,ispec

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  double precision, dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

! vector field in this element
  double precision, dimension(NDIM,NGLLX,NGLLX) :: vector_field_element

  logical, dimension(nspec) :: elastic
  double precision, dimension(npoin) :: potential_acoustic
  double precision, dimension(NDIM,npoin) :: veloc_elastic

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz

! local variables
  integer i,j,k,iglob

! space derivatives
  double precision tempx1l,tempx2l
  double precision hp1,hp2

! jacobian
  double precision xixl,xizl,gammaxl,gammazl

! simple copy of existing vector if elastic element
  if(elastic(ispec)) then

    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        vector_field_element(1,i,j) = veloc_elastic(1,iglob)
        vector_field_element(2,i,j) = veloc_elastic(2,iglob)
      enddo
    enddo

! compute gradient of potential to calculate vector if acoustic element
    else

! double loop over GLL points to compute and store gradients
    do j = 1,NGLLZ
      do i = 1,NGLLX

! derivative along x
        tempx1l = ZERO
        do k = 1,NGLLX
          hp1 = hprime_xx(k,i)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_acoustic(iglob)*hp1
        enddo

! derivative along z
        tempx2l = ZERO
        do k = 1,NGLLZ
          hp2 = hprime_zz(k,j)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_acoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

! derivatives of potential
        vector_field_element(1,i,j) = tempx1l*xixl + tempx2l*gammaxl
        vector_field_element(2,i,j) = tempx1l*xizl + tempx2l*gammazl

      enddo
    enddo

  endif ! end of test if acoustic or elastic element

  end subroutine compute_vector_one_element


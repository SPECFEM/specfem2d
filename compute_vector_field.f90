
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!  Main authors: Dimitri Komatitsch, Nicolas Le Goff and Roland Martin
!                 University of Pau and CNRS, France
!
!        (c) University of Pau and CNRS, France, November 2007
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

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

  logical, dimension(nspec) :: elastic
  real(kind=CUSTOM_REAL), dimension(npoin) :: potential_acoustic
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin) :: veloc_elastic
  double precision, dimension(NDIM,npoin) :: vector_field_display

! array with derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

! local variables
  integer i,j,ispec,iglob

! vector field in this element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLX) :: vector_field_element

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

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

! vector field in this element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLX) :: vector_field_element

  logical, dimension(nspec) :: elastic
  real(kind=CUSTOM_REAL), dimension(npoin) :: potential_acoustic
  real(kind=CUSTOM_REAL), dimension(NDIM,npoin) :: veloc_elastic

! array with derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz

! local variables
  integer i,j,k,iglob

! space derivatives
  real(kind=CUSTOM_REAL) tempx1l,tempx2l
  real(kind=CUSTOM_REAL) hp1,hp2

! jacobian
  real(kind=CUSTOM_REAL) xixl,xizl,gammaxl,gammazl

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
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_acoustic(iglob)*hp1
        enddo

! derivative along z
        tempx2l = ZERO
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
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


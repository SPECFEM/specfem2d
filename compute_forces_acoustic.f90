
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
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

  subroutine compute_forces_acoustic(npoin,nspec,nelemabs,numat, &
               anyabs,assign_external_model,ibool,kmato,numabs, &
               elastic,codeabs,potential_dot_dot_acoustic,potential_dot_acoustic, &
               potential_acoustic,density,elastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,rhoext,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right, &
               nspec_outer, we_are_in_phase_outer)

! compute forces for the acoustic elements

  implicit none

  include "constants.h"

  integer :: npoin,nspec,nelemabs,numat

  logical :: anyabs,assign_external_model

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato
  integer, dimension(nelemabs) :: numabs,ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right

  logical, dimension(nspec) :: elastic
  logical, dimension(4,nelemabs)  :: codeabs

  real(kind=CUSTOM_REAL), dimension(npoin) :: potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic
  double precision, dimension(numat) :: density
  double precision, dimension(4,numat) :: elastcoef
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext,rhoext

! derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

! Gauss-Lobatto-Legendre weights
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: wzgll

! for overlapping MPI communications with computation
  integer, intent(in)  :: nspec_outer
!!!!!!!  integer, dimension(max(1,nspec_outer)), intent(in)  :: ispec_inner_outer_to_glob
  logical, intent(in)  :: we_are_in_phase_outer

!---
!--- local variables
!---

  integer :: ispec,i,j,k,iglob,ispecabs,ibegin,iend,jbegin,jend

! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,dux_dxl,dux_dzl
  real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2

! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_relaxed,lambdal_relaxed,kappal,cpl,rhol

! loop over spectral elements
! do ispec_inner_outer = 1,nspec_outer

!   ispec = ispec_inner_outer_to_glob(ispec_inner_outer)

  integer :: ifirstelem,ilastelem

  if(we_are_in_phase_outer) then
    ifirstelem = 1
    ilastelem = nspec_outer
  else
    ifirstelem = nspec_outer + 1
    ilastelem = nspec
  endif

  do ispec = ifirstelem,ilastelem

!---
!--- acoustic spectral element
!---
    if(.not. elastic(ispec)) then

      rhol = density(kmato(ispec))

! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

! derivative along x and along z
          dux_dxi = ZERO
          dux_dgamma = ZERO

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

! derivatives of potential
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          jacobianl = jacobian(i,j,ispec)

! if external density model
        if(assign_external_model) rhol = rhoext(i,j,ispec)

! for acoustic medium
! also add GLL integration weights
          tempx1(i,j) = wzgll(j)*jacobianl*(xixl*dux_dxl + xizl*dux_dzl) / rhol
          tempx2(i,j) = wxgll(i)*jacobianl*(gammaxl*dux_dxl + gammazl*dux_dzl) / rhol

        enddo
      enddo

!
! second double-loop over GLL to compute all the terms
!
      do j = 1,NGLLZ
        do i = 1,NGLLX

          iglob = ibool(i,j,ispec)

! along x direction and z direction
! and assemble the contributions
          do k = 1,NGLLX
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
                           (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
          enddo

        enddo ! second loop over the GLL points
      enddo

    endif ! end of test if acoustic element

    enddo ! end of loop over all spectral elements

! only for the first call to compute_forces_acoustic (during computation on outer elements)
  if ( we_are_in_phase_outer ) then

!
!--- absorbing boundaries
!
  if(anyabs) then

    do ispecabs=1,nelemabs

      ispec = numabs(ispecabs)

! get elastic parameters of current spectral element
      lambdal_relaxed = elastcoef(1,kmato(ispec))
      mul_relaxed = elastcoef(2,kmato(ispec))
      kappal  = lambdal_relaxed + TWO*mul_relaxed/3._CUSTOM_REAL
      rhol = density(kmato(ispec))
      cpl = sqrt((kappal + 4._CUSTOM_REAL*mul_relaxed/3._CUSTOM_REAL)/rhol)

!--- left absorbing boundary
      if(codeabs(ILEFT,ispecabs)) then

        i = 1

        jbegin = jbegin_left(ispecabs)
        jend = jend_left(ispecabs)

        do j = jbegin,jend

          iglob = ibool(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)

          weight = jacobian1D * wzgll(j)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec)) &
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl/rhol

        enddo

      endif  !  end of left absorbing boundary

!--- right absorbing boundary
      if(codeabs(IRIGHT,ispecabs)) then

        i = NGLLX

        jbegin = jbegin_right(ispecabs)
        jend = jend_right(ispecabs)

        do j = jbegin,jend

          iglob = ibool(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)

          weight = jacobian1D * wzgll(j)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec)) &
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl/rhol

        enddo

      endif  !  end of right absorbing boundary

!--- bottom absorbing boundary
      if(codeabs(IBOTTOM,ispecabs)) then

        j = 1

        ibegin = ibegin_bottom(ispecabs)
        iend = iend_bottom(ispecabs)

! exclude corners to make sure there is no contradiction on the normal
        if(codeabs(ILEFT,ispecabs)) ibegin = 2
        if(codeabs(IRIGHT,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend

          iglob = ibool(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)

          weight = jacobian1D * wxgll(i)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec)) &
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl/rhol

        enddo

      endif  !  end of bottom absorbing boundary

!--- top absorbing boundary
      if(codeabs(ITOP,ispecabs)) then

        j = NGLLZ

        ibegin = ibegin_top(ispecabs)
        iend = iend_top(ispecabs)

! exclude corners to make sure there is no contradiction on the normal
        if(codeabs(ILEFT,ispecabs)) ibegin = 2
        if(codeabs(IRIGHT,ispecabs)) iend = NGLLX-1

        do i = ibegin,iend

          iglob = ibool(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)

          weight = jacobian1D * wxgll(i)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec)) &
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl/rhol

        enddo

      endif  !  end of top absorbing boundary

    enddo

  endif  ! end of absorbing boundaries

  endif ! end of computation that needs to be done only once, during the first call to compute_forces_acoustic

  end subroutine compute_forces_acoustic


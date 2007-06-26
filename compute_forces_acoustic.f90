
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

  subroutine compute_forces_acoustic(npoin,nspec,nelemabs,numat, &
               iglob_source,ispec_selected_source,is_proc_source,source_type,it,NSTEP,anyabs, &
               assign_external_model,initialfield,ibool,kmato,numabs, &
               elastic,codeabs,potential_dot_dot_acoustic,potential_dot_acoustic, &
               potential_acoustic,density,elastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,vsext,rhoext,source_time_function,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right)

! compute forces for the acoustic elements

  implicit none

  include "constants.h"

  integer :: npoin,nspec,nelemabs,numat,iglob_source,ispec_selected_source,is_proc_source,source_type,it,NSTEP

  logical :: anyabs,assign_external_model,initialfield

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato
  integer, dimension(nelemabs) :: numabs,ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right

  logical, dimension(nspec) :: elastic
  logical, dimension(4,nelemabs)  :: codeabs

  double precision, dimension(npoin) :: potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic
  double precision, dimension(numat) :: density
  double precision, dimension(4,numat) :: elastcoef
  double precision, dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext,vsext,rhoext
  double precision, dimension(NSTEP) :: source_time_function

! derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

! Gauss-Lobatto-Legendre weights
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: wzgll

!---
!--- local variables
!---

  integer :: ispec,i,j,k,iglob,ispecabs,ibegin,iend,jbegin,jend

! spatial derivatives
  double precision :: dux_dxi,dux_dgamma,dux_dxl,dux_dzl
  double precision :: nx,nz,rho_vp,rho_vs,weight,xxi,zxi,xgamma,zgamma,jacobian1D

  double precision, dimension(NGLLX,NGLLZ) :: tempx1,tempx2

! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl

! material properties of the elastic medium
  double precision :: mul_relaxed,lambdal_relaxed,kappal,cpl,csl,rhol

! loop over spectral elements
  do ispec = 1,nspec

!---
!--- acoustic spectral element
!---
    if(.not. elastic(ispec)) then

! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

! derivative along x and along z
          dux_dxi = ZERO
          dux_dgamma = ZERO

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + potential_acoustic(ibool(k,j,ispec))*hprime_xx(k,i)
            dux_dgamma = dux_dgamma + potential_acoustic(ibool(i,k,ispec))*hprime_zz(k,j)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

! derivatives of potential
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          jacobianl = jacobian(i,j,ispec)

! for acoustic medium
! also add GLL integration weights
          tempx1(i,j) = wzgll(j)*jacobianl*(xixl*dux_dxl + xizl*dux_dzl)
          tempx2(i,j) = wxgll(i)*jacobianl*(gammaxl*dux_dxl + gammazl*dux_dzl)

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
                           (tempx1(k,j)*hprimewgll_xx(i,k) + tempx2(i,k)*hprimewgll_zz(j,k))
          enddo

        enddo ! second loop over the GLL points
      enddo

    endif ! end of test if acoustic element

    enddo ! end of loop over all spectral elements

!
!--- absorbing boundaries
!
  if(anyabs) then

    do ispecabs=1,nelemabs

      ispec = numabs(ispecabs)

! get elastic parameters of current spectral element
      lambdal_relaxed = elastcoef(1,kmato(ispec))
      mul_relaxed = elastcoef(2,kmato(ispec))
      rhol  = density(kmato(ispec))
      kappal  = lambdal_relaxed + TWO*mul_relaxed/3.d0
      cpl = sqrt((kappal + 4.d0*mul_relaxed/3.d0)/rhol)
      csl = sqrt(mul_relaxed/rhol)


!--- left absorbing boundary
      if(codeabs(ILEFT,ispecabs)) then

        i = 1

        jbegin = jbegin_left(ispecabs)
        jend = jend_left(ispecabs)

        do j = jbegin,jend

          iglob = ibool(i,j,ispec)

          zgamma = xix(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec)) then
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl
          endif

        enddo

      endif  !  end of left absorbing boundary

!--- right absorbing boundary
      if(codeabs(IRIGHT,ispecabs)) then

        i = NGLLX

        jbegin = jbegin_right(ispecabs)
        jend = jend_right(ispecabs)

        do j = jbegin,jend

          iglob = ibool(i,j,ispec)

          zgamma = xix(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec)) then
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl
          endif

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

          xxi = gammaz(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = + zxi / jacobian1D
          nz = - xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec)) then
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl
          endif

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

          xxi = gammaz(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
          endif

          rho_vp = rhol*cpl
          rho_vs = rhol*csl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = + zxi / jacobian1D
          nz = - xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec)) then
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl
          endif

        enddo

      endif  !  end of top absorbing boundary

    enddo

  endif  ! end of absorbing boundaries


! --- add the source
  if(.not. initialfield) then

     if (is_proc_source == 1 ) then
! collocated force
! beware, for acoustic medium, source is a pressure source
        if(source_type == 1) then
           if(.not. elastic(ispec_selected_source)) then
              potential_dot_dot_acoustic(iglob_source) = potential_dot_dot_acoustic(iglob_source) + source_time_function(it)
           endif

! moment tensor
        else if(source_type == 2) then

           if(.not. elastic(ispec_selected_source)) stop 'cannot have moment tensor source in acoustic element'

        endif
     end if
  else
     stop 'wrong source type'
  endif

  end subroutine compute_forces_acoustic


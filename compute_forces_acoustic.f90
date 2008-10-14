
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.3
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour and CNRS, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT gps DOT caltech DOT edu
!               Jeroen Tromp, jtromp aT gps DOT caltech DOT edu
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

  subroutine compute_forces_acoustic(npoin,nspec,myrank,numat, &
               iglob_source,ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver,&
               source_type,it,NSTEP,anyabs, &
               assign_external_model,initialfield,ibool,kmato, &
               elastic,poroelastic,potential_dot_dot_acoustic,potential_dot_acoustic, &
               potential_acoustic,b_potential_dot_dot_acoustic,b_potential_acoustic,&
               density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,source_time_function,adj_sourcearrays,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right, &
               nspec_inner_outer, ispec_inner_outer_to_glob, num_phase_inner_outer, &
               nrec,isolver,save_forward,b_absorb_acoustic_left,&
               b_absorb_acoustic_right,b_absorb_acoustic_bottom,&
               b_absorb_acoustic_top,nspec_xmin,nspec_xmax,&
               nspec_zmin,nspec_zmax,ib_xmin,ib_xmax,ib_zmin,ib_zmax,kappa_ac_k)

! compute forces for the acoustic elements

  implicit none

  include "constants.h"

  integer :: npoin,nspec,myrank,numat,iglob_source,ispec_selected_source,is_proc_source,source_type,it,NSTEP
  integer :: nrec,isolver
  integer, dimension(nrec) :: ispec_selected_rec,which_proc_receiver
  integer :: nspec_xmin,nspec_xmax,nspec_zmin,nspec_zmax
  integer, dimension(nspec_xmin) :: ib_xmin,jbegin_left,jend_left
  integer, dimension(nspec_xmax) :: ib_xmax,jbegin_right,jend_right
  integer, dimension(nspec_zmin) :: ib_zmin,ibegin_bottom,iend_bottom
  integer, dimension(nspec_zmax) :: ib_zmax,ibegin_top,iend_top

  logical :: anyabs,assign_external_model,initialfield
  logical :: save_forward

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato

  logical, dimension(nspec) :: elastic,poroelastic

  real(kind=CUSTOM_REAL), dimension(npoin) :: potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic
  real(kind=CUSTOM_REAL), dimension(npoin) :: b_potential_dot_dot_acoustic,b_potential_acoustic
  double precision, dimension(2,numat) :: density
  double precision, dimension(4,3,numat) :: poroelastcoef
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext
  real(kind=CUSTOM_REAL), dimension(NSTEP) :: source_time_function

  real(kind=CUSTOM_REAL), dimension(nrec,NSTEP,NDIM,NGLLX,NGLLZ) :: adj_sourcearrays
  real(kind=CUSTOM_REAL), dimension(npoin) :: kappa_ac_k
  double precision, dimension(NGLLZ,nspec_xmin,NSTEP) :: b_absorb_acoustic_left
  double precision, dimension(NGLLZ,nspec_xmax,NSTEP) :: b_absorb_acoustic_right
  double precision, dimension(NGLLX,nspec_zmax,NSTEP) :: b_absorb_acoustic_top
  double precision, dimension(NGLLX,nspec_zmin,NSTEP) :: b_absorb_acoustic_bottom
! derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

! Gauss-Lobatto-Legendre weights
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: wzgll

! for overlapping MPI communications with computation
  integer, intent(in)  :: nspec_inner_outer
  integer, dimension(max(1,nspec_inner_outer)), intent(in)  :: ispec_inner_outer_to_glob
  logical, intent(in)  :: num_phase_inner_outer

!---
!--- local variables
!---

  integer :: ispec,ispec_inner_outer,i,j,k,iglob,ispecabs,ibegin,iend,jbegin,jend,irec,irec_local

! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,dux_dxl,dux_dzl
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_dux_dxl,b_dux_dzl
  real(kind=CUSTOM_REAL) :: weight,xxi,zxi,xgamma,zgamma,jacobian1D

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: b_tempx1,b_tempx2

! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl,nx,nz

! material properties of the acoustic medium
  real(kind=CUSTOM_REAL) :: kappal,cpl,rhol,rho_vp

! loop over spectral elements
  do ispec_inner_outer = 1,nspec_inner_outer

    ispec = ispec_inner_outer_to_glob(ispec_inner_outer)

!---
!--- acoustic spectral element
!---
    if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then

! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

! derivative along x and along z
          dux_dxi = ZERO
          dux_dgamma = ZERO

          b_dux_dxi = ZERO
          b_dux_dgamma = ZERO

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + potential_acoustic(ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + potential_acoustic(ibool(i,k,ispec))*hprime_zz(j,k)

            if(isolver == 2) then
            b_dux_dxi = b_dux_dxi + b_potential_acoustic(ibool(k,j,ispec))*hprime_xx(k,i)
            b_dux_dgamma = b_dux_dgamma + b_potential_acoustic(ibool(i,k,ispec))*hprime_zz(k,j)
            endif

          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

! derivatives of potential
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

            if(isolver == 2) then
          b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
          b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl
            endif

! kernels calculation
   if(isolver == 2) then
          iglob = ibool(i,j,ispec)
          kappa_ac_k(iglob) = dux_dxl *  b_dux_dxl
   endif

          jacobianl = jacobian(i,j,ispec)

! for acoustic medium
! also add GLL integration weights
          tempx1(i,j) = wzgll(j)*jacobianl*(xixl*dux_dxl + xizl*dux_dzl)
          tempx2(i,j) = wxgll(i)*jacobianl*(gammaxl*dux_dxl + gammazl*dux_dzl)

            if(isolver == 2) then
          b_tempx1(i,j) = wzgll(j)*jacobianl*(xixl*b_dux_dxl + xizl*b_dux_dzl)
          b_tempx2(i,j) = wxgll(i)*jacobianl*(gammaxl*b_dux_dxl + gammazl*b_dux_dzl)
            endif

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
            if(isolver == 2) then
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                           (b_tempx1(k,j)*hprimewgll_xx(i,k) + b_tempx2(i,k)*hprimewgll_zz(j,k))
            endif
          enddo

        enddo ! second loop over the GLL points
      enddo

    endif ! end of test if acoustic element

    enddo ! end of loop over all spectral elements

! only for the first call to compute_forces_acoustic (during computation on outer elements)
  if ( num_phase_inner_outer ) then

!
!--- absorbing boundaries
!
  if(anyabs) then

!--- left absorbing boundary
      if( nspec_xmin > 0 ) then

      do ispecabs = 1, nspec_xmin

      ispec = ib_xmin(ispecabs)

! get parameters of current spectral element
! acoustic (fluid) properties
    kappal = poroelastcoef(1,2,kmato(ispec)) 
    rhol = density(2,kmato(ispec))

      cpl = sqrt(kappal/rhol)

        i = 1

        jbegin = jbegin_left(ispecabs)
        jend = jend_left(ispecabs)

        do j = jbegin,jend

          iglob = ibool(i,j,ispec)

          zgamma = xix(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
          endif

          rho_vp = rhol*cpl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = - zgamma / jacobian1D
          nz = + xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl

             if(save_forward .and. isolver ==1) then
            b_absorb_acoustic_left(j,ispecabs,it) = potential_dot_acoustic(iglob)*weight/cpl
             elseif(isolver == 2) then
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                               b_absorb_acoustic_left(j,ispecabs,NSTEP-it+1)
             endif

          endif

         enddo

        enddo

      endif  !  end of left absorbing boundary

!--- right absorbing boundary
      if( nspec_xmax > 0 ) then
        
      do ispecabs = 1, nspec_xmax

      ispec = ib_xmax(ispecabs)

! get parameters of current spectral element
! acoustic (fluid) properties
    kappal = poroelastcoef(1,2,kmato(ispec)) 
    rhol = density(2,kmato(ispec))

      cpl = sqrt(kappal/rhol)

        i = NGLLX

        jbegin = jbegin_right(ispecabs)
        jend = jend_right(ispecabs)

        do j = jbegin,jend

          iglob = ibool(i,j,ispec)

          zgamma = xix(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
          endif

          rho_vp = rhol*cpl

          xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
          zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xgamma**2 + zgamma**2)
          nx = + zgamma / jacobian1D
          nz = - xgamma / jacobian1D

          weight = jacobian1D * wzgll(j)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl

             if(save_forward .and. isolver ==1) then
            b_absorb_acoustic_right(j,ispecabs,it) = potential_dot_acoustic(iglob)*weight/cpl
             elseif(isolver == 2) then
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                              b_absorb_acoustic_right(j,ispecabs,NSTEP-it+1)
             endif

          endif

         enddo

        enddo

      endif  !  end of right absorbing boundary

!--- bottom absorbing boundary
      if( nspec_zmin > 0) then

      do ispecabs = 1, nspec_zmin

      ispec = ib_zmin(ispecabs)

! get parameters of current spectral element
! acoustic (fluid) properties
    kappal = poroelastcoef(1,2,kmato(ispec)) 
    rhol = density(2,kmato(ispec))

      cpl = sqrt(kappal/rhol)

        j = 1

        ibegin = ibegin_bottom(ispecabs)
        iend = iend_bottom(ispecabs)

! exclude corners to make sure there is no contradiction on the normal
        if( nspec_xmin > 0 ) ibegin = 2
        if( nspec_xmax > 0 ) iend = NGLLX-1

        do i = ibegin,iend

          iglob = ibool(i,j,ispec)

          xxi = gammaz(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
          endif

          rho_vp = rhol*cpl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = + zxi / jacobian1D
          nz = - xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl

             if(save_forward .and. isolver ==1) then
            b_absorb_acoustic_bottom(j,ispecabs,it) = potential_dot_acoustic(iglob)*weight/cpl
             elseif(isolver == 2) then
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                               b_absorb_acoustic_bottom(j,ispecabs,NSTEP-it+1)
             endif

          endif

         enddo

        enddo

      endif  !  end of bottom absorbing boundary

!--- top absorbing boundary
      if( nspec_zmax > 0 ) then

      do ispecabs = 1, nspec_zmax

      ispec = ib_zmax(ispecabs)

! get parameters of current spectral element
! acoustic (fluid) properties
    kappal = poroelastcoef(1,2,kmato(ispec)) 
    rhol = density(2,kmato(ispec))

      cpl = sqrt(kappal/rhol)

        j = NGLLZ

        ibegin = ibegin_top(ispecabs)
        iend = iend_top(ispecabs)

! exclude corners to make sure there is no contradiction on the normal
        if( nspec_xmin > 0) ibegin = 2
        if( nspec_xmax > 0) iend = NGLLX-1

        do i = ibegin,iend

          iglob = ibool(i,j,ispec)

          xxi = gammaz(i,j,ispec) * jacobian(i,j,ispec)

! external velocity model
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
          endif

          rho_vp = rhol*cpl

          xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
          zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
          jacobian1D = sqrt(xxi**2 + zxi**2)
          nx = - zxi / jacobian1D
          nz = + xxi / jacobian1D

          weight = jacobian1D * wxgll(i)

! Sommerfeld condition if acoustic
          if(.not. elastic(ispec) .and. .not. poroelastic(ispec)) then
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_acoustic(iglob)*weight/cpl

             if(save_forward .and. isolver ==1) then
            b_absorb_acoustic_top(j,ispecabs,it) = potential_dot_acoustic(iglob)*weight/cpl
             elseif(isolver == 2) then
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - b_absorb_acoustic_top(j,ispecabs,NSTEP-it+1)
             endif

          endif

         enddo

        enddo

      endif  !  end of top absorbing boundary

  endif  ! end of absorbing boundaries

! --- add the source
  if(.not. initialfield) then

     if (is_proc_source == 1 ) then
! collocated force
! beware, for acoustic medium, source is a pressure source
        if(source_type == 1) then
           if(.not. elastic(ispec_selected_source) .and. .not. poroelastic(ispec_selected_source)) then

      if(isolver == 1) then  ! forward wavefield
      potential_dot_dot_acoustic(iglob_source) = potential_dot_dot_acoustic(iglob_source) + source_time_function(it)
      else                   ! backward wavefield
      b_potential_dot_dot_acoustic(iglob_source) = b_potential_dot_dot_acoustic(iglob_source) + source_time_function(NSTEP-it+1)
      endif

           endif

! moment tensor
        else if(source_type == 2) then

           if(.not. elastic(ispec_selected_source) .and. .not. poroelastic(ispec_selected_source)) then
              call exit_MPI('cannot have moment tensor source in acoustic element')
           endif
        endif
     endif

    if(isolver == 2) then   ! adjoint wavefield
      irec_local = 0
      do irec = 1,nrec
!   add the source (only if this proc carries the source)
      if(myrank == which_proc_receiver(irec)) then
           if(.not. elastic(ispec_selected_rec(irec)) .and. .not. poroelastic(ispec_selected_rec(irec))) then
      irec_local = irec_local + 1
! add source array
      do j=1,NGLLZ
        do i=1,NGLLX
      iglob = ibool(i,j,ispec_selected_rec(irec))
      potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + &
          adj_sourcearrays(irec_local,NSTEP-it+1,1,i,j)
        enddo
      enddo
            endif
      endif ! if this processor carries the adjoint source
      enddo ! irec = 1,nrec
    endif ! isolver == 2 adjoint wavefield

  else
     call exit_MPI('wrong source type')
  endif

  endif ! end of computation that needs to be done only once, during the first call to compute_forces_acoustic

  end subroutine compute_forces_acoustic


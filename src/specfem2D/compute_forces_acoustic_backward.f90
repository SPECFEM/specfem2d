
!========================================================================
!
!                  S P E C F E M 2 D  Version 7 . 0
!                  --------------------------------
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

  subroutine compute_forces_acoustic_backward(b_potential_dot_dot_acoustic,b_potential_acoustic)


  ! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use specfem_par, only: nglob,nspec,nelemabs,it,NSTEP, &
                         assign_external_model,ibool,kmato,numabs,acoustic, &
                         codeabs,codeabs_corner, &
                         density,poroelastcoef,xix,xiz,gammax,gammaz,jacobian, &
                         vpext,rhoext, &
                         hprime_xx,hprimewgll_xx, &
                         hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         AXISYM,coord, is_on_the_axis,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj, &
                         ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
                         ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2, &
                         ib_left,ib_right,ib_bottom,ib_top, &
                         b_absorb_acoustic_left,b_absorb_acoustic_right, &
                         b_absorb_acoustic_bottom,b_absorb_acoustic_top, &
                         STACEY_BOUNDARY_CONDITIONS

  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(nglob) :: b_potential_dot_dot_acoustic, b_potential_acoustic

  ! local parameters
  integer :: ispec,i,j,k,iglob,ispecabs,ibegin,iend,jbegin,jend
  integer :: ifirstelem,ilastelem

  ! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,dux_dxl,dux_dzl
  real(kind=CUSTOM_REAL) :: xxi

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  ! material properties of the acoustic medium
  real(kind=CUSTOM_REAL) :: mul_relaxed,lambdal_relaxed,kappal,cpl,rhol

  ifirstelem = 1
  ilastelem = nspec


! loop over spectral elements
  do ispec = ifirstelem,ilastelem
    ! acoustic spectral element
    if( acoustic(ispec) ) then
      rhol = density(1,kmato(ispec))

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! derivative along x and along z
          dux_dxi = 0._CUSTOM_REAL; dux_dgamma = 0._CUSTOM_REAL

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            if( AXISYM ) then
              if( is_on_the_axis(ispec) ) then
                dux_dxi = dux_dxi + b_potential_acoustic(ibool(k,j,ispec)) * hprimeBar_xx(i,k)
              else
                dux_dxi = dux_dxi + b_potential_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)
              endif
            else
              dux_dxi = dux_dxi + b_potential_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)
            endif
            dux_dgamma = dux_dgamma + b_potential_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)

          ! derivatives of potential
          dux_dxl = dux_dxi * xixl + dux_dgamma * gammaxl
          dux_dzl = dux_dxi * xizl + dux_dgamma * gammazl

          if( AXISYM .and. is_on_the_axis(ispec) .and. i == 1 ) then ! dchi/dr=rho * u_r=0 on the axis
            dux_dxl = ZERO
          endif

          jacobianl = jacobian(i,j,ispec)

          ! if external density model
          if( assign_external_model ) then
            rhol = rhoext(i,j,ispec)
          endif

          if( AXISYM ) then
            if( is_on_the_axis(ispec) .and. i == 1 ) then
              xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
              r_xiplus1(i,j) = xxi
            else if( is_on_the_axis(ispec) ) then
              r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i)+ONE)
            endif
          endif

          ! for acoustic medium also add integration weights
          if( AXISYM ) then
            if( is_on_the_axis(ispec) ) then
              tempx1(i,j) = wzgll(j) * r_xiplus1(i,j) * jacobianl * (xixl * dux_dxl + xizl * dux_dzl) / rhol
              tempx2(i,j) = wxglj(i) * r_xiplus1(i,j) * jacobianl * (gammaxl * dux_dxl + gammazl * dux_dzl) / rhol
            else
              tempx1(i,j) = wzgll(j) * coord(1,ibool(i,j,ispec)) * jacobianl * (xixl * dux_dxl + xizl * dux_dzl) / rhol
              tempx2(i,j) = wxgll(i) * coord(1,ibool(i,j,ispec)) * jacobianl * (gammaxl * dux_dxl + gammazl * dux_dzl) / rhol
            endif
          else
            tempx1(i,j) = wzgll(j) * jacobianl * (xixl * dux_dxl + xizl * dux_dzl) / rhol
            tempx2(i,j) = wxgll(i) * jacobianl * (gammaxl * dux_dxl + gammazl * dux_dzl) / rhol
          endif
        enddo
      enddo

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          if( assign_external_model ) then
            rhol = rhoext(i,j,ispec)
            cpl = vpext(i,j,ispec)
            !assuming that in fluid(acoustic) part input cpl is defined by sqrt(kappal/rhol), &
            !which is not the same as in cpl input in elastic part
            kappal = rhol * cpl * cpl
          else
            lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
            mul_relaxed = poroelastcoef(2,1,kmato(ispec))
            kappal  = lambdal_relaxed + TWO * mul_relaxed/3._CUSTOM_REAL
            rhol = density(1,kmato(ispec))
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
            if( AXISYM ) then
              if( is_on_the_axis(ispec) ) then
                b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                       (tempx1(k,j) * hprimeBarwglj_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
              else
                b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                         (tempx1(k,j) * hprimewgll_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
              endif
            else
              b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                       (tempx1(k,j) * hprimewgll_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
            endif
          enddo
        enddo ! second loop over the GLL points
      enddo

    endif ! end of test if acoustic element
  enddo ! end of loop over all spectral elements
!
!--- absorbing boundaries
!

! The outer boundary condition to use for PML elements in fluid layers is Neumann for the potential
! because we need Dirichlet conditions for the displacement vector, which means Neumann for the potential.
! Thus, there is nothing to enforce explicitly here.
! There is something to enforce explicitly only in the case of elastic elements, for which a Dirichlet
! condition is needed for the displacement vector, which is the vectorial unknown for these elements.

! for Stacey paraxial absorbing conditions (more precisely: Sommerfeld in the case of a fluid) we implement them here

  if( STACEY_BOUNDARY_CONDITIONS ) then
    do ispecabs=1,nelemabs
      ispec = numabs(ispecabs)
      ! Sommerfeld condition if acoustic
      if( acoustic(ispec) ) then
        !--- left absorbing boundary
        if( codeabs(IEDGE4,ispecabs) ) then
          i = 1
          jbegin = ibegin_edge4(ispecabs)
          jend = iend_edge4(ispecabs)
          do j = jbegin,jend
            iglob = ibool(i,j,ispec)
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                                  b_absorb_acoustic_left(j,ib_left(ispecabs),NSTEP-it+1)

          enddo
        endif  !  end of left absorbing boundary

        !--- right absorbing boundary
        if( codeabs(IEDGE2,ispecabs) ) then
          i = NGLLX
          jbegin = ibegin_edge2(ispecabs)
          jend = iend_edge2(ispecabs)
          do j = jbegin,jend
            iglob = ibool(i,j,ispec)
            ! adds (previously) stored contribution
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                                  b_absorb_acoustic_right(j,ib_right(ispecabs),NSTEP-it+1)
          enddo
        endif  !  end of right absorbing boundary

        !--- bottom absorbing boundary
        if( codeabs(IEDGE1,ispecabs) ) then
          j = 1
          ibegin = ibegin_edge1(ispecabs)
          iend = iend_edge1(ispecabs)
          ! exclude corners to make sure there is no contradiction on the normal
          if( codeabs_corner(1,ispecabs)) ibegin = 2
          if( codeabs_corner(2,ispecabs)) iend = NGLLX-1
          do i = ibegin,iend
            iglob = ibool(i,j,ispec)
            ! adds (previously) stored contribution
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                                  b_absorb_acoustic_bottom(i,ib_bottom(ispecabs),NSTEP-it+1)
          enddo
        endif  !  end of bottom absorbing boundary

        !--- top absorbing boundary
        if( codeabs(IEDGE3,ispecabs) ) then
          j = NGLLZ
          ibegin = ibegin_edge3(ispecabs)
          iend = iend_edge3(ispecabs)
          ! exclude corners to make sure there is no contradiction on the normal
          if( codeabs_corner(3,ispecabs)) ibegin = 2
          if( codeabs_corner(4,ispecabs)) iend = NGLLX-1
          do i = ibegin,iend
            iglob = ibool(i,j,ispec)
            b_potential_dot_dot_acoustic(iglob) = b_potential_dot_dot_acoustic(iglob) - &
                                                  b_absorb_acoustic_top(i,ib_top(ispecabs),NSTEP-it+1)
          enddo
        endif  !  end of top absorbing boundary
      endif ! acoustic ispec
    enddo
  endif  ! end of absorbing boundaries

 end subroutine compute_forces_acoustic_backward

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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

  subroutine compute_forces_acoustic(potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic, &
                                     PML_BOUNDARY_CONDITIONS,potential_acoustic_old,iphase)


! compute forces in the acoustic elements in forward simulation and in adjoint simulation in adjoint inversion

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,CPML_X_ONLY,CPML_Z_ONLY,IRIGHT,ILEFT,IBOTTOM,ITOP, &
    ZERO,ONE,TWO,TWO_THIRDS,IEDGE1,IEDGE2,IEDGE3,IEDGE4, &
    ALPHA_LDDRK,BETA_LDDRK,C_LDDRK

  use specfem_par, only: nglob, &
                         assign_external_model,ibool,kmato,ispec_is_acoustic, &
                         density,poroelastcoef,eta,rhoext,vpext,etaext, &
                         xix,xiz,gammax,gammaz,jacobian, &
                         hprime_xx,hprimewgll_xx, &
                         hprime_zz,hprimewgll_zz,wxgll,wzgll, &
                         AXISYM,is_on_the_axis,coord,hprimeBar_xx,hprimeBarwglj_xx,xiglj,wxglj

  ! overlapping communication
  use specfem_par, only: nspec_inner_acoustic,nspec_outer_acoustic,phase_ispec_inner_acoustic

  ! PML arrays
  use specfem_par, only: ispec_is_PML

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob),intent(inout) :: potential_dot_dot_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob),intent(in) :: potential_dot_acoustic,potential_acoustic

  logical,intent(in) :: PML_BOUNDARY_CONDITIONS
  real(kind=CUSTOM_REAL), dimension(nglob) :: potential_acoustic_old

  integer,intent(in) :: iphase

  ! local parameters
  integer :: ispec,i,j,k,iglob

  ! spatial derivatives
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dux_dxi,dux_dgamma
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dux_dxl,dux_dzl
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dvx_dxi,dvx_dgamma
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: dvx_dxl,dvx_dzl

  ! Bulk viscosity derivatives
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: deta_dxi,deta_dgamma
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: deta_dxl,deta_dzl

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: potential_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: d_potential_elem
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1_visco,tempx2_visco
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: eta_on_rhocsquare
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: temp3_visco

  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLLZ) :: r_xiplus1

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL), dimension(8,NGLLX,NGLLZ) :: deriv
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  real(kind=CUSTOM_REAL) :: rhol,cpl,etal,lambda_relaxed,mu_relaxed,kappal,fac
  real(kind=CUSTOM_REAL) :: temp1l,temp2l,temp1l_visco,temp2l_visco

  ! local PML parameters
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: potential_dot_dot_acoustic_PML

  integer :: num_elements,ispec_p

  ! choses inner/outer elements
  if (iphase == 1) then
    num_elements = nspec_outer_acoustic
  else
    num_elements = nspec_inner_acoustic
  endif

  ! loop over spectral elements
  do ispec_p = 1,num_elements

    ! returns element id from stored element list
    ispec = phase_ispec_inner_acoustic(ispec_p,iphase)

    ! only for acoustic spectral elements
    if (.not. ispec_is_acoustic(ispec)) cycle

    ! gets local potential for element
    rhol = density(1,kmato(ispec))
    etal = eta(1,kmato(ispec))

    lambda_relaxed = poroelastcoef(1,1,kmato(ispec))
    mu_relaxed = poroelastcoef(2,1,kmato(ispec))
    kappal = lambda_relaxed + mu_relaxed

    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        potential_elem(i,j) = potential_acoustic(iglob)
        d_potential_elem(i,j) = potential_dot_acoustic(iglob)

! Note from Quentin Brissaud, ISAE, France, April 2017 about viscoacoustic fluids:
! De la meme maniere que pour la partie viscoelastique, on realise le calul d'un facteur Q constant,
! en revanche je n'utilise que les solides correspondant a une geometrie isotropique
! (tout est donc sur la diagonale du tenseur des contraintes) et il n'y a que Qkappa
! car il n'y a pas d'ondes de cisaillement dans le fluide.

        ! stores local array for element xi/gamma/jacobian (for better performance)
        deriv(1,i,j) = xix(i,j,ispec)
        deriv(2,i,j) = xiz(i,j,ispec)
        deriv(3,i,j) = gammax(i,j,ispec)
        deriv(4,i,j) = gammaz(i,j,ispec)
        deriv(5,i,j) = jacobian(i,j,ispec)
        ! if external density model
        if (assign_external_model) then
          rhol = rhoext(i,j,ispec)
          etal = etaext(i,j,ispec)
          cpl  = vpext(i,j,ispec)
          kappal  = rhol*cpl**2
        endif
        ! viscosity
        etal = 0. !! DK DK and Quentin Brissaud: viscosity "eta" is ignored for now
        deriv(6,i,j) = jacobian(i,j,ispec) / rhol
        deriv(7,i,j) = jacobian(i,j,ispec) * etal/(rhol*kappal)
        ! deriv(8,i,j) = jacobian(i,j,ispec) * ( tau_eps_fluid(i,j,ispec) /  tau_sig_fluid(i,j,ispec) )

        eta_on_rhocsquare(i,j) = etal/kappal

! Note from Quentin Brissaud, ISAE, France, April 2017 about viscosity "eta" for fluids, which is ignored for now:
! DK DK the comment below is now a bit obsolete, the new implementation of viscoacoustic fluids based
! on N_SLS Zener bodies put in parallel, as in solid regions, being much better for most applications.
!
! Pour inclure la bulk viscosity pour les fluides, il faut regarder mon rapport de stage d'il y a quelques annees
! qui montre comment est implementee la viscosite (voir equation (41) du rapport).
! La viscosite telle qu'implementee peut varier selon z (profondeur ou altitude).
! J'ai rajoute deux variables globales: eta et etaext qui correspondent respectivement a la viscosite
! si aucun fichier externe de modele n'est fourni et a la viscosit& provenant d'un modele externe.
! Je n'ai pas fait les modifications afin de mettre la viscosite dans le fichier de "Par_file"
! (lorsqu'on definit les modeles a la main). J'ai simplement fait un test avec une viscosite constante
! que j'avais defini a la main dans le fichier "compute_forces_acoustic.f90" et dans le fichier "compute_coupling_acoustic_el.f90"
! et cela dampe bien l'amplitude (ainsi qu'implique un shift en frequence).
!
! Dans ce code officiel fusionn'e de SPECFEM2D, J'ai laisse la possibilite d'ajouter un jour la viscosite bulk (facteur "eta")
! mais je n'ai pas mis l'option dans le "Par_file", c'est pour l'instant mis a la main dans les fichiers
! "compute_forces_acoustic" et "coupling_acoustic_el". Il faudrait donc realiser trois etapes:
!   - Il faut creer une variable globale correspondant a la viscosite bulk et remplacer dans les fichiers
!       "compute_forces_acoustic" et "coupling_acoustic_el" les variables locales "eta"
!   - Il faut modifier la lecture de fichiers externes "define_external_model" afin de lire la viscosite bulk
!   - Pour l'ajouter en tant que parametre dans le "Par_file" il serait necessaire de modifier un peu la creation
!       de modele afin qu'on puisse indiquer la viscosite bulk (a la place de Qmu par exemple).

      enddo
    enddo

    ! first double loop over GLL points to compute and store gradients
    call mxm_2comp_singleA(dux_dxi,dux_dgamma,potential_elem,hprime_xx,hprime_zz)
    call mxm_2comp_singleA(dvx_dxi,dvx_dgamma,d_potential_elem,hprime_xx,hprime_zz)
    call mxm_2comp_singleA(deta_dxi,deta_dgamma,eta_on_rhocsquare,hprime_xx,hprime_zz)

    ! AXISYM case overwrites dux_dxi
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            ! derivative along x and along z
            dux_dxi(i,j) = 0._CUSTOM_REAL
            do k = 1,NGLLX
              dux_dxi(i,j) = dux_dxi(i,j) + potential_elem(k,j) * hprimeBar_xx(i,k)
            enddo
          enddo
        enddo
      endif
    endif

    ! gets derivatives of ux and uz with respect to x and z
    do j = 1,NGLLZ
      do i = 1,NGLLX
        xixl = deriv(1,i,j)
        xizl = deriv(2,i,j)
        gammaxl = deriv(3,i,j)
        gammazl = deriv(4,i,j)

        ! derivatives of potential
        dux_dxl(i,j) = dux_dxi(i,j) * xixl + dux_dgamma(i,j) * gammaxl
        dux_dzl(i,j) = dux_dxi(i,j) * xizl + dux_dgamma(i,j) * gammazl

        dvx_dxl(i,j) = dvx_dxi(i,j) * xixl + dvx_dgamma(i,j) * gammaxl
        dvx_dzl(i,j) = dvx_dxi(i,j) * xizl + dvx_dgamma(i,j) * gammazl

        deta_dxl(i,j) = deta_dxi(i,j) * xixl + deta_dgamma(i,j) * gammaxl
        deta_dzl(i,j) = deta_dxi(i,j) * xizl + deta_dgamma(i,j) * gammazl
      enddo
    enddo

    ! AXISYM case overwrite dux_dxl
    if (AXISYM) then
      if (is_on_the_axis(ispec)) then
        ! dchi/dr=rho * u_r=0 on the axis
        ! i == 1
        do j = 1,NGLLZ
          dux_dxl(1,j) = 0._CUSTOM_REAL
        enddo
      endif
    endif

    ! derivative along x and along zbb
    if (PML_BOUNDARY_CONDITIONS) then
      call pml_compute_memory_variables_acoustic(ispec,nglob,potential_acoustic_old,dux_dxl,dux_dzl)
    endif

    ! first double loop to compute gradient
    if (AXISYM) then
      ! AXISYM case
      if (is_on_the_axis(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            xixl = deriv(1,i,j)
            xizl = deriv(2,i,j)
            gammaxl = deriv(3,i,j)
            gammazl = deriv(4,i,j)
            jacobianl = deriv(5,i,j)
            fac = deriv(6,i,j) ! jacobian/rho

            if (i == 1) then
              ! dchi/dr=rho * u_r=0 on the axis
              dux_dxl(i,j) = 0._CUSTOM_REAL
              r_xiplus1(i,j) = gammaz(i,j,ispec) * jacobianl
            else
              r_xiplus1(i,j) = coord(1,ibool(i,j,ispec))/(xiglj(i) + ONE)
            endif
            tempx1(i,j) = r_xiplus1(i,j) * fac * (xixl * dux_dxl(i,j) + xizl * dux_dzl(i,j))
            tempx2(i,j) = r_xiplus1(i,j) * fac * (gammaxl * dux_dxl(i,j) + gammazl * dux_dzl(i,j))
          enddo
        enddo
      else
        do j = 1,NGLLZ
          do i = 1,NGLLX
            xixl = deriv(1,i,j)
            xizl = deriv(2,i,j)
            gammaxl = deriv(3,i,j)
            gammazl = deriv(4,i,j)
            jacobianl = deriv(5,i,j)
            fac = deriv(6,i,j) ! jacobian/rho

            tempx1(i,j) = coord(1,ibool(i,j,ispec)) * fac * (xixl * dux_dxl(i,j) + xizl * dux_dzl(i,j))
            tempx2(i,j) = coord(1,ibool(i,j,ispec)) * fac * (gammaxl * dux_dxl(i,j) + gammazl * dux_dzl(i,j))
          enddo
        enddo
      endif
    else
      ! default case
      do j = 1,NGLLZ
        do i = 1,NGLLX
          xixl = deriv(1,i,j)
          xizl = deriv(2,i,j)
          gammaxl = deriv(3,i,j)
          gammazl = deriv(4,i,j)
          jacobianl = deriv(5,i,j)
          fac = deriv(6,i,j) ! jacobian/rho

          tempx1(i,j) = fac * (xixl * dux_dxl(i,j) + xizl * dux_dzl(i,j))
          tempx2(i,j) = fac * (gammaxl * dux_dxl(i,j) + gammazl * dux_dzl(i,j))

          temp3_visco(i,j)  = fac * (dvx_dxl(i,j)*deta_dxl(i,j) + dvx_dzl(i,j)*deta_dzl(i,j))

          fac = deriv(7,i,j)
          tempx1_visco(i,j) = fac * (xixl * dvx_dxl(i,j) + xizl * dvx_dzl(i,j))
          tempx2_visco(i,j) = fac * (gammaxl * dvx_dxl(i,j) + gammazl * dvx_dzl(i,j))

        enddo
      enddo
    endif

    ! first double loop over GLL points to compute and store gradients
    if (PML_BOUNDARY_CONDITIONS) then
      ! calculates contribution from each C-PML element to update acceleration
      call pml_compute_accel_contribution_acoustic(ispec,nglob, &
                                                   potential_acoustic,potential_acoustic_old,potential_dot_acoustic, &
                                                   potential_dot_dot_acoustic_PML,r_xiplus1)
    endif

!
! second double-loop over GLL to compute all the terms
!

    ! along x direction and z direction
    ! and assemble the contributions
    if (AXISYM) then
      ! axisymmetric case
      if (is_on_the_axis(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            ! assembles the contributions
            temp1l = 0._CUSTOM_REAL
            temp2l = 0._CUSTOM_REAL
            do k = 1,NGLLX
              temp1l = temp1l + tempx1(k,j) * hprimeBarwglj_xx(k,i)
              temp2l = temp2l + tempx2(i,k) * hprimewgll_zz(k,j)
            enddo
            ! sums contributions from each element to the global values
            iglob = ibool(i,j,ispec)
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                - (wzgll(j) * temp1l + wxglj(i) * temp2l)
          enddo
        enddo
      else
        do j = 1,NGLLZ
          do i = 1,NGLLX
            ! assembles the contributions
            temp1l = 0._CUSTOM_REAL
            temp2l = 0._CUSTOM_REAL
            do k = 1,NGLLX
              temp1l = temp1l + tempx1(k,j) * hprimewgll_xx(k,i)
              temp2l = temp2l + tempx2(i,k) * hprimewgll_zz(k,j)
            enddo
            ! sums contributions from each element to the global values
            iglob = ibool(i,j,ispec)
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                                - (wzgll(j) * temp1l + wxgll(i) * temp2l)
          enddo
        enddo
      endif
    else
      ! default case
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! assembles the contributions
          temp1l = 0._CUSTOM_REAL
          temp2l = 0._CUSTOM_REAL
          temp1l_visco = 0._CUSTOM_REAL
          temp2l_visco = 0._CUSTOM_REAL
          do k = 1,NGLLX
            temp1l = temp1l + tempx1(k,j) * hprimewgll_xx(k,i)
            temp2l = temp2l + tempx2(i,k) * hprimewgll_zz(k,j)

            temp1l_visco = temp1l_visco + tempx1_visco(k,j) * hprimewgll_xx(k,i)
            temp2l_visco = temp2l_visco + tempx2_visco(i,k) * hprimewgll_zz(k,j)
          enddo
          ! sums contributions from each element to the global values
          iglob = ibool(i,j,ispec)
          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) &
                                              - (wzgll(j) * temp1l + wxgll(i) * temp2l) &
                                              - (wzgll(j) * temp1l_visco + wxgll(i) * temp2l_visco) &
                                              - wzgll(j) * wxgll(i) * temp3_visco(i,j)
        enddo
      enddo
    endif

    ! PML contribution
    if (PML_BOUNDARY_CONDITIONS) then
      if (ispec_is_PML(ispec)) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - potential_dot_dot_acoustic_PML(i,j)
          enddo
        enddo
      endif
    endif

  enddo ! end of loop over all spectral elements

  contains

!---------------------------------------------------------------------------------------

  subroutine mxm_2comp_singleA(x,z,A,B,C)

! matrix x matrix multiplication, merging 2 loops for x = A^t B^t and z = A C^t
!
! index notation:
! general matrix multiplication: uij = (A B)ij = Aik Bkj
!                          here: xij = (A^t B^t)ij = Akj Bik = (B A)ij
!                                zij = (A C^t)ij = Aik Cjk
!
! original loops:
!
!      do j = 1,NGLLZ
!        do i = 1,NGLLX
!          ! derivative along x and along z
!          dux_dxi(i,j) = 0._CUSTOM_REAL
!          dux_dgamma(i,j) = 0._CUSTOM_REAL
!
!          ! first double loop over GLL points to compute and store gradients
!          ! we can merge the two loops because NGLLX == NGLLZ
!          do k = 1,NGLLX
!            dux_dxi(i,j) = dux_dxi(i,j) + potential_elem(k,j) * hprime_xx(i,k)
!            dux_dgamma(i,j) = dux_dgamma(i,j) + potential_elem(i,k) * hprime_zz(j,k)
!          enddo
!        enddo
!      enddo

  use constants, only: NGLLX,NGLLZ,CUSTOM_REAL

  implicit none

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(out) :: x,z
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ),intent(in) :: A,B,C

  ! local parameters
  integer :: i,j,k

  select case(NGLLX)
  case (5)
    do j = 1,5
      do i = 1,5
        ! loop unrolling
        x(i,j) = A(1,j) * B(i,1) + A(2,j) * B(i,2) + A(3,j) * B(i,3) + A(4,j) * B(i,4) + A(5,j) * B(i,5)
        z(i,j) = A(i,1) * C(j,1) + A(i,2) * C(j,2) + A(i,3) * C(j,3) + A(i,4) * C(j,4) + A(i,5) * C(j,5)
      enddo
    enddo

  case default
    do j = 1,NGLLZ
      do i = 1,NGLLX
        x(i,j) = 0._CUSTOM_REAL
        z(i,j) = 0._CUSTOM_REAL
        do k = 1,NGLLX
          x(i,j) = x(i,j) + A(k,j) * B(i,k)
          z(i,j) = z(i,j) + A(i,k) * C(j,k)
        enddo
      enddo
    enddo
  end select

  end subroutine mxm_2comp_singleA

  end subroutine compute_forces_acoustic

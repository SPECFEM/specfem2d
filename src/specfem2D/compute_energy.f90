!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
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

  subroutine compute_and_output_energy()

  use constants, only: IOUT_ENERGY,CUSTOM_REAL

  use specfem_par, only: myrank,it,DT,kinetic_energy,potential_energy,t0,NTSTEP_BETWEEN_OUTPUT_ENERGY

  implicit none

  ! local parameters
  real(kind=CUSTOM_REAL) :: kinetic_energy_total,potential_energy_total

  ! checks if anything to do
  if (mod(it,NTSTEP_BETWEEN_OUTPUT_ENERGY) /= 0) return

  ! computes energy
  call compute_energy()

  ! computes total for all processes
  call sum_all_cr(kinetic_energy,kinetic_energy_total)
  call sum_all_cr(potential_energy,potential_energy_total)

  ! saves kinetic, potential and total energy for this time step in external file
  if (myrank == 0) then
    ! format: #time  #E_kin(kinetic energy)  #E_pot(potential energy)  #E_tot(total energy)
    write(IOUT_ENERGY,*) real(dble(it-1)*DT - t0,4),real(kinetic_energy_total,4), &
                         real(potential_energy_total,4),real(kinetic_energy_total + potential_energy_total,4)
  endif

  end subroutine compute_and_output_energy

!
!----------------------------------------------------------------------------------------
!

  subroutine compute_energy()

! compute kinetic and potential energy in the solid (acoustic elements are excluded)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NGLJ,NDIM,ZERO,TWO

  use specfem_par, only: AXISYM,is_on_the_axis,nspec,kinetic_energy,potential_energy, &
                         ibool,hprime_xx,hprime_zz,hprimeBar_xx,xix,xiz,gammax,gammaz,jacobian,wxgll,wzgll, &
                         mustore,rho_vpstore,rhostore, &
                         phistore,tortstore,kappaarraystore,mufr_store,rhoarraystore, &
                         ispec_is_poroelastic,ispec_is_elastic, &
                         P_SV,ispec_is_PML, &
                         GPU_MODE,any_acoustic,any_elastic,any_poroelastic

  ! wavefields
  use specfem_par, only: displ_elastic,veloc_elastic,accel_elastic, &
                         displs_poroelastic,displw_poroelastic,velocs_poroelastic,velocw_poroelastic, &
                         potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

  use specfem_par_gpu, only: Mesh_pointer,NGLOB_AB

  implicit none

! local variables
  integer :: i,j,k,ispec
  real(kind=CUSTOM_REAL) :: cpl
  real(kind=CUSTOM_REAL) :: mu_G,lambdal_G,lambdalplus2mul_G

  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl
  double precision :: rhol
  ! to evaluate cpI, cpII, and cs, and rI (poroelastic medium)
  double precision :: phi,tort,mu_fr,kappa_s,kappa_f,kappa_fr
  double precision :: rho_s,rho_f,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot
  double precision :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic

  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  double precision :: dwx_dxl,dwz_dxl,dwx_dzl,dwz_dzl
  ! vector field in an element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: vector_field_element
  ! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: pressure_element

  ! transfers fields
  if (GPU_MODE) then
    ! a simple workaround to avoid implementing the following routine in CUDA.
    ! be aware that these memory transfers will go through the bottleneck of memory bandwidth between CPU & GPU,
    ! thus slow down the simulation
    ! acoustic domains
    if (any_acoustic) then
      call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                          Mesh_pointer)
    endif
    ! elastic domains
    if (any_elastic) then
      call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ_elastic,veloc_elastic,accel_elastic,Mesh_pointer)
    endif
    ! poroelastic domains
    if (any_poroelastic) then
      stop 'Poroelastic domain transfers for GPU_MODE in compute_energy() not implemented yet'
    endif
  endif

  ! initialization
  ! please do *NOT* remove it
  kinetic_energy = 0._CUSTOM_REAL
  potential_energy = 0._CUSTOM_REAL

  ! loop over spectral elements
  do ispec = 1,nspec

    ! only compute energy in the main domain, not in the PMLs, in which it is not physical
    if (ispec_is_PML(ispec)) cycle

    !---
    !--- elastic spectral element
    !---
    if (ispec_is_elastic(ispec)) then

      ! double loop over GLL points
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! get elastic parameters of current grid point
          mul_unrelaxed_elastic = mustore(i,j,ispec)
          rhol = rhostore(i,j,ispec)
          cpl = rho_vpstore(i,j,ispec) / rhol

          lambdal_unrelaxed_elastic = rhol*cpl*cpl - TWO * mul_unrelaxed_elastic
          lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + TWO * mul_unrelaxed_elastic

          ! derivative along x and along z
          dux_dxi = 0._CUSTOM_REAL
          duz_dxi = 0._CUSTOM_REAL

          dux_dgamma = 0._CUSTOM_REAL
          duz_dgamma = 0._CUSTOM_REAL

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          if (AXISYM) then
            ! axisymmetric case
            if (is_on_the_axis(ispec)) then
              do k = 1,NGLJ
                dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
              enddo
            else
              do k = 1,NGLJ
                dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
                duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
              enddo
            endif
          else
            ! default, non-axisymmetric case
            do k = 1,NGLLX
              dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
              duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            enddo
          endif

          do k = 1,NGLLX
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          jacobianl = jacobian(i,j,ispec)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          ! elastic medium
          ! kinetic energy:
          !   E_kin = 1/2 * m * v**2
          !           -> 1/2 * rho * (v_x**2 + v_z**2) * w_i * w_j * J_ij
          !
          ! potential energy:
          !   E_pot = 1/2 * stress * strain
          !         = 1/2 * ( T_xx * eps_xx
          !                 + T_zz * eps_zz
          !                 + 2 T_xz * eps_xz)                  (since T_xz == T_zx and eps_xz == eps_zx)
          !         = 1/2 * ( [(lambda + 2 mu) eps_xx + lambda eps_zz] * eps_xx
          !                 + [(lambda + 2 mu) eps_zz + lambda eps_xx] * eps_zz
          !                 + 2 [mu 2 eps_xz] * eps_xz  )       (note: eps_xz = 1/2 (duz_dx + dux_dz))
          !         = 1/2 * ( (lambda + 2 mu)*dux_dx * dux_dx
          !                 + (lambda + 2 mu)*duz_dx * duz_dz
          !                 + 2 lambda * dux_dx * duz_dz
          !                 + 2 mu * (dux_dz + duz_dx) * 1/2 (dux_dz + duz_dx) )
          !
          if (P_SV) then
            ! P-SV waves
            ! compute kinetic energy
            ! TODO ABAB This integral over space should be adapted for axisym geometries if needed
            kinetic_energy = kinetic_energy  &
                + rhol * (veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) &
                * wxgll(i) * wzgll(j) * jacobianl / TWO

            ! compute potential energy
            ! TODO ABAB This integral over space should be adapted for axisym geometries if needed
            potential_energy = potential_energy &
                + ( lambdaplus2mu_unrelaxed_elastic * dux_dxl**2 &
                  + lambdaplus2mu_unrelaxed_elastic * duz_dzl**2 &
                  + TWO*lambdal_unrelaxed_elastic * dux_dxl*duz_dzl &
                  + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2 ) &
                * wxgll(i) * wzgll(j) * jacobianl / TWO
          else
            ! SH-waves
            ! compute kinetic energy
            ! note: only component array(1,..) is non-zero for SH waves and corresponds to y-direction displacement
            kinetic_energy = kinetic_energy  &
                + rhol * (veloc_elastic(1,ibool(i,j,ispec))**2 ) &
                * wxgll(i) * wzgll(j) * jacobianl / TWO

            ! compute potential energy
            ! note: T_xy = T_yx = mu 2 eps_yx = mu (duy_dx + dux_dy) = mu duy_dx     (sets dux_dy == 0)
            !       T_zy = T_yz = mu 2 eps_yz = mu (duy_dz + duz_dy) = mu duy_dz
            !
            !      -> E_pot = 1/2 * stress * strain
            !               = 1/2 * (T_yx * eps_yx + T_yz * eps_yz)
            !               = 1/2 * (mu duy_dx * 1/2 duy_dx + mu duy_dz * 1/2 duy_dz)
            !               = 1/2 * (mu (duy_dx**2 + duy_dz**2) / 2 )
            !
            ! reuses variable names from P-SV case, here: dux_dx == duy_dx and dux_dz == duy_dz
            potential_energy = potential_energy &
                + ( mul_unrelaxed_elastic * (dux_dxl**2 + dux_dzl**2) / TWO )&
                * wxgll(i) * wzgll(j) * jacobianl / TWO
          endif
        enddo
      enddo

    !---
    !--- poroelastic spectral element
    !---
    else if (ispec_is_poroelastic(ispec)) then

      ! get unrelaxed elastic parameters of current spectral element
      !for now replaced by solid, fluid, and frame parameters of current spectral element

      ! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! derivative along x and along z
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

          dwx_dxi = ZERO
          dwz_dxi = ZERO

          dwx_dgamma = ZERO
          dwz_dgamma = ZERO

          ! first double loop over GLL points to compute and store gradients
          ! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displs_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displs_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displs_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displs_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)


            dwx_dxi = dwx_dxi + displw_poroelastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            dwz_dxi = dwz_dxi + displw_poroelastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
            dwx_dgamma = dwx_dgamma + displw_poroelastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            dwz_dgamma = dwz_dgamma + displw_poroelastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
          enddo

          xixl = xix(i,j,ispec)
          xizl = xiz(i,j,ispec)
          gammaxl = gammax(i,j,ispec)
          gammazl = gammaz(i,j,ispec)
          jacobianl = jacobian(i,j,ispec)

          ! derivatives of displacement
          dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

          dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
          dwx_dzl = dwx_dxi*xizl + dwx_dgamma*gammazl

          dwz_dxl = dwz_dxi*xixl + dwz_dgamma*gammaxl
          dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

          ! gets poroelastic material
          phi = phistore(i,j,ispec)
          tort = tortstore(i,j,ispec)
          kappa_s = kappaarraystore(1,i,j,ispec)
          kappa_f = kappaarraystore(2,i,j,ispec)
          kappa_fr = kappaarraystore(3,i,j,ispec)
          mu_fr = mufr_store(i,j,ispec)

          ! Biot coefficients for the input phi
          call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

          !The RHS has the form : div T -phi/c div T_f + phi/ceta_fk^-1.partial t w
          !where T = G:grad u_s + C div w I
          !and T_f = C div u_s I + M div w I
          !we are expressing lambdaplus2mu, lambda, and mu for G, C, and M
          mu_G = mu_fr
          lambdal_G = H_biot - TWO*mu_fr
          lambdalplus2mul_G = lambdal_G + TWO*mu_G

          rho_s = rhoarraystore(1,i,j,ispec)
          rho_f = rhoarraystore(2,i,j,ispec)
          rho_bar = (1.d0 - phi)*rho_s + phi * rho_f

          ! compute potential energy
          potential_energy = potential_energy &
              + ( lambdalplus2mul_G*dux_dxl**2 &
              + lambdalplus2mul_G*duz_dzl**2 &
              + TWO*lambdal_G*dux_dxl*duz_dzl + mu_G*(dux_dzl + duz_dxl)**2 &
              + TWO*C_biot*dwx_dxl*dux_dxl + TWO*C_biot*dwz_dzl*duz_dzl &
              + TWO*C_biot*(dwx_dxl*duz_dzl + dwz_dzl*dux_dxl) &
              + M_biot*dwx_dxl**2 + M_biot*dwz_dzl**2 &
              + TWO*M_biot*dwx_dxl*dwz_dzl )*wxgll(i)*wzgll(j)*jacobianl / TWO

          ! compute kinetic energy
          if (phi > 0.0d0) then
            kinetic_energy = kinetic_energy &
              + ( rho_bar*(velocs_poroelastic(1,ibool(i,j,ispec))**2 &
              + velocs_poroelastic(2,ibool(i,j,ispec))**2) &
              + rho_f*tort/phi*(velocw_poroelastic(1,ibool(i,j,ispec))**2 &
              + velocw_poroelastic(2,ibool(i,j,ispec))**2) &
              + rho_f*(velocs_poroelastic(1,ibool(i,j,ispec))*velocw_poroelastic(1,ibool(i,j,ispec)) &
              + velocs_poroelastic(2,ibool(i,j,ispec))*velocw_poroelastic(2,ibool(i,j,ispec))) &
                 )*wxgll(i)*wzgll(j)*jacobianl / TWO
          else
            kinetic_energy = kinetic_energy  &
              + rho_s*(velocs_poroelastic(1,ibool(i,j,ispec))**2 &
              + velocs_poroelastic(2,ibool(i,j,ispec))**2)*wxgll(i)*wzgll(j)*jacobianl / TWO
          endif
        enddo
      enddo

    !---
    !--- acoustic spectral element
    !---
    else

      ! for the definition of potential energy in an acoustic fluid, see for instance
      ! equation (23) of M. Maess et al., Journal of Sound and Vibration 296 (2006) 264-276

      ! in case of an acoustic medium, a potential Chi of (density * displacement) is used as in Chaljub and Valette,
      ! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
      ! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
      ! This permits acoustic-elastic coupling based on a non-iterative time scheme.
      ! Displacement is then: u = grad(Chi) / rho
      ! Velocity is then: v = grad(Chi_dot) / rho (Chi_dot being the time derivative of Chi)
      ! and pressure is: p = - Chi_dot_dot  (Chi_dot_dot being the time second derivative of Chi).

      ! compute pressure in this element
      call compute_pressure_one_element(ispec,pressure_element,displ_elastic,displs_poroelastic,displw_poroelastic, &
                                        potential_dot_dot_acoustic,potential_acoustic)

      ! compute velocity vector field in this element
      call compute_vector_one_element(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,ispec,vector_field_element)

      ! double loop over GLL points
      do j = 1,NGLLZ
        do i = 1,NGLLX
          ! get elastic parameters of current grid point
          rhol = rhostore(i,j,ispec)
          cpl = rho_vpstore(i,j,ispec)/rhol

          jacobianl = jacobian(i,j,ispec)

          ! compute kinetic energy
          kinetic_energy = kinetic_energy &
              + rhol*(vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) *wxgll(i)*wzgll(j)*jacobianl / TWO

          ! compute potential energy
          potential_energy = potential_energy &
              + (pressure_element(i,j)**2)*wxgll(i)*wzgll(j)*jacobianl / (TWO * rhol * cpl**2)
        enddo
      enddo

    endif
  enddo

  end subroutine compute_energy

!
!----------------------------------------------------------------------------------------
!

  subroutine compute_integrated_energy_field_and_output()

  ! compute int_0^t v^2 dt and write it on file if needed

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,IIN,MAX_STRING_LEN,OUTPUT_FILES

  use specfem_par, only: myrank,it,coord,nspec,ibool,integrated_kinetic_energy_field,max_kinetic_energy_field, &
                         integrated_potential_energy_field,max_potential_energy_field,kinetic_effective_duration_field, &
                         potential_effective_duration_field,total_integrated_energy_field,max_total_energy_field, &
                         total_effective_duration_field,NTSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP

  implicit none

  ! local variables
  integer :: ispec,iglob
!!!! DK DK commenting this out for now because "call execute_command_line" is Fortran 2008
!!!! DK DK and some older compilers do not support it yet. We can probably put it back in a few years.
  !integer :: statMkdir
  character(len=MAX_STRING_LEN)  :: filename

  !! ABAB Uncomment to write the velocity profile in acoustic part
  !real(kind=CUSTOM_REAL) :: cpl
  !double precision :: rhol
  !double precision :: lambdal_unrelaxed_elastic
  !! ABAB

  ! computes maximum energy and integrated energy fields
  call compute_energy_fields()

  ! Create directories
  if (it == 1) then
    if (myrank == 0) then
!!!! DK DK commenting this out for now because "call execute_command_line" is Fortran 2008
!!!! DK DK and some older compilers do not support it yet. We can probably put it back in a few years.
     !call execute_command_line('mkdir -p '//trim(OUTPUT_FILES)//'energyFields',wait = .true.,cmdstat = statMkdir)
     !if (statMkdir /= 0) call exit_MPI(myrank,'Impossible to create '//trim(OUTPUT_FILES)//'energyFields')

     !call execute_command_line('mkdir -p '//trim(OUTPUT_FILES)//'energyFields/kinetic',wait = .true.,cmdstat = statMkdir)
     !if (statMkdir /= 0) call exit_MPI(myrank,'Impossible to create '//trim(OUTPUT_FILES)//'energyFields/kinetic')

     !call execute_command_line('mkdir -p '//trim(OUTPUT_FILES)//'energyFields/potential',wait = .true.,cmdstat = statMkdir)
     !if (statMkdir /= 0) call exit_MPI(myrank,'Impossible to create '//trim(OUTPUT_FILES)//'energyFields/potential')

     !call execute_command_line('mkdir -p '//trim(OUTPUT_FILES)//'energyFields/total',wait = .true.,cmdstat = statMkdir)
     !if (statMkdir /= 0) call exit_MPI(myrank,'Impossible to create '//trim(OUTPUT_FILES)//'energyFields/total')
    endif
    !call synchronize_all() ! Wait for first proc to create directories
  endif

  if (mod(it,NTSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
    ! write integrated kinetic energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/kinetic/integrated_kinetic_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(integrated_kinetic_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write max kinetic energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/kinetic/max_kinetic_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(max_kinetic_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write integrated potential energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/potential/integrated_potential_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(integrated_potential_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write max potential energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/potential/max_potential_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(max_potential_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write potential effective duration field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/potential/potential_effective_duration_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(potential_effective_duration_field(ispec),4)
    enddo
    close(IIN)

   ! write kinetic effective duration field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/kinetic/kinetic_effective_duration_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(kinetic_effective_duration_field(ispec),4)
    enddo
    close(IIN)

    ! write total integrated energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/total/total_integrated_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(total_integrated_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write max total energy field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/total/max_total_energy_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(max_total_energy_field(ispec),4)
    enddo
    close(IIN)

    ! write total effective duration field in external file
    write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'energyFields/total/total_effective_duration_field',myrank
    open(unit=IIN,file=trim(filename),status='unknown',action='write')
    ! loop over spectral elements
    do ispec = 1,nspec
      iglob = ibool(2,2,ispec)
      write(IIN,*) real(coord(1,iglob),4), &
                   real(coord(2,iglob),4),real(total_effective_duration_field(ispec),4)
    enddo
    close(IIN)
  endif

  ! ABAB Uncomment to write the velocity profile in the acoustic part in file
  !
  !  write(filename,"(a,i5.5)") trim(OUTPUT_FILES)//'velocities',myrank
  !  open(unit=IIN,file=trim(filename),status='unknown',action='write')
  !
  !  if (mod(it,NTSTEP_BETWEEN_OUTPUT_SEISMOS) == 0 .or. it == NSTEP) then
  !    ! loop over spectral elements
  !    do ispec = 1,nspec
  !      if (ispec_is_acoustic(ispec)) then
  !        ! get density of current grid point
  !        cpl = rho_vpstore(2,2,ispec)/rhostore(2,2,ispec)
  !        iglob = ibool(2,2,ispec)
  !        write(IIN,*) real(coord(2,iglob),4),cpl
  !      endif
  !    enddo
  !  endif
  !  close(IIN)

  end subroutine compute_integrated_energy_field_and_output

!
!----------------------------------------------------------------------------------------
!

  subroutine compute_energy_fields()

  ! computes maximum, integrated energy and duration fields

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM,TWO,ZERO

  use specfem_par, only: AXISYM,is_on_the_axis,nspec,ibool,deltat, &
                         ispec_is_elastic,ispec_is_poroelastic, &
                         integrated_kinetic_energy_field,max_kinetic_energy_field, &
                         integrated_potential_energy_field,max_potential_energy_field,kinetic_effective_duration_field, &
                         potential_effective_duration_field,total_integrated_energy_field,max_total_energy_field, &
                         total_effective_duration_field, &
                         mustore,rho_vpstore,rhostore, &
                         jacobian, &
                         hprime_xx,hprime_zz,hprimeBar_xx,xix,xiz,gammax,gammaz, &
                         GPU_MODE,any_acoustic,any_elastic,any_poroelastic, &
                         P_SV

  ! wavefields
  use specfem_par, only: displ_elastic,veloc_elastic,accel_elastic, &
                         displs_poroelastic,displw_poroelastic,velocs_poroelastic, &
                         potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic

  use specfem_par_gpu, only: Mesh_pointer,NGLOB_AB

  implicit none

  ! local variables
  integer :: ispec,i,j,k
  real(kind=CUSTOM_REAL) :: cpl

  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl
  double precision :: rhol
  double precision :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic
  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl

  ! vector field in an element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: vector_field_element
  ! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: pressure_element

  ! transfers fields
  if (GPU_MODE) then
    ! a simple workaround to avoid implementing the following routine in CUDA.
    ! be aware that these memory transfers will go through the bottleneck of memory bandwidth between CPU & GPU,
    ! thus slow down the simulation
    ! acoustic domains
    if (any_acoustic) then
      call transfer_fields_ac_from_device(NGLOB_AB,potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic, &
                                          Mesh_pointer)
    endif
    ! elastic domains
    if (any_elastic) then
      call transfer_fields_el_from_device(NDIM*NGLOB_AB,displ_elastic,veloc_elastic,accel_elastic,Mesh_pointer)
    endif
    ! poroelastic domains
    if (any_poroelastic) then
      stop 'Poroelastic domain transfers for GPU_MODE in compute_energy_fields() not implemented yet'
    endif
  endif

  ! We save the value at the GLL point:
  i=2
  j=2

  ! loop over spectral elements
  do ispec = 1,nspec

    !---
    !--- elastic spectral element
    !---
    if (ispec_is_elastic(ispec)) then

      ! get elastic parameters of current grid point
      mul_unrelaxed_elastic = mustore(i,j,ispec)
      rhol = rhostore(i,j,ispec)
      cpl = rho_vpstore(i,j,ispec)/rhol

      lambdal_unrelaxed_elastic = rhol*cpl*cpl - TWO * mul_unrelaxed_elastic
      lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + TWO * mul_unrelaxed_elastic

      ! derivative along x and along z
      dux_dxi = 0._CUSTOM_REAL
      duz_dxi = 0._CUSTOM_REAL

      dux_dgamma = 0._CUSTOM_REAL
      duz_dgamma = 0._CUSTOM_REAL

      ! first double loop over GLL points to compute and store gradients
      ! we can merge the two loops because NGLLX == NGLLZ
      if (AXISYM) then
        ! axisymmetric case
        if (is_on_the_axis(ispec)) then
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprimeBar_xx(i,k)
          enddo
        else
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
          enddo
        endif
      else
        ! default, non-axisymmetric case
        do k = 1,NGLLX
          dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
          duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
        enddo
      endif
      do k = 1,NGLLX
        dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
        duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
      enddo

      xixl = xix(i,j,ispec)
      xizl = xiz(i,j,ispec)
      gammaxl = gammax(i,j,ispec)
      gammazl = gammaz(i,j,ispec)
      jacobianl = jacobian(i,j,ispec)

      ! derivatives of displacement
      dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
      dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

      duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
      duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

      ! compute total integrated energy
      ! We record just one point per element (i=2, j=2)

      ! kinetic energy
      if (P_SV) then
        ! P-SV waves
        integrated_kinetic_energy_field(ispec) = integrated_kinetic_energy_field(ispec)  &
            +  rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) * deltat * 0.5_CUSTOM_REAL

        ! maximum value
        if (max_kinetic_energy_field(ispec) <&
              rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) * 0.5_CUSTOM_REAL) then
          max_kinetic_energy_field(ispec) = &
              rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) * 0.5_CUSTOM_REAL
        endif
      else
        ! SH waves
        integrated_kinetic_energy_field(ispec) = integrated_kinetic_energy_field(ispec)  &
            +  rhol*(veloc_elastic(1,ibool(i,j,ispec))**2) * deltat * 0.5_CUSTOM_REAL

        ! maximum value
        if (max_kinetic_energy_field(ispec) <&
              rhol*(veloc_elastic(1,ibool(i,j,ispec))**2) * 0.5_CUSTOM_REAL) then
          max_kinetic_energy_field(ispec) = &
              rhol*(veloc_elastic(1,ibool(i,j,ispec))**2) * 0.5_CUSTOM_REAL
        endif
      endif

      ! maximum E_kin
      if (max_kinetic_energy_field(ispec) > ZERO) then
        kinetic_effective_duration_field(ispec) = TWO*integrated_kinetic_energy_field(ispec)/max_kinetic_energy_field(ispec)
      endif

      ! potential energy
      if (P_SV) then
        ! P-SV waves
        integrated_potential_energy_field(ispec) = integrated_potential_energy_field(ispec) &
                + (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
                + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
                + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
                + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2) * deltat * 0.5_CUSTOM_REAL

        ! maximum value
        if (max_potential_energy_field(ispec) < (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
              + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
              + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
              + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2) / TWO) then
          max_potential_energy_field(ispec) = (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
            + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
            + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
            + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2) / TWO
        endif
      else
        ! SH waves
        integrated_potential_energy_field(ispec) = integrated_potential_energy_field(ispec) &
                + (mul_unrelaxed_elastic * (dux_dxl**2 + dux_dzl**2) / TWO ) * deltat * 0.5_CUSTOM_REAL

        ! maximum value
        if (max_potential_energy_field(ispec) < (mul_unrelaxed_elastic * (dux_dxl**2 + dux_dzl**2) / TWO ) / TWO ) then
          max_potential_energy_field(ispec) = (mul_unrelaxed_elastic * (dux_dxl**2 + dux_dzl**2) / TWO ) / TWO
        endif
      endif

      if (max_potential_energy_field(ispec) > ZERO) then
        potential_effective_duration_field(ispec) = TWO*integrated_potential_energy_field(ispec) / &
                                                    max_potential_energy_field(ispec)
      endif

      ! maximum total energy
      if (P_SV) then
        ! P-SV waves
        if (max_total_energy_field(ispec) < (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
              + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
              + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
              + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2) / TWO + &
              rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) / TWO) then
          max_total_energy_field(ispec) = (lambdaplus2mu_unrelaxed_elastic*dux_dxl**2 &
              + lambdaplus2mu_unrelaxed_elastic*duz_dzl**2 &
              + TWO*lambdal_unrelaxed_elastic*dux_dxl*duz_dzl &
              + mul_unrelaxed_elastic*(dux_dzl + duz_dxl)**2) / TWO + &
              rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) / TWO
        endif
      else
        ! SH waves
        if (max_total_energy_field(ispec) < (mul_unrelaxed_elastic * (dux_dxl**2 + dux_dzl**2) / TWO ) / TWO + &
                rhol*(veloc_elastic(1,ibool(i,j,ispec))**2) / TWO) then
            max_total_energy_field(ispec) = (mul_unrelaxed_elastic * (dux_dxl**2 + dux_dzl**2) / TWO ) / TWO + &
                rhol*(veloc_elastic(1,ibool(i,j,ispec))**2) / TWO
        endif
      endif

    !---
    !--- poroelastic spectral element
    !---
    else if (ispec_is_poroelastic(ispec)) then
       ! safety check
       call stop_the_code( &
'COMPUTE_INTEGRATED_ENERGY_FIELD is not available for poroelastic media yet (but it would be very easy to implement)')

    !---
    !--- acoustic spectral element
    !---
    else

      ! compute velocity vector field in this element
      call compute_vector_one_element(potential_dot_acoustic,veloc_elastic,velocs_poroelastic,ispec,vector_field_element)

      ! compute pressure in this element
      call compute_pressure_one_element(ispec,pressure_element,displ_elastic,displs_poroelastic,displw_poroelastic, &
                                        potential_dot_dot_acoustic,potential_acoustic)

      ! get density of current grid point
      rhol = rhostore(i,j,ispec)
      cpl = rho_vpstore(i,j,ispec) / rhol

      jacobianl = jacobian(i,j,ispec)

      ! compute total integrated energy ! = int_0^t v^2 dt
      integrated_kinetic_energy_field(ispec) = integrated_kinetic_energy_field(ispec)  &
           +  rhol * (vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) * deltat * 0.5_CUSTOM_REAL

      if (max_kinetic_energy_field(ispec) <&
            rhol * (vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) * 0.5_CUSTOM_REAL) then
        max_kinetic_energy_field(ispec) = &
            rhol * (vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) * 0.5_CUSTOM_REAL
      endif

      if (max_kinetic_energy_field(ispec) > ZERO) then
        kinetic_effective_duration_field(ispec) = TWO*integrated_kinetic_energy_field(ispec)/max_kinetic_energy_field(ispec)
      endif

      ! compute potential energy
      integrated_potential_energy_field(ispec) = integrated_potential_energy_field(ispec) + pressure_element(i,j)**2 * deltat &
                                                 / (TWO*rhol*cpl**2)
      if (max_potential_energy_field(ispec) < pressure_element(i,j)**2 / (TWO*rhol*cpl**2)) then
        max_potential_energy_field(ispec) = pressure_element(i,j)**2 / (TWO*rhol*cpl**2)
      endif

      if (max_potential_energy_field(ispec) > ZERO) then
        potential_effective_duration_field(ispec) = TWO*integrated_potential_energy_field(ispec) / &
                                                    max_potential_energy_field(ispec)
      endif

      if (max_total_energy_field(ispec) < pressure_element(i,j)**2 / (TWO*rhol*cpl**2) + &
            rhol * (vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) / TWO) then
        max_total_energy_field(ispec) = pressure_element(i,j)**2 / (TWO*rhol*cpl**2) + &
            rhol * (vector_field_element(1,i,j)**2 + vector_field_element(2,i,j)**2) / TWO
      endif

    endif

    total_integrated_energy_field(ispec) = integrated_kinetic_energy_field(ispec) + integrated_potential_energy_field(ispec)

    if (max_total_energy_field(ispec) > ZERO) then
      total_effective_duration_field(ispec) = TWO*total_integrated_energy_field(ispec) / max_total_energy_field(ispec)
    endif

  enddo

  end subroutine compute_energy_fields



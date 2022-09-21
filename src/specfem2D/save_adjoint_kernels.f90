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

  subroutine save_adjoint_kernels()

! saves adjoint sensitivity kernels to file

  use constants, only: IMAIN,SAVE_WEIGHTS

  use specfem_par, only: myrank, any_acoustic, any_elastic, any_poroelastic

  implicit none

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'Writing Kernels file'
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! elastic kernels
  if (any_elastic) call save_kernels_elastic()

  ! acoustic kernels
  if (any_acoustic) call save_kernels_acoustic()

  ! poroelastic kernels
  if (any_poroelastic) call save_kernels_poroelastic()

  ! save weights for volume integration,
  ! in order to benchmark the kernels with analytical expressions
  if (SAVE_WEIGHTS) call save_weights_kernel()

  end subroutine save_adjoint_kernels

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_elastic()

  use constants, only: NGLLX,NGLLZ,IMAIN,CUSTOM_REAL,FOUR_THIRDS,TWO_THIRDS,TWO,MAX_STRING_LEN, &
    OUTPUT_FILES

  use specfem_par, only: myrank, nspec, ibool, coord, save_ASCII_kernels, &
                         rho_kl, kappa_kl, mu_kl, rhop_kl, alpha_kl, beta_kl, &
                         bulk_c_kl, bulk_beta_kl, &
                         rhorho_el_Hessian_final1, rhorho_el_Hessian_final2, &
                         rhostore,mustore,kappastore, &
                         ispec_is_elastic, &
                         deltat,NTSTEP_BETWEEN_COMPUTE_KERNELS, &
                         APPROXIMATE_HESS_KL

  use specfem_par, only: ispec_is_anisotropic, c11_kl, c13_kl, c15_kl, c33_kl, c35_kl, c55_kl

  implicit none

  ! local parameters
  integer :: i, j, ispec, iglob, ier
  double precision :: xx, zz
  real(kind=CUSTOM_REAL) :: rhol,mul,kappal
  character(len=MAX_STRING_LEN) :: outputname

  ! elastic kernels
  ! only at the last time step we multiply by deltat and parameter value, it is not necessary to do it at each iteration
  ! note: acoustic kernels add these at each iteration step already

  do ispec = 1, nspec
    if (ispec_is_elastic(ispec)) then
      ! isotropic kernels
      do j = 1, NGLLZ
        do i = 1, NGLLX
          ! gets elastic moduli
          !
          ! daniel todo - please check:
          ! note: in case of attenuation, mu/kappa store and poroelastcoef parameters have been scaled
          !       to unrelaxed moduli in routine prepare_attenuation()
          !
          !       this should probably be corrected to have again the moduli at the (initial) reference frequency
          !       for calculating the kernel values below...
          rhol = rhostore(i,j,ispec)
          mul = mustore(i,j,ispec)
          kappal = kappastore(i,j,ispec)

          ! for parameterization (rho,mu,kappa): "primary" kernels
          !
          ! kernel expressions below are for relative perturbations
          ! for example using dln(rho) = drho/rho (and not for absolute perturbations drho)
          !
          ! dChi = \int_V [ K_rho dln(rho) + K_mu dln(mu) + K_kappa dln(kappa) ] d^3x
          !
          ! and K_rho   = - rho   int_0^T [ s^adj * \partial_t^2 s ] dt
          !     K_mu    = - 2 mu  int_0^T [ D^adj : D ] dt               (where D = 1/2[grad(s) + grad(s)^T] - 1/3 div(s) I )
          !     K_kappa = - kappa int_0^T [ div(s^adj) div(s) ] dt
          !
          ! (see e.g. Peter et al. 2011, eq. 15, 18, 19)

          ! density kernel
          rho_kl(i,j,ispec) =  - rhol * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * rho_kl(i,j,ispec)
          ! shear modulus kernel
          mu_kl(i,j,ispec) =  - TWO * mul * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * mu_kl(i,j,ispec)
          ! bulk modulus kernel
          kappa_kl(i,j,ispec) = - kappal * (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * kappa_kl(i,j,ispec)

          ! derived from "primary" kernels above...
          !

          ! daniel todo - please check:
          ! note: the definition of vp depends on the modulus definition of kappa (case AXISYM).
          !       this in turn can change the kernels formula...
          !
          !if (AXISYM) then ! CHECK kappa
          !  vp = sqrt((kappal + FOUR_THIRDS * mul)/rhol)
          !else
          !  vp = sqrt((kappal + mul)/rhol)
          !endif
          !vs = sqrt(mul/rhol)

          ! for parameterization (rho,beta,alpha):
          ! rho prime kernel
          rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + kappa_kl(i,j,ispec) + mu_kl(i,j,ispec)
          ! Vs kernel
          ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) CHECK Kappa
          beta_kl(i,j,ispec) = TWO * (mu_kl(i,j,ispec) - FOUR_THIRDS * mul/kappal * kappa_kl(i,j,ispec))
          ! Vp kernel
          ! ABAB !! Warning !! This is possibly false for plane strain (look for: bulk modulus plane strain) Check Kappa
          alpha_kl(i,j,ispec) = TWO * (1._CUSTOM_REAL + FOUR_THIRDS * mul/kappal) * kappa_kl(i,j,ispec)

          ! for bulk velocity c parameterization (rho,bulk_c,beta):
          bulk_c_kl(i,j,ispec) =  TWO * kappa_kl(i,j,ispec)
          bulk_beta_kl(i,j,ispec) =  TWO * mu_kl(i,j,ispec)
        enddo
      enddo

      ! Voigt kernels, e.g., see Sieminski, 2007a,b
      if (ispec_is_anisotropic(ispec)) then
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            ! "primary" kernels
            c11_kl(i,j,ispec) = - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * c11_kl(i,j,ispec)
            c13_kl(i,j,ispec) = - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * c13_kl(i,j,ispec)
            c15_kl(i,j,ispec) = - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * c15_kl(i,j,ispec)
            c33_kl(i,j,ispec) = - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * c33_kl(i,j,ispec)
            c35_kl(i,j,ispec) = - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * c35_kl(i,j,ispec)
            c55_kl(i,j,ispec) = - (deltat * NTSTEP_BETWEEN_COMPUTE_KERNELS) * c55_kl(i,j,ispec)
            ! rho prime kernel
            rhop_kl(i,j,ispec) = rho_kl(i,j,ispec) + c11_kl(i,j,ispec) + &
                                 c13_kl(i,j,ispec) + c15_kl(i,j,ispec) + c33_kl(i,j,ispec) + &
                                 c35_kl(i,j,ispec) + c55_kl(i,j,ispec)
          enddo
        enddo
      endif

    endif ! elastic element
  enddo !nspec loop

  ! saves elastic kernels to file
  ! elastic domains
  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Elastic kernels:'
    write(IMAIN,*) '  maximum value of rho  kernel      = ',maxval(rho_kl)
    if (count(ispec_is_anisotropic(:) .eqv. .true.) >= 1) then
      write(IMAIN,*) '  maximum value of c11 kernel       = ',maxval(c11_kl)
      write(IMAIN,*) '  maximum value of c13 kernel       = ',maxval(c13_kl)
      write(IMAIN,*) '  maximum value of c15 kernel       = ',maxval(c15_kl)
      write(IMAIN,*) '  maximum value of c33 kernel       = ',maxval(c33_kl)
      write(IMAIN,*) '  maximum value of c35 kernel       = ',maxval(c35_kl)
      write(IMAIN,*) '  maximum value of c55 kernel       = ',maxval(c55_kl)
    else
      write(IMAIN,*) '  maximum value of kappa kernel     = ',maxval(kappa_kl)
      write(IMAIN,*) '  maximum value of mu kernel        = ',maxval(mu_kl)
      write(IMAIN,*)
      write(IMAIN,*) '  maximum value of rho prime kernel = ',maxval(rhop_kl)
      write(IMAIN,*) '  maximum value of alpha kernel     = ',maxval(alpha_kl)
      write(IMAIN,*) '  maximum value of beta kernel      = ',maxval(beta_kl)
    endif
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! saves to files
  if (save_ASCII_kernels) then
    ! ascii format
    if (count(ispec_is_anisotropic(:) .eqv. .true.) >= 1) then
      ! anisotropic
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_cijkl_kernel.dat'
      open(unit = 97, file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            xx = coord(1,iglob)
            zz = coord(2,iglob)
            ! anisotropic kernel values
            write(97,'(9e15.5e4)') xx, zz, rho_kl(i,j,ispec), c11_kl(i,j,ispec), &
                                   c13_kl(i,j,ispec), c15_kl(i,j,ispec), c33_kl(i,j,ispec), c35_kl(i,j,ispec), &
                                   c55_kl(i,j,ispec)
          enddo
        enddo
      enddo

      close(97)

    else
      ! isotropic kernels
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_mu_kernel.dat'
      open(unit = 97, file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_alpha_beta_kernel.dat'
      open(unit = 98, file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            xx = coord(1,iglob)
            zz = coord(2,iglob)
            ! isotropic kernel values
            ! parameterization (rho,kappa,mu) "primary" kernels
            write(97,'(5e15.5e4)') xx,zz,rho_kl(i,j,ispec),kappa_kl(i,j,ispec),mu_kl(i,j,ispec)
            ! parameterization (rho,Vp,Vs)
            write(98,'(5e15.5e4)') xx,zz,rhop_kl(i,j,ispec),alpha_kl(i,j,ispec),beta_kl(i,j,ispec)
          enddo
        enddo
      enddo

      close(97)
      close(98)
    endif

    ! Hessian kernels
    if (APPROXIMATE_HESS_KL) then
      ! isotropic kernels
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian1_kernel.dat'
      open(unit = 97, file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian2_kernel.dat'
      open(unit = 98, file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            xx = coord(1,iglob)
            zz = coord(2,iglob)
            write(97,'(3e15.5e4)') xx,zz,rhorho_el_Hessian_final1(i,j,ispec)
            write(98,'(3e15.5e4)') xx,zz,rhorho_el_Hessian_final2(i,j,ispec)
          enddo
        enddo
      enddo

      close(97)
      close(98)
    endif

  else
    ! binary format
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kernel.bin'
    open(unit = 204, file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
    write(204) rho_kl
    close(204)

    if ((count(ispec_is_anisotropic(:)) >= 1) .eqv. .true.) then
      ! anisotropic
      write(outputname,'(a,i6.6,a)')'proc',myrank,'_c11_kernel.bin'
      open(unit = 205,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c13_kernel.bin'
      open(unit = 206,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c15_kernel.bin'
      open(unit = 207,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c33_kernel.bin'
      open(unit = 208,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c35_kernel.bin'
      open(unit = 209,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c55_kernel.bin'
      open(unit = 210,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

      write(205) c11_kl
      write(206) c13_kl
      write(207) c15_kl
      write(208) c33_kl
      write(209) c35_kl
      write(210) c55_kl

      close(205)
      close(206)
      close(207)
      close(208)
      close(209)
      close(210)

    else
      ! isotropic
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_kappa_kernel.bin'
      open(unit = 205, file =trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mu_kernel.bin'
      open(unit = 206, file =trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_kernel.bin'
      open(unit = 207, file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_alpha_kernel.bin'
      open(unit = 208, file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_beta_kernel.bin'
      open(unit = 209, file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted', iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_bulk_c_kernel.bin'
      open(unit = 210,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_bulk_beta_kernel.bin'
      open(unit = 211,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

      write(205) kappa_kl
      write(206) mu_kl
      write(207) rhop_kl
      write(208) alpha_kl
      write(209) beta_kl
      write(210) bulk_c_kl
      write(211) bulk_beta_kl

      close(205)
      close(206)
      close(207)
      close(208)
      close(209)
      close(210)
      close(211)
    endif

    ! Hessian kernels
    if (APPROXIMATE_HESS_KL) then
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian1_kernel.bin'
      open(unit = 214,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian2_kernel.bin'
      open(unit = 215,file=trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

      write(214) rhorho_el_Hessian_final1
      write(215) rhorho_el_Hessian_final2

      close(214)
      close(215)
    endif
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine save_kernels_elastic

!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_acoustic()

  use constants, only: NGLLX,NGLLZ,IMAIN,CUSTOM_REAL,MAX_STRING_LEN,OUTPUT_FILES

  use specfem_par, only: myrank, nspec, ibool, coord, save_ASCII_kernels, &
                         rho_ac_kl, kappa_ac_kl, alpha_ac_kl, rhop_ac_kl, &
                         rhorho_ac_Hessian_final1, rhorho_ac_Hessian_final2, &
                         APPROXIMATE_HESS_KL

  implicit none

  ! local parameters
  integer :: i, j, ispec, iglob, ier
  double precision :: xx, zz
  character(len=MAX_STRING_LEN) :: outputname

  ! acoustic domain
  ! note: acoustic kernels have already added deltat and parameter values, thus only file output needed
  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Acoustic kernels:'
    write(IMAIN,*) '  maximum value of rho kernel       = ',maxval(rho_ac_kl)
    write(IMAIN,*) '  maximum value of kappa kernel     = ',maxval(kappa_ac_kl)
    write(IMAIN,*)
    write(IMAIN,*) '  maximum value of rho prime kernel = ',maxval(rhop_ac_kl)
    write(IMAIN,*) '  maximum value of alpha kernel     = ',maxval(alpha_ac_kl)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! saves to file
  if (save_ASCII_kernels) then
    ! ascii format
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_kernel.dat'
    open(unit = 95, file = trim(OUTPUT_FILES)//trim(outputname),status ='unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_c_kernel.dat'
    open(unit = 96, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

    do ispec = 1, nspec
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          xx = coord(1,iglob)
          zz = coord(2,iglob)
          write(95,'(4e15.5e4)') xx,zz,rho_ac_kl(i,j,ispec),kappa_ac_kl(i,j,ispec)
          write(96,'(4e15.5e4)') xx,zz,rhop_ac_kl(i,j,ispec),alpha_ac_kl(i,j,ispec)
        enddo
      enddo
    enddo

    close(95)
    close(96)

    ! Hessian kernels
    if (APPROXIMATE_HESS_KL) then
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian1_acoustic_kernel.dat'
      open(unit = 95, file = trim(OUTPUT_FILES)//trim(outputname),status ='unknown',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian2_acoustic_kernel.dat'
      open(unit = 96, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            xx = coord(1,iglob)
            zz = coord(2,iglob)
            write(95,'(3e15.5e4)') xx,zz,rhorho_ac_Hessian_final1(i,j,ispec)
            write(96,'(3e15.5e4)') xx,zz,rhorho_ac_Hessian_final2(i,j,ispec)
          enddo
        enddo
      enddo

      close(95)
      close(96)
    endif

  else
    ! binary format
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_acoustic_kernel.bin'
    open(unit = 200, file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_kappa_acoustic_kernel.bin'
    open(unit = 201, file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_acoustic_kernel.bin'
    open(unit = 202, file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_c_acoustic_kernel.bin'
    open(unit = 203, file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted', iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
    if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

    write(200) rho_ac_kl
    write(201) kappa_ac_kl
    write(202) rhop_ac_kl
    write(203) alpha_ac_kl

    close(200)
    close(201)
    close(202)
    close(203)

    ! Hessian kernels
    if (APPROXIMATE_HESS_KL) then
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian1_acoustic_kernel.bin'
      open(unit=212,file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Hessian2_acoustic_kernel.bin'
      open(unit=213,file = trim(OUTPUT_FILES)//trim(outputname),status='unknown',action='write',form='unformatted',iostat=ier)
      if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
      if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
      write(212) rhorho_ac_Hessian_final1
      write(213) rhorho_ac_Hessian_final2
      close(212)
      close(213)
    endif
  endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine save_kernels_acoustic


!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_kernels_poroelastic()

  use constants, only: NGLLX,NGLLZ,IMAIN,CUSTOM_REAL,MAX_STRING_LEN,OUTPUT_FILES

  use specfem_par, only: myrank, nspec, ibool, coord, save_ASCII_kernels, &
                         rhot_kl, rhof_kl, sm_kl, eta_kl, mufr_kl, B_kl, &
                         C_kl, M_kl, rhob_kl, rhofb_kl, phi_kl, mufrb_kl, &
                         rhobb_kl, rhofbb_kl, phib_kl, cpI_kl, cpII_kl, cs_kl, ratio_kl, &
                         GPU_MODE

  implicit none

  ! local parameters
  integer :: i, j, ispec, iglob, ier
  double precision :: xx, zz
  character(len=MAX_STRING_LEN) :: outputname

  ! saves poroelastic kernels

  ! poro-elastic domains
  ! safety checks
  if (GPU_MODE) call stop_the_code('poroelastic kernel output not implemented on GPUs yet')
  if (.not. SAVE_ASCII_KERNELS) call stop_the_code('poroelastic simulations must use SAVE_ASCII_KERNELS')

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) 'Poroelastic kernels:'
    write(IMAIN,*) '  maximum value of rho kernel   = ',maxval(rhot_kl)
    write(IMAIN,*) '  maximum value of rho_f kernel = ',maxval(rhof_kl)
    write(IMAIN,*) '  maximum value of m kernel     = ',maxval(sm_kl)
    write(IMAIN,*) '  maximum value of M kernel     = ',maxval(M_kl)
    write(IMAIN,*) '  maximum value of B kernel     = ',maxval(B_kl)
    write(IMAIN,*) '  maximum value of C kernel     = ',maxval(C_kl)
    write(IMAIN,*) '  maximum value of mu_fr kernel = ',maxval(mufr_kl)
    write(IMAIN,*)
    write(IMAIN,*) '  maximum value of cpI kernel   = ',maxval(cpI_kl)
    write(IMAIN,*) '  maximum value of cpII kernel  = ',maxval(cpII_kl)
    write(IMAIN,*) '  maximum value of cs kernel    = ',maxval(cs_kl)
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  ! Primary kernels
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mu_B_C_kernel.dat'
  open(unit = 144, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
  if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
  if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_M_rho_rhof_kernel.dat'
  open(unit = 155, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
  if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
  if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_m_eta_kernel.dat'
  open(unit = 16, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
  if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
  if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
  ! Wavespeed kernels
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_cpI_cpII_cs_kernel.dat'
  open(unit = 20, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
  if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
  if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhobb_rhofbb_ratio_kernel.dat'
  open(unit = 21, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
  if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
  if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_phib_eta_kernel.dat'
  open(unit = 22, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
  if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
  if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
  ! Density normalized kernels
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mub_Bb_Cb_kernel.dat'
  open(unit = 17, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
  if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
  if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Mb_rhob_rhofb_kernel.dat'
  open(unit = 18, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
  if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
  if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)
  write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mb_etab_kernel.dat'
  open(unit = 19, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
  if (ier /= 0) call stop_the_code('Error writing kernel file to disk')
  if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

  ! ascii format
  do ispec = 1, nspec
    do j = 1, NGLLZ
      do i = 1, NGLLX
        iglob = ibool(i,j,ispec)
        xx = coord(1,iglob)
        zz = coord(2,iglob)
        ! primary kernels
        write(144,'(5e15.5e4)') xx,zz,mufr_kl(i,j,ispec),B_kl(i,j,ispec),C_kl(i,j,ispec)
        write(155,'(5e15.5e4)') xx,zz,M_kl(i,j,ispec),rhot_kl(i,j,ispec),rhof_kl(i,j,ispec)
        write(16,'(5e15.5e4)') xx,zz,sm_kl(i,j,ispec),eta_kl(i,j,ispec)
        ! density normalized kernels
        write(17,'(5e15.5e4)') xx,zz,mufrb_kl(i,j,ispec),B_kl(i,j,ispec),C_kl(i,j,ispec)
        write(18,'(5e15.5e4)') xx,zz,M_kl(i,j,ispec),rhob_kl(i,j,ispec),rhofb_kl(i,j,ispec)
        write(19,'(5e15.5e4)') xx,zz,phi_kl(i,j,ispec),eta_kl(i,j,ispec)
        ! wavespeed kernels
        write(20,'(5e15.5e4)') xx,zz,cpI_kl(i,j,ispec),cpII_kl(i,j,ispec),cs_kl(i,j,ispec)
        write(21,'(5e15.5e4)') xx,zz,rhobb_kl(i,j,ispec),rhofbb_kl(i,j,ispec),ratio_kl(i,j,ispec)
        write(22,'(5e15.5e4)') xx,zz,phib_kl(i,j,ispec),eta_kl(i,j,ispec)
      enddo
    enddo
  enddo

  close(144)
  close(155)
  close(16)
  close(17)
  close(18)
  close(19)
  close(20)
  close(21)
  close(22)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine save_kernels_poroelastic


!
!-------------------------------------------------------------------------------------------------
!

  subroutine save_weights_kernel()

!> Save weights for volume integration,
!! in order to benchmark the kernels with analytical expressions.

  use constants, only: NGLLX,NGLLZ,IMAIN,CUSTOM_REAL,MAX_STRING_LEN,OUTPUT_FILES

  use specfem_par, only: myrank, nspec, ibool, coord, save_ASCII_kernels, &
                         jacobian,wxgll,wzgll

  implicit none

  ! local parameters
  integer :: i, j, ispec, iglob, ier
  double precision :: xx, zz
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: weights_kernel
  real(kind=CUSTOM_REAL) :: jacobianl
  character(len=MAX_STRING_LEN) :: outputname

  allocate(weights_kernel(NGLLX,NGLLZ,nspec),stat=ier)
  if (ier /= 0) call stop_the_code('Error allocating weights')
  weights_kernel(:,:,:) = 0.0_CUSTOM_REAL

  do ispec = 1,nspec
    do j = 1, NGLLZ
      do i = 1, NGLLX
        jacobianl = jacobian(i,j,ispec)
        weights_kernel(i,j,ispec) = wxgll(i) * wzgll(j) * jacobianl
      enddo
    enddo
  enddo

  if (save_ASCII_kernels) then
    ! ascii format
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_weights_kernel.dat'
    open(unit = 144, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing weights file to disk')
    if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

    do ispec = 1, nspec
      do j = 1, NGLLZ
        do i = 1, NGLLX
          iglob = ibool(i,j,ispec)
          xx = coord(1,iglob)
          zz = coord(2,iglob)
          ! format: #x  #z  #weight
          write(144,'(5e15.5e4)') xx,zz,weights_kernel(i,j,ispec)
        enddo
      enddo
    enddo

    close(144)

  else
    ! binary format
    write(outputname,'(a,i6.6,a)') 'proc',myrank,'_weights_kernel.bin'
    open(unit = 144, file = trim(OUTPUT_FILES)//trim(outputname),status = 'unknown',form='unformatted',iostat=ier)
    if (ier /= 0) call stop_the_code('Error writing weights file to disk')
    if (myrank == 0) write(IMAIN,*) '  ',trim(OUTPUT_FILES)//trim(outputname)

    write(144) weights_kernel

    close(144)
  endif

  deallocate(weights_kernel)

  ! user output
  if (myrank == 0) then
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  end subroutine save_weights_kernel


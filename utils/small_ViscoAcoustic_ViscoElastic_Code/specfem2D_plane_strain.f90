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


! 21 juin 2016 introduction att calcul des varaibles tau_eps tau_sigma ...
! reste le calcul des variables a memoire et de la matrice de rigidite modifiee des variables a memoire
! sauvegarde avant de depoter le calcul des int de regidite dans un fichier a part

program serial_specfem2D

  use small_specfem_par

  implicit none

  logical ::  is_viscoelastic, is_viscoacoustic


!!!!!!!!!!
!!!!!!!!!! NGLLX and NGLLZ are set equal to 5,
!!!!!!!!!! therefore each element contains NGLLX * NGLLZ = 25 points.
!!!!!!!!!!

!!!!!!!!!!
!!!!!!!!!! All the calculations are done in single precision.
!!!!!!!!!! We do not need double precision in SPECFEM2D.
!!!!!!!!!!






  ! allocate Lagrange interpolators for receivers
  allocate(ispec_selected_rec(nrec),xi_receiver(nrec),gamma_receiver(nrec),st_xval(nrec),st_zval(nrec))
  allocate(x_final_receiver(nrec),z_final_receiver(nrec),seisNameX(nrec),seisNameZ(nrec))
  allocate(hxir_store(nrec,NGLLX),hgammar_store(nrec,NGLLZ),seismogramX(NSTEP,nrec),seismogramZ(NSTEP,nrec))
  allocate(hxir(NGLLX),hgammar(NGLLZ))
  allocate(hpxir(NGLLX),hpgammar(NGLLZ))
  allocate(forced(NGLOB))
  forced(:) = .false.

  seisNameX(1) = "X1.txt" ;  seisNameX(2) = "X2.txt" ;  seisNameX(3) = "X3.txt";  seisNameX(4) = "X4.txt";  seisNameX(5) = "X5.txt"
  seisNameX(6) = "X6.txt";  seisNameX(7) = "X7.txt";  seisNameZ(1) = "Z1.txt";  seisNameZ(2) = "Z2.txt";  seisNameZ(3) = "Z3.txt"
  seisNameZ(4) = "Z4.txt";  seisNameZ(5) = "Z5.txt";  seisNameZ(6) = "Z6.txt";  seisNameZ(7) = "Z7.txt";


  print *
  print *,'NSPEC = ',NSPEC
  print *,'NGLOB = ',NGLOB
  print *

  print *,'NSTEP = ',NSTEP
  print *,'deltat = ',deltat
  print *



  ! set up Gauss-Lobatto-Legendre points, weights and also derivation matrices

    call define_derivation_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz,NGLLX,NGLLZ)






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !--- definition of the mesh
  do iz = 0,nz
     do ix = 0,nx
        ! coordinates of the grid points (evenly spaced points along X and Z)
        xgrid(ix,iz) = xmin + (xmax - xmin) * dble(ix) / dble(nx)
        zgrid(ix,iz) = zmin + (zmax - zmin) * dble(iz) / dble(nz)
     enddo
  enddo

  ! create the coorg array
  do j = 0,nz
     do i = 0,nx
        ipoin = num(i,j,nx)
        coorg(1,ipoin) = xgrid(i,j)
        coorg(2,ipoin) = zgrid(i,j)
     enddo
  enddo

  ! create the knods array
  k = 0
  do j=0,nz-1
     do i=0,nx-1
        k = k + 1
        knods(1,k) = num(i,j,nx)
        knods(2,k) = num(i+1,j,nx)
        knods(3,k) = num(i+1,j+1,nx)
        knods(4,k) = num(i,j+1,nx)
     enddo
  enddo

  !
  !---- generate the global numbering
  !

  call createnum_slow(knods,ibool,nglob_to_compute,nspec,NGLLX,NGLLZ,ngnod) ! Create ibool and recompute nglob for checking
  if (nglob_to_compute /= NGLOB) stop 'error: incorrect total number of unique grid points found'
  if (minval(ibool) /= 1) stop 'error: incorrect minimum value of ibool'
  if (maxval(ibool) /= NGLOB) stop 'error: incorrect maximum value of ibool'

  !
  !----  set the coordinates of the points of the global grid
  !
  do ispec = 1,nspec
     do j = 1,NGLLZ
        do i = 1,NGLLX
           xi = xigll(i)
           gamma = zigll(j)

           call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl, &
                jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo,NDIM)
           if (jacobianl <= 0.d0) stop 'error: negative Jacobian found'

           coord(1,ibool(i,j,ispec)) = x
           coord(2,ibool(i,j,ispec)) = z

           xix(i,j,ispec) = xixl
           xiz(i,j,ispec) = xizl
           gammax(i,j,ispec) = gammaxl
           gammaz(i,j,ispec) = gammazl
           jacobian(i,j,ispec) = jacobianl

        enddo
     enddo
  enddo



    ! must be call here  because we need to compute the displacmetn gradient and we need : xix, xiz, gammax, gammaz, jacobian
    call  check_attenuation(is_viscoelastic,is_viscoacoustic) !check the type of simulation acoustic, viscoacoustic, elastic or viscoelastic ?
    if (is_viscoelastic .or. is_viscoacoustic) then
     call prepare_timerun_attenuation(is_viscoelastic,is_viscoacoustic)
    endif


  call build_forced(xforced,coord,NGLOB,forced,nb_forced)
  print *,'nb_forced', nb_forced

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! build the mass matrix
  rmass_inverse_elastic(:) = 0.
  rmass_inverse_acoustic(:) = 0.

  if (.not. is_acoustic) then
     do ispec = 1,NSPEC
        do j = 1,NGLLZ
           do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              rmass_inverse_elastic(iglob) = rmass_inverse_elastic(iglob) + wxgll(i)*wzgll(j)*rho*jacobian(i,j,ispec) !equivalent lines in the big code at lines 250-251 of Invert_mass_matrix
           enddo
        enddo
     enddo
  else
     do ispec = 1,NSPEC
        do j = 1,NGLLZ
           do i = 1,NGLLX
              iglob = ibool(i,j,ispec)
              rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                   + wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / kappal
           enddo
        enddo
     enddo
  endif



  ! we have built the real exactly diagonal mass matrix (not its inverse)
  ! therefore invert it here once and for all

  if (.not. is_acoustic) then
     do i = 1,NGLOB
        rmass_inverse_elastic(i) = 1. / rmass_inverse_elastic(i)
     enddo
  else
     do i = 1,NGLOB
        rmass_inverse_acoustic(i) = 1. / rmass_inverse_acoustic(i)
     enddo
  endif



  ! compute the position of the source and of the receiver
  x_source = coord(1,ibool(IGLL_SOURCE,JGLL_SOURCE,NSPEC_SOURCE))
  z_source = coord(2,ibool(IGLL_SOURCE,JGLL_SOURCE,NSPEC_SOURCE))
  x_receiver = coord(1,ibool(IGLL_RECEIVER,JGLL_RECEIVER,NSPEC_RECEIVER))
  z_receiver = coord(2,ibool(IGLL_RECEIVER,JGLL_RECEIVER,NSPEC_RECEIVER))

  call locate_receivers(ibool,coord,nspec,NGLOB,xigll,zigll, &
       nrec,st_xval,st_zval,ispec_selected_rec, &
       xi_receiver,gamma_receiver,x_source,z_source, &
       coorg,knods,ngnod,npgeo, &
       x_final_receiver, z_final_receiver,NDIM,NGLLX,NGLLZ)

  ! clear initial vectors before starting the time loop
  ! init elastic variables
  displ(:,:) = 0. !!!!!!!!!! VERYSMALLVAL
  veloc(:,:) = 0.
  accel(:,:) = 0.
  ! init acoutic variables
  potential_acoustic(:)= 0.
  potential_dot_acoustic(:)= 0.
  potential_dot_dot_acoustic(:)= 0.
  acoustic_displ(:,:) = 0.

  !  TODO a display subruntine
  ! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
  ! time_values(3): day of the month
  ! time_values(5): hour of the day
  ! time_values(6): minutes of the hour
  ! time_values(7): seconds of the minute
  ! time_values(8): milliseconds of the second
  ! this fails if we cross the end of the month
  time_start = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
       60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

  ! start of the time loop (which must remain serial obviously)
  do it = 1,NSTEP

     ! compute maximum of norm of displacement from time to time and display it
     ! in order to monitor the simulation
     if (mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 2 .or. it == NSTEP) then
        Usolidnorm = -1.
        do iglob = 1,NGLOB
           current_value = sqrt(displ(1,iglob)**2 + displ(2,iglob)**2)
           if (current_value > Usolidnorm) Usolidnorm = current_value
        enddo
        write(*,*) 'Time step # ',it,' out of ',NSTEP
        ! compute current time
        time = (it-1)*deltat
        write(*,*) 'Max norm displacement vector U in the solid (m) = ',Usolidnorm
        ! check stability of the code, exit if unstable
        if (Usolidnorm > STABILITY_THRESHOLD .or. Usolidnorm < 0) stop 'code became unstable and blew up'
        ! count elapsed wall-clock time
        call date_and_time(datein,timein,zone,time_values)
        time_end = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
             60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

        ! elapsed time since beginning of the simulation
        tCPU = time_end - time_start
        int_tCPU = int(tCPU)
        ihours = int_tCPU / 3600
        iminutes = (int_tCPU - 3600*ihours) / 60
        iseconds = int_tCPU - 3600*ihours - 60*iminutes
        write(*,*) 'Elapsed time in seconds = ',tCPU
        write(*,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
        write(*,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
        write(*,*)

        if (.not. is_acoustic) then !elastic
           ! draw the displacement vector field in a PostScript file
           call plot_post(displ,coord,ibool,NGLOB,NSPEC,x_source,z_source,st_xval,st_zval,it,deltat,NGLLX,NGLLZ,NDIM,nrec)
        else
           call plot_post(acoustic_displ, &
                coord,ibool,NGLOB,NSPEC,x_source,z_source,st_xval,st_zval,it,deltat,NGLLX,NGLLZ,NDIM,nrec)
        endif

     endif


     if (.not. is_acoustic) then !elastic
        do iglob = 1,NGLOB
           if (forced(iglob)) then
              call enforce_elastic_forcing(NGLOB,iglob,it,deltat,deltatover2,deltatsqover2,coord,displ,veloc,accel)
           else
              ! big loop over all the global points (not elements) in the mesh to update
              ! the displacement and velocity vectors and clear the acceleration vector
              displ(:,iglob) = displ(:,iglob) + deltat*veloc(:,iglob) + deltatsqover2*accel(:,iglob)
              veloc(:,iglob) = veloc(:,iglob) + deltatover2*accel(:,iglob) ! Predictor
              accel(:,iglob) = 0.
           endif
        enddo
     else !acoustic
        do iglob = 1,NGLOB
           if (forced(iglob)) then !acoustic and forced
              call enforce_acoustic_forcing(NGLOB,iglob,it,deltat,deltatover2,deltatsqover2,coord, &
                   potential_acoustic,potential_dot_acoustic,potential_dot_dot_acoustic)
           else !acoustic not forced
              !time sheme for acoustic media
              potential_acoustic(iglob) = potential_acoustic(iglob) + deltat * potential_dot_acoustic(iglob) &
                   + deltatsqover2 * potential_dot_dot_acoustic(iglob)
              potential_dot_acoustic(iglob) = potential_dot_acoustic(iglob) + deltatover2 * potential_dot_dot_acoustic(iglob) ! Predictor @ (it-1/2)*deltat
              potential_dot_dot_acoustic(iglob) = 0.

!              if (iglob==3000) then
!              print *,'********************************************************************************'
!              print *,'avant calcul force it = ' ,it , "potential_acoustic(3000)=  ", potential_acoustic(3000)
!              endif
           endif !end acoustic
        enddo
     endif


     ! big loop over all the elements in the mesh to localize data
     ! from the global vectors to the local mesh
     ! using indirect addressing (contained in array ibool)
     ! and then to compute the elemental contribution
     ! to the acceleration vector of each element of the finite-element mesh

     ! Start compute forces
     if (.not. is_acoustic) then !elastic
        ! get elastic parameters of current spectral element pas dasn la boucle car constant
        mu = rho*cs*cs
        lambda = rho*cp*cp - 2.*mu
        lambdaplus2mu = lambda + 2.*mu
        lambdaplusmu = lambda + mu

        if (is_viscoelastic) call update_memory_variable(is_viscoelastic,is_viscoacoustic)
        call compute_viscoelastic_forces(is_viscoelastic)

     else !end elastic then it an acoustic media

        if (is_viscoacoustic) call update_memory_variable(is_viscoelastic,is_viscoacoustic)
        call compute_viscoacoustic_forces(is_viscoacoustic)

     endif !end compute forces


     ! add the force source at a given grid point
     iglob_source = ibool(IGLL_SOURCE,JGLL_SOURCE,NSPEC_SOURCE) ! determine the global (grid) point at which the source is located
     time = (it-1)*deltat ! compute current time
     if (.not. is_acoustic) then !elastic
        accel(2,iglob_source) = accel(2,iglob_source) - factor_amplitude * (1.-2.*a*(time-t0)**2) * exp(-a*(time-t0)**2)
     else !acoustic
        potential_dot_dot_acoustic(iglob_source) = potential_dot_dot_acoustic(iglob_source) - &
             factor_amplitude * (1.-2.*a*(time-t0)**2) * exp(-a*(time-t0)**2)/kappal

        ! ligne 158-160 prepare_source_time_function
        !aval(i_source) = PI*PI*f0_source(i_source)*f0_source(i_source) ! set_source_parameters.f90:98:
        !source_time_function(i_source,it,i_stage) = - factor(i_source) * &
        !                        (ONE-TWO*aval(i_source)*t_used**2) * &
        !                        exp(-aval(i_source)*t_used**2)
        ! potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - & !lignes83-84 compute_add_sources_acoustic.f90
        !                                                 sourcearrays(i_source,1,i,j) * stf_used / kappastore(i,j,ispec)
     endif


     ! big loop over all the global points (not elements) in the mesh to update
     ! the acceleration and velocity vectors.
     ! To compute acceleration from the elastic forces we need to divide them by the mass matrix, i.e. multiply by its inverse
     if (.not. is_acoustic) then !elastic
        do iglob = 1,NGLOB
           ! big loop over all the global points (not elements) in the mesh to update
           ! the acceleration and velocity vectors.
           ! To compute acceleration from the elastic forces we need to divide them by the mass matrix, i.e. multiply by its inverse
           if (.not. forced(iglob)) then
              accel(1,iglob) = accel(1,iglob)*rmass_inverse_elastic(iglob)
              accel(2,iglob) = accel(2,iglob)*rmass_inverse_elastic(iglob)
              veloc(:,iglob) = veloc(:,iglob) + deltatover2*accel(:,iglob) !corrector
           endif
        enddo
     else ! acoustic and not forced
        do iglob = 1,NGLOB
           if (.not. forced(iglob)) then
              potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) * rmass_inverse_acoustic(iglob)
              potential_dot_acoustic(iglob) = potential_dot_acoustic(iglob) + deltatover2 * potential_dot_dot_acoustic(iglob) !corrector


           endif
        enddo
     endif


     if (.not. is_acoustic) then !elastic
        ! record a seismogram to check that the simulation went well
        do irec = 1,nrec
           ! initializes interpolated values
           dx = 0.d0
           dz = 0.d0
           ispec = ispec_selected_rec(irec)
           call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
           call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
           hxir_store(irec,:) = hxir(:)
           hgammar_store(irec,:) = hgammar(:)
           do j = 1,NGLLZ
              do i = 1,NGLLX
                 iglob = ibool(i,j,ispec)
                 hlagrange = hxir_store(irec,i)*hgammar_store(irec,j)
                 ! computes interpolated field
                 dx = dx + displ(1,iglob)*hlagrange
                 dz = dz + displ(2,iglob)*hlagrange
              enddo
           enddo
           seismogramX(it,irec) = dx
           seismogramZ(it,irec) = dz
        enddo
     else
        do irec = 1,nrec
           ! initializes interpolated values
           dx = 0.d0
           dz = 0.d0
           ispec = ispec_selected_rec(irec)
           call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
           call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
           hxir_store(irec,:) = hxir(:)
           hgammar_store(irec,:) = hgammar(:)
           do j = 1,NGLLZ
              do i = 1,NGLLX
                 iglob = ibool(i,j,ispec)
                 hlagrange = hxir_store(irec,i)*hgammar_store(irec,j)
                 ! computes interpolated field
                 dx = dx + acoustic_displ(1,iglob)*hlagrange
                 dz = dz + acoustic_displ(2,iglob)*hlagrange
              enddo
           enddo
           seismogramX(it,irec) = dx
           seismogramZ(it,irec) = dz
        enddo
     endif



     !seismogramX(it,1) = displ(1,50400)
     !seismogramX(it,2) = veloc(1,50400)
     !seismogramX(it,3) = accel(1,50400)


  enddo ! end of the serial time loop

  ! save the seismogram at the end of the run
  do irec = 1,nrec
     open(unit=IIN,file=seisNameX(irec),status='unknown')
     do it = 1,NSTEP
        write(IIN,*) (it-1)*deltat,seismogramX(it,irec)
     enddo
     close(IIN)
     open(unit=IIN,file=seisNameZ(irec),status='unknown')
     do it = 1,NSTEP
        write(IIN,*) (it-1)*deltat,seismogramZ(it,irec)
     enddo
     close(IIN)
  enddo


contains




    subroutine check_attenuation(is_viscoelastic,is_viscoacoustic)
    ! subroutine check_attenuation()
    ! check that attenuation values entered by the user make sense
    ! if is_acoustic=.true.  an Qkappa > 9999 then it is an acoustic simulation
    ! if is_acoustic=.true.  an Qkappa < 9999 then it is an viscoacoustic simulation
    ! if is_acoustic=.false. an Qkappa and Qmu are >9999 then is is an elastic simulation
    ! if is_acoustic=.false. an Qkappa and Qmu are both < 9999 then is is an viscoelastic simulation


    use small_specfem_par
    implicit none
    logical, intent(inout) :: is_viscoelastic,is_viscoacoustic

    !imput variable : is_acoustic
    !logical, intent(out)  :: is_viscoelastic, is_viscoacoustic

    if (.not. is_acoustic) then

       if (QKappa > 9998.999d0 .or. Qmu > 9998.999d0) then

          if ((QKappa <= 9998.999d0 .and. Qmu > 9998.999d0) &
               .or. (QKappa > 9998.999d0 .and. Qmu <= 9998.999d0)) then
             print *, 'QKappa and Qmu the both must be above or below 9999'
          endif

          is_viscoelastic=.false.
          print *,'Elastic simulation without attenuation&
               & : is_viscoelastic= ', is_viscoelastic
       else
          is_viscoelastic=.true.
          print *,'Elastic simulation with attenuation &
               &: is_viscoelastic= ' ,  is_viscoelastic
       endif

    else
       if (QKappa > 9998.999d0) then
          is_viscoacoustic=.false.
          print *,'Acoustic simulation without attenuation : &
               &is_viscoacoustic= ' ,is_viscoacoustic
       else
          is_viscoacoustic=.true.
          print *,'Acoustic simulation with attenuation &
               &: is_viscoacoustic= ', is_viscoacoustic
       endif
    endif

  end subroutine check_attenuation



  ! --------------------- --------------------- --------------------- ---------------------

  subroutine prepare_timerun_attenuation(is_viscoelastic,is_viscoacoustic )

    use small_specfem_par
    implicit none
    logical, intent(in) :: is_viscoelastic,is_viscoacoustic


    ! allocate memory variables for attenuation
    allocate(e1(NGLLX,NGLLZ,NSPEC,N_SLS), &
         e11(NGLLX,NGLLZ,NSPEC,N_SLS), &
         e13(NGLLX,NGLLZ,NSPEC,N_SLS))

    allocate(inv_tau_sigma_nu1(NGLLX,NGLLZ,NSPEC,N_SLS))
    allocate(inv_tau_sigma_nu2(NGLLX,NGLLZ,NSPEC,N_SLS))
    allocate(phi_nu1(NGLLX,NGLLZ,NSPEC,N_SLS))
    allocate(phi_nu2(NGLLX,NGLLZ,NSPEC,N_SLS))
    allocate(Mu_nu1(NGLLX,NGLLZ,NSPEC))
    allocate(Mu_nu2(NGLLX,NGLLZ,NSPEC))
    allocate(tau_epsilon_nu1(N_SLS))
    allocate(tau_epsilon_nu2(N_SLS))
    allocate(inv_tau_sigma_nu1_sent(N_SLS))
    allocate(inv_tau_sigma_nu2_sent(N_SLS))
    allocate(phi_nu1_sent(N_SLS))
    allocate(phi_nu2_sent(N_SLS))


    e1(:,:,:,:) = 0.0
    e11(:,:,:,:) = 0.0
    e13(:,:,:,:) = 0.0

    if (is_viscoelastic) then

        ! initialize to dummy values
        ! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
        inv_tau_sigma_nu1(:,:,:,:) = -1.0
        phi_nu1(:,:,:,:) = -1.0
        inv_tau_sigma_nu2(:,:,:,:) = -1.0
        phi_nu2(:,:,:,:) = -1.0
        Mu_nu1(:,:,:) = -1.0
        Mu_nu2(:,:,:) = -1.0

        do ispec = 1,nspec

            call viscoelastic_attenuation_model()
           ! si acosutique il faut appeler une autre fonction
           !   call attenuation_viscoacoustic_model ! todo

           do j = 1,NGLLZ
              do i = 1,NGLLX
                 inv_tau_sigma_nu1(i,j,ispec,:) = inv_tau_sigma_nu1_sent(:)
                 phi_nu1(i,j,ispec,:) = phi_nu1_sent(:)
                 inv_tau_sigma_nu2(i,j,ispec,:) = inv_tau_sigma_nu2_sent(:) !0 si acoustique
                 phi_nu2(i,j,ispec,:) = phi_nu2_sent(:)!0 si acoustique
                 Mu_nu1(i,j,ispec) = Mu_nu1_sent
                 Mu_nu2(i,j,ispec) = Mu_nu2_sent !0 si acoustique
    !             if (mod(ispec,2000)==0) then
    !                print *,'ipsec= ', ispec, 'Qmu= ',  Qmu
    !             endif
              enddo
           enddo
        enddo
    endif

    if (is_viscoacoustic) then
        print *,"prepare time_run ViscoAcoustic !"

    ! initialize to dummy values
        ! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
        inv_tau_sigma_nu1(:,:,:,:) = -1.0
        phi_nu1(:,:,:,:) = -1.0
!       inv_tau_sigma_nu2(:,:,:,:) = -1.0 ! not used in this case but variables are allocate
!       phi_nu2(:,:,:,:) = -1.0           ! not used in this case but variables are allocate
        Mu_nu1(:,:,:) = -1.0
!       Mu_nu2(:,:,:) = -1.0

        do ispec = 1,nspec

            call viscoelastic_attenuation_model()
           ! si acosutique il faut appeler une autre fonction
           !   call attenuation_viscoacoustic_model ! todo

           do j = 1,NGLLZ
              do i = 1,NGLLX
                 inv_tau_sigma_nu1(i,j,ispec,:) = inv_tau_sigma_nu1_sent(:)
                 phi_nu1(i,j,ispec,:) = phi_nu1_sent(:)
                 Mu_nu1(i,j,ispec) = Mu_nu1_sent

!                 if (mod(ispec,2000)==0) then
!                    print *,'ipsec= ', ispec, 'phi_nu1 = ', phi_nu1(1,1,ispec,1)
!                 endif
              enddo
           enddo
        enddo
    endif



    ! todo
    !   call shift_velocities_from_f0(vp,vs,rhol,mul,lambdal)

    ! enddo
  end subroutine prepare_timerun_attenuation

  !  --------------------- --------------------- --------------------- --------------------- ---------------------

end program serial_specfem2D


! ******************
! meshing subroutine
! ******************

!--- global node number

integer function num(i,j,nx)

  implicit none

  integer :: i,j,nx

  num = j*(nx+1) + i + 1

end function num




!--------------- --------------- --------------- --------------- --------------- ---------------
subroutine compute_viscoelastic_forces(is_viscoelastic)
! compute elastic and viscoelastic forces
! 1 compute the gradient of displacement
! 2 compute the constitutive law : elastic if is_viscoelastic=false and viscoelesatic if is_viscoelastic=true
! 3 compute the stiffness integral

  use small_specfem_par
  implicit none

  logical, intent(in) :: is_viscoelastic
  integer :: i_sls
  real(kind=4) :: e1_sum,e11_sum,e13_sum
!  ! for gradient
!  real(kind=4) dux_dxi,duz_dxi,dux_dgamma,duz_dgamma
!  real(kind=4) dux_dxl,dux_dzl,duz_dxl,duz_dzl

  do ispec = 1,NSPEC

           tempx1(:,:) = 0.
           tempz1(:,:) = 0.
           tempx2(:,:) = 0.
           tempz2(:,:) = 0.
  ! il vaudrait mieux declarer i et j k en local
  do j = 1,NGLLZ
     do i = 1,NGLLX

        ! compute the gradient of displacement
        ! derivative along x and along z
            dux_dxi = 0.
            duz_dxi = 0.
            dux_dgamma = 0.
            duz_dgamma = 0.

        ! first double loop over GLL points to compute and store gradients
        ! we can merge the two loops because NGLLX == NGLLZ
        do k = 1,NGLLX
           dux_dxi = dux_dxi + displ(1,ibool(k,j,ispec))*hprime_xx(i,k)
           duz_dxi = duz_dxi + displ(2,ibool(k,j,ispec))*hprime_xx(i,k)
           dux_dgamma = dux_dgamma + displ(1,ibool(i,k,ispec))*hprime_zz(j,k)
           duz_dgamma = duz_dgamma + displ(2,ibool(i,k,ispec))*hprime_zz(j,k)
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)


        ! derivatives of displacement
        dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

        duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
        duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

        ! end of computing the gradient of displacement

          ! compute the constitutive law

         if (is_viscoelastic) then

          ! When implementing viscoelasticity according to the Carcione 1993 paper, attenuation is
          ! non-causal rather than causal i.e. wave speed up instead of slowing down
          ! when attenuation is turned on. We fixed that issue (which is not incorrect but non traditional)
          ! by taking the unrelaxed state (infinite frequency) as a reference instead of the relaxed state (zero frequency)
          ! and also using equations in Carcione's 2007 book.
          ! See file doc/old_problem_attenuation_reference_Specfem2D_fixed_by_Xie_Zhinan.pdf
          ! and doc/how_we_modified_Carcione_1993_to_make_it_causal_and_include_the_missing_1_over_L_factor.pdf

          ! See also J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
          ! vol. 58(1), p. 110-120 (1993) for two memory-variable mechanisms (page 112).

          ! and J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation
          ! in a linear viscoelastic medium, Geophysical Journal International,
          ! vol. 95, p. 597-611 (1988) for two memory-variable mechanisms (page 604).

            sigma_xx = lambdaplus2mu*dux_dxl + lambda*duz_dzl
            sigma_xz = mu*(duz_dxl + dux_dzl)
            sigma_zz = lambdaplus2mu*duz_dzl + lambda*dux_dxl

!           elastic modulus must be unrelaxed
!           sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
!           sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
!           sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dx


            e1_sum = 0.; e11_sum = 0.; e13_sum = 0.
            do i_sls = 1,N_SLS
            ! sum on memory variables
            ! print *,"e1_sum = ", e1_sum
                e1_sum = e1_sum + e1(i,j,ispec,i_sls)
                e11_sum = e11_sum + e11(i,j,ispec,i_sls)
                e13_sum = e13_sum + e13(i,j,ispec,i_sls)
            enddo


           ! use the right formula with 1/N included
           ! i.e. use the unrelaxed moduli here (see Carcione's book, third edition, equation (3.189))
            sigma_xx = sigma_xx + lambdaplusmu * e1_sum + 2.0 * mu* e11_sum
            sigma_xz = sigma_xz + mu * e13_sum
            sigma_zz = sigma_zz + lambdaplusmu * e1_sum - 2.0 * mu * e11_sum
            sigma_zx = sigma_xz
        else
            ! no attenuation
            sigma_xx = lambdaplus2mu*dux_dxl + lambda*duz_dzl
            sigma_xz = mu*(duz_dxl + dux_dzl)
            sigma_zz = lambdaplus2mu*duz_dzl + lambda*dux_dxl
            sigma_zx = sigma_xz

        endif
        ! end of computing the constitutive law

        ! compute the stiffness integral and forces at each GLL points
        ! weak formulation term based on stress tensor (non-symmetric form)
        ! also add GLL integration weights
        jacobianl = jacobian(i,j,ispec)

        tempx1(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
        tempz1(i,j) = wzgll(j)*jacobianl*(sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z

        tempx2(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
        tempz2(i,j) = wxgll(i)*jacobianl*(sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z

     enddo
  enddo  ! end of the loops on the collocation points i,j

  !
  ! second double-loop over GLL to compute all the terms
  !
      do j = 1,NGLLZ
         do i = 1,NGLLX
            iglob = ibool(i,j,ispec)
            ! along x direction and z direction
            ! and assemble the contributions
            ! we can merge the two loops because NGLLX == NGLLZ
            do k = 1,NGLLX
               if (.not. forced(iglob)) then
                ! here accel is a force not an acceleration we did divided by the mass matrix yet
                  accel(1,iglob) = accel(1,iglob) - (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
                  accel(2,iglob) = accel(2,iglob) - (tempz1(k,j)*hprimewgll_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))
               endif
            enddo
         enddo
      enddo ! second loop over the GLL points

     enddo   ! end of main loop on all the elements


end subroutine compute_viscoelastic_forces
!--------------- --------------- --------------- --------------- --------------- ---------------

subroutine compute_viscoacoustic_forces(is_viscoacoustic)
    use small_specfem_par
    implicit none

    logical, intent(in) :: is_viscoacoustic
    integer :: i_sls
    real(kind=4) :: e1_sum


 do ispec = 1,NSPEC
         tempx1(:,:) = 0.
         tempx2(:,:) = 0.

         ! print *,"is_viscoacoustic = "  , is_viscoacoustic

       do j = 1,NGLLZ
          do i = 1,NGLLX

             dux_dxi = 0.   !ux here is chi !
             dux_dgamma = 0.


             ! first double loop over GLL points to compute and store gradients
             ! we can merge the two loops because NGLLX == NGLLZ
             do k = 1,NGLLX
                dux_dxi = dux_dxi + potential_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)       !dchi/dxi
                dux_dgamma = dux_dgamma + potential_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k) !dchi/dgamma
             enddo



             xixl = xix(i,j,ispec)
             xizl = xiz(i,j,ispec)
             gammaxl = gammax(i,j,ispec)
             gammazl = gammaz(i,j,ispec)




             if (is_viscoacoustic) then


                e1_sum = 0.
                do i_sls = 1,N_SLS
                ! sum on memory variables
!                if (mod(ispec,3000)==0) then
!                if (i==1 .and. j==1) print *,"! sum on memory variables e1_sum = ", e1_sum
!                if (i==1 .and. j==1) print *, "e1(i=1,j=1,ispec=3000,i_sls) =", e1(i,j,ispec,i_sls)
!                endif

                e1_sum = e1_sum + e1(i,j,ispec,i_sls)
                enddo


                   ! derivatives of potential gradient of chi
             dux_dxl = dux_dxi * xixl + dux_dgamma * gammaxl !dchi/dx
             dux_dzl = dux_dxi * xizl + dux_dgamma * gammazl !dchi/dz

             dux_dxl=dux_dxl+e1_sum
             dux_dzl=dux_dzl+e1_sum

             ! Acoustic displacements = gradient of chi / rho
             iglob = ibool(i,j,ispec)
             acoustic_displ(1,iglob) = dux_dxl / rho        !u_x
             acoustic_displ(2,iglob) = dux_dzl / rho        !u_z

!              if (mod(ispec,3000)==0) then
!             if (i==1 .and. j==1) print *,"  "
!             if (i==1 .and. j==1) print *,"la calcul visco forces"
!             if (i==1 .and. j==1) print *, "dux_dxi,dux_dgamma =", dux_dxi,dux_dgamma
!             if (i==1 .and. j==1) print *, "e1_sum =", e1_sum
!             if (i==1 .and. j==1) print *,"potential_acoustic(3000)", potential_acoustic (3000)
!             if (i==1 .and. j==1) print *,'**********************************s**********************************************'
!             if (i==1 .and. j==1) print *,"  "
!             endif

             else

             !print *,"acoustic no att"


             ! derivatives of potential gradient of chi
             dux_dxl = dux_dxi * xixl + dux_dgamma * gammaxl !dchi/dx
             dux_dzl = dux_dxi * xizl + dux_dgamma * gammazl !dchi/dz

             ! Acoustic displacements = gradient of chi / rho
             iglob = ibool(i,j,ispec)
             acoustic_displ(1,iglob) = dux_dxl / rho        !u_x
             acoustic_displ(2,iglob) = dux_dzl / rho        !u_z

             endif

             jacobianl = jacobian(i,j,ispec)
             ! for acoustic medium also add integration weights

                ! for acoustics grad chi = -p = sigma_xx thus memory variable modify grad(chi))
                !tempx1(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
                !tempx2(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x



             tempx1(i,j) = wzgll(j) * jacobianl * (xixl * dux_dxl + xizl * dux_dzl) / rho
             tempx2(i,j) = wxgll(i) * jacobianl * (gammaxl * dux_dxl + gammazl * dux_dzl) / rho
          enddo
       enddo ! end of the loops on the collocation points i,j


       do j = 1,NGLLZ
          do i = 1,NGLLX
             iglob = ibool(i,j,ispec)
             ! along x direction and z direction
             ! and assemble the contributions
             ! we can merge the two loops because NGLLX == NGLLZ
             do k = 1,NGLLX
                if (.not. forced(iglob)) then
                   potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
                    (tempx1(k,j) * hprimewgll_xx(k,i) + tempx2(i,k) * hprimewgll_zz(k,j))
                endif
             enddo !enddo on k
          enddo
       enddo ! second loop over the GLL points

 enddo !end ispec end acoustics

end subroutine compute_viscoacoustic_forces
!--------------- --------------- --------------- --------------- --------------- ---------------


subroutine update_memory_variable(is_viscoelastic,is_viscoacoustic)

    use small_specfem_par
    implicit none
    logical, intent(in) :: is_viscoelastic,is_viscoacoustic

  integer :: i_sls
  ! for attenuation
  real(kind=4) :: phinu1,phinu2,theta_n_u
  double precision :: tauinvnu1,tauinvnu2



  ! compute Grad(displ_elastic) at time step n for attenuation
  !call compute_gradient_attenuation(displ_elastic,dux_dxl_n,duz_dxl_n, &
  !    dux_dzl_n,duz_dzl_n,xix,xiz,gammax,gammaz,ibool,ispec_is_elastic,hprime_xx,hprime_zz,nspec,nglob)

  ! we need :
  !  displ_elastic_old(:,:) = displ_elastic(:,:) + deltatsquareover2/TWO * accel_elastic(:,:)
  ! if CONVOLUTION_MEMORY_VARIABLES is used !cf update_displacement scheme l=518

  !compute Grad(disp_elastic_old) at time step n-1 for attenuation only if CONVOLUTION_MEMORY_VARIABLES is used
  !call compute_gradient_attenuation(displ_elastic_old,dux_dxl_nsub1,duz_dxl_nsub1, &
  !    dux_dzl_nsub1,duz_dzl_nsub1,xix,xiz,gammax,gammaz,ibool,ispec_is_elastic,hprime_xx,hprime_zz,nspec,nglob)

  !update memory variable in viscoelastic simulation
  ! loop over spectral elements and GLL points
  if (is_viscoelastic) then
  do ispec = 1,nspec
     do j = 1,NGLLZ
        do i = 1,NGLLX

         ! derivative along x and along z
            dux_dxi = 0.
            duz_dxi = 0.
            dux_dgamma = 0.
            duz_dgamma = 0.

     ! first double loop over GLL points to compute and store gradients
     ! we can merge the two loops because NGLLX == NGLLZ
        do k = 1,NGLLX
           dux_dxi = dux_dxi + displ(1,ibool(k,j,ispec))*hprime_xx(i,k)
           duz_dxi = duz_dxi + displ(2,ibool(k,j,ispec))*hprime_xx(i,k)
           dux_dgamma = dux_dgamma + displ(1,ibool(i,k,ispec))*hprime_zz(j,k)
           duz_dgamma = duz_dgamma + displ(2,ibool(i,k,ispec))*hprime_zz(j,k)
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)


        ! print *, "ispec = ", ispec, 'i= ' ,i,"j= ",j, "xixl = ", xixl ,'xizl = ', xizl, "gammaxl = ",gammaxl, "gammazl = ", gammazl

        ! derivatives of displacement
        dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

        duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
        duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

        theta_n_u = dux_dxl + duz_dzl


        do i_sls = 1,N_SLS
            phinu1 = phi_nu1(i,j,ispec,i_sls)
            tauinvnu1 = inv_tau_sigma_nu1(i,j,ispec,i_sls)
            phinu2 = phi_nu2(i,j,ispec,i_sls)
            tauinvnu2 = inv_tau_sigma_nu2(i,j,ispec,i_sls)

             !print *,"isl = ", i_sls," theta_n_u = " , theta_n_u ,"e1(i,j,ispec,i_sls)= ", e1(i,j,ispec,i_sls)

             e1(i,j,ispec,i_sls) = e1(i,j,ispec,i_sls) + deltat * &
                (- e1(i,j,ispec,i_sls)*tauinvnu1 + phinu1 * theta_n_u)

             e11(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls) + deltat * &
                (- e11(i,j,ispec,i_sls)*tauinvnu2 + phinu2 * (dux_dxl-theta_n_u/2.0))

             e13(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls) + deltat * &
                (- e13(i,j,ispec,i_sls)*tauinvnu2 + phinu2 * (dux_dzl + duz_dxl))
        enddo
     enddo
    enddo
  enddo
  endif
  if (is_viscoacoustic) then

    do ispec = 1,nspec
     do j = 1,NGLLZ
        do i = 1,NGLLX


            dux_dxi = 0.   !ux here is chi !
            dux_dgamma = 0.


             ! first double loop over GLL points to compute and store gradients
             ! we can merge the two loops because NGLLX == NGLLZ
             do k = 1,NGLLX
                dux_dxi = dux_dxi + potential_acoustic(ibool(k,j,ispec)) * hprime_xx(i,k)       !dchi/dxi
                dux_dgamma = dux_dgamma + potential_acoustic(ibool(i,k,ispec)) * hprime_zz(j,k) !dchi/dgamma
             enddo

             xixl = xix(i,j,ispec)
             xizl = xiz(i,j,ispec)
             gammaxl = gammax(i,j,ispec)
             gammazl = gammaz(i,j,ispec)


             !print *, "xixl,xizl =", xixl,xizl,gammaxl,gammazl
             !print *, "dux_dxi,dux_dgamma =", dux_dxi,dux_dgamma

             !if (i==1 .and. j==1) print *,"potential_acoustic(ibool(k,j,ispec) = ", potential_acoustic(ibool(i,j,ispec))

             ! derivatives of potential gradient of chi
             dux_dxl = dux_dxi * xixl + dux_dgamma * gammaxl
             dux_dzl = dux_dxi * xizl + dux_dgamma * gammazl

            !a verifier
            theta_n_u = dux_dxl ! + dux_dzl !!! a verifier theta__u= dux_dxl
            !a verifier





        do i_sls = 1,N_SLS
            phinu1 = phi_nu1(i,j,ispec,i_sls)
            tauinvnu1 = inv_tau_sigma_nu1(i,j,ispec,i_sls)
            e1(i,j,ispec,i_sls) = e1(i,j,ispec,i_sls) + deltat * (- e1(i,j,ispec,i_sls)*tauinvnu1 + phinu1 * theta_n_u)
        enddo

!           if (mod(ispec,3000)==0) then
!            if (i==1 .and. j==1) print *,"  "
!            if (i==1 .and. j==1) print *,"ici update memory"
!            if (i==1 .and. j==1) print *, "theta_n_u", theta_n_u, 'it =', it, 'dux_dgamma = ', dux_dgamma, "dux_dxi= ", dux_dxi
!            if (i==1 .and. j==1) print *,"potential_acoustic(3000)", potential_acoustic (3000)
!            if (i==1 .and. j==1) print *,"e1(i=1,j=1,ispec=3000,i_sls = 1)", e1(i,j,ispec,1)
!            if (i==1 .and. j==1) print *,"e1(i=1,j=1,ispec=3000,i_sls = 2)", e1(i,j,ispec,2)
!            if (i==1 .and. j==1) print *,"e1(i=1,j=1,ispec=3000,i_sls = 3)", e1(i,j,ispec,3)
!            if (i==1 .and. j==1) print *,"  "
!           endif

     enddo
    enddo
  enddo




  endif



  end subroutine update_memory_variable ! end of update memory variable in viscoelastic simulation




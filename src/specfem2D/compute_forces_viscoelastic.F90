
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, INRIA and University of Pau, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
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

subroutine compute_forces_viscoelastic(p_sv,nglob,nspec,myrank,nelemabs,numat, &
     ispec_selected_source,ispec_selected_rec,is_proc_source,which_proc_receiver, &
     source_type,it,NSTEP,anyabs,assign_external_model, &
     initialfield,ATTENUATION_VISCOELASTIC_SOLID,anglesource,deltatcube, &
     deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,elastic,codeabs, &
     accel_elastic,veloc_elastic,displ_elastic,b_accel_elastic,b_displ_elastic, &
     density,poroelastcoef,xix,xiz,gammax,gammaz, &
     jacobian,vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,anisotropic,anisotropy, &
     source_time_function,sourcearray,adj_sourcearrays,e1,e11, &
     e13,dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
     dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1,hprime_xx,hprimewgll_xx, &
     hprime_zz,hprimewgll_zz,wxgll,wzgll,inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2,N_SLS, &
     deltat,coord,add_Bielak_conditions, &
     x0_source, z0_source, A_plane, B_plane, C_plane, anglesource_refl, c_inc, c_refl, time_offset,f0, &
     v0x_left,v0z_left,v0x_right,v0z_right,v0x_bot,v0z_bot,t0x_left,t0z_left,t0x_right,t0z_right,t0x_bot,t0z_bot,&
     nleft,nright,nbot,over_critical_angle,NSOURCES,nrec,SIMULATION_TYPE,SAVE_FORWARD,b_absorb_elastic_left,&
     b_absorb_elastic_right,b_absorb_elastic_bottom,b_absorb_elastic_top,nspec_left,nspec_right,&
     nspec_bottom,nspec_top,ib_left,ib_right,ib_bottom,ib_top,mu_k,kappa_k,&
     e1_LDDRK,e11_LDDRK,e13_LDDRK,alpha_LDDRK,beta_LDDRK, & 
     e1_initial_rk,e11_initial_rk,e13_initial_rk,e1_force_RK, e11_force_RK, e13_force_RK, &
     stage_time_scheme,i_stage,ADD_SPRING_TO_STACEY,x_center_spring,z_center_spring,nadj_rec_local, &
     is_PML,nspec_PML,spec_to_PML,region_CPML, &
     K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store, &
     rmemory_displ_elastic,rmemory_dux_dx,rmemory_dux_dz,rmemory_duz_dx,rmemory_duz_dz, &
     rmemory_dux_dx_prime,rmemory_dux_dz_prime,rmemory_duz_dx_prime,rmemory_duz_dz_prime, &
     rmemory_displ_elastic_LDDRK,rmemory_dux_dx_LDDRK,rmemory_dux_dz_LDDRK,& 
     rmemory_duz_dx_LDDRK,rmemory_duz_dz_LDDRK, &  
     PML_BOUNDARY_CONDITIONS,ROTATE_PML_ACTIVATE,ROTATE_PML_ANGLE)  


  ! compute forces for the elastic elements

  implicit none

  include "constants.h"

  logical :: p_sv
  integer :: NSOURCES, i_source
  integer :: nglob,nspec,myrank,nelemabs,numat,it,NSTEP
  integer, dimension(NSOURCES) :: ispec_selected_source,is_proc_source,source_type

  integer :: nrec,SIMULATION_TYPE
  integer, dimension(nrec) :: ispec_selected_rec,which_proc_receiver
  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top
  integer, dimension(nelemabs) :: ib_left
  integer, dimension(nelemabs) :: ib_right
  integer, dimension(nelemabs) :: ib_bottom
  integer, dimension(nelemabs) :: ib_top
  integer :: stage_time_scheme,i_stage,nadj_rec_local

  logical :: anyabs,assign_external_model,initialfield,ATTENUATION_VISCOELASTIC_SOLID,add_Bielak_conditions
  logical :: ADD_SPRING_TO_STACEY
  real(kind=CUSTOM_REAL) :: x_center_spring,z_center_spring

  logical :: SAVE_FORWARD

  double precision :: deltatcube,deltatfourth,twelvedeltat,fourdeltatsquare
  double precision, dimension(NSOURCES) :: anglesource

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato
  integer, dimension(nelemabs) :: numabs

  logical, dimension(nspec) :: elastic,anisotropic
  logical, dimension(4,nelemabs)  :: codeabs

  real(kind=CUSTOM_REAL), dimension(3,nglob) :: accel_elastic,veloc_elastic,displ_elastic
  real(kind=CUSTOM_REAL), dimension(3,nglob) :: temp_displ_elastic
  double precision, dimension(2,numat) :: density
  double precision, dimension(4,3,numat) :: poroelastcoef
  double precision, dimension(6,numat) :: anisotropy
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext,vsext,rhoext
  double precision, dimension(NGLLX,NGLLZ,nspec) ::  c11ext,c15ext,c13ext,c33ext,c35ext,c55ext

  real(kind=CUSTOM_REAL), dimension(NSOURCES,NSTEP,stage_time_scheme) :: source_time_function
  real(kind=CUSTOM_REAL), dimension(NSOURCES,NDIM,NGLLX,NGLLZ) :: sourcearray

  real(kind=CUSTOM_REAL), dimension(3,nglob) :: b_accel_elastic,b_displ_elastic
  real(kind=CUSTOM_REAL), dimension(nadj_rec_local,NSTEP,3,NGLLX,NGLLZ) :: adj_sourcearrays
  real(kind=CUSTOM_REAL), dimension(nglob) :: mu_k,kappa_k
  real(kind=CUSTOM_REAL), dimension(3,NGLLZ,nspec_left,NSTEP) :: b_absorb_elastic_left
  real(kind=CUSTOM_REAL), dimension(3,NGLLZ,nspec_right,NSTEP) :: b_absorb_elastic_right
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,nspec_top,NSTEP) :: b_absorb_elastic_top
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,nspec_bottom,NSTEP) :: b_absorb_elastic_bottom

  integer :: N_SLS
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec,N_SLS) :: e1,e11,e13
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec,N_SLS) :: e1_LDDRK,e11_LDDRK,e13_LDDRK
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec,N_SLS) :: e1_initial_rk,e11_initial_rk,e13_initial_rk
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec,N_SLS,stage_time_scheme) :: e1_force_RK, e11_force_RK, e13_force_RK
  double precision, dimension(NGLLX,NGLLZ,nspec,N_SLS) :: inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
  double precision, dimension(NGLLX,NGLLZ,nspec) :: Mu_nu1,Mu_nu2
  real(kind=CUSTOM_REAL) :: e1_sum,e11_sum,e13_sum
  integer :: i_sls

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: &
       dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n,dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1

  ! derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

  ! Gauss-Lobatto-Legendre weights
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  real(kind=CUSTOM_REAL), dimension(NGLLZ) :: wzgll

  ! Parameter for LDDRK time scheme
  real(kind=CUSTOM_REAL), dimension(Nstages) :: alpha_LDDRK,beta_LDDRK

  !temp variable
  real(kind=CUSTOM_REAL) :: temp_time_scheme,temper_time_scheme,weight_rk


  !---
  !--- local variables
  !---

  integer :: ispec,i,j,k,iglob,ispecabs,ibegin,iend,irec,irec_local

  ! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duy_dxi,duy_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duy_dxl,duz_dxl,dux_dzl,duy_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: dux_dxl_prime,duy_dxl_prime,duz_dxl_prime,dux_dzl_prime,duy_dzl_prime,duz_dzl_prime 
  real(kind=CUSTOM_REAL) :: theta,ct,st 
  real(kind=CUSTOM_REAL) :: b_dux_dxi,b_dux_dgamma,b_duy_dxi,b_duy_dgamma,b_duz_dxi,b_duz_dgamma
  real(kind=CUSTOM_REAL) :: b_dux_dxl,b_duy_dxl,b_duz_dxl,b_dux_dzl,b_duy_dzl,b_duz_dzl
  real(kind=CUSTOM_REAL) :: dsxx,dsxz,dszz
  real(kind=CUSTOM_REAL) :: b_dsxx,b_dsxz,b_dszz
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_xy,sigma_xz,sigma_zy,sigma_zz,sigma_zx
  real(kind=CUSTOM_REAL) :: sigma_xx_prime,sigma_xy_prime,sigma_xz_prime,sigma_zy_prime,sigma_zz_prime,sigma_zx_prime 
  real(kind=CUSTOM_REAL) :: b_sigma_xx,b_sigma_xy,b_sigma_xz,b_sigma_zy,b_sigma_zz,b_sigma_zx
  real(kind=CUSTOM_REAL) :: nx,nz,vx,vy,vz,vn,rho_vp,rho_vs,tx,ty,tz,weight,xxi,zxi,xgamma,zgamma,jacobian1D
  real(kind=CUSTOM_REAL) :: displx,disply,displz,displn,spring_position,displtx,displty,displtz

!! DK DK added this for Guenneau, March 2012
#ifdef USE_GUENNEAU
  integer :: kmato_ispec_outside_Guenneau
  real(kind=CUSTOM_REAL) :: ct, st, r, a, inva, lambda, mu, x, y, &
                            lambdaplus2mu, ct2 , ct3 , ct4 , twoct2 , st2 , st3 , st4
  real(kind=CUSTOM_REAL) :: epsilon_xx,epsilon_xz,epsilon_zx,epsilon_zz
  real(kind=CUSTOM_REAL) :: C1111, C1112, C1121, C1122, C1211, C1212, C1221, C1222, C2111, C2112, C2121, &
                            C2122, C2211, C2212, C2221, C2222
#endif
!! DK DK added this for Guenneau, March 2012

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: tempx1,tempx2,tempy1,tempy2,tempz1,tempz2
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ) :: b_tempx1,b_tempx2,b_tempy1,b_tempy2,b_tempz1,b_tempz2

  ! Jacobian matrix and determinant
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl,jacobianl

  ! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic, &
    lambdaplus2mu_unrelaxed_elastic,kappal,cpl,csl,rhol, &
    lambdal_relaxed_viscoelastic,mul_relaxed_viscoelastic,lambdalplus2mul_relaxed_viscoel

  ! for attenuation
  real(kind=CUSTOM_REAL) :: Un,Unp1,tauinv,Sn,Snp1,theta_n,theta_np1,tauinvsquare,tauinvcube,tauinvUn

  ! for anisotropy
  double precision ::  c11,c15,c13,c33,c35,c55

  ! for analytical initial plane wave for Bielak's conditions
  double precision :: veloc_horiz,veloc_vert,dxUx,dzUx,dxUz,dzUz,traction_x_t0,traction_z_t0,deltat
  double precision, dimension(NDIM,nglob), intent(in) :: coord
  double precision x0_source, z0_source, anglesource_refl, c_inc, c_refl, time_offset, f0
  double precision, dimension(NDIM) :: A_plane, B_plane, C_plane
  !over critical angle
  logical :: over_critical_angle
  integer :: nleft, nright, nbot
  double precision, dimension(nleft) :: v0x_left,v0z_left,t0x_left,t0z_left
  double precision, dimension(nright) :: v0x_right,v0z_right,t0x_right,t0z_right
  double precision, dimension(nbot) :: v0x_bot,v0z_bot,t0x_bot,t0z_bot
  integer count_left,count_right,count_bottom

  integer :: ifirstelem,ilastelem

!CPML coefficients and memory variables
  integer :: nspec_PML,ispec_PML
  integer, dimension(nspec) :: region_CPML
  logical, dimension(nspec) :: is_PML
  integer, dimension(nspec) :: spec_to_PML
  logical :: PML_BOUNDARY_CONDITIONS,ROTATE_PML_ACTIVATE 
  double precision ROTATE_PML_ANGLE                      

  real(kind=CUSTOM_REAL), dimension(2,3,NGLLX,NGLLZ,nspec_PML) :: rmemory_displ_elastic
  real(kind=CUSTOM_REAL), dimension(2,3,NGLLX,NGLLZ,nspec_PML) :: rmemory_displ_elastic_LDDRK
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML) :: &
    rmemory_dux_dx,rmemory_dux_dz,rmemory_duz_dx,rmemory_duz_dz
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML) :: &   
    rmemory_dux_dx_prime,rmemory_dux_dz_prime,rmemory_duz_dx_prime,rmemory_duz_dz_prime
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML) :: &
    rmemory_dux_dx_LDDRK,rmemory_dux_dz_LDDRK,rmemory_duz_dx_LDDRK,rmemory_duz_dz_LDDRK
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML) :: &
                  K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLZ,nspec_PML) :: accel_elastic_PML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML) ::PML_dux_dxl,PML_dux_dzl,PML_duz_dxl,PML_duz_dzl,&
                           PML_dux_dxl_new,PML_dux_dzl_new,PML_duz_dxl_new,PML_duz_dzl_new
  real(kind=CUSTOM_REAL) :: coef0, coef1, coef2,bb

  real(kind=CUSTOM_REAL) :: A0, A1, A2, A3, A4, A5, A6, A7, A8

  real(kind=CUSTOM_REAL) :: dux_dxi_new,dux_dgamma_new,duz_dxi_new,duz_dgamma_new
  real(kind=CUSTOM_REAL) :: dux_dxl_new,dux_dzl_new,duz_dxl_new,duz_dzl_new
  real(kind=CUSTOM_REAL), dimension(3,nglob) :: displ_elastic_new

! this to avoid a warning at execution time about an undefined variable being used
! for the SH component in the case of a P-SV calculation, and vice versa
  sigma_xx = 0
  sigma_xy = 0
  sigma_xz = 0
  sigma_zy = 0
  sigma_zz = 0
  sigma_zx = 0

  if( PML_BOUNDARY_CONDITIONS ) then
    accel_elastic_PML = 0._CUSTOM_REAL
    PML_dux_dxl = 0._CUSTOM_REAL
    PML_dux_dzl = 0._CUSTOM_REAL
    PML_duz_dxl = 0._CUSTOM_REAL
    PML_duz_dzl = 0._CUSTOM_REAL
    PML_dux_dxl_new = 0._CUSTOM_REAL
    PML_dux_dzl_new = 0._CUSTOM_REAL
    PML_duz_dxl_new = 0._CUSTOM_REAL
    PML_duz_dzl_new = 0._CUSTOM_REAL
    displ_elastic_new = displ_elastic + deltat * veloc_elastic
  endif


  ! compute Grad(displ_elastic) at time step n for attenuation
  if(ATTENUATION_VISCOELASTIC_SOLID) then
     temp_displ_elastic = displ_elastic + deltat * veloc_elastic + 0.5 * deltat**2 * accel_elastic
     call compute_gradient_attenuation(displ_elastic,dux_dxl_n,duz_dxl_n, &
          dux_dzl_n,duz_dzl_n,xix,xiz,gammax,gammaz,ibool,elastic,hprime_xx,hprime_zz,nspec,nglob)
  endif

  ifirstelem = 1
  ilastelem = nspec

  ! loop over spectral elements
  do ispec = ifirstelem,ilastelem

     tempx1(:,:) = ZERO
     tempy1(:,:) = ZERO
     tempz1(:,:) = ZERO
     tempx2(:,:) = ZERO
     tempy2(:,:) = ZERO
     tempz2(:,:) = ZERO
     if(SIMULATION_TYPE ==2)then
        b_tempx1(:,:) = ZERO
        b_tempy1(:,:) = ZERO
        b_tempz1(:,:) = ZERO
        b_tempx2(:,:) = ZERO
        b_tempy2(:,:) = ZERO
        b_tempz2(:,:) = ZERO
     endif

     !---
     !--- elastic spectral element
     !---
     if(elastic(ispec)) then

        ! get unrelaxed elastic parameters of current spectral element
        lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
        mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
        lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))

        ! first double loop over GLL points to compute and store gradients
        do j = 1,NGLLZ
           do i = 1,NGLLX

              !--- if external medium, get elastic parameters of current grid point
              if(assign_external_model) then
                 cpl = vpext(i,j,ispec)
                 csl = vsext(i,j,ispec)
                 rhol = rhoext(i,j,ispec)
                 mul_unrelaxed_elastic = rhol*csl*csl
                 lambdal_unrelaxed_elastic = rhol*cpl*cpl - TWO*mul_unrelaxed_elastic
                 lambdaplus2mu_unrelaxed_elastic = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic
              endif

              ! derivative along x and along z
              dux_dxi = ZERO
              duy_dxi = ZERO
              duz_dxi = ZERO

              dux_dgamma = ZERO
              duy_dgamma = ZERO
              duz_dgamma = ZERO

              if(SIMULATION_TYPE == 2) then ! Adjoint calculation, backward wavefield
                 b_dux_dxi = ZERO
                 b_duy_dxi = ZERO
                 b_duz_dxi = ZERO

                 b_dux_dgamma = ZERO
                 b_duy_dgamma = ZERO
                 b_duz_dgamma = ZERO
              endif

              ! first double loop over GLL points to compute and store gradients
              ! we can merge the two loops because NGLLX == NGLLZ
              do k = 1,NGLLX
                 dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
                 duy_dxi = duy_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
                 duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
                 dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
                 duy_dgamma = duy_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
                 duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)

                 if(SIMULATION_TYPE == 2) then ! Adjoint calculation, backward wavefield
                    b_dux_dxi = b_dux_dxi + b_displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
                    b_duy_dxi = b_duy_dxi + b_displ_elastic(2,ibool(k,j,ispec))*hprime_xx(i,k)
                    b_duz_dxi = b_duz_dxi + b_displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
                    b_dux_dgamma = b_dux_dgamma + b_displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
                    b_duy_dgamma = b_duy_dgamma + b_displ_elastic(2,ibool(i,k,ispec))*hprime_zz(j,k)
                    b_duz_dgamma = b_duz_dgamma + b_displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
                 endif
              enddo

              xixl = xix(i,j,ispec)
              xizl = xiz(i,j,ispec)
              gammaxl = gammax(i,j,ispec)
              gammazl = gammaz(i,j,ispec)

              ! derivatives of displacement
              dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
              dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

              duy_dxl = duy_dxi*xixl + duy_dgamma*gammaxl
              duy_dzl = duy_dxi*xizl + duy_dgamma*gammazl

              duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
              duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

              if(PML_BOUNDARY_CONDITIONS .and. is_PML(ispec)) then
                  ispec_PML=spec_to_PML(ispec)

                  PML_dux_dxl(i,j,ispec_PML) = dux_dxl
                  PML_dux_dzl(i,j,ispec_PML)=dux_dzl
                  PML_duz_dzl(i,j,ispec_PML)=duz_dzl
                  PML_duz_dxl(i,j,ispec_PML)=duz_dxl

                  ! derivative along x and along z
                  dux_dxi_new = ZERO
                  duz_dxi_new = ZERO
                  dux_dgamma_new = ZERO
                  duz_dgamma_new = ZERO

                  ! first double loop over GLL points to compute and store gradients
                  ! we can merge the two loops because NGLLX == NGLLZ
                  do k = 1,NGLLX
                    dux_dxi_new = dux_dxi_new + displ_elastic_new(1,ibool(k,j,ispec))*hprime_xx(i,k)
                    duz_dxi_new = duz_dxi_new + displ_elastic_new(3,ibool(k,j,ispec))*hprime_xx(i,k)
                    dux_dgamma_new = dux_dgamma_new + displ_elastic_new(1,ibool(i,k,ispec))*hprime_zz(j,k)
                    duz_dgamma_new = duz_dgamma_new + displ_elastic_new(3,ibool(i,k,ispec))*hprime_zz(j,k)
                  enddo

                  xixl = xix(i,j,ispec)
                  xizl = xiz(i,j,ispec)
                  gammaxl = gammax(i,j,ispec)
                  gammazl = gammaz(i,j,ispec)

                  ! derivatives of displacement
                  dux_dxl_new = dux_dxi_new*xixl + dux_dgamma_new*gammaxl
                  dux_dzl_new = dux_dxi_new*xizl + dux_dgamma_new*gammazl
                  duz_dxl_new = duz_dxi_new*xixl + duz_dgamma_new*gammaxl
                  duz_dzl_new = duz_dxi_new*xizl + duz_dgamma_new*gammazl

                  PML_dux_dxl_new(i,j,ispec_PML) = dux_dxl_new
                  PML_dux_dzl_new(i,j,ispec_PML)=dux_dzl_new
                  PML_duz_dzl_new(i,j,ispec_PML)=duz_dzl_new
                  PML_duz_dxl_new(i,j,ispec_PML)=duz_dxl_new
              endif


                if(is_PML(ispec) .and. PML_BOUNDARY_CONDITIONS) then
                  ispec_PML=spec_to_PML(ispec)
!------------------------------------------------------------------------------
!---------------------------- LEFT & RIGHT ------------------------------------
!------------------------------------------------------------------------------
                   if (region_CPML(ispec) == CPML_LEFT .or. region_CPML(ispec) == CPML_RIGHT) then

                    !---------------------- A8--------------------------
                    A8 = - d_x_store(i,j,ispec_PML) / (k_x_store(i,j,ispec_PML) ** 2)
                    bb = d_x_store(i,j,ispec_PML) / k_x_store(i,j,ispec_PML) + alpha_x_store(i,j,ispec_PML)

                    if(stage_time_scheme == 1) then

                    coef0 = exp(-bb * deltat)
                    if ( abs(bb) > 1.d-3 ) then
                      coef1 = (1.d0 - exp(-bb * deltat / 2.d0)) / bb
                      coef2 = (1.d0 - exp(-bb* deltat / 2.d0)) * exp(-bb * deltat / 2.d0)/ bb
                    else
                      coef1 = deltat / 2.0d0
                      coef2 = deltat / 2.0d0
                    end if

                    if(ROTATE_PML_ACTIVATE)then
                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dx(i,j,ispec_PML) = coef0 * rmemory_dux_dx(i,j,ispec_PML) &
                    + PML_dux_dxl_new(i,j,ispec_PML) * coef1 + PML_dux_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dz(i,j,ispec_PML) = coef0 * rmemory_dux_dz(i,j,ispec_PML) &
                    + PML_dux_dzl_new(i,j,ispec_PML) * coef1 + PML_dux_dzl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dx(i,j,ispec_PML) = coef0 * rmemory_duz_dx(i,j,ispec_PML) &
                    + PML_duz_dxl_new(i,j,ispec_PML) * coef1 + PML_duz_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dz(i,j,ispec_PML) = coef0 * rmemory_duz_dz(i,j,ispec_PML) &
                    + PML_duz_dzl_new(i,j,ispec_PML) * coef1 + PML_duz_dzl(i,j,ispec_PML) * coef2
                    else
                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dx(i,j,ispec_PML) = coef0 * rmemory_dux_dx(i,j,ispec_PML) &
                    + PML_dux_dxl_new(i,j,ispec_PML) * coef1 + PML_dux_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dx(i,j,ispec_PML) = coef0 * rmemory_duz_dx(i,j,ispec_PML) &
                    + PML_duz_dxl_new(i,j,ispec_PML) * coef1 + PML_duz_dxl(i,j,ispec_PML) * coef2
                    endif

                    endif 

                    if(stage_time_scheme == 6) then 

                     rmemory_dux_dx_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_dux_dx(i,j,ispec_PML) + PML_dux_dxl(i,j,ispec_PML))
                     rmemory_dux_dx(i,j,ispec_PML) = rmemory_dux_dx(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML)

                     rmemory_duz_dx_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_duz_dx(i,j,ispec_PML) + PML_duz_dxl(i,j,ispec_PML))
                     rmemory_duz_dx(i,j,ispec_PML) = rmemory_duz_dx(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML)

                    end if 
 
                    if(ROTATE_PML_ACTIVATE)then
                    dux_dxl = PML_dux_dxl(i,j,ispec_PML)  + A8 * rmemory_dux_dx(i,j,ispec_PML)
                    dux_dzl = PML_dux_dzl(i,j,ispec_PML)  + A8 * rmemory_dux_dz(i,j,ispec_PML)
                    duz_dxl = PML_duz_dxl(i,j,ispec_PML)  + A8 * rmemory_duz_dx(i,j,ispec_PML)
                    duz_dzl = PML_duz_dzl(i,j,ispec_PML)  + A8 * rmemory_duz_dz(i,j,ispec_PML)
                    else
                    dux_dxl = PML_dux_dxl(i,j,ispec_PML)  + A8 * rmemory_dux_dx(i,j,ispec_PML)
                    duz_dxl = PML_duz_dxl(i,j,ispec_PML)  + A8 * rmemory_duz_dx(i,j,ispec_PML)
                    endif



                    !---------------------- A5--------------------------
                    A5 = d_x_store(i,j,ispec_PML)
                    bb = alpha_x_store(i,j,ispec_PML)

                    if(stage_time_scheme == 1) then

                    coef0 = exp(- bb * deltat)

                    if ( abs( bb ) > 1.d-3) then
                      coef1 = (1.0d0 - exp(- bb * deltat / 2.0d0) ) / bb
                      coef2 = (1.0d0 - exp(- bb * deltat / 2.0d0) ) * exp(- bb * deltat / 2.0d0) / bb
                    else
                      coef1 = deltat / 2.0d0
                      coef2 = deltat / 2.0d0
                    end if

                    if(ROTATE_PML_ACTIVATE)then
                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dx_prime(i,j,ispec_PML) = coef0 * rmemory_dux_dx_prime(i,j,ispec_PML) &
                    + PML_dux_dxl_new(i,j,ispec_PML) *coef1 + PML_dux_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dz_prime(i,j,ispec_PML) = coef0 * rmemory_dux_dz_prime(i,j,ispec_PML) &
                    + PML_dux_dzl_new(i,j,ispec_PML) *coef1 + PML_dux_dzl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dx_prime(i,j,ispec_PML) = coef0 * rmemory_duz_dx_prime(i,j,ispec_PML) &
                    + PML_duz_dxl_new(i,j,ispec_PML) *coef1 + PML_duz_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dz_prime(i,j,ispec_PML) = coef0 * rmemory_duz_dz_prime(i,j,ispec_PML) &
                    + PML_duz_dzl_new(i,j,ispec_PML) *coef1 + PML_duz_dzl(i,j,ispec_PML) * coef2
                    else
                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dz(i,j,ispec_PML) = coef0 * rmemory_dux_dz(i,j,ispec_PML) &
                    + PML_dux_dzl_new(i,j,ispec_PML) *coef1 + PML_dux_dzl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dz(i,j,ispec_PML) = coef0 * rmemory_duz_dz(i,j,ispec_PML) &
                    + PML_duz_dzl_new(i,j,ispec_PML) *coef1 + PML_duz_dzl(i,j,ispec_PML) * coef2
                    endif

                    end if

                    if(stage_time_scheme == 6) then
                     rmemory_dux_dz_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_dux_dz(i,j,ispec_PML) + PML_dux_dzl(i,j,ispec_PML))
                     rmemory_dux_dz(i,j,ispec_PML) = rmemory_dux_dz(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML)

                     rmemory_duz_dz_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_duz_dz(i,j,ispec_PML) + PML_duz_dzl(i,j,ispec_PML))
                     rmemory_duz_dz(i,j,ispec_PML) = rmemory_duz_dz(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML)
                    end if 

                    if(ROTATE_PML_ACTIVATE)then
                    dux_dxl_prime = PML_dux_dxl(i,j,ispec_PML)  + A5 * rmemory_dux_dx_prime(i,j,ispec_PML)
                    dux_dzl_prime = PML_dux_dzl(i,j,ispec_PML)  + A5 * rmemory_dux_dz_prime(i,j,ispec_PML)
                    duz_dxl_prime = PML_duz_dxl(i,j,ispec_PML)  + A5 * rmemory_duz_dx_prime(i,j,ispec_PML)
                    duz_dzl_prime = PML_duz_dzl(i,j,ispec_PML)  + A5 * rmemory_duz_dz_prime(i,j,ispec_PML)
                    else
                    dux_dzl = PML_dux_dzl(i,j,ispec_PML)  + A5 * rmemory_dux_dz(i,j,ispec_PML)
                    duz_dzl = PML_duz_dzl(i,j,ispec_PML)  + A5 * rmemory_duz_dz(i,j,ispec_PML)
                    endif




!------------------------------------------------------------------------------
!---------------------------- CORNER ------------------------------------------
!------------------------------------------------------------------------------
                   elseif (region_CPML(ispec) == CPML_TOP_LEFT .or. region_CPML(ispec) == CPML_TOP_RIGHT .or. &
                           region_CPML(ispec) == CPML_BOTTOM_LEFT .or. region_CPML(ispec) == CPML_BOTTOM_RIGHT) then

                    !---------------------- A8--------------------------
                    A8 = (k_x_store(i,j,ispec_PML) * d_z_store(i,j,ispec_PML) &
                          - k_z_store(i,j,ispec_PML) * d_x_store(i,j,ispec_PML)) &
                         / (k_x_store(i,j,ispec_PML)**2)

                    bb = d_x_store(i,j,ispec_PML) / k_x_store(i,j,ispec_PML) + alpha_x_store(i,j,ispec_PML)

                    if(stage_time_scheme == 1) then 

                    coef0 = exp(- bb * deltat)

                    if ( abs(bb) > 1.d-3 ) then
                      coef1 = ( 1.d0 - exp(- bb * deltat / 2.d0) ) / bb
                      coef2 = ( 1.d0 - exp(- bb * deltat / 2.d0) ) * exp(- bb * deltat / 2.d0) / bb
                    else
                      coef1 = deltat / 2.0d0
                      coef2 = deltat / 2.0d0
                    end if

                    if(ROTATE_PML_ACTIVATE)then
                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dx(i,j,ispec_PML) = coef0 * rmemory_dux_dx(i,j,ispec_PML) &
                    + PML_dux_dxl_new(i,j,ispec_PML) * coef1 + PML_dux_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dz(i,j,ispec_PML) = coef0 * rmemory_dux_dz(i,j,ispec_PML) &
                    + PML_dux_dzl_new(i,j,ispec_PML) * coef1 + PML_dux_dzl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dx(i,j,ispec_PML) = coef0 * rmemory_duz_dx(i,j,ispec_PML) &
                    + PML_duz_dxl_new(i,j,ispec_PML) * coef1 + PML_duz_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dz(i,j,ispec_PML) = coef0 * rmemory_duz_dz(i,j,ispec_PML) &
                    + PML_duz_dzl_new(i,j,ispec_PML) * coef1 + PML_duz_dzl(i,j,ispec_PML) * coef2
                    else
                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dx(i,j,ispec_PML) = coef0*rmemory_dux_dx(i,j,ispec_PML) &
                    + PML_dux_dxl_new(i,j,ispec_PML) * coef1 + PML_dux_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dx(i,j,ispec_PML) = coef0*rmemory_duz_dx(i,j,ispec_PML) &
                    + PML_duz_dxl_new(i,j,ispec_PML) * coef1 + PML_duz_dxl(i,j,ispec_PML) * coef2
                    endif

                    end if

                    if(stage_time_scheme == 6) then 
                     rmemory_dux_dx_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_dux_dx(i,j,ispec_PML) + PML_dux_dxl(i,j,ispec_PML))
                     rmemory_dux_dx(i,j,ispec_PML) = rmemory_dux_dx(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML)

                     rmemory_duz_dx_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_duz_dx(i,j,ispec_PML) + PML_duz_dxl(i,j,ispec_PML))
                     rmemory_duz_dx(i,j,ispec_PML) = rmemory_duz_dx(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML)
                    end if  

                    if(ROTATE_PML_ACTIVATE)then
                    dux_dxl = PML_dux_dxl(i,j,ispec_PML)  + A8 * rmemory_dux_dx(i,j,ispec_PML)
                    dux_dzl = PML_dux_dzl(i,j,ispec_PML)  + A8 * rmemory_dux_dz(i,j,ispec_PML)
                    duz_dxl = PML_duz_dxl(i,j,ispec_PML)  + A8 * rmemory_duz_dx(i,j,ispec_PML)
                    duz_dzl = PML_duz_dzl(i,j,ispec_PML)  + A8 * rmemory_duz_dz(i,j,ispec_PML)
                    else
                    dux_dxl = PML_dux_dxl(i,j,ispec_PML)  + A8 * rmemory_dux_dx(i,j,ispec_PML)
                    duz_dxl = PML_duz_dxl(i,j,ispec_PML)  + A8 * rmemory_duz_dx(i,j,ispec_PML)
                    endif



                    !---------------------------- A5 ----------------------------
                    A5 =(k_z_store(i,j,ispec_PML) * d_x_store(i,j,ispec_PML) &
                         - k_x_store(i,j,ispec_PML) * d_z_store(i,j,ispec_PML)) &
                          / (k_z_store(i,j,ispec_PML)**2)

                    bb = d_z_store(i,j,ispec_PML) / k_z_store(i,j,ispec_PML) + alpha_z_store(i,j,ispec_PML)

                    if(stage_time_scheme == 1) then 

                    coef0 = exp(- bb * deltat)

                    if ( abs(bb) > 1.d-3 ) then
                      coef1 = ( 1.d0 - exp(- bb * deltat / 2.d0) ) / bb
                      coef2 = ( 1.d0 - exp(- bb * deltat / 2.d0) ) * exp(- bb * deltat / 2.d0) / bb
                    else
                      coef1 = deltat / 2.0d0
                      coef2 = deltat / 2.0d0
                    end if

                    if(ROTATE_PML_ACTIVATE)then

                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dx_prime(i,j,ispec_PML) = coef0 * rmemory_dux_dx_prime(i,j,ispec_PML) &
                    + PML_dux_dxl_new(i,j,ispec_PML) *coef1 + PML_dux_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dz_prime(i,j,ispec_PML) = coef0 * rmemory_dux_dz_prime(i,j,ispec_PML) &
                    + PML_dux_dzl_new(i,j,ispec_PML) *coef1 + PML_dux_dzl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dx_prime(i,j,ispec_PML) = coef0 * rmemory_duz_dx_prime(i,j,ispec_PML) &
                    + PML_duz_dxl_new(i,j,ispec_PML) *coef1 + PML_duz_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dz_prime(i,j,ispec_PML) = coef0 * rmemory_duz_dz_prime(i,j,ispec_PML) &
                    + PML_duz_dzl_new(i,j,ispec_PML) *coef1 + PML_duz_dzl(i,j,ispec_PML) * coef2

                    else

                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dz(i,j,ispec_PML) = coef0 * rmemory_dux_dz(i,j,ispec_PML) &
                    + PML_dux_dzl_new(i,j,ispec_PML) *coef1 + PML_dux_dzl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dz(i,j,ispec_PML) = coef0 * rmemory_duz_dz(i,j,ispec_PML) &
                    + PML_duz_dzl_new(i,j,ispec_PML) *coef1 + PML_duz_dzl(i,j,ispec_PML) * coef2
                    endif

                    end if

                    if(stage_time_scheme == 6) then  
                     rmemory_dux_dz_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_dux_dz(i,j,ispec_PML) + PML_dux_dzl(i,j,ispec_PML))
                     rmemory_dux_dz(i,j,ispec_PML) = rmemory_dux_dz(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML)

                     rmemory_duz_dz_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_duz_dz(i,j,ispec_PML) + PML_duz_dzl(i,j,ispec_PML))
                     rmemory_duz_dz(i,j,ispec_PML) = rmemory_duz_dz(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML)
                    end if   

                    if(ROTATE_PML_ACTIVATE)then
                    dux_dxl_prime = PML_dux_dxl(i,j,ispec_PML)  + A5 * rmemory_dux_dx_prime(i,j,ispec_PML)
                    dux_dzl_prime = PML_dux_dzl(i,j,ispec_PML)  + A5 * rmemory_dux_dz_prime(i,j,ispec_PML)
                    duz_dxl_prime = PML_duz_dxl(i,j,ispec_PML)  + A5 * rmemory_duz_dx_prime(i,j,ispec_PML)
                    duz_dzl_prime = PML_duz_dzl(i,j,ispec_PML)  + A5 * rmemory_duz_dz_prime(i,j,ispec_PML)
                    else
                    dux_dzl = PML_dux_dzl(i,j,ispec_PML)  + A5 * rmemory_dux_dz(i,j,ispec_PML)
                    duz_dzl = PML_duz_dzl(i,j,ispec_PML)  + A5 * rmemory_duz_dz(i,j,ispec_PML)
                    endif



                  elseif(region_CPML(ispec) == CPML_TOP .or. region_CPML(ispec) == CPML_BOTTOM) then
!------------------------------------------------------------------------------
!---------------------------- TOP & BOTTOM ------------------------------------
!------------------------------------------------------------------------------
                    !---------------------- A7 --------------------------
                    A7 = d_z_store(i,j,ispec_PML)
                    bb = alpha_z_store(i,j,ispec_PML)

                    if(stage_time_scheme == 1) then 
                    coef0 = exp(- bb * deltat)

                    if ( abs( bb ) > 1.d-3) then
                      coef1 = (1.0d0 - exp(- bb * deltat / 2.0d0) ) / bb
                      coef2 = (1.0d0 - exp(- bb * deltat / 2.0d0) ) * exp(- bb * deltat / 2.0d0) / bb
                    else
                      coef1 = deltat / 2.0d0
                      coef2 = deltat / 2.0d0
                    end if

                    if(ROTATE_PML_ACTIVATE)then
                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dx(i,j,ispec_PML) = coef0 * rmemory_dux_dx(i,j,ispec_PML) &
                    + PML_dux_dxl_new(i,j,ispec_PML) * coef1 + PML_dux_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dz(i,j,ispec_PML) = coef0 * rmemory_dux_dz(i,j,ispec_PML) &
                    + PML_dux_dzl_new(i,j,ispec_PML) * coef1 + PML_dux_dzl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dx(i,j,ispec_PML) = coef0 * rmemory_duz_dx(i,j,ispec_PML) &
                    + PML_duz_dxl_new(i,j,ispec_PML) * coef1 + PML_duz_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dz(i,j,ispec_PML) = coef0 * rmemory_duz_dz(i,j,ispec_PML) &
                    + PML_duz_dzl_new(i,j,ispec_PML) * coef1 + PML_duz_dzl(i,j,ispec_PML) * coef2
                    else
                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dx(i,j,ispec_PML) = coef0*rmemory_dux_dx(i,j,ispec_PML) &
                    + PML_dux_dxl_new(i,j,ispec_PML) * coef1 + PML_dux_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dx(i,j,ispec_PML) = coef0*rmemory_duz_dx(i,j,ispec_PML) &
                    + PML_duz_dxl_new(i,j,ispec_PML) * coef1 + PML_duz_dxl(i,j,ispec_PML) * coef2
                    endif

                    end if  

                    if(stage_time_scheme == 6) then   
                     rmemory_dux_dx_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_dux_dx(i,j,ispec_PML) + PML_dux_dxl(i,j,ispec_PML))
                     rmemory_dux_dx(i,j,ispec_PML) = rmemory_dux_dx(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_dux_dx_LDDRK(i,j,ispec_PML)

                     rmemory_duz_dx_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_duz_dx(i,j,ispec_PML) + PML_duz_dxl(i,j,ispec_PML))
                     rmemory_duz_dx(i,j,ispec_PML) = rmemory_duz_dx(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_duz_dx_LDDRK(i,j,ispec_PML)
                    end if  

                    if(ROTATE_PML_ACTIVATE)then
                    dux_dxl = PML_dux_dxl(i,j,ispec_PML)  + A7 * rmemory_dux_dx(i,j,ispec_PML)
                    dux_dzl = PML_dux_dzl(i,j,ispec_PML)  + A7 * rmemory_dux_dz(i,j,ispec_PML)
                    duz_dxl = PML_duz_dxl(i,j,ispec_PML)  + A7 * rmemory_duz_dx(i,j,ispec_PML)
                    duz_dzl = PML_duz_dzl(i,j,ispec_PML)  + A7 * rmemory_duz_dz(i,j,ispec_PML)
                    else
                    dux_dxl = dux_dxl  + A7 * rmemory_dux_dx(i,j,ispec_PML)
                    duz_dxl = duz_dxl  + A7 * rmemory_duz_dx(i,j,ispec_PML)
                    endif



                    !---------------------- A6 --------------------------
                    A6 = - d_z_store(i,j,ispec_PML) / ( k_z_store(i,j,ispec_PML) ** 2 )
                    bb = d_z_store(i,j,ispec_PML) / k_z_store(i,j,ispec_PML) + alpha_z_store(i,j,ispec_PML)

                    if(stage_time_scheme == 1) then 
                    coef0 = exp(-bb * deltat)

                    if ( abs(bb) > 1.d-3 ) then
                      coef1 = (1.d0 - exp(-bb * deltat / 2.d0)) / bb
                      coef2 = (1.d0 - exp(-bb* deltat / 2.d0)) * exp(-bb * deltat / 2.d0) / bb
                    else
                      coef1 = deltat / 2.0d0
                      coef2 = deltat / 2.0d0
                    end if

                    if(ROTATE_PML_ACTIVATE)then
                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dx_prime(i,j,ispec_PML) = coef0 * rmemory_dux_dx_prime(i,j,ispec_PML) &
                    + PML_dux_dxl_new(i,j,ispec_PML) *coef1 + PML_dux_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dz_prime(i,j,ispec_PML) = coef0 * rmemory_dux_dz_prime(i,j,ispec_PML) &
                    + PML_dux_dzl_new(i,j,ispec_PML) *coef1 + PML_dux_dzl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dx_prime(i,j,ispec_PML) = coef0 * rmemory_duz_dx_prime(i,j,ispec_PML) &
                    + PML_duz_dxl_new(i,j,ispec_PML) *coef1 + PML_duz_dxl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dz_prime(i,j,ispec_PML) = coef0 * rmemory_duz_dz_prime(i,j,ispec_PML) &
                    + PML_duz_dzl_new(i,j,ispec_PML) *coef1 + PML_duz_dzl(i,j,ispec_PML) * coef2
                    else
                    !! DK DK new from Wang eq (21)
                    rmemory_dux_dz(i,j,ispec_PML) = coef0 * rmemory_dux_dz(i,j,ispec_PML) &
                    + PML_dux_dzl_new(i,j,ispec_PML) *coef1 + PML_dux_dzl(i,j,ispec_PML) * coef2

                    !! DK DK new from Wang eq (21)
                    rmemory_duz_dz(i,j,ispec_PML) = coef0 * rmemory_duz_dz(i,j,ispec_PML) &
                    + PML_duz_dzl_new(i,j,ispec_PML) *coef1 + PML_duz_dzl(i,j,ispec_PML) * coef2
                    endif

                    end if

                    if(stage_time_scheme == 6) then 
                     rmemory_dux_dz_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_dux_dz(i,j,ispec_PML) + PML_dux_dzl(i,j,ispec_PML))
                     rmemory_dux_dz(i,j,ispec_PML) = rmemory_dux_dz(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_dux_dz_LDDRK(i,j,ispec_PML)

                     rmemory_duz_dz_LDDRK(i,j,ispec_PML) = alpha_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_duz_dz(i,j,ispec_PML) + PML_duz_dzl(i,j,ispec_PML))
                     rmemory_duz_dz(i,j,ispec_PML) = rmemory_duz_dz(i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_duz_dz_LDDRK(i,j,ispec_PML)
                    end if 

                    if(ROTATE_PML_ACTIVATE)then
                    dux_dxl_prime = PML_dux_dxl(i,j,ispec_PML)  + A6 * rmemory_dux_dx_prime(i,j,ispec_PML)
                    dux_dzl_prime = PML_dux_dzl(i,j,ispec_PML)  + A6 * rmemory_dux_dz_prime(i,j,ispec_PML)
                    duz_dxl_prime = PML_duz_dxl(i,j,ispec_PML)  + A6 * rmemory_duz_dx_prime(i,j,ispec_PML)
                    duz_dzl_prime = PML_duz_dzl(i,j,ispec_PML)  + A6 * rmemory_duz_dz_prime(i,j,ispec_PML)
                    else
                    dux_dzl = dux_dzl  + A6 * rmemory_dux_dz(i,j,ispec_PML)
                    duz_dzl = duz_dzl  + A6 * rmemory_duz_dz(i,j,ispec_PML)
                    endif



                  endif
              endif ! PML_BOUNDARY_CONDITIONS

              if(SIMULATION_TYPE == 2) then ! Adjoint calculation, backward wavefield
                 b_dux_dxl = b_dux_dxi*xixl + b_dux_dgamma*gammaxl
                 b_dux_dzl = b_dux_dxi*xizl + b_dux_dgamma*gammazl

                 b_duy_dxl = b_duy_dxi*xixl + b_duy_dgamma*gammaxl
                 b_duy_dzl = b_duy_dxi*xizl + b_duy_dgamma*gammazl

                 b_duz_dxl = b_duz_dxi*xixl + b_duz_dgamma*gammaxl
                 b_duz_dzl = b_duz_dxi*xizl + b_duz_dgamma*gammazl
              endif

              ! compute stress tensor (include attenuation or anisotropy if needed)

              if(ATTENUATION_VISCOELASTIC_SOLID) then

                 ! attenuation is implemented following the memory variable formulation of
                 ! J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
                 ! vol. 58(1), p. 110-120 (1993). More details can be found in
                 ! J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation in a linear
                 ! viscoelastic medium, Geophysical Journal International, vol. 95, p. 597-611 (1988).

                 ! When implementing viscoelasticity according to Carcione 1993 paper, the attenuation is
                 ! non-causal rather than causal. We fixed the problem by using equations in Carcione's
                 ! 2004 paper and his 2007 book.

                 !J. M. Carcione, H B. Helle, The physics and simulation of wave propagation at the ocean
                 ! bottom, Geophysics, vol. 69(3), p. 825-839, 2004
                 !J. M. Carcione, Wave fields in real media: wave propagation in anisotropic, anelastic
                 ! and porous media, Elsevier, p. 124-125, 2007

                 ! compute unrelaxed elastic coefficients from formulas in Carcione 2007 page 125
                 lambdal_relaxed_viscoelastic = (lambdal_unrelaxed_elastic + mul_unrelaxed_elastic) / Mu_nu1(i,j,ispec) &
                                                - mul_unrelaxed_elastic / Mu_nu2(i,j,ispec)
                 mul_relaxed_viscoelastic = mul_unrelaxed_elastic / Mu_nu2(i,j,ispec)
                 lambdalplus2mul_relaxed_viscoel = lambdal_relaxed_viscoelastic + TWO*mul_relaxed_viscoelastic

                 ! compute the stress using the unrelaxed Lame parameters (Carcione 2007 page 125)
                 sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
                 sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                 sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl

                 ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
                 ! beware: there is a bug in Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below
                 e1_sum = 0._CUSTOM_REAL
                 e11_sum = 0._CUSTOM_REAL
                 e13_sum = 0._CUSTOM_REAL

                 do i_sls = 1,N_SLS
                    e1_sum = e1_sum + e1(i,j,ispec,i_sls)
                    e11_sum = e11_sum + e11(i,j,ispec,i_sls)
                    e13_sum = e13_sum + e13(i,j,ispec,i_sls)
                 enddo

                 sigma_xx = sigma_xx + (lambdal_relaxed_viscoelastic + mul_relaxed_viscoelastic) * e1_sum &
                                     + TWO * mul_relaxed_viscoelastic * e11_sum
                 sigma_xz = sigma_xz + mul_relaxed_viscoelastic * e13_sum
                 sigma_zz = sigma_zz + (lambdal_relaxed_viscoelastic + mul_relaxed_viscoelastic) * e1_sum &
                                     - TWO * mul_relaxed_viscoelastic * e11_sum
                 sigma_zx = sigma_xz

                 if(PML_BOUNDARY_CONDITIONS .and. is_PML(ispec)) then
                     ispec_PML=spec_to_PML(ispec)
                     sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*PML_duz_dzl(i,j,ispec_PML)
                     sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*PML_dux_dxl(i,j,ispec_PML)
                     sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j,ispec_PML) + dux_dzl)
                     sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j,ispec_PML) + duz_dxl)
                 endif

              else

                 ! no attenuation
                 sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
                 sigma_xy = mul_unrelaxed_elastic*duy_dxl
                 sigma_xz = mul_unrelaxed_elastic*(duz_dxl + dux_dzl)
                 sigma_zy = mul_unrelaxed_elastic*duy_dzl
                 sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
                 sigma_zx = sigma_xz

!! DK DK added this for Guenneau, March 2012
#ifdef USE_GUENNEAU
  include "include_stiffness_Guenneau.f90"
#endif
!! DK DK added this for Guenneau, March 2012

                 if(PML_BOUNDARY_CONDITIONS .and. is_PML(ispec)) then
                     ispec_PML=spec_to_PML(ispec)
                     if(ROTATE_PML_ACTIVATE)then
                     theta = -ROTATE_PML_ANGLE/180.d0*Pi
                     if(it==1)write(*,*)theta,ROTATE_PML_ACTIVATE,cos(theta),sin(theta)
                     ct=cos(theta)
                     st=sin(theta)
                     sigma_xx_prime = lambdaplus2mu_unrelaxed_elastic*(ct**2*dux_dxl+ct*st*duz_dxl+ct*st*dux_dzl+st**2*duz_dzl) &
                                      + lambdal_unrelaxed_elastic*(st**2*PML_dux_dxl(i,j,ispec_PML)&
                                                                   -ct*st*PML_duz_dxl(i,j,ispec_PML)&
                                                                   -ct*st*PML_dux_dzl(i,j,ispec_PML)&
                                                                   +ct**2*PML_duz_dzl(i,j,ispec_PML))

                     sigma_xz_prime = mul_unrelaxed_elastic * (-ct*st*dux_dxl+ct**2*duz_dxl-st**2*dux_dzl+ct*st*duz_dzl) &
                                      +mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j,ispec_PML)&
                                                                   -st**2*PML_duz_dxl(i,j,ispec_PML)&
                                                                   +ct**2*PML_dux_dzl(i,j,ispec_PML)&
                                                                   +ct*st*PML_duz_dzl(i,j,ispec_PML))

                     sigma_zx_prime = mul_unrelaxed_elastic * (-ct*st*PML_dux_dxl(i,j,ispec_PML)&
                                                                   +ct**2*PML_duz_dxl(i,j,ispec_PML)&
                                                                   -st**2*PML_dux_dzl(i,j,ispec_PML)&
                                                                   +ct*st*PML_duz_dzl(i,j,ispec_PML)) &
                                      +mul_unrelaxed_elastic * (-ct*st*dux_dxl_prime-st**2*duz_dxl_prime &
                                                                 +ct**2*dux_dzl_prime+ct*st*duz_dzl_prime)

                     sigma_zz_prime = lambdaplus2mu_unrelaxed_elastic*(st**2*dux_dxl_prime-ct*st*duz_dxl_prime&
                                                                       -ct*st*dux_dzl_prime+ct**2*duz_dzl_prime) &
                                      + lambdal_unrelaxed_elastic*(ct**2*PML_dux_dxl(i,j,ispec_PML)&
                                                                   +ct*st*PML_duz_dxl(i,j,ispec_PML)&
                                                                   +ct*st*PML_dux_dzl(i,j,ispec_PML)&
                                                                   +st**2*PML_duz_dzl(i,j,ispec_PML))

                     sigma_xx = ct**2*sigma_xx_prime-ct*st*sigma_xz_prime-ct*st*sigma_zx_prime+st**2*sigma_zz_prime
                     sigma_xz = ct*st*sigma_xx_prime+ct**2*sigma_xz_prime-st**2*sigma_zx_prime-ct*st*sigma_zz_prime
                     sigma_zx = ct*st*sigma_xx_prime-st**2*sigma_xz_prime+ct**2*sigma_zx_prime-ct*st*sigma_zz_prime
                     sigma_zz = st**2*sigma_xx_prime+ct*st*sigma_xz_prime+ct*st*sigma_zx_prime+ct**2*sigma_zz_prime
                     else
                     sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*PML_duz_dzl(i,j,ispec_PML)
                     sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*PML_dux_dxl(i,j,ispec_PML)
                     sigma_zx = mul_unrelaxed_elastic * (PML_duz_dxl(i,j,ispec_PML) + dux_dzl)
                     sigma_xz = mul_unrelaxed_elastic * (PML_dux_dzl(i,j,ispec_PML) + duz_dxl)
                     endif
                 endif

                 if(SIMULATION_TYPE == 2) then ! Adjoint calculation, backward wavefield
                    b_sigma_xx = lambdaplus2mu_unrelaxed_elastic*b_dux_dxl + lambdal_unrelaxed_elastic*b_duz_dzl
                    b_sigma_xy = mul_unrelaxed_elastic*b_duy_dxl
                    b_sigma_xz = mul_unrelaxed_elastic*(b_duz_dxl + b_dux_dzl)
                    b_sigma_zy = mul_unrelaxed_elastic*b_duy_dzl
                    b_sigma_zz = lambdaplus2mu_unrelaxed_elastic*b_duz_dzl + lambdal_unrelaxed_elastic*b_dux_dxl

                    ! the stress tensor is symmetric
                    b_sigma_zx = b_sigma_xz
                 endif

              endif

              ! full anisotropy
              if(anisotropic(ispec)) then
                 if(assign_external_model) then
                    c11 = c11ext(i,j,ispec)
                    c13 = c13ext(i,j,ispec)
                    c15 = c15ext(i,j,ispec)
                    c33 = c33ext(i,j,ispec)
                    c35 = c35ext(i,j,ispec)
                    c55 = c55ext(i,j,ispec)
                 else
                    c11 = anisotropy(1,kmato(ispec))
                    c13 = anisotropy(2,kmato(ispec))
                    c15 = anisotropy(3,kmato(ispec))
                    c33 = anisotropy(4,kmato(ispec))
                    c35 = anisotropy(5,kmato(ispec))
                    c55 = anisotropy(6,kmato(ispec))
                 end if

                 ! implement anisotropy in 2D
                 sigma_xx = c11*dux_dxl + c15*(duz_dxl + dux_dzl) + c13*duz_dzl
                 sigma_zz = c13*dux_dxl + c35*(duz_dxl + dux_dzl) + c33*duz_dzl
                 sigma_xz = c15*dux_dxl + c55*(duz_dxl + dux_dzl) + c35*duz_dzl

              endif

              ! Pre-kernels calculation
              if(SIMULATION_TYPE == 2) then
                 iglob = ibool(i,j,ispec)
                 if(p_sv)then !P-SV waves
                    dsxx =  dux_dxl
                    dsxz = HALF * (duz_dxl + dux_dzl)
                    dszz =  duz_dzl

                    b_dsxx =  b_dux_dxl
                    b_dsxz = HALF * (b_duz_dxl + b_dux_dzl)
                    b_dszz =  b_duz_dzl

                    kappa_k(iglob) = (dux_dxl + duz_dzl) *  (b_dux_dxl + b_duz_dzl)
                    mu_k(iglob) = dsxx * b_dsxx + dszz * b_dszz + &
                         2._CUSTOM_REAL * dsxz * b_dsxz - 1._CUSTOM_REAL/3._CUSTOM_REAL * kappa_k(iglob)
                 else !SH (membrane) waves
                    mu_k(iglob) = duy_dxl * b_duy_dxl + duy_dzl * b_duy_dzl
                 endif
              endif

              jacobianl = jacobian(i,j,ispec)

! the stress tensor is symmetric by default, unless defined otherwise
#ifdef USE_GUENNEAU
              sigma_zx = sigma_xz
#endif

              ! weak formulation term based on stress tensor (non-symmetric form)
              ! also add GLL integration weights
              tempx1(i,j) = wzgll(j)*jacobianl*(sigma_xx*xixl+sigma_zx*xizl) ! this goes to accel_x
              tempy1(i,j) = wzgll(j)*jacobianl*(sigma_xy*xixl+sigma_zy*xizl) ! this goes to accel_y
              tempz1(i,j) = wzgll(j)*jacobianl*(sigma_xz*xixl+sigma_zz*xizl) ! this goes to accel_z

              tempx2(i,j) = wxgll(i)*jacobianl*(sigma_xx*gammaxl+sigma_zx*gammazl) ! this goes to accel_x
              tempy2(i,j) = wxgll(i)*jacobianl*(sigma_xy*gammaxl+sigma_zy*gammazl) ! this goes to accel_y
              tempz2(i,j) = wxgll(i)*jacobianl*(sigma_xz*gammaxl+sigma_zz*gammazl) ! this goes to accel_z

              if(SIMULATION_TYPE == 2) then ! Adjoint calculation, backward wavefield
                 b_tempx1(i,j) = wzgll(j)*jacobianl*(b_sigma_xx*xixl+b_sigma_zx*xizl) ! this goes to b_accel_x
                 b_tempy1(i,j) = wzgll(j)*jacobianl*(b_sigma_xy*xixl+b_sigma_zy*xizl) ! this goes to b_accel_y
                 b_tempz1(i,j) = wzgll(j)*jacobianl*(b_sigma_xz*xixl+b_sigma_zz*xizl) ! this goes to b_accel_z

                 b_tempx2(i,j) = wxgll(i)*jacobianl*(b_sigma_xx*gammaxl+b_sigma_zx*gammazl) ! this goes to b_accel_x
                 b_tempy2(i,j) = wxgll(i)*jacobianl*(b_sigma_xy*gammaxl+b_sigma_zy*gammazl) ! this goes to b_accel_y
                 b_tempz2(i,j) = wxgll(i)*jacobianl*(b_sigma_xz*gammaxl+b_sigma_zz*gammazl) ! this goes to b_accel_z
              endif

           enddo
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! PML
        if( is_PML(ispec) .and. PML_BOUNDARY_CONDITIONS ) then
            ispec_PML=spec_to_PML(ispec)
            do j = 1,NGLLZ
              do i = 1,NGLLX
                if ( assign_external_model) then
                 rhol = rhoext(i,j,ispec)
                else
                 rhol = density(1,kmato(ispec))
                endif
                    iglob=ibool(i,j,ispec)

!------------------------------------------------------------------------------
!---------------------------- LEFT & RIGHT ------------------------------------
!------------------------------------------------------------------------------
                   if (region_CPML(ispec) == CPML_LEFT .or. region_CPML(ispec) == CPML_RIGHT) then

                    bb = alpha_x_store(i,j,ispec_PML)
                    if(stage_time_scheme == 1) then  
                    coef0 = exp(- bb * deltat)

                    if ( abs( bb ) > 1.d-3) then
                       coef1 = (1.0d0 - exp(- bb * deltat / 2.0d0) ) / bb
                       coef2 = (1.0d0 - exp(- bb * deltat / 2.0d0) ) * exp(- bb * deltat / 2.0d0) / bb
                    else
                       coef1 = deltat / 2.0d0
                       coef2 = deltat / 2.0d0
                    end if

                    rmemory_displ_elastic(1,1,i,j,ispec_PML)=coef0 * rmemory_displ_elastic(1,1,i,j,ispec_PML) &
                    + displ_elastic_new(1,iglob) * coef1 + displ_elastic(1,iglob) * coef2
                    rmemory_displ_elastic(2,1,i,j,ispec_PML)=0.0

                    rmemory_displ_elastic(1,3,i,j,ispec_PML)=coef0 * rmemory_displ_elastic(1,3,i,j,ispec_PML) &
                    + displ_elastic_new(3,iglob) * coef1 + displ_elastic(3,iglob) * coef2
                    rmemory_displ_elastic(2,3,i,j,ispec_PML)=0.0
                    end if

!------------------------------------------------------------------------------
!-------------------------------- CORNER --------------------------------------
!------------------------------------------------------------------------------
                   elseif (region_CPML(ispec) == CPML_TOP_LEFT .or. region_CPML(ispec) == CPML_TOP_RIGHT .or. &
                           region_CPML(ispec) == CPML_BOTTOM_LEFT .or. region_CPML(ispec) == CPML_BOTTOM_RIGHT) then

                       !------------------------------------------------------------
                       bb = alpha_x_store(i,j,ispec_PML)

                      if(stage_time_scheme == 1) then 
                       coef0 = exp(- bb * deltat)

                       if ( abs(bb) > 1.d-3 ) then
                          coef1 = ( 1 - exp(- bb * deltat / 2) ) / bb
                          coef2 = ( 1 - exp(- bb * deltat / 2) ) * exp(- bb * deltat / 2) / bb
                       else
                          coef1 = deltat / 2.0d0
                          coef2 = deltat / 2.0d0
                       end if

                       rmemory_displ_elastic(1,1,i,j,ispec_PML)= &
                        coef0 * rmemory_displ_elastic(1,1,i,j,ispec_PML) &
                        + displ_elastic_new(1,iglob) * coef1 + displ_elastic(1,iglob) * coef2
                       rmemory_displ_elastic(1,3,i,j,ispec_PML)= &
                        coef0 * rmemory_displ_elastic(1,3,i,j,ispec_PML) &
                        + displ_elastic_new(3,iglob) * coef1 + displ_elastic(3,iglob) * coef2
                       end if

                       !------------------------------------------------------------
                       bb = alpha_z_store(i,j,ispec_PML)

                      if(stage_time_scheme == 1) then 
                       coef0 = exp(- bb * deltat)

                       if ( abs(bb) > 1.d-3 ) then
                          coef1 = ( 1 - exp(- bb * deltat / 2) ) / bb
                          coef2 = ( 1 - exp(- bb * deltat / 2) ) * exp(- bb * deltat / 2) / bb
                       else
                          coef1 = deltat / 2.0d0
                          coef2 = deltat / 2.0d0
                       end if

                       rmemory_displ_elastic(2,1,i,j,ispec_PML)=coef0 * rmemory_displ_elastic(2,1,i,j,ispec_PML) &
                        + displ_elastic_new(1,iglob)*(it+0.5)*deltat * coef1 &
                        + displ_elastic(1,iglob)*(it-0.5)*deltat * coef2
                       rmemory_displ_elastic(2,3,i,j,ispec_PML)=coef0 * rmemory_displ_elastic(2,3,i,j,ispec_PML) &
                        + displ_elastic_new(3,iglob)*(it+0.5)*deltat * coef1 &
                        + displ_elastic(3,iglob)*(it-0.5)*deltat * coef2
                       end if

                  elseif(region_CPML(ispec) == CPML_TOP .or. region_CPML(ispec) == CPML_BOTTOM) then
!------------------------------------------------------------------------------
!-------------------------------- TOP & BOTTOM --------------------------------
!------------------------------------------------------------------------------

                      !------------------------------------------------------------
                      bb = alpha_z_store(i,j,ispec_PML)

                      if(stage_time_scheme == 1) then
                      coef0 = exp(- bb * deltat)

                      if ( abs( bb ) > 1.d-3) then
                         coef1 = (1.0d0 - exp(- bb * deltat / 2.0d0) ) / bb
                         coef2 = (1.0d0 - exp(- bb * deltat / 2.0d0) ) * exp(- bb * deltat / 2.0d0) / bb
                      else
                         coef1 = deltat / 2.0d0
                         coef2 = deltat / 2.0d0
                      end if

                      rmemory_displ_elastic(1,1,i,j,ispec_PML)=0.d0
                      rmemory_displ_elastic(2,1,i,j,ispec_PML)=coef0 * rmemory_displ_elastic(2,1,i,j,ispec_PML) &
                      + displ_elastic_new(1,iglob) * coef1 + displ_elastic(1,iglob) * coef2

                      rmemory_displ_elastic(1,3,i,j,ispec_PML)=0.d0
                      rmemory_displ_elastic(2,3,i,j,ispec_PML)=coef0 * rmemory_displ_elastic(2,3,i,j,ispec_PML) &
                      + displ_elastic_new(3,iglob) * coef1 + displ_elastic(3,iglob) * coef2
                      end if

                endif


                   if (region_CPML(ispec) == CPML_LEFT .or. region_CPML(ispec) == CPML_RIGHT) then

                     A0 = - alpha_x_store(i,j,ispec_PML) * d_x_store(i,j,ispec_PML)
                     A1 = d_x_store(i,j,ispec_PML)
                     A2 = k_x_store(i,j,ispec_PML)
                     A3 = d_x_store(i,j,ispec_PML) * alpha_x_store(i,j,ispec_PML) ** 2
                     A4 = 0.d0

                    if(stage_time_scheme == 6) then

                     bb = alpha_x_store(i,j,ispec_PML)

                     rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML) = &
                     alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_displ_elastic(1,1,i,j,ispec_PML) + displ_elastic(1,iglob))
                     rmemory_displ_elastic(1,1,i,j,ispec_PML) = rmemory_displ_elastic(1,1,i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML)
                     rmemory_displ_elastic(2,1,i,j,ispec_PML) =0.d0

                     rmemory_displ_elastic_LDDRK(1,3,i,j,ispec_PML) = &
                     alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,3,i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_displ_elastic(1,3,i,j,ispec_PML) + displ_elastic(3,iglob))

                     rmemory_displ_elastic(1,3,i,j,ispec_PML) = rmemory_displ_elastic(1,3,i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,3,i,j,ispec_PML)
                     rmemory_displ_elastic(2,3,i,j,ispec_PML) =0.d0

                    end if 

                     accel_elastic_PML(1,i,j,ispec_PML)= wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * &
                          ( &
                          A0 * displ_elastic(1,iglob) + &
                          A1 *veloc_elastic(1,iglob)  + &
                          A3 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + &
                          A4 * rmemory_displ_elastic(2,1,i,j,ispec_PML)   &
                          )
                     accel_elastic_PML(3,i,j,ispec_PML)= wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * &
                          ( &
                          A0 * displ_elastic(3,iglob) + &
                          A1 * veloc_elastic(3,iglob)  + &
                          A3 * rmemory_displ_elastic(1,3,i,j,ispec_PML) + &
                          A4 * rmemory_displ_elastic(2,3,i,j,ispec_PML)   &
                          )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!corner
                   elseif (region_CPML(ispec) == CPML_TOP_LEFT .or. region_CPML(ispec) == CPML_TOP_RIGHT .or. &
                           region_CPML(ispec) == CPML_BOTTOM_LEFT .or. region_CPML(ispec) == CPML_BOTTOM_RIGHT) then

                     A0 = d_x_store(i,j,ispec_PML) * d_z_store(i,j,ispec_PML) &
                        - alpha_x_store(i,j,ispec_PML) * d_x_store(i,j,ispec_PML) * k_z_store(i,j,ispec_PML) &
                        - alpha_z_store(i,j,ispec_PML) * d_z_store(i,j,ispec_PML) * k_x_store(i,j,ispec_PML)

                     A1 = d_x_store(i,j,ispec_PML) * k_z_store(i,j,ispec_PML) + d_z_store(i,j,ispec_PML) * k_x_store(i,j,ispec_PML)

                     A2 = k_x_store(i,j,ispec_PML) * k_z_store(i,j,ispec_PML)

                     A3 = alpha_x_store(i,j,ispec_PML) ** 2*(d_x_store(i,j,ispec_PML) * k_z_store(i,j,ispec_PML)+ &
                            d_z_store(i,j,ispec_PML) * k_x_store(i,j,ispec_PML)) &
                            -2.d0 * alpha_x_store(i,j,ispec_PML)*d_x_store(i,j,ispec_PML)*d_z_store(i,j,ispec_PML)+ &
                            (it+0.5)*deltat*alpha_x_store(i,j,ispec_PML)**2*d_x_store(i,j,ispec_PML)*d_z_store(i,j,ispec_PML)
                     A4 = -alpha_x_store(i,j,ispec_PML) ** 2*d_x_store(i,j,ispec_PML)*d_z_store(i,j,ispec_PML)

                    if(stage_time_scheme == 6) then 
                     A3 = alpha_x_store(i,j,ispec_PML) ** 2*(d_x_store(i,j,ispec_PML) * k_z_store(i,j,ispec_PML)+ &
                            d_z_store(i,j,ispec_PML) * k_x_store(i,j,ispec_PML)) &
                            -2.d0 * alpha_x_store(i,j,ispec_PML)*d_x_store(i,j,ispec_PML)*d_z_store(i,j,ispec_PML)
                     A4 = alpha_x_store(i,j,ispec_PML) ** 2*d_x_store(i,j,ispec_PML)*d_z_store(i,j,ispec_PML)
                    end if   

                    if(stage_time_scheme == 6) then  

                     bb = alpha_z_store(i,j,ispec_PML)

                     rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML) = &
                     alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_displ_elastic(1,1,i,j,ispec_PML) + displ_elastic(1,iglob))
                     rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML) = &
                     alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_displ_elastic(2,1,i,j,ispec_PML) + rmemory_displ_elastic(1,1,i,j,ispec_PML))

                     rmemory_displ_elastic(1,1,i,j,ispec_PML) = rmemory_displ_elastic(1,1,i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,1,i,j,ispec_PML)
                     rmemory_displ_elastic(2,1,i,j,ispec_PML) = rmemory_displ_elastic(2,1,i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML)

                     rmemory_displ_elastic_LDDRK(1,3,i,j,ispec_PML) = &
                     alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,3,i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_displ_elastic(1,3,i,j,ispec_PML) + displ_elastic(3,iglob))
                     rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML) = &
                     alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_displ_elastic(2,3,i,j,ispec_PML) + rmemory_displ_elastic(1,3,i,j,ispec_PML))

                     rmemory_displ_elastic(1,3,i,j,ispec_PML) = rmemory_displ_elastic(1,3,i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(1,3,i,j,ispec_PML)
                     rmemory_displ_elastic(2,3,i,j,ispec_PML) = rmemory_displ_elastic(2,3,i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML)

                    end if

                     accel_elastic_PML(1,i,j,ispec_PML)= wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * &
                          ( &
                          A0 * displ_elastic(1,iglob) + &
                          A1 *veloc_elastic(1,iglob)  + &
                          A3 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + &
                          A4 * rmemory_displ_elastic(2,1,i,j,ispec_PML)   &
                           )
                     accel_elastic_PML(3,i,j,ispec_PML)= wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * &
                          ( &
                          A0 * displ_elastic(3,iglob) + &
                          A1 *veloc_elastic(3,iglob)  + &
                          A3 * rmemory_displ_elastic(1,3,i,j,ispec_PML) + &
                          A4 * rmemory_displ_elastic(2,3,i,j,ispec_PML)   &
                          )



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!corner
               elseif(region_CPML(ispec) == CPML_TOP .or. region_CPML(ispec) == CPML_BOTTOM) then

                     A0 = - alpha_z_store(i,j,ispec_PML) * d_z_store(i,j,ispec_PML)
                     A1 = d_z_store(i,j,ispec_PML)
                     A2 = k_z_store(i,j,ispec_PML)
                     A3 = 0.d0
                     A4 = d_z_store(i,j,ispec_PML) * alpha_z_store(i,j,ispec_PML) ** 2

                    if(stage_time_scheme == 6) then 

                     bb = alpha_z_store(i,j,ispec_PML)

                     rmemory_displ_elastic(1,1,i,j,ispec_PML) =0.d0
                     rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML) = &
                     alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_displ_elastic(2,1,i,j,ispec_PML) + displ_elastic(1,iglob))
                     rmemory_displ_elastic(2,1,i,j,ispec_PML) = rmemory_displ_elastic(2,1,i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,1,i,j,ispec_PML)

                     rmemory_displ_elastic(1,3,i,j,ispec_PML) =0.d0
                     rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML) = &
                     alpha_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML) &
                     + deltat * (-bb * rmemory_displ_elastic(2,3,i,j,ispec_PML) + displ_elastic(3,iglob))
                     rmemory_displ_elastic(2,3,i,j,ispec_PML) = rmemory_displ_elastic(2,3,i,j,ispec_PML) + &
                     beta_LDDRK(i_stage) * rmemory_displ_elastic_LDDRK(2,3,i,j,ispec_PML)

                    end if

                     accel_elastic_PML(1,i,j,ispec_PML)= wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * &
                          ( &
                          A0 * displ_elastic(1,iglob) + &
                          A1 *veloc_elastic(1,iglob)  + &
                          A3 * rmemory_displ_elastic(1,1,i,j,ispec_PML) + &
                          A4 * rmemory_displ_elastic(2,1,i,j,ispec_PML)   &
                          )
                     accel_elastic_PML(3,i,j,ispec_PML)= wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec) * &
                          ( &
                          A0 * displ_elastic(3,iglob) + &
                          A1 *veloc_elastic(3,iglob)  + &
                          A3 * rmemory_displ_elastic(1,3,i,j,ispec_PML) + &
                          A4 * rmemory_displ_elastic(2,3,i,j,ispec_PML)   &
                          )



               endif

           enddo
        enddo

     endif ! PML_BOUNDARY_CONDITIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
                 accel_elastic(1,iglob) = accel_elastic(1,iglob) &
                  - (tempx1(k,j)*hprimewgll_xx(k,i) + tempx2(i,k)*hprimewgll_zz(k,j))
                 accel_elastic(2,iglob) = accel_elastic(2,iglob) &
                  - (tempy1(k,j)*hprimewgll_xx(k,i) + tempy2(i,k)*hprimewgll_zz(k,j))
                 accel_elastic(3,iglob) = accel_elastic(3,iglob) &
                  - (tempz1(k,j)*hprimewgll_xx(k,i) + tempz2(i,k)*hprimewgll_zz(k,j))

                 if(SIMULATION_TYPE == 2) then ! Adjoint calculation, backward wavefield
                    b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - &
                         (b_tempx1(k,j)*hprimewgll_xx(k,i) + b_tempx2(i,k)*hprimewgll_zz(k,j))
                    b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - &
                         (b_tempy1(k,j)*hprimewgll_xx(k,i) + b_tempy2(i,k)*hprimewgll_zz(k,j))
                    b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) - &
                         (b_tempz1(k,j)*hprimewgll_xx(k,i) + b_tempz2(i,k)*hprimewgll_zz(k,j))
                 endif
              enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            if(is_PML(ispec) .and. PML_BOUNDARY_CONDITIONS)then
                ispec_PML=spec_to_PML(ispec)
                      accel_elastic(1,iglob) = accel_elastic(1,iglob) - accel_elastic_PML(1,i,j,ispec_PML)
                      accel_elastic(2,iglob) = accel_elastic(2,iglob)
                      accel_elastic(3,iglob) = accel_elastic(3,iglob) - accel_elastic_PML(3,i,j,ispec_PML)
            endif ! PML_BOUNDARY_CONDITIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           enddo ! second loop over the GLL points
        enddo

     endif ! end of test if elastic element

  enddo ! end of loop over all spectral elements


  !
  !--- absorbing boundaries
  !
  if(anyabs .and. .not. PML_BOUNDARY_CONDITIONS) then

     count_left=1
     count_right=1
     count_bottom=1

     do ispecabs = 1,nelemabs

        ispec = numabs(ispecabs)

        ! get elastic parameters of current spectral element
        lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
        mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
        rhol  = density(1,kmato(ispec))
        kappal  = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic/3._CUSTOM_REAL
        cpl = sqrt((kappal + 4._CUSTOM_REAL*mul_unrelaxed_elastic/3._CUSTOM_REAL)/rhol)
        csl = sqrt(mul_unrelaxed_elastic/rhol)

        !--- left absorbing boundary
        if(codeabs(IEDGE4,ispecabs)) then

           i = 1

           do j = 1,NGLLZ

              iglob = ibool(i,j,ispec)

              ! for analytical initial plane wave for Bielak's conditions
              ! left or right edge, horizontal normal vector
              if(add_Bielak_conditions .and. initialfield) then
                 if (.not.over_critical_angle) then
                    call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                         x0_source, z0_source, A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                         c_inc, c_refl, time_offset,f0)
                    traction_x_t0 = (lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic)*dxUx + lambdal_unrelaxed_elastic*dzUz
                    traction_z_t0 = mul_unrelaxed_elastic*(dxUz + dzUx)
                 else
                    veloc_horiz=v0x_left(count_left)
                    veloc_vert=v0z_left(count_left)
                    traction_x_t0=t0x_left(count_left)
                    traction_z_t0=t0z_left(count_left)
                    count_left=count_left+1
                 end if
              else
                 veloc_horiz = 0
                 veloc_vert = 0
                 traction_x_t0 = 0
                 traction_z_t0 = 0
              endif

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
              nx = - zgamma / jacobian1D
              nz = + xgamma / jacobian1D

              weight = jacobian1D * wzgll(j)

              ! Clayton-Engquist condition if elastic
              if(elastic(ispec)) then
                 vx = veloc_elastic(1,iglob) - veloc_horiz
                 vy = veloc_elastic(2,iglob)
                 vz = veloc_elastic(3,iglob) - veloc_vert

                 vn = nx*vx+nz*vz

                 tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
                 ty = rho_vs*vy
                 tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

                 displtx=0.0d0
                 displtz=0.0d0

                 if(ADD_SPRING_TO_STACEY)then

                 displx = displ_elastic(1,iglob)
                 disply = displ_elastic(2,iglob)
                 displz = displ_elastic(3,iglob)

                 spring_position=sqrt((coord(1,iglob)-x_center_spring)**2 +&
                                 (coord(2,iglob)-z_center_spring)**2)

                 displn = nx*displx+nz*displz

                 displtx = (lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic)/&
                            (2.0*spring_position)*displn*nx &
                           +mul_unrelaxed_elastic/(2.0*spring_position)*(displx-displn*nx)
                 displty = mul_unrelaxed_elastic*disply
                 displtz = (lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic)/&
                            (2.0*spring_position)*displn*nz &
                           +mul_unrelaxed_elastic/(2.0*spring_position)*(displz-displn*nz)

                 endif

                 accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + traction_x_t0+displtx)*weight
                 accel_elastic(2,iglob) = accel_elastic(2,iglob) - ty*weight
                 accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz + traction_z_t0+displtz)*weight

                 if(SAVE_FORWARD .and. SIMULATION_TYPE ==1) then
                    if(p_sv)then !P-SV waves
                       b_absorb_elastic_left(1,j,ib_left(ispecabs),it) = (tx + traction_x_t0)*weight
                       b_absorb_elastic_left(3,j,ib_left(ispecabs),it) = (tz + traction_z_t0)*weight
                    else !SH (membrane) waves
                       b_absorb_elastic_left(2,j,ib_left(ispecabs),it) = ty*weight
                    endif
                 elseif(SIMULATION_TYPE == 2) then
                    if(p_sv)then !P-SV waves
                       b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - &
                            b_absorb_elastic_left(1,j,ib_left(ispecabs),NSTEP-it+1)
                       b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) - &
                            b_absorb_elastic_left(3,j,ib_left(ispecabs),NSTEP-it+1)
                    else !SH (membrane) waves
                       b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - &
                            b_absorb_elastic_left(2,j,ib_left(ispecabs),NSTEP-it+1)
                    endif
                 endif

              endif

           enddo

        endif  !  end of left absorbing boundary

        !--- right absorbing boundary
        if(codeabs(IEDGE2,ispecabs)) then

           i = NGLLX

           do j = 1,NGLLZ

              iglob = ibool(i,j,ispec)

              ! for analytical initial plane wave for Bielak's conditions
              ! left or right edge, horizontal normal vector
              if(add_Bielak_conditions .and. initialfield) then
                 if (.not.over_critical_angle) then
                    call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                         x0_source, z0_source, A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                         c_inc, c_refl, time_offset,f0)
                    traction_x_t0 = (lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic)*dxUx + lambdal_unrelaxed_elastic*dzUz
                    traction_z_t0 = mul_unrelaxed_elastic*(dxUz + dzUx)
                 else
                    veloc_horiz=v0x_right(count_right)
                    veloc_vert=v0z_right(count_right)
                    traction_x_t0=t0x_right(count_right)
                    traction_z_t0=t0z_right(count_right)
                    count_right=count_right+1
                 end if
              else
                 veloc_horiz = 0
                 veloc_vert = 0
                 traction_x_t0 = 0
                 traction_z_t0 = 0
              endif

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

              ! Clayton-Engquist condition if elastic
              if(elastic(ispec)) then
                 vx = veloc_elastic(1,iglob) - veloc_horiz
                 vy = veloc_elastic(2,iglob)
                 vz = veloc_elastic(3,iglob) - veloc_vert

                 vn = nx*vx+nz*vz

                 tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
                 ty = rho_vs*vy
                 tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

                 displtx=0.0d0
                 displtz=0.0d0

                 if(ADD_SPRING_TO_STACEY)then

                 displx = displ_elastic(1,iglob)
                 disply = displ_elastic(2,iglob)
                 displz = displ_elastic(3,iglob)

                 spring_position=sqrt((coord(1,iglob)-x_center_spring)**2 +&
                                 (coord(2,iglob)-z_center_spring)**2)

                 displn = nx*displx+nz*displz

                 displtx = (lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic)/&
                            (2.0*spring_position)*displn*nx &
                           +mul_unrelaxed_elastic/(2.0*spring_position)*(displx-displn*nx)
                 displty = mul_unrelaxed_elastic*disply
                 displtz = (lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic)/&
                            (2.0*spring_position)*displn*nz &
                           +mul_unrelaxed_elastic/(2.0*spring_position)*(displz-displn*nz)

                 endif

                 accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx - traction_x_t0+displtx)*weight
                 accel_elastic(2,iglob) = accel_elastic(2,iglob) - ty*weight
                 accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz - traction_z_t0+displtz)*weight

                 if(SAVE_FORWARD .and. SIMULATION_TYPE ==1) then
                    if(p_sv)then !P-SV waves
                       b_absorb_elastic_right(1,j,ib_right(ispecabs),it) = (tx - traction_x_t0)*weight
                       b_absorb_elastic_right(3,j,ib_right(ispecabs),it) = (tz - traction_z_t0)*weight
                    else! SH (membrane) waves
                       b_absorb_elastic_right(2,j,ib_right(ispecabs),it) = ty*weight
                    endif
                 elseif(SIMULATION_TYPE == 2) then
                    if(p_sv)then !P-SV waves
                       b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - &
                            b_absorb_elastic_right(1,j,ib_right(ispecabs),NSTEP-it+1)
                       b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) - &
                            b_absorb_elastic_right(3,j,ib_right(ispecabs),NSTEP-it+1)
                    else! SH (membrane) waves
                       b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - &
                            b_absorb_elastic_right(2,j,ib_right(ispecabs),NSTEP-it+1)
                    endif
                 endif

              endif

           enddo

        endif  !  end of right absorbing boundary

        !--- bottom absorbing boundary
        if(codeabs(IEDGE1,ispecabs)) then

           j = 1

!! DK DK not needed           ! exclude corners to make sure there is no contradiction on the normal
           ibegin = 1
           iend = NGLLX
!! DK DK not needed           if(codeabs(IEDGE4,ispecabs)) ibegin = 2
!! DK DK not needed           if(codeabs(IEDGE2,ispecabs)) iend = NGLLX-1

           do i = ibegin,iend

              iglob = ibool(i,j,ispec)

              ! for analytical initial plane wave for Bielak's conditions
              ! top or bottom edge, vertical normal vector
              if(add_Bielak_conditions .and. initialfield) then
                 if (.not.over_critical_angle) then
                    call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                         x0_source, z0_source, A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                         c_inc, c_refl, time_offset,f0)
                    traction_x_t0 = mul_unrelaxed_elastic*(dxUz + dzUx)
                    traction_z_t0 = lambdal_unrelaxed_elastic*dxUx + (lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic)*dzUz
                 else
                    veloc_horiz=v0x_bot(count_bottom)
                    veloc_vert=v0z_bot(count_bottom)
                    traction_x_t0=t0x_bot(count_bottom)
                    traction_z_t0=t0z_bot(count_bottom)
                    count_bottom=count_bottom+1
                 end if
              else
                 veloc_horiz = 0
                 veloc_vert = 0
                 traction_x_t0 = 0
                 traction_z_t0 = 0
              endif

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

              ! Clayton-Engquist condition if elastic
              if(elastic(ispec)) then
                 vx = veloc_elastic(1,iglob) - veloc_horiz
                 vy = veloc_elastic(2,iglob)
                 vz = veloc_elastic(3,iglob) - veloc_vert

                 vn = nx*vx+nz*vz

                 tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
                 ty = rho_vs*vy
                 tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
                 if((codeabs(IEDGE4,ispecabs) .and. i == 1) .or. (codeabs(IEDGE2,ispecabs) .and. i == NGLLX)) then
                   tx = 0
                   ty = 0
                   tz = 0
                 endif

                 displtx=0.0d0
                 displtz=0.0d0

                 if(ADD_SPRING_TO_STACEY)then

                 displx = displ_elastic(1,iglob)
                 disply = displ_elastic(2,iglob)
                 displz = displ_elastic(3,iglob)

                 spring_position=sqrt((coord(1,iglob)-x_center_spring)**2 +&
                                 (coord(2,iglob)-z_center_spring)**2)

                 displn = nx*displx+nz*displz

                 displtx = (lambdal_unrelaxed_elastic+2.0*mul_unrelaxed_elastic)/&
                            (2.0*spring_position)*displn*nx &
                           +mul_unrelaxed_elastic/(2.0*spring_position)*(displx-displn*nx)
                 displty = mul_unrelaxed_elastic*disply
                 displtz = (lambdal_unrelaxed_elastic+2.0*mul_unrelaxed_elastic)/&
                            (2.0*spring_position)*displn*nz &
                           +mul_unrelaxed_elastic/(2.0*spring_position)*(displz-displn*nz)

                 if((codeabs(IEDGE4,ispecabs) .and. i == 1) .or. (codeabs(IEDGE2,ispecabs) .and. i == NGLLX)) then
                   displtx = 0
                   displty = 0
                   displtz = 0
                 endif


                 endif


                 accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + traction_x_t0+displtx)*weight
                 accel_elastic(2,iglob) = accel_elastic(2,iglob) - ty*weight
                 accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz + traction_z_t0+displtz)*weight

                 if(SAVE_FORWARD .and. SIMULATION_TYPE ==1) then
                    if(p_sv)then !P-SV waves
                       b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),it) = (tx + traction_x_t0)*weight
                       b_absorb_elastic_bottom(3,i,ib_bottom(ispecabs),it) = (tz + traction_z_t0)*weight
                    else!SH (membrane) waves
                       b_absorb_elastic_bottom(2,i,ib_bottom(ispecabs),it) = ty*weight
                    endif
                 elseif(SIMULATION_TYPE == 2) then
                    if(p_sv)then !P-SV waves
                       b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) - &
                            b_absorb_elastic_bottom(1,i,ib_bottom(ispecabs),NSTEP-it+1)
                       b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) - &
                            b_absorb_elastic_bottom(3,i,ib_bottom(ispecabs),NSTEP-it+1)
                    else!SH (membrane) waves
                       b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) - &
                            b_absorb_elastic_bottom(2,i,ib_bottom(ispecabs),NSTEP-it+1)
                    endif
                 endif

              endif

           enddo

        endif  !  end of bottom absorbing boundary

        !--- top absorbing boundary
        if(codeabs(IEDGE3,ispecabs)) then

           j = NGLLZ

!! DK DK not needed           ! exclude corners to make sure there is no contradiction on the normal
           ibegin = 1
           iend = NGLLX
!! DK DK not needed           if(codeabs(IEDGE4,ispecabs)) ibegin = 2
!! DK DK not needed           if(codeabs(IEDGE2,ispecabs)) iend = NGLLX-1

           do i = ibegin,iend

              iglob = ibool(i,j,ispec)

              ! for analytical initial plane wave for Bielak's conditions
              ! top or bottom edge, vertical normal vector
              if(add_Bielak_conditions .and. initialfield) then
                 call compute_Bielak_conditions(coord,iglob,nglob,it,deltat,dxUx,dxUz,dzUx,dzUz,veloc_horiz,veloc_vert, &
                      x0_source, z0_source, A_plane, B_plane, C_plane, anglesource(1), anglesource_refl, &
                      c_inc, c_refl, time_offset,f0)
                 traction_x_t0 = mul_unrelaxed_elastic*(dxUz + dzUx)
                 traction_z_t0 = lambdal_unrelaxed_elastic*dxUx + (lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic)*dzUz
              else
                 veloc_horiz = 0
                 veloc_vert = 0
                 traction_x_t0 = 0
                 traction_z_t0 = 0
              endif

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
              nx = - zxi / jacobian1D
              nz = + xxi / jacobian1D

              weight = jacobian1D * wxgll(i)

              ! Clayton-Engquist condition if elastic
              if(elastic(ispec)) then
                 vx = veloc_elastic(1,iglob) - veloc_horiz
                 vy = veloc_elastic(2,iglob)
                 vz = veloc_elastic(3,iglob) - veloc_vert

                 vn = nx*vx+nz*vz

                 tx = rho_vp*vn*nx+rho_vs*(vx-vn*nx)
                 ty = rho_vs*vy
                 tz = rho_vp*vn*nz+rho_vs*(vz-vn*nz)

! exclude corners to make sure there is no contradiction on the normal
! for Stacey absorbing conditions but not for incident plane waves;
! thus subtract nothing i.e. zero in that case
                 if((codeabs(IEDGE4,ispecabs) .and. i == 1) .or. (codeabs(IEDGE2,ispecabs) .and. i == NGLLX)) then
                   tx = 0
                   ty = 0
                   tz = 0
                 endif

                 displtx=0.0d0
                 displtz=0.0d0

                 if(ADD_SPRING_TO_STACEY)then

                 displx = displ_elastic(1,iglob)
                 disply = displ_elastic(2,iglob)
                 displz = displ_elastic(3,iglob)

                 spring_position=sqrt((coord(1,iglob)-x_center_spring)**2 +&
                                 (coord(2,iglob)-z_center_spring)**2)

                 displn = nx*displx+nz*displz

                 displtx = (lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic)/&
                            (2.0*spring_position)*displn*nx &
                           +mul_unrelaxed_elastic/(2.0*spring_position)*(displx-displn*nx)
                 displty = mul_unrelaxed_elastic*disply
                 displtz = (lambdal_unrelaxed_elastic+2*mul_unrelaxed_elastic)/&
                            (2.0*spring_position)*displn*nz &
                           +mul_unrelaxed_elastic/(2.0*spring_position)*(displz-displn*nz)

                 if((codeabs(IEDGE4,ispecabs) .and. i == 1) .or. (codeabs(IEDGE2,ispecabs) .and. i == NGLLX)) then
                   displtx = 0
                   displty = 0
                   displtz = 0
                 endif


                 endif


                 accel_elastic(1,iglob) = accel_elastic(1,iglob) - (tx + traction_x_t0+displtx)*weight
                 accel_elastic(2,iglob) = accel_elastic(2,iglob) - ty*weight
                 accel_elastic(3,iglob) = accel_elastic(3,iglob) - (tz + traction_z_t0+displtz)*weight

                 if(SAVE_FORWARD .and. SIMULATION_TYPE ==1) then
                    if(p_sv)then !P-SV waves
                       b_absorb_elastic_top(1,i,ib_top(ispecabs),it) = (tx- traction_x_t0)*weight
                       b_absorb_elastic_top(3,i,ib_top(ispecabs),it) = (tz- traction_z_t0)*weight
                    else!SH (membrane) waves
                       b_absorb_elastic_top(2,i,ib_top(ispecabs),it) = ty*weight
                    endif
                 elseif(SIMULATION_TYPE == 2) then
                    if(p_sv)then !P-SV waves
                       b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) &
                        - b_absorb_elastic_top(1,i,ib_top(ispecabs),NSTEP-it+1)
                       b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) &
                        - b_absorb_elastic_top(3,i,ib_top(ispecabs),NSTEP-it+1)
                    else!SH (membrane) waves
                       b_accel_elastic(2,iglob) = b_accel_elastic(2,iglob) &
                        - b_absorb_elastic_top(2,i,ib_top(ispecabs),NSTEP-it+1)
                    endif
                 endif

              endif

           enddo

        endif  !  end of top absorbing boundary

     enddo

  endif  ! end of absorbing boundaries

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !
  !--- absorbing boundaries
  !
  if( PML_BOUNDARY_CONDITIONS .and. anyabs) then

! we have to put Dirichlet on the boundary of the PML
     do ispecabs = 1,nelemabs

       ispec = numabs(ispecabs)

       if (is_PML(ispec)) then
      ispec_PML=spec_to_PML(ispec)
!--- left absorbing boundary
      if(codeabs(IEDGE4,ispecabs)) then
        i = 1
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          displ_elastic(:,iglob) = 0._CUSTOM_REAL
          veloc_elastic(:,iglob) = 0._CUSTOM_REAL
          accel_elastic(:,iglob) = 0._CUSTOM_REAL
       enddo
      endif
!--- right absorbing boundary
      if(codeabs(IEDGE2,ispecabs)) then
        i = NGLLX
        do j = 1,NGLLZ
          iglob = ibool(i,j,ispec)
          displ_elastic(:,iglob) = 0._CUSTOM_REAL
          veloc_elastic(:,iglob) = 0._CUSTOM_REAL
          accel_elastic(:,iglob) = 0._CUSTOM_REAL
        enddo

      endif
!--- bottom absorbing boundary
      if(codeabs(IEDGE1,ispecabs)) then
        j = 1
! exclude corners to make sure there is no contradiction on the normal
        ibegin = 1
        iend = NGLLX
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          displ_elastic(:,iglob) = 0._CUSTOM_REAL
          veloc_elastic(:,iglob) = 0._CUSTOM_REAL
          accel_elastic(:,iglob) = 0._CUSTOM_REAL
        enddo
      endif
!--- top absorbing boundary
      if(codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ
! exclude corners to make sure there is no contradiction on the normal
        ibegin = 1
        iend = NGLLX
        do i = ibegin,iend
          iglob = ibool(i,j,ispec)
          displ_elastic(:,iglob) = 0._CUSTOM_REAL
          veloc_elastic(:,iglob) = 0._CUSTOM_REAL
          accel_elastic(:,iglob) = 0._CUSTOM_REAL
        enddo
      endif  !  end of top absorbing boundary
      endif ! end of is_PML
    enddo ! end specabs loop
  endif  ! end of absorbing boundaries PML_BOUNDARY_CONDITIONS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


  ! --- add the source if it is a moment tensor
  if(.not. initialfield) then

     do i_source=1,NSOURCES
        ! if this processor core carries the source and the source element is elastic
        if (is_proc_source(i_source) == 1 .and. elastic(ispec_selected_source(i_source))) then

           ! moment tensor
           if(source_type(i_source) == 2) then

              if(.not.p_sv)  call exit_MPI('cannot have moment tensor source in SH (membrane) waves calculation')

              if(SIMULATION_TYPE == 1) then  ! forward wavefield
                 ! add source array
                 do j=1,NGLLZ
                    do i=1,NGLLX
                       iglob = ibool(i,j,ispec_selected_source(i_source))
                       accel_elastic(1,iglob) = accel_elastic(1,iglob) + &
                            sourcearray(i_source,1,i,j)*source_time_function(i_source,it,i_stage)
                       accel_elastic(3,iglob) = accel_elastic(3,iglob) + &
                            sourcearray(i_source,2,i,j)*source_time_function(i_source,it,i_stage)
                    enddo
                 enddo
              else                   ! backward wavefield
                 do j=1,NGLLZ
                    do i=1,NGLLX
                       iglob = ibool(i,j,ispec_selected_source(i_source))
                       b_accel_elastic(1,iglob) = b_accel_elastic(1,iglob) + &
                            sourcearray(i_source,1,i,j)*source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)
                       b_accel_elastic(3,iglob) = b_accel_elastic(3,iglob) + &
                            sourcearray(i_source,2,i,j)*source_time_function(i_source,NSTEP-it+1,stage_time_scheme-i_stage+1)
                    enddo
                 enddo
              endif  !endif SIMULATION_TYPE == 1

           endif !if(source_type(i_source) == 2)

        endif ! if this processor core carries the source and the source element is elastic
     enddo ! do i_source=1,NSOURCES

     if(SIMULATION_TYPE == 2) then   ! adjoint wavefield

        irec_local = 0
        do irec = 1,nrec
           !   add the source (only if this proc carries the source)
           if(myrank == which_proc_receiver(irec)) then

              irec_local = irec_local + 1
              if(elastic(ispec_selected_rec(irec))) then
                 ! add source array
                 do j=1,NGLLZ
                    do i=1,NGLLX
                       iglob = ibool(i,j,ispec_selected_rec(irec))
                       if(p_sv)then !P-SH waves
                          accel_elastic(1,iglob) = accel_elastic(1,iglob) + adj_sourcearrays(irec_local,NSTEP-it+1,1,i,j)
                          accel_elastic(3,iglob) = accel_elastic(3,iglob) + adj_sourcearrays(irec_local,NSTEP-it+1,3,i,j)
                       else !SH (membrane) waves
                          accel_elastic(2,iglob) = accel_elastic(2,iglob) + adj_sourcearrays(irec_local,NSTEP-it+1,2,i,j)
                       endif
                    enddo
                 enddo
              endif ! if element is elastic

           endif ! if this processor core carries the adjoint source and the source element is elastic
        enddo ! irec = 1,nrec

     endif ! if SIMULATION_TYPE == 2 adjoint wavefield

  endif ! if not using an initial field

  ! implement attenuation
  if(ATTENUATION_VISCOELASTIC_SOLID) then

     ! compute Grad(displ_elastic) at time step n+1 for attenuation
     call compute_gradient_attenuation(temp_displ_elastic,dux_dxl_np1,duz_dxl_np1, &
          dux_dzl_np1,duz_dzl_np1,xix,xiz,gammax,gammaz,ibool,elastic,hprime_xx,hprime_zz,nspec,nglob)

     ! update memory variables with fourth-order Runge-Kutta time scheme for attenuation
     ! loop over spectral elements
     do ispec = 1,nspec

        do j=1,NGLLZ
           do i=1,NGLLX

              theta_n   = dux_dxl_n(i,j,ispec) + duz_dzl_n(i,j,ispec)
              theta_np1 = dux_dxl_np1(i,j,ispec) + duz_dzl_np1(i,j,ispec)

              ! loop on all the standard linear solids
              do i_sls = 1,N_SLS

                 ! evolution e1
                 if(stage_time_scheme == 1) then
                    Un = e1(i,j,ispec,i_sls)
                    tauinv = - inv_tau_sigma_nu1(i,j,ispec,i_sls)
                    tauinvsquare = tauinv * tauinv
                    tauinvcube = tauinvsquare * tauinv
                    tauinvUn = tauinv * Un
                    Sn   = theta_n * phi_nu1(i,j,ispec,i_sls)
                    Snp1 = theta_np1 * phi_nu1(i,j,ispec,i_sls)
                    Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
                         twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
                         fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
                         deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
                    e1(i,j,ispec,i_sls) = Unp1
                 endif

                 if(stage_time_scheme == 6) then

                    tauinv = inv_tau_sigma_nu1(i,j,ispec,i_sls)
                    Un = e1(i,j,ispec,i_sls)
                    temp_time_scheme = phi_nu1(i,j,ispec,i_sls)
                    e1_LDDRK(i,j,ispec,i_sls) = alpha_LDDRK(i_stage) * e1_LDDRK(i,j,ispec,i_sls) + &
                                              deltat * (theta_n * temp_time_scheme - Un * tauinv)
                    e1(i,j,ispec,i_sls) = Un + beta_LDDRK(i_stage) * e1_LDDRK(i,j,ispec,i_sls)
                 endif

                 if(stage_time_scheme == 4) then

                    tauinv = inv_tau_sigma_nu1(i,j,ispec,i_sls)
                    Un = e1(i,j,ispec,i_sls)
                    temp_time_scheme = phi_nu1(i,j,ispec,i_sls)
                    e1_force_RK(i,j,ispec,i_sls,i_stage) = deltat * (theta_n * temp_time_scheme - Un * tauinv)

                    if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then
                       if(i_stage == 1)weight_rk = 0.5d0
                       if(i_stage == 2)weight_rk = 0.5d0
                       if(i_stage == 3)weight_rk = 1.0d0

                       if(i_stage==1)then
                       e1_initial_rk(i,j,ispec,i_sls) = e1(i,j,ispec,i_sls)
                       endif

                       e1(i,j,ispec,i_sls) = e1_initial_rk(i,j,ispec,i_sls) &
                        + weight_rk * e1_force_RK(i,j,ispec,i_sls,i_stage)

                    elseif(i_stage==4)then

                       e1(i,j,ispec,i_sls) = e1_initial_rk(i,j,ispec,i_sls) + 1.0d0 / 6.0d0 * &
                                             (e1_force_RK(i,j,ispec,i_sls,1) + 2.0d0 * e1_force_RK(i,j,ispec,i_sls,2) + &
                                             2.0d0 * e1_force_RK(i,j,ispec,i_sls,3) + e1_force_RK(i,j,ispec,i_sls,4))
                    endif
                 endif


                 ! evolution e11
                 if(stage_time_scheme == 1) then
                    Un = e11(i,j,ispec,i_sls)
                    tauinv = - inv_tau_sigma_nu2(i,j,ispec,i_sls)
                    tauinvsquare = tauinv * tauinv
                    tauinvcube = tauinvsquare * tauinv
                    tauinvUn = tauinv * Un
                    Sn   = (dux_dxl_n(i,j,ispec) - theta_n/TWO) * phi_nu2(i,j,ispec,i_sls)
                    Snp1 = (dux_dxl_np1(i,j,ispec) - theta_np1/TWO) * phi_nu2(i,j,ispec,i_sls)
                    Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
                         twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
                         fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
                         deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
                    e11(i,j,ispec,i_sls) = Unp1
                endif

                 if(stage_time_scheme == 6) then
                    temp_time_scheme = dux_dxl_n(i,j,ispec)-theta_n/TWO
                    temper_time_scheme = phi_nu2(i,j,ispec,i_sls)
                    e11_LDDRK(i,j,ispec,i_sls) = alpha_LDDRK(i_stage) * e11_LDDRK(i,j,ispec,i_sls) &
                                                 + deltat * (temp_time_scheme * temper_time_scheme)- &
                                                 deltat * (e11(i,j,ispec,i_sls) * inv_tau_sigma_nu2(i,j,ispec,i_sls))
                    e11(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls)+beta_LDDRK(i_stage)*e11_LDDRK(i,j,ispec,i_sls)
                 endif

                 if(stage_time_scheme == 4) then

                    temp_time_scheme = dux_dxl_n(i,j,ispec)-theta_n/TWO
                    temper_time_scheme = phi_nu2(i,j,ispec,i_sls)
                    e11_force_RK(i,j,ispec,i_sls,i_stage) = deltat * (temp_time_scheme * temper_time_scheme- &
                                               e11(i,j,ispec,i_sls) * inv_tau_sigma_nu2(i,j,ispec,i_sls))

                    if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then
                       if(i_stage == 1)weight_rk = 0.5d0
                       if(i_stage == 2)weight_rk = 0.5d0
                       if(i_stage == 3)weight_rk = 1.0d0
                       if(i_stage==1)then
                          e11_initial_rk(i,j,ispec,i_sls) = e11(i,j,ispec,i_sls)
                       endif
                       e11(i,j,ispec,i_sls) = e11_initial_rk(i,j,ispec,i_sls) &
                        + weight_rk * e11_force_RK(i,j,ispec,i_sls,i_stage)
                    elseif(i_stage==4)then
                      e11(i,j,ispec,i_sls) = e11_initial_rk(i,j,ispec,i_sls) + 1.0d0 / 6.0d0 * &
                                             (e11_force_RK(i,j,ispec,i_sls,1) + 2.0d0 * e11_force_RK(i,j,ispec,i_sls,2) + &
                                              2.0d0 * e11_force_RK(i,j,ispec,i_sls,3) + e11_force_RK(i,j,ispec,i_sls,4))
                    endif
                 endif

                 ! evolution e13
                 if(stage_time_scheme == 1) then
                    Un = e13(i,j,ispec,i_sls)
                    tauinv = - inv_tau_sigma_nu2(i,j,ispec,i_sls)
                    tauinvsquare = tauinv * tauinv
                    tauinvcube = tauinvsquare * tauinv
                    tauinvUn = tauinv * Un
                    Sn   = (dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)) * phi_nu2(i,j,ispec,i_sls)
                    Snp1 = (dux_dzl_np1(i,j,ispec) + duz_dxl_np1(i,j,ispec)) * phi_nu2(i,j,ispec,i_sls)
                    Unp1 = Un + (deltatfourth*tauinvcube*(Sn + tauinvUn) + &
                         twelvedeltat*(Sn + Snp1 + 2*tauinvUn) + &
                         fourdeltatsquare*tauinv*(2*Sn + Snp1 + 3*tauinvUn) + &
                         deltatcube*tauinvsquare*(3*Sn + Snp1 + 4*tauinvUn))* ONE_OVER_24
                    e13(i,j,ispec,i_sls) = Unp1
                 endif

                 if(stage_time_scheme == 6) then
                    temp_time_scheme=dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)
                    temper_time_scheme=phi_nu2(i,j,ispec,i_sls)
                    e13_LDDRK(i,j,ispec,i_sls) = alpha_LDDRK(i_stage) * e13_LDDRK(i,j,ispec,i_sls)+&
                                             deltat * (temp_time_scheme*temper_time_scheme)- &
                                            deltat * (e13(i,j,ispec,i_sls) * inv_tau_sigma_nu2(i,j,ispec,i_sls))
                    e13(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls)+beta_LDDRK(i_stage) * e13_LDDRK(i,j,ispec,i_sls)
                 endif

                 if(stage_time_scheme == 4) then
                    temp_time_scheme=dux_dzl_n(i,j,ispec) + duz_dxl_n(i,j,ispec)
                    temper_time_scheme=phi_nu2(i,j,ispec,i_sls)
                    e13_force_RK(i,j,ispec,i_sls,i_stage) = deltat * (temp_time_scheme*temper_time_scheme- &
                                            e13(i,j,ispec,i_sls) * inv_tau_sigma_nu2(i,j,ispec,i_sls))
                    if(i_stage==1 .or. i_stage==2 .or. i_stage==3)then
                       if(i_stage == 1)weight_rk = 0.5d0
                       if(i_stage == 2)weight_rk = 0.5d0
                       if(i_stage == 3)weight_rk = 1.0d0
                       if(i_stage==1)then
                          e13_initial_rk(i,j,ispec,i_sls) = e13(i,j,ispec,i_sls)
                       endif
                          e13(i,j,ispec,i_sls) = e13_initial_rk(i,j,ispec,i_sls) &
                            + weight_rk * e13_force_RK(i,j,ispec,i_sls,i_stage)
                    elseif(i_stage==4)then
                       e13(i,j,ispec,i_sls) = e13_initial_rk(i,j,ispec,i_sls) + 1.0d0 / 6.0d0 * &
                                              (e13_force_RK(i,j,ispec,i_sls,1) + 2.0d0 * e13_force_RK(i,j,ispec,i_sls,2) + &
                                              2.0d0 * e13_force_RK(i,j,ispec,i_sls,3) + e13_force_RK(i,j,ispec,i_sls,4))
                    endif
                 endif

              enddo

           enddo
        enddo
     enddo

  endif ! end of test on attenuation

end subroutine compute_forces_viscoelastic


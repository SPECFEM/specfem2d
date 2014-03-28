
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, Inria and University of Pau, France,
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

  subroutine compute_pressure_whole_medium(potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic,displ_elastic,&
                  displs_poroelastic,displw_poroelastic,acoustic,gravitoacoustic,elastic,poroelastic,vector_field_display, &
                  xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec, &
                  AXISYM,coord,jacobian,is_on_the_axis,hprimeBar_xx, &
                  nglob,nglob_acoustic,nglob_gravitoacoustic,nglob_elastic,nglob_poroelastic,assign_external_model, &
                  numat,kmato,density,porosity,tortuosity,poroelastcoef,vpext,vsext,rhoext, &
                  c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropic,anisotropy,e1,e11, &
                  ATTENUATION_VISCOELASTIC_SOLID,Mu_nu1,Mu_nu2,N_SLS)

! compute pressure in acoustic elements and in elastic elements

  implicit none

  include "constants.h"

  integer :: nspec,nglob,numat


  integer, dimension(nspec) :: kmato
  integer, dimension(NGLLX,NGLLX,nspec) :: ibool

  double precision, dimension(2,numat) :: density
  double precision, dimension(numat) :: porosity,tortuosity
  double precision, dimension(4,3,numat) :: poroelastcoef
  double precision, dimension(9,numat) :: anisotropy
  double precision, dimension(NGLLX,NGLLX,nspec) :: vpext,vsext,rhoext
  double precision, dimension(NGLLX,NGLLZ,nspec) ::  c11ext,c15ext,c13ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

  logical, dimension(nspec) :: acoustic,gravitoacoustic,elastic,poroelastic,anisotropic
  integer :: nglob_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: potential_dot_dot_acoustic
  integer :: nglob_gravitoacoustic
  real(kind=CUSTOM_REAL), dimension(nglob_gravitoacoustic) :: potential_dot_dot_gravitoacoustic
  integer :: nglob_elastic
  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: displ_elastic
  integer :: nglob_poroelastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic) :: displs_poroelastic,displw_poroelastic

  double precision, dimension(3,nglob) :: vector_field_display

! array with derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz


  logical :: AXISYM
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLJ) :: hprimeBar_xx
  logical, dimension(nspec) :: is_on_the_axis
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: jacobian
  double precision, dimension(NDIM,nglob), intent(in) :: coord

  logical :: assign_external_model,ATTENUATION_VISCOELASTIC_SOLID

  integer :: N_SLS
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec,N_SLS) :: e1,e11
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: Mu_nu1,Mu_nu2

! local variables
  integer :: i,j,ispec,iglob

! pressure in this element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: pressure_element

! loop over spectral elements
  do ispec = 1,nspec

! compute pressure in this element
    call compute_pressure_one_element(pressure_element,potential_dot_dot_acoustic,potential_dot_dot_gravitoacoustic,displ_elastic,&
         displs_poroelastic,displw_poroelastic,acoustic,gravitoacoustic,elastic,poroelastic,&
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec, &
         AXISYM,nglob,coord,jacobian,is_on_the_axis,hprimeBar_xx, &
         nglob_acoustic,nglob_gravitoacoustic,nglob_elastic,nglob_poroelastic,assign_external_model, &
         numat,kmato,density,porosity,tortuosity,poroelastcoef,vpext,vsext,rhoext, &
         c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropic,anisotropy,ispec,e1,e11, &
         ATTENUATION_VISCOELASTIC_SOLID,Mu_nu1,Mu_nu2,N_SLS)

! use vector_field_display as temporary storage, store pressure in its second component
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        vector_field_display(3,iglob) = pressure_element(i,j)
      enddo
    enddo

  enddo

  end subroutine compute_pressure_whole_medium

!
!=====================================================================
!

  subroutine compute_pressure_one_element(pressure_element,potential_dot_dot_acoustic, &
         potential_dot_dot_gravitoacoustic,displ_elastic,&
         displs_poroelastic,displw_poroelastic,acoustic,gravitoacoustic,elastic,poroelastic,&
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec, &
         AXISYM,nglob,coord,jacobian,is_on_the_axis,hprimeBar_xx, &
         nglob_acoustic,nglob_gravitoacoustic,nglob_elastic,nglob_poroelastic,assign_external_model, &
         numat,kmato,density,porosity,tortuosity,poroelastcoef,vpext,vsext,rhoext, &
         c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext,anisotropic,anisotropy,ispec,e1,e11, &
         ATTENUATION_VISCOELASTIC_SOLID,Mu_nu1,Mu_nu2,N_SLS)

! compute pressure in acoustic elements and in elastic elements

  implicit none

  include "constants.h"

  integer nspec,numat,ispec

  integer, dimension(nspec) :: kmato
  integer, dimension(NGLLX,NGLLX,nspec) :: ibool

  double precision, dimension(2,numat) :: density
  double precision, dimension(numat) :: porosity,tortuosity
  double precision, dimension(4,3,numat) :: poroelastcoef
  double precision, dimension(9,numat) :: anisotropy
  double precision, dimension(NGLLX,NGLLX,nspec) :: vpext,vsext,rhoext
  double precision, dimension(NGLLX,NGLLZ,nspec) ::  c11ext,c15ext,c13ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

! pressure in this element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: pressure_element

  logical, dimension(nspec) :: acoustic,gravitoacoustic,elastic,poroelastic,anisotropic
  integer :: nglob_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: potential_dot_dot_acoustic
  integer :: nglob_gravitoacoustic
  real(kind=CUSTOM_REAL), dimension(nglob_gravitoacoustic) :: potential_dot_dot_gravitoacoustic
  integer :: nglob_elastic
  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: displ_elastic
  integer :: nglob_poroelastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic) :: displs_poroelastic,displw_poroelastic

! array with derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz


  logical :: AXISYM
  integer :: nglob
  real(kind=CUSTOM_REAL), dimension(NGLJ,NGLJ) :: hprimeBar_xx
  logical, dimension(nspec) :: is_on_the_axis
  double precision, dimension(NDIM,nglob), intent(in) :: coord
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: jacobian

  logical :: assign_external_model,ATTENUATION_VISCOELASTIC_SOLID

  integer :: N_SLS
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec,N_SLS) :: e1,e11
  real(kind=CUSTOM_REAL) :: e1_sum,e11_sum
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: Mu_nu1,Mu_nu2
  integer :: i_sls

! local variables
  integer :: i,j,k,iglob

! jacobian
  real(kind=CUSTOM_REAL) :: xixl,xizl,gammaxl,gammazl

! spatial derivatives
  real(kind=CUSTOM_REAL) :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  real(kind=CUSTOM_REAL) :: dux_dxl,duz_dxl,dux_dzl,duz_dzl
  real(kind=CUSTOM_REAL) :: sigma_xx,sigma_yy,sigma_zz !! ,sigmap
  real(kind=CUSTOM_REAL) :: sigma_thetatheta,xxi
  real(kind=CUSTOM_REAL) :: dwx_dxi,dwx_dgamma,dwz_dxi,dwz_dgamma
  real(kind=CUSTOM_REAL) :: dwx_dxl,dwz_dzl

! material properties of the elastic medium
  real(kind=CUSTOM_REAL) :: mul_unrelaxed_elastic,lambdal_unrelaxed_elastic,lambdaplus2mu_unrelaxed_elastic,denst
  real(kind=CUSTOM_REAL) :: mul_relaxed_viscoelastic,lambdal_relaxed_viscoelastic,lambdalplus2mul_relaxed_viscoel,cpl,csl

  real(kind=CUSTOM_REAL) :: mul_s,kappal_s,rhol_s
  real(kind=CUSTOM_REAL) :: kappal_f,rhol_f
  real(kind=CUSTOM_REAL) :: mul_fr,kappal_fr,phil,tortl
  real(kind=CUSTOM_REAL) :: D_biot,H_biot,C_biot,M_biot,rhol_bar
  real(kind=CUSTOM_REAL) :: mul_G,lambdal_G,lambdalplus2mul_G

! for anisotropy
  double precision ::  c11,c15,c13,c33,c35,c55,c12,c23,c25

! if elastic element
!
! from L. S. Bennethum, Compressibility Moduli for Porous Materials Incorporating Volume Fraction,
! J. Engrg. Mech., vol. 132(11), p. 1205-1214 (2006), below equation (5):
! for a 3D isotropic solid, pressure is defined in terms of the trace of the stress tensor as
! p = -1/3 (t11 + t22 + t33) where t is the Cauchy stress tensor.

! to compute pressure in 3D in an elastic solid, one uses pressure = - trace(sigma) / 3
! sigma_ij = lambda delta_ij trace(epsilon) + 2 mu epsilon_ij
!          = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_ij
! sigma_xx = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_xx
! sigma_yy = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_yy
! sigma_zz = lambda (epsilon_xx + epsilon_yy + epsilon_zz) + 2 mu epsilon_zz
! pressure = - trace(sigma) / 3 = - (lambda + 2/3 mu) trace(epsilon) = - kappa * trace(epsilon)
!
! to compute pressure in 2D in an elastic solid in the plane strain convention i.e. in the P-SV case,
! one still uses pressure = - trace(sigma) / 3 but taking into account the fact
! that the off-plane strain epsilon_zz is zero by definition of the plane strain convention
! but thus the off-plane stress sigma_zz is not equal to zero,
! one has instead:  sigma_zz = lambda * (epsilon_xx + epsilon_yy), thus
! sigma_ij = lambda delta_ij trace(epsilon) + 2 mu epsilon_ij
!          = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_ij
! sigma_xx = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_xx
! sigma_yy = lambda (epsilon_xx + epsilon_yy) + 2 mu epsilon_yy
! sigma_zz = lambda * (epsilon_xx + epsilon_yy)
! pressure = - trace(sigma) / 3 = - (lambda + 2*mu/3) (epsilon_xx + epsilon_yy)

  if(elastic(ispec)) then

    ! get relaxed elastic parameters of current spectral element
    lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
    mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))
    lambdaplus2mu_unrelaxed_elastic = poroelastcoef(3,1,kmato(ispec))

    do j = 1,NGLLZ
      do i = 1,NGLLX

        !--- if external medium, get elastic parameters of current grid point
        if(assign_external_model) then
          cpl = vpext(i,j,ispec)
          csl = vsext(i,j,ispec)
          denst = rhoext(i,j,ispec)
          mul_unrelaxed_elastic = denst*csl*csl
          lambdal_unrelaxed_elastic = denst*cpl*cpl - TWO*mul_unrelaxed_elastic
        endif

        ! derivative along x and along z
        dux_dxi = ZERO
        duz_dxi = ZERO

        dux_dgamma = ZERO
        duz_dgamma = ZERO

        ! first double loop over GLL points to compute and store gradients
        ! we can merge the two loops because NGLLX == NGLLZ

        if(AXISYM) then
          if (is_on_the_axis(ispec)) then
            do k = 1,NGLJ
              dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
              duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprimeBar_xx(i,k)
              dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
              duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
            enddo
          else
            do k = 1,NGLJ
              dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
              duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
              dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
              duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
            enddo
          endif
        else
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(i,k)
            duz_dxi = duz_dxi + displ_elastic(3,ibool(k,j,ispec))*hprime_xx(i,k)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(j,k)
            duz_dgamma = duz_dgamma + displ_elastic(3,ibool(i,k,ispec))*hprime_zz(j,k)
          enddo
        endif

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! derivatives of displacement
        dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

        if (AXISYM .and. (abs(coord(1,ibool(i,j,ispec))) < TINYVAL)) then ! du_z/dr=0 on the axis
          duz_dxl = 0.d0
        endif

! compute diagonal components of the stress tensor (include attenuation or anisotropy if needed)

        if(ATTENUATION_VISCOELASTIC_SOLID) then

! attenuation is implemented following the memory variable formulation of
! J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993). More details can be found in
! J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation in a linear
! viscoelastic medium, Geophysical Journal International, vol. 95, p. 597-611 (1988).

! When implementing viscoelasticity according to Carcione 1993 paper, the attenuation is
! non-causal rather than causal. We fixed the problem by using equations in Carcione's
! 2004 paper and his 2007 book.

! J. M. Carcione, H. B. Helle, The physics and simulation of wave propagation at the ocean
!  bottom, Geophysics, vol. 69(3), p. 825-839, 2004
! J. M. Carcione, Wave fields in real media: wave propagation in anisotropic, anelastic
!  and porous media, Elsevier, p. 124-125, 2007

          ! compute unrelaxed elastic coefficients from formulas in Carcione 2007 page 125
          lambdal_relaxed_viscoelastic = (lambdal_unrelaxed_elastic + mul_unrelaxed_elastic) / Mu_nu1(i,j,ispec) &
                            - mul_unrelaxed_elastic / Mu_nu2(i,j,ispec)
          mul_relaxed_viscoelastic = mul_unrelaxed_elastic / Mu_nu2(i,j,ispec)
          lambdalplus2mul_relaxed_viscoel = lambdal_relaxed_viscoelastic + TWO*mul_relaxed_viscoelastic

          ! compute the stress using the unrelaxed Lame parameters (Carcione 2007 page 125)
          sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
          sigma_yy = lambdal_unrelaxed_elastic*(dux_dxl + duz_dzl)
          sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl

          ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
          ! beware: there is a bug in Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below
          e1_sum = 0._CUSTOM_REAL
          e11_sum = 0._CUSTOM_REAL

          do i_sls = 1,N_SLS
            e1_sum = e1_sum + e1(i,j,ispec,i_sls)
            e11_sum = e11_sum + e11(i,j,ispec,i_sls)
          enddo

          sigma_xx = sigma_xx + (lambdal_relaxed_viscoelastic + mul_relaxed_viscoelastic) * e1_sum &
                      + TWO * mul_unrelaxed_elastic * e11_sum
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
          sigma_yy = sigma_yy + (lambdal_relaxed_viscoelastic + mul_relaxed_viscoelastic) * e1_sum
          sigma_zz = sigma_zz + (lambdal_relaxed_viscoelastic + mul_relaxed_viscoelastic) * e1_sum &
                      - TWO * mul_relaxed_viscoelastic * e11_sum

        else

          ! no attenuation

          if(AXISYM) then
            if (is_on_the_axis(ispec)) then
              if (abs(coord(1,ibool(i,j,ispec))) < TINYVAL) then ! First GLJ point
                sigma_xx = 0._CUSTOM_REAL
                sigma_zz = 0._CUSTOM_REAL
                sigma_thetatheta = 0._CUSTOM_REAL
                xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
                do k = 1,NGLJ
                  sigma_xx = sigma_xx + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  sigma_zz = sigma_zz + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                  sigma_thetatheta = sigma_thetatheta + displ_elastic(1,ibool(k,j,ispec))*hprimeBar_xx(i,k)
                enddo
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                           + lambdal_unrelaxed_elastic*sigma_xx/xxi
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                           + lambdal_unrelaxed_elastic*sigma_zz/xxi
                sigma_thetatheta = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                   + lambdaplus2mu_unrelaxed_elastic*sigma_thetatheta/xxi
              else ! Not first GLJ point
                sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                           + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                           + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
                sigma_thetatheta = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                        + lambdaplus2mu_unrelaxed_elastic &
                                        * displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              endif
            else ! Not on the axis
              sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl &
                         + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                         + lambdal_unrelaxed_elastic*displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
              sigma_thetatheta = lambdal_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl &
                                      + lambdaplus2mu_unrelaxed_elastic &
                                      * displ_elastic(1,ibool(i,j,ispec))/coord(1,ibool(i,j,ispec))
            endif
            sigma_yy = sigma_thetatheta
          else ! Not axisym
            sigma_xx = lambdaplus2mu_unrelaxed_elastic*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
            ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
            sigma_yy = lambdal_unrelaxed_elastic*(dux_dxl + duz_dzl)
            sigma_zz = lambdaplus2mu_unrelaxed_elastic*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl
          endif

        endif

        ! full anisotropy
        if(anisotropic(ispec)) then
          if(assign_external_model) then
            c11 = c11ext(i,j,ispec)
            c15 = c15ext(i,j,ispec)
            c13 = c13ext(i,j,ispec)
            c33 = c33ext(i,j,ispec)
            c35 = c35ext(i,j,ispec)
            c55 = c55ext(i,j,ispec)
            c12 = c12ext(i,j,ispec)
            c23 = c23ext(i,j,ispec)
            c25 = c25ext(i,j,ispec)
          else
            c11 = anisotropy(1,kmato(ispec))
            c13 = anisotropy(2,kmato(ispec))
            c15 = anisotropy(3,kmato(ispec))
            c33 = anisotropy(4,kmato(ispec))
            c35 = anisotropy(5,kmato(ispec))
            c55 = anisotropy(6,kmato(ispec))
            c12 = anisotropy(7,kmato(ispec))
            c23 = anisotropy(8,kmato(ispec))
            c25 = anisotropy(9,kmato(ispec))
          endif

          duz_dxl = duz_dxi*xixl + duz_dgamma*gammaxl
          dux_dzl = dux_dxi*xizl + dux_dgamma*gammazl

          ! implement anisotropy in 2D
          sigma_xx = c11*dux_dxl + c13*duz_dzl + c15*(duz_dxl + dux_dzl)
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
          if(c12 < 1.e-7 .or. c23 < 1.e-7) stop 'cannot compute pressure for an anisotropic material if c12 or c23 are zero'
          sigma_yy = c12*dux_dxl + c23*duz_dzl + c25*(duz_dxl + dux_dzl)
          sigma_zz = c13*dux_dxl + c33*duz_dzl + c35*(duz_dxl + dux_dzl)

        endif

        ! store pressure
        ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
        pressure_element(i,j) = - (sigma_xx + sigma_yy + sigma_zz) / 3.d0

      enddo
    enddo

  else if(poroelastic(ispec)) then

    lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
    mul_unrelaxed_elastic = poroelastcoef(2,1,kmato(ispec))

    ! get poroelastic parameters of current spectral element
    phil = porosity(kmato(ispec))
    tortl = tortuosity(kmato(ispec))
    ! solid properties
    mul_s = poroelastcoef(2,1,kmato(ispec))
    kappal_s = poroelastcoef(3,1,kmato(ispec)) - FOUR_THIRDS*mul_s
    rhol_s = density(1,kmato(ispec))
    ! fluid properties
    kappal_f = poroelastcoef(1,2,kmato(ispec))
    rhol_f = density(2,kmato(ispec))
    ! frame properties
    mul_fr = poroelastcoef(2,3,kmato(ispec))
    kappal_fr = poroelastcoef(3,3,kmato(ispec)) - FOUR_THIRDS*mul_fr
    rhol_bar =  (1.d0 - phil)*rhol_s + phil*rhol_f
    ! Biot coefficients for the input phi
    D_biot = kappal_s*(1.d0 + phil*(kappal_s/kappal_f - 1.d0))
    H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) &
            + kappal_fr + FOUR_THIRDS*mul_fr
    C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
    M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
    ! where T = G:grad u_s + C div w I
    ! and T_f = C div u_s I + M div w I
    ! we are expressing lambdaplus2mu, lambda, and mu for G, C, and M
    mul_G = mul_fr
    lambdal_G = H_biot - TWO*mul_fr
    lambdalplus2mul_G = lambdal_G + TWO*mul_G

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

        ! derivatives of displacement
        dux_dxl = dux_dxi*xixl + dux_dgamma*gammaxl
        duz_dzl = duz_dxi*xizl + duz_dgamma*gammazl

        dwx_dxl = dwx_dxi*xixl + dwx_dgamma*gammaxl
        dwz_dzl = dwz_dxi*xizl + dwz_dgamma*gammazl

! compute diagonal components of the stress tensor (include attenuation if needed)

        if(ATTENUATION_VISCOELASTIC_SOLID) then
!-------------------- ATTENTION TO BE DEFINED ------------------------------!

! attenuation is implemented following the memory variable formulation of
! J. M. Carcione, Seismic modeling in viscoelastic media, Geophysics,
! vol. 58(1), p. 110-120 (1993). More details can be found in
! J. M. Carcione, D. Kosloff and R. Kosloff, Wave propagation simulation in a linear
! viscoelastic medium, Geophysical Journal International, vol. 95, p. 597-611 (1988).

! When implement viscoelasticity according to Carcione 1993 paper, the attenuation is
! non-causal rather than causal. We fixed the problem by using equations in Carcione's
! 2004 paper and his 2007 book.

! J. M. Carcione, H B. Helle, The physics and simulation of wave propagation at the ocean
!  bottom, Geophysics, vol. 69(3), p. 825-839, 2004
! J. M. Carcione, Wave fields in real media: wave propagation in anisotropic, anelastic
!  and porous media, Elsevier, p. 124-125, 2007

          ! compute relaxed elastic coefficients from formulas in Carcione 2007 page 125
          lambdal_relaxed_viscoelastic = (lambdal_unrelaxed_elastic + mul_unrelaxed_elastic) / Mu_nu1(i,j,ispec) &
                            - mul_unrelaxed_elastic / Mu_nu2(i,j,ispec)
          mul_relaxed_viscoelastic = mul_unrelaxed_elastic / Mu_nu2(i,j,ispec)
          lambdalplus2mul_relaxed_viscoel = lambdal_relaxed_viscoelastic + TWO*mul_relaxed_viscoelastic

          ! compute the stress using the unrelaxed Lame parameters (Carcione 2007 page 125)
          sigma_xx = (lambdal_unrelaxed_elastic + 2.0 * mul_unrelaxed_elastic)*dux_dxl + lambdal_unrelaxed_elastic*duz_dzl
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
!         sigma_yy = ...  ! it is not zero because of the plane strain formulation, thus it should be computed here
          stop 'pressure calculation not implemented for poroelastic media yet, you should compute sigma_yy here'
          sigma_zz = (lambdal_unrelaxed_elastic + 2.0 * mul_unrelaxed_elastic)*duz_dzl + lambdal_unrelaxed_elastic*dux_dxl

          ! add the memory variables using the relaxed parameters (Carcione 2007 page 125)
          ! beware: there is a bug in  Carcione's equation (2c) of his 1993 paper for sigma_zz, we fixed it in the code below
          e1_sum = 0._CUSTOM_REAL
          e11_sum = 0._CUSTOM_REAL

          do i_sls = 1,N_SLS
            e1_sum = e1_sum + e1(i,j,ispec,i_sls)
            e11_sum = e11_sum + e11(i,j,ispec,i_sls)
          enddo

          sigma_xx = sigma_xx + (lambdal_relaxed_viscoelastic + mul_relaxed_viscoelastic) * e1_sum &
                    + TWO * mul_relaxed_viscoelastic * e11_sum
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
!         sigma_yy = ...  ! it is not zero because of the plane strain formulation, thus it should be computed here
          stop 'pressure calculation not implemented for poroelastic media yet, you should compute sigma_yy here'
          sigma_zz = sigma_zz + (lambdal_relaxed_viscoelastic + mul_relaxed_viscoelastic) * e1_sum &
                    - TWO * mul_relaxed_viscoelastic * e11_sum

        else

          ! no attenuation
          sigma_xx = lambdalplus2mul_G*dux_dxl + lambdal_G*duz_dzl + C_biot*(dwx_dxl + dwz_dzl)
          ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
!         sigma_yy = ...  ! it is not zero because of the plane strain formulation, thus it should be computed here
          stop 'pressure calculation not implemented for poroelastic media yet, you should compute sigma_yy here'
          sigma_zz = lambdalplus2mul_G*duz_dzl + lambdal_G*dux_dxl + C_biot*(dwx_dxl + dwz_dzl)

!         sigmap = C_biot*(dux_dxl + duz_dzl) + M_biot*(dwx_dxl + dwz_dzl)

        endif

        ! store pressure
        ! sigma_yy is not equal to zero in a 2D medium because of the plane strain formulation
        pressure_element(i,j) = - (sigma_xx + sigma_yy + sigma_zz) / 3.d0
!       pressure_element2(i,j) = - sigmap
      enddo
    enddo

! pressure = - Chi_dot_dot if acoustic element
  else if(acoustic(ispec)) then

    do j = 1,NGLLZ
      do i = 1,NGLLX

        iglob = ibool(i,j,ispec)

        ! store pressure
        pressure_element(i,j) = - potential_dot_dot_acoustic(iglob)

      enddo
    enddo

  else if(gravitoacoustic(ispec)) then

    do j = 1,NGLLZ
      do i = 1,NGLLX

        iglob = ibool(i,j,ispec)

        ! store pressure
        pressure_element(i,j) = - potential_dot_dot_gravitoacoustic(iglob)

      enddo
    enddo

  endif ! end of test if acoustic or elastic or gravito element

  end subroutine compute_pressure_one_element


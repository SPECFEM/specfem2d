
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

  subroutine compute_elastic_energy(displ_elastic,veloc_elastic, &
         xix,xiz,gammax,gammaz,jacobian,ibool,elastic,hprime_xx,hprime_zz, &
         nspec,npoin,assign_external_model,it,deltat,t0,kmato,elastcoef,density, &
         vpext,vsext,rhoext,wxgll,wzgll,numat)

! compute kinetic and potential energy in the solid (acoustic elements are excluded)

  implicit none

  include "constants.h"

  integer :: nspec,npoin,numat

  integer :: it
  double precision :: t0,deltat

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  logical, dimension(nspec) :: elastic

  double precision, dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz,jacobian

  integer, dimension(nspec) :: kmato

  logical :: assign_external_model

  double precision, dimension(numat) :: density
  double precision, dimension(4,numat) :: elastcoef
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext,vsext,rhoext

  double precision, dimension(NDIM,npoin) :: displ_elastic,veloc_elastic

! Gauss-Lobatto-Legendre points and weights
  double precision, dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: wzgll

! array with derivatives of Lagrange polynomials
  double precision, dimension(NGLLX,NGLLX) :: hprime_xx
  double precision, dimension(NGLLZ,NGLLZ) :: hprime_zz

! local variables
  integer :: i,j,k,ispec

! spatial derivatives
  double precision :: dux_dxi,dux_dgamma,duz_dxi,duz_dgamma
  double precision :: dux_dxl,duz_dxl,dux_dzl,duz_dzl

! jacobian
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl

  double precision :: kinetic_energy,potential_energy
  double precision :: cpl,csl,rhol,mul_relaxed,lambdal_relaxed,lambdalplus2mul_relaxed

  kinetic_energy = ZERO
  potential_energy = ZERO

! loop over spectral elements
  do ispec = 1,nspec

!---
!--- elastic spectral element
!---
    if(elastic(ispec)) then

! get relaxed elastic parameters of current spectral element
      lambdal_relaxed = elastcoef(1,kmato(ispec))
      mul_relaxed = elastcoef(2,kmato(ispec))
      lambdalplus2mul_relaxed = elastcoef(3,kmato(ispec))
      rhol  = density(kmato(ispec))

! first double loop over GLL points to compute and store gradients
      do j = 1,NGLLZ
        do i = 1,NGLLX

!--- if external medium, get elastic parameters of current grid point
          if(assign_external_model) then
            cpl = vpext(i,j,ispec)
            csl = vsext(i,j,ispec)
            rhol = rhoext(i,j,ispec)
            mul_relaxed = rhol*csl*csl
            lambdal_relaxed = rhol*cpl*cpl - TWO*mul_relaxed
            lambdalplus2mul_relaxed = lambdal_relaxed + TWO*mul_relaxed
          endif

! derivative along x and along z
          dux_dxi = ZERO
          duz_dxi = ZERO

          dux_dgamma = ZERO
          duz_dgamma = ZERO

! first double loop over GLL points to compute and store gradients
! we can merge the two loops because NGLLX == NGLLZ
          do k = 1,NGLLX
            dux_dxi = dux_dxi + displ_elastic(1,ibool(k,j,ispec))*hprime_xx(k,i)
            duz_dxi = duz_dxi + displ_elastic(2,ibool(k,j,ispec))*hprime_xx(k,i)
            dux_dgamma = dux_dgamma + displ_elastic(1,ibool(i,k,ispec))*hprime_zz(k,j)
            duz_dgamma = duz_dgamma + displ_elastic(2,ibool(i,k,ispec))*hprime_zz(k,j)
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

! compute potential energy
          potential_energy = potential_energy + (lambdalplus2mul_relaxed*dux_dxl**2 &
              + lambdalplus2mul_relaxed*duz_dzl**2 &
              + two*lambdal_relaxed*dux_dxl*duz_dzl + mul_relaxed*(dux_dzl + duz_dxl)**2)*wxgll(i)*wzgll(j)*jacobianl

! compute kinetic energy
          kinetic_energy = kinetic_energy + &
              rhol*(veloc_elastic(1,ibool(i,j,ispec))**2 + veloc_elastic(2,ibool(i,j,ispec))**2) *wxgll(i)*wzgll(j)*jacobianl

        enddo
      enddo

    endif

  enddo

! do not forget to divide by two at the end
  kinetic_energy = kinetic_energy / TWO
  potential_energy = potential_energy / TWO

! save kinetic, potential and total energy for this time step in external file
  write(IENERGY,*) sngl(dble(it-1)*deltat - t0),sngl(kinetic_energy), &
                     sngl(potential_energy),sngl(kinetic_energy + potential_energy)

  end subroutine compute_elastic_energy


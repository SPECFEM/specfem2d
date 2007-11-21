
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!  Main authors: Dimitri Komatitsch, Nicolas Le Goff and Roland Martin
!                 University of Pau and CNRS, France
!
!                         (c) November 2007
!
!========================================================================

  subroutine compute_arrays_source(ispec_selected_source,xi_source,gamma_source,sourcearray, &
             Mxx,Mzz,Mxz,xix,xiz,gammax,gammaz,xigll,zigll,nspec)

  implicit none

  include "constants.h"

  integer ispec_selected_source
  integer nspec

  double precision xi_source,gamma_source
  double precision Mxx,Mzz,Mxz

  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: xix,xiz,gammax,gammaz

  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  double precision xixd,xizd,gammaxd,gammazd

! Gauss-Lobatto-Legendre points of integration and weights
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

! source arrays
  double precision, dimension(NGLLX,NGLLZ) :: G11,G13,G31,G33
  double precision, dimension(NGLLX) :: hxis,hpxis
  double precision, dimension(NGLLZ) :: hgammas,hpgammas

  integer k,m
  integer ir,iv

! calculate G_ij for general source location
! the source does not necessarily correspond to a Gauss-Lobatto point
  do m=1,NGLLZ
      do k=1,NGLLX

        xixd    = xix(k,m,ispec_selected_source)
        xizd    = xiz(k,m,ispec_selected_source)
        gammaxd = gammax(k,m,ispec_selected_source)
        gammazd = gammaz(k,m,ispec_selected_source)

        G11(k,m) = Mxx*xixd+Mxz*xizd
        G13(k,m) = Mxx*gammaxd+Mxz*gammazd
        G31(k,m) = Mxz*xixd+Mzz*xizd
        G33(k,m) = Mxz*gammaxd+Mzz*gammazd

!!!!        G21(k,m) = Mxy*xixd+Myz*xizd
!!!!        G23(k,m) = Mxy*gammaxd+Myz*gammazd

      enddo
  enddo

! compute Lagrange polynomials at the source location
  call lagrange_any(xi_source,NGLLX,xigll,hxis,hpxis)
  call lagrange_any(gamma_source,NGLLZ,zigll,hgammas,hpgammas)

! calculate source array
  do m=1,NGLLZ
    do k=1,NGLLX

      sourcearray(:,k,m) = ZERO

      do iv=1,NGLLZ
        do ir=1,NGLLX

          sourcearray(1,k,m) = sourcearray(1,k,m) + hxis(ir)*hgammas(iv) &
                                 *(G11(ir,iv)*hpxis(k)*hgammas(m) &
                                 +G13(ir,iv)*hxis(k)*hpgammas(m))

!        sourcearray(2,k,m) = sourcearray(2,k,m) + hxis(ir)*hgammas(iv) &
!                               *(G21(ir,iv)*hpxis(k)*hgammas(m) &
!                               +G23(ir,iv)*hxis(k)*hpgammas(m))

          sourcearray(2,k,m) = sourcearray(2,k,m) + hxis(ir)*hgammas(iv) &
                                 *(G31(ir,iv)*hpxis(k)*hgammas(m) &
                                 +G33(ir,iv)*hxis(k)*hpgammas(m))

        enddo
      enddo

    enddo
  enddo

  end subroutine compute_arrays_source


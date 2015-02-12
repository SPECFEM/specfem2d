
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

! ------------------------------------------------------------------------------------------------------


  subroutine compute_arrays_adj_source(xi_rec,gamma_rec)

 use constants
 use specfem_par

 implicit none



  double precision xi_rec, gamma_rec

  integer icomp, itime, i, k
  double precision :: junk
  character(len=3) :: comp(3)

  call lagrange_any(xi_rec,NGLLX,xigll,hxir,hpxir)
  call lagrange_any(gamma_rec,NGLLZ,zigll,hgammar,hpgammar)

  adj_sourcearray(:,:,:,:) = 0.


  if(seismotype == 1 .or. seismotype == 2 .or. seismotype == 3) then

  comp = (/"BXX","BXY","BXZ"/)

  do icomp = 1,3

    filename = 'SEM/'//trim(adj_source_file) // '.'// comp(icomp) // '.adj'
    open(unit = IIN, file = trim(filename), iostat = ios)
    if (ios /= 0) call exit_MPI(' file '//trim(filename)//'does not exist')

    do itime = 1, NSTEP
      read(IIN,*) junk, adj_src_s(itime,icomp)
    enddo
    close(IIN)

  enddo

  source_adjointe(irec_local,:,2) = adj_src_s(:,3)

  else if (seismotype == 4 ) then

    filename = 'SEM/'//trim(adj_source_file) // '.PRE.adj'
    open(unit = IIN, file = trim(filename), iostat = ios)
    if (ios /= 0) call exit_MPI(' file '//trim(filename)//'does not exist')

    do itime = 1, NSTEP
      read(IIN,*) junk, adj_src_s(itime,1)
    enddo
    close(IIN)

  else if (seismotype == 6 ) then

    filename = 'SEM/'//trim(adj_source_file) // '.POT.adj'
    open(unit = IIN, file = trim(filename), iostat = ios)
    if (ios /= 0) call exit_MPI(' file '//trim(filename)//'does not exist')

    do itime = 1, NSTEP
      read(IIN,*) junk, adj_src_s(itime,1)
    enddo
    close(IIN)

  endif

  source_adjointe(irec_local,:,1) = adj_src_s(:,1)

  do k = 1, NGLLZ
      do i = 1, NGLLX
        adj_sourcearray(:,:,i,k) = sngl(hxir(i) * hgammar(k)) * adj_src_s(:,:)
      enddo
  enddo


end subroutine compute_arrays_adj_source

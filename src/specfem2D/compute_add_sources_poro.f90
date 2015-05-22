
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
!=====================================================================

! for poro solver

  subroutine compute_add_sources_poro(accels_poroelastic,accelw_poroelastic,it,i_stage)

  use specfem_par, only: poroelastic,nglob_poroelastic, &
                         NSOURCES,source_type,anglesource,source_time_function, &
                         is_proc_source,ispec_selected_source, &
                         hxis_store,hgammas_store,ibool, &
                         porosity,tortuosity,density,kmato
  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic) :: accels_poroelastic,accelw_poroelastic

  integer :: it,i_stage

  !local variables
  integer :: i_source,i,j,iglob
  double precision :: hlagrange,phil,tortl,rhol_s,rhol_f,rhol_bar

  do i_source=1,NSOURCES
    ! if this processor core carries the source and the source element is elastic
    if (is_proc_source(i_source) == 1 .and. poroelastic(ispec_selected_source(i_source))) then
      phil = porosity(kmato(ispec_selected_source(i_source)))
      tortl = tortuosity(kmato(ispec_selected_source(i_source)))
      rhol_s = density(1,kmato(ispec_selected_source(i_source)))
      rhol_f = density(2,kmato(ispec_selected_source(i_source)))
      rhol_bar = (1._CUSTOM_REAL - phil)*rhol_s + phil*rhol_f

      ! collocated force
      if( source_type(i_source) == 1 ) then
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec_selected_source(i_source))
            hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
            ! s
            accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) - hlagrange * &
                   (1._CUSTOM_REAL - phil/tortl)*sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)
            accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + hlagrange * &
                   (1._CUSTOM_REAL - phil/tortl)*cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)
            ! w
            accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) - hlagrange * &
                   (1._CUSTOM_REAL - rhol_f/rhol_bar)*sin(anglesource(i_source))*source_time_function(i_source,it,i_stage)
            accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + hlagrange * &
                   (1._CUSTOM_REAL - rhol_f/rhol_bar)*cos(anglesource(i_source))*source_time_function(i_source,it,i_stage)
          enddo
        enddo
      endif
    endif ! if this processor core carries the source and the source element is elastic
  enddo ! do i_source=1,NSOURCES

 end subroutine compute_add_sources_poro


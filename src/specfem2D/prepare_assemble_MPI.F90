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

!
! This file contains subroutines related to assembling (of the mass matrix, potential_dot_dot and
! accel_elastic, accels_poroelastic, accelw_poroelastic).
! These subroutines are for the most part not used in the sequential version.
!

!-----------------------------------------------
! Determines the points that are on the interfaces with other partitions, to help
! build the communication buffers, and determines which elements are considered 'inner'
! (no points in common with other partitions) and 'outer' (at least one point in common
! with neighboring partitions).
! We have both acoustic and (poro)elastic buffers, for coupling between acoustic and (poro)elastic elements
! led us to have two sets of communications.
!-----------------------------------------------
  subroutine prepare_assemble_MPI()

  use constants, only: NGLLX,NGLLZ

  use specfem_par, only: ibool, knods, ngnod, nglob, &
    ispec_is_elastic, ispec_is_poroelastic, ispec_is_acoustic

  use specfem_par, only: ninterface, my_nelmnts_neighbors, my_interfaces, &
    nibool_interfaces_ext_mesh, ibool_interfaces_ext_mesh_init

  use specfem_par, only: NPROC, &
    ibool_interfaces_acoustic, ibool_interfaces_elastic, &
    ibool_interfaces_poroelastic, &
    nibool_interfaces_acoustic, nibool_interfaces_elastic, &
    nibool_interfaces_poroelastic, &
    inum_interfaces_acoustic, inum_interfaces_elastic, &
    inum_interfaces_poroelastic, &
    ninterface_acoustic, ninterface_elastic, ninterface_poroelastic

  implicit none

  ! local parameters
  integer  :: iinterface
  integer  :: ispec_interface

  logical, dimension(nglob)  :: mask_ibool_acoustic
  logical, dimension(nglob)  :: mask_ibool_elastic
  logical, dimension(nglob)  :: mask_ibool_poroelastic
  logical, dimension(nglob)  :: mask_ibool_ext_mesh

  integer  :: ixmin, ixmax, izmin, izmax, ix, iz
  integer, dimension(ngnod)  :: n
  integer  :: e1, e2, itype, ispec, k, sens, iglob
  integer  :: nglob_interface_acoustic
  integer  :: nglob_interface_elastic
  integer  :: nglob_interface_poroelastic
  integer :: npoin_interface_ext_mesh

  ! checks if anything to do
  if (NPROC <= 1) return

  ! initializes
  ! for all domains
  nibool_interfaces_ext_mesh(:) = 0
  ibool_interfaces_ext_mesh_init(:,:) = 0

  ! only specific domain interfaces
  ibool_interfaces_acoustic(:,:) = 0
  nibool_interfaces_acoustic(:) = 0

  ibool_interfaces_elastic(:,:) = 0
  nibool_interfaces_elastic(:) = 0

  ibool_interfaces_poroelastic(:,:) = 0
  nibool_interfaces_poroelastic(:) = 0

  do iinterface = 1, ninterface
    ! initializes interface point counters
    npoin_interface_ext_mesh = 0
    mask_ibool_ext_mesh(:) = .false.

    nglob_interface_acoustic = 0
    nglob_interface_elastic = 0
    nglob_interface_poroelastic = 0

    mask_ibool_acoustic(:) = .false.
    mask_ibool_elastic(:) = .false.
    mask_ibool_poroelastic(:) = .false.

    do ispec_interface = 1, my_nelmnts_neighbors(iinterface)
      ! element id
      ispec = my_interfaces(1,ispec_interface,iinterface)

      ! type of interface: 1 = common point, 2 = common edge
      itype = my_interfaces(2,ispec_interface,iinterface)

      ! element control node ids
      do k = 1, ngnod
        n(k) = knods(k,ispec)
      enddo

      ! common node ids
      e1 = my_interfaces(3,ispec_interface,iinterface)
      e2 = my_interfaces(4,ispec_interface,iinterface)

      call get_edge(ngnod, n, itype, e1, e2, ixmin, ixmax, izmin, izmax, sens)

      ! sets interface points (all material domains)
      do iz = izmin, izmax, sens
        do ix = ixmin, ixmax, sens
          ! global index
          iglob = ibool(ix,iz,ispec)

          if (.not. mask_ibool_ext_mesh(iglob)) then
            ! masks point as being accounted for
            mask_ibool_ext_mesh(iglob) = .true.
            ! adds point to interface
            npoin_interface_ext_mesh = npoin_interface_ext_mesh + 1
            ibool_interfaces_ext_mesh_init(npoin_interface_ext_mesh,iinterface) = iglob
          endif
        enddo
      enddo

      ! sets interface points for specific domains
      do iz = izmin, izmax, sens
        do ix = ixmin, ixmax, sens
          ! global index
          iglob = ibool(ix,iz,ispec)

          ! checks to which material this common interface belongs
          if (ispec_is_elastic(ispec)) then
            ! elastic element
            if (.not. mask_ibool_elastic(iglob)) then
              mask_ibool_elastic(iglob) = .true.
              nglob_interface_elastic = nglob_interface_elastic + 1
              ibool_interfaces_elastic(nglob_interface_elastic,iinterface) = iglob
            endif

          else if (ispec_is_poroelastic(ispec)) then
            ! poroelastic element
            if (.not. mask_ibool_poroelastic(iglob)) then
              mask_ibool_poroelastic(iglob) = .true.
              nglob_interface_poroelastic = nglob_interface_poroelastic + 1
              ibool_interfaces_poroelastic(nglob_interface_poroelastic,iinterface) = iglob
            endif

          else if (ispec_is_acoustic(ispec)) then
            ! acoustic element
            if (.not. mask_ibool_acoustic(iglob)) then
              mask_ibool_acoustic(iglob) = .true.
              nglob_interface_acoustic = nglob_interface_acoustic + 1
              ibool_interfaces_acoustic(nglob_interface_acoustic,iinterface) = iglob
            endif

          else
            call stop_the_code('Invalid element type found in prepare_assemble_MPI() routine')
          endif

        enddo
      enddo

    enddo

    ! stores total number of (global) points on this MPI interface
    nibool_interfaces_ext_mesh(iinterface) = npoin_interface_ext_mesh

    ! stores counters for interface points
    nibool_interfaces_acoustic(iinterface) = nglob_interface_acoustic
    nibool_interfaces_elastic(iinterface) = nglob_interface_elastic
    nibool_interfaces_poroelastic(iinterface) = nglob_interface_poroelastic

  enddo

  ! sets number of interfaces for each material domain
  ninterface_acoustic = 0
  ninterface_elastic =  0
  ninterface_poroelastic =  0

  ! loops over all MPI interfaces
  do iinterface = 1, ninterface
    ! sets acoustic MPI interface (local) indices in range [1,ninterface_acoustic]
    if (nibool_interfaces_acoustic(iinterface) > 0) then
      ninterface_acoustic = ninterface_acoustic + 1
      inum_interfaces_acoustic(ninterface_acoustic) = iinterface
    endif
    ! elastic
    if (nibool_interfaces_elastic(iinterface) > 0) then
      ninterface_elastic = ninterface_elastic + 1
      inum_interfaces_elastic(ninterface_elastic) = iinterface
    endif
    ! poroelastic
    if (nibool_interfaces_poroelastic(iinterface) > 0) then
      ninterface_poroelastic = ninterface_poroelastic + 1
      inum_interfaces_poroelastic(ninterface_poroelastic) = iinterface
    endif
  enddo

  end subroutine prepare_assemble_MPI


!-----------------------------------------------
! Get the points (ixmin, ixmax, izmin and izmax) on an node/edge for one element.
! 'sens' is used to have DO loops with increment equal to 'sens' (-/+1).
!-----------------------------------------------
  subroutine get_edge ( ngnod, n, itype, e1, e2, ixmin, ixmax, izmin, izmax, sens )

  use constants, only: NGLLX,NGLLZ

  implicit none

  integer, intent(in)  :: ngnod
  integer, dimension(ngnod), intent(in)  :: n
  integer, intent(in)  :: itype, e1, e2
  integer, intent(out)  :: ixmin, ixmax, izmin, izmax
  integer, intent(out)  :: sens

  if (itype == 1) then

    ! common single point

    ! checks which corner point is given
    if (e1 == n(1)) then
        ixmin = 1
        ixmax = 1
        izmin = 1
        izmax = 1
    endif
    if (e1 == n(2)) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = 1
        izmax = 1
    endif
    if (e1 == n(3)) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = NGLLZ
        izmax = NGLLZ
    endif
    if (e1 == n(4)) then
        ixmin = 1
        ixmax = 1
        izmin = NGLLZ
        izmax = NGLLZ
    endif
    sens = 1

  else if (itype == 2) then

    ! common edge

    ! checks which edge and corner points are given
    if (e1 == n(1)) then
        ixmin = 1
        izmin = 1
        if (e2 == n(2)) then
           ixmax = NGLLX
           izmax = 1
           sens = 1
        endif
        if (e2 == n(4)) then
           ixmax = 1
           izmax = NGLLZ
           sens = 1
        endif
     endif
     if (e1 == n(2)) then
        ixmin = NGLLX
        izmin = 1
        if (e2 == n(3)) then
           ixmax = NGLLX
           izmax = NGLLZ
           sens = 1
        endif
        if (e2 == n(1)) then
           ixmax = 1
           izmax = 1
           sens = -1
        endif
     endif
     if (e1 == n(3)) then
        ixmin = NGLLX
        izmin = NGLLZ
        if (e2 == n(4)) then
           ixmax = 1
           izmax = NGLLZ
           sens = -1
        endif
        if (e2 == n(2)) then
           ixmax = NGLLX
           izmax = 1
           sens = -1
        endif
     endif
     if (e1 == n(4)) then
        ixmin = 1
        izmin = NGLLZ
        if (e2 == n(1)) then
           ixmax = 1
           izmax = 1
           sens = -1
        endif
        if (e2 == n(3)) then
           ixmax = NGLLX
           izmax = NGLLZ
           sens = 1
        endif
     endif

  else

    call stop_the_code('ERROR get_edge unknown type')

  endif

  end subroutine get_edge

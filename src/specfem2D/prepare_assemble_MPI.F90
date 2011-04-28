
!========================================================================
!
!                   S P E C F E M 2 D  Version 6 . 2
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
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

!
! This file contains subroutines related to assembling (of the mass matrix, potential_dot_dot and
! accel_elastic, accels_poroelastic, accelw_poroelastic).
! These subroutines are for the most part not used in the sequential version.
!

#ifdef USE_MPI

!-----------------------------------------------
! Determines the points that are on the interfaces with other partitions, to help
! build the communication buffers, and determines which elements are considered 'inner'
! (no points in common with other partitions) and 'outer' (at least one point in common
! with neighbouring partitions).
! We have both acoustic and (poro)elastic buffers, for coupling between acoustic and (poro)elastic elements
! led us to have two sets of communications.
!-----------------------------------------------
  subroutine prepare_assemble_MPI(nspec,ibool,knods, ngnod,nglob, elastic, poroelastic, &
                                ninterface, max_interface_size, &
                                my_nelmnts_neighbours, my_interfaces, &
                                ibool_interfaces_acoustic, ibool_interfaces_elastic, &
                                ibool_interfaces_poroelastic, &
                                nibool_interfaces_acoustic, nibool_interfaces_elastic, &
                                nibool_interfaces_poroelastic, &
                                inum_interfaces_acoustic, inum_interfaces_elastic, &
                                inum_interfaces_poroelastic, &
                                ninterface_acoustic, ninterface_elastic, ninterface_poroelastic, &
                                mask_ispec_inner_outer )

  implicit none

  include 'constants.h'

  integer, intent(in)  :: nspec, nglob, ngnod
  logical, dimension(nspec), intent(in)  :: elastic, poroelastic
  integer, dimension(ngnod,nspec), intent(in)  :: knods
  integer, dimension(NGLLX,NGLLZ,nspec), intent(in)  :: ibool

  integer  :: ninterface
  integer  :: max_interface_size
  integer, dimension(ninterface)  :: my_nelmnts_neighbours
  integer, dimension(4,max_interface_size,ninterface)  :: my_interfaces
  integer, dimension(NGLLX*max_interface_size,ninterface)  :: &
       ibool_interfaces_acoustic,ibool_interfaces_elastic,ibool_interfaces_poroelastic
  integer, dimension(ninterface)  :: &
       nibool_interfaces_acoustic,nibool_interfaces_elastic,nibool_interfaces_poroelastic

  integer, dimension(ninterface), intent(out)  :: &
       inum_interfaces_acoustic, inum_interfaces_elastic, inum_interfaces_poroelastic
  integer, intent(out)  :: ninterface_acoustic, ninterface_elastic, ninterface_poroelastic

  logical, dimension(nspec), intent(inout)  :: mask_ispec_inner_outer

  ! local parameters
  integer  :: num_interface
  integer  :: ispec_interface
  logical, dimension(nglob)  :: mask_ibool_acoustic
  logical, dimension(nglob)  :: mask_ibool_elastic
  logical, dimension(nglob)  :: mask_ibool_poroelastic
  integer  :: ixmin, ixmax, izmin, izmax, ix, iz
  integer, dimension(ngnod)  :: n
  integer  :: e1, e2, itype, ispec, k, sens, iglob
  integer  :: nglob_interface_acoustic
  integer  :: nglob_interface_elastic
  integer  :: nglob_interface_poroelastic

  ! initializes
  ibool_interfaces_acoustic(:,:) = 0
  nibool_interfaces_acoustic(:) = 0
  ibool_interfaces_elastic(:,:) = 0
  nibool_interfaces_elastic(:) = 0
  ibool_interfaces_poroelastic(:,:) = 0
  nibool_interfaces_poroelastic(:) = 0

  do num_interface = 1, ninterface
    ! initializes interface point counters
    nglob_interface_acoustic = 0
    nglob_interface_elastic = 0
    nglob_interface_poroelastic = 0
    mask_ibool_acoustic(:) = .false.
    mask_ibool_elastic(:) = .false.
    mask_ibool_poroelastic(:) = .false.

    do ispec_interface = 1, my_nelmnts_neighbours(num_interface)
      ! element id
      ispec = my_interfaces(1,ispec_interface,num_interface)
      ! type of interface: 1 = common point, 2 = common edge
      itype = my_interfaces(2,ispec_interface,num_interface)
      ! element control node ids
      do k = 1, ngnod
        n(k) = knods(k,ispec)
      end do
      ! common node ids
      e1 = my_interfaces(3,ispec_interface,num_interface)
      e2 = my_interfaces(4,ispec_interface,num_interface)

      call get_edge(ngnod, n, itype, e1, e2, ixmin, ixmax, izmin, izmax, sens)

      do iz = izmin, izmax, sens
        do ix = ixmin, ixmax, sens
          ! global index
          iglob = ibool(ix,iz,ispec)

          ! checks to which material this common interface belongs
          if ( elastic(ispec) ) then
            ! elastic element
            if(.not. mask_ibool_elastic(iglob)) then
              mask_ibool_elastic(iglob) = .true.
              nglob_interface_elastic = nglob_interface_elastic + 1
              ibool_interfaces_elastic(nglob_interface_elastic,num_interface) = iglob
            end if
          else if ( poroelastic(ispec) ) then
            ! poroelastic element
            if(.not. mask_ibool_poroelastic(iglob)) then
              mask_ibool_poroelastic(iglob) = .true.
              nglob_interface_poroelastic = nglob_interface_poroelastic + 1
              ibool_interfaces_poroelastic(nglob_interface_poroelastic,num_interface) = iglob
            end if
          else
            ! acoustic element
            if(.not. mask_ibool_acoustic(iglob)) then
              mask_ibool_acoustic(iglob) = .true.
              nglob_interface_acoustic = nglob_interface_acoustic + 1
              ibool_interfaces_acoustic(nglob_interface_acoustic,num_interface) = iglob
            end if
          end if
        end do
      end do

    end do

    ! stores counters for interface points
    nibool_interfaces_acoustic(num_interface) = nglob_interface_acoustic
    nibool_interfaces_elastic(num_interface) = nglob_interface_elastic
    nibool_interfaces_poroelastic(num_interface) = nglob_interface_poroelastic

    ! sets inner/outer element flags
    do ispec = 1, nspec
      do iz = 1, NGLLZ
        do ix = 1, NGLLX
          if ( mask_ibool_acoustic(ibool(ix,iz,ispec)) &
            .or. mask_ibool_elastic(ibool(ix,iz,ispec)) &
            .or. mask_ibool_poroelastic(ibool(ix,iz,ispec)) ) then
               mask_ispec_inner_outer(ispec) = .true.
          endif

        enddo
      enddo
    enddo

  end do

  ! sets number of interfaces for each material domain
  ninterface_acoustic = 0
  ninterface_elastic =  0
  ninterface_poroelastic =  0

  ! loops over all MPI interfaces
  do num_interface = 1, ninterface
    ! sets acoustic MPI interface (local) indices in range [1,ninterface_acoustic]
    if ( nibool_interfaces_acoustic(num_interface) > 0 ) then
      ninterface_acoustic = ninterface_acoustic + 1
      inum_interfaces_acoustic(ninterface_acoustic) = num_interface
    end if
    ! elastic
    if ( nibool_interfaces_elastic(num_interface) > 0 ) then
      ninterface_elastic = ninterface_elastic + 1
      inum_interfaces_elastic(ninterface_elastic) = num_interface
    end if
    ! poroelastic
    if ( nibool_interfaces_poroelastic(num_interface) > 0 ) then
      ninterface_poroelastic = ninterface_poroelastic + 1
      inum_interfaces_poroelastic(ninterface_poroelastic) = num_interface
    end if
  end do

  end subroutine prepare_assemble_MPI


!-----------------------------------------------
! Get the points (ixmin, ixmax, izmin and izmax) on an node/edge for one element.
! 'sens' is used to have DO loops with increment equal to 'sens' (-/+1).
!-----------------------------------------------
  subroutine get_edge ( ngnod, n, itype, e1, e2, ixmin, ixmax, izmin, izmax, sens )

  implicit none

  include "constants.h"

  integer, intent(in)  :: ngnod
  integer, dimension(ngnod), intent(in)  :: n
  integer, intent(in)  :: itype, e1, e2
  integer, intent(out)  :: ixmin, ixmax, izmin, izmax
  integer, intent(out)  :: sens

  if ( itype == 1 ) then

    ! common single point

    ! checks which corner point is given
    if ( e1 == n(1) ) then
        ixmin = 1
        ixmax = 1
        izmin = 1
        izmax = 1
    end if
    if ( e1 == n(2) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = 1
        izmax = 1
    end if
    if ( e1 == n(3) ) then
        ixmin = NGLLX
        ixmax = NGLLX
        izmin = NGLLZ
        izmax = NGLLZ
    end if
    if ( e1 == n(4) ) then
        ixmin = 1
        ixmax = 1
        izmin = NGLLZ
        izmax = NGLLZ
    end if
    sens = 1

  else if( itype == 2 ) then

    ! common edge

    ! checks which edge and corner points are given
    if ( e1 ==  n(1) ) then
        ixmin = 1
        izmin = 1
        if ( e2 == n(2) ) then
           ixmax = NGLLX
           izmax = 1
           sens = 1
        end if
        if ( e2 == n(4) ) then
           ixmax = 1
           izmax = NGLLZ
           sens = 1
        end if
     end if
     if ( e1 == n(2) ) then
        ixmin = NGLLX
        izmin = 1
        if ( e2 == n(3) ) then
           ixmax = NGLLX
           izmax = NGLLZ
           sens = 1
        end if
        if ( e2 == n(1) ) then
           ixmax = 1
           izmax = 1
           sens = -1
        end if
     end if
     if ( e1 == n(3) ) then
        ixmin = NGLLX
        izmin = NGLLZ
        if ( e2 == n(4) ) then
           ixmax = 1
           izmax = NGLLZ
           sens = -1
        end if
        if ( e2 == n(2) ) then
           ixmax = NGLLX
           izmax = 1
           sens = -1
        end if
     end if
     if ( e1 == n(4) ) then
        ixmin = 1
        izmin = NGLLZ
        if ( e2 == n(1) ) then
           ixmax = 1
           izmax = 1
           sens = -1
        end if
        if ( e2 == n(3) ) then
           ixmax = NGLLX
           izmax = NGLLZ
           sens = 1
        end if
     end if

  else

    call exit_MPI('ERROR get_edge unknown type')

  end if

  end subroutine get_edge

#endif

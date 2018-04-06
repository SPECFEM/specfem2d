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

  subroutine get_coupling_edges()

! handles all domain coupling

  use constants, only: IMAIN,IBOTTOM,IRIGHT,ITOP,ILEFT,TINYVAL
  use specfem_par

  implicit none

  ! local parameters
  integer :: ispec_poroelastic,iedge_poroelastic
  integer :: ispec_acoustic,ispec_elastic
  integer :: iedge_acoustic,iedge_elastic
  integer :: inum,ipoin1D,iglob2
  integer :: i,j,iglob

  ! determine if coupled fluid-solid simulation
  coupled_acoustic_elastic = any_acoustic .and. any_elastic
  coupled_acoustic_poro = any_acoustic .and. any_poroelastic
  coupled_elastic_poro = any_elastic .and. any_poroelastic

  ! fluid/solid (elastic) edge detection
  ! the two elements (fluid and solid) forming an edge are already known (computed in meshfem2D),
  ! the common nodes forming the edge are computed here
  if (coupled_acoustic_elastic) then

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Mixed acoustic/elastic simulation'
      write(IMAIN,*)
      write(IMAIN,*) '  Beginning of fluid/solid edge detection'
      call flush_IMAIN()
    endif

    ! define the edges of a given element
    i_begin(IBOTTOM) = 1
    j_begin(IBOTTOM) = 1
    i_end(IBOTTOM) = NGLLX
    j_end(IBOTTOM) = 1

    i_begin(IRIGHT) = NGLLX
    j_begin(IRIGHT) = 1
    i_end(IRIGHT) = NGLLX
    j_end(IRIGHT) = NGLLZ

    i_begin(ITOP) = NGLLX
    j_begin(ITOP) = NGLLZ
    i_end(ITOP) = 1
    j_end(ITOP) = NGLLZ

    i_begin(ILEFT) = 1
    j_begin(ILEFT) = NGLLZ
    i_end(ILEFT) = 1
    j_end(ILEFT) = 1

    ! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = ipoin1D
      ivalue_inverse(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,IBOTTOM) = 1
      jvalue_inverse(ipoin1D,IBOTTOM) = 1

      ivalue(ipoin1D,IRIGHT) = NGLLX
      ivalue_inverse(ipoin1D,IRIGHT) = NGLLX
      jvalue(ipoin1D,IRIGHT) = ipoin1D
      jvalue_inverse(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1

      ivalue(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,ITOP) = ipoin1D
      jvalue(ipoin1D,ITOP) = NGLLZ
      jvalue_inverse(ipoin1D,ITOP) = NGLLZ

      ivalue(ipoin1D,ILEFT) = 1
      ivalue_inverse(ipoin1D,ILEFT) = 1
      jvalue(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,ILEFT) = ipoin1D

    enddo

    do inum = 1, num_fluid_solid_edges
       ispec_acoustic = fluid_solid_acoustic_ispec(inum)
       ispec_elastic = fluid_solid_elastic_ispec(inum)

        ! one element must be acoustic and the other must be elastic
        if (ispec_acoustic /= ispec_elastic .and. .not. ispec_is_elastic(ispec_acoustic) .and. &
             .not. ispec_is_poroelastic(ispec_acoustic) .and. ispec_is_elastic(ispec_elastic)) then

          ! loop on the four edges of the two elements
          do iedge_acoustic = 1,NEDGES
            do iedge_elastic = 1,NEDGES
              ! store the matching topology if the two edges match in inverse order
              if (ibool(i_begin(iedge_acoustic),j_begin(iedge_acoustic),ispec_acoustic) == &
                  ibool(i_end(iedge_elastic),j_end(iedge_elastic),ispec_elastic) .and. &
                  ibool(i_end(iedge_acoustic),j_end(iedge_acoustic),ispec_acoustic) == &
                  ibool(i_begin(iedge_elastic),j_begin(iedge_elastic),ispec_elastic)) then
                 fluid_solid_acoustic_iedge(inum) = iedge_acoustic
                 fluid_solid_elastic_iedge(inum) = iedge_elastic
!                  print *,'edge ',iedge_acoustic,' of acoustic element ',ispec_acoustic, &
!                          ' is in contact with edge ',iedge_elastic,' of elastic element ',ispec_elastic
              endif

            enddo
          enddo

       endif

    enddo

    ! make sure fluid/solid (elastic) matching has been perfectly detected: check that the grid points
    ! have the same physical coordinates
    ! loop on all the coupling edges

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  Checking fluid/solid edge topology...'
      call flush_IMAIN()
    endif

    do inum = 1,num_fluid_solid_edges

      ! get the edge of the acoustic element
      ispec_acoustic = fluid_solid_acoustic_ispec(inum)
      iedge_acoustic = fluid_solid_acoustic_iedge(inum)

      ! get the corresponding edge of the elastic element
      ispec_elastic = fluid_solid_elastic_ispec(inum)
      iedge_elastic = fluid_solid_elastic_iedge(inum)

      ! implement 1D coupling along the edge
      do ipoin1D = 1,NGLLX
        ! get point values for the elastic side, which matches our side in the inverse direction
        i = ivalue_inverse(ipoin1D,iedge_elastic)
        j = jvalue_inverse(ipoin1D,iedge_elastic)
        iglob = ibool(i,j,ispec_elastic)

        ! get point values for the acoustic side
        i = ivalue(ipoin1D,iedge_acoustic)
        j = jvalue(ipoin1D,iedge_acoustic)
        iglob2 = ibool(i,j,ispec_acoustic)

        ! if distance between the two points is not negligible, there is an error, since it should be zero
        if (sqrt((coord(1,iglob) - coord(1,iglob2))**2 + (coord(2,iglob) - coord(2,iglob2))**2) > TINYVAL) &
            call exit_MPI(myrank,'Error in fluid/solid coupling buffer')
      enddo

    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  End of fluid/solid edge detection'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  endif

  ! fluid/solid (poroelastic) edge detection
  ! the two elements (fluid and solid) forming an edge are already known (computed in meshfem2D),
  ! the common nodes forming the edge are computed here
  if (coupled_acoustic_poro) then

    ! checks
    if (ATTENUATION_VISCOACOUSTIC .or. ATTENUATION_PORO_FLUID_PART) &
      call stop_the_code('Attenuation not supported for mixed acoustic/poroelastic simulations')

    ! user output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Mixed acoustic/poroelastic simulation'
      write(IMAIN,*)
      write(IMAIN,*) '  Beginning of fluid/solid (poroelastic) edge detection'
      call flush_IMAIN()
    endif

    ! define the edges of a given element
    i_begin(IBOTTOM) = 1
    j_begin(IBOTTOM) = 1
    i_end(IBOTTOM) = NGLLX
    j_end(IBOTTOM) = 1

    i_begin(IRIGHT) = NGLLX
    j_begin(IRIGHT) = 1
    i_end(IRIGHT) = NGLLX
    j_end(IRIGHT) = NGLLZ

    i_begin(ITOP) = NGLLX
    j_begin(ITOP) = NGLLZ
    i_end(ITOP) = 1
    j_end(ITOP) = NGLLZ

    i_begin(ILEFT) = 1
    j_begin(ILEFT) = NGLLZ
    i_end(ILEFT) = 1
    j_end(ILEFT) = 1

    ! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = ipoin1D
      ivalue_inverse(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,IBOTTOM) = 1
      jvalue_inverse(ipoin1D,IBOTTOM) = 1

      ivalue(ipoin1D,IRIGHT) = NGLLX
      ivalue_inverse(ipoin1D,IRIGHT) = NGLLX
      jvalue(ipoin1D,IRIGHT) = ipoin1D
      jvalue_inverse(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1

      ivalue(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,ITOP) = ipoin1D
      jvalue(ipoin1D,ITOP) = NGLLZ
      jvalue_inverse(ipoin1D,ITOP) = NGLLZ

      ivalue(ipoin1D,ILEFT) = 1
      ivalue_inverse(ipoin1D,ILEFT) = 1
      jvalue(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,ILEFT) = ipoin1D

    enddo

    do inum = 1, num_fluid_poro_edges
       ispec_acoustic = fluid_poro_acoustic_ispec(inum)
       ispec_poroelastic = fluid_poro_poroelastic_ispec(inum)

        ! one element must be acoustic and the other must be poroelastic
        if (ispec_acoustic /= ispec_poroelastic .and. .not. ispec_is_poroelastic(ispec_acoustic) .and. &
                 .not. ispec_is_elastic(ispec_acoustic) .and. ispec_is_poroelastic(ispec_poroelastic)) then

          ! loop on the four edges of the two elements
          do iedge_acoustic = 1,NEDGES
            do iedge_poroelastic = 1,NEDGES

              ! store the matching topology if the two edges match in inverse order
              if (ibool(i_begin(iedge_acoustic),j_begin(iedge_acoustic),ispec_acoustic) == &
                   ibool(i_end(iedge_poroelastic),j_end(iedge_poroelastic),ispec_poroelastic) .and. &
                   ibool(i_end(iedge_acoustic),j_end(iedge_acoustic),ispec_acoustic) == &
                   ibool(i_begin(iedge_poroelastic),j_begin(iedge_poroelastic),ispec_poroelastic)) then
                 fluid_poro_acoustic_iedge(inum) = iedge_acoustic
                 fluid_poro_poroelastic_iedge(inum) = iedge_poroelastic
                endif

             enddo
          enddo

       endif

    enddo


    ! make sure fluid/solid (poroelastic) matching has been perfectly detected: check that the grid points
    ! have the same physical coordinates
    ! loop on all the coupling edges

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  Checking fluid/solid (poroelastic) edge topology'
      call flush_IMAIN()
    endif

    do inum = 1,num_fluid_poro_edges

      ! get the edge of the acoustic element
      ispec_acoustic = fluid_poro_acoustic_ispec(inum)
      iedge_acoustic = fluid_poro_acoustic_iedge(inum)

      ! get the corresponding edge of the poroelastic element
      ispec_poroelastic = fluid_poro_poroelastic_ispec(inum)
      iedge_poroelastic = fluid_poro_poroelastic_iedge(inum)

      ! implement 1D coupling along the edge
      do ipoin1D = 1,NGLLX

        ! get point values for the poroelastic side, which matches our side in the inverse direction
        i = ivalue_inverse(ipoin1D,iedge_poroelastic)
        j = jvalue_inverse(ipoin1D,iedge_poroelastic)
        iglob = ibool(i,j,ispec_poroelastic)

        ! get point values for the acoustic side
        i = ivalue(ipoin1D,iedge_acoustic)
        j = jvalue(ipoin1D,iedge_acoustic)
        iglob2 = ibool(i,j,ispec_acoustic)

        ! if distance between the two points is not negligible, there is an error, since it should be zero
        if (sqrt((coord(1,iglob) - coord(1,iglob2))**2 + (coord(2,iglob) - coord(2,iglob2))**2) > TINYVAL) &
            call exit_MPI(myrank,'Error in fluid/solid (poroelastic) coupling buffer')

      enddo

    enddo

    ! user output
    if (myrank == 0) then
      write(IMAIN,*) '  End of fluid/solid (poroelastic) edge detection'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  endif

  ! solid/porous edge detection
  ! the two elements forming an edge are already known (computed in meshfem2D),
  ! the common nodes forming the edge are computed here
  if (coupled_elastic_poro) then

    ! checks
    if (ATTENUATION_VISCOELASTIC .or. ATTENUATION_PORO_FLUID_PART) &
      call stop_the_code('Attenuation not supported for mixed elastic/poroelastic simulations')

    if (time_stepping_scheme == 2 .or. time_stepping_scheme == 3) &
      call stop_the_code('RK and LDDRK time scheme not supported for mixed elastic/poroelastic simulations')

    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) 'Mixed elastic/poroelastic simulation'
      write(IMAIN,*)
      write(IMAIN,*) '  Beginning of solid/porous edge detection'
      call flush_IMAIN()
    endif

    ! define the edges of a given element
    i_begin(IBOTTOM) = 1
    j_begin(IBOTTOM) = 1
    i_end(IBOTTOM) = NGLLX
    j_end(IBOTTOM) = 1

    i_begin(IRIGHT) = NGLLX
    j_begin(IRIGHT) = 1
    i_end(IRIGHT) = NGLLX
    j_end(IRIGHT) = NGLLZ

    i_begin(ITOP) = NGLLX
    j_begin(ITOP) = NGLLZ
    i_end(ITOP) = 1
    j_end(ITOP) = NGLLZ

    i_begin(ILEFT) = 1
    j_begin(ILEFT) = NGLLZ
    i_end(ILEFT) = 1
    j_end(ILEFT) = 1

    ! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = ipoin1D
      ivalue_inverse(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,IBOTTOM) = 1
      jvalue_inverse(ipoin1D,IBOTTOM) = 1

      ivalue(ipoin1D,IRIGHT) = NGLLX
      ivalue_inverse(ipoin1D,IRIGHT) = NGLLX
      jvalue(ipoin1D,IRIGHT) = ipoin1D
      jvalue_inverse(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1

      ivalue(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,ITOP) = ipoin1D
      jvalue(ipoin1D,ITOP) = NGLLZ
      jvalue_inverse(ipoin1D,ITOP) = NGLLZ

      ivalue(ipoin1D,ILEFT) = 1
      ivalue_inverse(ipoin1D,ILEFT) = 1
      jvalue(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,ILEFT) = ipoin1D

    enddo


    do inum = 1, num_solid_poro_edges
       ispec_elastic = solid_poro_elastic_ispec(inum)
       ispec_poroelastic = solid_poro_poroelastic_ispec(inum)

        ! one element must be elastic and the other must be poroelastic
        if (ispec_elastic /= ispec_poroelastic .and. ispec_is_elastic(ispec_elastic) .and. &
                 ispec_is_poroelastic(ispec_poroelastic)) then

          ! loop on the four edges of the two elements
          do iedge_poroelastic = 1,NEDGES
            do iedge_elastic = 1,NEDGES

              ! store the matching topology if the two edges match in inverse order
              if (ibool(i_begin(iedge_poroelastic),j_begin(iedge_poroelastic),ispec_poroelastic) == &
                   ibool(i_end(iedge_elastic),j_end(iedge_elastic),ispec_elastic) .and. &
                   ibool(i_end(iedge_poroelastic),j_end(iedge_poroelastic),ispec_poroelastic) == &
                   ibool(i_begin(iedge_elastic),j_begin(iedge_elastic),ispec_elastic)) then
                 solid_poro_elastic_iedge(inum) = iedge_elastic
                 solid_poro_poroelastic_iedge(inum) = iedge_poroelastic
                endif

             enddo
          enddo

       endif

    enddo

    ! make sure solid/porous matching has been perfectly detected: check that the grid points
    ! have the same physical coordinates
    ! loop on all the coupling edges

    if (myrank == 0) then
      write(IMAIN,*) '  Checking solid/porous edge topology...'
      call flush_IMAIN()
    endif

    do inum = 1,num_solid_poro_edges

      ! get the edge of the elastic element
      ispec_elastic = solid_poro_elastic_ispec(inum)
      iedge_elastic = solid_poro_elastic_iedge(inum)

      ! get the corresponding edge of the poroelastic element
      ispec_poroelastic = solid_poro_poroelastic_ispec(inum)
      iedge_poroelastic = solid_poro_poroelastic_iedge(inum)

      ! implement 1D coupling along the edge
      do ipoin1D = 1,NGLLX

        ! get point values for the poroelastic side, which matches our side in the inverse direction
        i = ivalue_inverse(ipoin1D,iedge_elastic)
        j = jvalue_inverse(ipoin1D,iedge_elastic)
        iglob = ibool(i,j,ispec_elastic)

        ! get point values for the elastic side
        i = ivalue(ipoin1D,iedge_poroelastic)
        j = jvalue(ipoin1D,iedge_poroelastic)
        iglob2 = ibool(i,j,ispec_poroelastic)

        ! if distance between the two points is not negligible, there is an error, since it should be zero
        if (sqrt((coord(1,iglob) - coord(1,iglob2))**2 + (coord(2,iglob) - coord(2,iglob2))**2) > TINYVAL) &
            call exit_MPI(myrank,'Error in solid/porous coupling buffer')

      enddo

    enddo

    if (myrank == 0) then
      write(IMAIN,*) '  End of solid/porous edge detection'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  endif

  ! synchronizes all processes
  call synchronize_all()

  ! exclude common points between acoustic absorbing edges and acoustic/elastic matching interfaces
  if (coupled_acoustic_elastic .and. anyabs) then

    if (myrank == 0) then
      write(IMAIN,*) 'excluding common points between acoustic absorbing edges '// &
                     'and acoustic/elastic matching interfaces, if any'
      call flush_IMAIN()
    endif
    ! excludes common points in acoustic domain
    call get_coupling_edges_exclude_common_points(nelemabs,numabs,num_fluid_solid_edges, &
                                                  fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                                                  ibegin_edge1,ibegin_edge2,ibegin_edge3,ibegin_edge4, &
                                                  iend_edge1,iend_edge2,iend_edge3,iend_edge4)
  endif

  ! exclude common points between acoustic absorbing edges and acoustic/poroelastic matching interfaces
  if (coupled_acoustic_poro .and. anyabs) then

    if (myrank == 0) then
      write(IMAIN,*) 'excluding common points between acoustic absorbing edges'// &
                     ' and acoustic/poroelastic matching interfaces, if any'
      call flush_IMAIN()
    endif
    ! excludes common points in acoustic domain
    call get_coupling_edges_exclude_common_points(nelemabs,numabs,num_fluid_poro_edges, &
                                                  fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge, &
                                                  ibegin_edge1,ibegin_edge2,ibegin_edge3,ibegin_edge4, &
                                                  iend_edge1,iend_edge2,iend_edge3,iend_edge4)
  endif

  ! exclude common points between poroelastic absorbing edges and elastic/poroelastic matching interfaces
  if (coupled_elastic_poro .and. anyabs) then

    if (myrank == 0) then
      write(IMAIN,*) 'excluding common points between poroelastic absorbing edges '// &
                     'and elastic/poroelastic matching interfaces, if any'
      call flush_IMAIN()
    endif
    ! excludes common points in poroelastic domain
    call get_coupling_edges_exclude_common_points(nelemabs,numabs,num_solid_poro_edges, &
                                                  solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge, &
                                                  ibegin_edge1_poro,ibegin_edge2_poro,ibegin_edge3_poro,ibegin_edge4_poro, &
                                                  iend_edge1_poro,iend_edge2_poro,iend_edge3_poro,iend_edge4_poro)
  endif

  end subroutine get_coupling_edges

!
!-------------------------------------------------------------------------------
!

  subroutine get_coupling_edges_exclude_common_points(nelemabs,numabs,num_domain_edges, &
                                                      domain_ispec,domain_iedge, &
                                                      ibegin_edge1,ibegin_edge2,ibegin_edge3,ibegin_edge4, &
                                                      iend_edge1,iend_edge2,iend_edge3,iend_edge4)

! excludes common GLL points (in one of the domains) between coupling domains to correct absorbing boundary

  use constants, only: IBOTTOM,ITOP,ILEFT,IRIGHT,NGLLX,NGLLZ

  implicit none

  integer,intent(in) :: nelemabs
  integer,dimension(nelemabs),intent(in) ::numabs

  integer,intent(in) :: num_domain_edges
  integer,dimension(num_domain_edges),intent(in) :: domain_ispec,domain_iedge

  integer,dimension(nelemabs),intent(inout) :: ibegin_edge1,ibegin_edge2,ibegin_edge3,ibegin_edge4
  integer,dimension(nelemabs),intent(inout) :: iend_edge1,iend_edge2,iend_edge3,iend_edge4

  ! local parameters
  integer :: ispec,ispecabs,inum
  integer :: ispec_domain,iedge_domain

  ! loop on all the absorbing elements
  do ispecabs = 1,nelemabs

    ispec = numabs(ispecabs)

    ! loop on all the coupling edges
    do inum = 1,num_domain_edges

      ! get the edge of the acoustic element
      ispec_domain = domain_ispec(inum)
      iedge_domain = domain_iedge(inum)

      ! if acoustic absorbing element and acoustic/poroelastic coupled element is the same
      if (ispec_domain == ispec) then

        if (iedge_domain == IBOTTOM) then
          ibegin_edge4(ispecabs) = 2
          ibegin_edge2(ispecabs) = 2
        endif

        if (iedge_domain == ITOP) then
          iend_edge4(ispecabs) = NGLLZ - 1
          iend_edge2(ispecabs) = NGLLZ - 1
        endif

        if (iedge_domain == ILEFT) then
          ibegin_edge1(ispecabs) = 2
          ibegin_edge3(ispecabs) = 2
        endif

        if (iedge_domain == IRIGHT) then
          iend_edge1(ispecabs) = NGLLX - 1
          iend_edge3(ispecabs) = NGLLX - 1
        endif

      endif

    enddo

  enddo

  end subroutine get_coupling_edges_exclude_common_points



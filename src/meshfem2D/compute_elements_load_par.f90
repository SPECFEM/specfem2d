!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                        Princeton University, USA
!                and CNRS / University of Marseille, France
!                 (there are currently many more authors!)
! (c) Princeton University and CNRS / University of Marseille, July 2012
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================

module compute_elements_load_par

  use part_unstruct_par, only: myrank,nelmnts,nb_edges,abs_surface,nelemabs,glob2loc_elmnts,abs_surface_merge,nelemabs_merge, &
                               nxread,nzread
  use shared_parameters, only: ATTENUATION_VISCOELASTIC,ATTENUATION_VISCOACOUSTIC,NELEM_PML_THICKNESS,PML_BOUNDARY_CONDITIONS, &
                               read_external_mesh,num_material,phi_read,QKappa,Qmu,absorbbottom,absorbtop,absorbleft,absorbright
  use constants, only: IMAIN,TINYVAL

  integer, dimension(:), allocatable  :: elmnts_load
  integer, dimension(:), allocatable  :: adjwgt

  logical, dimension(:), allocatable :: is_elastic,is_acoustic,is_viscoelastic,is_pml

  !! acoustic-elastic-poroelastic as well as CPML load balancing:
  !! we define here the relative cost of all types of spectral elements used in the code.
  integer, parameter :: ACOUSTIC_LOAD = 46
  integer, parameter :: ELASTIC_LOAD = 113
  integer, parameter :: VISCOACOUSTIC_LOAD = 99999 !! not implemented yet, but should be implemented for sure
  integer, parameter :: VISCOELASTIC_LOAD = 280

  integer, parameter :: ACOUSTIC_LOAD_PML = 790
  integer, parameter :: ELASTIC_LOAD_PML = 1049
  integer, parameter :: VISCOACOUSTIC_LOAD_PML = 99999 !! not implemented yet, but should be implemented for sure
  integer, parameter :: VISCOELASTIC_LOAD_PML = 1306

contains

!
!---------------------------------------------------------------------------------------
!

  subroutine compute_elements_load()
    !----------------------------------------------------------------------------------------------
    ! compute elements load for efficient partitioning
    ! (a PML element take more time to calculate than a fluid element for ex)
    ! Fill load_elmnts array contains the element weight to be considered in mesh decomposition
    ! in order to obtain a good load inbalance
    !----------------------------------------------------------------------------------------------

    implicit none

    ! local parameters
    integer :: ier,ispec,nelem_elastic,nelem_acoustic,nelem_viscoelastic,nelem_elastic_pml,nelem_acoustic_pml,nelem_viscoelastic_pml

    ! Allocations WARNING indices start at zero!
    allocate(elmnts_load(0:nelmnts-1), &
             adjwgt(0:nb_edges-1),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating arrays for weights')
    allocate(is_acoustic(0:nelmnts-1), &
             is_elastic(0:nelmnts-1), &
             is_viscoelastic(0:nelmnts-1),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating arrays for is_acoustic, is_elastic and is_viscoelastic')

    ! Initialize counters and arrays
    nelem_elastic = 0
    nelem_acoustic = 0
    nelem_viscoelastic = 0
    nelem_elastic_pml = 0
    nelem_acoustic_pml = 0
    nelem_viscoelastic_pml = 0
    elmnts_load(:) = ELASTIC_LOAD
    adjwgt(:) = ELASTIC_LOAD
    is_acoustic(:) = .false.
    is_elastic(:) = .false.
    is_viscoelastic(:) = .false.

    ! Loop on the elements to determine which element is elastic, acoustic or viscoelastic
    call determine_elements_type()

    ! Determine which element is PML and which is not when the mesh has been made by internal mesher.
    ! For external meshes the detection is made directly when reading the absorbing elements file in
    ! read_external_mesh_file.F90
    if (PML_BOUNDARY_CONDITIONS .and. (.not. read_external_mesh)) then
      call locate_pml_elements_internal()
    endif

    ! Loop on the elements to fill the array elmnts_load and count the elements (for user output)
    do ispec = 0,nelmnts-1
      if (is_elastic(ispec)) then
        if (is_pml(ispec)) then
          elmnts_load(ispec) = ELASTIC_LOAD_PML
          nelem_elastic_pml = nelem_elastic_pml + 1
        else
          nelem_elastic = nelem_elastic + 1
          !elmnts_load(ispec) = ELASTIC_LOAD ! Useless: it has been initialized as elastic
        endif
      else if (is_viscoelastic(ispec)) then
        if (is_pml(ispec)) then
          elmnts_load(ispec) = VISCOELASTIC_LOAD_PML
          nelem_viscoelastic_pml = nelem_viscoelastic_pml + 1
        else
          elmnts_load(ispec) = VISCOELASTIC_LOAD
          nelem_viscoelastic = nelem_viscoelastic + 1
        endif
      else if (is_acoustic(ispec)) then
        if (is_pml(ispec)) then
          elmnts_load(ispec) = ACOUSTIC_LOAD_PML
          nelem_acoustic_pml = nelem_acoustic_pml + 1
        else
          elmnts_load(ispec) = ACOUSTIC_LOAD
          nelem_acoustic = nelem_acoustic + 1
        endif
      endif
    enddo

    ! User output
    if (myrank == 0) then
      write(IMAIN,*)
      write(IMAIN,*) '************ Computing elements load ************'
      write(IMAIN,*) 'Number of elastic elements :',nelem_elastic
      write(IMAIN,*) 'Number of acoustic elements :',nelem_acoustic
      write(IMAIN,*) 'Number of viscoelastic elements :',nelem_viscoelastic
      write(IMAIN,*) 'Number of elastic PML elements :',nelem_elastic_pml
      write(IMAIN,*) 'Number of acoustic PML elements :',nelem_acoustic_pml
      write(IMAIN,*) 'Number of viscoelastic PML elements :',nelem_viscoelastic_pml
      write(IMAIN,*) '*************************************************'
      write(IMAIN,*)
      call flush_IMAIN()
    endif

  end subroutine compute_elements_load

!
!---------------------------------------------------------------------------------------
!

  subroutine determine_elements_type()
    ! Loop on the elements to determine which element is elastic, acoustic or viscoelastic

    implicit none

    integer :: ispec

    do ispec = 0,nelmnts-1
      if (phi_read(num_material(ispec+1)) < TINYVAL) then
        is_acoustic(ispec) = .false.
        is_elastic(ispec) = .true.
      else if (phi_read(num_material(ispec+1)) >= 1.d0) then
        is_acoustic(ispec) = .true.
        is_elastic(ispec) = .false.
      else
        is_acoustic(ispec) = .false.
        is_elastic(ispec) = .false.
      endif

      if (ATTENUATION_VISCOELASTIC) then
        if (((abs(Qkappa(num_material(ispec+1)) - 9999.0d0)) > TINYVAL) .or. &
            ((abs(Qmu(num_material(ispec+1)) - 9999.0d0)) > TINYVAL)) then
          is_viscoelastic(ispec) = .true.
          is_elastic(ispec) = .false.
          is_acoustic(ispec) = .false.
        endif
      endif

      if (ATTENUATION_VISCOACOUSTIC) then
        print *,'warning: there should be some code for ATTENUATION_VISCOACOUSTIC in subroutine determine_elements_type!'
      endif

    enddo

  end subroutine determine_elements_type

!
!---------------------------------------------------------------------------------------
!

  subroutine locate_pml_elements_internal()
    !----------------------------------------------------------------------------------------------
    ! Determine which element is PML and which is not when the mesh has been made by the internal
    ! mesher.
    ! The elements are numbered like this :
    !
    !                              nxread
    ! <----------------------------------------------------->
    !      ______________________________________________________
    !     |          |          |    |           | nxread*nzread | ^
    !     |__________|__________|____|___________|_______________| |
    !     |   ...    |    ...   |....|    ...    |      ...      | |
    !     |__________|__________|____|___________|_______________| |
    !     | nxread+1 | nxread+2 |....|2*nxread-1 |   2*nxread    | | nzread
    !     |__________|__________|____|___________|_______________| |
    !     |     1    |    2     |....| nxread-1  |    nxread     | |
    !     |__________|__________|____|___________|_______________| v
    !
    !     !WARNING! is_pml array starts from 0
    !
    !----------------------------------------------------------------------------------------------

    implicit none

    ! local parameters
    integer :: i,j

      ! PML Left
      if (absorbleft) then
        do i = 1,NELEM_PML_THICKNESS
          do j = 1,nzread
            is_pml((j-1)*nxread + (i-1)) = .true.
          enddo
        enddo
      endif
      ! PML Right
      if (absorbright) then
        do i = nxread-NELEM_PML_THICKNESS+1,nxread
          do j = 1,nzread
            is_pml((j-1)*nxread + (i-1)) = .true.
          enddo
        enddo
      endif
      ! PML Down
      if (absorbbottom) then
        do i = 1,nxread
          do j = 1,NELEM_PML_THICKNESS
            is_pml((j-1)*nxread + (i-1)) = .true.
          enddo
        enddo
      endif
      ! PML Up
      if (absorbtop) then
        do i = 1,nxread
          do j = nzread-NELEM_PML_THICKNESS+1,nzread
            is_pml((j-1)*nxread + (i-1)) = .true.
          enddo
        enddo
      endif

  end subroutine locate_pml_elements_internal

end module compute_elements_load_par


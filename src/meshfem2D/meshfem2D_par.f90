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


  module source_file_par

  use constants, only: MAX_STRING_LEN

  implicit none

  ! source type parameters
  integer, dimension(:),allocatable ::  source_type,time_function_type
  ! location
  double precision, dimension(:),allocatable :: xs,zs
  ! moment tensor
  double precision, dimension(:),allocatable :: Mxx,Mzz,Mxz
  ! source parameters
  double precision, dimension(:),allocatable :: f0_source,tshift_src,anglesource,factor,burst_band_width
  ! flag for fixation to surface
  logical, dimension(:),allocatable ::  source_surf
  ! File name can't exceed MAX_STRING_LEN characters
  character(len=MAX_STRING_LEN), dimension(:),allocatable :: name_of_source_file

  end module source_file_par

!
!---------------------------------------------------------------------------------------
!

  module decompose_par

  implicit none

  ! variables used for storing info about the mesh and partitions
  integer, dimension(:), allocatable  :: my_interfaces
  integer, dimension(:), allocatable  :: my_nb_interfaces

  end module decompose_par

!
!---------------------------------------------------------------------------------------
!

  module part_unstruct_par

! This module contains subroutines related to unstructured meshes and partitioning of the
! corresponding graphs.

  use shared_parameters, only: nelmnts,nxread,nzread,max_npoints_interface,number_of_interfaces, &
    nz_layer,number_of_layers,nx,nz,myrank

  implicit none

  integer, dimension(:), allocatable  :: elmnts
  integer, dimension(:), allocatable  :: elmnts_bis
  integer, dimension(:), allocatable  :: glob2loc_elmnts
  integer, dimension(:), allocatable  :: part

  integer :: nb_edges

  integer, dimension(:), allocatable  :: xadj_g
  integer, dimension(:), allocatable  :: adjncy_g

  integer :: nnodes
  double precision, dimension(:,:), allocatable  :: nodes_coords
  integer, dimension(:), allocatable  :: nnodes_elmnts
  integer, dimension(:), allocatable  :: nodes_elmnts
  integer, dimension(:), allocatable  :: glob2loc_nodes_nparts
  integer, dimension(:), allocatable  :: glob2loc_nodes_parts
  integer, dimension(:), allocatable  :: glob2loc_nodes

  ! interface data
  integer :: ninterfaces
  integer, dimension(:), allocatable  :: tab_size_interfaces, tab_interfaces

  integer :: nelem_acoustic_surface
  integer, dimension(:,:), allocatable  :: acoustic_surface
  integer :: nelem_acoustic_surface_loc

  integer :: nelem_on_the_axis
  integer, dimension(:), allocatable  :: ispec_of_axial_elements
  integer, dimension(:), allocatable  :: inode1_axial_elements, inode2_axial_elements
  integer :: nelem_on_the_axis_loc

  integer :: nelemabs
  integer, dimension(:,:), allocatable  :: abs_surface
  logical, dimension(:,:), allocatable  :: abs_surface_char
  integer, dimension(:), allocatable  :: abs_surface_merge,abs_surface_type
  integer :: nelemabs_loc

  integer :: nelemabs_merge
  integer, dimension(:), allocatable  :: ibegin_edge1,iend_edge1,ibegin_edge3,iend_edge3, &
       ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2

  ! for acoustic/elastic coupled elements
  integer :: nedges_coupled
  integer, dimension(:,:), allocatable  :: edges_coupled

  ! for acoustic/poroelastic coupled elements
  integer :: nedges_acporo_coupled
  integer, dimension(:,:), allocatable  :: edges_acporo_coupled

  ! for poroelastic/elastic coupled elements
  integer :: nedges_elporo_coupled
  integer, dimension(:,:), allocatable  :: edges_elporo_coupled

  ! for acoustic forcing elements
  integer :: nelemacforcing
  integer, dimension(:,:), allocatable :: acforcing_surface
  logical, dimension(:,:), allocatable  :: acforcing_surface_char
  integer, dimension(:), allocatable  :: acforcing_surface_merge,acforcing_surface_type
  integer :: nelemacforcing_loc

  integer :: nelemacforcing_merge
  integer, dimension(:), allocatable  :: ibegin_edge1_acforcing,iend_edge1_acforcing, &
       ibegin_edge3_acforcing,iend_edge3_acforcing,ibegin_edge4_acforcing,iend_edge4_acforcing, &
       ibegin_edge2_acforcing,iend_edge2_acforcing

  ! variables used for tangential detection
  integer ::  nnodes_tangential_curve
  double precision, dimension(:,:), allocatable  :: nodes_tangential_curve

  ! coordinates of the grid points of the mesh
  double precision, dimension(:,:), allocatable :: grid_point_x,grid_point_z

  integer :: npoints_interface_top
  double precision, dimension(:), allocatable :: xinterface_top,zinterface_top,coefs_interface_top

  ! local coupled edges
  integer :: nedges_coupled_loc
  integer :: nedges_acporo_coupled_loc
  integer :: nedges_elporo_coupled_loc

  ! to store the position of PML element in array region_pml_external_mesh
  ! this is only useful when using PML together with external mesh
  integer, dimension(:), allocatable :: region_pml_external_mesh

  integer :: remove_min_to_start_at_zero

  integer :: nspec
  integer :: npgeo
  integer :: iproc

  end module part_unstruct_par

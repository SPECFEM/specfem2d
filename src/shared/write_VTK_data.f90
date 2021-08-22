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


!------------------------------------------------------------------------------------
!
! VTK routines for meshfem3D
!
!------------------------------------------------------------------------------------

! for details:
!
! VTK type definition numbers:
!   https://vtk.org/doc/nightly/html/vtkCellType_8h_source.html
! ordering images see e.g.:
!   https://lorensen.github.io/VTKExamples/site/VTKBook/05Chapter5/


  subroutine write_VTK_data_ngnod_elem_i(nspec,nglob,NGNOD, &
                                         xstore_dummy,zstore_dummy,elmnts, &
                                         elem_flag,filename)

! stores element flags (integers) on 2D quad mesh

  use constants, only: IOUT_VTK,MAX_STRING_LEN

  implicit none

  integer :: nspec,nglob,NGNOD

  ! global coordinates
  integer, dimension(NGNOD,nspec) :: elmnts
  double precision, dimension(nglob) :: xstore_dummy,zstore_dummy

  ! element flag array
  integer, dimension(nspec) :: elem_flag

  ! file name
  character(len=MAX_STRING_LEN) :: filename

  ! local parameters
  integer :: ispec,i,itype

  ! safety check
  if (NGNOD == 4) then
    ! VTK_QUAD == 8 type
    itype = 8
  else if (NGNOD == 9) then
    ! VTK_BIQUADRATIC_QUAD = 28 type
    itype = 28
  else
    stop 'Error invalid NGNOD in write_VTK_data_ngnod_elem_i routine'
  endif

  ! debug
  !print *,'debug:',itype,nglob,nspec,xstore_dummy(1),zstore_dummy(1)
  !print *,'debug:',elmnts(:,1)  ! index values should start at 0

  ! write source and receiver VTK files for Paraview
  open(IOUT_VTK,file=trim(filename),status='unknown')
  write(IOUT_VTK,'(a)') '# vtk DataFile Version 3.1'
  write(IOUT_VTK,'(a)') 'material model VTK file'
  write(IOUT_VTK,'(a)') 'ASCII'
  write(IOUT_VTK,'(a)') 'DATASET UNSTRUCTURED_GRID'
  write(IOUT_VTK, '(a,i12,a)') 'POINTS ', nglob, ' float'
  do i = 1,nglob
    write(IOUT_VTK,'(3e18.6)') xstore_dummy(i),0.0,zstore_dummy(i)  ! assuming Y-coordinate at zero
  enddo
  write(IOUT_VTK,*) ""

  ! note: indices for vtk start at 0
  write(IOUT_VTK,'(a,i12,i12)') "CELLS ",nspec,nspec*(NGNOD+1)
  do ispec = 1,nspec
    if (NGNOD == 4) then
      ! quad4 element using an indexing (left,bottom),(right,bottom),(right,top),(left,top)
      write(IOUT_VTK,'(9i12)') NGNOD,elmnts(1,ispec)-1,elmnts(2,ispec)-1,elmnts(4,ispec)-1,elmnts(3,ispec)-1
    else
      ! quad9 element
      write(IOUT_VTK,'(9i12)') NGNOD,(elmnts(i,ispec)-1,i=1,NGNOD)
    endif
  enddo
  write(IOUT_VTK,*) ""

  ! type: hexahedrons
  write(IOUT_VTK,'(a,i12)') "CELL_TYPES ",nspec
  write(IOUT_VTK,'(6i12)') (itype,ispec=1,nspec)
  write(IOUT_VTK,*) ""

  write(IOUT_VTK,'(a,i12)') "CELL_DATA ",nspec
  write(IOUT_VTK,'(a)') "SCALARS elem_flag integer"
  write(IOUT_VTK,'(a)') "LOOKUP_TABLE default"
  do ispec = 1,nspec
    write(IOUT_VTK,*) elem_flag(ispec)
  enddo
  write(IOUT_VTK,*) ""
  close(IOUT_VTK)

  end subroutine write_VTK_data_ngnod_elem_i



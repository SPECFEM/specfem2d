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


subroutine read_model_nspec()

! reads in nspec from database

  use tomography_par,only: MAX_STRING_LEN,IIN,NSPEC,myrank

  implicit none

  ! local parameters
  integer :: ier
  character(len=MAX_STRING_LEN) :: m_file
  character(len=80) :: datlin

  ! opens database file
  write(m_file,'(a,i5.5)') './OUTPUT_FILES/'//'Database',myrank
  open(IIN,file=trim(m_file),status='old',action='read',iostat=ier)
  if (ier /= 0) then
    print *,'Error opening: ',trim(m_file)
    call exit_MPI(myrank,'file not found')
  endif

  ! reads lines until nspec
  read(IIN,"(a80)") datlin  ! header
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin  ! simulation_title
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin ! AXISYM
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin  ! SIMULATION_TYPE, ...
  read(IIN,"(a80)") datlin
  read(IIN,*) nspec

  close(IIN)

  if (myrank == 0) then
    print *,'number of spectral-elements: ',nspec
    print *,''
  endif

end subroutine read_model_nspec

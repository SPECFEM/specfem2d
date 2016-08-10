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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
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

  subroutine write_output_SU(x_source,z_source,irec,buffer_binary,number_of_components)

  use specfem_par, only: NSTEP,nrec,deltat,seismotype,st_xval, &
                          NSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current,P_SV, &
                          st_zval,subsamp_seismos

  implicit none

  double precision,intent(in) :: x_source,z_source
  integer,intent(in) :: irec,number_of_components

  ! to write seismograms in single precision SEP and double precision binary
  double precision, dimension(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrec,number_of_components),intent(in) :: buffer_binary

  ! local parameters
  integer :: deltat_int2
  integer :: isample

  ! scaling factor for Seismic Unix xsu dislay
  double precision, parameter :: FACTORXSU = 1.d0

  integer(kind=2) :: header2(2)

  if (seismo_offset == 0) then

     if (deltat*1.0d6 > 2**15) then
        deltat_int2 = 0
     else
        deltat_int2 = NINT(deltat*1.0d6, kind=2) ! deltat (unit: 10^{-6} second)
     endif

     ! write SU headers (refer to Seismic Unix for details)
     write(12,rec=(irec-1)*60+(irec-1)*NSTEP+1)  irec                          ! receiver ID
     write(12,rec=(irec-1)*60+(irec-1)*NSTEP+10) NINT(st_xval(irec)-x_source)  ! offset
     write(12,rec=(irec-1)*60+(irec-1)*NSTEP+19) NINT(x_source)                ! source location xs
     write(12,rec=(irec-1)*60+(irec-1)*NSTEP+20) NINT(z_source)                ! source location zs
     write(12,rec=(irec-1)*60+(irec-1)*NSTEP+21) NINT(st_xval(irec))           ! receiver location xr
     write(12,rec=(irec-1)*60+(irec-1)*NSTEP+22) NINT(st_zval(irec))           ! receiver location zr
     if (nrec > 1) write(12,rec=(irec-1)*60+(irec-1)*NSTEP+48) SNGL(st_xval(2)-st_xval(1)) ! receiver interval
     header2(1)=0  ! dummy
     header2(2)=int(NSTEP, kind=2)
     write(12,rec=(irec-1)*60+(irec-1)*NSTEP+29) header2
     header2(1)=deltat_int2
     header2(2)=0  ! dummy
     write(12,rec=(irec-1)*60+(irec-1)*NSTEP+30) header2
     if (seismotype /= 4 .and. seismotype /= 6 .and. P_SV) then
        ! headers
        if (seismo_offset == 0) then
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+1)  irec
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+10) NINT(st_xval(irec)-x_source)
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+19) NINT(x_source)
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+20) NINT(z_source)
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+21) NINT(st_xval(irec))
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+22) NINT(st_zval(irec))
           if (nrec > 1) write(14,rec=(irec-1)*60+(irec-1)*NSTEP+48) SNGL(st_xval(2)-st_xval(1))
           header2(1)=0  ! dummy
           header2(2)=int(NSTEP, kind=2)
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+29) header2
           header2(1)=deltat_int2
           header2(2)=0  ! dummy
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+30) header2
        endif
     endif

  endif

  ! the "60" in the following corresponds to 240 bytes header (note the reclength is 4 bytes)
  do isample = 1, seismo_current
     write(12,rec=irec*60+(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,irec,1))
     if (seismotype /= 4 .and. seismotype /= 6 .and. P_SV) then
        write(14,rec=irec*60+(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,irec,2))
     endif
  enddo

  end subroutine write_output_SU

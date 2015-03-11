
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

  subroutine write_output_SU(x_source,z_source,irec,buffer_binary,number_of_components)

  use specfem_par, only : NSTEP,nrec,deltat,seismotype,st_xval, &
                          NSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current,p_sv, &
                          st_zval,subsamp_seismos
  integer :: deltat_int2
  integer :: irec,isample,number_of_components

! to write seismograms in single precision SEP and double precision binary
! format
  double precision, dimension(NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrec,number_of_components) :: buffer_binary

! scaling factor for Seismic Unix xsu dislay
  double precision, parameter :: FACTORXSU = 1.d0


  double precision :: x_source,z_source
  integer(kind=2) :: header2(2)

  if (seismo_offset==0) then

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
     if (nrec>1) write(12,rec=(irec-1)*60+(irec-1)*NSTEP+48) SNGL(st_xval(2)-st_xval(1)) ! receiver interval
     header2(1)=0  ! dummy
     header2(2)=int(NSTEP, kind=2)
     write(12,rec=(irec-1)*60+(irec-1)*NSTEP+29) header2
     header2(1)=deltat_int2
     header2(2)=0  ! dummy
     write(12,rec=(irec-1)*60+(irec-1)*NSTEP+30) header2
     if ( seismotype /= 4 .and. seismotype /= 6 .and. p_sv) then
        ! headers
        if (seismo_offset==0) then
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+1)  irec
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+10) NINT(st_xval(irec)-x_source)
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+19) NINT(x_source)
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+20) NINT(z_source)
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+21) NINT(st_xval(irec))
           write(14,rec=(irec-1)*60+(irec-1)*NSTEP+22) NINT(st_zval(irec))
           if(nrec>1) write(14,rec=(irec-1)*60+(irec-1)*NSTEP+48) SNGL(st_xval(2)-st_xval(1))
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
     if ( seismotype /= 4 .and. seismotype /= 6 .and. p_sv) then
        write(14,rec=irec*60+(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(isample,irec,2))
     endif
  enddo

  end subroutine write_output_SU

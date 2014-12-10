  subroutine write_output_SU(x_source,z_source,irec,buffer_binary,number_of_components)

  use specfem_par, only : sisux,sisuz,siscurl,station_name,network_name, &
                          NSTEP,nrecloc,which_proc_receiver,nrec,myrank,deltat,seismotype,st_xval,t0, &
                          NSTEP_BETWEEN_OUTPUT_SEISMOS,seismo_offset,seismo_current,p_sv, &
                          st_zval,SU_FORMAT,save_ASCII_seismograms, &
                          save_binary_seismograms_double,subsamp_seismos
  integer :: deltat_int2
  integer :: irec,isample,number_of_components

! to write seismograms in single precision SEP and double precision binary
! format
  double precision, dimension(number_of_components,NSTEP_BETWEEN_OUTPUT_SEISMOS/subsamp_seismos,nrec) :: buffer_binary

! scaling factor for Seismic Unix xsu dislay
  double precision, parameter :: FACTORXSU = 1.d0


  double precision :: x_source,z_source
  integer(kind=2) :: header2(2)

  print*, 'x_source', x_source
  print*, 'z_source', z_source


          if (seismo_offset==0) then

             if (deltat*1.0d6 > 2**15) then
                deltat_int2 = 0
             else
                deltat_int2 = NINT(deltat*1.0d6, kind=2) ! deltat (unit: 10^{-6} second)
             endif

  print*, 'made it here-1'

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
             write(12,rec=irec*60+(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(1,isample,irec))
             if ( seismotype /= 4 .and. seismotype /= 6 .and. p_sv) then
                write(14,rec=irec*60+(irec-1)*NSTEP+seismo_offset+isample) sngl(buffer_binary(2,isample,irec))
             endif
          enddo

  end subroutine write_output_SU

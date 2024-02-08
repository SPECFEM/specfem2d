
  program compute_mean_seismograms

! compute the mean value of sigma_xz seismograms for LCND (Aix-en-Provence, France)

  implicit none

  integer, parameter :: nrec = 99
  integer, parameter :: NSTEP = 14000

  double precision, dimension(NSTEP,nrec) :: sisuz
  double precision, dimension(NSTEP) :: time

  double precision :: mean_value

  character(len=150) sisname

  integer irec,it

!----

! read seismograms in ASCII format
  do irec = 1,nrec
           print *,'reading seismogram ',irec,' out of ',nrec
           write(sisname,"('OUTPUT_FILES/S00',i2.2,'.AA.SXZ.semd')") irec
           open(unit=11,file=sisname,status='old')
           do it = 1,NSTEP
             read(11,*) time(it),sisuz(it,irec)
           enddo
           close(11)
  enddo

! write the mean seismogram
           print *,'writing the mean seismogram'
           open(unit=11,file='OUTPUT_FILES/S00_mean_value.AA.SXZ.semd',status='unknown')
           do it = 1,NSTEP
             mean_value = 0.d0
             do irec = 1,nrec
               mean_value = mean_value + sisuz(it,irec)
             enddo
             mean_value = mean_value / nrec
             write(11,*) sngl(time(it)),sngl(mean_value)
           enddo
           close(11)

  end program compute_mean_seismograms


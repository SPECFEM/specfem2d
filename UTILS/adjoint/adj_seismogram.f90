      program adj_seismogram

! This program cuts a certain portion of the seismograms and convert it
! into the adjoint source for generating banana-dougnut kernels

      implicit none

      integer, parameter :: NSTEP = 3000
      integer, parameter :: nrec = 1
      double precision, parameter :: t0 = 6d-2
      double precision, parameter :: deltat = 2d-4
      double precision, parameter :: EPS = 1.d-40
      integer :: itime,icomp,istart,iend,nlen,irec
      double precision :: time,tstart(nrec),tend(nrec)
      character(len=150), dimension(nrec) :: station_name
      double precision, dimension(NSTEP) :: time_window
      double precision :: seism(NSTEP,2),Nnorm,seism_win(NSTEP)
      double precision :: seism_veloc(NSTEP),seism_accel(NSTEP),ft_bar(NSTEP)
      character(len=3) :: comp(2)
      character(len=150) :: filename,filename2

      include 'constants.h'

      station_name(1) = 'S0001'
      tstart(1) = 0.031d0 + t0
      tend(1) = 0.121d0 + t0

      comp = (/"BHX","BHZ"/)
     
      do irec =1,nrec

        do icomp = 1, NDIM

      filename = 'OUTPUT_FILES/'//trim(station_name(irec))//'.AA.'// comp(icomp) // '.semd'
      open(unit = 10, file = trim(filename))

         do itime = 1,NSTEP
        read(10,*) time , seism(itime,icomp)
         enddo
 
        enddo

      close(10)


         istart = max(floor(tstart(irec)/deltat),1)
         iend = min(floor(tend(irec)/deltat),NSTEP)
         print*,'istart =',istart, 'iend =', iend
         print*,'tstart =',istart*deltat, 'tend =', iend*deltat
         if(istart >= iend) stop 'check istart,iend'
         nlen = iend - istart +1
         
       do icomp = 1, NDIM

      filename = 'OUTPUT_FILES/'//trim(station_name(irec))//'.AA.'// comp(icomp) // '.adj'
      open(unit = 11, file = trim(filename))

        time_window(:) = 0.d0
        seism_win(:) = seism(:,icomp)
        seism_veloc(:) = 0.d0
        seism_accel(:) = 0.d0

        do itime =istart,iend
!        time_window(itime) = 1.d0 - cos(pi*(itime-1)/NSTEP+1)**10   ! cosine window         
        time_window(itime) = 1.d0 - (2* (dble(itime) - istart)/(iend-istart) -1.d0)**2  ! Welch window         
        enddo

         do itime = 2,NSTEP-1
      seism_veloc(itime) = (seism_win(itime+1) - seism_win(itime-1))/(2*deltat)
         enddo
      seism_veloc(1) = (seism_win(2) - seism_win(1))/deltat
      seism_veloc(NSTEP) = (seism_win(NSTEP) - seism_win(NSTEP-1))/deltat
 
         do itime = 2,NSTEP-1
      seism_accel(itime) = (seism_veloc(itime+1) - seism_veloc(itime-1))/(2*deltat)
         enddo
      seism_accel(1) = (seism_veloc(2) - seism_veloc(1))/deltat
      seism_accel(NSTEP) = (seism_veloc(NSTEP) - seism_veloc(NSTEP-1))/deltat

      Nnorm = deltat * sum(time_window(:) * seism_win(:) * seism_accel(:))
!      Nnorm = deltat * sum(time_window(:) * seism_veloc(:) * seism_veloc(:))
! cross-correlation traveltime adjoint source
      if(abs(Nnorm) > EPS) then
!      ft_bar(:) = - seism_veloc(:) * time_window(:) / Nnorm
      ft_bar(:) = seism_veloc(:) * time_window(:) / Nnorm
      print*,'Norm =', Nnorm
      else
      print *, 'norm < EPS for file '
      print*,'Norm =', Nnorm
      ft_bar(:) = 0.d0
      endif

       do itime =1,NSTEP
        if(icomp == 1) then
      write(11,*) (itime-1)*deltat - t0, ft_bar(itime)
!      write(12,*) (itime-1)*deltat - t0, seism_veloc(itime)
        else
      write(11,*) (itime-1)*deltat - t0, 0.d0
        endif
       enddo

        enddo
      close(11)
!      close(12)

      enddo


      end program adj_seismogram

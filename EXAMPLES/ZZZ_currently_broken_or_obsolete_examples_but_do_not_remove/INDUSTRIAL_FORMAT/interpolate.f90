  program MAIN

  !Modified by Qiancheng Liu and Daniel Peter at Kaust

  implicit none
  !include "constants.h"
  character(len=512),parameter :: sep_directory='./SEG_2D_SALT/'
  character(len=512),parameter :: sep_header_file_vp='vp.H'
  character(len=512),parameter :: sep_header_file_vs='vs.H'   ! not used
  character(len=512),parameter :: sep_header_file_rho='rho.H' ! not used
  integer :: NX,NY,NZ,esize,nproc,scalar
  real :: OX,OY,OZ,DX,DY,DZ
  real,dimension(:,:),allocatable :: vp_SEP,vs_SEP,rho_SEP
  character(len=512) :: system_command,sep_file,data_format,mesh_file,wavespeed_file
  character(len=40) :: tmpstring
  character(len=20) :: numstring
  character(len=3) ::num
  logical closest_point  

  closest_point=.true. !Two options: closest_point ('closet_point=.true.') or bilinear interpolation ('closet_point=.false.')
  scalar=1
     ! creating mesh/grid
  print *, 'please input the number of processes'
  read(*,*) nproc

  call system('sed -i "s/SAVE_FORWARD                    = .true./SAVE_FORWARD                    = .false./g" ./DATA/Par_file')
  call system('sed -i "/NSTEP   /c NSTEP                              = 5          # total number of time steps" ./DATA/Par_file')
  call system('sed -i "/^MODEL   /c MODEL                              = default" ./DATA/Par_file')
  call system('sed -i "/^SAVE_MODEL   /c SAVE_MODEL                              = legacy" ./DATA/Par_file')

     ! reading SEP models
  print *, 'reading SEP models...'
     !!! VP model in SEP format
  call READ_SEP_HEADER(sep_directory,sep_header_file_vp, &
                          sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ, &
                          esize,data_format)
  allocate(vp_SEP(NX,NZ),vs_SEP(NX,NZ),rho_SEP(NX,NZ)) ! assuming that VP,VS,RHO models have the same dimension
  call READ_SEP_MODEL_2D(sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ,esize,data_format,vp_SEP)

  !call READ_SEP_HEADER(sep_directory,sep_header_file_vs, &
  !                     sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ, &
  !                     esize,data_format)
  !call READ_SEP_MODEL_2D(sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ,esize,data_format,vs_SEP)

  !call READ_SEP_HEADER(sep_directory,sep_header_file_rho, &
  !                     sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ, &
  !                     esize,data_format)
  !call READ_SEP_MODEL_2D(sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ,esize,data_format,rho_SEP)

  print *,'Writting interface data to file DATA/interface_industry.dat'

  open(unit=15,file='./DATA/interface_industry.dat',status='unknown')
  write(15,*)'#'
  write(15,*)'# number of interfaces'     
  write(15,*)'#'
  write(15,*)' 2'     
  write(15,*)'#'
  write(15,*)'# for each interface below, we give the number of points and then x,z for each point'
  write(15,*)'#'
  write(15,*)'#'
  write(15,*)'# interface number 1 (bottom of the mesh)'
  write(15,*)'#'
  write(15,*)' 2'
  write(15,*) OX,OZ
  write(15,*) OX+(NX-1)*DX,OZ
  write(15,*)'#'
  write(15,*)'# interface number 2 (topography, top of the mesh)'
  write(15,*)'#'
  write(15,*)' 2'
  write(15,*) OX,OZ+(NZ-1)*DZ
  write(15,*) OX+(NX-1)*DX,OZ+(NZ-1)*DZ
  write(15,*)'#'
  write(15,*)'# for each layer, we give the number of spectral elements in the vertical direction'
  write(15,*)'#'
  write(15,*)'#'
  write(15,*)'# layer number 1 (bottom layer)'
  write(15,*)'#'
  write(15,*) (NZ-1)/scalar
  close(15)

  print *,'OX=',OX
  write(numstring,*) OX
  write(tmpstring,*) 'xmin			=',numstring
  write(system_command,*) tmpstring
  write(system_command,*) 'sed -i "/xmin/c',tmpstring,'" ./DATA/Par_file'
  write (*,*) system_command
  call system(system_command)    

  write(numstring,*) OX+(NX-1)*DX
  write(tmpstring,*) 'xmax			=',numstring
  write(*,*) tmpstring
  write(system_command,*) 'sed -i "/xmax/c',tmpstring,'" ./DATA/Par_file'
  write (*,*) system_command
  call system(system_command)

  write(numstring,*) (NX-1)/scalar
  write(tmpstring,*) 'nx			 =',numstring
  write(*,*) tmpstring
  write(system_command,*) 'sed -i "/^nx/c',tmpstring,'" ./DATA/Par_file'
  write (*,*) system_command
  call system(system_command)

  call system('sed -i "$ d" ./DATA/Par_file')
  write(tmpstring,*) '1',(NX-1)/scalar,'1',(NZ-1)/scalar,'1'
  write (*,*) tmpstring
  write(system_command,*) 'sed -i "/nbregions   /a',tmpstring,'" ./DATA/Par_file'
     !write(system_command,*) 'echo ${tmpstring}'
  write (*,*) system_command
  call system(system_command)


  write(tmpstring,*) 'NPROC			 =',nproc
  write(*,*) tmpstring
  write(system_command,*) 'sed -i "/^NPROC/c',tmpstring,'" ./DATA/Par_file'
  write (*,*) system_command
  call system(system_command)

  call sleep(1)

  write(num,'(i2.2)') nproc
  print *,'call xmeshfem2d'
  call system('./xmeshfem2D >OUTPUT_FILES/output_mesher.txt') 

  print *,'call xspecfem2d'
  call system('mpirun -np ' //num//' ./xspecfem2D>OUTPUT_FILES/output_solver.txt')

  call system('sed -i "/^MODEL   /c MODEL                              = legacy" ./DATA/Par_file')
  call system('sed -i "/^SAVE_MODEL   /c SAVE_MODEL                              = default" ./DATA/Par_file')
  call system('sed -i "/NSTEP   /c NSTEP                            = 10000     # total number of time steps" ./DATA/Par_file')

  print *,'call xmeshfem2d'
  call system('./xmeshfem2D')

  !!! VS model, scaled or in SEP format
  !vs_SEP=0.0
  vs_SEP=vp_SEP/1.732
     !call READ_SEP_HEADER(sep_directory,sep_header_file_vs, &
     !                     sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ, &
     !                     esize,data_format)
     !call READ_SEP_MODEL_2D(sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ,esize,data_format,vs_SEP)

     !!! RHO model, scaled or in SEP format
  rho_SEP  =  1000.0
     !rho_SEP  =   (1.6612 * (vp_SEP/14500*4500 / 1000.0)     &
     !             -0.4720 * (vp_SEP/14500*4500 / 1000.0)**2  &
     !             +0.0671 * (vp_SEP/14500*4500 / 1000.0)**3  &
     !             -0.0043 * (vp_SEP/14500*4500 / 1000.0)**4  &
     !             +0.000106*(vp_SEP/14500*4500 / 1000.0)**5 )*1000.0
     !call READ_SEP_HEADER(sep_directory,sep_header_file_rho, &
     !                     sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ, &
     !                     esize,data_format)
     !call READ_SEP_MODEL_2D(sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ,esize,data_format,rho_SEP)

     ! interpolating models

  vp_SEP=1.0/vp_SEP
  vs_SEP=1.0/vs_SEP
  rho_SEP=1.0/rho_SEP

  call interpolate(vp_SEP,vs_SEP,rho_SEP,NX,NY,NZ,DX,DY,DZ,OX,OY,OZ,closest_point,nproc)

  print *,'call xspecfem2d'
  call system('mpirun -np ' //num//' ./xspecfem2D ')

  end program MAIN

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

  subroutine READ_SEP_HEADER(sep_directory,sep_header_file, &
                             sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ, &
                             esize,data_format)

  implicit none

  !!! input
  character(len=512),intent(in) :: sep_directory    ! where the sep header file and sep file are saved
  character(len=512),intent(in) :: sep_header_file  ! file name of the sep header file

  !!! output
  character(len=512),intent(out) :: sep_file    ! sep file specified in sep header file
  integer,intent(out) :: NX,NY,NZ               ! n1,n2,n3 in sep header file
  real,intent(out)    :: OX,OY,OZ               ! o1,o2,o3 in sep header file
  real,intent(out)    :: DX,DY,DZ               ! d1,d2,d3 in sep header file
  integer,intent(out) :: esize                  ! esize in sep header file
  character(len=512),intent(out) :: data_format ! data_format in sep header file

  !!! local
  integer :: ier
  character(len=512) :: junk,sep_header_file_complete

  sep_header_file_complete=trim(adjustl(sep_directory))//trim(adjustl(sep_header_file))

  open(unit=13,file=trim(adjustl(sep_header_file_complete)),status='old',iostat=ier)
  print *
  print *, '*******************************************************************************'
  print *, 'reading sep header file: '
  print *, trim(adjustl(sep_header_file_complete))
  if (ier /= 0) stop 'ERROR: cannot open sep header file'

  read(13,'(a3a)')     junk, sep_file
  read(13,'(a3i10)')     junk, NX
  read(13,'(a3i10)')     junk, NY
  read(13,'(a3i10)')     junk, NZ
  read(13,'(a3f20.0)') junk, OX
  read(13,'(a3f20.0)') junk, OY
  read(13,'(a3f20.0)') junk, OZ
  read(13,'(a3f20.0)') junk, DX
  read(13,'(a3f20.0)') junk, DY
  read(13,'(a3f20.0)') junk, DZ
  read(13,'(a6i10)')     junk, esize
  read(13,'(a13a)')    junk, data_format
  close(13)
  sep_file=trim(adjustl(sep_directory))//trim(adjustl(sep_file))
  data_format=data_format(1:len_trim(adjustl(data_format))-1)

  print *
  print *, 'sep file specified in the header file is: ', trim(adjustl(sep_file))
  print *, 'NX,NY,NZ = ', NX,NY,NZ
  print *, 'OX,OY,OZ = ', OX,OY,OZ
  print *, 'DX,DY,DZ = ', DX,DY,DZ
  print *, 'esize = ', esize
  print *, 'data_format = ', trim(adjustl(data_format))
  print *, '*******************************************************************************'
  print *

  end subroutine READ_SEP_HEADER

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

  subroutine READ_SEP_MODEL_2D(sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ,esize,data_format, &
                               model)

  implicit none

  !!! input
  character(len=512),intent(in) :: sep_file    ! sep file specified in sep header file
  integer,intent(in) :: NX,NY,NZ               ! n1,n2,n3 in sep header file
  real,intent(in)    :: OX,OY,OZ               ! o1,o2,o3 in sep header file
  real,intent(in)    :: DX,DY,DZ               ! d1,d2,d3 in sep header file
  integer,intent(in) :: esize                  ! esize in sep header file
  character(len=512),intent(in) :: data_format ! data_format in sep header file

  !!! output
  real,dimension(NX,NZ),intent(out) :: model

  !!! local
  integer :: ier

  ! check whether the model is 2D (NY==1)
  if (NY /= 1) stop 'ERROR: this only works for 2D problems (NY/n2 must be 1)'

  ! note that we keep NY as general in the following (for 3D problems in the future)
  open(unit=14,file=trim(adjustl(sep_file)),access='direct',status='old',recl=4*NX*NY*NZ,iostat=ier)
  print *, '*******************************************************************************'
  print *, 'reading sep file: '
  print *, trim(adjustl(sep_file))
  if (ier /= 0) stop 'ERROR: cannot open sep file'

  read(14,rec=1,iostat=ier) model(:,:)
  close(14)
  if (ier /= 0) stop 'ERROR: reading sep file'
  print *, 'done reading sucessfully'
  print *, '*******************************************************************************'
  print *

  end subroutine READ_SEP_MODEL_2D

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

  subroutine interpolate(vp,vs,rho,NX,NY,NZ,DX,DY,DZ,OX,OY,OZ,closest_point,nproc)

  implicit none

  !!! constants
  integer, parameter :: NGLLX=5,NGLLZ=5
  integer :: iproc
  double precision, parameter :: vp_water=5500 ! upper bound
  logical :: closest_point
  !!! input
  integer,intent(in) :: NX,NY,NZ,nproc
  real,intent(in)    :: OX,OY,OZ
  real,intent(in)    :: DX,DY,DZ
  real,dimension(NX,NZ),intent(in) :: vp,vs,rho

  character(len=512) :: input_file
  character(len=512) :: output_file
  character(len=512) :: system_command
  !!! output
  ! currently no output, all information is saved in wavespeed_file

  !!! local
  integer :: ier,ix,iz,ix1,iz1,i
  integer, dimension(NGLLX*NGLLZ) :: iglob
  double precision :: rho_temp,vp_temp,vs_temp,a1,b1,tmp1,tmp2
  double precision, dimension(NGLLX*NGLLZ) :: rho_new,vp_new,vs_new,x,z


  do iproc=1,nproc
    print *, 'iproc=',iproc
    write(input_file,'(a,i6.6,a)') 'DATA/proc',iproc-1,'_model_velocity.dat_input'
    write(output_file,'(a,i6.6,a)') 'DATA/proc',iproc-1,'_model_velocity.dat_tmp'
    open(unit=15,file=trim(adjustl(input_file)),status='old',iostat=ier)
    if (ier /= 0) stop 'Error opening DATA/proc*****_model_velocity.dat_input file.'
    open(unit=16,file=trim(adjustl(output_file)),status='unknown',iostat=ier)
    ier=0
    do while (ier == 0)
      do i=1,NGLLX*NGLLZ
        read(15,'(I10,5e15.5e4)',iostat=ier)  iglob(i),x(i),z(i),rho_temp,vp_temp,vs_temp
        if(closest_point) then       
          ix=NINT((x(i)-OX)/DX)+1
          iz=NINT((z(i)-OZ)/DZ)+1

          if (ix > NX) ix=NX
          if (ix < 1)  ix=1
          if (iz > NZ) iz=NZ
          if (iz < 1)  iz=1

          rho_new(i)=rho(ix,iz)
          vp_new(i)=vp(ix,iz)
          vs_new(i)=vs(ix,iz)
        else
          ix=INT((x(i)-OX)/DX)+1
          iz=INT((z(i)-OZ)/DZ)+1
	  ix1=ix+1
	  iz1=iz+1

          if (ix1 > NX) ix1=NX
          if (iz1 > NZ) iz1=NZ

	  a1=x(i)-OX-(ix-1)*DX
	  a1=a1/DX
	  b1=z(i)-OZ-(iz-1)*DZ
	  b1=b1/DZ

	  tmp1=(1-a1)*vp(ix,iz)+a1*vp(ix1,iz)
	  tmp2=(1-a1)*vp(ix,iz1)+a1*vp(ix1,iz1)
	  tmp1=(1-b1)*tmp1+b1*tmp2
	  vp_new(i)=tmp1

	  tmp1=(1-a1)/vs(ix,iz)+a1/vs(ix1,iz)
	  tmp2=(1-a1)/vs(ix,iz1)+a1/vs(ix1,iz1)
	  tmp1=(1-b1)*tmp1+b1*tmp2
	  vs_new(i)=1/tmp1

	  tmp1=(1-a1)*rho(ix,iz)+a1*rho(ix1,iz)
	  tmp2=(1-a1)*rho(ix,iz1)+a1*rho(ix1,iz1)
	  tmp1=(1-b1)*tmp1+b1*tmp2
      	  rho_new(i)=tmp1
        endif
      enddo
      if (ier == 0) then
        do i=1,NGLLX*NGLLZ
          write(16,'(I10,5e15.5e4)') iglob(i),x(i),z(i),1.0/rho_new(i),1.0/vp_new(i),1.0/vs_new(i)
        enddo
      endif
    enddo

    close(15)
    close(16)

    write(system_command,*) 'mv -f ',trim(adjustl(output_file)),' ',trim(adjustl(input_file))
    write(*,*) system_command
    call system(system_command)

  enddo

  end subroutine interpolate

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

  program MAIN

  implicit none
  !include "constants.h"
  character(len=512),parameter :: sep_directory='./EXAMPLES/INDUSTRIAL_FORMAT/SEG_2D_SALT/'
  character(len=512),parameter :: sep_header_file_vp='vp.H'
  character(len=512),parameter :: sep_header_file_vs='vs.H'   ! not used
  character(len=512),parameter :: sep_header_file_rho='rho.H' ! not used
  integer :: NX,NY,NZ,esize
  real :: OX,OY,OZ,DX,DY,DZ
  real,dimension(:,:),allocatable :: vp_SEP,vs_SEP,rho_SEP
  character(len=512) :: system_command,sep_file,data_format,mesh_file,wavespeed_file

     ! creating mesh/grid
     print *, 'generating mesh/grid...'
     write(system_command,"('sed -e ""s#^SIMULATION_TYPE.*#SIMULATION_TYPE = 1 #g""                     < ./DATA/Par_file > temp; mv temp ./DATA/Par_file')"); call system(system_command)
     write(system_command,"('sed -e ""s#^SAVE_FORWARD.*#SAVE_FORWARD = .false. #g""                     < ./DATA/Par_file > temp; mv temp ./DATA/Par_file')"); call system(system_command)
     write(system_command,"('sed -e ""s#^assign_external_model.*#assign_external_model = .false. #g""   < ./DATA/Par_file > temp; mv temp ./DATA/Par_file')"); call system(system_command)
     write(system_command,"('sed -e ""s#^READ_EXTERNAL_SEP_FILE.*#READ_EXTERNAL_SEP_FILE = .false. #g"" < ./DATA/Par_file > temp; mv temp ./DATA/Par_file')"); call system(system_command)
     write(system_command,"('sed -e ""s#^nt.*#nt = 5 #g""                                               < ./DATA/Par_file > temp; mv temp ./DATA/Par_file')"); call system(system_command)
     write(system_command,"('./xmeshfem2D')"); call system(system_command)
     write(system_command,"('./xspecfem2D')"); call system(system_command)

     ! reading SEP models
     print *, 'reading SEP models...'
     !!! VP model in SEP format
     call READ_SEP_HEADER(sep_directory,sep_header_file_vp, &
                          sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ, &
                          esize,data_format)
     allocate(vp_SEP(NX,NZ),vs_SEP(NX,NZ),rho_SEP(NX,NZ)) ! assuming that VP,VS,RHO models have the same dimension
     call READ_SEP_MODEL_2D(sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ,esize,data_format,vp_SEP)

     !!! VS model, scaled or in SEP format
     vs_SEP=0.0
     !vs_SEP=vp_SEP/1.732
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
     print *, 'interpolating models onto mesh/grid...'
     mesh_file="./DATA/model_velocity.dat_output"
     wavespeed_file="./DATA/model_velocity.dat_input"
     call interpolate(vp_SEP,vs_SEP,rho_SEP,NX,NY,NZ,DX,DY,DZ,OX,OY,OZ,mesh_file,wavespeed_file)
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
  print *, ''
  print *, '*******************************************************************************'
  print *, 'reading sep header file: '
  print *, trim(adjustl(sep_header_file_complete))
  if (ier/=0) stop 'ERROR: cannot open sep header file'

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

  print *, ''
  print *, 'sep file specified in the header file is: ', trim(adjustl(sep_file))
  print *, 'NX,NY,NZ = ', NX,NY,NZ
  print *, 'OX,OY,OZ = ', OX,OY,OZ
  print *, 'DX,DY,DZ = ', DX,DY,DZ
  print *, 'esize = ', esize
  print *, 'data_format = ', trim(adjustl(data_format))
  print *, '*******************************************************************************'
  print *, ''

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
  if (NY/=1) stop 'ERROR: this only works for 2D problems (NY/n2 must be 1)'

  ! note that we keep NY as general in the following (for 3D problems in the future)
  open(unit=14,file=trim(adjustl(sep_file)),access='direct',status='old',recl=4*NX*NY*NZ,iostat=ier)
  print *, '*******************************************************************************'
  print *, 'reading sep file: '
  print *, trim(adjustl(sep_file))
  if (ier/=0) stop 'ERROR: cannot open sep file'

  read(14,rec=1,iostat=ier) model(:,:)
  close(14)
  if (ier/=0) stop 'ERROR: reading sep file'
  print *, 'done reading sucessfully'
  print *, '*******************************************************************************'
  print *, ''

  end subroutine READ_SEP_MODEL_2D

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------

  subroutine interpolate(vp,vs,rho,NX,NY,NZ,DX,DY,DZ,OX,OY,OZ,mesh_file,wavespeed_file)

  implicit none

  !!! constants
  integer, parameter :: NGLLX=5,NGLLZ=5
  double precision, parameter :: vp_water=5500 ! upper bound

  !!! input
  integer,intent(in) :: NX,NY,NZ
  real,intent(in)    :: OX,OY,OZ
  real,intent(in)    :: DX,DY,DZ
  real,dimension(NX,NZ),intent(in) :: vp,vs,rho
  character(len=512),intent(in) :: mesh_file
  character(len=512),intent(in) :: wavespeed_file

  !!! output
  ! currently no output, all information is saved in wavespeed_file

  !!! local
  integer :: ier,ix,iz,i
  integer, dimension(NGLLX*NGLLZ) :: iglob
  double precision :: rho_temp,vp_temp,vs_temp
  double precision, dimension(NGLLX*NGLLZ) :: rho_new,vp_new,vs_new,x,z


  open(unit=15,file=trim(adjustl(mesh_file)),status='old')
  open(unit=16,file=trim(adjustl(wavespeed_file)),status='unknown')
  ier=0
  do while (ier==0)
     do i=1,NGLLX*NGLLZ
        read(15,'(I10, 5F13.4)',iostat=ier)  iglob(i),x(i),z(i),rho_temp,vp_temp,vs_temp
        ix=NINT((x(i)-OX)/DX)+1
        iz=NINT((z(i)-OZ)/DZ)+1
        if (ix>NX-2) ix=NX-2
        if (ix<1)  ix=1
        if (iz>NZ) iz=NZ
        if (iz<1)  iz=1
        rho_new(i)=rho(ix,iz)
        vp_new(i)=vp(ix,iz)
        vs_new(i)=vs(ix,iz)
        if(vs_temp<1000.0) vs_new(i)=0.0
        if(vs_temp<1000.0) rho_new(i)=1000.0
     enddo
     if (ier==0) then
     do i=1,NGLLX*NGLLZ
        write(16,'(I10, 5F13.4)') iglob(i),x(i),z(i),rho_new(i),vp_new(i),vs_new(i)
     enddo
     endif
  enddo

  close(15)
  close(16)

  end subroutine interpolate

!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------
!-----------------------------------------------------------------------------------


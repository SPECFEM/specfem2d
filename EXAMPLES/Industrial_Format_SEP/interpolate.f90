! this program creates first a specfem mesh based on the SEP-format mesh, then interpolates the wave speeds
! onto all GLL points and runs a new simulation with the interpolated mesh.

  program main_interpolate

  implicit none

!-----------------------------------------------------------------------

  ! SEP format model files
  character(len=512),parameter :: sep_directory = './SEG_2D_SALT/'
  character(len=512),parameter :: sep_header_file_vp  = 'vp.H'
  character(len=512),parameter :: sep_header_file_vs  = 'vs.H'  ! not used if scaling from vp
  character(len=512),parameter :: sep_header_file_rho = 'rho.H' ! not used if scaling from vp

  ! Option to scale VS & RHO from VP
  logical,parameter :: SCALE_FROM_VP = .true.

  ! Interpolation: closest_point ('closet_point=.true.') or bilinear interpolation ('closet_point=.false.')
  logical,parameter :: INTERPOLATE_FROM_CLOSEST_POINT = .true.

  ! Option for interpolation points: 1 == corners only / 2 == corner & midpoints used for interpolation
  integer,parameter :: INTERPOLATION_POINTS = 1

  ! lets user choose number of MPI processes
  logical,parameter :: USER_INPUT = .false.

!-----------------------------------------------------------------------

  integer :: NX,NY,NZ,esize,nproc,ier
  real :: OX,OY,OZ,DX,DY,DZ
  real,dimension(:,:),allocatable :: vp_SEP,vs_SEP,rho_SEP

  character(len=512) :: system_command,sep_file,data_format,line
  character(len=128) :: tmpstring
  character(len=128) :: numstring
  character(len=3) ::num

  ! user output
  print *
  print *, '***********************************************************'
  print *, 'Model Interpolation'
  print *, '***********************************************************'
  print *

  ! user input
  if (USER_INPUT) then
    print *, 'please input the number of processes:'
    read(*,*) nproc
    print *
  else
    ! reads from default Par_file
    nproc = 0
    open(unit=15,file='./DATA/Par_file',status='old',iostat=ier)
    if (ier /= 0) stop 'Error opening DATA/Par_file'
    do while (ier == 0)
      read(15,'(a512)',iostat=ier) line
      if (ier == 0) then
        ! left trim
        line = adjustl(line)
        ! skip comment lines
        if (line(1:1) == '#') cycle
        ! suppress trailing comment
        if (index(line,'#') > 0) line = line(1:index(line,'#')-1)
        ! checks parameter name
        !print *,'line: ',line(1:5)
        if (line(1:5) == 'NPROC') then
          read(line(index(line,'=')+1:len_trim(line)),'(i3)') nproc
          exit
        endif
      endif
    enddo
    close(15)
    if (nproc < 1) stop 'Error could not read nproc from Par_file'
  endif

  ! user output
  print *, 'setting:'
  print *, '  SEP model from directory : ',trim(sep_directory)
  print *, '  using number of processes: ',nproc
  if (INTERPOLATE_FROM_CLOSEST_POINT) then
    print *, '  interpolation type       : from closest point'
  else
    print *, '  interpolation type       : bilinear'
  endif
  print *, '  interpolation points set : ',INTERPOLATION_POINTS
  print *, '  using scaling from Vp    : ',SCALE_FROM_VP
  print *

  ! reading SEP models
  print *, 'reading SEP models...'

  ! VP model
  ! gets header info for VP model in SEP format
  call READ_SEP_HEADER(sep_directory,sep_header_file_vp, &
                       sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ, &
                       esize,data_format)

  ! assuming that VP,VS,RHO models have the same dimension
  allocate(vp_SEP(NX,NZ),vs_SEP(NX,NZ),rho_SEP(NX,NZ))

  ! reads in VP model
  call READ_SEP_MODEL_2D(sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ,esize,data_format,vp_SEP)

  ! VS model
  if (.not. SCALE_FROM_VP) then
    call READ_SEP_HEADER(sep_directory,sep_header_file_vs, &
                        sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ, &
                        esize,data_format)
    call READ_SEP_MODEL_2D(sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ,esize,data_format,vs_SEP)
  endif

  ! RHO model
  if (.not. SCALE_FROM_VP) then
    call READ_SEP_HEADER(sep_directory,sep_header_file_rho, &
                         sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ, &
                         esize,data_format)
    call READ_SEP_MODEL_2D(sep_file,NX,NY,NZ,OX,OY,OZ,DX,DY,DZ,esize,data_format,rho_SEP)
  endif

  ! creates interface file
  print *
  print *,'Writing interface data to file DATA/interface_industry.dat'
  print *

  open(unit=15,file='./DATA/interface_industry.dat',status='unknown')
  write(15,'(a)')'#'
  write(15,'(a)')'# number of interfaces'
  write(15,'(a)')'#'
  write(15,'(a)')' 2'
  write(15,'(a)')'#'
  write(15,'(a)')'# for each interface below, we give the number of points and then x,z for each point'
  write(15,'(a)')'#'
  write(15,'(a)')'#'
  write(15,'(a)')'# interface number 1 (bottom of the mesh)'
  write(15,'(a)')'#'
  write(15,'(a)')' 2'
  write(15,*) OX,OZ
  write(15,*) OX+(NX-1)*DX,OZ
  write(15,'(a)')'#'
  write(15,'(a)')'# interface number 2 (topography, top of the mesh)'
  write(15,'(a)')'#'
  write(15,'(a)')' 2'
  write(15,*) OX,OZ+(NZ-1)*DZ
  write(15,*) OX+(NX-1)*DX,OZ+(NZ-1)*DZ
  write(15,'(a)')'#'
  write(15,'(a)')'# for each layer, we give the number of spectral elements in the vertical direction'
  write(15,'(a)')'#'
  write(15,'(a)')'#'
  write(15,'(a)')'# layer number 1 (bottom layer)'
  write(15,'(a)')'#'
  write(15,*) (NZ-1)/INTERPOLATION_POINTS
  close(15)

  ! model geometry
  print *
  print *, '************ Preparing Par_file ************ '
  print *

  !print *, 'OX = ',OX
  write(numstring,*) OX
  write(tmpstring,*) 'xmin                            =',trim(numstring)
  tmpstring = adjustl(tmpstring)
  !write(*,*) trim(tmpstring)
  write(system_command,*) 'sed -i "s/^xmin .*/',trim(tmpstring),'/" ./DATA/Par_file'
  write (*,*) trim(system_command)
  call system(trim(system_command))

  write(numstring,*) OX+(NX-1)*DX
  write(tmpstring,*) 'xmax                            =',trim(numstring)
  tmpstring = adjustl(tmpstring)
  !write(*,*) trim(tmpstring)
  write(system_command,*) 'sed -i "s/^xmax .*/',trim(tmpstring),'/" ./DATA/Par_file'
  write (*,*) trim(system_command)
  call system(system_command)

  write(numstring,*) (NX-1)/INTERPOLATION_POINTS
  write(tmpstring,*) 'nx                              =',trim(numstring)
  tmpstring = adjustl(tmpstring)
  !write(*,*) trim(tmpstring)
  write(system_command,*) 'sed -i "s/^nx .*/',trim(tmpstring),'/" ./DATA/Par_file'
  write (*,*) trim(system_command)
  call system(trim(system_command))

  ! nbregions
  !call system('sed -i "$ d" ./DATA/Par_file')
  write(tmpstring,*) '1',(NX-1)/INTERPOLATION_POINTS,'1',(NZ-1)/INTERPOLATION_POINTS,'1'
  tmpstring = adjustl(tmpstring)
  !write (*,*) trim(tmpstring)
  !write(system_command,*) 'echo ${tmpstring}'
  !write(system_command,*) 'sed -i "/nbregions   /a',trim(tmpstring),'" ./DATA/Par_file'
  ! original Par_file line for regions to replace is: 1 20 1 20 1
  write(system_command,*) 'sed -i "s/^1 20 .*/',trim(tmpstring),'/" ./DATA/Par_file'
  write (*,*) trim(system_command)
  call system(trim(system_command))

  ! number of processes
  if (USER_INPUT) then
    write(tmpstring,*) 'NPROC                           =',nproc
    tmpstring = adjustl(tmpstring)
    !write(*,*) trim(tmpstring)
    write(system_command,*) 'sed -i "s/^NPROC .*/',trim(tmpstring),'/" ./DATA/Par_file'
    write (*,*) trim(system_command)
    call system(trim(system_command))
  endif
  print *

  ! backup Par_file
  call system('cp ./DATA/Par_file ./DATA/Par_file.org')

  call sleep(1)

  ! note: first simulation is used to create a mesh and initial model files,
  !       which will be read in a second run for interpolation

  ! user output
  print *
  print *, '************ Setting up first simulation ************ '
  print *

  ! Par_file flags to create new mesh files
  call system('cp ./DATA/Par_file.org ./DATA/Par_file')
  call system('sed -i "s/^MODEL .*/MODEL                            = default/" ./DATA/Par_file')
  call system('sed -i "s/^SAVE_MODEL .*/SAVE_MODEL                  = legacy/" ./DATA/Par_file')
  call system('sed -i "s/^NSTEP .*/NSTEP                            = 5 /" ./DATA/Par_file')
  call system('cp ./DATA/Par_file ./OUTPUT_FILES/Par_file.step_1')

  ! runs mesher
  print *,'call xmeshfem2d'
  call system('./bin/xmeshfem2D > OUTPUT_FILES/output_mesher.step_1.txt')

  ! runs first forward simulation to generate new specfem files
  print *,'call xspecfem2d'
  write(num,'(i2.2)') nproc
  call system('mpirun -np ' //num//' ./bin/xspecfem2D > OUTPUT_FILES/output_solver.step_1.txt')

  print *
  print *, '************ Setting up new simulation ************ '
  print *

  ! now uses default specfem file format
  ! backup Par_file
  call system('cp ./DATA/Par_file.org ./DATA/Par_file')
  call system('sed -i "s/^MODEL .*/MODEL                            = legacy/" ./DATA/Par_file')
  call system('sed -i "s/^SAVE_MODEL .*/SAVE_MODEL                  = default/" ./DATA/Par_file')
  call system('cp ./DATA/Par_file ./OUTPUT_FILES/Par_file.step_2')

  ! runs mesher
  print *,'call xmeshfem2d'
  call system('./bin/xmeshfem2D')

  ! scaling model values
  if (SCALE_FROM_VP) then
    ! VS model, scaled or in SEP format
    vs_SEP(:,:) = vp_SEP(:,:) / 1.732

    ! RHO model, scaled or in SEP format
    ! note: by default, just assumes constant density
    ! constant density
    rho_SEP(:,:) =  1000.0

    ! scaling rule
    !rho_SEP  =   (1.6612 * (vp_SEP/14500*4500 / 1000.0)     &
    !             -0.4720 * (vp_SEP/14500*4500 / 1000.0)**2  &
    !             +0.0671 * (vp_SEP/14500*4500 / 1000.0)**3  &
    !             -0.0043 * (vp_SEP/14500*4500 / 1000.0)**4  &
    !             +0.000106*(vp_SEP/14500*4500 / 1000.0)**5 )*1000.0
  endif

  ! user output
  print *
  print *, '************ Interpolating model ************'
  print *

  print *, 'input model:'
  print *, '  Vp min/max  = ',minval(vp_SEP(:,:)),maxval(vp_SEP(:,:))
  print *, '  Vs min/max  = ',minval(vs_SEP(:,:)),maxval(vs_SEP(:,:))
  print *, '  rho min/max = ',minval(rho_SEP(:,:)),maxval(rho_SEP(:,:))

  call interpolate_slowness_gll(vp_SEP,vs_SEP,rho_SEP,NX,NY,NZ,DX,DY,DZ,OX,OY,OZ,INTERPOLATE_FROM_CLOSEST_POINT,nproc)

  ! user output
  print *
  print *, 'Interpolation done'
  print *

  ! runs forward simulation
  print *,'call xspecfem2d'
  call system('mpirun -np ' //num//' ./bin/xspecfem2D ')

  end program main_interpolate

!
!-----------------------------------------------------------------------------------
!

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

!
!-----------------------------------------------------------------------------------
!

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
  integer :: ier,idummy
  character(len=1) :: strdummy

  ! to avoid compiler warning
  idummy = int(OX * OY * OZ)
  idummy = int(DX * DY * DZ)
  idummy = esize
  strdummy(1:1) = data_format(1:1)

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

!
!-----------------------------------------------------------------------------------
!

  subroutine interpolate_slowness_gll(vp,vs,rho,NX,NY,NZ,DX,DY,DZ,OX,OY,OZ,closest_point,nproc)

  implicit none

  !!! constants
  integer, parameter :: NGLLX = 5, NGLLZ = 5
  !double precision, parameter :: vp_water = 5500 ! upper bound

  !!! input
  integer,intent(in) :: NX,NY,NZ
  real,intent(in)    :: OX,OY,OZ
  real,intent(in)    :: DX,DY,DZ
  real,dimension(NX,NZ),intent(in) :: vp,vs,rho
  logical,intent(in) :: closest_point
  integer,intent(in) :: nproc

  !!! output
  ! currently no output, all information is saved in output file

  ! local parameters
  integer :: iproc
  integer :: ier,ix,iz,ix1,iz1,i,npoints
  integer, dimension(NGLLX*NGLLZ) :: iglob
  double precision, dimension(NGLLX*NGLLZ) :: rho_new,vp_new,vs_new,x,z
  double precision :: rho_temp,vp_temp,vs_temp,a1,b1,tmp1,tmp2,tmp3
  character(len=512) :: input_file,output_file
  character(len=512) :: system_command
  double precision :: vp_min,vp_max,vs_min,vs_max,rho_min,rho_max
  real :: rdummy

  ! to avoid compiler warning
  rdummy = OY * DY * NY

  print *
  print *, 'Interpolating GLL model values'
  print *

  ! interpolating model values onto all GLL points
  vp_min = 1.d30
  vp_max = -1.d30
  vs_min = 1.d30
  vs_max = - 1.d30
  rho_min = 1.d30
  rho_max = - 1.d30

  ! for each process
  do iproc = 1,nproc

    ! user output
    print *, 'iproc = ',iproc

    ! input file
    write(input_file,'(a,i6.6,a)') 'DATA/proc',iproc-1,'_model_velocity.dat_input'
    open(unit=15,file=trim(adjustl(input_file)),status='old',iostat=ier)
    if (ier /= 0) stop 'Error opening DATA/proc*****_model_velocity.dat_input file.'

    ! counts number of points
    npoints = 0
    do while (ier == 0)
      ! loops over GLL points of every element
      do i = 1,NGLLX*NGLLZ
          ! reads model point
          read(15,'(I10,5e15.5e4)',iostat=ier)  iglob(i),x(i),z(i),rho_temp,vp_temp,vs_temp
          if (ier /= 0) exit

          npoints = npoints + 1
      enddo
    enddo
    rewind(15)
    print *, 'number of points = ',npoints

    ! number of lines must be a multiple of NGLLX * NGLLZ
    if (mod(npoints,(NGLLX*NGLLZ)) /= 0) then
      print *,'Error: invalid number of lines in input file '//trim(input_file)
      print *,'  number of points = ',npoints,' should be a multiple of NGLLX * NGLLZ = ',NGLLX,' * ',NGLLZ
      stop 'Error invalid number of lines in input file'
    endif

    ! output file
    write(output_file,'(a,i6.6,a)') 'DATA/proc',iproc-1,'_model_velocity.dat_tmp'
    open(unit=16,file=trim(adjustl(output_file)),status='unknown',iostat=ier)
    if (ier /= 0) stop 'Error opening output file DATA/proc****_model_velocity.dat_tmp for interpolation.'

    ! interpolation
    if (closest_point) then

      ! takes value from closest point
      do while (ier == 0)
        ! loops over GLL points of every element
        do i = 1,NGLLX*NGLLZ
          ! reads model point
          read(15,'(I10,5e15.5e4)',iostat=ier)  iglob(i),x(i),z(i),rho_temp,vp_temp,vs_temp
          if (ier /= 0) exit

          ix = NINT((x(i)-OX)/DX) + 1
          iz = NINT((z(i)-OZ)/DZ) + 1

          if (ix > NX) ix = NX
          if (ix < 1)  ix = 1
          if (iz > NZ) iz = NZ
          if (iz < 1)  iz = 1

          ! closest point values
          rho_new(i) = rho(ix,iz)
          vp_new(i) = vp(ix,iz)
          vs_new(i) = vs(ix,iz)

          ! statistics
          if (rho_new(i) < rho_min) rho_min = rho_new(i)
          if (rho_new(i) > rho_max) rho_max = rho_new(i)

          if (vp_new(i) < vp_min) vp_min = vp_new(i)
          if (vp_new(i) > vp_max) vp_max = vp_new(i)

          if (vs_new(i) < vs_min) vs_min = vs_new(i)
          if (vs_new(i) > vs_max) vs_max = vs_new(i)
        enddo

        ! writes out interpolated values
        ! converting from slowness to wave speed values again & buoyancy to density
        if (ier == 0) then
          do i = 1,NGLLX*NGLLZ
            write(16,'(I10,5e15.5e4)') iglob(i),x(i),z(i), rho_new(i), vp_new(i), vs_new(i)
          enddo
        endif
      enddo

    else

      ! bilinear interpolation
      !
      !  (ix,iz)    (ix+1,iz)
      !   x---------- x
      !   |   *       |
      !   |           |
      !   |           |
      !   x-----------x
      !  (ix,iz+1)   (ix+1,iz+1)
      !
      do while (ier == 0)
        ! loops over GLL points of every element
        do i = 1,NGLLX*NGLLZ
          ! reads model point
          read(15,'(I10,5e15.5e4)',iostat=ier)  iglob(i),x(i),z(i),rho_temp,vp_temp,vs_temp
          if (ier /= 0) exit

          ix = INT((x(i)-OX)/DX) + 1
          iz = INT((z(i)-OZ)/DZ) + 1
          ix1 = ix+1
          iz1 = iz+1

          if (ix1 > NX) ix1 = NX
          if (iz1 > NZ) iz1 = NZ

          a1 = x(i)-OX-(ix-1)*DX
          a1 = a1/DX
          b1 = z(i)-OZ-(iz-1)*DZ
          b1 = b1/DZ

          ! note:
          ! - Vp: we use slowness values for interpolation. this garantees better interpolation results for traveltimes.
          ! - Vs: shear speeds might be zero in a coupled fluid-elastic model, so it is inverted directly.
          ! - density: it will be inverted to have "lightness" or "buoyancy" for interpolation.
          !            this is also done in a staggered-grid approach, see e.g.,
          !            Virieux 1987, "P-SV wave propagation in heterogeneous media: Velocity-stress finite-difference method",
          !            Geophysics, 51, p. 889 - 901.

          ! slowness
          tmp1 = (1.d0-a1)*1.d0/vp(ix,iz) + a1*1.d0/vp(ix1,iz)
          tmp2 = (1.d0-a1)*1.d0/vp(ix,iz1) + a1*1.d0/vp(ix1,iz1)
          tmp3 = (1.d0-b1)*tmp1 + b1*tmp2
          ! converts back to Vp
          vp_new(i) = 1.d0/tmp3

          ! Vs
          tmp1 = (1.d0-a1)*vs(ix,iz) + a1*vs(ix1,iz)
          tmp2 = (1.d0-a1)*vs(ix,iz1) + a1*vs(ix1,iz1)
          tmp3 = (1.d0-b1)*tmp1 + b1*tmp2
          vs_new(i) = tmp3

          ! buoyancy
          tmp1 = (1.d0-a1)*1.d0/rho(ix,iz) + a1*1.d0/rho(ix1,iz)
          tmp2 = (1.d0-a1)*1.d0/rho(ix,iz1) + a1*1.d0/rho(ix1,iz1)
          tmp3 = (1.d0-b1)*tmp1 + b1*tmp2
          ! converts back to density
          rho_new(i) = 1.d0/tmp3

          ! statistics
          if (rho_new(i) < rho_min) rho_min = rho_new(i)
          if (rho_new(i) > rho_max) rho_max = rho_new(i)

          if (vp_new(i) < vp_min) vp_min = vp_new(i)
          if (vp_new(i) > vp_max) vp_max = vp_new(i)

          if (vs_new(i) < vs_min) vs_min = vs_new(i)
          if (vs_new(i) > vs_max) vs_max = vs_new(i)
        enddo

        ! writes out interpolated values
        if (ier == 0) then
          do i = 1,NGLLX*NGLLZ
            write(16,'(I10,5e15.5e4)') iglob(i),x(i),z(i), rho_new(i), vp_new(i), vs_new(i)
          enddo
        endif
      enddo
    endif

    ! closes files
    close(15)
    close(16)

    ! replaces model file
    print *, 'replacing model file'
    write(system_command,*) 'mv -f ',trim(adjustl(input_file)),' ',trim(adjustl(input_file))//'.org'
    write(*,*) trim(system_command)
    call system(trim(system_command))
    write(system_command,*) 'mv -f ',trim(adjustl(output_file)),' ',trim(adjustl(input_file))
    write(*,*) trim(system_command)
    call system(trim(system_command))

  enddo ! iproc

  ! user output
  print *
  print *, 'interpolated output model:'
  print *, '  Vp min/max  = ',sngl(vp_min),sngl(vp_max)
  print *, '  Vs min/max  = ',sngl(vs_min),sngl(vs_max)
  print *, '  rho min/max = ',sngl(rho_min),sngl(rho_max)
  print *

  end subroutine interpolate_slowness_gll


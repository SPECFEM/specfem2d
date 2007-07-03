
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
!                         Dimitri Komatitsch
!                     University of Pau, France
!
!                          (c) April 2007
!
!========================================================================

!========================================================================
!
!  Basic mesh generator for SPECFEM2D
!
!========================================================================

! If you use this code for your own research, please cite:
!
! @ARTICLE{KoTr99,
! author={D. Komatitsch and J. Tromp},
! title={Introduction to the spectral-element method for 3-{D} seismic wave propagation},
! journal={Geophys. J. Int.},
! year=1999,
! volume=139,
! number=3,
! pages={806-822},
! doi={10.1046/j.1365-246x.1999.00967.x}}
!
! @ARTICLE{KoVi98,
! author={D. Komatitsch and J. P. Vilotte},
! title={The spectral-element method: an efficient tool to simulate the seismic response of 2{D} and 3{D} geological structures},
! journal={Bull. Seismol. Soc. Am.},
! year=1998,
! volume=88,
! number=2,
! pages={368-392}}

program meshfem2D

  use part_unstruct
  implicit none

  include "constants.h"

! coordinates of the grid points of the mesh
  double precision, dimension(:,:), allocatable :: x,z

! to compute the coordinate transformation
  integer :: ioffset
  double precision :: gamma,absx,a00,a01,bot0,top0

! to store density and velocity model
  double precision, dimension(:), allocatable :: rho,cp,cs,aniso3,aniso4
  integer, dimension(:), allocatable :: icodemat
  integer, dimension(:), allocatable :: num_material

! interface data
  integer interface_current,ipoint_current,number_of_interfaces,npoints_interface_bottom,npoints_interface_top
  integer ilayer,number_of_layers,max_npoints_interface
  double precision xinterface_dummy,zinterface_dummy,xinterface_dummy_previous
  integer, dimension(:), allocatable :: nz_layer
  double precision, dimension(:), allocatable :: &
         xinterface_bottom,zinterface_bottom,coefs_interface_bottom, &
         xinterface_top,zinterface_top,coefs_interface_top

! for the source and receivers
  integer source_type,time_function_type,nrec_total,irec_global_number
  double precision xs,zs,f0,t0,angleforce,Mxx,Mzz,Mxz,factor,xrec,zrec

  character(len=50) interfacesfile,title

  integer imaterial_number,inumelem
  integer nelemabs,nelem_acoustic_surface,npgeo,nspec
  integer ix,iz,irec,i,j
  integer icol, ili, istepx, istepz
  integer ixdebregion,ixfinregion,izdebregion,izfinregion
  integer iregion,imaterial,nbregion,nb_materials
  integer NTSTEP_BETWEEN_OUTPUT_INFO,pointsdisp,subsamp,seismotype,imagetype
  logical generate_STATIONS 
  integer ngnod,nt,nx,nz,nxread,nzread,icodematread,ireceiverlines,nreceiverlines

  integer, dimension(:), allocatable :: nrec

  logical output_postscript_snapshot,output_color_image,plot_lowerleft_corner_only

  double precision tang1,tangN,vpregion,vsregion,poisson_ratio
  double precision cutsnaps,sizemax_arrows,anglerec,xmin,xmax,deltat
  double precision rhoread,cpread,csread,aniso3read,aniso4read

  double precision, dimension(:), allocatable :: xdeb,zdeb,xfin,zfin

  logical interpol,gnuplot,assign_external_model,outputgrid
  logical abstop,absbottom,absleft,absright,any_abs
  logical source_surf,meshvect,initialfield,modelvect,boundvect
  logical TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON

  logical, dimension(:), allocatable :: enreg_surf

  integer, external :: num_4, num_9
  double precision, external :: value_spline

! flag to indicate an anisotropic material
  integer, parameter :: ANISOTROPIC_MATERIAL = 2

! file number for interface file
  integer, parameter :: IIN_INTERFACES = 15

! ignore variable name field (junk) at the beginning of each input line
  logical, parameter :: IGNORE_JUNK = .true.,DONT_IGNORE_JUNK = .false.

! parameters for external mesh
  logical  :: read_external_mesh
  character(len=256)  :: mesh_file, nodes_coords_file, materials_file, free_surface_file, absorbing_surface_file
  
! variables used for storing info about the mesh and partitions    
  integer, dimension(:), pointer  :: elmnts
  integer, dimension(:), pointer  :: elmnts_bis
  double precision, dimension(:,:), pointer  :: nodes_coords
  integer, dimension(:,:), pointer  :: acoustic_surface
  integer, dimension(:,:), pointer  :: abs_surface
  
  integer, dimension(:), pointer  :: xadj
  integer, dimension(:), pointer  :: adjncy
  integer, dimension(:), pointer  :: nnodes_elmnts
  integer, dimension(:), pointer  :: nodes_elmnts
  
  integer, dimension(:), pointer  :: vwgt
  integer, dimension(:), pointer  :: adjwgt
  integer, dimension(:), pointer  :: part

  integer, dimension(:), pointer  :: glob2loc_elmnts
  integer, dimension(:), pointer  :: glob2loc_nodes_nparts
  integer, dimension(:), pointer  :: glob2loc_nodes_parts
  integer, dimension(:), pointer  :: glob2loc_nodes
  integer, dimension(:), pointer  :: tab_size_interfaces, tab_interfaces
  integer, dimension(:), allocatable  :: my_interfaces
  integer, dimension(:), allocatable  :: my_nb_interfaces
    
  integer  :: nedges_coupled, nedges_coupled_loc
  integer, dimension(:,:), pointer  :: edges_coupled

  integer  :: num_start
  integer  :: nelmnts
  integer  :: num_elmnt
  integer  :: nnodes
  integer  :: num_node
  integer  :: nb_edges
  integer  ::  ninterfaces
  integer  :: my_ninterface
  integer  :: nelem_acoustic_surface_loc, nelemabs_loc
  logical, dimension(:,:), pointer  :: abs_surface_char
  integer, dimension(:), pointer  :: abs_surface_merge
  integer  :: nelemabs_merge
  integer, dimension(:), pointer  :: ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
       jbegin_left,jend_left,jbegin_right,jend_right

! variables used for partitionning
  integer  :: nproc
  integer  :: partitionning_method
  character(len=256)  :: partitionning_strategy
  character(len=256)  :: scotch_strategy
  integer, dimension(0:4)  :: metis_options
  character(len=256)  :: prname

#if defined USE_METIS || defined USE_SCOTCH
  integer  :: edgecut
#endif

  integer  :: iproc



! ***
! *** read the parameter file
! ***

  print *,'Reading the parameter file ... '
  print *

  open(unit=IIN,file='DATA/Par_file',status='old')

! read file names and path for output
  call read_value_string(IIN,IGNORE_JUNK,title)
  call read_value_string(IIN,IGNORE_JUNK,interfacesfile)

  write(*,*) 'Title of the simulation'
  write(*,*) title
  print *

! read info about external mesh
  call read_value_logical(IIN,IGNORE_JUNK,read_external_mesh)
  call read_value_string(IIN,IGNORE_JUNK,mesh_file)
  call read_value_string(IIN,IGNORE_JUNK,nodes_coords_file)
  call read_value_string(IIN,IGNORE_JUNK,materials_file)
  call read_value_string(IIN,IGNORE_JUNK,free_surface_file)
  call read_value_string(IIN,IGNORE_JUNK,absorbing_surface_file)

! read info about partitionning
  call read_value_integer(IIN,IGNORE_JUNK,nproc)
  if ( nproc <= 0 ) then
     print *, 'Number of processes (nproc) must be greater than or equal to one.' 
     stop 
  end if
  
#ifndef USE_MPI
  if ( nproc > 1 ) then
     print *, 'Number of processes (nproc) must be equal to one when not using MPI.' 
     print *, 'Please recompile with -DUSE_MPI in order to enable use of MPI.'
     stop 
  end if
  
#endif
  
  call read_value_integer(IIN,IGNORE_JUNK,partitionning_method)
  call read_value_string(IIN,IGNORE_JUNK,partitionning_strategy)
  select case(partitionning_method) 
  case(1)
  case(2)
     partitionning_strategy = trim(partitionning_strategy)
     if ( partitionning_strategy(1:1) == '0' ) then
        metis_options(0) = 0
     else
        do i = 1, 5
           metis_options = iachar(partitionning_strategy(i:i)) - iachar('0')
        end do
     end if
     
  case(3)
     scotch_strategy = trim(partitionning_strategy)
     
  case default 
     print *, 'Invalid partionning method number.'
     print *, 'Partionning method', partitionning_method, 'was requested, but is not available.'
     stop
  end select
  
 
! read grid parameters
  call read_value_double_precision(IIN,IGNORE_JUNK,xmin)
  call read_value_double_precision(IIN,IGNORE_JUNK,xmax)
  call read_value_integer(IIN,IGNORE_JUNK,nx)
  call read_value_integer(IIN,IGNORE_JUNK,ngnod)
  if ( ngnod == 9 .and. read_external_mesh ) then
     print *, 'Number of control nodes must be equal to four when reading from external mesh.'
     print *, 'ngnod = 9 is not yet supported.'
     stop 
  end if
  
  call read_value_logical(IIN,IGNORE_JUNK,initialfield)
  call read_value_logical(IIN,IGNORE_JUNK,assign_external_model)
  call read_value_logical(IIN,IGNORE_JUNK,TURN_ANISOTROPY_ON)
  call read_value_logical(IIN,IGNORE_JUNK,TURN_ATTENUATION_ON)
  
  if ( read_external_mesh ) then
     call read_mesh(mesh_file, nelmnts, elmnts, nnodes, num_start)
     
  else
     ! get interface data from external file to count the spectral elements along Z
     print *,'Reading interface data from file DATA/',interfacesfile(1:len_trim(interfacesfile)),' to count the spectral elements'
     open(unit=IIN_INTERFACES,file='DATA/'//interfacesfile,status='old')
     
     max_npoints_interface = -1
     
     ! read number of interfaces
     call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,number_of_interfaces)
     if(number_of_interfaces < 2) stop 'not enough interfaces (minimum is 2)'

     ! loop on all the interfaces
     do interface_current = 1,number_of_interfaces
        
        call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_bottom)
        if(npoints_interface_bottom < 2) stop 'not enough interface points (minimum is 2)'
        max_npoints_interface = max(npoints_interface_bottom,max_npoints_interface)
        print *,'Reading ',npoints_interface_bottom,' points for interface ',interface_current
        
        ! loop on all the points describing this interface
        do ipoint_current = 1,npoints_interface_bottom
           call read_two_interface_points(IIN_INTERFACES,DONT_IGNORE_JUNK,xinterface_dummy,zinterface_dummy)
           if(ipoint_current > 1 .and. xinterface_dummy <= xinterface_dummy_previous) &
                stop 'interface points must be sorted in increasing X'
           xinterface_dummy_previous = xinterface_dummy
        enddo
        
     enddo
     
     ! define number of layers
     number_of_layers = number_of_interfaces - 1
     
     allocate(nz_layer(number_of_layers))
     
     ! loop on all the layers
     do ilayer = 1,number_of_layers
        
        ! read number of spectral elements in vertical direction in this layer
        call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,nz_layer(ilayer))
        if(nz_layer(ilayer) < 1) stop 'not enough spectral elements along Z in layer (minimum is 1)'
        print *,'There are ',nz_layer(ilayer),' spectral elements along Z in layer ',ilayer
        
     enddo

     close(IIN_INTERFACES)
     
     ! compute total number of spectral elements in vertical direction
     nz = sum(nz_layer)
     
     print *
     print *,'Total number of spectral elements along Z = ',nz
     print *
     
     nxread = nx
     nzread = nz
     
     ! multiply by 2 if elements have 9 nodes
     if(ngnod == 9) then
        nx = nx * 2
        nz = nz * 2
        nz_layer(:) = nz_layer(:) * 2
     endif
     
     nelmnts = nxread * nzread
     allocate(elmnts(0:ngnod*nelmnts-1))
     if ( ngnod == 4 ) then
        num_elmnt = 0
        do j = 1, nzread
           do i = 1, nxread
              elmnts(num_elmnt*ngnod)   = (j-1)*(nxread+1) + (i-1)
              elmnts(num_elmnt*ngnod+1) = (j-1)*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*ngnod+2) = j*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*ngnod+3) = j*(nxread+1) + (i-1)
              num_elmnt = num_elmnt + 1
           end do
        end do
      else
        num_elmnt = 0
        do j = 1, nzread
           do i = 1, nxread
              elmnts(num_elmnt*ngnod)   = (j-1)*(nxread+1) + (i-1)
              elmnts(num_elmnt*ngnod+1) = (j-1)*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*ngnod+2) = j*(nxread+1) + (i-1) + 1
              elmnts(num_elmnt*ngnod+3) = j*(nxread+1) + (i-1) 
              elmnts(num_elmnt*ngnod+4) = (nxread+1)*(nzread+1) + (j-1)*nxread + (i-1)
              elmnts(num_elmnt*ngnod+5) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2 + 2
              elmnts(num_elmnt*ngnod+6) = (nxread+1)*(nzread+1) + j*nxread + (i-1)
              elmnts(num_elmnt*ngnod+7) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2
              elmnts(num_elmnt*ngnod+8) = (nxread+1)*(nzread+1) + nxread*(nzread+1) + (j-1)*(nxread*2+1) + (i-1)*2 + 1
              num_elmnt = num_elmnt + 1
           end do
        end do
        
     end if
  end if

! read absorbing boundaries parameters
  call read_value_logical(IIN,IGNORE_JUNK,any_abs)
  call read_value_logical(IIN,IGNORE_JUNK,absbottom)
  call read_value_logical(IIN,IGNORE_JUNK,absright)
  call read_value_logical(IIN,IGNORE_JUNK,abstop)
  call read_value_logical(IIN,IGNORE_JUNK,absleft)
  if ( .not. any_abs ) then
     absbottom = .false.
     absright = .false.
     abstop = .false.
     absleft = .false.
  end if

! read time step parameters
  call read_value_integer(IIN,IGNORE_JUNK,nt)
  call read_value_double_precision(IIN,IGNORE_JUNK,deltat)

! read source parameters
  call read_value_logical(IIN,IGNORE_JUNK,source_surf)
  call read_value_double_precision(IIN,IGNORE_JUNK,xs)
  call read_value_double_precision(IIN,IGNORE_JUNK,zs)
  call read_value_integer(IIN,IGNORE_JUNK,source_type)
  call read_value_integer(IIN,IGNORE_JUNK,time_function_type)
  call read_value_double_precision(IIN,IGNORE_JUNK,f0)
  call read_value_double_precision(IIN,IGNORE_JUNK,angleforce)
  call read_value_double_precision(IIN,IGNORE_JUNK,Mxx)
  call read_value_double_precision(IIN,IGNORE_JUNK,Mzz)
  call read_value_double_precision(IIN,IGNORE_JUNK,Mxz)
  call read_value_double_precision(IIN,IGNORE_JUNK,factor)

! if Dirac source time function, use a very thin Gaussian instead
! if Heaviside source time function, use a very thin error function instead
  if(time_function_type == 4 .or. time_function_type == 5) f0 = 1.d0 / (10.d0 * deltat)

! time delay of the source in seconds, use a 20 % security margin (use 2 / f0 if error function)
  if(time_function_type == 5) then
    t0 = 2.0d0 / f0
  else
    t0 = 1.20d0 / f0
  endif

  print *
  print *,'Source:'
  print *,'Position xs, zs = ',xs,zs
  print *,'Frequency, delay = ',f0,t0
  print *,'Source type (1=force, 2=explosion): ',source_type
  print *,'Time function type (1=Ricker, 2=First derivative, 3=Gaussian, 4=Dirac, 5=Heaviside): ',time_function_type
  print *,'Angle of the source if force = ',angleforce
  print *,'Mxx of the source if moment tensor = ',Mxx
  print *,'Mzz of the source if moment tensor = ',Mzz
  print *,'Mxz of the source if moment tensor = ',Mxz
  print *,'Multiplying factor = ',factor

! read receiver line parameters
  call read_value_integer(IIN,IGNORE_JUNK,seismotype)
  call read_value_logical(IIN,IGNORE_JUNK,generate_STATIONS)
  call read_value_integer(IIN,IGNORE_JUNK,nreceiverlines)
  call read_value_double_precision(IIN,IGNORE_JUNK,anglerec)

  if(nreceiverlines < 1) stop 'number of receiver lines must be greater than 1'

! allocate receiver line arrays
  allocate(nrec(nreceiverlines))
  allocate(xdeb(nreceiverlines))
  allocate(zdeb(nreceiverlines))
  allocate(xfin(nreceiverlines))
  allocate(zfin(nreceiverlines))
  allocate(enreg_surf(nreceiverlines))

! loop on all the receiver lines
  do ireceiverlines = 1,nreceiverlines
    call read_value_integer(IIN,IGNORE_JUNK,nrec(ireceiverlines))
    call read_value_double_precision(IIN,IGNORE_JUNK,xdeb(ireceiverlines))
    call read_value_double_precision(IIN,IGNORE_JUNK,zdeb(ireceiverlines))
    call read_value_double_precision(IIN,IGNORE_JUNK,xfin(ireceiverlines))
    call read_value_double_precision(IIN,IGNORE_JUNK,zfin(ireceiverlines))
    call read_value_logical(IIN,IGNORE_JUNK,enreg_surf(ireceiverlines))
  enddo

! read display parameters
  call read_value_integer(IIN,IGNORE_JUNK,NTSTEP_BETWEEN_OUTPUT_INFO)
  call read_value_logical(IIN,IGNORE_JUNK,output_postscript_snapshot)
  call read_value_logical(IIN,IGNORE_JUNK,output_color_image)
  call read_value_integer(IIN,IGNORE_JUNK,imagetype)
  call read_value_double_precision(IIN,IGNORE_JUNK,cutsnaps)
  call read_value_logical(IIN,IGNORE_JUNK,meshvect)
  call read_value_logical(IIN,IGNORE_JUNK,modelvect)
  call read_value_logical(IIN,IGNORE_JUNK,boundvect)
  call read_value_logical(IIN,IGNORE_JUNK,interpol)
  call read_value_integer(IIN,IGNORE_JUNK,pointsdisp)
  call read_value_integer(IIN,IGNORE_JUNK,subsamp)
  call read_value_double_precision(IIN,IGNORE_JUNK,sizemax_arrows)
  call read_value_logical(IIN,IGNORE_JUNK,gnuplot)
  call read_value_logical(IIN,IGNORE_JUNK,outputgrid)

! can use only one point to display lower-left corner only for interpolated snapshot
  if(pointsdisp < 3) then
    pointsdisp = 3
    plot_lowerleft_corner_only = .true.
  else
    plot_lowerleft_corner_only = .false.
  endif

! read the different material materials
  call read_value_integer(IIN,IGNORE_JUNK,nb_materials)
  if(nb_materials <= 0) stop 'Negative number of materials not allowed!'

  allocate(icodemat(nb_materials))
  allocate(rho(nb_materials))
  allocate(cp(nb_materials))
  allocate(cs(nb_materials))
  allocate(aniso3(nb_materials))
  allocate(aniso4(nb_materials))
  allocate(num_material(nelmnts))

  icodemat(:) = 0
  rho(:) = 0.d0
  cp(:) = 0.d0
  cs(:) = 0.d0
  aniso3(:) = 0.d0
  aniso4(:) = 0.d0
  num_material(:) = 0

  do imaterial=1,nb_materials
    call read_material_parameters(IIN,DONT_IGNORE_JUNK,i,icodematread,rhoread,cpread,csread,aniso3read,aniso4read)
    if(i < 1 .or. i > nb_materials) stop 'Wrong material number!'
    icodemat(i) = icodematread
    rho(i) = rhoread
    cp(i) = cpread
    cs(i) = csread

    if(rho(i) <= 0.d0 .or. cp(i) <= 0.d0 .or. cs(i) < 0.d0) stop 'negative value of velocity or density'

    aniso3(i) = aniso3read
    aniso4(i) = aniso4read
  enddo
 
  print *
  print *, 'Nb of solid or fluid materials = ',nb_materials
  print *
  do i=1,nb_materials
    if(icodemat(i) /= ANISOTROPIC_MATERIAL) then
      print *,'Material #',i,' isotropic'
      print *,'rho,cp,cs = ',rho(i),cp(i),cs(i)
      if(cs(i) < TINYVAL) then
        print *,'Material is fluid'
      else
        print *,'Material is solid'
      endif
    else
      print *,'Material #',i,' anisotropic'
      print *,'rho,c11,c13,c33,c44 = ',rho(i),cp(i),cs(i),aniso3(i),aniso4(i)
    endif
  print *
  enddo

  
  if ( read_external_mesh ) then
     call read_mat(materials_file, nelmnts, num_material)
  else
     ! read the material numbers for each region
     call read_value_integer(IIN,IGNORE_JUNK,nbregion)

     if(nbregion <= 0) stop 'Negative number of regions not allowed!'

     print *
     print *, 'Nb of regions in the mesh = ',nbregion
     print *

     do iregion = 1,nbregion

        call read_region_coordinates(IIN,DONT_IGNORE_JUNK,ixdebregion,ixfinregion,izdebregion,izfinregion,imaterial_number)

        if(imaterial_number < 1) stop 'Negative material number not allowed!'
        if(ixdebregion < 1) stop 'Left coordinate of region negative!'
        if(ixfinregion > nxread) stop 'Right coordinate of region too high!'
        if(izdebregion < 1) stop 'Bottom coordinate of region negative!'
        if(izfinregion > nzread) stop 'Top coordinate of region too high!'

        print *,'Region ',iregion
        print *,'IX from ',ixdebregion,' to ',ixfinregion
        print *,'IZ from ',izdebregion,' to ',izfinregion

        if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL) then
           vpregion = cp(imaterial_number)
           vsregion = cs(imaterial_number)
           print *,'Material # ',imaterial_number,' isotropic'
           if(vsregion < TINYVAL) then
              print *,'Material is fluid'
           else
              print *,'Material is solid'
           endif
           print *,'vp = ',vpregion
           print *,'vs = ',vsregion
           print *,'rho = ',rho(imaterial_number)
           poisson_ratio = 0.5d0*(vpregion*vpregion-2.d0*vsregion*vsregion) / (vpregion*vpregion-vsregion*vsregion)
           print *,'Poisson''s ratio = ',poisson_ratio
           if(poisson_ratio <= -1.00001d0 .or. poisson_ratio >= 0.50001d0) stop 'incorrect value of Poisson''s ratio'
        else
           print *,'Material # ',imaterial_number,' anisotropic'
           print *,'c11 = ',cp(imaterial_number)
           print *,'c13 = ',cs(imaterial_number)
           print *,'c33 = ',aniso3(imaterial_number)
           print *,'c44 = ',aniso4(imaterial_number)
           print *,'rho = ',rho(imaterial_number)
        endif
        print *,' -----'

        ! store density and velocity model
        do i = ixdebregion,ixfinregion
           do j = izdebregion,izfinregion
              num_material((j-1)*nxread+i) = imaterial_number
           enddo
        enddo
        
     enddo

     if(minval(num_material) <= 0) stop 'Velocity model not entirely set...'

  end if

  close(IIN)

  print *
  print *,'Parameter file successfully read... '

!---

  if(ngnod /= 4 .and. ngnod /= 9) stop 'ngnod different from 4 or 9!'

  print *
  print *,'The mesh contains ',nelmnts,' elements'
  print *
  print *,'Control elements have ',ngnod,' nodes'
  print *

!---

  if ( .not. read_external_mesh ) then
     ! allocate arrays for the grid
     allocate(x(0:nx,0:nz))
     allocate(z(0:nx,0:nz))

     x(:,:) = 0.d0
     z(:,:) = 0.d0

     ! get interface data from external file
     print *,'Reading interface data from file DATA/',interfacesfile(1:len_trim(interfacesfile))
     open(unit=IIN_INTERFACES,file='DATA/'//interfacesfile,status='old')

     allocate(xinterface_bottom(max_npoints_interface))
     allocate(zinterface_bottom(max_npoints_interface))
     allocate(coefs_interface_bottom(max_npoints_interface))

     allocate(xinterface_top(max_npoints_interface))
     allocate(zinterface_top(max_npoints_interface))
     allocate(coefs_interface_top(max_npoints_interface))

     ! read number of interfaces
     call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,number_of_interfaces)

     ! read bottom interface
     call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_bottom)

     ! loop on all the points describing this interface
     do ipoint_current = 1,npoints_interface_bottom
        call read_two_interface_points(IIN_INTERFACES,DONT_IGNORE_JUNK, &
             xinterface_bottom(ipoint_current),zinterface_bottom(ipoint_current))
     enddo

     ! loop on all the layers
     do ilayer = 1,number_of_layers

        ! read top interface
        call read_value_integer(IIN_INTERFACES,DONT_IGNORE_JUNK,npoints_interface_top)

        ! loop on all the points describing this interface
        do ipoint_current = 1,npoints_interface_top
           call read_two_interface_points(IIN_INTERFACES,DONT_IGNORE_JUNK, &
                xinterface_top(ipoint_current),zinterface_top(ipoint_current))
        enddo

        ! compute the spline for the bottom interface, impose the tangent on both edges
        tang1 = (zinterface_bottom(2)-zinterface_bottom(1)) / (xinterface_bottom(2)-xinterface_bottom(1))
        tangN = (zinterface_bottom(npoints_interface_bottom)-zinterface_bottom(npoints_interface_bottom-1)) / &
             (xinterface_bottom(npoints_interface_bottom)-xinterface_bottom(npoints_interface_bottom-1))
        call spline(xinterface_bottom,zinterface_bottom,npoints_interface_bottom,tang1,tangN,coefs_interface_bottom)

        ! compute the spline for the top interface, impose the tangent on both edges
        tang1 = (zinterface_top(2)-zinterface_top(1)) / (xinterface_top(2)-xinterface_top(1))
        tangN = (zinterface_top(npoints_interface_top)-zinterface_top(npoints_interface_top-1)) / &
             (xinterface_top(npoints_interface_top)-xinterface_top(npoints_interface_top-1))
        call spline(xinterface_top,zinterface_top,npoints_interface_top,tang1,tangN,coefs_interface_top)

        ! check if we are in the last layer, which contains topography,
        ! and modify the position of the source accordingly if it is located exactly at the surface
        if(source_surf .and. ilayer == number_of_layers) &
             zs = value_spline(xs,xinterface_top,zinterface_top,coefs_interface_top,npoints_interface_top)

        ! compute the offset of this layer in terms of number of spectral elements below along Z
        if(ilayer > 1) then
           ioffset = sum(nz_layer(1:ilayer-1))
        else
           ioffset = 0
        endif

        !--- definition of the mesh

        do ix = 0,nx

           ! evenly spaced points along X
           absx = xmin + (xmax - xmin) * dble(ix) / dble(nx)

           ! value of the bottom and top splines
           bot0 = value_spline(absx,xinterface_bottom,zinterface_bottom,coefs_interface_bottom,npoints_interface_bottom)
           top0 = value_spline(absx,xinterface_top,zinterface_top,coefs_interface_top,npoints_interface_top)

           do iz = 0,nz_layer(ilayer)

              ! linear interpolation between bottom and top
              gamma = dble(iz) / dble(nz_layer(ilayer))
              a00 = 1.d0 - gamma
              a01 = gamma

              ! coordinates of the grid points
              x(ix,iz + ioffset) = absx
              z(ix,iz + ioffset) = a00*bot0 + a01*top0

           enddo

        enddo

        ! the top interface becomes the bottom interface before switching to the next layer
        npoints_interface_bottom = npoints_interface_top
        xinterface_bottom(:) = xinterface_top(:)
        zinterface_bottom(:) = zinterface_top(:)

     enddo

     close(IIN_INTERFACES)
     
     nnodes = (nz+1)*(nx+1)
     allocate(nodes_coords(2,nnodes))
     if ( ngnod == 4 ) then
        do j = 0, nz
           do i = 0, nx
              num_node = num_4(i,j,nxread)
              nodes_coords(1, num_node) = x(i,j)
              nodes_coords(2, num_node) = z(i,j)
              
           end do
        end do
        
     else
        do j = 0, nz
           do i = 0, nx
              num_node = num_9(i,j,nxread,nzread)
              nodes_coords(1, num_node) = x(i,j)
              nodes_coords(2, num_node) = z(i,j)
              
           end do
        end do
        
     end if
  else
     call read_nodes_coords(nodes_coords_file, nnodes, nodes_coords)
  end if
  

  if ( read_external_mesh ) then
     call read_acoustic_surface(free_surface_file, nelem_acoustic_surface, acoustic_surface, &
          nelmnts, num_material, ANISOTROPIC_MATERIAL, nb_materials, icodemat, cs, num_start)
     
     if ( any_abs ) then
        call read_abs_surface(absorbing_surface_file, nelemabs, abs_surface, num_start)
     end if
     
  else
     ! count the number of acoustic free-surface elements
     nelem_acoustic_surface = 0
     j = nzread
     do i = 1,nxread
        imaterial_number = num_material((j-1)*nxread+i)
        if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. cs(imaterial_number) < TINYVAL ) then
           nelem_acoustic_surface = nelem_acoustic_surface + 1
        end if
     enddo
     
     allocate(acoustic_surface(4,nelem_acoustic_surface))
     
     nelem_acoustic_surface = 0
     j = nzread
     do i = 1,nxread
        imaterial_number = num_material((j-1)*nxread+i)
        if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. cs(imaterial_number) < TINYVAL ) then
           nelem_acoustic_surface = nelem_acoustic_surface + 1
           acoustic_surface(1,nelem_acoustic_surface) = (j-1)*nxread + (i-1)
           acoustic_surface(2,nelem_acoustic_surface) = 2
           acoustic_surface(3,nelem_acoustic_surface) = elmnts(3+ngnod*((j-1)*nxread+i-1))
           acoustic_surface(4,nelem_acoustic_surface) = elmnts(2+ngnod*((j-1)*nxread+i-1))          
        end if
     end do
     
     !
     !--- definition of absorbing boundaries
     !
     nelemabs = 0
     if(absbottom) nelemabs = nelemabs + nxread
     if(abstop) nelemabs = nelemabs + nxread
     if(absleft) nelemabs = nelemabs + nzread
     if(absright) nelemabs = nelemabs + nzread
     
     allocate(abs_surface(4,nelemabs))
     
     ! generate the list of absorbing elements
     if(nelemabs > 0) then
        nelemabs = 0
        do iz = 1,nzread
           do ix = 1,nxread
              inumelem = (iz-1)*nxread + ix
              if(absbottom    .and. iz == 1) then
                 nelemabs = nelemabs + 1
                 abs_surface(1,nelemabs) = inumelem-1
                 abs_surface(2,nelemabs) = 2
                 abs_surface(3,nelemabs) = elmnts(0+ngnod*(inumelem-1))
                 abs_surface(4,nelemabs) = elmnts(1+ngnod*(inumelem-1))
              end if
              if(absright .and. ix == nxread) then
                 nelemabs = nelemabs + 1
                 abs_surface(1,nelemabs) = inumelem-1
                 abs_surface(2,nelemabs) = 2
                 abs_surface(3,nelemabs) = elmnts(1+ngnod*(inumelem-1))
                 abs_surface(4,nelemabs) = elmnts(2+ngnod*(inumelem-1))
              end if
              if(abstop   .and. iz == nzread) then
                 nelemabs = nelemabs + 1
                 abs_surface(1,nelemabs) = inumelem-1
                 abs_surface(2,nelemabs) = 2
                 abs_surface(3,nelemabs) = elmnts(3+ngnod*(inumelem-1)) 
                 abs_surface(4,nelemabs) = elmnts(2+ngnod*(inumelem-1)) 
              end if
              if(absleft .and. ix == 1) then
                 nelemabs = nelemabs + 1
                 abs_surface(1,nelemabs) = inumelem-1
                 abs_surface(2,nelemabs) = 2
                 abs_surface(3,nelemabs) = elmnts(0+ngnod*(inumelem-1)) 
                 abs_surface(4,nelemabs) = elmnts(3+ngnod*(inumelem-1)) 
              end if
           end do
        end do
     end if
     
  end if
  
  
! compute min and max of X and Z in the grid
  print *
  print *,'Min and max value of X in the grid = ',minval(nodes_coords(1,:)),maxval(nodes_coords(1,:))
  print *,'Min and max value of Z in the grid = ',minval(nodes_coords(2,:)),maxval(nodes_coords(2,:))
  print *
  
  
! ***
! *** create a Gnuplot file that displays the grid
! ***

  if ( .not. read_external_mesh ) then
  print *
  print *,'Saving the grid in Gnuplot format...'

  open(unit=20,file='OUTPUT_FILES/gridfile.gnu',status='unknown')

! draw horizontal lines of the grid
  print *,'drawing horizontal lines of the grid'
  istepx = 1
  if(ngnod == 4) then
    istepz = 1
  else
    istepz = 2
  endif
  do ili=0,nz,istepz
    do icol=0,nx-istepx,istepx
      write(20,*) sngl(x(icol,ili)),sngl(z(icol,ili))
      write(20,*) sngl(x(icol+istepx,ili)),sngl(z(icol+istepx,ili))
      write(20,10)
    enddo
  enddo

! draw vertical lines of the grid
  print *,'drawing vertical lines of the grid'
  if(ngnod == 4) then
    istepx = 1
  else
    istepx = 2
  endif
  istepz = 1
  do icol=0,nx,istepx
    do ili=0,nz-istepz,istepz
      write(20,*) sngl(x(icol,ili)),sngl(z(icol,ili))
      write(20,*) sngl(x(icol,ili+istepz)),sngl(z(icol,ili+istepz))
      write(20,10)
    enddo
  enddo

 10 format('')

  close(20)

! create a Gnuplot script to display the grid
  open(unit=20,file='OUTPUT_FILES/plotgnu',status='unknown')
  write(20,*) '#set term X11'
  write(20,*) 'set term postscript landscape monochrome solid "Helvetica" 22'
  write(20,*) 'set output "grid.ps"'
  write(20,*) '#set xrange [',sngl(minval(x)),':',sngl(maxval(x)),']'
  write(20,*) '#set yrange [',sngl(minval(z)),':',sngl(maxval(z)),']'
! use same unit length on both X and Y axes
  write(20,*) 'set size ratio -1'
  write(20,*) 'plot "gridfile.gnu" title "Macrobloc mesh" w l'
  write(20,*) 'pause -1 "Hit any key..."'
  close(20)

  print *,'Grid saved in Gnuplot format...'
  print *
  end if
     

  !*****************************
  ! Partitionning
  !*****************************
  allocate(part(0:nelmnts-1))

! if ngnod == 9, we work on a subarray of elmnts, which represents the elements with for nodes only
! construction of the graph
  if ( ngnod == 9 ) then
     allocate(elmnts_bis(0:ESIZE*nelmnts-1))
     do i = 0, nelmnts-1
        elmnts_bis(i*esize:i*esize+esize-1) = elmnts(i*ngnod:i*ngnod+esize-1)
     end do
        
     call mesh2dual_ncommonnodes(nelmnts, (nxread+1)*(nzread+1), elmnts_bis, xadj, adjncy, nnodes_elmnts, nodes_elmnts,1)
     
  else
     call mesh2dual_ncommonnodes(nelmnts, nnodes, elmnts, xadj, adjncy, nnodes_elmnts, nodes_elmnts,1)
     
  end if
     
  nb_edges = xadj(nelmnts)
 
! giving weight to edges and vertices. Currently not used.
  call read_weights(nelmnts, vwgt, nb_edges, adjwgt)
     
  if ( nproc == 1 ) then
      part(:) = 0
  else

! partitioning
     select case (partitionning_method)
     case(1)
        do iproc = 0, nproc-2
           part(iproc*floor(real(nelmnts)/real(nproc)):(iproc+1)*floor(real(nelmnts)/real(nproc))-1) = iproc
        end do
        part(floor(real(nelmnts)/real(nproc))*(nproc-1):nelmnts-1) = nproc - 1
        
     case(2)
#ifdef USE_METIS
        call Part_metis(nelmnts, xadj, adjncy, vwgt, adjwgt, nproc, nb_edges, edgecut, part, metis_options)
#else
        print *, 'This version of SPECFEM was not compiled with support of METIS.'
        print *, 'Please recompile with -DUSE_METIS in order to enable use of METIS.'
        stop
#endif
        
     case(3)
#ifdef USE_SCOTCH
        call Part_scotch(nelmnts, xadj, adjncy, vwgt, adjwgt, nproc, nb_edges, edgecut, part, scotch_strategy)
#else
        print *, 'This version of SPECFEM was not compiled with support of SCOTCH.'
        print *, 'Please recompile with -DUSE_SCOTCH in order to enable use of SCOTCH.'
        stop
#endif
        
     end select
 
  end if
  
! beware of fluid solid edges : coupled elements are transfered to the same partition
  if ( ngnod == 9 ) then
     call acoustic_elastic_repartitioning (nelmnts, nnodes, elmnts_bis, nb_materials, cs, num_material, &
          nproc, part, nedges_coupled, edges_coupled)
  else
     call acoustic_elastic_repartitioning (nelmnts, nnodes, elmnts, nb_materials, cs, num_material, &
          nproc, part, nedges_coupled, edges_coupled)
  end if
  
! local number of each element for each partition
  call Construct_glob2loc_elmnts(nelmnts, part, nproc, glob2loc_elmnts)
  
  if ( ngnod == 9 ) then
     deallocate(nnodes_elmnts)
     deallocate(nodes_elmnts)
     allocate(nnodes_elmnts(0:nnodes-1))
     allocate(nodes_elmnts(0:nsize*nnodes-1))
     nnodes_elmnts(:) = 0
     nodes_elmnts(:) = 0
     do i = 0, ngnod*nelmnts-1
        nodes_elmnts(elmnts(i)*nsize+nnodes_elmnts(elmnts(i))) = i/ngnod
        nnodes_elmnts(elmnts(i)) = nnodes_elmnts(elmnts(i)) + 1
        
     end do
  end if
  
! local number of each node for each partition
  call Construct_glob2loc_nodes(nelmnts, nnodes, nnodes_elmnts, nodes_elmnts, part, nproc, &
       glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes)

! construct the interfaces between partitions (used for MPI assembly)
  if ( nproc /= 1 ) then
     if ( ngnod == 9 ) then
        call Construct_interfaces(nelmnts, nproc, part, elmnts_bis, xadj, adjncy, tab_interfaces, &
             tab_size_interfaces, ninterfaces, nb_materials, cs, num_material)
     else
        call Construct_interfaces(nelmnts, nproc, part, elmnts, xadj, adjncy, tab_interfaces, &
             tab_size_interfaces, ninterfaces, nb_materials, cs, num_material)
     end if
     print *, '04'
     allocate(my_interfaces(0:ninterfaces-1))
     allocate(my_nb_interfaces(0:ninterfaces-1))
     print *, '05'
  end if
  
! setting absorbing boundaries by elements instead of edges
  if ( any_abs ) then
     call merge_abs_boundaries(nelemabs, nelemabs_merge, abs_surface, abs_surface_char, abs_surface_merge, &
          ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
          jbegin_left,jend_left,jbegin_right,jend_right, &
          nedges_coupled, edges_coupled, nb_materials, cs, num_material, &
          nelmnts, &
          elmnts, ngnod)
     print *, 'nelemabs_merge', nelemabs_merge
  end if

! *** generate the databases for the solver

  do iproc = 0, nproc-1
     
     write(prname, "('/Database',i5.5)") iproc
     open(unit=15,file='./OUTPUT_FILES'//prname,status='unknown')

     write(15,*) '#'
     write(15,*) '# Database for SPECFEM2D'
     write(15,*) '# Dimitri Komatitsch, (c) University of Pau, France'
     write(15,*) '#'
  
     write(15,*) 'Title of the simulation'
     write(15,"(a50)") title
     

     call write_glob2loc_nodes_database(15, iproc, npgeo, nodes_coords, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
          glob2loc_nodes, nnodes, 1)
     

     call write_partition_database(15, iproc, nspec, nelmnts, elmnts, glob2loc_elmnts, glob2loc_nodes_nparts, &
          glob2loc_nodes_parts, glob2loc_nodes, part, num_material, ngnod, 1)
     

     write(15,*) 'npgeo'
     write(15,*) npgeo
     
     write(15,*) 'gnuplot interpol'
     write(15,*) gnuplot,interpol
     
     write(15,*) 'NTSTEP_BETWEEN_OUTPUT_INFO'
     write(15,*) NTSTEP_BETWEEN_OUTPUT_INFO
     
     write(15,*) 'output_postscript_snapshot output_color_image colors numbers'
     write(15,*) output_postscript_snapshot,output_color_image,' 1 0'

     write(15,*) 'meshvect modelvect boundvect cutsnaps subsamp sizemax_arrows'
     write(15,*) meshvect,modelvect,boundvect,cutsnaps,subsamp,sizemax_arrows
     
     write(15,*) 'anglerec'
     write(15,*) anglerec
     
     write(15,*) 'initialfield'
     write(15,*) initialfield

     write(15,*) 'seismotype imagetype'
     write(15,*) seismotype,imagetype

     write(15,*) 'assign_external_model outputgrid TURN_ANISOTROPY_ON TURN_ATTENUATION_ON'
     write(15,*) assign_external_model,outputgrid,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON
     
     write(15,*) 'nt deltat'
     write(15,*) nt,deltat

     write(15,*) 'source'
     write(15,*) source_type,time_function_type,xs,zs,f0,t0,factor,angleforce,Mxx,Mzz,Mxz
     
     write(15,*) 'Coordinates of macrobloc mesh (coorg):'
     

     call write_glob2loc_nodes_database(15, iproc, npgeo, nodes_coords, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
          glob2loc_nodes, nnodes, 2)
     

     write(15,*) 'numat ngnod nspec pointsdisp plot_lowerleft_corner_only'
     write(15,*) nb_materials,ngnod,nspec,pointsdisp,plot_lowerleft_corner_only

     
     if ( any_abs ) then
        call write_abs_merge_database(15, nelemabs_merge, nelemabs_loc, &
             abs_surface_char, abs_surface_merge, &
             ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
             jbegin_left,jend_left,jbegin_right,jend_right, &
             glob2loc_elmnts, part, iproc, 1)
     else
        nelemabs_loc = 0
     end if
          
     call Write_surface_database(15, nelem_acoustic_surface, acoustic_surface, nelem_acoustic_surface_loc, &
          iproc, glob2loc_elmnts, &
          glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes, part, 1)
     

     call write_fluidsolid_edges_database(15, nedges_coupled, nedges_coupled_loc, &
          edges_coupled, glob2loc_elmnts, part, iproc, 1)
     
     write(15,*) 'nelemabs nelem_acoustic_surface num_fluid_solid_edges'
     write(15,*) nelemabs_loc,nelem_acoustic_surface_loc,nedges_coupled_loc
     
     
     write(15,*) 'Material sets (num 1 rho vp vs 0 0) or (num 2 rho c11 c13 c33 c44)'
     do i=1,nb_materials
        write(15,*) i,icodemat(i),rho(i),cp(i),cs(i),aniso3(i),aniso4(i)
     enddo
  
     write(15,*) 'Arrays kmato and knods for each bloc:'

     
     call write_partition_database(15, iproc, nspec, nelmnts, elmnts, glob2loc_elmnts, glob2loc_nodes_nparts, &
          glob2loc_nodes_parts, glob2loc_nodes, part, num_material, ngnod, 2)
     
     if ( nproc /= 1 ) then 
        call Write_interfaces_database(15, tab_interfaces, tab_size_interfaces, nproc, iproc, ninterfaces, &
             my_ninterface, my_interfaces, my_nb_interfaces, glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
             glob2loc_nodes, 1)
        
        write(15,*) 'Interfaces:'
        write(15,*) my_ninterface, maxval(my_nb_interfaces)
        
        call Write_interfaces_database(15, tab_interfaces, tab_size_interfaces, nproc, iproc, ninterfaces, &
             my_ninterface, my_interfaces, my_nb_interfaces, glob2loc_elmnts, glob2loc_nodes_nparts, glob2loc_nodes_parts, &
             glob2loc_nodes, 2)

     else
        write(15,*) 'Interfaces:'
        write(15,*) 0, 0
        
     end if
     

     write(15,*) 'List of absorbing elements (bottom right top left):'
     if ( any_abs ) then
        call write_abs_merge_database(15, nelemabs_merge, nelemabs_loc, &
             abs_surface_char, abs_surface_merge, &
             ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
             jbegin_left,jend_left,jbegin_right,jend_right, &
             glob2loc_elmnts, part, iproc, 2)
     end if
     
     write(15,*) 'List of acoustic free-surface elements:'
     call Write_surface_database(15, nelem_acoustic_surface, acoustic_surface, nelem_acoustic_surface_loc, &
          iproc, glob2loc_elmnts, &
          glob2loc_nodes_nparts, glob2loc_nodes_parts, glob2loc_nodes, part, 2)
     

     write(15,*) 'List of acoustic elastic coupled edges:'
     call write_fluidsolid_edges_database(15, nedges_coupled, nedges_coupled_loc, &
          edges_coupled, glob2loc_elmnts, part, iproc, 2)
  end do
  
  
! print position of the source
  print *
  print *,'Position (x,z) of the source = ',xs,zs
  print *

!--- compute position of the receivers and write the STATIONS file

  if ( generate_STATIONS ) then
  print *
  print *,'writing the DATA/STATIONS file'
  print *

! total number of receivers in all the receiver lines
  nrec_total = sum(nrec)

  print *
  print *,'There are ',nrec_total,' receivers'

  print *
  print *,'Position (x,z) of the ',nrec_total,' receivers'
  print *

  open(unit=15,file='DATA/STATIONS',status='unknown')
  write(15,*) nrec_total

  irec_global_number = 0

! loop on all the receiver lines
  do ireceiverlines = 1,nreceiverlines

! loop on all the receivers of this receiver line
    do irec = 1,nrec(ireceiverlines)

! compute global receiver number
      irec_global_number = irec_global_number + 1

! compute coordinates of the receiver
      if(nrec(ireceiverlines) > 1) then
        xrec = xdeb(ireceiverlines) + dble(irec-1)*(xfin(ireceiverlines)-xdeb(ireceiverlines))/dble(nrec(ireceiverlines)-1)
        zrec = zdeb(ireceiverlines) + dble(irec-1)*(zfin(ireceiverlines)-zdeb(ireceiverlines))/dble(nrec(ireceiverlines)-1)
      else
        xrec = xdeb(ireceiverlines)
        zrec = zdeb(ireceiverlines)
      endif

! modify position of receiver if we must record exactly at the surface
      if(enreg_surf(ireceiverlines)) &
        zrec = value_spline(xrec,xinterface_top,zinterface_top,coefs_interface_top,npoints_interface_top)

! display position of the receiver
      print *,'Receiver ',irec_global_number,' = ',xrec,zrec

      write(15,"('S',i4.4,'    AA ',f20.7,1x,f20.7,'       0.0         0.0')") irec_global_number,xrec,zrec

    enddo
  enddo

  close(15)
  end if

  print *
  

  end program meshfem2D

! *******************
! meshing subroutines
! *******************


!--- global node number

  integer function num(i,j,nx)

    implicit none
    
    integer i,j,nx
    
    num = j*(nx+1) + i + 1
    
  end function num
  

 !---  global node number (when ngnod==4).
  integer function num_4(i,j,nx)

    implicit none
    
    integer i,j,nx
    
    num_4 = j*(nx+1) + i + 1
    
  end function num_4

  
 !---  global node number (when ngnod==9).
  integer function num_9(i,j,nx,nz)
    
    implicit none
    
    integer i,j,nx,nz
    
    
    if ( (mod(i,2) == 0) .and. (mod(j,2) == 0) ) then
       num_9 = j/2 * (nx+1) + i/2 + 1
    else 
       if ( mod(j,2) == 0 ) then
          num_9 = (nx+1)*(nz+1) + j/2 * nx + ceiling(real(i)/real(2))  
       else 
          num_9 = (nx+1)*(nz+1) + nx*(nz+1) + floor(real(j)/real(2))*(nx*2+1) + i + 1
          
       end if
    end if
    
    
  end function num_9

  


!--- spline to describe the interfaces

  double precision function value_spline(x,xinterface,zinterface,coefs_interface,npoints_interface)

  implicit none

  integer npoints_interface
  double precision x,xp
  double precision, dimension(npoints_interface) :: xinterface,zinterface,coefs_interface

  value_spline = 0.d0

  xp = x

! assign the value on the edge if point is outside the model
  if(xp < xinterface(1)) xp = xinterface(1)
  if(xp > xinterface(npoints_interface)) xp = xinterface(npoints_interface)

  call splint(xinterface,zinterface,coefs_interface,npoints_interface,xp,value_spline)

  end function value_spline

! --------------------------------------

! compute spline coefficients (Numerical Recipes)
! modified to use dynamic allocation

  subroutine spline(x,y,n,yp1,ypn,y2)

  implicit none

  integer n
  double precision, dimension(n) :: x,y,y2
  double precision, dimension(:), allocatable :: u
  double precision yp1,ypn

  integer i,k
  double precision sig,p,qn,un

  allocate(u(n))

  y2(1)=-0.5d0
  u(1)=(3.d0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)

  do i=2,n-1
    sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
    p=sig*y2(i-1)+2.d0
    y2(i)=(sig-1.d0)/p
    u(i)=(6.d0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
  enddo

  qn=0.5d0
  un=(3.d0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.d0)

  do k=n-1,1,-1
    y2(k)=y2(k)*y2(k+1)+u(k)
  enddo

  deallocate(u)

  end subroutine spline

! --------------

! evaluate spline (adapted from Numerical Recipes)

  subroutine splint(xa,ya,y2a,n,x,y)

  implicit none

  integer n
  double precision, dimension(n) :: XA,YA,Y2A
  double precision x,y

  integer k,klo,khi
  double precision h,a,b

  KLO = 1
  KHI = N

  do while (KHI-KLO > 1)
    K=(KHI+KLO)/2
    if(XA(K) > X) then
      KHI=K
    else
      KLO=K
    endif
  enddo

  H = XA(KHI) - XA(KLO)
  IF (H == 0.d0) stop 'bad input in spline evaluation'

  A = (XA(KHI)-X) / H
  B = (X-XA(KLO)) / H

  Y = A*YA(KLO) + B*YA(KHI) + ((A**3-A)*Y2A(KLO) + (B**3-B)*Y2A(KHI))*(H**2)/6.d0

  end subroutine splint


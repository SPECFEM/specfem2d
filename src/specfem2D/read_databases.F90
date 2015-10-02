
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
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
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

  subroutine read_databases_init()

! starts reading in parameters from input Database file

  use specfem_par

  implicit none

  ! local parameters

  !integer :: ier,nproc
  character(len=80) :: datlin
  character(len=256)  :: prname

  ! opens Database file
  write(prname,230) myrank
  open(unit=IIN,file=prname,status='old',action='read',iostat=ier)
  if( ier /= 0 ) call exit_MPI('error opening file OUTPUT/Database***')

  !---  read job title and skip remaining titles of the input file
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a50)") simulation_title

  !---- print the date, time and start-up banner
  if (myrank == 0) call datim(simulation_title)

  if (myrank == 0) then
    write(IOUT,*)
    write(IOUT,*)
    write(IOUT,*) '*********************'
    write(IOUT,*) '****             ****'
    write(IOUT,*) '****  SPECFEM2D  ****'
    write(IOUT,*) '****             ****'
    write(IOUT,*) '*********************'
  endif

  !---- read parameters from input file
  read(IIN,"(a80)") datlin
  read(IIN,*) AXISYM
  if (myrank == 0 .and. AXISYM) then
    write(IOUT,*)
    write(IOUT,*)
    write(IOUT,*) '====================================================='
    write(IOUT,*) '=== A x i s y m m e t r i c   S i m u l a t i o n ==='
    write(IOUT,*) '====================================================='
  endif

  read(IIN,"(a80)") datlin
  read(IIN,*) SIMULATION_TYPE, NOISE_TOMOGRAPHY, SAVE_FORWARD, UNDO_ATTENUATION

  read(IIN,"(a80)") datlin
  read(IIN,*) npgeo,nproc

  read(IIN,"(a80)") datlin
  read(IIN,*) output_grid_Gnuplot,interpol

  read(IIN,"(a80)") datlin
  read(IIN,*) NSTEP_BETWEEN_OUTPUT_INFO

  read(IIN,"(a80)") datlin
  read(IIN,*) NSTEP_BETWEEN_OUTPUT_SEISMOS

  read(IIN,"(a80)") datlin
  read(IIN,*) NSTEP_BETWEEN_OUTPUT_IMAGES

  read(IIN,"(a80)") datlin
  read(IIN,*) PML_BOUNDARY_CONDITIONS

  read(IIN,"(a80)") datlin
  read(IIN,*) ROTATE_PML_ACTIVATE

  read(IIN,"(a80)") datlin
  read(IIN,*) ROTATE_PML_ANGLE

  read(IIN,"(a80)") datlin
  read(IIN,*) read_external_mesh

  read(IIN,"(a80)") datlin
  read(IIN,*) NELEM_PML_THICKNESS

  read(IIN,"(a80)") datlin
  read(IIN,*) NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS

  read(IIN,"(a80)") datlin
  read(IIN,*) subsamp_seismos,imagetype_JPEG,imagetype_wavefield_dumps

  read(IIN,"(a80)") datlin
  read(IIN,*) output_postscript_snapshot,output_color_image,colors,numbers

  read(IIN,"(a80)") datlin
  read(IIN,*) meshvect,modelvect,boundvect,cutsnaps,subsamp_postscript,sizemax_arrows
  cutsnaps = cutsnaps / 100.d0

  read(IIN,"(a80)") datlin
  read(IIN,*) anglerec

  read(IIN,"(a80)") datlin
  read(IIN,*) initialfield,add_Bielak_conditions
  if(add_Bielak_conditions .and. .not. initialfield) &
    stop 'need to have an initial field to add Bielak plane wave conditions'

  read(IIN,"(a80)") datlin
  read(IIN,*) seismotype,imagetype_postscript
  if(seismotype < 1 .or. seismotype > 6) call exit_MPI('Wrong type for seismogram output')
  if(imagetype_postscript < 1 .or. imagetype_postscript > 4) call exit_MPI('Wrong type for PostScript snapshots')

  if(SAVE_FORWARD .and. (seismotype /= 1 .and. seismotype /= 6)) then
    print *, '***** WARNING *****'
    print *, 'seismotype =',seismotype
    print *, 'Save forward wavefield => seismogram must be in displacement for (poro)elastic or potential for acoustic'
    print *, 'Seismotype must be changed to 1 (elastic/poroelastic adjoint source) or 6 (acoustic adjoint source)'
  !  stop
  endif

  read(IIN,"(a80)") datlin
  read(IIN,'(a100)') MODEL

  read(IIN,"(a80)") datlin
  read(IIN,'(a100)') SAVE_MODEL

  read(IIN,"(a80)") datlin
  read(IIN,'(a100)') TOMOGRAPHY_FILE

  read(IIN,"(a80)") datlin
  read(IIN,*) output_grid_ASCII,output_energy,output_wavefield_dumps

  read(IIN,"(a80)") datlin
  read(IIN,*) use_binary_for_wavefield_dumps

  read(IIN,"(a80)") datlin
  read(IIN,*) ATTENUATION_VISCOELASTIC_SOLID,ATTENUATION_PORO_FLUID_PART

  read(IIN,"(a80)") datlin
  read(IIN,*) save_ASCII_seismograms

  read(IIN,"(a80)") datlin
  read(IIN,*) save_binary_seismograms_single,save_binary_seismograms_double

  read(IIN,"(a80)") datlin
  read(IIN,*) USE_TRICK_FOR_BETTER_PRESSURE

  read(IIN,"(a80)") datlin
  read(IIN,*) save_ASCII_kernels

  read(IIN,"(a80)") datlin
  read(IIN,*) DRAW_SOURCES_AND_RECEIVERS

  read(IIN,"(a80)") datlin
  read(IIN,*) Q0,freq0

  read(IIN,"(a80)") datlin
  read(IIN,*) p_sv

  read(IIN,"(a80)") datlin
  read(IIN,*) factor_subsample_image

  read(IIN,"(a80)") datlin
  read(IIN,*) USE_CONSTANT_MAX_AMPLITUDE

  read(IIN,"(a80)") datlin
  read(IIN,*) CONSTANT_MAX_AMPLITUDE_TO_USE

  read(IIN,"(a80)") datlin
  read(IIN,*) USE_SNAPSHOT_NUMBER_IN_FILENAME

  read(IIN,"(a80)") datlin
  read(IIN,*) DRAW_WATER_IN_BLUE

  read(IIN,"(a80)") datlin
  read(IIN,*) US_LETTER

  read(IIN,"(a80)") datlin
  read(IIN,*) POWER_DISPLAY_COLOR

  read(IIN,"(a80)") datlin
  read(IIN,*) SU_FORMAT

  read(IIN,"(a80)") datlin
  read(IIN,*) USER_T0

  read(IIN,"(a80)") datlin
  read(IIN,*) time_stepping_scheme

  read(IIN,"(a80)") datlin
  read(IIN,*) ADD_SPRING_TO_STACEY

  read(IIN,"(a80)") datlin
  read(IIN,*) ADD_PERIODIC_CONDITIONS

  read(IIN,"(a80)") datlin
  read(IIN,*) PERIODIC_HORIZ_DIST

  read(IIN,"(a80)") datlin
  read(IIN,*) GPU_MODE

  !---- check parameters read
  if (myrank == 0) then
    write(IOUT,200) npgeo,NDIM
    write(IOUT,600) NSTEP_BETWEEN_OUTPUT_INFO,colors,numbers
    write(IOUT,700) seismotype,anglerec
    write(IOUT,750) initialfield,add_Bielak_conditions,&
                    ATTENUATION_VISCOELASTIC_SOLID, &
                    output_grid_ASCII,output_energy
    write(IOUT,800) imagetype_postscript,100.d0*cutsnaps,subsamp_postscript
  endif


  !---- check model
  if (trim(MODEL) == 'default') then
      assign_external_model = .false.
  else
      assign_external_model = .true.
  endif

  !---- read time step
  read(IIN,"(a80)") datlin
  read(IIN,*) NSTEP,deltat

  read(IIN,"(a80)") datlin
  read(IIN,*) NT_DUMP_ATTENUATION
  if (myrank == 0) write(IOUT,703) NSTEP,deltat,NSTEP*deltat

  if(SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. &
     (ATTENUATION_PORO_FLUID_PART)) then
    print *, '*************** error ***************'
    stop 'Anisotropy & Viscous damping are not presently implemented for adjoint calculations'
  endif

  if(SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. &
     ATTENUATION_PORO_FLUID_PART .and. (.not. UNDO_ATTENUATION)) then
    print *, '*************** error ***************'
    stop 'attenuation is only implemented for adjoint calculations with UNDO_ATTENUATION'
  endif

! make sure NSTEP_BETWEEN_OUTPUT_SEISMOS is a multiple of subsamp_seismos
  if(mod(NSTEP_BETWEEN_OUTPUT_SEISMOS,subsamp_seismos) /= 0) &
    stop 'error: NSTEP_BETWEEN_OUTPUT_SEISMOS must be a multiple of subsamp_seismos'

! make sure NSTEP is a multiple of subsamp_seismos
! if not, increase it a little bit, to the next multiple
  if(mod(NSTEP,subsamp_seismos) /= 0) then
    if (myrank == 0) then
      print *,'NSTEP is not a multiple of subsamp_seismos'
      NSTEP = (NSTEP/subsamp_seismos + 1)*subsamp_seismos
      print *,'thus increasing it automatically to the next multiple, which is ',NSTEP
      print *
    endif
  endif

  nproc_read_from_database=nproc

! output seismograms at least once at the end of the simulation
  NSTEP_BETWEEN_OUTPUT_SEISMOS = min(NSTEP,NSTEP_BETWEEN_OUTPUT_SEISMOS)

! read the ACOUSTIC_FORCING flag
  read(IIN,"(a80)") datlin
  read(IIN,*) ACOUSTIC_FORCING

  !----  read source information
  read(IIN,"(a80)") datlin
  read(IIN,*) NSOURCES

  ! output formats
230 format('./OUTPUT_FILES/Database',i5.5)

200 format(//1x,'C o n t r o l',/1x,13('='),//5x,&
  'Number of spectral element control nodes. . .(npgeo) =',i8/5x, &
  'Number of space dimensions. . . . . . . . . . (NDIM) =',i8)

600 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Display frequency . . . .(NSTEP_BETWEEN_OUTPUT_INFO) = ',i6/ 5x, &
  'Color display . . . . . . . . . . . . . . . (colors) = ',i6/ 5x, &
  '        ==  0     black and white display              ',  / 5x, &
  '        ==  1     color display                        ',  /5x, &
  'Numbered mesh . . . . . . . . . . . . . . .(numbers) = ',i6/ 5x, &
  '        ==  0     do not number the mesh               ',  /5x, &
  '        ==  1     number the mesh                      ')

700 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Seismograms recording type . . . . . . .(seismotype) = ',i6/5x, &
  'Angle for first line of receivers. . . . .(anglerec) = ',f6.2)

750 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Read external initial field. . . . . .(initialfield) = ',l6/5x, &
  'Add Bielak conditions . . . .(add_Bielak_conditions) = ',l6/5x, &
  'Attenuation on/off .(ATTENUATION_VISCOELASTIC_SOLID) = ',l6/5x, &
  'Save grid in ASCII file or not . (output_grid_ASCII) = ',l6/5x, &
  'Save a file with total energy or not.(output_energy) = ',l6)

800 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Vector display type . . . . . . . . .(imagetype_postscript) = ',i6/5x, &
  'Percentage of cut for vector plots. . . . . . . .(cutsnaps) = ',f6.2/5x, &
  'Subsampling of velocity model display. (subsamp_postscript) = ',i6)

703 format(//' I t e r a t i o n s '/1x,19('='),//5x, &
      'Number of time iterations . . . . .(NSTEP) =',i8,/5x, &
      'Time step increment. . . . . . . .(deltat) =',1pe15.6,/5x, &
      'Total simulation duration . . . . . (ttot) =',1pe15.6)

  end subroutine read_databases_init

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_sources()

! reads source parameters
  use specfem_par, only : NSOURCES,source_type,time_function_type,name_of_source_file,burst_band_width, &
                          x_source,z_source,Mxx,Mzz,Mxz,f0,tshift_src,factor,anglesource
  implicit none
  include "constants.h"


  ! local parameters
  integer :: i_source
  character(len=80) :: datlin

  ! initializes
  source_type(:) = 0
  time_function_type(:) = 0
  x_source(:) = 0.d0
  z_source(:) = 0.d0
  Mxx(:) = 0.d0
  Mzz(:) = 0.d0
  Mxz(:) = 0.d0
  f0(:) = 0.d0
  tshift_src(:) = 0.d0
  factor(:) = 0.d0
  anglesource(:) = 0.d0

  ! reads in source info from Database file
  do i_source=1,NSOURCES
     read(IIN,"(a80)") datlin
     read(IIN,*) source_type(i_source),time_function_type(i_source)
     read(IIN,"(a100)") name_of_source_file(i_source)
     read(IIN,*) burst_band_width(i_source),x_source(i_source),z_source(i_source),f0(i_source), &
                 tshift_src(i_source),factor(i_source),anglesource(i_source),Mxx(i_source),Mzz(i_source),Mxz(i_source)
  enddo

  end subroutine read_databases_sources

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_atten()

  use specfem_par, only : N_SLS,f0_attenuation,READ_VELOCITIES_AT_f0

! reads attenuation information

  implicit none
  include "constants.h"

  ! local parameters
  character(len=80) :: datlin

  read(IIN,"(a80)") datlin
  read(IIN,*) N_SLS, f0_attenuation, READ_VELOCITIES_AT_f0
  if(N_SLS < 2) stop 'must have N_SLS >= 2 even if attenuation if off because it is used to assign some arrays'

  end subroutine read_databases_atten

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_coorg_elem()

! reads the spectral macrobloc nodal coordinates

  use specfem_par, only : myrank,npgeo,coorg,numat,ngnod,nspec, &
                          pointsdisp,plot_lowerleft_corner_only, &
                          nelemabs,nelem_acforcing,nelem_acoustic_surface, &
                          num_fluid_solid_edges,num_fluid_poro_edges, &
                          num_solid_poro_edges,nnodes_tangential_curve, &
                          nelem_on_the_axis


  implicit none
  include "constants.h"

  ! local parameters
  integer :: ipoin,ip,id
  double precision, dimension(:), allocatable :: coorgread
  character(len=80) :: datlin

  ! initializes
  coorg(:,:) = 0.d0

  ! reads the spectral macrobloc nodal coordinates
  read(IIN,"(a80)") datlin

  ! reads in values
  ipoin = 0
  allocate(coorgread(NDIM))
  do ip = 1,npgeo
    ! reads coordinates
    read(IIN,*) ipoin,(coorgread(id),id =1,NDIM)

    if(ipoin<1 .or. ipoin>npgeo) call exit_MPI('Wrong control point number')

    ! saves coordinate array
    coorg(:,ipoin) = coorgread

  enddo
  deallocate(coorgread)

  !---- read the basic properties of the spectral elements
  read(IIN,"(a80)") datlin
  read(IIN,*) numat,ngnod,nspec,pointsdisp,plot_lowerleft_corner_only

  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,*) nelemabs,nelem_acforcing,nelem_acoustic_surface,num_fluid_solid_edges, &
              num_fluid_poro_edges,num_solid_poro_edges,nnodes_tangential_curve, &
              nelem_on_the_axis
  !---- print element group main parameters
  if (myrank == 0 ) then
    write(IOUT,107)
    write(IOUT,207) nspec,ngnod,NGLLX,NGLLZ,NGLLX*NGLLZ,pointsdisp,numat,nelem_acforcing
  endif

  ! output formats
107 format(/5x,'--> Spectral Elements (for mesh slice 0 only if using MPI runs) <--',//)

207 format(5x,'Number of spectral elements . . . . . . . . .  (nspec) =',i7,/5x, &
               'Number of control nodes per element . . . . . (ngnod) =',i7,/5x, &
               'Number of points in X-direction . . . . . . . (NGLLX) =',i7,/5x, &
               'Number of points in Y-direction . . . . . . . (NGLLZ) =',i7,/5x, &
               'Number of points per element. . . . . . (NGLLX*NGLLZ) =',i7,/5x, &
               'Number of points for display . . . . . . (pointsdisp) =',i7,/5x, &
               'Number of element material sets . . . . . . . (numat) =',i7,/5x, &
               'Number of acoustic forcing elements (nelem_acforcing) =',i7)

  end subroutine read_databases_coorg_elem


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_mato()

! reads spectral macrobloc data

  use specfem_par, only : nspec,ngnod,kmato,knods,region_CPML
  implicit none
  include "constants.h"

  ! local parameters
  integer :: n,k,ispec,kmato_read,pml_read
  integer, dimension(:), allocatable :: knods_read
  character(len=80) :: datlin

  ! initializes
  kmato(:) = 0
  knods(:,:) = 0
  region_CPML(:) = 0

  ! reads spectral macrobloc data
  read(IIN,"(a80)") datlin

  ! reads in values
  allocate(knods_read(ngnod))
  n = 0
  do ispec = 1,nspec
    ! format: #element_id  #material_id #node_id1 #node_id2 #...
    read(IIN,*) n,kmato_read,(knods_read(k), k=1,ngnod),pml_read
      ! material association
      kmato(n) = kmato_read
      region_CPML(n) = pml_read
      ! element control node indices
      knods(:,n)= knods_read(:)
  enddo
  deallocate(knods_read)


  end subroutine read_databases_mato

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_ninterface()

! reads in interface dimensions

  use specfem_par, only : ninterface,max_interface_size

  implicit none
  include "constants.h"

  ! local parameters
  character(len=80) :: datlin

  read(IIN,"(a80)") datlin
  read(IIN,*) ninterface, max_interface_size

  end subroutine read_databases_ninterface

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_interfaces()

! reads in interfaces

  use specfem_par, only : ninterface,my_neighbours,my_nelmnts_neighbours,my_interfaces
  implicit none
  include "constants.h"

  ! local parameters
  integer :: num_interface,ie,my_interfaces_read

  ! initializes
  my_neighbours(:) = -1
  my_nelmnts_neighbours(:) = 0
  my_interfaces(:,:,:) = -1

  ! reads in interfaces
  do num_interface = 1, ninterface
    ! format: #process_interface_id  #number_of_elements_on_interface
    ! where
    !     process_interface_id = rank of (neighbor) process to share MPI interface with
    !     number_of_elements_on_interface = number of interface elements
    read(IIN,*) my_neighbours(num_interface), my_nelmnts_neighbours(num_interface)

    ! loops over interface elements
    do ie = 1, my_nelmnts_neighbours(num_interface)
      ! format: #(1)spectral_element_id  #(2)interface_type  #(3)node_id1  #(4)node_id2
      !
      ! interface types:
      !     1  -  corner point only
      !     2  -  element edge
      read(IIN,*) my_interfaces_read, my_interfaces(2,ie,num_interface), &
              my_interfaces(3,ie,num_interface), my_interfaces(4,ie,num_interface)

      my_interfaces(1,ie,num_interface) = my_interfaces_read

    enddo
  enddo

  end subroutine read_databases_interfaces


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_absorbing()

! reads in absorbing edges

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par, only : myrank,nelemabs,nspec,anyabs, &
                          ibegin_edge1,iend_edge1,ibegin_edge2,iend_edge2, &
                          ibegin_edge3,iend_edge3,ibegin_edge4,iend_edge4, &
                          numabs,codeabs,codeabs_corner,typeabs, &
                          nspec_left,nspec_right,nspec_bottom,nspec_top, &
                          ib_right,ib_left,ib_bottom,ib_top,PML_BOUNDARY_CONDITIONS

  implicit none
  include "constants.h"

  ! local parameters
  integer :: inum,inum_duplicate,numabsread,typeabsread
  logical :: codeabsread(4)
  character(len=80) :: datlin

  integer :: nelemabs_tot,nspec_left_tot,nspec_right_tot,nspec_bottom_tot,nspec_top_tot

#ifdef USE_MPI
  integer :: ier
#endif

  ! initializes
  codeabs(:,:) = .false.
  codeabs_corner(:,:) = .false.
  typeabs(:) = 0

  ibegin_edge1(:) = 0
  iend_edge1(:) = 0
  ibegin_edge3(:) = 0
  iend_edge3(:) = 0

  ibegin_edge4(:) = 0
  iend_edge4(:) = 0
  ibegin_edge2(:) = 0
  iend_edge2(:) = 0

  nspec_left = 0
  nspec_right = 0
  nspec_bottom = 0
  nspec_top = 0

  ib_right(:) = 0
  ib_left(:) = 0
  ib_bottom(:) = 0
  ib_top(:) = 0

  ! reads in absorbing edges
  read(IIN,"(a80)") datlin

  ! reads in values
  if( anyabs ) then

    ! reads absorbing boundaries
    do inum = 1,nelemabs

! beware here and below that external meshes (for instance coming from CUBIT or Gmsh)
! may have rotated elements and thus edge 1 may not correspond to the bottom,
! edge 2 may not correspond to the right, edge 3 may not correspond to the top,
! and edge 4 may not correspond to the left.
      read(IIN,*) numabsread,codeabsread(1),codeabsread(2),codeabsread(3),&
                  codeabsread(4), typeabsread, ibegin_edge1(inum), iend_edge1(inum), &
                  ibegin_edge2(inum), iend_edge2(inum), ibegin_edge3(inum), &
                  iend_edge3(inum), ibegin_edge4(inum), iend_edge4(inum)

      if(numabsread < 1 .or. numabsread > nspec) &
        call exit_MPI('Wrong absorbing element number')

      numabs(inum) = numabsread

      codeabs(IEDGE1,inum) = codeabsread(1)
      codeabs(IEDGE2,inum) = codeabsread(2)
      codeabs(IEDGE3,inum) = codeabsread(3)
      codeabs(IEDGE4,inum) = codeabsread(4)

      typeabs(inum) = typeabsread

! check that a single edge is defined for each element cited
! (since elements with two absorbing edges MUST be cited twice, each time with a different "typeabs()" code
      if(count(codeabs(:,inum) .eqv. .true.) /= 1) then
        print *,'error for absorbing element inum = ',inum
        stop 'must have one and only one absorbing edge per absorbing line cited'
      endif

    enddo

! detection of the corner element
    do inum = 1,nelemabs

      if (codeabs(IEDGE1,inum)) then
        do inum_duplicate = 1,nelemabs
           if (inum == inum_duplicate) then
             ! left for blank, since no operation is needed.
           else
             if (numabs(inum) == numabs(inum_duplicate)) then
                if (codeabs(IEDGE4,inum_duplicate)) then
                   codeabs_corner(1,inum) = .true.
                endif

                if (codeabs(IEDGE2,inum_duplicate)) then
                   codeabs_corner(2,inum) = .true.
                endif

             endif
           endif
        enddo
      endif

      if (codeabs(IEDGE3,inum)) then
        do inum_duplicate = 1,nelemabs
           if (inum == inum_duplicate) then
             ! left for blank, since no operation is needed.
           else
             if (numabs(inum) == numabs(inum_duplicate)) then
                if (codeabs(IEDGE4,inum_duplicate)) then
                   codeabs_corner(3,inum) = .true.
                endif

                if (codeabs(IEDGE2,inum_duplicate)) then
                   codeabs_corner(4,inum) = .true.
                endif

             endif
           endif
        enddo
      endif

    enddo

! detection of the corner element

    ! boundary element numbering
    do inum = 1,nelemabs

      if (codeabs(IEDGE1,inum)) then
        nspec_bottom = nspec_bottom + 1
        ib_bottom(inum) =  nspec_bottom

      else if (codeabs(IEDGE2,inum)) then
        nspec_right = nspec_right + 1
        ib_right(inum) =  nspec_right

      else if (codeabs(IEDGE3,inum)) then
        nspec_top = nspec_top + 1
        ib_top(inum) = nspec_top

      else if (codeabs(IEDGE4,inum)) then
        nspec_left = nspec_left + 1
        ib_left(inum) =  nspec_left

      else
        stop 'incorrect absorbing boundary element type read'
      endif
    enddo

  else ! if this MPI slice has no absorbing element at all

      nelemabs = 0
      nspec_left = 0
      nspec_right = 0
      nspec_bottom = 0
      nspec_top = 0

  endif

#ifdef USE_MPI
      call MPI_REDUCE(nelemabs, nelemabs_tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
      call MPI_REDUCE(nspec_left, nspec_left_tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
      call MPI_REDUCE(nspec_right, nspec_right_tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
      call MPI_REDUCE(nspec_bottom, nspec_bottom_tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
      call MPI_REDUCE(nspec_top, nspec_top_tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
#else
      nelemabs_tot = nelemabs
      nspec_left_tot = nspec_left
      nspec_right_tot = nspec_right
      nspec_bottom_tot = nspec_bottom
      nspec_top_tot = nspec_top
#endif

    if (myrank == 0 .and. .not. PML_BOUNDARY_CONDITIONS) then
      write(IOUT,*)
      write(IOUT,*) 'Number of absorbing elements: ',nelemabs_tot
      write(IOUT,*) '  nspec_left = ',nspec_left_tot
      write(IOUT,*) '  nspec_right = ',nspec_right_tot
      write(IOUT,*) '  nspec_bottom = ',nspec_bottom_tot
      write(IOUT,*) '  nspec_top = ',nspec_top_tot
      write(IOUT,*)
    endif

  end subroutine read_databases_absorbing

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_acoustic_forcing()

! reads in absorbing edges

  use specfem_par, only : myrank,nelem_acforcing,nspec,ACOUSTIC_FORCING, &
                          ibegin_edge1_acforcing,iend_edge1_acforcing,ibegin_edge2_acforcing,iend_edge2_acforcing, &
                          ibegin_edge3_acforcing,iend_edge3_acforcing,ibegin_edge4_acforcing,iend_edge4_acforcing, &
                          numacforcing,codeacforcing,typeacforcing, &
                          nspec_left_acforcing,nspec_right_acforcing,nspec_bottom_acforcing,nspec_top_acforcing, &
                          ib_right_acforcing,ib_left_acforcing,ib_bottom_acforcing,ib_top_acforcing

  implicit none
  include "constants.h"

  ! local parameters
  integer :: inum,numacforcingread,typeacforcingread
  logical :: codeacforcingread(4)
  character(len=80) :: datlin

  ! initializes
  codeacforcing(:,:) = .false.
  typeacforcing(:) = 0

  ibegin_edge1_acforcing(:) = 0
  iend_edge1_acforcing(:) = 0
  ibegin_edge3_acforcing(:) = 0
  iend_edge3_acforcing(:) = 0

  ibegin_edge4_acforcing(:) = 0
  iend_edge4_acforcing(:) = 0
  ibegin_edge2_acforcing(:) = 0
  iend_edge2_acforcing(:) = 0

  nspec_left_acforcing = 0
  nspec_right_acforcing = 0
  nspec_bottom_acforcing = 0
  nspec_top_acforcing = 0

  ib_right_acforcing(:) = 0
  ib_left_acforcing(:) = 0
  ib_bottom_acforcing(:) = 0
  ib_top_acforcing(:) = 0

  ! reads in absorbing edges
  read(IIN,"(a80)") datlin

  ! reads in values
  if( ACOUSTIC_FORCING ) then

    ! reads absorbing boundaries
    do inum = 1,nelem_acforcing

! beware here and below that external meshes (for instance coming from CUBIT or Gmsh)
! may have rotated elements and thus edge 1 may not correspond to the bottom,
! edge 2 may not correspond to the right, edge 3 may not correspond to the top,
! and edge 4 may not correspond to the left.
      read(IIN,*) numacforcingread,codeacforcingread(1),codeacforcingread(2),codeacforcingread(3),&
                  codeacforcingread(4), typeacforcingread, ibegin_edge1_acforcing(inum), iend_edge1_acforcing(inum), &
                  ibegin_edge2_acforcing(inum), iend_edge2_acforcing(inum), ibegin_edge3_acforcing(inum), &
                  iend_edge3_acforcing(inum), ibegin_edge4_acforcing(inum), iend_edge4_acforcing(inum)

      if(numacforcingread < 1 .or. numacforcingread > nspec) &
        call exit_MPI('Wrong absorbing element number')


        numacforcing(inum) = numacforcingread

      codeacforcing(IEDGE1,inum) = codeacforcingread(1)
      codeacforcing(IEDGE2,inum) = codeacforcingread(2)
      codeacforcing(IEDGE3,inum) = codeacforcingread(3)
      codeacforcing(IEDGE4,inum) = codeacforcingread(4)

      typeacforcing(inum) = typeacforcingread

! check that a single edge is defined for each element cited
! (since elements with two absorbing edges MUST be cited twice, each time with a different "typeacforcing()" code
      if(count(codeacforcing(:,inum) .eqv. .true.) /= 1) then
        print *,'error for absorbing element inum = ',inum
        stop 'must have one and only one absorbing edge per absorbing line cited'
      endif

    enddo

    ! boundary element numbering
    do inum = 1,nelem_acforcing
      if (typeacforcing(inum) == IBOTTOM) then
        nspec_bottom_acforcing = nspec_bottom_acforcing + 1
        ib_bottom_acforcing(inum) =  nspec_bottom_acforcing

      else if (typeacforcing(inum) == IRIGHT) then
        nspec_right_acforcing = nspec_right_acforcing + 1
        ib_right_acforcing(inum) =  nspec_right_acforcing

      else if (typeacforcing(inum) == ITOP) then
        nspec_top_acforcing = nspec_top_acforcing + 1
        ib_top_acforcing(inum) = nspec_top_acforcing

      else if (typeacforcing(inum) == ILEFT) then
        nspec_left_acforcing = nspec_left_acforcing + 1
        ib_left_acforcing(inum) =  nspec_left_acforcing

      else
        stop 'incorrect absorbing boundary element type read'
      endif
    enddo

    if (myrank == 0) then
      write(IOUT,*)
      write(IOUT,*) 'Number of acoustic forcing elements: ',nelem_acforcing
      write(IOUT,*) '  nspec_left_acforcing = ',nspec_left_acforcing
      write(IOUT,*) '  nspec_right_acforcing = ',nspec_right_acforcing
      write(IOUT,*) '  nspec_bottom_acforcing = ',nspec_bottom_acforcing
      write(IOUT,*) '  nspec_top_acforcing = ',nspec_top_acforcing
      write(IOUT,*)
    endif

  endif

  end subroutine read_databases_acoustic_forcing

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_free_surf()

! reads acoustic free surface data
  use specfem_par, only : nelem_acoustic_surface,acoustic_edges,any_acoustic_edges
  implicit none
  include "constants.h"

  ! local parameters
  integer :: inum,acoustic_edges_read
  character(len=80) :: datlin

  ! initializes
  acoustic_edges(:,:) = 0

  ! reads in any possible free surface edges
  read(IIN,"(a80)") datlin

  if( any_acoustic_edges ) then
    do inum = 1,nelem_acoustic_surface
      read(IIN,*) acoustic_edges_read, acoustic_edges(2,inum), acoustic_edges(3,inum), &
           acoustic_edges(4,inum)

      acoustic_edges(1,inum) = acoustic_edges_read

    enddo

  endif

  end subroutine read_databases_free_surf

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_coupled()

! reads acoustic elastic coupled edges
! reads acoustic poroelastic coupled edges
! reads poroelastic elastic coupled edges
  use specfem_par, only : num_fluid_solid_edges,any_fluid_solid_edges, &
                          fluid_solid_acoustic_ispec,fluid_solid_elastic_ispec, &
                          num_fluid_poro_edges,any_fluid_poro_edges, &
                          fluid_poro_acoustic_ispec,fluid_poro_poroelastic_ispec, &
                          num_solid_poro_edges,any_solid_poro_edges, &
                          solid_poro_elastic_ispec,solid_poro_poroelastic_ispec

  implicit none
  include "constants.h"

  ! local parameters
  integer :: inum
  integer :: fluid_solid_acoustic_ispec_read,fluid_solid_elastic_ispec_read, &
    fluid_poro_acoustic_ispec_read,fluid_poro_poro_ispec_read, &
    solid_poro_poro_ispec_read,solid_poro_elastic_ispec_read
  character(len=80) :: datlin

  ! initializes
  fluid_solid_acoustic_ispec(:) = 0
  fluid_solid_elastic_ispec(:) = 0
  fluid_poro_acoustic_ispec(:) = 0
  fluid_poro_poroelastic_ispec(:) = 0
  solid_poro_elastic_ispec(:) = 0
  solid_poro_poroelastic_ispec(:) = 0

  ! reads acoustic elastic coupled edges
  read(IIN,"(a80)") datlin

  if ( any_fluid_solid_edges ) then
    do inum = 1, num_fluid_solid_edges
      read(IIN,*) fluid_solid_acoustic_ispec_read,fluid_solid_elastic_ispec_read

      fluid_solid_acoustic_ispec(inum) = fluid_solid_acoustic_ispec_read
      fluid_solid_elastic_ispec(inum) = fluid_solid_elastic_ispec_read
    enddo
  endif

  ! reads acoustic poroelastic coupled edges
  read(IIN,"(a80)") datlin

  if ( any_fluid_poro_edges ) then
    do inum = 1, num_fluid_poro_edges
      read(IIN,*) fluid_poro_acoustic_ispec_read,fluid_poro_poro_ispec_read

      fluid_poro_acoustic_ispec(inum) = fluid_poro_acoustic_ispec_read
      fluid_poro_poroelastic_ispec(inum) = fluid_poro_poro_ispec_read
    enddo
  endif

  ! reads poroelastic elastic coupled edges
  read(IIN,"(a80)") datlin

  if ( any_solid_poro_edges ) then
    do inum = 1, num_solid_poro_edges
      read(IIN,*) solid_poro_poro_ispec_read,solid_poro_elastic_ispec_read

      solid_poro_elastic_ispec(inum) = solid_poro_elastic_ispec_read
      solid_poro_poroelastic_ispec(inum) = solid_poro_poro_ispec_read
    enddo
  endif

  end subroutine read_databases_coupled

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_tangential_detection_curve()

! reads tangential detection curve
  use specfem_par, only : nnodes_tangential_curve,nodes_tangential_curve, &
                          force_normal_to_surface,rec_normal_to_surface, &
                          any_tangential_curve

  implicit none
  include "constants.h"

  ! local parameters
  integer :: i
  character(len=80) :: datlin

  ! initializes
  nodes_tangential_curve(:,:) = 0.d0

  ! reads tangential detection curve
  read(IIN,"(a80)") datlin
  read(IIN,*) force_normal_to_surface,rec_normal_to_surface

  if( any_tangential_curve ) then
    do i = 1, nnodes_tangential_curve
      read(IIN,*) nodes_tangential_curve(1,i),nodes_tangential_curve(2,i)
    enddo
  else
    force_normal_to_surface = .false.
    rec_normal_to_surface = .false.
  endif

  end subroutine read_tangential_detection_curve

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_axial_elements()

! reads axial elements data
  use specfem_par, only : nelem_on_the_axis,ispec_of_axial_elements
  implicit none
  include "constants.h"

  ! local parameters
  integer :: i
  character(len=80) :: datlin

  ! initializes
  ispec_of_axial_elements(:) = 0

  ! reads in any possible axial elements
  read(IIN,"(a80)") datlin

  do i = 1,nelem_on_the_axis
    read(IIN,*) ispec_of_axial_elements(i)
  enddo

  end subroutine read_databases_axial_elements



!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, INRIA and University of Pau, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
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

  subroutine read_databases_init(myrank,ipass, &
                  simulation_title,SIMULATION_TYPE,NOISE_TOMOGRAPHY,SAVE_FORWARD,npgeo,nproc, &
                  output_grid_Gnuplot,interpol,NSTEP_BETWEEN_OUTPUT_INFO,NSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP_BETWEEN_OUTPUT_IMAGES, &
                  PML_BOUNDARY_CONDITIONS,NELEM_PML_THICKNESS, &
                  NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS,subsamp_seismos,imagetype_JPEG,imagetype_wavefield_dumps, &
                  output_postscript_snapshot,output_color_image,colors,numbers, &
                  meshvect,modelvect,boundvect,cutsnaps,subsamp_postscript,sizemax_arrows, &
                  anglerec,initialfield,add_Bielak_conditions, &
                  seismotype,imagetype_postscript,assign_external_model,READ_EXTERNAL_SEP_FILE, &
                  output_grid_ASCII,output_energy,output_wavefield_dumps,use_binary_for_wavefield_dumps, &
                  ATTENUATION_VISCOELASTIC_SOLID,ATTENUATION_PORO_FLUID_PART,save_ASCII_seismograms, &
                  save_binary_seismograms_single,save_binary_seismograms_double, &
                  DRAW_SOURCES_AND_RECEIVERS,Q0,freq0,p_sv,NSTEP,deltat,NSOURCES, &
                  factor_subsample_image,USE_SNAPSHOT_NUMBER_IN_FILENAME,DRAW_WATER_IN_BLUE,US_LETTER, &
                  POWER_DISPLAY_COLOR,PERFORM_CUTHILL_MCKEE,SU_FORMAT,USER_T0,time_stepping_scheme,&
                  ADD_SPRING_TO_STACEY,ADD_PERIODIC_CONDITIONS,PERIODIC_horiz_dist,PERIODIC_DETECT_TOL,&
                  read_external_mesh)

! starts reading in parameters from input Database file

  implicit none
  include "constants.h"

  integer :: myrank,ipass
  character(len=60) simulation_title
  integer :: SIMULATION_TYPE,NOISE_TOMOGRAPHY,npgeo,nproc
  integer :: colors,numbers,subsamp_postscript,seismotype,imagetype_postscript
  logical :: SAVE_FORWARD,output_grid_Gnuplot,interpol,output_postscript_snapshot, &
    output_color_image
  logical :: meshvect,modelvect,boundvect,initialfield,add_Bielak_conditions, &
    assign_external_model,READ_EXTERNAL_SEP_FILE, &
    output_grid_ASCII,output_energy,output_wavefield_dumps,p_sv,use_binary_for_wavefield_dumps
  logical :: ATTENUATION_VISCOELASTIC_SOLID,ATTENUATION_PORO_FLUID_PART,PML_BOUNDARY_CONDITIONS, &
             save_ASCII_seismograms,save_binary_seismograms_single,save_binary_seismograms_double,DRAW_SOURCES_AND_RECEIVERS

  double precision :: cutsnaps,sizemax_arrows,anglerec
  double precision :: Q0,freq0
  double precision :: deltat

  integer :: NSTEP,NSOURCES
  integer :: NSTEP_BETWEEN_OUTPUT_INFO,NSTEP_BETWEEN_OUTPUT_SEISMOS,NSTEP_BETWEEN_OUTPUT_IMAGES,NSTEP_BETWEEN_OUTPUT_WAVE_DUMPS, &
             subsamp_seismos,imagetype_JPEG,imagetype_wavefield_dumps,NELEM_PML_THICKNESS

! factor to subsample color images output by the code (useful for very large models)
  integer :: factor_subsample_image

! use snapshot number in the file name of JPG color snapshots instead of the time step
  logical :: USE_SNAPSHOT_NUMBER_IN_FILENAME

! display acoustic layers as constant blue, because they likely correspond to water in the case of ocean acoustics
! or in the case of offshore oil industry experiments.
! (if off, display them as greyscale, as for elastic or poroelastic elements)
  logical :: DRAW_WATER_IN_BLUE

! US letter paper or European A4
  logical :: US_LETTER

! non linear display to enhance small amplitudes in color images
  double precision :: POWER_DISPLAY_COLOR

! perform inverse Cuthill-McKee (1969) permutation for mesh numbering
  logical :: PERFORM_CUTHILL_MCKEE

! output seismograms in Seismic Unix format (adjoint traces will be read in the same format)
  logical :: SU_FORMAT

! use this t0 as earliest starting time rather than the automatically calculated one
! (must be positive and bigger than the automatically one to be effective;
!  simulation will start at t = - t0)
  double precision :: USER_T0

! value of time_stepping_scheme to decide which time scheme will be used
! # 1 = Newmark (2nd order), 2 = LDDRK4-6 (4th-order 6-stage low storage Runge-Kutta)
! 3 = classical 4th-order 4-stage Runge-Kutta
  integer :: time_stepping_scheme

!! DK DK for add spring to stacey absorbing boundary condition
  logical :: ADD_SPRING_TO_STACEY

!! DK DK for horizontal periodic conditions: detect common points between left and right edges
  logical :: ADD_PERIODIC_CONDITIONS

!! DK DK horizontal periodicity distance for periodic conditions
  double precision :: PERIODIC_horiz_dist

!! DK DK grid point detection tolerance for periodic conditions
  double precision :: PERIODIC_DETECT_TOL

!!! DK DK for CPML_element_file
   logical :: read_external_mesh

  ! local parameters
  integer :: ier
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
  if (myrank == 0 .and. ipass == 1) call datim(simulation_title)

  if (myrank == 0 .and. ipass == 1) then
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
  read(IIN,*) SIMULATION_TYPE, NOISE_TOMOGRAPHY, SAVE_FORWARD

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
    print*, '***** WARNING *****'
    print*, 'seismotype =',seismotype
    print*, 'Save forward wavefield => seismogram must be in displacement for (poro)elastic or potential for acoustic'
    print*, 'Seismotype must be changed to 1 (elastic/poroelastic adjoint source) or 6 (acoustic adjoint source)'
    stop
  endif

  read(IIN,"(a80)") datlin
  read(IIN,*) assign_external_model,READ_EXTERNAL_SEP_FILE

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
  read(IIN,*) DRAW_SOURCES_AND_RECEIVERS

  read(IIN,"(a80)") datlin
  read(IIN,*) Q0,freq0

  read(IIN,"(a80)") datlin
  read(IIN,*) p_sv

  read(IIN,"(a80)") datlin
  read(IIN,*) factor_subsample_image

  read(IIN,"(a80)") datlin
  read(IIN,*) USE_SNAPSHOT_NUMBER_IN_FILENAME

  read(IIN,"(a80)") datlin
  read(IIN,*) DRAW_WATER_IN_BLUE

  read(IIN,"(a80)") datlin
  read(IIN,*) US_LETTER

  read(IIN,"(a80)") datlin
  read(IIN,*) POWER_DISPLAY_COLOR

  read(IIN,"(a80)") datlin
  read(IIN,*) PERFORM_CUTHILL_MCKEE

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
  read(IIN,*) PERIODIC_horiz_dist

  read(IIN,"(a80)") datlin
  read(IIN,*) PERIODIC_DETECT_TOL

  !---- check parameters read
  if (myrank == 0 .and. ipass == 1) then
    write(IOUT,200) npgeo,NDIM
    write(IOUT,600) NSTEP_BETWEEN_OUTPUT_INFO,colors,numbers
    write(IOUT,700) seismotype,anglerec
    write(IOUT,750) initialfield,add_Bielak_conditions,assign_external_model,&
                    READ_EXTERNAL_SEP_FILE,ATTENUATION_VISCOELASTIC_SOLID, &
                    output_grid_ASCII,output_energy
    write(IOUT,800) imagetype_postscript,100.d0*cutsnaps,subsamp_postscript
  endif

  !---- read time step
  read(IIN,"(a80)") datlin
  read(IIN,*) NSTEP,deltat
  if (myrank == 0 .and. ipass == 1) write(IOUT,703) NSTEP,deltat,NSTEP*deltat

  if(SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. &
     (ATTENUATION_VISCOELASTIC_SOLID .or. ATTENUATION_PORO_FLUID_PART)) then
    print *, '*************** error ***************'
    stop 'Anisotropy & Attenuation & Viscous damping are not presently implemented for adjoint calculations'
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

! output seismograms at least once at the end of the simulation
  NSTEP_BETWEEN_OUTPUT_SEISMOS = min(NSTEP,NSTEP_BETWEEN_OUTPUT_SEISMOS)

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
  'Assign external model . . . .(assign_external_model) = ',l6/5x, &
  'Read external SEP file . . .(READ_EXTERNAL_SEP_FILE) = ',l6/5x, &
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

  subroutine read_databases_sources(NSOURCES,source_type,time_function_type, &
                      x_source,z_source,Mxx,Mzz,Mxz,f0,tshift_src,factor,anglesource)

! reads source parameters

  implicit none
  include "constants.h"

  integer :: NSOURCES
  integer, dimension(NSOURCES) :: source_type,time_function_type
  double precision, dimension(NSOURCES) :: x_source,z_source, &
    Mxx,Mzz,Mxz,f0,tshift_src,factor,anglesource

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
     read(IIN,*) source_type(i_source),time_function_type(i_source), &
                 x_source(i_source),z_source(i_source),f0(i_source),tshift_src(i_source), &
                 factor(i_source),anglesource(i_source),Mxx(i_source),Mzz(i_source),Mxz(i_source)
  enddo

  end subroutine read_databases_sources

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_atten(N_SLS,f0_attenuation)

! reads attenuation information

  implicit none
  include "constants.h"

  integer :: N_SLS
  double precision :: f0_attenuation

  ! local parameters
  character(len=80) :: datlin

  read(IIN,"(a80)") datlin
  read(IIN,*) N_SLS, f0_attenuation
  if(N_SLS < 1) stop 'must have N_SLS >= 1 even if attenuation if off because it is used to assign some arrays'

  end subroutine read_databases_atten

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_coorg_elem(myrank,ipass,npgeo,coorg,numat,ngnod,nspec, &
                              pointsdisp,plot_lowerleft_corner_only, &
                              nelemabs,nelem_acoustic_surface, &
                              num_fluid_solid_edges,num_fluid_poro_edges, &
                              num_solid_poro_edges,nnodes_tangential_curve)

! reads the spectral macrobloc nodal coordinates

  implicit none
  include "constants.h"

  integer :: myrank,ipass,npgeo
  double precision, dimension(NDIM,npgeo) :: coorg

  integer :: numat,ngnod,nspec
  integer :: pointsdisp
  logical :: plot_lowerleft_corner_only
  integer :: nelemabs,nelem_acoustic_surface, &
    num_fluid_solid_edges,num_fluid_poro_edges, &
    num_solid_poro_edges,nnodes_tangential_curve

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
  read(IIN,*) nelemabs,nelem_acoustic_surface,num_fluid_solid_edges,num_fluid_poro_edges,&
              num_solid_poro_edges,nnodes_tangential_curve

  !---- print element group main parameters
  if (myrank == 0 .and. ipass == 1) then
    write(IOUT,107)
    write(IOUT,207) nspec,ngnod,NGLLX,NGLLZ,NGLLX*NGLLZ,pointsdisp,numat,nelemabs
  endif

  ! output formats
107 format(/5x,'--> Spectral Elements <--',//)

207 format(5x,'Number of spectral elements . . . . .  (nspec) =',i7,/5x, &
               'Number of control nodes per element .  (ngnod) =',i7,/5x, &
               'Number of points in X-direction . . .  (NGLLX) =',i7,/5x, &
               'Number of points in Y-direction . . .  (NGLLZ) =',i7,/5x, &
               'Number of points per element. . .(NGLLX*NGLLZ) =',i7,/5x, &
               'Number of points for display . . .(pointsdisp) =',i7,/5x, &
               'Number of element material sets . . .  (numat) =',i7,/5x, &
               'Number of absorbing elements . . . .(nelemabs) =',i7)

  end subroutine read_databases_coorg_elem


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_mato(ipass,nspec,ngnod,kmato,knods, &
                                perm,antecedent_list,region_CPML)

! reads spectral macrobloc data

  implicit none
  include "constants.h"

  integer :: ipass,ngnod,nspec
  integer, dimension(nspec) :: kmato
  integer, dimension(nspec) :: region_CPML
  integer, dimension(ngnod,nspec) :: knods

  integer, dimension(nspec) :: perm,antecedent_list

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
    if(ipass == 1) then
      ! material association
      kmato(n) = kmato_read
      region_CPML(n) = pml_read
      ! element control node indices
      knods(:,n)= knods_read(:)
    else if(ipass == 2) then
      kmato(perm(antecedent_list(n))) = kmato_read
      region_CPML(perm(antecedent_list(n))) = pml_read
      knods(:,perm(antecedent_list(n)))= knods_read(:)
    else
      call exit_MPI('error: maximum is 2 passes')
    endif
  enddo
  deallocate(knods_read)


  end subroutine read_databases_mato

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_ninterface(ninterface,max_interface_size)

! reads in interface dimensions

  implicit none
  include "constants.h"

  integer :: ninterface,max_interface_size

  ! local parameters
  character(len=80) :: datlin

  read(IIN,"(a80)") datlin
  read(IIN,*) ninterface, max_interface_size

  end subroutine read_databases_ninterface

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_interfaces(ipass,ninterface,nspec,max_interface_size, &
                              my_neighbours,my_nelmnts_neighbours,my_interfaces, &
                              perm,antecedent_list)

! reads in interfaces

  implicit none
  include "constants.h"

  integer :: ipass,nspec
  integer :: ninterface,max_interface_size
  integer, dimension(ninterface) :: my_neighbours,my_nelmnts_neighbours
  integer, dimension(4,max_interface_size,ninterface) :: my_interfaces

  integer, dimension(nspec) :: perm,antecedent_list

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

      if(ipass == 1) then
        my_interfaces(1,ie,num_interface) = my_interfaces_read
      else if(ipass == 2) then
        my_interfaces(1,ie,num_interface) = perm(antecedent_list(my_interfaces_read))
      else
        call exit_MPI('error: maximum number of passes is 2')
      endif

    enddo
  enddo

  end subroutine read_databases_interfaces


!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_absorbing(myrank,ipass,nelemabs,nspec,anyabs, &
                            ibegin_edge1,iend_edge1,ibegin_edge2,iend_edge2, &
                            ibegin_edge3,iend_edge3,ibegin_edge4,iend_edge4, &
                            numabs,codeabs,typeabs,perm,antecedent_list, &
                            nspec_left,nspec_right,nspec_bottom,nspec_top, &
                            ib_right,ib_left,ib_bottom,ib_top)

! reads in absorbing edges

  implicit none
  include "constants.h"

  integer :: myrank,ipass,nspec
  integer :: nelemabs
  integer, dimension(nelemabs) :: numabs,ibegin_edge1,iend_edge1, &
    ibegin_edge3,iend_edge3,ibegin_edge4,iend_edge4,ibegin_edge2,iend_edge2
  logical, dimension(4,nelemabs) :: codeabs
  integer, dimension(nelemabs) :: typeabs
  logical :: anyabs
  integer, dimension(nspec) :: perm,antecedent_list
  integer :: nspec_left,nspec_right,nspec_bottom,nspec_top

  integer, dimension(nelemabs) :: ib_right,ib_left,ib_bottom,ib_top

  ! local parameters
  integer :: inum,numabsread,typeabsread
  logical :: codeabsread(4)
  character(len=80) :: datlin

  ! initializes
  codeabs(:,:) = .false.
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

      if(ipass == 1) then
        numabs(inum) = numabsread
      else if(ipass == 2) then
        numabs(inum) = perm(antecedent_list(numabsread))
      else
        call exit_MPI('error: maximum number of passes is 2')
      endif

      codeabs(IEDGE1,inum) = codeabsread(1)
      codeabs(IEDGE2,inum) = codeabsread(2)
      codeabs(IEDGE3,inum) = codeabsread(3)
      codeabs(IEDGE4,inum) = codeabsread(4)

      typeabs(inum) = typeabsread

! check that a single edge is defined for each element cited
! (since elements with two absorbing edges MUST be cited twice, each time with a different "typeabs()" code
      if(count(codeabs(:,inum) .eqv. .true.) /= 1) stop 'must have one and only one absorbing edge per absorbing line cited'

    enddo

    ! boundary element numbering
    do inum = 1,nelemabs
      if (typeabs(inum) == IBOTTOM) then
        nspec_bottom = nspec_bottom + 1
        ib_bottom(inum) =  nspec_bottom

      else if (typeabs(inum) == IRIGHT) then
        nspec_right = nspec_right + 1
        ib_right(inum) =  nspec_right

      else if (typeabs(inum) == ITOP) then
        nspec_top = nspec_top + 1
        ib_top(inum) = nspec_top

      else if (typeabs(inum) == ILEFT) then
        nspec_left = nspec_left + 1
        ib_left(inum) =  nspec_left

      else
        stop 'incorrect absorbing boundary element type read'
      endif
    enddo

    if (myrank == 0 .and. ipass == 1) then
      write(IOUT,*)
      write(IOUT,*) 'Number of absorbing elements: ',nelemabs
      write(IOUT,*) '  nspec_left = ',nspec_left
      write(IOUT,*) '  nspec_right = ',nspec_right
      write(IOUT,*) '  nspec_bottom = ',nspec_bottom
      write(IOUT,*) '  nspec_top = ',nspec_top
      write(IOUT,*)
    endif

  endif

  end subroutine read_databases_absorbing

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_free_surf(ipass,nelem_acoustic_surface,nspec, &
                            acoustic_edges,perm,antecedent_list,any_acoustic_edges)

! reads acoustic free surface data

  implicit none
  include "constants.h"

  integer :: ipass,nspec
  integer :: nelem_acoustic_surface
  integer, dimension(4,nelem_acoustic_surface) :: acoustic_edges
  logical :: any_acoustic_edges

  integer, dimension(nspec) :: perm,antecedent_list

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

      if(ipass == 1) then
        acoustic_edges(1,inum) = acoustic_edges_read
      else if(ipass == 2) then
        acoustic_edges(1,inum) = perm(antecedent_list(acoustic_edges_read))
      else
        call exit_MPI('error: maximum number of passes is 2')
      endif

    enddo

  endif

  end subroutine read_databases_free_surf

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_coupled(ipass,nspec,num_fluid_solid_edges,any_fluid_solid_edges, &
                            fluid_solid_acoustic_ispec,fluid_solid_elastic_ispec, &
                            num_fluid_poro_edges,any_fluid_poro_edges, &
                            fluid_poro_acoustic_ispec,fluid_poro_poroelastic_ispec, &
                            num_solid_poro_edges,any_solid_poro_edges, &
                            solid_poro_elastic_ispec,solid_poro_poroelastic_ispec, &
                            perm,antecedent_list)

! reads acoustic elastic coupled edges
! reads acoustic poroelastic coupled edges
! reads poroelastic elastic coupled edges

  implicit none
  include "constants.h"

  integer :: ipass,nspec

  integer :: num_fluid_solid_edges
  logical :: any_fluid_solid_edges
  integer, dimension(num_fluid_solid_edges) :: fluid_solid_acoustic_ispec,fluid_solid_elastic_ispec

  integer :: num_fluid_poro_edges
  logical :: any_fluid_poro_edges
  integer, dimension(num_fluid_poro_edges) :: fluid_poro_acoustic_ispec,fluid_poro_poroelastic_ispec

  integer :: num_solid_poro_edges
  logical :: any_solid_poro_edges
  integer, dimension(num_solid_poro_edges) :: solid_poro_elastic_ispec,solid_poro_poroelastic_ispec

  integer, dimension(nspec) :: perm,antecedent_list

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

      if(ipass == 1) then
        fluid_solid_acoustic_ispec(inum) = fluid_solid_acoustic_ispec_read
        fluid_solid_elastic_ispec(inum) = fluid_solid_elastic_ispec_read
      else if(ipass == 2) then
        fluid_solid_acoustic_ispec(inum) = perm(antecedent_list(fluid_solid_acoustic_ispec_read))
        fluid_solid_elastic_ispec(inum) = perm(antecedent_list(fluid_solid_elastic_ispec_read))
      else
        call exit_MPI('error: maximum number of passes is 2')
      endif
    enddo
  endif

  ! reads acoustic poroelastic coupled edges
  read(IIN,"(a80)") datlin

  if ( any_fluid_poro_edges ) then
    do inum = 1, num_fluid_poro_edges
      read(IIN,*) fluid_poro_acoustic_ispec_read,fluid_poro_poro_ispec_read

      if(ipass == 1) then
        fluid_poro_acoustic_ispec(inum) = fluid_poro_acoustic_ispec_read
        fluid_poro_poroelastic_ispec(inum) = fluid_poro_poro_ispec_read
      else if(ipass == 2) then
        fluid_poro_acoustic_ispec(inum) = perm(antecedent_list(fluid_poro_acoustic_ispec_read))
        fluid_poro_poroelastic_ispec(inum) = perm(antecedent_list(fluid_poro_poro_ispec_read))
      else
        call exit_MPI('error: maximum number of passes is 2')
      endif
    enddo
  endif

  ! reads poroelastic elastic coupled edges
  read(IIN,"(a80)") datlin

  if ( any_solid_poro_edges ) then
    do inum = 1, num_solid_poro_edges
      read(IIN,*) solid_poro_poro_ispec_read,solid_poro_elastic_ispec_read

      if(ipass == 1) then
        solid_poro_elastic_ispec(inum) = solid_poro_elastic_ispec_read
        solid_poro_poroelastic_ispec(inum) = solid_poro_poro_ispec_read
      else if(ipass == 2) then
        solid_poro_elastic_ispec(inum) = perm(antecedent_list(solid_poro_elastic_ispec_read))
        solid_poro_poroelastic_ispec(inum) = perm(antecedent_list(solid_poro_poro_ispec_read))
      else
        call exit_MPI('error: maximum number of passes is 2')
      endif
    enddo
  endif

  end subroutine read_databases_coupled

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_final(nnodes_tangential_curve,nodes_tangential_curve, &
                                force_normal_to_surface,rec_normal_to_surface, &
                                any_tangential_curve )

! reads tangential detection curve
! and closes Database file

  implicit none
  include "constants.h"

  integer :: nnodes_tangential_curve
  logical :: any_tangential_curve
  double precision, dimension(2,nnodes_tangential_curve) :: nodes_tangential_curve

  logical :: force_normal_to_surface,rec_normal_to_surface

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

  ! closes input Database file
  close(IIN)

  end subroutine read_databases_final



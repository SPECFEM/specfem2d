
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
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
                  simulation_title,SIMULATION_TYPE,SAVE_FORWARD,npgeo, &
                  gnuplot,interpol,NTSTEP_BETWEEN_OUTPUT_INFO, &
                  output_postscript_snapshot,output_color_image,colors,numbers, &
                  meshvect,modelvect,boundvect,cutsnaps,subsamp,sizemax_arrows, &
                  anglerec,initialfield,add_Bielak_conditions, &
                  seismotype,imagetype,assign_external_model,READ_EXTERNAL_SEP_FILE, &
                  outputgrid,OUTPUT_ENERGY,TURN_ATTENUATION_ON, &
                  TURN_VISCATTENUATION_ON,Q0,freq0,p_sv, &
                  NSTEP,deltat,NTSTEP_BETWEEN_OUTPUT_SEISMO,NSOURCE)

! starts reading in parameters from input Database file

  implicit none
  include "constants.h"

  integer :: myrank,ipass
  character(len=60) simulation_title
  integer :: SIMULATION_TYPE,npgeo
  integer :: colors,numbers,subsamp,seismotype,imagetype
  logical :: SAVE_FORWARD,gnuplot,interpol,output_postscript_snapshot, &
    output_color_image
  logical :: meshvect,modelvect,boundvect,initialfield,add_Bielak_conditions, &
    assign_external_model,READ_EXTERNAL_SEP_FILE, &
    outputgrid,OUTPUT_ENERGY,p_sv
  logical :: TURN_ATTENUATION_ON,TURN_VISCATTENUATION_ON

  double precision :: cutsnaps,sizemax_arrows,anglerec
  double precision :: Q0,freq0
  double precision :: deltat

  integer :: NSTEP,NSOURCE
  integer :: NTSTEP_BETWEEN_OUTPUT_INFO,NTSTEP_BETWEEN_OUTPUT_SEISMO

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
  read(IIN,*) SIMULATION_TYPE, SAVE_FORWARD

  read(IIN,"(a80)") datlin
  read(IIN,*) npgeo

  read(IIN,"(a80)") datlin
  read(IIN,*) gnuplot,interpol

  read(IIN,"(a80)") datlin
  read(IIN,*) NTSTEP_BETWEEN_OUTPUT_INFO

  read(IIN,"(a80)") datlin
  read(IIN,*) output_postscript_snapshot,output_color_image,colors,numbers

  read(IIN,"(a80)") datlin
  read(IIN,*) meshvect,modelvect,boundvect,cutsnaps,subsamp,sizemax_arrows
  cutsnaps = cutsnaps / 100.d0

  read(IIN,"(a80)") datlin
  read(IIN,*) anglerec

  read(IIN,"(a80)") datlin
  read(IIN,*) initialfield,add_Bielak_conditions
  if(add_Bielak_conditions .and. .not. initialfield) &
    stop 'need to have an initial field to add Bielak plane wave conditions'

  read(IIN,"(a80)") datlin
  read(IIN,*) seismotype,imagetype
  if(seismotype < 1 .or. seismotype > 6) call exit_MPI('Wrong type for seismogram output')
  if(imagetype < 1 .or. imagetype > 4) call exit_MPI('Wrong type for snapshots')

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
  read(IIN,*) outputgrid,OUTPUT_ENERGY,TURN_ATTENUATION_ON

  read(IIN,"(a80)") datlin
  read(IIN,*) TURN_VISCATTENUATION_ON,Q0,freq0

  read(IIN,"(a80)") datlin
  read(IIN,*) p_sv

  !---- check parameters read
  if (myrank == 0 .and. ipass == 1) then
    write(IOUT,200) npgeo,NDIM
    write(IOUT,600) NTSTEP_BETWEEN_OUTPUT_INFO,colors,numbers
    write(IOUT,700) seismotype,anglerec
    write(IOUT,750) initialfield,add_Bielak_conditions,assign_external_model,&
                    READ_EXTERNAL_SEP_FILE,TURN_ATTENUATION_ON, &
                    outputgrid,OUTPUT_ENERGY
    write(IOUT,800) imagetype,100.d0*cutsnaps,subsamp
  endif

  !---- read time step
  read(IIN,"(a80)") datlin
  read(IIN,*) NSTEP,deltat
  if (myrank == 0 .and. ipass == 1) write(IOUT,703) NSTEP,deltat,NSTEP*deltat

  if( SIMULATION_TYPE == 1 .and. SAVE_FORWARD .and. &
    (TURN_ATTENUATION_ON .or. TURN_VISCATTENUATION_ON) ) then
    print*, '*************** WARNING ***************'
    print*, 'Anisotropy & Attenuation & Viscous damping are not presently implemented for adjoint calculations'
    stop
  endif

  NTSTEP_BETWEEN_OUTPUT_SEISMO = min(NSTEP,NTSTEP_BETWEEN_OUTPUT_INFO)

  !----  read source information
  read(IIN,"(a80)") datlin
  read(IIN,*) NSOURCE

  ! output formats
230 format('./OUTPUT_FILES/Database',i5.5)
  
200 format(//1x,'C o n t r o l',/1x,13('='),//5x,&
  'Number of spectral element control nodes. . .(npgeo) =',i8/5x, &
  'Number of space dimensions. . . . . . . . . . (NDIM) =',i8)

600 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Display frequency . . . (NTSTEP_BETWEEN_OUTPUT_INFO) = ',i6/ 5x, &
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
  'Turn attenuation on or off. . .(TURN_ATTENUATION_ON) = ',l6/5x, &
  'Save grid in external file or not. . . .(outputgrid) = ',l6/5x, &
  'Save a file with total energy or not.(OUTPUT_ENERGY) = ',l6)

800 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Vector display type . . . . . . . . . . .(imagetype) = ',i6/5x, &
  'Percentage of cut for vector plots . . . .(cutsnaps) = ',f6.2/5x, &
  'Subsampling for velocity model display. . .(subsamp) = ',i6)

703 format(//' I t e r a t i o n s '/1x,19('='),//5x, &
      'Number of time iterations . . . . .(NSTEP) =',i8,/5x, &
      'Time step increment. . . . . . . .(deltat) =',1pe15.6,/5x, &
      'Total simulation duration . . . . . (ttot) =',1pe15.6)

  end subroutine read_databases_init

!
!-------------------------------------------------------------------------------------------------
!

  subroutine read_databases_sources(NSOURCE,source_type,time_function_type, &
                      x_source,z_source,Mxx,Mzz,Mxz,f0,tshift_src,factor,angleforce)

! reads source parameters

  implicit none
  include "constants.h"

  integer :: NSOURCE
  integer, dimension(NSOURCE) :: source_type,time_function_type
  double precision, dimension(NSOURCE) :: x_source,z_source, &
    Mxx,Mzz,Mxz,f0,tshift_src,factor,angleforce

  ! local parameters
  integer :: i_source
  character(len=80) :: datlin

  ! reads in source info from Database file
  do i_source=1,NSOURCE
     read(IIN,"(a80)") datlin
     read(IIN,*) source_type(i_source),time_function_type(i_source), &
                 x_source(i_source),z_source(i_source),f0(i_source),tshift_src(i_source), &
                 factor(i_source),angleforce(i_source),Mxx(i_source),Mzz(i_source),Mxz(i_source)
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

  ! reads the spectral macrobloc nodal coordinates
  ipoin = 0
  read(IIN,"(a80)") datlin

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
107 format(/5x,'--> Isoparametric Spectral Elements <--',//)

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
                                perm,antecedent_list)

! reads spectral macrobloc data

  implicit none
  include "constants.h"

  integer :: ipass,ngnod,nspec
  integer, dimension(nspec) :: kmato
  integer, dimension(ngnod,nspec) :: knods
  
  integer, dimension(nspec) :: perm,antecedent_list

  ! local parameters
  integer :: n,k,ispec,kmato_read
  integer, dimension(:), allocatable :: knods_read  
  character(len=80) :: datlin

  ! reads spectral macrobloc data
  n = 0
  read(IIN,"(a80)") datlin
  allocate(knods_read(ngnod))
  do ispec = 1,nspec
    read(IIN,*) n,kmato_read,(knods_read(k), k=1,ngnod)
    if(ipass == 1) then
      kmato(n) = kmato_read
      knods(:,n)= knods_read(:)
    else if(ipass == 2) then
      kmato(perm(antecedent_list(n))) = kmato_read
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

  ! reads in interfaces
  do num_interface = 1, ninterface
    read(IIN,*) my_neighbours(num_interface), my_nelmnts_neighbours(num_interface)
    do ie = 1, my_nelmnts_neighbours(num_interface)
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
                            ibegin_bottom,iend_bottom,jbegin_right,jend_right, &
                            ibegin_top,iend_top,jbegin_left,jend_left, &
                            numabs,codeabs,perm,antecedent_list)

! reads in absorbing edges

  implicit none
  include "constants.h"

  integer :: myrank,ipass,nspec
  integer :: nelemabs
  integer, dimension(nelemabs) :: numabs,ibegin_bottom,iend_bottom, &
    ibegin_top,iend_top,jbegin_left,jend_left,jbegin_right,jend_right
  logical, dimension(4,nelemabs) :: codeabs
  logical :: anyabs
  integer, dimension(nspec) :: perm,antecedent_list

  ! local parameters
  integer :: inum,numabsread
  logical :: codeabsread(4)
  character(len=80) :: datlin

  ! reads in absorbing edges
  read(IIN,"(a80)") datlin

  if( anyabs ) then
    do inum = 1,nelemabs
      read(IIN,*) numabsread,codeabsread(1),codeabsread(2),codeabsread(3),&
                  codeabsread(4), ibegin_bottom(inum), iend_bottom(inum), &
                  jbegin_right(inum), jend_right(inum), ibegin_top(inum), &
                  iend_top(inum), jbegin_left(inum), jend_left(inum)

      if(numabsread < 1 .or. numabsread > nspec) &
        call exit_MPI('Wrong absorbing element number')

      if(ipass == 1) then
        numabs(inum) = numabsread
      else if(ipass == 2) then
        numabs(inum) = perm(antecedent_list(numabsread))
      else
        call exit_MPI('error: maximum number of passes is 2')
      endif

      codeabs(IBOTTOM,inum) = codeabsread(1)
      codeabs(IRIGHT,inum) = codeabsread(2)
      codeabs(ITOP,inum) = codeabsread(3)
      codeabs(ILEFT,inum) = codeabsread(4)
    enddo

    if (myrank == 0 .and. ipass == 1) then
      write(IOUT,*)
      write(IOUT,*) 'Number of absorbing elements: ',nelemabs
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

  
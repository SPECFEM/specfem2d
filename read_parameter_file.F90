
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

module parameter_file

  ! note: we use this module definition only to be able to allocate
  !          arrays for receiverlines and materials in this subroutine rather than in the main
  !          routine in meshfem2D.F90
  
  ! note 2: the filename ending is .F90 to have pre-compilation with pragmas 
  !            (like #ifndef USE_MPI) working properly
  
  implicit none
  character(len=100) :: interfacesfile,title
  
  integer :: SIMULATION_TYPE
  logical :: SAVE_FORWARD,read_external_mesh

  character(len=256) :: mesh_file, nodes_coords_file, materials_file, &
                        free_surface_file, absorbing_surface_file
  character(len=256)  :: tangential_detection_curve_file

  ! variables used for partitioning
  integer :: nproc,partitioning_method

  double precision :: xmin,xmax
  integer :: nx,ngnod

  logical :: initialfield,add_Bielak_conditions,assign_external_model, &
            READ_EXTERNAL_SEP_FILE,TURN_ATTENUATION_ON,TURN_VISCATTENUATION_ON

  double precision :: Q0,freq0

  logical :: p_sv
  logical :: any_abs,absbottom,absright,abstop,absleft
  
  integer :: nt
  double precision :: deltat

  integer :: NSOURCE
  logical :: force_normal_to_surface

  ! variables used for attenuation
  integer  :: N_SLS
  double precision  :: f0_attenuation

  integer :: seismotype
  logical :: generate_STATIONS
  
  integer :: nreceiverlines
  double precision :: anglerec  
  logical :: rec_normal_to_surface

  integer, dimension(:), pointer :: nrec
  double precision, dimension(:), pointer :: xdeb,zdeb,xfin,zfin
  logical, dimension(:), pointer :: enreg_surf
  
  integer :: NTSTEP_BETWEEN_OUTPUT_INFO
  logical :: output_postscript_snapshot,output_color_image
  integer :: imagetype
  double precision :: cutsnaps
  logical :: meshvect,modelvect,boundvect,interpol
  integer :: pointsdisp,subsamp
  double precision :: sizemax_arrows
  logical :: gnuplot,outputgrid,OUTPUT_ENERGY
  logical :: plot_lowerleft_corner_only

  ! to store density and velocity model
  integer :: nb_materials
  integer, dimension(:),pointer :: icodemat
  double precision, dimension(:),pointer :: rho_s,cp,cs, &
    aniso3,aniso4,aniso5,aniso6,aniso7,aniso8,Qp,Qs
  double precision, dimension(:),pointer :: rho_f,phi,tortuosity,permxx,permxz,&
       permzz,kappa_s,kappa_f,kappa_fr,eta_f,mu_fr

contains

  subroutine read_parameter_file()
  
! reads in DATA/Par_file
  
  implicit none
  include "constants.h"
  
  ! local parameters
  integer :: ios,ireceiverlines

  ! read file names and path for output
  call read_value_string(IIN,IGNORE_JUNK,title)
  call read_value_string(IIN,IGNORE_JUNK,interfacesfile)

  write(*,*) 'Title of the simulation'
  write(*,*) title
  print *

  ! read type of simulation
  call read_value_integer(IIN,IGNORE_JUNK,SIMULATION_TYPE)
  call read_value_logical(IIN,IGNORE_JUNK,SAVE_FORWARD)

  ! read info about external mesh
  call read_value_logical(IIN,IGNORE_JUNK,read_external_mesh)
  call read_value_string(IIN,IGNORE_JUNK,mesh_file)
  call read_value_string(IIN,IGNORE_JUNK,nodes_coords_file)
  call read_value_string(IIN,IGNORE_JUNK,materials_file)
  call read_value_string(IIN,IGNORE_JUNK,free_surface_file)
  call read_value_string(IIN,IGNORE_JUNK,absorbing_surface_file)
  call read_value_string(IIN,IGNORE_JUNK,tangential_detection_curve_file)

  ! read info about partitioning
  call read_value_integer(IIN,IGNORE_JUNK,nproc)
  call read_value_integer(IIN,IGNORE_JUNK,partitioning_method)

  ! read grid parameters
  call read_value_double_precision(IIN,IGNORE_JUNK,xmin)
  call read_value_double_precision(IIN,IGNORE_JUNK,xmax)
  call read_value_integer(IIN,IGNORE_JUNK,nx)
  call read_value_integer(IIN,IGNORE_JUNK,ngnod)
  call read_value_logical(IIN,IGNORE_JUNK,initialfield)
  call read_value_logical(IIN,IGNORE_JUNK,add_Bielak_conditions)
  call read_value_logical(IIN,IGNORE_JUNK,assign_external_model)
  call read_value_logical(IIN,IGNORE_JUNK,READ_EXTERNAL_SEP_FILE)
  call read_value_logical(IIN,IGNORE_JUNK,TURN_ATTENUATION_ON)
  ! read viscous attenuation parameters (poroelastic media)
  call read_value_logical(IIN,IGNORE_JUNK,TURN_VISCATTENUATION_ON)
  call read_value_double_precision(IIN,IGNORE_JUNK,Q0)
  call read_value_double_precision(IIN,IGNORE_JUNK,freq0)
  ! determine if body or surface (membrane) waves calculation
  call read_value_logical(IIN,IGNORE_JUNK,p_sv)

  ! read absorbing boundaries parameters
  call read_value_logical(IIN,IGNORE_JUNK,any_abs)
  call read_value_logical(IIN,IGNORE_JUNK,absbottom)
  call read_value_logical(IIN,IGNORE_JUNK,absright)
  call read_value_logical(IIN,IGNORE_JUNK,abstop)
  call read_value_logical(IIN,IGNORE_JUNK,absleft)

  ! read time step parameters
  call read_value_integer(IIN,IGNORE_JUNK,nt)
  call read_value_double_precision(IIN,IGNORE_JUNK,deltat)
  
  ! read source infos
  call read_value_integer(IIN,IGNORE_JUNK,NSOURCE)
  call read_value_logical(IIN,IGNORE_JUNK,force_normal_to_surface)

  ! read constants for attenuation
  call read_value_integer(IIN,IGNORE_JUNK,N_SLS)
  call read_value_double_precision(IIN,IGNORE_JUNK,f0_attenuation)

  ! read receiver line parameters
  call read_value_integer(IIN,IGNORE_JUNK,seismotype)
  call read_value_logical(IIN,IGNORE_JUNK,generate_STATIONS)
  call read_value_integer(IIN,IGNORE_JUNK,nreceiverlines)
  call read_value_double_precision(IIN,IGNORE_JUNK,anglerec)
  call read_value_logical(IIN,IGNORE_JUNK,rec_normal_to_surface)

  if(nreceiverlines < 1) stop 'number of receiver lines must be greater than 1'

  ! allocate receiver line arrays
  allocate(nrec(nreceiverlines))
  allocate(xdeb(nreceiverlines))
  allocate(zdeb(nreceiverlines))
  allocate(xfin(nreceiverlines))
  allocate(zfin(nreceiverlines))
  allocate(enreg_surf(nreceiverlines),stat=ios)
  if( ios /= 0 ) stop 'error allocating receiver lines'

  ! loop on all the receiver lines
  do ireceiverlines = 1,nreceiverlines
     call read_value_integer(IIN,IGNORE_JUNK,nrec(ireceiverlines))
     call read_value_double_precision(IIN,IGNORE_JUNK,xdeb(ireceiverlines))
     call read_value_double_precision(IIN,IGNORE_JUNK,zdeb(ireceiverlines))
     call read_value_double_precision(IIN,IGNORE_JUNK,xfin(ireceiverlines))
     call read_value_double_precision(IIN,IGNORE_JUNK,zfin(ireceiverlines))
     call read_value_logical(IIN,IGNORE_JUNK,enreg_surf(ireceiverlines))
     if (read_external_mesh .and. enreg_surf(ireceiverlines)) then
        stop 'Cannot use enreg_surf with external meshes!'
     endif
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
  call read_value_logical(IIN,IGNORE_JUNK,OUTPUT_ENERGY)


  ! read the different material materials
  call read_value_integer(IIN,IGNORE_JUNK,nb_materials)
  if(nb_materials <= 0) stop 'Negative number of materials not allowed!'

  allocate(icodemat(nb_materials))
  allocate(cp(nb_materials))
  allocate(cs(nb_materials))
  allocate(aniso3(nb_materials))
  allocate(aniso4(nb_materials))
  allocate(aniso5(nb_materials))
  allocate(aniso6(nb_materials))
  allocate(aniso7(nb_materials))
  allocate(aniso8(nb_materials))
  allocate(Qp(nb_materials))
  allocate(Qs(nb_materials))
  allocate(rho_s(nb_materials))
  allocate(rho_f(nb_materials))
  allocate(phi(nb_materials))
  allocate(tortuosity(nb_materials))
  allocate(permxx(nb_materials))
  allocate(permxz(nb_materials))
  allocate(permzz(nb_materials))
  allocate(kappa_s(nb_materials))
  allocate(kappa_f(nb_materials))
  allocate(kappa_fr(nb_materials))
  allocate(eta_f(nb_materials))
  allocate(mu_fr(nb_materials))

  call read_materials(nb_materials,icodemat,cp,cs, &
                      aniso3,aniso4,aniso5,aniso6,aniso7,aniso8, &
                      Qp,Qs,rho_s,rho_f,phi,tortuosity, &
                      permxx,permxz,permzz,kappa_s,kappa_f,kappa_fr, &
                      eta_f,mu_fr)
  

  ! checks input parameters
  call check_parameters()
    
  end subroutine read_parameter_file
  
!
!-------------------------------------------------------------------------------------------------
!
  
  subroutine check_parameters()
  
  implicit none
  
  ! checks partitioning
  if ( nproc <= 0 ) then
     print *, 'Number of processes (nproc) must be greater than or equal to one.'
     stop
  endif

#ifndef USE_MPI
  if ( nproc > 1 ) then
     print *, 'Number of processes (nproc) must be equal to one when not using MPI.'
     print *, 'Please recompile with -DUSE_MPI in order to enable use of MPI.'
     stop
  endif
#endif

  if(partitioning_method /= 1 .and. partitioning_method /= 3) then
     print *, 'Invalid partitioning method number.'
     print *, 'Partitioning method ',partitioning_method,' was requested, but is not available.'
     print *, 'Support for the METIS graph partitioner has been discontinued, please use SCOTCH (option 3) instead.'
     stop
  endif

  ! checks absorbing boundaries
  if ( .not. any_abs ) then
     absbottom = .false.
     absright = .false.
     abstop = .false.
     absleft = .false.
  endif

  ! can use only one point to display lower-left corner only for interpolated snapshot
  if(pointsdisp < 3) then
     pointsdisp = 3
     plot_lowerleft_corner_only = .true.
  else
     plot_lowerleft_corner_only = .false.
  endif
  
  end subroutine check_parameters
  
end module parameter_file
  

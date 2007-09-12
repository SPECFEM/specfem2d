
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

!====================================================================================
!
! An explicit 2D spectral element solver for the anelastic anisotropic wave equation
!
!====================================================================================

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

!
! version 5.2, Dimitri Komatitsch, April 2007:
!               - general fluid/solid implementation with any number, shape and orientation of
!                 matching edges
!               - absorbing edges with any normal vector
!               - general numbering of absorbing and acoustic free surface edges
!               - cleaned implementation of attenuation as in Carcione (1993)
!               - merged loops in the solver for efficiency
!               - simplified input of external model
!               - added CPU time information
!               - translated many comments from French to English
!
! version 5.1, Dimitri Komatitsch, January 2005:
!               - more general mesher with any number of curved layers
!               - Dirac and Gaussian time sources and corresponding convolution routine
!               - option for acoustic medium instead of elastic
!               - receivers at any location, not only grid points
!               - moment-tensor source at any location, not only a grid point
!               - color snapshots
!               - more flexible DATA/Par_file with any number of comment lines
!               - Xsu scripts for seismograms
!               - subtract t0 from seismograms
!               - seismograms and snapshots in pressure in addition to vector field
!
! version 5.0, Dimitri Komatitsch, May 2004:
!               - got rid of useless routines, suppressed commons etc.
!               - weak formulation based explicitly on stress tensor
!               - implementation of full anisotropy
!               - implementation of attenuation based on memory variables
!
! based on SPECFEM2D version 4.2, June 1998
! (c) by Dimitri Komatitsch, Harvard University, USA
! and Jean-Pierre Vilotte, Institut de Physique du Globe de Paris, France
!
! itself based on SPECFEM2D version 1.0, 1995
! (c) by Dimitri Komatitsch and Jean-Pierre Vilotte,
! Institut de Physique du Globe de Paris, France
!

! in case of an acoustic medium, a displacement potential Chi is used as in Chaljub and Valette,
! Geophysical Journal International, vol. 158, p. 131-141 (2004) and *NOT* a velocity potential
! as in Komatitsch and Tromp, Geophysical Journal International, vol. 150, p. 303-318 (2002).
! This permits acoustic-elastic coupling based on a non-iterative time scheme.
! Displacement is then: u = grad(Chi)
! Velocity is then: v = grad(Chi_dot)       (Chi_dot being the time derivative of Chi)
! and pressure is: p = - rho * Chi_dot_dot  (Chi_dot_dot being the time second derivative of Chi).
! The source in an acoustic element is a pressure source.

  program specfem2D

  implicit none

  include "constants.h"
#ifdef USE_MPI
  include "mpif.h"
#endif

  character(len=80) datlin

  integer :: source_type,time_function_type
  double precision :: x_source,z_source,xi_source,gamma_source,Mxx,Mzz,Mxz,f0,t0,factor,angleforce,hdur,hdur_gauss
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLZ) :: sourcearray

  double precision, dimension(:,:), allocatable :: coorg
  double precision, dimension(:), allocatable :: coorgread

! receiver information
  integer, dimension(:), allocatable :: ispec_selected_rec
  double precision, dimension(:), allocatable :: xi_receiver,gamma_receiver,st_xval,st_zval

! for seismograms
  double precision, dimension(:,:), allocatable :: sisux,sisuz

! vector field in an element
  real(kind=CUSTOM_REAL), dimension(NDIM,NGLLX,NGLLX) :: vector_field_element

! pressure in an element
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: pressure_element

  integer :: i,j,k,l,it,irec,ipoin,ip,id,nbpoin,inump,n,ispec,npoin,npgeo,iglob
  logical :: anyabs
  double precision :: dxd,dzd,valux,valuz,hlagrange,rhol,cosrot,sinrot,xi,gamma,x,z

! coefficients of the explicit Newmark time scheme
  integer NSTEP
  double precision deltatover2,deltatsquareover2,time,deltat

! Gauss-Lobatto-Legendre points and weights
  double precision, dimension(NGLLX) :: xigll
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  double precision, dimension(NGLLZ) :: zigll
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wzgll

! derivatives of Lagrange polynomials
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLX) :: hprime_xx,hprimewgll_xx
  real(kind=CUSTOM_REAL), dimension(NGLLZ,NGLLZ) :: hprime_zz,hprimewgll_zz

! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl,jacobianl

! material properties of the elastic medium
  double precision :: mul_relaxed,lambdal_relaxed,cpsquare

  real(kind=CUSTOM_REAL), dimension(:,:), allocatable :: accel_elastic,veloc_elastic,displ_elastic
  double precision, dimension(:,:), allocatable :: coord, flagrange,xinterp,zinterp,Uxinterp,Uzinterp,elastcoef,vector_field_display

! for acoustic medium
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: potential_dot_dot_acoustic,potential_dot_acoustic,potential_acoustic
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: rmass_inverse_elastic,rmass_inverse_acoustic
  double precision, dimension(:), allocatable :: density,displread,velocread,accelread

  double precision, dimension(:), allocatable :: vp_display

  double precision, dimension(:,:,:), allocatable :: vpext,vsext,rhoext
  double precision :: previous_vsext

  double precision, dimension(:,:,:), allocatable :: shape2D,shape2D_display
  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable  :: xix,xiz,gammax,gammaz,jacobian

  double precision, dimension(:,:,:,:), allocatable :: dershape2D,dershape2D_display

  integer, dimension(:,:,:), allocatable :: ibool
  integer, dimension(:,:), allocatable  :: knods
  integer, dimension(:), allocatable :: kmato,numabs, &
     ibegin_bottom,iend_bottom,ibegin_top,iend_top,jbegin_left,jend_left,jbegin_right,jend_right

  integer ispec_selected_source,iglob_source,ix_source,iz_source,is_proc_source,nb_proc_source
  double precision a,displnorm_all
  real(kind=CUSTOM_REAL), dimension(:), allocatable :: source_time_function
  double precision, external :: erf

  double precision :: vpmin,vpmax

  integer :: colors,numbers,subsamp,imagetype,NTSTEP_BETWEEN_OUTPUT_INFO,nrec,seismotype
  integer :: numat,ngnod,nspec,pointsdisp,nelemabs,nelem_acoustic_surface,ispecabs

  logical interpol,meshvect,modelvect,boundvect,assign_external_model,initialfield, &
    outputgrid,gnuplot,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON,output_postscript_snapshot,output_color_image, &
    plot_lowerleft_corner_only

  double precision :: cutsnaps,sizemax_arrows,anglerec,xirec,gammarec

! for absorbing and acoustic free surface conditions
  integer :: ispec_acoustic_surface,inum,numabsread
  logical :: codeabsread(4)
  real(kind=CUSTOM_REAL) :: nx,nz,weight,xxi,zgamma

  logical, dimension(:,:), allocatable  :: codeabs

! for attenuation
  double precision  :: Qp_attenuation
  double precision  :: Qs_attenuation
  double precision  :: f0_attenuation
  integer nspec_allocate
  double precision :: deltatsquare,deltatcube,deltatfourth,twelvedeltat,fourdeltatsquare

  real(kind=CUSTOM_REAL), dimension(:,:,:,:), allocatable :: e1,e11,e13
  double precision, dimension(N_SLS) :: inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
  double precision :: Mu_nu1,Mu_nu2

  real(kind=CUSTOM_REAL), dimension(:,:,:), allocatable :: &
    dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n,dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1

! for fluid/solid coupling and edge detection
  logical, dimension(:), allocatable :: elastic
  integer, dimension(NEDGES) :: i_begin,j_begin,i_end,j_end
  integer, dimension(NGLLX,NEDGES) :: ivalue,jvalue,ivalue_inverse,jvalue_inverse
  integer, dimension(:), allocatable :: fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge, &
                                        fluid_solid_elastic_ispec,fluid_solid_elastic_iedge
  integer :: num_fluid_solid_edges,ispec_acoustic,ispec_elastic, &
             iedge_acoustic,iedge_elastic,ipoin1D,iglob2
  logical :: any_acoustic,any_elastic,any_elastic_glob,coupled_acoustic_elastic
  real(kind=CUSTOM_REAL) :: displ_x,displ_z,displ_n,zxi,xgamma,jacobian1D,pressure

! for color images
  integer :: NX_IMAGE_color,NZ_IMAGE_color
  integer  :: npgeo_glob
  double precision :: xmin_color_image,xmax_color_image, &
    zmin_color_image,zmax_color_image,size_pixel_horizontal,size_pixel_vertical
  integer, dimension(:,:), allocatable :: iglob_image_color,copy_iglob_image_color
  double precision, dimension(:,:), allocatable :: image_color_data
  double precision, dimension(:,:), allocatable :: image_color_vp_display

  double precision  :: xmin_color_image_loc, xmax_color_image_loc, zmin_color_image_loc, &
       zmax_color_image_loc
  integer  :: min_i, min_j, max_i, max_j
  integer  :: nb_pixel_loc
  integer, dimension(:), allocatable  :: nb_pixel_per_proc
  double precision  :: i_coord, j_coord
  double precision, dimension(2,4)  :: elmnt_coords
  integer, dimension(:), allocatable  :: num_pixel_loc
  integer, dimension(:,:), allocatable  :: num_pixel_recv
  double precision, dimension(:), allocatable  :: data_pixel_recv
  double precision, dimension(:), allocatable  :: data_pixel_send
  logical  :: pixel_is_in
  double precision  :: dist_pixel, dist_min_pixel
#ifdef USE_MPI
  integer, dimension(MPI_STATUS_SIZE)  :: request_mpi_status
#endif

! timing information for the stations
  character(len=MAX_LENGTH_STATION_NAME), allocatable, dimension(:) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), allocatable, dimension(:) :: network_name

! title of the plot
  character(len=60) simulation_title

! Lagrange interpolators at receivers
  double precision, dimension(:), allocatable :: hxir,hgammar,hpxir,hpgammar
  double precision, dimension(:,:), allocatable :: hxir_store,hgammar_store

! for Lagrange interpolants
  double precision, external :: hgll

! timer to count elapsed time
  character(len=8) datein
  character(len=10) timein
  character(len=5)  :: zone
  integer, dimension(8) :: time_values
  integer ihours,iminutes,iseconds,int_tCPU
  double precision :: time_start,time_end,tCPU

! for MPI and partitionning
  integer  :: ier
  integer  :: nproc
  integer  :: myrank
  integer  :: iproc
  character(len=256)  :: prname

  integer  :: ninterface
  integer  :: max_interface_size
  integer, dimension(:), allocatable  :: my_neighbours
  integer, dimension(:), allocatable  :: my_nelmnts_neighbours
  integer, dimension(:,:,:), allocatable  :: my_interfaces
  integer, dimension(:,:), allocatable  :: ibool_interfaces_acoustic,ibool_interfaces_elastic
  integer, dimension(:), allocatable  :: nibool_interfaces_acoustic,nibool_interfaces_elastic

  integer  :: ninterface_acoustic, ninterface_elastic
  integer, dimension(:), allocatable  :: inum_interfaces_acoustic, inum_interfaces_elastic

#ifdef USE_MPI
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_ac
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_ac
  integer, dimension(:), allocatable  :: tab_requests_send_recv_acoustic
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_send_faces_vector_el
  real(kind=CUSTOM_REAL), dimension(:,:), allocatable  :: buffer_recv_faces_vector_el
  integer, dimension(:), allocatable  :: tab_requests_send_recv_elastic
  integer  :: max_ibool_interfaces_size_ac, max_ibool_interfaces_size_el
#endif

! for overlapping MPI communications with computation
  integer  :: nspec_outer, nspec_inner, num_ispec_outer, num_ispec_inner
  integer, dimension(:), allocatable  :: ispec_outer_to_glob, ispec_inner_to_glob
  logical, dimension(:), allocatable  :: mask_ispec_inner_outer

  integer, dimension(:,:), allocatable  :: acoustic_surface
  integer, dimension(:,:), allocatable  :: acoustic_edges

  integer  :: ixmin, ixmax, izmin, izmax

  integer  :: ie, num_interface

  integer  :: nrecloc, irecloc
  integer, dimension(:), allocatable :: recloc, which_proc_receiver

! mask to sort ibool
  integer, dimension(:), allocatable :: mask_ibool
  integer, dimension(:,:,:), allocatable :: copy_ibool_ori
  integer :: inumber

!***********************************************************************
!
!             i n i t i a l i z a t i o n    p h a s e
!
!***********************************************************************

#ifdef USE_MPI
  call MPI_INIT(ier)
  call MPI_COMM_SIZE(MPI_COMM_WORLD,nproc,ier)
  call MPI_COMM_RANK(MPI_COMM_WORLD,myrank,ier)

#else
  nproc = 1
  myrank = 0
  ier = 0
  ninterface_acoustic = 0
  ninterface_elastic = 0
  iproc = 0

#endif

  write(prname,230)myrank
230   format('./OUTPUT_FILES/Database',i5.5)

  open(unit=IIN,file=prname,status='old',action='read')


! determine if we write to file instead of standard output
  if(IOUT /= ISTANDARD_OUTPUT) open(IOUT,file='simulation_results.txt',status='unknown')

!
!---  read job title and skip remaining titles of the input file
!
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a80)") datlin
  read(IIN,"(a50)") simulation_title

!
!---- print the date, time and start-up banner
!
  call datim(simulation_title)

  write(IOUT,*)
  write(IOUT,*)
  write(IOUT,*) '*********************'
  write(IOUT,*) '****             ****'
  write(IOUT,*) '****  SPECFEM2D  ****'
  write(IOUT,*) '****             ****'
  write(IOUT,*) '*********************'

!
!---- read parameters from input file
!

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
  read(IIN,*) initialfield

  read(IIN,"(a80)") datlin
  read(IIN,*) seismotype,imagetype
  if(seismotype < 1 .or. seismotype > 4) call exit_MPI('Wrong type for seismogram output')
  if(imagetype < 1 .or. imagetype > 4) call exit_MPI('Wrong type for snapshots')

  read(IIN,"(a80)") datlin
  read(IIN,*) assign_external_model,outputgrid,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON

!---- check parameters read
  write(IOUT,200) npgeo,NDIM
  write(IOUT,600) NTSTEP_BETWEEN_OUTPUT_INFO,colors,numbers
  write(IOUT,700) seismotype,anglerec
  write(IOUT,750) initialfield,assign_external_model,TURN_ANISOTROPY_ON,TURN_ATTENUATION_ON,outputgrid
  write(IOUT,800) imagetype,100.d0*cutsnaps,subsamp

!---- read time step
  read(IIN,"(a80)") datlin
  read(IIN,*) NSTEP,deltat
  write(IOUT,703) NSTEP,deltat,NSTEP*deltat

!
!----  read source information
!
  read(IIN,"(a80)") datlin
  read(IIN,*) source_type,time_function_type,x_source,z_source,f0,t0,factor,angleforce,Mxx,Mzz,Mxz

!
!----  read attenuation information
!
  read(IIN,"(a80)") datlin
  read(IIN,*) Qp_attenuation, Qs_attenuation, f0_attenuation

!
!-----  check the input
!
 if(.not. initialfield) then
   if (source_type == 1) then
     write(IOUT,212) x_source,z_source,f0,t0,factor,angleforce
   else if(source_type == 2) then
     write(IOUT,222) x_source,z_source,f0,t0,factor,Mxx,Mzz,Mxz
   else
     call exit_MPI('Unknown source type number !')
   endif
 endif

! for the source time function
  a = pi*pi*f0*f0

!-----  convert angle from degrees to radians
  angleforce = angleforce * pi / 180.d0

!
!---- read the spectral macrobloc nodal coordinates
!
  allocate(coorg(NDIM,npgeo))

  ipoin = 0
  read(IIN,"(a80)") datlin
  allocate(coorgread(NDIM))
  do ip = 1,npgeo
   read(IIN,*) ipoin,(coorgread(id),id =1,NDIM)
   if(ipoin<1 .or. ipoin>npgeo) call exit_MPI('Wrong control point number')
   coorg(:,ipoin) = coorgread
  enddo
  deallocate(coorgread)

!
!---- read the basic properties of the spectral elements
!
  read(IIN,"(a80)") datlin
  read(IIN,*) numat,ngnod,nspec,pointsdisp,plot_lowerleft_corner_only
  read(IIN,"(a80)") datlin
  read(IIN,*) nelemabs,nelem_acoustic_surface,num_fluid_solid_edges

!
!---- allocate arrays
!
  allocate(shape2D(ngnod,NGLLX,NGLLZ))
  allocate(dershape2D(NDIM,ngnod,NGLLX,NGLLZ))
  allocate(shape2D_display(ngnod,pointsdisp,pointsdisp))
  allocate(dershape2D_display(NDIM,ngnod,pointsdisp,pointsdisp))
  allocate(xix(NGLLX,NGLLZ,nspec))
  allocate(xiz(NGLLX,NGLLZ,nspec))
  allocate(gammax(NGLLX,NGLLZ,nspec))
  allocate(gammaz(NGLLX,NGLLZ,nspec))
  allocate(jacobian(NGLLX,NGLLZ,nspec))
  allocate(flagrange(NGLLX,pointsdisp))
  allocate(xinterp(pointsdisp,pointsdisp))
  allocate(zinterp(pointsdisp,pointsdisp))
  allocate(Uxinterp(pointsdisp,pointsdisp))
  allocate(Uzinterp(pointsdisp,pointsdisp))
  allocate(density(numat))
  allocate(elastcoef(4,numat))
  allocate(kmato(nspec))
  allocate(knods(ngnod,nspec))
  allocate(ibool(NGLLX,NGLLZ,nspec))
  allocate(elastic(nspec))

! --- allocate arrays for absorbing boundary conditions
  if(nelemabs <= 0) then
    nelemabs = 1
    anyabs = .false.
  else
    anyabs = .true.
  endif
  allocate(numabs(nelemabs))
  allocate(codeabs(4,nelemabs))

  allocate(ibegin_bottom(nelemabs))
  allocate(iend_bottom(nelemabs))
  allocate(ibegin_top(nelemabs))
  allocate(iend_top(nelemabs))

  allocate(jbegin_left(nelemabs))
  allocate(jend_left(nelemabs))
  allocate(jbegin_right(nelemabs))
  allocate(jend_right(nelemabs))

!
!---- print element group main parameters
!
  write(IOUT,107)
  write(IOUT,207) nspec,ngnod,NGLLX,NGLLZ,NGLLX*NGLLZ,pointsdisp,numat,nelemabs

! set up Gauss-Lobatto-Legendre derivation matrices
  call define_derivation_matrices(xigll,zigll,wxgll,wzgll,hprime_xx,hprime_zz,hprimewgll_xx,hprimewgll_zz)

!
!---- read the material properties
!
  call gmat01(density,elastcoef,numat)

!
!----  read spectral macrobloc data
!
  n = 0
  read(IIN,"(a80)") datlin
  do ispec = 1,nspec
    read(IIN,*) n,kmato(n),(knods(k,n), k=1,ngnod)
  enddo

!
!----  determine if each spectral element is elastic or acoustic
!
  any_acoustic = .false.
  any_elastic = .false.
  do ispec = 1,nspec
    mul_relaxed = elastcoef(2,kmato(ispec))
    if(mul_relaxed < TINYVAL) then
      elastic(ispec) = .false.
      any_acoustic = .true.
    else
      elastic(ispec) = .true.
      any_elastic = .true.
    endif
  enddo

  if(TURN_ATTENUATION_ON) then
    nspec_allocate = nspec
  else
    nspec_allocate = 1
  endif

! allocate memory variables for attenuation
  allocate(e1(NGLLX,NGLLZ,nspec_allocate,N_SLS))
  allocate(e11(NGLLX,NGLLZ,nspec_allocate,N_SLS))
  allocate(e13(NGLLX,NGLLZ,nspec_allocate,N_SLS))
  e1(:,:,:,:) = 0._CUSTOM_REAL
  e11(:,:,:,:) = 0._CUSTOM_REAL
  e13(:,:,:,:) = 0._CUSTOM_REAL

  allocate(dux_dxl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(duz_dzl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(duz_dxl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(dux_dzl_n(NGLLX,NGLLZ,nspec_allocate))
  allocate(dux_dxl_np1(NGLLX,NGLLZ,nspec_allocate))
  allocate(duz_dzl_np1(NGLLX,NGLLZ,nspec_allocate))
  allocate(duz_dxl_np1(NGLLX,NGLLZ,nspec_allocate))
  allocate(dux_dzl_np1(NGLLX,NGLLZ,nspec_allocate))

! define the attenuation constants
  call attenuation_model(Qp_attenuation,Qs_attenuation,f0_attenuation, &
      inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2)

!
!----  read interfaces data
!
  print *, 'read the interfaces', myrank
  read(IIN,"(a80)") datlin
  read(IIN,*) ninterface, max_interface_size
  if ( ninterface == 0 ) then
     !allocate(my_neighbours(1))
     !allocate(my_nelmnts_neighbours(1))
     !allocate(my_interfaces(4,1,1))
     !allocate(ibool_interfaces(NGLLX*1,1,1))
     !allocate(nibool_interfaces(1,1))

  else
     allocate(my_neighbours(ninterface))
     allocate(my_nelmnts_neighbours(ninterface))
     allocate(my_interfaces(4,max_interface_size,ninterface))
     allocate(ibool_interfaces_acoustic(NGLLX*max_interface_size,ninterface))
     allocate(ibool_interfaces_elastic(NGLLX*max_interface_size,ninterface))
     allocate(nibool_interfaces_acoustic(ninterface))
     allocate(nibool_interfaces_elastic(ninterface))
     allocate(inum_interfaces_acoustic(ninterface))
     allocate(inum_interfaces_elastic(ninterface))

     do num_interface = 1, ninterface
        read(IIN,*) my_neighbours(num_interface), my_nelmnts_neighbours(num_interface)
        do ie = 1, my_nelmnts_neighbours(num_interface)
           read(IIN,*) my_interfaces(1,ie,num_interface), my_interfaces(2,ie,num_interface), &
                my_interfaces(3,ie,num_interface), my_interfaces(4,ie,num_interface)

        end do
     end do
     print *, 'end read the interfaces', myrank

  end if


!
!----  read absorbing boundary data
!
  read(IIN,"(a80)") datlin
  if(anyabs) then
     do inum = 1,nelemabs
      read(IIN,*) numabsread,codeabsread(1),codeabsread(2),codeabsread(3),codeabsread(4), ibegin_bottom(inum), iend_bottom(inum), &
           jbegin_right(inum), jend_right(inum), ibegin_top(inum), iend_top(inum), jbegin_left(inum), jend_left(inum)
      if(numabsread < 1 .or. numabsread > nspec) call exit_MPI('Wrong absorbing element number')
      numabs(inum) = numabsread
      codeabs(IBOTTOM,inum) = codeabsread(1)
      codeabs(IRIGHT,inum) = codeabsread(2)
      codeabs(ITOP,inum) = codeabsread(3)
      codeabs(ILEFT,inum) = codeabsread(4)
    enddo
    write(IOUT,*)
    write(IOUT,*) 'Number of absorbing elements: ',nelemabs
  endif

!
!----  read acoustic free surface data
!
  read(IIN,"(a80)") datlin
  if(nelem_acoustic_surface > 0) then
     allocate(acoustic_edges(4,nelem_acoustic_surface))
      do inum = 1,nelem_acoustic_surface
        read(IIN,*) acoustic_edges(1,inum), acoustic_edges(2,inum), acoustic_edges(3,inum), &
             acoustic_edges(4,inum)
     end do
     allocate(acoustic_surface(5,nelem_acoustic_surface))
     call construct_acoustic_surface ( nspec, ngnod, knods, nelem_acoustic_surface, &
          acoustic_edges, acoustic_surface)
    write(IOUT,*)
    write(IOUT,*) 'Number of free surface elements: ',nelem_acoustic_surface
  endif

!
!---- read acoustic elastic coupled edges
!
  read(IIN,"(a80)") datlin
  if ( num_fluid_solid_edges > 0 ) then
     allocate(fluid_solid_acoustic_ispec(num_fluid_solid_edges))
     allocate(fluid_solid_acoustic_iedge(num_fluid_solid_edges))
     allocate(fluid_solid_elastic_ispec(num_fluid_solid_edges))
     allocate(fluid_solid_elastic_iedge(num_fluid_solid_edges))
     do inum = 1, num_fluid_solid_edges
        read(IIN,*) fluid_solid_acoustic_ispec(inum), fluid_solid_elastic_ispec(inum)
     end do
  else
     allocate(fluid_solid_acoustic_ispec(1))
     allocate(fluid_solid_acoustic_iedge(1))
     allocate(fluid_solid_elastic_ispec(1))
     allocate(fluid_solid_elastic_iedge(1))

  end if


!
!---- close input file
!
  close(IIN)

!
!---- compute shape functions and their derivatives for SEM grid
!
  do j = 1,NGLLZ
    do i = 1,NGLLX
      call define_shape_functions(shape2D(:,i,j),dershape2D(:,:,i,j),xigll(i),zigll(j),ngnod)
    enddo
  enddo

!
!---- generate the global numbering
!

! "slow and clean" or "quick and dirty" version
  if(FAST_NUMBERING) then
    call createnum_fast(knods,ibool,shape2D,coorg,npoin,npgeo,nspec,ngnod)
  else
    call createnum_slow(knods,ibool,npoin,nspec,ngnod)
  endif

! create a new indirect addressing array instead, to reduce cache misses
! in memory access in the solver
  allocate(copy_ibool_ori(NGLLX,NGLLZ,nspec))
  allocate(mask_ibool(npoin))
  mask_ibool(:) = -1
  copy_ibool_ori(:,:,:) = ibool(:,:,:)

  inumber = 0
  do ispec=1,nspec
    do j=1,NGLLZ
      do i=1,NGLLX
        if(mask_ibool(copy_ibool_ori(i,j,ispec)) == -1) then
! create a new point
          inumber = inumber + 1
          ibool(i,j,ispec) = inumber
          mask_ibool(copy_ibool_ori(i,j,ispec)) = inumber
        else
! use an existing point created previously
          ibool(i,j,ispec) = mask_ibool(copy_ibool_ori(i,j,ispec))
        endif
      enddo
    enddo
  enddo
  deallocate(copy_ibool_ori)
  deallocate(mask_ibool)

!---- compute shape functions and their derivatives for regular interpolated display grid
  do j = 1,pointsdisp
    do i = 1,pointsdisp
      xirec  = 2.d0*dble(i-1)/dble(pointsdisp-1) - 1.d0
      gammarec  = 2.d0*dble(j-1)/dble(pointsdisp-1) - 1.d0
      call define_shape_functions(shape2D_display(:,i,j),dershape2D_display(:,:,i,j),xirec,gammarec,ngnod)
    enddo
  enddo

!---- compute Lagrange interpolants on a regular interpolated grid in (xi,gamma)
!---- for display (assumes NGLLX = NGLLZ)
  do j=1,NGLLX
    do i=1,pointsdisp
      xirec  = 2.d0*dble(i-1)/dble(pointsdisp-1) - 1.d0
      flagrange(j,i) = hgll(j-1,xirec,xigll,NGLLX)
    enddo
  enddo

! read total number of receivers
  open(unit=IIN,file='DATA/STATIONS',status='old')
  read(IIN,*) nrec
  close(IIN)

  write(IOUT,*)
  write(IOUT,*) 'Total number of receivers = ',nrec
  write(IOUT,*)

  if(nrec < 1) call exit_MPI('need at least one receiver')

! receiver information
  allocate(ispec_selected_rec(nrec))
  allocate(st_xval(nrec))
  allocate(st_zval(nrec))
  allocate(xi_receiver(nrec))
  allocate(gamma_receiver(nrec))
  allocate(station_name(nrec))
  allocate(network_name(nrec))
  allocate(recloc(nrec))
  allocate(which_proc_receiver(nrec))

! allocate 1-D Lagrange interpolators and derivatives
  allocate(hxir(NGLLX))
  allocate(hpxir(NGLLX))
  allocate(hgammar(NGLLZ))
  allocate(hpgammar(NGLLZ))

! allocate Lagrange interpolators for receivers
  allocate(hxir_store(nrec,NGLLX))
  allocate(hgammar_store(nrec,NGLLZ))

! allocate other global arrays
  allocate(coord(NDIM,npoin))

! to display acoustic elements
  allocate(vector_field_display(NDIM,npoin))

  if(assign_external_model) then
    allocate(vpext(NGLLX,NGLLZ,nspec))
    allocate(vsext(NGLLX,NGLLZ,nspec))
    allocate(rhoext(NGLLX,NGLLZ,nspec))
  else
    allocate(vpext(1,1,1))
    allocate(vsext(1,1,1))
    allocate(rhoext(1,1,1))
  endif

!
!----  set the coordinates of the points of the global grid
!
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX

        xi = xigll(i)
        gamma = zigll(j)

        call recompute_jacobian(xi,gamma,x,z,xixl,xizl,gammaxl,gammazl,jacobianl,coorg,knods,ispec,ngnod,nspec,npgeo)

        coord(1,ibool(i,j,ispec)) = x
        coord(2,ibool(i,j,ispec)) = z

        xix(i,j,ispec) = xixl
        xiz(i,j,ispec) = xizl
        gammax(i,j,ispec) = gammaxl
        gammaz(i,j,ispec) = gammazl
        jacobian(i,j,ispec) = jacobianl

      enddo
    enddo
  enddo

!
!--- save the grid of points in a file
!
  if(outputgrid) then
    write(IOUT,*)
    write(IOUT,*) 'Saving the grid in a text file...'
    write(IOUT,*)
    open(unit=55,file='OUTPUT_FILES/grid_points_and_model.txt',status='unknown')
    write(55,*) npoin
    do n = 1,npoin
      write(55,*) (coord(i,n), i=1,NDIM)
    enddo
    close(55)
  endif

!
!-----   plot the GLL mesh in a Gnuplot file
!
  if(gnuplot) call plotgll(knods,ibool,coorg,coord,npoin,npgeo,ngnod,nspec)

!
!----  assign external velocity and density model if needed
!
  if(assign_external_model) then
    write(IOUT,*)
    write(IOUT,*) 'Assigning external velocity and density model...'
    write(IOUT,*)
    if(TURN_ANISOTROPY_ON .or. TURN_ATTENUATION_ON) &
         call exit_MPI('cannot have anisotropy nor attenuation if external model in current version')
    any_acoustic = .false.
    any_elastic = .false.
    do ispec = 1,nspec
      previous_vsext = -1.d0
      do j = 1,NGLLZ
        do i = 1,NGLLX
          iglob = ibool(i,j,ispec)
          call define_external_model(coord(1,iglob),coord(2,iglob),kmato(ispec), &
                                         rhoext(i,j,ispec),vpext(i,j,ispec),vsext(i,j,ispec),myrank)
! stop if the same element is assigned both acoustic and elastic points in external model
          if(.not. (i == 1 .and. j == 1) .and. &
            ((vsext(i,j,ispec) >= TINYVAL .and. previous_vsext < TINYVAL) .or. &
             (vsext(i,j,ispec) < TINYVAL .and. previous_vsext >= TINYVAL)))  &
                call exit_MPI('external velocity model cannot be both fluid and solid inside the same spectral element')
          if(vsext(i,j,ispec) < TINYVAL) then
            elastic(ispec) = .false.
            any_acoustic = .true.
          else
            elastic(ispec) = .true.
            any_elastic = .true.
          endif
          previous_vsext = vsext(i,j,ispec)
        enddo
      enddo
    enddo
  endif

!
!----  perform basic checks on parameters read
!
any_elastic_glob = any_elastic
#ifdef USE_MPI
  call MPI_ALLREDUCE(any_elastic, any_elastic_glob, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ier)
#endif

! for acoustic
  if(TURN_ANISOTROPY_ON .and. .not. any_elastic_glob) &
    call exit_MPI('cannot have anisotropy if acoustic simulation only')

  if(TURN_ATTENUATION_ON .and. .not. any_elastic_glob) &
    call exit_MPI('currently cannot have attenuation if acoustic simulation only')

! for attenuation
  if(TURN_ANISOTROPY_ON .and. TURN_ATTENUATION_ON) then
    call exit_MPI('cannot have anisotropy and attenuation both turned on in current version')
  end if
!
!----   define coefficients of the Newmark time scheme
!
  deltatover2 = HALF*deltat
  deltatsquareover2 = HALF*deltat*deltat

!---- define actual location of source and receivers
  if(source_type == 1) then
! collocated force source
    call locate_source_force(coord,ibool,npoin,nspec,x_source,z_source,source_type, &
      ix_source,iz_source,ispec_selected_source,iglob_source,is_proc_source,nb_proc_source)

! check that acoustic source is not exactly on the free surface because pressure is zero there
    if ( is_proc_source == 1 ) then
       do ispec_acoustic_surface = 1,nelem_acoustic_surface
          ispec = acoustic_surface(1,ispec_acoustic_surface)
          if( .not. elastic(ispec) .and. ispec == ispec_selected_source ) then
             do j = acoustic_surface(4,ispec_acoustic_surface), acoustic_surface(5,ispec_acoustic_surface)
                do i = acoustic_surface(2,ispec_acoustic_surface), acoustic_surface(3,ispec_acoustic_surface)
                   iglob = ibool(i,j,ispec)
                   if ( iglob_source == iglob ) then
 call exit_MPI('an acoustic source cannot be located exactly on the free surface because pressure is zero there')
                   end if
                end do
             end do
          endif
       enddo
    end if

  else if(source_type == 2) then
! moment-tensor source
     call locate_source_moment_tensor(ibool,coord,nspec,npoin,xigll,zigll,x_source,z_source, &
          ispec_selected_source,is_proc_source,nb_proc_source,nproc,myrank,xi_source,gamma_source,coorg,knods,ngnod,npgeo)

! compute source array for moment-tensor source
    call compute_arrays_source(ispec_selected_source,xi_source,gamma_source,sourcearray, &
               Mxx,Mzz,Mxz,xix,xiz,gammax,gammaz,xigll,zigll,nspec)

  else
    call exit_MPI('incorrect source type')
  endif


! locate receivers in the mesh
  call locate_receivers(ibool,coord,nspec,npoin,xigll,zigll,nrec,nrecloc,recloc,which_proc_receiver,nproc,myrank,&
       st_xval,st_zval,ispec_selected_rec, &
       xi_receiver,gamma_receiver,station_name,network_name,x_source,z_source,coorg,knods,ngnod,npgeo)

! allocate seismogram arrays
  allocate(sisux(NSTEP,nrecloc))
  allocate(sisuz(NSTEP,nrecloc))

! check if acoustic receiver is exactly on the free surface because pressure is zero there
  do ispec_acoustic_surface = 1,nelem_acoustic_surface
     ispec = acoustic_surface(1,ispec_acoustic_surface)
     ixmin = acoustic_surface(2,ispec_acoustic_surface)
     ixmax = acoustic_surface(3,ispec_acoustic_surface)
     izmin = acoustic_surface(4,ispec_acoustic_surface)
     izmax = acoustic_surface(5,ispec_acoustic_surface)
     do irecloc = 1,nrecloc
        irec = recloc(irecloc)
        if(.not. elastic(ispec) .and. ispec == ispec_selected_rec(irec)) then
           if ( (izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==NGLLX .and. &
                gamma_receiver(irec) < -0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==NGLLX .and. &
                gamma_receiver(irec) > 0.99d0) .or.&
                (izmin==1 .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
                xi_receiver(irec) < -0.99d0) .or.&
                (izmin==1 .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                xi_receiver(irec) > 0.99d0) .or.&
                (izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==1 .and. &
                gamma_receiver(irec) < -0.99d0 .and. xi_receiver(irec) < -0.99d0) .or.&
                (izmin==1 .and. izmax==1 .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                gamma_receiver(irec) < -0.99d0 .and. xi_receiver(irec) > 0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
                gamma_receiver(irec) > 0.99d0 .and. xi_receiver(irec) < -0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                gamma_receiver(irec) > 0.99d0 .and. xi_receiver(irec) > 0.99d0) ) then
              if(seismotype == 4) then
call exit_MPI('an acoustic pressure receiver cannot be located exactly on the free surface because pressure is zero there')
              else
                 print *, '**********************************************************************'
                 print *, '*** Warning: acoustic receiver located exactly on the free surface ***'
                 print *, '*** Warning: tangential component will be zero there               ***'
                 print *, '**********************************************************************'
                 print *
              endif
           endif
        endif
     enddo
  enddo

! define and store Lagrange interpolators at all the receivers
  do irec = 1,nrec
    call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
    call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
    hxir_store(irec,:) = hxir(:)
    hgammar_store(irec,:) = hgammar(:)
  enddo

! displacement, velocity, acceleration and inverse of the mass matrix for elastic elements
  if(any_elastic) then
    allocate(displ_elastic(NDIM,npoin))
    allocate(veloc_elastic(NDIM,npoin))
    allocate(accel_elastic(NDIM,npoin))
    allocate(rmass_inverse_elastic(npoin))
  else
! allocate unused arrays with fictitious size
    allocate(displ_elastic(1,1))
    allocate(veloc_elastic(1,1))
    allocate(accel_elastic(1,1))
    allocate(rmass_inverse_elastic(1))
  endif

! potential, its first and second derivative, and inverse of the mass matrix for acoustic elements
  if(any_acoustic) then
    allocate(potential_acoustic(npoin))
    allocate(potential_dot_acoustic(npoin))
    allocate(potential_dot_dot_acoustic(npoin))
    allocate(rmass_inverse_acoustic(npoin))
  else
! allocate unused arrays with fictitious size
    allocate(potential_acoustic(1))
    allocate(potential_dot_acoustic(1))
    allocate(potential_dot_dot_acoustic(1))
    allocate(rmass_inverse_acoustic(1))
  endif

!
!---- build the global mass matrix and invert it once and for all
!
  if(any_elastic) rmass_inverse_elastic(:) = ZERO
  if(any_acoustic) rmass_inverse_acoustic(:) = ZERO
  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
! if external density model
        if(assign_external_model) then
          rhol = rhoext(i,j,ispec)
          cpsquare = vpext(i,j,ispec)**2
        else
          rhol = density(kmato(ispec))
          lambdal_relaxed = elastcoef(1,kmato(ispec))
          mul_relaxed = elastcoef(2,kmato(ispec))
          cpsquare = (lambdal_relaxed + 2.d0*mul_relaxed) / rhol
        endif
! for acoustic medium
        if(elastic(ispec)) then
          rmass_inverse_elastic(iglob) = rmass_inverse_elastic(iglob) + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)
        else
          rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) + wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / cpsquare
        endif
      enddo
    enddo
  enddo

#ifdef USE_MPI
  if ( nproc > 1 ) then
! preparing for MPI communications
    allocate(mask_ispec_inner_outer(nspec))    
    mask_ispec_inner_outer(:) = .false.

    call prepare_assemble_MPI (nspec,ibool, &
          knods, ngnod, &
          npoin, elastic, &
          ninterface, max_interface_size, &
          my_nelmnts_neighbours, my_interfaces, &
          ibool_interfaces_acoustic, ibool_interfaces_elastic, &
          nibool_interfaces_acoustic, nibool_interfaces_elastic, &
          inum_interfaces_acoustic, inum_interfaces_elastic, &
          ninterface_acoustic, ninterface_elastic, &
          mask_ispec_inner_outer &
          )

    nspec_outer = count(mask_ispec_inner_outer)     
    nspec_inner = nspec - nspec_outer          

    allocate(ispec_outer_to_glob(nspec_outer))
    allocate(ispec_inner_to_glob(nspec_inner))

! building of corresponding arrays between inner/outer elements and their global number
    num_ispec_outer = 0
    num_ispec_inner = 0
    do ispec = 1, nspec
      if ( mask_ispec_inner_outer(ispec) ) then
        num_ispec_outer = num_ispec_outer + 1
        ispec_outer_to_glob(num_ispec_outer) = ispec
      else
        num_ispec_inner = num_ispec_inner + 1
        ispec_inner_to_glob(num_ispec_inner) = ispec
           
      endif 
    enddo

  max_ibool_interfaces_size_ac = maxval(nibool_interfaces_acoustic(:))
  max_ibool_interfaces_size_el = NDIM*maxval(nibool_interfaces_elastic(:))
  allocate(tab_requests_send_recv_acoustic(ninterface_acoustic*2))
  allocate(buffer_send_faces_vector_ac(max_ibool_interfaces_size_ac,ninterface_acoustic))
  allocate(buffer_recv_faces_vector_ac(max_ibool_interfaces_size_ac,ninterface_acoustic))
  allocate(tab_requests_send_recv_elastic(ninterface_elastic*2))
  allocate(buffer_send_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic))
  allocate(buffer_recv_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic))

! creating mpi non-blocking persistent communications for acoustic elements
  call create_MPI_req_SEND_RECV_ac( &
     ninterface, ninterface_acoustic, &
     nibool_interfaces_acoustic, &
     my_neighbours, &
     max_ibool_interfaces_size_ac, &
     buffer_send_faces_vector_ac, &
     buffer_recv_faces_vector_ac, &
     tab_requests_send_recv_acoustic, &
     inum_interfaces_acoustic &
     )

! creating mpi non-blocking persistent communications for elastic elements
  call create_MPI_req_SEND_RECV_el( &
     ninterface, ninterface_elastic, &
     nibool_interfaces_elastic, &
     my_neighbours, &
     max_ibool_interfaces_size_el, &
     buffer_send_faces_vector_el, &
     buffer_recv_faces_vector_el, &
     tab_requests_send_recv_elastic, &
     inum_interfaces_elastic &
     )

! assembling the mass matrix
  call assemble_MPI_scalar(myrank,rmass_inverse_acoustic, rmass_inverse_elastic,npoin, &
     ninterface, max_interface_size, max_ibool_interfaces_size_ac, max_ibool_interfaces_size_el, &
     ibool_interfaces_acoustic,ibool_interfaces_elastic, nibool_interfaces_acoustic,nibool_interfaces_elastic, my_neighbours)

  else 
    ninterface_acoustic = 0
    ninterface_elastic = 0
    
    nspec_outer = 0
    nspec_inner = nspec
     
    allocate(ispec_inner_to_glob(nspec_inner))
    do ispec = 1, nspec
      ispec_inner_to_glob(ispec) = ispec
    enddo

  end if ! end of test on wether there is more than one process ( nproc>1 )

#else
  nspec_outer = 0
  nspec_inner = nspec
     
  allocate(ispec_inner_to_glob(nspec_inner))
  do ispec = 1, nspec
     ispec_inner_to_glob(ispec) = ispec
  enddo

#endif


! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
  if(any_elastic) where(rmass_inverse_elastic <= 0._CUSTOM_REAL) rmass_inverse_elastic = 1._CUSTOM_REAL
  if(any_acoustic) where(rmass_inverse_acoustic <= 0._CUSTOM_REAL) rmass_inverse_acoustic = 1._CUSTOM_REAL

! compute the inverse of the mass matrix
  if(any_elastic) rmass_inverse_elastic(:) = 1._CUSTOM_REAL / rmass_inverse_elastic(:)
  if(any_acoustic) rmass_inverse_acoustic(:) = 1._CUSTOM_REAL / rmass_inverse_acoustic(:)

! check the mesh, stability and number of points per wavelength
  call checkgrid(vpext,vsext,rhoext,density,elastcoef,ibool,kmato,coord,npoin,vpmin,vpmax, &
                 assign_external_model,nspec,numat,deltat,f0,t0,initialfield,time_function_type, &
                 coorg,xinterp,zinterp,shape2D_display,knods,simulation_title,npgeo,pointsdisp,ngnod,any_elastic,myrank,nproc)

! convert receiver angle to radians
  anglerec = anglerec * pi / 180.d0



!
!---- for color images
!

  if(output_color_image) then

! horizontal size of the image
  xmin_color_image_loc = minval(coord(1,:))
  xmax_color_image_loc = maxval(coord(1,:))

! vertical size of the image, slightly increase it to go beyond maximum topography
  zmin_color_image_loc = minval(coord(2,:))
  zmax_color_image_loc = maxval(coord(2,:))

! global values
  xmin_color_image = xmin_color_image_loc
  xmax_color_image = xmax_color_image_loc
  zmin_color_image = zmin_color_image_loc
  zmax_color_image = zmax_color_image_loc
  npgeo_glob = npgeo

#ifdef USE_MPI
  call MPI_ALLREDUCE(xmin_color_image_loc, xmin_color_image, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(xmax_color_image_loc, xmax_color_image, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(zmin_color_image_loc, zmin_color_image, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(zmax_color_image_loc, zmax_color_image, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE(npgeo, npgeo_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ier)

#endif

  zmax_color_image = zmin_color_image + 1.05d0 * (zmax_color_image - zmin_color_image)

! compute number of pixels in the horizontal direction based on typical number
! of spectral elements in a given direction (may give bad results for very elongated models)
  NX_IMAGE_color = nint(sqrt(dble(npgeo_glob))) * (NGLLX-1) + 1

! compute number of pixels in the vertical direction based on ratio of sizes
  NZ_IMAGE_color = nint(NX_IMAGE_color * (zmax_color_image - zmin_color_image) / (xmax_color_image - xmin_color_image))

! convert pixel sizes to even numbers because easier to reduce size, create MPEG movies in postprocessing
  NX_IMAGE_color = 2 * (NX_IMAGE_color / 2)
  NZ_IMAGE_color = 2 * (NZ_IMAGE_color / 2)

! check that image size is not too big
  if ( NX_IMAGE_color > 99999 ) then
      call exit_MPI('output image too big : NX_IMAGE_color > 99999.')
  end if
  if ( NZ_IMAGE_color > 99999 ) then
      call exit_MPI('output image too big : NZ_IMAGE_color > 99999.')
  end if

! allocate an array for image data
  allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color))
  allocate(image_color_vp_display(NX_IMAGE_color,NZ_IMAGE_color))

! allocate an array for the grid point that corresponds to a given image data point
  allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color))
  allocate(copy_iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color))

! create all the pixels
  write(IOUT,*)
  write(IOUT,*) 'locating all the pixels of color images'

  size_pixel_horizontal = (xmax_color_image - xmin_color_image) / dble(NX_IMAGE_color-1)
  size_pixel_vertical = (zmax_color_image - zmin_color_image) / dble(NZ_IMAGE_color-1)

  iglob_image_color(:,:) = -1

! cheking which pixels are inside each elements
  nb_pixel_loc = 0
  do ispec = 1, nspec

     do k = 1, 4
        elmnt_coords(1,k) = coorg(1,knods(k,ispec))
        elmnt_coords(2,k) = coorg(2,knods(k,ispec))

     end do

! avoid working on the whole pixel grid
     min_i = floor(minval((elmnt_coords(1,:) - xmin_color_image))/size_pixel_horizontal) + 1
     max_i = ceiling(maxval((elmnt_coords(1,:) - xmin_color_image))/size_pixel_horizontal) + 1
     min_j = floor(minval((elmnt_coords(2,:) - zmin_color_image))/size_pixel_vertical) + 1
     max_j = ceiling(maxval((elmnt_coords(2,:) - zmin_color_image))/size_pixel_vertical) + 1

     do j = min_j, max_j
        do i = min_i, max_i
           i_coord = (i-1)*size_pixel_horizontal + xmin_color_image
           j_coord = (j-1)*size_pixel_vertical + zmin_color_image

! checking if the pixel is inside the element (must be a convex quadrilateral)
           call is_in_convex_quadrilateral( elmnt_coords, i_coord, j_coord, pixel_is_in)

! if inside, getting the nearest point inside the element
           if ( pixel_is_in ) then
              dist_min_pixel = HUGEVAL
              do k = 1, NGLLX
                 do l = 1, NGLLZ
                    iglob = ibool(k,l,ispec)
                    dist_pixel = (coord(1,iglob)-i_coord)**2 + (coord(2,iglob)-j_coord)**2
                    if (dist_pixel < dist_min_pixel) then
                       dist_min_pixel = dist_pixel
                       iglob_image_color(i,j) = iglob

                    end if

                 end do
              end do
              if ( dist_min_pixel >= HUGEVAL ) then
                 call exit_MPI('Error in detecting pixel for color image')

              end if
              nb_pixel_loc = nb_pixel_loc + 1

           end if

        end do
     end do
  end do

! creating and filling array num_pixel_loc with the positions of each colored
! pixel owned by the local process (useful for parallel jobs)
  allocate(num_pixel_loc(nb_pixel_loc))

  nb_pixel_loc = 0
  do i = 1, NX_IMAGE_color
     do j = 1, NZ_IMAGE_color
        if ( iglob_image_color(i,j) /= -1 ) then
           nb_pixel_loc = nb_pixel_loc + 1
           num_pixel_loc(nb_pixel_loc) = (j-1)*NX_IMAGE_color + i

        end if

     end do
  end do



! filling array iglob_image_color, containing info on which process owns which pixels.
#ifdef USE_MPI
  allocate(nb_pixel_per_proc(nproc))

  call MPI_GATHER( nb_pixel_loc, 1, MPI_INTEGER, nb_pixel_per_proc(1), 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)

  if ( myrank == 0 ) then
     allocate(num_pixel_recv(maxval(nb_pixel_per_proc(:)),nproc))
     allocate(data_pixel_recv(maxval(nb_pixel_per_proc(:))))
  end if

  allocate(data_pixel_send(nb_pixel_loc))
  if ( nproc > 1 ) then
     if ( myrank == 0 ) then

        do iproc = 1, nproc-1

           call MPI_RECV(num_pixel_recv(1,iproc+1),nb_pixel_per_proc(iproc+1), MPI_INTEGER, &
                iproc, 42, MPI_COMM_WORLD, request_mpi_status, ier)
           do k = 1, nb_pixel_per_proc(iproc+1)
              j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
              i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color
              iglob_image_color(i,j) = iproc

           end do
        end do

     else
        call MPI_SEND(num_pixel_loc(1),nb_pixel_loc,MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)

     end if

  end if
#else
   allocate(nb_pixel_per_proc(1))
   deallocate(nb_pixel_per_proc)
   allocate(num_pixel_recv(1,1))
   deallocate(num_pixel_recv)
   allocate(data_pixel_recv(1))
   deallocate(data_pixel_recv)
   allocate(data_pixel_send(1))
   deallocate(data_pixel_send)
#endif

  write(IOUT,*) 'done locating all the pixels of color images'

  endif

!
!---- initialize seismograms
!
  sisux = ZERO
  sisuz = ZERO

  cosrot = cos(anglerec)
  sinrot = sin(anglerec)

! initialize arrays to zero
  displ_elastic = ZERO
  veloc_elastic = ZERO
  accel_elastic = ZERO

  potential_acoustic = ZERO
  potential_dot_acoustic = ZERO
  potential_dot_dot_acoustic = ZERO

!
!----  read initial fields from external file if needed
!
  if(initialfield) then
    write(IOUT,*)
    write(IOUT,*) 'Reading initial fields from external file...'
    write(IOUT,*)
    if(any_acoustic) call exit_MPI('initial field currently implemented for purely elastic simulation only')
    open(unit=55,file='OUTPUT_FILES/wavefields.txt',status='unknown')
    read(55,*) nbpoin
    if(nbpoin /= npoin) call exit_MPI('Wrong number of points in input file')
    allocate(displread(NDIM))
    allocate(velocread(NDIM))
    allocate(accelread(NDIM))
    do n = 1,npoin
      read(55,*) inump, (displread(i), i=1,NDIM), &
          (velocread(i), i=1,NDIM), (accelread(i), i=1,NDIM)
      if(inump<1 .or. inump>npoin) call exit_MPI('Wrong point number')
      displ_elastic(:,inump) = displread
      veloc_elastic(:,inump) = velocread
      accel_elastic(:,inump) = accelread
    enddo
    deallocate(displread)
    deallocate(velocread)
    deallocate(accelread)
    close(55)
    write(IOUT,*) 'Max norm of initial elastic displacement = ',maxval(sqrt(displ_elastic(1,:)**2 + displ_elastic(2,:)**2))
  endif

  deltatsquare = deltat * deltat
  deltatcube = deltatsquare * deltat
  deltatfourth = deltatsquare * deltatsquare

  twelvedeltat = 12.d0 * deltat
  fourdeltatsquare = 4.d0 * deltatsquare

! compute the source time function and store it in a text file
  if(.not. initialfield) then

    allocate(source_time_function(NSTEP))

      if ( myrank == 0 ) then
    write(IOUT,*)
    write(IOUT,*) 'Saving the source time function in a text file...'
    write(IOUT,*)
    open(unit=55,file='OUTPUT_FILES/source.txt',status='unknown')
      endif

! loop on all the time steps
    do it = 1,NSTEP

! compute current time
      time = (it-1)*deltat

! Ricker (second derivative of a Gaussian) source time function
      if(time_function_type == 1) then
        source_time_function(it) = - factor * (ONE-TWO*a*(time-t0)**2) * exp(-a*(time-t0)**2)

! first derivative of a Gaussian source time function
      else if(time_function_type == 2) then
        source_time_function(it) = - factor * TWO*a*(time-t0) * exp(-a*(time-t0)**2)

! Gaussian or Dirac (we use a very thin Gaussian instead) source time function
      else if(time_function_type == 3 .or. time_function_type == 4) then
        source_time_function(it) = factor * exp(-a*(time-t0)**2)

! Heaviside source time function (we use a very thin error function instead)
      else if(time_function_type == 5) then
        hdur = 1.d0 / f0
        hdur_gauss = hdur * 5.d0 / 3.d0
        source_time_function(it) = factor * 0.5d0*(1.0d0 + erf(SOURCE_DECAY_MIMIC_TRIANGLE*(time-t0)/hdur_gauss))

      else
        call exit_MPI('unknown source time function')
      endif

! output absolute time in third column, in case user wants to check it as well
      if ( myrank == 0 ) then
      write(55,*) sngl(time),real(source_time_function(it),4),sngl(time-t0)
      endif
   enddo

      if ( myrank == 0 ) then
    close(55)
      endif

! nb_proc_source is the number of processes that own the source (the nearest point). It can be greater 
! than one if the nearest point is on the interface between several partitions with an explosive source.
! since source contribution is linear, the source_time_function is cut down by that number (it would have been similar 
! if we just had elected one of those processes).
    source_time_function(:) = source_time_function(:) / nb_proc_source

  else

    allocate(source_time_function(1))

 endif


! determine if coupled fluid-solid simulation
  coupled_acoustic_elastic = any_acoustic .and. any_elastic

! fluid/solid edge detection
! the two elements (fluid and solid) forming an edge are already known (computed in meshfem2D),
! the common nodes forming the edge are computed here
  if(coupled_acoustic_elastic) then
    print *
    print *,'Mixed acoustic/elastic simulation'
    print *
    print *,'Beginning of fluid/solid edge detection (slow algorithm for now, will be improved later)'

! define the edges of a given element
    i_begin(IBOTTOM) = 1
    j_begin(IBOTTOM) = 1
    i_end(IBOTTOM) = NGLLX
    j_end(IBOTTOM) = 1

    i_begin(IRIGHT) = NGLLX
    j_begin(IRIGHT) = 1
    i_end(IRIGHT) = NGLLX
    j_end(IRIGHT) = NGLLZ

    i_begin(ITOP) = NGLLX
    j_begin(ITOP) = NGLLZ
    i_end(ITOP) = 1
    j_end(ITOP) = NGLLZ

    i_begin(ILEFT) = 1
    j_begin(ILEFT) = NGLLZ
    i_end(ILEFT) = 1
    j_end(ILEFT) = 1

! define i and j points for each edge
    do ipoin1D = 1,NGLLX

      ivalue(ipoin1D,IBOTTOM) = ipoin1D
      ivalue_inverse(ipoin1D,IBOTTOM) = NGLLX - ipoin1D + 1
      jvalue(ipoin1D,IBOTTOM) = 1
      jvalue_inverse(ipoin1D,IBOTTOM) = 1

      ivalue(ipoin1D,IRIGHT) = NGLLX
      ivalue_inverse(ipoin1D,IRIGHT) = NGLLX
      jvalue(ipoin1D,IRIGHT) = ipoin1D
      jvalue_inverse(ipoin1D,IRIGHT) = NGLLZ - ipoin1D + 1

      ivalue(ipoin1D,ITOP) = NGLLX - ipoin1D + 1
      ivalue_inverse(ipoin1D,ITOP) = ipoin1D
      jvalue(ipoin1D,ITOP) = NGLLZ
      jvalue_inverse(ipoin1D,ITOP) = NGLLZ

      ivalue(ipoin1D,ILEFT) = 1
      ivalue_inverse(ipoin1D,ILEFT) = 1
      jvalue(ipoin1D,ILEFT) = NGLLZ - ipoin1D + 1
      jvalue_inverse(ipoin1D,ILEFT) = ipoin1D

    enddo


    do inum = 1, num_fluid_solid_edges
       ispec_acoustic =  fluid_solid_acoustic_ispec(inum)
       ispec_elastic =  fluid_solid_elastic_ispec(inum)

! one element must be acoustic and the other must be elastic
        if(ispec_acoustic /= ispec_elastic .and. .not. elastic(ispec_acoustic) .and. elastic(ispec_elastic)) then

! loop on the four edges of the two elements
          do iedge_acoustic = 1,NEDGES
            do iedge_elastic = 1,NEDGES

! store the matching topology if the two edges match in inverse order
              if(ibool(i_begin(iedge_acoustic),j_begin(iedge_acoustic),ispec_acoustic) == &
                   ibool(i_end(iedge_elastic),j_end(iedge_elastic),ispec_elastic) .and. &
                   ibool(i_end(iedge_acoustic),j_end(iedge_acoustic),ispec_acoustic) == &
                   ibool(i_begin(iedge_elastic),j_begin(iedge_elastic),ispec_elastic)) then
                 fluid_solid_acoustic_iedge(inum) = iedge_acoustic
                 fluid_solid_elastic_iedge(inum) = iedge_elastic
!                  print *,'edge ',iedge_acoustic,' of acoustic element ',ispec_acoustic, &
!                          ' is in contact with edge ',iedge_elastic,' of elastic element ',ispec_elastic
                endif

             enddo
          enddo

       endif

    enddo


! make sure fluid/solid matching has been perfectly detected: check that the grid points
! have the same physical coordinates
! loop on all the coupling edges

    print *,'Checking fluid/solid edge topology...'

    do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
      ispec_acoustic = fluid_solid_acoustic_ispec(inum)
      iedge_acoustic = fluid_solid_acoustic_iedge(inum)

! get the corresponding edge of the elastic element
      ispec_elastic = fluid_solid_elastic_ispec(inum)
      iedge_elastic = fluid_solid_elastic_iedge(inum)

! implement 1D coupling along the edge
      do ipoin1D = 1,NGLLX

! get point values for the elastic side, which matches our side in the inverse direction
        i = ivalue_inverse(ipoin1D,iedge_elastic)
        j = jvalue_inverse(ipoin1D,iedge_elastic)
        iglob = ibool(i,j,ispec_elastic)

! get point values for the acoustic side
        i = ivalue(ipoin1D,iedge_acoustic)
        j = jvalue(ipoin1D,iedge_acoustic)
        iglob2 = ibool(i,j,ispec_acoustic)

! if distance between the two points is not negligible, there is an error, since it should be zero
        if(sqrt((coord(1,iglob) - coord(1,iglob2))**2 + (coord(2,iglob) - coord(2,iglob2))**2) > TINYVAL) &
            call exit_MPI( 'error in fluid/solid coupling buffer')

      enddo

    enddo

    print *,'End of fluid/solid edge detection'
    print *

  else



  endif

! exclude common points between acoustic absorbing edges and acoustic/elastic matching interfaces
  if(coupled_acoustic_elastic .and. anyabs) then

    print *,'excluding common points between acoustic absorbing edges and acoustic/elastic matching interfaces, if any'

! loop on all the absorbing elements
    do ispecabs = 1,nelemabs

      ispec = numabs(ispecabs)

! loop on all the coupling edges
      do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_solid_acoustic_ispec(inum)
        iedge_acoustic = fluid_solid_acoustic_iedge(inum)

! if acoustic absorbing element and acoustic/elastic coupled element is the same
        if(ispec_acoustic == ispec) then

           if(iedge_acoustic == IBOTTOM) then
            jbegin_left(ispecabs) = 2
            jbegin_right(ispecabs) = 2
          endif

          if(iedge_acoustic == ITOP) then
            jend_left(ispecabs) = NGLLZ - 1
            jend_right(ispecabs) = NGLLZ - 1
          endif

          if(iedge_acoustic == ILEFT) then
            ibegin_bottom(ispecabs) = 2
            ibegin_top(ispecabs) = 2
          endif

          if(iedge_acoustic == IRIGHT) then
            iend_bottom(ispecabs) = NGLLX - 1
            iend_top(ispecabs) = NGLLX - 1
          endif

        endif

      enddo

    enddo

  endif

#ifdef USE_MPI
  if(OUTPUT_ENERGY) stop 'energy calculation only serial right now, should add an MPI_REDUCE in parallel'
#endif
! open the file in which we will store the energy curve
  if(OUTPUT_ENERGY) open(unit=IENERGY,file='energy.gnu',status='unknown')

!
!----          s t a r t   t i m e   i t e r a t i o n s
!

  write(IOUT,400)

! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_start = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
               60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

! to display the P-velocity model in background on color images
  allocate(vp_display(npoin))
  do ispec = 1,nspec
! get relaxed elastic parameters of current spectral element
    rhol = density(kmato(ispec))
    lambdal_relaxed = elastcoef(1,kmato(ispec))
    mul_relaxed = elastcoef(2,kmato(ispec))
    do j = 1,NGLLZ
      do i = 1,NGLLX
!--- if external medium, get elastic parameters of current grid point
          if(assign_external_model) then
            vp_display(ibool(i,j,ispec)) = vpext(i,j,ispec)
          else
            cpsquare = (lambdal_relaxed + 2.d0*mul_relaxed) / rhol
            vp_display(ibool(i,j,ispec)) = sqrt(cpsquare)
          endif
      enddo
    enddo
  enddo

! getting velocity for each local pixels
  image_color_vp_display(:,:) = 0.d0
      
  do k = 1, nb_pixel_loc
    j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
    i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color
    image_color_vp_display(i,j) = vp_display(iglob_image_color(i,j))

  enddo

! assembling array image_color_vp_display on process zero for color output
#ifdef USE_MPI
  if ( nproc > 1 ) then
    if ( myrank == 0 ) then
      do iproc = 1, nproc-1
        call MPI_RECV(data_pixel_recv(1),nb_pixel_per_proc(iproc+1), MPI_DOUBLE_PRECISION, &
                iproc, 43, MPI_COMM_WORLD, request_mpi_status, ier)

        do k = 1, nb_pixel_per_proc(iproc+1)
          j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
          i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color
          image_color_vp_display(i,j) = data_pixel_recv(k)

        enddo
      enddo

    else
      do k = 1, nb_pixel_loc
        j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
        i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color
        data_pixel_send(k) = vp_display(iglob_image_color(i,j))

      enddo

      call MPI_SEND(data_pixel_send(1),nb_pixel_loc,MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)

    endif
  endif

#endif

! *********************************************************
! ************* MAIN LOOP OVER THE TIME STEPS *************
! *********************************************************

  do it = 1,NSTEP

! compute current time
    time = (it-1)*deltat

! update displacement using finite-difference time scheme (Newmark)
    if(any_elastic) then
      displ_elastic = displ_elastic + deltat*veloc_elastic + deltatsquareover2*accel_elastic
      veloc_elastic = veloc_elastic + deltatover2*accel_elastic
      accel_elastic = ZERO
    endif

    if(any_acoustic) then

      potential_acoustic = potential_acoustic + deltat*potential_dot_acoustic + deltatsquareover2*potential_dot_dot_acoustic
      potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic
      potential_dot_dot_acoustic = ZERO

! free surface for an acoustic medium
      call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
           potential_acoustic,acoustic_surface, &
           ibool,nelem_acoustic_surface,npoin,nspec)

! *********************************************************
! ************* compute forces for the acoustic elements
! *********************************************************

! first call, computation on outer elements, absorbing conditions and source
    call compute_forces_acoustic(npoin,nspec,nelemabs,numat, &
               iglob_source,ispec_selected_source,is_proc_source,source_type,it,NSTEP,anyabs, &
               assign_external_model,initialfield,ibool,kmato,numabs, &
               elastic,codeabs,potential_dot_dot_acoustic,potential_dot_acoustic, &
               potential_acoustic,density,elastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,source_time_function,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right, &
               nspec_outer, ispec_outer_to_glob, .true. &
               )

 endif ! end of test if any acoustic element

! *********************************************************
! ************* add coupling with the elastic side
! *********************************************************

  if(coupled_acoustic_elastic) then

! loop on all the coupling edges
      do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_solid_acoustic_ispec(inum)
        iedge_acoustic = fluid_solid_acoustic_iedge(inum)

! get the corresponding edge of the elastic element
        ispec_elastic = fluid_solid_elastic_ispec(inum)
        iedge_elastic = fluid_solid_elastic_iedge(inum)

! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

! get point values for the elastic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_elastic)
          j = jvalue_inverse(ipoin1D,iedge_elastic)
          iglob = ibool(i,j,ispec_elastic)

          displ_x = displ_elastic(1,iglob)
          displ_z = displ_elastic(2,iglob)

! get point values for the acoustic side
          i = ivalue(ipoin1D,iedge_acoustic)
          j = jvalue(ipoin1D,iedge_acoustic)
          iglob = ibool(i,j,ispec_acoustic)

! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
! Blackwell Science, page 110, equation (4.60).
          if(iedge_acoustic == IBOTTOM .or. iedge_acoustic == ITOP) then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
          else
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
          endif

! compute dot product
          displ_n = displ_x*nx + displ_z*nz

! formulation with generalized potential
          weight = jacobian1D * wxgll(i)

          potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

        enddo

      enddo

   endif

! assembling potential_dot_dot for acoustic elements (send)
#ifdef USE_MPI
  if ( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_ac_start(potential_dot_dot_acoustic,npoin, &
           ninterface, ninterface_acoustic, &
           inum_interfaces_acoustic, &
           max_interface_size, max_ibool_interfaces_size_ac,&
           ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
           tab_requests_send_recv_acoustic, &
           buffer_send_faces_vector_ac &
           )
  endif
#endif

! second call, computation on inner elements
  if(any_acoustic) then
    call compute_forces_acoustic(npoin,nspec,nelemabs,numat, &
               iglob_source,ispec_selected_source,is_proc_source,source_type,it,NSTEP,anyabs, &
               assign_external_model,initialfield,ibool,kmato,numabs, &
               elastic,codeabs,potential_dot_dot_acoustic,potential_dot_acoustic, &
               potential_acoustic,density,elastcoef,xix,xiz,gammax,gammaz,jacobian, &
               vpext,source_time_function,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll, &
               ibegin_bottom,iend_bottom,ibegin_top,iend_top, &
               jbegin_left,jend_left,jbegin_right,jend_right, &
	       nspec_inner, ispec_inner_to_glob, .false. &	      
	       )
   endif

! assembling potential_dot_dot for acoustic elements (receive)
#ifdef USE_MPI
  if ( nproc > 1 .and. any_acoustic .and. ninterface_acoustic > 0) then
    call assemble_MPI_vector_ac_wait(potential_dot_dot_acoustic,npoin, &
           ninterface, ninterface_acoustic, &
           inum_interfaces_acoustic, &
           max_interface_size, max_ibool_interfaces_size_ac,&
           ibool_interfaces_acoustic, nibool_interfaces_acoustic, &
           tab_requests_send_recv_acoustic, &
           buffer_recv_faces_vector_ac &
           )
  endif
#endif


! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

  if(any_acoustic) then

    potential_dot_dot_acoustic = potential_dot_dot_acoustic * rmass_inverse_acoustic
    potential_dot_acoustic = potential_dot_acoustic + deltatover2*potential_dot_dot_acoustic

! free surface for an acoustic medium
    call enforce_acoustic_free_surface(potential_dot_dot_acoustic,potential_dot_acoustic, &
                potential_acoustic,acoustic_surface, &
                ibool,nelem_acoustic_surface,npoin,nspec)
  endif

! *********************************************************
! ************* main solver for the elastic elements
! *********************************************************

! first call, computation on outer elements, absorbing conditions and source
 if(any_elastic) &
    call compute_forces_elastic(npoin,nspec,nelemabs,numat,iglob_source, &
               ispec_selected_source,is_proc_source,source_type,it,NSTEP,anyabs,assign_external_model, &
               initialfield,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,angleforce,deltatcube, &
               deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,elastic,codeabs, &
               accel_elastic,veloc_elastic,displ_elastic,density,elastcoef,xix,xiz,gammax,gammaz, &
               jacobian,vpext,vsext,rhoext,source_time_function,sourcearray, &
               e1,e11,e13,dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
               dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll,inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2, &
               nspec_outer, ispec_outer_to_glob, .true. &
               )

! *********************************************************
! ************* add coupling with the acoustic side
! *********************************************************

    if(coupled_acoustic_elastic) then

! loop on all the coupling edges
      do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
        ispec_acoustic = fluid_solid_acoustic_ispec(inum)
        iedge_acoustic = fluid_solid_acoustic_iedge(inum)

! get the corresponding edge of the elastic element
        ispec_elastic = fluid_solid_elastic_ispec(inum)
        iedge_elastic = fluid_solid_elastic_iedge(inum)

! implement 1D coupling along the edge
        do ipoin1D = 1,NGLLX

! get point values for the acoustic side, which matches our side in the inverse direction
          i = ivalue_inverse(ipoin1D,iedge_acoustic)
          j = jvalue_inverse(ipoin1D,iedge_acoustic)
          iglob = ibool(i,j,ispec_acoustic)

! get density of the fluid, depending if external density model
          if(assign_external_model) then
            rhol = rhoext(i,j,ispec_acoustic)
          else
            rhol = density(kmato(ispec_acoustic))
          endif

! compute pressure on the fluid/solid edge
          pressure = - rhol * potential_dot_dot_acoustic(iglob)

! get point values for the elastic side
          i = ivalue(ipoin1D,iedge_elastic)
          j = jvalue(ipoin1D,iedge_elastic)
          iglob = ibool(i,j,ispec_elastic)

! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
! Blackwell Science, page 110, equation (4.60).
          if(iedge_acoustic == IBOTTOM .or. iedge_acoustic == ITOP) then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            nx = + zxi / jacobian1D
            nz = - xxi / jacobian1D
          else
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            nx = + zgamma / jacobian1D
            nz = - xgamma / jacobian1D
          endif

! formulation with generalized potential
          weight = jacobian1D * wxgll(i)

          accel_elastic(1,iglob) = accel_elastic(1,iglob) + weight*nx*pressure
          accel_elastic(2,iglob) = accel_elastic(2,iglob) + weight*nz*pressure

        enddo

      enddo

    endif

! assembling accel_elastic for elastic elements (send)
#ifdef USE_MPI
 if ( nproc > 1 .and. any_elastic .and. ninterface_elastic > 0) then
    call assemble_MPI_vector_el_start(accel_elastic,npoin, &
     ninterface, ninterface_elastic, &
     inum_interfaces_elastic, &
     max_interface_size, max_ibool_interfaces_size_el,&
     ibool_interfaces_elastic, nibool_interfaces_elastic, &
     tab_requests_send_recv_elastic, &
     buffer_send_faces_vector_el &
     )
  endif
#endif

! second call, computation on inner elements and update of 
  if(any_elastic) &
    call compute_forces_elastic(npoin,nspec,nelemabs,numat,iglob_source, &
               ispec_selected_source,is_proc_source,source_type,it,NSTEP,anyabs,assign_external_model, &
               initialfield,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,angleforce,deltatcube, &
               deltatfourth,twelvedeltat,fourdeltatsquare,ibool,kmato,numabs,elastic,codeabs, &
               accel_elastic,veloc_elastic,displ_elastic,density,elastcoef,xix,xiz,gammax,gammaz, &
               jacobian,vpext,vsext,rhoext,source_time_function,sourcearray, &
               e1,e11,e13,dux_dxl_n,duz_dzl_n,duz_dxl_n,dux_dzl_n, &
               dux_dxl_np1,duz_dzl_np1,duz_dxl_np1,dux_dzl_np1,hprime_xx,hprimewgll_xx, &
               hprime_zz,hprimewgll_zz,wxgll,wzgll,inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2,Mu_nu1,Mu_nu2, &
	       nspec_inner, ispec_inner_to_glob, .false. &
	       )

! assembling accel_elastic for elastic elements (receive)
#ifdef USE_MPI
 if ( nproc > 1 .and. any_elastic .and. ninterface_elastic > 0) then
    call assemble_MPI_vector_el_wait(accel_elastic,npoin, &
     ninterface, ninterface_elastic, &
     inum_interfaces_elastic, &
     max_interface_size, max_ibool_interfaces_size_el,&
     ibool_interfaces_elastic, nibool_interfaces_elastic, &
     tab_requests_send_recv_elastic, &
     buffer_recv_faces_vector_el &
     )
  end if
#endif


! ************************************************************************************
! ************* multiply by the inverse of the mass matrix and update velocity
! ************************************************************************************

  if(any_elastic) then
    accel_elastic(1,:) = accel_elastic(1,:) * rmass_inverse_elastic
    accel_elastic(2,:) = accel_elastic(2,:) * rmass_inverse_elastic
    veloc_elastic = veloc_elastic + deltatover2*accel_elastic
  endif

!----  compute kinetic and potential energy
!
  if(OUTPUT_ENERGY) &
     call compute_energy(displ_elastic,veloc_elastic, &
         xix,xiz,gammax,gammaz,jacobian,ibool,elastic,hprime_xx,hprime_zz, &
         nspec,npoin,assign_external_model,it,deltat,t0,kmato,elastcoef,density, &
         vpext,vsext,rhoext,wxgll,wzgll,numat, &
         pressure_element,vector_field_element,e1,e11, &
         potential_dot_acoustic,potential_dot_dot_acoustic,TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,Mu_nu1,Mu_nu2)

!----  display time step and max of norm of displacement
  if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5) then

    write(IOUT,*)
    if(time >= 1.d-3 .and. time < 1000.d0) then
      write(IOUT,"('Time step number ',i7,'   t = ',f9.4,' s')") it,time
    else
      write(IOUT,"('Time step number ',i7,'   t = ',1pe12.6,' s')") it,time
    endif

    if(any_elastic) then
      displnorm_all = maxval(sqrt(displ_elastic(1,:)**2 + displ_elastic(2,:)**2))
      write(IOUT,*) 'Max norm of vector field in solid = ',displnorm_all
! check stability of the code in solid, exit if unstable
      if(displnorm_all > STABILITY_THRESHOLD) call exit_MPI('code became unstable and blew up in solid')
    endif

    if(any_acoustic) then
      displnorm_all = maxval(abs(potential_acoustic(:)))
      write(IOUT,*) 'Max absolute value of scalar field in fluid = ',displnorm_all
! check stability of the code in fluid, exit if unstable
      if(displnorm_all > STABILITY_THRESHOLD) call exit_MPI('code became unstable and blew up in fluid')
    endif

    write(IOUT,*)
  endif

! loop on all the receivers to compute and store the seismograms
  do irecloc = 1,nrecloc

     irec = recloc(irecloc)

    ispec = ispec_selected_rec(irec)

! compute pressure in this element if needed
    if(seismotype == 4) then

      call compute_pressure_one_element(pressure_element,potential_dot_dot_acoustic,displ_elastic,elastic, &
            xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
            numat,kmato,density,elastcoef,vpext,vsext,rhoext,ispec,e1,e11, &
            TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,Mu_nu1,Mu_nu2)

    else if(.not. elastic(ispec)) then

! for acoustic medium, compute vector field from gradient of potential for seismograms
      if(seismotype == 1) then
        call compute_vector_one_element(vector_field_element,potential_acoustic,displ_elastic,elastic, &
               xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,ispec)
      else if(seismotype == 2) then
        call compute_vector_one_element(vector_field_element,potential_dot_acoustic,veloc_elastic,elastic, &
               xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,ispec)
      else if(seismotype == 3) then
        call compute_vector_one_element(vector_field_element,potential_dot_dot_acoustic,accel_elastic,elastic, &
               xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,ispec)
      endif

    endif

! perform the general interpolation using Lagrange polynomials
    valux = ZERO
    valuz = ZERO

    do j = 1,NGLLZ
      do i = 1,NGLLX

        iglob = ibool(i,j,ispec)

        hlagrange = hxir_store(irec,i)*hgammar_store(irec,j)

        if(seismotype == 4) then

          dxd = pressure_element(i,j)
          dzd = ZERO

        else if(.not. elastic(ispec)) then

          dxd = vector_field_element(1,i,j)
          dzd = vector_field_element(2,i,j)

        else if(seismotype == 1) then

          dxd = displ_elastic(1,iglob)
          dzd = displ_elastic(2,iglob)

        else if(seismotype == 2) then

          dxd = veloc_elastic(1,iglob)
          dzd = veloc_elastic(2,iglob)

        else if(seismotype == 3) then

          dxd = accel_elastic(1,iglob)
          dzd = accel_elastic(2,iglob)

        endif

! compute interpolated field
        valux = valux + dxd*hlagrange
        valuz = valuz + dzd*hlagrange

      enddo
    enddo

! rotate seismogram components if needed, except if recording pressure, which is a scalar
    if(seismotype /= 4) then
      sisux(it,irecloc) =   cosrot*valux + sinrot*valuz
      sisuz(it,irecloc) = - sinrot*valux + cosrot*valuz
    else
      sisux(it,irecloc) = valux
      sisuz(it,irecloc) = ZERO
    endif

  enddo

!
!----  display results at given time steps
!
  if(mod(it,NTSTEP_BETWEEN_OUTPUT_INFO) == 0 .or. it == 5 .or. it == NSTEP) then

!
!----  PostScript display
!
  if(output_postscript_snapshot) then

  write(IOUT,*) 'Writing PostScript file'

  if(imagetype == 1) then

    write(IOUT,*) 'drawing displacement vector as small arrows...'

    call compute_vector_whole_medium(potential_acoustic,displ_elastic,elastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin)

    call plotpost(vector_field_display,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs, &
          nelem_acoustic_surface, acoustic_edges, &
          simulation_title,npoin,npgeo,vpmin,vpmax,nrec, &
          colors,numbers,subsamp,imagetype,interpol,meshvect,modelvect, &
          boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
          nspec,ngnod,coupled_acoustic_elastic,any_acoustic,plot_lowerleft_corner_only, &
          fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges, &
          myrank, nproc)

  else if(imagetype == 2) then

    write(IOUT,*) 'drawing velocity vector as small arrows...'

    call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,elastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin)

    call plotpost(vector_field_display,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs, &
          nelem_acoustic_surface, acoustic_edges, &
          simulation_title,npoin,npgeo,vpmin,vpmax,nrec, &
          colors,numbers,subsamp,imagetype,interpol,meshvect,modelvect, &
          boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
          nspec,ngnod,coupled_acoustic_elastic,any_acoustic,plot_lowerleft_corner_only, &
          fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges, &
          myrank, nproc)

  else if(imagetype == 3) then

    write(IOUT,*) 'drawing acceleration vector as small arrows...'

    call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,elastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin)

    call plotpost(vector_field_display,coord,vpext,x_source,z_source,st_xval,st_zval, &
          it,deltat,coorg,xinterp,zinterp,shape2D_display, &
          Uxinterp,Uzinterp,flagrange,density,elastcoef,knods,kmato,ibool, &
          numabs,codeabs,anyabs, &
          nelem_acoustic_surface, acoustic_edges, &
          simulation_title,npoin,npgeo,vpmin,vpmax,nrec, &
          colors,numbers,subsamp,imagetype,interpol,meshvect,modelvect, &
          boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,numat,pointsdisp, &
          nspec,ngnod,coupled_acoustic_elastic,any_acoustic,plot_lowerleft_corner_only, &
          fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges, &
          myrank, nproc)

  else if(imagetype == 4) then

    write(IOUT,*) 'cannot draw scalar pressure field as a vector plot, skipping...'

  else
    call exit_MPI('wrong type for snapshots')
  endif

  if(imagetype /= 4) write(IOUT,*) 'PostScript file written'

  endif

!
!----  display color image
!
  if(output_color_image) then

  write(IOUT,*) 'Creating color image of size ',NX_IMAGE_color,' x ',NZ_IMAGE_color

  if(imagetype == 1) then

    write(IOUT,*) 'drawing image of vertical component of displacement vector...'

    call compute_vector_whole_medium(potential_acoustic,displ_elastic,elastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin)

  else if(imagetype == 2) then

    write(IOUT,*) 'drawing image of vertical component of velocity vector...'

    call compute_vector_whole_medium(potential_dot_acoustic,veloc_elastic,elastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin)

  else if(imagetype == 3) then

    write(IOUT,*) 'drawing image of vertical component of acceleration vector...'

    call compute_vector_whole_medium(potential_dot_dot_acoustic,accel_elastic,elastic,vector_field_display, &
          xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin)

  else if(imagetype == 4) then

    write(IOUT,*) 'drawing image of pressure field...'

    call compute_pressure_whole_medium(potential_dot_dot_acoustic,displ_elastic,elastic,vector_field_display, &
         xix,xiz,gammax,gammaz,ibool,hprime_xx,hprime_zz,nspec,npoin,assign_external_model, &
         numat,kmato,density,elastcoef,vpext,vsext,rhoext,e1,e11, &
         TURN_ATTENUATION_ON,TURN_ANISOTROPY_ON,Mu_nu1,Mu_nu2)

  else
    call exit_MPI('wrong type for snapshots')
  endif

  image_color_data(:,:) = 0.d0
      
  do k = 1, nb_pixel_loc
     j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
     i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color
     image_color_data(i,j) = vector_field_display(2,iglob_image_color(i,j))

  end do

! assembling array image_color_data on process zero for color output
#ifdef USE_MPI
  if ( nproc > 1 ) then
     if ( myrank == 0 ) then

        do iproc = 1, nproc-1

           call MPI_RECV(data_pixel_recv(1),nb_pixel_per_proc(iproc+1), MPI_DOUBLE_PRECISION, &
                iproc, 43, MPI_COMM_WORLD, request_mpi_status, ier)

           do k = 1, nb_pixel_per_proc(iproc+1)
              j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
              i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color
              image_color_data(i,j) = data_pixel_recv(k)

           end do
        end do

     else
        do k = 1, nb_pixel_loc
           j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
           i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color
           data_pixel_send(k) = vector_field_display(2,iglob_image_color(i,j))

        end do

        call MPI_SEND(data_pixel_send(1),nb_pixel_loc,MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)

     end if
  end if

#endif


  if ( myrank == 0 ) then
     call create_color_image(image_color_data,iglob_image_color,NX_IMAGE_color,NZ_IMAGE_color,it,cutsnaps,image_color_vp_display)
     write(IOUT,*) 'Color image created'
  endif

  endif

!----  save temporary or final seismograms
  call write_seismograms(sisux,sisuz,station_name,network_name,NSTEP, &
        nrecloc,which_proc_receiver,nrec,myrank,deltat,seismotype,st_xval,it,t0)

! count elapsed wall-clock time
  call date_and_time(datein,timein,zone,time_values)
! time_values(3): day of the month
! time_values(5): hour of the day
! time_values(6): minutes of the hour
! time_values(7): seconds of the minute
! time_values(8): milliseconds of the second
! this fails if we cross the end of the month
  time_end = 86400.d0*time_values(3) + 3600.d0*time_values(5) + &
             60.d0*time_values(6) + time_values(7) + time_values(8) / 1000.d0

! elapsed time since beginning of the simulation
  tCPU = time_end - time_start
  int_tCPU = int(tCPU)
  ihours = int_tCPU / 3600
  iminutes = (int_tCPU - 3600*ihours) / 60
  iseconds = int_tCPU - 3600*ihours - 60*iminutes
  write(*,*) 'Elapsed time in seconds = ',tCPU
  write(*,"(' Elapsed time in hh:mm:ss = ',i4,' h ',i2.2,' m ',i2.2,' s')") ihours,iminutes,iseconds
  write(*,*) 'Mean elapsed time per time step in seconds = ',tCPU/dble(it)
  write(*,*)

  endif

  enddo ! end of the main time loop

!----  close energy file and create a gnuplot script to display it
  if(OUTPUT_ENERGY) then
    close(IENERGY)
    open(unit=IENERGY,file='plotenergy',status='unknown')
    write(IENERGY,*) 'set term postscript landscape color solid "Helvetica" 22'
    write(IENERGY,*) 'set output "energy.ps"'
    write(IENERGY,*) 'set xlabel "Time (s)"'
    write(IENERGY,*) 'set ylabel "Energy (J)"'
    write(IENERGY,*) 'plot "energy.gnu" us 1:4 t ''Total Energy'' w l 1, "energy.gnu" us 1:3 t ''Potential Energy'' w l 2'
    close(IENERGY)
  endif

! print exit banner
  call datim(simulation_title)

!
!----  close output file
!
  if(IOUT /= ISTANDARD_OUTPUT) close(IOUT)

!
!----  formats
!

 400 format(/1x,41('=')/,' =  T i m e  e v o l u t i o n  l o o p  ='/1x,41('=')/)

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
  'Assign external model . . . .(assign_external_model) = ',l6/5x, &
  'Turn anisotropy on or off. . . .(TURN_ANISOTROPY_ON) = ',l6/5x, &
  'Turn attenuation on or off. . .(TURN_ATTENUATION_ON) = ',l6/5x, &
  'Save grid in external file or not. . . .(outputgrid) = ',l6)

 800 format(//1x,'C o n t r o l',/1x,13('='),//5x, &
  'Vector display type . . . . . . . . . . .(imagetype) = ',i6/5x, &
  'Percentage of cut for vector plots . . . .(cutsnaps) = ',f6.2/5x, &
  'Subsampling for velocity model display. . .(subsamp) = ',i6)

 703 format(//' I t e r a t i o n s '/1x,19('='),//5x, &
      'Number of time iterations . . . . .(NSTEP) =',i8,/5x, &
      'Time step increment. . . . . . . .(deltat) =',1pe15.6,/5x, &
      'Total simulation duration . . . . . (ttot) =',1pe15.6)

 107 format(/5x,'--> Isoparametric Spectral Elements <--',//)

 207 format(5x,'Number of spectral elements . . . . .  (nspec) =',i7,/5x, &
               'Number of control nodes per element .  (ngnod) =',i7,/5x, &
               'Number of points in X-direction . . .  (NGLLX) =',i7,/5x, &
               'Number of points in Y-direction . . .  (NGLLZ) =',i7,/5x, &
               'Number of points per element. . .(NGLLX*NGLLZ) =',i7,/5x, &
               'Number of points for display . . .(pointsdisp) =',i7,/5x, &
               'Number of element material sets . . .  (numat) =',i7,/5x, &
               'Number of absorbing elements . . . .(nelemabs) =',i7)

 212 format(//,5x,'Source Type. . . . . . . . . . . . . . = Collocated Force',/5x, &
                  'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
                  'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
                  'Angle from vertical direction (deg). . =',1pe20.10,/5x)

 222 format(//,5x,'Source Type. . . . . . . . . . . . . . = Moment-tensor',/5x, &
                  'X-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Y-position (meters). . . . . . . . . . =',1pe20.10,/5x, &
                  'Fundamental frequency (Hz) . . . . . . =',1pe20.10,/5x, &
                  'Time delay (s) . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Multiplying factor . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mxx. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mzz. . . . . . . . . . . . . . . . . . =',1pe20.10,/5x, &
                  'Mxz. . . . . . . . . . . . . . . . . . =',1pe20.10)

end program specfem2D



subroutine is_in_convex_quadrilateral ( elmnt_coords, x_coord, z_coord, is_in)

  implicit none

  double precision, dimension(2,4)  :: elmnt_coords
  double precision, intent(in)  :: x_coord, z_coord
  logical, intent(out)  :: is_in

  real :: x1, x2, x3, x4, z1, z2, z3, z4
  real  :: normal1, normal2, normal3, normal4


  x1 = elmnt_coords(1,1)
  x2 = elmnt_coords(1,2)
  x3 = elmnt_coords(1,3)
  x4 = elmnt_coords(1,4)
  z1 = elmnt_coords(2,1)
  z2 = elmnt_coords(2,2)
  z3 = elmnt_coords(2,3)
  z4 = elmnt_coords(2,4)

  normal1 = (z_coord-z1) * (x2-x1) - (x_coord-x1) * (z2-z1)
  normal2 = (z_coord-z2) * (x3-x2) - (x_coord-x2) * (z3-z2)
  normal3 = (z_coord-z3) * (x4-x3) - (x_coord-x3) * (z4-z3)
  normal4 = (z_coord-z4) * (x1-x4) - (x_coord-x4) * (z1-z4)

  if ( (normal1 < 0) .or. (normal2 < 0) .or. (normal3 < 0) .or. (normal4 < 0)  ) then
     is_in = .false.
  else
     is_in = .true.
  end if



end subroutine is_in_convex_quadrilateral

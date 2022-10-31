
  program generate_mesh

! Dimitri Komatitsch, Harvard University, USA, around 1998: 2D mesh generator for the global Earth

  implicit none

! dans le code SPECFEM3D_GLOBE on fait :
  ! un doubling sous le Moho
  ! un doubling un peu plus bas que sous la d670
  ! un doubling un peu au dessus du milieu du outer core, entre la CMB et l'ICB

! ici on fait un doubling de moins, et l'un des doublings n'est pas tout a fait au meme endroit :
  ! un doubling juste sous la d670
  ! un doubling un peu au dessus du milieu du outer core, entre la CMB et l'ICB

! in the case of very large meshes, this option can be useful to switch from ASCII to binary for the mesh files
!!#define USE_BINARY_FOR_EXTERNAL_MESH_DATABASE

! 1 = generate the whole mesh   2 = generate only half the mesh
  integer, parameter :: factor_divide_mesh = 2

! multiplication factor for the mesh to be accurate at 1 second of minimum seismic period
!  1 is accurate down to periods of about 45 seconds
!  2 is accurate down to periods of about 22.5 seconds
! 22 is accurate down to periods of about 2 seconds
! 45 is accurate down to periods of about 1 second
! 64 is accurate down to periods of about 0.7 second
! 90 is accurate down to periods of about 0.5 second
  integer, parameter :: multiplication_factor = 2 ! 1 ! 22 ! 45

! output a Gnuplot grid or not (set to false when creating a very large mesh for very short periods)
  logical, parameter :: output_gnuplot_grid = .false. ! .true.

!  The routine uses 9 control nodes defined as follows:
!
!                               4 . . . . 7 . . . . 3
!                               .                   .
!                               .         y         .
!                               .                   .
!                               8         9  x      6
!                               .                   .
!                               .                   .
!                               .                   .
!                               1 . . . . 5 . . . . 2

! nombre de dimensions du probleme physique
  integer, parameter :: NDIM = 2

! nombre de noeuds de controle par element
  integer, parameter :: ngnod = 9

! double du nombre d'elements a la surface
! double du nombre d'elements spectraux dans la direction radiale
  integer, parameter :: nspec_surf_whole_circle      = 512 * multiplication_factor
  integer, parameter :: nspec_rad_670_surf           = 10 * multiplication_factor
  integer, parameter :: nspec_rad_CMB_670            = 56 * multiplication_factor
  integer, parameter :: nspec_rad_doubling_OC_to_CMB = 24 * multiplication_factor
  integer, parameter :: nspec_rad_ICB_to_doubling_OC = 40 * multiplication_factor
  integer, parameter :: nspec_rad_Cube_ICB           = 16 * multiplication_factor

! taille du cube a l'interieur de la graine
  double precision, parameter :: radius_cube = 650000.d0

! to inflate the central cube (set to 0.d0 for a non-inflated cube)
  double precision, parameter :: CENTRAL_CUBE_INFLATE_FACTOR = 0.41d0

! flags for material properties
  integer, parameter :: IREGION_MANTLE_CRUST_ABOVE_d670 = 1
  integer, parameter :: IREGION_MANTLE_BELOW_d670 = 2
  integer, parameter :: IREGION_OUTER_CORE = 3
  integer, parameter :: IREGION_INNER_CORE = 4

! flags for symmetry-axis edges
  integer, parameter :: IBOTTOM = 1
  integer, parameter :: IRIGHT = 2
  integer, parameter :: ITOP = 3
  integer, parameter :: ILEFT = 4

  double precision :: ratio_x, ratio_y, fact_x, fact_y, xi, eta

! prevision exacte du nombre max d'elements et prevision grossiere du nombre max de points
  integer, parameter :: nspec_max_670_surf = (nspec_surf_whole_circle/2)*(nspec_rad_670_surf/2)
  integer, parameter :: nspec_max_CMB_670  = (nspec_surf_whole_circle/4)*(nspec_rad_CMB_670/4 - 1) + &
              (nspec_surf_whole_circle/4)*3
  integer, parameter :: nspec_max_doubling_OC_CMB  = (nspec_surf_whole_circle/4)*(nspec_rad_doubling_OC_to_CMB/4)
  integer, parameter :: nspec_max_ICB_doubling_OC  = (nspec_surf_whole_circle/8)*(nspec_rad_ICB_to_doubling_OC/4 - 1) + &
              (nspec_surf_whole_circle/8)*3
  integer, parameter :: nspec_max_Cube_ICB = (nspec_surf_whole_circle/8)*(nspec_rad_Cube_ICB/4)
  integer, parameter :: nspec_max_inside_cube = (nspec_surf_whole_circle/32)*(nspec_surf_whole_circle/32)

  integer, parameter :: nspec_exact = (nspec_max_670_surf + nspec_max_CMB_670 + nspec_max_doubling_OC_CMB + &
        nspec_max_ICB_doubling_OC + nspec_max_Cube_ICB + nspec_max_inside_cube) / factor_divide_mesh

! prevision surestimee du nombre max de points pour routine Fischer
  integer, parameter :: npoin_max = nspec_exact * ngnod

! topologie des elements du maillage et coordonnees des points de controle
  integer, dimension(ngnod,nspec_exact) :: ibool
  integer, dimension(nspec_exact) :: idoubling
  double precision, dimension(ngnod,nspec_exact) :: xcoord,ycoord

  double precision, dimension(npoin_max) :: xp,yp,work
  integer, dimension(npoin_max) :: iwork
  equivalence(work,iwork)
  integer, dimension(npoin_max) :: loc,ind,ninseg,iglob
  logical, dimension(npoin_max) :: ifseg

! save memory in the case of very large meshes to generate by using two different shapes for the same array
  equivalence(ibool,iglob)

  double precision, parameter :: R_EARTH = 6371000.d0

! values for AK135F_NO_MUD (modified version of AK135)
  double precision, parameter :: R670   = 5711000.d0
  double precision, parameter :: RCMB   = 3479500.d0
  double precision, parameter :: RICB   = 1217500.d0

! values for PREM
!  double precision, parameter :: R670   = 5701000.d0  ! at 670 km depth
!  double precision, parameter :: RCMB   = 3480000.d0  !   2891 km depth
!  double precision, parameter :: RICB   = 1221500.d0

  double precision, parameter :: R_DOUBLING_OUTER_CORE   = RICB + 0.56*(RCMB - RICB)

  double precision, parameter :: PI = 3.141592653589793d0
  double precision, parameter :: PI_OVER_TWO = PI / 2.d0

  double precision zero,one
  parameter(zero=0.d0, one=1.d0)

! to be able to estimate the shortest seismic period that can be resolved by the mesh
  double precision, parameter :: cs_min_in_crust_in_km_per_s = 3.460000d0

  logical :: generate_only_half_the_mesh
  double precision :: delta_theta,size_of_a_surface_element_in_km,radius_interf
  integer :: nspec,ispec,ipoin,ia1,ia2,icentral_cube1,icentral_cube2,ispec_count

! %%%%% pour le maillage elements spectraux %%%%%

!---- zone d670 -> surface

  double precision x1surf(0:nspec_surf_whole_circle)
  double precision y1surf(0:nspec_surf_whole_circle)

  double precision x1bot(0:nspec_surf_whole_circle)
  double precision y1bot(0:nspec_surf_whole_circle)

  double precision x1vol(0:nspec_surf_whole_circle,0:nspec_rad_670_surf)
  double precision y1vol(0:nspec_surf_whole_circle,0:nspec_rad_670_surf)

!---- zone CMB -> d670

  double precision x2surf(0:nspec_surf_whole_circle)
  double precision y2surf(0:nspec_surf_whole_circle)

  double precision x2bot(0:nspec_surf_whole_circle)
  double precision y2bot(0:nspec_surf_whole_circle)

  double precision x2vol(0:nspec_surf_whole_circle,0:nspec_rad_CMB_670)
  double precision y2vol(0:nspec_surf_whole_circle,0:nspec_rad_CMB_670)

!---- zone ICB -> CMB

  double precision x3surf(0:nspec_surf_whole_circle)
  double precision y3surf(0:nspec_surf_whole_circle)

  double precision x3bot(0:nspec_surf_whole_circle)
  double precision y3bot(0:nspec_surf_whole_circle)

  double precision x3vol(0:nspec_surf_whole_circle,0:nspec_rad_ICB_to_doubling_OC)
  double precision y3vol(0:nspec_surf_whole_circle,0:nspec_rad_ICB_to_doubling_OC)

!---- zone Cube -> ICB

  double precision x4surf(0:nspec_surf_whole_circle)
  double precision y4surf(0:nspec_surf_whole_circle)

  double precision x4bot(0:nspec_surf_whole_circle)
  double precision y4bot(0:nspec_surf_whole_circle)

  double precision x4vol(0:nspec_surf_whole_circle,0:nspec_rad_Cube_ICB)
  double precision y4vol(0:nspec_surf_whole_circle,0:nspec_rad_Cube_ICB)

!---- zone interieur du Cube

  double precision x5vol(0:nspec_surf_whole_circle/16,0:nspec_surf_whole_circle/16)
  double precision y5vol(0:nspec_surf_whole_circle/16,0:nspec_surf_whole_circle/16)

! save a lot of memory in the case of a very dense mesh by using the same memory block for the five arrays
  equivalence(x1vol,x2vol,x3vol,x4vol,x5vol)
  equivalence(y1vol,y2vol,y3vol,y4vol,y5vol)

  integer ix,irad,ia
  double precision xicoord,radcoord,thetacoord,xlincoord

! %%%%% pour le maillage elements spectraux %%%%%

! %%%%% pour routine numbering Fischer %%%%%
! pour numerotation locale -> globale
  integer ie,nseg,ioff,iseg,ig,ieoff
  integer ilocnum,i,j,npoin

! mesh tolerance for fast global numbering
  double precision, parameter :: SMALLVALTOL = 1.d-10

! very large and very small values
  double precision, parameter :: HUGEVAL = 1.d+30

  double precision xmaxval,xminval,ymaxval,yminval,xtol,xtypdist

! elements 9 noeuds donc on a seulement la moitie des elements
  delta_theta = 2. * pi / dble(nspec_surf_whole_circle/2)
  size_of_a_surface_element_in_km = delta_theta*R_EARTH/1000.

  print *,'Generate mesh'
  print *
  print *,'Number of elements at the surface of the whole circle = ',nspec_surf_whole_circle
  print *
  print *,'Size of a spectral element at the surface = ',size_of_a_surface_element_in_km,' km'
  print *
  print *,'Since minimum S velocity in the crust is about ',cs_min_in_crust_in_km_per_s,' km/s'
  print *,'and since one spectral element is needed to resolve the shortest wavelength,'
  print *,'this mesh is accurate approximately down to a minimum period of ', &
                   size_of_a_surface_element_in_km/cs_min_in_crust_in_km_per_s,' seconds'
  print *

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!
!---- generation de la grille pour SPECFEM90
!

  print *,'Generating the grid for SPECFEM2D...'
  print *

! generate only half the mesh or the whole mesh
  if (factor_divide_mesh < 1 .or. factor_divide_mesh > 2) stop 'incorrect value of factor_divide_mesh'

  if (factor_divide_mesh == 1) then
    generate_only_half_the_mesh = .false.
  else
    generate_only_half_the_mesh = .true.
  endif

  if (generate_only_half_the_mesh) then
    print *,'generating only half the mesh'
  else
    print *,'generating the whole mesh'
  endif
  print *

! verification de la coherence des parametres entres

! --- pour pouvoir generer des elements 9 noeuds
  if (mod(nspec_rad_670_surf,2) /= 0) stop 'nspec_rad_670_surf should be a multiple of 2'
  if (mod(nspec_rad_CMB_670,4) /= 0) stop 'nspec_rad_CMB_670 should be a multiple of 4'
  if (mod(nspec_rad_ICB_to_doubling_OC,4) /= 0) stop 'nspec_rad_ICB_to_doubling_OC should be a multiple of 4'
  if (mod(nspec_rad_Cube_ICB,4) /= 0) stop 'nspec_rad_Cube_ICB should be a multiple of 4'

! --- pour permettre le deraffinement
  if (mod(nspec_surf_whole_circle,16) /= 0) stop 'nspec_surf_whole_circle should be a multiple of 16'
  if (nspec_rad_CMB_670 < 6) stop 'nspec_rad_CMB_670 should be greater than 6'
  if (nspec_rad_ICB_to_doubling_OC < 6) stop 'nspec_rad_ICB_to_doubling_OC should be greater than 6'

! %%%% zone d670 -> surface %%%%

! generation maillage de la surface

  do ix=0,nspec_surf_whole_circle

  xicoord  = dble(ix)/dble(nspec_surf_whole_circle)

! parcours de l'angle dans le sens negatif
  thetacoord = +pi/2 - 2*pi*xicoord

! coordonnees cartesiennes correspondantes

!---- surface

  radius_interf = R_EARTH

  x1surf(ix) = radius_interf * dcos(thetacoord)
  y1surf(ix) = radius_interf * dsin(thetacoord)

!---- d670

  radius_interf = R670

  x1bot(ix) = radius_interf * dcos(thetacoord)
  y1bot(ix) = radius_interf * dsin(thetacoord)

!---- volume

  do irad=0,nspec_rad_670_surf
      radcoord  = dble(irad)/dble(nspec_rad_670_surf)
      x1vol(ix,irad) = x1surf(ix) * radcoord + x1bot(ix) * (one - radcoord)
      y1vol(ix,irad) = y1surf(ix) * radcoord + y1bot(ix) * (one - radcoord)
  enddo

  enddo

! %%% ecrire la grille dans un tableau pour routine de Fischer
  ispec = 0

! %%% bloc d670 -> surface
  do irad=0,nspec_rad_670_surf-2,2
  do ix=0,nspec_surf_whole_circle/factor_divide_mesh-2,2

      ispec = ispec + 1

      idoubling(ispec) = IREGION_MANTLE_CRUST_ABOVE_d670

      xcoord(1,ispec) = x1vol(ix  ,irad)
      xcoord(5,ispec) = x1vol(ix+1,irad)
      xcoord(2,ispec) = x1vol(ix+2,irad)
      xcoord(8,ispec) = x1vol(ix  ,irad+1)
      xcoord(9,ispec) = x1vol(ix+1,irad+1)
      xcoord(6,ispec) = x1vol(ix+2,irad+1)
      xcoord(4,ispec) = x1vol(ix  ,irad+2)
      xcoord(7,ispec) = x1vol(ix+1,irad+2)
      xcoord(3,ispec) = x1vol(ix+2,irad+2)

      ycoord(1,ispec) = y1vol(ix  ,irad)
      ycoord(5,ispec) = y1vol(ix+1,irad)
      ycoord(2,ispec) = y1vol(ix+2,irad)
      ycoord(8,ispec) = y1vol(ix  ,irad+1)
      ycoord(9,ispec) = y1vol(ix+1,irad+1)
      ycoord(6,ispec) = y1vol(ix+2,irad+1)
      ycoord(4,ispec) = y1vol(ix  ,irad+2)
      ycoord(7,ispec) = y1vol(ix+1,irad+2)
      ycoord(3,ispec) = y1vol(ix+2,irad+2)
  enddo
  enddo

! %%%% zone CMB -> d670 %%%%

! generation maillage de la surface

  do ix=0,nspec_surf_whole_circle

  xicoord  = dble(ix)/dble(nspec_surf_whole_circle)

! parcours de l'angle dans le sens negatif
  thetacoord = +pi/2 - 2*pi*xicoord

! coordonnees cartesiennes correspondantes

!---- d670

  radius_interf = R670

  x2surf(ix) = radius_interf * dcos(thetacoord)
  y2surf(ix) = radius_interf * dsin(thetacoord)

!---- CMB

  radius_interf = RCMB

  x2bot(ix) = radius_interf * dcos(thetacoord)
  y2bot(ix) = radius_interf * dsin(thetacoord)

!---- volume

  do irad=0,nspec_rad_CMB_670
      radcoord  = dble(irad)/dble(nspec_rad_CMB_670)
      x2vol(ix,irad) = x2surf(ix) * radcoord + x2bot(ix) * (one - radcoord)
      y2vol(ix,irad) = y2surf(ix) * radcoord + y2bot(ix) * (one - radcoord)
  enddo

  enddo

! --- bloc principal
  do irad=0,nspec_rad_CMB_670-8,4
  do ix=0,nspec_surf_whole_circle/factor_divide_mesh-4,4

      ispec = ispec + 1

      idoubling(ispec) = IREGION_MANTLE_BELOW_d670

      xcoord(1,ispec) = x2vol(ix  ,irad)
      xcoord(5,ispec) = x2vol(ix+2,irad)
      xcoord(2,ispec) = x2vol(ix+4,irad)
      xcoord(8,ispec) = x2vol(ix  ,irad+2)
      xcoord(9,ispec) = x2vol(ix+2,irad+2)
      xcoord(6,ispec) = x2vol(ix+4,irad+2)
      xcoord(4,ispec) = x2vol(ix  ,irad+4)
      xcoord(7,ispec) = x2vol(ix+2,irad+4)
      xcoord(3,ispec) = x2vol(ix+4,irad+4)

      ycoord(1,ispec) = y2vol(ix  ,irad)
      ycoord(5,ispec) = y2vol(ix+2,irad)
      ycoord(2,ispec) = y2vol(ix+4,irad)
      ycoord(8,ispec) = y2vol(ix  ,irad+2)
      ycoord(9,ispec) = y2vol(ix+2,irad+2)
      ycoord(6,ispec) = y2vol(ix+4,irad+2)
      ycoord(4,ispec) = y2vol(ix  ,irad+4)
      ycoord(7,ispec) = y2vol(ix+2,irad+4)
      ycoord(3,ispec) = y2vol(ix+4,irad+4)
  enddo
  enddo

! --- zone de raccord geometrique conforme
  irad=nspec_rad_CMB_670-4
  do ix=0,nspec_surf_whole_circle/factor_divide_mesh-8,8

      ispec = ispec + 1

      idoubling(ispec) = IREGION_MANTLE_BELOW_d670

      xcoord(1,ispec) = x2vol(ix  ,irad+2)
      xcoord(5,ispec) = x2vol(ix+1,irad+2)
      xcoord(2,ispec) = x2vol(ix+2,irad+2)
      xcoord(8,ispec) = x2vol(ix  ,irad+3)
      xcoord(9,ispec) = x2vol(ix+1,irad+3)
      xcoord(6,ispec) = x2vol(ix+2,irad+3)
      xcoord(4,ispec) = x2vol(ix  ,irad+4)
      xcoord(7,ispec) = x2vol(ix+1,irad+4)
      xcoord(3,ispec) = x2vol(ix+2,irad+4)

      ycoord(1,ispec) = y2vol(ix  ,irad+2)
      ycoord(5,ispec) = y2vol(ix+1,irad+2)
      ycoord(2,ispec) = y2vol(ix+2,irad+2)
      ycoord(8,ispec) = y2vol(ix  ,irad+3)
      ycoord(9,ispec) = y2vol(ix+1,irad+3)
      ycoord(6,ispec) = y2vol(ix+2,irad+3)
      ycoord(4,ispec) = y2vol(ix  ,irad+4)
      ycoord(7,ispec) = y2vol(ix+1,irad+4)
      ycoord(3,ispec) = y2vol(ix+2,irad+4)

      ispec = ispec + 1

      idoubling(ispec) = IREGION_MANTLE_BELOW_d670

      xcoord(1,ispec) = x2vol(ix  ,irad)
      xcoord(5,ispec) = x2vol(ix+2,irad)
      xcoord(2,ispec) = x2vol(ix+4,irad)
      xcoord(8,ispec) = x2vol(ix  ,irad+1)
      xcoord(9,ispec) = (x2vol(ix+1,irad+1) + x2vol(ix+2,irad+1))/2.
      xcoord(6,ispec) = x2vol(ix+3,irad+1)
      xcoord(4,ispec) = x2vol(ix  ,irad+2)
      xcoord(7,ispec) = x2vol(ix+1,irad+2)
      xcoord(3,ispec) = x2vol(ix+2,irad+2)

      ycoord(1,ispec) = y2vol(ix  ,irad)
      ycoord(5,ispec) = y2vol(ix+2,irad)
      ycoord(2,ispec) = y2vol(ix+4,irad)
      ycoord(8,ispec) = y2vol(ix  ,irad+1)
      ycoord(9,ispec) = (y2vol(ix+1,irad+1) + y2vol(ix+2,irad+1))/2.
      ycoord(6,ispec) = y2vol(ix+3,irad+1)
      ycoord(4,ispec) = y2vol(ix  ,irad+2)
      ycoord(7,ispec) = y2vol(ix+1,irad+2)
      ycoord(3,ispec) = y2vol(ix+2,irad+2)

      ispec = ispec + 1

      idoubling(ispec) = IREGION_MANTLE_BELOW_d670

      xcoord(1,ispec) = x2vol(ix+2,irad+2)
      xcoord(5,ispec) = x2vol(ix+3,irad+1)
      xcoord(2,ispec) = x2vol(ix+4,irad)
      xcoord(8,ispec) = x2vol(ix+2,irad+3)
      xcoord(9,ispec) = (x2vol(ix+3,irad+2) + x2vol(ix+3,irad+3))/2.
      xcoord(6,ispec) = x2vol(ix+4,irad+2)
      xcoord(4,ispec) = x2vol(ix+2,irad+4)
      xcoord(7,ispec) = x2vol(ix+3,irad+4)
      xcoord(3,ispec) = x2vol(ix+4,irad+4)

      ycoord(1,ispec) = y2vol(ix+2,irad+2)
      ycoord(5,ispec) = y2vol(ix+3,irad+1)
      ycoord(2,ispec) = y2vol(ix+4,irad)
      ycoord(8,ispec) = y2vol(ix+2,irad+3)
      ycoord(9,ispec) = (y2vol(ix+3,irad+2) + y2vol(ix+3,irad+3))/2.
      ycoord(6,ispec) = y2vol(ix+4,irad+2)
      ycoord(4,ispec) = y2vol(ix+2,irad+4)
      ycoord(7,ispec) = y2vol(ix+3,irad+4)
      ycoord(3,ispec) = y2vol(ix+4,irad+4)

  enddo

! --- zone de raccord geometrique conforme inverse
  irad=nspec_rad_CMB_670-4
  do ix=4,nspec_surf_whole_circle/factor_divide_mesh-4,8

      ispec = ispec + 1

      idoubling(ispec) = IREGION_MANTLE_BELOW_d670

      xcoord(1,ispec) = x2vol(ix+2,irad+2)
      xcoord(5,ispec) = x2vol(ix+3,irad+2)
      xcoord(2,ispec) = x2vol(ix+4,irad+2)
      xcoord(8,ispec) = x2vol(ix+2,irad+3)
      xcoord(9,ispec) = x2vol(ix+3,irad+3)
      xcoord(6,ispec) = x2vol(ix+4,irad+3)
      xcoord(4,ispec) = x2vol(ix+2,irad+4)
      xcoord(7,ispec) = x2vol(ix+3,irad+4)
      xcoord(3,ispec) = x2vol(ix+4,irad+4)

      ycoord(1,ispec) = y2vol(ix+2,irad+2)
      ycoord(5,ispec) = y2vol(ix+3,irad+2)
      ycoord(2,ispec) = y2vol(ix+4,irad+2)
      ycoord(8,ispec) = y2vol(ix+2,irad+3)
      ycoord(9,ispec) = y2vol(ix+3,irad+3)
      ycoord(6,ispec) = y2vol(ix+4,irad+3)
      ycoord(4,ispec) = y2vol(ix+2,irad+4)
      ycoord(7,ispec) = y2vol(ix+3,irad+4)
      ycoord(3,ispec) = y2vol(ix+4,irad+4)

      ispec = ispec + 1

      idoubling(ispec) = IREGION_MANTLE_BELOW_d670

      xcoord(1,ispec) = x2vol(ix  ,irad)
      xcoord(5,ispec) = x2vol(ix+2,irad)
      xcoord(2,ispec) = x2vol(ix+4,irad)
      xcoord(8,ispec) = x2vol(ix+1,irad+1)
      xcoord(9,ispec) = (x2vol(ix+2,irad+1) + x2vol(ix+3,irad+1))/2.
      xcoord(6,ispec) = x2vol(ix+4,irad+1)
      xcoord(4,ispec) = x2vol(ix+2,irad+2)
      xcoord(7,ispec) = x2vol(ix+3,irad+2)
      xcoord(3,ispec) = x2vol(ix+4,irad+2)

      ycoord(1,ispec) = y2vol(ix  ,irad)
      ycoord(5,ispec) = y2vol(ix+2,irad)
      ycoord(2,ispec) = y2vol(ix+4,irad)
      ycoord(8,ispec) = y2vol(ix+1,irad+1)
      ycoord(9,ispec) = (y2vol(ix+2,irad+1) + y2vol(ix+3,irad+1))/2.
      ycoord(6,ispec) = y2vol(ix+4,irad+1)
      ycoord(4,ispec) = y2vol(ix+2,irad+2)
      ycoord(7,ispec) = y2vol(ix+3,irad+2)
      ycoord(3,ispec) = y2vol(ix+4,irad+2)

      ispec = ispec + 1

      idoubling(ispec) = IREGION_MANTLE_BELOW_d670

      xcoord(1,ispec) = x2vol(ix,irad)
      xcoord(5,ispec) = x2vol(ix+1,irad+1)
      xcoord(2,ispec) = x2vol(ix+2,irad+2)
      xcoord(8,ispec) = x2vol(ix,irad+2)
      xcoord(9,ispec) = (x2vol(ix+1,irad+2) + x2vol(ix+1,irad+3))/2.
      xcoord(6,ispec) = x2vol(ix+2,irad+3)
      xcoord(4,ispec) = x2vol(ix,irad+4)
      xcoord(7,ispec) = x2vol(ix+1,irad+4)
      xcoord(3,ispec) = x2vol(ix+2,irad+4)

      ycoord(1,ispec) = y2vol(ix,irad)
      ycoord(5,ispec) = y2vol(ix+1,irad+1)
      ycoord(2,ispec) = y2vol(ix+2,irad+2)
      ycoord(8,ispec) = y2vol(ix,irad+2)
      ycoord(9,ispec) = (y2vol(ix+1,irad+2) + y2vol(ix+1,irad+3))/2.
      ycoord(6,ispec) = y2vol(ix+2,irad+3)
      ycoord(4,ispec) = y2vol(ix,irad+4)
      ycoord(7,ispec) = y2vol(ix+1,irad+4)
      ycoord(3,ispec) = y2vol(ix+2,irad+4)

  enddo

! %%%% zone doubling in the CMB -> CMB %%%%

! generation maillage de la surface

  do ix=0,nspec_surf_whole_circle

  xicoord  = dble(ix)/dble(nspec_surf_whole_circle)

! parcours de l'angle dans le sens negatif
  thetacoord = +pi/2 - 2*pi*xicoord

! coordonnees cartesiennes correspondantes

!---- CMB

  radius_interf = RCMB

  x3surf(ix) = radius_interf * dcos(thetacoord)
  y3surf(ix) = radius_interf * dsin(thetacoord)

!---- doubling in the outer core

  radius_interf = R_DOUBLING_OUTER_CORE

  x3bot(ix) = radius_interf * dcos(thetacoord)
  y3bot(ix) = radius_interf * dsin(thetacoord)

!---- volume

  do irad=0,nspec_rad_doubling_OC_to_CMB
      radcoord  = dble(irad)/dble(nspec_rad_doubling_OC_to_CMB)
      x3vol(ix,irad) = x3surf(ix) * radcoord + x3bot(ix) * (one - radcoord)
      y3vol(ix,irad) = y3surf(ix) * radcoord + y3bot(ix) * (one - radcoord)
  enddo

  enddo

! --- bloc principal
  do irad=0,nspec_rad_doubling_OC_to_CMB-4,4
  do ix=0,nspec_surf_whole_circle/factor_divide_mesh-4,4

      ispec = ispec + 1

      idoubling(ispec) = IREGION_OUTER_CORE

      xcoord(1,ispec) = x3vol(ix  ,irad)
      xcoord(5,ispec) = x3vol(ix+2,irad)
      xcoord(2,ispec) = x3vol(ix+4,irad)
      xcoord(8,ispec) = x3vol(ix  ,irad+2)
      xcoord(9,ispec) = x3vol(ix+2,irad+2)
      xcoord(6,ispec) = x3vol(ix+4,irad+2)
      xcoord(4,ispec) = x3vol(ix  ,irad+4)
      xcoord(7,ispec) = x3vol(ix+2,irad+4)
      xcoord(3,ispec) = x3vol(ix+4,irad+4)

      ycoord(1,ispec) = y3vol(ix  ,irad)
      ycoord(5,ispec) = y3vol(ix+2,irad)
      ycoord(2,ispec) = y3vol(ix+4,irad)
      ycoord(8,ispec) = y3vol(ix  ,irad+2)
      ycoord(9,ispec) = y3vol(ix+2,irad+2)
      ycoord(6,ispec) = y3vol(ix+4,irad+2)
      ycoord(4,ispec) = y3vol(ix  ,irad+4)
      ycoord(7,ispec) = y3vol(ix+2,irad+4)
      ycoord(3,ispec) = y3vol(ix+4,irad+4)
  enddo
  enddo

! %%%% zone ICB -> doubling in the CMB %%%%

! generation maillage de la surface

  do ix=0,nspec_surf_whole_circle

  xicoord  = dble(ix)/dble(nspec_surf_whole_circle)

! parcours de l'angle dans le sens negatif
  thetacoord = +pi/2 - 2*pi*xicoord

! coordonnees cartesiennes correspondantes

!---- doubling in the outer core

  radius_interf = R_DOUBLING_OUTER_CORE

  x3surf(ix) = radius_interf * dcos(thetacoord)
  y3surf(ix) = radius_interf * dsin(thetacoord)

!---- ICB

  radius_interf = RICB

  x3bot(ix) = radius_interf * dcos(thetacoord)
  y3bot(ix) = radius_interf * dsin(thetacoord)

!---- volume

  do irad=0,nspec_rad_ICB_to_doubling_OC
      radcoord  = dble(irad)/dble(nspec_rad_ICB_to_doubling_OC)
      x3vol(ix,irad) = x3surf(ix) * radcoord + x3bot(ix) * (one - radcoord)
      y3vol(ix,irad) = y3surf(ix) * radcoord + y3bot(ix) * (one - radcoord)
  enddo

  enddo

! --- bloc principal
  do irad=0,nspec_rad_ICB_to_doubling_OC-8,4
  do ix=0,nspec_surf_whole_circle/factor_divide_mesh-8,8

      ispec = ispec + 1

      idoubling(ispec) = IREGION_OUTER_CORE

      xcoord(1,ispec) = x3vol(ix  ,irad)
      xcoord(5,ispec) = x3vol(ix+4,irad)
      xcoord(2,ispec) = x3vol(ix+8,irad)
      xcoord(8,ispec) = x3vol(ix  ,irad+2)
      xcoord(9,ispec) = x3vol(ix+4,irad+2)
      xcoord(6,ispec) = x3vol(ix+8,irad+2)
      xcoord(4,ispec) = x3vol(ix  ,irad+4)
      xcoord(7,ispec) = x3vol(ix+4,irad+4)
      xcoord(3,ispec) = x3vol(ix+8,irad+4)

      ycoord(1,ispec) = y3vol(ix  ,irad)
      ycoord(5,ispec) = y3vol(ix+4,irad)
      ycoord(2,ispec) = y3vol(ix+8,irad)
      ycoord(8,ispec) = y3vol(ix  ,irad+2)
      ycoord(9,ispec) = y3vol(ix+4,irad+2)
      ycoord(6,ispec) = y3vol(ix+8,irad+2)
      ycoord(4,ispec) = y3vol(ix  ,irad+4)
      ycoord(7,ispec) = y3vol(ix+4,irad+4)
      ycoord(3,ispec) = y3vol(ix+8,irad+4)
  enddo
  enddo

! --- zone de raccord geometrique conforme
  irad=nspec_rad_ICB_to_doubling_OC-4
  do ix=0,nspec_surf_whole_circle/factor_divide_mesh-16,16

      ispec = ispec + 1

      idoubling(ispec) = IREGION_OUTER_CORE

      xcoord(1,ispec) = x3vol(ix  ,irad+2)
      xcoord(5,ispec) = x3vol(ix+2,irad+2)
      xcoord(2,ispec) = x3vol(ix+4,irad+2)
      xcoord(8,ispec) = x3vol(ix  ,irad+3)
      xcoord(9,ispec) = x3vol(ix+2,irad+3)
      xcoord(6,ispec) = x3vol(ix+4,irad+3)
      xcoord(4,ispec) = x3vol(ix  ,irad+4)
      xcoord(7,ispec) = x3vol(ix+2,irad+4)
      xcoord(3,ispec) = x3vol(ix+4,irad+4)

      ycoord(1,ispec) = y3vol(ix  ,irad+2)
      ycoord(5,ispec) = y3vol(ix+2,irad+2)
      ycoord(2,ispec) = y3vol(ix+4,irad+2)
      ycoord(8,ispec) = y3vol(ix  ,irad+3)
      ycoord(9,ispec) = y3vol(ix+2,irad+3)
      ycoord(6,ispec) = y3vol(ix+4,irad+3)
      ycoord(4,ispec) = y3vol(ix  ,irad+4)
      ycoord(7,ispec) = y3vol(ix+2,irad+4)
      ycoord(3,ispec) = y3vol(ix+4,irad+4)

      ispec = ispec + 1

      idoubling(ispec) = IREGION_OUTER_CORE

      xcoord(1,ispec) = x3vol(ix  ,irad)
      xcoord(5,ispec) = x3vol(ix+4,irad)
      xcoord(2,ispec) = x3vol(ix+8,irad)
      xcoord(8,ispec) = x3vol(ix  ,irad+1)
      xcoord(9,ispec) = (x3vol(ix+2,irad+1) + x3vol(ix+4,irad+1))/2.
      xcoord(6,ispec) = x3vol(ix+6,irad+1)
      xcoord(4,ispec) = x3vol(ix  ,irad+2)
      xcoord(7,ispec) = x3vol(ix+2,irad+2)
      xcoord(3,ispec) = x3vol(ix+4,irad+2)

      ycoord(1,ispec) = y3vol(ix  ,irad)
      ycoord(5,ispec) = y3vol(ix+4,irad)
      ycoord(2,ispec) = y3vol(ix+8,irad)
      ycoord(8,ispec) = y3vol(ix  ,irad+1)
      ycoord(9,ispec) = (y3vol(ix+2,irad+1) + y3vol(ix+4,irad+1))/2.
      ycoord(6,ispec) = y3vol(ix+6,irad+1)
      ycoord(4,ispec) = y3vol(ix  ,irad+2)
      ycoord(7,ispec) = y3vol(ix+2,irad+2)
      ycoord(3,ispec) = y3vol(ix+4,irad+2)

      ispec = ispec + 1

      idoubling(ispec) = IREGION_OUTER_CORE

      xcoord(1,ispec) = x3vol(ix+4,irad+2)
      xcoord(5,ispec) = x3vol(ix+6,irad+1)
      xcoord(2,ispec) = x3vol(ix+8,irad)
      xcoord(8,ispec) = x3vol(ix+4,irad+3)
      xcoord(9,ispec) = (x3vol(ix+6,irad+2) + x3vol(ix+6,irad+3))/2.
      xcoord(6,ispec) = x3vol(ix+8,irad+2)
      xcoord(4,ispec) = x3vol(ix+4,irad+4)
      xcoord(7,ispec) = x3vol(ix+6,irad+4)
      xcoord(3,ispec) = x3vol(ix+8,irad+4)

      ycoord(1,ispec) = y3vol(ix+4,irad+2)
      ycoord(5,ispec) = y3vol(ix+6,irad+1)
      ycoord(2,ispec) = y3vol(ix+8,irad)
      ycoord(8,ispec) = y3vol(ix+4,irad+3)
      ycoord(9,ispec) = (y3vol(ix+6,irad+2) + y3vol(ix+6,irad+3))/2.
      ycoord(6,ispec) = y3vol(ix+8,irad+2)
      ycoord(4,ispec) = y3vol(ix+4,irad+4)
      ycoord(7,ispec) = y3vol(ix+6,irad+4)
      ycoord(3,ispec) = y3vol(ix+8,irad+4)

  enddo

! --- zone de raccord geometrique conforme inverse
  irad=nspec_rad_ICB_to_doubling_OC-4
  do ix=8,nspec_surf_whole_circle/factor_divide_mesh-8,16

      ispec = ispec + 1

      idoubling(ispec) = IREGION_OUTER_CORE

      xcoord(1,ispec) = x3vol(ix+4,irad+2)
      xcoord(5,ispec) = x3vol(ix+6,irad+2)
      xcoord(2,ispec) = x3vol(ix+8,irad+2)
      xcoord(8,ispec) = x3vol(ix+4,irad+3)
      xcoord(9,ispec) = x3vol(ix+6,irad+3)
      xcoord(6,ispec) = x3vol(ix+8,irad+3)
      xcoord(4,ispec) = x3vol(ix+4,irad+4)
      xcoord(7,ispec) = x3vol(ix+6,irad+4)
      xcoord(3,ispec) = x3vol(ix+8,irad+4)

      ycoord(1,ispec) = y3vol(ix+4,irad+2)
      ycoord(5,ispec) = y3vol(ix+6,irad+2)
      ycoord(2,ispec) = y3vol(ix+8,irad+2)
      ycoord(8,ispec) = y3vol(ix+4,irad+3)
      ycoord(9,ispec) = y3vol(ix+6,irad+3)
      ycoord(6,ispec) = y3vol(ix+8,irad+3)
      ycoord(4,ispec) = y3vol(ix+4,irad+4)
      ycoord(7,ispec) = y3vol(ix+6,irad+4)
      ycoord(3,ispec) = y3vol(ix+8,irad+4)

      ispec = ispec + 1

      idoubling(ispec) = IREGION_OUTER_CORE

      xcoord(1,ispec) = x3vol(ix  ,irad)
      xcoord(5,ispec) = x3vol(ix+4,irad)
      xcoord(2,ispec) = x3vol(ix+8,irad)
      xcoord(8,ispec) = x3vol(ix+2,irad+1)
      xcoord(9,ispec) = (x3vol(ix+4,irad+1) + x3vol(ix+6,irad+1))/2.
      xcoord(6,ispec) = x3vol(ix+8,irad+1)
      xcoord(4,ispec) = x3vol(ix+4,irad+2)
      xcoord(7,ispec) = x3vol(ix+6,irad+2)
      xcoord(3,ispec) = x3vol(ix+8,irad+2)

      ycoord(1,ispec) = y3vol(ix  ,irad)
      ycoord(5,ispec) = y3vol(ix+4,irad)
      ycoord(2,ispec) = y3vol(ix+8,irad)
      ycoord(8,ispec) = y3vol(ix+2,irad+1)
      ycoord(9,ispec) = (y3vol(ix+4,irad+1) + y3vol(ix+6,irad+1))/2.
      ycoord(6,ispec) = y3vol(ix+8,irad+1)
      ycoord(4,ispec) = y3vol(ix+4,irad+2)
      ycoord(7,ispec) = y3vol(ix+6,irad+2)
      ycoord(3,ispec) = y3vol(ix+8,irad+2)

      ispec = ispec + 1

      idoubling(ispec) = IREGION_OUTER_CORE

      xcoord(1,ispec) = x3vol(ix,irad)
      xcoord(5,ispec) = x3vol(ix+2,irad+1)
      xcoord(2,ispec) = x3vol(ix+4,irad+2)
      xcoord(8,ispec) = x3vol(ix,irad+2)
      xcoord(9,ispec) = (x3vol(ix+2,irad+2) + x3vol(ix+2,irad+3))/2.
      xcoord(6,ispec) = x3vol(ix+4,irad+3)
      xcoord(4,ispec) = x3vol(ix,irad+4)
      xcoord(7,ispec) = x3vol(ix+2,irad+4)
      xcoord(3,ispec) = x3vol(ix+4,irad+4)

      ycoord(1,ispec) = y3vol(ix,irad)
      ycoord(5,ispec) = y3vol(ix+2,irad+1)
      ycoord(2,ispec) = y3vol(ix+4,irad+2)
      ycoord(8,ispec) = y3vol(ix,irad+2)
      ycoord(9,ispec) = (y3vol(ix+2,irad+2) + y3vol(ix+2,irad+3))/2.
      ycoord(6,ispec) = y3vol(ix+4,irad+3)
      ycoord(4,ispec) = y3vol(ix,irad+4)
      ycoord(7,ispec) = y3vol(ix+2,irad+4)
      ycoord(3,ispec) = y3vol(ix+4,irad+4)

  enddo

! %%%% zone Cube -> ICB %%%%

! generation maillage de la surface

  do ix=0,nspec_surf_whole_circle

  xicoord  = dble(ix)/dble(nspec_surf_whole_circle)

! parcours de l'angle dans le sens negatif with a different starting point
  thetacoord = +pi/2 - 2*pi*xicoord

! coordonnees cartesiennes correspondantes

!---- ICB

  radius_interf = RICB

  x4surf(ix) = radius_interf * dcos(thetacoord)
  y4surf(ix) = radius_interf * dsin(thetacoord)

!---- Cube

  if (thetacoord > pi/4) thetacoord = thetacoord - 2*pi

! Xmax face (right)
  if (thetacoord >= - pi/4 .and. thetacoord <= + pi/4) then
      ratio_x = +1
      ratio_y = (thetacoord - (-pi/4)) / (pi/2)

! Ymin face (bottom)
  else if (thetacoord >= - 3*pi/4 .and. thetacoord <= - pi/4) then
      ratio_x = (thetacoord - (-3*pi/4)) / (pi/2)
      ratio_y = 0

! Xmin face (left)
  else if (thetacoord >= - 5*pi/4 .and. thetacoord <= - 3*pi/4) then
      ratio_x = 0
      ratio_y = 1 - (thetacoord - (-5*pi/4)) / (pi/2)

! Ymax face (top)
  else
      ratio_x = 1 - (thetacoord - (-7*pi/4)) / (pi/2)
      ratio_y = +1

  endif

!---- volume

! use a "flat" cubed sphere to create the central cube

! map ratio to [-1,1] and then map to real radius
! then add deformation
      fact_x = 2.d0*ratio_x-1.d0
      fact_y = 2.d0*ratio_y-1.d0

      xi = PI_OVER_TWO*fact_x
      eta = PI_OVER_TWO*fact_y

      x4bot(ix) = radius_cube * fact_x * (1 + cos(eta)*CENTRAL_CUBE_INFLATE_FACTOR / PI)
      y4bot(ix) = radius_cube * fact_y * (1 + cos(xi)*CENTRAL_CUBE_INFLATE_FACTOR / PI)

!---- volume

  do irad=0,nspec_rad_Cube_ICB
      radcoord  = dble(irad)/dble(nspec_rad_Cube_ICB)
      x4vol(ix,irad) = x4surf(ix) * radcoord + x4bot(ix) * (one - radcoord)
      y4vol(ix,irad) = y4surf(ix) * radcoord + y4bot(ix) * (one - radcoord)
  enddo

  enddo

  do irad=0,nspec_rad_Cube_ICB-4,4
  do ix=0,nspec_surf_whole_circle/factor_divide_mesh-8,8

      ispec = ispec + 1

      idoubling(ispec) = IREGION_INNER_CORE

      xcoord(1,ispec) = x4vol(ix  ,irad)
      xcoord(5,ispec) = x4vol(ix+4,irad)
      xcoord(2,ispec) = x4vol(ix+8,irad)
      xcoord(8,ispec) = x4vol(ix  ,irad+2)
      xcoord(9,ispec) = x4vol(ix+4,irad+2)
      xcoord(6,ispec) = x4vol(ix+8,irad+2)
      xcoord(4,ispec) = x4vol(ix  ,irad+4)
      xcoord(7,ispec) = x4vol(ix+4,irad+4)
      xcoord(3,ispec) = x4vol(ix+8,irad+4)

      ycoord(1,ispec) = y4vol(ix  ,irad)
      ycoord(5,ispec) = y4vol(ix+4,irad)
      ycoord(2,ispec) = y4vol(ix+8,irad)
      ycoord(8,ispec) = y4vol(ix  ,irad+2)
      ycoord(9,ispec) = y4vol(ix+4,irad+2)
      ycoord(6,ispec) = y4vol(ix+8,irad+2)
      ycoord(4,ispec) = y4vol(ix  ,irad+4)
      ycoord(7,ispec) = y4vol(ix+4,irad+4)
      ycoord(3,ispec) = y4vol(ix+8,irad+4)

  enddo
  enddo

!----
!---- generer l'interieur du cube
!----

  do ix=0,nspec_surf_whole_circle/16

      xlincoord  = dble(ix)/dble(nspec_surf_whole_circle/16)

!---- volume

  do irad=0,nspec_surf_whole_circle/16
      radcoord  = dble(irad)/dble(nspec_surf_whole_circle/16)

! use a "flat" cubed sphere to create the central cube
      ratio_x = xlincoord
      ratio_y = radcoord

! map ratio to [-1,1] and then map to real radius
! then add deformation
      fact_x = 2.d0*ratio_x-1.d0
      fact_y = 2.d0*ratio_y-1.d0

      xi = PI_OVER_TWO*fact_x
      eta = PI_OVER_TWO*fact_y

      x5vol(ix,irad) = radius_cube * fact_x * (1 + cos(eta)*CENTRAL_CUBE_INFLATE_FACTOR / PI)
      y5vol(ix,irad) = radius_cube * fact_y * (1 + cos(xi)*CENTRAL_CUBE_INFLATE_FACTOR / PI)

  enddo

  enddo

  if (generate_only_half_the_mesh) then
    icentral_cube1 = nspec_surf_whole_circle/factor_divide_mesh/16
    icentral_cube2 = nspec_surf_whole_circle/16-2
  else
    icentral_cube1 = 0
    icentral_cube2 = nspec_surf_whole_circle/16-2
  endif

  do irad=0,nspec_surf_whole_circle/16-2,2
  do ix=icentral_cube1,icentral_cube2,2

      ispec = ispec + 1

      idoubling(ispec) = IREGION_INNER_CORE

      xcoord(1,ispec) = x5vol(ix  ,irad)
      xcoord(5,ispec) = x5vol(ix+1,irad)
      xcoord(2,ispec) = x5vol(ix+2,irad)
      xcoord(8,ispec) = x5vol(ix  ,irad+1)
      xcoord(9,ispec) = x5vol(ix+1,irad+1)
      xcoord(6,ispec) = x5vol(ix+2,irad+1)
      xcoord(4,ispec) = x5vol(ix  ,irad+2)
      xcoord(7,ispec) = x5vol(ix+1,irad+2)
      xcoord(3,ispec) = x5vol(ix+2,irad+2)

      ycoord(1,ispec) = y5vol(ix  ,irad)
      ycoord(5,ispec) = y5vol(ix+1,irad)
      ycoord(2,ispec) = y5vol(ix+2,irad)
      ycoord(8,ispec) = y5vol(ix  ,irad+1)
      ycoord(9,ispec) = y5vol(ix+1,irad+1)
      ycoord(6,ispec) = y5vol(ix+2,irad+1)
      ycoord(4,ispec) = y5vol(ix  ,irad+2)
      ycoord(7,ispec) = y5vol(ix+1,irad+2)
      ycoord(3,ispec) = y5vol(ix+2,irad+2)

  enddo
  enddo

! stocker le nombre total d'elements spectraux generes
  nspec = ispec

  print *
  print *,'Total number of spectral elements = ',nspec
  print *

  if (nspec /= nspec_exact) stop 'incorrect number of spectral elements generated'

!
!---- generation de la numerotation pour SPECFEM90
!

  print *,'Generating the global numbering...'

! get coordinates of the grid points
  xp(:) = 0.d0
  yp(:) = 0.d0
  do ispec=1,nspec
   ieoff = ngnod*(ispec - 1)
   ilocnum = 0

  do ia = 1,ngnod
    ilocnum = ilocnum + 1
    xp(ilocnum + ieoff) = xcoord(ia,ispec)
    yp(ilocnum + ieoff) = ycoord(ia,ispec)
  enddo

  enddo

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

!  Establish initial pointers
  do ie=1,nspec
   ieoff = ngnod*(ie -1)
   do ix=1,ngnod
      loc (ix+ieoff) = ix+ieoff
   enddo
  enddo

! set up local geometric tolerances

  xtypdist=+HUGEVAL

  do ispec=1,nspec

  xminval=+HUGEVAL
  yminval=+HUGEVAL
  xmaxval=-HUGEVAL
  ymaxval=-HUGEVAL
  ieoff=ngnod*(ispec-1)
  do ilocnum=1,ngnod
    xmaxval=max(xp(ieoff+ilocnum),xmaxval)
    xminval=min(xp(ieoff+ilocnum),xminval)
    ymaxval=max(yp(ieoff+ilocnum),ymaxval)
    yminval=min(yp(ieoff+ilocnum),yminval)
  enddo

! compute the minimum typical "size" of an element in the mesh
  xtypdist = min(xtypdist,xmaxval-xminval)
  xtypdist = min(xtypdist,ymaxval-yminval)

  enddo

! define a tolerance, small with respect to the minimum size
  xtol = SMALLVALTOL * xtypdist

  ifseg(:) = .false.
  nseg        = 1
  ifseg(1)    = .true.
  ninseg(1)   = npoin_max

  do j=1,NDIM
!  Sort within each segment
   ioff=1
   do iseg=1,nseg
      if (j == 1) then
         call rank(xp(ioff),ind,ninseg(iseg))
      else
         call rank(yp(ioff),ind,ninseg(iseg))
      endif
      call swap(xp(ioff),work,ind,ninseg(iseg))
      call swap(yp(ioff),work,ind,ninseg(iseg))
      call iswap(loc(ioff),iwork,ind,ninseg(iseg))
      ioff=ioff+ninseg(iseg)
   enddo
!  Check for jumps in current coordinate
   if (j == 1) then
     do i=2,npoin_max
       if (abs(xp(i)-xp(i-1)) > xtol) ifseg(i)=.true.
     enddo
   else
     do i=2,npoin_max
       if (abs(yp(i)-yp(i-1)) > xtol) ifseg(i)=.true.
     enddo
   endif
!  Count up number of different segments
   nseg = 0
   do i=1,npoin_max
      if (ifseg(i)) then
         nseg = nseg+1
         ninseg(nseg) = 1
      else
         ninseg(nseg) = ninseg(nseg) + 1
      endif
   enddo
  enddo
!
!  Assign global node numbers (now sorted lexicographically!)
!
  ig = 0
  do i=1,npoin_max
   if (ifseg(i)) ig=ig+1
   iglob(loc(i)) = ig
  enddo

  npoin = ig

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! verifier la coherence du nombre de points generes
  if (npoin > npoin_max) stop 'Too many points generated'

! verification de la coherence de la numerotation generee
! conversion from iglob() to ibool() is handled automatically by the "equivalence" statement
  if (minval(ibool) /= 1 .or. maxval(ibool) /= npoin) stop 'Error when generating global numbering'

  print *,'Total number of points of the global mesh: ',npoin
  print *

! generer les coordonnees des points du maillage global
! in global numbering
  print *,'Generating the coordinates of the points of the global mesh...'
  do ispec=1,nspec
  do ia = 1,ngnod
      xp(ibool(ia,ispec)) = xcoord(ia,ispec)
      yp(ibool(ia,ispec)) = ycoord(ia,ispec)
  enddo
  enddo
  print *

! ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

!
!---- write the external mesh database for SPECFEM2D
!

  print *,'Writing the external mesh database for SPECFEM2D...'

#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  print *,'using binary file format to store the mesh databases'
#else
  print *,'using ASCII file format to store the mesh databases'
#endif

! write the mesh points
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  open(unit=22,file='DATA/Nodes_AK135F_NO_MUD',form='unformatted',status='unknown',action='write')
  write(22) npoin
#else
  open(unit=22,file='DATA/Nodes_AK135F_NO_MUD',form='formatted',status='unknown',action='write')
  write(22,*) npoin
#endif
  do ipoin = 1,npoin
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
    write(22) xp(ipoin),yp(ipoin)
#else
    write(22,*) xp(ipoin),yp(ipoin)
#endif
  enddo
  close(22)

! write the mesh elements
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  open(unit=22,file='DATA/Mesh_AK135F_NO_MUD',form='unformatted',status='unknown',action='write')
  write(22) nspec
#else
  open(unit=22,file='DATA/Mesh_AK135F_NO_MUD',form='formatted',status='unknown',action='write')
  write(22,*) nspec
#endif
  do ispec = 1,nspec
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
    write(22) (ibool(ia,ispec), ia=1,ngnod)
#else
    write(22,*) (ibool(ia,ispec), ia=1,ngnod)
#endif
  enddo
  close(22)

! write the material properties
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  open(unit=22,file='DATA/Material_AK135F_NO_MUD',form='unformatted',status='unknown',action='write')
#else
  open(unit=22,file='DATA/Material_AK135F_NO_MUD',form='formatted',status='unknown',action='write')
#endif
  do ispec = 1,nspec
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
    write(22) idoubling(ispec)
#else
    write(22,*) idoubling(ispec)
#endif
  enddo
  close(22)

! write empty file for the absorbing boundary
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  open(unit=22,file='DATA/Surf_abs_AK135F_NO_MUD',form='unformatted',status='unknown',action='write')
  write(22) 0
#else
  open(unit=22,file='DATA/Surf_abs_AK135F_NO_MUD',form='formatted',status='unknown',action='write')
  write(22,*) 0
#endif
  close(22)

! write empty file for the acoustic free surface boundary (since there is no acoustic free surface in the 1D Earth without oceans)
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  open(unit=22,file='DATA/Surf_free_AK135F_NO_MUD',form='unformatted',status='unknown',action='write')
  write(22) 0
#else
  open(unit=22,file='DATA/Surf_free_AK135F_NO_MUD',form='formatted',status='unknown',action='write')
  write(22,*) 0
#endif
  close(22)

  print *,'Done writing the external mesh database for SPECFEM2D'
  print *

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! list all the elements that are in contact with the symmetry axis, and by what edge

  ispec_count = 0

! count the number of elements that are in contact with the symmetry axis
  do ispec=1,nspec
    i = 0

    if (xcoord(1,ispec) < 0.001d0) i = i + 1
    if (xcoord(2,ispec) < 0.001d0) i = i + 1
    if (xcoord(3,ispec) < 0.001d0) i = i + 1
    if (xcoord(4,ispec) < 0.001d0) i = i + 1

    if (i > 0) then
! we know in advance from the way the mesh is designed that single points can never be detected, only a full edge,
! and a single edge can of course be detected otherwise that element would have a 180-degree angle
      if (i /= 2) stop 'error: element in contact with the symmetry axis by incorrect number of points'
      ispec_count = ispec_count + 1
    endif
  enddo

  print *,'number of elements in contact with the symmetry axis = ',ispec_count

! then save them to a file and also save the (only) edge that is in contact
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
  open(unit=22,file='DATA/Symmetry_axis_elements_AK135F_NO_MUD',form='unformatted',status='unknown',action='write')
  write(22) ispec_count
#else
  open(unit=22,file='DATA/Symmetry_axis_elements_AK135F_NO_MUD',form='formatted',status='unknown',action='write')
  write(22,*) ispec_count
#endif

  do ispec=1,nspec
#ifdef USE_BINARY_FOR_EXTERNAL_MESH_DATABASE
    if (xcoord(1,ispec) < 0.001d0 .and. xcoord(2,ispec) < 0.001d0) write(22) ispec,' 2 ',ibool(1,ispec),ibool(2,ispec),IBOTTOM
    if (xcoord(2,ispec) < 0.001d0 .and. xcoord(3,ispec) < 0.001d0) write(22) ispec,' 2 ',ibool(2,ispec),ibool(3,ispec),IRIGHT
    if (xcoord(3,ispec) < 0.001d0 .and. xcoord(4,ispec) < 0.001d0) write(22) ispec,' 2 ',ibool(3,ispec),ibool(4,ispec),ITOP
    if (xcoord(4,ispec) < 0.001d0 .and. xcoord(1,ispec) < 0.001d0) write(22) ispec,' 2 ',ibool(4,ispec),ibool(1,ispec),ILEFT
#else
    if (xcoord(1,ispec) < 0.001d0 .and. xcoord(2,ispec) < 0.001d0) write(22,*) ispec,' 2 ',ibool(1,ispec),ibool(2,ispec),IBOTTOM
    if (xcoord(2,ispec) < 0.001d0 .and. xcoord(3,ispec) < 0.001d0) write(22,*) ispec,' 2 ',ibool(2,ispec),ibool(3,ispec),IRIGHT
    if (xcoord(3,ispec) < 0.001d0 .and. xcoord(4,ispec) < 0.001d0) write(22,*) ispec,' 2 ',ibool(3,ispec),ibool(4,ispec),ITOP
    if (xcoord(4,ispec) < 0.001d0 .and. xcoord(1,ispec) < 0.001d0) write(22,*) ispec,' 2 ',ibool(4,ispec),ibool(1,ispec),ILEFT
#endif
  enddo

  close(22)

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! ***
! *** generer un fichier 'GNUPLOT' pour le controle de la grille ***
! ***

  if (output_gnuplot_grid) then

    print *
    print *,'Writing the grid in GNUPLOT format...'

    open(unit=20,file='gridfile.gnu',status='unknown')

    do ispec=1,nspec

  ! draw the four edges of each element (using straight lines to simplify)
      ia1 = 1
      ia2 = 2
      write(20,15) sngl(xcoord(ia1,ispec)),sngl(ycoord(ia1,ispec))
      write(20,15) sngl(xcoord(ia2,ispec)),sngl(ycoord(ia2,ispec))
      write(20,10)

      ia1 = 2
      ia2 = 3
      write(20,15) sngl(xcoord(ia1,ispec)),sngl(ycoord(ia1,ispec))
      write(20,15) sngl(xcoord(ia2,ispec)),sngl(ycoord(ia2,ispec))
      write(20,10)

      ia1 = 3
      ia2 = 4
      write(20,15) sngl(xcoord(ia1,ispec)),sngl(ycoord(ia1,ispec))
      write(20,15) sngl(xcoord(ia2,ispec)),sngl(ycoord(ia2,ispec))
      write(20,10)

      ia1 = 4
      ia2 = 1
      write(20,15) sngl(xcoord(ia1,ispec)),sngl(ycoord(ia1,ispec))
      write(20,15) sngl(xcoord(ia2,ispec)),sngl(ycoord(ia2,ispec))
      write(20,10)

    enddo

    close(20)

  ! cree le script de dessin pour gnuplot
    open(unit=20,file='plotgrid.gnu',status='unknown')
    write(20,*) '#set term postscript landscape monochrome solid "Helvetica" 22'
    write(20,*) '#set output "grille.ps"'
    write(20,*) 'set term x11'
    write(20,*) 'set size ratio -1'
    write(20,*) 'plot "gridfile.gnu" title "Macroblocs mesh" w l'
    write(20,*) 'pause -1 "Hit any key..."'
    close(20)

    print *,'Done writing the grid in GNUPLOT format'
    print *

 10   format('')
 15   format(e12.5,1x,e12.5)

  endif

  print *
  print *,'All Done'
  print *

  end

!-----------------------------------------------------------------------

  subroutine rank(A,IND,N)
!
! Use Heap Sort (p 233 Numerical Recipes)
!
  implicit none

  integer N
  double precision A(N)
  integer IND(N)

  integer i,j,l,ir,indx
  double precision q

  do J=1,N
   IND(j)=j
  enddo

  if (n == 1) return
  L=n/2+1
  ir=n
  100 continue
   if (l > 1) then
     l=l-1
     indx=ind(l)
     q=a(indx)
   ELSE
     indx=ind(ir)
     q=a(indx)
     ind(ir)=ind(1)
     ir=ir-1
     if (ir == 1) then
       ind(1)=indx
       return
     endif
   endif
   i=l
   j=l+l
  200 continue
   if (J <= IR) then
      if (J < IR) then
         if (A(IND(j)) < A(IND(j+1))) j=j+1
      endif
      if (q < A(IND(j))) then
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      endif
   goto 200
   endif
   IND(I)=INDX
  goto 100

  end subroutine rank

!-----------------------------------------------------------------------

  subroutine swap(a,w,ind,n)
!
! Use IND to sort array A (p 233 Numerical Recipes)
!
  implicit none

  integer n
  double precision A(N),W(N)
  integer IND(N)

  integer j

  W(:) = A(:)

  do J=1,N
    A(j) = W(ind(j))
  enddo

  end subroutine swap

!-----------------------------------------------------------------------

  subroutine iswap(a,w,ind,n)
!
! Use IND to sort array A
!
  implicit none

  integer n
  integer A(N),W(N),IND(N)

  integer j

  W(:) = A(:)

  do J=1,N
    A(j) = W(ind(j))
  enddo

  end subroutine iswap


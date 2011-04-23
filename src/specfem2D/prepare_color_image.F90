
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


  subroutine prepare_color_image_init(NX_IMAGE_color,NZ_IMAGE_color, &
                            xmin_color_image,xmax_color_image, &
                            zmin_color_image,zmax_color_image, &
                            coord,nglob,npgeo)

  implicit none
  include "constants.h"
#ifdef USE_MPI
  include "mpif.h"
#endif

  integer :: NX_IMAGE_color,NZ_IMAGE_color

  integer :: nglob,npgeo
  double precision, dimension(NDIM,nglob) :: coord

  double precision :: xmin_color_image,xmax_color_image, &
    zmin_color_image,zmax_color_image

  ! local parameters
  integer  :: npgeo_glob
  double precision  :: xmin_color_image_loc, xmax_color_image_loc, &
      zmin_color_image_loc,zmax_color_image_loc
#ifdef USE_MPI
  integer :: ier
#endif

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
  NZ_IMAGE_color = nint(NX_IMAGE_color * (zmax_color_image - zmin_color_image) &
                                      / (xmax_color_image - xmin_color_image))

  ! convert pixel sizes to even numbers because easier to reduce size,
  ! create MPEG movies in postprocessing
  NX_IMAGE_color = 2 * (NX_IMAGE_color / 2)
  NZ_IMAGE_color = 2 * (NZ_IMAGE_color / 2)

  ! check that image size is not too big
  if (NX_IMAGE_color > 99999) call exit_MPI('output image too big : NX_IMAGE_color > 99999.')
  if (NZ_IMAGE_color > 99999) call exit_MPI('output image too big : NZ_IMAGE_color > 99999.')

  end subroutine prepare_color_image_init


!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_color_image_pixels(myrank,NX_IMAGE_color,NZ_IMAGE_color, &
                            xmin_color_image,xmax_color_image, &
                            zmin_color_image,zmax_color_image, &
                            coord,nglob,coorg,npgeo,nspec,ngnod,knods,ibool, &
                            nb_pixel_loc,iglob_image_color)

  implicit none
  include "constants.h"

  integer :: myrank
  integer :: NX_IMAGE_color,NZ_IMAGE_color
  double precision :: xmin_color_image,xmax_color_image, &
    zmin_color_image,zmax_color_image

  integer :: nglob,nspec,npgeo,ngnod
  double precision, dimension(NDIM,nglob) :: coord
  double precision, dimension(NDIM,npgeo) :: coorg

  integer, dimension(ngnod,nspec) :: knods
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool

  integer :: nb_pixel_loc
  integer, dimension(NX_IMAGE_color,NZ_IMAGE_color) :: iglob_image_color

  ! local parameters
  double precision  :: size_pixel_horizontal,size_pixel_vertical
  double precision, dimension(2,4)  :: elmnt_coords
  double precision  :: i_coord, j_coord
  double precision  :: dist_pixel, dist_min_pixel
  integer  :: min_i, min_j, max_i, max_j
  integer  :: ispec,i,j,k,l,iglob
  logical  :: pixel_is_in

  ! create all the pixels
  if (myrank == 0) then
    write(IOUT,*)
    write(IOUT,*) 'locating all the pixels of color images'
  endif

  size_pixel_horizontal = (xmax_color_image - xmin_color_image) / dble(NX_IMAGE_color-1)
  size_pixel_vertical = (zmax_color_image - zmin_color_image) / dble(NZ_IMAGE_color-1)

  iglob_image_color(:,:) = -1

  ! checking which pixels are inside each elements

  nb_pixel_loc = 0
  do ispec = 1, nspec

    do k = 1, 4
      elmnt_coords(1,k) = coorg(1,knods(k,ispec))
      elmnt_coords(2,k) = coorg(2,knods(k,ispec))
    enddo

    ! avoid working on the whole pixel grid
    min_i = floor(minval((elmnt_coords(1,:) - xmin_color_image))/size_pixel_horizontal) + 1
    max_i = ceiling(maxval((elmnt_coords(1,:) - xmin_color_image))/size_pixel_horizontal) + 1
    min_j = floor(minval((elmnt_coords(2,:) - zmin_color_image))/size_pixel_vertical) + 1
    max_j = ceiling(maxval((elmnt_coords(2,:) - zmin_color_image))/size_pixel_vertical) + 1

    ! avoid edge effects
    if(min_i < 1) min_i = 1
    if(min_j < 1) min_j = 1

    if(max_i > NX_IMAGE_color) max_i = NX_IMAGE_color
    if(max_j > NZ_IMAGE_color) max_j = NZ_IMAGE_color

     do j = min_j, max_j
        do i = min_i, max_i
           i_coord = (i-1)*size_pixel_horizontal + xmin_color_image
           j_coord = (j-1)*size_pixel_vertical + zmin_color_image

           ! checking if the pixel is inside the element (must be a convex quadrilateral)
           call is_in_convex_quadrilateral( elmnt_coords, i_coord, j_coord, pixel_is_in)

           ! if inside, getting the nearest point inside the element!
           if ( pixel_is_in ) then
              dist_min_pixel = HUGEVAL
              do k = 1, NGLLX
                 do l = 1, NGLLZ
                    iglob = ibool(k,l,ispec)
                    dist_pixel = (coord(1,iglob)-i_coord)**2 + (coord(2,iglob)-j_coord)**2
                    if (dist_pixel < dist_min_pixel) then
                       dist_min_pixel = dist_pixel
                       iglob_image_color(i,j) = iglob

                    endif

                 enddo
              enddo
              if ( dist_min_pixel >= HUGEVAL ) then
                 call exit_MPI('Error in detecting pixel for color image')

              endif
              nb_pixel_loc = nb_pixel_loc + 1

           endif

        enddo
     enddo
  enddo

  end subroutine prepare_color_image_pixels


!
!-------------------------------------------------------------------------------------------------
!


  subroutine prepare_color_image_vp(nglob,image_color_vp_display,iglob_image_color, &
                            NX_IMAGE_color,NZ_IMAGE_color,nb_pixel_loc, &
                            num_pixel_loc,nspec,poroelastic,ibool,kmato, &
                            numat,density,poroelastcoef,porosity,tortuosity, &
                            nproc,myrank,assign_external_model,vpext)

! stores P-velocity model in image_color_vp_display

  implicit none
  include "constants.h"
#ifdef USE_MPI
  include "mpif.h"
#endif

  integer :: nglob,nspec
  integer :: NX_IMAGE_color,NZ_IMAGE_color
  double precision, dimension(NX_IMAGE_color,NZ_IMAGE_color) :: image_color_vp_display
  integer, dimension(NX_IMAGE_color,NZ_IMAGE_color) :: iglob_image_color

  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato

  logical, dimension(nspec) :: poroelastic

  integer :: nb_pixel_loc
  integer, dimension(nb_pixel_loc) :: num_pixel_loc

  logical :: assign_external_model
  integer :: nproc,myrank
  integer :: numat
  double precision, dimension(2,numat) :: density
  double precision, dimension(4,3,numat) :: poroelastcoef
  double precision, dimension(numat) :: porosity,tortuosity
  double precision, dimension(NGLLX,NGLLX,nspec) :: vpext

  ! local parameters
  double precision, dimension(:), allocatable :: vp_display
  double precision :: rhol,mul_relaxed,lambdal_relaxed
  double precision :: rhol_s,rhol_f,rhol_bar,phil,tortl,mul_s,kappal_s,kappal_f, &
    mul_fr,kappal_fr
  double precision :: afactor,bfactor,cfactor,D_biot,H_biot,C_biot,&
    M_biot,B_biot,cpIsquare,cpIIsquare,cssquare
  double precision :: gamma1,gamma2,gamma3,gamma4,ratio
  integer  :: i,j,k,ispec
#ifdef USE_MPI
  double precision, dimension(:), allocatable  :: data_pixel_recv
  double precision, dimension(:), allocatable  :: data_pixel_send
  integer, dimension(:,:), allocatable  :: num_pixel_recv
  integer, dimension(:), allocatable  :: nb_pixel_per_proc
  integer, dimension(MPI_STATUS_SIZE)  :: request_mpi_status
  integer :: ier,iproc
#else
  integer :: dummy
#endif

  ! to display the P-velocity model in background on color images
  allocate(vp_display(nglob))

  do ispec = 1,nspec

    if(poroelastic(ispec)) then
      !get parameters of current spectral element
      phil = porosity(kmato(ispec))
      tortl = tortuosity(kmato(ispec))
      !solid properties
      mul_s = poroelastcoef(2,1,kmato(ispec))
      kappal_s = poroelastcoef(3,1,kmato(ispec)) - 4.d0*mul_s/3.d0
      rhol_s = density(1,kmato(ispec))
      !fluid properties
      kappal_f = poroelastcoef(1,2,kmato(ispec))
      rhol_f = density(2,kmato(ispec))
      !frame properties
      mul_fr = poroelastcoef(2,3,kmato(ispec))
      kappal_fr = poroelastcoef(3,3,kmato(ispec)) - 4.d0*mul_fr/3.d0
      rhol_bar =  (1.d0 - phil)*rhol_s + phil*rhol_f
      !Biot coefficients for the input phi
      D_biot = kappal_s*(1.d0 + phil*(kappal_s/kappal_f - 1.d0))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) &
              + kappal_fr + 4.d0*mul_fr/3.d0
      C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
      M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
      B_biot = H_biot - 4.d0*mul_fr/3.d0
      ! Approximated velocities (no viscous dissipation)
      afactor = rhol_bar - phil/tortl*rhol_f
      bfactor = H_biot + phil*rhol_bar/(tortl*rhol_f)*M_biot - TWO*phil/tortl*C_biot
      cfactor = phil/(tortl*rhol_f)*(H_biot*M_biot - C_biot*C_biot)
      cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
      cpIIsquare = (bfactor - sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
      cssquare = mul_fr/afactor

      ! Approximated ratio r = amplitude "w" field/amplitude "s" field (no viscous dissipation)
      ! used later for wavespeed kernels calculation, which are presently implemented for inviscid case,
      ! contrary to primary and density-normalized kernels, which are consistent with viscous fluid case.
      gamma1 = H_biot - phil/tortl*C_biot
      gamma2 = C_biot - phil/tortl*M_biot
      gamma3 = phil/tortl*( M_biot*(afactor/rhol_f + phil/tortl) - C_biot)
      gamma4 = phil/tortl*( C_biot*(afactor/rhol_f + phil/tortl) - H_biot)
      ratio = HALF*(gamma1 - gamma3)/gamma4 &
            + HALF*sqrt((gamma1-gamma3)**2/gamma4**2 &
            + 4.d0 * gamma2/gamma4)

      do j = 1,NGLLZ
        do i = 1,NGLLX
          vp_display(ibool(i,j,ispec)) = sqrt(cpIsquare)
        enddo
      enddo

    else
      ! get relaxed elastic parameters of current spectral element
      rhol = density(1,kmato(ispec))
      lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
      mul_relaxed = poroelastcoef(2,1,kmato(ispec))
      do j = 1,NGLLZ
        do i = 1,NGLLX
          !--- if external medium, get elastic parameters of current grid point
          if(assign_external_model) then
            vp_display(ibool(i,j,ispec)) = vpext(i,j,ispec)
          else
            vp_display(ibool(i,j,ispec)) = sqrt((lambdal_relaxed + 2.d0*mul_relaxed) / rhol)
          endif
        enddo
      enddo
    endif !if(poroelastic(ispec)) then
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

  allocate(nb_pixel_per_proc(nproc))
  nb_pixel_per_proc(:) = 0
  call MPI_GATHER( nb_pixel_loc, 1, MPI_INTEGER, nb_pixel_per_proc(1), &
                  1, MPI_INTEGER, 0, MPI_COMM_WORLD, ier)


  if ( myrank == 0 ) then
     allocate(num_pixel_recv(maxval(nb_pixel_per_proc(:)),nproc))
     allocate(data_pixel_recv(maxval(nb_pixel_per_proc(:))))
  endif
  allocate(data_pixel_send(nb_pixel_loc))

  if (nproc > 1) then
    if (myrank == 0) then
      do iproc = 1, nproc-1

        call MPI_RECV(num_pixel_recv(1,iproc+1),nb_pixel_per_proc(iproc+1), MPI_INTEGER, &
                iproc, 42, MPI_COMM_WORLD, request_mpi_status, ier)

        call MPI_RECV(data_pixel_recv(1),nb_pixel_per_proc(iproc+1), MPI_DOUBLE_PRECISION, &
                iproc, 43, MPI_COMM_WORLD, request_mpi_status, ier)

        do k = 1, nb_pixel_per_proc(iproc+1)
          j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
          i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

          ! checks bounds
          if( i < 1 .or. i > NX_IMAGE_color .or. j < 1 .or. j > NZ_IMAGE_color ) then
            print*,'image vp bounds:',myrank,iproc,k, &
              num_pixel_recv(k,iproc+1),nb_pixel_per_proc(iproc+1)
            print*,'  i: ',i,NX_IMAGE_color
            print*,'  j: ',j,NZ_IMAGE_color
          endif

          image_color_vp_display(i,j) = data_pixel_recv(k)
        enddo
      enddo

    else
      do k = 1, nb_pixel_loc
        j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
        i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color
        data_pixel_send(k) = vp_display(iglob_image_color(i,j))
      enddo

      call MPI_SEND(num_pixel_loc(1),nb_pixel_loc,MPI_INTEGER, &
              0, 42, MPI_COMM_WORLD, ier)

      call MPI_SEND(data_pixel_send(1),nb_pixel_loc,MPI_DOUBLE_PRECISION, &
              0, 43, MPI_COMM_WORLD, ier)

    endif
  endif

  deallocate(nb_pixel_per_proc)
  deallocate(data_pixel_send)
  if( myrank == 0 ) then
    deallocate(num_pixel_recv)
    deallocate(data_pixel_recv)
  endif
#else
  ! to avoid compiler warnings
  dummy = myrank
  dummy = nproc
#endif

  end subroutine prepare_color_image_vp

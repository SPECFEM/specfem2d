!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

  subroutine prepare_color_image_init()

  use constants, only: NX_NZ_IMAGE_MAX,NGLLX

  use specfem_par, only: coord,npgeo,myrank

  use specfem_par_movie, only: NX_IMAGE_color,NZ_IMAGE_color, &
                          xmin_color_image,xmax_color_image, &
                          zmin_color_image,zmax_color_image,factor_subsample_image
  implicit none

  ! local parameters
  integer  :: npgeo_glob
  double precision  :: xmin_color_image_loc, xmax_color_image_loc, &
                       zmin_color_image_loc,zmax_color_image_loc

  ! horizontal min and max coordinates of the image
  xmin_color_image_loc = minval(coord(1,:))
  xmax_color_image_loc = maxval(coord(1,:))

  ! vertical min and max coordinates of the image, slightly increase it to go beyond maximum topography
  zmin_color_image_loc = minval(coord(2,:))
  zmax_color_image_loc = maxval(coord(2,:))

! global values
  xmin_color_image = xmin_color_image_loc
  xmax_color_image = xmax_color_image_loc
  zmin_color_image = zmin_color_image_loc
  zmax_color_image = zmax_color_image_loc
  npgeo_glob = npgeo

  call min_all_all_dp(xmin_color_image_loc, xmin_color_image)
  call max_all_all_dp(xmax_color_image_loc, xmax_color_image)
  call min_all_all_dp(zmin_color_image_loc, zmin_color_image)
  call max_all_all_dp(zmax_color_image_loc, zmax_color_image)

  ! collects total on all processes
  call sum_all_all_i(npgeo, npgeo_glob)

  zmax_color_image = zmin_color_image + 1.05d0 * (zmax_color_image - zmin_color_image)

  ! compute number of pixels in the horizontal direction based on typical number
  ! of spectral elements in a given direction (may give bad results for very elongated models)
  NX_IMAGE_color = nint(sqrt(dble(npgeo_glob))) * (NGLLX-1) + 1

  ! compute number of pixels in the vertical direction based on ratio of sizes
  NZ_IMAGE_color = nint(NX_IMAGE_color * (zmax_color_image - zmin_color_image) &
                                      / (xmax_color_image - xmin_color_image))

  ! convert pixel sizes to even numbers because easier to reduce size,
  ! create MPEG movies in postprocessing
  NX_IMAGE_color = 2 * nint((NX_IMAGE_color / 2 + 1) / factor_subsample_image)
  NZ_IMAGE_color = 2 * nint((NZ_IMAGE_color / 2 + 1) / factor_subsample_image)

  ! check that image size is not too big
! because from http://www.jpegcameras.com/libjpeg/libjpeg-2.html
! we know that the size limit of the image in each dimension is 65535:
! "JPEG supports image dimensions of 1 to 64K pixels in either direction".
  if (NX_IMAGE_color < 4) call exit_MPI(myrank,'output image too small: NX_IMAGE_color < 4.')
  if (NZ_IMAGE_color < 4) call exit_MPI(myrank,'output image too small: NZ_IMAGE_color < 4.')

  if (NX_IMAGE_color > 65534) &
    call exit_MPI(myrank,'output image too big: NX_IMAGE_color > 65534; increase factor_subsample_image in DATA/Par_file.')
  if (NZ_IMAGE_color > 65534) &
    call exit_MPI(myrank,'output image too big: NZ_IMAGE_color > 65534; increase factor_subsample_image in DATA/Par_file.')

  if (NX_IMAGE_color > NX_NZ_IMAGE_MAX) then
    print *,'NX_IMAGE_color,NX_NZ_IMAGE_MAX = ',NX_IMAGE_color,NX_NZ_IMAGE_MAX
    call exit_MPI(myrank, &
      'output image too big: NX_IMAGE_color > NX_NZ_IMAGE_MAX; increase factor_subsample_image or change NX_NZ_IMAGE_MAX.')
  endif
  if (NZ_IMAGE_color > NX_NZ_IMAGE_MAX) then
    print *,'NZ_IMAGE_color,NX_NZ_IMAGE_MAX = ',NZ_IMAGE_color,NX_NZ_IMAGE_MAX
    call exit_MPI(myrank, &
      'output image too big: NZ_IMAGE_color > NX_NZ_IMAGE_MAX; increase factor_subsample_image or change NX_NZ_IMAGE_MAX.')
  endif

  end subroutine prepare_color_image_init


!
!-------------------------------------------------------------------------------------------------
!

  subroutine prepare_color_image_pixels()

  use constants, only: IMAIN,HUGEVAL,NGLLX,NGLLZ

  use specfem_par, only: myrank,coord,coorg,nspec,knods,ibool, &
                          NSOURCES,nrec,x_source,z_source,st_xval,st_zval


  use specfem_par_movie, only: NX_IMAGE_color,NZ_IMAGE_color, &
                              xmin_color_image,xmax_color_image, &
                              zmin_color_image,zmax_color_image, &
                              nb_pixel_loc,iglob_image_color, &
                              DRAW_SOURCES_AND_RECEIVERS, &
                              ix_image_color_source,iy_image_color_source,ix_image_color_receiver,iy_image_color_receiver

  implicit none

  ! local parameters
  double precision  :: size_pixel_horizontal,size_pixel_vertical
  double precision, dimension(2,4)  :: elmnt_coords
  double precision  :: i_coord, j_coord
  double precision  :: dist_pixel, dist_min_pixel
  integer  :: min_i, min_j, max_i, max_j
  integer  :: ispec,i,j,k,l,iglob,ier,pixel_k,pixel_l
  integer  :: nb_pixel_total
  logical  :: pixel_is_in

  ! create all the pixels
  if (myrank == 0) then
    write(IMAIN,*) '  locating all the pixels of color images'
    call flush_IMAIN()
  endif

  size_pixel_horizontal = (xmax_color_image - xmin_color_image) / dble(NX_IMAGE_color-1)
  size_pixel_vertical = (zmax_color_image - zmin_color_image) / dble(NZ_IMAGE_color-1)

  ! initializes
  iglob_image_color(:,:) = -1

  ! checking which pixels are inside each element
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
    if (min_i < 1) min_i = 1
    if (min_j < 1) min_j = 1

    if (max_i > NX_IMAGE_color) max_i = NX_IMAGE_color
    if (max_j > NZ_IMAGE_color) max_j = NZ_IMAGE_color

    do j = min_j, max_j
      do i = min_i, max_i
        i_coord = (i-1)*size_pixel_horizontal + xmin_color_image
        j_coord = (j-1)*size_pixel_vertical + zmin_color_image

        ! checking if the pixel is inside the element (must be a convex quadrilateral)
        call is_in_convex_quadrilateral( elmnt_coords, i_coord, j_coord, pixel_is_in)

        ! if inside, getting the nearest point inside the element!
        if (pixel_is_in) then
          dist_min_pixel = HUGEVAL
          do k = 1, NGLLX
            do l = 1, NGLLZ
              iglob = ibool(k,l,ispec)
              dist_pixel = (coord(1,iglob)-i_coord)**2 + (coord(2,iglob)-j_coord)**2
              if (dist_pixel < dist_min_pixel) then
                dist_min_pixel = dist_pixel
                pixel_l = l
                pixel_k = k
              endif
            enddo
          enddo
          ! checks if pixel found
          if (dist_min_pixel >= HUGEVAL) &
            call exit_MPI(myrank,'Error in detecting pixel for color image')

          ! sets closest GLL point for pixel location
          if (iglob_image_color(i,j) == -1) then
            ! sets new pixel
            iglob = ibool(pixel_k,pixel_l,ispec)
            iglob_image_color(i,j) = iglob
            nb_pixel_loc = nb_pixel_loc + 1
          endif
        endif
      enddo
    enddo
  enddo

  ! user output
  call sum_all_i(nb_pixel_loc,nb_pixel_total)
  if (myrank == 0) then
    write(IMAIN,*) '  total number of image pixels = ',nb_pixel_total
    call flush_IMAIN()
  endif


!
!----  find pixel position of the sources and receivers
!
  if (DRAW_SOURCES_AND_RECEIVERS .and. myrank == 0) then

    ! find pixel position of the sources with orange crosses
    allocate(ix_image_color_source(NSOURCES), &
             iy_image_color_source(NSOURCES),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating image source arrays')

    do i = 1,NSOURCES
      ix_image_color_source(i) = int((x_source(i) - xmin_color_image) / size_pixel_horizontal) + 1
      iy_image_color_source(i) = int((z_source(i) - zmin_color_image) / size_pixel_vertical) + 1

      ! avoid edge effects
      if (ix_image_color_source(i) < 1) ix_image_color_source(i) = 1
      if (iy_image_color_source(i) < 1) iy_image_color_source(i) = 1

      if (ix_image_color_source(i) > NX_IMAGE_color) ix_image_color_source(i) = NX_IMAGE_color
      if (iy_image_color_source(i) > NZ_IMAGE_color) iy_image_color_source(i) = NZ_IMAGE_color
    enddo

    ! find pixel position of the receivers with green squares
    allocate(ix_image_color_receiver(nrec), &
             iy_image_color_receiver(nrec),stat=ier)
    if (ier /= 0) call stop_the_code('Error allocating image receiver arrays')

    do i = 1,nrec
      ix_image_color_receiver(i) = int((st_xval(i) - xmin_color_image) / size_pixel_horizontal) + 1
      iy_image_color_receiver(i) = int((st_zval(i) - zmin_color_image) / size_pixel_vertical) + 1

      ! avoid edge effects
      if (ix_image_color_receiver(i) < 1) ix_image_color_receiver(i) = 1
      if (iy_image_color_receiver(i) < 1) iy_image_color_receiver(i) = 1

      if (ix_image_color_receiver(i) > NX_IMAGE_color) ix_image_color_receiver(i) = NX_IMAGE_color
      if (iy_image_color_receiver(i) > NZ_IMAGE_color) iy_image_color_receiver(i) = NZ_IMAGE_color
    enddo

  endif

  end subroutine prepare_color_image_pixels


!
!-------------------------------------------------------------------------------------------------
!


  subroutine prepare_color_image_vp()

! stores P-velocity model in image_color_vp_display

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: NGLLX,NGLLZ,HALF,TWO,IMAIN

  use specfem_par, only: nglob,nspec,ispec_is_elastic,ispec_is_poroelastic,ibool,kmato, &
    density,poroelastcoef,NPROC,myrank,assign_external_model,vpext

  use specfem_par_movie, only: image_color_vp_display,iglob_image_color, &
    NX_IMAGE_color,NZ_IMAGE_color,nb_pixel_loc,num_pixel_loc,DRAW_WATER_IN_BLUE

  implicit none

  ! local parameters
  double precision :: vp_of_the_model
  double precision, dimension(:), allocatable :: vp_display
  double precision :: rhol,mul_relaxed,lambdal_relaxed

  double precision :: phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot

  double precision :: afactor,bfactor,cfactor
  double precision :: cpIsquare

  integer  :: i,j,k,ispec
#ifdef USE_MPI
  double precision, dimension(:), allocatable  :: data_pixel_recv
  double precision, dimension(:), allocatable  :: data_pixel_send
  integer, dimension(:,:), allocatable  :: num_pixel_recv
  integer, dimension(:), allocatable  :: tmp_nb_pixel_per_proc
  integer :: iproc
#else
  integer :: dummy
#endif

  ! user output
  if (myrank == 0) then
    write(IMAIN,*) '  coloring image background based on vp'
    call flush_IMAIN()
  endif

  ! to display the P-velocity model in background on color images
  allocate(vp_display(nglob))

  do ispec = 1,nspec

    if (ispec_is_poroelastic(ispec)) then
      !get parameters of current spectral element
      call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

      ! Biot coefficients for the input phi
      call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

      ! Approximated velocities (no viscous dissipation)
      afactor = rho_bar - phi/tort*rho_f
      bfactor = H_biot + phi*rho_bar/(tort*rho_f)*M_biot - TWO*phi/tort*C_biot
      cfactor = phi/(tort*rho_f)*(H_biot*M_biot - C_biot*C_biot)

      cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)

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
          if (assign_external_model) then
            vp_display(ibool(i,j,ispec)) = vpext(i,j,ispec)
          else
            vp_display(ibool(i,j,ispec)) = sqrt((lambdal_relaxed + 2.d0*mul_relaxed) / rhol)
          endif
        enddo
      enddo
    endif ! of if (ispec_is_poroelastic(ispec)) then

! now display acoustic layers as constant blue, because they likely correspond to water in the case of ocean acoustics
! or in the case of offshore oil industry experiments.
! (simply turn DRAW_WATER_IN_BLUE off in DATA/Par_file if you want to suppress this (for instance when running
!  a purely acoustic simulation with different acoustic media for the oil industry, one then wants to see the different
!  acoustic wave speeds displayed as a grey scale).
! For now, in this routine, use -1 as a flag to label such acoustic points
    if (DRAW_WATER_IN_BLUE .and. .not. ispec_is_elastic(ispec) .and. .not. ispec_is_poroelastic(ispec)) then
      do j = 1,NGLLZ
        do i = 1,NGLLX

          !--- if external medium, get elastic parameters of current grid point
          if (assign_external_model) then
            vp_of_the_model = vpext(i,j,ispec)
          else
            vp_of_the_model = sqrt((lambdal_relaxed + 2.d0*mul_relaxed) / rhol)
          endif

! test that water is indeed water and not an acoustic version of a sediment for instance
! thus check that Vp is the typical Vp of water
          if (abs(vp_of_the_model - 1480.d0) <= 50.d0) vp_display(ibool(i,j,ispec)) = -1

        enddo
      enddo
    endif

  enddo

  ! getting velocity for each local pixels
  image_color_vp_display(:,:) = 0.d0

  do k = 1, nb_pixel_loc
    j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
    i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

! avoid edge effects
    if (i < 1) i = 1
    if (i > NX_IMAGE_color) i = NX_IMAGE_color

    if (j < 1) j = 1
    if (j > NZ_IMAGE_color) j = NZ_IMAGE_color

    if (iglob_image_color(i,j) /= -1) image_color_vp_display(i,j) = vp_display(iglob_image_color(i,j))
  enddo

! assembling array image_color_vp_display on process zero for color output
#ifdef USE_MPI

  allocate(tmp_nb_pixel_per_proc(0:NPROC-1))
  tmp_nb_pixel_per_proc(:) = 0
  call gather_all_singlei(nb_pixel_loc,tmp_nb_pixel_per_proc,NPROC)

  if (myrank == 0) then
     allocate(num_pixel_recv(maxval(tmp_nb_pixel_per_proc(:)),NPROC))
     allocate(data_pixel_recv(maxval(tmp_nb_pixel_per_proc(:))))
  endif
  allocate(data_pixel_send(nb_pixel_loc))

  if (NPROC > 1) then
    if (myrank == 0) then
      do iproc = 1, NPROC-1

        call recv_i(num_pixel_recv(1,iproc+1), tmp_nb_pixel_per_proc(iproc), iproc, 42)
        call recv_dp(data_pixel_recv(1),tmp_nb_pixel_per_proc(iproc), iproc, 43)

        do k = 1, tmp_nb_pixel_per_proc(iproc)
          j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
          i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

! avoid edge effects
          if (i < 1) i = 1
          if (i > NX_IMAGE_color) i = NX_IMAGE_color

          if (j < 1) j = 1
          if (j > NZ_IMAGE_color) j = NZ_IMAGE_color

          image_color_vp_display(i,j) = data_pixel_recv(k)
        enddo
      enddo

    else
      do k = 1, nb_pixel_loc
        j = ceiling(real(num_pixel_loc(k)) / real(NX_IMAGE_color))
        i = num_pixel_loc(k) - (j-1)*NX_IMAGE_color

! avoid edge effects
        if (i < 1) i = 1
        if (i > NX_IMAGE_color) i = NX_IMAGE_color

        if (j < 1) j = 1
        if (j > NZ_IMAGE_color) j = NZ_IMAGE_color

        if (iglob_image_color(i,j) /= -1) data_pixel_send(k) = vp_display(iglob_image_color(i,j))
      enddo

      call send_i(num_pixel_loc(1), nb_pixel_loc, 0, 42)
      call send_dp(data_pixel_send(1), nb_pixel_loc, 0, 43)

    endif
  endif

  deallocate(tmp_nb_pixel_per_proc)
  deallocate(data_pixel_send)
  if (myrank == 0) then
    deallocate(num_pixel_recv)
    deallocate(data_pixel_recv)
  endif
#else
  ! to avoid compiler warnings
  dummy = myrank
  dummy = NPROC
#endif

  end subroutine prepare_color_image_vp


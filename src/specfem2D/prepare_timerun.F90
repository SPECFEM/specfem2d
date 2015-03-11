subroutine prepare_timerun_mass_matrix()


#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

integer ispec
! inner/outer elements in the case of an MPI simulation
  integer :: ispec_inner,ispec_outer

#ifdef USE_MPI
  include "precision.h"
#endif


  !
  !---- build the global mass matrix
  !
  call invert_mass_matrix_init()

#ifdef USE_MPI
  if ( nproc > 1 ) then

    ! preparing for MPI communications
    allocate(mask_ispec_inner_outer(nspec))
    mask_ispec_inner_outer(:) = .false.

    call get_MPI()

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

    ! buffers for MPI communications
    max_ibool_interfaces_size_ac = maxval(nibool_interfaces_acoustic(:))
    max_ibool_interfaces_size_el = 3*maxval(nibool_interfaces_elastic(:))
    max_ibool_interfaces_size_po = NDIM*maxval(nibool_interfaces_poroelastic(:))
    max_nibool_interfaces_ext_mesh = maxval(nibool_interfaces_ext_mesh(:))
      allocate(tab_requests_send_recv_acoustic(ninterface_acoustic*2))
      allocate(buffer_send_faces_vector_ac(max_ibool_interfaces_size_ac,ninterface_acoustic))
      allocate(buffer_recv_faces_vector_ac(max_ibool_interfaces_size_ac,ninterface_acoustic))
      allocate(tab_requests_send_recv_elastic(ninterface_elastic*2))
      allocate(buffer_send_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic))
      allocate(buffer_recv_faces_vector_el(max_ibool_interfaces_size_el,ninterface_elastic))
      allocate(tab_requests_send_recv_poro(ninterface_poroelastic*4))
      allocate(buffer_send_faces_vector_pos(max_ibool_interfaces_size_po,ninterface_poroelastic))
      allocate(buffer_recv_faces_vector_pos(max_ibool_interfaces_size_po,ninterface_poroelastic))
      allocate(buffer_send_faces_vector_pow(max_ibool_interfaces_size_po,ninterface_poroelastic))
      allocate(buffer_recv_faces_vector_pow(max_ibool_interfaces_size_po,ninterface_poroelastic))

! assembling the mass matrix
    call assemble_MPI_scalar(rmass_inverse_acoustic,nglob_acoustic, &
                            rmass_inverse_elastic_one,rmass_inverse_elastic_three,nglob_elastic, &
                            rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic,nglob_poroelastic)

  else
    ninterface_acoustic = 0
    ninterface_elastic = 0
    ninterface_poroelastic = 0

    num_ispec_outer = 0
    num_ispec_inner = 0
    allocate(mask_ispec_inner_outer(1))

    nspec_outer = 0
    nspec_inner = nspec

    allocate(ispec_inner_to_glob(nspec_inner))
    do ispec = 1, nspec
      ispec_inner_to_glob(ispec) = ispec
    enddo

  endif ! end of test on whether there is more than one process (nproc > 1)

#else
  num_ispec_outer = 0
  num_ispec_inner = 0
  allocate(mask_ispec_inner_outer(1))

  nspec_outer = 0
  nspec_inner = nspec

  allocate(ispec_outer_to_glob(1))
  allocate(ispec_inner_to_glob(nspec_inner))
  do ispec = 1, nspec
     ispec_inner_to_glob(ispec) = ispec
  enddo

#endif

    ! loop over spectral elements
    do ispec_outer = 1,nspec_outer
    ! get global numbering for inner or outer elements
      ispec = ispec_outer_to_glob(ispec_outer)
    enddo

    ! loop over spectral elements
    do ispec_inner = 1,nspec_inner
    ! get global numbering for inner or outer elements
      ispec = ispec_inner_to_glob(ispec_inner)
    enddo

    allocate(ibool_outer(NGLLX,NGLLZ,nspec_outer))
    allocate(ibool_inner(NGLLX,NGLLZ,nspec_inner))
    allocate(ispec_is_inner(nspec))
    ispec_is_inner(:) = .false.

    ! loop over spectral elements
    do ispec_outer = 1,nspec_outer
    ! get global numbering for inner or outer elements
      ispec = ispec_outer_to_glob(ispec_outer)
      ibool_outer(:,:,ispec_outer) = ibool(:,:,ispec)
    enddo

    ! loop over spectral elements
    do ispec_inner = 1,nspec_inner
    ! get global numbering for inner or outer elements
      ispec = ispec_inner_to_glob(ispec_inner)
      ibool_inner(:,:,ispec_inner) = ibool(:,:,ispec)
      ispec_is_inner(ispec) = .true.
    enddo

    ! reduces cache misses for outer elements
    call get_global_indirect_addressing(nspec_outer,nglob,ibool_outer,copy_ibool_ori,integer_mask_ibool)

    ! the total number of points without multiples in this region is now known
    nglob_outer = maxval(ibool_outer)

    ! reduces cache misses for inner elements
    call get_global_indirect_addressing(nspec_inner,nglob,ibool_inner,copy_ibool_ori,integer_mask_ibool)

    ! the total number of points without multiples in this region is now known
    nglob_inner = maxval(ibool_inner)

  call invert_mass_matrix()

! check the mesh, stability and number of points per wavelength
  if(DISPLAY_SUBSET_OPTION == 1) then
    UPPER_LIMIT_DISPLAY = nspec
  else if(DISPLAY_SUBSET_OPTION == 2) then
    UPPER_LIMIT_DISPLAY = nspec_inner
  else if(DISPLAY_SUBSET_OPTION == 3) then
    UPPER_LIMIT_DISPLAY = nspec_outer
  else if(DISPLAY_SUBSET_OPTION == 4) then
    UPPER_LIMIT_DISPLAY = NSPEC_DISPLAY_SUBSET
  else
    stop 'incorrect value of DISPLAY_SUBSET_OPTION'
  endif
  call checkgrid()

! convert receiver angle to radians
  anglerec = anglerec * pi / 180.d0

end subroutine prepare_timerun_mass_matrix





subroutine prepare_timerun_image_coloring()


#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

  integer i,j

#ifdef USE_MPI
  include "precision.h"
  integer k
#endif


!
!---- for color images
!

  if(output_color_image) then
    ! prepares dimension of image
    call prepare_color_image_init()

    ! allocate an array for image data
    allocate(image_color_data(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 1'
    allocate(image_color_vp_display(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 2'

    ! allocate an array for the grid point that corresponds to a given image data point
    allocate(iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 3'
    allocate(copy_iglob_image_color(NX_IMAGE_color,NZ_IMAGE_color),stat=ier); if(ier /= 0) stop 'error in an allocate statement 4'

    ! creates pixels indexing
    call prepare_color_image_pixels()

    ! creating and filling array num_pixel_loc with the positions of each colored
    ! pixel owned by the local process (useful for parallel jobs)
    allocate(num_pixel_loc(nb_pixel_loc))

    nb_pixel_loc = 0
    do i = 1, NX_IMAGE_color
       do j = 1, NZ_IMAGE_color
          if ( iglob_image_color(i,j) /= -1 ) then
             nb_pixel_loc = nb_pixel_loc + 1
             num_pixel_loc(nb_pixel_loc) = (j-1)*NX_IMAGE_color + i
          endif
       enddo
    enddo

! filling array iglob_image_color, containing info on which process owns which pixels.
#ifdef USE_MPI
    allocate(nb_pixel_per_proc(nproc))

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
                  iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
             do k = 1, nb_pixel_per_proc(iproc+1)
                j = ceiling(real(num_pixel_recv(k,iproc+1)) / real(NX_IMAGE_color))
                i = num_pixel_recv(k,iproc+1) - (j-1)*NX_IMAGE_color

                ! avoid edge effects
                if(i < 1) i = 1
                if(j < 1) j = 1

                if(i > NX_IMAGE_color) i = NX_IMAGE_color
                if(j > NZ_IMAGE_color) j = NZ_IMAGE_color

                iglob_image_color(i,j) = iproc

             enddo
          enddo

       else
          call MPI_SEND(num_pixel_loc(1),nb_pixel_loc,MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
       endif
    endif
#endif

    if (myrank == 0) write(IOUT,*) 'done locating all the pixels of color images'

  endif ! color_image


end subroutine prepare_timerun_image_coloring











subroutine prepare_timerun_kernel()


#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif

!
!----- Allocate sensitivity kernel arrays
!

  if(SIMULATION_TYPE == 3) then

    if(any_elastic) then

      if(.not. save_ASCII_kernels) then
        if (NEW_BINARY_FORMAT) then
          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kernel.bin'
          open(unit = 200, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ios)
          if (ios /= 0) stop 'Error writing kernel file to disk'

          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_kappa_kernel.bin'
          open(unit = 201, file ='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ios)
          if (ios /= 0) stop 'Error writing kernel file to disk'

          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mu_kernel.bin'
          open(unit = 202, file ='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ios)
          if (ios /= 0) stop 'Error writing kernel file to disk'

          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_kernel.bin'
          open(unit = 203, file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ios)
          if (ios /= 0) stop 'Error writing kernel file to disk'

          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_alpha_kernel.bin'
          open(unit = 204, file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ios)
          if (ios /= 0) stop 'Error writing kernel file to disk'

          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_beta_kernel.bin'
          open(unit = 205, file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ios)
          if (ios /= 0) stop 'Error writing kernel file to disk'

        else
          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_alpha_beta_kernel.bin'
          open(unit = 97, file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',iostat=ios)
          if (ios /= 0) stop 'Error writing kernel file to disk'

          write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_mu_kernel.bin'
          open(unit = 98, file='OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted', iostat=ios)
          if (ios /= 0) stop 'Error writing kernel file to disk'
        endif

      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_alpha_beta_kernel.dat'
        open(unit = 97, file = 'OUTPUT_FILES/'//outputname,status='unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'

        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_mu_kernel.dat'
        open(unit = 98, file = 'OUTPUT_FILES/'//outputname,status='unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      endif

      rho_kl(:,:,:) = 0._CUSTOM_REAL
      mu_kl(:,:,:) = 0._CUSTOM_REAL
      kappa_kl(:,:,:) = 0._CUSTOM_REAL

      rhop_kl(:,:,:) = 0._CUSTOM_REAL
      beta_kl(:,:,:) = 0._CUSTOM_REAL
      alpha_kl(:,:,:) = 0._CUSTOM_REAL
      rhorho_el_hessian_final2(:,:,:) = 0._CUSTOM_REAL
      rhorho_el_hessian_temp2(:) = 0._CUSTOM_REAL
      rhorho_el_hessian_final1(:,:,:) = 0._CUSTOM_REAL
      rhorho_el_hessian_temp1(:) = 0._CUSTOM_REAL
    endif

    if(any_poroelastic) then

      ! Primary kernels
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mu_B_C_kernel.dat'
      open(unit = 144, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_M_rho_rhof_kernel.dat'
      open(unit = 155, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_m_eta_kernel.dat'
      open(unit = 16, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      ! Wavespeed kernels
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_cpI_cpII_cs_kernel.dat'
      open(unit = 20, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhobb_rhofbb_ratio_kernel.dat'
      open(unit = 21, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_phib_eta_kernel.dat'
      open(unit = 22, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      ! Density normalized kernels
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mub_Bb_Cb_kernel.dat'
      open(unit = 17, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_Mb_rhob_rhofb_kernel.dat'
      open(unit = 18, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'
      write(outputname,'(a,i6.6,a)') 'proc',myrank,'_mb_etab_kernel.dat'
      open(unit = 19, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
      if (ios /= 0) stop 'Error writing kernel file to disk'

      rhot_kl(:,:,:) = 0._CUSTOM_REAL
      rhof_kl(:,:,:) = 0._CUSTOM_REAL
      eta_kl(:,:,:) = 0._CUSTOM_REAL
      sm_kl(:,:,:) = 0._CUSTOM_REAL
      mufr_kl(:,:,:) = 0._CUSTOM_REAL
      B_kl(:,:,:) = 0._CUSTOM_REAL
      C_kl(:,:,:) = 0._CUSTOM_REAL
      M_kl(:,:,:) = 0._CUSTOM_REAL

      rhob_kl(:,:,:) = 0._CUSTOM_REAL
      rhofb_kl(:,:,:) = 0._CUSTOM_REAL
      phi_kl(:,:,:) = 0._CUSTOM_REAL
      mufrb_kl(:,:,:) = 0._CUSTOM_REAL
      Bb_kl(:,:,:) = 0._CUSTOM_REAL
      Cb_kl(:,:,:) = 0._CUSTOM_REAL
      Mb_kl(:,:,:) = 0._CUSTOM_REAL

      rhobb_kl(:,:,:) = 0._CUSTOM_REAL
      rhofbb_kl(:,:,:) = 0._CUSTOM_REAL
      phib_kl(:,:,:) = 0._CUSTOM_REAL
      cs_kl(:,:,:) = 0._CUSTOM_REAL
      cpI_kl(:,:,:) = 0._CUSTOM_REAL
      cpII_kl(:,:,:) = 0._CUSTOM_REAL
      ratio_kl(:,:,:) = 0._CUSTOM_REAL
    endif

    if(any_acoustic) then

      if(.not. save_ASCII_kernels)then
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_kernel.bin'
        open(unit = 95, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',&
            iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rho_kappa_kernel.dat'
        open(unit = 95, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      endif

      if(.not. save_ASCII_kernels)then
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_c_kernel.bin'
        open(unit = 96, file = 'OUTPUT_FILES/'//outputname,status='unknown',action='write',form='unformatted',&
            iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      else
        write(outputname,'(a,i6.6,a)') 'proc',myrank,'_rhop_c_kernel.dat'
        open(unit = 96, file = 'OUTPUT_FILES/'//outputname,status = 'unknown',iostat=ios)
        if (ios /= 0) stop 'Error writing kernel file to disk'
      endif

      rho_ac_kl(:,:,:) = 0._CUSTOM_REAL
      kappa_ac_kl(:,:,:) = 0._CUSTOM_REAL

      rhop_ac_kl(:,:,:) = 0._CUSTOM_REAL
      alpha_ac_kl(:,:,:) = 0._CUSTOM_REAL
      rhorho_ac_hessian_final2(:,:,:) = 0._CUSTOM_REAL
      rhorho_ac_hessian_final1(:,:,:) = 0._CUSTOM_REAL
    endif

  endif ! if(SIMULATION_TYPE == 3)


end subroutine prepare_timerun_kernel









subroutine prepare_timerun_pml()


#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

  integer i

#ifdef USE_MPI
  include "precision.h"
#endif

if (GPU_MODE .and. PML_BOUNDARY_CONDITIONS ) stop 'error : PML not implemented on GPU mode. Please use Stacey instead'

  ! PML absorbing conditions
    anyabs_glob=anyabs
#ifdef USE_MPI
    call MPI_ALLREDUCE(anyabs, anyabs_glob, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ier)
#endif

    ! allocate this in all cases, even if PML is not set, because we use it for PostScript display as well
    allocate(is_PML(nspec),stat=ier)
    if(ier /= 0) stop 'error: not enough memory to allocate array is_PML'
    is_PML(:) = .false.

    if(PML_BOUNDARY_CONDITIONS .and. anyabs_glob ) then
      allocate(spec_to_PML(nspec),stat=ier)
      if(ier /= 0) stop 'error: not enough memory to allocate array spec_to_PML'

      allocate(which_PML_elem(4,nspec),stat=ier)
      if(ier /= 0) stop 'error: not enough memory to allocate array which_PML_elem'
      which_PML_elem(:,:) = .false.

      if(SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))then
        allocate(PML_interior_interface(4,nspec),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array PML_interior_interface'
        PML_interior_interface = .false.
      else
        allocate(PML_interior_interface(4,1))
      endif

! add support for using PML in MPI mode with external mesh
      if(read_external_mesh)then
        allocate(mask_ibool_pml(nglob))
      else
        allocate(mask_ibool_pml(1))
      endif

      call pml_init()

      if((SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD)) .and. PML_BOUNDARY_CONDITIONS)then

        if(nglob_interface > 0) then
          allocate(point_interface(nglob_interface),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array point_interface'
        endif

        if(any_elastic .and. nglob_interface > 0)then
          allocate(pml_interface_history_displ(3,nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_displ'
          allocate(pml_interface_history_veloc(3,nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_veloc'
          allocate(pml_interface_history_accel(3,nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_accel'
        endif

        if(any_acoustic .and. nglob_interface > 0)then
          allocate(pml_interface_history_potential(nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_potential'
          allocate(pml_interface_history_potential_dot(nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_potential_dot'
          allocate(pml_interface_history_potential_dot_dot(nglob_interface,NSTEP),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array pml_interface_history_potential_dot_dot'
        endif

        if(nglob_interface > 0) then
          call determin_interface_pml_interior()
          deallocate(PML_interior_interface)
          deallocate(mask_ibool_pml)
        endif

        if(any_elastic .and. nglob_interface > 0)then
          write(outputname,'(a,i6.6,a)') 'pml_interface_elastic',myrank,'.bin'
          open(unit=71,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
        endif

        if(any_acoustic .and. nglob_interface > 0)then
          write(outputname,'(a,i6.6,a)') 'pml_interface_acoustic',myrank,'.bin'
          open(unit=72,file='OUTPUT_FILES/'//outputname,status='unknown',form='unformatted')
        endif
      else
        allocate(point_interface(1))
        allocate(pml_interface_history_displ(3,1,1))
        allocate(pml_interface_history_veloc(3,1,1))
        allocate(pml_interface_history_accel(3,1,1))
        allocate(pml_interface_history_potential(1,1))
        allocate(pml_interface_history_potential_dot(1,1))
        allocate(pml_interface_history_potential_dot_dot(1,1))
      endif

      if(SIMULATION_TYPE == 3 .and. PML_BOUNDARY_CONDITIONS)then

        if(any_elastic .and. nglob_interface > 0)then
          do it = 1,NSTEP
            do i = 1, nglob_interface
              read(71)pml_interface_history_accel(1,i,it),pml_interface_history_accel(2,i,it),&
                      pml_interface_history_accel(3,i,it),&
                      pml_interface_history_veloc(1,i,it),pml_interface_history_veloc(2,i,it),&
                      pml_interface_history_veloc(3,i,it),&
                      pml_interface_history_displ(1,i,it),pml_interface_history_displ(2,i,it),&
                      pml_interface_history_displ(3,i,it)
          enddo
         enddo
       endif

       if(any_acoustic .and. nglob_interface > 0)then
         do it = 1,NSTEP
           do i = 1, nglob_interface
             read(72)pml_interface_history_potential_dot_dot(i,it),pml_interface_history_potential_dot(i,it),&
                     pml_interface_history_potential(i,it)
           enddo
         enddo
       endif
     endif

      deallocate(which_PML_elem)

      if (nspec_PML==0) nspec_PML=1 ! DK DK added this

      if (nspec_PML > 0) then
        allocate(K_x_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array K_x_store'
        allocate(K_z_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array K_z_store'
        allocate(d_x_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array d_x_store'
        allocate(d_z_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array d_z_store'
        allocate(alpha_x_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array alpha_x_store'
        allocate(alpha_z_store(NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array alpha_z_store'
        K_x_store(:,:,:) = ZERO
        K_z_store(:,:,:) = ZERO
        d_x_store(:,:,:) = ZERO
        d_z_store(:,:,:) = ZERO
        alpha_x_store(:,:,:) = ZERO
        alpha_z_store(:,:,:) = ZERO
        call define_PML_coefficients()
      else
        allocate(K_x_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array K_x_store'
        allocate(K_z_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array K_z_store'
        allocate(d_x_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array d_x_store'
        allocate(d_z_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array d_z_store'
        allocate(alpha_x_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array alpha_x_store'
        allocate(alpha_z_store(NGLLX,NGLLZ,1),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array alpha_z_store'
        K_x_store(:,:,:) = ZERO
        K_z_store(:,:,:) = ZERO
        d_x_store(:,:,:) = ZERO
        d_z_store(:,:,:) = ZERO
        alpha_x_store(:,:,:) = ZERO
        alpha_z_store(:,:,:) = ZERO
      endif

      ! elastic PML memory variables
      if (any_elastic .and. nspec_PML>0) then
        allocate(rmemory_displ_elastic(2,3,NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_displ_elastic'
        allocate(rmemory_dux_dx(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx'
        allocate(rmemory_dux_dz(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz'
        allocate(rmemory_duz_dx(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx'
        allocate(rmemory_duz_dz(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz'
        if(any_acoustic .and. num_fluid_solid_edges > 0)then
          allocate(rmemory_fsb_displ_elastic(1,3,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_fsb_displ_elastic'
          allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_sfb_potential_ddot_acoustic'
        endif

        if(ROTATE_PML_ACTIVATE)then
          allocate(rmemory_dux_dx_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx_prime'
          allocate(rmemory_dux_dz_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz_prime'
          allocate(rmemory_duz_dx_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx_prime'
          allocate(rmemory_duz_dz_prime(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz_prime'
        else
          allocate(rmemory_dux_dx_prime(1,1,1,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx_prime'
          allocate(rmemory_dux_dz_prime(1,1,1,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz_prime'
          allocate(rmemory_duz_dx_prime(1,1,1,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx_prime'
          allocate(rmemory_duz_dz_prime(1,1,1,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz_prime'
        endif

        if(time_stepping_scheme == 2)then
          allocate(rmemory_displ_elastic_LDDRK(2,3,NGLLX,NGLLZ,nspec_PML),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_displ_elastic'
          allocate(rmemory_dux_dx_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dx'
          allocate(rmemory_dux_dz_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_dux_dz'
          allocate(rmemory_duz_dx_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dx'
          allocate(rmemory_duz_dz_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_duz_dz'
          if(any_acoustic .and. num_fluid_solid_edges > 0)then
            allocate(rmemory_fsb_displ_elastic_LDDRK(1,3,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
            if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_fsb_displ_elastic'
            allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,num_fluid_solid_edges),stat=ier)
            if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_sfb_potential_ddot_acoustic'
          endif
        else
          allocate(rmemory_displ_elastic_LDDRK(1,1,1,1,1),stat=ier)
          allocate(rmemory_dux_dx_LDDRK(1,1,1,2),stat=ier)
          allocate(rmemory_dux_dz_LDDRK(1,1,1,2),stat=ier)
          allocate(rmemory_duz_dx_LDDRK(1,1,1,2),stat=ier)
          allocate(rmemory_duz_dz_LDDRK(1,1,1,2),stat=ier)
          if(any_acoustic .and. num_fluid_solid_edges > 0)then
            allocate(rmemory_fsb_displ_elastic_LDDRK(1,3,NGLLX,NGLLZ,1),stat=ier)
            if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_fsb_displ_elastic'
            allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1),stat=ier)
            if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_sfb_potential_ddot_acoustic'
          endif
        endif

        rmemory_displ_elastic(:,:,:,:,:) = ZERO
        rmemory_dux_dx(:,:,:,:) = ZERO
        rmemory_dux_dz(:,:,:,:) = ZERO
        rmemory_duz_dx(:,:,:,:) = ZERO
        rmemory_duz_dz(:,:,:,:) = ZERO

        if(any_acoustic .and. num_fluid_solid_edges > 0)then
          rmemory_fsb_displ_elastic(:,:,:,:,:) = ZERO
          rmemory_sfb_potential_ddot_acoustic(:,:,:,:) = ZERO
        endif

        if(ROTATE_PML_ACTIVATE)then
          rmemory_dux_dx_prime(:,:,:,:) = ZERO
          rmemory_dux_dz_prime(:,:,:,:) = ZERO
          rmemory_duz_dx_prime(:,:,:,:) = ZERO
          rmemory_duz_dz_prime(:,:,:,:) = ZERO
        endif

        if(time_stepping_scheme == 2)then
          rmemory_displ_elastic_LDDRK(:,:,:,:,:) = ZERO
          rmemory_dux_dx_LDDRK(:,:,:,:) = ZERO
          rmemory_dux_dz_LDDRK(:,:,:,:) = ZERO
          rmemory_duz_dx_LDDRK(:,:,:,:) = ZERO
          rmemory_duz_dz_LDDRK(:,:,:,:) = ZERO
          if(any_acoustic .and. num_fluid_solid_edges > 0)then
            rmemory_fsb_displ_elastic_LDDRK(:,:,:,:,:) = ZERO
            rmemory_sfb_potential_ddot_acoustic_LDDRK(:,:,:,:) = ZERO
          endif
        endif

      else

        allocate(rmemory_displ_elastic(1,1,1,1,1))
        allocate(rmemory_dux_dx(1,1,1,1))
        allocate(rmemory_dux_dz(1,1,1,1))
        allocate(rmemory_duz_dx(1,1,1,1))
        allocate(rmemory_duz_dz(1,1,1,1))
        if(any_acoustic .and. num_fluid_solid_edges > 0)then
          allocate(rmemory_fsb_displ_elastic(1,3,NGLLX,NGLLZ,1))
          allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,1))
          allocate(rmemory_fsb_displ_elastic_LDDRK(1,3,NGLLX,NGLLZ,1))
          allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1))
        endif

        allocate(rmemory_dux_dx_prime(1,1,1,1))
        allocate(rmemory_dux_dz_prime(1,1,1,1))
        allocate(rmemory_duz_dx_prime(1,1,1,1))
        allocate(rmemory_duz_dz_prime(1,1,1,1))

        allocate(rmemory_displ_elastic_LDDRK(1,1,1,1,1))
        allocate(rmemory_dux_dx_LDDRK(1,1,1,1))
        allocate(rmemory_dux_dz_LDDRK(1,1,1,1))
        allocate(rmemory_duz_dx_LDDRK(1,1,1,1))
        allocate(rmemory_duz_dz_LDDRK(1,1,1,1))
      endif

      if (any_acoustic .and. nspec_PML>0) then
        allocate(rmemory_potential_acoustic(2,NGLLX,NGLLZ,nspec_PML),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_potential_acoustic'
        allocate(rmemory_acoustic_dux_dx(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dx'
        allocate(rmemory_acoustic_dux_dz(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
        if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dz'

        rmemory_potential_acoustic = ZERO
        rmemory_acoustic_dux_dx = ZERO
        rmemory_acoustic_dux_dz = ZERO

        if(time_stepping_scheme == 2)then
          allocate(rmemory_potential_acoustic_LDDRK(2,NGLLX,NGLLZ,nspec_PML),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_potential_acoustic'
          allocate(rmemory_acoustic_dux_dx_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dx'
          allocate(rmemory_acoustic_dux_dz_LDDRK(NGLLX,NGLLZ,nspec_PML,2),stat=ier)
          if(ier /= 0) stop 'error: not enough memory to allocate array rmemory_acoustic_dux_dz'
        else
          allocate(rmemory_potential_acoustic_LDDRK(1,1,1,1),stat=ier)
          allocate(rmemory_acoustic_dux_dx_LDDRK(1,1,1,1),stat=ier)
          allocate(rmemory_acoustic_dux_dz_LDDRK(1,1,1,1),stat=ier)
        endif

        rmemory_potential_acoustic_LDDRK = ZERO
        rmemory_acoustic_dux_dx_LDDRK = ZERO
        rmemory_acoustic_dux_dz_LDDRK = ZERO

      else
        allocate(rmemory_potential_acoustic(1,1,1,1))
        allocate(rmemory_acoustic_dux_dx(1,1,1,1))
        allocate(rmemory_acoustic_dux_dz(1,1,1,1))
      endif

    else
      allocate(rmemory_dux_dx(1,1,1,1))
      allocate(rmemory_dux_dz(1,1,1,1))
      allocate(rmemory_duz_dx(1,1,1,1))
      allocate(rmemory_duz_dz(1,1,1,1))
      allocate(rmemory_fsb_displ_elastic(1,3,NGLLX,NGLLZ,1))
      allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,1))
      allocate(rmemory_fsb_displ_elastic_LDDRK(1,3,NGLLX,NGLLZ,1))
      allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1))

      allocate(rmemory_dux_dx_prime(1,1,1,1))
      allocate(rmemory_dux_dz_prime(1,1,1,1))
      allocate(rmemory_duz_dx_prime(1,1,1,1))
      allocate(rmemory_duz_dz_prime(1,1,1,1))

      allocate(rmemory_displ_elastic(1,1,1,1,1))

      allocate(rmemory_displ_elastic_LDDRK(1,1,1,1,1))
      allocate(rmemory_dux_dx_LDDRK(1,1,1,1))
      allocate(rmemory_dux_dz_LDDRK(1,1,1,1))
      allocate(rmemory_duz_dx_LDDRK(1,1,1,1))
      allocate(rmemory_duz_dz_LDDRK(1,1,1,1))

      allocate(rmemory_potential_acoustic(1,1,1,1))
      allocate(rmemory_acoustic_dux_dx(1,1,1,1))
      allocate(rmemory_acoustic_dux_dz(1,1,1,1))

      allocate(rmemory_potential_acoustic_LDDRK(1,1,1,1))
      allocate(rmemory_acoustic_dux_dx_LDDRK(1,1,1,1))
      allocate(rmemory_acoustic_dux_dz_LDDRK(1,1,1,1))

      allocate(spec_to_PML(1))

      allocate(K_x_store(1,1,1))
      allocate(K_z_store(1,1,1))
      allocate(d_x_store(1,1,1))
      allocate(d_z_store(1,1,1))
      allocate(alpha_x_store(1,1,1))
      allocate(alpha_z_store(1,1,1))
    endif ! PML_BOUNDARY_CONDITIONS


  ! avoid a potential side effect owing to the "if" statements above: this array may be unallocated,
  ! if so we need to allocate a dummy version in order to be able to use that array as an argument
  ! in some subroutine calls below
  if(.not. allocated(rmemory_fsb_displ_elastic)) allocate(rmemory_fsb_displ_elastic(1,3,NGLLX,NGLLZ,1))
  if(.not. allocated(rmemory_sfb_potential_ddot_acoustic)) allocate(rmemory_sfb_potential_ddot_acoustic(1,NGLLX,NGLLZ,1))
  if(.not. allocated(rmemory_fsb_displ_elastic_LDDRK)) then
    allocate(rmemory_fsb_displ_elastic_LDDRK(1,3,NGLLX,NGLLZ,1))
  endif
  if(.not. allocated(rmemory_sfb_potential_ddot_acoustic_LDDRK)) then
    allocate(rmemory_sfb_potential_ddot_acoustic_LDDRK(1,NGLLX,NGLLZ,1))
  endif


end subroutine prepare_timerun_pml









subroutine prepare_timerun_read()


#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

  integer i,ispec,ispec2,j

#ifdef USE_MPI
  include "precision.h"
#endif




  ! starts reading in Database file
  call read_databases_init()

  if(nproc_read_from_database < 1) stop 'should have nproc_read_from_database >= 1'
  if(SIMULATION_TYPE == 3 .and.(time_stepping_scheme == 2 .or. time_stepping_scheme == 3)) &
                                  stop 'RK and LDDRK time scheme not supported for adjoint inversion'
  if(nproc /= nproc_read_from_database) stop 'must always have nproc == nproc_read_from_database'

! add a small crack (discontinuity) in the medium manually
  npgeo_ori = npgeo
  if(ADD_A_SMALL_CRACK_IN_THE_MEDIUM) npgeo = npgeo + NB_POINTS_TO_ADD_TO_NPGEO

  !
  !--- source information
  !
    allocate( source_type(NSOURCES) )
    allocate( time_function_type(NSOURCES) )
    allocate( x_source(NSOURCES) )
    allocate( z_source(NSOURCES) )
    allocate( ix_image_color_source(NSOURCES) )
    allocate( iy_image_color_source(NSOURCES) )
    allocate( f0(NSOURCES) )
    allocate( tshift_src(NSOURCES) )
    allocate( factor(NSOURCES) )
    allocate( anglesource(NSOURCES) )
    allocate( Mxx(NSOURCES) )
    allocate( Mxz(NSOURCES) )
    allocate( Mzz(NSOURCES) )
    allocate( aval(NSOURCES) )
    allocate( ispec_selected_source(NSOURCES) )
    allocate( iglob_source(NSOURCES) )
    allocate( source_courbe_eros(NSOURCES) )
    allocate( xi_source(NSOURCES) )
    allocate( gamma_source(NSOURCES) )
    allocate( is_proc_source(NSOURCES) )
    allocate( nb_proc_source(NSOURCES) )
    allocate( sourcearray(NSOURCES,NDIM,NGLLX,NGLLZ) )

  ! reads in source infos
  call read_databases_sources()

  !if(AXISYM) factor = factor/(TWO*PI)   !!!!!axisym TODO verify

  ! sets source parameters
  call set_sources()

  !----  define time stepping scheme
  if(time_stepping_scheme == 1)then
    stage_time_scheme=1
  else if(time_stepping_scheme == 2)then
    stage_time_scheme=Nstages
  else if(time_stepping_scheme == 3)then
    stage_time_scheme=4
  endif

  !----  read attenuation information
  call read_databases_atten()

  ! if source is not a Dirac or Heavyside then f0_attenuation is f0 of the first source
  if(.not. (time_function_type(1) == 4 .or. time_function_type(1) == 5)) then
    f0_attenuation = f0(1)
  endif

  !---- read the spectral macrobloc nodal coordinates
  allocate(coorg(NDIM,npgeo))

  ! reads the spectral macrobloc nodal coordinates
  ! and basic properties of the spectral elements
  !! DK DK  call read_databases_coorg_elem(myrank,npgeo,coorg,numat,ngnod,nspec, &
  !! DK DK  added a crack manually
  call read_databases_coorg_elem()

  !---- allocate arrays
    allocate(shape2D(ngnod,NGLLX,NGLLZ))
    allocate(dershape2D(NDIM,ngnod,NGLLX,NGLLZ))
    allocate(shape2D_display(ngnod,pointsdisp,pointsdisp))
    allocate(dershape2D_display(NDIM,ngnod,pointsdisp,pointsdisp))
    if(AXISYM) then
      allocate(flagrange_GLJ(NGLJ,pointsdisp))
    else
      allocate(flagrange_GLJ(1,1))
    endif
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
    allocate(density(2,numat))
    allocate(anisotropy(9,numat))
    allocate(porosity(numat))
    allocate(tortuosity(numat))
    allocate(permeability(3,numat))
    allocate(poroelastcoef(4,3,numat))
    allocate(already_shifted_velocity(numat))
    allocate(QKappa_attenuation(numat))
    allocate(Qmu_attenuation(numat))
    allocate(kmato(nspec))
    allocate(knods(ngnod,nspec))
    allocate(ibool(NGLLX,NGLLZ,nspec))
    allocate(elastic(nspec))
    allocate(acoustic(nspec))
    allocate(gravitoacoustic(nspec))
    allocate(poroelastic(nspec))
    allocate(anisotropic(nspec))
    allocate(inv_tau_sigma_nu1(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(inv_tau_sigma_nu2(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(phi_nu1(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(phi_nu2(NGLLX,NGLLZ,nspec,N_SLS))
    allocate(tau_epsilon_nu1(N_SLS))
    allocate(tau_epsilon_nu2(N_SLS))
    allocate(inv_tau_sigma_nu1_sent(N_SLS))
    allocate(inv_tau_sigma_nu2_sent(N_SLS))
    allocate(phi_nu1_sent(N_SLS))
    allocate(phi_nu2_sent(N_SLS))

    already_shifted_velocity(:) = .false.

  !
  !---- read the material properties
  !
  call gmat01(f0(1))
  !
  !----  read spectral macrobloc data
  !

! add support for using PML in MPI mode with external mesh
  allocate(region_CPML(nspec))
  call read_databases_mato()

! add a small crack (discontinuity) in the medium manually
  if(ADD_A_SMALL_CRACK_IN_THE_MEDIUM) then

#ifdef USE_MPI
  stop 'currently only serial runs are handled when adding a crack manually'
#endif
!! DK DK material number 2 indicates the spectral elements that form the left vertical side of the crack
  check_nb_points_to_add_to_npgeo = count(kmato == 2)
  print *
  print *,'adding a crack manually'
  print *,'need to add ',nb_points_to_add_to_npgeo,' npgeo mesh points to do that'

  if(check_nb_points_to_add_to_npgeo /= NB_POINTS_TO_ADD_TO_NPGEO) &
    stop 'must have check_nb_points_to_add_to_npgeo == NB_POINTS_TO_ADD_TO_NPGEO when adding a crack manually'

  if(ngnod /= 4) stop 'must currently have ngnod == 4 when adding a crack manually'

  if(FAST_NUMBERING) stop 'must not have FAST_NUMBERING when adding a crack manually'

!! DK DK modify arrays "knods" and "coorg" to introduce the crack manually by duplicating and splitting the nodes
  already_found_a_crack_element = .false.
  current_last_point = npgeo_ori

  do ispec = 1,nspec-1
!! DK DK my convention is to introduce a vertical crack between two elements with material numbers 2 and 3
    if(kmato(ispec) == 2 .and. kmato(ispec+1) == 3) then

      print *,'adding a crack between elements ',ispec,' and ',ispec+1

!! DK DK duplicate and split the lower-right corner of this element,
!! DK DK except if it is the first crack element found, because then it is the crack
!! DK DK tip and thus it should be assembled rather than split.
!! DK DK Lower-right corner of an element is local npgeo point #2
      if(already_found_a_crack_element .and. knods(2,ispec) <= npgeo_ori) then
        current_last_point = current_last_point + 1
        original_value = knods(2,ispec)
!! DK DK split this point number in all the elements in which it appears
        do ispec2 = 1,nspec
! do this only for elements that define the left vertical edge of the crack
          if(kmato(ispec2) /= 2) cycle
          do ignod = 1,ngnod
            if(knods(ignod,ispec2) == original_value) then
              knods(ignod,ispec2) = current_last_point
              coorg(:,current_last_point) = coorg(:,original_value)
            endif
          enddo
        enddo
      endif

!! DK DK duplicate and split the upper-right corner of this element
      already_found_a_crack_element = .true.

!! DK DK Upper-right corner of an element is local npgeo point #3
      if(knods(3,ispec) <= npgeo_ori) then

        current_last_point = current_last_point + 1
        original_value = knods(3,ispec)
!! DK DK split this point number in all the elements in which it appears
        do ispec2 = 1,nspec
! do this only for elements that define the left vertical edge of the crack
          if(kmato(ispec2) /= 2) cycle
          do ignod = 1,ngnod
            if(knods(ignod,ispec2) == original_value) then
              knods(ignod,ispec2) = current_last_point
              coorg(:,current_last_point) = coorg(:,original_value)
            endif
          enddo
        enddo
      endif

    endif ! of if(kmato(ispec) == 2 .and. kmato(ispec+1) == 3)

  enddo

  if(current_last_point /= npgeo) then
    print *,'current_last_point = ',current_last_point
    print *,'npgeo_new = ',npgeo
    stop 'did not find the right total number of points, should have current_last_point == npgeo_new'
  endif

  endif ! of if(ADD_A_SMALL_CRACK_IN_THE_MEDIUM) then

!-------------------------------------------------------------------------------
!----  determine if each spectral element is elastic, poroelastic, or acoustic
!-------------------------------------------------------------------------------
  call initialize_simulation_domains()

  if(PML_BOUNDARY_CONDITIONS .and. any_poroelastic) then
    stop 'PML boundary conditions not implemented for poroelastic simulations yet'
  endif

  if(PML_BOUNDARY_CONDITIONS .and. any_elastic .and. (.not. p_sv)) then
    stop 'PML boundary conditions not implemented for SH simulations yet'
  endif

  if(PML_BOUNDARY_CONDITIONS .and. time_stepping_scheme == 3) then
    stop 'PML boundary conditions not implemented with standard Runge Kutta scheme'
  endif

#ifdef USE_MPI
  if(myrank == 0)then
   if(time_stepping_scheme == 3) then
    stop 'MPI support for standard Runge-Kutta scheme is not implemented'
   endif
  endif
#endif

  ! allocate memory variables for attenuation
    allocate(e1(NGLLX,NGLLZ,nspec_allocate,N_SLS))
    allocate(e11(NGLLX,NGLLZ,nspec_allocate,N_SLS))
    allocate(e13(NGLLX,NGLLZ,nspec_allocate,N_SLS))

    e1(:,:,:,:) = 0._CUSTOM_REAL
    e11(:,:,:,:) = 0._CUSTOM_REAL
    e13(:,:,:,:) = 0._CUSTOM_REAL

    if(time_stepping_scheme == 2)then
      allocate(e1_LDDRK(NGLLX,NGLLZ,nspec_allocate,N_SLS))
      allocate(e11_LDDRK(NGLLX,NGLLZ,nspec_allocate,N_SLS))
      allocate(e13_LDDRK(NGLLX,NGLLZ,nspec_allocate,N_SLS))
    else
      allocate(e1_LDDRK(1,1,1,1))
      allocate(e11_LDDRK(1,1,1,1))
      allocate(e13_LDDRK(1,1,1,1))
    endif
    e1_LDDRK(:,:,:,:) = 0._CUSTOM_REAL
    e11_LDDRK(:,:,:,:) = 0._CUSTOM_REAL
    e13_LDDRK(:,:,:,:) = 0._CUSTOM_REAL

    if(time_stepping_scheme == 3)then
      allocate(e1_initial_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS))
      allocate(e11_initial_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS))
      allocate(e13_initial_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS))
      allocate(e1_force_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS,stage_time_scheme))
      allocate(e11_force_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS,stage_time_scheme))
      allocate(e13_force_rk(NGLLX,NGLLZ,nspec_allocate,N_SLS,stage_time_scheme))
    else
      allocate(e1_initial_rk(1,1,1,1))
      allocate(e11_initial_rk(1,1,1,1))
      allocate(e13_initial_rk(1,1,1,1))
      allocate(e1_force_rk(1,1,1,1,1))
      allocate(e11_force_rk(1,1,1,1,1))
      allocate(e13_force_rk(1,1,1,1,1))
    endif
    e1_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
    e11_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
    e13_initial_rk(:,:,:,:) = 0._CUSTOM_REAL
    e1_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
    e11_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
    e13_force_rk(:,:,:,:,:) = 0._CUSTOM_REAL
    allocate(Mu_nu1(NGLLX,NGLLZ,nspec))
    allocate(Mu_nu2(NGLLX,NGLLZ,nspec))

! initialize to dummy values
! convention to indicate that Q = 9999 in that element i.e. that there is no viscoelasticity in that element
  inv_tau_sigma_nu1(:,:,:,:) = -1._CUSTOM_REAL
  phi_nu1(:,:,:,:) = -1._CUSTOM_REAL
  inv_tau_sigma_nu2(:,:,:,:) = -1._CUSTOM_REAL
  phi_nu2(:,:,:,:) = -1._CUSTOM_REAL
  Mu_nu1(:,:,:) = -1._CUSTOM_REAL
  Mu_nu2(:,:,:) = -1._CUSTOM_REAL

! define the attenuation quality factors.
! they can be different for each element.
!! DK DK if needed in the future, here the quality factor could be different for each point
  do ispec = 1,nspec

!   attenuation is not implemented in acoustic (i.e. fluid) media for now, only in viscoelastic (i.e. solid) media
    if(acoustic(ispec)) cycle

!   check that attenuation values entered by the user make sense
    if((QKappa_attenuation(kmato(ispec)) <= 9998.999d0 .and. Qmu_attenuation(kmato(ispec)) >  9998.999d0) .or. &
       (QKappa_attenuation(kmato(ispec)) >  9998.999d0 .and. Qmu_attenuation(kmato(ispec)) <= 9998.999d0)) stop &
     'need to have Qkappa and Qmu both above or both below 9999 for a given material; trick: use 9998 if you want to turn off one'

!   if no attenuation in that elastic element
    if(QKappa_attenuation(kmato(ispec)) > 9998.999d0) cycle

    call attenuation_model(QKappa_attenuation(kmato(ispec)),Qmu_attenuation(kmato(ispec)))

    do j = 1,NGLLZ
      do i = 1,NGLLX
        inv_tau_sigma_nu1(i,j,ispec,:) = inv_tau_sigma_nu1_sent(:)
        phi_nu1(i,j,ispec,:) = phi_nu1_sent(:)
        inv_tau_sigma_nu2(i,j,ispec,:) = inv_tau_sigma_nu2_sent(:)
        phi_nu2(i,j,ispec,:) = phi_nu2_sent(:)
        Mu_nu1(i,j,ispec) = Mu_nu1_sent
        Mu_nu2(i,j,ispec) = Mu_nu2_sent
      enddo
    enddo

    if(ATTENUATION_VISCOELASTIC_SOLID .and. READ_VELOCITIES_AT_F0 .and. .not. assign_external_model) then
      if(anisotropic(ispec) .or. poroelastic(ispec) .or. gravitoacoustic(ispec)) &
         stop 'READ_VELOCITIES_AT_F0 only implemented for non anisotropic, non poroelastic, non gravitoacoustic materials for now'
      n = kmato(ispec)
      if(.not. already_shifted_velocity(n)) then
        rho = density(1,n)
        lambda = poroelastcoef(1,1,n)
        mu = poroelastcoef(2,1,n)
        vp = dsqrt((lambda + TWO * mu) / rho)
        vs = dsqrt(mu / rho)
        call shift_velocities_from_f0(vp,vs,rho,mu,lambda)
        poroelastcoef(1,1,n) = lambda
        poroelastcoef(2,1,n) = mu
        poroelastcoef(3,1,n) = lambda + TWO*mu
        already_shifted_velocity(n) = .true.
      endif
    endif

 enddo

! allocate memory variables for viscous attenuation (poroelastic media)
    if(ATTENUATION_PORO_FLUID_PART) then
      allocate(rx_viscous(NGLLX,NGLLZ,nspec))
      allocate(rz_viscous(NGLLX,NGLLZ,nspec))
      allocate(viscox(NGLLX,NGLLZ,nspec))
      allocate(viscoz(NGLLX,NGLLZ,nspec))

      if(time_stepping_scheme == 2) then
      allocate(rx_viscous_LDDRK(NGLLX,NGLLZ,nspec))
      allocate(rz_viscous_LDDRK(NGLLX,NGLLZ,nspec))
      endif

      if(time_stepping_scheme == 3) then
      allocate(rx_viscous_initial_rk(NGLLX,NGLLZ,nspec))
      allocate(rz_viscous_initial_rk(NGLLX,NGLLZ,nspec))
      allocate(rx_viscous_force_RK(NGLLX,NGLLZ,nspec,stage_time_scheme))
      allocate(rz_viscous_force_RK(NGLLX,NGLLZ,nspec,stage_time_scheme))
      endif

    else
      allocate(rx_viscous(NGLLX,NGLLZ,1))
      allocate(rz_viscous(NGLLX,NGLLZ,1))
      allocate(viscox(NGLLX,NGLLZ,1))
      allocate(viscoz(NGLLX,NGLLZ,1))
    endif

  !
  !----  read interfaces data
  !
  call read_databases_ninterface()
  if ( ninterface > 0 ) then
       allocate(my_neighbours(ninterface))
       allocate(my_nelmnts_neighbours(ninterface))
       allocate(my_interfaces(4,max_interface_size,ninterface))
       allocate(ibool_interfaces_acoustic(NGLLX*max_interface_size,ninterface))
       allocate(ibool_interfaces_elastic(NGLLX*max_interface_size,ninterface))
       allocate(ibool_interfaces_poroelastic(NGLLX*max_interface_size,ninterface))
       allocate(ibool_interfaces_ext_mesh_init(NGLLX*max_interface_size,ninterface))
       allocate(nibool_interfaces_acoustic(ninterface))
       allocate(nibool_interfaces_elastic(ninterface))
       allocate(nibool_interfaces_poroelastic(ninterface))
       allocate(nibool_interfaces_ext_mesh(ninterface))
       allocate(inum_interfaces_acoustic(ninterface))
       allocate(inum_interfaces_elastic(ninterface))
       allocate(inum_interfaces_poroelastic(ninterface))
   call read_databases_interfaces()

  else
       allocate(my_neighbours(1))
       allocate(my_nelmnts_neighbours(1))
       allocate(my_interfaces(1,1,1))
       allocate(ibool_interfaces_acoustic(1,1))
       allocate(ibool_interfaces_elastic(1,1))
       allocate(ibool_interfaces_poroelastic(1,1))
       allocate(ibool_interfaces_ext_mesh_init(1,1))
       allocate(nibool_interfaces_acoustic(1))
       allocate(nibool_interfaces_elastic(1))
       allocate(nibool_interfaces_poroelastic(1))
       allocate(nibool_interfaces_ext_mesh(1))
       allocate(inum_interfaces_acoustic(1))
       allocate(inum_interfaces_elastic(1))
       allocate(inum_interfaces_poroelastic(1))
  endif


! --- allocate arrays for absorbing boundary conditions

  if(nelemabs <= 0) then
    nelemabs = 1
    anyabs = .false.
  else
    anyabs = .true.
  endif

    allocate(numabs(nelemabs))
    allocate(codeabs(4,nelemabs))

!---codeabs_corner(1,nelemabs) denotes whether element is on bottom-left corner of absorbing boundary or not
!---codeabs_corner(2,nelemabs) denotes whether element is on bottom-right corner of absorbing boundary or not
!---codeabs_corner(3,nelemabs) denotes whether element is on top-left corner of absorbing boundary or not
!---codeabs_corner(4,nelemabs) denotes whether element is on top-right corner of absorbing boundary or not
    allocate(codeabs_corner(4,nelemabs))
    allocate(typeabs(nelemabs))

    allocate(ibegin_edge1(nelemabs))
    allocate(iend_edge1(nelemabs))
    allocate(ibegin_edge3(nelemabs))
    allocate(iend_edge3(nelemabs))

    allocate(ibegin_edge4(nelemabs))
    allocate(iend_edge4(nelemabs))
    allocate(ibegin_edge2(nelemabs))
    allocate(iend_edge2(nelemabs))

    allocate(ibegin_edge1_poro(nelemabs))
    allocate(iend_edge1_poro(nelemabs))
    allocate(ibegin_edge3_poro(nelemabs))
    allocate(iend_edge3_poro(nelemabs))

    allocate(ibegin_edge4_poro(nelemabs))
    allocate(iend_edge4_poro(nelemabs))
    allocate(ibegin_edge2_poro(nelemabs))
    allocate(iend_edge2_poro(nelemabs))

    allocate(ib_left(nelemabs))
    allocate(ib_right(nelemabs))
    allocate(ib_bottom(nelemabs))
    allocate(ib_top(nelemabs))

  !
  !----  read absorbing boundary data
  !
  call read_databases_absorbing()

  if(anyabs .and. (.not. PML_BOUNDARY_CONDITIONS))then
    STACEY_BOUNDARY_CONDITIONS = .true.
  else
    STACEY_BOUNDARY_CONDITIONS = .false.
  endif


  if( anyabs ) then
    ! files to save absorbed waves needed to reconstruct backward wavefield for adjoint method
      if(any_elastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3).and. (.not. PML_BOUNDARY_CONDITIONS)) then
        allocate(b_absorb_elastic_left(3,NGLLZ,nspec_left,NSTEP))
        allocate(b_absorb_elastic_right(3,NGLLZ,nspec_right,NSTEP))
        allocate(b_absorb_elastic_bottom(3,NGLLX,nspec_bottom,NSTEP))
        allocate(b_absorb_elastic_top(3,NGLLX,nspec_top,NSTEP))
      else
        allocate(b_absorb_elastic_left(1,1,1,1))
        allocate(b_absorb_elastic_right(1,1,1,1))
        allocate(b_absorb_elastic_bottom(1,1,1,1))
        allocate(b_absorb_elastic_top(1,1,1,1))
      endif
      if(any_poroelastic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3).and. (.not. PML_BOUNDARY_CONDITIONS)) then
        allocate(b_absorb_poro_s_left(NDIM,NGLLZ,nspec_left,NSTEP))
        allocate(b_absorb_poro_s_right(NDIM,NGLLZ,nspec_right,NSTEP))
        allocate(b_absorb_poro_s_bottom(NDIM,NGLLX,nspec_bottom,NSTEP))
        allocate(b_absorb_poro_s_top(NDIM,NGLLX,nspec_top,NSTEP))
        allocate(b_absorb_poro_w_left(NDIM,NGLLZ,nspec_left,NSTEP))
        allocate(b_absorb_poro_w_right(NDIM,NGLLZ,nspec_right,NSTEP))
        allocate(b_absorb_poro_w_bottom(NDIM,NGLLX,nspec_bottom,NSTEP))
        allocate(b_absorb_poro_w_top(NDIM,NGLLX,nspec_top,NSTEP))
      else
        allocate(b_absorb_poro_s_left(1,1,1,1))
        allocate(b_absorb_poro_s_right(1,1,1,1))
        allocate(b_absorb_poro_s_bottom(1,1,1,1))
        allocate(b_absorb_poro_s_top(1,1,1,1))
        allocate(b_absorb_poro_w_left(1,1,1,1))
        allocate(b_absorb_poro_w_right(1,1,1,1))
        allocate(b_absorb_poro_w_bottom(1,1,1,1))
        allocate(b_absorb_poro_w_top(1,1,1,1))
      endif
      if(any_acoustic .and. (SAVE_FORWARD .or. SIMULATION_TYPE == 3) .and. (.not. PML_BOUNDARY_CONDITIONS)) then
        allocate(b_absorb_acoustic_left(NGLLZ,nspec_left,NSTEP))
        allocate(b_absorb_acoustic_right(NGLLZ,nspec_right,NSTEP))
        allocate(b_absorb_acoustic_bottom(NGLLX,nspec_bottom,NSTEP))
        allocate(b_absorb_acoustic_top(NGLLX,nspec_top,NSTEP))
      else
        allocate(b_absorb_acoustic_left(1,1,1))
        allocate(b_absorb_acoustic_right(1,1,1))
        allocate(b_absorb_acoustic_bottom(1,1,1))
        allocate(b_absorb_acoustic_top(1,1,1))
      endif

  else

    if(.not. allocated(b_absorb_elastic_left)) then
      allocate(b_absorb_elastic_left(1,1,1,1))
      allocate(b_absorb_elastic_right(1,1,1,1))
      allocate(b_absorb_elastic_bottom(1,1,1,1))
      allocate(b_absorb_elastic_top(1,1,1,1))
    endif

    if(.not. allocated(b_absorb_poro_s_left)) then
      allocate(b_absorb_poro_s_left(1,1,1,1))
      allocate(b_absorb_poro_s_right(1,1,1,1))
      allocate(b_absorb_poro_s_bottom(1,1,1,1))
      allocate(b_absorb_poro_s_top(1,1,1,1))
      allocate(b_absorb_poro_w_left(1,1,1,1))
      allocate(b_absorb_poro_w_right(1,1,1,1))
      allocate(b_absorb_poro_w_bottom(1,1,1,1))
      allocate(b_absorb_poro_w_top(1,1,1,1))
    endif

    if(.not. allocated(b_absorb_acoustic_left)) then
      allocate(b_absorb_acoustic_left(1,1,1))
      allocate(b_absorb_acoustic_right(1,1,1))
      allocate(b_absorb_acoustic_bottom(1,1,1))
      allocate(b_absorb_acoustic_top(1,1,1))
    endif

  endif

! --- allocate arrays for acoustic forcing boundary conditions

  if(.not. ACOUSTIC_FORCING) then
    nelem_acforcing = 1
  endif

    allocate(numacforcing(nelem_acforcing))
    allocate(codeacforcing(4,nelem_acforcing))
    allocate(typeacforcing(nelem_acforcing))

    allocate(ibegin_edge1_acforcing(nelem_acforcing))
    allocate(iend_edge1_acforcing(nelem_acforcing))
    allocate(ibegin_edge3_acforcing(nelem_acforcing))
    allocate(iend_edge3_acforcing(nelem_acforcing))

    allocate(ibegin_edge4_acforcing(nelem_acforcing))
    allocate(iend_edge4_acforcing(nelem_acforcing))
    allocate(ibegin_edge2_acforcing(nelem_acforcing))
    allocate(iend_edge2_acforcing(nelem_acforcing))

    allocate(ib_left_acforcing(nelem_acforcing))
    allocate(ib_right_acforcing(nelem_acforcing))
    allocate(ib_bottom_acforcing(nelem_acforcing))
    allocate(ib_top_acforcing(nelem_acforcing))

  !
  !----  read acoustic forcing boundary data
  !
  call read_databases_acoustic_forcing()


!
!----  read acoustic free surface data
!
  if(nelem_acoustic_surface > 0) then
    any_acoustic_edges = .true.
  else
    any_acoustic_edges = .false.
    nelem_acoustic_surface = 1
  endif
  allocate(acoustic_edges(4,nelem_acoustic_surface))
  allocate(acoustic_surface(5,nelem_acoustic_surface))
  call read_databases_free_surf()
  ! resets nelem_acoustic_surface
  if( any_acoustic_edges .eqv. .false. ) nelem_acoustic_surface = 0

  ! constructs acoustic surface
  if(nelem_acoustic_surface > 0) then
    call construct_acoustic_surface ()
    if (myrank == 0) then
      write(IOUT,*)
      write(IOUT,*) 'Number of free surface elements: ',nelem_acoustic_surface
    endif
  endif


  !
  !---- read coupled edges
  !
  if( num_fluid_solid_edges > 0 ) then
    any_fluid_solid_edges = .true.
  else
    any_fluid_solid_edges = .false.
    num_fluid_solid_edges = 1
  endif
  allocate(fluid_solid_acoustic_ispec(num_fluid_solid_edges))
  allocate(fluid_solid_acoustic_iedge(num_fluid_solid_edges))
  allocate(fluid_solid_elastic_ispec(num_fluid_solid_edges))
  allocate(fluid_solid_elastic_iedge(num_fluid_solid_edges))
  if( num_fluid_poro_edges > 0 ) then
    any_fluid_poro_edges = .true.
  else
    any_fluid_poro_edges = .false.
    num_fluid_poro_edges = 1
  endif
  allocate(fluid_poro_acoustic_ispec(num_fluid_poro_edges))
  allocate(fluid_poro_acoustic_iedge(num_fluid_poro_edges))
  allocate(fluid_poro_poroelastic_ispec(num_fluid_poro_edges))
  allocate(fluid_poro_poroelastic_iedge(num_fluid_poro_edges))
  if ( num_solid_poro_edges > 0 ) then
    any_solid_poro_edges = .true.
  else
    any_solid_poro_edges = .false.
    num_solid_poro_edges = 1
  endif
  allocate(solid_poro_elastic_ispec(num_solid_poro_edges))
  allocate(solid_poro_elastic_iedge(num_solid_poro_edges))
  allocate(solid_poro_poroelastic_ispec(num_solid_poro_edges))
  allocate(solid_poro_poroelastic_iedge(num_solid_poro_edges))

  call read_databases_coupled()

  ! resets counters
  if( any_fluid_solid_edges .eqv. .false. ) num_fluid_solid_edges = 0
  if( any_fluid_poro_edges .eqv. .false. ) num_fluid_poro_edges = 0
  if( any_solid_poro_edges .eqv. .false. ) num_solid_poro_edges = 0


  !
  !---- read tangential detection curve
  !      and close Database file
  !
  if (nnodes_tangential_curve > 0) then
    any_tangential_curve = .true.
  else
    any_tangential_curve = .false.
    nnodes_tangential_curve = 1
  endif
  allocate(nodes_tangential_curve(2,nnodes_tangential_curve))
  allocate(dist_tangential_detection_curve(nnodes_tangential_curve))
  call read_tangential_detection_curve()
  ! resets nnode_tangential_curve
  if( any_tangential_curve .eqv. .false. ) nnodes_tangential_curve = 0

!
!----  read axial elements data
!
  allocate(is_on_the_axis(nspec),stat=ier)
  if(ier /= 0) stop 'error: not enough memory to allocate array is_on_the_axis'
  is_on_the_axis(:) = .false.
  if(nelem_on_the_axis == 0) then
    allocate(ispec_of_axial_elements(1))
  else
    allocate(ispec_of_axial_elements(nelem_on_the_axis))
    call read_databases_axial_elements()
    call build_is_on_the_axis()
  endif
  if (myrank == 0 .and. AXISYM) then
    write(IOUT,*)
    write(IOUT,*) 'Number of elements on the axis (for process 0): ',nelem_on_the_axis
  endif

!
!----  end of reading
!

! closes input Database file
 close(IIN)

end subroutine prepare_timerun_read







subroutine prepare_timerun_noise()


#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

  integer i,j,iglob,ispec

#ifdef USE_MPI
  include "precision.h"
#endif


!<NOISE_TOMOGRAPHY

  if (NOISE_TOMOGRAPHY /= 0) then

    !allocate arrays for noise tomography
    allocate(time_function_noise(NSTEP))
    allocate(source_array_noise(3,NGLLX,NGLLZ,NSTEP))
    allocate(mask_noise(nglob))
    allocate(surface_movie_x_noise(nglob))
    allocate(surface_movie_y_noise(nglob))
    allocate(surface_movie_z_noise(nglob))

    !read in parameters for noise tomography
    call read_parameters_noise()

  endif ! NOISE_TOMOGRAPHY /= 0


  if (NOISE_TOMOGRAPHY == 1) then
    call compute_source_array_noise()

    !write out coordinates of mesh
    open(unit=504,file='OUTPUT_FILES/mesh_spec',status='unknown',action='write')
      do ispec = 1, nspec
        do j = 1, NGLLZ
          do i = 1, NGLLX
            iglob = ibool(i,j,ispec)
            write(504,'(1pe11.3,1pe11.3,2i3,i7)') coord(1,iglob), coord(2,iglob), i, j, ispec
         enddo
        enddo
      enddo
    close(504)

    open(unit=504,file='OUTPUT_FILES/mesh_glob',status='unknown',action='write')
      do iglob = 1, nglob
        write(504,'(1pe11.3,1pe11.3,i7)') coord(1,iglob), coord(2,iglob), iglob
      enddo
    close(504)

    !write out spatial distribution of noise sources
    call create_mask_noise()
    open(unit=504,file='OUTPUT_FILES/mask_noise',status='unknown',action='write')
      do iglob = 1, nglob
            write(504,'(1pe11.3,1pe11.3,1pe11.3)') coord(1,iglob), coord(2,iglob), mask_noise(iglob)
      enddo
    close(504)

    !write out velocity model
    if(assign_external_model) then
      open(unit=504,file='OUTPUT_FILES/model_rho_vp_vs',status='unknown',action='write')
        do ispec = 1, nspec
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)
              write(504,'(1pe11.3,1pe11.3,1pe11.3,1pe11.3,1pe11.3)') &
                coord(1,iglob), coord(2,iglob), rhoext(i,j,ispec), vpext(i,j,ispec), vsext(i,j,ispec)
            enddo
          enddo
        enddo
      close(504)
    else
      open(unit=504,file='OUTPUT_FILES/model_rho_kappa_mu',status='unknown',action='write')
        do ispec = 1, nspec
          do j = 1, NGLLZ
            do i = 1, NGLLX
              iglob = ibool(i,j,ispec)
              write(504,'(1pe11.3,1pe11.3,1pe11.3,1pe11.3,1pe11.3)') &
                coord(1,iglob), coord(2,iglob), density(1,kmato(ispec)), &
                poroelastcoef(1,1,kmato(ispec)) + 2.d0/3.d0*poroelastcoef(2,1,kmato(ispec)), &
                poroelastcoef(2,1,kmato(ispec))

            enddo
          enddo
        enddo
      close(504)
    endif

  else if (NOISE_TOMOGRAPHY == 2) then
    call create_mask_noise()

  else if (NOISE_TOMOGRAPHY == 3) then

    if (output_wavefields_noise) then
      call create_mask_noise()

      !prepare array that will hold wavefield snapshots
      noise_output_ncol = 5
      allocate(noise_output_array(noise_output_ncol,nglob))
      allocate(noise_output_rhokl(nglob))
    endif

  endif

!>NOISE_TOMOGRAPHY

end subroutine prepare_timerun_noise








subroutine prepare_timerun_attenuation()


#ifdef USE_MPI
  use mpi
#endif

  use specfem_par

  implicit none

#ifdef USE_MPI
  include "precision.h"
#endif


! Precompute Runge Kutta coefficients if viscous attenuation
  if(ATTENUATION_PORO_FLUID_PART) then
! viscous attenuation is implemented following the memory variable formulation of
! J. M. Carcione Wave fields in real media: wave propagation in anisotropic,
! anelastic and porous media, Elsevier, p. 304-305, 2007
    theta_e = (sqrt(Q0**2+1.d0) +1.d0)/(2.d0*pi*freq0*Q0)
    theta_s = (sqrt(Q0**2+1.d0) -1.d0)/(2.d0*pi*freq0*Q0)

    thetainv = - 1.d0 / theta_s
    alphaval = 1.d0 + deltat*thetainv + deltat**2*thetainv**2 / 2.d0 + &
      deltat**3*thetainv**3 / 6.d0 + deltat**4*thetainv**4 / 24.d0
    betaval = deltat / 2.d0 + deltat**2*thetainv / 3.d0 + deltat**3*thetainv**2 / 8.d0 + deltat**4*thetainv**3 / 24.d0
    gammaval = deltat / 2.d0 + deltat**2*thetainv / 6.d0 + deltat**3*thetainv**2 / 24.d0
   print*,'************************************************************'
   print*,'****** Visco attenuation coefficients (poroelastic) ********'
   print*,'theta_e = ', theta_e
   print*,'theta_s = ', theta_s
   print*,'alpha = ', alphaval
   print*,'beta = ', betaval
   print*,'gamma = ', gammaval
   print*,'************************************************************'

! initialize memory variables for attenuation
    viscox(:,:,:) = 0.d0
    viscoz(:,:,:) = 0.d0
    rx_viscous(:,:,:) = 0.d0
    rz_viscous(:,:,:) = 0.d0
    if(time_stepping_scheme == 2) then
     rx_viscous_LDDRK = 0.d0
     rz_viscous_LDDRK = 0.d0
    endif

    if(time_stepping_scheme == 3) then
     rx_viscous_initial_rk = 0.d0
     rz_viscous_initial_rk = 0.d0
     rx_viscous_force_RK = 0.d0
     rz_viscous_force_RK = 0.d0
    endif

  endif

! allocate arrays for postscript output
#ifdef USE_MPI
  if(modelvect) then
  d1_coorg_recv_ps_velocity_model=2
  call mpi_allreduce(nspec,d2_coorg_recv_ps_velocity_model,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  d2_coorg_recv_ps_velocity_model=d2_coorg_recv_ps_velocity_model*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
       ((NGLLX-subsamp_postscript)/subsamp_postscript)*4
  d1_RGB_recv_ps_velocity_model=1
  call mpi_allreduce(nspec,d2_RGB_recv_ps_velocity_model,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  d2_RGB_recv_ps_velocity_model=d2_RGB_recv_ps_velocity_model*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
       ((NGLLX-subsamp_postscript)/subsamp_postscript)*4
  else
  d1_coorg_recv_ps_velocity_model=1
  d2_coorg_recv_ps_velocity_model=1
  d1_RGB_recv_ps_velocity_model=1
  d2_RGB_recv_ps_velocity_model=1
  endif

  d1_coorg_send_ps_element_mesh=2
  if ( ngnod == 4 ) then
    if ( numbers == 1 ) then
      d2_coorg_send_ps_element_mesh=nspec*5
      if ( colors == 1 ) then
        d1_color_send_ps_element_mesh=2*nspec
      else
        d1_color_send_ps_element_mesh=1*nspec
      endif
    else
      d2_coorg_send_ps_element_mesh=nspec*6
      if ( colors == 1 ) then
        d1_color_send_ps_element_mesh=1*nspec
      endif
    endif
  else
    if ( numbers == 1 ) then
      d2_coorg_send_ps_element_mesh=nspec*((pointsdisp-1)*3+max(0,pointsdisp-2)+1+1)
      if ( colors == 1 ) then
        d1_color_send_ps_element_mesh=2*nspec
      else
        d1_color_send_ps_element_mesh=1*nspec
      endif
    else
      d2_coorg_send_ps_element_mesh=nspec*((pointsdisp-1)*3+max(0,pointsdisp-2)+1)
      if ( colors == 1 ) then
        d1_color_send_ps_element_mesh=1*nspec
      endif
    endif
  endif

  call mpi_allreduce(d1_coorg_send_ps_element_mesh,d1_coorg_recv_ps_element_mesh,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  call mpi_allreduce(d2_coorg_send_ps_element_mesh,d2_coorg_recv_ps_element_mesh,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  call mpi_allreduce(d1_color_send_ps_element_mesh,d1_color_recv_ps_element_mesh,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  d1_coorg_send_ps_abs=4
  d2_coorg_send_ps_abs=4*nelemabs
  call mpi_allreduce(d1_coorg_send_ps_abs,d1_coorg_recv_ps_abs,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  call mpi_allreduce(d2_coorg_send_ps_abs,d2_coorg_recv_ps_abs,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  d1_coorg_send_ps_free_surface=4
  d2_coorg_send_ps_free_surface=4*nelem_acoustic_surface
  call mpi_allreduce(d1_coorg_send_ps_free_surface,d1_coorg_recv_ps_free_surface,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  call mpi_allreduce(d2_coorg_send_ps_free_surface,d2_coorg_recv_ps_free_surface,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)

  d1_coorg_send_ps_vector_field=8
  if(interpol) then
    if(plot_lowerleft_corner_only) then
      d2_coorg_send_ps_vector_field=nspec*1*1
    else
      d2_coorg_send_ps_vector_field=nspec*pointsdisp*pointsdisp
    endif
  else
    d2_coorg_send_ps_vector_field=nglob
  endif
  call mpi_allreduce(d1_coorg_send_ps_vector_field,d1_coorg_recv_ps_vector_field,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)
  call mpi_allreduce(d2_coorg_send_ps_vector_field,d2_coorg_recv_ps_vector_field,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ier)


#else
  d1_coorg_recv_ps_velocity_model=1
  d2_coorg_recv_ps_velocity_model=1
  d1_RGB_recv_ps_velocity_model=1
  d2_RGB_recv_ps_velocity_model=1

  d1_coorg_send_ps_element_mesh=1
  d2_coorg_send_ps_element_mesh=1
  d1_coorg_recv_ps_element_mesh=1
  d2_coorg_recv_ps_element_mesh=1
  d1_color_send_ps_element_mesh=1
  d1_color_recv_ps_element_mesh=1

  d1_coorg_send_ps_abs=1
  d2_coorg_send_ps_abs=1
  d1_coorg_recv_ps_abs=1
  d2_coorg_recv_ps_abs=1
  d1_coorg_send_ps_free_surface=1
  d2_coorg_send_ps_free_surface=1
  d1_coorg_recv_ps_free_surface=1
  d2_coorg_recv_ps_free_surface=1

  d1_coorg_send_ps_vector_field=1
  d2_coorg_send_ps_vector_field=1
  d1_coorg_recv_ps_vector_field=1
  d2_coorg_recv_ps_vector_field=1

#endif
  d1_coorg_send_ps_velocity_model=2
  d2_coorg_send_ps_velocity_model=nspec*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
                                        ((NGLLX-subsamp_postscript)/subsamp_postscript)*4
  d1_RGB_send_ps_velocity_model=1
  d2_RGB_send_ps_velocity_model=nspec*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
                                      ((NGLLX-subsamp_postscript)/subsamp_postscript)

  allocate(coorg_send_ps_velocity_model(d1_coorg_send_ps_velocity_model,d2_coorg_send_ps_velocity_model))
  allocate(RGB_send_ps_velocity_model(d1_RGB_send_ps_velocity_model,d2_RGB_send_ps_velocity_model))

  allocate(coorg_recv_ps_velocity_model(d1_coorg_recv_ps_velocity_model,d2_coorg_recv_ps_velocity_model))
  allocate(RGB_recv_ps_velocity_model(d1_RGB_recv_ps_velocity_model,d2_RGB_recv_ps_velocity_model))

  allocate(coorg_send_ps_element_mesh(d1_coorg_send_ps_element_mesh,d2_coorg_send_ps_element_mesh))
  allocate(coorg_recv_ps_element_mesh(d1_coorg_recv_ps_element_mesh,d2_coorg_recv_ps_element_mesh))
  allocate(color_send_ps_element_mesh(d1_color_send_ps_element_mesh))
  allocate(color_recv_ps_element_mesh(d1_color_recv_ps_element_mesh))

  allocate(coorg_send_ps_abs(d1_coorg_send_ps_abs,d2_coorg_send_ps_abs))
  allocate(coorg_recv_ps_abs(d1_coorg_recv_ps_abs,d2_coorg_recv_ps_abs))

  allocate(coorg_send_ps_free_surface(d1_coorg_send_ps_free_surface,d2_coorg_send_ps_free_surface))
  allocate(coorg_recv_ps_free_surface(d1_coorg_recv_ps_free_surface,d2_coorg_recv_ps_free_surface))

  allocate(coorg_send_ps_vector_field(d1_coorg_send_ps_vector_field,d2_coorg_send_ps_vector_field))
  allocate(coorg_recv_ps_vector_field(d1_coorg_recv_ps_vector_field,d2_coorg_recv_ps_vector_field))

! to dump the wave field
  this_is_the_first_time_we_dump = .true.

end subroutine prepare_timerun_attenuation

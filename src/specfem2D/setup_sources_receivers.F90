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


  subroutine setup_sources_receivers()

  use constants,only: NGLLX,NGLLZ,NDIM,IMAIN

  use specfem_par, only: NSOURCES,initialfield,source_type, &
                         coord,ibool,nglob,nspec,nelem_acoustic_surface,acoustic_surface, &
                         ispec_is_elastic,ispec_is_poroelastic, &
                         x_source,z_source,ispec_selected_source,ispec_selected_rec, &
                         is_proc_source,nb_proc_source, &
                         sourcearray,Mxx,Mzz,Mxz,xix,xiz,gammax,gammaz,xigll,zigll,npgeo, &
                         nproc,myrank,xi_source,gamma_source,coorg,knods,ngnod, &
                         nrec,nrecloc,recloc,which_proc_receiver,st_xval,st_zval, &
                         xi_receiver,gamma_receiver,station_name,network_name,x_final_receiver,&
                         z_final_receiver,iglob_source
  implicit none

  ! Local variables
  integer ispec_acoustic_surface
  integer  :: ixmin, ixmax, izmin, izmax,i_source,ispec
#ifndef USE_MPI
  integer irec
#endif

  do i_source= 1,NSOURCES

    if (source_type(i_source) == 1) then

      ! collocated force source
      call locate_source_force(ibool,coord,nspec,nglob,xigll,zigll,x_source(i_source),z_source(i_source), &
          ispec_selected_source(i_source),is_proc_source(i_source),nb_proc_source(i_source), &
          nproc,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo, &
          iglob_source(i_source))

      ! check that acoustic source is not exactly on the free surface because pressure is zero there
      if (is_proc_source(i_source) == 1) then
        do ispec_acoustic_surface = 1,nelem_acoustic_surface
          ispec = acoustic_surface(1,ispec_acoustic_surface)
          ixmin = acoustic_surface(2,ispec_acoustic_surface)
          ixmax = acoustic_surface(3,ispec_acoustic_surface)
          izmin = acoustic_surface(4,ispec_acoustic_surface)
          izmax = acoustic_surface(5,ispec_acoustic_surface)
          if (.not. ispec_is_elastic(ispec) .and. .not. ispec_is_poroelastic(ispec) .and. &
            ispec == ispec_selected_source(i_source)) then
            if ((izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==NGLLX .and. &
                gamma_source(i_source) < -0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==NGLLX .and. &
                gamma_source(i_source) > 0.99d0) .or.&
                (izmin==1 .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
                xi_source(i_source) < -0.99d0) .or.&
                (izmin==1 .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                xi_source(i_source) > 0.99d0) .or.&
                (izmin==1 .and. izmax==1 .and. ixmin==1 .and. ixmax==1 .and. &
                gamma_source(i_source) < -0.99d0 .and. xi_source(i_source) < -0.99d0) .or.&
                (izmin==1 .and. izmax==1 .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                gamma_source(i_source) < -0.99d0 .and. xi_source(i_source) > 0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==1 .and. ixmax==1 .and. &
                gamma_source(i_source) > 0.99d0 .and. xi_source(i_source) < -0.99d0) .or.&
                (izmin==NGLLZ .and. izmax==NGLLZ .and. ixmin==NGLLX .and. ixmax==NGLLX .and. &
                gamma_source(i_source) > 0.99d0 .and. xi_source(i_source) > 0.99d0)) then
              call exit_MPI(myrank,'an acoustic source cannot be located exactly '// &
                            'on the free surface because pressure is zero there')
            endif
          endif
        enddo
      endif

    else if (source_type(i_source) == 2) then
      ! moment-tensor source
      call locate_source_moment_tensor(ibool,coord,nspec,nglob,xigll,zigll,x_source(i_source),z_source(i_source), &
             ispec_selected_source(i_source),is_proc_source(i_source),nb_proc_source(i_source), &
             nproc,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo)

      ! compute source array for moment-tensor source
      call compute_arrays_source(ispec_selected_source(i_source),xi_source(i_source),gamma_source(i_source),&
             sourcearray(i_source,1,1,1), &
             Mxx(i_source),Mzz(i_source),Mxz(i_source),xix,xiz,gammax,gammaz,xigll,zigll,nspec)

    else if (.not.initialfield) then

      call exit_MPI(myrank,'incorrect source type')

    endif

  enddo ! do i_source= 1,NSOURCES

  ! locate receivers in the mesh
  call locate_receivers(ibool,coord,nspec,nglob,xigll,zigll, &
                      nrec,nrecloc,recloc,which_proc_receiver,nproc,myrank, &
                      st_xval,st_zval,ispec_selected_rec, &
                      xi_receiver,gamma_receiver,station_name,network_name, &
                      x_source(1),z_source(1), &
                      coorg,knods,ngnod,npgeo, &
                      x_final_receiver,z_final_receiver)

!! DK DK this below not supported in the case of MPI yet, we should do a MPI_GATHER() of the values
!! DK DK and use "if (myrank == which_proc_receiver(irec)) then" to display the right sources
!! DK DK and receivers carried by each mesh slice, and not fictitious values coming from other slices
#ifndef USE_MPI
  if (myrank == 0) then

     ! write actual source locations to file
     ! note that these may differ from input values, especially if source_surf = .true. in SOURCE
     ! note that the exact source locations are determined from (ispec,xi,gamma) values
     open(unit=14,file='OUTPUT_FILES/for_information_SOURCE_actually_used',status='unknown')
     do i_source= 1,NSOURCES
        write(14,*) x_source(i_source), z_source(i_source)
     enddo
     close(14)

     ! write out actual station locations (compare with STATIONS from meshfem2D)
     ! NOTE: this will be written out even if use_existing_STATIONS = .true.
     open(unit=15,file='OUTPUT_FILES/for_information_STATIONS_actually_used',status='unknown')
     do irec = 1,nrec
        write(15,"('S',i4.4,'    AA ',f20.7,1x,f20.7,'       0.0         0.0')") &
             irec,x_final_receiver(irec),z_final_receiver(irec)
     enddo
     close(15)

  endif
#endif

  end subroutine setup_sources_receivers



! =====

  subroutine add_adjoint_sources_SU()

  use specfem_par, only: AXISYM,xiglj,is_on_the_axis,myrank, NSTEP, nrec, xi_receiver, gamma_receiver, which_proc_receiver, &
                         xigll,zigll,hxir,hgammar,hpxir,hpgammar, &
                         adj_sourcearray, adj_sourcearrays, &
                         r4head, header2, source_adjointe, GPU_MODE, ispec_selected_rec

  use constants,only: CUSTOM_REAL,MAX_STRING_LEN,NGLLX,NGLLZ,NGLJ
  implicit none

  ! local parameters
  integer :: i, k, irec, irec_local, ier
  real(kind=CUSTOM_REAL), allocatable, dimension(:,:) :: adj_src_s
  character(len=MAX_STRING_LEN) :: filename

  irec_local = 0
  write(filename, "('./SEM/Ux_file_single.su.adj')")
  open(111,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
  if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

  write(filename, "('./SEM/Uy_file_single.su.adj')")
  open(112,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
  if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

  write(filename, "('./SEM/Uz_file_single.su.adj')")
  open(113,file=trim(filename),access='direct',recl=240+4*NSTEP,iostat = ier)
  if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' does not exist')

  allocate(adj_src_s(NSTEP,3))
  adj_src_s(:,:) = 0.

  do irec = 1, nrec
    if (myrank == which_proc_receiver(irec)) then
      irec_local = irec_local + 1
      adj_sourcearray(:,:,:,:) = 0.0
      read(111,rec=irec,iostat=ier) r4head, adj_src_s(:,1)
      if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')

      read(112,rec=irec,iostat=ier) r4head, adj_src_s(:,2)
      if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')

      read(113,rec=irec,iostat=ier) r4head, adj_src_s(:,3)
      if (ier /= 0) call exit_MPI(myrank,'file '//trim(filename)//' read error')

      header2=int(r4head(29), kind=2)
      if (irec==1) print *, r4head(1),r4head(19),r4head(20),r4head(21),r4head(22),header2(2)

      if (AXISYM) then
        if (is_on_the_axis(ispec_selected_rec(irec))) then ! TODO verify ispec_selected_rec and not ispec_selected_source
          call lagrange_any(xi_receiver(irec),NGLJ,xiglj,hxir,hpxir)
        else
          call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
        endif
      else
        call lagrange_any(xi_receiver(irec),NGLLX,xigll,hxir,hpxir)
      endif

      call lagrange_any(gamma_receiver(irec),NGLLZ,zigll,hgammar,hpgammar)
      source_adjointe(irec_local,:,1) = adj_src_s(:,1)
      source_adjointe(irec_local,:,2) = adj_src_s(:,3)

      if (.not. GPU_MODE) then
        do k = 1, NGLLZ
          do i = 1, NGLLX
            adj_sourcearray(:,:,i,k) = hxir(i) * hgammar(k) * adj_src_s(:,:)
          enddo
        enddo
        adj_sourcearrays(irec_local,:,:,:,:) = adj_sourcearray(:,:,:,:)
      endif

    endif !  if (myrank == which_proc_receiver(irec))
  enddo ! irec
  close(111)
  close(112)
  close(113)
  deallocate(adj_src_s)

  end subroutine add_adjoint_sources_SU


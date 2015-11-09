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



subroutine init_host_to_dev_variable()

  use specfem_par
  implicit none

  integer :: i_spec_free, ipoint1D, i, j, k, ispec, ispecabs, i_source, ispec_inner, ispec_outer

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Initialisation variables pour routine prepare_constants_device
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



deltatf=sngl(deltat)
deltatover2f=sngl(deltatover2)
deltatsquareover2f=sngl(deltatsquareover2)
b_deltatf=sngl(b_deltat)
b_deltatover2f=sngl(b_deltatover2)
b_deltatsquareover2f=sngl(b_deltatsquareover2)




NSPEC_AB = nspec
NGLOB_AB = nglob




allocate(ibool_interfaces_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))

ibool_interfaces_ext_mesh(:,:)=0
do j=1,ninterface
do i=1,nibool_interfaces_ext_mesh(j)
ibool_interfaces_ext_mesh(i,j)=ibool_interfaces_ext_mesh_init(i,j)
enddo
enddo






!!!
allocate(free_ac_ispec(nelem_acoustic_surface))

free_ac_ispec(:)=acoustic_surface(1,:)



  ! checks
  if( nelemabs < 0 ) then
    print *,'host_to_dev: reading in negative nelemabs ',nelemabs,'...resetting to zero'
    nelemabs = 0
  endif

  allocate(abs_boundary_ij(2,NGLLX,nelemabs), &
           abs_boundary_jacobian1Dw(NGLLX,nelemabs), &
           abs_boundary_normal(NDIM,NGLLX,nelemabs),&
           cote_abs(nelemabs),stat=ier)
  if( ier /= 0 ) stop 'error allocating array abs_boundary_ispec etc.'



  if(STACEY_BOUNDARY_CONDITIONS) then

    do ispecabs=1,nelemabs
      ispec = numabs(ispecabs)
        !--- left absorbing boundary
        if(codeabs(IEDGE4,ispecabs)) then
        i = 1
        do j = 1,NGLLZ

            abs_boundary_ij(1,j,ispecabs) = i
            abs_boundary_ij(2,j,ispecabs) = j

            xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)

            abs_boundary_normal(1,j,ispecabs) = - zgamma / jacobian1D
            abs_boundary_normal(2,j,ispecabs) = + xgamma / jacobian1D

            abs_boundary_jacobian1Dw(j,ispecabs) = jacobian1D * wzgll(j)

            cote_abs(ispecabs) = 4

          enddo

        !--- right absorbing boundary
        else if(codeabs(IEDGE2,ispecabs)) then
          i = NGLLX
           do j = 1,NGLLZ

            abs_boundary_ij(1,j,ispecabs) = i
            abs_boundary_ij(2,j,ispecabs) = j

            xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
            zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)

            abs_boundary_normal(1,j,ispecabs) = + zgamma / jacobian1D
            abs_boundary_normal(2,j,ispecabs) = - xgamma / jacobian1D

            abs_boundary_jacobian1Dw(j,ispecabs) = jacobian1D * wzgll(j)

            cote_abs(ispecabs) = 2

          enddo

           !--- bottom absorbing boundary
      else if(codeabs(IEDGE1,ispecabs)) then
        j = 1
         do i = 1,NGLLX

            abs_boundary_ij(1,i,ispecabs) = i
            abs_boundary_ij(2,i,ispecabs) = j

              xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
              zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
              jacobian1D = sqrt(xxi**2 + zxi**2)

              abs_boundary_normal(1,i,ispecabs) = + zxi / jacobian1D
              abs_boundary_normal(2,i,ispecabs) = - xxi / jacobian1D

              abs_boundary_jacobian1Dw(i,ispecabs) = jacobian1D * wxgll(i)

              cote_abs(ispecabs) = 1

          enddo

              !--- top absorbing boundary
      else if(codeabs(IEDGE3,ispecabs)) then
        j = NGLLZ
        do i = 1,NGLLX

            abs_boundary_ij(1,i,ispecabs) = i
            abs_boundary_ij(2,i,ispecabs) = j

              xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
              zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
              jacobian1D = sqrt(xxi**2 + zxi**2)

              abs_boundary_normal(1,i,ispecabs) = - zxi / jacobian1D
              abs_boundary_normal(2,i,ispecabs) = + xxi / jacobian1D

              abs_boundary_jacobian1Dw(i,ispecabs) = jacobian1D * wxgll(i)

              cote_abs(ispecabs) = 3

         enddo

        endif
      enddo
   endif



   nsources_local=0

   do i = 1, NSOURCES

      if (is_proc_source(i) == 1) then
      nsources_local = nsources_local +1
      endif
   enddo


  allocate(source_time_function_loc(nsources_local,NSTEP))
  allocate(ispec_selected_source_loc(nsources_local))

  j=0
  do i = 1, NSOURCES
      if (is_proc_source(i) == 1) then
        if (j>nsources_local) stop 'error with the number of local sources'
        j=j+1
        source_time_function_loc(j,:) = source_time_function(i,:,1)
        ispec_selected_source_loc(j)  = ispec_selected_source(i)
      endif
  enddo


  if ( nsources_local > 0 ) then
    allocate(sourcearray_loc(nsources_local,NDIM,NGLLX,NGLLX))
  else
    allocate(sourcearray_loc(1,1,1,1))
  endif

  k=0
  do i_source=1,NSOURCES

    if (is_proc_source(i_source) == 1) then

    k = k + 1

    if(source_type(i_source) == 1) then

      if( acoustic(ispec_selected_source(i_source)) ) then

        do j = 1,NGLLZ
                    do i = 1,NGLLX

        sourcearray_loc(k,1,i,j) = sngl(hxis_store(i_source,i) * hgammas_store(i_source,j))

                    enddo
        enddo

      else if ( elastic(ispec_selected_source(i_source)) ) then
        do j = 1,NGLLZ
                    do i = 1,NGLLX

        sourcearray_loc(k,1,i,j) = - sngl(sin(anglesource(i_source)) * hxis_store(i_source,i) * hgammas_store(i_source,j))
        sourcearray_loc(k,2,i,j) = sngl(cos(anglesource(i_source)) * hxis_store(i_source,i) * hgammas_store(i_source,j))

                    enddo
        enddo
      endif ! is elastic

     else

         sourcearray_loc(k,:,:,:) = sourcearray(i_source,:,:,:)
         sourcearray_loc(k,:,:,:) = sourcearray(i_source,:,:,:)

     endif ! Source_type

     endif ! is_proc_source

  enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!! Init pour prepare acoustique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






  ! sets up elements for loops in acoustic simulations
  nspec_inner_acoustic = 0
  nspec_outer_acoustic = 0
  if( any_acoustic ) then
    ! counts inner and outer elements
    do ispec = 1, nspec
      if(acoustic(ispec) ) then
        if( ispec_is_inner(ispec) .eqv. .true. ) then
          nspec_inner_acoustic = nspec_inner_acoustic + 1
        else
          nspec_outer_acoustic = nspec_outer_acoustic + 1
        endif
      endif
    enddo

    ! stores indices of inner and outer elements for faster(?) computation
    num_phase_ispec_acoustic = max(nspec_inner_acoustic,nspec_outer_acoustic)
    if( num_phase_ispec_acoustic < 0 ) stop 'error acoustic simulation: num_phase_ispec_acoustic is < zero'

    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if( ier /= 0 ) stop 'error allocating array phase_ispec_inner_acoustic'
    phase_ispec_inner_acoustic(:,:) = 0

    ispec_inner = 0
    ispec_outer = 0
    do ispec = 1, nspec
      if(acoustic(ispec) ) then
        if( ispec_is_inner(ispec) .eqv. .true. ) then
          ispec_inner = ispec_inner + 1
          phase_ispec_inner_acoustic(ispec_inner,2) = ispec
        else
          ispec_outer = ispec_outer + 1
          phase_ispec_inner_acoustic(ispec_outer,1) = ispec
        endif
      endif
    enddo
  else
    ! allocates dummy array
    num_phase_ispec_acoustic = 0
    allocate( phase_ispec_inner_acoustic(num_phase_ispec_acoustic,2),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array phase_ispec_inner_acoustic'
    phase_ispec_inner_acoustic(:,:) = 0
  endif

!

  allocate(free_surface_ij(2,NGLLX,nelem_acoustic_surface))

  do i_spec_free = 1, nelem_acoustic_surface


if (acoustic_surface(2,i_spec_free) ==acoustic_surface(3,i_spec_free)) then

     do j=1,5
      free_surface_ij(1,j,i_spec_free) = acoustic_surface(2,i_spec_free)
     enddo


else

      j=1

      do i = acoustic_surface(2,i_spec_free), acoustic_surface(3,i_spec_free)

      free_surface_ij(1,j,i_spec_free) = i

      j=j+1

      enddo
endif


if (acoustic_surface(4,i_spec_free) ==acoustic_surface(5,i_spec_free)) then

     do j=1,5
      free_surface_ij(2,j,i_spec_free) = acoustic_surface(4,i_spec_free)
     enddo


else

    j=1

      do i = acoustic_surface(4,i_spec_free), acoustic_surface(5,i_spec_free)

      free_surface_ij(2,j,i_spec_free) = i

      j=j+1

      enddo
endif

  enddo


!

 allocate(coupling_ac_el_ispec(num_fluid_solid_edges))
 allocate(coupling_ac_el_ij(2,NGLLX,num_fluid_solid_edges))
 allocate(coupling_ac_el_normal(2,NGLLX,num_fluid_solid_edges))
 allocate(coupling_ac_el_jacobian1Dw(NGLLX,num_fluid_solid_edges))
      do inum = 1,num_fluid_solid_edges

        ! get the edge of the acoustic element
        ispec_acoustic = fluid_solid_acoustic_ispec(inum)
        iedge_acoustic = fluid_solid_acoustic_iedge(inum)
        coupling_ac_el_ispec(inum) = ispec_acoustic

        ! get the corresponding edge of the elastic element
        ispec_elastic = fluid_solid_elastic_ispec(inum)
        iedge_elastic = fluid_solid_elastic_iedge(inum)

        ! implement 1D coupling along the edge
        do ipoint1D = 1,NGLLX

          ! get point values for the elastic side, which matches our side in the inverse direction
        coupling_ac_el_ij(1,ipoint1D,inum) = ivalue(ipoint1D,iedge_acoustic)
        coupling_ac_el_ij(2,ipoint1D,inum) = jvalue(ipoint1D,iedge_acoustic)

        i = ivalue(ipoint1D,iedge_acoustic)
        j = jvalue(ipoint1D,iedge_acoustic)




          if(iedge_acoustic == ITOP)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            coupling_ac_el_normal(1,ipoint1D,inum) = - zxi / jacobian1D
            coupling_ac_el_normal(2,ipoint1D,inum) = + xxi / jacobian1D
            coupling_ac_el_jacobian1Dw(ipoint1D,inum) = jacobian1D * wxgll(i)

          else if(iedge_acoustic == IBOTTOM)then
            xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xxi**2 + zxi**2)
            coupling_ac_el_normal(1,ipoint1D,inum) = + zxi / jacobian1D
            coupling_ac_el_normal(2,ipoint1D,inum) = - xxi / jacobian1D
            coupling_ac_el_jacobian1Dw(ipoint1D,inum) = jacobian1D * wxgll(i)

          else if(iedge_acoustic ==ILEFT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            coupling_ac_el_normal(1,ipoint1D,inum) = - zgamma / jacobian1D
            coupling_ac_el_normal(2,ipoint1D,inum) = + xgamma / jacobian1D
            coupling_ac_el_jacobian1Dw(ipoint1D,inum) = jacobian1D * wzgll(j)

          else if(iedge_acoustic ==IRIGHT)then
            xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
            jacobian1D = sqrt(xgamma**2 + zgamma**2)
            coupling_ac_el_normal(1,ipoint1D,inum) = + zgamma / jacobian1D
            coupling_ac_el_normal(2,ipoint1D,inum) = - xgamma / jacobian1D
            coupling_ac_el_jacobian1Dw(ipoint1D,inum) = jacobian1D * wzgll(j)
          endif


        enddo

      enddo


!!


num_colors_outer_acoustic = 0
num_colors_inner_acoustic = 0
allocate(num_elem_colors_acoustic(1))
num_elem_colors_acoustic(1)=0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Initialisation parametres pour simulation elastique
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  ! sets up elements for loops in acoustic simulations
  nspec_inner_elastic = 0
  nspec_outer_elastic = 0
  if( any_elastic ) then
    ! counts inner and outer elements
    do ispec = 1, nspec
      if(elastic(ispec) ) then
        if( ispec_is_inner(ispec) .eqv. .true. ) then
          nspec_inner_elastic = nspec_inner_elastic + 1
        else
          nspec_outer_elastic = nspec_outer_elastic + 1
        endif
      endif
    enddo

    ! stores indices of inner and outer elements for faster(?) computation
    num_phase_ispec_elastic = max(nspec_inner_elastic,nspec_outer_acoustic)
    if( num_phase_ispec_elastic < 0 ) stop 'error elastic simulation: num_phase_ispec_elastic is < zero'

    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if( ier /= 0 ) stop 'error allocating array phase_ispec_inner_elastic'
    phase_ispec_inner_elastic(:,:) = 0

    ispec_inner = 0
    ispec_outer = 0
    do ispec = 1, nspec
      if( elastic(ispec) ) then
        if( ispec_is_inner(ispec) .eqv. .true. ) then
          ispec_inner = ispec_inner + 1
          phase_ispec_inner_elastic(ispec_inner,2) = ispec
        else
          ispec_outer = ispec_outer + 1
          phase_ispec_inner_elastic(ispec_outer,1) = ispec
        endif
      endif
    enddo
  else
    ! allocates dummy array
    num_phase_ispec_elastic = 0
    allocate( phase_ispec_inner_elastic(num_phase_ispec_elastic,2),stat=ier)
    if( ier /= 0 ) stop 'error allocating dummy array phase_ispec_inner_elastic'
    phase_ispec_inner_elastic(:,:) = 0
  endif

!


num_colors_outer_elastic = 0
num_colors_inner_elastic = 0
allocate(num_elem_colors_elastic(1))
num_elem_colors_elastic(1)=0

!

ANY_ANISOTROPY= .false.

    do ispec = 1, nspec
        if (anisotropic(ispec) ) ANY_ANISOTROPY=.true.
    enddo

if (ANY_ANISOTROPY) then

allocate(c11store(NGLLX,NGLLZ,NSPEC))
allocate(c13store(NGLLX,NGLLZ,NSPEC))
allocate(c15store(NGLLX,NGLLZ,NSPEC))
allocate(c33store(NGLLX,NGLLZ,NSPEC))
allocate(c35store(NGLLX,NGLLZ,NSPEC))
allocate(c55store(NGLLX,NGLLZ,NSPEC))
allocate(c12store(NGLLX,NGLLZ,NSPEC))
allocate(c23store(NGLLX,NGLLZ,NSPEC))
allocate(c25store(NGLLX,NGLLZ,NSPEC))

if (assign_external_model) then

do ispec=1,nspec
        do j = 1,NGLLZ
             do i = 1,NGLLX

                c11store(i,j,ispec) = c11ext(i,j,ispec)
                c13store(i,j,ispec) = c13ext(i,j,ispec)
                c15store(i,j,ispec) = c15ext(i,j,ispec)
                c33store(i,j,ispec) = c33ext(i,j,ispec)
                c35store(i,j,ispec) = c35ext(i,j,ispec)
                c55store(i,j,ispec) = c55ext(i,j,ispec)
                c12store(i,j,ispec) = c12ext(i,j,ispec)
                c23store(i,j,ispec) = c23ext(i,j,ispec)
                c25store(i,j,ispec) = c25ext(i,j,ispec)

               enddo
       enddo
enddo

else

do ispec=1,nspec
        do j = 1,NGLLZ
             do i = 1,NGLLX

                    c11store(i,j,ispec) = sngl(anisotropy(1,kmato(ispec)))
                    c13store(i,j,ispec) = sngl(anisotropy(2,kmato(ispec)))
                    c15store(i,j,ispec) = sngl(anisotropy(3,kmato(ispec)))
                    c33store(i,j,ispec) = sngl(anisotropy(4,kmato(ispec)))
                    c35store(i,j,ispec) = sngl(anisotropy(5,kmato(ispec)))
                    c55store(i,j,ispec) = sngl(anisotropy(6,kmato(ispec)))
                    c12store(i,j,ispec) = sngl(anisotropy(7,kmato(ispec)))
                    c23store(i,j,ispec) = sngl(anisotropy(8,kmato(ispec)))
                    c25store(i,j,ispec) = sngl(anisotropy(9,kmato(ispec)))

               enddo
       enddo
enddo

endif

else

allocate(c11store(1,1,1))
allocate(c13store(1,1,1))
allocate(c15store(1,1,1))
allocate(c33store(1,1,1))
allocate(c35store(1,1,1))
allocate(c55store(1,1,1))
allocate(c12store(1,1,1))
allocate(c23store(1,1,1))
allocate(c25store(1,1,1))

endif



!!!!


allocate(displ_2D(2,nglob_elastic))
allocate(veloc_2D(2,nglob_elastic))
allocate(accel_2D(2,nglob_elastic))
displ_2D(1,:)=displ_elastic(1,:)
displ_2D(2,:)=displ_elastic(3,:)

veloc_2D(1,:)=veloc_elastic(1,:)
veloc_2D(2,:)=veloc_elastic(3,:)

accel_2D(1,:)=accel_elastic(1,:)
accel_2D(2,:)=accel_elastic(3,:)

if(SIMULATION_TYPE == 3 .and. any_elastic) then
  allocate(b_displ_2D(2,nglob))
  allocate(b_veloc_2D(2,nglob))
  allocate(b_accel_2D(2,nglob))


b_displ_2D(1,:)=b_displ_elastic(1,:)
b_displ_2D(2,:)=b_displ_elastic(3,:)

b_veloc_2D(1,:)=b_veloc_elastic(1,:)
b_veloc_2D(2,:)=b_veloc_elastic(3,:)

b_accel_2D(1,:)=b_accel_elastic(1,:)
b_accel_2D(2,:)=b_accel_elastic(3,:)

endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate(tab_requests_send_recv_scalar(2*ninterface))
allocate(b_tab_requests_send_recv_scalar(2*ninterface))
allocate(tab_requests_send_recv_vector(2*ninterface))
allocate(b_tab_requests_send_recv_vector(2*ninterface))


allocate(buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
allocate(b_buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
allocate(buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
allocate(b_buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
allocate(buffer_send_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))
allocate(b_buffer_send_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))
allocate(buffer_recv_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))
allocate(b_buffer_recv_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))



allocate(cosrot_irecf(nrecloc))
allocate(sinrot_irecf(nrecloc))


do i=1,nrecloc
cosrot_irecf(i)=sngl(cosrot_irec(i))
sinrot_irecf(i)=sngl(sinrot_irec(i))
enddo


end subroutine init_host_to_dev_variable

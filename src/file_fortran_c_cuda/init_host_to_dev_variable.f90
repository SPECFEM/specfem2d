 subroutine init_host_to_dev_variable(nspec,nglob,poroelastcoef,numat,ninterface,nelemabs,numabs,STACEY_BOUNDARY_CONDITIONS, &
       codeabs,is_proc_source,density,acoustic,elastic,kmato,anisotropy&
         ,which_proc_receiver,fluid_solid_acoustic_ispec, fluid_solid_acoustic_iedge,&
          fluid_solid_elastic_ispec, fluid_solid_elastic_iedge,anisotropic,nrecloc,cosrot_irec,sinrot_irec,displ,veloc,accel,&
          b_displ,b_veloc,b_accel,nglob_elastic,assign_external_model,&
          vpext,vsext,rhoext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext)

  use specfem_par
  implicit none

  integer :: nspec, nglob, numat, ninterface, nelemabs,ispecabs,nrecloc,nglob_elastic
  integer, dimension(nelemabs) :: numabs
  double precision, dimension(4,3,numat) :: poroelastcoef
  logical, dimension(4,nelemabs) :: codeabs
  integer, dimension(NSOURCES) :: is_proc_source
  logical :: STACEY_BOUNDARY_CONDITIONS
  double precision, dimension(2,numat) :: density
  logical, dimension(nspec) :: acoustic, elastic
  integer, dimension(nspec):: kmato
  integer, dimension(nrec):: which_proc_receiver
  double precision, dimension(9,numat) :: anisotropy
  integer, dimension(num_fluid_solid_edges) :: fluid_solid_acoustic_ispec, fluid_solid_acoustic_iedge
  integer, dimension(num_fluid_solid_edges) :: fluid_solid_elastic_ispec, fluid_solid_elastic_iedge
  logical, dimension(nspec) :: anisotropic
  double precision, dimension(nrecloc) :: cosrot_irec,sinrot_irec
  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: displ,veloc,accel
  real(kind=CUSTOM_REAL), dimension(3,nglob) :: b_displ,b_veloc,b_accel
  logical :: assign_external_model
  double precision, dimension(NGLLX,NGLLZ,nspec) :: vpext,vsext,rhoext
  double precision, dimension(NGLLX,NGLLZ,nspec) ::  c11ext,c15ext,c13ext,c33ext,c35ext,c55ext,c12ext,c23ext,c25ext
  real(kind=CUSTOM_REAL) temp

!local
  integer :: i, j, ispec, ier, ispec_inner, ispec_outer,i_spec_free, ispec_elastic, iedge_elastic, ispec_acoustic,&
             iedge_acoustic, inum, iglob, ipoint1d
  real(kind=CUSTOM_REAL) :: xgamma, zgamma, xxi, zxi, jacobian1D
  real(kind=CUSTOM_REAL) :: lambdal_unrelaxed_elastic

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


  allocate(kappastore(NGLLX,NGLLZ,NSPEC_AB))
  allocate(mustore(NGLLX,NGLLZ,NSPEC_AB))
  allocate(rhostore(NGLLX,NGLLZ,NSPEC_AB))
  allocate(rho_vp(NGLLX,NGLLZ,NSPEC_AB))
  allocate(rho_vs(NGLLX,NGLLZ,NSPEC_AB))


if (assign_external_model) then

do ispec=1,nspec
          do j = 1,NGLLZ
              do i = 1,NGLLX

            if(CUSTOM_REAL == SIZE_REAL) then

              rhostore(i,j,ispec)    = sngl (rhoext(i,j,ispec))
              rho_vp(i,j,ispec)      = rhostore(i,j,ispec) * sngl (vpext(i,j,ispec))
              rho_vs(i,j,ispec)      = rhostore(i,j,ispec) * sngl (vsext(i,j,ispec))
              mustore(i,j,ispec)     = rho_vs(i,j,ispec) * sngl (vsext(i,j,ispec))
              kappastore(i,j,ispec)  = rho_vp(i,j,ispec) * sngl (vpext(i,j,ispec))-TWO*TWO*mustore(i,j,ispec)/3._CUSTOM_REAL

            else

              rhostore(i,j,ispec)    = rhoext(i,j,ispec)
              rho_vp(i,j,ispec)      = rhostore(i,j,ispec) * vpext(i,j,ispec)
              rho_vs(i,j,ispec)      = rhostore(i,j,ispec) * vsext(i,j,ispec)
              mustore(i,j,ispec)     = rho_vs(i,j,ispec) * vsext(i,j,ispec)
              kappastore(i,j,ispec)  = rho_vp(i,j,ispec) * vpext(i,j,ispec)-TWO*TWO*mustore(i,j,ispec)/3._CUSTOM_REAL

            endif

                enddo
        enddo
enddo

else ! Internal rho vp vs model

do ispec=1,nspec
          do j = 1,NGLLZ
              do i = 1,NGLLX

            if(CUSTOM_REAL == SIZE_REAL) then
            
              rhostore(i,j,ispec)       = sngl(density(1,kmato(ispec)))
              lambdal_unrelaxed_elastic = sngl(poroelastcoef(1,1,kmato(ispec)))
              mul_unrelaxed_elastic     = sngl(poroelastcoef(2,1,kmato(ispec)))
              mustore(i,j,ispec)        = mul_unrelaxed_elastic
              kappastore(i,j,ispec)     = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic/3._CUSTOM_REAL
              rho_vp(i,j,ispec)         = sngl(density(1,kmato(ispec)) * sqrt((kappastore(i,j,ispec) + &
                                          4._CUSTOM_REAL*mul_unrelaxed_elastic/ &
                                          3._CUSTOM_REAL)/density(1,kmato(ispec))))
              rho_vs(i,j,ispec)         = sngl(density(1,kmato(ispec)) * sqrt(mul_unrelaxed_elastic/density(1,kmato(ispec))))
            else
              rhostore(i,j,ispec)       = density(1,kmato(ispec))
              lambdal_unrelaxed_elastic = poroelastcoef(1,1,kmato(ispec))
              mul_unrelaxed_elastic     = poroelastcoef(2,1,kmato(ispec))
              mustore(i,j,ispec)        = mul_unrelaxed_elastic   
              kappastore(i,j,ispec)     = lambdal_unrelaxed_elastic + TWO*mul_unrelaxed_elastic/3._CUSTOM_REAL
              rho_vp(i,j,ispec)         = density(1,kmato(ispec)) * sqrt((kappastore(i,j,ispec) + &
                                          4._CUSTOM_REAL*mul_unrelaxed_elastic/ &
                                          3._CUSTOM_REAL)/density(1,kmato(ispec)))
              rho_vs(i,j,ispec)         = density(1,kmato(ispec)) * sqrt(mul_unrelaxed_elastic/density(1,kmato(ispec)))
            endif    


                enddo
        enddo
enddo

endif ! Internal/External model





!!!





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
    print*,'host_to_dev: reading in negative nelemabs ',nelemabs,'...resetting to zero'
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



   allocate(islice_selected_source(NSOURCES))
   nsources_local=0

   do i = 1, NSOURCES

      if (is_proc_source(i) == 1) then 
      nsources_local = nsources_local +1 
      islice_selected_source(i) = myrank
      endif
   enddo


  allocate(source_time_function_loc(nsources_local,NSTEP))


  allocate(num_src_loc(NSOURCES))


  j=0
  do i = 1, NSOURCES
      if (is_proc_source(i) == 1) then
        j=j+1
        source_time_function_loc(j,:) = source_time_function(i,:,1)

        num_src_loc(i)=j
      endif
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
!!! Initialisation paramètres pour simulation elastique
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

                c11store(i,j,ispec) = sngl(c11ext(i,j,ispec))
                c13store(i,j,ispec) = sngl(c13ext(i,j,ispec))
                c15store(i,j,ispec) = sngl(c15ext(i,j,ispec))
                c33store(i,j,ispec) = sngl(c33ext(i,j,ispec))
                c35store(i,j,ispec) = sngl(c35ext(i,j,ispec))
                c55store(i,j,ispec) = sngl(c55ext(i,j,ispec))
                c12store(i,j,ispec) = sngl(c12ext(i,j,ispec))
                c23store(i,j,ispec) = sngl(c23ext(i,j,ispec))
                c25store(i,j,ispec) = sngl(c25ext(i,j,ispec))

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
displ_2D(1,:)=displ(1,:)
displ_2D(2,:)=displ(3,:)

veloc_2D(1,:)=veloc(1,:)
veloc_2D(2,:)=veloc(3,:)

accel_2D(1,:)=accel(1,:)
accel_2D(2,:)=accel(3,:)

if(SIMULATION_TYPE == 3 .and. any_elastic) then
  allocate(b_displ_2D(2,nglob))
  allocate(b_veloc_2D(2,nglob))
  allocate(b_accel_2D(2,nglob))


b_displ_2D(1,:)=b_displ(1,:)
b_displ_2D(2,:)=b_displ(3,:)

b_veloc_2D(1,:)=b_veloc(1,:)
b_veloc_2D(2,:)=b_veloc(3,:)

b_accel_2D(1,:)=b_accel(1,:)
b_accel_2D(2,:)=b_accel(3,:)

endif



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate(tab_requests_send_recv_scalar(2*ninterface))
allocate(b_tab_requests_send_recv_scalar(2*ninterface))
allocate(tab_requests_send_recv_vector(2*ninterface))
allocate(b_tab_requests_send_recv_vector(2*ninterface))


allocate(tab_requests_send_recv_vector_c(2*ninterface))

if (.not. CUDA_AWARE_MPI ) then


allocate(buffer_send_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
allocate(buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
allocate(b_buffer_recv_scalar_ext_mesh(max_nibool_interfaces_ext_mesh,ninterface))
allocate(buffer_send_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))
allocate(b_buffer_send_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))
allocate(buffer_recv_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))
allocate(b_buffer_recv_vector_ext_mesh(2,max_nibool_interfaces_ext_mesh,ninterface))


else

allocate(buffer_send_scalar_ext_mesh(1,1))
allocate(buffer_recv_scalar_ext_mesh(1,1))
allocate(b_buffer_recv_scalar_ext_mesh(1,1))
allocate(buffer_send_vector_ext_mesh(1,1,1))
allocate(b_buffer_send_vector_ext_mesh(1,1,1))
allocate(buffer_recv_vector_ext_mesh(1,1,1))
allocate(b_buffer_recv_vector_ext_mesh(1,1,1))

endif



!!

allocate(cosrot_irecf(nrecloc))
allocate(sinrot_irecf(nrecloc))


do i=1,nrecloc
cosrot_irecf(i)=sngl(cosrot_irec(i))
sinrot_irecf(i)=sngl(sinrot_irec(i))
enddo


end subroutine init_host_to_dev_variable

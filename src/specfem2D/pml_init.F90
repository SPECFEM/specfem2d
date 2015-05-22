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

  subroutine pml_init()


! Rough strategy to select the PML damping parameters:

! In order to achieve absorbing for gradient incident wave for zPML:
! we need to set K>1 to bending the wavefront aside from its tangential propagation
! direction.
! Then a bigger damping factor d and alpha could help damping out the
! bended wave along its propagation close to normal direction.

! Keeping in mind, with K>1, we decrease the absorbing efficiency a little
! bit in normal direction.
! Though not exact, you can thought like this: for example with K_x>1, you
! scale the element size in x direction with a big factor.
! Then the resolution decrease in PML but not in interior, thus you can
! see reflections due to different resolutions in two sides of PML interface.
! However if the resolution is already very accurate, that is you have a
! relatively fine grid, then the above effect will be small.


#ifdef USE_MPI
  use mpi
#endif

use specfem_par, only: myrank,SIMULATION_TYPE,SAVE_FORWARD,nspec,nglob,ibool,anyabs,nelemabs,codeabs,numabs,&
                       NELEM_PML_THICKNESS,nspec_PML,is_PML,which_PML_elem,spec_to_PML,region_CPML,&
#ifdef USE_MPI
  ier, &
#endif
                       PML_interior_interface,nglob_interface,mask_ibool,read_external_mesh

  implicit none
  include 'constants.h'

  integer, dimension(nglob) ::   icorner_iglob

  integer :: nspec_PML_tot,ibound,ispecabs,ncorner,i_coef,i,j,k,ispec,iglob

  nspec_PML = 0

  ! detection of PML elements
  if(.not. read_external_mesh) then

    ! ibound is the side we are looking (bottom, right, top or left)
    do ibound=1,4
      icorner_iglob = 0
      ncorner=0

      if(anyabs) then
        ! mark any elements on the boundary as PML and list their corners
        do ispecabs = 1,nelemabs
          ispec = numabs(ispecabs)
          !array to know which PML it is
          which_PML_elem(ibound,ispec)=codeabs(ibound,ispecabs)
          if(codeabs(ibound,ispecabs)) then ! we are on the good absorbing boundary
            do j=1,NGLLZ,NGLLZ-1; do i=1,NGLLX,NGLLX-1
              iglob=ibool(i,j,ispec)
              k=1
              do while(k<=ncorner .and. icorner_iglob(k)/=iglob)
                k=k+1
              enddo
              ncorner=ncorner+1
              icorner_iglob(ncorner) = iglob
            enddo; enddo
          endif ! we are on the good absorbing boundary
        enddo
      endif

     !find elements stuck to boundary elements to define the 4 elements PML thickness
     !we take 4 elements for the PML thickness
     do i_coef=2,NELEM_PML_THICKNESS

       do ispec=1,nspec
         if(.not. which_PML_elem(ibound,ispec)) then
           do j=1,NGLLZ,NGLLZ-1; do i=1,NGLLX,NGLLX-1
             iglob=ibool(i,j,ispec)
             do k=1,ncorner
               if(iglob==icorner_iglob(k)) which_PML_elem(ibound,ispec) = .true.
             enddo
           enddo; enddo
         endif
       enddo

       ! list every corner of each PML element detected
       ncorner=0
       icorner_iglob=0
       nspec_PML=0
       do ispec=1,nspec
          if(which_PML_elem(ibound,ispec)) then
            is_PML(ispec)=.true.
            do j=1,NGLLZ,NGLLZ-1; do i=1,NGLLX,NGLLX-1
              iglob=ibool(i,j,ispec)
              k=1
              do while(k<=ncorner .and. icorner_iglob(k)/=iglob)
                k=k+1
              enddo
              ncorner=ncorner+1
              icorner_iglob(ncorner) = iglob
            enddo; enddo
            nspec_PML=nspec_PML+1
          endif
       enddo

     enddo !end nelem_thickness loop

     if(SIMULATION_TYPE == 3 .or.  (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))then

       do i_coef=NELEM_PML_THICKNESS,NELEM_PML_THICKNESS+1
         do ispec=1,nspec
           if(.not. which_PML_elem(ibound,ispec)) then
             do j=1,NGLLZ,NGLLZ-1; do i=1,NGLLX,NGLLX-1
               iglob=ibool(i,j,ispec)
               do k=1,ncorner
                 if(iglob==icorner_iglob(k)) PML_interior_interface(ibound,ispec) = .true.
               enddo
             enddo; enddo
           endif
         enddo
       enddo !end nelem_thickness loop

     endif !end of SIMULATION_TYPE == 3

   enddo ! end loop on the four boundaries

     if(SIMULATION_TYPE == 3 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))then
       nglob_interface = 0
       do ispec = 1,nspec
         if(PML_interior_interface(IBOTTOM,ispec) .and. (.not. PML_interior_interface(IRIGHT,ispec)) .and. &
            (.not. PML_interior_interface(ILEFT,ispec)) .and. (.not. which_PML_elem(IRIGHT,ispec)) .and. &
            (.not. which_PML_elem(ILEFT,ispec)))then
            nglob_interface = nglob_interface + 5
         else if(PML_interior_interface(ITOP,ispec) .and. (.not. PML_interior_interface(IRIGHT,ispec)) .and. &
                 (.not. PML_interior_interface(ILEFT,ispec)) .and. (.not. which_PML_elem(IRIGHT,ispec)) .and. &
                 (.not. which_PML_elem(ILEFT,ispec)))then
            nglob_interface = nglob_interface + 5
         else if(PML_interior_interface(IRIGHT,ispec) .and. (.not. PML_interior_interface(IBOTTOM,ispec)) .and. &
                 (.not. PML_interior_interface(ITOP,ispec)) .and. (.not. which_PML_elem(IBOTTOM,ispec)) .and. &
                 (.not. which_PML_elem(ITOP,ispec)))then
            nglob_interface = nglob_interface + 5
         else if(PML_interior_interface(ILEFT,ispec) .and. (.not. PML_interior_interface(IBOTTOM,ispec)) .and. &
                 (.not. PML_interior_interface(ITOP,ispec)) .and. (.not. which_PML_elem(IBOTTOM,ispec)) .and. &
                 (.not. which_PML_elem(ITOP,ispec)))then
            nglob_interface = nglob_interface + 5
         else if(PML_interior_interface(ILEFT,ispec) .and. PML_interior_interface(IBOTTOM,ispec))then
            nglob_interface = nglob_interface + 10
         else if(PML_interior_interface(IRIGHT,ispec) .and. PML_interior_interface(IBOTTOM,ispec))then
            nglob_interface = nglob_interface + 10
         else if(PML_interior_interface(ILEFT,ispec) .and. PML_interior_interface(ITOP,ispec))then
            nglob_interface = nglob_interface + 10
         else if(PML_interior_interface(IRIGHT,ispec) .and. PML_interior_interface(ITOP,ispec))then
            nglob_interface = nglob_interface + 10
         endif
       enddo
    endif

   do ispec=1,nspec
     if(is_PML(ispec)) then
! element is in the left cpml layer
       if((which_PML_elem(ILEFT,ispec).eqv. .true.)   .and. (which_PML_elem(IRIGHT,ispec)  .eqv. .false.) .and. &
          (which_PML_elem(ITOP,ispec)  .eqv. .false.) .and. (which_PML_elem(IBOTTOM,ispec).eqv. .false.)) then
         region_CPML(ispec) = CPML_X_ONLY
! element is in the right cpml layer
       else if((which_PML_elem(ILEFT,ispec).eqv. .false.) .and. (which_PML_elem(IRIGHT,ispec)  .eqv. .true.) .and. &
               (which_PML_elem(ITOP,ispec) .eqv. .false.) .and. (which_PML_elem(IBOTTOM,ispec) .eqv. .false.)) then
         region_CPML(ispec) = CPML_X_ONLY
! element is in the top cpml layer
       else if((which_PML_elem(ILEFT,ispec).eqv. .false.) .and. (which_PML_elem(IRIGHT,ispec)  .eqv. .false.) .and. &
               (which_PML_elem(ITOP,ispec) .eqv. .true. ) .and. (which_PML_elem(IBOTTOM,ispec) .eqv. .false.)) then
         region_CPML(ispec) = CPML_Z_ONLY
! element is in the bottom cpml layer
       else if((which_PML_elem(ILEFT,ispec).eqv. .false.) .and. (which_PML_elem(IRIGHT,ispec)  .eqv. .false.) .and. &
               (which_PML_elem(ITOP,ispec) .eqv. .false.) .and. (which_PML_elem(IBOTTOM,ispec) .eqv. .true. )) then
         region_CPML(ispec) = CPML_Z_ONLY
! element is in the left-top cpml corner
       else if((which_PML_elem(ILEFT,ispec).eqv. .true. ) .and. (which_PML_elem(IRIGHT,ispec)  .eqv. .false.).and. &
               (which_PML_elem(ITOP,ispec) .eqv. .true. ) .and. (which_PML_elem(IBOTTOM,ispec) .eqv. .false.)) then
         region_CPML(ispec) = CPML_XZ_ONLY
! element is in the right-top cpml corner
       else if((which_PML_elem(ILEFT,ispec).eqv. .false. ).and. (which_PML_elem(IRIGHT,ispec)  .eqv. .true. ).and. &
               (which_PML_elem(ITOP,ispec) .eqv. .true.  ).and. (which_PML_elem(IBOTTOM,ispec) .eqv. .false.)) then
         region_CPML(ispec) = CPML_XZ_ONLY
! element is in the left-bottom cpml corner
       else if((which_PML_elem(ILEFT,ispec).eqv. .true.  ).and. (which_PML_elem(IRIGHT,ispec)  .eqv. .false.).and. &
               (which_PML_elem(ITOP,ispec) .eqv. .false. ).and. (which_PML_elem(IBOTTOM,ispec) .eqv. .true. )) then
         region_CPML(ispec) = CPML_XZ_ONLY
! element is in the right-bottom cpml corner
       else if((which_PML_elem(ILEFT,ispec).eqv. .false. ).and. (which_PML_elem(IRIGHT,ispec)  .eqv. .true.).and. &
               (which_PML_elem(ITOP,ispec) .eqv. .false. ).and. (which_PML_elem(IBOTTOM,ispec) .eqv. .true.)) then
         region_CPML(ispec) = CPML_XZ_ONLY
       else
         region_CPML(ispec) = 0
       endif
     endif
   enddo

   !construction of table to use less memory for absorbing coefficients
     spec_to_PML=0
     nspec_PML=0
     do ispec=1,nspec
        if (is_PML(ispec)) then
           nspec_PML=nspec_PML+1
           spec_to_PML(ispec)=nspec_PML
        endif
     enddo

  endif !end of detection of element inside PML layer for inner mesher

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if(read_external_mesh) then

    if(.not. allocated(mask_ibool)) allocate(mask_ibool(nglob))

    is_PML(:) = .false.
    which_PML_elem(:,:) = .false.
    nspec_PML = 0
    spec_to_PML=0
    mask_ibool(:) = .false.
    do ispec=1,nspec
      if(region_CPML(ispec) /= 0) then
        nspec_PML = nspec_PML + 1
        is_PML(ispec)=.true.
        spec_to_PML(ispec)=nspec_PML
      endif

      if(SIMULATION_TYPE == 3 .or.  (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))then
        if(region_CPML(ispec) == 0) then
          do i = 1, NGLLX;  do j = 1, NGLLZ
            iglob = ibool(i,j,ispec)
            mask_ibool(iglob) = .true.
          enddo; enddo
        endif
      endif
    enddo

    nglob_interface = 0
    if(SIMULATION_TYPE == 3 .or.  (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))then
      do ispec=1,nspec
        if(region_CPML(ispec) /= 0) then
          do i = 1, NGLLX; do j = 1, NGLLZ
            iglob = ibool(i,j,ispec)
            if(mask_ibool(iglob))nglob_interface = nglob_interface + 1
          enddo; enddo
        endif
      enddo
    endif

  endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#ifdef USE_MPI
  call MPI_REDUCE(nspec_PML, nspec_PML_tot, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
#else
  nspec_PML_tot = nspec_PML
#endif

  if(myrank == 0) then
    write(IOUT,*) "Total number of PML spectral elements: ", nspec_PML_tot
    write(IOUT,*)
  endif

end subroutine pml_init

!
!-------------------------------------------------------------------------------------------------
!
 subroutine determin_interface_pml_interior()

  use specfem_par, only: nglob_interface,nspec,ibool,PML_interior_interface,&
                         which_PML_elem,point_interface,read_external_mesh,mask_ibool,region_CPML
  implicit none

  integer i,j,iglob,ispec
  include 'constants.h'

  nglob_interface = 0

  if(.not. read_external_mesh) then
       do ispec = 1,nspec
         if(PML_interior_interface(IBOTTOM,ispec) .and. (.not. PML_interior_interface(IRIGHT,ispec)) .and. &
            (.not. PML_interior_interface(ILEFT,ispec)) .and. (.not. which_PML_elem(IRIGHT,ispec)) .and. &
            (.not. which_PML_elem(ILEFT,ispec)))then
            point_interface(nglob_interface + 1) = ibool(1,1,ispec)
            point_interface(nglob_interface + 2) = ibool(2,1,ispec)
            point_interface(nglob_interface + 3) = ibool(3,1,ispec)
            point_interface(nglob_interface + 4) = ibool(4,1,ispec)
            point_interface(nglob_interface + 5) = ibool(5,1,ispec)
            nglob_interface = nglob_interface + 5
         else if(PML_interior_interface(ITOP,ispec) .and. (.not. PML_interior_interface(IRIGHT,ispec)) .and. &
                 (.not. PML_interior_interface(ILEFT,ispec)) .and. (.not. which_PML_elem(IRIGHT,ispec)) .and. &
                 (.not. which_PML_elem(ILEFT,ispec)))then
            point_interface(nglob_interface + 1) = ibool(1,NGLLZ,ispec)
            point_interface(nglob_interface + 2) = ibool(2,NGLLZ,ispec)
            point_interface(nglob_interface + 3) = ibool(3,NGLLZ,ispec)
            point_interface(nglob_interface + 4) = ibool(4,NGLLZ,ispec)
            point_interface(nglob_interface + 5) = ibool(5,NGLLZ,ispec)
            nglob_interface = nglob_interface + 5
         else if(PML_interior_interface(IRIGHT,ispec) .and. (.not. PML_interior_interface(IBOTTOM,ispec)) .and. &
                 (.not. PML_interior_interface(ITOP,ispec)) .and. (.not. which_PML_elem(IBOTTOM,ispec)) .and. &
                 (.not. which_PML_elem(ITOP,ispec)))then
            point_interface(nglob_interface + 1) = ibool(NGLLX,1,ispec)
            point_interface(nglob_interface + 2) = ibool(NGLLX,2,ispec)
            point_interface(nglob_interface + 3) = ibool(NGLLX,3,ispec)
            point_interface(nglob_interface + 4) = ibool(NGLLX,4,ispec)
            point_interface(nglob_interface + 5) = ibool(NGLLX,5,ispec)
            nglob_interface = nglob_interface + 5
         else if(PML_interior_interface(ILEFT,ispec) .and. (.not. PML_interior_interface(IBOTTOM,ispec)) .and. &
                 (.not. PML_interior_interface(ITOP,ispec)) .and. (.not. which_PML_elem(IBOTTOM,ispec)) .and. &
                 (.not. which_PML_elem(ITOP,ispec)))then
            point_interface(nglob_interface + 1) = ibool(1,1,ispec)
            point_interface(nglob_interface + 2) = ibool(1,2,ispec)
            point_interface(nglob_interface + 3) = ibool(1,3,ispec)
            point_interface(nglob_interface + 4) = ibool(1,4,ispec)
            point_interface(nglob_interface + 5) = ibool(1,5,ispec)
            nglob_interface = nglob_interface + 5
         else if(PML_interior_interface(ILEFT,ispec) .and. PML_interior_interface(IBOTTOM,ispec))then
            point_interface(nglob_interface + 1) = ibool(1,1,ispec)
            point_interface(nglob_interface + 2) = ibool(1,2,ispec)
            point_interface(nglob_interface + 3) = ibool(1,3,ispec)
            point_interface(nglob_interface + 4) = ibool(1,4,ispec)
            point_interface(nglob_interface + 5) = ibool(1,5,ispec)
            point_interface(nglob_interface + 6) = ibool(1,1,ispec)
            point_interface(nglob_interface + 7) = ibool(2,1,ispec)
            point_interface(nglob_interface + 8) = ibool(3,1,ispec)
            point_interface(nglob_interface + 9) = ibool(4,1,ispec)
            point_interface(nglob_interface + 10)= ibool(5,1,ispec)
            nglob_interface = nglob_interface + 10
         else if(PML_interior_interface(IRIGHT,ispec) .and. PML_interior_interface(IBOTTOM,ispec))then
            point_interface(nglob_interface + 1) = ibool(NGLLX,1,ispec)
            point_interface(nglob_interface + 2) = ibool(NGLLX,2,ispec)
            point_interface(nglob_interface + 3) = ibool(NGLLX,3,ispec)
            point_interface(nglob_interface + 4) = ibool(NGLLX,4,ispec)
            point_interface(nglob_interface + 5) = ibool(NGLLX,5,ispec)
            point_interface(nglob_interface + 6) = ibool(1,1,ispec)
            point_interface(nglob_interface + 7) = ibool(2,1,ispec)
            point_interface(nglob_interface + 8) = ibool(3,1,ispec)
            point_interface(nglob_interface + 9) = ibool(4,1,ispec)
            point_interface(nglob_interface + 10)= ibool(5,1,ispec)
            nglob_interface = nglob_interface + 10
         else if(PML_interior_interface(ILEFT,ispec) .and. PML_interior_interface(ITOP,ispec))then
            point_interface(nglob_interface + 1) = ibool(1,1,ispec)
            point_interface(nglob_interface + 2) = ibool(1,2,ispec)
            point_interface(nglob_interface + 3) = ibool(1,3,ispec)
            point_interface(nglob_interface + 4) = ibool(1,4,ispec)
            point_interface(nglob_interface + 5) = ibool(1,5,ispec)
            point_interface(nglob_interface + 6) = ibool(1,NGLLZ,ispec)
            point_interface(nglob_interface + 7) = ibool(2,NGLLZ,ispec)
            point_interface(nglob_interface + 8) = ibool(3,NGLLZ,ispec)
            point_interface(nglob_interface + 9) = ibool(4,NGLLZ,ispec)
            point_interface(nglob_interface + 10)= ibool(5,NGLLZ,ispec)
            nglob_interface = nglob_interface + 10
         else if(PML_interior_interface(IRIGHT,ispec) .and. PML_interior_interface(ITOP,ispec))then
            point_interface(nglob_interface + 1) = ibool(NGLLX,1,ispec)
            point_interface(nglob_interface + 2) = ibool(NGLLX,2,ispec)
            point_interface(nglob_interface + 3) = ibool(NGLLX,3,ispec)
            point_interface(nglob_interface + 4) = ibool(NGLLX,4,ispec)
            point_interface(nglob_interface + 5) = ibool(NGLLX,5,ispec)
            point_interface(nglob_interface + 6) = ibool(1,NGLLZ,ispec)
            point_interface(nglob_interface + 7) = ibool(2,NGLLZ,ispec)
            point_interface(nglob_interface + 8) = ibool(3,NGLLZ,ispec)
            point_interface(nglob_interface + 9) = ibool(4,NGLLZ,ispec)
            point_interface(nglob_interface + 10)= ibool(5,NGLLZ,ispec)
            nglob_interface = nglob_interface + 10
         endif
       enddo
  endif

  if(read_external_mesh) then
    nglob_interface = 0
    do ispec=1,nspec
      if(region_CPML(ispec) /= 0) then
        do i = 1, NGLLX; do j = 1, NGLLZ
          iglob = ibool(i,j,ispec)
          if(mask_ibool(iglob))then
            nglob_interface = nglob_interface + 1
            point_interface(nglob_interface)= iglob
          endif
        enddo; enddo
      endif
    enddo
  endif

 end subroutine determin_interface_pml_interior
!
!-------------------------------------------------------------------------------------------------
!

 subroutine define_PML_coefficients()

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par, only: f0,elastic,acoustic,&
                         NSOURCES,ispec_selected_source,&
                         nspec,kmato,density,poroelastcoef,&
                         ibool,coord,is_PML,region_CPML,spec_to_PML,&
                         K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  implicit none

  include "constants.h"

!ZN  an automatic adjustment of PML parameter for elongated model.
!ZN  The goal here is to improve the absorbing efficiency of PML for waves with large incident angle,
!ZN  while keep the layer thickness of PML constant.
!ZN  The idea behind this optimization is illustrated at the top the this file.

!ZN  we start this simple automatic adjustment only when you allow to do it, that is
!ZN  The flag PML_parameter_adjustment_for_elongated_model = .true.
!ZN  Since an bigger K only works when the mesh resolution in PML are fine enough,
!ZN  we require an fine mesh to do PML_parameter_adjustment, that is the maximum elment size in PML
!ZN  should be K_MIN_PML times less than the minimum wave velocity divided by dominate frequency.

!ZN  we do not proved and we can not ensure that such an adjustment is optimum, since we do the adjustment based on
!ZN  our own experience in numerical simulation.


!ZN  double precision :: f0_temp
  double precision :: f0_max

! PML fixed parameters to compute parameter in PML
  double precision, parameter :: NPOWER = 2.d0
  double precision, parameter :: Rcoef = 0.001d0

! PML flexible parameters to compute parameter in PML
  double precision :: damping_modified_factor_acoustic,damping_modified_factor_elastic
  double precision :: K_MAX_PML,K_MIN_PML
  double precision :: ALPHA_MAX_PML

! material properties of the elastic medium
  integer i,j,ispec,iglob,ispec_PML,i_source
  double precision :: lambdalplus2mul_relaxed,rhol
  double precision :: d_x, d_z, K_x, K_z, alpha_x, alpha_z
! define an alias for y and z variable names (which are the same)
  double precision :: d0_z_bottom_acoustic, d0_x_right_acoustic, d0_z_top_acoustic, d0_x_left_acoustic
  double precision :: d0_z_bottom_elastic, d0_x_right_elastic, d0_z_top_elastic, d0_x_left_elastic
  double precision :: abscissa_in_PML, abscissa_normalized

  double precision :: thickness_PML_z_min_bottom,thickness_PML_z_max_bottom,&
                      thickness_PML_x_min_right,thickness_PML_x_max_right,&
                      thickness_PML_z_min_top,thickness_PML_z_max_top,&
                      thickness_PML_x_min_left,thickness_PML_x_max_left,&
                      thickness_PML_z_bottom,thickness_PML_x_right,&
                      thickness_PML_z_top,thickness_PML_x_left

  double precision :: xmin, xmax, zmin, zmax, xorigin, zorigin, xval, zval
  double precision :: vpmax_acoustic, vpmax_elastic
  double precision :: xoriginleft, xoriginright, zorigintop, zoriginbottom

  integer :: NSOURCES_glob
  double precision :: averagex_source, averagez_source, averagex_source_sum, averagez_source_sum
  double precision :: rough_estimate_incident_angle

#ifdef USE_MPI
! for MPI and partitioning
  integer :: ier
  double precision :: f0_max_glob
  double precision :: thickness_PML_z_min_bottom_glob,thickness_PML_z_max_bottom_glob,&
                      thickness_PML_x_min_right_glob,thickness_PML_x_max_right_glob,&
                      thickness_PML_z_min_top_glob,thickness_PML_z_max_top_glob,&
                      thickness_PML_x_min_left_glob,thickness_PML_x_max_left_glob
  double precision :: xmin_glob, xmax_glob, zmin_glob, zmax_glob
  double precision :: vpmax_glob_acoustic, vpmax_glob_elastic
#endif

! compute the maximum dominant frequency of all sources
  f0_max = maxval(f0(:))
#ifdef USE_MPI
  call MPI_ALLREDUCE (f0_max, f0_max_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  f0_max = f0_max_glob
#endif
! finish the computation of the maximum dominant frequency of all sources

! reflection coefficient (Inria report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  ALPHA_MAX_PML = PI*f0_max ! from Festa and Vilotte
  K_MAX_PML = 1.0d0 ! from Gedney page 8.11
  K_MIN_PML = 1.0d0
!Due to experience, the d parameter defition according to Festa and Vilotte is small, thus we use damping_modified_factor_acoustic
!to increase the d parameter for PML implementation for acoustic simulation.
  damping_modified_factor_acoustic = 0.5d0
  damping_modified_factor_elastic = 1.0d0

! check that NPOWER is okay
  if(NPOWER < 1) stop 'NPOWER must be greater than 1'

! get minimum and maximum values of mesh coordinates
  xmin = minval(coord(1,:))
  zmin = minval(coord(2,:))
  xmax = maxval(coord(1,:))
  zmax = maxval(coord(2,:))
  if(zmax - zmin < 0.0 .or. xmax - xmin < 0.0) stop 'there are errors in the mesh'
#ifdef USE_MPI
  call MPI_ALLREDUCE (xmin, xmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (zmin, zmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (xmax, xmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (zmax, zmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  xmin = xmin_glob; zmin = zmin_glob
  xmax = xmax_glob; zmax = zmax_glob
#endif

! get the center(origin) of mesh coordinates
  xorigin = xmin + (xmax - xmin)/2.d0
  zorigin = zmin + (zmax - zmin)/2.d0

! determinate the thickness of PML on each side on which PML are activated
  thickness_PML_z_min_bottom=1.d30
  thickness_PML_z_max_bottom=-1.d30

  thickness_PML_x_min_right=1.d30
  thickness_PML_x_max_right=-1.d30

  thickness_PML_z_min_top=1.d30
  thickness_PML_z_max_top=-1.d30

  thickness_PML_x_min_left=1.d30
  thickness_PML_x_max_left=-1.d30

  do ispec=1,nspec
     if(is_PML(ispec)) then
       do j=1,NGLLZ; do i=1,NGLLX
!!!bottom_case
         if(coord(2,ibool(i,j,ispec)) < zorigin) then
           if(region_CPML(ispec) == CPML_Z_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
             thickness_PML_z_max_bottom=max(coord(2,ibool(i,j,ispec)),thickness_PML_z_max_bottom)
             thickness_PML_z_min_bottom=min(coord(2,ibool(i,j,ispec)),thickness_PML_z_min_bottom)
           endif
         endif
!!!right case
         if(coord(1,ibool(i,j,ispec)) > xorigin) then
           if(region_CPML(ispec) == CPML_X_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
             thickness_PML_x_max_right=max(coord(1,ibool(i,j,ispec)),thickness_PML_x_max_right)
             thickness_PML_x_min_right=min(coord(1,ibool(i,j,ispec)),thickness_PML_x_min_right)
           endif
         endif
!!!top case
         if(coord(2,ibool(i,j,ispec)) > zorigin) then
           if(region_CPML(ispec) == CPML_Z_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
             thickness_PML_z_max_top=max(coord(2,ibool(i,j,ispec)),thickness_PML_z_max_top)
             thickness_PML_z_min_top=min(coord(2,ibool(i,j,ispec)),thickness_PML_z_min_top)
           endif
         endif
!!!left case
         if(coord(1,ibool(i,j,ispec)) < xorigin) then
           if(region_CPML(ispec) == CPML_X_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
             thickness_PML_x_max_left=max(coord(1,ibool(i,j,ispec)),thickness_PML_x_max_left)
             thickness_PML_x_min_left=min(coord(1,ibool(i,j,ispec)),thickness_PML_x_min_left)
           endif
         endif
       enddo; enddo
     endif
  enddo

#ifdef USE_MPI
!!!bottom_case
  call MPI_ALLREDUCE (thickness_PML_z_max_bottom, thickness_PML_z_max_bottom_glob, &
                      1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (thickness_PML_z_min_bottom, thickness_PML_z_min_bottom_glob, &
                      1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  thickness_PML_z_max_bottom=thickness_PML_z_max_bottom_glob
  thickness_PML_z_min_bottom=thickness_PML_z_min_bottom_glob

!!!right_case
  call MPI_ALLREDUCE (thickness_PML_x_max_right, thickness_PML_x_max_right_glob, &
                      1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (thickness_PML_x_min_right, thickness_PML_x_min_right_glob, &
                      1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  thickness_PML_x_max_right=thickness_PML_x_max_right_glob
  thickness_PML_x_min_right=thickness_PML_x_min_right_glob

!!!top_case
  call MPI_ALLREDUCE (thickness_PML_z_max_top, thickness_PML_z_max_top_glob, &
                      1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (thickness_PML_z_min_top, thickness_PML_z_min_top_glob, &
                      1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  thickness_PML_z_max_top=thickness_PML_z_max_top_glob
  thickness_PML_z_min_top=thickness_PML_z_min_top_glob

!!!left_case
  call MPI_ALLREDUCE (thickness_PML_x_max_left, thickness_PML_x_max_left_glob, &
                      1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (thickness_PML_x_min_left, thickness_PML_x_min_left_glob, &
                      1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  thickness_PML_x_max_left=thickness_PML_x_max_left_glob
  thickness_PML_x_min_left=thickness_PML_x_min_left_glob
#endif

  thickness_PML_x_left= thickness_PML_x_max_left - thickness_PML_x_min_left
  thickness_PML_x_right= thickness_PML_x_max_right - thickness_PML_x_min_right
  thickness_PML_z_bottom= thickness_PML_z_max_bottom - thickness_PML_z_min_bottom
  thickness_PML_z_top= thickness_PML_z_max_top - thickness_PML_z_min_top

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = thickness_PML_x_left+xmin
  xoriginright = xmax - thickness_PML_x_right
  zoriginbottom = thickness_PML_z_bottom + zmin
  zorigintop = zmax-thickness_PML_z_top

! compute d0 from Inria report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  vpmax_acoustic = 0.0d0
  vpmax_elastic = 0.0d0
  do ispec = 1,nspec
    if(is_PML(ispec)) then
      if(acoustic(ispec))then
! From gmat01.f90 we know, in acoustic region
! lambdalplus2mul_relaxed = kappal  = poroelastcoef(3,1,kmato(ispec)) = rhol * vp_acoustic * vp_acoustic
        lambdalplus2mul_relaxed = poroelastcoef(3,1,kmato(ispec))
        rhol = density(1,kmato(ispec))
        vpmax_acoustic=max(vpmax_acoustic,sqrt(lambdalplus2mul_relaxed/rhol))
      else if(elastic(ispec))then
        ! get relaxed elastic parameters of current spectral element
        lambdalplus2mul_relaxed = poroelastcoef(3,1,kmato(ispec))
        rhol = density(1,kmato(ispec))
        vpmax_elastic=max(vpmax_elastic,sqrt(lambdalplus2mul_relaxed/rhol))
      else
        stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
      endif
    endif
  enddo

! compute the average position of all sources (not the plane wave incident)
  if (PML_parameter_adjustment) then
    averagex_source = 0.0d0
    averagez_source = 0.0d0
    do i_source = 1, NSOURCES
      ispec = ispec_selected_source(i_source)
      do i = 1, NGLLX
        do j = 1, NGLLZ
           averagex_source = averagex_source + coord(1,ibool(i,j,ispec))
           averagez_source = averagez_source + coord(2,ibool(i,j,ispec))
        enddo
      enddo
    enddo
  endif

  averagex_source_sum = 0.0d0
  averagez_source_sum = 0.0d0
  NSOURCES_glob = 0
#ifdef USE_MPI
  if (PML_parameter_adjustment) then
    call MPI_ALLREDUCE (averagex_source, averagex_source_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
    call MPI_ALLREDUCE (averagez_source, averagez_source_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ier)
    call MPI_ALLREDUCE (NSOURCES, NSOURCES_glob, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ier)
  endif
#endif

  if (PML_parameter_adjustment) then
    averagex_source = averagex_source / (NSOURCES*NGLLX*NGLLZ)
    averagez_source = averagez_source / (NSOURCES*NGLLX*NGLLZ)
  endif

#ifdef USE_MPI
  if (PML_parameter_adjustment) then
    averagex_source = averagex_source_sum / (NSOURCES_glob*NGLLX*NGLLZ)
    averagez_source = averagez_source_sum / (NSOURCES_glob*NGLLX*NGLLZ)
  endif
#endif
! finish the computation of the average position of all sources (not the plane wave incident)

#ifdef USE_MPI
  call MPI_ALLREDUCE (vpmax_acoustic, vpmax_glob_acoustic, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (vpmax_elastic, vpmax_glob_elastic, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  vpmax_acoustic=vpmax_glob_acoustic
  vpmax_elastic=vpmax_glob_elastic
#endif

  d0_x_left_acoustic = - (NPOWER + 1) * vpmax_acoustic * log(Rcoef) / (2.d0 * thickness_PML_x_left)
  d0_x_right_acoustic = - (NPOWER + 1) * vpmax_acoustic * log(Rcoef) / (2.d0 * thickness_PML_x_right)
  d0_z_bottom_acoustic = - (NPOWER + 1) * vpmax_acoustic * log(Rcoef) / (2.d0 * thickness_PML_z_bottom)
  d0_z_top_acoustic = - (NPOWER + 1) * vpmax_acoustic * log(Rcoef) / (2.d0 * thickness_PML_z_top)

  d0_x_left_elastic = - (NPOWER + 1) * vpmax_elastic * log(Rcoef) / (2.d0 * thickness_PML_x_left)
  d0_x_right_elastic = - (NPOWER + 1) * vpmax_elastic * log(Rcoef) / (2.d0 * thickness_PML_x_right)
  d0_z_bottom_elastic = - (NPOWER + 1) * vpmax_elastic * log(Rcoef) / (2.d0 * thickness_PML_z_bottom)
  d0_z_top_elastic = - (NPOWER + 1) * vpmax_elastic * log(Rcoef) / (2.d0 * thickness_PML_z_top)

  d_x_store = 0.d0
  d_z_store = 0.d0
  K_x_store = 1.0d0
  K_z_store = 1.0d0
  alpha_x_store = 0.d0
  alpha_z_store = 0.d0

!!!----------------------------------------------------------------------------
!!! without PML_parameter_adjustment
!!!----------------------------------------------------------------------------
  ! define damping profile at the grid points
  if (.not. PML_parameter_adjustment) then
    do ispec = 1,nspec
      ispec_PML = spec_to_PML(ispec)
      if(is_PML(ispec)) then
        do j=1,NGLLZ; do i=1,NGLLX
          d_x = 0.d0; d_z = 0.d0
          K_x = 1.0d0; K_z = 1.0d0
          alpha_x = 0.d0; alpha_z = 0.d0

          iglob=ibool(i,j,ispec)
          ! abscissa of current grid point along the damping profile
          xval = coord(1,iglob)
          zval = coord(2,iglob)

!!!! ---------- bottom edge
          if(zval < zorigin) then
            if(region_CPML(ispec) == CPML_Z_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
              abscissa_in_PML = zoriginbottom - zval
              if(abscissa_in_PML >= 0.d0) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_z_bottom
!ZN                d_z = d0_z_bottom / damping_modified_factor * abscissa_normalized**NPOWER
                if(acoustic(ispec))then
                  d_z = d0_z_bottom_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                else if(elastic(ispec))then
                  d_z = d0_z_bottom_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                else
                  stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                endif
                K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
              else
                d_z = 0.d0; K_z = 1.0d0; alpha_z = 0.d0
              endif

              if(region_CPML(ispec) == CPML_Z_ONLY)then
                d_x = 0.d0; K_x = 1.0d0; alpha_x = 0.d0
              endif
            endif
          endif

!!!! ---------- top edge
          if(zval > zorigin) then
            if(region_CPML(ispec) == CPML_Z_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
              abscissa_in_PML = zval - zorigintop
              if(abscissa_in_PML >= 0.d0) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_z_top
!ZN                d_z = d0_z_top / damping_modified_factor * abscissa_normalized**NPOWER
                if(acoustic(ispec))then
                  d_z = d0_z_top_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                else if(elastic(ispec))then
                  d_z = d0_z_top_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                else
                  stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                endif
                K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
              else
                d_z = 0.d0; K_z = 1.0d0; alpha_z = 0.d0
              endif

              if(region_CPML(ispec) == CPML_Z_ONLY)then
                d_x = 0.d0; K_x = 1.0d0; alpha_x = 0.d0
              endif
            endif
          endif

!!!! ---------- right edge
          if(xval > xorigin) then
            if(region_CPML(ispec) == CPML_X_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
            ! define damping profile at the grid points
              abscissa_in_PML = xval - xoriginright
              if(abscissa_in_PML >= 0.d0) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x_right
!ZN                d_x = d0_x_right / damping_modified_factor * abscissa_normalized**NPOWER
                if(acoustic(ispec))then
                  d_x = d0_x_right_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                else if(elastic(ispec))then
                  d_x = d0_x_right_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                else
                  stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                endif
                K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
              else
                d_x = 0.d0; K_x = 1.0d0; alpha_x = 0.d0
              endif

              if(region_CPML(ispec) == CPML_X_ONLY)then
                d_z = 0.d0; K_z = 1.0d0; alpha_z = 0.d0
              endif
            endif
          endif

!!!! ---------- left edge
          if(xval < xorigin) then
            if(region_CPML(ispec) == CPML_X_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
              abscissa_in_PML = xoriginleft - xval
              if(abscissa_in_PML >= 0.d0) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x_left
!ZN                d_x = d0_x_left / damping_modified_factor * abscissa_normalized**NPOWER
                if(acoustic(ispec))then
                  d_x = d0_x_left_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                else if(elastic(ispec))then
                  d_x = d0_x_left_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                else
                  stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                endif
                K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
              else
                d_x = 0.d0; K_x = 1.0d0; alpha_x = 0.d0
              endif

              if(region_CPML(ispec) == CPML_X_ONLY)then
                 d_z = 0.d0; K_z = 1.0d0; alpha_z = 0.d0
              endif
            endif
          endif

          d_x_store(i,j,ispec_PML) = d_x
          d_z_store(i,j,ispec_PML) = d_z
          K_x_store(i,j,ispec_PML) = K_x
          K_z_store(i,j,ispec_PML) = K_z
          alpha_x_store(i,j,ispec_PML) = alpha_x
          alpha_z_store(i,j,ispec_PML) = alpha_z

        enddo; enddo
      endif
    enddo
  endif
!!!----------------------------------------------------------------------------
!!! without PML_parameter_adjustment
!!!----------------------------------------------------------------------------

!!!----------------------------------------------------------------------------
!!! with PML_parameter_adjustment
!!!----------------------------------------------------------------------------
  ! define damping profile at the grid points
  if (PML_parameter_adjustment) then
    do ispec = 1,nspec
      ispec_PML = spec_to_PML(ispec)
      if(is_PML(ispec)) then

        do j=1,NGLLZ; do i=1,NGLLX
          d_x = 0.d0; d_z = 0.d0
          K_x = 1.0d0; K_z = 1.0d0
          alpha_x = 0.d0; alpha_z = 0.d0

          iglob=ibool(i,j,ispec)
          ! abscissa of current grid point along the damping profile
          xval = coord(1,iglob)
          zval = coord(2,iglob)

!!!! ---------- bottom edge
          if(zval < zorigin) then
            if(region_CPML(ispec) == CPML_Z_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
              abscissa_in_PML = zoriginbottom - zval
              if(abscissa_in_PML >= 0.d0) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_z_bottom
!ZN                d_z = d0_z_bottom / damping_modified_factor * abscissa_normalized**NPOWER
                rough_estimate_incident_angle = &
                      max(abs(xmin-averagex_source)/abs(averagez_source-zoriginbottom),&
                          abs(xmax-averagex_source)/abs(averagez_source-zoriginbottom))
                if(rough_estimate_incident_angle <= 1.0d0)then
                  damping_modified_factor_acoustic = 0.5d0
                  damping_modified_factor_elastic = 1.0d0
                  K_MAX_PML = 1.0d0
                  K_MIN_PML = 1.0d0
                  ALPHA_MAX_PML = 1.0d0
                  if(acoustic(ispec))then
                    d_z = d0_z_bottom_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_z = d0_z_bottom_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif
                else if(rough_estimate_incident_angle > 1.0d0 .and. &
                       rough_estimate_incident_angle <= 6.0d0 )then

                  damping_modified_factor_acoustic = 0.5d0
                  damping_modified_factor_elastic = 0.8d0
                  K_MAX_PML = 2.0d0
                  K_MIN_PML = 1.0d0
                  ALPHA_MAX_PML = 2.5d0
                  if(acoustic(ispec))then
                    d_z = d0_z_bottom_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_z = d0_z_bottom_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif

                else if(rough_estimate_incident_angle > 6.0d0 )then

                  damping_modified_factor_acoustic = 0.4d0
                  damping_modified_factor_elastic = 0.7d0
                  K_MAX_PML = 3.5d0
                  K_MIN_PML = 1.5d0
                  ALPHA_MAX_PML = 4.0d0
                  if(acoustic(ispec))then
                    d_z = d0_z_bottom_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_z = d0_z_bottom_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif

                endif
              else
                d_z = 0.d0; K_z = 1.0d0; alpha_z = 0.d0
              endif

              if(region_CPML(ispec) == CPML_Z_ONLY)then
                d_x = 0.d0; K_x = 1.0d0; alpha_x = 0.d0
              endif
            endif
          endif

!!!! ---------- top edge
          if(zval > zorigin) then
            if(region_CPML(ispec) == CPML_Z_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
              abscissa_in_PML = zval - zorigintop
              if(abscissa_in_PML >= 0.d0) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_z_top
!ZN                d_z = d0_z_top / damping_modified_factor * abscissa_normalized**NPOWER
                rough_estimate_incident_angle =  &
                      max(abs(xmin-averagex_source)/abs(averagez_source-zorigintop),&
                          abs(xmax-averagex_source)/abs(averagez_source-zorigintop))
                if(rough_estimate_incident_angle <= 1.0d0)then
                  damping_modified_factor_acoustic = 0.5d0
                  damping_modified_factor_elastic = 1.0d0
                  K_MAX_PML = 1.0d0
                  K_MIN_PML = 1.0d0
                  ALPHA_MAX_PML = 1.0d0
                  if(acoustic(ispec))then
                    d_z = d0_z_top_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_z = d0_z_top_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif
                else if(rough_estimate_incident_angle > 1.0d0 .and. &
                       rough_estimate_incident_angle <= 6.0d0 )then

                  damping_modified_factor_acoustic = 0.5d0
                  damping_modified_factor_elastic = 0.8d0
                  K_MAX_PML = 2.0d0
                  K_MIN_PML = 1.0d0
                  ALPHA_MAX_PML = 2.5d0
                  if(acoustic(ispec))then
                    d_z = d0_z_top_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_z = d0_z_top_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif

                else if(rough_estimate_incident_angle > 6.0d0 )then

                  damping_modified_factor_acoustic = 0.4d0
                  damping_modified_factor_elastic = 0.7d0
                  K_MAX_PML = 3.5d0
                  K_MIN_PML = 1.5d0
                  ALPHA_MAX_PML = 4.0d0
                  if(acoustic(ispec))then
                    d_z = d0_z_top_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_z = d0_z_top_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_z = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_z = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif
                endif
              else
                d_z = 0.d0; K_z = 1.0d0; alpha_z = 0.d0
              endif

              if(region_CPML(ispec) == CPML_Z_ONLY)then
                d_x = 0.d0; K_x = 1.0d0; alpha_x = 0.d0
              endif
            endif
          endif

!!!! ---------- right edge
          if(xval > xorigin) then
            if(region_CPML(ispec) == CPML_X_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
            ! define damping profile at the grid points
              abscissa_in_PML = xval - xoriginright
              if(abscissa_in_PML >= 0.d0) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x_right
!ZN                d_x = d0_x_right / damping_modified_factor * abscissa_normalized**NPOWER
                rough_estimate_incident_angle = &
                      max(abs(zmin-averagez_source)/abs(averagex_source-xoriginright),&
                          abs(zmax-averagez_source)/abs(averagex_source-xoriginright))
                if(rough_estimate_incident_angle <= 1.0d0)then
                  damping_modified_factor_acoustic = 0.5d0
                  damping_modified_factor_elastic = 1.0d0
                  K_MAX_PML = 1.0d0
                  K_MIN_PML = 1.0d0
                  ALPHA_MAX_PML = 1.0d0
                  if(acoustic(ispec))then
                    d_x = d0_x_right_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_x = d0_x_right_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif
                else if(rough_estimate_incident_angle > 1.0d0 .and. &
                       rough_estimate_incident_angle <= 6.0d0 )then

                  damping_modified_factor_acoustic = 0.5d0
                  damping_modified_factor_elastic = 0.8d0
                  K_MAX_PML = 2.0d0
                  K_MIN_PML = 1.0d0
                  ALPHA_MAX_PML = 2.5d0
                  if(acoustic(ispec))then
                    d_x = d0_x_right_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_x = d0_x_right_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif

                else if(rough_estimate_incident_angle > 6.0d0 )then

                  damping_modified_factor_acoustic = 0.4d0
                  damping_modified_factor_elastic = 0.7d0
                  K_MAX_PML = 3.5d0
                  K_MIN_PML = 1.5d0
                  ALPHA_MAX_PML = 4.0d0
                  if(acoustic(ispec))then
                    d_x = d0_x_right_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_x = d0_x_right_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif
                endif

              else
                d_x = 0.d0; K_x = 1.0d0; alpha_x = 0.d0
              endif

              if(region_CPML(ispec) == CPML_X_ONLY)then
                d_z = 0.d0; K_z = 1.0d0; alpha_z = 0.d0
              endif
            endif
          endif

!!!! ---------- left edge
          if(xval < xorigin) then
            if(region_CPML(ispec) == CPML_X_ONLY  .or. region_CPML(ispec) == CPML_XZ_ONLY) then
              abscissa_in_PML = xoriginleft - xval
              if(abscissa_in_PML >= 0.d0) then
                abscissa_normalized = abscissa_in_PML / thickness_PML_x_left
!ZN                d_x = d0_x_left / damping_modified_factor * abscissa_normalized**NPOWER
                rough_estimate_incident_angle =  &
                      max(abs(zmin-averagez_source)/abs(averagex_source-xoriginleft),&
                          abs(zmax-averagez_source)/abs(averagex_source-xoriginleft))
                if(rough_estimate_incident_angle <= 1.0d0)then
                  damping_modified_factor_acoustic = 0.5d0
                  damping_modified_factor_elastic = 1.0d0
                  K_MAX_PML = 1.0d0
                  K_MIN_PML = 1.0d0
                  ALPHA_MAX_PML = 1.0d0
                  if(acoustic(ispec))then
                    d_x = d0_x_left_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_x = d0_x_left_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif
                else if(rough_estimate_incident_angle > 1.0d0 .and. &
                       rough_estimate_incident_angle <= 6.0d0 )then

                  damping_modified_factor_acoustic = 0.5d0
                  damping_modified_factor_elastic = 0.8d0
                  K_MAX_PML = 2.0d0
                  K_MIN_PML = 1.0d0
                  ALPHA_MAX_PML = 2.5d0
                  if(acoustic(ispec))then
                    d_x = d0_x_left_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_x = d0_x_left_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif

                else if(rough_estimate_incident_angle > 6.0d0 )then

                  damping_modified_factor_acoustic = 0.4d0
                  damping_modified_factor_elastic = 0.7d0
                  K_MAX_PML = 3.5d0
                  K_MIN_PML = 1.5d0
                  ALPHA_MAX_PML = 4.0d0
                  if(acoustic(ispec))then
                    d_x = d0_x_left_acoustic / damping_modified_factor_acoustic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else if(elastic(ispec))then
                    d_x = d0_x_left_elastic / damping_modified_factor_elastic * abscissa_normalized**NPOWER
                    K_x = K_MIN_PML + (K_MAX_PML - 1.0d0) * abscissa_normalized**NPOWER
                    alpha_x = ALPHA_MAX_PML * (1.0d0 - abscissa_normalized)
                  else
                    stop 'PML only implemented for purely elastic or purely acoustic or acoustic/elastic simulation'
                  endif
                endif

              else
                d_x = 0.d0; K_x = 1.0d0; alpha_x = 0.d0
              endif

              if(region_CPML(ispec) == CPML_X_ONLY)then
                 d_z = 0.d0; K_z = 1.0d0; alpha_z = 0.d0
              endif
            endif
          endif

          d_x_store(i,j,ispec_PML) = d_x
          d_z_store(i,j,ispec_PML) = d_z
          K_x_store(i,j,ispec_PML) = K_x
          K_z_store(i,j,ispec_PML) = K_z
          alpha_x_store(i,j,ispec_PML) = alpha_x
          alpha_z_store(i,j,ispec_PML) = alpha_z

        enddo; enddo
      endif
    enddo
  endif
!!!----------------------------------------------------------------------------
!!! with PML_parameter_adjustment
!!!----------------------------------------------------------------------------

 end subroutine define_PML_coefficients


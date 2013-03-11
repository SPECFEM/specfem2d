!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
! Copyright CNRS, INRIA and University of Pau, France,
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

  subroutine pml_init(nspec,nglob,anyabs,ibool,nelemabs,codeabs,numabs,&
                    nspec_PML,is_PML,which_PML_elem,spec_to_PML, &
                    icorner_iglob,NELEM_PML_THICKNESS,&
                    read_external_mesh,region_CPML,&
                    SIMULATION_TYPE,PML_interior_interface,nglob_interface,SAVE_FORWARD)


  implicit none
  include 'constants.h'

#ifdef USE_MPI
  include 'mpif.h'
#endif

  integer ::  SIMULATION_TYPE,nglob_interface

  integer :: nspec,nglob,nelemabs,nspec_PML,NELEM_PML_THICKNESS
  logical :: anyabs
  logical, dimension(4,nspec) :: PML_interior_interface

  integer :: ibound,ispecabs,ncorner,ispec,iglob
  integer :: i,j,k,i_coef

  logical, dimension(4,nspec) :: which_PML_elem
  integer, dimension(nglob) ::   icorner_iglob
  integer, dimension(nelemabs) :: numabs
  logical, dimension(4,nelemabs) :: codeabs
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: spec_to_PML
  logical, dimension(nspec) :: is_PML

!! DK DK for CPML_element_file
  logical :: read_external_mesh
  integer, dimension(nspec) :: region_CPML
  logical :: SAVE_FORWARD

  !!!detection of PML elements

  if(.not. read_external_mesh) then
  nspec_PML = 0

     !ibound is the side we are looking (bottom, right, top or left)
     do ibound=1,4

     icorner_iglob = ZERO
     ncorner=0

     if (anyabs) then
     !mark any elements on the boundary as PML and list their corners
     do ispecabs = 1,nelemabs
        ispec = numabs(ispecabs)

        !array to know which PML it is
        which_PML_elem(ibound,ispec)=codeabs(ibound,ispecabs)

        if(codeabs(ibound,ispecabs)) then ! we are on the good absorbing boundary
        do j=1,NGLLZ,NGLLZ-1
           do i=1,NGLLX,NGLLX-1
              iglob=ibool(i,j,ispec)
              k=1
              do while(k<=ncorner .and. icorner_iglob(k)/=iglob)
                 k=k+1
              end do
              ncorner=ncorner+1
              icorner_iglob(ncorner) = iglob

           enddo
        enddo
        endif ! we are on the good absorbing boundary

     enddo
     end if

     !find elements stuck to boundary elements to define the 4 elements PML thickness
     !we take 4 elements for the PML thickness
     do i_coef=2,NELEM_PML_THICKNESS

        do ispec=1,nspec
           if (.not. which_PML_elem(ibound,ispec)) then
              do j=1,NGLLZ,NGLLZ-1
                 do i=1,NGLLX,NGLLX-1
                    iglob=ibool(i,j,ispec)
                    do k=1,ncorner
                       if (iglob==icorner_iglob(k)) then
                          which_PML_elem(ibound,ispec) = .true.
                       end if
                    end do
                 end do
              end do
           end if
        end do

        !list every corner of each PML elements detected
        ncorner=0
        icorner_iglob=ZERO
        nspec_PML=0
        do ispec=1,nspec
           if (which_PML_elem(ibound,ispec)) then
              is_PML(ispec)=.true.
              do j=1,NGLLZ,NGLLZ-1
                 do i=1,NGLLX,NGLLX-1
                    iglob=ibool(i,j,ispec)
                    k=1
                    do while(k<=ncorner .and. icorner_iglob(k)/=iglob)
                       k=k+1
                    end do
                    ncorner=ncorner+1
                    icorner_iglob(ncorner) = iglob
                 end do
              end do
              nspec_PML=nspec_PML+1
           end if
        end do

     end do !end nelem_thickness loop

     if(SIMULATION_TYPE == 2 .or.  (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))then

      do i_coef=NELEM_PML_THICKNESS,NELEM_PML_THICKNESS+1
        do ispec=1,nspec
           if (.not. which_PML_elem(ibound,ispec)) then
              do j=1,NGLLZ,NGLLZ-1
                 do i=1,NGLLX,NGLLX-1
                    iglob=ibool(i,j,ispec)
                    do k=1,ncorner
                       if (iglob==icorner_iglob(k)) then
                          PML_interior_interface(ibound,ispec) = .true.
                       end if
                    end do
                 end do
              end do
           end if
        end do

      end do !end nelem_thickness loop

     endif !end of SIMULATION_TYPE == 2

     write(IOUT,*) "number of PML spectral elements on side ", ibound,":", nspec_PML

     enddo ! end loop on the 4 boundaries

 if(SIMULATION_TYPE == 2 .or. (SIMULATION_TYPE == 1 .and. SAVE_FORWARD))then
       nglob_interface = 0
       do ispec = 1,nspec
         if(PML_interior_interface(IBOTTOM,ispec) &
            .and. (.not. PML_interior_interface(IRIGHT,ispec)) &
            .and. (.not. PML_interior_interface(ILEFT,ispec))  &
            .and. (.not. which_PML_elem(IRIGHT,ispec)) &
            .and. (.not. which_PML_elem(ILEFT,ispec)))then
            nglob_interface = nglob_interface + 5
         elseif(PML_interior_interface(ITOP,ispec) &
            .and. (.not. PML_interior_interface(IRIGHT,ispec)) &
            .and. (.not. PML_interior_interface(ILEFT,ispec))  &
            .and. (.not. which_PML_elem(IRIGHT,ispec)) &
            .and. (.not. which_PML_elem(ILEFT,ispec)))then
            nglob_interface = nglob_interface + 5
         elseif(PML_interior_interface(IRIGHT,ispec) &
            .and. (.not. PML_interior_interface(IBOTTOM,ispec)) &
            .and. (.not. PML_interior_interface(ITOP,ispec))    &
            .and. (.not. which_PML_elem(IBOTTOM,ispec)) &
            .and. (.not. which_PML_elem(ITOP,ispec)))then
            nglob_interface = nglob_interface + 5
         elseif(PML_interior_interface(ILEFT,ispec) &
            .and. (.not. PML_interior_interface(IBOTTOM,ispec)) &
            .and. (.not. PML_interior_interface(ITOP,ispec))    &
            .and. (.not. which_PML_elem(IBOTTOM,ispec)) &
            .and. (.not. which_PML_elem(ITOP,ispec)))then
            nglob_interface = nglob_interface + 5
         elseif(PML_interior_interface(ILEFT,ispec) &
                .and. PML_interior_interface(IBOTTOM,ispec))then
            nglob_interface = nglob_interface + 10
         elseif(PML_interior_interface(IRIGHT,ispec) &
                .and. PML_interior_interface(IBOTTOM,ispec))then
            nglob_interface = nglob_interface + 10
         elseif(PML_interior_interface(ILEFT,ispec) &
                .and. PML_interior_interface(ITOP,ispec))then
            nglob_interface = nglob_interface + 10
         elseif(PML_interior_interface(IRIGHT,ispec) &
                .and. PML_interior_interface(ITOP,ispec))then
            nglob_interface = nglob_interface + 10
         endif
       enddo
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do ispec=1,nspec
     if(is_PML(ispec)) then

! element is in the left cpml layer
       if( &
         (which_PML_elem(ILEFT,ispec)   .eqv. .true.)  .and. &
         (which_PML_elem(IRIGHT,ispec)  .eqv. .false.) .and. &
         (which_PML_elem(ITOP,ispec)    .eqv. .false.) .and. &
         (which_PML_elem(IBOTTOM,ispec) .eqv. .false.) ) then
         region_CPML(ispec) = CPML_LEFT
! element is in the right cpml layer
       else if( &
         (which_PML_elem(ILEFT,ispec)   .eqv. .false.) .and. &
         (which_PML_elem(IRIGHT,ispec)  .eqv. .true.)  .and. &
         (which_PML_elem(ITOP,ispec)    .eqv. .false.) .and. &
         (which_PML_elem(IBOTTOM,ispec) .eqv. .false.) ) then
         region_CPML(ispec) = CPML_RIGHT
! element is in the top cpml layer
       else if( &
         (which_PML_elem(ILEFT,ispec)   .eqv. .false.) .and. &
         (which_PML_elem(IRIGHT,ispec)  .eqv. .false.) .and. &
         (which_PML_elem(ITOP,ispec)    .eqv. .true. ) .and. &
         (which_PML_elem(IBOTTOM,ispec) .eqv. .false.) ) then
         region_CPML(ispec) = CPML_TOP
! element is in the bottom cpml layer
       else if( &
         (which_PML_elem(ILEFT,ispec)   .eqv. .false.) .and. &
         (which_PML_elem(IRIGHT,ispec)  .eqv. .false.) .and. &
         (which_PML_elem(ITOP,ispec)    .eqv. .false.) .and. &
         (which_PML_elem(IBOTTOM,ispec) .eqv. .true. ) ) then
         region_CPML(ispec) = CPML_BOTTOM
! element is in the left-top cpml corner
       else if( &
         (which_PML_elem(ILEFT,ispec)   .eqv. .true.  ).and. &
         (which_PML_elem(IRIGHT,ispec)  .eqv. .false. ).and. &
         (which_PML_elem(ITOP,ispec)    .eqv. .true.  ).and. &
         (which_PML_elem(IBOTTOM,ispec) .eqv. .false. )) then
         region_CPML(ispec) = CPML_TOP_LEFT
! element is in the right-top cpml corner
       else if( &
         (which_PML_elem(ILEFT,ispec)   .eqv. .false. ).and. &
         (which_PML_elem(IRIGHT,ispec)  .eqv. .true.  ).and. &
         (which_PML_elem(ITOP,ispec)    .eqv. .true.  ).and. &
         (which_PML_elem(IBOTTOM,ispec) .eqv. .false. )) then
         region_CPML(ispec) = CPML_TOP_RIGHT
! element is in the left-bottom cpml corner
       else if( &
         (which_PML_elem(ILEFT,ispec)   .eqv. .true.  ).and. &
         (which_PML_elem(IRIGHT,ispec)  .eqv. .false. ).and. &
         (which_PML_elem(ITOP,ispec)    .eqv. .false. ).and. &
         (which_PML_elem(IBOTTOM,ispec) .eqv. .true.  )) then
         region_CPML(ispec) = CPML_BOTTOM_LEFT
! element is in the right-bottom cpml corner
       else if( &
         (which_PML_elem(ILEFT,ispec)   .eqv. .false. ).and. &
         (which_PML_elem(IRIGHT,ispec)  .eqv. .true.  ).and. &
         (which_PML_elem(ITOP,ispec)    .eqv. .false. ).and. &
         (which_PML_elem(IBOTTOM,ispec) .eqv. .true.  )) then
         region_CPML(ispec) = CPML_BOTTOM_RIGHT
       else

         region_CPML(ispec) = 0

       endif
     endif
  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

     !construction of table to use less memory for absorbing coefficients
     spec_to_PML=0
     nspec_PML=0
     do ispec=1,nspec
        if (is_PML(ispec)) then
           nspec_PML=nspec_PML+1
           spec_to_PML(ispec)=nspec_PML
        end if
     enddo

     write(IOUT,*) "number of PML spectral elements :", nspec_PML

     endif

  if(read_external_mesh) then
  is_PML(:) = .false.
  which_PML_elem(:,:) = .false.

  nspec_PML = 0
  do ispec=1,nspec
    if(region_CPML(ispec) /= 0) then
      nspec_PML = nspec_PML + 1
      is_PML(ispec)=.true.
      spec_to_PML(ispec)=nspec_PML
    endif
  enddo

  write(IOUT,*) "number of PML spectral elements :", nspec_PML

  endif

  end subroutine pml_init
!
!-------------------------------------------------------------------------------------------------
!
 subroutine determin_interface_pml_interior(nglob_interface,nspec,ibool,PML_interior_interface,&
                                            which_PML_elem,point_interface)

  implicit none
  include 'constants.h'

  integer :: nglob_interface, nspec
  logical, dimension(4,nspec) :: PML_interior_interface
  logical, dimension(4,nspec) :: which_PML_elem
  integer :: ispec
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nglob_interface) :: point_interface

  nglob_interface = 0
       do ispec = 1,nspec
         if(PML_interior_interface(IBOTTOM,ispec) &
            .and. (.not. PML_interior_interface(IRIGHT,ispec)) &
            .and. (.not. PML_interior_interface(ILEFT,ispec))  &
            .and. (.not. which_PML_elem(IRIGHT,ispec)) &
            .and. (.not. which_PML_elem(ILEFT,ispec)))then
            point_interface(nglob_interface + 1) = ibool(1,1,ispec)
            point_interface(nglob_interface + 2) = ibool(2,1,ispec)
            point_interface(nglob_interface + 3) = ibool(3,1,ispec)
            point_interface(nglob_interface + 4) = ibool(4,1,ispec)
            point_interface(nglob_interface + 5) = ibool(5,1,ispec)
            nglob_interface = nglob_interface + 5
         elseif(PML_interior_interface(ITOP,ispec) &
            .and. (.not. PML_interior_interface(IRIGHT,ispec)) &
            .and. (.not. PML_interior_interface(ILEFT,ispec))  &
            .and. (.not. which_PML_elem(IRIGHT,ispec)) &
            .and. (.not. which_PML_elem(ILEFT,ispec)))then
            point_interface(nglob_interface + 1) = ibool(1,NGLLZ,ispec)
            point_interface(nglob_interface + 2) = ibool(2,NGLLZ,ispec)
            point_interface(nglob_interface + 3) = ibool(3,NGLLZ,ispec)
            point_interface(nglob_interface + 4) = ibool(4,NGLLZ,ispec)
            point_interface(nglob_interface + 5) = ibool(5,NGLLZ,ispec)
            nglob_interface = nglob_interface + 5
         elseif(PML_interior_interface(IRIGHT,ispec) &
            .and. (.not. PML_interior_interface(IBOTTOM,ispec)) &
            .and. (.not. PML_interior_interface(ITOP,ispec))    &
            .and. (.not. which_PML_elem(IBOTTOM,ispec)) &
            .and. (.not. which_PML_elem(ITOP,ispec)))then
            point_interface(nglob_interface + 1) = ibool(NGLLX,1,ispec)
            point_interface(nglob_interface + 2) = ibool(NGLLX,2,ispec)
            point_interface(nglob_interface + 3) = ibool(NGLLX,3,ispec)
            point_interface(nglob_interface + 4) = ibool(NGLLX,4,ispec)
            point_interface(nglob_interface + 5) = ibool(NGLLX,5,ispec)
            nglob_interface = nglob_interface + 5
         elseif(PML_interior_interface(ILEFT,ispec) &
            .and. (.not. PML_interior_interface(IBOTTOM,ispec)) &
            .and. (.not. PML_interior_interface(ITOP,ispec))    &
            .and. (.not. which_PML_elem(IBOTTOM,ispec)) &
            .and. (.not. which_PML_elem(ITOP,ispec)))then
            point_interface(nglob_interface + 1) = ibool(1,1,ispec)
            point_interface(nglob_interface + 2) = ibool(1,2,ispec)
            point_interface(nglob_interface + 3) = ibool(1,3,ispec)
            point_interface(nglob_interface + 4) = ibool(1,4,ispec)
            point_interface(nglob_interface + 5) = ibool(1,5,ispec)
            nglob_interface = nglob_interface + 5
         elseif(PML_interior_interface(ILEFT,ispec) &
                .and. PML_interior_interface(IBOTTOM,ispec))then
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
         elseif(PML_interior_interface(IRIGHT,ispec) &
                .and. PML_interior_interface(IBOTTOM,ispec))then
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
         elseif(PML_interior_interface(ILEFT,ispec) &
                .and. PML_interior_interface(ITOP,ispec))then
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
         elseif(PML_interior_interface(IRIGHT,ispec) &
                .and. PML_interior_interface(ITOP,ispec))then
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

 end subroutine determin_interface_pml_interior
!
!-------------------------------------------------------------------------------------------------
!

 subroutine define_PML_coefficients(npoin,nspec,is_PML,ibool,coord,&
          region_CPML,kmato,density,poroelastcoef,numat,f0_temp,&
          myrank,&
          K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store,&
          nspec_PML,spec_to_PML)

  implicit none

  include "constants.h"

#ifdef USE_MPI
  include "mpif.h"
#endif

  integer nspec, npoin, i, j,numat, ispec,iglob,nspec_PML
  double precision :: f0_temp

  logical, dimension(nspec) :: is_PML
  integer, dimension(nspec) :: region_CPML
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec_PML) ::  &
                    K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  double precision, dimension(NDIM,npoin) ::  coord
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  double precision, dimension(2,numat) ::  density
  double precision, dimension(4,3,numat) ::  poroelastcoef
  integer, dimension(nspec) :: kmato

  double precision :: d_x, d_z, K_x, K_z, alpha_x, alpha_z
! define an alias for y and z variable names (which are the same)
  double precision :: d0_z_bottom, d0_x_right, d0_z_top, d0_x_left
  double precision :: abscissa_in_PML, abscissa_normalized

  double precision :: thickness_PML_z_min_bottom,thickness_PML_z_max_bottom,&
       thickness_PML_x_min_right,thickness_PML_x_max_right,&
       thickness_PML_z_min_top,thickness_PML_z_max_top,&
       thickness_PML_x_min_left,thickness_PML_x_max_left,&
       thickness_PML_z_bottom,thickness_PML_x_right,&
       thickness_PML_z_top,thickness_PML_x_left

  double precision :: xmin, xmax, zmin, zmax, vpmax, xval, zval
  double precision :: xoriginleft, xoriginright, zorigintop, zoriginbottom

#ifdef USE_MPI
! for MPI and partitioning
  integer  :: ier

  double precision :: thickness_PML_z_min_bottom_glob,thickness_PML_z_max_bottom_glob,&
       thickness_PML_x_min_right_glob,thickness_PML_x_max_right_glob,&
       thickness_PML_z_min_top_glob,thickness_PML_z_max_top_glob,&
       thickness_PML_x_min_left_glob,thickness_PML_x_max_left_glob
  double precision :: xmin_glob, xmax_glob, zmin_glob, zmax_glob, vpmax_glob
#endif

  integer :: myrank

  !PML fixed parameters

! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0

  double precision, parameter :: K_MAX_PML = 1.d0 ! from Gedney page 8.11

  double precision :: ALPHA_MAX_PML

  double precision :: Rcoef

! material properties of the elastic medium
  double precision :: lambdalplus2mul_relaxed

!DK,DK
  integer, dimension(nspec) :: spec_to_PML
  integer :: ispec_PML
  double precision, parameter :: damping_modified_factor = 1.2d0

! reflection coefficient (INRIA report section 6.1) http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  Rcoef = 0.001d0

  ALPHA_MAX_PML = 2.d0*PI*(f0_temp/2.d0) ! from Festa and Vilotte

! check that NPOWER is okay
  if(NPOWER < 1) stop 'NPOWER must be greater than 1'

  ! get minimum and maximum values of mesh coordinates
  xmin = minval(coord(1,:))
  zmin = minval(coord(2,:))
  xmax = maxval(coord(1,:))
  zmax = maxval(coord(2,:))

#ifdef USE_MPI
  call MPI_ALLREDUCE (xmin, xmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (zmin, zmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (xmax, xmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (zmax, zmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  xmin = xmin_glob
  zmin = zmin_glob
  xmax = xmax_glob
  zmax = zmax_glob
#endif

   thickness_PML_z_min_bottom=1.d30
   thickness_PML_z_max_bottom=-1.d30

   thickness_PML_x_min_right=1.d30
   thickness_PML_x_max_right=-1.d30

   thickness_PML_z_min_top=1.d30
   thickness_PML_z_max_top=-1.d30

   thickness_PML_x_min_left=1.d30
   thickness_PML_x_max_left=-1.d30

   do ispec=1,nspec
      if (is_PML(ispec)) then
         do j=1,NGLLZ
            do i=1,NGLLX

!!!bottom_case
               if (region_CPML(ispec) == CPML_BOTTOM_RIGHT .or. region_CPML(ispec) == CPML_BOTTOM_LEFT .or. &
                   region_CPML(ispec) == CPML_BOTTOM) then
                  thickness_PML_z_max_bottom=max(coord(2,ibool(i,j,ispec)),thickness_PML_z_max_bottom)
                  thickness_PML_z_min_bottom=min(coord(2,ibool(i,j,ispec)),thickness_PML_z_min_bottom)
               endif

!!!right case
               if (region_CPML(ispec) == CPML_BOTTOM_RIGHT .or. region_CPML(ispec) == CPML_TOP_RIGHT .or. &
                   region_CPML(ispec) == CPML_RIGHT) then
                  thickness_PML_x_max_right=max(coord(1,ibool(i,j,ispec)),thickness_PML_x_max_right)
                  thickness_PML_x_min_right=min(coord(1,ibool(i,j,ispec)),thickness_PML_x_min_right)
               endif

!!!top case
               if (region_CPML(ispec) == CPML_TOP_RIGHT .or. region_CPML(ispec) == CPML_TOP_LEFT .or. &
                   region_CPML(ispec) == CPML_TOP) then
                  thickness_PML_z_max_top=max(coord(2,ibool(i,j,ispec)),thickness_PML_z_max_top)
                  thickness_PML_z_min_top=min(coord(2,ibool(i,j,ispec)),thickness_PML_z_min_top)
               endif

!!!left case
               if (region_CPML(ispec) == CPML_BOTTOM_LEFT .or. region_CPML(ispec) == CPML_TOP_LEFT .or. &
                   region_CPML(ispec) == CPML_LEFT) then
                  thickness_PML_x_max_left=max(coord(1,ibool(i,j,ispec)),thickness_PML_x_max_left)
                  thickness_PML_x_min_left=min(coord(1,ibool(i,j,ispec)),thickness_PML_x_min_left)
               endif

            enddo
         enddo
      endif
   enddo

#ifdef USE_MPI

!!!bottom
  call MPI_ALLREDUCE (thickness_PML_z_max_bottom, thickness_PML_z_max_bottom_glob, &
       1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (thickness_PML_z_min_bottom, thickness_PML_z_min_bottom_glob, &
       1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  thickness_PML_z_max_bottom=thickness_PML_z_max_bottom_glob
  thickness_PML_z_min_bottom=thickness_PML_z_min_bottom_glob

!!!right
  call MPI_ALLREDUCE (thickness_PML_x_max_right, thickness_PML_x_max_right_glob, &
       1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (thickness_PML_x_min_right, thickness_PML_x_min_right_glob, &
       1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  thickness_PML_x_max_right=thickness_PML_x_max_right_glob
  thickness_PML_x_min_right=thickness_PML_x_min_right_glob

!!!top
  call MPI_ALLREDUCE (thickness_PML_z_max_top, thickness_PML_z_max_top_glob, &
       1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (thickness_PML_z_min_top, thickness_PML_z_min_top_glob, &
       1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  thickness_PML_z_max_top=thickness_PML_z_max_top_glob
  thickness_PML_z_min_top=thickness_PML_z_min_top_glob

!!!left
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

   vpmax = 0
   do ispec = 1,nspec
      if (is_PML(ispec)) then
        ! get relaxed elastic parameters of current spectral element
        lambdalplus2mul_relaxed = poroelastcoef(3,1,kmato(ispec))
        vpmax=max(vpmax,sqrt(lambdalplus2mul_relaxed/density(1,kmato(ispec))))
      endif
   enddo

#ifdef USE_MPI
   call MPI_ALLREDUCE (vpmax, vpmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
   vpmax=vpmax_glob
#endif

! compute d0 from INRIA report section 6.1 http://hal.inria.fr/docs/00/07/32/19/PDF/RR-3471.pdf
  d0_x_left = - (NPOWER + 1) * vpmax * log(Rcoef) / (2.d0 * thickness_PML_x_left)
  d0_x_right = - (NPOWER + 1) * vpmax * log(Rcoef) / (2.d0 * thickness_PML_x_right)
  d0_z_bottom = - (NPOWER + 1) * vpmax * log(Rcoef) / (2.d0 * thickness_PML_z_bottom)
  d0_z_top = - (NPOWER + 1) * vpmax * log(Rcoef) / (2.d0 * thickness_PML_z_top)

!   if (myrank == 0) then
      write(IOUT,*)
      write(IOUT,*) 'PML properties -------',myrank,'myrank'
      write(IOUT,*) '     Vpmax=', vpmax
      write(IOUT,*) '     log(Rcoef)=',log(Rcoef)
      write(IOUT,*) '     thickness_PML_z_bottom =',thickness_PML_z_bottom
      write(IOUT,*) '     thickness_PML_x_right  =',thickness_PML_x_right
      write(IOUT,*) '     thickness_PML_z_top    =',thickness_PML_z_top
      write(IOUT,*) '     thickness_PML_x_left   =',thickness_PML_x_left
      write(IOUT,*) '     d0_bottom       =', d0_z_bottom
      write(IOUT,*) '     d0_right        =', d0_x_right
      write(IOUT,*) '     d0_top          =', d0_z_top
      write(IOUT,*) '     d0_left         =', d0_x_left
!   endif

   d_x = 0.d0
   d_z = 0.d0
   K_x = 0.d0
   K_z = 0.d0
   alpha_x = 0.d0
   alpha_z = 0.d0

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = thickness_PML_x_left+xmin
  xoriginright = xmax - thickness_PML_x_right
  zoriginbottom = thickness_PML_z_bottom + zmin
  zorigintop = zmax-thickness_PML_z_top


 do ispec = 1,nspec
  ispec_PML=spec_to_PML(ispec)
    if (is_PML(ispec)) then
       do j=1,NGLLZ
          do i=1,NGLLX
          K_x_store(i,j,ispec_PML) = 0._CUSTOM_REAL
          K_z_store(i,j,ispec_PML) = 0._CUSTOM_REAL
          d_x_store(i,j,ispec_PML) = 0._CUSTOM_REAL
          d_z_store(i,j,ispec_PML) = 0._CUSTOM_REAL
          alpha_x_store(i,j,ispec_PML) = 0._CUSTOM_REAL
          alpha_z_store(i,j,ispec_PML) = 0._CUSTOM_REAL

             iglob=ibool(i,j,ispec)
             ! abscissa of current grid point along the damping profile
             xval = coord(1,iglob)
             zval = coord(2,iglob)

!!!! ---------- bottom edge
               if (region_CPML(ispec) == CPML_BOTTOM_RIGHT .or. region_CPML(ispec) == CPML_BOTTOM_LEFT .or. &
                   region_CPML(ispec) == CPML_BOTTOM) then
                ! define damping profile at the grid points
                abscissa_in_PML = zoriginbottom - zval
                if(abscissa_in_PML >= 0.d0) then
                   abscissa_normalized = abscissa_in_PML / thickness_PML_z_bottom

                   d_z = d0_z_bottom / damping_modified_factor * abscissa_normalized**NPOWER

!DK DK we keep the equation to define an nonconstant alpha_z in case user needed,
!However for users who want to use nonconstant alpha_z, you also need to change
!routines for CMPL computation. For example in compute_forces_viscoelastic.f90

!                   alpha_z = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) &
!                   + ALPHA_MAX_PML / 2.d0

!DK DK Here we set alpha_z=alpha_x=const where alpha_z or alpha_x is nonzero

                   alpha_z = ALPHA_MAX_PML / 2.d0

                   K_z = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
                else
                   d_z = 0.d0
                   alpha_z = 0.d0
                   K_z = 1.d0
                endif
             endif

!!!! ---------- top edge
               if (region_CPML(ispec) == CPML_TOP_RIGHT .or. region_CPML(ispec) == CPML_TOP_LEFT .or. &
                   region_CPML(ispec) == CPML_TOP) then
                ! define damping profile at the grid points
                abscissa_in_PML = zval - zorigintop
                if(abscissa_in_PML >= 0.d0) then
                   abscissa_normalized = abscissa_in_PML / thickness_PML_z_top

                   d_z = d0_z_top / damping_modified_factor * abscissa_normalized**NPOWER

!                   alpha_z = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) &
!                   + ALPHA_MAX_PML / 2.d0

                   alpha_z = ALPHA_MAX_PML / 2.d0

                   K_z = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
                else
                   d_z = 0.d0
                   alpha_z = 0.d0
                   K_z = 1.d0
                endif
             endif

!!!! ---------- right edge
               if (region_CPML(ispec) == CPML_BOTTOM_RIGHT .or. region_CPML(ispec) == CPML_TOP_RIGHT .or. &
                   region_CPML(ispec) == CPML_RIGHT) then
                ! define damping profile at the grid points
                abscissa_in_PML = xval - xoriginright
                if(abscissa_in_PML >= 0.d0) then
                   abscissa_normalized = abscissa_in_PML / thickness_PML_x_right

                   d_x = d0_x_right / damping_modified_factor * abscissa_normalized**NPOWER

!                   alpha_x = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) &
!                   + ALPHA_MAX_PML / 2.d0

                   alpha_x = ALPHA_MAX_PML / 2.d0

                   K_x = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
                else
                   d_x = 0.d0
                   alpha_x = 0.d0
                   K_x = 1.d0
                endif
             endif

!!!! ---------- left edge
               if (region_CPML(ispec) == CPML_BOTTOM_LEFT .or. region_CPML(ispec) == CPML_TOP_LEFT .or. &
                   region_CPML(ispec) == CPML_LEFT) then
                ! define damping profile at the grid points
                abscissa_in_PML = xoriginleft - xval
                if(abscissa_in_PML >= 0.d0) then
                   abscissa_normalized = abscissa_in_PML / thickness_PML_x_left

                   d_x = d0_x_left / damping_modified_factor * abscissa_normalized**NPOWER

!                   alpha_x = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) &
!                   + ALPHA_MAX_PML / 2.d0

                   alpha_x = ALPHA_MAX_PML / 2.d0

                   K_x = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
                else
                   d_x = 0.d0
                   alpha_x = 0.d0
                   K_x = 1.d0
                endif
             endif

       if(CUSTOM_REAL == SIZE_REAL) then
          K_x_store(i,j,ispec_PML) = sngl(K_x)
          K_z_store(i,j,ispec_PML) = sngl(K_z)
          d_x_store(i,j,ispec_PML) = sngl(d_x)
          d_z_store(i,j,ispec_PML) = sngl(d_z)
          alpha_x_store(i,j,ispec_PML) = sngl(alpha_x)
          alpha_z_store(i,j,ispec_PML) = sngl(alpha_z)
       else
          K_x_store(i,j,ispec_PML) = K_x
          K_z_store(i,j,ispec_PML) = K_z
          d_x_store(i,j,ispec_PML) = d_x
          d_z_store(i,j,ispec_PML) = d_z
          alpha_x_store(i,j,ispec_PML) = alpha_x
          alpha_z_store(i,j,ispec_PML) = alpha_z
       endif

       enddo
     enddo
    endif
 enddo

  end subroutine define_PML_coefficients


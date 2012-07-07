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
                    nspec_PML,is_PML,which_PML_elem,which_PML_poin,spec_to_PML,ibool_PML, &
                    npoin_PML,icorner_iglob,NELEM_PML_THICKNESS)

  implicit none
  include 'constants.h'

#ifdef USE_MPI
  include 'mpif.h'
#endif

  integer :: nspec,nglob,nelemabs,nspec_PML,npoin_PML,NELEM_PML_THICKNESS
  logical :: anyabs

  integer :: ibound,ispecabs,ncorner,ispec,iglob
  integer :: i,j,k,i_coef

  logical, dimension(4,nspec) :: which_PML_elem
  logical, dimension(4,nglob) :: which_PML_poin
  integer, dimension(nglob) ::   icorner_iglob
  integer, dimension(nelemabs) :: numabs
  logical, dimension(4,nelemabs) :: codeabs
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: spec_to_PML
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool_PML
  integer, dimension(:), allocatable :: iPML_to_iglob
  logical, dimension(nspec) :: is_PML

  !!!detection of PML elements

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

              !array to know which PML it is
              which_PML_poin(ibound,iglob) = codeabs(ibound,ispecabs)

           enddo
        enddo
        endif ! we are on the good absorbing boundary

     enddo
     end if

     !find elements stuck to boundary elements to define the 4 elements PML thickness
     !we take 4 elements for the PML thickness
     do i_coef=2,NELEM_PML_THICKNESS

        do ispec=1,nspec
           if (.not.which_PML_elem(ibound,ispec)) then
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
                    which_PML_poin(ibound,iglob) = .true.
                 end do
              end do
              nspec_PML=nspec_PML+1
           end if
        end do

     end do !end nelem_thickness loop

     write(IOUT,*) "number of PML spectral elements on side ", ibound,":", nspec_PML

     enddo ! end loop on the 4 boundaries

     !construction of table to use less memory for absorbing coefficients
  !     allocate(spec_to_PML(nspec))
     spec_to_PML=0
     nspec_PML=0
     do ispec=1,nspec
        if (is_PML(ispec)) then
           nspec_PML=nspec_PML+1
           spec_to_PML(ispec)=nspec_PML
        end if
     enddo

  !     allocate(ibool_PML(NGLLX,NGLLZ,nspec))
     if (nspec_PML > 0) then
        allocate(iPML_to_iglob(nspec_PML*NGLLX*NGLLZ))
     else
        allocate(iPML_to_iglob(1))
     end if

     iPML_to_iglob(:)=0
     ibool_PML=0

     npoin_PML=0
     do ispec=1,nspec
        if (is_PML(ispec)) then
           do j=1,NGLLZ
              do i=1,NGLLX
                 iglob=ibool(i,j,ispec)
                 k=1
                 do while (iglob/=iPML_to_iglob(k).and.iPML_to_iglob(k)/=0)
                    k=k+1
                 end do
                 ibool_PML(i,j,ispec)=k
                 if (k>npoin_PML) then
                    npoin_PML=npoin_PML+1
                    iPML_to_iglob(k)=iglob
                 end if
              end do
           end do
        end if
     end do

     deallocate(iPML_to_iglob)

     write(IOUT,*) "number of PML spectral elements :", nspec_PML
     write(IOUT,*) "number of PML spectral points   :", npoin_PML

  end subroutine pml_init

!
!-------------------------------------------------------------------------------------------------
!

 subroutine define_PML_coefficients(npoin,nspec,is_PML,ibool,coord,&
          which_PML_elem,kmato,density,poroelastcoef,numat,f0_temp,npoin_PML,&
          ibool_PML,myrank,&
            K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store)

  implicit none

  include "constants.h"

#ifdef USE_MPI
  include "mpif.h"
#endif

  integer nspec, npoin, i, j,numat, ispec,iglob,npoin_PML,iPML
  double precision :: f0_temp

  logical, dimension(nspec) :: is_PML
  logical, dimension(4,nspec) :: which_PML_elem
  real(kind=CUSTOM_REAL), dimension(npoin_PML) ::  &
                    K_x_store,K_z_store,d_x_store,d_z_store,alpha_x_store,alpha_z_store

  real(kind=CUSTOM_REAL), dimension(NDIM,npoin) ::  coord
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool, ibool_PML
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
!       thickness_PML_z_bottom_glob,thickness_PML_x_right_glob,&
!       thickness_PML_z_top_glob,thickness_PML_x_left_glob

  double precision :: xmin_glob, xmax_glob, zmin_glob, zmax_glob, vpmax_glob
#endif

  integer :: myrank

  !PML fixed parameters

! power to compute d0 profile
  double precision, parameter :: NPOWER = 2.d0

  double precision, parameter :: K_MAX_PML = 1.d0 ! from Gedney page 8.11

  double precision :: ALPHA_MAX_PML

  double precision :: Rcoef

  double precision :: factorx, factorz

! material properties of the elastic medium
  double precision :: lambdalplus2mul_relaxed

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
               if (which_PML_elem(IBOTTOM,ispec)) then
                  thickness_PML_z_max_bottom=max(coord(2,ibool(i,j,ispec)),thickness_PML_z_max_bottom)
                  thickness_PML_z_min_bottom=min(coord(2,ibool(i,j,ispec)),thickness_PML_z_min_bottom)
               endif

!!!right case
               if (which_PML_elem(IRIGHT,ispec)) then
                  thickness_PML_x_max_right=max(coord(1,ibool(i,j,ispec)),thickness_PML_x_max_right)
                  thickness_PML_x_min_right=min(coord(1,ibool(i,j,ispec)),thickness_PML_x_min_right)
               endif

!!!top case
               if (which_PML_elem(ITOP,ispec)) then
                  thickness_PML_z_max_top=max(coord(2,ibool(i,j,ispec)),thickness_PML_z_max_top)
                  thickness_PML_z_min_top=min(coord(2,ibool(i,j,ispec)),thickness_PML_z_min_top)
               endif

!!!left case
               if (which_PML_elem(ILEFT,ispec)) then
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

   if (myrank == 0) then
      write(IOUT,*)
      write(IOUT,*) 'PML properties -------'
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
   endif

   d_x = ZERO
   d_z = ZERO
   K_x = ZERO
   K_z = ZERO
   alpha_x = ZERO
   alpha_z = ZERO

! origin of the PML layer (position of right edge minus thickness, in meters)
  xoriginleft = thickness_PML_x_left+xmin
  xoriginright = xmax - thickness_PML_x_right
  zoriginbottom = thickness_PML_z_bottom + zmin
  zorigintop = zmax-thickness_PML_z_top


 do ispec = 1,nspec

    if (is_PML(ispec)) then

       do j=1,NGLLZ
          do i=1,NGLLX
             iglob=ibool(i,j,ispec)
             iPML=ibool_PML(i,j,ispec)
             if(iPML < 1) stop 'error: iPML < 1 in a PML element'
             ! abscissa of current grid point along the damping profile
             xval = coord(1,iglob)
             zval = coord(2,iglob)

!!!! ---------- bottom edge
             if (which_PML_elem(IBOTTOM,ispec)) then
                ! define damping profile at the grid points
                abscissa_in_PML = zoriginbottom - zval
                if(abscissa_in_PML >= 0.d0) then
                   abscissa_normalized = abscissa_in_PML / thickness_PML_z_bottom
                   d_z = d0_z_bottom / 0.6d0 * abscissa_normalized**NPOWER

                   alpha_z = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) &
                   + ALPHA_MAX_PML / 2.d0

                   K_z = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
                else
                   d_z = 0.d0
                   alpha_z = 0.d0
                   K_z = 1.d0
                endif
                factorz = 1.d0
                if (which_PML_elem(IRIGHT,ispec) .or. which_PML_elem(ILEFT,ispec)) then
                   factorx = 1.d0
                else
                   factorx = 0.d0
                endif
             endif

!!!! ---------- top edge
             if (which_PML_elem(ITOP,ispec)) then
                ! define damping profile at the grid points
                abscissa_in_PML = zval - zorigintop
                if(abscissa_in_PML >= 0.d0) then
                   abscissa_normalized = abscissa_in_PML / thickness_PML_z_top
                   d_z = d0_z_top / 0.6d0 * abscissa_normalized**NPOWER
                   alpha_z = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) &
                   + ALPHA_MAX_PML / 2.d0

                   K_z = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
                else
                   d_z = 0.d0
                   alpha_z = 0.d0
                   K_z = 1.d0
                endif
                factorz = 1.d0
                if (which_PML_elem(IRIGHT,ispec) .or. which_PML_elem(ILEFT,ispec)) then
                   factorx = 1.d0
                else
                   factorx = 0.d0
                endif
             endif

!!!! ---------- right edge
             if (which_PML_elem(IRIGHT,ispec)) then
                ! define damping profile at the grid points
                abscissa_in_PML = xval - xoriginright
                if(abscissa_in_PML >= 0.d0) then
                   abscissa_normalized = abscissa_in_PML / thickness_PML_x_right
                   d_x = d0_x_right / 0.6d0 * abscissa_normalized**NPOWER
                   alpha_x = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) &
                   + ALPHA_MAX_PML / 2.d0

                   K_x = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
                else
                   d_x = 0.d0
                   alpha_x = 0.d0
                   K_x = 1.d0
                endif
                factorx = 1.d0
                if (which_PML_elem(IBOTTOM,ispec) .or. which_PML_elem(ITOP,ispec)) then
                   factorz = 1.d0
                else
                   factorz = 0.d0
                endif
             endif

!!!! ---------- left edge
             if (which_PML_elem(ILEFT,ispec)) then
                ! define damping profile at the grid points
                abscissa_in_PML = xoriginleft - xval
                if(abscissa_in_PML >= 0.d0) then
                   abscissa_normalized = abscissa_in_PML / thickness_PML_x_left
                   d_x = d0_x_left / 0.6d0 * abscissa_normalized**NPOWER
                   alpha_x = ALPHA_MAX_PML * (1.d0 - abscissa_normalized) &
                   + ALPHA_MAX_PML / 2.d0

                   K_x = 1.d0 + (K_MAX_PML - 1.d0) * abscissa_normalized**NPOWER
                else
                   d_x = 0.d0
                   alpha_x = 0.d0
                   K_x = 1.d0
                endif
                factorx = 1.d0
                if (which_PML_elem(IBOTTOM,ispec) .or. which_PML_elem(ITOP,ispec)) then
                   factorz = 1.d0
                else
                   factorz = 0.d0
                endif
             endif

          K_x_store(iPML) = K_x
          K_z_store(iPML) = K_z
          d_x_store(iPML) = d_x
          d_z_store(iPML) = d_z
          alpha_x_store(iPML) = alpha_x
          alpha_z_store(iPML) = alpha_z

       enddo
     enddo
    endif
 enddo

  end subroutine define_PML_coefficients



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

  subroutine checkgrid()

! check the mesh, stability and number of points per wavelength

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par, only: vpext,vsext,rhoext,density,poroelastcoef, &
                      porosity,tortuosity,permeability,ibool,kmato, &
                      coord,vpImin,vpImax,vpIImin,vpIImax, &
                      assign_external_model,nspec,UPPER_LIMIT_DISPLAY,deltat, &
                      f0,initialfield,time_function_type, &
                      coorg,xinterp,zinterp,shape2D_display,knods,simulation_title, &
                      pointsdisp,ngnod,any_elastic,any_poroelastic,all_anisotropic, &
                      myrank,nproc,NSOURCES,poroelastic, &
                      freq0,Q0,ATTENUATION_PORO_FLUID_PART,US_LETTER,output_postscript_snapshot
  implicit none

  include "constants.h"

! option to display only part of the mesh and not the whole mesh,
! for instance to analyze Cuthill-McKee mesh partitioning etc.

  ! local parameters
  double precision vpIImax_local,vpIImin_local
  double precision vsmin,vsmax,densmin,densmax,vpImax_local,vpImin_local,vsmin_local
  double precision kappa_s,kappa_f,kappa_fr,mu_s,mu_fr,denst_s,denst_f,denst,phi,tort,cpIloc,cpIIloc,csloc
  double precision D_biot,H_biot,C_biot,M_biot,cpIsquare,cpIIsquare,cssquare
  double precision f0max,w_c,eta_f,perm
  double precision lambdaplus2mu,mu
  double precision distance_min,distance_max,distance_min_local,distance_max_local,lambdaS_local,lambdaPI_local
  double precision courant_stability_number_max,lambdaPImin,lambdaPImax,lambdaPIImin,lambdaPIImax,lambdaSmin,lambdaSmax
  double precision distance_1,distance_2,distance_3,distance_4

! for the stability condition
! maximum polynomial degree for which we can compute the stability condition
  integer, parameter :: NGLLX_MAX_STABILITY = 15
  double precision :: percent_GLL(NGLLX_MAX_STABILITY)

! color palette
  integer, parameter :: NUM_COLORS = 236
  double precision, dimension(NUM_COLORS) :: red,green,blue

  double precision :: xmax,zmax,height,usoffset,sizex,sizez,courant_stability_number
  double precision :: x1,z1,x2,z2,ratio_page,xmin,zmin

#ifdef USE_MPI
  integer  :: icol
  double precision  :: vpImin_glob,vpImax_glob,vsmin_glob,vsmax_glob,densmin_glob,densmax_glob
  double precision  :: vpIImin_glob,vpIImax_glob
  double precision  :: distance_min_glob,distance_max_glob
  double precision  :: courant_stability_max_glob,lambdaPImin_glob,lambdaPImax_glob,&
                       lambdaPIImin_glob,lambdaPIImax_glob,lambdaSmin_glob,lambdaSmax_glob, &
                       lambdaPmin_in_fluid_histo_glob,lambdaPmax_in_fluid_histo_glob
  double precision  :: xmin_glob, xmax_glob, zmin_glob, zmax_glob
#endif

  logical  :: any_elastic_glob,any_poroelastic_glob
  double precision, dimension(2,nspec*5)  :: coorg_send
  double precision, dimension(:,:), allocatable  :: coorg_recv
  integer, dimension(nspec)  :: RGB_send
  integer, dimension(:), allocatable  :: RGB_recv
  real, dimension(nspec)  :: greyscale_send
  real, dimension(:), allocatable  :: greyscale_recv
  integer :: nspec_recv
  integer :: num_ispec
  integer :: iproc
  integer :: ier
  integer :: i,j,ispec,material
  integer :: is,ir,in,nnum

! for histogram of number of points per wavelength
  double precision :: lambdaPmin_in_fluid_histo,lambdaPmax_in_fluid_histo,lambdaSmin_histo,lambdaSmax_histo, &
                   min_nb_of_points_per_wavelength,max_nb_of_points_per_wavelength,nb_of_points_per_wavelength, &
                   scaling_factor,scaling_factor_S,scaling_factor_P
  logical :: any_fluid_histo,any_fluid_histo_glob
  integer, parameter :: NCLASSES = 20
  integer, dimension(0:NCLASSES-1) :: classes_wavelength,classes_wavelength_all
  integer :: iclass,nspec_all,ipass,nspec_counted,nspec_counted_all,nspec_counted_all_solid,nspec_counted_all_fluid
  logical :: create_wavelength_histogram
  double precision :: current_percent,total_percent

  ! check
  if(UPPER_LIMIT_DISPLAY > nspec) &
    call exit_MPI('cannot have UPPER_LIMIT_DISPLAY > nspec in checkgrid.F90')

  if (myrank == 0) then
    write(IOUT,*) "  checking mesh and stability"
    call flush_IOUT()
  endif

#ifndef USE_MPI
  allocate(coorg_recv(1,1))
  allocate(RGB_recv(1))
  allocate(greyscale_recv(1))
  nspec_recv = 0
  ier = 0
  iproc = nproc
  deallocate(coorg_recv)
  deallocate(RGB_recv)
  deallocate(greyscale_recv)
#endif


! define percentage of smallest distance between GLL points for NGLLX points
! percentages were computed by calling the GLL points routine for each degree
  call checkgrid_setup_GLLper(percent_GLL,NGLLX_MAX_STABILITY)

! define color palette in random order
  call checkgrid_setup_colorp(red,green,blue,NUM_COLORS)

!---- compute parameters for the spectral elements

  vpImin = HUGEVAL
  vpImax = -HUGEVAL

  if(any_elastic .or. any_poroelastic) then
    vsmin = HUGEVAL
    vsmax = -HUGEVAL
  else
    vsmin = 0
    vsmax = 0
  endif

  if(any_poroelastic) then
    vpIImin = HUGEVAL
    vpIImax = -HUGEVAL
  else
    vpIImin = 0
    vpIImax = 0
  endif

  densmin = HUGEVAL
  densmax = -HUGEVAL

  distance_min = HUGEVAL
  distance_max = -HUGEVAL

  courant_stability_number_max = -HUGEVAL

  lambdaPImin = HUGEVAL
  lambdaPImax = -HUGEVAL

  if(any_elastic .or. any_poroelastic) then
    lambdaSmin = HUGEVAL
    lambdaSmax = -HUGEVAL
  else
    lambdaSmin = 0
    lambdaSmax = 0
  endif

  if(any_poroelastic) then
    lambdaPIImin = HUGEVAL
    lambdaPIImax = -HUGEVAL
  else
    lambdaPIImin = 0
    lambdaPIImax = 0
  endif

  lambdaPmin_in_fluid_histo = HUGEVAL
  lambdaPmax_in_fluid_histo = -HUGEVAL

  any_fluid_histo = .false.

  do ispec=1,nspec

    material = kmato(ispec)

    if(poroelastic(ispec)) then

      ! poroelastic material

      phi = porosity(material)
      tort = tortuosity(material)
      perm = permeability(1,material)
      ! solid properties
      mu_s = poroelastcoef(2,1,material)
      kappa_s = poroelastcoef(3,1,material) - FOUR_THIRDS*mu_s
      denst_s = density(1,material)
      denst = denst_s
      ! fluid properties
      kappa_f = poroelastcoef(1,2,material)
      denst_f = density(2,material)
      eta_f = poroelastcoef(2,2,material)
      ! frame properties
      mu_fr = poroelastcoef(2,3,material)
      kappa_fr = poroelastcoef(3,3,material) - FOUR_THIRDS*mu_fr
      ! Biot coefficients for the input phi
      D_biot = kappa_s*(1.d0 + phi*(kappa_s/kappa_f - 1.d0))
      H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr + FOUR_THIRDS*mu_fr
      C_biot = kappa_s*(kappa_s - kappa_fr)/(D_biot - kappa_fr)
      M_biot = kappa_s*kappa_s/(D_biot - kappa_fr)

      call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
             tort,denst_s,denst_f,eta_f,perm,f0(1),freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

      cpIloc = sqrt(cpIsquare)
      cpIIloc = sqrt(cpIIsquare)
      csloc = sqrt(cssquare)
    else
      mu = poroelastcoef(2,1,material)
      lambdaplus2mu  = poroelastcoef(3,1,material)
      denst = density(1,material)

      cpIloc = sqrt(lambdaplus2mu/denst)
      cpIIloc = 0.d0
      csloc = sqrt(mu/denst)
    endif

    vpImax_local = -HUGEVAL
    vpImin_local = HUGEVAL
    vpIImax_local = -HUGEVAL
    vpIImin_local = HUGEVAL
    vsmin_local = HUGEVAL

    distance_min_local = HUGEVAL
    distance_max_local = -HUGEVAL

    do j=1,NGLLZ
      do i=1,NGLLX

!--- if heterogeneous formulation with external velocity model
        if(assign_external_model) then
          cpIloc = vpext(i,j,ispec)
          csloc = vsext(i,j,ispec)
          denst = rhoext(i,j,ispec)
        endif

!--- compute min and max of velocity and density models
        vpImin = min(vpImin,cpIloc)
        vpImax = max(vpImax,cpIloc)

! ignore acoustic and elastic regions with cpII = 0
        if(cpIIloc > 1.d-20) vpIImin = min(vpIImin,cpIIloc)
        vpIImax = max(vpIImax,cpIIloc)

! ignore fluid regions with Vs = 0
        if(csloc > 1.d-20) vsmin = min(vsmin,csloc)
        vsmax = max(vsmax,csloc)

        densmin = min(densmin,denst)
        densmax = max(densmax,denst)

        vpImax_local = max(vpImax_local,cpIloc)
        vpImin_local = min(vpImin_local,cpIloc)
        vpIImax_local = max(vpIImax_local,cpIIloc)
        vpIImin_local = min(vpIImin_local,cpIIloc)
        vsmin_local = min(vsmin_local,csloc)

      enddo
    enddo

! compute minimum and maximum size of edges of this grid cell
    distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
               (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

    distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

    distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

    distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
               (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

    distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
    distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

    distance_min = min(distance_min,distance_min_local)
    distance_max = max(distance_max,distance_max_local)

    courant_stability_number_max = max(courant_stability_number_max, &
                vpImax_local * deltat / (distance_min_local * percent_GLL(NGLLX)))

! check if fluid region with Vs = 0
    if(vsmin_local > 1.d-20) then
      lambdaSmin = min(lambdaSmin,vsmin_local / (distance_max_local / (NGLLX - 1)))
      lambdaSmax = max(lambdaSmax,vsmin_local / (distance_max_local / (NGLLX - 1)))
    else
      any_fluid_histo = .true.
      lambdaPmin_in_fluid_histo = min(lambdaPmin_in_fluid_histo,vpImin_local / (distance_max_local / (NGLLX - 1)))
      lambdaPmax_in_fluid_histo = max(lambdaPmax_in_fluid_histo,vpImin_local / (distance_max_local / (NGLLX - 1)))
    endif

    lambdaPImin = min(lambdaPImin,vpImin_local / (distance_max_local / (NGLLX - 1)))
    lambdaPImax = max(lambdaPImax,vpImin_local / (distance_max_local / (NGLLX - 1)))

    if(cpIIloc > 1.d-20) then
      lambdaPIImin = min(lambdaPIImin,vpIImin_local / (distance_max_local / (NGLLX - 1)))
      lambdaPIImax = max(lambdaPIImax,vpIImin_local / (distance_max_local / (NGLLX - 1)))
    endif

  enddo

#ifdef USE_MPI
  call MPI_ALLREDUCE (vpImin, vpImin_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (vpImax, vpImax_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (vpIImin, vpIImin_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (vpIImax, vpIImax_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (vsmin, vsmin_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (vsmax, vsmax_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (densmin, densmin_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (densmax, densmax_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (distance_min, distance_min_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (distance_max, distance_max_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (courant_stability_number_max, courant_stability_max_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (lambdaPImin, lambdaPImin_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (lambdaPImax, lambdaPImax_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (lambdaPIImin, lambdaPIImin_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (lambdaPIImax, lambdaPIImax_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (lambdaSmin, lambdaSmin_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (lambdaSmax, lambdaSmax_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (lambdaPmin_in_fluid_histo, lambdaPmin_in_fluid_histo_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (lambdaPmax_in_fluid_histo, lambdaPmax_in_fluid_histo_glob, 1, MPI_DOUBLE_PRECISION, &
                    MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (any_elastic, any_elastic_glob, 1, MPI_LOGICAL, &
                    MPI_LOR, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (any_poroelastic, any_poroelastic_glob, 1, MPI_LOGICAL, &
                    MPI_LOR, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (any_fluid_histo, any_fluid_histo_glob, 1, MPI_LOGICAL, &
                    MPI_LOR, MPI_COMM_WORLD, ier)
  vpImin = vpImin_glob
  vpImax = vpImax_glob
  vpIImin = vpIImin_glob
  vpIImax = vpIImax_glob
  vsmin = vsmin_glob
  vsmax = vsmax_glob
  densmin = densmin_glob
  densmax = densmax_glob
  distance_min = distance_min_glob
  distance_max = distance_max_glob
  courant_stability_number_max = courant_stability_max_glob
  lambdaPImin = lambdaPImin_glob
  lambdaPImax = lambdaPImax_glob
  lambdaPIImin = lambdaPIImin_glob
  lambdaPIImax = lambdaPIImax_glob
  lambdaSmin = lambdaSmin_glob
  lambdaSmax = lambdaSmax_glob
  lambdaPmin_in_fluid_histo = lambdaPmin_in_fluid_histo_glob
  lambdaPmax_in_fluid_histo = lambdaPmax_in_fluid_histo_glob
#else
  any_elastic_glob = any_elastic
  any_poroelastic_glob = any_poroelastic
  any_fluid_histo_glob = any_fluid_histo
#endif

  if ( myrank == 0 ) then
    if(.not. all_anisotropic) then
      write(IOUT,*)
      write(IOUT,*) '********'
      write(IOUT,*) 'Model: P (or PI) velocity min,max = ',vpImin,vpImax
      write(IOUT,*) 'Model: PII velocity min,max = ',vpIImin,vpIImax
      write(IOUT,*) 'Model: S velocity min,max = ',vsmin,vsmax
      write(IOUT,*) 'Model: density min,max = ',densmin,densmax
      write(IOUT,*) '********'
      write(IOUT,*)

      write(IOUT,*)
      write(IOUT,*) '*********************************************'
      write(IOUT,*) '*** Verification of simulation parameters ***'
      write(IOUT,*) '*********************************************'
      write(IOUT,*)
      write(IOUT,*) '*** Max grid size = ',distance_max
      write(IOUT,*) '*** Min grid size = ',distance_min
      write(IOUT,*) '*** Max/min ratio = ',distance_max / distance_min
      write(IOUT,*)
      write(IOUT,*) '*** Max CFL stability condition of the time scheme &
                         &based on P wave velocity (must be below about 0.50 or so) = ',courant_stability_number_max
      write(IOUT,*)
    endif

    create_wavelength_histogram = .false.

! only if time source is not a Dirac or Heaviside (otherwise maximum frequency of spectrum undefined)
! and if source is not an initial field, for the same reason
    if(.not. initialfield) then
      f0max = -HUGEVAL

      do i = 1,NSOURCES

        ! excludes Dirac and Heaviside sources
        if(time_function_type(i) /= 4 .and. time_function_type(i) /= 5) then

          ! sets min/max frequency
          if(f0(i) > f0max) f0max = f0(i)

          if( i == NSOURCES ) then
            write(IOUT,*) '----'
            write(IOUT,*) ' Nb pts / lambdaPI_fmax min = ',lambdaPImin/(2.5d0*f0max)
            write(IOUT,*) ' Nb pts / lambdaPI_fmax max = ',lambdaPImax/(2.5d0*f0max)
            write(IOUT,*) '----'
            write(IOUT,*) ' Nb pts / lambdaPII_fmax min = ',lambdaPIImin/(2.5d0*f0max)
            write(IOUT,*) ' Nb pts / lambdaPII_fmax max = ',lambdaPIImax/(2.5d0*f0max)
            write(IOUT,*) '----'
            write(IOUT,*) ' Nb pts / lambdaS_fmax min = ',lambdaSmin/(2.5d0*f0max)
            write(IOUT,*) ' Nb pts / lambdaS_fmax max = ',lambdaSmax/(2.5d0*f0max)
            write(IOUT,*) '----'

! for histogram
            lambdaPmin_in_fluid_histo = lambdaPmin_in_fluid_histo/(2.5d0*f0max)
            lambdaPmax_in_fluid_histo = lambdaPmax_in_fluid_histo/(2.5d0*f0max)

            lambdaSmin_histo = lambdaSmin/(2.5d0*f0max)
            lambdaSmax_histo = lambdaSmax/(2.5d0*f0max)

            create_wavelength_histogram = .true.

          endif

        endif
      enddo
    endif
  endif

#ifdef USE_MPI
  call MPI_BCAST(lambdaPmin_in_fluid_histo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(lambdaPmax_in_fluid_histo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(lambdaSmin_histo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(lambdaSmax_histo,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(f0max,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ier)
  call MPI_BCAST(create_wavelength_histogram,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ier)
#endif

!!!!!!!!!!!!!!!! DK DK: added histogram of minimum number of points per wavelength

!! DK DK take into account the fact that there is no S velocity in the fluid
!! DK DK in this case, for the fluid, use the P wave data

  if(create_wavelength_histogram) then

! create statistics about mesh sampling (number of points per wavelength)

  nspec_counted_all_solid = 0
  nspec_counted_all_fluid = 0

! first pass is for S wave sampling in solid, second pass is for P wave sampling in fluid
  do ipass = 1,2

  nspec_counted = 0

  if(ipass == 1) then
    min_nb_of_points_per_wavelength = lambdaSmin_histo
    max_nb_of_points_per_wavelength = lambdaSmax_histo
! do not create this histogram if the model is entirely fluid
    if(.not. any_elastic_glob .and. .not. any_poroelastic_glob) cycle
  else
! do not create this histogram if the model is entirely solid
    if(.not. any_fluid_histo_glob) cycle
    min_nb_of_points_per_wavelength = lambdaPmin_in_fluid_histo
    max_nb_of_points_per_wavelength = lambdaPmax_in_fluid_histo
  endif

! when the grid is regular and the medium is homogeneous, the minimum and the maximum are equal
! and thus we cannot create an histogram; in such a case, let us artificially create a non-empty range
  if(abs(max_nb_of_points_per_wavelength - min_nb_of_points_per_wavelength) < 1.d-10) then
    min_nb_of_points_per_wavelength = min_nb_of_points_per_wavelength * 0.99d0
    max_nb_of_points_per_wavelength = max_nb_of_points_per_wavelength * 1.01d0
  endif

! erase histogram of wavelength
  classes_wavelength(:) = 0

! loop on all the elements
  do ispec = 1,nspec

    material = kmato(ispec)

    if(poroelastic(ispec)) then

      ! poroelastic material

      phi = porosity(material)
      tort = tortuosity(material)
      perm = permeability(1,material)
      ! solid properties
      mu_s = poroelastcoef(2,1,material)
      kappa_s = poroelastcoef(3,1,material) - FOUR_THIRDS*mu_s
      denst_s = density(1,material)
      denst = denst_s
      ! fluid properties
      kappa_f = poroelastcoef(1,2,material)
      denst_f = density(2,material)
      eta_f = poroelastcoef(2,2,material)
      ! frame properties
      mu_fr = poroelastcoef(2,3,material)
      kappa_fr = poroelastcoef(3,3,material) - FOUR_THIRDS*mu_fr
      ! Biot coefficients for the input phi
      D_biot = kappa_s*(1.d0 + phi*(kappa_s/kappa_f - 1.d0))
      H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr + FOUR_THIRDS*mu_fr
      C_biot = kappa_s*(kappa_s - kappa_fr)/(D_biot - kappa_fr)
      M_biot = kappa_s*kappa_s/(D_biot - kappa_fr)

      call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
             tort,denst_s,denst_f,eta_f,perm,f0(1),freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

      cpIloc = sqrt(cpIsquare)
      cpIIloc = sqrt(cpIIsquare)
      csloc = sqrt(cssquare)
    else
      mu = poroelastcoef(2,1,material)
      lambdaplus2mu  = poroelastcoef(3,1,material)
      denst = density(1,material)

      cpIloc = sqrt(lambdaplus2mu/denst)
      cpIIloc = 0.d0
      csloc = sqrt(mu/denst)
    endif

    vpImin_local = HUGEVAL
    vpIImin_local = HUGEVAL
    vsmin_local = HUGEVAL

    distance_min_local = HUGEVAL
    distance_max_local = -HUGEVAL

    do j=1,NGLLZ
      do i=1,NGLLX

!--- if heterogeneous formulation with external velocity model
        if(assign_external_model) then
          cpIloc = vpext(i,j,ispec)
          csloc = vsext(i,j,ispec)
        endif

        vpImin_local = min(vpImin_local,cpIloc)
        vpIImin_local = min(vpIImin_local,cpIIloc)
        vsmin_local = min(vsmin_local,csloc)

      enddo
    enddo

! compute minimum and maximum size of edges of this grid cell
    distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
               (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

    distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

    distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

    distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
               (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

    distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
    distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

  if(ipass == 1) then

! in first pass, only solid regions, thus ignore fluid regions with Vs = 0
    if(vsmin_local > 1.d-20) then
      nb_of_points_per_wavelength = vsmin_local / (distance_max_local / (NGLLX - 1))

      nspec_counted = nspec_counted + 1

      nb_of_points_per_wavelength = nb_of_points_per_wavelength/(2.5d0*f0max)

! store number of points per wavelength in histogram
      iclass = int((nb_of_points_per_wavelength - min_nb_of_points_per_wavelength) / &
                   (max_nb_of_points_per_wavelength - min_nb_of_points_per_wavelength) * dble(NCLASSES))
      if(iclass < 0) iclass = 0
      if(iclass > NCLASSES-1) iclass = NCLASSES-1
      classes_wavelength(iclass) = classes_wavelength(iclass) + 1

    endif

  else

! in second pass, only fluid regions, thus ignore solid regions with Vs > 0
    if(abs(vsmin_local) < 1.d-20) then
      if(vpIImin_local <= ZERO) then
        nb_of_points_per_wavelength = vpImin_local / (distance_max_local / (NGLLX - 1))
      else
        nb_of_points_per_wavelength = min(vpImin_local,vpIImin_local) / (distance_max_local / (NGLLX - 1))
      endif

      nspec_counted = nspec_counted + 1

      nb_of_points_per_wavelength = nb_of_points_per_wavelength/(2.5d0*f0max)

! store number of points per wavelength in histogram
      iclass = int((nb_of_points_per_wavelength - min_nb_of_points_per_wavelength) / &
                   (max_nb_of_points_per_wavelength - min_nb_of_points_per_wavelength) * dble(NCLASSES))
      if(iclass < 0) iclass = 0
      if(iclass > NCLASSES-1) iclass = NCLASSES-1
      classes_wavelength(iclass) = classes_wavelength(iclass) + 1
    endif

  endif

  enddo

#ifdef USE_MPI
  call MPI_REDUCE(classes_wavelength, classes_wavelength_all, NCLASSES, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
#else
  classes_wavelength_all(:) = classes_wavelength(:)
#endif

#ifdef USE_MPI
  call MPI_REDUCE(nspec, nspec_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
  call MPI_REDUCE(nspec_counted, nspec_counted_all, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ier)
#else
  nspec_counted_all = nspec_counted
  nspec_all = nspec
#endif

  if(ipass == 1) then
    nspec_counted_all_solid = nspec_counted_all
  else
    nspec_counted_all_fluid = nspec_counted_all
  endif

! create histogram of wavelength and save in Gnuplot file
  if (myrank == 0) then

  write(IOUT,*)
  write(IOUT,*) '-----------------------------------------'
  write(IOUT,*)
  if(ipass == 1) then
    write(IOUT,*) 'histogram of min number of points per S wavelength in solid regions:'
    write(IOUT,*)
    write(IOUT,*) 'there are ',nspec_counted_all,' elements out of ',nspec_all,' in solid regions'
  else
    write(IOUT,*) 'histogram of min number of points per P wavelength in fluid regions:'
    write(IOUT,*)
    write(IOUT,*) 'there are ',nspec_counted_all,' elements out of ',nspec_all,' in fluid regions'
  endif
  write(IOUT,*) '  (i.e., ',sngl(100.d0*nspec_counted_all/dble(nspec_all)),'% of the total)'
  write(IOUT,*)
  write(IOUT,*) '(too small = poor resolution of calculations -'
  write(IOUT,*) ' too big = wasting memory and CPU time)'
  write(IOUT,*) '(threshold value is around 4.5 points per S wavelength'
  write(IOUT,*) ' in elastic regions and 5.5 per P wavelength in fluid regions):'
  write(IOUT,*)

  total_percent = 0.
  scaling_factor = max_nb_of_points_per_wavelength - min_nb_of_points_per_wavelength

  if(ipass == 1) then
    open(unit=14,file='points_per_wavelength_histogram_S_in_solid.txt',status='unknown')
    scaling_factor_S = scaling_factor
  else
    open(unit=14,file='points_per_wavelength_histogram_P_in_fluid.txt',status='unknown')
    scaling_factor_P = scaling_factor
  endif
  do iclass = 0,NCLASSES-1
    current_percent = 100.*dble(classes_wavelength_all(iclass))/dble(nspec_counted_all)
    total_percent = total_percent + current_percent
    write(IOUT,*) sngl(min_nb_of_points_per_wavelength + scaling_factor*iclass/dble(NCLASSES)),' - ', &
                  sngl(min_nb_of_points_per_wavelength + scaling_factor*(iclass+1)/dble(NCLASSES)),classes_wavelength_all(iclass), &
                  ' ',sngl(current_percent),' %'
    write(14,*) 0.5*(sngl(min_nb_of_points_per_wavelength + scaling_factor*iclass/dble(NCLASSES)) + &
                     sngl(min_nb_of_points_per_wavelength + scaling_factor*(iclass+1)/dble(NCLASSES))),' ',sngl(current_percent)
  enddo
  close(14)

  if(total_percent < 99.9d0 .or. total_percent > 100.1d0) then
    write(IOUT,*) 'total percentage = ',total_percent,' %'
    stop 'total percentage should be 100%'
  else
    write(IOUT,*)
    write(IOUT,*) 'total percentage = ',total_percent,' %'
  endif

  endif ! of if myrank == 0

  enddo ! end of the two passes on S wavelength data and P wavelength data

! create script for Gnuplot histogram file
  if (myrank == 0) then

  open(unit=14,file='plot_points_per_wavelength_histogram.gnu',status='unknown')
  write(14,*) 'set term wxt'

  if(nspec_counted_all_solid > 0) then
    write(14,*) '#set term gif'
    write(14,*) '#set output "points_per_wavelength_histogram_S_in_solid.gif"'
    write(14,*)
    write(14,*) 'set boxwidth ',real(scaling_factor_S/NCLASSES)
    write(14,*) 'set xlabel "Range of min number of points per S wavelength in solid"'
    write(14,*) 'set ylabel "Percentage of elements (%)"'
    write(14,*) 'plot "points_per_wavelength_histogram_S_in_solid.txt" with boxes'
    write(14,*) 'pause -1 "hit any key..."'
  endif

  if(nspec_counted_all_fluid > 0) then
    write(14,*) '#set term gif'
    write(14,*) '#set output "points_per_wavelength_histogram_P_in_fluid.gif"'
    write(14,*)
    write(14,*) 'set boxwidth ',real(scaling_factor_P/NCLASSES)
    write(14,*) 'set xlabel "Range of min number of points per P wavelength in fluid"'
    write(14,*) 'set ylabel "Percentage of elements (%)"'
    write(14,*) 'plot "points_per_wavelength_histogram_P_in_fluid.txt" with boxes'
    write(14,*) 'pause -1 "hit any key..."'
  endif

  close(14)

  write(IOUT,*)
  write(IOUT,*)
  write(IOUT,*) 'total number of elements in fluid and solid regions = ',nspec_all
  write(IOUT,*)

  endif ! of if myrank == 0

  endif ! of if create_wavelength_histogram

!!!!!!!!!!!!!!!! DK DK: added histogram of minimum number of points per wavelength

!
!--------------------------------------------------------------------------------
!

! do not create the PostScript file below if the postscript flag is off
  if(.not. output_postscript_snapshot) return

! A4 or US letter paper
  if(US_LETTER) then
    usoffset = 1.75d0
    sizex = 27.94d0
    sizez = 21.59d0
  else
    usoffset = 0.d0
    sizex = 29.7d0
    sizez = 21.d0
  endif

! height of domain numbers in centimeters
  height = 0.25d0

! get minimum and maximum values of mesh coordinates
  xmin = minval(coord(1,:))
  zmin = minval(coord(2,:))
  xmax = maxval(coord(1,:))
  zmax = maxval(coord(2,:))

#ifdef USE_MPI
  call MPI_ALLREDUCE (xmin, xmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (xmax, xmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (zmin, zmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (zmax, zmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  xmin = xmin_glob
  xmax = xmax_glob
  zmin = zmin_glob
  zmax = zmax_glob

#endif

! ratio of physical page size/size of the domain meshed
  ratio_page = min(rpercentz*sizez/(zmax-zmin),rpercentx*sizex/(xmax-xmin)) / 100.d0


  if (myrank == 0) then

    write(IOUT,*)
    write(IOUT,*) 'Creating PostScript file with stability condition'

!
!---- open PostScript file
!
    open(unit=24,file='OUTPUT_FILES/mesh_stability.ps',status='unknown')

!
!---- write PostScript header
!
    write(24,10) simulation_title
    write(24,*) '/CM {28.5 mul} def'
    write(24,*) '/LR {rlineto} def'
    write(24,*) '/LT {lineto} def'
    write(24,*) '/L {lineto} def'
    write(24,*) '/MR {rmoveto} def'
    write(24,*) '/MV {moveto} def'
    write(24,*) '/M {moveto} def'
    write(24,*) '/ST {stroke} def'
    write(24,*) '/CP {closepath} def'
    write(24,*) '/RG {setrgbcolor} def'
    write(24,*) '/GF {gsave fill grestore} def'
    write(24,*) '% different useful symbols'
    write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
    write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
    write(24,*) 'CP fill} def'
    write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
    write(24,*) 'CP fill} def'
    write(24,*) '/Cross {gsave 0.05 CM setlinewidth'
    write(24,*) 'gsave 3 3 MR -6. -6. LR ST grestore'
    write(24,*) 'gsave 3 -3 MR -6. 6. LR ST grestore'
    write(24,*) '0.01 CM setlinewidth} def'
    write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
    write(24,*) '/Diamond {gsave 0.05 CM setlinewidth 0 4.2 MR'
    write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
    write(24,*) 'grestore 0.01 CM setlinewidth} def'
    write(24,*) '%'
    write(24,*) '% macro to draw the contour of the elements'
    write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
    write(24,*) '%'
    write(24,*) '.01 CM setlinewidth'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.35 CM scalefont setfont'
    write(24,*) '%'
    write(24,*) '/vshift ',-height/2,' CM def'
    write(24,*) '/Rshow { currentpoint stroke MV'
    write(24,*) 'dup stringwidth pop neg vshift MR show } def'
    write(24,*) '/Cshow { currentpoint stroke MV'
    write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
    write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM scalefont setfont} def'
    write(24,*) '%'
    write(24,*) 'gsave newpath 90 rotate'
    write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
    write(24,*) '%'

    !
    !--- write captions of PostScript figure
    !
    write(24,*) '0 setgray'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.5 CM scalefont setfont'

    write(24,*) '%'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.6 CM scalefont setfont'
    write(24,*) '.4 .9 .9 setrgbcolor'
    write(24,*) '11 CM 1.1 CM MV'
    write(24,*) '(X axis) show'
    write(24,*) '%'
    write(24,*) '1.4 CM 9.5 CM MV'
    write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
    write(24,*) '(Z axis) show'
    write(24,*) 'grestore'
    write(24,*) '%'
    write(24,*) '/Times-Roman findfont'
    write(24,*) '.7 CM scalefont setfont'
    write(24,*) '.8 0 .8 setrgbcolor'
    write(24,*) '24.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(Mesh stability condition \(red = bad\)) show'
    write(24,*) 'grestore'
    write(24,*) '25.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(',simulation_title,') show'
    write(24,*) 'grestore'
    write(24,*) '26.45 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,*) '(2D Spectral Element Method) show'
    write(24,*) 'grestore'

    write(24,*) '%'
    write(24,*) '1 1 scale'
    write(24,*) '%'

    !
    !---- draw the spectral element mesh
    !
    write(24,*) '%'
    write(24,*) '% spectral element mesh'
    write(24,*) '%'
    write(24,*) '0 setgray'

    num_ispec = 0
  endif

  do ispec = 1, nspec
    if ( myrank == 0 ) then
      num_ispec = num_ispec + 1
      write(24,*) '% elem ',num_ispec
    endif

    do i=1,pointsdisp
      do j=1,pointsdisp
        xinterp(i,j) = 0.d0
        zinterp(i,j) = 0.d0
        do in = 1,ngnod
          nnum = knods(in,ispec)
          xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
          zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
        enddo
      enddo
    enddo

    is = 1
    ir = 1
    x1 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
    z1 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
    x1 = x1 * centim
    z1 = z1 * centim
    if ( myrank == 0 ) then
      write(24,*) 'mark'
      write(24,681) x1,z1
    else
      coorg_send(1,(ispec-1)*5+1) = x1
      coorg_send(2,(ispec-1)*5+1) = z1
    endif

    ! draw straight lines if elements have 4 nodes

    ir=pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    if ( myrank == 0 ) then
      write(24,681) x2,z2
    else
      coorg_send(1,(ispec-1)*5+2) = x2
      coorg_send(2,(ispec-1)*5+2) = z2
    endif

    ir=pointsdisp
    is=pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    if ( myrank == 0 ) then
      write(24,681) x2,z2
    else
      coorg_send(1,(ispec-1)*5+3) = x2
      coorg_send(2,(ispec-1)*5+3) = z2
    endif

    is=pointsdisp
    ir=1
    x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    if ( myrank == 0 ) then
      write(24,681) x2,z2
    else
      coorg_send(1,(ispec-1)*5+4) = x2
      coorg_send(2,(ispec-1)*5+4) = z2
    endif

    ir=1
    is=2
    x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    if ( myrank == 0 ) then
      write(24,681) x2,z2
      write(24,*) 'CO'
    else
      coorg_send(1,(ispec-1)*5+5) = x2
      coorg_send(2,(ispec-1)*5+5) = z2
    endif

    material = kmato(ispec)

    if(poroelastic(ispec)) then

      ! poroelastic material

      phi=porosity(material)
      tort=tortuosity(material)
      perm=permeability(1,material)
      ! solid properties
      mu_s = poroelastcoef(2,1,material)
      kappa_s = poroelastcoef(3,1,material) - FOUR_THIRDS*mu_s
      denst_s = density(1,material)
      denst = denst_s
      ! fluid properties
      kappa_f = poroelastcoef(1,2,material)
      denst_f = density(2,material)
      eta_f = poroelastcoef(2,2,material)
      ! frame properties
      mu_fr = poroelastcoef(2,3,material)
      kappa_fr = poroelastcoef(3,3,material) - FOUR_THIRDS*mu_fr
      ! Biot coefficients for the input phi
      D_biot = kappa_s*(1.d0 + phi*(kappa_s/kappa_f - 1.d0))
      H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr + FOUR_THIRDS*mu_fr
      C_biot = kappa_s*(kappa_s - kappa_fr)/(D_biot - kappa_fr)
      M_biot = kappa_s*kappa_s/(D_biot - kappa_fr)

      call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
           tort,denst_s,denst_f,eta_f,perm,f0(1),freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

      cpIloc = sqrt(cpIsquare)
    else
      lambdaplus2mu  = poroelastcoef(3,1,material)
      denst = density(1,material)

      cpIloc = sqrt(lambdaplus2mu/denst)
    endif

    vpImax_local = -HUGEVAL

    distance_min_local = HUGEVAL
    distance_max_local = -HUGEVAL

    do j=1,NGLLZ
      do i=1,NGLLX

        !--- if heterogeneous formulation with external velocity model
        if(assign_external_model) then
          cpIloc = vpext(i,j,ispec)
          denst = rhoext(i,j,ispec)
        endif

        vpImax_local = max(vpImax_local,cpIloc)

      enddo
    enddo

! compute minimum and maximum size of edges of this grid cell
    distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
             (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

    distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
             (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

    distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
             (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

    distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
             (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

    distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
    distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

    distance_min = min(distance_min,distance_min_local)
    distance_max = max(distance_max,distance_max_local)

    courant_stability_number = vpImax_local * deltat / (distance_min_local * percent_GLL(NGLLX))

! display bad elements that are above the threshold
    if(courant_stability_number >= THRESHOLD_POSTSCRIPT * courant_stability_number_max) then
      if ( myrank == 0 ) then
        write(24,*) '1 0 0 RG GF 0 setgray ST'
      else
        RGB_send(ispec) = 1
      endif
    else
! do not color the elements if below the threshold
      if ( myrank == 0 ) then
        write(24,*) 'ST'
      else
        RGB_send(ispec) = 0
      endif
    endif

  enddo ! end of loop on all the spectral elements

#ifdef USE_MPI
  if (myrank == 0 ) then

    do iproc = 1, nproc-1
      call MPI_RECV (nspec_recv, 1, MPI_INTEGER, &
              iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      allocate(coorg_recv(2,nspec_recv*5))
      allocate(RGB_recv(nspec_recv))
      call MPI_RECV (coorg_recv(1,1), nspec_recv*5*2, MPI_DOUBLE_PRECISION, &
              iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
      call MPI_RECV (RGB_recv(1), nspec_recv, MPI_INTEGER, &
              iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

      do ispec = 1, nspec_recv
        num_ispec = num_ispec + 1
        write(24,*) '% elem ',num_ispec
        write(24,*) 'mark'
        write(24,681) coorg_recv(1,(ispec-1)*5+1), coorg_recv(2,(ispec-1)*5+1)
        write(24,681) coorg_recv(1,(ispec-1)*5+2), coorg_recv(2,(ispec-1)*5+2)
        write(24,681) coorg_recv(1,(ispec-1)*5+3), coorg_recv(2,(ispec-1)*5+3)
        write(24,681) coorg_recv(1,(ispec-1)*5+4), coorg_recv(2,(ispec-1)*5+4)
        write(24,681) coorg_recv(1,(ispec-1)*5+5), coorg_recv(2,(ispec-1)*5+5)
        write(24,*) 'CO'
        if ( RGB_recv(ispec)  == 1) then
          write(24,*) '1 0 0 RG GF 0 setgray ST'
        else
          write(24,*) 'ST'
        endif
      enddo
      deallocate(coorg_recv)
      deallocate(RGB_recv)

    enddo

  else
    call MPI_SEND (nspec, 1, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
    call MPI_SEND (coorg_send, nspec*5*2, MPI_DOUBLE_PRECISION, 0, 42, MPI_COMM_WORLD, ier)
    call MPI_SEND (RGB_send, nspec, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
  endif

#endif

  if ( myrank == 0 ) then
    write(24,*) '%'
    write(24,*) 'grestore'
    write(24,*) 'showpage'

    close(24)

    write(IOUT,*) 'End of creation of PostScript file with stability condition'
  endif

!
!--------------------------------------------------------------------------------
!

  if (myrank == 0) then

    write(IOUT,*)
    write(IOUT,*) 'Creating PostScript file with mesh dispersion'

!
!---- open PostScript file
!
    if(any_elastic_glob .or. any_poroelastic_glob) then
      open(unit=24,file='OUTPUT_FILES/mesh_S_wave_dispersion.ps',status='unknown')
    else
      open(unit=24,file='OUTPUT_FILES/mesh_P_wave_dispersion.ps',status='unknown')
    endif

!
!---- write PostScript header
!
  write(24,10) simulation_title
  write(24,*) '/CM {28.5 mul} def'
  write(24,*) '/LR {rlineto} def'
  write(24,*) '/LT {lineto} def'
  write(24,*) '/L {lineto} def'
  write(24,*) '/MR {rmoveto} def'
  write(24,*) '/MV {moveto} def'
  write(24,*) '/M {moveto} def'
  write(24,*) '/ST {stroke} def'
  write(24,*) '/CP {closepath} def'
  write(24,*) '/RG {setrgbcolor} def'
  write(24,*) '/GF {gsave fill grestore} def'
  write(24,*) '% different useful symbols'
  write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
  write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/Cross {gsave 0.05 CM setlinewidth'
  write(24,*) 'gsave 3 3 MR -6. -6. LR ST grestore'
  write(24,*) 'gsave 3 -3 MR -6. 6. LR ST grestore'
  write(24,*) '0.01 CM setlinewidth} def'
  write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
  write(24,*) '/Diamond {gsave 0.05 CM setlinewidth 0 4.2 MR'
  write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
  write(24,*) 'grestore 0.01 CM setlinewidth} def'
  write(24,*) '%'
  write(24,*) '% macro to draw the contour of the elements'
  write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
  write(24,*) '%'
  write(24,*) '.01 CM setlinewidth'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.35 CM scalefont setfont'
  write(24,*) '%'
  write(24,*) '/vshift ',-height/2,' CM def'
  write(24,*) '/Rshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop neg vshift MR show } def'
  write(24,*) '/Cshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
  write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM scalefont setfont} def'
  write(24,*) '%'
  write(24,*) 'gsave newpath 90 rotate'
  write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
  write(24,*) '%'

!
!--- write captions of PostScript figure
!
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.5 CM scalefont setfont'

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM scalefont setfont'
  write(24,*) '.4 .9 .9 setrgbcolor'
  write(24,*) '11 CM 1.1 CM MV'
  write(24,*) '(X axis) show'
  write(24,*) '%'
  write(24,*) '1.4 CM 9.5 CM MV'
  write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
  write(24,*) '(Z axis) show'
  write(24,*) 'grestore'
  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.7 CM scalefont setfont'
  write(24,*) '.8 0 .8 setrgbcolor'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  if(any_elastic_glob) then
    write(24,*) '(Mesh elastic S-wave dispersion \(red = good, blue = bad\)) show'
  else
    write(24,*) '(Mesh acoustic P-wave dispersion \(red = good, blue = bad\)) show'
  endif
  write(24,*) 'grestore'
  write(24,*) '25.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(',simulation_title,') show'
  write(24,*) 'grestore'
  write(24,*) '26.45 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(2D Spectral Element Method) show'
  write(24,*) 'grestore'

  write(24,*) '%'
  write(24,*) '1 1 scale'
  write(24,*) '%'

!
!---- draw the spectral element mesh
!
  write(24,*) '%'
  write(24,*) '% spectral element mesh'
  write(24,*) '%'
  write(24,*) '0 setgray'

  num_ispec = 0
  endif

  do ispec = 1, nspec
     if ( myrank == 0 ) then
        num_ispec = num_ispec + 1
        write(24,*) '% elem ',num_ispec
     endif

  do i=1,pointsdisp
  do j=1,pointsdisp
  xinterp(i,j) = 0.d0
  zinterp(i,j) = 0.d0
  do in = 1,ngnod
    nnum = knods(in,ispec)
      xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
      zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
  enddo
  enddo
  enddo

  is = 1
  ir = 1
  x1 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z1 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  if ( myrank == 0 ) then
     write(24,*) 'mark'
     write(24,681) x1,z1
  else
     coorg_send(1,(ispec-1)*5+1) = x1
     coorg_send(2,(ispec-1)*5+1) = z1
  endif

! draw straight lines if elements have 4 nodes

  ir=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     coorg_send(1,(ispec-1)*5+2) = x2
     coorg_send(2,(ispec-1)*5+2) = z2
  endif

  ir=pointsdisp
  is=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     coorg_send(1,(ispec-1)*5+3) = x2
     coorg_send(2,(ispec-1)*5+3) = z2
  endif

  is=pointsdisp
  ir=1
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     coorg_send(1,(ispec-1)*5+4) = x2
     coorg_send(2,(ispec-1)*5+4) = z2
  endif

  ir=1
  is=2
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
     write(24,*) 'CO'
  else
     coorg_send(1,(ispec-1)*5+5) = x2
     coorg_send(2,(ispec-1)*5+5) = z2
  endif

    material = kmato(ispec)

   if(poroelastic(ispec)) then
    phi = porosity(material)
    tort = tortuosity(material)
    perm = permeability(1,material)
! solid properties
    mu_s = poroelastcoef(2,1,material)
    kappa_s = poroelastcoef(3,1,material) - FOUR_THIRDS*mu_s
    denst_s = density(1,material)
    denst = denst_s
! fluid properties
    kappa_f = poroelastcoef(1,2,material)
    denst_f = density(2,material)
    eta_f = poroelastcoef(2,2,material)
! frame properties
    mu_fr = poroelastcoef(2,3,material)
    kappa_fr = poroelastcoef(3,3,material) - FOUR_THIRDS*mu_fr
! Biot coefficients for the input phi
      D_biot = kappa_s*(1.d0 + phi*(kappa_s/kappa_f - 1.d0))
      H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr + FOUR_THIRDS*mu_fr
      C_biot = kappa_s*(kappa_s - kappa_fr)/(D_biot - kappa_fr)
      M_biot = kappa_s*kappa_s/(D_biot - kappa_fr)

    call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
             tort,denst_s,denst_f,eta_f,perm,f0(1),freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

    cpIloc = sqrt(cpIsquare)
    csloc = sqrt(cssquare)
   else
    mu = poroelastcoef(2,1,material)
    lambdaplus2mu  = poroelastcoef(3,1,material)
    denst = density(1,material)

    cpIloc = sqrt(lambdaplus2mu/denst)
    csloc = sqrt(mu/denst)
   endif

  vpImax_local = -HUGEVAL
  vpImin_local = HUGEVAL
  vsmin_local = HUGEVAL

  distance_min_local = HUGEVAL
  distance_max_local = -HUGEVAL

  do j=1,NGLLZ
    do i=1,NGLLX

!--- if heterogeneous formulation with external velocity model
    if(assign_external_model) then
      cpIloc = vpext(i,j,ispec)
      csloc = vsext(i,j,ispec)
      denst = rhoext(i,j,ispec)
    endif

    vpImax_local = max(vpImax_local,cpIloc)
    vpImin_local = min(vpImin_local,cpIloc)
    vsmin_local = min(vsmin_local,csloc)

    enddo
  enddo

! compute minimum and maximum size of edges of this grid cell
  distance_1 = sqrt((coord(1,ibool(1,1,ispec)) - coord(1,ibool(NGLLX,1,ispec)))**2 + &
               (coord(2,ibool(1,1,ispec)) - coord(2,ibool(NGLLX,1,ispec)))**2)

  distance_2 = sqrt((coord(1,ibool(NGLLX,1,ispec)) - coord(1,ibool(NGLLX,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,1,ispec)) - coord(2,ibool(NGLLX,NGLLZ,ispec)))**2)

  distance_3 = sqrt((coord(1,ibool(NGLLX,NGLLZ,ispec)) - coord(1,ibool(1,NGLLZ,ispec)))**2 + &
               (coord(2,ibool(NGLLX,NGLLZ,ispec)) - coord(2,ibool(1,NGLLZ,ispec)))**2)

  distance_4 = sqrt((coord(1,ibool(1,NGLLZ,ispec)) - coord(1,ibool(1,1,ispec)))**2 + &
               (coord(2,ibool(1,NGLLZ,ispec)) - coord(2,ibool(1,1,ispec)))**2)

  distance_min_local = min(distance_1,distance_2,distance_3,distance_4)
  distance_max_local = max(distance_1,distance_2,distance_3,distance_4)

  distance_min = min(distance_min,distance_min_local)
  distance_max = max(distance_max,distance_max_local)

! display mesh dispersion for S waves if there is at least one elastic element in the mesh
  if(any_elastic_glob .or. any_poroelastic_glob) then

! ignore fluid regions with Vs = 0
  if(vsmin_local > 1.d-20) then

    lambdaS_local = vsmin_local / (distance_max_local / (NGLLX - 1))

! display very good elements that are above the threshold in red
    if(lambdaS_local >= THRESHOLD_POSTSCRIPT * lambdaSmax) then
       if ( myrank == 0 ) then
          write(24,*) '1 0 0 RG GF 0 setgray ST'
       else
          RGB_send(ispec) = 1
       endif

! display bad elements that are below the threshold in blue
    else if(lambdaS_local <= (1. + (1. - THRESHOLD_POSTSCRIPT)) * lambdaSmin) then
       if ( myrank == 0 ) then
          write(24,*) '0 0 1 RG GF 0 setgray ST'
       else
          RGB_send(ispec) = 3
       endif

    else
! do not color the elements if not close to the threshold
       if ( myrank == 0 ) then
          write(24,*) 'ST'
       else
          RGB_send(ispec) = 0
       endif
    endif

  else
! do not color the elements if S-wave velocity undefined
     if ( myrank == 0 ) then
        write(24,*) 'ST'
     else
        RGB_send(ispec) = 0
     endif
  endif

! display mesh dispersion for P waves if there is no elastic element in the mesh
  else

    lambdaPI_local = vpImin_local / (distance_max_local / (NGLLX - 1))

! display very good elements that are above the threshold in red
    if(lambdaPI_local >= THRESHOLD_POSTSCRIPT * lambdaPImax) then
       if ( myrank == 0 ) then
          write(24,*) '1 0 0 RG GF 0 setgray ST'
       else
          RGB_send(ispec) = 1
       endif

! display bad elements that are below the threshold in blue
    else if(lambdaPI_local <= (1. + (1. - THRESHOLD_POSTSCRIPT)) * lambdaPImin) then
       if ( myrank == 0 ) then
          write(24,*) '0 0 1 RG GF 0 setgray ST'
       else
          RGB_send(ispec) = 3
       endif

    else
! do not color the elements if not close to the threshold
       if ( myrank == 0 ) then
          write(24,*) 'ST'
       else
          RGB_send(ispec) = 0
       endif
    endif

  endif

  enddo ! end of loop on all the spectral elements

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        allocate(coorg_recv(2,nspec_recv*5))
        allocate(RGB_recv(nspec_recv))
        call MPI_RECV (coorg_recv(1,1), nspec_recv*5*2, MPI_DOUBLE_PRECISION, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        call MPI_RECV (RGB_recv(1), nspec_recv, MPI_INTEGER, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        do ispec = 1, nspec_recv
           num_ispec = num_ispec + 1
           write(24,*) '% elem ',num_ispec
           write(24,*) 'mark'
           write(24,681) coorg_recv(1,(ispec-1)*5+1), coorg_recv(2,(ispec-1)*5+1)
           write(24,681) coorg_recv(1,(ispec-1)*5+2), coorg_recv(2,(ispec-1)*5+2)
           write(24,681) coorg_recv(1,(ispec-1)*5+3), coorg_recv(2,(ispec-1)*5+3)
           write(24,681) coorg_recv(1,(ispec-1)*5+4), coorg_recv(2,(ispec-1)*5+4)
           write(24,681) coorg_recv(1,(ispec-1)*5+5), coorg_recv(2,(ispec-1)*5+5)
           write(24,*) 'CO'
           if ( RGB_recv(ispec)  == 1) then
              write(24,*) '1 0 0 RG GF 0 setgray ST'
           endif
           if ( RGB_recv(ispec)  == 3) then
              write(24,*) '0 0 1 RG GF 0 setgray ST'
           endif
           if ( RGB_recv(ispec)  == 0) then
              write(24,*) 'ST'
           endif

        enddo
        deallocate(coorg_recv)
        deallocate(RGB_recv)

     enddo

  else
     call MPI_SEND (nspec, 1, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
     call MPI_SEND (coorg_send, nspec*5*2, MPI_DOUBLE_PRECISION, 0, 42, MPI_COMM_WORLD, ier)
     call MPI_SEND (RGB_send, nspec, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)

  endif
#endif

  if ( myrank == 0 ) then
     write(24,*) '%'
     write(24,*) 'grestore'
     write(24,*) 'showpage'

     close(24)

     write(IOUT,*) 'End of creation of PostScript file with mesh dispersion'

  endif

!
!--------------------------------------------------------------------------------
!

  if (myrank == 0) then

    write(IOUT,*)
    write(IOUT,*) 'Creating PostScript file with velocity model'

!
!---- open PostScript file
!
  open(unit=24,file='OUTPUT_FILES/P_velocity_model.ps',status='unknown')

!
!---- write PostScript header
!
  write(24,10) simulation_title
  write(24,*) '/CM {28.5 mul} def'
  write(24,*) '/LR {rlineto} def'
  write(24,*) '/LT {lineto} def'
  write(24,*) '/L {lineto} def'
  write(24,*) '/MR {rmoveto} def'
  write(24,*) '/MV {moveto} def'
  write(24,*) '/M {moveto} def'
  write(24,*) '/ST {stroke} def'
  write(24,*) '/CP {closepath} def'
  write(24,*) '/RG {setrgbcolor} def'
  write(24,*) '/GF {gsave fill grestore} def'
  write(24,*) '% different useful symbols'
  write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
  write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/Cross {gsave 0.05 CM setlinewidth'
  write(24,*) 'gsave 3 3 MR -6. -6. LR ST grestore'
  write(24,*) 'gsave 3 -3 MR -6. 6. LR ST grestore'
  write(24,*) '0.01 CM setlinewidth} def'
  write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
  write(24,*) '/Diamond {gsave 0.05 CM setlinewidth 0 4.2 MR'
  write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
  write(24,*) 'grestore 0.01 CM setlinewidth} def'
  write(24,*) '%'
  write(24,*) '% macro to draw the contour of the elements'
  write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
  write(24,*) '%'
  write(24,*) '.01 CM setlinewidth'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.35 CM scalefont setfont'
  write(24,*) '%'
  write(24,*) '/vshift ',-height/2,' CM def'
  write(24,*) '/Rshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop neg vshift MR show } def'
  write(24,*) '/Cshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
  write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM scalefont setfont} def'
  write(24,*) '%'
  write(24,*) 'gsave newpath 90 rotate'
  write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
  write(24,*) '%'

!
!--- write captions of PostScript figure
!
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.5 CM scalefont setfont'

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM scalefont setfont'
  write(24,*) '.4 .9 .9 setrgbcolor'
  write(24,*) '11 CM 1.1 CM MV'
  write(24,*) '(X axis) show'
  write(24,*) '%'
  write(24,*) '1.4 CM 9.5 CM MV'
  write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
  write(24,*) '(Z axis) show'
  write(24,*) 'grestore'
  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.7 CM scalefont setfont'
  write(24,*) '.8 0 .8 setrgbcolor'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(P-velocity model \(dark = fast, light = slow\)) show'
  write(24,*) 'grestore'
  write(24,*) '25.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(',simulation_title,') show'
  write(24,*) 'grestore'
  write(24,*) '26.45 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(2D Spectral Element Method) show'
  write(24,*) 'grestore'

  write(24,*) '%'
  write(24,*) '1 1 scale'
  write(24,*) '%'

!
!---- draw the spectral element mesh
!
  write(24,*) '%'
  write(24,*) '% spectral element mesh'
  write(24,*) '%'
  write(24,*) '0 setgray'

  num_ispec = 0
endif

  do ispec = 1, UPPER_LIMIT_DISPLAY
     if ( myrank == 0 ) then
        num_ispec = num_ispec + 1
        write(24,*) '% elem ',num_ispec
     endif
  do i=1,pointsdisp
  do j=1,pointsdisp
  xinterp(i,j) = 0.d0
  zinterp(i,j) = 0.d0
  do in = 1,ngnod
    nnum = knods(in,ispec)
      xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
      zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
  enddo
  enddo
  enddo

  is = 1
  ir = 1
  x1 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z1 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  if ( myrank == 0 ) then
     write(24,*) 'mark'
     write(24,681) x1,z1
  else
     coorg_send(1,(ispec-1)*5+1) = x1
     coorg_send(2,(ispec-1)*5+1) = z1
  endif

! draw straight lines if elements have 4 nodes

  ir=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     coorg_send(1,(ispec-1)*5+2) = x2
     coorg_send(2,(ispec-1)*5+2) = z2
  endif

  ir=pointsdisp
  is=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     coorg_send(1,(ispec-1)*5+3) = x2
     coorg_send(2,(ispec-1)*5+3) = z2
  endif

  is=pointsdisp
  ir=1
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     coorg_send(1,(ispec-1)*5+4) = x2
     coorg_send(2,(ispec-1)*5+4) = z2
  endif

  ir=1
  is=2
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
     write(24,*) 'CO'
  else
     coorg_send(1,(ispec-1)*5+5) = x2
     coorg_send(2,(ispec-1)*5+5) = z2
  endif

  if((vpImax-vpImin)/vpImin > 0.02d0) then
  if(assign_external_model) then
! use lower-left corner
    x1 = (vpext(1,1,ispec)-vpImin) / (vpImax-vpImin)
  else
    material = kmato(ispec)
   if(poroelastic(ispec)) then
    phi = porosity(material)
    tort = tortuosity(material)
    perm = permeability(1,material)
!solid properties
    mu_s = poroelastcoef(2,1,material)
    kappa_s = poroelastcoef(3,1,material) - FOUR_THIRDS*mu_s
    denst_s = density(1,material)
!fluid properties
    kappa_f = poroelastcoef(1,2,material)
    denst_f = density(2,material)
    eta_f = poroelastcoef(2,2,material)
!frame properties
    mu_fr = poroelastcoef(2,3,material)
    kappa_fr = poroelastcoef(3,3,material) - FOUR_THIRDS*mu_fr
!Biot coefficients for the input phi
      D_biot = kappa_s*(1.d0 + phi*(kappa_s/kappa_f - 1.d0))
      H_biot = (kappa_s - kappa_fr)*(kappa_s - kappa_fr)/(D_biot - kappa_fr) + kappa_fr + FOUR_THIRDS*mu_fr
      C_biot = kappa_s*(kappa_s - kappa_fr)/(D_biot - kappa_fr)
      M_biot = kappa_s*kappa_s/(D_biot - kappa_fr)

    call get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mu_fr,phi, &
             tort,denst_s,denst_f,eta_f,perm,f0(1),freq0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

    cpIloc = sqrt(cpIsquare)
   else
    lambdaplus2mu  = poroelastcoef(3,1,material)
    denst = density(1,material)
    cpIloc = sqrt(lambdaplus2mu/denst)
   endif
    x1 = (cpIloc-vpImin)/(vpImax-vpImin)
  endif
  else
    x1 = 0.5d0
 endif

! rescale to avoid very dark gray levels
  x1 = x1*0.7 + 0.2
  if(x1 > 1.d0) x1=1.d0

! invert scale: white = vpmin, dark gray = vpmax
  x1 = 1.d0 - x1

! display P-velocity model using gray levels
  if ( myrank == 0 ) then
     write(24,*) sngl(x1),' setgray GF 0 setgray ST'
  else
     greyscale_send(ispec) = sngl(x1)
  endif
  enddo ! end of loop on all the spectral elements

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        allocate(coorg_recv(2,nspec_recv*5))
        allocate(greyscale_recv(nspec_recv))
        call MPI_RECV (coorg_recv(1,1), nspec_recv*5*2, MPI_DOUBLE_PRECISION, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        call MPI_RECV (greyscale_recv(1), nspec_recv, MPI_REAL, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        do ispec = 1, nspec_recv
           num_ispec = num_ispec + 1
           write(24,*) '% elem ',num_ispec
           write(24,*) 'mark'
           write(24,681) coorg_recv(1,(ispec-1)*5+1), coorg_recv(2,(ispec-1)*5+1)
           write(24,681) coorg_recv(1,(ispec-1)*5+2), coorg_recv(2,(ispec-1)*5+2)
           write(24,681) coorg_recv(1,(ispec-1)*5+3), coorg_recv(2,(ispec-1)*5+3)
           write(24,681) coorg_recv(1,(ispec-1)*5+4), coorg_recv(2,(ispec-1)*5+4)
           write(24,681) coorg_recv(1,(ispec-1)*5+5), coorg_recv(2,(ispec-1)*5+5)
           write(24,*) 'CO'
           write(24,*) greyscale_recv(ispec), ' setgray GF 0 setgray ST'

        enddo
        deallocate(coorg_recv)
        deallocate(greyscale_recv)

     enddo

  else
     call MPI_SEND (UPPER_LIMIT_DISPLAY, 1, MPI_INTEGER, &
        0, 42, MPI_COMM_WORLD, ier)
     call MPI_SEND (coorg_send, UPPER_LIMIT_DISPLAY*5*2, MPI_DOUBLE_PRECISION, &
        0, 42, MPI_COMM_WORLD, ier)
     call MPI_SEND (greyscale_send, UPPER_LIMIT_DISPLAY, MPI_INTEGER, &
        0, 42, MPI_COMM_WORLD, ier)
  endif
#endif

  if (myrank == 0) then

     write(24,*) '%'
     write(24,*) 'grestore'
     write(24,*) 'showpage'

     close(24)

     write(IOUT,*) 'End of creation of PostScript file with velocity model'

  endif

  if (myrank == 0) then

    write(IOUT,*)
    write(IOUT,*) 'Creating PostScript file with mesh partitioning'

!
!---- open PostScript file
!
  open(unit=24,file='OUTPUT_FILES/mesh_partitioning.ps',status='unknown')

!
!---- write PostScript header
!
  write(24,10) simulation_title
  write(24,*) '/CM {28.5 mul} def'
  write(24,*) '/LR {rlineto} def'
  write(24,*) '/LT {lineto} def'
  write(24,*) '/L {lineto} def'
  write(24,*) '/MR {rmoveto} def'
  write(24,*) '/MV {moveto} def'
  write(24,*) '/M {moveto} def'
  write(24,*) '/ST {stroke} def'
  write(24,*) '/CP {closepath} def'
  write(24,*) '/RG {setrgbcolor} def'
  write(24,*) '/GF {gsave fill grestore} def'
  write(24,*) '% different useful symbols'
  write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
  write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/Cross {gsave 0.05 CM setlinewidth'
  write(24,*) 'gsave 3 3 MR -6. -6. LR ST grestore'
  write(24,*) 'gsave 3 -3 MR -6. 6. LR ST grestore'
  write(24,*) '0.01 CM setlinewidth} def'
  write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
  write(24,*) '/Diamond {gsave 0.05 CM setlinewidth 0 4.2 MR'
  write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
  write(24,*) 'grestore 0.01 CM setlinewidth} def'
  write(24,*) '%'
  write(24,*) '% macro to draw the contour of the elements'
  write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
  write(24,*) '%'
  write(24,*) '.01 CM setlinewidth'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.35 CM scalefont setfont'
  write(24,*) '%'
  write(24,*) '/vshift ',-height/2,' CM def'
  write(24,*) '/Rshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop neg vshift MR show } def'
  write(24,*) '/Cshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
  write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM scalefont setfont} def'
  write(24,*) '%'
  write(24,*) 'gsave newpath 90 rotate'
  write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
  write(24,*) '%'

!
!--- write captions of PostScript figure
!
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.5 CM scalefont setfont'

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM scalefont setfont'
  write(24,*) '.4 .9 .9 setrgbcolor'
  write(24,*) '11 CM 1.1 CM MV'
  write(24,*) '(X axis) show'
  write(24,*) '%'
  write(24,*) '1.4 CM 9.5 CM MV'
  write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
  write(24,*) '(Z axis) show'
  write(24,*) 'grestore'
  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.7 CM scalefont setfont'
  write(24,*) '.8 0 .8 setrgbcolor'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(Mesh partitioning) show'
  write(24,*) 'grestore'
  write(24,*) '25.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(',simulation_title,') show'
  write(24,*) 'grestore'
  write(24,*) '26.45 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(2D Spectral Element Method) show'
  write(24,*) 'grestore'

  write(24,*) '%'
  write(24,*) '1 1 scale'
  write(24,*) '%'

!
!---- draw the spectral element mesh
!
  write(24,*) '%'
  write(24,*) '% spectral element mesh'
  write(24,*) '%'
  write(24,*) '0 setgray'

  num_ispec = 0
  endif

  do ispec = 1, UPPER_LIMIT_DISPLAY

     if ( myrank == 0 ) then
        num_ispec = num_ispec + 1
        write(24,*) '% elem ',num_ispec
     endif

  do i=1,pointsdisp
  do j=1,pointsdisp
  xinterp(i,j) = 0.d0
  zinterp(i,j) = 0.d0
  do in = 1,ngnod
    nnum = knods(in,ispec)
      xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
      zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
  enddo
  enddo
  enddo

  is = 1
  ir = 1
  x1 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z1 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  if ( myrank == 0 ) then
     write(24,*) 'mark'
     write(24,681) x1,z1
  else
     coorg_send(1,(ispec-1)*5+1) = x1
     coorg_send(2,(ispec-1)*5+1) = z1
  endif

! draw straight lines if elements have 4 nodes

  ir=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     coorg_send(1,(ispec-1)*5+2) = x2
     coorg_send(2,(ispec-1)*5+2) = z2
  endif

  ir=pointsdisp
  is=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     coorg_send(1,(ispec-1)*5+3) = x2
     coorg_send(2,(ispec-1)*5+3) = z2
  endif

  is=pointsdisp
  ir=1
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     coorg_send(1,(ispec-1)*5+4) = x2
     coorg_send(2,(ispec-1)*5+4) = z2
  endif

  ir=1
  is=2
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
     write(24,*) 'CO'
  else
     coorg_send(1,(ispec-1)*5+5) = x2
     coorg_send(2,(ispec-1)*5+5) = z2
  endif

  if ( myrank == 0 ) then
        write(24,*) red(1), green(1), blue(1), 'RG GF 0 setgray ST'
     endif

  enddo ! end of loop on all the spectral elements

#ifdef USE_MPI
  if (myrank == 0 ) then

      do iproc = 1, nproc-1

! use a different color for each material set
        icol = mod(iproc, NUM_COLORS) + 1

        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        allocate(coorg_recv(2,nspec_recv*5))
        call MPI_RECV (coorg_recv(1,1), nspec_recv*5*2, MPI_DOUBLE_PRECISION, &
            iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        do ispec = 1, nspec_recv
           num_ispec = num_ispec + 1
           write(24,*) '% elem ',num_ispec
           write(24,*) 'mark'
           write(24,681) coorg_recv(1,(ispec-1)*5+1), coorg_recv(2,(ispec-1)*5+1)
           write(24,681) coorg_recv(1,(ispec-1)*5+2), coorg_recv(2,(ispec-1)*5+2)
           write(24,681) coorg_recv(1,(ispec-1)*5+3), coorg_recv(2,(ispec-1)*5+3)
           write(24,681) coorg_recv(1,(ispec-1)*5+4), coorg_recv(2,(ispec-1)*5+4)
           write(24,681) coorg_recv(1,(ispec-1)*5+5), coorg_recv(2,(ispec-1)*5+5)
           write(24,*) 'CO'

           write(24,*) red(icol), green(icol), blue(icol), ' RG GF 0 setgray ST'

        enddo
        deallocate(coorg_recv)

     enddo

  else
     call MPI_SEND (UPPER_LIMIT_DISPLAY, 1, MPI_INTEGER, &
        0, 42, MPI_COMM_WORLD, ier)
     call MPI_SEND (coorg_send, UPPER_LIMIT_DISPLAY*5*2, MPI_DOUBLE_PRECISION, &
        0, 42, MPI_COMM_WORLD, ier)

  endif
#endif

 if (myrank == 0) then
   write(24,*) '%'
   write(24,*) 'grestore'
   write(24,*) 'showpage'

   close(24)

   write(IOUT,*) 'End of creation of PostScript file with partitioning'
   write(IOUT,*)
 endif

 10  format('%!PS-Adobe-2.0',/,'%%',/,'%% Title: ',a100,/,'%% Created by: Specfem2D',/,'%% Author: Dimitri Komatitsch',/,'%%')

 681 format(f6.2,1x,f6.2)

  end subroutine checkgrid


!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkgrid_setup_GLLper(percent_GLL,NGLLX_MAX_STABILITY)

  implicit none
  include "constants.h"

  integer :: NGLLX_MAX_STABILITY
  double precision :: percent_GLL(NGLLX_MAX_STABILITY)

  if( NGLLX_MAX_STABILITY /= 15 ) call exit_MPI('check NGLLX_MAX_STABILITY in checkgrid.f90')

! define percentage of smallest distance between GLL points for NGLLX points
! percentages were computed by calling the GLL points routine for each degree
  percent_GLL(:) = 100.d0

  percent_GLL(2) = 100.d0
  percent_GLL(3) = 50.d0
  percent_GLL(4) = 27.639320225002102d0
  percent_GLL(5) = 17.267316464601141d0
  percent_GLL(6) = 11.747233803526763d0
  percent_GLL(7) = 8.4888051860716516d0
  percent_GLL(8) = 6.4129925745196719d0
  percent_GLL(9) = 5.0121002294269914d0
  percent_GLL(10) = 4.0233045916770571d0
  percent_GLL(11) = 3.2999284795970416d0
  percent_GLL(12) = 2.7550363888558858d0
  percent_GLL(13) = 2.3345076678918053d0
  percent_GLL(14) = 2.0032477366369594d0
  percent_GLL(15) = 1.7377036748080721d0

! convert to real percentage
  percent_GLL(:) = percent_GLL(:) / 100.d0

  if(NGLLX > NGLLX_MAX_STABILITY) then
    call exit_MPI('cannot estimate the stability condition for that degree')
  endif

  end subroutine checkgrid_setup_GLLper


!
!-------------------------------------------------------------------------------------------------
!

  subroutine checkgrid_setup_colorp(red,green,blue,NUM_COLORS)

! color palette

  implicit none
  integer :: NUM_COLORS
  double precision, dimension(NUM_COLORS) :: red,green,blue

  if( NUM_COLORS /= 236 ) call exit_MPI('check NUM_COLORS in checkgrid.f90')

! red
  red(1) = 1.00000000000000
  green(1) = 0.000000000000000E+000
  blue(1) = 0.000000000000000E+000

! DodgerBlue2
  red(2) = 0.109803921568627
  green(2) = 0.525490196078431
  blue(2) = 0.933333333333333

! gold
  red(3) = 1.00000000000000
  green(3) = 0.840000000000000
  blue(3) = 0.000000000000000E+000

! springgreen
  red(4) = 0.000000000000000E+000
  green(4) = 1.00000000000000
  blue(4) = 0.500000000000000

! NavajoWhite
  red(5) = 1.00000000000000
  green(5) = 0.870588235294118
  blue(5) = 0.678431372549020

! SteelBlue3
  red(6) = 0.309803921568627
  green(6) = 0.580392156862745
  blue(6) = 0.803921568627451

! Ivory3
  red(7) = 0.803921568627451
  green(7) = 0.803921568627451
  blue(7) = 0.756862745098039

! SkyBlue4
  red(8) = 0.290196078431373
  green(8) = 0.439215686274510
  blue(8) = 0.545098039215686

! Snow
  red(9) = 0.980392156862745
  green(9) = 0.980392156862745
  blue(9) = 0.980392156862745

! SteelBlue
  red(10) = 0.274509803921569
  green(10) = 0.509803921568627
  blue(10) = 0.705882352941177

! Bisque3
  red(11) = 0.803921568627451
  green(11) = 0.717647058823529
  blue(11) = 0.619607843137255

! Salmon
  red(12) = 0.980392156862745
  green(12) = 0.501960784313725
  blue(12) = 0.447058823529412

! SlateBlue2
  red(13) = 0.478431372549020
  green(13) = 0.403921568627451
  blue(13) = 0.933333333333333

! NavajoWhite2
  red(14) = 0.933333333333333
  green(14) = 0.811764705882353
  blue(14) = 0.631372549019608

! MediumBlue
  red(15) = 0.000000000000000E+000
  green(15) = 0.000000000000000E+000
  blue(15) = 0.803921568627451

! LightCoral
  red(16) = 0.941176470588235
  green(16) = 0.501960784313725
  blue(16) = 0.501960784313725

! FloralWhite
  red(17) = 1.00000000000000
  green(17) = 0.980392156862745
  blue(17) = 0.941176470588235

! Cornsilk3
  red(18) = 0.803921568627451
  green(18) = 0.784313725490196
  blue(18) = 0.694117647058824

! GhostWhite
  red(19) = 0.972549019607843
  green(19) = 0.972549019607843
  blue(19) = 1.00000000000000

! blue
  red(20) = 0.000000000000000E+000
  green(20) = 0.000000000000000E+000
  blue(20) = 1.00000000000000

! Linen
  red(21) = 0.980392156862745
  green(21) = 0.941176470588235
  blue(21) = 0.901960784313726

! peachpuff
  red(22) = 1.00000000000000
  green(22) = 0.850000000000000
  blue(22) = 0.730000000000000

! Cornsilk1
  red(23) = 1.00000000000000
  green(23) = 0.972549019607843
  blue(23) = 0.862745098039216

! LightSalmon
  red(24) = 1.00000000000000
  green(24) = 0.627450980392157
  blue(24) = 0.478431372549020

! DeepSkyBlue1
  red(25) = 0.000000000000000E+000
  green(25) = 0.749019607843137
  blue(25) = 1.00000000000000

! LemonChiffon4
  red(26) = 0.545098039215686
  green(26) = 0.537254901960784
  blue(26) = 0.439215686274510

! PeachPuff1
  red(27) = 1.00000000000000
  green(27) = 0.854901960784314
  blue(27) = 0.725490196078431

! BlanchedAlmond
  red(28) = 1.00000000000000
  green(28) = 0.921568627450980
  blue(28) = 0.803921568627451

! SlateBlue3
  red(29) = 0.411764705882353
  green(29) = 0.349019607843137
  blue(29) = 0.803921568627451

! LightSkyBlue1
  red(30) = 0.690196078431373
  green(30) = 0.886274509803922
  blue(30) = 1.00000000000000

! DarkViolet
  red(31) = 0.580392156862745
  green(31) = 0.000000000000000E+000
  blue(31) = 0.827450980392157

! Azure3
  red(32) = 0.756862745098039
  green(32) = 0.803921568627451
  blue(32) = 0.803921568627451

! LavenderBlush3
  red(33) = 0.803921568627451
  green(33) = 0.756862745098039
  blue(33) = 0.772549019607843

! Honeydew1
  red(34) = 0.941176470588235
  green(34) = 1.00000000000000
  blue(34) = 0.941176470588235

! Ivory2
  red(35) = 0.933333333333333
  green(35) = 0.933333333333333
  blue(35) = 0.878431372549020

! RosyBrown
  red(36) = 0.737254901960784
  green(36) = 0.560784313725490
  blue(36) = 0.560784313725490

! Thistle
  red(37) = 0.847058823529412
  green(37) = 0.749019607843137
  blue(37) = 0.847058823529412

! Orange
  red(38) = 1.00000000000000
  green(38) = 0.647058823529412
  blue(38) = 0.000000000000000E+000

! DarkSeaGreen
  red(39) = 0.560784313725490
  green(39) = 0.737254901960784
  blue(39) = 0.560784313725490

! Moccasin
  red(40) = 1.00000000000000
  green(40) = 0.894117647058824
  blue(40) = 0.709803921568627

! DeepSkyBlue2
  red(41) = 0.000000000000000E+000
  green(41) = 0.698039215686274
  blue(41) = 0.933333333333333

! SlateGray4
  red(42) = 0.423529411764706
  green(42) = 0.482352941176471
  blue(42) = 0.545098039215686

! Beige
  red(43) = 0.960784313725490
  green(43) = 0.960784313725490
  blue(43) = 0.862745098039216

! Gold
  red(44) = 1.00000000000000
  green(44) = 0.843137254901961
  blue(44) = 0.000000000000000E+000

! SlateBlue
  red(45) = 0.415686274509804
  green(45) = 0.352941176470588
  blue(45) = 0.803921568627451

! SteelBlue1
  red(46) = 0.388235294117647
  green(46) = 0.721568627450980
  blue(46) = 1.00000000000000

! SaddleBrown
  red(47) = 0.545098039215686
  green(47) = 0.270588235294118
  blue(47) = 7.450980392156863E-002

! Pink
  red(48) = 1.00000000000000
  green(48) = 0.752941176470588
  blue(48) = 0.796078431372549

! Black
  red(49) = 0.000000000000000E+000
  green(49) = 0.000000000000000E+000
  blue(49) = 0.000000000000000E+000

! SlateGrey
  red(50) = 0.439215686274510
  green(50) = 0.501960784313725
  blue(50) = 0.564705882352941

! Ivory
  red(51) = 1.00000000000000
  green(51) = 1.00000000000000
  blue(51) = 0.941176470588235

! OliveDrab
  red(52) = 0.419607843137255
  green(52) = 0.556862745098039
  blue(52) = 0.137254901960784

! Ivory1
  red(53) = 1.00000000000000
  green(53) = 1.00000000000000
  blue(53) = 0.941176470588235

! SkyBlue
  red(54) = 0.529411764705882
  green(54) = 0.807843137254902
  blue(54) = 0.921568627450980

! MistyRose3
  red(55) = 0.803921568627451
  green(55) = 0.717647058823529
  blue(55) = 0.709803921568627

! LimeGreen
  red(56) = 0.196078431372549
  green(56) = 0.803921568627451
  blue(56) = 0.196078431372549

! Purple
  red(57) = 0.627450980392157
  green(57) = 0.125490196078431
  blue(57) = 0.941176470588235

! SkyBlue2
  red(58) = 0.494117647058824
  green(58) = 0.752941176470588
  blue(58) = 0.933333333333333

! Red
  red(59) = 1.00000000000000
  green(59) = 0.000000000000000E+000
  blue(59) = 0.000000000000000E+000

! DarkKhaki
  red(60) = 0.741176470588235
  green(60) = 0.717647058823529
  blue(60) = 0.419607843137255

! MediumTurquoise
  red(61) = 0.282352941176471
  green(61) = 0.819607843137255
  blue(61) = 0.800000000000000

! Grey
  red(62) = 0.745098039215686
  green(62) = 0.745098039215686
  blue(62) = 0.745098039215686

! Coral
  red(63) = 1.00000000000000
  green(63) = 0.498039215686275
  blue(63) = 0.313725490196078

! NavajoWhite4
  red(64) = 0.545098039215686
  green(64) = 0.474509803921569
  blue(64) = 0.368627450980392

! SlateBlue4
  red(65) = 0.278431372549020
  green(65) = 0.235294117647059
  blue(65) = 0.545098039215686

! RoyalBlue4
  red(66) = 0.152941176470588
  green(66) = 0.250980392156863
  blue(66) = 0.545098039215686

! YellowGreen
  red(67) = 0.603921568627451
  green(67) = 0.803921568627451
  blue(67) = 0.196078431372549

! DeepSkyBlue3
  red(68) = 0.000000000000000E+000
  green(68) = 0.603921568627451
  blue(68) = 0.803921568627451

! goldenrod
  red(69) = 0.854901960784314
  green(69) = 0.647058823529412
  blue(69) = 0.125490196078431

! AntiqueWhite4
  red(70) = 0.545098039215686
  green(70) = 0.513725490196078
  blue(70) = 0.470588235294118

! lemonchiffon
  red(71) = 1.00000000000000
  green(71) = 0.980000000000000
  blue(71) = 0.800000000000000

! GreenYellow
  red(72) = 0.678431372549020
  green(72) = 1.00000000000000
  blue(72) = 0.184313725490196

! LightSlateGray
  red(73) = 0.466666666666667
  green(73) = 0.533333333333333
  blue(73) = 0.600000000000000

! RoyalBlue
  red(74) = 0.254901960784314
  green(74) = 0.411764705882353
  blue(74) = 0.882352941176471

! DarkGreen
  red(75) = 0.000000000000000E+000
  green(75) = 0.392156862745098
  blue(75) = 0.000000000000000E+000

! NavajoWhite3
  red(76) = 0.803921568627451
  green(76) = 0.701960784313725
  blue(76) = 0.545098039215686

! Azure1
  red(77) = 0.941176470588235
  green(77) = 1.00000000000000
  blue(77) = 1.00000000000000

! PowderBlue
  red(78) = 0.690196078431373
  green(78) = 0.878431372549020
  blue(78) = 0.901960784313726

! slateblue
  red(79) = 0.420000000000000
  green(79) = 0.350000000000000
  blue(79) = 0.800000000000000

! MediumOrchid
  red(80) = 0.729411764705882
  green(80) = 0.333333333333333
  blue(80) = 0.827450980392157

! turquoise
  red(81) = 0.250000000000000
  green(81) = 0.880000000000000
  blue(81) = 0.820000000000000

! Snow1
  red(82) = 1.00000000000000
  green(82) = 0.980392156862745
  blue(82) = 0.980392156862745

! violet
  red(83) = 0.930000000000000
  green(83) = 0.510000000000000
  blue(83) = 0.930000000000000

! DeepPink
  red(84) = 1.00000000000000
  green(84) = 7.843137254901961E-002
  blue(84) = 0.576470588235294

! MistyRose4
  red(85) = 0.545098039215686
  green(85) = 0.490196078431373
  blue(85) = 0.482352941176471

! PeachPuff3
  red(86) = 0.803921568627451
  green(86) = 0.686274509803922
  blue(86) = 0.584313725490196

! MediumSeaGreen
  red(87) = 0.235294117647059
  green(87) = 0.701960784313725
  blue(87) = 0.443137254901961

! Honeydew4
  red(88) = 0.513725490196078
  green(88) = 0.545098039215686
  blue(88) = 0.513725490196078

! Tan
  red(89) = 0.823529411764706
  green(89) = 0.705882352941177
  blue(89) = 0.549019607843137

! DarkGoldenrod
  red(90) = 0.721568627450980
  green(90) = 0.525490196078431
  blue(90) = 4.313725490196078E-002

! Blue2
  red(91) = 0.000000000000000E+000
  green(91) = 0.000000000000000E+000
  blue(91) = 0.933333333333333

! Maroon
  red(92) = 0.690196078431373
  green(92) = 0.188235294117647
  blue(92) = 0.376470588235294

! LightSkyBlue3
  red(93) = 0.552941176470588
  green(93) = 0.713725490196078
  blue(93) = 0.803921568627451

! LemonChiffon2
  red(94) = 0.933333333333333
  green(94) = 0.913725490196078
  blue(94) = 0.749019607843137

! Snow3
  red(95) = 0.803921568627451
  green(95) = 0.788235294117647
  blue(95) = 0.788235294117647

! Ivory4
  red(96) = 0.545098039215686
  green(96) = 0.545098039215686
  blue(96) = 0.513725490196078

! AntiqueWhite3
  red(97) = 0.803921568627451
  green(97) = 0.752941176470588
  blue(97) = 0.690196078431373

! Bisque4
  red(98) = 0.545098039215686
  green(98) = 0.490196078431373
  blue(98) = 0.419607843137255

! Snow2
  red(99) = 0.933333333333333
  green(99) = 0.913725490196078
  blue(99) = 0.913725490196078

! SlateGray1
  red(100) = 0.776470588235294
  green(100) = 0.886274509803922
  blue(100) = 1.00000000000000

! Seashell2
  red(101) = 0.933333333333333
  green(101) = 0.898039215686275
  blue(101) = 0.870588235294118

! Aquamarine
  red(102) = 0.498039215686275
  green(102) = 1.00000000000000
  blue(102) = 0.831372549019608

! SlateGray2
  red(103) = 0.725490196078431
  green(103) = 0.827450980392157
  blue(103) = 0.933333333333333

! White
  red(104) = 1.00000000000000
  green(104) = 1.00000000000000
  blue(104) = 1.00000000000000

! LavenderBlush
  red(105) = 1.00000000000000
  green(105) = 0.941176470588235
  blue(105) = 0.960784313725490

! DodgerBlue3
  red(106) = 9.411764705882353E-002
  green(106) = 0.454901960784314
  blue(106) = 0.803921568627451

! RoyalBlue3
  red(107) = 0.227450980392157
  green(107) = 0.372549019607843
  blue(107) = 0.803921568627451

! LightYellow
  red(108) = 1.00000000000000
  green(108) = 1.00000000000000
  blue(108) = 0.878431372549020

! DeepSkyBlue
  red(109) = 0.000000000000000E+000
  green(109) = 0.749019607843137
  blue(109) = 1.00000000000000

! AntiqueWhite2
  red(110) = 0.933333333333333
  green(110) = 0.874509803921569
  blue(110) = 0.800000000000000

! CornflowerBlue
  red(111) = 0.392156862745098
  green(111) = 0.584313725490196
  blue(111) = 0.929411764705882

! PeachPuff4
  red(112) = 0.545098039215686
  green(112) = 0.466666666666667
  blue(112) = 0.396078431372549

! SpringGreen
  red(113) = 0.000000000000000E+000
  green(113) = 1.00000000000000
  blue(113) = 0.498039215686275

! Honeydew
  red(114) = 0.941176470588235
  green(114) = 1.00000000000000
  blue(114) = 0.941176470588235

! Honeydew2
  red(115) = 0.878431372549020
  green(115) = 0.933333333333333
  blue(115) = 0.878431372549020

! LightSeaGreen
  red(116) = 0.125490196078431
  green(116) = 0.698039215686274
  blue(116) = 0.666666666666667

! NavyBlue
  red(117) = 0.000000000000000E+000
  green(117) = 0.000000000000000E+000
  blue(117) = 0.501960784313725

! Azure4
  red(118) = 0.513725490196078
  green(118) = 0.545098039215686
  blue(118) = 0.545098039215686

! MediumAquamarine
  red(119) = 0.400000000000000
  green(119) = 0.803921568627451
  blue(119) = 0.666666666666667

! SkyBlue3
  red(120) = 0.423529411764706
  green(120) = 0.650980392156863
  blue(120) = 0.803921568627451

! LavenderBlush2
  red(121) = 0.933333333333333
  green(121) = 0.878431372549020
  blue(121) = 0.898039215686275

! Bisque1
  red(122) = 1.00000000000000
  green(122) = 0.894117647058824
  blue(122) = 0.768627450980392

! DarkOrange
  red(123) = 1.00000000000000
  green(123) = 0.549019607843137
  blue(123) = 0.000000000000000E+000

! LightSteelBlue
  red(124) = 0.690196078431373
  green(124) = 0.768627450980392
  blue(124) = 0.870588235294118

! SteelBlue2
  red(125) = 0.360784313725490
  green(125) = 0.674509803921569
  blue(125) = 0.933333333333333

! LemonChiffon3
  red(126) = 0.803921568627451
  green(126) = 0.788235294117647
  blue(126) = 0.647058823529412

! DarkSlateBlue
  red(127) = 0.282352941176471
  green(127) = 0.239215686274510
  blue(127) = 0.545098039215686

! Seashell
  red(128) = 1.00000000000000
  green(128) = 0.960784313725490
  blue(128) = 0.933333333333333

! Firebrick
  red(129) = 0.698039215686274
  green(129) = 0.133333333333333
  blue(129) = 0.133333333333333

! LightGray
  red(130) = 0.827450980392157
  green(130) = 0.827450980392157
  blue(130) = 0.827450980392157

! Blue
  red(131) = 0.000000000000000E+000
  green(131) = 0.000000000000000E+000
  blue(131) = 1.00000000000000

! Bisque2
  red(132) = 0.933333333333333
  green(132) = 0.835294117647059
  blue(132) = 0.717647058823529

! WhiteSmoke
  red(133) = 0.960784313725490
  green(133) = 0.960784313725490
  blue(133) = 0.960784313725490

! SeaGreen
  red(134) = 0.180392156862745
  green(134) = 0.545098039215686
  blue(134) = 0.341176470588235

! Burlywood
  red(135) = 0.870588235294118
  green(135) = 0.721568627450980
  blue(135) = 0.529411764705882

! RoyalBlue2
  red(136) = 0.262745098039216
  green(136) = 0.431372549019608
  blue(136) = 0.933333333333333

! RoyalBlue1
  red(137) = 0.282352941176471
  green(137) = 0.462745098039216
  blue(137) = 1.00000000000000

! SteelBlue4
  red(138) = 0.211764705882353
  green(138) = 0.392156862745098
  blue(138) = 0.545098039215686

! AliceBlue
  red(139) = 0.941176470588235
  green(139) = 0.972549019607843
  blue(139) = 1.00000000000000

! LightSlateBlue
  red(140) = 0.517647058823529
  green(140) = 0.439215686274510
  blue(140) = 1.00000000000000

! MistyRose1
  red(141) = 1.00000000000000
  green(141) = 0.894117647058824
  blue(141) = 0.882352941176471

! SandyBrown
  red(142) = 0.956862745098039
  green(142) = 0.643137254901961
  blue(142) = 0.376470588235294

! DarkOliveGreen
  red(143) = 0.333333333333333
  green(143) = 0.419607843137255
  blue(143) = 0.184313725490196

! Yellow
  red(144) = 1.00000000000000
  green(144) = 1.00000000000000
  blue(144) = 0.000000000000000E+000

! SlateGray3
  red(145) = 0.623529411764706
  green(145) = 0.713725490196078
  blue(145) = 0.803921568627451

! HotPink
  red(146) = 1.00000000000000
  green(146) = 0.411764705882353
  blue(146) = 0.705882352941177

! Violet
  red(147) = 0.933333333333333
  green(147) = 0.509803921568627
  blue(147) = 0.933333333333333

! LightSkyBlue
  red(148) = 0.529411764705882
  green(148) = 0.807843137254902
  blue(148) = 0.980392156862745

! Cornsilk2
  red(149) = 0.933333333333333
  green(149) = 0.909803921568627
  blue(149) = 0.803921568627451

! MidnightBlue
  red(150) = 9.803921568627451E-002
  green(150) = 9.803921568627451E-002
  blue(150) = 0.439215686274510

! AntiqueWhite
  red(151) = 0.980392156862745
  green(151) = 0.921568627450980
  blue(151) = 0.843137254901961

! PaleGreen
  red(152) = 0.596078431372549
  green(152) = 0.984313725490196
  blue(152) = 0.596078431372549

! MedSpringGreen
  red(153) = 0.000000000000000E+000
  green(153) = 0.980392156862745
  blue(153) = 0.603921568627451

! DodgerBlue1
  red(154) = 0.117647058823529
  green(154) = 0.564705882352941
  blue(154) = 1.00000000000000

! Blue3
  red(155) = 0.000000000000000E+000
  green(155) = 0.000000000000000E+000
  blue(155) = 0.803921568627451

! Cyan
  red(156) = 0.000000000000000E+000
  green(156) = 1.00000000000000
  blue(156) = 1.00000000000000

! LemonChiffon
  red(157) = 1.00000000000000
  green(157) = 0.980392156862745
  blue(157) = 0.803921568627451

! mediumorchid
  red(158) = 0.730000000000000
  green(158) = 0.330000000000000
  blue(158) = 0.830000000000000

! Turquoise
  red(159) = 0.250980392156863
  green(159) = 0.878431372549020
  blue(159) = 0.815686274509804

! IndianRed
  red(160) = 0.803921568627451
  green(160) = 0.360784313725490
  blue(160) = 0.360784313725490

! DodgerBlue
  red(161) = 0.117647058823529
  green(161) = 0.564705882352941
  blue(161) = 1.00000000000000

! Seashell3
  red(162) = 0.803921568627451
  green(162) = 0.772549019607843
  blue(162) = 0.749019607843137

! BlueViolet
  red(163) = 0.541176470588235
  green(163) = 0.168627450980392
  blue(163) = 0.886274509803922

! DeepSkyBlue4
  red(164) = 0.000000000000000E+000
  green(164) = 0.407843137254902
  blue(164) = 0.545098039215686

! PaleVioletRed
  red(165) = 0.858823529411765
  green(165) = 0.439215686274510
  blue(165) = 0.576470588235294

! Azure2
  red(166) = 0.878431372549020
  green(166) = 0.933333333333333
  blue(166) = 0.933333333333333

! greenyellow
  red(167) = 0.680000000000000
  green(167) = 1.00000000000000
  blue(167) = 0.180000000000000

! LightGoldenrod
  red(168) = 0.933333333333333
  green(168) = 0.866666666666667
  blue(168) = 0.509803921568627

! MistyRose
  red(169) = 1.00000000000000
  green(169) = 0.894117647058824
  blue(169) = 0.882352941176471

! LightSkyBlue4
  red(170) = 0.376470588235294
  green(170) = 0.482352941176471
  blue(170) = 0.545098039215686

! OrangeRed
  red(171) = 1.00000000000000
  green(171) = 0.270588235294118
  blue(171) = 0.000000000000000E+000

! DimGrey
  red(172) = 0.411764705882353
  green(172) = 0.411764705882353
  blue(172) = 0.411764705882353

! MediumVioletRed
  red(173) = 0.780392156862745
  green(173) = 8.235294117647059E-002
  blue(173) = 0.521568627450980

! DarkSlateGray
  red(174) = 0.184313725490196
  green(174) = 0.309803921568627
  blue(174) = 0.309803921568627

! yellow
  red(175) = 1.00000000000000
  green(175) = 1.00000000000000
  blue(175) = 0.000000000000000E+000

! Plum
  red(176) = 0.866666666666667
  green(176) = 0.627450980392157
  blue(176) = 0.866666666666667

! DarkTurquoise
  red(177) = 0.000000000000000E+000
  green(177) = 0.807843137254902
  blue(177) = 0.819607843137255

! DodgerBlue4
  red(178) = 6.274509803921569E-002
  green(178) = 0.305882352941176
  blue(178) = 0.545098039215686

! Cornsilk
  red(179) = 1.00000000000000
  green(179) = 0.972549019607843
  blue(179) = 0.862745098039216

! SkyBlue1
  red(180) = 0.529411764705882
  green(180) = 0.807843137254902
  blue(180) = 1.00000000000000

! Seashell1
  red(181) = 1.00000000000000
  green(181) = 0.960784313725490
  blue(181) = 0.933333333333333

! lavender
  red(182) = 0.901960784313726
  green(182) = 0.901960784313726
  blue(182) = 0.980392156862745

! Snow4
  red(183) = 0.545098039215686
  green(183) = 0.537254901960784
  blue(183) = 0.537254901960784

! Peru
  red(184) = 0.803921568627451
  green(184) = 0.521568627450980
  blue(184) = 0.247058823529412

! PeachPuff
  red(185) = 1.00000000000000
  green(185) = 0.854901960784314
  blue(185) = 0.725490196078431

! Green
  red(186) = 0.000000000000000E+000
  green(186) = 1.00000000000000
  blue(186) = 0.000000000000000E+000

! Blue1
  red(187) = 0.000000000000000E+000
  green(187) = 0.000000000000000E+000
  blue(187) = 1.00000000000000

! Seashell4
  red(188) = 0.545098039215686
  green(188) = 0.525490196078431
  blue(188) = 0.509803921568627

! dodgerblue
  red(189) = 0.120000000000000
  green(189) = 0.560000000000000
  blue(189) = 1.00000000000000

! MistyRose2
  red(190) = 0.933333333333333
  green(190) = 0.835294117647059
  blue(190) = 0.823529411764706

! Tomato
  red(191) = 1.00000000000000
  green(191) = 0.388235294117647
  blue(191) = 0.278431372549020

! Wheat
  red(192) = 0.960784313725490
  green(192) = 0.870588235294118
  blue(192) = 0.701960784313725

! LightBlue
  red(193) = 0.678431372549020
  green(193) = 0.847058823529412
  blue(193) = 0.901960784313726

! Chocolate
  red(194) = 0.823529411764706
  green(194) = 0.411764705882353
  blue(194) = 0.117647058823529

! Blue4
  red(195) = 0.000000000000000E+000
  green(195) = 0.000000000000000E+000
  blue(195) = 0.545098039215686

! LavenderBlush1
  red(196) = 1.00000000000000
  green(196) = 0.941176470588235
  blue(196) = 0.960784313725490

! Magenta
  red(197) = 1.00000000000000
  green(197) = 0.000000000000000E+000
  blue(197) = 1.00000000000000

! darkturquoise
  red(198) = 0.000000000000000E+000
  green(198) = 0.810000000000000
  blue(198) = 0.820000000000000

! blueviolet
  red(199) = 0.540000000000000
  green(199) = 0.170000000000000
  blue(199) = 0.890000000000000

! MintCream
  red(200) = 0.960784313725490
  green(200) = 1.00000000000000
  blue(200) = 0.980392156862745

! PaleGoldenrod
  red(201) = 0.933333333333333
  green(201) = 0.909803921568627
  blue(201) = 0.666666666666667

! MediumPurple
  red(202) = 0.576470588235294
  green(202) = 0.439215686274510
  blue(202) = 0.858823529411765

! PapayaWhip
  red(203) = 1.00000000000000
  green(203) = 0.937254901960784
  blue(203) = 0.835294117647059

! LavenderBlush4
  red(204) = 0.545098039215686
  green(204) = 0.513725490196078
  blue(204) = 0.525490196078431

! Cornsilk4
  red(205) = 0.545098039215686
  green(205) = 0.533333333333333
  blue(205) = 0.470588235294118

! LtGoldenrodYello
  red(206) = 0.980392156862745
  green(206) = 0.980392156862745
  blue(206) = 0.823529411764706

! limegreen
  red(207) = 0.200000000000000
  green(207) = 0.800000000000000
  blue(207) = 0.200000000000000

! LemonChiffon1
  red(208) = 1.00000000000000
  green(208) = 0.980392156862745
  blue(208) = 0.803921568627451

! DarkOrchid
  red(209) = 0.600000000000000
  green(209) = 0.196078431372549
  blue(209) = 0.800000000000000

! SlateBlue1
  red(210) = 0.513725490196078
  green(210) = 0.435294117647059
  blue(210) = 1.00000000000000

! chartreuse
  red(211) = 0.500000000000000
  green(211) = 1.00000000000000
  blue(211) = 0.000000000000000E+000

! PaleTurquoise
  red(212) = 0.686274509803922
  green(212) = 0.933333333333333
  blue(212) = 0.933333333333333

! NavajoWhite1
  red(213) = 1.00000000000000
  green(213) = 0.870588235294118
  blue(213) = 0.678431372549020

! LightSkyBlue2
  red(214) = 0.643137254901961
  green(214) = 0.827450980392157
  blue(214) = 0.933333333333333

! VioletRed
  red(215) = 0.815686274509804
  green(215) = 0.125490196078431
  blue(215) = 0.564705882352941

! mocassin
  red(216) = 1.00000000000000
  green(216) = 0.890000000000000
  blue(216) = 0.710000000000000

! OldLace
  red(217) = 0.992156862745098
  green(217) = 0.960784313725490
  blue(217) = 0.901960784313726

! deeppink
  red(218) = 1.00000000000000
  green(218) = 8.000000000000000E-002
  blue(218) = 0.580000000000000

! Honeydew3
  red(219) = 0.756862745098039
  green(219) = 0.803921568627451
  blue(219) = 0.756862745098039

! Gainsboro
  red(220) = 0.862745098039216
  green(220) = 0.862745098039216
  blue(220) = 0.862745098039216

! DarkSalmon
  red(221) = 0.913725490196078
  green(221) = 0.588235294117647
  blue(221) = 0.478431372549020

! AntiqueWhite1
  red(222) = 1.00000000000000
  green(222) = 0.937254901960784
  blue(222) = 0.858823529411765

! LightCyan
  red(223) = 0.878431372549020
  green(223) = 1.00000000000000
  blue(223) = 1.00000000000000

! ForestGreen
  red(224) = 0.133333333333333
  green(224) = 0.545098039215686
  blue(224) = 0.133333333333333

! Orchid
  red(225) = 0.854901960784314
  green(225) = 0.439215686274510
  blue(225) = 0.839215686274510

! PeachPuff2
  red(226) = 0.933333333333333
  green(226) = 0.796078431372549
  blue(226) = 0.678431372549020

! LightPink
  red(227) = 1.00000000000000
  green(227) = 0.713725490196078
  blue(227) = 0.756862745098039

! Sienna
  red(228) = 0.627450980392157
  green(228) = 0.321568627450980
  blue(228) = 0.176470588235294

! darkorchid
  red(229) = 0.600000000000000
  green(229) = 0.200000000000000
  blue(229) = 0.800000000000000

! MediumSlateBlue
  red(230) = 0.482352941176471
  green(230) = 0.407843137254902
  blue(230) = 0.933333333333333

! CadetBlue
  red(231) = 0.372549019607843
  green(231) = 0.619607843137255
  blue(231) = 0.627450980392157

! LawnGreen
  red(232) = 0.486274509803922
  green(232) = 0.988235294117647
  blue(232) = 0.000000000000000E+000

! Chartreuse
  red(233) = 0.498039215686275
  green(233) = 1.00000000000000
  blue(233) = 0.000000000000000E+000

! Brown
  red(234) = 0.647058823529412
  green(234) = 0.164705882352941
  blue(234) = 0.164705882352941

! Azure
  red(235) = 0.941176470588235
  green(235) = 1.00000000000000
  blue(235) = 1.00000000000000

! Bisque
  red(236) = 1.00000000000000
  green(236) = 0.894117647058824
  blue(236) = 0.768627450980392

  end subroutine checkgrid_setup_colorp



!========================================================================
!
!                   S P E C F E M 2 D  Version 6.0
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France,
! and Princeton University, USA.
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

subroutine setup_sources_receivers(NSOURCE,assign_external_model,initialfield,numat,source_type,&
     coord,ibool,kmato,npoin,nspec,nelem_acoustic_surface,acoustic_surface,elastic,poroelastic, &
     x_source,z_source,ix_source,iz_source,ispec_selected_source,ispec_selected_rec,iglob_source, &
     is_proc_source,nb_proc_source,rho_at_source_location,ipass,&
     sourcearray,Mxx,Mzz,Mxz,xix,xiz,gammax,gammaz,xigll,zigll,npgeo,density,rhoext,&
     nproc,myrank,xi_source,gamma_source,coorg,knods,ngnod, &
     nrec,nrecloc,recloc,which_proc_receiver,st_xval,st_zval, &
     xi_receiver,gamma_receiver,station_name,network_name,x_final_receiver,z_final_receiver)

  implicit none

  include "constants.h"

  logical :: assign_external_model,initialfield
  integer :: NSOURCE, numat
  integer :: npgeo,ngnod,myrank,ipass,nproc
  integer :: npoin,nspec,nelem_acoustic_surface

  ! Gauss-Lobatto-Legendre points 
  double precision, dimension(NGLLX) :: xigll
  double precision, dimension(NGLLZ) :: zigll

  ! for receivers 
  integer  :: nrec,nrecloc
  integer, dimension(nrec) :: recloc, which_proc_receiver
  integer, dimension(nrec) :: ispec_selected_rec
  double precision, dimension(nrec) :: xi_receiver,gamma_receiver,st_xval,st_zval
  double precision, dimension(nrec) :: x_final_receiver, z_final_receiver

  ! timing information for the stations
  character(len=MAX_LENGTH_STATION_NAME), dimension(nrec) :: station_name
  character(len=MAX_LENGTH_NETWORK_NAME), dimension(nrec) :: network_name

  ! for sources
  double precision :: rho_at_source_location
  integer, dimension(NSOURCE) :: source_type
  integer, dimension(NSOURCE) :: ispec_selected_source,iglob_source,ix_source,iz_source,is_proc_source,nb_proc_source  
  real(kind=CUSTOM_REAL), dimension(NSOURCE,NDIM,NGLLX,NGLLZ) :: sourcearray
  double precision, dimension(NSOURCE) :: x_source,z_source,xi_source,gamma_source,Mxx,Mzz,Mxz

  logical, dimension(nspec) :: elastic,poroelastic
  integer, dimension(nspec) :: kmato
  integer, dimension(ngnod,nspec) :: knods
  integer, dimension(5,nelem_acoustic_surface) :: acoustic_surface
  integer, dimension(NGLLX,NGLLZ,nspec)  :: ibool
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec)  :: xix,xiz,gammax,gammaz
  double precision, dimension(NDIM,npgeo) :: coorg
  double precision, dimension(2,numat) :: density
  double precision, dimension(NGLLX,NGLLZ,nspec) :: rhoext
  double precision, dimension(NDIM,npoin) :: coord

  ! Local variables
  integer i_source,ispec,ispec_acoustic_surface,i,j,iglob

  do i_source=1,NSOURCE

     if(source_type(i_source) == 1) then

        ! collocated force source
        call locate_source_force(coord,ibool,npoin,nspec,x_source(i_source),z_source(i_source), &
             ix_source(i_source),iz_source(i_source),ispec_selected_source(i_source),iglob_source(i_source), &
             is_proc_source(i_source),nb_proc_source(i_source),ipass)

        ! get density at the source in order to implement collocated force with the right
        ! amplitude later
        if(is_proc_source(i_source) == 1) then
           rho_at_source_location  = density(1,kmato(ispec_selected_source(i_source)))
           ! external velocity model
           if(assign_external_model) rho_at_source_location = &
                rhoext(ix_source(i_source),iz_source(i_source),ispec_selected_source(i_source))
        endif

        ! check that acoustic source is not exactly on the free surface because pressure is zero there
        if(is_proc_source(i_source) == 1) then
           do ispec_acoustic_surface = 1,nelem_acoustic_surface
              ispec = acoustic_surface(1,ispec_acoustic_surface)
              if( .not. elastic(ispec) .and. .not. poroelastic(ispec) .and. ispec == ispec_selected_source(i_source) ) then
                 do j = acoustic_surface(4,ispec_acoustic_surface), acoustic_surface(5,ispec_acoustic_surface)
                    do i = acoustic_surface(2,ispec_acoustic_surface), acoustic_surface(3,ispec_acoustic_surface)
                       iglob = ibool(i,j,ispec)
                       if ( iglob_source(i_source) == iglob ) then

            call exit_MPI('an acoustic source cannot be located exactly on the free surface because pressure is zero there')

                       endif
                    enddo
                 enddo
              endif
           enddo
        endif

     else if(source_type(i_source) == 2) then
        ! moment-tensor source
        call locate_source_moment_tensor(ibool,coord,nspec,npoin,xigll,zigll,x_source(i_source),z_source(i_source), &
             ispec_selected_source(i_source),is_proc_source(i_source),nb_proc_source(i_source),&
             nproc,myrank,xi_source(i_source),gamma_source(i_source),coorg,knods,ngnod,npgeo,ipass)

        ! compute source array for moment-tensor source
        call compute_arrays_source(ispec_selected_source(i_source),xi_source(i_source),gamma_source(i_source),&
             sourcearray(i_source,:,:,:), &
             Mxx(i_source),Mzz(i_source),Mxz(i_source),xix,xiz,gammax,gammaz,xigll,zigll,nspec)

     else if(.not.initialfield) then
        call exit_MPI('incorrect source type')
     endif


     ! locate receivers in the mesh
     call locate_receivers(ibool,coord,nspec,npoin,xigll,zigll,nrec,nrecloc,recloc,which_proc_receiver,nproc,myrank,&
          st_xval,st_zval,ispec_selected_rec, &
          xi_receiver,gamma_receiver,station_name,network_name,x_source(i_source),z_source(i_source),&
          coorg,knods,ngnod,npgeo,ipass, &
          x_final_receiver, z_final_receiver)

  enddo ! do i_source=1,NSOURCE


end subroutine setup_sources_receivers


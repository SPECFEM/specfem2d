
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
!=====================================================================

! for acoustic solver

  subroutine compute_add_sources_acoustic(potential_dot_dot_acoustic,it,i_stage)

  use specfem_par, only: acoustic,nglob_acoustic,&
                         NSOURCES,source_type,source_time_function,&
                         is_proc_source,ispec_selected_source,&
                         hxis_store,hgammas_store,ibool,kappastore
  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: potential_dot_dot_acoustic
  integer :: it,i_stage

  !local variables
  integer :: i_source,i,j,iglob
  double precision :: hlagrange

  do i_source=1,NSOURCES
    ! if this processor core carries the source and the source element is acoustic
    if (is_proc_source(i_source) == 1 .and. acoustic(ispec_selected_source(i_source))) then
      ! collocated force
      ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
      ! the sign is negative because pressure p = - Chi_dot_dot therefore we need
      ! to add minus the source to Chi_dot_dot to get plus the source in pressure
      if( source_type(i_source) == 1 ) then
        ! forward wavefield
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec_selected_source(i_source))
            hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) - &
                                                source_time_function(i_source,it,i_stage)*hlagrange &
                                                !ZN becareful the following line is new added, thus when do comparison
                                                !ZN of the new code with the old code, you will have big difference if you
                                                !ZN do not tune the source
                                                / kappastore(i,j,ispec_selected_source(i_source))
          enddo
        enddo
      ! moment tensor
      else if(source_type(i_source) == 2) then
         call exit_MPI('cannot have moment tensor source in acoustic element')
      endif
    endif ! if this processor core carries the source and the source element is acoustic
  enddo ! do i_source=1,NSOURCES

  end subroutine compute_add_sources_acoustic
!
!=====================================================================
! for acoustic solver for adjoint propagation wave field

  subroutine compute_add_sources_acoustic_adjoint()

  use specfem_par, only: myrank,potential_dot_dot_acoustic,acoustic,NSTEP,it,&
                         nrec,which_proc_receiver,ispec_selected_rec,adj_sourcearrays,&
                         ibool,kappastore
  implicit none
  include "constants.h"

  !local variables
  integer :: irec_local,irec,i,j,iglob

  irec_local = 0
  do irec = 1,nrec
    ! add the source (only if this proc carries the source)
    if( myrank == which_proc_receiver(irec) ) then
      irec_local = irec_local + 1
      if( acoustic(ispec_selected_rec(irec)) ) then
        ! add source array
        do j=1,NGLLZ
          do i=1,NGLLX
            iglob = ibool(i,j,ispec_selected_rec(irec))
            potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + &
                                                adj_sourcearrays(irec_local,NSTEP-it+1,1,i,j) &
                                                !ZN becareful the following line is new added, thus when do comparison
                                                !ZN of the new code with the old code, you will have big difference if you
                                                !ZN do not tune the source
                                                / kappastore(i,j,ispec_selected_rec(irec))
          enddo
        enddo
      endif ! if element acoustic
    endif ! if this processor core carries the adjoint source
  enddo ! irec = 1,nrec

  end subroutine compute_add_sources_acoustic_adjoint
!=====================================================================

  subroutine add_acoustic_forcing_at_rigid_boundary(potential_dot_dot_acoustic)

  use specfem_par, only: nglob_acoustic,nelem_acforcing,codeacforcing,numacforcing,acoustic,&
                         ibool,xix,xiz,jacobian,gammax,gammaz,wxgll,wzgll,&
                         PML_BOUNDARY_CONDITIONS,is_PML

  implicit none
  include "constants.h"

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: potential_dot_dot_acoustic

  !local variables
  integer :: inum,ispec,i,j,iglob
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight,&
                            displ_x,displ_z,displ_n

  ! loop on all the forced edges
  do inum = 1,nelem_acforcing

    ispec = numacforcing(inum)
    if(.not. acoustic(ispec)) cycle ! acoustic spectral element

    !--- left acoustic forcing boundary
    if( codeacforcing(IEDGE4,inum) ) then
      i = 1
      do j = 1,NGLLZ
        iglob = ibool(i,j,ispec)
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)

        ! define displacement components which will force the boundary
        if( PML_BOUNDARY_CONDITIONS ) then
          if( is_PML(ispec) ) then
            displ_x = 0
            displ_z = 0
          else
            call acoustic_forcing_boundary(iglob)
          endif
        else
          call acoustic_forcing_boundary(iglob)
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
      enddo
    endif  !  end of left acoustic forcing boundary

    !--- right acoustic forcing boundary
    if(codeacforcing(IEDGE2,inum)) then
      i = NGLLX
      do j = 1,NGLLZ
        iglob = ibool(i,j,ispec)
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)

        ! define displacement components which will force the boundary
        if( PML_BOUNDARY_CONDITIONS ) then
          if( is_PML(ispec) ) then
            displ_x = 0
            displ_z = 0
          else
            call acoustic_forcing_boundary(iglob)
          endif
        else
          call acoustic_forcing_boundary(iglob)
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
      enddo
    endif  !  end of right acoustic forcing boundary

    !--- bottom acoustic forcing boundary
    if( codeacforcing(IEDGE1,inum) ) then
      j = 1
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        weight = jacobian1D * wxgll(i)

        ! define displacement components which will force the boundary
        if( PML_BOUNDARY_CONDITIONS ) then
          if( is_PML(ispec) ) then
            displ_x = 0
            displ_z = 0
          else
            call acoustic_forcing_boundary(iglob)
          endif
        else
          call acoustic_forcing_boundary(iglob)
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
      enddo
    endif  !  end of bottom acoustic forcing boundary

    !--- top acoustic forcing boundary
    if( codeacforcing(IEDGE3,inum) ) then
      j = NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        weight = jacobian1D * wxgll(i)

        ! define displacement components which will force the boundary
        if( PML_BOUNDARY_CONDITIONS ) then
          if( is_PML(ispec) ) then
            displ_x = 0
            displ_z = 0
          else
            call acoustic_forcing_boundary(iglob)
          endif
        else
          call acoustic_forcing_boundary(iglob)
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
      enddo
    endif  !  end of top acoustic forcing boundary
  enddo

  end subroutine add_acoustic_forcing_at_rigid_boundary

!=====================================================================
! *********************************************************
! ** impose displacement from acoustic forcing at a rigid boundary
! ** force potential_dot_dot_gravito by displacement
! *********************************************************
  subroutine add_acoustic_forcing_at_rigid_boundary_gravitoacoustic()

  use specfem_par, only: nelem_acforcing,codeacforcing,numacforcing, &
                         gravitoacoustic,potential_dot_dot_gravito, &
                         potential_gravitoacoustic,potential_gravito, &
                         it,ibool,xix,xiz,jacobian,gammax,gammaz,wxgll,wzgll,hprime_xx,hprime_zz, &
                         iglobzero,assign_external_model,rhoext,gravityext,Nsqext, &
                         PML_BOUNDARY_CONDITIONS,is_PML

  implicit none
  include "constants.h"

  !local variables
  integer :: inum,ispec,i,j,k,iglob
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight,&
                            tempx1l,hp1,tempx2l,hp2,xixl,xizl,gammaxl,gammazl, &
                            rhol,gravityl,Nsql,displ_x,displ_z,displ_n


  ! loop on all the forced edges
  do inum = 1,nelem_acforcing

    ispec = numacforcing(inum)
    ! gravito spectral element
    if( .not. gravitoacoustic(ispec) ) cycle

    !--- left acoustic forcing boundary
    if( codeacforcing(IEDGE4,inum) ) then
      i = 1
      do j = 1,NGLLZ
        iglob = ibool(i,j,ispec)
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)

        ! define displacement components which will force the boundary
        if( PML_BOUNDARY_CONDITIONS ) then
          if( is_PML(ispec) ) then
            displ_x = 0
            displ_z = 0
          else
            call acoustic_forcing_boundary(iglob)
          endif
        else
          call acoustic_forcing_boundary(iglob)
        endif

        ! compute displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! if external density model
        if( assign_external_model ) then
          rhol = rhoext(i,j,ispec)
          gravityl = gravityext(i,j,ispec)
        endif

        ! impose potential_gravito in order to have z displacement equal to forced value
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if( abs(nz) > TINYVAL ) then
          potential_gravito(iglob) = ( rhol*displ_n -(tempx1l*xizl + tempx2l*gammazl)*nz - &
                                       (tempx1l*xixl + tempx2l*gammaxl)*nx ) / (0._CUSTOM_REAL - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n
      enddo
    endif  !  end of left acoustic forcing boundary

    !--- right acoustic forcing boundary
    if( codeacforcing(IEDGE2,inum) ) then
      i = NGLLX
      do j = 1,NGLLZ
        iglob = ibool(i,j,ispec)
        xgamma = - xiz(i,j,ispec) * jacobian(i,j,ispec)
        zgamma = + xix(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)

        ! define displacement components which will force the boundary
        if( PML_BOUNDARY_CONDITIONS ) then
          if( is_PML(ispec) ) then
            displ_x = 0
            displ_z = 0
          else
            call acoustic_forcing_boundary(iglob)
          endif
        else
          call acoustic_forcing_boundary(iglob)
        endif

        ! compute displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! if external density model
        if( assign_external_model ) then
          rhol = rhoext(i,j,ispec)
          gravityl = gravityext(i,j,ispec)
        endif

        ! impose potential_gravito in order to have z displacement equal to forced value
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if( abs(nz) > TINYVAL ) then
          potential_gravito(iglob) = ( rhol*displ_n - (tempx1l*xizl + tempx2l*gammazl)*nz - &
                                       (tempx1l*xixl + tempx2l*gammaxl)*nx ) / (0._CUSTOM_REAL - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n

      enddo
    endif  !  end of right acoustic forcing boundary

    !--- bottom acoustic forcing boundary
    if( codeacforcing(IEDGE1,inum) ) then
      j = 1
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        weight = jacobian1D * wxgll(i)

        ! define displacement components which will force the boundary
        if( PML_BOUNDARY_CONDITIONS ) then
          if( is_PML(ispec) ) then
            displ_x = 0
            displ_z = 0
          else
            call acoustic_forcing_boundary(iglob)
          endif
        else
          call acoustic_forcing_boundary(iglob)
        endif

        ! compute displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! if external density model
        if( assign_external_model ) then
           rhol = rhoext(i,j,ispec)
           gravityl = gravityext(i,j,ispec)
        endif

        ! impose potential_gravito in order to have z displacement equal to forced value
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if( abs(nz) > TINYVAL) then
          potential_gravito(iglob) = ( rhol*displ_n - (tempx1l*xizl + tempx2l*gammazl)*nz - &
                                       (tempx1l*xixl + tempx2l*gammaxl)*nx ) / (0._CUSTOM_REAL - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n
      enddo
    endif  !  end of bottom acoustic forcing boundary

    !--- top acoustic forcing boundary
    if( codeacforcing(IEDGE3,inum) ) then
      j = NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        xxi = + gammaz(i,j,ispec) * jacobian(i,j,ispec)
        zxi = - gammax(i,j,ispec) * jacobian(i,j,ispec)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        weight = jacobian1D * wxgll(i)

        ! define displacement components which will force the boundary
        if( PML_BOUNDARY_CONDITIONS ) then
          if( is_PML(ispec) ) then
            displ_x = 0
            displ_z = 0
          else
            call acoustic_forcing_boundary(iglob)
          endif
        else
          call acoustic_forcing_boundary(iglob)
        endif

        ! compute z displacement at this point
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + potential_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + potential_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! if external density model
        if( assign_external_model ) then
          rhol = rhoext(i,j,ispec)
          gravityl = gravityext(i,j,ispec)
          Nsql = Nsqext(i,j,ispec)
        endif

        ! impose potential_gravito in order to have z displacement equal to forced value on the boundary
        !!!! Passe deux fois sur le meme iglob
        !!!! Mais vrai pour tous les points partages entre deux elements
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if( abs(nz) > TINYVAL ) then
          potential_gravito(iglob) = ( rhol*displ_n - (tempx1l*xizl + tempx2l*gammazl)*nz - &
                                      (tempx1l*xixl + tempx2l*gammaxl)*nx ) / (0._CUSTOM_REAL - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_gravito(iglob) = potential_dot_dot_gravito(iglob) - rhol*weight*displ_n
      enddo

      !write(*,*) 'ispec detection =',ispec
      !if ((ispec==2000).and.(mod(it,100)==0)) then
      if( (ispec==800) .and. (mod(it,100)==0) ) then
      !if ((ispec==800)) then
        iglobzero=iglob
        write(*,*) ispec,it,Nsql,rhol,displ_n, &
                   maxval(potential_dot_dot_gravito),potential_dot_dot_gravito(iglob), &
                   maxval(potential_gravitoacoustic),potential_gravitoacoustic(iglob), &
                   maxval(potential_gravito),potential_gravito(iglob)
      endif
    endif  !  end of top acoustic forcing boundary
  enddo

  end subroutine add_acoustic_forcing_at_rigid_boundary_gravitoacoustic


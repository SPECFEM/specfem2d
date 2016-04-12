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
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!=====================================================================

! for acoustic solver

  subroutine compute_add_sources_acoustic(minus_pressure_acoustic,it,i_stage)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ

  use specfem_par, only: ispec_is_acoustic,nglob_acoustic,&
                         NSOURCES,source_type,source_time_function,&
                         is_proc_source,ispec_selected_source,&
                         hxis_store,hgammas_store,ibool,kappastore,myrank
  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic),intent(inout) :: minus_pressure_acoustic
  integer,intent(in) :: it,i_stage

  !local variables
  integer :: i_source,i,j,iglob
  double precision :: hlagrange

  do i_source= 1,NSOURCES
    ! if this processor core carries the source and the source element is acoustic
    if (is_proc_source(i_source) == 1 .and. ispec_is_acoustic(ispec_selected_source(i_source))) then
      ! collocated force
      ! beware, for acoustic medium, source is: pressure divided by Kappa of the fluid
      ! the sign is negative because pressure p = - minus_pressure therefore we need
      ! to add minus the source to minus_pressure to get plus the source in pressure
      if (source_type(i_source) == 1) then
        ! forward wavefield
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec_selected_source(i_source))
            hlagrange = hxis_store(i_source,i) * hgammas_store(i_source,j)
            minus_pressure_acoustic(iglob) = minus_pressure_acoustic(iglob) - &
                                                source_time_function(i_source,it,i_stage)*hlagrange &
                                                !ZN becareful the following line is new added, thus when do comparison
                                                !ZN of the new code with the old code, you will have big difference if you
                                                !ZN do not tune the source
                                                / kappastore(i,j,ispec_selected_source(i_source))
          enddo
        enddo
      ! moment tensor
      else if (source_type(i_source) == 2) then
         call exit_MPI(myrank,'cannot have moment tensor source in acoustic element')
      endif
    endif ! if this processor core carries the source and the source element is acoustic
  enddo ! do i_source= 1,NSOURCES

  end subroutine compute_add_sources_acoustic

!
!=====================================================================
!

! for acoustic solver for adjoint propagation wave field

  subroutine compute_add_sources_acoustic_adjoint()

  use constants,only: NGLLX,NGLLZ

  use specfem_par, only: myrank,minus_pressure_acoustic,ispec_is_acoustic,NSTEP,it,&
                         nrec,which_proc_receiver,ispec_selected_rec,adj_sourcearrays,&
                         ibool,kappastore
  implicit none

  !local variables
  integer :: irec_local,irec,i,j,iglob
  integer :: it_tmp

  ! time step index
  it_tmp = NSTEP - it + 1

  irec_local = 0
  do irec = 1,nrec
    ! add the source (only if this proc carries the source)
    if (myrank == which_proc_receiver(irec)) then
      irec_local = irec_local + 1
      if (ispec_is_acoustic(ispec_selected_rec(irec))) then
        ! add source array
        do j = 1,NGLLZ
          do i = 1,NGLLX
            iglob = ibool(i,j,ispec_selected_rec(irec))
            minus_pressure_acoustic(iglob) = minus_pressure_acoustic(iglob) + &
                                                adj_sourcearrays(irec_local,it_tmp,1,i,j) &
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

!
!=====================================================================
!

  subroutine add_acoustic_forcing_at_rigid_boundary(minus_pressure_acoustic)

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob_acoustic,nelem_acforcing,codeacforcing,numacforcing,ispec_is_acoustic,&
                         ibool,xix,xiz,jacobian,gammax,gammaz,wxgll,wzgll

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: minus_pressure_acoustic

  !local variables
  integer :: inum,ispec,i,j,iglob
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight,&
                            displ_x,displ_z,displ_n

  ! loop on all the forced edges
  do inum = 1,nelem_acforcing

    ispec = numacforcing(inum)
    if (.not. ispec_is_acoustic(ispec)) cycle ! acoustic spectral element

    !--- left acoustic forcing boundary
    if (codeacforcing(IEDGE4,inum)) then
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
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
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
        minus_pressure_acoustic(iglob) = minus_pressure_acoustic(iglob) + weight*displ_n
      enddo
    endif  !  end of left acoustic forcing boundary

    !--- right acoustic forcing boundary
    if (codeacforcing(IEDGE2,inum)) then
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
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
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
        minus_pressure_acoustic(iglob) = minus_pressure_acoustic(iglob) + weight*displ_n
      enddo
    endif  !  end of right acoustic forcing boundary

    !--- bottom acoustic forcing boundary
    if (codeacforcing(IEDGE1,inum)) then
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
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
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
        minus_pressure_acoustic(iglob) = minus_pressure_acoustic(iglob) + weight*displ_n
      enddo
    endif  !  end of bottom acoustic forcing boundary

    !--- top acoustic forcing boundary
    if (codeacforcing(IEDGE3,inum)) then
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
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
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
        minus_pressure_acoustic(iglob) = minus_pressure_acoustic(iglob) + weight*displ_n
      enddo
    endif  !  end of top acoustic forcing boundary
  enddo

  end subroutine add_acoustic_forcing_at_rigid_boundary

!
!=====================================================================
!

! *********************************************************
! ** impose displacement from acoustic forcing at a rigid boundary
! ** force minus_pressure_gravito by displacement
! *********************************************************

  subroutine add_acoustic_forcing_at_rigid_boundary_gravitoacoustic()

  use constants,only: CUSTOM_REAL,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4,TINYVAL

  use specfem_par, only: nelem_acforcing,codeacforcing,numacforcing, &
                         ispec_is_gravitoacoustic,minus_pressure_gravito, &
                         minus_int_int_pressure_gravitoacoustic,minus_int_int_pressure_gravito, &
                         it,ibool,xix,xiz,jacobian,gammax,gammaz,wxgll,wzgll,hprime_xx,hprime_zz, &
                         iglobzero,assign_external_model,rhoext,gravityext,Nsqext

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML

  implicit none

  !local variables
  integer :: inum,ispec,i,j,k,iglob
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight,&
                            tempx1l,hp1,tempx2l,hp2,xixl,xizl,gammaxl,gammazl, &
                            rhol,gravityl,Nsql,displ_x,displ_z,displ_n


  ! loop on all the forced edges
  do inum = 1,nelem_acforcing

    ispec = numacforcing(inum)
    ! gravito spectral element
    if (.not. ispec_is_gravitoacoustic(ispec) ) cycle

    !--- left acoustic forcing boundary
    if (codeacforcing(IEDGE4,inum)) then
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
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
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
          tempx1l = tempx1l + minus_int_int_pressure_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + minus_int_int_pressure_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! if external density model
        if (assign_external_model) then
          rhol = rhoext(i,j,ispec)
          gravityl = gravityext(i,j,ispec)
        endif

        ! impose minus_int_int_pressure_gravito in order to have z displacement equal to forced value
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          minus_int_int_pressure_gravito(iglob) = ( rhol*displ_n -(tempx1l*xizl + tempx2l*gammazl)*nz - &
                                       (tempx1l*xixl + tempx2l*gammaxl)*nx ) / (0._CUSTOM_REAL - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        minus_pressure_gravito(iglob) = minus_pressure_gravito(iglob) - rhol*weight*displ_n
      enddo
    endif  !  end of left acoustic forcing boundary

    !--- right acoustic forcing boundary
    if (codeacforcing(IEDGE2,inum)) then
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
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
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
          tempx1l = tempx1l + minus_int_int_pressure_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + minus_int_int_pressure_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! if external density model
        if (assign_external_model) then
          rhol = rhoext(i,j,ispec)
          gravityl = gravityext(i,j,ispec)
        endif

        ! impose minus_int_int_pressure_gravito in order to have z displacement equal to forced value
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          minus_int_int_pressure_gravito(iglob) = ( rhol*displ_n - (tempx1l*xizl + tempx2l*gammazl)*nz - &
                                       (tempx1l*xixl + tempx2l*gammaxl)*nx ) / (0._CUSTOM_REAL - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        minus_pressure_gravito(iglob) = minus_pressure_gravito(iglob) - rhol*weight*displ_n

      enddo
    endif  !  end of right acoustic forcing boundary

    !--- bottom acoustic forcing boundary
    if (codeacforcing(IEDGE1,inum)) then
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
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
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
          tempx1l = tempx1l + minus_int_int_pressure_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + minus_int_int_pressure_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! if external density model
        if (assign_external_model) then
           rhol = rhoext(i,j,ispec)
           gravityl = gravityext(i,j,ispec)
        endif

        ! impose minus_int_int_pressure_gravito in order to have z displacement equal to forced value
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          minus_int_int_pressure_gravito(iglob) = ( rhol*displ_n - (tempx1l*xizl + tempx2l*gammazl)*nz - &
                                       (tempx1l*xixl + tempx2l*gammaxl)*nx ) / (0._CUSTOM_REAL - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        minus_pressure_gravito(iglob) = minus_pressure_gravito(iglob) - rhol*weight*displ_n
      enddo
    endif  !  end of bottom acoustic forcing boundary

    !--- top acoustic forcing boundary
    if (codeacforcing(IEDGE3,inum)) then
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
        if (PML_BOUNDARY_CONDITIONS) then
          if (ispec_is_PML(ispec)) then
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
          tempx1l = tempx1l + minus_int_int_pressure_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + minus_int_int_pressure_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        ! if external density model
        if (assign_external_model) then
          rhol = rhoext(i,j,ispec)
          gravityl = gravityext(i,j,ispec)
          Nsql = Nsqext(i,j,ispec)
        endif

        ! impose minus_int_int_pressure_gravito in order to have z displacement equal to forced value on the boundary
        !!!! Passe deux fois sur le meme iglob
        !!!! Mais vrai pour tous les points partages entre deux elements
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
          minus_int_int_pressure_gravito(iglob) = ( rhol*displ_n - (tempx1l*xizl + tempx2l*gammazl)*nz - &
                                      (tempx1l*xixl + tempx2l*gammaxl)*nx ) / (0._CUSTOM_REAL - gravityl*nz)
        else
          write(*,*) 'STOP : forcing surface element along z',i,j,ispec,iglob,nx,nz
          stop
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        minus_pressure_gravito(iglob) = minus_pressure_gravito(iglob) - rhol*weight*displ_n
      enddo

      ! debugging
      !write(*,*) 'ispec detection =',ispec
      !if ((ispec==2000).and.(mod(it,100)==0)) then
      if ((ispec==800) .and. (mod(it,100)==0)) then
      !if ((ispec==800)) then
        iglobzero=iglob
        write(*,*) ispec,it,Nsql,rhol,displ_n, &
                   maxval(minus_pressure_gravito),minus_pressure_gravito(iglob), &
                   maxval(minus_int_int_pressure_gravitoacoustic),minus_int_int_pressure_gravitoacoustic(iglob), &
                   maxval(minus_int_int_pressure_gravito),minus_int_int_pressure_gravito(iglob)
      endif
    endif  !  end of top acoustic forcing boundary
  enddo

  end subroutine add_acoustic_forcing_at_rigid_boundary_gravitoacoustic


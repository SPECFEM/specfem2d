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

  subroutine add_acoustic_forcing_at_rigid_boundary(potential_dot_dot_acoustic)

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4

  use specfem_par, only: nglob_acoustic,nelem_acforcing,codeacforcing,numacforcing,ispec_is_acoustic, &
                         ibool,xix,xiz,jacobian,gammax,gammaz,wxgll,wzgll

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: potential_dot_dot_acoustic

  !local variables
  integer :: inum,ispec,i,j,iglob
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight, &
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
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
        endif
        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
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
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n

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
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
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
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
        endif

        ! compute dot product
        displ_n = displ_x*nx + displ_z*nz
        potential_dot_dot_acoustic(iglob) = potential_dot_dot_acoustic(iglob) + weight*displ_n
      enddo
    endif  !  end of top acoustic forcing boundary
  enddo

  end subroutine add_acoustic_forcing_at_rigid_boundary

!
!=====================================================================
!

! *********************************************************
! ** impose displacement from acoustic forcing at a rigid boundary
! ** force potential_dot_dot_gravito by displacement
! *********************************************************

  subroutine add_acoustic_forcing_at_rigid_boundary_gravitoacoustic()

  use constants, only: CUSTOM_REAL,NGLLX,NGLLZ,IEDGE1,IEDGE2,IEDGE3,IEDGE4,TINYVAL

  use specfem_par, only: nelem_acforcing,codeacforcing,numacforcing, &
                         ispec_is_gravitoacoustic,potential_dot_dot_gravito, &
                         potential_gravitoacoustic,potential_gravito, &
                         it,ibool,xix,xiz,jacobian,gammax,gammaz,wxgll,wzgll,hprime_xx,hprime_zz, &
                         iglobzero,assign_external_model,rhoext,gravityext,Nsqext

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML

  implicit none

  !local variables
  integer :: inum,ispec,i,j,k,iglob
  real(kind=CUSTOM_REAL) :: xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight, &
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
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
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
        if (assign_external_model) then
          rhol = rhoext(i,j,ispec)
          gravityl = gravityext(i,j,ispec)
        endif

        ! impose potential_gravito in order to have z displacement equal to forced value
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
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
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
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
        if (assign_external_model) then
          rhol = rhoext(i,j,ispec)
          gravityl = gravityext(i,j,ispec)
        endif

        ! impose potential_gravito in order to have z displacement equal to forced value
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
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
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
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
        if (assign_external_model) then
           rhol = rhoext(i,j,ispec)
           gravityl = gravityext(i,j,ispec)
        endif

        ! impose potential_gravito in order to have z displacement equal to forced value
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
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
            displ_x = 0.0d0
            displ_z = 0.0d0
          else
            call acoustic_forcing_boundary(iglob,displ_x,displ_z)
          endif
        else
          call acoustic_forcing_boundary(iglob,displ_x,displ_z)
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
        if (assign_external_model) then
          rhol = rhoext(i,j,ispec)
          gravityl = gravityext(i,j,ispec)
          Nsql = Nsqext(i,j,ispec)
        endif

        ! impose potential_gravito in order to have z displacement equal to forced value on the boundary
        !!!! Passe deux fois sur le meme iglob
        !!!! Mais vrai pour tous les points partages entre deux elements
        iglob = ibool(i,j,ispec)
        displ_n = displ_x*nx + displ_z*nz
        if (abs(nz) > TINYVAL) then
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

      ! debugging
      !write(*,*) 'ispec detection =',ispec
      !if ((ispec==2000) .and. (mod(it,100)==0)) then
      if ((ispec == 800) .and. (mod(it,100) == 0)) then
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


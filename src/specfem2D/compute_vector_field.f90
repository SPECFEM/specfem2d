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

  subroutine compute_vector_whole_medium(field_acoustic,field_gravitoacoustic, &
                                         field_gravito,field_elastic,fields_poroelastic)

! compute Grad(potential) in acoustic elements
! and combine with existing velocity vector field in elastic elements

  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM, &
    nspec,ibool, &
    nglob_acoustic,nglob_elastic,nglob_gravitoacoustic,nglob_poroelastic, &
    vector_field_display

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: field_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob_gravitoacoustic) :: field_gravitoacoustic
  real(kind=CUSTOM_REAL), dimension(nglob_gravitoacoustic) :: field_gravito
  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: field_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic) :: fields_poroelastic
  ! vector field in an element
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLZ) :: vector_field_element

  ! local parameters
  integer :: i,j,ispec,iglob

  ! loop over spectral elements
  do ispec = 1,nspec
    ! computes vector field in this element
    call compute_vector_one_element(field_acoustic,field_gravitoacoustic, &
                                    field_gravito,field_elastic,fields_poroelastic,ispec,vector_field_element)

    ! stores the result on global nodes
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        vector_field_display(:,iglob) = vector_field_element(:,i,j)
      enddo
    enddo
  enddo

  end subroutine compute_vector_whole_medium

!
!=====================================================================
!

  subroutine compute_vector_one_element(field_acoustic,field_gravitoacoustic, &
                                        field_gravito,field_elastic,fields_poroelastic,ispec,vector_field_element)

! compute Grad(potential) if acoustic element or copy existing vector if elastic element

  use specfem_par,only: CUSTOM_REAL,NGLLX,NGLLZ,NDIM, &
    nglob_acoustic,nglob_elastic,nglob_gravitoacoustic,nglob_poroelastic, &
    assign_external_model,density,kmato,gravityext,rhoext, &
    hprimeBar_xx,hprime_xx,hprime_zz, &
    xix,xiz,gammax,gammaz,ibool, &
    ispec_is_elastic,ispec_is_poroelastic,ispec_is_acoustic,ispec_is_gravitoacoustic, &
    AXISYM,is_on_the_axis

  implicit none

  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: field_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob_gravitoacoustic) :: field_gravitoacoustic
  real(kind=CUSTOM_REAL), dimension(nglob_gravitoacoustic) :: field_gravito
  real(kind=CUSTOM_REAL), dimension(3,nglob_elastic) :: field_elastic
  real(kind=CUSTOM_REAL), dimension(NDIM,nglob_poroelastic) :: fields_poroelastic

  integer,intent(in) :: ispec

  ! vector field in an element
  real(kind=CUSTOM_REAL), dimension(3,NGLLX,NGLLZ),intent(out) :: vector_field_element

  ! local variables
  integer i,j,k,iglob
  ! Jacobian matrix and determinant
  double precision :: xixl,xizl,gammaxl,gammazl
  double precision :: gravityl,hp1,hp2
  double precision :: rhol
  double precision :: tempx1l,tempx2l

  ! initializes
  vector_field_element(:,:,:) = 0._CUSTOM_REAL

  ! determines vector field
  if (ispec_is_elastic(ispec)) then
    ! elastic element
    ! simple copy of existing vector if elastic element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        vector_field_element(1,i,j) = field_elastic(1,iglob)
        vector_field_element(2,i,j) = field_elastic(2,iglob)
        vector_field_element(3,i,j) = field_elastic(3,iglob)
      enddo
    enddo

  else if (ispec_is_poroelastic(ispec)) then
    ! poro-elastic element
    ! simple copy of existing vector if poroelastic element
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)
        vector_field_element(1,i,j) = fields_poroelastic(1,iglob)
        vector_field_element(2,i,j) = 0._CUSTOM_REAL
        vector_field_element(3,i,j) = fields_poroelastic(2,iglob)
      enddo
    enddo

  else if (ispec_is_acoustic(ispec)) then
    ! acoustic element

    ! compute gradient of potential to calculate vector if acoustic element
    ! we then need to divide by density because the potential is a potential of (density * displacement)
    rhol = density(1,kmato(ispec))

    ! double loop over GLL points to compute and store gradients
    do j = 1,NGLLZ
      do i = 1,NGLLX
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        if (AXISYM) then
          if (is_on_the_axis(ispec)) then
            do k = 1,NGLLX
              hp1 = hprimeBar_xx(i,k)
              iglob = ibool(k,j,ispec)
              tempx1l = tempx1l + field_acoustic(iglob)*hp1
            enddo
          else
            !AXISYM but not on the axis
            do k = 1,NGLLX
              hp1 = hprime_xx(i,k)
              iglob = ibool(k,j,ispec)
              tempx1l = tempx1l + field_acoustic(iglob)*hp1
            enddo
          endif
        else
          !not AXISYM
          do k = 1,NGLLX
            hp1 = hprime_xx(i,k)
            iglob = ibool(k,j,ispec)
            tempx1l = tempx1l + field_acoustic(iglob)*hp1
          enddo
        endif

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + field_acoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        if (assign_external_model) rhol = rhoext(i,j,ispec)

        ! derivatives of potential
        vector_field_element(1,i,j) = (tempx1l*xixl + tempx2l*gammaxl) / rhol        !u_x
        vector_field_element(2,i,j) = 0._CUSTOM_REAL
        vector_field_element(3,i,j) = (tempx1l*xizl + tempx2l*gammazl) / rhol        !u_z
      enddo
    enddo

  else if (ispec_is_gravitoacoustic(ispec)) then
    ! gravito-acoustic element

    ! compute gradient of potential to calculate vector if gravitoacoustic element
    ! we then need to divide by density because the potential is a potential of (density * displacement)
    rhol = density(1,kmato(ispec))

    ! double loop over GLL points to compute and store gradients
    do j = 1,NGLLZ
      do i = 1,NGLLX
        ! derivative along x
        tempx1l = 0._CUSTOM_REAL
        do k = 1,NGLLX
          hp1 = hprime_xx(i,k)
          iglob = ibool(k,j,ispec)
          tempx1l = tempx1l + field_gravitoacoustic(iglob)*hp1
        enddo

        ! derivative along z
        tempx2l = 0._CUSTOM_REAL
        do k = 1,NGLLZ
          hp2 = hprime_zz(j,k)
          iglob = ibool(i,k,ispec)
          tempx2l = tempx2l + field_gravitoacoustic(iglob)*hp2
        enddo

        xixl = xix(i,j,ispec)
        xizl = xiz(i,j,ispec)
        gammaxl = gammax(i,j,ispec)
        gammazl = gammaz(i,j,ispec)

        if (assign_external_model) then
          rhol = rhoext(i,j,ispec)
          gravityl = gravityext(i,j,ispec)
        endif

        ! derivatives of potential
        vector_field_element(1,i,j) = (tempx1l*xixl + tempx2l*gammaxl) / rhol
        vector_field_element(2,i,j) = 0._CUSTOM_REAL
        vector_field_element(3,i,j) = (tempx1l*xizl + tempx2l*gammazl) / rhol

        ! add the gravito potential along the z component
        iglob = ibool(i,j,ispec)
        ! remove gravito contribution
        ! sign gravito correction
        vector_field_element(3,i,j) = vector_field_element(3,i,j) - (field_gravito(iglob)*gravityl) / rhol

      enddo
    enddo

  endif ! end of test if acoustic or elastic element

  end subroutine compute_vector_one_element


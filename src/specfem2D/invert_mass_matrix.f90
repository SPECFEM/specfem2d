
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
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

  subroutine invert_mass_matrix_init(any_elastic,any_acoustic,any_poroelastic, &
                                rmass_inverse_elastic,nglob_elastic, &
                                rmass_inverse_acoustic,nglob_acoustic, &
                                rmass_s_inverse_poroelastic, &
                                rmass_w_inverse_poroelastic,nglob_poroelastic, &
                                nspec,ibool,kmato,wxgll,wzgll,jacobian, &
                                elastic,poroelastic, &
                                assign_external_model,numat, &
                                density,poroelastcoef,porosity,tortuosity, &
                                vpext,rhoext)

!  builds the global mass matrix

  implicit none
  include 'constants.h'

  logical any_elastic,any_acoustic,any_poroelastic

  ! inverse mass matrices
  integer :: nglob_elastic
  real(kind=CUSTOM_REAL), dimension(nglob_elastic) :: rmass_inverse_elastic

  integer :: nglob_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: rmass_inverse_acoustic

  integer :: nglob_poroelastic
  real(kind=CUSTOM_REAL), dimension(nglob_poroelastic) :: &
    rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic

  integer :: nspec
  integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
  integer, dimension(nspec) :: kmato
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wxgll
  real(kind=CUSTOM_REAL), dimension(NGLLX) :: wzgll
  real(kind=CUSTOM_REAL), dimension(NGLLX,NGLLZ,nspec) :: jacobian

  logical,dimension(nspec) :: elastic,poroelastic

  logical :: assign_external_model
  integer :: numat
  double precision, dimension(2,numat) :: density
  double precision, dimension(4,3,numat) :: poroelastcoef
  double precision, dimension(numat) :: porosity,tortuosity
  double precision, dimension(NGLLX,NGLLX,nspec) :: vpext,rhoext

  ! local parameters
  integer :: ispec,i,j,iglob
  double precision :: rhol,kappal,mul_relaxed,lambdal_relaxed
  double precision :: rhol_s,rhol_f,rhol_bar,phil,tortl

  ! initializes mass matrix
  if(any_elastic) rmass_inverse_elastic(:) = 0._CUSTOM_REAL
  if(any_poroelastic) rmass_s_inverse_poroelastic(:) = 0._CUSTOM_REAL
  if(any_poroelastic) rmass_w_inverse_poroelastic(:) = 0._CUSTOM_REAL
  if(any_acoustic) rmass_inverse_acoustic(:) = 0._CUSTOM_REAL

  do ispec = 1,nspec
    do j = 1,NGLLZ
      do i = 1,NGLLX
        iglob = ibool(i,j,ispec)

        ! if external density model (elastic or acoustic)
        if(assign_external_model) then
          rhol = rhoext(i,j,ispec)
          kappal = rhol * vpext(i,j,ispec)**2
        else
          rhol = density(1,kmato(ispec))
          lambdal_relaxed = poroelastcoef(1,1,kmato(ispec))
          mul_relaxed = poroelastcoef(2,1,kmato(ispec))
          kappal = lambdal_relaxed + 2.d0/3.d0*mul_relaxed
        endif

        if( poroelastic(ispec) ) then

          ! material is poroelastic

          rhol_s = density(1,kmato(ispec))
          rhol_f = density(2,kmato(ispec))
          phil = porosity(kmato(ispec))
          tortl = tortuosity(kmato(ispec))
          rhol_bar = (1.d0-phil)*rhol_s + phil*rhol_f

          ! for the solid mass matrix
          rmass_s_inverse_poroelastic(iglob) = rmass_s_inverse_poroelastic(iglob)  &
                  + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)*(rhol_bar - phil*rhol_f/tortl)
          ! for the fluid mass matrix
          rmass_w_inverse_poroelastic(iglob) = rmass_w_inverse_poroelastic(iglob) &
                  + wxgll(i)*wzgll(j)*jacobian(i,j,ispec)*(rhol_bar*rhol_f*tortl  &
                  - phil*rhol_f*rhol_f)/(rhol_bar*phil)

        elseif( elastic(ispec) ) then

          ! for elastic medium

          rmass_inverse_elastic(iglob) = rmass_inverse_elastic(iglob)  &
                  + wxgll(i)*wzgll(j)*rhol*jacobian(i,j,ispec)

        else

          ! for acoustic medium

          rmass_inverse_acoustic(iglob) = rmass_inverse_acoustic(iglob) &
                  + wxgll(i)*wzgll(j)*jacobian(i,j,ispec) / kappal

        endif

      enddo
    enddo
  enddo ! do ispec = 1,nspec

  end subroutine invert_mass_matrix_init
!
!-------------------------------------------------------------------------------------------------
!

  subroutine invert_mass_matrix(any_elastic,any_acoustic,any_poroelastic, &
                                rmass_inverse_elastic,nglob_elastic, &
                                rmass_inverse_acoustic,nglob_acoustic, &
                                rmass_s_inverse_poroelastic, &
                                rmass_w_inverse_poroelastic,nglob_poroelastic)

! inverts the global mass matrix

  implicit none
  include 'constants.h'

  logical any_elastic,any_acoustic,any_poroelastic

! inverse mass matrices
  integer :: nglob_elastic
  real(kind=CUSTOM_REAL), dimension(nglob_elastic) :: rmass_inverse_elastic

  integer :: nglob_acoustic
  real(kind=CUSTOM_REAL), dimension(nglob_acoustic) :: rmass_inverse_acoustic

  integer :: nglob_poroelastic
  real(kind=CUSTOM_REAL), dimension(nglob_poroelastic) :: &
    rmass_s_inverse_poroelastic,rmass_w_inverse_poroelastic


! fill mass matrix with fictitious non-zero values to make sure it can be inverted globally
  if(any_elastic) &
    where(rmass_inverse_elastic <= 0._CUSTOM_REAL) rmass_inverse_elastic = 1._CUSTOM_REAL
  if(any_poroelastic) &
    where(rmass_s_inverse_poroelastic <= 0._CUSTOM_REAL) rmass_s_inverse_poroelastic = 1._CUSTOM_REAL
  if(any_poroelastic) &
    where(rmass_w_inverse_poroelastic <= 0._CUSTOM_REAL) rmass_w_inverse_poroelastic = 1._CUSTOM_REAL
  if(any_acoustic) &
    where(rmass_inverse_acoustic <= 0._CUSTOM_REAL) rmass_inverse_acoustic = 1._CUSTOM_REAL

! compute the inverse of the mass matrix
  if(any_elastic) &
    rmass_inverse_elastic(:) = 1._CUSTOM_REAL / rmass_inverse_elastic(:)
  if(any_poroelastic) &
    rmass_s_inverse_poroelastic(:) = 1._CUSTOM_REAL / rmass_s_inverse_poroelastic(:)
  if(any_poroelastic) &
    rmass_w_inverse_poroelastic(:) = 1._CUSTOM_REAL / rmass_w_inverse_poroelastic(:)
  if(any_acoustic) &
    rmass_inverse_acoustic(:) = 1._CUSTOM_REAL / rmass_inverse_acoustic(:)

  end subroutine invert_mass_matrix

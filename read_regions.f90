
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

  subroutine read_regions(nbregion,nb_materials,icodemat,cp,cs, &
                          rho_s,Qp,Qs,aniso3,aniso4,aniso5,aniso6,aniso7,aniso8, &
                          nelmnts,num_material,nxread,nzread)

! reads in material definitions in DATA/Par_file

  implicit none
  include "constants.h"

  integer :: nbregion,nb_materials
  integer, dimension(nb_materials) :: icodemat
  double precision, dimension(nb_materials) :: rho_s,cp,cs, &
    aniso3,aniso4,aniso5,aniso6,aniso7,aniso8,Qp,Qs

  integer :: nelmnts
  integer,dimension(nelmnts) :: num_material
  integer :: nxread,nzread

  ! local parameters
  integer :: iregion,ixdebregion,ixfinregion,izdebregion,izfinregion,imaterial_number
  integer :: i,j
  double precision :: vpregion,vsregion,poisson_ratio

  ! read the material numbers for each region
  call read_value_integer(IIN,IGNORE_JUNK,nbregion)

  if(nbregion <= 0) stop 'Negative number of regions not allowed!'

  print *
  print *, 'Nb of regions in the mesh = ',nbregion
  print *

  do iregion = 1,nbregion

    call read_region_coordinates(IIN,DONT_IGNORE_JUNK,ixdebregion,ixfinregion, &
                                izdebregion,izfinregion,imaterial_number)

    if(imaterial_number < 1) stop 'Negative material number not allowed!'
    if(ixdebregion < 1) stop 'Left coordinate of region negative!'
    if(ixfinregion > nxread) stop 'Right coordinate of region too high!'
    if(izdebregion < 1) stop 'Bottom coordinate of region negative!'
    if(izfinregion > nzread) stop 'Top coordinate of region too high!'

    print *,'Region ',iregion
    print *,'IX from ',ixdebregion,' to ',ixfinregion
    print *,'IZ from ',izdebregion,' to ',izfinregion

    if(icodemat(imaterial_number) /= ANISOTROPIC_MATERIAL .and. icodemat(imaterial_number) /= POROELASTIC_MATERIAL) then

       ! isotropic material
       vpregion = cp(imaterial_number)
       vsregion = cs(imaterial_number)
       print *,'Material # ',imaterial_number,' isotropic'
       if(vsregion < TINYVAL) then
          print *,'Material is fluid'
       else
          print *,'Material is solid'
       endif
       print *,'vp = ',vpregion
       print *,'vs = ',vsregion
       print *,'rho = ',rho_s(imaterial_number)
       poisson_ratio = 0.5d0*(vpregion*vpregion-2.d0*vsregion*vsregion) / (vpregion*vpregion-vsregion*vsregion)
       print *,'Poisson''s ratio = ',poisson_ratio
       if(poisson_ratio <= -1.00001d0 .or. poisson_ratio >= 0.50001d0) stop 'incorrect value of Poisson''s ratio'
       print *,'Qp = ',Qp(imaterial_number)
       print *,'Qs = ',Qs(imaterial_number)
    elseif(icodemat(imaterial_number) == POROELASTIC_MATERIAL) then

       ! poroelastic material
       print *,'Material # ',imaterial_number,' isotropic'
       print *,'Material is poroelastic'
    else

       ! anisotropic material
       print *,'Material # ',imaterial_number,' anisotropic'
       print *,'cp = ',cp(imaterial_number)
       print *,'cs = ',cs(imaterial_number)
       print *,'c11 = ',aniso3(imaterial_number)
       print *,'c13 = ',aniso4(imaterial_number)
       print *,'c15 = ',aniso5(imaterial_number)
       print *,'c33 = ',aniso6(imaterial_number)
       print *,'c35 = ',aniso7(imaterial_number)
       print *,'c55 = ',aniso8(imaterial_number)
       print *,'rho = ',rho_s(imaterial_number)
       print *,'Qp = ',Qp(imaterial_number)
       print *,'Qs = ',Qs(imaterial_number)
    endif
    print *,' -----'

    ! store density and velocity model
    do i = ixdebregion,ixfinregion
       do j = izdebregion,izfinregion
          num_material((j-1)*nxread+i) = imaterial_number
       enddo
    enddo

  enddo

  if(minval(num_material) <= 0) stop 'Velocity model not entirely set...'

  end subroutine read_regions

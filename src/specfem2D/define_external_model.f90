
!========================================================================
!
!                   S P E C F E M 2 D  Version 6 . 2
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


  subroutine define_external_model(x,y,iflag_element,myrank,rho,vp,vs,Qp_attenuation,&
       Qs_attenuation,c11,c13,c15,c33,c35,c55 )

  implicit none

  include "constants.h"

! user can modify this routine to assign any different external Earth model (rho, vp, vs)
! based on the x and y coordinates of that grid point and the flag of the region it belongs to

  integer, intent(in) :: iflag_element,myrank

  double precision, intent(in) :: x,y

  double precision, intent(out) :: rho,vp,vs
  double precision, intent(out) :: Qp_attenuation,Qs_attenuation
  double precision, intent(out) :: c11,c15,c13,c33,c35,c55

! dummy routine here, just to demonstrate how the model can be assigned
   if(myrank == 0 .and. iflag_element == 1 .or. x < 1700.d0 .or. y >= 2300.d0) then
     rho = 2000.d0
     vp = 3000.d0
     vs = vp / sqrt(3.d0)
     Qp_attenuation = 0
     Qs_attenuation = 0
     c11 = 169.d9
     c13 = 122.d9
     c15 = 0.d0
     c33 = c11
     c35 = 0.d0
     c55 = 75.3d9
   else
     rho = 2500.d0
     vp = 3600.d0
     vs = vp / 2.d0
     Qp_attenuation = 60
     Qs_attenuation = 60
     c11 = 0.d0
     c13 = 0.d0
     c15 = 0.d0
     c33 = 0.d0
     c35 = 0.d0
     c55 = 0.d0
   endif

  end subroutine define_external_model

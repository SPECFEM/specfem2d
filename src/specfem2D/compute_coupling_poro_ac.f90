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
! for poro solver

 subroutine compute_coupling_poro_ac()

  use specfem_par, only: num_fluid_poro_edges,ibool,wxgll,wzgll,xix,xiz,&
                         gammax,gammaz,jacobian,ivalue,jvalue,ivalue_inverse,jvalue_inverse,&
                         fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge, &
                         fluid_poro_poroelastic_ispec,fluid_poro_poroelastic_iedge, &
                         porosity,tortuosity,density,kmato, &
                         potential_dot_dot_acoustic,b_potential_dot_dot_acoustic,potential_acoustic_adj_coupling, &
                         accels_poroelastic,accelw_poroelastic,b_accels_poroelastic,b_accelw_poroelastic, &
                         SIMULATION_TYPE

  implicit none
  include "constants.h"

  !local variable
  integer :: inum,ispec_acoustic,iedge_acoustic,ispec_poroelastic,iedge_poroelastic, &
             ipoin1D,i,j,ii2,jj2,iglob
  real(kind=CUSTOM_REAL) :: pressure,b_pressure,&
                            xxi,zxi,xgamma,zgamma,jacobian1D,nx,nz,weight
  double precision :: phil,tortl,rhol_f,rhol_s,rhol_bar

  ! loop on all the coupling edges
  do inum = 1,num_fluid_poro_edges
    ! get the edge of the acoustic element
    ispec_acoustic = fluid_poro_acoustic_ispec(inum)
    iedge_acoustic = fluid_poro_acoustic_iedge(inum)
    ! get the corresponding edge of the poroelastic element
    ispec_poroelastic = fluid_poro_poroelastic_ispec(inum)
    iedge_poroelastic = fluid_poro_poroelastic_iedge(inum)

    ! implement 1D coupling along the edge
    do ipoin1D = 1,NGLLX
      ! get point values for the acoustic side, which matches our side in the inverse direction
      i = ivalue_inverse(ipoin1D,iedge_acoustic)
      j = jvalue_inverse(ipoin1D,iedge_acoustic)
      iglob = ibool(i,j,ispec_acoustic)

      ! get poroelastic parameters
      phil = porosity(kmato(ispec_poroelastic))
      tortl = tortuosity(kmato(ispec_poroelastic))
      rhol_f = density(2,kmato(ispec_poroelastic))
      rhol_s = density(1,kmato(ispec_poroelastic))
      rhol_bar = (1._CUSTOM_REAL-phil)*rhol_s + phil*rhol_f

      ! compute pressure on the fluid/porous medium edge
      pressure = - potential_dot_dot_acoustic(iglob)
      if(SIMULATION_TYPE == 3) then
        b_pressure = - b_potential_dot_dot_acoustic(iglob)
        ! new definition of adjoint displacement and adjoint potential
        pressure = potential_acoustic_adj_coupling(iglob)
      endif

      ! get point values for the poroelastic side
      ii2 = ivalue(ipoin1D,iedge_poroelastic)
      jj2 = jvalue(ipoin1D,iedge_poroelastic)
      iglob = ibool(ii2,jj2,ispec_poroelastic)

      ! compute the 1D Jacobian and the normal to the edge: for their expression see for instance
      ! O. C. Zienkiewicz and R. L. Taylor, The Finite Element Method for Solid and Structural Mechanics,
      ! Sixth Edition, electronic version, www.amazon.com, p. 204 and Figure 7.7(a),
      ! or Y. K. Cheung, S. H. Lo and A. Y. T. Leung, Finite Element Implementation,
      ! Blackwell Science, page 110, equation (4.60).
      if( iedge_acoustic == ITOP ) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = - zxi / jacobian1D
        nz = + xxi / jacobian1D
        weight = jacobian1D * wxgll(i)
      else if( iedge_acoustic == IBOTTOM ) then
        xxi = + gammaz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zxi = - gammax(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xxi**2 + zxi**2)
        nx = + zxi / jacobian1D
        nz = - xxi / jacobian1D
        weight = jacobian1D * wxgll(i)
      else if( iedge_acoustic ==ILEFT ) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = - zgamma / jacobian1D
        nz = + xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)
      else if( iedge_acoustic ==IRIGHT ) then
        xgamma = - xiz(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        zgamma = + xix(i,j,ispec_acoustic) * jacobian(i,j,ispec_acoustic)
        jacobian1D = sqrt(xgamma**2 + zgamma**2)
        nx = + zgamma / jacobian1D
        nz = - xgamma / jacobian1D
        weight = jacobian1D * wzgll(j)
      endif

      ! contribution to the solid phase
      accels_poroelastic(1,iglob) = accels_poroelastic(1,iglob) + weight*nx*pressure*(1._CUSTOM_REAL-phil/tortl)
      accels_poroelastic(2,iglob) = accels_poroelastic(2,iglob) + weight*nz*pressure*(1._CUSTOM_REAL-phil/tortl)
      ! contribution to the fluid phase
      accelw_poroelastic(1,iglob) = accelw_poroelastic(1,iglob) + weight*nx*pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
      accelw_poroelastic(2,iglob) = accelw_poroelastic(2,iglob) + weight*nz*pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)

      if( SIMULATION_TYPE == 3 ) then
        ! contribution to the solid phase
        b_accels_poroelastic(1,iglob) = b_accels_poroelastic(1,iglob) + weight*nx*b_pressure*(1._CUSTOM_REAL-phil/tortl)
        b_accels_poroelastic(2,iglob) = b_accels_poroelastic(2,iglob) + weight*nz*b_pressure*(1._CUSTOM_REAL-phil/tortl)
        ! contribution to the fluid phase
        b_accelw_poroelastic(1,iglob) = b_accelw_poroelastic(1,iglob) + weight*nx*b_pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
        b_accelw_poroelastic(2,iglob) = b_accelw_poroelastic(2,iglob) + weight*nz*b_pressure*(1._CUSTOM_REAL-rhol_f/rhol_bar)
      endif
    enddo
  enddo

 end subroutine compute_coupling_poro_ac

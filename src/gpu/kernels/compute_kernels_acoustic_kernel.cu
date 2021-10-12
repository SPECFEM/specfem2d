/*
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
! the Free Software Foundation; either version 3 of the License, or
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
!========================================================================
*/


__global__ void compute_kernels_acoustic_kernel(int* ispec_is_acoustic,
                                                int* d_ibool,
                                                realw* rhostore,
                                                realw* kappastore,
                                                realw* d_hprime_xx,
                                                realw* d_xix,realw* d_xiz,
                                                realw* d_gammax,realw* d_gammaz,
                                                realw* potential_acoustic,
                                                realw* potential_dot_dot_acoustic,
                                                realw* b_potential_acoustic,
                                                realw* rho_ac_kl,
                                                realw* kappa_ac_kl,
                                                int NSPEC_AB,
                                                realw deltat) {

  int ispec = blockIdx.x + blockIdx.y*gridDim.x;
  int ij = threadIdx.x;

  // local and global indices
  int ij_ispec = ij + NGLL2*ispec;
  int ij_ispec_padded = ij + NGLL2_PADDED*ispec;
  int iglob;

  // shared memory between all threads within this block
  __shared__ realw scalar_field_displ[NGLL2];
  __shared__ realw scalar_field_accel[NGLL2];

  int active = 0;

  // calcul the displacement by computing the gradient of potential / rho
  // and calcul the acceleration by computing the gradient of potential_dot_dot / rho
  //
  // note: we use the displ potential_acoustic for the adjoint wavefield, and not the accel potential_dot_dot_acoustic,
  //       to match the rho-kernel expressions from Luo et al. (2013)

  // handles case when there is 1 extra block (due to rectangular grid)
  if (ispec < NSPEC_AB) {
    // acoustic elements only
    if (ispec_is_acoustic[ispec]) {
      active = 1;

      // copy field values
      iglob = d_ibool[ij_ispec_padded] - 1;
      scalar_field_displ[ij] = b_potential_acoustic[iglob];
      // old: scalar_field_accel[ij] = potential_dot_dot_acoustic[iglob];
      scalar_field_accel[ij] = potential_acoustic[iglob];
    }
  }

  // synchronizes threads
  __syncthreads();

  if (active) {
    realw accel_elm[2];
    realw b_displ_elm[2];
    realw rhol,kappal;

    // gets material parameter
    rhol = rhostore[ij_ispec_padded];

    // displacement vector from backward field
    compute_gradient_kernel(ij,ispec,scalar_field_displ,b_displ_elm,
                            d_hprime_xx,
                            d_xix,d_xiz,d_gammax,d_gammaz,
                            rhol);

    // acceleration vector
    compute_gradient_kernel(ij,ispec,scalar_field_accel,accel_elm,
                            d_hprime_xx,
                            d_xix,d_xiz,d_gammax,d_gammaz,
                            rhol);

    // acoustic kernel integration

    // old expression (from elastic kernels)
    //rho_ac_kl(i,j,ispec) = rho_ac_kl(i,j,ispec) - rhol  * &
    //      dot_product(accel_loc(:),b_displ_loc(:)) * deltat
    //kappa_ac_kl(i,j,ispec) = kappa_ac_kl(i,j,ispec) - kappal * &
    //      potential_dot_dot_acoustic(iglob) / kappal * &
    //      b_potential_dot_dot_acoustic(iglob) / kappal
    //      * deltat

    // new expression (from PDE-constrained optimization, coupling terms changed as well)
    // note: from Luo et al. (2013), kernels are given for absolute perturbations d(1/rho) and d(1/kappa)
    //       and only change from Peter et al. (2011) for acoustic kernels:
    //         K_kappa = - int_0^T [ phi^adj \partial_t^2 phi ] dt     see (A-27)
    //         K_rho   = - int_0^T [ grad(phi^adj) * grad(phi) ] dt        (A-28)
    //
    //       since we output relative perturbation kernels for elastic domains, we make here also use of the variation
    //         d(1/rho) = - 1/rho^2 d(rho) = - 1/rho d(rho)/rho = - 1/rho dln(rho)
    //
    //       to obtain relative kernels, we start from eq. (A-24)
    //         dChi = int_V [ d(1/kappa) K_kappa + d(1/rho) K_rho ] d^3x (+ elastic kernels)
    //              = int_V [ (-1/kappa K_kappa) dln(kappa) + (- 1/rho K_rho) dln(rho)
    //
    //       and see that relative perturbation kernels are given by
    //          \tilde{K}_kappa = - 1/kappa K_kappa
    //                          = + 1/kappa int_0^T [ phi^adj \partial_t^2 phi ] dt
    //                          = + 1/kappa int_0^T [ \partial_t^2 phi^adj phi ] dt              (equivalent)
    //
    //          \tilde{K}_rho   = - 1/rho   K_rho
    //                          = + 1/rho   int_0^T [ grad(phi^adj) * grad(phi) ] dt
    //                          = + rho     int_0^T [ 1/rho grad(phi^adj) * 1/rho grad(phi) ] dt   (equivalent)
    //
    // density kernel
    rho_ac_kl[ij_ispec] +=  rhol * (accel_elm[0]*b_displ_elm[0] + accel_elm[1]*b_displ_elm[1] ) * deltat;

    // bulk modulus kernel
    kappal = kappastore[ij_ispec];
    kappa_ac_kl[ij_ispec] += potential_dot_dot_acoustic[iglob] * b_potential_acoustic[iglob] * deltat / kappal ;

  } // active
}


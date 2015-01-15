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

*/

#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"


/* ----------------------------------------------------------------------------------------------- */

// ASSEMBLY - mpi data transfer between CPU-GPU

/* ----------------------------------------------------------------------------------------------- */

// prepares a device array with with all inter-element edge-nodes -- this
// is followed by a memcpy and MPI operations
__global__ void prepare_boundary_potential_on_device(realw* d_potential_dot_dot_acoustic,
                                                     realw* d_send_potential_dot_dot_buffer,
                                                     const int ninterface_ac,
                                                     const int max_nibool_interfaces_ext_mesh,
                                                     const int* d_nibool_interfaces_ext_mesh,
                                                     const int* d_ibool_interfaces_ext_mesh,
                                                     const int* inum_inter_acoustic) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int ientry,iglob,num_int;

  for(int iinterface=0; iinterface < ninterface_ac; iinterface++) {

   num_int=inum_inter_acoustic[iinterface]-1;

    if(id<d_nibool_interfaces_ext_mesh[num_int]) {

      // entry in interface array
      ientry = id + max_nibool_interfaces_ext_mesh*num_int;
      // global index in wavefield
      iglob = d_ibool_interfaces_ext_mesh[ientry] - 1;

      d_send_potential_dot_dot_buffer[ientry] = d_potential_dot_dot_acoustic[iglob];
    }
  }

}


/* ----------------------------------------------------------------------------------------------- */

// prepares and transfers the inter-element edge-nodes to the host to be MPI'd
extern "C"
void FC_FUNC_(transfer_boun_pot_from_device,
              TRANSFER_BOUN_POT_FROM_DEVICE)(long* Mesh_pointer,
                                             realw* send_potential_dot_dot_buffer,
                                             const int* FORWARD_OR_ADJOINT){

TRACE("transfer_boun_pot_from_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // checks if anything to do
  if( mp->size_mpi_buffer_potential > 0 ){

    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)(mp->max_nibool_interfaces_ext_mesh))/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);

    if(*FORWARD_OR_ADJOINT == 1) {

     prepare_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                                   mp->d_send_potential_dot_dot_buffer,
                                                                                   mp->ninterface_acoustic,
                                                                                   mp->max_nibool_interfaces_ext_mesh,
                                                                                   mp->d_nibool_interfaces_ext_mesh,
                                                                                   mp->d_ibool_interfaces_ext_mesh,
                                                                                   mp->d_inum_interfaces_acoustic);

      // synchronizes
      //synchronize_cuda();
      // explicitly waits until previous compute stream finishes
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      print_CUDA_error_if_any(cudaStreamSynchronize(mp->compute_stream),95);

      print_CUDA_error_if_any(cudaMemcpy(send_potential_dot_dot_buffer,mp->d_send_potential_dot_dot_buffer,
                                         mp->size_mpi_buffer_potential*sizeof(realw),cudaMemcpyDeviceToHost),98000);
    }
    else if(*FORWARD_OR_ADJOINT == 3) {
      // backward/reconstructed wavefield buffer
      prepare_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_dot_dot_acoustic,
                                                                                   mp->d_b_send_potential_dot_dot_buffer,
                                                                                   mp->ninterface_acoustic,
                                                                                   mp->max_nibool_interfaces_ext_mesh,
                                                                                   mp->d_nibool_interfaces_ext_mesh,
                                                                                   mp->d_ibool_interfaces_ext_mesh,
                                                                                   mp->d_inum_interfaces_acoustic);

      // synchronizes
      //synchronize_cuda();
      // explicitly waits until previous compute stream finishes
      // (cudaMemcpy implicitly synchronizes all other cuda operations)
      cudaStreamSynchronize(mp->compute_stream);

      print_CUDA_error_if_any(cudaMemcpy(send_potential_dot_dot_buffer,mp->d_b_send_potential_dot_dot_buffer,
                                         mp->size_mpi_buffer_potential*sizeof(realw),cudaMemcpyDeviceToHost),98001);
    }
  }

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("after prepare_boundary_potential_on_device");
#endif


  // finish timing of kernel+memcpy
  // cudaEventRecord( stop, 0 );
  // cudaEventSynchronize( stop );
  // cudaEventElapsedTime( &time, start, stop );
  // cudaEventDestroy( start );
  // cudaEventDestroy( stop );
  // printf("boundary xfer d->h Time: %f ms\n",time);
#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_boun_pot_from_device");
#endif
}

/* ----------------------------------------------------------------------------------------------- */

// Assembly

/* ----------------------------------------------------------------------------------------------- */


__global__ void assemble_boundary_potential_on_device(realw* d_potential_dot_dot_acoustic,
                                                      realw* d_send_potential_dot_dot_buffer,
                                                      const int ninterface_ac,
                                                      const int max_nibool_interfaces_ext_mesh,
                                                      const int* d_nibool_interfaces_ext_mesh,
                                                      const int* d_ibool_interfaces_ext_mesh,
                                                      const int* inum_inter_acoustic) {

  int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
  int ientry,iglob,num_int;

  for( int iinterface=0; iinterface < ninterface_ac; iinterface++) {
   num_int=inum_inter_acoustic[iinterface]-1;


    if(id<d_nibool_interfaces_ext_mesh[num_int]) {

      // entry in interface array
      ientry = id + max_nibool_interfaces_ext_mesh*num_int;
      // global index in wavefield
      iglob = d_ibool_interfaces_ext_mesh[ientry] - 1;

      // for testing atomic operations against not atomic operations (0.1ms vs. 0.04 ms)
      // d_potential_dot_dot_acoustic[3*(d_ibool_interfaces_ext_mesh[id+max_nibool_interfaces_ext_mesh*iinterface]-1)] +=
      // d_send_potential_dot_dot_buffer[3*(id + max_nibool_interfaces_ext_mesh*iinterface)];
    atomicAdd(&d_potential_dot_dot_acoustic[iglob],d_send_potential_dot_dot_buffer[ientry]);
    }
  }
  // ! This step is done via previous function transfer_and_assemble...
  // ! do iinterface = 1, num_interfaces_ext_mesh
  // !   do ipoin = 1, nibool_interfaces_ext_mesh(iinterface)
  // !     array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) = &
  // !          array_val(:,ibool_interfaces_ext_mesh(ipoin,iinterface)) + buffer_recv_vector_ext_mesh(:,ipoin,iinterface)
  // !   enddo
  // ! enddo
}


/* ----------------------------------------------------------------------------------------------- */

extern "C"
void FC_FUNC_(transfer_asmbl_pot_to_device,
              TRANSFER_ASMBL_POT_TO_DEVICE)(long* Mesh_pointer,
                                            realw* buffer_recv_scalar_ext_mesh,
                                            const int* FORWARD_OR_ADJOINT) {

TRACE("transfer_asmbl_pot_to_device");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  // Cuda timing
  //cudaEvent_t start, stop;
  //start_timing_cuda(&start,&stop);

  // checks if anything to do
  if( mp->size_mpi_buffer_potential > 0 ){


    // assembles on GPU
    int blocksize = BLOCKSIZE_TRANSFER;
    int size_padded = ((int)ceil(((double)mp->max_nibool_interfaces_ext_mesh)/((double)blocksize)))*blocksize;

    int num_blocks_x, num_blocks_y;
    get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

    dim3 grid(num_blocks_x,num_blocks_y);
    dim3 threads(blocksize,1,1);


    // synchronizes
    synchronize_cuda();

    if(*FORWARD_OR_ADJOINT == 1) {
      // copies buffer onto GPU
      print_CUDA_error_if_any(cudaMemcpy(mp->d_send_potential_dot_dot_buffer, buffer_recv_scalar_ext_mesh,
                                         mp->size_mpi_buffer_potential*sizeof(realw), cudaMemcpyHostToDevice),98010);

      //assemble forward field
      assemble_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,
                                                                                    mp->d_send_potential_dot_dot_buffer,
                                                                                    mp->ninterface_acoustic,
                                                                                    mp->max_nibool_interfaces_ext_mesh,
                                                                                    mp->d_nibool_interfaces_ext_mesh,
                                                                                    mp->d_ibool_interfaces_ext_mesh,
                                                                                   mp->d_inum_interfaces_acoustic);


    }
    else if(*FORWARD_OR_ADJOINT == 3) {
      // copies buffer onto GPU
      print_CUDA_error_if_any(cudaMemcpy(mp->d_b_send_potential_dot_dot_buffer, buffer_recv_scalar_ext_mesh,
                                         mp->size_mpi_buffer_potential*sizeof(realw), cudaMemcpyHostToDevice),98011);

      //assemble reconstructed/backward field
      assemble_boundary_potential_on_device<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_dot_dot_acoustic,
                                                                                    mp->d_b_send_potential_dot_dot_buffer,
                                                                                    mp->ninterface_acoustic,
                                                                                    mp->max_nibool_interfaces_ext_mesh,
                                                                                    mp->d_nibool_interfaces_ext_mesh,
                                                                                    mp->d_ibool_interfaces_ext_mesh,
                                                                                   mp->d_inum_interfaces_acoustic);
    }
  }

  // Cuda timing
  //stop_timing_cuda(&start,&stop,"assemble_boundary_potential_on_device");

#ifdef ENABLE_VERY_SLOW_ERROR_CHECKING
  exit_on_cuda_error("transfer_asmbl_pot_to_device");
#endif
}


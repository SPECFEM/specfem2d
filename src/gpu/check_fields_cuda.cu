/*
!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
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

#include "mesh_constants_cuda.h"

/* ----------------------------------------------------------------------------------------------- */

// Helper functions

/* ----------------------------------------------------------------------------------------------- */

void cudaMemoryTest(int posId)
{
    // call this function at different places to locate a bug
    const unsigned int N = 1048576;
    const unsigned int bytes = N * sizeof(int);
    int *h_a = (int*)malloc(bytes);
    int *d_a;
    cudaSafeCall(cudaMalloc((int**)&d_a, bytes), posId);
    memset(h_a, 0, bytes);
    cudaSafeCall(cudaMemcpy(d_a, h_a, bytes, cudaMemcpyHostToDevice), posId);
    cudaSafeCall(cudaMemcpy(h_a, d_a, bytes, cudaMemcpyDeviceToHost), posId);
}

/* ----------------------------------------------------------------------------------------------- */

/*
 __global__ void check_field(int* ibool,int* nibool, int max_nibool,int num_interfaces_ext_mesh)
 {


 int id = threadIdx.x + blockIdx.x*blockDim.x + blockIdx.y*gridDim.x*blockDim.x;
 int ientry,iglob;

 for(int iinterface=0; iinterface < num_interfaces_ext_mesh; iinterface++) {
 if (id<nibool[iinterface]) {

 // entry in interface array
 ientry = id + max_nibool*iinterface;
 // global index in wavefield
 iglob = ibool[ientry] - 1;

 cuPrintf("valeurs de iglob %d, de l'indice d'entree %d,du nombre de points dans l'interface %d : %d\n", iglob,ientry,iinterface,nibool[iinterface]);
 }}
 }*/


/* ----------------------------------------------------------------------------------------------- */

// GPU device memory functions

/* ----------------------------------------------------------------------------------------------- */

void get_free_memory(double* free_db, double* used_db, double* total_db) {

  TRACE("get_free_memory");

  // gets memory usage in byte
  size_t free_byte ;
  size_t total_byte ;
  cudaError_t cuda_status = cudaMemGetInfo( &free_byte, &total_byte ) ;
  if (cudaSuccess != cuda_status) {
    printf("Error: cudaMemGetInfo fails, %s \n", cudaGetErrorString(cuda_status) );
    exit(EXIT_FAILURE);
  }

  *free_db = (double)free_byte ;
  *total_db = (double)total_byte ;
  *used_db = *total_db - *free_db ;
  return;
}

/* ----------------------------------------------------------------------------------------------- */

// Saves GPU memory usage to file
void output_free_memory(int myrank,char* info_str) {

  TRACE("output_free_memory");

  FILE* fp;
  char filename[BUFSIZ];
  double free_db,used_db,total_db;
  int do_output_info;

  // by default, only main process outputs device infos to avoid file cluttering
  do_output_info = 0;
  if (myrank == 0) {
    do_output_info = 1;
    sprintf(filename,OUTPUT_FILES_PATH"/gpu_device_mem_usage.txt");
  }
  // debugging
  if (DEBUG) {
    do_output_info = 1;
    sprintf(filename,OUTPUT_FILES_PATH"/gpu_device_mem_usage_proc_%06d.txt",myrank);
  }

  // outputs to file
  if (do_output_info) {

    // gets memory usage
    get_free_memory(&free_db,&used_db,&total_db);

    // file output
    fp = fopen(filename,"a+");
    if (fp != NULL) {
      fprintf(fp,"%d: @%s GPU memory usage: used = %f MB, free = %f MB, total = %f MB\n", myrank, info_str,
              used_db/1024.0/1024.0, free_db/1024.0/1024.0, total_db/1024.0/1024.0);
      fclose(fp);
    }
  }
}

/* ----------------------------------------------------------------------------------------------- */

// Fortran-callable version of above method
extern "C"
void FC_FUNC_(output_free_device_memory,
              OUTPUT_FREE_DEVICE_MEMORY)(int* myrank_f) {
  TRACE("output_free_device_memory");

  char info[64];
  int myrank = *myrank_f;

  sprintf(info,"f %d:",myrank);
  output_free_memory(myrank,info);
}


/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(get_free_device_memory,
              get_FREE_DEVICE_MEMORY)(realw* free, realw* used, realw* total) {
  TRACE("get_free_device_memory");

  double free_db,used_db,total_db;

  get_free_memory(&free_db,&used_db,&total_db);

  // converts to MB
  *free = (realw) free_db/1024.0/1024.0;
  *used = (realw) used_db/1024.0/1024.0;
  *total = (realw) total_db/1024.0/1024.0;
  return;
}



/* ----------------------------------------------------------------------------------------------- */

// Auxiliary functions

/* ----------------------------------------------------------------------------------------------- */

/*
__global__ void memset_to_realw_kernel(realw* array, int size, realw value){

  unsigned int tid = threadIdx.x;
  unsigned int bx = blockIdx.y*gridDim.x+blockIdx.x;
  unsigned int i = tid + bx*blockDim.x;

  if (i < size) {
    array[i] = *value;
  }
}
*/

/* ----------------------------------------------------------------------------------------------- */

realw get_device_array_maximum_value(realw* array, int size){

// get maximum of array on GPU by copying over to CPU and handle it there

  realw max = 0.0f;

  // checks if anything to do
  if (size > 0) {
    realw* h_array;

    // explicitly wait for cuda kernels to finish
    // (cudaMemcpy implicitly synchronizes all other cuda operations)
    synchronize_cuda();

    h_array = (realw*)calloc(size,sizeof(realw));
    print_CUDA_error_if_any(cudaMemcpy(h_array,array,sizeof(realw)*size,cudaMemcpyDeviceToHost),33001);

    // finds maximum value in array
    max = h_array[0];
    for( int i=1; i < size; i++){
      if (abs(h_array[i]) > max ) max = abs(h_array[i]);
    }
    free(h_array);
  }
  return max;
}



/* ----------------------------------------------------------------------------------------------- */

// ACOUSTIC simulations

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(get_norm_acoustic_from_device,
              GET_NORM_ACOUSTIC_FROM_DEVICE)(realw* norm,long* Mesh_pointer,const int* FORWARD_OR_ADJOINT) {

  TRACE("get_norm_acoustic_from_device");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  realw max1 = 0.0;
  //realw max2 = 0.0;
  //realw max3 = 0.0;
  realw *d_max1;
  //realw *d_max2;
  //realw *d_max3;

  //initializes
  *norm = 0.0f;

  // way 2 b: timing Elapsed time: 1.236916e-03
  // launch simple reduction kernel
  int blocksize = BLOCKSIZE_TRANSFER;

  int size = mp->NGLOB_AB;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // on host (allocates & initializes to zero)
  realw* h_max1 = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));
  //realw* h_max2 = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));
  //realw* h_max3 = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));

  // allocates memory on device
  print_CUDA_error_if_any(cudaMalloc((void**)&d_max1,num_blocks_x*num_blocks_y*sizeof(realw)),78001);
  // initializes values to zero
  print_CUDA_error_if_any(cudaMemset(d_max1,0,num_blocks_x*num_blocks_y*sizeof(realw)),77002);

  //print_CUDA_error_if_any(cudaMalloc((void**)&d_max2,num_blocks_x*num_blocks_y*sizeof(realw)),78001);
  //print_CUDA_error_if_any(cudaMemset(d_max2,0,num_blocks_x*num_blocks_y*sizeof(realw)),77002);

  //print_CUDA_error_if_any(cudaMalloc((void**)&d_max3,num_blocks_x*num_blocks_y*sizeof(realw)),78001);
  //print_CUDA_error_if_any(cudaMemset(d_max3,0,num_blocks_x*num_blocks_y*sizeof(realw)),77002);

  if (*FORWARD_OR_ADJOINT == 1) {
    get_maximum_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_acoustic,size,d_max1);
    //get_maximum_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_acoustic,size,d_max2);
    //get_maximum_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_potential_dot_dot_acoustic,size,d_max3);
  }else if (*FORWARD_OR_ADJOINT == 3) {
    get_maximum_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_acoustic,size,d_max1);
    //get_maximum_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_dot_acoustic,size,d_max2);
    //get_maximum_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_potential_dot_dot_acoustic,size,d_ma3);
  }

  GPU_ERROR_CHECKING ("get_maximum_kernel");

  // synchronizes
  //synchronize_cuda();
  // explicitly waits for stream to finish
  // (cudaMemcpy implicitly synchronizes all other cuda operations)
  cudaStreamSynchronize(mp->compute_stream);

  print_CUDA_error_if_any(cudaMemcpy(h_max1,d_max1,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),222);

  //print_CUDA_error_if_any(cudaMemcpy(h_max2,d_max2,num_blocks_x*num_blocks_y*sizeof(realw),
  //                                   cudaMemcpyDeviceToHost),222);

  //print_CUDA_error_if_any(cudaMemcpy(h_max3,d_max3,num_blocks_x*num_blocks_y*sizeof(realw),
  //                                   cudaMemcpyDeviceToHost),222);

  // determines max for all blocks
  max1 = h_max1[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if (max1 < h_max1[i]) max1 = h_max1[i];
  }
  /*
  max2 = h_max2[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if (max2 < h_max2[i]) max2 = h_max2[i];
  }
  */
  /*
  max3 = h_max3[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if (max3 < h_max3[i]) max3 = h_max3[i];
  }
  */
  cudaFree(d_max1);
  free(h_max1);
  //cudaFree(d_max2);
  //free(h_max2);
  //cudaFree(d_max3);
  //free(h_max3);

  // return result
  //norm[0] = max1;
  //norm[1] = max2;
  //norm[2] = max3;
  // only potential_acoustic for now:
  *norm = max1;

  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);

  GPU_ERROR_CHECKING ("get_norm_acoustic_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

// ELASTIC simulations

/* ----------------------------------------------------------------------------------------------- */


extern "C"
void FC_FUNC_(get_norm_elastic_from_device,
              GET_NORM_ELASTIC_FROM_DEVICE)(realw* norm,long* Mesh_pointer,const int* FORWARD_OR_ADJOINT) {

  TRACE("\tget_norm_elastic_from_device");
  //double start_time = get_time();

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container
  realw max1; //,max2,max3;
  realw res1; //,res2,res3;
  realw *d_max1;
  //realw *d_max2;
  //realw *d_max3;

  //initializes
  *norm = 0.0f;

  // launch simple reduction kernel
  int blocksize = BLOCKSIZE_TRANSFER;

  int size = mp->NGLOB_AB;
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

  // on host (allocates & initializes to zero)
  realw* h_max1 = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));
  //realw* h_max2 = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));
  //realw* h_max3 = (realw*) calloc(num_blocks_x*num_blocks_y,sizeof(realw));

  // allocates memory on device
  // and initializes values to zero
  print_CUDA_error_if_any(cudaMalloc((void**)&d_max1,num_blocks_x*num_blocks_y*sizeof(realw)),77001);
  print_CUDA_error_if_any(cudaMemset(d_max1,0,num_blocks_x*num_blocks_y*sizeof(realw)),77002);

  //print_CUDA_error_if_any(cudaMalloc((void**)&d_max2,num_blocks_x*num_blocks_y*sizeof(realw)),77001);
  //print_CUDA_error_if_any(cudaMemset(d_max2,0,num_blocks_x*num_blocks_y*sizeof(realw)),77002);

  //print_CUDA_error_if_any(cudaMalloc((void**)&d_max3,num_blocks_x*num_blocks_y*sizeof(realw)),77001);
  //print_CUDA_error_if_any(cudaMemset(d_max3,0,num_blocks_x*num_blocks_y*sizeof(realw)),77002);

  if (*FORWARD_OR_ADJOINT == 1) {
    get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_displ,size,d_max1);
    //get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_veloc,size,d_max2);
    //get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_accel,size,d_max3);
  }else if (*FORWARD_OR_ADJOINT == 3) {
    get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_displ,size,d_max1);
    //get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_veloc,size,d_max2);
    //get_maximum_vector_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->d_b_accel,size,d_max3);
  }

  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);

  GPU_ERROR_CHECKING ("get_norm_elastic_from_device");

  // synchronizes
  //synchronize_cuda();
  // explicitly waits for stream to finish
  // (cudaMemcpy implicitly synchronizes all other cuda operations)
  cudaStreamSynchronize(mp->compute_stream);

  // copies reduction array back to CPU
  print_CUDA_error_if_any(cudaMemcpy(h_max1,d_max1,num_blocks_x*num_blocks_y*sizeof(realw),
                                     cudaMemcpyDeviceToHost),222);
  //print_CUDA_error_if_any(cudaMemcpy(h_max2,d_max2,num_blocks_x*num_blocks_y*sizeof(realw),
  //                                   cudaMemcpyDeviceToHost),222);
  //print_CUDA_error_if_any(cudaMemcpy(h_max3,d_max3,num_blocks_x*num_blocks_y*sizeof(realw),
  //                                   cudaMemcpyDeviceToHost),222);

  // determines max for all blocks
  max1 = h_max1[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if (max1 < h_max1[i]) max1 = h_max1[i];
  }
  /*
  max2 = h_max2[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if (max2 < h_max2[i]) max2 = h_max2[i];
  }
  */
  /*
  max3 = h_max3[0];
  for(int i=1;i<num_blocks_x*num_blocks_y;i++) {
    if (max3 < h_max3[i]) max3 = h_max3[i];
  }
  */
  res1 = sqrt(max1);
  //res2 = sqrt(max2);
  //res3 = sqrt(max3);
  // debug
  // int it = 0; // not available.
  //printf("iteration %d  Valeur max de displ : %.12f  processus %d \n",it,res1,mp->myrank);
  //printf("iteration %d  Valeur max de veloc : %.12f  processus %d \n",it,res2,mp->myrank);
  //printf("iteration %d  Valeur max de accel : %.12f  processus %d \n",it,res3,mp->myrank);

  // return result
  *norm = res1;

  // debug
  //printf("rank % d - type: %d norm: %e \n",mp->myrank,*type,*norm);

  cudaFree(d_max1);
  free(h_max1);
  //cudaFree(d_max2);
  //free(h_max2);
  //cudaFree(d_max3);
  //free(h_max3);

  //double end_time = get_time();
  //printf("Elapsed time: %e\n",end_time-start_time);

  GPU_ERROR_CHECKING ("get_norm_elastic_from_device");
}

/* ----------------------------------------------------------------------------------------------- */

// unused ...

/* ----------------------------------------------------------------------------------------------- */

/*
extern "C"
void FC_FUNC_(get_max_accel,
              GET_MAX_ACCEL)(int* itf,int* sizef,long* Mesh_pointer) {

TRACE("get_max_accel");

  Mesh* mp = (Mesh*)(*Mesh_pointer);
  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  int size = *sizef;
  int it = *itf;
  realw* accel_cpy = (realw*)malloc(size*sizeof(realw));
  cudaMemcpy(accel_cpy,mp->d_accel,size*sizeof(realw),cudaMemcpyDeviceToHost);
  realw maxval=0;
  for(int i=0;i<size;++i) {
    maxval = MAX(maxval,accel_cpy[i]);
  }
  printf("%d/%d: max=%e\n",it,procid,maxval);
  free(accel_cpy);
}
*/


/* ----------------------------------------------------------------------------------------------- */
//daniel: helper function
/*
 __global__ void check_phase_ispec_kernel(int num_phase_ispec,
 int* phase_ispec,
 int NSPEC_AB,
 int* ier) {

 int i,ispec,iphase,count0,count1;
 *ier = 0;

 for(iphase=0; iphase < 2; iphase++){
 count0 = 0;
 count1 = 0;

 for(i=0; i < num_phase_ispec; i++){
 ispec = phase_ispec[iphase*num_phase_ispec + i] - 1;
 if (ispec < -1 || ispec >= NSPEC_AB) {
 printf("Error in d_phase_ispec_inner_elastic %d %d\n",i,ispec);
 *ier = 1;
 return;
 }
 if (ispec >= 0) { count0++;}
 if (ispec < 0) { count1++;}
 }

 printf("check_phase_ispec done: phase %d, count = %d %d \n",iphase,count0,count1);

 }
 }

 void check_phase_ispec(long* Mesh_pointer,int type){

 Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

 printf("check phase_ispec for type=%d\n",type);

 dim3 grid(1,1);
 dim3 threads(1,1,1);

 int* h_debug = (int*) calloc(1,sizeof(int));
 int* d_debug;
 cudaMalloc((void**)&d_debug,sizeof(int));

 if (type == 1) {
 check_phase_ispec_kernel<<<grid,threads>>>(mp->num_phase_ispec_elastic,
 mp->d_phase_ispec_inner_elastic,
 mp->NSPEC_AB,
 d_debug);
 }else if (type == 2) {
 check_phase_ispec_kernel<<<grid,threads>>>(mp->num_phase_ispec_acoustic,
 mp->d_phase_ispec_inner_acoustic,
 mp->NSPEC_AB,
 d_debug);
 }

 cudaMemcpy(h_debug,d_debug,1*sizeof(int),cudaMemcpyDeviceToHost);
 cudaFree(d_debug);
 if (*h_debug != 0) {printf("error for type=%d\n",type); exit(1);}
 free(h_debug);
 fflush(stdout);

 GPU_ERROR_CHECKING ("check_phase_ispec");
 }
*/

/* ----------------------------------------------------------------------------------------------- */
//daniel: helper function
/*
 __global__ void check_ispec_is_kernel(int NSPEC_AB,
 int* ispec_is,
 int* ier) {

 int ispec,count0,count1;

 *ier = 0;
 count0 = 0;
 count1 = 0;
 for(ispec=0; ispec < NSPEC_AB; ispec++){
 if (ispec_is[ispec] < -1 || ispec_is[ispec] > 1) {
 printf("Error in ispec_is %d %d\n",ispec,ispec_is[ispec]);
 *ier = 1;
 return;
 //exit(1);
 }
 if (ispec_is[ispec] == 0) {count0++;}
 if (ispec_is[ispec] != 0) {count1++;}
 }
 printf("check_ispec_is done: count = %d %d\n",count0,count1);
 }

 void check_ispec_is(long* Mesh_pointer,int type){

 Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

 printf("check ispec_is for type=%d\n",type);

 dim3 grid(1,1);
 dim3 threads(1,1,1);

 int* h_debug = (int*) calloc(1,sizeof(int));
 int* d_debug;
 cudaMalloc((void**)&d_debug,sizeof(int));

 if (type == 0) {
 check_ispec_is_kernel<<<grid,threads>>>(mp->NSPEC_AB,
 mp->d_ispec_is_inner,
 d_debug);
 }else if (type == 1) {
 check_ispec_is_kernel<<<grid,threads>>>(mp->NSPEC_AB,
 mp->d_ispec_is_elastic,
 d_debug);
 }else if (type == 2) {
 check_ispec_is_kernel<<<grid,threads>>>(mp->NSPEC_AB,
 mp->d_ispec_is_acoustic,
 d_debug);
 }

 cudaMemcpy(h_debug,d_debug,1*sizeof(int),cudaMemcpyDeviceToHost);
 cudaFree(d_debug);
 if (*h_debug != 0) {printf("error for type=%d\n",type); exit(1);}
 free(h_debug);
 fflush(stdout);

 GPU_ERROR_CHECKING ("check_ispec_is");
 }
*/
/* ----------------------------------------------------------------------------------------------- */
//daniel: helper function
/*
 __global__ void check_array_ispec_kernel(int num_array_ispec,
 int* array_ispec,
 int NSPEC_AB,
 int* ier) {

 int i,ispec,count0,count1;

 *ier = 0;
 count0 = 0;
 count1 = 0;

 for(i=0; i < num_array_ispec; i++){
 ispec = array_ispec[i] - 1;
 if (ispec < -1 || ispec >= NSPEC_AB) {
 printf("Error in d_array_ispec %d %d\n",i,ispec);
 *ier = 1;
 return;
 }
 if (ispec >= 0) { count0++;}
 if (ispec < 0) { count1++;}
 }

 printf("check_array_ispec done: count = %d %d \n",count0,count1);
 }

 void check_array_ispec(long* Mesh_pointer,int type){

 Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

 printf("check array_ispec for type=%d\n",type);

 dim3 grid(1,1);
 dim3 threads(1,1,1);

 int* h_debug = (int*) calloc(1,sizeof(int));
 int* d_debug;
 cudaMalloc((void**)&d_debug,sizeof(int));

 if (type == 1) {
 check_array_ispec_kernel<<<grid,threads>>>(mp->d_num_abs_boundary_faces,
 mp->d_abs_boundary_ispec,
 mp->NSPEC_AB,
 d_debug);
 }

 cudaMemcpy(h_debug,d_debug,1*sizeof(int),cudaMemcpyDeviceToHost);
 cudaFree(d_debug);
 if (*h_debug != 0) {printf("error for type=%d\n",type); exit(1);}
 free(h_debug);
 fflush(stdout);

 GPU_ERROR_CHECKING ("check_array_ispec");
 }
*/

/* ----------------------------------------------------------------------------------------------- */

// Check functions

/* ----------------------------------------------------------------------------------------------- */

//max: helper functions

/*
extern "C"
void FC_FUNC_(check_max_norm_displ_gpu,
              CHECK_MAX_NORM_DISPL_GPU)(int* size, realw* displ,long* Mesh_pointer,int* announceID) {

  TRACE("check_max_norm_displ_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  cudaMemcpy(displ, mp->d_displ,*size*sizeof(realw),cudaMemcpyDeviceToHost);
  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(displ[i]));
  }
  printf("%d: maxnorm of forward displ = %e\n",*announceID,maxnorm);
}
*/

/* ----------------------------------------------------------------------------------------------- */
/*
extern "C"
void FC_FUNC_(check_max_norm_vector,
              CHECK_MAX_NORM_VECTOR)(int* size, realw* vector1, int* announceID) {

TRACE("check_max_norm_vector");

  int procid;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD,&procid);
#else
  procid = 0;
#endif
  realw maxnorm=0;
  int maxloc;
  for(int i=0;i<*size;i++) {
    if (maxnorm<fabsf(vector1[i])) {
      maxnorm = vector1[i];
      maxloc = i;
    }
  }
  printf("%d:maxnorm of vector %d [%d] = %e\n",procid,*announceID,maxloc,maxnorm);
}
*/

/* ----------------------------------------------------------------------------------------------- */

/*
extern "C"
void FC_FUNC_(check_max_norm_displ,
              CHECK_MAX_NORM_DISPL)(int* size, realw* displ, int* announceID) {

TRACE("check_max_norm_displ");

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(displ[i]));
  }
  printf("%d: maxnorm of forward displ = %e\n",*announceID,maxnorm);
}
*/

/* ----------------------------------------------------------------------------------------------- */

/*
extern "C"
void FC_FUNC_(check_max_norm_b_displ_gpu,
              CHECK_MAX_NORM_B_DISPL_GPU)(int* size, realw* b_displ,long* Mesh_pointer,int* announceID) {

  TRACE("check_max_norm_b_displ_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  realw* b_accel = (realw*)malloc(*size*sizeof(realw));

  cudaMemcpy(b_displ, mp->d_b_displ,*size*sizeof(realw),cudaMemcpyDeviceToHost);
  cudaMemcpy(b_accel, mp->d_b_accel,*size*sizeof(realw),cudaMemcpyDeviceToHost);

  realw maxnorm=0;
  realw maxnorm_accel=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_displ[i]));
    maxnorm_accel = MAX(maxnorm,fabsf(b_accel[i]));
  }
  free(b_accel);
  printf("%d: maxnorm of backward displ = %e\n",*announceID,maxnorm);
  printf("%d: maxnorm of backward accel = %e\n",*announceID,maxnorm_accel);
}
*/

/* ----------------------------------------------------------------------------------------------- */

/*
extern "C"
void FC_FUNC_(check_max_norm_b_accel_gpu,
              CHECK_MAX_NORM_B_ACCEL_GPU)(int* size, realw* b_accel,long* Mesh_pointer,int* announceID) {

  TRACE("check_max_norm_b_accel_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  cudaMemcpy(b_accel, mp->d_b_accel,*size*sizeof(realw),cudaMemcpyDeviceToHost);

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_accel[i]));
  }
  printf("%d: maxnorm of backward accel = %e\n",*announceID,maxnorm);
}
*/

/* ----------------------------------------------------------------------------------------------- */

/*
extern "C"
void FC_FUNC_(check_max_norm_b_veloc_gpu,
              CHECK_MAX_NORM_B_VELOC_GPU)(int* size, realw* b_veloc,long* Mesh_pointer,int* announceID) {

  TRACE("check_max_norm_b_veloc_gpu");

  Mesh* mp = (Mesh*)(*Mesh_pointer); //get mesh pointer out of fortran integer container

  cudaMemcpy(b_veloc, mp->d_b_veloc,*size*sizeof(realw),cudaMemcpyDeviceToHost);

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_veloc[i]));
  }
  printf("%d: maxnorm of backward veloc = %e\n",*announceID,maxnorm);
}
*/

/* ----------------------------------------------------------------------------------------------- */

/*
extern "C"
void FC_FUNC_(check_max_norm_b_displ,
              CHECK_MAX_NORM_B_DISPL)(int* size, realw* b_displ,int* announceID) {

TRACE("check_max_norm_b_displ");

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_displ[i]));
  }
  printf("%d:maxnorm of backward displ = %e\n",*announceID,maxnorm);
}
*/

/* ----------------------------------------------------------------------------------------------- */

/*
extern "C"
void FC_FUNC_(check_max_norm_b_accel,
              CHECK_MAX_NORM_B_ACCEL)(int* size, realw* b_accel,int* announceID) {

TRACE("check_max_norm_b_accel");

  realw maxnorm=0;

  for(int i=0;i<*size;i++) {
    maxnorm = MAX(maxnorm,fabsf(b_accel[i]));
  }
  printf("%d:maxnorm of backward accel = %e\n",*announceID,maxnorm);
}
*/

/* ----------------------------------------------------------------------------------------------- */

/*
extern "C"
void FC_FUNC_(check_error_vectors,
              CHECK_ERROR_VECTORS)(int* sizef, realw* vector1,realw* vector2) {

TRACE("check_error_vectors");

  int size = *sizef;

  double diff2 = 0;
  double sum = 0;
  double temp;
  double maxerr=0;
  int maxerrorloc;

  for(int i=0;i<size;++i) {
    temp = vector1[i]-vector2[i];
    diff2 += temp*temp;
    sum += vector1[i]*vector1[i];
    if (maxerr < fabsf(temp)) {
      maxerr = abs(temp);
      maxerrorloc = i;
    }
  }

  printf("rel error = %f, maxerr = %e @ %d\n",diff2/sum,maxerr,maxerrorloc);
  int myrank;
#ifdef WITH_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
#else
  myrank = 0;
#endif
  if (myrank == 0) {
    for(int i=maxerrorloc;i>maxerrorloc-5;i--) {
      printf("[%d]: %e vs. %e\n",i,vector1[i],vector2[i]);
    }
  }

}
*/


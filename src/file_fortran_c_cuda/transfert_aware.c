#include <stdio.h>
#include <cuda.h>
#include <cublas.h>

#include <sys/time.h>
#include <sys/resource.h>

#include "config.h"
#include "mesh_constants_cuda.h"

#include <mpi.h>
#include "precisionc.h"


extern
void FC_FUNC_(transfer_accel_aware,
              TRANSFER_ACCEL_AWARE)(long * Mesh_pointer,int* NPROC, int * max_nibool_interfaces_ext_meshf,
                                    int * nibool_interfaces_ext_mesh,int * my_neighbours, MPI_Request * tab_requests_send_recv_vector,
                                    int * inum_interfaces_elastic) {



Mesh* mp = (Mesh*)(*Mesh_pointer);
int iinterface,num_int;
int max_nibool_interfaces_ext_mesh=* max_nibool_interfaces_ext_meshf;



MPI_Status status[mp->ninterface_elastic];


if(*NPROC > 1) {

    for(iinterface = 0;iinterface < mp->ninterface_elastic;iinterface++)
    {

     num_int = inum_interfaces_elastic[iinterface]-1;

MPI_Isend(mp->d_send_accel_buffer + 2*max_nibool_interfaces_ext_mesh*num_int,2*nibool_interfaces_ext_mesh[num_int], CUSTOM_MPI_TYPE,
         my_neighbours[num_int],0,MPI_COMM_WORLD,&tab_requests_send_recv_vector[num_int]);


MPI_Irecv(mp->d_recv_accel_buffer + 2*max_nibool_interfaces_ext_mesh*iinterface,2*nibool_interfaces_ext_mesh[num_int], CUSTOM_MPI_TYPE,
          my_neighbours[num_int],0,MPI_COMM_WORLD,&tab_requests_send_recv_vector[num_int+ mp->ninterface_elastic ]);

/*
MPI_Sendrecv(mp->d_send_accel_buffer + 2*max_nibool_interfaces_ext_mesh*iinterface,2*nibool_interfaces_ext_mesh[iinterface], CUSTOM_MPI_TYPE,
         my_neighbours[iinterface],0,mp->d_recv_accel_buffer + 2*max_nibool_interfaces_ext_mesh*iinterface,2*nibool_interfaces_ext_mesh[iinterface], CUSTOM_MPI_TYPE,
          my_neighbours[iinterface],0,MPI_COMM_WORLD,&status[iinterface]);*/


   }

}

}


void FC_FUNC_(receive_accel_aware,
              RECEIVE_ACCEL_AWARE)(MPI_Request * tab_requests_send_recv_vector,int* ninterface_elastic,int *NPROC,int * inum_interfaces_elastic){


MPI_Status status[ *ninterface_elastic];
int iinterface,num_int;  


if(*NPROC > 1) {

    for(iinterface = 0;iinterface < *ninterface_elastic;iinterface++)

{   num_int=inum_interfaces_elastic[iinterface]-1;
    
    MPI_Wait(&tab_requests_send_recv_vector[num_int+*ninterface_elastic],&status[num_int]);


}
}
}

void FC_FUNC_(wait_sent_request,
              WAIT_SENT_REQUEST)(MPI_Request * tab_requests_send_recv_vector,int* ninterface_elastic,int *NPROC,int * inum_interfaces_elastic){


MPI_Status status[ *ninterface_elastic];
int iinterface,num_int;  


if(*NPROC > 1) {

    for(iinterface = 0;iinterface < *ninterface_elastic;iinterface++)

{   num_int=inum_interfaces_elastic[iinterface]-1;
    
    MPI_Wait(&tab_requests_send_recv_vector[num_int],&status[num_int]);


}
}
}









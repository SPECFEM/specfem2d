/*
**************************

check_opencl_device utility

**************************

this utility program will output GPU OpenCL device informations helpful for debugging OpenCL.


for compilation, see the command-line examples given here:

example on OsX:

clang -framework OpenCL -o check_opencl_device check_opencl_device.c

or on linux:

gcc -lOpenCL -o check_opencl_device check_opencl_device.c

execute by:

./check_opencl_device

*/


// includes

#include <stdio.h>
#include <stdlib.h>
#ifdef __APPLE__
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

/*----------------------------------------------------------------------------------------*/

//helper functions

cl_int mocl_errcode;

/* The OpenCL Extension Wrangler Library
 * https://code.google.com/p/clew/
 * MIT License
 * */

const char* clewErrorString (cl_int error) {
  static const char* strings[] = {
    // Error Codes
    "CL_SUCCESS"                                  //   0
    , "CL_DEVICE_NOT_FOUND"                         //  -1
    , "CL_DEVICE_NOT_AVAILABLE"                     //  -2
    , "CL_COMPILER_NOT_AVAILABLE"                   //  -3
    , "CL_MEM_OBJECT_ALLOCATION_FAILURE"            //  -4
    , "CL_OUT_OF_RESOURCES"                         //  -5
    , "CL_OUT_OF_HOST_MEMORY"                       //  -6
    , "CL_PROFILING_INFO_NOT_AVAILABLE"             //  -7
    , "CL_MEM_COPY_OVERLAP"                         //  -8
    , "CL_IMAGE_FORMAT_MISMATCH"                    //  -9
    , "CL_IMAGE_FORMAT_NOT_SUPPORTED"               //  -10
    , "CL_BUILD_PROGRAM_FAILURE"                    //  -11
    , "CL_MAP_FAILURE"                              //  -12

    , ""    //  -13
    , ""    //  -14
    , ""    //  -15
    , ""    //  -16
    , ""    //  -17
    , ""    //  -18
    , ""    //  -19

    , ""    //  -20
    , ""    //  -21
    , ""    //  -22
    , ""    //  -23
    , ""    //  -24
    , ""    //  -25
    , ""    //  -26
    , ""    //  -27
    , ""    //  -28
    , ""    //  -29

    , "CL_INVALID_VALUE"                            //  -30
    , "CL_INVALID_DEVICE_TYPE"                      //  -31
    , "CL_INVALID_PLATFORM"                         //  -32
    , "CL_INVALID_DEVICE"                           //  -33
    , "CL_INVALID_CONTEXT"                          //  -34
    , "CL_INVALID_QUEUE_PROPERTIES"                 //  -35
    , "CL_INVALID_COMMAND_QUEUE"                    //  -36
    , "CL_INVALID_HOST_PTR"                         //  -37
    , "CL_INVALID_MEM_OBJECT"                       //  -38
    , "CL_INVALID_IMAGE_FORMAT_DESCRIPTOR"          //  -39
    , "CL_INVALID_IMAGE_SIZE"                       //  -40
    , "CL_INVALID_SAMPLER"                          //  -41
    , "CL_INVALID_BINARY"                           //  -42
    , "CL_INVALID_BUILD_OPTIONS"                    //  -43
    , "CL_INVALID_PROGRAM"                          //  -44
    , "CL_INVALID_PROGRAM_EXECUTABLE"               //  -45
    , "CL_INVALID_KERNEL_NAME"                      //  -46
    , "CL_INVALID_KERNEL_DEFINITION"                //  -47
    , "CL_INVALID_KERNEL"                           //  -48
    , "CL_INVALID_ARG_INDEX"                        //  -49
    , "CL_INVALID_ARG_VALUE"                        //  -50
    , "CL_INVALID_ARG_SIZE"                         //  -51
    , "CL_INVALID_KERNEL_ARGS"                      //  -52
    , "CL_INVALID_WORK_DIMENSION"                   //  -53
    , "CL_INVALID_WORK_GROUP_SIZE"                  //  -54
    , "CL_INVALID_WORK_ITEM_SIZE"                   //  -55
    , "CL_INVALID_GLOBAL_OFFSET"                    //  -56
    , "CL_INVALID_EVENT_WAIT_LIST"                  //  -57
    , "CL_INVALID_EVENT"                            //  -58
    , "CL_INVALID_OPERATION"                        //  -59
    , "CL_INVALID_GL_OBJECT"                        //  -60
    , "CL_INVALID_BUFFER_SIZE"                      //  -61
    , "CL_INVALID_MIP_LEVEL"                        //  -62
    , "CL_INVALID_GLOBAL_WORK_SIZE"                 //  -63
    , "CL_UNKNOWN_ERROR_CODE"
  };

  if (error >= -63 && error <= 0) {
    return strings[-error];
  } else {
    return strings[64];
  }
}


static inline cl_int _clCheck(cl_int errcode, const char *file, int line, const char *func) {
  mocl_errcode = errcode;
  if (mocl_errcode != CL_SUCCESS) {
    fprintf (stderr, "OpenCL Error %d/%s at %s:%d %s\n", mocl_errcode,
             clewErrorString(mocl_errcode),
             file, line, func);
    fflush(NULL);
    exit(1);
  }
  return errcode;
}

#define clCheck(to_check) _clCheck(to_check,__FILE__, __LINE__,  __func__)

const int BUFFER_LENGTH = 1024;

/*----------------------------------------------------------------------------------------*/

int main(int argc, char* const argv[]) {

  // platform infos
  cl_int errcode = CL_SUCCESS;
  cl_platform_id *platform_ids;
  cl_uint num_platforms;
  size_t info_length;
  int i,j;
  char info[BUFFER_LENGTH];
  // device infos
  cl_uint num_devices;
  cl_device_id *devices;
  cl_device_type device_type;
  size_t max_work_group_size;
  cl_ulong mem_size;
  cl_uint units;
  char name[BUFFER_LENGTH];
  size_t image2d_max_size[2];

  // gets platform info

  // first OpenCL call
  // only gets number of platforms
  clCheck( clGetPlatformIDs(0, NULL, &num_platforms) );

  // checks if OpenCL platforms available
  if (num_platforms == 0) {
    fprintf(stderr,"OpenCL error: No OpenCL platform available!\n");
    exit(1);
  }

  platform_ids = (cl_platform_id *) malloc(num_platforms * sizeof(cl_platform_id));

  // gets platform infos
  clCheck( clGetPlatformIDs(num_platforms, platform_ids, NULL));

  // loops over all platforms
  for (j = 0; j < num_platforms; j++) {
    fprintf (stdout, "-----------\n");
    fprintf (stdout, "Platform %i:\n",j);
    fprintf (stdout, "-----------\n");

    // checks vendor and platform names
    // gets property info length
    clCheck( clGetPlatformInfo(platform_ids[j], CL_PLATFORM_NAME, 0, NULL, &info_length));
    // checks info
    if( info_length > BUFFER_LENGTH ){ fprintf(stderr,"OpenCL error: OpenCL platform info length invalid for name!\n"); exit(1);}
    // gets info string
    clCheck( clGetPlatformInfo(platform_ids[j], CL_PLATFORM_NAME, info_length, info, NULL));
    fprintf (stdout, "  Platform Name   = %s\n", info);

    // gets property info length
    clCheck( clGetPlatformInfo(platform_ids[j], CL_PLATFORM_VENDOR, 0, NULL, &info_length));
    // checks info
    if( info_length > BUFFER_LENGTH ){ fprintf(stderr,"OpenCL error: OpenCL platform info length invalid for vendor!\n"); exit(1);}
    // gets info string
    clCheck( clGetPlatformInfo(platform_ids[j], CL_PLATFORM_VENDOR, info_length, info, NULL));
    fprintf (stdout, "  Platform Vendor = %s\n", info);

    // devices

    // only gets number of devices for this platform
    clCheck( clGetDeviceIDs(platform_ids[j], CL_DEVICE_TYPE_ALL, 0, NULL, &num_devices));

    // checks
    if (num_devices > 0) {
      // fills device infos
      devices = (cl_device_id *) malloc(num_devices * sizeof(cl_device_id));

      // gets device infos
      clCheck( clGetDeviceIDs(platform_ids[j], CL_DEVICE_TYPE_ALL, num_devices, devices, NULL));

      // loops over all devices
      for (i = 0; i < num_devices; i++) {
        fprintf (stdout, "  ----------\n");
        fprintf (stdout, "  Device %i:\n",i);
        fprintf (stdout, "  ----------\n");

        // display device properties
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_NAME, sizeof(name), name, NULL));
        fprintf (stdout, "  Device Name = %s\n", name);
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_VENDOR, sizeof(name), name, NULL));
        fprintf (stdout, "  Device Vendor = %s\n", name);
        fprintf (stdout, "  Memory:\n");
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_LOCAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL));
        fprintf (stdout, "    local_mem_size (in KB) : %f\n", mem_size / 1024.f);
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_GLOBAL_MEM_SIZE, sizeof(mem_size), &mem_size, NULL));
        fprintf (stdout, "    global_mem_size (in MB): %f\n", mem_size / (1024.f * 1024.f));
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_IMAGE2D_MAX_WIDTH, sizeof(size_t), &image2d_max_size[0], NULL));
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_IMAGE2D_MAX_HEIGHT, sizeof(size_t), &image2d_max_size[1], NULL));
        fprintf (stdout, "    image2d_max_size: %zu x %zu\n", image2d_max_size[0], image2d_max_size[1]);
        fprintf(stdout,"  blocks:\n");
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof(units), &units, NULL));
        fprintf (stdout, "    max_compute_units: %u\n", units);
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof(max_work_group_size), &max_work_group_size, NULL));
        fprintf (stdout, "    max_work_group_size: %lu\n", max_work_group_size);
        fprintf(stdout,"  features:\n");
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_VERSION, sizeof(name), name, NULL));
        fprintf (stdout, "    device version : %s\n", name);
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_TYPE, sizeof(device_type), &device_type, NULL));
        fprintf (stdout, "    device type: %d\n", (int) device_type);
        clCheck( clGetDeviceInfo(devices[i], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof(units), &units, NULL));
        fprintf (stdout, "    device max_clock_frequency: %u\n", units);
        clCheck( clGetDeviceInfo(devices[i], CL_DRIVER_VERSION, sizeof(name), name, NULL));
        fprintf (stdout, "    driver version : %s\n", name);
      }
      // frees temporary array
      free(devices);
    } else {
      fprintf(stdout,"  has no OpenCL device of type %d!\n", (int) CL_DEVICE_TYPE_ALL);
    }
  }
}

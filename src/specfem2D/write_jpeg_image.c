#include <stdio.h>
#include "libjpeg/jpeglib.h"
#include <stdlib.h>

// DK DK sample code written by Junaed Sattar, Oct 2005  http://www.cim.mcgill.ca/~junaed/libjpeg.php
// DK DK modified by Dimitri Komatitsch, CNRS Marseille, France, Oct 2011

/*
 * write_jpeg_image() writes the raw image data stored in the raw_image buffer
 * to a JPEG image with default compression and smoothing options in the file
 * specified by *filename.
 *
 * returns positive integer if successful, -1 otherwise
 *
 * parameter *filename char string specifying the file name to save to
 *
 * "raw_image" points to the raw, uncompressed image
 */
int write_jpeg_image_( unsigned char *raw_image, int *width_in, int *height_in, char *filename )
{

  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;

  /* dimensions of the image we want to write */
  int width = *width_in;
  int height = *height_in;

  int bytes_per_pixel = 3;   /* or 1 for GRAYSCALE images */
  int color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

  /* this is a pointer to one row of image data */
  JSAMPROW row_pointer[1];
  FILE *outfile = fopen( filename, "wb" );

  if ( !outfile )
  {
    printf("Error opening output jpeg file %s\n!", filename );
    return -1;
  }

  /* here we set up the standard libjpeg error handler */
  cinfo.err = jpeg_std_error( &jerr );

  /* setup compression process and destination, then write JPEG header */
  jpeg_create_compress(&cinfo);

  /* this makes the library write to outfile */
  jpeg_stdio_dest(&cinfo, outfile);

  /* Setting the parameters of the output file here */
  cinfo.image_width = width;
  cinfo.image_height = height;
  cinfo.input_components = bytes_per_pixel;
  cinfo.in_color_space = color_space;

  /* default compression parameters, we should not be worried about these */
  jpeg_set_defaults( &cinfo );

  /* Now you can set any non-default parameters you wish to.
     Here we just illustrate the use of quality (quantization table) scaling */
  int quality = 92; // DK DK can be between 0 and 100
  jpeg_set_quality(&cinfo, quality, TRUE /* limit to baseline-JPEG values */);

  /* Now do the compression */
  jpeg_start_compress( &cinfo, TRUE );
  /* write one row at a time */
  while( cinfo.next_scanline < cinfo.image_height )
  {
    row_pointer[0] = &raw_image[ cinfo.next_scanline * cinfo.image_width *  cinfo.input_components];
    jpeg_write_scanlines( &cinfo, row_pointer, 1 );
  }

  /* clean up after we are done compressing */
  jpeg_finish_compress( &cinfo );
  jpeg_destroy_compress( &cinfo );
  fclose( outfile );

  /* success code is 1 */
  return 1;
}


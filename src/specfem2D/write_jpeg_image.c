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

#include "config.h"
#include <stdio.h>
#include "libjpeg/jpeglib.h"
#include <stdlib.h>
#include <ctype.h>

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
int
FC_FUNC_(write_jpeg_image,WRITE_JPEG_IMAGE)(
  unsigned char *raw_image,
  int *width_in,
  int *height_in,
  char *filename)
{

  struct jpeg_compress_struct cinfo;
  struct jpeg_error_mgr jerr;

  /* dimensions of the image we want to write */
  int width = *width_in;
  int height = *height_in;

  int bytes_per_pixel = 3;   /* or 1 for GRAYSCALE images */
  J_COLOR_SPACE color_space = JCS_RGB; /* or JCS_GRAYSCALE for grayscale images */

  /* this is a pointer to one row of image data */
  JSAMPROW row_pointer[1];

// DK DK trim white spaces in filename coming from Fortran call
// DK DK taken from http://stackoverflow.com/questions/122616/how-do-i-trim-leading-trailing-whitespace-in-a-standard-way
  char *end;

  // Trim leading white spaces
  while(isspace(*filename)) filename++;

  // find first white space after the end of the file name
  end = filename;
  while(!isspace(*end)) end++;

  // write null terminator to keep the useful part of the file name only
  *end = 0;

  FILE *outfile = fopen( filename, "wb" );

  if (!outfile )
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


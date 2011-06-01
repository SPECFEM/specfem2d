/*
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
 */

/* Dimitri Komatitsch, University of Toulouse, May 2011: */

/* added code to force Flush-To-Zero (FTZ) on Intel processors */
/* otherwise Gradual Underflow can be extremely slow. With Intel */
/* ifort one can use compiler option -ftz, but no such option exists */
/* in gcc and gfortran therefore we call an assembler routine directly here. */
/* Very precise description available at */
/* http://software.intel.com/en-us/articles/x87-and-sse-floating-point-assists-in-ia-32-flush-to-zero-ftz-and-denormals-are-zero-daz/ */
/* and at http://www.rz.uni-karlsruhe.de/rz/docs/VTune/reference/vc148.htm */
/* from http://software.intel.com/en-us/articles/x87-and-sse-floating-point-assists-in-ia-32-flush-to-zero-ftz-and-denormals-are-zero-daz : */

/* Flush-To-Zero (FTZ) mode */

/* FTZ is a method of bypassing IEEE 754 methods of dealing with */
/* invalid floating-point numbers due to underflows. This mode is much */
/* faster. Two conditions must be met for FTZ processing to occur: */

/*   * The FTZ bit (bit 15) in the MXCSR register must be masked (value = 1). */
/*   * The underflow exception (bit 11) needs to be masked (value = 1). */

#include "config.h"

#define FTZ_BIT 15
#define UNDERFLOW_EXCEPTION_MASK 11

#ifdef HAVE_XMMINTRIN
  #define FORCE_FTZ
  #include <xmmintrin.h>
#elif HAVE_EMMINTRIN
  #include <emmintrin.h>
  #define FORCE_FTZ
#endif

void
FC_FUNC_(force_ftz,FORCE_FTZ)()
{
#ifdef FORCE_FTZ
  unsigned int x;

  /* force FTZ by setting bits 11 and 15 to one */
  x = _mm_getcsr();
  x |= (1 << FTZ_BIT);
  x |= (1 << UNDERFLOW_EXCEPTION_MASK);
  _mm_setcsr(x);
#endif
}

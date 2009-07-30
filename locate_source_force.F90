
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.2
!                   ------------------------------
!
! Copyright Universite de Pau et des Pays de l'Adour, CNRS and INRIA, France.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic wave equation
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

  subroutine locate_source_force(coord,ibool,npoin,nspec,x_source,z_source,ix_source,iz_source, &
     ispec_source,iglob_source,is_proc_source,nb_proc_source,ipass)

!
!----- calculer la position reelle de la source
!

  implicit none

  include "constants.h"
#ifdef USE_MPI
  include "mpif.h"
#endif

  integer npoin,nspec,ipass
  integer ibool(NGLLX,NGLLZ,nspec)

  double precision x_source,z_source
  double precision coord(NDIM,npoin)

  integer, intent(inout)  :: is_proc_source, nb_proc_source

  integer ip,ix,iz,numelem,ilowx,ilowz,ihighx,ihighz,ix_source,iz_source,ispec_source,iglob_source

  double precision distminmax,distmin,xp,zp,dist, dist_glob

  integer  :: ierror

  ierror = 0
  is_proc_source = 0

  distminmax = -HUGEVAL

      distmin = +HUGEVAL

      ilowx = 1
      ilowz = 1
      ihighx = NGLLX
      ihighz = NGLLZ

! look for the closest grid point
      do numelem = 1,nspec

      do ix = ilowx,ihighx
      do iz = ilowz,ihighz

! global point number
        ip = ibool(ix,iz,numelem)

! coordinates of this grid point
            xp = coord(1,ip)
            zp = coord(2,ip)

            dist = sqrt((xp-x_source)**2 + (zp-z_source)**2)

! keep the point for which distance is minimum
            if(dist < distmin) then
              distmin = dist
              iglob_source = ip
              ix_source = ix
              iz_source = iz
              ispec_source = numelem
            endif

      enddo
      enddo

      enddo

  distminmax = max(distmin,distminmax)

#ifdef USE_MPI
! global minimum distance computed over all processes
  call MPI_ALLREDUCE (distminmax, dist_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)

#else
  dist_glob = distminmax

#endif

! check if this process contains the source
  if (dist_glob == distminmax) is_proc_source = 1

#ifdef USE_MPI
! determining the number of processes that contain the source (useful when the source is located on an interface)
  call MPI_ALLREDUCE (is_proc_source, nb_proc_source, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)

#else
  nb_proc_source = is_proc_source
#endif

  if (nb_proc_source < 1) call exit_MPI('error locating force source')

  if (is_proc_source == 1 .and. ipass == 1) then
     write(IOUT,200)
     write(IOUT,"(1x,f12.3,1x,f12.3,1x,f12.3,1x,f12.3,f12.3,1x,i5.5)") x_source,z_source, &
          coord(1,iglob_source),coord(2,iglob_source),distmin,nb_proc_source
     write(IOUT,*)
     write(IOUT,*)
     write(IOUT,"('Maximum distance between asked and real =',f12.3)") distminmax
  endif

#ifdef USE_MPI
  call MPI_BARRIER(MPI_COMM_WORLD,ierror)
#endif

 200 format(//1x,48('=')/,' =  S o u r c e s  ', &
  'r e a l  p o s i t i o n s  ='/1x,48('=')// &
  '    Source    x-asked      z-asked     x-obtain     z-obtain       dist    nb_proc_source'/)

  end subroutine locate_source_force


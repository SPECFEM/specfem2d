
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) January 2005
!
!========================================================================

! modify an external grid file (list of points and coordinates) to include the
! velocity model (rho, vp, vs, in this order)

  program create_earth_model

  implicit none

  integer ipoin,npoin

  double precision rho,vp,vs

  double precision, dimension(:), allocatable :: xgrid,zgrid

  include "constants.h"

! read the grid from an existing text file
  print *
  print *,'Reading the grid from an existing text file...'
  print *

  open(unit=55,file='OUTPUT_FILES/grid_points_and_model.txt',status='old')

  read(55,*) npoin

  print *,'There are ',npoin,' grid points'

  allocate(xgrid(npoin))
  allocate(zgrid(npoin))

  do ipoin = 1,npoin
    read(55,*) xgrid(ipoin),zgrid(ipoin)
  enddo

  close(55)

! write the velocity model to the same text file
  print *
  print *,'Saving the grid and the velocity model in the same text file...'
  print *

  open(unit=55,file='OUTPUT_FILES/grid_points_and_model.txt',status='unknown')

  write(55,*) npoin

  do ipoin = 1,npoin

! user should change this to assign these values depending on the position of the grid point
    rho = 2200.d0
    vp = 3000.d0
    vs = 1732.d0

    write(55,*) xgrid(ipoin),zgrid(ipoin),rho,vp,vs

  enddo

  close(55)

  end program create_earth_model


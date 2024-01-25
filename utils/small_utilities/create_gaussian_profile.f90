
! create a Gaussian profile for bathymetry

  program create_Gaussian_profile

  implicit none

  integer, parameter :: N = 200

  double precision, parameter :: horiz_size = 20000.d0, x0 = horiz_size / 2.d0, sigma = horiz_size / 40.d0, height = 300.d0

  double precision :: deltax,x

  integer :: i

  deltax = horiz_size / dble(N-1)

  do i = 1,N
    x = (i-1) * deltax
    print *,x, -3000.d0 + height * exp(-(x-x0)**2/(2*sigma**2))
  enddo

 end program create_Gaussian_profile


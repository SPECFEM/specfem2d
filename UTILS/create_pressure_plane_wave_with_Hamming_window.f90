
  program create_Hamming_window

! create a "pseudo plane wave" source file using tapering based on a Hamming window

! Dimitri Komatitsch, CNRS Marseille, France, February 2014,
! based on discussions with Paul Cristini, also from CNRS Marseille, France.

  implicit none

! pi
  double precision, parameter :: PI = 3.141592653589793d0

  integer, parameter :: NSOURCES = 1000

! the plane wave will extend in the vertical direction from zmin to zmax
! and its amplitude in between will vary as a Hamming apodization window
  double precision, parameter :: xs_center = -0.1d0
  double precision, parameter :: zs_center = 0.07d0
  double precision, parameter :: zs_size   = 0.43d0
  double precision, parameter :: zs_min = zs_center - zs_size/2.d0
  double precision, parameter :: zs_max = zs_center + zs_size/2.d0

! angle in degrees by which we rotate the plane wave
  double precision, parameter :: angle_rotate_source = -45.d0 * (PI / 180.d0)

  double precision, parameter :: factor_max = 1.d10

  integer :: isource

  double precision :: x,hamming,xs,zs,xval,zval,xprime,zprime

  do isource = 1,NSOURCES

! Hamming apodization window
! see e.g. http://docs.scipy.org/doc/numpy/reference/generated/numpy.hamming.html
! and http://www.mathworks.fr/fr/help/signal/ref/hamming.html
    x = dble(isource - 1) / dble(NSOURCES - 1)
    hamming = 0.54d0 - 0.46d0*cos(2*PI*x)

    xs = xs_center
    zs = zs_min + x * zs_size

! subtract the position of the center of rotation
    xval = xs - xs_center
    zval = zs - zs_center

! rotate it
    xprime = xval*cos(angle_rotate_source) - zval*sin(angle_rotate_source)
    zprime = xval*sin(angle_rotate_source) + zval*cos(angle_rotate_source)

! add the position of the center of rotation back
    xprime = xprime + xs_center
    zprime = zprime + zs_center

    write(*,*) '# source ',isource
    write(*,*) 'source_surf                     = .false.'
    write(*,*) 'xs                              = ',xprime
    write(*,*) 'zs                              = ',zprime
    write(*,*) 'source_type                     = 1'
    write(*,*) 'time_function_type              = 1'
    write(*,*) 'f0                              = 1.d6'
    write(*,*) 't0                              = 0.d0'
    write(*,*) 'angleforce                      = 0.0'
    write(*,*) 'Mxx                             = 1.'
    write(*,*) 'Mzz                             = 1.'
    write(*,*) 'Mxz                             = 0.'
    write(*,*) 'factor                          = ',factor_max * hamming

  enddo

  end program create_Hamming_window


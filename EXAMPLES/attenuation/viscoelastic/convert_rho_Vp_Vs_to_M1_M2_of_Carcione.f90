
  program convert_rho_Vp_Vs_to_M1_M2

!! DK DK, June 2017

  implicit none

! unrelaxed (f = +infinity) values
! double precision, parameter :: Vp = 3297.849d0
! double precision, parameter :: Vs = 2222.536d0

! relaxed (f = 0) values
  double precision, parameter :: Vp = 3000.d0
  double precision, parameter :: Vs = 2000.d0

  double precision, parameter :: rho = 2000.d0

  double precision :: M1,M2

  M2 = Vs**2 * 2.d0 * rho

  M1 = 2 * Vp**2 * rho - M2

  print *,'M1 = ',M1
  print *,'M2 = ',M2

  end program convert_rho_Vp_Vs_to_M1_M2


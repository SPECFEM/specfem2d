
  program generate

  implicit none

! DK DK DK pour l'integrale de couplage fluide/solide
      double precision xmax,zmax
      parameter(xmax = 6400.d0)
      parameter(zmax = 4800.d0)
      integer Narch
      parameter(Narch = 6)
      double precision hauteurarch
      parameter(hauteurarch = 180.d0)
      integer ntopo
      parameter(ntopo = 800)
      integer nx
      parameter(nx = 120)
! DK DK DK pour l'integrale de couplage fluide/solide

  double precision pi
  parameter(pi = 3.14159265d0)

  integer i

  double precision :: xpoint,factorarch,topo

! print *,ntopo

  do i=1,ntopo

      xpoint = dble(i-1)*xmax/dble(ntopo - 1)

      factorarch = 2.d0*pi*dble(Narch)/xmax
      topo = zmax/2.d0 + hauteurarch*dsin(factorarch*xpoint)

      print *,xpoint,topo

  enddo

  end program generate


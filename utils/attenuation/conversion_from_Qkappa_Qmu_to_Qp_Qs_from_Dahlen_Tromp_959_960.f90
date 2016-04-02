
  program conversion

! Dimitri Komatitsch, CNRS Marseille, France, October 2015

! see formulas 9.59 and 9.60 in the book of Dahlen and Tromp, 1998
! (in that book, P is called alpha and S is called beta)

  implicit none

  double precision :: Qkappa,Qmu,Qp,Qs,inverse_of_Qp,cp,cs

!!! this is for the Carcione et al. 1988 example

! enter your Qkappa and Qmu here
  Qkappa = 40.d0
  Qmu = 20.d0

! enter the cp and cs velocities of the medium here, at the frequency at which you want this conversion to be performed
  cp = 3000.d0
  cs = 2000.d0

!!! this is for the Carcione 1993 example

! enter your Qkappa and Qmu here
! Qkappa = 20.d0
! Qmu = 10.d0

! enter the cp and cs velocities of the medium here, at the frequency at which you want this conversion to be performed
! cp = 3249.d0
! cs = 2235.d0

! Qs is the same as Qmu
  Qs = Qmu

! for Qp the formula is more complex
  inverse_of_Qp = (1.d0 - (4.d0/3.d0)*(cs**2)/(cp**2))/Qkappa + (4.d0/3.d0)*(cs**2)/(cp**2)/Qmu
  Qp = 1.d0/inverse_of_Qp

! In 2D plane strain, one spatial dimension is much greater than the others
! (see for example: http://www.engineering.ucsb.edu/~hpscicom/projects/stress/introge.pdf)
! and thus kappa = lambda + mu in 2D plane strain (instead of kappa = lambda + 2/3 mu in 3D).
! See for example equation 6 in http://cherrypit.princeton.edu/papers/paper-99.pdf.
! In 2D axisymmetric I think the 2/3 coefficient is OK, but it would be worth doublechecking.
  stop 'error: this code is currently wrong because it uses the 3D formulas instead of 2D plane strain formulas for Kappa. &
    & See comments about this in the source code of this program as well as in the users manual.'

! print the result
  print *,'Qp = ',Qp
  print *,'Qs = ',Qs

  end program conversion


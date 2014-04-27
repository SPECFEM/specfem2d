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

!----
!---- subroutine to compute poroelastic velocities cpI, cpII, & cs as a function of the dominant frequency
!----

  subroutine get_poroelastic_velocities(cpIsquare,cpIIsquare,cssquare,H_biot,C_biot,M_biot,mul_fr,phil, &
             tortl,rhol_s,rhol_f,etal_f,perm,fi,f0,Q0,w_c,ATTENUATION_PORO_FLUID_PART)

  implicit none

  include "constants.h"

  double precision :: f0,w0il
  double precision :: H_biot,C_biot,M_biot
  double precision :: cpIsquare,cpIIsquare
  double precision :: cssquare,att_I,att_II
  double precision :: etal_f,rhol_f,rhol_s,rhol_bar,perm
  double precision :: mul_fr,phil,tortl

  double precision :: a_r,a_i,b_r,b_i,cc,alpha,aa1,aa2
  double precision :: xx,yy, gXI, gYI,gXII,gYII,w_c,f_c
  double precision :: wi,fi,taus,taue,Q0,bbr,bbi

  double precision :: gA,gB,sa,sb,xxs,yys
  logical :: ATTENUATION_PORO_FLUID_PART

    rhol_bar =  (1.d0 - phil)*rhol_s + phil*rhol_f

    w_c = etal_f*phil/(tortl*rhol_f*perm)
    f_c = w_c/(2*pi)

    wi=2.d0*pi*fi

    alpha=10.d0**dlog10(wi)
    w0il =  2.d0*pi*f0
    taue = (sqrt(Q0*Q0+1) +1)/(w0il*Q0)
    taus = (sqrt(Q0*Q0+1) -1)/(w0il*Q0)

     if(ATTENUATION_PORO_FLUID_PART) then
! high frequency, with memory variables
    bbr = etal_f/perm*(1.d0+alpha*alpha*taus*taue)/(1.d0 + alpha*alpha*taus*taus)
    bbi = etal_f/perm*alpha*(taue-taus)/(1.d0 + alpha*alpha*taus*taus)
     else
! low frequency
    bbr = etal_f/perm
    bbi = 0.d0
     endif

! cs
     gA = (rhol_f*tortl*rhol_bar-phil*rhol_f**2)**2/(phil*rhol_bar)**2 - (bbr**2-bbi**2)/alpha**2*&
          (phil*rhol_f/(rhol_bar*tortl) -1.d0) - bbi/alpha*phil*rhol_f/(rhol_bar*tortl)*&
          (rhol_f*tortl*rhol_bar-phil*rhol_f**2)/(phil*rhol_bar)
     gB = -2.d0*bbr*bbi/alpha**2*(phil*rhol_f/(rhol_bar*tortl) -1.d0) + bbr/alpha*phil*rhol_f/&
          (rhol_bar*tortl)*(rhol_f*tortl*rhol_bar-phil*rhol_f**2)/(phil*rhol_bar)
!
     sa = (rhol_f*tortl*rhol_bar-phil*rhol_f**2)**2/(phil*rhol_bar)**2 + (bbr**2-bbi**2)/alpha**2
     sb = 2.d0*bbr*bbi/alpha**2
!
     xxs = sa*gA + sb*gB
     yys = gA*sb - sa*gB

     cssquare = mul_fr/(rhol_bar-phil*rhol_f/tortl) * 2.d0*(gA**2+gB**2)/(sqrt(xxs**2+yys**2)+xxs)


! cpI & cpII
      a_r = rhol_bar - phil*rhol_f/tortl - phil*rhol_bar/(tortl*rhol_f)*bbi/alpha
      a_i = phil*rhol_bar/(tortl*rhol_f)*bbr
      b_r = H_biot + M_biot*phil*rhol_bar/(tortl*rhol_f) - 2.d0*phil*C_biot/tortl - &
          phil*H_biot/(tortl*rhol_f)*bbi/alpha
      b_i = phil*H_biot/(tortl*rhol_f)*bbr
      cc = phil/(tortl*rhol_f)*(H_biot*M_biot - C_biot*C_biot)
!
      xx = b_r*b_r - b_i*b_i/(alpha*alpha) - 4.d0*a_r*cc
      yy = 2.d0*b_r*b_i/alpha - 4.d0*a_i/alpha*cc
!
      gXI = a_r*(b_r + sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx)) + &
            a_i/alpha*(b_i/alpha + sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))

      gYI = a_i/alpha*(b_r + sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx)) - &
            a_r*(b_i/alpha + sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))
      gYI = -gYI

      gXII = a_r*(b_r - sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx)) + &
            a_i/alpha*(b_i/alpha - sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))

      gYII = a_i/alpha*(b_r - sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx)) - &
            a_r*(b_i/alpha - sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))
      gYII = -gYII
!
!
!
      cpIsquare = ((b_r + sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx))**2 + &
                  (b_i/alpha + sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))**2)/&
                  (sqrt(gXI**2+gYI**2) + gXI)

      cpIIsquare = ((b_r - sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx))**2 + &
                  (b_i/alpha - sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))**2)/&
                  (sqrt(gXII**2+gYII**2) + gXII)

! attenuation factors
      att_I = -alpha*sign(1.d0,yy)*sqrt(sqrt(gXI**2+gYI**2)-gXI) / &
               sqrt((b_r + sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx))**2+&
                   (b_i/alpha + sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))**2)
      att_II = -alpha*sign(1.d0,yy)*sqrt(sqrt(gXII**2+gYII**2)-gXII) / &
               sqrt((b_r - sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)+xx))**2+&
                   (b_i/alpha - sign(1.d0,yy)*sqrt(0.5)*sqrt(sqrt(xx**2+yy**2)-xx))**2)

! inverse quality factors
        aa1 = -gYI/gXI
        aa2 = -gYII/gXII

   end subroutine get_poroelastic_velocities


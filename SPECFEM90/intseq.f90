!=====================================================================
!
!                 S p e c f e m  V e r s i o n  4 . 2
!                 -----------------------------------
!
!                         Dimitri Komatitsch
!    Department of Earth and Planetary Sciences - Harvard University
!                         Jean-Pierre Vilotte
!                 Departement de Sismologie - IPGP - Paris
!                           (c) June 1998
!
!=====================================================================

  subroutine intseq
!
!=======================================================================
!
!     "i n t s e q" : Read time iteration and non linear iteration
!                     parameters
!
!=======================================================================
!
!     output variables :
!     ----------------
!               .ncycl :  number of time steps
!               .niter :  Number of non linear or corrector iterations
!               .deltat :  Time step increment
!
!               .nftfl :  Load time function for collocated nodal forces
!               .nftfk :  Load time function for kinematic constrains
!
!=======================================================================
!
  use loadft
  use iounit
  use infos
  use timeparams

  implicit none

  character(len=80) datlin

!
!-----------------------------------------------------------------------
!
!---- read first sequence parameters for dynamic analysis
!
  read(iin ,40) datlin
  read(iin ,*) ncycl,deltat,niter

!
!---- read load time functions parameters
!
  read(iin ,40) datlin
  read(iin ,*) nltfl

  if(iecho  /=  0) then
!
!----    print requested output
!
   write(iout,100) ncycl,deltat,ncycl*deltat,niter

  endif

  return
!
!---- formats
!
  40    format(a80)
  100   format(//' I t e r a t i o n   i n f o s '/1x,29('='),//5x, &
      'Number of time iterations . . . . .(ncycl) =',i8,/5x, &
      'Time step increment . . . . . . . .(deltat) =',1pe15.6,/5x, &
      'Total simulation duration . . . . . (ttot) =',1pe15.6,/5x, &
      'Number of corrector iterations. . .(niter) =',i8)

  end subroutine intseq


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

  subroutine createnum_fast()

! same as subroutine "createnum_slow" but with a faster algorithm

  use specfem_par, only: knods,ibool,shape2D,coorg,nglob,nspec,ngnod,myrank
  implicit none

  include "constants.h"

! additional arrays needed for this fast version
  integer, dimension(:), allocatable :: locval,ind,ninseg,iglob,iwork
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,work

  integer :: ispec,nseg,ioff,iseg,ig,i,j
  integer :: nxyz,ntot,ieoff,ilocnum,iy,ix,in,nnum

  double precision :: xmaxval,xminval,ymaxval,yminval,xtol,xtypdist
  double precision :: xcor,ycor


!----  create global mesh numbering
  if(myrank == 0) then
    write(IOUT,*)
    write(IOUT,*)
    write(IOUT,*) 'Generating global mesh numbering (fast version)...'
    write(IOUT,*)
  endif

  nxyz = NGLLX*NGLLZ
  ntot = nxyz*nspec

  allocate(locval(ntot))
  allocate(ind(ntot))
  allocate(ninseg(ntot))
  allocate(iglob(ntot))
  allocate(ifseg(ntot))
  allocate(xp(ntot))
  allocate(yp(ntot))
  allocate(work(ntot))
  allocate(iwork(ntot))

! compute coordinates of the grid points
  do ispec = 1,nspec
   ieoff = nxyz*(ispec - 1)
   ilocnum = 0

  do iy = 1,NGLLX
  do ix = 1,NGLLX

    ilocnum = ilocnum + 1

    xcor = zero
    ycor = zero
    do in = 1,ngnod
        nnum = knods(in,ispec)
        xcor = xcor + shape2D(in,ix,iy)*coorg(1,nnum)
        ycor = ycor + shape2D(in,ix,iy)*coorg(2,nnum)
    enddo

    xp(ilocnum + ieoff) = xcor
    yp(ilocnum + ieoff) = ycor

  enddo
  enddo

  enddo

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Establish initial pointers
  do ispec = 1,nspec
   ieoff = nxyz*(ispec -1)
   do ix = 1,nxyz
      locval (ix+ieoff) = ix+ieoff
   enddo
  enddo

! set up a local geometric tolerance

  xtypdist = +HUGEVAL

  do ispec = 1,nspec

  xminval = +HUGEVAL
  yminval = +HUGEVAL
  xmaxval = -HUGEVAL
  ymaxval = -HUGEVAL
  ieoff = nxyz*(ispec-1)
  do ilocnum = 1,nxyz
    xmaxval = max(xp(ieoff+ilocnum),xmaxval)
    xminval = min(xp(ieoff+ilocnum),xminval)
    ymaxval = max(yp(ieoff+ilocnum),ymaxval)
    yminval = min(yp(ieoff+ilocnum),yminval)
  enddo

! compute the minimum typical "size" of an element in the mesh
  xtypdist = min(xtypdist,xmaxval-xminval)
  xtypdist = min(xtypdist,ymaxval-yminval)

  enddo

! define a tolerance, small with respect to the minimum size
  xtol = SMALLVALTOL * xtypdist

  ifseg(:) = .false.
  nseg = 1
  ifseg(1) = .true.
  ninseg(1) = ntot

  do j=1,NDIM
!  Sort within each segment
   ioff=1
   do iseg=1,nseg
      if(j == 1) then
        call rank (xp(ioff),ind,ninseg(iseg))
      else
        call rank (yp(ioff),ind,ninseg(iseg))
      endif
      call swap(xp(ioff),work,ind,ninseg(iseg))
      call swap(yp(ioff),work,ind,ninseg(iseg))
      call iswap(locval(ioff),iwork,ind,ninseg(iseg))
      ioff=ioff+ninseg(iseg)
   enddo
!  Check for jumps in current coordinate
   if (j == 1) then
     do i=2,ntot
     if (abs(xp(i)-xp(i-1)) > xtol) ifseg(i)=.true.
     enddo
   else
     do i=2,ntot
     if (abs(yp(i)-yp(i-1)) > xtol) ifseg(i)=.true.
     enddo
   endif
!  Count up number of different segments
   nseg = 0
   do i=1,ntot
      if (ifseg(i)) then
         nseg = nseg+1
         ninseg(nseg) = 1
      else
         ninseg(nseg) = ninseg(nseg) + 1
      endif
   enddo
  enddo
!
!  Assign global node numbers (now sorted lexicographically!)
!
  ig = 0
  do i=1,ntot
   if (ifseg(i)) ig=ig+1
   iglob(locval(i)) = ig
  enddo

  nglob = ig

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! get result in my format
  do ispec=1,nspec
   ieoff = nxyz*(ispec - 1)
   ilocnum = 0
  do iy = 1,NGLLX
  do ix = 1,NGLLX
      ilocnum = ilocnum + 1
      ibool(ix,iy,ispec) = iglob(ilocnum + ieoff)
  enddo
  enddo
  enddo

  deallocate(locval)
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iglob)
  deallocate(ifseg)
  deallocate(xp)
  deallocate(yp)
  deallocate(work)
  deallocate(iwork)

! check the numbering obtained
  if(minval(ibool) /= 1 .or. maxval(ibool) /= nglob) call exit_MPI('Error while generating global numbering')

  end subroutine createnum_fast


!-----------------------------------------------------------------------

  subroutine rank(A,IND,N)
!
! Use Heap Sort (p 233 Numerical Recipes)
!
  implicit none

  integer N
  double precision A(N)
  integer IND(N)

  integer i,j,l,ir,indx
  double precision q

  do J=1,N
   IND(j)=j
  enddo

  if(n == 1) return
  L=n/2+1
  ir=n
  100 continue
   IF(l > 1) THEN
     l=l-1
     indx=ind(l)
     q=a(indx)
   ELSE
     indx=ind(ir)
     q=a(indx)
     ind(ir)=ind(1)
     ir=ir-1
     if(ir == 1) then
       ind(1)=indx
       return
     endif
   endif
   i=l
   j=l+l
  200 continue
   IF(J <= IR) THEN
      IF(J < IR) THEN
         IF(A(IND(j)) < A(IND(j+1))) j=j+1
      endif
      IF(q < A(IND(j))) THEN
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      endif
   GOTO 200
   endif
   IND(I)=INDX
  GOTO 100

  end subroutine rank

!-----------------------------------------------------------------------

  subroutine swap(a,w,ind,n)
!
! Use IND to sort array A (p 233 Numerical Recipes)
!
  implicit none

  integer n
  double precision A(N),W(N)
  integer IND(N)

  integer j

  W(:) = A(:)

  do J=1,N
    A(j) = W(ind(j))
  enddo

  end subroutine swap

!-----------------------------------------------------------------------

  subroutine iswap(a,w,ind,n)
!
! Use IND to sort array A
!
  implicit none

  integer n
  integer A(N),W(N),IND(N)

  integer j

  W(:) = A(:)

  do J=1,N
    A(j) = W(ind(j))
  enddo

  end subroutine iswap


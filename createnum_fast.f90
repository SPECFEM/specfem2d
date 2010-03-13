
!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
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

  subroutine createnum_fast(knods,ibool,shape,coorg,npoin,npgeo,nspec,ngnod,myrank,ipass)

! equivalent de la routine "createnum_slow" mais algorithme plus rapide

  implicit none

  include "constants.h"

  integer npoin,npgeo,nspec,ngnod,myrank,ipass
  integer knods(ngnod,nspec),ibool(NGLLX,NGLLZ,nspec)
  double precision shape(ngnod,NGLLX,NGLLX)
  double precision coorg(NDIM,npgeo)

  integer i,j

! tableaux supplementaires pour cette version rapide
  integer, dimension(:), allocatable :: loc,ind,ninseg,iglob,iwork
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,work

  integer ie,nseg,ioff,iseg,ig
  integer nxyz,ntot,ispec,ieoff,ilocnum,iy,ix,in,nnum

  double precision xmaxval,xminval,ymaxval,yminval,xtol,xtypdist
  double precision xcor,ycor


!----  create global mesh numbering
  if(myrank == 0 .and. ipass == 1) then
    write(IOUT,*)
    write(IOUT,*)
    write(IOUT,*) 'Generating global mesh numbering (fast version)...'
    write(IOUT,*)
  endif

  nxyz = NGLLX*NGLLZ
  ntot = nxyz*nspec

  allocate(loc(ntot))
  allocate(ind(ntot))
  allocate(ninseg(ntot))
  allocate(iglob(ntot))
  allocate(ifseg(ntot))
  allocate(xp(ntot))
  allocate(yp(ntot))
  allocate(work(ntot))
  allocate(iwork(ntot))

! compute coordinates of the grid points
  do ispec=1,nspec
   ieoff = nxyz*(ispec - 1)
   ilocnum = 0

  do iy = 1,NGLLX
  do ix = 1,NGLLX

    ilocnum = ilocnum + 1

    xcor = zero
    ycor = zero
    do in = 1,ngnod
        nnum = knods(in,ispec)
        xcor = xcor + shape(in,ix,iy)*coorg(1,nnum)
        ycor = ycor + shape(in,ix,iy)*coorg(2,nnum)
    enddo

    xp(ilocnum + ieoff) = xcor
    yp(ilocnum + ieoff) = ycor

  enddo
  enddo

  enddo

! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

! Establish initial pointers
  do ie=1,nspec
   ieoff = nxyz*(ie -1)
   do ix=1,nxyz
      loc (ix+ieoff) = ix+ieoff
   enddo
  enddo

! set up local geometric tolerances

  xtypdist=+HUGEVAL

  do ie=1,nspec

  xminval=+HUGEVAL
  yminval=+HUGEVAL
  xmaxval=-HUGEVAL
  ymaxval=-HUGEVAL
  ieoff=nxyz*(ie-1)
  do ilocnum=1,nxyz
    xmaxval=max(xp(ieoff+ilocnum),xmaxval)
    xminval=min(xp(ieoff+ilocnum),xminval)
    ymaxval=max(yp(ieoff+ilocnum),ymaxval)
    yminval=min(yp(ieoff+ilocnum),yminval)
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
      call iswap(loc(ioff),iwork,ind,ninseg(iseg))
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
   iglob(loc(i)) = ig
  enddo

  npoin = ig

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

  deallocate(loc)
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iglob)
  deallocate(ifseg)
  deallocate(xp)
  deallocate(yp)
  deallocate(work)
  deallocate(iwork)

! check the numbering obtained
  if(minval(ibool) /= 1 .or. maxval(ibool) /= npoin) call exit_MPI('Error while generating global numbering')

  if(myrank == 0 .and. ipass == 1) then
    write(IOUT,*)
    write(IOUT,*) 'Total number of points of the global mesh: ',npoin
    write(IOUT,*)
  endif

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
   ENDIF
   i=l
   j=l+l
  200 continue
   IF(J <= IR) THEN
      IF(J < IR) THEN
         IF(A(IND(j)) < A(IND(j+1))) j=j+1
      ENDIF
      IF(q < A(IND(j))) THEN
         IND(I)=IND(J)
         I=J
         J=J+J
      ELSE
         J=IR+1
      ENDIF
   GOTO 200
   ENDIF
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


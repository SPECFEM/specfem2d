
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

  subroutine createnum_fast(knods,ibool,shape,coorg,npoin,ndime,npgeo)
!
!=======================================================================
!
!  "c r e a t e n u m _ f a s t": Equivalent de la routine "createnum_slow"
!           mais avec un algorithme "sale mais tres rapide"
!
!  Cette version rapide necessite l'allocation de tableaux supplementaires
!
!  Cette version rapide n'accepte pas les conditions periodiques
!  En cas de conditions periodiques, utiliser la version lente
!
!=======================================================================
!

  use iounit
  use infos
  use spela202

  implicit none

  integer npoin,ndime,npgeo
  integer knods(ngnod,nspec),ibool(nxgll,nygll,nspec)
  double precision shape(ngnod,nxgll,nxgll)
  double precision coorg(ndime,npgeo)

  integer i,j

  double precision, parameter :: smallvaltol = 0.000001d0
  double precision, parameter :: HUGEVAL=1.0d+30

  double precision, parameter :: zero = 0.d0

! tableaux supplementaires pour cette version rapide

  integer, dimension(:), allocatable :: loc,ind,ninseg,iglob,iwork
  logical, dimension(:), allocatable :: ifseg
  double precision, dimension(:), allocatable :: xp,yp,work

  integer ie,nseg,ioff,iseg,ig
  integer nxyz,ntot,ispec,ieoff,ilocnum,iy,ix,in,nnum

  double precision xmaxval,xminval,ymaxval,yminval,xtol,xtypdist
  double precision xcor,ycor

!
!-----------------------------------------------------------------------
!

!
!----  create global numbering from mesh structure
!
  print *
  print *
  print *,'Generating global numbering from mesh structure (fast version)...'
  print *

  nxyz   = nxgll*nygll
  ntot   = nxyz*nspec

  print *,'Allocating a few more arrays for the fast version'

  allocate(loc(ntot))
  allocate(ind(ntot))
  allocate(ninseg(ntot))
  allocate(iglob(ntot))
  allocate(ifseg(ntot))
  allocate(xp(ntot))
  allocate(yp(ntot))
  allocate(work(ntot))
  allocate(iwork(ntot))

  print *,'Generating the numbering'

! compute coordinates of the grid points
  do ispec=1,nspec
   ieoff = nxyz*(ispec - 1)
   ilocnum = 0

  do iy = 1,nxgll
  do ix = 1,nxgll

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
    xmaxval=dmax1(xp(ieoff+ilocnum),xmaxval)
    xminval=dmin1(xp(ieoff+ilocnum),xminval)
    ymaxval=dmax1(yp(ieoff+ilocnum),ymaxval)
    yminval=dmin1(yp(ieoff+ilocnum),yminval)
  enddo

! compute the minimum typical "size" of an element in the mesh
  xtypdist = dmin1(xtypdist,xmaxval-xminval)
  xtypdist = dmin1(xtypdist,ymaxval-yminval)

  enddo

! define a tolerance, small with respect to the minimum size
  xtol=smallvaltol*xtypdist

  ifseg(:) = .false.
  nseg = 1
  ifseg(1) = .true.
  ninseg(1) = ntot

  do j=1,ndime
!  Sort within each segment
   ioff=1
   do iseg=1,nseg
      if  (j == 1) then
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

! recuperer resultat a mon format
  do ispec=1,nspec
   ieoff = nxyz*(ispec - 1)
   ilocnum = 0
  do iy = 1,nxgll
  do ix = 1,nxgll
      ilocnum = ilocnum + 1
      ibool(ix,iy,ispec) = iglob(ilocnum + ieoff)
  enddo
  enddo
  enddo

  print *,'Deallocating the arrays'

  deallocate(loc)
  deallocate(ind)
  deallocate(ninseg)
  deallocate(iglob)
  deallocate(ifseg)
  deallocate(xp)
  deallocate(yp)
  deallocate(work)
  deallocate(iwork)

! verification de la coherence de la numerotation generee
  if(minval(ibool) /= 1 .or. maxval(ibool) /= npoin) &
          stop 'Error while generating global numbering'

  print *
  print *,'Total number of points of the global mesh: ',npoin
  print *

  return
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

  if (n == 1) return
  L=n/2+1
  ir=n
  100 CONTINUE
   IF (l > 1) THEN
      l=l-1
      indx=ind(l)
      q=a(indx)
   ELSE
      indx=ind(ir)
      q=a(indx)
      ind(ir)=ind(1)
      ir=ir-1
      if (ir == 1) then
         ind(1)=indx
         return
      endif
   ENDIF
   i=l
   j=l+l
  200    CONTINUE
   IF (J <= IR) THEN
      IF (J < IR) THEN
         IF ( A(IND(j)) < A(IND(j+1)) ) j=j+1
      ENDIF
      IF (q < A(IND(j))) THEN
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

  do J=1,N
   W(j)=A(j)
  enddo

  do J=1,N
   A(j)=W(ind(j))
  enddo

  RETURN
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

  do J=1,N
   W(j)=A(j)
  enddo

  do J=1,N
   A(j)=W(ind(j))
  enddo

  RETURN
  end subroutine iswap


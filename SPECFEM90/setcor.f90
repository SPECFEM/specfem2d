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

  subroutine setcor(coord,npoin,ndime,knods,shape,ibool,coorg, &
                       nxgll,nygll,nspec,npgeo,ngnod,ioutputgrid)
!
!=======================================================================
!
!     "s e t c o r" : set the global nodal coordinates
!
!=======================================================================
!
  use iounit
  use infos
  use label1

  implicit none

  integer npoin,ndime,nxgll,nygll,nspec,npgeo,ngnod

  integer knods(ngnod,nspec),ibool(nxgll,nygll,nspec)
  double precision coord(ndime,npoin),coorg(ndime,npgeo)
  double precision shape(ngnod,nxgll,nygll)

  logical ioutputgrid

  integer n,i,ip1,ip2,ispel,in,nnum
  double precision xcor,zcor

  double precision, parameter :: zero = 0.d0

!
!----  initialisation des labels
!
  labelc(1) = '   x1'
  labelc(2) = '   x2'
  labelc(3) = '   x3'

!
!----  generation des coordonnees physiques des points globaux
!
  do ispel = 1,nspec
  do ip1 = 1,nxgll
  do ip2 = 1,nygll

      xcor = zero
      zcor = zero
      do in = 1,ngnod
          nnum = knods(in,ispel)
          xcor = xcor + shape(in,ip1,ip2)*coorg(1,nnum)
          zcor = zcor + shape(in,ip1,ip2)*coorg(2,nnum)
      enddo

      coord(1,ibool(ip1,ip2,ispel)) = xcor
      coord(2,ibool(ip1,ip2,ispel)) = zcor

   enddo
   enddo
   enddo

!
!----  check the input
!
  if(iecho  ==  2) then
   do n = 1,npoin
      if(mod(n,50)  ==  1) write(iout,100) (labelc(i),i=1,ndime)
      write(iout,200) n, (coord(i,n), i=1,ndime)
   enddo
  endif

!
!----  sauvegarde de la grille de points dans un fichier
!
  if(ioutputgrid) then
    print *
    print *,'Saving the grid in a text file...'
    print *
    open(unit=55,file='gridpoints.txt',status='unknown')
    write(55,*) npoin
    do n = 1,npoin
      write(55,*) n, (coord(i,n), i=1,ndime)
    enddo
    close(55)
  endif

  return
!
!---- formats
!
  100   format(///' n o d a l   c o o r d i n a t e    d a t a'/1x, &
         42('=')///,4x,' node number ',10x,2(a5,12x))
  200   format(4x,i7,10x,3(1pe15.8,2x))

  end subroutine setcor

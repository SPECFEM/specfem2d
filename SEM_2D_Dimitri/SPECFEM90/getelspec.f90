
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

  subroutine getelspec(knods,kmato,numabs,codeabs,codeperio,anyabs,anyperio)
!
!=======================================================================
!
!    "g e t e l s p e c": Read elements topology and material set for
!           spectral elements bloc
!
!=======================================================================
!

  use iounit
  use infos
  use spela202
  use codebord

  implicit none

  character(len=80) datlin

  integer knods(ngnod,nspec),kmato(nspec)
  integer numabs(nelemabs),codeabs(4,nelemabs)
  integer codeperio(4,nelemperio)
  logical anyabs,anyperio

  integer ie,n,i,k
  integer inum,itourne,ntourne,nperio,idummy,numabsread

  integer codeabsread(4),codeperioread(4)

!
!-----------------------------------------------------------------------
!

!
!----  read spectral macroblocs data
!
  n = 0
  read(iin,40) datlin
  do ie = 1,nspec
  read(iin,*) n,kmato(n),(knods(k,n), k=1,ngnod)
  enddo

!
!----  check the input
!
  if(iecho  ==  2) then
 do ie = 1,nspec
    if(mod(ie,50)  ==  1) write(iout,150) (i, i=1,ngnod)
    write(iout,200) ie,kmato(ie),(knods(k,ie), k=1,ngnod)
 enddo
  endif

!
!----  lire bords absorbants et bords periodiques
!
  if(anyperio) then
  read(iin ,40) datlin
  do n=1,nelemperio
    read(iin ,*) inum,codeperioread(1), &
           codeperioread(2),codeperioread(3),codeperioread(4)
    if(inum < 1 .or. inum > nelemperio) stop 'Wrong periodic element number'
    codeperio(1,inum) = codeperioread(1)
    codeperio(2,inum) = codeperioread(2)
    codeperio(3,inum) = codeperioread(3)
    codeperio(4,inum) = codeperioread(4)
  enddo
  write(*,*)
  write(*,*) 'Number of periodic elements : ',nelemperio
  endif

  if(anyabs) then
  read(iin ,40) datlin
  do n=1,nelemabs
    read(iin ,*) inum,numabsread,codeabsread(1), &
           codeabsread(2),codeabsread(3),codeabsread(4)
    if(inum < 1 .or. inum > nelemabs) stop 'Wrong absorbing element number'
    numabs(inum) = numabsread
    codeabs(ihaut,inum)   = codeabsread(1)
    codeabs(ibas,inum)    = codeabsread(2)
    codeabs(igauche,inum) = codeabsread(3)
    codeabs(idroite,inum) = codeabsread(4)

!----  eventuellement tourner element counterclockwise si condition absorbante

     if(codeabs(ibas,inum)     ==  iaretebas    .or. &
            codeabs(ihaut,inum)    ==  iaretehaut   .or. &
            codeabs(igauche,inum)  ==  iaretegauche .or. &
            codeabs(idroite,inum)  ==  iaretedroite) then
           ntourne = 0

  else if(codeabs(ibas,inum)     ==  iaretegauche .or. &
            codeabs(ihaut,inum)    ==  iaretedroite .or. &
            codeabs(igauche,inum)  ==  iaretehaut   .or. &
            codeabs(idroite,inum)  ==  iaretebas) then
           ntourne = 3

  else if(codeabs(ibas,inum)     ==  iaretehaut   .or. &
            codeabs(ihaut,inum)    ==  iaretebas    .or. &
            codeabs(igauche,inum)  ==  iaretedroite .or. &
            codeabs(idroite,inum)  ==  iaretegauche) then
           ntourne = 2

  else if(codeabs(ibas,inum)     ==  iaretedroite .or. &
            codeabs(ihaut,inum)    ==  iaretegauche .or. &
            codeabs(igauche,inum)  ==  iaretebas    .or. &
            codeabs(idroite,inum)  ==  iaretehaut) then
           ntourne = 1
  else
     stop 'Error in absorbing conditions numbering'

  endif

!----  rotate element counterclockwise
  if(ntourne  /=  0) then

  do itourne = 1,ntourne

      idummy = knods(1,numabs(inum))
      knods(1,numabs(inum)) = knods(2,numabs(inum))
      knods(2,numabs(inum)) = knods(3,numabs(inum))
      knods(3,numabs(inum)) = knods(4,numabs(inum))
      knods(4,numabs(inum)) = idummy

    if(ngnod  ==  9) then
      idummy = knods(5,numabs(inum))
      knods(5,numabs(inum)) = knods(6,numabs(inum))
      knods(6,numabs(inum)) = knods(7,numabs(inum))
      knods(7,numabs(inum)) = knods(8,numabs(inum))
      knods(8,numabs(inum)) = idummy
    endif

!----  tourner aussi le numero d'arete pour condition periodique si necessaire
  if(anyperio) then
  do nperio=1,nelemperio
    if(codeperio(1,nperio)  ==  numabs(inum)) then
        codeperio(2,nperio) = codeperio(2,nperio) - 1
        if(codeperio(2,nperio)  ==  0) codeperio(2,nperio) = 4
    endif
    if(codeperio(3,nperio)  ==  numabs(inum)) then
        codeperio(4,nperio) = codeperio(4,nperio) - 1
        if(codeperio(4,nperio)  ==  0) codeperio(4,nperio) = 4
    endif
  enddo
  endif

  enddo

  endif

  enddo
  write(*,*)
  write(*,*) 'Number of absorbing elements : ',nelemabs
  endif

  return
!
!---- formats
!
  40  format(a80)
  150 format(///' S p e c t r a l   m a c r o b l o c s   t o p o l o g y'/1x, &
            55('='),//5x,'macrobloc  material',9('  node ',i1,:,2x),/5x, &
           'number      number',//)
  200 format(4x,i7,9(3x,i7))

  end subroutine getelspec

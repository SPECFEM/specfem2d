
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.0
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) May 2004
!
!========================================================================

  subroutine createnum_slow(knods,ibool,npoin,nspec,ngnod)

! generate the global numbering

  implicit none

  include "constants.h"

  integer npoin,nspec,ngnod

  integer knods(ngnod,nspec),ibool(NGLLX,NGLLY,nspec)

  integer i,j,num2,i2,j2,ipos,ipos2,iloc,jloc,kloc
  integer ngnodloc,ngnodother,nedgeloc,nedgeother,npedge,numelem,npcorn

  logical alreadyexist

  integer ngnoddeb(4),ngnodfin(4)

!----  create global mesh numbering
  print *
  print *,'Generating global mesh numbering (slow version)...'
  print *

  npoin = 0
  npedge = 0
  npcorn = 0

! definition des aretes par rapport aux quatre points de controle

! --- arete 1 relie point 1 a point 2
  ngnoddeb(1)= 1
  ngnodfin(1)= 2

! --- arete 2 relie point 2 a point 3
  ngnoddeb(2)= 2
  ngnodfin(2)= 3

! --- arete 3 relie point 3 a point 4
  ngnoddeb(3)= 3
  ngnodfin(3)= 4

! --- arete 4 relie point 4 a point 1
  ngnoddeb(4)= 4
  ngnodfin(4)= 1

! initialisation du tableau de numerotation globale
  ibool(:,:,:) = 0

  do numelem = 1,nspec
  do i=1,NGLLX
    do j=1,NGLLY

! verifier que le point n'a pas deja ete genere

  if(ibool(i,j,numelem) == 0) then

!
!---- point interieur a un element, donc forcement unique
!
  if(i /= 1 .and. i /= NGLLX .and. j /= 1 .and. j /= NGLLY) then

    npoin = npoin + 1
    ibool(i,j,numelem) = npoin

!
!---- point au coin d'un element, rechercher les coins des autres elements
!
  else if((i == 1 .and. j == 1) .or. (i == 1 .and. j == NGLLY) .or. &
          (i == NGLLX .and. j == 1) .or. (i == NGLLX .and. j == NGLLY)) then

! trouver numero local du coin
  if(i == 1 .and. j == 1) then
    ngnodloc = 1
  else if(i == NGLLX .and. j == 1) then
    ngnodloc = 2
  else if(i == NGLLX .and. j == NGLLY) then
    ngnodloc = 3
  else if(i == 1 .and. j == NGLLY) then
    ngnodloc = 4
  endif

! rechercher si existe deja, forcement dans un element precedent

  alreadyexist = .false.

  if(numelem > 1) then

  do num2=1,numelem-1

! ne rechercher que sur les 4 premiers points de controle et non sur ngnod
    do ngnodother=1,4

! voir si ce coin a deja ete genere
      if(knods(ngnodother,num2) == knods(ngnodloc,numelem)) then
        alreadyexist = .true.

! obtenir la numerotation dans l'autre element
          if(ngnodother == 1) then
            i2 = 1
            j2 = 1
          else if(ngnodother == 2) then
            i2 = NGLLX
            j2 = 1
          else if(ngnodother == 3) then
            i2 = NGLLX
            j2 = NGLLY
          else if(ngnodother == 4) then
            i2 = 1
            j2 = NGLLY
            else
                  stop 'bad corner'
            endif

! affecter le meme numero
          ibool(i,j,numelem) = ibool(i2,j2,num2)

! sortir de la recherche
          goto 134

      endif
    enddo
  enddo

 134  continue

  endif

! si un ancien point n'a pas ete trouve, en generer un nouveau
  if(.not. alreadyexist) then
    npcorn = npcorn + 1
    npoin = npoin + 1
    ibool(i,j,numelem) = npoin
  endif

!
!---- point a l'interieur d'une arete, rechercher si autre arete correspondante
!
  else

! trouver numero local de l'arete
  if(j == 1) then
    nedgeloc = 1
  else if(i == NGLLX) then
    nedgeloc = 2
  else if(j == NGLLY) then
    nedgeloc = 3
  else if(i == 1) then
    nedgeloc = 4
  endif

! rechercher si existe deja, forcement dans un element precedent

  alreadyexist = .false.

  if(numelem > 1) then

  do num2=1,numelem-1

! rechercher sur les 4 aretes
    do nedgeother=1,4

!--- detecter un eventuel defaut dans la structure topologique du maillage

  if((knods(ngnoddeb(nedgeother),num2) == knods(ngnoddeb(nedgeloc),numelem)) &
       .and. &
    (knods(ngnodfin(nedgeother),num2) == knods(ngnodfin(nedgeloc),numelem))) then
  stop 'Improper topology of the input mesh detected'

!--- sinon voir si cette arete a deja ete generee

  else if((knods(ngnoddeb(nedgeother),num2) == knods(ngnodfin(nedgeloc),numelem)) &
       .and. &
    (knods(ngnodfin(nedgeother),num2) == knods(ngnoddeb(nedgeloc),numelem))) then

        alreadyexist = .true.

! obtenir la numerotation dans l'autre element
! maillage conforme donc on doit supposer que NGLLX == NGLLY

! generer toute l'arete pour eviter des recherches superflues
  do kloc = 2,NGLLX-1

! calculer l'abscisse le long de l'arete de depart
          if(nedgeloc == 1) then
            iloc = kloc
            jloc = 1
            ipos = iloc
          else if(nedgeloc == 2) then
            iloc = NGLLX
            jloc = kloc
            ipos = jloc
          else if(nedgeloc == 3) then
            iloc = kloc
            jloc = NGLLY
            ipos = NGLLX - iloc + 1
          else if(nedgeloc == 4) then
            iloc = 1
            jloc = kloc
            ipos = NGLLY - jloc + 1
            else
                  stop 'bad nedgeloc'
            endif

! calculer l'abscisse le long de l'arete d'arrivee
! topologie du maillage coherente, donc sens de parcours des aretes opposes

        ipos2 = NGLLX - ipos + 1

! calculer les coordonnees reelles dans l'element d'arrivee
          if(nedgeother == 1) then
            i2 = ipos2
            j2 = 1
          else if(nedgeother == 2) then
            i2 = NGLLX
            j2 = ipos2
          else if(nedgeother == 3) then
            i2 = NGLLX - ipos2 + 1
            j2 = NGLLY
          else if(nedgeother == 4) then
            i2 = 1
            j2 = NGLLY - ipos2 + 1
            else
                  stop 'bad nedgeother'
            endif

! verifier que le point de depart n'existe pas deja
      if(ibool(iloc,jloc,numelem) /= 0) stop 'point genere deux fois'

! verifier que le point d'arrivee existe bien deja
      if(ibool(i2,j2,num2) == 0) stop 'point inconnu dans le maillage'

! affecter le meme numero
      ibool(iloc,jloc,numelem) = ibool(i2,j2,num2)

  enddo

! sortir de la recherche
          goto 135

      endif
    enddo
  enddo

 135  continue

  endif

! si un ancien point n'a pas ete trouve, en generer un nouveau
  if(.not. alreadyexist) then
    npedge = npedge + 1
    npoin = npoin + 1
    ibool(i,j,numelem) = npoin
  endif

  endif

  endif

    enddo
  enddo
  enddo

! verification de la coherence de la numerotation generee
  if(minval(ibool) /= 1 .or. maxval(ibool) /= npoin) stop 'Error while generating global numbering'

  print *,'Total number of points of the global mesh: ',npoin
  print *,'distributed as follows:'
  print *
  print *,'Number of interior points: ',npoin-npedge-npcorn
  print *,'Number of edge points (without corners): ',npedge
  print *,'Number of corner points: ',npcorn
  print *

  end subroutine createnum_slow


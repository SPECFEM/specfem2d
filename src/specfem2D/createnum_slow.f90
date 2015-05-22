
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

  subroutine createnum_slow()

! generate the global numbering

  use specfem_par, only : knods,ibool,nglob,nspec,myrank

  implicit none

  include "constants.h"

  integer i,j,num2,i2,j2,ipos,ipos2,iloc,jloc,kloc
  integer ngnodloc,ngnodother,nedgeloc,nedgeother,npedge,numelem,npcorn

  logical alreadyexist

  integer, dimension(NEDGES) :: ngnod_begin,ngnod_end


!----  create global mesh numbering
  if(myrank == 0) then
    write(IOUT,*)
    write(IOUT,*) 'Generating global mesh numbering (slow version)...'
    write(IOUT,*)
  endif

  nglob = 0
  npedge = 0
  npcorn = 0

! define edges from the four control points

! --- edge 1 linking point 1 to point 2
  ngnod_begin(1)= 1
  ngnod_end(1)= 2

! --- edge 2 linking point 2 to point 3
  ngnod_begin(2)= 2
  ngnod_end(2)= 3

! --- edge 3 linking point 3 to point 4
  ngnod_begin(3)= 3
  ngnod_end(3)= 4

! --- edge 4 linking point 4 to point 1
  ngnod_begin(4)= 4
  ngnod_end(4)= 1

! initialisation du tableau de numerotation globale
  ibool(:,:,:) = 0

  do numelem = 1,nspec
  do i=1,NGLLX
    do j=1,NGLLZ

! verifier que le point n'a pas deja ete genere

  if(ibool(i,j,numelem) == 0) then

!
!---- point interieur a un element, donc forcement unique
!
  if(i /= 1 .and. i /= NGLLX .and. j /= 1 .and. j /= NGLLZ) then

    nglob = nglob + 1
    ibool(i,j,numelem) = nglob

!
!---- point au coin d'un element, rechercher les coins des autres elements
!
  else if((i == 1 .and. j == 1) .or. (i == 1 .and. j == NGLLZ) .or. &
          (i == NGLLX .and. j == 1) .or. (i == NGLLX .and. j == NGLLZ)) then

! trouver numero local du coin
  if(i == 1 .and. j == 1) then
    ngnodloc = 1
  else if(i == NGLLX .and. j == 1) then
    ngnodloc = 2
  else if(i == NGLLX .and. j == NGLLZ) then
    ngnodloc = 3
  else if(i == 1 .and. j == NGLLZ) then
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
            j2 = NGLLZ
          else if(ngnodother == 4) then
            i2 = 1
            j2 = NGLLZ
          else
             call exit_MPI('bad corner')
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
    nglob = nglob + 1
    ibool(i,j,numelem) = nglob
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
  else if(j == NGLLZ) then
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

  if((knods(ngnod_begin(nedgeother),num2) == knods(ngnod_begin(nedgeloc),numelem)) &
       .and. &
    (knods(ngnod_end(nedgeother),num2) == knods(ngnod_end(nedgeloc),numelem))) then
     call exit_MPI('Improper topology of the input mesh detected')

!--- sinon voir si cette arete a deja ete generee

  else if((knods(ngnod_begin(nedgeother),num2) == knods(ngnod_end(nedgeloc),numelem)) &
       .and. &
    (knods(ngnod_end(nedgeother),num2) == knods(ngnod_begin(nedgeloc),numelem))) then

        alreadyexist = .true.

! obtenir la numerotation dans l'autre element
! maillage conforme donc on doit supposer que NGLLX == NGLLZ

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
            jloc = NGLLZ
            ipos = NGLLX - iloc + 1
          else if(nedgeloc == 4) then
            iloc = 1
            jloc = kloc
            ipos = NGLLZ - jloc + 1
            else
               call exit_MPI('bad nedgeloc')
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
            j2 = NGLLZ
          else if(nedgeother == 4) then
            i2 = 1
            j2 = NGLLZ - ipos2 + 1
            else
               call exit_MPI('bad nedgeother')
            endif

! verifier que le point de depart n'existe pas deja
      if(ibool(iloc,jloc,numelem) /= 0) call exit_MPI('point generated twice')

! verifier que le point d'arrivee existe bien deja
      if(ibool(i2,j2,num2) == 0) call exit_MPI('unknown point in the mesh')

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
    nglob = nglob + 1
    ibool(i,j,numelem) = nglob
  endif

  endif

  endif

    enddo
  enddo
  enddo

! verification de la coherence de la numerotation generee
  if(minval(ibool) /= 1 .or. maxval(ibool) /= nglob) call exit_MPI('Error while generating global numbering')

  end subroutine createnum_slow


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

  subroutine modifperio(ibool,iboolori,codeperio)
!
!=======================================================================
!
!    "m o d i f p e r i o": Modify the numbering to take periodic
!                           boundary conditions into account
!
!=======================================================================
!
  use spela202
  use codebord

  implicit none

  integer ibool(nxgll,nygll,nspec),iboolori(nxgll,nygll,nspec)
  integer codeperio(4,nelemperio)

  integer n,num2,i2,j2,ipos,ipos2,iloc,jloc,kloc
  integer iloc1,jloc1,iloc2,jloc2
  integer iother1,jother1,iother2,jother2,iboolother1,iboolother2
  integer ispec,ix,iy,nedgeloc,nedgeother,numelem

!
!-----------------------------------------------------------------------
!

  print *
  print *
  print *,'Modifying global numbering for periodic boundaries...'
  print *

! sauvegarder ancienne numerotation pour representations graphiques
  iboolori = ibool

    do n=1,nelemperio

      numelem    = codeperio(1,n)
      nedgeloc   = codeperio(2,n)
      num2       = codeperio(3,n)
      nedgeother = codeperio(4,n)

!
!---- point a l'interieur d'une arete, modifier dans arete correspondante
!

! obtenir la numerotation dans l'autre element
! maillage conforme donc on doit supposer que nxgll == nygll

! modifier tout l'interieur de l'arete
  do kloc = 2,nxgll-1

! calculer l'abscisse le long de l'arete de depart
          select case (nedgeloc)
            case(1)
              iloc = kloc
              jloc = 1
              ipos = iloc
            case(2)
              iloc = nxgll
              jloc = kloc
              ipos = jloc
            case(3)
              iloc = kloc
              jloc = nygll
              ipos = nxgll - iloc + 1
            case(4)
              iloc = 1
              jloc = kloc
              ipos = nygll - jloc + 1
          end select

! calculer l'abscisse le long de l'arete d'arrivee
! topologie du maillage coherente, donc sens de parcours des aretes opposes

          ipos2 = nxgll - ipos + 1

! calculer les coordonnees reelles dans l'element d'arrivee
          select case (nedgeother)
            case(1)
              i2 = ipos2
              j2 = 1
            case(2)
              i2 = nxgll
              j2 = ipos2
            case(3)
              i2 = nxgll - ipos2 + 1
              j2 = nygll
            case(4)
              i2 = 1
              j2 = nygll - ipos2 + 1
          end select

! implementer la periodicite en affectant le meme numero global
          ibool(i2,j2,num2) = ibool(iloc,jloc,numelem)

  enddo

!
!---- cas particulier des coins, recherche sur tous les coins du maillage
!

! determiner les deux coins delimitant l'arete de depart
          select case (nedgeloc)
            case(1)
              iloc1 = 1
              jloc1 = 1
              iloc2 = nxgll
              jloc2 = 1
            case(2)
              iloc1 = nxgll
              jloc1 = 1
              iloc2 = nxgll
              jloc2 = nygll
            case(3)
              iloc1 = nxgll
              jloc1 = nygll
              iloc2 = 1
              jloc2 = nygll
            case(4)
              iloc1 = 1
              jloc1 = nygll
              iloc2 = 1
              jloc2 = 1
          end select

! determiner les deux coins delimitant l'arete d'arrivee
          select case (nedgeother)
            case(1)
              iother1 = 1
              jother1 = 1
              iother2 = nxgll
              jother2 = 1
            case(2)
              iother1 = nxgll
              jother1 = 1
              iother2 = nxgll
              jother2 = nygll
            case(3)
              iother1 = nxgll
              jother1 = nygll
              iother2 = 1
              jother2 = nygll
            case(4)
              iother1 = 1
              jother1 = nygll
              iother2 = 1
              jother2 = 1
          end select

  iboolother1 = ibool(iother1,jother1,num2)
  iboolother2 = ibool(iother2,jother2,num2)

! rechercher correspondants de ces deux coins parmi autres coins du maillage
  do ispec = 1,nspec

  if(ispec /= numelem) then
    do ix = 1,nxgll,nxgll-1
      do iy = 1,nygll,nygll-1

! affecter le meme numero global en tenant compte du sens inverse des aretes
      if(ibool(ix,iy,ispec) == iboolother2) &
          ibool(ix,iy,ispec) = ibool(iloc1,jloc1,numelem)

      if(ibool(ix,iy,ispec) == iboolother1) &
          ibool(ix,iy,ispec) = ibool(iloc2,jloc2,numelem)

      enddo
    enddo
  endif

  enddo

  enddo

  return
  end subroutine modifperio

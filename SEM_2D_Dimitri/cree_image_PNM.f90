
!========================================================================
!
!                   S P E C F E M 2 D  Version 5.1
!                   ------------------------------
!
!                         Dimitri Komatitsch
!          Universite de Pau et des Pays de l'Adour, France
!
!                          (c) January 2005
!
!========================================================================

  subroutine cree_image_PNM(donnees_image_PNM_2D,iglob_image_PNM_2D,NX,NY,it,cutvect)

! routine d'affichage du deplacement sous forme d'image en couleurs

! pour voir les snapshots : display image*.pnm
! pour les convertir en autre format : convert image0001.pnm image0001.jpg

  implicit none

  include "constants.h"

  integer NX,NY,it

  double precision cutvect

  integer, dimension(NX,NY) :: iglob_image_PNM_2D

  double precision, dimension(NX,NY) :: donnees_image_PNM_2D

  integer ix,iy

  double precision amplitude_max

  character(len=100) nom_fichier,system_command

  double precision valeur_normalisee
  integer :: R, G, B

! ouverture du fichier image
  write(nom_fichier,"('OUTPUT_FILES/image',i6.6,'.pnm')") it

! ouvrir le fichier
  open(unit=27, file=nom_fichier, status='unknown')

  write(27,"('P3')") ! ecrire P3 = format d'image PNM

  write(27,*) NX,NY ! ecrire la taille
  write(27,*) '255' ! nombre de nuances

! calculer l'amplitude maximum
  amplitude_max = maxval(abs(donnees_image_PNM_2D))

! supprimer les petites amplitudes considerees comme du bruit
  where(abs(donnees_image_PNM_2D) < amplitude_max * cutvect) donnees_image_PNM_2D = 0.d0

! dans le format PNM, l'image commence par le coin en haut a gauche
  do iy=NY,1,-1
    do ix=1,NX

! regarder si le pixel est defini ou non (au dessus de la topographie par exemple)
      if(iglob_image_PNM_2D(ix,iy) == -1) then

! utiliser couleur bleu ciel pour afficher les zones non definies situees au dessus de la topo
        R = 204
        G = 255
        B = 255

      else

! definir les donnees comme etant le deplacement normalise entre [-1:1]
! et converti a l'entier le plus proche
! en se rappelant que l'amplitude peut etre negative
        valeur_normalisee = donnees_image_PNM_2D(ix,iy) / amplitude_max

! supprimer valeurs en dehors de [-1:+1]
        if(valeur_normalisee < -1.d0) valeur_normalisee = -1.d0
        if(valeur_normalisee > 1.d0) valeur_normalisee = 1.d0

! utiliser rouge si deplacement positif, bleu si negatif, pas de vert
        if(valeur_normalisee >= 0.d0) then
          R = nint(255.d0*valeur_normalisee**POWER_DISPLAY_PNM)
          G = 0
          B = 0
        else
          R = 0
          G = 0
          B = nint(255.d0*abs(valeur_normalisee)**POWER_DISPLAY_PNM)
        endif

      endif

! ecrire l'image en couleur
      write(27,"(i3,' ',i3,' ',i3)") R,G,B

    enddo
  enddo

! fermer le fichier
  close(27)

! open image file and create system command to convert image to more convenient format
  write(system_command,"('cd OUTPUT_FILES ; convert image',i6.6, &
          & '.pnm image',i6.6,'.gif ; rm image',i6.6,'.pnm')") it,it,it

! call the system to convert image to GIF
  call system(system_command)

  end subroutine cree_image_PNM


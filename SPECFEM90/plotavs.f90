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

  subroutine plotavs(displ,coord,kmato,ibool,it)

!
! routine sauvegarde fichier AVS
!

  use constspec
  use mesh01
  use spela202

  implicit none

  integer kmato(nspec)
  integer ibool(nxgll,nygll,nspec)
  double precision displ(ndime,npoin),coord(ndime,npoin)
  integer it

  integer icell,i,j,ispec,iavsfile,ip
  double precision rmaxnorm
  character(len=40) name

  print *,'Entering AVS file generation...'

! file number for AVS output
  iavsfile = 34

!---- ouverture du fichier AVS
  write(name,222) it
  open(unit=iavsfile,file=name,status='unknown')
 222 format('avs',i5.5,'.inp')

! nb de noeuds, de cellules, de donnees par cellule
  write(iavsfile,180) npoin,nspec*(nxgll-1)*(nygll-1),1,0,0

! numero et coordonnees des points du maillage (3D fictif avec coordonnee nulle)
  do ip=1,npoin
    write(iavsfile,200) ip,coord(1,ip),coord(2,ip)
  enddo

! numero et topologie des cellules
  icell = 0
      do ispec=1,nspec
            do i=1,nxgll-1
            do j=1,nxgll-1

      icell = icell + 1
      write(iavsfile,210) icell,kmato(ispec),ibool(i,j+1,ispec), &
               ibool(i,j,ispec),ibool(i+1,j,ispec),ibool(i+1,j+1,ispec)

            enddo
            enddo
      enddo

! structure data vector et labels bidons
      write(iavsfile,*) ' 1 1'
      write(iavsfile,*) ' Label1, Label2'

! donnees aux noeuds (norme du vecteur deplacement, normalisee a 1)
  rmaxnorm = maxval(sqrt(displ(1,:)**2 + displ(2,:)**2))
  do ip=1,npoin
    write(iavsfile,205) ip,sqrt(displ(1,ip)**2 + displ(2,ip)**2)/rmaxnorm
  enddo

  print *,'Max norme dans output AVS = ',rmaxnorm

  close(iavsfile)

  print *,'End of AVS file generation...'

180 format(i6,1x,i6,1x,i3,1x,i3,1x,i3)
200 format(i6,1x,e12.5,' 0. ',e12.5)
205 format(i6,1x,e12.5)
210 format(i6,1x,i4,' quad ',i6,1x,i6,1x,i6,1x,i6)

  end subroutine plotavs

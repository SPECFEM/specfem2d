!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module arraydir
!
!=======================================================================
!
!     "arraydir" : for directory of dynamically allocated arrays
!      ----------
!
!=======================================================================
!
  implicit none

  integer, parameter :: iinteg = 1, isngl = 2, idouble = 3
  integer, parameter :: maxnbarrays = 250
  integer, dimension(maxnbarrays), save :: arraysizes,arraytypes
  character(len=12), dimension(maxnbarrays), save :: arraynames
  integer, save :: nbarrays

  end module arraydir
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module captio
!
!=======================================================================
!     "captio" :
!      ----------
!
!=======================================================================
!
  implicit none

  character(len=50), save :: stitle
  character(len=80), save :: jtitle

  end module captio
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module codebord
!
!=======================================================================
!
!  "codebord" : Code bords absorbants et periodiques
!   --------
!
!=======================================================================
!
  implicit none

  integer, parameter :: ihaut   = 1
  integer, parameter :: ibas    = 2
  integer, parameter :: igauche = 3
  integer, parameter :: idroite = 4

! --- code des numeros d'aretes pour les bords absorbants
  integer, parameter :: iaretebas    = 1
  integer, parameter :: iaretedroite = 2
  integer, parameter :: iaretehaut   = 3
  integer, parameter :: iaretegauche = 4

  end module codebord
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module constspec
!
!=======================================================================
!     "constspec" :
!      ----------
!
!=======================================================================
!
  implicit none

  logical, save :: display,ignuplot,interpol,sismos, &
      ivectplot,imeshvect,isymbols,simuornot, &
      imodelvect,iboundvect,initialfield,usletter, &
      ireadmodel,ioutputgrid,iavs

  integer, save :: nrec,isamp,itaff,itfirstaff, &
      icolor,inumber,isubsamp,nrec1,nrec2,irepr,n1ana,n2ana, &
      isismostype,ivecttype,iaffinfo

  double precision, save :: anglerec,anglerec2, &
      cutvect,scalex,scalez,angle,rapport, &
      sizex,sizez,orig_x,orig_z,rapp_page, &
      sizemax,dispmax,factorana,factorxsu,xmin,zmin

  double precision, parameter :: centim = 28.5d0

  end module constspec
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module defpi
!
!=======================================================================
!     "defpi" : Define the constant number pi
!      -----
!
!=======================================================================
!
  implicit none

  double precision, parameter :: pi = 3.141592653589793d0

  end module defpi
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module energie
!
!=======================================================================
!     "energie" :
!      ----------
!
!=======================================================================
!
  implicit none

  integer, save :: ienergy
  logical, save :: compenergy

  end module energie
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module infos
!
!=======================================================================
!     "infos" :
!      ------
!
!=======================================================================
!
  implicit none

  integer, save :: iecho,iexec

  end module infos
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module iounit
!
!=======================================================================
!     "iounit" :
!      ----------
!
!=======================================================================
!
  implicit none

  integer, save :: iin, iout

  end module iounit
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module label1
!
!=======================================================================
!
!     "label1" : Coordinate labels
!      --------
!
!=======================================================================
!
  implicit none

  character(len=5), dimension(3), save :: labelc

  end module label1
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module loadft
!
!=======================================================================
!     "loadft" :
!      ----------
!
!=======================================================================
!
  implicit none

  integer, save :: nltfl

  end module loadft
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module mesh01
!
!=======================================================================
!     "mesh01" :
!      ----------
!
!=======================================================================
!
  implicit none

  integer, save :: npoin,ndofn,ndime,npgeo

  end module mesh01
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module palette
!
!=======================================================================
!     "palette" :
!      ----------
!
!=======================================================================
!
  implicit none

  integer, parameter :: maxcolors = 128

  double precision, dimension(maxcolors), save :: red,green,blue

  end module palette
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module spela202
!
!=======================================================================
!     "spela202" :
!      ----------
!
!=======================================================================
!
  implicit none

  integer, save :: numat,ngnod,nxgll,nygll,nspec,iptsdisp,nelemabs,nelemperio

  end module spela202
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module timeparams
!
!=======================================================================
!     "timeparams" :
!      ----------
!
!=======================================================================
!
  implicit none

  double precision, save :: deltat,time
  integer, save :: ncycl,niter

  end module timeparams
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module verifs
!
!=======================================================================
!     "verifs" :
!      ----------
!
!=======================================================================
!
  implicit none

  double precision, save :: rsizemin,rsizemax,cpoverdxmin,cpoverdxmax, &
  rlamdaSmin,rlamdaSmax,rlamdaPmin,rlamdaPmax,valseuil,freqmaxrep,vpmin,vpmax

  end module verifs
!
!=====================================================================
!
!                  S p e c f e m  V e r s i o n  3 . 0
!                  -----------------------------------
!
!               Dimitri Komatitsch and Jean-Pierre Vilotte
!                        Departement de Sismologie
!       (c) Institut de Physique du Globe de Paris, Octobre 1997
!
!=====================================================================
!
  module vparams
!
!=======================================================================
!     "vparams" :
!      --------
!
!=======================================================================
!
  implicit none

  double precision, save :: cp1,cs1,rho1,cp2,cs2,rho2,xt1,zt1,xt2,zt2

  end module vparams

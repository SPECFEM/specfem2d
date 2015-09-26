
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

  subroutine plotpost()

!
! PostScript display routine
!

#ifdef USE_MPI
  use mpi
#endif

  use specfem_par, only: vector_field_display,coord,vpext,x_source,z_source,st_xval,st_zval,it,deltat,coorg, &
                         xinterp,zinterp,shape2D_display,Uxinterp,Uzinterp,flagrange,density,porosity,tortuosity,&
                         AXISYM,is_on_the_axis,flagrange_GLJ, &
                         poroelastcoef,knods,kmato,ibool, &
                         numabs,codeabs,typeabs,anyabs,nelem_acoustic_surface, acoustic_edges, &
                         simulation_title,nglob,vpImin,vpImax,nrec,NSOURCES, &
                         colors,numbers,subsamp_postscript,imagetype_postscript,interpol,meshvect,modelvect, &
                         boundvect,assign_external_model,cutsnaps,sizemax_arrows,nelemabs,pointsdisp, &
                         nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poro,coupled_elastic_poro, &
                         any_acoustic,any_poroelastic,plot_lowerleft_corner_only, &
                         fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges, &
                         fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge,num_fluid_poro_edges, &
                         solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,num_solid_poro_edges, &
                         poroelastic,myrank,nproc,ier, &
                         coorg_send_ps_velocity_model,RGB_send_ps_velocity_model, &
                         coorg_recv_ps_velocity_model,RGB_recv_ps_velocity_model,&
                         coorg_send_ps_element_mesh,color_send_ps_element_mesh, &
                         coorg_recv_ps_element_mesh,color_recv_ps_element_mesh, &
                         coorg_send_ps_abs,coorg_recv_ps_abs, &
                         coorg_send_ps_free_surface,coorg_recv_ps_free_surface, &
                         coorg_send_ps_vector_field,coorg_recv_ps_vector_field,US_LETTER,is_PML

  implicit none

  include "constants.h"


! color palette
  integer, parameter :: NUM_COLORS = 236
  double precision, dimension(NUM_COLORS) :: red,green,blue

  double precision, dimension(:,:), allocatable  :: coorg_send
  double precision, dimension(:,:), allocatable  :: coorg_recv
  integer k,j,ispec,material,is,ir,imat,icol,l,line_length
  integer index_char,ii,ipoin,in,nnum,inum,ideb,ifin,iedge
! for the file name
  character(len=100) :: file_name
  integer  :: buffer_offset, RGB_offset
  double precision convert,x1,cpIloc,xa,za,xb,zb,lambdaplus2mu,denst
  double precision z1,x2,z2,d,d1,d2,dummy,theta,thetaup,thetadown
  double precision :: cpIsquare
  double precision :: mul_s,kappal_s,rhol_s
  double precision :: kappal_f,rhol_f
  double precision :: mul_fr,kappal_fr,phil,tortl
  double precision ratio_page,dispmax,xmin,zmin
  logical :: anyabs_glob, coupled_acoustic_elastic_glob, coupled_acoustic_poro_glob, &
             coupled_elastic_poro_glob

  integer  :: nspec_recv,num_spec,pointsdisp_loop
  integer  :: nb_coorg_per_elem, nb_color_per_elem,i,iproc
! to suppress useless white spaces in postscript lines
  character(len=100) :: postscript_line
  character(len=1), dimension(100) :: ch1,ch2
  equivalence (postscript_line,ch1)
  logical :: first

  double precision :: afactor,bfactor,cfactor,D_biot,H_biot,C_biot,M_biot,rhol_bar
  double precision xmax,zmax,height,xw,zw,usoffset,sizex,sizez,timeval
#ifdef USE_MPI
  double precision  :: xmin_glob, xmax_glob, zmin_glob, zmax_glob
  double precision  :: dispmax_glob
#endif

#ifndef USE_MPI
! this to avoid warnings by the compiler about unused variables in the case
! of a serial code, therefore use them once and do nothing: just set them to zero
  nspec_recv = 0
  nb_coorg_per_elem = 0
  nb_color_per_elem = 0
  ier = 0
  num_spec = 0
  iproc = nproc
  coorg_recv_ps_velocity_model = 0
  RGB_recv_ps_velocity_model = 0
  coorg_recv_ps_element_mesh = 0
  color_recv_ps_element_mesh = 0
  coorg_recv_ps_abs = 0
  coorg_recv_ps_free_surface = 0
  coorg_recv_ps_vector_field = 0
  allocate(coorg_recv(1,1))
  deallocate(coorg_recv)
#endif

! A4 or US letter paper
  if(US_LETTER) then
    usoffset = 1.75d0
    sizex = 27.94d0
    sizez = 21.59d0
  else
    usoffset = 0.d0
    sizex = 29.7d0
    sizez = 21.d0
  endif

! height of domain numbers in centimeters
  height = 0.05d0

! define color palette in random order

! red
  red(1) = 1.00000000000000
  green(1) = 0.000000000000000E+000
  blue(1) = 0.000000000000000E+000

! DodgerBlue2
  red(2) = 0.109803921568627
  green(2) = 0.525490196078431
  blue(2) = 0.933333333333333

! gold
  red(3) = 1.00000000000000
  green(3) = 0.840000000000000
  blue(3) = 0.000000000000000E+000

! springgreen
  red(4) = 0.000000000000000E+000
  green(4) = 1.00000000000000
  blue(4) = 0.500000000000000

! NavajoWhite
  red(5) = 1.00000000000000
  green(5) = 0.870588235294118
  blue(5) = 0.678431372549020

! SteelBlue3
  red(6) = 0.309803921568627
  green(6) = 0.580392156862745
  blue(6) = 0.803921568627451

! Ivory3
  red(7) = 0.803921568627451
  green(7) = 0.803921568627451
  blue(7) = 0.756862745098039

! SkyBlue4
  red(8) = 0.290196078431373
  green(8) = 0.439215686274510
  blue(8) = 0.545098039215686

! Snow
  red(9) = 0.980392156862745
  green(9) = 0.980392156862745
  blue(9) = 0.980392156862745

! SteelBlue
  red(10) = 0.274509803921569
  green(10) = 0.509803921568627
  blue(10) = 0.705882352941177

! Bisque3
  red(11) = 0.803921568627451
  green(11) = 0.717647058823529
  blue(11) = 0.619607843137255

! Salmon
  red(12) = 0.980392156862745
  green(12) = 0.501960784313725
  blue(12) = 0.447058823529412

! SlateBlue2
  red(13) = 0.478431372549020
  green(13) = 0.403921568627451
  blue(13) = 0.933333333333333

! NavajoWhite2
  red(14) = 0.933333333333333
  green(14) = 0.811764705882353
  blue(14) = 0.631372549019608

! MediumBlue
  red(15) = 0.000000000000000E+000
  green(15) = 0.000000000000000E+000
  blue(15) = 0.803921568627451

! LightCoral
  red(16) = 0.941176470588235
  green(16) = 0.501960784313725
  blue(16) = 0.501960784313725

! FloralWhite
  red(17) = 1.00000000000000
  green(17) = 0.980392156862745
  blue(17) = 0.941176470588235

! Cornsilk3
  red(18) = 0.803921568627451
  green(18) = 0.784313725490196
  blue(18) = 0.694117647058824

! GhostWhite
  red(19) = 0.972549019607843
  green(19) = 0.972549019607843
  blue(19) = 1.00000000000000

! blue
  red(20) = 0.000000000000000E+000
  green(20) = 0.000000000000000E+000
  blue(20) = 1.00000000000000

! Linen
  red(21) = 0.980392156862745
  green(21) = 0.941176470588235
  blue(21) = 0.901960784313726

! peachpuff
  red(22) = 1.00000000000000
  green(22) = 0.850000000000000
  blue(22) = 0.730000000000000

! Cornsilk1
  red(23) = 1.00000000000000
  green(23) = 0.972549019607843
  blue(23) = 0.862745098039216

! LightSalmon
  red(24) = 1.00000000000000
  green(24) = 0.627450980392157
  blue(24) = 0.478431372549020

! DeepSkyBlue1
  red(25) = 0.000000000000000E+000
  green(25) = 0.749019607843137
  blue(25) = 1.00000000000000

! LemonChiffon4
  red(26) = 0.545098039215686
  green(26) = 0.537254901960784
  blue(26) = 0.439215686274510

! PeachPuff1
  red(27) = 1.00000000000000
  green(27) = 0.854901960784314
  blue(27) = 0.725490196078431

! BlanchedAlmond
  red(28) = 1.00000000000000
  green(28) = 0.921568627450980
  blue(28) = 0.803921568627451

! SlateBlue3
  red(29) = 0.411764705882353
  green(29) = 0.349019607843137
  blue(29) = 0.803921568627451

! LightSkyBlue1
  red(30) = 0.690196078431373
  green(30) = 0.886274509803922
  blue(30) = 1.00000000000000

! DarkViolet
  red(31) = 0.580392156862745
  green(31) = 0.000000000000000E+000
  blue(31) = 0.827450980392157

! Azure3
  red(32) = 0.756862745098039
  green(32) = 0.803921568627451
  blue(32) = 0.803921568627451

! LavenderBlush3
  red(33) = 0.803921568627451
  green(33) = 0.756862745098039
  blue(33) = 0.772549019607843

! Honeydew1
  red(34) = 0.941176470588235
  green(34) = 1.00000000000000
  blue(34) = 0.941176470588235

! Ivory2
  red(35) = 0.933333333333333
  green(35) = 0.933333333333333
  blue(35) = 0.878431372549020

! RosyBrown
  red(36) = 0.737254901960784
  green(36) = 0.560784313725490
  blue(36) = 0.560784313725490

! Thistle
  red(37) = 0.847058823529412
  green(37) = 0.749019607843137
  blue(37) = 0.847058823529412

! Orange
  red(38) = 1.00000000000000
  green(38) = 0.647058823529412
  blue(38) = 0.000000000000000E+000

! DarkSeaGreen
  red(39) = 0.560784313725490
  green(39) = 0.737254901960784
  blue(39) = 0.560784313725490

! Moccasin
  red(40) = 1.00000000000000
  green(40) = 0.894117647058824
  blue(40) = 0.709803921568627

! DeepSkyBlue2
  red(41) = 0.000000000000000E+000
  green(41) = 0.698039215686274
  blue(41) = 0.933333333333333

! SlateGray4
  red(42) = 0.423529411764706
  green(42) = 0.482352941176471
  blue(42) = 0.545098039215686

! Beige
  red(43) = 0.960784313725490
  green(43) = 0.960784313725490
  blue(43) = 0.862745098039216

! Gold
  red(44) = 1.00000000000000
  green(44) = 0.843137254901961
  blue(44) = 0.000000000000000E+000

! SlateBlue
  red(45) = 0.415686274509804
  green(45) = 0.352941176470588
  blue(45) = 0.803921568627451

! SteelBlue1
  red(46) = 0.388235294117647
  green(46) = 0.721568627450980
  blue(46) = 1.00000000000000

! SaddleBrown
  red(47) = 0.545098039215686
  green(47) = 0.270588235294118
  blue(47) = 7.450980392156863E-002

! Pink
  red(48) = 1.00000000000000
  green(48) = 0.752941176470588
  blue(48) = 0.796078431372549

! Black
  red(49) = 0.000000000000000E+000
  green(49) = 0.000000000000000E+000
  blue(49) = 0.000000000000000E+000

! SlateGrey
  red(50) = 0.439215686274510
  green(50) = 0.501960784313725
  blue(50) = 0.564705882352941

! Ivory
  red(51) = 1.00000000000000
  green(51) = 1.00000000000000
  blue(51) = 0.941176470588235

! OliveDrab
  red(52) = 0.419607843137255
  green(52) = 0.556862745098039
  blue(52) = 0.137254901960784

! Ivory1
  red(53) = 1.00000000000000
  green(53) = 1.00000000000000
  blue(53) = 0.941176470588235

! SkyBlue
  red(54) = 0.529411764705882
  green(54) = 0.807843137254902
  blue(54) = 0.921568627450980

! MistyRose3
  red(55) = 0.803921568627451
  green(55) = 0.717647058823529
  blue(55) = 0.709803921568627

! LimeGreen
  red(56) = 0.196078431372549
  green(56) = 0.803921568627451
  blue(56) = 0.196078431372549

! Purple
  red(57) = 0.627450980392157
  green(57) = 0.125490196078431
  blue(57) = 0.941176470588235

! SkyBlue2
  red(58) = 0.494117647058824
  green(58) = 0.752941176470588
  blue(58) = 0.933333333333333

! Red
  red(59) = 1.00000000000000
  green(59) = 0.000000000000000E+000
  blue(59) = 0.000000000000000E+000

! DarkKhaki
  red(60) = 0.741176470588235
  green(60) = 0.717647058823529
  blue(60) = 0.419607843137255

! MediumTurquoise
  red(61) = 0.282352941176471
  green(61) = 0.819607843137255
  blue(61) = 0.800000000000000

! Grey
  red(62) = 0.745098039215686
  green(62) = 0.745098039215686
  blue(62) = 0.745098039215686

! Coral
  red(63) = 1.00000000000000
  green(63) = 0.498039215686275
  blue(63) = 0.313725490196078

! NavajoWhite4
  red(64) = 0.545098039215686
  green(64) = 0.474509803921569
  blue(64) = 0.368627450980392

! SlateBlue4
  red(65) = 0.278431372549020
  green(65) = 0.235294117647059
  blue(65) = 0.545098039215686

! RoyalBlue4
  red(66) = 0.152941176470588
  green(66) = 0.250980392156863
  blue(66) = 0.545098039215686

! YellowGreen
  red(67) = 0.603921568627451
  green(67) = 0.803921568627451
  blue(67) = 0.196078431372549

! DeepSkyBlue3
  red(68) = 0.000000000000000E+000
  green(68) = 0.603921568627451
  blue(68) = 0.803921568627451

! goldenrod
  red(69) = 0.854901960784314
  green(69) = 0.647058823529412
  blue(69) = 0.125490196078431

! AntiqueWhite4
  red(70) = 0.545098039215686
  green(70) = 0.513725490196078
  blue(70) = 0.470588235294118

! lemonchiffon
  red(71) = 1.00000000000000
  green(71) = 0.980000000000000
  blue(71) = 0.800000000000000

! GreenYellow
  red(72) = 0.678431372549020
  green(72) = 1.00000000000000
  blue(72) = 0.184313725490196

! LightSlateGray
  red(73) = 0.466666666666667
  green(73) = 0.533333333333333
  blue(73) = 0.600000000000000

! RoyalBlue
  red(74) = 0.254901960784314
  green(74) = 0.411764705882353
  blue(74) = 0.882352941176471

! DarkGreen
  red(75) = 0.000000000000000E+000
  green(75) = 0.392156862745098
  blue(75) = 0.000000000000000E+000

! NavajoWhite3
  red(76) = 0.803921568627451
  green(76) = 0.701960784313725
  blue(76) = 0.545098039215686

! Azure1
  red(77) = 0.941176470588235
  green(77) = 1.00000000000000
  blue(77) = 1.00000000000000

! PowderBlue
  red(78) = 0.690196078431373
  green(78) = 0.878431372549020
  blue(78) = 0.901960784313726

! slateblue
  red(79) = 0.420000000000000
  green(79) = 0.350000000000000
  blue(79) = 0.800000000000000

! MediumOrchid
  red(80) = 0.729411764705882
  green(80) = 0.333333333333333
  blue(80) = 0.827450980392157

! turquoise
  red(81) = 0.250000000000000
  green(81) = 0.880000000000000
  blue(81) = 0.820000000000000

! Snow1
  red(82) = 1.00000000000000
  green(82) = 0.980392156862745
  blue(82) = 0.980392156862745

! violet
  red(83) = 0.930000000000000
  green(83) = 0.510000000000000
  blue(83) = 0.930000000000000

! DeepPink
  red(84) = 1.00000000000000
  green(84) = 7.843137254901961E-002
  blue(84) = 0.576470588235294

! MistyRose4
  red(85) = 0.545098039215686
  green(85) = 0.490196078431373
  blue(85) = 0.482352941176471

! PeachPuff3
  red(86) = 0.803921568627451
  green(86) = 0.686274509803922
  blue(86) = 0.584313725490196

! MediumSeaGreen
  red(87) = 0.235294117647059
  green(87) = 0.701960784313725
  blue(87) = 0.443137254901961

! Honeydew4
  red(88) = 0.513725490196078
  green(88) = 0.545098039215686
  blue(88) = 0.513725490196078

! Tan
  red(89) = 0.823529411764706
  green(89) = 0.705882352941177
  blue(89) = 0.549019607843137

! DarkGoldenrod
  red(90) = 0.721568627450980
  green(90) = 0.525490196078431
  blue(90) = 4.313725490196078E-002

! Blue2
  red(91) = 0.000000000000000E+000
  green(91) = 0.000000000000000E+000
  blue(91) = 0.933333333333333

! Maroon
  red(92) = 0.690196078431373
  green(92) = 0.188235294117647
  blue(92) = 0.376470588235294

! LightSkyBlue3
  red(93) = 0.552941176470588
  green(93) = 0.713725490196078
  blue(93) = 0.803921568627451

! LemonChiffon2
  red(94) = 0.933333333333333
  green(94) = 0.913725490196078
  blue(94) = 0.749019607843137

! Snow3
  red(95) = 0.803921568627451
  green(95) = 0.788235294117647
  blue(95) = 0.788235294117647

! Ivory4
  red(96) = 0.545098039215686
  green(96) = 0.545098039215686
  blue(96) = 0.513725490196078

! AntiqueWhite3
  red(97) = 0.803921568627451
  green(97) = 0.752941176470588
  blue(97) = 0.690196078431373

! Bisque4
  red(98) = 0.545098039215686
  green(98) = 0.490196078431373
  blue(98) = 0.419607843137255

! Snow2
  red(99) = 0.933333333333333
  green(99) = 0.913725490196078
  blue(99) = 0.913725490196078

! SlateGray1
  red(100) = 0.776470588235294
  green(100) = 0.886274509803922
  blue(100) = 1.00000000000000

! Seashell2
  red(101) = 0.933333333333333
  green(101) = 0.898039215686275
  blue(101) = 0.870588235294118

! Aquamarine
  red(102) = 0.498039215686275
  green(102) = 1.00000000000000
  blue(102) = 0.831372549019608

! SlateGray2
  red(103) = 0.725490196078431
  green(103) = 0.827450980392157
  blue(103) = 0.933333333333333

! White
  red(104) = 1.00000000000000
  green(104) = 1.00000000000000
  blue(104) = 1.00000000000000

! LavenderBlush
  red(105) = 1.00000000000000
  green(105) = 0.941176470588235
  blue(105) = 0.960784313725490

! DodgerBlue3
  red(106) = 9.411764705882353E-002
  green(106) = 0.454901960784314
  blue(106) = 0.803921568627451

! RoyalBlue3
  red(107) = 0.227450980392157
  green(107) = 0.372549019607843
  blue(107) = 0.803921568627451

! LightYellow
  red(108) = 1.00000000000000
  green(108) = 1.00000000000000
  blue(108) = 0.878431372549020

! DeepSkyBlue
  red(109) = 0.000000000000000E+000
  green(109) = 0.749019607843137
  blue(109) = 1.00000000000000

! AntiqueWhite2
  red(110) = 0.933333333333333
  green(110) = 0.874509803921569
  blue(110) = 0.800000000000000

! CornflowerBlue
  red(111) = 0.392156862745098
  green(111) = 0.584313725490196
  blue(111) = 0.929411764705882

! PeachPuff4
  red(112) = 0.545098039215686
  green(112) = 0.466666666666667
  blue(112) = 0.396078431372549

! SpringGreen
  red(113) = 0.000000000000000E+000
  green(113) = 1.00000000000000
  blue(113) = 0.498039215686275

! Honeydew
  red(114) = 0.941176470588235
  green(114) = 1.00000000000000
  blue(114) = 0.941176470588235

! Honeydew2
  red(115) = 0.878431372549020
  green(115) = 0.933333333333333
  blue(115) = 0.878431372549020

! LightSeaGreen
  red(116) = 0.125490196078431
  green(116) = 0.698039215686274
  blue(116) = 0.666666666666667

! NavyBlue
  red(117) = 0.000000000000000E+000
  green(117) = 0.000000000000000E+000
  blue(117) = 0.501960784313725

! Azure4
  red(118) = 0.513725490196078
  green(118) = 0.545098039215686
  blue(118) = 0.545098039215686

! MediumAquamarine
  red(119) = 0.400000000000000
  green(119) = 0.803921568627451
  blue(119) = 0.666666666666667

! SkyBlue3
  red(120) = 0.423529411764706
  green(120) = 0.650980392156863
  blue(120) = 0.803921568627451

! LavenderBlush2
  red(121) = 0.933333333333333
  green(121) = 0.878431372549020
  blue(121) = 0.898039215686275

! Bisque1
  red(122) = 1.00000000000000
  green(122) = 0.894117647058824
  blue(122) = 0.768627450980392

! DarkOrange
  red(123) = 1.00000000000000
  green(123) = 0.549019607843137
  blue(123) = 0.000000000000000E+000

! LightSteelBlue
  red(124) = 0.690196078431373
  green(124) = 0.768627450980392
  blue(124) = 0.870588235294118

! SteelBlue2
  red(125) = 0.360784313725490
  green(125) = 0.674509803921569
  blue(125) = 0.933333333333333

! LemonChiffon3
  red(126) = 0.803921568627451
  green(126) = 0.788235294117647
  blue(126) = 0.647058823529412

! DarkSlateBlue
  red(127) = 0.282352941176471
  green(127) = 0.239215686274510
  blue(127) = 0.545098039215686

! Seashell
  red(128) = 1.00000000000000
  green(128) = 0.960784313725490
  blue(128) = 0.933333333333333

! Firebrick
  red(129) = 0.698039215686274
  green(129) = 0.133333333333333
  blue(129) = 0.133333333333333

! LightGray
  red(130) = 0.827450980392157
  green(130) = 0.827450980392157
  blue(130) = 0.827450980392157

! Blue
  red(131) = 0.000000000000000E+000
  green(131) = 0.000000000000000E+000
  blue(131) = 1.00000000000000

! Bisque2
  red(132) = 0.933333333333333
  green(132) = 0.835294117647059
  blue(132) = 0.717647058823529

! WhiteSmoke
  red(133) = 0.960784313725490
  green(133) = 0.960784313725490
  blue(133) = 0.960784313725490

! SeaGreen
  red(134) = 0.180392156862745
  green(134) = 0.545098039215686
  blue(134) = 0.341176470588235

! Burlywood
  red(135) = 0.870588235294118
  green(135) = 0.721568627450980
  blue(135) = 0.529411764705882

! RoyalBlue2
  red(136) = 0.262745098039216
  green(136) = 0.431372549019608
  blue(136) = 0.933333333333333

! RoyalBlue1
  red(137) = 0.282352941176471
  green(137) = 0.462745098039216
  blue(137) = 1.00000000000000

! SteelBlue4
  red(138) = 0.211764705882353
  green(138) = 0.392156862745098
  blue(138) = 0.545098039215686

! AliceBlue
  red(139) = 0.941176470588235
  green(139) = 0.972549019607843
  blue(139) = 1.00000000000000

! LightSlateBlue
  red(140) = 0.517647058823529
  green(140) = 0.439215686274510
  blue(140) = 1.00000000000000

! MistyRose1
  red(141) = 1.00000000000000
  green(141) = 0.894117647058824
  blue(141) = 0.882352941176471

! SandyBrown
  red(142) = 0.956862745098039
  green(142) = 0.643137254901961
  blue(142) = 0.376470588235294

! DarkOliveGreen
  red(143) = 0.333333333333333
  green(143) = 0.419607843137255
  blue(143) = 0.184313725490196

! Yellow
  red(144) = 1.00000000000000
  green(144) = 1.00000000000000
  blue(144) = 0.000000000000000E+000

! SlateGray3
  red(145) = 0.623529411764706
  green(145) = 0.713725490196078
  blue(145) = 0.803921568627451

! HotPink
  red(146) = 1.00000000000000
  green(146) = 0.411764705882353
  blue(146) = 0.705882352941177

! Violet
  red(147) = 0.933333333333333
  green(147) = 0.509803921568627
  blue(147) = 0.933333333333333

! LightSkyBlue
  red(148) = 0.529411764705882
  green(148) = 0.807843137254902
  blue(148) = 0.980392156862745

! Cornsilk2
  red(149) = 0.933333333333333
  green(149) = 0.909803921568627
  blue(149) = 0.803921568627451

! MidnightBlue
  red(150) = 9.803921568627451E-002
  green(150) = 9.803921568627451E-002
  blue(150) = 0.439215686274510

! AntiqueWhite
  red(151) = 0.980392156862745
  green(151) = 0.921568627450980
  blue(151) = 0.843137254901961

! PaleGreen
  red(152) = 0.596078431372549
  green(152) = 0.984313725490196
  blue(152) = 0.596078431372549

! MedSpringGreen
  red(153) = 0.000000000000000E+000
  green(153) = 0.980392156862745
  blue(153) = 0.603921568627451

! DodgerBlue1
  red(154) = 0.117647058823529
  green(154) = 0.564705882352941
  blue(154) = 1.00000000000000

! Blue3
  red(155) = 0.000000000000000E+000
  green(155) = 0.000000000000000E+000
  blue(155) = 0.803921568627451

! Cyan
  red(156) = 0.000000000000000E+000
  green(156) = 1.00000000000000
  blue(156) = 1.00000000000000

! LemonChiffon
  red(157) = 1.00000000000000
  green(157) = 0.980392156862745
  blue(157) = 0.803921568627451

! mediumorchid
  red(158) = 0.730000000000000
  green(158) = 0.330000000000000
  blue(158) = 0.830000000000000

! Turquoise
  red(159) = 0.250980392156863
  green(159) = 0.878431372549020
  blue(159) = 0.815686274509804

! IndianRed
  red(160) = 0.803921568627451
  green(160) = 0.360784313725490
  blue(160) = 0.360784313725490

! DodgerBlue
  red(161) = 0.117647058823529
  green(161) = 0.564705882352941
  blue(161) = 1.00000000000000

! Seashell3
  red(162) = 0.803921568627451
  green(162) = 0.772549019607843
  blue(162) = 0.749019607843137

! BlueViolet
  red(163) = 0.541176470588235
  green(163) = 0.168627450980392
  blue(163) = 0.886274509803922

! DeepSkyBlue4
  red(164) = 0.000000000000000E+000
  green(164) = 0.407843137254902
  blue(164) = 0.545098039215686

! PaleVioletRed
  red(165) = 0.858823529411765
  green(165) = 0.439215686274510
  blue(165) = 0.576470588235294

! Azure2
  red(166) = 0.878431372549020
  green(166) = 0.933333333333333
  blue(166) = 0.933333333333333

! greenyellow
  red(167) = 0.680000000000000
  green(167) = 1.00000000000000
  blue(167) = 0.180000000000000

! LightGoldenrod
  red(168) = 0.933333333333333
  green(168) = 0.866666666666667
  blue(168) = 0.509803921568627

! MistyRose
  red(169) = 1.00000000000000
  green(169) = 0.894117647058824
  blue(169) = 0.882352941176471

! LightSkyBlue4
  red(170) = 0.376470588235294
  green(170) = 0.482352941176471
  blue(170) = 0.545098039215686

! OrangeRed
  red(171) = 1.00000000000000
  green(171) = 0.270588235294118
  blue(171) = 0.000000000000000E+000

! DimGrey
  red(172) = 0.411764705882353
  green(172) = 0.411764705882353
  blue(172) = 0.411764705882353

! MediumVioletRed
  red(173) = 0.780392156862745
  green(173) = 8.235294117647059E-002
  blue(173) = 0.521568627450980

! DarkSlateGray
  red(174) = 0.184313725490196
  green(174) = 0.309803921568627
  blue(174) = 0.309803921568627

! yellow
  red(175) = 1.00000000000000
  green(175) = 1.00000000000000
  blue(175) = 0.000000000000000E+000

! Plum
  red(176) = 0.866666666666667
  green(176) = 0.627450980392157
  blue(176) = 0.866666666666667

! DarkTurquoise
  red(177) = 0.000000000000000E+000
  green(177) = 0.807843137254902
  blue(177) = 0.819607843137255

! DodgerBlue4
  red(178) = 6.274509803921569E-002
  green(178) = 0.305882352941176
  blue(178) = 0.545098039215686

! Cornsilk
  red(179) = 1.00000000000000
  green(179) = 0.972549019607843
  blue(179) = 0.862745098039216

! SkyBlue1
  red(180) = 0.529411764705882
  green(180) = 0.807843137254902
  blue(180) = 1.00000000000000

! Seashell1
  red(181) = 1.00000000000000
  green(181) = 0.960784313725490
  blue(181) = 0.933333333333333

! lavender
  red(182) = 0.901960784313726
  green(182) = 0.901960784313726
  blue(182) = 0.980392156862745

! Snow4
  red(183) = 0.545098039215686
  green(183) = 0.537254901960784
  blue(183) = 0.537254901960784

! Peru
  red(184) = 0.803921568627451
  green(184) = 0.521568627450980
  blue(184) = 0.247058823529412

! PeachPuff
  red(185) = 1.00000000000000
  green(185) = 0.854901960784314
  blue(185) = 0.725490196078431

! Green
  red(186) = 0.000000000000000E+000
  green(186) = 1.00000000000000
  blue(186) = 0.000000000000000E+000

! Blue1
  red(187) = 0.000000000000000E+000
  green(187) = 0.000000000000000E+000
  blue(187) = 1.00000000000000

! Seashell4
  red(188) = 0.545098039215686
  green(188) = 0.525490196078431
  blue(188) = 0.509803921568627

! dodgerblue
  red(189) = 0.120000000000000
  green(189) = 0.560000000000000
  blue(189) = 1.00000000000000

! MistyRose2
  red(190) = 0.933333333333333
  green(190) = 0.835294117647059
  blue(190) = 0.823529411764706

! Tomato
  red(191) = 1.00000000000000
  green(191) = 0.388235294117647
  blue(191) = 0.278431372549020

! Wheat
  red(192) = 0.960784313725490
  green(192) = 0.870588235294118
  blue(192) = 0.701960784313725

! LightBlue
  red(193) = 0.678431372549020
  green(193) = 0.847058823529412
  blue(193) = 0.901960784313726

! Chocolate
  red(194) = 0.823529411764706
  green(194) = 0.411764705882353
  blue(194) = 0.117647058823529

! Blue4
  red(195) = 0.000000000000000E+000
  green(195) = 0.000000000000000E+000
  blue(195) = 0.545098039215686

! LavenderBlush1
  red(196) = 1.00000000000000
  green(196) = 0.941176470588235
  blue(196) = 0.960784313725490

! Magenta
  red(197) = 1.00000000000000
  green(197) = 0.000000000000000E+000
  blue(197) = 1.00000000000000

! darkturquoise
  red(198) = 0.000000000000000E+000
  green(198) = 0.810000000000000
  blue(198) = 0.820000000000000

! blueviolet
  red(199) = 0.540000000000000
  green(199) = 0.170000000000000
  blue(199) = 0.890000000000000

! MintCream
  red(200) = 0.960784313725490
  green(200) = 1.00000000000000
  blue(200) = 0.980392156862745

! PaleGoldenrod
  red(201) = 0.933333333333333
  green(201) = 0.909803921568627
  blue(201) = 0.666666666666667

! MediumPurple
  red(202) = 0.576470588235294
  green(202) = 0.439215686274510
  blue(202) = 0.858823529411765

! PapayaWhip
  red(203) = 1.00000000000000
  green(203) = 0.937254901960784
  blue(203) = 0.835294117647059

! LavenderBlush4
  red(204) = 0.545098039215686
  green(204) = 0.513725490196078
  blue(204) = 0.525490196078431

! Cornsilk4
  red(205) = 0.545098039215686
  green(205) = 0.533333333333333
  blue(205) = 0.470588235294118

! LtGoldenrodYello
  red(206) = 0.980392156862745
  green(206) = 0.980392156862745
  blue(206) = 0.823529411764706

! limegreen
  red(207) = 0.200000000000000
  green(207) = 0.800000000000000
  blue(207) = 0.200000000000000

! LemonChiffon1
  red(208) = 1.00000000000000
  green(208) = 0.980392156862745
  blue(208) = 0.803921568627451

! DarkOrchid
  red(209) = 0.600000000000000
  green(209) = 0.196078431372549
  blue(209) = 0.800000000000000

! SlateBlue1
  red(210) = 0.513725490196078
  green(210) = 0.435294117647059
  blue(210) = 1.00000000000000

! chartreuse
  red(211) = 0.500000000000000
  green(211) = 1.00000000000000
  blue(211) = 0.000000000000000E+000

! PaleTurquoise
  red(212) = 0.686274509803922
  green(212) = 0.933333333333333
  blue(212) = 0.933333333333333

! NavajoWhite1
  red(213) = 1.00000000000000
  green(213) = 0.870588235294118
  blue(213) = 0.678431372549020

! LightSkyBlue2
  red(214) = 0.643137254901961
  green(214) = 0.827450980392157
  blue(214) = 0.933333333333333

! VioletRed
  red(215) = 0.815686274509804
  green(215) = 0.125490196078431
  blue(215) = 0.564705882352941

! mocassin
  red(216) = 1.00000000000000
  green(216) = 0.890000000000000
  blue(216) = 0.710000000000000

! OldLace
  red(217) = 0.992156862745098
  green(217) = 0.960784313725490
  blue(217) = 0.901960784313726

! deeppink
  red(218) = 1.00000000000000
  green(218) = 8.000000000000000E-002
  blue(218) = 0.580000000000000

! Honeydew3
  red(219) = 0.756862745098039
  green(219) = 0.803921568627451
  blue(219) = 0.756862745098039

! Gainsboro
  red(220) = 0.862745098039216
  green(220) = 0.862745098039216
  blue(220) = 0.862745098039216

! DarkSalmon
  red(221) = 0.913725490196078
  green(221) = 0.588235294117647
  blue(221) = 0.478431372549020

! AntiqueWhite1
  red(222) = 1.00000000000000
  green(222) = 0.937254901960784
  blue(222) = 0.858823529411765

! LightCyan
  red(223) = 0.878431372549020
  green(223) = 1.00000000000000
  blue(223) = 1.00000000000000

! ForestGreen
  red(224) = 0.133333333333333
  green(224) = 0.545098039215686
  blue(224) = 0.133333333333333

! Orchid
  red(225) = 0.854901960784314
  green(225) = 0.439215686274510
  blue(225) = 0.839215686274510

! PeachPuff2
  red(226) = 0.933333333333333
  green(226) = 0.796078431372549
  blue(226) = 0.678431372549020

! LightPink
  red(227) = 1.00000000000000
  green(227) = 0.713725490196078
  blue(227) = 0.756862745098039

! Sienna
  red(228) = 0.627450980392157
  green(228) = 0.321568627450980
  blue(228) = 0.176470588235294

! darkorchid
  red(229) = 0.600000000000000
  green(229) = 0.200000000000000
  blue(229) = 0.800000000000000

! MediumSlateBlue
  red(230) = 0.482352941176471
  green(230) = 0.407843137254902
  blue(230) = 0.933333333333333

! CadetBlue
  red(231) = 0.372549019607843
  green(231) = 0.619607843137255
  blue(231) = 0.627450980392157

! LawnGreen
  red(232) = 0.486274509803922
  green(232) = 0.988235294117647
  blue(232) = 0.000000000000000E+000

! Chartreuse
  red(233) = 0.498039215686275
  green(233) = 1.00000000000000
  blue(233) = 0.000000000000000E+000

! Brown
  red(234) = 0.647058823529412
  green(234) = 0.164705882352941
  blue(234) = 0.164705882352941

! Azure
  red(235) = 0.941176470588235
  green(235) = 1.00000000000000
  blue(235) = 1.00000000000000

! Bisque
  red(236) = 1.00000000000000
  green(236) = 0.894117647058824
  blue(236) = 0.768627450980392

! get minimum and maximum values of mesh coordinates
  xmin = minval(coord(1,:))
  zmin = minval(coord(2,:))
  xmax = maxval(coord(1,:))
  zmax = maxval(coord(2,:))

#ifdef USE_MPI
  call MPI_ALLREDUCE (xmin, xmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (zmin, zmin_glob, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (xmax, xmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  call MPI_ALLREDUCE (zmax, zmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  xmin = xmin_glob
  zmin = zmin_glob
  xmax = xmax_glob
  zmax = zmax_glob
#endif

  if ( myrank == 0 ) then
     write(IOUT,*) 'X min, max = ',xmin,xmax
     write(IOUT,*) 'Z min, max = ',zmin,zmax
  endif

! ratio of physical page size/size of the domain meshed
  ratio_page = min(rpercentz*sizez/(zmax-zmin),rpercentx*sizex/(xmax-xmin)) / 100.d0

! compute the maximum of the norm of the vector
  dispmax = maxval(sqrt(vector_field_display(1,:)**2 + vector_field_display(3,:)**2))
#ifdef USE_MPI
  call MPI_ALLREDUCE (dispmax, dispmax_glob, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ier)
  dispmax = dispmax_glob
#endif
  if ( myrank == 0 ) then
     write(IOUT,*) 'Max norm = ',dispmax
  endif

!
!---- open PostScript file
!
  if ( myrank == 0 ) then
  write(file_name,"('OUTPUT_FILES/vect',i7.7,'.ps')") it
  open(unit=24,file=file_name,status='unknown')

!
!---- write PostScript header
!
  write(24,10) simulation_title
  write(24,*) '/CM {28.5 mul} def'
  write(24,*) '/LR {rlineto} def'
  write(24,*) '/LT {lineto} def'
  write(24,*) '/L {lineto} def'
  write(24,*) '/MR {rmoveto} def'
  write(24,*) '/MV {moveto} def'
  write(24,*) '/M {moveto} def'
  write(24,*) '/ST {stroke} def'
  write(24,*) '/CP {closepath} def'
  write(24,*) '/RG {setrgbcolor} def'
  write(24,*) '/GF {gsave fill grestore} def'
  write(24,*) '% different useful symbols'
  write(24,*) '/Point {2 0 360 arc CP 0 setgray fill} def'
  write(24,*) '/VDot {-0.75 -1.5 MR 1.5 0 LR 0 3. LR -1.5 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/HDot {-1.5 -0.75 MR 3. 0 LR 0 1.5 LR -3. 0 LR'
  write(24,*) 'CP fill} def'
  write(24,*) '/Cross {gsave 0.05 CM setlinewidth'
  write(24,*) 'gsave 3 3 MR -6. -6. LR ST grestore'
  write(24,*) 'gsave 3 -3 MR -6. 6. LR ST grestore'
  write(24,*) '0.01 CM setlinewidth} def'
  write(24,*) '/SmallLine {MV 0.07 CM 0 rlineto} def'
  write(24,*) '/Diamond {gsave 0.05 CM setlinewidth 0 4.2 MR'
  write(24,*) '-3 -4.2 LR 3 -4.2 LR 3 4.2 LR CP ST'
  write(24,*) 'grestore 0.01 CM setlinewidth} def'
  write(24,*) '%'
  write(24,*) '% gray levels for the velocity model'
  write(24,*) '/BK {setgray fill} def'
  write(24,*) '% black and white version'
  write(24,*) '%/BK {pop 1 setgray fill} def'
  write(24,*) '%'
  write(24,*) '% magenta for vectors'
  write(24,*) '/Colvects {0 setlinewidth 1. 0. 1. RG} def'
  write(24,*) '% black and white version'
  write(24,*) '%/Colvects {0 setlinewidth 0. setgray} def'
  write(24,*) '%'
  write(24,*) '% chartreuse for macrobloc mesh'
  write(24,*) '/Colmesh {0 setlinewidth 0.5 1. 0. RG} def'
  write(24,*) '% black and white version'
  write(24,*) '%/Colmesh {0 setlinewidth 0. setgray} def'
  write(24,*) '%'
  write(24,*) '% cyan for sources and receivers'
  write(24,*) '/Colreceiv {0. 1. 1. RG} def'
  write(24,*) '% black and white version'
  write(24,*) '%/Colreceiv {0. setgray} def'
  write(24,*) '%'
  write(24,*) '% macro to draw an arrow'
  write(24,*) '/F {MV LR gsave LR ST grestore LR ST} def'
  write(24,*) '% macro to draw the contour of the elements'
  write(24,*) '/CO {M counttomark 2 idiv {L} repeat cleartomark CP} def'
  write(24,*) '%'
  write(24,*) '0 setlinewidth'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.35 CM scalefont setfont'
  write(24,*) '%'
  write(24,*) '/vshift ',-height/2,' CM def'
  write(24,*) '/Rshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop neg vshift MR show } def'
  write(24,*) '/Cshow { currentpoint stroke MV'
  write(24,*) 'dup stringwidth pop -2 div vshift MR show } def'
  write(24,*) '/fN {/Helvetica-Bold findfont ',height,' CM scalefont setfont} def'
  write(24,*) '%'
  write(24,*) 'gsave newpath 90 rotate'
  write(24,*) '0 ',-sizez,' CM translate 1. 1. scale'
  write(24,*) '% uncomment this to zoom on parts of the mesh'
  write(24,*) '% -32 CM -21 CM translate 3. 3. scale'
  write(24,*) '% -52 CM -24 CM translate 4. 4. scale'
  write(24,*) '%'

!
!--- write captions of PostScript figure
!
  write(24,*) '0 setgray'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.5 CM scalefont setfont'

  write(24,*) '24. CM 1.2 CM MV'
  write(24,610) usoffset,it
  write(24,*) '%'

  write(24,*) '24. CM 1.95 CM MV'
  timeval = it*deltat
  if(timeval >= 1.d-3 .and. timeval < 1000.d0) then
    write(24,600) usoffset,timeval
  else
    write(24,601) usoffset,timeval
  endif
  write(24,*) '%'
  write(24,*) '24. CM 2.7 CM MV'
  write(24,640) usoffset,dispmax
  write(24,*) '%'
  write(24,*) '24. CM 3.45 CM MV'
  write(24,620) usoffset,cutsnaps*100.d0

  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.6 CM scalefont setfont'
  if(colors == 1) write(24,*) '.4 .9 .9 setrgbcolor'
  write(24,*) '11 CM 1.1 CM MV'
  write(24,*) '(X axis) show'
  write(24,*) '%'
  write(24,*) '1.4 CM 9.5 CM MV'
  write(24,*) 'currentpoint gsave translate 90 rotate 0 0 moveto'
  write(24,*) '(Z axis) show'
  write(24,*) 'grestore'
  write(24,*) '%'
  write(24,*) '/Times-Roman findfont'
  write(24,*) '.7 CM scalefont setfont'
  if(colors == 1) write(24,*) '.8 0 .8 setrgbcolor'
  write(24,*) '24.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  if(imagetype_postscript == 1) then
    write(24,*) '(Displacement vector field) show'
  else if(imagetype_postscript == 2) then
    write(24,*) '(Velocity vector field) show'
  else if(imagetype_postscript == 3) then
    write(24,*) '(Acceleration vector field) show'
  else
    call exit_MPI('Bad field code in PostScript display')
  endif
  write(24,*) 'grestore'
  write(24,*) '25.35 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
  write(24,*) '(',simulation_title,') show'
  write(24,*) 'grestore'
  write(24,*) '26.45 CM 18.9 CM MV'
  write(24,*) usoffset,' CM 2 div neg 0 MR'
  write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'

  if(coupled_acoustic_elastic) then
    write(24,*) '(Coupled Acoustic/Elastic Wave 2D - SEM) show'
  else if(coupled_acoustic_poro) then
    write(24,*) '(Coupled Acoustic/Poroelastic Wave 2D - SEM) show'
  else if(coupled_elastic_poro) then
    write(24,*) '(Coupled Elastic/Poroelastic Wave 2D - SEM) show'
  else if(any_acoustic) then
    write(24,*) '(Acoustic Wave 2D - Spectral Element Method) show'
  else if(any_poroelastic) then
    write(24,*) '(Poroelastic Wave 2D - Spectral Element Method) show'
  else
    write(24,*) '(Elastic Wave 2D - Spectral Element Method) show'
  endif

  write(24,*) 'grestore'

  write(24,*) '%'
  write(24,*) '1 1 scale'
  write(24,*) '%'

!
!---- print the spectral elements mesh in PostScript
!

  endif


  convert = PI / 180.d0

!
!----  draw the velocity model in background
!
  if(modelvect) then

  buffer_offset = 0
  RGB_offset = 0

  do ispec=1,nspec
    do i=1,NGLLX-subsamp_postscript,subsamp_postscript
          do j=1,NGLLX-subsamp_postscript,subsamp_postscript

  if((vpImax-vpImin)/vpImin > 0.02d0) then

  if(assign_external_model) then

    x1 = (vpext(i,j,ispec)-vpImin) / (vpImax-vpImin)

  else

    material = kmato(ispec)

    if(poroelastic(ispec)) then

      ! poroelastic material

! get elastic parameters of current spectral element
    phil = porosity(kmato(ispec))
    tortl = tortuosity(kmato(ispec))
!solid properties
    mul_s = poroelastcoef(2,1,kmato(ispec))
    kappal_s = poroelastcoef(3,1,kmato(ispec)) - FOUR_THIRDS*mul_s
    rhol_s = density(1,kmato(ispec))
!fluid properties
    kappal_f = poroelastcoef(1,2,kmato(ispec))
    rhol_f = density(2,kmato(ispec))
!frame properties
    mul_fr = poroelastcoef(2,3,kmato(ispec))
    kappal_fr = poroelastcoef(3,3,kmato(ispec)) - FOUR_THIRDS*mul_fr
    rhol_bar =  (1.d0 - phil)*rhol_s + phil*rhol_f
!Biot coefficients for the input phi
      D_biot = kappal_s*(1.d0 + phil*(kappal_s/kappal_f - 1.d0))
      H_biot = (kappal_s - kappal_fr)*(kappal_s - kappal_fr)/(D_biot - kappal_fr) + kappal_fr + FOUR_THIRDS*mul_fr
      C_biot = kappal_s*(kappal_s - kappal_fr)/(D_biot - kappal_fr)
      M_biot = kappal_s*kappal_s/(D_biot - kappal_fr)
! Approximated velocities (no viscous dissipation)
      afactor = rhol_bar - phil/tortl*rhol_f
      bfactor = H_biot + phil*rhol_bar/(tortl*rhol_f)*M_biot - 2.d0*phil/tortl*C_biot
      cfactor = phil/(tortl*rhol_f)*(H_biot*M_biot - C_biot*C_biot)
      cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
      cpIloc = sqrt(cpIsquare)

    else

      lambdaplus2mu  = poroelastcoef(3,1,material)
      denst = density(1,material)
      cpIloc = sqrt(lambdaplus2mu/denst)

    endif

    x1 = (cpIloc-vpImin)/(vpImax-vpImin)

  endif

  else
    x1 = 0.5d0
  endif

! rescale to avoid very dark gray levels
  x1 = x1*0.7 + 0.2
  if(x1 > 1.d0) x1=1.d0

! invert scale: white = vpImin, dark gray = vpImax
  x1 = 1.d0 - x1

  xw = coord(1,ibool(i,j,ispec))
  zw = coord(2,ibool(i,j,ispec))
  xw = (xw-xmin)*ratio_page + orig_x
  zw = (zw-zmin)*ratio_page + orig_z
  xw = xw * centim
  zw = zw * centim
  if ( myrank == 0 ) then
     write(24,500) xw,zw
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_velocity_model(1,buffer_offset) = xw
     coorg_send_ps_velocity_model(2,buffer_offset) = zw
  endif

  xw = coord(1,ibool(i+subsamp_postscript,j,ispec))
  zw = coord(2,ibool(i+subsamp_postscript,j,ispec))
  xw = (xw-xmin)*ratio_page + orig_x
  zw = (zw-zmin)*ratio_page + orig_z
  xw = xw * centim
  zw = zw * centim
  if ( myrank == 0 ) then
     write(24,499) xw,zw
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_velocity_model(1,buffer_offset) = xw
     coorg_send_ps_velocity_model(2,buffer_offset) = zw
  endif

  xw = coord(1,ibool(i+subsamp_postscript,j+subsamp_postscript,ispec))
  zw = coord(2,ibool(i+subsamp_postscript,j+subsamp_postscript,ispec))
  xw = (xw-xmin)*ratio_page + orig_x
  zw = (zw-zmin)*ratio_page + orig_z
  xw = xw * centim
  zw = zw * centim
  if ( myrank == 0 ) then
     write(24,499) xw,zw
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_velocity_model(1,buffer_offset) = xw
     coorg_send_ps_velocity_model(2,buffer_offset) = zw
  endif

  xw = coord(1,ibool(i,j+subsamp_postscript,ispec))
  zw = coord(2,ibool(i,j+subsamp_postscript,ispec))
  xw = (xw-xmin)*ratio_page + orig_x
  zw = (zw-zmin)*ratio_page + orig_z
  xw = xw * centim
  zw = zw * centim
  if ( myrank == 0 ) then
     write(24,499) xw,zw
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_velocity_model(1,buffer_offset) = xw
     coorg_send_ps_velocity_model(2,buffer_offset) = zw
  endif

! display P-velocity model using gray levels
  if ( myrank == 0 ) then
     write(24,604) x1
  else
     RGB_offset = RGB_offset + 1
     RGB_send_ps_velocity_model(1,RGB_offset) = x1
  endif

          enddo
    enddo
  enddo

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        call MPI_RECV (coorg_recv_ps_velocity_model(1,1), &
             2*nspec_recv*((NGLLX-subsamp_postscript)/subsamp_postscript)*((NGLLX-subsamp_postscript)/subsamp_postscript)*4, &
             MPI_DOUBLE_PRECISION, iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        call MPI_RECV (RGB_recv_ps_velocity_model(1,1), nspec_recv*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
             ((NGLLX-subsamp_postscript)/subsamp_postscript), &
             MPI_DOUBLE_PRECISION, iproc, 42, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        buffer_offset = 0
        RGB_offset = 0
        do ispec = 1, nspec_recv
           do i=1,NGLLX-subsamp_postscript,subsamp_postscript
              do j=1,NGLLX-subsamp_postscript,subsamp_postscript
                 buffer_offset = buffer_offset + 1
                 write(24,500) coorg_recv_ps_velocity_model(1,buffer_offset), &
                               coorg_recv_ps_velocity_model(2,buffer_offset)
                 buffer_offset = buffer_offset + 1
                 write(24,499) coorg_recv_ps_velocity_model(1,buffer_offset), &
                               coorg_recv_ps_velocity_model(2,buffer_offset)
                 buffer_offset = buffer_offset + 1
                 write(24,499) coorg_recv_ps_velocity_model(1,buffer_offset), &
                               coorg_recv_ps_velocity_model(2,buffer_offset)
                 buffer_offset = buffer_offset + 1
                 write(24,499) coorg_recv_ps_velocity_model(1,buffer_offset), &
                               coorg_recv_ps_velocity_model(2,buffer_offset)
                 RGB_offset = RGB_offset + 1
                 write(24,604) RGB_recv_ps_velocity_model(1,RGB_offset)
              enddo
           enddo
        enddo

     enddo
  else
     call MPI_SEND (nspec, 1, MPI_INTEGER, 0, 42, MPI_COMM_WORLD, ier)
     call MPI_SEND (coorg_send_ps_velocity_model(1,1), 2*nspec*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
          ((NGLLX-subsamp_postscript)/subsamp_postscript)*4, &
          MPI_DOUBLE_PRECISION, 0, 42, MPI_COMM_WORLD, ier)
     call MPI_SEND (RGB_send_ps_velocity_model(1,1), nspec*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
          ((NGLLX-subsamp_postscript)/subsamp_postscript), &
          MPI_DOUBLE_PRECISION, 0, 42, MPI_COMM_WORLD, ier)
  endif


#endif


  endif

!
!---- draw the spectral element mesh
!

  if ( myrank == 0 ) then
     write(24,*) '%'
     write(24,*) '% spectral element mesh'
     write(24,*) '%'
  endif

  buffer_offset = 0
  RGB_offset = 0

  do ispec=1,nspec

  if ( myrank == 0 ) write(24,*) '% elem ',ispec

  do i=1,pointsdisp
  do j=1,pointsdisp
  xinterp(i,j) = 0.d0
  zinterp(i,j) = 0.d0
  do in = 1,ngnod
    nnum = knods(in,ispec)
      xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
      zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
  enddo
  enddo
  enddo

  is = 1
  ir = 1
  x1 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z1 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  if ( myrank == 0 ) then
     write(24,*) 'mark'
     write(24,681) x1,z1
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_element_mesh(1,buffer_offset) = x1
     coorg_send_ps_element_mesh(2,buffer_offset) = z1
  endif

  if(ngnod == 4) then

! draw straight lines if elements have 4 nodes

  ir=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_element_mesh(1,buffer_offset) = x2
     coorg_send_ps_element_mesh(2,buffer_offset) = z2
  endif

  ir=pointsdisp
  is=pointsdisp
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_element_mesh(1,buffer_offset) = x2
     coorg_send_ps_element_mesh(2,buffer_offset) = z2
  endif

  is=pointsdisp
  ir=1
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_element_mesh(1,buffer_offset) = x2
     coorg_send_ps_element_mesh(2,buffer_offset) = z2
  endif

  ir=1
  is=2
  x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
  z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,681) x2,z2
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_element_mesh(1,buffer_offset) = x2
     coorg_send_ps_element_mesh(2,buffer_offset) = z2
  endif

  else

! draw curved lines if elements have 9 nodes
  do ir=2,pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    if ( myrank == 0 ) then
       write(24,681) x2,z2
    else
       buffer_offset = buffer_offset + 1
       coorg_send_ps_element_mesh(1,buffer_offset) = x2
       coorg_send_ps_element_mesh(2,buffer_offset) = z2
    endif
  enddo

  ir=pointsdisp
  do is=2,pointsdisp
    x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    if ( myrank == 0 ) then
       write(24,681) x2,z2
    else
       buffer_offset = buffer_offset + 1
       coorg_send_ps_element_mesh(1,buffer_offset) = x2
       coorg_send_ps_element_mesh(2,buffer_offset) = z2
    endif
  enddo

  is=pointsdisp
  do ir=pointsdisp-1,1,-1
    x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    if ( myrank == 0 ) then
       write(24,681) x2,z2
    else
       buffer_offset = buffer_offset + 1
       coorg_send_ps_element_mesh(1,buffer_offset) = x2
       coorg_send_ps_element_mesh(2,buffer_offset) = z2
    endif
  enddo

  ir=1
  do is=pointsdisp-1,2,-1
    x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
    z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
    x2 = x2 * centim
    z2 = z2 * centim
    if ( myrank == 0 ) then
       write(24,681) x2,z2
    else
       buffer_offset = buffer_offset + 1
       coorg_send_ps_element_mesh(1,buffer_offset) = x2
       coorg_send_ps_element_mesh(2,buffer_offset) = z2
    endif
  enddo

  endif

  if ( myrank == 0 ) then
     write(24,*) 'CO'
  endif

  if(colors == 1) then

! use a different color for each material set
  imat = kmato(ispec)
  icol = mod(imat - 1,NUM_COLORS) + 1

! display all the PML layers in a different (constant) color if needed
  if(DISPLAY_PML_IN_DIFFERENT_COLOR .and. is_PML(ispec)) then
    icol = ICOLOR_FOR_PML_DISPLAY
    ! make sure that number exists
    if(icol > NUM_COLORS) icol = NUM_COLORS
  endif

  if (  myrank == 0 ) then
    if(meshvect) then
      write(24,680) red(icol),green(icol),blue(icol)
    else
      write(24,679) red(icol),green(icol),blue(icol)
    endif
  else
     RGB_offset = RGB_offset + 1
     color_send_ps_element_mesh(RGB_offset) = icol
  endif

  endif

  if ( myrank == 0 ) then
  if(meshvect) then
    if(modelvect) then
      write(24,*) 'Colmesh ST'
    else
      write(24,*) '0 setgray ST'
    endif
  endif
  endif

! write the element number, the group number and the material number inside the element
  if(numbers == 1) then

  xw = (coorg(1,knods(1,ispec)) + coorg(1,knods(2,ispec)) + coorg(1,knods(3,ispec)) + coorg(1,knods(4,ispec))) / 4.d0
  zw = (coorg(2,knods(1,ispec)) + coorg(2,knods(2,ispec)) + coorg(2,knods(3,ispec)) + coorg(2,knods(4,ispec))) / 4.d0
  xw = (xw-xmin)*ratio_page + orig_x
  zw = (zw-zmin)*ratio_page + orig_z
  xw = xw * centim
  zw = zw * centim

  if ( myrank == 0 ) then
  if(colors == 1) write(24,*) '1 setgray'
  endif

  if ( myrank == 0 ) then
     write(24,500) xw,zw
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_element_mesh(1,buffer_offset) = x2
     coorg_send_ps_element_mesh(2,buffer_offset) = z2
  endif

! write spectral element number
  if ( myrank == 0 ) then
     write(24,502) ispec
  else
     RGB_offset = RGB_offset + 1
     color_send_ps_element_mesh(RGB_offset) = ispec
  endif

  endif

  enddo

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, iproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        nb_coorg_per_elem = 1
        if ( numbers == 1 ) then
           nb_coorg_per_elem = nb_coorg_per_elem + 1
        endif
        if ( ngnod == 4 ) then
           nb_coorg_per_elem = nb_coorg_per_elem + 4
        else
           nb_coorg_per_elem = nb_coorg_per_elem + 3*(pointsdisp-1)+(pointsdisp-2)
        endif
        nb_color_per_elem = 0
        if ( colors == 1 ) then
           nb_color_per_elem = nb_color_per_elem + 1
        endif
        if ( numbers == 1 ) then
           nb_color_per_elem = nb_color_per_elem + 1
        endif

        call MPI_RECV (coorg_recv_ps_element_mesh(1,1), 2*nspec_recv*nb_coorg_per_elem, &
             MPI_DOUBLE_PRECISION, iproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        call MPI_RECV (color_recv_ps_element_mesh(1), nspec_recv*nb_coorg_per_elem, &
             MPI_INTEGER, iproc, 43, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        buffer_offset = 0
        RGB_offset = 0
        num_spec = nspec
        do ispec = 1, nspec_recv
           num_spec = num_spec + 1
           write(24,*) '% elem ',num_spec
           buffer_offset = buffer_offset + 1
           write(24,*) 'mark'
           write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
           if ( ngnod == 4 ) then
              buffer_offset = buffer_offset + 1
              write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
              buffer_offset = buffer_offset + 1
              write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
              buffer_offset = buffer_offset + 1
              write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
              buffer_offset = buffer_offset + 1
              write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)

           else
              do ir=2,pointsdisp
                 buffer_offset = buffer_offset + 1
                 write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
              enddo
              do is=2,pointsdisp
                 buffer_offset = buffer_offset + 1
                 write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
              enddo
              do ir=pointsdisp-1,1,-1
                 buffer_offset = buffer_offset + 1
                 write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
              enddo
              do is=pointsdisp-1,2,-1
                 buffer_offset = buffer_offset + 1
                 write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
              enddo

           endif

           write(24,*) 'CO'
           if ( colors == 1 ) then
              if(meshvect) then
                 RGB_offset = RGB_offset + 1
                 write(24,680) red(color_recv_ps_element_mesh(RGB_offset)),&
                               green(color_recv_ps_element_mesh(RGB_offset)),&
                               blue(color_recv_ps_element_mesh(RGB_offset))
              else
                 RGB_offset = RGB_offset + 1
                 write(24,679) red(color_recv_ps_element_mesh(RGB_offset)),&
                               green(color_recv_ps_element_mesh(RGB_offset)),&
                               blue(color_recv_ps_element_mesh(RGB_offset))
              endif
           endif
           if(meshvect) then
              if(modelvect) then
                 write(24,*) 'Colmesh ST'
              else
                 write(24,*) '0 setgray ST'
              endif
           endif
           if(numbers == 1) then
              if(colors == 1) write(24,*) '1 setgray'
              buffer_offset = buffer_offset + 1
              write(24,500) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
              RGB_offset = RGB_offset + 1
              write(24,502) color_recv_ps_element_mesh(RGB_offset)
           endif

        enddo

     enddo
  else
     call MPI_SEND (nspec, 1, MPI_INTEGER, 0, 43, MPI_COMM_WORLD, ier)
     nb_coorg_per_elem = 1
     if ( numbers == 1 ) then
        nb_coorg_per_elem = nb_coorg_per_elem + 1
     endif
     if ( ngnod == 4 ) then
        nb_coorg_per_elem = nb_coorg_per_elem + 4
     else
        nb_coorg_per_elem = nb_coorg_per_elem + 3*(pointsdisp-1)+(pointsdisp-2)
     endif
     nb_color_per_elem = 0
     if ( colors == 1 ) then
        nb_color_per_elem = nb_color_per_elem + 1
     endif
     if ( numbers == 1 ) then
        nb_color_per_elem = nb_color_per_elem + 1
     endif
     call MPI_SEND (coorg_send_ps_element_mesh(1,1), 2*nspec*nb_coorg_per_elem, &
          MPI_DOUBLE_PRECISION, 0, 43, MPI_COMM_WORLD, ier)
     if ( nb_color_per_elem > 0 ) then
        call MPI_SEND (color_send_ps_element_mesh(1), nspec*nb_color_per_elem, &
             MPI_INTEGER, 0, 43, MPI_COMM_WORLD, ier)
     endif

  endif

#endif

!
!--- draw absorbing boundaries with a thick color line
!
  anyabs_glob = anyabs
#ifdef USE_MPI
  call MPI_ALLREDUCE(anyabs, anyabs_glob, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  if(anyabs_glob .and. boundvect) then
  if ( myrank == 0 ) then
  write(24,*) '%'
  write(24,*) '% boundary conditions on the mesh'
  write(24,*) '%'

! use green color
  write(24,*) '0 1 0 RG'

  write(24,*) '0.02 CM setlinewidth'
  endif

  buffer_offset = 0

  if ( anyabs ) then
  do inum = 1,nelemabs
  ispec = numabs(inum)

  do iedge = 1,4

  if(codeabs(iedge,inum)) then ! codeabs(:,:) is defined as "logical" in MAIN program

  if(iedge == IEDGE1) then
    ideb = 1
    ifin = 2
  else if(iedge == IEDGE2) then
    ideb = 2
    ifin = 3
  else if(iedge == IEDGE3) then
    ideb = 3
    ifin = 4
  else if(iedge == IEDGE4) then
    ideb = 4
    ifin = 1
  else
    call exit_MPI('Wrong codeabs() absorbing boundary code')
  endif

! draw the Stacey absorbing boundary line segment in different colors depending on its type
  if ( myrank == 0 ) then
  if(typeabs(inum) == IBOTTOM) then
    write(24,*) '0 1 0 RG'  ! Green
  else if(typeabs(inum) == IRIGHT) then
    write(24,*) '0 0 1 RG'  ! Blue
  else if(typeabs(inum) == ITOP) then
    write(24,*) '1 0.7529 0.7960 RG' ! Pink
  else if(typeabs(inum) == ILEFT) then
    write(24,*) '1 0.6470 0 RG' ! Orange
  else
    call exit_MPI('Wrong typeabs() absorbing boundary code')
  endif
  endif

  x1 = (coorg(1,knods(ideb,ispec))-xmin)*ratio_page + orig_x
  z1 = (coorg(2,knods(ideb,ispec))-zmin)*ratio_page + orig_z
  x2 = (coorg(1,knods(ifin,ispec))-xmin)*ratio_page + orig_x
  z2 = (coorg(2,knods(ifin,ispec))-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,602) x1,z1,x2,z2
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_abs(1,buffer_offset) = x1
     coorg_send_ps_abs(2,buffer_offset) = z1
     coorg_send_ps_abs(3,buffer_offset) = x2
     coorg_send_ps_abs(4,buffer_offset) = z2
  endif

  endif ! of if(codeabs(iedge,inum))
  enddo ! of do iedge = 1,4

  enddo
  endif

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, iproc, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        if ( nspec_recv > 0 ) then
        call MPI_RECV (coorg_recv_ps_abs(1,1), 4*nspec_recv, &
             MPI_DOUBLE_PRECISION, iproc, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        buffer_offset = 0
        do ispec = 1, nspec_recv
           buffer_offset = buffer_offset + 1
           write(24,602) coorg_recv_ps_abs(1,buffer_offset), coorg_recv_ps_abs(2,buffer_offset), &
                coorg_recv_ps_abs(3,buffer_offset), coorg_recv_ps_abs(4,buffer_offset)
        enddo
        endif
     enddo
  else
     call MPI_SEND (buffer_offset, 1, MPI_INTEGER, 0, 44, MPI_COMM_WORLD, ier)
     if ( buffer_offset > 0 ) then
     call MPI_SEND (coorg_send_ps_abs(1,1), 4*buffer_offset, &
          MPI_DOUBLE_PRECISION, 0, 44, MPI_COMM_WORLD, ier)
     endif

  endif

#endif

  if ( myrank == 0 ) then
    write(24,*) '0 setgray'
    write(24,*) '0 setlinewidth'
  endif

  endif

!
!--- draw free surface with a thick color line
!

  if ( myrank == 0 ) then
  write(24,*) '%'
  write(24,*) '% free surface on the mesh'
  write(24,*) '%'

! use orange color
  write(24,*) '1 0.66 0 RG'

  write(24,*) '0.02 CM setlinewidth'
  endif

  buffer_offset = 0

  if ( nelem_acoustic_surface > 0 ) then
  do inum = 1,nelem_acoustic_surface
  ispec = acoustic_edges(1,inum)

  x1 = (coorg(1,acoustic_edges(3,inum))-xmin)*ratio_page + orig_x
  z1 = (coorg(2,acoustic_edges(3,inum))-zmin)*ratio_page + orig_z
  x2 = (coorg(1,acoustic_edges(4,inum))-xmin)*ratio_page + orig_x
  z2 = (coorg(2,acoustic_edges(4,inum))-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,602) x1,z1,x2,z2
  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_free_surface(1,buffer_offset) = x1
     coorg_send_ps_free_surface(2,buffer_offset) = z1
     coorg_send_ps_free_surface(3,buffer_offset) = x2
     coorg_send_ps_free_surface(4,buffer_offset) = z2
  endif

  enddo
  endif

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, iproc, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        if ( nspec_recv > 0 ) then
        call MPI_RECV (coorg_recv_ps_free_surface(1,1), 4*nspec_recv, &
             MPI_DOUBLE_PRECISION, iproc, 44, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        buffer_offset = 0
        do ispec = 1, nspec_recv
           buffer_offset = buffer_offset + 1
           write(24,602) coorg_recv_ps_free_surface(1,buffer_offset), coorg_recv_ps_free_surface(2,buffer_offset), &
                coorg_recv_ps_free_surface(3,buffer_offset), coorg_recv_ps_free_surface(4,buffer_offset)
        enddo
        endif
     enddo
  else
     call MPI_SEND (buffer_offset, 1, MPI_INTEGER, 0, 44, MPI_COMM_WORLD, ier)
     if ( buffer_offset > 0 ) then
     call MPI_SEND (coorg_send_ps_free_surface(1,1), 4*buffer_offset, &
          MPI_DOUBLE_PRECISION, 0, 44, MPI_COMM_WORLD, ier)
     endif

  endif

#endif

  if ( myrank == 0 ) then
    write(24,*) '0 setgray'
    write(24,*) '0 setlinewidth'
  endif

!
!----  draw the fluid-solid coupling edges with a thick color line
!
  coupled_acoustic_elastic_glob = coupled_acoustic_elastic
#ifdef USE_MPI
  call MPI_ALLREDUCE(coupled_acoustic_elastic, coupled_acoustic_elastic_glob, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  if(coupled_acoustic_elastic_glob .and. boundvect) then

  if ( myrank == 0 ) then
  write(24,*) '%'
  write(24,*) '% fluid-solid coupling edges in the mesh'
  write(24,*) '%'

  write(24,*) '0.02 CM setlinewidth'
  endif

  if ( myrank /= 0 .and. num_fluid_solid_edges > 0 ) allocate(coorg_send(4,num_fluid_solid_edges))
  buffer_offset = 0

! loop on all the coupling edges
  do inum = 1,num_fluid_solid_edges

! get the edge of the acoustic element
   ispec = fluid_solid_acoustic_ispec(inum)
   iedge = fluid_solid_acoustic_iedge(inum)

! use pink color
  if ( myrank == 0 ) write(24,*) '1 0.75 0.8 RG'

  if(iedge == ITOP) then
    ideb = 3
    ifin = 4
  else if(iedge == IBOTTOM) then
    ideb = 1
    ifin = 2
  else if(iedge == ILEFT) then
    ideb = 4
    ifin = 1
  else if(iedge == IRIGHT) then
    ideb = 2
    ifin = 3
  else
    call exit_MPI('Wrong fluid-solid coupling edge code')
  endif

  x1 = (coorg(1,knods(ideb,ispec))-xmin)*ratio_page + orig_x
  z1 = (coorg(2,knods(ideb,ispec))-zmin)*ratio_page + orig_z
  x2 = (coorg(1,knods(ifin,ispec))-xmin)*ratio_page + orig_x
  z2 = (coorg(2,knods(ifin,ispec))-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,602) x1,z1,x2,z2
  else
     buffer_offset = buffer_offset + 1
     coorg_send(1,buffer_offset) = x1
     coorg_send(2,buffer_offset) = z1
     coorg_send(3,buffer_offset) = x2
     coorg_send(4,buffer_offset) = z2
  endif

  enddo

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, iproc, 45, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        if ( nspec_recv > 0 ) then
        allocate(coorg_recv(4,nspec_recv))
        call MPI_RECV (coorg_recv(1,1), 4*nspec_recv, &
             MPI_DOUBLE_PRECISION, iproc, 45, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        buffer_offset = 0
        do ispec = 1, nspec_recv
           buffer_offset = buffer_offset + 1
           write(24,*) '1 0.75 0.8 RG'
           write(24,602) coorg_recv(1,buffer_offset), coorg_recv(2,buffer_offset), &
                coorg_recv(3,buffer_offset), coorg_recv(4,buffer_offset)
        enddo
        deallocate(coorg_recv)
        endif
     enddo
  else
     call MPI_SEND (buffer_offset, 1, MPI_INTEGER, 0, 45, MPI_COMM_WORLD, ier)
     if ( buffer_offset > 0 ) then
     call MPI_SEND (coorg_send(1,1), 4*buffer_offset, &
          MPI_DOUBLE_PRECISION, 0, 45, MPI_COMM_WORLD, ier)
     deallocate(coorg_send)
     endif
  endif

#endif

  if ( myrank == 0 ) then
    write(24,*) '0 setgray'
    write(24,*) '0 setlinewidth'
  endif

  endif

!
!----  draw the fluid-porous coupling edges with a thick color line
!
  coupled_acoustic_poro_glob = coupled_acoustic_poro
#ifdef USE_MPI
  call MPI_ALLREDUCE(coupled_acoustic_poro, coupled_acoustic_poro_glob, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  if(coupled_acoustic_poro_glob .and. boundvect) then

  if ( myrank == 0 ) then
  write(24,*) '%'
  write(24,*) '% fluid-porous coupling edges in the mesh'
  write(24,*) '%'

  write(24,*) '0.02 CM setlinewidth'
  endif

  if ( myrank /= 0 .and. num_fluid_poro_edges > 0 ) allocate(coorg_send(4,num_fluid_poro_edges))
  buffer_offset = 0

! loop on all the coupling edges
  do inum = 1,num_fluid_poro_edges

! get the edge of the acoustic element
   ispec = fluid_poro_acoustic_ispec(inum)
   iedge = fluid_poro_acoustic_iedge(inum)

! use pink color
  if ( myrank == 0 ) write(24,*) '1 0.75 0.8 RG'

  if(iedge == ITOP) then
    ideb = 3
    ifin = 4
  else if(iedge == IBOTTOM) then
    ideb = 1
    ifin = 2
  else if(iedge == ILEFT) then
    ideb = 4
    ifin = 1
  else if(iedge == IRIGHT) then
    ideb = 2
    ifin = 3
  else
    call exit_MPI('Wrong fluid-solid coupling edge code')
  endif

  x1 = (coorg(1,knods(ideb,ispec))-xmin)*ratio_page + orig_x
  z1 = (coorg(2,knods(ideb,ispec))-zmin)*ratio_page + orig_z
  x2 = (coorg(1,knods(ifin,ispec))-xmin)*ratio_page + orig_x
  z2 = (coorg(2,knods(ifin,ispec))-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,602) x1,z1,x2,z2
  else
     buffer_offset = buffer_offset + 1
     coorg_send(1,buffer_offset) = x1
     coorg_send(2,buffer_offset) = z1
     coorg_send(3,buffer_offset) = x2
     coorg_send(4,buffer_offset) = z2
  endif

  enddo

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, iproc, 45, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        if ( nspec_recv > 0 ) then
        allocate(coorg_recv(4,nspec_recv))
        call MPI_RECV (coorg_recv(1,1), 4*nspec_recv, &
             MPI_DOUBLE_PRECISION, iproc, 45, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        buffer_offset = 0
        do ispec = 1, nspec_recv
           buffer_offset = buffer_offset + 1
           write(24,*) '1 0.75 0.8 RG'
           write(24,602) coorg_recv(1,buffer_offset), coorg_recv(2,buffer_offset), &
                coorg_recv(3,buffer_offset), coorg_recv(4,buffer_offset)
        enddo
        deallocate(coorg_recv)
        endif
     enddo
  else
     call MPI_SEND (buffer_offset, 1, MPI_INTEGER, 0, 45, MPI_COMM_WORLD, ier)
     if ( buffer_offset > 0 ) then
     call MPI_SEND (coorg_send(1,1), 4*buffer_offset, &
          MPI_DOUBLE_PRECISION, 0, 45, MPI_COMM_WORLD, ier)
     deallocate(coorg_send)
     endif
  endif

#endif

  if ( myrank == 0 ) then
    write(24,*) '0 setgray'
    write(24,*) '0 setlinewidth'
  endif

  endif

!
!----  draw the solid-porous coupling edges with a thick color line
!
  coupled_elastic_poro_glob = coupled_elastic_poro
#ifdef USE_MPI
  call MPI_ALLREDUCE(coupled_elastic_poro, coupled_elastic_poro_glob, 1, MPI_LOGICAL, MPI_LOR, MPI_COMM_WORLD, ier)
#endif

  if(coupled_elastic_poro_glob .and. boundvect) then

  if ( myrank == 0 ) then
  write(24,*) '%'
  write(24,*) '% solid-porous coupling edges in the mesh'
  write(24,*) '%'

  write(24,*) '0.02 CM setlinewidth'
  endif

  if ( myrank /= 0 .and. num_solid_poro_edges > 0 ) allocate(coorg_send(4,num_solid_poro_edges))
  buffer_offset = 0

! loop on all the coupling edges
  do inum = 1,num_solid_poro_edges

! get the edge of the poroelastic element
   ispec = solid_poro_poroelastic_ispec(inum)
   iedge = solid_poro_poroelastic_iedge(inum)

! use pink color
  if ( myrank == 0 ) write(24,*) '1 0.75 0.8 RG'

  if(iedge == ITOP) then
    ideb = 3
    ifin = 4
  else if(iedge == IBOTTOM) then
    ideb = 1
    ifin = 2
  else if(iedge == ILEFT) then
    ideb = 4
    ifin = 1
  else if(iedge == IRIGHT) then
    ideb = 2
    ifin = 3
  else
    call exit_MPI('Wrong fluid-solid coupling edge code')
  endif

  x1 = (coorg(1,knods(ideb,ispec))-xmin)*ratio_page + orig_x
  z1 = (coorg(2,knods(ideb,ispec))-zmin)*ratio_page + orig_z
  x2 = (coorg(1,knods(ifin,ispec))-xmin)*ratio_page + orig_x
  z2 = (coorg(2,knods(ifin,ispec))-zmin)*ratio_page + orig_z
  x1 = x1 * centim
  z1 = z1 * centim
  x2 = x2 * centim
  z2 = z2 * centim
  if ( myrank == 0 ) then
     write(24,602) x1,z1,x2,z2
  else
     buffer_offset = buffer_offset + 1
     coorg_send(1,buffer_offset) = x1
     coorg_send(2,buffer_offset) = z1
     coorg_send(3,buffer_offset) = x2
     coorg_send(4,buffer_offset) = z2
  endif

  enddo

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, iproc, 45, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        if ( nspec_recv > 0 ) then
        allocate(coorg_recv(4,nspec_recv))
        call MPI_RECV (coorg_recv(1,1), 4*nspec_recv, &
             MPI_DOUBLE_PRECISION, iproc, 45, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        buffer_offset = 0
        do ispec = 1, nspec_recv
           buffer_offset = buffer_offset + 1
           write(24,*) '1 0.75 0.8 RG'
           write(24,602) coorg_recv(1,buffer_offset), coorg_recv(2,buffer_offset), &
                coorg_recv(3,buffer_offset), coorg_recv(4,buffer_offset)
        enddo
        deallocate(coorg_recv)
        endif
     enddo
  else
     call MPI_SEND (buffer_offset, 1, MPI_INTEGER, 0, 45, MPI_COMM_WORLD, ier)
     if ( buffer_offset > 0 ) then
     call MPI_SEND (coorg_send(1,1), 4*buffer_offset, &
          MPI_DOUBLE_PRECISION, 0, 45, MPI_COMM_WORLD, ier)
     deallocate(coorg_send)
     endif
  endif

#endif

  if ( myrank == 0 ) then
    write(24,*) '0 setgray'
    write(24,*) '0 setlinewidth'
  endif

  endif

!
!----  draw the normalized vector field
!

  if ( myrank == 0 ) then
! return if the maximum vector equals zero (no source)
  if(dispmax == 0.d0) then
    write(IOUT,*) 'null vector: returning!'
    return
  endif

  write(24,*) '%'
  write(24,*) '% vector field'
  write(24,*) '%'

! color arrows if we draw the velocity model in the background
  if(modelvect) then
        write(24,*) 'Colvects'
  else
        write(24,*) '0 setgray'
  endif
  endif

  if(interpol) then

  if (myrank == 0) write(IOUT,*) 'Interpolating the vector field...'

! option to plot only lowerleft corner value to avoid very large files if dense meshes
  if(plot_lowerleft_corner_only) then
    pointsdisp_loop = 1
  else
    pointsdisp_loop = pointsdisp
  endif

  buffer_offset = 0

  do ispec=1,nspec

! interpolation on a uniform grid
#ifdef USE_MPI
  if(myrank == 0 .and. mod(ispec,1000) == 0) write(IOUT,*) 'Interpolation uniform grid element ',ispec,' on processor core 0'
#else
  if(mod(ispec,1000) == 0) write(IOUT,*) 'Interpolation uniform grid element ',ispec
#endif

  do i=1,pointsdisp_loop
  do j=1,pointsdisp_loop

  xinterp(i,j) = 0.d0
  zinterp(i,j) = 0.d0
  do in = 1,ngnod
    nnum = knods(in,ispec)
      xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
      zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
  enddo

  Uxinterp(i,j) = 0.d0
  Uzinterp(i,j) = 0.d0

  do k=1,NGLLX
  do l=1,NGLLX
    if(AXISYM) then
      if(is_on_the_axis(ispec)) then
        Uxinterp(i,j) = Uxinterp(i,j) + vector_field_display(1,ibool(k,l,ispec))*flagrange_GLJ(k,i)*flagrange_GLJ(l,j)
        Uzinterp(i,j) = Uzinterp(i,j) + vector_field_display(3,ibool(k,l,ispec))*flagrange_GLJ(k,i)*flagrange_GLJ(l,j)
      else
        Uxinterp(i,j) = Uxinterp(i,j) + vector_field_display(1,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
        Uzinterp(i,j) = Uzinterp(i,j) + vector_field_display(3,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
      endif
    else
      Uxinterp(i,j) = Uxinterp(i,j) + vector_field_display(1,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
      Uzinterp(i,j) = Uzinterp(i,j) + vector_field_display(3,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
    endif
  enddo
  enddo

  x1 =(xinterp(i,j)-xmin)*ratio_page
  z1 =(zinterp(i,j)-zmin)*ratio_page
  x2 = Uxinterp(i,j)*sizemax_arrows/dispmax
  z2 = Uzinterp(i,j)*sizemax_arrows/dispmax

  d = sqrt(x2**2 + z2**2)

! ignore if vector is too small
  if(d > cutsnaps*sizemax_arrows) then

  d1 = d * ARROW_RATIO
  d2 = d1 * cos(ARROW_ANGLE*convert)

  dummy = x2/d
  if(dummy > 0.9999d0) dummy = 0.9999d0
  if(dummy < -0.9999d0) dummy = -0.9999d0
  theta = acos(dummy)

  if(z2 < 0.d0) theta = 360.d0*convert - theta
  thetaup = theta - ARROW_ANGLE*convert
  thetadown = theta + ARROW_ANGLE*convert

! draw the vector
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  x2 = x2 * centim
  z2 = z2 * centim
  xa = -d2*cos(thetaup)
  za = -d2*sin(thetaup)
  xa = xa * centim
  za = za * centim
  xb = -d2*cos(thetadown)
  zb = -d2*sin(thetadown)
  xb = xb * centim
  zb = zb * centim
  if ( myrank == 0 ) then
  write(postscript_line,700) xb,zb,xa,za,x2,z2,x1,z1
! suppress useless white spaces to make PostScript file smaller
! suppress leading white spaces again, if any
  postscript_line = adjustl(postscript_line)

  line_length = len_trim(postscript_line)
  index_char = 1
  first = .false.
  do ii = 1,line_length-1
    if(ch1(ii) /= ' ' .or. first) then
      if(ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
        ch2(index_char) = ch1(ii)
        index_char = index_char + 1
        first = .true.
      endif
    endif
  enddo
  ch2(index_char) = ch1(line_length)
  write(24,"(100(a1))") (ch2(ii), ii=1,index_char)

  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_vector_field(1,buffer_offset) = xb
     coorg_send_ps_vector_field(2,buffer_offset) = zb
     coorg_send_ps_vector_field(3,buffer_offset) = xa
     coorg_send_ps_vector_field(4,buffer_offset) = za
     coorg_send_ps_vector_field(5,buffer_offset) = x2
     coorg_send_ps_vector_field(6,buffer_offset) = z2
     coorg_send_ps_vector_field(7,buffer_offset) = x1
     coorg_send_ps_vector_field(8,buffer_offset) = z1
  endif

  endif

  enddo
  enddo
  enddo

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, iproc, 46, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        if ( nspec_recv > 0 ) then
        call MPI_RECV (coorg_recv_ps_vector_field(1,1), 8*nspec_recv, &
             MPI_DOUBLE_PRECISION, iproc, 46, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        buffer_offset = 0
        do ispec = 1, nspec_recv
           buffer_offset = buffer_offset + 1
             write(postscript_line,700) coorg_recv_ps_vector_field(1,buffer_offset), &
                  coorg_recv_ps_vector_field(2,buffer_offset), &
                  coorg_recv_ps_vector_field(3,buffer_offset), coorg_recv_ps_vector_field(4,buffer_offset), &
                  coorg_recv_ps_vector_field(5,buffer_offset), coorg_recv_ps_vector_field(6,buffer_offset), &
                  coorg_recv_ps_vector_field(7,buffer_offset), coorg_recv_ps_vector_field(8,buffer_offset)

             ! suppress useless white spaces to make PostScript file smaller
             ! suppress leading white spaces again, if any
             postscript_line = adjustl(postscript_line)

             line_length = len_trim(postscript_line)
             index_char = 1
             first = .false.
             do ii = 1,line_length-1
                if(ch1(ii) /= ' ' .or. first) then
                   if(ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
                      ch2(index_char) = ch1(ii)
                      index_char = index_char + 1
                      first = .true.
                   endif
                endif
             enddo
             ch2(index_char) = ch1(line_length)
             write(24,"(100(a1))") (ch2(ii), ii=1,index_char)
          enddo
          endif
       enddo
    else
       call MPI_SEND (buffer_offset, 1, MPI_INTEGER, 0, 46, MPI_COMM_WORLD, ier)
       if ( buffer_offset > 0 ) then
       call MPI_SEND (coorg_send_ps_vector_field(1,1), 8*buffer_offset, &
            MPI_DOUBLE_PRECISION, 0, 46, MPI_COMM_WORLD, ier)
       endif

  endif

#endif


! draw the vectors at the nodes of the mesh if we do not interpolate the display on a regular grid
  else

  buffer_offset = 0

  do ipoin=1,nglob

  x1 =(coord(1,ipoin)-xmin)*ratio_page
  z1 =(coord(2,ipoin)-zmin)*ratio_page

  x2 = vector_field_display(1,ipoin)*sizemax_arrows/dispmax
  z2 = vector_field_display(3,ipoin)*sizemax_arrows/dispmax

  d = sqrt(x2**2 + z2**2)

! ignore if vector is too small
  if(d > cutsnaps*sizemax_arrows) then

  d1 = d * ARROW_RATIO
  d2 = d1 * cos(ARROW_ANGLE*convert)

  dummy = x2/d
  if(dummy > 0.9999d0) dummy = 0.9999d0
  if(dummy < -0.9999d0) dummy = -0.9999d0
  theta = acos(dummy)

  if(z2 < 0.d0) theta = 360.d0*convert - theta
  thetaup = theta - ARROW_ANGLE*convert
  thetadown = theta + ARROW_ANGLE*convert

! draw the vector
  x1 = (orig_x+x1) * centim
  z1 = (orig_z+z1) * centim
  x2 = x2 * centim
  z2 = z2 * centim
  xa = -d2*cos(thetaup)
  za = -d2*sin(thetaup)
  xa = xa * centim
  za = za * centim
  xb = -d2*cos(thetadown)
  zb = -d2*sin(thetadown)
  xb = xb * centim
  zb = zb * centim
  if ( myrank == 0 ) then
  write(postscript_line,700) xb,zb,xa,za,x2,z2,x1,z1

! suppress useless white spaces to make PostScript file smaller
! suppress leading white spaces again, if any
  postscript_line = adjustl(postscript_line)

  line_length = len_trim(postscript_line)
  index_char = 1
  first = .false.
  do ii = 1,line_length-1
    if(ch1(ii) /= ' ' .or. first) then
      if(ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
        ch2(index_char) = ch1(ii)
        index_char = index_char + 1
        first = .true.
      endif
    endif
  enddo
  ch2(index_char) = ch1(line_length)
  write(24,"(100(a1))") (ch2(ii), ii=1,index_char)

  else
     buffer_offset = buffer_offset + 1
     coorg_send_ps_vector_field(1,buffer_offset) = xb
     coorg_send_ps_vector_field(2,buffer_offset) = zb
     coorg_send_ps_vector_field(3,buffer_offset) = xa
     coorg_send_ps_vector_field(4,buffer_offset) = za
     coorg_send_ps_vector_field(5,buffer_offset) = x2
     coorg_send_ps_vector_field(6,buffer_offset) = z2
     coorg_send_ps_vector_field(7,buffer_offset) = x1
     coorg_send_ps_vector_field(8,buffer_offset) = z1
  endif
  endif

  enddo

#ifdef USE_MPI
  if (myrank == 0 ) then

     do iproc = 1, nproc-1
        call MPI_RECV (nspec_recv, 1, MPI_INTEGER, iproc, 47, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)
        if ( nspec_recv > 0 ) then
        call MPI_RECV (coorg_recv_ps_vector_field(1,1), 8*nspec_recv, &
             MPI_DOUBLE_PRECISION, iproc, 47, MPI_COMM_WORLD, MPI_STATUS_IGNORE, ier)

        buffer_offset = 0
        do ispec = 1, nspec_recv
           buffer_offset = buffer_offset + 1
             write(postscript_line,700) coorg_recv_ps_vector_field(1,buffer_offset), &
                  coorg_recv_ps_vector_field(2,buffer_offset), &
                  coorg_recv_ps_vector_field(3,buffer_offset), coorg_recv_ps_vector_field(4,buffer_offset), &
                  coorg_recv_ps_vector_field(5,buffer_offset), coorg_recv_ps_vector_field(6,buffer_offset), &
                  coorg_recv_ps_vector_field(7,buffer_offset), coorg_recv_ps_vector_field(8,buffer_offset)

             ! suppress useless white spaces to make PostScript file smaller
             ! suppress leading white spaces again, if any
             postscript_line = adjustl(postscript_line)

             line_length = len_trim(postscript_line)
             index_char = 1
             first = .false.
             do ii = 1,line_length-1
                if(ch1(ii) /= ' ' .or. first) then
                   if(ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
                      ch2(index_char) = ch1(ii)
                      index_char = index_char + 1
                      first = .true.
                   endif
                endif
             enddo
             ch2(index_char) = ch1(line_length)
             write(24,"(100(a1))") (ch2(ii), ii=1,index_char)
          enddo
          endif
       enddo
    else
       call MPI_SEND (buffer_offset, 1, MPI_INTEGER, 0, 47, MPI_COMM_WORLD, ier)
       if ( buffer_offset > 0 ) then
       call MPI_SEND (coorg_send_ps_vector_field(1,1), 8*buffer_offset, &
            MPI_DOUBLE_PRECISION, 0, 47, MPI_COMM_WORLD, ier)
       endif
  endif

#endif

  endif

  if ( myrank == 0 ) then
  write(24,*) '0 setgray'

! sources and receivers in color if velocity model
  if(modelvect) then
    write(24,*) 'Colreceiv'
  else
    write(24,*) '0 setgray'
  endif

!
!----  write position of the source
!
  do i=1,NSOURCES
    if(i == 1) write(24,*) '% beginning of source line'
    if(i == NSOURCES) write(24,*) '% end of source line'
  xw = x_source(i)
  zw = z_source(i)
  xw = (xw-xmin)*ratio_page + orig_x
  zw = (zw-zmin)*ratio_page + orig_z
  xw = xw * centim
  zw = zw * centim
  write(24,500) xw,zw
  write(24,*) 'Cross'
  enddo

!
!----  write position of the receivers
!
  do i=1,nrec
    if(i == 1) write(24,*) '% beginning of receiver line'
    if(i == nrec) write(24,*) '% end of receiver line'

    xw = st_xval(i)
    zw = st_zval(i)

    xw = (xw-xmin)*ratio_page + orig_x
    zw = (zw-zmin)*ratio_page + orig_z
    xw = xw * centim
    zw = zw * centim
    write(24,500) xw,zw
    write(24,*) 'Diamond'
  enddo

  write(24,*) '%'
  write(24,*) 'grestore'
  write(24,*) 'showpage'

  close(24)
  endif

 10  format('%!PS-Adobe-2.0',/,'%%',/,'%% Title: ',a100,/,'%% Created by: Specfem2D',/,'%% Author: Dimitri Komatitsch',/,'%%')
 600 format(f6.3,' neg CM 0 MR (Time =',f8.3,' s) show')
 601 format(f6.3,' neg CM 0 MR (Time =',1pe12.3,' s) show')
 610 format(f6.3,' neg CM 0 MR (Time step = ',i7,') show')
 620 format(f6.3,' neg CM 0 MR (Cut =',f5.2,' \%) show')
 640 format(f6.3,' neg CM 0 MR (Max norm =',1pe12.3,') show')

 499 format(f8.3,1x,f8.3,' L')
 500 format(f8.3,1x,f8.3,' M')
 502 format('fN (',i7,') Cshow')
 679 format(f12.6,1x,f12.6,1x,f12.6,' RG fill stroke')
 680 format(f12.6,1x,f12.6,1x,f12.6,' RG GF')
 681 format(f6.2,1x,f6.2)
 602 format(f6.2,1x,f6.2,' M ',f6.2,1x,f6.2,' L ST')
 604 format('CP ',f12.6,' BK')
 700 format(8(f6.2,1x),'F')

  end subroutine plotpost


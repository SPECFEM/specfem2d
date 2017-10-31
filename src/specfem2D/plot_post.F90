!========================================================================
!
!                   S P E C F E M 2 D  Version 7 . 0
!                   --------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
! The full text of the license is available in file "LICENSE".
!
!========================================================================

  subroutine plot_post()

!
! PostScript display routine
!

#ifdef USE_MPI
  use mpi
#endif

  use constants, only: IMAIN,NGLLX,NGLLZ,ARROW_ANGLE,ARROW_RATIO,CENTIM, &
    DISPLAY_PML_IN_DIFFERENT_COLOR,ICOLOR_FOR_PML_DISPLAY, &
    IEDGE1,IEDGE2,IEDGE3,IEDGE4, &
    IRIGHT,ILEFT,IBOTTOM,ITOP, &
    ORIG_X,ORIG_Z,PI,RPERCENTX,RPERCENTZ,STABILITY_THRESHOLD, &
    DISPLAY_COLORS,DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT,OUTPUT_FILES

  use specfem_par, only: coord,vpext,x_source,z_source,st_xval,st_zval,it,deltat,coorg,density, &
                         AXISYM,is_on_the_axis,flagrange_GLJ, &
                         poroelastcoef,knods,kmato,ibool, &
                         numabs,codeabs,typeabs,anyabs,nelem_acoustic_surface, acoustic_edges, &
                         nglob,nrec,NSOURCES, &
                         assign_external_model,nelemabs,pointsdisp, &
                         nspec,ngnod,coupled_acoustic_elastic,coupled_acoustic_poro,coupled_elastic_poro, &
                         any_acoustic,any_poroelastic, &
                         fluid_solid_acoustic_ispec,fluid_solid_acoustic_iedge,num_fluid_solid_edges, &
                         fluid_poro_acoustic_ispec,fluid_poro_acoustic_iedge,num_fluid_poro_edges, &
                         solid_poro_poroelastic_ispec,solid_poro_poroelastic_iedge,num_solid_poro_edges, &
                         ispec_is_poroelastic,myrank,NPROC

  ! PML arrays
  use specfem_par, only: PML_BOUNDARY_CONDITIONS,ispec_is_PML

  ! movie images
  use specfem_par_movie, only: vector_field_display,simulation_title, &
    xinterp,zinterp,Uxinterp,Uzinterp,flagrange,shape2D_display, &
    subsamp_postscript,imagetype_postscript,interpol,meshvect,modelvect, &
    cutsnaps,sizemax_arrows,boundvect,plot_lowerleft_corner_only, &
    vpImin,vpImax, &
    coorg_send_ps_velocity_model,RGB_send_ps_velocity_model, &
    coorg_recv_ps_velocity_model,RGB_recv_ps_velocity_model, &
    coorg_send_ps_element_mesh,color_send_ps_element_mesh, &
    coorg_recv_ps_element_mesh,color_recv_ps_element_mesh, &
    coorg_send_ps_abs,coorg_recv_ps_abs, &
    coorg_send_ps_free_surface,coorg_recv_ps_free_surface, &
    coorg_send_ps_vector_field,coorg_recv_ps_vector_field,US_LETTER

  implicit none

  ! flag to display a deformed mesh instead of the displacement vector, can be useful for instance
  ! in non-destructive testing, for Lamb waves and so on.
  logical, parameter :: DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR = .false.

  ! to compute the deformed mesh we take the original mesh and add the displacement of the current point
  ! to the coordinates of that original (undeformed) mesh point, with this scaling factor that the user
  ! can make different from 1 if he/she wants in order to reduce or increase the visual effect of the deformation
  ! for display purposes
  double precision, parameter :: SCALING_FACTOR_FOR_DEFORMED_MESH = 1.d0  !  0.05d0 !  1000.d0

  ! suppress the display of the receivers with diamonds or not if drawing the deformed mesh
  ! (can be useful e.g. when showing Lamb waves i.e. when the mesh is very elongated and a lot of diamonds
  !  on the display would make it hard to see the deformed mesh behind)
  logical, parameter :: SUPPRESS_DISPLAY_OF_RECEIVERS_IF_DEFORMED_MESH = .true.

  ! color palette
  integer, parameter :: NUM_COLORS = 236
  double precision, dimension(NUM_COLORS) :: red,green,blue

  double precision, dimension(:,:), allocatable  :: coorg_send
  double precision, dimension(:,:), allocatable  :: coorg_recv
  integer :: k,j,ispec,material,is,ir,imat,icol,l,line_length
  integer :: index_char,ii,ipoin,in,nnum,inum,ideb,ifin,iedge
  integer :: ier

  ! for the file name
  character(len=100) :: file_name
  integer  :: buffer_offset, RGB_offset
  double precision convert,x1,cpIloc,xa,za,xb,zb,lambdaplus2mu,denst
  double precision z1,x2,z2,d,d1,d2,dummy,theta,thetaup,thetadown
  double precision :: cpIsquare

  double precision :: phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar
  double precision :: D_biot,H_biot,C_biot,M_biot

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

  double precision :: afactor,bfactor,cfactor
  double precision :: xmax,zmax,height,xw,zw,usoffset,sizex,sizez,timeval
  ! for MPI collection
  double precision :: xmin_glob, xmax_glob, zmin_glob, zmax_glob
  double precision :: dispmax_glob

#ifndef USE_MPI
! this to avoid warnings by the compiler about unused variables in the case
! of a serial code, therefore use them once and do nothing: just set them to zero
  nspec_recv = 0
  nb_coorg_per_elem = 0
  nb_color_per_elem = 0
  ier = 0
  num_spec = 0
  iproc = NPROC
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
  if (US_LETTER) then
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
  call set_color_palette(red,green,blue,NUM_COLORS)

  ! get minimum and maximum values of mesh coordinates
  xmin = minval(coord(1,:))
  zmin = minval(coord(2,:))
  xmax = maxval(coord(1,:))
  zmax = maxval(coord(2,:))

  call min_all_all_dp(xmin, xmin_glob)
  call min_all_all_dp(zmin, zmin_glob)
  call max_all_all_dp(xmax, xmax_glob)
  call max_all_all_dp(zmax, zmax_glob)
  xmin = xmin_glob
  zmin = zmin_glob
  xmax = xmax_glob
  zmax = zmax_glob

  if (myrank == 0) then
     write(IMAIN,*) '  X min, max = ',xmin,xmax
     write(IMAIN,*) '  Z min, max = ',zmin,zmax
  endif

  ! ratio of physical page size/size of the domain meshed
  ratio_page = min(rpercentz*sizez/(zmax-zmin),rpercentx*sizex/(xmax-xmin)) / 100.d0

  ! compute the maximum of the norm of the vector
  dispmax = maxval(sqrt(vector_field_display(1,:)**2 + vector_field_display(2,:)**2))

  call max_all_all_dp(dispmax, dispmax_glob)
  dispmax = dispmax_glob

  ! checks maximum value
! this trick checks for NaN (Not a Number), which is not even equal to itself
  if (dispmax > STABILITY_THRESHOLD .or. dispmax < 0 .or. dispmax /= dispmax) then
    print *,'Warning: failed creating postscript image, maximum value of display is invalid'
    print *,'display max = ',dispmax,' with threshold at ', STABILITY_THRESHOLD
    print *,'Please check your simulation setup...'
    call exit_MPI(myrank,'Code became unstable and blew up ( vector_field_display array in routine plot_post() )')
  endif

  ! user output
  if (myrank == 0) then
     write(IMAIN,*) '  Max norm = ',dispmax
     call flush_IMAIN()
  endif
  call synchronize_all()

!
!---- open PostScript file
!
  if (myrank == 0) then
    write(file_name,"(a,i7.7,a)") trim(OUTPUT_FILES)//'vect',it,'.ps'
    open(unit=24,file=file_name,status='unknown',iostat=ier)
    if (ier /= 0) call exit_MPI(myrank,'Error opening postscript file for image output')

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
    if (timeval >= 1.d-3 .and. timeval < 1000.d0) then
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
    if (DISPLAY_COLORS == 1) write(24,*) '.4 .9 .9 setrgbcolor'
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
    if (DISPLAY_COLORS == 1) write(24,*) '.8 0 .8 setrgbcolor'
    write(24,*) '24.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'

  if (DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then
    write(24,*) '(Deformed mesh) show'
  else
    if (imagetype_postscript == 1) then
      write(24,*) '(Displacement vector field) show'
    else if (imagetype_postscript == 2) then
      write(24,*) '(Velocity vector field) show'
    else if (imagetype_postscript == 3) then
      write(24,*) '(Acceleration vector field) show'
    else
      call exit_MPI(myrank,'Bad field code in PostScript display')
    endif
  endif

    write(24,*) 'grestore'
    write(24,*) '25.35 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'
    write(24,15) simulation_title
    write(24,*) 'grestore'
    write(24,*) '26.45 CM 18.9 CM MV'
    write(24,*) usoffset,' CM 2 div neg 0 MR'
    write(24,*) 'currentpoint gsave translate -90 rotate 0 0 moveto'

    if (coupled_acoustic_elastic) then
      write(24,*) '(Coupled Acoustic/Elastic Wave 2D - SEM) show'
    else if (coupled_acoustic_poro) then
      write(24,*) '(Coupled Acoustic/Poroelastic Wave 2D - SEM) show'
    else if (coupled_elastic_poro) then
      write(24,*) '(Coupled Elastic/Poroelastic Wave 2D - SEM) show'
    else if (any_acoustic) then
      write(24,*) '(Acoustic Wave 2D - Spectral Element Method) show'
    else if (any_poroelastic) then
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
  if (modelvect) then

    buffer_offset = 0
    RGB_offset = 0

    do ispec= 1,nspec
      do i = 1,NGLLX-subsamp_postscript,subsamp_postscript
        do j = 1,NGLLX-subsamp_postscript,subsamp_postscript

          if ((vpImax-vpImin)/vpImin > 0.02d0) then

            if (assign_external_model) then

              x1 = (vpext(i,j,ispec)-vpImin) / (vpImax-vpImin)

            else

              material = kmato(ispec)

              if (ispec_is_poroelastic(ispec)) then

                ! poroelastic material

                ! get elastic parameters of current spectral element
                call get_poroelastic_material(ispec,phi,tort,mu_s,kappa_s,rho_s,kappa_f,rho_f,eta_f,mu_fr,kappa_fr,rho_bar)

                ! Biot coefficients for the input phi
                call get_poroelastic_Biot_coeff(phi,kappa_s,kappa_f,kappa_fr,mu_fr,D_biot,H_biot,C_biot,M_biot)

                ! Approximated velocities (no viscous dissipation)
                afactor = rho_bar - phi/tort*rho_f
                bfactor = H_biot + phi*rho_bar/(tort*rho_f)*M_biot - 2.d0*phi/tort*C_biot
                cfactor = phi/(tort*rho_f)*(H_biot*M_biot - C_biot*C_biot)
                cpIsquare = (bfactor + sqrt(bfactor*bfactor - 4.d0*afactor*cfactor))/(2.d0*afactor)
                cpIloc = sqrt(cpIsquare)

              else

                lambdaplus2mu  = poroelastcoef(3,1,material)
                denst = density(1,material)
                cpIloc = sqrt(lambdaplus2mu/denst)

              endif

              x1 = (cpIloc-vpImin) / (vpImax-vpImin)

            endif

          else
            x1 = 0.5d0
          endif

          ! rescale to avoid very dark gray levels
          x1 = x1*0.7 + 0.2
          if (x1 > 1.d0) x1 = 1.d0

          ! invert scale: white = vpImin, dark gray = vpImax
          x1 = 1.d0 - x1

          xw = coord(1,ibool(i,j,ispec))
          zw = coord(2,ibool(i,j,ispec))
          xw = (xw-xmin)*ratio_page + orig_x
          zw = (zw-zmin)*ratio_page + orig_z
          xw = xw * centim
          zw = zw * centim
          if (myrank == 0) then
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
          if (myrank == 0) then
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
          if (myrank == 0) then
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
          if (myrank == 0) then
             write(24,499) xw,zw
          else
             buffer_offset = buffer_offset + 1
             coorg_send_ps_velocity_model(1,buffer_offset) = xw
             coorg_send_ps_velocity_model(2,buffer_offset) = zw
          endif

          ! display P-velocity model using gray levels
          if (myrank == 0) then
             write(24,604) x1
          else
             RGB_offset = RGB_offset + 1
             RGB_send_ps_velocity_model(1,RGB_offset) = x1
          endif

        enddo
      enddo
    enddo

#ifdef USE_MPI
    if (NPROC > 1) then
      if (myrank == 0) then
        ! master collects
        do iproc = 1, NPROC-1

          call recv_singlei(nspec_recv, iproc, 42)
          call recv_dp(coorg_recv_ps_velocity_model(1,1), &
               2*nspec_recv*((NGLLX-subsamp_postscript)/subsamp_postscript)*((NGLLX-subsamp_postscript)/subsamp_postscript)*4, &
               iproc, 42)
          call recv_dp(RGB_recv_ps_velocity_model(1,1), nspec_recv*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
               ((NGLLX-subsamp_postscript)/subsamp_postscript), &
               iproc, 42)

          buffer_offset = 0
          RGB_offset = 0
          do ispec = 1, nspec_recv
             do i = 1,NGLLX-subsamp_postscript,subsamp_postscript
                do j = 1,NGLLX-subsamp_postscript,subsamp_postscript
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
        call send_singlei(nspec, 0, 42)
        call send_dp(coorg_send_ps_velocity_model(1,1), 2*nspec*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
                     ((NGLLX-subsamp_postscript)/subsamp_postscript)*4, &
                     0, 42)
        call send_dp(RGB_send_ps_velocity_model(1,1), nspec*((NGLLX-subsamp_postscript)/subsamp_postscript)* &
                     ((NGLLX-subsamp_postscript)/subsamp_postscript), &
                     0, 42)
      endif
    endif
    call synchronize_all()
#endif

  endif ! modelvect

!
!---- draw the spectral element mesh
!

  if (myrank == 0) then
     write(24,*) '%'
     write(24,*) '% spectral element mesh'
     write(24,*) '%'
  endif

  buffer_offset = 0
  RGB_offset = 0
  color_send_ps_element_mesh(:) = 0

  do ispec = 1,nspec

    if (myrank == 0) write(24,*) '% elem ',ispec

    do i = 1,pointsdisp
      do j = 1,pointsdisp

        xinterp(i,j) = 0.d0
        zinterp(i,j) = 0.d0
        do in = 1,ngnod
          nnum = knods(in,ispec)
          xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
          zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
        enddo

        if (DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then

          Uxinterp(i,j) = 0.d0
          Uzinterp(i,j) = 0.d0

          do k = 1,NGLLX
            do l= 1,NGLLX
              if (AXISYM) then
                if (is_on_the_axis(ispec)) then
                  Uxinterp(i,j) = Uxinterp(i,j) + vector_field_display(1,ibool(k,l,ispec))*flagrange_GLJ(k,i)*flagrange_GLJ(l,j)
                  Uzinterp(i,j) = Uzinterp(i,j) + vector_field_display(2,ibool(k,l,ispec))*flagrange_GLJ(k,i)*flagrange_GLJ(l,j)
                else
                  Uxinterp(i,j) = Uxinterp(i,j) + vector_field_display(1,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
                  Uzinterp(i,j) = Uzinterp(i,j) + vector_field_display(2,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
                endif
              else
                Uxinterp(i,j) = Uxinterp(i,j) + vector_field_display(1,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
                Uzinterp(i,j) = Uzinterp(i,j) + vector_field_display(2,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
              endif
            enddo
          enddo

          ! to compute the deformed mesh we take the original mesh and add the displacement of the current point
          ! to the coordinates of that original (undeformed) mesh point, with a scaling factor that the user
          ! can make different from 1 if he/she wants in order to reduce or increase the visual effect of the deformation
          ! for display purposes
          xinterp(i,j) = xinterp(i,j) + Uxinterp(i,j) * SCALING_FACTOR_FOR_DEFORMED_MESH
          zinterp(i,j) = zinterp(i,j) + Uzinterp(i,j) * SCALING_FACTOR_FOR_DEFORMED_MESH

        endif ! of if (DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then

      enddo
    enddo

    is = 1
    ir = 1
    x1 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
    z1 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
    x1 = x1 * centim
    z1 = z1 * centim
    if (myrank == 0) then
      write(24,*) 'mark'
      write(24,681) x1,z1
    else
      buffer_offset = buffer_offset + 1
      coorg_send_ps_element_mesh(1,buffer_offset) = x1
      coorg_send_ps_element_mesh(2,buffer_offset) = z1
    endif

    if (ngnod == 4) then

      ! draw straight lines if elements have 4 nodes

      ir=pointsdisp
      x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
      z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
      x2 = x2 * centim
      z2 = z2 * centim
      if (myrank == 0) then
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
      if (myrank == 0) then
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
      if (myrank == 0) then
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
      if (myrank == 0) then
         write(24,681) x2,z2
      else
         buffer_offset = buffer_offset + 1
         coorg_send_ps_element_mesh(1,buffer_offset) = x2
         coorg_send_ps_element_mesh(2,buffer_offset) = z2
      endif

    else

      ! draw curved lines if elements have 9 nodes
      do ir = 2,pointsdisp
        x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
        z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
        x2 = x2 * centim
        z2 = z2 * centim
        if (myrank == 0) then
           write(24,681) x2,z2
        else
           buffer_offset = buffer_offset + 1
           coorg_send_ps_element_mesh(1,buffer_offset) = x2
           coorg_send_ps_element_mesh(2,buffer_offset) = z2
        endif
      enddo

      ir=pointsdisp
      do is= 2,pointsdisp
        x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
        z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
        x2 = x2 * centim
        z2 = z2 * centim
        if (myrank == 0) then
           write(24,681) x2,z2
        else
           buffer_offset = buffer_offset + 1
           coorg_send_ps_element_mesh(1,buffer_offset) = x2
           coorg_send_ps_element_mesh(2,buffer_offset) = z2
        endif
      enddo

      is=pointsdisp
      do ir =pointsdisp-1,1,-1
        x2 = (xinterp(ir,is)-xmin)*ratio_page + orig_x
        z2 = (zinterp(ir,is)-zmin)*ratio_page + orig_z
        x2 = x2 * centim
        z2 = z2 * centim
        if (myrank == 0) then
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
        if (myrank == 0) then
           write(24,681) x2,z2
        else
           buffer_offset = buffer_offset + 1
           coorg_send_ps_element_mesh(1,buffer_offset) = x2
           coorg_send_ps_element_mesh(2,buffer_offset) = z2
        endif
      enddo

    endif

    if (myrank == 0) then
       write(24,*) 'CO'
    endif

    if (DISPLAY_COLORS == 1) then

      ! use a different color for each material set
      imat = kmato(ispec)
      icol = mod(imat - 1,NUM_COLORS) + 1

      ! display all the PML layers in a different (constant) color if needed
      if (PML_BOUNDARY_CONDITIONS .and. ispec_is_PML(ispec)) then
        if (DISPLAY_PML_IN_DIFFERENT_COLOR) then
          icol = ICOLOR_FOR_PML_DISPLAY
          ! make sure that number exists
          if (icol > NUM_COLORS) icol = NUM_COLORS
        endif
      endif

      if ( myrank == 0) then
        if (meshvect) then
          write(24,680) red(icol),green(icol),blue(icol)
        else
          write(24,679) red(icol),green(icol),blue(icol)
        endif
      else
        RGB_offset = RGB_offset + 1
        color_send_ps_element_mesh(RGB_offset) = icol
      endif
    endif

    if (myrank == 0) then
      if (meshvect) then
        if (modelvect) then
          write(24,*) 'Colmesh ST'
        else
          write(24,*) '0 setgray ST'
        endif
      endif
    endif

    ! write the element number, the group number and the material number inside the element
    if (DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT == 1) then

      xw = (coorg(1,knods(1,ispec)) + coorg(1,knods(2,ispec)) + coorg(1,knods(3,ispec)) + coorg(1,knods(4,ispec))) / 4.d0
      zw = (coorg(2,knods(1,ispec)) + coorg(2,knods(2,ispec)) + coorg(2,knods(3,ispec)) + coorg(2,knods(4,ispec))) / 4.d0
      xw = (xw-xmin)*ratio_page + orig_x
      zw = (zw-zmin)*ratio_page + orig_z
      xw = xw * centim
      zw = zw * centim

      if (myrank == 0) then
        if (DISPLAY_COLORS == 1) write(24,*) '1 setgray'
      endif

      if (myrank == 0) then
        write(24,500) xw,zw
      else
        buffer_offset = buffer_offset + 1
        coorg_send_ps_element_mesh(1,buffer_offset) = x2
        coorg_send_ps_element_mesh(2,buffer_offset) = z2
      endif

      ! write spectral element number
      if (myrank == 0) then
        write(24,502) ispec
      else
        RGB_offset = RGB_offset + 1
        color_send_ps_element_mesh(RGB_offset) = ispec
      endif

    endif

  enddo ! ispec

#ifdef USE_MPI
  if (NPROC > 1) then
    if (myrank == 0) then
      ! master collects
      do iproc = 1, NPROC-1
        call recv_singlei(nspec_recv,iproc,43)

        nb_coorg_per_elem = 1
        if (DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT == 1) then
          nb_coorg_per_elem = nb_coorg_per_elem + 1
        endif
        if (ngnod == 4) then
          nb_coorg_per_elem = nb_coorg_per_elem + 4
        else
          nb_coorg_per_elem = nb_coorg_per_elem + 3*(pointsdisp-1)+(pointsdisp-2)
        endif
        call recv_dp(coorg_recv_ps_element_mesh(1,1), 2*nspec_recv*nb_coorg_per_elem, iproc, 43)

        nb_color_per_elem = 0
        if (DISPLAY_COLORS == 1) then
          nb_color_per_elem = nb_color_per_elem + 1
        endif
        if (DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT == 1) then
          nb_color_per_elem = nb_color_per_elem + 1
        endif
        if (nb_color_per_elem > 0) then
          call recv_i(color_recv_ps_element_mesh(1), nspec_recv*nb_color_per_elem, iproc, 43)
        endif

        buffer_offset = 0
        RGB_offset = 0
        num_spec = nspec
        do ispec = 1, nspec_recv
          num_spec = num_spec + 1
          write(24,*) '% elem ',num_spec
          buffer_offset = buffer_offset + 1
          write(24,*) 'mark'
          write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
          if (ngnod == 4) then
            buffer_offset = buffer_offset + 1
            write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
            buffer_offset = buffer_offset + 1
            write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
            buffer_offset = buffer_offset + 1
            write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
            buffer_offset = buffer_offset + 1
            write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
          else
            do ir = 2,pointsdisp
              buffer_offset = buffer_offset + 1
              write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
            enddo
            do is= 2,pointsdisp
              buffer_offset = buffer_offset + 1
              write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
            enddo
            do ir =pointsdisp-1,1,-1
              buffer_offset = buffer_offset + 1
              write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
            enddo
            do is=pointsdisp-1,2,-1
              buffer_offset = buffer_offset + 1
              write(24,681) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
            enddo
          endif

          write(24,*) 'CO'
          if (DISPLAY_COLORS == 1) then
            if (meshvect) then
              RGB_offset = RGB_offset + 1
              write(24,680) red(color_recv_ps_element_mesh(RGB_offset)), &
                            green(color_recv_ps_element_mesh(RGB_offset)), &
                            blue(color_recv_ps_element_mesh(RGB_offset))
            else
              RGB_offset = RGB_offset + 1
              write(24,679) red(color_recv_ps_element_mesh(RGB_offset)), &
                            green(color_recv_ps_element_mesh(RGB_offset)), &
                            blue(color_recv_ps_element_mesh(RGB_offset))
            endif
          endif
          if (meshvect) then
            if (modelvect) then
              write(24,*) 'Colmesh ST'
            else
              write(24,*) '0 setgray ST'
            endif
          endif
          if (DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT == 1) then
            if (DISPLAY_COLORS == 1) write(24,*) '1 setgray'
            buffer_offset = buffer_offset + 1
            write(24,500) coorg_recv_ps_element_mesh(1,buffer_offset), coorg_recv_ps_element_mesh(2,buffer_offset)
            RGB_offset = RGB_offset + 1
            write(24,502) color_recv_ps_element_mesh(RGB_offset)
          endif
        enddo
      enddo
    else
      call send_singlei(nspec, 0, 43)

      nb_coorg_per_elem = 1
      if (DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT == 1) then
        nb_coorg_per_elem = nb_coorg_per_elem + 1
      endif
      if (ngnod == 4) then
        nb_coorg_per_elem = nb_coorg_per_elem + 4
      else
        nb_coorg_per_elem = nb_coorg_per_elem + 3*(pointsdisp-1)+(pointsdisp-2)
      endif
      call send_dp(coorg_send_ps_element_mesh(1,1), 2*nspec*nb_coorg_per_elem, 0, 43)

      nb_color_per_elem = 0
      if (DISPLAY_COLORS == 1) then
        nb_color_per_elem = nb_color_per_elem + 1
      endif
      if (DISPLAY_ELEMENT_NUMBERS_POSTSCRIPT == 1) then
        nb_color_per_elem = nb_color_per_elem + 1
      endif
      if (nb_color_per_elem > 0) then
        call send_i(color_send_ps_element_mesh(1), nspec*nb_color_per_elem, 0, 43)
      endif
    endif
  endif
  call synchronize_all()
#endif

!
!--- draw absorbing boundaries with a thick color line
!
  ! sets global flag for all slices
  call any_all_l(anyabs, anyabs_glob)

  if (anyabs_glob .and. boundvect .and. .not. DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then

    if (myrank == 0) then
      write(24,*) '%'
      write(24,*) '% boundary conditions on the mesh'
      write(24,*) '%'

      write(24,*) '0.05 CM setlinewidth'
    endif

    buffer_offset = 0

    if (anyabs) then
      do inum = 1,nelemabs
        ispec = numabs(inum)

        do iedge = 1,4

          if (codeabs(iedge,inum)) then ! codeabs(:,:) is defined as "logical" in MAIN program

            if (iedge == IEDGE1) then
              ! bottom
              ideb = 1
              ifin = 2
            else if (iedge == IEDGE2) then
              ! right
              ideb = 2
              ifin = 3
            else if (iedge == IEDGE3) then
              ! top
              ideb = 3
              ifin = 4
            else if (iedge == IEDGE4) then
              ! left
              ideb = 4
              ifin = 1
            else
              call exit_MPI(myrank,'Wrong absorbing boundary code')
            endif

            x1 = (coorg(1,knods(ideb,ispec))-xmin)*ratio_page + orig_x
            z1 = (coorg(2,knods(ideb,ispec))-zmin)*ratio_page + orig_z
            x2 = (coorg(1,knods(ifin,ispec))-xmin)*ratio_page + orig_x
            z2 = (coorg(2,knods(ifin,ispec))-zmin)*ratio_page + orig_z
            x1 = x1 * centim
            z1 = z1 * centim
            x2 = x2 * centim
            z2 = z2 * centim

            if (myrank == 0) then
            ! draw the Stacey absorbing boundary line segment in different colors depending on its type
              if (typeabs(inum) == IBOTTOM) then
                write(24,*) '0 1 0 RG'  ! green
              else if (typeabs(inum) == IRIGHT) then
                write(24,*) '0 0 1 RG'  ! blue
              else if (typeabs(inum) == ITOP) then
                write(24,*) '1 0.7529 0.7960 RG' ! pink
              else if (typeabs(inum) == ILEFT) then
                write(24,*) '1 0.6470 0 RG' ! orange
              else
                call exit_MPI(myrank,'Wrong absorbing boundary code')
              endif
              write(24,602) x1,z1,x2,z2
            else
              buffer_offset = buffer_offset + 1
              coorg_send_ps_abs(1,buffer_offset) = x1
              coorg_send_ps_abs(2,buffer_offset) = z1
              coorg_send_ps_abs(3,buffer_offset) = x2
              coorg_send_ps_abs(4,buffer_offset) = z2
              coorg_send_ps_abs(5,buffer_offset) = typeabs(inum)
            endif

          endif ! of if (codeabs(iedge,inum))
        enddo ! of do iedge = 1,4

      enddo
    endif ! anyabs

#ifdef USE_MPI
    if (NPROC > 1) then
      if (myrank == 0) then
        ! master collects
        do iproc = 1, NPROC-1
          call recv_singlei(nspec_recv,iproc, 44)
          if (nspec_recv > 0) then
            call recv_dp(coorg_recv_ps_abs(1,1), 5*nspec_recv, iproc, 44)

            buffer_offset = 0
            do ispec = 1, nspec_recv
              buffer_offset = buffer_offset + 1
            ! draw the Stacey absorbing boundary line segment in different colors depending on its type
              if (coorg_recv_ps_abs(5,buffer_offset) == IBOTTOM) then
                write(24,*) '0 1 0 RG'  ! green
              else if (coorg_recv_ps_abs(5,buffer_offset) == IRIGHT) then
                write(24,*) '0 0 1 RG'  ! blue
              else if (coorg_recv_ps_abs(5,buffer_offset) == ITOP) then
                write(24,*) '1 0.7529 0.7960 RG' ! pink
              else if (coorg_recv_ps_abs(5,buffer_offset) == ILEFT) then
                write(24,*) '1 0.6470 0 RG' ! orange
              else
                call exit_MPI(myrank,'Wrong absorbing boundary code')
              endif
              write(24,602) coorg_recv_ps_abs(1,buffer_offset), coorg_recv_ps_abs(2,buffer_offset), &
                            coorg_recv_ps_abs(3,buffer_offset), coorg_recv_ps_abs(4,buffer_offset)
            enddo
          endif
        enddo
      else
        call send_singlei(buffer_offset, 0, 44)
        if (buffer_offset > 0) then
          call send_dp(coorg_send_ps_abs(1,1), 5*buffer_offset, 0, 44)
        endif
      endif
    endif
    call synchronize_all()
#endif

    if (myrank == 0) then
      write(24,*) '0 setgray'
      write(24,*) '0 setlinewidth'
    endif

  endif ! anyabs_glob .and. boundvect

!
!--- draw free surface with a thick color line
!

  if (.not. DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then

  if (myrank == 0) then
    write(24,*) '%'
    write(24,*) '% free surface on the mesh'
    write(24,*) '%'

    ! use orange color
    write(24,*) '1 0.66 0 RG'
    write(24,*) '0.05 CM setlinewidth'
  endif

  buffer_offset = 0

  if (nelem_acoustic_surface > 0) then
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
      if (myrank == 0) then
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
  if (NPROC > 1) then
    if (myrank == 0) then
      ! master collects
      do iproc = 1, NPROC-1
        call recv_singlei(nspec_recv, iproc, 44)
        if (nspec_recv > 0) then
          call recv_dp(coorg_recv_ps_free_surface(1,1), 4*nspec_recv,iproc, 44)

          buffer_offset = 0
          do ispec = 1, nspec_recv
            buffer_offset = buffer_offset + 1
            write(24,602) coorg_recv_ps_free_surface(1,buffer_offset), coorg_recv_ps_free_surface(2,buffer_offset), &
                          coorg_recv_ps_free_surface(3,buffer_offset), coorg_recv_ps_free_surface(4,buffer_offset)
          enddo
        endif
      enddo
    else
      call send_singlei(buffer_offset, 0, 44)
      if (buffer_offset > 0) then
        call send_dp(coorg_send_ps_free_surface(1,1), 4*buffer_offset, 0, 44)
      endif
    endif
  endif
  call synchronize_all()
#endif

  if (myrank == 0) then
    write(24,*) '0 setgray'
    write(24,*) '0 setlinewidth'
  endif

  endif ! of if (.not. DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then

!
!----  draw the fluid-solid coupling edges with a thick color line
!
  ! sets global flag for all slices
  call any_all_l(coupled_acoustic_elastic, coupled_acoustic_elastic_glob)

  if (coupled_acoustic_elastic_glob .and. boundvect .and. .not. DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then

    if (myrank == 0) then
      write(24,*) '%'
      write(24,*) '% fluid-solid coupling edges in the mesh'
      write(24,*) '%'

      ! use grey color
      write(24,*) '0.65 0.65 0.65 RG'
      write(24,*) '0.05 CM setlinewidth'
    endif

    if (myrank /= 0 .and. num_fluid_solid_edges > 0 ) allocate(coorg_send(4,num_fluid_solid_edges))
    buffer_offset = 0

    ! loop on all the coupling edges
    do inum = 1,num_fluid_solid_edges

      ! get the edge of the acoustic element
      ispec = fluid_solid_acoustic_ispec(inum)
      iedge = fluid_solid_acoustic_iedge(inum)

      if (iedge == ITOP) then
        ideb = 3
        ifin = 4
      else if (iedge == IBOTTOM) then
        ideb = 1
        ifin = 2
      else if (iedge == ILEFT) then
        ideb = 4
        ifin = 1
      else if (iedge == IRIGHT) then
        ideb = 2
        ifin = 3
      else
        call exit_MPI(myrank,'Wrong fluid-solid coupling edge code')
      endif

      x1 = (coorg(1,knods(ideb,ispec))-xmin)*ratio_page + orig_x
      z1 = (coorg(2,knods(ideb,ispec))-zmin)*ratio_page + orig_z
      x2 = (coorg(1,knods(ifin,ispec))-xmin)*ratio_page + orig_x
      z2 = (coorg(2,knods(ifin,ispec))-zmin)*ratio_page + orig_z
      x1 = x1 * centim
      z1 = z1 * centim
      x2 = x2 * centim
      z2 = z2 * centim
      if (myrank == 0) then
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
    if (NPROC > 1) then
      if (myrank == 0) then
        ! master collects
        do iproc = 1, NPROC-1
          call recv_singlei(nspec_recv, iproc, 45)
          if (nspec_recv > 0) then
            allocate(coorg_recv(4,nspec_recv))
            call recv_dp(coorg_recv(1,1), 4*nspec_recv, iproc, 45)

            buffer_offset = 0
            do ispec = 1, nspec_recv
              buffer_offset = buffer_offset + 1
              write(24,602) coorg_recv(1,buffer_offset), coorg_recv(2,buffer_offset), &
                            coorg_recv(3,buffer_offset), coorg_recv(4,buffer_offset)
            enddo
            deallocate(coorg_recv)
          endif
        enddo
      else
        call send_singlei(buffer_offset, 0, 45)
        if (buffer_offset > 0) then
          call send_dp(coorg_send(1,1), 4*buffer_offset, 0, 45)
          deallocate(coorg_send)
        endif
      endif
    endif
    call synchronize_all()
#endif

    if (myrank == 0) then
      write(24,*) '0 setgray'
      write(24,*) '0 setlinewidth'
    endif

  endif ! coupled_acoustic_elastic_glob .and. boundvect

!
!----  draw the fluid-porous coupling edges with a thick color line
!
  ! sets global flag for all slices
  call any_all_l(coupled_acoustic_poro, coupled_acoustic_poro_glob)

  if (coupled_acoustic_poro_glob .and. boundvect .and. .not. DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then

    if (myrank == 0) then
      write(24,*) '%'
      write(24,*) '% fluid-porous coupling edges in the mesh'
      write(24,*) '%'

      ! use grey color
      write(24,*) '0.65 0.65 0.65 RG'
      write(24,*) '0.05 CM setlinewidth'
    endif

    if (myrank /= 0 .and. num_fluid_poro_edges > 0 ) allocate(coorg_send(4,num_fluid_poro_edges))
    buffer_offset = 0

    ! loop on all the coupling edges
    do inum = 1,num_fluid_poro_edges

      ! get the edge of the acoustic element
      ispec = fluid_poro_acoustic_ispec(inum)
      iedge = fluid_poro_acoustic_iedge(inum)

      if (iedge == ITOP) then
        ideb = 3
        ifin = 4
      else if (iedge == IBOTTOM) then
        ideb = 1
        ifin = 2
      else if (iedge == ILEFT) then
        ideb = 4
        ifin = 1
      else if (iedge == IRIGHT) then
        ideb = 2
        ifin = 3
      else
        call exit_MPI(myrank,'Wrong fluid-solid coupling edge code')
      endif

      x1 = (coorg(1,knods(ideb,ispec))-xmin)*ratio_page + orig_x
      z1 = (coorg(2,knods(ideb,ispec))-zmin)*ratio_page + orig_z
      x2 = (coorg(1,knods(ifin,ispec))-xmin)*ratio_page + orig_x
      z2 = (coorg(2,knods(ifin,ispec))-zmin)*ratio_page + orig_z
      x1 = x1 * centim
      z1 = z1 * centim
      x2 = x2 * centim
      z2 = z2 * centim
      if (myrank == 0) then
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
    if (NPROC > 1) then
      if (myrank == 0) then
        ! master collects
        do iproc = 1, NPROC-1
          call recv_singlei(nspec_recv, iproc, 45)
          if (nspec_recv > 0) then
            allocate(coorg_recv(4,nspec_recv))
            call recv_dp(coorg_recv(1,1), 4*nspec_recv, iproc, 45)

            buffer_offset = 0
            do ispec = 1, nspec_recv
              buffer_offset = buffer_offset + 1
              write(24,602) coorg_recv(1,buffer_offset), coorg_recv(2,buffer_offset), &
                            coorg_recv(3,buffer_offset), coorg_recv(4,buffer_offset)
            enddo
            deallocate(coorg_recv)
          endif
        enddo
      else
        call send_singlei(buffer_offset, 0, 45)
        if (buffer_offset > 0) then
          call send_dp(coorg_send(1,1), 4*buffer_offset, 0, 45)
          deallocate(coorg_send)
        endif
      endif
    endif
    call synchronize_all()
#endif

    if (myrank == 0) then
      write(24,*) '0 setgray'
      write(24,*) '0 setlinewidth'
    endif

  endif ! coupled_acoustic_poro_glob .and. boundvect

!
!----  draw the solid-porous coupling edges with a thick color line
!
  ! sets global flag for all slices
  call any_all_l(coupled_elastic_poro, coupled_elastic_poro_glob)


  if (coupled_elastic_poro_glob .and. boundvect .and. .not. DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then

    if (myrank == 0) then
      write(24,*) '%'
      write(24,*) '% solid-porous coupling edges in the mesh'
      write(24,*) '%'

      ! use grey color
      write(24,*) '0.65 0.65 0.65 RG'
      write(24,*) '0.05 CM setlinewidth'
    endif

    if (myrank /= 0 .and. num_solid_poro_edges > 0 ) allocate(coorg_send(4,num_solid_poro_edges))
    buffer_offset = 0

    ! loop on all the coupling edges
    do inum = 1,num_solid_poro_edges

      ! get the edge of the poroelastic element
      ispec = solid_poro_poroelastic_ispec(inum)
      iedge = solid_poro_poroelastic_iedge(inum)

      if (iedge == ITOP) then
        ideb = 3
        ifin = 4
      else if (iedge == IBOTTOM) then
        ideb = 1
        ifin = 2
      else if (iedge == ILEFT) then
        ideb = 4
        ifin = 1
      else if (iedge == IRIGHT) then
        ideb = 2
        ifin = 3
      else
        call exit_MPI(myrank,'Wrong fluid-solid coupling edge code')
      endif

      x1 = (coorg(1,knods(ideb,ispec))-xmin)*ratio_page + orig_x
      z1 = (coorg(2,knods(ideb,ispec))-zmin)*ratio_page + orig_z
      x2 = (coorg(1,knods(ifin,ispec))-xmin)*ratio_page + orig_x
      z2 = (coorg(2,knods(ifin,ispec))-zmin)*ratio_page + orig_z
      x1 = x1 * centim
      z1 = z1 * centim
      x2 = x2 * centim
      z2 = z2 * centim
      if (myrank == 0) then
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
    if (NPROC > 1) then
      if (myrank == 0) then
        ! master collects
        do iproc = 1, NPROC-1
          call recv_singlei(nspec_recv, iproc, 45)
          if (nspec_recv > 0) then
            allocate(coorg_recv(4,nspec_recv))
            call recv_dp(coorg_recv(1,1), 4*nspec_recv, iproc, 45)

            buffer_offset = 0
            do ispec = 1, nspec_recv
              buffer_offset = buffer_offset + 1
              write(24,602) coorg_recv(1,buffer_offset), coorg_recv(2,buffer_offset), &
                            coorg_recv(3,buffer_offset), coorg_recv(4,buffer_offset)
            enddo
            deallocate(coorg_recv)
          endif
        enddo
      else
        call send_singlei(buffer_offset, 0, 45)
        if (buffer_offset > 0) then
          call send_dp(coorg_send(1,1), 4*buffer_offset, 0, 45)
          deallocate(coorg_send)
        endif
      endif
    endif
    call synchronize_all()
#endif

    if (myrank == 0) then
      write(24,*) '0 setgray'
      write(24,*) '0 setlinewidth'
    endif

  endif ! coupled_elastic_poro_glob .and. boundvect

!
!----  draw the normalized vector field
!

  if (.not. DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then

  ! warning if the maximum vector equals zero (no source)
  if (dispmax == 0.d0) then
    if (myrank == 0) then
      write(IMAIN,*) 'null vector: maximum vector length is zero!'
      call flush_IMAIN()
    endif
  endif

  if (myrank == 0) then
    write(24,*) '%'
    write(24,*) '% vector field'
    write(24,*) '%'

    ! color arrows if we draw the velocity model in the background
    if (modelvect) then
      write(24,*) 'Colvects'
    else
      write(24,*) '0 setgray'
    endif
  endif

  if (interpol) then

    if (myrank == 0) then
      write(IMAIN,*) '  Interpolating the vector field...'
      call flush_IMAIN()
    endif

    ! option to plot only lowerleft corner value to avoid very large files if dense meshes
    if (plot_lowerleft_corner_only) then
      pointsdisp_loop = 1
    else
      pointsdisp_loop = pointsdisp
    endif

    buffer_offset = 0

    do ispec = 1,nspec

! interpolation on a uniform grid
#ifdef USE_MPI
      if (myrank == 0 .and. mod(ispec,1000) == 0) &
        write(IMAIN,*) '  Interpolation uniform grid element ',ispec,' on processor core 0'
#else
      if (mod(ispec,1000) == 0) &
        write(IMAIN,*) '  Interpolation uniform grid element ',ispec
#endif

      do i = 1,pointsdisp_loop
        do j = 1,pointsdisp_loop

          xinterp(i,j) = 0.d0
          zinterp(i,j) = 0.d0
          do in = 1,ngnod
            nnum = knods(in,ispec)
            xinterp(i,j) = xinterp(i,j) + shape2D_display(in,i,j)*coorg(1,nnum)
            zinterp(i,j) = zinterp(i,j) + shape2D_display(in,i,j)*coorg(2,nnum)
          enddo

          Uxinterp(i,j) = 0.d0
          Uzinterp(i,j) = 0.d0

          do k = 1,NGLLX
            do l= 1,NGLLX
              if (AXISYM) then
                if (is_on_the_axis(ispec)) then
                  Uxinterp(i,j) = Uxinterp(i,j) + vector_field_display(1,ibool(k,l,ispec))*flagrange_GLJ(k,i)*flagrange_GLJ(l,j)
                  Uzinterp(i,j) = Uzinterp(i,j) + vector_field_display(2,ibool(k,l,ispec))*flagrange_GLJ(k,i)*flagrange_GLJ(l,j)
                else
                  Uxinterp(i,j) = Uxinterp(i,j) + vector_field_display(1,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
                  Uzinterp(i,j) = Uzinterp(i,j) + vector_field_display(2,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
                endif
              else
                Uxinterp(i,j) = Uxinterp(i,j) + vector_field_display(1,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
                Uzinterp(i,j) = Uzinterp(i,j) + vector_field_display(2,ibool(k,l,ispec))*flagrange(k,i)*flagrange(l,j)
              endif
            enddo
          enddo

          x1 =(xinterp(i,j)-xmin)*ratio_page
          z1 =(zinterp(i,j)-zmin)*ratio_page

          if (dispmax > 0.d0) then
            x2 = Uxinterp(i,j)*sizemax_arrows/dispmax
            z2 = Uzinterp(i,j)*sizemax_arrows/dispmax
          else
            x2 = 0.d0
            z2 = 0.d0
          endif
          d = sqrt(x2**2 + z2**2)

          ! ignore if vector is too small
          if (d > cutsnaps*sizemax_arrows) then

            d1 = d * ARROW_RATIO
            d2 = d1 * cos(ARROW_ANGLE*convert)

            dummy = x2/d
            if (dummy > 0.9999d0) dummy = 0.9999d0
            if (dummy < -0.9999d0) dummy = -0.9999d0
            theta = acos(dummy)

            if (z2 < 0.d0) theta = 360.d0*convert - theta
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
            if (myrank == 0) then
              write(postscript_line,700) xb,zb,xa,za,x2,z2,x1,z1
              ! suppress useless white spaces to make PostScript file smaller
              ! suppress leading white spaces again, if any
              postscript_line = adjustl(postscript_line)

              line_length = len_trim(postscript_line)
              index_char = 1
              first = .false.
              do ii = 1,line_length-1
                if (ch1(ii) /= ' ' .or. first) then
                  if (ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
                    ch2(index_char) = ch1(ii)
                    index_char = index_char + 1
                    first = .true.
                  endif
                endif
              enddo
              ch2(index_char) = ch1(line_length)
              write(24,"(100(a1))") (ch2(ii), ii= 1,index_char)
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

    enddo ! ispec

#ifdef USE_MPI
    if (myrank == 0) then
      ! master collects
      do iproc = 1, NPROC-1
        call recv_singlei(nspec_recv, iproc, 46)
        if (nspec_recv > 0) then
          call recv_dp(coorg_recv_ps_vector_field(1,1), 8*nspec_recv, iproc, 46)

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
              if (ch1(ii) /= ' ' .or. first) then
                if (ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
                  ch2(index_char) = ch1(ii)
                  index_char = index_char + 1
                  first = .true.
                endif
              endif
            enddo
            ch2(index_char) = ch1(line_length)
            write(24,"(100(a1))") (ch2(ii), ii= 1,index_char)
          enddo
        endif
      enddo
    else
      call send_singlei(buffer_offset, 0, 46)

      if (buffer_offset > 0) then
        call send_dp(coorg_send_ps_vector_field(1,1), 8*buffer_offset, 0, 46)

      endif
    endif
    call synchronize_all()
#endif

! draw the vectors at the nodes of the mesh if we do not interpolate the display on a regular grid
  else

    buffer_offset = 0

    do ipoin= 1,nglob

      x1 =(coord(1,ipoin)-xmin)*ratio_page
      z1 =(coord(2,ipoin)-zmin)*ratio_page

      if (dispmax > 0.d0) then
        x2 = vector_field_display(1,ipoin)*sizemax_arrows/dispmax
        z2 = vector_field_display(2,ipoin)*sizemax_arrows/dispmax
      else
        x2 = 0.d0
        z2 = 0.d0
      endif
      d = sqrt(x2**2 + z2**2)

      ! ignore if vector is too small
      if (d > cutsnaps*sizemax_arrows) then

        d1 = d * ARROW_RATIO
        d2 = d1 * cos(ARROW_ANGLE*convert)

        dummy = x2/d
        if (dummy > 0.9999d0) dummy = 0.9999d0
        if (dummy < -0.9999d0) dummy = -0.9999d0
        theta = acos(dummy)

        if (z2 < 0.d0) theta = 360.d0*convert - theta
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
        if (myrank == 0) then
          write(postscript_line,700) xb,zb,xa,za,x2,z2,x1,z1

          ! suppress useless white spaces to make PostScript file smaller
          ! suppress leading white spaces again, if any
          postscript_line = adjustl(postscript_line)

          line_length = len_trim(postscript_line)
          index_char = 1
          first = .false.
          do ii = 1,line_length-1
            if (ch1(ii) /= ' ' .or. first) then
              if (ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
                ch2(index_char) = ch1(ii)
                index_char = index_char + 1
                first = .true.
              endif
            endif
          enddo
          ch2(index_char) = ch1(line_length)
          write(24,"(100(a1))") (ch2(ii), ii= 1,index_char)
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
    if (myrank == 0) then
      ! master collects
      do iproc = 1, NPROC-1
        call recv_singlei(nspec_recv, iproc, 47)

        if (nspec_recv > 0) then
          call recv_dp(coorg_recv_ps_vector_field(1,1), 8*nspec_recv, iproc, 47)

          buffer_offset = 0
          do ispec = 1, nspec_recv
            buffer_offset = buffer_offset + 1
            write(postscript_line,700) coorg_recv_ps_vector_field(1,buffer_offset),coorg_recv_ps_vector_field(2,buffer_offset), &
                                       coorg_recv_ps_vector_field(3,buffer_offset),coorg_recv_ps_vector_field(4,buffer_offset), &
                                       coorg_recv_ps_vector_field(5,buffer_offset),coorg_recv_ps_vector_field(6,buffer_offset), &
                                       coorg_recv_ps_vector_field(7,buffer_offset),coorg_recv_ps_vector_field(8,buffer_offset)

            ! suppress useless white spaces to make PostScript file smaller
            ! suppress leading white spaces again, if any
            postscript_line = adjustl(postscript_line)

            line_length = len_trim(postscript_line)
            index_char = 1
            first = .false.
            do ii = 1,line_length-1
              if (ch1(ii) /= ' ' .or. first) then
                if (ch1(ii) /= ' ' .or. ch1(ii+1) /= ' ') then
                  ch2(index_char) = ch1(ii)
                  index_char = index_char + 1
                  first = .true.
                endif
              endif
            enddo
            ch2(index_char) = ch1(line_length)
            write(24,"(100(a1))") (ch2(ii), ii= 1,index_char)
          enddo
        endif
      enddo
    else
      call send_singlei(buffer_offset, 0, 47)
      if (buffer_offset > 0) then
        call send_dp(coorg_send_ps_vector_field(1,1), 8*buffer_offset, 0, 47)
      endif
    endif
    call synchronize_all()
#endif

  endif ! of interpolated values or values at GLL points

  endif ! of if (.not. DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR) then

  ! source/receiver marks
  if (myrank == 0) then
    write(24,*) '0 setgray'

    ! sources and receivers in color if velocity model
    if (modelvect) then
      write(24,*) 'Colreceiv'
    else
      write(24,*) '0 setgray'
    endif

    !
    !----  write position of the source
    !
    do i = 1,NSOURCES
      if (i == 1) write(24,*) '% beginning of source line'
      if (i == NSOURCES) write(24,*) '% end of source line'
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

  if (.not. (DISPLAY_DEFORMED_MESH_INSTEAD_OF_DISPLACEMENT_VECTOR .and. SUPPRESS_DISPLAY_OF_RECEIVERS_IF_DEFORMED_MESH)) then
    do i = 1,nrec
      if (i == 1) write(24,*) '% beginning of receiver line'
      if (i == nrec) write(24,*) '% end of receiver line'

      xw = st_xval(i)
      zw = st_zval(i)

      xw = (xw-xmin)*ratio_page + orig_x
      zw = (zw-zmin)*ratio_page + orig_z
      xw = xw * centim
      zw = zw * centim
      write(24,500) xw,zw
      write(24,*) 'Diamond'
    enddo
  endif

    write(24,*) '%'
    write(24,*) 'grestore'
    write(24,*) 'showpage'

    close(24)

  endif

  call synchronize_all()

 10  format('%!PS-Adobe-2.0',/,'%%',/,'%% Title: ',a100,/,'%% Created by: Specfem2D',/,'%% Author: Dimitri Komatitsch',/,'%%')
 15  format('(',a100,') show')
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

  end subroutine plot_post


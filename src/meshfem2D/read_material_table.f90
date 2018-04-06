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

  subroutine read_material_table()

! reads in material definitions in DATA/Par_file

  use constants, only: IMAIN,TINYVAL,ISOTROPIC_MATERIAL,ANISOTROPIC_MATERIAL,POROELASTIC_MATERIAL

  use shared_parameters, only: AXISYM,nbmodels,icodemat,cp,cs, &
                              aniso3,aniso4,aniso5,aniso6,aniso7,aniso8,aniso9,aniso10,aniso11,aniso12, &
                              comp_g,QKappa,Qmu, &
                              rho_s_read,rho_f_read, &
                              phi_read,tortuosity_read, &
                              permxx_read,permxz_read,permzz_read,kappa_s_read,kappa_f_read,kappa_fr_read, &
                              eta_f_read,mu_fr_read

  implicit none

  ! local parameters
  integer :: imaterial,i,icodematread,number_of_materials_defined_by_tomo_file
  integer :: reread_nbmodels
  double precision :: val0read,val1read,val2read,val3read,val4read, &
                      val5read,val6read,val7read,val8read,val9read,val10read,val11read,val12read

  integer,external :: err_occurred


  ! re-read number of models to reposition read-header
  call read_value_integer_p(reread_nbmodels, 'mesher.nbmodels')
  if (err_occurred() /= 0) call stop_the_code('error reading parameter nbmodels in Par_file')

  ! check
  if (reread_nbmodels /= nbmodels) call stop_the_code('Invalid reread value of nbmodels in reading material table')

  ! safety check
  if (nbmodels <= 0) call stop_the_code('Non-positive number of materials not allowed!')

  ! allocates material tables
  allocate(icodemat(nbmodels))
  allocate(cp(nbmodels))
  allocate(cs(nbmodels))

  allocate(aniso3(nbmodels))
  allocate(aniso4(nbmodels))
  allocate(aniso5(nbmodels))
  allocate(aniso6(nbmodels))
  allocate(aniso7(nbmodels))
  allocate(aniso8(nbmodels))
  allocate(aniso9(nbmodels))
  allocate(aniso10(nbmodels))
  allocate(aniso11(nbmodels))
  allocate(aniso12(nbmodels))

  allocate(comp_g(nbmodels))
  allocate(QKappa(nbmodels))
  allocate(Qmu(nbmodels))

  allocate(rho_s_read(nbmodels))
  allocate(rho_f_read(nbmodels))

  allocate( phi_read(nbmodels), &
            tortuosity_read(nbmodels), &
            permxx_read(nbmodels), &
            permxz_read(nbmodels), &
            permzz_read(nbmodels), &
            kappa_s_read(nbmodels), &
            kappa_f_read(nbmodels), &
            kappa_fr_read(nbmodels), &
            eta_f_read(nbmodels), &
            mu_fr_read(nbmodels))

  ! initializes material properties
  icodemat(:) = 0

  cp(:) = 0.d0
  cs(:) = 0.d0

  aniso3(:) = 0.d0
  aniso4(:) = 0.d0
  aniso5(:) = 0.d0
  aniso6(:) = 0.d0
  aniso7(:) = 0.d0
  aniso8(:) = 0.d0
  aniso9(:) = 0.d0
  aniso10(:) = 0.d0
  aniso11(:) = 0.d0

  comp_g(:) = 0.0d0
  QKappa(:) = 9999.d0
  Qmu(:) = 9999.d0

  rho_s_read(:) = 0.d0
  rho_f_read(:) = 0.d0

  phi_read(:) = 0.d0
  tortuosity_read(:) = 0.d0
  permxx_read(:) = 0.d0
  permxz_read(:) = 0.d0
  permzz_read(:) = 0.d0
  kappa_s_read(:) = 0.d0
  kappa_f_read(:) = 0.d0
  kappa_fr_read(:) = 0.d0
  eta_f_read(:) = 0.d0
  mu_fr_read(:) = 0.d0

  number_of_materials_defined_by_tomo_file = 0

  ! reads in material definitions
  do imaterial = 1,nbmodels

    call read_material_parameters_p(i,icodematread, &
                                    val0read,val1read,val2read,val3read, &
                                    val4read,val5read,val6read,val7read, &
                                    val8read,val9read,val10read,val11read,val12read)

    ! checks material id
    if (i < 1 .or. i > nbmodels) call stop_the_code('Wrong material number!')
    icodemat(i) = icodematread

    ! sets material properties
    if (icodemat(i) == ISOTROPIC_MATERIAL) then
      ! isotropic materials
      rho_s_read(i) = val0read
      cp(i) = val1read
      cs(i) = val2read
      comp_g(i) = val3read
      QKappa(i) = val5read
      Qmu(i) = val6read

      ! for Cs we use a less restrictive test because acoustic media have Cs exactly equal to zero
      if (rho_s_read(i) <= 0.00000001d0 .or. cp(i) <= 0.00000001d0 .or. cs(i) < 0.d0) then
        write(*,*) 'rho,cp,cs=',rho_s_read(i),cp(i),cs(i)
        call stop_the_code('negative value of velocity or density')
      endif
      if (QKappa(i) <= 0.00000001d0 .or. Qmu(i) <= 0.00000001d0) then
        call stop_the_code('non-positive value of QKappa or Qmu')
      endif

      aniso3(i) = val3read
      aniso4(i) = val4read
      if (abs(cs(i)) > TINYVAL) then
        phi_read(i) = 0.d0           ! elastic
      else
        phi_read(i) = 1.d0           ! acoustic
      endif

    else if (icodemat(i) == ANISOTROPIC_MATERIAL) then
      ! anisotropic materials
      rho_s_read(i) = val0read
      aniso3(i) = val1read
      aniso4(i) = val2read
      aniso5(i) = val3read
      aniso6(i) = val4read
      aniso7(i) = val5read
      aniso8(i) = val6read
      aniso9(i) = val7read
      aniso10(i) = val8read
      aniso11(i) = val9read
      aniso12(i) = val10read ! This value will be used only in AXISYM

    else if (icodemat(i) == POROELASTIC_MATERIAL) then
      ! poroelastic materials
      rho_s_read(i) = val0read
      rho_f_read(i) = val1read
      phi_read(i) = val2read
      tortuosity_read(i) = val3read
      permxx_read(i) = val4read
      permxz_read(i) = val5read
      permzz_read(i) = val6read
      kappa_s_read(i) = val7read
      kappa_f_read(i) = val8read
      kappa_fr_read(i) = val9read
      eta_f_read(i) = val10read
      mu_fr_read(i) = val11read
      Qmu(i) = val12read

      if (rho_s_read(i) <= 0.d0 .or. rho_f_read(i) <= 0.d0) &
        call stop_the_code('non-positive value of density')
      if (phi_read(i) <= 0.d0 .or. tortuosity_read(i) <= 0.d0) &
        call stop_the_code('non-positive value of porosity or tortuosity')
      if (kappa_s_read(i) <= 0.d0 .or. kappa_f_read(i) <= 0.d0 .or. kappa_fr_read(i) <= 0.d0 .or. mu_fr_read(i) <= 0.d0) &
        call stop_the_code('non-positive value of modulus')
      if (Qmu(i) <= 0.00000001d0) &
        call stop_the_code('non-positive value of Qmu')

    else if (icodemat(i) <= 0) then
      ! tomographic material
      number_of_materials_defined_by_tomo_file = number_of_materials_defined_by_tomo_file + 1

      if (number_of_materials_defined_by_tomo_file > 1) &
        call stop_the_code( &
'Just one material can be defined by a tomo file for now (we would need to write a nummaterial_velocity_file)')

      ! Assign dummy values for now (they will be read by the solver). Vs must be == 0 for acoustic media anyway
      rho_s_read(i) = -1.0d0
      cp(i) = -1.0d0
      cs(i) = val2read
      QKappa(i) = -1.0d0
      Qmu(i) = -1.0d0
      aniso3(i) = -1.0d0
      aniso4(i) = -1.0d0
      if (abs(cs(i)) > TINYVAL) then
        phi_read(i) = 0.d0           ! elastic
      else
        phi_read(i) = 1.d0           ! acoustic
      endif

    else
      ! default
      call stop_the_code('Unknown material code')

    endif

  enddo ! nbmodels

  ! user output
  write(IMAIN,*) 'Materials:'
  write(IMAIN,*) '  Nb of solid, fluid or porous materials = ',nbmodels
  write(IMAIN,*)

  do i = 1,nbmodels
     if (i == 1) write(IMAIN,*) '--------'
     if (icodemat(i) == ISOTROPIC_MATERIAL) then
        write(IMAIN,*) 'Material #',i,' isotropic'
        write(IMAIN,*) 'rho,cp,cs   = ',rho_s_read(i),cp(i),cs(i)
        write(IMAIN,*) 'Qkappa, Qmu = ',QKappa(i),Qmu(i)
        if (cs(i) < TINYVAL) then
           write(IMAIN,*) 'Material is fluid'
        else
           write(IMAIN,*) 'Material is solid'
        endif
     else if (icodemat(i) == ANISOTROPIC_MATERIAL) then
        write(IMAIN,*) 'Material #',i,' anisotropic'
        write(IMAIN,*) 'rho,cp,cs    = ',rho_s_read(i),cp(i),cs(i)
        if (AXISYM) then
          write(IMAIN,*) 'c11,c13,c15,c33,c35,c55,c12,c23,c25,c22 = ', &
                          aniso3(i),aniso4(i),aniso5(i),aniso6(i),aniso7(i),aniso8(i), &
                          aniso9(i),aniso10(i),aniso11(i),aniso12(i)
        else
          write(IMAIN,*) 'c11,c13,c15,c33,c35,c55,c12,c23,c25 = ', &
                          aniso3(i),aniso4(i),aniso5(i),aniso6(i),aniso7(i),aniso8(i), &
                          aniso9(i),aniso10(i),aniso11(i)
          write(IMAIN,*) 'QKappa,Qmu = ',QKappa(i),Qmu(i)
        endif
     else if (icodemat(i) == POROELASTIC_MATERIAL) then
        write(IMAIN,*) 'Material #',i,' isotropic'
        write(IMAIN,*) 'rho_s, kappa_s         = ',rho_s_read(i),kappa_s_read(i)
        write(IMAIN,*) 'rho_f, kappa_f, eta_f  = ',rho_f_read(i),kappa_f_read(i),eta_f_read(i)
        write(IMAIN,*) 'phi, tortuosity        = ',phi_read(i),tortuosity_read(i)
        write(IMAIN,*) 'permxx, permxz, permzz = ',permxx_read(i),permxz_read(i),permzz_read(i)
        write(IMAIN,*) 'kappa_fr, mu_fr, Qmu   = ',kappa_fr_read(i),mu_fr_read(i),Qmu(i)
        write(IMAIN,*) 'Material is porous'
     else if (icodemat(i) <= 0) then
        write(IMAIN,*) 'Material #',i,' will be read in an external tomography file (TOMOGRAPHY_FILE in Par_file)'
     else
        call stop_the_code('Unknown material code')
     endif
     write(IMAIN,*) '--------'
     call flush_IMAIN()
  enddo

  end subroutine read_material_table

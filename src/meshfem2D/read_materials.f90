
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

  subroutine read_materials(nb_materials,icodemat,cp,cs, &
                            aniso3,aniso4,aniso5,aniso6,aniso7,aniso8,aniso9,aniso10,aniso11, &
                            QKappa,Qmu,rho_s,rho_f,phi,tortuosity, &
                            permxx,permxz,permzz,kappa_s,kappa_f,kappa_fr, &
                            eta_f,mu_fr)

! reads in material definitions in DATA/Par_file

  implicit none
  include "constants.h"

  integer :: nb_materials

  integer, dimension(nb_materials) :: icodemat

  double precision, dimension(nb_materials) :: rho_s,cp,cs, &
    aniso3,aniso4,aniso5,aniso6,aniso7,aniso8,aniso9,aniso10,aniso11,QKappa,Qmu
  double precision, dimension(nb_materials) :: rho_f,phi,tortuosity,permxx,permxz,&
       permzz,kappa_s,kappa_f,kappa_fr,eta_f,mu_fr

  ! local parameters
  integer :: imaterial,i,icodematread
  double precision :: val0read,val1read,val2read,val3read,val4read, &
       val5read,val6read,val7read,val8read,val9read,val10read,val11read,val12read

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
  QKappa(:) = 9999.d0
  Qmu(:) = 9999.d0
  rho_s(:) = 0.d0
  rho_f(:) = 0.d0
  phi(:) = 0.d0
  tortuosity(:) = 0.d0
  permxx(:) = 0.d0
  permxz(:) = 0.d0
  permzz(:) = 0.d0
  kappa_s(:) = 0.d0
  kappa_f(:) = 0.d0
  kappa_fr(:) = 0.d0
  eta_f(:) = 0.d0
  mu_fr(:) = 0.d0

  ! reads in material parameters
  do imaterial=1,nb_materials

     call read_material_parameters_p(i,icodematread, &
                              val0read,val1read,val2read,val3read, &
                              val4read,val5read,val6read,val7read, &
                              val8read,val9read,val10read,val11read,val12read)

     ! checks material id
     if(i < 1 .or. i > nb_materials) stop 'Wrong material number!'
     icodemat(i) = icodematread


     ! sets material properties
     if(icodemat(i) == ISOTROPIC_MATERIAL) then

        ! isotropic materials

        rho_s(i) = val0read
        cp(i) = val1read
        cs(i) = val2read
        QKappa(i) = val5read
        Qmu(i) = val6read

! for Cs we use a less restrictive test because acoustic media have Cs exactly equal to zero
        if(rho_s(i) <= 0.00000001d0 .or. cp(i) <= 0.00000001d0 .or. cs(i) < 0.d0) &
            stop 'negative value of velocity or density'
        if(QKappa(i) <= 0.00000001d0 .or. Qmu(i) <= 0.00000001d0) stop 'non-positive value of QKappa or Qmu'

        aniso3(i) = val3read
        aniso4(i) = val4read
        if(cs(i) /= 0.d0) then
           phi(i) = 0.d0           ! elastic
        else
           phi(i) = 1.d0           ! acoustic
        endif
     else if (icodemat(i) == ANISOTROPIC_MATERIAL) then

        ! anisotropic materials

        rho_s(i) = val0read
        aniso3(i) = val1read
        aniso4(i) = val2read
        aniso5(i) = val3read
        aniso6(i) = val4read
        aniso7(i) = val5read
        aniso8(i) = val6read
        aniso9(i) = val7read
        aniso10(i) = val8read
        aniso11(i) = val9read

     else

        ! poroelastic materials

        rho_s(i) = val0read
        rho_f(i) = val1read
        phi(i) = val2read
        tortuosity(i) = val3read
        permxx(i) = val4read
        permxz(i) = val5read
        permzz(i) = val6read
        kappa_s(i) = val7read
        kappa_f(i) = val8read
        kappa_fr(i) = val9read
        eta_f(i) = val10read
        mu_fr(i) = val11read
        Qmu(i) = val12read

        if(rho_s(i) <= 0.d0 .or. rho_f(i) <= 0.d0) stop 'non-positive value of density'
        if(phi(i) <= 0.d0 .or. tortuosity(i) <= 0.d0) stop 'non-positive value of porosity or tortuosity'
        if(kappa_s(i) <= 0.d0 .or. kappa_f(i) <= 0.d0 .or. kappa_fr(i) <= 0.d0 .or. mu_fr(i) <= 0.d0) then
           stop 'non-positive value of modulus'
        endif
        if(Qmu(i) <= 0.00000001d0) stop 'non-positive value of Qmu'

     endif
  enddo

  ! user output
  print *
  print *, 'Nb of solid, fluid or porous materials = ',nb_materials
  print *
  do i=1,nb_materials
     if(icodemat(i) /= ANISOTROPIC_MATERIAL .and. icodemat(i) /= POROELASTIC_MATERIAL) then
        print *,'Material #',i,' isotropic'
        print *,'rho,cp,cs = ',rho_s(i),cp(i),cs(i),QKappa(i),Qmu(i)
        if(cs(i) < TINYVAL) then
           print *,'Material is fluid'
        else
           print *,'Material is solid'
        endif
     else if(icodemat(i) == POROELASTIC_MATERIAL) then
        print *,'Material #',i,' isotropic'
        print *,'rho_s, kappa_s= ',rho_s(i),kappa_s(i)
        print *,'rho_f, kappa_f, eta_f= ',rho_f(i),kappa_f(i),eta_f(i)
        print *,'phi, tortuosity, permxx, permxz, permzz= ',phi(i),tortuosity(i),permxx(i),permxz(i),permzz(i)
        print *,'kappa_fr, mu_fr, Qmu= ',kappa_fr(i),mu_fr(i),Qmu(i)
        print *,'Material is porous'
     else
        print *,'Material #',i,' anisotropic'
        print *,'rho,cp,cs = ',rho_s(i),cp(i),cs(i)
        print*,'c11,c13,c15,c33,c35,c55,c12,c23,c25 = ',aniso3(i),aniso4(i),aniso5(i),aniso6(i),aniso7(i),aniso8(i), &
                                                        aniso9(i),aniso10(i),aniso11(i)
        print *,'QKappa,Qmu = ',QKappa(i),Qmu(i)
     endif
     print *
  enddo

  end subroutine read_materials


!========================================================================
!
!                   S P E C F E M 2 D  Version 6.1
!                   ------------------------------
!
! Copyright Universite de Pau, CNRS and INRIA, France,
! and Princeton University / California Institute of Technology, USA.
! Contributors: Dimitri Komatitsch, dimitri DOT komatitsch aT univ-pau DOT fr
!               Nicolas Le Goff, nicolas DOT legoff aT univ-pau DOT fr
!               Roland Martin, roland DOT martin aT univ-pau DOT fr
!               Christina Morency, cmorency aT princeton DOT edu
!               Pieyre Le Loher, pieyre DOT le-loher aT inria.fr
!
! This software is a computer program whose purpose is to solve
! the two-dimensional viscoelastic anisotropic or poroelastic wave equation
! using a spectral-element method (SEM).
!
! This software is governed by the CeCILL license under French law and
! abiding by the rules of distribution of free software. You can use,
! modify and/or redistribute the software under the terms of the CeCILL
! license as circulated by CEA, CNRS and INRIA at the following URL
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


subroutine read_external_model(any_acoustic,any_elastic,any_poroelastic, &
                elastic,poroelastic,anisotropic,nspec,npoin,N_SLS,ibool, &
                f0_attenuation,inv_tau_sigma_nu1_sent,phi_nu1_sent, &
                inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu1_sent,Mu_nu2_sent, &
                inv_tau_sigma_nu1,inv_tau_sigma_nu2,phi_nu1,phi_nu2,Mu_nu1,Mu_nu2,&
                coord,kmato,myrank,rhoext,vpext,vsext, &
                Qp_attenuationext,Qs_attenuationext,c11ext,c13ext,c15ext,c33ext,c35ext,c55ext,READ_EXTERNAL_SEP_FILE)

  implicit none
  include "constants.h"

integer :: nspec,myrank,npoin
double precision  :: f0_attenuation

! Mesh
integer, dimension(NGLLX,NGLLZ,nspec) :: ibool
double precision, dimension(NDIM,npoin) :: coord

! Material properties
logical :: any_acoustic,any_elastic,any_poroelastic,READ_EXTERNAL_SEP_FILE
integer, dimension(nspec) :: kmato
logical, dimension(nspec) :: elastic,poroelastic
double precision, dimension(NGLLX,NGLLZ,nspec) :: rhoext,vpext,vsext

! for attenuation
integer :: N_SLS
double precision :: Mu_nu1_sent,Mu_nu2_sent
double precision, dimension(N_SLS) :: inv_tau_sigma_nu1_sent,phi_nu1_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent
double precision, dimension(NGLLX,NGLLZ,nspec,N_SLS) :: inv_tau_sigma_nu1,phi_nu1,inv_tau_sigma_nu2,phi_nu2
double precision, dimension(NGLLX,NGLLZ,nspec) :: Mu_nu1,Mu_nu2
double precision, dimension(NGLLX,NGLLZ,nspec) :: Qp_attenuationext,Qs_attenuationext

! for anisotropy
logical, dimension(nspec) :: anisotropic
double precision, dimension(NGLLX,NGLLZ,nspec) :: c11ext,c13ext,c15ext,c33ext,c35ext,c55ext

! Local variables
integer :: i,j,ispec,iglob
double precision :: previous_vsext
double precision :: tmp1, tmp2,tmp3

if(READ_EXTERNAL_SEP_FILE) then
        write(IOUT,*)
        write(IOUT,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
        write(IOUT,*) 'Assigning external velocity and density model (elastic and/or acoustic)...'
        write(IOUT,*) 'Read outside SEG model...'
        write(IOUT,*) '!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'

        open(unit=1001,file='DATA/model_velocity.dat_input',status='unknown')
        do ispec = 1,nspec
           do j = 1,NGLLZ
              do i = 1,NGLLX
                 iglob = ibool(i,j,ispec)
                 read(1001,*) tmp1,tmp2,tmp3,rhoext(i,j,ispec),vpext(i,j,ispec),vsext(i,j,ispec)
!     vsext(i,j,ispec)=0.0
              end do
           end do
        end do
        close(1001)

    else
       do ispec = 1,nspec
          do j = 1,NGLLZ
             do i = 1,NGLLX

                iglob = ibool(i,j,ispec)
                call define_external_model(coord(1,iglob),coord(2,iglob),kmato(ispec),myrank,&
                rhoext(i,j,ispec),vpext(i,j,ispec),vsext(i,j,ispec), &
                Qp_attenuationext(i,j,ispec),Qs_attenuationext(i,j,ispec),&
                c11ext(i,j,ispec),c13ext(i,j,ispec),c15ext(i,j,ispec),c33ext(i,j,ispec),c35ext(i,j,ispec),c55ext(i,j,ispec))
                if((c11ext(i,j,ispec) /= 0) .or. (c13ext(i,j,ispec) /= 0) .or. (c15ext(i,j,ispec) /= 0) .or. &
                       (c33ext(i,j,ispec) /= 0) .or. (c35ext(i,j,ispec) /= 0) .or. (c55ext(i,j,ispec) /= 0)) then
                   ! vp, vs : dummy values, trick to avoid floating point errors
                   vpext(i,j,ispec) = 20.d0
                   vsext(i,j,ispec) = 10.d0
                end if
             end do
          end do
       end do
    end if

      any_acoustic = .false.
      any_elastic = .false.
      any_poroelastic = .false.
      do ispec = 1,nspec
         previous_vsext = -1.d0
         do j = 1,NGLLZ
            do i = 1,NGLLX
               iglob = ibool(i,j,ispec)
               if(.not. (i == 1 .and. j == 1) .and. &
               ((vsext(i,j,ispec) >= TINYVAL .and. previous_vsext < TINYVAL) .or. &
               (vsext(i,j,ispec) < TINYVAL .and. previous_vsext >= TINYVAL)))  &
               call exit_MPI('external velocity model cannot be both fluid and solid inside the same spectral element')
               if((c11ext(i,j,ispec) /= 0) .or. (c13ext(i,j,ispec) /= 0) .or. (c15ext(i,j,ispec) /= 0) .or. &
                       (c33ext(i,j,ispec) /= 0) .or. (c35ext(i,j,ispec) /= 0) .or. (c55ext(i,j,ispec) /= 0)) then
                  anisotropic(ispec) = .true.
                  poroelastic(ispec) = .false.
                  elastic(ispec) = .true.
                  any_elastic = .true.
                  Qp_attenuationext(i,j,ispec) = 10.d0
                  Qs_attenuationext(i,j,ispec) = 10.d0
               elseif(vsext(i,j,ispec) < TINYVAL) then
                  elastic(ispec) = .false.
                  poroelastic(ispec) = .false.
                  any_acoustic = .true.
               else
                  poroelastic(ispec) = .false.
                  elastic(ispec) = .true.
                  any_elastic = .true.
               endif
               call attenuation_model(N_SLS,Qp_attenuationext(i,j,ispec),Qs_attenuationext(i,j,ispec), &
                    f0_attenuation,inv_tau_sigma_nu1_sent,phi_nu1_sent,inv_tau_sigma_nu2_sent,phi_nu2_sent,Mu_nu1_sent,Mu_nu2_sent)
               inv_tau_sigma_nu1(i,j,ispec,:) = inv_tau_sigma_nu1_sent(:)
               phi_nu1(i,j,ispec,:) = phi_nu1_sent(:)
               inv_tau_sigma_nu2(i,j,ispec,:) = inv_tau_sigma_nu2_sent(:)
               phi_nu2(i,j,ispec,:) = phi_nu2_sent(:)
               Mu_nu1(i,j,ispec) = Mu_nu1_sent
               Mu_nu2(i,j,ispec) = Mu_nu2_sent
               previous_vsext = vsext(i,j,ispec)
            enddo
         enddo
      enddo

end subroutine read_external_model

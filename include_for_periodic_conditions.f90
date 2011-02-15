
!----------------------------------------------------------------------
          do ispecperio2 = 1,NSPEC_PERIO

            ispec2 = numperio_right(ispecperio2)

            if(codeabs_perio_right(ILEFT,ispecperio2)) then
               i2 = 1
               do j2 = 1,NGLLZ
                  iglob2 = ibool(i2,j2,ispec2)
                  if(sqrt(abs(coord(2,iglob) - coord(2,iglob2))**2 + &
                     (abs(coord(1,iglob) - coord(1,iglob2)) - PERIODIC_horiz_dist)**2) < PERIODIC_DETECT_TOL) then
                    print *,iglob,' and ',iglob2,' are the same periodic point, merging them'
!                   print *,'horiz dist is = ',abs(coord(1,iglob) - coord(1,iglob2))
!                   print *,ispec,i,j,ispec2,i2,j2
!                   ibool(i2,j2,ispec2) = ibool(i,j,ispec)
!--------------------------------------------------------------------------------
                    iglob_target_to_replace = ibool(i2,j2,ispec2)
                    do ispec3 = 1,nspec
                      do j3 = 1,NGLLZ
                        do i3 = 1,NGLLX
                          if(ibool(i3,j3,ispec3) == iglob_target_to_replace) ibool(i3,j3,ispec3) = ibool(i,j,ispec)
                        enddo
                      enddo
                    enddo
!--------------------------------------------------------------------------------
                  endif
               enddo
            endif

            if(codeabs_perio_right(IRIGHT,ispecperio2)) then
               i2 = NGLLX
               do j2 = 1,NGLLZ
                  iglob2 = ibool(i2,j2,ispec2)
                  if(sqrt(abs(coord(2,iglob) - coord(2,iglob2))**2 + &
                     (abs(coord(1,iglob) - coord(1,iglob2)) - PERIODIC_horiz_dist)**2) < PERIODIC_DETECT_TOL) then
                    print *,iglob,' and ',iglob2,' are the same periodic point, merging them'
!                   print *,'horiz dist is = ',abs(coord(1,iglob) - coord(1,iglob2))
!                   print *,ispec,i,j,ispec2,i2,j2
!                   ibool(i2,j2,ispec2) = ibool(i,j,ispec)
!--------------------------------------------------------------------------------
                    iglob_target_to_replace = ibool(i2,j2,ispec2)
                    do ispec3 = 1,nspec
                      do j3 = 1,NGLLZ
                        do i3 = 1,NGLLX
                          if(ibool(i3,j3,ispec3) == iglob_target_to_replace) ibool(i3,j3,ispec3) = ibool(i,j,ispec)
                        enddo
                      enddo
                    enddo
!--------------------------------------------------------------------------------
                  endif
               enddo
            endif

            if(codeabs_perio_right(IBOTTOM,ispecperio2)) then
               j2 = 1
               do i2 = 1,NGLLX
                  iglob2 = ibool(i2,j2,ispec2)
                  if(sqrt(abs(coord(2,iglob) - coord(2,iglob2))**2 + &
                     (abs(coord(1,iglob) - coord(1,iglob2)) - PERIODIC_horiz_dist)**2) < PERIODIC_DETECT_TOL) then
                    print *,iglob,' and ',iglob2,' are the same periodic point, merging them'
!                   print *,'horiz dist is = ',abs(coord(1,iglob) - coord(1,iglob2))
!                   print *,ispec,i,j,ispec2,i2,j2
!                   ibool(i2,j2,ispec2) = ibool(i,j,ispec)
!--------------------------------------------------------------------------------
                    iglob_target_to_replace = ibool(i2,j2,ispec2)
                    do ispec3 = 1,nspec
                      do j3 = 1,NGLLZ
                        do i3 = 1,NGLLX
                          if(ibool(i3,j3,ispec3) == iglob_target_to_replace) ibool(i3,j3,ispec3) = ibool(i,j,ispec)
                        enddo
                      enddo
                    enddo
!--------------------------------------------------------------------------------
                  endif
               enddo
            endif

            if(codeabs_perio_right(ITOP,ispecperio2)) then
               j2 = NGLLZ
               do i2 = 1,NGLLX
                  iglob2 = ibool(i2,j2,ispec2)
                  if(sqrt(abs(coord(2,iglob) - coord(2,iglob2))**2 + &
                     (abs(coord(1,iglob) - coord(1,iglob2)) - PERIODIC_horiz_dist)**2) < PERIODIC_DETECT_TOL) then
                    print *,iglob,' and ',iglob2,' are the same periodic point, merging them'
!                   print *,'horiz dist is = ',abs(coord(1,iglob) - coord(1,iglob2))
!                   print *,ispec,i,j,ispec2,i2,j2
!                   ibool(i2,j2,ispec2) = ibool(i,j,ispec)
!--------------------------------------------------------------------------------
                    iglob_target_to_replace = ibool(i2,j2,ispec2)
                    do ispec3 = 1,nspec
                      do j3 = 1,NGLLZ
                        do i3 = 1,NGLLX
                          if(ibool(i3,j3,ispec3) == iglob_target_to_replace) ibool(i3,j3,ispec3) = ibool(i,j,ispec)
                        enddo
                      enddo
                    enddo
!--------------------------------------------------------------------------------
                  endif
               enddo
            endif

          enddo
!----------------------------------------------------------------------


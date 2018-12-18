
! extract thermocline from PNM image

! use "convert -compress none" to convert to ASCII PNM image
! and then "tr '\040' '\012'<thermocline.pnm > thermocline2.pnm" to suppress white spaces

  program extract_thermocline

  implicit none

  integer :: ix,iy,R,G,B,NX,NY

  logical :: first_red_pixel_in_this_line

  double precision :: cp,depth

! skip header
  read(*,*)
  read(*,*) NX,NY
  read(*,*)

! in the PNM format, the image starts in the upper-left corner
  do iy=NY,1,-1
    first_red_pixel_in_this_line = .true.
    do ix=1,NX
      read(*,*) R
      read(*,*) G
      read(*,*) B
! if pixel is red
      if ((R > 240 .and. G < 180 .and. B < 180 .and. first_red_pixel_in_this_line) .or. &
         (ix == NX .and. first_red_pixel_in_this_line)) then

        first_red_pixel_in_this_line = .false.

        cp = 1472.5d0 + (1495.d0 - 1472.5d0) * dble(ix-1)/ dble(NX-1)
        depth = -2000 + 2000 * dble(iy-1)/ dble(NY-1)

        if (ix < NX) then
          print *, depth,cp
        else
          print *, depth,'UNDEF'
        endif
      endif
    enddo
  enddo

 end program extract_thermocline



! permute color palette to have list of colors in random order

program permute_color_palette

implicit none

integer, parameter :: N = 236
character(len=50) nom(N)
double precision, dimension(N) :: r,g,b

integer, dimension(N) :: done,perm

integer i,irandom_color

real random_val

done(:) = -1

do i=1,N

777 continue
  call random_number(random_val)
  irandom_color = nint(random_val * (N+3)) - 1
  if(irandom_color < 1) irandom_color = 1
  if(irandom_color > N) irandom_color = N
if(done(irandom_color) /= -1) goto 777

perm(i) = irandom_color
done(irandom_color) = 100

enddo

!print *,'random done'
!do i=1,N
!write(*,*) perm(i)
!enddo

do i=1,N
read(*,*) nom(perm(i))
read(*,*) r(perm(i))
read(*,*) g(perm(i))
read(*,*) b(perm(i))
enddo

do i=1,N
write(*,*) '!#',nom(i)
write(*,*) '##red(',i,')#=#',r(i)
write(*,*) '##green(',i,')#=#',g(i)
write(*,*) '##blue(',i,')#=#',b(i)
write(*,*)
enddo

end program permute_color_palette


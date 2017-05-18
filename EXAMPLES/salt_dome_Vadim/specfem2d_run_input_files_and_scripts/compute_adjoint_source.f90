  program adjoit_source

    character(len=150) :: dir_data,   dir_synth, dir_adj
    character(len=150) :: file_data, file_synth, file_adj
    character(len=150) :: arg(3)
    character(len=150) :: name_in

    integer            :: i,nt
    double precision   :: t, x, z
    double precision, dimension(:), allocatable :: time, pressure_synth, pressure_data, adjoint_source

    do i = 1, 3
       call getarg(i, arg(i))
       write (*,*) arg(i)
    enddo
    dir_data=arg(1)
    dir_synth=arg(2)
    dir_adj=arg(3)


    open(10, file='list_stations_to_read.txt')

    do
       read(10,'(a)',end=99) name_in

       file_data  =  trim(dir_data)//trim(name_in)
       file_synth = trim(dir_synth)//trim(name_in)
       file_adj   = trim(dir_adj)//name_in(1:8)//'.POT.adj'

       open(11,file=trim(file_data))
       open(12,file=trim(file_synth))
       open(13,file=trim(file_adj))
       !count number of lines
       nt=0
       do
          read(11,*,end=97) t,z
          nt = nt + 1
       enddo
       97 close(11)

       allocate(time(nt), pressure_synth(nt), pressure_data(nt), adjoint_source(nt))
       open(11,file=trim(file_data))
       do ii=1,nt
          read(11,*) t,z
          read(12,*) t,x
          time(i)=t
          pressure_synth(i)=z
          pressure_data(i)=x
          !! need to compute 2nd time derivative of pressure residual

          !write(13,*)       t,z-x
       enddo
       98 continue
       close(11)
       close(12)

       adjoint_source(:)= pressure_synth(:)-pressure_data(:)
       dt_square = (time(2) - time(1))**2
       write(13,*) time(1),0.
       do it=2,nt-1
          write(13,*) time(it), adjoint_source(it-1) + adjoint_source(it+1) - 2.*adjoint_source(it) / dt_square
       enddo
       write(13,*) time(nt),0.
       close(13)
       deallocate(time, adjoint_source, pressure_synth, pressure_data)
    enddo
    99 close(10)


  end program adjoit_source

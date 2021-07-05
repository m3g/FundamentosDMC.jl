!
! Computes the g(r) up to distance = 10.
!

program radialdist

  use ff
  implicit none
  integer, parameter :: n = 100
  integer, parameter :: ndist = 100
  integer :: nin, i, j, idist, ioerr, nframes
  character(len=200) :: title
  character(len=3) :: element
  double precision :: x(n,2), gr(0:ndist), dist, density, xj, yj
  double precision, parameter :: pi = 4.d0*datan(1.d0)
  double precision :: image

  write(*,*) ' Program g(r) '

  ! Reseting histogram

  do i = 0, ndist
    gr(i) = 0.d0
  end do

  ! Reading coordinates, computing distances, and adding to histogram

  open(10,file='trajectory.xyz',action='read',status='old')
  nframes = 0
  do 

    read(10,*,iostat=ioerr) nin
    if ( ioerr /= 0 ) exit
    nframes = nframes + 1
    read(10,*) title
    do i = 1, n
      read(10,*) element, x(i,1), x(i,2)
    end do

     do i = 1, n-1
       do j = i+1, n
         xj = x(j,1)-x(i,1)
         yj = x(j,2)-x(i,2)
         xj = image(xj)
         yj = image(yj)
         dist = dsqrt(xj**2 + yj**2)
         idist = int(dist*10)
         if ( idist <= ndist ) gr(idist) = gr(idist) + 2.d0
       end do
     end do

  end do
  close(10)

  ! Normalizing by the number of frames and number of atoms

  do i = 0, ndist
    gr(i) = gr(i) / ( nframes*(n-1) )
  end do

  ! Normalizing by the number of atoms expected in an area of a hollow 
  ! disk of radius "dist" and width 0.1 

  density = n / side**2
  do i = 0, ndist
    dist = dble(i)/10.d0
    gr(i) = gr(i) / (density*pi*((dist+0.1)**2-dist**2))
  end do
 
  write(*,*) ' Number of frames read: ', nframes

  ! Print results

  open(10,file='gr.dat')
  do i = 0, ndist
    dist = dble(i)/10.d0 + 0.05
    write(10,*) dist, gr(i)
  end do
  close(10)

  write(*,*) ' Wrote file: gr.dat '

end program radialdist















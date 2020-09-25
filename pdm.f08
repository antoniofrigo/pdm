program pdm
  use stellingwerf
  implicit none
  character(len=100)             :: filename ! Filename
  integer                        :: n_lines = 0 ! Lines in file
  real, dimension(:),allocatable :: x_values ! Time values
  real, dimension(:),allocatable :: y_values ! Magnitudes
  integer                        :: io

  integer :: i ! dummy variable
  
  real :: start, finish ! CPU Time
  filename = "OGLE-BLG-CEP-027.dat"

  ! Finds length of file
  open(unit = 1, file = filename)
  do
    read(1,*,iostat=io)
    if (io.ne.0) exit
    n_lines = n_lines + 1
  end do
  close(1)

  ! Reads in the data
  open(unit = 1, file = filename)
  allocate(x_values(n_lines))
  allocate(y_values(n_lines))
  do i = 1, n_lines
    read(1,*) x_values(i), y_values(i)
  end do
  close(1)

  call cpu_time(start)
  call stellingwerf_period(x_values, y_values)
  call cpu_time(finish)
  write(*,*) "Total Time:", finish - start
  stop 0
end program pdm

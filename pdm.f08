program pdm
  use stellingwerf, only: find_segments, largest_segment, &
    welford_modified, pooled_variance, variance, output_results, &
    pooled_variance_seg
  implicit none
  character(len=100) :: data_file, output_file ! Filenames
  character(len=100) :: item ! Commandline arguments
  character(len=1)   :: creturn ! Carriage return
  integer            :: LUN1, LUN2

  ! Number of bins, separation of segments in stdevs, maximum frequency
  integer            :: num_bins, stdev_sep, max_freq

  integer :: io, k, i, j ! Dummy variables
  integer :: num_line = 0 ! Number of lines in data file
  integer :: num_segs, total  ! Number of segments, total trials
  real    :: start, finish ! Runtime variables

  real(kind=kind(1.0d0)), dimension(:,:), allocatable :: observs ! Observations
  real(kind=kind(1.0d0)), dimension(:,:), allocatable :: theta ! (period, theta)
  real(kind=kind(1.0d0)), dimension(:), allocatable   :: x_values ! Time obs.
  real(kind=kind(1.0d0)), dimension(:), allocatable   :: y_values ! Mag. obs.

  ! Arrays for storing variance information
  ! var_array (mean, mean_old, msq, num_points)
  ! var_array_seg(num_points, var_p)
  real(kind=kind(1.0d0)), dimension(:,:), allocatable :: var_array
  real(kind=kind(1.0d0)), dimension(:,:), allocatable :: var_array_seg

  integer, dimension(:), allocatable :: segments ! Indices of segments

  ! Largest segment (in time), frequency step, period
  real(kind=kind(1.0d0)) :: max_seg, freq_step, period
  ! Variance of entire data set, pooled variance, first scan's best period
  real(kind=kind(1.0d0)) :: var_o, var_p, first_period

  ! Commandline arguments in the form 
  ! data_file num_bins stdev_sep max_freq
  call get_command_argument(1,item)
  read(item, *) data_file
  call get_command_argument(2,item)
  read(item, *) num_bins
  call get_command_argument(3,item)
  read(item, *) stdev_sep
  call get_command_argument(4,item)
  read(item, *) max_freq

  ! Reading in data to observ
  output_file = "output_theta.dat"
  open(newunit=LUN1, file=data_file,status='OLD',action='READ',iostat=io)
  do 
    read(LUN1,*,iostat=io) 
    if (io.ne.0) exit
    num_line = num_line + 1
  end do
  close(LUN1)

  allocate(observs(2,num_line))

  open(newunit=LUN1, file=data_file,status='OLD',action='READ',iostat=io)
  do i = 1, num_line
    read(LUN1,*,iostat=io) observs(1,i), observs(2,i)
  end do
  close(LUN1)

  ! Allocates all the arrays
  allocate(x_values(num_line))
  allocate(y_values(num_line))
  allocate(var_array(4,num_bins))
  allocate(var_array_seg(4,num_bins))
  var_array = 0.0d0
  var_array_seg = 0.0d0

  ! FULL PDM PROCEDURE
  call cpu_time(start)
  ! Finds largest segment
  call find_segments(observs, segments, num_line, stdev_sep)
  call largest_segment(max_seg, segments)

  ! Finds total trials and allocates memory accordingly
  total = ceiling(2.0d0 * max_seg * max_freq)
  allocate(theta(2,total+2000))
  theta = 1.0d0
  var_o = variance(observs(2,:))

  x_values = observs(1,:)
  y_values = observs(2,:)
  num_segs = size(segments) - 1

  creturn = achar(13)

  ! The primary scan locates the approximate location of the minima, while
  ! the secondary scan checks the neighborhood for any other smaller theta.
  ! This means that we can still have high precision, so long as the primary
  ! scan gets sufficiently close to the true minima.

  ! To properly utilize the segments, we find the variance for each bin 
  ! within a given segment which are stored in var_array. These are then pooled
  ! together for each segment and are stored in var_array_seg. We can then pool
  ! these together again and find the pooled variance for the entire set.
  ! This allows us to use substantially smaller frequency steps.

  ! We use the method to compute the overall variance of the data set,
  ! except that we do not take time modulo the period and we do not have 
  ! any bins. This is necessary to keep the s^2/sigma^2 ratio unity.
  do j = 1, num_segs ! Segment loop
    write(*,"(A,A5,I15,A4,I15)", advance='no') creturn,"Trial",i," of ",total
    var_array_seg(1,j) = segments(j+1) - segments(j) ! Num. of pts in seg.
    var_array_seg(2,j) = variance(y_values(segments(j):segments(j+1))) ! Var. in segment
  end do
  var_o = 0.0d0
  call pooled_variance_seg(var_o, var_array_seg) ! Overall variance of the set

  ! PRIMARY SCAN
  do i = 1, total ! Time Loop
    do j = 1, num_segs ! Segment loop
      write(*,"(A,A5,I15,A4,I15)", advance='no') creturn,"Trial",i," of ",total
      period = 2.0d0 * max_seg / (1.0d0*i)
      var_array = 0.0d0
      var_p = 0.0d0

      ! Finds the variance and stores in in var_array within seg.
      call welford_modified(var_array, &
        mod(x_values(segments(j):segments(j+1)),period)/period,&
        y_values(segments(j):segments(j+1)),&
        num_bins) 
      ! Finds pooled var. of all bins
      call pooled_variance(var_p, var_array) 

      var_array_seg(1,j) = segments(j+1) - segments(j) ! Num. of pts in seg.
      var_array_seg(2,j) = var_p ! Pooled var. in segment
    end do
    ! Determines Theta for each time step
    var_p = 0.0d0
    call pooled_variance_seg(var_p, var_array_seg)
    theta(1,i) = period
    theta(2,i) = var_p/var_o
  end do
  write(*,"(A)") creturn

  ! Current best-estimated period
  first_period = theta(1,minloc(theta(2,:),1))

  ! SECONDARY SCAN (refines the estimate of the period)
  do i = 1, 2000 ! Time Loop
    var_array_seg = 0.0d0
    do j = 1, num_segs ! Segment loop
      write(*,"(A,A5,I15,A4,I15)", advance='no') creturn,"Trial",i," of ",2000
      period = 1/(1/first_period - 1/(1.0d0*max_seg) + i/(1000.d0*max_seg))
      var_array = 0.0d0
      var_p = 0.0d0

      ! Finds the variance and stores in in var_array within seg.
      call welford_modified(var_array, &
        mod(x_values(segments(j):segments(j+1)),period)/period,&
        y_values(segments(j):segments(j+1)),&
        num_bins) 
      ! Finds pooled var. of all bins
      call pooled_variance(var_p, var_array) 

      var_array_seg(1,j) = segments(j+1) - segments(j) ! Num. of pts in seg.
      var_array_seg(2,j) = var_p ! Pooled var. in segment
    end do
    ! Determines Theta for each time step
    var_p = 0.0d0
    call pooled_variance_seg(var_p, var_array_seg)
    theta(1,total + i) = period
    theta(2,total + i) = var_p/var_o
  end do
  write(*,"(A)") creturn

  ! Writes output to stdout and to file
  call cpu_time(finish)
  call output_results(theta, observs,output_file, total)
  write(*,*) "Time elapsed: ", finish - start
  deallocate(observs, theta, segments, var_array, x_values, y_values)
  stop 0
end program pdm

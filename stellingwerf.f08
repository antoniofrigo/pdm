module stellingwerf
  contains
    subroutine stellingwerf_period(x_values, y_values,m,max_freq)
      ! computes the best period using stellingwerf's algorithm
      implicit none
      real, dimension(:) :: x_values, y_values ! input points
      character(len=100) :: output_form
      integer :: total

      integer, optional, value :: m ! num bins
      real, optional, value    :: max_freq ! maximum_frequency

      real, allocatable, dimension(:,:) :: var_array
      real, allocatable, dimension(:,:) :: sig_array
      ! mean, mean_old, variance * (count-1), count
      ! stores information of array_y in bins

      real :: var_o ! variance of the original array
      real :: mean ! mean of the original array
      real :: var_p ! variance of the bins (pooled)
      real :: period ! guessed period
      real :: frequency ! guessed frequency
      real :: interval = 0! period interval, assumes x_values sorted
      real :: msq_sum ! sum of msqs across all
      real :: n_sum ! sum of sum(n) - sizeof(n) across all

      integer :: i,j  ! dummy variable
      integer :: s = 1! dummy variable
      integer :: test_interval ! dummy variable

      ! segmentation 
      integer, dimension(10000)          :: seg_index_o ! indices of breakpoints (with 0s)
      integer, allocatable, dimension(:) :: seg_index ! indices of breakpoints
      real, dimension(:), allocatable    :: diff ! difference array
      real                               :: mean_d, msq_d ! mean difference, msq of difference
      integer                            :: counter = 0 ! number of breakpoints

      if(.not. present(m)) m = 5
      if(.not. present(max_freq)) max_freq = 20

      allocate(diff(size(x_values)-1))
      seg_index_o = 0

      do i = 1, (size(x_values) - 1)
        diff(i) = x_values(i+1) - x_values(i)
      end do
      
      call welford(diff, msq_d, mean_d)
      do i = 1, size(diff)-1, 1
        if (diff(i).ge.(40)) then
          counter = counter + 1
          seg_index_o(counter) = i
        end if
      end do

      ! Set of break points [)
      allocate(seg_index(count(seg_index_o/=0)))
      seg_index = pack(seg_index_o, seg_index_o /=0)

      do i = 1, size(seg_index) - 1, 1
        if (i.eq.1) then
          test_interval = x_values(seg_index(1)-1) - x_values(1)
        else if (i.eq.(size(seg_index) - 1)) then
          test_interval = x_values(size(x_values)) - x_values(seg_index(i))
        else
          test_interval = x_values(seg_index(i+1)-1) - x_values(seg_index(i))
        end if
        if (test_interval.ge.interval) then
          interval = test_interval
        end if
      end do

      total = ceiling(2*interval*max_freq)
      ! End segmentation

      allocate(var_array(4,m * size(seg_index)))
      allocate(sig_array(total,2))
      sig_array = 1.0

      call welford(y_values, var_o, mean)

      write(*,*) "Total cases:", total
      do i = 1, total ! period loop
        do j = 1, size(seg_index)-1, 1 ! seg loop
          var_array = 0.0
          var_p = 0.0
          frequency = 1/(2 * interval) * i
          period = 1 / frequency
          call welford_modified(var_array, &
            mod(x_values(seg_index(j):seg_index(j+1)),period) / period, &
            y_values(seg_index(j):seg_index(j+1)), m)
        end do
        call pooled_variance(var_p, var_array(3,:), var_array(4,:))
        sig_array(i,1) = var_p / var_o
        sig_array(i,2) = period
      end do

      output_form = "(I3, F12.6, F12.6)"

      write(*,"(A3, A12, A12)") "\n#", "Theta", "Period (d)"
      s=minloc(sig_array(:,1), 1)
      write(*,output_form) 1, sig_array(s,1), sig_array(s,2)
      sig_array(s,1) = 1.0
      s=minloc(sig_array(:,1), 1)
      write(*,output_form)2, sig_array(s,1), sig_array(s,2)
      sig_array(s,1) = 1.0
      s=minloc(sig_array(:,1), 1)
      write(*,output_form) 3, sig_array(s,1), sig_array(s,2)
    end subroutine

    subroutine welford_modified(var_array, array_x, array_y, m)
      ! modified version of welford's online
      ! algorithm to compute the variance for each bin
      ! expects input of two arrays
      implicit none
      real, dimension(:) :: array_x
      real, dimension(:) :: array_y
      integer            :: m

      integer :: i ! dummy variable
      integer :: j ! dummy variable

      real, dimension(:,:) :: var_array
      ! mean, mean_old, variance * (count-1), count
      ! stores information of array_y in bins

      do i = 1,size(array_y)
        j = int(ceiling(m * array_x(i)))
        if (j.eq.0) then
          j = 1
        else if (j.eq.6) then
          j = 5
        end if
        var_array(4,j) = var_array(4,j) + 1 ! Ups count
        var_array(2,j) = var_array(1,j) ! Sets mean_old = mean
        
        ! mean = ((count - 1) * mean_old + y)/count
        var_array(1,j) = (var_array(4,j) - 1) * var_array(2,j) + array_y(i)
        var_array(1,j) = var_array(1,j) / var_array(4,j)

        ! msq = msq + (y - mean)(y - mean_old)
        var_array(3,j) = var_array(3,j) + &
          (array_y(i) - var_array(1,j)) * (array_y(i) - var_array(2,j))
      end do
    end subroutine

    subroutine welford(array, var_o, mean)
      ! computes variance of an array using welford's method
      ! in the standard way
      implicit none
      real, dimension(:) :: array
      integer            :: i
      real               :: var_o

      real :: mean, mean_old
      real :: msq

      do i = 1, size(array)
        mean_old = mean
        mean = (( i - 1.0 ) * mean_old + array(i))/i
        msq = msq + (array(i) - mean) * (array(i) - mean_old)
      end do
      var_o = msq / (size(array) - 1.0)
    end subroutine

    subroutine pooled_variance(var_p, msq_array, counts)
      ! computes pooled variance of all the segments
      implicit none
      real, dimension(:) :: msq_array ! array of (x- mean)^2 values per bin
      real, dimension(:) :: counts ! counts per bin
      real               :: var_p

      var_p = real(sum(msq_array))/real((sum(counts) - size(counts)))
    end subroutine
end module

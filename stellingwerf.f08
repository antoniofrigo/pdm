module stellingwerf
  contains
    subroutine find_segments(observs, segments, num, stdevs)
      ! Finds the locations of the segments
      implicit none
      real(kind=kind(1.0d0)),dimension(:,:), intent(in) :: observs
      integer,dimension(:), allocatable, intent(inout)  :: segments

      integer, intent(in) :: num, stdevs

      real(kind=kind(1.0d0)), dimension(:), allocatable :: diffs

      real(kind=kind(1.0d0)) :: disp
      integer                :: i, counter=1
      allocate(diffs(num-1))

      ! Creates finite difference array
      do i = 1, num-1
        diffs(i) = observs(1,i+1) - observs(1,i)
      end do

      ! Finds number of start indices
      disp = mean(diffs) + stdevs * sqrt(variance(diffs))
      do i = 1, num -1
        if (diffs(i).ge.disp) then
          counter = counter + 1
        end if
      end do

      ! Allocates memory for segment indices and sets them
      allocate(segments(counter+1))
      segments = 0
      counter = 2
      do i = 1, num-1
        if (diffs(i).ge.disp) then
          segments(counter) = i+1
          counter = counter + 1
        end if
      end do
      segments(1) = 1
      segments(counter) = num + 1

      deallocate(diffs)
    end subroutine

    subroutine largest_segment(max_seg, segments)
      ! Finds largest segment
      implicit none
      real(kind=kind(1.0d0)), intent(inout) :: max_seg
      integer,dimension(:), intent(in)   :: segments

      integer :: i, max_size

      max_seg = 0
      max_size = size(segments) - 1

      do i = 1, max_size
        if ((segments(i+1) - segments(i)).gt.max_seg) then
          max_seg = segments(i+1) - segments(i)
        end if
      end do
    end subroutine

    real(kind=kind(1.0d0)) function mean(array)
      implicit none
      real(kind=kind(1.0d0)), dimension(:), intent(in) :: array

      mean = 0.0d0
      mean = sum(array)/(size(array)*1.0d0)
      return
    end function

    real(kind=kind(1.0d0)) function variance(array)
      ! Standard Welford implementation
      implicit none
      real(kind=kind(1.0d0)), dimension(:), intent(in) :: array

      real(kind=kind(1.0d0)) :: msq, mu, mu_old

      integer :: n, i
      msq = 0.0d0
      mu = 0.0d0
      mu_old = 0.0d0
      n = size(array)

      do i = 1, n
        mu_old = mu
        mu = mu + array(i)
        if (i.eq.1) cycle
        msq = msq + (array(i) - mu_old/(i-1.0d0)) * (array(i) - mu/(1.0d0*i))
      end do

      variance = msq/( n - 1.0d0)
    end function

    subroutine welford_modified(var_array, x_val, y_val, num_bins)
      ! Calcualtes 'variance' for each bin
      implicit none
      integer                                :: num_bins
      real(kind=kind(1.0d0)), dimension(:)   :: x_val
      real(kind=kind(1.0d0)), dimension(:)   :: y_val
      real(kind=kind(1.0d0)), dimension(:,:) :: var_array

      integer :: i, j

      do i = 1,size(x_val)
        j = int(ceiling(num_bins * x_val(i)))
        if (j.eq.0) then
          j = 1
        else if (j.eq.(num_bins+1)) then
          j = num_bins
        end if
        var_array(4,j) = var_array(4,j) + 1 ! Ups count
        var_array(2,j) = var_array(1,j) ! Sets mean_old = mean
        
        var_array(1,j) = var_array(2,j) + y_val(i)

        ! msq = msq + (y - mean)(y - mean_old)
        if (var_array(4,j).eq.1.0) then
          var_array(3,j) = 0.0
        else
          var_array(3,j) = var_array(3,j) + &
            (y_val(i) - var_array(1,j)/var_array(4,j)) * &
            (y_val(i) - var_array(2,j)/(var_array(4,j) - 1))
        end if
      end do
    end subroutine

    subroutine pooled_variance(var_p, var_array)
      ! Calculates pooled variance
      implicit none
      real(kind=kind(1.0d0)), dimension(:,:), intent(in) :: var_array
      real(kind=kind(1.0d0)), intent (inout)             :: var_p
      var_p = sum(var_array(3,:))/(sum(var_array(4,:)) - size(var_array(4,:)))
    end subroutine

    subroutine pooled_variance_seg(var_p, var_array_seg)
      ! Calculates pooled variance
      implicit none
      real(kind=kind(1.0d0)), dimension(:,:), intent(in) :: var_array_seg
      real(kind=kind(1.0d0)), intent (inout)             :: var_p
      var_p = sum(var_array_seg(2,:))/(sum(var_array_seg(1,:)) - &
        size(var_array_seg(1,:)))
    end subroutine


    subroutine output_results(theta, observs, output_file, total)
      ! Outputs results to console and to file
      implicit none
      character(len=100)                                 :: output_file
      character(len=100)                                 :: fmt1
      integer, intent(in)                                :: total
      real(kind=kind(1.0d0)), dimension(:,:), intent(in) :: theta
      real(kind=kind(1.0d0)), dimension(:,:), intent(in) :: observs

      integer :: LUN, i
      integer :: k1, k2, k3

      fmt1 = "(I3, F15.10, F15.10)"
      write(*,"(A3, A15, A15)") "  #", "Period (d)", "Theta"
      k1 = minloc(theta(2,:),1)
      write(*,fmt1) 1, theta(1,k1), theta(2,k1)

      open(newunit=LUN, file = output_file)
      do i = 1, size(theta(1,:))
        write(LUN,*) 1.0d0/theta(1,i), theta(2,i)
      end do
      close(unit=LUN)
    end subroutine

end module

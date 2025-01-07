!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function factorial(n) result(res)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    integer, intent(in) :: n
    real*8 :: res
    integer :: i

    ! Initialize result to 1 for multiplication
    res = 1.0d0
    do i = n, 1, -1
        res = res * real(i, kind=8)
    end do
end function factorial

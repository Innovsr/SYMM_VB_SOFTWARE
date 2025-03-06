!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function factorial(n) result(res)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    implicit none
    integer, intent(in) :: n
    integer :: i, res

    ! Initialize result to 1 for multiplication
    res = 1
    do i = n, 1, -1
        res = res * i
    end do
end function factorial

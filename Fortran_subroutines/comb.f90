
!!!! function for calculation combination set (mcn) 
function comb(m, n)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
implicit none
integer:: m, n, comb, factorial

comb = factorial(m)/(factorial(m - n) * factorial(n))

return
end function comb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

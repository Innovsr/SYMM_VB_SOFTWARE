
!!!! function for calculation combination set (mcn) 
function comb(m, n)
integer:: m, n, comb, factorial

comb = factorial(m)/(factorial(m - n) * factorial(n))

return
end function comb


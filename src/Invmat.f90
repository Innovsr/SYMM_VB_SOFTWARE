!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! this subroutine do the matrix inversion required for diagonalisation !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      SUBROUTINE Invmat(Ndim,Ifail)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
use commondat
use commondat1
      Implicit DOUBLE PRECISION(A-H,O-Z)
      Dimension A(1000,1000),B(1000,1000)
      integer::ndim,ifail



!print*,'enter Invmat'
do i=1,ndim
  do j=1,ndim
    A(i,j)=ind_mat(i,j)
  enddo
enddo

111 format(a,30I5)

      N=Ndim
      Ifail = 0
      Xdet = 1.D0
      Do I = 1,N
       Do J = 1,N
        B(J,I) = 0.D0
       Enddo
       B(I,I) = 1.D0
      Enddo

      Do I = 1,N-1
!--Find the Max
       Amax = Dabs(A(I,I))
       Nmax = I

       Do J = I+1,N
        Aij = Dabs(A(I,J))
        If(Aij.Gt.Amax) Then
         Amax = Aij
         Nmax = J
        Endif
       Enddo

       If(Amax.Lt.1.D-11) Then
        Ifail = 1
        Xdet = 0.D0
        Return
       Endif

       If(Nmax.Ne.I) Then
        Xdet = -Xdet
        Do J = 1,N
         Atmp = A(J,I)
         A(J,I) = A(J,Nmax)
         A(J,Nmax) = Atmp
         Atmp = B(J,I)
         B(J,I) = B(J,Nmax)
         B(J,Nmax) = Atmp
        Enddo
       Endif

       Do J = I+1,N
        Atmp = A(I,J)/A(I,I)
        Do K = 1,N
         A(K,J) = A(K,J)-Atmp*A(K,I)
         B(K,J) = B(K,J)-Atmp*B(K,I)
        Enddo
       Enddo
        Do K = 1,N
        Enddo
      Enddo

      If(Dabs(A(N,N)).Lt.1.D-11) Then
       Ifail = 1
       Xdet = 0.D0
       Return
      Endif

      Do I = 1,N
       Xdet = Xdet*A(I,I)
      Enddo

!-- Vanish Upper
      Do I = N,2,-1
       Do J = I-1,1,-1
        Atmp = A(I,J)/A(I,I)
        Do K = 1,N
         A(K,J) = A(K,J)-Atmp*A(K,I)
         B(K,J) = B(K,J)-Atmp*B(K,I)
        Enddo
       Enddo
      Enddo

!-- Normalized
      Do I = 1,N
       Atmp = A(I,I)
       Do J = 1,N
        B(J,I) = B(J,I)/Atmp
       Enddo
      Enddo
      End


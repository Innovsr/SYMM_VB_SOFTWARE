!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Calculate transpose of matrix
!-------------------------------
      SUBROUTINE MatTran(A,N)
!-------------------------------
!
!
!
        implicit none
        integer N,I,J
        Double Precision A,x
        Dimension A(N,N)

        Do I=1,N-1
          Do J=I+1,N
            x=A(J,I)
            A(J,I)=A(I,J)
            A(I,J)=x
          Enddo
        Enddo

        Return
      End


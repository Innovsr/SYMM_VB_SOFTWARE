!--------------------------------------------
      SUBROUTINE MatLDR(name,col,N,D)
!--------------------------------------------
!
!     A = U * D * VT**T
!     by two-sides Jacobi's method
!     U,VT::: left and right eigen value of A
!     right eigenvectors are stored in transpose form
!     D::: diagonal matrix stored in one dimension array
!     N::: the length of row/column of matrix A
!     Nlt::: is the Nulity of matrix A, which is equal to the number of zero eigenvalues 
!     Up::: upper main tridiagonal
!     Low::: lowwer main tridiagonal
!     Amax::: maximum off-diagonal matrix elements
!
      USE commondat1
      Implicit none
      Double Precision A(N,N),U(N,N),D(N),VT(N,N),bond_vec_mat(150,150)
      Double Precision Amax,Anzo,Aii,Ajj,Aij,Aji,Aod,Adif,tmp
      integer N,Nlt,col(1000)
      integer I,J,K,L,Iplus,faiil
      integer iter,Itmax,Dotri
      Double Precision Zer,One,Eps,ovlp
      INTEGER IW,IER
      character(len=3)::name
common/fail/faiil
common /bnd_mat/bond_vec_mat

      Data Zer,One,Eps/0.D0,1.D0,1.0D-20/


if(name.eq.'str')then
  do i=1,N
    do j=1,N
      A(i,j)=ovlp_mat_norm(col(i),col(j))
    enddo
  enddo
endif

if(name.eq.'ind')then
  do i=1,N
    do j=1,N
      A(i,j)=ind_mat(i,j)
    enddo
  enddo
endif

if(name.eq.'bnv')then
  do i=1,N
    do j=1,N
      A(i,j)=bond_vec_mat(i,j)
    enddo
  enddo
endif

      do i = 1, N
        D(i) = Zer
        do j = 1, N
          U(j,i) = Zer
          VT(j,i) = Zer
        enddo
        U(i,i) = One
        VT(i,i) = One
      enddo
      Nlt = 0
      if(N.le.0) Return
      if(N.eq.1) then
          D(1) = A(1,1)
          if(dabs(D(1)) .lt. 1.d-12 ) Nlt = 1
          Return
      endif
      itmax = 200
      iter = 0
      Amax = One
      do while(Amax .gt. Eps)
          iter = iter + 1
          Amax = Zer
          do i = 2, N
          do j = 1, i - 1
            Aii = A(i,i)
            Ajj = A(j,j)
            Aij = A(i,j)
            Aji = A(j,i)
            if(dabs(Aij) .gt. Amax) Amax = dabs(Aij)
            if(dabs(Aji) .gt. Amax) Amax = dabs(Aji)
            if(dabs(Aij) .lt. Eps .and. dabs(Aji) .lt. Eps) cycle
            Call TwoSideJacobi(i,j,Aii,Ajj,Aij,Aji,U,VT,A,N)
          end do
          end do
      end do

      do i = 1,N
        D(i) = A(i,i)
      End do
      if((iter / N) .eq. itmax) then
        write(*,*)'fail in MatLDR'
        write(*,*)'Average number of sweep:::',iter / N
        Stop
      endif

      Do I = 1, N-1
        tmp = dabs(D(I))
        K = I
        Do J = I+1, N
          If(dabs(D(J)) .GT. tmp) then
            tmp = dabs(D(J))
            K = J
          Endif
        Enddo
        If( K .NE. I) Then
          tmp = D(I)
          D(I) = D(K)
          D(K) = tmp
          Do J = 1, N
            tmp = U(J,I)
            U(J,I) = U(J,K)
            U(J,K) = tmp
            tmp = VT(J,I)
            VT(J,I) = VT(J,K)
            VT(J,K) = tmp
          Enddo
        Endif
      Enddo
      Do I=N,1,-1
        If(Dabs(D(I)).LT.1.d-12) Then
          Nlt=Nlt+1
        Else  
          Exit
        Endif
      Enddo
      Call MatTran(VT,N)

! diadonal elementsare stored in D

      ovlp=1.0
      do i=1,N
        ovlp=ovlp*D(i)
      enddo
      Return
      End

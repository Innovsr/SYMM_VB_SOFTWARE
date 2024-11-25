
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! The two-sided Jacobi method iteratively reduces off-diagonal elements of A to zero while maintaining its symmetry. !
!!The eigenvalues accumulate on the diagonal of A, and the eigenvectors are built in U and VT.                        !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!----------------------------------------------------------------
      SUBROUTINE TwoSideJacobi(i,j,Aii,Ajj,Aij,Aji,U,VT,A,N)
!----------------------------------------------------------------
!
!     This SUBROUTINE calculate two sides Jacobi rotation to following
!     matrix
!     | Aii Aij |       | Aii~  0  |
!     |         |  -->  |          |
!     | Aji Ajj |       |  0  Ajj~ |
!
      Implicit none
      Double Precision A(N,N),U(N,N),D(N),VT(N,N)
      Double Precision t1,c1,s1
      Double Precision t2,c2,s2
      Double Precision tp,cp,sp
      Double Precision tm,cm,sm
      Double Precision Aii,Ajj,Aij,Aji,tt,Aod,Adif,tmp
      integer N
      integer I,J,k
      Double Precision Zer,One,Eps
      Data Zer,One,Eps/0.D0,1.D0,1.0D-20/
      Adif = Aii - Ajj
      Aod = Aij + Aji
      if(dabs(Adif) .lt. dabs(Aod)) then
        tt = Adif / Aod
        tp = dsqrt(One + tt*tt)
        if(tt .lt. Zer) tp = -tp
        tp = tp - tt
      else if (dabs(Adif).lt.Eps) then
        tp = Zer
      else
        tt = Aod / Adif
        if(dabs(tt) .gt. 1.d-7) then
          tp = (dsqrt(One + tt*tt) - One) / tt
        else
          tp = 0.5d0*tt - 0.125d0*tt*tt + 0.0625d0*tt*tt*tt 
        endif
      endif
      cp = One / dsqrt(One + tp*tp)
      sp = tp * cp
      Adif = Aii + Ajj
      Aod = Aji - Aij
      if(dabs(Adif) .lt. dabs(Aod)) then
        tt = Adif / Aod
        tm = dsqrt(One + tt*tt)
        if(tt .lt. Zer) tm = -tm
        tm = tm - tt
      else if (dabs(Adif).lt.Eps) then
        tm = Zer
      else
        tt = Aod / Adif
        if(dabs(tt) .gt. 1.d-7) then
          tm = (dsqrt(One + tt*tt) - One) / tt
        else
          tm = 0.5d0*tt - 0.125d0*tt*tt + 0.0625d0*tt*tt*tt 
        endif
      endif
      cm = One / dsqrt(One + tm*tm)
      sm = tm * cm
      c1 = cp*cm - sp*sm
      s1 = sp*cm + cp*sm
      c2 = cp*cm + sp*sm
      s2 = sp*cm - cp*sm
      Do k = 1, N
        tmp=c1*U(k,j)-s1*U(k,i)
        U(k,i)=s1*U(k,j)+c1*U(k,i)
        U(k,j)=tmp
        tmp=c1*A(j,k)-s1*A(i,k)
        A(i,k)=s1*A(j,k)+c1*A(i,k)
        A(j,k)=tmp
      enddo
!     update VT and A
      Do k = 1, N
        tmp=c2*VT(k,j)-s2*VT(k,i)
        VT(k,i)=s2*VT(k,j)+c2*VT(k,i)
        VT(k,j)=tmp
        tmp=c2*A(k,j)-s2*A(k,i)
        A(k,i)=s2*A(k,j)+c2*A(k,i)
        A(k,j)=tmp
      enddo

      A(j,i) = Zer
      A(i,j) = Zer

      return
      End
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

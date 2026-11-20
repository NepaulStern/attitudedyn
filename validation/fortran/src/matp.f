      subroutine matp (a,n,nadim,k,b,nbdim,l,c,ncdim)
      implicit real*8 (a-h,o-z)
      dimension a(nadim,k),b(nbdim,l),c(ncdim,l)
      data zero/0.0/
      do 50 j=1,l
      do 50 i=1,n
      hq=zero
      do 49 ml6=1,k
 49   hq=hq+a(i,ml6)*b(ml6,j)
 50   c(i,j)=hq
      return
      end
      subroutine matpt(a,n,nadim,k,b,nbdim,l,c,ncdim)
      implicit real*8 (a-h,o-z)
      dimension a(nadim,k),b(nbdim,l),c(ncdim,l)
      data zero/0.0/
      do 50 i=1,k
      do 50 j=1,l
      hq=zero
      do 49 ml7=1,n
 49   hq=hq+a(ml7,i)*b(ml7,j)
 50   c(i,j)=hq
      return
      end
      DOUBLE PRECISION FUNCTION DOTV(X,N,Y)
      DOUBLE PRECISION X(N),Y(N)
C      DOTV WILL FORM THE DOT PRODUCT OF X AND Y IN DOUBLE PRECISON
      DOTV=0.
      DO 10 I=1,N
 10   DOTV=DOTV+X(I)*Y(I)
      RETURN
      END
      SUBROUTINE VEC(P,X,N,Q,Y)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION P(N),X(N),Y(N)
      DO 10 I=1,N
 10   Y(I)=P(I)+Q*X(I)
      RETURN
      END
      SUBROUTINE MOVE(X1,X2,M)
      DOUBLE PRECISION X1(M),X2(M)
      DO 100 I=1,M
 100  X2(I)=X1(I)
      RETURN
      END
      FUNCTION ZZD(N1,IA,IB,I1,I2,J1,J2,K1,K2,L1,L2,M1,M2,A,B)
C
C     THE TWO VECTORS TO BE MULTIPLIED ARE ELEMENTS OF THE MATRICES
C     A AND B. THE ELEMENTS TO BE USED MUST HAVE CONSTANT ADDRESS
C     DIFFERENCES (BUT MAY BE DIFFERENT BETWEEN A AND B). THE FIRST,
C     SECOND, AND LAST ELEMENTS OF A TO BE USED HAVE SUBSCRIPTS (I1,I2),
C     (J1,J2),(K1,K2). THE FIRST AND SECOND ELEMENTS OF B HAVE SUBSCRIPT
C     (L1,L2) AND (M1,M2).
C
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION A(IA,N1),B(IB,N1)
      SUM=0.
      IAD1=(I1-1)*IA+I2
      IAD2=(J1-1)*IA+J2
      IAD3=(K1-1)*IA+K2
      INC1=IAD2-IAD1
      INCT=IAD3-IAD1
      ICT=INCT/INC1+1
      IA1=J1-I1
      IA2=J2-I2
      IB1=M1-L1
      IB2=M2-L2
      DO 100 I=1,ICT
100   SUM=SUM+A(I1+(I-1)*IA1,I2+(I-1)*IA2)*B(L1+(I-1)*IB1,L2+(I-1)*IB2)
      ZZD=SUM
      RETURN
      END

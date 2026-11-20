c
c
c   imsl routine name   - ehobks
c 
c---------------------------------------------------------------------- 
c 
c   computer            - cdcft5/single 
c 
c   latest revision     - june 1, 1982
c 
c   purpose             - nucleus called only by imsl routine eigrs 
c 
c   precision/hardware  - single and double/h32 
c                       - single/h36,h48,h60
c 
c   reqd. imsl routines - none required 
c 
c   notation            - information on special notation and 
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c 
c   copyright           - 1982 by imsl, inc. all rights reserved. 
c 
c   warranty            - imsl warrants only that imsl testing has been 
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c 
c---------------------------------------------------------------------- 
c 
      subroutine ehobks (a,n,m1,m2,z,iz)
      implicit real*8(a-h,o-z)
c 
      dimension          a(31375),z(iz,iz) 
c                                  first executable statement 
      if (n .eq. 1) go to 30
      do 25 i=2,n 
         l = i-1
         ia = (i*l)/2 
         h = a(ia+i)
         if (h.eq.0.) go to 25
c                                  derives eigenvectors m1 to m2 of 
c                                  the original matrix from eigenvector 
c                                  m1 to m2 of the symmetric
c                                  tridiagonal matrix 
         do 20 j = m1,m2
            s = 0.0 
            do 10 k = 1,l 
               s = s+a(ia+k)*z(k,j) 
   10       continue
            s = s/h 
            do 15 k=1,l 
               z(k,j) = z(k,j)-s*a(ia+k)
   15       continue
   20    continue 
   25 continue
   30 return
      end 
c 
c 
c   imsl routine name   - ehouss
c 
c---------------------------------------------------------------------- 
c 
c   computer            - cdcft5/single 
c 
c   latest revision     - november 1, 1984
c 
c   purpose             - nucleus called only by imsl routine eigrs 
c 
c   precision/hardware  - single and double/h32 
c                       - single/h36,h48,h60
c 
c   reqd. imsl routines - none required 
c 
c   notation            - information on special notation and 
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c 
c   copyright           - 1982 by imsl, inc. all rights reserved. 
c 
c   warranty            - imsl warrants only that imsl testing has been 
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c 
c---------------------------------------------------------------------- 
c 
      subroutine ehouss (a,n,d,e,e2)
      implicit real*8(a-h,o-z)
c 
      dimension          a(31375),d(n),e(n),e2(n) 
      real*8             a,d,e,e2,zero,h,scale,f,g,hh 
      data               zero/0.0/
c                                  first executable statement 
      np1 = n+1 
      nn = (n*np1)/2-1
      nbeg = nn+1-n 
      do 70 ii = 1,n
         i = np1-ii 
         l = i-1
         h = zero 
         scale = zero 
         if (l .lt. 1) go to 10 
c                                  scale row (algol tol then not needed 
         nk = nn
         do 5 k = 1,l 
            scale = scale+abs(a(nk))
            nk = nk-1 
    5    continue 
         if (scale .ne. zero) go to 15
   10    e(i) = zero
         e2(i) = zero 
         go to 65 
   15    nk = nn
         do 20 k = 1,l
            a(nk) = a(nk)/scale 
            h = h+a(nk)*a(nk) 
            nk = nk-1 
   20    continue 
         e2(i) = scale*scale*h
         f = a(nn)
         g = -sign(sqrt(h),f) 
         e(i) = scale*g 
         h = h-f*g
         a(nn) = f-g
         if (l .eq. 1) go to 55 
         f = zero 
         jk1 = 1
         do 40 j = 1,l
            g = zero
            ik = nbeg+1 
            jk = jk1
c                                  form element of a*u
            do 25 k = 1,j 
               g = g+a(jk)*a(ik)
               jk = jk+1
               ik = ik+1
   25       continue
            jp1 = j+1 
            if (l .lt. jp1) go to 35
            jk = jk+j-1 
            do 30 k = jp1,l 
               g = g+a(jk)*a(ik)
               jk = jk+k
               ik = ik+1
   30       continue
c                                  form element of p
   35       e(j) = g/h
            f = f+e(j)*a(nbeg+j)
            jk1 = jk1+j 
   40    continue 
         hh = f/(h+h) 
c                                  form reduced a 
         jk = 1 
         do 50 j = 1,l
            f = a(nbeg+j) 
            g = e(j)-hh*f 
            e(j) = g
            do 45 k = 1,j 
               a(jk) = a(jk)-f*e(k)-g*a(nbeg+k) 
               jk = jk+1
   45       continue
   50    continue 
   55    do 60 k = 1,l
            a(nbeg+k) = scale*a(nbeg+k) 
   60    continue 
   65    d(i) = a(nbeg+i) 
         a(nbeg+i) = h*scale*scale
         nbeg = nbeg-i+1
         nn = nn-i
   70 continue
      return
      end 
c 
c 
c   imsl routine name   - eigrs 
c 
c---------------------------------------------------------------------- 
c 
c   computer            - cdcft5/single 
c 
c   latest revision     - june 1, 1980
c 
c   purpose             - eigenvalues and (optionally) eigenvectors of
c                           a real symmetric matrix 
c 
c   usage               - call eigrs (a,n,jobn,d,z,iz,wk,ier) 
c 
c   arguments    a      - input real symmetric matrix of order n, 
c                           whose eigenvalues and eigenvectors
c                           are to be computed. input a is
c                           destroyed if ijob is equal to 0 or 1. 
c                n      - input order of the matrix a.
c                jobn   - input option parameter.  if jobn.ge.10
c                         a is assumed to be in full storage mode 
c                         (in this case, a must be dimensioned exactly
c                         n by n in the calling program). 
c                         if jobn.lt.10 then a is assumed to be in
c                         symmetric storage mode.  define 
c                         ijob=mod(jobn,10).  then when 
c                           ijob = 0, compute eigenvalues only
c                           ijob = 1, compute eigenvalues and eigen-
c                             vectors.
c                           ijob = 2, compute eigenvalues, eigenvectors 
c                             and performance index.
c                           ijob = 3, compute performance index only. 
c                           if the performance index is computed, it is 
c                           returned in wk(1). the routines have
c                           performed (well, satisfactorily, poorly) if 
c                           wk(1) is (less than 1, between 1 and 100, 
c                           greater than 100).
c                d      - output vector of length n,
c                           containing the eigenvalues of a.
c                z      - output n by n matrix containing 
c                           the eigenvectors of a.
c                           the eigenvector in column j of z corres-
c                           ponds to the eigenvalue d(j). 
c                           if ijob = 0, z is not used. 
c                iz     - input row dimension of matrix z exactly as
c                           specified in the dimension statement in the 
c                           calling program.
c                wk     - work area, the length of wk depends 
c                           on the value of ijob, when
c                           ijob = 0, the length of wk is at least n. 
c                           ijob = 1, the length of wk is at least n. 
c                           ijob = 2, the length of wk is at least
c                             n(n+1)/2+n. 
c                           ijob = 3, the length of wk is at least 1. 
c                ier    - error parameter (output)
c                         terminal error
c                           ier = 128+j, indicates that eqrt2s failed 
c                             to converge on eigenvalue j. eigenvalues
c                             and eigenvectors 1,...,j-1 have been
c                             computed correctly, but the eigenvalues 
c                             are unordered. the performance index
c                             is set to 1000.0
c                         warning error (with fix)
c                           in the following, ijob = mod(jobn,10).
c                           ier = 66, indicates ijob is less than 0 or
c                             ijob is greater than 3. ijob set to 1.
c                           ier = 67, indicates ijob is not equal to
c                             zero, and iz is less than the order of
c                             matrix a. ijob is set to zero.
c 
c   precision/hardware  - single and double/h32 
c                       - single/h36,h48,h60
c 
c   reqd. imsl routines - ehobks,ehouss,eqrt2s
c 
c   notation            - information on special notation and 
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c 
c   copyright           - 1980 by imsl, inc. all rights reserved. 
c 
c   warranty            - imsl warrants only that imsl testing has been 
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c 
c---------------------------------------------------------------------- 
c 
      subroutine eigrs  (a,n,jobn,d,z,iz,wk,ier)
      implicit real*8(a-h,o-z)
c                                  specifications for arguments 
      integer            n,jobn,iz,ier
      real*8             a(31375),d(250),wk(500),z(iz,iz)
c                                  specifications for local variables 
      integer            ijob,ir,jr,ij,ji,np1 
      integer            jer,na,nd,iiz,ibeg,il,kk,lk,i,j,k,l
      real*8             anorm,asum,pi,sumz,sumr,an,s,ten,rdelp,zero, 
     1                   one,thous
      data               rdelp/.710543e-14/ 
      data               zero,one/0.0,1.0/,ten/10.0/,thous/1000.0/
c                                  initialize error parameters
c                                  first executable statement 
      ier = 0 
      jer = 0 
      if (jobn.lt.10) go to 15
c                                  convert to symmetric storage mode
      k = 1 
      ji = n-1
      ij = 1
      do 10 j=1,n 
         do 5 i=1,j 
            a(k) = a(ij)
            ij = ij+1 
            k = k+1 
    5    continue 
         ij = ij + ji 
         ji = ji - 1
   10 continue
   15 ijob = mod(jobn,10) 
      if (ijob.ge.0.and.ijob.le.3) go to 20 
c                                  warning error - ijob is not in the 
c                                    range
      ier = 66
      ijob = 1
      go to 25
   20 if (ijob.eq.0) go to 35 
   25 if (iz.ge.n) go to 30 
c                                  warning error - iz is less than n
c                                    eigenvectors can not be computed,
c                                    ijob set to zero 
      ier = 67
      ijob = 0
   30 if (ijob.eq.3) go to 75 
   35 na = (n*(n+1))/2
      if (ijob.ne.2) go to 45 
      do 40 i=1,na
         wk(i) = a(i) 
   40 continue
c                                  save input a if ijob = 2 
   45 nd = 1
      if (ijob.eq.2) nd = na+1
c                                  reduce a to symmetric tridiagonal
c                                    form 
      call ehouss (a,n,d,wk(nd),wk(nd)) 
      iiz = 1 
      if (ijob.eq.0) go to 60 
      iiz = iz
c                                  set z to the identity matrix 
      do 55 i=1,n 
         do 50 j=1,n
            z(i,j) = zero 
   50    continue 
         z(i,i) = one 
   55 continue
c                                  compute eigenvalues and eigenvectors 
   60 call eqrt2s (d,wk(nd),n,z,iiz,jer)
      if (ijob.eq.0) go to 9000 
      if (jer.gt.128) go to 65
c                                  back transform eigenvectors
      call ehobks (a,n,1,n,z,iz)
   65 if (ijob.le.1) go to 9000 
c                                  move input matrix back to a
      do 70 i=1,na
         a(i) = wk(i) 
   70 continue
      wk(1) = thous 
      if (jer.ne.0) go to 9000
c                                  compute 1 - norm of a
   75 anorm = zero
      ibeg = 1
      do 85 i=1,n 
         asum = zero
         il = ibeg
         kk = 1 
         do 80 l=1,n
            asum = asum+abs(a(il))
            if (l.ge.i) kk = l
            il = il+kk
   80    continue 
         if(anorm.lt.asum)anorm=asum
c        anorm = amax1(anorm,asum)
         ibeg = ibeg+i
   85 continue
      if (anorm.eq.zero) anorm = one
c                                  compute performance index
      pi = zero 
      do 100 i=1,n
         ibeg = 1 
         s = zero 
         sumz = zero
         do 95 l=1,n
            lk = ibeg 
            kk = 1
            sumz = sumz+abs(z(l,i)) 
            sumr = -d(i)*z(l,i) 
            do 90 k=1,n 
               sumr = sumr+a(lk)*z(k,i) 
               if (k.ge.l) kk = k 
               lk = lk+kk 
   90       continue
            s = s+abs(sumr) 
            ibeg = ibeg+l 
   95    continue 
         if (sumz.eq.zero) go to 100
         ssumz=s/sumz
         if(pi.lt.ssumz)pi=ssumz
c        pi = amax1(pi,s/sumz)
  100 continue
      an = n
      pi = pi/(anorm*ten*an*rdelp)
      wk(1) = pi
      if (jobn.lt.10) go to 9000
c                                  convert back to full storage mode
      np1 = n+1 
      ij = (n-1)*np1 + 2
      k = (n*(np1))/2 
      do 110 jr=1,n 
         j = np1-jr 
         do 105 ir=1,j
            ij = ij-1 
            a(ij) = a(k)
            k = k-1 
  105    continue 
         ij = ij-jr 
  110 continue
      ji = 0
      k = n-1 
      do 120 i=1,n
         ij = i-n 
         do 115 j=1,i 
            ij = ij+n 
            ji = ji+1 
            a(ij) = a(ji) 
  115    continue 
         ji = ji + k
         k = k-1
  120 continue
 9000 continue
      if (jer.eq.0) go to 9005
      ier = jer 
 9005 return
      end 
c
c
c   imsl routine name   - eqrt2s
c 
c---------------------------------------------------------------------- 
c 
c   computer            - cdcft5/single 
c 
c   latest revision     - november 1, 1984
c 
c   purpose             - eigenvalues and (optionally) eigenvectors of
c                           a symmetric tridiagonal matrix using the
c                           ql method.
c 
c   usage               - call eqrt2s (d,e,n,z,iz,ier)
c 
c   arguments    d      - on input, the vector d of length n contains 
c                           the diagonal elements of the symmetric
c                           tridiagonal matrix t. 
c                           on output, d contains the eigenvalues of
c                           t in ascending order. 
c                e      - on input, the vector e of length n contains 
c                           the sub-diagonal elements of t in position
c                           2,...,n. on output, e is destroyed. 
c                n      - order of tridiagonal matrix t.(input) 
c                z      - on input, z contains the identity matrix of 
c                           order n.
c                           on output, z contains the eigenvectors
c                           of t. the eigenvector in column j of z
c                           corresponds to the eigenvalue d(j). 
c                iz     - input row dimension of matrix z exactly as
c                           specified in the dimension statement in the 
c                           calling program. if iz is less than n, the
c                           eigenvectors are not computed. in this case 
c                           z is not used.
c                ier    - error parameter 
c                         terminal error
c                           ier = 128+j, indicates that eqrt2s failed 
c                             to converge on eigenvalue j. eigenvalues
c                             and eigenvectors 1,...,j-1 have been
c                             computed correctly, but the eigenvalues 
c                             are unordered.
c 
c   precision/hardware  - single and double/h32 
c                       - single/h36,h48,h60
c 
c   notation            - information on special notation and 
c                           conventions is available in the manual
c                           introduction or through imsl routine uhelp
c 
c   copyright           - 1978 by imsl, inc. all rights reserved. 
c 
c   warranty            - imsl warrants only that imsl testing has been 
c                           applied to this code. no other warranty,
c                           expressed or implied, is applicable.
c 
c---------------------------------------------------------------------- 
c 
      subroutine eqrt2s (d,e,n,z,iz,ier)
      implicit real*8(a-h,o-z)
c 
      dimension          d(250),e(250),z(iz,iz)
      data               rdelp/.710543e-14/ 
      data               zero,one/0.0,1.0/
c                                  move the last n-1 elements 
c                                  of e into the first n-1 locations
c                                  first executable statement 
      ier  = 0
      if (n .eq. 1) go to 9005
      do 5  i=2,n 
         e(i-1) = e(i)
    5 continue
      e(n) = zero 
      b = zero
      f = zero
      do  60  l=1,n 
         j = 0
         h = rdelp*(abs(d(l))+abs(e(l)))
         if (b.lt.h) b = h
c                                  look for small sub-diagonal element
         do 10  m=l,n 
            k=m 
            if (abs(e(k)) .le. b) go to 15
   10    continue 
   15    m = k
         if (m.eq.l) go to 55 
   20    if (j .eq. 30) go to 85
         j = j+1
         l1 = l+1 
         g = d(l) 
         p = (d(l1)-g)/(e(l)+e(l))
         r = abs(p) 
         if (rdelp*abs(p) .lt. 1.0) r = sqrt(p*p+one) 
         d(l) = e(l)/(p+sign(r,p))
         h = g-d(l) 
         do 25 i = l1,n 
            d(i) = d(i)-h 
   25    continue 
         f = f+h
c                                  ql transformation
         p = d(m) 
         c = one
         s = zero 
         mm1 = m-1
         mm1pl = mm1+l
         if (l.gt.mm1) go to 50 
         do 45 ii=l,mm1 
            i = mm1pl-ii
            g = c*e(i)
            h = c*p 
            if (abs(p).lt.abs(e(i))) go to 30 
            c = e(i)/p
            r = sqrt(c*c+one) 
            e(i+1) = s*p*r
            s = c/r 
            c = one/r 
            go to 35
   30       c = p/e(i)
            r = sqrt(c*c+one) 
            e(i+1) = s*e(i)*r 
            s = one/r 
            c = c*s 
   35       p = c*d(i)-s*g
            d(i+1) = h+s*(c*g+s*d(i)) 
            if (iz .lt. n) go to 45 
c                                  form vector
            do 40 k=1,n 
               h = z(k,i+1) 
               z(k,i+1) = s*z(k,i)+c*h
               z(k,i) = c*z(k,i)-s*h
   40       continue
   45    continue 
   50    e(l) = s*p 
         d(l) = c*p 
         if ( abs(e(l)) .gt.b) go to 20 
   55    d(l) = d(l) + f
   60 continue
c                                  order eigenvalues and eigenvectors 
      do  80  i=1,n 
         k = i
         p = d(i) 
         ip1 = i+1
         if (ip1.gt.n) go to 70 
         do 65  j=ip1,n 
            if (d(j) .ge. p) go to 65 
            k = j 
            p = d(j)
   65    continue 
   70    if (k.eq.i) go to 80 
         d(k) = d(i)
         d(i) = p 
         if (iz .lt. n) go to 80
         do 75 j = 1,n
            p = z(j,i)
            z(j,i) = z(j,k) 
            z(j,k) = p
   75    continue 
   80 continue
      go to 9005
   85 ier = 128+l 
 9000 continue
 9005 return
      end 

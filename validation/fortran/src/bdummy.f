
cvms  *************************
c     to bring the data block information in file bdummy into the gqopt
c     linking, the data block in file bdummy is given the name bdummy
c     block data
      block data bdummy
cvms  *************************
      implicit real*8(a-h,o-z)
      character*1 inlist,newpag,nlist,clist
      common/bdum/i
      common/binput/inflg
      common/bspd/ispd
      common/bgrdx1/ntrm,nbhhh,lbhhh,kbhhh
      common/bnobs/nobss
      common/bopt/ix,lt,ifp,isp,nloop,ist,iloop
c      common/bstak/nq,ntop
      common/bstak/nq
      common/bfidif/fdfrac,fdmin
      common/bgrdx2/rtm,eig,rtmult
      common/bgrdx3/acti,actw,zful,betm
      common/bdfp/stpmin,fmin
      common/bpat/stpsiz,shrink
      common/bprint/ipt,nfile,ndig,npunch,jpt,mfile
      common/bread/nread
      common/bread1/ifile
      common/blnsr/step1,stpacc,nlnsr
      common/bnel/stp1,var,konvge,nrst
      common/bnlys2/inlist,newpag,nlist,clist
      common/bstop/nvar1,istop(3)
      common/bprec/rsmall,absmal
      common/btrat/itrflg
      common/bopt5/pm2,jthet,jphi,jran,jsave
      common/bopt6/deltaz,sdelta,rho
      common/bopt7/alpha,beta,omega,epstil,rh ,deltil,sigma,bigl,akappa
      common/bopt8/iv8,iv3
      common/bnwrp1/ivflim,invlim
      common/bnwrp2/epcor,epsvrt
      common/bssvm1/thetaa,phii,shrip1,shrip2,epss
      common/bpnlt1/shrk,xl1,xl2,ivers1
      common/bprob/iprob,jprob
      common/bprob2/kprob
      common/bstop1/jstop(2)
      common/btob1/itob
cvms  *************************
c     these are the initial declarations used for the common
c     block used by the random number generator ggnml (from
c     the imsl package) found in subroutine rxx.
      common /vmsran/ vmseed,vmsnr,vmsr
      real*8  vmseed
      integer vmsnr
      real    vmsr(1)
      data vmsnr /0/
cvms  *************************
      data ivflim,invlim/7,7/
      data epcor,epsvrt/0.1d0,1.d-15/
      data rsmall,absmal/1.d-14,1.d-35/
      data inlist,newpag,nlist,clist/'t','t','t','t'/
      data inflg/0/,nbhhh/0/
      data ispd /1/
      data nobss/0/
      data stp1,konvge,nrst/5.d-2,-1,0/
      data nvar1,istop/-1,0,0,0/
      data ix,lt,ifp,isp,nloop,ist/1,0,3,2,7,1/
      data ipt,nfile,ndig,npunch,nread/2,6,8,9,9/
      data jpt,mfile/0,6/
      data ifile/10/
      data stpmin,fmin/1.d-10,0.d0/
c      data nq,ntop/1,0/
      data nq/1/
      data rtm,eig,rtmult,fdfrac,fdmin/.25d0,.1d-07,.115d01,.1d-3,.1d-5/
      data thetaa,phii,shrip1,shrip2,epss/.25d0,.5d0,.2d0,2.d0,.5d-1/
      data acti,actw,zful,betm/.2d0,.2d0,.2d01,.5d0/
      data stpsiz,shrink/1.d0,.5d0/
      data step1,stpacc,nlnsr/1.d0,-1.d0,20/
      data iprob/0/,jprob/1/,kprob/0/
      data itob/0/
      data jstop/0,0/
      data deltaz,sdelta,rho/4.d0,1.d-4,5.d-1/
      data shrk/0.65d0/,xl1/0.1d0/,xl2/0.1d0/
      data ivers1/0/
      data iv8/2/,iv3/0/
      data alpha,beta,omega,epstil,rh,deltil,sigma,bigl,akappa/0.1,0.6,
     1  0.3,1.d-3,0.3,0.1,0.5,1000.0,30.0/
c   stpacc will be reset by opt for default value of acc
      data itrflg/0/
      end

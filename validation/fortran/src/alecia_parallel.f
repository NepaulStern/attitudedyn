c
c       this program estimates the alecia model.
c
c                       main ->      covar,estim1,estim2,initia,
c                                    likely11,likely12,lmtests,
c                                    olsresids,parset0,vardecomp,
c                                    vartest1,vartest2,waldtests
c normadj1,normadj2, -> anormc
c probmake1,probmake2
c normadj2,probmake2 -> anormd
c estim1,likely11, ->   bestupd 
c likely12,vartest1
c main,estim1,estim2 -> covar ->     eigrs,likobs11,likobs12,linv3f,
c                                    obsmake1,obsmake2,parset0,parset1,
c                                    parset211,parset212
c estim2 ->             cpu_time
c datget1,lmtests ->    datextraget
c initia ->             datget1 ->   datextraget
c initia ->             datget2
c initia ->             datget3
c covar,lmcovar, ->     eigrs
c vartest1
c likobs11 ->           errget1
c likobs12, ->          errget2
c likobs12t,likobs2
c initia ->             errmake ->   mdnris,randfn,utrunc
c main ->               estim1 ->    bestupd,covar,likely11,likely12,
c                                    opt,parset0,vardecomp
c main ->               estim2 ->    covar,cpu_time,likely2,linv3f,
c                                    parset0
c likobs12,likobs12t -> ghk1 ->      linv3f
c likobs2 ->            ghk2 ->      linv3f
c main ->               initia ->    datget1,datget2,datget3,errmake,
c                                    parget
c main,estim1, ->       likely11 ->  bestupd,likobs11,obsmake1,parset1,
c nmsimp,vartest1                    parset211
c main,estim1, ->       likely12 ->  bestupd,likobs12,obsmake2,parset1,
c nmsimp,vartest1                    parset212,probmake1
c estim2 ->             likely2 ->   likobs2,obsmake2,parset1,parset22,
c                                    probmake2
c covar,likely11, ->    likobs11 ->  errget1,normadj1,preplik11,
c vartest1                           probmake1
c covar,likely12, ->    likobs12 ->  errget2,ghk1,normadj1,preplik12,
c vartest1,waldtests                 probmake1
c lmcovar ->            likobs12t -> errget2,ghk1,normadj1,preplik12t,
c                                    probmake1
c likely2,vartest2 ->   likobs2 ->   errget2,ghk2,normadj2,preplik2
c covar,estim2,ghk1, -> linv3f ->    ludatn,luelmn
c ghk2,olsresids,
c waldtests
c lmtests ->            lmcovar ->   likobs12t,lmcovar,obsmake2,
c                                    parset0t,parset1t,parset212t
c main ->               lmtests ->   datextraget,lmcovar
c linv3f ->             ludatn
c linv3f ->             luelmn
c errmake,normadj1, ->  mdnris ->    merfi
c normadj2
c mdnris ->             merfi
c opt ->                nmsimp ->    likely11,likely12
c likobs11,likobs12, -> normadj1 ->  anormc,mdnris,utrunc
c likobs12t
c likobs2               normadj2 ->  anormc,anormd,mdnris,utrunc
c covar,likely11, ->    obsmake1
c vartest1
c covar,likely12, ->    obsmake2
c likely2,lmcovar,
c vardecomp,vartest2,
c waldtests
c main ->               olsresids -> linv3f
c parset212, ->         omegamake1
c parset212t
c parset22 ->           omegamake2
c estim1 ->             opt ->       nmsimp
c initia ->             parget
c main,covar,estim1, -> parset0
c estim2,vartest1,
c vartest2,waldtests
c lmcovar ->            parset0t
c covar,likely11, ->    parset1
c likely12,likely2,
c vartest1,vartest2,
c waldtests
c lmcovar ->            parset1t
c covar,likely11, ->    parset211
c vartest1
c covar,likely12, ->    parset212 ->  omegamake1
c vardecomp,vartest1,
c waldtests
c lmcovar ->            parset212t -> omegamake1
c likely2,vartest2 ->   parset22 ->   omegamake2
c likobs11 ->           preplik11
c likobs12,vardecomp -> preplik12
c likobs12t ->          preplik12t
c likobs2 ->            preplik2
c likobs11, ->          probmake1 ->  anormc
c likobs12t,likobs12
c likobs2 ->            probmake2 ->  anormc,anormd
c errmake ->            randfn
c errmake,normadj1, ->  utrunc
c normadj2
c main,estim1 ->        vardecomp ->  obsmake2,parset212,preplik12
c main ->               vartest1 ->   bestupd,eigrs,likely11,likely12,
c                                     likobs11,likobs12,obsmake1,
c                                     parset0,parset1,parset211,
c                                     parset212
c main ->               vartest2 ->   likobs2,obsmake2,parset0,parset1,
c                                     parset22
c main ->               waldtests ->  likobs12,linv3f,obsmake2,parset0,
c                                     parset1,parset212
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+(nys*nsegm2s)+
     & (nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension parm(ncoefs)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
c
c       initialize.
c
      call initia
c
      if(method.eq.0)then
         call parset0(parm,nparm)
         if(model.eq.1)call likely11(parm,nparm,obj,*3)
         if(model.eq.2)call likely12(parm,nparm,obj,*3)
   3     continue
         stop
         endif
      if(method.eq.1)call estim1
      if(method.eq.2)call estim2
      if(method.eq.3)call covar
      if(method.eq.4)call vardecomp
      if(method.eq.5)call lmtests
      if(method.eq.6)call waldtests
      if(method.eq.7)call olsresids
      if(method.eq.8)call vartest1
      if(method.eq.9)call vartest2
      stop
      end
c
c
      subroutine anormc(arg,val)
c 
c       this subroutine evaluates the standard normal distribution
c       function.
c 
      real xt,valt
      real*8 arg,val
c 
c       standardize the argument of the distribution function 
c 
      xt=sngl(arg)*.707107
c 
c       evaluate using the error function 
c 
      if(xt.lt.0.)valt=1.-(.5*erfc(xt)) 
      if(xt.ge.0.)valt=.5*erfc(-xt) 
      val=dble(valt)
      return
      end 
c
c
      subroutine anormd(arg,farg)
      implicit real*8(a-h,o-z)
      data piinv/0.398942d0/
      farg=piinv*exp(-0.5d0*arg*arg)
      return
      end
c
c
      subroutine bestupd(icase,f)
c
c       this subroutine keeps track of and outputs best f values.
c
      implicit real*8(a-h,o-z)
      common/best/fbest(3)
      if(icase.eq.0)then
         do 3 i=1,3
   3        fbest(i)=-1.d20
         return
         endif
      if(icase.eq.1)then
         if(f.gt.fbest(3))then
            if(f.gt.fbest(2))then
               fbest(3)=fbest(2)
               if(f.gt.fbest(1))then
                  fbest(2)=fbest(1)
                  fbest(1)=f
                  goto 2
                  endif
               fbest(2)=f
               goto 2
               endif
            fbest(3)=f
            endif
c   2     write(6,100)f,(fbest(i),i=1,3)
c 100     format(1x,'f = ',g15.8,2x,'fbest =',3(2x,g15.8))
         endif
   2   return
      end
c
c
      subroutine covar
c
c       this subroutine computes the asymptotic covariance matrix and 
c       reports the statistical significance of the parameter 
c       estimates.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2),ncoefss=
     & ncoefs*(ncoefs+1)/2,ncoef2s=ncoefs*2)
      implicit real*8(a-h,o-z)
      character*2 start,star(3)
      character*8 alabel
      real xextra,xv
      dimension amean(ncoefs),cov(ncoefs,ncoefs),covs(ncoefss),
     & dft(ncoefs),eigval(ncoefs),eigvec(ncoefs,ncoefs),parm(ncoefs),
     & work1(ncoefs),work2(ncoef2s)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/labels/alabel(ncoefs)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      data eigtol/0.001d0/
      data star/'  ','* ','**'/
c
c       initialize.
c
      do 13 i=1,ncoef
  13     amean(i)=0.d0
      do 6 i=1,ncoef
         do 7 j=1,i
   7        cov(i,j)=0.d0
   6     continue
      call parset0(parm,nparm)
c
c       loop over observations.
c
      anobs=0.d0
      do 2 i=1,nobs
         if(((i/1000)*1000).eq.i)write(6,111)i
 111     format(1x,'i = ',i5)
         if(incluv(i).eq.1)then
            anobs=anobs+1.d0
c
c       evaluate at the parameter estimate.
c
            call parset1(parm)
            if(model.eq.1)then
               call parset211
               call obsmake1(i)
               call likobs11(f0)
               endif
            if(model.eq.2)then
               call parset212(penalty)
               call obsmake2(i)
               call likobs12(f0)
               f0=f0+penalty
               endif
c
c       loop over parameters and compute numerical derivatives.
c
            iparm=0
            do 3 j=1,ncoef
               if(icoef(j).eq.0)then
                  iparm=iparm+1
                  delta=coef(j)/100.d0
                  if(delta.gt.1.d0)delta=1.d0
                  if(delta.lt.0.01d0)delta=0.01d0
                  coef(j)=coef(j)+delta
                  if(model.eq.1)then
                     call parset211
                     call likobs11(f1)
                     endif
                  if(model.eq.2)then
                     call parset212(penalty)
                     call likobs12(f1)
                     f1=f1+penalty
                     endif
                  dft(iparm)=(f1-f0)/delta
                  coef(j)=coef(j)-delta              
                  endif
   3           continue
c
c       update mean and covariance matrix.
c
            do 4 j=1,iparm
               amean(j)=amean(j)+dft(j)
               do 5 k=1,j
   5              cov(j,k)=cov(j,k)+(dft(j)*dft(k))
   4           continue
            endif
   2     continue
c
c       update and output mean derivatives.
c
      write(6,100)
 100  format(1x,50('='))
      do 8 i=1,iparm
   8     amean(i)=amean(i)/anobs
      ij=0
      do 9 i=1,iparm
         do 10 j=1,i
            cov(i,j)=(cov(i,j)/anobs)-(amean(i)*amean(j))
            cov(j,i)=cov(i,j)
            ij=ij+1
  10        covs(ij)=cov(i,j)
   9     continue
      write(6,109)
 109  format(1x,'mean derivatives')
      write(6,101)
 101  format(6x,'variable',2x,'mean',8x,'std err',3x,'t-stat')
      do 11 i=1,iparm
         stderr=(cov(i,i)/anobs)**.5d0
         tstat=amean(i)/stderr
         atstat=abs(tstat)
         if(atstat.lt.1.53d0)start=star(1)
         if((atstat.ge.1.53d0).and.(atstat.le.1.96d0))start=star(2)
         if(atstat.gt.1.96d0)start=star(3)
  11     write(6,102)i,alabel(i),amean(i),start,stderr,tstat
 102  format(1x,i3,2x,a8,2x,f8.3,a2,2(2x,f8.3))
c
c       compute eigenvalues.
c
      iflag=0
      call eigrs(covs,iparm,1,eigval,eigvec,ncoefs,work1,ier)
      do 12 i=1,iparm
         if(eigval(i).lt.eigtol)then
            if(iflag.eq.0)then
               iflag=1
               write(6,103)
 103           format(1x,'problem eigenvalues and eigenvectors',/,1x,
     &          6x,'eigenvalue',7x,'eigenvector')
               endif
            write(6,104)i,eigval(i),(eigvec(j,i),j=1,iparm)
 104        format(1x,i3,2x,g15.8,5(2x,f7.4),50(/,21x,5(2x,f7.4)))
            endif
  12     continue
c
c       invert.
c
      d1=-1.d0
      call linv3f(cov,work1,1,iparm,ncoefs,d1,d2,work2,ier)
      if(ier.ne.0)write(6,106)ier
 106  format(1x,'ier = ',i3)
      write(6,100)
      write(6,112)
 112  format(1x,'asymptotic covariance matrix')
      write(6,113)
 113  format(5x,'covariance')
      do 14 i=1,iparm
  14     write(6,114)i,(cov(i,j),j=1,i)
 114  format(1x,i3,5(2x,f8.3),50(/,4x,5(2x,f8.3)))
      write(6,100)
      write(6,110)
 110  format(1x,'estimates and significance')
      write(6,105)
 105  format(6x,'variable',4x,'estimate',2x,'std err',3x,'t-stat')
      do 15 i=1,iparm
         stderr=(cov(i,i)/anobs)**.5d0
         tstat=parm(i)/stderr
         atstat=abs(tstat)
         if(atstat.lt.1.53d0)start=star(1)
         if((atstat.ge.1.53d0).and.(atstat.le.1.96d0))start=star(2)
         if(atstat.gt.1.96d0)start=star(3)
  15     write(6,108)i,alabel(i),parm(i),start,stderr,tstat
 108  format(1x,i3,2x,a8,2x,f8.3,a2,2(2x,f8.3))
      return
      end
c
c
      subroutine datextraget(icase)
c
c       this subroutine gets extra data for a better age variable and
c       one of the lagrange multiplier tests.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nrs=10,nsegs=8,nsegm2s=
     & nsegs-2,ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*
     & nbetas)+(nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      real xextra,xextratt,xv
      dimension iagex(nobss,ntimes),xexmom(5),xmom(5)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      if(icase.eq.1)then
         open(unit=8,file='data/besextradata.prn')
         nextra=57000
         do 10 i=1,nobs
            do 11 j=1,ntime
  11           xextra(i,j)=-99.
  10        continue
         do 12 i=1,3
  12        xexmom(i)=0.d0
         xexmom(4)=1.d20
         xexmom(5)=-1.d20
         do 13 i=1,nextra
            read(8,*)idextra,itimeextra,iageextra
            if(iageextra.lt.0)goto 13   
            do 14 j=1,nobs
               if(idextra.eq.idv(j))then
                  iagex(j,itimeextra)=iageextra
                  if(incluv(i).eq.1)then
                     if(model.eq.1)then
                        if(iageextra.lt.0)iaget=8
                        if(iageextra.gt.68)iaget=4
                        if((iageextra.le.68).and.(iageextra.ge.50))
     &                   iaget=5
                        if((iageextra.le.51).and.(iageextra.ge.34))
     &                   iaget=6
                        if((iageextra.le.33).and.(iageextra.ge.18))
     &                   iaget=7
                        endif
                     if(model.eq.2)then
                        iaget=0
                        if(iageextra.lt.0)iaget=9
                        if(iageextra.gt.68)iaget=6
                        if((iageextra.le.68).and.(iageextra.ge.50))
     &                   iaget=7
                        if((iageextra.lt.50).and.(iageextra.ge.34))
     &                   iaget=8
                        do 15 k=4,8
  15                       xv(j,itimeextra,k)=0.
                        if(iaget.gt.0)xv(j,itimeextra,iaget)=1.
                        endif
                     endif
                  endif
  14           continue
  13        continue
         close(8)
         endif
      if(icase.eq.2)then
         open(unit=8,file='data/besextradata.prn')
         nextra=57000
         do 4 i=1,nobs
            do 5 j=1,ntime
   5           xextra(i,j)=-99.
   4        continue
         do 6 i=1,3
   6        xmom(i)=0.d0
         xmom(4)=1.d20
         xmom(5)=-1.d20
         do 2 i=1,nextra
            read(8,*)idextra,itimeextra,iageextra
            if(iageextra.lt.0)goto 2
            do 3 j=1,nobs
               if(idextra.eq.idv(j))then
                  iagex(j,itimeextra)=iageextra
                  xextra(j,itimeextra)=0.
                  if((iageextra.ge.16).and.(iageextra.le.25))
     &             xextra(j,itimeextra)=1.
                  xextrat=dble(xextra(j,itimeextra))
                  xmom(1)=xmom(1)+1.d0
                  xmom(2)=xmom(2)+xextrat
                  xmom(3)=xmom(3)+(xextrat*xextrat)
                  if(xmom(4).gt.xextrat)xmom(4)=xextrat
                  if(xmom(5).lt.xextrat)xmom(5)=xextrat
                  goto 2
                  endif
   3           continue
   2        continue
         close(8)
         do 7 i=1,nobs
            do 8 j=1,ntime
               if(xextra(i,j).lt.0.)then
                  do 9 k=1,ntime
                     if(xextra(i,k).gt.0.)then
                        adjust=(time(j)-time(k))/365.d0
                        xextratt=float(iagex(i,k))+sngl(adjust)
                        xextra(i,j)=0.
                        if((xextratt.ge.16.).and.(xextratt.le.25.))
     &                   xextra(i,j)=1.d0
                        xextrat=dble(xextra(i,j))
                        xmom(1)=xmom(1)+1.d0
                        xmom(2)=xmom(2)+xextrat
                        xmom(3)=xmom(3)+(xextrat*xextrat)
                        if(xmom(4).gt.xextrat)xmom(4)=xextrat
                        if(xmom(5).lt.xextrat)xmom(5)=xextrat
                        goto 8
                        endif
   9                 continue
                  endif
   8           continue
   7        continue
         xmom(2)=xmom(2)/xmom(1)
         xmom(3)=((xmom(3)/xmom(1))-(xmom(2)*xmom(2)))**.5d0
         write(6,102)(xmom(i),i=1,5)
 102     format(1x,'# extra observations = ',f8.0,/,1x,'mean = ',f6.3,
     &    /,1x,'std dev = ',f6.3,/,1x,'minimum = ',f6.3,/,1x,
     &    'maximum = ',f6.3)
         endif
      return
      end
c
c
      subroutine datget1
c
c       this subroutine gets data from the british election panel 
c       study.
c
      parameter(nbetas=20,nbounds=200,nfacs=2,nobss=96529,nsegs=8,
     & nsegm2s=nsegs-2,nsegp1s=nsegs+1,ntimes=34,ntimep2s=ntimes+2,
     & nxs=nbetas,ncauses=nxs+1,nx2s=nxs*2,nys=6,nyp1s=nys+1,
     & ny3s=nys*3,ncoefs=(nys*ntimes)+(nys*nbetas)+(nys*nsegm2s)+
     & (nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      character*2 yesno(2)
      character*8 aunitan(4),ax1(nx2s),ax2(nx2s),ay2(ny3s)
      real xextra,xinc(15),xv
      dimension ibound(nbounds,3),icause(ntimes,ncauses),
     & icodens(500,4),icodent(3),icountry(ntimes),id(nobss),
     & ideppat(200,nyp1s),iduvot(ntimes),ieduc(ntimes),
     & ieducbdens(100,3),ieducdens(30000,7),ieducdent(7),
     & ieducr(ntimes),ieffic(ntimes,2),income(ntimes),
     & interden(2,3,ntimep2s),inumdens(ntimes),ipartidst(ntimes),
     & iparty(ntimes),ipartysq(ntimes),irelig(ntimes),isateu(ntimes),
     & isatuk(ntimes),itimet(ntimes),itrustmp(ntimes),iwavet(ntimes),
     & iwindex(39),iydens(nys,nsegp1s),iyvb(nys),timet(ntimes,3),
     & xmom(nxs,6),xmt(nxs,3),ymom(nys,6),ymt(nys,3)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      data aunitan/'observat','ion','period',' '/
      data ax1/'constant',' ','female',' ','educatio','n','coh: < 1',
     & '946','coh: 194','6-1964','coh: 196','5-1980','coh: 198',
     & '1-1996','age miss','ing','income/1','00k','income m','issing',
     & 'religiou','s','conserva','tive','labor',' ','other pa','rty',
     & 'party st','rength','party va','rs missg','satisfie','d uk dem',
     & 'satisfie','d eu dem','mp trust',' ',' ',' '/
      data ax2/'constant',' ','female',' ','hs diplo','ma',
     & 'college ','degree','post-gra','d degree','coh: < 1','946',
     & 'coh: 194','6-1964','coh: 196','4-1980','age miss','ing',
     & 'income <',' 25K','25K< inc','ome <50K','income m','issing',
     & 'religiou','s','conserva','tive','labor',' ','other pa','rty',
     & 'prty str','ngth = 1','prty str','ngth = 2','scotland',' ',
     & 'wales',' '/
      data ay2/'duty to ','vote',' ','satisfie','d w/ uk ','democracy',
     & 'trust mp','s',' ','efficacy',': unders','tand','efficacy',
     & ': pol ca','re',' ',' ',' '/
      data iwindex/5,9,10,11,13,16,19,24,28,29,30,31,33,37,43,46,47,
     & 48,49,51,52,62,66,67,68,69,70,71,75,76,81,85,86,87,88,89,90,94,
     & 95/
      data xinc/0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,
     & 0.475,0.55,0.65,0.85,1.25,2.0/
      data yesno/'  ','no'/
c
c       initialize.
c
      do 5 i=1,ntime
         timet(i,1)=0.d0
         timet(i,2)=0.d0
   5     timet(i,3)=0.d0
c
c       variables if model = 1:
c       dependent variable (iyv): duty to vote (ordered discrete) 
c       explanatory variables:
c       1) constant
c       2) female
c       3) education dummy (at least college)
c       4) cohort: < born 1946
c       5) cohort: born 1946-1964
c       6) cohort: born 1965-1980
c       7) cohort: born 1981-1996
c       8) age missing
c       9) income/100k
c      10) income missing
c      11) religious
c      12) party: conservative
c      13) party: labor
c      14) other party
c      15) party strength > 0
c      16) party strength missing
c      17) satisfied w/ uk democracy
c      18) satisfied w/ eu democracy
c      19) trust mps
c
c       variables if model = 2:
c       dependent variables (iyv):
c       1) duty to vote (ordered discrete w/ 5 levels)
c       2) satisfaction with uk democracy (ordered discrete w/ 4 
c          levels)
c       3) trust mps (ordered discrete w/ 7 levels)
c       4) efficacy: understand (ordered discrete w/ 5 levels)
c       5) efficacy: pol care (ordered discrete w/ 5 levels)
c
c       explanatory variables:
c       1) constant
c       2) female
c       3) hs diploma
c       4) college degree
c       5) post-grad degree
c       6) cohort: born < 1946
c       7) cohort: born 1946-1964
c       8) cohort: born 1965-1980
c       9) age missing
c      10) income < 25K
c      11) 25K < income < 50K
c      12) income missing
c      13) religious
c      14) party: conservative
c      15) party: labor
c      16) other party
c      17) party strength = 1
c      18) party strength = 2
c      19) scotland
c      20) wales
c
c       read data from aleciadata1.
c
      write(6,125)
 125  format(1x,'read aleciadata1c.prn')
      open(unit=8,file='data/aleciadata1c.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                     4    8    3x,i5 
c         2) iwavet(1)       dummy for wave 1        16   16    7x,i1
c         3) iwavet(2)       dummy for wave 2        24   24    7x,i1
c         4) iwavet(3)       dummy for wave 3        32   32    7x,i1
c         5) iwavet(4)       dummy for wave 4        40   40    7x,i1
c         6) iwavet(5)       dummy for wave 5        48   48    7x,i1
c         7) iwavet(6)       dummy for wave 6        56   56    7x,i1
c         8) iwavet(7)       dummy for wave 7        64   64    7x,i1
c         9) iwavet(8)       dummy for wave 8        72   72    7x,i1
c        10) iwavet(9)       dummy for wave 9        80   80    7x,i1
c        11) iwavet(10)      dummy for wave 10       88   88    7x,i1
c        12) iwavet(11)      dummy for wave 11       96   96    7x,i1
c        13) iwavet(12)      dummy for wave 12      104  104    7x,i1
c        14) iwavet(13)      dummy for wave 13      112  112    7x,i1
c        15) iwavet(14)      dummy for wave 14      120  120    7x,i1
c        16) iwavet(15)      dummy for wave 15      128  128    7x,i1
c        17) iwavet(16)      dummy for wave 16      136  136    7x,i1
c        18) iwavet(17)      dummy for wave 17      144  144    7x,i1
c        19) iwavet(18)      dummy for wave 18      152  152    7x,i1
c        20) iwavet(19)      dummy for wave 19      160  160    7x,i1
c        21) itimet(1)       date of wave 1         164  168    3x,i5
c        22) itimet(2)       date of wave 2         172  176    3x,i5
c        23) itimet(3)       date of wave 3         180  184    3x,i5
c        24) itimet(4)       date of wave 4         188  192    3x,i5
c        25) itimet(5)       date of wave 5         196  200    3x,i5
c        26) itimet(6)       date of wave 6         204  208    3x,i5
c        27) itimet(7)       date of wave 7         212  216    3x,i5
c        28) itimet(8)       date of wave 8         220  224    3x,i5
c        29) itimet(9)       date of wave 9         228  232    3x,i5
c        30) itimet(10)      date of wave 10        236  240    3x,i5
c        31) itimet(11)      date of wave 11        244  248    3x,i5
c        32) itimet(12)      date of wave 12        252  256    3x,i5
c        33) itimet(13)      date of wave 13        260  264    3x,i5
c        34) itimet(14)      date of wave 14        268  272    3x,i5
c        35) itimet(15)      date of wave 15        276  280    3x,i5
c        36) itimet(16)      date of wave 16        284  288    3x,i5
c        37) itimet(17)      date of wave 17        292  296    3x,i5
c        38) itimet(18)      date of wave 18        300  304    3x,i5
c        39) itimet(19)      date of wave 19        308  312    3x,i5
c                                                             153x
      ico1=0
      ico2=0
      do 2 i=1,nobs
         read(8,100)idt,(iwavet(j),j=1,ntime),(itimet(j),j=1,19)
 100     format(3x,i5,19(7x,i1),19(3x,i5),153x,11(3x,i5))
c
c       construct id, iobsrvv, istartv, and iendv, and update timet.
c
         iflag=0
         do 3 j=1,ntime
            iobsrvv(i,j)=iwavet(j)
            if((j.eq.5).or.((j.ge.9).and.(j.le.11)))then
               iobsrvv(i,j)=0
               iwavet(j)=0
               endif
            if((iwavet(j).ne.1).and.(j.ne.5).and.((j.lt.9).or.
     &       (j.gt.11)))then
               if(itimet(j).ne.-999)then
                  ico1=ico1+1
c                 write(6,101)i,j,iwavet(j),itimet(j),ico1
                  endif
               endif
            if(iwavet(j).eq.1)then
               if(itimet(j).lt.0)then
                  write(6,101)i,j,iwavet(j),itimet(j)
 101              format(1x,'mismatch for iwavet and itimet: i,j,',
     &             'iwavet,itimet: ',i5,2(2x,i2),2(2x,i5))
                  stop
                  endif
               if(iflag.eq.0)then
                  iflag=1
                  id(i)=idt
                  istartv(i)=j
                  endif
               xv(i,j,1)=1.
               timet(j,1)=timet(j,1)+1.d0
               timet(j,2)=timet(j,2)+dble(float(itimet(j)))
               timet(j,3)=timet(j,3)+dble(float(itimet(j)*itimet(j)))
               endif
   3        continue
         if(iflag.eq.0)then
            id(i)=idt
            ico2=ico2+1
            incluv(i)=0
            goto 2
            endif
         incluv(i)=1
         istart=istartv(i)
         jt=ntime+1
         do 4 j=istart,ntime
            jt=jt-1
            if(iobsrvv(i,jt).eq.1)then
               iendv(i)=jt
               goto 2
               endif
   4        continue
   2     continue
c
c       adjust and output timet.
c
      write(6,103)
 103  format(1x,'average times')
      write(6,104)
 104  format(1x,'wave',2x,'# obs',5x,'mean',6x,'std dev')
      do 7 i=1,ntime
         if(timet(i,1).le.0.d0)then
            if(i.eq.5)timet(i,2)=20191.31d0
            if(i.eq.9)timet(i,2)=20693.23d0
            if(i.eq.10)timet(i,2)=20782.39d0
            if(i.eq.11)timet(i,2)=20871.55d0
            timet(i,3)=0.d0
            write(6,156)i,(timet(i,j),j=1,3)
 156        format(1x,i2,4x,f8.0,2x,f8.1,2x,f8.3,5x,'interpolated')
            if(icoef(i).eq.0)then
               write(6,157)i
 157           format(1x,'icoef(',i2,') = 0')
               stop
               endif
            endif
         if(timet(i,1).gt.0.d0)then
            timet(i,2)=timet(i,2)/timet(i,1)
            timet(i,3)=((timet(i,3)/timet(i,1))-(timet(i,2)*
     &       timet(i,2)))**.5d0
            write(6,105)i,(timet(i,j),j=1,3)
 105        format(1x,i2,4x,f8.0,2x,f8.1,2x,f8.3)
            endif
   7     time(i)=timet(i,2)
      if((ico1.gt.0).or.(ico2.gt.0))write(6,124)ico1,ico2
 124  format(1x,'# cases w/ iwave eq 0 and itime ne 0: ',i5,/,1x,
     & '# cases w/ no useful waves: ',i5)
      close(8)
c
c       read data from aleciadata2.
c
      write(6,126)
 126  format(1x,'read aleciadata2c.prn')
      open(unit=8,file='data/aleciadata2c.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                     4    8    3x,i5 
c         2) ipartidst(1)    party id strength       12   16    3x,i5
c                            in wave 1
c         3) iparty(1)       party id in wave 1      20   24    3x,i5
c         4) ipartysq(1)     party id squeeze in     28   32    3x,i5
c                            wave 1
c         5) ipartidst(2)    party id strength       36   40    3x,i5
c                            in wave 2
c         6) iparty(2)       party id in wave 2      44   48    3x,i5
c         7) ipartysq(2)     party id squeeze in     52   56    3x,i5
c                            wave 2
c         8) ipartidst(3)    party id strength       60   64    3x,i5       
c                            in wave 3
c         9) iparty(3)       party id in wave 3      68   72    3x,i5
c        10) ipartysq(3)     party id squeeze in     76   80    3x,i5
c                            wave 3
c        11) ipartidst(4)    party id strength       84   88    3x,i5       
c                            in wave 4
c        12) iparty(4)       party id in wave 4      92   96    3x,i5
c        13) ipartysq(4)     party id squeeze in    100  104    3x,i5
c                            wave 4
c        14) ipartidst(6)    party id strength      108  112    3x,i5
c                            in wave 6
c        15) iparty(6)       party id in wave 6     116  120    3x,i5
c        16) ipartysq(6)     party id squeeze in    124  128    3x,i5
c                            wave 6
c        17) ipartidst(7)    party id strength      132  136    3x,i5   
c                            in wave 7
c        18) iparty(7)       party id in wave 7     140  144    3x,i5
c        19) ipartysq(7)     party id squeeze in    148  152    3x,i5
c        20) ipartidst(8)    party id strength      156  160    3x,i5     
c                            in wave 8
c        21) iparty(8)       party id in wave 8     164  168    3x,i5
c        22) ipartysq(8)     party id squeeze in    172  176    3x,i5
c                            wave 8
c        23) ipartidst(9)    party id strength      180  184    3x,i5
c                            in wave 9
c        24) iparty(9)       party id in wave 9     188  192    3x,i5
c        25) ipartysq(9)     party id squeeze in    196  200    3x,i5
c                            wave 9
c        26) ipartidst(10)   party id strength      204  208    3x,i5  
c                            in wave 10
c        27) iparty(10)      party id in wave 10    212  216    3x,i5
c        28) ipartysq(10)    party id squeeze in    220  224    3x,i5
c                            wave 10
c        29) ipartidst(11)   party id strength      228  232    3x,i5   
c                            in wave 11
c        30) iparty(11)      party id in wave 11    236  240    3x,i5
c        31) ipartysq(11)    party id squeeze in    244  248    3x,i5
c                            wave 11
c        32) ipartidst(12)   party id strength      252  256    3x,i5  
c                            in wave 12
c        33) iparty(12)      party id in wave 12    260  264    3x,i5
c        34) ipartysq(12)    party id squeeze in    268  272    3x,i5
c                            wave 12
c        35) ipartidst(13)   party id strength      276  280    3x,i5  
c                            in wave 13
c        36) iparty(13)      party id in wave 13    284  288    3x,i5
c        37) ipartysq(13)    party id squeeze in    292  296    3x,i5
c                            wave 13
c        38) ipartidst(14)   party id strength      300  304    3x,i5  
c                            in wave 14
c        39) iparty(14)      party id in wave 14    308  312    3x,i5
c        40) ipartysq(14)    party id squeeze in    316  320    3x,i5
c                            wave 14
c        41) ipartidst(15)   party id strength      324  328    3x,i5  
c                            in wave 15
c        42) iparty(15)      party id in wave 15    332  336    3x,i5
c        43) ipartysq(15)    party id squeeze in    340  344    3x,i5
c                            wave 15
c        44) ipartidst(16)   party id strength      348  352    3x,i5  
c                            in wave 16
c        45) iparty(16)      party id in wave 16    356  360    3x,i5
c        46) ipartysq(16)    party id squeeze in    364  368    3x,i5
c                            wave 16
c        47) ipartidst(17)   party id strength      372  376    3x,i5  
c                            in wave 17
c        48) iparty(17)      party id in wave 17    380  384    3x,i5
c        49) ipartysq(17)    party id squeeze in    388  392    3x,i5
c                            wave 17
c        50) ipartidst(18)   party id strength      396  400    3x,i5  
c                            in wave 18
c        51) iparty(18)      party id in wave 18    404  408    3x,i5
c        52) ipartysq(18)    party id squeeze in    412  416    3x,i5
c                            wave 18
c        53) ipartidst(19)   party id strength      420  424    3x,i5  
c                            in wave 19
c        54) iparty(19)      party id in wave 19    428  432    3x,i5
c        55) ipartysq(19)    party id squeeze in    436  440    3x,i5
c                            wave 19
c
      do 8 i=1,nobs
         read(8,106)idt,(ipartidst(j),iparty(j),ipartysq(j),j=1,4),
     &    (ipartidst(j),iparty(j),ipartysq(j),j=6,19)
 106     format(55(3x,i5))
         if(idt.ne.id(i))then
            write(6,107)i,idt,id(i)
 107        format(1x,'mismatch for i,idt,id(i):',3(2x,i5))
            stop
            endif
         if(incluv(i).eq.1)then
            do 9 j=1,ntime
               if(iobsrvv(i,j).eq.1)then
                  if(model.eq.1)then
                     do 10 k=12,16
  10                    xv(i,j,k)=0.
                     if(j.ne.5)then
                        if(iparty(j).eq.1)xv(i,j,12)=1.
                        if(iparty(j).eq.2)xv(i,j,13)=1.
                        if((iparty(j).gt.2).and.(iparty(j).lt.10))
     &                   xv(i,j,14)=1.
                        if(iparty(j).eq.9999)then
                           if(ipartysq(j).eq.1)xv(i,j,12)=1.
                           if(ipartysq(j).eq.2)xv(i,j,13)=1.
                           if((ipartysq(j).gt.2).and.
     &                      (ipartysq(j).lt.10))xv(i,j,14)=1.
                           endif
                        if(ipartidst(j).eq.1)xv(i,j,15)=1.
                        endif
                     if(j.eq.5)xv(i,j,16)=1.
                     endif
                  if(model.eq.2)then
                     do 79 k=14,18
  79                    xv(i,j,k)=0.
                     if(j.ne.5)then
                        if(iparty(j).eq.1)xv(i,j,14)=1.
                        if(iparty(j).eq.2)xv(i,j,15)=1.
                        if((iparty(j).gt.2).and.(iparty(j).lt.10))
     &                   xv(i,j,16)=1.
                        if(iparty(j).eq.9999)then
                           if(ipartysq(j).eq.1)xv(i,j,14)=1.
                           if(ipartysq(j).eq.2)xv(i,j,15)=1.
                           if((ipartysq(j).gt.2).and.
     &                      (ipartysq(j).lt.10))xv(i,j,16)=1.
                           endif
                        if(ipartidst(j).eq.1)xv(i,j,17)=1.
                        if(ipartidst(j).eq.2)xv(i,j,18)=1.
                        endif
                     if(j.eq.5)xv(i,j,18)=1.
                     endif
                  endif
   9           continue
            endif
   8     continue
      close(8)
c
c       read data from aleciadata3.
c
      write(6,127)
 127  format(1x,'read aleciadata3c.prn')
      open(unit=8,file='data/aleciadata3c.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                     4    8    3x,i5 
c         2) itrustmp(1)     trust mps in wave 1     12   16    3x,i5
c         3) iduvot(1)       duty to vote in         44   48   27x,i5
c                            wave 1
c         4) ieffic(1,1)     efficacy: understand    60   64   11x,i5
c                            in wave 1
c         5) ieffic(1,2)     efficacy: pol care      84   88   19x,i5
c                            in wave 1
c         6) isatuk(1)       satisfied w/ uk        108  112   19x,i5
c                            democracy in wave 1
c         7) isateu(1)       satisfied w/ eu        140  144   27x,i5
c                            democracy in wave 1
c         8) itrustmp(2)     trust mps in wave 2    156  160   11x,i5
c         9) ieffic(2,1)     efficacy: understand   188  192   27x,i5
c                            in wave 2
c        10) ieffic(2,2)     efficacy: pol care     212  216   19x,i5
c                            in wave 2
c        11) iduvot(2)       duty to vote in        244  248   27x,i5
c                            wave 2
c        12) isatuk(2)       satisfied w/ uk        268  272   19x,i5
c                            democracy in wave 2
c        13) isateu(2)       satisfied w/ eu        300  304   27x,i5
c                            democracy in wave 2
c        14) itrustmp(3)     trust mps in wave 3    316  320   11x,i5
c        15) ieffic(3,1)     efficacy: understand   332  336   11x,i5
c                            in wave 3
c        16) ieffic(3,2)     efficacy: pol care     356  360   19x,i5
c                            in wave 3
c        17) iduvot(3)       duty to vote in        380  384   19x,i5
c                            wave 3
c        19) isatuk(3)       satisfied w/ uk        404  408   19x,i5
c                            democracy in wave 3
c        20) isateu(3)       satisfied w/ eu        436  440   27x,i5
c                            democracy in wave 3
c        21) itrustmp(4)     trust mps in wave 4    452  456   11x,i5
c        22) ieffic(4,1)     efficacy: understand   468  472   11x,i5
c                            in wave 4
c        23) ieffic(4,2)     efficacy: pol care     492  496   19x,i5
c                            in wave 4
c        24) iduvot(4)       duty to vote in        548  552   51x,i5
c                            wave 4
c
      do 19 i=1,nobs
         do 27 j=1,4
  27        iduvot(j)=-99
         read(8,111)idt,(itrustmp(j),(ieffic(j,k),k=1,2),iduvot(j),
     &    isatuk(j),isateu(j),j=1,3),itrustmp(4),(ieffic(4,j),j=1,2),
     &    iduvot(4)
 111     format(2(3x,i5),27x,i5,11x,i5,2(19x,i5),27x,i5,11x,i5,
     &    2(27x,i5,19x,i5),27x,i5,2(11x,i5),3(19x,i5),27x,i5,2(11x,i5),
     &    19x,i5,51x,i5)
         if(idt.ne.id(i))then
            write(6,107)i,idt,id(i)
            stop
            endif
         if(incluv(i).eq.1)then
            if(model.eq.1)then
               do 20 j=1,4
                  do 21 k=17,19
  21                 xv(i,j,k)=0.
                  if(iobsrvv(i,j).eq.1)then
                     if(j.lt.4)then
                        if((isatuk(j).eq.3).or.(isatuk(j).eq.4))
     &                   xv(i,j,17)=1.
                        if((isateu(j).eq.3).or.(isateu(j).eq.4))
     &                   xv(i,j,18)=1.
                        if(isatuk(j).lt.0)xv(i,j,17)=-99.
                        if(isateu(j).lt.0)xv(i,j,18)=-99.
                        endif
                     if((itrustmp(j).ge.4).and.(itrustmp(j).le.7))
     &                xv(i,j,19)=1.
                     iyv(i,1,j)=iduvot(j)
                     endif
  20              continue
               endif
            if(model.eq.2)then
               do 80 j=1,4
                  iyv(i,2,j)=0
                  if(iobsrvv(i,j).eq.1)then
                     iyv(i,1,j)=iduvot(j)
                     if(j.lt.4)then
                        if((isatuk(j).ge.1).and.(isatuk(j).le.4))then
                           iyv(i,2,j)=isatuk(j)
                           endif
                        if((isatuk(j).lt.0).or.(isatuk(j).gt.4))
     &                   iyv(i,2,j)=-99
                        endif
                     if((itrustmp(j).ge.1).and.(itrustmp(j).le.7))then
                        iyv(i,3,j)=itrustmp(j)
                        endif
                     iyv(i,1,j)=iduvot(j)
                     do 81 k=1,2
                        kp3=k+3
                        iyv(i,kp3,j)=0
                        if((ieffic(j,k).ge.1).and.(ieffic(j,k).le.5))
     &                   then
                           if(k.eq.1)iyv(i,kp3,j)=ieffic(j,k)
                           if(k.eq.2)iyv(i,kp3,j)=6-ieffic(j,k)
                           endif
                        if((ieffic(j,k).lt.1).or.(ieffic(j,k).gt.5))
     &                   iyv(i,kp3,j)=-99
  81                    continue
                     endif
  80              continue
               endif
            endif
  19     continue
      close(8)
c
c       read data from aleciadata4.
c
      write(6,128)
 128  format(1x,'read aleciadata4c.prn')
      open(unit=8,file='data/aleciadata4c.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                     4    8    3x,i5 
c                                                             512x
c         2) isatuk(4)       satisfied w/ uk        524  528    3x,i5
c                            democracy in wave 4
c         3) isateu(4)       satisfied w/ eu        556  560   27x,i5
c                            democracy in wave 4
c
      do 22 i=1,nobs
         read(8,112)idt,isatuk(4),isateu(4)
 112     format(3x,i5,515x,i5,27x,i5)
         if(idt.ne.id(i))then
            write(6,107)i,idt,id(i)
            stop
            endif
         if(incluv(i).eq.1)then
            if(iobsrvv(i,4).eq.1)then
               if(model.eq.1)then
                  xv(i,4,17)=0.
                  xv(i,4,18)=0.
                  if((isatuk(4).eq.3).or.(isatuk(4).eq.4))xv(i,4,17)=1.
                  if((isateu(4).eq.3).or.(isateu(4).eq.4))xv(i,4,18)=1.
                  if(isatuk(4).lt.0)xv(i,4,17)=-99.
                  if(isateu(4).lt.0)xv(i,4,18)=-99.
                  endif
               if(model.eq.2)then
                  iyv(i,2,4)=0
                  if((isatuk(4).ge.1).and.(isatuk(4).le.4))then
                     iyv(i,2,4)=isatuk(4)
                     endif
                  if((isatuk(4).lt.1).or.(isatuk(4).gt.4))iyv(i,2,4)=
     &             -99
                  endif
               endif
            endif
  22     continue
      close(8)
c
c       read data from aleciadata5.
c
      write(6,129)
 129  format(1x,'read aleciadata5c.prn')
      open(unit=8,file='data/aleciadata5c.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                     4    8    3x,i5 
c         2) itrustmp(6)     trust mps in wave 6     36   40   27x,i5
c         3) ieffic(6,1)     efficacy: understand    44   48    3x,i5
c                            in wave 6
c         4) ieffic(6,2)     efficacy: pol care      68   72   19x,i5
c                            in wave 6
c         5) iduvot(6)       duty to vote in        116  120   43x,i5
c                            wave 6
c         6) isatuk(6)       satisfied w/ uk        140  144   19x,i5
c                            democracy in wave 6
c         7) isateu(6)       satisfied w/ eu        172  176   27x,i5
c                            democracy in wave 6
c         8) itrustmp(7)     trust mps in wave 7    204  208   27x,i5
c         9) ieffic(7,1)     efficacy: understand   220  224   11x,i5
c                            in wave 7
c        10) ieffic(7,2)     efficacy: pol care     244  248   19x,i5
c                            in wave 7
c        11) isatuk(7)       satisfied w/ uk        276  280   27x,i5
c                            democracy in wave 7
c        12) isateu(7)       satisfied w/ eu        308  312   27x,i5
c                            democracy in wave 7
c        13) iduvot(7)       duty to vote in        316  320    3x,i5
c                            wave 7
c        14) ieffic(8,1)     efficacy: understand   412  416   91x,i5
c                            in wave 8
c        15) ieffic(8,2)     efficacy: pol care     420  424    3x,i5
c                            in wave 8
c        16) isatuk(8)       satisfied w/ uk        436  440   11x,i5
c                            democracy in wave 8
c        17) isateu(8)       satisfied w/ eu        468  472   27x,i5
c                            democracy in wave 8
c        18) iduvot(8)       duty to vote in        508  512   35x,i5
c                            wave 8
c        19) itrustmp(8)     trust mps in wave 8    532  536   27x,i5
c
      do 23 i=1,nobs
         do 28 j=6,8
  28        iduvot(j)=-99
         read(8,113)idt,itrustmp(6),(ieffic(6,j),j=1,2),iduvot(6),
     &    isatuk(6),isateu(6),itrustmp(7),(ieffic(7,j),j=1,2),
     &    isatuk(7),isateu(7),iduvot(7),(ieffic(8,j),j=1,2),isatuk(8),
     &    isateu(8),iduvot(8),itrustmp(8)
 113     format(3x,i5,27x,i5,3x,i5,19x,i5,43x,i5,19x,i5,2(27x,i5),
     &    11x,i5,19x,i5,2(27x,i5),3x,i5,91x,i5,3x,i5,11x,i5,27x,i5,
     &    35x,i5,27x,i5)
         if(idt.ne.id(i))then
            write(6,107)i,idt,id(i)
            stop
            endif
         if(incluv(i).eq.1)then
            do 24 j=6,8
               if(iobsrvv(i,j).eq.1)then
                  if(model.eq.1)then
                     xv(i,j,19)=-99.
                     if(itrustmp(j).gt.0)then
                        xv(i,j,19)=0.
                        if((itrustmp(j).ge.4).and.(itrustmp(j).le.7))
     &                   xv(i,j,19)=1.
                        endif
                     xv(i,j,17)=0.
                     xv(i,j,18)=0.
                     if((isatuk(j).eq.3).or.(isatuk(j).eq.4))
     &                xv(i,j,17)=1.
                     if((isateu(j).eq.3).or.(isateu(j).eq.4))
     &                xv(i,j,18)=1.
                     if(isatuk(j).lt.0)xv(i,j,17)=-99.
                     if(isateu(j).lt.0)xv(i,j,18)=-99.
                     iyv(i,1,j)=iduvot(j)
                     endif
                  if(model.eq.2)then
                     iyv(i,1,j)=iduvot(j)
                     iyv(i,2,j)=-99
                     if((isatuk(j).ge.1).and.(isatuk(j).le.4))
     &                iyv(i,2,j)=isatuk(j)
                     iyv(i,3,j)=-99
                     if((itrustmp(j).ge.1).and.(itrustmp(j).le.7))
     &                iyv(i,3,j)=itrustmp(j)
                     do 82 k=1,2
                        kp3=k+3
                        iyv(i,kp3,j)=-99
                        if((ieffic(j,k).ge.1).and.(ieffic(j,k).le.5))
     &                   then
                            if(k.eq.1)iyv(i,kp3,j)=ieffic(j,k)
                            if(k.eq.2)iyv(i,kp3,j)=6-ieffic(j,k)
                            endif
  82                    continue
                     endif
                  endif
  24           continue
            endif
  23     continue
      close(8)
c
c       read data from aleciadata6.
c
      write(6,130)
 130  format(1x,'read aleciadata6c.prn')
      open(unit=8,file='data/aleciadata6c.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                     4    8    3x,i5
c         2) itrustmp(9)     trust mps in wave 9     36   40   27x,i5
c         3) isatuk(9)       satisfied w/ uk         44   48    3x,i5
c                            democracy in wave 9
c         4) isateu(9)       satisfied w/ eu         76   80   27x,i5
c                            democracy in wave 9
c         5) ieffic(9,1)     efficacy: understand    84   88    3x,i5
c                            in wave 9
c         6) ieffic(9,2)     efficacy: pol care     108  112   19x,i5
c                            in wave 9
c         7) itrustmp(10)    trust mps in wave 10   148  152   35x,i5
c         8) ieffic(10,1)    efficacy: understand   156  160    3x,i5
c                            in wave 10
c         9) ieffic(10,2)    efficacy: pol care     180  184   19x,i5
c                            in wave 10
c        10) isatuk(10)      satisfied w/ uk        196  200   11x,i5
c                            democracy in wave 10
c        11) isateu(10)      satisfied w/ eu        228  232   27x,i5
c                            democracy in wave 10
c        12) ieffic(11,1)    efficacy: understand   276  280   43x,i5
c                            in wave 11
c        13) ieffic(11,2)    efficacy: pol care     300  304   19x,i5
c                            in wave 11
c        14) isatuk(11)      satisfied w/ uk        324  328   19x,i5
c                            democracy in wave 11
c        15) itrustmp(12)    trust mps in wave 12   372  376   43x,i5
c        16) iduvot(12)      duty to vote in        404  408   27x,i5
c                            wave 12
c        17) isatuk(13)      satisfied w/ uk        444  448   59x,i5
c                            democracy in wave 13
c        18) iduvot(14)      duty to vote in        548  552   99x,i5
c                            wave 14
c
      do 25 i=1,nobs
         do 29 j=9,14
            iduvot(j)=-99
            itrustmp(j)=-99
            isatuk(j)=-99
  29        isateu(j)=-99
         read(8,114)idt,itrustmp(9),isatuk(9),isateu(9),
     &    (ieffic(9,j),j=1,2),itrustmp(10),(ieffic(10,j),j=1,2),
     &    isatuk(10),isateu(10),(ieffic(11,j),j=1,2),isatuk(11),
     &    itrustmp(12),iduvot(12),isatuk(13),iduvot(14)
 114     format(3x,i5,27x,i5,3x,i5,27x,i5,3x,i5,19x,i5,35x,i5,3x,i5,
     &    19x,i5,11x,i5,27x,i5,43x,i5,2(19x,i5),43x,i5,27x,i5,59x,i5,
     &    99x,i5)
         if(idt.ne.id(i))then
            write(6,107)i,idt,id(i)
            stop
            endif
         if(incluv(i).eq.1)then
            do 26 j=9,14
               if(iobsrvv(i,j).eq.1)then
                  if(model.eq.1)then
                     xv(i,j,19)=-99.
                     if(itrustmp(j).gt.0)then
                        xv(i,j,19)=0.
                        if((itrustmp(j).ge.4).and.(itrustmp(j).le.7))
     &                   xv(i,j,19)=1.
                        endif
                     xv(i,j,17)=0.
                     xv(i,j,18)=0.
                     if((isatuk(j).eq.3).or.(isatuk(j).eq.4))
     &                xv(i,j,17)=1.
                     if((isateu(j).eq.3).or.(isateu(j).eq.4))
     &                xv(i,j,18)=1.
                     if(isatuk(j).lt.0)xv(i,j,17)=-99.
                     if(isateu(j).lt.0)xv(i,j,18)=-99.
                     iyv(i,1,j)=iduvot(j)
                     endif
                  if(model.eq.2)then
                     iyv(i,1,j)=iduvot(j)
                     do 83 k=2,5
  83                    iyv(i,k,j)=-99
                     if((isatuk(j).ge.1).and.(isatuk(j).le.4))
     &                iyv(i,2,j)=isatuk(j)
                     if((itrustmp(j).ge.1).and.(itrustmp(j).le.7))
     &                iyv(i,3,j)=itrustmp(j)
                     do 84 k=1,2
                        if((ieffic(j,k).ge.1).and.(ieffic(j,k).le.5))
     &                   then
                           if(k.eq.1)iyv(i,k+3,j)=ieffic(j,k)
                           if(k.eq.2)iyv(i,k+3,j)=6-ieffic(j,k)
                           endif
  84                    continue
                     endif
                  endif
  26           continue
            endif
  25     continue
      close(8)
c
c       read data from aleciadata7.
c
      write(6,131)
 131  format(1x,'read aleciadata7c.prn')
      open(unit=8,file='data/aleciadata7c.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                     4    8    3x,i5 
c         2) itrustmp(15)    trust mps in wave 15    52   56   43x,i5
c         3) isatuk(15)      satisfied w/ uk         76   80   19x,i5 
c                            democracy in wave 15
c         4) iduvot(15)      duty to vote in        108  112   27x,i5
c                            wave 15
c         5) ieffic(15,1)    efficacy: understand   124  128   11x,i5
c                            in wave 15
c         6) ieffic(15,2)    efficacy: pol care     148  152   19x,i5
c                            in wave 15
c         7) isatuk(16)      satisfied w/ uk        244  248   91x,i5
c                            democracy in wave 16
c         8) ieffic(16,1)    efficacy: understand   276  280   27x,i5
c                            in wave 16
c         9) ieffic(16,2)    efficacy: pol care     300  304   19x,i5
c                            in wave 16
c        10) itrustmp(16)    trust mps in wave 16   316  320   11x,i5
c        11) itrustmp(17)    trust mps in wave 17   380  384   59x,i5
c        12) ieffic(17,1)    efficacy: understand   412  416   27x,i5
c                            in wave 17
c        13) ieffic(17,2)    efficacy: pol care     436  440   19x,i5
c                            in wave 17
c        14) iduvot(17)      duty to vote in        452  456   11x,i5  
c                            wave 17
c        15) isatuk(17)      satisfied w/ uk        468  472   11x,i5
c                            democracy in wave 17
c        16) itrustmp(18)    trust mps in wave 18   532  536   59x,i5
c        17) iduvot(18)      duty to vote in        556  560   19x,i5
c                            wave 18
c        18) isatuk(19)      satisfied w/ uk        612  616   51x,i5
c                            democracy in wave 19
c        19) itrustmp(19)    trust mps in wave 19   644  648   27x,i5
c
      do 30 i=1,nobs
         do 31 j=15,19
            iduvot(j)=-99
            itrustmp(j)=-99
            isatuk(j)=-99
  31        isateu(j)=-99
         read(8,115)idt,itrustmp(15),isatuk(15),iduvot(15),
     &    (ieffic(15,j),j=1,2),isatuk(16),(ieffic(16,j),j=1,2),
     &    (itrustmp(j),j=16,17),(ieffic(17,j),j=1,2),iduvot(17),
     &    isatuk(17),itrustmp(18),iduvot(18),isatuk(19),
     &    itrustmp(19)
 115     format(3x,i5,43x,i5,19x,i5,27x,i5,11x,i5,19x,i5,91x,i5,27x,i5,
     &    19x,i5,11x,i5,59x,i5,27x,i5,19x,i5,2(11x,i5),59x,i5,19x,i5,
     &    51x,i5,27x,i5)
         if(idt.ne.id(i))then
            write(6,107)i,idt,id(i)
            stop
            endif
         if(incluv(i).eq.1)then
            do 32 j=15,19
               if(iobsrvv(i,j).eq.1)then
                  if(model.eq.1)then
                     xv(i,j,19)=-99.
                     if(itrustmp(j).gt.0)then
                        xv(i,j,19)=0.
                        if((itrustmp(j).ge.4).and.(itrustmp(j).le.7))
     &                   xv(i,j,19)=1.
                        endif
                     xv(i,j,17)=0.
                     xv(i,j,18)=0.
                     if((isatuk(j).eq.3).or.(isatuk(j).eq.4))
     &                xv(i,j,17)=1.
                     if((isateu(j).eq.3).or.(isateu(j).eq.4))
     &                xv(i,j,18)=1.
                     if(isatuk(j).lt.0)xv(i,j,17)=-99.
                     if(isateu(j).lt.0)xv(i,j,18)=-99.
                     iyv(i,1,j)=iduvot(j)
                     endif
                  if(model.eq.2)then
                     iyv(i,1,j)=iduvot(j)
                     do 85 k=2,5
  85                    iyv(i,k,j)=-99
                     if((isatuk(j).ge.1).and.(isatuk(j).le.4))
     &                iyv(i,2,j)=isatuk(j)
                     if((itrustmp(j).ge.1).and.(itrustmp(j).le.7))
     &                iyv(i,3,j)=itrustmp(j)
                     do 86 k=1,2
                        if((ieffic(j,k).ge.1).and.(ieffic(j,k).le.5))
     &                   then
                           if(k.eq.1)iyv(i,k+3,j)=ieffic(j,k)
                           if(k.eq.2)iyv(i,k+3,j)=6-ieffic(j,k)
                           endif
  86                    continue
                     endif
                  endif
  32           continue
            endif
  30     continue
      close(8)
c
c       exclude observations with not enough observed dependent 
c       variables.
c
      ideldep=0
      do 78 i=1,ntime
  78     inumdens(i)=0
      do 75 i=1,nobs
         if(incluv(i).eq.1)then
            istart=istartv(i)
            iend=iendv(i)
            inum=0
            do 76 j=istart,iend
               if((iobsrvv(i,j).eq.1).and.(iyv(i,1,j).ge.0))inum=inum+1
  76           continue
            if(inum.le.ifcond1)then
               ideldep=ideldep+1
               incluv(i)=0
               endif
            if(inum.gt.ifcond1)inumdens(inum)=inumdens(inum)+1
            endif
  75     continue
      write(6,152)ideldep
 152  format(1x,'# obs deleted for no dependent variables: ',i5)
      write(6,153)
 153  format(1x,'frequency of # dependent variables')
      write(6,154)
 154  format(1x,'# dep vars',2x,'frequency')
      do 77 i=1,ntime
         if(inumdens(i).gt.0)write(6,155)i,inumdens(i)
 155     format(1x,i2,10x,i5)
  77     continue
c
c       read data from aleciadata8.
c
      write(6,132)
 132  format(1x,'read aleciadata8c.prn')
      open(unit=8,file='data/aleciadata8c.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                     4    8    3x,i5 
c         2) icountry(1)     country in wave 1       14   16    5x,i3
c         3) icountry(2)     country in wave 2       22   24    5x,i3
c         4) icountry(3)     country in wave 3       30   32    5x,i3
c         5) icountry(4)     country in wave 4       38   40    5x,i3
c         6) icountry(5)     country in wave 5       46   48    5x,i3
c         7) icountry(6)     country in wave 6       54   56    5x,i3
c         8) icountry(7)     country in wave 7       62   64    5x,i3
c         9) icountry(8)     country in wave 8       70   72    5x,i3
c        10) icountry(9)     country in wave 9       78   80    5x,i3
c        11) icountry(10)    country in wave 10      86   88    5x,i3
c        12) icountry(11)    country in wave 11      94   96    5x,i3
c        13) icountry(12)    country in wave 12     102  104    5x,i3
c        14) icountry(13)    country in wave 13     110  112    5x,i3
c        15) icountry(14)    country in wave 14     118  120    5x,i3
c        16) icountry(15)    country in wave 15     126  128    5x,i3
c        17) icountry(16)    country in wave 16     134  136    5x,i3
c        18) icountry(17)    country in wave 17     134  136    5x,i3
c        19) icountry(18)    country in wave 18     142  144    5x,i3
c        20) icountry(19)    country in wave 19     150  152    5x,i3
c                                                               8x
c        21) ifemale         female + 1             174  176    5x,i3
c        22) iage            age                    182  184    5x,i3  
c        23) income(1)       bracketed income in    198  200   13x,i3
c                            wave 1
c        24) irelig(1)       religion in wave 1     230  232   29x,i3
c        25) income(2)       bracketed income in    262  264   29x,i3
c                            wave 2
c        26) irelig(2)       religion in wave 2     294  296   29x,i3
c        27) income(3)       bracketed income in    326  328   29x,i3
c                            wave 3
c        28) irelig(3)       religion in wave 3     358  360   29x,i3
c        29) income(4)       bracketed income in    390  392   29x,i3
c                            wave 4
c        30) irelig(4)       religion in wave 4     422  424   29x,i3
c        31) income(5)       bracketed income in    454  456   29x,i3
c                            wave 5
c        32) irelig(5)       religion in wave 5     486  488   29x,i3
c        33) income(6)       bracketed income in    518  520   29x,i3
c                            wave 6
c        34) irelig(6)       religion in wave 6     550  552   29x,i3
c        35) income(7)       bracketed income in    582  584   29x,i3
c                            wave 7
c        36) irelig(7)       religion in wave 7     614  616   29x,i3
c
      ncause=nx+1
      do 39 i=1,ntime
         do 40 j=1,ncause
  40        icause(i,j)=0
  39     continue
      nicodens=0
      do 11 i=1,nobs
         read(8,108)idt,(icountry(j),j=1,19),ifemale,iage,(income(j),
     &    irelig(j),j=1,7)
 108     format(3x,i5,19(5x,i3),8x,2(5x,i3),13x,i3,13(29x,i3))
         if(idt.ne.id(i))then
            write(6,107)i,idt,id(i)
            stop
            endif
         if(model.eq.2)then
            do 88 j=1,3
  88           icodent(j)=0
            do 87 j=1,ntime
               if((icountry(j).ge.1).and.(icountry(j).le.3))
     &          icodent(icountry(j))=icodent(icountry(j))+1
  87           continue
            icomax=-5
            icon=-5
            do 89 j=1,3
               if(icodent(j).gt.icomax)then
                  icon=j
                  icomax=icodent(j)
                  endif
  89           continue
            do 90 j=1,ntime
               do 91 k=19,20
  91              xv(i,j,k)=0.
               if(icon.lt.1)then
                  write(6,102)i
 102              format(1x,'no country for i = ',i8)
                  stop
                  endif
               if((icon.ge.2).and.(icon.le.3))xv(i,j,icon+17)=1.
  90           continue
            if(nicodens.gt.0)then
               do 92 j=1,nicodens
                  do 93 k=1,3
                     if(icodent(k).ne.icodens(j,k))goto 92
  93                 continue
                  icodens(j,4)=icodens(j,4)+1
                  goto 94
  92              continue
               endif
            nicodens=nicodens+1
            do 95 j=1,3
  95           icodens(nicodens,j)=icodent(j)
            icodens(nicodens,4)=1
  94        continue
            endif
         if(ifemale.lt.0)then
            incluv(i)=0
            do 54 j=1,19
               if(iobsrvv(i,j).eq.1)icause(j,3)=icause(j,3)+1
  54           continue
            endif
         if(incluv(i).eq.1)then
            ifemale=ifemale-1
            if(model.eq.1)then
               if(iage.lt.0)iaget=8
               if(iage.gt.68)iaget=4
               if((iage.le.68).and.(iage.ge.50))iaget=5
               if((iage.le.51).and.(iage.ge.34))iaget=6
               if((iage.le.33).and.(iage.ge.18))iaget=7
               endif
            if(model.eq.2)then
               iaget=0
               if(iage.lt.0)iaget=9
               if(iage.gt.68)iaget=6
               if((iage.le.68).and.(iage.ge.50))iaget=7
               if((iage.le.51).and.(iage.ge.34))iaget=8
               endif
            do 55 j=1,ntime
               if(iobsrvv(i,j).eq.1)then
                  if(model.eq.1)then
                     xv(i,j,2)=ifemale
                     do 13 k=4,8
  13                    xv(i,j,k)=0.
                     xv(i,j,iaget)=1.
                     endif
                  if(model.eq.2)then
                     xv(i,j,2)=ifemale
                     do 96 k=6,8
  96                    xv(i,j,k)=0.
                     if(iaget.gt.0)xv(i,j,iaget)=1.
                     endif
                  endif
  55           continue
            do 12 j=1,7
               if(iobsrvv(i,j).eq.1)then
                  if(model.eq.1)then
                     do 56 k=9,11
  56                    xv(i,j,k)=0.
                     if((income(j).ge.1).and.(income(j).le.15))
     &                xv(i,j,9)=xinc(income(j))
                     if((income(j).lt.1).or.(income(j).gt.15))
     &                xv(i,j,10)=1.
                     if((irelig(j).gt.1).and.(irelig(j).le.15))
     &                xv(i,j,11)=1.
                      endif
                  if(model.eq.2)then
                     do 97 k=10,12
  97                    xv(i,j,k)=0.
                     if((income(j).ge.1).and.(income(j).le.5))
     &                xv(i,j,10)=1.
                     if((income(j).ge.6).and.(income(j).le.9))
     &                xv(i,j,11)=1.
                     if((income(j).lt.1).or.(income(j).gt.15))
     &                xv(i,j,12)=1.
                     if((irelig(j).gt.1).and.(irelig(j).le.15))
     &                xv(i,j,13)=1.
                     endif
                  endif
  12           continue
            endif
  11     continue
      write(6,117)
      write(6,119)
 119  format(1x,'distribution of county combinations')
      write(6,158)
 158  format(1x,'# obs',5x,'england',2x,'scotland',2x,'wales')
      nicodup=0
      do 98 i=1,nicodens
         nicon=-5
         icomax=-5
         do 99 j=1,nicodens
            if(icodens(j,4).gt.icomax)then
               nicon=j
               icomax=icodens(j,4)
               endif
  99        continue
         if(icomax.gt.0)then
            write(6,159)icodens(nicon,4),(icodens(nicon,j),j=1,3)
 159        format(1x,i8,2x,i2,7x,i2,8x,i2)
            idup=0
            do 220 j=1,3
               if(icodens(nicon,j).gt.0)idup=idup+1
 220           continue
            if(idup.gt.1)nicodup=nicodup+1
            icodens(nicon,4)=-10
            goto 98
            endif
  98     continue
      write(6,169)nicodup
 169  format(1x,'# obs w/ > 1 country: ',i5)
      close(8)
c
c       read better age data.
c
      call datextraget(1)
c
c       read data from aleciadata9.
c
      write(6,133)
 133  format(1x,'read aleciadata9c.prn')
      open(unit=8,file='data/aleciadata9c.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                     4    8    3x,i5 
c         2) income(8)       bracketed income in     22   24   13x,i3
c                            wave 8
c         3) irelig(8)       religion in wave 8      54   56   29x,i3
c         4) income(9)       bracketed income in     86   88   29x,i3
c                            wave 9
c         5) irelig(9)       religion in wave 9     118  120   29x,i3
c         6) income(10)      bracketed income in    150  152   29x,i3
c                            wave 10
c         7) irelig(10)      religion in wave 10    182  184   29x,i3
c         8) income(11)      bracketed income in    214  216   29x,i3
c                            wave 11
c         9) irelig(11)      religion in wave 11    246  248   29x,i3
c        10) income(12)      bracketed income in    278  280   29x,i3
c                            wave 12
c        11) irelig(12)      religion in wave 12    310  312   29x,i3
c        12) income(13)      bracketed income in    342  344   29x,i3
c                            wave 13
c        13) irelig(13)      religion in wave 13    374  376   29x,i3
c        14) income(14)      bracketed income in    406  408   29x,i3
c                            wave 14
c        15) irelig(14)      religion in wave 14    438  440   29x,i3
c        16) income(15)      bracketed income in    470  472   29x,i3
c                            wave 15
c        17) irelig(15)      religion in wave 15    502  504   29x,i3
c        18) income(16)      bracketed income in    534  536   29x,i3
c                            wave 16
c        19) irelig(16)      religion in wave 16    566  568   29x,i3
c
      do 14 i=1,nobs
         read(8,109)idt,(income(j),irelig(j),j=8,16)
 109     format(3x,i5,13x,i3,17(29x,i3))
         if(idt.ne.id(i))then
            write(6,107)i,idt,id(i)
            stop
            endif
         if(incluv(i).eq.1)then
            do 15 j=8,16
               if(iobsrvv(i,j).eq.1)then
                  if(model.eq.1)then
                     do 57 k=9,11
  57                    xv(i,j,k)=0.
                     if((income(j).ge.1).and.(income(j).le.15))
     &                xv(i,j,9)=xinc(income(j))
                     if((income(j).lt.1).or.(income(j).gt.15))
     &                xv(i,j,10)=1.
                     if((irelig(j).gt.1).and.(irelig(j).le.15))
     &                xv(i,j,11)=1.
                     endif
                  if(model.eq.2)then
                     do 200 k=10,12
 200                    xv(i,j,k)=0.
                     if((income(j).ge.1).and.(income(j).le.5))
     &                xv(i,j,10)=1.
                     if((income(j).ge.6).and.(income(j).le.9))
     &                xv(i,j,11)=1.
                     if((income(j).lt.1).or.(income(j).gt.15))
     &                xv(i,j,12)=1.
                     if((irelig(j).gt.1).and.(irelig(j).le.15))
     &                xv(i,j,13)=1.
                     endif
                  endif
  15           continue
            endif
  14     continue
      close(8)
c
c       read data from aleciadata10.
c
      write(6,134)
 134  format(1x,'read aleciadata10c.prn')
      open(unit=8,file='data/aleciadata10c.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                     4    8    3x,i5 
c         2) income(17)      bracketed income in     22   24   13x,i3
c                            wave 17
c         3) irelig(17)      religion in wave 17     54   56   29x,i3
c         4) income(18)      bracketed income in     86   88   13x,i3
c                            wave 18
c         5) irelig(18)      religion in wave 18    118  120   29x,i3
c         6) income(19)      bracketed income in    150  152   13x,i3
c                            wave 19
c         7) irelig(19)      religion in wave 19    182  184   29x,i3
c                                                              24x
c         8) ieduc(1)        education in wave 1    213  216    4x,i4
c         9) ieduc(2)        education in wave 2    221  224    4x,i4
c        10) ieduc(3)        education in wave 3    229  232    4x,i4
c        11) ieduc(4)        education in wave 4    237  240    4x,i4
c        12) ieduc(5)        education in wave 5    245  248    4x,i4
c        13) ieduc(6)        education in wave 6    253  256    4x,i4
c        14) ieduc(7)        education in wave 7    261  264    4x,i4
c        15) ieduc(8)        education in wave 8    269  272    4x,i4
c        16) ieduc(9)        education in wave 9    277  280    4x,i4
c        17) ieduc(10)       education in wave 10   285  288    4x,i4
c        18) ieduc(11)       education in wave 11   293  296    4x,i4
c        19) ieduc(12)       education in wave 12   301  304    4x,i4
c        20) ieduc(13)       education in wave 13   309  312    4x,i4
c        21) ieduc(14)       education in wave 14   317  320    4x,i4
c        22) ieduc(15)       education in wave 15   325  328    4x,i4
c        23) ieduc(16)       education in wave 16   333  336    4x,i4
c        24) ieduc(17)       education in wave 17   341  344    4x,i4
c        25) ieduc(18)       education in wave 18   349  352    4x,i4
c        26) ieduc(19)       education in wave 19   357  360    4x,i4
c
      do 16 i=1,nobs
         read(8,110)idt,(income(j),irelig(j),j=17,19),(ieduc(j),j=1,19)
 110     format(3x,i5,13x,i3,5(29x,i3),24x,19(4x,i4))
         if(idt.ne.id(i))then
            write(6,107)i,idt,id(i)
            stop
            endif
         if(incluv(i).eq.1)then
            do 17 j=17,19
               if(iobsrvv(i,j).eq.1)then
                  if(model.eq.1)then
                     do 58 k=9,11
  58                    xv(i,j,k)=0.
                     if((income(j).ge.1).and.(income(j).le.15))
     &                xv(i,j,9)=xinc(income(j))
                     if((income(j).lt.1).or.(income(j).gt.15))
     &                xv(i,j,10)=1.
                     if((irelig(j).gt.1).and.(irelig(j).le.15))
     &                xv(i,j,11)=1.
                     endif
                  if(model.eq.2)then
                     do 201 k=10,12
 201                    xv(i,j,k)=0.
                     if((income(j).ge.1).and.(income(j).le.5))
     &                xv(i,j,10)=1.
                     if((income(j).ge.6).and.(income(j).le.9))
     &                xv(i,j,11)=1.
                     if((income(j).lt.1).or.(income(j).gt.15))
     &                xv(i,j,12)=1.
                     if((irelig(j).gt.1).and.(irelig(j).le.15))
     &                xv(i,j,13)=1.
                     endif
                  endif
  17           continue
            if(model.eq.1)then
               do 18 j=1,ntime
                  xv(i,j,3)=0.
                  if(iobsrvv(i,j).eq.1)then
                     if((ieduc(j).eq.4).or.(ieduc(j).eq.5))xv(i,j,3)=1.
                     endif
  18              continue
               endif
            endif
  16     continue
      close(8)
c
c       read data from aleciaeducationc.
c
      if(model.eq.2)then
         write(6,161)
 161     format(1x,'read aleciaeducationc.prn')
         open(unit=8,file='data/aleciaeducationc.prn')
c
c       read relevant data:
c            variable        definition           start  end    format    
c         1) idt             id #                    1     5    i5
c                                                               152x
c         2) ieducr(1)       education recode in   163   165    5x,i3
c                            wave 1
c         3) ieducr(2)       education recode in   171   173    5x,i3
c                            wave 2
c         4) ieducr(3)       education recode in   179   181    5x,i3
c                            wave 3
c         5) ieducr(4)       education recode in   187   189    5x,i3
c                            wave 4
c         6) ieducr(5)       education recode in   195   197    5x,i3
c                            wave 5
c         7) ieducr(6)       education recode in   203   205    5x,i3
c                            wave 6
c         8) ieducr(7)       education recode in   211   213    5x,i3
c                            wave 7
c         9) ieducr(8)       education recode in   219   221    5x,i3
c                            wave 8
c        10) ieducr(9)       education recode in   227   229    5x,i3
c                            wave 9
c        11) ieducr(10)      education recode in   235   237    5x,i3
c                            wave 10
c        12) ieducr(11)      education recode in   243   245    5x,i3
c                            wave 11
c        13) ieducr(12)      education recode in   251   253    5x,i3
c                            wave 12
c        14) ieducr(13)      education recode in   259   261    5x,i3
c                            wave 13
c        15) ieducr(14)      education recode in   267   269    5x,i3
c                            wave 14
c        16) ieducr(15)      education recode in   275   277    5x,i3
c                            wave 15
c        17) ieducr(16)      education recode in   283   285    5x,i3
c                            wave 16
c        18) ieducr(17)      education recode in   291   293    5x,i3
c                            wave 17
c        19) ieducr(18)      education recode in   299   301    5x,i3
c                            wave 18
c        20) ieducr(19)      education recode in   307   309    5x,i3
c                            wave 19
c
         ieducbord=0
         neducdent=0
         neducbdens=0
         do 202 i=1,nobs
            read(8,162)idt,(ieducr(j),j=1,19)
 162        format(i5,152x,19(5x,i3))
            do 206 j=1,6
 206           ieducdent(j)=0
            ieducmax=-5
            iflag=0
            do 203 j=1,ntime
               if((ieducr(j).ge.0).and.(ieducr(j).le.5))then
                  ieducrjp1=ieducr(j)+1
                  ieducdent(ieducrjp1)=ieducdent(ieducrjp1)+1
                  if(ieducr(j).lt.ieducmax)then
                     if(iflag.eq.0)then
                        iflag=1
                        ieducbord=ieducbord+1
                        if(neducbdens.gt.0)then
                           do 228 k=1,neducbdens
                              if((ieducbdens(k,1).eq.ieducmax).and.
     &                         (ieducbdens(k,2).eq.ieducr(j)))then
                                 ieducbdens(k,3)=ieducbdens(k,3)+1
                                 goto 229
                                 endif
 228                          continue
                           endif
                        neducbdens=neducbdens+1
                        ieducbdens(neducbdens,1)=ieducmax
                        ieducbdens(neducbdens,2)=ieducr(j)
                        ieducbdens(neducbdens,3)=1
 229                    continue
                        endif
                     endif
                  if(ieducr(j).gt.ieducmax)ieducmax=ieducr(j)
                  endif
 203           continue
            if(neducdent.gt.0)then
               do 207 j=1,neducdent
                  do 208 k=1,6
                     if(ieducdent(k).ne.ieducdens(j,k))goto 207
 208                 continue
                  ieducdens(j,7)=ieducdens(j,7)+1
                  goto 209     
 207              continue
               endif
            neducdent=neducdent+1
            do 210 j=1,6
 210           ieducdens(neducdent,j)=ieducdent(j)
            ieducdens(neducdent,7)=1
 209        continue
            do 205 j=1,ntime
               do 204 k=3,5
 204              xv(i,j,k)=0.
               if((ieducmax.eq.2).or.(ieducmax.eq.3))xv(i,j,3)=1.
               if(ieducmax.eq.4)xv(i,j,4)=1.
               if(ieducmax.eq.5)xv(i,j,5)=1.
 205           continue
 202        continue
         write(6,117)
         write(6,163)ieducbord
 163     format(1x,'# obs w/ decreasing changes in education: ',i5)
         if(ieducbord.gt.0)then
            write(6,174)
 174        format(1x,'decreasing patterns')
            write(6,175)
 175        format(1x,'rank',2x,'frequency',2x,'first',2x,'second')
            do 230 i=1,neducbdens
               ieducbmax=-5
               iconeducb=-5
               do 231 j=1,neducbdens
                  if(ieducbdens(j,3).gt.ieducbmax)then
                     ieducbmax=ieducbdens(j,3)
                     iconeducb=j
                     endif
 231              continue
               if(ieducbmax.gt.0)then
                  write(6,173)i,ieducbmax,(ieducbdens(iconeducb,j),
     &             j=1,2)
 173              format(1x,i4,2x,i4,7x,i2,5x,i2)
                  ieducbdens(iconeducb,3)=-10
                  endif
 230           continue
            endif
         write(6,164)
 164     format(1x,'patterns of education')
         write(6,165)(i,i=0,5)
 165     format(18x,'pattern',/,1x,'rank',2x,'frequency',6(2x,i3))
         do 211 i=1,neducdent
            neducmax=-5
            neduc=-5
            do 212 j=1,neducdent
               if(ieducdens(j,7).gt.neducmax)then
                  neduc=j
                  neducmax=ieducdens(j,7)
                  endif
 212           continue
            if(neducmax.gt.0)then
               write(6,166)i,neducmax,(ieducdens(neduc,j),j=1,6)
 166           format(1x,i4,2x,i5,4x,6(2x,i3))
               ieducdens(neduc,7)=-10
               endif
 211        continue
         close(8)
         endif
c
c       operate on missing democracy satisfaction and trust variables.
c
      if((model.eq.1).and.(ifsattr.ge.1))then
         if((ifsattr.ge.1).and.(ifsattr.lt.4))then
            write(6,144)ifsattr
 144        format(1x,'code not written yet for ifsattr = ',i1)
            stop
            endif
         if(ifsattr.eq.4)then
            ntimep1=ntime+1
            ntimep2=ntime+2
            do 71 i=1,2
               do 67 j=1,3
                  do 68 k=1,ntimep2
  68                 interden(i,j,k)=0
  67              continue
  71           continue
            do 59 i=1,nobs
               if(incluv(i).eq.1)then
                  istart=istartv(i)
                  iend=iendv(i)
                  iflag=0
                  do 60 j=istart,iend
                     if(iobsrvv(i,j).eq.1)then
                        do 61 k=1,3
                           kp16=k+16
                           if(xv(i,j,kp16).lt.0.)then
                              if(j.eq.istart)then
                                 interpt=99
                                 goto 63
                                 endif
                              if(j.gt.istart)then
                                 ltt=j
                                 do 62 l=1,ntime
                                    ltt=ltt-1
                                    if(ltt.ge.istart)then
                                       if((xv(i,ltt,kp16).eq.0.).or.
     &                                  (xv(i,ltt,kp16).eq.1.))then
                                          interpt=l
                                          goto 63
                                          endif
                                       endif
  62                                continue
                                 interpt=99
                                 endif
  63                          continue
                              if(j.eq.iend)then
                                 if(interpt.eq.99)interpt=199
                                 if(interpt.lt.ntime)interpt=99
                                 endif
                              if(j.lt.iend)then
                                 ltt=j
                                 do 64 l=1,ntime
                                    ltt=ltt+1
                                    if(ltt.gt.iend)goto 70
                                    if(ltt.le.iend)then
                                       if((xv(i,ltt,kp16).eq.0.).or.
     &                                  (xv(i,ltt,kp16).eq.1.))then
                                          if(interpt.lt.ntime)interpt=
     &                                     interpt+l
                                          goto 65
                                          endif
                                       endif
  64                                continue
  70                             continue
                                 if(interpt.eq.99)interpt=199
                                 if(interpt.lt.ntime)interpt=99
  65                             continue
                                 endif
                              if(iflag.eq.0)then
                                 if(interpt.lt.99)
     &                            interden(1,k,interpt)=
     &                            interden(1,k,interpt)+1
                                 if(interpt.eq.99)
     &                            interden(1,k,ntimep1)=
     &                            interden(1,k,ntimep1)+1
                                 if(interpt.eq.199)
     &                            interden(1,k,ntimep2)=
     &                            interden(1,k,ntimep2)+1
                                 iflag=1
                                 endif
                              if(interpt.lt.99)interden(2,k,interpt)=
     &                         interden(2,k,interpt)+1
                              if(interpt.eq.99)interden(2,k,ntimep1)=
     &                         interden(2,k,ntimep1)+1
                              if(interpt.eq.199)interden(2,k,ntimep2)=
     &                         interden(2,k,ntimep2)+1
                              endif
  61                       continue
                        endif
  60                 continue
                  endif
  59           continue
            write(6,117)
            write(6,145)
 145        format(1x,'interpolation analysis for democracy',
     &       ' satisfaction and trust variables')
            do 72 i=1,2
               ie=i*2
               ib=ie-1
               write(6,149)(aunitan(j),j=ib,ie)
 149           format(1x,'analysis by ',2a8)
               do 69 j=1,3
                  je=(j+16)*2
                  jb=je-1
                  write(6,151)
 151              format(1x,50('.'))
                  write(6,146)(ax1(k),k=jb,je)
 146              format(1x,'analysis for ',2a8)
                  do 66 k=1,ntimep2
                     if(interden(i,j,k).gt.0)then
                        if(k.le.ntime)write(6,147)k,interden(i,j,k)
 147                    format(1x,'interpolation from ',i2,' periods',
     &                   ' apart',3x,i6)
                        if(k.gt.ntime)then
                           kt=k-ntime
                           write(6,148)yesno(kt),interden(i,j,k)
 148                       format(1x,a2,' extrapolation',22x,i6)
                           endif
                        endif
  66                 continue
  69              continue
               write(6,150)
 150           format(1x,50('-'))
  72           continue
            stop
            endif
         endif
c
c       compute descriptive statistics.
c
      do 49 i=1,nxs
         do 33 j=1,4
  33        xmom(i,j)=0.d0
         xmom(i,5)=1.d20
  49     xmom(i,6)=-1.d20
      if(model.eq.1)nyt=1
      if(model.eq.2)nyt=ny
      do 214 i=1,nyt
         do 42 j=1,4
  42        ymom(i,j)=0.d0
         ymom(i,5)=1.d20
 214     ymom(i,6)=-1.d20
      do 216 i=1,nyt
         do 52 j=1,nsegp1s
  52        iydens(i,j)=0
 216     continue
      nbound=0
      do 34 i=1,nobs
         if(incluv(i).eq.1)then
            istart=istartv(i)
            iend=iendv(i)
            if(nbound.gt.0)then
               do 36 j=1,nbound
                  if((istart.eq.ibound(j,1)).and.(iend.eq.ibound(j,2)))
     &             then
                     ibound(j,3)=ibound(j,3)+1
                     goto 37
                     endif
  36              continue
               endif
            nbound=nbound+1
            ibound(nbound,1)=istart
            ibound(nbound,2)=iend
            ibound(nbound,3)=1
  37        continue
            do 43 k=1,nx    
               do 45 l=1,3
  45              xmt(k,l)=0.d0
  43           continue
            do 215 k=1,nyt
               do 6 l=1,3
   6              ymt(k,l)=0.d0
 215           continue
            do 35 j=istart,iend
               if(iobsrvv(i,j).eq.1)then
                  do 38 k=1,nx
                     xt=dble(xv(i,j,k))
                     xmom(k,1)=xmom(k,1)+1.d0
                     xmom(k,2)=xmom(k,2)+xt
                     xmom(k,3)=xmom(k,3)+(xt*xt)
                     if(xmom(k,5).gt.xt)xmom(k,5)=xt
                     if(xmom(k,6).lt.xt)xmom(k,6)=xt
                     xmt(k,1)=xmt(k,1)+1.d0
                     xmt(k,2)=xmt(k,2)+xt
                     xmt(k,3)=xmt(k,3)+(xt*xt)
  38                 continue
                  do 213 k=1,nyt
                     if((iyv(i,k,j).ge.0).and.(iyv(i,k,j).lt.10))then
                        yt=dble(float(iyv(i,k,j)))
                        ymom(k,1)=ymom(k,1)+1.d0
                        ymom(k,2)=ymom(k,2)+yt
                        ymom(k,3)=ymom(k,3)+(yt*yt)
                        if(ymom(k,5).gt.yt)ymom(k,5)=yt
                        if(ymom(k,6).lt.yt)ymom(k,6)=yt
                        ymt(k,1)=ymt(k,1)+1.d0
                        ymt(k,2)=ymt(k,2)+yt
                        ymt(k,3)=ymt(k,3)+(yt*yt)
                        iyvt=iyv(i,k,j)+1
                        iydens(k,iyvt)=iydens(k,iyvt)+1
                        endif
                     if((iyv(i,k,j).lt.0).or.(iyv(i,k,j).gt.9))
     &                icause(j,1)=icause(j,1)+1
 213                 continue
                  do 41 k=1,nx
                     xt=dble(xv(i,j,k))
                     if(xt.lt.0.d0)icause(j,k+1)=icause(j,k+1)+1
  41                 continue
                  endif
  35           continue
            do 217 j=1,nyt
               if(ymt(j,1).gt.0.d0)then
                  ymt(j,2)=ymt(j,2)/ymt(j,1)
                  ymt(j,3)=(ymt(j,3)/ymt(j,1))-(ymt(j,2)*ymt(j,2)) 
                  ymom(j,4)=ymom(j,4)+ymt(j,3)
                  endif
 217           continue
            do 44 j=1,nx
               if(xmt(j,1).gt.0.d0)then
                  xmt(j,2)=xmt(j,2)/xmt(j,1)
                  xmom(j,4)=xmom(j,4)+(xmt(j,1)*((xmt(j,3)/xmt(j,1))-
     &             (xmt(j,2)*xmt(j,2))))
                  endif
  44           continue
            endif
  34     continue
c
c       output descriptive statistics.
c
      write(6,117)
 117  format(1x,50('='))
      write(6,118)
 118  format(1x,'ibound statistics')
      write(6,160)
 160  format(1x,'#',6x,1x,'frequency',2x,'start',2x,'end')
      do 46 i=1,nbound
         ibmax=-89
         do 47 j=1,nbound
            if(ibound(j,3).gt.ibmax)then
               ibmax=ibound(j,3)
               ibi=j
               endif
  47        continue
         write(6,116)i,ibound(ibi,3),ibound(ibi,1),ibound(ibi,2)
 116     format(1x,i5,2x,i6,5x,i2,5x,i2)
  46     ibound(ibi,3)=-99
      write(6,117)
      do 219 i=1,nyt
         ymom(i,2)=ymom(i,2)/ymom(i,1)
         ymom(i,3)=((ymom(i,3)/ymom(i,1))-(ymom(i,2)*ymom(i,2)))**.5d0
         ymom(i,4)=ymom(i,4)/ymom(i,1)
         if(model.eq.1)write(6,120)(ymom(1,j),j=1,6)
 120     format(1x,'descriptive statistics for dependent variable',/,
     &    1x,'# obs: ',f8.0,/,1x,'mean: ',f8.3,/,1x,'std dev: ',f8.3,
     &    /,1x,'within std dev: ',f8.3,/,1x,'minimum: ',f8.3,/,1x,
     &    'maximum: ',f8.3)
         if(model.eq.2)then
            ie=i*3
            ib=ie-2
            write(6,167)(ay2(j),j=ib,ie),(ymom(i,j),j=1,6)
 167        format(1x,'descriptive statistics for ',3a8,/,1x,'# obs: ',
     &       f8.0,/,1x,'mean: ',f8.3,/,1x,'std dev: ',f8.3,/,1x,
     &       'within std dev: ',f8.3,/,1x,'minimum: ',f8.3,/,1x,
     &       'maximum: ',f8.3)
            endif
 219     continue
      write(6,135)
 135  format(1x,'frequency for dependent variable')
      if(model.eq.1)write(6,137)
 137  format(1x,'iy',2x,'frequency')
      if(model.eq.2)write(6,168)(i,i=1,nyt)
 168  format(5x,'frequency',/,1x,'iy',2x,5(i3,5x))
      do 50 i=1,nyp1s
  50     write(6,136)i,(iydens(j,i),j=1,nyt)
 136  format(1x,i2,5(2x,i6))
      if(model.eq.2)then
         write(6,117)
         write(6,170)
 170     format(1x,'dependent variable patterns (by observation-year)')
         nyp1=ny+1
         ndeppat=0
         do 223 i=1,nobs
            if(incluv(i).eq.1)then
               istart=istartv(i)
               iend=iendv(i)
               do 224 j=istart,iend
                  if(iobsrvv(i,j).eq.1)then
                     if(ndeppat.gt.0)then
                        do 225 k=1,ndeppat
                           do 226 l=1,ny
                              iyvb(l)=0
                              if((iyv(i,l,j).ge.0).and.
     &                         (iyv(i,l,j).lt.10))iyvb(l)=1
                              if(iyvb(l).ne.ideppat(k,l))goto 225
 226                          continue
                           ideppat(k,nyp1)=ideppat(k,nyp1)+1
                           goto 224
 225                       continue
                        endif
                     ndeppat=ndeppat+1
                     do 227 k=1,ny
 227                    ideppat(ndeppat,k)=iyvb(k)
                     ideppat(ndeppat,nyp1)=1
                     endif
 224              continue
               endif
 223        continue
         write(6,171)(i,i=1,ny)
 171     format(20x,'dependent variable',/,1x,'rank',2x,'frequency',
     &    2x,6(2x,i1))
         do 221 i=1,ndeppat
            idepmax=-5
            idep=-5
            do 222 j=1,ndeppat
               if(ideppat(j,nyp1).gt.idepmax)then
                  idepmax=ideppat(j,nyp1)
                  idep=j
                  endif
 222           continue
            write(6,172)i,idepmax,(ideppat(idep,j),j=1,ny)
 172        format(1x,i4,2x,i6,5x,6(2x,i1))
            ideppat(idep,nyp1)=-10
 221        continue
         endif
      write(6,117)
      write(6,121)
 121  format(1x,'descriptive statistics for explanatory variables')
      write(6,122)
 122  format(49x,'within',/,1x,'variable',10x,'# obs',5x,'mean',6x,
     & 'std dev',3x,'std dev',3x,'minimum',3x,'maximum')
      do 48 i=1,nx
         ie=i*2
         ib=ie-1
         xmom(i,2)=xmom(i,2)/xmom(i,1)
         xmom(i,3)=((xmom(i,3)/xmom(i,1))-(xmom(i,2)*xmom(i,2)))**.5d0
         xmom(i,4)=xmom(i,4)/xmom(i,1)
         if(model.eq.1)write(6,123)(ax1(j),j=ib,ie),(xmom(i,j),j=1,6)
         if(model.eq.2)write(6,123)(ax2(j),j=ib,ie),(xmom(i,j),j=1,6)
 123     format(1x,2a8,2x,f8.0,5(2x,f8.3))
  48     continue
      write(6,117)
      write(6,138)
 138  format(1x,'missing variable analysis')
      write(6,139)
 139  format(1x,'variable',10x,'frequency')
      do 53 i=1,ncause
         do 73 j=1,ntime
            if(icause(j,i).gt.0)goto 74
  73        continue
         goto 53
  74     continue
         if(i.eq.1)write(6,142)
 142     format(1x,'dependent variable')
         if(i.gt.1)then
            ie=(i-1)*2
            ib=ie-1
            if(model.eq.1)write(6,143)(ax1(j),j=ib,ie)
            if(model.eq.2)write(6,143)(ax2(j),j=ib,ie)
 143        format(1x,2a8)
            endif    
         write(6,141)
 141     format(1x,'year',2x,'frequency',2x,'probability')
         do 51 j=1,ntime
            if(icause(j,i).gt.0)then
               prob=dble(float(icause(j,i)))/timet(j,1)
               write(6,140)j,icause(j,i),prob
 140           format(1x,i2,4x,i6,5x,f6.4)
               endif
  51        continue
  53     continue
      return
      end
c
c
      subroutine datget2
c
c       this subroutine reads data from generalized social survey.
c
      parameter(nbetas=20,nmisss=4,nobss=96529,npatts=1000,ntimes=34,
     & ntime12s=ntimes*12,nxs=nbetas,nys=6,ny3s=nys*3)
      implicit real*8(a-h,o-z)
      character*8 amiss(nmisss),axvar(nxs),ayvar(ny3s)
      real xextra,xv
      dimension adjtime(ntimes),ameanlag(nys,2,2),autocov(nys,3),
     & cmom1(nys,2),cmom2(nys,nys,2),corr(nys,nys),corrlag(nys,nys),
     & cov(nys,nys),covlag(nys,nys,2),dens(2,ntime12s),idens(2),
     & idepfreq1(nys,8),idepfreq2(nys,ntimes,8),idatev(nobss,ntimes),
     & ijunk(18),imiss(nmisss),ipatt(npatts,ntimes),iwave(ntimes),
     & jpatt(npatts),timet(2,2),varlag(nys,2),xmom(nxs,5),ymom(nys,8),
     & ymomt(nys,3)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      data adjtime/2,3,4,5,6,7,8,10,12,13,14,15,16,17,18,19,20,21,23,
     & 24,26,28,30,32,34,36,38,40,42,44,46,48,51,52/
      data amiss/'dep var','age','educ','female'/
      data axvar/'constant','age','educ','female','white','black',
     & 'income',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' ',' '/
      data ayvar/'party id',' ',' ','help sic','k',' ','help bla',
     & 'cks',' ','equal we','alth',' ','help poo','r',' ','help not',
     & ' ',' '/
      data time/2,3,4,5,6,7,8,10,12,13,14,15,16,17,18,19,20,21,23,24,
     & 26,28,30,32,34,36,38,40,42,44,46,48,51,52/
c
c       initialize.
c
      do 9 i=1,ny
         do 10 j=1,8
  10        idepfreq1(i,j)=0
   9     continue
      do 54 i=1,ny
         do 55 j=1,ntime
            do 56 k=1,8
  56           idepfreq2(i,j,k)=0
  55        continue
  54     continue
      do 12 i=1,nx
         do 13 j=1,3
  13        xmom(i,j)=0.d0
         xmom(i,4)=1.d20
  12     xmom(i,5)=-1.d20
      do 65 i=1,ny
         do 66 j=1,6
  66        ymom(i,j)=0.d0
         ymom(i,7)=1.d20
  65     ymom(i,8)=-1.d20
      do 76 i=1,ny
         do 77 j=1,i
            cmom2(i,j,1)=0.d0
  77        cmom2(i,j,2)=0.d0
  76     continue
      do 80 i=1,ny
         do 81 j=1,3
  81        autocov(i,j)=0.d0
  80     continue
      do 90 i=1,ny
         do 91 j=1,ny
            covlag(i,j,1)=0.d0
  91        covlag(i,j,2)=0.d0
  90     continue
      do 92 i=1,ny
         do 95 j=1,2
            ameanlag(i,j,1)=0.d0
  95        ameanlag(i,j,2)=0.d0
  92     continue
      do 97 i=1,ny
         varlag(i,1)=0.d0
  97     varlag(i,2)=0.d0
      do 67 i=1,ntime
  67     iwave(i)=0
      nmiss=4
      do 29 i=1,nmiss
  29     imiss(i)=0
      iobs=0
      open(unit=8,file='data/gss.prn')
      do 2 i=1,nobs
c
c       read variables:
c        1) idt               id #
c        2) iyear             year of interview
c        3) idate             month and date of interview
c        4) iballot           ballot used for interview
c        5) ipartyid          party id (dependent variable)
c        6) ihelpsick         government should pay for medical care
c                             (dependent variable)
c        7) ihlpblacks        government should aid blacks (dependent 
c                             variable)
c        8) ieqwealth         government should educe income 
c                             differences (dependent variable)
c        9) ihelppoor         government should improve standard of   
c                             living (dependent variable)
c       10) ihelpnot          help not (dependent variable)
c       10) iage              age (explanatory variable)
c       11) ieduc             education (explanatory variable)
c       12) isex              gender (explanatory variable)
c       13) irace             race (explanatory variable)
c       14) income            income brackets (explanatory variable) 
c
         read(8,*)idt,iyear,idate,intv1,intv2,iballot,iage,ieduc,isex,
     &    irace,income,ipartyid,ieqwealth,(ijunk(j),j=1,18),ihelppoor,
     &    ihelpnot,ihelpsick,ihlpblacks
c
c       create adjusted year.
c
         if(iyear.le.1990)then
            if(iyear.le.1982)then
               if(iyear.le.1976)then
                  if(iyear.eq.1972)iyearadj=1
                  if(iyear.eq.1973)iyearadj=2
                  if(iyear.eq.1974)iyearadj=3
                  if(iyear.eq.1975)iyearadj=4
                  if(iyear.eq.1976)iyearadj=5
                  goto 5
                  endif
               if(iyear.gt.1976)then
                  if(iyear.eq.1977)iyearadj=6
                  if(iyear.eq.1978)iyearadj=7
                  if(iyear.eq.1980)iyearadj=8
                  if(iyear.eq.1982)iyearadj=9
                  goto 5
                  endif
               endif
            if(iyear.gt.1982)then
               if(iyear.le.1986)then
                  if(iyear.eq.1983)iyearadj=10
                  if(iyear.eq.1984)iyearadj=11
                  if(iyear.eq.1985)iyearadj=12
                  if(iyear.eq.1986)iyearadj=13
                  goto 5
                  endif
               if(iyear.gt.1976)then
                  if(iyear.eq.1987)iyearadj=14
                  if(iyear.eq.1988)iyearadj=15
                  if(iyear.eq.1989)iyearadj=16
                  if(iyear.eq.1990)iyearadj=17
                  goto 5
                  endif
               endif
            endif
         if(iyear.gt.1990)then
            if(iyear.le.2006)then
               if(iyear.le.1998)then
                  if(iyear.eq.1991)iyearadj=18
                  if(iyear.eq.1993)iyearadj=19
                  if(iyear.eq.1994)iyearadj=20
                  if(iyear.eq.1996)iyearadj=21
                  if(iyear.eq.1998)iyearadj=22
                  goto 5
                  endif
               if(iyear.gt.1998)then
                  if(iyear.eq.2000)iyearadj=23
                  if(iyear.eq.2002)iyearadj=24
                  if(iyear.eq.2004)iyearadj=25
                  if(iyear.eq.2006)iyearadj=26
                  goto 5
                  endif
               endif
            if(iyear.gt.2006)then
               if(iyear.le.2014)then
                  if(iyear.eq.2008)iyearadj=27
                  if(iyear.eq.2010)iyearadj=28
                  if(iyear.eq.2012)iyearadj=29
                  if(iyear.eq.2014)iyearadj=30
                  goto 5
                  endif
               if(iyear.gt.2014)then
                  if(iyear.eq.2016)iyearadj=31
                  if(iyear.eq.2018)iyearadj=32
                  if(iyear.eq.2021)iyearadj=33
                  if(iyear.eq.2022)iyearadj=34
                  endif
               endif
            endif
   5     continue
c
c       check for existing person.
c
         if(iobs.gt.0)then
            do 3 j=1,iobs
               if(idt.eq.idv(j))then
                  jt=j
                  goto 4
                  endif
   3           continue
            endif
         iobs=iobs+1
         idv(iobs)=idt
         idatev(iobs,iyearadj)=idate
         jt=iobs
c
c       construct dependent variables.
c        1) party id
c        2) help sick
c        3) help blacks
c        4) equal wealth
c        5) help poor 
c        6) help not
c
   4     iyv(jt,1,iyearadj)=ipartyid+1
         iyv(jt,2,iyearadj)=ihelpsick
         iyv(jt,3,iyearadj)=ihlpblacks
         iyv(jt,4,iyearadj)=ieqwealth
         iyv(jt,5,iyearadj)=ihelppoor
         iyv(jt,6,iyearadj)=ihelpnot 
c
c       construct explanatory variables.
c        1) constant
c        2) age
c        3) education
c        4) female
c        5) white
c        6) black
c        7) income
c
         xv(jt,iyearadj,1)=1.
         xv(jt,iyearadj,2)=float(iage)
         xv(jt,iyearadj,3)=float(ieduc)
         xv(jt,iyearadj,4)=float(isex-1)
         xv(jt,iyearadj,5)=0.
         xv(jt,iyearadj,6)=0.
         if(irace.eq.1)xv(jt,iyearadj,5)=1.
         if(irace.eq.2)xv(jt,iyearadj,6)=1.
         if(income.le.4)then
            if(income.eq.1)xv(jt,iyearadj,7)=.5
            if(income.eq.2)xv(jt,iyearadj,7)=2.
            if(income.eq.3)xv(jt,iyearadj,7)=3.5
            if(income.eq.4)xv(jt,iyearadj,7)=4.5
            goto 2
            endif
         if((income.gt.4).and.(income.lt.8))then
            if(income.eq.5)xv(jt,iyearadj,7)=5.5
            if(income.eq.6)xv(jt,iyearadj,7)=6.5
            if(income.eq.7)xv(jt,iyearadj,7)=7.5
            if(income.eq.8)xv(jt,iyearadj,7)=9.
            goto 2 
            endif
         if(income.ge.9)then
            if(income.eq.9)xv(jt,iyearadj,7)=12.5
            if(income.eq.10)xv(jt,iyearadj,7)=17.5
            if(income.eq.11)xv(jt,iyearadj,7)=22.5
            if(income.eq.12)xv(jt,iyearadj,7)=40.
            endif
   2     continue
      close(8)
      nobs=iobs
      write(6,100)
 100  format(1x,50('='))
c
c       construct incluv for dependent variables.
c
      ncluv=nobs
      do 22 i=1,nobs
         do 23 j=1,ntime
            do 24 k=1,ny
               if((iyv(i,k,j).ge.1).and.(iyv(i,k,j).le.nseg(k)))then
                  incluv(i)=1
                  goto 22
                  endif
  24           continue
  23        continue
         incluv(i)=0
         ncluv=ncluv-1
         imiss(1)=imiss(1)+1
  22     continue
c
c       adjust incluv for age.
c
      do 32 i=1,nobs
         do 33 j=1,ntime
            if((xv(i,j,2).ge.0.).and.(xv(i,j,2).le.16))then
               do 34 k=1,ntime
                  if(xv(i,k,2).lt.0.)then
                     xv(i,k,2)=xv(i,j,2)+float(k-j)
                     if(xv(i,k,2).lt.0.)xv(i,k,2)=0.
                     endif
  34              continue
               goto 32
               endif
  33        continue
         if(incluv(i).eq.1)imiss(2)=imiss(2)+1
         incluv(i)=0
         ncluv=ncluv-1
         imiss(2)=imiss(2)+1
  32     continue
c
c       adjust incluv for education.
c
      do 25 i=1,nobs
         do 26 j=1,ntime
            if((xv(i,j,3).ge.0.).and.(xv(i,j,3).le.16))then
               do 27 k=1,ntime
                  if(xv(i,k,3).lt.0.)xv(i,k,3)=xv(i,j,3)
  27              continue
               goto 25
               endif
  26        continue
         if(incluv(i).eq.1)imiss(2)=imiss(2)+1
         incluv(i)=0
         ncluv=ncluv-1
         imiss(3)=imiss(3)+1
  25     continue
c
c       adjust incluv for female.
c
      do 35 i=1,nobs
         do 36 j=1,ntime
            if((xv(i,j,4).ge.0.).and.(xv(i,j,4).le.16))then
               do 37 k=1,ntime
                  if(xv(i,k,4).lt.0.)xv(i,k,4)=xv(i,j,4)
  37              continue
               goto 35
               endif
  36        continue
         if(incluv(i).eq.1)imiss(4)=imiss(4)+1
         incluv(i)=0
         ncluv=ncluv-1
         imiss(4)=imiss(4)+1
  35     continue
      nobs=ncluv
c
c       compute istartv, iendv, and iobsrvv.
c
      do 46 i=1,nobs
         iflag=0
         do 49 j=1,ntime
  49        iobsrvv(i,j)=0
         do 47 j=1,ntime
            if(incluv(i).eq.1)then
               if(iflag.eq.0)then
                  istartv(i)=j
                  iwave(j)=iwave(j)+1
                  iendv(i)=j
                  iobsrvv(i,j)=1
                  iflag=1
                  goto 47   
                  endif
               if(iflag.eq.1)then
                  iwave(j)=iwave(j)+1
                  iendv(i)=j
                  iobsrvv(i,j)=1
                  endif
               endif
  47        continue
         if(istartv(i).lt.13)then
            if(iendv(i).ge.13)istartv(i)=13
            if(iendv(i).lt.13)incluv(i)=0
            endif
  46     continue
c
c       compute descriptive statistics for dependent variables.
c
      do 68 i=1,nobs
         if(incluv(i).eq.0)goto 68
         do 72 j=1,ny
            do 73 k=1,3
  73           ymomt(j,k)=0.d0
  72        continue
         istart=istartv(i)
         iend=iendv(i)
         do 69 j=istart,iend
            do 70 k=1,ny
               if((iyv(i,k,j).ge.1).and.((k.ne.1).or.
     &          (iyv(i,k,j).ne.8)))then
                  yt=dble(float(iyv(i,k,j)))
                  ymom(k,1)=ymom(k,1)+1.d0
                  ymom(k,2)=ymom(k,2)+yt
                  ymom(k,3)=ymom(k,3)+(yt*yt)
                  if(ymom(k,7).gt.yt)ymom(k,7)=yt
                  if(ymom(k,8).lt.yt)ymom(k,8)=yt
                  ymomt(k,1)=ymomt(k,1)+1.d0
                  ymomt(k,2)=ymomt(k,2)+yt
                  ymomt(k,3)=ymomt(k,3)+(yt*yt)
                  endif
  70           continue
            do 71 k=1,ny
               if(ymomt(k,1).gt.0.d0)then
                  ymom(k,4)=ymom(k,4)+ymomt(k,1)
                  ymom(k,5)=ymom(k,5)+(ymomt(k,2)/ymomt(k,1))
                  endif
  71           continue
            do 74 k=1,ny
               if((iyv(i,k,j).ge.1).and.((k.ne.1).or.
     &          (iyv(i,k,j).ne.8)))then
                  ytk=dble(float(iyv(i,k,j)))
                  cmom1(k,1)=cmom1(k,1)+1.d0
                  cmom1(k,2)=cmom1(k,2)+ytk
                  do 75 l=1,k
                     if(iyv(i,l,j).ge.0)then
                        ytl=dble(float(iyv(i,l,j)))
                        cmom2(k,l,1)=cmom2(k,l,1)+1.d0
                        cmom2(k,l,2)=cmom2(k,l,2)+(ytk*ytl)
                        endif
  75                 continue
                  endif
  74           continue
  69        continue
         istartt=istartv(i)
         if(istartt.lt.2)istartt=2
         do 79 j=istartt,iend
            jm1=j-1
            do 82 k=1,ny
               if((iyv(i,k,j).ge.1).and.(iyv(i,k,jm1).ge.1).and.
     &          ((k.ne.1).or.(iyv(i,k,j).ne.8)))then
                  autocov(k,1)=autocov(k,1)+1.d0
                  autocov(k,2)=autocov(k,2)+dble(float(iyv(i,k,jm1)))
                  autocov(k,3)=autocov(k,3)+dble(float(iyv(i,k,j)*
     &             iyv(i,k,jm1)))
                  endif
  82           continue
  79        continue
         do 87 j=istartt,iend
            jm1=j-1
            do 99 k=1,ny
               if((iyv(i,k,j).ge.1).and.((k.ne.1).or.
     &          (iyv(i,k,j).ne.8)))then
                  yvtk=dble(float(iyv(i,k,j)))
                  ameanlag(k,1,1)=ameanlag(k,1,1)+1.d0
                  ameanlag(k,1,2)=ameanlag(k,1,2)+yvtk
                  varlag(k,1)=varlag(k,1)+(yvtk*yvtk)
                  endif
               if((iyv(i,k,jm1).ge.1).and.((k.ne.1).or.
     &          (iyv(i,k,jm1).ne.8)))then
                  ameanlag(k,2,1)=ameanlag(k,2,1)+1.d0
                  yvtkm=dble(float(iyv(i,k,jm1)))
                  ameanlag(k,2,2)=ameanlag(k,2,2)+yvtkm
                  varlag(k,2)=varlag(k,2)+(yvtkm*yvtkm)
                  endif
  99           continue
            do 88 k=1,ny
               if((iyv(i,k,j).ge.1).and.((k.ne.1).or.
     &          (iyv(i,k,j).ne.8)))then
                  yvtk=dble(float(iyv(i,k,j)))
                  do 89 l=1,ny
                     if((iyv(i,l,jm1).ge.1).and.((l.ne.1).or.
     &                (iyv(i,l,jm1).ne.8)))then
                        yvtlm=dble(float(iyv(i,l,jm1)))
                        covlag(k,l,1)=covlag(k,l,1)+1.d0
                        covlag(k,l,2)=covlag(k,l,2)+(yvtk*yvtlm)
                        endif
  89                 continue
                  endif
  88           continue
  87        continue
  68     continue
c
c       compute descriptive statistics for explanatory variables.
c
      do 30 i=1,nobs
         if(incluv(i).eq.1)then 
            do 31 j=1,ntime
               if(xv(i,j,1).eq.1.)then
                  do 11 k=1,7
                     xvt=dble(xv(i,j,k))
                     xmom(k,1)=xmom(k,1)+1.d0
                     xmom(k,2)=xmom(k,2)+xvt
                     xmom(k,3)=xmom(k,3)+(xvt*xvt)
                     if(xmom(k,4).gt.xvt)xmom(k,4)=xvt
                     if(xmom(k,5).lt.xvt)xmom(k,5)=xvt
  11                 continue
                  endif
  31           continue
            endif
  30     continue
c
c       output missing value analysis.
c
      write(6,110)
 110  format(1x,'missing value analysis')
      write(6,112)
 112  format(1x,'#',3x,'variable',2x,'number missing')
      do 28 i=1,nmiss
  28     write(6,111)i,amiss(i),imiss(i)
 111  format(1x,i2,2x,a8,2x,i8)
      write(6,113)nobs
 113  format(1x,'nobs = ',i5)
c
c       output frequency of dependent variables.
c
      write(6,103)
 103  format(1x,'frequency of dependent variables')
      write(6,102)(i,i=1,8)
 102  format(1x,'variable',16x,8(2x,i1,4x))
      do 50 i=1,nobs
         if(incluv(i).eq.1)then
            istart=istartv(i)
            iend=iendv(i)
            do 8 j=1,ny
               do 51 k=istart,iend
                  iyvt=iyv(i,j,k)
                  if((iyvt.ge.1).and.(iyvt.le.nseg(j)))then
                     idepfreq1(j,iyvt)=idepfreq1(j,iyvt)+1
                     idepfreq2(j,k,iyvt)=idepfreq2(j,k,iyvt)+1
                     endif
  51              continue
   8           continue
            endif
  50     continue
      do 21 i=1,ny
         ie=i*3
         ib=ie-2
  21     write(6,101)(ayvar(j),j=ib,ie),(idepfreq1(i,j),j=1,nseg(i))
 101  format(1x,3a8,8(2x,i5))
      write(6,117)
 117  format(1x,'frequency of dependent variable disaggregated by',
     & ' year')
      do 52 i=1,ny
         ie=i*3
         ib=ie-2
         write(6,118)(ayvar(j),j=ib,ie)
 118     format(1x,'statistics for ',3a8)
         write(6,119)(j,j=1,8)
 119     format(1x,'year',8(2x,i2,2x))
         do 53 j=1,ntime
  53        write(6,120)j,(idepfreq2(i,j,k),k=1,nseg(i))
 120     format(1x,i2,2x,8(2x,i4))
  52     continue
c
c       output descriptive statistics of explanatory variables.
c
      write(6,100)
      write(6,104)
 104  format(1x,'descriptive statistics of explanatory variables')
      write(6,105)
 105  format(1x,'variable',2x,'# obs',5x,'mean',6x,'std dev',3x,
     & 'minimum',3x,'maximum')
      do 6 i=1,nx
         xmom(i,2)=xmom(i,2)/xmom(i,1)
         xmom(i,3)=((xmom(i,3)/xmom(i,1))-(xmom(i,2)*xmom(i,2)))**.5d0
   6     write(6,106)axvar(i),(xmom(i,j),j=1,5)
 106  format(1x,a8,2x,f8.0,4(2x,f6.3,2x))
c
c       output frequency of patterns of data by waves.
c
      write(6,100)
      write(6,107)
 107  format(1x,'frequency of patterns of data by waves')
      npatt=0
      do 14 i=1,nobs
         if(npatt.gt.0)then
            do 15 j=1,npatt
               do 16 k=1,ntime
                  if((ipatt(j,k).eq.1).and.((iyv(i,1,k).lt.1).or.
     &             (iyv(i,1,k).gt.nseg(1))))goto 15
                  if((ipatt(j,k).eq.0).and.(iyv(i,1,k).ge.1).and.
     &             (iyv(i,1,k).le.nseg(1)))goto 15
  16              continue
               jpatt(j)=jpatt(j)+1
               goto 14
  15           continue
            endif
         npatt=npatt+1
         jpatt(npatt)=1
         do 19 j=1,ntime
            do 17 k=1,ny
               if((iyv(i,k,j).ge.1).and.(iyv(i,k,j).le.nseg(k)))then
                  ipatt(npatt,j)=1        
                  goto 19
                  endif
  17           continue
            ipatt(npatt,j)=0
  19        continue
  14     continue
      write(6,108)
 108  format(1x,'rank',3x,'# obs',2x,'pattern')
      do 7 i=1,npatt
         jpatmax=0
         do 40 j=1,npatt
            if(jpatt(j).gt.jpatmax)then
               jpatmax=jpatt(j)
               jpati=j
               endif
  40        continue
         if(jpatmax.ge.5)then
            write(6,109)i,jpatt(jpati),(ipatt(jpati,j),j=1,ntime)
 109        format(1x,i5,2x,i5,34(2x,i2))
            jpatt(jpati)=0
            endif
   7     continue
c
c       compute the density of times between interviews with and 
c       without using month and date.
c
      do 48 i=1,2
         do 41 j=1,ntime
  41        dens(i,j)=0.d0
  48     continue
      do 38 i=1,nobs
         if(incluv(i).eq.1)then
            iflag=0
            do 39 j=1,ntime
               if(xv(i,j,1).eq.1.)then
                  if(iflag.eq.0)then
                     timet(1,1)=adjtime(j)
                     month=idatev(i,j)/100
                     amonth=dble(float(month))/12.d0
                     timet(2,1)=adjtime(j)+amonth
                     iflag=1
                     goto 39
                     endif
                  if(iflag.eq.1)then
                     timet(1,2)=adjtime(j)
                     month=idatev(i,j)/100
                     amonth=dble(float(month))/12.d0
                     timet(2,2)=adjtime(j)+amonth
                     dtime=timet(1,2)-timet(1,1)
                     itime=ifix(sngl(dtime*12.d0))
                     dens(1,itime)=dens(1,itime)+1.d0
                     dtime=timet(2,2)-timet(2,1)
                     itime=ifix(sngl(dtime*12.d0))
                     dens(2,itime)=dens(2,itime)+1.d0
                     timet(1,1)=timet(1,2)
                     timet(2,1)=timet(2,2)
                     endif
                  endif
  39           continue
            endif
  38     continue
      write(6,100)
      write(6,114)
 114  format(1x,'density of times between interviews with and'
     & ' without month and date information')
      write(6,115)
 115  format(15x,'density',2x,'density',/,15x,'without',2x,'with',/,
     & 1x,'years',2x,'month',2x,'months',3x,'months')
      ndens=ntime*12
      do 43 i=1,2
         denstot=0.d0
         do 44 j=1,ndens
  44        denstot=denstot+dens(i,j)
         do 45 j=1,ndens
  45        dens(i,j)=dens(i,j)/denstot
  43     continue
      do 42 i=1,ndens
         idens(1)=i/12
         idens(2)=i-(idens(1)*12)
  42     write(6,116)(idens(j),j=1,2),(dens(j,i),j=1,2)
 116  format(1x,i3,4x,i2,5x,f5.3,2x,f5.3)
c
c       output number of observations for each wave.
c
      write(6,100)
      write(6,129)
 129  format(1x,'number of observations for each wave')
      write(6,132)
 132  format(1x,'wave',2x,'# obs')
      do 57 i=1,ntime
  57     write(6,133)i,iwave(i)
 133  format(1x,i2,4x,i5)
c
c       output descriptive statistics for each dependent variable
c       including person-specific standard deviation and idiosyncratic
c       standard deviation.
c
      write(6,100)
      write(6,130)
 130  format(1x,'descriptive statistics for dependent variables')
      write(6,128)
 128  format(67x,'person-',3x,'idiosyn-',/,67x,'specific',2x,'cratic',
     & /,1x,'variable',18x,'# obs',5x,'mean',6x,'std dev',13x,
     & 'std dev',3x,'std dev',3x,'minimum',3x,'maximum')
      do 58 i=1,ny
         ie=i*3
         ib=ie-2
         ymom(i,2)=ymom(i,2)/ymom(i,1)
         ymom(i,3)=((ymom(i,3)/ymom(i,1))-(ymom(i,2)*ymom(i,2)))**.5d0
         ymom(i,5)=(ymom(i,5)/ymom(i,4))**.5d0
         ymom(i,6)=((ymom(i,3)*ymom(i,3))-(ymom(i,5)*ymom(i,5)))**.5d0
  58     write(6,121)(ayvar(j),j=ib,ie),(ymom(i,j),j=1,8)
 121  format(1x,3a8,2x,f8.0,2(2x,f8.3),2x,f8.0,2(2x,f8.3),2(2x,f8.1))
c
c       output correlation matrix of contemporaneous dependent 
c       variables.
c
      write(6,100)
      write(6,124)(i,i=1,ny)
 124  format(2x,6(2x,i1,5x))
      do 83 i=1,ny
  83     cmom1(i,2)=cmom1(i,2)/cmom1(i,1)
      do 59 i=1,ny
         do 60 j=1,i
            cov(i,j)=(cmom2(i,j,2)/cmom2(i,j,1))-(cmom1(i,2)*
     &       cmom1(j,2))
  60        cov(j,i)=cov(i,j)
  59     continue
         do 61 i=1,ny
            do 62 j=1,i
  62           corr(i,j)=cov(i,j)/((cov(i,i)*cov(j,j))**.5d0)
  61     write(6,122)i,(corr(i,j),j=1,i)
 122  format(1x,i1,6(2x,f6.3))
c
c       output first autocovariance for each dependent variable.
c
      write(6,100)
      write(6,125)
 125  format(1x,'first autocovariances for each dependent variable')
      write(6,127)
 127  format(1x,'variable',18x,'autocovariance')
      do 63 i=1,ny
         ie=i*3
         ib=ie-2
         autocov(i,2)=autocov(i,2)/autocov(i,1)
         autocov(i,3)=(autocov(i,3)/autocov(i,1))-(autocov(i,2)*
     &    autocov(i,2))
  63     write(6,126)(ayvar(j),j=ib,ie),autocov(i,3)
 126  format(1x,3a8,2x,f8.3)
c
c       output lagged covariance table.
c
      write(6,100)
      write(6,135)
 135  format(1x,'lagged covariance matrix')
      write(6,136)(i,i=1,ny)
 136  format(11x,'lagged',/,11x,6(i2,8x))
      do 84 i=1,ny
         ameanlag(i,1,2)=ameanlag(i,1,2)/ameanlag(i,1,1)
  84     ameanlag(i,2,2)=ameanlag(i,2,2)/ameanlag(i,2,1)
      do 85 i=1,ny
         do 86 j=1,ny
  86        covlag(i,j,2)=(covlag(i,j,2)/covlag(i,j,1))-
     &       (ameanlag(i,1,2)*ameanlag(j,2,2))
  85     write(6,134)i,(covlag(i,j,2),j=1,ny)
 134  format(1x,i3,7x,6(f8.3,2x))
c
c       output lagged correlation table.
c
      write(6,100)
      write(6,137)
 137  format(1x,'lagged correlation matrix')
      write(6,136)(i,i=1,ny)
      do 94 i=1,ny
         do 98 j=1,2
  98        varlag(i,j)=(varlag(i,j)/ameanlag(i,j,1))-(ameanlag(i,j,2)
     &       *ameanlag(i,j,2))
  94     continue
      do 200 i=1,ny
         do 93 j=1,ny
  93        corrlag(i,j)=covlag(i,j,2)/((varlag(i,1)*varlag(j,2))
     &       **.5d0)
 200     write(6,134)i,(corrlag(i,j),j=1,ny)
      return
      end
c
c
      subroutine datget3
c
c       this subroutine reads data from generalized social survey,
c       second try.
c
      parameter(nbetas=20,nmisss=4,nobss=96529,npatts=1000,ntimes=34,
     & ntime12s=ntimes*12,nxs=nbetas,nys=6,ny3s=nys*3)
      implicit real*8(a-h,o-z)
      character*8 amiss(nmisss),axvar(nxs),ayvar(ny3s)
      real xextra,xv
      dimension adjtime(ntimes),ameanlag(nys,2,2),autocov(nys,3),
     & cmom1(nys,2),cmom2(nys,nys,2),corr(nys,nys),corrlag(nys,nys),
     & cov(nys,nys),covlag(nys,nys,2),dens(2,ntime12s),
     & idatev(nobss,ntimes),idens(2),idepfreq1(nys,8),
     & idepfreq2(nys,ntimes,8),imiss(nmisss),ipatt(npatts,ntimes),
     & iwave(ntimes),jpatt(npatts),timet(2,2),varlag(nys,2),
     & xmom(nxs,5),ymom(nys,8),ymomt(nys,3)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      data adjtime/2,3,4,5,6,7,8,10,12,13,14,15,16,17,18,19,20,21,23,
     & 24,26,28,30,32,34,36,38,40,42,44,46,48,51,52/
      data amiss/'dep var','age','educ','female'/
      data axvar/'constant','age','hs dipl','college','adv deg',
     & 'female','white','black','inc mid','inc high',' ',' ',' ',' ',
     & ' ',' ',' ',' ',' ',' '/
      data ayvar/'polit id','eology',' ','military',' spendin','g',
     & 'welfare',' ',' ','immigran','ts',' ','party id',' ',' ',
     & 'environm','ent',' '/
      data time/2,3,4,5,6,7,8,10,12,13,14,15,16,17,18,19,20,21,23,24,
     & 26,28,30,32,34,36,38,40,42,44,46,48,51,52/
c
c       initialize.
c
      do 9 i=1,ny
         do 10 j=1,8
  10        idepfreq1(i,j)=0
   9     continue
      do 54 i=1,ny
         do 55 j=1,ntime
            do 56 k=1,8
  56           idepfreq2(i,j,k)=0
  55        continue
  54     continue
      do 12 i=1,nx
         do 13 j=1,3
  13        xmom(i,j)=0.d0
         xmom(i,4)=1.d20
  12     xmom(i,5)=-1.d20
      do 65 i=1,ny
         do 66 j=1,6
  66        ymom(i,j)=0.d0
         ymom(i,7)=1.d20
  65     ymom(i,8)=-1.d20
      do 201 i=1,ny
         cmom1(i,1)=0.d0
 201     cmom1(i,2)=0.d0
      do 76 i=1,ny
         do 77 j=1,i
            cmom2(i,j,1)=0.d0
  77        cmom2(i,j,2)=0.d0
  76     continue
      do 80 i=1,ny
         do 81 j=1,3
  81        autocov(i,j)=0.d0
  80     continue
      do 90 i=1,ny
         do 91 j=1,ny
            covlag(i,j,1)=0.d0
  91        covlag(i,j,2)=0.d0
  90     continue
      do 92 i=1,ny
         do 95 j=1,2
            ameanlag(i,j,1)=0.d0
  95        ameanlag(i,j,2)=0.d0
  92     continue
      do 97 i=1,ny
         varlag(i,1)=0.d0
  97     varlag(i,2)=0.d0
      do 67 i=1,ntime
  67     iwave(i)=0
      nmiss=4
      do 29 i=1,nmiss
  29     imiss(i)=0
      iobs=0
      open(unit=8,file='data/gss2.prn')
      do 2 i=1,nobs
c
c       read variables:
c        1) idt               id #
c        2) iyear             year of interview
c        3) iage              age (explanatory variable)
c        4) ieduc             education (explanatory variable)
c        5) ipolyv            political ideology (dependent variable)
c        6) irace             race (explanatory variable)
c        7) isex              gender (explanatory variable)
c        8) ipanel            panel #
c        9) idate             date interviewed
c       10) imonth            month interviewed
c       11) itrvdate          days from 1/1/60 to date interviewed  
c       12) imilitary         government should spend on military
c                             (dependent variable)
c       13) ijunk1            government should spend on law 
c                             enforcement (extra dependent variable)
c       14) ienvir            government should spend more on         
c                             environment (dependent variable)
c       15) iwelfare          government should spend more on welfare
c                             (dependent variable)
c       16) ijunk3            government should spend more on improving
c                             the condition of blacks (extra dependent 
c                             variable)
c       17) ijunk4            government should spend more on 
c                             supporting scientific research (extra 
c                             dependent variable)
c       18) immigr            do you think the number of immigrants 
c                             being let in is ... (dependent variable)
c       19) ipartyid          party id (dependent variable)
c       20) income            family income (explanatory variable)
c
         read(8,*)idt,iyear,iage,ieduc,ipolyv,irace,isex,ipanel,idate,
     &    imonth,itrvdate,imilitary,ijunk1,ienvir,iwelfare,ijunk3,
     &    ijunk4,immigr,ipartyid,income
c
c       create adjusted year.
c
         iyearadj=((iyear-2006)/2)+1
c
c       check for existing person.
c
         if(iobs.gt.0)then
            do 3 j=1,iobs
               if(idt.eq.idv(j))then
                  jt=j
                  goto 4
                  endif
   3           continue
            endif
         iobs=iobs+1
         idv(iobs)=idt
         idatev(iobs,iyearadj)=((imonth-1)*30)+idate
         jt=iobs
c
c       construct dependent variables.
c        1) political ideology (ipolyv)
c        2) military spending (imilitary)
c        3) welfare (iwelfare)
c        4) immigrant (immigr)
c        5) party id (ipartyid)
c        6) environment (ienvir)
c
   4     iyv(jt,1,iyearadj)=ipolyv
         iyv(jt,2,iyearadj)=imilitary
         iyv(jt,3,iyearadj)=iwelfare
         iyv(jt,4,iyearadj)=immigr
         iyv(jt,5,iyearadj)=ipartyid+1
         iyv(jt,6,iyearadj)=ienvir   
c
c       construct explanatory variables.
c        1) constant
c        2) age
c        3) high school diploma
c        4) college degree
c        5) advanced degree
c        6) female
c        7) white
c        8) black
c        9) income 15k - 25k
c       10) income 25k or more
c
         xv(jt,iyearadj,1)=1.
         xv(jt,iyearadj,2)=float(iage)
         do 202 j=3,5
 202        xv(jt,iyearadj,j)=0.
         if((ieduc.eq.1).or.(ieduc.eq.2))xv(jt,iyearadj,3)=1.
         if(ieduc.eq.3)xv(jt,iyearadj,4)=1.
         if(ieduc.eq.4)xv(jt,iyearadj,5)=1.
         xv(jt,iyearadj,6)=float(isex-1)
         xv(jt,iyearadj,7)=0.
         xv(jt,iyearadj,8)=0.
         if(irace.eq.1)xv(jt,iyearadj,7)=1.
         if(irace.eq.2)xv(jt,iyearadj,8)=1.
         xv(jt,iyearadj,9)=0.
         xv(jt,iyearadj,10)=0.
         if(income.eq.2)xv(jt,iyearadj,9)=1.
         if(income.eq.3)xv(jt,iyearadj,10)=1.
   2     continue
      close(8)
      nobs=iobs
      write(6,100)
 100  format(1x,50('='))
c
c       construct incluv for dependent variables.
c
      do 22 i=1,nobs
         do 23 j=1,ntime
            do 24 k=1,ny
               if((iyv(i,k,j).ge.1).and.(iyv(i,k,j).le.nseg(k)))then
                  incluv(i)=1
                  goto 22
                  endif
  24           continue
  23        continue
         incluv(i)=0
         imiss(1)=imiss(1)+1
  22     continue
c
c       adjust incluv for age.
c
      do 32 i=1,nobs
         do 33 j=1,ntime
            if(xv(i,j,2).ge.16.)then
               do 34 k=1,ntime
                  if(xv(i,k,2).lt.16.)then
                     xv(i,k,2)=xv(i,j,2)+float(2*(k-j))
                     endif
  34              continue
               goto 32
               endif
  33        continue
         if(incluv(i).eq.1)imiss(2)=imiss(2)+1
         incluv(i)=0
         imiss(2)=imiss(2)+1
  32     continue
c
c       adjust incluv for female.
c
      do 35 i=1,nobs
         do 36 j=1,ntime
            if((xv(i,j,6).ge.0.).and.(xv(i,j,6).le.1.))then
               do 37 k=1,ntime
                  if(xv(i,k,6).lt.0.)xv(i,k,6)=xv(i,j,6)
  37              continue
               goto 35
               endif
  36        continue
         if(incluv(i).eq.1)imiss(4)=imiss(4)+1
         incluv(i)=0
         imiss(4)=imiss(4)+1
  35     continue
c
c       compute istartv, iendv, and iobsrvv.
c
      do 46 i=1,nobs
         iflag=0
         do 49 j=1,ntime
  49        iobsrvv(i,j)=0
         do 47 j=1,ntime
            if(incluv(i).eq.1)then
               if(iflag.eq.0)then
                  do 204 k=1,ny
                     if(iyv(i,k,j).gt.0)then
                        istartv(i)=j
                        iwave(j)=iwave(j)+1
                        iendv(i)=j
                        iobsrvv(i,j)=1
                        iflag=1
                        goto 47   
                        endif
 204                 continue
                  endif
               if(iflag.eq.1)then
                  do 203 k=1,ny
                     if(iyv(i,k,j).gt.0)then
                        iwave(j)=iwave(j)+1
                        iendv(i)=j
                        iobsrvv(i,j)=1
                        endif
 203                 continue
                  endif
               endif
  47        continue
  46     continue
c
c       compute descriptive statistics for dependent variables.
c
      do 68 i=1,nobs
         if(incluv(i).eq.0)goto 68
         do 72 j=1,ny
            do 73 k=1,3
  73           ymomt(j,k)=0.d0
  72        continue
         istart=istartv(i)
         iend=iendv(i)
         do 69 j=istart,iend
            do 70 k=1,ny
               if(iyv(i,k,j).ge.1)then
                  yt=dble(float(iyv(i,k,j)))
                  ymom(k,1)=ymom(k,1)+1.d0
                  ymom(k,2)=ymom(k,2)+yt
                  ymom(k,3)=ymom(k,3)+(yt*yt)
                  if(ymom(k,7).gt.yt)ymom(k,7)=yt
                  if(ymom(k,8).lt.yt)ymom(k,8)=yt
                  ymomt(k,1)=ymomt(k,1)+1.d0
                  ymomt(k,2)=ymomt(k,2)+yt
                  ymomt(k,3)=ymomt(k,3)+(yt*yt)
                  endif
  70           continue
            do 71 k=1,ny
               if(ymomt(k,1).gt.0.d0)then
                  ymom(k,4)=ymom(k,4)+ymomt(k,1)
                  ymom(k,5)=ymom(k,5)+(ymomt(k,2)/ymomt(k,1))
                  endif
  71           continue
            do 74 k=1,ny
               if(iyv(i,k,j).ge.1)then 
                  ytk=dble(float(iyv(i,k,j)))
                  cmom1(k,1)=cmom1(k,1)+1.d0
                  cmom1(k,2)=cmom1(k,2)+ytk
                  do 75 l=1,k
                     if(iyv(i,l,j).ge.0)then
                        ytl=dble(float(iyv(i,l,j)))
                        cmom2(k,l,1)=cmom2(k,l,1)+1.d0
                        cmom2(k,l,2)=cmom2(k,l,2)+(ytk*ytl)
                        endif
  75                 continue
                  endif
  74           continue
  69        continue
         istartt=istartv(i)
         if(istartt.lt.2)istartt=2
         do 79 j=istartt,iend
            jm1=j-1
            do 82 k=1,ny
               if((iyv(i,k,j).ge.1).and.(iyv(i,k,jm1).ge.1))then
                  autocov(k,1)=autocov(k,1)+1.d0
                  autocov(k,2)=autocov(k,2)+dble(float(iyv(i,k,jm1)))
                  autocov(k,3)=autocov(k,3)+dble(float(iyv(i,k,j)*
     &             iyv(i,k,jm1)))
                  endif
  82           continue
  79        continue
         do 87 j=istartt,iend
            jm1=j-1
            do 99 k=1,ny
               if((iyv(i,k,j).ge.1).and.((k.ne.1).or.
     &          (iyv(i,k,j).ne.8)))then
                  yvtk=dble(float(iyv(i,k,j)))
                  ameanlag(k,1,1)=ameanlag(k,1,1)+1.d0
                  ameanlag(k,1,2)=ameanlag(k,1,2)+yvtk
                  varlag(k,1)=varlag(k,1)+(yvtk*yvtk)
                  endif
               if(iyv(i,k,jm1).ge.1)then
                  ameanlag(k,2,1)=ameanlag(k,2,1)+1.d0
                  yvtkm=dble(float(iyv(i,k,jm1)))
                  ameanlag(k,2,2)=ameanlag(k,2,2)+yvtkm
                  varlag(k,2)=varlag(k,2)+(yvtkm*yvtkm)
                  endif
  99           continue
            do 88 k=1,ny
               if((iyv(i,k,j).ge.1).and.((k.ne.1).or.
     &          (iyv(i,k,j).ne.8)))then
                  yvtk=dble(float(iyv(i,k,j)))
                  do 89 l=1,ny
                     if((iyv(i,l,jm1).ge.1).and.((l.ne.1).or.
     &                (iyv(i,l,jm1).ne.8)))then
                        yvtlm=dble(float(iyv(i,l,jm1)))
                        covlag(k,l,1)=covlag(k,l,1)+1.d0
                        covlag(k,l,2)=covlag(k,l,2)+(yvtk*yvtlm)
                        endif
  89                 continue
                  endif
  88           continue
  87        continue
  68     continue
c
c       compute descriptive statistics for explanatory variables.
c
      do 30 i=1,nobs
         if(incluv(i).eq.1)then 
            do 31 j=1,ntime
               if(xv(i,j,1).eq.1.)then
                  do 11 k=1,nx
                     xvt=dble(xv(i,j,k))
                     xmom(k,1)=xmom(k,1)+1.d0
                     xmom(k,2)=xmom(k,2)+xvt
                     xmom(k,3)=xmom(k,3)+(xvt*xvt)
                     if(xmom(k,4).gt.xvt)xmom(k,4)=xvt
                     if(xmom(k,5).lt.xvt)xmom(k,5)=xvt
  11                 continue
                  endif
  31           continue
            endif
  30     continue
c
c       output missing value analysis.
c
      write(6,110)
 110  format(1x,'missing value analysis')
      write(6,112)
 112  format(1x,'#',3x,'variable',2x,'number missing')
      do 28 i=1,nmiss
  28     write(6,111)i,amiss(i),imiss(i)
 111  format(1x,i2,2x,a8,2x,i8)
      write(6,113)nobs
 113  format(1x,'nobs = ',i5)
c
c       output frequency of dependent variables.
c
      write(6,103)
 103  format(1x,'frequency of dependent variables')
      write(6,102)(i,i=1,7)
 102  format(1x,'variable',16x,8(2x,i1,4x))
      do 50 i=1,nobs
         if(incluv(i).eq.1)then
            istart=istartv(i)
            iend=iendv(i)
            do 8 j=1,ny
               do 51 k=istart,iend
                  iyvt=iyv(i,j,k)
                  if((iyvt.ge.1).and.(iyvt.le.nseg(j)))then
                     idepfreq1(j,iyvt)=idepfreq1(j,iyvt)+1
                     idepfreq2(j,k,iyvt)=idepfreq2(j,k,iyvt)+1
                     endif
  51              continue
   8           continue
            endif
  50     continue
      do 21 i=1,ny
         ie=i*3
         ib=ie-2
  21     write(6,101)(ayvar(j),j=ib,ie),(idepfreq1(i,j),j=1,nseg(i))
 101  format(1x,3a8,8(2x,i5))
      write(6,117)
 117  format(1x,'frequency of dependent variable disaggregated by',
     & ' year')
      do 52 i=1,ny
         ie=i*3
         ib=ie-2
         write(6,118)(ayvar(j),j=ib,ie)
 118     format(1x,'statistics for ',3a8)
         write(6,119)(j,j=1,8)
 119     format(1x,'year',8(2x,i2,2x))
         do 53 j=1,ntime
  53        write(6,120)j,(idepfreq2(i,j,k),k=1,nseg(i))
 120     format(1x,i2,2x,8(2x,i4))
  52     continue
c
c       output descriptive statistics of explanatory variables.
c
      write(6,100)
      write(6,104)
 104  format(1x,'descriptive statistics of explanatory variables')
      write(6,105)
 105  format(1x,'variable',2x,'# obs',5x,'mean',6x,'std dev',3x,
     & 'minimum',3x,'maximum')
      do 6 i=1,nx
         xmom(i,2)=xmom(i,2)/xmom(i,1)
         xmom(i,3)=((xmom(i,3)/xmom(i,1))-(xmom(i,2)*xmom(i,2)))**.5d0
   6     write(6,106)axvar(i),(xmom(i,j),j=1,5)
 106  format(1x,a8,2x,f8.0,4(2x,f6.3,2x))
c
c       output frequency of patterns of data by waves.
c
      write(6,100)
      write(6,107)
 107  format(1x,'frequency of patterns of data by waves')
      npatt=0
      do 14 i=1,nobs
         if(npatt.gt.0)then
            do 15 j=1,npatt
               do 16 k=1,ntime
                  if((ipatt(j,k).eq.1).and.((iyv(i,1,k).lt.1).or.
     &             (iyv(i,1,k).gt.nseg(1))))goto 15
                  if((ipatt(j,k).eq.0).and.(iyv(i,1,k).ge.1).and.
     &             (iyv(i,1,k).le.nseg(1)))goto 15
  16              continue
               jpatt(j)=jpatt(j)+1
               goto 14
  15           continue
            endif
         npatt=npatt+1
         jpatt(npatt)=1
         do 19 j=1,ntime
            do 17 k=1,ny
               if((iyv(i,k,j).ge.1).and.(iyv(i,k,j).le.nseg(k)))then
                  ipatt(npatt,j)=1        
                  goto 19
                  endif
  17           continue
            ipatt(npatt,j)=0
  19        continue
  14     continue
      write(6,108)
 108  format(1x,'rank',3x,'# obs',2x,'pattern')
      do 7 i=1,npatt
         jpatmax=0
         do 40 j=1,npatt
            if(jpatt(j).gt.jpatmax)then
               jpatmax=jpatt(j)
               jpati=j
               endif
  40        continue
         if(jpatmax.ge.5)then
            write(6,109)i,jpatt(jpati),(ipatt(jpati,j),j=1,ntime)
 109        format(1x,i5,2x,i5,34(2x,i2))
            jpatt(jpati)=0
            endif
   7     continue
c
c       compute the density of times between interviews with and 
c       without using month and date.
c
      do 48 i=1,2
         do 41 j=1,ntime
  41        dens(i,j)=0.d0
  48     continue
      do 38 i=1,nobs
         if(incluv(i).eq.1)then
            iflag=0
            do 39 j=1,ntime
               if(xv(i,j,1).eq.1.)then
                  if(iflag.eq.0)then
                     timet(1,1)=adjtime(j)
                     month=idatev(i,j)/100
                     amonth=dble(float(month))/12.d0
                     timet(2,1)=adjtime(j)+amonth
                     iflag=1
                     goto 39
                     endif
                  if(iflag.eq.1)then
                     timet(1,2)=adjtime(j)
                     month=idatev(i,j)/100
                     amonth=dble(float(month))/12.d0
                     timet(2,2)=adjtime(j)+amonth
                     dtime=timet(1,2)-timet(1,1)
                     itime=ifix(sngl(dtime*12.d0))
                     dens(1,itime)=dens(1,itime)+1.d0
                     dtime=timet(2,2)-timet(2,1)
                     itime=ifix(sngl(dtime*12.d0))
                     dens(2,itime)=dens(2,itime)+1.d0
                     timet(1,1)=timet(1,2)
                     timet(2,1)=timet(2,2)
                     endif
                  endif
  39           continue
            endif
  38     continue
      write(6,100)
      write(6,114)
 114  format(1x,'density of times between interviews with and'
     & ' without month and date information')
      write(6,115)
 115  format(15x,'density',2x,'density',/,15x,'without',2x,'with',/,
     & 1x,'years',2x,'month',2x,'months',3x,'months')
      ndens=ntime*12
      do 43 i=1,2
         denstot=0.d0
         do 44 j=1,ndens
  44        denstot=denstot+dens(i,j)
         do 45 j=1,ndens
  45        dens(i,j)=dens(i,j)/denstot
  43     continue
      do 42 i=1,ndens
         idens(1)=i/12
         idens(2)=i-(idens(1)*12)
  42     write(6,116)(idens(j),j=1,2),(dens(j,i),j=1,2)
 116  format(1x,i3,4x,i2,5x,f5.3,2x,f5.3)
c
c       output number of observations for each wave.
c
      write(6,100)
      write(6,129)
 129  format(1x,'number of observations for each wave')
      write(6,132)
 132  format(1x,'wave',2x,'# obs')
      do 57 i=1,ntime
  57     write(6,133)i,iwave(i)
 133  format(1x,i2,4x,i5)
c
c       output descriptive statistics for each dependent variable
c       including person-specific standard deviation and idiosyncratic
c       standard deviation.
c
      write(6,100)
      write(6,130)
 130  format(1x,'descriptive statistics for dependent variables')
      write(6,128)
 128  format(67x,'person-',3x,'idiosyn-',/,67x,'specific',2x,'cratic',
     & /,1x,'variable',18x,'# obs',5x,'mean',6x,'std dev',13x,
     & 'std dev',3x,'std dev',3x,'minimum',3x,'maximum')
      do 58 i=1,ny
         ie=i*3
         ib=ie-2
         ymom(i,2)=ymom(i,2)/ymom(i,1)
         ymom(i,3)=((ymom(i,3)/ymom(i,1))-(ymom(i,2)*ymom(i,2)))**.5d0
         ymom(i,5)=(ymom(i,5)/ymom(i,4))**.5d0
         ymom(i,6)=((ymom(i,3)*ymom(i,3))-(ymom(i,5)*ymom(i,5)))**.5d0
  58     write(6,121)(ayvar(j),j=ib,ie),(ymom(i,j),j=1,8)
 121  format(1x,3a8,2x,f8.0,2(2x,f8.3),2x,f8.0,2(2x,f8.3),2(2x,f8.1))
c
c       output correlation matrix of contemporaneous dependent 
c       variables.
c
      write(6,100)
      write(6,124)(i,i=1,ny)
 124  format(2x,6(2x,i1,5x))
      do 83 i=1,ny
  83     cmom1(i,2)=cmom1(i,2)/cmom1(i,1)
      do 59 i=1,ny
         do 60 j=1,i
            cov(i,j)=(cmom2(i,j,2)/cmom2(i,j,1))-(cmom1(i,2)*
     &       cmom1(j,2))
  60        cov(j,i)=cov(i,j)
  59     continue
         do 61 i=1,ny
            do 62 j=1,i
  62           corr(i,j)=cov(i,j)/((cov(i,i)*cov(j,j))**.5d0)
  61     write(6,122)i,(corr(i,j),j=1,i)
 122  format(1x,i1,6(2x,f6.3))
c
c       output first autocovariance for each dependent variable.
c
      write(6,100)
      write(6,125)
 125  format(1x,'first autocovariances for each dependent variable')
      write(6,127)
 127  format(1x,'variable',18x,'autocovariance')
      do 63 i=1,ny
         ie=i*3
         ib=ie-2
         autocov(i,2)=autocov(i,2)/autocov(i,1)
         autocov(i,3)=(autocov(i,3)/autocov(i,1))-(autocov(i,2)*
     &    autocov(i,2))
  63     write(6,126)(ayvar(j),j=ib,ie),autocov(i,3)
 126  format(1x,3a8,2x,f8.3)
c
c       output lagged covariance table.
c
      write(6,100)
      write(6,135)
 135  format(1x,'lagged covariance matrix')
      write(6,136)(i,i=1,ny)
 136  format(11x,'lagged',/,11x,6(i2,8x))
      do 84 i=1,ny
         ameanlag(i,1,2)=ameanlag(i,1,2)/ameanlag(i,1,1)
  84     ameanlag(i,2,2)=ameanlag(i,2,2)/ameanlag(i,2,1)
      do 85 i=1,ny
         do 86 j=1,ny
  86        covlag(i,j,2)=(covlag(i,j,2)/covlag(i,j,1))-
     &       (ameanlag(i,1,2)*ameanlag(j,2,2))
  85     write(6,134)i,(covlag(i,j,2),j=1,ny)
 134  format(1x,i3,7x,6(f8.3,2x))
c
c       output lagged correlation table.
c
      write(6,100)
      write(6,137)
 137  format(1x,'lagged correlation matrix')
      write(6,136)(i,i=1,ny)
      do 94 i=1,ny
         do 98 j=1,2
  98        varlag(i,j)=(varlag(i,j)/ameanlag(i,j,1))-(ameanlag(i,j,2)
     &       *ameanlag(i,j,2))
  94     continue
      do 200 i=1,ny
         do 93 j=1,ny
  93        corrlag(i,j)=covlag(i,j,2)/((varlag(i,1)*varlag(j,2))
     &       **.5d0)
 200     write(6,134)i,(corrlag(i,j),j=1,ny)
      return
      end
c
c
      subroutine errget1(ir,ia,e,u)
c
c       this subroutine gets errors for the simple model.
c
      parameter(nbetas=20,nfacs=2,nrs=10,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension u(ntimes)
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
!$omp threadprivate(/datat/)
      if(ia.eq.1)then
         e=er(ir)
         do 2 i=istart,iend 
   2        u(i)=ur(ir,1,i)
         endif
      if(ia.eq.2)then
         e=-1.d0*e
         do 3 i=istart,iend 
   3        u(i)=1.d0-u(i)
         endif
      return
      end
c
c
      subroutine errget2(ir,ia,e,u)
c
c       this subroutine gets errors for the complex model.
c
      parameter(nbetas=20,nfacs=2,nrs=10,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension e(nfacs),u(nys,ntimes)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
!$omp threadprivate(/datat/)
      if(ia.eq.1)then
         do 5 i=1,ny
            do 6 j=istart,iend
   6           u(i,j)=ur(ir,i,j)
   5        continue
         do 2 i=1,nfac
   2        e(i)=er2(ir,i)
         endif
      if(ia.eq.2)then
         do 3 i=1,nfac
   3        e(i)=-1.d0*e(i)
         do 7 i=1,ny
            do 8 j=istart,iend
   8           u(i,j)=-1.d0*u(i,j)
   7        continue
         endif
      return
      end
c
c
      subroutine errmake
c
c       this subroutine draws errors.
c
      parameter(nfacs=2,nobss=96529,nrs=10,ntimes=34,nys=6)
      implicit real*8(a-h,o-z)
      character*8 adist(2)
      real e,erv,randfn,u,urv
      dimension usum(2,5)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/errors/erv(nobss,nfacs,nrs),urv(nobss,nrs,nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      data adist/'uniform','normal'/
      iseed=1
      u=randfn(iseed)
      iseed=0
      do 8 i=1,2
         do 9 j=1,3
   9        usum(i,j)=0.d0
         usum(i,4)=1.d20
   8     usum(i,5)=-1.d20
c
c       draw erv.
c
      if(model.eq.1)nfact=1
      if(model.eq.2)nfact=nfac
      write(6,103)
 103  format(1x,'draws in errmake:')
      write(6,104)
 104  format(1x,'draws of erv:')
      do 2 i=1,nobs
         if(((i/10000)*10000).eq.i)write(6,105)i
 105     format(1x,'i = ',i8)
         do 3 j=1,nr
            do 10 k=1,nfact
               u=randfn(iseed)
               iseed=0
               call utrunc(u)
               ut=dble(u)
               usum(1,1)=usum(1,1)+1.d0
               usum(1,2)=usum(1,2)+ut
               usum(1,3)=usum(1,3)+(ut*ut)
               if(usum(1,4).gt.ut)usum(1,4)=ut
               if(usum(1,5).lt.ut)usum(1,5)=ut
               call mdnris(u,e,ier)
               erv(i,k,j)=e
               et=dble(e)
               usum(2,1)=usum(2,1)+1.d0
               usum(2,2)=usum(2,2)+et
               usum(2,3)=usum(2,3)+(et*et)
               if(usum(2,4).gt.et)usum(2,4)=et
               if(usum(2,5).lt.et)usum(2,5)=et
  10        continue
   3        continue
   2     continue
c
c       draw urv.
c
      write(6,106)
 106  format(1x,'draws of urv')
      do 4 i=1,nobs
         if(((i/10000)*10000).eq.i)write(6,105)i
         do 5 j=1,nr
            do 6 k=1,ntime
               do 11 l=1,ny
                  u=randfn(iseed)
                  iseed=0
                  call utrunc(u)
                  ut=dble(u)
                  usum(1,1)=usum(1,1)+1.d0
                  usum(1,2)=usum(1,2)+ut
                  usum(1,3)=usum(1,3)+(ut*ut)
                  if(usum(1,4).gt.ut)usum(1,4)=ut
                  if(usum(1,5).lt.ut)usum(1,5)=ut
                  call mdnris(u,e,ier)
                  urv(i,j,l,k)=e
                  et=dble(e)
                  usum(2,1)=usum(2,1)+1.d0
                  usum(2,2)=usum(2,2)+et
                  usum(2,3)=usum(2,3)+(et*et)
                  if(usum(2,4).gt.et)usum(2,4)=et
                  if(usum(2,5).lt.et)usum(2,5)=et
  11              continue
   6           continue
   5        continue
   4     continue
c
c       adjust and output moments.
c
      write(6,100)
 100  format(1x,'moments of draws')
      do 7 i=1,2
         write(6,101)adist(i)
 101     format(1x,'moments of ',a8,' draws')
         usum(i,2)=usum(i,2)/usum(i,1)
         usum(i,3)=((usum(i,3)/usum(i,1))-(usum(i,2)*usum(i,2)))**.5d0
         write(6,102)(usum(i,j),j=1,5)
 102     format(1x,'# obs: ',f10.0,/,1x,'mean: ',f8.4,/,1x,'std dev: ',
     &    f8.4,/,1x,'minimum: ',f8.4,/,1x,'maximum: ',f8.4)
   7     continue
      return
      end
c
c
      subroutine estim1
c
c       this subroutine is the driver for estimating parameters using
c       gqopt.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+(nys*nsegm2s)+
     & (nys*(nfacs+nys))+(nys*(nys-1)/2),naints=(ncoefs*ncoefs)+
     & (4*ncoefs)+1)
      implicit real*8(a-h,o-z)
      character*8 alabel
      dimension parm(ncoefs),parmold(ncoefs)
      common/bnel/stp1,var,konvge,nrst
      common/bstack/aint(naints)
      common/bstak/nq
      common/labels/alabel(ncoefs)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      external likely11,likely12,nmsimp
      data dobjtol/0.1d-3/
      data dparm2tol/0.1d-3/
c
c       initialize.
c
      objold=0.d0
      call bestupd(0,objold)
      call parset0(parm,nparm)
      do 2 i=1,nparm
   2     parmold(i)=parm(i)
      if(model.eq.1)call likely11(parmold,nparm,objold,*3)
      if(model.eq.2)call likely12(parmold,nparm,objold,*3)
   3  continue
      anparm=dble(float(nparm))
      mx=1    
      acc=.1d-4
      stp1h=stp1
      nq=(nparm*nparm)+(4*nparm)+1
      do 4 i=1,8
         stp1=stp1h
c
c       call gqopt.
c
         if(model.eq.1)call opt(parm,nparm,obj,nmsimp,niter,mx,ier,acc,
     &    likely11,alabel)
         if(model.eq.2)call opt(parm,nparm,obj,nmsimp,niter,mx,ier,acc,
     &    likely12,alabel)
         dobj=log(objold/obj)
         objold=obj
         dparm2=0.d0
         do 5 j=1,nparm
            dparmt=parm(j)-parmold(j)
            dparm2=dparm2+(dparmt*dparmt)
   5        parmold(j)=parm(j)
         dparm2=(dparm2/anparm)**.5d0
         write(6,100)i,obj,dobj,dparm2
 100     format(1x,'iteration = ',i1,2x,'obj = ',g15.8,2x,'dobj = ',
     &    g15.8,2x,'dparm2 = ',g15.8)
         if((dobj.lt.dobjtol).and.(dparm2.lt.dparm2tol))then
            write(6,101)
 101        format(1x,'convergence reached')
            call covar
            call vardecomp
            return
            endif
   4     continue
c      write(6,102)
c 102  format(1x,'no convergence')
      call covar
      call vardecomp
      return
      end
c
c
      subroutine estim2
c
c       this subroutine estimates the parameters of the model using
c       modified bhhh.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+(nys*nsegm2s)+(nys*
     & (nfacs+nys))+(nys*(nys-1)/2),ncoef2s=2*ncoefs)
      implicit real*8(a-h,o-z)
      dimension delta(ncoefs),df(ncoefs),df2(ncoefs,ncoefs),
     & parm(ncoefs),parmold(ncoefs),work(ncoefs),work2(ncoef2s)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      data dobjtol/.1d-2/
      data dparm2tol/.1d-2/
      call cpu_time(time0)
      iter=1
   9  call parset0(parm,nparm)
      anparm=dble(float(nparm))
      call likely2(parm,nparm,f,df,df2,1)
      do 4 i=1,nparm
   4     parmold(i)=parm(i)
      fold=f
      d1=-1.d0
      call linv3f(df2,work,1,nparm,ncoefs,d1,d2,work2,ier)
      do 2 i=1,nparm
         parm(i)=0.d0
         do 3 j=1,nparm
   3        parm(i)=parm(i)+(df2(i,j)*df(j))
   2     delta(i)=parm(i)-parmold(i)
      iterg=-1
      iterb=-1
   6  continue
      if((iterg.gt.niterg).or.(iterb.gt.niterb))goto 7
      call likely2(parm,nparm,f,df,df2,0)
      if(f.gt.fold)then
         do 5 i=1,nparm
            parmold(i)=parm(i)
            delta(i)=delta(i)/2.d0
   5        parm(i)=parm(i)+delta(i)
         if(iterb.gt.0)goto 7
         iterg=iterg+1
         go to 6
         endif
      if(f.le.fold)then
         if(iterg.gt.0)goto 7
         do 8 i=1,nparm
            delta(i)=delta(i)/2.d0
   8        parm(i)=parmold(i)+delta(i)
         iterb=iterb+1
         goto 6
         endif
   7  continue
      if(iterb.gt.niterb)then
         write(6,100)niterb,f,fold
 100     format(1x,'no improvement after ',i2,' half-steps: f = ',
     &    g15.8,2x,'fold = ',g15.8)
         call cpu_time(time1)
         dtime=time1-time0
         write(6,103)dtime
 103     format(1x,'time = 'g15.8)
         call covar
         stop
         endif
      dobj=log(fold/f)
      dparm2=0.d0
      do 10 i=1,nparm
         dparm2t=parm(i)-parmold(i)
  10     dparm2=dparm2+(dparm2t*dparm2t)
      dparm2=(dparm2/anparm)**.5d0
      if((iter.gt.niter).or.((dobj.lt.dobjtol).and.
     & (dparm2.lt.dparm2tol)))then
         if(iter.gt.niter)write(6,101)
 101     format(1x,'hit maximum # iterations')
         if((dobj.lt.dobjtol).and.(dparm2.lt.dparm2tol))write(6,102)
 102     format(1x,'convergence')
         call cpu_time(time1)
         dtime=time1-time0
         write(6,103)dtime
         call covar
         endif
      if((iter.le.niter).and.((dobj.ge.dobjtol).or.
     & (dparm2.gt.dparm2tol)))then
         iter=iter+1
         goto 9
         endif
      return
      end
c
c
      subroutine ghk1(iyc,omega,ystar,adjmean,adjstdev)
c
c       this subroutine provides ghk step information for likobs12.
c
      parameter(nys=6,ny2s=nys*2)
      implicit real*8(a-h,o-z)
      dimension a(nys,nys),aone(nys),omega(nys,nys),work(nys),
     & work2(ny2s),ystar(nys)
      if(iyc.eq.1)then
c
c       initialize for first element.
c
         adjmean=0.d0
         adjstdev=omega(1,1)**.5d0
         return
         endif
c
c       otherwise, invert omega.
c
      do 2 i=1,iyc
         do 3 j=1,iyc
   3        a(i,j)=omega(i,j)
   2     continue
      d1=-1.d0
      call linv3f(a,work,1,iyc,nys,d1,d2,work2,ier)
c
c       multiply a11*a21T.
c
      iycm1=iyc-1
      do 4 i=1,iycm1     
         do 5 j=1,iycm1
   5        aone(i)=a(1,1)*a(j+1,1)
   4     continue
      adjmean=0.d0
      do 6 i=1,iycm1
   6     adjmean=adjmean-(aone(i)*ystar(i))
      adjstdev=a(1,1)**(-.5d0)
      return
      end
c
c
      subroutine ghk2(iyc,omega,domega,ystar,dystar,adjmean,dadjmean,
     & adjstdev,dadjstdev,nparm,ifder)
c
c       this subroutine provides ghk step information for likobs2.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ny2s=nys*2,ncoefs=(nys*ntimes)+(nys*nbetas)+(nys*nsegm2s)+
     & (nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension a(nys,nys),aone(nys),da(nys,nys,ncoefs),
     & dadjmean(ncoefs),dadjstdev(ncoefs),daone(nys,ncoefs),
     & dat(nys,nys,ncoefs),domega(nys,nys,ncoefs),
     & dystar(nys,ncoefs),omega(nys,nys),work(nys),work2(ny2s),
     & ystar(nys)
      if(iyc.eq.1)then
c
c       initialize for first element.
c
         adjmean=0.d0
         adjstdev=omega(1,1)**.5d0
         return
         endif
c
c       otherwise, invert omega.
c
      do 2 i=1,iyc
         do 3 j=1,iyc
   3        a(i,j)=omega(i,j)
   2     continue
      if(ifder.eq.1)then
         do 11 i=1,iyc
            do 12 j=1,iyc
               do 13 k=1,nparm
  13              dat(i,j,k)=domega(i,j,k)
  12           continue
  11        continue
         endif
      d1=-1.d0
      call linv3f(a,work,1,iyc,nys,d1,d2,work2,ier)
      if(ifder.eq.1)then
         do 14 i=1,nparm
            do 15 j=1,iyc
               do 16 k=1,iyc
                  da(j,k,i)=0.d0
                  do 17 l=1,iyc
                     do 18 m=1,iyc
  18                    da(j,k,i)=da(j,k,i)+(a(j,l)*dat(l,m,i)*a(k,m))
  17                 continue
  16              continue
  15           continue
  14        continue
         endif
c
c       multiply a11*a21T.
c
      iycm1=iyc-1
      do 4 i=1,iycm1     
         do 5 j=1,iycm1
   5        aone(i)=a(1,1)*a(j+1,1)
   4     continue
      if(ifder.eq.1)then
         do 19 i=1,nparm
            do 20 j=1,iycm1
               do 21 k=1,iycm1
  21              daone(j,i)=(da(1,1,i)*a(k+1,1))+(a(1,1)*da(k+1,1,i))
  20           continue
  19        continue
         endif
      adjmean=0.d0
      do 6 i=1,iycm1
   6     adjmean=adjmean-(aone(i)*ystar(i))
      adjstdev=a(1,1)**(-.5d0)
      if(ifder.eq.1)then
         do 10 i=1,nparm
  10        dadjmean(i)=0.d0
         do 7 i=1,iycm1
            do 8 j=1,nparm
   8           dadjmean(j)=dadjmean(j)-(daone(i,j)*ystar(i))+
     &          (aone(i)*dystar(i,j))
   7        continue
         do 9 i=1,nparm
   9        dadjstdev(i)=-.5d0*da(1,1,i)/(a(1,1)**.5d0)
         endif
      return
      end
c
c
      subroutine initia
c
c       this subroutine is the driver for initialization.
c
      implicit real*8(a-h,o-z)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
c
c       get parameters.
c
      call parget
c
c       get data.
c
      if(idata.eq.1)call datget1
      if(idata.eq.2)call datget2
      if(idata.eq.3)call datget3
c
c       draw errors.
c
      call errmake
      return
      end
c
c
      subroutine likely11(parm,nparm,f,*)
c
c       this subroutine simulates the log likelihood function for the 
c       simple model.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      real xextra,xv
      dimension parm(ncoefs)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      call parset1(parm)
      call parset211
      f=0.d0
      do 2 i=1,nobs
         if(incluv(i).eq.1)then
            call obsmake1(i)
            call likobs11(ft)
            f=f+ft
            endif
   2     continue
      call bestupd(1,f)
      return
      end
c
c
      subroutine likely12(parm,nparm,f,*)
c
c       this subroutine simulates the log likelihood function for the 
c       complex model.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      real xextra,xv
      dimension parm(ncoefs)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      call parset1(parm)
      call parset212(penalty)
      f=0.d0
!$omp parallel do private (i,ft) reduction (+:f)
      do 2 i=1,nobs
         if(incluv(i).eq.1)then
            call obsmake2(i)
            call likobs12(ft)
            f=f+ft
            endif
   2     continue
!$omp end parallel do
      f=f-penalty
      call bestupd(1,f)
      return
      end
c
c
      subroutine likely2(parm,nparm,f,df,df2,ifder)
c
c       this subroutine simulates the log likelihood function for the 
c       complex model along with first partial derivatives.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      real xextra,xv
      dimension df(ncoefs),dft(ncoefs),df2(ncoefs,ncoefs),
     & dpenalty(ncoefs),parm(ncoefs)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      call parset1(parm)
      call parset22(penalty,dpenalty,nparm,ifder)
      f=0.d0
      if(ifder.eq.1)then
         do 4 i=1,nparm
   4        df(i)=0.d0
         do 8 i=1,nparm
            do 9 j=1,nparm
   9           df2(i,j)=0.d0
   8        continue
         endif
      do 2 i=1,nobs
         if(incluv(i).eq.1)then
            call obsmake2(i)
            call likobs2(ft,dft,nparm,ifder)
            f=f+ft
            if(ifder.eq.1)then
               do 5 j=1,nparm
   5              df(j)=df(j)+dft(j)
               do 6 j=1,nparm
                  do 7 k=1,nparm
   7                 df2(j,k)=df2(j,k)+(dft(j)*dft(k))
   6              continue
               endif
            endif
   2     continue
      f=f-penalty
      if(ifder.eq.1)then
         do 3 i=1,nparm
   3        df(i)=df(i)-dpenalty(i)
         endif
      return
      end
c
c
      subroutine likobs11(f)
c
c       this subroutine simulates the log likelihood contribution for
c       observation iobs for the simple model.
c
      parameter(nbetas=20,nfacs=2,nrs=10,nsegs=8,nsegm2s=nsegs-2,
     & nsegp1s=nsegs+1,ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+
     & (nys*nbetas)+(nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension ftt(nrs,2),u(ntimes)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
c
c       prepare non-random terms.
c
      call preplik11
c
c       loop over draws.
c
      do 2 i=1,nr
         do 3 j=1,2
            call errget1(i,j,e,u)
            ftt(i,j)=0.d0
            iflag=0
            do 4 k=istart,iend 
               if((iobsrv(k).eq.1).and.(iy(1,k).ge.0).and.
     &          (iy(1,k).le.nseg(1)))then
                  if(iflag.eq.0)then
                     if(ifcond1.eq.0)then
                        ystar(1,k)=xb(1,k)+(sige*e)
                        stdev=((1.d0-(rhop(1,1,k)*rhop(1,1,k)))**.5d0)
                        arg1=(thresh(1,iy(1,k))-ystar(1,k))/stdev
                        arg2=(thresh(1,iy(1,k)+1)-ystar(1,k))/stdev
                        call probmake1(arg1,arg2,fttt)
                        ftt(i,j)=ftt(i,j)+fttt
                        iyk=iy(1,k)
                        ystark=ystar(1,k)
                        uk=u(k)
                        call normadj1(iyk,1,ystark,uk,utc)
                        ystar(1,k)=ystar(1,k)+utc
                        endif
                     if(ifcond1.eq.1)then
                        ystark=xb(1,k)+(sige*e)
                        uk=u(k)
                        iyk=iy(1,k)
                        call normadj1(iyk,1,ystark,uk,utc)
                        ystar(1,k)=ystark+utc
                        iflag=1
                        goto 4
                        endif
                     endif
                  if(iflag.eq.1)then
                     ystar(1,k)=(rhop(1,1,k)*ystar(1,k-1))+xb(1,k)+
     &                (sige*e)
                     arg1=thresh(1,iy(1,k))-ystar(1,k)
                     arg2=thresh(1,iy(1,k)+1)-ystar(1,k)
                     call probmake1(arg1,arg2,fttt)
                     ftt(i,j)=ftt(i,j)+fttt
                     iyk=iy(1,k)
                     ystark=ystar(1,k)
                     uk=u(k)
                     call normadj1(iyk,1,ystark,uk,utc)
                     ystar(1,k)=ystar(1,k)+utc
                     endif
                  endif
               if((iobsrv(k).eq.0).or.(iy(1,k).lt.0).or.
     &          (iy(1,k).gt.nseg(1)))then
                  if(iflag.eq.1)then
                     ystar(1,k)=(rhop(1,1,k)*ystar(1,k-1))+xb(1,k)+
     &                (sige*e)+u(k)
                     endif
                  endif
   4           continue
   3        continue
   2     continue
c
c       find maximum.
c
      fttmax=-1.d20
      do 7 i=1,nr
         do 8 j=1,2
            if(ftt(i,j).gt.fttmax)fttmax=ftt(i,j)
   8        continue
   7     continue
      f=0.d0
      do 9 i=1,nr
         do 10 j=1,2
  10        f=f+exp(ftt(i,j)-fttmax)
   9     continue
      f=fttmax+log(f)-log(dble(float(nr*2)))
      return
      end
c
c
      subroutine likobs12(f)
c
c       this subroutine simulates the log likelihood contribution for
c       observation iobs for the complex model.
c
      parameter(nbetas=20,nfacs=2,nrs=10,nsegs=8,nsegm2s=nsegs-2,
     & nsegp1s=nsegs+1,ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+
     & (nys*nbetas)+(nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension e(nfacs),ftt(nrs,2),u(nys,ntimes),ystart(nys)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
!$omp threadprivate(/datat/)
c
c       prepare non-random terms.
c
      call preplik12
c
c       loop over draws.
c
      do 2 i=1,nr
         do 3 j=1,2
            call errget2(i,j,e,u)
            ftt(i,j)=0.d0
            do 4 k=istart,iend 
               km1=k-1
               do 11 l=1,ny
                  if((iobsrv(k).eq.1).and.(iy(l,k).ge.1).and.
     &             (iy(l,k).le.nseg(l)))then
                     ltt=l
                     mt=l+1
                     do 12 m=1,l
                        mt=mt-1
  12                    ystart(m)=ystar(mt,k)
                     call ghk1(ltt,omega,ystart,adjmean,adjstdev)
                     ystar(l,k)=xb(l,k)
                     if(km1.ge.istart)then
                        do 14 m=1,ny
  14                       ystar(l,k)=ystar(l,k)+(rhop(l,m,k)*
     &                      ystar(m,km1))
                        endif
                     do 13 m=1,nfac
  13                    ystar(l,k)=ystar(l,k)+(alamb(l,m)*e(m))
                     arg1=(thresh(l,iy(l,k))-ystar(l,k)-adjmean)/
     &                adjstdev
                     arg2=(thresh(l,iy(l,k)+1)-ystar(l,k)-adjmean)/
     &                adjstdev
                     call probmake1(arg1,arg2,fttt)
                     ftt(i,j)=ftt(i,j)+fttt
                     iyk=iy(l,k)
                     ystark=ystar(l,k)
                     uk=u(l,k)
                     call normadj1(iyk,l,ystark,uk,utc)
                     ystar(l,k)=ystar(l,k)+utc
                     endif
                  if((iobsrv(k).eq.0).or.(iy(l,k).lt.1).or.
     &             (iy(l,k).gt.nseg(l)))then
                     ystar(l,k)=xb(l,k)+u(l,k)
                     if(km1.ge.istart)then
                        do 15 m=1,ny
  15                       ystar(l,k)=ystar(l,k)+(rhop(l,m,k)*
     &                      ystar(m,km1))
                        endif
                     do 16 m=1,nfac
  16                    ystar(l,k)=ystar(l,k)+(alamb(l,m)*e(m))
                     endif
  11              continue
   4           continue
   3        continue
   2     continue
c
c       find maximum.
c
      fttmax=-1.d20
      do 7 i=1,nr
         do 8 j=1,2
            if(ftt(i,j).gt.fttmax)fttmax=ftt(i,j)
   8        continue
   7     continue
      f=0.d0
      do 9 i=1,nr
         do 10 j=1,2
  10        f=f+exp(ftt(i,j)-fttmax)
   9     continue
      f=fttmax+log(f)-log(dble(float(nr*2)))
      return
      end
c
c
      subroutine likobs12t(f)
c
c       this subroutine simulates the log likelihood contribution for
c       observation iobs for the complex model with extra parameters
c       for lagrange multiplier tests.
c
      parameter(nbetas=20,nfacs=2,nrs=10,nsegs=8,nsegm2s=nsegs-2,
     & nsegp1s=nsegs+1,ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+
     & (nys*nbetas)+(nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension e(nfacs),ftt(nrs,2),u(nys,ntimes),ystart(nys)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/lmparms/alambt(nys,nfacs),omegat(nys,nys),
     & rhoptt(nys,nys,ntimes),rhot(nys,nys),iswitch(ntimes),lmflag
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
c
c       prepare non-random terms.
c
      call preplik12t
c
c       loop over draws.
c
      do 2 i=1,nr
         do 3 j=1,2
            call errget2(i,j,e,u)
            ftt(i,j)=0.d0
            do 4 k=istart,iend 
               km1=k-1
               do 11 l=1,ny
                  if((iobsrv(k).eq.1).and.(iy(l,k).ge.1).and.
     &             (iy(l,k).le.nseg(l)))then
                     ltt=l
                     mt=l+1
                     do 12 m=1,l
                        mt=mt-1
  12                    ystart(m)=ystar(mt,k)
                     if(iswitch(k).eq.0)call ghk1(ltt,omega,ystart,
     &                adjmean,adjstdev)
                     if(iswitch(k).eq.1)call ghk1(ltt,omegat,ystart,
     &                adjmean,adjstdev)
                     ystar(l,k)=xb(l,k)
                     if(km1.ge.istart)then
                        if(iswitch(k).eq.0)then
                           do 14 m=1,ny
  14                          ystar(l,k)=ystar(l,k)+(rhop(l,m,k)*
     &                         ystar(m,km1))
                           endif
                        if(iswitch(k).eq.1)then
                           do 5 m=1,ny
   5                          ystar(l,k)=ystar(l,k)+(rhoptt(l,m,k)*
     &                         ystar(m,km1))
                           endif
                        endif
                     if(iswitch(k).eq.0)then
                        do 13 m=1,nfac
  13                       ystar(l,k)=ystar(l,k)+(alamb(l,m)*e(m))
                        endif
                     if(iswitch(k).eq.1)then
                        do 6 m=1,nfac
   6                       ystar(l,k)=ystar(l,k)+(alambt(l,m)*e(m))
                        endif
                     arg1=(thresh(l,iy(l,k))-ystar(l,k)-adjmean)/
     &                adjstdev
                     arg2=(thresh(l,iy(l,k)+1)-ystar(l,k)-adjmean)/
     &                adjstdev
                     call probmake1(arg1,arg2,fttt)
                     ftt(i,j)=ftt(i,j)+fttt
                     iyk=iy(l,k)
                     ystark=ystar(l,k)
                     uk=u(l,k)
                     call normadj1(iyk,l,ystark,uk,utc)
                     ystar(l,k)=ystar(l,k)+utc
                     endif
                  if((iobsrv(k).eq.0).or.(iy(l,k).lt.1).or.
     &             (iy(l,k).gt.nseg(l)))then
                     ystar(l,k)=xb(l,k)+u(l,k)
                     if(km1.ge.istart)then
                        if(iswitch(k).eq.0)then
                           do 15 m=1,ny
  15                          ystar(l,k)=ystar(l,k)+(rhop(l,m,k)*
     &                         ystar(m,km1))
                           endif
                        if(iswitch(k).eq.1)then
                           do 17 m=1,ny
  17                          ystar(l,k)=ystar(l,k)+(rhoptt(l,m,k)*
     &                         ystar(m,km1))
                           endif
                        endif
                     if(iswitch(k).eq.0)then
                        do 16 m=1,nfac
  16                       ystar(l,k)=ystar(l,k)+(alamb(l,m)*e(m))
                        endif
                     if(iswitch(k).eq.1)then
                        do 18 m=1,nfac
  18                       ystar(l,k)=ystar(l,k)+(alambt(l,m)*e(m))
                        endif
                     endif
  11              continue
   4           continue
   3        continue
   2     continue
c
c       find maximum.
c
      fttmax=-1.d20
      do 7 i=1,nr
         do 8 j=1,2
            if(ftt(i,j).gt.fttmax)fttmax=ftt(i,j)
   8        continue
   7     continue
      f=0.d0
      do 9 i=1,nr
         do 10 j=1,2
  10        f=f+exp(ftt(i,j)-fttmax)
   9     continue
      f=fttmax+log(f)-log(dble(float(nr*2)))
      return
      end
c
c
      subroutine likobs2(f,df,nparm,ifder)
c
c       this subroutine simulates the log likelihood contribution and
c       analytical derivatives for observation iobs for the complex 
c       model.
c
      parameter(nbetas=20,nfacs=2,nrs=10,nsegs=8,nsegm2s=nsegs-2,
     & nsegp1s=nsegs+1,ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+
     & (nys*nbetas)+(nys*nsegm2s)+nys*(nfacs+nys)+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension dadjmean(ncoefs),dadjstdev(ncoefs),darg1(ncoefs),
     & darg2(ncoefs),df(ncoefs),dfac(ncoefs),dftt(nrs,2,ncoefs),
     & dfttmax(ncoefs),dfttt(ncoefs),dutc(ncoefs),
     & dystar(nys,ntimes,ncoefs),dystark(ncoefs),
     & dystart(nys,ncoefs),e(nfacs),ftt(nrs,2),ifttmax(2),
     & u(nys,ntimes),ystart(nys),ystartt(nys,ntimes)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
c
c       prepare non-random terms.
c
      call preplik2(nparm,ifder)
c
c       loop over draws.
c
      do 2 i=1,nr
         do 3 j=1,2
            call errget2(i,j,e,u)
            ftt(i,j)=0.d0
            if(ifder.eq.1)then
               do 22 k=1,nparm
  22              dftt(i,j,k)=0.d0
               endif
            do 4 k=istart,iend
               km1=k-1
               do 11 l=1,ny
                  if((iobsrv(k).eq.1).and.(iy(l,k).ge.1).and.
     &             (iy(l,k).le.nseg(l)))then
                     ltt=l
                     mt=l+1
                     do 12 m=1,l
                        mt=mt-1
  12                    ystart(m)=ystar(mt,k)
                     if(ifder.eq.1)then
                        do 29 m=1,l
                           do 30 n=1,nparm
  30                          dystart(m,n)=dystar(mt,k,n)
  29                       continue
                        endif
                     call ghk2(ltt,omega,domega,ystart,dystart,adjmean,
     &                dadjmean,adjstdev,dadjstdev,nparm,ifder)
                     ystar(l,k)=xb(l,k)
                     if(ifder.eq.1)then
                        do 28 m=1,nparm
  28                       dystar(l,k,m)=dxb(l,k,m)
                        endif
                     if(km1.ge.istart)then
                        if(ifder.eq.1)then
                           do 26 m=1,ny
                              do 27 n=1,nparm
  27                             dystar(l,k,n)=dystar(l,k,n)+
     &                            (drhop(l,m,k,n)*ystar(m,km1))+
     &                            (rhop(l,m,k)*dystar(m,km1,n))
  26                          continue
                           endif
                        do 14 m=1,ny
  14                       ystar(l,k)=ystar(l,k)+(rhop(l,m,k)*
     &                      ystar(m,km1))
                        endif
                     do 13 m=1,nfac
                        ystar(l,k)=ystar(l,k)+(alamb(l,m)*e(m))
                        if(ifder.eq.1)then
                           do 25 n=1,nparm
  25                          dystar(l,k,n)=dystar(l,k,n)+
     &                         (dalamb(l,m,n)*e(m))
                           endif
  13                    continue
                     arg1=(thresh(l,iy(l,k))-ystar(l,k)-adjmean)/
     &                adjstdev
                     arg2=(thresh(l,iy(l,k)+1)-ystar(l,k)-adjmean)/
     &                adjstdev
                     if(ifder.eq.1)then
                        do 23 m=1,nparm
  23                       darg1(m)=((dthresh(l,iy(l,k),m)-
     &                      dystar(l,k,m)-dadjmean(m))/adjstdev)-
     &                      (arg1*dadjstdev(m)/adjstdev)
                        do 24 m=1,nparm
  24                       darg2(m)=((dthresh(l,iy(l,k)+1,m)-
     &                      dystar(l,k,m)-dadjmean(m))/adjstdev)-
     &                      (arg2*dadjstdev(m)/adjstdev)
                        endif
                     call probmake2(arg1,darg1,arg2,darg2,fttt,dfttt,
     &                nparm,ifder)
                     ftt(i,j)=ftt(i,j)+fttt
                     if(ifder.eq.1)then
                        do 21 m=1,nparm
  21                       dftt(i,j,m)=dftt(i,j,m)+dfttt(m)
                        endif
                     iyk=iy(l,k)
                     ystark=ystar(l,k)
                     if(ifder.eq.1)then
                        do 31 m=1,nparm
  31                       dystark(m)=dystar(l,k,m)
                        endif
                     uk=u(l,k)
                     call normadj2(iyk,l,ystark,dystark,uk,utc,dutc,
     &                nparm,ifder)
                     ystar(l,k)=ystar(l,k)+utc
                     if(ifder.eq.1)then
                        do 38 m=1,nparm
  38                       dystar(l,k,m)=dystar(l,k,m)+dutc(m)
                        endif
                     endif
                  if((iobsrv(k).eq.0).or.(iy(l,k).lt.1).or.
     &             (iy(l,k).gt.nseg(l)))then
                     ystar(l,k)=xb(l,k)+u(l,k)
                     if(ifder.eq.1)then
                        do 32 m=1,nparm
  32                       dystar(l,k,m)=dxb(l,k,m)
                        endif
                     if(km1.ge.istart)then
                        do 15 m=1,ny
                           ystartt(l,k)=ystar(l,k)
  15                       ystar(l,k)=ystar(l,k)+(rhop(l,m,k)*
     &                      ystartt(m,km1))
                        if(ifder.eq.1)then
                           do 33 m=1,nparm
                              do 34 n=1,ny
  34                             dystar(l,k,m)=dystar(l,k,m)+
     &                            (rhop(l,n,k)*dystar(l,km1,m))+
     &                            (drhop(l,n,k,m)*ystartt(n,km1))
  33                          continue
                           endif 
                        endif
                     do 16 m=1,nfac
  16                    ystar(l,k)=ystar(l,k)+(alamb(l,m)*e(m))
                     if(ifder.eq.1)then
                        do 35 m=1,nparm
                           do 36 n=1,nfac
  36                          dystar(l,k,m)=dystar(l,k,m)+
     &                         (dalamb(l,n,m)*e(n))
  35                       continue
                        endif
                     endif
  11              continue
   4           continue
   3        continue
   2     continue
c
c       find maximum.
c
      fttmax=-1.d20
      do 7 i=1,nr
         do 8 j=1,2
            if(ftt(i,j).gt.fttmax)then
               ifttmax(1)=i
               ifttmax(2)=j
               fttmax=ftt(i,j)
               endif
   8        continue
   7     continue
      if(ifder.eq.1)then
         do 20 i=1,nparm
  20        dfttmax(i)=dftt(ifttmax(1),ifttmax(2),i)
         do 37 i=1,nparm
  37        df(i)=0.d0
         endif
      f=0.d0
      do 9 i=1,nr
         do 10 j=1,2
            fac=exp(ftt(i,j)-fttmax)
            if(ifder.eq.1)then
               do 19 k=1,nparm
  19              dfac(k)=fac*(dftt(i,j,k)-dfttmax(k))
               endif
            f=f+fac
            if(ifder.eq.1)then
               do 18 k=1,nparm
  18              df(k)=df(k)+dfac(k)
               endif
  10        continue
   9     continue
      if(ifder.eq.1)then
         do 17 i=1,nparm
  17        df(i)=dfttmax(i)+(df(i)/f)
         endif
      f=fttmax+log(f)-log(dble(float(nr*2)))
      return
      end
c
c
      subroutine linv3f (a,b,ijob,n,ia,d1,d2,wkarea,ier)                
c
c   imsl routine name   - linv3f                                        
c                                                                       
c-----------------------------------------------------------------------
c                                                                       
c   computer            - prim77/single                                 
c                                                                       
c   latest revision     - june 1, 1982                                  
c                                                                       
c   purpose             - in place inverse, equation solution, and/or   
c                           determinant evaluation - full storage mode  
c                                                                       
c   usage               - call linv3f (a,b,ijob,n,ia,d1,d2,wkarea,ier)  
c                                                                       
c   arguments    a      - input/output matrix of dimension n by n. see  
c                           parameter ijob.                             
c                b      - input/output vector of length n when ijob =   
c                           2 or 3. otherwise, b is not used.           
c                         on input, b contains the right hand side of   
c                           of the equation ax = b.                     
c                         on output, the solution x replaces b.         
c                ijob   - input option parameter. ijob = i implies:     
c                           i = 1, invert matrix a. a is replaced by    
c                             its inverse.                              
c                           i = 2, solve the equation ax = b. a is      
c                             replaced by the lu decomposition of a     
c                             rowwise permutation of a, where u is      
c                             upper triangular and l is lower           
c                             triangular with unit diagonal.            
c                             the unit diagonal of l is not stored.     
c                           i = 3, solve ax = b and invert matrix a.    
c                             a is replaced by its inverse.             
c                           i = 4, compute the determinant of a.        
c                             a is replaced by the lu decomposition     
c                             of a rowwise permutation of a.            
c                n      - order of a. (input)                           
c                ia     - row dimension of matrix a exactly as          
c                           specified in the dimension statement in the 
c                           calling program. (input)                    
c                d1     - input/output. if the d1 and d2 components     
c                d2         of determinant(a) = d1*2**d2 are desired,   
c                           input d1.ge.0. otherwise, input d1.lt.0.    
c                           d2 is never input.                          
c                wkarea - work area of length at least 2*n for ijob=1   
c                           or ijob=3.                                  
c                         work area of length at least n for ijob=2     
c                           or ijob=4.                                  
c                ier    - error parameter. (output)                     
c                         warning with fix                              
c                           ier = 65 indicates that ijob was less than  
c                             1 or greater than 4. ijob is assumed to   
c                             be 4.                                     
c                         terminal error                                
c                           ier = 130 indicates that matrix a is        
c                             algorithmically singular. (see the        
c                             chapter l prelude).                       
c                                                                       
c   precision/hardware  - single and double/h32                         
c                       - single/h36,h48,h60                            
c                                                                       
c   reqd. imsl routines - ludatn,luelmn
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
c-----------------------------------------------------------------------
c                                                                       
c                                                                       
      implicit real*8 (a-h,o-z)
c     implicit real (a-h,o-z)
      real*8             a(ia,n),b(n),wkarea(1000)  
c     real*4             a(ia,n),b(n),wkarea(400)                       
      data               zero/0.d0/,one/1.d0/                           
c                                  first executable statement           
c                                  lu decomposition of a                
      call ludatn (a,ia,n,a,ia,0,c1,c2,wkarea,wkarea,wa,ier)            
      if (d1 .lt. zero .and. ijob .ge. 1 .and. ijob .lt. 4) go to 5     
      d1 = c1                                                           
      d2 = c2                                                           
    5 if (ier .ge. 128) go to 60                                        
      if (ijob .le. 0 .or. ijob .gt. 4) go to 55                        
c                                  solve ax = b                         
      if (ijob .eq. 2 .or . ijob .eq. 3) call luelmn (a,ia,n,b,wkarea,b)
      if (ijob .ne. 1 .and. ijob .ne. 3) go to 9005                     
c                                  matrix inversion                     
      a(n,n) = one/a(n,n)                                               
      nm1 = n-1                                                         
      if (nm1 .lt. 1) go to 9005                                        
      do 40 ii=1,nm1                                                    
         l = n-ii                                                       
         m = l+1                                                        
         do 15 i=m,n                                                    
            sum = zero                                                  
            do 10 k=m,n                                                 
               sum = sum-a(i,k)*a(k,l)                                  
   10       continue                                                    
            wkarea(n+i) = sum                                           
   15    continue                                                       
         do 20 i=m,n                                                    
            a(i,l) = wkarea(n+i)                                        
   20    continue                                                       
         do 30 j=l,n                                                    
            sum = zero                                                  
            if (j .eq. l) sum = one                                     
            do 25 k=m,n                                                 
               sum = sum-a(l,k)*a(k,j)                                  
   25       continue                                                    
            wkarea(n+j) = sum/a(l,l)                                    
   30    continue                                                       
         do 35 j=l,n                                                    
            a(l,j) = wkarea(n+j)                                        
   35    continue                                                       
   40 continue                                                          
c                                  permute columns of a inverse         
      do 50 i=1,n                                                       
         j = n-i+1                                                      
         jp = wkarea(j)                                                 
         if (j .eq. jp) go to 50                                        
         do 45 k=1,n                                                    
            c = a(k,jp)                                                 
            a(k,jp) = a(k,j)                                            
            a(k,j) = c                                                  
   45    continue                                                       
   50 continue                                                          
      go to 9005                                                        
   55 continue                                                          
c                                  warning with fix - ijob was set      
c                                  incorrectly                          
      ier = 65                                                          
      go to 9000                                                        
c                                  terminal error - matrix a is         
c                                  algorithmically singular             
   60 ier = 130                                                         
 9000 continue                                                          
c     call uertst(ier,'linv3f')                                         
 9005 return                                                            
      end                                                               
c
c
      subroutine lmcovar
c
c       this subroutine computes the terms necessary to do each of the
c       lagrange multiplier tests.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2),nextras=
     & (nys*nfacs)+(nys*nys)+(nys*(nys-1)/2),ncoefts=ncoefs+nextras,
     & ncoeftss=ncoefts*(ncoefts+1)/2,ncoeft2s=2*ncoefts)
      implicit real*8(a-h,o-z)
      character*2 start,star(3)
      character*8 alabel
      real xextra,xv
      dimension amean(ncoefts),cov(ncoefts,ncoefts),covs(ncoeftss),
     & dft(ncoefs),eigval(ncoefts),eigvec(ncoefts,ncoefts),
     & parm(ncoefts),work1(ncoefts),work2(ncoeft2s)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/coefts/coeft(ncoefts),icoeft(ncoefts),ncoeft
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/labels/alabel(ncoefs)
      data eigtol/0.001d0/
      data star/'  ','* ','**'/
c
c       initialize.
c
      nextra=(ny*nfac)+(ny*ny)+(ny*(ny-1)/2)
      ncoeft=ncoef+nextra
      do 13 i=1,ncoeft
  13     amean(i)=0.d0
      do 6 i=1,ncoeft
         do 7 j=1,i
   7        cov(i,j)=0.d0
   6     continue
      call parset0t(parm,nparm)
c
c       loop over observations.
c
      anobs=0.d0
      do 2 i=1,nobs
         if(((i/100)*100).eq.i)write(6,111)i
 111     format(1x,'i = ',i5)
         if(incluv(i).eq.1)then
            anobs=anobs+1.d0
c
c       evaluate at the parameter estimate.
c
            call parset1t(parm)
            call parset212t(penalty)
            call obsmake2(i)
            call likobs12t(f0)
            f0=f0+penalty
c
c       loop over parameters and compute numerical derivatives.
c
            iparm=0
            do 3 j=1,ncoeft
               if(icoeft(j).eq.0)then
                  iparm=iparm+1
                  delta=coeft(j)/100.d0
                  if(delta.gt.1.d0)delta=1.d0
                  if(delta.lt.0.01d0)delta=0.01d0
                  coeft(j)=coeft(j)+delta
                  call parset212t(penalty)
                  call likobs12t(f1)
                  f1=f1+penalty
                  dft(iparm)=(f1-f0)/delta
                  coeft(j)=coeft(j)-delta              
                  endif
   3           continue
c
c       update mean and covariance matrix.
c
            do 4 j=1,iparm
               amean(j)=amean(j)+dft(j)
               do 5 k=1,j
   5              cov(j,k)=cov(j,k)+(dft(j)*dft(k))
   4           continue
            endif
   2     continue
c
c       update and output mean derivatives.
c
      write(6,100)
 100  format(1x,50('='))
      do 8 i=1,iparm
   8     amean(i)=amean(i)/anobs
      ij=0
      do 9 i=1,iparm
         do 10 j=1,i
            cov(i,j)=(cov(i,j)/anobs)-(amean(i)*amean(j))
            cov(j,i)=cov(i,j)
            ij=ij+1
  10        covs(ij)=cov(i,j)
   9     continue
      write(6,109)
 109  format(1x,'mean derivatives')
      write(6,101)
 101  format(6x,'variable',2x,'mean',8x,'std err',3x,'t-stat')
      do 11 i=1,iparm
         stderr=(cov(i,i)/anobs)**.5d0
         tstat=amean(i)/stderr
         atstat=abs(tstat)
         if(atstat.lt.1.53d0)start=star(1)
         if((atstat.ge.1.53d0).and.(atstat.le.1.96d0))start=star(2)
         if(atstat.gt.1.96d0)start=star(3)
  11     write(6,102)i,alabel(i),amean(i),start,stderr,tstat
 102  format(1x,i3,2x,a8,2x,f8.3,a2,2(2x,f8.3))
c
c       compute eigenvalues.
c
      iflag=0
      call eigrs(covs,iparm,1,eigval,eigvec,ncoefts,work1,ier)
      do 12 i=1,iparm
         if(eigval(i).lt.eigtol)then
            if(iflag.eq.0)then
               iflag=1
               write(6,103)
 103           format(1x,'problem eigenvalues and eigenvectors',/,1x,
     &          6x,'eigenvalue',7x,'eigenvector')
               endif
            write(6,104)i,eigval(i),(eigvec(j,i),j=1,iparm)
 104        format(1x,i3,2x,g15.8,5(2x,f7.4),50(/,21x,5(2x,f7.4)))
            endif
  12     continue
c
c       adjust covariance matrix for bad eigenvalues.
c
      if(iflag.eq.1)then
         do 17 i=1,iparm
            if(eigval(i).lt.eigtol)eigval(i)=eigtol
  17        continue
         do 18 i=1,iparm
            do 19 j=1,i
               cov(i,j)=0.d0
               do 20 k=1,iparm
  20              cov(i,j)=cov(i,j)+(eigvec(i,k)*eigval(k)*eigvec(j,k))
  19           cov(j,i)=cov(i,j)
  18        continue
         endif 
c
c       invert.
c
      d1=-1.d0
      call linv3f(cov,work1,1,iparm,ncoefts,d1,d2,work2,ier)
      if(ier.ne.0)then
         write(6,106)ier
 106     format(1x,'linv3f: ier = ',i3)
         stop
         endif
c
c       construct lm test statistic.
c
      iparm=0
      do 14 i=1,ncoef
         if(icoeft(i).eq.0)then
            iparm=iparm+1
            amean(iparm)=0.d0
            endif
  14     continue
      iparm=0
      chi=0.d0
      do 15 i=1,ncoeft
         if(icoeft(i).eq.0)then
            iparm=iparm+1
            jparm=0
            do 16 j=1,ncoeft
               if(icoeft(j).eq.0)then
                  jparm=jparm+1
                  chi=chi+(amean(iparm)*cov(iparm,jparm)*amean(jparm))
                  endif
  16           continue
            endif
  15     continue
      chi=chi*anobs
      ndf=iparm
      andf=dble(float(ndf))
      adjust=(chi-andf)/((2.d0*andf)**.5d0)
      write(6,112)chi,ndf,adjust
 112  format(1x,'chi sq stat = ',g15.8,/,1x,'degrees of freedom = ',
     & i3,/,1x,'adjusted test stat = ',g15.8)
      return
      end
c
c
      subroutine lmtests
c
c       this subroutine performs four lagrange multiplier tests.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+(nys*nsegm2s)+nys*
     & (nfacs+nys)+(nys*(nys-1)/2),nextras=(nys*nfacs)+(nys*nys)+(nys*
     & (nys-1)/2),ncoefts=ncoefs+nextras)
      implicit real*8(a-h,o-z)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/coefts/coeft(ncoefts),icoeft(ncoefts),ncoeft
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/lmparms/alambt(nys,nfacs),omegat(nys,nys),
     & rhoptt(nys,nys,ntimes),rhot(nys,nys),iswitch(ntimes),lmflag
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
c
c       construct icoeft.
c
      do 2 i=1,ncoef
   2     icoeft(i)=icoef(i)
      nextra=(ny*nfac)+(ny*ny)+(ny*(ny-1)/2)
      ij1=ncoef
      ij2=(ny*ntime*ifwvdum)+(ny*nx)
      do 3 i=1,nextra
         ij1=ij1+1
         ij2=ij2+1
         coeft(ij1)=coef(ij2)
   3     icoeft(ij1)=icoef(ij2)
c
c       perform lagrange multiplier test for covariance parameters 
c       varying by education.
c
      lmflag=1
      call lmcovar
c
c       perform lagrange multiplier test for covariance parameters
c       varying by age.
c
      lmflag=2
      call datextraget(2)
      call lmcovar
c
c       perform lagrange multiplier test for covariance parameters
c       varying between pre- and post-brexit.
c
      lmflag=3
      call lmcovar
c
c       perform lagrange multiplier test for covariance parameters 
c       varying between election cycles and not election cycles.
c
      lmflag=4
      call lmcovar
      return
      end
c
c
      subroutine ludatn (a,ia,n,xlu,ilu,idgt,d1,d2,apvt,equil,wa,ier)   
c
c   imsl routine name   - ludatn                                        
c                                                                       
c-----------------------------------------------------------------------
c                                                                       
c   computer            - prim77/single                                 
c                                                                       
c   latest revision     - june 1, 1982                                  
c                                                                       
c   purpose             - nucleus called only by imsl subroutine leqt2f 
c                                                                       
c   precision/hardware  - single and double/h32                         
c                       - single/h36,h48,h60                            
c                                                                       
c   reqd. imsl routines -       
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
c-----------------------------------------------------------------------
c                                                                       
      implicit real*8 (a-h,o-z)
c     implicit real*4 (a-h,o-z)
c                                                                       
      real*8          a(ia,n),xlu(ilu,ilu),apvt(ilu),equil(n)           
c     real*4          a(ia,n),xlu(ilu,ilu),apvt(ilu),equil(n)           
      data  zero,one,four,sixtn,sixth/0.d0,1.d0,4.d0,16.d0,.0625d0/     
c                                  first executable statement           
c                                  initialization                       
      ier = 0                                                           
      rn = n                                                            
      wrel = zero                                                       
      d1 = one                                                          
      d2 = zero                                                         
      biga = zero                                                       
      do 10 i=1,n                                                       
         big = zero                                                     
         do 5 j=1,n                                                     
            p = a(i,j)                                                  
            xlu(i,j) = p                                                
c           p = dabs(p)                                                 
            p = abs(p)                                                 
            if (p .gt. big) big = p                                     
    5    continue                                                       
         if (big .gt. biga) biga = big                                  
         if (big .eq. zero) go to 110                                   
         equil(i) = one/big                                             
   10 continue                                                          
      do 105 j=1,n                                                      
         jm1 = j-1                                                      
         if (jm1 .lt. 1) go to 40                                       
c                                  compute u(i,j), i=1,...,j-1          
         do 35 i=1,jm1                                                  
            sum = xlu(i,j)                                              
            im1 = i-1                                                   
            if (idgt .eq. 0) go to 25                                   
c                                  with accuracy test                   
c           ai = dabs(sum)                                              
            ai = abs(sum)                                              
            wi = zero                                                   
            if (im1 .lt. 1) go to 20                                    
            do 15 k=1,im1                                               
               t = xlu(i,k)*xlu(k,j)                                    
               sum = sum-t                                              
c              wi = wi+dabs(t)                                          
               wi = wi+abs(t)                                          
   15       continue                                                    
            xlu(i,j) = sum                                              
c  20       wi = wi+dabs(sum)                                           
   20       wi = wi+abs(sum)                                           
            if (ai .eq. zero) ai = biga                                 
            test = wi/ai                                                
            if (test .gt. wrel) wrel = test                             
            go to 35                                                    
c                                  without accuracy                     
   25       if (im1 .lt. 1) go to 35                                    
            do 30 k=1,im1                                               
               sum = sum-xlu(i,k)*xlu(k,j)                              
   30       continue                                                    
            xlu(i,j) = sum                                              
   35    continue                                                       
   40    p = zero                                                       
c                                  compute u(j,j) and l(i,j), i=j+1,...,
         do 70 i=j,n                                                    
            sum = xlu(i,j)                                              
            if (idgt .eq. 0) go to 55                                   
c                                  with accuracy test                   
c           ai = dabs(sum)                                              
            ai = abs(sum)                                              
            wi = zero                                                   
            if (jm1 .lt. 1) go to 50                                    
            do 45 k=1,jm1                                               
               t = xlu(i,k)*xlu(k,j)                                    
               sum = sum-t                                              
c              wi = wi+dabs(t)                                          
               wi = wi+abs(t)                                          
   45       continue                                                    
            xlu(i,j) = sum                                              
c  50       wi = wi+dabs(sum)                                           
   50       wi = wi+abs(sum)                                           
            if (ai .eq. zero) ai = biga                                 
            test = wi/ai                                                
            if (test .gt. wrel) wrel = test                             
            go to 65                                                    
c                                  without accuracy test                
   55       if (jm1 .lt. 1) go to 65                                    
            do 60 k=1,jm1                                               
               sum = sum-xlu(i,k)*xlu(k,j)                              
   60       continue                                                    
            xlu(i,j) = sum                                              
c  65       q = equil(i)*dabs(sum)                                      
   65       q = equil(i)*abs(sum)                                      
            if (p .ge. q) go to 70                                      
            p = q                                                       
            imax = i                                                    
   70    continue                                                       
c                                  test for algorithmic singularity     
         q = rn+p                                                       
   71    if (q .eq. rn) go to 110                                       
         if (j .eq. imax) go to 80                                      
c                                  interchange rows j and imax          
         d1 = -d1                                                       
         do 75 k=1,n                                                    
            p = xlu(imax,k)                                             
            xlu(imax,k) = xlu(j,k)                                      
            xlu(j,k) = p                                                
   75    continue                                                       
         equil(imax) = equil(j)                                         
   80    apvt(j) = imax                                                 
         d1 = d1*xlu(j,j)                                               
c  85    if (dabs(d1) .le. one) go to 90                                
   85    if (abs(d1) .le. one) go to 90                                
         d1 = d1*sixth                                                  
         d2 = d2+four                                                   
         go to 85                                                       
c  90    if (dabs(d1) .ge. sixth) go to 95                              
   90    if (abs(d1) .ge. sixth) go to 95                              
         d1 = d1*sixtn                                                  
         d2 = d2-four                                                   
         go to 90                                                       
   95    continue                                                       
         jp1 = j+1                                                      
         if (jp1 .gt. n) go to 105                                      
c                                  divide by pivot element u(j,j)       
         p = xlu(j,j)                                                   
         do 100 i=jp1,n                                                 
            xlu(i,j) = xlu(i,j)/p                                       
  100    continue                                                       
  105 continue                                                          
c                                  perform accuracy test                
      if (idgt .eq. 0) go to 9005                                       
      p = 3*n+3                                                         
      wa = p*wrel                                                       
      q = wa+10.0**(-idgt)                                              
  106 if (q .ne. wa) go to 9005                                         
      ier = 34                                                          
      go to 9000                                                        
c                                  algorithmic singularity              
  110 ier = 129                                                         
      d1 = zero                                                         
      d2 = zero                                                         
 9000 continue                                                          
c                                  print error                          
c     call uertst(ier,'ludatn')                                         
 9005 return                                                            
      end                                                               
c
c
      subroutine luelmn (a,ia,n,b,apvt,x)                               
c
c   imsl routine name   - luelmn                                        
c                                                                       
c-----------------------------------------------------------------------
c                                                                       
c   computer            - prim77/single                                 
c                                                                       
c   latest revision     - june 1, 1982                                  
c                                                                       
c   purpose             - nucleus called only by imsl subroutine leqt2f 
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
c-----------------------------------------------------------------------
c                                                                       
      implicit real*8 (a-h,o-z)
c     implicit real*4 (a-h,o-z)
c                                                                       
      real*8          a(ia,n),b(ia),apvt(ia),x(ia)                      
c     real*4          a(ia,n),b(ia),apvt(ia),x(ia)                      
c                                  first executable statement           
c                                  solve ly = b for y                   
      do 5 i=1,n                                                        
    5 x(i) = b(i)                                                       
      iw = 0                                                            
      do 20 i=1,n                                                       
         ip = apvt(i)                                                   
         sum = x(ip)                                                    
         x(ip) = x(i)                                                   
         if (iw .eq. 0) go to 15                                        
         im1 = i-1                                                      
         do 10 j=iw,im1                                                 
            sum = sum-a(i,j)*x(j)                                       
   10    continue                                                       
         go to 20                                                       
   15    if (sum .ne. 0.) iw = i                                        
   20 x(i) = sum                                                        
c                                  solve ux = y for x                   
      do 30 ib=1,n                                                      
         i = n+1-ib                                                     
         ip1 = i+1                                                      
         sum = x(i)                                                     
         if (ip1 .gt. n) go to 30                                       
         do 25 j=ip1,n                                                  
            sum = sum-a(i,j)*x(j)                                       
   25   continue                                                        
   30 x(i) = sum/a(i,i)                                                 
      return                                                            
      end                                                               
c
c
c   imsl routine name   - mdnris
c 
c-----------------------------------------------------------------------
c 
c   computer            - cdcft5/single 
c 
c   latest revision     - november 1, 1984
c 
c   purpose             - inverse standard normal (gaussian)
c                           probability distribution function 
c 
c   usage               - call mdnris (p,y,ier) 
c 
c   arguments    p      - input value in the exclusive range (0.0,1.0)
c                y      - output value of the inverse normal (0,1)
c                           probability distribution function 
c                ier    - error parameter (output)
c                         terminal error
c                           ier = 129 indicates p lies outside the legal
c                             range. plus or minus machine infinity is
c                             given as the result (sign is the sign of
c                             the function value of the nearest legal 
c                             argument).
c 
c   precision/hardware  - single/all
c 
c   reqd. imsl routines - merfi 
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
c-----------------------------------------------------------------------
c 
      subroutine mdnris (p,y,ier) 
c                                  specifications for arguments 
      real               p,y
      integer            ier
c                                  specifications for local variables 
      real               eps,g0,g1,g2,g3,h0,h1,h2,a,w,wi,sn,sd
      real               sigma,sqrt2,x,xinf 
c     data               xinf/.12650140831e+323/
      data               sqrt2/1.4142135623731/ 
      data               eps/.710543e-14/ 
      data               g0/.18511591e-3/,g1/-.20281520e-2/ 
      data               g2/-.14983844/,g3/.10786386e-1/
      data               h0/.99529751e-1/,h1/.52117329/ 
      data               h2/-.68883009e-1/
c                                  first executable statement 
      xinf=huge(sigma)
      ier = 0 
      if (p .gt. 0.0 .and. p .lt. 1.0) go to 5
      ier = 129 
      sigma = sign(1.0,p) 
      y = sigma * xinf
      go to 9000
    5 if(p.le.eps) go to 10 
      x = 1.0 -(p + p)
      call merfi (x,y,ier)
      y = -sqrt2 * y
      go to 9005
c                                  p too small, compute y directly
   10 a = p+p 
      w = (-log(a+(a-a*a)))**.5
c                                  use a rational function in 1./w
      wi = 1./w 
      sn = ((g3*wi+g2)*wi+g1)*wi
      sd = ((wi+h2)*wi+h1)*wi+h0
      y = w + w*(g0+sn/sd)
      y = -y*sqrt2
      go to 9005
 9000 continue
 9005 return
      end 
c
c
      subroutine merfi (p,y,ier)
c
c   imsl routine name   - merfi 
c 
c---------------------------------------------------------------------- 
c 
c   computer            - cdcft5/single 
c 
c   latest revision     - january 1, 1978 
c 
c   purpose             - inverse error function
c 
c   usage               - call merfi (p,y,ier)
c 
c   arguments    p      - input value in the exclusive range (-1.0,1.0) 
c                y      - output value of the inverse error function
c                ier    - error parameter (output)
c                         terminal error
c                           ier = 129 indicates p lies outside the lega 
c                             range. plus or minus machine infinity is
c                             given as the result (sign is the sign of
c                             the function value of the nearest legal 
c                             argument).
c 
c   precision/hardware  - single/all
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
      implicit real(a-h,o-z)
c                                  specifications for arguments 
      real               p,y
      integer            ier
c                                  specifications for local variables 
      real               a(65),h1,h2,h3,h4,xinf,x,sigma,z,w,x3,x4,x5, 
     *                   x6,b 
      integer            n,ipp,l,lb2
      data               a(1),a(2),a(3),a(4),a(5),a(6),a(7),a(8),a(9),
     *                   a(10),a(11),a(12),a(13),a(14),a(15),a(16), 
     *                   a(17),a(18),a(19),a(20),a(21),a(22),a(23)
     *                   /.99288537662,.12046751614,
     *                   .16078199342e-01,.26867044372e-02, 
     *                   .49963473024e-03,.98898218599e-04, 
     *                   .20391812764e-04,.4327271618e-05,
     *                   .938081413e-06,.206734721e-06, 
     *                   .46159699e-07,.10416680e-07, 
     *                   .2371501e-08,.543928e-09,
     *                   .125549e-09,.29138e-10,
     *                   .6795e-11,.1591e-11, 
     *                   .374e-12,.88e-13,
     *                   .21e-13,.5e-14,
     *                   .1e-14/
      data               a(24),a(25),a(26),a(27),a(28),a(29),a(30), 
     *                   a(31),a(32),a(33),a(34),a(35),a(36),a(37), 
     *                   a(38),a(39),a(40)
     *                   /.91215880342e00,-.16266281868e-01,
     *                   .43355647295e-03,.21443857007e-03, 
     *                   .2625751076e-05,-.302109105e-05, 
     *                   -.12406062e-07,.62406609e-07,
     *                   -.540125e-09,-.142328e-08, 
     *                   .34384e-10,.33585e-10, 
     *                   -.1458e-11,-.81e-12, 
     *                   .53e-13,.2e-13,
     *                   -.2e-14/ 
      data               a(41),a(42),a(43),a(44),a(45),a(46),a(47), 
     *                   a(48),a(49),a(50),a(51),a(52),a(53),a(54), 
     *                   a(55),a(56),a(57),a(58),a(59),a(60),a(61), 
     *                   a(62),a(63),a(64),a(65)
     *                   /.95667970902,-.023107004309,
     *                   -.43742360975e-02,-.57650342265e-03, 
     *                   -.10961022307e-04,.25108547025e-04,
     *                   .10562336068e-04,.275441233e-05, 
     *                   .432484498e-06,-.20530337e-07, 
     *                   -.43891537e-07,-.1768401e-07,
     *                   -.3991289e-08,-.186932e-09,
     *                   .272923e-09,.132817e-09, 
     *                   .31834e-10,.1670e-11,
     *                   -.2036e-11,-.965e-12,
     *                   -.22e-12,-.1e-13,
     *                   .13e-13,.6e-14,
     *                   .1e-14/
      data               h1,h2,h3,h4/-1.5488130424, 
     *                   2.5654901231,-.55945763133,
     *                   2.2879157163/
c     data               xinf/.12650140831e+323/
c                                  first executable statement 
      xinf=huge(sigma)
      x = p 
      ier = 0 
      sigma = sign(1.,x)
      if (.not.(x.gt.-1..and.x.lt.1.)) go to 35 
      z = abs(x)
      if(z.gt..8) go to 20
      w = z*z/.32-1.
      n = 22
      ipp = 1 
      l = 1 
    5 lb2 = 1 
      x3 = 1. 
      x4 = w
      x6 = a(ipp) 
   10 x6 = x6 + a(ipp+lb2) * x4 
      x5 = x4 * w * 2.-x3 
      x3 = x4 
      x4 = x5 
      lb2 = lb2 + 1 
      if (lb2 .le. n) go to 10
      go to (15,30),l 
   15 y = z * x6 * sigma
      go to 9005
   20 b = (-log(1.-z*z))**.5 
      if (z .gt. .9975) go to 25
      w = h1*b+h2 
      ipp = 24
      l = 2 
      n = 16
      go to 5 
   25 w = h3 * b + h4 
      ipp = 41
      n = 24
      l = 2 
      go to 5 
   30 y = b * x6 * sigma
      go to 9005
   35 y = sigma*xinf
      ier = 129 
 9000 continue
 9005 return
      end 
c
c
      subroutine normadj1(iyk,iyl,ystark,uk,utc)
c
c       this subroutine adjusts the standard normal to be consistent
c       with the observed choice.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,nsegp1s=
     & nsegs+1,ntimes=34,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+nys*(nfacs+nys)+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      real etcr,utcr
      call anormc(uk,fuk)
      arg1=thresh(iyl,iyk)-ystark
      call anormc(arg1,pfac1)
      arg2=thresh(iyl,iyk+1)-ystark
      call anormc(arg2,pfac2)
      fac1=pfac2-pfac1
      utc=(fac1*fuk)+pfac1
      utcr=sngl(utc)
      call utrunc(utcr)
      call mdnris(utcr,etcr,ier)
      utc=dble(etcr)
      return
      end
c
c
      subroutine normadj2(iyk,iyl,ystark,dystark,uk,utc,dutc,nparm,
     & ifder)
c
c       this subroutine adjusts the standard normal to be consistent
c       with the observed choice.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,nsegp1s=
     & nsegs+1,ntimes=34,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+nys*(nfacs+nys)+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension darg1(ncoefs),darg2(ncoefs),dfac1(ncoefs),
     & dpfac1(ncoefs),dpfac2(ncoefs),dutc(ncoefs),dystark(ncoefs)
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      real etcr,utcr
      call anormc(uk,fuk)
      arg1=thresh(iyl,iyk)-ystark
      call anormc(arg1,pfac1)
      arg2=thresh(iyl,iyk+1)-ystark
      call anormc(arg2,pfac2)
      fac1=pfac2-pfac1
      utc=(fac1*fuk)+pfac1
      utcr=sngl(utc)
      call utrunc(utcr)
      call mdnris(utcr,etcr,ier)
      utc=dble(etcr)
      if(ifder.eq.1)then
         call anormd(arg1,pfac1t)
         call anormd(arg2,pfac2t)
         call anormd(utc,etcrt)
         do 2 i=1,nparm
            darg1(i)=dthresh(iyl,iyk,i)-dystark(i)
            darg2(i)=dthresh(iyl,iyk+1,i)-dystark(i)
            dpfac1(i)=pfac1t*darg1(i)
            dpfac2(i)=pfac2t*darg2(i)
            dfac1(i)=dpfac2(i)-dpfac1(i)
            dutc(i)=(fuk*dfac1(i))+dpfac1(i)
   2        dutc(i)=dutc(i)/etcrt
         endif
      return
      end
c
c
      subroutine obsmake1(iobs)
c
c       this subroutine pulls out the data relevant to iobs for the 
c       simple model.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nrs=10,nsegs=8,nsegm2s=
     & nsegs-2,ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+
     & (nys*nbetas)+(nys*nsegm2s)+nys*(nfacs+nys)+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      real erv,urv,xextra,xv
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/errors/erv(nobss,nfacs,nrs),urv(nobss,nrs,nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
!$omp threadprivate(/datat/)
      istart=istartv(iobs)
      iend=iendv(iobs)
      do 7 i=istart,iend 
         do 8 j=1,nx
   8        x(i,j)=dble(xv(iobs,i,j))
   7     continue
      do 5 i=istart,iend 
   5     iy(1,i)=iyv(iobs,1,i)
      do 6 i=istart,iend 
   6     iobsrv(i)=iobsrvv(iobs,i)
      do 2 i=1,nr
   2     er(i)=dble(erv(iobs,1,i))
      do 3 i=1,nr
         do 4 j=istart,iend 
   4        ur(i,1,j)=dble(urv(iobs,i,1,j))
   3     continue
      return
      end
c
c
      subroutine obsmake2(iobs)
c
c       this subroutine pulls out the data relevant to iobs for the 
c       complex model.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nrs=10,nsegs=8,nsegm2s=
     & nsegs-2,ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+
     & (nys*nbetas)+(nys*nsegm2s)+nys*(nfacs+nys)+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      real erv,urv,xextra,xv
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/errors/erv(nobss,nfacs,nrs),urv(nobss,nrs,nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
!$omp threadprivate(/datat/)
      istart=istartv(iobs)
      iend=iendv(iobs)
      do 7 i=istart,iend 
         do 8 j=1,nx
   8        x(i,j)=dble(xv(iobs,i,j))
   7     continue
      do 5 i=istart,iend 
         do 9 j=1,ny
   9        iy(j,i)=iyv(iobs,j,i)
   5     continue
      do 6 i=istart,iend 
   6     iobsrv(i)=iobsrvv(iobs,i)
      do 2 i=1,nr
         do 10 j=1,nfac
  10        er2(i,j)=dble(erv(iobs,j,i))
   2     continue
      do 3 i=1,nr
         do 4 j=istart,iend 
            do 11 k=1,ny
  11           ur(i,k,j)=dble(urv(iobs,i,k,j))
   4        continue
   3     continue
      if(method.eq.5)then
         do 12 i=1,ntime
  12        xextrat(i)=dble(xextra(iobs,i))
         endif
      return
      end
c
c
      subroutine olsresids
c
c       this subroutine performs ols on each dependent variable and then
c       looks at the covariance characteristics of the residuals.
c
      parameter(nobss=96529,ntimes=34,nxs=20,nx2s=nxs*2,nys=6)
      implicit real*8(a-h,o-z)
      real xextra,xv
      dimension b(nxs),corr(nys,nys),corr1(nys,nys),cov(nys,nys),
     & cov1(nys,nys),resid(nys,nobss,ntimes),work1(nxs),work2(nx2s),
     & xx(nxs,nxs),xy(nxs)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      data tolres/1.d-3/
c
c       do ols on each dependent variable and compute residuals.
c
      do 2 i=1,ny
         do 3 j=1,nx
            xy(j)=0.d0
            anobs=0.d0
            do 4 k=1,nobs
               if(incluv(k).eq.1)then
                  istart=istartv(k)
                  iend=iendv(k)
                  do 21 l=istart,iend
                     if((iyv(k,i,l).gt.0).and.(iyv(k,i,l).lt.10).and.
     &                (xv(k,l,1).eq.1.))then
                        anobs=anobs+1.d0
                        xy(j)=xy(j)+(dble(xv(k,l,j))*
     &                   dble(float(iyv(k,i,l)))) 
                        endif
  21                 continue
                  endif
   4           continue
   3        xy(j)=xy(j)/anobs
         do 5 j=1,nx
            do 6 k=1,j
               xx(j,k)=0.d0
               do 7 l=1,nobs
                  if(incluv(l).eq.1)then
                     istart=istartv(l)
                     iend=iendv(l)
                     do 22 m=istart,iend
                        if((iyv(l,i,m).gt.0).and.(iyv(l,i,m).lt.10).and.
     &                   (xv(l,m,1).eq.1.))xx(j,k)=xx(j,k)+
     &                   (dble(xv(l,m,j)*xv(l,m,k)))
  22                    continue
                     endif
   7              continue
               xx(j,k)=xx(j,k)/anobs
   6           xx(k,j)=xx(j,k)
   5        continue
         d1=-1.d0
         call linv3f(xx,work1,1,nx,nxs,d1,d2,work2,ier)
         do 8 j=1,nx
            b(j)=0.d0
            do 9 k=1,nx
   9           b(j)=b(j)+(xx(j,k)*xy(k))
   8        continue
         residtot=0.d0
         do 10 j=1,nobs
            if(incluv(j).eq.1)then
               istart=istartv(j)
               iend=iendv(j)
               do 23 k=istart,iend
                  if((iyv(j,i,k).gt.0).and.(iyv(j,i,k).lt.10).and.
     &             (xv(j,k,1).eq.1.))then
                     resid(i,j,k)=dble(float(iyv(j,i,k)))
                     do 11 l=1,nx
  11                    resid(i,j,k)=resid(i,j,k)-(b(l)*dble(xv(j,k,l)))
                     residtot=residtot+resid(i,j,k)
                     endif
  23              continue
               endif
  10        continue
         residtot=residtot/anobs
         write(6,100)i,residtot
 100     format(1x,'check: resid total: ',i1,2x,g15.8)
         if(abs(residtot).gt.tolres)stop
   2     continue
c
c       compute contemporaneous covariance matrix of residuals.
c
      write(6,101)
 101  format(1x,50('='))
      write(6,102)
 102  format(1x,'contemporaneous correlation matrix')
      write(6,103)(i,i=1,ny)
 103  format(1x,'variable',2x,6(i1,9x))
      do 12 i=1,ny
         do 13 j=1,i
            cov(i,j)=0.d0
            anobs=0.d0
            do 14 k=1,nobs
               if(incluv(k).eq.1)then
                  istart=istartv(k)
                  iend=iendv(k)
                  do 24 l=istart,iend
                     if((iyv(k,i,l).gt.0).and.(iyv(k,i,l).lt.10).and.
     &                (xv(k,l,1).eq.1.))then
                        anobs=anobs+1.d0
                        cov(i,j)=cov(i,j)+(resid(i,k,l)*resid(j,k,l))
                        endif
  24                 continue
                  endif
  14           continue
  13        cov(i,j)=cov(i,j)/anobs
  12     continue
      do 15 i=1,ny
         do 16 j=1,i
  16        corr(i,j)=cov(i,j)/((cov(i,i)*cov(j,j))**.5d0)
  15     write(6,104)i,(corr(i,j),j=1,i)
 104  format(1x,i1,7x,6(2x,f6.3))
c
c       compute correlation of each variable with lags of each
c        variable.
c
      write(6,101)
      write(6,105)
 105  format(1x,'correlation of each variable with each lagged',
     & ' variable')
      write(6,106)(i,i=1,ny)
 106  format(1x,'variable',2x,6(i1,7x))
      do 17 i=1,ny
         do 18 j=1,ny
            anobsc=0.d0
            cov1(i,j)=0.d0
            do 19 k=2,nobs
               if(incluv(k).eq.1)then
                  istartp1=istartv(k)+1
                  iend=iendv(k)
                  do 25 l=istartp1,iend
                     lm1=l-1
                     if((iyv(k,i,l).gt.0).and.(iyv(k,i,l).lt.10).and.
     &                (xv(k,l,1).eq.1.).and.(iyv(k,j,lm1).gt.0).and.
     &                (iyv(k,j,lm1).lt.10).and.(xv(k,lm1,1).eq.1.))then
                        anobsc=anobsc+1.d0
                        cov1(i,j)=cov1(i,j)+(resid(i,k,l)*
     &                   resid(j,k,lm1))
                        endif
  25                 continue
                  endif
  19           continue
  18        continue
         do 20 j=1,ny
  20        corr1(i,j)=(cov1(i,j)/anobsc)/((cov(i,i)*cov(j,j))**.5d0)
  17     write(6,107)i,(corr1(i,j),j=1,ny)
 107  format(1x,i1,7x,6(2x,f6.3))
c
c       compute idiosyncratic and person-specific standard deviations.
c
      write(6,101)
      write(6,109)
 109  format(1x,'person-specific and idiosyncratic standard',
     & ' deviations for each set of residuals')
      do 27 i=1,ny
         sig2i=0.d0
         sig2p=0.d0
         anobsp=0.d0
         anobsi=0.d0
         do 28 j=1,nobs
            if(incluv(j).eq.1)then
               istart=istartv(j)
               iend=iendv(j)
               amean=0.d0
               anobsa=0.d0
               do 29 k=istart,iend
                  if((iyv(j,i,k).gt.0).and.(iyv(j,i,k).lt.10).and.
     &             (xv(j,k,1).eq.1.))then
                     sig2i=sig2i+(resid(i,j,k)*resid(i,j,k))
                     anobsi=anobsi+1.d0
                     amean=amean+resid(i,j,k)
                     anobsa=anobsa+1.d0
                     endif
  29              continue
               if(anobsa.gt.0.d0)then
                  amean=amean/anobsa
                  anobsp=anobsp+1.d0
                  sig2p=sig2p+(amean*amean)
                  endif
               endif
  28        continue
         sig2p=(sig2p/anobsp)**.5d0
         sig2i=(sig2i/anobsi)**.5d0
  27     write(6,108)i,sig2p,sig2i
 108  format(1x,'dep var: ',i1,2x,'pers-spec std: ',g15.8,2x,
     & 'idiosyn std: ',g15.8) 
      return
      end
c
c
      subroutine omegamake1(chol,omega,penalty)
c
c       this subroutine constructs omega from chol.  the diagonal 
c       elements of omega must be one.  if the off-diagonal elements 
c       are too large for that to happen, then omega is adjusted, and a
c       penalty is imposed proportional to the square of the 
c       adjustment.
c
      parameter(ntimes=34,nys=6)
      implicit real*8(a-h,o-z)
      dimension chol(nys,nys),omega(nys,nys)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      data penrat/1.d2/
      data tol/.1d-1/
      penalty=0.d0
      omega(1,1)=1.d0
      do 2 i=2,ny
         im1=i-1
         do 3 j=1,im1
            omega(i,j)=0.d0
            do 4 k=1,j
   4           omega(i,j)=omega(i,j)+(chol(i,k)*chol(j,k))
   3        omega(j,i)=omega(i,j)
         omegsum=1.d0
         do 5 j=1,im1
   5        omegsum=omegsum-(omega(i,j)*omega(i,j))
         if(omegsum.ge.tol)omega(i,i)=1.d0
         if(omegsum.lt.tol)then
            omega(i,i)=1.d0
            dev=tol-omegsum
            penalty=penalty+(penrat*dev*dev)
            ratio=(1.d0-tol)/(1.d0-omegsum)
            do 6 j=1,im1
   6           omega(i,j)=omega(i,j)*ratio
            endif
   2     continue
      return
      end
c
c
      subroutine omegamake2(chol,dchol,omega,domega,penalty,dpenalty,
     & nparm,ifder)
c
c       this subroutine constructs omega from chol and its partial
c       derivatives.  the diagonal elements of omega must be one.  if 
c       the off-diagonal elements   are too large for that to happen, 
c       then omega is adjusted, and apenalty is imposed proportional to
c       the square of the adjustment.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+(nys*nsegm2s)+
     & (nys*(nfacs+nys)+(nys*(nys-1)/2)))
      implicit real*8(a-h,o-z)
      dimension chol(nys,nys),dchol(nys,nys,ncoefs),ddev(nys,ncoefs),
     & dev(nys),domega(nys,nys,ncoefs),domegsum(ncoefs),
     & dpenalty(ncoefs),dratio(ncoefs),omega(nys,nys),omegsum(nys)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      data penrat/1.d2/
      data tol/.1d-1/
      penalty=0.d0
      omega(1,1)=1.d0
      do 2 i=2,ny
         im1=i-1
         do 3 j=1,im1
            omega(i,j)=0.d0
            do 4 k=1,j
   4           omega(i,j)=omega(i,j)+(chol(i,k)*chol(j,k))
   3        omega(j,i)=omega(i,j)
         omegsum(i)=1.d0
         do 5 j=1,im1
   5        omegsum(i)=omegsum(i)-(omega(i,j)*omega(i,j))
         if(omegsum(i).ge.tol)omega(i,i)=1.d0
         if(omegsum(i).lt.tol)then
            omega(i,i)=1.d0
            dev(i)=tol-omegsum(i)
            penalty=penalty+(penrat*dev(i)*dev(i))
            ratio=(1.d0-tol)/(1.d0-omegsum(i))
            do 6 j=1,im1
   6           omega(i,j)=omega(i,j)*ratio
            endif
   2     continue
      if(ifder.eq.1)then
         do 7 i=1,nparm
            dpenalty(i)=0.d0
            domega(1,1,i)=0.d0
            do 8 j=2,ny
               jm1=j-1
               do 9 k=1,jm1
                  domega(j,k,i)=0.d0
                  do 10 l=1,k
  10                 domega(j,k,i)=domega(j,k,i)+(chol(j,l)*
     &                dchol(k,l,i))+(chol(k,l)*dchol(j,l,i))
   9              domega(k,j,i)=domega(j,k,i)
               domegsum(i)=0.d0
               do 11 k=1,jm1
  11              domegsum(i)=domegsum(i)-(2.d0*omega(j,k)*
     &             domega(j,k,i))
               if(omegsum(j).ge.tol)domega(j,j,i)=0.d0
               if(omegsum(j).lt.tol)then
                  domega(j,j,i)=0.d0
                  ddev(j,i)=-domega(j,j,i)
                  dpenalty(i)=dpenalty(i)+(2.d0*penrat*dev(j)*
     &             ddev(j,i))
                  dratio(i)=-1.d0*domegsum(i)*(1.d0-tol)/
     &             ((1.d0-omegsum(j))*(1.d0-omegsum(j)))
                  do 12 k=1,im1
  12                 domega(j,k,i)=(omega(j,k)*dratio(i))+(ratio*
     &                domega(j,k,i))
                  endif
   8           continue
   7        continue
         endif
      return
      end
c
c
      subroutine parget
c
c       this subroutine gets parameters.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+(nys*nsegm2s)+nys*
     & (nfacs+nys)+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      character*8 ajunk,alabel
      dimension nsegm2(nys)
      common/bnel/stp1,var,konvge,nrst
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/labels/alabel(ncoefs)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      data niterg/3/
      data niterb/6/
      open(unit=8,file='alecia.p')
      read(8,*)method,model,idata,nobs,ny,nx,nr,ntime,nfac,irhobessel,
     & ifcond1,niter,stp1,nrst,ifsattr,ifwvdum
      if(model.eq.1)then
         if(ny.gt.1)then
            write(6,103)ny
 103        format(1x,'model = 1 and ny = ',i2)
            stop
            endif
         if(nfac.gt.1)then
            write(6,104)nfac
 104        format(1x,'model = 1 and nfac = ',i2)
            stop
            endif
         endif
      if((model.eq.2).and.(ifcond1.eq.1))then
         write(6,102)
 102     format(1x,'model = 2 and ifcond1 = 1')
         stop
         endif
      if(idata.eq.1)then
         if(nobs.ne.96529)then
            write(6,105)idata,nobs
 105        format(1x,'idata = ',i1,2x,'nobs = ',i5)
            stop
            endif
         if(ny.ne.5)then
            write(6,106)idata,ny
 106       format(1x,'idata = ',i1,2x,'ny = ',i2)
            stop
            endif
         if(nx.ne.20)then
            write(6,107)idata,nx
 107        format(1x,'idata = ',i1,2x,'nx = ',i2)
            stop
            endif
         if(ntime.ne.19)then
            write(6,108)idata,ntime
 108        format(1x,'idata = ',i1,2x,'ntime = ',i3)
            stop
            endif
         endif
      if(idata.eq.2)then
         if(nobs.ne.72390)then
            write(6,105)idata,nobs
            stop
            endif
         if(ny.ne.6)then
            write(6,106)idata,ny
            stop
            endif
         if(nx.ne.7)then
            write(6,107)idata,nx
            stop
            endif
         if(ntime.ne.34)then
            write(6,108)idata,ntime
            stop
            endif
         endif
      if(idata.eq.3)then
         if(nobs.ne.21648)then
            write(6,105)idata,nobs
            stop
            endif
         if(ny.ne.6)then
            write(6,106)idata,ny
            stop
            endif
         if(nx.ne.10)then
            write(6,107)idata,nx
            stop
            endif
         if(ntime.ne.8)then
            write(6,108)idata,ntime
            stop
            endif
         endif
      read(8,*)(nseg(i),i=1,ny)
      do 5 i=1,ny
   5     nsegm2(i)=nseg(i)-2
      nbeta=nx
      if(model.eq.1)ncoef=ntime+nbeta+nsegm2(1)+2
      if(model.eq.2)then
         ncoef=(ny*ntime)+(ny*nbeta)+(ny*(nfac+ny))+(ny*(ny-1)/2)
         do 4 i=1,ny
   4        ncoef=ncoef+nsegm2(i)
         endif
      read(8,*)(coef(i),i=1,ncoef)
      read(8,*)(icoef(i),i=1,ncoef)
      if((model.eq.1).and.(ifsattr.eq.0))then
         do 3 i=17,19
            if(icoef(i).eq.0)then
               write(6,101)i
 101           format(1x,'ifsattr = 1 and icoef(',i2,') = 0')
               stop
               endif
   3        continue
         endif
      read(8,100)ajunk
 100  format(a8)
      iparm=0
      do 2 i=1,ncoef
         read(8,100)ajunk
         if(icoef(i).eq.0)then
            iparm=iparm+1
            alabel(iparm)=ajunk
            endif
   2     continue
      close(8)
      return
      end
c
c
      subroutine parset0(parm,nparm)
c
c       this subroutine moves coef to parm.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ncoefs=(ntimes*nys)+(nbetas*nys)+(nys*nsegm2s)+
     & (nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension parm(ncoefs)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      nparm=0
      do 2 i=1,ncoef
         if(icoef(i).eq.0)then
            nparm=nparm+1
            parm(nparm)=coef(i)
            endif
   2     continue
      return
      end
c
c
      subroutine parset0t(parm,nparm)
c
c       this subroutine moves coef to parm for lagrange multiplier
c       tests.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ncoefs=(ntimes*nys)+(nbetas*nys)+(nys*nsegm2s)+
     & (nys*(nfacs+nys))+(nys*(nys-1)/2),nextras=(nys*nfacs)+(nys*nys)+
     & (nys*(nys-1)/2),ncoefts=ncoefs+nextras)
      implicit real*8(a-h,o-z)
      dimension parm(ncoefts)
      common/coefts/coeft(ncoefts),icoeft(ncoefts),ncoeft
      nparm=0
      do 2 i=1,ncoeft
         if(icoeft(i).eq.0)then
            nparm=nparm+1
            parm(nparm)=coeft(i)
            endif
   2     continue
      return
      end
c
c
      subroutine parset1(parm)
c
c       this subroutine moves parm into coef.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ncoefs=(ntimes*nys)+(nbetas*nys)+(nys*nsegm2s)+
     & (nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension parm(ncoefs)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      iparm=0
      do 2 i=1,ncoef
         if(icoef(i).eq.0)then
            iparm=iparm+1
            coef(i)=parm(iparm)
            endif
   2     continue
      return
      end
c
c
      subroutine parset1t(parm)
c
c       this subroutine moves parm into coeft for lagrange multiplier
c       tests.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ncoefs=(ntimes*nys)+(nbetas*nys)+(nys*nsegm2s)+
     & (nys*(nfacs+nys))+(nys*(nys-1)/2),nextras=(nys*nfacs)+(nys*nys)+
     & (nys*(nys-1)/2),ncoefts=ncoefs+nextras)
      implicit real*8(a-h,o-z)
      dimension parm(ncoefts)
      common/coefts/coeft(ncoefts),icoeft(ncoefts),ncoeft
      iparm=0
      do 2 i=1,ncoeft
         if(icoeft(i).eq.0)then
            iparm=iparm+1
            coeft(i)=parm(iparm)
            endif
   2     continue
      return
      end
c
c
      subroutine parset211
c
c       this subroutine puts coef into theta for the simple model.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,nsegp1s=
     & nsegs+1,ntimes=34,nys=6,ncoefs=(ntimes*nys)+(nbetas*nys)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      ij=0
      do 3 i=1,ntime
         ij=ij+1
   3     fixeff(1,i)=coef(ij)
      do 2 i=1,nx
         ij=ij+1
   2     beta(1,i)=coef(ij)
      ij=ij+1
      rho(1,1)=exp(coef(ij))
      rho(1,1)=(2.d0*rho(1,1)/(1.d0+rho(1,1)))-1.d0
      ij=ij+1
      sige=exp(coef(ij))
      thresh(1,1)=-100.d0
      thresh(1,2)=0.d0
      ij=ij+1
      thresh(1,3)=exp(coef(ij))
      ij=ij+1
      thresh(1,4)=thresh(1,3)+exp(coef(ij))
      ij=ij+1
      thresh(1,5)=thresh(1,4)+exp(coef(ij))
      thresh(1,6)=100.d0
      return
      end
c
c
      subroutine parset212(penalty)
c
c       this subroutine puts coef into theta for the complex model.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,nsegp1s=
     & nsegs+1,ntimes=34,nys=6,ncoefs=(ntimes*nys)+(nbetas*nys)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension chol(nys,nys)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      common/wald/icowald(2,ncoefs)
      ij=0
      if(ifwvdum.eq.1)then
         do 11 i=1,ny
            do 12 j=1,ntime
               ij=ij+1
  12           fixeff(i,j)=coef(ij)
  11        continue
         endif
      do 3 i=1,ny
         do 2 j=1,nx
            ij=ij+1
   2        beta(i,j)=coef(ij)
   3     continue
      if(method.eq.6)then
         do 14 i=1,ij
            icowald(1,i)=1
  14        icowald(2,i)=1
         endif
      ijt=ij
      do 4 i=1,ny
         do 5 j=1,ny
            ij=ij+1
            rho(i,j)=exp(coef(ij))
   5        rho(i,j)=(2.d0*rho(i,j)/(1.d0+rho(i,j)))-1.d0
   4     continue
      if(method.eq.6)then
         do 15 i=1,ny
            do 16 j=1,ny
               ijt=ijt+1
               icowald(1,ijt)=0
               icowald(2,ijt)=0
               if(i.eq.j)icowald(2,ijt)=1
  16           continue
  15        continue
         endif
      do 6 i=1,nfac
         do 7 j=1,ny 
            ij=ij+1
   7        alamb(j,i)=coef(ij)
   6     continue
      if(method.eq.6)then
         do 17 i=1,ny
            ijt=ijt+1
            icowald(1,ijt)=0
  17        icowald(2,ijt)=1
         endif
      do 8 i=2,ny
         im1=i-1
         do 9 j=1,im1
            ij=ij+1
   9        chol(i,j)=coef(ij)
   8     continue
      if(method.eq.6)then
         do 18 i=1,ny
            im1=i-1
            do 19 j=1,im1
               ijt=ijt+1
               icowald(1,ijt)=0
  19           icowald(2,ijt)=0
  18        continue
         endif
      call omegamake1(chol,omega,penalty)
      do 10 i=1,ny
         thresh(i,1)=-100.d0
         thresh(i,2)=0.d0
         ij=ij+1
         thresh(i,3)=exp(coef(ij))
         nsegm2=nseg(i)-1
         if(nsegm2.ge.3)then
            do 13 j=3,nsegm2
               ij=ij+1
  13           thresh(i,j+1)=thresh(i,j)+exp(coef(ij))
            endif
  10     thresh(i,nseg(i)+1)=100.d0
      if(method.eq.6)then
         do 20 i=1,ny
            do 21 j=1,nsegm2
               ijt=ijt+1
               icowald(1,ijt)=1
  21           icowald(2,ijt)=1
  20        continue
         endif
      return
      end
c
c
      subroutine parset212t(penalty)
c
c       this subroutine puts coeft into theta for the complex model for
c       lagrange multiplier tests.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,nsegp1s=
     & nsegs+1,ntimes=34,nys=6,ncoefs=(ntimes*nys)+(nbetas*nys)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2),nextras=
     & (nys*nfacs)+(nys*nys)+(nys*(nys-1)/2),ncoefts=ncoefs+nextras)
      implicit real*8(a-h,o-z)
      dimension chol(nys,nys),cholt(nys,nys)
      common/coefts/coeft(ncoefts),icoeft(ncoefts),ncoeft
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/lmparms/alambt(nys,nfacs),omegat(nys,nys),
     & rhoptt(nys,nys,ntimes),rhot(nys,nys),iswitch(ntimes),lmflag
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      ij=0
      if(ifwvdum.eq.1)then
         do 11 i=1,ny
            do 12 j=1,ntime
               ij=ij+1
  12           fixeff(i,j)=coeft(ij)
  11        continue
         endif
      do 3 i=1,ny
         do 2 j=1,nx
            ij=ij+1
   2        beta(i,j)=coeft(ij)
   3     continue
      do 4 i=1,ny
         do 5 j=1,ny
            ij=ij+1
            rho(i,j)=exp(coeft(ij))
   5        rho(i,j)=(2.d0*rho(i,j)/(1.d0+rho(i,j)))-1.d0
   4     continue
      do 6 i=1,nfac
         do 7 j=1,ny 
            ij=ij+1
   7        alamb(j,i)=coeft(ij)
   6     continue
      do 8 i=2,ny
         im1=i-1
         do 9 j=1,im1
            ij=ij+1
   9        chol(i,j)=coeft(ij)
   8     continue
      call omegamake1(chol,omega,penalty)
      do 10 i=1,ny
         thresh(i,1)=-100.d0
         thresh(i,2)=0.d0
         ij=ij+1
         thresh(i,3)=exp(coeft(ij))
         nsegm2=nseg(i)-1
         if(nsegm2.ge.3)then
            do 13 j=3,nsegm2
               ij=ij+1
  13           thresh(i,j+1)=thresh(i,j)+exp(coeft(ij))
            endif
  10     thresh(i,nseg(i)+1)=100.d0
      do 14 i=1,ny
         do 15 j=1,ny
            ij=ij+1
            rhot(i,j)=exp(coeft(ij))
  15        rhot(i,j)=(2.d0*rhot(i,j)/(1.d0+rhot(i,j)))-1.d0
  14     continue
      do 16 i=1,nfac
         do 17 j=1,ny 
            ij=ij+1
  17        alambt(j,i)=coeft(ij)
  16     continue
      do 18 i=2,ny
         im1=i-1
         do 19 j=1,im1
            ij=ij+1
  19        cholt(i,j)=coeft(ij)
  18     continue
      call omegamake1(cholt,omegat,penalty)
      return
      end
c
c
      subroutine parset22(penalty,dpenalty,nparm,ifder)
c
c       this subroutine puts coef into theta for the complex model.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,nsegp1s=
     & nsegs+1,ntimes=34,nys=6,ncoefs=(ntimes*nys)+(nbetas*nys)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension chol(nys,nys),dchol(nys,nys,ncoefs),dpenalty(ncoefs)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      ij=0
      if(ifwvdum.eq.1)then
         do 11 i=1,ny
            do 12 j=1,ntime
               ij=ij+1
  12           fixeff(i,j)=coef(ij)
  11        continue
         endif
      do 3 i=1,ny
         do 2 j=1,nx
            ij=ij+1
   2        beta(i,j)=coef(ij)
   3     continue
      do 4 i=1,ny
         do 5 j=1,ny
            ij=ij+1
            rho(i,j)=exp(coef(ij))
   5        rho(i,j)=(2.d0*rho(i,j)/(1.d0+rho(i,j)))-1.d0
   4     continue
      do 6 i=1,nfac
         do 7 j=1,ny 
            ij=ij+1
   7        alamb(j,i)=coef(ij)
   6     continue
      do 8 i=2,ny
         im1=i-1
         do 9 j=1,im1
            ij=ij+1
   9        chol(i,j)=coef(ij)
   8     continue
      do 10 i=1,ny
         thresh(i,1)=-100.d0
         thresh(i,2)=0.d0
         ij=ij+1
         thresh(i,3)=exp(coef(ij))
         nsegm2=nseg(i)-1
         if(nsegm2.ge.3)then
            do 13 j=3,nsegm2
               ij=ij+1
  13           thresh(i,j+1)=thresh(i,j)+exp(coef(ij))
            endif
  10     thresh(i,nseg(i)+1)=100.d0
      if(ifder.eq.1)then
         do 15 i=1,ncoef
            do 16 j=1,ny
               do 17 k=1,ny
  17              dfixeff(j,k,i)=0.d0
  16           continue
            do 18 j=1,ny
               do 19 k=1,nx
  19              dbeta(j,k,i)=0.d0
  18           continue
            do 20 j=1,ny
               do 21 k=1,ny
  21              drho(j,k,i)=0.d0
  20           continue 
            do 22 j=1,nfac
               do 23 k=1,ny
  23              dalamb(k,j,i)=0.d0
  22           continue
            do 24 j=1,ny
               do 25 k=1,ny
  25              dchol(j,k,i)=0.d0
  24           continue
            do 26 j=1,ny
               nsegj=nseg(j)
               do 27 k=1,nsegj
  27              dthresh(j,k,i)=0.d0 
  26           continue
  15        continue
         if(ifwvdum.eq.1)then
            n1=ny*ntime
            n2=n1+(ny*nx)
            n3=n2+(ny*ny)
            n4=n3+(nfac*ny)
            n5=n4+((ny*(ny-1))/2)
            iparm=0
            ij=0
            do 14 i=1,ncoef
               ij=ij+1
               if(icoef(ij).eq.0)then
                  i1=1
                  i2=0
                  iparm=iparm+1
                  if(ij.le.n1)then
                     i2=i2+1
                     if(i2.gt.ny)then
                        i1=i1+1
                        i2=1
                        endif
                     dfixeff(i1,i2,iparm)=1.d0
                     endif
                  if((ij.gt.n1).and.(ij.le.n2))then
                     i2=i2+1
                     if(i2.gt.nx)then
                        i1=i1+1
                        i2=1
                        endif
                     dbeta(i2,i1,iparm)=1.d0
                     endif
                  if((ij.gt.n2).and.(ij.le.n3))then
                     i2=i2+1
                     if(i2.gt.ny)then
                        i1=i1+1
                        i2=1
                        endif
                     erho=exp(coef(ij))
                     drho(i1,i2,iparm)=((2.d0*erho)/(1.d0+erho)-(2.d0*
     &                ((erho/(1.d0+erho))**.5d0)))
                     endif
                  if((ij.gt.n3).and.(ij.le.n4))then
                     i2=i2+1
                     if(i2.gt.nfac)then
                        i1=i1+1
                        i2=1
                        endif
                     dalamb(i1,i2,iparm)=1.d0
                     endif
                  if((ij.gt.n4).and.(ij.le.n5))then
                     i2=i2+1
                     if(i2.gt.(i2-1))then
                        i1=i1+1
                        i2=1
                        endif
                     dchol(i1,i2,iparm)=1.d0
                     endif
                  if(ij.gt.n5)then
                     do 28 j=1,ny
                        do 29 k=3,nsegm2
                           do 30 l=k,nsegm2 
  30                          dthresh(j,l,iparm)=exp(coef(ij))
  29                       continue
  28                    continue
                     endif
                  endif
  14           continue
            endif
         if(ifwvdum.eq.0)then
            n1=(ny*nx)
            n2=n1+(ny*ny)
            n3=n2+(nfac*ny)
            n4=n3+((ny*(ny-1))/2)
            iparm=0
            ij=0
            do 31 i=1,ncoef
               ij=ij+1
               if(icoef(ij).eq.0)then
                  i1=1
                  i2=0
                  iparm=iparm+1
                  if(ij.le.n1)then
                     i2=i2+1
                     if(i2.gt.nx)then
                        i1=i1+1
                        i2=1
                        endif
                     dbeta(i2,i1,iparm)=1.d0
                     endif
                  if((ij.gt.n1).and.(ij.le.n2))then
                     i2=i2+1
                     if(i2.gt.ny)then
                        i1=i1+1
                        i2=1
                        endif
                     erho=exp(coef(ij))
                     drho(i1,i2,iparm)=((2.d0*erho)/(1.d0+erho)-(2.d0*
     &                ((erho/(1.d0+erho))**.5d0)))
                     endif
                  if((ij.gt.n2).and.(ij.le.n3))then
                     i2=i2+1
                     if(i2.gt.nfac)then
                        i1=i1+1
                        i2=1
                        endif
                     dalamb(i1,i2,iparm)=1.d0
                     endif
                  if((ij.gt.n3).and.(ij.le.n4))then
                     i2=i2+1
                     if(i2.gt.(i2-1))then
                        i1=i1+1
                        i2=1
                        endif
                     dchol(i1,i2,iparm)=1.d0
                     endif
                  if(ij.gt.n4)then
                     do 33 j=1,ny
                        do 34 k=3,nsegm2 
                           do 32 l=k,nsegm2 
  32                          dthresh(j,l,iparm)=exp(coef(ij))
  34                       continue
  33                    continue
                     endif
                  endif
  31           continue
            endif
         endif
      call omegamake2(chol,dchol,omega,domega,penalty,dpenalty,nparm,
     & ifder)
      return
      end
c
c
      subroutine preplik11
c
c       this subroutine constructs the nonrandom part of the dependent 
c       variable for the simple model.
c
      parameter(nbetas=20,nfacs=2,nrs=10,nsegs=8,nsegm2s=nsegs-2,
     & nsegp1s=nsegs+1,ntimes=34,nxs=nbetas,nys=6,ncoefs=(ntimes*nys)+
     & (nbetas*nys)+(nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      do 2 i=istart,iend 
         xb(1,i)=0.d0
         if(ifwvdum.eq.1)xb(1,i)=xb(1,i)+fixeff(1,i)
         do 3 j=1,nx
   3        xb(1,i)=xb(1,i)+(beta(1,j)*x(i,j))
   2     continue
      istartp1=istart+1
      if(iend.ge.istartp1)then
         do 4 i=istartp1,iend 
   4        rhop(1,1,i)=rho(1,1)**((time(i)-time(i-1))/365.d0)
         endif
      return
      end
c
c
      subroutine preplik12
c
c       this subroutine constructs the nonrandom part of the dependent 
c       variable for the complex model.
c
      parameter(nbetas=20,nfacs=2,nrs=10,nsegs=8,nsegm2s=nsegs-2,
     & nsegp1s=nsegs+1,ntimes=34,nxs=nbetas,nys=6,ncoefs=(ntimes*nys)+
     & (nbetas*nys)+(nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      data pi/3.1415926535898d0/
!$omp threadprivate(/datat/)
      do 2 i=istart,iend 
         do 5 j=1,ny
            xb(j,i)=0.d0
            if(ifwvdum.eq.1)xb(j,i)=xb(j,i)+fixeff(j,i)
            do 3 k=1,nx
   3           xb(j,i)=xb(j,i)+(beta(j,k)*x(i,k))
   5        continue
   2     continue
      istartp1=istart+1
      if(iend.ge.istartp1)then
         do 4 i=istartp1,iend 
            do 6 j=1,ny
               do 7 k=1,ny
                  if(rho(j,k).ne.0.d0)then
                     if((irhobessel.eq.0).or.(rho(j,k).gt.0.d0))then
                        fac=rho(j,k)
                        afac=abs(fac)
                        rfac=fac/afac
                        if(idata.eq.1)rhop(j,k,i)=rfac*(afac**(12.d0*
     &                   (time(i)-time(i-1))/365.d0))
                        if((idata.eq.2).or.(idata.eq.3))rhop(j,k,i)=
     &                   rfac*(afac**(time(i)-time(i-1)))
                        endif
                     if((irhobessel.eq.1).and.(rho(j,k).lt.0.d0))then
                        fac=rho(j,k)
                        afac=abs(fac)
                        if(idata.eq.1)fact=12.d0*(time(i)-time(i-1))/
     &                   365.d0
                        if((idata.eq.2).or.(idata.eq.3))fact=time(i)-
     &                   time(i-1)
                        rhop(j,k,i)=cos(pi*fact/2.d0)*
     &                   (afac**fact)
                        endif
                     endif
                  if(rho(j,k).eq.0.d0)rhop(j,k,i)=0.d0
   7              continue
   6           continue
   4        continue
         endif
      return
      end
c
c
      subroutine preplik12t
c
c       this subroutine constructs the nonrandom part of the dependent 
c       variable for the complex model for lagrange multiplier tests.
c
      parameter(nbetas=20,nfacs=2,nrs=10,nsegs=8,nsegm2s=nsegs-2,
     & nsegp1s=nsegs+1,ntimes=34,nxs=nbetas,nys=6,ncoefs=(ntimes*nys)+
     & (nbetas*nys)+(nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/lmparms/alambt(nys,nfacs),omegat(nys,nys),
     & rhoptt(nys,nys,ntimes),rhot(nys,nys),iswitch(ntimes),lmflag
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      data pi/3.1415926535898d0/
      do 2 i=istart,iend 
         do 5 j=1,ny
            xb(j,i)=0.d0
            if(ifwvdum.eq.1)xb(j,i)=xb(j,i)+fixeff(j,i)
            do 3 k=1,nx
   3           xb(j,i)=xb(j,i)+(beta(j,k)*x(i,k))
   5        continue
   2     continue
      do 11 i=istart,iend
  11     iswitch(i)=0
      if((lmflag.eq.1).and.((x(istart,4).eq.1.d0).or.
     & (x(istart,5).eq.1.d0)))then
         do 12 i=istart,iend
  12        iswitch(i)=1
         endif
      if(lmflag.eq.2)then
         do 13 i=istart,iend
            if(xextrat(i).lt.0.d0)iobsrv(i)=0
            if(xextrat(i).eq.0.d0)iswitch(i)=0
            if(xextrat(i).eq.1.d0)iswitch(i)=1
  13        continue
         endif
      istartp1=istart+1
      if(iend.ge.istartp1)then
         do 4 i=istartp1,iend 
            if((lmflag.eq.3).and.(i.gt.8))iswitch(i)=1
            if(lmflag.eq.4)then
               if((i.ge.4).and.(i.le.6))iswitch(i)=1
               if((i.eq.12).or.(i.eq.13))iswitch(i)=1
               if((i.eq.18).or.(i.eq.19))iswitch(i)=1
               endif
            if(iswitch(i).eq.0)then
               do 6 j=1,ny
                  do 7 k=1,ny
                     if(rho(j,k).ne.0.d0)then
                        if((irhobessel.eq.0).or.(rho(j,k).gt.0.d0))then
                           fac=rho(j,k)
                           afac=abs(fac)
                           rfac=fac/afac
                           rhop(j,k,i)=rfac*(afac**(12.d0*(time(i)-
     &                      time(i-1))/365.d0))
                           endif
                        if((irhobessel.eq.1).and.(rho(j,k).lt.0.d0))
     &                   then
                           fac=rho(j,k)
                           afac=abs(fac)
                           fact=12.d0*(time(i)-time(i-1))/365.d0
                           rhop(j,k,i)=cos(pi*fact/2.d0)*(afac**fact)
                           endif
                        endif
                     if(rho(j,k).eq.0.d0)rhop(j,k,i)=0.d0
   7                 continue
   6              continue
               endif
            if(iswitch(i).eq.1)then
               do 9 j=1,ny
                  do 10 k=1,ny
                     if(rhot(j,k).ne.0.d0)then
                        if((irhobessel.eq.0).or.(rhot(j,k).gt.0.d0))
     &                   then
                           fac=rhot(j,k)
                           afac=abs(fac)
                           rfac=fac/afac
                           rhoptt(j,k,i)=rfac*(afac**(12.d0*(time(i)-
     &                      time(i-1))/365.d0))
                           endif
                        if((irhobessel.eq.1).and.(rhot(j,k).lt.0.d0))
     &                   then
                           fac=rhot(j,k)
                           afac=abs(fac)
                           fact=12.d0*(time(i)-time(i-1))/365.d0
                           rhoptt(j,k,i)=cos(pi*fact/2.d0)*(afac**fact)
                           endif
                        endif
                     if(rhot(j,k).eq.0.d0)rhoptt(j,k,i)=0.d0
  10                 continue
   9              continue
               endif
   4        continue
         endif
      return
      end
c
c
      subroutine preplik2(nparm,ifder)
c
c       this subroutine constructs the nonrandom part of the dependent 
c       variable for the complex model with partial derivatives.
c
      parameter(nbetas=20,nfacs=2,nrs=10,nsegs=8,nsegm2s=nsegs-2,
     & nsegp1s=nsegs+1,ntimes=34,nxs=nbetas,nys=6,ncoefs=(ntimes*nys)+
     & (nbetas*nys)+(nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      data pi/3.1415926535898d0/
      do 2 i=istart,iend 
         do 5 j=1,ny
            xb(j,i)=0.d0
            if(ifwvdum.eq.1)xb(j,i)=xb(j,i)+fixeff(j,i)
            do 3 k=1,nx
   3           xb(j,i)=xb(j,i)+(beta(j,k)*x(i,k))
   5        continue
   2     continue
      if(ifder.eq.1)then
         do 14 i=istart,iend
            do 11 j=1,ny
               do 12 k=1,nparm
                  dxb(j,i,k)=0.d0
                  if(ifwvdum.eq.1)dxb(j,i,k)=dxb(j,i,k)+dfixeff(j,i,k)
                  do 13 l=1,nx
  13                 dxb(j,i,k)=dxb(j,i,k)+(dbeta(j,l,k)*x(i,l))
  12              continue
  11           continue
  14        continue
         endif
      istartp1=istart+1
      if(iend.ge.istartp1)then
         do 4 i=istartp1,iend 
            do 6 j=1,ny
               do 7 k=1,ny
                  if(rho(j,k).ne.0.d0)then
                     if((irhobessel.eq.0).or.(rho(j,k).gt.0.d0))then
                        fac=rho(j,k)
                        afac=abs(fac)
                        rfac=fac/afac
                        rhop(j,k,i)=rfac*(afac**(12.d0*(time(i)-
     &                   time(i-1))/365.d0))
                        endif
                     if((irhobessel.eq.1).and.(rho(j,k).lt.0.d0))then
                        fac=rho(j,k)
                        afac=abs(fac)
                        fact=12.d0*(time(i)-time(i-1))/365.d0
                        rhop(j,k,i)=cos(pi*fact/2.d0)*(afac**fact)
                        endif
                     if(ifder.eq.1)then
                        if((irhobessel.eq.0).or.(rho(j,k).gt.0.d0))then
                           do 9 l=1,nparm
   9                          drhop(j,k,i,l)=rhop(j,k,i)*drho(j,k,l)
                           endif
                        if((irhobessel.eq.1).and.(rho(j,k).lt.0.d0))
     &                   then
                           do 10 l=1,nparm
  10                          drhop(j,k,i,l)=cos(pi*fact/2.d0)*
     &                         rhop(j,k,i)*drho(j,k,l)
                           endif
                        endif
                     endif
                  if(rho(j,k).eq.0.d0)then
                     rhop(j,k,i)=0.d0
                     if(ifder.eq.1)then
                        do 8 l=1,nparm
   8                       drhop(j,k,i,l)=0.d0
                        endif
                     endif
   7              continue
   6           continue
   4        continue
         endif
      return
      end
c
c
      subroutine probmake1(arg1,arg2,f)
c
c       this subroutine computes the log ordered probit term at the 
c       heart of the log likelihood contribution.
c
      implicit real*8(a-h,o-z)
      data big1/4.d0/
      data pi/3.1415d0/
      data small1/-4.d0/
      data small2/-5.d0/
      if((arg2.ge.small1).and.(arg1.le.big1))then
c
c       compute for case not in extremes.
c
         call anormc(arg2,farg2)
         call anormc(arg1,farg1)
         f=log(farg2-farg1)
         return
         endif
      if(arg2.lt.small1)then
         delta=arg2-arg1
         fac1=delta*((2.d0*arg2)-delta)/2.d0
         fac2=log(-1.d0*arg2)
         fac3=arg2*arg2/2.d0
         fac4=log(2.d0*pi)/2.d0
         if(arg2.ge.small2)then
            f1=-1.d0*(exp(fac1)+fac2+fac3+fac4)
            call anormc(arg2,farg2)
            call anormc(arg1,farg1)
            f2=log(farg2-farg1)
            weight=abs(arg2-small2)
            f=(weight*f1)+((1.d0-weight)*f2)
            return
            endif
         if(arg2.lt.small2)then
            f=-1.d0*(exp(fac1)+fac2+fac3+fac4)
            return
            endif
         endif
      if(arg1.gt.big1)then
         arg2t=-1.d0*arg1
         arg1t=-1.d0*arg2
         delta=arg2t-arg1t
         fac1=delta*((2.d0*arg2t)-delta)/2.d0
         fac2=log(-1.d0*arg2t)
         fac3=arg2t*arg2t/2.d0
         fac4=log(2.d0*pi)/2.d0
         if(arg2t.ge.small2)then
            f1=-1.d0*(exp(fac1)+fac2+fac3+fac4)
            call anormc(arg2t,farg2)
            call anormc(arg1t,farg1)
            f2=log(farg2-farg1)
            weight=abs(arg2t-small2)
            f=(weight*f1)+((1.d0-weight)*f2)
            return
            endif
         if(arg2t.lt.small2)then
            f=-1.d0*(exp(fac1)+fac2+fac3+fac4)
            return
            endif
         endif
      return
      end
c
c
      subroutine probmake2(arg1,darg1,arg2,darg2,f,df,nparm,ifder)
c
c       this subroutine computes the log ordered probit term and
c       partial derivatives at the heart of the log likelihood 
c       contribution.
c
      parameter(nbetas=20,nfacs=2,nsegs=8,nsegm2s=nsegs-2,ntimes=34,
     & nys=6,ncoefs=(ntimes*nys)+(nbetas*nys)+(nys*nsegm2s)+(nys*
     & (nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      dimension darg1(ncoefs),darg1t(ncoefs),darg2(ncoefs),
     & darg2t(ncoefs),ddelta(ncoefs),df(ncoefs),dfac1(ncoefs),
     & dfac2(ncoefs),dfac3(ncoefs),dfarg1(ncoefs),dfarg2(ncoefs),
     & df1(ncoefs),df2(ncoefs),dweight(ncoefs)
      data big1/4.d0/
      data pi/3.1415d0/
      data small1/-4.d0/
      data small2/-5.d0/
      if((arg2.ge.small1).and.(arg1.le.big1))then
c
c       compute for case not in extremes.
c
         call anormc(arg2,farg2)
         call anormc(arg1,farg1)
         if(ifder.eq.1)then
            call anormd(arg2,dfarg2t)
            call anormd(arg1,dfarg1t)
            do 2 i=1,nparm
               dfarg2(i)=darg2(i)*dfarg2t
   2           dfarg1(i)=darg1(i)*dfarg1t
            endif
         f=log(farg2-farg1)
         if(ifder.eq.1)then
            fac=farg2-farg1
            do 3 i=1,nparm
   3           df(i)=(dfarg2(i)-dfarg1(i))/fac
            endif
         return
         endif
      if(arg2.lt.small1)then
         delta=arg2-arg1
         if(ifder.eq.1)then
            do 4 i=1,nparm
   4           ddelta(i)=darg2(i)-darg1(i)
            endif
         fac1a=((2.d0*arg2)-delta)/2.d0
         fac1=delta*fac1a
         fac2=log(-1.d0*arg2)
         fac3=arg2*arg2/2.d0
         fac4=log(2.d0*pi)/2.d0
         if(ifder.eq.1)then
            do 5 i=1,nparm
               dfac1(i)=(ddelta(i)*fac1a)+(delta*(darg2(i)-ddelta(i)))
               dfac2(i)=-1.d0*darg2(i)/arg2
   5           dfac3(i)=arg2*darg2(i)
            endif
         if(arg2.ge.small2)then
            f1=-1.d0*(exp(fac1)+fac2+fac3+fac4)
            call anormc(arg2,farg2)
            call anormc(arg1,farg1)
            f2=log(farg2-farg1)
            weight=abs(arg2-small2)
            f=(weight*f1)+((1.d0-weight)*f2)
            if(ifder.eq.1)then
               call anormd(arg2,dfarg2t)
               call anormd(arg1,dfarg1t)
               do 6 i=1,nparm
                  df1(i)=-1.d0*((exp(fac1)*dfac1(i))+dfac2(i)+dfac3(i))
                  df2(i)=((dfarg2t*darg2(i))-(dfarg1t*darg1(i)))/
     &             (farg2-farg1)
                  if((arg2-small2).gt.0.d0)dweight(i)=darg2(i)
                  if((arg2-small2).lt.0.d0)dweight(i)=-1.d0*darg2(i)
   6              df(i)=(weight*df1(i))+(dweight(i)*f1)-(weight*
     &             df2(i))-(f2*dweight(i))
               endif
            return
            endif
         if(arg2.lt.small2)then
            f=-1.d0*(exp(fac1)+fac2+fac3+fac4)
            if(ifder.eq.1)then
               do 7 i=1,nparm
   7              df(i)=-1.d0*((exp(fac1)*dfac1(i))+dfac2(i)+dfac3(i))
               endif
            return
            endif
         endif
      if(arg1.gt.big1)then
         arg2t=-1.d0*arg1
         arg1t=-1.d0*arg2
         delta=arg2t-arg1t
         fac1a=((2.d0*arg2t)-delta)/2.d0
         fac1=delta*fac1a
         fac2=log(-1.d0*arg2t)
         fac3=arg2t*arg2t/2.d0
         fac4=log(2.d0*pi)/2.d0
         if(ifder.eq.1)then
            do 8 i=1,nparm
               darg2t(i)=-1.d0*darg1(i)
               darg1t(i)=-1.d0*darg2(i)
               ddelta(i)=darg2t(i)-darg1t(i)
               dfac1(i)=(ddelta(i)*fac1a)+(delta*(darg2(i)-ddelta(i)))
               dfac2(i)=-1.d0*darg2(i)/arg2
   8           dfac3(i)=arg2*darg2(i)
            endif
         if(arg2t.ge.small2)then
            f1=-1.d0*(exp(fac1)+fac2+fac3+fac4)
            call anormc(arg2t,farg2)
            call anormc(arg1t,farg1)
            f2=log(farg2-farg1)
            weight=abs(arg2t-small2)
            f=(weight*f1)+((1.d0-weight)*f2)
            if(ifder.eq.1)then
               call anormd(arg2t,dfarg2t)
               call anormd(arg1t,dfarg1t)
               do 9 i=1,nparm
                  df1(i)=-1.d0*(exp(fac1)*dfac1(i)+dfac2(i)+dfac3(i))
                  df2(i)=((dfarg2t*darg2t(i))-(dfarg1t*darg1t(i)))/
     &             (farg2-farg1)
                  if(arg2t.ge.small2)dweight(i)=darg2t(i)
                  if(arg2t.lt.small2)dweight(i)=-1.d0*darg2t(i)
   9              df(i)=(f1*dweight(i)+(weight*df1(i))-(f2*dweight(i))-
     &             (weight*df2(i)))
               endif
            return
            endif
         if(arg2t.lt.small2)then
            f=-1.d0*(exp(fac1)+fac2+fac3+fac4)
            if(ifder.eq.1)then
               do 10 i=1,nparm
  10              df(i)=-1.d0*((exp(fac1)*dfac1(i))+dfac2(i)+dfac3(i))
               endif
            return
            endif
         endif
      return
      end
c
c
      function randfn(idum)
      dimension v(97)
      save y,v
      if (idum .ne. 0) then
        do k=1,97
          call random_number(v(k))
        enddo
          call random_number(y)
      endif
      j=1+int(97*y)
      y = v(j)
      randfn=v(j)
      call random_number(v(j))
      return
      end
c
c
      subroutine utrunc(u)
c
c       this subroutine censors uniform random variables at both ends.
c
      implicit real*8(a-h,o-z)
      real big,small,u
      data big/0.9999/
      data small/0.0001/
      if(u.gt.big)u=big
      if(u.lt.small)u=small
      return
      end
c
c
      subroutine vardecomp
c
c       this subroutine performs a variance decomposition based on
c       parameter estimates.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nrs=10,nsegs=8,
     & nsegm2s=nsegs-2,nsegp1s=nsegs+1,ntimes=34,nxs=nbetas,nys=6,
     & ny3s=nys*3,ncoefs=(ntimes*nys)+(nbetas*nys)+(nys*nsegm2s)+
     & (nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      character*8 adecomp(4),ay2(ny3s)
      real xextra,xv
      dimension alambvar(nys),deltabar(nys,3),deltavar(nys),
     & omegacon(nys,2),omegavar(nys),xbbar(nys,3),xbvar(nys)
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datat/drhop(nys,nys,ntimes,ncoefs),dxb(nys,ntimes,ncoefs),
     & er(nrs),er2(nrs,nfacs),rhop(nys,nys,ntimes),ur(nrs,nys,ntimes),
     & x(ntimes,nxs),xb(nys,ntimes),xextrat(ntimes),ystar(nys,ntimes),
     & iend,iobsrv(ntimes),istart,iy(nys,ntimes)
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/theta/alamb(nys,nfacs),beta(nys,nbetas),
     & dalamb(nys,nfacs,ncoefs),dbeta(nys,nbetas,ncoefs),
     & dfixeff(nys,ntimes,ncoefs),domega(nys,nys,ncoefs),
     & drho(nys,nys,ncoefs),dthresh(nys,nsegp1s,ncoefs),
     & fixeff(nys,ntimes),omega(nys,nys),rho(nys,nys),sige,
     & thresh(nys,nsegp1s)
      data adecomp/'x-beta','delta','lambda','omega'/
      data ay2/'duty to ','vote',' ','satisfie','d w/ uk ','democracy',
     & 'trust mp','s',' ','efficacy',': unders','tand','efficacy',
     & ': pol ca','re',' ',' ',' '/
c
c       initialize.
c
      do 6 i=1,ny
         do 7 j=1,3
   7        xbbar(i,j)=0.d0
   6     continue
      do 9 i=1,ny
         do 10 j=1,3
  10        deltabar(i,j)=0.d0
   9     continue
      do 13 i=1,ny
  13     alambvar(i)=0.d0
      do 21 i=1,ny
         do 22 j=1,2
  22        omegacon(i,j)=0.d0
  21     continue
      call parset212(penalty)
c
c       compute alamb-hat var.
c
      do 11 i=1,ny
         do 12 j=1,nfac
  12        alambvar(i)=alambvar(i)+(alamb(i,j)*alamb(i,j))
  11     continue
c
c       loop over observations.
c
      do 2 i=1,nobs
         if(((i/10000)*10000).eq.i)write(6,106)i
 106     format(1x,'i = ',i8)
         if(incluv(i).eq.1)then
            call obsmake2(i)
            call preplik12
            istartp1=istart+1
c
c       update omega contribution.
c
            do 17 j=1,ny
               omegacont=1.d0
               do 25 k=istartp1,iend
                  km1=k-1
                  do 18 l=1,km1
                     do 19 m=1,ny
                        rhopt=1.d0
                        do 20 n=l,km1
  20                       rhopt=rhopt*rhop(j,m,n)
  19                    omegacont=omegacont+(rhopt*rhopt)
  18                 continue
  25              continue
               omegacon(j,1)=omegacon(j,1)+1.d0
  17           omegacon(j,2)=omegacon(j,2)+omegacont
            do 3 j=istartp1,iend
c
c       update delta-hat bar.
c
               do 8 k=1,ny
                  deltabar(k,1)=deltabar(k,1)+1.d0
                  deltabar(k,2)=deltabar(k,2)+fixeff(k,j)
   8              deltabar(k,3)=deltabar(k,3)+(fixeff(k,j)*
     &             fixeff(k,j))
c
c       update x-beta-hat bar.
c
               do 4 k=1,ny
                  xb(k,j)=fixeff(k,j)
                  do 5 l=1,nx
   5                 xb(k,j)=xb(k,j)+(beta(k,l)*x(j,l))
                  xbbar(k,1)=xbbar(k,1)+1.d0
                  xbbar(k,2)=xbbar(k,2)+xb(k,j)
   4              xbbar(k,3)=xbbar(k,3)+(xb(k,j)*xb(k,j))
   3           continue
            endif
   2     continue
c
c       compute delta var.
c
      do 15 i=1,ny
         deltabar(i,2)=deltabar(i,2)/deltabar(i,1)
  15     deltavar(i)=(deltabar(i,3)/deltabar(i,1))-(deltabar(i,2)*
     &    deltabar(i,2))
c
c       compute x-beta-hat var.
c
      do 16 i=1,ny
         xbbar(i,2)=xbbar(i,2)/xbbar(i,1)
  16     xbvar(i)=(xbbar(i,3)/xbbar(i,1))-(xbbar(i,2)*xbbar(i,2))
c
c       compute omega var.
c
      do 23 i=1,ny
  23     omegavar(i)=omegacon(i,2)/omegacon(i,1)
c
c       output.
c
      write(6,100)
 100  format(1x,50('='))
      write(6,101)
 101  format(1x,'variance decomposition for each dependent variable')
      do 24 i=1,ny
         ie=i*3
         ib=ie-2
         write(6,102)(ay2(j),j=ib,ie)
 102     format(1x,'variance decomposition for ',3a8)
         write(6,103)
 103     format(1x,'source of',/,1x,'variation',2x,'size',6x,
     &    'proportion')
         vartot=xbvar(i)+deltavar(i)+alambvar(i)+omegavar(i)
         prop=xbvar(i)/vartot
         write(6,104)adecomp(1),xbvar(i),prop
         prop=deltavar(i)/vartot
         write(6,104)adecomp(2),deltavar(i),prop
         prop=alambvar(i)/vartot
         write(6,104)adecomp(3),alambvar(i),prop
         prop=omegavar(i)/vartot
         write(6,104)adecomp(4),omegavar(i),prop
 104     format(1x,a8,3x,f8.3,2x,f8.4)
         write(6,105)vartot
 105     format(1x,'total:',5x,f8.3)
  24     continue
      return
      end
c
c
      subroutine vartest1
c
c       this subroutine tests to see if the log likelihood is sensitive
c       to each of the paramters and if the outer product of 
c       derivatives is positive definite.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2),ncoefss=
     & ncoefs*(ncoefs+1)/2)
      implicit real*8(a-h,o-z)
      real xextra,xv
      dimension amean(ncoefs),cov(ncoefs,ncoefs),covs(ncoefss),
     & dfv(ncoefs),eigval(ncoefs),eigvec(ncoefs,ncoefs),parm(ncoefs),
     & work1(ncoefs)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,method,
     & model,niter,niterb,niterg,nr
      data eigtol/0.1d-3/
c
c       initialize.
c
      f0=0.d0
      call bestupd(0,f0)
      call parset0(parm,nparm)
      if(model.eq.1)then
         call parset211
         call likely11(parm,nparm,f0,*3)
         endif
      if(model.eq.2)then
         call parset212(penalty)
         call likely12(parm,nparm,f0,*3)
         f0=f0+penalty
         endif
   3  continue
c
c       look for parameters with no change in log likelihood function.
c
      iparm=0
      do 4 i=1,ncoef
         if(icoef(i).eq.0)then
            iparm=iparm+1
            delta=abs(parm(iparm))/1.d2
            if(delta.gt.1.d0)delta=1.d0
            if(delta.lt.0.01d0)delta=0.01d0
            parm(iparm)=parm(iparm)+delta
            call parset1(parm)
            if(model.eq.1)then
               call parset211
               call likely11(parm,nparm,f1,*3)
               endif
            if(model.eq.2)then
               call parset212(penalty)
               call likely12(parm,nparm,f1,*3)
               f1=f1+penalty
               endif
            df=f1-f0
            if(abs(df).gt.0.d0)then
               parm(iparm)=parm(iparm)-delta
               goto 4
               endif
            parm(iparm)=parm(iparm)-delta+1.d0
            call parset1(parm)
            if(model.eq.1)then
               call parset211
               call likely11(parm,nparm,f1,*3)
               endif
            if(model.eq.2)then
               call parset212(penalty)
               call likely12(parm,nparm,f1,*3)
               f1=f1+penalty
               endif
            df=f1-f0
            if(abs(df).gt.0.d0)then
               write(6,100)i,iparm
 100           format(1x,'change of 1 needed for coef(',i3,') and',
     &          ' parm(',i3,')')
               parm(iparm)=parm(iparm)-1.d0
               goto 4
               endif
            write(6,101)i,iparm
 101        format(1x,'no change for coef(',i3,') and parm(',i3,')')
            endif
   4     continue
c
c       construct outer product of derivatives and check for positive
c       definiteness.
c
      anobs=dble(float(nobs))
      do 16 i=1,nparm
  16     amean(i)=0.d0
      do 13 i=1,nparm
         do 14 j=1,i
  14        cov(i,j)=0.d0
  13     continue
      do 5 i=1,nobs
         if(incluv(i).eq.1)then
            call parset1(parm)
            if(model.eq.1)then
               call parset211
               call obsmake1(i)
               call likobs11(f0)
               endif
            if(model.eq.2)then
               call parset212(penalty)
               call obsmake2(i)
               call likobs12(f0)
               f0=f0+penalty
               endif
            iparm=0
            do 2 j=1,ncoef
               if(icoef(j).eq.0)then
                  iparm=iparm+1
                  delta=abs(parm(iparm))/1.d2
                  if(delta.gt.1.d0)delta=1.d0
                  if(delta.lt.0.01d0)delta=0.01d0
                  parm(iparm)=parm(iparm)+delta
                  call parset1(parm)
                  if(model.eq.1)then
                     call parset211
                     call likobs11(f1)
                     endif
                  if(model.eq.2)then
                     call parset212(penalty)
                     call likobs12(f1)
                     f1=f1+penalty
                     endif
                  dfv(iparm)=(f1-f0)/delta
                  parm(iparm)=parm(iparm)-delta
                  endif
   2           continue
            do 15 j=1,iparm
  15           amean(j)=amean(j)+dfv(j)
            do 6 j=1,iparm
               do 7 k=1,j
   7              cov(j,k)=cov(j,k)+(dfv(j)*dfv(k))
   6           continue
            endif
   5     continue
c
c       adjust amean and cov.
c
      do 8 i=1,iparm
   8     amean(i)=amean(i)/anobs
      ij=0
      do 9 i=1,iparm
         do 10 j=1,i
            cov(i,j)=(cov(i,j)/anobs)-(amean(i)*amean(j))
            ij=ij+1
  10        covs(ij)=cov(i,j)
   9     continue
c
c       compute eigenvalues and associated eigenvectors.
c
      call eigrs(covs,iparm,1,eigval,eigvec,ncoefs,work1,ier)
      iflag=0
      do 11 i=1,iparm
         if(eigval(i).ge.eigtol)goto 12
         if(eigval(i).lt.eigtol)then
            if(iflag.eq.0)then
               iflag=1
               write(6,103)
 103           format(1x,'problem eigenvalues and associated',
     &          ' eigenvectors')
               write(6,104)
 104           format(6x,'eigenvalue',7x,'eigenvector')
               endif
            write(6,105)i,eigval(i),(eigvec(j,i),j=1,iparm)
 105        format(1x,i3,2x,g15.8,5(2x,f7.4),50(/,21x,5(2x,f7.4)))
            endif
  11     continue
  12  continue
      if(iflag.eq.0)write(6,106)
 106  format(1x,'no problem eigenvalues')
      return
      end
c
c
      subroutine vartest2
c
c       this subroutine tests to see if analytical partial derivatives 
c       are close to numerical partial derivatives for each parameter
c       and each observation.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      real xextra,xv
      dimension dpenalty(ncoefs),df0(ncoefs),df1(ncoefs),parm(ncoefs)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      data tol/.01d0/
      call parset0(parm,nparm)
      do 2 i=1,nobs
         write(6,101)i
 101     format(1x,'obs #: ',i8)
         if(incluv(i).eq.1)then
            call parset1(parm)
            call parset22(penalty0,dpenalty,nparm,1)
            call obsmake2(i)
            call likobs2(f0,df0,nparm,1)
            f0=f0-penalty0
            do 4 j=1,nparm
   4           df0(j)=df0(j)-dpenalty(j)
            iparm=0
            do 3 j=1,ncoef
               if(icoef(j).eq.0)then
                  iparm=iparm+1
                  delta=abs(parm(iparm))/1.d2
                  if(delta.gt.1.d0)delta=1.d0
                  if(delta.lt.0.01d0)delta=0.01d0
                  parm(iparm)=parm(iparm)+delta
                  call parset1(parm)
                  call parset22(penalty1,dpenalty,nparm,1)
                  call likobs2(f1,df1,nparm,1)
                  f1=f1-penalty1
                  do 5 k=1,nparm
   5                 df1(k)=df1(k)-dpenalty(k)
                  dfvn=(f1-f0)/delta
                  parm(iparm)=parm(iparm)-delta
                  endif
               if(df0(j).gt.df1(j))then
                  if(((df0(j)+tol).gt.dfvn).and.(dfvn.gt.(df1(j)-tol)))
     &             goto 2
                  endif
               if(df0(j).le.df1(j))then
                  if(((df0(j)-tol).le.dfvn).and.(dfvn.le.(df1(j)+tol)))
     &             goto 2
                  endif
               write(6,100)i,j,df0(j),df1(j),dfvn
 100           format(1x,'bad derivative: i = ',i8,2x,'j = ',i3,2x,
     &          'df0 and df1 = ',2(g15.8,2x),'dfvn = ',g15.8)
               stop
   3           continue
            endif
   2     continue
      return
      end
c
c
      subroutine waldtests
c
c       this subroutine performs wald tests for all of the covariance
c       terms equal to zero and all of the covariance terms except for
c       diagonal terms equal to zero.
c
      parameter(nbetas=20,nfacs=2,nobss=96529,nsegs=8,nsegm2s=nsegs-2,
     & ntimes=34,nxs=nbetas,nys=6,ncoefs=(nys*ntimes)+(nys*nbetas)+
     & (nys*nsegm2s)+(nys*(nfacs+nys))+(nys*(nys-1)/2),ncoef2s=
     & ncoefs*2)
      implicit real*8(a-h,o-z)
      character*8 awald(8)
      real xextra,xv
      dimension amean(ncoefs),ameanw(2,ncoefs),cov(ncoefs,ncoefs),
     & covw(ncoefs,ncoefs),dft(ncoefs),parm(ncoefs),work1(ncoefs),
     & work2(ncoef2s)
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/wald/icowald(2,ncoefs)
      data awald/'no covar','iance te','rms',' ','only var','iable-sp',
     & 'ecific t','erms'/
c
c       initialize.
c
      do 13 i=1,ncoef
  13     amean(i)=0.d0
      do 6 i=1,ncoef
         do 7 j=1,i
   7        cov(i,j)=0.d0
   6     continue
      call parset0(parm,nparm)
c
c       loop over observations.
c
      anobs=0.d0
      do 2 i=1,nobs
         if(((i/1000)*1000).eq.i)write(6,101)i
 101     format(1x,'i = ',i5)
         if(incluv(i).eq.1)then
            anobs=anobs+1.d0
c
c       evaluate at the parameter estimate.
c
            call parset1(parm)
            call parset212(penalty)
            call obsmake2(i)
            call likobs12(f0)
            f0=f0+penalty
c
c       loop over parameters and compute numerical derivatives.
c
            iparm=0
            do 3 j=1,ncoef
               if((icoef(j).eq.0).and.(icowald(1,j).eq.0).or.
     &          (icowald(2,j).eq.0))then
                  iparm=iparm+1
                  delta=coef(j)/100.d0
                  if(delta.gt.1.d0)delta=1.d0
                  if(delta.lt.0.01d0)delta=0.01d0
                  coef(j)=coef(j)+delta
                  call parset212(penalty)
                  call likobs12(f1)
                  f1=f1+penalty
                  dft(iparm)=(f1-f0)/delta
                  coef(j)=coef(j)-delta              
                  endif
   3           continue
c
c       update means and covariance matrices.
c
            do 4 j=1,iparm
               amean(j)=amean(j)+dft(j)
               do 5 k=1,j
   5              cov(j,k)=cov(j,k)+(dft(j)*dft(k))
   4           continue
            endif
   2     continue
      do 11 i=1,2
c
c       do for each wald test.
c
         ie=i*4
         ib=ie-3
         write(6,102)(awald(j),j=ib,ie)
 102     format(1x,'test for ',4a8)
c
c       update and output mean derivatives.
c
         write(6,100)
 100     format(1x,50('='))
         do 14 j=1,ncoef
  14        ameanw(i,j)=0.d0
         iwo=0
         jj=0
         do 8 j=1,ncoef
            if((icoef(j).eq.0).and.(icowald(i,j).eq.0))then
               jj=jj+1
               iwo=iwo+1
               ameanw(i,iwo)=amean(jj)/anobs
               endif
   8        continue
         jj=0
         do 9 j=1,ncoef
            if((icoef(j).eq.0).and.(icowald(i,j).eq.0))then
               jj=jj+1
               kk=0
               do 10 k=1,j
                  if((icoef(k).eq.0).and.(icowald(i,k).eq.0))then
                     kk=kk+1
                     covw(jj,kk)=(cov(jj,kk)/anobs)-(ameanw(i,jj)*
     &                ameanw(i,kk))
                     covw(kk,jj)=covw(jj,kk)
                     endif
  10              continue
               endif
   9        continue
c
c       invert.
c
         d1=-1.d0
         call linv3f(covw,work1,1,iwo,ncoefs,d1,d2,work2,ier)
         if(ier.ne.0)write(6,106)ier
 106     format(1x,'ier = ',i3)
         jt=0
         chi=0.d0
         do 15 j=1,ncoef 
            if(icowald(i,j).eq.0)then
               jt=jt+1
               kt=0
               do 16 k=1,ncoef
                  if(icowald(i,k).eq.0)then
                     kt=kt+1
                     chi=chi+(coef(j)*covw(jt,kt)*coef(k))
                     endif
  16              continue
               endif
  15        continue
         chi=chi*anobs
         ndeg=jt
         andeg=dble(float(ndeg))
         anorm=(chi-andeg)/((2.d0*andeg)**.5d0)
         write(6,103)chi,ndeg,anorm
 103     format(1x,'chi2 = ',g15.8,2x,'df = ',i2,2x,'normal = ',g15.8)
  11     continue
      return
      end

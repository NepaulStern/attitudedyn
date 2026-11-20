c
c       helpers for the one-evaluation oracle. not part of the r package.
c
      subroutine oracle_shared(ishared)
      integer ishared,ifirst,icache
      character*32 flag
      save ifirst,icache
      data ifirst/0/
      if(ifirst.eq.0)then
         ifirst=1
         icache=0
         flag=' '
         call getenv('ATTDYN_SHARED_MATH',flag)
         if(flag(1:1).eq.'1')icache=1
         endif
      ishared=icache
      return
      end
c
c
      subroutine oracle_dump
      parameter(nbetas=20,nfacs=2,nobss=96529,nrs=10,nsegs=8,
     & nsegm2s=nsegs-2,ntimes=34,nxs=nbetas,nys=6,ncoefs=
     & (nys*ntimes)+(nys*nbetas)+(nys*nsegm2s)+(nys*(nfacs+nys))+
     & (nys*(nys-1)/2))
      implicit real*8(a-h,o-z)
      real erv,urv,xextra,xv
      integer*4 magic,i4,nseg4(6)
      real*4 r4
      common/coefs/coef(ncoefs),icoef(ncoefs),ncoef
      common/datan/time(ntimes),nfac,nobs,nseg(nys),ntime,nx,ny
      common/datav/xextra(nobss,ntimes),xv(nobss,ntimes,nxs),
     & idv(nobss),iendv(nobss),incluv(nobss),iobsrvv(nobss,ntimes),
     & istartv(nobss),iyv(nobss,nys,ntimes)
      common/errors/erv(nobss,nfacs,nrs),urv(nobss,nrs,nys,ntimes)
      common/oracleft/ftall(nobss)
      common/progparms/idata,ifcond1,ifsattr,ifwvdum,irhobessel,
     & method,model,niter,niterb,niterg,nr
      open(21,file='oracle_body.bin',access='stream',
     & form='unformatted',status='replace')
      magic=20261120
      write(21)magic
      i4=nobs
      write(21)i4
      i4=ny
      write(21)i4
      i4=nx
      write(21)i4
      i4=ntime
      write(21)i4
      i4=nfac
      write(21)i4
      i4=nr
      write(21)i4
      i4=ncoef
      write(21)i4
      do 10 i=1,ny
   10    nseg4(i)=nseg(i)
      write(21)(nseg4(i),i=1,ny)
      write(21)(time(i),i=1,ntime)
      write(21)(coef(i),i=1,ncoef)
      do 11 i=1,ncoef
         i4=icoef(i)
         write(21)i4
   11 continue
      write(21)(ftall(i),i=1,nobs)
      do 20 i=1,nobs
         i4=incluv(i)
         write(21)i4
         i4=istartv(i)
         write(21)i4
         i4=iendv(i)
         write(21)i4
         do 21 j=1,ntime
            i4=iobsrvv(i,j)
            write(21)i4
   21    continue
         do 22 j=1,ny
            do 23 k=1,ntime
               i4=iyv(i,j,k)
               write(21)i4
   23       continue
   22    continue
         do 24 j=1,ntime
            do 25 k=1,nx
               r4=xv(i,j,k)
               write(21)r4
   25       continue
   24    continue
   20 continue
      do 30 i=1,nobs
         do 31 j=1,nr
            do 32 k=1,nfac
               r4=erv(i,k,j)
               write(21)r4
   32       continue
   31    continue
   30 continue
      do 40 i=1,nobs
         do 41 j=1,nr
            do 42 l=1,ny
               do 43 k=1,ntime
                  r4=urv(i,j,l,k)
                  write(21)r4
   43          continue
   42       continue
   41    continue
   40 continue
      close(21)
      return
      end

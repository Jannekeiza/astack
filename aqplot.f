      program aqplot
c
c     This program is designed to plot seismic traces
c     stored in .aq format.
c
      implicit none
      integer i,j,ltpb,ns,year,month,day,hour,minute
      integer npoints,ral
      real chpb,second,sampr,tshft,lft,rgt,mn,mx
      real aver,maxd,amp(10000),ptim(10000)
      real pbw,pbh,elat,elon,depth,x,z,swpol
      character*20 aqfile,odev,sname
      character label*100
      open(unit=10,file='aqplot.in',status='old')
      read(10,1)aqfile
      read(10,1)odev
      read(10,*)pbw,pbh
      read(10,*)chpb
      read(10,*)ltpb
      read(10,*)ral
1     format(a20)
      close(10)
c
c     Read in header information of .aq file
c
      open(unit=20,file=aqfile,status='old')
      read(20,*)ns
      read(20,*)elat,elon,depth
      read(20,*)year,month,day
      read(20,*)hour,minute,second
      read(20,*)sampr
      write(label,2)aqfile,elat,elon,depth,day,month,year,
     *hour,minute,second
2     format(a15,2f7.2,f6.1,'km  ',i2,'/',i2,'/',i4,' ',i2,
     *':',i2,':',f5.2)
c
c     Read in header of first station
c
      read(20,*)swpol,npoints,tshft,sname
c
c     Initiate plotting routine
c
      lft=0.0
      rgt=npoints*sampr     
      mn=0.0
      mx=ns+1
      call pgbegin(0,odev,1,1)
      pbh=pbh/pbw
      call pgpap(pbw,pbh)
      call pgsch(chpb)
      call pgslw(ltpb)
      call pgenv(lft,rgt,mn,mx,0,0)
      if(ral.eq.1)then
         call pglabel('time(s)','station number',label)
      else
         call pglabel('time(s)','station number',' ')
      endif
c
c     Now start a loop to plot out all traces
c
      do i=1,ns
c
c        Read in data
c
         if(i.ne.1)then
            read(20,*)swpol,npoints,tshft,sname
         endif
         read(20,*)(amp(j),j=1,npoints)
c
c        Now plot trace using normalization; start by removing average.
c    
         aver=0.0
         do j=1,npoints
            aver=aver+amp(j)
         enddo
         aver=aver/npoints
         do j=1,npoints
            amp(j)=amp(j)-aver
         enddo
c
c        Now find maximum deviation from average
c
         maxd=abs(amp(1))
         do j=1,npoints
            if(maxd.lt.abs(amp(j)))maxd=abs(amp(j))
         enddo
c
c        Apply normalization and DC offset for plotting
c
         if(maxd.gt.1.0e-6)then
            do j=1,npoints
               amp(j)=swpol*amp(j)/(1.333*maxd)+i
            enddo
         else
            do j=1,npoints
               amp(j)=swpol*amp(j)+i
            enddo
         endif
c
c        Apply time shift
c
         do j=1,npoints
            ptim(j)=(j-1)*sampr+tshft
         enddo
c
c        Plot trace
c
         call pgline(npoints,ptim,amp)
c
c        Apply labels to each trace
c
         x=0.02*npoints*sampr
         z=i+0.1
         call pgtext(x,z,sname)
      enddo
      call pgend
      close(20)
      stop
      end

c*********************************************************************
c
c     program tcas
c
c     Seismology Group, Research School of Earth Sciences
c     The Australian National University
c
c     Non-linear inversion scheme for differential travel times 
c
c     Uses adaptive stacking to determine optimum trace
c     alignment given some initial approximate alignment. If
c     initial alignment is model predicted (e.g. AK135), then
c     model traveltime residuals will be obtained.
c
c
c     AUTHOR:  Brian L.N. Kennett (RSES,ANU)
c              November 2002
c
c     MODIFIED: Nick Rawlinson (RSES,ANU)
c               August 2003
c
c*********************************************************************
      Program tcas
      implicit none
c*********************************************************************
c
c---------------------------------------------------------------------
c     SET UNIT NUMBERS FOR FILE INPUT AND OUTPUT
c---------------------------------------------------------------------
      integer  EVT, RUN, INI, RES, FIN
      parameter (RUN=15, INI=16, EVT=17, RES=18, FIN=19)
c
c     EVT = unit number of input seismogram file
c     RUN = unit number of file containing run information
c     INI = unit number of initial stack output file
c     FIN = unit number of final stack output file
c     RES = unit number of file containing residuals and errors
c
c----------------------------------------------------------------------
c     COMMON BLOCK DECLARATIONS
c----------------------------------------------------------------------
c
      integer MAXST,MAXD
      parameter (MAXST=100,MAXD=5000)
      integer npts,nstkwb,nstkwl,nsta,loff,nstk
      integer lim1,lim2
      real zssl,zscp,sampin
      real zu,dtcs,dtmo,stkwb,stkwl
c
      common /IST1/ npts(MAXD),nstkwb,nstkwl,nsta,loff,nstk
      common /IST2/ lim1,lim2
      common /RST1/ zu(MAXST,MAXD),zssl(MAXD),zscp(MAXD)
      common /RST2/ dtcs(MAXST),dtmo(MAXST),sampin,stkwb,stkwl
c
c     MAXST = maximum number of stations used (plus 2)
c     MAXD = maximum datalength
c     npts = number of data points per trace
c     zssl = linear trace stack
c     zscp = quadratic trace stack
c     sampin = sample interval
c     zu = Array containing all normalized data
c     dtcs = local time shift from stacking iteration
c     dtmo = initial (model) moveout correction
c     nstkwb = number of samples to start of stack window
c     nstkwl = number of samples in stack window
c     tshify = Time shift for stacked traces
c     nsta = number of stations
c     stkwb,stkwl = start,length of stack window
c     loff = reference sample of stacked trace
c     nstk = maximum number of samples in stacked trace
c     lim1,lim2 = Limits of stack window
c
c----------------------------------------------------------------------
c     OTHER VARIABLES
c----------------------------------------------------------------------
c
      integer i,j,l,m,j1,j2
      integer jim1,jim2,imo,lu,js,jm,jmi
      integer nsi,ist,icl,swl,swr,w(MAXST)
      integer iyr,imon,iday,ihr,imin,nstb,npmx
      real ttto,den,err(MAXST),wsp(5000),emax,emin
      real evlat,evlon,evdep,sec,scu,pstakn
      character*30 infile
c      character aqini*20,aqfini*20,akph*8,ttres*20
      character akph*8
      character*8 csta(MAXST)
      character(len=256) :: directory, efile, aqini, aqfini, ttres
      real vu(2000),erl,errl,errr,zv(MAXST,MAXD),tshift(MAXST)
      real dtcmin,dtcmax,pjgl,wsl,dtmx
      real wm,ws,tshify
c
c
c     nsi = Number of stacking iterations
c     ttto = Time to trace origin
c     efile = Event file name
c     vu = Data buffer
c     aqini = Name of initial stack aq file
c     aqfini = Name of final stack aq file
c     akph = AK135 phase name for alignment
c     ttres = Name of traveltime residual output
c     ist,icl = counters for error determination
c     err = pick error estimated from trace power
c     wsp = power of weighted stack
c     erl = error factor
c     swl,swr = switches for half width error determination
c     errl,errr = Half width errors (left and right)
c     iyr,imon,iday = Year,month,day of trace start time
c     ihr,imin,sec = Hour,minute,second of trace start time
c     evlat,evlon,evdep = event latitude, longitude, depth
c     zv = Array containing all data in raw form
c     tshift = time shift to align each trace
c     w = weight applied to each trace (usually 0 or 1)
c     dtcmin,dtcmax = min,max bounds on differential time search
c     pjgl = Stack index (the Lp norm used)
c     csta = station name
c     scu = denominator for trace normalization
c     nstb=nsta+2
c     wsl = weight applied to stacked trace
c     tshify = Time shift for stacked traces
c     pstakn = L2 measure of trace misfit
c     dtmx = maximum model predicted time shift
c     npmx = trace with maximum number of samples
c     emax,emin = Maximum and minimum errors (ms)
c
c----------------------------------------------------------------------
c
c     Open file containing input parameters
c
      open(RUN,file='tcas.cmd',status='old')
c
c     Read in number of stacking iterations
c
      read(RUN,*)nsi
c
c     Read in stack index   
c
      read(RUN,*) pjgl
c
c     Read in error factor
c
      read(RUN,*)erl
c
c     Read in error limits
c
      read(RUN,*)emin,emax
      emin=emin/1000.0
      emax=emax/1000.0
c
c     Read in stack window   
c
      read(RUN,*) stkwb,stkwl
c
c     Read in name of event file
c
      read(RUN,fmt='(A)') infile
      print*, infile


c
c     Read in bounds on differential time search
c
      read(RUN,*) dtcmin
      read(RUN,*) dtcmax
      read(RUN,fmt='(A)') directory
      efile=trim(directory)//"/Input_data/"//trim(infile)
      print*, trim(efile)
      open(UNIT=EVT, FILE=trim(efile), STATUS='old')
      rewind(EVT)
      close(RUN)
c
c     Construct initial and final stack output
c     filenames.
c
      aqini=trim(directory)//"/Output_data/asi"//infile(1:24)//".aq"
      aqfini=trim(directory)//"/Output_data/asf"//infile(1:24)//".aq"
c
c     Construct file name for traveltime
c     residual output.
c
      ttres=trim(directory)//"/Output_data/ts"//infile(1:24)//".ttr"

      print*, trim(aqini)
      print*, trim(aqfini)
      print*, trim(ttres)
c
c----------------------------------------------------------------------
c     OBTAIN EVENT  DATA
c----------------------------------------------------------------------
c
c     Read in AQ formatted event data
c
      read(EVT,*) nsta
      read(EVT,*) evlat, evlon, evdep
      read(EVT,*) iyr,imon,iday
      read(EVT,*) ihr,imin,sec,ttto
      read(EVT,*) sampin,akph
      DO i=1,nsta
         read(EVT,*)  w(i),npts(i),tshift(i),csta(i)
         read(EVT,*)  (vu(j),j=1,npts(i))
         scu = 0.
         err(i)=0.0
         DO j=1,npts(i)
            zu(i,j) = w(i)*vu(j)
            if(abs(zu(i,j)).gt.scu) scu = abs(zu(i,j))
            zv(i,j) = vu(j)
         ENDDO
c
c        Normalise working copy zu
c
         IF(scu.gt.0.0) THEN
            DO j=1,npts(i)
               zu(i,j) = zu(i,j)/scu
            ENDDO
         ENDIF
      ENDDO
      CLOSE(EVT)
      nstkwb = int(stkwb/sampin)
      nstkwl = int(stkwl/sampin)
c                    
c     Apply moveout correction (probably from AK135)
c
      DO i =1,nsta
        dtmo(i) = -tshift(i)
        dtcs(i) = 0.
      ENDDO
c
c     Calculate array bounds for stacked trace
c
      dtmx=0.
      npmx=0
      DO l=1,nsta
         if(abs(dtmo(l)).gt.dtmx)dtmx=abs(dtmo(l))
         if(npts(l).gt.npmx)npmx=npts(l)
      ENDDO
      loff = abs(dtcmin)
      if(abs(dtcmax).gt.abs(dtcmin))loff=abs(dtcmax)
      loff = nint((dtmx+loff)/sampin)+1
      nstk=2*loff+npmx
      lim1 = loff+nstkwb
      lim2 = loff+nstkwb+nstkwl
c
c     Stack all preliminarily aligned traces
c
      call pstack(pstakn)
c
c     Write preliminarily aligned data to file
c
      nstb = nsta+2
      OPEN(INI,file=trim(aqini),status='unknown')
      write(INI,*) nstb
      write(INI,*) evlat, evlon, evdep
      write(INI,*) iyr,imon,iday
      write(INI,*) ihr,imin,sec
      write(INI,*) sampin
      DO i=1,nsta
        write(INI,*)  w(i),npts(i),-dtmo(i),"  ",csta(i)
        write(INI,*)  (zv(i,j),j=1,npts(i))
      ENDDO  
      wsl = 1.
      j1 = loff
      j2 = loff+nstk-1
      tshify = 0.
      write(INI,*) wsl,nstk,tshify," zssl"
      write(INI,*) (zssl(j),j=j1,j2)
      write(INI,*) wsl,nstk,tshify," zscp"
      write(INI,*) (zscp(j),j=j1,j2)      
      CLOSE(INI)    
c----------------------------------------------------------------------
c     Start the adaptive stacking procedure
c----------------------------------------------------------------------
      jim1 = nint(dtcmin/sampin)
      jim2 = nint(dtcmax/sampin) 
      DO m = 1,nsi
         DO i = 1,nsta  
            imo = nint(dtmo(i)/sampin)      
            IF(w(i).eq.0) THEN
               dtcs(i) = 0.
               GOTO 19
            ENDIF
            wm = 1.e6
            DO js = jim1,jim2
               ws = 0.
               DO l = lim1,lim2
                  lu = l-loff+imo+js
                  ws = ws + abs(zssl(l)-zu(i,lu))**pjgl
               ENDDO
               ws = ws/stkwl
               wsp(js+1-jim1)=ws
               IF(ws.lt.wm) THEN
                  wm = ws
                  jm = js
               ENDIF
            ENDDO 
c
c           Below, we estimate error from the power
c           of the stacked trace.
c
            if(m.eq.nsi)then
c
c              Determine the width of the stack
c              and tag the location of the minimum
c
               ist=jim2-jim1+1
               jmi=jm+1-jim1
c
c              Calculate the cross-over point for <jmi
c            
               if(jmi.eq.1)then
                  swl=2
               else
                  l=jmi-1
                  swl=0
                  do while(swl.eq.0)
                     if(wsp(l).ge.wsp(jmi)*erl)then
                        icl=l
                        swl=1
                     else
                        l=l-1
                        if(l.eq.0)swl=2
                     endif
                  enddo
                  if(swl.eq.1)then
                     den=wsp(icl+1)-wsp(icl)
                     if(abs(den).gt.1.0e-5)then
                        errl=sampin*(wsp(jmi)*erl-wsp(icl))/den
                        errl=(jmi-icl)*sampin-errl
                     else
                        errl=(jmi-icl)*sampin
                     endif
                  endif
               endif
c
c              Calculate cross-over point for >jmi
c
               if(jmi.eq.ist)then
                  swr=2
               else
                  l=jmi+1
                  swr=0
                  do while(swr.eq.0)
                     if(wsp(l).ge.wsp(jmi)*erl)then
                        icl=l
                        swr=1
                     else
                        l=l+1
                        if(l.eq.ist)swr=2
                     endif
                  enddo
                  if(swr.eq.1)then
                     den=wsp(icl)-wsp(icl-1)
                     if(abs(den).gt.1.0e-5)then
                        errr=sampin*(wsp(jmi)*erl-wsp(icl-1))/den
                        errr=(icl-jmi-1)*sampin+errr
                     else
                        errr=(icl-jmi)*sampin
                     endif
                  endif
               endif
c
c              Take average of errr and errl for actual error
c
               if(swr.eq.1.and.swl.eq.1)then
                  err(i)=(errr+errl)/2.0
               else if(swl.eq.1)then
                  err(i)=errl
               else if(swr.eq.1)then
                  err(i)=errr
               else
                  err(i)=emax
               endif
c
c              Constrain error limits
c
               if(err(i).lt.emin)then
                  err(i)=emin
               else if(err(i).gt.emax)then
                  err(i)=emax
               endif
            endif
c
c           Determine the time shift
c
            dtcs(i) = -float(jm)*sampin        
 19         CONTINUE
         ENDDO
         call pstack(pstakn)
         write(6,*) "pstakn = ", pstakn
      ENDDO
c
c     Write residuals and associated errors to file
c
      open(RES,file=trim(ttres),status='unknown')
      write(RES,*)nsta
      write(RES,'(f9.4,2X,f9.4,2X,f8.3)')evlat,evlon,evdep
      write(RES,'(i6,i5,i5)')iyr,imon,iday
      write(RES,'(i5,i5,2X,f7.4)')ihr,imin,sec
      write(RES,'(f9.4)')ttto
      write(RES,'(f7.4)')sampin
      write(RES,*)akph
      write(RES,*)nsi
      DO i=1,nsta
         write(RES,'(i4,2X,a5,2X,f7.4,2X,f7.4,2X,f8.4,i4)') 
     *         i,csta(i),dtcs(i),err(i),tshift(i),w(i)
        tshift(i) = -dtmo(i)+dtcs(i)
      ENDDO
      close(RES)
c
c     Write final trace stack to file
c
      nstb = nsta+2
      open(FIN,file=trim(aqfini),status='unknown')
      write(FIN,*) nstb
      write(FIN,*) evlat, evlon, evdep
      write(FIN,*) iyr,imon,iday
      write(FIN,*) ihr,imin,sec
      write(FIN,*) sampin
      DO i=1,nsta
         write(FIN,*)  w(i),npts(i),tshift(i),"  ",csta(i)
         write(FIN,*)  (zv(i,j),j=1,npts(i))
      ENDDO  
      wsl = 1.
      j1 = loff
      j2 = loff+nstk-1
      tshify = 0.
      write(FIN,*) wsl,nstk,tshify,"  zssl"
      write(FIN,*) (zssl(j),j=j1,j2)
      write(FIN,*) wsl,nstk,tshify,"  zscp"
      write(FIN,*) (zscp(j),j=j1,j2)        
      CLOSE(FIN)   
      stop
      end

*********************************************************************
      subroutine pstack(pstakn)	
      implicit none
*********************************************************************
c
c     Perform a linear and quadratic stack of all traces
c
c--------------------------------------------------------------------
c     COMMON BLOCK DECLARATIONS
c--------------------------------------------------------------------
c
      integer MAXST,MAXD
      parameter (MAXST=100,MAXD=5000)
      integer npts,nstkwb,nstkwl,nsta,loff,nstk
      integer lim1,lim2
      real zssl,zscp,sampin
      real zu,dtcs,dtmo,stkwb,stkwl
c
      common /IST1/ npts(MAXD),nstkwb,nstkwl,nsta,loff,nstk
      common /IST2/ lim1,lim2
      common /RST1/ zu(MAXST,MAXD),zssl(MAXD),zscp(MAXD)
      common /RST2/ dtcs(MAXST),dtmo(MAXST),sampin,stkwb,stkwl
c
c     MAXST = maximum number of stations used (plus 2)
c     MAXD = maximum datalength
c     npts = number of data points per trace
c     zssl = linear trace stack
c     zscp = quadratic trace stack
c     sampin = sample interval
c     zu = Array containing all normalized data
c     dtcs = local time shift from stacking iteration
c     dtmo = initial (model) moveout correction
c     nstkwb = number of samples to start of stack window
c     nstkwl = number of samples in stack window
c     tshify = Time shift for stacked traces
c     nsta = number of stations
c     stkwb,stkwl = start,length of stack window
c     loff = reference sample of stacked trace
c     nstk = maximum number of samples in stacked trace
c     lim1,lim2 = Limits of stack window
c
c---------------------------------------------------------------------
c     OTHER VARIABLES
c---------------------------------------------------------------------
c
      integer i,l,lmn
      real zcu(MAXD),scz,pstakn
c     pstakn = L2 measure of trace misfit
c---------------------------------------------------------------------
      DO l=1,nstk
         zssl(l) = 0.
         zscp(l) = 0.
      ENDDO
c
c     loop on stations
c
      DO i = 1,nsta
         DO l=1,nstk
            zcu(l) = 0.
         ENDDO
         lmn = nint((-dtmo(i)+dtcs(i))/sampin)
         DO l=1,npts(i)
            zcu(loff+lmn+l) = zu(i,l)
         ENDDO    
c
c        linear and quadratic stack
c
         DO l=lim1,lim2
            zssl(l) = zssl(l) + zcu(l)
            zscp(l) = zscp(l) + zcu(l)*zcu(l)
         ENDDO
      ENDDO
      pstakn = 0.
      DO l=lim1,lim2
         pstakn = pstakn + abs(zscp(l))
      ENDDO
      pstakn = pstakn/(real(nsta)*stkwl)
      scz = 0.
      DO l=lim1,lim2
         IF(abs(zssl(l)).gt.scz) scz = abs(zssl(l))
      ENDDO
      IF(scz.gt.0.0) THEN
         DO l=lim1,lim2
            zssl(l) = zssl(l)/scz
         ENDDO
      ENDIF
      return
      end

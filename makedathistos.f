	subroutine makedathistos(runno,dummyflag,posiflag,counts,charge,charge_elcl,b4cflag, ca48flag, doing_shms, doing_hms)

	implicit none
	save

	include 'hbook.inc'
	INCLUDE 'ntuple.cmn'
	include 'cuts.inc'
	real*8 counts,charge,charge_elcl
	integer*4	lrecl, istat
	integer*4	runno

	character*80 infile
	logical dummyflag,posiflag,b4cflag,ca48flag,doing_shms,doing_hms,do_skim

	integer GetNumBranches
	external GetNumBranches


	do_skim=.false.
	write(33,*) '------------- Analyzing run number ',runno,' -------------'
	if(do_skim) then
	   write(33,*) '-------------------- Using skim files --------------------'
	endif

	lrecl = 4096
	if(do_skim) then
	   if(doing_shms) then
	      write(infile,'("skimfiles/shms_",i4,"_skim.txt")') runno
	   endif

	   if(doing_hms) then
	      write(infile,'("skimfiles/hms_",i4,"_skim.txt")') runno
	   endif
	else ! using root trees
	   if(doing_shms) then
	      write(infile,'("rootfiles/shms_replay_production_",i4,"_-1.root")') runno
	   endif

	   if(doing_hms) then
	      write(infile,'("rootfiles/hms_replay_production_",i4,"_-1.root")') runno
	   endif
	endif

	write(33,*) infile

	if(do_skim) then
	   open(unit=1,file=infile,status='old')
	else
	   call InitRootNT(infile,'READ');
	   print *, GetNumBranches(),' branches'
	endif

	call readdat(runno,dummyflag,posiflag,counts,charge,charge_elcl, b4cflag, ca48flag, doing_shms,doing_hms,do_skim)

	call hrend ('ZXC')
	close (1)
	close (10)
	end



	subroutine readdat(runno,dummyflag,posiflag,counts,charge,charge_elcl,b4cflag, ca48flag, doing_shms,doing_hms,do_skim)

	implicit none

	INCLUDE 'ntuple.cmn'
	include 'hbook.inc'
	include 'cuts.inc'
	include 'kinematics.cmn'
c	include 'hcer_eff.cmn'
	include 'extcor.cmn'
	include 'cryocor.cmn'

	integer*4 ierr, nevt, ievt, runno,neloss,istat
c	integer*4 i
	integer runcount,runid,kind(32),loc
	integer usflag
c	integer nvar
c	real*8 rhigh(49),rlow(49)
	real*8 hsxfp,hsyfp,hsxpfp,hsypfp
	real*8 hsytar,hsxptar,hsyptar,hszbeam,hsdelta
	real*8 hsztar,ntrack,gtrindex
	real*8 hcer_npe,hsprtrk,hsshtrk,hsshsum,hsbeta,hsdedx
	real*8 hsphi,hsp,nu,y_scale
	real*8 dy_tmp,cor
	real*8 x_bj_ntup,Q2_ntup,w_ntup,xi_ntup,hstheta_ntup
	real*8 x_bj,Q2,xi,hstheta,w
	real*8 eventID,ev_type
	real*8 hselreal,hselclean,hspipre,hselhi,hsello
	real*8 hsprhi,hsprlo,hsshlo,hstart,ntup_charge
	real*8 gfrx_raw,gfry_raw,gfrx,gfry,bcm1bcm2,bcm4a,bcm4c
C New pass 3 variables
	real*8 hsshtrkp,hsshtrkn,hsshsuma,hsshsumb,hsshsumc,hsshsumd
	real*8 hsscin,hsstof
C cer eff variables
	real*8 hxcer,hcereff,hmscereff_xem2,xlo,xhi
	real*8 posicor,pcor1,pcor2
	real*8 pioncor,picor1,picor2,picor3,picor4,picor5,picor6
	real*8 eprime,eloss_tmp, w2e
	real*8 dxp,ddel
	real*8 xifunc,pmin,pmax
	real*8 xmin,xmax,xbjfunc
	real*8 temp,ext_weight
	real*8 ext_func
	real*8 c0_new, c0_old
	real*8 c0_yp
	real*8 dzero
	real*4 zero
	real*8 calcor,cercor,deltatmp,deltacor,ajicer
	real*8 thetamin,thetamax,epmin,epmax,deltamin,deltamax
	real*8 c1,c2
	real*8 xcer,ycer,cer_eff, caleff
	real*8 boilslope,cur_flag,index
	real*8 ytarweight
	real*8 abishekcereff,shmscereff,hmscereff,ybender,xbender
	real*8 he3_boilslope,current
	logical shms_dipole_exit, inside_dipole
	


	
	real*8 charge,charge_elcl

	integer*4 npass
	real*8 counts
	integer idbase
	integer hbookc(51)
	common/hcbook/hbookc

	integer GetNtEntries
	external GetNtEntries 


	logical dummyflag,posiflag,upflag,first,b4cflag,ca48flag,doing_shms,doing_hms,do_skim
	logical loop

	data runcount /0/
	data first /.true./
	save   !remember it all

	runcount=runcount+1
	runid = runcount*100


	zero =0.0
	dzero=0.0
c	write(6,*) 'runcount and runid',runcount,runid

C Read in kinematics.
	call getkine(runno,dummyflag,posiflag,charge,b4cflag,ca48flag,current,boilslope,doing_shms,doing_hms)

	write(6,*) 'targ_A', targ_A
	if(dummyflag) then
	   write(6,*) 'dummy run: B4C or Ca48?',b4cflag,ca48flag
	   if(b4cflag) then
	      if(targ_A.eq.10) then !B4C
		 charge = charge*4.319654
	      elseif(targ_A.eq.11) then !11B4C
		 charge = charge*4.231550
	      else
		 write(6,*) 'Not a B4C target - Danger Will Robinson'
c		 stop
	      endif
	   elseif(ca48flag) then
	      charge = charge*7.49901
	   else
	      if(targ_A.eq.1.0 .or. targ_A.eq.4.0) then !loop 2
		 charge = charge*3.5889
	      elseif (targ_A.eq.2.0 .or. targ_A.eq.3.0) then !loop 1
		 charge = charge*4.5733
	      else
		 write(6,*) 'Doing dummy on a solid target?'
		 write(6,*) 'Screw you guys - Im going home.'
		 stop
	      endif
	   endif
	endif
C Initialize some arrays and counters.
	npass=0


c read in matrix elements
c	call h_targ_trans_init

	
	if(.not.do_skim) then
	   nevt = GetNtEntries()
	   write(6,*) 'number of entries: ',nevt
	   if(.not.posiflag) then
	      if(nevt.eq.0) then
		 write(6,*) 'No entries in tree for this run - stopping analysis'
		 stop
	      endif
	   endif
	   call TagNtBranch(hsdelta,'H.gtr.dp ') 
	   call TagNtBranch(hsxptar,'H.gtr.th ') 
	   call TagNtBranch(hsyptar,'H.gtr.ph ') 
	   call TagNtBranch(hsytar,'H.gtr.y ') 
	   call TagNtBranch(hsxfp,'H.dc.x_fp ') 
	   call TagNtBranch(hsyfp,'H.dc.y_fp ') 
	   call TagNtBranch(hsxpfp,'H.dc.xp_fp ') 
	   call TagNtBranch(hsypfp,'H.dc.yp_fp ') 
	   call TagNtBranch(hcer_npe,'H.cer.npeSum ') 
	   call TagNtBranch(hsshsum,'H.cal.etottracknorm ') 
c	   call TagNtBranch(hsbeta,'H.gtr.beta ') 
c	   call TagNtBranch(ntrack,'H.dc.ntrack ') 
c	   call TagNtBranch(gtrindex,'H.gtr.index ') 
c	   call TagNtBranch(gfry_raw,'H.rb.raster.fryaRawAdc') ! not even calibrated raster in tree
c	   call TagNtBranch(bcm4c,'H.bcm.bcm4c.AvgCurrent ') 
cc	   call TagNtBranch(bcm4a,'H.bcm.bcm4a.AvgCurrent ') 
	   call TagNtBranch(cur_flag,'H.bcm.CurrentFlag ') 
	   call TagNtBranch(index, 'H.gtr.index ')
c	   call TagNtBranch(hselhi, 'T.shms.pEL_HI_tdcTime ')
c	   call TagNtBranch(hselclean, 'T.shms.pEL_CLEAN_tdcTime ')
	endif


	ievt =0
	loop=.true.
        do while (loop)
	   if(do_skim) then
	      if((abs(th_deg-13.0).lt.1.)) then
		 read(1,*,end=99) hsdelta,hsxptar,hsyptar,hsytar,hsxfp,hsyfp,hsxpfp,hsypfp,
     >            hcer_npe,hsshsum,hsshtrk,hsbeta
	      else
		 read(1,*,end=99) hsdelta,hsxptar,hsyptar,hsytar,hsxfp,hsyfp,hsxpfp,hsypfp,
     >            hcer_npe,hsshsum,hsshtrk,hsbeta,bcm4a,bcm4c
c		 read(1,*,end=99) hsdelta,hsxptar,hsyptar,hsytar,hsxfp,hsyfp,hsxpfp,hsypfp,
c     >            hcer_npe,hsshsum,hsshtrk,hsbeta
	      endif
	   else
	      if(ievt.gt.nevt-2) then
		 loop=.false.
	      endif
	      call ReadNtBranch(ievt)
	   endif
	   ievt = ievt+1

c	   if(targ_A.ge.1 .and. targ_A.le.4) then
c	      cor = 0.991117*hsdelta-0.00190166*hsdelta**2+0.000295015*hsdelta**3+
c     >              2.4572E-6*hsdelta**4-2.68198E-6*hsdelta**5+
c     >              1.77248E-7*hsdelta**6-hsdelta
c	      hsdelta = hsdelta-cor
c	   else
c	      cor = 0.991333*hsdelta-0.0026865*hsdelta**2+0.000233225*hsdelta**3+
c     >              0.00001006*hsdelta**4-1.96515E-6*hsdelta**5+
c     >              8.42355E-8*hsdelta**6-hsdelta
c	      hsdelta = hsdelta-cor
c	   endif

c	   call h_targ_trans(hsxfp,hsxpfp,hsyfp,hsypfp,gfry,hsdelta,
c     >           hsxptar,hsyptar,hsytar,istat)

C PID cuts
	  if (hcer_npe .le. 2.0) goto 888
	  if(hsshsum .le. 0.7) goto 888
C Now do cuts on reconstructed quantities 
C HMS delta cut
	  if (hsdelta.ge.hdeltacuthi) goto 888
	  if (hsdelta.le.hdeltacutlo) goto 888
	  if (isnan(hsdelta)) goto 888
c	  if(abs(hsytar).gt.200.0) goto 888
c	  if(isnan(hsytar)) goto 888
	  if (abs(hsxptar).ge.0.085)goto 888
	  if (abs(hsyptar).ge.0.032)goto 888
	  if(cur_flag.ne.1) goto 888
cc	  if(index.le.-1.0) goto 888


	  hsztar=hsytar/sin(th_rad)

	  if(ytar.gt.0.5) then
	     if(doing_shms) then
		if(hsytar.lt.0.115) goto 888
	     else
		if(hsytar.gt.0.0) goto 888
	     endif
	  endif

	  if(ytar.lt.-0.5) then
	     if(doing_shms) then
		if(hsytar.gt.0.115) goto 888
	     else
		if(hsytar.lt.0.0) goto 888
	     endif
	  endif


C use eloss extracted above in physics recon.
	  call get_physics(dzero,hsdelta,hsxptar,hsyptar,W,Q2,x_bj,xi,hstheta,eprime,doing_shms,doing_hms)

c	  if(targ_A.lt.1.5) then ! W cut for hydrogen
cdg	     if(w**2.lt.2.0) goto 888
c	  endif

	  if(targ_A.gt.2 .and. targ_A.lt.4) then
	     boilslope=abs(he3_boilslope(hsztar))/100.0/100.0
	  endif
	  currentcor = 1.0/(1.0-current*boilslope)
	  if(dummyflag) then
	     currentcor =1.0
	  endif
	  caleff=1.0
	  if(doing_shms) then
	     caleff = -0.000005759 + 0.9984 - exp(-1.98 - 0.2356 * eprime**3)
	  endif

	  npass=npass+1
C If all cuts are passed, fill histos
	  if(dummyflag) then
	     if(posiflag) then
		idbase = 9
	     else
		idbase = 6
	     endif
	  else 
	     if(posiflag) then
		idbase = 4
	     else
		idbase = 1
	     endif
	  endif


	  ytarweight=1.0
c	  if(dummyflag) then !extra correction for difference in ext. radcor
c	     if(b4cflag) then
c		ytarweight=1.0
c	     else	      
c		usflag=0
c		if(targ_A.eq.1.0) then !loop 2
c		   if(hsytar>0) then
c		      if(doing_shms) then
c			 usflag=1
c			 ytarweight =3.789357/4.30842
c		      else
c			 ytarweight=3.789357/3.381714
c		      endif
c		   else
c		      if(doing_hms) then
c			 usflag=1
c 			 ytarweight=3.789357/4.30842
c		      else			 
c			 ytarweight=3.789357/3.381714
c		      endif
c		   endif
c		elseif(targ_A.eq.2.0) then ! loop 3
c		   if(hsytar>0) then
c		      if(doing_shms) then
c			 usflag=1
c			 ytarweight =4.06343/4.971257
c		      else
c			 ytarweight=4.06343/3.43567
c		      endif
c		   else
c		      if(doing_hms) then
c			 usflag=1
c			 ytarweight =4.06343/4.971257
c		      else
c			 ytarweight=4.06343/3.43567
c		      endif
c		   endif
c		endif
c		if(th_deg.gt.14.0) then
c		   ytarweight = ytarweight*ext_func(eprime,hstheta,usflag)
c		endif
c		write(6,*) 'ytarweight',ytarweight,ext_func(eprime,hstheta,usflag)
c	     endif
c	  endif


c get cerenkov efficiency
	  cer_eff=hmscereff_xem2(hsdelta)
c	  cer_eff=1.0
	  
c add positron subtraction correction 
	  posicor=1.0
	  pcor1=-1.0
	  pcor2=-1.0
	  if(doing_hms) then
	     if(abs(thdegcent-26.0)<0.1) then
		if(targ_A.eq.1.0) then
		   pcor1=3.200
		   pcor2=-3.067
c		   if(dummyflag) then
c		      pcor1=1.7244
c		      pcor2=-2.1505
c		   endif
		elseif(targ_A.eq.2.0) then
		   pcor1=3.238
		   pcor2=-3.018
c		   if(dummyflag) then
c		      pcor1=1.7244
c		      pcor2=-2.1505
c		   endif
		elseif(targ_A.eq.12.0) then
		   pcor1=2.2685
		   pcor2=-2.6566
		elseif(targ_A.eq.40.0) then
		   pcor1=2.3573
		   pcor2=-2.3468
		endif
		if(eprime.lt.3.0745) then
		   posicor = 1.0-exp(pcor1+eprime*pcor2)
		else
		   posicor=1.0
		endif
	     endif
	     if(abs(thdegcent-20.0)<0.1) then
		if(targ_A.eq.1.0) then
		   pcor1=2.8300
		   pcor2=-2.4893
		   if(dummyflag) then
		      pcor1=1.7244
		      pcor2=-2.1505
		   endif
		elseif(targ_A.eq.2.0) then
		   pcor1=2.8778
		   pcor2=-2.4620
		   if(dummyflag) then
		      pcor1=1.7244
		      pcor2=-2.1505
		   endif
		elseif(targ_A.eq.3.0) then
		   pcor1=2.4665
		   pcor2=-2.4223
		   if(dummyflag) then
		      pcor1=1.7244
		      pcor2=-2.1505
		   endif
		elseif(targ_A.eq.4.0) then
		   pcor1=2.3352
		   pcor2=-2.2914
		   if(dummyflag) then
		      pcor1=1.7244
		      pcor2=-2.1505
		   endif
		elseif(targ_A.eq.6.0) then
		   pcor1=1.7863
		   pcor2=-2.2679
		elseif(targ_A.eq.7.0) then
		   pcor1=1.9246
		   pcor2=-2.3439
		elseif(targ_A.eq.9.0) then
		   pcor1=2.0396
		   pcor2=-2.2490
		elseif(targ_A.eq.10.0) then
		   pcor1=1.9480
		   pcor2=-2.2513
		elseif(targ_A.eq.11.0) then
		   pcor1=2.0291
		   pcor2=-2.2794
		elseif(targ_A.eq.12.0) then
		   pcor1=1.9864
		   pcor2=-2.2508
		elseif(targ_A.eq.27.0) then
		   pcor1=1.7341
		   pcor2=-2.1419
		elseif(targ_A.eq.40.0) then ! calcium 40
		   pcor1=2.2831
		   pcor2=-2.1151
		elseif(targ_A.eq.48.0) then ! calcium 48
		   pcor1=2.3785
		   pcor2=-2.1161
		elseif(targ_A.eq.47.9183) then ! titanium
		   pcor1=1.7862
		   pcor2=-2.1606
		elseif(targ_A.eq.54.0) then ! Fe 54
		   pcor1=1.9325
		   pcor2=-2.1383
		elseif(targ_A.eq.58.0) then ! nickel 58
		   pcor1=1.5876
		   pcor2=-2.0857
		elseif(targ_A.eq.64.0) then !nickel 64
		   pcor1=1.5647
		   pcor2=-2.0778
		elseif(targ_A.eq.64.234) then !copper
		   pcor1=2.1489
		   pcor2=-1.9655
		elseif(targ_A.eq.107.96322) then !silver
		   pcor1=2.2068
		   pcor2=-2.0492
		elseif(targ_A.eq.118.7705) then !tin
		   pcor1=1.2090
		   pcor2=-1.6651
		elseif(targ_A.eq.197.0) then !gold
		   pcor1=1.7988
		   pcor2=-1.9448
		elseif(targ_A.eq.232.0) then ! thorium
		   pcor1=2.1016
		   pcor2=-2.0198
		endif
		if(eprime.lt.3.65) then
		   posicor = 1.0-exp(pcor1+eprime*pcor2)
		else
		   posicor=1.0
		endif
	     endif
	  endif
	  if(pcor1.lt.0) then
	     posicor=1.0
	  endif
c add pion subtraction correction if not using cerenkov- this only works at 21 degrees
	  pioncor=1.0
c	  posicor=1.0


	  ext_weight = ytarweight*denscor*currentcor*posicor*pioncor/(cer_eff*caleff)
c	  if(ievt.lt.1000) then
c	     write(6,*) 'Total weight: ', ext_weight, counts
c	     write(6,*) 'ytarweight, boilslope',ytarweight, boilslope
c	     write(6,*) 'Denscor, posicor, pioncor: ', denscor,posicor,pioncor, currentcor
c	     write(6,*) 'cer_eff, caleff: ', cer_eff, caleff
c	  endif
	  counts = counts+1.0*ext_weight



	  call hfill(3000+idbase,sngl(hsdelta),zero,sngl(ext_weight))
	  call hfill(3100+idbase,sngl(1000.0*hsxptar),zero,sngl(ext_weight))
	  call hfill(3200+idbase,sngl(1000.0*hsyptar),zero,sngl(ext_weight))
	  call hfill(3300+idbase,sngl(hsytar),zero,sngl(ext_weight))
	  call hfill(3400+idbase,sngl(hsxfp),zero,sngl(ext_weight))
	  call hfill(3500+idbase,sngl(hsxpfp),zero,sngl(ext_weight))
	  call hfill(3600+idbase,sngl(hsyfp),zero,sngl(ext_weight))
	  call hfill(3700+idbase,sngl(hsypfp),zero,sngl(ext_weight))
	  call hfill(3800+idbase,sngl(W),zero,sngl(ext_weight))
	  call hfill(3900+idbase,sngl(Q2),zero,sngl(ext_weight))
	  call hfill(4200+idbase,sngl(hstheta-th_rad),zero,sngl(ext_weight))
	  call hfill(4300+idbase,sngl(eprime),zero,sngl(ext_weight))
	  call hfill(4400+idbase,sngl(hsdelta),sngl(hstheta-th_rad),sngl(ext_weight))
	  call hfill(4100+idbase,sngl(xi),zero,sngl(ext_weight))
	  call hfill(4000+idbase,sngl(x_bj),zero,sngl(ext_weight))
	  

	  if( (.not.posiflag).and.(hselclean.gt.0.0)) then !fill "elclean" histos
	     idbase = idbase+2
	     call hfill(3000+idbase,sngl(hsdelta),zero,sngl(ext_weight))
	     call hfill(3100+idbase,sngl(hsxptar),zero,sngl(ext_weight))
	     call hfill(3200+idbase,sngl(hsyptar),zero,sngl(ext_weight))
	     call hfill(3300+idbase,sngl(hsytar),zero,sngl(ext_weight))
	     call hfill(3400+idbase,sngl(hsxfp),zero,sngl(ext_weight))
	     call hfill(3500+idbase,sngl(hsxpfp),zero,sngl(ext_weight))
	     call hfill(3600+idbase,sngl(hsyfp),zero,sngl(ext_weight))
	     call hfill(3700+idbase,sngl(hsypfp),zero,sngl(ext_weight))
	     call hfill(3800+idbase,sngl(W),zero,sngl(ext_weight))
	     call hfill(3900+idbase,sngl(Q2),zero,sngl(ext_weight))
	     call hfill(4200+idbase,sngl(hstheta-th_rad),zero,sngl(ext_weight))
	     call hfill(4300+idbase,sngl(eprime),zero,sngl(ext_weight))
	     call hfill(4400+idbase,sngl(hsdelta),sngl(hstheta-th_rad),sngl(ext_weight))
	     call hfill(4100+idbase,sngl(xi),zero,sngl(ext_weight))
	     call hfill(4000+idbase,sngl(x_bj),zero,sngl(ext_weight))
	  endif

	  
 888	  continue

	enddo

	if(.not.do_skim) then
	   call RootNTOutp();
	endif

 99	continue

	write(33,*) 'Number of events passing cuts is',npass
	write(6,*) 'Number of events passing cuts is',npass
	write(33,*) 'Number of weighted events passing cuts is',counts
	write(33,*) 'fractional error= (%)',100.0*sqrt(dble(npass))/npass
	if(posiflag) then
	   write(33,*) 'positron counts per charge=',counts/charge
	   write(33,*) 'yield uncertainty=',sqrt(counts)/charge
	else
	   write(33,*) 'electron counts per charge=',counts/charge
	   write(33,*) 'yield uncertainty=',sqrt(counts)/charge
	endif



	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	real*8 function xifunc(eb,ep,th)

	real*8 eb,ep,th
	real*8 x,nu,Q2,mp

	parameter(mp=0.93827231)

	nu = eb-ep
	Q2=4.*eb*ep*sin(th/2.)**2
	x=Q2/2./mp/nu

	xifunc = 2.*x/(1.0+sqrt(1.+4.*mp**2*x**2/Q2))

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	real*8 function xbjfunc(eb,ep,th)

	real*8 eb,ep,th
	real*8 x,nu,Q2,mp

	parameter(mp=0.93827231)

	nu = eb-ep
	Q2=4.*eb*ep*sin(th/2.)**2
	xbjfunc=Q2/2./mp/nu

	return
	end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	real*8 function ext_func(ep,theta,usflag)
	
	implicit none

	real*8 ep,theta,thdeg
	real*8 thhi,thlo,delta_th
	real*8 ephi,eplo,delta_ep
	real*8 A1,A2,A3,A4,A12,A34,A1234

	integer thcount,epcount,usflag
	logical upflag

	include 'extcor.cmn'


	thdeg = theta*180.0/3.141592654

	if(thdeg.lt.thext(1).or.thdeg.gt.thext(50)) then
	   ext_func=1.0
	   write(6,*) 'Warning: theta beyond external extrat table limits!',thdeg,thext(1),thext(50)
	   write(6,*) 'using thmin/thmax for now...'
	   if(thdeg.gt.thext(50)) thdeg=thext(50)
	   if(thdeg.lt.thext(1)) thdeg=thext(1)
c	   return
	endif

	if(ep.lt.epext(1).or.ep.gt.epext(npmax)) then
	   ext_func=1.0
	   write(6,*) 'Warning: eprime beyond radiation table limits! (in ext_func)',ep
c	   return
	endif

	do thcount=1,50-1
	   thhi=0.0
	   thlo=0.0
	   ephi=0.0
	   eplo=0.0
	   if( (thdeg.gt.thext(thcount)) .and. (thdeg.le.thext(thcount+1)) ) then
	      thhi=thext(thcount+1)
	      thlo=thext(thcount)
	      delta_th=thhi-thlo
	      do epcount=1,npmax
		 if( (ep.gt.epext(epcount)) .and. (ep.le.epext(epcount+1)) )then
		    ephi=epext(epcount+1)
		    eplo=epext(epcount)
		    delta_ep=ephi-eplo

		    if(usflag>0) then
		       A1 = extratus(thcount,epcount)
		       A2 = extratus(thcount+1,epcount)
		       A3 = extratus(thcount,epcount+1)
		       A4 = extratus(thcount+1,epcount+1)
		    else
		       A1 = extratds(thcount,epcount)
		       A2 = extratds(thcount+1,epcount)
		       A3 = extratds(thcount,epcount+1)
		       A4 = extratds(thcount+1,epcount+1)
		    endif


		    if(A1.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A1'
		       if(A2.ne.10.0) then
			  A1 = A2
		       endif
		    endif

		    if(A2.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A1'
		       if(A1.ne.10.0) then
			  A2 = A1
		       endif
		    endif

		    if(A3.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A1'
		       if(A4.ne.10.0) then
			  A3 = A4
		       endif
		    endif

		    if(A4.eq.10.0) then
		       write(6,*) 'Warning: unphysical table value for A1'
		       if(A3.ne.10.0) then
			  A4 = A3
		       endif
		    endif

		    if((A1.ne.10.0).and.(A2.ne.10.0).and.(A3.ne.10.0).and.(A4.ne.10)) then
		       A12 = (A1*(thhi-thdeg) + A2*(thdeg-thlo))/delta_th
		       A34 = (A3*(thhi-thdeg) + A4*(thdeg-thlo))/delta_th
		       A1234 = (A12*(ephi-ep) + A34*(ep-eplo))/delta_ep
		    elseif ( (A1.eq.10.0).or.(A2.eq.10.0)) then
		       A34 = (A3*(thhi-thdeg) + A4*(thdeg-thlo))/delta_th
		       A1234 = A34
		       write(6,*) 'Radiative corrections not defined for A12 theta, using A34'

		    elseif ( (A3.eq.10.0).or.(A4.eq.10.0)) then
		       A12 = (A1*(thhi-thdeg) + A2*(thdeg-thlo))/delta_th
		       A1234 = A12
		       write(6,*) 'Radiative corrections not defined for A34 theta, using A12'
		    endif

		 endif   !ep check
	      enddo      !loop over ep
	   endif    !theta check
	enddo       !loop over theta

	ext_func=A1234

	return
	end


*************************************************************************************
	real*8 function hmscereff_xem2(delta)

	implicit none

	real*8 pcent,delta
	real*8 a1,b1,b2,c1
	real*8 eff


	parameter(a1=0.9994)

        parameter(b1= 0.99727)
	parameter(b2= -0.00213)

	parameter(c1=0.9962)

* -------low delta parametrisation (delta < 1%)-------------------------------------

	if(delta.le.-1.0) then
	   eff=a1
* -------middle delta parametrisation (1% to 0.5%)-------------------------------------
	elseif(delta.gt.-1.0 .and. delta.lt.0.5) then
	   eff=b1+b2*delta
* -----hi delta  (>0.5%)---------------------------------------------------
	else
	   eff=c1
	endif

	hmscereff_xem2=eff

	return
	end



	
*************************************************************************************
	real*8 function shmscereff(delta,A)

	implicit none

	real*8 delta,A
	real*8 P1,P2,P3,P4,P5,P6,P7
	real*8 d1,d2

	if(A.gt.3.0) then
	   P1=0.99253
	   P2=0.99857
	   P3=-0.25180E-03
	   P4 =-0.33961E-04
	   P5=1.0108   
	   P6=-0.25214E-01 
	   P7=0.40108E-02
	   d1=5.5
	   d2=0.4
	else
	   P1=0.99233      
	   P2=0.99922      
	   P3=0.61793E-04  
	   P4=-0.10144E-04 
	   P5=0.99617      
	   P6=-0.67952E-02 
	   P7=0.10429E-02  
	   d1=6.0
	   d2=-0.5
	endif




	if(delta.gt.d1) then
	   shmscereff=p1
	elseif(delta.lt.d2) then
	   shmscereff=p2+p3*delta+p4*delta**2
	else
         shmscereff=p5+p6*delta+p7*delta**2
	endif

	return
	end

	real*8 function abishekcereff(delta,A)
	real*8 delta,A,x,eff

	p0=2.512
	p1=2.396
	p2=1.058
	p3=-0.02447
	p4=0.9958
	p5=-0.000122
	p6=0.9995

	x=delta
	if(a.gt.9.5) then
	   p0=2.512
	   p1=2.396
	   p2=1.058
	   p3=-0.02447
	   p4=0.9958
	   p5=-0.000122
	   p6=0.9995
	elseif(a.lt.3) then
	   p0=1.679
	   p1=3.806
	   p2=0.8675
	   p3=-0.01079
	   p4=0.9973
	   p5=-0.0002058
	   p6=0.9994
	elseif(a.gt.8.5 .and. a.lt.9.5) then
	   p0=2.283
	   p1=2.305
	   p2=0.9276
	   p3=-0.01993
	   p4=0.9941
	   p5=-4.483E-05
	   p6=0.9991
	endif

	if(x.le.p0) then
	   eff = p6+p3*exp(-((x-p0)/p2)**2)
	else
	   eff = p3*exp( -((x-p0)/p1)**2) 
     >         + (p4+x*p5)
	endif

	abishekcereff = eff
c      write(6,*) eff

	return
	end

	logical function shms_dipole_exit(xfp,xpfp,yfp,ypfp)

	real*8 xfp,xpfp,yfp,ypfp
	real*8 xdip,ydip
	real*8 crad,hwid,voffset
	real*8 zdip
	logical inside

	parameter(zdip=-307.0)
	parameter(crad=23.81)
	hwid=11.549/2.0

	xdip=xfp+zdip*xpfp
	ydip=yfp+zdip*ypfp

	inside=.true.
	if(abs(ydip).lt.hwid) then
	   if(abs(xdip).gt.(crad+voffset)) then
	      inside=.false.
	   endif
	else
	   if(ydip.ge.hwid) then
	      if( (xdip-voffset)**2+(ydip-hwid)**2 .gt. crad**2 ) then
		 inside=.false.
	      endif
	   endif

	   if(ydip.le.-hwid) then
	       if( (xdip-voffset)**2+(ydip+hwid)**2 .gt. crad**2 ) then
		 inside=.false.
	      endif
	   endif
	endif

	shms_dipole_exit=inside

	return
	end



PRO scaling,check=check,quick=quick,load=load
	
	if keyword_set(load)then goto,skipall
	
	interelm= 1
	PsepAbel= [ 4.31     , 4.05     , 5.80     , 5.77     ,$
	            8.99     , 14.21    , 5.91     , 7.67     ,$
		    5.67     , 6.97     ]
	shot_db = [ 34108    , 34358	, 30776    , 30306    ,$
	            35158    , 30506    , 34973    , 35356    ,$
		    30298    , 30307    ] 
	twindow = [ [3.6,4.6], [4.0,6.0], [3.7,3.9], [3.8,4.1],$
	            [2.8,3.1], [2.8,3.1], [3.4,3.9], [2.5,3.5],$
		    [3.8,4.2], [3.8,4.1]]
	los     = [ 'ROV-14' , 'ROV-13' , 'ROV014' , 'ROV014' ,$
	            'ROV-14' , 'ROV014' , 'ROV-13' , 'ROV-14' ,$
		    'ROV014' , 'ROV014'] 
	tran    = [ 1.0/0.69 , 1.0/0.67 , 1.0	   , 1.0      ,$
	            1.0      , 1.0	, 1.0/0.67 , 1.00     ,$
		    1.0      , 1.0     ]
	chan    = [ 0	     , 19	, 0	   , 0        ,$
	            22       , 0	, 19	   , 22       ,$
		    0        , 0       ]
	diag    = [ 'FVL'    , 'EVL'	, 'EVL'    , 'EVL'    ,$
	            'FVL'    , 'EVL'    , 'EVL'    , 'FVL'    ,$
		    'EVL'    ,'EVL'   ] 
	useevl  = [ 0	     , 1	, 0	   , 0        ,$
	            1	     , 0	, 1	   , 1        ,$
		    0        , 0       ]
	tmax    = [4.0       , 4.0      , 4.0      , 4.0      ,$
	           4.0       , 4.0      , 4.0      , 4.0      ,$
		   4.0       , 4.0    ]
	tmin    = [0.0       , 0.0      , 0.0      , 0.0      ,$
	           0.0       , 0.0      , 0.0      , 0.0      ,$
		   0.0       , 0.0    ]
	temin   = [3.3       , 3.3      , 3.2      , 3.4      ,$
	           3.5       , 3.1      , 3.1      , 3.3      ,$
		   3.3       , 3.4     ]
	temax   = [3.8       , 3.8      , 3.5      , 3.7      ,$
	           3.6       , 3.3      , 3.6      , 3.8      ,$
		   3.5       , 3.6     ]



	a = fltarr(n_elements(shot_db))
	k = fltarr(n_elements(shot_db))
	ip = fltarr(n_elements(shot_db))
	psep = fltarr(n_elements(shot_db))
	nesep = fltarr(n_elements(shot_db))
	cn_mean = fltarr(n_elements(shot_db))
	cn_err = fltarr(n_elements(shot_db))
	for i=0,n_elements(shot_db)-1 do begin
		if keyword_set(check)then begin
			debug=1
			if shot_db[i] ne check then goto,skip
		endif

		t0 = twindow[0,i]
		t1 = twindow[1,i]
		; print kappa and amin
		read_signal_mrm,0L,shot_db[i],'TOT','kappa',time1,kappa,1
		id = where(time1 ge t0 and time1 le t1)
		k[i] = mean(kappa[id])

		read_signal_mrm,0L,shot_db[i],'TOT','ahor',time1,amin,1
		id = where(time1 ge t0 and time1 le t1)
		a[i] = mean(amin[id])

		read_signal_mrm,0L,shot_db[i],'TOT','P_TOT',time1,Ptot,1
		id  = where(time1 ge t0 and time1 le t1)
		pin = mean(ptot[id])

		if max(time1) lt t1 then begin
			read_signal_mrm,0L,shot_db[i],'NIS','PNI',time1,PNBI,1
			id  = where(time1 ge t0 and time1 le t1)
			pni = mean(pnbi[id])
		
			read_signal_mrm,0L,shot_db[i],'ICP','PICRN',time1,PICRN,1
			id  = where(time1 ge t0 and time1 le t1)
			pic = mean(picrn[id])

			read_signal_mrm,0L,shot_db[i],'ECS','PECRH',time1,PECRH,1
			id  = where(time1 ge t0 and time1 le t1)
			pec = mean(pecrh[id])

			pin = pni + pic + pec
		endif
		read_signal_mrm,0L,shot_db[i],'BPD','Prad',time1,Prad,1
		id  = where(time1 ge t0 and time1 le t1)
		pr  = mean(prad[id])-1.0
		
		psep[i] = pin - pr

		read_signal_mrm,0L,shot_db[i],'MAG','Ipa',time1,curr,1
		id = where(time1 ge t0 and time1 le t1)
		ip[i] = mean(curr[id])
		
		augped,shot_db[i],dens
		nesep[i]=dens
			
		
		if keyword_set(medpow)then begin
			if psep[i] ge 4.0e6 and psep[i] le 6.0e6 then begin
				print,'Calculating'
			endif else goto,skip
		endif
		if keyword_set(debug)then begin
			print,'Pinj: ',pin
			print,'Prad: ',pr
			print,'Psep: ',psep[i]
			print,'Ip  : ',ip[i]
			print,'nspx: ',nesep[i]
			bp = ((4*!pi*0.0000001) * ip[i])/(2.0*!pi * a[i] * sqrt(0.5*(1+k[i]^2)))
			ngw= (nesep[i]/1e20) / ((ip[i]/1e6) / !pi / a[i]^2)
			print,'1/(bp*(1+k^2)^1.5 * fgw^2): ',1/(bp * (1+k[i]^2)^(1.5) * ngw^2) 
			print,'Gold: ',4 * ((Psep[i]/1e6) / (bp * (1+k[i]^2)^(1.5) * ngw^2))/18.3
		endif
		if ~keyword_set(quick)then begin
		analysis,shot_db[i],time=time,cn_up=cn_up,cn_low=cn_low,sv=10,lowerte=temin[i],$
				    interelm=interelm,/no402,use_evl=useevl[i],los=los[i],trans=tran[i],$
				    channel=chan[i],diag=diag[i],tdiv=tdiv,xr=[t0,t1],upperte=temax[i]
		tdiv_range = [tmin[i],tmax[i]]
		id = where(tdiv ge tmin[i] and tdiv le tmax[i])
		stats = moment([cn_up[id]*100,cn_low[id]*100])
		cn_mean[i] = stats[0]
		cn_err[i]  = sqrt(stats[1])
		print,'Spec: ',cn_mean[i]
		print,'Err : ',cn_err[i]
		
		if keyword_set(debug)then begin
			setgraphics,xs=xs,ys=ys,ncol=1,nrow=1,psplot=psplot,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list	
			user_psym,4,/fill
			cn_range = [3,5]
			plot,tdiv,cn_up*100,psym=8,col=colors.black,back=colors.white,xr=[0,20],yr=[0,50],$
			xtitle='T!ldiv!n [eV]',ytitle='c!lN!n [%]'
			oplot,tdiv,cn_low*100,psym=8,col=colors.black
			oplot,[tdiv_range[0],tdiv_range[0]],[0,100],linest=5,col=colors.black
			oplot,[tdiv_range[1],tdiv_range[1]],[0,100],linest=5,col=colors.black
			oplot,[mean(tdiv_range),mean(tdiv_range)],[stats[0],1000],psym=8,col=colors.red
			err_plot,[mean(tdiv_range),mean(tdiv_range)],[stats[0],1000],[sqrt(stats[1]),0.1],col=colors.red
			cursor,x,y,/up
			
			if !MOUSE.button eq 4 then stop
		endif
		endif 
		if keyword_set(debug)then stop
		skip:
	endfor
	for i=0,n_elements(psep)-1 do begin
		print,'---------------------------'
		print,'Shot: ',shot_db[i]
		print,'Psep: ',psep[i]
		print,'Ip  : ',ip[i]
		print,'k   : ',k[i]
		print,'a   : ',a[i]
		print,'nspx: ',nesep[i]
		bp = ((4*!pi*0.0000001) * ip[i])/(2.0*!pi * a[i] * sqrt(0.5*(1+k[i]^2)))
		ngw= (nesep[i]/1e20) / ((ip[i]/1e6) / !pi / a[i]^2)
		print,'Bp   : ',bp
		print,'fnsep,GW: ',ngw
		print,'1/(bp*(1+k^2)^1.5 * fgw^2): ',1/(bp * (1+k[i]^2)^(1.5) * ngw^2) 
		print,'Gold: ',4 * ((Psep[i]/1e6) / (bp * (1+k[i]^2)^(1.5) * ngw^2))/18.3
		print,'cN spec: ',cn_mean[i],' +/- ',cn_err[i]
	endfor
	if keyword_set(quick)then stop


	; Medium power scaling	
	psplot=0
	rep:
	
	bp = ((4*!pi*0.0000001) * ip)/(2.0*!pi * a * sqrt(0.5*(1+k^2)))
	ngw = (nesep/1e20) / ((ip/1e6) / !pi / a^2)
	if keyword_set(psplot)then begin
		xs = 6 & ys =8
	endif else begin
		xs = 600 & ys =800
	end
		
	skipall:
	if keyword_set(load)then begin
		restore,'save/cn_database.sav'	
		psep = database.psep
		bp=database.bp
		cn_mean=database.cn_mean
		cn_err = database.cn_err
		ngw = database.ngw
		k = database.k
		a = database.a
		
	endif
	if keyword_set(medpow)then id = where(Psep ge 4.0e6 and Psep le 6.0e6) else id = where(Psep ge 0)
	xr = [0,50]
	xpow = 1.0
	xprc = 4.0
	gd_scal = (10.7^(xpow))/(0.34 * (1+1.63^2)^1.5 * 0.496^2)
	gd_perc = xprc
	xaxis = (Psep[id]/1e6)^(xpow)/ (bp[id] * (1+k[id]^2)^(1.5) * ngw[id]^2) * (gd_perc/gd_scal)
	yaxis = cn_mean[id]
	err   = cn_err[id]
	fit   = linfit(xaxis,yaxis,measure_errors=err)
	xarr  = findgen(1000)/999.0*(max(xaxis)-min(xaxis))+min(xaxis)
	
	setgraphics,xs=xs,ys=ys,ncol=2,nrow=1,psplot=psplot,filename='figures/scaling.ps',colors=colors,colpick=colpick,collabel=collabel,full_list=full_list	
	user_psym,1,/fill
	plot,xaxis,yaxis,col=colors.black,back=colors.white,xr=xr,yr=xr,psym=8,xtitle='Y * P!lsep!uX!n/(B!lp!n(1+k!u2!n)!u1.5!nf!u2!lGW!n)',ytitle='c!lN!n [%]'
	err_plot,xaxis,yaxis,err,col=colors.black
	oplot,[0,100],[0,100],linest=5,col=colors.black
	xv = (max(xr)-min(xr))*0.1 + min(xr)
	yv1 = (max(xr)-min(xr))*0.8 + min(xr)
	yv2 = (max(xr)-min(xr))*0.7 + min(xr)
	xyouts,[xv,xv],[yv1,yv2],['X= 1.0',string(xprc,format='("Y= ",D3.1," %")')],col=[colors.black,colors.black]	
	xpow = 0.5
	xprc = 1.2
	gd_scal = (10.7^(xpow))/(0.34 * (1+1.63^2)^1.5 * 0.496^2)
	gd_perc = xprc
	xaxis = (Psep[id]/1e6)^(xpow) / (bp[id] * (1+k[id]^2)^(1.5) * ngw[id]^2) * (gd_perc/gd_scal)
	yaxis = cn_mean[id]
	err   = cn_err[id]
	fit   = linfit(xaxis,yaxis,measure_errors=err)
	xarr  = findgen(1000)/999.0*(max(xaxis)-min(xaxis))+min(xaxis)
	xr = [0,20]
	plot,xaxis,yaxis,xr=xr,yr=xr,col=colors.black,back=colors.white,psym=8,xtitle='Y * P!lsep!uX!n /(B!lp!n(1+k!u2!n)!u1.5!nf!u2!lGW!n)',ytitle='c!lN!n [%]'
	err_plot,xaxis,yaxis,err,col=colors.black
	xv = (max(xr)-min(xr))*0.1 + min(xr)
	yv1 = (max(xr)-min(xr))*0.8 + min(xr)
	yv2 = (max(xr)-min(xr))*0.7 + min(xr)
	xyouts,[xv,xv],[yv1,yv2],[string(xpow,format='("X= ",D3.1)'),string(xprc,format='("Y= ",D3.1," %")')],col=[colors.black,colors.black]
	oplot,[0,100],[0,100],linest=5,col=colors.black
	
	if ~keyword_set(psplot)then begin
		ans=''
		read,ans,prompt='Make postscript (y/n): '
		if ans eq 'y' or ans eq 'Y' then begin
			psplot=1
			goto,rep
		endif
	endif else setgraphics,psplot=psplot,/close
	Stop
END  	

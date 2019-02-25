PRO pscan,debug=debug,check=check,quick=quick
	
	interelm= 1
	
	shot_db = [ 34108    , 34358	, 30776    , 30306    , 33032	 , 32950    , 32244    , 32932    , 30296    ,$
	            35158    , 30505	, 30506    , 32952    , 34971	 , 34973    , 34275    , 34981    , 35356    ,$
		    30298    , 30308    , 30309    , 33254    , 30295 ] 
	twindow = [ [3.6,4.6], [4.0,6.0], [3.6,4.0], [3.0,4.2], [3.0,7.0], [4.0,6.0], [2.2,3.8], [4.5,6.0], [4.0,5.5],$
	            [1.7,2.5], [3.5,5.0], [2.4,3.0], [4.5,6.0], [3.0,4.0], [3.0,4.0], [3.5,4.0], [5.6,6.6], [2.5,4.0],$
		    [3.8,4.2], [3.3,4.2], [3.6,4.0], [2.5,3.3], [3.0,4.5]]
	los     = [ 'ROV-14' , 'ROV-13' , 'ROV014' , 'ROV014' , 'ROV012' , 'ROV014' , 'ROV014' , 'ROV012' , 'ROV014' ,$
	            'ROV-14' , 'ROV014' , 'ROV014' , 'ROV012' , 'ROV-13' , 'ROV-13' , 'ROV-13' , 'ROV-13' , 'ROV-14' ,$
		    'ROV014' , 'ROV014' , 'ROV014' , 'ROV013' , 'ROV014' ] 
	tran    = [ 1.0/0.69 , 1.0/0.67 , 1.0	   , 1.0      , 1.0	 , 1.0      , 1.0      , 1.0      , 1.0      ,$
	            1.0      , 1.0	, 1.0	   , 1.0      , 1.0/0.67 , 1.0/0.67 , 1.0/0.67 , 1.0/0.67 , 1.00     ,$
		    1.0      , 1.0      , 1.0      , 1.0      , 1.0    ]
	chan    = [ 0	     , 19	, 0	   , 0        , 0	 , 0	    , 9        , 0        , 0        ,$
	            22       , 0	, 0	   , 0        , 19	 , 19	    , 0        , 19       , 0        ,$
		    0        , 0        , 0        , 0        , 0       ]
	diag    = [ 'FVL'    , 'EVL'	, 'EVL'    , 'EVL'    , 'EVL'	 , 'EVL'    , 'EVL'    , 'EVL'    , 'EVL'    ,$
	            'FVL'    , 'EVL'	, 'EVL'    , 'EVL'    , 'EVL'	 , 'EVL'    , 'EVL'    , 'EVL'    , 'EVL'    ,$
		    'EVL'    , 'EVL'    , 'EVL'    , 'EVL'    , 'EVL'   ] 
	useevl  = [ 0	     , 1	, 0	   , 0        , 0	 , 0	    , 1        , 0        , 0        ,$
	            1	     , 0	, 0	   , 0        , 1	 , 1	    , 0        , 1        , 0        ,$
		    0        , 0        , 0        , 0        , 0     ]
	tmax    = [4.0       , 4.0      , 4.0      , 4.0      , 4.0      , 2.0      , 4.0      , 3.0      , 3.0      ,$
	           4.0       , 4.0      , 4.0      , 6.0      , 4.0      , 4.0      , 6.0      , 4.0      , 4.0      ,$
		   4.0       , 4.0      , 4.0      , 4.0      , 4.0  ]
	tmin    = [0.0       , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      ,$
	           0.0       , 0.0      , 0.0      , 5.0      , 0.0      , 0.0      , 5.0      , 0.0      , 0.0      ,$
		   0.0       , 0.0      , 0.0      , 0.0      , 0.0  ]
	temin   = [3.3       , 3.3      , 3.3      , 3.3      , 3.3      , 3.3      , 3.3      , 3.3      , 3.3      ,$
	           3.5       , 3.3      , 3.2      , 3.8      , 3.3      , 3.3      , 3.3      , 3.3      , 3.3      ,$
		   3.3       , 3.3      , 3.3      , 3.3      , 3.3 ]

	setgraphics,xs=xs,ys=ys,ncol=1,nrow=1,psplot=psplot,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list	


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
			
		
		
		if ~keyword_set(quick)then begin
			analysis,shot_db[i],time=time,cn_up=cn_up,cn_low=cn_low,sv=10,lowerte=temin[i],$
				    interelm=interelm,/no402,use_evl=useevl[i],los=los[i],trans=tran[i],$
				    channel=chan[i],diag=diag[i],tdiv=tdiv,xr=[t0,t1]
			tdiv_range = [tmin[i],tmax[i]]
			id = where(tdiv ge tmin[i] and tdiv le tmax[i])
			stats = moment([cn_up[id]*100,cn_low[id]*100])
			cn_mean[i] = stats[0]
			cn_err[i]  = sqrt(stats[1])
		
			if keyword_set(debug)then begin
				user_psym,4,/fill
				cn_range = [3,5]
				plot,tdiv,cn_up*100,psym=8,col=colors.black,back=colors.white,xr=[0,20],yr=[0,50],$
				xtitle='T!ldiv!n [eV]',ytitle='c!lN!n [%]'
				oplot,tdiv,cn_low*100,psym=8,col=colors.black
				oplot,[tdiv_range[0],tdiv_range[0]],[0,100],linest=5,col=colors.black
				oplot,[tdiv_range[1],tdiv_range[1]],[0,100],linest=5,col=colors.black
				oplot,[mean(tdiv_range),mean(tdiv_range)],[stats[0],1000],psym=8,col=colors.red
				err_plot,[mean(tdiv_range),mean(tdiv_range)],[stats[0],1000],[sqrt(stats[1]),0.1],col=colors.red
			endif
		endif	
		print,'Pinj: ',pin
		print,'Prad: ',pr
		print,'Psep: ',psep[i]
		print,'Ip  : ',ip[i]
		print,'nspx: ',nesep[i]
		bp = ((4*!pi*0.0000001) * ip[i])/(2.0*!pi * a[i] * sqrt(0.5*(1+k[i]^2)))
		ngw= (nesep[i]/1e20) / ((ip[i]/1e6) / !pi / a[i]^2)
		print,'1/fsp/bp: ',1/ngw^2/bp
		print,'Gold: ',4 * ((Psep[i]/1e6) / (bp * (1+k[i]^2)^(1.5) * ngw^2))/18.3
		if ~keyword_set(quick)then begin
			print,'Spec: ',cn_mean[i]
			print,'Err : ',cn_err[i]
			if keyword_set(debug)then begin
				cursor,x,y,/up	
				if !MOUSE.button eq 4 then stop				
			endif
		endif
		skip:
	endfor
	
	bp = ((4*!pi*0.0000001) * ip)/(2.0*!pi * a * sqrt(0.5*(1+k^2)))
	ngw = (nesep/1e20) / ((ip/1e6) / !pi / a^2)
	scaler = (Psep/1e6) / (bp * (1+k^2)^(1.5) * ngw^2)
	xarr = findgen(1000)/999.0*(max(scaler)-min(scaler))+min(scaler)
	scaler1 = (Psep/1e6) / ngw^2
	xarr1 = findgen(1000)/999.0*(max(scaler1)-min(scaler1))+min(scaler1)
	scaler2 = (Psep/1e6) / bp / ngw^2
	xarr2 = findgen(1000)/999.0*(max(scaler2)-min(scaler2))+min(scaler2)
	scaler3 = 4.0 * scaler / 18.3
	xarr3 = findgen(1000)/999.0*(max(scaler3)-min(scaler3))+min(scaler3)

	result = linfit(scaler,cn_mean,measure_errors=cn_err)
	result1 = linfit(scaler1,cn_mean,measure_errors=cn_err)
	result2 = linfit(scaler2,cn_mean,measure_errors=cn_err)
	result3 = linfit(scaler3,cn_mean,measure_errors=cn_err)
	
	fit =  result[1] * xarr + result[0]
	fit1 =  result1[1] * xarr1 + result1[0]
	fit2 =  result2[1] * xarr2 + result2[0]
	fit3 =  result3[1] * xarr3 + result3[0]
	user_psym,3,/fill
	
	setgraphics,xs=1000,ys=800,ncol=2,nrow=2,psplot=psplot,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list	

	plot,scaler,cn_mean,col=colors.black,back=colors.white,psym=8,yr=[0,40],xtitle='P!lsep!n/(B!lp!n(1+k!u2!n)!u1.5!nf!uGW!lne,sep!u2!n)',ytitle='c!lN!n [%]'
	err_plot,scaler,cn_mean,cn_err,col=colors.black
	oplot,xarr,fit,linest=5,col=colors.black

	plot,scaler1,cn_mean,col=colors.black,back=colors.white,psym=8,yr=[0,40],xtitle='P!lsep!n/(f!uGW!le,sep!u2!n)',ytitle='c!lN!n [%]'
	err_plot,scaler1,cn_mean,cn_err,col=colors.black
	oplot,xarr1,fit1,linest=5,col=colors.black

	plot,scaler2,cn_mean,col=colors.black,back=colors.white,psym=8,yr=[0,40],xtitle='P!lsep!n/(B!lp!nf!uGW!lne,sep!u2!n)',ytitle='c!lN!n [%]'
	err_plot,scaler2,cn_mean,cn_err,col=colors.black
	oplot,xarr2,fit2,linest=5,col=colors.black

	plot,scaler3,cn_mean,col=colors.black,back=colors.white,psym=8,xr=[0,100],yr=[0,100],xtitle='4%/18.3 * P!lsep!n/(B!lp!n(1+k!u2!n)!u1.5!nf!uGW!lne,sep!u2!n)',ytitle='c!lN!n [%]'
	err_plot,scaler3,cn_mean,cn_err,col=colors.black
	oplot,xarr3,fit3,linest=5,col=colors.black
	oplot,[0,100],[0,100],col=colors.red,linest=5
	for i=0,n_elements(shot_db)-1 do print,string(shot_db[i],cn_mean[i],k[i],a[i],psep[i]/1e6,ngw[i],scaler[i],format='("Shot: ",i5,"; cN= ",d7.2,"; kappa= ",d7.2,"; a= ",d7.2,"; Psep= ",d7.2,"; nesep,GW =",d7.2,"; scaler=",d7.2)')
	
	STOP
END  	

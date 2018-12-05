FUNCTION fetch_data,shot, sig,tr=tr,$
                          backc=backc,$
			  fitback=fitback,$
			  mdlfile=mdlfile,$
			  psplot=psplot,$
			  calwave=calwave,$
			  debug=debug,$
			  wshift=wshift,$
			  diag=diag,$
			  resolution=resolution,$
			  nomodel=nomodel,$
			  use_tau=use_tau, $
                          keep_file = keep_file,$
			  overview=overview,$
			  interelm=interelm,$
			  elmcond=elmcond,$
			  machine=machine,$
			  rrange=rrange,$
			  jet_coord=jet_coord,$
			  fitall=fitall,$
			  kt3a=kt3a,neon=neon,$
			  spec=spec,use_bart=use_bart
			  
			  
    	!QUIET=1
	!EXCEPT=0
    	
    	IF ~KEYWORD_SET(machine)THEN machine = 'JET'
	instr_lorz  = 0.1
	
	IF KEYWORD_SET(calwave) THEN BEGIN
	    nomodel=1
	    calwave=STRING(shot,FORMAT='(I5)')
	ENDIF    

;	**********************************
;	**** Read data ****
;	**********************************
	data = get_data(shot,machine,diag=diag)
	
	
	
;	**********************************
;	**** Setup spectral line fitting *
;	**********************************

    	gauss = gauss_lines(species=['He','C','O','N','W'],instr_func=data.instr_func,diag=diag)
	voigt = voigt_lines(species=['D'],instr_func=data.instr_func,diag=diag)
	
;	**********************************
;	**** Analyse time resolution  ****
;	**********************************

	FOR i=0,N_ELEMENTS(red_los)-1 DO PRINT,STRING(i,FORMAT='("LOS(",i2,"):")')+red_los(i)
	PRINT,'Time resolution: ',data.time(1)-data.time(0)

;	**********************************
;	**** Pick and sort the sightlines*
;	**********************************

	id_sig   = WHERE(STRPOS(red_los,sig) NE -1)	
	IF keyword_set(rrange) and machine eq 'JET' then id_sig=where(jet_coord ge rrange[0] and jet_coord le rrange[1])
    	jet_coord = jet_coord[id_sig]
	nsig      = N_ELEMENTS(id_sig)
	arr       = 0
	IF machine EQ 'AUG' THEN BEGIN
	    FOR i=0,nsig-1 DO arr=[arr,LONG(STRMID(red_los(id_sig[i]),4,2))] & arr = arr[1:*]
	ENDIF ELSE BEGIN
	    FOR i=0,nsig-1 DO arr=[arr,LONG(STRMID(red_los(id_sig[i]),2,2))] & arr = arr[1:*]
	ENDELSE
	PRINT,'Using sightlines: '
	FOR i=0,nsig-1 DO print,red_los(id_sig[i])
	srt    = SORT(arr) & id_sig = id_sig[srt] & arr = arr[srt]
	!P.CHARSIZE=2.0
	IF id_sig[0] NE -1 THEN BEGIN 
	    	IF ~KEYWORD_SET(tr)THEN tr=[MIN(data.time),MAX(data.time)]
		idt = WHERE(data.time GE tr[0] AND data.time LE tr[1])
		IF idt[0] EQ -1 THEN BEGIN
		    PRINT,'Warning: Specified time range outside experiment range'
		    PRINT,'Defaulting to full diagnostic range.'
		    idt = INDGEN(N_ELEMENTS(data.time))
		ENDIF    
		IF ~KEYWORD_SET(resolution) THEN resolution   = 0.0001
		act_resol    = data.time(1)-data.time(0)
		nires        = resolution/act_resol
		nires        = FLOOR(nires)>1
		PRINT,'Binning time slices every: ',nires
		icount_total = nsig * N_ELEMENTS(idt)/nires &	icount = 0.0
;	*************************************************************
;	**** Array details                                       ****
;	**** H I(0) :  n=7-2                                     ****
;	**** H I(1) :  n=6-2                                     ****
;	**** N II(0):  399.55 nm 2s2 2p 3p 1D --> 2s2 2p 3s 1P   ****
;	**** N II(1):  404.30 nm 2s2 2p 4f 3G --> 2s2 2p 3d 3F   **** 
;	**** N II(2):  402.65 nm 2s2 2p 4f 1G --> 2s2 2p 3d 3F   ****
;	**** N II(3):  407.40 nm 2s2 2p 4f 3F --> 2s2 2p 3d 3F   ****
;	**** N II(4):  408.85 nm 2s2 2p 4f 1F --> 2s2 2p 3d 3F   ****
;	**** N III(0): 409.70 nm 2s2 3p 2P    --> 2s2 3s 2S      ****
;	**** N III(1): 400.36 nm 2s2 5f 2F    --> 2s2 4d 2D      ****
;	**** N IV(0):  405.70 nm 1s2 2s 3d 1D --> 1s2 2s 3p 1P   **** 
;   	*************************************************************

		ntime          = N_ELEMENTS(data.time)
		wi_time        = FLTARR(ntime,nsig,3)
		wi_time_err    = FLTARR(ntime,nsig,3)
		nii_time       = FLTARR(ntime,nsig,5)
		nii_time_err   = FLTARR(ntime,nsig,5)
		niii_time      = FLTARR(ntime,nsig,2)
		niii_time_err  = FLTARR(ntime,nsig,2)
		niv_time       = FLTARR(ntime,nsig,1)
		niv_time_err   = FLTARR(ntime,nsig,1)
		nv_time        = FLTARR(ntime,nsig,1)
		nv_time_err    = FLTARR(ntime,nsig,1)
		hi_time        = FLTARR(ntime,nsig,2)
		hi_time_err    = FLTARR(ntime,nsig,2)
		te_time        = FLTARR(ntime,nsig)
		ne_time        = FLTARR(ntime,nsig)
		te_err_time    = FLTARR(ntime,nsig)
		ne_err_time    = FLTARR(ntime,nsig)
		ne_balmer      = FLTARR(ntime,nsig)
		ne_balmer_err  = FLTARR(ntime,nsig)
		n1_time        = FLTARR(ntime,nsig)
		n1_time_err_up = FLTARR(ntime,nsig)
		n1_time_err_lw = FLTARR(ntime,nsig)

;	**********************************
;	**** Calibrate wavelength     ****
;	**********************************
    	    	if spec eq 'kt3a' then wcal_file='Save/wcal_kt3a'+DIAG+STRING(shot,format='(I5)')+'.sav' else $
		                       wcal_file='Save/wcal'+DIAG+STRING(shot,format='(I5)')+'.sav'
	        
		IF FILE_TEST(wcal_file)THEN BEGIN
	            RESTORE,wcal_file,/verb 
	            PRINT,'Restored wavelength calibration file'
	        ENDIF ELSE wcal=0.0

		IF KEYWORD_SET(calwave)THEN wcal=0.0

		FOR i=0,nsig-1 DO BEGIN

;	**********************************
;	**** Retrieve data for sightline**
;	**********************************

		    IF (size(data.emiss))[0] GT 2 THEN BEGIN
    		    	wavelength = REFORM(data.wvlngth(*,id_sig[i]))
		   	emiss      = ABS(REFORM(data.emiss(*,id_sig[i],*)))
	    	   	yerr       = ABS(REFORM(data.error(*,id_sig[i],*)))
			time       = data.time
		    ENDIF ELSE BEGIN
    		    	wavelength = data.wvlngth
		   	emiss      = data.emiss
	    	   	yerr       = data.error
			time       = data.time
		    ENDELSE	
		    time       = data.time
		    IF KEYWORD_SET(overview)THEN BEGIN
			window,0,ys=900,xs=1200
			!p.multi=[0,1,2]
			tsel = 50
			tfind:
			id = where(abs(time-tsel) eq min(abs(time-tsel)))
			print,'Time: ',tsel
			plot,wavelength,emiss[*,id[0]],/ylog,yr=[1e16,1e20],xs=1
			CONTOUR,emiss>1e15<(mean(emiss)*4),wavelength,time,nle=50,/fill,xs=1
			oplot,[0,1000],[tsel,tsel],linest=5
			cursor,x,y,/up
			tsel = y
			goto,tfind
			STOP
		    ENDIF
		    IF ~KEYWORD_SET(calwave)THEN wavelength = wavelength+ wcal

;	**********************************
;	**** Store emiss traces       ****
;	**********************************
		    ires = 0
		    emissivity =0
		    yerror =0
		    FOR j=idt[0],MAX(idt) DO BEGIN
			    ires = ires + 1
			    emissivity = emissivity + emiss(*,j)
			    yerror     = yerror + yerr(*,j)
			    IF ires EQ nires THEN BEGIN
				emissivity = emissivity / nires
				yerror     = yerror / nires
				IF KEYWORD_SET(DEBUG)THEN PRINT,'Time / s: ',time(j)
				idnan = WHERE(FINITE(emissivity) EQ 0 OR emissivity EQ 0,complement=idreal)
				IF idnan[0] NE -1 THEN emissivity[idnan]=INTERPOL(emissivity[idreal],idreal,idnan)
;	*************************************************************
;	**** Assess background                                   ****
;	*************************************************************

	    	    	    	IF KEYWORD_SET(backc)THEN useback=backc

;	*************************************************************
;	**** Fit spectra over 398 < lambda < 409 nm              ****
;	*************************************************************

				
				if keyword_set(kt3a) then nii_range=[363,390] else nii_range=[399,409] 
				if keyword_set(calwave) or keyword_set(fitall) then nii_range=[0,1000] 
				id = WHERE(wavelength GT nii_range[0] and wavelength lt nii_range[1] )

				IF keyword_set(calwave)then test=cal_wav(wavelength[id],emissivity[id],shotstr,diag,spec)
				id_reduce_gauss=WHERE(gauss.pos LE MAX(wavelength[id])and gauss.pos GE MIN(wavelength[id]))
				id_reduce_voigt=WHERE(voigt.pos LE MAX(wavelength[id])and voigt.pos GE MIN(wavelength[id]))
				gauss_new = {pos:gauss.pos[id_reduce_gauss],$
				             ion:gauss.ion[id_reduce_gauss],$
					     couple:gauss.couple[id_reduce_gauss],$
					     trn:gauss.trn[id_reduce_gauss],$
					     fwhm:gauss.fwhm,$
					     erc:gauss.erc,diag:gauss.diag}
				gauss_use_1=gauss_new
				if ~keyword_set(fitback)then fixback=0 else fixback=abs(1-fitback)
				if id_reduce_voigt[0] ne -1 then begin
					voigt_new = {pos:voigt.pos[id_reduce_voigt],$
				             fwhml:voigt.fwhml[id_reduce_voigt],$
				             fwhmg:voigt.fwhmg,$
					     couple:voigt.couple[id_reduce_voigt],$
					     nbalmer:voigt.nbalmer[id_reduce_voigt]}
					voigt_use_1=voigt_new
					params_1= run_ffs_fit(wavelength[id],emissivity[id],yerr=yerror[id],$
		                                      mdlfile=mdlfile,background=useback,fixback=fixback,$
						      instr_func=data.instr_func,debug=debug,voigt=voigt_use_1,$
				                      gauss=gauss_use_1,psplot=psplot,/nomodel,use_tau=use_tau)
				endif else $
					params_1= run_ffs_fit(wavelength[id],emissivity[id],yerr=yerror[id],mdlfile=mdlfile,$
								 background=useback,fixback=fixback,instr_func=data.instr_func,debug=debug,$
				                                 gauss=gauss_use_1,psplot=psplot,/nomodel,use_tau=use_tau)

;	*************************************************************
;	**** Store N II 402 and 404 nm lines                     ****
;	*************************************************************
				
				nii_time(j,i,1)     = params_1.n404  
				nii_time_err(j,i,1) = params_1.n404_err  
				nii_time(j,i,2)     = params_1.n402  
				nii_time_err(j,i,2) = params_1.n402_err  
				wi_time(j,i,0)      = params_1.nwi
				niv_time(j,i,0)     = params_1.nvi
				wi_time_err(j,i,0)  = params_1.nwi_err
				niv_time_err(j,i,0) = params_1.nvi_err
				nv_time(j,i,0)      = params_1.nv
				nv_time_err(j,i,0)  = params_1.nv_err
				IF mean(wavelength) gt 400 and mean(wavelength) lt 410 then begin				
								      
;	*************************************************************
;	**** Fit spectra over 396 < lambda < 400.5 nm            ****
;	*************************************************************
				id = WHERE(wavelength gt 396.5 and wavelength lt 400.5 ) 
				id_reduce_gauss=WHERE(gauss.pos LE MAX(wavelength[id])and gauss.pos GE MIN(wavelength[id]))
				gauss_new = {pos:gauss.pos[id_reduce_gauss],$
				             ion:gauss.ion[id_reduce_gauss],$
					     couple:gauss.couple[id_reduce_gauss],$
					     fwhm:gauss.fwhm,$
					     trn:gauss.trn[id_reduce_gauss],$
					     erc:gauss.erc,diag:gauss.diag}
				id_reduce_voigt=WHERE(voigt.pos LE MAX(wavelength[id])and voigt.pos GE MIN(wavelength[id]))
				
				voigt_new = {pos:voigt.pos[id_reduce_voigt],$
				             fwhml:voigt.fwhml[id_reduce_voigt],$
				             fwhmg:voigt.fwhmg,$
					     couple:voigt.couple[id_reduce_voigt],$
					     nbalmer:voigt.nbalmer[id_reduce_voigt]}
				if id_reduce_gauss[0] ne -1 then gauss_use_2=gauss_new
				if id_reduce_voigt[0] ne -1 then voigt_use_2=voigt_new
				newnorm         = max(emissivity[id])
				params_2        = run_ffs_fit(wavelength[id],emissivity[id],yerr=yerror[id],$
				                              mdlfile=mdlfile,instr_func=data.instr_func,debug=debug,$
							      gauss=gauss_use_2,voigt=voigt_use_2,use_tau=use_tau,$
				                              background=(params_1.background[0]*params_1.norm)/newnorm,$
							      fixback=1,psplot=psplot,/nomodel)
;	*************************************************************
;	**** Store N II 399 and Stark density                    ****
;	*************************************************************
				
				nii_time(j,i,0)      = params_2.n399  
				nii_time_err(j,i,0)  = params_2.n399_err  
				hi_time(j,i,0)       = params_2.h72
				hi_time_err(j,i,0)   = params_2.h72_err
				niii_time(j,i,1)     = params_2.niii2  
				niii_time_err(j,i,1) = params_2.niii2_err  
				
;	*************************************************************
;	**** Fit spectra over 406 < lambda < 411 nm              ****
;	*************************************************************

				id = WHERE(wavelength gt 408 and wavelength lt 411 ) 
				id_reduce_gauss=WHERE(gauss.pos LE MAX(wavelength[id])and gauss.pos GE MIN(wavelength[id]))
				gauss_new = {pos:gauss.pos[id_reduce_gauss],$
				             ion:gauss.ion[id_reduce_gauss],$
					     trn:gauss.trn[id_reduce_gauss],$
					     couple:gauss.couple[id_reduce_gauss],$
					     fwhm:gauss.fwhm,$
					     erc:gauss.erc,diag:gauss.diag}
				id_reduce_voigt=WHERE(voigt.pos LE MAX(wavelength[id])and voigt.pos GE MIN(wavelength[id]))
				
				voigt_new = {pos:voigt.pos[id_reduce_voigt],$
				             fwhml:voigt.fwhml[id_reduce_voigt],$
				             fwhmg:voigt.fwhmg,$
					     couple:voigt.couple[id_reduce_voigt],$
					     nbalmer:voigt.nbalmer[id_reduce_voigt]}
				if id_reduce_gauss[0] ne -1 then gauss_use_3=gauss_new
				if id_reduce_voigt[0] ne -1 then voigt_use_3=voigt_new
				
				newnorm         = max(emissivity[id])
				params_3        = run_ffs_fit(wavelength[id],emissivity[id],yerr=yerror[id],$
				                              mdlfile=mdlfile,instr_func=data.instr_func,debug=debug,$
							      gauss=gauss_use_3,voigt=voigt_use_3,$
				                              background=(params_1.background[0]*params_1.norm)/newnorm,fixback=1,$
						              psplot=psplot,/nomodel,use_tau=use_tau)

;	*************************************************************
;	**** Store N II 408 and 409 lines                        ****
;	*************************************************************
				
				hi_time(j,i,1)       = params_3.h62
				hi_time_err(j,i,1)   = params_3.h62_err
				niii_time(j,i,0)     = params_3.niii  
				niii_time_err(j,i,0) = params_3.niii_err  
				nii_time(j,i,3)      = params_1.n408    
				nii_time_err(j,i,3)  = params_1.n408_err    
				nii_time(j,i,4)      = params_3.n409   
				nii_time_err(j,i,4)  = params_3.n409_err       
				ne_balmer(j,i)       = params_2.balmer_ne
	    	    	    	ne_balmer_err(j,i)   = params_2.balmer_ne_err

    	    	    	    	IF KEYWORD_SET(debug)THEN BEGIN
				    IF KEYWORD_SET(psplot)THEN makeps,file='overview_spectra.ps',xs=8,ys=5
				    PRINT,'Ratio 404/399:', nii_time(j,i,1)/nii_time(j,i,0)
				    PRINT,'Ratio 404/402:', nii_time(j,i,1)/nii_time(j,i,2)
				    PRINT,'Ratio 408/399:', nii_time(j,i,3)/nii_time(j,i,0)
				    PRINT,'Ratio 404/408:', nii_time(j,i,1)/nii_time(j,i,3)
				    PRINT,'Ratio 404/409:', nii_time(j,i,1)/nii_time(j,i,4)
				    PRINT,'Ratio 402/399:', nii_time(j,i,2)/nii_time(j,i,0)
				    PRINT,'Ratio 400/410:', niii_time(j,i,1)/niii_time(j,i,0)
				    adas_colors,colors=colors
    	    	    	    	    plot,wavelength,emissivity,/ylog,xs=1,BACK=255,COL=0 ,xr=[396,MAX(wavelength)], YR=[1e16,1e20]    
				    id = WHERE(params_1.wavelength ge MAX(params_2.wavelength) AND params_1.wavelength LE MIN(params_3.wavelength))
				    oplot,params_1.wavelength[id],params_1.gaussians[id]*params_1.norm,col=colors.red
    	    	    	    	    oplot,params_3.wavelength,params_3.gaussians*params_3.norm,col=colors.red
    	    	    	    	    oplot,params_2.wavelength,params_2.gaussians*params_2.norm,col=colors.red
    	    	    	    	    oplot,params_2.wavelength,params_2.balmerfit*params_2.norm,col=colors.blue
    	    	    	    	    oplot,params_3.wavelength,params_3.balmerfit*params_3.norm,col=colors.blue
    	    	    	    	    oplot,params_1.wavelength[id],params_1.fitted[id]*params_1.norm,col=colors.green
   	    	    	    	    oplot,params_2.wavelength,params_2.fitted*params_2.norm,col=colors.green
    	    	    	    	    oplot,params_3.wavelength,params_3.fitted*params_3.norm,col=colors.green
				    stop
				ENDIF

;	*************************************************************
;	**** Fit full N II spectra to determine Te and ne        ****
;	*************************************************************
				
			    	endif
				ires       = 0
			    	emissivity = 0.0
			    	IF icount EQ 0 THEN progress,0.0,/reset,label='Progress (%)',frequency=1.0
		    	    	icount     = icount + 1
		    	    	progress,100.0*icount/icount_total
			    ENDIF
		    ENDFOR
		ENDFOR
		progress,100.0,/last
	ENDIF ELSE STOP,'Error: No sightlines exist for this signal ...'

;	*************************************************************
;	**** Prepare data for output structure                   ****
;	*************************************************************
	id_use         = FINDGEN(N_ELEMENTS(idt)/nires)*nires+idt[0]+nires-1
	wi_time        = wi_time[id_use,*,*] 
	wi_time_err    = wi_time_err[id_use,*,*] 
	te_time        = te_time[id_use,*,*] 
	ne_time        = ne_time[id_use,*,*] 
	te_err_time    = te_err_time[id_use,*,*] 
	ne_err_time    = ne_err_time[id_use,*,*] 
	ne_balmer      = ne_balmer[id_use,*,*] 
	ne_balmer_err  = ne_balmer_err[id_use,*,*] 
	n1_time        = n1_time[id_use,*,*] 
	n1_time_err_up = n1_time_err_up[id_use,*,*] 
	n1_time_err_lw = n1_time_err_lw[id_use,*,*] 
	nii_time       = nii_time[id_use,*,*] 
	niii_time      = niii_time[id_use,*,*] 
	niv_time       = niv_time[id_use,*,*] 
	nii_time_err   = nii_time_err[id_use,*,*] 
	niii_time_err  = niii_time_err[id_use,*,*] 
	niv_time_err   = niv_time_err[id_use,*,*] 
	balmer         = hi_time[id_use,*,*] 
	balmer_err     = hi_time_err[id_use,*,*] 
	time           = time[id_use,*,*] 
	output = { te:te_time ,$
	           dens:ne_time,$
		   te_err:te_err_time,$
		   ne_err:ne_err_time,$
		   dens_balmer:ne_balmer,$
		   dens_balmer_err:ne_balmer_err,$
		   nconc:n1_time,$
		   nconc_err_up:n1_time_err_up,$
		   nconc_err_lw:n1_time_err_lw,$		   
		   wi:wi_time,$
		   wi_err:wi_time_err,$
		   nii:nii_time,$
		   nii_err:nii_time_err,$
		   niii:niii_time,$
		   niii_err:niii_time_err,$
		   niv:niv_time,$
		   niv_err:niv_time_err,$
		   nv:nv_time,$
		   nv_err:nv_time_err,$
		   balmer:balmer,$
		   balmer_err:balmer_err,$
		   time:time,$
		   wavelength:wavelength,$
		   los_names:red_los[id_sig],$
		   rvals:jet_coord}
	RETURN,output
END

























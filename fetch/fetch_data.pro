FUNCTION fetch_data,shot, sig,tr=tr,$
                          backc=backc,$
			  fitback=fitback,$
			  species=species,$
			  mdlfile=mdlfile,$
			  psplot=psplot,$
			  range=range,$
			  calwave=calwave,$
			  debug=debug,$
			  wshift=wshift,$
			  diag=diag,$
			  resolution=resolution,$
			  use_tau=use_tau, $
                          keep_file = keep_file,$
			  overview=overview,$
			  interelm=interelm,$
			  machine=machine,$
			  rrange=rrange,$
			  fitall=fitall,$
			  kt3a=kt3a,$
			  save=save,$
			  load=load,$
			  append=append,$
			  quick=quick,$
			  stark=stark,$
			  pixelview=pixelview
			  
			  
    	!QUIET=1
	!EXCEPT=0
    	
    	IF ~KEYWORD_SET(machine)THEN machine = 'JET'
	instr_lorz  = 0.1
	
	IF KEYWORD_SET(calwave) THEN BEGIN
	    nomodel=1
	    calwave=STRING(shot,FORMAT='(I5)')
	ENDIF    

	shotstr = string(shot,format='(I5)')

	if keyword_set(load)then begin
		if ~keyword_set(append)then append='data'
		restore,'save/'+shotstr+'/'+sig+'-'+append+'.idl',/verb
		return,output
	endif				  

	if ~keyword_set(mdlfile)then mdlfile = 'tmp/'+shotstr + sig + '.mdl'

;	**********************************
;	**** Read data ****
;	**********************************

	data = get_data(shot,machine,diag=diag,interelm=interelm,debug=debug,sig=sig,jet_coord=jet_coord)
		
;	**********************************
;	**** Setup spectral line fitting *
;	**********************************
	if ~keyword_set(species)then species=['C','O','N','W']
    	gauss = gauss_lines(species=species,instr_func=data.instr_func,diag=diag,/unknown)
	voigt = voigt_lines(species=['D'],instr_func=data.instr_func,diag=diag)
	
;	**********************************
;	**** Analyse time resolution  ****
;	**********************************

	FOR i=0,N_ELEMENTS(data.los)-1 DO PRINT,STRING(i,FORMAT='("LOS(",i2,"):")')+data.los(i)
	PRINT,'Time resolution: ',data.time(1)-data.time(0)

;	**********************************
;	**** Pick and sort the sightlines*
;	**********************************

	id_sig   = WHERE(STRPOS(data.los,sig) NE -1)	
	nsig      = N_ELEMENTS(id_sig)
	arr       = 0
	IF machine EQ 'AUG' THEN BEGIN
	    FOR i=0,nsig-1 DO arr=[arr,LONG(STRMID(data.los(id_sig[i]),4,2))] & arr = arr[1:*]
	ENDIF ELSE BEGIN
	    FOR i=0,nsig-1 DO arr=[arr,LONG(STRMID(data.los(id_sig[i]),2,2))] & arr = arr[1:*]
	ENDELSE
	PRINT,'Using sightlines: '
	FOR i=0,nsig-1 DO print,data.los(id_sig[i])
	srt    = SORT(arr) & id_sig = id_sig[srt] & arr = arr[srt]
	!P.CHARSIZE=2.0
	
	if keyword_set(pixelview)then pixelview,data,id_sig=id_sig,xr=tr
	
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
		ari_time       = FLTARR(ntime,nsig,2)
		ari_time_err   = FLTARR(ntime,nsig,2)
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
		wcal_file='tmp/wcal'+DIAG+shotstr+'.sav'
	        IF FILE_TEST(wcal_file)THEN BEGIN
	            RESTORE,wcal_file,/verb 
	            PRINT,'Restored wavelength calibration file'
	        ENDIF ELSE BEGIN
		    wcal=0.0
		    IF ~KEYWORD_SET(calwave)THEN BEGIN
		    	cont='N'
		    	read,cont,prompt='Wavelength calibration not available, continue (Y/N)?'
		    	if strlowcase(cont) ne 'y' then return,-1
		    ENDIF
		END
		
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
			window,0
			CONTOUR,emiss>1e15<(mean(emiss)*7),wavelength,time,nle=50,/fill,yr=tr,xtitle='Wavelength [nm]',ytitle='Time [s]'
			id_nii   = where(wavelength ge 399.3 and wavelength le 399.7)
			id_nii2  = where(wavelength ge 403.0 and wavelength le 405.0)
			id_niii  = where(wavelength ge 400.4 and wavelength le 400.5)
			id_niv   = where(wavelength ge 405.6 and wavelength le 405.7)
			id_di    = where(wavelength ge 396.5 and wavelength le 397.2)
			id_brem  = where(wavelength ge 401.0 and wavelength le 402.0)
			nii_avr  = time
			nii2_avr  = time
			niii_avr = time
			niv_avr  = time
			di_avr   = time
			brem_avr   = time
			window,2
			!p.multi=0
			spec_avr = wavelength
			for i=0,n_elements(wavelength)-1 do spec_avr[i]=int_tabulated(time,emiss[i,*])
			plot,wavelength,spec_avr,xs=1,xtitle='Wavelength [nm]',ytitle='Time integrated spectra'
;			for i=0,n_elements(time)-1 do nii_avr[i]  = int_tabulated(wavelength[id_nii],emiss[id_nii,i])
;			for i=0,n_elements(time)-1 do nii_avr[i]  = int_tabulated(wavelength[id_nii],emiss[id_nii,i])
;			for i=0,n_elements(time)-1 do nii2_avr[i]  = int_tabulated(wavelength[id_nii2],emiss[id_nii2,i])
;			for i=0,n_elements(time)-1 do niii_avr[i] = int_tabulated(wavelength[id_niii],emiss[id_niii,i])
;			for i=0,n_elements(time)-1 do niv_avr[i] = int_tabulated(wavelength[id_niv],emiss[id_niv,i])
;			for i=0,n_elements(time)-1 do di_avr[i]   = int_tabulated(wavelength[id_di],emiss[id_di,i])
;			for i=0,n_elements(time)-1 do brem_avr[i]   = int_tabulated(wavelength[id_brem],emiss[id_brem,i])
			window,1,xs=1000,ys=900
			!p.multi=[0,1,4]
			!p.charsize=2.0
			adas_colors,colors=colors
			plot,time,nii_avr,xtitle='Time [s]',ytitle='N II @ 399.5',xr=tr
			oplot,time,smooth(nii_avr,15),col=colors.red
			plot,time,nii2_avr/nii_avr,xtitle='Time [s]',ytitle='N II ratio',xr=tr,yr=[0,0.5]
			oplot,time,smooth(nii2_avr,15)/smooth(nii_avr,15),col=colors.red
			plot,time,niii_avr,xtitle='Time [s]',ytitle='N III @ 400.5',xr=tr
			oplot,time,smooth(niii_avr,15),col=colors.red
			plot,time,di_avr,xtitle='Time [s]',ytitle='D I @ 397',xr=tr
			oplot,time,smooth(di_avr,15),col=colors.red
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
				timeavr = (time[j] + time[j]+resolution)/2.0
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

				nii_range=[395,406.5] 
				if keyword_set(quick)then nii_range=[398,406] 
				if keyword_set(stark)then nii_range=[395,398] 
				if keyword_set(calwave)then nii_range=[0,1000] 
				id = WHERE(wavelength GT nii_range[0] and wavelength lt nii_range[1] )
				IF keyword_set(calwave)then begin
				    test=cal_wav(wavelength,emissivity,shotstr,diag)
				    return,-1
				endif
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
				
				hi_time(j,i,1)       = params_1.h62
				hi_time_err(j,i,1)   = params_1.h62_err
				ne_balmer(j,i)       = params_1.balmer_ne
	    	    	    	ne_balmer_err(j,i)   = params_1.balmer_ne_err
				niii_time(j,i,0)     = params_1.niii  
				niii_time_err(j,i,0) = params_1.niii_err  
				nii_time(j,i,3)      = params_1.n408    
				nii_time_err(j,i,3)  = params_1.n408_err    
				nii_time(j,i,4)      = params_1.n409   
				nii_time_err(j,i,4)  = params_1.n409_err 
				nii_time(j,i,0)      = params_1.n399  
				nii_time_err(j,i,0)  = params_1.n399_err  
				hi_time(j,i,0)       = params_1.h72
				hi_time_err(j,i,0)   = params_1.h72_err
				niii_time(j,i,1)     = params_1.niii2  
				niii_time_err(j,i,1) = params_1.niii2_err  
				nii_time(j,i,1)      = params_1.n404  
				nii_time_err(j,i,1)  = params_1.n404_err  
				nii_time(j,i,2)      = params_1.n402  
				nii_time_err(j,i,2)  = params_1.n402_err  
				nii_time(j,i,3)      = params_1.n395  
				nii_time_err(j,i,3)  = params_1.n395_err  
				wi_time(j,i,0)       = params_1.nwi
				niv_time(j,i,0)      = params_1.nvi
				wi_time_err(j,i,0)   = params_1.nwi_err
				niv_time_err(j,i,0)  = params_1.nvi_err
				nv_time(j,i,0)       = params_1.nv
				nv_time_err(j,i,0)   = params_1.nv_err
								      

				skipnext:
;	*************************************************************
;	**** Fit full N II spectra to determine Te and ne        ****
;	*************************************************************
				
				ires       = 0
			    	emissivity = 0.0
			    	yerror     = 0.0
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
	ari_time       = ari_time[id_use,*,*] 
	ari_time_err   = ari_time_err[id_use,*,*] 
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
		   ari:ari_time,$
		   ari_err:ari_time_err,$
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
		   rvals:data.rvals[id_sig],$
		   los_names:data.los[id_sig]}
	if keyword_set(save)then begin
		cmd = 'mkdir -p save/'+shotstr
		spawn,cmd
		if ~keyword_set(append)then append='data'
		file = 'save/'+shotstr+'/'+sig+'-'+append+'.idl'
		save,file=file,output
	endif
	RETURN,output
END

























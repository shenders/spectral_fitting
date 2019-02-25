PRO calc,shot,diag=diag,channel=channel,te_limits=te_limits,interelm=interelm,deltaL=deltaL,sv=sv,xr=xr,timearr_elm=timearr_elm,sig3995_elm=sig3995_elm,ratio_elm=ratio_elm,cn_up_elm=cn_up_elm,cn_low_elm=cn_low_elm,debug=debug

; Set defaults
	if ~keyword_set(diag)then diag='EVL'
	if ~keyword_set(elmdiag)then elmdiag=diag
	if ~keyword_set(channel)then channel=0
	if ~keyword_set(elmchannel)then elmchannel=channel
	if ~keyword_set(sv)then sv=4
	if keyword_set(debug)then begin	
		adas_colors,colors=colors
		window,0,xs=1400,ys=1000
		!p.multi=[0,2,2]
		!p.charsize=2.0
	endif
; Get 399.5 line
	read_signal_mrm,0L,shot,diag,'N_1_3995',timearr,sig3995,2		
	sig3995 = sig3995[*,channel]

; Get 404.1 line
	read_signal_mrm,0L,shot,diag,'N_1_4041',timearr,sig4041,2,text=text
	sig4041 = sig4041[*,channel]

; Get interELM signal
	if keyword_set(xr)then begin
		id = where(timearr ge xr[0] and timearr le xr[1])
		timearr = timearr[id]
		sig3995 = sig3995[id]
		sig4041 = sig4041[id]
	endif
	if keyword_set(interelm)then begin	
		telm       = find_elm(shot,timearr)
		if ~keyword_set(elmcond)then elmcond=4.5
		id     = where(telm ge elmcond)
		if id[0] ne -1 then begin
			timearr_elm = timearr[id]
			sig3995_elm = sig3995[id]
			sig4041_elm = sig4041[id]
		endif
	endif else begin
		timearr_elm = timearr
		sig3995_elm = sig3995
		sig4041_elm = sig4041
	end
	
	
	
	ratio     = sig4041/sig3995
	ratio_elm = smooth(sig4041_elm,sv)/smooth(sig3995_elm,sv)


; Get 4041/3995 intensity ratio

	if keyword_set(debug)then begin	
		plot,timearr,ratio,yr=[0,1.0],xr=xr,xs=1,ytitle='N II 4041/3995'
		oplot,timearr_elm,ratio_elm,col=colors.red,thick=2.0
	endif 
; Temperature upper and lower limits
	if keyword_set(te_limits)then begin
		te_lower = te_limits[0]
		te_upper = te_limits[1]
	endif else begin
		te_lower = 3.5
		te_upper = 4.0
	end
; Density function based on temperature limits
	basic_fit,te_val=te_lower,func=func_te_lower,dens=dens_te_lower
	basic_fit,te_val=te_upper,func=func_te_upper,dens=dens_te_upper

; Guess transmission
	transmission = 1.0/0.64

; Get nitrogen concentrations from spectroscopy

	densfit     = interpol(dens_te_lower,func_te_lower,smooth(ratio,sv))
	densfit_elm = interpol(dens_te_lower,func_te_lower,smooth(ratio_elm,sv))
	if keyword_set(debug)then begin	
		plot,timearr,densfit/1e14,xr=xr,xs=1,ytitle='n!le!n [10!u20!n m!u-3!n]',/nodata
		oplot,timearr_elm,densfit_elm/1e14,col=colors.red,thick=2.0
	endif
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_lower+fltarr(n_elements(timearr)),$
		dens=densfit,data=exc3995,block=15
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_lower+fltarr(n_elements(timearr)),$
		dens=densfit,data=rec3995,block=65
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_lower+fltarr(n_elements(timearr_elm)),$
		dens=densfit_elm,data=exc3995_elm,block=15
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_lower+fltarr(n_elements(timearr_elm)),$
		dens=densfit_elm,data=rec3995_elm,block=65
			      
	run_adas405,elem='n',year='96',uid='adas',te=te_lower+fltarr(n_elements(timearr)),$
		dens=densfit,frac=frac
	run_adas405,elem='n',year='96',uid='adas',te=te_lower+fltarr(n_elements(timearr_elm)),$
		dens=densfit_elm,frac=frac_elm

	pec3995 = frac.ion[*,1] * exc3995 + frac.ion[*,2] * rec3995

	pec3995_elm = frac_elm.ion[*,1] * exc3995_elm + frac_elm.ion[*,2] * rec3995_elm

	cn_low = 3.14 * 4.0 * smooth(sig3995,sv) * transmission / (densfit * pec3995) / DeltaL / (densfit * 1e6)

	cn_low_elm = 3.14 * 4.0 * smooth(sig3995_elm,sv) * transmission / (densfit_elm * pec3995_elm) / DeltaL / (densfit_elm * 1e6)

	densfit     = interpol(dens_te_upper,func_te_upper,smooth(ratio,sv))
	densfit_elm = interpol(dens_te_upper,func_te_upper,smooth(ratio_elm,sv))
	if keyword_set(debug)then oplot,timearr_elm,densfit_elm/1e14,col=colors.red,linest=5,thick=2.0

	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_upper+fltarr(n_elements(timearr)),$
		dens=densfit,data=exc3995,block=15
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_upper+fltarr(n_elements(timearr)),$
		dens=densfit,data=rec3995,block=65
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_upper+fltarr(n_elements(timearr_elm)),$
		dens=densfit_elm,data=exc3995_elm,block=15
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_upper+fltarr(n_elements(timearr_elm)),$
		dens=densfit_elm,data=rec3995_elm,block=65
			      
	run_adas405,elem='n',year='96',uid='adas',te=te_upper+fltarr(n_elements(timearr)),$
		dens=densfit,frac=frac
	run_adas405,elem='n',year='96',uid='adas',te=te_upper+fltarr(n_elements(timearr_elm)),$
		dens=densfit_elm,frac=frac_elm

	pec3995 = frac.ion[*,1] * exc3995 + frac.ion[*,2] * rec3995

	pec3995_elm = frac_elm.ion[*,1] * exc3995_elm + frac_elm.ion[*,2] * rec3995_elm

	cn_up = 3.14 * 4.0 * smooth(sig3995,sv) * transmission / (densfit * pec3995) / DeltaL / (densfit * 1e6)

	cn_up_elm = 3.14 * 4.0 * smooth(sig3995_elm,sv) * transmission / (densfit_elm * pec3995_elm) / DeltaL / (densfit_elm * 1e6)
	if keyword_set(debug)then begin
		plot,timearr,sig3995/1e19,xr=xr,xs=1,ytitle='N II 399.5 [10!u19!n ph/s/m!u2!n/sr]'
		oplot,timearr_elm,smooth(sig3995_elm/1e19,sv),col=colors.red,thick=2.0
	
		plot,timearr,cn_up*100,/nodata,xr=xr,yr=[0,40],xs=1,xtitle='Time [s]',ytitle='c!lN!n [%]'	
		oband,timearr_elm,cn_up_elm * 100,cn_low_elm * 100,/norm,col=colors.red
		
		stop
	endif
END
Pro cn,shot,diag=diag,channel=channel,interelm=interelm,sv=sv,xr=xr,debug=debug,te_limits=te_limits


; Set defaults
	if ~keyword_set(sv)then sv=4
	if shot ge 35156 and shot le 35167 then begin		
		diag    = ['EVL','EVL','EVL','EVL','FVL']
		channel = [21   , 20  , 22  , 23  , 22  ]
		deltaL  = [0.05  , 0.05 , 0.05 , 0.05 , 0.05 ]
	endif	
;	diag='FVL'
;	channel=22
;	debug=1
	istore = -1
	conc_lw= -1
	conc_up= -1
	ratio  = -1
	time   = -1
	intens = -1
	for i=0,n_elements(diag)-1 do begin
		calc,shot,diag=diag[i],channel=channel[i],deltaL=deltaL[i],te_limits=te_limits,interelm=interelm,sv=sv,xr=xr,timearr_elm=timearr_elm,sig3995_elm=sig3995_elm,ratio_elm=ratio_elm,cn_up_elm=cn_up_elm,cn_low_elm=cn_low_elm,debug=debug
		istore = [istore,i+fltarr(n_elements(cn_up_elm))]
		conc_up = [conc_up,cn_up_elm*100]
		conc_lw = [conc_lw,cn_low_elm*100]
		ratio   = [ratio,ratio_elm]
		intens  = [intens,sig3995_elm]
		time    = [time,timearr_elm]
	endfor	
	istore = istore[1:*] 
	conc_lw= conc_lw[1:*] 
	conc_up= conc_up[1:*] 
	ratio  = ratio[1:*]  
	time   = time[1:*]
	intens = intens[1:*]
	adas_colors,colors=colors
	window,/free,xs=600,ys=900
	!p.multi=[0,1,3]
	!p.charsize=2.0
	!p.thick=2.0
	colpick = [colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.orange,colors.red]
	for i=0,n_elements(diag)-1 do begin & $
		id = where(istore eq i) & $
		if i eq 0 then plot,time[id],ratio[id],/nodata,xr=xr,xs=1,yr=[0,max(ratio)],ys=1,back=colors.white,col=colors.black & $
		oplot,time[id],ratio[id],col=colpick[i]	& $	
	endfor
	for i=0,n_elements(diag)-1 do begin & $
		id = where(istore eq i) & $
		if i eq 0 then plot,time[id],intens[id]/1e19,/nodata,xr=xr,xs=1,yr=[0,max(intens)/1e19],ys=1,back=colors.white,col=colors.black & $
		oplot,time[id],intens[id]/1e19,col=colpick[i]	& $	
	endfor
	for i=0,n_elements(diag)-1 do begin & $
		id = where(istore eq i) & $
		if i eq 0 then plot,time[id],conc_up[id],/nodata,xr=xr,xs=1,yr=[0,40],col=colors.black,back=colors.white,ytitle='cN [%]',xtitle='Time [s]' & $
		oband,time[id],conc_up[id],conc_lw[id],/norm,col=colpick[i] & $	
	endfor
stop
END

Function get_kt3,pulse,debug=debug,spec=spec,rrange=rrange,psplot=psplot,sig=sig,nored=nored
    if ~keyword_set(spec)then spec='kt3b'
    data = agm_readspec(pulse,spec=strlowcase(spec))
    agm_process_data,data,/wavecal,/radcal,/calibrate,/cal_si
    wave = data.data.wave_cal
    wave = REVERSE(wave,1)
    emis = data.data.data
    emis = REVERSE(emis,1)
    time = data.data.time
    los  = 'SP'+STRCOMPRESS(STRING(data.data.track,FORMAT='(I2)'),/remove_all)
    rval = data.data.coord
    yerr = data.data.toterror
    yerr = REVERSE(yerr,1)
; Load in vessel details
    openr,unit,'/home/shenders/jams/jet_vessel_new.txt',/get_lun
    xwl = -1
    ywl = -1
    while ~eof(unit)do begin
    	readf,unit,xjt,yjt
	xwl = [xwl,xjt]
	ywl = [ywl,yjt]
    end
    xwl = xwl[1:*]
    ywl = ywl[1:*]
    close,unit
    free_lun,unit
    if keyword_set(rrange)then begin
	id = where(rval ge rrange[0] and rval le rrange[1])
    	if id[0] ne -1 then begin
	    emis = emis[*,id,*]
	    wave = wave[*,id]
	    los  = los[id]
	    yerr = yerr[*,id,*]
	    rval = rval[id]
	endif
    endif else begin
    	; skip first few channels
	nored=1
	if ~keyword_set(nored)then begin
    	    skip_chan = 3
	    emis = emis[*,skip_chan:*,*]
	    wave = wave[*,skip_chan:*]
	    los  = los[skip_chan:*]
	    yerr = yerr[*,skip_chan:*,*]
	    rval = rval[skip_chan:*]
    	endif
    end
    for i=0,n_elements(los)-1 do print,los(i),rval(i)
    if keyword_set(debug)then begin
    	if keyword_set(psplot)then makeps,file='kt3_geom.ps',xs=8,ys=5
	window,/free
	plot,xwl,ywl,xr=[2.26,3.06],xs=1,yr=[-1.75,-1.3],ys=1,/iso
	
	if keyword_set(rrange)then begin
		id = where(data.data.coord ge rrange[0] and data.data.coord le rrange[1])
		for i=id[0],max(id) do oplot,[data.data.geo.origin[0],data.data.geo.origin[0]+6*data.data.geo.vector[0,i]],$
	    	    	    [data.data.geo.origin[1],data.data.geo.origin[1]+6*data.data.geo.vector[1,i]],thick=5
    	endif else begin
	    	for i=0,21 do oplot,[data.data.geo.origin[0],data.data.geo.origin[0]+6*data.data.geo.vector[0,i]],$
	    	    	    [data.data.geo.origin[1],data.data.geo.origin[1]+6*data.data.geo.vector[1,i]]
	
	endelse
    endif
    strc = {phflx    : emis , $
    	    lamgrid  : wave , $
	    losnames : los  , $
	    time     : time , $
	    yerr     : yerr , $
	    rval     : rval   }
    
    RETURN,strc
END    	    
	    

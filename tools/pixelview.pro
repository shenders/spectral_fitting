Pro pixelview,data,xr=xr,psplot=psplot,id_sig=id_sig


    wave = data.wvlngth[*,0] ; function of wave, los
    emis = data.emiss[*,*,*] ; function of wave, los, time
    time = data.time         ; time
    rval = data.rvals
    
    if keyword_set(xr)then begin
    	id_time = where(time ge xr[0] and time le xr[1])
	time = time[id_time]
	emis = emis[*,*,id_time]
    endif
    miss_pix = 0
    emis = emis[*,miss_pix:*,*]
    los  = data.los[miss_pix:*]
    rval = rval[miss_pix:*]
    id_sig=id_sig[miss_pix:*]
    if keyword_set(id_sig)then begin
       emis = emis[*,id_sig,*]
       los  = los[id_sig]
       rval = rval[id_sig]
    endif    
    tint = total(emis,3)     ; sum dimension in time
    lint = total(tint,2)     ; sum dimension in los
    
    rep:
    setgraphics,xs=800,ys=600,colors=colors
    plot,wave,lint,/ylog,back=colors.white,col=colors.black
    cursor,x,y,/up
    
    id_pix = where(abs(wave-x) eq min(abs(wave-x)))
    intens = emis[id_pix[0],*,*]/1e18
    
    radial = fltarr(n_elements(emis[0,*,0]))
    for i=0,n_elements(radial)-1 do radial[i]=float(strmid(los[i],2,2))
    
    setgraphics,xs=800,ys=600,title='Left click for individual spectra; middle click to quit; right click to repeat',psplot=psplot,file='figures/pixelview.ps'
    shimage,reform(intens),radia,time,xtitle='R [m]',ytitle='Time [s]',title='N II [10!u18!n ph/s/m!u2!n/sr/nm]'
    cursor,x,y,/up
    if !mouse.button eq 4 then begin
    	wdelete,!window
	wdelete,!window
	goto,rep
    endif else begin
    	if !mouse.button eq 2 then begin
	    wdelete,!window
	    wdelete,!window
	    goto,fin
	endif else begin
    	    id_time = where(abs(time-y) eq min(abs(time-y)))
	    id_los  = where(abs(radial-x) eq min(abs(radial-x)))
	    oplot,[x,x],[0,1000],linest=5,col=colors.white
	    oplot,[0,1000],[y,y],linest=5,col=colors.white
    	    setgraphics,xs=800,ys=600,title='Left click to quit; right click to repeat'
	    plot,wave,emis[*,id_los[0],id_time[0]],/ylog,col=colors.black,back=colors.white,xtitle='Wavelength [nm]',ytitle='ph/s/m2/sr/nm'
    	    cursor,x,y,/up
	    if !mouse.button eq 4 then begin
    	    	wdelete,!window
    	    	wdelete,!window
    	    	wdelete,!window
 	    	goto,rep
	    endif
	end
    end
    fin:
Stop	
End

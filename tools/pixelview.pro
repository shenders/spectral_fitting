Pro pixelview,data,xr=xr


    wave = data.wvlngth[*,0] ; function of wave, los
    emis = data.emiss[*,*,*] ; function of wave, los, time
    time = data.time         ; time
    
    if keyword_set(xr)then begin
    	id_time = where(time ge xr[0] and time le xr[1])
	time = time[id_time]
	emis = emis[*,*,id_time]
    endif
    miss_pix = 0
    emis = emis[*,miss_pix:*,*]
    
    tint = total(emis,3)     ; sum dimension in time
    lint = total(tint,2)     ; sum dimension in los
    
    rep:
    setgraphics,xs=800,ys=600,colors=colors
    plot,wave,lint,/ylog,back=colors.white,col=colors.black
    cursor,x,y,/up
    
    id_pix = where(abs(wave-x) eq min(abs(wave-x)))
    intens = emis[id_pix[0],*,*]/1e18
    
    radial = fltarr(n_elements(emis[0,*,0]))
    for i=0,n_elements(radial)-1 do radial[i]=float(strmid(data.los[i],2,2))
    
    setgraphics,xs=800,ys=600,title='Left click for individual spectra; middle click to quit; right click to repeat'
    shimage,reform(intens),radial,time,xtitle='LOS',ytitle='Time [s]'
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
	    oplot,[x,x],[0,1000],linest=5,col=colors.white
	    oplot,[0,1000],[y,y],linest=5,col=colors.white
    	    setgraphics,xs=800,ys=600,title='Left click to quit; right click to repeat'
	    plot,wave,emis[*,x-miss_pix,id_time[0]],/ylog,col=colors.black,back=colors.white,xtitle='Wavelength [nm]',ytitle='ph/s/m2/sr/nm'
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

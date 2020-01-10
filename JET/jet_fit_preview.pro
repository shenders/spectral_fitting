Pro jet_fit_preview,shot,xr=xr,debug=debug,usesave=usesave
    
    rval = [2.52331,	 2.54117,      2.55904,	   2.57690,	2.59476,      2.61263,	  2.63049,      2.64836,$
            2.66622,	 2.68408,      2.70195,	   2.71981,	2.73767,      2.75554,	  2.77340,      2.79127,$
            2.80913,	 2.82699,      2.84486,	   2.86272,	2.88059,      2.89845]
    id_r = where(rval gt 2.55 and rval lt 2.86)
    rval = rval[id_r]
    shotstr = string(shot,format='(i5)')


    trace='save/'+shotstr+'/SP-data.idl'
    restore,trace[0]
    time    = output.time
    nii_3995 = output.nii[*,id_r,0]
    nii_4041 = output.nii[*,id_r,1]
    nii_4026 = output.nii[*,id_r,2]
    nii_3995_err = output.nii_err[*,id_r,0]
    nii_4041_err = output.nii_err[*,id_r,1]
    nii_4026_err = output.nii_err[*,id_r,2]

    if keyword_set(xr)then begin
    	    id = where(time ge xr[0] and time le xr[1])
    	    time = time[id]
    	    nii_3995 = nii_3995[id,*]
    	    nii_4041 = nii_4041[id,*]
    	    nii_4026 = nii_4026[id,*]
    	    nii_3995_err = nii_3995_err[id,*]
    	    nii_4041_err = nii_4041_err[id,*]
    	    nii_4026_err = nii_4026_err[id,*]
    endif
    
    exp_ratio1 = smooth(nii_4041,10,/edge_tr)/smooth(nii_3995,10,/edge_tr)
    exp_ratio2 = smooth(nii_4041,10,/edge_tr)/smooth(nii_4026,10,/edge_tr)
    num    = 40
    te_arr = [2.0,3.0,4.0,5.0,6.0]
    dens   = adas_vector(low=1.1e13,high=1e16,num=num)
    replot:
    setgraphics,colors=colors,colpick=colpick,/full,psplot=psplot,xs=800,ys=600,file='figures/ratios_db.ps'
    user_psym,1
    plot,[6,10],[0,0.5],/nodata,col=colors.black,back=colors.white,$
    xtitle='N II 404.1/402.6 [-]',ytitle='N II 404.1/399.5 [-]'
    for i=0,n_elements(te_arr)-1 do begin
	te = fltarr(num)+te_arr[i]
	if ~keyword_set(usesave)then begin
	    atomdb,te,dens,tec3995=tec3995,$
		           tec4041=tec4041,$
	       		   tec4026=tec4026
	    ratio1 = tec4041/tec3995
	    ratio2 = tec4041/tec4026
	    save,file='tmp/te_'+string(te[0],format='(d3.1)')+'.sav',ratio1,ratio2 
	endif else restore,'tmp/te_'+string(te[0],format='(d3.1)')+'.sav'
	oplot,ratio2,ratio1,col=colpick[i]
	
    endfor
    colpick = [colpick,colpick]
    if keyword_set(x0)then begin
    	id_time = where(time ge x0 and time le x1)
	id_spce = where(rval ge y0 and rval le y1)
    	for i=0,n_elements(id_spce)-1 do oplot,exp_ratio2[id_time,id_spce[i]],exp_ratio1[id_time,id_spce[i]],psym=8,col=colpick[i]	
    endif else begin
    	for i=0,n_elements(rval)-1 do oplot,exp_ratio2[*,i],exp_ratio1[*,i],psym=8,col=colpick[i]
    end
    
    setgraphics,xs=600,ys=400,colors=colors,/full,colpick=colpick
    rat=retrieve_ratio(nii_4041,output.los_names,shot)
    plot,time,rat,col=colors.black,back=colors.white,xtitle='Time [s]',ytitle='X/S-point ratio'

    if keyword_set(debug)then begin
	setgraphics,xs=600,ys=400,colors=colors,/full,colpick=colpick
    	shimage,nii_4041/1e17,time,rval,xtitle='Time [s]',ytitle='Radius [m]',/nobar
	if keyword_set(x0) then oplot,[x0,x0,x1,x1,x0],[y0,y1,y1,y0,y0],col=colors.white,linest=5
    	cursor,x0,y0,/up
    	cursor,x1,y1,/up
	oplot,[x0,x0,x1,x1,x0],[y0,y1,y1,y0,y0],col=colors.white,linest=5
	wdelete,!window
	wdelete,!window
	usesave=1
	goto,replot
    endif

Stop
End

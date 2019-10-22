Pro front,shot,xr=xr,sm=sm,psplot=psplot

	read_signal_mrm,0L,shot,'EVL','N_1_3995',time,nii_3995,2,exp=exp
	read_signal_mrm,0L,shot,'EVL','N_1_4041',time,nii_4041,2,exp=exp
	read_signal_mrm,0L,shot,'EVL','D_0_4101',time,di_4101,2,exp=exp
	read_signal_mrm,0L,shot,'EVL','N_2_4099',time,niii_4099,2,exp=exp
	read_signal_mrm,0L,shot,'MAC','Tdiv',tdiv_time,tdiv,2,exp=exp
	read_signal_mrm,0L,shot,'TOT','P_TOT',ptime,ptot,2,exp=exp
	
	if keyword_Set(psplot)then begin
		xs = 10 & ys = 7
	endif else begin
		xs = 1200 & ys = 900	
	end
	setgraphics,nrow=2,ncol=2,colpick=colpick,colors=colors,xs=xs,ys=ys,/full_list,psplot=psplot,filename='figures/shotdetails.ps',/landscape

	; Define channels

	channel = 9
	if shot eq 35844 then begin
		period = 0.505
		tstart = 2.1
		windows = [tstart,tstart+3*period,tstart+6*period]
	endif
	if shot eq 35848 then begin
		period = 0.505
		tstart = 2.6
		windows = [tstart,tstart+3*period,tstart+6*period]
	endif
	if shot eq 35846 or shot eq 35850 then begin
		period = 0.433
		tstart = 2.1
		windows = [tstart,tstart+3*period,tstart+6*period]	
	endif
	
	lowerte = 3.3
	upperte = 3.6
	
	nii_store = -1
	nii_rat_store = -1
	cn_store = -1
	cn_err_store = -1
	tdiv_store = -1
	power_store = -1
	for j=0,2 do begin
		if j eq 0 then plot,tdiv_time,tdiv,title='Tdiv',xtitle='Time [s]',ytitle='[eV]',xs=1,xr=xr,back=colors.white,/nodata,col=colors.black,yr=[0,40]
		if keyword_set(sm)then oplot,tdiv_time,smooth(tdiv,sm),col=colors.black else oplot,tdiv_time,tdiv,col=colors.black 
	endfor
	if shot eq 35844 or shot eq 35848 then begin
	for j=0,2 do begin
			t1 = windows[j]
			t2 = windows[j]+3*period
			bindata,tdiv_time,tdiv,[t1,t2],x1,x1_err
			oplot,[t1,t2],[x1,x1],col=colors.red,linest=5
			oplot,[t1,t1],[0,100],col=colors.black,linest=5
			oplot,[t2,t2],[0,100],col=colors.black,linest=5
			tdiv_store = [tdiv_store,x1]   
	endfor	
	endif
	if shot eq 35846 or shot eq 35850 then begin
	for j=0,2 do begin
		for i=0,2 do begin
			t1 = windows[j]+i*period
			t2 = windows[j]+i*period+period
			bindata,tdiv_time,tdiv,[t1,t2],x1,x1_err
			oplot,[t1,t2],[x1,x1],col=colors.red,linest=5
			oplot,[t1,t1],[0,100],col=colors.black,linest=5
			oplot,[t2,t2],[0,100],col=colors.black,linest=5
			tdiv_store = [tdiv_store,x1]   
		endfor	
	endfor	
	
	endif

	for j=0,2 do begin
		full_trange = [windows[j],windows[j]+3*period]
		if j eq 0 then plot,time,nii_3995[*,channel]/1e18,title='N II @ 399.5 nm',yr=[0,10],xs=1,xr=xr,xtitle='Time [s]',ytitle='[ph/s/m!u2!n/sr]',back=colors.white,/nodata,col=colors.black
		if keyword_set(sm)then oplot,time,smooth(nii_3995[*,channel]/1e18,sm),col=colors.black else oplot,time,nii_3995[*,channel]/1e18,col=colors.black
		
	endfor
	for j=0,2 do begin
		for i=0,2 do begin
			t1 = windows[j]+i*period
			t2 = windows[j]+i*period+period
			bindata,time,nii_3995[*,channel]/1e18,[t1,t2],x1,x1_err
			oplot,[t1,t2],[x1,x1],col=colors.red,linest=5
			oplot,[t1,t1],[0,100],col=colors.black,linest=5
			oplot,[t2,t2],[0,100],col=colors.black,linest=5
			nii_store = [nii_store,x1]   
		endfor
	endfor	
	if shot eq 35844 or shot eq 35848 then begin
	for j=0,2 do begin
		if j eq 0 then plot,ptime,ptot/1e6,xs=1,xr=xr,title='Total power',xtitle='Time [s]',ytitle='[MW]',back=colors.white,/nodata,col=colors.black
		if keyword_set(sm)then oplot,ptime,smooth(ptot/1e6,sm),col=colors.black else oplot,ptime,ptot/1e6,col=colors.black 
	endfor
	for j=0,2 do begin
		for i=0,2 do begin
			t1 = windows[j]+i*period
			t2 = windows[j]+i*period+period
			bindata,ptime,ptot/1e6,[t1,t2],x1,x1_err
			oplot,[t1,t2],[x1,x1],col=colors.red,linest=5
			oplot,[t1,t1],[0,100],col=colors.black,linest=5
			oplot,[t2,t2],[0,100],col=colors.black,linest=5
			power_store = [power_store,x1]   
		endfor
		;cursor,x,y,/down
	endfor
	endif
	if shot eq 35846 or shot eq 35850 then begin
	for j=0,2 do begin
		if j eq 0 then plot,ptime,ptot/1e6,xs=1,xr=xr,title='Total power',xtitle='Time [s]',ytitle='[MW]',back=colors.white,/nodata,col=colors.black
		if keyword_set(sm)then oplot,ptime,smooth(ptot/1e6,sm),col=colors.black else oplot,ptime,ptot/1e6,col=colors.black 
	endfor
	for j=0,2 do begin
			t1 = windows[j]
			t2 = windows[j]+3*period
			bindata,ptime,ptot/1e6,[t1,t2],x1,x1_err
			oplot,[t1,t2],[x1,x1],col=colors.red,linest=5
			oplot,[t1,t1],[0,100],col=colors.black,linest=5
			oplot,[t2,t2],[0,100],col=colors.black,linest=5
			power_store = [power_store,x1]   
		;cursor,x,y,/down
	endfor
	endif
	
	for j=0,2 do begin
		full_trange = [windows[j],windows[j]+3*period]
		if j eq 0 then plot,time,nii_4041[*,channel]/nii_3995[*,channel],title='N II 4041/3995 ratio',yr=[0,0.4],xs=1,xr=xr,xtitle='Time [s]',ytitle='[-]',back=colors.white,/nodata,col=colors.black
		if keyword_set(sm)then oplot,time,smooth(nii_4041[*,channel],sm)/smooth(nii_3995[*,channel],sm),col=colors.black else oplot,time,nii_4041[*,channel]/nii_3995[*,channel],col=colors.black
		
	endfor
	for j=0,2 do begin
		for i=0,2 do begin
			t1 = windows[j]+i*period
			t2 = windows[j]+i*period+period
			bindata,time,(nii_4041[*,channel]/nii_3995[*,channel])<0.4,[t1,t2],x1,x1_err
			oplot,[t1,t2],[x1,x1],col=colors.red,linest=5
			oplot,[t1,t1],[0,100],col=colors.black,linest=5
			oplot,[t2,t2],[0,100],col=colors.black,linest=5
			nii_rat_store = [nii_rat_store,x1] 
		endfor
	endfor	
	setgraphics,psplot=psplot,/close
	
	analysis,shot,/use_evl,channel=channel,xr=xr,cn_mean=cn_mean,cn_err=cn_err,time=cn_time,lowerte=lowerte,upperte=upperte

        for j=0,2 do begin
	       for i=0,2 do begin
		       t1 = windows[j]+i*period
		       t2 = windows[j]+i*period+period
		       bindata,cn_time,cn_mean*100,[t1,t2],x1,x1_err
		       cn_store = [cn_store,x1]   
		       bindata,cn_time,cn_err*100,[t1,t2],x1,x1_err
		       cn_err_store = [cn_err_store,x1]   
	       endfor
        endfor  
	nii_store   = nii_store[1:*]  
	nii_rat_store= nii_rat_store[1:*]  
	cn_store    = cn_store[1:*]    
	cn_err_store= cn_err_store[1:*]    
	tdiv_store  = tdiv_store[1:*]  
	power_store = power_store[1:*]  
	deltaL,findgen(100)*20/99.0,dl_used,/rov014
	deltaL,tdiv_store,dl,/rov014
	if shot eq 35844 or shot eq 35848 then begin
	if keyword_Set(psplot)then begin
		xs = 8 & ys = 10
	endif else begin
		xs = 800 & ys = 1000	
	end
	setgraphics,nrow=1,ncol=2,xs=xs,ys=ys,colpick=colpick,psplot=psplot,filename='figures/concentrations.ps',/full,/portrait
	user_psym,1,/fill
	if shot eq 35844 then jrange = [0,2] else jrange = [1,1]
	for j=jrange[0],jrange[1] do begin
		t1 = [3*j,3*(j+1)-1]
		if j eq jrange[0] then plot,power_store[t1[0]:t1[1]],cn_store[t1[0]:t1[1]],xtitle='Power [MW]',ytitle='cN [%]',yr=[0,6],back=colors.white,/nodata,col=colors.black
		oplot,power_store[t1[0]:t1[1]],cn_store[t1[0]:t1[1]],col=colpick[j],psym=-8,symsize=2.0
		err_plot,power_store[t1[0]:t1[1]],cn_store[t1[0]:t1[1]],cn_err_store[t1[0]:t1[1]],col=colpick[j]
	endfor
	for j=jrange[0],jrange[1] do begin
		t1 = 3*j
		p0 = power_store[t1]
		p1 = power_store[t1+1]
		p2 = power_store[t1+2]
		c0 = cn_store[t1]
		c1 = cn_store[t1+1]
		c2 = cn_store[t1+2]
		c0e = cn_err_store[t1]
		c1e = cn_err_store[t1+1]
		c2e = cn_err_store[t1+2]
		rat0  = c0/c0
		rat0e = rat0 * sqrt((c0e/c0)^2+(c0e/c0)^2)
		rat1  = c1/c0
		rat1e = rat1 * sqrt((c1e/c1)^2+(c0e/c0)^2)
		rat2  = c2/c0
		rat2e = rat2 * sqrt((c2e/c2)^2+(c0e/c0)^2)
		print,rat0,rat0e
		print,rat1,rat1e
		print,rat2,rat2e 
		if j eq jrange[0] then plot,[p0/p0,p1/p0,p2/p0]*100-100,[rat0,rat1,rat2]*100-100,xtitle='Power increase [%]',ytitle='cN increase [%]',xr=[0,100],yr=[0,100],back=colors.white,/nodata,col=colors.black
		oplot,[p0/p0,p1/p0,p2/p0]*100-100,[rat0,rat1,rat2]*100-100,col=colpick[j],psym=-8,symsize=2.0
		;err_plot,[p0/p0,p1/p0,p2/p0]*100-100,[rat0,rat1,rat2]*100-100,[rat0e,rat1e,rat2e]*100,col=colpick[j]
	endfor
	oplot,[0,100],[0,100],linest=5,col=colors.black
		rat1  = c1/c0
		rat1e = rat1 * sqrt((c1e/c1)^2+(c0e/c0)^2)
	setgraphics,psplot=psplot,/close
	endif

	if shot eq 35846 or shot eq 35850 then begin
	if keyword_Set(psplot)then begin
		xs = 8 & ys = 6
	endif else begin
		xs = 800 & ys = 600	
	end
	setgraphics,nrow=1,ncol=1,xs=xs,ys=ys,colpick=colpick,psplot=psplot,filename='figures/concentrations.ps',/full,/portrait
	user_psym,1,/fill
	for j=0,2 do begin
		t1 = [3*j,3*(j+1)-1]
		if j eq 0 then plot,tdiv_store[t1[0]:t1[1]],cn_store[t1[0]:t1[1]],xtitle='Tdiv [eV]',ytitle='cN [%]',xr=[0,30],yr=[0,10],back=colors.white,/nodata,col=colors.black
		oplot,tdiv_store[t1[0]:t1[1]],cn_store[t1[0]:t1[1]],col=colpick[j],psym=8,symsize=2.0
		err_plot,tdiv_store[t1[0]:t1[1]],cn_store[t1[0]:t1[1]],cn_err_store[t1[0]:t1[1]],col=colpick[j]
	endfor
	
	setgraphics,psplot=psplot,/close
	endif
	
Stop
End

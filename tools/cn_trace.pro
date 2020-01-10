Pro cn_trace,shot,$
             los=los,$
	     append=append,$
	     xr=xr,$
	     use_evl=use_evl,$
	     diag=diag,$
	     channel=channel,$
	     interelm=interelm,$
	     transmission=transmission,$
	     full=full,$
	     sm=sm



	data  =  get_nii(shot,$
		 los=los,$
		 xr=xr,$
		 interelm=interelm,$
 		 channel=channel,$
		 use_evl=use_evl,$
		 append=append,$
		 diag=diag)
		 
	if ~keyword_set(transmission)then transmission=1.0
	if ~keyword_set(sm)then sm=20.0
 
       	; get cn, Te and dens from line ratios

        if keyword_set(full)then begin
		plasma = optimize_ratios(data,shot,los,transmission,sm=sm,debug=debug)
       	endif else begin
	   	plasma = estimate_ratios(data,shot,los,transmission,sm=sm,debug=debug,lowerte=3.1,upperte=3.5)
       	end
	
	nrow = 2
	ncol = 1
	xspc = 0.1
	yspc = 0.1
	setgraphics,nrow=nrow,ncol=ncol,colors=colors,xs=800,ys=900
	
	plot,data.time,plasma.te_upper,col=colors.black,back=colors.white,$
	xr=xr,xs=1,ys=9,yr=[0,6],ytitle='T!le !n[eV]',xtitle='Time [s]',pos=graphpos(0,1,nrow,ncol,xspc=xspc,yspc=yspc)
	
	oplot,data.time,plasma.te_lower,col=colors.black
	axis,yaxis=1,/save,yr=[0,3],ytitle='n!le !n[10!u20!n m!u-3!n]',col=colors.blue
	oplot,data.time,plasma.dens_upper/1e14,col=colors.blue
	oplot,data.time,plasma.dens_lower/1e14,col=colors.blue
	
	plot,data.time,plasma.cn_upper*100,col=colors.black,back=colors.white,$
	xr=xr,xs=1,ys=1,yr=[0,15],ytitle='c!lN !n[%]',xtitle='Time [s]',pos=graphpos(0,0,nrow,ncol,xspc=xspc,yspc=yspc)
	oplot,data.time,plasma.cn_lower*100,col=colors.black

	stop
End

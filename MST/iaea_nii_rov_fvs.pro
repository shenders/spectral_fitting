Pro legend,text,xr,yr,col
	y0 = yr[0]+(yr[1] - yr[0])*0.8
	x0 = xr[0]+(xr[1] - xr[0])*0.7
	xyouts,x0,y0,text,col=col,charsize=1.4
End
PRO iaea_nii_rov_fvs,psplot=psplot,shot=shot,rawelms=rawelms,xr=xr,yr=yr,elmcond=elmcond,dynamic=dynamic,use_rov8=use_rov8,interelm=interelm

	if ~keyword_Set(shot)then shot=30776
	time=-1
	read_signal_mrm,0L,shot,'MAC','Tdiv',time,tdiv,2,exp=exp
	;restore,'save/'+string(shot,format='(i5)')+'/ROV008-data.idl'
	;n08  = output.nii[*,0,0]/1e18
	;t08  = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV-10-data.idl'
	n10  = output.nii[*,0,0]/1e17
	t10  = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV-12-data.idl'
	n12  = output.nii[*,0,0]/1e17
	t12  = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV-14-data.idl'
	n14  = output.nii[*,0,0]/1e17
	t14  = output.time	
	n08 = n10
	t08 = t10
	if ~keyword_set(elmcond)then elmcond = 6.5
	if keyword_set(rawelms)then telm    = find_elms(shot,t10,use_rov8=use_rov8,dynamic=dynamic) else telm= find_elm(shot,t10)
	if keyword_set(interelm)then id      = where(telm ge elmcond) else id=findgen(n_elements(t10))
	xspc    = 0.05
	yspc    = 0.00
	
	nrow    = 5
	ncol    = 1
	if ~keyword_set(xr)then xr=[min(t10),max(t10)]
	setgraphics,nrow=nrow,ncol=ncol,xs=900,ys=800,colors=colors,psplot=psplot,file='strikepoint_scan.ps',/portrait
	plot,time,tdiv,col=colors.black,backg=colors.white,position=graphpos(0,4,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,xs=1,yr=[-10,30],$
	xtitle='Time [s]',ytitle='[eV]' & legend,'(e) Tdiv',xr,[-10,30],colors.black
	
	yr = [0,40]
	y10 = interpol(n10[id],t10[id],t10)
	y12 = interpol(n12[id],t12[id],t10)
	y14 = interpol(n14[id],t14[id],t10)
	sv = 15
	plot,t08,smooth(y10,sv,/edge_trun),col=colors.black,xs=9,backg=colors.white,yr=yr,position=graphpos(0,3,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	oplot,t08,smooth(y12,sv,/edge_trun),col=colors.blue 
	oplot,t08,smooth(y14,sv,/edge_trun),col=colors.red 
	legend,'(d) N II sightlines',xr,yr,colors.black
	yr = [0,40]
	user_psym,5 & plot,t10,n10,psym=8,col=colors.black,xs=9,backg=colors.white,yr=yr,position=graphpos(0,2,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,t10[id],n10[id],psym=8,col=colors.black & legend,'(c) ROV010',xr,yr,colors.black

	yr = [0,40]
	user_psym,5 & plot,t12,n12,psym=8,col=colors.black,xs=9,backg=colors.white,yr=yr,position=graphpos(0,1,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,t12[id],n12[id],psym=8,col=colors.blue & legend,'(b) ROV012',xr,yr,colors.black

	yr = [0,30]
	user_psym,5 & plot,t14,n14,psym=8,col=colors.black,xs=1,backg=colors.white,$
	yr=yr,position=graphpos(0,0,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,t14[id],n14[id],psym=8,col=colors.red & legend,'(a) ROV014',xr,yr,colors.black

	if keyword_set(psplot)then setgraphics,/close
STOP
End

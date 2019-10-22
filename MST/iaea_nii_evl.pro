PRO iaea_nii_evl,psplot=psplot,shot=shot,rawelms=rawelms,xr=xr,yr=yr,elmcond=elmcond,dynamic=dynamic,use_rov8=use_rov8,$
                 channels=channels,diag=diag,interelm=interelm

	if ~keyword_Set(shot)then shot=30776
	time=-1
	read_signal_mrm,0L,shot,'MAC','Tdiv',time,tdiv,2,exp=exp

	read_signal_mrm,0L,shot,diag[0],'N_1_3995',t10,n10,2,exp=exp
	read_signal_mrm,0L,shot,diag[1],'N_1_3995',t12,n12,2,exp=exp
	read_signal_mrm,0L,shot,diag[2],'N_1_3995',t14,n14,2,exp=exp


	if ~keyword_set(xr)then xr=[2,6]

	id = where(t10 ge xr[0] and t10 le xr[1])
	
	n10  = n10[id,channels[0]]/1e17
	t10  = t10[id]
	n12  = n12[id,channels[1]]/1e17
	t12  = t12[id]
	n14  = n14[id,channels[2]]/1e17
	t14  = t14[id]



	if ~keyword_set(elmcond)then elmcond = 3.0
	if keyword_set(rawelms)then telm    = find_elms(shot,t10,use_rov8=use_rov8,dynamic=dynamic) else telm= find_elm(shot,t10)
	if keyword_set(interelm)then id      = where(telm ge elmcond) else id=findgen(n_elements(t14))
	xspc    = 0.05
	yspc    = 0.00
	
	nrow    = 5
	ncol    = 1
	if ~keyword_set(xr)then xr=[min(t10),max(t10)]
	setgraphics,nrow=nrow,ncol=ncol,xs=900,ys=800,colors=colors,psplot=psplot,file='strikepoint_scan.ps',/portrait
	plot,time,tdiv,col=colors.black,backg=colors.white,position=graphpos(0,4,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,xs=1,yr=[0,30],$
	xtitle='Time [s]',ytitle='[eV]' & legend,'(e) Tdiv',colors.black
	
	
	y10 = interpol(n10[id],t10[id],t10)
	y12 = interpol(n12[id],t12[id],t10)
	y14 = interpol(n14[id],t14[id],t10)
	sv = 15
	plot,t10,smooth(y10,sv,/edge_trun),col=colors.black,xs=9,backg=colors.white,yr=[0,max(n14)*1.5],position=graphpos(0,3,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,$
	ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	oplot,t10,smooth(y12,sv,/edge_trun),col=colors.blue 
	oplot,t10,smooth(y14,sv,/edge_trun),col=colors.red 
	legend,'(d) N II sightlines',colors.black

	user_psym,5 & plot,t10,n10,psym=8,col=colors.black,xs=9,backg=colors.white,yr=[0,max(n10)*1.5],position=graphpos(0,2,nrow,ncol,xspc=xspc,yspc=yspc),$
	xr=xr,ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,t10[id],n10[id],psym=8,col=colors.black & legend,'(c) ROV010',colors.black

	
	user_psym,5 & plot,t12,n12,psym=8,col=colors.black,xs=9,backg=colors.white,yr=[0,max(n12)*1.5],position=graphpos(0,1,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,t12[id],n12[id],psym=8,col=colors.blue & legend,'(b) ROV012',colors.black

	user_psym,5 & plot,t14,n14,psym=8,col=colors.black,backg=colors.white,xs=1,$
	yr=[0,max(n14)*1.5],position=graphpos(0,0,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,t14[id],n14[id],psym=8,col=colors.red & legend,'(a) ROV014',colors.black

	if keyword_set(psplot)then setgraphics,/close
STOP
End

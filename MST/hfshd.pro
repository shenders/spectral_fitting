Pro HFSHD,shot,xr=xr


	read_signal_mrm,0L,shot,'EVL','Ne',time,ne_evl,2,exp=exp

	ziv4 = ne_evl[*,14]
	ziv6 = ne_evl[*,15]
	
	zon1 = ne_evl[*,17]
	zon5 = ne_evl[*,19]
	
	rtn3 = ne_evl[*,0]
	rtn5 = ne_evl[*,1]
	
	adas_colors,colors=colors
	window,0,xs=1000,ys=1000
	!p.multi=[0,1,3]
	!p.charsize=3.0
	
	plot,time,ziv4,xr=xr,xs=1,yr=[0,2e21],col=0,back=255
	oplot,time,ziv6,col=colors.red
	plot,time,rtn3,xr=xr,xs=1,yr=[0,2e21],col=0,back=255
	oplot,time,rtn5,col=colors.red
	plot,time,zon1,xr=xr,xs=1,yr=[0,2e21],col=0,back=255
	oplot,time,zon5,col=colors.red
stop
End

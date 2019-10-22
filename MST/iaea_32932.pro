Pro legend,text,xr,yr,col
	y0 = yr[0]+(yr[1] - yr[0])*0.8
	x0 = xr[0]+(xr[1] - xr[0])*0.7
	xyouts,x0,y0,text,col=col,charsize=1.4
End
PRO iaea_32932,psplot=psplot

	shot=32932
	time=-1
	read_signal_mrm,0L,shot,'MAC','Tdiv',time,tdiv,2,exp=exp
	restore,'save/'+string(shot,format='(i5)')+'/ROV008-data.idl'
	n08  = output.nii[*,0,0]/1e18
	t08  = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV010-data.idl'
	n10  = output.nii[*,0,0]/1e18
	t10  = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV012-data.idl'
	n12  = output.nii[*,0,0]/1e18
	t12  = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV014-data.idl'
	n14  = output.nii[*,0,0]/1e18
	t14  = output.time	
	
	elmcond = 4.5
	telm    = find_elm(shot,t08)
	id      = where(telm ge elmcond)
	xspc    = 0.05
	yspc    = 0.00
	
	xr      = [0,30]
	yr      = [0,8]
	
	nrow    = 4
	ncol    = 2
	td = interpol(tdiv,time,t08)
	setgraphics,nrow=nrow,ncol=ncol,xs=900,ys=800,colors=colors,psplot=psplot,file='strikepoint_scan.ps',/portrait
	user_psym,5 & plot,td,n08,psym=8,col=colors.black,xs=8,backg=colors.white,position=graphpos(0,3,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,yr=[0,12],$
	ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtitle='Tdiv [eV]'
	user_psym,5,/fill & oplot,td[id],n08[id],psym=8,col=colors.red & legend,'(d) ROV008',xr,[0,12],colors.black
	x = td[id] & y=n08[id] & id_sort = sort(x) & x = x[id_sort] & y = y[id_sort]
	xx = svdfit(x,y,3,yfit=yfit) & x01 = x & y01=yfit
	oplot,x,yfit,col=colors.green,thick=2

	user_psym,5 & plot,td,n10,psym=8,col=colors.black,xs=8,backg=colors.white,position=graphpos(0,2,nrow,ncol,xspc=xspc,yspc=yspc),$
	xr=xr,yr=[0,8],ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,td[id],n10[id],psym=8,col=colors.red & legend,'(c) ROV010',xr,[0,8],colors.black
	x = td[id] & y=n10[id] & id_sort = sort(x) & x = x[id_sort] & y = y[id_sort]
	xx = svdfit(x,y,4,yfit=yfit) & x11 = x & y11=yfit
	oplot,x,yfit,col=colors.green,thick=2

	user_psym,5 & plot,td,n12,psym=8,col=colors.black,xs=8,backg=colors.white,position=graphpos(0,1,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,yr=[0,8],ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,td[id],n12[id],psym=8,col=colors.red & legend,'(b) ROV012',xr,[0,8],colors.black
	x = td[id] & y=n12[id] & id_sort = sort(x) & x = x[id_sort] & y = y[id_sort]
	xx = svdfit(x,y,3,yfit=yfit) & x21 = x & y21=yfit
	oplot,x,yfit,col=colors.green,thick=2

	user_psym,5 & plot,td,n14,psym=8,col=colors.black,xs=0,backg=colors.white,position=graphpos(0,0,nrow,ncol,xspc=xspc,yspc=yspc),title='AUG #32932',xr=xr,yr=[0,5],ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,td[id],n14[id],psym=8,col=colors.red & legend,'(a) ROV014',xr,[0,5],colors.black
	x = td[id] & y=n14[id] & id_sort = sort(x) & x = x[id_sort] & y = y[id_sort]
	xx = svdfit(x,y,3,yfit=yfit) & x31 = x & y31=yfit
	oplot,x,yfit,col=colors.green,thick=2


	shot=30623
	time=-1
	read_signal_mrm,0L,shot,'MAC','Tdiv',time,tdiv,2,exp=exp
	restore,'save/'+string(shot,format='(i5)')+'/ROV008-data.idl'
	n08  = output.nii[*,0,0]/1e18
	t08  = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV010-data.idl'
	n10  = output.nii[*,0,0]/1e18
	t10  = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV012-data.idl'
	n12  = output.nii[*,0,0]/1e18
	t12  = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV014-data.idl'
	n14  = output.nii[*,0,0]/1e18
	t14  = output.time	
	telm    = find_elm(shot,t08)
	id      = where(telm ge elmcond)
	td = interpol(tdiv,time,t08)

	user_psym,5 & plot,td,n08,psym=8,col=colors.black,xs=8,backg=colors.white,position=graphpos(1,3,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,yr=[0,60],$
	ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtitle='Tdiv [eV]'
	user_psym,5,/fill & oplot,td[id],n08[id],psym=8,col=colors.red & legend,'(i) ROV008',xr,[0,60],colors.black
	x = td[id] & y=n08[id] & id_sort = sort(x) & x = x[id_sort] & y = y[id_sort]
	xx = svdfit(x,y,3,yfit=yfit) & x02 = x & y02=yfit
	oplot,x,yfit,col=colors.green,thick=2

	user_psym,5 & plot,td,n10,psym=8,col=colors.black,xs=8,backg=colors.white,position=graphpos(1,2,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,yr=[0,30],ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,td[id],n10[id],psym=8,col=colors.red & legend,'(h) ROV010',xr,[0,30],colors.black
	x = td[id] & y=n10[id] & id_sort = sort(x) & x = x[id_sort] & y = y[id_sort]
	xx = svdfit(x,y,4,yfit=yfit) & x12 = x & y12=yfit
	oplot,x,yfit,col=colors.green,thick=2

	user_psym,5 & plot,td,n12,psym=8,col=colors.black,xs=8,backg=colors.white,position=graphpos(1,1,nrow,ncol,xspc=xspc,yspc=yspc),xr=xr,yr=[0,40],ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,td[id],n12[id],psym=8,col=colors.red & legend,'(g) ROV012',xr,[0,40],colors.black
	x = td[id] & y=n12[id] & id_sort = sort(x) & x = x[id_sort] & y = y[id_sort]
	xx = svdfit(x,y,3,yfit=yfit) & x22 = x & y22=yfit
	oplot,x,yfit,col=colors.green,thick=2

	user_psym,5 & plot,td,n14,psym=8,col=colors.black,xs=0,backg=colors.white,position=graphpos(1,0,nrow,ncol,xspc=xspc,yspc=yspc),title='AUG #30623',xr=xr,yr=[0,20],ytitle='[10!u18!n ph/s/m!u2!n/sr]',xtickname=replicate(' ',30)
	user_psym,5,/fill & oplot,td[id],n14[id],psym=8,col=colors.red & legend,'(f) ROV014',xr,[0,20],colors.black
	x = td[id] & y=n14[id] & id_sort = sort(x) & x = x[id_sort] & y = y[id_sort]
	xx = svdfit(x,y,3,yfit=yfit) & x32 = x & y32=yfit
	oplot,x,yfit,col=colors.green,thick=2
	
	
	
	setgraphics,/close
STOP
End

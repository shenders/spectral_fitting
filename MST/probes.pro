Pro gfunct, X, A, F, pder
; Eich function 
fx  = 4.5
q0  = a[0]
S   = a[1]
lq  = 3.7e-3;a[2]
s0  = a[2]
qBG = 0.25
sbar= x-s0
f1  = S/2.0/lq/fx
f2  = sbar/lq/fx
if s < 0 then f = fltarr(n_elements(sbar)) - 1.0 else F   = (q0/2.0) * exp(f1^2 - f2) * erfc( f1 - f2 ) + qBG

End

Function binner,time,data,dt
	
	bin     = 0.0
	iter    = 0.0
	tsum    = 0.0
	t0      = time[0]
	bindata = -1.0
	bintime = -1.0
	for i=0,n_elements(time)-1 do begin
		tsum    = time[i] - t0
		if tsum gt dt then begin
			bindata = [bindata,bin/iter]
			bintime = [bintime,(t0+time[i-1])/2.0]
			bin     = 0.0
			iter    = 0.0
			t0      = time[i] 
		endif else begin
			bin     = bin + data[i]
			iter    = iter + 1.0
		end	
	endfor
	bindata = bindata[1:*]
	bintime = bintime[1:*]
	avrdata = mean(bindata[where(finite(bindata))])
	return,{data:bindata,time:bintime,avr:avrdata}
End
Function probe_map,shot

	if shot eq 30505 then begin
		; Mapping function	
		number  = [ 2     , 3     , 1     , 5     , 6     , 2     , 6     , 2     , 4     , 5     , 2     , 4     , 6    ]
		channel = [ 21    , 25    , 25    , 24    , 21    , 22    , 22    , 23    , 23    , 23    , 24    , 24    , 24   ]
		names   = [ 'ua1' , 'res1', 'res2', 'res3', 'ua3' , 'ua4' , 'ua6' , 'ua7' , 'ua8' , 'ua9' , 'uaa' , 'uab' , 'uac']
		spos    = [ 0.996 , 0.996 , 1.021 , 1.046 , 1.046 , 1.071 , 1.121 , 1.146 , 1.171 , 1.221 , 1.271 , 1.302 , 1.331]
		ds      = (spos-1.0)*100 ; in cm
		area    = 0.005 * 0.025 +fltarr(n_elements(number))
	endif
	
	return,{number:number,$
	        channel:channel,$
		ds:ds,$
		area:area,$
		names:names}
End
PRO rollover,shot,names,trange=trange,psplot=psplot,dt=dt,elmcond=elmcond

	; timing
	if ~keyword_set(trange)then trange=[2.0,5.8]
	if ~keyword_set(dt)then dt= 0.0025
	if ~keyword_set(elmcond)then elmcond = 5.0

	; plotting
	nrow=2
	ncol=1
	xspc    = 0.05
	yspc    = 0.00
	setgraphics,nrow=nrow,ncol=ncol,xs=800,ys=800,colors=colors,psplot=psplot,file='figures/rollover.ps',/portrait,colpick=colpick

	; get details
	
	details = probe_map(shot)
	
	colpick = [colors.black,colors.blue,colors.red,colors.cyan,colors.orange,colors.aqua]
	for i=0,n_elements(names)-1 do begin
		id = where(details.names eq names[i])
		print,'Fetching '+details.names[id[0]]
		diag='CH'+string(details.channel[id[0]],format='(i2)')
		time = -1 & trace = -1
		read_signal_mrm,0L,shot,'LSF',diag,time,trace,2,exp=exp
		x       = time
		y       = trace[*,details.number[id[0]]-1]/details.area[id[0]]/1.6e-19/1e23 
		telm    = find_elm(shot,x)
		idelm   = where(telm ge elmcond)
		bins    = binner(x[idelm],y[idelm],dt)
		if i eq 0 then begin
			plot,bins.time,bins.data,col=colors.black,back=colors.white,yr=[0,max(bins.data)*1.2],$
		  	ytitle='J!lsat!n [10!u23!n m!u-2!ns!u-1!n]',xr=trange,/nodata,$
			position=graphpos(0,0,nrow,ncol,xspc=xspc,yspc=yspc),xs=1,xtickname=replicate(' ',30)
			legend,string(shot,format='("#AUG ",i5)'),colors.black,yshift=0.1,xshift=-0.65
		endif
		oplot,bins.time,bins.data,col=colpick[i]
		legend,string(details.ds[id[0]]-3,format='("dS = ",f5.1," [cm]")'),colpick[i],yshift=-i*0.1,xshift=-0.65
	endfor
	x = -1 & y=-1
	read_signal_mrm,0L,shot,'MAC','Tdiv',x,y,2,exp=exp
	plot,x,y,xr=trange,yr=[-10,30],col=colors.black,back=colors.white,xtitle='Time [s]',ytitle='Tdiv [eV]',$
			position=graphpos(0,1,nrow,ncol,xspc=xspc,yspc=yspc),xs=1
	oplot,[-10,100],[0,0],linest=5,col=colors.black
End	

PRO heatflux,shot,trange,ch=ch,num=num,q0=q0,lq=lq,spread=spread,psplot=psplot,dt=dt,elmcond=elmcond

	; timing
	
	if ~keyword_set(dt)then dt= 0.0025
	if ~keyword_set(elmcond)then elmcond = 5.0

	; plotting
	
	setgraphics,nrow=nrow,ncol=ncol,xs=500,ys=500,colors=colors,psplot=psplot,file='figures/'+string(trange[0]*1000,trange[1]*1000,format='(i4,"-",i4,"-strikepoint_scan.ps")'),/portrait
	
	; get details
	
	details = probe_map(shot)	
	nprobes = n_elements(details.number)
	bindata = -1
	user_psym,5,/fill
	
	
	for i=0,nprobes-1 do begin
		print,'Fetching '+details.names[i]
		diag='CH'+string(details.channel[i],format='(i2)')
		time = -1 & trace = -1
	      	read_signal_mrm,0L,shot,'LSF',diag,time,trace,2,exp=exp
		id   = where(time ge trange[0] and time le trange[1])
		if id[0] ne -1 then begin
			x       = time[id]
			y       = trace[id,details.number[i]-1]/details.area[i]/1.6e-19/1e23 
			telm    = find_elm(shot,x)
			idelm   = where(telm ge elmcond)
			bins    = binner(x[idelm],y[idelm],dt)
			z       = details.ds[i] + fltarr(n_elements(bins.data))
			bindata = [bindata,bins.avr]
			if i eq 0 then plot,z,bins.data,back=colors.white,col=colors.black,/nodata,yr=[0,5],xr=[-5,15],xs=1,$
			      ytitle='J!lsat!n [10!u23!n m!u-2!ns!u-1!n]',xtitle='DS [cm]'			
			oplot,z-3,bins.data,col=colors.black,psym=8
		endif
	endfor
	
	bindata = bindata[1:*]
	;oplot,details.ds-3,bindata,col=colors.red
	weights = 1/bindata
	a = [10,2.79E-3,0.00]
	yfit = curvefit((details.ds-3.0)/100.0,bindata,weights,a,sigma,function_name='gfunct',/noderivative)
	xx = findgen(1000)*40/999.0-10.0
        gfunct, XX/100.0, A, yy, pder
	oplot,xx,yy,col=colors.red
	lq = a[2]
	spread = a[1]
	q0 = a[0]
End
Pro probes,psplot=psplot
	!quiet=1
	dt= 0.01
	elmcond = 5.0
	rollover,30505, ['res2','ua3','ua4'],psplot=psplot,dt=dt,elmcond=elmcond
	heatflux,30505, [3.8,4.1],q0=q01,lq=lq1,spread=s1,psplot=psplot,dt=dt,elmcond=elmcond
	heatflux,30505, [4.4,4.7],q0=q02,lq=lq2,spread=s2,psplot=psplot,dt=dt,elmcond=elmcond
	heatflux,30505, [4.7,5.0],q0=q03,lq=lq3,spread=s3,psplot=psplot,dt=dt,elmcond=elmcond
	print,string(lq1,s1,q01,format='("Lambda_q = ",E9.2,"; S = ",E9.2,"; q0 = ",E9.2)')
	print,string(lq2,s2,q02,format='("Lambda_q = ",E9.2,"; S = ",E9.2,"; q0 = ",E9.2)')
	print,string(lq3,s3,q03,format='("Lambda_q = ",E9.2,"; S = ",E9.2,"; q0 = ",E9.2)')

End

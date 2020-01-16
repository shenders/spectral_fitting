Function h98,shot
    ppfread,shot=shot,dda='KG1L',dtype='LAD3',data=N,t=t_N
    ppfread,shot=shot,dda='EFIT',dtype='WDIA',data=Wmhd,t=t_wmhd
    ppfread,shot=shot,dda='EFIT',dtype='ELON',data=kappa,t=t_kappa
    ppfread,shot=shot,dda='EFIT',dtype='POHM',data=Pohm,t=t_pohm
    ppfread,shot=shot,dda='EFIT',dtype='RMAG',data=Rmaj,t=t_rmaj
    ppfread,shot=shot,dda='NBI',dtype='NBLM',data=Pnbi,t=t_pnbi 
    ppfread,shot=shot,dda='ICRH',dtype='PTOT',data=Picrh,t=t_picrh
    ppfread,shot=shot,dda='MAGN',dtype='IPLA',data=Ip,t=t_ip 
    ppfread,shot=shot,dda='MAGN',dtype='BVAC',data=Bvac,t=t_bvac 
         

    Return,hfac
End

PRO jet_control_room,pulse,tr=tr,average=average
    !quiet=1
    if ~keyword_set(tr)then tr=[40,60]
    for i=0,n_elements(pulse)-1 do begin
    	shot=pulse[i]
    	if keyword_set(average)then begin
	    print,' '
	    print,'-----------------------'
	    print,'Average values for t=',string(average-0.2,average+0.2,shot,format='(D5.2,"-",D5.2, " in JPN",I5)')
	endif
    	setgraphics,xs=1000,ys=1000,psplot=psplot,file='jet_trace.ps',colors=colors,nrow=3,ncol=3,title=string(shot,format='(i5)')
; Total input and radiated power
    	ppfread,shot=shot,dda='NBI',dtype='NBLM',data=data,t=t & plot,t,data/1e6,col=colors.black,back=colors.white,xtit='Time [s]',ytit='Input power [MW]',yr=[0,30],/nodata,xr=tr,xs=1 & oplot,t,data/1e6,col=colors.blue
    	ppfread,shot=shot,dda='ICRH',dtype='PTOT',data=data2,t=t2 & oplot,t2,data2/1e6,col=colors.red
    	ppfread,shot=shot,dda='EFIT',dtype='POHM',data=Pohm,t=t_pohm
    	ppfread,shot=shot,dda='SCAL',dtype='PSOL',data=data3,t=t3 & if keyword_set(t3)then oplot,t3,data3/1e6,col=colors.green
    	pin = (data + interpol(data2,t2,t) + interpol(Pohm,t_pohm,t))/1e6
	oplot,t,pin,col=colors.black
    	if keyword_set(average)then begin
	    id=where(t ge average-0.2 and t le average+0.2)
	    print,'Total input power (MW): ',mean(pin[id])
	endif
    	if keyword_set(average)then begin
	    id=where(t3 ge average-0.2 and t3 le average+0.2)
	    print,'Psep (MW): ',mean(data3[id]/1e6)
	endif
    	
; Radiated power   
    	ppfread,shot=shot,dda='BOLO',dtype='TOPI',data=data,t=t & plot,t,data/1e6>0,col=colors.black,back=colors.white,xtit='Time [s]',ytit='Total radiated power [MW]',xr=tr,xs=1,yr=[0,mean(data/1e6)*5.0]
; Gas valves
    	ppfread,shot=shot,dda='GASM',dtype='MAJR',data=data,t=t & plot,t,data/1e21,col=colors.black,back=colors.white,xtit='Time [s]',ytit='Gas valve, cZ [10!u21!n #/s, %]',/nodata,xr=tr,xs=1 & oplot,t,data/1e21,col=colors.red
    	ppfread,shot=shot,dda='GASM',dtype='MN1R',data=data2,t=t2 & oplot,t2,data2/1e21,col=colors.blue
    	cn = (interpol(data2/7,t2,t) / (data+interpol(data2/7,t2,t)))*100
	;oplot,t,cn,col=colors.black
    	if keyword_set(average)then begin
	    id=where(t ge average-0.2 and t le average+0.2)
	    print,'Valve flux ratio (cN) [%]: ',mean(cn[id])
	endif
; Diamagnetic energy
    	ppfread,shot=shot,dda='EFIT',dtype='WDIA',data=data,t=t & plot,t,data/1e6,col=colors.black,yr=[0,6],back=colors.white,xtit='Time [s]',ytit='Stored energy, core Te, Ip [MJ,keV,MA]',/nodata,xr=tr,xs=1 & oplot,t,data/1e6,col=colors.red
    	if keyword_set(average)then begin
	    id=where(t ge average-0.2 and t le average+0.2)
	    print,'Diamagnetic energy [MJ]: ',mean(data[id]/1e6)
	endif
    	ppfread,shot=shot,dda='HRTX',dtype='TE0',data=data,t=t & oplot,t,data/1e3,col=colors.blue
    	if keyword_set(average)then begin
	    id=where(t ge average-0.2 and t le average+0.2)
	    print,'Core Te [keV]: ',mean(data[id]/1e3)
	endif
    	ppfread,shot=shot,dda='MAGN',dtype='IPLA',data=data,t=t & oplot,t,-data/1e6,col=colors.orange
    	if keyword_set(average)then begin
	    id=where(t ge average-0.2 and t le average+0.2)
	    print,'Ip [MA]: ',mean(data[id]/1e6)
	endif
; Neutrons
    	ppfread,shot=shot,dda='TIN',dtype='RNT',data=data,t=t & plot,t,data/1e15,col=colors.black,back=colors.white,yr=[0,5],xtit='Time [s]',ytit='Neutron rates [10!u15!n #/s]',/nodata,xr=tr,xs=1 & oplot,t,data/1e15,col=colors.red
    	if keyword_set(average)then begin
	    id=where(t ge average-0.2 and t le average+0.2)
	    print,'Neutron rates [10^15]: ',mean(data[id]/1e15)
	endif
; Core and edge density
    	ppfread,shot=shot,dda='KG1V',dtype='LID3',data=data,t=t & plot,t,data/1e19,col=colors.black,back=colors.white,yr=[0,30],xtit='Time [s]',ytit='Line integrated density [10!u19!n m!u-2!n]',/nodata,xr=tr,xs=1 & oplot,t,data/1e19,col=colors.red
    	ppfread,shot=shot,dda='KG1V',dtype='LID4',data=data,t=t & oplot,t,data/1e19,col=colors.blue
    	if keyword_set(average)then begin
	    id=where(t ge average-0.2 and t le average+0.2)
	    print,'Line integrated edge density [10^19]: ',mean(data[id]/1e19)
	endif
   
; Triangularity 
    	ppfread,shot=shot,dda='EFIT',dtype='TRIL',data=data,t=t & plot,t,data,col=colors.black,back=colors.white,xtit='Time [s]',ytit='Triangularity',/nodata,xr=tr,xs=1,yr=[0,0.5] & oplot,t,data,col=colors.blue
    	ppfread,shot=shot,dda='EFIT',dtype='TRIU',data=data1,t=t1 & oplot,t1,data1,col=colors.red
    	tri = (data + data1)/2.0
	oplot,t,tri,col=colors.black
    	if keyword_set(average)then begin
	    id=where(t ge average-0.2 and t le average+0.2)
	    print,'Triangularity: ',mean(tri[id])
	endif
; Zeff
    	ppfread,shot=shot,dda='KS3',dtype='ZEFV',data=data,t=t & plot,t,data,col=colors.black,back=colors.white,xtit='Time [s]',ytit='ZEFF',/nodata,xr=tr,xs=1,yr=[0,5] & oplot,t,data,col=colors.red

; H-factor    
    	ppfread,shot=shot,dda='SCAL',dtype='H98Y',data=data,t=t & if n_elements(t) gt 1 then begin plot,t,data,col=colors.black,back=colors.white,xtit='Time [s]',ytit='H98Y',/nodata,xr=tr,xs=1,yr=[0,1.5] & oplot,t,data,col=colors.red & endif
    	if keyword_set(average)then begin
	    id=where(t ge average-0.2 and t le average+0.2)
	    print,'H98: ',mean(data[id])
	endif

    endfor
stop
End    


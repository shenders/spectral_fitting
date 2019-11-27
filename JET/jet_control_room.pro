PRO jet_control_room,shot,tr=tr
    
    if ~keyword_set(tr)then tr=[40,60]

    setgraphics,xs=1000,ys=1000,psplot=psplot,file='jet_trace.ps',colors=colors,nrow=3,ncol=2
; Total input and radiated power
    ppfread,shot=shot,dda='NBI',dtype='NBLM',data=data,t=t & plot,t,data/1e6,col=colors.black,back=colors.white,xtit='Time [s]',ytit='Input power [MW]',yr=[0,30],/nodata,xr=tr,xs=1 & oplot,t,data/1e6,col=colors.blue
    ppfread,shot=shot,dda='ICRH',dtype='PTOT',data=data2,t=t2 & oplot,t2,data2/1e6,col=colors.red
    oplot,t,(data + interpol(data2,t2,t))/1e6,col=colors.black
; Radiated power   
    ppfread,shot=shot,dda='BOLO',dtype='TOPI',data=data,t=t & plot,t,data/1e6>0,col=colors.black,back=colors.white,xtit='Time [s]',ytit='Total radiated power [MW]',xr=tr,xs=1,yr=[0,mean(data/1e6)*5.0]
; Gas valves
    ppfread,shot=shot,dda='GASM',dtype='MAJR',data=data,t=t & plot,t,data/1e21,col=colors.black,back=colors.white,xtit='Time [s]',ytit='Gas valve, cZ [10!u21!n #/s, %]',/nodata,xr=tr,xs=1 & oplot,t,data/1e21,col=colors.red
    ppfread,shot=shot,dda='GASM',dtype='MN1R',data=data2,t=t2 & oplot,t2,data2/1e21,col=colors.blue
    oplot,t,(interpol(data2/7,t2,t) / (data+interpol(data2/7,t2,t)))*100,col=colors.black
; Diamagnetic energy
    ppfread,shot=shot,dda='EFIT',dtype='WDIA',data=data,t=t & plot,t,data/1e6,col=colors.black,yr=[0,6],back=colors.white,xtit='Time [s]',ytit='Stored energy, core Te, Ip [MJ,keV,MA]',/nodata,xr=tr,xs=1 & oplot,t,data/1e6,col=colors.red
    ppfread,shot=shot,dda='HRTX',dtype='TE0',data=data,t=t & oplot,t,data/1e3,col=colors.blue
    ppfread,shot=shot,dda='MAGN',dtype='IPLA',data=data,t=t & oplot,t,-data/1e6,col=colors.orange
; Neutrons
    ppfread,shot=shot,dda='TIN',dtype='RNT',data=data,t=t & plot,t,data/1e15,col=colors.black,back=colors.white,xtit='Time [s]',ytit='Neutron rates [10!u15!n #/s]',/nodata,xr=tr,xs=1 & oplot,t,data/1e15,col=colors.red
; Core and edge density
    ppfread,shot=shot,dda='KG1V',dtype='LID3',data=data,t=t & plot,t,data/1e19,col=colors.black,back=colors.white,xtit='Time [s]',ytit='Line integrated density [10!u19!n m!u-2!n]',/nodata,xr=tr,xs=1 & oplot,t,data/1e19,col=colors.red
    ppfread,shot=shot,dda='KG1V',dtype='LID4',data=data,t=t & oplot,t,data/1e19,col=colors.blue
   
    
stop
End    


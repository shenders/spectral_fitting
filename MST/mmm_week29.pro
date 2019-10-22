PRO basic_fit,te_val=te_val,func=func,dens=dens,setdens=setdens,tec3995=tec3995,exc3995=exc3995,rec3995=rec3995,use402=use402
	
	if keyword_set(setdens)then begin
		te=adas_vector(high=10,low=2,num=100)
		te_val=te
		dens = fltarr(100)+setdens
	endif else begin
		dens=adas_vector(high=1e15,low=1e13,num=100)
		if ~keyword_set(te_val)then te_val=3.5
		te = te_val+fltarr(100)
	end
; 399.5 line
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=exc3995,block=15
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=rec3995,block=65
; 404.2 line	
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=exc4042,block=21
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=rec4042,block=71
; Ionisation balance
    	run_adas405, uid='adas', year='96', elem='n', te=te, dens=dens, frac=frac
; Line ratio at fixed temperature	
	tec3995 = (frac.ion[*,1] * exc3995 + frac.ion[*,2] * rec3995) 
	tec4042 = (frac.ion[*,1] * exc4042 + frac.ion[*,2] * rec4042) 

	if keyword_set(use402)then begin
	; 402.6 line
		read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=exc4026,block=19
		read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=rec4026,block=69
		tec3995 = (frac.ion[*,1] * exc4026 + frac.ion[*,2] * rec4026) 
	endif
	func    = tec4042 / tec3995
	
END
PRO mmm_week29,psplot=psplot

x=fetch_data(36702,['ROV-09','ROV-10','ROV-12','ROV-14','ZON-01','ZON-03','ZON-05'],tr=[2.0,4.0],machine='AUG',diag='EVS',res=0.02,species=['O','C','He','W'],/quick,/load)

setgraphics,nrow=2,ncol=1,xs=600,ys=1000,psplot=psplot,file='HeI36702.ps',colors=colors

plot,x.time,x.hei[*,2,0]/1e17,title='ZON sightlines',back=colors.white,xtitle='Time [s]',ytitle='He I [10!u17!n ph/s/m!u2!n/sr]',col=colors.black 
oplot,x.time,x.hei[*,1,0]/1e17,col=colors.blue & oplot,x.time,x.hei[*,0,0]/1e17,col=colors.cyan


plot,x.time,x.hei[*,6,0]/1e17,title='ROV sightlines',col=colors.black,back=colors.white,yr=[0,3.5],xtitle='Time [s]',ytitle='He I [10!u17!n ph/s/m!u2!n/sr]' 
oplot,x.time,x.hei[*,5,0]/1e17,col=colors.blue & oplot,x.time,x.hei[*,4,0]/1e17,col=colors.cyan & oplot,x.time,x.hei[*,3,0]/1e17,col=colors.green

x=fetch_data(36703,['ROV-09','ROV-10','ROV-12','ROV-14','ZON-01','ZON-03','ZON-05'],tr=[2.0,4.0],machine='AUG',diag='EVS',res=0.02,species=['O','C','He','W'],/quick,/load)

setgraphics,nrow=2,ncol=2,xs=1200,ys=800,psplot=psplot,file='NII36703.ps',colors=colors,title='36703'

plot,x.time,x.nii[*,2,0]/1e17,title='ZON sightlines',back=colors.white,xtitle='Time [s]',ytitle='N II [10!u17!n ph/s/m!u2!n/sr]',col=colors.black 
oplot,x.time,x.nii[*,1,0]/1e17,col=colors.blue & oplot,x.time,x.nii[*,0,0]/1e17,col=colors.cyan


plot,x.time,x.nii[*,6,0]/1e17,title='ROV sightlines',col=colors.black,back=colors.white,yr=[0,40],xtitle='Time [s]',ytitle='N II [10!u17!n ph/s/m!u2!n/sr]' 
oplot,x.time,x.nii[*,5,0]/1e17,col=colors.blue & oplot,x.time,x.nii[*,4,0]/1e17,col=colors.cyan & oplot,x.time,x.nii[*,3,0]/1e17,col=colors.green


plot,x.time,x.nii[*,2,1]/x.nii[*,2,0],title='ZON sightlines',back=colors.white,yr=[0,0.7],xtitle='Time [s]',ytitle='N II 4041 / 3995 ratio',col=colors.black 
oplot,x.time,x.nii[*,1,1]/x.nii[*,1,0],col=colors.blue & oplot,x.time,x.nii[*,0,1]/x.nii[*,0,0],col=colors.cyan


plot,x.time,x.nii[*,6,1]/x.nii[*,6,0],title='ROV sightlines',col=colors.black,back=colors.white,yr=[0,0.7],xtitle='Time [s]',ytitle='N II 4041 / 3995 ratio' 
oplot,x.time,x.nii[*,5,1]/x.nii[*,5,0],col=colors.blue & oplot,x.time,x.nii[*,4,1]/x.nii[*,4,0],col=colors.cyan & oplot,x.time,x.nii[*,3,1]/x.nii[*,3,0],col=colors.green

x=fetch_data(36709,['ROV-09','ROV-10','ROV-12','ROV-14','ZON-01','ZON-03','ZON-05'],tr=[2.0,4.0],machine='AUG',diag='EVS',res=0.02,species=['O','C','He','W'],/quick,/load)
y=fetch_data(36709,['ROV-09','ROV-10','ROV-12','ROV-14','ZON-01','ZON-03','ZON-05'],tr=[2.0,4.0],machine='AUG',diag='EVS',res=0.02,species=['O','C','He','W'],/quick,/load,append='stark')

setgraphics,nrow=2,ncol=2,xs=1000,ys=800,psplot=psplot,file='36709.ps',colors=colors,title='36709',/landscape

plot,x.time,x.nii[*,6,0]/1e17,title='Nitrogen',col=colors.black,xr=[2.0,3.3],xs=1,back=colors.white,yr=[0,40],xtitle='Time [s]',ytitle='N II [10!u17!n ph/s/m!u2!n/sr]' 
oplot,x.time,x.nii[*,5,0]/1e17,col=colors.blue & oplot,x.time,x.nii[*,4,0]/1e17,col=colors.cyan & oplot,x.time,x.nii[*,3,0]/1e17,col=colors.green

plot,y.time,y.balmer[*,3,0]/1e16,yr=[0,20],title='Deuterium',xr=[2.0,3.3],xs=1,back=colors.white,xtitle='Time [s]',ytitle='D I [10!u16!n ph/s/m!u2!n/sr]',col=colors.black,/nodata 
oplot,y.time,y.balmer[*,3,0]/1e16,col=colors.green

te_lower = 3.5
basic_fit,te_val=te_lower,func=func_te_lower,dens=dens_te_lower

plot,x.time,interpol(dens_te_lower,func_te_lower,x.nii[*,6,1]/x.nii[*,6,0])/1e14,xr=[2.0,3.3],xs=1,title='N II ratio analysis',yr=[0,20],col=colors.black,back=colors.white,xtitle='Time [s]',ytitle='n!le!n [10!u20!n m!u-3!n]' 
oplot,x.time,interpol(dens_te_lower,func_te_lower,x.nii[*,5,1]/x.nii[*,5,0])/1e14,col=colors.blue 
oplot,x.time,interpol(dens_te_lower,func_te_lower,x.nii[*,4,1]/x.nii[*,4,0])/1e14,col=colors.cyan 
oplot,x.time,interpol(dens_te_lower,func_te_lower,x.nii[*,3,1]/x.nii[*,3,0])/1e14,col=colors.green

sm=10

plot,y.time,smooth(y.dens_balmer[*,3,0]/1e20,sm,/edge_truncate),xr=[2.0,3.3],xs=1,title='Stark broadening',back=colors.white,xtitle='Time [s]',ytitle='n!le!n [10!u20!n m!u-3!n]',col=colors.black,/nodata 
oplot,y.time,smooth(y.dens_balmer[*,3,0]/1e20,sm,/edge_truncate),col=colors.green

setgraphics,nrow=2,ncol=2,xs=1000,ys=800,psplot=psplot,file='test.ps',colors=colors,title='36709',/landscape

END

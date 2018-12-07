PRO basic_fit,te_val=te_val,func=func,dens=dens,tec3995=tec3995,exc3995=exc3995,rec3995=rec3995,use402=use402
	dens=adas_vector(high=1e15,low=1e13,num=100)
	if ~keyword_set(te_val)then te_val=3.5
	te = te_val+fltarr(100)
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
Pro analysis,shot,trace,sv=sv,tr=tr,sightline=sightline,dl=dl,psplot=psplot

shotstr  = string(shot,format='(i5)') 
trace    = ['ROV-11-interelm-25ms',$
	    'ROV-12-interelm-25ms',$
	    'ROV-14-interelm-25ms']

restore,'save/'+shotstr+'/'+trace[0]+'.idl'

time    = output.time
nii3995 = fltarr(n_elements(time),n_elements(trace))
nii4041 = fltarr(n_elements(time),n_elements(trace))
density = fltarr(n_elements(time),n_elements(trace))

for i=0,n_elements(trace)-1 do begin
	restore,'save/'+shotstr+'/'+trace[i]+'.idl'
	nii3995[*,i] = output.nii[*,0,0]
	nii4041[*,i] = output.nii[*,0,1]
	density[*,i] = output.dens_balmer[*,0]
endfor	    

if keyword_set(psplot)then begin
	makeps,file='figures/analysis.ps',xs=10,ys=8,/landscape
endif else begin
	window,0,xsize=1200,ysize=800
end
adas_colors,colors=colors
colpick = [colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.orange,colors.red]
!p.multi=[0,2,2]
!p.thick=2.0
if ~keyword_set(sv)then sv = 6

plot,time,nii3995[*,0],xr=tr,back=colors.white,xtitle='Time [s]',ytitle='ph/s/nm/sr',col=colors.black,/nodata,title='N II emission @399.5 nm)',yr=[0,max(nii3995)]
for i=0,n_elements(trace)-1 do oplot,time,smooth(nii3995[*,i],sv),col=colpick[i]

plot,time,smooth(density[*,0],sv),xr=tr,/nodata,xtitle='Time [s]',ytitle='m-3',back=colors.white,col=colors.black,title='Density (Stark broadening)'
for i=0,n_elements(trace)-1 do oplot,time,smooth(density[*,i],sv),col=colpick[i]

; Temperature upper and lower limits
te_lower = 3.5
te_upper = 4.0

; Density function based on temperature limits
basic_fit,te_val=te_lower,func=func_te_lower,dens=dens_te_lower
basic_fit,te_val=te_upper,func=func_te_upper,dens=dens_te_upper

plot,time,fltarr(n_elements(time)),xr=tr,/nodata,xtitle='Time [s]',ytitle='m-3',yr=[0,3e20],title='Density (N II line ratio)',col=colors.black,back=colors.white
densfitu = fltarr(n_elements(time),n_elements(trace))
densfitl = fltarr(n_elements(time),n_elements(trace))
for i=0,n_elements(trace)-1 do begin
	densfitu[*,i] = interpol(dens_te_upper,func_te_upper,smooth(nii4041[*,i],sv)/smooth(nii3995[*,i],sv))
	densfitl[*,i] = interpol(dens_te_lower,func_te_lower,smooth(nii4041[*,i],sv)/smooth(nii3995[*,i],sv))
	oband,time,densfitu[*,i]*1e6,densfitl[*,i]*1e6,col=colpick[i],/norm
endfor

plot,time,fltarr(n_elements(time)),xr=tr,/nodata,xtitle='Time [s]',ytitle='%',back=colors.white,col=colors.black,yr=[0,50],title='Nitrogen concentration'
if ~keyword_set(sightline)then sightline = [1,2,4]
if ~keyword_set(dl)then DL = [0.055,0.075,0.1]
for i=0,n_elements(sightline)-1 do begin
	; Guess transmission
	transmission = 1.0/0.64
	; Guess Delta L
	DeltaL = DL[i]
		
	; Get nitrogen concentrations from spectroscopy

	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_lower+fltarr(n_elements(time)),$
		dens=densfitl[*,sightline[i]],data=exc3995,block=15
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_lower+fltarr(n_elements(time)),$
		dens=densfitl[*,sightline[i]],data=rec3995,block=65
			      
	run_adas405,elem='n',year='96',uid='adas',te=te_lower+fltarr(n_elements(time)),$
		dens=densfitl[*,sightline[i]],frac=frac

	pec3995 = frac.ion[*,1] * exc3995 + frac.ion[*,2] * rec3995

	cn_low = 3.14 * 4.0 * smooth(nii3995[*,sightline[i]],sv) * transmission / (densfitl[*,sightline[i]] * pec3995) / DeltaL / (densfitl[*,sightline[i]]*1e6)

	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_upper+fltarr(n_elements(time)),$
		dens=densfitu[*,sightline[i]],data=exc3995,block=15
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te_upper+fltarr(n_elements(time)),$
		dens=densfitu[*,sightline[i]],data=rec3995,block=65
			      
	run_adas405,elem='n',year='96',uid='adas',te=te_upper+fltarr(n_elements(time)),$
		dens=densfitu[*,sightline[i]],frac=frac

	pec3995 = frac.ion[*,1] * exc3995 + frac.ion[*,2] * rec3995

	cn_up = 3.14 * 4.0 * smooth(nii3995[*,sightline[i]],sv) * transmission / (densfitu[*,sightline[i]] * pec3995) / DeltaL / (densfitu[*,sightline[i]]*1e6)

	oband,time,cn_up * 100,cn_low * 100,/norm,col=colpick[sightline[i]]
endfor
read_data,'Npuff',shot,output_n2,time_n2
read_data,'Dpuff',shot,output_d2,time_d2
oplot,time_n2,(output_n2/7.0)/(output_n2/7.0 + output_d2)*100.0,col=colors.black,linest=5
end
	
	    
	    
	    

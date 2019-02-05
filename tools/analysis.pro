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
Pro analysis,shot,los=los,append=append,sv=sv,sm=sm,xr=xr,$
		  interelm=interelm,sightline=sightline,dl=dl,$
		  psplot=psplot,trans=trans,upperte=upperte,lowerte=lowerte,$
		  channel=channel,use_evl=use_evl,setdens=setdens,ratio402=ratio402,$
		  diag=diag,sig3995=sig3995,sig4041=sig4041,exp=exp,no402=no402,png=png,$
		  time=time,cn_up=cn_up,cn_low=cn_low,debug=debug,tdiv=tdiv


shotstr  = string(shot,format='(i5)') 
if ~keyword_set(append)then append='data'
if ~keyword_set(los)then los='ROV014' 
sig=los
if ~keyword_set(diag)then diag='EVL' 
if ~keyword_set(sig3995)then sig3995='N_1_3995' 
if ~keyword_set(sig4041)then sig4041='N_1_4041' 
if ~keyword_set(sv)then sv = 15
if ~keyword_set(sm)then sv402 = sv else sv402=sm
if ~keyword_set(trace)then trace='save/'+shotstr+'/'+sig+'-'+append+'.idl'
if ~keyword_Set(channel)then channel=22

if ~keyword_set(use_evl)then begin
	restore,trace[0]
	time    = output.time

endif else begin
	read_signal_mrm,0L,shot,diag,sig3995,time,nii_3995,2,exp=exp
	read_signal_mrm,0L,shot,diag,sig4041,time,nii_4041,2,exp=exp
end
if ~keyword_set(xr)then xr = [min(time),max(time)]
if keyword_set(interelm)then begin
	telm    = find_elm(shot,time)
	if ~keyword_set(elmcond)then elmcond=4.5
	idelm     = where(telm ge elmcond)
	if idelm[0] eq -1 then idelm = findgen(n_elements(time))
endif else idelm = findgen(n_elements(time))

nii3995 = fltarr(n_elements(time),n_elements(trace))
nii4041 = fltarr(n_elements(time),n_elements(trace))
nii4026 = fltarr(n_elements(time),n_elements(trace))
nii4091 = fltarr(n_elements(time),n_elements(trace))
density = fltarr(n_elements(time),n_elements(trace))

for i=0,n_elements(trace)-1 do begin
	if ~keyword_set(use_evl) then begin
		restore,trace[i]
		nii3995[*,i] = output.nii[*,0,0]
		nii4041[*,i] = output.nii[*,0,1]
		nii4026[*,i] = output.nii[*,0,2]
		nii4091[*,i] = output.nii[*,0,4]
		density[*,i] = output.dens_balmer[*,0]
	endif else begin
		nii3995[*,i] = nii_3995[*,channel]
		nii4041[*,i] = nii_4041[*,channel]
		nii4026[*,i] = nii_4041[*,channel]
		nii4091[*,i] = nii_4041[*,channel]
		density[*,i] = fltarr(n_elements(time))+1e20
	endelse		
endfor	    

if keyword_set(debug)then begin
	if keyword_set(psplot)then begin
		xs= 10 & ys=8
		filename='figures/'+shotstr+'_NII_analysis.ps'
		landscape=1
	endif else begin
		xs=1200 & ys=800
	end
	setgraphics,xs=xs,ys=ys,ncol=2,nrow=2,psplot=psplot,filename=filename,portrait=portrait,landscape=landscape,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list	

	colpick = [colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.orange,colors.red]
	plot,time,nii3995[*,0]/1e19,xr=xr,xs=1,back=colors.white,xtitle='Time [s]',ytitle='Radiance [10!u19!n ph/s/m!u2!n/sr]',col=colors.black,/nodata,title='N II @ 399.5 nm',yr=[0,max(nii3995)/1e19]
	user_psym,1
	for i=0,n_elements(trace)-1 do oplot,time,nii3995[*,i]/1e19,col=colors.black,psym=8
	if keyword_set(interelm)then begin
		user_psym,1,/fill
		for i=0,n_elements(trace)-1 do oplot,time[idelm],nii3995[idelm,i]/1e19,col=colors.red,psym=8
	endif
endif
nii3995=interpol(nii3995[idelm,*],time[idelm],time)
nii4041=interpol(nii4041[idelm,*],time[idelm],time)
nii4026=interpol(nii4026[idelm,*],time[idelm],time)
density=interpol(density[idelm,*],time[idelm],time)

idchk = where(nii3995 gt 1e17)
nii3995=nii3995[idchk]
nii4041=nii4041[idchk]
nii4026=nii4026[idchk]
density=density[idchk]
time = time[idchk]

if keyword_set(ratio402)then rat402 = fltarr(n_elements(time))+ratio402
if keyword_set(debug)then begin

	oplot,time,smooth(nii3995,sv)/1e19,col=colpick[0]

	plot,time,smooth(nii4041,sv)/smooth(nii3995,sv),xr=xr,xs=1,/nodata,xtitle='Time [s]',ytitle='[-]',$
	  back=colors.white,col=colors.black,title='N II line ratios',yr=[0,0.4]
	
	oplot,time,smooth(nii4041,sv)/smooth(nii3995,sv),col=colpick[0]
	xyouts,[min(xr)+0.1],[0.28],['404.1/399.5'],col=colors.black
	if ~keyword_set(no402)then begin
		oplot,time,(smooth(nii4041,sv402)/smooth(nii4026,sv402))/100,col=colpick[2]
		xyouts,[min(xr)+0.1],[0.11],['404.1/402.6 /100'],col=colors.black
		if keyword_set(ratio402)then begin
			oplot,time,rat402/100,col=colors.red
		endif
	endif
endif
; Temperature upper and lower limits
if ~keyword_set(lowerte)then te_lower = 3.3 else te_lower = lowerte
if ~keyword_set(upperte)then te_upper = te_lower+0.5 else te_upper = upperte

; Density function based on temperature limits
basic_fit,te_val=te_lower,func=func_te_lower,dens=dens_te_lower
basic_fit,te_val=te_upper,func=func_te_upper,dens=dens_te_upper
basic_fit,te_val=te_lower,func=func_te_lower_402,dens=dens_te_lower_402,/use402
basic_fit,te_val=te_upper,func=func_te_upper_402,dens=dens_te_upper_402,/use402
if keyword_set(debug)then begin
	plot,time,fltarr(n_elements(time)),xr=xr,xs=1,/nodata,xtitle='Time [s]',ytitle='n!le,N II!n [10!u20!n m!u-3!n]',yr=[0,4],title='Electron density',col=colors.black,back=colors.white
endif
densfitu = fltarr(n_elements(time))
densfitl = fltarr(n_elements(time))
densfitu_402 = fltarr(n_elements(time))
densfitl_402 = fltarr(n_elements(time))

densfitu = interpol(dens_te_upper,func_te_upper,smooth(nii4041,sv)/smooth(nii3995,sv))
densfitl = interpol(dens_te_lower,func_te_lower,smooth(nii4041,sv)/smooth(nii3995,sv))

if keyword_set(debug)then begin
	oband,time,densfitu*1e6/1e20,densfitl*1e6/1e20,col=colpick[0],/norm
	xyouts,[mean([xr[0],mean(xr)])],[3],[string(te_lower,te_upper,format='(D3.1," eV < T!le!n < ",D3.1," eV")')],col=colors.black
endif
if ~keyword_set(no402)then begin
	if ~keyword_set(ratio402)then rat402 = (smooth(nii4041,sv402)/smooth(nii4026,sv402))<max(func_te_upper_402)
	densfitu_402 = interpol(dens_te_upper_402,func_te_upper_402,rat402)
	densfitl_402 = interpol(dens_te_lower_402,func_te_lower_402,rat402)
	if keyword_set(debug)then begin
		oplot,time,densfitu_402*1e6/1e20,col=colpick[2]
		oplot,time,densfitl_402*1e6/1e20,col=colpick[2]
	endif		
endif

if keyword_set(setdens)then begin
	densfitu = densfitu_402>5e19
	densfitl = densfitl_402>5e19
	te_lower = fltarr(n_elements(time))
	te_upper = fltarr(n_elements(time))
	ratio = smooth(nii4041,sv)/smooth(nii3995,sv)
	for j=0,n_elements(densfitu[*,0])-1 do begin
		basic_fit,te_val=tel,func=func_te_lower,setdens=densfitl[j]
		basic_fit,te_val=teu,func=func_te_upper,setdens=densfitu[j]
		te_lower[j] = interpol(tel,func_te_lower,ratio[j])
		te_upper[j] = interpol(teu,func_te_upper,ratio[j])
	endfor
endif	

if keyword_set(debug)then begin
	plot,time,fltarr(n_elements(time)),xr=xr,xs=1,/nodata,xtitle='Time [s]',ytitle='c!lN!n [%]',back=colors.white,col=colors.black,yr=[0,30],title='Nitrogen concentration'
endif

read_data,'Tdiv',shot,output_tdiv,time_tdiv
id =where(time_tdiv ge 1 and time_tdiv le 8)
tdiv = interpol(output_tdiv[id,0],time_tdiv[id],time)
losnum = fix(strmid(los,4,2))
if losnum ge 13 then rov014 = 1 else rov014 = 0
deltaL,tdiv,dl,rov014=rov014

; Guess transmission
if ~keyword_set(trans)then transmission=1.0 else transmission=trans
;transmission = 1.0;/0.64
; Guess Delta L
DeltaL = dl
		
; Get nitrogen concentrations from spectroscopy
if ~keyword_set(setdens)then begin
	telow = te_lower + fltarr(n_elements(time))
	teupp = te_upper + fltarr(n_elements(time))
endif else begin
	telow = te_lower
	teupp = te_upper
end

read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=telow,$
		dens=densfitl,data=exc3995,block=15
read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=telow,$
		dens=densfitl,data=rec3995,block=65
			      
run_adas405,elem='n',year='96',uid='adas',te=telow,$
		dens=densfitl,frac=frac

pec3995 = frac.ion[*,1] * exc3995 + frac.ion[*,2] * rec3995

cn_low = 3.14 * 4.0 * smooth(nii3995,sv) * transmission / (densfitl * pec3995) / DeltaL / (densfitl*1e6)

read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=teupp,$
		dens=densfitu,data=exc3995,block=15
read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=teupp,$
		dens=densfitu,data=rec3995,block=65
			      
run_adas405,elem='n',year='96',uid='adas',te=teupp,$
		dens=densfitu,frac=frac

pec3995 = frac.ion[*,1] * exc3995 + frac.ion[*,2] * rec3995

cn_up = 3.14 * 4.0 * smooth(nii3995,sv) * transmission / (densfitu * pec3995) / DeltaL / (densfitu*1e6)
if keyword_set(debug)then begin
	oband,time,(cn_up * 100)>0.1,(cn_low * 100)>0.1,/norm,col=colpick[0]
endif
if keyword_set(debug)then begin
	id = where(time_tdiv ge xr[0] and time_tdiv le xr[1])
	oplot,time_tdiv[id],output_tdiv[id],col=colors.red
	xyouts,[mean(xr)],[25],['Tdiv'],col=colors.red
	xyouts,[mean(xr)],[22],['cN [N II]'],col=colpick[0]
	xyouts,[mean(xr)],[19],['cN [Flux]'],col=colors.black
	read_signal_mrm,0L,shot,'UVS','N_tot',time_n2,output_n2,1
	read_signal_mrm,0L,shot,'UVS','D_tot',time_d2,output_d2,1
	oplot,time_n2,(output_n2/7.0)/(output_n2/7.0 + output_d2)*100.0,col=colors.black,linest=5
	if keyword_set(psplot)then setgraphics,/close,psplot=psplot
	if keyword_set(png)then write_png,'figures/'+string(shot,format='(i5)')+sig+'.png',tvrd(/true)
	stop
endif
if keyword_set(xr)then begin
	id = where(time ge xr[0] and time le xr[1])
	time = time[id]
	cn_up= cn_up[id]
	cn_low=cn_low[id]
	tdiv=tdiv[id]
endif

end
	
	    
	    
	    

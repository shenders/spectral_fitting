PRO aug_div_pressure,shot,trange,press,error,scl=scl,debug=debug,nosave=nosave,restore=restore

; Set constants
	if ~keyword_set(scl)then scl=2.5
	kb           = 1.38064E-23
	flux_to_dens = 4.0 / 1240.0
	dens_to_pres = 293.15 * kb
	torr_to_Pa   = 133.32
	
; real-time estimate

	read_signal_mrm,0L,shot,'DDS','nDivIst',time_dds,data_dds,1
	
	Pdiv_dds = data_dds * dens_to_pres

; calibrated measurements

	read_signal_mrm,0L,shot,'IOC','F03',time_ioc,data_ioc_f03,1
	id_f03 = where(finite(data_ioc_f03))	
	Pdiv_ioc_f03 = (data_ioc_f03 * flux_to_dens) * dens_to_pres

	read_signal_mrm,0L,shot,'IOC','F01',time_ioc,data_ioc_f01,1
	id_f01 = where(finite(data_ioc_f01))	
	Pdiv_ioc_f01 = (data_ioc_f01 * flux_to_dens) * dens_to_pres

	read_signal_mrm,0L,shot,'IOC','F09',time_ioc,data_ioc_f09,1	
	id_f09 = where(finite(data_ioc_f09))	
	Pdiv_ioc_f09 = (data_ioc_f09 * flux_to_dens) * dens_to_pres

	read_signal_mrm,0L,shot,'IOC','F04',time_ioc,data_ioc_f04,1	
	id_f04 = where(finite(data_ioc_f04))	
	Pdiv_ioc_f04 = (data_ioc_f04 * flux_to_dens) * dens_to_pres

; reliable measurements

	read_signal_mrm,0L,shot,'MSP','B25_08Fu',time_msp,data_msp,1	
	Pdiv_msp = (data_msp );* torr_to_Pa)
	Pdiv_msp = 100 * Pdiv_msp /(0.0337 * alog(Pdiv_msp*100)+0.7304)

; plot 
	iter=0

	if keyword_set(debug)then setgraphics,colors=colors,xs=600,ys=400
	rep:
	tmp = moment(Pdiv_msp[where(time_msp ge trange[0] and time_msp le trange[1])])
	press = tmp[0]
	error = sqrt(tmp[1])
	if iter eq 1 then begin
		press = y
		error = y*0.2
	endif
	if keyword_set(debug)then begin
		iter=0
		plot,time_dds,Pdiv_dds*scl,col=colors.black,title=string(shot,format='(i5)'),/nodata,back=colors.white,xr=[1.0,7.0],yr=[0,max(Pdiv_ioc_f09[id_f09])>max(Pdiv_ioc_f03[id_f03])>max(Pdiv_ioc_f04[id_f04])>max(Pdiv_ioc_f01[id_f01])]
		oplot,time_ioc[id_f09],Pdiv_ioc_f09[id_f09],col=colors.red
		oplot,time_ioc[id_f03],Pdiv_ioc_f03[id_f03],col=colors.cyan
		oplot,time_ioc[id_f04],Pdiv_ioc_f04[id_f04],col=colors.green
		oplot,time_ioc[id_f01],Pdiv_ioc_f01[id_f01],col=colors.orange
		oplot,time_msp,Pdiv_msp,col=colors.blue
		oplot,time_dds,Pdiv_dds*scl,col=colors.black
	
		oplot,[trange[0],trange[0]],[0,5000],col=colors.black,linest=5
		oplot,[trange[1],trange[1]],[0,5000],col=colors.black,linest=5

		user_psym,5,/fill & oplot,[mean(trange),-10],[press,-10],col=colors.red,psym=8
		errors,[mean(trange),-10],[press,-10],ystd=[error,-10],col=colors.red
		iter = 1
		if keyword_set(restore)then restore,'save/'+string(shot,format='(I5)')+'/div_pressure.sav'
		if keyword_set(restore)then errors,[mean(trange),-10],[press,-10],ystd=[error,-10],col=colors.green
		cursor,x,y,/up
		if !mouse.button eq 4 then goto,rep

	endif
	cmd='mkdir -p save/'+string(shot,format='(I5)')
	spawn,cmd
	cmd = 'mv save/'+string(shot,format='(I5)')+'/div_pressure.sav save/'+string(shot,format='(I5)')+'/div_pressure2.sav' 
	if ~keyword_set(nosave)then spawn,cmd
	if ~keyword_set(nosave)then save,file='save/'+string(shot,format='(I5)')+'/div_pressure.sav',press,error
END

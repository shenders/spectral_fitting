Pro graphs,shot,psplot=psplot

	
	read_signal_mrm,0L,shot,'EVL','N_1_3995',time,nii_3995,2,exp=exp
	read_signal_mrm,0L,shot,'EVL','N_1_4041',time,nii_4041,2,exp=exp
	if keyword_set(psplot)then begin
		xs = 10 & ys = 8
	endif else begin
		xs = 1000 & ys = 800
	end
	setgraphics,colors=colors,nrow=2,ncol=2,xs=xs,ys=ys,psplot=psplot,/landscape,filename='figures/35819NII.ps'
	colpick = [colors.black,colors.blue,colors.red]
	sm = 30
	chn = [9 , 17, 19]
	delta = [0.05,0.20,0.10]
	lowerte = [3.0, 3.4, 3.8]
	user_psym,1
	for i = 0,2 do begin
		if i eq 0 then plot,time,smooth(nii_3995[*,chn[i]],sm)/1e18,ytitle='[10!u18!n ph/s/m!u2!n/sr]',$
		               xr=[4.0,7.0],yr=[0,4],back=colors.white,col=colors.black,/nodata,xtitle='Time [s]'
		oplot,time,nii_3995[*,chn[i]]/1e18,col=colpick[i],psym=8
		oplot,time,smooth(nii_3995[*,chn[i]],sm)/1e18,col=colpick[i]
	endfor
	xyouts,[4.05],[3.5],['(a) N II @ 399.5 nm'],col=colors.black

	for i = 0,2 do begin
		if i eq 0 then plot,time,smooth(nii_4041[*,chn[i]],sm)/smooth(nii_3995[*,chn[i]],sm),ytitle='[-]',$
		               xr=[4.0,7.0],yr=[0,0.4],back=colors.white,col=colors.black,/nodata,xtitle='Time [s]'
		oplot,time,smooth(nii_4041[*,chn[i]],sm)/smooth(nii_3995[*,chn[i]],sm),col=colpick[i]
	endfor
	xyouts,[4.05],[0.35],['(b) N II 404.1 / 399.5 nm'],col=colors.black
	
	if keyword_set(psplot)then setgraphics,psplot=psplot,/close
	
	;setgraphics,colors=colors,nrow=1,ncol=2,xs=xs,ys=ys,psplot=psplot,/portrait,filename='figures/35819cN.ps'
	for i = 0,2 do begin
		analysis,shot,los=los,xr=[4.0,7.0],$
		  channel=chn[i],/use_evl,/no402,time=time,lowerte=lowerte[i],$
		  cn_up=cn_up,cn_low=cn_low,cn_mean=cn_mean,cn_err=cn_err,tdiv=tdiv,del=delta[i],$
		  tflux=tflux,fluxratio=fluxratio,/onlydens,udensfit=densfitu,ldensfit=densfitl
		if i eq 0 then plot,time,densfitu * 1e6 / 1e20,ytitle='[10!u20!n m!u-3!n]',$
		               xr=[4.0,7.0],yr=[0,4],back=colors.white,col=colors.black,/nodata,xtitle='Time [s]'
		oband,time,densfitu * 1e6 / 1e20,densfitl * 1e6 / 1e20,col=colpick[i],/norm
	endfor 
	xyouts,[4.05],[3.5],['(c) n!le,N II!u'],col=colors.black
	for i = 0,2 do begin
		analysis,shot,los=los,xr=[4.0,7.0],$
		  channel=chn[i],/use_evl,/no402,time=time,$
		  cn_up=cn_up,cn_low=cn_low,cn_mean=cn_mean,cn_err=cn_err,tdiv=tdiv,del=delta[i],$
		  tflux=tflux,fluxratio=fluxratio,lowerte=lowerte[i]
		if i eq 0 then plot,time,cn_mean * 100,ytitle='[%]',$
		               xr=[4.0,7.0],yr=[0,6],back=colors.white,col=colors.black,/nodata,xtitle='Time [s]'
		oband,time,cn_mean * 100,cn_err * 100,col=colpick[i]
	endfor 
	xyouts,[4.05],[(3.5/4.0)*6],['(d) c!lN!u'],col=colors.black
	oplot,tflux,fluxratio,col=colors.black,linest=5

	if keyword_set(psplot)then setgraphics,psplot=psplot,/close
	
	
stop
End

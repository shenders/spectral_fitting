Pro get_cn,shot,los,$
                debug=debug,$
		interelm=interelm,$
		use_evl=use_evl,$
		channel=channel,$
		diag=diag,$
		full=full,$
		xr=xr,$
		rawelm=rawelm,$
		elmcond=elmcond,$
		transmission=transmission,$
		sm=sm,$
		lowerte=lowerte,$
		upperte=upperte,$
		use_rov8=use_rov8,$
		dynamic=dynamic,$
		preset=preset,psplot=psplot
	!quiet=1

	if ~keyword_set(sm)then sm = 20
	if ~keyword_set(upperte)then upperte = 3.8
	if ~keyword_set(lowerte)then lowerte = 3.3

	; Get N II 3995, 4041, and 4042 lines

	data  =  get_nii(shot,$
		 los=los,$
		 xr=xr,$
		 interelm=interelm,$
 		 channel=channel,$
		 use_evl=use_evl,$
		 append=append,$
		 diag=diag,$
		 sig3995=sig3995,$
		 rawelm=rawelm,$
		 elmcond=elmcond,$
		 no404=no404,$
		 use_rov8=use_rov8,$
		 dynamic=dynamic,$
		 preset=preset)

	if ~keyword_set(transmission)then transmission=1.0

	; get cn, Te and dens from line ratios

	if keyword_set(full)then begin
		plasma = optimize_ratios(data,shot,los,transmission,sm=sm,debug=debug)
	endif else begin
		plasma = estimate_ratios(data,shot,los,transmission,sm=sm,debug=debug,lowerte=lowerte,upperte=upperte)
	end
	
	x1=-1 & y1=-1 & x2=-1 & y2=-1
	read_signal_mrm,0L,shot,'UVS','D_tot',x1,y1,2,exp=exp
	read_signal_mrm,0L,shot,'UVS','N_tot',x2,y2,2,exp=exp
	cn_flux = (y2/7)/(y2/7+y1)
	id = where(x1 ge min(data.time) and x1 le max(data.time))
	cn_flux = cn_flux[id]
	cn_time = x1[id]
	cn_tdiv = interpol(data.tdiv,data.time,cn_time)
	cn_idsort = sort(cn_tdiv)
	!mouse.button=0
	iter = 0
	x1   = 3.0
	x2   = 5.0
	nrow = 3
	ncol = 1
	if keyword_set(debug)then setgraphics,nrow=nrow,ncol=ncol,colors=colors,xs=600,ys=1000,psplot=psplot,file='figures/cn_analysis.ps'
	while !mouse.button ne 4 do begin
		ii   = 0
		exp_ratio1 = smooth(data.nii4041,sm,/edge_truncate) / smooth(data.nii3995,sm,/edge_truncate)
		exp_ratio2 = smooth(data.nii4041,sm,/edge_truncate) / smooth(data.nii4026,sm,/edge_truncate)
		rawtdiv    = interpol(data.tdiv,data.time,data.rawtime)
		idsort     = sort(data.tdiv)
		idsortraw  = sort(rawtdiv)
		xr = [min(data.time)*1.05,max(data.time)*0.95]
		xr = [0,15]
		if keyword_set(debug)then begin
			plot,data.tdiv[idsort],exp_ratio1[idsort],col=colors.black,back=colors.white,yr=[0.,0.3],/nodata,xr=xr,xs=1,$
			ytitle='N II line ratios [-]',position=graphpos(0,ii,nrow,ncol,xspc=xspc,yspc=yspc),xtickname=replicate(' ',30)
			oband,data.tdiv[idsort],plasma.ratio1_lower[idsort],plasma.ratio1_upper[idsort],/norm,col=colors.black
			oband,data.tdiv[idsort],plasma.ratio2_lower[idsort]/70,plasma.ratio2_upper[idsort]/70,/norm,col=colors.blue
			oplot,data.tdiv[idsort],exp_ratio1[idsort],col=colors.red
			oplot,data.tdiv[idsort],exp_ratio2[idsort]/70,col=colors.green
			legend,'404.1/399.5',colors.black,yshift=-0.07
			legend,'404.1/402.6/70',colors.black,yshift=-0.5
			ii = ii+1

			plot,data.tdiv[idsort],plasma.dens_upper[idsort]*1e6/1e20,col=colors.black,back=colors.white,$
			ytitle='[10!u20!n m-3; eV]',/nodata,xr=xr,xs=1,yr=[0,5],position=graphpos(0,ii,nrow,ncol,xspc=xspc,yspc=yspc),xtickname=replicate(' ',30)
			oband,data.tdiv[idsort],plasma.dens_lower[idsort]*1e6/1e20,plasma.dens_upper[idsort]*1e6/1e20,/norm,col=colors.black
			oband,data.tdiv[idsort],plasma.te_lower[idsort],plasma.te_upper[idsort],/norm,col=colors.red
			legend,'n!le,N II!n',colors.black,yshift=-0.6
			legend,'T!le,N II!n',colors.black
			ii = ii+1

			plot,data.tdiv[idsort],plasma.cn_upper[idsort] * 100.0,col=colors.black,back=colors.white,$
			xtitle='Tdiv [eV]',xr=xr,xs=1,position=graphpos(0,ii,nrow,ncol,xspc=xspc,yspc=yspc),$
			ytitle='c!lN!n [%]',/nodata,yr=[0,20]
			oband,data.tdiv[idsort],plasma.cn_lower[idsort]*100,plasma.cn_upper[idsort]*100,/norm,col=colors.black
			oplot,cn_tdiv[cn_idsort],cn_flux[cn_idsort]*100,col=colors.blue
			if iter ne 0 then cursor,x1,y1,/up
			if !mouse.button eq 2 then goto,quit
			oplot,[x1,x1],[0,40],col=colors.red
			if iter ne 0 then cursor,x2,y2,/up
			oplot,[x2,x2],[0,40],col=colors.red
		endif
		id      = where(data.tdiv ge x1<x2 and data.tdiv le x1>x2)
		cn_mean = (mean(plasma.cn_upper[id]*100) + mean(plasma.cn_lower[id]*100))/2.0
		cn_err  = abs(cn_mean-mean(plasma.cn_upper[id]*100))
		
		cn_dens = (mean(plasma.dens_upper[id]) + mean(plasma.dens_lower[id]))/2.0
		cn_dens_err  = abs(cn_dens-mean(plasma.dens_upper[id]))

		cn_te = (mean(plasma.te_upper[id]) + mean(plasma.te_lower[id]))/2.0
		cn_te_err  = abs(cn_te-mean(plasma.te_upper[id]))

		if keyword_set(debug)then begin
			user_psym,5,/fill
			oplot,[mean([x1,x2]),1000],[cn_mean,1000],col=colors.red,psym=8
			err_plot,[mean([x1,x2]),1000],[cn_mean,1000],[cn_err,0.1],col=colors.red		
			legend,'c!lN!n',colors.black,yshift=-0.2
			legend,'T!ldiv!n',colors.black
			if ~keyword_set(psplot) then cursor,x,y,/up else !mouse.button = 4
			if !mouse.button eq 2 then goto,quit
			iter = iter + 1
		endif else !mouse.button=4
	end
	if keyword_set(psplot)then setgraphics,/close
	if ~keyword_set(nosave)then begin
		conc = cn_mean
		err  = cn_err
		te_nii = cn_te
		te_nii_err = cn_te_err
		dens_nii = cn_dens
		dens_nii_err = cn_dens_err
		cmd = 'mkdir -p save/cn/'+string(shot,format='(i5)')
		spawn,cmd
		file = 'save/cn/'+string(shot,format='(i5)')+'/cn_database_2.sav'
		save,file=file,conc,err,te_nii,te_nii_err,dens_nii,dens_nii_err
		print,string(conc,err,format='("cN: ",f5.2," +/- ",f4.2)')		
		print,'Save cN to file: ',file
	endif
	quit:
End

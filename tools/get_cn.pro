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
		nocn=nocn,$
		preset=preset,psplot=psplot,$
		t0=t0,$
		load=load,$
		fluxonly=fluxonly,$
		switchgas=switchgas,$
		overplot=overplot,$
		colpick=colpick
	!quiet=1
    	
	if ~keyword_set(t0)then t0 = 40.0
	if ~keyword_set(sm)then sm = 10
	if ~keyword_set(upperte)then upperte = 3.8
	if ~keyword_set(lowerte)then lowerte = 3.3

	adas_colors,colors=colors
	if keyword_set(colpick)then begin
	    if strlowcase(colpick) eq 'black' then coluse=colors.black
	    if strlowcase(colpick) eq 'blue'  then coluse=colors.blue
	    if strlowcase(colpick) eq 'green' then coluse=colors.green
	    if strlowcase(colpick) eq 'orange'then coluse=colors.orange
	    if strlowcase(colpick) eq 'cyan'  then coluse=colors.cyan
	    if strlowcase(colpick) eq 'red'   then coluse=colors.red
	endif else coluse=colors.black   
	shotstr = string(shot,format='(I5)')
    	filename = 'save/'+shotstr+'/'+los+'get_cn.sav'
	if keyword_set(aug)then begin
    	    read_signal_mrm,0L,shot,'UVS','D_tot',x1,y1,2,exp=exp
	    read_signal_mrm,0L,shot,'UVS','N_tot',x2,y2,2,exp=exp
	    cn_flux = (y2/7)/(y2/7+y1)
	    id = where(x1 ge min(data.time) and x1 le max(data.time))
	    cn_flux = cn_flux[id]
	    cn_time = x1[id]
	    cn_tdiv = interpol(data.tdiv,data.time,cn_time)
	    cn_idsort = sort(cn_tdiv)
	endif else begin
    	    if keyword_set(switchgas)then begin
    	    	ppfread,shot=shot,dda='GASM',dtype='MN1R',data=d_flux,t=time_d
    	    	ppfread,shot=shot,dda='GASM',dtype='MAJR',data=n_flux,t=time_n
	    endif else begin
    	    	ppfread,shot=shot,dda='GASM',dtype='MAJR',data=d_flux,t=time_d
    	    	ppfread,shot=shot,dda='GASM',dtype='MN1R',data=n_flux,t=time_n	    
	    end
	    cn_flux =( interpol(n_flux/7,time_n,time_d) / (d_flux+interpol(n_flux/7,time_n,time_d)))
	    cn_time =time_d
	    cn_tdiv = -1+fltarr(n_elements(cn_time))
	end   
    	if keyword_set(fluxonly)then begin
	    if ~keyword_set(overplot)then begin
	    	setgraphics,xs=800,ys=600,psplot=psplot,file='JET_1d_conc.ps' ,colors=colors & plot,cn_time-t0,cn_flux*100,xr=[0,2.5],back=colors.white,xtitle='Time - t!lseeding!n [s]',ytitle='c!lN!n [%]',col=colors.black,/nodata,yr=[0,50]
	    endif
	    oplot,cn_time-t0,cn_flux*100,col=coluse
	    stop
	end
	if keyword_set(load) and file_test(filename) then begin
	    restore,filename
	endif else begin
       	    load_jet_paths
	    
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
		plasma = optimize_ratios(data,shot,los,transmission,sm=sm,debug=debug,nocn=nocn)
	    endif else begin
		plasma = estimate_ratios(data,shot,los,transmission,sm=sm,debug=debug,lowerte=lowerte,upperte=upperte)
	    end
	    save,file=filename,plasma,data
	end
	x1=-1 & y1=-1 & x2=-1 & y2=-1
	!mouse.button=0
	iter = 0
	nrow = 3
	ncol = 1
	if n_elements(plasma.cn_upper[0,*]) gt 1 then begin
	    x1   = 1.0
	    x2   = 1.5
	    dens2d = plasma.cn_upper
	    temp2d = plasma.cn_upper
	    conc2d = plasma.cn_upper
	    for i=0,n_elements(plasma.cn_upper[0,*])-1 do begin
	    	for j=0,n_elements(plasma.cn_upper[*,0])-1 do begin
		    dens2d[j,i]=(plasma.dens_upper[j,i]+plasma.dens_lower[j,i])/2.0
		    temp2d[j,i]=(plasma.te_upper[j,i]+plasma.te_lower[j,i])/2.0
		    conc2d[j,i]=(plasma.cn_upper[j,i]+plasma.cn_lower[j,i])/2.0		    
		endfor
	    endfor
    	    rval = data.rvals    	    
	    dens2d = transpose(dens2d)/1e14 & conc2d = transpose(conc2d)*100 & emiss2d = transpose(data.nii4041)/1e18 & temp2d = transpose(temp2d)
	    time = data.time
	    if keyword_set(debug)then begin
    	    	setgraphics,xs=800,ys=600,psplot=psplot,file='JET_dens.ps' & shimage,dens2d,range=[0,5],rval,time,title='Density [10!u20!n m!u-3!n]',ytitle='Time [s]',xtitle='R [m]',xr=[2.6,2.9]
	    	setgraphics,xs=800,ys=600,psplot=psplot,file='JET_temp.ps' & shimage,temp2d,rval,time,range=[2,7],title='Temperature [eV]',ytitle='Time [s]',xtitle='R [m]',xr=[2.6,2.9]
	    	if ~keyword_set(nocn)then begin
	    	    setgraphics,xs=800,ys=600,psplot=psplot,file='JET_conc.ps' & shimage,conc2d,rval,time,range=[0,40],title='Concentration [%]',ytitle='Time [s]',xtitle='R [m]',xr=[2.6,2.9]
	    	endif
	    	setgraphics,xs=800,ys=600,psplot=psplot,file='JET_emis.ps' & shimage,emiss2d,rval,time,range=[0,max(emiss2d[5:*,*])*0.8],title='N II @ 404.1 nm [10!u18!n ph/s/m!u2!n/sr]',ytitle='Time [s]',xtitle='R [m]',xr=[2.6,2.9]
	    endif
	    id_range = where(rval ge 2.73 and rval le 2.78)
	    cn_data  = fltarr(n_elements(time))
	    cn_derr  = cn_data
	    dens  = fltarr(n_elements(time))
	    ne_err = cn_data
	    temp  = fltarr(n_elements(time))
	    te_err  = cn_data
	    
	    for i=0,n_elements(time)-1 do begin stat=moment(conc2d[id_range,i]) & cn_data[i] = stat[0] & cn_derr[i]=sqrt(stat[1]) & endfor
	    for i=0,n_elements(time)-1 do begin stat=moment(dens2d[id_range,i]) & dens[i] = stat[0] & ne_err[i]=sqrt(stat[1]) & endfor
	    for i=0,n_elements(time)-1 do begin stat=moment(temp2d[id_range,i]) & temp[i] = stat[0] & te_err[i]=sqrt(stat[1]) & endfor
	    if ~keyword_set(overplot)then setgraphics,xs=800,ys=600,psplot=psplot,file='JET_1d_conc.ps' ,colors=colors
	    while !mouse.button ne 4 do begin
    	    	plot,time-t0,cn_data,xr=[0,3],back=colors.white,xtitle='Time - t!lseeding!n [s]',ytitle='c!lN!n [%]',col=colors.black,/nodata,yr=[0,50] 
	    	oplot,time-t0,cn_data,col=coluse & errors,time-t0,cn_data,ystd=cn_derr,col=coluse
	    	oplot,time-t0,data.tdiv*10,col=colors.red
	    	if keyword_set(overplot)then yshift=-0.1*overplot
	    	legend,shotstr,coluse,yshift=yshift
	    	if iter ne 0 then cursor,x1,y1,/up
	    	if !mouse.button eq 2 then goto,quit
	    	oplot,[x1,x1],[0,100],col=colors.red
	    	if iter ne 0 then cursor,x2,y2,/up
	    	oplot,[x2,x2],[0,100],col=colors.red
	    	id      = where(time-t0 ge x1<x2 and time-t0 le x1>x2)	    	    
	    	stats   = moment(cn_data[id])	    
	    	cn_mean = mean(cn_data[id])
	     	cn_err  = mean(cn_derr[id])
	    
	    	cn_dens = mean(dens[id])
	    	cn_dens_err  = mean(ne_err[id])

	    	cn_te = mean(temp[id])
	    	cn_te_err  = mean(te_err[id])
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
	endif else begin
	    x1   = 3.0
	    x2   = 5.0
	    if keyword_set(debug)then setgraphics,nrow=nrow,ncol=ncol,colors=colors,xs=600,ys=1000,psplot=psplot,file='figures/cn_analysis.ps'
	
	    while !mouse.button ne 4 do begin
		ii   = 0
		exp_ratio1 = smooth(data.nii4041,sm,/edge_truncate) / smooth(data.nii3995,sm,/edge_truncate)
		exp_ratio2 = smooth(data.nii4041,sm,/edge_truncate) / smooth(data.nii4026,sm,/edge_truncate)
		rawtdiv    = interpol(data.tdiv,data.time,data.rawtime)
		idsort     = sort(data.tdiv)
		idsortraw  = sort(rawtdiv)
		xr = [min(data.time)*1.05,max(data.time)*0.95]
		xr = [0,25]
		if keyword_set(debug)then begin
			plot,data.tdiv[idsort],exp_ratio1[idsort],col=colors.black,back=colors.white,yr=[0.,0.3],/nodata,xr=xr,xs=1,$
			ytitle='N II line ratios [-]',position=graphpos(0,ii,nrow,ncol,xspc=xspc,yspc=yspc),xtickname=replicate(' ',30)
			oband,data.tdiv[idsort],plasma.ratio1_lower[idsort],plasma.ratio1_upper[idsort],/norm,col=colors.black
			oband,data.tdiv[idsort],plasma.ratio2_lower[idsort]/70,plasma.ratio2_upper[idsort]/70,/norm,col=colors.blue
			oplot,data.tdiv[idsort],exp_ratio1[idsort],col=colors.red
			oplot,data.tdiv[idsort],exp_ratio2[idsort]/70,col=colors.green
			legend,'404.1/399.5',colors.black,yshift=-0.07
			legend,'404.1/402.6/70',colors.black,yshift=-0.5
		    	dl = length(data.tdiv[idsort],shot,los,upperdl=upperdl,lowerdl=lowerdl)
		        ii = ii +1
			plot,data.tdiv[idsort],dl*100,col=colors.black,back=colors.white,yr=[0,10],$
			   position=graphpos(0,ii,nrow,ncol,xspc=xspc,yspc=yspc),xtickname=replicate(' ',30),xr=xr,xs=1,$
			   ytitle='Delta L [cm]'

		        oplot,[x1,x1],[0,10],col=colors.red,linest=5
		        oplot,[x2,x2],[0,10],col=colors.red,linest=5
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
			if iter ne 0 then cursor,x1,y1,/up
			if !mouse.button eq 2 then goto,quit
			oplot,[x1,x1],[0,40],col=colors.red,linest=5
			if iter ne 0 then cursor,x2,y2,/up
			oplot,[x2,x2],[0,40],col=colors.red,linest=5
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

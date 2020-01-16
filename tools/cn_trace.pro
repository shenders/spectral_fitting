Pro cn_trace,shot,$
             los=los,$
	     append=append,$
	     xr=xr,$
	     use_evl=use_evl,$
	     diag=diag,$
	     channel=channel,$
	     interelm=interelm,$
	     transmission=transmission,$
	     full=full,$
	     sm=sm,$
	     load=load,$
	     machine=machine

	if ~keyword_set(transmission)then transmission=1.0
	if ~keyword_set(sm)then sm=20.0
	if ~keyword_set(machine)then machine='AUG'
    	
	shotstr = string(shot,format='(I5)')
	filename = 'save/'+shotstr+'/'+los+'get_cn.sav'

	if keyword_set(load) and file_test(filename) then begin
	    restore,filename
	endif else begin

	    data  =  get_nii(shot,$
		 los=los,$
		 xr=xr,$
		 interelm=interelm,$
 		 channel=channel,$
		 use_evl=use_evl,$
		 append=append,$
		 diag=diag)
		 
 
       	; get cn, Te and dens from line ratios

            if keyword_set(full)then begin
		plasma = optimize_ratios(data,shot,los,transmission,sm=sm,debug=debug)
       	    endif else begin
	       	plasma = estimate_ratios(data,shot,los,transmission,sm=sm,debug=debug,lowerte=3.1,upperte=3.5)
       	    end
	end
	
	if machine eq 'JET' then begin
	    
    	    rval = data.rvals    	    
	    if ~keyword_set(horizontal)then id_range = where(rval ge 2.73 and rval le 2.78) else id_range = where(rval ge 2.65 and rval le 2.68)
	    dens2d = plasma.cn_upper
	    temp2d = plasma.cn_upper
	    conc2d = plasma.cn_upper
    	    emiss2d = transpose(data.nii4041)/1e18	    
	    for i=0,n_elements(plasma.cn_upper[0,*])-1 do begin
	    	for j=0,n_elements(plasma.cn_upper[*,0])-1 do begin
		    dens2d[j,i]=(plasma.dens_upper[j,i]+plasma.dens_lower[j,i])/2.0
		    temp2d[j,i]=(plasma.te_upper[j,i]+plasma.te_lower[j,i])/2.0
		    conc2d[j,i]=(plasma.cn_upper[j,i]+plasma.cn_lower[j,i])/2.0		    
		endfor
	    endfor
	    dens2d = transpose(dens2d)/1e14 & conc2d = transpose(conc2d)*100 & emiss2d = transpose(data.nii4041)/1e18 & temp2d = transpose(temp2d)
	    time = data.time	    
	    cn_mean  = fltarr(n_elements(time))
	    temp = fltarr(n_elements(time)) & dens = fltarr(n_elements(time)) & cn_mean = fltarr(n_elements(time))
	    te_err = fltarr(n_elements(time)) & ne_err = fltarr(n_elements(time)) & cn_err = fltarr(n_elements(time))
	    
	    for i=0,n_elements(time)-1 do begin stat=moment(conc2d[id_range,i]) & cn_mean[i] = stat[0] & cn_err[i]=sqrt(stat[1]) & endfor
	    for i=0,n_elements(time)-1 do begin stat=moment(dens2d[id_range,i]) & dens[i] = stat[0] & ne_err[i]=sqrt(stat[1]) & endfor
	    for i=0,n_elements(time)-1 do begin stat=moment(temp2d[id_range,i]) & temp[i] = stat[0] & te_err[i]=sqrt(stat[1]) & endfor
	    
	endif else begin
	    
	    time = data.time
	    temp = fltarr(n_elements(time)) & dens = fltarr(n_elements(time)) & cn_mean = fltarr(n_elements(time))
	    te_err = fltarr(n_elements(time)) & ne_err = fltarr(n_elements(time)) & cn_err = fltarr(n_elements(time))
	    for i=0,n_elements(time)-1 do temp[i] = (plasma.te_upper[i]+plasma.te_lower[i])/2.0
	    for i=0,n_elements(time)-1 do te_err[i] = temp[i]-plasma.te_lower[i]
	    for i=0,n_elements(time)-1 do dens[i] = (plasma.dens_upper[i]+plasma.dens_lower[i])/2.0
	    for i=0,n_elements(time)-1 do ne_err[i] = dens[i]-plasma.dens_lower[i]
	    for i=0,n_elements(time)-1 do cn_mean[i] = (plasma.cn_upper[i]+plasma.cn_lower[i])/2.0
	    for i=0,n_elements(time)-1 do cn_err[i] = cn_mean[i]-plasma.cn_lower[i]
	    dens = dens / 1e14	  
	    cn_mean = cn_mean * 100
	    cn_err = cn_err * 100  
	end
	nrow = 2
	ncol = 1
	xspc = 0.1
	yspc = 0.1
	setgraphics,nrow=nrow,ncol=ncol,colors=colors,xs=800,ys=900
	
	plot,time,temp,col=colors.black,back=colors.white,$
	xr=xr,xs=1,ys=9,yr=[0,6],ytitle='T!le !n[eV]',xtitle='Time [s]',pos=graphpos(0,1,nrow,ncol,xspc=xspc,yspc=yspc)
	errors,time,temp,ystd=te_err,col=colors.black
	axis,yaxis=1,/save,yr=[0,3],ytitle='n!le !n[10!u20!n m!u-3!n]',col=colors.blue
	oplot,time,dens,col=colors.blue
	errors,time,dens,ystd=ne_err,col=colors.blue
	
	plot,time,cn_mean,col=colors.black,back=colors.white,$
	xr=xr,xs=1,ys=1,yr=[0,max(cn_mean)],ytitle='c!lN !n[%]',xtitle='Time [s]',pos=graphpos(0,0,nrow,ncol,xspc=xspc,yspc=yspc)
    	errors,time,cn_mean,ystd=cn_err,col=colors.black

	stop
End

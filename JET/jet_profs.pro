Pro jet_profs,shot,dda=dda,tr=tr,xr=xr,shift=shift,$
    	      tesep = tesep, nesep=nesep, err=err,$
	      debug=debug

    if ~keyword_set(shift)then begin
    	shift = 0.02
	if shot eq 85262 then shift = 0.023
	if shot eq 85264 then shift = 0.026
	if shot eq 85265 then shift = 0.026
	if shot eq 85266 then shift = 0.027
	if shot eq 85267 then shift = 0.027
	if shot eq 85268 then shift = 0.027
	if shot eq 85270 then shift = 0.027
	if shot eq 85272 then shift = 0.025
	if shot eq 85274 then shift = 0.025
	if shot eq 85276 then shift = 0.026
	if shot eq 96075 then shift = 0.0475
	if shot eq 85423 then shift = 0.0395
	if shot eq 85417 then shift = 0.041
	if shot eq 96227 then begin
	    shift = 0.030
	    dda = 'EFTP'
	endif
    endif
    
    if ~keyword_set(tr)then tr=[40,60]
    if ~keyword_set(xr)then xr=[0,1.1]

    ; read in equilibrium
    if ~keyword_set(dda)then dda='EFIT'
    ppfread,shot=shot,dda=dda,dtype='PSI',data=psirz,t=time,x=rad
    ppfread,shot=shot,dda=dda,dtype='PSIR',data=psir
    ppfread,shot=shot,dda=dda,dtype='PSIZ',data=psiz
    ppfread,shot=shot,dda=dda,dtype='RMAG',data=rmag,t=t1
    ppfread,shot=shot,dda=dda,dtype='ZMAG',data=zmag,t=t2
    ppfread,shot=shot,dda=dda,dtype='FAXS',data=axes,t=t3
    ppfread,shot=shot,dda=dda,dtype='FBND',data=bndy,t=t4
    nr  = n_elements(psir)
    nz  = n_elements(psiz)
    nt  = n_elements(time)
    psi = fltarr(nr,nz,nt)
    for t=0,nt-1 do begin
        for i=0,n_elements(psiz)-1 do begin
    		psi[*,i,t] = psirz[i*nr:(i+1)*nr-1,t]
    	endfor
    endfor
    
    midplane = fltarr(nt)
    psin     = fltarr(nr,nt)
    for t=0,nt-1 do begin
    	midplane  = where(abs(psiz-zmag[t]) eq min(abs(psiz-zmag[t])))
        psin[*,t] = (psi[*,midplane,t] - axes[t])/(bndy[t] - axes[t])
    endfor 
    
    ; read in temperature
    ppfread,shot=shot,dda='HRTS',dtype='TE',data=temp,t=t,x=rad

    ; read in density
    ppfread,shot=shot,dda='HRTS',dtype='NE',data=dens,t=t,x=rad
    
    if keyword_set(shift)then rad = rad + shift

    ; convert to psiN
    psi_hrts = fltarr(n_elements(rad),n_elements(t))
    for i=0,n_elements(t)-1 do begin
    	id = where(abs(time-t[i]) eq min(abs(time-t[i])))
    	psi_hrts[*,i] = interpol(psin[*,id[0]],psir,rad)
    endfor
    
    if keyword_set(debug)then setgraphics,nrow=2,ncol=1,xs=800,ys=800,colors=colors
    lims = [0.99,1.01]
    id = where(t ge tr[0] and t le tr[1])
    ne_mean = fltarr(n_elements(id))
    te_mean = fltarr(n_elements(id))
    if keyword_set(debug)then plot,psi_hrts[*,0],temp[*,0],/nodata,col=colors.black,back=colors.white,xr=xr,yr=[0,1500],xs=1,$
    xtitle='Psi!lN!n',ytitle='T!le!n [eV]',title=string(shot,format='("JPN #",i5)')
    for i=min(id),max(id) do begin
    	if keyword_set(debug)then oplot,psi_hrts[*,i],temp[*,i],col=colors.black
	vals = interpol(temp[*,i],psi_hrts[*,i],findgen(10)*(lims[1]-lims[0])/9.0+lims[0])
	stat = moment(vals)
    	te_mean[i-min(id)]=stat[0]
    endfor
    if keyword_set(debug)then oplot,[0,100],[mean(te_mean),mean(te_mean)],col=colors.red,linest=5
    if keyword_set(debug)then oplot,[0,100],[100,100],col=colors.green,linest=5
    if keyword_set(debug)then oplot,[1,1],[0,100000],col=colors.red,linest=5
    if keyword_set(debug)then plot,psi_hrts[*,0],dens[*,0],/nodata,col=colors.black,back=colors.white,xr=xr,yr=[0,10],xs=1,$
    xtitle='Psi!lN!n',ytitle='n!le!n [10!u19!nm!u-3!n]'
    for i=min(id),max(id) do begin
    	if keyword_set(debug)then oplot,psi_hrts[*,i],dens[*,i]/1e19,col=colors.black
	vals = interpol(dens[*,i]/1e19,psi_hrts[*,i],findgen(10)*(lims[1]-lims[0])/9.0+lims[0])
	stat = moment(vals)
    	ne_mean[i-min(id)]=stat[0]
    endfor
    
    if keyword_set(debug)then oplot,[0,100],[mean(ne_mean),mean(ne_mean)],col=colors.red,linest=5
    if keyword_set(debug)then oplot,[1,1],[0,100000],col=colors.red,linest=5
    if keyword_set(debug)then print,'Average Te [eV]        = ',mean(te_mean)
    if keyword_set(debug)then print,'Average ne [10^19 m-3] = ',mean(ne_mean)
    
    
    tesep = mean(te_mean)
    nesep = mean(ne_mean)
    
    ; Use standard error of mean
    
    err   = (1.0/sqrt(n_elements(ne_mean)))* sqrt((moment(ne_mean))[1])

end

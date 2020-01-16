Pro jet_get_flux,shots,twn1,twn2,flux_cn,flux_cn_err
    
    flux_cn = fltarr(n_elements(shots))
    flux_cn_err = flux_cn
    for i=0,n_elements(shots)-1 do begin
    	shot = shots[i]
    	if keyword_set(switchgas)then begin
    	    ppfread,shot=shot,dda='GASM',dtype='MN1R',data=d_flux,t=time_d
    	    ppfread,shot=shot,dda='GASM',dtype='MAJR',data=n_flux,t=time_n
	endif else begin
    	    ppfread,shot=shot,dda='GASM',dtype='MAJR',data=d_flux,t=time_d
    	    ppfread,shot=shot,dda='GASM',dtype='MN1R',data=n_flux,t=time_n	
	end
	conc =( interpol(n_flux/7,time_n,time_d) / (d_flux+interpol(n_flux/7,time_n,time_d)))
	time =time_d
    	id = where(time ge twn1[i] and time le twn2[i])
	flux_cn[i] = mean(conc[id]) * 100.0
	flux_cn_err[i] = sqrt((moment(conc[id]))[1]) * 100.0
	debug=0
	if keyword_set(debug) then begin
	    plot,time,conc*100,yr=[0,50],xr=[45,57]
	    oplot,[twn1[i],twn1[i]],[0,100],linest=5
	    oplot,[twn2[i],twn2[i]],[0,100],linest=5
	    cursor,x,y,/up
	endif    
    endfor
End


Pro jet_db,psplot=psplot

    ; list shots
    shots = [85417, 85264, 85265, 85266, 85267, 85268, 85270, 85272, 85274, 85276, 96075, 85423, 96227] 
    twn1  = [53.0 , 53.0 , 53.0 , 53.0 , 53.0 , 53.0 , 53.0 , 53.0 , 53.0 , 53.0 , 49.0 , 50.0 , 50.0 ]
    twn2  = [55.0 , 55.0 , 55.0 , 55.0 , 55.0 , 55.0 , 55.0 , 55.0 , 55.0 , 55.0 , 50.5 , 52.0 , 52.0 ]
    psep  = -1
    tsep  = -1
    pinj  = -1
    nsep  = -1
    icur  = -1
    kapp  = -1
    btor  = -1
    amin  = -1
    cn_mean = -1
    cn_err  = -1
    psep_err  = -1
    icur_err  = -1
    nsep_err  = -1
    pinj_err  = -1

   
    ; fetch machine data
    for i=0,n_elements(shots)-1 do begin
    	t=-1 & data=-1 & ppfread,shot=shots[i],dda='MAGN',dtype='IPLA',data=data,t=t & id = where(t ge twn1[i] and t le twn2[i])
    	icur = [icur,mean(data[id])]
    	icur_err = [icur_err,sqrt((moment(data[id]))[1])]
        t=-1 & data=-1 & ppfread,shot=shots[i],dda='EFIT',dtype='ELON',data=data,t=t & id = where(t ge twn1[i] and t le twn2[i])
	kapp = [kapp,mean(data[id])]
        t=-1 & data=-1 & ppfread,shot=shots[i],dda='SCAL',dtype='BT',data=data,t=t & id = where(t ge twn1[i] and t le twn2[i])
	btor = [btor,mean(data[id])]
        t=-1 & data=-1 & ppfread,shot=shots[i],dda='SCAL',dtype='AMIN',data=data,t=t & id = where(t ge twn1[i] and t le twn2[i])
	amin = [amin,mean(data[id])]
        t=-1 & data=-1 & ppfread,shot=shots[i],dda='SCAL',dtype='PIN',data=data,t=t & id = where(t ge twn1[i] and t le twn2[i]) 
	pinj = [pinj,mean(data[id])]
	pinj_err = [pinj_err,sqrt((moment(data[id]))[1])]
    endfor

    btor = btor[1:*]
    kapp = kapp[1:*]
    amin = amin[1:*]
    icur = icur[1:*]
    pinj = pinj[1:*]
    pinj_err = pinj_err[1:*]
    
    ; fetch Psep
    for i=0,n_elements(shots)-1 do begin
     	jet_radiated_power,shots[i],time,pr_main=prad,pr_div=pdiv,pr_tot=ptot
	stats = moment(prad[where(time ge twn1[i] and time le twn2[i])])
	psep = [psep,pinj[i] - stats[0]]
	psep_err = [psep_err,sqrt(stats[1] + pinj_err[i]^2)]
    endfor
    psep = psep[1:*]/2.6
    psep_err = psep_err[1:*]/2.6

    ; get own estimate of ne,sep
    for i=0,n_elements(shots)-1 do begin
    	jet_profs,shots[i],nesep=nesep,err=err,tesep=tesep
	nsep = [nsep,nesep]
	tsep = [tsep,tesep]
	nsep_err = [nsep_err,err]
    endfor
    tsep = tsep[1:*]
    nsep = nsep[1:*]
    nsep_err = nsep_err[1:*]

    ; fetch concentrations
    for i=0,n_elements(shots)-1 do begin
    	restore,string(shots[i],format='("save/cn/",i5,"/cn_database_2.sav")')
    	cn_mean = [cn_mean,conc]
	cn_err  = [cn_err,err]
    endfor
    cn_mean = cn_mean[1:*]
    cn_err  = cn_err[1:*]

    ; get flux concentrations
    jet_get_flux,shots,twn1,twn2,flux_cn,flux_cn_err
     
    ;fac = 5.92
    ;if keyword_set(redpsep)then fac=fac; /10;2
    amin_alpha = -2.4
    scaling = (psep/1e6)^(0.999) * (abs(icur)/1e6)^(1.04) * nsep^(-2.63) * (1.0+kapp^2)^(-1.0) * amin^(amin_alpha)
    scaling_err = scaling * sqrt( (0.999 * Psep_err / Psep)^2 + (2.63 * nsep_err / nsep)^2)
    goldston = (psep/1e6) * (abs(icur)/1e6) * nsep^(-2) * (1.0+kapp^2)^(-1.0) * amin^(amin_alpha)
    id = where(cn_mean gt 20.0)
    fac = cn_mean[id[0]]/ scaling[id[0]]
    fac = 15.5
    scaling = scaling * fac
    scaling_err = scaling_err * fac

    setgraphics,colors=colors,xs=800,ys=600,psplot=psplot,file='cn_flux.ps'
    fac_flux = 0.5
    plot,flux_cn * fac_flux,cn_mean,col=colors.black,back=colors.white,psym=5,xr=[0,30],yr=[0,30],$
    xtitle=string(fac_flux,format='("Flux c!lN!n x",f3.1," [%]")'),ytitle='Spectroscopy c!lN!n [%]'
    errors,flux_cn * fac_flux,cn_mean,ystd=cn_err,xstd=flux_cn_err * fac_flux,col=colors.black
    oplot,[0,100],[0,100],linest=5,col=colors.gray

    setgraphics,colors=colors,xs=800,ys=600,psplot=psplot,file='figures/cn_scal.ps'
    
    user_psym,1,/fill & plot,scaling,cn_mean,col=colors.black,back=colors.white,psym=5,xr=[0,30],yr=[0,30],$
    xtitle='Scaling c!lN!n [%]',ytitle='Spectroscopy c!lN!n [%]'
    oplot,[0,100],[0,100],linest=5,col=colors.gray
    errors,scaling,cn_mean,ystd=cn_err,xstd=scaling_err,col=colors.black
    print,'       Shot       cN           scaling      psep         nsep         Tsep         Ip           flux cN     '   
    for i=0,n_elements(shots)-1 do print,shots[i],cn_mean[i],scaling[i],psep[i]/1e6,nsep[i],tsep[i],abs(icur[i])/1e6,flux_cn[i]
    
    cn_mean_jet = cn_mean
    regress_jet = scaling
    regress_jet_err = scaling_err

    restore,'save/database_scaling_aug.sav'
    user_psym,5,/fill
    cn_mean_aug = cn_mean
    regress_aug     = fac * (Psep/1e6)^(0.999)  * (dens/1e19)^(-2.63) * (pcur/1e6)^(1.04) * (1+k^2)^(-1.0) * aminor^(amin_alpha)
    regress_err_aug = regress_aug * sqrt( (0.999 * Psep_err / Psep)^2 + (2.63 * dens_err / dens)^2)

    oplot,regress_aug,cn_mean_aug,col=colors.red,psym=8
    errors,regress_aug,cn_mean_aug,ystd=cn_err,xstd=regress_err,col=colors.red
    regress = [regress_jet,regress_aug]
    cn_mean = [cn_mean_jet,cn_mean_aug]
    rsqr = r2(regress,cn_mean,b=b,m=m)
    legend,'R!u2!n='+string(rsqr,format='(d4.2)'),colors.black,yshift=-0.1,xshift=-0.6
    legend,string(b,m,format='("y=",d5.2,"+",d4.2,"x")'),colors.black,yshift=0.0,xshift=-0.6

    ; recompute Goldston's table
    psep_device = [3.83, 10.7, 14.0, 100.0] / 2.5
    ip_device   = [0.82, 1.2 , 2.5 , 15.0 ]
    ngw_device  = [5.39, 1.44, 0.98, 1.19 ] * 10
    k95_device  = [1.51, 1.63, 1.73, 1.80 ]
    a_device    = [0.22, 0.52, 0.90, 2.00 ]
    device      = ['CMOD','AUG','JET','ITER']
    cn_device   = 15.5 * psep_device * (0.5 * ngw_device)^(-2.63) * ip_device * (1+k95_device^2)^(-1.0) * a_device^(-2.4)
    print,device
    
    print,'cN: ',cn_device
    
    
    cz = findgen(100)*3.0/99.0
    yaxis_n = ((2-cz * (7-1.0))/(2.0*(1.0-7*cz)+14.0*cz))^(0.5)
    yaxis_ne = ((2-cz * (10-1.0))/(2.0*(1.0-10*cz)+20.0*cz))^(0.5)
    
    
    
    stop
End

Function retrieve,shot,t0=t0,cont=cont,debug=debug

    shotstr = string(shot,format='(i5)')
    trace='save/'+shotstr+'/SP-data.idl'
    restore,trace[0]
    time = output.time
    if keyword_set(t0)then time = time-t0
    nii  = output.nii[*,*,1]
    sp   = fltarr(n_elements(output.los_names)) & for i=0,n_elements(output.los_names)-1 do sp[i]  = float(strmid(output.los_names[i],2,2))

    ; Average strike-point and x-point areas
    ; note there is a difference in fibre location between 85### and 96###
    if shot gt 96000 then begin
        ids    = [17,18,19]
    	idx    = [13,14,15]
    	idxx   = [7,8,9]
    endif else begin
        ids    = [16,17,18]
    	idx    = [12,13,14]    
    	idxx   = [6,7,8]    
    end
    id_sp  = fltarr(n_elements(ids))
    id_xp  = fltarr(n_elements(idx))
    id_xxp  = fltarr(n_elements(idx))
    avr_sp = fltarr(n_elements(time))
    avr_xp = fltarr(n_elements(time))
    avr_xxp = fltarr(n_elements(time))
    
    for i=0,n_elements(ids)-1  do id_sp[i]  = where(sp eq ids[i])
    for i=0,n_elements(idx)-1  do id_xp[i]  = where(sp eq idx[i])
    for i=0,n_elements(idxx)-1  do id_xxp[i]  = where(sp eq idxx[i])
    for i=0,n_elements(time)-1 do avr_sp[i] = mean(nii[i,id_sp])     
    for i=0,n_elements(time)-1 do avr_xp[i] = mean(nii[i,id_xp])     
    for i=0,n_elements(time)-1 do avr_xxp[i] = mean(nii[i,id_xxp])     
    ratio  = avr_sp / avr_xp
    
    if keyword_set(debug)then begin
    	ppfread,shot=shot,dda='GASM',dtype='MAJR',data=data2,t=t2 
	n2 = interpol(data2,t2-t0,time)
	setgraphics,colors=colors,colpick=colpick,/full,psplot=psplot,xs=800,ys=600    
    	id = where(time le 2.7)
	plot,time[id],smooth(avr_xp[id],4,/edge_trun)/max(avr_xp[id])/0.9,/nodata,xr=[0,2.7],xs=1,background=colors.white,col=colors.black,yr=[0,1.1],ys=1
	user_psym,3,/fill
	oplot,time[id],avr_xp[id]/max(avr_xp[id])/0.9,col=colors.black
	oplot,time[id],avr_sp[id]/max(avr_sp[id])/0.9,col=colors.red
	oplot,time[id],avr_xxp[id]/max(avr_xxp[id]),col=colors.blue
	oplot,time[id],n2[id]/max(n2[id]),col=colors.green
	stop
    endif
    
    
    if keyword_set(cont)then begin
    	
	setgraphics,colors=colors,colpick=colpick,/full,psplot=psplot,xs=800,ys=600
	shimage,nii[*,2:*],time,sp[2:*],xtitle='Time [s]',ytitle='SP##'
    endif
    return,{time:time,$
            nii:nii,$
    	    ratio:ratio}
End
Pro get_jet_deltal,shots,psplot=psplot,cont=cont,debug=debug
    
    x0 = retrieve(85417,t0=52)
    x1 = retrieve(85264,t0=52)
    x2 = retrieve(85265,t0=52)
    x3 = retrieve(85266,t0=52)
    x4 = retrieve(85267,t0=52)
    x5 = retrieve(85268,t0=52)
    x6 = retrieve(85270,t0=52)
    x7 = retrieve(85272,t0=52)
    x8 = retrieve(85274,t0=52)
    x9 = retrieve(85276,t0=52)
    x10= retrieve(96075,t0=49)
    x11= retrieve(96227,t0=49)
    setgraphics,colors=colors,colpick=colpick,/full,psplot=psplot,xs=800,ys=600,file='figures/ratios_detachment.ps'

    plot,x1.time,smooth(x1.ratio,2),yr=[0,2],col=colors.black,back=colors.white,xr=[0,2],$
           xtitle='Time -t!lseeding!n [s]',ytitle='N II: Strike-point/X-point',/nodata
    oplot,x0.time,smooth(x0.ratio,2),col=colors.pink
    oplot,x1.time,smooth(x1.ratio,2),col=colors.red
    oplot,x2.time,smooth(x2.ratio,2),col=colors.blue
    oplot,x3.time,smooth(x3.ratio,2),col=colors.orange
    oplot,x4.time,smooth(x4.ratio,2),col=colors.green
    oplot,x5.time,smooth(x5.ratio,2),col=colors.cyan
    oplot,x6.time,smooth(x6.ratio,2),col=colors.magenta
    oplot,x7.time,smooth(x7.ratio,2),col=colors.sky
    oplot,x8.time,smooth(x8.ratio,2),col=colors.gray
    oplot,x9.time,smooth(x9.ratio,2),col=colors.navy
    oplot,x10.time,smooth(x10.ratio,2),col=colors.gold
    oplot,x11.time,smooth(x11.ratio,2),col=colors.aqua


End

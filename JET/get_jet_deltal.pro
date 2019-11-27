
Pro get_jet_deltal
	
    ; highly detached plasma
    shot=96228
    shotstr = string(shot,format='(i5)')
    trace='save/'+shotstr+'/SP-data.idl'
    restore,trace[0]
    time1     = output.time
    nii1 = output.nii[*,*,1]
    sp1 = fltarr(n_elements(output.los_names)) & for i=0,n_elements(output.los_names)-1 do sp1[i]  = float(strmid(output.los_names[i],2,2))
    rat1 = nii1[*,where(sp1 eq 15)]/nii1[*,where(sp1 eq 8)]
    t0_1 = 49.5
    
    shot=96081
    shotstr = string(shot,format='(i5)')
    trace='save/'+shotstr+'/SP-data.idl'
    restore,trace[0]
    time2     = output.time
    nii2 = output.nii[*,*,1]
    sp2 = fltarr(n_elements(output.los_names)) & for i=0,n_elements(output.los_names)-1 do sp2[i]  = float(strmid(output.los_names[i],2,2))
    rat2 = nii2[*,where(sp2 eq 15)]/nii2[*,where(sp2 eq 8)]
    t0_2 = 48.0
    
    ; partially detached plasma
    shot=96073
    shotstr = string(shot,format='(i5)')
    trace='save/'+shotstr+'/SP-data.idl'
    restore,trace[0]
    time3     = output.time
    nii3 = output.nii[*,*,1]
    sp3 = fltarr(n_elements(output.los_names)) & for i=0,n_elements(output.los_names)-1 do sp3[i]  = float(strmid(output.los_names[i],2,2))
    rat3 = nii3[*,where(sp3 eq 15)]/nii3[*,where(sp3 eq 8)]
    t0_3 = 49.5

    ; partially detached plasma
    shot=96227
    shotstr = string(shot,format='(i5)')
    trace='save/'+shotstr+'/SP-data.idl'
    restore,trace[0]
    time4     = output.time
    nii4 = output.nii[*,*,1]
    sp4 = fltarr(n_elements(output.los_names)) & for i=0,n_elements(output.los_names)-1 do sp4[i]  = float(strmid(output.los_names[i],2,2))
    rat4 = nii4[*,where(sp4 eq 15)]/nii4[*,where(sp4 eq 8)]
    t0_4 = 50.5

    setgraphics,colors=colors,colpick=colpick,/full,psplot=psplot,xs=800,ys=600,file='figures/ratios_db.ps'

    plot,time1-t0_1,smooth(rat1,4),yr=[0,2.5],col=colors.black,back=colors.white,xr=[0,3]
    oplot,time2-t0_2,smooth(rat2,4),col=colors.blue
    oplot,time3-t0_3,smooth(rat3,4),col=colors.red
    oplot,time4-t0_4,smooth(rat4,4),col=colors.orange

stop

End

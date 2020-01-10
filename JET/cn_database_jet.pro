Pro info,debug=debug,pickshot=pickshot,psplot=psplot,flux=flux

    	shot_db = [ 85417    , 85418    , 85265    , 85270    , 85272    , 85274            ]
	
	twindow = [ [4.0,4.25],[3.5,4.2], [3.5,4.4], [3.5,5.5], [4.0,6.0], [3.0,5.2]        ]

		    
	tran    = [ 1.0      , 1.0      , 1.0      ,  1.0     , 1.0      , 1.0              ]
	
	full    = [ 1.0      , 1.0      , 1.0      , 1.0      , 1.0      , 1.0              ]	    
	
    	interelm= 0
	nshots  = n_elements(shot_db)
    	
	upperte = fltarr(nshots)+3.5
	lowerte = fltarr(nshots)+3.1
		
	sm      = 20	

	for i=0,nshots-1 do begin
		if keyword_Set(pickshot)then begin
			if pickshot ne shot_db[i] then goto,skip
		endif

		xr = [twindow[0,i],twindow[1,i]]
		if keyword_set(flux)then begin
			save,file='save/'+string(shot_db[i],format='(i5)')+'/flux_cn.sav',flux,flux_err,denom2,denom2_err,denom,denom_err
			goto,skip
		endif
		get_cn,	shot_db[i],los[i],$
                	debug=debug,$
			interelm=interelm[i],$
			use_evl=use_evl[i],$
			channel=channel[i],$
			diag=diag[i],$
			full=full[i],$
			xr=xr,$
			elmcond=elmcond[i],$
			transmission=tran[i],$
			sm=sm,$
			lowerte=lowerte[i],$
			upperte=upperte[i],psplot=psplot
		skip:
	endfor
	
End
Pro cn_database_jet,debug=debug,shot=shot,flux=flux,psplot=psplot

    	shots = [85417, 85418, 85265, 85270, 85272, 85274, 96227]

	if keyword_set(shot)then shots = shot 
	for i=0,n_elements(shots)-1 do info,debug=debug,psplot=psplot,pickshot=shots[i],flux=flux

End

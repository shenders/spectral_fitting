Pro info,debug=debug,pickshot=pickshot,psplot=psplot,flux=flux

	shot_db = [ 30298    , 30306    , 30307    , 30554    , 30623    , 30505    , 30506 ,$
	            30776    , 34108    , 34358    , 34359    , 34368    , 34971    , 32273 , $
		    34973    , 35158    , 35356    , 32244    , 34281    , 35846    , 32932 ]
	
	twindow = [ [4.0,4.25],[3.5,4.2], [3.5,4.4], [3.5,5.5], [4.0,6.0], [3.0,5.2], [2.2,3.1],$
	            [2.5,4.0], [3.5,4.6], [3.5,6.0], [3.5,5.1], [3.5,5.7], [3.0,3.9], [3.0,4.5],$
		    [3.0,3.8], [1.9,3.0], [2.0,3.5], [2.0,3.8], [3.5,4.4], [3.5,6.0], [5.0,6.0] ]

	los     = [ 'ROV014' , 'ROV014' , 'ROV014' , 'ROV014' , 'ROV014' , 'ROV014' , 'ROV014',$
	            'ROV014' , 'ROV-12' , 'ROV-13' , 'ROV-13' , 'ROV-13' , 'ROV-13' , 'ROV014',$
		    'ROV-13' , 'ROV-14' , 'ROV-14' , 'ROV014' , 'ROV-13' , 'ROV-14' , 'ROV014']
		    
	use_evl = [ 0        , 0        , 0        , 0        , 0        , 0        , 0   ,$
	            0        , 0        , 1        , 1        , 1        , 1        , 0   ,$
		    1        , 1        , 1        , 0        , 1        , 1        , 0 ]	    
	
	channel = [ -1       , -1       ,  -1      , -1       , -1       , -1       , -1  ,$
	            -1       , -1       , 19       , 19       , 19       , 19       , -1  ,$
		    19       , 22       , 22       , -1       , 19       , 9        , -1]
		    
	diag    = [ 'EVL'    , 'EVL'    , 'EVL'    , 'EVL'    , 'EVL'    , 'EVL'    , 'EVL', $
	            'EVL'    , 'FVL'    , 'EVL'    , 'EVL'    , 'EVL'    , 'EVL'    , 'EVL', $
		    'EVL'    , 'FVL'    , 'FVL'    , 'EVL'    , 'EVL'    , 'EVL'    , 'EVL' ]
		    
	tran    = [ 1.0      , 1.0      , 1.0      ,  1.0     , 1.0      , 1.0      , 1.0 ,$
	            1.0      , 1.0/0.64 , 1.0/0.67 , 1.0/0.67 , 1.0/0.67 , 1.0/0.67 , 1.0 ,$
		    1.0/0.67 , 1.0      , 1.0      , 1.0      , 1.0/0.67 , 1.0      , 1.0]
	
	full    = [ 1.0      , 1.0      , 1.0      , 1.0      , 1.0      , 1.0      , 1.0 ,$
	            1.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 1.0 ,$
		    0.0      , 0.0      , 0.0      , 1.0      , 0.0      , 0.0      , 1.0]	    
	
	elmcond = [ 4.0      , 4.0      , 4.0      , 3.0      , 4.0      , 3.0      , 3.0 ,$
	            4.0      , 4.5      , 4.5      , 2.0      , 2.0      , 3.0      , 3.0 ,$ 
	            2.0      , 3.0      , 4.5      , 4.5      , 3.0      , 3.0      , 3.0]
		  
	interelm= [ 0.0      , 0.0      , 0.0      , 1.0      , 1.0      , 1.0      , 1.0 ,$
	            1.0      , 0.0      , 0.0      , 0.0      , 0.0      , 0.0      , 1.0 ,$ 
	            0.0      , 1.0      , 1.0      , 1.0      , 1.0      , 1.0      , 1.0]

	nshots  = n_elements(shot_db)

	upperte = fltarr(nshots)+3.5
	lowerte = fltarr(nshots)+3.1
	
	;upperte[where(shot_db eq 35158)] = 3.1
	;lowerte[where(shot_db eq 35158)] = 3.5
	
	sm      = 20	

	for i=0,nshots-1 do begin
		if keyword_Set(pickshot)then begin
			if pickshot ne shot_db[i] then goto,skip
		endif

		xr = [twindow[0,i],twindow[1,i]]
		if keyword_set(flux)then begin
			time1= -1
			d_tot  = -1
			n_tot  = -1
			read_signal_mrm,0L,shot_db[i],'UVS','D_tot',time1,d_tot,2
			read_signal_mrm,0L,shot_db[i],'UVS','N_tot',time1,n_tot,2
			id = where(time1 ge xr[0] and time1 le xr[1])
			flux = (n_tot/7.0) *100.0/ (d_tot + n_tot/7.0)
			flux_dens =  (d_tot + n_tot/7.0)/1e22
			flux_dens2 =  (d_tot + n_tot)/1e22
			tmp = moment(flux[id])
			flux = tmp[0]
			flux_err = sqrt(tmp[1])	
			tmp = moment(flux_dens[id])
			denom = tmp[0]
			denom_err =sqrt(tmp[1])
			tmp = moment(flux_dens2[id])
			denom2 = tmp[0]
			denom2_err =sqrt(tmp[1])

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
Pro cn_database,debug=debug,shot=shot,flux=flux,psplot=psplot

	shots= [32244	   , 34368   , 30554	, 30776    ,  34358    ,30307	   ,30306   ,30298    ,34971	  , 34973    ,35158	, 34359    , 30505    , 30506	    ]
	if keyword_set(shot)then shots = shot 
	for i=0,n_elements(shots)-1 do info,debug=debug,psplot=psplot,pickshot=shots[i],flux=flux

End

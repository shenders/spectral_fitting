Pro iaea_nii

	shot_db = [ 30623    , 34108    , 34358	   , 30776    , 30306    ,$
	            35158    , 34973    , 35356    , 30554    ,$
		    30298    , 30307    ] 
	twindow = [ [4.9,6.0], [3.5,4.6], [4.0,6.0], [3.1,4.1], [3.5,4.5],$
	            [1.9,2.5], [3.4,3.9], [2.0,3.5], [3.6,5.5],$
		    [3.5,4.5], [3.5,4.5]]
	interelm= [ 1.0      , 1.0      , 1.0      , 1.0      , 1.0      ,$
	            1.0      , 1.0      , 1.00     , 1.0      ,$
		    1.0      , 1.0      ]
	useevl  = [ 0        , 1        , 1	   , 0	      , 0        ,$
	            1	     , 1	, 1        , 0        ,$
		    0        , 0        ]

	los12   = [ 'ROV012' , 'ROV-12' , 'ROV-11' , 'ROV012' , 'ROV012' ,$
	            'ROV-12' , 'ROV-11' , 'ROV-12' , 'ROV012' ,$
		    'ROV012' , 'ROV012' ] 
	tran12  = [ 1.0      , 1.0/0.64 , 1.0/0.61 , 1.0      , 1.0      ,$
	            1.0      , 1.0/0.61 , 1.00     , 1.0      ,$
		    1.0      , 1.0      ]
	chan12  = [ 0        , 4        , 22	   , 0	      , 0        ,$
	            22       , 22	, 22       , 0        ,$
		    0        , 0        ]
	diag12  = [ 'EVL'    ,'FVL'     , 'EVL'	   , 'EVL'    , 'EVL'    ,$
	            'EVL'    ,'EVL'     , 'EVL'    , 'EVL'    ,$
		    'EVL'    ,'EVL'     ] 

	los14   = [ 'ROV014' , 'ROV-14' , 'ROV-13' , 'ROV014' , 'ROV014' ,$
	            'ROV-14' , 'ROV-13' , 'ROV-14' , 'ROV014' ,$
		    'ROV014' , 'ROV014' ] 
	tran14  = [ 1.0      , 1.0/0.70 , 1.0/0.67 , 1.0      , 1.0      ,$
	            1.0      , 1.0/0.67 , 1.00     , 1.0      ,$
		    1.0      , 1.0      ]
	chan14  = [ 0        , 5        , 19	   , 0	      , 0        ,$
	            22       , 19       , 22       , 0        ,$
		    0        , 0        ]
	diag14  = [ 'EVL'    ,'FVL'     , 'EVL'	   , 'EVL'    , 'EVL'    ,$
	            'FVL'    ,'EVL'     , 'FVL'    , 'EVL'    ,$
		    'EVL'    ,'EVL'     ] 

	rawelm = fltarr(n_elements(shot_db))+0.0 & rawelm[where(shot_db eq 30306 or shot_db eq 30307 or shot_db eq 35356 or shot_db eq 35158)]=1 
	elmcond= fltarr(n_elements(shot_db))+4.5 & elmcond[where(shot_db eq 30306 or shot_db eq 30307 or shot_db eq 35158)]=2 & elmcond[where(shot_db eq 35356)]=3

	setgraphics,colors=colors
	for i=0,n_elements(shot_db)-1 do begin
		if keyword_set(pickshot)then begin
			if shot_db[i] ne pickshot then goto,skip
		endif
		t0 = twindow[0,i]
		t1 = twindow[1,i]
		data12 = get_nii(shot_db[i],diag=diag12[i],los=los12[i],use_evl=useevl[i],channel=chan12[i],xr=[twindow[0,i],twindow[1,i]],rawelm=rawelm[i],elmcond=elmcond[i],interelm=interelm[i])
		data14 = get_nii(shot_db[i],diag=diag14[i],los=los14[i],use_evl=useevl[i],channel=chan14[i],xr=[twindow[0,i],twindow[1,i]],rawelm=rawelm[i],elmcond=elmcond[i],interelm=interelm[i])
		plot,data12.tdiv,smooth(data12.nii3995_ie/1e19,10,/edge_truncate),back=colors.white,col=colors.black,$
			title=string(shot_db[i],format='("AUG #",i5)'),xtitle='Tdiv [eV]',ytitle='N II [10!u19!n ph/s/m!u2!n/sr]'
		oplot,data14.tdiv,smooth(data14.nii3995_ie/1e19,10,/edge_truncate),col=colors.red
		cursor,x,y,/up
		skip:
	endfor

End

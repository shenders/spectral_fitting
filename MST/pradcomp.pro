Pro pradcomp
	tdiv=8
	if tdiv eq 8 then begin
		; Set shots for Te=8 eV
		shotdb= [35157, 35167, 35358]
		shots = [35157 , 35157 , 35167 , 35167 , 35167 , 35167 , 35358 , 35358 , 35358 , 35381 , 35381]
		nshot = 4
		; Set initial time for 200 ms window
		times = [2.8   , 4.8   , 2.8   , 3.6   , 4.8   , 5.5   , 3.3   , 4.8   , 5.3   , 2.3   , 2.9  ]
		twin  = [0.2   , 0.17  , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2  ]
	endif
	if tdiv eq 0 then begin
		; Set shots for Te=0 eV
		shotdb= [35158, 35165, 35381]
		shots = [35158 , 35158 , 35158 , 35158 , 35165 , 35165 , 35381 , 35381         ]
		nshot = 3
		; Set initial time for 200 ms window
		times = [2.8   , 3.8   , 4.8   , 5.8   , 2.8   , 3.8   , 2.8   , 4.3           ]
		twin  = [0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2           ]
	endif
	if tdiv eq 5 then begin
		; Set shots for Te=5 eV
		shotdb= [35157, 35167, 35358]
		shots = [33257 , 33257 , 33257 , 33257     ]
		nshot = 1
		; Set initial time for 200 ms window
		times = [2.5   , 3.4   , 4.2   , 5.2       ]
		twin  = [0.1   , 0.1   , 0.1   , 0.1       ]
	endif
	adas_colors, colors=colors
	for i=0,n_elements(shotdb)-1 do begin
		read_signal_mrm,0L,shotdb[i],'BPD','Prad',tline,prad,2
		id = where(shots eq shotdb[i])
		window,0
		!p.multi=0
		plot,tline,prad,xr=[1,6],yr=[0,1e7]
		for j=0,n_elements(id)-1 do begin
			x=bolom(shotdb[i],[times[id[j]],times[id[j]]+twin[id[j]]],/noplot)
			oplot,[times[id[j]],times[id[j]]+twin[id[j]]],[x,x]*1e6,col=colors.red
		endfor
		stop
	endfor
End

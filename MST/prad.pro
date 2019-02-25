Pro prad,debug=debug,tdiv=tdiv,traces=traces
	
;	tdiv = 8

	if ~keyword_set(tdiv)then tdiv=0

	if tdiv eq 8 then begin
		; Set shots for Te=8 eV
		shots = [35157 , 35157 , 35167 , 35167 , 35167 , 35167 , 35358 , 35358 , 35358 , 35381 , 35381]
		nshot = 4
		; Set initial time for 200 ms window
		times = [2.8   , 4.8   , 2.8   , 3.6   , 4.8   , 5.5   , 3.3   , 4.8   , 5.3   , 2.3   , 2.9  ]
		twin  = [0.2   , 0.17  , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2  ]
	endif
	if tdiv eq 0 then begin
		; Set shots for Te=0 eV
		shots = [35158 , 35158 , 35158 , 35158 , 35165 , 35165 , 35381 , 35381         ]
		nshot = 3
		; Set initial time for 200 ms window
		times = [2.8   , 3.8   , 4.8   , 5.8   , 2.8   , 3.8   , 2.8   , 4.3           ]
		twin  = [0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2   , 0.2           ]
	endif
	if tdiv eq 5 then begin
		; Set shots for Te=5 eV
		shots = [33257 , 33257 , 33257 , 33257     ]
		nshot = 1
		; Set initial time for 200 ms window
		times = [2.5   , 3.4   , 4.2   , 5.2       ]
		twin  = [0.1   , 0.1   , 0.1   , 0.1       ]
	endif
	

	; Initialise arrays
	ar_arr = fltarr(n_elements(shots))
	ar_err = fltarr(n_elements(shots))
	pr_arr = fltarr(n_elements(shots))
	pr_err = fltarr(n_elements(shots))
	pt_arr = fltarr(n_elements(shots))
	pt_err = fltarr(n_elements(shots))
	pd_arr = fltarr(n_elements(shots))
	pd_err = fltarr(n_elements(shots))
	
	if keyword_set(traces)then begin
		; Loop around database
		ref = 0
		setgraphics,xs=1400,ys=800,ncol=3,nrow=nshot,psplot=psplot,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list
		xrange=[2.0,6.0]
		for i=0,n_elements(shots)-1 do begin
			; Initialise shots
			if ref eq shots[i] then begin
				swch = 0 
			endif else begin
				ref = shots[i]
				swch = 1
			endelse
			if swch eq 1 then begin
				if shots[i] eq 35167 then begin		
					read_signal_mrm,0L,35158,'UVS','CFA03A',tline,ar,2
					ar = ar/1.1e21 * 0.82e21
				endif else begin			
					read_signal_mrm,0L,shots[i],'UVS','CFA03A',tline,ar,2
				endelse
				ymax = 3
				plot, tline, ar/1e21, xr=xrange, xs=1, yr=[0,ymax],/nodata, col=colors.black, back=colors.white,$
				title=string(shots[i],format='("Ar seeding rate: ",I5)'),xtitle='Time [s]', ytitle='[10!u21!n #/s]'
				id = where(shots eq shots[i])
				for j=0,n_elements(id)-1 do begin
					xstart = [times[id[j]],times[id[j]]]
					obandx,xstart,[0,ymax],xstart + twin[id[j]],/norm,col=colors.red
				endfor
				oplot,tline, ar/1e21, col=colors.black
			endif
		endfor
		ref = 0
		for i=0,n_elements(shots)-1 do begin
			; Initialise shots
			if ref eq shots[i] then begin
				swch = 0 
			endif else begin
				ref = shots[i]
				swch = 1
			endelse
			if swch eq 1 then begin
		
				read_signal_mrm,0L,shots[i],'BPD','Prad',tline,prad,2
				read_signal_mrm,0L,shots[i],'BPD','Pradtot',tline,pradtot,2
				
				ymax = 15
				plot, tline, prad/1e6, xr=xrange, xs=1, yr=[0,ymax],/nodata, col=colors.black, back=colors.white,$
				title=string(shots[i],format='("Prad[black]; Pradtot[blue]: ",I5)'),xtitle='Time [s]', ytitle='[MW]'

				id = where(shots eq shots[i])
				
				for j=0,n_elements(id)-1 do begin
					xstart = [times[id[j]],times[id[j]]]
					obandx,xstart,[0,ymax],xstart + twin[id[j]],/norm,col=colors.red
				endfor
				
				oplot,tline, prad/1e6, col=colors.black
				oplot,tline, pradtot/1e6, col=colors.blue
			endif
		endfor
		ref = 0
		for i=0,n_elements(shots)-1 do begin
			; Initialise shots
			if ref eq shots[i] then begin
				swch = 0 
			endif else begin
				ref = shots[i]
				swch = 1
			endelse
			if swch eq 1 then begin
		
				read_signal_mrm,0L,shots[i],'MAC','Tdiv',tline,tdv,2
				ymax = 20
				plot, tline, tdv, xr=xrange, xs=1, yr=[-10,ymax],/nodata, col=colors.black, back=colors.white,$
				title=string(shots[i],format='("Tdiv: ",I5)'),xtitle='Time [s]', ytitle='[eV]'
				id = where(shots eq shots[i])
				
				for j=0,n_elements(id)-1 do begin
					xstart = [times[id[j]],times[id[j]]]
					obandx,xstart,[-10,ymax],xstart + twin[id[j]],/norm,col=colors.red
				endfor
				
				oplot,tline, tdv, col=colors.black
				
			endif
			
			
		endfor
	
	endif
	; Loop around database
	ref = 0
	for i=0,n_elements(shots)-1 do begin

		if keyword_set(debug)then begin
			print,'--------------------'
			print,'Shot: ',shots[i]
		endif 
		
		; Initialise shots
		if ref eq shots[i] then begin
			swch = 0 
		endif else begin
			ref = shots[i]
			swch = 1
		endelse

		; Retrieve Ar seeding rate
		if shots[i] eq 35167 then begin		
			read_signal_mrm,0L,35158,'UVS','CFA03A',tline,ar,2
			ar = ar/1.1e21 * 0.82e21
		endif else begin
			if shots[i] eq 35381 and swch eq 1 then begin	
				read_signal_mrm,0L,35158,'UVS','CFA03A',tline,ar,2
			endif else begin
				read_signal_mrm,0L,shots[i],'UVS','CFA03A',tline,ar,2
			end	
		end
		bindata,tline,ar/1e21,[times[i],times[i]+ twin[i]],avr_ar,err_ar

		if keyword_set(debug)then begin
			print,'Ar seeding rate: ',avr_ar
		endif 

		ar_arr[i] = avr_ar
		ar_err[i] = err_ar

		; Retrieve Prad 
		
		if shots[i] eq 35381 and swch eq 1 then begin
			read_signal_mrm,0L,35158,'BPD','Prad',tline,prad,2
		endif else begin
			read_signal_mrm,0L,shots[i],'BPD','Prad',tline,prad,2
		end
		bindata,tline,prad,[times[i],times[i]+ twin[i]],avr_pr,err_pr
		if swch eq 1 then begin
			ref_pr_avr = avr_pr
			ref_pr_err = err_pr
		endif

		if keyword_set(debug)then begin
			print,'Prad: ',avr_pr
		endif 

		pr_arr[i] = avr_pr/ref_pr_avr
		pr_err[i] = pr_arr[i] * sqrt((err_pr/avr_pr)^2 + (ref_pr_err/ref_pr_avr)^2)
		pr_arr[i] = (pr_arr[i] - 1.0) * 100
		pr_err[i] = (pr_err[i]) * 100
		; Retrieve Pradtot

		if shots[i] eq 35381 and swch eq 1 then begin
			read_signal_mrm,0L,35158,'BPD','Pradtot',tline,pradtot,2
		endif else begin
			read_signal_mrm,0L,shots[i],'BPD','Pradtot',tline,pradtot,2
		end

		bindata,tline,pradtot,[times[i],times[i]+ twin[i]],avr_pt,err_pt
		if swch eq 1 then begin
			ref_pt_avr = avr_pt
			ref_pt_err = err_pt
		endif

		if keyword_set(debug)then begin
			print,'Pradtot: ',avr_pt
		endif 

		pt_arr[i] = avr_pt/ref_pt_avr
		pt_err[i] = pt_arr[i] * sqrt((err_pt/avr_pt)^2 + (ref_pt_err/ref_pt_avr)^2)
		pt_arr[i] = (pt_arr[i] - 1.0) * 100
		pt_err[i] = (pt_err[i]) * 100

		pd_arr[i] = (avr_pt-avr_pr)/(ref_pt_avr - ref_pr_avr)
		pd_arr[i] = (pd_arr[i] - 1.0) * 100
	endfor
	
	if keyword_set(debug)then begin
		print,'-----------------------------------'
	endif 
	
	
	
	
	setgraphics,xs=xs,ys=ys,ncol=ncol,nrow=nrow,psplot=psplot,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list
	user_psym,1,/fill
	
	plot  , ar_arr, pr_arr, psym=5, yr=[0,100], xr=[0,2], col=colors.black, back=colors.white, /nodata,$
		xtitle = 'Ar [10!u21!n #/s]', ytitle= '% Radiation Increase', title=string(tdiv,format='("Tdiv = ",D3.1," eV")')
	xyouts, [0.3,0.3],[60,55],['BPD: Prad','BPD: Pradtot'],col=[colors.blue,colors.red]
	oplot , ar_arr, pr_arr, psym=8, col=colors.blue, symsize=1.5
	errors, ar_arr, pr_arr, ystd = pr_err, col=colors.blue, xmax=2, xmin=0, ymax=100.0, ymin=0
	oplot , ar_arr, pt_arr, psym=8, col=colors.red , symsize=1.5
	errors, ar_arr, pt_arr, ystd = pt_err, col=colors.red, xmax=2, xmin=0, ymax=100.0, ymin=0  

	;oplot , ar_arr, pd_arr, psym=8, col=colors.green , symsize=1.5
	
	if tdiv eq 8 then begin
		id = where(ar_arr le 0.9)
		pr_fit = linfit(ar_arr[id],pr_arr[id],measure_errors = pr_err[id])
		xarr = findgen(100)*0.8/99.0
		yarr = interpol(pr_fit[0] + pr_fit[1] * ar_arr, ar_arr, xarr)
		oplot,xarr,yarr,col=colors.blue
		id = where(ar_arr le 1.0)
		xarr = findgen(100)*0.95/99.0
		pt_fit = linfit(ar_arr[id],pt_arr[id],measure_errors = pt_err[id])
		yarr = interpol(pt_fit[0] + pt_fit[1] * ar_arr, ar_arr, xarr)
		oplot,xarr,yarr,col=colors.red

		id = where(ar_arr ge 0.8)
		pr_fit = linfit(ar_arr[id],pr_arr[id],measure_errors = pr_err[id])
		xarr = findgen(100)*1.25/99.0 + 0.85
		yarr = interpol(pr_fit[0] + pr_fit[1] * ar_arr, ar_arr, xarr)
		oplot,xarr,yarr,col=colors.blue,linest=5
		id = where(ar_arr ge 0.9)
		xarr = findgen(100)*1.25/99.0 + 1.0
		pt_fit = linfit(ar_arr[id],pt_arr[id],measure_errors = pt_err[id])
		yarr = interpol(pt_fit[0] + pt_fit[1] * ar_arr, ar_arr, xarr)
		oplot,xarr,yarr,col=colors.red,linest=5
	
		oplot,[0.9,0.9],[0,100],col=colors.black, linest=5
	endif
	if tdiv lt 7 then begin
		id = where(ar_arr le 1.2)
		pr_fit = linfit(ar_arr[id],pr_arr[id],measure_errors = pr_err[id])
		xarr = findgen(100)*1.1/99.0
		yarr = interpol(pr_fit[0] + pr_fit[1] * ar_arr, ar_arr, xarr)
		oplot,xarr,yarr,col=colors.blue
		pt_fit = linfit(ar_arr[id],pt_arr[id],measure_errors = pt_err[id])
		yarr = interpol(pt_fit[0] + pt_fit[1] * ar_arr, ar_arr, xarr)
		oplot,xarr,yarr,col=colors.red

		id = where(ar_arr ge 1.0)
		pr_fit = linfit(ar_arr[id],pr_arr[id],measure_errors = pr_err[id])
		xarr = findgen(100)*0.9/99.0 + 1.1
		yarr = interpol(pr_fit[0] + pr_fit[1] * ar_arr, ar_arr, xarr)
		oplot,xarr,yarr,col=colors.blue,linest=5
		pt_fit = linfit(ar_arr[id],pt_arr[id],measure_errors = pt_err[id])
		yarr = interpol(pt_fit[0] + pt_fit[1] * ar_arr, ar_arr, xarr)
		oplot,xarr,yarr,col=colors.red,linest=5
	
		oplot,[1.2,1.2],[0,100],col=colors.black, linest=5
	endif
stop
End

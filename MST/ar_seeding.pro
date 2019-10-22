Function fittype,type=type
	if ~keyword_set(type)then type='quadratic'
	if type eq 'linear' then a=[1,1] else a=[1,1,1]
	return,a
end
Pro ar_seeding,abel=abel,traces=traces,shots=shots,$
               times=times,twin=twin,psplot=psplot,$
	       power=power,diag=diag,shot_ignore=shot_ignore
	
	if ~keyword_set(tdiv)then tdiv=0
	if ~keyword_set(power)then power='med'
	if power eq 'med' then begin
		shots = [35157 , 35157 , 35158 , 35158 , 35158 , 35158 , 35167 , 35167 , 35167 , 35167 , 35358 , 35358 , 35358 , 35381 , 35381 , 35381]
		; Set initial time for 200 ms window
		times = [2.5   , 4.5   , 2.5   , 3.5   , 4.1   , 5.5   , 2.5   , 3.55  , 4.5   , 5.5   , 3.0   , 4.0   , 5.0   , 2.0   , 3.0   , 4.0]
		ref   = [1     , 0     , 1     , 0     , 0     , 0     , 1     , 0     , 0     , 0     , 1     , 0     , 0     , 1.0   , 0     , 0  ]
		twin  = fltarr(n_elements(times))+0.5
		tdivs = [8     , 8     , 0     , 0     , 0     , 0     , 8     , 8     , 8     , 8     , 8     , 8     , 8     , 8     , 8     , 0]
		title = 'Medium power'
		file  = 'figures/ar_frac_medpow.ps'
		tdiv  = -5
	endif
	if power eq 'high' then begin
		shots = [33257 , 33257 , 33257 , 33257 , 33258 , 33258 , 33258 , 33258]
		times = [2.4   , 3.0   , 4.0   , 5.0   , 2.4   , 3.0   , 3.85  , 5.0  ]
		ref   = [1     , 0     , 0     , 0     , 1     , 0     , 0     , 0    ]
		twin  = [0.3   , 0.3   , 0.3   , 0.3   , 0.3   , 0.3   , 0.3   , 0.3  ]
		tdivs = fltarr(n_elements(twin))+5
		file  = 'figures/ar_frac_highpow.ps'
		tdiv  = -5
		title = 'High power'
	endif
	
	if ~keyword_set(diag)then diag = 'BPD'
	if diag eq 'BPD' then begin
		shot_ignore = [-1]
		trace1 = 'Prad'
		trace2 = 'Pradtot'
		exp    = 'AUGD'
	endif else begin
		diag   = 'BPT'
		shot_ignore = [-1] 
		trace1 = 'Pr_main'
		trace2 = 'Pr_tot'
		exp    = 'davidp' 
	end
	arr_keep=-1
	for i=0,n_elements(shots)-1 do begin
		id_shot_check = where(shot_ignore eq shots[i])
		if id_shot_check[0] eq -1 then arr_keep=[arr_keep,i] 
	endfor
	arr_keep=arr_keep[1:*]
	shots=shots[arr_keep]
	times=times[arr_keep]
	ref  =ref[arr_keep]
	twin =twin[arr_keep]
	tdivs=tdivs[arr_keep]
	
	nshot = 0 & for i=0,n_elements(shots)-1 do begin & if ref[i] eq 1 then begin nshot = nshot + 1 & endif & endfor

	if keyword_set(traces)then begin
		; Loop around database

		setgraphics,xs=1400,ys=800,ncol=4,nrow=nshot,psplot=psplot,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list
		xrange=[1.0,7.0]
		for i=0,n_elements(shots)-1 do begin
			; Initialise shots
			if ref[i] eq 1 then begin
				read_signal_mrm,0L,shots[i],'TOT','P_TOT',tline,p_tot,2
				ymax = 25
				plot, tline, p_tot/1e6, xr=xrange, xs=1, yr=[0,ymax],/nodata, col=colors.black, back=colors.white,$
				title=string(shots[i],format='("Total power: ",I5)'),xtitle='Time [s]', ytitle='[MW]'
				id = where(shots eq shots[i] and ref ne 2)
				for j=0,n_elements(id)-1 do begin
					xstart = [times[id[j]],times[id[j]]]
					obandx,xstart,[0,ymax],xstart + twin[id[j]],/norm,col=colors.red
				endfor
				oplot,tline, p_tot/1e6, col=colors.black
			endif
		endfor

		for i=0,n_elements(shots)-1 do begin
			; Initialise shots
			if ref[i] eq 1 then begin
				if shots[i] eq 35167 then begin		
					read_signal_mrm,0L,35158,'UVS','CFA03A',tline,ar,2
					ar = ar/1.1e21 * 0.82e21
				endif else begin			
					read_signal_mrm,0L,shots[i],'UVS','CFA03A',tline,ar,2
				endelse
				ymax = 3
				plot, tline, ar/1e21, xr=xrange, xs=1, yr=[0,ymax],/nodata, col=colors.black, back=colors.white,$
				title=string(shots[i],format='("Ar seeding rate: ",I5)'),xtitle='Time [s]', ytitle='[10!u21!n #/s]'
				id = where(shots eq shots[i] and ref ne 2)
				for j=0,n_elements(id)-1 do begin
					xstart = [times[id[j]],times[id[j]]]
					obandx,xstart,[0,ymax],xstart + twin[id[j]],/norm,col=colors.red
				endfor
				oplot,tline, ar/1e21, col=colors.black
			endif
		endfor

		for i=0,n_elements(shots)-1 do begin
			; Initialise shots
			if ref[i] eq 1 then begin
				if keyword_set(abel)then begin
					runabel,shots=shots[i],/noplot,/load,time=tline2,prad=prad
					read_signal_mrm,0L,shots[i],diag,trace2,tline,pradtot,2,exp=exp
					prad = prad * 1e6
				endif else begin
					read_signal_mrm,0L,shots[i],diag,trace1,tline2,prad,2,exp=exp
					read_signal_mrm,0L,shots[i],diag,trace2,tline,pradtot,2,exp=exp
				end		
				if power eq 'high' then ymax = 20 else ymax = 15
				plot, tline2, prad/1e6, xr=xrange, xs=1, yr=[0,ymax],/nodata, col=colors.black, back=colors.white,$
				title=string(shots[i],format='("Prad[blue]; Pradtot[black]: ",I5)'),xtitle='Time [s]', ytitle='[MW]'

				id = where(shots eq shots[i] and ref ne 2)
				
				for j=0,n_elements(id)-1 do begin
					xstart = [times[id[j]],times[id[j]]]
					obandx,xstart,[0,ymax],xstart + twin[id[j]],/norm,col=colors.red
				endfor
				
				oplot,tline2, prad/1e6, col=colors.blue
				oplot,tline, pradtot/1e6, col=colors.black
			endif
		endfor

		for i=0,n_elements(shots)-1 do begin
			; Initialise shots

			if ref[i] eq 1 then begin
				read_signal_mrm,0L,shots[i],'MAC','Tdiv',tline,tdv,2
				read_signal_mrm,0L,shots[i],'UVS','N_tot',tline_n2,rates_n2,2
				
				ymax = 20
				plot, tline, tdv, xr=xrange, xs=1, yr=[-10,ymax],/nodata, col=colors.black, back=colors.white,$
				title=string(shots[i],format='("Tdiv [black]| N_tot [blue] | cN [green]: ",I5)'),xtitle='Time [s]', ytitle='[eV,E21]'
				id = where(shots eq shots[i] and ref ne 2)
				for j=0,n_elements(id)-1 do begin
					xstart = [times[id[j]],times[id[j]]]
					obandx,xstart,[-10,ymax],xstart + twin[id[j]],/norm,col=colors.red
				endfor
				
				oplot,tline, tdv, col=colors.black
				oplot,tline_n2,rates_n2/1e21,col=colors.blue
				if power ne 'high'then begin
					runcn,shots=shots[i],/load,time=cn_time,cn_mean=cn_mean
					oplot,cn_time,cn_mean*100,col=colors.green
				endif
			endif
			
			
		endfor
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
	n2_arr = fltarr(n_elements(shots))
	n2_err = fltarr(n_elements(shots))
	cn_arr = fltarr(n_elements(shots))
	cn_err = fltarr(n_elements(shots))
	td_arr = fltarr(n_elements(shots))
	td_err = fltarr(n_elements(shots))
	

	; Loop around database

	for i=0,n_elements(shots)-1 do begin

		if keyword_set(debug)then begin
			print,'--------------------'
			print,'Shot: ',shots[i]
		endif 
		
		; Retrieve Ar seeding rate
		if shots[i] eq 35167 then begin		
			read_signal_mrm,0L,35158,'UVS','CFA03A',tline,ar,2
			ar = ar/1.1e21 * 0.82e21
		endif else begin
			if shots[i] eq 35381 and ref[i] gt 1 then begin	
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
		if shots[i] eq 35381 and ref[i] gt 1 then begin
			if keyword_set(abel)then begin
				runabel,shots=35158,/load,/noplot,time=tline,prad=prad
				prad = prad * 1e6
			endif else begin
				read_signal_mrm,0L,35158,diag,trace1,tline,prad,2,exp=exp
			end
		endif else begin
			if keyword_set(abel)then begin
				runabel,shots=shots[i],/load,/noplot,time=tline,prad=prad
				prad = prad * 1e6
			endif else begin
				read_signal_mrm,0L,shots[i],diag,trace1,tline,prad,2,exp=exp
			end
		end
		bindata,tline,prad/1e6,[times[i],times[i]+ twin[i]],avr_pr,err_pr
		if ref[i] gt 0 then begin
			ref_pr_avr = avr_pr
			ref_pr_err = err_pr
		endif

		if keyword_set(debug)then begin
			print,'Prad: ',avr_pr
		endif 
		if power eq 'high' then norm = 7 else norm = 2.0
		pr_arr[i] = avr_pr* norm/ref_pr_avr
		;pr_err[i] = pr_arr[i] * sqrt((err_pr/avr_pr)^2 + (ref_pr_err/ref_pr_avr)^2)
		pr_err[i] = err_pr
		;pr_arr[i] = (pr_arr[i] - 1.0) * 100
		;pr_err[i] = (pr_err[i]) * 100
		; Retrieve Pradtot

		if shots[i] eq 35381 and ref[i] gt 1 then begin
			read_signal_mrm,0L,35158,diag,trace2,tline,pradtot,2,exp=exp
		endif else begin
			read_signal_mrm,0L,shots[i],diag,trace2,tline,pradtot,2,exp=exp
		end

		bindata,tline,pradtot/1e6,[times[i],times[i]+ twin[i]],avr_pt,err_pt
		if ref[i] gt 0 then begin
			ref_pt_avr = avr_pt
			ref_pt_err = err_pt
		endif

		if keyword_set(debug)then begin
			print,'Pradtot: ',avr_pt
		endif 

		if power eq 'high' then norm = 15 else norm = 6.0
		pt_arr[i] = avr_pt*norm/ref_pt_avr
		;pt_err[i] = pt_arr[i] * sqrt((err_pt/avr_pt)^2 + (ref_pt_err/ref_pt_avr)^2)
		pt_err[i] = err_pt
		;pt_arr[i] = (pt_arr[i] - 1.0) * 100
		;pt_err[i] = (pt_err[i]) * 100

		pd_arr[i] = (avr_pt-avr_pr)/(ref_pt_avr - ref_pr_avr)
		pd_arr[i] = (pd_arr[i] - 1.0) * 100
		; Retrieve N2 
		if shots[i] eq 35381 and ref[i] gt 1 then begin
			read_signal_mrm,0L,35158,'UVS','N_tot',N2_time,N2_rate,2
			read_signal_mrm,0L,35158,'UVS','D_tot',D2_time,D2_rate,2
			n2_rate = n2_rate/1e21
			d2_rate = d2_rate/1e21
		endif else begin
			read_signal_mrm,0L,shots[i],'UVS','N_tot',N2_time,N2_rate,2
			read_signal_mrm,0L,shots[i],'UVS','D_tot',D2_time,D2_rate,2
			n2_rate = n2_rate/1e21
			d2_rate = d2_rate/1e21
		end
		bindata,N2_time,(N2_rate/7)/(D2_rate+N2_rate/7),[times[i],times[i]+ twin[i]],avr_n2,err_n2

		if ref[i] gt 0 then begin
			ref_n2_avr = avr_n2
			ref_n2_err = err_n2
		endif

		if keyword_set(debug)then begin
			print,'N2: ',avr_n2
		endif 

		n2_arr[i] = avr_n2
		n2_err[i] = err_n2
		n2_arr[i] = n2_arr[i] * 100
		n2_err[i] = n2_err[i] * 100

		; Retrieve cN
		if power ne 'high' then begin 
		if shots[i] eq 35381 and ref[i] gt 1 then begin
			runcn,shots=35158,/load,cn_mean=cn_mean
		endif else begin
			runcn,shots=shots[i],/load,time=cn_time,cn_mean=cn_mean
		end
		bindata,cn_time,cn_mean,[times[i],times[i]+ twin[i]],avr_cn,err_cn
		if ref[i] gt 0 then begin
			ref_cn_avr = avr_cn
			ref_cn_err = err_cn
		endif

		if keyword_set(debug)then begin
			print,'cN: ',avr_cn
		endif 

		;cn_arr[i] = avr_cn/ref_cn_avr
		cn_arr[i] = avr_cn
		cn_err[i] = err_cn
		cn_arr[i] = cn_arr[i] * 100
		cn_err[i] = cn_err[i] * 100

		endif else begin
			cn_arr[i] = -1.0
			cn_err[i] = -0.1
		end
		; Retrieve Tdiv 
		if shots[i] eq 35381 and ref[i] gt 1 then begin
			read_signal_mrm,0L,35158,'MAC','Tdiv',Tdiv_time,Tdiv_rate,2
		endif else begin
			read_signal_mrm,0L,shots[i],'MAC','Tdiv',Tdiv_time,Tdiv_rate,2
		end
		bindata,Tdiv_time,Tdiv_rate,[times[i],times[i]+ twin[i]],avr_td,err_td

		if keyword_set(debug)then begin
			print,'Tdiv: ',avr_td
		endif 
		td_arr[i] = avr_td
		td_err[i] = err_td
	endfor
	
	if keyword_set(debug)then begin
		print,'-----------------------------------'
	endif 
	
	
	
	
	setgraphics,xs=800,ys=1000,ncol=3,nrow=1,psplot=psplot,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list,file=file
	user_psym,1,/fill
	xr=[-0.05,2.2]
	idsort = sort(ar_arr)
	if power eq 'high'then yr=[0,25] else yr=[0,10]
	plot  , ar_arr, pr_arr, psym=5, yr=yr, xr=xr,xs=1, col=colors.black, back=colors.white, /nodata,$
		ytitle= 'Radiated power [MW]', title=title, xtitle = 'Ar [10!u21!n #/s]'
	xyouts, [1.7,1.7],[9,7],['P!lcore!n','P!ltotal!n'],col=[colors.blue,colors.red]
	oplot , ar_arr, pr_arr, psym=8, col=colors.blue, symsize=1.5
	errors, ar_arr, pr_arr, ystd = pr_err, col=colors.blue, xmax=xr[1], xmin=xr[0], ymax=120.0, ymin=0
	oplot , ar_arr, pt_arr, psym=8, col=colors.red , symsize=1.5
	errors, ar_arr, pt_arr, ystd = pt_err, col=colors.red, xmax=xr[1], xmin=xr[0], ymax=120.0, ymin=0  
	res=svdfit(ar_arr[idsort],pt_arr[idsort],yfit=yfit,measure_errors=pt_err[idsort],a=fittype(type='quad'))
	xaxis = findgen(100)*20/99.0
	yaxis = fltarr(100)
	for i=0,n_elements(res)-1 do yaxis = yaxis + res[i] * xaxis^(i)
	;oplot,xaxis,yaxis,col=colors.black
	
	res=svdfit(ar_arr[idsort],pr_arr[idsort],yfit=yfit,measure_errors=pr_err[idsort],a=fittype(type='quad'))	
	xaxis1 = findgen(100)*20/99.0
	yaxis1 = fltarr(100)
	for i=0,n_elements(res)-1 do yaxis1 = yaxis1 + res[i] * xaxis1^(i)
	;oplot,xaxis1,yaxis1,col=colors.black
	if power eq 'med'then begin
		oplot,[-1,10],[2,2],linest=5,col=colors.blue
		oplot,[-1,10],[6,6],linest=5,col=colors.red
	endif else begin
		oplot,[-1,10],[7,7],linest=5,col=colors.blue
		oplot,[-1,10],[15,15],linest=5,col=colors.red
	
	end
	plot  , ar_arr, n2_arr, psym=5, yr=[0,15], xr=xr,xs=1, col=colors.black, back=colors.white, /nodata,$
		ytitle= 'Nitrogen concentration [%]'
	oplot , ar_arr, n2_arr, psym=8, col=colors.black, symsize=2
	errors, ar_arr, n2_arr, ystd = n2_err, col=colors.black, xmax=xr[1], xmin=xr[0], ymax=120.0, ymin=0
	xyouts, [1.5,1.5],[12,10],['N!l2!n valve flux','c!lN!n div. spec.'],col=[colors.black,colors.aqua]
	user_psym,2,/fill
	if power eq 'high' then begin
		res=svdfit(ar_arr[idsort],n2_arr[idsort],yfit=yfit,measure_errors=n2_err[idsort],a=fittype(type='quad'))
	endif else begin
		newxarr = [ar_arr,ar_arr]
		newyarr = [cn_arr,n2_arr]
		newearr = [cn_err,n2_err]
		idsort  = sort(newxarr)
		res=svdfit(newxarr[idsort],newyarr[idsort],yfit=yfit,measure_errors=newearr[idsort],a=fittype(type='quad'))
		oplot , ar_arr, cn_arr, psym=8, col=colors.aqua, symsize=2
		errors, ar_arr, cn_arr, ystd = cn_err, col=colors.aqua, xmax=xr[1], xmin=xr[0], ymax=120.0, ymin=0
	end
	xaxis = findgen(100)*20/99.0
	yaxis = fltarr(100)
	for i=0,n_elements(res)-1 do yaxis = yaxis + res[i] * xaxis^(i)
	;oplot,xaxis,yaxis,col=colors.black

	plot  , ar_arr, td_arr, psym=5, yr=[-10,15], xr=xr,xs=1, col=colors.black, back=colors.white, /nodata,$
		xtitle = 'Ar [10!u21!n #/s]', ytitle= 'Tdiv [eV]'
	
	id = where(tdivs eq 8)
	if id[0] ne -1 then begin
		oplot , ar_arr[id], td_arr[id], psym=8, col=colors.black, symsize=1.5
		errors, ar_arr[id], td_arr[id], ystd = td_err, col=colors.black, xmax=xr[1], xmin=xr[0], ymax=15, ymin=-10
		oplot,xr,[8,8],col=colors.black,linest=5
	endif
	id = where(tdivs eq 5)
	if id[0] ne -1 then begin
		oplot , ar_arr[id], td_arr[id], psym=8, col=colors.black, symsize=1.5
		errors, ar_arr[id], td_arr[id], ystd = td_err, col=colors.black, xmax=xr[1], xmin=xr[0], ymax=15, ymin=-10
		oplot,xr,[5,5],col=colors.black,linest=5
	endif
	id = where(tdivs eq 0)
	if id[0] ne -1 then begin
		oplot , ar_arr[id], td_arr[id], psym=8, col=colors.black, symsize=1.5
		errors, ar_arr[id], td_arr[id], ystd = td_err, col=colors.black, xmax=xr[1], xmin=xr[0], ymax=15, ymin=-10
		oplot,xr,[0,0],col=colors.black,linest=5
	endif

	setgraphics,psplot=psplot,/close
stop
End

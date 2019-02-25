Pro regress_cn,abel=abel

	restore,'save/cn_database.sav'

	
	; loop over all values
	
	
	PsepAbel= [ 4.31     , 4.05     , 5.80     , 5.77     ,$
	            8.99     , 14.21    , 5.91     , 7.67     ,$
		    5.67     , 6.97     ] * 1e6
	
	idsort = sort(database.cn_mean)
	if keyword_set(abel)then psep   = PsepAbel[idsort] else psep   = database.psep[idsort]
	yvals  = database.cn_mean[idsort]
	cn_err = database.cn_err[idsort]
	bp     = database.bp[idsort]
	k      = database.k[idsort]
	ngw    = database.ngw[idsort]
	num    = 20
	psep_vals = adas_vector(high=2.0,low=0.1,num=num,/linear)
	bp_vals   = adas_vector(high=-0.5,low=-3.0,num=num,/linear)	
	ngw_vals  = adas_vector(high=-0.5,low=-3.0,num=num,/linear)
	k_vals    = adas_vector(high=-1.0,low=-2.0,num=num,/linear)
	cn_vals   = adas_vector(high=3.0,low=4.0,num=num,/linear)
	
	r2    = fltarr(num,num,num,num,num)
	r2max = 0
	for i=0,n_elements(psep_vals)-1 do begin
		for j=0,n_elements(bp_vals)-1 do begin
			for l=0,n_elements(ngw_vals)-1 do begin
				for m=0,n_elements(k_vals)-1 do begin
					for n=0,n_elements(cn_vals)-1 do begin
						ybar   = mean(yvals)
						gldstn = (10.7^psep_vals[i]) * 0.34^bp_vals[j] * 0.496^ngw_vals[l] * (1.0+1.63^2)^k_vals[m]
						f_i    = (psep/1e6)^psep_vals[i] * bp^bp_vals[j] * ngw^ngw_vals[l] * (1.0+k^2)^k_vals[m] * (cn_vals[n]/gldstn)
						ss_tot = 0 & for ii = 0,n_elements(yvals)-1 do ss_tot = ss_tot + (yvals[ii]-ybar)^2
						ss_res = 0 & for ii = 0,n_elements(yvals)-1 do ss_res = ss_res + (yvals[ii]-f_i[ii])^2
						R2[i,j,l,m,n] = 1 - ss_res/ss_tot
						if r2[i,j,l,m,n] gt r2max then begin
							istore = i
							jstore = j
							lstore = l
							mstore = m
							nstore = n
							gold   = gldstn
							r2max  = r2[i,j,l,m,n]
						endif
					endfor
				endfor
			endfor
		endfor
	endfor


	xvals = (psep/1e6)^psep_vals[istore] * bp^bp_vals[jstore] * ngw^ngw_vals[lstore] * (1.0+k^2)^k_vals[mstore] * (cn_vals[nstore]/gold)
	
	user_psym,1,/fill
	setgraphics,xs=xs,ys=ys,ncol=1,nrow=1,psplot=psplot,filename='figures/scaling.ps',colors=colors,colpick=colpick,collabel=collabel,full_list=full_list	
	
	plot,xvals,yvals,psym=8,yr=[0,20],xr=[0,20],back=colors.white,col=colors.black
	err_plot,xvals,yvals,cn_err,col=colors.black

	fit   = linfit(xvals,yvals,measure_errors=cn_err)
	xarr  = findgen(1000)/999.0*100
	yarr  = fit[0] + xarr * fit[1]
	
	oplot,[0,100],[0,100],col=colors.black

	print,'Psep^[x]= ',psep_vals[istore]
	print,'Bp^[x]= ',bp_vals[jstore]
	print,'ngw^[x]= ',ngw_vals[lstore]
	print,'(1+k^2)^[x] = ',k_vals[mstore]
	print,'cn^[x] = ',cn_vals[nstore]
stop


End

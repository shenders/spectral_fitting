Function calc_profiles, data, dens_arr, te_arr, ratio1, ratio2, exp_ratio1, exp_ratio2 ,debug=debug

	num        = n_elements(dens_arr)
	idx        = fltarr(num)	
	min2       = fltarr(num)
	te         = fltarr(n_elements(exp_ratio1))
	dens       = fltarr(n_elements(exp_ratio1))
	sim_ratio1 = fltarr(n_elements(exp_ratio1))
	sim_ratio2 = fltarr(n_elements(exp_ratio1))
	
	for i=0,n_elements(exp_ratio1)-1 do begin
		for j=0,num-1 do begin
			min1   = abs((ratio1[*,j] - exp_ratio1[i])/exp_ratio1[i]) * $
	         		 abs((ratio2[*,j] - exp_ratio2[i])/exp_ratio2[i])
			idx[j] = where(min1 eq min(min1))
			min2[j]= min1[idx[j]]
		endfor
		idy     = where(min2 eq min(min2))
        	; check if reasonable solution found
		sim_ratio1[i] = ratio1[idx[idy[0]],idy[0]]
		sim_ratio2[i] = ratio2[idx[idy[0]],idy[0]]
		if abs(sim_ratio1[i]-exp_ratio1[i])/exp_ratio1[i] lt 0.05 and $
		   abs(sim_ratio2[i]-exp_ratio2[i])/exp_ratio2[i] lt 0.05 then begin
			dens[i] = dens_arr[idx[idy[0]]]
			te[i]   = te_arr[idy[0]]
		endif else begin
			dens[i]       = -1
			te[i]         = -1	
			sim_ratio1[i] = -1
			sim_ratio2[i] = -1
		end
	endfor
	id = where(dens ne -1,complement=iy)
	if id[0] ne -1 then begin
		te[iy]   = interpol(te[id],data.time[id],data.time[iy])
		dens[iy] = interpol(dens[id],data.time[id],data.time[iy])
		sim_ratio1[iy] = interpol(sim_ratio1[id],data.time[id],data.time[iy])
		sim_ratio2[iy] = interpol(sim_ratio2[id],data.time[id],data.time[iy])
	endif	

	return,{sim_ratio1: sim_ratio1,$
        	sim_ratio2: sim_ratio2,$
		dens      : dens      ,$
		te        : te         }
end

Function optimize_ratios,data,shot,los,transmission,sm=sm,debug=debug

	num      = 40
	te_arr   = adas_vector(high=6,low=2.,num=num)
	dens_arr = adas_vector(high=1e15,low=1e13,num=num)
	ratio1   = fltarr(num,num)
	ratio2   = fltarr(num,num)

	exp_ratio1 = smooth(data.nii4041,sm,/edge_truncate) / smooth(data.nii3995,sm,/edge_truncate)
	exp_ratio2 = smooth(data.nii4041,sm,/edge_truncate) / smooth(data.nii4026,sm,/edge_truncate)
	
	err_ratio1 = exp_ratio1 * sqrt((data.nii3995_err/data.nii3995)^2+(data.nii4041_err/data.nii4041)^2)
	err_ratio2 = exp_ratio2 * sqrt((data.nii4026_err/data.nii4026)^2+(data.nii4041_err/data.nii4041)^2)
	
	exp_ratio1_upper = exp_ratio1 + err_ratio1/2.0
	exp_ratio1_lower = exp_ratio1 - err_ratio1/2.0

	exp_ratio2_upper = exp_ratio2 + err_ratio2/2.0
	exp_ratio2_lower = exp_ratio2 - err_ratio2/2.0

	print,'Fetching atomic data...'
	for i=0,num-1 do begin
		atomdb,te_arr[i]+fltarr(num),dens_arr,$
		               tec3995=tec3995,$
		               tec4041=tec4041,$
		               tec4026=tec4026
		ratio1[*,i] = tec4041 / tec3995
		ratio2[*,i] = tec4041 / tec4026
	endfor

	; find best solution for 4 different cases

	if keyword_set(debug)then begin
		setgraphics,colors=colors,colpick=colpick,/full
		plot,[6,10],[0,0.5],/nodata,col=colors.black,back=colors.white	
		for i=0,num-1 do oplot,ratio2[*,i],ratio1[*,i],col=colpick[i]
		user_psym,5,/fill & oplot,exp_ratio2_upper,exp_ratio1_upper,psym=8,col=colors.black
		user_psym,5,/fill & oplot,exp_ratio2_upper,exp_ratio1_lower,psym=8,col=colors.black
		user_psym,5,/fill & oplot,exp_ratio2_lower,exp_ratio1_upper,psym=8,col=colors.black
		user_psym,5,/fill & oplot,exp_ratio2_lower,exp_ratio1_lower,psym=8,col=colors.black
	endif

        ; case 1: ratio1 + err with ratio2 + err
        res1 = calc_profiles(data, dens_arr, te_arr, ratio1, ratio2, exp_ratio1_upper, exp_ratio2_upper,debug=debug)

        ; case 2: ratio1 + err with ratio2 - err
        res2 = calc_profiles(data, dens_arr, te_arr, ratio1, ratio2, exp_ratio1_upper, exp_ratio2_lower,debug=debug)

        ; case 3: ratio1 - err with ratio2 + err
        res3 = calc_profiles(data, dens_arr, te_arr, ratio1, ratio2, exp_ratio1_lower, exp_ratio2_upper,debug=debug)

        ; case 3: ratio1 - err with ratio2 - err
        res4 = calc_profiles(data, dens_arr, te_arr, ratio1, ratio2, exp_ratio1_lower, exp_ratio2_lower,debug=debug)

	; find maximum and minimum values for ratios
	te_lower     = fltarr(n_elements(exp_ratio1))
	te_upper     = fltarr(n_elements(exp_ratio1))
	dens_lower   = fltarr(n_elements(exp_ratio1))
	dens_upper   = fltarr(n_elements(exp_ratio1))
	ratio1_lower = fltarr(n_elements(exp_ratio1))
	ratio2_lower = fltarr(n_elements(exp_ratio1))
	ratio1_upper = fltarr(n_elements(exp_ratio1))
	ratio2_upper = fltarr(n_elements(exp_ratio1))
	
	for i=0,n_elements(data.time)-1 do begin
		te_lower[i]     = res1.te[i]<res2.te[i]<res3.te[i]<res4.te[i]
		te_upper[i]     = res1.te[i]>res2.te[i]>res3.te[i]>res4.te[i]
		dens_lower[i]   = res1.dens[i]<res2.dens[i]<res3.dens[i]<res4.dens[i]
		dens_upper[i]   = res1.dens[i]>res2.dens[i]>res3.dens[i]>res4.dens[i]
		ratio1_lower[i] = res1.sim_ratio1[i]<res2.sim_ratio1[i]<res3.sim_ratio1[i]<res4.sim_ratio1[i]
		ratio1_upper[i] = res1.sim_ratio1[i]>res2.sim_ratio1[i]>res3.sim_ratio1[i]>res4.sim_ratio1[i]
		ratio2_lower[i] = res1.sim_ratio2[i]<res2.sim_ratio2[i]<res3.sim_ratio2[i]<res4.sim_ratio2[i]
		ratio2_upper[i] = res1.sim_ratio2[i]>res2.sim_ratio2[i]>res3.sim_ratio2[i]>res4.sim_ratio2[i]
	endfor

	; find line-integrated nii
        cn_1 = calc_cn(shot,los,transmission,te_lower,dens_lower,data,sm,err=err_1)
        cn_2 = calc_cn(shot,los,transmission,te_lower,dens_upper,data,sm,err=err_2)
        cn_3 = calc_cn(shot,los,transmission,te_upper,dens_lower,data,sm,err=err_3)
        cn_4 = calc_cn(shot,los,transmission,te_upper,dens_upper,data,sm,err=err_4)
		
	cn_lower = fltarr(n_elements(data.time))
	cn_upper = fltarr(n_elements(data.time))

	for i=0,n_elements(data.time)-1 do begin
		cn_lower[i] = (cn_1[i])<(cn_2[i])<(cn_3[i])<(cn_4[i])
		cn_upper[i] = (cn_1[i])>(cn_2[i])>(cn_3[i])>(cn_4[i])
	endfor
	return,{ratio1_upper:ratio1_upper,$
		ratio2_upper:ratio2_upper,$
		ratio1_lower:ratio1_lower,$
		ratio2_lower:ratio2_lower,$
		cn_upper:cn_upper,$
		cn_lower:cn_lower,$
		dens_upper:dens_upper,$
		dens_lower:dens_lower,$
		te_lower:te_lower,$
		te_upper:te_upper}
End

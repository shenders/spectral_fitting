Function estimate_ratios,data,shot,los,transmission,$
                         sm=sm,full=full,debug=debug,$
			 lowerte=lowerte,upperte=upperte

	
	num      = 40
	dens_arr = adas_vector(high=5e14,low=5e13,num=num)
	exp_ratio1 = smooth(data.nii4041,sm,/edge_truncate) / smooth(data.nii3995,sm,/edge_truncate)	
	err_ratio1 = exp_ratio1 * sqrt((data.nii3995_err/data.nii3995)^2+(data.nii4041_err/data.nii4041)^2)	
	exp_ratio1_upper = exp_ratio1 + err_ratio1/2.0
	exp_ratio1_lower = exp_ratio1 - err_ratio1/2.0

	te_arr   = fltarr(num)+lowerte
	atomdb,te_arr,dens_arr,tec3995=tec3995,$
		               tec4041=tec4041,$
		               tec4026=tec4026
	dens1 = interpol(dens_arr,tec4041/tec3995,exp_ratio1_lower)
	dens2 = interpol(dens_arr,tec4041/tec3995,exp_ratio1_upper)

	te_arr   = fltarr(num)+upperte
	atomdb,te_arr,dens_arr,tec3995=tec3995,$
		               tec4041=tec4041,$
		               tec4026=tec4026
	dens3 = interpol(dens_arr,tec4041/tec3995,exp_ratio1_lower)
	dens4 = interpol(dens_arr,tec4041/tec3995,exp_ratio1_upper)

	dens_lower = fltarr(n_elements(data.time))
	dens_upper = fltarr(n_elements(data.time))
	te_lower   = fltarr(n_elements(data.time))+lowerte
	te_upper   = fltarr(n_elements(data.time))+upperte

	for i=0,n_elements(data.time)-1 do begin
		dens_lower[i] = dens3[i]<dens4[i]
		dens_upper[i] = dens1[i]>dens2[i]		
	endfor

	atomdb,te_lower,dens_upper,tec3995=tec3995,$
		                   tec4041=tec4041,$
		                   tec4026=tec4026
	ratio1_upper = tec4041/tec3995
	ratio2_upper = tec4041/tec4026

	atomdb,te_upper,dens_lower,tec3995=tec3995,$
		                   tec4041=tec4041,$
		                   tec4026=tec4026
	ratio1_lower = tec4041/tec3995
	ratio2_lower = tec4041/tec4026

	; find concentration
        cn_lower = calc_cn(shot,los,transmission,te_lower,dens_upper,data,sm)
        cn_upper = calc_cn(shot,los,transmission,te_upper,dens_lower,data,sm)
	
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

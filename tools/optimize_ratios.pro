Function calc_profiles, data, dens_arr, te_arr, ratio1, ratio2, exp_ratio1, exp_ratio2 ,debug=debug

	num        = n_elements(dens_arr)
	idx        = fltarr(num)	
	min2       = fltarr(num)
	te         = exp_ratio1
	dens       = exp_ratio1
	sim_ratio1 = exp_ratio1
	sim_ratio2 = exp_ratio1
	for k=0,n_elements(exp_ratio1[0,*])-1 do begin
    	    for i=0,n_elements(exp_ratio1[*,0])-1 do begin
		for j=0,num-1 do begin
			min1   = abs((ratio1[*,j] - exp_ratio1[i,k])/exp_ratio1[i,k]) * $
	         		 abs((ratio2[*,j] - exp_ratio2[i,k])/exp_ratio2[i,k])
			idx[j] = where(min1 eq min(min1))
			min2[j]= min1[idx[j]]
		endfor
		idy     = where(min2 eq min(min2))
        	; check if reasonable solution found
		sim_ratio1[i,k] = ratio1[idx[idy[0]],idy[0]]
		sim_ratio2[i,k] = ratio2[idx[idy[0]],idy[0]]
		if abs(sim_ratio1[i,k]-exp_ratio1[i,k])/exp_ratio1[i,k] lt 0.05 and $
		   abs(sim_ratio2[i,k]-exp_ratio2[i,k])/exp_ratio2[i,k] lt 0.05 then begin
			dens[i,k] = dens_arr[idx[idy[0]]]
			te[i,k]   = te_arr[idy[0]]
		endif else begin
			dens[i,k]       = -1
			te[i,k]         = -1	
			sim_ratio1[i,k] = -1
			sim_ratio2[i,k] = -1
		end
	    endfor
	    id = where(dens[*,k] ne -1,complement=iy)
	    if id[0] ne -1 then begin
		te[iy,k]   = interpol(te[id,k],data.time[id],data.time[iy])
		dens[iy,k] = interpol(dens[id,k],data.time[id],data.time[iy])
		sim_ratio1[iy,k] = interpol(sim_ratio1[id,k],data.time[id],data.time[iy])
		sim_ratio2[iy,k] = interpol(sim_ratio2[id,k],data.time[id],data.time[iy])
	    endif	
	endfor

	return,{sim_ratio1: sim_ratio1,$
        	sim_ratio2: sim_ratio2,$
		dens      : dens      ,$
		te        : te         }
end

Function optimize_ratios,data,shot,los,transmission,sm=sm,debug=debug,nocn=nocn

	num      = 40
	te_arr   = adas_vector(high=6,low=2.,num=num)
	dens_arr = adas_vector(high=1e15,low=1e13,num=num)
	ratio1   = fltarr(num,num)
	ratio2   = fltarr(num,num)
	exp_ratio1 = data.nii4041
	exp_ratio2 = data.nii4041
	err_ratio1 = data.nii4041
	err_ratio2 = data.nii4041
	exp_ratio1_upper = data.nii4041
	exp_ratio2_lower = data.nii4041
	exp_ratio2_upper = data.nii4041
	exp_ratio1_lower = data.nii4041
    	for i=0,n_elements(data.nii4041[0,*])-1 do begin
    	    exp_ratio1[*,i] = smooth(data.nii4041[*,i],sm,/edge_truncate) / smooth(data.nii3995[*,i],sm,/edge_truncate)
	    exp_ratio2[*,i] = smooth(data.nii4041[*,i],sm,/edge_truncate) / smooth(data.nii4026[*,i],sm,/edge_truncate)
	    
	    err_ratio1[*,i] = exp_ratio1[*,i] * sqrt((data.nii3995_err[*,i]/data.nii3995[*,i])^2+(data.nii4041_err[*,i]/data.nii4041[*,i])^2)
	    err_ratio2[*,i] = exp_ratio2[*,i] * sqrt((data.nii4026_err[*,i]/data.nii4026[*,i])^2+(data.nii4041_err[*,i]/data.nii4041[*,i])^2)
	
	    exp_ratio1_upper[*,i] = exp_ratio1[*,i] + err_ratio1[*,i]/2.0
	    exp_ratio1_lower[*,i] = exp_ratio1[*,i] - err_ratio1[*,i]/2.0

	    exp_ratio2_upper[*,i] = exp_ratio2[*,i] + err_ratio2[*,i]/2.0
	    exp_ratio2_lower[*,i] = exp_ratio2[*,i] - err_ratio2[*,i]/2.0
    	endfor

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

        ; case 1: ratio1 + err with ratio2 + err
        res1 = calc_profiles(data, dens_arr, te_arr, ratio1, ratio2, exp_ratio1_upper, exp_ratio2_upper,debug=debug)

        ; case 2: ratio1 + err with ratio2 - err
        res2 = calc_profiles(data, dens_arr, te_arr, ratio1, ratio2, exp_ratio1_upper, exp_ratio2_lower,debug=debug)

        ; case 3: ratio1 - err with ratio2 + err
        res3 = calc_profiles(data, dens_arr, te_arr, ratio1, ratio2, exp_ratio1_lower, exp_ratio2_upper,debug=debug)

        ; case 3: ratio1 - err with ratio2 - err
        res4 = calc_profiles(data, dens_arr, te_arr, ratio1, ratio2, exp_ratio1_lower, exp_ratio2_lower,debug=debug)

	; find maximum and minimum values for ratios
	te_lower     = exp_ratio1
	te_upper     = exp_ratio1
	dens_lower   = exp_ratio1
	dens_upper   = exp_ratio1
	ratio1_lower = exp_ratio1
	ratio2_lower = exp_ratio1
	ratio1_upper = exp_ratio1
	ratio2_upper = exp_ratio1
	for k=0,n_elements(exp_ratio1[0,*])-1 do begin
    	    for i=0,n_elements(data.time)-1 do begin
		te_lower[i,k]     = res1.te[i,k]<res2.te[i,k]<res3.te[i,k]<res4.te[i,k]
		te_upper[i,k]     = res1.te[i,k]>res2.te[i,k]>res3.te[i,k]>res4.te[i,k]
		dens_lower[i,k]   = res1.dens[i,k]<res2.dens[i,k]<res3.dens[i,k]<res4.dens[i,k]
		dens_upper[i,k]   = res1.dens[i,k]>res2.dens[i,k]>res3.dens[i,k]>res4.dens[i,k]
		ratio1_lower[i,k] = res1.sim_ratio1[i,k]<res2.sim_ratio1[i,k]<res3.sim_ratio1[i,k]<res4.sim_ratio1[i,k]
		ratio1_upper[i,k] = res1.sim_ratio1[i,k]>res2.sim_ratio1[i,k]>res3.sim_ratio1[i,k]>res4.sim_ratio1[i,k]
		ratio2_lower[i,k] = res1.sim_ratio2[i,k]<res2.sim_ratio2[i,k]<res3.sim_ratio2[i,k]<res4.sim_ratio2[i,k]
		ratio2_upper[i,k] = res1.sim_ratio2[i,k]>res2.sim_ratio2[i,k]>res3.sim_ratio2[i,k]>res4.sim_ratio2[i,k]
	    endfor
	endfor

	; find line-integrated nii
	cn_lower = exp_ratio1
	cn_upper = exp_ratio1

        if ~keyword_set(nocn)then begin
    	    cn_1 = calc_cn(shot,los,transmission,te_lower,dens_lower,data,sm,err=err_1,/jet)
            cn_2 = calc_cn(shot,los,transmission,te_lower,dens_upper,data,sm,err=err_2,/jet)
            cn_3 = calc_cn(shot,los,transmission,te_upper,dens_lower,data,sm,err=err_3,/jet)
            cn_4 = calc_cn(shot,los,transmission,te_upper,dens_upper,data,sm,err=err_4,/jet)
	    for k=0,n_elements(exp_ratio1[0,*])-1 do begin
	    	for i=0,n_elements(data.time)-1 do begin
		    cn_lower[i,k] = (cn_1[i,k])<(cn_2[i,k])<(cn_3[i,k])<(cn_4[i,k])
		    cn_upper[i,k] = (cn_1[i,k])>(cn_2[i,k])>(cn_3[i,k])>(cn_4[i,k])
	    	endfor
	    endfor
	endif 
	
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

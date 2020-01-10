Pro calc_profs, nii3995, nii4041, nii4026 , te=te, dens=dens, calc=calc
	num        = 40
	exp_ratio1 = nii4041/nii3995
	exp_ratio2 = nii4041/nii4026
	
	te_arr   = adas_vector(high=6,low=1.,num=num)
	dens_arr = adas_vector(high=1e15,low=1e13,num=num)
	ratio1   = fltarr(num,num)
	ratio2   = fltarr(num,num)
	print,'Fetching atomic data...'
	if keyword_set(calc)then begin
		for i=0,num-1 do begin
			atomdb,te_arr[i]+fltarr(num),dens_arr,$
			               	tec3995=tec3995,$
		        	       	tec4041=tec4041,$
		               		tec4026=tec4026
			ratio1[*,i] = tec4041 / tec3995
			ratio2[*,i] = tec4041 / tec4026
		endfor
		save,file='tmp/ratio_analysis.sav',ratio1,ratio2,dens_arr,te_arr
	endif else restore,'tmp/ratio_analysis.sav'

	tol = 0.05
	num        = n_elements(dens_arr)
	idx        = fltarr(num)	
	min2       = fltarr(num)
	for j=0,num-1 do begin
		min1   = abs((ratio1[*,j] - exp_ratio1)/exp_ratio1) * $
	 		 abs((ratio2[*,j] - exp_ratio2)/exp_ratio2)
		idx[j] = where(min1 eq min(min1))
		min2[j]= min1[idx[j]]
	endfor

	; check for multiple solutions and choose solution that matches best ratio
	;
	
	order = sort(min2)
	chi = fltarr(n_elements(order))
	for k=0,n_elements(order)-1 do begin
		chi[k] = (exp_ratio1 - ratio1[idx[order[k]],order[k]])^2/exp_ratio1
		chi[k] = chi[k] +  (exp_ratio2 - ratio2[idx[order[k]],order[k]])^2/exp_ratio2
	endfor
	
	idz = where(chi eq min(chi))
	idy = order[idz[0]]
        ; check if reasonable solution found
	
	sim_ratio1 = ratio1[idx[idy[0]],idy[0]]
	sim_ratio2 = ratio2[idx[idy[0]],idy[0]]
	
	if abs(sim_ratio1-exp_ratio1)/exp_ratio1 lt tol and $
	   abs(sim_ratio2-exp_ratio2)/exp_ratio2 lt tol then begin
		dens = dens_arr[idx[idy[0]]]
		te   = te_arr[idy[0]]
		atomdb,te,dens,$
	               tec3995=tec3995,$
		       tec4041=tec4041,$
	       	       tec4026=tec4026
			
		dl = 0.05
		cn = nii3995 * 4.00 * !pi  /( tec3995 * dens ) / (dens * 1e6) / dl
		
		print,'----------------------------------------'
		print,'Input ratio1: ',exp_ratio1
		print,'Input ratio2: ',exp_ratio2
		print,'Te [eV] = ',te
		print,'dens [cm-3] = ',dens
		print,'dL = 5 cm'
		print,'cN [%] = ',cn *100
		print,'Input ratio1/Sim ratio1 = ',exp_Ratio1/sim_ratio1
		print,'Input ratio2/Sim ratio2 = ',exp_Ratio2/sim_ratio2		
		print,'----------------------------------------'
		
	endif else begin
		print,'----------------------------------------'
		print,'No solution found :( '
		print,'----------------------------------------'
	end
	
end

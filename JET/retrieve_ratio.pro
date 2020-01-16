Function retrieve_ratio,nii,rvals,horizontal=horizontal,data=data

    ; nii = fltarr(time,los)
    ; rvals = fltarr(los)
    
    if keyword_set(data)then begin
    	nii = data.nii[*,*,1]
	rvals = data.rvals
    endif
    
    nt   = n_elements(nii[*,0])

    ; Average strike-point and x-point areas
    
    if ~keyword_set(horizontal)then begin
        id_sp = where(rvals ge 2.82 and rvals le 2.87)
    	id_xp = where(rvals ge 2.75 and rvals le 2.80)
    endif else begin
        id_sp = where(rvals ge 2.68 and rvals le 2.72)
    	id_xp = where(rvals ge 2.64 and rvals le 2.67)        
    end

    avr_sp = fltarr(nt)
    avr_xp = fltarr(nt)
        
    for i=0,nt-1  do avr_sp[i] = mean(nii[i,id_sp])     
    for i=0,nt-1  do avr_xp[i] = mean(nii[i,id_xp])     

    return, smooth(avr_sp / avr_xp,2)
End

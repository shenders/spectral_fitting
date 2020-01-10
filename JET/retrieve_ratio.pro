Function retrieve_ratio,nii,los,shot

    ; nii = fltarr(time,los)
    ; los = ['SP1','SP2',etc]

    sp   = fltarr(n_elements(los)) & for i=0,n_elements(los)-1 do sp[i]  = float(strmid(los[i],2,2))
    nt   = n_elements(nii[*,0])

    ; Average strike-point and x-point areas
    ; note there is a difference in fibre location between 85### and 96###

    if shot gt 96000 then begin
        ids    = [17,18,19]
    	idx    = [13,14,15]
    endif else begin
        ids    = [16,17,18]
    	idx    = [12,13,14]    
    end
    id_sp  = fltarr(n_elements(ids))
    id_xp  = fltarr(n_elements(idx))
    avr_sp = fltarr(nt)
    avr_xp = fltarr(nt)
    
    for i=0,n_elements(ids)-1 do id_sp[i]  = where(sp eq ids[i])
    for i=0,n_elements(idx)-1 do id_xp[i]  = where(sp eq idx[i])
    for i=0,nt-1  do avr_sp[i] = mean(nii[i,id_sp])     
    for i=0,nt-1  do avr_xp[i] = mean(nii[i,id_xp])     

    return, smooth(avr_sp / avr_xp,2)
End

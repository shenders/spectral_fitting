Pro combine_stark,shot,tr=tr,los,load=load,rescale=rescale
    
    if ~keyword_set(load)then save=1

    x=fetch_data(shot,los,load=load,save=save,tr=tr,append='stark',machine='JET',diag='KT3B',/stark)
    y=fetch_data(shot,los,load=load,save=save,tr=tr,machine='JET',diag='KT3B',/quick)
    
    y.nii[*,*,3]     = x.nii[*,*,3]
    y.nii_err[*,*,3] = x.nii_err[*,*,3]
    
    if keyword_set(rescale)then begin
    	; first store saturated 399.5nm line as array 4
    	y.nii[*,*,4]     = y.nii[*,*,0]
    	y.nii_err[*,*,4] = y.nii_err[*,*,0]
    	; use branching ratio of 395.5 line to set 399.5 line
    	branching_ratio  = 1.0/0.0937
	y.nii[*,*,0]     = x.nii[*,*,3] * branching_ratio
    	y.nii_err[*,*,0] = x.nii_err[*,*,3] * branching_ratio
    endif
    
    y.dens_balmer    = x.dens_balmer
    y.dens_balmer_err= x.dens_balmer_err

    shotstr = string(shot,format='(i5)')
    append = 'data'
    output = y
    file = 'save/'+shotstr+'/'+los+'-'+append+'.idl'
    save,file=file,output

    stop
end

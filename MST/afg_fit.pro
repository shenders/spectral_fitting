FUNCTION afg_fit,data,use_tau=use_tau,$
		      debug=debug,$
		      tr=tr,$
		      smoothfunc=smoothfunc,$
		      mdlfile=mdlfile,$
		      noprog=noprog
	
	
	
	
	if ~keyword_set(use_tau)then use_tau = 200
	!quiet  = 1
	!except = 0
	time    = data.time
	ntime   = n_elements(time)
	if (size(data.nii))[0] gt 2 then begin
		nlines  = n_elements(data.nii(0,0,*))
		nii     = reform(data.nii[*,0,*],ntime,nlines)
		err     = reform(data.nii_err[*,0,*],ntime,nlines)
	endif else begin
		nlines  = n_elements(data.nii(0,*))
		nii     = data.nii[*,*]
		err     = data.nii_err[*,*]
	end

	if keyword_set(smoothfunc)then begin
		for i=0,n_elements(nii(0,*))-1 do nii(*,i)=smooth(nii(*,i),smoothfunc)
	endif
	if ~keyword_Set(smoothfunc)then smoothfunc = -1

	if keyword_set(tr)then begin
		id_reduce = where(time ge tr[0] and time le tr[1])
		time      = time[id_reduce]
		nii       = nii[id_reduce,*]
	endif	
	
	ntime        = n_elements(time)
	dens         = fltarr(ntime)
	temp         = fltarr(ntime)
	nconc        = fltarr(ntime)
	dens_err     = fltarr(ntime)
	temp_err     = fltarr(ntime)
	nconc_err_lw = fltarr(ntime)
	nconc_err_up = fltarr(ntime)
	vac_to_air   = 0.113
	nist_399     = [399.613                                         ]-vac_to_air
	nist_402     = [402.722 , 404.048                               ]-vac_to_air
	nist_404     = [403.622 , 404.245 , 404.467 , 404.592 , 405.805 ]-vac_to_air
	nist_407     = [407.420 , 407.808 , 408.342 , 409.706 ]-vac_to_air
	nist_408     = [408.846 ]-vac_to_air
    	wblock       = [ 4074.20  , 4078.08  , 4083.42  , 4097.06  ]
	norm         = [ 0.387051 , 0.147080 , 0.328994 , 0.136876 ]
	fwhm         = 0.085
	sigma        = fwhm / 2.35482
	factor       = sqrt(2.0 * 3.141)
	wavelength   = findgen(1024) * (407.0 - 399.0) / 1023.0 + 399.0
	skip399      = 0
	skip402      = 0
	skip404      = 0
	skip407      = 1
	skip408      = 1
	icount_total = ntime*1.0 
	for i=0,ntime-1 do begin
		spec         = 0
		spec_err     = 0
    		if skip399 eq 0 then begin
			norm     = [ 1.0 ]
			scal     = nii(i,0) * norm / sigma / factor
			scal_err = err(i,0) * norm / sigma / factor
			for j=0,n_elements(norm)-1 do spec     = spec + scal[j] * $
							exp(-(wavelength - nist_399[j])^2/2.0/sigma^2)
			for j=0,n_elements(norm)-1 do spec_err = spec_err+ scal_err[j] * $
							exp(-(wavelength - nist_399[j])^2/2.0/sigma^2)
		endif	
    		if skip402 eq 0 then begin
			norm     = [ 0.92675 , 0.0732504 ]
			scal     = nii(i,2) * norm / sigma / factor
			scal_err = err(i,2) * norm / sigma / factor
			for j=0,n_elements(norm)-1 do spec     = spec + scal[j] * $
							exp(-(wavelength - nist_402[j])^2/2.0/sigma^2)
			for j=0,n_elements(norm)-1 do spec_err = spec_err+ scal_err[j] * $
							exp(-(wavelength - nist_402[j])^2/2.0/sigma^2)
		endif	
    		if skip404 eq 0 then begin
			norm     = [ 0.207273 , 0.549047 , 0.178974 , 0.0322398 , 0.0324669]
			scal     = nii(i,1) * norm / sigma / factor
			scal_err = err(i,1) * norm / sigma / factor
			for j=0,n_elements(norm)-1 do spec     = spec + scal[j] * $
							exp(-(wavelength - nist_404[j])^2/2.0/sigma^2)
			for j=0,n_elements(norm)-1 do spec_err = spec_err+ scal_err[j] * $
							exp(-(wavelength - nist_404[j])^2/2.0/sigma^2)
		endif	
    		if skip407 eq 0 then begin
			norm     = [ 0.387051 , 0.147080 , 0.328994 , 0.136876 ]
			scal     = nii(i,3) * norm / sigma / factor
			scal_err = err(i,1) * norm / sigma / factor
			for j=0,n_elements(norm)-1 do spec     = spec + scal[j] * $
							exp(-(wavelength - nist_407[j])^2/2.0/sigma^2)
			for j=0,n_elements(norm)-1 do spec_err = spec_err+ scal_err[j] * $
							exp(-(wavelength - nist_407[j])^2/2.0/sigma^2)
		endif	
    		if skip408 eq 0 then begin
			norm     = [ 1.0]
			scal     = nii(i,4) * norm / sigma / factor
			scal_err = err(i,1) * norm / sigma / factor
			for j=0,n_elements(norm)-1 do spec     = spec + scal[j] * $
							exp(-(wavelength - nist_408[j])^2/2.0/sigma^2)
			for j=0,n_elements(norm)-1 do spec_err = spec_err+ scal_err[j] * $
							exp(-(wavelength - nist_408[j])^2/2.0/sigma^2)
		endif	
		bconst     = 0.001
		background = bconst * max(spec)
		spec       = spec   + background
		yerr       = spec_err + background*0.1
		params     = run_ffs_fit(wavelength,spec,yerr=yerr,$
			     mdlfile=mdlfile,instr_func=fwhm,$
    			     background=bconst,fixback=1,debug=debug,$
			     psplot=psplot,use_tau=use_tau)
		dens[i]         = params.dens
		temp[i]         = params.te
		nconc[i]        = params.nconc
		dens_err[i]     = params.ne_err
		temp_err[i]     = params.te_err
		nconc_err_up[i] = params.nconc_err_up 
		nconc_err_lw[i] = params.nconc_err_lw
		if ~keyword_set(noprog)then begin
			IF i EQ 0 THEN progress,0.0,/reset,label='Progress (%)',frequency=1.0
			icount = i*1.0
			progress,100.0*icount/icount_total
		endif	
	endfor
	if ~keyword_set(noprog)then progress,100.0,/last
	
	output = {	tau:use_tau,$
			nii:nii,$
			err:err,$
			time:time,$
			te:temp,$
			te_err:temp_err,$
			dens:dens,$
			dens_err:dens_err,$
			nconc:nconc,$
			nconc_err_up:nconc_err_up,$
			nconc_err_lw:nconc_err_lw,$
			smoothfunc:smoothfunc}
			
RETURN,output
END

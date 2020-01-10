Function run_ffs_fit, x, spectrum,fixback=fixback,yerr=yerr,$
                      mdlfile=mdlfile,dens=dens,nomodel=nomodel,$
		      psplot=psplot,debug=debug,$
		      instr_func=instr_func,gauss=gauss,voigt=voigt,$
		      background=background,$
		      backwave=backwave,use_tau=use_tau

;	************************************
;	**** Initialise data            ****
;	************************************
	
	id   = where(spectrum gt 0 and finite(spectrum) eq 1)
	y    = spectrum[id]
    	x    = x[id]
	IF KEYWORD_SET(yerr)then yerr = yerr[id]	
	ynorm          = where(x lt 400)
	if ynorm[0] eq -1 then maxy=max(y) else maxy = max(y[ynorm])
	y              = y/maxy
	if ~keyword_set(background)then begin
		nbins          = 10
        	binarr         = FINDGEN(nbins)*(MAX(ALOG10(y)-MIN(ALOG10(y))))/(nbins-1)+MIN(ALOG10(y))
		hist           = HISTOGRAM(ALOG10(y),nbins=nbins)
		aval           = [MAX(hist),binarr[WHERE(hist EQ MAX(hist))],0.5]
		ga_back        = GAUSSFIT(binarr,hist,aval,nterms=3)
		hist_back      = 10^(aval[1])
		IF hist_back LT min(y) THEN BEGIN
		    miny = MIN(y)
		    id = where(y lt 3.0*miny)
		    backc=mean(y[id])
		ENDIF ELSE backc          = hist_back
		backc          = hist_back
	endif else backc=background	
	if backc lt 0.5 * min (y) then backc = min(y)
	IF ~KEYWORD_SET(fixback)THEN fixback  = 0

	instr_func_str = STRTRIM(STRING(instr_func,FORMAT='(D6.3)'),2)
    	background_str = STRTRIM(STRING(backc,FORMAT='(D10.7)'),2)
    	background_up  = STRTRIM(STRING(backc*1.25,FORMAT='(D10.7)'),2)
    	background_lw  = STRTRIM(STRING(backc*0.8,FORMAT='(D10.7)'),2)
	IF ~KEYWORD_SET(yerr)then err  = FLTARR(N_ELEMENTS(y))+backc*0.1 ELSE err  = yerr/maxy
;	************************************
;	**** Setup MDL file             ****
;	************************************
	for i=100,128 do free_lun,i
	GET_LUN,unit_write & IF KEYWORD_SET(mdlfile)THEN mdl_file=mdlfile ELSE mdl_file = 'tmp/default_model.mdl'
	IF KEYWORD_SET(nowrite)THEN OPENW,unit_write,'temp.mdl' ELSE OPENW,unit_write,mdl_file
    	junk= ''
    	PRINTF,unit_write,'(model'
    	PRINTF,unit_write,'    (+'

;	************************************
;	**** AFG details                ****
;	************************************

	IF KEYWORD_SET(use_tau)THEN afg_file = 'adasn1_t' else afg_file = 'adasn1_r' 	
	IF ~KEYWORD_SET(nomodel)THEN PRINTF,unit_write,'       (* (broaden-gaussian (adas-'+afg_file+' nspec) nbroad ) mul)'
    	PRINTF,unit_write,'       (background-linear backg)'
	
;	************************************
;	**** Get info                   ****
;	************************************

	
	IF N_ELEMENTS(x) LT 2 THEN GOTO,nofit
	IF ~KEYWORD_SET(gauss) THEN ngauss = 0 ELSE BEGIN
	    ngauss     = N_ELEMENTS(gauss.pos)
    	    gauss_str  = REPLICATE({pos:'',area:'',areaup:'',arealw:'',posup:'',poslw:'',fwhm:'',name:''},ngauss)
	ENDELSE
	IF ~KEYWORD_SET(voigt) THEN nvoigt = 0 ELSE BEGIN
	    nvoigt = N_ELEMENTS(voigt.pos)
    	    voigt_str  = REPLICATE({pos:'',area:'',areaup:'',arealw:'',posup:'',poslw:'',fwhmg:'',fwhml:''},nvoigt)
	ENDELSE

	FOR i=0,ngauss-1 DO BEGIN
	    gauss_str(i).pos   = STRTRIM(STRING(gauss.pos(i),FORMAT='(D10.5)'),2)    
	    gauss_str(i).name  = STRCOMPRESS(STRLOWCASE(gauss.ion(i)),/REMOVE_ALL)
	    gauss_str(i).name  = gauss_str(i).name+STRING(gauss.pos(i)*100,FORMAT='(I5)')+gauss.trn(i)
	    
	    if keyword_Set(calwave) then begin
	    	gauss_str(i).posup = STRTRIM(STRING(gauss.pos(i)+0.1,FORMAT='(D10.5)'),2)    
	    	gauss_str(i).poslw = STRTRIM(STRING(gauss.pos(i)-0.1,FORMAT='(D10.5)'),2)    
	    endif else begin
	    	gauss_str(i).posup = STRTRIM(STRING(gauss.pos(i)+0.05,FORMAT='(D10.5)'),2)    
	    	gauss_str(i).poslw = STRTRIM(STRING(gauss.pos(i)-0.05,FORMAT='(D10.5)'),2)    
    	    endelse
;	************************************
;	**** Estimate area guess        ****
;	************************************
    	    inten      = MAX((INTERPOL(y,x,FINDGEN(10)*0.2/9.0+gauss.pos(i)-0.005)))-backc
	    sigma      = gauss.fwhm / 2.35482
	    area_guess = inten * sigma * SQRT(2.0*3.141) 
	    area_guess = area_guess > 1E-5
	    gauss_str(i).area  = STRTRIM(STRING(area_guess,FORMAT='(D10.5)'),2)    
	    gauss_str(i).areaup= STRTRIM(STRING(area_guess*2.5,FORMAT='(D10.5)'),2)    
	    gauss_str(i).arealw= STRTRIM(STRING(area_guess*0.1,FORMAT='(D10.5)'),2)    
	    gauss_str(i).fwhm  = STRTRIM(STRING(gauss.fwhm,FORMAT='(D10.5)'),2)
	    str                = STRTRIM(STRING(i+1,FORMAT='(I4)'),2)
	    PRINTF,unit_write,'       (gaussian '+gauss_str(i).name+')' 	     
	    
	ENDFOR

	FOR i=0,nvoigt-1 DO BEGIN
	    voigt_str(i).pos   = STRTRIM(STRING(voigt.pos(i),FORMAT='(D10.5)'),2)    
	    voigt_str(i).posup = STRTRIM(STRING(voigt.pos(i)+0.03,FORMAT='(D10.5)'),2)    
	    voigt_str(i).poslw = STRTRIM(STRING(voigt.pos(i)-0.03,FORMAT='(D10.5)'),2)    
	    inten      = MAX((INTERPOL(y,x,FINDGEN(10)*0.2/9.0+voigt.pos(i)-0.1)))-backc
	    sigma      = voigt.fwhmg / 2.35482
	    area_guess = inten * sigma * SQRT(2.0*3.141) * 3.0 
	    area_guess = area_guess > 0.0
	    voigt_str(i).area  = STRTRIM(STRING(area_guess,FORMAT='(D10.5)'),2)    
	    voigt_str(i).areaup= STRTRIM(STRING(area_guess*2.5,FORMAT='(D10.5)'),2)    
	    voigt_str(i).arealw= STRTRIM(STRING(area_guess*0.01,FORMAT='(D10.5)'),2)    
	    voigt_str(i).fwhml = STRTRIM(STRING(voigt.fwhml(i),FORMAT='(D10.2)'),2) 
	    voigt_str(i).fwhmg = STRTRIM(STRING(voigt.fwhmg,FORMAT='(D10.5)'),2)    
	    str                = STRTRIM(STRING(i+1,FORMAT='(I4)'),2)
	    strname            = STRING(voigt.pos(i)*10,FORMAT='(I4)')+STRING(voigt.nbalmer[i],FORMAT='(i1)')+'2'  
	    PRINTF,unit_write,'       (broaden_lorentz52 (gaussian hi'+strname+') linebl'+str+')'
	ENDFOR
	

    	PRINTF,unit_write,'     )'
	PRINTF,unit_write,'adas_fit)'
    	PRINTF,unit_write,'(setval backg.m 0.0000)' 
    	PRINTF,unit_write,'(fixed backg.m)' 
    	PRINTF,unit_write,'(setval backg.c '+background_str+')'
        IF KEYWORD_SET(fixback)THEN PRINTF,unit_write,'(fixed backg.c)' ELSE $
	                            PRINTF,unit_write,'(setlimits backg.c '+background_lw+' '+background_up+')'
	IF ~KEYWORD_SET(nomodel)THEN BEGIN
		PRINTF,unit_write,'(setval nspec.te 4.5)(setlimits nspec.te 2.0 10.0)'
		PRINTF,unit_write,'(setval nspec.dens 1.00E14)(setlimits nspec.dens 1e13 1e15)' 
		PRINTF,unit_write,'(setval nspec.norm 1)'
		PRINTF,unit_write,'(setval nspec.ionbal 1)'
		IF KEYWORD_SET(use_tau)THEN PRINTF,unit_write,'(setval nspec.tau '+STRING(use_tau,format='(e8.2)')+')(fixed nspec.tau)' 
		PRINTF,unit_write,'(setval nspec.nii 1e13)(fixed nspec.nii)' 
		PRINTF,unit_write,'(setval nspec.niii 1e13)(fixed nspec.niii)' 
		PRINTF,unit_write,'(setval nbroad.fwhm '+instr_func_str+')(fixed nbroad.fwhm)'
		PRINTF,unit_write,'(setval mul.factor 0.091) (setlimits mul.factor 0.001 1.0)'
	ENDIF
; Test Stark broadening measurement
        IF nvoigt GT 0 THEN BEGIN 
	PRINTF,unit_write,'(dummy d)'
        IF KEYWORD_SET(dens)THEN BEGIN
	        dens_str=STRTRIM(STRING(dens,FORMAT='(E8.2)'),2)
	        PRINTF,unit_write,'(setval d.dens '+dens_str+')'
        	PRINTF,unit_write,'(fixed d.dens)'
        ENDIF ELSE BEGIN
	    	PRINTF,unit_write,'(setval d.dens 3.0e20)'
        	PRINTF,unit_write,'(setlimits d.dens 1.0e19 1.0e21)'
	ENDELSE	
    	ENDIF	
; End update
	FOR i=0,ngauss-1 DO BEGIN
	    IF gauss.erc(i)NE -1 THEN BEGIN
	    str                = STRTRIM(STRING(i+1,FORMAT='(I4)'),2)
    	    PRINTF,unit_write,'(setval '+gauss_str(i).name+'.pos '+gauss_str(i).pos+')'
	    if strmid(gauss_str(i).name,0,2) eq 'xx' then begin 
		    PRINTF,unit_write,'(fixed '+gauss_str(i).name+'.pos )'
	    endif else begin
		    PRINTF,unit_write,'(setlimits '+gauss_str(i).name+'.pos '+gauss_str(i).poslw+' '+gauss_str(i).posup+')'
	    end
	    PRINTF,unit_write,'(setval '+gauss_str(i).name+'.fwhm '+gauss_str(i).fwhm+')' 
	    fix_fwhm=1
	    IF fix_fwhm THEN PRINTF,unit_write,'(fixed '+gauss_str(i).name+'.fwhm )' $
	                ELSE PRINTF,unit_write,'(setlimits '+gauss_str(i).name+'.fwhm 0.08 0.45)'
    	    
	    FOR ii=1,MAX(gauss.couple) DO BEGIN
	    	IF (gauss.couple(i) GT ii-1 and gauss.couple(i) LT ii)THEN BEGIN
		    idcoupled  = WHERE(gauss.couple EQ ii)
 		    str_couple = STRTRIM(STRING(idcoupled[0]+1,FORMAT='(I4)'),2)
		    mul_couple = STRTRIM(STRING(gauss.couple(i)-(ii-1),FORMAT='(D10.4)'),2)
;	    	    PRINTF,unit_write,'(coupled '+gauss_str(i).name+'.area (* '+gauss_str(i).name+str_couple+'.area '+mul_couple+'))'
	    	    PRINTF,unit_write,'(coupled '+gauss_str(i).name+'.area (* '+gauss_str(idcoupled[0]).name+'.area '+mul_couple+'))'
		    GOTO,coupling_complete
		ENDIF
	    ENDFOR
	    PRINTF,unit_write,'(setval '+gauss_str(i).name+'.area '+gauss_str(i).area+')'
	    PRINTF,unit_write,'(setlimits '+gauss_str(i).name+'.area 0.00005 '+gauss_str(i).areaup+')'
	    coupling_complete:
	    PRINTF,unit_write,''
	    ENDIF
	ENDFOR
	FOR i=0,nvoigt-1 DO BEGIN
	    str     = STRTRIM(STRING(i+1,FORMAT='(I4)'),2)
	    strname = STRING(voigt.pos(i)*10,FORMAT='(I4)')+STRING(voigt.nbalmer[i],FORMAT='(i1)')+'2'  
; Test Stark broadening measurement
	    IF voigt.nbalmer[i] EQ 6 THEN parms=['0.7149','3.954e-16']
	    IF voigt.nbalmer[i] EQ 7 THEN parms=['0.7120','6.258e-16']
	    IF voigt.nbalmer[i] EQ 8 THEN parms=['0.7159','7.378e-16']
	    IF voigt.nbalmer[i] EQ 9 THEN parms=['0.7177','8.947e-16']
	    IF voigt.nbalmer[i] GT 9 THEN parms=['0.7190','9.500e-16']
	    PRINTF,unit_write,'(couple linebl'+str+'.fwhm (* (^ d.dens '+parms[0]+') '+parms[1]+') )'	    
    	    PRINTF,unit_write,'(setval hi'+strname+'.pos '+voigt_str(i).pos+')'
	    PRINTF,unit_write,'(setlimits hi'+strname+'.pos '+voigt_str(i).poslw+' '+voigt_str(i).posup+' )'
    	    IF voigt.couple(i) NE 0 and voigt.couple(i) LT 1.0 THEN BEGIN
	    	idcoupled  = WHERE(voigt.couple EQ 1.0)
		str_couple = STRING(voigt.pos(idcoupled[0])*10,FORMAT='(I4)')+STRING(voigt.nbalmer[idcoupled[0]],FORMAT='(i1)')+'2'  
		mul_couple = STRTRIM(STRING(voigt.couple(i),FORMAT='(D10.4)'),2)
	    	PRINTF,unit_write,'(coupled hi'+strname+'.area (* hi'+str_couple+'.area '+mul_couple+'))'
	    ENDIF ELSE BEGIN
	    	PRINTF,unit_write,'(setval hi'+strname+'.area '+voigt_str(i).area+')'
		PRINTF,unit_write,'(setlimits hi'+strname+'.area 0.00005 '+voigt_str(i).areaup+')'
	    ENDELSE	    	
	    PRINTF,unit_write,'(setval hi'+strname+'.fwhm '+voigt_str(i).fwhmg+')'    	
	    PRINTF,unit_write,'(fixed hi'+strname+'.fwhm'+')'    
; End update
	    PRINTF,unit_write,''
	ENDFOR
	
	CLOSE,unit_write & FREE_LUN,unit_write
;	************************************
;	**** Setup FFS model            ****
;	************************************
        modelname    = STRMID(mdl_file,0,STRPOS(mdl_file,'.mdl'))
	model	     = obj_new('ffs_model', modelname=modelname)
        parser	     = obj_new('ffs_parser', file=mdl_file)    

;	************************************
;	**** Fit the data               ****
;	************************************
    	if x[0] gt x[1] then stop,'Error: wavelength axis must be in ascending order'
	tmp           = parser->apply(model)
    	tmp           = model->setxdata(x)
    	tmp           = model->pdsetup()
    	init          = (*(model->getresult())).intensity
    	xset              = (*(model->getresult())).wavelength
	orig_parvals      = model->getparvals()
	orig_parvalsfree  = model->getparvals(/free)
	IF KEYWORD_SET(debug)THEN BEGIN
	    adas_colors,colors=colors 
	    wvstr = STRTRIM(STRING(x[0],x[N_ELEMENTS(x)-1],FORMAT='("wv_",I3,"-",I3)'),2)   	    	    
	    !p.thick=1.0
	    IF KEYWORD_SET(psplot)THEN makeps,XS=8,YS=5,file='debug_'+wvstr+'.ps' ELSE WINDOW,0
	    PLOT,x,y*maxy,XTITLE='Wavelength / nm',YTITLE='ph/s/m2/sr/nm',yr=[MIN(y*maxy),MAX(y*maxy)],/YLOG,XSTY=1,/NODATA  ;& OPLOT,xset,init*maxy,linest=5,thick=3.0
	    user_psym,1 & OPLOT,x,y*maxy;,PSYM=8
	ENDIF
	fit           = obj_new('ffs_fit',debug=debug)
    	data          = {x:x, y:y, error:err}
	tmp           = fit->apply(model, data)
        parnames_full = model->getparnames(/full)
        parnames_err  = model->getparnames(/full,/free)
        parvals       = model->getparvals()
    	parvalserr    = fit->geterrors() 
        parvalsfree   = fit->getp()
	if n_elements(parnames_err) gt n_elements(parvalserr) then parnames_err=parnames_err[1:*]
		
;	************************************
;	**** Return the background      ****
;	************************************
    	mback = parvals(WHERE(parnames_full eq 'backg.m'))
    	cback = parvals(WHERE(parnames_full eq 'backg.c'))
    	cbackerr = parvalserr(WHERE(parnames_err eq 'backg.c'))
	IF ~KEYWORD_SET(backwave)THEN backwave=(MIN(x)+MAX(x))/2.0
	backc = INTERPOL((cback[0] + mback[0] * (x - MIN(x))) * maxy,x,backwave)
	continuum  = (cback[0] + mback[0] * (x - MIN(x)))

;	************************************
;	**** Return the  parameters     ****

	res  = *(model->getresult())

	continuum = (cback[0] + mback[0] * (x - MIN(x)))
	IF KEYWORD_SET(debug)THEN BEGIN
    	    !p.thick=2.0
	    oplot, res.wavelength, res.intensity*maxy
	    oplot, x , continuum*maxy,col=colors.green
    	ENDIF

	func = 0.0
	FOR i=0,ngauss-1 DO BEGIN
	    IF gauss.erc(i)NE -1 THEN BEGIN
	    	str   = STRTRIM(STRING(i+1,FORMAT='(I4)'),2)
	    	fwhm  = parvals(WHERE(parnames_full eq gauss_str(i).name+'.fwhm'))
	    	sigma = fwhm[0] / 2.35482
	    	area  = parvals(WHERE(parnames_full eq gauss_str(i).name+'.area'))
	    	inten = area[0] / sigma[0] / SQRT(2.0 * 3.141)
	    	pos   = parvals(WHERE(parnames_full eq gauss_str(i).name+'.pos'))
	    	func  = func + inten[0] * EXP(-(x-pos[0])^2/2.0/sigma[0]^2)
	    	IF KEYWORD_SET(debug)THEN BEGIN
	    		OPLOT,x,inten[0]*maxy* EXP(-(x-pos[0])^2/2.0/sigma[0]^2),col=colors.red
	    	ENDIF
	    ENDIF
	ENDFOR
	resid = res.intensity - func - continuum
	IF KEYWORD_SET(debug)THEN BEGIN
	    OPLOT,x,resid*maxy,col=colors.blue
	ENDIF	
	IF KEYWORD_SET(debug)THEN BEGIN
	    IF KEYWORD_SET(psplot)THEN makeps,XS=8,YS=3,file='debug_'+wvstr+'_resid.ps' ELSE WINDOW,1
	    fitted=(*(model->getresult())).intensity	    
	    user_psym,1,/fill & PLOT,x,(y - fitted)/y,yr=[-1,1],xs=1,psym=8
	    oplot,x,res.intensity,col=colors.green,thick=1.0
	    oplot,x,-res.intensity,col=colors.green,thick=1.0
	ENDIF
	if ~keyword_set(nomodel)then begin
	    te_n1         = (parvals(WHERE(parnames_full EQ 'nspec.te')))[0]
	    ne_n1         = (parvals(WHERE(parnames_full EQ 'nspec.dens')))[0]
	    te_n1_err     = (parvalserr(WHERE(parnames_err EQ 'nspec.te')))[0]
	    ne_n1_err     = (parvalserr(WHERE(parnames_err EQ 'nspec.dens')))[0]
            obj           = obj_new('afg_'+afg_file)
            pars          = obj->getpars()
            pars.te       = te_n1
            pars.dens     = ne_n1
            pars.norm     = 0
            ok            = obj->setpars(pars=pars)
            result        = obj->getcalc()

            pars.te       = te_n1+te_n1_err
            pars.dens     = ne_n1+ne_n1_err
            pars.norm     = 0
            ok            = obj->setpars(pars=pars)
            err_up        = obj->getcalc()

            pars.te       = te_n1-te_n1_err
            pars.dens     = ne_n1-ne_n1_err
            pars.norm     = 0
            ok            = obj->setpars(pars=pars)
            err_lw        = obj->getcalc()
    	    
	    
    	    const         = 2.0*SQRT(2.0*alog(2.0))
	    sigma_sum     = (instr_func / const)^2
	    norm	  = 1.0/(SQRT(2*!PI*sigma_sum))
	    spec          = FLTARR(N_ELEMENTS(x))
	    
	    FOR i=0,N_ELEMENTS(result.wv)-1 DO spec[*]  = spec[*]+result.intensity(i)  *  norm   * EXP(-((x-result.wv(i)/10.0)^2/(2.0*sigma_sum)))
	    
	    id        = WHERE(x GE 399.0 AND x LE 400.0)
	    
	    IF id[0] ne -1 THEN BEGIN
	    	int_399      = INT_TABULATED(x[id],y[id]*maxy-backc)
	    	nconc        = int_399 * 4.0 * 3.141 / result.intensity[0] ; m-2
    	    	nconc_err_up = int_399 * 4.0 * 3.141 / err_up.intensity[0] ; m-2
    	    	nconc_err_lw = int_399 * 4.0 * 3.141 / err_lw.intensity[0] ; m-2
    	    ENDIF ELSE BEGIN
	    	nconc        = -1
    	    	nconc_err_up = -1
    	    	nconc_err_lw = -1
    	    ENDELSE	    	
	    
	    n_395_int = -1
	    n_399_int = -1
	    n_wi_int  = -1
	    n_402_int = -1
	    n_404_int = -1
	    n_408_int = -1
	    n_409_int = -1
	    Ar_750_int = -1
	    Ar_751_int = -1
	    n_nvi_int = -1
	    n_niii_int= -1
	    n_niii2_int= -1
	    h_72_int  = -1
	    h_62_int  = -1
	    n_399_err = -1
	    n_395_err = -1
	    n_wi_err  = -1
	    n_402_err = -1
	    n_404_err = -1
	    n_408_err = -1
	    n_409_err = -1
	    Ar_750_err = -1
	    Ar_751_err = -1
	    n_nvi_err = -1
	    n_niii_err= -1
	    n_niii2_err= -1
	    h_72_err = -1
	    h_62_err = -1
	    n_v_int   = -1
	    n_v_err   = -1
	    ne_ii_int = -1
	    ne_ii_err = -1
	    ne_iv_int = -1
	    ne_iv_err = -1
	    
    	endif else begin

	    wave      = 395.58
	    n_id      = where(strpos(parnames_full,'.pos') ne -1 and strpos(parnames_full,'n2') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    n_395_int = parvals[where(parnames_full EQ str[0])]*maxy
	    n_395_err = parvalserr[where(parnames_err EQ str[0])]*maxy
	    	    
	    wave      = 399.5
	    n_id      = where(strpos(parnames_full,'.pos') ne -1 and strpos(parnames_full,'n2') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    n_399_int = parvals[where(parnames_full EQ str[0])]*maxy
	    n_399_err = parvalserr[where(parnames_err EQ str[0])]*maxy

	    wave      = 400.88
	    n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    n_wi_int  = parvals[where(parnames_full EQ str[0])]*maxy
	    n_wi_err  = parvalserr[where(parnames_err EQ str[0])]*maxy

	    wave      = [402.6,403.93] & n_402_int=0.0 
	    cple      = [1    ,0] 
	    for ii=0,n_elements(wave)-1 do begin
	    	n_id      = where(strpos(parnames_full,'.pos') ne -1 and strpos(parnames_full,'n2') ne -1)
		id_pos    = where(abs(parvals(n_id)-wave(ii)) eq min(abs(parvals(n_id)-wave(ii))))
	    	str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    	n_402_int = n_402_int+parvals[where(parnames_full EQ str[0])]*maxy
	    	if cple(ii) eq 1 then n_402_error = parvalserr[where(parnames_err EQ str[0])]/parvals[where(parnames_full EQ str[0])]
	    endfor
	    n_402_err = n_402_int * n_402_error[0]

	    wave      = [403.5,404.13,404.35,404.47,405.68] & n_404_int=0.0 
	    cple      = [0    ,1     ,0     ,0     ,0     ] 
	    for ii=0,n_elements(wave)-1 do begin
	    	n_id      = where(strpos(parnames_full,'.pos') ne -1 and strpos(parnames_full,'n2') ne -1)
	    	id_pos    = where(abs(parvals(n_id)-wave(ii)) eq min(abs(parvals(n_id)-wave(ii))))
	    	str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    	n_404_int = n_404_int+parvals[where(parnames_full EQ str[0])]*maxy
	    	if cple(ii) eq 1 then n_404_error = parvalserr[where(parnames_err EQ str[0])]/parvals[where(parnames_full EQ str[0])]
	    endfor
	    n_404_err = n_404_int * n_404_error[0]

	    wave      = 405.77
	    n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    n_nvi_int = parvals[where(parnames_full EQ str[0])]*maxy
	    n_nvi_err = parvalserr[where(parnames_err EQ str[0])]*maxy
	    
	    wave      = [407.30,407.69,408.21] & n_408_int=0.0 
	    cple      = [1     ,0     ,0     ] 
	    for ii=0,n_elements(wave)-1 do begin
	    	n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    	id_pos    = where(abs(parvals(n_id)-wave(ii)) eq min(abs(parvals(n_id)-wave(ii))))
	    	str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    	n_408_int = n_408_int+parvals[where(parnames_full EQ str[0])]*maxy
	    	if cple(ii) eq 1 then n_408_error = parvalserr[where(parnames_err EQ str[0])]/parvals[where(parnames_full EQ str[0])]
	    endfor
	    n_408_err = n_408_int * n_408_error[0]
	    
	    wave      = 408.7
	    n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    n_409_int = parvals[where(parnames_full EQ str[0])]*maxy
	    n_409_err = parvalserr[where(parnames_err EQ str[0])]*maxy
	    
	    wave      = 460.15
	    n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    n_v_int   = parvals[where(parnames_full EQ str[0])]*maxy
	    n_v_err   = parvalserr[where(parnames_err EQ str[0])]*maxy

	    wave      = 750.38
	    n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    Ar_750_int = parvals[where(parnames_full EQ str[0])]*maxy
	    Ar_750_err = parvalserr[where(parnames_err EQ str[0])]*maxy

	    wave      = 751.46
	    n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    Ar_751_int = parvals[where(parnames_full EQ str[0])]*maxy
	    Ar_751_err = parvalserr[where(parnames_err EQ str[0])]*maxy

	    wave      = 369.4
	    n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    ne_ii_int = parvals[where(parnames_full EQ str[0])]*maxy
	    ne_ii_err = parvalserr[where(parnames_err EQ str[0])]*maxy

	    wave      = 381.3
	    n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    ne_iv_int = parvals[where(parnames_full EQ str[0])]*maxy
	    ne_iv_err = parvalserr[where(parnames_err EQ str[0])]*maxy

	    wave      = [409.735,410.34] & n_niii_int=0.0 
	    cple      = [1      ,0     ] 
	    for ii=0,n_elements(wave)-1 do begin
	    	n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    	id_pos    = where(abs(parvals(n_id)-wave(ii)) eq min(abs(parvals(n_id)-wave(ii))))
	    	str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    	n_niii_int = n_niii_int+parvals[where(parnames_full EQ str[0])]*maxy
	    	if cple(ii) eq 1 then n_niii_error = parvalserr[where(parnames_err EQ str[0])]/parvals[where(parnames_full EQ str[0])]
	    endfor
	    n_niii_err = n_niii_int * n_niii_error[0]
	    
	    wave      = [399.86 ,400.35] & n_niii2_int=0.0 
	    cple      = [0      ,1     ] 
	    for ii=0,n_elements(wave)-1 do begin
	    	n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    	id_pos    = where(abs(parvals(n_id)-wave(ii)) eq min(abs(parvals(n_id)-wave(ii))))
	    	str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    	n_niii2_int = n_niii2_int+parvals[where(parnames_full EQ str[0])]*maxy
	    	if cple(ii) eq 1 then n_niii2_error = parvalserr[where(parnames_err EQ str[0])]/parvals[where(parnames_full EQ str[0])]
	    endfor
	    n_niii2_err = n_niii2_int * n_niii2_error[0]

	    wave      = 396.89
	    n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    h_72_int  = parvals[where(parnames_full EQ str[0])]*maxy
	    h_72_err  = parvalserr[where(parnames_err EQ str[0])]*maxy

	    wave      = 410.04
	    n_id      = where(strpos(parnames_full,'.pos') ne -1)
	    id_pos    = where(abs(parvals(n_id)-wave) eq min(abs(parvals(n_id)-wave)))
	    str       = strmid(parnames_full[n_id[id_pos]],0,strpos(parnames_full[n_id[id_pos]],'.pos'))+'.area'
	    h_62_int  = parvals[where(parnames_full EQ str[0])]*maxy
	    h_62_err  = parvalserr[where(parnames_err EQ str[0])]*maxy

	    te_n1=-1
	    ne_n1=-1
	    te_n1_err=-1
	    ne_n1_err=-1
	    nconc	 = -1
    	    nconc_err_up = -1
    	    nconc_err_lw = -1
	    
	endelse     
	IF KEYWORD_SET(debug)THEN BEGIN
    	    for ii=0,n_elements(parvalsfree)-1 do print,parnames_err(ii),parvalsfree(ii),orig_parvalsfree(ii)
	    if keyword_set(nomodel)then begin
		print,'N II line ratio details:'
		print,'4041/3995: ',n_404_int/n_399_int
		print,'4041/4026: ',n_404_int/n_402_int
		print,'395.5 nm:  ',n_395_int
		calc_profs, n_399_int[0], n_404_int[0],n_402_int[0] 
		IF nvoigt NE 0 THEN BEGIN
		    PRINT,'Balmer Feature Details:'
	    	    id_check=where(parnames_err EQ 'd.dens')
		    if id_check[0] ne -1 then PRINT,'ne / m-3    : ',parvals[where(parnames_full EQ 'd.dens')],' +/- ',parvalserr[where(parnames_err EQ 'd.dens')] else $
		    			      PRINT,'ne / m-3    : ',parvals[where(parnames_full EQ 'd.dens')]
	    	ENDIF
	    endif else begin
    	    	print,'Normalised by: ',maxy   
	    	PRINT,'N II Feature Details:'
	    	PRINT,'Te / eV     : ',te_n1,' +/- ',te_n1_err
	    	PRINT,'ne / m-3    : ',ne_n1*1e6,' +/- ',ne_n1_err*1e6
	    	PRINT,'dLnN1 (+err) / m-2 : ',MAX([nconc_err_up,nconc_err_lw])
	    	PRINT,'dLnN1 / m-2 : ',nconc
	    	PRINT,'dLnN1 (-err) / m-2 : ',MIN([nconc_err_up,nconc_err_lw])
	    	IF nvoigt NE 0 THEN BEGIN
    	    	    PRINT,'Balmer Feature Details:'
	    	    PRINT,'ne / m-3    : ',parvals[where(parnames_full EQ 'd.dens')],' +/- ',parvalserr[where(parnames_err EQ 'd.dens')]
	    	ENDIF
        	endelse
	stop
    	ENDIF

	return,{norm:maxy, te:te_n1, $
	                   te_err:te_n1_err,$
			   ne_err:ne_n1_err,$
			   dens:ne_n1,$
			   ari:ar_750_int,$
			   ari2:ar_751_int,$
			   ari_err:ar_750_err,$
			   ari2_err:ar_751_err,$
			   niii:n_niii_int,$
			   niii2:n_niii2_int,$
			   nvi:n_nvi_int,$
			   neii:ne_ii_int,$
			   neiv:ne_iv_int,$
			   neii_err:ne_ii_err,$
			   neiv_err:ne_iv_err,$
			   nwi:n_wi_int,$
			   niii_err:n_niii_err,$
			   niii2_err:n_niii2_err,$
			   nvi_err:n_nvi_err,$
			   nwi_err:n_wi_err,$
			   nconc:nconc,$
			   nconc_err_up:nconc_err_up,$
			   nconc_err_lw:nconc_err_lw,$
	                   h72:h_72_int,$
			   h62:h_62_int,$
	                   h72_err:h_72_err,$
			   h62_err:h_62_err,$
			   n395:n_395_int,$
			   n399:n_399_int,$
			   n402:n_402_int,$
	                   n395_err:n_395_err,$
	                   n399_err:n_399_err,$
			   n402_err:n_402_err,$
			   balmer_ne_err:parvalserr[where(parnames_full EQ 'd.dens')],$
			   balmer_ne:parvals[where(parnames_full EQ 'd.dens')],$
			   n404:n_404_int,$
			   n408:n_408_int,$
			   n409:n_409_int,$
			   nv:n_v_int,$
			   nv_err:n_v_err,$
			   n404_err:n_404_err,$
			   n408_err:n_408_err,$
			   n409_err:n_409_err,$
			   wavelength:x,$
			   balmerfit:resid,$
			   fitted:res.intensity,$
			   gaussians:func,$
			   background:cback, $
			   backgrounderr:cbackerr, $
			   backm:mback}
	nofit:
	return,{norm:maxy, te:-1, $
	                   te_err:-1,$
			   ne_err:-1,$
			   dens:-1,$
			   niii:-1,$
			   ari:-1,$
			   ari2:-1,$
			   ari_err:-1,$
			   ari2_err:-1,$
			   niii2:-1,$
			   nvi:-1,$
			   nwi:-1,$
			   niii_err:-1,$
			   niii2_err:-1,$
			   nvi_err:-1,$
			   nwi_err:-1,$
	                   h72:-1,$
			   h62:-1,$
	                   h72_err:-1,$
			   h62_err:-1,$
			   nv:-1,$
			   nv_err:-1,$
			   nconc:-1,$
			   nconc_err_up:-1,$
			   nconc_err_lw:-1,$
			   n395:-1,$
			   n399:-1,$
			   n402:-1,$
	                   n395_err:-1,$
	                   n399_err:-1,$
			   n402_err:-1,$
			   balmer_ne_err:-1,$
			   balmer_ne:-1,$
			   neii:-1,$
			   neiv:-1,$
			   neii_err:-1,$
			   neiv_err:-1,$
			   n404:-1,$
			   n408:-1,$
			   n409:-1,$
			   n404_err:-1,$
			   n408_err:-1,$
			   n409_err:-1,$
			   wavelength:-1,$
			   balmerfit:-1,$
			   fitted:-1,$
			   gaussians:-1,$
			   background:-1, $
			   backgrounderr:-1, $
			   backm:-1}

END

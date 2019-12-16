Function get_data,shot,machine,diag=diag,interelm=interelm,debug=debug,sig=sig,jet_coord=jet_coord

IF machine eq 'AUG' THEN BEGIN
;    !path = expand_path('+/afs/ipp/u/mcavedon/VS/lib:')+':'+!path
    IF ~KEYWORD_SET(diag)THEN diag       = 'EVS'
    IF ~KEYWORD_SET(wshift)THEN wshift   = 0.0
    shotstr = STRING(shot,FORMAT='(I5)')
	IF shot lt 33725 THEN BEGIN
		read_spec,shot,diag,data 
		red_emiss  = data.phflx
		red_wavel  = data.lamgrid
		red_los    = data.losnames
		red_time   = data.time
		red_wavel  = (red_wavel+wshift)
		red_error  = data.sy
		instr_func = data.FWHM_NM
		red_time   = red_time[0:N_ELEMENTS(red_emiss(0,0,*))-1]
		PRINT,'Time resolution: ',red_time(1)-red_time(0)
		IF ~KEYWORD_SET(resolution) THEN resolution   = 0.0001
		act_resol    = red_time(1)-red_time(0)

		if keyword_set(interelm)then begin
			telm       = find_elm(shot,red_time,red_emiss(*,0,0))
			if ~keyword_set(elmcond)then elmcond=4.5
			idback     = where(telm ge elmcond)
			if idback[0] ne -1 then begin
				red_time   = red_time[idback]
				red_emiss  = red_emiss[*,*,idback]
				red_error  = red_error[*,*,idback]
			endif
		endif
	ENDIF ELSE BEGIN
		read_xvs_diag,shot, diag, $
                    dwdp, exptime, sad2, ctsph, $
                    time, lam, offset, sens, spec, $
                    wlen, wslit, gratcons, op_ang, pixw,$
                    magnification, fwhm_pix, $
                    neon_done, neon, lambda_neon, $
                    r1,phi1,z1, r2,phi2,z2, $
                    error=error, $
                    los_name= los_name, $
                    dat = dat, $ 
		    exp_name = exp_name,$
                    no_copy=no_copy, $
                    no_smear_cor = no_smear_cor, $
                    read_again = read_again
		
		red_emiss  = spec
		IF diag NE 'DVS' THEN BEGIN
			red_emiss  = TRANSPOSE(red_emiss,[0,2,1])
			FOR i=0,n_elements(time)-1 DO $
		   		red_emiss[*,*,i] = red_emiss[*,*,i] / sens[*,*] / exptime / dwdp 
		ENDIF ELSE BEGIN
			FOR i=0,n_elements(time)-1 DO $
		   		red_emiss[*,i] = red_emiss[*,i] / sens[*] / exptime / dwdp 
		ENDELSE		 		
		red_wavel  = lam
		red_los    = los_name
		red_time   = time
		red_wavel  = (red_wavel+wshift)
		red_error  = red_emiss * 0.1
		IF diag EQ 'FVS' THEN instr_func = 0.08
		instr_func = 0.08
		PRINT,'Time resolution: ',red_time(1)-red_time(0)
		IF ~KEYWORD_SET(resolution) THEN resolution   = 0.0001
		act_resol    = red_time(1)-red_time(0)
		if keyword_set(interelm)then begin
			telm       = find_elm(shot,red_time,red_emiss(*,0,0))
			if ~keyword_set(elmcond)then elmcond=5
			idback     = where(telm ge elmcond)
			if idback[0] ne -1 then begin
				red_time   = red_time[idback]
				red_emiss  = red_emiss[*,*,idback]
				red_error  = red_error[*,*,idback]
			endif
		endif
	ENDELSE	
	rvals = -1
ENDIF
IF machine EQ 'JET'THEN BEGIN
        !PATH=!PATH + ':' + $
        expand_path( '+~cxs/utilities' ) + ':' + $
        expand_path( '+~cxs/calibration' ) + ':' + $
        expand_path( '+~cxs/instrument_data' )
;	 + ':' + $
;         expand_path( '+~cxs/idl/ks457_0/programs') 
;        expand_path( '+~cxs/ks6read/' ) + ':' + $
;        expand_path( '+~cxs/ktread/' ) + ':' + $
;        expand_path( '+~cxs/kx1read/' ) + ':' + $

    IF ~KEYWORD_SET(diag)THEN diag        = 'KT3B'
    IF ~KEYWORD_SET(wshift)THEN wshift      = 0.0
    instr_lorz  = 0.1

    IF KEYWORD_SET(calwave) THEN BEGIN
    	nomodel=1
    	calwave=STRING(shot,FORMAT='(I5)')
    ENDIF    

;**********************************
;**** Read data from read_spec ****
;**********************************

    shotstr = STRING(shot,FORMAT='(I5)')
;    if keyword_set(kt3a)then $
;	  data       = get_kt3a(shot,spec=spec,debug=debug,rrange=rrange,psplot=psplot,sig=sig) else $
;	  data       = get_kt3b(shot,spec=spec,debug=debug,rrange=rrange,psplot=psplot,sig=sig)
    data = get_kt3(shot,spec=diag,debug=debug,rrange=rrange,psplot=psplot,sig=sig)
    red_emiss  = data.phflx
    red_wavel  = data.lamgrid
    red_los    = data.losnames
    red_time   = data.time
    red_wavel  = (red_wavel+wshift)
    red_error  = data.phflx*0.01
    rvals      = data.rval
    IF strlowcase(diag) eq 'kt3b' then instr_func = 0.076 else instr_func = 0.21
ENDIF

data = {emiss:red_emiss,wvlngth:red_wavel,time:red_time,error:red_error,los:red_los,instr_func:instr_func,rvals:rvals}

Return, data
End

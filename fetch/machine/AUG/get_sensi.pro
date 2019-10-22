pro get_sensi, shot, diag, lam0, slit, blende, em_gain, ad_gain, $
               gratcons, losnames, npix, nchan, setup_num, camera_num, $
               verbose=verbose, no_calib=no_calib,$
               sens


pixel = findgen(npix)+1
sens = fltarr(npix,nchan)+1.
no_calib = 1

CASE diag of
    'DVS':begin
        blende_cal = 2.8
        if em_gain eq 1 then gain = ad_gain*0.35
        if em_gain eq 2 then gain = ad_gain
        if em_gain gt 2 then gain = ad_gain*em_gain*0.51
    end
    'EVS':begin
        blende_cal = 4.0
        if shot gt 27500 or shot lt 10000 then begin
            gain = 2.^(ad_gain-1)*em_gain
        endif else begin
            if em_gain eq 1 then gain = ad_gain*0.35
            if em_gain eq 2 then gain = ad_gain
            if em_gain gt 2 then gain = ad_gain*em_gain*0.51
        endelse
    end
    'FVS':begin
        blende_cal = 4.0
        if shot gt 26500 or shot lt 10000 then begin
            gain = 2.^(ad_gain-1)*em_gain
        endif else begin
            if em_gain eq 1 then gain = ad_gain*0.35
            if em_gain eq 2 then gain = ad_gain
            if em_gain gt 2 then gain = ad_gain*em_gain*0.51
        endelse
    end
    'HVS':begin
        blende_cal = 1.0
        gain = 2.^(ad_gain-1)*em_gain
    end
    'CFR':begin
        blende_cal = 2.8
        gain = ad_gain*em_gain
    end 
    'COR':begin
        blende_cal = 2.8
        gain = ad_gain*em_gain
    end 
    'BES':begin
        blende_cal = 10.0
        if shot gt 29000 and shot le 30015 then setup_num = 1L
        gain = ad_gain
    end
    'BEP':begin
        blende_cal = blende
        if em_gain eq 0 then begin
            if ad_gain eq 1 then gain = 0.487585
            if ad_gain eq 2 then gain = 1.0
            if ad_gain eq 3 then gain = 1.895
        endif else begin
            if ad_gain eq 1 then gain = em_gain*0.151
            if ad_gain eq 2 then gain = em_gain*0.255
            if ad_gain eq 3 then gain = em_gain*0.558
        endelse
    end
    'CDL':begin
        blende_cal = 1.
        gain = em_gain*ad_gain
    end
    'LVS':begin
        blende_cal = 2.8
        if em_gain eq 0 then begin
            if ad_gain eq 1 then gain = 1.0
            if ad_gain eq 2 then gain = 2.0
            if ad_gain eq 3 then gain = 4.0
        endif
        if em_gain eq 3000 and ad_gain eq 4 then gain = 19.02
        if em_gain eq 3500 and ad_gain eq 4 then gain = 99.02
    end             
    ELSE:begin
        print,' '
        print,' Diag '+diag+' not included in get_sensi!'
        print,' '
        return
    end
ENDCASE
gain = gain * (blende_cal/blende)^2.

if keyword_set(verbose) then begin
    print,' '
    print,' get_sensi: Blende: ', blende
    print,' get_sensi: AD-Gain:', ad_gain
    print,' get_sensi: EM-Gain:', em_gain
    print,' get_sensi: Gain:   ', gain
    print,' get_sensi: gratcons:   ', gratcons
endif 

if diag eq 'HVS' and lam0 eq 750. and em_gain ne 15 then begin
    if keyword_set(verbose) then $
      print,' get_sensi: reading private calibration!'
    path = '/afs/ipp/u/sfp/ausw/CALIB/LOS_SPEC_CALIB/SENS/'
    for ichan = 0,nchan-1 do begin
        filename = path+strtrim(losnames(ichan))+'_'+diag+'_'+'Chan'+$
                   string(ichan+1,f='(i2.2)')+'_'+string(lam0,f='(f6.2)')+$
                   'nm_'+string(gratcons,f='(i4.4)')+'Lmm.sav'
        exists = file_search(filename)
        if exists[0] ne '' then begin
            restore,filename
            sens[*,ichan] = sensitivity * (slit/slit_cal) * (gain/gain_cal)
            no_calib = 0
            if keyword_set(verbose) then $
              print,' get_sensi: Calib found for ',losnames[ichan]
        endif else begin
            if keyword_set(verbose) then $
              print,' get_sensi: No calib for ',losnames[ichan]
        endelse 
    endfor 
endif else begin       
    dum_sens = fltarr(npix)
    status = 0L
    ldiag = lonarr(3)
    for i = 0,2 do ldiag[i] = long(byte(strmid(diag,i,1)))
    if (!VERSION.MEMORY_BITS eq 32) then begin 
        lib = '/afs/ipp/u/sprd/CALIB/calcsens32.so'
    endif else begin   
        lib = '/afs/ipp/u/sprd/CALIB/calcsens64.so'
    endelse 
    for ichan = 0,nchan-1 do begin
        llosn = lonarr(8)
        n = strlen(losnames[ichan])
        for i = 0,n-1 do llosn[i] = long(byte(strmid(losnames[ichan],i,1)))
        for i = n,7 do llosn[i] = long(byte(' '))
        
        s = Call_External(lib,'calcsens_', $
                          long(shot), ldiag, gratcons, $
                          long(camera_num), long(setup_num), $
                          llosn, long(ichan+1), slit, gain, $
                          lam0, pixel, long(npix), dum_sens, status)
        if status ne 0 then begin
            print,' get_sensi: No calib for ',losnames[ichan] 
        endif else begin
            no_calib = 0
        endelse 
        sens[*,ichan] = dum_sens
    endfor
endelse



END

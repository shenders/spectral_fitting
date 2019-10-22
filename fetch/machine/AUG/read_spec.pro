pro read_spec,shot, diag, $
              exp=exp, path=path, error=error, $
              save_data=save_data, netcdf=netcdf, $
              verbose=verbose, $
              data

; needed subroutines
;@/afs/ipp/u/sfp/idlpro/input/magnification
;@/afs/ipp/u/sfp/idlpro/input/dispersion
;@/afs/ipp/u/sfp/idlpro/input/read_los_coord
;@/afs/ipp/u/sfp/idlpro/input/get_sensi
;@/afs/ipp/u/sfp/idlpro/input/write_netcdf
;@/afs/ipp/u/sfp/idlpro/input/read_netcdf


if n_params() eq 0 then begin
    print, ' '
    print, 'USAGE: '
    print, 'read_spec,shot, diag, '
    print, '          exp=exp, path=path, error=error, '
    print, '          save_data=save_data, netcdf=netcdf, '
    print, '          verbose=verbose, '
    print, '          data'
    print, ' '
    return
endif

libddww='/afs/ipp/aug/ads/lib64/@sys/libddww8.so' 

filename = getenv(diag+'_PATH')+'all_spectra_'+diag+'_'+string(shot,f='(i5.5)')
if keyword_set(netcdf) then filename += '.ncdf' else filename += '.sav'
if keyword_set(path) then filename = path+filename
exist = file_search(filename)

IF exist eq '' THEN BEGIN
    
  if not keyword_set(exp) then exp='AUGD'
  error     = 0L
  shot      = long(shot)
  edt       = 0L
  dia_ref   = 0L
  date      = 'YY:MMM:DD;hh:mm:ss'
  phys_unit = 0L
  dim       = 0L
  stop      = 0
  
;; open shotfile
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddopen', $
                    error, exp, diag, shot, edt, dia_ref, date)
  
  if error ne 0 then begin
      message = exp+' '+diag+' '+string(shot,f='(i5.5)')+$
                '('+string(edt,f='(i2)')+')'
      ddww_error,error,message
      return
  endif
  if keyword_set(verbose) then $
    print,' Date of Edititon: ', date


;; ---------------------------------------------------------------------------------------
;:  READ SPECTROMETER SETTINGS
;; ---------------------------------------------------------------------------------------
  pname = 'PARAM' 
  if diag eq 'LVS' then begin
      sfh_type = 0L
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, pname, 'SFH_TYPE', $
                        1L, 1L, sfh_type, phys_unit)
      if sfh_type ne 5 then begin
          print,' '
          print,' ERROR: '
          print,' Old LVS shotfile header!'
          print,' '
          stop = 1
          goto,close_shotfile
      endif
  endif 
  slit = 0.
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'SLIT', $
                    2L, 1L, slit, phys_unit)
  
; repair wrong slit in L0 shotfile
  if diag eq 'DVS' then begin
      if shot le 24268 and shot ge 24189 then slit = 150.
  endif 
  if diag eq 'EVS' then begin
      if shot ge 23144 and shot le 23153 then slit = 200.
      if shot ge 24801 and shot le 24815 then slit = 150.
      if shot eq 26346 or (shot ge 23873 and shot le 23927) or $
        (shot ge 29852 and shot le 29973) then slit = 100.
  endif
  if diag eq 'LVS' then begin
      if shot ge 26156 and shot le 26395 then slit = 200.
  endif
  if diag eq 'CFR' then begin
      if shot eq 28731 or shot eq 28732 or shot eq 29731 then slit = 150.
  endif 
  if diag eq 'BES' then begin
      if shot ge 28525 and shot le 29390 then slit = 200.
      if shot ge 29879 and shot le 29892 then slit = 60.
  endif

  op_ang = 0.
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'OP_ANG', $
                    2L, 1L, op_ang, phys_unit)

  foc_len = 0.
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'FOC_LEN', $
                    2L, 1L, foc_len, phys_unit)
  if error ne 0 then begin
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, pname, 'FOC_LEN1', $
                        2L, 1L, foc_len, phys_unit)
      if error ne 0 then ddww_error, error, 'FOC_LEN'
  endif
  if diag eq 'DVS' and shot gt 30150 then foc_len = 0.175165
  if diag eq 'BES' and shot gt 30150 then foc_len = 1.002402

  pixw = 13.
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'PIXW', $
                    2L, 1L, pixw, phys_unit)
  if diag eq 'DVS' and shot le 30150 then pixw = 7.735
  if diag eq 'HVS' and pixw eq 0 then pixw = 13
  if diag eq 'BES' then pixh = pixw*200./150. else pixh = pixw

  em_gain = 1L
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'EM-GAIN', $
                    1L, 1L, em_gain, phys_unit)

  ad_gain = 1L
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'AD-GAIN', $
                    1L, 1L, ad_gain, phys_unit)

  dyn_range = 16L
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'DynRange', $
                    1L, 1L, dyn_range, phys_unit)

  delta_h = 0.
  if shot gt 99999 then $
    s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                      error, dia_ref, pname, 'DELTA_H', $
                      2L, 1L, delta_h, phys_unit)
  if diag eq 'BES' and (shot lt 15000 or shot gt 27401) then $
    delta_h = -10.8

  blende = 1.
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'BLENDE_1', $
                    2L, 1L, blende, phys_unit)
  
; repair wrong blende in L0 shotfile 
  if diag eq 'EVS' then begin
      if shot ge 23144 and shot le 23153 then blende = 4.0
  endif
  if diag eq 'FVS' then begin
      if shot ge 22825 and shot le 22860 then blende = 8.0
      if shot ge 31470 and shot le 31606 then blende = 11.0
  endif

  ctsph = 0.
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'CTS_PHOT', $
                    2L, 1L, ctsph, phys_unit)
  sad2 = 0.
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'SDEV_AD', $
                    2L, 1L, sad2, phys_unit)
  if diag eq 'LVS' and em_gain eq 0 then begin
      ctsph = 0.32
      if ad_gain eq 1 then sad2 = 5.6
      if ad_gain eq 2 then sad2 = 9.9
      if ad_gain eq 3 then sad2 = 13.9      
  endif
  if diag eq 'CDL' then begin
      sad2 = 7.
      ctsph = 0.7
  endif 

  mode = 2
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'MODE', $
                    11L, 1L, mode, phys_unit)
  
  if diag eq 'GVS' or diag eq 'HVS' then begin
      pname = 'CCD'
      sname = 'SpecGroo'
  endif else begin
      pname = 'PARAM'
      sname = 'GRATCONS'
  endelse 
  gratcons = 0.
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, sname, $
                    2L, 1L, gratcons, phys_unit)

  pname = 'WL_MOTOR'
  sname = 'WAVE'
  if diag eq 'GVS' or diag eq 'HVS' then begin
      pname = 'CCD'
      sname = 'SpecCent'
  endif 
  if diag eq 'CDL' or diag eq 'LVS' then begin
      pname = 'PARAM'
      sname = 'WAVELEN'
  endif 
  lam0 = 0.
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, sname, $
                    2L, 1L, lam0, phys_unit)
  if error ne 0 then begin
      ddww_error,error,' WAVELEN'
      stop = 1
      goto,close_shotfile
  endif 
  

;; ---------------------------------------------------------------------------------------
;;  READ # of TRIGGERS, EXPTIME and TIMEBASE
;; ---------------------------------------------------------------------------------------
  if diag eq 'GVS' or diag eq 'HVS' or diag eq 'CDL' then tbname = 'T-Base' $
  else tbname = 'TIME'
  
  tbeg = 0.0
  tend = 0.0
  ntrig = 0L
  npretrig = 0L
  error_tb = 0L
  s = CALL_EXTERNAL(libddww,'ddgetaug','ddtrange', $
                    error_tb, dia_ref, tbname, tbeg, tend, ntrig, npretrig)

  if error_tb ne 0 then begin
      ddww_error,error_tb,'ddtrange'
      print,' #triggers according to PPG:',ntrig
  endif else begin  
      time = fltarr(ntrig)
      ndataread = 0L
      s = CALL_EXTERNAL(libddww,'ddgetaug','ddtbase',$
                        error, dia_ref, tbname, 1L, ntrig, $
                        2L, ntrig, time, ndataread)
  endelse  
      
  waittime = 0.
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, 'PARAM', 'WAITTIME', $
                    2L, 1L, waittime, phys_unit)
    
; read CAMAC PPG device
  if diag eq 'GVS' or diag eq 'HVS' or diag eq 'CDL' or diag eq 'LVS' then begin
      pulsa = intarr(16)
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, 'PPG', 'PULSES', $
                        11L, 16L,pulsa, phys_unit)
      if error ne 0 then begin
          ddww_error,error,'PPG PULSES'
          goto, close_shotfile
      endif 
      
      reptime = intarr(16)
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, 'PPG', 'RESOLUT', $
                        11L, 16L, reptime, phys_unit)
      if error ne 0 then begin
          ddww_error,error,'PPG RESOLUT'
          goto, close_shotfile
      endif 
      
      reptfac = intarr(16)
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, 'PPG', 'RESFACT', $
                        11L, 16L, reptfac, phys_unit)
      if error ne 0 then begin
          ddww_error,error,'PPG RESFACT'
          goto, close_shotfile
      endif 
      exptime = reptime[0]*10.^(reptfac[0]-6)    ;[s]

      if pulsa[1] gt 0 then begin
          bgr_fra = 1
          n_offset = pulsa[1]-3
          offset_index = indgen(n_offset)+pulsa[0]+3
      endif else begin
          bgr_fra = 0
          n_offset = 25
          offset_index = indgen(n_offset)+pulsa[0]-n_offset
      endelse 
      
      
; read UDC-PPG device
  endif else begin
      pulsa = lonarr(7)     
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, 'PPG', 'PULS', $
                        1L, 7L,pulsa, phys_unit)
      if error ne 0 then begin
          ddww_error,error,'PPG PULS'
          goto, close_shotfile
      endif
      
      reptime = lonarr(7)
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, 'PPG', 'FREQ', $
                        1L, 7L, reptime, phys_unit)
      if error ne 0 then begin
          ddww_error,error,'PPG RESOLUT'
          goto, close_shotfile
      endif
      
      reptfac = intarr(7)
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, 'PPG', 'RES', $
                        11L, 7L, reptfac, phys_unit)
      if error ne 0 then begin
          ddww_error,error,'PPG RES'
          goto, close_shotfile
      endif
      exptime = reptime[0]*10.^(reptfac[0]*3-9)  ;[s]      

      output = intarr(7)
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, 'PPG', 'OUT', $
                        11L, 7L, output, phys_unit)
      if error ne 0 then begin
          ddww_error,error,'PPG OUT'
          goto, close_shotfile
      endif
      output_index = where(output eq 1,noutidx)
      ww = where(pulsa[output_index] gt 0,noutput)
      if noutput gt 1 then begin
          bgr_fra = 1
          n_offset = pulsa[output_index[1]]-3
          offset_index = indgen(n_offset)+pulsa[0]+3
      endif else begin
          bgr_fra = 0
          n_offset = 25
          offset_index = indgen(n_offset)+pulsa[0]-n_offset
      endelse 
  endelse 
     
; Frame Transfer Mode
  if mode eq 2 then begin
      if keyword_set(verbose) then begin
          print, ' '
          print,' Frame Transfer Mode'
      endif 
      ntime = pulsa[0]
      ; built time vector manually if there was a problem with the timebase
      if error_tb ne 0 then begin
          outind = where(output eq 1,nout)
          ntrig = total(pulsa[outind])
          time = findgen(ntime)*exptime-0.5*exptime+waittime
          if keyword_set(verbose) then print, ' Built time vector manually!'
      endif else begin
          time = time-0.5*exptime+waittime
      endelse 
      if ntime ne ntrig and keyword_set(verbose)then begin
          print,' '
          print,' WARNING: '
          print,'    ntime             ',ntime
          print,'    not equal ntrig:',ntrig
          print,' '
          time = temporary(time[0:ntime-1])
      endif
  endif 
  
  ; TDI Mode
  if mode eq 1 then begin
      if keyword_set(verbose) then begin
          print,' '
          print,' TDI Mode'
      endif 
      exptime = 0.
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, 'PARAM', 'EXPTIME', $
                        2L, 1L, exptime, phys_unit)
      ntime = ntrig
      npuls = pulsa[0]
      mean_dtm = 0.
      n_mean = 0.

      for itime = 2L,ntime-1 do begin
          dtm = time[itime]-time[itime-1]
          if abs(dtm-exptime)/exptime lt 0.1 then begin
              mean_dtm = mean_dtm+dtm
              n_mean = n_mean+1
          endif 
      endfor
      exptime = mean_dtm/n_mean
      dt_clear = 5.3e-6
      dt_ppgw = 2.0e-6
      nt = 0L
      for ipuls = 0,npuls-1 do begin
          for i = 0L,(ntime/npuls)-1 do begin
              time[nt] = ipuls*reptime[0]+(float(i+1)-0.5)*exptime+dt_clear+dt_ppgw
              nt += 1L
          endfor
      endfor
      bgr_fra = 0
      n_offset = 20
      offset_index = indgen(n_offset)+20
  endif 
  
;; ---------------------------------------------------------------------------------------
;;  READ CCD DATA
;; ---------------------------------------------------------------------------------------
  if diag eq 'GVS' or diag eq 'HVS' or diag eq 'CDL' then begin
      pname = 'CCD'
      sname = 'CCDSig'
  endif else begin
      pname = 'READ_CCD'
      sname = 'CCD_DATA'
  endelse

  nchan = 0
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'ydim', $
                    11L, 1L, nchan, phys_unit)
  if error ne 0 then begin
      ddww_error,error,pname+' ydim'
      goto, close_shotfile
  endif 
  
  xdim = 0
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'xdim', $
                    11L, 1L, xdim, phys_unit)

  startx = 0
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'startx', $
                    11L, 1L, startx, phys_unit)
  if diag eq 'CDL' and shot lt 32460 then startx += 1

  npix = xdim +startx-1
  if error ne 0 then begin
      ddww_error,error,pname+' xdim'
      goto, close_shotfile
  endif 
  
  length  = long(xdim*nchan)
  rawdata = lonarr(length,ntrig)
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddsgroup', $
                    error, dia_ref, sname, 1L, length, 1L, length, $
                    rawdata, dim)
  if error ne 0 then begin
      ddww_error,error,sname
      goto, close_shotfile
  endif

  counts = float(reform(rawdata,xdim,nchan,ntrig,/overwrite))
  if startx ne 1 then counts = temporary(counts[abs(startx)+1:xdim-1,*,*])

; check for saturation
  saturated = where(counts ge (2.^dyn_range-50.), nsaturated)

  npixy = 0
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'yDimDet', $
                    11L, 1L, npixy, phys_unit)
  if error ne 0 then begin
      ddww_error,error,pname+' yDimDet'
      goto, close_shotfile
  endif

  starty = intarr(nchan)
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'starty', $
                    11L, long(nchan), starty, phys_unit)
  if error ne 0 then begin
      if diag eq 'HVS' and shot le 31776 then begin
          starty = [  1,  65, 129, 193, 257, 321, 385, 449, $
                    513, 577, 641, 705, 769, 833, 897, 961]
      endif else begin
          ddww_error,error,pname+' starty'
      endelse
  endif 
  
  groupy = intarr(nchan)+1
  s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                    error, dia_ref, pname, 'groupy', $
                    11L, long(nchan), groupy, phys_unit)
  if error ne 0 then begin
      if diag eq 'HVS' and shot le 31776 then begin
          groupy = intarr(nchan)+64
      endif else begin
          ddww_error,error,pname+' groupy'
      endelse
  endif 

;; ---------------------------------------------------------------------------------------
;;  READ LOSNAMES AND COORDINATES
;; ---------------------------------------------------------------------------------------
  pname    = 'PARAM'
  losdum   = '12345678'
  losnames = strarr(nchan)
  
  for ichan = 0, nchan-1 do begin
      dumstr = 'CHAN_'+string(ichan+1,f='(I2.2)')
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, pname, dumstr, $
                        6L,8L, losdum, phys_unit)
      if error ne 0 then begin
          ddww_error,error,' reading losnames '
          goto, close_shotfile
      endif
      losnames[ichan] = string(byte(losdum))
  endfor
  ; repair wrong losnames in L0 file
  if shot eq 31312 and diag eq 'EVS' then losnames([2,3]) = losnames([3,2])
  if shot eq 32201 and diag eq 'EVS' then losnames([0])  = 'RH2-1   '
  if shot eq 32201 and diag eq 'HVS' then losnames([13]) = 'ROV007  '
  
  read_los_coord,shot,losnames,verbose=verbose,$
                 r1,z1,f1,r2,z2,f2

  
;; ---------------------------------------------------------------------------------------
;;  READ WAVELENGTH POLYNOM COEFFICIENTS of INTERNAL
;;  SPECTROMETER CALIBRATION for GVS and HVS
;; ---------------------------------------------------------------------------------------
  if diag eq 'GVS' or diag eq 'HVS' then begin
      pname = 'CCD'
      p_order = 'a'
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, pname, 'p_order', $
                        6L, 1L, p_order, phys_unit)
      p_order = strtrim(p_order)
      
      p_coeff = dblarr(6)
      s = CALL_EXTERNAL(libddww, 'ddgetaug', 'ddparm', $
                        error, dia_ref, pname, 'p_coeff', $
                        3L, 6L, p_coeff, phys_unit)
      
     ; repair wrong dispersion in L0 file 
      if diag eq 'HVS' and shot eq 30624 then $
        p_coeff = [450.84816d, 0.019872289d, -1.2454068d-07]
  endif
  
;; ---------------------------------------------------------------------------------------
;;  CLOSE SHOTFILE
;; ---------------------------------------------------------------------------------------
close_shotfile:
  error_2 = 0L
  s = CALL_EXTERNAL(libddww,'ddgetaug','ddclose', $
                    error_2, dia_ref)
  if error_2 ne 0 then begin
      ddww_error,error_2,exp+' '+diag+' '+string(shot,f='(i5.5)')
      return
  endif
  if stop eq 1 then stop

;; ---------------------------------------------------------------------------------------
;;  WAVELENGTH CALIBRATION
;; ---------------------------------------------------------------------------------------
  if diag eq 'GVS' or diag eq 'HVS' then begin
      dwdp   = p_coeff[1]
      d2wdp2 = p_coeff[2]
      
      xx  = dindgen(npix)
      wavel = dblarr(npix)
      for i = 0,long(p_order) do wavel[*] += p_coeff[i]*xx[*]^(i)
  endif else begin
      dispersion, gratcons, op_ang, foc_len, lam0, pixw, $
                  dwdp,d2wdp2

      pixarr = findgen(npix)-((npix-1.)/2.)
      wavel = lam0 + dwdp*pixarr + 0.5*d2wdp2*pixarr^2.
  endelse 

  
; Calculate parabula shift
  if total(groupy) eq nchan then begin
      d_z = groupy*0.
      wvl_shiftv = d_z
  endif else begin
      d_z = ((2.*starty+groupy-1-npixy)*0.5-delta_h)*pixh
      wvl_shiftv = 0.5*lam0*(d_z/foc_len/1e6)^2.
  endelse 


; determine magnification and slit, dfibre, fwhm in nm
  if diag eq 'EVS' or diag eq 'FVS' then $
    f2_f1 = foc_len/0.28 else f2_f1 = 1.
  fwhm_gauss = 10.
  if diag eq 'FVS' and shot ge 25903 then fwhm_gauss = 38.4 ;55.6
  if diag eq 'EVS' and shot ge 29852 then fwhm_gauss = 46.  ;35.6
  if diag eq 'CDL' then fwhm_gauss = 51.
  if diag eq 'DVS' then fwhm_gauss = 23.
  if diag eq 'CFR' then fwhm_gauss = 41.7
 
  mag = magnification(gratcons,op_ang,lam0,f2_f1)
  slit_nm = slit*mag*dwdp/pixw
  fibre_nm = 400.*mag*dwdp/pixw
  fwhm_nm = fwhm_gauss*dwdp/pixw

  if keyword_set(verbose) then begin
      print,' exptime [ms]:            ', exptime/1e-3
      print,' # of frames:             ', ntime
      print,' # of channels:           ', nchan
      print,' # of pixels:             ', npix
      print,' Lambda at central pixel: '+string(lam0,f='(f7.3)')
      print,' Delta Lambda per Pixel:  '+string(dwdp,f='(f6.4)')
      print,' Slitwidth [mu]:          '+string(slit,f='(f5.1)')
      print,' '
      print,' Delta_h [mu]:',delta_h
      print,' Vertical shift [nm]:'
      for ichan=0,nchan-1 do print,ichan+1,d_z[ichan],wvl_shiftv[ichan]
  endif

  lamgrid = fltarr(npix,nchan)
  for ichan = 0,nchan-1 do lamgrid[*,ichan] = wavel-wvl_shiftv[ichan]

;; ---------------------------------------------------------------------------------------
;;  REMOVE OFFSET
;; ---------------------------------------------------------------------------------------
  offset_error = 0
  if max(counts[*,*,offset_index]) ge 2000. or $
    max(counts[*,*,offset_index]) eq 0. then begin
      print,' '
      print,' ERROR in offset determination! '
      print,' Intensity of last 25 frames either too high or 0.!'
      print,' Max/Mean offset [counts]: '+string(max(counts[*,*,offset_index]),f='(i5)')+$
            '/'+string(mean(counts[*,*,offset_index]),f='(i5)')
      print,' Try first 25 frames instead!'
      n_offset = 25
      offset_index = indgen(n_offset)+20
      if max(counts[*,*,offset_index]) ge 2000. or $
        max(counts[*,*,offset_index]) eq 0. then begin
          print,' '
          print,' Also intensity of first 25 frames either too high or 0.!'
          print,' Max/Mean offset [counts]: '+$
                string(max(counts[*,*,offset_index]),f='(i5)')+'/'+$
                string(mean(counts[*,*,offset_index]),f='(i5)')
          print,' No Background substraction applied!'
          print,' '
          offset_error = 1
      endif 
  endif
  
  offset = fltarr(npix,nchan)
  if offset_error eq 0 then begin
      for ichan = 0,nchan-1 do $
        offset[*,ichan] = total(counts[*,ichan,offset_index],3)/n_offset
      for itime = 0L,ntime-1 do counts[*,*,itime] = (counts[*,*,itime]-offset[*,*]) > 0.
      if nsaturated gt 0 then counts[saturated] = !Values.f_NaN
  endif

  if diag eq 'FVS' and (shot eq 31178 or shot eq 31180) then begin
      if keyword_set(verbose) then print,' Apply dedicated Background substraction!'
      offset = fltarr(npix,nchan)
      for ichan = 0,nchan-1 do begin
          offset[*,ichan] = total(counts[*,ichan,offset_index],3)/n_offset
          offset[253:290,ichan] = 0.
          wo = where(offset[*,ichan] ne 0)
          ydum = reform(offset[wo,ichan])
          xdum = reform(lamgrid[wo,ichan])
          offset[*,ichan] = interpol(ydum,xdum,lamgrid[*,ichan])
      endfor
      for itime = 0L,ntime-1 do counts[*,*,itime] = (counts[*,*,itime]-offset[*,*]) > 0.
      if nsaturated gt 0 then counts[saturated] = !Values.f_NaN
  endif 
                  
  if keyword_set(verbose) then begin
      print,' '
      print,' Dynamic range:        ',dyn_range
      print,' # of saturated frames:',nsaturated
      if bgr_fra eq 1 then print,' Dedicated Background Frames found!'
      if offset_error eq 0 then  $
        print,' Max/Mean offset [counts]: '+string(max(offset),f='(i5)')+$
              '/'+string(mean(offset),f='(i5)')
  endif
  
  if ntime ne ntrig and mode eq 2 then begin
      if keyword_set(verbose) then print,' Reform counts from ntrig->ntime'
      counts = temporary(counts[*,*,0:ntime-1])
  endif
  

  
;; ---------------------------------------------------------------------------------------
;;  REMOVE SMEAR EFFECT
;;
;;  ---------------------------------------------------------------------------------------
  bgr_index = where(losnames eq 'BACKGRND',nwo)
  if nwo gt 0 and total(groupy) ne nchan and diag ne 'GVS' then begin
      bgr_index = bgr_index[0]
      smear = reform(counts[*,bgr_index,*]); > 0.

      for ichan = 0,nchan-1 do $
        counts[*,ichan,*] -= smear[*,*]/groupy[bgr_index]*groupy[ichan]

      not_saturated = where(finite(counts))
      counts[not_saturated] = floor(counts[not_saturated]) > 0.
      
      counts[*,bgr_index,*] = smear[*,*]
      losnames(bgr_index) = 'SMEAR'
      if keyword_set(verbose) then print,' Remove smear'
  endif
  wo = where(losnames ne 'SMEAR',nwo)


  
;; ---------------------------------------------------------------------------------------
;;  APPLY INTENSITY CALIBRATION
;; ---------------------------------------------------------------------------------------
  camera_num = 0L
  setup_num = 0L
  get_sensi, shot, diag, lam0, slit, blende, em_gain, ad_gain, $
             gratcons, losnames, npix, nchan, setup_num, camera_num, $
             verbose=verbose, no_calib=no_calib,$
             sens

  phflx = counts
  sy    = counts
  if no_calib eq 0 then begin
      wo = where(losnames ne 'SMEAR')
      for itime = 0L,ntime-1 do begin  
          phflx[*,wo,itime] = counts[*,wo,itime]/sens[*,wo]/exptime/dwdp
          sy[*,wo,itime] = sqrt(counts[*,wo,itime]*ctsph+sad2)/sens[*,wo]/exptime/dwdp
      endfor
  endif else begin
      wo = where(losnames ne 'SMEAR')
      for itime = 0L,ntime-1 do begin  
          sy[*,wo,itime] = sqrt(counts[*,wo,itime]*ctsph+sad2);/exptime/dwdp
      endfor
  endelse 
      
  
  ; DEFINE OUTPUT STRUCTURE
  data = {time:     time, $
          rawwvl:   wavel,$
          counts:   counts, $
          lamgrid:  lamgrid, $
          phflx:    phflx, $
          sy:       sy, $
          sens:     sens, $
          offset:   offset, $
          exptime:  exptime, $
          lam0:     lam0, $
          dwdp:     dwdp, $
          d2wdp2:   d2wdp2, $
          pixw:     pixw, $
          slit:     slit, $
          slit_nm:  slit_nm, $
          gratcons: gratcons, $
          op_ang:   op_ang, $
          foc_len:  foc_len, $
          mag:      mag, $
          fibre_nm: fibre_nm, $
          fwhm_nm:  fwhm_nm, $
          no_calib: no_calib, $
          losnames: losnames, $
          r1:       r1, $
          z1:       z1, $
          f1:       f1, $
          r2:       r2, $
          z2:       z2, $
          f2:       f2, $
          sad2:     sad2, $
          ctsph:    ctsph, $
          em_gain:  em_gain, $
          ad_gain:  ad_gain, $
          shot:     shot, $
          diag:     diag $
         } 
         
  ; write netCDF file
  if keyword_set(save_data) then begin
      if keyword_set(netcdf) then begin
          write_netcdf,data,filename,error_netcdf
          if error_netcdf ne 0 then print,' ERROR in write netCDF file!'
      endif else begin
          save, data, file=filename
      endelse 
      if keyword_set(verbose) then begin
          print,' '
          print,' Data saved to: '+filename
          print,' '
      endif
  endif 
  
ENDIF ELSE BEGIN

  if keyword_set(verbose) then begin
      print,' '
      print,' Read Data from: '+filename
      print,' '
  endif
  error = 0
  ; read shot and check the type of netcdf
  if keyword_set(netcdf) then begin
      read_netcdf,filename,'shot',shot,/gatr
      if shot gt 33725 then begin
          read_netcdf,filename,'spectra',counts
          read_netcdf,filename,'offset',offset
          read_netcdf,filename,'neon',neon
          read_netcdf,filename,'lambda',lamgrid
          read_netcdf,filename,'sens',sens
          read_netcdf,filename,'time',time
          read_netcdf,filename,'R_LOS_1',r1
          read_netcdf,filename,'Z_LOS_1',z1
          read_netcdf,filename,'PHI_LOS_1',f1
          read_netcdf,filename,'R_LOS_2',r2
          read_netcdf,filename,'Z_LOS_2',z2
          read_netcdf,filename,'PHI_LOS_2',f2
          read_netcdf,filename,'fwhm_pix',fwhm_pix
          read_netcdf,filename,'losnam',losnames
          read_netcdf,filename,'exptime',exptime,/gatr
          read_netcdf,filename,'gratcons',gratcons,/gatr
          read_netcdf,filename,'op_ang',op_ang,/gatr
          read_netcdf,filename,'foc_len',foc_len,/gatr
          read_netcdf,filename,'dwdp',dwdp,/gatr
          read_netcdf,filename,'magnification',mag,/gatr
          read_netcdf,filename,'sad2',sad2,/gatr
          read_netcdf,filename,'ctsph',ctsph,/gatr
          read_netcdf,filename,'wlen',lam0,/gatr
          read_netcdf,filename,'slit',slit,/gatr
          slit_nm = slit*1e-3
          read_netcdf,filename,'pixw',pixw,/gatr
          read_netcdf,filename,'lambda_neon',lambda_neon,/gatr
          read_netcdf,filename,'neon_done',neon_done,/gatr
          fwhm_nm = fwhm_pix*dwdp
          fibre_nm = 400.*mag*dwdp/pixw
          ; Calculate calibrated data and errobar
          phflx = counts
          sy    = counts
          for itime = 0L,n_elements(time)-1 do begin  
              phflx[*,itime,*] = counts[*,itime,*]/sens[*,*]/exptime/dwdp
              sy[*,itime,*] = sqrt(counts[*,itime,*]*ctsph+sad2)/sens[*,*]/exptime/dwdp
          endfor
          ; DEFINE OUTPUT STRUCTURE
          data = {time:     time, $
                  lamgrid:  lamgrid, $
                  counts:   counts, $
                  phflx:    phflx, $
                  sy:       sy, $
                  sens:     sens, $
                  offset:   offset, $
                  exptime:  exptime, $
                  lam0:     lam0, $
                  dwdp:     dwdp, $
                  pixw:     pixw, $
                  slit:     slit, $
                  slit_nm:  slit_nm, $
                  gratcons: gratcons, $
                  op_ang:   op_ang, $
                  foc_len:  foc_len, $
                  mag:      mag, $
                  fwhm_nm:  fwhm_nm, $
                  losnames: losnames, $
                  fibre_nm: fibre_nm,$
                  r1:       r1, $
                  z1:       z1, $
                  f1:       f1, $
                  r2:       r2, $
                  z2:       z2, $
                  f2:       f2, $
                  sad2:     sad2, $
                  ctsph:    ctsph, $
                  shot:     shot $
                 }
      endif else begin
          read_netcdf,filename,'TIME',time
          read_netcdf,filename,'LAMGRID',lamgrid
          read_netcdf,filename,'COUNTS',counts
          read_netcdf,filename,'PHFLX',phflx
          read_netcdf,filename,'SY',sy
          read_netcdf,filename,'SENS',sens
          read_netcdf,filename,'OFFSET',offset
          read_netcdf,filename,'EXPTIME',exptime
          read_netcdf,filename,'LAM0',lam0
          read_netcdf,filename,'DWDP',dwdp
          read_netcdf,filename,'D2WDP2',d2wdp2
          read_netcdf,filename,'PIXW',pixw
          read_netcdf,filename,'SLIT',slit
          read_netcdf,filename,'SLIT_NM',slit_nm
          read_netcdf,filename,'GRATCONS',gratcons
          read_netcdf,filename,'OP_ANG',op_ang
          read_netcdf,filename,'FOC_LEN',foc_len
          read_netcdf,filename,'MAG',mag
          read_netcdf,filename,'FIBRE_NM',fibre_nm
          read_netcdf,filename,'FWHM_NM',fwhm_nm
          read_netcdf,filename,'NO_CALIB',no_calib
          read_netcdf,filename,'LOSNAMES',losnames
          read_netcdf,filename,'R1',r1
          read_netcdf,filename,'Z1',z1
          read_netcdf,filename,'F1',f1
          read_netcdf,filename,'R2',r2
          read_netcdf,filename,'Z2',z2
          read_netcdf,filename,'F2',f2
          read_netcdf,filename,'SAD2',sad2
          read_netcdf,filename,'CTSPH',ctsph
          read_netcdf,filename,'DIAG',diag
          ; DEFINE OUTPUT STRUCTURE
          data = {time:     time, $
                  lamgrid:  lamgrid, $
                  counts:   counts, $
                  phflx:    phflx, $
                  sy:       sy, $
                  sens:     sens, $
                  offset:   offset, $
                  exptime:  exptime, $
                  lam0:     lam0, $
                  dwdp:     dwdp, $
                  d2wdp2:   d2wdp2, $
                  pixw:     pixw, $
                  slit:     slit, $
                  slit_nm:  slit_nm, $
                  gratcons: gratcons, $
                  op_ang:   op_ang, $
                  foc_len:  foc_len, $
                  mag:      mag, $
                  fibre_nm: fibre_nm, $
                  fwhm_nm:  fwhm_nm, $
                  no_calib: no_calib, $
                  losnames: losnames, $
                  r1:       r1, $
                  z1:       z1, $
                  f1:       f1, $
                  r2:       r2, $
                  z2:       z2, $
                  f2:       f2, $
                  sad2:     sad2, $
                  ctsph:    ctsph, $
                  shot:     shot, $
                  diag:     diag $
                 }
      endelse
  endif else begin
      restore,filename,verbose=verbose
  endelse
  
ENDELSE

end_pro:       
END

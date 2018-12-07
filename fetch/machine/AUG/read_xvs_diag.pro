@/afs/ipp/u/rld/idlpro/input/read_netcdf.pro
@/afs/ipp/u/sprd/idlpro/fit_spec/shot_to_year.pro

pro  read_xvs_diag, shot, diag, $
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

diag = strupcase(diag)
error = 0L
ishot = 0 
new_read = 0
lshot = long(shot)
year = shot_to_year(shot)
add = ''

if not keyword_set(no_copy) then begin
   file = '/tmp/all_spectra_'+diag+'_'+string(shot,f='(i5.5)')
   ;file = 'all_spectra_'+diag
   exists = file_search(file)
   if exists eq '' then new_read=1
endif else begin
   file = '/tmp/all_spectra_'+diag
   exists = file_search(file)
   if exists eq '' then begin 
	new_read=1
   endif else begin
        read_netcdf, file, 'shot', ishot, /gatr
        ff = abs(float(ishot-shot))
        if ff gt 0.5 then new_read = 1
   endelse
endelse

if keyword_set(read_again) then new_read = 1

if new_read eq 1 then begin
    expmt = 'AUGD'
    if keyword_set(exp_name) then expmt = exp_name
    if diag eq 'DVS' then begin
       if year ge 2017 then add='_2017' 
       line = '/afs/ipp/u/sprd/DVL/dvsget'+add+' '+' '''+expmt+''' '+string(shot,f='(i5.5)')
    endif
    if diag eq 'EVS' then begin
       if year ge 2012 then add='_2017' 
       line = '/afs/ipp/u/sprd/EVL/evsget'+add+' '+' '''+expmt+''' '+string(shot,f='(i5.5)')
    endif
    if diag eq 'FVS' then begin
       if year ge 2012 then add='_2017' 
       line = '/afs/ipp/u/sprd/FVL/fvsget'+add+' '+' '''+expmt+''' '+string(shot,f='(i5.5)')
    endif
    if diag eq 'GVS' then begin
       line = 'cd /tmp; /afs/ipp/u/sprd/GVL/gvsget'+add+' '+' '''+expmt+''' '+string(shot,f='(i5.5)')
    endif   
    if diag eq 'LVS' then begin
       if year ge 2014 then add='_2014'
       line = '/afs/ipp/u/sprd/LVL/lvsget'+add+' '+' '''+expmt+''' '+string(shot,f='(i5.5)')
    endif
    if diag eq 'BES' then begin
       if year ge 2014 then add='_2014'
       line = '/afs/ipp/u/cxrs/BEL/besget'+add+' '+' '''+expmt+''' '+string(shot,f='(i5.5)')
     endif
    if diag eq 'BEP' then begin
       add=' '
       line = '/afs/ipp/u/cxrs/BPZ/bepget'+add+' '+' '''+expmt+''' '+string(shot,f='(i5.5)')
    endif
    if diag eq 'CDL' then begin
       if year ge 2017 then add='_2017'
       line = 'cd /tmp; /afs/ipp/u/disp/CDM/cdlget'+add+' '+' '''+expmt+''' '+string(shot,f='(i5.5)')
    endif
    if diag eq 'CFR' then $
       line = 'cd /tmp; /afs/ipp/u/cxrs/idl/CFR_ausw/cfrget '+' '''+expmt+''' '+string(shot,f='(i5.5)')

    if keyword_set(no_smear_cor) then line = line +' 0 ' else line = line + ' '
    
    spawn, line

    if not keyword_set(no_copy) then begin
        line = string('mv /tmp/all_spectra_'+diag+' /tmp/all_spectra_'+diag+'_'+$
                     string(shot,f='(i5.5)'))
       spawn, line
    endif

    exists = file_search(file)
    if exists eq '' then begin 
       error=1L
       message, file+' does not exist!',/info
       return
    endif else begin
       read_netcdf, file, 'shot', ishot, /gatr
       ff = abs(float(ishot-shot))
       if ff gt 0.5 then begin 
          message,'no '+file+' with  shot='+string(shot,f='(i5.5)')+ 'found!',/info
          message,expmt+'/'+diag+'/'+string(shot,f='(i5.5)')+' does not exist!',/info
          error = 1L
          return
       endif
    endelse
endif

read_netcdf, file, 'dwdp', dwdp, /gatr
read_netcdf, file, 'exptime', exptime, /gatr
read_netcdf, file, 'sad2', sad2, /gatr
read_netcdf, file, 'ctsph', ctsph, /gatr
read_netcdf, file, 'slit', wslit, /gatr
read_netcdf, file, 'wlen', wlen, /gatr
read_netcdf, file, 'pixw', pixw, /gatr
read_netcdf, file, 'foc_len', foc_len, /gatr

read_netcdf, file, 'time',  time, /all
read_netcdf, file, 'lambda',  lam, /all
read_netcdf, file, 'offset',  offset, /all
read_netcdf, file, 'sens',  sens, /all
read_netcdf, file, 'spectra', spec, /all
read_netcdf, file, 'losnam', los_name, /all
read_netcdf, file, 'R_LOS_1', r1, /all
read_netcdf, file, 'R_LOS_2', r2, /all
read_netcdf, file, 'Z_LOS_1', z1, /all
read_netcdf, file, 'Z_LOS_2', z2, /all
read_netcdf, file, 'PHI_LOS_1', phi1, /all
read_netcdf, file, 'PHI_LOS_2', phi2, /all
phi1 = phi1*!Pi/180.
phi2 = phi2*!Pi/180.

nlos = n_elements(r1)
fwhm_pix = dblarr(nlos) + 1.5
wslit_image = dblarr(nlos)
neon_done = 0
neon=1.

read_netcdf, file, 'fwhm_pix', fwhm_pix, /all
read_netcdf, file, 'wslit_shape', wslit_image, /all
read_netcdf, file, 'gratcons', gratcons, /gatr
read_netcdf, file, 'op_ang', op_ang, /gatr
read_netcdf, file, 'magnification', magnification, /gatr
read_netcdf, file, 'neon_done', neon_done, /gatr
read_netcdf, file, 'neon', neon, /all
if neon_done eq 2 then read_netcdf, file, 'neonfit', neonfit, /all else neonfit = neon

ichan_neon = 1
read_netcdf, file, 'ichan_neon', ichan_neon, /gatr
read_netcdf, file, 'lambda_neon', lambda_neon, /gatr

dat = {counts: spec, $
       lam:  lam, $
       time: time, $
       offset: offset, $
       sens: sens, $
       dwdp: dwdp, $
       exptime: exptime, $
       ctsph: ctsph, $
       sad2: sad2, $
       wslit: wslit, $
       wslit_image: wslit_image, $
       wlen: wlen, $
       pixw: pixw, $
       magnification: magnification, $
       fwhm_pix: fwhm_pix, $
       gratcons: gratcons, $
       op_ang: op_ang, $
       foc_len: foc_len, $
       neon: neon, $
       neonfit: neonfit, $
       neon_done: neon_done, $
       ichan_neon: ichan_neon, $
       lambda_neon: lambda_neon, $
       losnames: los_name, $
       r1: r1, $
       z1: z1, $
       phi1: phi1, $
       r2: r2, $
       z2: z2, $
       phi2: phi2 }
       

return
end

Function get_nii,shot,$
		 los=los,$
		 xr=xr,$
		 interelm=interelm,$
 		 channel=channel,$
		 use_evl=use_evl,$
		 append=append,$
		 diag=diag,$
		 sig3995=sig3995,$
		 rawelm=rawelm,$
		 elmcond=elmcond,$
		 no404=no404,$
		 use_rov8=use_rov8,$
		 dynamic=dynamic,$
		 preset=preset

shotstr  = string(shot,format='(i5)') 
if ~keyword_set(append)then append='data'
if ~keyword_set(los)then los='ROV-12' 
if ~keyword_set(diag)then diag='EVL' 
if ~keyword_set(sig3995)then sig3995='N_1_3995' 
if ~keyword_set(sig4041)then sig4041='N_1_4041' 
if ~keyword_Set(channel)then channel=22
if ~keyword_set(trace)then trace='save/'+shotstr+'/'+los+'-'+append+'.idl'


if ~keyword_set(use_evl)then begin
	restore,trace[0]
	time    = output.time
	nii_3995 = output.nii[*,0,0]
	nii_4041 = output.nii[*,0,1]
	nii_4026 = output.nii[*,0,2]
	nii_3995_err = output.nii_err[*,0,0]
	nii_4041_err = output.nii_err[*,0,1]
	nii_4026_err = output.nii_err[*,0,2]
	if keyword_set(xr)then begin
		id = where(time ge xr[0] and time le xr[1])
		time = time[id]
		nii_3995 = nii_3995[id,*]
		nii_4041 = nii_4041[id,*]
		nii_4026 = nii_4026[id,*]
		nii_3995_err = nii_3995_err[id,*]
		nii_4041_err = nii_4041_err[id,*]
		nii_4026_err = nii_4026_err[id,*]
	endif

endif else begin
	read_signal_mrm,0L,shot,diag,sig3995,time,nii_3995,2,exp=exp
	if ~keyword_set(no404)then read_signal_mrm,0L,shot,diag,sig4041,time,nii_4041,2,exp=exp
	if keyword_set(xr)then begin
		id = where(time ge xr[0] and time le xr[1])
	endif else id = findgen(n_elements(time))
	time = time[id]
	nii_3995 = nii_3995[id,channel]
	if ~keyword_set(no404)then nii_4041 = nii_4041[id,channel] else nii_4041 = nii_3995[id,channel]
	nii_4026 = nii_4041
	nii_3995_err = nii_3995*0.01
	nii_4041_err = nii_4041*0.01
	nii_4026_err = nii_4026*0.01
end
if ~keyword_set(xr)then xr = [min(time),max(time)]
if keyword_set(interelm)then begin
	print,'Using interELM ...'
	telm = find_elm(shot,time)
	if ~keyword_set(elmcond)then elmcond=4.5
	idelm     = where(telm ge elmcond)
	if idelm[0] eq -1 then idelm = findgen(n_elements(time))
endif else idelm = findgen(n_elements(time))

rawtime    = time
time       = time[idelm]
rawNii3995 = nii_3995
nii3995    = nii_3995[idelm,*]
rawNii4041 = nii_4041
nii4041    = nii_4041[idelm,*]
rawNii4026 = nii_4026
nii4026    = nii_4026[idelm,*]

rawNii3995_err = nii_3995_err
nii3995_err    = nii_3995_err[idelm,*]
rawNii4041_err = nii_4041_err
nii4041_err    = nii_4041_err[idelm,*]
rawNii4026_err = nii_4026_err
nii4026_err    = nii_4026_err[idelm,*]

read_signal_mrm,0L,shot,'MAC','Tdiv',x,y,2,exp=exp
id         = where(x ge min(time) and x le max(time))
output_tdiv= y[id]
time_tdiv  = x[id]
tdiv       = interpol(smooth(output_tdiv,40,/edge_truncate),time_tdiv,time)

return,{ rawnii3995:rawnii3995,$
         rawnii4041:rawnii4041,$
         rawnii4026:rawnii4026,$
         nii3995:nii3995,$
         nii4041:nii4041,$
         nii4026:nii4026,$
         rawnii3995_err:rawnii3995_err,$
         rawnii4041_err:rawnii4041_err,$
         rawnii4026_err:rawnii4026_err,$
         nii3995_err:nii3995_err,$
         nii4041_err:nii4041_err,$
         nii4026_err:nii4026_err,$
	 tdiv:tdiv,$
	 rawtime:rawtime,$
	 time:time}
End
	 
	 
	 
	 
	 

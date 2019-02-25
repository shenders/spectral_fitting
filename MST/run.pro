Pro run,interelm=interelm

; Setup the run sequence

shots    = [35157      , 35158     , 35165     , 35167]
diags    = ['GVL'      , 'FVL'     , 'UVS'     , 'UVS'    , 'MAC' ]
sigs     = ['Ar0_7504' , 'N_1_3995', 'N_tot   ', 'CFA03A' , 'Tdiv']
LOS_idx  = [0          , 22	   , 0         , 0        , 0     ]
yscale   = [1e17       , 1e19      , 1e22      , 1e21     , 1     ] 
ymax     = [1e20       , 1e20      , 1e20      , 1e20     , 20    ] 
tres     = 4E-3
trange   = [2.0,6.5]
npts     = FLOOR((trange[1]-trange[0])/tres)
timebase = findgen(npts)*(trange[1]-trange[0])/(npts-1)+trange[0]

intens   = fltarr(n_elements(shots),n_elements(diags),n_elements(timebase))

for s = 0,n_elements(shots)-1 do begin
	shot = shots[s]
	for i=0,n_elements(diags)-1 do begin
		read_signal_mrm,0L,shot,diags[i],sigs[i],time,sig,2		
		if keyword_set(interelm)then begin	
			telm       = find_elm(shot,time)
			if ~keyword_set(elmcond)then elmcond=4.5
			idback     = where(telm ge elmcond)
			if idback[0] ne -1 then begin
				time = time[idback]
				sig  = sig[idback,*]
			endif
		endif
		intens[s,i,*] = interpol(sig[*,LOS_idx[i]]/yscale[i],time,timebase)
	endfor
endfor

ncols=n_elements(shots)
nrows=n_elements(diags)
!p.multi=[0,ncols,nrows] 
adas_colors,colors=colors
window,/free,xs=1400,ys=1000
!p.charsize=2.1
for j=0,n_elements(diags)-1 do begin
	for i=0,n_elements(shots)-1 do begin
		plot,timebase,intens[i,j,*]>0,xtitle='Time [s]',ytitle='ph/s/m!u2!n/sr',yr=[0,ymax[j]<max(intens[*,j,*])],$
		title=string(shots[i],format='(i5)')+': '+diags[j]+': '+sigs[j],col=colors.black,back=colors.white
		if diags[j] ne 'UVS' then oplot,timebase,smooth(intens[i,j,*]>0,15),col=colors.red,thick=2.0
	endfor
endfor
stop
END

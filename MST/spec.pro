Pro spec,ar=ar,n=n,dens=dens,interelm=interelm,combine=combine,sf=sf,pngplot=pngplot,tr=tr

if ~keyword_set(sf)then sf=20
if keyword_set(ar)then begin
	shots     = [35157      , 35158   , 35165    , 35167          ]
	diags     = ['GVL'	,'UVS'    , 'UVS'    , 'MAC' , 'ELM'  ]
	sigs      = ['Ar0_7504' ,'N_tot'  , 'CFA03A' , 'Tdiv', 'f_ELM']
	LOS_idx   = [0  	,0	  , 0	     , 0     , 0      ]
	LOS_str   = ['ROV-09'	,'NULL'   , 'NULL'   , 'NULL', 'NULL' ]
	idmark    = [1  	,2	  , 3	     , 4     , 5      ]
	yscale    = [1e17	,1e22	  , 1e21     , 1     , 1      ]
	yscalestr = ['x1e17'	,'x1e22'  , 'x1e21'  , ''    , ''     ]
	ymax      = [0.6  	,3	  , 3	     , 20    , 250    ]
	ytit      = ['ph/s/m!u2!n/sr', '#/s' , '#/s' , 'eV'  , 'Hz'   ]
endif

if keyword_set(n)then begin
	shots    = [35157      , 35158     , 35165     , 35167                 ]
	diags    = ['EVL'      , 'EVL'     , 'EVL'     , 'EVL'     , 'FVL'     ]
	sigs     = ['N_1_3995' , 'N_1_3995', 'N_1_3995', 'N_1_3995', 'N_1_3995']
	LOS_idx  = [21         , 20	   , 22        , 23        , 22        ]
	LOS_str  = ['ROV-10'   , 'ROV-11'  , 'ROV-12'  , 'ROV-13'  , 'ROV-14'  ]
	idmark   = [0          , 1	   , 2         , 3         , 4         ]
	yscale   = [1e19       , 1e19	   , 1e19      , 1e19      , 1e19      ]
	yscalestr= ['x1e19'    , 'x1e19'   , 'x1e19'   , 'x1e19'   , 'x1e19'    ]
	ymax     = [1.         , 1. 	   , 1.        , 1.        , 1.        ]* 0.8 
	ytit     = ['ph/s/m!u2!n/sr', 'ph/s/m!u2!n/sr','ph/s/m!u2!n/sr','ph/s/m!u2!n/sr','ph/s/m!u2!n/sr'   ]
endif
if keyword_set(dens)then begin
	shots    = [35157      , 35158     , 35165     , 35167                  	   ]
	diags    = ['EVL'      , 'EVL'     , 'EVL'     , 'EVL'     , 'FVL'     , 'FVL'     ]
	sigs     = ['Ne1'      , 'Ne1'     , 'Ne1'     , 'Ne1'     , 'Ne'      , 'Ne'      ]
	LOS_idx  = [21         , 20	   , 22        , 23        , 22        , 8	   ]
	LOS_str  = ['ROV-10'   , 'ROV-11'  , 'ROV-12'  , 'ROV-13'  , 'ROV-14'  , 'ZON-01'  ]
	idmark   = [0          , 1	   , 2         , 3         , 4         , 5	   ]
	yscale   = [1e21       , 1e21	   , 1e21      , 1e21      , 1e21      , 1e19	   ]
	ymax     = [1.         , 1. 	   , 1.        , 1.        , 1.        , 1.0	   ]
endif

tres      = 4E-3
if keyword_set(tr)then trange=tr else trange    = [2.0,6.2]
npts      = FLOOR((trange[1]-trange[0])/tres)

intens      = -1
intens_elm  = -1
intens_elm_raw  = -1
time        = -1
time_elm    = -1
shotarr     = -1
diagarr     = -1
shotarr_elm = -1
diagarr_elm = -1
adas_colors,colors=colors
colpick     = [colors.navy,colors.blue,colors.aqua,colors.magenta,colors.red,colors.orange,colors.gold]

for s = 0,n_elements(shots)-1 do begin
	shot = shots[s]
	for i=0,n_elements(diags)-1 do begin
		read_signal_mrm,0L,shot,diags[i],sigs[i],timearr,sig,2		
		if keyword_set(interelm)then begin	
			telm       = find_elm(shot,timearr)
			if ~keyword_set(elmcond)then elmcond=4.5
			idback     = where(telm ge elmcond)
			if idback[0] ne -1 then begin
				timearr_elm = timearr[idback]
				sig_elm     = sig[idback,*]
			endif
		endif else begin
			timearr_elm = timearr
			sig_elm  = sig
		end
		intens     = [intens,sig[*,LOS_idx[i]]/yscale[i]]
		intens_elm = [intens_elm,smooth(sig_elm[*,LOS_idx[i]]/yscale[i],sf)]
		intens_elm_raw = [intens_elm_raw,sig_elm[*,LOS_idx[i]]/yscale[i]]
		time       = [time,timearr]
		time_elm   = [time_elm,timearr_elm]
		shotarr    = [shotarr,fltarr(n_elements(timearr))+shot]
		diagarr    = [diagarr,fltarr(n_elements(timearr))+idmark[i]]
		shotarr_elm= [shotarr_elm,fltarr(n_elements(timearr_elm))+shot]
		diagarr_elm= [diagarr_elm,fltarr(n_elements(timearr_elm))+idmark[i]]
	endfor
endfor
intens      = intens[1:*]	 
intens_elm  = intens_elm[1:*] 
intens_elm_raw  = intens_elm_raw[1:*] 
time        = time[1:*]	 
time_elm    = time_elm[1:*]   
shotarr     = shotarr[1:*]
diagarr     = diagarr[1:*]
shotarr_elm = shotarr_elm[1:*]
diagarr_elm = diagarr_elm[1:*]

if keyword_set(combine)then begin
	nrows = 1 
	window,/free,xs=700,ys=1000
endif else begin 
	nrows = n_elements(diags)
	window,/free,xs=1400,ys=1000
end
ncols=n_elements(shots)

if keyword_set(combine)then !p.multi=[0,nrows,ncols] else !p.multi=[0,ncols,nrows] 

!p.charsize=2.5

if keyword_set(combine)then begin
for j=0,n_elements(shots)-1 do begin
	for i=0,n_elements(diags)-1 do begin
		id     = where(shotarr eq shots[j] and diagarr eq idmark[i])
		id_elm = where(shotarr_elm eq shots[j] and diagarr_elm eq idmark[i])
		if i eq 0 then begin
			plot,time_elm[id_elm],intens_elm[id_elm]>0,xtitle='Time [s]',xs=1,ys=1,/nodata,ytitle=ytit[i],xr=trange,yr=[0,ymax[i]<max(intens)],$
			title=string(shots[j],format='(i5)')+': '+diags[i]+': '+sigs[i]+' '+yscalestr[i],col=colors.black,back=colors.white
		endif 
		oplot,time_elm[id_elm],intens_elm_raw[id_elm]>0,col=colpick[i],thick=0.5
		oplot,time_elm[id_elm],intens_elm[id_elm]>0,col=colpick[i],thick=2.0
	endfor
endfor

endif else begin
for i=0,n_elements(diags)-1 do begin
	for j=0,n_elements(shots)-1 do begin
		id     = where(shotarr eq shots[j] and diagarr eq idmark[i])
		id_elm = where(shotarr_elm eq shots[j] and diagarr_elm eq idmark[i])
		plot,time[id],intens[id]>0,xtitle='Time [s]',ytitle=ytit[i],xs=1,ys=1,xr=trange,yr=[0,ymax[i]<max(intens)],$
		title=string(shots[j],format='("#",i5)')+': '+diags[i]+': '+LOS_str[i]+': '+sigs[i]+' '+yscalestr[i],col=colors.black,back=colors.white
		if keyword_set(interelm)then begin
			if diags[i] ne 'UVS' then oplot,time_elm[id_elm],intens_elm[id_elm]>0,col=colors.red
		endif
	endfor
endfor
end
if keyword_set(pngplot)then write_png,'out.png',tvrd(/true)

if keyword_set(combine)then begin
nrows = 1 
window,/free,xs=700,ys=1000
yarr = fltarr(n_elements(diags))
for j=0,n_elements(shots)-1 do begin
	timewindows = [3.5,4.5,5.5]
	for ii = 0, n_elements(timewindows)-1 do begin	
		for i=0,n_elements(diags)-1 do begin
			id_elm = where(shotarr_elm eq shots[j] and diagarr_elm eq idmark[i])		
			xtemp = time_elm[id_elm]
			ytemp = intens_elm[id_elm]		
			id = where(abs(xtemp - timewindows[ii]) eq min(abs(xtemp-timewindows[ii])))
			yarr[i] = ytemp[id]
		endfor
		
		xax    = [10,11,12,13,14]
		if ii eq 0 then begin
			plot,xax,yarr>0,xtitle='ROV-##',xs=1,ys=1,/nodata,ytitle=ytit[0],yr=[0,ymax[0]<max(intens)],$
			title=string(shots[j],format='(i5)')+' '+yscalestr[0],col=colors.black,back=colors.white
		endif 
		user_psym,1,/fill
		oplot,xax,yarr,col=colpick[ii+1],thick=2.0,psym=-8
	endfor
endfor
if keyword_set(pngplot)then write_png,'out2.png',tvrd(/true)
endif

stop
End

Function gauss_lines,species=species,instr_func=instr_func,diag=diag,unknown=unknown

if ~keyword_set(instr_func)then instr_func=0.08
if ~keyword_set(diag)then diag='EVS'
if ~keyword_set(species)then species = 'NONE'

; Helium lines
;He_lines = [396.50  ,402.62  ]
;He_ions  = ['He1'   , 'He1'  ]
;He_trn   = ['4p2s'  ,'4p2s' ]
;He_cpl   = [0.6     , 1.00]

; Carbon lines
C_lines  = [406.817 , 406.914 , 407.013 , 407.219, 407.588 ]
C_ions   = ['C3'    , 'C3'    , 'C3'    , 'C2'   , 'C2'    ]
C_trn    = ['5g4f'  , '5g4f'  , '5g4f'  , 'xxxx' , 'xxxx'  ]
C_cpl    = [0.470   , 0.537   , 1.00    , 0.00   , 0.0     ]


; Nitrogen lines
N_lines  = [399.50  , 399.863 , 400.358 , 402.62  , 403.510 , 403.935 , 404.15  , $
            404.37  , 404.490 , 405.68  , 405.776 , 407.31  , 407.69  , 408.15  , $
            408.21  , 408.73  , 409.73  , 410.34  , 395.585 , 744.229 , 746.831   ]
N_ions   = ['N2'    , 'N3'    , 'N3'	, 'N2'    , 'N2'    , 'N2'    , 'N2'	, $
            'N2'    , 'N2'    , 'N2'	, 'N4'    , 'N2'    , 'N2'    , 'xx'    , $ 
            'N2'    , 'N2'    , 'N3'	, 'N3'    , 'N2'    , 'N1'    , 'N1'      ]
N_trn    = ['3p3s'  , '5f4d'  , '5f4d'  , '4f3d'  , '4f3d'  , '4f3d'  , '4f3d'  , $
            '4f3d'  , '4f3d'  , '4f3d'  , '3d3p'  , '4f3d'  , '4f3d'  , 'xxxx'  , $
            '4f3d'  , '4f3d'  , '3p3s'  , '3p3s'  , '3p3s'  , '3p3s'  , '3p3s'    ]
N_cpl    = [0.00    , 0.60    , 1.00	, 2.00    , 2.378   , 1.070   , 3.00    , $
            2.326   , 2.0587  , 2.05913 , 0.00    , 4.00    , 3.38    , 0.00	, $
            3.85    , 0.00    , 5.0	, 4.50    , 0.0     , 0.0     , 0.0       ]

; Oxygen lines
O_lines  = [397.30  , 398.27  ,408.92  , 409.29]
O_ions   = ['O1'    , 'O1'    ,'O2'    , 'O2'  ]
O_trn    = ['3p3s'  , '3d3p'  ,'4f3d'  , '4f3d']
O_cpl    = [0.00    , 0.00    ,0.00    , 0.00  ]

; Neon lines
Ne_lines = [ 375.120 , 370.960 , 377.710 , 373.490 , 366.410 , 376.630 , 369.420 , $
             372.710 , 364.390 , 371.310 , 381.360 , 371.350 ]
Ne_trn	 = [ '3p3s'  , '3p3s'  , '3p3s'  , '3p3s'  , '3p3s'  , '3p3s'  , '3p3s'  , $
             '3p3s'  , '3p3s'  , '3p3s'  , '3p3s'  , '3p3s'  ]
Ne_ions  = [ 'Ne2'   , 'Ne2'   , 'Ne2'   , 'Ne2'   , 'Ne2'   , 'Ne2'   , 'Ne2'   , $
             'Ne2'   , 'Ne2'   , 'Ne2'   , 'Ne4'   , 'Ne4'   ]
Ne_cpl	 = [ 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , 0.00    , $
             0.00    , 0.00    , 0.00    , 0.00    , 0.00    ]
; Ar lines
Ar_lines = [ 750.387 , 750.515 , 751.465 ]
Ar_trn	 = [ '4p4s'  , '4d4p'  , '4p4s'  ]
Ar_ions	 = [ 'Ar1'   , 'Ar2'   , 'Ar1'   ]
Ar_cpl	 = [ 0.00    , 0.00    , 0.00    ]

; Tungsten lines
W_lines  = [400.876 , 407.449]
W_ions   = ['W1'    , 'W1'   ]
W_trn    = ['6p6s'  , '6p6s' ]
W_cpl    = [0.00    , 0.00   ]

; Unknown lines
; Put He I lines into unknown so that wavelength positions are fixed [396.50 and 402.62]
extra_line = 0
if keyword_set(extra_line)then begin
	id = where(species eq 'He')
	if id[0] ne -1 then begin	
		X_lines  = [396.50  ,402.62, 399.29 , 399.35 , 399.57, 396.20  , 396.89  , 399.7   , 401.38  , 401.493 , $
            	402.478 , 402.533 , 402.6, 402.7, 407.88  , 408.06  , 408.82  ,409.55 , 410.6 , 750.0 , 750.9, 751.99,409.87,410.47, 395.42]
		X_ions   = strarr(n_elements(X_lines))+'xx'
		X_trn    = strarr(n_elements(X_lines))+'xxxx'
		X_cpl    = strarr(n_elements(X_lines))+0.0
		X_cpl[3:4]   = [1.0, 0.65]
		X_cpl[10:13] = [2.0, 1.55,1.3,1.05]
		X_cpl[0:1]   = [2.6, 3.0] ; He I lines
		X_ions[0:1]  = ['He1','He1'] ; He I lines
		X_trn[0:1]  = ['4p2s'  ,'4p2s'] ; He I lines
	endif else begin
		X_lines  = [399.29 , 399.35 , 399.57, 396.20  , 396.89  , 399.7   , 401.38  , 401.493 , $
            	402.478 , 402.533 , 402.6, 402.7,407.88  , 408.06  , 408.82  ,409.55 , 410.6 , 750.0 , 750.9, 751.99,409.87,410.47, 395.42]
		X_ions   = strarr(n_elements(X_lines))+'xx'
		X_trn    = strarr(n_elements(X_lines))+'xxxx'
		X_cpl    = strarr(n_elements(X_lines))+0.0
		X_cpl[1:2]   = [1.0, 0.65]
		X_cpl[8:11] = [2.0, 1.55,1.3,1.05]
	end
endif else begin
	id = where(species eq 'He')
	if id[0] ne -1 then begin	
		X_lines  = [396.472  ,402.62, 399.29 , 399.35 , 399.57, 396.20  , 396.89  , 399.7   , 401.38  , 401.493 , $
            	402.478 , 402.533 , 407.88  , 408.06  , 408.82  ,409.57 , 410.6 , 750.0 , 750.9, 751.99,409.83,410.45, 395.42]
		X_ions   = strarr(n_elements(X_lines))+'xx'
		X_trn    = strarr(n_elements(X_lines))+'xxxx'
		X_cpl    = strarr(n_elements(X_lines))+0.0
		X_cpl[3:4]   = [1.0, 0.65]
		X_cpl[10:11] = [2.0, 1.69]
		X_cpl[0:1]   = [2.6, 3.0] ; He I lines
		X_ions[0:1]  = ['He1','He1'] ; He I lines
		X_trn[0:1]  = ['4p2s'  ,'4p2s'] ; He I lines
	endif else begin
		X_lines  = [399.29 , 399.35 , 399.57, 396.20  , 396.89  , 399.7   , 401.38  , 401.493 , $
            	402.478 , 402.533 , 407.88  , 408.06  , 408.82  ,409.57 , 410.6 , 750.0 , 750.9, 751.99,409.83,410.45, 395.42]
		X_ions   = strarr(n_elements(X_lines))+'xx'
		X_trn    = strarr(n_elements(X_lines))+'xxxx'
		X_cpl    = strarr(n_elements(X_lines))+0.0
		X_cpl[1:2] = [1.0, 0.65]
		X_cpl[8:9] = [2.0, 1.69]
	end
end
pos = 0.0
ion = '-1'
trn = '-1'		     
cpl = 0.0

;id = where(species eq 'He')
;if id[0] ne -1 then begin
;    pos = [pos,He_lines]
;    ion = [ion,He_ions]
;    trn = [trn,He_trn]
;    max_cpl = max(cpl)
;    id = where(He_cpl > 0)
;    if id[0] ne -1 then He_cpl[id] = He_cpl[id] + max_cpl
;    cpl = [cpl,He_cpl]
;endif
id = where(species eq 'C')
if id[0] ne -1 then begin
    pos = [pos,C_lines]
    ion = [ion,C_ions]
    trn = [trn,C_trn]
    max_cpl = max(cpl)
    id = where(C_cpl > 0)
    if id[0] ne -1 then C_cpl[id] = C_cpl[id] + max_cpl
    cpl = [cpl,C_cpl]
endif
id = where(species eq 'N')
if id[0] ne -1 then begin
    pos = [pos,N_lines]
    ion = [ion,N_ions]
    trn = [trn,N_trn]
    max_cpl = max(cpl)
    id = where(N_cpl > 0)
    if id[0] ne -1 then N_cpl[id] = N_cpl[id] + max_cpl
    cpl = [cpl,N_cpl]
endif
id = where(species eq 'O')
if id[0] ne -1 then begin
    pos = [pos,O_lines]
    ion = [ion,O_ions]
    trn = [trn,O_trn]
    max_cpl = max(cpl)
    id = where(O_cpl > 0)
    if id[0] ne -1 then O_cpl[id] = O_cpl[id] + max_cpl
    cpl = [cpl,O_cpl]
endif
id = where(species eq 'Ne')
if id[0] ne -1 then begin
    pos = [pos,Ne_lines]
    ion = [ion,Ne_ions]
    trn = [trn,Ne_trn]
    max_cpl = max(cpl)
    id = where(Ne_cpl > 0)
    if id[0] ne -1 then Ne_cpl[id] = Ne_cpl[id] + max_cpl
    cpl = [cpl,Ne_cpl]
endif
id = where(species eq 'Ar')
if id[0] ne -1 then begin
    pos = [pos,Ar_lines]
    ion = [ion,Ar_ions]
    trn = [trn,Ar_trn]
    max_cpl = max(cpl)
    id = where(Ar_cpl > 0)
    if id[0] ne -1 then Ar_cpl[id] = Ar_cpl[id] + max_cpl
    cpl = [cpl,Ar_cpl]
endif
id = where(species eq 'W')
if id[0] ne -1 then begin
    pos = [pos,W_lines]
    ion = [ion,W_ions]
    trn = [trn,W_trn]
    max_cpl = max(cpl)
    id = where(W_cpl > 0)
    if id[0] ne -1 then W_cpl[id] = W_cpl[id] + max_cpl
    cpl = [cpl,W_cpl]
endif

if keyword_set(unknown) then begin
    pos = [pos,X_lines]
    ion = [ion,X_ions]
    trn = [trn,X_trn]
    max_cpl = max(cpl)
    id = where(X_cpl > 0)
    if id[0] ne -1 then X_cpl[id] = X_cpl[id] + max_cpl
    cpl = [cpl,X_cpl]
endif

if n_elements(pos) > 1 then begin
    pos = pos[1:*]
    ion = ion[1:*]
    trn = trn[1:*]
    cpl = cpl[1:*]
endif

gauss        = {pos:pos, ion:ion, trn:trn, couple:cpl, fwhm: instr_func, erc: fltarr(1000)+0.0, diag:diag}
idsort	     = sort(gauss.pos)
gauss.pos    = gauss.pos[idsort]
gauss.ion    = gauss.ion[idsort]
gauss.couple = gauss.couple[idsort]
gauss.trn    = gauss.trn[idsort]
return, gauss

End

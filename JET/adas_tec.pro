;-----------------------------------------------------------------------------
;+
;
; NAME     : adas_tec
;
; PURPOSE  : 
;   	    Function to calculate the dominant total emission coefficients (TEC) for 
;   	    a particular element over a pre-defined spectrometer range and resolution.
;   	    These TECs are a combination of the ionisation balance and PEC.
;
; ROUTINES: 
;           Name                Description
;   	    ----  	    	-----------
;
;   	    xxeiz0  	    	ADAS: Converts element symbol to nuclear charge
;   	    xxeiam  	    	ADAS: Returns element mass
;   	    adas_colors	    	ADAS: Returns color codes for plotting
;   	    run_adas405    	ADAS: Calculate ionisation balance
;   	    read_adf15    	ADAS: Reads adf15 files
;   	    split_multiplet_ssh   	ADAS: Statistically split term resolved cross-sections  	       
;
; PARAMETERS:
;   	    Inputs   	    	Description


;   	    ------  	    	-----------
;	    	       	   
;   	    temp    	*** 	Electron temperature / eV
;   	    dens    	*** 	Electron density / cm-3
;   	    elem(,)    	*** 	Element symbol
;   	    wmin    	*** 	Minimum wavelength / Ang
;   	    wmax    	*** 	Maximum wavelength / Ang
;   	    resol   	*** 	Spectrometer resolution / Ang
;   	    inst_fwhm   *** 	Spectrometer instrument function / Ang
;   
;   	    Keywords   	    	Description
;   	    -------- 	    	-----------
;
;   	    verbose 	*** 	Print code progress to screen
;   	    doplot  	*** 	Plot tec in IDL plot window
;   	    adashome  	*** 	Directory of adas folder (default './')
;   	    adf15_code 	*** 	Code naming of adf04 files (default 'pecssh')
;   	    year 	*** 	Year of data for ionisation balance (default '96')
;
;   	    Outputs 	    	Description   
;   	    -------  	    	-----------
;
;   	    wavelength	*** 	Wavelength array / Ang
;   	    tec     	*** 	Total emission coefficient array / s-1
;   	    err     	*** 	Error status: -1 failed, 0 success
;   	    message 	*** 	Error message
;
; NOTES:
;   	    Uses the SPAWN command therefore IDL v5.3 cannot be used. 
;   	
;   	    This routine is designed to integrate into the python
;   	    program [] and therefore ADF04 and ADF15 files have been
;   	    provided for this particular use. However these files can  
;   	    be aligned to the ADAS database files using appropriate keywords. 
;
;   	    This code does use ADAS system codes, therefore IDL pathways
;   	    must include adashome/idl/*. This can be setup in bash using:
;           . /home/adas/setup/adas_setup.ksh

;   	    Currently, the code assumes the adf04 and adf15 naming structure:
;   	    adf15 : adashome/adas/adf15/[adf15_code]#[z0]/[adf15_code]#[z0]_[llu]#[ion].dat
;   	
;   	    Doppler broadening assumes that Timp = Te
;
;   	    A typical command line call:
;   	
;   	    IDL> temp   = 10.0   ; eV
;   	    IDL> dens   = 1e13   ; cm-3
;   	    IDL> elem   = 'ne'   ; neon
;   	    IDL> wmin   = 3600.0 ; 3600 Ang minimum
;   	    IDL> wmin   = 4000.0 ; 4000 Ang maximum
;   	    IDL> resol  = 0.5    ; spectrometer resolution Ang 
;   	    IDL> inst_fwhm  = 1.0    ; spectrometer instrument function Ang 
;   	    IDL> result = adas_tec(temp,dens,elem,wmin,wmax,resol,inst_fwhm)
;   	    IDL> plot,result.wavelength,result.tec
;
; CATEGORY:
;           Adas system.
;
; WRITTEN:
;           S. S. Henderson, University of Strathclyde.
;
; MODIFIED:
;           1.1   S. Henderson - initial version    
;
; VERSION:
;           1.1   15-10-15
;-
;-----------------------------------------------------------------------------
PRO psplot,file=file
	IF ~KEYWORD_SET(file)THEN file='plot.ps'
	!P.FONT=0
	!X.THICK=5
	!Y.THICK=5
	!Z.THICK=5
	!P.THICK=5
	SET_PLOT,'ps'
	DEVICE,COLOR=1,XSIZE=8,YSIZE=6,/INCHES,BITS_PER_PIXEL=64,FILE=FILE,$
	FONT_SIZE=9,YOFFSET=0.1
END
FUNCTION norm_split_multiplet,l_lw,s_up,l_up
    IF l_lw EQ l_up   THEN norm = 1.0 / ((2.0*s_up+1.0)*(2.0*l_up+1.0)*l_up*(l_up+1.0))           
    IF l_lw EQ l_up+1 THEN norm = 1.0 / ((2.0*s_up+1.0)*(2.0*l_up+1.0)*(l_up+1.0)*(2.0*l_up+3.0))  
    IF l_lw EQ l_up+2 THEN norm = 1.0 / ((2.0*s_up+1.0)*(2.0*l_up+1.0)*(l_up+1.0)*(2.0*l_up+3.0))  
    IF l_lw EQ l_up-1 THEN norm = 1.0 / ((2.0*s_up+1.0)*(2.0*l_up+1.0)*l_up*(2.0*l_up-1.0))       
    IF l_lw EQ l_up-2 THEN norm = 1.0 / ((2.0*s_up+1.0)*(2.0*l_up+1.0)*l_up*(2.0*l_up-1.0))       
    RETURN,norm
END    	
FUNCTION adas_tec, temp , dens , elem , wmin , wmax , resol , inst_fwhm , $              
		   verbose    = verbose    , doplot     = doplot        , $
		   adashome   = adashome   , adf15_code = adf15_code    , $
		   year       = year       , yscale     = yscale        , $
		   psplot     = psplot     , elem_mul   = elem_mul

;   *********************************************************
;   **** Step 0.0: Initialise inputs             	 ****
;   *********************************************************

    IF ~KEYWORD_SET(adashome)THEN adashome = '/work/shenders/MST2/'
    IF ~KEYWORD_SET(adf15_code)THEN adf15_code = 'pec96'
    IF ~KEYWORD_SET(year)THEN year = '96'
    IF ~KEYWORD_SET(elem_mul)THEN elem_mul = FLTARR(N_ELEMENTS(elem))+1.0

    wv  = 0.0
    gf  = 0.0 
    ab  = 0.0 
    el  = ''
    js  = 0.0
    
    FOR ii = 0, N_ELEMENTS(elem)-1 DO BEGIN

;   *********************************************************
;   **** Step 1.0: Calculate ionisation balance     	 ****
;   *********************************************************

    IF KEYWORD_SET(verbose) THEN PRINT,'Running ADAS405 to calculate ionisation balance ...'
    
    
    z0_nuc      = xxeiz0(elem(ii))
    z0_nuc_str  = STRCOMPRESS(STRING(z0_nuc,FORMAT='(I2)'),/REMOVE_ALL)
    run_adas405,    uid  = 'adas'  , $
    		    year = year	   , $
    		    elem = elem(ii), $
    		    te   = temp    , $
    		    dens = dens    , $
    		    frac = frac

;   *********************************************************
;   **** Step 2.0: NIST+ADAS Database                 	 ****
;   *********************************************************

    adf15_dir   = adashome+'adas/adf15/'+adf15_code+'#'+z0_nuc_str+'/'

    IF z0_nuc EQ 10 THEN BEGIN
;	********************************	
;	**** NeII Transition: 	    ****	
;	**** 2s2 2p4 3p 4P  -->     ****
;	**** 2s2 2p4 3s 4P          ****
;	********************************	
	
	j_tran        = [[0.5,0.5],[0.5,1.5],[1.5,0.5],[1.5,1.5],[1.5,2.5],[2.5,1.5],[2.5,2.5]]		 
	NIST_wv       = [3751.246 ,3709.622 ,3777.136 ,3734.939 ,3664.074 ,3766.259 ,3694.215 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#ne_pju#ne1.dat'
	read_adf15,file=file,block=17,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=70,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.5,1,1.5,1,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,1.5,1.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** NeII Transition: 	    ****	
;	**** 2s2 2p4 3p 2D  -->     ****
;	**** 2s2 2p4 3s 2P          ****
;	********************************	
	
	j_tran        = [[1.5,0.5],[1.5,1.5],[2.5,1.5]]		 
	NIST_wv       = [3727.108 ,3643.929 ,3713.083 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#ne_pju#ne1.dat'
	read_adf15,file=file,block=24,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=77,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.5,1.0,0.5,2.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,0.5,2.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** NeIV Transition: 	    ****	
;	**** 2s2 2p2 3p 2S -->      ****
;	**** 2s2 2p2 3s 2P          ****
;	********************************	
	j_tran        = [[0.5,0.5],[1.5,0.5]]		 
	NIST_wv       = [3713.500 ,3813.670 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 3.0
	file     = adf15_dir + 'pec96#ne_pju#ne3.dat'
	read_adf15,file=file,block=8,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=53,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.5,0.0,0.5,1.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(0.0,0.5,1.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF
    ENDIF
    IF z0_nuc EQ 7 THEN BEGIN
;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 3p 1D -->       ****
;	**** 2s2 2p 3s 1P           ****
;	********************************	
	j_tran        = [[2.0,1.0]]		 
	NIST_wv       = [3995.00 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=2,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=52,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.0,1.0,0.0,2.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,0.0,2.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 3p 1S -->       ****
;	**** 2s2 2p 3s 1P           ****
;	********************************	
	j_tran        = [[1.0,0.0]]		 
	NIST_wv       = [3437.150 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=7,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=57,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.0,1.0,0.0,0.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,0.0,0.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 3d 1P -->       ****
;	**** 2s2 2p 3p 1P           ****
;	********************************	
	j_tran        = [[1.0,1.0]]		 
	NIST_wv       = [3919.000 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=40,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=90,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.0,1.0,0.0,1.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,0.0,1.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 4s 3P -->       ****
;	**** 2s2 2p 3p 3P           ****
;	********************************	
	j_tran        = [[0.0,1.0],[1.0,0.0],[1.0,1.0],[1.0,2.0],[2.0,1.0],[2.0,2.0]]		 
	NIST_wv       = [3855.096 ,3842.187 ,3847.997 ,3856.062 ,3829.795 ,3838.374 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=15,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=65,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.0,1.0,1.0,1.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,1.0,1.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 3p 1D -->       ****
;	**** 2s2 2p 3s 3P           ****
;	********************************	

	j_tran        = [[2.0,1.0]]		 
	NIST_wv       = [3955.85 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=19,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=69,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.0,1.0,0.0,2.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,0.0,2.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 4f 3G -->       ****
;	**** 2s2 2p 3d 3F           ****
;	********************************	

	j_tran        = [[3.0,2.0],[4.0,3.0],[5.0,4.0],[3.0,3.0],[4.0,4.0]]		 
	NIST_wv       = [4035.080 ,4043.530 ,4041.310 ,4045.920 ,4058.050 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=8,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=58,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.0,3.0,1.0,4.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(3.0,1.0,4.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 4f 1G -->       ****
;	**** 2s2 2p 3d 3F           ****
;	********************************	
	j_tran        = [[4.0,3.0]]		 
	NIST_wv       = [4027.220]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=36,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=86,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.0,3.0,0.0,4.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(3.0,0.0,4.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 4f 1F -->       ****
;	**** 2s2 2p 3d 1D           ****
;	********************************	

	j_tran        = [[3.0,2.0]]		 
	NIST_wv       = [4177.340 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=28,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=78,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.0,2.0,0.0,3.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(2.0,0.0,3.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 4f 3F -->       ****
;	**** 2s2 2p 3d 1D           ****
;	********************************	

	j_tran        = [[2.0,2.0],[3.0,2.0]]		 
	NIST_wv       = [4176.840 ,4172.770 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=27,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=77,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.0,2.0,1.0,3.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(2.0,1.0,3.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 4s 1P -->       ****
;	**** 2s2 2p 3p 1D           ****
;	********************************	

	j_tran        = [[2.0,1.0]]		 
	NIST_wv       = [4227.740 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=42,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=92,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.0,2.0,0.0,1.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(2.0,0.0,1.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 4f 3F -->       ****
;	**** 2s2 2p 3d 3D           ****
;   	**** Based on NIST and Cowan****
;	********************************	

	j_tran        = [[2.0,1.0],[2.0,2.0],[2.0,3.0],[3.0,2.0],[3.0,3.0],[4.0,3.0]]		 
	NIST_wv       = [4237.050 ,4242.430 ,4247.900 ,4236.910 ,4243.700 ,4241.780 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=10,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=60,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.0,2.0,1.0,3.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(2.0,1.0,3.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 4f 1F -->       ****
;	**** 2s2 2p 3d 3D           ****
;	********************************	

	j_tran        = [[3.0,2.0]]		 
	NIST_wv       = [4242.950];[4239.910]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=39,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=89,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.0,2.0,0.0,3.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(2.0,0.0,3.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 4f 3D -->       ****
;	**** 2s2 2p 3d 3D           ****
;	********************************	

	j_tran        = [[1.0,1.0],[1.0,2.0],[2.0,1.0],[2.0,2.0],[2.0,3.0],[3.0,2.0],[3.0,3.0]]		 
	NIST_wv       = [4158.170 ,4162.330 ,4170.560 ,4174.740 ,4180.030 ,4186.930 ,4180.850 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=33,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=83,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.0,2.0,1.0,2.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(2.0,1.0,2.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N II Transition: 	    ****	
;	**** 2s2 2p 4f 3G -->       ****
;	**** 2s2 2p 3d 3D           ****
;	********************************	

	j_tran        = [[4.0,3.0],[3.0,2.0],[3.0,3.0]]		 
	NIST_wv       = [4201.160 ,4197.150 ,4202.510 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 1.0
	file     = adf15_dir + 'pec96#n_vsu#n1.dat'
	read_adf15,file=file,block=30,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=80,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.0,2.0,1.0,4.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(2.0,1.0,4.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N III Transition: 	    ****	
;	**** 2s2 3p 2P -->          ****
;	**** 2s2 3s 2S              ****
;	********************************	

	j_tran        = [[0.5,0.5],[1.5,0.5]]		 
	NIST_wv       = [4103.430 ,4097.330 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 2.0
	file     = adf15_dir + 'pec96#n_vsu#n2.dat'
	read_adf15,file=file,block=1,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=51,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.5,0.0,0.5,1.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(0.0,0.5,1.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
	ENDIF

;	********************************	
;	**** N III Transition: 	    ****	
;	**** 2s2p3d 2D -->          ****
;	**** 2s2p3p 2P              ****
;	********************************	

	j_tran        = [[1.5,0.5],[1.5,1.5],[2.5,1.5]]		 
	NIST_wv       = [3934.500 ,3942.880 ,3938.520]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 2.0
	file     = adf15_dir + 'pec96#n_vsu#n2.dat'
	read_adf15,file=file,block=11,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=61,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.5,1.0,0.5,2.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,0.5,2.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N III Transition: 	    ****	
;	**** 2s2p3p 4S -->          ****
;	**** 2s2p3s 4P              ****
;	********************************	

	j_tran        = [[0.5,1.5],[1.5,1.5],[2.5,1.5]]		 
	NIST_wv       = [3745.950 ,3754.670 ,3771.050]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 2.0
	file     = adf15_dir + 'pec96#n_vsu#n2.dat'
	read_adf15,file=file,block=7,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=57,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.5,1.0,1.5,0.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,1.5,0.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N III Transition: 	    ****	
;	**** 2s2p3d 4P -->          ****
;	**** 2s2p3p 4D              ****
;	********************************	

	j_tran        = [[0.5,0.5],[0.5,1.5],[1.5,0.5],[1.5,1.5],[1.5,2.5],[2.5,1.5],[2.5,2.5]]		 
	NIST_wv       = [3752.630 ,3757.650 ,3757.570 ,3762.600 ,3771.360 ,3770.360 ,3779.160 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 2.0
	file     = adf15_dir + 'pec96#n_vsu#n2.dat'
	read_adf15,file=file,block=26,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=76,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.5,2.0,1.5,1.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(2.0,1.5,1.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N III Transition: 	    ****	
;	**** 2s2p3p 2D -->          ****
;	**** 2s2p3s 2P              ****
;	********************************	

	j_tran        = [[1.5,0.5],[1.5,1.5],[2.5,1.5]]		 
	NIST_wv       = [4195.760 ,4215.770 ,4200.100]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 2.0
	file     = adf15_dir + 'pec96#n_vsu#n2.dat'
	read_adf15,file=file,block=9,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=59,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.5,1.0,0.5,2.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,0.5,2.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N IV Transition: 	    ****	
;	**** 1s2 2s 3d 1D -->       ****
;	**** 1s2 2s 3p 1P           ****
;	**** 4057.76 Ang 	    ****
;	********************************	

	j_tran        = [[2.0,1.0]]		 
	NIST_wv       = [4057.760 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 3.0
	file     = adf15_dir + 'pec96#n_pju#n3.dat'
	read_adf15,file=file,block=28,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=78,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.0,1.0,0.0,2.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,0.0,2.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N IV Transition: 	    ****	
;	**** 1s2 2s 3d 3P -->       ****
;	**** 1s2 2s 3p 3S           ****
;	********************************	

	j_tran        = [[0.0,1.0],[1.0,1.0],[2.0,1.0]]		 
	NIST_wv       = [3484.960 ,3482.990 ,3478.710 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 3.0
	file     = adf15_dir + 'pec96#n_pju#n3.dat'
	read_adf15,file=file,block=13,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=63,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.0,0.0,1.0,1.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(0.0,1.0,1.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N IV Transition: 	    ****	
;	**** 1s2 2s 3d 3P -->       ****
;	**** 1s2 2s 3p 3P           ****
;	********************************	

	j_tran        = [[0.0,1.0],[1.0,0.0],[1.0,1.0],[2.0,1.0],[1.0,2.0],[2.0,2.0]]		 
	NIST_wv       = [3461.360 ,3461.360 ,3454.650 ,3443.610 ,3474.530 ,3463.370 ]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN
	ionstage = 3.0
	file     = adf15_dir + 'pec96#n_pju#n3.dat'
	read_adf15,file=file,block=12,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=62,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(1.0,1.0,1.0,1.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,1.0,1.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF

;	********************************	
;	**** N IV Transition: 	    ****	
;	**** 1s2 2s 3d 1D -->       ****
;	**** 1s2 2s 3p 1P           ****
;	********************************	

	j_tran        = [[2.0,1.0]]		 
	NIST_wv       = [3747.540]		 
	IF MIN(NIST_wv) LE wmax AND MAX(NIST_wv) GE wmin THEN BEGIN 
	ionstage = 3.0
	file     = adf15_dir + 'pec96#n_pju#n3.dat'
	read_adf15,file=file,block=21,te=temp,dens=dens,data=data1	; Excitation
	read_adf15,file=file,block=71,te=temp,dens=dens,data=data2	; Recombination
	frac_split    = split_multiplet_ssh(0.0,1.0,0.0,2.0,j_low=j1,j_up=j2)	; Multiplet split
	norm          = norm_split_multiplet(1.0,0.0,2.0)
	data          = data1 * frac.ion(*,ionstage) + data2 * frac.ion(*,ionstage+1)
	FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
	    FOR i=0,N_ELEMENTS(NIST_wv)-1 DO BEGIN
		gf  = [gf,data(j) * elem_mul(ii) * frac_split(WHERE(j_tran(0,i) EQ j2 AND j_tran(1,i) EQ j1)) * norm ]
		js  = [js,j]
		wv  = [wv,NIST_wv(i)]
		ab  = [ab,ionstage]
		el  = [el,elem(ii)]
	    ENDFOR
	ENDFOR
    	ENDIF
    ENDIF
    
    ENDFOR
    wv_all  = wv[1:*]
    gf_all  = gf[1:*] 
    ab_all  = ab[1:*] 
    el_all  = el[1:*] 
    js_all  = js[1:*]
    xxeiam,elem,amss 
    npix       = 1024
    wavelength = wmin+((wmax-wmin)/(npix-1))*FINDGEN(npix)
    tec        = FLTARR(npix,n_elements(temp))
    
    FOR j=0,N_ELEMENTS(temp)-1 DO BEGIN
    	id_j  = WHERE(js_all EQ j)
    	wv_j  = wv_all[id_j]
    	gf_j  = gf_all[id_j] 
    	ab_j  = ab_all[id_j] 
    	el_j  = el_all[id_j] 
    	id  = WHERE(wv_j ge wmin and wv_j le wmax)
    	IF id[0] EQ -1 THEN GOTO,nolines
    	wv  = wv_j[id]
    	gf  = gf_j[id] 
    	ab  = ab_j[id] 
    	el  = el_j[id] 
    
;   *********************************************************
;   **** Step 3.0: Doppler broaden total emission   	 ****
;   ****           coefficients                          ****
;   *********************************************************
    
    	IF KEYWORD_SET(verbose)THEN PRINT,'Doppler broadening lines ... '
    	tion       = temp(j)
    	id_sort    = SORT(wv)
    	wv         = wv[id_sort]
    	ab         = ab[id_sort]
    	el         = el[id_sort]
    	gf         = gf[id_sort]>1e-30    
    	const      = 2.0*SQRT(2.0*alog(2.0))
    	FOR i=0,N_ELEMENTS(wv)-1 DO BEGIN
    	    doppler_fwhm = SQRT(tion / amss) * (1.09E-05 * wv(i))
	    sigma_dopp   = doppler_fwhm / const 
	    sigma_inst   = inst_fwhm / const
	    sigma_sum    = sigma_dopp^2+sigma_inst^2 
	    norm         = 1.0/(SQRT(2*!PI*sigma_sum))
	    tec[*,j]          = tec[*,j]+gf[i] * norm * EXP(-((wavelength-wv(i))^2/(2.0*sigma_sum)))
    	ENDFOR

    	gf  = gf * dens(j)
    	tec[*,j] = tec[*,j] * dens(j) ; Making units of /s

        IF KEYWORD_SET(doplot) AND j EQ 0 THEN BEGIN
    	    IF KEYWORD_SET(psplot)THEN PSPLOT,file='plot.ps' ELSE WINDOW,0,XS=700,YS=700
	    adas_colors,colors=colors
    	    !P.MULTI=[0,1,2]
	    IF ~KEYWORD_SET(yscale)THEN yscale = MAX(tec)
	    PLOT,wavelength,tec,YRANGE=[1e-10,YSCALE],XS=1,POS=[0.15,0.5,0.9,0.95],YLOG=YLOG,/NODATA,BACKGROUND=colors.white,COL=colors.black,CHARSIZE=2.0,YTITLE='Total Emission Coefficient / s!u-1!n'
	    OPLOT,wavelength,tec,COL=colors.black
	    PLOT,wavelength,tec,/NODATA,YR=[1E-6,1e2],XTITLE='Wavelength / Ang',/YLOG,BACKGROUND=colors.white,COL=colors.black,CHARSIZE=2.0,POS=[0.15,0.15,0.9,0.45],XSTYLE=9
	    coltable = [colors.red,colors.blue,colors.green,colors.gold,colors.orange,colors.aqua]
	    PRINT,'Plot legend: Z=0 - red; Z=1 - blue; Z=2 - green; Z=3 - gold; Z=4 - orange; Z>4 - aqua;' 
	    PRINT,'Plot legend: Solid lines 1st, 3rd, 5th etc. impurity; dashed lines 2nd, 4th, 6th etc. impurity;' 
	    FOR k=0,N_ELEMENTS(elem)-1 DO BEGIN
	    FOR j=0,MAX(ab) DO BEGIN
	    	id = WHERE(ab EQ j AND el EQ elem(k))
	    	IF (-1.0)^k EQ -1.0 THEN LINEST=5 ELSE LINEST=0
	    	IF j LT N_ELEMENTS(coltable)THEN ab_col = coltable(j) ELSE ab_col = coltable(N_ELEMENTS(coltable)-1)
	    	IF id[0] NE -1 THEN FOR i=0,N_ELEMENTS(wv[id])-1 DO OPLOT,[wv(id[i]),wv(id[i])],[1E-6,gf(id[i])],COL=ab_col,LINEST=linest
	    ENDFOR
	    ENDFOR
	    IF KEYWORD_SET(psplot)THEN BEGIN
		DEVICE,/CLOSE
		SET_PLOT,'X'
		!P.FONT=-1
		!P.THICK=1.5
		!X.THICK=1.5
		!Y.THICK=1.5
		!Z.THICK=1.5
	    ENDIF	
    	ENDIF
    ENDFOR
    tec = REFORM(tec)
    RETURN,{wavelength:wavelength, tec:tec,err:0,message:'Calculation complete.',wv:wv,gf:gf,ab:ab,el:el}
    NOLINES:
    npix = ROUND((wmax-wmin)/resol)    
    RETURN,{wavelength:wmin+((wmax-wmin)/(npix-1))*FINDGEN(npix), tec:FLTARR(npix),err:-1,message:'No lines in this spectral region.',wv:wv,gf:gf,ab:ab,el:el}
    
END


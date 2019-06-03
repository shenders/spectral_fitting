PRO user_psym,psym,fill=fill
;   ****************************
;   **** PSYM=1 : Circle    ****
;   **** PSYM=2 : Square    ****
;   **** PSYM=3 : Triangle  ****
;   **** PSYM=4 : Star	    ****
;   **** PSYM=5 : Diamond    ****
;   ****************************
    CASE psym OF
    	1 : Begin
    	    a = FINDGEN(17) * (!PI*2.0/16.0)
    	    USERSYM, cos(A), sin(A), FILL=FILL
        End
        2 : Begin
            USERSYM, [-1,-1,1, 1,-1], [-1, 1,1,-1,-1], FILL=FILL
        End
        3 : Begin
            USERSYM, [-1,0,1,-1], [-1,1,-1,-1], FILL=FILL
        End
        4 : Begin
    	    ang = (360. / 10 * findgen(11) + 90) / !RADEG  ;star angles every 36 deg
    	    r = ang*0
    	    r[2*indgen(6)] = 1.
    	    cp5 = cos(!pi/5.)
    	    r1 = 2. * cp5 - 1. / cp5
    	    r[2*indgen(5)+1] = r1
    	    r = r  / sqrt(!pi/4.) * 2. / (1.+r1)
    	    xarr = r * cos(ang)   &   yarr = r * sin(ang)
    	    USERSYM, xarr, yarr, FILL=FILL
        End
        5 : Begin
            USERSYM, [-1,0,1,0,-1], [0,1,0,-1,0], FILL=FILL
        End
    ENDCASE
    	
END
;+
; PROJECT:
;       ADAS
;
; NAME:
;       ADAS_COLORS
;
; PURPOSE: Another general purpose colour naming scheme.
;
;
; EXPLANATION: Returns a structure with colour names which should
;              work for all output devices, postscript and screens
;              of different colour depths.
;
;
; EXAMPLE:
;        adas_colors, colors=colors
;        plot, [3,4], [3,4], color=colors.red
;
; INPUTS:
;        None.
;
; OPTIONAL INPUTS:
;       None
;
; OUTPUTS:
;       The function returns a structure of name of colours. Currently,
;                  WHITE         DARKGRAY      LIGHTGRAY
;                  DIMGRAY       GRAY          ORANGE
;                  WHEAT         HOTPINK       GOLD
;                  PEACHPUFF     CHARCOAL      BEIGE
;                  SKY           ORCHID        AQUA
;                  PINK          NAVY          BLUE
;                  RED           GREEN         MAGENTA
;                  YELLOW        CYAN          BLACK
;
; OPTIONAL OUTPUTS:
;       None
;
; KEYWORD PARAMETERS:
;       help - returns this header as help
;
;
; WRITTEN:
;       Allan Whiterford
;
; MODIFIED:
;       1.1     Martin O'Mullane
;                - First release
;
; VERSION:
;       1.1     21-11-2007
;
;-
;-----------------------------------------------------------------------------

PRO adas_colors, colors = colors


usetrue=0

if !d.name eq 'PS' then begin
        usetrue=0
endif

if !d.name eq 'X' then begin
        device,get_decomposed=usetrue
        device,get_visual_depth=depth

        if depth lt 24 then begin
                device,decomposed=0
                device,get_decomposed=usetrue
        endif
endif

; RGB values from X11 rgb.txt.
; what about other devices? probably a z-buffer is the one to catch here

names  = ['Black']
rvalue = [  0    ]
gvalue = [  0    ]
bvalue = [  0    ]

names  = [names,   'Magenta', 'Cyan', 'Yellow', 'Green']
rvalue = [rvalue,     255,       0,      255,       0  ]
gvalue = [gvalue,       0,     255,      255,     255  ]
bvalue = [bvalue,     255,     255,        0,       0  ]

names  = [names,  'Red', 'Blue', 'Navy', 'Pink', 'Aqua']
rvalue = [rvalue,  255,     0,      0,    255,    112]
gvalue = [gvalue,    0,     0,      0,    127,    219]
bvalue = [bvalue,    0,   255,    115,    127,    147]

names  = [names,  'Orchid', 'Sky', 'Beige', 'Charcoal']
rvalue = [rvalue,   219,      0,     255,       80    ]
gvalue = [gvalue,   112,    163,     171,       80    ]
bvalue = [bvalue,   219,    255,     127,       80    ]

names  = [names,  'PeachPuff', 'Gold', 'HotPink', 'Wheat', 'Orange' ]
rvalue = [rvalue,   255,       255,     255,       255,      255    ]
gvalue = [gvalue,   218,       215,     105,       231,      165    ]
bvalue = [bvalue,   185,         0,     180,       186,        0    ]

names  = [names,   'Gray', 'DimGray', 'LightGray', 'DarkGray' ]
rvalue = [rvalue,    135,   105,       211,         169       ]
gvalue = [gvalue,    135,   105,       211,         169       ]
bvalue = [bvalue,    135,   105,       211,         169       ]

names  = [names,   'White']
rvalue = [rvalue,    255  ]
gvalue = [gvalue,    255  ]
bvalue = [bvalue,    255  ]


if usetrue eq 1 then begin
        value=create_struct('n_colours',n_elements(names))
        for i=0,n_elements(names)-1 do begin
                value=create_struct(names[i],rvalue[i]+256l*gvalue[i]+256l*256l*bvalue[i],value)
        end
endif else begin
        value=create_struct('n_colours',n_elements(names))
        offset=!d.table_size-n_elements(names)-2
        tvlct,rvalue,gvalue,bvalue,offset
        for i=0,n_elements(names)-1 do begin
                value=create_struct(names[i],i+offset,value)
        end
endelse

if arg_present(colors) then colors=value

END
Pro read_analysis,file,debug=debug

; Use ADAS routine to set the color

	adas_colors,colors=colors

; Open file

	openr,unit,file,/get_lun


; Read opening strings

	str     = ''
	for i=0,4 do readf,unit,str
	readf,unit,ntime
	time    = fltarr(ntime)

; Set arrays according to number of time points

	rawNii3995 = fltarr(ntime)
	rawNii4041 = fltarr(ntime)
	anlNii3995 = fltarr(ntime)
	anlNii4041 = fltarr(ntime)
	cn_mean = fltarr(ntime)
	cn_err  = fltarr(ntime)
	deltal  = fltarr(ntime)
	ne_mean = fltarr(ntime)
	ne_err  = fltarr(ntime)
	te_mean = fltarr(ntime)
	te_err  = fltarr(ntime)
	cN_flux = fltarr(ntime)
	tdiv    = fltarr(ntime)

; Read data

	readf,unit,str
	readf,unit,time
	readf,unit,str
	readf,unit,rawNii3995
	readf,unit,str
	readf,unit,rawNii4041
	readf,unit,str
	readf,unit,anlNii3995
	readf,unit,str
	readf,unit,anlNii4041
	readf,unit,str
	readf,unit,cn_flux
	readf,unit,str
	readf,unit,cn_mean
	readf,unit,cn_err
	readf,unit,str
	readf,unit,deltal
	readf,unit,str
	readf,unit,ne_mean
	readf,unit,ne_err
	readf,unit,str
	readf,unit,te_mean
	readf,unit,te_err
	readf,unit,str
	readf,unit,tdiv

; Plot data if requested

	if keyword_set(debug)then begin
		window,0,xs=1400,ys=900
		!p.multi=[0,3,2]
		!p.charsize=3
		user_psym,1,/fill & plot,time,rawNii3995>0,psym=8,col=colors.black,back=colors.white,xtitle='Time [s]',ytitle='N II radiance @ 399.5 nm [10!u19!n ph/s/m!u2!n/sr]',xr=xr,xs=1
		oplot,time,anlNii3995,col=colors.red
		plot,time,anlNii4041/anlNii3995,col=colors.black,back=colors.white,xtitle='Time [s]',ytitle='N II ratio 404.1/399.5 [-]',xr=xr,xs=1,yr=[0,0.4]
		plot,time,cn_mean,back=colors.white,col=colors.black,xtitle='Time [s]',ytitle='c!lN!n [%]; c!ln,flux-ratio!n [%]',xr=xr,xs=1,yr=[0,30]
		err_plot,time,cn_mean,cn_err,col=colors.black
		oplot,time,cn_flux,col=colors.red
		plot,time,deltal*100,back=colors.white,col=colors.black,xtitle='Time [s]',ytitle='Delta L [cm]; Tdiv [eV]',xr=xr,xs=1,yr=[0,40]
		oplot,time,tdiv,col=colors.red
		plot,time,ne_mean,back=colors.white,col=colors.black,xtitle='Time [s]',ytitle='n!le,N II!n [10!u20!n m!u-3!n]',yr=[0,3],xr=xr,xs=1
		err_plot,time,ne_mean,ne_err,col=colors.black
		plot,time,te_mean,back=colors.white,col=colors.black,xtitle='Time [s]',ytitle='T!le!n [eV]',xr=xr,xs=1,yr=[0,5]
		err_plot,time,te_mean,te_err,col=colors.black
	endif
End

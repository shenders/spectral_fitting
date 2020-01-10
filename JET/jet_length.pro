Function jet_length,tdiv,los,upperdl=upperdl,lowerdl=lowerdl,doplot=doplot,psplot=psplot

; For JET, tdiv refers to the ratio of N II emission at the strike-point to x-point
; JPN 85423
;SP1      2.54359
;SP2      2.56143
;SP3      2.57927
;SP4      2.59711
;SP5      2.61495
;SP6      2.63279
;SP7      2.65064
;SP8      2.66848
;SP9      2.68632
;SP10      2.70416
;SP11      2.72200
;SP12      2.73984
;SP13      2.75768
;SP14      2.77552
;SP15      2.79336
;SP16      2.81120
;SP17      2.82904
;SP18      2.84688
;SP19      2.86472
;SP20      2.88256
;SP21      2.90040
;SP22      2.91824

if keyword_set(doplot)then begin
	tdiv = findgen(100)*(2.0-0.5)/99.0+0.5
	nrow = 3
	ncol = 1
	setgraphics,xs=800,ys=1000,nrow=nrow,ncol=ncol,colors=colors,psplot=psplot,file='figures/dl_model_jet.ps'

	los  = 'SP10' 
	x = jet_length(tdiv,los,upperdl=upperdl,lowerdl=lowerdl)
	yr=[0,(1/0.85) * max(upperdl)*100]
	plot,tdiv,upperdl*100,/nodata,back=colors.white,col=colors.black,yr=yr,ys=1,ytitle='dL [cm]',position=graphpos(0,0,nrow,ncol,xspc=xspc,yspc=yspc),xtickname=replicate(' ',30)
	oband,tdiv,lowerdl*100,upperdl*100,/norm,col=colors.black
	oplot,tdiv,x*100,col=colors.red
    	legend,'R=2.70 m',colors.black,yshift=-0.6
	
	los  = 'SP14' 
	x = jet_length(tdiv,los,upperdl=upperdl,lowerdl=lowerdl)	
	plot,tdiv,upperdl*100,/nodata,back=colors.white,col=colors.black,yr=yr,ys=1,ytitle='dL [cm]',position=graphpos(0,1,nrow,ncol,xspc=xspc,yspc=yspc),xtickname=replicate(' ',30)
	oband,tdiv,lowerdl*100,upperdl*100,/norm,col=colors.black
	oplot,tdiv,x*100,col=colors.red
    	legend,'R=2.77 m',colors.black,yshift=-0.6

	los  = 'SP17' 
	x = jet_length(tdiv,los,upperdl=upperdl,lowerdl=lowerdl)
	plot,tdiv,upperdl*100,/nodata,back=colors.white,col=colors.black,yr=yr,ys=1,xtitle='N II X/S-point ratio [-]',ytitle='dL [cm]',position=graphpos(0,2,nrow,ncol,xspc=xspc,yspc=yspc)
	oband,tdiv,lowerdl*100,upperdl*100,/norm,col=colors.black
	oplot,tdiv,x*100,col=colors.red
    	legend,'R=2.82 m',colors.black,yshift=-0.6

	stop
endif
; Close to X-point
if los eq 'SP11' or los eq 'SP10' or los eq 'SP9' or los eq 'SP8' or $
   los eq 'SP7'  or los eq 'SP6'  or los eq 'SP5' or los eq 'SP4' or $
   los eq 'SP3'  or los eq 'SP2'  or los eq 'SP1' then begin
	Tdiv_ref = [0.5,0.6,1.0,1.5,2.0]
	func_up  = [8.0,6.5,5.5,5.5,5.5]/100
	func_dn  = [7.0,5.5,4.5,4.5,4.5]/100
	func_dl  = func_dn
	for i=0,n_elements(func_dl)-1 do func_dl[i] = (func_up[i]+func_dn[i])/2.0
	upperdl  = interpol(func_up ,tdiv_ref,tdiv)
	lowerdl  = interpol(func_dn ,tdiv_ref,tdiv)
	dl       = interpol(func_dl ,tdiv_ref,tdiv)
	return,dl
endif

; Mid way between x-point and strike-point
if los eq 'SP15' or los eq 'SP14' or los eq 'SP13' or los eq 'SP12' then begin
	Tdiv_ref = [0.5,0.6,1.0,1.5,2.0]
	func_up  = [6.5,8.0,6.5,5.5,5.5]/100
	func_dn  = [4.5,7.0,5.5,4.5,4.5]/100
	func_dl  = func_dn
	for i=0,n_elements(func_dl)-1 do func_dl[i] = (func_up[i]+func_dn[i])/2.0
	upperdl  = interpol(func_up ,tdiv_ref,tdiv)
	lowerdl  = interpol(func_dn ,tdiv_ref,tdiv)
	dl       = interpol(func_dl ,tdiv_ref,tdiv)
	return,dl
endif

; Close to strike-point
if los eq 'SP19' or los eq 'SP18' or los eq 'SP17' or los eq 'SP16' then begin
	Tdiv_ref = [0.5,0.6,1.0,1.5,2.0]
	func_up  = [2.5,5.0,6.5,8.0,8.0]/100
	func_dn  = [1.5,4.0,5.5,7.0,7.0]/100
	func_dl  = func_dn
	for i=0,n_elements(func_dl)-1 do func_dl[i] = (func_up[i]+func_dn[i])/2.0
	upperdl  = interpol(func_up ,tdiv_ref,tdiv)
	lowerdl  = interpol(func_dn ,tdiv_ref,tdiv)
	dl       = interpol(func_dl ,tdiv_ref,tdiv)
	return,dl
endif
return,-1
End

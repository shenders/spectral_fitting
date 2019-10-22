Function length,tdiv,shot,los,upperdl=upperdl,lowerdl=lowerdl,doplot=doplot,psplot=psplot

if keyword_set(doplot)then begin
	tdiv = findgen(100)*15/99.0
	shot = 30554
	los  = 'ROV014' 
	x = length(tdiv,shot,los,upperdl=upperdl,lowerdl=lowerdl)
	setgraphics,xs=800,ys=600,colors=colors,psplot=psplot,file='figures/dl_model.ps'
	yr=[0,(1/0.85) * max(upperdl)*100]
	
	plot,tdiv,upperdl*100,/nodata,back=colors.white,col=colors.black,yr=yr,ys=1,xtitle='Tdiv [eV]',ytitle='dL [cm]'
	oband,tdiv,lowerdl*100,upperdl*100,/norm,col=colors.black
	oplot,tdiv,x*100,col=colors.red
	stop
endif
if shot lt 33804 then begin
	if los eq 'ROV010' then los = 'ROV010'
	if los eq 'ROV011' then los = 'ROV011'
	if los eq 'ROV012' then los = 'ROV012'
	if los eq 'ROV013' then los = 'ROV013'
	if los eq 'ROV014' then los = 'ROV014'
endif
if shot lt 35056 and shot gt 33804 then begin
	if los eq 'ROV-09' then los = 'ROV010'
	if los eq 'ROV-10' then los = 'ROV011'
	if los eq 'ROV-11' then los = 'ROV012'
	if los eq 'ROV-12' then los = 'ROV013'
	if los eq 'ROV-13' then los = 'ROV014'
	if los eq 'ROV-14' then los = 'ROV015'
endif
if shot gt 35056 then begin
	if los eq 'ROV-10' then los = 'ROV010'
	if los eq 'ROV-11' then los = 'ROV011'
	if los eq 'ROV-12' then los = 'ROV012'
	if los eq 'ROV-13' then los = 'ROV013'
	if los eq 'ROV-14' then los = 'ROV014'
endif

if los eq 'ROV014' or los eq 'ROV013' then begin
; Only assessed for ROV014, ROV-13, ROV-014
	Tdiv_ref = [-5,0,3,12,20]
	if shot eq 30505 or shot eq 30506 then begin
		tdiv_ref = tdiv_ref -2.0	
	endif
	func_up  = [0.65,1.02,1.02,0.72,0.72]/10
	func_dn  = [0.48,0.85,0.85,0.55,0.55]/10
	func_dl  = func_dn
	for i=0,n_elements(func_dl)-1 do func_dl[i] = (func_up[i]+func_dn[i])/2.0
	upperdl  = interpol(func_up ,tdiv_ref,tdiv)
	lowerdl  = interpol(func_dn ,tdiv_ref,tdiv)
	dl       = interpol(func_dl ,tdiv_ref,tdiv)
	norm     = 1.0 / max(dl) * 0.08
	if shot eq 30505 or shot eq 30506 then begin
		
	endif
	dl       = dl*norm
	upperdl  = upperdl * norm
	lowerdl  = lowerdl * norm
	return,dl
endif
if los eq 'ROV012' then begin
; Only assessed for ROV014, ROV-13, ROV-014
	Tdiv_ref = [-5,1,5,20]
	func_up  = [0.2,0.5,1.05,0.85]
	norm     = func_up * 0.107 / max(func_up)
	upperdl  = interpol(func_up * norm,tdiv_ref,tdiv)
	Tdiv_ref = [-5,1,6,20]
	func_dn  = [0.1,0.2,0.9,0.7]
	lowerdl  = interpol(func_dn * norm,tdiv_ref,tdiv)
	dl       = fltarr(n_elements(tdiv))
	for i=0,n_elements(tdiv)-1 do dl[i] = (upperdl[i]+lowerdl[i])/2.0
	return,dl
endif
return,-1
End

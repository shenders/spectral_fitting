Function get_ratios,shot,los=los
	if ~keyword_set(los)then los='ROV014'
	restore,'save/'+string(shot,format='(i5)')+'/'+los+'-data.idl'
	read_signal_mrm,0L,shot,'MAC','Tdiv',x,y,2,exp=exp
	telm   = find_elm(shot,output.time)
	id     = where(telm ge 3.0)	
	nii3995= output.nii[id,0,0]
	nii4041= output.nii[id,0,1]
	nii4026= output.nii[id,0,2]
	time   = output.time[id]
	ratio1 = smooth(nii4041,10,/edge_tr)/smooth(nii3995,10,/edge_tr)
	ratio2 = smooth(nii4041,10,/edge_tr)/smooth(nii4026,10,/edge_tr)
	tdiv   = interpol(y,x,time)
	id     = where(tdiv gt 0 and tdiv le 20 and nii3995 gt 3e18)
	
return,{ratio1:ratio1[id],ratio2:ratio2[id]}
End

Pro ratios_db,los=los,showdens=showdens,psplot=psplot
!quiet=1
num    = 100
te_arr = [2.0,3.0,4.0,5.0,6.0]
te     = adas_vector(low=2.0,high=6,num=num)
dens_arr= [5e13,1e14,1.5e14,2.5e14,5e14]
dens   = adas_vector(low=1e13,high=1e16,num=num)
setgraphics,colors=colors,colpick=colpick,psplot=psplot,xs=800,ys=600,file='figures/ratios_db.ps'
user_psym,1,/fill
plot,[6,10],[0,0.5],/nodata,col=colors.black,back=colors.white,$
  xtitle='N II 404.1/402.6 [-]',ytitle='N II 404.1/399.5 [-]'

if keyword_set(showdens)then begin
	for i=0,n_elements(dens_arr)-1 do begin
		dens = fltarr(num)+dens_arr[i]
		atomdb,te,dens,tec3995=tec3995,$
		           tec4041=tec4041,$
	       		   tec4026=tec4026
		ratio1 = tec4041/tec3995
		ratio2 = tec4041/tec4026
		oplot,ratio2,ratio1,col=colpick[i]
	endfor
endif else begin
	for i=0,n_elements(te_arr)-1 do begin
		te = fltarr(num)+te_arr[i]
		atomdb,te,dens,tec3995=tec3995,$
		           tec4041=tec4041,$
	       		   tec4026=tec4026
		ratio1 = tec4041/tec3995
		ratio2 = tec4041/tec4026
		oplot,ratio2,ratio1,col=colpick[i]
	endfor
end
shot_db = [30506,30505,30776,30298,30306,30307,30554,32244,30623]
if ~keyword_set(los)then los = ['ROV014','ROV012','ROV010']
collos  = [colors.blue,colors.green,colors.red]
for j=0,n_elements(los)-1 do begin
	for i=0,n_elements(shot_db)-1 do begin
		x=get_ratios(shot_db[i],los=los[j])
		oplot,x.ratio2,x.ratio1,psym=8,col=colors.black
	endfor
end
for j=0,n_elements(los)-1 do begin
	r1 = -1.0
	r2 = -1.0
	for i=0,n_elements(shot_db)-1 do begin
		x=get_ratios(shot_db[i],los=los[j])
		r1 = [r1,x.ratio1]
		r2 = [r2,x.ratio2]
	endfor
	r1 = r1[1:*] & r2 = r2[1:*]
	m1 = moment(r1) & m2 = moment(r2)
	oplot,[-10,m2[0]],[-10,m1[0]],psym=8,col=collos[j]
	errors,[-10,m2[0]],[-10,m1[0]],xstd=[-1.0,sqrt(m2[1])],ystd=[-1.0,sqrt(m1[1])],col=collos[j]
	print,'Ratio 1: ',m1[0],' +/- ',sqrt(m1[1]) & print,'Ratio 2: ',m2[0],' +/- ',sqrt(m2[1])
end
stop

End

Pro flux_comparison,val,diff,debug=debug,extra=extra
val = -1
diff= -1

if keyword_set(extra)then begin
	
	a = get_nii(30623,$
	  los='ROV014',$
	  xr=[0.0,4.0],$
	  /interelm,$
	  append='start',$
	  sig3995=sig3995)
	diff = [diff, 12.49/12.3]
	t0_a = 3.5
	id   = where(a.rawtime ge t0_a-0.5 and a.rawtime le t0_a)
	val  = [val,mean(a.rawnii3995[id])]
	
	a = get_nii(32932,$
	  los='ROV014',$
	  xr=[0.0,4.0],$
	  /interelm,$
	  append='start',$
	  sig3995=sig3995)
	diff = [diff, 2.26/13.0]
	t0_a = 3.0
	id   = where(a.rawtime ge t0_a-0.5 and a.rawtime le t0_a)
	val  = [val,mean(a.rawnii3995[id])]

	a = get_nii(30506,$
	  los='ROV014',$
	  xr=[0.0,4.0],$
	  /interelm,$
	  append='start',$
	  sig3995=sig3995)
	diff = [diff, 6.17/7.3]
	t0_a = 2.15
	id   = where(a.rawtime ge t0_a-0.5 and a.rawtime le t0_a)
	val  = [val,mean(a.rawnii3995[id])]

	a = get_nii(32273,$
	  los='ROV014',$
	  xr=[0.0,4.0],$
	  /interelm,$
	  append='start',$
	  sig3995=sig3995)
	diff = [diff, 3.57/7.3]
	t0_a = 2.10
	id   = where(a.rawtime ge t0_a-0.5 and a.rawtime le t0_a)
	val  = [val,mean(a.rawnii3995[id])]
	if ~keyword_set(debug)then begin
		val=val[1:*]
		diff=diff[1:*]
		return
	endif
endif

a = get_nii(32244,$
	  los='ROV014',$
	  xr=[0.0,4.0],$
	  /interelm,$
	  append='start',$
	  sig3995=sig3995)
diff = [diff, 5.78/10.3]
t0_a = 2.08
id   = where(a.rawtime ge t0_a-0.5 and a.rawtime le t0_a)
val  = [val,mean(a.rawnii3995[id])]

i = get_nii(30306,$
	  los='ROV014',$
	  xr=[0.0,4.0],$
	  /interelm,$
	  append='start',$
	  sig3995=sig3995)
diff = [diff, 5.78/10.3]
t0_i = 3.5
id   = where(i.rawtime ge t0_i-0.5 and i.rawtime le t0_i)
val  = [val,mean(i.rawnii3995[id])]

v = get_nii(30554,$
	  los='ROV014',$
	  xr=[0.0,4.0],$
	  /interelm,$
	  append='start',$
	  sig3995=sig3995)
diff = [diff, 7.76 / 9.76]
t0_v = 3.48
id   = where(v.rawtime ge t0_v-0.5 and v.rawtime le t0_v)
val  = [val,mean(v.rawnii3995[id])]

b = get_nii(30776,$
	  los='ROV014',$
	  xr=[0.0,4.0],$
	  /interelm,$
	  append='start',$
	  sig3995=sig3995)
diff = [diff, 9.3/12.46]
t0_b = 1.92
id   = where(b.rawtime ge t0_b-0.5 and b.rawtime le t0_b)
val  = [val,mean(b.rawnii3995[id])]

z = get_nii(30505,$
	  los='ROV014',$
	  xr=[0.0,4.0],$
	  /interelm,$
 	  channel=22,$
	  append='start',$
	  sig3995=sig3995)
diff = [diff, 4.78/4.43]
t0_z = 2.01
id   = where(z.rawtime ge t0_z-0.5 and z.rawtime le t0_z)
val  = [val,mean(z.rawnii3995[id])]

x = get_nii(35158,$
	  los='ROV-14',$
	  xr=[0.0,4.0],$
	  /interelm,$
 	  channel=22,$
	  /use_evl,$
	  append=append,$
	  diag='FVL',$
	  sig3995=sig3995)
diff = [diff,9.22/12.29]
t0_x = 1.62
id   = where(x.rawtime ge t0_x-0.5 and x.rawtime le t0_x)
val  = [val,mean(x.rawnii3995[id])]

y = get_nii(34971,$
	  los='ROV-13',$
	  xr=[0.0,4.0],$
	  /interelm,$
 	  channel=19,$
	  /use_evl,$
	  append=append,$
	  diag='EVL',$
	  sig3995=sig3995)
diff = [diff,1.571/13.93]
t0_y = 3.01
id   = where(y.rawtime ge t0_y-0.5 and y.rawtime le t0_y)
val  = [val,mean(y.rawnii3995[id])]

u = get_nii(34973,$
	  los='ROV-13',$
	  xr=[0.0,4.0],$
	  /interelm,$
 	  channel=19,$
	  /use_evl,$
	  append=append,$
	  diag='EVL',$
	  sig3995=sig3995)
diff = [diff,2.27/15.10]
t0_u = 3.01
id   = where(u.rawtime ge t0_u-0.5 and u.rawtime le t0_u)
val  = [val,mean(u.rawnii3995[id])]

n = get_nii(34358,$
	  los='ROV-13',$
	  xr=[0.0,4.0],$
	  /interelm,$
 	  channel=19,$
	  /use_evl,$
	  append=append,$
	  diag='EVL',$
	  sig3995=sig3995)
diff = [diff,1.35/7.35]
t0_n = 3.5
id   = where(n.rawtime ge t0_n-0.5 and n.rawtime le t0_n)
val  = [val,mean(n.rawnii3995[id])]

c = get_nii(34359,$
	  los='ROV-13',$
	  xr=[0.0,4.0],$
	  /interelm,$
 	  channel=19,$
	  /use_evl,$
	  append=append,$
	  diag='EVL',$
	  sig3995=sig3995)
diff = [diff,2.02/7.5]
t0_c = 3.5
id   = where(c.rawtime ge t0_c-0.5 and c.rawtime le t0_c)
val  = [val,mean(c.rawnii3995[id])]

m = get_nii(34368,$
	  los='ROV-12',$
	  xr=[0.0,4.0],$
	  /interelm,$
 	  channel=19,$
	  /use_evl,$
	  append=append,$
	  diag='EVL',$
	  sig3995=sig3995)
diff = [diff,2.15/7.31]
t0_m = 3.49
id   = where(m.rawtime ge t0_m-0.5 and m.rawtime le t0_m)
val  = [val,mean(m.rawnii3995[id])]


val = val[1:*]
diff= diff[1:*]
if keyword_set(debug)then begin
	setgraphics,xs=600,ys=400,colors=colors,nrow=1,ncol=1
	user_psym,5,/fill & plot,val/1e17,diff,col=colors.black,back=colors.white,psym=8,/nodata,$
	xtitle='N II @399.5nm [10!u17!n ph/s/m!u2!n/sr]',ytitle='c!lN!n spec/flux'
	if keyword_set(extra)then begin
		xval = findgen(100)/99.0*20.0
		yval = 0.1+0.09*xval
		oband,xval,fltarr(100)+0.8,fltarr(100)+1.2,col=colors.blue,/norm
		oplot,xval,yval,linest=5,col=colors.black
		oplot,[10,10],[0,10],linest=5,col=colors.black
	endif
	oplot,val/1e17,diff,col=colors.black,psym=8,symsize=2.0
endif
end

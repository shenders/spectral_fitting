Pro los,tau=tau

num     = 100
dist    = findgen(num)
setgraphics,colors=colors,xs=800,ys=900,nrow=3,ncol=2
plot,[5,11],[0,1.2],/nodata,col=colors.black,back=colors.white,xs=1,ys=1,xtitle='N II: 404.1/402.6 nm',ytitle='N II: 404.1/399.5 nm'
for i=2,8 do begin
	te   = i + fltarr(num)
	dens = adas_vector(high=1e15,low=5e13,num=num)
	atomdb,te,dens,tec3995=tec3995,$
		  tec4041=tec4041,$
		  tec4026=tec4026,$
		  tau=tau
	oplot,tec4041/tec4026,tec4041/tec3995,col=colors.black
endfor

te_grad = adas_vector(high=10,low=1,num=num)
dens    = fltarr(num)+1e14

atomdb,te_Grad,dens,tec3995=tec3995,$
		   tec4041=tec4041,$
		   tec4026=tec4026,$
		   tau=tau
plot,te_Grad,tec3995,col=colors.black,back=colors.white,title='Exc + rec',yr=[0,max(tec3995)>max(tec4041)],xtitle='T!le!n [eV]',ytitle='TEC [cm!u3!n/s]'
oplot,te_Grad,tec4041,col=colors.red

plot,[5,11],[0,0.8],/nodata,col=colors.black,back=colors.white,xs=1,ys=1,xtitle='N II: 404.1/402.6 nm',ytitle='N II: 404.1/399.5 nm'
for i=2,8 do begin
	te   = i + fltarr(num)
	dens = adas_vector(high=1e15,low=5e13,num=num)
	atomdb,te,dens,tec3995=tec3995,$
		  tec4041=tec4041,$
		  tec4026=tec4026,$
		  tau=tau,/exc_only
	oplot,tec4041/tec4026,tec4041/tec3995,col=colors.black
endfor

te_grad = adas_vector(high=10,low=1,num=num)
dens    = fltarr(num)+1e14
atomdb,te_Grad,dens,tec3995=tec3995,$
		   tec4041=tec4041,$
		   tec4026=tec4026,$
		   tau=tau,/exc_only

plot,te_grad,tec3995,col=colors.black,back=colors.white,title='Exc only',yr=[0,max(tec3995)>max(tec4041)],xtitle='T!le!n [eV]',ytitle='TEC [cm!u3!n/s]'
oplot,te_grad,tec4041,col=colors.red
plot,[5,11],[0,4],/nodata,col=colors.black,back=colors.white,xs=1,ys=1,xtitle='N II: 404.1/402.6 nm',ytitle='N II: 404.1/399.5 nm'
for i=2,8 do begin
	te   = i + fltarr(num)
	dens = adas_vector(high=1e15,low=5e13,num=num)
	atomdb,te,dens,tec3995=tec3995,$
		  tec4041=tec4041,$
		  tec4026=tec4026,$
		  tau=tau,/rec_only
	oplot,tec4041/tec4026,tec4041/tec3995,col=colors.black
endfor
atomdb,te_Grad,dens,tec3995=tec3995,$
		   tec4041=tec4041,$
		   tec4026=tec4026,$
		   tau=tau,/rec_only
te_grad = adas_vector(high=10,low=1,num=num)
dens    = fltarr(num)+1e14
plot,te_grad,tec3995,col=colors.black,back=colors.white,title='Rec only',yr=[0,max(tec3995)>max(tec4041)],xtitle='T!le!n [eV]',ytitle='TEC [cm!u3!n/s]'
oplot,te_grad,tec4041,col=colors.red


stop
end

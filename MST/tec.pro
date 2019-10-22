Pro tec,psplot=psplot

num  = 100
te   = adas_vector(low=1,high=15,num=num,/linear)
dens = fltarr(num)+1e14

atomdb,te,dens,tec3995=tec3995,$
	       tec4041=tec4041,$
	       tec4026=tec4026

atomdb,te,dens,tec3995=tec3995_tr,$
	       tec4041=tec4041_tr,$
	       tec4026=tec4026_tr,tau=0.1e-3

setgraphics,colors=colors,psplot=psplot,file='tec.ps',xs=800,ys=600
plot,te,tec3995/max(tec3995),background=colors.white,col=colors.black,$
 xtitle='T!le!n [eV]',ytitle='N II TEC [-]'
oplot,te,tec3995_tr/max(tec3995_tr),col=colors.red
xyouts,[11,11],[0.8,0.7],['tau=inf','tau=0.1 ms'],col=[colors.black,colors.red]
stop
end


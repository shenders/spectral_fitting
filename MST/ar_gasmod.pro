Pro ar_gasmod,psplot=psplot

channel = 0

read_signal_mrm,0L,35157,'GVL','Ar0_7504',t_ar1,ar1,2
read_signal_mrm,0L,35158,'GVL','Ar0_7504',t_ar2,ar2,2
read_signal_mrm,0L,35157,'UVS','CFA03A',t_arvalve1,arvalve1,2
read_signal_mrm,0L,35158,'UVS','CFA03A',t_arvalve2,arvalve2,2

arvalve1 = arvalve1 / 1e21
arvalve2 = arvalve2 / 1e21
ar1 = ar1 / 1e16
ar2 = ar2 / 1e16

setgraphics,xs=800,ys=1000,nrow=1,ncol=2,colors=colors,psplot=psplot,file='figures/ar_valve.ps'
xr = [2.0,7.0]
yr = [0.0,5.0]
plot,t_arvalve1,arvalve1,col=colors.black,back=colors.white,xr=xr,xs=1,/nodata,yr=yr,ytitle='[1E21 #/s, 1E16 ph/s/m!u2!n/sr]'
oplot,t_arvalve1,arvalve1,col=colors.black
oplot,t_ar1,ar1,col=colors.blue
xyouts,[2.3,2.3,2.3],[4.5,4.1,3.7],['AUG #35157','Ar CFA03A valve rate','Divertor Ar I 750.4 nm'],col=[colors.black,colors.black,colors.black]
oplot,[2.1,2.2],[4.15,4.15],col=colors.black
oplot,[2.1,2.2],[3.75,3.75],col=colors.blue

plot,t_arvalve2,arvalve2,col=colors.black,back=colors.white,xr=xr,xs=1,/nodata,yr=yr,ytitle='[1E21 #/s, 1E16 ph/s/m!u2!n/sr]',xtitle='Time [s]'
oplot,t_arvalve2,arvalve2,col=colors.black
oplot,t_ar2,ar2,col=colors.blue
xyouts,[2.3,2.3,2.3],[4.5,4.1,3.7],['AUG #35158','Ar CFA03A valve rate','Divertor Ar I 750.4 nm'],col=[colors.black,colors.black,colors.black]
oplot,[2.1,2.2],[4.15,4.15],col=colors.black
oplot,[2.1,2.2],[3.75,3.75],col=colors.blue
setgraphics,/close,psplot=psplot

Stop
End

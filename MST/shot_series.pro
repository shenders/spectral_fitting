Pro shot_series	

channel = 19
use_evl = 1
elmcond = 3.0
sv      = 20
lowerte = 3.1
upperte = 3.6
debug   =1
x280 = get_nii(34280,los='ROV-13',$
	       xr=[3.5,5],$
	       /interelm,$
 	       channel=channel,$
	       use_evl=use_evl,$
	       append=append,$
	       diag=diag,$
	       elmcond=elmcond)

x281 = get_nii(34281,los='ROV-13',$
	       xr=[3.5,5],$
	       /interelm,$
 	       channel=channel,$
	       use_evl=use_evl,$
	       append=append,$
	       diag=diag,$
	       elmcond=elmcond)

cn280 = estimate_ratios(x280,34280,'ROV-13',1.0/0.67,sm=sv,debug=debug,lowerte=lowerte,upperte=upperte)
cn281 = estimate_ratios(x281,34281,'ROV-13',1.0/0.67,sm=sv,debug=debug,lowerte=lowerte,upperte=upperte)

setgraphics,nrow=2,ncol=1,xs=600,ys=1000,colors=colors 

user_psym,5 & plot,x280.rawtime,x280.rawnii3995/1e17,back=colors.white,col=colors.black,psym=8,yr=[0,15] 
user_psym,5,/fill & oplot,x280.time,x280.nii3995/1e17,col=colors.black,psym=8
oplot,x280.time,smooth(x280.nii3995/1e17,sv,/edge_truncate),col=colors.red
oplot,x280.time,x280.tdiv,col=colors.red

user_psym,5 & plot,x281.rawtime,x281.rawnii3995/1e17,back=colors.white,col=colors.black,psym=8,yr=[0,40] 
user_psym,5,/fill & oplot,x281.time,x281.nii3995/1e17,col=colors.black,psym=8
oplot,x281.time,smooth(x281.nii3995/1e17,sv,/edge_truncate),col=colors.red
oplot,x281.time,x281.tdiv,col=colors.red

setgraphics,nrow=2,ncol=1,xs=600,ys=1000,colors=colors 

plot,x280.tdiv,cn280.cn_upper*100,/nodata,back=colors.white,col=colors.black,yr=[0,5]
oband,x280.tdiv,cn280.cn_upper*100,cn280.cn_lower*100,/norm,col=colors.black
 
plot,x281.tdiv,cn281.cn_upper*100,/nodata,back=colors.white,col=colors.black,yr=[0,5]
oband,x281.tdiv,cn281.cn_upper*100,cn281.cn_lower*100,/norm,col=colors.black
stop

End

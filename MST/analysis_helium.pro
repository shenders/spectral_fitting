Pro analysis_helium,shot,psplot=psplot
if shot eq 36708 then begin

read_signal_mrm,0L,shot,'EVL','N_1_3995',t1,td1,2
read_signal_mrm,0L,shot,'EVL','Ne',t1,stark,2
time = t1
x14=td1[*,8]
x12=td1[*,21]
x10=td1[*,20]
x9=td1[*,12]
s14=stark[*,8]
s12=stark[*,21]
s10=stark[*,20]
s9=stark[*,12]

setgraphics,colors=colors,xs=800,ys=900,nrow=2,ncol=1,psplot=psplot,file='He_Hmode_36708.ps'
plot,time,x12/1e18,yr=[0,10],/nodata,xs=1,col=colors.black,back=colors.white,$
 xtitle='Time [s]',ytitle='N II [10!u18!n ph/s/m!u2!n/sr]',xr=[2.0,5.1]
oplot,time,x9/1e18,col=colors.orange
oplot,time,x10/1e18,col=colors.red
oplot,time,x14/1e18,col=colors.green
oplot,time,x12/1e18,col=colors.blue
x12=fetch_data(36708,'ROV-12',machine='AUG',tr=[2.7,3.5],diag='EVS',species=['He','C','W','O','N','D'],res=0.1)
x10=fetch_data(36708,'ROV-10',machine='AUG',tr=[2.7,3.5],diag='EVS',species=['He','C','W','O','N','D'],res=0.1)
x9=fetch_data(36708,'ROV-09',machine='AUG',tr=[2.7,3.5],diag='EVS',species=['He','C','W','O','N','D'],res=0.1)
oplot,[2.8,2.8],[0,1e20],linest=5,col=colors.black
oplot,[2.9,2.9],[0,1e20],linest=5,col=colors.black
oplot,[3.4,3.4],[0,1e20],linest=5,col=colors.black
oplot,[3.5,3.5],[0,1e20],linest=5,col=colors.black

sm=4


ratio_x9_1 = smooth(x9.nii[*,0,1],sm,/edge_trun)/smooth(x9.nii[*,0,2],sm,/edge_trun)
ratio_x9_2 = smooth(x9.nii[*,0,1],sm,/edge_trun)/smooth(x9.nii[*,0,0],sm,/edge_trun)
ratio_x10_1 = smooth(x10.nii[*,0,1],sm,/edge_trun)/smooth(x10.nii[*,0,2],sm,/edge_trun)
ratio_x10_2 = smooth(x10.nii[*,0,1],sm,/edge_trun)/smooth(x10.nii[*,0,0],sm,/edge_trun)
ratio_x12_1 = smooth(x12.nii[*,0,1],sm,/edge_trun)/smooth(x12.nii[*,0,2],sm,/edge_trun)
ratio_x12_2 = smooth(x12.nii[*,0,1],sm,/edge_trun)/smooth(x12.nii[*,0,0],sm,/edge_trun)


num    = 100
te_arr = [1.0,2.0,3.0,4.0,5.0,6.0]
te     = adas_vector(low=1.0,high=6,num=num)
dens_arr= [5e13,1e14,1.5e14,2.5e14,5e14,1e15]
dens   = adas_vector(low=1e13,high=1e16,num=num)

user_psym,1,/fill
plot,[5,11],[0,0.7],/nodata,col=colors.black,back=colors.white,$
  xtitle='N II 404.1/402.6 [-]',ytitle='N II 404.1/399.5 [-]',xs=1,ys=1

for i=0,n_elements(te_arr)-1 do begin
	te = fltarr(num)+te_arr[i]
	atomdb,te,dens,tec3995=tec3995,$
	           tec4041=tec4041,$
       		   tec4026=tec4026
	ratio1 = tec4041/tec3995
	ratio2 = tec4041/tec4026
	id = where(ratio2 le 10)
	oplot,ratio2[id],ratio1[id],col=colors.black
endfor
num    = 100
te_arr = [1.0,2.0,3.0,4.0,5.0,6.0]
te     = adas_vector(low=1.0,high=6,num=num)
dens_arr= [5e13,1e14,2e14,1e15]
dens   = adas_vector(low=1e13,high=1e16,num=num)
for i=0,n_elements(dens_arr)-1 do begin
	dens = fltarr(num)+dens_arr[i]
	atomdb,te,dens,tec3995=tec3995,$
	           tec4041=tec4041,$
       		   tec4026=tec4026
	ratio1 = tec4041/tec3995
	ratio2 = tec4041/tec4026
	id = where(ratio2 le 10)
	oplot,ratio2[id],ratio1[id],col=colors.black,linest=5
endfor


id = where(x12.time gt 2.8 and x12.time le 2.9)

x12_mean_x1 = mean(ratio_x12_1[id])<10
x12_mean_y1 = mean(ratio_x12_2[id])
x10_mean_x1 = mean(ratio_x10_1[id])<10
x10_mean_y1 = mean(ratio_x10_2[id])
x9_mean_x1 = mean(ratio_x9_1[id])<10
x9_mean_y1 = mean(ratio_x9_2[id])
id = where(x12.time gt 3.4 and x12.time le 3.5)
x12_mean_x2 = mean(ratio_x12_1[id])<10
x12_mean_y2 = mean(ratio_x12_2[id])
x10_mean_x2 = mean(ratio_x10_1[id])<10
x10_mean_y2 = mean(ratio_x10_2[id])
x9_mean_x2 = mean(ratio_x9_1[id])<10
x9_mean_y2 = mean(ratio_x9_2[id])

user_psym,4,/fill & oplot,[x12_mean_x1,-1],[x12_mean_y1,-1],psym=8,col=colors.blue,symsize=2.5
user_psym,4,/fill & oplot,[-1,x12_mean_x2],[-1,x12_mean_y2],psym=8,col=colors.blue,symsize=2.5
user_psym,4 & oplot,[-1,x12_mean_x2],[-1,x12_mean_y2],psym=8,col=colors.black,symsize=2.5
;user_psym,4,/fill & oplot,[x10_mean_x1,-1],[x10_mean_y1,-1],psym=8,col=colors.red,symsize=2.5
;user_psym,4,/fill & oplot,[-1,x10_mean_x2],[-1,x10_mean_y2],psym=8,col=colors.red,symsize=2.5
;user_psym,4 & oplot,[-1,x10_mean_x2],[-1,x10_mean_y2],psym=8,col=colors.black,symsize=2.5
user_psym,4,/fill & oplot,[x9_mean_x1,-1],[x9_mean_y1,-1],psym=8,col=colors.orange,symsize=2.5
user_psym,4,/fill & oplot,[-1,x9_mean_x2],[-1,x9_mean_y2],psym=8,col=colors.orange,symsize=2.5
user_psym,4 & oplot,[-1,x9_mean_x2],[-1,x9_mean_y2],psym=8,col=colors.black,symsize=2.5

endif
if shot eq 36702 then begin
read_signal_mrm,0L,shot,'EVL','He0_4026',t1,td1,2
read_signal_mrm,0L,shot,'EVL','He0_3965',t1,td2,2
time = t1

x14=td1[*,8]
x12=td1[*,21]
x10=td1[*,20]
x9=td1[*,12]
y14=td2[*,8]
y12=td2[*,21]
y10=td2[*,20]
y9=td2[*,12]
setgraphics,colors=colors,xs=800,ys=900,nrow=2,ncol=1,psplot=psplot,file='He_Lmode_36702.ps'
plot,time,x12/1e17,yr=[0,4],/nodata,xs=1,col=colors.black,back=colors.white,$
 xtitle='Time [s]',ytitle='He I [10!u17!n ph/s/m!u2!n/sr]',xr=[2.5,4.2]
oplot,time,x14/1e17,col=colors.green
oplot,time,x12/1e17,col=colors.blue
oplot,time,x10/1e17,col=colors.red
oplot,time,x9/1e17,col=colors.orange
sm=4
plot,time,x12/y12,/nodata,xs=1,col=colors.black,back=colors.white,$
 xtitle='Time [s]',ytitle='He I 402/396',xr=[2.5,4.2],yr=[0,6]
oplot,time,smooth(x14,sm,/edge_trun)/smooth(y14,sm,/edge_trun),col=colors.green
oplot,time,smooth(x12,sm,/edge_trun)/smooth(y12,sm,/edge_trun),col=colors.blue
oplot,time,smooth(x10,sm,/edge_trun)/smooth(y10,sm,/edge_trun),col=colors.red
oplot,time,smooth(x9,sm,/edge_trun)/smooth(y9,sm,/edge_trun),col=colors.orange
stop
endif

if shot eq 36703 then begin

read_signal_mrm,0L,shot,'EVL','N_1_3995',t1,td1,2
time = t1
x14=td1[*,8]
x12=td1[*,21]
x10=td1[*,20]
x9=td1[*,12]

setgraphics,colors=colors,xs=800,ys=900,nrow=2,ncol=1,psplot=psplot,file='He_Lmode_36703.ps'
plot,time,x12/1e18,yr=[0,4],/nodata,xs=1,col=colors.black,back=colors.white,$
 xtitle='Time [s]',ytitle='N II [10!u18!n ph/s/m!u2!n/sr]',xr=[3.5,5.1]
oplot,time,x14/1e18,col=colors.green
oplot,time,x12/1e18,col=colors.blue
oplot,time,x10/1e18,col=colors.red
oplot,time,x9/1e18,col=colors.orange

oplot,[3.7,3.7],[0,1e20],linest=5,col=colors.black
oplot,[3.8,3.8],[0,1e20],linest=5,col=colors.black
oplot,[4.4,4.4],[0,1e20],linest=5,col=colors.black
oplot,[4.5,4.5],[0,1e20],linest=5,col=colors.black
oplot,[4.9,4.9],[0,1e20],linest=5,col=colors.black
oplot,[5.0,5.0],[0,1e20],linest=5,col=colors.black

sm=2
x10=fetch_data(36703,'ROV-10',tr=[3.6,5.1],machine='AUG',diag='EVS',species=['He','C','W','O','N'],res=0.1,/save)
x9=fetch_data(36703,'ROV-09',tr=[3.6,5.1],machine='AUG',diag='EVS',species=['He','C','W','O','N'],res=0.1,/save)
ratio_x9_1 = smooth(x9.nii[*,0,1],sm,/edge_trun)/smooth(x9.nii[*,0,2],sm,/edge_trun)
ratio_x9_2 = smooth(x9.nii[*,0,1],sm,/edge_trun)/smooth(x9.nii[*,0,0],sm,/edge_trun)
ratio_x10_1 = smooth(x10.nii[*,0,1],sm,/edge_trun)/smooth(x10.nii[*,0,2],sm,/edge_trun)
ratio_x10_2 = smooth(x10.nii[*,0,1],sm,/edge_trun)/smooth(x10.nii[*,0,0],sm,/edge_trun)

num    = 100
te_arr = [1.0,2.0,3.0,4.0,5.0,6.0]
te     = adas_vector(low=1.0,high=6,num=num)
dens_arr= [5e13,1e14,1.5e14,2.5e14,5e14,1e15]
dens   = adas_vector(low=1e13,high=1e16,num=num)

user_psym,1,/fill
plot,[5,11],[0,0.8],/nodata,col=colors.black,back=colors.white,$
  xtitle='N II 404.1/402.6 [-]',ytitle='N II 404.1/399.5 [-]',xs=1,ys=1

for i=0,n_elements(te_arr)-1 do begin
	te = fltarr(num)+te_arr[i]
	atomdb,te,dens,tec3995=tec3995,$
	           tec4041=tec4041,$
       		   tec4026=tec4026
	ratio1 = tec4041/tec3995
	ratio2 = tec4041/tec4026
	id = where(ratio2 le 10)
	oplot,ratio2[id],ratio1[id],col=colors.black
endfor
num    = 100
te_arr = [1.0,2.0,3.0,4.0,5.0,6.0]
te     = adas_vector(low=1.0,high=6,num=num)
dens_arr= [5e13,1e14,2e14,1e15]
dens   = adas_vector(low=1e13,high=1e16,num=num)
for i=0,n_elements(dens_arr)-1 do begin
	dens = fltarr(num)+dens_arr[i]
	atomdb,te,dens,tec3995=tec3995,$
	           tec4041=tec4041,$
       		   tec4026=tec4026
	ratio1 = tec4041/tec3995
	ratio2 = tec4041/tec4026
	id = where(ratio2 le 10)
	oplot,ratio2[id],ratio1[id],col=colors.black,linest=5
endfor


id = where(x10.time gt 3.7 and x10.time le 3.8)

x9_mean_x1 = mean(ratio_x9_1[id])<10
x9_mean_y1 = mean(ratio_x9_2[id])
x10_mean_x1 = mean(ratio_x10_1[id])<10
x10_mean_y1 = mean(ratio_x10_2[id])
id = where(x10.time gt 4.4 and x10.time le 4.5)
x9_mean_x2 = mean(ratio_x9_1[id])<10
x9_mean_y2 = mean(ratio_x9_2[id])
x10_mean_x2 = mean(ratio_x10_1[id])<10
x10_mean_y2 = mean(ratio_x10_2[id])
id = where(x10.time gt 4.7 and x10.time le 4.8)
x9_mean_x3 = mean(ratio_x9_1[id])<10
x9_mean_y3 = mean(ratio_x9_2[id])
x10_mean_x3 = mean(ratio_x10_1[id])<10
x10_mean_y3 = mean(ratio_x10_2[id])

user_psym,4,/fill & oplot,[x9_mean_x1,x9_mean_x2,x9_mean_x3],[x9_mean_y1,x9_mean_y2,x9_mean_y3],psym=8,col=colors.orange,symsize=2.5
user_psym,4,/fill & oplot,[x10_mean_x1,x10_mean_x2,x10_mean_x3],[x10_mean_y1,x10_mean_y2,x10_mean_y3],psym=8,col=colors.red,symsize=2.5

endif
Stop
End

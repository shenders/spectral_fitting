Pro plot_legend,text,col=col
	y0 = !y.crange[0] + (!y.crange[1]-!y.crange[0]) *0.8
	y1 = !y.crange[0] + (!y.crange[1]-!y.crange[0]) *0.98
	x0 = !x.crange[0] + (!x.crange[1]-!x.crange[0]) *0.7
	x1 = !x.crange[0] + (!x.crange[1]-!x.crange[0]) *0.98
	xtext = x0 + (x1-x0)*0.05
	ytext = y0 + (y1-y0)*0.4
	oplot,[x0,x1,x1,x0,x0],[y0,y0,y1,y1,y0],col=col,thick=!x.thick
	xyouts,[xtext],[ytext],[text],col=col,charsize=!x.charsize
End

Pro power_comp,psplot=psplot,power=power

if ~keyword_set(power)then power='med'
setgraphics,xs=1200,ys=1000,ncol=2,nrow=3,psplot=psplot,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list

if power eq 'high' then goto,high_power

t0_ar1 = 3.0
t0_ar2 = 3.0 
t0_ar3 = 3.0 
t0_ar4 = 3.5
t0_ar5 = 2.5
t0_ar6 = 3.0


read_signal_mrm,0L,35157,'UVS','CFA03A',tline1,ar1,2
read_signal_mrm,0L,35158,'UVS','CFA03A',tline2,ar2,2
read_signal_mrm,0L,35158,'UVS','CFA03A',tline3,ar3,2
ar3 = ar3/1.1e21 * 0.82e21
read_signal_mrm,0L,35358,'UVS','CFA03A',tline4,ar4,2
read_signal_mrm,0L,35381,'UVS','CFA03A',tline5,ar5,2
read_signal_mrm,0L,35165,'UVS','CFA03A',tline6,ar6,2

;window,0,xs=1500,ys=1200
;!p.multi=[0,2,3]
;!p.thick=2.0
;!p.charsize=2.0

plot ,tline1-t0_ar1,smooth(ar1,300)/1e21,xr=[-0.5,3.],xs=1,col=0,back=255,yr=[0,3],$
      /nodata,xtitle='Time - t!lAr!n [s]',ytitle='[10!u21!n #/s]'
oplot,tline2-t0_ar2,smooth(ar2,300)/1e21,col=colors.red
oplot,tline3-t0_ar3,smooth(ar3,300)/1e21,col=colors.blue
oplot,tline4-t0_ar4,ar4/1e21,col=colors.aqua
oplot,tline5-t0_ar5,ar5/1e21,col=colors.gold
id = where(tline6 le t0_ar6+1.0)
;oplot,tline6[id]-t0_ar6,(smooth(ar6,300))[id]/1e21,col=colors.green
plot_legend,'(a) Ar seeding',col=colors.black

read_signal_mrm,0L,35157,'MAC','Tdiv',t1,tdiv1,2
read_signal_mrm,0L,35158,'MAC','Tdiv',t2,tdiv2,2
read_signal_mrm,0L,35167,'MAC','Tdiv',t3,tdiv3,2
read_signal_mrm,0L,35358,'MAC','Tdiv',t4,tdiv4,2
read_signal_mrm,0L,35381,'MAC','Tdiv',t5,tdiv5,2
read_signal_mrm,0L,35165,'MAC','Tdiv',t6,tdiv6,2

plot, t1-t0_ar1,smooth(tdiv1,40),xr=[-0.5,3.],xs=1,yr=[-5,15],col=0,back=255,$
      /nodata,xtitle='Time - t!lAr!n [s]',ytitle='[eV]'
oplot,t2-t0_ar2,smooth(tdiv2,40),col=colors.red
oplot,t3-t0_ar3,smooth(tdiv3,40),col=colors.blue
oplot,t4-t0_ar4,smooth(tdiv4,40),col=colors.aqua
oplot,t5-t0_ar5,smooth(tdiv5,40),col=colors.gold
id = where(t6 le t0_ar6+1.0)
;oplot,t6[id]-t0_ar6,(smooth(tdiv6,40))[id],col=colors.green
plot_legend,'(d) Tdiv',col=colors.black


diag   = 'BPT'
trace1 = 'Pr_main'
trace2 = 'Pr_tot'
exp    = 'davidp' 

read_signal_mrm,0L,35157,diag,trace1,tline1,prad1,2,exp=exp
read_signal_mrm,0L,35157,diag,trace2,tline1,pradtot1,2,exp=exp
read_signal_mrm,0L,35158,diag,trace1,tline2,prad2,2,exp=exp
read_signal_mrm,0L,35158,diag,trace2,tline2,pradtot2,2,exp=exp
read_signal_mrm,0L,35167,diag,trace1,tline3,prad3,2,exp=exp
read_signal_mrm,0L,35167,diag,trace2,tline3,pradtot3,2,exp=exp
read_signal_mrm,0L,35358,diag,trace1,tline4,prad4,2,exp=exp
read_signal_mrm,0L,35358,diag,trace2,tline4,pradtot4,2,exp=exp
read_signal_mrm,0L,35381,diag,trace1,tline5,prad5,2,exp=exp
read_signal_mrm,0L,35381,diag,trace2,tline5,pradtot5,2,exp=exp
read_signal_mrm,0L,35165,diag,trace1,tline6,prad6,2,exp=exp
read_signal_mrm,0L,35165,diag,trace2,tline6,pradtot6,2,exp=exp

plot ,tline1-t0_ar1,smooth(prad1,3)/1e6,xr=[-0.5,3.],xs=1,yr=[0,15],col=0,back=255,$
      /nodata,xtitle='Time - t!lAr!n [s]',ytitle='[MW]'
oplot,tline2-t0_ar2,smooth(prad2,3)/1e6,col=colors.red
oplot,tline3-t0_ar3,smooth(prad3,8)/1e6,col=colors.blue
oplot,tline4-t0_ar4,smooth(prad4,3)/1e6,col=colors.aqua
oplot,tline5-t0_ar5,smooth(prad5,3)/1e6,col=colors.gold
id = where(tline6 le t0_ar6+1.0)
;oplot,tline6[id]-t0_ar6,(smooth(prad6,3))[id]/1e6,col=colors.green
plot_legend,'(b) Pcore',col=colors.black

read_signal_mrm,0L,35157,'UVS','N_tot',tline_n1,rates_n1,2
read_signal_mrm,0L,35158,'UVS','N_tot',tline_n2,rates_n2,2
read_signal_mrm,0L,35167,'UVS','N_tot',tline_n3,rates_n3,2
read_signal_mrm,0L,35358,'UVS','N_tot',tline_n4,rates_n4,2
read_signal_mrm,0L,35381,'UVS','N_tot',tline_n5,rates_n5,2
read_signal_mrm,0L,35165,'UVS','N_tot',tline_n6,rates_n6,2
read_signal_mrm,0L,35157,'UVS','D_tot',tline_d1,rates_d1,2
read_signal_mrm,0L,35158,'UVS','D_tot',tline_d2,rates_d2,2
read_signal_mrm,0L,35167,'UVS','D_tot',tline_d3,rates_d3,2
read_signal_mrm,0L,35358,'UVS','D_tot',tline_d4,rates_d4,2
read_signal_mrm,0L,35381,'UVS','D_tot',tline_d5,rates_d5,2
read_signal_mrm,0L,35165,'UVS','D_tot',tline_d6,rates_d6,2

cn_flux1 = (rates_n1/7) / (rates_d1 + rates_n1/7)
cn_flux2 = (rates_n2/7) / (rates_d2 + rates_n2/7)
cn_flux3 = (rates_n3/7) / (rates_d3 + rates_n3/7)
cn_flux4 = (rates_n4/7) / (rates_d4 + rates_n4/7)
cn_flux5 = (rates_n5/7) / (rates_d5 + rates_n5/7)
cn_flux6 = (rates_n6/7) / (rates_d6 + rates_n6/7)

plot ,tline_n1-t0_ar1,smooth(cn_flux1*100,40),xr=[-0.5,3.],xs=1,yr=[0,15],col=0,back=255,$
      /nodata,xtitle='Time - t!lAr!n [s]',ytitle='[%]'
oplot,tline_n2-t0_ar2,smooth(cn_flux2*100,40),col=colors.red
oplot,tline_n3-t0_ar3,smooth(cn_flux3*100,40),col=colors.blue
oplot,tline_n4-t0_ar4,smooth(cn_flux4*100,40),col=colors.aqua
oplot,tline_n5-t0_ar5,smooth(cn_flux5*100,40),col=colors.gold
id = where(tline_n6 le t0_ar6+1.0)
;oplot,tline_n6[id]-t0_ar6,(smooth(cn_flux6*100,40))[id],col=colors.green
plot_legend,'(e) c!lN!n [Flux]',col=colors.black


plot,tline1-t0_ar1,smooth(pradtot1,3)/1e6,xr=[-0.5,3.],xs=1,yr=[0,15],col=0,back=255,$
      /nodata,xtitle='Time - t!lAr!n [s]',ytitle='[MW]'
oplot,tline2-t0_ar2,smooth(pradtot2,3)/1e6,col=colors.red 
oplot,tline3-t0_ar3,smooth(pradtot3,8)/1e6,col=colors.blue
oplot,tline4-t0_ar4,smooth(pradtot4,3)/1e6,col=colors.aqua
oplot,tline5-t0_ar5,smooth(pradtot5,3)/1e6,col=colors.gold
id = where(tline6 le t0_ar6+1.0)
;oplot,tline6[id]-t0_ar6,(smooth(pradtot6,3))[id]/1e6,col=colors.green
plot_legend,'(c) Ptot',col=colors.black


runcn,shots=35157,/load,time=cn_time1,cn_mean=cn_mean1
runcn,shots=35158,/load,time=cn_time2,cn_mean=cn_mean2
runcn,shots=35167,/load,time=cn_time3,cn_mean=cn_mean3
runcn,shots=35358,/load,time=cn_time4,cn_mean=cn_mean4
runcn,shots=35381,/load,time=cn_time5,cn_mean=cn_mean5
runcn,shots=35165,/load,time=cn_time6,cn_mean=cn_mean6


plot ,cn_time1-t0_ar1,smooth(cn_mean1*100,5),xr=[-0.5,3.],xs=1,yr=[0,15],col=0,back=255,$
      /nodata,xtitle='Time - t!lAr!n [s]',ytitle='[%]'
oplot,cn_time2-t0_ar2,smooth(cn_mean2*100,5),col=colors.red
oplot,cn_time3-t0_ar3,smooth(cn_mean3*100,5),col=colors.blue
oplot,cn_time4-t0_ar4,smooth(cn_mean4*100,5),col=colors.aqua
oplot,cn_time5-t0_ar5,smooth(cn_mean5*100,5),col=colors.gold
id = where(cn_time6 le t0_ar6+1.0)
;oplot,cn_time6[id]-t0_ar6,(smooth(cn_mean6*100,5))[id],col=colors.green
plot_legend,'(f) c!lN!n [Spec]',col=colors.black

goto, finish

high_power:

t0_ar1 = 2.6
t0_ar2 = 2.6


read_signal_mrm,0L,33257,'UVS','CFA03A',tline1,ar1,2
read_signal_mrm,0L,33258,'UVS','CFA03A',tline2,ar2,2


plot ,tline1-t0_ar1,smooth(ar1,1)/1e21,xr=[-0.5,3.],xs=1,col=0,back=255,yr=[0,3],$
      xtitle='Time - t!lAr!n [s]',ytitle='[10!u21!n #/s]'
oplot,tline2-t0_ar2,smooth(ar2,1)/1e21,col=colors.red
plot_legend,'(a) Ar seeding',col=colors.black

read_signal_mrm,0L,33257,'MAC','Tdiv',t1,tdiv1,2
read_signal_mrm,0L,33258,'MAC','Tdiv',t2,tdiv2,2

plot, t1-t0_ar1,smooth(tdiv1,40),xr=[-0.5,3.],xs=1,yr=[-5,15],col=0,back=255,$
      xtitle='Time - t!lAr!n [s]',ytitle='[eV]'
oplot,t2-t0_ar2,smooth(tdiv2,40),col=colors.red

plot_legend,'(d) Tdiv',col=colors.black


diag   = 'BPT'
trace1 = 'Pr_main'
trace2 = 'Pr_tot'
exp    = 'davidp' 

read_signal_mrm,0L,33257,diag,trace1,tline1,prad1,2,exp=exp
read_signal_mrm,0L,33257,diag,trace2,tline1,pradtot1,2,exp=exp
read_signal_mrm,0L,33258,diag,trace1,tline2,prad2,2,exp=exp
read_signal_mrm,0L,33258,diag,trace2,tline2,pradtot2,2,exp=exp

plot ,tline1-t0_ar1,smooth(prad1,3)/1e6,xr=[-0.5,3.],xs=1,yr=[0,15],col=0,back=255,$
      xtitle='Time - t!lAr!n [s]',ytitle='[MW]'
oplot,tline2-t0_ar2,smooth(prad2,3)/1e6,col=colors.red
plot_legend,'(b) Pcore',col=colors.black

read_signal_mrm,0L,33257,'UVS','N_tot',tline_n1,rates_n1,2
read_signal_mrm,0L,33258,'UVS','N_tot',tline_n2,rates_n2,2
read_signal_mrm,0L,33257,'UVS','D_tot',tline_d1,rates_d1,2
read_signal_mrm,0L,33258,'UVS','D_tot',tline_d2,rates_d2,2

cn_flux1 = (rates_n1/7) / (rates_d1 + rates_n1/7)
cn_flux2 = (rates_n2/7) / (rates_d2 + rates_n2/7)

plot ,tline_n1-t0_ar1,smooth(cn_flux1*100,40),xr=[-0.5,3.],xs=1,yr=[0,15],col=0,back=255,$
      xtitle='Time - t!lAr!n [s]',ytitle='[%]'
oplot,tline_n2-t0_ar2,smooth(cn_flux2*100,40),col=colors.red
plot_legend,'(e) c!lN!n [Flux]',col=colors.black


plot,tline1-t0_ar1,smooth(pradtot1,3)/1e6,xr=[-0.5,3.],xs=1,yr=[0,25],col=0,back=255,$
      xtitle='Time - t!lAr!n [s]',ytitle='[MW]'
oplot,tline2-t0_ar2,smooth(pradtot2,3)/1e6,col=colors.red 
plot_legend,'(c) Ptot',col=colors.black



finish:
stop

End

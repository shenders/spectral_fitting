Pro iaea_30554,debug=debug

shot=30554
setgraphics,nrow=3,ncol=1,colors=colors
xr=[3.0,5.5]
x=-1 & y=-1 & read_signal_mrm,0L,shot,'MAC','Tdiv',x,y,2,exp=exp

plot,x,y,xr=xr,back=colors.white,col=colors.black

x=-1 & y=-1 & read_signal_mrm,0L,shot,'UVS','N_tot',x,y,2,exp=exp

plot,x,y,xr=xr,back=colors.white,col=colors.black

x=-1 & y=-1 & read_signal_mrm,0L,shot,'UVS','D_tot',x,y,2,exp=exp

oplot,x,y,col=colors.red

x=-1 & y=-1 & read_signal_mrm,0L,shot,'LSD','te-ua4',x,y,2,exp=exp
id = where(x ge xr[0] and x le xr[1]) & x=x[id] & y=y[id]
telm = find_elms(shot,x) 
idelm = where(telm ge 3.0)
xx=x[idelm]
yy=y[idelm]

plot,x,y,xr=xr,back=colors.white,col=colors.black,yr=[0,50]
oplot,xx,yy,col=colors.red

stop

End

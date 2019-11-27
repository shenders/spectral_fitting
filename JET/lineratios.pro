Pro lineratios,psplot=psplot

    npts = 20
    temp = [4.0,6.0,8.0,10]
    dens = adas_vector(high=1e15,low=1e13,num=npts)
    tec1 = fltarr(npts,npts)
    tec2 = fltarr(npts,npts)

    for i=0,n_elements(temp)-1 do begin
	
	te = fltarr(npts)+temp[i]
	
	read_adf15,file='/home/adas/adas/adf15/pec96#ne/pec96#ne_pju#ne1.dat',block=17,te=te,dens=dens,data=data1_ex
	read_adf15,file='/home/adas/adas/adf15/pec96#ne/pec96#ne_pju#ne1.dat',block=70,te=te,dens=dens,data=data1_rc

	read_adf15,file='/home/adas/adas/adf15/pec96#ne/pec96#ne_pju#ne1.dat',block=24,te=te,dens=dens,data=data2_ex
	read_adf15,file='/home/adas/adas/adf15/pec96#ne/pec96#ne_pju#ne1.dat',block=77,te=te,dens=dens,data=data2_rc
    	
	run_adas405,year='96',uid='adas',elem='ne',te=te,dens=dens,frac=frac
	
	tec1[*,i] = data1_ex * frac.ion[*,1] + data1_rc * frac.ion[*,2]
	tec2[*,i] = data2_ex * frac.ion[*,1] + data2_rc * frac.ion[*,2]

    endfor
    
    setgraphics,nrow=1,ncol=2,xs=800,ys=1000,psplot=psplot,filename=filename,colors=colors,/portrait

    plot,dens*1e6,tec1[*,0]/tec2[*,0],/xlog,back=255,col=0,ytitle='!u4!nP-!u4!nP / !u2!nD-!u2!nP',xtitle='n!le!n [m!u-3!n]'
    oplot,dens*1e6,tec1[*,1]/tec2[*,1],col=colors.blue
    oplot,dens*1e6,tec1[*,2]/tec2[*,2],col=colors.green
    oplot,dens*1e6,tec1[*,3]/tec2[*,3],col=colors.red
    oplot,dens*1e6,tec1[*,4]/tec2[*,4],col=colors.orange
    oplot,dens*1e6,tec1[*,5]/tec2[*,5],col=colors.beige
    
    te = adas_vector(high=15,low=1.0,num=npts)
    run_adas405,year='96',uid='adas',elem='ne',te=te,dens=1e13+fltarr(npts),frac=frac
    read_adf15,file='/home/adas/adas/adf15/pec96#ne/pec96#ne_pju#ne1.dat',block=17,te=te,dens=1e14+fltarr(npts),data=data1_ex
    read_adf15,file='/home/adas/adas/adf15/pec96#ne/pec96#ne_pju#ne1.dat',block=70,te=te,dens=1e14+fltarr(npts),data=data1_rc
    tec1  = data1_ex * frac.ion[*,1] + data1_rc * frac.ion[*,2]
    plot,te,tec1/max(tec1),col=0,back=255,ytitle='Total emission coefficient',xtitle='T!le!n [eV]'
    

End

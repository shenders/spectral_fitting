Pro scaling_mr,psplot=psplot,tt=tt


	Psep = [5.80 , 4.05 , 4.31 , 6.97 , 5.77 , 5.67 , 5.91 , 7.67 , 8.99 , 14.21]
	nsep = [1.65 , 3.04 , 2.57 , 1.78 , 1.80 , 2.04 , 2.80 , 2.76 , 1.77 , 2.41 ] 
	Ip   = [1.04 , 0.84 , 0.84 , 0.84 , 0.84 , 0.84 , 0.84 , 1.00 , 1.00 , 1.20 ]
	amin = [0.49 , 0.50 , 0.50 , 0.51 , 0.51 , 0.50 , 0.49 , 0.49 , 0.50 , 0.50 ]
	kappa= [1.76 , 1.68 , 1.69 , 1.72 , 1.73 , 1.74 , 1.73 , 1.70 , 1.69 , 1.75 ]
	Bt   = [2.50 , 2.50 , 2.50 , 2.50 , 2.50 , 2.50 , 2.50 , 2.50 , 2.50 , 2.50 ]
	spec = [12.3 , 2.14 , 3.00 , 8.91 , 9.32 , 6.53 , 2.65 , 4.38 , 14.6 , 7.96 ]
	err  = [2.68 , 0.76 , 1.04 , 1.69 , 2.17 , 1.42 , 0.82 , 1.46 , 2.89 , 1.09 ]
	cz   = fltarr(n_elements(Psep))
	for i=0,n_elements(psep)-1 do begin	
		x=lengfunc(scen=6,psep=psep[i],$
		                  nsep=nsep[i]*1e19,$
				  ip=ip[i]    ,$
				  amin=amin[i],$
				  kappa=kappa[i],$
				  tt=tt,$
				  bt=bt[i],/arne)
				
		cZ[i] = x
	endfor
	
	user_psym,1,/fill
	setgraphics,xs=xs,ys=ys,ncol=1,nrow=1,psplot=psplot,filename='figures/scaling_mreinke.ps',colors=colors,colpick=colpick,collabel=collabel,full_list=full_list	
	
	plot,cz*100/10,spec,psym=8,xr=[0,20],yr=[0,20],col=colors.black,back=colors.white,$
	ytitle='Spectroscopy c!lN!n [%]',xtitle='MReinke2017 c!lN!n [%]'
	err_plot,cz*100/10,spec,err,col=colors.black
	oplot,[0,100],[0,100],linest=5,col=colors.black
stop
End

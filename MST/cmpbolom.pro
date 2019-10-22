Pro cmpbolom,shot

	runabel,shots=shot,/noplot,/load,time=time,prad=pradabel

	read_signal_mrm,0L,shot,'BPT','Pr_tot',t1,pr_tot,2,exp='davidp'
	read_signal_mrm,0L,shot,'BPT','Pr_main',t2,pr_main,2,exp='davidp'
	read_signal_mrm,0L,shot,'BPT','Pr_div',t3,pr_div,2,exp='davidp'

	read_signal_mrm,0L,shot,'BPD','Pradtot',tline,pradtot,2
	read_signal_mrm,0L,shot,'BPD','Prad',tline,prad,2

	setgraphics,colors=colors	
	plot,t1,pr_tot,/nodata,xr=[2.0,6.0],yr=[0,1e7],back=colors.white,col=colors.black,xtitle='Time [s]',ytitle='[MW]'
	
	oplot,t1,pr_tot,col=colors.black
	oplot,t2,pr_main,col=colors.red
		
	oplot,tline,pradtot,col=colors.black,linest=5
	oplot,tline,prad,col=colors.red,linest=5
	
	user_psym,1,/fill
	oplot,time,pradabel*1e6,col=colors.red,psym=-8
stop
end

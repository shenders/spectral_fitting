Pro jet_radiated_power,shot,time,pr_tot=pr_tot,pr_main=pr_main,pr_div=pr_div,sm=sm

    	ppfread,shot=shot,dda='BOLP',dtype='TOBP',data=pr_main,t=tline 
    	ppfread,shot=shot,dda='BOLP',dtype='TODP',data=pr_div,t=tline 
    	ppfread,shot=shot,dda='BOLP',dtype='TOPO',data=pr_tot,t=tline 
	
	if keyword_set(sm)then begin
	    pr_main = smooth(pr_main,sm)
	    pr_div  = smooth(pr_div,sm)
	    pr_tot  = smooth(pr_tot,sm)
	endif
	
	if keyword_set(time)then begin
	    pr_main = interpol(pr_main,tline,time)
	    pr_div  = interpol(pr_div,tline,time)
	    pr_tot  = interpol(pr_tot,tline,time)
    	endif else time = tline

End

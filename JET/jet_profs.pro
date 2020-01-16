Pro jet_profs,shot,$
    	      tesep = tesep,$
	      nesep=nesep,$
	      err=nesep_err,$
	      debug=debug
	      
    if shot eq 96075 then seq = 182
    if shot eq 85417 then seq = 842
    if shot eq 85264 then seq = 391
    if shot eq 85265 then seq = 419
    if shot eq 85266 then seq = 902
    if shot eq 85267 then seq = 1061
    if shot eq 85268 then seq = 413

    if shot eq 85270 then seq = 467
    if shot eq 85272 then seq = 482
    if shot eq 85274 then seq = 1128
    if shot eq 85276 then seq = 445
    if shot eq 85423 then seq = 455
    if shot eq 96227 then seq = 261

        lims = [0.99,1.01]
    
    ppfread,shot=shot,dda='PED',dtype='TEFP',ppfuid='shenders',seq=seq,data=data,x=x    
    te_vals = interpol(data,x,findgen(10)*(lims[1]-lims[0])/9.0+lims[0])
    stats = moment(te_vals)
    tesep = stats[0]
    tesep_err = sqrt(stats[1])
    if keyword_set(debug)then begin
    	setgraphics,colors=colors 	
    	plot,x,data,xr=[0.8,1.2],/nodata,back=colors.white,col=colors.black,ys=9
	oplot,x,data,col=colors.black
	oplot,[0,2],[100,100],col=colors.black,linest=5
    endif
    ppfread,shot=shot,dda='PED',dtype='NEFP',ppfuid='shenders',seq=seq,data=data,x=x    
    ne_vals = interpol(data/1e19,x,findgen(10)*(lims[1]-lims[0])/9.0+lims[0])
    stats = moment(ne_vals)
    nesep = stats[0]
    nesep_err = 1.0/sqrt(n_elements(ne_vals)) * sqrt(stats[1])
    
    if keyword_set(debug)then begin 	
    	print,nesep
    	print,tesep
    	axis,yaxis=1,/save,yr=[0,max(data)],col=colors.blue
	oplot,x,data,col=colors.blue
    stop
    endif
end












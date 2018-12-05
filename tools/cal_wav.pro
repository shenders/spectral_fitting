FUNCTION cal_wav,x,y,shot,diag,spec


	    nist_vac = [373.791,375.574,376.094,379.234,388.909,397.01,399.613,404.245,410.174]
	    nist_air = nist_vac-0.113
	    id = where(nist_air ge min(x) and nist_air le max(x))
	    nist_air=nist_air[id] 
	    
	    plot,x,y,/ylog,xr=[min(x),max(x)]
	    xmax=-1
	    print,'Click on centre wavelength for line: '
	    for i=0,n_elements(nist_air)-1 do begin
	    	print,nist_air(i)
		oplot,[nist_air(i),nist_air(i)],[1e10,1e22],linest=5
		cursor,x1,y1,/up
		idx  = where(x ge x1-0.1 and x le x1+0.1)
		idm  = where(y[idx] eq max(y[idx]))
		xmax = [xmax,x[idx[idm[0]]]]
	    endfor
	    xmax=xmax[1:*]
	    wcal = interpol(nist_air - xmax,nist_air,x)
	    print,nist_air - xmax
	    print,nist_air
	    print,xmax
	    plot,x,wcal
	    wcal_file='Save/wcal'+DIAG+STRING(shot,format='(I5)')+'.sav'
	    save,file=wcal_file,wcal
	    STOP
return,wcal
END	    


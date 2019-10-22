Pro straight,x , a, f, pder

	f    = a[0] * x 
	pder = [ [ x ]  ]	

end

Pro regress,xdata,ydata,err,xfit,yfit,chisq,a,func_type=func_type,debug=debug

	if ~keyword_set(func_type)then func_type='straight'

	id = sort(xdata)
	xx = xdata(id)
	yy = ydata(id)
	ee = err(id)
	a  = [1.0]
	x=curvefit(xx,yy,ee,a,sig,function_name='straight',itmax=100,chisq=chisq)
	xfit = findgen(100)*(max(xx)-min(xx))/99.0 + min(xx)
	yfit = a[0] * xfit 
	
	if keyword_set(debug)then begin
		plot,xx,yy,psym=4,yr=[0,max(yy)*1.2],xr=[0,max(xx)]
		oplot,xfit,yfit,linest=5
		errors,xx,yy,ystd=ee
		print,chisq
		stop
	endif
	
End




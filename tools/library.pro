Function constants,symbol,help=help,split=split

    if keyword_set(help)then begin
    	print,'Constants in database'
	print,'---------------------'
	print,'c   : speed of light'
	print,'e   : electron charge'
	print,'me  : mass of electron'
	print,'mi  : mass of proton'
	print,"h   : planck's constant"
	print,'hc  : h * c'
	print,"kb  : boltzmann's constant"
	print,'ry  : rydberg energy'
	print,'hbar: h/pi'
	print,'pi  : pi'
	print,'a0  : bohr radius'
	print,'mu0 : vacuum permeability'
	print,'epsilon0: vacuum permittivity'
	print,'alpha: fine structure constant'
	print,'ev_to_k: convert eV to kelvin'
	print,'recip_4_pi: 1/(4pi)'
	return,-1
    endif
    if ~keyword_set(split)then begin
        if symbol eq 'pi' then return,3.14159265359
    	if symbol eq 'c' then return,2.9979E+8
    	if symbol eq 'e' then return,1.60217662E-19
    	if symbol eq 'me' then return,9.10938356E-31
    	if symbol eq 'h' then return,6.62607004E-34
    	if symbol eq 'hc' then return,1.98644582E-25
    	if symbol eq 'kb' then return,1.38064852E-23
    	if symbol eq 'hbar' then return,1.054571817E-34
    	if symbol eq 'mp' then return,1.6726219E-27
    	if symbol eq 'ry' then return,13.605693
    	if symbol eq 'ev_to_k' then return,1.1604518E4
    	if symbol eq 'alpha' then return,7.2974E-03
    	if symbol eq 'mu0' then return,1.256637062E-6
    	if symbol eq 'epsilon0' then return,8.8541878128E-12
    	if symbol eq 'a0' then return,5.2917721067E-11
    	if symbol eq 'recip_4_pi'then return,7.957747E-2
    endif else begin
        if symbol eq 'pi' then return,[3.14159265359,0]
    	if symbol eq 'c' then return,[2.9979,8]
    	if symbol eq 'e' then return,[1.60217662,-19]
    	if symbol eq 'me' then return,[9.10938356,-31]
    	if symbol eq 'h' then return,[6.62607004,-34]
    	if symbol eq 'hc' then return,[1.98644582,-25]
    	if symbol eq 'kb' then return,[1.38064852,-23]
    	if symbol eq 'hbar' then return,[1.054571817,-34]
    	if symbol eq 'mp' then return,[1.6726219,-27]
    	if symbol eq 'ry' then return,[1.3605693,1]
    	if symbol eq 'ev_to_k' then return,[1.1604518,4]
    	if symbol eq 'alpha' then return,[7.2974,-03]
    	if symbol eq 'mu0' then return,[1.256637062,-6]
    	if symbol eq 'epsilon0' then return,[8.8541878128,-12]
    	if symbol eq 'a0' then return,[5.2917721067,-11]
    	if symbol eq 'recip_4_pi'then return,[7.957747,-2]
	    
    end
    print,'Symbol not recognised'
    return,-1
end
    


Pro load_jet_paths
    !PATH=!PATH + ':' + $
    expand_path( '+~cxs/utilities' ) + ':' + $
    expand_path( '+~cxs/calibration' ) + ':' + $
    expand_path( '+~cxs/instrument_data' )
End
Function graphpos,xidx,yidx,rows,cols,xspc=xspc,yspc=yspc,xlim=xlim,ylim=ylim

    if ~keyword_set(xlim)then xlim = [0.1,0.9]
    if ~keyword_set(ylim)then ylim = [0.15,0.9]
    if ~keyword_set(xspc)then xspc = 0.00
    if ~keyword_set(yspc)then yspc = 0.00
    xdif = xlim[1]-xlim[0]-xspc*(cols-1)-xspc
    ydif = ylim[1]-ylim[0]-yspc*(rows-1)-yspc
    tgrp = rows * cols
    xdiv = xdif / cols
    ydiv = ydif / rows
    x0   = xlim[0] + xdiv * xidx + xspc * xidx + xspc/2.0
    x1   = x0 + xdiv
    y1   = 1.0 - (ylim[0] + ydiv * yidx + yspc * yidx)+yspc/2.0
    y0   = y1 - ydiv
return,[x0,y0,x1,y1]
end
Pro legend,text,col,yshift=yshift,xshift=xshift,ylog=ylog,xlog=xlog
	if ~keyword_set(yshift)then yshift=0.0
	if ~keyword_set(xshift)then xshift=0.0
	xr = !x.crange
	yr = !y.crange	
	y0 = yr[0]+(yr[1] - yr[0])*(0.8+yshift)
	x0 = xr[0]+(xr[1] - xr[0])*(0.7+xshift)
	if keyword_set(xlog)then begin
    	    x0 = 10^(x0)
	endif 	
	if keyword_set(ylog)then begin
	    y0 = 10^(y0)
	endif 

    	xyouts,x0,y0,text,col=col,charsize=1.4

End
Function calc_cn,shot,los,transmission,te,dens,data,sm,err=err_1,jet=jet,horizontal=horizontal

	cn = data.nii3995
	for i=0,n_elements(cn[0,*])-1 do begin
    	    id_nan = where(finite(dens[*,i]) eq 0)
	    if id_nan[0] eq -1 then begin
	    	atomdb,te[*,i],dens[*,i],tec3995=tec3995
	    	if keyword_set(jet)then begin
    		    dl = jet_length(data.tdiv,data.los_names[i],upperdl=upperdl,lowerdl=lowerdl,doplot=doplot,psplot=psplot,horizontal=horizontal)
	    	endif else begin
	    	    dl = length(data.tdiv,shot,los,upperdl=upperdl,lowerdl=lowerdl)
		end
		cn[*,i] = smooth(data.nii3995[*,i],sm,/edge_truncate) * 4.0 * !pi * transmission /( tec3995 * dens[*,i] ) / (dens[*,i] * 1e6) / dl
	    endif else cn[*,i]=-1
	endfor
	return,cn
	
End

Function r2,x,y,b=b,m=m	
	xsum = total(x)
	ysum = total(y)
	xysum=0.0 & for i=0,n_elements(x)-1 do xysum = xysum + x[i]*y[i]
	xsqr=0.0 & for i=0,n_elements(x)-1 do xsqr = xsqr + x[i]^2
	ysqr=0.0 & for i=0,n_elements(x)-1 do ysqr = ysqr + y[i]^2
	N  = n_elements(x)
	m  = (N*xysum -xsum*ysum)/(N*xsqr-xsum*xsum)
	b  = (xsqr*ysum - xsum*xysum)/(N*xsqr-xsum*xsum)
	r2 = (N*xysum - xsum*ysum)^2/( (N*xsqr - xsum*xsum)*(N*ysqr-ysum*ysum))
	return,r2
End

Pro bindata,time,yvals,trange,avr,err

	x   = moment(yvals[where(time ge trange[0] and time le trange[1])])
	avr = x[0]
	err = sqrt(x[1])

End

PRO errors,x,y,xstd=xstd,ystd=ystd,col=col,xupper=xupper,xlower=xlower

    	ymax=!y.crange(1)
    	ymin=!y.crange(0)

    	xmax=!x.crange(1)
    	xmin=!x.crange(0)

	if keyword_set(xupper)then begin
		for i=0,n_elements(x)-1 do oplot,[xlower[i],xupper[i]],[y[i],y[i]],col=col
		y2 = (ymax - ymin) * 0.01
		for i=0,n_elements(x)-1 do oplot,[xupper[i],xupper[i]],[y[i]-y2,y[i]+y2],col=col
		for i=0,n_elements(x)-1 do oplot,[xlower[i],xlower[i]],[y[i]-y2,y[i]+y2],col=col
	endif

	if keyword_set(xstd)then begin
		for i=0,n_elements(x)-1 do oplot,[x[i]-xstd[i],x[i]+xstd[i]],[y[i],y[i]],col=col
		y2 = (ymax - ymin) * 0.01
		for i=0,n_elements(x)-1 do oplot,[x[i]+xstd[i],x[i]+xstd[i]],[y[i]-y2,y[i]+y2],col=col
		for i=0,n_elements(x)-1 do oplot,[x[i]-xstd[i],x[i]-xstd[i]],[y[i]-y2,y[i]+y2],col=col
	endif
	if keyword_set(ystd)then begin
		for i=0,n_elements(x)-1 do oplot,[x[i],x[i]],[y[i]+ystd[i],y[i]-ystd[i]],col=col
		x2 = (xmax - xmin) * 0.01
		for i=0,n_elements(x)-1 do oplot,[x[i]-x2,x[i]+x2],[y[i]+ystd[i],y[i]+ystd[i]],col=col
		for i=0,n_elements(x)-1 do oplot,[x[i]-x2,x[i]+x2],[y[i]-ystd[i],y[i]-ystd[i]],col=col
	endif
END	
Pro augped,shot,nesep,tesep,nesep_err,tesep_err,debug=debug,xmin=xmin,xmax=xmax

	if ~keyword_set(xmin)then xmin=0.999
	if ~keyword_set(xmax)then xmax=1.001
	dir = '/afs/ipp/home/s/shenders/augped/output/'
	shotstr = string(shot,format='(i5)')
	spawn,'ls '+dir,files
	id  = where(strpos(files,shotstr) ne -1)
	if id[0] eq -1 then stop,'No files exist with this shot number'
	if n_elements(id) gt 1 then begin
		for i=0,n_elements(id)-1 do print,i,': ',files[id[i]]
		read,val,prompt='More than 1 file exists, which file (e.g. 0,1...): '
		id = id[val] 
	endif else val=0
	file = files[id]
	print,'Reading file: ',file
	openr,unit,dir+file,/get_lun
	txt = ''
	itest = 0
	x = -1.0
	dens = -1.0
	te = -1.0
	dens_err = -1.0
	te_err = -1.0
	while ~eof(unit) do begin
		readf,unit,txt
		check='Fit to: Te'
		if strpos(txt,check) ne -1 then begin
			for i=0,2 do readf,unit,txt
			fin = 0
			while fin ne 1 do begin
				readf,unit,rmaj,rho,psi,fit,dfit
				x = [x,rho]
				te = [te,fit]
				te_err = [te_err,dfit]
				if rho eq 1.2 then fin = 1
			end
		endif		
		check='Fit to: ne'
		if strpos(txt,check) ne -1 then begin
			for i=0,2 do readf,unit,txt
			fin = 0
			while fin ne 1 do begin
				readf,unit,rmaj,rho,psi,fit,dfit
				dens = [dens,fit]
				dens_err = [dens_err,dfit]
				if rho eq 1.2 then fin = 1
			end
		endif		

	end
	close,unit
	free_lun,unit
	x = x[1:*]
	dens = dens[1:*]
	te = te[1:*]
	dens_err = dens_err[1:*]
	te_err = te_err[1:*]
	id = where(x ge xmin and x le xmax)
	print,string(mean(dens[id]),format='("ne,sep= ",E8.2)')
	if keyword_set(debug)then begin
		plot,x,dens>0
		oplot,[xmin,xmin],[-5e20,5e20],linest=5
		oplot,[xmax,xmax],[-5e20,5e20],linest=5
		oplot,[0,2],[mean(dens[id]),mean(dens[id])],linest=5
		oplot,x[id],dens[id],thick=4,psym=5
		stop 
	endif
	nesep = mean(dens[id])
	Tesep = mean(te[id])
	nesep_err = mean(dens_err[id])
	Tesep_err = mean(te_err[id])
END

Pro deltaL,tdiv,dl,upperdl=upperdl,lowerdl=lowerdl,machine=machine,los=los,rov014=rov014

	if ~keyword_Set(machine)then machine='AUG'
	if ~keyword_Set(los)then los = '-'
	if machine eq 'AUG'then begin
		teprof = findgen(100)/99.0 * 99.0 + 1.0
		if los eq 'ROV-14' or los eq 'ROV014' then begin
			set_te = [-2.0,5,40] 
			del_l_u = interpol([0.069,0.107,0.035],$
	 		                    set_te,teprof)
		endif else begin
			if los eq 'ZON-01' then begin
				set_te = [-2.0,5,40] 
				del_l_u = interpol([0.069,0.107,0.035]+0.05,$
	 		                            set_te,teprof)
			
			endif else begin
				if los eq 'ZON-05' then begin
					set_te = [-2.0,5,40] 
					del_l_u = interpol([0.069,0.107,0.035],$
	 		                                    set_te,teprof)
				
				endif else begin
					if keyword_set(rov014) then begin
						set_te = [-2.0,5,40] 
						del_l_u = interpol([0.069,0.107,0.035],$
	 		                    			   set_te,teprof)					
					endif else begin
						set_te = [1.0,7.0,100]
						del_l_u = interpol([0.069,0.107,0.035],$
	 		   		                  set_te,teprof)				
					end
				end
			end
		end
		
		del_l_l = del_l_u * 0.8
	endif
	
	upperdl = interpol(del_l_u,teprof,tdiv)
	lowerdl = interpol(del_l_l,teprof,tdiv)
	dl = upperdl
	for i=0,n_elements(tdiv)-1 do dl[i]=mean([upperdl[i],lowerdl[i]])
END
	
PRO setgraphics,xs=xs,ys=ys,ncol=ncol,nrow=nrow,psplot=psplot,landscape=landscape,portrait=portrait,close=close,filename=filename,$
                colors=colors,colpick=colpick,collabel=collabel,full_list=full_list,title=title

	adas_colors,colors=colors
	
	colpick  = [colors.blue, colors.green,colors.orange,colors.red,colors.magenta]
	collabel = ['blue', 'green','orange','red','magenta']
	
	if keyword_set(full_list)then begin
		colpick = [colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink,$		   
		   colors.navy,colors.blue,colors.sky,colors.cyan,colors.aqua,colors.yellow,colors.gold,$
	           colors.orange,colors.peachpuff,colors.red,colors.magenta,colors.pink]
		collabel = ['navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink',$		    
		    'navy','blue','sky','cyan','aqua','yellow','gold',$
	            'orange','peachpuff','red','magenta','pink']	
	endif
	if ~keyword_set(ncol)then ncol=0
	if ~keyword_set(nrow)then nrow=0
	if ~keyword_set(psplot) and ~keyword_set(close) then begin
		window,/free,xs=xs,ys=ys,title=title
		!p.multi=[0,ncol,nrow]
		!p.charsize=2.0
		!p.thick=2.0
	endif else begin
		if keyword_set(close) and keyword_set(psplot) then device,/close else begin
			if keyword_set(psplot)then begin
				if ~keyword_set(filename)then filename='output.ps'
				!p.multi=[0,nrow,ncol]
				!p.charsize=2.0
				!p.font=0
				!p.thick=6.0
				!x.thick=6.0
				!y.thick=6.0
				!z.thick=6.0
				set_plot,'ps'
				if xs gt 100 then xs=xs/100
				if ys gt 100 then ys=ys/100
				device,color=1,xsize=xs,ysize=ys,/inches,bits_per_pixel=64,file=filename,$
				font_size=11,landscape=landscape,portrait=portrait,/encapsulated
			endif
		end
	end	

END


PRO basic_fit,te_val=te_val,func=func,dens=dens,tec3995=tec3995,exc3995=exc3995,rec3995=rec3995,use402=use402
	dens=adas_vector(high=1e15,low=1e13,num=100)
	if ~keyword_set(te_val)then te_val=3.5
	te = te_val+fltarr(100)
; 399.5 line
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=exc3995,block=15
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=rec3995,block=65
; 404.2 line	
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=exc4042,block=21
	read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=rec4042,block=71
; Ionisation balance
    	run_adas405, uid='adas', year='96', elem='n', te=te, dens=dens, frac=frac
; Line ratio at fixed temperature	
	tec3995 = (frac.ion[*,1] * exc3995 + frac.ion[*,2] * rec3995) 
	tec4042 = (frac.ion[*,1] * exc4042 + frac.ion[*,2] * rec4042) 

	if keyword_set(use402)then begin
	; 402.6 line
		read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=exc4026,block=19
		read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=rec4026,block=69
		tec3995 = (frac.ion[*,1] * exc4026 + frac.ion[*,2] * rec4026) 
	endif
	func    = tec4042 / tec3995
	
END
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
	    wcal_file='tmp/wcal'+DIAG+STRING(shot,format='(I5)')+'.sav'
	    save,file=wcal_file,wcal
return,wcal
END	    
FUNCTION find_elm,shot,tline,emiss,psplot=psplot,plot_stats=plot_stats,binsize=binsize,experiment=experiment,trange=trange
	
	
	IF keyword_set(trange)then begin
		id=where(tline ge trange[0] and tline le trange[1])
		tline=tline[id]
		emiss=emiss[id]
	endif
	
	tr=[min(tline),max(tline)]
	res = tline[1]-tline[0]
	tline=tline+res/2.0
;	read_data,'ELM',shot,f_elm,time,experiment=experiment
;    	diagnostic = 'ELM'
;	signalname = 'f_ELM'
	file = '~/data/elms/'+string(shot,format='(i5)')
	if file_test(file)then begin
		restore,file,/verb
		elmtime=telm
	endif else begin	
		read_signal_mrm,0L,shot,'ELM','f_ELM',time,f_elm,1	
		idgtzero = WHERE(time GT 0)
		time     = time[idgtzero]
		f_elm    = f_elm[idgtzero]
		idelms           = WHERE(time GE tr[0] AND time LE tr[1])
		time             = time[idelms]
		elmtime          = time
	end
	telm             = FLTARR(n_elements(tline))
	time_to_next_elm = FLTARR(n_elements(tline))
	FOR i=0,N_ELEMENTS(tline)-1 DO BEGIN
		id2      = WHERE(elmtime LE tline(i))
		id2      = MAX(id2)
		IF id2 EQ -1 THEN telm(i)=-1 ELSE telm(i)  = tline(i) - elmtime[id2] 
		IF id2 EQ N_ELEMENTS(elmtime)-1 THEN time_to_next_elm(i)=-1 ELSE time_to_next_elm(i)  = elmtime[id2+1]-tline(i)
	ENDFOR
	telm = telm*1e3
	time_to_next_elm=time_to_next_elm*1e3

	IF KEYWORD_SET(plot_stats)THEN BEGIN

		bins = [0.0001,0.0005,0.0010,0.0015,0.0020,0.0025,0.0030,0.0035,0.0040,0.0045,0.005,$
       			0.0055,0.0060,0.0065,0.0070,0.0075,0.0080,0.0085,0.0090,0.0095,0.010]*1e3	
		bin       = 0.0
		data1     = 0.0
		telm_plot = telm
		norm      = MAX(emiss)
		nii       = emiss/norm
		FOR i=0,N_ELEMENTS(bins)-2 DO BEGIN
       			id   = where(telm_plot ge bins(i) and telm_plot lt bins(i+1))
       			IF id[0] ne -1 THEN BEGIN
	     			bin   = [bin,(bins(i)+bins(i+1))/2.0]    
	     			ndata = N_ELEMENTS(id)
	     			data1 = [data1,(MOMENT(nii(id)))[0]] 
       			ENDIF
		ENDFOR
		bin   = bin[1:*]
		data1 = data1[1:*]*norm
		!p.charsize=2.0
		user_psym,1
		IF KEYWORD_SET(psplot)THEN makeps,file=title+'elm_stats_hi.ps',xs=8,ys=5 ELSE BEGIN
			window,0,xs=700,ys=600 
		ENDELSE	

		IF KEYWORD_SET(psplot)THEN makeps,file=title+'elm_stats_nii.ps',xs=8,ys=5
		plot,telm_plot,nii*norm,psym=8,xr=[0,10],/ylog 
		oplot,bin,data1
		oplot,[res,res]*1e3,[0,1e23],linest=5
	ENDIF
	RETURN,telm
END
pro read_signal_mrm, ier,       $
                 shot,      $			;<-
                 diag_name, $			;<-
                 sig_name,  $			;<-
                 time,      $			;->
                 sig,       $			;-> Daten
                 phys_dim,  $ 			;-> Physikalische Einheit
                 edition = edition, $
                 time_ft = time_ft, $
                 indices = indices,  $
                 exp=exp,   $
		 text=text, $
		 debug=debug
;		 format_ext=format_ext


 if (n_params() le 6) then begin
	print,'usage:    read_signal, ier,shot,diag_name,sig_name,time,sig,phys_dim,edition=edition,time_ft=time_ft,indices=indices,exp=exp'
	return
 endif




libddww =  '/afs/ipp/aug/ads/lib64/@sys/libddww8.so'



 ;
 ; General information about data types:
 ;
 ; text description:
 format_txt = ['BYTE','CHAR','SHORT_INT','INTEGER','IEEE_FLOAT','IEEE_DOUBLE', $
               'LOGICAL','CHAR_REAL','U_SHORT','IBM_REAL','IBM_DOUBLE', $
               'CHAR_8','CHAR_16','CHAR_32','CHAR_48','CHAR_64','CHAR_72']
 ; DDWW array type at return:
 format_val = [1,2,3,4,5,6,7,8,9,10,11,1794,3842,7938,12034,16130,18178] 
 ; type which is used for reading in IDL and translated to DDWW:
 ;type_val   = [2,6,2,2,2,2,1,2,2, 2, 3,   6,   6,   6,    6,    6,    6]
 type_val   = [2,6,2,2,2,3,1,2,2, 2, 3,   6,   6,   6,    6,    6,    6]
 ; IDL array type:
 ;array_type = [4,7,4,4,4,4,3,4,4, 4, 5,   7,   7,   7,    7,    7,    7]       
 array_type = [4,7,4,4,4,5,3,4,4, 4, 5,   7,   7,   7,    7,    7,    7]       




 if (keyword_set(edition) eq 0) then ed = 0l else ed = long(edition)
  
 shot      = long(shot)
 ier       = 0L
 dia_ref   = 0L
 k1        = 0L
 k2        = 0L
 dim       = 0L
 ctrl      = 3L
 ncal      = 0L
 ntval     = 0L
 npretrig  = 0L
 tbeg      = 0.0 
 tend      = 0.0 
 if KEYWORD_SET(exp) then begin
    exp    = STRUPCASE(STRTRIM(exp, 2))
 endif else begin
    exp    = 'AUGD'
 endelse
 datum     = string('',format='(a18)')
 phys_dim  = string('',format='(a12)')

;_______________________________________________________________________________

;OPEN SHOTFILE

 s = call_external(libddww,'ddgetaug','ddopen',         $
                   ier,exp,diag_name,shot,ed,dia_ref,datum)
 if (ier gt 0) then begin
    s = call_external(libddww,'ddgetaug','xxerror',ier,3L,'ddopen:       ')
    return
 endif

    
;GET INFORMATION

 text	= '                                                                 '
 obuf	= lonarr(26)
 ier	= 0L
 s = CALL_EXTERNAL(libddww,'ddgetaug','ddobjhdr', $
                   ier,dia_ref,sig_name,obuf,text)
 if (ier ne 0) then begin
    s = CALL_EXTERNAL(libddww,'ddgetaug','xxerror',ier, ctrl, ' ')
    goto, close
 endif  
  
 obj_type = 0l
 s = CALL_EXTERNAL(libddww,'ddgetaug','ddobjval', $
                   ier,dia_ref,sig_name,'objtype',obj_type)
 if (ier ne 0) then begin
    s = CALL_EXTERNAL(libddww,'ddgetaug','xxerror',ier, ctrl, ' ')
    goto, close
 endif  
 
 obj_size = 0l
 s = CALL_EXTERNAL(libddww,'ddgetaug','ddobjval', $
                   ier,dia_ref,sig_name,'size',obj_size)
 if (ier ne 0) then $
    s = CALL_EXTERNAL(libddww,'ddgetaug','xxerror',ier, ctrl, ' ')
 k1 = 1l
 k2 = long(obj_size)

 dim_arr = lonarr(3)   
 s = CALL_EXTERNAL (libddww,'ddgetaug','ddobjval', $
                    ier,dia_ref,sig_name,'indices',dim_arr)
 if (ier ne 0) then $
    s = CALL_EXTERNAL(libddww,'ddgetaug','xxerror',ier, ctrl, ' ')

 relations = lonarr(8)
 s = CALL_EXTERNAL(libddww,'ddgetaug','ddobjval', $
                   ier,dia_ref,sig_name,'relations',relations)
 if (ier ne 0) then $
    s = CALL_EXTERNAL(libddww,'ddgetaug','xxerror',ier, ctrl, ' ')
 tmp = where(relations gt 0,count)

 if (count gt 0) then begin

   ;READ TIME-BASE

    ;
    ; read format of timebase:
    ;
    ier = 0L
    typ = 0L
    tname = '        '
    ind = lonarr(4)
    ; get timebase name for signal:
    s = CALL_EXTERNAL(libddww,'ddgetaug','ddsinfo', $
                   ier,dia_ref,sig_name,typ,tname,ind)
    if (ier ne 0) then s=CALL_EXTERNAL(libddww,'ddgetaug','xxerror',ier, ctrl, ' ')
    if(keyword_set(debug)) then print,'ddsinfo = ',s,sig_name,typ,tname,ind

    ier = 0L
    time_format = lonarr(3)
    ; get format for timebase itself:
    s = CALL_EXTERNAL(libddww,'ddgetaug','ddobjval', $
                   ier,dia_ref,tname,'format',time_format)
    if (ier ne 0) then s=CALL_EXTERNAL(libddww,'ddgetaug','xxerror',ier, ctrl, ' ')
    time_ind = type_val(where(time_format(0) eq format_val) > 0)
    time_type = (long(type_val(where(time_format(0) eq format_val) > 0)))(0)
    time_arr_typ = (array_type(where(time_format(0) eq format_val) > 0))(0)
    time_format_txt = (format_txt(where(time_format(0) eq format_val) > 0))(0)

    if(keyword_set(debug)) then begin
        print,'time_format = ',time_format
        print, 'time_ind = ',time_ind
        print, 'time_typ = ',time_type
        print, 'time_arr_typ = ',time_arr_typ
        print, 'time_format_txt = ',time_format_txt
    endif

    ier = 0L
    _tbeg = 0.0
    _tend = 0.0
    _ntval = 0L
    _npretrig = 0L
    s = CALL_EXTERNAL (libddww,'ddgetaug','ddtrange', $
                       ier,dia_ref,tname,_tbeg,_tend,_ntval,_npretrig)
    if (ier ne 0) then begin
	s = CALL_EXTERNAL (libddww,'ddgetaug','xxerror',ier,ctrl,' ')
	goto, close
    endif
    if(keyword_set(debug)) then begin
	print,'t=',_tbeg,_tend,'size=',_ntval,_npretrig
    endif

    if (keyword_set(debug)) then begin
      _k1 = 0L
      _k2 = 0L
      ier = 0L
      s = CALL_EXTERNAL (libddww,'ddgetaug','ddtindex', $
                       ier,dia_ref,tname,_tbeg,_tend,_k1,_k2)
;      ntval = long(k2-k1+1)
      print,'size = ',_ntval,_k1,_k2,obj_size,_tbeg,_tend
    endif
;    if (keyword_set(format_ext) ) then begin
;	case format_ext of
;	    'float': begin
;		time_type = 2L
;		time_si = long([1,1,4,1])
;		end
;	    'double': begin
;		time_type = 3L
;		time_si = long([1,1,5,1])
;		end
;	    else: begin
;		time_type = 2L
;		time_si = long([1,1,4,1])
;		end
;	endcase
;    endif else begin
;	time_type = 2L
;	time_si = long([1,1,4,1])
;    endelse

;;    max_ntval = 5000000l
;    ier = 0L
;    max_ntval = _ntval
;    time_si = [1,max_ntval,time_arr_typ,max_ntval]
;    time = make_array(size=time_si)
;    s = CALL_EXTERNAL (libddww,'ddgetaug','ddtbase',  $
;                       ier,dia_ref,tname,1l,max_ntval,time_type,max_ntval,time,dim)
;    if (ier ne 0) then s = CALL_EXTERNAL (libddww,'ddgetaug','xxerror',ier,ctrl,' ')
;    k1 = 1l
;    tbeg = time(k1-1)
;    tmp = max(time(0:dim-1),ntval)
;    tend = time(ntval)
;    ntval = long(ntval+1)
;    k2 = ntval
    k1		= 1L
    k2		= long(_ntval)
    ntval	= long(_ntval)
    tbeg	= _tbeg
    tend	= _tend
    tb_index = (where(ntval eq [obj_size,dim_arr]))(0)
    if (tb_index ne 0) then message,'Zeitbasis nicht als erster Index abgespeichert',/inform
    if(keyword_set(debug)) then begin
	print,'size = ',_ntval,obj_size,ntval
	print, 'tb_index = ',where(ntval eq [obj_size,dim_arr])
	print,'ind=',k1, k2, ntval, 'size=',obj_size, dim_arr, 'tb_index = ',tb_index
	print, obuf
    endif

    if (keyword_set(time_ft) and (obj_size eq ntval)) then begin
       time_ft = float(time_ft)
       tbeg = time_ft(0) > tbeg
       tend = tbeg > time_ft(1) < tend
;       k1 = long((where(time ge tbeg))(0))+1l
;       k2 = long((where(time ge tend))(0))+1l
;       obj_size = ntval  
       if(keyword_set(debug)) then print,'before ddtindex: t=',tbeg,tend,',ind=',k1,k2,ntval
       s = CALL_EXTERNAL (libddww,'ddgetaug','ddtindex', $
                          ier,dia_ref,tname,tbeg,tend,k1,k2)
       if (ier ne 0) then s = CALL_EXTERNAL (libddww,'ddgetaug','xxerror',ier,ctrl,' ')
       ntval = long(k2-k1+1)
       if(keyword_set(debug)) then print,'after ddtindex: t=',tbeg,tend,',ind=',k1,k2,ntval
    endif
;    else begin
;	print, 'NOT reducing time ......'
;	k1 = 1l
;	k2 = obj_size
;    endelse

    time_si = [1,ntval,time_arr_typ,ntval]
    time = make_array(size=time_si)
    s = CALL_EXTERNAL (libddww,'ddgetaug','ddtbase',  $
                       ier,dia_ref,tname,k1,k2,time_type,ntval,time,dim)
    if (ier ne 0) then begin
       s = CALL_EXTERNAL(libddww,'ddgetaug','xxerror',ier, ctrl, ' ')
       goto, close
    endif   
 endif else begin
    time = indgen(obj_size)+1
    ntval = n_elements(time)
 endelse


 ;
 ; Read format type of data array
 ; for signal or signalgroup:
 ;
 format = lonarr(3)
 s = CALL_EXTERNAL(libddww,'ddgetaug','ddobjval', $
                   ier,dia_ref,sig_name,'format',format)
 if (ier ne 0) then s = CALL_EXTERNAL(libddww,'ddgetaug','xxerror',ier, ctrl, ' ')

 ; DDWW data type:
 type = (long(type_val(where(format(0) eq format_val) > 0)))(0)
 ; IDL data type from / for size:
 arr_typ = (array_type(where(format(0) eq format_val) > 0))(0)
 ; text describtion:
 data_format_txt = (format_txt(where(format(0) eq format_val) > 0))(0)

 if(keyword_set(debug)) then begin
	print, format
	print, type
	print, arr_typ
	print, data_format_txt
 endif


 ;
 ;READ DATA
 ;
 case obj_type of
    
    6 : begin  ;signalgroup
        dim_arr = lonarr(3)   
        s = CALL_EXTERNAL (libddww,'ddgetaug','ddobjval', $
                           ier,dia_ref,sig_name,'indices',dim_arr)

        if (keyword_set(indices)) then begin
           indices = long(1 > [indices,1,1,1] < dim_arr)  
;           si = [1,obj_size,arr_typ,obj_size]
           si = [1,ntval,arr_typ,ntval]
           sig = make_array(size=si)
;           s = CALL_EXTERNAL (libddww,'ddgetaug','ddcxsig',    $
;                    ier,dia_ref,sig_name,k1,k2,indices,type,obj_size,sig,dim, $
;                    ncal,phys_dim)
           if (tb_index eq 0) then begin
              ;s = CALL_EXTERNAL (libddww,'ddgetaug','ddccxsig',    $
              ;                   ier,dia_ref,sig_name,k1,k2,indices,type,ntval,sig,dim, $
              ;                   ncal,phys_dim)
              s = CALL_EXTERNAL (libddww,'ddgetaug','ddcxsig',    $
                                 ier,dia_ref,sig_name,k1,k2,indices,type,ntval,sig,dim, $
                                 ncal,phys_dim)
           endif else begin
              tmp = make_array(size=[1,1,arr_typ,1])
              k11 = indices(tb_index-1)
              indi = indices
              for i=k1,k2 do begin
                 indi(tb_index-1) = i                 
                 ;s = CALL_EXTERNAL (libddww,'ddgetaug','ddccxsig',    $
                 ;                   ier,dia_ref,sig_name,k11,k11,indi,type,1L,tmp,dim, $
                 ;                   ncal,phys_dim)
                 s = CALL_EXTERNAL (libddww,'ddgetaug','ddcxsig',    $
                                    ier,dia_ref,sig_name,k11,k11,indi,type,1L,tmp,dim, $
                                    ncal,phys_dim)
                 sig(i-k1) = tmp(0)
              endfor
           endelse
           
        endif else begin   
           n_el_ind = 1
           for i=0,n_elements(dim_arr)-1 do n_el_ind = n_el_ind*dim_arr(i)
           si=[1+n_elements(dim_arr),ntval,dim_arr,arr_typ,ntval*n_el_ind]
           sig = make_array(size=si)
           s = CALL_EXTERNAL (libddww,'ddgetaug','ddcsgrp',ier,    $
           dia_ref,sig_name,k1,k2,type,ntval,sig,dim,ncal,phys_dim)

;           s = CALL_EXTERNAL (libddww,'ddgetaug','ddccsgrp',ier,    $
;                    dia_ref,sig_name,k1,k2,type,ntval,sig,dim,ncal,phys_dim)
        endelse
                    
        end
    
    7 : begin  ;signal
        si = [1,ntval,arr_typ,ntval]
;        si = [1,obj_size,arr_typ,obj_size]
        sig = make_array(size=si)
;        s = CALL_EXTERNAL (libddww,'ddgetaug','ddcsgnl',ier,dia_ref, $
;
;sig_name,k1,k2,type,obj_size,sig,dim,ncal,phys_dim)
        s = CALL_EXTERNAL (libddww,'ddgetaug','ddccsgnl',ier,dia_ref, $
                           sig_name,k1,k2,type,ntval,sig,dim,ncal,phys_dim)
;	my_local_len = 0l
;        s = CALL_EXTERNAL (libddww,'ddgetaug','ddsignal',ier,dia_ref, $
;                           sig_name,k1,k2,type,ntval,sig,my_local_len)
        end                    
    
    else: begin
	message, ' Invalid data type', /continue
	end

 endcase                         

 if (ier eq 555352323) then ier=0  ;warning: no calibration
 if (ier eq 555483394) then ier=0  ;warning: no PARAM_SET found
 if (ier eq 555417858) then ier=0  ;warning: no PARAM_SET found
 if (ier ne 0) then $
    s = CALL_EXTERNAL (libddww,'ddgetaug','xxerror',ier, ctrl, ' ')

close:

;CLOSE SHOTFILE

 ierc = 0l
 s = CALL_EXTERNAL (libddww,'ddgetaug','ddclose',ierc, dia_ref)                   
 if (ierc ne 0) then $
    s = CALL_EXTERNAL (libddww,'ddgetaug','xxerror',ier, ctrl, ' ')
    
 edition=ed   


return
end
;   read_signal.pro - read any shot file signal/-group and its timebase
;
;
;	Wolfgang Suttrop, wls@slcwls
;	E1, Tel. 1466, L6, Zi. 156
;
;	27-Sep-94  wls  extracted from read_cec
;	09-Dec-94  wls	time can be any dimension (follow relations)
;	21-Dec-94  wls	time and area base read optional,
;			making use of DDAINFO and DDAGROUP
;	03-Apr-98  wls  allow for time dependent area base
;	16-Mar-99  wls	caution multiple time bases, absent signal

PRO read_signal,error,diaref,name,t1,t2,data,$
    time=time,area=area,raw=raw

;	error	(<=)	libddww error code
;	diaref	(=>)	ddww file handle, shot file must be open!
;	name	(=>)	name of signal (group) of interest
;	t1,t2	(=>)	desired time range
;	data	(<=)	signal (group) data, memory allocated by 'read_signal'.
;			'data' assumes the dimensions given in shot file
;			exept for the time index which is sized according to
;			the t1, t2 range specified
;	time	(<=)	optional time base, memory allocated inside function
;	area	(<=)	optional area base, memory allocated inside function
;	raw	(=>)	optional keyword: if set (raw=1 or /raw) raw data is
;			read, i.e. no calibration steps are performed. 
;			By default, calibration is performed.

; Shot file access library
;  defsysv, '!libddww', exists=exists
;  if exists eq 0 then defsysv, '!libddww', '/usr/ads/lib/libddww.so', 1

; make sure parameters are of correct type
    if (!VERSION.MEMORY_BITS eq 32) then begin              
	_libddww = '/afs/ipp/aug/ads/lib/@sys/libddww.so'  
    endif else begin                                        
	_libddww = '/afs/ipp/aug/ads/lib64/@sys/libddww8.so'
    endelse


  error = LONG(error)
  diaref = LONG(diaref)

  ctrl = 3L
  k1=0L
  k2=0L
  rleng = 0L
  type = 2L
  ncal = 0L

  rname = string('',format='(a8)')
  buf=lonarr(26)
  rbuf=lonarr(26)
  text=string('',format='(a80)')
  physdim  = string('',format='(a12)')

; get time interval index range

; -------------------------------------------------------
; this doen't work if time is not first dimension !!!
;
;  s = CALL_EXTERNAL (!libddww,'ddgetaug','ddtindex',$
;	error, diaref, name, FLOAT(t1), FLOAT(t2), k1, k2)
;  IF error NE 0 THEN BEGIN
;    print, 'cannot find time interval in signal ', name
;    error = -2L
;    RETURN
;  ENDIF
; -------------------------------------------------------



; get object header of this object for size & type information

  s = CALL_EXTERNAL (!libddww,'ddgetaug','ddobjhdr',$
                      error, diaref, name, buf, text)
  IF error NE 0 THEN GOTO,quit
;debug
;  print,'buf  ', buf

  IF buf(2) EQ -1 THEN BEGIN
    print,'%READ_SIGNAL: signal '+name+' not found'
    error = -1L
    RETURN
  ENDIF

  IF buf(16) LT 1 OR buf(16) GT 4 THEN BEGIN
    print,'%READ_SIGNAL: unsupported dimension ',buf(16),' encountered'
    error = -1L
    RETURN
  ENDIF

  IF buf(0) LT 6 OR buf(0) GT 7 THEN BEGIN
    print,'%READ_SIGNAL: '+name+' must be signal or signal group'
    error = -1L
    RETURN
  ENDIF

  fi=[0L,0L,0L,0L]			; default index range: from start ...
  li=LONG(reverse(buf(18:21)-1))	; ... to actual length in each dim.
  ld=li-fi+1				; length in dimensions
; dim=(SIZE(WHERE(ld GT 1)))(1)		; number of dimensions


; ----------------------------------------------------------------------
; follow all relations and read time base

  ant=1   ; read area base only at one point if no time base found
  ak1=1   ; start with 1st index
  dim=0	  ; dimension teller

  tdim = -1 ; time base dimension
  adim = -1 ; are base dimension

  FOR i=4,11 DO BEGIN
    IF buf(i) GT 0 THEN BEGIN	; is a defined relation
      s = CALL_EXTERNAL (_libddww,'ddgetaug','ddobjname',$
                      error,diaref,buf(i),rname)
      IF error NE 0 THEN GOTO,quit
      s = CALL_EXTERNAL (_libddww,'ddgetaug','ddobjhdr',$
                      error, diaref, rname, rbuf, text)
      IF error NE 0 THEN GOTO,quit
;debug
;      print,'rbuf ', rbuf

; evaluate relations

      CASE rbuf(0) OF
        8: BEGIN		

; read time base array

	     IF tdim LT 0 THEN BEGIN
	       IF(ld(dim) NE rbuf(21)) THEN BEGIN
		  print,'Time Base Dimensionen passen nicht zum Signal'
		  print,'oder das Programm hat noch Macken'
		  error = -2
                  GOTO,quit
	       ENDIF

	       length = LONG (rbuf(21))
	       k1=1L  &  k2= LONG(length)
               time = fltarr(length)
               s = CALL_EXTERNAL (_libddww,'ddgetaug','ddtbase', $
                       error,diaref,rname,k1,k2,type,length,time,rleng)
               IF error NE 0 THEN GOTO,quit

; determine desired subrange

	       wh=WHERE(time GE t1 AND time LE t2)
	       IF wh(0) EQ -1 THEN BEGIN
    	         print,'%READ_SIGNAL: signal ',name,' has time range t=',$
		  time(0),' ... ',time(rleng-1),' s.'
                 print,'There are NO DATA POINTS between ',t1,' s and ',t2,' s'
	         error= -2
	         GOTO,quit	
	       ENDIF

	       IF N_ELEMENTS(wh) NE rleng THEN BEGIN   ; adjust to subrange
		  time = TEMPORARY(time(wh))
	          fi(dim) = MIN(wh)
	          li(dim) = MAX(wh)
;	          ld(dim)=li(dim)-fi(dim)+1
	       ENDIF
	       ant = li(dim)-fi(dim)+1  ; # time indices for area base
	       ak1 = fi(dim)+1          ; start time index for area base
	       tdim = dim
	       dim=dim+1
	    ENDIF
	  END   ; 8=time base

       13: BEGIN

; check area base dimensions
  	     adim=0	 ; area base dimension
  	     FOR j=0,2 DO IF rbuf(18+j) GT 0 THEN BEGIN
		IF ld(dim) NE rbuf(18+j) THEN BEGIN
		  print,'Area Base Dimensionen passen nicht zum Signal'
		  print,'oder das Programm hat noch Macken'
                  error= -2
		  GOTO,quit
		ENDIF
		adim=adim+1
		dim=dim+1
	     ENDIF
	   END   ; 13=area base

        ELSE:  ; do nothing (yet)

      ENDCASE
    ENDIF
  ENDFOR


; ----------------------------------------------------------------------
; follow all relations and read area base

  FOR i=4,11 DO BEGIN
    IF buf(i) GT 0 THEN BEGIN	; is a defined relation
      s = CALL_EXTERNAL (_libddww,'ddgetaug','ddobjname',$
                      error,diaref,buf(i),rname)
      IF error NE 0 THEN GOTO,quit
      s = CALL_EXTERNAL (_libddww,'ddgetaug','ddobjhdr',$
                      error, diaref, rname, rbuf, text)
      IF error NE 0 THEN GOTO,quit

; evaluate relations

      CASE rbuf(0) OF
        8: BEGIN		; ignore time base here
	   END

       13: BEGIN

; read area base
	     IF KEYWORD_SET(area) THEN BEGIN
             print,'read area base'
             
	       IF rbuf(21) LT ant THEN ant = 1
               CASE adim OF 
                 1: BEGIN 
	              areabuf=fltarr(rbuf(18))
	              IF ant GT 1 THEN area=fltarr(ant, rbuf(18))
	            END
		 2: BEGIN
	              areabuf=fltarr(rbuf(18), rbuf(19))
	              IF ant GT 1 THEN area=fltarr(ant, rbuf(18), rbuf(19))
		    END
		 3: BEGIN
	              areabuf=fltarr(rbuf(18), rbuf(19), rbuf(20))
	              IF ant GT 1 THEN area=fltarr(ant, rbuf(18), rbuf(19), rbuf(20))
	 	    END
	       ENDCASE
	       sa=SIZE(areabuf)
	       alength=LONG(sa((SIZE(sa))(1)-1))

               FOR j=0,ant-1 DO BEGIN	; for every time index
	         IF ant EQ 1 THEN ak=1L ELSE ak=LONG(j+ak1)
		 s = CALL_EXTERNAL (_libddww,'ddgetaug','ddagroup', error,$
              	      diaref, name,ak,ak,type,alength,areabuf,rleng)
                 IF error NE 0 THEN GOTO,quit

	         IF ant EQ 1 THEN area=areabuf $
	         ELSE BEGIN
                   CASE adim OF 
                     1: area(j,*)=areabuf
		     2: area(j,*,*)=areabuf
		     3: area(j,*,*,*)=areabuf
	           ENDCASE
	         ENDELSE
	       ENDFOR
	     ENDIF
	   END
        ELSE:  ; do nothing
      ENDCASE
    ENDIF
  ENDFOR



; ------------------------------------------------
; allocate data storage space


  ld(0)=li(0)-fi(0)+1

;  IF KEYWORD_SET(raw) AND $
;       ((buf(0) EQ 7 AND buf(14) EQ 3)   OR   (buf(0) EQ 6 AND buf(4) EQ 0)) $
;       THEN BEGIN
  IF KEYWORD_SET(raw) THEN BEGIN
    CASE buf[14] OF
      2: BEGIN  ; character*1
        type=6L
        string_length = 1L ; one byte long
        CASE buf(16) OF
          1: data=strarr(ld(0))
          2: data=strarr(ld(0), ld(1))
          3: data=strarr(ld(0), ld(1), ld(2))
          4: data=strarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
        data = data + ' '
      END
      3: BEGIN  ; short integer
        type=0L
        CASE buf(16) OF
          1: data=intarr(ld(0))
          2: data=intarr(ld(0), ld(1))
          3: data=intarr(ld(0), ld(1), ld(2))
          4: data=intarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
      END
      4: BEGIN  ; long integer
        type=1L
        CASE buf(16) OF
          1: data=lonarr(ld(0))
          2: data=lonarr(ld(0), ld(1))
          3: data=lonarr(ld(0), ld(1), ld(2))
          4: data=lonarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
      END
      5: BEGIN  ; float
        type=2L
        CASE buf(16) OF
          1: data=fltarr(ld(0))
          2: data=fltarr(ld(0), ld(1))
          3: data=fltarr(ld(0), ld(1), ld(2))
          4: data=fltarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
      END
      6: BEGIN  ; double
        type=3L
        CASE buf(16) OF
          1: data=dblarr(ld(0))
          2: data=dblarr(ld(0), ld(1))
          3: data=dblarr(ld(0), ld(1), ld(2))
          4: data=dblarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
      END
      1794: BEGIN  ; character*8
        type=6L
        string_length = 8L ; eight bytes long
        CASE buf(16) OF
          1: data=strarr(ld(0))
          2: data=strarr(ld(0), ld(1))
          3: data=strarr(ld(0), ld(1), ld(2))
          4: data=strarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
        data = data + STRING('',FORMAT='(a8)')
      END
      3842: BEGIN  ; character*16
        type=6L
        string_length = 16L ; sixteen bytes long
        CASE buf(16) OF
          1: data=strarr(ld(0))
          2: data=strarr(ld(0), ld(1))
          3: data=strarr(ld(0), ld(1), ld(2))
          4: data=strarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
        data = data + STRING('',FORMAT='(a16)')
      END
      7938: BEGIN  ; character*32
        type=6L
        string_length = 32L ; thirty two bytes long
        CASE buf(16) OF
          1: data=strarr(ld(0))
          2: data=strarr(ld(0), ld(1))
          3: data=strarr(ld(0), ld(1), ld(2))
          4: data=strarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
        data = data + STRING('',FORMAT='(a32)')
      END
      12034: BEGIN  ; character*48
        type=6L
        string_length = 48L ; forty eight bytes long
        CASE buf(16) OF
          1: data=strarr(ld(0))
          2: data=strarr(ld(0), ld(1))
          3: data=strarr(ld(0), ld(1), ld(2))
          4: data=strarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
        data = data + STRING('',FORMAT='(a48)')
      END
      16130: BEGIN  ; character*64
        type=6L
        string_length = 64L ; sixty four bytes long
        CASE buf(16) OF
          1: data=strarr(ld(0))
          2: data=strarr(ld(0), ld(1))
          3: data=strarr(ld(0), ld(1), ld(2))
          4: data=strarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
        data = data + STRING('',FORMAT='(a64)')
      END
      18178: BEGIN  ; character*72
        type=6L
        string_length = 72L ; seventy two bytes long
        CASE buf(16) OF
          1: data=strarr(ld(0))
          2: data=strarr(ld(0), ld(1))
          3: data=strarr(ld(0), ld(1), ld(2))
          4: data=strarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
        data = data + STRING('',FORMAT='(a72)')
      END
    ENDCASE
  ENDIF ELSE BEGIN
    CASE buf[14] OF
      6: BEGIN  ; double
        type=3L
        CASE buf(16) OF
          1: data=dblarr(ld(0))
          2: data=dblarr(ld(0), ld(1))
          3: data=dblarr(ld(0), ld(1), ld(2))
          4: data=dblarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
      END
      ELSE: BEGIN  ; float
        type=2L
        CASE buf(16) OF
          1: data=fltarr(ld(0))
          2: data=fltarr(ld(0), ld(1))
          3: data=fltarr(ld(0), ld(1), ld(2))
          4: data=fltarr(ld(0), ld(1), ld(2), ld(3))
        ENDCASE
      END
    ENDCASE
  ENDELSE

; read signal (group) array

  k1=LONG(fi(0)+1)
  k2=LONG(li(0)+1)
  length=LONG(ld(0))
  IF type EQ 6L THEN length = length*string_length ; now in bytes

; Here, I'm having trouble with reading the entire signal group
; ==> try DDxtrsignal (as in ISIS) and ask Annedore!
  IF buf(0) EQ 6 THEN   $      ; signal group
    IF KEYWORD_SET(raw) THEN $
      s = CALL_EXTERNAL (_libddww,'ddgetaug','ddsgroup',$
                     error,diaref,name,k1,k2,type,length,data,rleng)  $
    ELSE  $
      s = CALL_EXTERNAL (_libddww,'ddgetaug','ddcsgrp',$
            error,diaref,name,k1,k2,type,length,data,rleng,ncal,physdim) $

;  IF buf[0] EQ 6 THEN BEGIN    ; signal group
;    temp = data[*,0,0,0]
;    FOR k = 0L, ld[3]-1 DO BEGIN
;      FOR j = 0L, ld[2]-1 DO BEGIN
;        FOR i = 0L, ld[1]-1 DO BEGIN
;          indices = [i+1,j+1,k+1]
;          IF KEYWORD_SET(raw) THEN $
;            s = CALL_EXTERNAL (!libddww,'ddgetaug','ddxtrsignal',$
;                               error,diaref,name,k1,k2,indices, $
;                               type,length,temp,rleng)  $
;          ELSE  $
;            s = CALL_EXTERNAL (!libddww,'ddgetaug','ddcxsig',$
;                               error,diaref,name,k1,k2,indices, $
;                               type,length,temp,rleng,ncal,physdim)
;          data[*,i,j,k] = temp
;        ENDFOR
;      ENDFOR
;    ENDFOR
;  ENDIF $

  ELSE IF buf(0) EQ 7 THEN  $  ; signal
    IF KEYWORD_SET(raw) THEN $
      s = CALL_EXTERNAL (_libddww,'ddgetaug','ddsignal',$
                     error,diaref,name,k1,k2,type,length,data,rleng)  $
    ELSE  $
      s = CALL_EXTERNAL (!libddww,'ddgetaug','ddcsgnl',$
            error,diaref,name,k1,k2,type,length,data,rleng,ncal,physdim)

; Ignore warnings
    IF BYTE(error,2) EQ 1 THEN error = 0L
;  IF    BYTE(error,0) EQ 33 $ ; DD-library
;    AND BYTE(error,1) EQ 27 $ ; ddcsgrp
;    AND BYTE(error,2) EQ  1 $ ; warning
;    AND BYTE(error,3) EQ  2 $ ; "No PARAM_SET found"
;  THEN error = 0L

  IF error NE 0L THEN GOTO,quit

  CASE buf(16) OF
    1: ;
    2: data=TEMPORARY(data(*,fi(1):li(1)))
    3: data=TEMPORARY(data(*,fi(1):li(1),fi(2):li(2)))
    4: data=TEMPORARY(data(*,fi(1):li(1),fi(2):li(2),fi(3):li(3)))
 ENDCASE

quit:
  IF error GT 0 THEN $
    s = CALL_EXTERNAL (_libddww,'ddgetaug','xxerror', error, ctrl, ' ')

END
pro read_data,sig,shot,output,time,trange=trange,diag=diag,signal=signal,ascii=ascii,doplot=doplot,res=res,experiment = experiment,edition=edition

;    if shot eq 32244 and sig eq 'ELM' then begin
;    	experiment = 'MBERN'
;	edition = 2L
;    endif		
    if shot eq 33258 and sig eq 'ELM' then begin
    	experiment = 'SHENDERS'
	edition = 0L
    endif		
    if shot eq 33268 and sig eq 'ELM' then begin
    	experiment = 'SHENDERS'
	edition = 0L
    endif		
    
    ; Set the correct version of libddww

    if (!VERSION.MEMORY_BITS eq 32) then begin              
	_libddww = '/afs/ipp/aug/ads/lib/@sys/libddww.so'  
    endif else begin                                        
	_libddww = '/afs/ipp/aug/ads/lib64/@sys/libddww8.so'
    endelse

    ; Choose the right shotfile

    if ~keyword_set(experiment)then experiment='AUGD'
    if sig eq 'Nimp' then begin
    	diagnostic = 'CES'
	signalname = 'nimp'
    endif	
    if sig eq 'ipsa' then begin
    	diagnostic = 'MAC'
	signalname = 'Ipolsola'
    endif	
    if sig eq 'ipsi' then begin
    	diagnostic = 'MAC'
	signalname = 'Ipolsoli'
    endif	
    if sig eq 'ELM' then begin
    	diagnostic = 'ELM'
	signalname = 'f_ELM'
    endif	
    if sig eq 'ELMi' then begin
    	diagnostic = 'POT'
	signalname = 'ELMi-Ha'
    endif	
    if sig eq 'ELMa' then begin
    	diagnostic = 'POT'
	signalname = 'ELMa-Ha'
    endif	
    if sig eq 'Npuff' then begin
    	diagnostic = 'UVS'
	signalname = 'N_tot'
    endif	
    if sig eq 'Tdiv' then begin
    	diagnostic = 'DDS'
	signalname = 'Tdiv'
    endif	
    if sig eq 'Wmhd' then begin
    	diagnostic = 'GQI'
	signalname = 'Wmhd'
    endif	
    if sig eq 'Dpuff' then begin
    	diagnostic = 'UVS'
	signalname = 'D_tot'
    endif	
    if sig eq 'Nrad' then begin
    	diagnostic = 'GVL'
	signalname = 'Ar0_7504'
    endif	
     if sig eq 'LSD' then begin
    	diagnostic = 'LSD'
	signalname = 'ne-ua3'
    endif	
    
    if keyword_set(diag)then diagnostic=diag
    if keyword_set(signal)then signalname=signal 

    if ~keyword_set(edition)then edition=1L				; Last closed edition

    ; Open the shotfile
    error=0L
    diaref=0L
    date='123456789012345668'		; String with 18 characters
    shot=LONG(shot)
    result=call_external(_libddww,'ddgetaug','ddopen', $
           error,experiment,diagnostic,shot,edition,diaref,date)
    if (error ne 0) then begin		; Print error message
	result=call_external(_libddww,'ddgetaug','xxerror', $
	   error,3L,diagnostic)
	return
    endif

    ; Read the timebase

    type=2L				; real values
    k1=1L				; Read from index 1 to 35000
    k2=35000L				
    time=fltarr(k2)
    leng=0L
    result=call_external(_libddww,'ddgetaug','ddtbase', $
           error,diaref,signalname,k1,k2,type,k2,time,leng)
    if (error ne 0) then begin
	result=call_external(_libddww,'ddgetaug','xxerror', $
	   error,3L,signalname)
	return
    endif
    nt=leng				; Number of values returned

    ; Read the plasmacurrent

    output=fltarr(k2,k2)
    result=call_external(_libddww,'ddgetaug','ddsignal', $
           error,diaref,signalname,k1,k2,type,k2,output,leng)
    if (error ne 0) then begin
	result=call_external(_libddww,'ddgetaug','xxerror', $
	   error,3L,signalname)
	return
    endif

    ; Close the shotfile

    result=call_external(_libddww,'ddgetaug','ddclose', $
           error,diaref)
	IF KEYWORD_SET(trange)THEN BEGIN
		id   = WHERE(time GE trange[0] and time LE trange[1])
		IF id[0] NE -1 THEN BEGIN
			time   = time[id]
			output = output[id]
		ENDIF
	ENDIF		 
	IF KEYWORD_SET(res)THEN BEGIN
		ntime  = (MAX(time)-MIN(time))/res
		tline  = FINDGEN(ntime)*(MAX(time)-MIN(time))/(ntime-1.0)+MIN(time)
		out    = INTERPOL(output,time,tline)
		time   = tline
		output = out
	ENDIF		 
	IF KEYWORD_SET(ascii)THEN BEGIN
		GET_LUN,unit_write & file = 'cview_out.txt'
		OPENW,unit_write,file
		PRINTF,unit_write,'  Time / s    Value'
		FOR i=0,N_ELEMENTS(time)-1 DO PRINTF,unit_write,STRING(time(i),output(i),FORMAT='(E12.4,E12.4)')
		CLOSE,unit_write & FREE_LUN,unit_write
	
	ENDIF	   
	IF KEYWORD_SET(doplot)THEN BEGIN
		PLOT,time,output,XTIT='Time / s',YTIT='Value'
	ENDIF	   
end
PRO user_psym,psym,fill=fill
;   ****************************
;   **** PSYM=1 : Circle    ****
;   **** PSYM=2 : Square    ****
;   **** PSYM=3 : Triangle  ****
;   **** PSYM=4 : Star	    ****
;   **** PSYM=5 : Diamond    ****
;   ****************************
    CASE psym OF
    	1 : Begin
    	    a = FINDGEN(17) * (!PI*2.0/16.0)
    	    USERSYM, cos(A), sin(A), FILL=FILL
        End
        2 : Begin
            USERSYM, [-1,-1,1, 1,-1], [-1, 1,1,-1,-1], FILL=FILL
        End
        3 : Begin
            USERSYM, [-1,0,1,-1], [-1,1,-1,-1], FILL=FILL
        End
        4 : Begin
    	    ang = (360. / 10 * findgen(11) + 90) / !RADEG  ;star angles every 36 deg
    	    r = ang*0
    	    r[2*indgen(6)] = 1.
    	    cp5 = cos(!pi/5.)
    	    r1 = 2. * cp5 - 1. / cp5
    	    r[2*indgen(5)+1] = r1
    	    r = r  / sqrt(!pi/4.) * 2. / (1.+r1)
    	    xarr = r * cos(ang)   &   yarr = r * sin(ang)
    	    USERSYM, xarr, yarr, FILL=FILL
        End
        5 : Begin
            USERSYM, [-1,0,1,0,-1], [0,1,0,-1,0], FILL=FILL
        End
    ENDCASE
    	
END
pro win,num=num,free=free

if ~keyword_set(num)then num=3
if keyword_set(free)then window,/free,xsize=700,ysize=500 else window,num,xsize=700,ysize=500

end
function pos, dim , graph , data , spos=spos

; dim = [3 x 1] - 3 rows of data
; dim = [1 x 3] - 3 columns of data
; dim = [3 x 2] - 3 rows by 2 columns of data
; graph 1 starts from the top left and follows right and then down
; coded so far - [1x1], [2x1], [3x1], [4x1], [5x1], [1x2], [2x2]
if ~keyword_set(spos)then begin
left   = 0.15
right  = 0.95
bottom = 0.15
top    = 0.95
endif else begin
left   = spos[0]
right  = spos[1]
bottom = spos[2]
top    = spos[3]
endelse
   
if dim[1] eq 1 then begin
  if dim[0] eq 1 then begin
    vint = (top - bottom ) / 1.0
    hint = (right - left ) / 1.0
    if graph eq 1 then begin
      position=[left,bottom,right,top]
    endif
  endif
  if dim[0] eq 2 then begin
    vint = (top - bottom ) / 2.0
    hint = (right - left ) / 1.0
    if graph eq 1 then begin
      position=[left,bottom+vint,right,top]
    endif
    if graph eq 2 then begin
      position=[left,bottom,right,bottom+vint]
    endif  
  endif
  if dim[0] eq 3 then begin
    vint = (top - bottom ) / 3.0
    hint = (right - left ) / 1.0
    if graph eq 1 then begin
      position=[left,bottom+2*vint,right,top]
    endif
    if graph eq 2 then begin
      position=[left,bottom+vint,right,bottom+2*vint]
    endif  
    if graph eq 3 then begin
      position=[left,bottom,right,bottom+vint]
    endif  
  endif
  if dim[0] eq 4 then begin
    vint = (top - bottom ) / 4.0
    hint = (right - left ) / 1.0
    if graph eq 1 then begin
      position=[left,bottom+3*vint,right,top]
    endif
    if graph eq 2 then begin
      position=[left,bottom+2*vint,right,bottom+3*vint]
    endif  
    if graph eq 3 then begin
      position=[left,bottom+vint,right,bottom+2*vint]
    endif  
    if graph eq 4 then begin
      position=[left,bottom,right,bottom+vint]
    endif  
  endif
  if dim[0] eq 5 then begin
    vint = (top - bottom ) / 5.0
    hint = (right - left ) / 1.0
    if graph eq 1 then begin
      position=[left,bottom+4*vint,right,top]
    endif
    if graph eq 2 then begin
      position=[left,bottom+3*vint,right,bottom+4*vint]
    endif  
    if graph eq 3 then begin
      position=[left,bottom+2*vint,right,bottom+3*vint]
    endif  
    if graph eq 4 then begin
      position=[left,bottom+vint,right,bottom+2*vint]
    endif  
    if graph eq 5 then begin
      position=[left,bottom,right,bottom+vint]
    endif  
  endif
endif

if dim[1] eq 2 then begin
  if dim[0] eq 1 then begin
    vint = (top - bottom ) / 1.0
    hint = (right - left ) / 2.0
    if graph eq 1 then begin
      position=[left,bottom,left+hint,top]
    endif
    if graph eq 2 then begin
      position=[left+hint,bottom,right,top]
    endif  
  endif
  if dim[0] eq 2 then begin
    vint = (top - bottom ) / 2.0
    hint = (right - left ) / 2.0
    if graph eq 1 then begin
      position=[left,bottom+vint,left+hint,top]
    endif
    if graph eq 2 then begin
      position=[left+hint,bottom+vint,right,top]
    endif  
    if graph eq 3 then begin
      position=[left,bottom,left+hint,bottom+vint]
    endif  
    if graph eq 4 then begin
      position=[left+hint,bottom,right,bottom+vint]
    endif  
  endif
endif



return,position
end

function xz
return,systime(1)
end
 
FUNCTION CROSS,VECTOR1,VECTOR2
  MAG=VECTOR1
  MAG[0]=  VECTOR1[1]*VECTOR2[2]-VECTOR1[2]*VECTOR2[1]
  MAG[1]=-(VECTOR1[0]*VECTOR2[2]-VECTOR1[2]*VECTOR2[0])
  MAG[2]=  VECTOR1[0]*VECTOR2[1]-VECTOR1[1]*VECTOR2[0]
  RETURN,MAG
END  
FUNCTION DOT,VECTOR1,VECTOR2
  IF N_ELEMENTS(VECTOR1) NE N_ELEMENTS(VECTOR2) THEN RETURN,-1 ELSE BEGIN
  MAG=0 & FOR I=0,N_ELEMENTS(VECTOR1)-1 DO MAG=MAG+VECTOR1[I]*VECTOR2[I] 
  RETURN,MAG
  ENDELSE
END
FUNCTION NORM,VECTOR
  SUM=0.0 & FOR I=0,N_ELEMENTS(VECTOR)-1 DO SUM=SUM+VECTOR[I]^2
  RETURN,SQRT(SUM)
END
FUNCTION VECANGLE,VECTOR1,VECTOR2
  ANGLE=ACOS(DOT(VECTOR1,VECTOR2)/(ABS(VECTOR1)*ABS(VECTOR2)))
  RETURN,ANGLE
END  
FUNCTION BINOMIALVEC,VECTOR1,VECTOR2
  VECTOR3=VECTOR1+VECTOR2
  SQUARE=DOT(VECTOR3,VECTOR3)
  EXPANSION=(ABSVEC(VECTOR1))^2+2*(DOT(VECTOR1,VECTOR2)+ABSVEC(CROSS(VECTOR1,VECTOR2)))+(ABSVEC(VECTOR2))^2
  PRINT,SQUARE,EXPANSION
  RETURN,SQUARE
END 

PRO SHIMAGE,DATA,X,Y,RANGE=RANGE,TITLE=TITLE,XTITLE=XTITLE,YTITLE=YTITLE,XRANGE=XRANGE,YRANGE=YRANGE,CHARSIZE=CHARSIZE,NOBAR=NOBAR,XSTYLE=XSTYLE,YSTYLE=YSTYLE
;============================
; ENSURE X-AXIS IS LINEAR
;============================
res=x[1]-x[0]
adas_colors,colors=colors
loadct,5
if ~keyword_set(range)then range=[min(data),max(data)]
IF ~KEYWORD_SET(NOBAR)THEN BEGIN
  res=x[1]-x[0]
  IF ~KEYWORD_SET(RANGE)THEN RANGE=[MIN(DATA),MAX(DATA)]  
  PLOTIMAGE,DATA,RANGE=RANGE,$
  IMGXRANGE=[MIN(X)-res/2.0,MAX(X)+res/2.0],IMGYRANGE=[MIN(Y),MAX(Y)],XTITLE=XTITLE,YTITLE=YTITLE,CHARSIZE=CHARSIZE,$
  NCOLORS=NCOLORS,POSITION=[0.2,0.2,0.9,0.75],XRANGE=XRANGE,YRANGE=YRANGE
  COLORBAR,DIVISIONS=5,RANGE=RANGE,FORMAT='(d5.1)',POSITION=[0.85, 0.2, 0.90, 0.9],CHARSIZE=CHARSIZE,TITLE=TITLE,NCOLORS=NCOLORS
ENDIF ELSE BEGIN
  IF ~KEYWORD_SET(RANGE)THEN RANGE=[MIN(DATA),MAX(DATA)]  
  PLOTIMAGE,DATA,RANGE=RANGE,$
  IMGXRANGE=[MIN(X)-res/2.0,MAX(X)+res/2.0],IMGYRANGE=[MIN(Y),MAX(Y)],TITLE=TITLE,XTITLE=XTITLE,YTITLE=YTITLE,CHARSIZE=CHARSIZE,$
  NCOLORS=NCOLORS,POSITION=[0.2,0.2,0.9,0.9],XRANGE=XRANGE,YRANGE=YRANGE
ENDELSE
END


PRO SHCONTOUR,DATA,X,Y,RANGE=RANGE,XRANGE=XRANGE,YRANGE=YRANGE,XTITLE=XTITLE,YTITLE=YTITLE,CHARSIZE=CHARSIZE,YLOG=YLOG,XLOG=XLOG,XSTYLE=XSTYLE,YSTYLE=YSTYLE

LOADCT,5
IF ~KEYWORD_SET(NLE)THEN NLE=50
adas_colors,colors=colors
CONTOUR,DATA,X,Y,XR=XRANGE,YR=YRANGE,/FILL,NLE=NLE,POSITION=[0.2,0.1,0.9,0.8],YLOG=YLOG,XLOG=XLOG,XSTYLE=XSTYLE,YSTYLE=YSTYLE,XTITLE=XTITLE,YTITLE=YTITLE,back=colors.white,col=colors.black
COLORBAR,DIVISIONS=5,RANGE=[MIN(DATA),MAX(DATA)],FORMAT='(d6.2)',POSITION=[0.9, 0.2, 0.97, 0.9],CHARSIZE=CHARSIZE,col=colors.black

END


FUNCTION MYGAUSS, X, M
RETURN,Y=M[0]*EXP(-((X-M[1])/M[2])^2/2.0)
END

FUNCTION IDSELECT,X,Y
ID=WHERE(ABS(X-Y) EQ MIN(ABS(X-Y)))
RETURN,ID
END
;********************************************************************************
;********************************************************************************
PRO ERR_PLOT, X, Y, DELY, COLPICK=COLPICK, THICK=THICK,YLOG=YLOG
X=REFORM(X) & Y=REFORM(Y) & DELY=REFORM(DELY)
IF KEYWORD_SET(YLOG)THEN BEGIN
 FOR I=0,N_ELEMENTS(X)-1 DO BEGIN
   OPLOT,[X[I],X[I]],[(Y[I]-DELY[I])>1E-10,(Y[I]+DELY[I])>1E-10],THICK=THICK,COL=COLPICK
   DIFF = 0.005*(!X.crange[1]-!X.crange[0])
   ;IF I EQ 0 THEN DIFF=X[I+1]-X[I] ELSE DIFF=X[I]-X[I-1]
   ;DIFF=DIFF*0.2
   OPLOT,[X[I]-DIFF,X[I]+DIFF],[Y[I]-DELY, Y[I]-DELY],THICK=THICK,COL=COLPICK
   OPLOT,[X[I]-DIFF,X[I]+DIFF],[Y[I]+DELY, Y[I]+DELY],THICK=THICK,COL=COLPICK
 ENDFOR
ENDIF ELSE BEGIN
 FOR I=0,N_ELEMENTS(X)-1 DO BEGIN
   OPLOT,[X[I],X[I]],[Y[I]-DELY[I],Y[I]+DELY[I]],THICK=THICK,COL=COLPICK
   ;IF I EQ 0 THEN DIFF=X[I+1]-X[I] ELSE DIFF=X[I]-X[I-1]
   ;DIFF=DIFF*0.1
   DIFF = 0.005*(!X.crange[1]-!X.crange[0])
   OPLOT,[X[I]-DIFF,X[I]+DIFF],[Y[I]-DELY[I], Y[I]-DELY[I]],THICK=THICK,COL=COLPICK
   OPLOT,[X[I]-DIFF,X[I]+DIFF],[Y[I]+DELY[I], Y[I]+DELY[I]],THICK=THICK,COL=COLPICK
 ENDFOR
ENDELSE 
END

;********************************************************************************
FUNCTION STRAIGHTFIT, X, Y, W=W, M=M, C=C, DELC=DELC, DELM=DELM
;================================================================================
; Purpose: Perform a straight line fit to
;          data
;
; Equations taken from: 
; G. L. Squires, Practical Physics, 4th Ed, Cambridge University Press, pp 40
;
; Author  : Stuart Henderson
; Date    : Sept 2013
; Version : 1.0
;================================================================================
N=N_ELEMENTS(X)                        
IF ~KEYWORD_SET(W)THEN W=FLTARR(N)+1.0 
D=TOTAL(W*X^2)-(1.0/TOTAL(W))*TOTAL(W*X)^2
E=TOTAL(W*X*Y)-(1.0/TOTAL(W))*TOTAL(W*X)*TOTAL(W*Y)
F=TOTAL(W*Y^2)-(1.0/TOTAL(W))*TOTAL(W*Y)^2
M=E/D
C=MEAN(Y)-M*MEAN(X)
DELM2=(1./(N-2))*(D*F-E^2)/D^2 & DELM=SQRT(DELM2)
DELC2=(1./(N-2))*(D/TOTAL(W)+(MEAN(X))^2)*(D*F-E^2)/D^2 & DELC=SQRT(DELC2)
RETURN,M*X+C
END
;********************************************************************************
FUNCTION MAP3D, X, Y, Z, DATA, NEWX, NEWY, NEWZ
;================================================================================
; Purpose: Perform 3D interpolation
;
; Author  : Stuart Henderson
; Date    : Sept 2013
; Version : 1.0
;================================================================================ 
X_ID = (NEWX-MIN(X))/(MAX(X)-MIN(X)) * (N_ELEMENTS(X)-1) 
Y_ID = (NEWY-MIN(Y))/(MAX(Y)-MIN(Y)) * (N_ELEMENTS(Y)-1) 
Z_ID = (NEWZ-MIN(Z))/(MAX(Z)-MIN(Z)) * (N_ELEMENTS(Z)-1) 
RETURN,REFORM(INTERPOLATE(DATA,X_ID,Y_ID,Z_ID,/GRID))
END
FUNCTION MAP_2D, X, Y,DATA, NEWX, NEWY
;================================================================================
; Purpose: Perform 3D interpolation
;
; Author  : Stuart Henderson
; Date    : Sept 2013
; Version : 1.0
;================================================================================ 
X_ID = (NEWX-MIN(X))/(MAX(X)-MIN(X)) * (N_ELEMENTS(X)-1) 
Y_ID = (NEWY-MIN(Y))/(MAX(Y)-MIN(Y)) * (N_ELEMENTS(Y)-1) 
RETURN,REFORM(INTERPOLATE(DATA,X_ID,Y_ID,/GRID))
END
;********************************************************************************
PRO MAKEPS,FILENAME=FILENAME,CLOSE=CLOSE,AR=AR,XSIZE=XSIZE,YSIZE=YSIZE,LANDSCAPE=LANDSCAPE,PORTRAIT=PORTRAIT
;================================================================================
; Purpose: Produce a postscript file
;
; Author  : Stuart Henderson
; Date    : Sept 2013
; Version : 1.0
;================================================================================ 
IF ~KEYWORD_SET(CLOSE)THEN BEGIN
    MUL=!P.MULTI
    IF ~KEYWORD_SET(FILENAME)THEN FILENAME='Plot.ps'
    IF ~KEYWORD_SET(AR)THEN AR=1.5
    IF ~KEYWORD_SET(XSIZE)THEN XSIZE=7
    IF ~KEYWORD_SET(YSIZE)THEN YSIZE=6
    !P.FONT=0
    !P.THICK=6
    !X.THICK=6
    !Y.THICK=6
    !Z.THICK=6
    !P.CHARSIZE=1.5
    SET_PLOT,'ps'
    DEVICE,COLOR=1,XSIZE=XSIZE,YSIZE=YSIZE,/INCHES,BITS_PER_PIXEL=64,FILE=FILENAME,$
           FONT_SIZE=11,LANDSCAPE=LANDSCAPE,PORTRAIT=PORTRAIT,/ENCAPSULATED
ENDIF ELSE BEGIN
    DEVICE,/CLOSE_FILE
    SET_PLOT,'X'
    !P.FONT=-1
    !P.THICK=1.5
    !X.THICK=1.5
    !Y.THICK=1.5
    !Z.THICK=1.5
    !P.CHARSIZE=1.5
ENDELSE    
END
;********************************************************************************
PRO COLORBAR, BOTTOM=bottom, CHARSIZE=charsize, COLOR=color, DIVISIONS=divisions, $
   FORMAT=format, POSITION=position, MAXRANGE=maxrange, MINRANGE=minrange, NCOLORS=ncolors, $
   TITLE=title, VERTICAL=vertical, TOP=top, RIGHT=right, MINOR=minor, $
   RANGE=range, FONT=font, TICKLEN=ticklen, _EXTRA=extra, INVERTCOLORS=invertcolors, $
   TICKNAMES=ticknames, REVERSE=reverse, ANNOTATECOLOR=annotatecolor, XLOG=xlog, YLOG=ylog, $
   NODISPLAY=nodisplay

    compile_opt idl2

    ; Return to caller on error.
    On_Error, 2

    ; Save the current plot state.
    bang_p = !P
    bang_x = !X
    bang_Y = !Y
    bang_Z = !Z
    bang_Map = !Map

    ; Are scalable pixels available on the device?
    IF (!D.Flags AND 1) NE 0 THEN scalablePixels = 1 ELSE scalablePixels = 0

    ; Which release of IDL is this?
    thisRelease = Float(!Version.Release)

    ; Check and define keywords.
    IF N_ELEMENTS(ncolors) EQ 0 THEN BEGIN

       ; Most display devices to not use the 256 colors available to
       ; the PostScript device. This presents a problem when writing
       ; general-purpose programs that can be output to the display or
       ; to the PostScript device. This problem is especially bothersome
       ; if you don't specify the number of colors you are using in the
       ; program. One way to work around this problem is to make the
       ; default number of colors the same for the display device and for
       ; the PostScript device. Then, the colors you see in PostScript are
       ; identical to the colors you see on your display. Here is one way to
       ; do it.

       IF scalablePixels THEN BEGIN
          oldDevice = !D.NAME

             ; What kind of computer are we using? SET_PLOT to appropriate
             ; display device.

          thisOS = !VERSION.OS_FAMILY
          thisOS = STRMID(thisOS, 0, 3)
          thisOS = STRUPCASE(thisOS)
          CASE thisOS of
             'MAC': SET_PLOT, thisOS
             'WIN': SET_PLOT, thisOS
             ELSE: SET_PLOT, 'X'
          ENDCASE

          ; Here is how many colors we should use.
          ncolors = !D.TABLE_SIZE
          SET_PLOT, oldDevice
        ENDIF ELSE ncolors = !D.TABLE_SIZE
    ENDIF
    IF N_ELEMENTS(bottom) EQ 0 THEN bottom = 0B
    IF N_ELEMENTS(charsize) EQ 0 THEN charsize = !P.Charsize
    IF N_ELEMENTS(format) EQ 0 THEN format = '(I0)'
    IF N_ELEMENTS(color) EQ 0 THEN color = !P.Color
    minrange = (N_ELEMENTS(minrange) EQ 0) ? 0. : Float(minrange)
    maxrange = (N_ELEMENTS(maxrange) EQ 0) ? Float(ncolors) : Float(maxrange)
    IF N_ELEMENTS(ticklen) EQ 0 THEN ticklen = 0.2
    IF N_ELEMENTS(minor) EQ 0 THEN minor = 2
    IF N_ELEMENTS(range) NE 0 THEN BEGIN
       minrange = Float(range[0])
       maxrange = Float(range[1])
    ENDIF
    IF N_ELEMENTS(divisions) EQ 0 THEN divisions = 6
    IF N_ELEMENTS(font) EQ 0 THEN font = !P.Font
    IF N_ELEMENTS(title) EQ 0 THEN title = ''
    xlog = Keyword_Set(xlog)
    ylog = Keyword_Set(ylog)

    ; You can't have a format set *and* use ticknames.
    IF N_ELEMENTS(ticknames) NE 0 THEN format = ""

    ; If the format is NOT null, then format the ticknames yourself.
    ; Can't assume minrange is less than maxrange.
    IF (xlog XOR ylog) EQ 0 THEN BEGIN
        IF format NE "" THEN BEGIN
           IF minrange LT maxrange THEN BEGIN
               step = (maxrange - minrange) / divisions
               levels = minrange > (Indgen(divisions+1) * step + minrange) < maxrange
               IF StrPos(StrLowCase(format), 'i') NE -1 THEN levels = Round(levels)
               ticknames = String(levels, Format=format)
               format = "" ; No formats allowed in PLOT call now that we have ticknames.
	       
           ENDIF ELSE BEGIN
               step = (minrange - maxrange) / divisions
               levels = maxrange > (Indgen(divisions+1) * step + maxrange) < minrange
               levels = Reverse(levels)
               IF StrPos(StrLowCase(format), 'i') NE -1 THEN levels = Round(levels)
               ticknames = String(levels, Format=format)
               format = "" ; No formats allowed in PLOT call now that we have ticknames.
           ENDELSE
        ENDIF
    ENDIF

    IF KEYWORD_SET(vertical) THEN BEGIN
       bar = REPLICATE(1B,20) # BINDGEN(ncolors)
       IF Keyword_Set(invertcolors) THEN bar = Reverse(bar, 2)
       IF N_ELEMENTS(position) EQ 0 THEN BEGIN
          position = [0.88, 0.1, 0.95, 0.9]
       ENDIF ELSE BEGIN
          IF position[2]-position[0] GT position[3]-position[1] THEN BEGIN
             position = [position[1], position[0], position[3], position[2]]
          ENDIF
          IF position[0] GE position[2] THEN Message, "Position coordinates can't be reconciled."
          IF position[1] GE position[3] THEN Message, "Position coordinates can't be reconciled."
       ENDELSE
    ENDIF ELSE BEGIN
       bar = BINDGEN(ncolors) # REPLICATE(1B, 20)
       IF Keyword_Set(invertcolors) THEN bar = Reverse(bar, 1)
       IF N_ELEMENTS(position) EQ 0 THEN BEGIN
          position = [0.1, 0.88, 0.9, 0.95]
       ENDIF ELSE BEGIN
          IF position[3]-position[1] GT position[2]-position[0] THEN BEGIN
             position = [position[1], position[0], position[3], position[2]]
          ENDIF
          IF position[0] GE position[2] THEN Message, "Position coordinates can't be reconciled."
          IF position[1] GE position[3] THEN Message, "Position coordinates can't be reconciled."
       ENDELSE
    ENDELSE

    ; Scale the color bar.
     bar = BYTSCL(bar, TOP=(ncolors-1) < (255-bottom)) + bottom

     IF Keyword_Set(reverse) THEN BEGIN
       IF Keyword_Set(vertical) THEN bar = Reverse(bar,2) ELSE bar = Reverse(bar,1)
     ENDIF

    ; Get starting locations in NORMAL coordinates.
    xstart = position[0]
    ystart = position[1]

    ; Get the size of the bar in NORMAL coordinates.
    xsize = (position[2] - position[0])
    ysize = (position[3] - position[1])

    ; Display the color bar in the window. Sizing is
    ; different for PostScript and regular display.
    IF scalablePixels THEN BEGIN

       TV, bar, xstart, ystart, XSIZE=xsize, YSIZE=ysize, /Normal

    ENDIF ELSE BEGIN

       bar = CONGRID(bar, CEIL(xsize*!D.X_VSize), CEIL(ysize*!D.Y_VSize))

       ; Decomposed color off if device supports it.
       CASE  StrUpCase(!D.NAME) OF
            'X': BEGIN
                IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
                Device, Decomposed=0
                ENDCASE
            'WIN': BEGIN
                IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
                Device, Decomposed=0
                ENDCASE
            'MAC': BEGIN
                IF thisRelease GE 5.2 THEN Device, Get_Decomposed=thisDecomposed
                Device, Decomposed=0
                ENDCASE
            ELSE:
       ENDCASE

       TV, bar, xstart, ystart, /Normal

       ; Restore Decomposed state if necessary.
       CASE StrUpCase(!D.NAME) OF
          'X': BEGIN
             IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
             ENDCASE
          'WIN': BEGIN
             IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
             ENDCASE
          'MAC': BEGIN
             IF thisRelease GE 5.2 THEN Device, Decomposed=thisDecomposed
             ENDCASE
          ELSE:
       ENDCASE

    ENDELSE

    ; Annotate the color bar.
    IF N_Elements(annotateColor) NE 0 THEN $
      color = FSC_Color(annotateColor, color, NODISPLAY=Keyword_Set(nodisplay))

    IF KEYWORD_SET(vertical) THEN BEGIN

       IF KEYWORD_SET(right) THEN BEGIN

          PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=1, $
             YTICKS=divisions, XSTYLE=1, YSTYLE=9, $
             POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
             XTICKFORMAT='(A1)', YTICKFORMAT='(A1)', YMINOR=minor, _EXTRA=extra, $
             YTICKNAME=ticknames, FONT=font, YLOG=ylog

          AXIS, YAXIS=1, YRANGE=[minrange, maxrange], YTICKFORMAT=format, YTICKS=divisions, $
             YTICKLEN=ticklen, YSTYLE=1, COLOR=color, CHARSIZE=charsize, $
             FONT=font, YTITLE=title, _EXTRA=extra, YMINOR=minor, YTICKNAME=ticknames, YLOG=ylog

       ENDIF ELSE BEGIN

          PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=1,  $
             YTICKS=divisions, YSTYLE=1, XSTYLE=1, YTITLE=title, $
             POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
             XTICKFORMAT='(A1)', YTICKFORMAT=format, YMinor=minor, _EXTRA=extra, $
             YTICKNAME=ticknames, YLOG=ylog, YTICKLEN=ticklen, FONT=font

       ENDELSE

    ENDIF ELSE BEGIN

       IF KEYWORD_SET(top) THEN BEGIN

          PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=divisions, $
             YTICKS=1, XSTYLE=9, YSTYLE=1, $
             POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
             YTICKFORMAT='(A1)', XTICKFORMAT='(A1)', XTICKLEN=ticklen, $
             XRANGE=[minrange, maxrange], FONT=font, XMINOR=minor,_EXTRA=extra, $
             XTICKNAME=ticknames, XLOG=xlog

          AXIS, XTICKS=divisions, XSTYLE=1, COLOR=color, CHARSIZE=charsize, $
             XTICKFORMAT=format, XTICKLEN=ticklen, XRANGE=[minrange, maxrange], XAXIS=1, $
             FONT=font, XTITLE=title, _EXTRA=extra, XCHARSIZE=charsize, XMINOR=minor, $
             XTICKNAME=ticknames, XLOG=xlog

       ENDIF ELSE BEGIN

          PLOT, [minrange,maxrange], [minrange,maxrange], /NODATA, XTICKS=divisions, $
             YTICKS=1, XSTYLE=1, YSTYLE=1, TITLE=title, $
             POSITION=position, COLOR=color, CHARSIZE=charsize, /NOERASE, $
             YTICKFORMAT='(A1)', XTICKFORMAT=format, XTICKLEN=ticklen, $
             XRANGE=[minrange, maxrange], FONT=font, XMinor=minor, _EXTRA=extra, $
             XTICKNAME=ticknames, XLOG=xlog

        ENDELSE

    ENDELSE

    ; Restore the previous plot and map system variables.
    !P = bang_p
    !X = bang_x
    !Y = bang_y
    !Z = bang_z
    !Map = bang_map

END
;+
; NAME:
;   PLOTIMAGE
;
; AUTHOR:
;   Craig B. Markwardt, NASA/GSFC Code 662, Greenbelt, MD 20770
;   craigm@lheamail.gsfc.nasa.gov
;
; PURPOSE:
;   Displays an image via a "PLOT"-like interface.
;
; CALLING SEQUENCE:
;   PLOTIMAGE, img, [xrange=xrange,] [yrange=yrange,] ...
;
; DESCRIPTION:
;
;   PLOTIMAGE displays an image (or slice of an image) on the current
;   graphics device.  The syntax is very similar to the PLOT command,
;   in the sense that an XRANGE and YRANGE for the plot can be
;   specified.  
;
;   PLOTIMAGE keeps separate the notions of the image coordinate
;   system and the displayed coordinate system, which allows any input
;   image to be "cropped," "zoomed," or "flipped."
;
;   PLOTIMAGE allows the user to express image extents in physical
;   units rather than pixel units.
;
;   The image coordinate system specifies the physical coordinates of
;   original image data, IMG.  The image is considered to be a 2D
;   array (IMG = ARRAY(NX,NY)), where the values are attached to the
;   midpoint of each geometric pixel.  The image has NX columns and NY
;   rows.  Physical coordinates are attached to each pixel by using
;   the IMGXRANGE and IMGYRANGE keywords.  The IMGXRANGE keyword is a
;   two-element array specifying the "left" and "right" boundaries of
;   the image pixels in physical units; the IMGYRANGE keyword
;   specifies the "top" and "bottom" boundaries of the image.  This is
;   illustrated in Figure 1 for a simplified case.
;
;                                   ___
;         +-----------+-----------+  ^  IMGYRANGE[1]
;         |           |           |  |
;         | IMG[0,1]  | IMG[1,1]  |  |
;         |     +     |     +     |  |
;         |           |           |  |
;         |           |           |  |
;         +-----------+-----------+  |
;         |           |           |  |
;         | IMG[0,0]  | IMG[1,0]  |  |
;         |     +     |     +     |  |
;         |           |           |  |
;         |           |           |  v
;         +-----------+-----------+ ___ IMGYRANGE[0]
;        |                         |
;        |<----------------------->|
;        IMGXRANGE[0]   IMGXRANGE[1]
;
;      Figure 1.  Simplified example of a 2x2 input image,
;      demonstrating that IMG[*,*] values refer to the pixel
;      mid-points, and that IMGXRANGE and IMGYRANGE ranges specify the
;      physical coordinates of the outer edges of the image extent in
;      X and Y, respectively.
; 
;
;   The displayed plot coordinate system is entirely independent of
;   the native image coordinates.  Users can set up the plot scale
;   using any combination of {X,Y}RANGE, {X,Y}STYLE and/or {X,Y}LOG,
;   as they would for any IDL plot, using physical units.  The input
;   image will then be overlayed on this coordinate system.
;
;   If the displayed plot coordinates are narrower than the native
;   image coordinates, then the displayed portion of the image will be
;   cropped to fit.  If the displayed coordinates are wider than the
;   native image coordinates, then the image will be displayed with
;   blank spaces on either side (see Figure 2).  A mirror "flip" is
;   also possible in X and/or Y, if XRANGE or YRANGE are specified in
;   reverse order.
;                                                 ___
;      +---------------------------------------+   ^
;      |            ___                        |   |
;      |             ^  +---------------+      |   |
;      |             |  |               |      |   |
;      |             |  |               |      |   |
;      |    IMGYRANGE|  |     IMG       |      |   | YRANGE
;      |             |  |               |      |   |
;      |             v  |               |      |   |
;      |            ___ +---------------+      |   |
;      |               |<-- IMGXRANGE -->|     |   |
;      |                                       |   v
;      +---------------------------------------+  ___
;     |<-------------   XRANGE   -------------->|
;
;     Figure 2.  Example of an image whose native image coordinates
;     are embedded in a wider plot display range.
;
;   The standard [XY]STYLE keywords can be used to style either axis.
;   However at the very least [XY]STYLE=1 is always implied, i.e. the
;   plot limits exactly obey the [XY]RANGE keyword values.
;
;   If XLOG or YLOG are set, then the image is assumed to be sampled
;   on a logarithmic grid, and logarithmic axes are displayed
;   accordingly.  PLOTIMAGE does not attempt to resample the image
;   from linear scale to logarithmic scale, or reverse.
;
;   Psuedocolor images may be of any type, but must rescaled to a byte
;   range by using the RANGE keyword.  By default the color range used
;   in the rescaling operation is 0 to !D.N_COLORS - 3B.  The extra
;   two color values are reserved for the background and default pen
;   colors.  This behavior can be adjusted by specifying the BOTTOM
;   and/or NCOLORS keywords.
;
;   Truecolor images must always be of type BYTE and one of their
;   dimensions must have 3 elements, corresponding to the three color
;   planes of the image.
;
;
; INPUTS:
;
;   IMG - Array to be displayed.  For single-plane images (i.e.,
;         pseudocolor), the image must be two dimensional and of any
;         real numeric type.  For images that are not of BYTE type,
;         the RANGE keyword must be supplied, and then PLOTIMAGE will
;         rescale the image values to a byte range.
;
;         An image declared as ARRAY(NX,NY) will be NX pixels in the
;         x-direction and NY pixels in the y-direction.  The image is
;         resampled to fill the desired display region (and optionally
;         smoothed).
;
;         For three-plane images (i.e., truecolor) the image must be
;         of type BYTE.  One of the dimensions of the array must have
;         three elements.  Hence it must be one of BYTARR(NX, NY, 3),
;         BYTARR(NX, 3, NY) or BYTARR(3, NX, NY).  The 3-element
;         dimension is recognized automatically.
;
; OPTIONAL INPUTS:
;   NONE
;
; INPUT KEYWORD PARAMETERS:
;
;   IMGXRANGE, IMGYRANGE - Each is a two component vector that
;                          describes the X and Y position of the outer
;                          edges of the first and last pixels.
;                          Default: IMGXRANGE = [0,NX]
;                                   IMGYRANGE = [0,NY]
;
;   XRANGE, YRANGE - Each is a two component vector that specifies the
;                    X and Y plot ranges, respectively.  These values
;                    are not required to coincide with IMG[XY]RANGE.
;                    Default: XRANGE=IMGXRANGE
;                             YRANGE=IMGYRANGE
;
;   POSITION - Position of the inner plot window in the standard
;              graphics keyword format.  Overrides PANEL and SUBPANEL.
;
;   INTERP - if set, interpolate (smooth) the image before displaying.
;            This keyword applies to the screen displays.  For printed
;            images that are coarser than MIN_DPI, the image is
;            implicitly interpolated regardless of INTERP.
;
;   PRESERVE_ASPECT - if set, preserve the aspect ratio of the
;                     original image (in pixels).  The result will be
;                     the largest image that fits in the display
;                     region while maintaining approximately square
;                     pixels.  However, PIXEL_ASPECT_RATIO overrides
;                     PRESERVE_ASPECT.  The POSITION keyword will be
;                     reset upon output to the ultimate image
;                     position.
;                     DEFAULT: not set (image will fill POSITION rectangle)
;
;   PIXEL_ASPECT_RATIO - The ratio of width to height for each pixel.
;                        If specified, then the image will be scaled
;                        so that each pixel has the specified aspect
;                        ratio.  If not specified, then the image will
;                        be scaled independently in X and Y in order
;                        to fill the POSITION rectangle.  NOTE: If you
;                        want to change the overall image aspect
;                        ratio, then use the POSITION keyword.
;                  DEFAULT: undefined (image will fill POSITION rectangle)
;
;   MIN_DPI - if printing, the minimum dot-per-inch pixel resolution
;             for the resulting image.  Output images that would be
;             coarser than this value are resampled to have a
;             resolution of at least MIN_DPI, and smoothed.  Some
;             common resolutions are: screen, 90 dpi; dot matrix, 72
;             dpi; laser printer 300-600 dpi.  Note that large values
;             of MIN_DPI will produce very large output files.
;             Default: 0 (i.e., the output image will not be smoothed)
;
;   RANGE - a two element vector.  If the image is single plane (i.e.,
;           pseudocolor) the input image can be of any real numeric
;           type, and then must be rescaled into byte range with this
;           keyword.  In contrast, truecolor images must always be of
;           type BYTE.  Values are scaled into byte range with the
;           following statement:
;              RESULT = BYTSCL(INPUT, MIN=RANGE(0), MAX=RANGE(1), $
;                              TOP=NCOLORS-1) + BOTTOM
;           so that pixels with an intensity RANGE(0) are set to
;           BOTTOM; those with RANGE(1) are set to the maximum color.
;           Default: no range scaling occurs (and the image must hence
;                    be of type BYTE -- otherwise an error occurs)
;
;   NCOLORS - number of color table values be used in the byte
;             rescaling operation.
;             Default: !D.N_COLORS - BOTTOM - 1 (for default pen color)
;
;   BOTTOM - bottom-most value of the color table to be used in the
;            byte rescaling operation.
;            Default: 1 (for default background color)
;
;   NOERASE - If set, the display is not erased before graphics
;             operations.
;
;   NODATA - If set, the image is not actually displayed, but
;            coordinate axes may be drawn.
;
;   NOAXES - An attempt is made to render the image without coordinate
;            axes.  However, it's usually more straightforward to set 
;            XSTYLE=4 or YSTYLE=4, which is the standard IDL way to
;            disable coordinate axes.
;
;   ORDER - same interpretation as the !ORDER system variable; 
;           if ORDER=0, then the first pixel is drawn in the lower
;           left corner; if ORDER=1, then the first pixel is drawn in
;           the upper left corner.
;           Default: 0
;
;
;   PANEL, SUBPANEL - An alternate way to more precisely specify the
;                     plot and annotation positions.  See SUBCELL.
;
;   PLOTIMAGE will pass other keywords directly to the PLOT command
;   used for generating the plot axes.  XSTYLE=1 and YSTYLE=1 are
;   enforced.
;
; OUTPUTS:
;   NONE
;
; PROCEDURE:
;
; EXAMPLE:
;
;   This example constructs an image whose values are found by
;       z(x,y) = cos(x) * sin(y)
;   and x and y are in the range [-2,2] and [4,8], respectively.
;   The image is then plotted, with appropriate axes.
;
;   x = findgen(20)/5. - 2. + .1   ; 0.1 = half-pixel
;   y = findgen(20)/5. + 4. + .1
;   zz = cos(x) # sin(y)
;   imgxrange = [-2.,2.]           ; extend to pixel edges
;   imgyrange = [4.,8.]
;   plotimage, bytscl(zz), imgxrange=imgxrange, imgyrange=imgyrange
;
;   This second example plots the same image, but with a plot range
;   much larger than the image's.
;
;   xr=[-10.,10]
;   yr=[-10.,10]
;   plotimage, bytscl(zz), imgxrange=imgxrange, imgyrange=imgyrange, $
;      xrange=xr, yrange=yr
;
; SEE ALSO:
;
;   OPLOTIMAGE, BYTSCL
;
; EXTERNAL SUBROUTINES:
;
;   SUBCELL, DEFSUBCELL
;
; MODIFICATION HISTORY:
;   Written, CM, 1997
;   Correct various one-off problems, 02 Feb 1999, CM
;   Made self-contained with some pre-processing, 17 Oct 1999, CM
;   Corrected bug in newly introduced CONGRID functions, 18 Oct 1999, CM
;   Correct behavior with no POSITION keyword, 17 Nov 1999, CM
;   Simplified axis plotting, 17 Nov 1999, CM
;   Use _EXTRA keyword in first PLOT, but with blank TITLEs, 11 Jan
;     2000, CM
;   Correct implementation of X/YSTYLE in first PLOT, 11 Feb 2000, CM
;   Correct CONGRID implementation (small effect when enlarging most
;     images), 14 Feb 2000, CM
;   Major changes: 19 Apr 2000
;      - now handle decomposed color, automatic color mapping via
;        RANGE, and 24-bit multiplane images
;      - new PRESERVE_ASPECT keyword to keep square pixels
;      - removed legacy TVIMAGE code
;      - smoothing is more configurable, esp. for printers, but is not
;        done by default; more printers are supported
;   Corrected INTERPOLATE behavior (thanks to Liam Gumley
;     <Liam.Gumley@ssec.wisc.edu>), other minor tweaks, CM 20 Apr 2000
;   Added ability to use PRESERVE_ASPECT with POSITION, PANEL or
;     SUBPANEL keywords CM 20 Oct 2000
;   Oops, a typo is now fixed, CM 23 Oct 2000
;   Add fix for MacIntoshes and DECOMPOSED color, Tupper, 02 Aug 2001
;   Better behavior with fractional pixels (ie, when the image pixels
;     are very large compared to the screen), 23 Aug 2001
;   Add support for Z buffer, CM, 20 Oct 2002
;   Memory conservation: use REVERSE() to reverse IMG; rewrote
;     PLOTIMAGE_RESAMP to rescale entire image instead of each color plane
;     separately.  Jeff Guerber, 30 July 2003
;   Add PIXEL_ASPECT_RATIO keyword, 22-25 Nov 2005
;   Check for the case of an 1xNXxNY 3D image and treat it as a 2D
;     image.  The "1" dimension can be anywhere, CM, 03 Sep 2006
;   Add the ORDER keyword parameter, CM, 20 Mar 2007
;   Enable XLOG and YLOG keywords, for logarithmic axes;
;     doesn't actually resample the image from linear<->log, CM
;     21 Jan 2009
;   Documentation, CM, 21 Jan 2009
;   Allow reverse color scale, CM, 13 Nov 2010
;
;   $Id: plotimage.pro,v 1.15 2010/11/13 09:53:39 cmarkwar Exp $
;
;-
; Copyright (C) 1997-2001,2003,2005,2006,2007,2009,2010 Craig Markwardt
; This software is provided as is without any warranty whatsoever.
; Permission to use, copy, modify, and distribute modified or
; unmodified copies is granted, provided this copyright and disclaimer
; are included unchanged.
;-
;%insert HERE
;%include subcell.pro
;%include defsubcell.pro

; Utility routine to resample an image
;
;  IMAGE - image data ARRAY(NX,NY,BDEPTH)
;  NX,NY - original X,Y image size
;  BDEPTH- original image depth
;  NEWX, NEWY- desired X,Y image size
;  INTERP - if set, then use bilinear interpolation, otherwise nearest neighbor
function defsubcell, default

  if n_elements(default) EQ 0 then default = [-1.,-1,-1,-1]
  mysubcell = default
  defaultsubpos = [ 0.08, 0.08, 0.95, 0.95 ]

  iwh = where(mysubcell LT 0, ict)
  if ict GT 0 then $
    mysubcell(iwh) = defaultsubpos(iwh)

  return, mysubcell
end

function subcell, subpos, position, margin=margin

  ;; Default value for subposition
  if n_elements(subpos) EQ 0 then mysubpos = [-1.,-1,-1,-1] $
  else mysubpos = subpos

  ;; Default value for position - full screen
  if n_elements(position) EQ 0 then position = [0.,0.,1.,1.]

  ;; Get margins if necessary
  if keyword_set(margin) EQ 1 OR n_elements(subpos) EQ 0 then $
    mysubpos = defsubcell(mysubpos)

  ;; Compute new window position
  x0 = position(0)
  y0 = position(1)
  dx = position(2)-position(0)
  dy = position(3)-position(1)

  newsubpos = reform(mysubpos * 0, 4)
  newsubpos([0,2]) = x0 + dx * mysubpos([0,2])
  newsubpos([1,3]) = y0 + dy * mysubpos([1,3])

  return, newsubpos
end


function plotimage_resamp, image, nx, ny, bdepth, newx, newy, interp=interp

  ;; Sometimes the final dimension is lost.  Put it back
  image = reform(image, nx, ny, bdepth, /overwrite)

  ;; Correct interpolation
  srx = float(nx)/newx * findgen(newx) - 0.5 + 0.5*(float(nx)/newx)
  sry = float(ny)/newy * findgen(newy) - 0.5 + 0.5*(float(ny)/newy)
  srz = indgen(bdepth)
  if keyword_set(interp) then $
    return, interpolate(image, srx, sry, srz, /grid)

  ;; Simple nearest neighbor interpolation
  return, interpolate(image, round(srx), round(sry), srz, /grid)
end

pro plotimage_pos, xrange0, imgxrange0, imgxsize, xreverse, srcxpix, imgxpanel, $
                   logscale=logscale, $
                   quiet=quiet, status=status, pixtolerance=pixtolerance

  if keyword_set(logscale) then begin
     if min(xrange0) LE 0 OR min(imgxrange0) LE 0 then $
        message, ('ERROR: if XLOG or YLOG is set, then the image boundary cannot '+$
                  'cross or touch zero.  Did you forget to set IMGXRANGE or IMGYRANGE?')
     xrange    = alog10(xrange0)
     imgxrange = alog10(imgxrange0)
  endif else begin
     xrange    = xrange0
     imgxrange = imgxrange0
  endelse

  if n_elements(pixtolerance) EQ 0 then pixtolerance = 1.e-2
  status = 0
  ;; Decide if image must be reversed
  xreverse = 0
  if double(xrange(1)-xrange(0))*(imgxrange(1)-imgxrange(0)) LT 0 then begin
      xreverse = 1
      imgxrange = [imgxrange(1), imgxrange(0)]
  endif

  srcxpix  = [ 0L, imgxsize-1 ]
  ;; Size of one x pix
  dx = double(imgxrange(1) - imgxrange(0)) / imgxsize

  if min(xrange) GE max(imgxrange) OR max(xrange) LE min(imgxrange) then begin
      message, 'WARNING: No image data in specified plot RANGE.', /info, $
        noprint=keyword_set(quiet)
      return
  endif

  ;; Case where xrange cuts off image at left
  if (xrange(0) - imgxrange(0))/dx GT 0 then begin
      offset = double(xrange(0)-imgxrange(0))/dx
      if abs(offset-round(offset)) LT pixtolerance then $
        offset = round(offset)
      srcxpix(0) = floor(offset)
      froffset = offset - floor(offset)
      if abs(froffset) GT pixtolerance then begin
          xrange = double(xrange)
          xrange(0) = imgxrange(0) +dx*srcxpix(0)
      endif
  endif

  ;; Case where xrange cuts off image at right
  if (xrange(1) - imgxrange(1))/dx LT 0 then begin
      offset = double(xrange(1)-imgxrange(0))/dx
      if abs(offset-round(offset)) LT pixtolerance then $
        offset = round(offset)
      srcxpix(1) = ceil(offset) - 1
      froffset = offset - ceil(offset)
      if abs(froffset) GT pixtolerance then begin
          xrange = double(xrange)
          srcxpix(1) = srcxpix(1) < (imgxsize-1)
          xrange(1) = imgxrange(0) + dx*(srcxpix(1)+1)
      endif
  endif

  imgxpanel = [0., 1.]
  if (xrange(0) - imgxrange(0))/dx LT 0 then $
    imgxpanel(0) = (imgxrange(0) - xrange(0))/(xrange(1)-xrange(0))
  if (xrange(1) - imgxrange(1))/dx GT 0 then $
    imgxpanel(1) = (imgxrange(1) - xrange(0))/(xrange(1)-xrange(0))

  status = 1
  return
end

;; Main program
pro plotimage, img0, xrange=xrange0, yrange=yrange0, $
               imgxrange=imgxrange0, imgyrange=imgyrange0, $
               xlog=xlog, ylog=ylog, $
               position=position, panel=panel, subpanel=subpanel, $
               xstyle=xstyle, ystyle=ystyle, title=title, $
               interp=interp0, quiet=quiet, dither=dither, $
               preserve_aspect=paspect, pixel_aspect_ratio=asprat, $
               min_dpi=min_dpi, order=order, $
               ncolors=ncolors0, bottom=bottom0, range=range, $
               noerase=noerase0, nodata=nodata, noaxes=noaxes, $
               pixtolerance=pixtolerance, _EXTRA=extra

  ;; Return to user when an error is encountered
  on_error, 2

  ;; Usage message
  if n_params() EQ 0 then begin
      message, 'PLOTIMAGE, image, xrange=, yrange=, imgxrange=, imgyrange=,..', $
               /info
      return
  endif

  ;; Must have a byte-scaled image already
  imgsize  = size(img0)

  ;; Make sure windowing exists (borrowed from IMDISP)
  if ((!d.flags and 256) ne 0) and (!d.window lt 0) then begin
      window, /free, /pixmap
      wdelete, !d.window
  endif

  ;; Parameter checking
  if n_elements(ystyle) EQ 0 then ystyle = 0L
  if n_elements(xstyle) EQ 0 then xstyle = 0L
  if keyword_set(nodata) then mynodata = 1 else mynodata = 0
  if n_elements(pixtolerance) EQ 0 then pixtolerance = 1.e-2
  if n_elements(title) EQ 0 then title = ''
  if n_elements(min_dpi) EQ 0 then min_dpi = 0
  interp = keyword_set(interp0)
  noerase = keyword_set(noerase0)
  imgpanel = [0., 0., 1., 1.]

  ;; Default handling of color table stuff
  if n_elements(bottom0) EQ 0 then bottom0 = 1B
  bottom = byte(bottom0(0)) < 255B
  dncolors = min([!d.n_colors, !d.table_size, 256])
  if n_elements(ncolors0) EQ 0 then ncolors0 = dncolors - 1 - bottom
  ;; Make sure color table values are in bounds
  ncolors = floor(ncolors0(0)) < 256
  if bottom + ncolors GT 256 then ncolors = 256 - bottom

  ;; Image size and dimensions
  nimgdims  = imgsize(0)
  imgtype   = imgsize(nimgdims+1)
  if nimgdims LT 2 OR nimgdims GT 3 then begin
      message, 'ERROR: image must have 2 or 3 dimensions'
  endif

  if nimgdims EQ 2 then begin
      ;; Two dimensional image is pseudo color
      img = img0

      ONE_CHANNEL_IMAGE:
      imgxsize = imgsize(1)
      imgysize = imgsize(2)
      bdepth = 1

      if imgtype NE 1 then begin
          if n_elements(range) LT 2 then $
            message, 'ERROR: non-byte image must be scaled with RANGE keyword'
          if range(0) LE range(1) then begin
              img = bytscl(img, min=range(0), max=range(1), top=ncolors-1B) $
                + bottom
          endif else begin
              ;; Reverse color scheme
              img = bytscl(img, min=range(1), max=range(0), top=ncolors-1B)
              img = ncolors-1B-img + bottom
          endelse
      endif
      img = reform(img, imgxsize, imgysize, bdepth, /overwrite)
  endif else begin
      wh = where(imgsize(1:3) EQ 1, ct)
      if ct GT 0 then begin
          imgxsize = 1
          imgysize = 1
          j = 0
          for i = 1, 3 do if imgsize(i) NE 1 then begin
              if j EQ 0 then imgxsize = imgsize(i) else imgysize = imgsize(i)
              j = j + 1
          endif
          img = reform(img0, imgxsize, imgysize)
          imgsize = size(img)

          goto, ONE_CHANNEL_IMAGE
      endif else begin
          ;; Three dimensional image has three planes
          wh = where(imgsize(1:3) EQ 3, ct)
          if imgtype NE 1 then $
            message, 'ERROR: true color image must of type byte'
          if ct EQ 0 then $
            message, ('ERROR: True color image must have 3 elements '+$
                      'in one of its dimensions')
          truedim = wh(0)
          
          ;; Shuffle the data so planes are interleaved ...
          case truedim of
              0: img = transpose(img0, [1,2,0]) ;; ... from pixels interleaved
              1: img = transpose(img0, [0,2,1]) ;; ... from rows interleaved
              2: img = img0                 ;; ... by straight copying
          end

          imgsize = size(img)
          imgxsize = imgsize(1)
          imgysize = imgsize(2)
          bdepth = imgsize(3)

      endelse

  endelse

  ;; By default, we have no info about the image, and display the
  ;; whole thing
  if n_elements(imgxrange0) LT 2 then imgxrange = [ 0., imgxsize ] $
  else imgxrange = 0. + imgxrange0(0:1)
  if n_elements(xrange0) LT 2 then xrange = imgxrange $
  else xrange = 0. + xrange0(0:1)

  status = 0
  plotimage_pos, xrange, imgxrange, imgxsize, xreverse, srcxpix, imgxpanel, $
    quiet=keyword_set(quiet), status=status, pixtolerance=pixtolerance, $
    logscale=xlog
  if status EQ 0 then mynodata = 1 $
  else imgpanel([0,2]) = imgxpanel

  ;; By default, we have no info about the image, and display the
  ;; whole thing
  if n_elements(imgyrange0) LT 2 then imgyrange = [ 0., imgysize ] $
  else imgyrange = 0. + imgyrange0(0:1)
  if n_elements(yrange0) LT 2 then yrange = imgyrange $
  else yrange = 0. + yrange0(0:1)
  if keyword_set(order) then yrange = [yrange(1), yrange(0)]

  status = 0
  plotimage_pos, yrange, imgyrange, imgysize, yreverse, srcypix, imgypanel, $
    quiet=keyword_set(quiet), status=status, pixtolerance=pixtolerance, $
    logscale=ylog
  if status EQ 0 then mynodata = 1 $
  else imgpanel([1,3]) = imgypanel

  ;; Dimensions of output image in pixels
  nx = srcxpix(1)-srcxpix(0)+1
  ny = srcypix(1)-srcypix(0)+1

  ;; Create a coordinate system by plotting with no data or axes
  if n_elements(position) EQ 0 AND n_elements(panel) EQ 0 AND $
    n_elements(subpanel) EQ 0 then begin

      ;; If PANEL/SUBPANEL is not given, then plot once to set up
      ;; axes, despite NOAXES
      plot, xrange, yrange, noerase=noerase, /nodata, $
        xstyle=xstyle OR 5, ystyle=xstyle OR 5, xlog=xlog, ylog=ylog, $
        xrange=xrange, yrange=yrange, xtitle='', ytitle='', title='', $
        _EXTRA=extra

      ;; Retrieve axis settings
      xwindow = !x.window
      ywindow = !y.window

      subpanel1 = [xwindow(0), ywindow(0), xwindow(1), ywindow(1)]
      imgposition = subcell(imgpanel, subpanel1)
      position = subpanel1

  endif else begin

      ;; Construct the plot size from panel info.  Default is full-screen
      if NOT keyword_set(noerase) then erase
      if n_elements(position) GE 4 then begin
          imgposition = subcell(imgpanel, position)
      endif else begin
          if n_elements(panel) LT 4 then panel = [0.0,0.0,1.0,1.0]
          if n_elements(subpanel) LT 4 then subpanel = [-1., -1, -1, -1]
          subpanel = defsubcell(subpanel)

          imgposition = subcell(subcell(imgpanel, subpanel), panel)
          position = subcell(subpanel, panel)
      endelse

      xwindow = position([0,2])
      ywindow = position([1,3])

  endelse

  ;; If the aspect is to be preserved then we need to recompute the
  ;; position after considering the image size.  Since we have already
  ;; computed the outer envelope of the image from either the POSITION
  ;; or PANEL, or from the plot window itself, we can now go to the
  ;; logic which estimates the aspect-corrected size.

  if (keyword_set(paspect) OR n_elements(asprat) GT 0) AND $
    nx GT 0 AND ny GT 0 then begin

      if n_elements(asprat) EQ 0 then asprat1 = 1.0 $
      else asprat1 = asprat(0) + 0.

      ;; If we are preserving the aspect, then re-plot after scaling
      ;; the POSITION

      imgaspect = float(ny)/float(nx)/asprat1
      dispaspect = (ywindow(1)-ywindow(0))*!d.y_vsize $
        / ((xwindow(1)-xwindow(0))*!d.x_vsize)

      ;; Compute the new image dimensions
      if imgaspect GT dispaspect then begin
          x0 = total(xwindow)/2
          dx = (ywindow(1)-ywindow(0))*!d.y_vsize/(imgaspect*!d.x_vsize)
          xwindow = x0 + dx*[-0.5,0.5]
      endif else begin
          y0 = total(ywindow)/2
          dy = (xwindow(1)-xwindow(0))*!d.x_vsize*imgaspect/!d.y_vsize
          ywindow = y0 + dy*[-0.5,0.5]
      endelse

      subpanel1 = [xwindow(0), ywindow(0), xwindow(1), ywindow(1)]
      imgposition = subcell(imgpanel, subpanel1)
      position = subpanel1

      ;; Replot to regain coordinate system
      plot, xrange, yrange, /noerase, /nodata, $
        xstyle=xstyle OR 5, ystyle=xstyle OR 5, xlog=xlog, ylog=ylog, $
        xrange=xrange, yrange=yrange, xtitle='', ytitle='', title='', $
        position=position, _EXTRA=extra

  endif

  ;; Draw the image data
  if NOT keyword_set(mynodata) then begin

      ;; Reverse X- or Y- directions if necessary
      if xreverse then $
        srcxpix = imgxsize - 1 - [srcxpix(1), srcxpix(0)]
      if yreverse then $
        srcypix = imgysize - 1 - [srcypix(1), srcypix(0)]

      ;; Extract relevant image elements
      img = (temporary(img))(srcxpix(0):srcxpix(1), srcypix(0):srcypix(1),*)
      img = reform(img, nx, ny, bdepth, /overwrite)

      ;; Complete the extraction, if reversed
      if xreverse then begin
          img = reverse(img, 1, /overwrite)
          img = reform(img, nx, ny, bdepth, /overwrite)
      endif
      if yreverse then begin
          img = reverse(img, 2, /overwrite)
          img = reform(img, nx, ny, bdepth, /overwrite)
      endif

      ;; Compute the image position on screen in pixels
      x0 = round(imgposition(0) * !d.x_vsize)
      y0 = round(imgposition(1) * !d.y_vsize)
      dx = round((imgposition(2) - imgposition(0)) * !d.x_vsize) > 1
      dy = round((imgposition(3) - imgposition(1)) * !d.y_vsize) > 1

      ;; Decide which output type
      windowing = (!d.name EQ 'WIN') OR (!d.name EQ 'MAC') OR (!d.name EQ 'X')
      printing = (!d.name EQ 'PRINTER') OR (!d.flags AND 1) NE 0

      ;; Decide whether to resample the image
      rescaling = (windowing OR (!d.name EQ 'Z')) $
        AND ((dx NE nx) OR (dy NE ny))

      ;; If printing, and the printed resolution of the image will be
      ;; too coarse, then we should resample and interpolate
      dpi = min([nx*!d.x_px_cm/dx, ny*!d.y_px_cm/dy]*2.54) ; d.p.i. of image
      dxsize = dx & dysize = dy
      if printing AND (dpi LT min_dpi(0)) then begin
          dx = round(min_dpi(0)*dx/(2.54*!d.x_px_cm)) > nx
          dy = round(min_dpi(0)*dy/(2.54*!d.y_px_cm)) > ny
          interp = 1
          rescaling = 1
      endif

      ;; Rescale the image if needed
      if rescaling then begin
          img = plotimage_resamp(temporary(img), nx, ny, bdepth, $
            dx, dy, interp=interp)
          img = reform(img, dx, dy, bdepth, /overwrite)
      endif

      ;; Generic printer device
      if !d.name EQ 'PRINTER' then begin
          if bdepth EQ 3 then begin
              device, /true_color
              tv, img, x0, y0, xsize=dxsize, ysize=dysize, true=3
          endif else begin
              device, /index_color
              tv, img, x0, y0, xsize=dxsize, ysize=dysize
          endelse
          goto, DONE_IMG
      endif

      ;; Devices with scalable pixels
      if (!d.flags AND 1) NE 0 then begin
          if bdepth EQ 3 then begin
              tvlct, r, g, b, /get
              loadct, 0, /silent
              tv, img, x0, y0, xsize=dxsize, ysize=dysize, true=3
              tvlct, r, g, b
          endif else begin
              tv, img, x0, y0, xsize=dxsize, ysize=dysize
          endelse
          goto, DONE_IMG
      endif

      ;; Get visual depth (in bytes) and decomposed state
      decomposed0 = 0
      vdepth = 1
      version = float(!version.release)
      if windowing then begin

          ;; Visual depth
          if version GE 5.1 then begin
              device, get_visual_depth=vdepth
              vdepth = vdepth / 8
          endif else begin
              if !d.n_colors GT 256 then vdepth = 3
          endelse

          ;; Decomposed state
          if vdepth GT 1 then begin
              if version GE 5.2 then device, get_decomposed=decomposed0
              if bdepth EQ 3 then    device, decomposed=1 $
              else                   device, decomposed=0
          endif
      endif

      ;; If visual is 8-bit but image is 24-bit, then quantize
      if vdepth LE 1 AND bdepth EQ 3 then begin
          img = color_quan(temporary(img), 3, r, g, b, colors=ncolors-1, $
                           dither=keyword_set(dither)) + bottom
          tvlct, r, g, b, bottom
          bdepth = 1
      endif

      ;; Put the image
      if bdepth EQ 3 then tv, img, x0, y0, true=3 $
      else                tv, img, x0, y0

      ;; Restore the decomposed state
      if windowing then begin
          if vdepth GT 1 then device, decomposed=decomposed0
          ;; Tupper supplies following work-around for MacIntoshes
          if (!d.name EQ 'MAC') then tv, [0], -1, -1
      endif
  endif

  ;; Plot the axes if requested
  DONE_IMG:
  if NOT keyword_set(noaxes) then begin
      if n_elements(xrange) EQ 0 then begin
          if n_elements(imgxrange) GT 1 then xrange=imgxrange $
          else xrange = [0L, imgxsize]
      endif
      if n_elements(yrange) EQ 0 then begin
          if n_elements(imgyrange) GT 1 then yrange=imgyrange $
          else yrange = [0L, imgysize]
      endif

      plot, xrange, yrange, /noerase, /nodata, /normal, $
        xrange=xrange, yrange=yrange, xlog=xlog, ylog=ylog, $
        xstyle=xstyle OR 1, ystyle=ystyle OR 1, title=title, $
        position=position, _EXTRA=extra
  endif

  return
end
pro oband,x,y,err,border=border,_extra=_extra,xlog=xlog,ylog=ylog,norm=norm,transparent=transparent
;+
; ROUTINE:   oband
;
; USEAGE:    oband,x,y1,y2,border=border,color=color,spacing=spacing,$
;                  fill_pattern=fill_pattern,line_fill=line_fill,$
;                  pattern=pattern,orientation=orientation
;                  
;
; PURPOSE:   Over plot shaded band on an x-y plot.  Shaded region is
;            between y1 and y2 as a function of x.  This is useful for
;            indicating an error bounds on an x-y plot.  
;
; INPUTS:
;
; x          vector of x values (data coordinates)
; y1         vector of lower y values (data coordinates)
; y2         vector of upper y values (data coordinates)
; 
;            NOTE: y1 need not be smaller than y2. The shaded region
;                  always extends from y1 to y2 no matter which one is
;                  greater.
;
;            NOTE: If y1 or y2 is a scalor it will be used internally
;                  as if it were a constant value vector of the same
;                  length as x.
;
; OPTIONAL KEYWORDS:
;
;   border
;        if set, draw a border around the filled region, the numerical
;        value of BORDER is the color index used to draw the border
;    
;
; Keywords recognized by the POLYFILL procedure:
;
;   COLOR
;        color index used to fill region
;
;   LINE_FILL
;        Set this keyword to indicate that polygons are 
;        to be filled with parallel lines, rather than using 
;        solid or patterned filling methods.When using the 
;        line-drawing method of filling, the thickness, line-
;        style, orientation, and spacing of the lines may be 
;        specified with keywords.
;        
;   PATTERN
;        A rectangular array of pixels giving the fill 
;        pattern. If this keyword parameter is omitted, POLY-
;        FILL fills the area with a solid color. The pattern 
;        array may be of any size; if it is smaller than the 
;        filled area the pattern array is cyclically 
;        repeated.
;
;   SPACING
;        The spacing, in centimeters, between the parallel
;        lines used to fill polygons.
;        
;   ORIENTATION
;        Orientation angle of lines used to fill region.
;        
;Graphics Keywords Accepted
;
;        See Chapter 2, Graphics Keywords and System Variables, 
;        for the description of graphics and ploting keywords not 
;        listed above. CHANNEL CLIP DATA DEVICE LINE NOCLIP NORMAL 
;        T3D THICK Z.
; 
; EXAMPLE:
;
;        x=indgen(200)
;        y=randf(200,3)
;        y1=y+.1*(1.+randf(200,2)^2)
;        y2=y-.1*(1.+randf(200,2)^2)
;        plot,x,y
;        oband,x,y1,y2,color=100
;        oplot,x,y
;
;        plot,x,y
;        oband,x,y,0,/line_fill,orien=45,border=100
;
;
; DISCUSSION: If y1 and y2 are more than 1 element, make sure their
;             array lengths are the same as x.
;     
;
;  author:  Paul Ricchiazzi                            12apr93
;           Institute for Computational Earth System Science
;           University of California, Santa Barbara
;-

if keyword_set(norm)then begin
    y1=y
    y2=err
endif else begin   
    y1=y+err
    y2=y-err
end
ny1=n_elements(y1)
ny2=n_elements(y2)
nx=n_elements(x)

if nx eq 0 then begin & xhelp,'oband' & return & end

if ny1 ne nx and ny1 ne 1 then message,'array length mismatch: y1, x'
if ny2 ne nx and ny2 ne 1 then message,'array length mismatch: y2, x'

if ny1 eq 1 then yy1=replicate(y1,nx) else yy1=y1
if ny2 eq 1 then yy2=replicate(y2,nx) else yy2=y2

xxx=[x,reverse(x)]
yyy=[yy1,reverse(yy2)]

if keyword_set(ylog)then begin
    ymax=10^!y.crange(1)
    ymin=10^!y.crange(0)
endif else begin
    ymax=!y.crange(1)
    ymin=!y.crange(0)
endelse
if keyword_set(xlog)then begin
    xmax=10^!x.crange(1)
    xmin=10^!x.crange(0)
endif else begin
    xmax=!x.crange(1)
    xmin=!x.crange(0)
endelse
yyy=(ymin >= yyy) < ymax
xxx=(xmin >= xxx) < xmax
polyfill,xxx,yyy,_extra=_extra,transparent=transparent

if n_elements(border) ne 0 then oplot,xxx,yyy,color=border

end


pro obandx,x,y,err,border=border,_extra=_extra,xlog=xlog,ylog=ylog,norm=norm,transparent=transparent
;+
; ROUTINE:   obandx
;-

if keyword_set(norm)then begin
    x1=x
    x2=err
endif else begin   
    x1=x+err
    x2=x-err
end
nx1=n_elements(x1)
nx2=n_elements(x2)
ny=n_elements(y)

if ny eq 0 then begin & xhelp,'oband' & return & end

if nx1 ne ny and nx1 ne 1 then message,'array length mismatch: x1, y'
if nx2 ne ny and nx2 ne 1 then message,'array length mismatch: x2, y'

if nx1 eq 1 then xx1=replicate(x1,ny) else xx1=x1
if nx2 eq 1 then xx2=replicate(x2,ny) else xx2=x2

yyy=[y,reverse(y)]
xxx=[xx1,reverse(xx2)]

if keyword_set(ylog)then begin
    ymax=10^!y.crange(1)
    ymin=10^!y.crange(0)
endif else begin
    ymax=!y.crange(1)
    ymin=!y.crange(0)
endelse
if keyword_set(xlog)then begin
    xmax=10^!x.crange(1)
    xmin=10^!x.crange(0)
endif else begin
    xmax=!x.crange(1)
    xmin=!x.crange(0)
endelse
yyy=(ymin >= yyy) < ymax
xxx=(xmin >= xxx) < xmax
polyfill,xxx,yyy,_extra=_extra,transparent=transparent

if n_elements(border) ne 0 then oplot,xxx,yyy,color=border

end


pro progress,percentin,bigstep=bigstepin,smallstep=smallstepin, $
                 reset=reset,last=last,label=label,noeta=noeta, $
                 frequency=frequency,progstring=pstring, $
                 noprint=noprint,norecurse=norecurse, $
                 minval=minvalin,maxval=maxvalin

;+
;NAME:
;     PROGRESS
;PURPOSE:
;     Prints a progress summary in percent done and time remaining as
;     a process runs.  Call this routine multiple times from the
;     running process while updating the percent done parameter.
;CATEGORY:
;CALLING SEQUENCE:
;     progress,percent
;INPUTS:
;     percent = percent finished (in range 0 to 100, unless minval and
;               maxval are set)
;OPTIONAL INPUT PARAMETERS:
;KEYWORD PARAMETERS
;     /reset = this must be set for the first call to set up variables.
;     /last = the process is done and this is the last call.  Optional
;             for the last call.  Makes sure the printout goes all the
;             way to 100%, does a final carriage return, and resets
;             some variables.
;     label = string to print at the front of the progress line.  Only
;             used on the first call when /reset is set.
;     bigstep = percentage multiple to print the percent done as a
;               number (def = 25).  Only used on the first call when
;               /reset is set.  Integer.
;     smallstep = percentage multiple to print a dot (def = 5).  Only
;                 used on the first call when /reset is set.  Integer.
;     minval = percent value at process start (default = 0.0).  Only
;                 used on the first call when /reset is set.  Integer.
;     maxval = percent value at process end (default = 100.0).  Only
;                 used on the first call when /reset is set.  Integer.
;     /noeta = Do not print the estimated time to completion.  This
;              feature depends on your terminal accepting 8b as a
;              backspace character. If this does not work, the
;              formatting will be messed up.  So, if your formatting
;              is messed up, set this keyword to turn the feature off.
;     frequency = If set, update the estimated time to completion
;                 if at least 'frequency' seconds have passed since
;                 the last update.  The time is always printed when a
;                 dot or number is printed, as well. If set to 2, for
;                 example, update approximately every two seconds, etc. 
;                 If set to 0, update on every call to progress.
;                 The default (not set) is to update the time only
;                 when a dot or number is printed. 
;OUTPUTS:
;COMMON BLOCKS:
;SIDE EFFECTS:
;RESTRICTIONS:
;     Assumes that nothing else is printed between calls.  If it is,
;     the formatting will be messed up.  If you change the size of the
;     terminal during the process, the formatting may also be messed
;     up if it might have gone over the length of a line.  The
;     counting is done with integers, so you can't count, for example,
;     from 0. to 1. with dots printed at intervals of 0.1, regardless
;     of how you set minval and maxval.
;PROCEDURE:
;     Call the routine multiple times as progress is made.  It will
;     print numbers and dots to indicate the progress, e.g.
;     Label:  0 .... 25 .... 50 .... 75 .... 100  | Time=00:00:00
;     Also prints an estimated time to completion, HH:MM:SS
;EXAMPLE:
;     progress,0.0,/reset,label='Progress (%)'
;     for i=0,n-1 do begin
;        ; your processing goes here
;        progress,100.*float(i+1.0)/n
;     endfor
;     progress,100.,/last
;     Progress (%):  0 .... 25 .... 50 .... 75 .... 100  | Time=00:00:00
;
;     Or, you can use the minval and maxval keywords to change the
;     range of the counting:
;
;     progress,0.,/reset,label='test',bigstep=128,smallstep=32,maxval=512
;     for i=0,511 do begin
;        wait,0.1
;        progress,i+1
;     endfor
;     progress,512.,/last
;     test:  0 ... 128 ... 256 ... 384 ... 512  | Time=00:00:51
;MODIFICATION HISTORY:
;     T. Metcalf 2005-Jan-06
;     2005-Jan-10 Added frequency keyword.
;     2005-Jan-12 Move time to end of final string so that it does not
;                 move.
;     2005-Jan-13 If /last is set, the time printed is the total elapsed
;                 time.
;     2005-Jan-19 Added minval and maxval keywords.
;     2005-Feb-01 Check for rpercent eq 0 in eta calculation.
;     2005-Feb-11 Remove strcompress around user label.
;-

common progress_private,bigstep,smallstep,lastpercent,starttime,lasttime, $
                        ncharacters,progstring,ttysize,minval,maxval

if n_elements(frequency) GE 1 then freq = float(frequency[0])>0.0 else freq = 0.0

if keyword_set(reset) OR not keyword_set(bigstep) or not keyword_set(smallstep) or $
   n_elements(lastpercent) LE 0 then $
   doreset = 1 else doreset = 0
if n_elements(lastpercent) GT 0 and not keyword_set(doreset) then $
   if lastpercent LT minval then doreset = 1

if keyword_set(doreset)  then begin
   if not keyword_set(bigstepin) then bigstep=25 $
   else bigstep = round(bigstepin)
   if not keyword_set(smallstepin) then smallstep=5 $
   else smallstep = round(smallstepin)
   if n_elements(minvalin) LE 0 then minval = 0.0 else minval = float(minvalin)
   if n_elements(maxvalin) LE 0 then maxval = 100.0 else maxval = float(maxvalin)
   if smallstep GT bigstep then smallstep=bigstep

   if NOT keyword_set(norecurse) then begin ; Get full output string for counting
      progress,100.0,/reset,/last,/noprint,progstring=ps,label=label, $
               /norecurse,/noeta,bigstep=bigstepin,smallstep=smallstepin, $
               minval=minvalin, maxval=maxvalin
      ncharacters = strlen(ps) ; the length of the full progress string
      test = ''
      ; maxch is the longest possible string computed from
      ; ncharacters and the longest time string we're likely need
      maxch =  ncharacters + strlen(' | Time=00000:00:00') 
      ttysize = maxch  ; Maximum required tty width.
      for i=1,maxch-1 do begin 
         ; Figure out how wide the tty is.  There has got to be
         ; a better way to do this!!!
         ; String switches to an array when the tty line would be full.
         test = string('a',test)  ; add one character
         if (size(test))[0] then begin ; scalar string or string array?
            ttysize = i ; if we never get here, the tty is bigger than we need
            break
         endif
      endfor
      ; Better way(?), but works only for some unix flavors,
      ; so do don't want to do this.
      ; The IDL ioctl command could also be used, but it is at least
      ; as system dependent.
      ; spawn,['stty','size'],result,/noshell 
      ; ttysize=splitstr(result,' ')
      ; ttysize = ttysize[n_elements(ttysize)-1]
   endif

   lastpercent = round(minval)-1.0   ; must be exactly round(minval)-1 at the start
   starttime = systime(1)
   lasttime = starttime
   progstring = ''
   if keyword_set(label) then begin
      ;s = strcompress(label+': ')
      s = label+': '
      if NOT keyword_set(noprint) then print,s,format='(a,$)'
      progstring = progstring + s
   endif
endif   ; end of initialization

percent = (float(percentin) < maxval) > minval

if percent LT lastpercent then return 

if keyword_set(last) then npercent = round(maxval) else npercent = round(percent)
nlastpercent = round(lastpercent)

; Print whatever dots and numbers we need to

nprinted = 0
if npercent GT nlastpercent then begin
   for ipct=nlastpercent+1,npercent do begin
      if ipct MOD bigstep EQ 0 then begin
         s = ' '+strcompress(string(ipct),/remove_all)+' '
         if NOT keyword_set(noprint) then print,s,format='(a,$)' 
         progstring = progstring + s
         nprinted = nprinted + strlen(s)
      endif else if ipct MOD smallstep EQ 0 then begin
         s = '.'
         if NOT keyword_set(noprint) then print,s,format='(a,$)'
         progstring = progstring + s
         nprinted = nprinted + strlen(s)
      endif
   endfor
endif

; Print the estimated time left

if (keyword_set(nprinted) OR $
    n_elements(frequency) GE 1 OR $
    keyword_set(last)) AND $
   NOT keyword_set(doreset) AND $ ; no time has passed, can't get eta
   NOT keyword_set(noeta) then begin
   ; Compute & print estimated time to completion
   thistime = systime(1)
   if (thistime - lasttime) GE freq OR $
      keyword_set(nprinted) OR $
      keyword_set(last) then begin
      rpercent = float(percent-minval)/float(maxval-minval)
      if keyword_set(last) OR rpercent EQ 0.0 then $
         eta = (thistime-starttime) $ ; /last -> total elapsed time
      else $
         eta = (1.0-rpercent)*(thistime-starttime)/rpercent
      h = fix(eta / 3600.0)
      eta = eta - h*3600.0
      m = fix(eta/60.0)
      eta = eta - m*60.0
      s = fix(eta)
      if h LT 10 then h = '0'+string(h) else h = string(h)
      if m LT 10 then m = '0'+string(m) else m = string(m)
      if s LT 10 then s = '0'+string(s) else s = string(s)
      eta = ' | Time='+strcompress(h+':'+m+':'+s,/remove_all)
      ; Get number of spaces needed to push time string out to
      ; where the end of the final string will be so it will
      ; not appear to move.
      if keyword_set(ncharacters) then $
         nspace = ncharacters - strlen(progstring) $
      else nspace=0
      if strlen(progstring)+nspace+strlen(eta) GE ttysize then begin
         ; We will eventually go over the end of the line so 
         ; reduce number of spaces.
         nspace = (ttysize-strlen(progstring)-strlen(eta))
      endif
      if nspace LT 0 then begin
         ; We will go over the end of the line with this call,
         ; so we need to erase the last time and move to a new 
         ; line.
         if NOT keyword_set(noprint) then begin
            ; erase the last time and write
            ; enough to get to the next line
            nerase = strlen(eta)+nspace
            if nerase GT 0 then $
              print,string(replicate(32b,nerase)),format='(a,$)'
            progstring = ''  ; starting a new line
         endif
      endif
      if nspace GT 0 then spacestr = string(replicate(32b,nspace)) $
      else spacestr = ''
      ; 8b is a backspace so the eta will be overwritten on next call
      eta = spacestr + eta + $
            string(replicate(8b,strlen(eta)+strlen(spacestr)))
      if NOT keyword_set(noprint) then print,eta,format='(a,$)'
      lasttime = thistime  ; the last time the ETA was printed
   endif
endif

if keyword_set(last) then begin
   if NOT keyword_set(noprint) then print
   lastpercent = minval-1.0  ; force a reset on the next call
endif

pstring = progstring
lastpercent = percent

end


PRO LIBRARY
END

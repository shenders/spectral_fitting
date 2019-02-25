Pro findshots, year=year, overwrite=overwrite,$
			  ipval=ipval,$
			  p_range=p_range,$
			  d_range=d_range,$
			  seeding=seeding
			  

; Set shot range
	shotnums   = [27567 , 29157 , 30198 , 31833 , 32820 , 33804 , 35056 , 35266]
	shotyears  = [2012  , 2013  , 2014  , 2015  , 2016  , 2017  , 2018  , 2019]
	if ~keyword_set(year) then year_choice= 2012 else year_choice = year
	id = where(shotyears eq year_choice)
	shot_range = [shotnums[id[0]],shotnums[id[0]+1]]
; Build database
	file = string(year_choice,format='("save/database/shotlog",i4,".idl")')
	if ~file_test(file)or keyword_set(overwrite) then begin
		ip_data   = -1
		ip_shot   = -1
		pinj_data = -1
		pinj_time = -1
		pinj_shot = -1
		tdiv_data = -1
		tdiv_time = -1
		tdiv_shot = -1
		dinj_data = -1
		dinj_time = -1
		dinj_shot = -1
		ninj_data = -1
		ninj_time = -1
		ninj_shot = -1
		edgn_data = -1
		edgn_time = -1
		edgn_shot = -1		
		for i=shot_range[0],shot_range[1] do begin
		; Log plasma current
			read_signal_mrm,0L,i,'MAG','Ipa',time,ip,1
			id = where(time gt 2.0 and time le 8)
			id2 = where(ip[id] ge 300E3)
		; If shot lasts longer than 2s, then log traces
			if id[0] ne -1 and id2[0] ne -1 then begin
				ip = ip[id]
				time=time[id]
				id = where(ip gt 0.5*max(ip))
				shot_time = time[id]
			; Set reduced time array to save disk space
				resol   = 0.1
				tr      = [min(shot_time),max(shot_time)<8.0]
				npts    = floor((tr[1]-tr[0])/resol)				
				if npts eq 0 then goto,skipshot
				ip_data = [ip_data,mean(ip[id])]					
				ip_shot = [ip_shot,i]
				timearr = findgen(npts)*(tr[1]-tr[0])/(npts-1)+tr[0]
			; Log total power		
				;read_signal_mrm,0L,i,'TOT','P_TOT',time1,ptot,2
				;pinj_data = [pinj_data,interpol(ptot,time1,timearr)]
				;pinj_time = [pinj_time,timearr]
				;pinj_shot = [pinj_shot,fltarr(npts)+i]
			; Log Tdiv		
				read_signal_mrm,0L,i,'MAC','Tdiv',time11,tdiv,2
				tdiv_data = [tdiv_data,interpol(tdiv,time11,timearr)]
				tdiv_time = [tdiv_time,timearr]
				tdiv_shot = [tdiv_shot,fltarr(npts)+i]
			; Log fuelling rate
				;read_signal_mrm,0L,i,'UVS','D_tot',time2,dtot,2
				;dinj_data = [dinj_data,interpol(dtot,time2,timearr)]
				;dinj_time = [dinj_time,timearr]
				;dinj_shot = [dinj_shot,fltarr(npts)+i]
			; Log N seeding rate
				read_signal_mrm,0L,i,'UVS','N_tot',time3,ntot,2
				ninj_data = [ninj_data,interpol(ntot,time3,timearr)]
				ninj_time = [ninj_time,timearr]
				ninj_shot = [ninj_shot,fltarr(npts)+i]				
			; Log edge density 
				;read_signal_mrm,0L,i,'DCN','H-5',time4,edens,2
				;edgn_data = [edgn_data,interpol(edens,time4,timearr)]
				;edgn_time = [edgn_time,timearr]
				;edgn_shot = [edgn_shot,fltarr(npts)+i]
				skipshot:
			endif
		endfor
		ip_data   = ip_data[1:*]
		ip_shot   = ip_shot[1:*]
		;pinj_data = pinj_data[1:*] 
		;pinj_time = pinj_time[1:*] 
		;pinj_shot = pinj_shot[1:*]
		tdiv_data = tdiv_data[1:*] 
		tdiv_time = tdiv_time[1:*] 
		tdiv_shot = tdiv_shot[1:*]
		;dinj_data = dinj_data[1:*] 
		;dinj_time = dinj_time[1:*] 
		;dinj_shot = dinj_shot[1:*]
		ninj_data = ninj_data[1:*] 
		ninj_time = ninj_time[1:*] 
		ninj_shot = ninj_shot[1:*]
		;edgn_data = edgn_data[1:*] 
		;edgn_time = edgn_time[1:*] 
		;edgn_shot = edgn_shot[1:*]
		save,file=file,ip_data,ip_shot,$
			       pinj_data,pinj_time,pinj_shot,$
			       dinj_data,dinj_time,dinj_shot,$
			       ninj_data,ninj_time,ninj_shot,$
			       edgn_data,edgn_time,edgn_shot,$
			       tdiv_data,tdiv_time,tdiv_shot			       
	endif else begin
		print,'Restoring database'
		restore,file
	endelse
; Filter Ip range
	if keyword_set(ipval)then begin
		id = where(ip_data ge ipval*1e6-0.12e6 and ip_data lt ipval*1e6+0.05e6)
		if id[0] eq -1 then stop,'No shots in this range'
		ip_data = ip_data[id]
		ip_shot = ip_shot[id]
		idkeep_pinj = -1
		idkeep_dinj = -1
		idkeep_ninj = -1
		idkeep_edgn = -1
		for i=0,n_elements(ip_shot)-1 do begin
			id = where(pinj_shot eq ip_shot[i])
			if id[0] ne -1 then idkeep_pinj = [idkeep_pinj,id]
			id = where(dinj_shot eq ip_shot[i])
			if id[0] ne -1 then idkeep_dinj = [idkeep_dinj,id]
			id = where(ninj_shot eq ip_shot[i])
			if id[0] ne -1 then idkeep_ninj = [idkeep_ninj,id]
			id = where(edgn_shot eq ip_shot[i])
			if id[0] ne -1 then idkeep_edgn = [idkeep_edgn,id]
		endfor
		if n_elements(idkeep_pinj) eq 1 then stop,'No shots in this range'
		idkeep_pinj = idkeep_pinj[1:*]
		idkeep_dinj = idkeep_dinj[1:*]
		idkeep_ninj = idkeep_ninj[1:*]
		idkeep_edgn = idkeep_edgn[1:*]
		pinj_data   = pinj_data[idkeep_pinj]
		pinj_time   = pinj_time[idkeep_pinj]
		pinj_shot   = pinj_shot[idkeep_pinj]
		dinj_data   = dinj_data[idkeep_dinj]
		dinj_time   = dinj_time[idkeep_dinj]
		dinj_shot   = dinj_shot[idkeep_dinj]
		ninj_data   = ninj_data[idkeep_ninj]
		ninj_time   = ninj_time[idkeep_ninj]
		ninj_shot   = ninj_shot[idkeep_ninj]
	endif

; Filter power range
	if keyword_set(p_range)then begin
		id = where(pinj_data ge p_range[0]*1e6 and pinj_data le p_range[1]*1e6)
		if id[0] eq -1 then stop,'No shots in this range'
		pinj_data = pinj_data[id]
		pinj_time = pinj_time[id]
		pinj_shot = pinj_shot[id]
		; Check more than one point exists per shot
		idkeep_pinj = -1
		shotpick    = 0
		for i=0,n_elements(pinj_shot)-1 do begin
			if pinj_shot[i] ne shotpick then begin
				shotpick = pinj_shot[i]
				id = where(pinj_shot eq pinj_shot[i])
				if n_elements(id) gt 1 then idkeep_pinj = [idkeep_pinj,id]
			endif
		endfor
		if n_elements(idkeep_pinj) eq 1 then stop,'No shots in this range'
		idkeep_pinj = idkeep_pinj[1:*]
		pinj_data   = pinj_data[idkeep_pinj]
		pinj_time   = pinj_time[idkeep_pinj]
		pinj_shot   = pinj_shot[idkeep_pinj]			
		idkeep_ip   = -1
		idkeep_dinj = -1
		idkeep_ninj = -1
		idkeep_edgn = -1
		shotpick    = 0
		for i=0,n_elements(pinj_shot)-1 do begin
			if pinj_shot[i] ne shotpick then begin
				shotpick = pinj_shot[i]
				id = where(ip_shot eq pinj_shot[i])
				if id[0] ne -1 then idkeep_ip   = [idkeep_ip,id]
				id = where(dinj_shot eq pinj_shot[i])
				if id[0] ne -1 then idkeep_dinj = [idkeep_dinj,id]
				id = where(ninj_shot eq pinj_shot[i])
				if id[0] ne -1 then idkeep_ninj = [idkeep_ninj,id]
				id = where(edgn_shot eq pinj_shot[i])
				if id[0] ne -1 then idkeep_edgn = [idkeep_edgn,id]
			endif
		endfor
		if n_elements(idkeep_ip) eq 1 then stop,'No shots in this range'
		idkeep_ip   = idkeep_ip[1:*]
		idkeep_dinj = idkeep_dinj[1:*]
		idkeep_ninj = idkeep_ninj[1:*]
		idkeep_edgn = idkeep_edgn[1:*]
		ip_data     = ip_data[idkeep_ip]
		ip_shot     = ip_shot[idkeep_ip]
		dinj_data   = dinj_data[idkeep_dinj]
		dinj_time   = dinj_time[idkeep_dinj]
		dinj_shot   = dinj_shot[idkeep_dinj]
		ninj_data   = ninj_data[idkeep_ninj]
		ninj_time   = ninj_time[idkeep_ninj]
		ninj_shot   = ninj_shot[idkeep_ninj]
	endif

; Filter fuelling range
	if keyword_set(d_range)then begin
		id = where(dinj_data ge d_range[0] and dinj_data le d_range[1])
		if id[0] eq -1 then stop,'No shots in this range'
		dinj_data = dinj_data[id]
		dinj_time = dinj_time[id]
		dinj_shot = dinj_shot[id]
		idkeep_ip   = -1
		idkeep_pinj = -1
		idkeep_ninj = -1
		idkeep_edgn = -1
		shotpick    = 0
		for i=0,n_elements(dinj_shot)-1 do begin
			if dinj_shot[i] ne shotpick then begin
				shotpick = dinj_shot[i]
				id = where(ip_shot eq dinj_shot[i])
				if id[0] ne -1 then idkeep_ip   = [idkeep_ip,id]
				id = where(pinj_shot eq dinj_shot[i])
				if id[0] ne -1 then idkeep_pinj = [idkeep_pinj,id]
				id = where(ninj_shot eq dinj_shot[i])
				if id[0] ne -1 then idkeep_ninj = [idkeep_ninj,id]
				id = where(edgn_shot eq dinj_shot[i])
				if id[0] ne -1 then idkeep_edgn = [idkeep_edgn,id]
			endif
		endfor
		if n_elements(idkeep_ip) eq 1 then stop,'No shots in this range'
		idkeep_ip   = idkeep_ip[1:*]
		idkeep_pinj = idkeep_pinj[1:*]
		idkeep_ninj = idkeep_ninj[1:*]
		idkeep_edgn = idkeep_edgn[1:*]
		ip_data     = ip_data[idkeep_ip]
		ip_shot     = ip_shot[idkeep_ip]
		pinj_data   = pinj_data[idkeep_pinj]
		pinj_time   = pinj_time[idkeep_pinj]
		pinj_shot   = pinj_shot[idkeep_pinj]
		ninj_data   = ninj_data[idkeep_ninj]
		ninj_time   = ninj_time[idkeep_ninj]
		ninj_shot   = ninj_shot[idkeep_ninj]
	endif

; Filter seeding range
	if keyword_set(seeding)then begin
		id = where(ninj_data ge 1e21)
		if id[0] eq -1 then stop,'No shots in this range'
		ninj_data = ninj_data[id]
		ninj_time = ninj_time[id]
		ninj_shot = ninj_shot[id]
		idkeep_ip   = -1
		idkeep_pinj = -1
		idkeep_dinj = -1
		idkeep_tdiv = -1
		idkeep_edgn = -1
		shotpick    = 0
		for i=0,n_elements(ninj_shot)-1 do begin
			if ninj_shot[i] ne shotpick then begin
				shotpick = ninj_shot[i]
				id = where(ip_shot eq ninj_shot[i])
				if id[0] ne -1 then idkeep_ip   = [idkeep_ip,id]
				id = where(pinj_shot eq ninj_shot[i])
				if id[0] ne -1 then idkeep_pinj = [idkeep_pinj,id]
				id = where(dinj_shot eq ninj_shot[i])
				if id[0] ne -1 then idkeep_dinj = [idkeep_dinj,id]
				id = where(edgn_shot eq ninj_shot[i])
				if id[0] ne -1 then idkeep_edgn = [idkeep_edgn,id]
				id = where(tdiv_shot eq ninj_shot[i])
				if id[0] ne -1 then idkeep_tdiv = [idkeep_tdiv,id]
			endif
		endfor
		if n_elements(idkeep_ip) eq 1 then stop,'No shots in this range'
		idkeep_ip   = idkeep_ip[1:*]
		;idkeep_pinj = idkeep_pinj[1:*]
		;idkeep_dinj = idkeep_dinj[1:*]
		;idkeep_edgn = idkeep_edgn[1:*]
		idkeep_tdiv = idkeep_tdiv[1:*]
		ip_data     = ip_data[idkeep_ip]
		ip_shot     = ip_shot[idkeep_ip]
		;pinj_data   = pinj_data[idkeep_pinj]
		;pinj_time   = pinj_time[idkeep_pinj]
		;pinj_shot   = pinj_shot[idkeep_pinj]
		;dinj_data   = dinj_data[idkeep_dinj]
		;dinj_time   = dinj_time[idkeep_dinj]
		;dinj_shot   = dinj_shot[idkeep_dinj]
		tdiv_data   = tdiv_data[idkeep_tdiv]
		tdiv_time   = tdiv_time[idkeep_tdiv]
		tdiv_shot   = tdiv_shot[idkeep_tdiv]
	endif


	setgraphics,xs=600,ys=800,ncol=2,nrow=1,colors=colors,colpick=colpick,collabel=collabel

	Print,'Shots available'
	Print,'---------------'
	imul   = 0
	istart = 0
	iiend  = 0
	repeatgraphs:
	istart = istart + 5.0*imul
	iend   = istart + 4
	if iend gt n_elements(ip_shot)-1 then begin
		iiend = 1 
		iend = n_elements(ip_shot)-1
	endif
	for i=istart,iend do print,string(ip_shot[i],ip_data[i],format='("#",i5," Ip (A) = ",e9.2)')," ",collabel[i-istart]
	
	plot,[0,8],[0,max(tdiv_data)<20],/nodata,xs=1,backgr=colors.white,col=colors.black,xtitle='Time [s]',ytitle='T!ldiv!n [eV]'
	for i=istart,iend do begin
		id = where(tdiv_shot eq ip_shot[i])
		oplot,tdiv_time[id],tdiv_data[id],col=colpick[i-istart]
	endfor
	
	

	plot,[0,8],[0,max(ninj_data)],/nodata,xs=1,backgr=colors.white,col=colors.black,xtitle='Time [s]',ytitle='N!l2!n [#/s]'
	for i=istart,iend do begin
		id = where(ninj_shot eq ip_shot[i])
		oplot,ninj_time[id],ninj_data[id],col=colpick[i-istart]
	endfor
	cursor,x,y,/up
	imul = imul + 1
	Print,'---------------'
	
	if iiend ne 1 then goto,repeatgraphs
Stop
End 	

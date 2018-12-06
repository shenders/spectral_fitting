pro read_data,sig,shot,output,time,trange=trange,ascii=ascii,doplot=doplot,res=res,experiment = experiment,edition=edition

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
	_libddww = '/afs/ipp/aug/ads/lib64/@sys/libddww.so'
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
	signalname = 'Ne'
    endif	
     if sig eq 'LSD' then begin
    	diagnostic = 'LSD'
	signalname = 'ne-ua3'
    endif	
   if ~keyword_set(edition)then edition=0L				; Last closed edition

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
    k2=200000L				
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

    output=fltarr(k2)
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

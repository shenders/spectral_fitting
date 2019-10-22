Pro runabel,load=load,shots=shots,tres=tres,twin=twin,time=time,prad=prad,nosave=nosave,noplot=noplot,flag=flag,debug=debug

	if ~keyword_set(shots)then shots = [35157 , 35167 , 35358 , 35381 ]
	if ~keyword_set(tres)then tres  = 0.1
	if ~keyword_set(twin)then twin  = [2.0,6.5]
	tim   = twin[1]-twin[0] 
	narr  = tim/tres
	time  = findgen(narr)*tim/(narr-1)+twin[0]
	prad  = fltarr(narr)
	for i=0,n_elements(shots)-1 do begin
		file = 'save/abel/'+string(shots[i],format='(i5)')+'/abel.sav'
		if keyword_set(load)then begin
			restore,file
		endif else begin
			cmd = 'mkdir -p save/abel/'+string(shots[i],format='(i5)')
			spawn,cmd
			for j=0,narr-1 do prad[j]=bolom(shots[i],[time[j],time[j]+tres],noplot=noplot,flag=flag)
			if ~keyword_set(nosave)then save,file=file,time,prad
		end
		
		if keyword_set(debug)then begin
			read_signal_mrm,0L,shots[i],'BPD','Prad',tline,pradbpd,2
			read_signal_mrm,0L,shots[i],'BPD','Pradtot',tline,pradtot,2
			setgraphics,colors=colors,ncol=1,nrow=1
			plot,tline,pradbpd,xr=twin,xs=1,back=colors.white,yr=[0,2e7],col=colors.black,title=string(shots[i],format='(i5)')
			oplot,tline,pradtot,col=colors.blue
			if n_elements(time) gt 1 then begin
				oplot,time,prad*1e6,col=colors.red
			endif else begin
				oplot,[time,time+tres],[prad,prad]*1e6,col=colors.red
			end
		endif
	endfor
End

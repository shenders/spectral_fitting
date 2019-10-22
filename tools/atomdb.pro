Pro atomdb,te_arr,dens_arr,tec3995=tec3995,$
		   tec4041=tec4041,$
		   tec4026=tec4026,$
		   tau=tau,rec_only=rec_only,$
		   exc_only=exc_only
		   
	num = 40.0
	div = ceil(n_elements(te_arr)/num)
	tec3995 = -1
	tec4041 = -1
	tec4026 = -1
	for i=0,div-1 do begin
		if i*num+num-1 gt n_elements(te_arr)-1 then begin
			te   = abs(te_arr[i*num:*])>0.5
			dens = abs(dens_arr[i*num:*])
		endif else begin
			te   = abs(te_arr[i*num:i*num+num-1])>0.5
			dens = abs(dens_arr[i*num:i*num+num-1])
		end
		; 399.5 line
		read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=exc3995,block=15
		read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=rec3995,block=65
		; 404.2 line	
		read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=exc4041,block=21
		read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=rec4041,block=71
		; 402.6 line	
		read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=exc4026,block=19
		read_adf15,file='model/atomic_data/pec98#n_ssh_pju#n1.dat',te=te,dens=dens,data=rec4026,block=69
		
		if keyword_set(exc_only)then begin
			rec3995 = 0.0
			rec4041 = 0.0
			rec4026 = 0.0
		endif
		if keyword_set(rec_only)then begin
			exc3995 = 0.0
			exc4041 = 0.0
			exc4026 = 0.0
		endif
		
		if keyword_set(tau)then begin
			run_ionbal,pop=frac,tau=tau,te=te,elem='n',dens=dens,ionbal=ionbal
			tec3995 = [tec3995,(reform(frac[1,*]) * exc3995 + reform(frac[2,*]) * rec3995) ] 
			tec4041 = [tec4041,(reform(frac[1,*]) * exc4041 + reform(frac[2,*]) * rec4041) ]
			tec4026 = [tec4026,(reform(frac[1,*]) * exc4026 + reform(frac[2,*]) * rec4026) ]
		endif else begin
			run_adas405, defyear = '89', uid='adas', year='96', elem='n', te=te, dens=dens, frac=frac
			tec3995 = [tec3995,(frac.ion[*,1] * exc3995 + frac.ion[*,2] * rec3995) ] 
			tec4041 = [tec4041,(frac.ion[*,1] * exc4041 + frac.ion[*,2] * rec4041) ]
			tec4026 = [tec4026,(frac.ion[*,1] * exc4026 + frac.ion[*,2] * rec4026) ]
		end
	endfor
	tec3995 = tec3995[1:*]
	tec4041 = tec4041[1:*]
	tec4026 = tec4026[1:*]
End

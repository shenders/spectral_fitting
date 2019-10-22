Pro runcn,load=load,shots=shots,time=time,cn_mean=cn_mean,cn_err=cn_err,nosave=nosave,noplot=noplot,debug=debug,diag=diag,use_evl=use_evl,channel=channel

	if ~keyword_set(shots)then shots = [35157 , 35167 , 35358 , 35381 ]
	channel = 22
	diag='FVL'
	upper=3.7
	lower=3.5
	for i=0,n_elements(shots)-1 do begin
		file = 'save/cn/'+string(shots[i],format='(i5)')+'/cn.sav'
		if keyword_set(load)then begin
			restore,file
		endif else begin
			cmd = 'mkdir -p save/cn/'+string(shots[i],format='(i5)')
			spawn,cmd
			analysis,shots[i],/interelm,/use_evl,debug=debug,diag=diag,channel=channel,time=time,cn_mean=cn_mean,xr=[2.0,6.0],upper=upper,lower=lower,cn_err=cn_err
			if ~keyword_set(nosave)then save,file=file,time,cn_mean,cn_err
		end
		if keyword_set(debug)then begin
			setgraphics,colors=colors,ncol=1,nrow=1
			plot,time,cn_mean*100,xr=twin,xs=1,back=colors.white,yr=[0,20],col=colors.black,title=string(shots[i],format='(i5)')
			err_plot,time,cn_mean*100,cn_err*100,col=colors.black
		endif
	endfor
End

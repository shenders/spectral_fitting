Pro nii_frac,tdiv=tdiv

; Plot N2 and Ar seeding rates together

if ~keyword_set(tr)then tr=[1.5,6]
if ~keyword_set(tdiv)then tdiv=0

if tdiv eq 8 then begin
	shots    = [ 35157              , 35167            , 35358             , 35381              ]
	times    = [[2.7,-1.0,4.7,-1.0] , [2.7,3.7,4.7,5.5], [3.2,4.5,5.5,-1.0], [2.2,3.0,-1.0,-1.0]]
	diag     = 'EVL'
endif
if tdiv eq 0 then begin
	shots    = [ 35158              , 35165              , 35381		  ]
	times    = [[2.7,3.7,4.7,5.7]   , [2.7,3.7,-1.0,-1.0], [2.7,4.0,-1.0,-1.0]]
	diag     = 'FVL'
endif
ncols = 1
nrows = n_elements(shots)

fraction = -1
err      = -1
ar       = -1
for i=0,n_elements(shots)-1 do begin

	analysis,shots[i],/interelm,/use_evl,channel=22,time=time,cn_mean=cn,xr=[2.0,6.0],upper=3.7,lower=3.5,cn_err=cn_err

	if shots[i] eq 35167 then begin		
		read_signal_mrm,0L,35158,'UVS','CFA03A',tline,ar_rate,2
		ar = ar/1.1e21 * 0.82e21
	endif else begin			
		read_signal_mrm,0L,shots[i],'UVS','CFA03A',tline,ar_rate,2
	endelse

	for j=0,n_elements(times[*,i])-1 do begin
		if times[j,i] eq -1 then goto,skip
		id = where(time ge times[j,i] and time le times[j,i]+0.3)
		if j eq 0 then begin
			if shots[i] eq 35381 and tdiv eq 0 then begin
				analysis,35158,/interelm,/use_evl,channel=22,diag=diag,time=time,cn_mean=cn1,xr=[2.0,6.0],upper=3.7,lower=3.5,cn_err=cn_err1
				initial = mean(cn1[id])
				init_err = mean(cn_err1[id])			
			endif else begin
				initial = mean(cn[id])
				init_err = mean(cn_err[id])
			end
		endif
		if shots[i] eq 35381 then begin
			if j eq 0 then begin
				fraction = [fraction,1.0] 
				err      = [err,0.1]
			endif else begin
				fraction = [fraction,mean(cn[id])/initial]
				err      = [err,mean(cn[id])/initial * sqrt((mean(cn_err[id])/mean(cn[id]))^2 )]			
			end
		endif else begin
			fraction = [fraction,mean(cn[id])/initial]
			err      = [err,mean(cn[id])/initial * sqrt((mean(cn_err[id])/mean(cn[id]))^2 )]
		end
		id = where(tline ge times[j,i] and tline le times[j,i]+0.3)
		if fraction[-1] eq 1.0 then ar=[ar,0] else ar = [ar,mean(ar_rate[id]/1e20)]

		skip:
		
	endfor	
endfor
fraction = fraction[1:*]*100
err      = err[1:*]*100
ar       = ar[1:*]
setgraphics,xs=800,ys=600,ncol=1,nrow=1,psplot=psplot,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list	

user_psym,1,/fill
plot,ar,fraction,psym=8,yr=[0,150],col=colors.black,back=colors.white,$
 xtitle='Ar valve rate [10!u20!n #/s]',ytitle='Divertor c!lN!n fraction [%]'
err_plot,ar,fraction,err,col=colors.black
res=svdfit(ar,fraction,yfit=yfit,measure_errors=err,a=[1,1,1])
xaxis = findgen(100)*20/99.0
yaxis = fltarr(100)
for i=0,n_elements(res)-1 do yaxis = yaxis + res[i] * xaxis^(i)
oplot,xaxis,yaxis,col=colors.red

stop
End

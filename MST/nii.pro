Pro nii,pngplot=pngplot,tr=tr,tdiv=tdiv,ratio=ratio

; Plot N2 and Ar seeding rates together

if ~keyword_set(tr)then tr=[1.5,6]
if ~keyword_set(tdiv)then tdiv=0

if tdiv eq 8 then begin
	shots     = [ 35157   , 35167  , 35358 , 35381 ]
	arsteps   = [[3,4,5] , [3,4,5], [-10,3.5,5.0],[2.5,3.5,4.5]]
endif
if tdiv eq 0 then begin
	shots     = [ 35158   , 35165  , 35356 , 35381 ]
	arsteps   = [[3,4,5] , [3,4,5], [-10,3.5,5.0],[2.5,3.5,4.5]]
endif
ncols = 5
nrows =n_elements(shots)
setgraphics,xs=1400,ys=1000,ncol=ncols,nrow=nrows,psplot=psplot,colors=colors,colpick=colpick,collabel=collabel,full_list=full_list	
for i=0,n_elements(shots)-1 do begin
	read_signal_mrm,0L,shots[i],'UVS','N_tot',t_n2,n2,2
	if shots[i] eq 35167 then begin
		read_signal_mrm,0L,35158,'UVS','CFA03A',t_ar,ar,2
		ar = ar*0.8
	endif else read_signal_mrm,0L,shots[i],'UVS','CFA03A',t_ar,ar,2			
	plot,t_ar,ar/1e21,back=colors.white,col=colors.black,$
		ytitle='[1e21,1e22 #/s]',xr=tr,xs=1,yr=[0,3],ys=1,xtitle='Time [s]',$
		title=string(shots[i],format='(i5,": ")')+'UVS: N_tot/CFA03A'
	oplot,t_n2,n2/1e22,col=colors.blue
	oplot,[arsteps[0,i],arsteps[0,i],arsteps[1,i],arsteps[1,i],arsteps[2,i],arsteps[2,i]],[-50,1e22,1e22,-50,-50,1e22],thick=1.0,lines=5,col=colors.black
endfor

for i=0,n_elements(shots)-1 do begin

	channel=0
	if shots[i] eq 35382 then begin
		read_signal_mrm,0L,shots[i],'GVL','Ar3_4430',t_ar1,ar1,2
		ar = ar1
		t_ar = t_ar1	
		channel = 5	
	endif else begin
		read_signal_mrm,0L,shots[i],'GVL','Ar0_7504',t_ar,ar,2
	end
	plot,t_ar,ar[*,channel]/1e17,back=colors.white,col=colors.black,$
		ytitle='[1e17 ph/s/m!u2!n/sr]',xr=tr,xs=1,yr=[0,0.6],ys=1,xtitle='Time [s]',$
		title=string(shots[i],format='(i5,": ")')+'GVL: Ar0_7504',/nodata

	if shots[i] eq 35356 or shots[i] eq 35358 then begin
		print,'no gvs'
	endif else oplot,t_ar,ar[*,0]/1e17,col=colors.blue
	
	if shots[i] eq 35167 then begin
		read_signal_mrm,0L,35158,'UVS','CFA03A',t_ar,ar,2
		ar = ar*0.8
	endif else read_signal_mrm,0L,shots[i],'UVS','CFA03A',t_ar,ar,2			
	oplot,t_ar,ar/1e22,col=colors.black,thick=0.5
	oplot,[arsteps[0,i],arsteps[0,i],arsteps[1,i],arsteps[1,i],arsteps[2,i],arsteps[2,i]],[-50,1e22,1e22,-50,-50,1e22],thick=1.0,lines=5,col=colors.black
endfor
if keyword_set(ratio)then begin
for i=0,n_elements(shots)-1 do begin
	read_signal_mrm,0L,shots[i],'FVL','N_1_3995',t_nii,nii399,2
	read_signal_mrm,0L,shots[i],'FVL','N_1_4041',t_nii,nii404,2
	chn=22 & nii = nii404[*,chn]/nii399[*,chn]


	plot,t_nii,nii,back=colors.white,col=colors.black,$
		ytitle='[ph/s/m2/sr]',xr=tr,xs=1,yr=[0,0.4],ys=1,xtitle='Time [s]',$
		title=string(shots[i],format='(i5,": ")')+'FVL: N_1_3995',/nodata


	read_signal_mrm,0L,shots[i],'EVL','N_1_3995',t_nii2,nii399,2
	read_signal_mrm,0L,shots[i],'EVL','N_1_4041',t_nii2,nii404,2
	chn=22 & nii2 = nii404[*,chn]/nii399[*,chn]

	read_signal_mrm,0L,shots[i],'EVL','N_1_3995',t_nii3,nii399,2
	read_signal_mrm,0L,shots[i],'EVL','N_1_4041',t_nii3,nii404,2
	chn=16 & nii3 = nii404[*,chn]/nii399[*,chn]

	read_signal_mrm,0L,shots[i],'EVL','N_1_3995',t_nii4,nii399,2
	read_signal_mrm,0L,shots[i],'EVL','N_1_4041',t_nii4,nii404,2
	chn=11 & nii4 = nii404[*,chn]/nii399[*,chn]

	telm = find_elm(shots[i],t_nii)
	elmcond=4.5
	idelm     = where(telm ge elmcond)
	if idelm[0] eq -1 then idelm = findgen(n_elements(time))
	nii  = nii[idelm]  & t_nii  = t_nii[idelm]
	nii2 = nii2[idelm] & t_nii2 = t_nii2[idelm]
	nii3 = nii3[idelm] & t_nii3 = t_nii3[idelm]
	nii4 = nii4[idelm] & t_nii4 = t_nii4[idelm]
	user_psym,1,/fill
	oplot,t_nii,smooth(nii,10),col=colors.red
	user_psym,3,/fill
	oplot,t_nii2,smooth(nii2,10),col=colors.blue
	oplot,t_nii3,smooth(nii3,10),col=colors.green
	oplot,t_nii4,smooth(nii4,10),col=colors.orange
	
	oplot,[arsteps[0,i],arsteps[0,i],arsteps[1,i],arsteps[1,i],arsteps[2,i],arsteps[2,i]],[-50,1e22,1e22,-50,-50,1e22],thick=1.0,lines=5,col=colors.black
	oplot,[0,10],[8e18,8e18],linest=5,col=colors.black,thick=1.0
	oplot,[0,10],[6e18,6e18],linest=5,col=colors.black,thick=1.0
	oplot,[0,10],[5e18,5e18],linest=5,col=colors.black,thick=1.0
	oplot,[0,10],[4e18,4e18],linest=5,col=colors.black,thick=1.0
	oplot,[0,10],[3e18,3e18],linest=5,col=colors.black,thick=1.0
	oplot,[0,10],[2e18,2e18],linest=5,col=colors.black,thick=1.0
endfor



endif else begin
for i=0,n_elements(shots)-1 do begin
	if shots[i] eq 33257 then begin
		restore,'save/33257/ROV013-data.idl'
		nii = output.nii[*,0,0] * 3
		t_nii=output.time
	endif else begin
		read_signal_mrm,0L,shots[i],'FVL','N_1_3995',t_nii,nii,2
		nii = nii[*,22]

	end
	plot,t_nii,nii,back=colors.white,col=colors.black,$
		ytitle='[ph/s/m2/sr]',xr=tr,xs=1,yr=[0,1e19],ys=1,xtitle='Time [s]',$
		title=string(shots[i],format='(i5,": ")')+'FVL: N_1_3995',/nodata
	read_signal_mrm,0L,shots[i],'EVL','N_1_3995',t_nii2,nii2,2
	nii2 = nii2[*,22]

	read_signal_mrm,0L,shots[i],'EVL','N_1_3995',t_nii3,nii3,2
	nii3 = nii3[*,16]

	read_signal_mrm,0L,shots[i],'EVL','N_1_3995',t_nii4,nii4,2
	nii4 = nii4[*,11]

	telm = find_elm(shots[i],t_nii)
	elmcond=4.5
	idelm     = where(telm ge elmcond)
	if idelm[0] eq -1 then idelm = findgen(n_elements(time))
	nii  = nii[idelm]  & t_nii  = t_nii[idelm]
	nii2 = nii2[idelm] & t_nii2 = t_nii2[idelm]
	nii3 = nii3[idelm] & t_nii3 = t_nii3[idelm]
	nii4 = nii4[idelm] & t_nii4 = t_nii4[idelm]
	user_psym,1,/fill
	oplot,t_nii,smooth(nii,10),col=colors.red
	user_psym,3,/fill
	oplot,t_nii2,smooth(nii2,10),col=colors.blue
	oplot,t_nii3,smooth(nii3,10),col=colors.green
	oplot,t_nii4,smooth(nii4,10),col=colors.orange
	
	oplot,[arsteps[0,i],arsteps[0,i],arsteps[1,i],arsteps[1,i],arsteps[2,i],arsteps[2,i]],[-50,1e22,1e22,-50,-50,1e22],thick=1.0,lines=5,col=colors.black
	oplot,[0,10],[8e18,8e18],linest=5,col=colors.black,thick=1.0
	oplot,[0,10],[6e18,6e18],linest=5,col=colors.black,thick=1.0
	oplot,[0,10],[5e18,5e18],linest=5,col=colors.black,thick=1.0
	oplot,[0,10],[4e18,4e18],linest=5,col=colors.black,thick=1.0
	oplot,[0,10],[3e18,3e18],linest=5,col=colors.black,thick=1.0
	oplot,[0,10],[2e18,2e18],linest=5,col=colors.black,thick=1.0
endfor
endelse
for i=0,n_elements(shots)-1 do begin
	if shots[i] eq 33257 then begin
		restore,'save/33257/ROV013-data.idl'
		nii = output.wi[*,0,0] * 3
		t_nii=output.time
	endif else begin
		read_signal_mrm,0L,shots[i],'FVL','W_0_4009',t_nii,nii,2
		nii = nii[*,22]

	end
	plot,t_nii,nii,back=colors.white,col=colors.black,$
		ytitle='[ph/s/m2/sr]',xr=tr,xs=1,yr=[0,2e17],ys=1,xtitle='Time [s]',$
		title=string(shots[i],format='(i5,": ")')+'EVL: W_0_4009',/nodata
	;oplot,t_nii,nii,col=colors.red
	if shots[i] eq 33257 then begin
		restore,'save/33257/ROV013-data.idl'
		nii2 = output.wi[*,0,0] * 3.0
		t_nii2=output.time
	endif else begin
		read_signal_mrm,0L,shots[i],'EVL','W_0_4009',t_nii2,nii2,2
		nii2 = nii2[*,21]
	end
	;oplot,t_nii2,nii2,col=colors.blue

	user_psym,3,/fill
	oplot,t_nii2,smooth(nii2,10),col=colors.blue
	oplot,[arsteps[0,i],arsteps[0,i],arsteps[1,i],arsteps[1,i],arsteps[2,i],arsteps[2,i]],[-50,1e22,1e22,-50,-50,1e22],thick=1.0,lines=5,col=colors.black
endfor

for i=0,n_elements(shots)-1 do begin
	read_signal_mrm,0L,shots[i],'MAC','Tdiv',t_div,tdiv,2
	read_signal_mrm,0L,shots[i],'INJ','TdivSoll',t_inj,inj,2
	plot,t_div,tdiv,back=colors.white,col=colors.black,$
		ytitle='[eV]',xr=tr,xs=1,yr=[-5,30],ys=1,xtitle='Time [s]',$
		title=string(shots[i],format='(i5,": ")')+'MAC: Tdiv',/nodata
	oplot,t_div,tdiv,col=colors.red
	oplot,t_inj,inj,col=colors.black
	oplot,[arsteps[0,i],arsteps[0,i],arsteps[1,i],arsteps[1,i],arsteps[2,i],arsteps[2,i]],[-50,1e22,1e22,-50,-50,1e22],thick=1.0,lines=5,col=colors.black
endfor

;for i=0,n_elements(shots)-1 do begin
;	read_signal_mrm,0L,shots[i],'ELM','f_ELM',t_elm,elm,2
;	plot,t_elm,elm,back=colors.white,col=colors.black,$
;		ytitle='[Hz]',xr=tr,xs=1,yr=[0,250],ys=1,xtitle='Time [s]',$
;		title=string(shots[i],format='(i5,": ")')+'ELM: f_ELM'
;	oplot,[arsteps[0,i],arsteps[0,i],arsteps[1,i],arsteps[1,i],arsteps[2,i],arsteps[2,i]],[-50,1e22,1e22,-50,-50,1e22],thick=1.0,lines=5,col=colors.black
;endfor

if keyword_set(pngplot)then write_png,'out.png',tvrd(/true)


stop
end

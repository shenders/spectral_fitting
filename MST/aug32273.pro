pro plot_profs,output,idelm,colors=colors,sig=sig,single=single	
	user_psym,5,/fill
	if keyword_set(single)then begin
		plot,output.time[idelm],output.nii[idelm,0,0]/1e18,back=colors.white,col=colors.black,title=sig,psym=8,xtitle='Time [s]',ytitle='N II [10!u18!n ph/s/m!u2!n/sr]'	
		oplot,output.time[idelm],output.niii[idelm,0,1]/1e17,col=colors.blue,psym=8	
		oplot,output.time[idelm],output.niv[idelm]/1e16,col=colors.red,psym=8		
	endif else begin
		plot,output.time,output.nii[*,0,0]/1e18,back=colors.white,col=colors.black,title=sig,psym=8,xtitle='Time [s]',ytitle='N II [10!u18!n ph/s/m!u2!n/sr]'	
		oplot,output.time[idelm],output.nii[idelm,0,0]/1e18,col=colors.red,psym=8	
		plot,output.time,output.niii[*,0,1]/1e17,back=colors.white,col=colors.black,title=sig,psym=8,xtitle='Time [s]',ytitle='N III [10!u17!n ph/s/m!u2!n/sr]'	
		oplot,output.time[idelm],output.niii[idelm,0,1]/1e17,col=colors.red,psym=8	
		plot,output.time,output.niv/1e16,back=colors.white,col=colors.black,psym=8,title=sig,xtitle='Time [s]',ytitle='N IV [10!u16!n ph/s/m!u2!n/sr]'	
		oplot,output.time[idelm],output.niv[idelm]/1e16,col=colors.red,psym=8	
	end
end
Pro aug32273,single=single,psplot=psplot

shot=32273 & shotstr='32273'
restore,'save/'+shotstr+'/ROV014-data.idl'
print,'Using interELM ...'
telm = find_elm(shot,output.time)
idelm = where(telm ge 4.5)
if keyword_set(single)then setgraphics,ncol=3,nrow=1,xs=1000,ys=500,psplot=psplot,colors=colors else setgraphics,nrow=3,psplot=psplot,file='N_intens.ps',ncol=3,/landscape,xs=1000,ys=800,colors=colors
plot_profs,output,idelm,colors=colors,sig='ROV014',single=single
restore,'save/'+shotstr+'/ROV012-data.idl' & plot_profs,output,idelm,colors=colors,sig='ROV012',single=single
restore,'save/'+shotstr+'/ROV010-data.idl' & plot_profs,output,idelm,colors=colors,sig='ROV010',single=single
if keyword_set(psplot)then setgraphics,/close
Stop

END

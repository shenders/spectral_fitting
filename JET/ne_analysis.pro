Pro ne_analysis,psplot=psplot

; Restore save files containing time averaged spectrum from KT3A SP-8:
; w_ne.sav - 87522 t=45-55s
; wo_ne.sav - 87182 t=45-55s
;
restore,'JET/w_ne.sav'
w1=wavelength-0.3
s1=spec_avr
restore,'JET/wo_ne.sav'
w2=wavelength-0.3
s2=spec_avr
setgraphics,nrow=2,ncol=1,xs=800,ys=1200,psplot=psplot,filename=filename,colors=colors,/portrait

plot,w1,s1/10/1e19,xr=[360,385],col=colors.black,back=colors.white,$
     yr=[0,1],pos=[0.15,0.5,0.9,0.85]

oplot,w2,s2/10/1e19,col=colors.red

ne_1=[375.12,370.96,377.71,373.49,366.41,376.63,369.42]                                                                                                                    
ne_2=[372.71,364.39,371.31]                                                                                                                                                
ne_3=[371.35,381.36]                                                                                                                                                       

for i=0,n_elements(ne_1)-1 do oplot,[ne_1[i],ne_1[i]],[0,5],col=colors.blue,linest=5                                                                                       
for i=0,n_elements(ne_2)-1 do oplot,[ne_2[i],ne_2[i]],[0,5],col=colors.green,linest=5                                                                                      
for i=0,n_elements(ne_3)-1 do oplot,[ne_3[i],ne_3[i]],[0,5],col=colors.orange,linest=5    

dens_arr = [7e19 , 8e19 , 9e19 , 1e20]
te_arr   = [10.0 , 6.0  , 4.0  , 2.0 ]
c_ne     = fltarr(n_elements(dens_arr))+0.000002
dist     = [0.3  , 0.1  , 0.05 , 0.01]
tec      = 0.0
for i=0,n_elements(te_arr)-1 do begin
    x=adas_tec(te_arr[i],dens_arr[i],'ne',3600.0,3850.0,0.0,1.6)
    tec = tec + x.tec * c_ne[i]*dens_arr[i] * dist[i]
endfor

plot,x.wavelength/10,tec/4/!pi/1e19,xr=[360,385],pos=[0.15,0.15,0.9,0.45],xtitle='Wavelength [nm]',xs=9,col=colors.black,back=colors.white    

Stop
End                                                                                 

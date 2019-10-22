Pro gfunct0,x , a, f, pder

f    = a[0] / x
pder = [ [ 1.0/x ]  ]	
end
Pro gfunct1,x , a, f, pder

f    = a[0] * x^(-2.0)
pder = [ [ x^(-2.0) ]  ]	
end
Pro gfunct2,x , a, f, pder

f    = a[0] * x^(-3.0)
pder = [ [ x^(-3.0)]  ]	
end
Pro gfunct3,x , a, f, pder

f    = a[0] * x^(-2.5)
pder = [ [ x^(-2.5) ] ]	
end

Pro pfunct1,x , a, f, pder

f    = a[0] * x
pder = [ [ x ]  ]	
end

Pro pfunct2,x , a, f, pder

f    = a[0] * x^(0.8)
pder = [ [ x^(0.8) ]  ]	
end
Pro pfunct3,x , a, f, pder

f    = a[0] * x^(0.6)
pder = [ [ x^(0.6)]  ]	
end

Pro augped_plot,full=full,psplot=psplot,extended=extended,shots=shots,$
  trange=trange,skip=skip,debug=debug,define_dens=define_dens,use_flux=use_flux
  
if keyword_set(extended)then begin
	shots= [ 36655    , 30506      ,30623	   , 32273	]
	trange=[ [3.4,3.6], [2.65,2.75],[4.97,5.04], [3.84,4.05]] 
endif else begin
	shots= [32244	   , 34368   , 30554	, 30776    ,  34358    ,30307	   ,30306   ,30298    ,34971	  , 34973    ,35158	, 34359     , 30505    ]
	trange=[[2.45,2.75],[4.0,4.2],[3.9,4.1] , [3.1,3.3], [4.3,5.0] ,[3.55,3.65],[3.5,3.7],[3.8,4.1],[3.35,3.43], [3.2,3.3],[2.1,2.28], [4.0,4.3], [4.3,4.8]] 
end
if keyword_set(skip)then return

if keyword_set(full)then begin
setgraphics,ncol=3,nrow=2,xs=1500,ys=900,colpick=colpick,colors=colors,collabel=collabel ,/full_list
user_psym,1
plot,[0,15],[0,10],/nodata,back=colors.white,col=colors.black,xtitle='T!ldiv!n [eV]',ytitle='H-5 n!le!n [10!u19!n m!u-3!n]'

for i=0,n_elements(shots)-1 do begin
	read_signal_mrm,0L,shots[i],'DCN','H-5',t,h5,1
	id = where(t ge trange[0,i] and t le trange[1,i])
	read_signal_mrm,0L,shots[i],'MAC','Tdiv',t1,td1,2
	idd = where(t1 ge 1.0 and t1 le 7.0)
	td2 = smooth(td1,40)
	tdiv = interpol(td2[idd],t1[idd],t[id])
	oplot,tdiv,h5[id]/1e19,psym=8,col=colpick[i]
endfor
for i=0,n_elements(shots)-1 do print,collabel[i],shots[i]
endif 

;------------------------------------------------------------------
; Get ne,sep and Te,sep from saved AUGPED files
;------------------------------------------------------------------

te=-1 & dens=-1 & te_err=-1 & dens_err=-1
for i=0,n_elements(shots)-1 do begin
	augped,shots[i],nesep,tesep,nesep_err,tesep_err
	te=[te,tesep]
	dens=[dens,nesep]
	te_err=[te_err,tesep_err]
	dens_err=[dens_err,nesep_err]
endfor
te=te[1:*]
dens=dens[1:*]
te_err=te_err[1:*]
dens_err=dens_err[1:*]


;------------------------------------------------------------------
; Get Psep from Pierre's inverted bolometry calculations
;------------------------------------------------------------------

pmain=-1 & pmain_err = -1
diag = 'BPT'
exp='davidp'
for i=0,n_elements(shots)-1 do begin
	time=-1
	pcore=-1
	read_signal_mrm,0L,shots[i],diag,'Pr_main',time,Pcore,1,exp=exp
	
	id = where(time ge trange[0,i] and time le trange[1,i])		
	tmp = moment(Pcore[id])
	Pmain = [pmain,tmp[0]]
	Pmain_err = [pmain_err,sqrt(tmp[1])]
endfor

pmain=pmain[1:*]
pmain_err=pmain_err[1:*]

;------------------------------------------------------------------
; Get Pdiv from calibrated DDS nDivIst measurements
;------------------------------------------------------------------

Pdiv=-1.0 & Pdiv_err=-1
for i=0,n_elements(shots)-1 do begin
	aug_div_pressure,shots[i],trange[*,i]
	restore,'save/'+string(shots[i],format='(I5)')+'/div_pressure.sav'
	Pdiv = [Pdiv,press]
	Pdiv_err = [Pdiv_err,error]
endfor
Pdiv = Pdiv[1:*]
Pdiv_err = Pdiv_err[1:*]


;------------------------------------------------------------------
; Collect remaining parameters
;------------------------------------------------------------------

user_psym,1,/fill
Psep=-1.0 & Psep_err=-1.0
bp=-1.0
ngw=-1.0 & ngw_err=-1.0
aminor=-1.0
tdiv=-1.0
strikep = -1
k=-1.0
cn_mean=-1.0
cn_err=-1.0
cn_dens=-1.0
cn_dens_err=-1.0
cn_te=-1.0
cn_te_err=-1.0
cn_flux = -1.0
cn_flux_err = -1.0
flux_ne = -1.0
flux_ne_err = -1.0
flux_e = -1.0
flux_e_err = -1.0
pcur=-1.0
ndiv=-1.0
Pinj=-1.0
ne_core=-1.0
for i=0,n_elements(shots)-1 do begin

;------------------------------------------------------------------
; Psep parameter
;------------------------------------------------------------------

	time1 = -1
	ptot = -1
	read_signal_mrm,0L,shots[i],'TOT','P_TOT',time1,Ptot,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	tmp = moment(ptot[id])
	pin = tmp[0]
	pin_err = sqrt(tmp[1])
	if max(time1) lt trange[1,i] then begin
		time1 = -1
		pnbi =-1
		read_signal_mrm,0L,shots[i],'NIS','PNI',time1,PNBI,1
		id = where(time1 ge trange[0,i] and time1 le trange[1,i])
		tmp = moment(ptot[id])
		pni = tmp[0]
		pni_err = sqrt(tmp[1])
		
		time1=-1
		picrn=-1
		read_signal_mrm,0L,shots[i],'ICP','PICRN',time1,PICRN,1
		id = where(time1 ge trange[0,i] and time1 le trange[1,i])
		tmp = moment(ptot[id])
		pic = tmp[0]
		pic_err = sqrt(tmp[1])
		
		time1 = -1
		pecrh = -1
		read_signal_mrm,0L,shots[i],'ECS','PECRH',time1,PECRH,1
		id = where(time1 ge trange[0,i] and time1 le trange[1,i])
		tmp = moment(ptot[id])
		pec = tmp[0]
		pec_err = sqrt(tmp[1])

		pin = pni + pic + pec
		pin_err = sqrt(pni_err^2 + pic_err^2 + pec_err^2)
	endif
	Pinj = [Pinj,pin]
	Pup  = (Pin - Pmain[i])/2.37
	Pup_err = Pup * sqrt((sqrt(Pin_err^2 + Pmain_err[i]^2)/Pup)^2)
	Psep = [Psep,Pup]
	Psep_err = [Psep_err,Pup_err]

	; Plasma kappa, minor radius, current
	time1 = -1
	kappa = -1
	read_signal_mrm,0L,shots[i],'TOT','kappa',time1,kappa,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	k = [k,mean(kappa[id])]
	kap = mean(kappa[id])
	
	time1 = -1
	amin = -1
	read_signal_mrm,0L,shots[i],'TOT','ahor',time1,amin,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	a = mean(amin[id])
	aminor = [aminor,a]
	
	time1=-1
	curr=-1
	read_signal_mrm,0L,shots[i],'MAG','Ipa',time1,curr,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	ip = mean(curr[id])
	pcur = [pcur,ip]

	time1=-1
	curr=-1
	read_signal_mrm,0L,shots[i],'MAC','Tdiv',time1,tdiv1,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	tdiv = [tdiv,mean(tdiv1[id])]

	time1= -1
	ne0  = -1
	read_signal_mrm,0L,shots[i],'TOT','n/nGW',time1,ne0,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	
	GW = ((ip/1e6) / !pi / a^2) 
	ne_core = [ne_core,mean(ne0[id])*GW* 1e20]

; Derivation of Goldston's scaling law
;-------------------------------------
;	cz = Psep / Bp / (1+k^2)^2 / fGW
;	   = Psep / (Ip/a/sqrt(0.5*(1+k^2)^2)) / (1+k^2)^2 / (ne,sep/(Ip/a^2))^2
;	   = Psep * a * sqrt(0.5(1+k^2)^2) Ip^2 * a^4 / Ip / (1+k^2)^2/ ne,sep^2
;	   = Psep * a^5 * (0.5(1+k^2)^2)^0.5 Ip / (1+k^2)^2 / ne,sep^2
;	   = Psep * a^5 * Ip / (1+k^2) / ne,sep^2   
	
	bp = [bp,((4*!pi*0.0000001) * ip)/(2.0*!pi * a * sqrt(0.5*(1+mean(kap)^2)))]
	ngw= [ngw,(dens[i]/1e20) / GW]
	ngw_err = [ngw_err,(dens_err[i]/1e20) / GW]
		
	file = 'save/cn/'+string(shots[i],format='(i5)')+'/cn_database_2.sav'
	restore,file
	cn_mean = [cn_mean,conc]
	cn_err = [cn_err,err]
	cn_dens = [cn_dens,dens_nii]
	cn_dens_err = [cn_dens_err,dens_nii_err]
	cn_te = [cn_te,te_nii]
	cn_te_err = [cn_te_err,te_nii_err]

	file = 'save/'+string(shots[i],format='(i5)')+'/flux_cn.sav'
	restore,file
	cn_flux = [cn_flux,flux]
	cn_flux_err = [cn_flux_err,flux_err]
	flux_ne = [flux_ne,denom]
	flux_ne_err = [flux_ne_err,denom_err]
	flux_e = [flux_e,denom2]
	flux_e_err = [flux_e_err,denom2_err]

	time1= -1
	spnt  = -1
	read_signal_mrm,0L,shots[i],'FPG','Suna2b',time1,spnt,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	strikep = [strikep,mean(spnt[id])]

End
Pcur = Pcur[1:*]
Psep = Psep[1:*] & Psep_err=Psep_err[1:*]
Pinj = Pinj[1:*]
Tdiv = Tdiv[1:*]
aminor = aminor[1:*]
cn_flux = cn_flux[1:*]
cn_flux_err = cn_flux_err[1:*]
flux_ne = flux_ne[1:*]
flux_ne_err = flux_ne_err[1:*]
flux_e = flux_e[1:*]
flux_e_err = flux_e_err[1:*]

strikep = strikep[1:*]
ev2k = 1.1604E4
boltz= 1.38064852E-23
ne_core=ne_core[1:*]
bp=bp[1:*]
ngw=ngw[1:*]
ngw_err = ngw_err[1:*]
k=k[1:*]
cn_mean=cn_mean[1:*]
cn_err=cn_err[1:*]
cn_dens=cn_dens[1:*]
cn_dens_err=cn_dens_err[1:*]
cn_te=cn_te[1:*]
cn_te_err=cn_te_err[1:*]
scal = 1.3

if keyword_set(define_dens)then begin
	fac  = dens / (Pdiv^0.31 * 2.65 * 1e19)
	dens = Pdiv^0.31 * 2.65 * 1e19
	ngw  = ngw / fac
endif
if keyword_set(use_flux)then begin
	cn_mean = cn_flux
	cn_err = cn_flux_err
endif
;-------------------------------------------------------
; Plot the nesep vs. Pdiv and ne,sep vs Psep
;-------------------------------------------------------

setgraphics,ncol=1,nrow=2,xs=600,ys=800,colpick=colpick,colors=colors,collabel=collabel ,/full_list,psplot=psplot,file='database_parms.eps'

plot,ngw,Psep/1e6,psym=8,yr=[0,6],xr=[0,0.4],$
col=colors.black,back=colors.white,ytitle='P!lsep!n [MW]',xtitle='f!lGW,sep!n=n!le,sep!n / n!lGW!n'
errors,ngw,Psep/1e6,ystd=Psep_err/1e6,xstd=ngw_err,col=colors.black

xa  = (findgen(100)*(10-0.1)/99.0+0.1)^0.31 * 2.65
xae = xa * 0.31 * (Pdiv_err/Pdiv)  
plot,Pdiv,dens/1e19,back=colors.white,col=colors.black,psym=8,/xlog,/ylog,$
 xtitle='p!ldiv!n [Pa]',ytitle='n!le,sep!n [10!u19!n m!u-3!n]',yr=[1,7],ys=1,xr=[0.1,10]
errors,Pdiv,dens/1e19,ystd=dens_err/1e19,xstd=Pdiv_err,col=colors.black
oplot,(findgen(100)*(10-0.1)/99.0+0.1),xa,col=colors.black,linest=5

;-------------------------------------------------------
; Plot the ne,NII vs. Pdiv, Psep v cn and ngw vs cn
;-------------------------------------------------------

setgraphics,ncol=1,nrow=2,xs=600,ys=800,colpick=colpick,colors=colors,collabel=collabel ,/full_list,psplot=psplot,file='ne_nii_correlation.eps'
plot,Pdiv,cn_dens/1e14,back=colors.white,col=colors.black,psym=8,xtitle='p!ldiv!n [Pa]',ytitle='n!le,N II!n [10!u20!n m!u-3!n]',xr=[0,5],yr=[0,3]
errors,Pdiv,cn_dens/1e14,ystd=cn_dens_err/1e14,xstd=Pdiv_err,col=colors.black
ii = sort(Pdiv)
xx = Pdiv[ii]
yy = cn_dens[ii]/1e14
rsqr = r2(xx,yy,b=b,m=m)
xfit = findgen(100)/99.0 * 20.0 
yfit = b + m * xfit
oplot,xfit,yfit,col=colors.black,linest=5
legend,'R!u2!n='+string(rsqr,format='(d4.2)'),colors.black,yshift=-0.1,xshift=-0.6
legend,string(b,m,format='("y=",d5.2,"+",d4.2,"x")'),colors.black,yshift=0.0,xshift=-0.6

ngw_scan = [3.5,5.5]*1e6/2.37

id = where(Psep gt ngw_scan[0] and psep le ngw_scan[1])
print,'ngw: ',ngw[id]
print,'Bp for ne,sep scan: ',bp[id]
print,'Kappa for ne,sep scan: ',k[id]

plot,ngw[id],cn_mean[id],back=colors.white,col=colors.black,psym=8,$
  yr=[0,15],ytitle='Spectroscopy c!lN !n[%]',xtitle='f!lGW,sep!n=n!le,sep!n / n!lGW!n'
errors,ngw[id],cn_mean[id],ystd=cn_err[id],xstd=ngw_err[id],col=colors.black

xval = ngw[id]
yval = cn_mean[id]
err  = cn_err[id]
xs   = sort(xval)
xval = xval[xs]
yval = yval[xs]
err  = err[xs]
a=[60.0]
xnew = findgen(100)/99.0 * (max(xval)-min(xval)) + min(xval)
x=curvefit(xval,yval,err,a,sig,function_name='gfunct3',itmax=100)
ynew = a[0] * (xnew^(-2.5))
oplot,xnew,ynew,linest=5,col=colors.green

a=[40.0]
xnew = findgen(100)/99.0 * (max(xval)-min(xval)) + min(xval)
x=curvefit(xval,yval,err,a,sig,function_name='gfunct2',itmax=100)
ynew = a[0] * (xnew^(-3.0))
oplot,xnew,ynew,linest=5,col=colors.red


a=[30.0]
xnew = findgen(100)/99.0 * (max(xval)-min(xval)) + min(xval)
x=curvefit(xval,yval,err,a,sig,function_name='gfunct1',itmax=100)
ynew = a[0] * (xnew^(-2.0))
oplot,xnew,ynew,linest=5,col=colors.blue
legend,'1/f!lGW,sep!u2.0',colors.blue,yshift=0.0
legend,'1/f!lGW,sep!u2.5',colors.green,yshift=-0.1
legend,'1/f!lGW,sep!u3.0',colors.red,yshift=-0.2


;-------------------------------------------------------
; Plot the scaling law and Pdiv comparison to spectroscopy 
;-------------------------------------------------------
rep:
strscl = string(scal,format='(d3.1)')
;-------------------------------------------------------
; Calculate Goldston scaling factor
;-------------------------------------------------------

gold     = (4.0 / 18.3) * ((Psep/1e6)^(1.0) / (bp * (1+k^2)^(1.5) * ngw^(2.)))
gold_err = gold * sqrt( (Psep_err/Psep)^2 + ( 2.0 * ngw_err / ngw)^2 )

;-------------------------------------------------------
; Calculate Kallenbach scaling factor
;-------------------------------------------------------

R       = 1.65
eps     = aminor/R
Pup     = Psep/1e6*0.8
lambdaq = 1.35E-3 * Pup^(-0.02) * R^(0.04) * Bp^(-0.92) * eps^(0.42)
lambda_int = 0.005
kall    = ((Pup/R/Pdiv) * 1.3  - 1) * 1/18.0 

;-------------------------------------------------------
; Print out database values
;-------------------------------------------------------
print,'        Shot       cn_flux         Spec cN     Goldston         Psep          Pin     f_GW,sep       ne,sep      ne,core     div pres         tdiv      Kappa   Bp  kappa ip cz_red'
for i=0,n_elements(shots)-1 do print,shots[i],',',cn_flux[i],',',cn_mean[i],',',gold[i],',',psep[i]/1e6,',',Pinj[i]/1e6,',',ngw[i],',',dens[i]/1e19,',',ne_core[i]/1e19,',',Pdiv[i],',',tdiv[i],',',k[i],',',bp[i],',',k[i],',',pcur[i]/1e6,cn_mean[i]*aminor[i]^3

;-------------------------------------------------------
; Compare with flux ratio
;-------------------------------------------------------

setgraphics,ncol=1,nrow=2,xs=600,ys=800,colpick=colpick,colors=colors,collabel=collabel ,/full_list,psplot=psplot,file='flux_ratio.eps'

plot,cn_flux,cn_mean,psym=8,yr=[0,20],xr=[0,20],$
 col=colors.black,back=colors.white,ytitle='Spectroscopy c!lN !n[%]',xtitle='Flux ratio c!lN !n[%]'
errors,cn_flux,cn_mean,xstd=cn_flux_err,ystd=cn_err,col=colors.black
oplot,[0,20],[0,20],linest=5,col=colors.black

flux_comparison,val,diff
flux_comparison,val_ext,diff_ext,/extra
user_psym,1,/fill & plot,val/1e17,diff,col=colors.black,back=colors.white,psym=8,symsize=2.0,$
xtitle='N II @399.5nm [10!u17!n ph/s/m!u2!n/sr]',ytitle='c!lN!n spec/flux',xr=[0,20]
rsqr = r2(val/1e17,diff,b=b,m=m)
xfit = findgen(100)/99.0 * 20.0 
yfit = b + m * xfit
oplot,xfit,yfit,col=colors.black,linest=5
legend,'R!u2!n='+string(rsqr,format='(d4.2)'),colors.black,yshift=-0.1,xshift=-0.6
legend,string(b,m,format='("y=",d5.2,"+",d4.2,"x")'),colors.black,yshift=0.0,xshift=-0.6
legend,'N II before seeding',colors.black,yshift=0.1,xshift=-0.6
oplot,val_ext/1e17,diff_ext,psym=8,col=colors.red,symsize=2.0

;-------------------------------------------------------
; Compare with Goldston and regression
;-------------------------------------------------------

setgraphics,ncol=1,nrow=2,xs=600,ys=800,colpick=colpick,colors=colors,collabel=collabel ,/full_list,psplot=psplot,file='cn_scaling.eps'

plot,gold,cn_mean,back=colors.white,col=colors.black,psym=8,yr=[0,max(cn_mean+cn_err)>max(gold)],xr=[0,max(cn_mean+cn_err)>max(gold)],$
    xtitle='0.219 P!lsep!n / <B!lp!n> (1+k!u2!n)!u3/2!n f!lGW,sep!u2!n',ytitle='Spectroscopy c!lN!n [%]'
 
errors,gold,cn_mean,ystd=cn_err,xstd=gold_err,col=colors.black

rsqr = r2(gold,cn_mean,b=b,m=m)
xfit = findgen(100)/99.0 * 20.0 
yfit = b + m * xfit
oplot,xfit,yfit,col=colors.black,linest=5
legend,'R!u2!n='+string(rsqr,format='(d4.2)'),colors.black,yshift=-0.1,xshift=-0.6
legend,string(b,m,format='("y=",d5.2,"+",d4.2,"x")'),colors.black,yshift=0.0,xshift=-0.6
legend,'(a) Goldston scaling comparison',colors.black,yshift=0.1,xshift=-0.6


;------------------------
; Regression results
;------------------------
regress     = 5.92 * (Psep/1e6)^(0.999) * (dens/1e19)^(-2.63) * (pcur/1e6)^(1.04) * (1+kappa^2)^(-1.0) * aminor^(-3.0)
regress_err = regress * sqrt( (0.999 * Psep_err / Psep)^2 + (2.63 * dens_err / dens)^2)
rsqr = r2(regress,cn_mean,b=b,m=m)

plot,regress,cn_mean,back=colors.white,col=colors.black,psym=8,yr=[0,max(cn_mean+cn_err)>max(gold)],xr=[0,max(cn_mean+cn_err)>max(gold)],$
    xtitle='5.92 P!lsep!u0.999!n I!lp!u1.04!n n!le,sep!u-2.63!n (1+k!u2!n)!u-1!n a!u-3!n',ytitle='Spectroscopy c!lN!n [%]'

errors,regress,cn_mean,ystd=cn_err,xstd=regress_err,col=colors.black
xfit = findgen(100)/99.0 * 20.0 
yfit = b + m * xfit
oplot,xfit,yfit,col=colors.black,linest=5
legend,'R!u2!n='+string(rsqr,format='(d4.2)'),colors.black,yshift=-0.1,xshift=-0.6
legend,string(b,m,format='("y=",d5.2,"+",d4.2,"x")'),colors.black,yshift=0.0,xshift=-0.6
legend,'(b) Regression',colors.black,yshift=0.1,xshift=-0.6

read,scal,prompt='Set goldston scaling factor: '
if scal ne -1 then begin
	goto,rep
endif

stop
end

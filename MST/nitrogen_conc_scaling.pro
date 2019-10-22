; Code to produce plots for cn scaling comparisons
; requires cn_database.pro to save concentrations

Pro nitrogen_conc_scaling,all=all

shots= [30776    ,34358    ,34108    ,30307    ,30306    ,30298    ,34973    ,35356    ,35158    ,30506 ]
trange=[[3.1,3.3],[4.3,5.0],[3.8,4.0],[3.4,3.6],[3.5,3.7],[3.8,4.1],[3.2,3.3],[2.0,2.5],[1.95,2.3],[2.65,2.75]] 
if keyword_set(all)then begin
	setgraphics,ncol=3,nrow=1,xs=1500,ys=500,colpick=colpick,colors=colors,collabel=collabel ,/full_list
endif else begin
	setgraphics,ncol=2,nrow=1,xs=1500,ys=500,colpick=colpick,colors=colors,collabel=collabel ,/full_list
end

if keyword_set(all)then begin
user_psym,1
plot,[0,15],[0,8],/nodata,back=colors.white,col=colors.black,xtitle='T!ldiv!n [eV]',ytitle='H-5 n!le!n [10!u19!n m!u-3!n]'
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
te=-1 & dens=-1
for i=0,n_elements(shots)-1 do begin
	augped,shots[i],nesep,tesep
	te=[te,tesep]
	dens=[dens,nesep]
endfor
te=te[1:*]
dens=dens[1:*]
user_psym,1,/fill
diag = 'BPT'
exp='davidp'
Psep=-1.0
bp=-1.0
ngw=-1.0
k=-1.0
cn_mean=-1.0
cn_err=-1.0
for i=0,n_elements(shots)-1 do begin
	; Radiated power and Psep
	read_signal_mrm,0L,shots[i],'TOT','P_TOT',time1,Ptot,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	pin = mean(ptot[id])
	if max(time1) lt trange[1,i] then begin
		read_signal_mrm,0L,shots[i],'NIS','PNI',time1,PNBI,1
		id = where(time1 ge trange[0,i] and time1 le trange[1,i])
		pni = mean(pnbi[id])
	
		read_signal_mrm,0L,shots[i],'ICP','PICRN',time1,PICRN,1
		id = where(time1 ge trange[0,i] and time1 le trange[1,i])
		pic = mean(picrn[id])

		read_signal_mrm,0L,shots[i],'ECS','PECRH',time1,PECRH,1
		id = where(time1 ge trange[0,i] and time1 le trange[1,i])
		pec = mean(pecrh[id])

		pin = pni + pic + pec
	endif
	read_signal_mrm,0L,shots[i],diag,'Pr_tot',time,Ptot,1,exp=exp
	read_signal_mrm,0L,shots[i],diag,'Pr_main',time,Pmain,1,exp=exp
	id = where(time ge trange[0,i] and time le trange[1,i])
	Pcore = mean(Pmain[id])
	Psep = [Psep,Pin - Pcore]

	; Plasma kappa, minor radius, current
	read_signal_mrm,0L,shots[i],'TOT','kappa',time1,kappa,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	k = [k,mean(kappa[id])]
	kap = mean(kappa[id])
	
	read_signal_mrm,0L,shots[i],'TOT','ahor',time1,amin,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	a = mean(amin[id])
	
	read_signal_mrm,0L,shots[i],'MAG','Ipa',time1,curr,1
	id = where(time1 ge trange[0,i] and time1 le trange[1,i])
	ip = mean(curr[id])

	bp = [bp,((4*!pi*0.0000001) * ip)/(2.0*!pi * a * sqrt(0.5*(1+kap^2)))]
	ngw= [ngw,(dens/1e20) / ((ip/1e6) / !pi / a^2)]
	
	file = 'save/cn/'+string(shots[i],format='(i5)')+'/cn_database.sav'
	restore,file
	cn_mean = [cn_mean,conc]
	cn_err = [cn_err,err]
End
Psep = Psep[1:*]
bp=bp[1:*]
ngw=ngw[1:*]
k=k[1:*]
cn_mean=cn_mean[1:*]
cn_err=cn_err[1:*]
scal=1.2
rep:
gold = scal * ((Psep/1e6) / (bp * (1+k^2)^(1.5) * ngw^2))/18.3
plot,Psep/1e6,ngw,psym=8,xr=[0,max(psep/1e6)],yr=[0,0.3],col=colors.black,back=colors.white,xtitle='P!lsep!n [MW]',ytitle='n!le,sep!n / n!lGW!n'; [10!u19!n m!u-3!n]'
plot,gold,cn_mean,back=colors.white,col=colors.black,psym=8,xr=[0,20],yr=[0,20]
err_plot,gold,cn_mean,cn_err,col=colors.black
oplot,[0,100],[0,100],col=colors.black,linest=5
for i=0,n_elements(shots)-1 do print,shots[i],cn_mean[i],gold[i],psep[i]/1e6,ngw[i],dens[i],bp[i]

read,scal,prompt='Set Goldston scaler: '
if scal ne '-1' then goto,rep

stop
end

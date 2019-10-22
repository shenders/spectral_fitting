Pro workhorse, shot, n12, n14, t1, t2, xr=xr, debug=debug, $
               combine=combine, nofit=nofit, colors=colors,$
	       xyy=xyy, choice=choice, graphtype=graphtype,$
	       showrov12=showrov12, interelm=interelm

	elmcond   = 4.5
	time0 = -1
	time1 = -1 
	read_signal_mrm,0L,shot,'MAC','Tdiv',time0,tdiv,2,exp=exp
	read_signal_mrm,0L,shot,'FPG','Suna2b',time1,spnt,2,exp=exp

	telm1  = find_elm(shot,t1)
	id1    = where(telm1 ge elmcond)
	telm2  = find_elm(shot,t2)
	id2    = where(telm2 ge elmcond)
	
	if ~keyword_set(interelm)then begin
		id1 = findgen(n_elements(t1))
		id2 = findgen(n_elements(t2))
	endif
	
	if ~keyword_set(xr)then xr=[min(t1),max(t1)]
		
	tdiv_int1 = interpol(tdiv,time0,t1[id1])	
	tdiv_int2 = interpol(tdiv,time0,t2[id2])
	
	id3 = where(time1 ge xr[0] and time1 le xr[1])

	sp_int1 = interpol(spnt,time1,t1[id1])
	sp_int2 = interpol(spnt,time1,t2[id2])
	spnt    = mean(spnt[id3])

	if ~keyword_set(nofit)then begin
		yfit12 = smooth(n12[id1],20,/edge_truncate)
		yfit14 = smooth(n14[id2],20,/edge_truncate)
		yspnt1 = smooth(sp_int1,1,/edge_truncate)
		yspnt2 = smooth(sp_int2,1,/edge_truncate)
		ytdiv1 = smooth(tdiv_int1,20,/edge_truncate)
		ytdiv2 = smooth(tdiv_int2,20,/edge_truncate)		
	endif else begin
		yfit12 = smooth(n12[id1],1,/edge_truncate)
		yfit14 = smooth(n14[id2],1,/edge_truncate)
		yspnt1 = smooth(sp_int1,1,/edge_truncate)
		yspnt2 = smooth(sp_int2,1,/edge_truncate)
		ytdiv1 = smooth(tdiv_int1,1,/edge_truncate)
		ytdiv2 = smooth(tdiv_int2,1,/edge_truncate)		
	end	
	
	if keyword_set(debug)then begin
		plot,t1,n12,back=colors.white,col=colors.black,thick=0.1
		oplot,t2,n14,col=colors.red,thick=0.1
		oplot,t1[id1],yfit12,linest=5,col=colors.blue,thick=6
		oplot,t2[id2],yfit14,linest=5,col=colors.green,thick=6
		plot,t1[id1],tdiv_int1,back=colors.white,col=colors.black,psym=8,/nodata
		oplot,t2[id2],tdiv_int2,col=colors.red,psym=8
		oplot,t2[id2],ytdiv2,col=colors.red,linest=5,thick=4
	end
	user_psym,1,/fill
	
	title=string(shot,format='("#",i5)')+' SP: '+string(spnt,format='(d7.4)')
	
	if graphtype eq 1 then begin
		title=''
		if keyword_set(combine)then oplot,ytdiv1,(yfit12/yfit14)-1,psym=8,col=choice else $
			plot,ytdiv1,(yfit12/yfit14)-1,$
			ytitle='N II: (ROV-12/ROV-14)-1',xtitle='Tdiv [eV]',yr=[-1.0,1.0],ys=1,$
				title=title,back=colors.white,col=colors.black,psym=8,xr=[0,15]
	
		if keyword_set(combine)then xyouts,[6],[xyy],[string(shot,format='(i5)')+' SP: '+string(spnt,format='(d7.4)')],$
		    col=[choice] else xyouts,[6],[xyy],[string(shot,format='(i5)')+' SP: '+string(spnt,format='(d7.4)')],$
	    		col=[choice]
		
		oplot,[0,20],[0,0],linest=5,col=colors.black
	endif
	if graphtype eq 2 then begin
		
		read_signal_mrm,0L,shot,'UVS','N_tot',time2,ntot,2,exp=exp
		
		ntot_int1 = interpol(ntot,time2,t1[id1])       
		ntot_int2 = interpol(ntot,time2,t2[id2])
		
		yfit12 = yfit12 / ntot_int1
		yfit14 = yfit14 / ntot_int2
		if ~keyword_set(combine)then begin		
		plot,ytdiv2,yfit14/max(yfit14),$
			ytitle='N II: ROV014/max(ROV014)',xtitle='Tdiv [eV]',yr=[0,1.2],ys=1,$
				back=colors.white,col=colors.black,psym=8,xr=[0,15],/nodata
		oband,[0,3,12,20],[1.02,1.02,0.72,0.72],[0.85,0.85,0.55,0.55],/norm,col=colors.black
		endif
		oplot,ytdiv2,yfit14/max(yfit14),col=[choice];,psym=-8

		xyouts,[1],[xyy],[string(shot,format='("AUG #",i5)')+' S: '+string(spnt,format='(d6.4)')],$
		    col=[choice] 
		if keyword_set(showrov12)then begin
			oplot,ytdiv1,yfit12/max(yfit12),col=[choice],psym=8
			oplot,[1,5,20],[0.5,1.05,0.85],col=colors.black,thick=4
			oplot,[1,6,20],[0.2,0.9,0.7],col=colors.black,thick=4
		endif
	endif
	if graphtype eq 3 then begin
		title = string(shot,format='(i5)')+' Tdiv: '+string(mean(tdiv_int1),format='(d7.4)')
		
		plot,yspnt1,yfit12,yr=[0,1.0e19],$
		ytitle='N II',xtitle='Strike point position [S]',$
		title=title,back=colors.white,col=colors.black,psym=8,xr=[1.02,1.08]
		oplot,yspnt2,yfit14,col=colors.red,psym=8
		
	endif
	if keyword_set(debug)then stop

End

Pro strikepoint,debug=debug,nofit=nofit,showrov12=showrov12,single=single,interelm=interelm,psplot=psplot

combine = 1

setgraphics,colors=colors,nrow=1,ncol=1,ys=600,xs=800,psplot=psplot,file='figures/strikepoint_scan.ps' 

if keyword_set(single)then begin
	restore,'save/'+string(single,format='(i5)')+'/ROV012-data.idl'
	n12  = output.nii[*,0,0]
	t1   = output.time
	restore,'save/'+string(single,format='(i5)')+'/ROV014-data.idl'
	n14 = output.nii[*,0,0]
	t2  = output.time
	workhorse, single, n12, n14, t1, t2, graphtype = 1,showrov12=showrov12,interelm=interelm, debug=debug, nofit=nofit, colors=colors, xyy=-0.3, choice=colors.black
	stop	
endif
interelm=1
xyy=[0.5,0.42,0.34,0.26,0.18,0.1,0.02]
colpick=[colors.navy,colors.blue,colors.cyan,colors.green,colors.red,colors.orange,colors.gold]
i=2
ii=0

shot = 30505
	time=-1
	time1=-1
	read_signal_mrm,0L,shot,'MAC','Tdiv',time,tdiv,2,exp=exp
	restore,'save/'+string(shot,format='(i5)')+'/ROV012-data.idl'
	n12  = output.nii[*,0,0]
	t1   = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV014-data.idl'
	n14 = output.nii[*,0,0]
	t2  = output.time	
	id = where(t1 ge 3.8 and t1 le 5.0)
	t1 = t1[id]
	t2 = t2[id]
	n12 = n12[id]
	n14 = n14[id]	
	workhorse, shot, n12, n14, t1, t2, graphtype = i,showrov12=showrov12,interelm=interelm,debug=debug, nofit=nofit, colors=colors, xyy=xyy[ii], choice=colpick[ii]
	ii = ii+1
shot = 32932

	time=-1
	time1=-1
	read_signal_mrm,0L,shot,'MAC','Tdiv',time,tdiv,2,exp=exp
	restore,'save/'+string(shot,format='(i5)')+'/ROV012-data.idl'
	n12  = output.nii[*,0,0]
	t1   = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV014-data.idl'
	n14 = output.nii[*,0,0]
	t2  = output.time	
	id = where(t1 ge 5.2 and t1 le 5.9)
	t1 = t1[id]
	t2 = t2[id]
	n12 = n12[id]
	n14 = n14[id]	
	workhorse, shot, n12, n14, t1, t2, graphtype = i,showrov12=showrov12,interelm=interelm, combine=combine,debug=debug, nofit=nofit, colors=colors, xyy=xyy[ii], choice=colpick[ii]
	ii = ii+1

shot = 30776

	restore,'save/'+string(shot,format='(i5)')+'/ROV012-data.idl'
	n12  = output.nii[*,0,0]
	t1   = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV014-data.idl'
	n14 = output.nii[*,0,0]
	t2  = output.time
	id = where(t1 ge 2.8 and t1 le 4.0)
	t1 = t1[id]
	t2 = t2[id]
	n12 = n12[id]
	n14 = n14[id]	
	workhorse, shot, n12, n14, t1, t2, graphtype = i,showrov12=showrov12,interelm=interelm, combine=combine, debug=debug, nofit=nofit, colors=colors, xyy=xyy[ii], choice=colpick[ii]
	ii = ii+1

	
shot = 30623

	time=-1
	time1=-1
	read_signal_mrm,0L,shot,'MAC','Tdiv',time,tdiv,2,exp=exp
	restore,'save/'+string(shot,format='(i5)')+'/ROV012-data.idl'
	n12  = output.nii[*,0,0]
	t1   = output.time
	restore,'save/'+string(shot,format='(i5)')+'/ROV014-data.idl'
	n14 = output.nii[*,0,0]
	t2  = output.time	
	id = where(t1 ge 4.3 and t1 le 5.7)
	t1 = t1[id]
	t2 = t2[id]
	n12 = n12[id]
	n14 = n14[id]	
	workhorse, shot, n12, n14, t1, t2, graphtype = i,showrov12=showrov12,interelm=interelm, debug=debug, combine=combine, nofit=nofit, colors=colors, xyy=xyy[ii], choice=colpick[ii]
	ii = ii+1

shot = 34358
	ch2 = 22 ; ROV-12
	ch4 = 19 ; ROV-14
	xr  = [3.75,5.8]
	t1  = -1
	read_signal_mrm,0L,shot,'EVL','N_1_3995',t1,nii,2,exp=exp
	
	id  = where(t1 ge xr[0] and t1 le xr[1])
	t1  = t1[id]
	t2  = t1
	n12 = nii[id,ch2]
	n14 = nii[id,ch4]

	workhorse, shot, n12, n14, t1, t2, graphtype = i,showrov12=showrov12,interelm=interelm, debug=debug, combine=combine, nofit=nofit, colors=colors, xyy=xyy[ii], choice=colpick[ii]
	ii = ii+1

shot = 34280
	ch2 = 22 ; ROV-12
	ch4 = 19 ; ROV-14
	xr  = [3.5,5.0]
	t1  = -1
	read_signal_mrm,0L,shot,'EVL','N_1_3995',t1,nii,2,exp=exp
	
	id  = where(t1 ge xr[0] and t1 le xr[1])
	t1  = t1[id]
	t2  = t1
	n12 = nii[id,ch2]
	n14 = nii[id,ch4]

	workhorse, shot, n12, n14, t1, t2, graphtype = i,showrov12=showrov12,interelm=interelm, debug=debug, combine=combine, nofit=nofit, colors=colors, xyy=xyy[ii], choice=colpick[ii]
	ii = ii+1


stop




















setgraphics,colors=colors,nrow=3,ncol=1,ys=1000,xs=400
shot = 35846
	ch2 = 22 ; ROV-12
	ch4 = 9 ; ROV-14
	xr  = [5.2,5.75]
	t1  = -1
	read_signal_mrm,0L,shot,'EVL','N_1_3995',t1,nii,2,exp=exp
	
	id  = where(t1 ge xr[0] and t1 le xr[1])
	t1  = t1[id]
	t2  = t1
	nii = nii[id,*]
	n12 = nii[*,ch2]
	n14 = nii[*,ch4]

	workhorse, shot, n12, n14, t1, t2, graphtype = 3,showrov12=showrov12,interelm=interelm, debug=debug, combine=combine, nofit=nofit, colors=colors, xyy=-0.5, choice=colors.blue
	
shot = 35846
	ch2 = 22 ; ROV-12
	ch4 = 9 ; ROV-14
	xr  = [4,4.7]
	t1  = -1
	read_signal_mrm,0L,shot,'EVL','N_1_3995',t1,nii,2,exp=exp
	
	id  = where(t1 ge xr[0] and t1 le xr[1])
	t1  = t1[id]
	t2  = t1
	nii = nii[id,*]
	n12 = nii[*,ch2]
	n14 = nii[*,ch4]
	workhorse, shot, n12, n14, t1, t2, graphtype = 3,showrov12=showrov12,interelm=interelm, debug=debug, combine=combine, nofit=nofit, colors=colors, xyy=-0.5, choice=colors.blue

shot = 35846
	ch2 = 22 ; ROV-12
	ch4 = 9 ; ROV-14
	xr  = [2.6,3.6]
	t1  = -1
	read_signal_mrm,0L,shot,'EVL','N_1_3995',t1,nii,2,exp=exp
	id  = where(t1 ge xr[0] and t1 le xr[1])
	t1  = t1[id]
	t2  = t1
	nii = nii[id,*]
	n12 = nii[*,ch2]
	n14 = nii[*,ch4]
	workhorse, shot, n12, n14, t1, t2, graphtype = 3,showrov12=showrov12,interelm=interelm, debug=debug, combine=combine, nofit=nofit, colors=colors, xyy=-0.5, choice=colors.blue



End

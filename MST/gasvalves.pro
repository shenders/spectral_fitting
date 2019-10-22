Pro gasvalves


shots= [30776    ,34358    ,34108    ,30307    ,30306    ,30298    ,34973    ,35356     ,35158    ,30506 ]
adas_colors,colors=colors
v_store = -1
for i=0,n_elements(shots)-1 do begin
	read_signal_mrm,0L,shots[i],'UVS','CFDu01X',t,v1,1
	read_signal_mrm,0L,shots[i],'UVS','CFDu05X',t,v2,1
	read_signal_mrm,0L,shots[i],'UVS','CFDu09X',t,v3,1
	read_signal_mrm,0L,shots[i],'UVS','CFDu13X',t,v4,1
;	plot,t,v1,col=colors.black,back=colors.white,yr=[0,max(v1)>max(v2)>max(v3)>max(v4)]
;	oplot,t,v2,col=colors.red
;	oplot,t,v3,col=colors.blue
;	oplot,t,v4,col=colors.green
	id1 = where(v1 ne 0)
	id2 = where(v2 ne 0)
	id3 = where(v3 ne 0)
	id4 = where(v4 ne 0)
	v_active = 0
	if id1[0] ne -1 then v_active = v_active + 1
	if id2[0] ne -1 then v_active = v_active + 1
	if id3[0] ne -1 then v_active = v_active + 1
	if id4[0] ne -1 then v_active = v_active + 1
	print,'Valves active: ',v_active 
	v_store = [v_store,v_active]

endfor
v_store = v_store[1:*]

for i=0,n_elements(shots)-1 do print,shots[i],v_store[i]
stop
End

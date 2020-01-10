Function read_jet_database,pulse

    openr,unit,'JET/JET_database.txt',/get_lun
    
    shot = -1
    psep = -1
    nsep = -1
    csep = -1
    twn1  = -1
    twn2  = -1
    
    ; read headers
    str = ''
    readf,unit,str
    
    while ~eof(unit)do begin
    	readf,unit,st,tw1,tw2,nseed,pt,pt_i,ps,dwdt,nes,nep,tep,cnp,ze,zv,bo,tt,ts
	shot = [shot,st]
	psep = [psep,ps]
	nsep = [nsep,nes]
	csep = [csep,cnp]
	twn1 = [twn1,tw1]
	twn2 = [twn2,tw2]
    end
    
    shot = shot[1:*]
    psep = psep[1:*]
    nsep = nsep[1:*]
    csep = csep[1:*]
    twn1 = twn1[1:*]
    twn2 = twn2[1:*]
    
    id = where(shot eq pulse)
    
    return,{shot:shot[id],$
            twn1:twn1[id],$
	    twn2:twn2[id],$
	    nsep:nsep[id],$
	    psep:psep[id],$
	    csep:csep[id]}

End

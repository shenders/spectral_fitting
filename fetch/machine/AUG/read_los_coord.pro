pro read_los_coord, shot, losnames, $
                    year=year,coord_file=coord_file,verbose=verbose, $
                    r1, z1, f1, r2, z2, f2

nlos = n_elements(losnames)
if not keyword_set(coord_file) then begin
    if keyword_set(year) then begin
        syear = string(year,f='(i4)')
    endif else begin
        if (shot ge 19584 and shot lt 20660) then syear =  '2005'
        if (shot ge 20660 and shot lt 21482) then syear =  '2006'
        if (shot ge 21482 and shot lt 22586) then syear =  '2007'
        if (shot ge 22586 and shot lt 24190) then syear =  '2008'
        if (shot ge 24190 and shot lt 25891) then syear =  '2009'
        if (shot ge 25891 and shot lt 27401) then syear =  '2010'
        if (shot ge 27401 and shot lt 28524) then syear =  '2012'
        if (shot ge 28524 and shot lt 30135) then syear =  '2013'
        if (shot ge 30135 and shot lt 31777) then syear =  '2014'
        if (shot ge 31777 or  shot lt 10000) then syear =  '2015'
    endelse
    coord_file = '/afs/ipp/u/sprd/loscoord/LOS_COORD_'+syear
endif 
exist = file_search(coord_file)
if exist[0] eq '' then begin
    print,' ERROR:  LOS_COORD file ',coord_file,' not found!'
    stop
endif


;read names and coordinates of all LOS
openr,lun,coord_file,/get_lun
cstr = ' '

;read header
readf,lun,cstr
readf,lun,cstr

r1 = fltarr(nlos)
f1 = fltarr(nlos)
z1 = fltarr(nlos)
r2 = fltarr(nlos)
f2 = fltarr(nlos)
z2 = fltarr(nlos)
found = intarr(nlos)
while not EOF(lun) do begin
    readf,lun,cstr
    astr = strsplit(cstr,'''',/extract)
    if n_elements(astr) gt 1 then begin
        idx = where(strtrim(losnames,2) eq strtrim(astr[0],2), ifound)
        if ifound eq 1 then begin
            bstr =  strsplit(astr[1],' ',/extract)
            r1[idx] = float(bstr[0])
            f1[idx] = float(bstr[1])
            z1[idx] = float(bstr[2])
            r2[idx] = float(bstr[3])
            f2[idx] = float(bstr[4])
            z2[idx] = float(bstr[5])
            found[idx] = 1
        endif 
    endif 
endwhile
close, lun
free_lun, lun

if total(found) ne nlos then begin
    wo = where(found eq 0)
    print,'Could not find ',losnames(wo)
    print,'in file: ',coord_file
endif 

end

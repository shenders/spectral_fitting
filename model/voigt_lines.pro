Function voigt_lines,species=species,instr_func=instr_func,diag=diag

if ~keyword_set(instr_func)then instr_func=0.08
if ~keyword_set(diag)then diag='EVS'
if ~keyword_set(species)then species = 'D'


; Deuterium lines
D_lines  = [396.89 , 410.07]
D_trn    = [7      , 6     ]

pos = 0.0
trn = 0.0		     

id = where(species eq 'D')
if id[0] ne -1 then begin
    pos = [pos,D_lines]
    trn = [trn,D_trn]
endif

if n_elements(pos) > 1 then begin
    pos = pos[1:*]
    trn = trn[1:*]
endif

voigt = {pos : pos , fwhml : fltarr(n_elements(pos))+0.1, fwhmg : instr_func  ,couple:fltarr(n_elements(pos))+0.0,nbalmer:trn}				

return, voigt
End

Function calc_qpll,qpllu,trange,Lzrange,nu,kappa0,cz,tu,cumtrapz=cumtrapz

    if ~keyword_set(cumtrapz) then begin
	v0   = (2 * alog10(qpllu)) - alog10(1e20)
    	v1   = alog10(int_tabulated(Trange,Lzrange * 1e20 * sqrt(Trange)) * 1e-20) + (2 * alog10(nu)) + alog10(2*kappa0*cz) + (2 * alog10(tu)) - alog10(1e20)
    	v2   = 10^(v0) - 10^(v1)
    	if 10^(v0) lt 10^(v1) then begin
		print,'Numerical error: SQRT(-ve)'
		return,-1
	endif
    endif else begin
    	v2    = fltarr(n_elements(trange))
    	v2[0] = sqrt(qpllu^2)
    	for i=1,n_elements(trange)-1 do begin
    		v0 = (2 * alog10(qpllu)) - alog10(1e20)
		v1 = alog10(int_tabulated(Trange[0:i],Lzrange[0:i] * 1e20 * sqrt(Trange[0:i])) * 1e-20) + (2 * alog10(nu)) + alog10(2*kappa0*cz) + (2 * alog10(tu)) - alog10(1e20)
		v2[i] = 10^(v0) + 10^(v1)
	endfor
    end
    return, sqrt(v2) * 1e10	
End

Function LengFunc,scen=scen,$
                  tau=tau,$
		  nu=nu,$
		  Psep=Psep,$
		  Ip=Ip,$
		  amin=amin,$
		  Bt=Bt,$
		  kappa=kappa,$
		  lhat=lhat,$
		  qstar=qstar,$
		  nsep=nsep,$
		  year=year,$
		  Tt=Tt,$
		  arne=arne

; Set input defaults and machine parameters

if ~keyword_set(tau) then  tau     = 0.1E-3
if ~keyword_set(scen)then  scen    = 1
if ~keyword_set(R)then     R       = 1.65
if ~keyword_set(eps)then   eps     = 0.33
if ~keyword_set(amin)then  amin    = 0.54
if ~keyword_set(kappa)then kappa   = 1.7
if ~keyword_set(nu)then    nu      = 3.04E19
if ~keyword_set(Psep)then  Psep    = 3.6E6
if ~keyword_set(Ip)then    Ip      = 0.84
if ~keyword_set(lhat)then  lhat    = 4.33
if ~keyword_set(qstar)then qstar   = 3.00
if ~keyword_set(fsep)then  fsep    = 1.00
if ~keyword_set(fGW)then   fGW     = 0.80
if ~keyword_set(fLH)then   fLH     = 1.00
if ~keyword_set(Psep)then  Psep    = 1.00
if ~keyword_set(Psep)then  nsep    = 1.00

; Set constants 
mu0     = 4.0 * !pi * 1e-7
kappa0  = 2000
echarge = 1.60217662e-19
mi      = 2 * 1.672619E-27

; Setup cases matching table 1 of Reinke Nucl. Fusion 57 (2017) 034004
case scen of
	1: begin
		Bt      = 1.0
	end
	2: begin
		Bt      = 3.0
	end
	3: begin
		Bt      = 6.0
	end
	4: begin
		Bt      = 9.0
	end
	5: begin
		Bt      = 12.0
	end
	Else: begin
		;Bt      = 2.5
		;amin    = 0.51
		;Ip      = 0.84
		;nsep    = 3e19
		;Psep    = 6.0
		fsep    = 1.0		
		eps     = 0.51 / R
		qstar   = !pi * amin^2 * (1+kappa^2) * Bt / mu0 / (Ip*1e6) / R
		nvol    = nsep / fsep
		S       = 2 * sqrt(2) * !pi^2 * R^2 * eps * (1+kappa^2)^(0.5)
		PLH     = 0.049 * Bt^(0.80)*(nvol/1e20)^(0.72)*S^(0.94)
		fLh     = (Psep/2.37) / PLH
		nGW     = Ip / !pi / amin^2
		fGW     = (nvol/1e20) / nGW
	end
endcase

l       = lhat * qstar * R
amin    = eps * R
Ip      = 1e-6 * !pi * amin^2 * (1+kappa^2) * Bt / mu0 / qstar / R
S       = 2 * sqrt(2) * !pi^2 * R^2 * eps * (1+kappa^2)^(0.5)
nvol    = fGW * Ip /!pi / amin^2
PLH     = 1E6 * 0.049 * Bt^(0.80)*nvol^(0.72)*S^(0.94)
Psol    = fLH * PLH 
Bp      = 2E-1 * Ip  / (amin * sqrt(0.5 * (1 + kappa^2)))
B       = SQRT(Bt^2 + Bp^2)
lambdaq = 1.35E-3 * Psol^(-0.02) * R^(0.04) * Bp^(-0.92) * eps^(0.42)
qpllu_c = 117.9 * Psol * B / R / eps^(0.42) ; Psol * B / 2 / !pi / R / lambdaq / Bp
qpllu   = 1E9 * 0.112 * Bt^(2.52) * fLH * (fGW / qstar)^(0.72) * R^(0.16) * eps^(0.52) * (1 + kappa^2)^(1.19)
qpllu   = Psol * B / 2.0 / !pi / R / lambdaq / Bp
nu      = 1E20 * 0.795 * fsep * fGW * (1 + kappa^2) * Bt / R / qstar
; Set target Te and initial cZ
if ~keyword_set(Tt)then Tt      = 5.0 ; eV    
czinit  = 0.01
qacc    = 0.01
Tuacc   = 0.01
itmax   = 1e7
; First guesses at Tu and cz:
TuDTC   = ((7.0/2.0)*qpllu*L/kappa0)^(2.0/7.0) 
Tu      = TuDTC
cz      = czinit
error1  = Tuacc+10
iter1   = 0

; Retrieve Lz data with ne * tau 
if keyword_set(arne)then begin
	openr,unit,'MST/Lz_Kallenbach_N',/get_lun
	unt = 0
	Lz_func=fltarr(200,5)
	Te_func=fltarr(200)
	while ~eof(unit) do begin
    	Lz_AK = fltarr(5)
    	readf,unit,Te_AK,Lz_AK
    	LZ_func[unt,*]=Lz_AK
    	Te_func[unt]=Te_AK
    	unt = unt+1
	end
	Lz      = [[te_func],[Lz_func[*,0]]]

endif else begin
	num     = 200
	te      = adas_vector(high=1e4,low=1,num=num)
	dens    = 1e14 + fltarr(num)
	if ~keyword_set(year)then year = '96'
	run_ionbal,prad=power,tau=tau,te=te,elem='n',dens=dens, files = {acd : '/afs/ipp/home/a/adas/adas/adf11/acd'+year+'/acd'+year+'_n.dat', $
            	    	    	    	    	    	    	 scd : '/afs/ipp/home/a/adas/adas/adf11/scd'+year+'/scd'+year+'_n.dat', $
            	    	    	    	    	    	    	 plt : '/afs/ipp/home/a/adas/adas/adf11/plt'+year+'/plt'+year+'_n.dat', $
            	    	    	    	    	    	    	 prb : '/afs/ipp/home/a/adas/adas/adf11/prb'+year+'/prb'+year+'_n.dat'}
	Lz      = [[te],[power]]
end
; Lz function must be 2D array of data vs Te

while abs(error1) gt Tuacc and iter1 le itmax do begin
    ; Pick out the particular part of the Lz function needed:
    inds    = where(Lz[*,0] ge Tt and Lz[*,0] le Tu)
    Trange  = Lz[inds,0]
    Lzrange = 10^(alog10(Lz[inds,1]))    
    ; Sheath q||, should be equal to qpllt from LM:
    qpllsheath = (7.0/2.0)*nu*Tu*echarge*sqrt(2*Tt*echarge/mi);
    ; qpllt from LM:
    qpllt = calc_qpll(qpllu,trange,Lzrange,nu,kappa0,cz,tu)    
    ; Iterate until qpllt=qpllsheath (within specified accuracy):
    error2 = qacc+10
    iter2  = 0
    while abs(error2) ge qacc and iter2 le itmax do begin
	if qpllt lt qpllsheath then begin
            cz = 0.999*cz;
        endif else begin
            cz = 1.001*cz;
        end
        qpllt = calc_qpll(qpllu,trange,Lzrange,nu,kappa0,cz,tu)
        if qpllt eq -1 then return,-1
	error2 = (qpllt-qpllsheath)/qpllsheath;
	iter2 = iter2 + 1
    end

    ; Recalculate Tu:
    qpll   = calc_qpll(qpllt,trange,Lzrange,nu,kappa0,cz,tu,/cumtrapz)
    s      = fltarr(n_elements(trange))
    for i=1,n_elements(trange)-1 do s[i] = kappa0*int_tabulated(Trange[0:i],Trange[0:i]^(5.0/2.0)/qpll[0:i]);
    Tuz    = interpol(Trange,s,L)
    Tu     = 0.8*Tu+0.2*Tuz;
    error1 = (Tu-Tuz)/Tu;
    iter1  = iter1 + 1
end
print,'cZ[%]        : ',cZ*100   
print,'nu [1E20 m-3]: ',nu
print,'Tu [eV]      : ',Tu
print,'TuDTC [eV]   : ',TuDTC
print,'q||u [GW/m2] : ',qpllu/1e9
print,'q||t [GW/m2] : ',qpllt/1e9
return,cZ
END

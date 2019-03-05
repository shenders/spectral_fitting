Function LengFunc,debug=debug,mreinke=mreinke,tau=tau,year=year,arne=arne,ylog=ylog
;qpllu,nu,Tt,kappa0,L,mi,echarge,Lz

kappa0  = 2000
kappa   = 1.7
R       = 1.65
l       = 4.33 * 3.0 * R
epsilon = 0.33
echarge = 1.60217662e-19
mi      = 2 * 1.672619E-27

if keyword_set(mreinke)then begin
    qpllu   = 2.01E9
    nu      = 1.5E20
    if ~keyword_set(tau)then tau = 0.1E-3
    print,'Using MReinke table values'
endif else begin
    Bt      = 2.5
    nu      = 2e19
    Psep    = 6E6
    Ip      = 1.0
    amin    = 0.50
    epsilon = amin / R
    Bp      = ((4.0 * !pi * 1E-7) * (Ip * 1E6)) / (2.0 * !pi * amin * sqrt(0.5 * (1 + kappa0^2)))
    B       = SQRT(Bt^2 + Bp^2)
    lambdaq = 1.35 * Psep^(-0.02)*R^(0.04)*Bp^(-0.92)*epsilon^(0.42) * 1E-3
    qpllu   = Psep * B / 2 / !pi / R / lambdaq / bp
end
Tt      = 5.0 ; eV    
czinit  = 0.01
qacc    = 0.001
Tuacc   = 0.001

; First guesses at Tu and cz:
TuDTC   = ((7.0/2.0)*qpllu*L/kappa0)^(2.0/7.0) 
Tu      = TuDTC
cz      = czinit
error1  = Tuacc+10

num     = 100
te      = adas_vector(high=1e4,low=5,num=num)
dens    = 1e14 + fltarr(num)

openr,unit,'Lz_Kallenbach_N',/get_lun
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
Print,'Using Tau = ',tau
if ~keyword_set(year)then year = '96'
run_ionbal,prad=power,tau=tau,te=te,elem='n',dens=dens, files = {acd : '/home/adas/adas/adf11/acd'+year+'/acd'+year+'_n.dat', $
            	    	    	    	    	    	    	 scd : '/home/adas/adas/adf11/scd'+year+'/scd'+year+'_n.dat', $
            	    	    	    	    	    	    	 plt : '/home/adas/adas/adf11/plt'+year+'/plt'+year+'_n.dat', $
            	    	    	    	    	    	    	 prb : '/home/adas/adas/adf11/prb'+year+'/prb'+year+'_n.dat', $
								 ccd : '/home/adas/adas/adf11/ccd'+year+'/ccd'+year+'_n.dat'  }

if keyword_set(arne)then Lz = [[te_func],[lz_func]] else Lz      = [[te],[power]]

;plot,te,energy.total*1e-6,/xlog,/ylog   
plot,te,power,/xlog ,/ylog,yr=[1e-34,1e-30],xs=1
oplot,te_func,lz_func[*,0],psym=5
; Lz function must be 2D array of data vs Te

while abs(error1) gt Tuacc do begin
    ; Pick out the particular part of the Lz function needed:
    inds    = where(Lz[*,0] ge Tt and Lz[*,0] le Tu)
    Trange  = [tt,Lz[inds,0],tu]
    Lzrange = [10^interpol(alog10(Lz[*,1]),alog10(Lz[*,0]),alog10(tt)),10^(alog10(Lz[inds,1])),10^interpol(alog10(Lz[*,1]),alog10(Lz[*,0]),alog10(tu))]
    
    ; Sheath q||, should be equal to qpllt from LM:
    qpllsheath = (7.0/2.0)*nu*Tu*echarge*sqrt(2*Tt*echarge/mi);
    ; qpllt from LM:
    qpllt      = sqrt(qpllu^2-int_tabulated(Trange,Lzrange*sqrt(Trange))*2*kappa0*nu*nu*Tu^2*cz);
    ; Iterate until qpllt=qpllsheath (within specified accuracy):
    error2 = qacc+10;
    if keyword_set(debug)then begin
        print,'Debug ...'
    	print,'---------'
    	print,'Tu   = ',Tu
    	print,'Tt   = ',Tt
    	print,'q||u = ',qpllu
    	print,'q||t = ',qpllt
	print,'q||s = ',qpllsheath
    	print,'cZ   = ',cz
	res = '' & read,res,prompt='Continue (y/n)? '
    	if res ne 'y' then stop
    endif
    while abs(error2) ge qacc do begin
        qpllt_prev = qpllt;
        qpllt = sqrt(qpllu^2-int_tabulated(Trange,Lzrange*sqrt(Trange))*2*kappa0*nu*nu*Tu^2*cz);
        if qpllt lt qpllsheath then begin
            cz = 0.999*cz;
        endif else begin
            cz = 1.001*cz;
        end
        error2 = (qpllt-qpllsheath)/qpllsheath;
	if keyword_set(debug)then print,error2,cz
    end

    ; Recalculate Tu:
    s    = fltarr(n_elements(trange))
    qpll = fltarr(n_elements(trange))
    qpll[0] = sqrt(qpllt^2)
    for i=1,n_elements(trange)-1 do qpll[i] = sqrt(qpllt^2+2*kappa0*nu*Tu^2*cz*int_tabulated(Trange[0:i],Lzrange[0:i]*sqrt(Trange[0:i]))*nu);
    for i=1,n_elements(trange)-1 do s[i]    = kappa0*int_tabulated(Trange[0:i],Trange[0:i]^(5.0/2.0)/qpll[0:i]);
    Tuz  = interpol(Trange,s,L)
    Tu = 0.8*Tu+0.2*Tuz;
    error1 = (Tu-Tuz)/Tu;
end
print,'cZ: ',cZ   
print,'Tu: ',Tu
print,'TuDTC: ',TuDTC
print,'q||t:',qpllt
Stop
END

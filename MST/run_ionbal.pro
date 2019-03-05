PRO run_ionbal,pop=pop,prad=prad,p_rad=p_rad,tau=tau,te=te,elem=elem,dens=dens,files=files,ionbal=ionbal


if ~keyword_set(te)then itval=200 else itval=n_elements(te)
if ~keyword_set(tau)then tau=1.0000e-01  / 1e3 
if ~keyword_set(te)then te    = adas_vector(low=1,high=1000,num=itval)
if ~keyword_set(dens)then dens  = fltarr(itval) + 1e14
if ~keyword_set(elem)then elem  = 'n'
if ~keyword_set(files)then files = {acd : '/home/adas/adas/adf11/acd96/acd96_n.dat', $
            	    	    	    scd : '/home/adas/adas/adf11/scd96/scd96_n.dat', $
            	    	    	    plt : '/home/adas/adas/adf11/plt96/plt96_n.dat', $
            	    	    	    prb : '/home/adas/adas/adf11/prb96/prb96_n.dat'  }

iz0   = xxeiz0(elem)

if keyword_set(ionbal)then begin
    run_adas405, uid='adas', year='96', elem=elem, te=te, dens=dens, $
             power=power, frac=frac, files=files
    pop=transpose(frac.ion)	     
    prad=power.total/1e6
    p_rad=power.ion/1e6
endif else begin

    prad   = fltarr(itval)
    pop   = fltarr(iz0+1, itval)
    prad  = fltarr(itval)
    p_rad = fltarr(itval, iz0)
    p_plt = fltarr(itval, iz0)
    time  = adas_vector(low=0, high=100, num=2, /linear)
    ionbal_td, files=files, te=te, dens=dens, tconf=tau, time=time, $
    	  pop=fraction, prad=power,p_rad=p_rad
    pop   = reform(fraction(*,*,1))
    prad  = reform(power[*,1]) / 1e6
    p_rad = reform(p_rad[*,*,1]) / 1e6
endelse
END

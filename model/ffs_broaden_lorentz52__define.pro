;+
; PROJECT:
;   FFS - Framework for Feature Synthesis (Object oriented model/fitting).
;
; NAME:
;   FFS_BROADEN_LORENTZ52
;
; PURPOSE:
;   An FFS operator element - takes output of other ffs elements and performs
;   lorentzian broadening.
;
; EXPLANATION:
;   Performs broadening operation on other FFS element results. 
;   If the results are gridded, then the input profile is broadened on each
;   pixel.
;
; USE:
;   Example of stand alone use to follow at later date. Currently only used 
;   indirectly by FFS system parser.
;
; INITIALISATION SYNTAX:
;   broaden_obj = obj_new('ffs_broaden_lorentz', fwhm=fwhm, $
;       trap=trap, debug=debug)
;   PARAMETERS:
;     unusedparserinput -   input required by the ffs_parser, but unused by
;                           this element
;   KEYWORDS:
;     fwhm	-   full width at half maximum value
;     trap	-   parameter controlling point at which evaluation is
;                 truncated [abs(x-x0) = trap * fwhm].
;     debug     -   set this keyword to enable debug information to be printed 
;                   to the terminal.
;
; PUBLIC ROUTINES:
;   In addtion to the methods listed below, this object inherits methods from
;   ffs_element - refer to this class' documentation for more details.
;
;   [calculate]
;     PURPOSE:
;       Performs lorentzian broadening operation on operand elements.
;     INPUTS:
;       in   -  An array of other FFS element evaluation output structures i.e.
;               a structure with the fields:
;                   wavelength  -   double array of the wavelength grid.
;	    	    intensity   -   double array of profile intensity values.
;	    	    gridded     -   0 or 1 to signify if intensities mapped to 
;   	    	    	            wavelength grid.	
;     OUTPUTS:
;	Returns a structure with the following fields:
;	wavelength  -	double array of the wavelength grid.	 
;	intensity   -	double array of profile intensity values.
;	gridded     -	0 or 1 to signify if intensities mapped to 
;   	    	    	wavelength grid. Always 1 in this case.
;     SIDE EFFECTS:
;	None.
;
; AUTHOR:
;   Christopher Nicholas, University of Strathclyde
;   (Adaptation of 'ffs_broaden' which was based on Andrew Meigs' GLV).
;
; VERSION HISTORY:
;   1.1  CHN 14/08/2008
;        * Initial commit to CVS.
;   1.2  CHN 01/03/2010
;        * 'calculate' now expands evaluation range over a larger wavelength
;           interval to compensate for inaccuracy where child result has been
;           truncated at the wavelength bounds. After calculation, data is then
;           trimmed to orignal size.
;   1.3 AGM 28/04/2016
;        * Copied (a year ago) from ffs_broaden_lorentz and modified
;        for lorentz with power 2 replaced by 5/2 (see ffs_lorentzian52__define.pro)
;-
function ffs_broaden_lorentz52::calculate, in   
    if n_elements(in) eq 0 then return, $
        self->seterrmsg('no operand elements for operator element to work on')
    x0 = (*in[0]).wavelength
    intensity = (*in[0]).intensity

    x = self -> getxdata()
    pars = self -> getpars()
    result = 0.0d
    fwhm = pars[0]->getvalue()
    trap = pars[1]->getvalue() 
    norm = (5.0/4.0)*((fwhm/2.0)^1.5)/gamma(2.0/5.0)/gamma(3.0/5.0)
    alp = 5.0/2.0

    padded = 0
    if (*in[0]).gridded eq 1 then begin
        expanded = self->ffs_element::expand_wavegrid(x0, intensity, fwhm)
        x = expanded.x
        x0 = expanded.x0
        intensity = expanded.intensity
        ; midpoint riemann sum
        area = intensity * [x0[1]-x0[0],0.5*(x0[2:*] - x0[0:n_elements(x0)-2]),x0[n_elements(x0)-1]-x0[n_elements(x0)-2]]
        padded = 1
    endif else area = intensity
   
 
   if fwhm ne 0.0d then begin
     for i=0, n_elements(x0)-1 do begin
      ;redefine trap in terms of full width (chn):
      result = result + (abs(x-x0[i]) le trap*fwhm) * (area[i]*norm)/( (abs(x-x0[i]))^alp + (fwhm/2.0d)^alp)
      ;result = result + (abs(x-x0[i]) le trap*(fwhm/2.0d)) * (area[i]*(fwhm/2.0d)/((x-x0[i])^2 + (fwhm/2.0d)^2)/!dpi)
     endfor
   endif else begin
      result=area*((0.*x)*(x ne x0)+(x eq x0)) ; if fwhm eq 0 return delta function
   endelse
   
    if padded eq 1 then begin
        trimmed = self->ffs_element::trim_wavegrid(x, result)
        x = trimmed.wavelength
        result = trimmed.intensity
    endif
   
    return,{wavelength:x, intensity:result, gridded:1}
end
    
pro ffs_broaden_lorentz52::cleanup
    self->ffs_element::cleanup
end

function ffs_broaden_lorentz52::init, $
    unusedparserinput, $
    fwhm=fwhm, $
    trap=trap, $
    debug=debug
        
    deftrap=40.0d0
    
    if self->ffs_element::init(debug=debug) ne 1 then return, 0
    
    if self->addpar(parname='fwhm') eq 0 then return, 0
    if self->addprop(name='trap') eq 0 then return, 0
      
    if $
        self->setparhardlimits(parname='fwhm', [0.0d,!values.d_infinity]) eq 0 $
        then return, 0

    if n_elements(fwhm) gt 0 then $
        if self->setparvals(fwhm, parname='fwhm') eq 0 then return, 0
        
    if n_elements(trap) gt 0 then begin
      if self->setparvals(trap, name='trap') eq 0 then return, 0
    endif else if self->setparvals(deftrap, parname='trap') eq 0 then return, 0
    
    self.minchildren = 1
    self.maxchildren = 32767
    return, 1
end


pro ffs_broaden_lorentz52__define
    self = {ffs_broaden_lorentz52, $
            inherits ffs_element}
end

function magnification, lines_mm, opening, cwav, image_ratio, $
                        order=order 
;     
;     calculates magnification factor of slit width
;     
;     INPUT: 
;     lines_mm : lines per mm
;     opening  : opening angle, 
;                angle between incoming and outgoing light ray [degree] 
;     cwav     : wavelength at central pixel [nm]
;     f1       : Brennweite Linse 1
;     f2       : Brennweite Linse 2
;
;     OUTPUT:
;     mag      : magnification factor of slit width in dispersion direction
;
;     g   = number of lines per nm
;     phi = opening in rad
;

m=-1
if keyword_set(order) then m=order

g = 1.e-6*lines_mm
phi_half = .5* opening/ 180. * !Pi

theta = asin(.5*m*g*cwav/cos(phi_half))
beta  =  phi_half-theta
alpha = -phi_half-theta
mag = cos(alpha)/cos(beta)*image_ratio

return, mag

end


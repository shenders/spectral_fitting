pro dispersion, lines_mm, op_ang, foc_len, cwav, pixw, $
                order=order, $
                dwdp,d2wdp2 

;     
;     calculates dispersion dWavelength/dPixw [nm/mu]
;     
;     input: 
;     lines_mm:  lines per mm
;     op_ang:    op_ang angle, 
;                angle between incoming and outgoing light ray [degree]
;     foc_len:   focal length [m]
;     cwav:      wavelength at central pixel [nm]
;     pixw:      pixel width [mu]
;     
;     output: dispersion (1st and 2nd order)
;
;
      m = -1.
      if keyword_set(order) then m = order
      pix_order = 1

      gr = lines_mm*1.e3                     ;[lines/mm -> lines/m]
      cw = cwav*1.e-9                        ;[nm -> m]
      pw = pixw*1.e-6                        ;[mu -> m]
      phi_half = 0.5*op_ang/180.*!Pi         ;[rad]
      theta = asin(0.5*m*gr*cw/cos(phi_half))
      beta = theta-phi_half

      dbeta_dpix = pix_order*(pw/foc_len)
      dwdp = cos(beta)/gr*dbeta_dpix*1e9           ;[m/pix -> nm/pix]
      d2wdp2 = -sin(beta)/m/gr	*dbeta_dpix^2.*1e9 ;[m/pix -> nm/pix]
end

Pro seeding_scan

    !PATH=!PATH + ':' + $
    expand_path( '+~cxs/utilities' ) + ':' + $
    expand_path( '+~cxs/calibration' ) + ':' + $
    expand_path( '+~cxs/instrument_data' )

    s1=get_kt3(85273,spec='kt3b')
    s2=get_kt3(85272,spec='kt3b')
    s3=get_kt3(85266,spec='kt3b')
    s4=get_kt3(85270,spec='kt3b')
    s5=get_kt3(85274,spec='kt3b')
    s6=get_kt3(85276,spec='kt3b')
    s7=get_kt3(85277,spec='kt3b')

    id_t1 = where(abs(s1.time-54.0) eq min(abs(s1.time-54.0)))
    id_t2 = where(abs(s2.time-54.0) eq min(abs(s2.time-54.0)))
    id_t3 = where(abs(s3.time-54.0) eq min(abs(s3.time-54.0)))
    id_t4 = where(abs(s4.time-54.0) eq min(abs(s4.time-54.0)))
    id_t5 = where(abs(s5.time-54.0) eq min(abs(s5.time-54.0)))
    id_t6 = where(abs(s6.time-54.0) eq min(abs(s6.time-54.0)))
    id_t7 = where(abs(s7.time-54.0) eq min(abs(s7.time-54.0)))
    
    setgraphics,xs=800,ys=600,colors=colors
    plot,s1.lamgrid,s1.phflx[*,10,id_t1[0]]/max(s1.phflx[*,10,id_t1[0]]),/ylog,back=colors.white,col=colors.black
    oplot,s2.lamgrid,s2.phflx[*,10,id_t2[0]]/max(s2.phflx[*,10,id_t2[0]]),col=colors.red
    oplot,s3.lamgrid,s3.phflx[*,10,id_t3[0]]/max(s3.phflx[*,10,id_t3[0]]),col=colors.blue
    oplot,s4.lamgrid,s4.phflx[*,10,id_t4[0]]/max(s4.phflx[*,10,id_t4[0]]),col=colors.orange
    
    id_pix = 672; 404.1 pixel    

    oplot,[s1.lamgrid[id_pix],s1.lamgrid[id_pix]],[1e-3,1e25],linest=5,col=colors.black

    setgraphics,xs=800,ys=600,colors=colors
    miss_pix = 4
    rval = s1.rval
    
    plot,rval[miss_pix:*],s1.phflx[id_pix,miss_pix:*,id_t1[0]]/max(s1.phflx[id_pix,miss_pix:*,id_t1[0]]),$
                                       xtitle='R [m]',ytitle='N II normalised',$
                                       xs=1,xr=[2.6,2.8],back=colors.white,col=colors.black
    oplot,rval[miss_pix:*],s2.phflx[id_pix,miss_pix:*,id_t2[0]]/max(s2.phflx[id_pix,miss_pix:*,id_t2[0]]),col=colors.red
    oplot,rval[miss_pix:*],s3.phflx[id_pix,miss_pix:*,id_t3[0]]/max(s3.phflx[id_pix,miss_pix:*,id_t3[0]]),col=colors.blue
    oplot,rval[miss_pix:*],s4.phflx[id_pix,miss_pix:*,id_t4[0]]/max(s4.phflx[id_pix,miss_pix:*,id_t4[0]]),col=colors.orange
    oplot,rval[miss_pix:*],s5.phflx[id_pix,miss_pix:*,id_t5[0]]/max(s5.phflx[id_pix,miss_pix:*,id_t5[0]]),col=colors.aqua
    oplot,rval[miss_pix:*],s6.phflx[id_pix,miss_pix:*,id_t6[0]]/max(s6.phflx[id_pix,miss_pix:*,id_t6[0]]),col=colors.magenta
    oplot,rval[miss_pix:*],s7.phflx[id_pix,miss_pix:*,id_t7[0]]/max(s7.phflx[id_pix,miss_pix:*,id_t7[0]]),col=colors.cyan
stop
end    
    

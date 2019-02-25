Function h98 , shot , debug=debug

read_signal_mrm,0L,shot,'GQH','Rgeo',t1,R,2
read_signal_mrm,0L,shot,'GQH','ahor',t2,a,2
read_signal_mrm,0L,shot,'GQH','k',t3,k,2
read_signal_mrm,0L,shot,'GQH','Circumf',t4,C,2
S   = R * 2 * 3.141 * c
read_signal_mrm,0L,shot,'TOT','H-0_corr',t5,N1d,2 & N1 = interpol(N1d[where(finite(N1d))],t5[where(finite(N1d))],t1)
read_signal_mrm,0L,shot,'DCN','H-0',t6,H1d,2 & H1 = interpol(H1d[where(finite(H1d))],t6[where(finite(H1d))],t1)
read_signal_mrm,0L,shot,'GQH','lenH-1',t7,lh1d,2 & lh1 = interpol(lh1d[where(finite(lh1d))],t7[where(finite(lh1d))],t1)
N   = H1/lh1
read_signal_mrm,0L,shot,'MAI','BTF',t8,Btd,2 & bt = interpol(btd[where(finite(btd))],t8[where(finite(btd))],t1)
Lh  = 0.0488 * 1e6 * ((-Bt)^0.803) * ((N/1e20)^0.717) * (S^0.941)
read_signal_mrm,0L,shot,'TTH','L2H_SCAL',t9,Lh2d,2 & lh2 = interpol(lh2d[where(finite(lh2d))],t9[where(finite(lh2d))],t1)
read_signal_mrm,0L,shot,'NIS','PNI',t10,NId,2 & ni = interpol(nid[where(finite(nid))],t10[where(finite(nid))],t1)
read_signal_mrm,0L,shot,'ICP','PICRN',t11,ICd,2 & ic = interpol(icd[where(finite(icd))],t11[where(finite(icd))],t1)
read_signal_mrm,0L,shot,'ECS','PECRH',t12,ECd,2 & ec = interpol(ecd[where(finite(ecd))],t12[where(finite(ecd))],t1)
read_signal_mrm,0L,shot,'TOT','P_OH',t13,OHd,2 & oh = interpol(ohd[where(finite(ohd))],t13[where(finite(ohd))],t1)
Pt1 = NI + IC + EC ;+ OH
Pt2 = NI + IC + EC ;+ OH
read_signal_mrm,0L,shot,'FPC','IpiFP',t14,Ipd,2 & ip = interpol(ipd[where(finite(ipd))],t14[where(finite(ipd))],t1)
hsc = 2.0^0.19 * 0.0562 * (Ip/1e6)^0.93 * (Pt2/1e6)^(-0.69) * R^1.39 * a^0.58 * k^0.78 * (N/1e19)^0.41 * (-Bt)^0.15
read_signal_mrm,0L,shot,'GQH','Wmhd',t15,WMd,2 & WM = interpol(wmd[where(finite(wmd))],t15[where(finite(wmd))],t1)
ta  = WM/Pt2
H98 = ta / hsc
time=t1
if keyword_set(debug)then begin
	window,/free
	plot,time,h98,xtitle='Time [s]',ytitle='H98',yr=[0,1.5]
endif
Return, {H98: H98,$
	 LH:LH,$
	 time:t1}
          

End

Pro stark

	read_signal_mrm,0L,36321,'EVL','Ne',time,ne0,2
	dens=ne0[*,22]
	plot,time,dens,yr=[0,2e20],xr=[3.5,5.0]
	stop
end

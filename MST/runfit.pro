PRO runfit,block

if block eq 1 then begin
	x=fetch_data(34108,'ROV-12',tr=[3.6,4.6],machine='AUG',diag='FVS',/save)
	x=fetch_data(34108,'ROV-14',tr=[3.6,4.6],machine='AUG',diag='FVS',/save)
	x=fetch_data(30776,'ROV012',tr=[2.5,4.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(30776,'ROV014',tr=[2.5,4.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(30306,'ROV012',tr=[3.0,4.5],machine='AUG',diag='EVS',/save)
	x=fetch_data(30306,'ROV014',tr=[3.0,4.5],machine='AUG',diag='EVS',/save)	
endif
if block eq 2 then begin
	x=fetch_data(33032,'ROV012',tr=[3.0,7.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(32932,'ROV012',tr=[4.5,6.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(32932,'ROV014',tr=[4.5,6.0],machine='AUG',diag='EVS',/save)
endif
if block eq 3 then begin
	x=fetch_data(30505,'ROV012',tr=[3.5,5.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(30505,'ROV014',tr=[3.5,5.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(30506,'ROV012',tr=[2.4,3.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(30506,'ROV014',tr=[2.4,3.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(32952,'ROV012',tr=[4.5,6.1],machine='AUG',diag='EVS',/save)
	x=fetch_data(32952,'ROV014',tr=[4.5,6.1],machine='AUG',diag='EVS',/save)
endif
if block eq 4 then begin
	x=fetch_data(32950,'ROV012',tr=[4.0,6.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(32950,'ROV014',tr=[4.0,6.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(34275,'ROV-11',tr=[3.5,4.0],machine='AUG',diag='EVS',/save)
	x=fetch_data(34275,'ROV-13',tr=[3.5,4.0],machine='AUG',diag='EVS',/save)
endif
if block eq 5 then begin
	x=fetch_data(32952,'ROV012',tr=[4.5,6.2],machine='AUG',diag='EVS',/save)
	x=fetch_data(32952,'ROV014',tr=[4.5,6.2],machine='AUG',diag='EVS',/save)
endif
end

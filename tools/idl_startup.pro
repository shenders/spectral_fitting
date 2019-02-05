print,'************************************************'
print,'------------------------------------------------'
print,'         Spectral fitting with FFS'
print,'Version : 1.0'
print,'Author  : Stuart Henderson '
print,'Contact : stuart.henderson@ukaea.uk'
print,'Year    : 2018'
print,'------------------------------------------------'
device, true_color=24
device, decomposed=0
device, retain=2
LOADCT,5
!path='fetch:'+!path
!path='model:'+!path
!path='tools:'+!path
!path='MST:'+!path
!path='fetch/machine:'+!path
!path='fetch/machine/JET:'+!path
!path='fetch/machine/AUG:'+!path
.r library.pro
print,'************************************************'
print,'Start by running eg. AUG:'
print,"IDL> x=fetch_data(35158,'ROV-12',tr=[2.0,6.2],machine='AUG',diag='EVS')"
print,' For more info readme.txt'
print,'************************************************'

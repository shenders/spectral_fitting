;$Id: //depot/Release/ENVI51_IDL83/idl/idldir/lib/eigenvec.pro#1 $
;
; Copyright (c) 1994-2013, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;       EIGENVEC
;
; PURPOSE:
;       This function computes the eigenvectors of an N by N real, non-
;       symmetric array using inverse subspace iteration. The result is 
;       a complex array with a column dimension equal to N and a row 
;       dimension equal to the number of eigenvalues.
;
; CATEGORY:
;       Linear Algebra / Eigensystems
;
; CALLING SEQUENCE:
;       Result = Eigenvec(A, Eval)
;
; INPUTS:
;       A:    An N by N nonsymmetric array of type float or double.
;
;    EVAL:    An N-element complex vector of eigenvalues.
;
; KEYWORD PARAMETERS:
;       DOUBLE:  If set to a non-zero value, computations are done in
;                double precision arithmetic.
;
;        ITMAX:  The number of iterations performed in the computation
;                of each eigenvector. The default value is 4.
;
;     RESIDUAL:  Use this keyword to specify a named variable which returns
;                the residuals for each eigenvalue/eigenvector(lambda/x) pair.
;                The residual is based on the definition Ax - (lambda)x = 0
;                and is an array of the same size and type as RESULT. The rows
;                this array correspond to the residuals for each eigenvalue/
;                eigenvector pair. 
;
; EXAMPLE:
;       Define an N by N real, nonsymmetric array.
;         a = [[1.0, -2.0, -4.0,  1.0], $
;              [0.0, -2.0,  3.0,  4.0], $
;              [2.0, -6.0, -1.0,  4.0], $
;              [3.0, -3.0,  1.0, -2.0]]
;
;       Compute the eigenvalues of A using double-precision complex arithmetic.
;         eval = HQR(ELMHES(a), /double)
;
;       Print the eigenvalues. The correct solution should be:
;       (0.26366259, -6.1925899), (0.26366259, 6.1925899), $
;       (-4.9384492,  0.0000000), (0.41112397, 0.0000000)
;         print, eval
;
;       Compute the eigenvectors of A. The eigenvectors are returned in the 
;       rows of EVEC.
;         result = EIGENVEC(a, eval, residual = residual)
;
;       Print the eigenvectors.
;         print, evec(*,0), evec(*,1), evec(*,2), evec(*,3)
;
;       The accuracy of each eigenvalue/eigenvector (lamda/x) 
;       pair may be checked by printing the residual array. This array is the
;       same size and type as RESULT and returns the residuals as its rows.
;       The residual is based on the mathematical definition of an eigenvector,
;       Ax - (lambda)x = 0.
;
; PROCEDURE:
;       EIGENVEC computes the set of eigenvectors that correspond to a given 
;       set of eigenvalues using Inverse Subspace Iteration. The eigenvectors 
;       are computed up to a scale factor and are of Euclidean length. The
;       existence and uniqueness of eigenvectors are not guaranteed.
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, December 1994
;       Modified:    GGS, RSI, April 1996
;                    Modified keyword checking and use of double precision. 
;-

FUNCTION EigenVec_2, A, Eval, Double = Double, ItMax = ItMax, $
                                             Residual = Residual

  ON_ERROR, 2  ;Return to caller if error occurs.

  if N_PARAMS() ne 2 then $
    MESSAGE, "Incorrect number of input arguments."
    
  TypeA = SIZE(A)
  TypeEval = SIZE(Eval)

  if TypeA[1] ne TypeA[2] then $
    MESSAGE, "Input array must be square."

  if TypeA[3] ne 4 and TypeA[3] ne 5 then $
    MESSAGE, "Input array must be float or double."

  if TypeEval[TypeEval[0]+1] ne 6 and TypeEval[TypeEval[0]+1] ne 9 then $
    MESSAGE, "Eigenvalues must be complex or double-complex."

  Enum = TypeEval[TypeEval[0]+2] ;Number of eigenvalues.
  if TypeA[2] ne Enum then $
    MESSAGE, "Input array and eigenvalues are of incompatible size."

  ;If the DOUBLE keyword is not set then the internal precision and
  ;result are determined by the type of input.
  if N_ELEMENTS(Double) eq 0 then $
    Double = (TypeA[TypeA[0]+1] eq 5 or TypeEval[TypeEval[0]+1] eq 9)

  if N_ELEMENTS(ItMax) eq 0 then ItMax = 4

  Diag = LINDGEN(TypeA[1]) * (TypeA[1]+1) ;Diagonal indices.

  ;Double Precision.
  if Double ne 0 then begin
    Evec = DCOMPLEXARR(TypeA[1], Enum) ;Eigenvector storage array with number
                                    ;of rows equal to number of eigenvalues.
    for k = 0, Enum - 1 do begin
      Alud = A  ;Create a copy of the array for next eigenvalue computation.
      if IMAGINARY(Eval[k]) ne 0 then begin ;Complex eigenvalue.
        Alud = DCOMPLEX(Alud)
        Alud[Diag] = Alud[Diag] - Eval[k]
        ;Avoid intermediate variables. re = DOUBLE(Alud) im = IMAGINARY(Alud)
        Comp = [[DOUBLE(Alud), -IMAGINARY(Alud)], $
                [IMAGINARY(Alud), DOUBLE(Alud)]]
        ;Initial eigenvector.
        B = REPLICATE(1.0d, 2*TypeA[1]) / SQRT(2.0d * TypeA[1])
        LUDC, Comp, Index, DOUBLE = DOUBLE
        it = 0
        while it lt ItMax do begin ;Iteratively compute the eigenvector.
          X = LUSOL(Comp, Index, B, DOUBLE = DOUBLE)
          B = X / SQRT(TOTAL(X^2, 1, DOUBLE = DOUBLE)) ;Normalize eigenvector.
          it = it + 1
        endwhile
        ;Row vector storage.
        Evec[*, k] = DCOMPLEX(B[0:TypeA[1]-1], B[TypeA[1]:*])
      endif else begin ;Real eigenvalue
        Alud[Diag] = Alud[Diag] - DOUBLE(Eval[k])
	B = REPLICATE(1.0d, TypeA[1]) / SQRT(TypeA[1]+0.0d)
        LUDC, Alud, Index, DOUBLE = DOUBLE,interchanges=test
        print,index,test
	it = 0
        while it lt ItMax do begin
          X = LUSOL(Alud, Index, B, DOUBLE = DOUBLE)
          B = X / SQRT(TOTAL(X^2, 1, DOUBLE = DOUBLE)) ;Normalize eigenvector.
          it = it + 1
        endwhile
        Evec[*, k] = DCOMPLEX(B, 0.0d0) ;Row vector storage.
      endelse
    endfor
    
    if ARG_PRESENT(Residual) then begin ;Compute eigenvalue/vector residuals.
      Residual = DCOMPLEXARR(TypeA[1], Enum) ;Dimensioned the same as Evec.
	for k = 0, Enum - 1 do $
          Residual[*,k] = (A##Evec[*,k]) - (Eval[k] * Evec[*,k])
    endif
  endif else begin ;Single Precision.
    Evec = COMPLEXARR(TypeA[1], Enum) ;Eigenvector storage array.
    for k = 0, Enum - 1 do begin
      Alud = A  ;Create a copy of the array for next eigenvalue computation.
      if IMAGINARY(Eval[k]) ne 0 then begin ;Complex eigenvalue.
        Alud = COMPLEX(Alud)
        Alud[Diag] = Alud[Diag] - Eval[k]
        ;Avoid intermediate variables. re = FLOAT(Alud) im = IMAGINARY(Alud)
        Comp = [[FLOAT(Alud), -IMAGINARY(Alud)], $
                [IMAGINARY(Alud), FLOAT(Alud)]]
        ;Initial eigenvector.
        B = REPLICATE(1.0, 2*TypeA[1]) / SQRT(2.0 * TypeA[1])
        LUDC, Comp, Index, DOUBLE = DOUBLE
        it = 0
        while it lt ItMax do begin ;Iteratively compute the eigenvector. 
          X = LUSOL(Comp, Index, B, DOUBLE = DOUBLE)
          B = X / SQRT(TOTAL(X^2, 1)) ;Normalize eigenvector.
          it = it + 1
        endwhile
        ;Row vector storage.
        Evec[*, k] = COMPLEX(B[0:TypeA[1]-1], B[TypeA[1]:*])
      endif else begin ;Real eigenvalue 
        Alud[Diag] = Alud[Diag] - FLOAT(Eval[k])
        B = REPLICATE(1.0, TypeA[1]) / SQRT(TypeA[1])
        LUDC, Alud, Index, DOUBLE = DOUBLE
        it = 0
        while it lt ItMax do begin
          X = LUSOL(Alud, Index, B, DOUBLE = DOUBLE)
          B = X / SQRT(TOTAL(X^2, 1))  ;Normalize eigenvector.
          it = it + 1
        endwhile
        Evec[*, k] = COMPLEX(B, 0.0) ;Row vector storage.
      endelse
    endfor
    if ARG_PRESENT(Residual) then begin ;Compute eigenvalue/vector residuals.
      Residual = COMPLEXARR(TypeA[1], Enum) ;Dimensioned the same as Evec.
	for k = 0, Enum - 1 do $
          Residual[*,k] = (A##Evec[*,k]) - (Eval[k] * Evec[*,k])
    endif
  endelse
  
  if Double eq 0 then RETURN, COMPLEX(Evec) else RETURN, Evec

END
;$Id: //depot/Release/ENVI50_IDL82/idl/idldir/lib/cramer.pro#1 $
;
; Copyright (c) 1994-2012, Exelis Visual Information Solutions, Inc. All
;       rights reserved. Unauthorized reproduction is prohibited.
;+
; NAME:
;       CRAMER
;
; PURPOSE:
;       This function solves an n by n linear system of equations 
;       using Cramer's rule.
;
; CATEGORY:
;       Linear Algebra.
;
; CALLING SEQUENCE:
;       Result = CRAMER(A, B)
;
; INPUTS:
;       A:      An N by N array of type: float, or double.
;
;       B:      An N-element vector of type: float, or double.
;
; KEYWORD PARAMETERS:
;       DOUBLE: If set to a non-zero value, computations are done in
;               double precision arithmetic.
;
;       ZERO:   Use this keyword to set the value of floating-point
;               zero. A floating-point zero on the main diagonal of
;               a triangular matrix results in a zero determinant.
;               A zero determinant results in a 'Singular matrix'
;               error and stops the execution of CRAMER.PRO.
;               For single-precision inputs, the default value is 
;               1.0e-6. For double-precision inputs, the default value 
;               is 1.0e-12.
;
; EXAMPLE:
;       Define an array (a).
;         a = [[ 2.0,  1.0,  1.0], $
;              [ 4.0, -6.0,  0.0], $
;              [-2.0,  7.0,  2.0]]
;
;       And right-side vector (b).
;         b = [3.0, 10.0, -5.0]
;
;       Compute the solution of the system, ax = b.
;         result = cramer(a, b)
;
; PROCEDURE:
;       CRAMER.PRO uses ratios of column-wise permutations of the array (a)
;       to calculate the solution vector (x) of the linear system, ax = b.
;
; REFERENCE:
;       ADVANCED ENGINEERING MATHEMATICS (seventh edition)
;       Erwin Kreyszig
;       ISBN 0-471-55380-8
;
; MODIFICATION HISTORY:
;       Written by:  GGS, RSI, February 1994
;       Modified:    GGS, RSI, November 1994
;                    Added support for double precision results.
;       Modified:    GGS, RSI, April 1996
;                    Modified keyword checking and use of double precision.
;-

FUNCTION Cramer_2, A, B, Double = Double, Zero = Zero

  ;ON_ERROR, 2  ;Return to caller if error occurs.

  TypeA = SIZE(A)
  TypeB = SIZE(B)

  if TypeA[1] ne TypeA[2] then $
    MESSAGE, "Input array must be square."

  if TypeA[3] ne 4 and TypeA[3] ne 5 then $
    MESSAGE, "Input array must be float or double."

  ;If the DOUBLE keyword is not set then the internal precision and
  ;result are identical to the type of input.
  if N_ELEMENTS(Double) eq 0 then $
    Double = (TypeA[TypeA[0]+1] eq 5 or TypeB[TypeB[0]+1] eq 5) 

  if N_ELEMENTS(Zero) eq 0 and Double eq 0 then $
    Zero = 1.0e-6  ;Single-precision zero.
  if N_ELEMENTS(Zero) eq 0 and Double ne 0 then $
    Zero = 1.0d-12 ;Double-precision zero.

  DetermA = DETERM(A, Double = Double, Zero = Zero)
  if DetermA eq 0 then MESSAGE, "Input array is singular."

  if Double eq 0 then xOut = FLTARR(TypeA[1]) $ 
  else xOut = DBLARR(TypeA[1])

  for k = 0, TypeA[1]-1 do begin
    ColumnK = A[k,*] ;Save the Kth column of a.
    a[k,*] = B ;Permute the Kth column of A with B.
               ;Solve for the Kth component of the solution xOut
    xOut[k] = DETERM(A, Double = Double, Zero = Zero) / DetermA
    a[k,*] = ColumnK ;Restore A to its original state.
  endfor

  RETURN, xOut

END
;--------------------------------------------------------------------------------
; Time dependent ionisation balance (eigen-value method)
;
; INPUT:
;
; files  - structure with acd and scd adf11 dataset names
;          metastable resolved excluded (no qcd or xcd)
; te     - vector of temperature (eV)
; dens   - vector of electron density (cm-3)
; time   - vector of times (sec)
; meta   - initial distribution of stages
; tconf  - confinement time (sec) if not set defaults to infinity
;
;
; OUTPUT:
;
; pop    - population fractions for each time.
; eig    - structure with reaction matrix, eigenvalues with eigenvectors and
;          simultaneous equations solution.
; pdiff  - structure of time integrals of population difference from equilibrium
;          to last time-point in 'time' (pd_tmax as used in 406 for gcf) and
;          to infinity (pd_inf).
;
; NOTES:
;
; The functioning and argument list are modelled on adas405 but the files
; must be given as a files structure and all four must be supplied so
; calculating just the ionisation balance is not possible.
;
; EXAMPLE:
;
; files = { acd : '/home/adas/adas/adf11/acd96/acd96_c.dat', $
;           scd : '/home/adas/adas/adf11/scd96/scd96_c.dat', $
;           plt : '/home/adas/adas/adf11/plt96/plt96_c.dat', $
;           prb : '/home/adas/adas/adf11/prb96/prb96_c.dat'  }
; te    = [30.0, 100.0, 300.0]
; dens  = fltarr(3) + 1e12
; meta  = [0,0,1.0,0,0,0,0]
; time  = adas_vector(low=0,high=2.0e-2,num=41,/linear)
;
; ionbal_td, files=files, te=te, dens=dens, meta=meta, $
;            time=time, pop=pop, pdiff=pd, eigen=eigen
;
;
; WRITTEN  : Martin O'Mullane
;
; MODIFIED:
;       1.1     Martin O'Mullane
;                - First version with Te/dens vectors.
;       1.2     Martin O'Mullane
;                - Do not use Cramer method for Z0 above 10.
;       1.3     Martin O'Mullane
;                - Documentation fixes.
; VERSION:
;       1.1     17-04-2011
;       1.2     13-07-2016
;       1.3     24-04-2017
;
;--------------------------------------------------------------------------------

PRO ionbal_td, files = files, $
               te    = te,    $
               dens  = dens,  $
               time  = time,  $
               meta  = meta,  $
               pop   = pop,   $
               prad  = prad,  $
               p_rad = p_rad, $
               p_plt = p_plt, $
               pdiff = pdiff, $
               eigen = eigen, $
               tconf = tconf

; Extract info from inputs and setup variables

a11_acd = files.acd
a11_scd = files.scd
a11_plt = files.plt
a11_prb = files.prb

xxdata_11, file=a11_acd, class='acd', fulldata=all

iz0     = all.iz0
n_lev   = iz0 + 1           ; neutral to fully stripped
n_times = n_elements(time)

; If no metastables chosen assume all in the neutral

if n_elements(meta) EQ 0 then begin

   meta = dblarr(iz0+1)
   meta[0] = 1.0
endif


; Number of Te/dens pairs

itval = n_elements(te)
if itval NE n_elements(dens) then message, 'te and dens must be same size'

pop_o    = dblarr(n_lev, itval)
pd_max_o = dblarr(n_lev, itval)
pd_inf_o = dblarr(n_lev, itval)


; Get ionisation and recombination data from adf11 datasets


acd = dblarr(n_lev, itval)
scd = dblarr(n_lev, itval)
plt = dblarr(n_lev, itval)
prb = dblarr(n_lev, itval)

for j = 0, n_lev-2 do begin

   iz1 = j + 1
   read_adf11, file=a11_acd, class='acd', iz0=iz0, iz1=iz1, $
               te=te, dens=dens, data=tmp
   acd[j,*] = tmp

   read_adf11, file=a11_scd, class='scd', iz0=iz0, iz1=iz1, $
               te=te, dens=dens, data=tmp
   scd[j,*] = tmp

   read_adf11, file=a11_plt, class='plt', iz0=iz0, iz1=iz1, $
               te=te, dens=dens, data=tmp
   plt[j,*] = tmp

   read_adf11, file=a11_prb, class='prb', iz0=iz0, iz1=iz1, $
               te=te, dens=dens, data=tmp
   prb[j,*] = tmp

endfor


; Extra confinement loss time

if n_elements(tconf) GT 0 then begin
   aconf = tconf
   for i=0,n_elements(tconf)-1 do begin
    	if tconf(i) GE 100.0 then aconf(i) = 0.0 else aconf(i) = 1.0 / tconf(i)
   endfor	
endif else aconf = 0.0
if n_elements(aconf) ne n_lev then aconf = fltarr(n_lev) + aconf

; Evaluate for each time

pop     = dblarr(n_lev, itval, n_times)
prad    = dblarr(itval, n_times)
p_rad   = dblarr(itval, n_lev, n_times)
p_plt   = dblarr(itval, n_lev, n_times)

for ite = 0, itval-1 do begin

   ; Setup rate equations; dn/dt = An (n vector of stage populations, A nxn matrix)
   ; Note that IDL array index starts at 0

   a = dblarr(n_lev, n_lev)

   a[1,0] = dens[ite]*scd[0,ite]
   a[0,1] = dens[ite]*acd[0,ite]
   a[0,0] = -a[1,0] - aconf(0)

   for i = 1, n_lev-2 do begin
     a[i+1,i] = dens[ite]*scd[i,ite]
     a[i,i+1] = dens[ite]*acd[i,ite]
     a[i,i]   = -a[i+1,i] - a[i-1,i] - aconf(i)
   endfor

   a[n_lev-1, n_lev-1] = -a[n_lev-2,n_lev-1] - aconf(n_lev-1)

   for i = 0, n_lev-1 do a[0,i] = a[0,i] + aconf(i)

   ; The matrix 'a' was setup above following mathematical convention of a[row, column]
   ; IDL assumes a[column, row] so transpose 'a' to keep IDL happy.

   a = transpose(a)

   ; Compute eigenvalues (eval) and eigenvectors (evec) with built-in IDL routines

   eval = hqr(elmhes(a, /double), /double)
   if ite eq 3 then begin
;   print,a
   evec = eigenvec_2(a, eval,double=1,residual=resid)
   endif
   evec = eigenvec(a, eval,double=1,residual=resid)
;   if ite eq 3 then begin
;    print,resid
;    print,evec
;    print,'-------'
;   endif

   ; Our equation should not give complex results so take real part only.
   ; It may be a good idea to check if this is true!

   eval = real_part(eval)
   evec = real_part(evec)


   ; It is worthwile checking the results using the eigenvalue equation:
   ; ie does a*z = lambda*z ?

   maxerr = 0.0D0
   for i = 0, n_lev - 1 do begin
      alhs = a ## evec[*,i]
      arhs = eval[i] * evec[*,i]
      maxerr = maxerr > max(abs(alhs - arhs))
   endfor
   if maxerr GT 1.0D-6 then begin
      message, 'Solution not valid for index' + string(ite), /continue
      message, 'maxerr is ' + string(maxerr), /continue
   endif


   ; Now solve the simultaneous equations with our 'meta' starting fraction
   ; Use built-in IDL Cramer rule procedure - again evec needs to be rotated.


   ; Take 'num' elements either side of main eigenvector - set to 80
   ; effectively turning off this option.

   if n_lev GT 10 then begin

      num = 80

      res      = max(eval, dmax)
      i1       = max([0,dmax - num])
      i2       = min([dmax+num, n_lev-1])
      nn       = i2 - i1 + 1
      c        = dblarr(n_lev)
      a_mat    = transpose(evec[i1:i2, i1:i2])
      b_vec    = [1.0D0, dblarr(i2-i1)]

      c1       = invert(a_mat)##b_vec
      c[i1:i2] = c1

   endif else begin
      rhs = double(meta)
      c   = cramer_2(transpose(evec), rhs, double=1)
   endelse

   ; Now calculate population of the stages at each time

   for it = 0, n_times-1 do begin

      tm = time[it]

      for j = 0, n_lev-1 do begin
         pop[*, ite, it] = pop[*, ite, it] + c[j] * evec[*, j] * exp(eval[j] * tm)
      endfor

   endfor

   pop = pop > 1e-45

   ; Radiated power

   for j = 0, iz0-1 do prad[ite,*] = prad[ite,*] + $
                                     pop[j, ite,*] * plt[j, ite] + $
                                     pop[j+1, ite, *] * prb[j, ite]
    
   ; Radiated power - ion resolved

   for j = 0, iz0-1 do p_rad[ite,j,*] = pop[j, ite,*] * plt[j, ite] 

   
   ; Finally the integrals

   pd_max = 0.0D0
   pd_inf = 0.0D0
   tm_max = time[n_times-1]

   for j = 0, n_lev-1 do begin
      pd_max = pd_max + (c[j] * evec[*, j] / eval[j] * (exp(eval[j] * tm_max) - 1.0D0))
      pd_inf = pd_inf - (c[j] * evec[*, j] / eval[j])
   endfor

  ; Store in output arrays

  pd_max_o[*,ite] = pd_max
  pd_inf_o[*,ite] = pd_inf

endfor

; Output integrals

if arg_present(pop) then pop = pop

pdiff = { pd_max : pd_max_o, $
          pd_inf : pd_inf_o  }


; Return the eigenvectors for completness

eigen = { a   : a,    $
          vec : evec, $
          val : eval, $
          c   : c     }

END

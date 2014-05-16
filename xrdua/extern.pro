;    This file is part of XRDUA.
;
;    XRDUA is free software: you can redistribute it and/or modify
;    it under the terms of the GNU General Public License as published by
;    the Free Software Foundation, either version 3 of the License, or
;    (at your option) any later version.
;
;    XRDUA is distributed in the hope that it will be useful,
;    but WITHOUT ANY WARRANTY; without even the implied warranty of
;    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
;    GNU General Public License for more details.
;
;    You should have received a copy of the GNU General Public License
;    along with XRDUA.  If not, see <http://www.gnu.org/licenses/>.
;
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mrandomn, seed, covar, nrand

if n_params() lt 2 then begin
    print, 'Syntax- Result = mrandomn( seed, covar, [nrand] )'
    return, 0
endif

;check inputs and set up defaults
if n_elements(nrand) eq 0 then nrand = 1
if size(covar, /n_dim) ne 2 then begin
    print, 'COVAR must be a matrix.'
    return, 0
endif

np = (size(covar))[1]
if (size(covar))[2] ne np then begin
    print, 'COVAR must be a square matrix.'
    return, 0
endif

diag = lindgen(np) * (np + 1L)
epsilon = randomn(seed, nrand, np) ;standard normal random deviates (NP x NRAND matrix)

A = covar  ;store covariance into dummy variable for input into TRIRED

choldc, A, P, /double           ;do Cholesky decomposition

for j = 0, np - 1 do for i = j, np - 1 do A[i,j] = 0d

A[diag] = P

;transform standard normal deviates so they have covariance matrix COVAR
epsilon = A ## epsilon

return, epsilon
end;function mrandomn
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hist_nd,V,bs,MIN=mn,MAX=mx,NBINS=nbins,REVERSE_INDICES=ri
  s=size(V,/DIMENSIONS)
  if n_elements(s) ne 2 then message,'Input must be N (dimensions) x P (points)'
  if s[0] gt 8 then message, 'Only up to 8 dimensions allowed'

  imx=max(V,DIMENSION=2,MIN=imn)

  if n_elements(mx) eq 0 then mx=imx
  if n_elements(mn) eq 0 then mn=imn

  if s[0] gt 1 then begin
     if n_elements(mn)    eq 1 then mn=replicate(mn,s[0])
     if n_elements(mx)    eq 1 then mx=replicate(mx,s[0])
     if n_elements(bs)    eq 1 then bs=replicate(bs,s[0])
     if n_elements(nbins) eq 1 then nbins=replicate(nbins,s[0])
  endif

  if ~array_equal(mn le mx,1b) then $
     message,'Min must be less than or equal to max.'

  if n_elements(bs) eq 0 then begin
     if n_elements(nbins) ne 0 then begin
        nbins=long(nbins)       ;No fractional bins, please
        bs=float(mx-mn)/nbins   ;a correct formulation
     endif else message,'Must pass either binsize or NBINS'
  endif else nbins=long((mx-mn)/bs+1)

  total_bins=product(nbins,/PRESERVE_TYPE) ;Total number of bins
  h=long((V[s[0]-1,*]-mn[s[0]-1])/bs[s[0]-1])
  ;; The scaled indices, s[n]+N[n-1]*(s[n-1]+N[n-2]*(s[n-2]+...
  for i=s[0]-2,0,-1 do h=nbins[i]*h + long((V[i,*]-mn[i])/bs[i])

  out_of_range=[~array_equal(mn le imn,1b),~array_equal(mx ge imx,1b)]
  if ~array_equal(out_of_range,0b) then begin
     in_range=1
     if out_of_range[0] then $  ;out of range low
        in_range=total(V ge rebin(mn,s,/SAMP),1,/PRESERVE_TYPE) eq s[0]
     if out_of_range[1] then $  ;out of range high
        in_range AND= total(V le rebin(mx,s,/SAMP),1,/PRESERVE_TYPE) eq s[0]
     h=(temporary(h) + 1L)*temporary(in_range) - 1L
  endif

  ret=make_array(TYPE=3,DIMENSION=nbins,/NOZERO)
  if arg_present(ri) then $
     ret[0]=histogram(h,MIN=0L,MAX=total_bins-1L,REVERSE_INDICES=ri) $
  else $
     ret[0]=histogram(h,MIN=0L,MAX=total_bins-1L)
  return,ret
end;function hist_nd
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Bsort, Array, Asort, REVERSE = rev

        N = N_elements( Array )
        if N lt 1 then begin
                return, [0L]
           endif

        if N lt 2 then begin
            asort = array       ;MDM added 24-Sep-91
            return,[0L]    ;Only 1 element
        end
;
; sort array (in descending order if REVERSE keyword specified )
;
        subs = sort( Array )
        if keyword_set( REV ) then subs = rotate(subs,5)
        Asort = Array[subs]
;
; now sort subscripts into ascending order
; when more than one Asort has same value
;
             weq = where( (shift( Asort, -1 ) eq Asort) , Neq )

        if (Neq EQ n) then return,lindgen(n) ;Array is degenerate equal values

        if (Neq GT 0) then begin

                if (Neq GT 1) then begin              ;find clumps of equality

                        wclump = where( (shift( weq, -1 ) - weq) GT 1, Nclump )
                        Nclump = Nclump + 1

                  endif else Nclump = 1

                if (Nclump LE 1) then begin
                        Clump_Beg = 0
                        Clump_End = Neq-1
                  endif else begin
                        Clump_Beg = [0,wclump+1]
                        Clump_End = [wclump,Neq-1]
                   endelse

                weq_Beg = weq[ Clump_Beg ]              ;subscript ranges
                weq_End = weq[ Clump_End ] + 1          ; of Asort equalities.

                for ic = 0L, Nclump-1 do begin          ;sort each clump.

                        subic = subs[ weq_Beg[ic] : weq_End[ic] ]
                        subs[ weq_Beg[ic] ] = subic[ sort( subic ) ]
                  endfor

                if N_params() GE 2 then Asort = Array[subs]     ;resort array.
           endif

return, subs
end;function Bsort
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frebin,image,nsout,nlout,total=total
;+
; NAME:
;   FREBIN
;
; PURPOSE:
;   Shrink or expand the size of an array an arbitary amount using interpolation
;
; EXPLANATION: 
;   FREBIN is an alternative to CONGRID or REBIN.    Like CONGRID it
;   allows expansion or contraction by an arbitary amount. ( REBIN requires 
;   integral factors of the original image size.)    Like REBIN it conserves 
;   flux by ensuring that each input pixel is equally represented in the output
;   array.       
;
; CALLING SEQUENCE:
;   result = FREBIN( image, nsout, nlout, [ /TOTAL] )
;
; INPUTS:
;    image - input image, 1-d or 2-d numeric array
;    nsout - number of samples in the output image, numeric scalar
;
; OPTIONAL INPUT:
;    nlout - number of lines in the output image, numeric scalar
;            If not supplied, then set equal to 1
;
; OPTIONAL KEYWORD INPUTS:
;   /total - if set, the output pixels will be the sum of pixels within
;          the appropriate box of the input image.  Otherwise they will
;          be the average.    Use of the /TOTAL keyword conserves surface flux.
; 
; OUTPUTS:
;    The resized image is returned as the function result.    If the input
;    image is of type DOUBLE or FLOAT then the resized image is of the same
;    type.     If the input image is BYTE, INTEGER or LONG then the output
;    image is usually of type FLOAT.   The one exception is expansion by
;    integral amount (pixel duplication), when the output image is the same
;    type as the input image.  
;     
; EXAMPLE:
;     Suppose one has an 800 x 800 image array, im, that must be expanded to
;     a size 850 x 900 while conserving surface flux:
;
;     IDL> im1 = frebin(im,850,900,/total) 
;
;     im1 will be a 850 x 900 array, and total(im1) = total(im)
; NOTES:
;    If the input image sizes are a multiple of the output image sizes
;    then FREBIN is equivalent to the IDL REBIN function for compression,
;    and simple pixel duplication on expansion.
;
;    If the number of output pixels are not integers, the output image
;    size will be truncated to an integer.  The platescale, however, will
;    reflect the non-integer number of pixels.  For example, if you want to
;    bin a 100 x 100 integer image such that each output pixel is 3.1
;    input pixels in each direction use:
;           n = 100/3.1   ; 32.2581
;          image_out = frebin(image,n,n)
;
;     The output image will be 32 x 32 and a small portion at the trailing
;     edges of the input image will be ignored.
; 
; PROCEDURE CALLS:
;    None.
; HISTORY:
;    Adapted from May 1998 STIS  version, written D. Lindler, ACC
;    Added /NOZERO, use INTERPOLATE instead of CONGRID, June 98 W. Landsman  
;    Fixed for nsout non-integral but a multiple of image size  Aug 98 D.Lindler
;    DJL, Oct 20, 1998, Modified to work for floating point image sizes when
;        expanding the image. 
;    Improve speed by addressing arrays in memory order W.Landsman Dec/Jan 2001
;-
;----------------------------------------------------------------------------
      if N_params() LT 1 then begin
           print,'Syntax = newimage = FREBIN(image, nsout, nlout, [/TOTAL])'  
           return,-1
       endif

       if n_elements(nlout) eq 0 then nlout=1
;
; determine size of input image
;
    ns = n_elements(image[*,0])
    nl = n_elements(image)/ns
;
; determine if we can use the standard rebin function
;
        dtype = size(image,/TNAME)
    if dtype EQ 'DOUBLE' then begin
        sbox = ns/double(nsout) 
        lbox = nl/double(nlout)
       end else begin
        sbox = ns/float(nsout) 
        lbox = nl/float(nlout)
    end    

; Contraction by an integral amount 

    if (nsout eq long(nsout)) and (nlout eq long(nlout)) then begin
    if ((ns mod nsout) EQ 0) and ((nl mod nlout) EQ 0) then $
                if (dtype EQ 'DOUBLE') or (dtype EQ 'FLOAT') then begin
            if keyword_set(total) then $
           return,rebin(image,nsout,nlout)*sbox*lbox else $
           return,rebin(image,nsout,nlout) 
                endif else begin 
            if keyword_set(total) then $
           return,rebin(float(image),nsout,nlout)*sbox*lbox else $
           return,rebin(float(image),nsout,nlout)
                endelse 


; Expansion by an integral amount
    if ((nsout mod ns) EQ 0) and ((nlout mod nl) EQ 0) then begin
                xindex = long(lindgen(nsout)/(nsout/ns))
                if nl EQ 1 then begin
         if keyword_set(total) then $
        return,interpolate(image,xindex)*sbox else $        
        return,interpolate(image,xindex)  
                endif
                yindex = long(lindgen(nlout)/(nlout/nl))
         if keyword_set(total) then $
        return,interpolate(image,xindex,yindex,/grid)*sbox*lbox else $
        return,interpolate(image,xindex,yindex,/grid)  
    endif
   endif
        ns1 = ns-1
        nl1 = nl-1

; Do 1-d case separately

  if nl EQ 1 then begin
           if dtype eq 'DOUBLE' then result = dblarr(nsout,/NOZERO) $
                    else result = fltarr(nsout,/NOZERO)
        for i=0L,nsout-1 do begin
                rstart = i*sbox           ;starting position for each box
                istart = long(rstart)
                rstop = rstart + sbox      ;ending position for each box
                istop = long(rstop)<ns1
                frac1 = rstart-istart
                frac2 = 1.0 - (rstop-istop)
;
; add pixel values from istart to istop and  subtract fraction pixel 
; from istart to rstart and fraction pixel from rstop to istop
;
                result[i] = total(image[istart:istop]) $
                   - frac1 * image[istart]  $
                   - frac2 * image[istop] 
        endfor
         if keyword_set(total) then return,result $
                      else return,temporary(result)/(sbox*lbox)
 endif 

; Now do 2-d case
; First, bin in second dimension
;
        if dtype eq 'DOUBLE' then temp = dblarr(ns,nlout, /NOZERO) $
                     else temp = fltarr(ns,nlout, /NOZERO)

; loop on output image lines
;
        for i=0L,nlout-1 do begin
                rstart = i*lbox        ;starting position for each box
                istart = long(rstart)
                rstop = rstart + lbox    ;ending position for each box
                istop = long(rstop)<nl1
                frac1 = rstart-istart
                frac2 = 1.0 - (rstop-istop)
;
; add pixel values from istart to istop and  subtract fraction pixel 
; from istart to rstart and fraction pixel from rstop to istop
;

                     if istart EQ istop then $
                  temp[0,i] = (1.0 - frac1 - frac2)*image[*,istart] $
                       else $
                  temp[0,i] = total(image[*,istart:istop],2) $
                   - frac1 * image[*,istart]  $
                   - frac2 * image[*,istop] 
        endfor
           temp = transpose(temp)
;
; bin in first dimension
;
        if dtype eq 'DOUBLE' then result = dblarr(nlout,nsout,/NOZERO) $
                     else result = fltarr(nlout,nsout,/NOZERO)

;
; loop on output image samples
;
        for i=0L,nsout-1 do begin
                rstart = i*sbox           ;starting position for each box
                istart = long(rstart)
                rstop = rstart + sbox      ;ending position for each box
                istop = long(rstop)<ns1
                frac1 = rstart-istart
                frac2 = 1.0 - (rstop-istop)
;
; add pixel values from istart to istop and  subtract fraction pixel 
; from istart to rstart and fraction pixel from rstop to istop
;

            if istart eq istop then $
                        result[0,i] = (1.-frac1-frac2)*temp[*,istart] else $
                result[0,i] = total(temp[*,istart:istop],2)   $
                    - frac1 * temp[*,istart]  $
                    - frac2 * temp[*,istop]
        end

;            
        if keyword_set(total) then $
                        return, transpose(result) $
               else return, transpose(result)/(sbox*lbox)
                      
end;function frebin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PRO ploterror, x, y, xerr, yerr, NOHAT=hat, HATLENGTH=hln, ERRTHICK=eth, $
      ERRSTYLE=est, TYPE=itype, XRANGE = xrange, XLOG=xlog, YLOG=ylog, $
      NSKIP = nskip, NOCLIP = noclip, ERRCOLOR = ecol, YRANGE = yrange, $
      NSUM = nsum, _EXTRA = pkey

;+
; NAME:
;     PLOTERROR
; PURPOSE:
;     Plot data points with accompanying X or Y error bars.
; EXPLANATION:
;     This is a greatly enhanced version of the standard IDL Library routine
;     PLOTERR
;
; CALLING SEQUENCE:
;     ploterror, [ x,]  y, [xerr], yerr [, TYPE=, /NOHAT, HATLENGTH= , NSUM =
;                  ERRTHICK=, ERRSTYLE=, ERRCOLOR=, NSKIP=, .. PLOT keywords]
;
; INPUTS:
;     X = array of abcissae.
;     Y = array of Y values.
;     XERR = array of error bar values (along X)
;     YERR = array of error bar values (along Y)
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;     TYPE = type of plot produced.  The possible types are:
;              TYPE = 0 :       X Linear - Y Linear  (default)
;              TYPE = 1 :       X Linear - Y Log
;              TYPE = 2 :       X Log    - Y Linear
;              TYPE = 3 :       X Log    - Y Log
;              Actually, if 0 is specified, the XLOG and YLOG keywords
;              are used.  If these aren't specified, then a linear-linear
;              plot is produced.  This keyword is available to maintain
;              compatibility with the previous version of PLOTERROR.
;     /NOHAT     = if specified and non-zero, the error bars are drawn
;              without hats.
;     HATLENGTH = the length of the hat lines in device units used to cap the 
;              error bars.   Defaults to !D.X_VSIZE / 100).
;     ERRTHICK  = the thickness of the error bar lines.  Defaults to the
;              THICK plotting keyword.
;     ERRSTYLE  = the line style to use when drawing the error bars.  Uses
;              the same codes as LINESTYLE.
;     ERRCOLOR =  scalar integer (0 - !D.N_TABLE) specifying the color to
;              use for the error bars
;     NSKIP = Integer specifying the error bars to be plotted.   For example,
;              if NSKIP = 2 then every other error bar is plotted; if NSKIP=3
;              then every third error bar is plotted.   Default is to plot
;              every error bar (NSKIP = 1)
;     NSUM =  Number of points to average over before plotting, default=!P.NSUM
;             The errors are also averaged, and then divided by sqrt(NSUM).   
;             This  approximation is meaningful only when the neighboring error
;             bars have similar sizes.    PLOTERROR does not pass the NSUM 
;             keyword to the PLOT command, but rather computes the binning 
;             itself using the  FREBIN function.
;
;     Any valid keywords to the PLOT command (e.g. PSYM, YRANGE) are also 
;     accepted by PLOTERROR via the _EXTRA facility.
;
; RESTRICTIONS:
;       Arrays must not be of type string.  There must be enough points to plot.
;       If only three parameters are input, they will be taken as X, Y and
;       YERR respectively.
;
;       PLOTERROR cannot be used for asymmetric error bars.   Instead use
;       OPLOTERROR with the /LOBAR and /HIBAR keywords.
;
;       Any data points with NAN values in the X, Y, or error vectors are 
;       ignored.
; EXAMPLE:
;       Suppose one has X and Y vectors with associated errors XERR and YERR
;
;       (1) Plot Y vs. X with both X and Y errors and no lines connecting
;           the points
;                  IDL> ploterror, x, y, xerr, yerr, psym=3
;
;       (2) Like (1) but plot only the Y errors bars and omits "hats"
;                  IDL> ploterror, x, y, yerr, psym=3, /NOHAT
;
; WARNING:
;       This an enhanced version of the procedure PLOTERR in the standard IDL
;       distribution.    It was renamed from PLOTERR to PLOTERROR in June 1998
;       in the IDL Astronomy Library to avoid conflict with the RSI procedure.
;
; PROCEDURE:
;       A plot of X versus Y with error bars drawn from Y - YERR to Y + YERR
;       and optionally from X - XERR to X + XERR is written to the output device
;
; PROCEDURE CALLS:
;     FREBIN - used to compute binning if NSUM keyword is present
; MODIFICATION HISTORY:
;     William Thompson        Applied Research Corporation  July, 1986
;     DMS, April, 1989        Modified for Unix
;     Michael R. Greason      ST Systems
;     May, 1991               Added most of the plotting keywords, put hats
;                               on the error bars.
;     K. Venkatakrishna       Added option to plot xerr, May, 1992
;     Michael R. Greason      Corrected handling of reversed axes.  Aug. 1992
;     W. Landsman             Use _EXTRA keyword                    July 1995
;     W. Landsman             Plot more than 32767 points           Feb 1996
;     W. Landsman     Fix Y scaling when only XRANGE supplied       Nov 1996
;     W. Landsman     Added NSKIP keyword                           Dec 1996
;     W. Landsman     Use XLOG, YLOG instead of XTYPE, YTYPE        Jan 1998
;     W. Landsman     Rename to PLOTERROR, OPLOTERROR               Jun 1998
;     W. Landsman     Convert to IDL V5.0                           Jun 1998
;     W. Landsman  Better default scaling when NSKIP supplied       Oct 1998 
;     W. Landsman  Ignore !P.PSYM when drawing error bars           Jan 1999
;     W. Landsman  Handle NSUM keyword correctly                    Aug 1999
;     W. Landsman  Fix case of /XLOG but no X error bars            Oct 1999
;     W. Landsman  Work in the presence of NAN values               Nov 2000
;     W. Landsman  Improve logic when NSUM or !P.NSUM is set        Jan 2001
;     W. Landsman  Only draw error bars with in XRANGE (for speed)  Jan 2002
;     W. Landsman  Fix Jan 2002 update to work with log plots       Jun 2002
;     W. Landsman  Added _STRICT_EXTRA                              Jul 2005
;-
;                       Check the parameters.
 On_error, 2
 compile_opt idl2

 np = N_params()
 IF (np LT 2) THEN BEGIN
        print, "PLOTERROR must be called with at least two parameters."
        print, "Syntax: ploterror, [x,] y, [xerr], yerr"
        RETURN
 ENDIF

; Error bar keywords (except for HATLENGTH; this one will be taken care of 
; later, when it is time to deal with the error bar hats).

 IF (keyword_set(hat)) THEN hat = 0 ELSE hat = 1
 if (n_elements(eth) EQ 0) THEN eth = !P.THICK
 IF (n_elements(est) EQ 0) THEN est = 0
 IF (n_elements(ecol) EQ 0) THEN ecol = !P.COLOR
 if N_elements( NOCLIP ) EQ 0 then noclip = 0
 if not keyword_set(NSKIP) then nskip = 1
 if N_elements(nsum) EQ 0 then nsum = !P.NSUM

;                Other keywords.

 IF (keyword_set(itype)) THEN BEGIN
    CASE (itype) OF
           1 :  ylog = 1    ; X linear, Y log
           2 :  xlog = 1    ; X log, Y linear
           3 :  BEGIN        ; X log, Y log
            xlog = 1
            ylog = 1
            END
        ELSE : 
    ENDCASE
 ENDIF
 if not keyword_set(XLOG) then xlog = 0
 if not keyword_set(YLOG) then ylog = 0
;            If no x array has been supplied, create one.  Make
;            sure the rest of the procedure can know which parameter
;            is which.

 IF np EQ 2 THEN BEGIN            ; Only Y and YERR passed.
    yerr = y
    yy = x
    xx = lindgen(n_elements(yy))
        xerr = make_array(size=size(xx))

 ENDIF ELSE IF np EQ 3 THEN BEGIN     ; X, Y, and YERR passed.
        yerr = xerr
        yy = y
        xx = x

 ENDIF ELSE BEGIN                        ; X, Y, XERR and YERR passed.
    yy = y
        g = where(finite(xerr))
        xerr[g] = abs(xerr[g])
    xx = x
 ENDELSE

 g = where(finite(yerr))               ;Don't take absolute value of NAN values
 yerr[g] = abs(yerr[g])

;            Determine the number of points being plotted.  This
;            is the size of the smallest of the three arrays
;            passed to the procedure.  Truncate any overlong arrays.

 n = N_elements(xx) < N_elements(yy)

 IF np GT 2 then n = n < N_elements(yerr)   
 IF np EQ 4 then n = n < N_elements(xerr)

 IF n LT 2 THEN $
    message,'Not enough points to plot.'

 xx = xx[0:n-1]
 yy = yy[0:n-1]
 yerr = yerr[0:n-1]
 IF np EQ 4 then xerr = xerr[0:n-1]

; If NSUM is greater than one, then we need to smooth ourselves (using FREBIN)

 if nsum GT 1 then begin
      n1 = float(n) / nsum
      n  = long(n1)
      xx = frebin(xx, n1)
      yy = frebin(yy, n1)
      yerror = frebin(yerr,n1)/sqrt(nsum)
      if NP EQ 4 then xerror = frebin(xerr,n1)/sqrt(nsum)
  endif else begin
      yerror = yerr
      if NP EQ 4 then xerror = xerr
  endelse


; If no y-range was passed via keyword or system variable, force one large 
; enough to display all the data and the entire error bars.     
; If a reversed y-range was passed, switch ylo and yhi.

 ylo = yy - yerror
 yhi = yy + yerror

 if not keyword_set( YRANGE ) then yrange = !Y.RANGE
 IF yrange[0] EQ yrange[1] THEN BEGIN
    if keyword_set( XRANGE ) then  begin
        good = where( (xx GT min(xrange)) and (xx LT max(xrange)) )
        yrange = [min(ylo[good],/NAN), max(yhi[good], /NAN)]
    endif else yrange = [min(ylo,/NAN), max(yhi, /NAN)]
 ENDIF ELSE IF yrange[0] GT yrange[1] THEN BEGIN
    ylo = yy + yerror
    yhi = yy - yerror
 ENDIF

;        Similarly for x-range

 if not keyword_set( XRANGE ) then xrange = !X.RANGE
 if NP EQ 4 then begin
   xlo = xx - xerror
   xhi = xx + xerror
   IF xrange[0] EQ xrange[1] THEN xrange = [min(xlo,/NAN), max(xhi,/NAN)]
   IF xrange[0] GT xrange[1] THEN BEGIN
      xlo = xx + xerror
      xhi = xx - xerror
   ENDIF
 endif

; Plot the positions.    Always set NSUM = 1 since we already took care of 
; smoothing with FREBIN

 plot, xx, yy, XRANGE = xrange, YRANGE = yrange, XLOG = xlog, YLOG = ylog, $
         _STRICT_EXTRA = pkey, NOCLIP = noclip, NSum= 1

;    Plot the error bars.   Compute the hat length in device coordinates
;       so that it remains fixed even when doing logarithmic plots.

    data_low = convert_coord(xx,ylo,/TO_DEVICE)
    data_hi = convert_coord(xx,yhi,/TO_DEVICE)
    if NP EQ 4 then begin
       x_low = convert_coord(xlo,yy,/TO_DEVICE)
       x_hi = convert_coord(xhi,yy,/TO_DEVICE)
    endif
    ycrange = !Y.CRANGE   &  xcrange = !X.CRANGE
    sv_psym = !P.PSYM & !P.PSYM = 0
    if ylog EQ 1 then ylo = ylo > 10^ycrange[0]
    if (xlog EQ 1) and (np EQ 4) then xlo = xlo > 10^xcrange[0]
; Only draw error bars for X values within XCRANGE
    if xlog EQ 1 then xcrange = 10^xcrange
    g = where((xx GT xcrange[0]) and (xx LE xcrange[1]), Ng)

    if (Ng GT 0) and (Ng NE n) then begin  
          istart = min(g, max = iend)  
    endif else begin
          istart = 0L & iend = n-1
    endelse

 FOR i = istart, iend, Nskip DO BEGIN     

    plots, [xx[i],xx[i]], [ylo[i],yhi[i]], LINESTYLE=est,THICK=eth,  $
        NOCLIP = noclip, COLOR = ecol
;                                                         Plot X-error bars 
    if np EQ 4 then plots, [xlo[i],xhi[i]],[yy[i],yy[i]],LINESTYLE=est, $
        THICK=eth, COLOR = ecol, NOCLIP = noclip
    IF (hat NE 0) THEN BEGIN
        IF (N_elements(hln) EQ 0) THEN hln = !D.X_VSIZE/100. 
        exx1 = data_low[0,i] - hln/2.
        exx2 = exx1 + hln
        plots, [exx1,exx2], [data_low[1,i],data_low[1,i]],COLOR=ecol, $
                      LINESTYLE=est,THICK=eth,/DEVICE, noclip = noclip
        plots, [exx1,exx2], [data_hi[1,i],data_hi[1,i]], COLOR = ecol,$
                       LINESTYLE=est,THICK=eth,/DEVICE, noclip = noclip

;                                                        Plot Y-error bars

                IF np EQ 4 THEN BEGIN
                   IF (N_elements(hln) EQ 0) THEN hln = !D.Y_VSIZE/100.
                   eyy1 = x_low[1,i] - hln/2.
                   eyy2 = eyy1 + hln
                   plots, [x_low[0,i],x_low[0,i]], [eyy1,eyy2],COLOR = ecol, $
                         LINESTYLE=est,THICK=eth,/DEVICE, NOCLIP = noclip
                   plots, [x_hi[0,i],x_hi[0,i]], [eyy1,eyy2],COLOR = ecol, $
                         LINESTYLE=est,THICK=eth,/DEVICE, NOCLIP = noclip
                ENDIF
    ENDIF
    NOPLOT:
 ENDFOR
 !P.PSYM = sv_psym
;
 RETURN
END;PRO ploterror
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
PRO  oploterror, x, y, xerr, yerr, NOHAT=hat, HATLENGTH=hln, ERRTHICK=eth, $
      ERRSTYLE=est, THICK = thick, NOCLIP=noclip, ERRCOLOR = ecol, Nsum = nsum,$
      NSKIP=nskip, LOBAR=lobar, HIBAR=hibar,_EXTRA = pkey
;+
; NAME:
;      OPLOTERROR
; PURPOSE:
;      Over-plot data points with accompanying X or Y error bars.
; EXPLANATION:
;      For use instead of PLOTERROR when the plotting system has already been
;      defined. 
;
; CALLING SEQUENCE:
;      oploterror, [ x,]  y, [xerr], yerr,   
;            [ /NOHAT, HATLENGTH= , ERRTHICK =, ERRSTYLE=, ERRCOLOR =, 
;              /LOBAR, /HIBAR, NSKIP = , NSUM = , ... OPLOT keywords ]
; INPUTS:
;      X = array of abcissae, any datatype except string
;      Y = array of Y values, any datatype except string
;      XERR = array of error bar values (along X)
;      YERR = array of error bar values (along Y)
;
; OPTIONAL INPUT KEYWORD PARAMETERS:
;      /NOHAT     = if specified and non-zero, the error bars are drawn
;                  without hats.
;      HATLENGTH = the length of the hat lines used to cap the error bars.
;                  Defaults to !D.X_VSIZE / 100).
;      ERRTHICK  = the thickness of the error bar lines.  Defaults to the
;                  THICK plotting keyword.
;      ERRSTYLE  = the line style to use when drawing the error bars.  Uses
;                  the same codes as LINESTYLE.
;      ERRCOLOR =  scalar integer (0 - !D.N_TABLE) specifying the color to
;                  use for the error bars
;      NSKIP = Positive Integer specifying the error bars to be plotted.   
;            For example, if NSKIP = 2 then every other error bar is 
;            plotted; if NSKIP=3 then every third error bar is plotted.   
;            Default is to plot every error bar (NSKIP = 1)
;      NSUM =  Number of points to average over before plotting, default = 
;             !P.NSUM  The errors are also averaged, and then divided by 
;             sqrt(NSUM).   This approximation is meaningful only when the 
;             neighboring error bars have similar sizes.
; 
;      /LOBAR = if specified and non-zero, will draw only the -ERR error bars.
;      /HIBAR = if specified and non-zero, will draw only the +ERR error bars.
;                  If neither LOBAR or HIBAR are set _or_ if both are set,
;                  you will get both error bars.  Just specify one if you
;                  only want one set.
;     Any valid keywords to the OPLOT command (e.g. PSYM, YRANGE) are also 
;     accepted by OPLOTERROR via the _EXTRA facility.
;
; NOTES:
;     If only two parameters are input, they are taken as Y and YERR.  If only
;     three parameters are input, they will be taken as X, Y and YERR, 
;     respectively.
;
; EXAMPLE:
;      Suppose one has X and Y vectors with associated errors XERR and YERR
;      and that a plotting system has already been defined:
;
;       (1) Overplot Y vs. X with both X and Y errors and no lines connecting
;           the points
;                  IDL> oploterror, x, y, xerr, yerr, psym=3
;
;       (2) Like (1) but overplot only the Y errors bars and omits "hats"
;                  IDL> oploterror, x, y, yerr, psym=3, /NOHAT
;
;       (3) Like (2) but suppose one has a positive error vector YERR1, and 
;               a negative error vector YERR2 (asymmetric error bars)
;                  IDL> oploterror, x, y, yerr1, psym=3, /NOHAT,/HIBAR
;                  IDL> oploterror, x, y, yerr2, psym=3, /NOHAT,/LOBAR
;
; PROCEDURE:
;      A plot of X versus Y with error bars drawn from Y - YERR to Y + YERR
;      and optionally from X - XERR to X + XERR is written to the output device
;
; WARNING:
;      This an enhanced version of the procedure OPLOTERR in the standard RSI
;      library.    It was renamed to OPLOTERROR in June 1998 in the IDL 
;      Astronomy library.
;
; MODIFICATION HISTORY:
;      Adapted from the most recent version of PLOTERR.  M. R. Greason,
;            Hughes STX, 11 August 1992.
;      Added COLOR keyword option to error bars W. Landsman   November 1993
;      Add ERRCOLOR, use _EXTRA keyword,           W. Landsman, July 1995
;      Remove spurious call to PLOT_KEYWORDS     W. Landsman, August 1995
;      OPLOT more than 32767 error bars          W. Landsman, Feb 1996
;      Added NSKIP keyword                       W. Landsman, Dec 1996
;      Added HIBAR and LOBAR keywords, M. Buie, Lowell Obs., Feb 1998
;      Rename to OPLOTERROR    W. Landsman    June 1998
;      Converted to IDL V5.0   W. Landsman    June 1998
;      Ignore !P.PSYM when drawing error bars   W. Landsman   Jan 1999
;      Handle NSUM keyword correctly           W. Landsman    Aug 1999
;      Check limits for logarithmic axes       W. Landsman    Nov. 1999
;      Work in the presence of  NAN values     W. Landsman    Dec 2000
;      Improve logic when NSUM or !P.NSUM is set  W. Landsman      Jan 2001
;      Remove NSUM keyword from PLOTS call    W. Landsman      March 2001
;      Only draw error bars with in XRANGE (for speed)  W. Landsman Jan 2002
;      Fix Jan 2002 update to work with log plots  W. Landsman Jun 2002
;      Added STRICT_EXTRA keyword   W. Landsman     July 2005
;-
;                  Check the parameters.
;
 On_error, 2
 compile_opt idl2
 np = N_params()
 IF (np LT 2) THEN BEGIN
      print, "OPLOTERR must be called with at least two parameters."
      print, "Syntax: oploterr, [x,] y, [xerr], yerr, [..oplot keywords... "
      print,'     /NOHAT, HATLENGTH = , ERRTHICK=, ERRSTLYE=, ERRCOLOR='
      print,'     /LOBAR, /HIBAR, NSKIP= ]'
      RETURN
 ENDIF

; Error bar keywords (except for HATLENGTH; this one will be taken care of 
; later, when it is time to deal with the error bar hats).

 IF (keyword_set(hat)) THEN hat = 0 ELSE hat = 1
 if not keyword_set(THICK) then thick = !P.THICK
 IF (n_elements(eth) EQ 0) THEN eth = thick
 IF (n_elements(est) EQ 0) THEN est = 0
 IF (n_elements(ecol) EQ 0) THEN ecol = !P.COLOR
 if N_elements( NOCLIP ) EQ 0 THEN noclip = 0
 if not keyword_set(NSKIP) then nskip = 1
 if N_elements(nsum) EQ 0 then nsum = !P.NSUM
 if not keyword_set(lobar) and not keyword_set(hibar) then begin
      lobar=1
      hibar=1
 endif else if keyword_set(lobar) and keyword_set(hibar) then begin
      lobar=1
      hibar=1
 endif else if keyword_set(lobar) then begin
      lobar=1
      hibar=0
 endif else begin
      lobar=0
      hibar=1
 endelse
;
; If no X array has been supplied, create one.  Make sure the rest of the 
; procedure can know which parameter is which.
;
 IF np EQ 2 THEN BEGIN                  ; Only Y and YERR passed.
      yerr = y
      yy = x
      xx = indgen(n_elements(yy))
      xerr = make_array(size=size(xx))

 ENDIF ELSE IF np EQ 3 THEN BEGIN       ; X, Y, and YERR passed.
        yerr = xerr
        yy = y
        xx = x

 ENDIF ELSE BEGIN                        ; X, Y, XERR and YERR passed.
      yy = y
      g = where(finite(xerr))
      xerr[g] = abs(xerr[g])
      xx = x
 ENDELSE

 g = where(finite(yerr))
 yerr[g] = abs(yerr[g])

;
;                  Determine the number of points being plotted.  This
;                  is the size of the smallest of the three arrays
;                  passed to the procedure.  Truncate any overlong arrays.
;

 n = N_elements(xx) < N_elements(yy)

 IF np GT 2 then n = n < N_elements(yerr)   
 IF np EQ 4 then n = n < N_elements(xerr)

 xx = xx[0:n-1]
 yy = yy[0:n-1]
 yerr = yerr[0:n-1]
 IF np EQ 4 then xerr = xerr[0:n-1]

; If NSUM is greater than one, then we need to smooth ourselves (using FREBIN)

 if NSum GT 1 then begin
      n1 = float(n) / nsum
      n  = long(n1)
      xx = frebin(xx, n1)
      yy = frebin(yy, n1)
      yerror = frebin(yerr,n1)/sqrt(nsum)
      if NP EQ 4 then xerror = frebin(xerr,n1)/sqrt(nsum)
  endif else begin
      yerror = yerr
      if NP EQ 4 then xerror = xerr
  endelse

 ylo = yy - yerror*lobar
 yhi = yy + yerror*hibar

 if Np EQ 4 then begin
     xlo = xx - xerror*lobar
     xhi = xx + xerror*hibar
 endif
;
;                  Plot the positions.
;
 if n NE 1 then begin
     oplot, xx, yy, NOCLIP=noclip,THICK = thick,_STRICT_EXTRA = pkey 
 endif else begin 
     plots, xx, yy, NOCLIP=noclip,THICK = thick,_STRICT_EXTRA = pkey
 endelse
;
; Plot the error bars.   Compute the hat length in device coordinates
; so that it remains fixed even when doing logarithmic plots.
;
 data_low = convert_coord(xx,ylo,/TO_DEVICE)
 data_hi = convert_coord(xx,yhi,/TO_DEVICE)
 if NP EQ 4 then begin
    x_low = convert_coord(xlo,yy,/TO_DEVICE)
    x_hi = convert_coord(xhi,yy,/TO_DEVICE)
 endif
 ycrange = !Y.CRANGE   &  xcrange = !X.CRANGE
    if !Y.type EQ 1 then ylo = ylo > 10^ycrange[0]
    if (!X.type EQ 1) and (np EQ 4) then xlo = xlo > 10^xcrange[0]
 sv_psym = !P.PSYM & !P.PSYM = 0     ;Turn off !P.PSYM for error bars
; Only draw error bars for X values within XCRANGE
    if !X.TYPE EQ 1 then xcrange = 10^xcrange
    g = where((xx GT xcrange[0]) and (xx LE xcrange[1]), Ng)
    if (Ng GT 0) and (Ng NE n) then begin  
          istart = min(g, max = iend)  
    endif else begin
          istart = 0L & iend = n-1
    endelse
    
 FOR i = istart, iend, Nskip DO BEGIN

    plots, [xx[i],xx[i]], [ylo[i],yhi[i]], LINESTYLE=est,THICK=eth,  $
           NOCLIP = noclip, COLOR = ecol

    ; Plot X-error bars 
    ;
    if np EQ 4 then $
       plots, [xlo[i],xhi[i]],[yy[i],yy[i]],LINESTYLE=est, $
              THICK=eth, COLOR = ecol, NOCLIP = noclip

    IF (hat NE 0) THEN BEGIN
       IF (N_elements(hln) EQ 0) THEN hln = !D.X_VSIZE/100. 
       exx1 = data_low[0,i] - hln/2.
       exx2 = exx1 + hln
       if lobar then $
          plots, [exx1,exx2], [data_low[1,i],data_low[1,i]],COLOR=ecol, $
                 LINESTYLE=est,THICK=eth,/DEVICE, noclip = noclip
       if hibar then $
          plots, [exx1,exx2], [data_hi[1,i],data_hi[1,i]], COLOR = ecol,$
                 LINESTYLE=est,THICK=eth,/DEVICE, noclip = noclip
;                                          
       IF np EQ 4 THEN BEGIN
          IF (N_elements(hln) EQ 0) THEN hln = !D.Y_VSIZE/100.
             eyy1 = x_low[1,i] - hln/2.
             eyy2 = eyy1 + hln
             if lobar then $
                plots, [x_low[0,i],x_low[0,i]], [eyy1,eyy2],COLOR = ecol, $
                       LINESTYLE=est,THICK=eth,/DEVICE, NOCLIP = noclip
             if hibar then $
                plots, [x_hi[0,i],x_hi[0,i]], [eyy1,eyy2],COLOR = ecol, $
                       LINESTYLE=est,THICK=eth,/DEVICE, NOCLIP = noclip
          ENDIF
       ENDIF
    NOPLOT:
ENDFOR
 !P.PSYM = sv_psym 
;
RETURN
END; PRO  oploterror
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
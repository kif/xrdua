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

function HXIE_inittransfo,shiftx=shiftx,shifty=shifty,angle=angle,scale=scale,X0=X0,Y0=Y0,size=size
if n_elements(shiftx) ne 1 then shiftx=0.
if n_elements(shifty) ne 1 then shifty=0.
if n_elements(angle) ne 1 then angle=0.
if n_elements(scale) ne 1 then scale=1.
if n_elements(size) eq 2 then begin
    ; center of image with dimension s[0] x s[1]
    X0=(size[0]-1)/2.
    Y0=(size[1]-1)/2.
endif
if n_elements(X0) ne 1 then X0=0.
if n_elements(Y0) ne 1 then Y0=0.
return,{shiftx:shiftx,shifty:shifty,angle:angle,scale:scale,X0:X0,Y0:Y0}
end;function HXIE_inittransfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro HXie_encrect,transfo,sz,dimx,dimy,XL,YL
; Concider the transformation of an image as the transformation of
; a rectangle. This procedure calculates the enclosing box of
; this transformed rectangle. The enclosing box is defined by
; its left-bottom coordinate in the image coordinate system
; and its size.

; Input: 
; transfo: transformation info
; sz: size of transformed rectangle
;
; Output:
; dimx,dimy: size of enclosing box
; XL,YL: left-bottom corner of enclosing box

; Forward transformation of corners (cfr. HXie_transformation)
theta = (transfo.angle mod 360) * !dpi/180
c = cos(theta)*transfo.scale
s = sin(theta)*transfo.scale
XC = transfo.X0+transfo.shiftx
YC = transfo.Y0+transfo.shifty
P = [XC-c*transfo.X0+s*transfo.Y0,-s,c,0]
Q = [YC-s*transfo.X0-c*transfo.Y0,c,s,0]
x = [0,0,sz[0],sz[0]]-0.5  ; boarder of outer pixels
y = [0,sz[1],0,sz[1]]-0.5 ; boarder of outer pixels
xx = P[0]+P[1]*y+P[2]*x
yy = Q[0]+Q[1]*y+Q[2]*x
    
; Enclosing box
XL=round(min(xx,max=XR))
XR=round(XR)
dimx=XR-XL+1
YL=round(min(yy,max=YR))
YR=round(YR)
dimy=YR-YL+1

end;pro HXie_encrect
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXIE_inversetransfo,transfo,sz=sz
; Get the inverse transformation (cfr. HXie_transformation)

; X' = M3(0,1,-XL,-YL).M2(angle,scale,X0+shiftx,Y0+shifty).M1(0,1,-X0,-Y0).X       (M.X=A.X+T)
;    = A2.(X+T1)+T2+T3
; X  = A2^(-1).(X'-T2-T3)-T1
;     = A2'.(X+T1')+T2'+T3'
;     
; A2' = A2^(-1)        => angle'=-angle, scale'=1/scale
; -T1' = T2+T3        => [X0',Y0'] = [X0+shiftx-XL,Y0+shifty-YL]
; T2' + T3' = -T1    => [X0'+shiftx'-XL',Y0'+shifty'-YL'] = [X0,Y0]
;                    <=> [X0+shiftx-XL+shiftx'-XL',Y0+shifty-YL+shifty'-YL'] = [X0,Y0]
;                    => shiftx' = -shiftx+XL+XL'
;                       shifty' = -shifty+YL+YL'

if n_elements(sz) ne 2 then begin
    XL=0
    YL=0
endif else HXie_encrect,transfo,sz,dimx,dimy,XL,YL

X0 = transfo.X0+transfo.shiftx-XL
Y0 = transfo.Y0+transfo.shifty-YL

return,HXIE_inittransfo(shiftx=-transfo.shiftx,shifty=-transfo.shifty,$
        angle=-transfo.angle,scale=1./transfo.scale,X0=X0,Y0=Y0)

end;function HXIE_inversetransfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXIE_combinetransfo,transfo1,transfo2,sz=sz
; Combine two transformations (cfr. HXie_transformation)

; transformation1: 
; X' = M3(0,1,-XL,-YL).M2(angle,scale,X0+shiftx,Y0+shifty).M1(0,1,-X0,-Y0).X       (M.X=A.X+T)
;    = A2.(X+T1)+T2+T3
;    
; transformation2:
; X''= M3'(0,1,-XL',-YL').M2'(angle',scale',X0'+shiftx',Y0'+shifty').M1'(0,1,-X0',-Y0').X'
;    = A2'.(X'+T1')+T2'+T3'
;     = A2'.(A2.(X+T1)+T2+T3+T1')+T2'+T3'
;     = A2'.A2.(X+T1)+A2'.(T2+T3+T1')+T2'+T3'
;     = A2''.(X+T1'')+T2''+T3''
;
; T1'' = T1                        => [X0'',Y0''] = [X0,Y0]
; A2'' = A2'.A2                    => angle''=angle+angle', scale''=scale'.scale
; T2'' = A2'.(T2+T3+T1')+T2'     => [X0''+shiftx'',Y0''+shifty''] = [X0+shiftx'',Y0+shifty''] = [shiftx'',shifty'']-T1 = A2'(T2+T3+T1')+T2'
;                                <=> [shiftx'',shifty''] = A2'.(T2+T3+T1')+T2'+T1
;                                <=>                     = A2'.[X0+shiftx-XL-X0',Y0+shifty-YL-Y0'] + [X0'+shiftx'-X0,Y0'+shifty'-Y0]
; T3'' = T3'                     => XL'' = XL'

theta = (transfo2.angle mod 360) * !dpi/180
c = cos(theta)*transfo2.scale
s = sin(theta)*transfo2.scale

if n_elements(sz) ne 2 then begin
    XL=0
    YL=0
endif else HXie_encrect,transfo1,sz,dimx,dimy,XL,YL
XC = transfo1.X0+transfo1.shiftx-XL-transfo2.X0
YC = transfo1.Y0+transfo1.shifty-YL-transfo2.Y0
XT = transfo2.X0+transfo2.shiftx-transfo1.X0
YT = transfo2.Y0+transfo2.shifty-transfo1.Y0
shift=[[c,-s],[s,c]]##[[XC],[YC]]+[[XT],[YT]]

return,HXIE_inittransfo(shiftx=shift[0],shifty=shift[1],$
        angle=transfo1.angle+transfo2.angle,scale=transfo1.scale*transfo2.scale,X0=transfo1.X0,Y0=transfo1.Y0)
        
end;function HXIE_combinetransfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXie_transformation,source,transfo,resize=resize, MISSING = MISSING, INTERP = INTERP, $
                XC=XC,YC=YC,XL=XL,YL=YL

; source: 2D array
; transfo: rotate, scale and translate information
;              shiftx: pixels in source
;              shifty: pixels in source
;              angle: rotation in degrees around [X0,Y0]
;              scale: scale in %
;              X0,Y0: center of rotation is source
;              
; XC,YC: center of rotation in output image
; resize: make sure the transformed source is completely visible (output not the same size as source)
;
; Rotate, scale and translate:
;    [x',y',1]^T = [[scale.cos(angle),-scale.sin(angle),shiftx],[scale.sin(angle),scale.cos(angle),shifty],[0,0,1]].[x,y,1]^T
;                = [[c,-s,shiftx],[s,c,shifty],[0,0,1]].[x,y,1]^T
;                = M(angle,scale,shiftx,shifty).[x,y,1]^T
;
; Forward transformation:
; origin = center of the left-bottom corner pixel in image, i.e. the first pixel in the matrix
; 1. Shift center of rotation (X0,Y0) to origin: M1(0,1,-X0,-Y0)
; 2. Scale, rotate and shift origin to (XC,YC)=(X0+shiftx,Y0+shifty): M2(angle,scale,XC,YC)
; 3. Optional: Shift (XL,YL) to origin: M3(0,1,-XL,-YL)  (this can be included in XC,YC)
;    =>  X' = (M3.)M2.M1.X
;    <=> x' = XC-cX0+sY0 - sy + cx       (c=scale.cos(angle), s=scale.sin(angle))
;        y' = YC-sX0-cY0 + cy + sx
;        
; Remark: XL and YL are the left-bottom coordinates of the enclosing box of the image
;         after applying M1 and M2
;
; Backward transformation:
; 1. Optional: Shift origin to (XL,YL): M4(0,1,XL,YL)  (this can be included in XC,YC)
; 2. Shift (XC,YC) to origin: M5(0,1,-XC,-YC)
; 3. Scale, rotate and shift origin to (X0,Y0): M6(-angle,1/scale,X0,Y0)
;    =>  X = M6.M5.(M4.)X'
;    <=> x = X0-cXC+sYC - sy' + cx'       (c=cos(-angle)/scale, s=sin(-angle)/scale)
;        y = Y0-sXC-cYC + cy' + sx'
;
; Backward (polynomial) warping:
;     f: X -> Y: Y = f(X)
;     M: X -> X': X' = M(X)
;     g: X' -> Y: Y = g(X') = g(M(X))
;     =>  g(X') = f(M^(-1)(X'))
;     <=> dest[x',y'] = source[x,y]
;                     
;             X' -----g(X')-----> Y
;             |                 /
;             |               /
;             |             /
;            M(X)         /
;             |        f(X) = g(M(X))
;             |       /
;             |     /
;             |   /
;             | /
;             X
;     
;     M(X): 2D polynomial
;     => g(x',y') = f(x,y) = f(a(x',y'),b(x',y'))
;     polynomials in x and y of degree N:
;    a(x,y) = sum_ij P_ij.x^j.y^i
;    b(x,y) = sum_ij Q_ij.x^j.y^i
;    for degree 1:
;    a(x,y) = P_00 + P_10.y + P_01.x + P_11.x.y
;    b(x,y) = Q_00 + Q_10.y + Q_01.x + Q_11.x.y
;
; Difference with IDL's ROT:
;    1. rotation in other direction
;    2. adds translation
;    3. resize possiblity

; center of rotation
X0 = transfo.X0
Y0 = transfo.Y0

; new position of center of rotation (i.e. after transformation)
if keyword_set(resize) then begin
    HXie_encrect,transfo,size(source,/dim),dimx,dimy,XL,YL
    XC = X0+transfo.shiftx-XL
    YC = Y0+transfo.shifty-YL
endif else begin
    XC = X0+transfo.shiftx
    YC = Y0+transfo.shifty
    XL = 0L
    YL = 0L
endelse

; Coefficients for backward polynomial warping
theta = (-transfo.angle mod 360) * !dpi/180
c = cos(theta)/transfo.scale
s = sin(theta)/transfo.scale
P = [X0-c*XC+s*YC,-s,c,0]
Q = [Y0-s*XC-c*YC,c,s,0]

if n_elements(MISSING) eq 0 then MISSING=0
if n_elements(INTERP) eq 0 then INTERP=2
if INTERP eq 2 then CUBIC=-0.5
if keyword_set(resize) then $
    return,poly_2d(source, P, Q, INTERP, dimx, dimy, MISSING = MISSING, CUBIC=CUBIC) $ ; cubic interpolation
else $
    return,poly_2d(source, P, Q, INTERP, MISSING = MISSING, CUBIC=CUBIC)
end;function HXie_transformation
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro HXie_transform,source,transfo,error=error,_ref_extra=extra
; HXie_transformation for RGB or grayscale images
error=1b
s=size(source,/dim)
case n_elements(s) of
2:  begin
    source=HXie_transformation(source,transfo,_EXTRA=extra)
    endcase
3:    begin
    ind=where(s eq 3,ct)
    if ct ne 1 then return
    case ind[0] of
    0:    begin
        im1=HXie_transformation(reform(source[0,*,*]),transfo,_EXTRA=extra)
        im2=HXie_transformation(reform(source[1,*,*]),transfo,_EXTRA=extra)
        im3=HXie_transformation(reform(source[2,*,*]),transfo,_EXTRA=extra)
        source=transpose([[[temporary(im1)]],[[temporary(im2)]],[[temporary(im3)]]],[2,0,1])
        endcase
    1:    begin
        im1=HXie_transformation(reform(source[*,0,*]),transfo,_EXTRA=extra)
        im2=HXie_transformation(reform(source[*,1,*]),transfo,_EXTRA=extra)
        im3=HXie_transformation(reform(source[*,2,*]),transfo,_EXTRA=extra)
        source=transpose([[[temporary(im1)]],[[temporary(im2)]],[[temporary(im3)]]],[0,2,1])
        endcase
    2:    begin
        im1=HXie_transformation(reform(source[*,*,0]),transfo,_EXTRA=extra)
        im2=HXie_transformation(reform(source[*,*,1]),transfo,_EXTRA=extra)
        im3=HXie_transformation(reform(source[*,*,2]),transfo,_EXTRA=extra)
        source=[[[temporary(im1)]],[[temporary(im2)]],[[temporary(im3)]]]
        endcase
    else: return
    endcase
    
    endcase
else: return
endcase

error=0b
end;pro HXie_transform
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXie_grayscale,img,error=error
; Make grayscale from RGB image

error=1b

s=size(img,/dim)
case n_elements(s) of
2:  begin
    error=0b
    return,img
    endcase
3:    begin
    ind=where(s eq 3,ct)
    if ct ne 1 then return,0
    case ind[0] of
    0: img2=reform(0.299 * img[0,*,*] + 0.587 * img[1,*,*] + 0.114 * img[2,*,*])
    1: img2=reform(0.299 * img[*,0,*] + 0.587 * img[*,1,*] + 0.114 * img[*,2,*])
    2: img2=0.299 * img[*,*,0] + 0.587 * img[*,*,1] + 0.114 * img[*,*,2]
    endcase
    error=0b
    return,img2
    endcase
else: return,0
endcase

end;function HXie_grayscale
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXie_zerof,s
; Position of [u,v]=[0,0] in the Fourier domain of an s[0] x s[1] array
return,[s[0]/2-1+(s[0] mod 2),s[1]/2-1+(s[1] mod 2)]
end;function HXie_zerof
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro HXie_HighPassCenter,ARR,CEN=CEN
; ARR in frequency domain
s=size(ARR,/dim)

; u,v in frequency domain
piu=rebin((lindgen(s[0])-(s[0]-2+(s[0] mod 2))/2.)*(!pi/s[0]),s,/sample)
piv=rebin((lindgen(1,s[1])-(s[1]-2+(s[1] mod 2))/2.)*(!pi/s[1]),s,/sample)
h=cos(piu)*cos(piv)

; Sort ARR so that it corresponds with u and v (i.e. center FFT)
CEN=HXie_zerof(s)
ARR=shift(ARR,CEN[0],CEN[1])

; Filter
ARR*=(1-h)*(2-h)
end;pro HXie_HighPassCenter
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXie_LogPolarBase,s,CEN
; base so that the whole image-range in transformed to log-polar coordinates (i.e. r not too small and not too large)
rmax=s-CEN-1
rmax[0]>=CEN[0]
rmax=sqrt(rmax[0]^2.+rmax[1]^2.)
return,10^(alog10(rmax)/(s[0]-1))
end;function HXie_LogPolarBase
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro HXie_LogPolarTrans, RECT, CEN
; RECT: Fourier domain with [u,v]=[0,0] at pixel CEN
; 
; Convert from Cartesian coordinates to log-polar
;    Cartesian to log-polar (forward warping)
;        R = logb(r)   r=sqrt(x^2+y^2)
;        t = atan(y/x)
;    log-polar to Cartesian (backward warping)
;        x = b^R cos(t)
;        y = b^R sin(t)
; 
; Theta: only upper 2 quadrants (lower 2 are a copy because |F(u,v)|=|F(-u,-v)|)
; t range of log-polar plot (must be equal spaced): [0,!pi] = [0,...,s[1]-1]*!pi/(s[1]-1)) 
; R range of log-polar plot (must be equal spaced): [0,...,s[0]-1]
; r range of log-polar plot: b^[0,...,s[0]-1] 
; 
; Choose b so that:     b^(s[0]-1) = rmax
;                    <=> b = 10 ^ (alog10(rmax)/(s[0]-1))

s=size(RECT,/dim)
dTheta=!pi/(s[1]-1)
b=HXie_LogPolarBase(s,CEN)

r=rebin(b^lindgen(s[0]),s,/sample)
theta=rebin(dTheta*lindgen(1,s[1]),s,/sample)
x=r*cos(theta)+CEN[0]
y=temporary(r)*sin(temporary(theta))+CEN[1]

RECT=interpolate(RECT,x,y,missing=0) ; bilinear interpolation
end;pro HXie_LogPolarTrans
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro HXie_ScaleRot2shift, Im, CEN=CEN
; moduli of (Forward) Fourier transform
Im=FFT(Im,-1,/overwrite)
Im=abs(Im)

; Highpass filter + move [u,v]=[0,0] to CEN (~middle of the image)
HXie_HighPassCenter, Im, CEN=CEN

; log-polar coordinates
HXie_LogPolarTrans, Im, CEN
end;pro HXie_ScaleRot2shift
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXie_FourierFreq,n
return,(lindgen(n)-(n-2+(n mod 2))/2.)/n
end;function HXie_FourierFreq
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXie_FourierShift,Im1,Im2,logpolar=logpolar,$
                        maxshiftxrange=maxshiftxrange,maxshiftyrange=maxshiftyrange,debug=debug

; Fourier transform
Im1=FFT(Im1,-1,/overwrite)
Im2=FFT(Im2,-1,/overwrite)

; Cross-power spectrum of Im1 and Im2
R=Im1*conj(Im2)/(abs(Im1)*abs(Im2))

; Inverse cross-power spectrum
R = FFT(R, 1, /overwrite)
R = abs(R)
; Center [u,v]=[0,0]
s=size(R,/dim)
CEN=HXie_zerof(s)
R=shift(R,CEN[0],CEN[1])

; Get maximal position
maxshiftxrange=[0,s[0]-1]-CEN[0]
maxshiftyrange=[0,s[1]-1]-CEN[1]
maxn = MAX(R, position,/nan)
x1=(position[0] MOD s[0]) - CEN[0]
y1=(position[0] / s[0]) - CEN[1]

if keyword_set(logpolar) and keyword_set(debug) then begin
    window,3 & tvscl,r<(max(r)*0.1)
    plots,CEN[[0,0]],[0,s[1]-1],/device
    plots,[0,s[0]-1],CEN[[1,1]],/device
endif

return,[x1,y1]

; Get shift with subpixel precision: (C. G. Schiek 2007)

; Highest surrounding pixel (only diagonal)
x2=x1+[1,1,-1,-1]
y2=y1+[1,-1,1,-1]
maxn = MAX(R[x2,y2], position,/nan)
x2=x2[position[0]]
y2=y2[position[0]]

if keyword_set(logpolar) then alpha=1.55 else alpha=0.65
w1=R[x1,y1]^alpha+R[x1,y2]^alpha
w2=R[x2,y1]^alpha+R[x2,y2]^alpha
v1=R[x1,y1]^alpha+R[x2,y1]^alpha
v2=R[x1,y2]^alpha+R[x2,y2]^alpha

x=(w1*x1+w2*x2)/(w1+w2)
y=(v1*y1+v2*y2)/(v1+v2)
return,[x,y]

end;function HXie_FourierShift
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXie_equalsize,Im1,Im2,transfo=transfo
; Make pictures of the same size and return the transformation Im2 -> Im1

s1=size(Im1,/dim)
s2=size(Im2,/dim)
if array_equal(s1,s2) then return,0

if n_elements(Im1) gt n_elements(Im2) then begin
    ; Make Im1 smaller
    scale=s2/float(s1)
    scale=scale[0]>scale[1]
    
    transfo=HXIE_inittransfo(scale=scale)
    Im1=HXie_transformation(Im1,transfo)
    
    s=floor(s1*scale)<s2
    Im1=Im1[0:s[0]-1,0:s[1]-1]
    Im2=Im2[0:s[0]-1,0:s[1]-1]
    
    transfo=HXIE_inittransfo(scale=1./scale)
    return,-1
endif else begin
    ; Make Im2 smaller
    scale=s1/float(s2)
    scale=scale[0]>scale[1]
    
    transfo=HXIE_inittransfo(scale=scale)
    Im2=HXie_transformation(Im2,transfo)
    
    s=floor(s2*scale)<s1
    Im1=Im1[0:s[0]-1,0:s[1]-1]
    Im2=Im2[0:s[0]-1,0:s[1]-1]
    return,1
endelse

end;function HXie_equalsize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXie_register,source,dest,error=error,maxscalerange=maxscalerange,maxanglerange=maxanglerange,$
                                    maxshiftxrange=maxshiftxrange,maxshiftyrange=maxshiftyrange,debug=debug
;---------------------------------------------------------------------------------
;  Based on:
;  
;  Authors:
;  Hongjie Xie,  Nigel Hicks,  G. Randy Keller, Haitao Huang, Vladik Kreinovich
;
;  Title of the paper:
;  An IDL/ENVI implementation of the FFT Based Algorithm for Automatic Image Registration
;  
;---------------------------------------------------------------------------------
;
;  INPUT:
;  source (RGB or grayscale): this is a rotated, scaled and shifted version of dest
;  dest (RGB or grayscale): destination image
;  
;  RESTRICTIONS:
;    1. rotation angle must be in ]-90,90[
;    2. scale and shift are limited by the dimensions of the images
;
;  ISSUES: 
;    1. Rotation can't be derived very well from rectangular images
;    
;---------------------------------------------------------------------------------
; 
;  AFFINE TRANSFORMATION IN FOURIER DOMAIN:
;  
;  Function over domain X'
;             g: X' -> Y: Y = g(X')
;  Coordinate transformation (affine):  linear transformation and translation
;             M: X -> X': X' = M(X) = A.X + X0
;                 A = [[s.cos(theta),-s.sin(theta)],[s.sin(angle),s.cos(angle)]]
;                 x' = x.s.cos(theta) - y.s.sin(theta) + x0
;                 y' = x.s.sin(theta) + y.s.cos(theta) + y0
;             M^(-1): X' -> X: X = M^(-1)(X') = A^(-1).(X'-X0)
;             
;             2D transformation and A only scaling/rotating (i.r. rigid):
;                 A = s.Rot(theta)
;                 |A| = s^2
;                 A^T = s.Rot(-theta)
;                 A^(-1) = Rot(-theta)/s
;                 A^(-T) = Rot(theta)/s
;  Same function but now over domain X
;             f: X -> Y: Y = f(X) = g(M(X))
;             g: X' -> Y: Y = g(X') = f(M^(-1)(X'))
;                     
;             X' -----g(X')-----> Y
;             |                 /
;             |               /
;             |             /
;            M(X)         /
;             |        f(X) = g(M(X))
;             |       /
;             |     /
;             |   /
;             | /
;             X
;             
;  Integrate function g over its domain:
;             /int g(X') dX' = /int g(M(X)) |JM(X)| dX        (subsitute X' by X)
;                            = /int g(M(X)) |A| dX            (Jacobian determinant = determinant of A)
;                            = /int f(X) |A| dX
;                            
;  Fourier transformation:
;             G: U' -> G(U'): G(U') = /int g(X') exp[-2.pi.i.<U',X'>] dX'                                     (<.,.> inner product: U'^T.X')
;                                  = /int g(M(X)) exp[-2.pi.i.<U',M(X)>] |JM(X)| dX                       (substitution X'=M(X))
;                                  = |A| . exp[-2.pi.i.U'^T.X0] . /int f(X) exp[-2.pi.i.U'^T.A.X] dX      (matrix representation of M(X) = A.X + X0)
;                                  = |A| . exp[-2.pi.i.U'^T.X0] . /int f(X) exp[-2.pi.i.U^T.X] dX         (U=A^T.U')
;                                  = |A| . exp[-2.pi.i.U'^T.X0] . F(U)                                    (Fourier transform of f)
;                                  = s^2 . exp[-2.pi.i.<U',X0>] . F(s.Rot(-theta).U')                     (A= 2D scaling + rotation)
;                                  
;             F: U -> F(U): F(U) = /int f(X) exp[-2.pi.i.<U,X>] dX                                                      (<.,.> inner product: U'^T.X)
;                               = /int f(M^(-1)(X')) exp[-2.pi.i.<U,M^(-1)(X')>] |JM^(-1)(X')| dX'                      (X = M^(-1)(X'))
;                               = |A^(-1)| . exp[2.pi.i.U^T.A^(-1).X0] . /int g(X') exp[-2.pi.i.U^T.A^(-1).X'] dX      (matrix representation of M^(-1)(X))
;                               = 1/|A| . exp[2.pi.i.U^T.A^(-1).X0] . /int g(X') exp[-2.pi.i.U'^T.X'] dX               (U'=A^(-T).U)
;                               = 1/|A| . exp[2.pi.i.U^T.A^(-1).X0] . G(U')                                            (Fourier transform of g)
;                               = 1/s^2 . exp[2.pi.i.<Rot(theta).U/s,X0>] . G(Rot(theta).U/s)                          (A= 2D scaling + rotation)
;                               
;              As a result we see that:
;                  X' = s   . Rot(theta).X + X0
;                  U' = 1/s . Rot(theta).U
;              This means that the fourier transform is
;                  1. shift invariant
;                  2. preserves rotation
;                  3. inverses scaling
;  
;  Applied on images:
;      Im2(f,source) ---(scale, rotate,shift)---> Im1(g,dest):
;       => g(X') = f(M^(-1)(X'))
;      <=> Im1[X'] = Im2[M^(-1)(X')]
;  
;  Pure translation (HXie_FourierShift):
;                  g(X') = f(M(-1)(X')) = f(X-X0)
;              <=>    G(U) = exp[2.pi.i.<U,-X0>] . F(U)  (U'=U)
;                     This is called the fourier shift property and often written as
;                     FFT(f(X-X0)) = FFT(f(X)).exp[2.pi.i.<U,-X0>]
;              <=> R(U) = G(U).F(U)*/(|G|(U).|F|(U)) = exp[2.pi.i.<U,-X0>] 
;              <=> R(U) = FFT(delta(X-X0))(U)                                    (Fourier transform of the delta function)
;              <=> IFFT(R(U))(X) = IR(X) = delta(X-X0)                            (Inverse fourier transform)
;              <=> maximum of |IR|(X) when X=X0                                (X0 is returned by HXie_FourierShift)
;              
;  Transform images so M becomes a pure translation X0 (HXie_ScaleRot2shift):
;                  g(X') = f(M(-1)(X'))
;              <=> |G|(U') = s^2 |F|(U)
;                      U = s.Rot(-theta).U'
;                       u = u'.cos(theta).s + v'.sin(theta).s
;                       v = -u'.sin(theta).s + v'.cos(theta).s
;                     Log-polar transformation:
;                     <=> R = logb(sqrt(u^2+v^2)) = logb(s) + logb(sqrt(u'^2+v'^2)) = R' + logb(s)
;                     <=> tan(t) = v/u = tan(t'-theta) <=> t = t' - theta 
;             <=>    |G|(R',t') = s^2 |F|(R,t)
;             <=> |G|(R',t') = s^2 |F|(R'+logb(s),t-theta)
;             
;             When considering |G| and |F| in log-polar coordinates as new functions
;             then |F| can be transformed to |G| with a pure translation X0=[-logb(s),theta].
;             Don't forget the minus sign.
;              
;---------------------------------------------------------------------------------

; Initialize source and destination
Im1=HXie_grayscale(dest,error=error)
if error then return,HXIE_inittransfo()
Im2=HXie_grayscale(source,error=error)
error or= n_elements(Im1) eq 1 or n_elements(Im2) eq 1
if error then return,HXIE_inittransfo()

; Make same size
bresize=HXie_equalsize(Im1,Im2,transfo=transfo0)

; Make floating point [0,1]
Im1/=255.
Im2/=255.

; Keep source and estination
Im1_keep=Im1
Im2_keep=Im2

; Convert Im1 and Im2 so they only differ by a shift
HXie_ScaleRot2shift,Im1
HXie_ScaleRot2shift,Im2

; Retrieve the shift difference to transform Im2 -> Im1
; shift = [-logb(scale),theta] in [logb(pixels),pixels]
shift=HXie_FourierShift(Im1,Im2,/logpolar,maxshiftxrange=maxscalerange,maxshiftyrange=maxanglerange,debug=debug)

; Extract scale and rotation angle from shift
s=size(Im2,/dim)
CEN=HXie_zerof(s)
b=HXie_LogPolarBase(s,CEN)
scale = b ^ (-shift[0]) ; R = shift[0] (cfr. HXie_LogPolarTrans)
maxscalerange = b ^ (-maxscalerange)
angle = shift[1]*180./(s[1]-1.) ; t = shift[1]*!pi/(s[1]-1)) + rad -> degrees (cfr. HXie_LogPolarTrans)
maxanglerange *= 180./(s[1]-1.)
transfo=HXIE_inittransfo(angle=angle,scale=scale,X0=CEN[0],Y0=CEN[1])

if keyword_set(debug) then begin
    tmp=indgen(11)-5
    print,'Adjacent angles:',(shift[1]+tmp)*180./(s[1]-1.)
    print,'Adjacent scales:',b ^ (-(shift[0]+tmp))
endif

; Remove rotation and scaling from the source image
Im1=temporary(Im1_keep)
if shift[0] eq 0 and shift[1] eq 0 then Im2=temporary(Im2_keep) $
else Im2=HXie_transformation(temporary(Im2_keep),transfo)

; Retrieve the shift difference to transform Im2 -> Im1
; shift = [shiftx,shifty] in [pixels,pixels]
shift=HXie_FourierShift(Im1,Im2,maxshiftxrange=maxshiftxrange,maxshiftyrange=maxshiftyrange)
transfo.shiftx=shift[0]
transfo.shifty=shift[1]

; Combine with initial transformation
if bresize ne 0 then $
    if bresize lt 0 then transfo=HXIE_combinetransfo(transfo,transfo0) $
    else transfo=HXIE_combinetransfo(transfo0,transfo)

if keyword_set(debug) then begin
    print,'Experimental transformation:'
    help,transfo,/struc
endif

; Return transformation Im2 -> Im1
error=0b
return,transfo
end;function HXie_register
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro HXie_resizebeforeinsert,Im,XL,YL,sdest
; Im: grayscale image
; XL,YL: pixel in destination where [0,0] goes
; sdest: size of destination

if XL lt 0 then begin
    Im=Im[-XL:*,*]
    XL=0
endif
if YL lt 0 then begin
    Im=Im[*,-YL:*]
    YL=0
endif
smod=sdest-round([XL,YL])
sz=size(Im,/dim)
if sz[0] gt smod[0] then Im=Im[0:smod[0]-1,*]
if sz[1] gt smod[1] then Im=Im[*,0:smod[1]-1]
end;pro HXie_resizebeforeinsert
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro HXie_blend_nan,source,dest2,transparancy
ind=where(~finite(source),ct)
source*=(1-transparancy)
source+=dest2*transparancy
if ct ne 0 then source[ind]=dest2[ind]
end;pro HXie_blend_nan
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXie_blend,source,dest,XL,YL,transparancy,error=error

; Make blend

sdest=size(dest,/dim)
dest2=dest
case n_elements(sdest) of
2:    begin
    ; Make transformed source grayscale
    if n_elements(size(source,/dim)) eq 3 then source=HXie_grayscale(source)
    
    ; Make sure a blend can be make
    HXie_resizebeforeinsert,source,XL,YL,sdest
    
    ; Blend
    if n_elements(transparancy) ne 0 then begin
        transparancy=0>transparancy<1
        s=size(source,/dim)
        HXie_blend_nan,source,dest2[XL:XL+s[0]-1,YL:YL+s[1]-1],transparancy
        dest2[XL,YL]=temporary(source)
    endif else begin
        ind=where(~finite(source),ct)
        if ct ne 0 then source[ind]=0
        dest2*=0
        dest2[XL,YL]=temporary(source)
        dest2=[[[dest]],[[dest2]],[[dest2*0]]]
    endelse
    
    endcase
3:    begin
    
    if n_elements(transparancy) eq 0 then transparancy=0.5 $
    else transparancy=0>transparancy<1
    
    indRGB=where(sdest eq 3,comp=tmp)
    sdest=sdest[tmp]
    
    s=size(source,/dim)
    if n_elements(s) eq 3 then begin
        ; source: RGB
        ; Dest: RGB
        
        ; Make sure that channels are at the same dimensions
        indRGB2=where(s eq 3)
        if indRGB2 ne indRGB then begin
            case indRGB+indRGB2 of
            1: tr=[1,0,2]
            2: tr=(indRGB2 eq 2)?[2,0,1]:[1,2,0]
            3: tr=[0,2,1]
            endcase
            source=transpose(source,tr)
        endif
        
        ; Blend channels
        case indRGB of
        0:    begin
            tmp0=reform(source[0,*,*])
            tmp1=reform(source[1,*,*])
            tmp2=reform(source[2,*,*])
            HXie_resizebeforeinsert,tmp0,XL[0],YL[0],sdest
            HXie_resizebeforeinsert,tmp1,XL[0],YL[0],sdest
            HXie_resizebeforeinsert,tmp2,XL,YL,sdest
            source=transpose([[[temporary(tmp0)]],[[temporary(tmp1)]],[[temporary(tmp2)]]],[2,0,1])
            s=size(source,/dim)
            HXie_blend_nan,source,dest2[*,XL:XL+s[1]-1,YL:YL+s[2]-1],transparancy
            dest2[0,XL,YL]=temporary(source)
            endcase
        1:    begin
            tmp0=reform(source[*,0,*])
            tmp1=reform(source[*,1,*])
            tmp2=reform(source[*,2,*])
            HXie_resizebeforeinsert,tmp0,XL[0],YL[0],sdest
            HXie_resizebeforeinsert,tmp1,XL[0],YL[0],sdest
            HXie_resizebeforeinsert,tmp2,XL,YL,sdest
            source=transpose([[[temporary(tmp0)]],[[temporary(tmp1)]],[[temporary(tmp2)]]],[0,2,1])
            s=size(source,/dim)
            HXie_blend_nan,source,dest2[XL:XL+s[0]-1,*,YL:YL+s[2]-1],transparancy
            dest2[XL,0,YL]=temporary(source)
            endcase
        2:    begin
            tmp0=reform(source[*,*,0])
            tmp1=reform(source[*,*,1])
            tmp2=reform(source[*,*,2])
            HXie_resizebeforeinsert,tmp0,XL[0],YL[0],sdest
            HXie_resizebeforeinsert,tmp1,XL[0],YL[0],sdest
            HXie_resizebeforeinsert,tmp2,XL,YL,sdest
            source=[[[temporary(tmp0)]],[[temporary(tmp1)]],[[temporary(tmp2)]]]
            s=size(source,/dim)
            HXie_blend_nan,source,dest2[XL:XL+s[0]-1,YL:YL+s[1]-1,*],transparancy
            dest2[XL,YL,0]=temporary(source)
            endcase
        else: return,0
        endcase
        
    endif else begin
        ; source: grayscale
        ; Dest: RGB
        
        HXie_resizebeforeinsert,source,XL,YL,sdest
        s=size(source,/dim)
        
        ; Blend each channel
        case indRGB of
        0:    begin
            source=reform(source,1,s[0],s[1],/overwrite)
            source3=source
            HXie_blend_nan,source3,dest2[0,XL:XL+s[0]-1,YL:YL+s[1]-1],transparancy
            dest2[0,XL,YL]=temporary(source3)
            source3=source
            HXie_blend_nan,source3,dest2[1,XL:XL+s[0]-1,YL:YL+s[1]-1],transparancy
            dest2[1,XL,YL]=temporary(source3)
            source3=source
            HXie_blend_nan,source3,dest2[2,XL:XL+s[0]-1,YL:YL+s[1]-1],transparancy
            dest2[2,XL,YL]=temporary(source3)
            endcase
        1:    begin
            source=reform(source,s[0],1,s[1],/overwrite)
            source3=source
            HXie_blend_nan,source3,dest2[XL:XL+s[0]-1,0,YL:YL+s[1]-1],transparancy
            dest2[XL,0,YL]=temporary(source3)
            source3=source
            HXie_blend_nan,source3,dest2[XL:XL+s[0]-1,1,YL:YL+s[1]-1],transparancy
            dest2[XL,1,YL]=temporary(source3)
            source3=source
            HXie_blend_nan,source3,dest2[XL:XL+s[0]-1,2,YL:YL+s[1]-1],transparancy
            dest2[XL,2,YL]=temporary(source3)
            endcase
        2:    begin
            source=reform(source,s[0],s[1],1,/overwrite)
            source3=source
            HXie_blend_nan,source3,dest2[XL:XL+s[0]-1,YL:YL+s[1]-1,0],transparancy
            dest2[XL,YL,0]=temporary(source3)
            source3=source
            HXie_blend_nan,source3,dest2[XL:XL+s[0]-1,YL:YL+s[1]-1,1],transparancy
            dest2[XL,YL,1]=temporary(source3)
            source3=source
            HXie_blend_nan,source3,dest2[XL:XL+s[0]-1,YL:YL+s[1]-1,2],transparancy
            dest2[XL,YL,2]=temporary(source3)
            endcase
        endcase
    endelse

    endcase
endcase

error=0b
return,dest2
end;function HXie_blend
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HXie_regandblend, source,dest,transparancy,error=error,debug=debug
error=1b

; Find transformation
transfo=HXie_register(source,dest,error=error,debug=debug)
if error then return,dest

; Transform
source2=float(source) ; because it needs to contain NaNs
HXie_transform,source2,transfo,error=error,XL=XL,YL=YL,interp=0,missing=!values.F_NAN,/resize
if error then return,dest

; Blend
return,HXie_blend(source2,dest,XL,YL,transparancy,error=error)
end;function HXie_regandblend
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro HXie_tvsclblend,i,img

screen=GET_SCREEN_SIZE()/2
s=size(img,/dim)

if n_elements(s) eq 3 then begin
    ind=where(s eq 3,comp=ind2)
    
    mag=screen/float(s[ind2])
    mag=mag[0]<mag[1]<1

    window,i,xsize=s[ind2[0]]*mag,ysize=s[ind2[1]]*mag
    case ind[0] of
    0: tvscl,congrid(img,s[0],mag*s[1],mag*s[2]),true=1
    1: tvscl,congrid(img,mag*s[0],s[1],mag*s[2]),true=2
    2: tvscl,congrid(img,mag*s[0],mag*s[1],s[2]),true=3
    endcase
endif else begin
    mag=screen/float(s)
    mag=mag[0]<mag[1]<1
    window,i,xsize=s[0]*mag,ysize=s[1]*mag
    tvscl,congrid(img,mag*s[0],mag*s[1])
endelse
end;pro HXie_tvsclblend
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
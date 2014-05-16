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

function gfunctest_wnoise,x,function_name,params,gaussian=gaussian
call_procedure, function_name, x, params, y

if keyword_set(gaussian) then y += randomn(seed, n_elements(x)) * sqrt(y) $
else for i=0l,n_elements(x)-1 do y[i] += (y[i] eq 0)?0:randomn(seed,poisson=y[i])
return,y
end;function gfunctest_wnoise
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gausstest_wnoise,b,e,n,area,position,sigma,x=x,gaussian=gaussian
inc=(e-b)/(n-1.)
x=b+inc*lindgen(n)

y=x*0
for i=0l,n_elements(area)-1 do $
    y+=gfunctest_wnoise(x,'gfunc_gauss',[area[i],position[i],sigma[i]],gaussian=gaussian)

return,y
end;function gausstest_wnoise
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc_stepfunc,X,A,F,PDER, struc=struc
; F = step1(X,A[0:2]) (+ A[3] + A[4]*X + A[5]*X^2)
; A[0]: area of first stepfunction
; A[1]: position of first stepfunction
; A[2]: sigma of first stepfunction
; A[3:5]: optional: constant, linear or parabolic

n = n_elements(a)
nx = N_ELEMENTS(x)

; First edge
if A[2] ne 0.0 then begin
    Z1 = (X-A[1])/A[2]
    EZ1 = ERF(Z1/sqrt(2))
endif else begin
    ; Z=-Inf or Inf => EZ=0 or 1
    Z1 = REPLICATE(FIX(100+A[1], TYPE=SIZE(x,/TYPE)), nx)
    EZ1 = Z1*0
    ind=where(X gt A[1],ct)
    if ct ne 0 then EZ1[ind]++
endelse
EZN1=(1+EZ1)/2.

; Step1+a+bx+cx^2
case n of
    3:  F = A[0]*EZN1
    4:  F = A[0]*EZN1 + A[3]
    5:  F = A[0]*EZN1 + A[3] + A[4]*X
    6:  F = A[0]*EZN1 + A[3] + A[4]*X + A[5]*X^2
ENDCASE

IF N_PARAMS(0) LE 3 THEN RETURN ;NEED PARTIAL?

PDER = FLTARR(nx, n)
PDER[*,0] = EZN1
if A[2] ne 0. then begin
    PDER[*,1] = -A[0]/(sqrt(2*!pi)*A[2])*exp(-Z1^2./2)
    PDER[*,2] = PDER[*,1] * Z1
endif
if n gt 3 then PDER[*,3] = 1.
if n gt 4 then PDER[*,4] = X
if n gt 5 then PDER[*,5] = X^2
end;pro gfunc_stepfunc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc_bellfunc,X,A,F,PDER, struc=struc
; F = step1(X,A[0:2]) - step2(X,A[3:5]) (+ A[6] + A[7]*X + A[8]*X^2)
; A[0]: area of first stepfunction
; A[1]: position of first stepfunction
; A[2]: sigma of first stepfunction
; A[3]: area of second stepfunction
; A[4]: position of second stepfunction
; A[5]: sigma of second stepfunction
; A[6:8]: optional: constant, linear or parabolic

n = n_elements(a)
nx = N_ELEMENTS(x)

; First edge
if A[2] ne 0.0 then begin
    Z1 = (X-A[1])/A[2]
    EZ1 = ERF(Z1/sqrt(2))
endif else begin
    ; Z=-Inf or Inf => EZ=0 or 1
    Z1 = REPLICATE(FIX(100+A[1], TYPE=SIZE(x,/TYPE)), nx)
    EZ1 = Z1*0
    ind=where(X gt A[1],ct)
    if ct ne 0 then EZ1[ind]++
endelse
EZN1=(1+EZ1)/2.

; Second edge
if A[5] ne 0.0 then begin
    Z2 = (X-A[4])/A[5]
    EZ2 = ERF(Z2/sqrt(2))
endif else begin
    ; Z=-Inf or Inf => EZ=0 or 1
    Z2 = REPLICATE(FIX(100+A[4], TYPE=SIZE(x,/TYPE)), nx)
    EZ2 = Z2*0
    ind=where(X gt A[4],ct)
    if ct ne 0 then EZ2[ind]++
endelse
EZN2=(1+EZ2)/2.

; Step1-Step2+ax+b
case n of
    6:  F = A[0]*EZN1 - A[3]*EZN2
    7:  F = A[0]*EZN1 - A[3]*EZN2 + A[6]
    8:  F = A[0]*EZN1 - A[3]*EZN2 + A[6] + A[7]*X
    9:  F = A[0]*EZN1 - A[3]*EZN2 + A[6] + A[7]*X + A[8]*X^2
ENDCASE

IF N_PARAMS(0) LE 3 THEN RETURN ;NEED PARTIAL?

PDER = FLTARR(nx, n)
PDER[*,0] = EZN1
if A[2] ne 0. then begin
    PDER[*,1] = -A[0]/(sqrt(2*!pi)*A[2])*exp(-Z1^2./2)
    PDER[*,2] = PDER[*,1] * Z1
endif
PDER[*,3] = -EZN2
if A[5] ne 0. then begin
    PDER[*,4] = A[3]/(sqrt(2*!pi)*A[5])*exp(-Z2^2./2)
    PDER[*,5] = PDER[*,4] * Z2
endif
if n gt 6 then PDER[*,6] = 1.
if n gt 7 then PDER[*,7] = X
if n gt 8 then PDER[*,8] = X^2
end;pro gfunc_bellfunc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PRO gfunc_PearsonVII,X,A,F,PDER, struc=struc
; F = PearsonVII + (+ A[4] + A[5]*X + A[6]*X^2)
n = n_elements(a)
nx = N_ELEMENTS(x)

xu=A[1]
W=A[2]
R=A[3]

c=2.^(1./R)-1
gR=exp(lngamma(R)-lngamma(R-0.5))
dgR=digamma(R)-digamma(R-0.5)
AA=2*gR*sqrt(c)/sqrt(!pi)

Barg=1+4*c*((X-xu)/W)^2
B=Barg^(-R)
multpderI=AA/W*B

case n of
    4:  F = A[0]*multpderI
    5:  F = A[0]*multpderI + A[4]
    6:  F = A[0]*multpderI + A[4] + A[5]*X
    7:  F = A[0]*multpderI + A[4] + A[5]*X + A[6]*X^2
ENDCASE

IF N_PARAMS(0) LE 3 THEN RETURN ;NEED PARTIAL?

PDER = FLTARR(nx, n)
PDER[*,0] = multpderI
PDER[*,1] = (AA*A[0]/W^3*8*R*c)*Barg^(-R-1)*(X-xu)
PDER[*,2] = (AA*A[0]/W^4*8*R*c)*Barg^(-R-1)*(X-xu)^2. - (AA*A[0]/W^2)*B
PDER[*,3] = (A[0]*AA/W)*B*(dgR-alog(Barg)+(4*2^(1./R)*alog(2)/R/W^2)*(X-xu)^2./Barg)- $
            (gR*A[0]*2^(1./R)*alog(2)/sqrt(!pi*c)/W/R^2)*B

if n gt 4 then PDER[*,4] = 1.
if n gt 5 then PDER[*,5] = X
if n gt 6 then PDER[*,6] = X^2
END;PRO gfunc_PearsonVII
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PRO gfunc_gauss,X,A,F,PDER, struc=struc
; F = gaussian + (+ A[3] + A[4]*X + A[5]*X^2)
n = n_elements(a)
nx = N_ELEMENTS(x)

if a[2] ne 0.0 then begin
    Z = (X-A[1])/A[2]   ;GET Z
    EZ = EXP(-Z^2/2.)   ;GAUSSIAN PART
    if keyword_set(struc) then EZ/=sqrt(2*!pi)*A[2] ; NORMALIZE
endif else begin
    z = REPLICATE(FIX(100, TYPE=SIZE(x,/TYPE)), nx)
    ez = z*0
endelse

case n of
    3:  F = A[0]*EZ
    4:  F = A[0]*EZ + A[3]
    5:  F = A[0]*EZ + A[3] + A[4]*X
    6:  F = A[0]*EZ + A[3] + A[4]*X + A[5]*X^2 ;FUNCTIONS.
ENDCASE

IF N_PARAMS(0) LE 3 THEN RETURN ;NEED PARTIAL?

PDER = FLTARR(nx, n) ;YES, MAKE ARRAY.
PDER[*,0] = EZ      ;COMPUTE PARTIALS
if a[2] ne 0. then begin
    PDER[*,1] = A[0] * EZ * Z/A[2]
    if keyword_set(struc) then PDER[*,2] = A[0] * EZ * (Z^2-1)/A[2] $
    else PDER[*,2] = PDER[*,1] * Z
endif
if n gt 3 then PDER[*,3] = 1.
if n gt 4 then PDER[*,4] = X
if n gt 5 then PDER[*,5] = X^2
END;PRO gfunc_gauss
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PVconstr,x,a,m,s,n
return,a*((1-n)/sqrt(2*!dpi)/s*exp(-(x-m)^2./2/s^2.)+$
         n*s/(2*!dpi)/((x-m)^2.+(s/2)^2.))
end;function PVconstr
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncBackVar, X, A, F,pder
area=A[0]
position=A[1]
sigma=A[2]
F=area/(sqrt(2*!dpi)*sigma)*exp(-(X-position)^2/(2*sigma^2))
if N_PARAMS() GE 4 then begin
pder=dblarr(n_elements(x),n_elements(A))
;each row is a partial derivative with respect to one of the parameter
pder[*,0]=1/(sqrt(2*!dpi)*sigma)*exp(-(X-position)^2/(2*sigma^2))
pder[*,1]=area*pder[*,0]*(X-position)/sigma^2
pder[*,2]=area*pder[*,0]*(-1/sigma+(X-position)^2/sigma^3)
endif
end;pro gfuncBackVar
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncSin, X, A, F,pder,struc=struc

iA=0
iAfix=0
i=0
AA=(struc.bfix[i++]?struc.Afix[iAfix++]:A[iA++])
BB=(struc.bfix[i++]?struc.Afix[iAfix++]:A[iA++])
CC=(struc.bfix[i++]?struc.Afix[iAfix++]:A[iA++])
DD=(struc.bfix[i++]?struc.Afix[iAfix++]:A[iA++])
EE=(struc.bfix[i++]?struc.Afix[iAfix++]:A[iA++])
FF=(struc.bfix[i++]?struc.Afix[iAfix++]:A[iA++])

N=DD*DD*sin(X-FF)^2.+EE*EE*cos(X-FF)^2.
if (DD eq 0) and (EE eq 0) then N[*]=1
R=DD*EE/sqrt(N)
F = AA + (BB+R)*cos(X + CC)

if N_PARAMS() GE 4 then begin
pder=fltarr(n_elements(X),n_elements(A))

i=0
j=0
if ~struc.bfix[i++] then pder[*,j++]=1
if ~struc.bfix[i++] then pder[*,j++]=cos(X + CC)
if ~struc.bfix[i++] then pder[*,j++]=-(BB+R)*sin(X + CC)
if ~struc.bfix[i++] then pder[*,j++]=R*(1./DD-DD*sin(X-FF)^2./N)*cos(X + CC)
if ~struc.bfix[i++] then pder[*,j++]=R*(1./EE-EE*cos(X-FF)^2./N)*cos(X + CC)
if ~struc.bfix[i++] then pder[*,j++]=R/N*sin(X-FF)*cos(X-FF)*(DD*DD-EE*EE)*cos(X + CC)

endif

end;pro gfuncSin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc1, X, A, F,struc,pder
pol=struc.(0)
GR=struc.(1)
apeaks=(n_elements(A)-GR-1)/3
indar=indgen(apeaks)*3
indpos=indar+1
indsigma=indar+2
indcoeff=indgen(GR+1)+3*apeaks
area=A[indar]
position=A[indpos]
sigma=A[indsigma]
coeff=A[indcoeff]
F=X*0.
for k=0l,apeaks-1 do F=F+area[k]/(sqrt(2*!dpi)*sigma[k])*exp(-(X-position[k])^2/(2*sigma[k]^2))
for k=0l,GR do F=F+coeff[k]*pol[*,k]
if N_PARAMS() GE 5 then begin
pder=dblarr(n_elements(x),n_elements(A))
;each row is a partial derivative with respect to one of the parameter
for k=0l,apeaks-1 do begin
    j=3*k
    pder[*,j]=1/(sqrt(2*!dpi)*sigma[k])*exp(-(X-position[k])^2/(2*sigma[k]^2))
    pder[*,j+1]=area[k]*pder[*,j]*(X-position[k])/(sigma[k]^2)
    pder[*,j+2]=area[k]*pder[*,j]*(-1/sigma[k]+(X-position[k])^2/sigma[k]^3)
endfor
for k=0l,GR do begin
    pder[*,k+apeaks*3]=pol[*,k]
endfor
endif
end;pro gfunc1
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc2, X, A, F,pder,dummy=dummy
area=A[0]
position=A[1]
sigma=A[2]
F=area/(sqrt(2*!pi)*sigma)*exp(-(X-position)^2/(2*sigma^2))+A[3]+A[4]*X
if N_PARAMS() GE 4 then begin
pder=fltarr(n_elements(x),n_elements(A))
;each row is a partial derivative with respect to one of the parameter
pder[*,0]=1/(sqrt(2*!dpi)*sigma)*exp(-(X-position)^2/(2*sigma^2))
pder[*,1]=area*pder[*,0]*(X-position)/sigma^2
pder[*,2]=area*pder[*,0]*(-1/sigma+(X-position)^2/sigma^3)
pder[*,3]=0.*X+1.
pder[*,4]=X
endif
end;pro gfunc2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc3, X, A, F,pder
apeaks=(n_elements(A))/3
indar=indgen(apeaks)*3
indpos=indar+1
indsigma=indar+2
area=A[indar]
position=A[indpos]
sigma=A[indsigma]
F=X*0.
for k=0l,apeaks-1 do F=F+area[k]/(sqrt(2*!dpi)*sigma[k])*exp(-(X-position[k])^2/(2*sigma[k]^2))
if N_PARAMS() GE 4 then begin
pder=dblarr(n_elements(x),n_elements(A))
;each row is a partial derivative with respect to one of the parameter
for k=0l,apeaks-1 do begin
    j=3*k
    pder[*,j]=1/(sqrt(2*!dpi)*sigma[k])*exp(-(X-position[k])^2/(2*sigma[k]^2))
    pder[*,j+1]=area[k]*pder[*,j]*(X-position[k])/(sigma[k]^2)
    pder[*,j+2]=area[k]*pder[*,j]*(-1/sigma[k]+(X-position[k])^2/sigma[k]^3)
endfor
endif
end;pro gfunc3
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc4, X, A, F,pder
; Pseudo Voigt: PV=(1-n)G+nL
apeaks=(n_elements(A)-1)/3
indar=indgen(apeaks)*3
indpos=indar+1
indsigma=indar+2
area=A[indar]
position=A[indpos]
sigma=A[indsigma]
n=A[n_elements(A)-1]
F=X*0.
for k=0l,apeaks-1 do F=F+area[k]*((1-n)/(sqrt(2*!dpi)*sigma[k])*exp(-(X-position[k])^2/(2*sigma[k]^2))+$
                     n*sigma[k]/(2*!dpi)/((X-position[k])^2+(sigma[k]/2)^2) )
if N_PARAMS() GE 4 then begin
pder=dblarr(n_elements(x),n_elements(A))
;each row is a partial derivative with respect to one of the parameter
for k=0l,apeaks-1 do begin
    j=3*k
    G=1/(sqrt(2*!dpi)*sigma[k])*exp(-(X-position[k])^2/(2*sigma[k]^2))
    L=sigma[k]/(2*!dpi)/((X-position[k])^2+(sigma[k]/2)^2)
    pder[*,j]=(1-n)*G+n*L
    pder[*,j+1]=area[k]*((1-n)*G*(X-position[k])/(sigma[k]^2)+$
                  n*L*2*(X-position[k])/((X-position[k])^2+(sigma[k]/2)^2) )
    pder[*,j+2]=area[k]*((1-n)*G*(-1/sigma[k]+(X-position[k])^2/sigma[k]^3)+$
                  n*L*(1/sigma[k]+!dpi*L)  )
    pder[*,n_elements(A)-1]=pder[*,n_elements(A)-1]+area[k]*(L-G)
endfor
endif
end;pro gfunc4
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc6, X, A, F,pder
; Pseudo Voigt: PV=(1-n)G+nL
apeaks=(n_elements(A))/4
indar=indgen(apeaks)*3
indpos=indar+1
indsigma=indar+2
area=A[indar]
position=A[indpos]
sigma=A[indsigma]
n=A[3*apeaks:*]
F=X*0.
for k=0l,apeaks-1 do F=F+area[k]*((1-n[k])/(sqrt(2*!dpi)*sigma[k])*exp(-(X-position[k])^2/(2*sigma[k]^2))+$
                     n[k]*sigma[k]/(2*!dpi)/((X-position[k])^2+(sigma[k]/2)^2) )
if N_PARAMS() GE 4 then begin
pder=dblarr(n_elements(x),n_elements(A))
;each row is a partial derivative with respect to one of the parameter
for k=0l,apeaks-1 do begin
    j=3*k
    G=1/(sqrt(2*!dpi)*sigma[k])*exp(-(X-position[k])^2/(2*sigma[k]^2))
    L=sigma[k]/(2*!dpi)/((X-position[k])^2+(sigma[k]/2)^2)
    pder[*,j]=(1-n[k])*G+n[k]*L
    pder[*,j+1]=area[k]*((1-n[k])*G*(X-position[k])/(sigma[k]^2)+$
                  n[k]*L*2*(X-position[k])/((X-position[k])^2+(sigma[k]/2)^2) )
    pder[*,j+2]=area[k]*((1-n[k])*G*(-1/sigma[k]+(X-position[k])^2/sigma[k]^3)+$
                  n[k]*L*(1/sigma[k]+!dpi*L)  )
    pder[*,k+3*apeaks]=area[k]*(L-G)
endfor
endif
end;pro gfunc6
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc6c, X, A, F,struc,pder
; Pseudo Voigt: PV=(1-n)G+nL
; with constraints on all param

apeaks=(n_elements(A))/4
cj=struc.j
constrar=struc.constrar

cte=2*sqrt(2*alog(2))
indar=indgen(apeaks)*4
indpos=indar+1
indsigma=indar+2
indn=indar+3

; calculate physical meaningful parameters AA
; A is fit parameter
AA=reform(A,4,apeaks)
AAA=AA;=> fit parameters
constrar=struc.constrar
for i=0l,n_elements(cj)-1 do $
AA[cj[i],*]=constrar[2*cj[i],*]+2.*constrar[2*cj[i]+1,*]/!dpi*atan(AA[cj[i],*])
AA=reform(AA,4*apeaks);=> physical meaningful parameters
AA[indsigma]=AA[indsigma]

area=AA[indar]
position=AA[indpos]
sigma=AA[indsigma]
n=AA[indn]

F=X*0.
for k=0l,apeaks-1 do F=F+area[k]*((1-n[k])/(sqrt(2*!dpi)*sigma[k])*exp(-(X-position[k])^2/(2*sigma[k]^2))+$
                     n[k]*sigma[k]/(2*!dpi)/((X-position[k])^2+(sigma[k]/2)^2) )
if N_PARAMS() GE 5 then begin
pder=dblarr(n_elements(x),n_elements(A))
;each row is a partial derivative with respect to one of the parameter
for k=0l,apeaks-1 do begin
    j=4*k
    dAAdA=fltarr(4)+1
    dAAdA[cj]=2/!dpi*constrar[cj*2+1,k]/(1+AAA[cj,k]^2.)

    G=1/(sqrt(2*!dpi)*sigma[k])*exp(-(X-position[k])^2/(2*sigma[k]^2))
    L=sigma[k]/(2*!dpi)/((X-position[k])^2+(sigma[k]/2)^2)
    pder[*,j]=dAAdA[0]*((1-n[k])*G+n[k]*L)
    pder[*,j+1]=dAAdA[1]*area[k]*((1-n[k])*G*(X-position[k])/(sigma[k]^2)+$
                  n[k]*L*2*(X-position[k])/((X-position[k])^2+(sigma[k]/2)^2) )
    pder[*,j+2]=dAAdA[2]*area[k]*((1-n[k])*G*(-1/sigma[k]+(X-position[k])^2/sigma[k]^3)+$
                  n[k]*L*(1/sigma[k]+!dpi*L)  )
    pder[*,j+3]=dAAdA[3]*area[k]*(L-G)
endfor
endif
end;pro gfunc6c
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncTilt0, Xin, AA, F, pder, struc=struc
; optional fit parameters: a, b, center, dist

; ----Parameters (fixed and variable)----
nAA=n_elements(AA)
offset1=0
offset2=0
if struc.refbool[0] then begin
    a=AA[offset1]
    offset1++
endif else begin
    a=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[1] then begin
    b=AA[offset1]
    offset1++
endif else begin
    b=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[2] then begin
    xc=AA[offset1]
    offset1++
    yc=AA[offset1]
    offset1++
endif else begin
    xc=struc.Afix[offset2]
    offset2++
    yc=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[3] then begin
    di=AA[offset1]
    offset1++
endif else begin
    di=struc.Afix[offset2]
    offset2++
endelse

xx=Xin[*,0]
yy=Xin[*,1]
nx=n_elements(xx)

xd=xx-xc
yd=yy-yc
xt=xd*cos(b)-yd*sin(b)
yt=xd*sin(b)+yd*cos(b)

T=xt^2.+yt^2.*cos(a)^2.
N=di-sin(a)*yt
F=atan(sqrt(T),N)

if N_PARAMS() GE 4 then begin
pder=dblarr(nx,nAA)
;each row is a partial derivative with respect to one of the parameter

P=1./(1+T/N^2.)
c1=-sqrt(T)/N^2.*P
c2=P/(2*N*sqrt(T))
offset=0

if struc.refbool[0] then begin
    Ta=-2.*yt^2.*cos(a)*sin(a)
    Na=-yt*cos(a)
    pder[*,offset]=(c1*Na + c2*Ta)
    offset++
endif

if struc.refbool[1] then begin
    Tb=-2.*xt*yt*sin(a)^2.
    Nb=-xt*sin(a)
    pder[*,offset]=(c1*Nb + c2*Tb)
    offset++
endif

if struc.refbool[2] then begin
    Txc=-2*xt*cos(b)-2*yt*sin(b)*cos(a)^2.
    Nxc=sin(a)*sin(b)
    pder[*,offset]=(c1*Nxc + c2*Txc)
    offset++
    Tyc=2*xt*sin(b)-2*yt*cos(b)*cos(a)^2.
    Nyc=sin(a)*cos(b)
    pder[*,offset]=(c1*Nyc + c2*Tyc)
    offset++
endif

if struc.refbool[3] then begin
    pder[*,offset]=c1
    offset++
endif


endif

end;pro gfuncTilt0
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncTilt1, Xin, AA, F, pder, struc=struc
; fit parameters: d-spacings
; optional fit parameters: a, b, center, dist

; ----Parameters (fixed and variable)----
nAA=n_elements(AA)
offset1=0
offset2=0
if struc.refbool[0] then begin
    a=AA[offset1]
    offset1++
endif else begin
    a=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[1] then begin
    b=AA[offset1]
    offset1++
endif else begin
    b=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[2] then begin
    xc=AA[offset1]
    offset1++
    yc=AA[offset1]
    offset1++
endif else begin
    xc=struc.Afix[offset2]
    offset2++
    yc=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[3] then begin
    di=AA[offset1]
    offset1++
endif else begin
    di=struc.Afix[offset2]
    offset2++
endelse

d=AA[offset1:nAA-1]

xx=Xin[*,0]
yy=Xin[*,1]
nx=n_elements(xx)
nj=struc.nj
M=struc.M

xd=xx-xc
yd=yy-yc
xt=xd*cos(b)-yd*sin(b)
yt=xd*sin(b)+yd*cos(b)

T=xt^2.+yt^2.*cos(a)^2.
N=di-sin(a)*yt
X=0.5*atan(sqrt(T),N)
pderd=sin(X)

F=pderd
ind=0
d2=fltarr(nx)
for j=0l,M-1 do begin
    d2[ind:ind+nj[j]-1]=d[j]
    ind+=nj[j]
endfor
F*=d2

if N_PARAMS() GE 4 then begin
pder=dblarr(nx,nAA)
;each row is a partial derivative with respect to one of the parameter

P=0.5*d2*cos(X)/(1+T/N^2.)
c1=-sqrt(T)/N^2.*P
c2=P/(2*N*sqrt(T))
offset=0

if struc.refbool[0] then begin
    Ta=-2.*yt^2.*cos(a)*sin(a)
    Na=-yt*cos(a)
    pder[*,offset]=(c1*Na + c2*Ta)
    offset++
endif

if struc.refbool[1] then begin
    Tb=-2.*xt*yt*sin(a)^2.
    Nb=-xt*sin(a)
    pder[*,offset]=(c1*Nb + c2*Tb)
    offset++
endif

if struc.refbool[2] then begin
    Txc=-2*xt*cos(b)-2*yt*sin(b)*cos(a)^2.
    Nxc=sin(a)*sin(b)
    pder[*,offset]=(c1*Nxc + c2*Txc)
    offset++
    Tyc=2*xt*sin(b)-2*yt*cos(b)*cos(a)^2.
    Nyc=sin(a)*cos(b)
    pder[*,offset]=(c1*Nyc + c2*Tyc)
    offset++
endif

if struc.refbool[3] then begin
    pder[*,offset]=c1
    offset++
endif

ind=0
for j=0l,M-1 do begin
    pder[ind:ind+nj[j]-1,j+offset]=pderd[ind:ind+nj[j]-1]
    ind+=nj[j]
endfor

endif

end;pro gfuncTilt1
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncTilt2, Xin, AA, F, pder, struc=struc
; fit parameters: wavelength
; optional fit parameters: a, b, center, dist

; ----Parameters (fixed and variable)----
nAA=n_elements(AA)
offset1=0
offset2=0
if struc.refbool[0] then begin
    a=AA[offset1]
    offset1++
endif else begin
    a=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[1] then begin
    b=AA[offset1]
    offset1++
endif else begin
    b=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[2] then begin
    xc=AA[offset1]
    offset1++
    yc=AA[offset1]
    offset1++
endif else begin
    xc=struc.Afix[offset2]
    offset2++
    yc=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[3] then begin
    di=AA[offset1]
    offset1++
endif else begin
    di=struc.Afix[offset2]
    offset2++
endelse

lambda=AA[offset1]

xx=Xin[*,0]
yy=Xin[*,1]
nx=n_elements(xx)

xd=xx-xc
yd=yy-yc
xt=xd*cos(b)-yd*sin(b)
yt=xd*sin(b)+yd*cos(b)

T=xt^2.+yt^2.*cos(a)^2.
N=di-sin(a)*yt
X=0.5*atan(sqrt(T),N)
pderl=1./sin(X)

F=lambda*pderl

if N_PARAMS() GE 4 then begin
pder=dblarr(nx,nAA)
;each row is a partial derivative with respect to one of the parameter

P=-0.5*lambda*cos(X)/sin(X)^2./(1+T/N^2.)
c1=-sqrt(T)/N^2.*P
c2=P/(2*N*sqrt(T))
offset=0

if struc.refbool[0] then begin
    Ta=-2.*yt^2.*cos(a)*sin(a)
    Na=-yt*cos(a)
    pder[*,offset]=(c1*Na + c2*Ta)
    offset++
endif

if struc.refbool[1] then begin
    Tb=-2.*xt*yt*sin(a)^2.
    Nb=-xt*sin(a)
    pder[*,offset]=(c1*Nb + c2*Tb)
    offset++
endif

if struc.refbool[2] then begin
    Txc=-2*xt*cos(b)-2*yt*sin(b)*cos(a)^2.
    Nxc=sin(a)*sin(b)
    pder[*,offset]=(c1*Nxc + c2*Txc)
    offset++
    Tyc=2*xt*sin(b)-2*yt*cos(b)*cos(a)^2.
    Nyc=sin(a)*cos(b)
    pder[*,offset]=(c1*Nyc + c2*Tyc)
    offset++
endif

if struc.refbool[3] then begin
    pder[*,offset]=c1
    offset++
endif

pder[*,offset]=pderl

endif

end;pro gfuncTilt2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncTilt3, Xin, AA, F, pder, struc=struc
; fit parameters: wavelength and d-spacing
; optional fit parameters: a, b, center, dist

; ----Parameters (fixed and variable)----
nAA=n_elements(AA)
offset1=0
offset2=0
if struc.refbool[0] then begin
    a=AA[offset1]
    offset1++
endif else begin
    a=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[1] then begin
    b=AA[offset1]
    offset1++
endif else begin
    b=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[2] then begin
    xc=AA[offset1]
    offset1++
    yc=AA[offset1]
    offset1++
endif else begin
    xc=struc.Afix[offset2]
    offset2++
    yc=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[3] then begin
    di=AA[offset1]
    offset1++
endif else begin
    di=struc.Afix[offset2]
    offset2++
endelse

lambda=AA[offset1]
d=AA[offset1+1:nAA-1]

xx=Xin[*,0]
yy=Xin[*,1]
nx=n_elements(xx)
nj=struc.nj
M=struc.M

xd=xx-xc
yd=yy-yc
xt=xd*cos(b)-yd*sin(b)
yt=xd*sin(b)+yd*cos(b)

T=xt^2.+yt^2.*cos(a)^2.
N=di-sin(a)*yt
X=0.5*atan(sqrt(T),N)

pderd=sin(X)/lambda

F=pderd
ind=0
d2=fltarr(nx)
for j=0l,M-1 do begin
    d2[ind:ind+nj[j]-1]=d[j]
    ind=ind+nj[j]
endfor
F=F*d2

if N_PARAMS() GE 4 then begin
pder=dblarr(nx,nAA)
;each row is a partial derivative with respect to one of the parameter

P=0.5*d2*cos(X)/lambda/(1+T/N^2.)
c1=-sqrt(T)/N^2.*P
c2=P/(2*N*sqrt(T))
offset=0

if struc.refbool[0] then begin
    Ta=-2.*yt^2.*cos(a)*sin(a)
    Na=-yt*cos(a)
    pder[*,offset]=(c1*Na + c2*Ta)
    offset++
endif

if struc.refbool[1] then begin
    Tb=-2.*xt*yt*sin(a)^2.
    Nb=-xt*sin(a)
    pder[*,offset]=(c1*Nb + c2*Tb)
    offset++
endif

if struc.refbool[2] then begin
    Txc=-2*xt*cos(b)-2*yt*sin(b)*cos(a)^2.
    Nxc=sin(a)*sin(b)
    pder[*,offset]=(c1*Nxc + c2*Txc)
    offset++
    Tyc=2*xt*sin(b)-2*yt*cos(b)*cos(a)^2.
    Nyc=sin(a)*cos(b)
    pder[*,offset]=(c1*Nyc + c2*Tyc)
    offset++
endif

if struc.refbool[3] then begin
    pder[*,offset]=c1
    offset++
endif

pder[*,offset]=-d2*pderd/lambda
offset++

ind=0
for j=0l,M-1 do begin
    pder[ind:ind+nj[j]-1,j+offset]=pderd[ind:ind+nj[j]-1]
    ind+=nj[j]
endfor

endif

end;pro gfuncTilt3
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncTilt4, Xin, AA, F, pder, struc=struc
; optional fit parameters: a, b, center, dist

; ----Parameters (fixed and variable)----
nAA=n_elements(AA)
offset1=0
offset2=0
if struc.refbool[0] then begin
    a=AA[offset1]
    offset1++
endif else begin
    a=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[1] then begin
    b=AA[offset1]
    offset1++
endif else begin
    b=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[2] then begin
    xc=AA[offset1]
    offset1++
    yc=AA[offset1]
    offset1++
endif else begin
    xc=struc.Afix[offset2]
    offset2++
    yc=struc.Afix[offset2]
    offset2++
endelse

if struc.refbool[3] then begin
    di=AA[offset1]
    offset1++
endif else begin
    di=struc.Afix[offset2]
    offset2++
endelse

tt=AA[offset1:nAA-1]

xx=Xin[*,0]
yy=Xin[*,1]
nx=n_elements(xx)
nj=struc.nj
M=struc.M

xd=xx-xc
yd=yy-yc
xt=xd*cos(b)-yd*sin(b)
yt=xd*sin(b)+yd*cos(b)

T=xt^2.+yt^2.*cos(a)^2.
N=di-sin(a)*yt
F=atan(sqrt(T),N)

ind=0
for j=0l,M-1 do begin
    F[ind:ind+nj[j]-1]=F[ind:ind+nj[j]-1]-tt[j]
    ind=ind+nj[j]
endfor

if N_PARAMS() GE 4 then begin
pder=dblarr(nx,nAA)
;each row is a partial derivative with respect to one of the parameter

P=1./(1+T/N^2.)
c1=-sqrt(T)/N^2.*P
c2=P/(2*N*sqrt(T))
offset=0

if struc.refbool[0] then begin
    Ta=-2.*yt^2.*cos(a)*sin(a)
    Na=-yt*cos(a)
    pder[*,offset]=(c1*Na + c2*Ta)
    offset++
endif

if struc.refbool[1] then begin
    Tb=-2.*xt*yt*sin(a)^2.
    Nb=-xt*sin(a)
    pder[*,offset]=(c1*Nb + c2*Tb)
    offset++
endif

if struc.refbool[2] then begin
    Txc=-2*xt*cos(b)-2*yt*sin(b)*cos(a)^2.
    Nxc=sin(a)*sin(b)
    pder[*,offset]=(c1*Nxc + c2*Txc)
    offset++
    Tyc=2*xt*sin(b)-2*yt*cos(b)*cos(a)^2.
    Nyc=sin(a)*cos(b)
    pder[*,offset]=(c1*Nyc + c2*Tyc)
    offset++
endif

if struc.refbool[3] then begin
    pder[*,offset]=c1
    offset++
endif

ind=0
for j=0l,M-1 do begin
    pder[ind:ind+nj[j]-1,j+offset]=-1
    ind+=nj[j]
endfor

endif

end;pro gfuncTilt4
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc2dgbis, X, A, F,pder,struc=struc
ar=A[0]
mx=A[1]
my=A[2]
sx=A[3]
sy=A[4]
sxy=A[5]
bkg0=A[6]
corr=sxy/(sx*sy)
correx=1.-corr^2.
xmap=X[*,0]
ymap=X[*,1]
z=(xmap-mx)^2./sx^2.-2.*corr*(xmap-mx)*(ymap-my)/sx/sy+(ymap-my)^2./sy^2.
E=exp(-z/(2.*correx))
F=ar/(2*!pi*sx*sy*sqrt(correx))*E+bkg0


if N_PARAMS() GE 4 then begin
pder=dblarr(n_elements(F),n_elements(A))
;each row is a partial derivative with respect to one of the parameter

pder[*,0]= E/(2*!pi*sx*sy*sqrt(correx))

pder[*,1]= ar*E/(2*!pi*sy*sx^2.*correx^1.5)* $
( (xmap-mx)/sx - corr*(ymap-my)/sy )

pder[*,2]= ar*E/(2*!pi*sx*sy^2.*correx^1.5)* $
( (ymap-my)/sy - corr*(xmap-mx)/sx )

pder[*,3]= ar*E/(2*!pi*sy*sx^2*correx^1.5)* $
(-1. + (xmap-mx)^2./sx^2. - (xmap-mx)*(ymap-my)*corr/(sy*sx) + z*corr^2./correx )

pder[*,4]= ar*E/(2*!pi*sx*sy^2*correx^1.5)* $
(-1. + (ymap-my)^2./sy^2. - (ymap-my)*(xmap-mx)*corr/(sx*sy) + z*corr^2./correx )

pder[*,5]= ar*E/(2*!pi*sy^2.*sx^2.*correx^1.5)* $
( corr + (xmap-mx)*(ymap-my)/(sx*sy) - z*corr/correx )

pder[*,6]= 1
endif

end;pro gfunc2dgbis
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc2dg, X, A, F,pder,dummy=dummy
ar=A[0]
mx=A[1]
my=A[2]
sx=A[3]
sy=A[4]
sxy=A[5]
bkg0=A[6]

corr=sxy/(sx*sy)
correx=1.-corr^2.
xx=X.x
yy=X.y
aybin=n_elements(yy)
axbin=n_elements(xx)
xmap=(xx)#replicate(1,aybin)
ymap=replicate(1,axbin)#(yy)
xmap=reform(xmap,axbin*aybin,/overwrite)
ymap=reform(ymap,axbin*aybin,/overwrite)
z=(xmap-mx)^2./sx^2.-2.*corr*(xmap-mx)*(ymap-my)/sx/sy+(ymap-my)^2./sy^2.
E=exp(-z/(2.*correx))
F=ar/(2*!pi*sx*sy*sqrt(correx))*E+bkg0

if N_PARAMS() GE 4 then begin
pder=dblarr(axbin*aybin,n_elements(A))
;each row is a partial derivative with respect to one of the parameter

pder[*,0]= E/(2*!pi*sx*sy*sqrt(correx))

pder[*,1]= ar*E/(2*!pi*sy*sx^2.*correx^1.5)* $
( (xmap-mx)/sx - corr*(ymap-my)/sy )

pder[*,2]= ar*E/(2*!pi*sx*sy^2.*correx^1.5)* $
( (ymap-my)/sy - corr*(xmap-mx)/sx )

pder[*,3]= ar*E/(2*!pi*sy*sx^2*correx^1.5)* $
(-1. + (xmap-mx)^2./sx^2. - (xmap-mx)*(ymap-my)*corr/(sy*sx) + z*corr^2./correx )

pder[*,4]= ar*E/(2*!pi*sx*sy^2*correx^1.5)* $
(-1. + (ymap-my)^2./sy^2. - (ymap-my)*(xmap-mx)*corr/(sx*sy) + z*corr^2./correx )

pder[*,5]= ar*E/(2*!pi*sy^2.*sx^2.*correx^1.5)* $
( corr + (xmap-mx)*(ymap-my)/(sx*sy) - z*corr/correx )

pder[*,6]= 1
endif

end;pro gfunc2dg
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc2dPVc, X, A, F, struc ,pder
apeaks=n_elements(a)/7
cj=struc.j

; calculate physical meaningful parameters AA
; A is fit parameter
AA=reform(A,7,apeaks)
AAA=AA;=> fit parameters
;ar0,Car,mx0,Cmx,my0,Cmy,sx0,Csx,sy0,Csy,sxy0,Csxy,n0,Cn
constrar=struc.constrar
for i=0l,n_elements(cj)-1 do $
AA[cj[i],*]=constrar[2*cj[i],*]+2.*constrar[2*cj[i]+1,*]/!dpi*atan(AA[cj[i],*])
AA=reform(AA,7*apeaks);=> physical meaningful parameters

indar=indgen(apeaks)*7
indmx=indar+1
indmy=indar+2
indsx=indar+3
indsy=indar+4
indsxy=indar+5
indn=indar+6
ar=AA[indar]
mx=AA[indmx]
my=AA[indmy]
sx=AA[indsx]
sy=AA[indsy]
sxy=AA[indsxy]
n=AA[indn]

corr=sxy/(sx*sy)
correx=1.-corr^2.
xx=X.x
yy=X.y
aybin=n_elements(yy)
axbin=n_elements(xx)
xmap=(xx)#replicate(1,aybin)
ymap=replicate(1,axbin)#(yy)
xmap=reform(xmap,axbin*aybin,/overwrite)
ymap=reform(ymap,axbin*aybin,/overwrite)

F=xmap*0.
if N_PARAMS() GE 5 then $
pder=dblarr(axbin*aybin,n_elements(A))
;each row is a partial derivative with respect to one of the parameter

for k=0l,apeaks-1 do begin
z=(xmap-mx[k])^2./sx[k]^2.-2.*corr[k]*(xmap-mx[k])*(ymap-my[k])/sx[k]/sy[k]+(ymap-my[k])^2./sy[k]^2.
E=exp(-z/(2.*correx[k]))
EE=1/(z/correx[k]+0.25)
commF=ar[k]/(2*!dpi*sx[k]*sy[k]*sqrt(correx[k]))
G=E
L=EE/(2*!dpi)
F=F+((1-n[k])*G+n[k]*L)*commF

if N_PARAMS() GE 5 then begin
j=7*k

dAAdA=fltarr(7)+1
dAAdA[cj]=2/!dpi*constrar[cj*2+1,k]/(1+AAA[cj,k]^2.)

;A
commpder=1./(2*!dpi*sx[k]*sy[k]*sqrt(correx[k]))
Gpder=E
Lpder=EE/(2*!dpi)
pder[*,j]=dAAdA[0]*commpder*((1-n[k])*Gpder+n[k]*Lpder)

;mx
commpder=ar[k]/(2*!dpi*sy[k]*sx[k]^2.*correx[k]^1.5)* $
( (xmap-mx[k])/sx[k] - corr[k]*(ymap-my[k])/sy[k] )
Gpder=E
Lpder=EE^2./!dpi
pder[*,j+1]=dAAdA[1]*commpder*((1-n[k])*Gpder+n[k]*Lpder)

;my
commpder=ar[k]/(2*!dpi*sx[k]*sy[k]^2.*correx[k]^1.5)* $
( (ymap-my[k])/sy[k] - corr[k]*(xmap-mx[k])/sx[k] )
Gpder=E
Lpder=EE^2./!dpi
pder[*,j+2]=dAAdA[2]*commpder*((1-n[k])*Gpder+n[k]*Lpder)

;sx
commpder=ar[k]/(2*!dpi*sy[k]*sx[k]^2*correx[k]^1.5)
comm=(xmap-mx[k])^2./sx[k]^2. - (xmap-mx[k])*(ymap-my[k])*corr[k]/(sy[k]*sx[k]) + z*corr[k]^2./correx[k]
Gpder=E* (-1. + comm )
Lpder=EE/(2*!dpi)* (-1. + 2*EE*comm)
pder[*,j+3]=dAAdA[3]*commpder*((1-n[k])*Gpder+n[k]*Lpder)

;sy
commpder=ar[k]/(2*!dpi*sx[k]*sy[k]^2*correx[k]^1.5)
comm=(ymap-my[k])^2./sy[k]^2. - (ymap-my[k])*(xmap-mx[k])*corr[k]/(sx[k]*sy[k]) + z*corr[k]^2./correx[k]
Gpder=E* (-1. + comm )
Lpder=EE/(2*!dpi)* (-1. + 2*EE*comm )
pder[*,j+4]=dAAdA[4]*commpder*((1-n[k])*Gpder+n[k]*Lpder)

;sxy
commpder=ar[k]/(2*!dpi*sy[k]^2.*sx[k]^2.*correx[k]^1.5)
comm=(xmap-mx[k])*(ymap-my[k])/(sx[k]*sy[k]) - z*corr[k]/correx[k]
Gpder=E* ( corr[k] + comm )
Lpder=EE/(2*!dpi)* ( corr[k] + 2*EE*comm )
pder[*,j+5]=dAAdA[5]*commpder*((1-n[k])*Gpder+n[k]*Lpder)

;n
pder[*,j+6]=dAAdA[6]*(L-G)*commF
endif

endfor

end; pro gfunc2dPVc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfunc2dPV, X, A, F,pder
apeaks=n_elements(a)/7
indar=indgen(apeaks)*7
indmx=indar+1
indmy=indar+2
indsx=indar+3
indsy=indar+4
indsxy=indar+5
indn=indar+6
ar=A[indar]
mx=A[indmx]
my=A[indmy]
sx=A[indsx]
sy=A[indsy]
sxy=A[indsxy]
n=A[indn]

corr=sxy/(sx*sy)
correx=1.-corr^2.
xx=X.x
yy=X.y
aybin=n_elements(yy)
axbin=n_elements(xx)
xmap=(xx)#replicate(1,aybin)
ymap=replicate(1,axbin)#(yy)
xmap=reform(xmap,axbin*aybin,/overwrite)
ymap=reform(ymap,axbin*aybin,/overwrite)

F=xmap*0.
if N_PARAMS() GE 4 then $
pder=dblarr(axbin*aybin,n_elements(A))
;each row is a partial derivative with respect to one of the parameter

for k=0l,apeaks-1 do begin
z=(xmap-mx[k])^2./sx[k]^2.-2.*corr[k]*(xmap-mx[k])*(ymap-my[k])/sx[k]/sy[k]+(ymap-my[k])^2./sy[k]^2.
E=exp(-z/(2.*correx[k]))
EE=1/(z/correx[k]+0.25)
commF=ar[k]/(2*!dpi*sx[k]*sy[k]*sqrt(correx[k]))
G=E
L=EE/(2*!dpi)
F=F+((1-n[k])*G+n[k]*L)*commF

if N_PARAMS() GE 4 then begin
j=7*k

commpder=1./(2*!dpi*sx[k]*sy[k]*sqrt(correx[k]))
Gpder=E
Lpder=EE/(2*!dpi)
pder[*,j]=commpder*((1-n[k])*Gpder+n[k]*Lpder)

commpder=ar[k]/(2*!dpi*sy[k]*sx[k]^2.*correx[k]^1.5)* $
( (xmap-mx[k])/sx[k] - corr[k]*(ymap-my[k])/sy[k] )
Gpder=E
Lpder=EE^2./!dpi
pder[*,j+1]=commpder*((1-n[k])*Gpder+n[k]*Lpder)

commpder=ar[k]/(2*!dpi*sx[k]*sy[k]^2.*correx[k]^1.5)* $
( (ymap-my[k])/sy[k] - corr[k]*(xmap-mx[k])/sx[k] )
Gpder=E
Lpder=EE^2./!dpi
pder[*,j+2]=commpder*((1-n[k])*Gpder+n[k]*Lpder)

commpder=ar[k]/(2*!dpi*sy[k]*sx[k]^2*correx[k]^1.5)
comm=(xmap-mx[k])^2./sx[k]^2. - (xmap-mx[k])*(ymap-my[k])*corr[k]/(sy[k]*sx[k]) + z*corr[k]^2./correx[k]
Gpder=E* (-1. + comm )
Lpder=EE/(2*!dpi)* (-1. + 2*EE*comm)
pder[*,j+3]=commpder*((1-n[k])*Gpder+n[k]*Lpder)

commpder=ar[k]/(2*!dpi*sx[k]*sy[k]^2*correx[k]^1.5)
comm=(ymap-my[k])^2./sy[k]^2. - (ymap-my[k])*(xmap-mx[k])*corr[k]/(sx[k]*sy[k]) + z*corr[k]^2./correx[k]
Gpder=E* (-1. + comm )
Lpder=EE/(2*!dpi)* (-1. + 2*EE*comm )
pder[*,j+4]=commpder*((1-n[k])*Gpder+n[k]*Lpder)

commpder=ar[k]/(2*!dpi*sy[k]^2.*sx[k]^2.*correx[k]^1.5)
comm=(xmap-mx[k])*(ymap-my[k])/(sx[k]*sy[k]) - z*corr[k]/correx[k]
Gpder=E* ( corr[k] + comm )
Lpder=EE/(2*!dpi)* ( corr[k] + 2*EE*comm )
pder[*,j+5]=commpder*((1-n[k])*Gpder+n[k]*Lpder)

pder[*,j+6]=(L-G)*commF
endif

endfor

end; pro gfunc2dPV
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function marquardt,x, y, weights, a, sigma, FUNCTION_NAME = Function_Name, $
                        ITMAX=itmax, ITER=iter, TOL=tol, $
                        CHISQ=chisq,STRUC=STRUC,nfree=nfree,silent=silent,CError=CError
;Modified curvefit

;                Input                                Output
; Variables
; x:             un depending variable                no change
; y:             depending variable                    no change
; weights:         usually 1/y                            no change
; a:             initial guess of parameters         fitted parameters
; sigma:        none                                standard deviations
; Keywords
; FUNCTION_NAME    name procedure (calc pder and yfit)    no change
; ITMAX            max number of iterations allowed    no change
; ITER            none                                number of iter performed
; TOL            Convergence tolerance                no change
; CHISQ            none                                Reduced Chisq of fit
; STRUC            structure of gfunc input            no change
; NFREE            none                                degrees of freedom
; SILENT        set/not set                            none
; CERROR        none                                convertion error code

; ----Initialize----
if n_elements(function_name) LE 0 then function_name = "FUNCT";procedure that calculates yfit and pder
if n_elements(tol) EQ 0 then tol = 1.e-3      ;Convergence tolerance
if n_elements(itmax) EQ 0 then itmax = 20     ;Maximum # iterations
nterms = n_elements(a)         ; # of parameters
nfree = n_elements(y) - nterms ; Degrees of freedom
if (nfree LE 0) and not(keyword_set(silent)) then message, 'Marquardt - not enough data points.'
flambda = 0.001                   ;Initial lambda
diag = lindgen(nterms)*(nterms+1) ; Subscripts of diagonal elements
CError=0B

; ----Outher loop----
FOR iter = 1, itmax DO begin
    ; ----From parameters -> yfit+pder -> chisq,alpha,beta,sigma----
    if n_elements(struc) ne 0 then call_procedure, function_name, x, a, yfit,struc, pder else $
    call_procedure, function_name, x, a, yfit, pder
    if nterms EQ 1 then pder = reform(pder, n_elements(y), 1);eg Array[5]->Array[5,1]
    beta = (y-yfit)*Weights # pder
    alpha = transpose(pder) # (Weights # (fltarr(nterms)+1)*pder)
    sigmaold = sqrt( 1.0 / alpha[diag] )
    sigma  = sigmaold
    chisqold = total(Weights*(y-yfit)^2)/nfree
    chisq = chisqold
    yfitold = yfit
    done_early = chisqold LT total(abs(y))/1e7/NFREE
    if done_early then GOTO, done
    c = sqrt(alpha[diag])
    c = c # c
    lambdaCount = 0

    ; ----Inner loop----
    REPEAT begin
        ; ----Increase lambda and calculate parameters----
        ; ----Do this until chisq start decreasing----
        lambdaCount = lambdaCount + 1
        array = alpha / c ;Normalize alpha to have unit diagonal
        array[diag] = array[diag]*(1.+flambda)
        if n_elements(array) EQ 1 then array = (1.0 / array) $
        else array = invert(array)
        b = a + array/c # transpose(beta);parameter increments
        if n_elements(struc) ne 0 then call_procedure, function_name, x, b, yfit,struc else $  ; Evaluate function,with new param
        call_procedure, function_name, x, b, yfit
        chisq = total(Weights*(y-yfit)^2)/nfree    ; New chisq
        sigma = sqrt(array[diag]/alpha[diag])      ; New sigma
        if (finite(chisq) EQ 0) OR $
           (lambdaCount GT 30 AND chisq GE chisqold) then begin
           ; Reject changes made this iteration, use old values.
           yfit  = yfitold
           sigma = sigmaold
           chisq = chisqold
           CError=1B
           if not(keyword_set(silent)) then message, 'Failed to converge while increasing lambda.', /INFORMATIONAL
           GOTO, done
        endif
        flambda = flambda*10. ; Assume fit got worse
    ENDREP UNTIL chisq LE chisqold ; Inner loop

    flambda = flambda/100.
    a=b                                    ; Save new parameter estimate.
    if ((chisqold-chisq)/chisqold) LE tol then GOTO,done   ;Finished?
ENDFOR ; Outer loop

CError=1B
if not(keyword_set(silent)) then MESSAGE,'Failed to converge while decreasing lambda.', /INFORMATIONAL
done:  if done_early then iter = iter - 1
       return,yfit          ; return result
END ;function marquardt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MarquardtC,x, y, weights, a, sigma, FUNCTION_NAME = Function_Name, $
                        ITMAX=itmax, ITER=iter, TOL=tol, $
                        CHISQ=chisq,STRUC=STRUC,nfree=nfree,constr=constr,$
                        silent=silent,CError=CError,NUMDER=numder, CHECKROUTINE = checkroutine
;Modified curvefit

;                Input                                Output
; Variables
; x:             un depending variable                no change
; y:             depending variable                    no change
; weights:         usually 1/y                            no change
; a:             initial guess of parameters         fitted parameters
; sigma:        none                                standard deviations
; Keywords
; FUNCTION_NAME    name procedure (calc pder and yfit)    no change
; ITMAX            max number of iterations allowed    no change
; ITER            none                                number of iter performed
; TOL            Convergence tolerance                no change
; CHISQ            none                                Reduced Chisq of fit
; STRUC            structure of gfunc input            no change
; NFREE            none                                degrees of freedom
; CONSTR        restraints (no restr: 0)            no change
; SILENT        none                                none
; CERROR        none                                convertion error code
; NUMDER        set or not                            none

; ----Initialize----
CError=0B
if n_elements(function_name) LE 0 then function_name = "FUNCT" ;procedure that calculates yfit and pder
if n_elements(tol) EQ 0 then tol = 1.e-3      ;Convergence tolerance
if n_elements(itmax) EQ 0 then itmax = 20     ;Maximum # iterations
nterms = n_elements(a)         ; # of parameters
nY=n_elements(y)
nfree = nY - nterms ; Degrees of freedom
if (nfree LE 0) and not(keyword_set(silent)) then begin
    message, 'Marquardt - not enough data points.', /INFORMATIONAL
    return,-1
endif
if n_elements(constr) EQ 0 then constr=a*0
flambda = 0.001                   ;Initial lambda
diag = lindgen(nterms)*(nterms+1) ; Subscripts of diagonal elements
bnumder=keyword_set(numder)
if bnumder then begin
    res = machar(DOUBLE=double)
    eps = sqrt(res.eps)
    pder = fltarr(nY, nterms)
endif

; ----Prepare for restraints----
ao=a
ac=constr
aci=1./ac^2
ind=where(constr eq 0,ct)
if ct ne 0 then aci[ind]=0

; ----Outher loop----
FOR iter = 1, itmax DO begin
    ; ----From parameters -> yfit+pder -> chisq,alpha,beta,sigma----
    if bnumder then begin
        if n_elements(struc) ne 0 then call_procedure, function_name, x, a, yfit,struc,pdercalc else $ ;,pdercalc
        call_procedure, function_name, x, a, yfit
        yfit=float(yfit)

        FOR term=0, nterms-1 DO begin
              p = a       ; Copy current parameters

              ; Increment size for forward difference derivative
              inc = eps * abs(p[term])
              if (inc eq 0) then inc = eps
              p[term] += inc

              if n_elements(struc) ne 0 then call_procedure, function_name, x, p, yfit1,struc else $
              call_procedure, function_name, x, p, yfit1

              pder[*,term] = (float(yfit1)-yfit)/inc
        ENDFOR

        ; Compare calculated and numeric pder (use pdercalc in first call_procedure!!!)
;;        xrange=[2200,2400]
;        window
;        for term=0,(size(pder,/dimensions))[1]-1 do begin
;            plot,pder[*,term],xrange=xrange,psym=4 & oplot,pdercalc[*,term]
;            ;plot,pder[*,term]-pdercalc[*,term],xrange=xrange
;            wait,0.5
;        endfor
;        stop

    endif else begin
        if n_elements(struc) ne 0 then call_procedure, function_name, x, a, yfit,struc, pder else $
        call_procedure, function_name, x, a, yfit, pder
        yfit=float(yfit)
        pder=float(pder)
    endelse

    if nterms EQ 1 then pder = reform(pder, nY, 1);eg Array[1,5]->Array[5,1]
    beta = (y-yfit)*Weights # pder + ((ao-a)*aci)
    alpha = transpose(pder) # (Weights # (fltarr(nterms)+1)*pder)
    alpha[diag]+=aci
    sigmaold = sqrt( 1.0 / alpha[diag] )
    sigma  = sigmaold
    chisqold = (total(Weights*(y-yfit)^2)+ total((ao-a)^2*aci))/nfree
    chisq = chisqold
    yfitold = yfit
    done_early = chisqold LT total(abs(y))/1e7/NFREE
    if done_early then GOTO, done
    c = sqrt(alpha[diag])
    c = c # c
    lambdaCount = 0

    ; ----Inner loop----
    REPEAT begin
        ; ----Increase lambda and calculate parameters----
        ; ----Do this until chisq start decreasing----
        lambdaCount++
        array = alpha / c ;Normalize alpha to have unit diagonal
        array[diag] *= (1.+flambda)
        if n_elements(array) EQ 1 then array = (1.0 / array) $
        else array = invert(array)
        b = a + array/c # transpose(beta);parameter increments
        if n_elements(struc) ne 0 then call_procedure, function_name, x, b, yfit,struc else $  ; Evaluate function,with new param
        call_procedure, function_name, x, b, yfit
        yfit=float(yfit)
        chisq = (total(Weights*(y-yfit)^2) +total((ao-b)^2*aci))/nfree ; New chisq
        sigma = sqrt(array[diag]/alpha[diag])                          ; New sigma
        if (finite(chisq) EQ 0) OR $
           (lambdaCount GT 30 AND chisq GE chisqold) then begin
           ; Reject changes made this iteration, use old values.
           yfit  = yfitold
           sigma = sigmaold
           chisq = chisqold
           CError=1B
           if not(keyword_set(silent)) then message, 'Failed to converge while increasing lambda.', /INFORMATIONAL
           GOTO, done
        endif
        flambda *= 10. ; Assume fit got worse
    ENDREP UNTIL chisq LE chisqold ; Inner loop

    flambda /= 100.
    a=b                                    ; Save new parameter estimate.
    if ((chisqold-chisq)/chisqold) LE tol then GOTO,done   ;Finished?
ENDFOR ; Outer loop

if not(keyword_set(silent)) then MESSAGE, 'Failed to converge while decreasing lambda.', /INFORMATIONAL
done:  if done_early then iter--

CError=1B
return,yfit          ; return result
END ;function MarquardtC
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;pro chiplot,ev
;ID=widget_info(ev.top,FIND_BY_UNAME='draw')
;widget_control,ID,get_uvalue=temp
;zs=temp.zs
;widget_control,ev.top,get_uvalue=list
;
;var1=list.var1
;var2=list.var2
;z=list.z
;img=list.img
;x=list.x
;y=list.y
;img=list.img
;
;erase
;ny=n_elements(y)
;nx=n_elements(x)
;x=x#(intarr(ny)+1)
;y=(intarr(nx)+1)#y
;shade_surf,img,x,y,/t3d,$
;    xtitle='VAR1',ytitle='VAR2',ztitle='Chisq.Red.',$
;    zrange=[min(img),(zs>0.1)*max(img)]
;plots,var1,var2,z,/data,/t3d,psym=-2
;loadct,3,/silent
;plots,var1[n_elements(var1)-1],var2[n_elements(var2)-1],z[n_elements(z)-1],/data,/t3d,psym=2,color=16*14
;loadct,0,/silent
;end;pro chiplot
;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ChiConstruct,x,y,weights,nfree,a,ao,aci,struc,function_name,nvar,varnr,var1in,zin,var2in
if n_elements(zin) le 3 then return
; ----Display chisq search----
case nvar of
1:    begin
    var1=var1in[1:*]
    z=zin[1:*]
    ; Make Plot
    xmap=MakeRange(min(var1),max(var1),incper=1.,rangeext=50.)
    nx=n_elements(xmap)
    img=fltarr(nx)
    for i=0l,nx-1 do begin
        var=a
        var[varnr]=xmap[i]
        if n_elements(struc) ne 0 then call_procedure, function_name, x, var, yfit2,struc else $
           call_procedure, function_name, x, var, yfit2
           img[i] = total(Weights*(y-yfit2)^2)/nfree +total((ao-var)^2*aci)
    endfor
    window
    plot,xmap,img
    plots,var1,z,psym=-2
    loadct,3,/silent
    plots,var1[n_elements(var1)-1],z[n_elements(z)-1],psym=-2,color=10*16
    loadct,0,/silent
    endcase
2:    begin
    var1=var1in[1:*]
    var2=var2in[1:*]
    z=zin[1:*]
    ; Make Map
    vmin=min(var1)
    vmax=max(var1)
    if float(vmin) eq float(vmax) then begin
        vmin=vmin-0.01*vmin
        vmax=vmax+0.01*vmax
    endif
    xmap=MakeRange(vmin,vmax,incper=2.,rangeext=50.)
    vmin=min(var2)
    vmax=max(var2)
    if float(vmin) eq float(vmax) then begin
        vmin=vmin-0.01*vmin
        vmax=vmax+0.01*vmax
    endif
    ymap=MakeRange(vmin,vmax,incper=2.,rangeext=50.)
    nx=n_elements(xmap)
    ny=n_elements(ymap)
    img=fltarr(nx,ny)
    for i=0l,nx-1 do $
        for j=0l,ny-1 do begin
            var=a
            var[varnr]=[xmap[i],ymap[j]]
            if n_elements(struc) ne 0 then call_procedure, function_name, x, var, yfit2,struc else $
            call_procedure, function_name, x, var, yfit2
            img[i,j] = total(Weights*(y-yfit2)^2)/nfree +total((ao-var)^2*aci)
        endfor
;    struct={var1:var1,var2:var2,z:z,img:img,x:xmap,y:ymap,proplot:'chiplot'}
;    flplot,struct
    ;@@@@@@
    ; Plot sections of chisq-surface with a plane trough optimal minimum
    img=fltarr(nx)
    for j=0l,nx-1 do begin
        var=a
        var[varnr[0]]=xmap[j]
        if n_elements(struc) ne 0 then call_procedure, function_name, x, var, yfit2,struc else $
           call_procedure, function_name, x, var, yfit2
           img[j] = total(Weights*(y-yfit2)^2)/nfree +total((ao-var)^2*aci)
    endfor
    window,0
    plot,xmap,img,xtitle='Var1',ytitle='chisqr',ystyle=1
    plots,[a[varnr[0]],a[varnr[0]]],[!Y.CRange[0], !Y.CRange[1]],linestyle=1,NOCLIP = 0
    img=fltarr(ny)
    for j=0l,ny-1 do begin
        var=a
        var[varnr[1]]=ymap[j]
        if n_elements(struc) ne 0 then call_procedure, function_name, x, var, yfit2,struc else $
           call_procedure, function_name, x, var, yfit2
           img[j] = total(Weights*(y-yfit2)^2)/nfree +total((ao-var)^2*aci)
    endfor
    window,1
    plot,ymap,img,xtitle='Var2',ytitle='chisqr',ystyle=1
    plots,[a[varnr[1]],a[varnr[1]]],[!Y.CRange[0], !Y.CRange[1]],linestyle=1,NOCLIP = 0
    temp=''
    ; Print absolute difference in terms when delta(b)
    ;call_procedure, function_name, x, [a[0],a[1]], yfit2,struc,out=out1
    ;call_procedure, function_name, x, [a[0],a[1]+0.2], yfit2,struc,out=out2
    ;print,out2-out1
    ; Display partial derivative
    ;tiff=fltarr(max(x[*,0])+1,max(x[*,1])+1)
    ;call_procedure, function_name, x, ao, yfit2,struc,pder,out=out
    ;tiff[x[*,0],x[*,1]]=abs(pder[*,0])
    ;window
    ;shade_surf,tiff
    ;@@@@@@
    endcase
endcase
end;pro ChiConstruct
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MarquardtMap,x, y, weights, a, sigma, FUNCTION_NAME = Function_Name, $
                        ITMAX=itmax, ITER=iter, TOL=tol, CHISQ=chisq,$
                        STRUC=STRUC,nfree=nfree,constr=constr,varnr=varnr,silent=silent
;Modified curvefit
; Purpose: Construct 'searchpath' on Chisq-surface
; Remark: This representation will only be correct when
; one or two parameters are present. The Chisq-surface
; is constructed in an area where the targetted parameters
; stayed in during search. The not-targeted parameters
; are used also for Chisq-surface construction: we will
; use there last refined values before 'ChiConstruct' call.
; This way the 'searchpath' may well be NOT on the
; Chisq surface (except for the last point of course)!!!!!!
; (unless only the non-targetted parameters don't change a lot)
; Also the minimum in the displayed Chisq-surface may
; not be the minimum of the fit when more then two parameters
; are present.

;                Input                                Output
; Variables
; X:             un depending variable                no change
; Y:             depending variable                    no change
; weights:         usually 1/y                            no change
; A:             initial guess of parameters         fitted parameters
; Keywords
; FUNCTION_NAME    name procedure (calc pder and yfit)    no change
; ITMAX            max number of iterations allowed    no change
; ITER            none                                number of iter performed
; TOL            Convergence tolerance                no change
; CHISQ            none                                Chisq of fit
; STRUC            structure of gfunc input            no change
; NFREE            none                                degrees of freedom
; CONSTR        restraints (no constr: -1)            no change
; VARNR            Variables to watch                    no change
; SILENT        none                                none

; ----Initialize----
if n_elements(function_name) LE 0 then function_name = "FUNCT" ;procedure that calculates yfit and pder
if n_elements(tol) EQ 0 then tol = 1.e-3      ;Convergence tolerance
if n_elements(itmax) EQ 0 then itmax = 20     ;Maximum # iterations
nterms = n_elements(a)         ; # of parameters
nfree = n_elements(y) - nterms ; Degrees of freedom
if (nfree LE 0) and not(keyword_set(silent)) then message, 'Marquardt - not enough data points.'
flambda = 0.001                   ;Initial lambda
diag = lindgen(nterms)*(nterms+1) ; Subscripts of diagonal elements
if n_elements(constr) eq 0 then constr=a*0.-1 ;Fit without restraints
nvar=n_elements(varnr) ;Prepare storage of Search path in Chisq-space
case nvar of
1:     var1=-1.
2:    begin
    var1=-1.
    if nterms le 1 then begin
        varnr=varnr[0]
        nvar=1
    endif else var2=-1.
    endcase
else: begin
    if nterms le 1 then begin
        varnr=[0]
        var1=-1.
        nvar=1
    endif else begin
        varnr=[0,1]
        var1=-1.
        var2=-1.
        nvar=2
    endelse
    endcase
endcase
z=-1.


; ----Prepare for restraints----
ao=a
ac=constr
aci=1./ac^2
ind=where(constr eq 0,ct)
if (ct ne 0) and not(keyword_set(silent)) then message, 'Marquardt - restraints can`t be zero.'
ind=where(constr eq -1,ct)
if ct ne 0 then aci[ind]=0

; ----Outher loop----
FOR iter = 1, itmax DO begin
    ; ----From parameters -> yfit -> chisq,alpha,beta,sigma----
    if n_elements(struc) ne 0 then call_procedure, function_name, x, a, yfit,struc, pder else $
    call_procedure, function_name, x, a, yfit, pder

    var1=[var1,a[varnr[0]]];@@@@@@@@@@@@@@@
    if nvar eq 2 then var2=[var2,a[varnr[1]]];@@@@@@@@@@@@@@@

    if nterms EQ 1 then pder = reform(pder, n_elements(y), 1);eg Array[1,5]->Array[5,1]
    beta = (y-yfit)*Weights # pder + ((ao-a)*aci)
    alpha = transpose(pder) # (Weights # (fltarr(nterms)+1)*pder)
    alpha[diag]=alpha[diag]+aci
    sigmaold = sqrt( 1.0 / alpha[diag] )
    sigma  = sigmaold
    chisqold = total(Weights*(y-yfit)^2)/nfree + total((ao-a)^2*aci)
    chisq = chisqold

    z=[z,chisq];@@@@@@@@@@@@@@@
    ;ChiConstruct,x,y,weights,nfree,a,ao,aci,struc,function_name,nvar,varnr,var1,z,var2

    yfitold = yfit
    done_early = chisqold LT total(abs(y))/1e7/NFREE
    if done_early then GOTO, done
    c = sqrt(alpha[diag])
    c = c # c
    lambdaCount = 0

    ; ----Inner loop----
    REPEAT begin
        ; ----Increase lambda and calculate parameters----
        ; ----Do this until chisq starts decreasing----
        lambdaCount = lambdaCount + 1
        array = alpha / c ;Normalize alpha to have unit diagonal
        array[diag] = array[diag]*(1.+flambda)
        if n_elements(array) EQ 1 then array = (1.0 / array) $
        else array = invert(array)
        b = a + array/c # transpose(beta);parameter increments

        var1=[var1,b[varnr[0]]];@@@@@@@@@@@@@@@
        if nvar eq 2 then var2=[var2,b[varnr[1]]];@@@@@@@@@@@@@@@

        if n_elements(struc) ne 0 then call_procedure, function_name, x, b, yfit,struc else $  ; Evaluate function,with new param
        call_procedure, function_name, x, b, yfit
        chisq = total(Weights*(y-yfit)^2)/nfree +total((ao-b)^2*aci) ; New chisq

        z=[z,chisq];@@@@@@@@@@@@@@@
        ;ChiConstruct,x,y,weights,nfree,b,ao,aci,struc,function_name,nvar,varnr,var1,z,var2

        sigma = sqrt(array[diag]/alpha[diag])                          ; New sigma
        if (finite(chisq) EQ 0) OR $
           (lambdaCount GT 30 AND chisq GE chisqold) then begin
           ; Reject changes made this iteration, use old values.
           yfit  = yfitold
           sigma = sigmaold
           chisq = chisqold
           if not(keyword_set(silent)) then message, 'Failed to converge while increasing lambda.', /INFORMATIONAL
           GOTO, done
        endif
        flambda = flambda*10. ; Assume fit got worse
    ENDREP UNTIL chisq LE chisqold ; Inner loop

    flambda = flambda/100.
    a=b                                    ; Save new parameter estimate.
    if ((chisqold-chisq)/chisqold) LE tol then GOTO,done   ;Finished?
ENDFOR ; Outer loop

if not(keyword_set(silent)) then MESSAGE, 'Failed to converge while decreasing lambda.', /INFORMATIONAL
done:      if done_early then iter = iter - 1
        ChiConstruct,x,y,weights,nfree,a,ao,aci,struc,function_name,nvar,varnr,var1,z,var2
        return,yfit          ; return result
END ;function MarquardtMap
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NLLS_norm,vec,machvals

nvec=n_elements(vec)
agiant = machvals.rgiant / nvec
adwarf = machvals.rdwarf * nvec

mx = max(vec, min=mn)
mx = max(abs([mx,mn]))
if mx EQ 0 then return, vec[0]*0.

if mx GT agiant OR mx LT adwarf then return,mx * sqrt(total((vec/mx)^2))$
else                                 return,sqrt( total(vec^2) )
  
end;function NLLS_norm
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NLLS_LMpar, r, ipvt, diag, qtb, delta, x, sdiag, machvals, par=par
; Solve this:
;    R.x = Q^T.wResid AND sqrt(par).D.x = 0
;    
; NLLS_LMpar input
;     r: upper triangle of R
;     ipvt: pivoting array
;     diag: diagonal terms of D
;     qtb: Q^T.wResid (first n terms)
;     delta: upperbound on dxnorm
;     par: initial estimate
; 
; NLLS_LMpar output
;     P^T(R^T.R+par.D.D).P=S^T.S
;     
;     r: lower triangle contains strict upper triangle of S
;     x: least squares solution x of the outer+inner system of equations
;     sdiag: diagonal terms of S (not used)
;     par: final estimate
        
sz = size(r)
ncol = sz[1]
nrow = sz[2]
delm = lindgen(nrow) * (ncol+1) ;; Diagonal elements of r

; Compute and store in x the gauss-newton direction.  If the
; jacobian is rank-deficient, obtain a least-squares solution
nsing = nrow
wa1 = qtb
wh = where(r[delm] EQ 0, ct)
if ct GT 0 then begin
    nsing = wh[0]
    wa1[wh[0]:*] = 0
endif

if nsing GE 1 then begin
    for j=nsing-1,0,-1 do begin
        wa1[j] /= r[j,j]
        if (j-1 GE 0) then $
            wa1[0:(j-1)] -= r[0:(j-1),j]*wa1[j]
    endfor
endif

; Note: ipvt here is a permutation array
x[ipvt] = wa1

; Initialize the iteration counter.  Evaluate the function at the
; origin, and test for acceptance of the gauss-newton direction
iter = 0L
wa2 = diag * x
dxnorm = NLLS_norm(wa2,machvals)
fp = dxnorm - delta
if fp gt 0.1*delta then begin

    ; If the jacobian is not rank deficient, the newton step provides a
    ; lower bound, parl, for the zero of the function.  Otherwise set
    ; this bound to zero.

    zero = wa2[0]*0.
    parl = zero
    if nsing GE nrow then begin
        wa1 = diag[ipvt]*wa2[ipvt]/dxnorm

        wa1[0] /= r[0,0] ;; Degenerate case
        for j=1L, nrow-1 do begin   ;; Note "1" here, not zero
            sum = total(r[0:(j-1),j]*wa1[0:(j-1)])
            wa1[j] = (wa1[j] - sum)/r[j,j]
        endfor

        temp = NLLS_norm(wa1,machvals)
        parl = ((fp/delta)/temp)/temp
    endif

    ; Calculate an upper bound, paru, for the zero of the function
    for j=0l, nrow-1 do begin
          sum = total(r[0:j,j]*qtb[0:j])
          wa1[j] = sum/diag[ipvt[j]]
    endfor
    gnorm = NLLS_norm(wa1,machvals)
    paru  = gnorm/delta
    if paru EQ 0 then paru = machvals.DWARF/(delta<0.1)

    ; If the input par lies outside of the interval (parl,paru), set
    ; par to the closer endpoint
    par= parl>par<paru
    if par EQ 0 then par = gnorm/dxnorm

    ; Beginning of an interation
    repeat begin
        iter++

        ; Evaluate the function at the current value of par
        if par EQ 0 then par = max([machvals.DWARF, paru*0.001])
        temp = sqrt(par)
        wa1 = temp * diag
        QRsolve, r, ipvt, wa1, qtb, x, sdiag
        wa2 = diag*x
        dxnorm = NLLS_norm(wa2,machvals)
        temp = fp
        fp = dxnorm - delta

        bterminate=(abs(fp) LE 0.1D*delta) $
            OR ((parl EQ 0) AND (fp LE temp) AND (temp LT 0)) $
            OR (iter EQ 10)

        if ~bterminate then begin
            ; Compute the newton correction
            wa1 = diag[ipvt]*wa2[ipvt]/dxnorm

            for j=0l,nrow-2 do begin
                wa1[j] = wa1[j]/sdiag[j]
                wa1[j+1:nrow-1] = wa1[j+1:nrow-1] - r[j+1:nrow-1,j]*wa1[j]
            endfor
            wa1[nrow-1] = wa1[nrow-1]/sdiag[nrow-1] ;; Degenerate case

            temp = NLLS_norm(wa1,machvals)
            parc = ((fp/delta)/temp)/temp

            ; Depending on the sign of the function, update parl or paru
            if fp GT 0 then parl = max([parl,par])
            if fp LT 0 then paru = min([paru,par])

            ; Compute an improved estimate for par
            par = max([parl, par+parc])
        endif
    endrep until bterminate
endif

if iter EQ 0 then return, par[0]*0.
return, par
end;function NLLS_LMpar
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro NLLS_tie,A,indtied,tied
A[indtied]=A[tied]
end;pro NLLS_tie
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakeAInfo,n
; limited: [low,high]
;    1 -> limit
;    2 -> don't limit during fit but reject afterwards
return,replicate({fixed:0b,nonls:0b,limited:bytarr(2),limits:fltarr(2),$
    limitedcombo:bytarr(2),comboID:0UL,combocoeff:0.,limitscombo:fltarr(2),$
    tied:-1L},n)
end;function MakeAInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NLLSPoissonW,y,type
ny=n_elements(y)

case type of
0:     return,make_array(ny,value=y[0]*0.+1.) ; no weighting
1:     return,sqrt(y>1) ; counts cut off at 1
2:     begin ; counts cut off at the lowest positive number
    ySTDEV=sqrt(y>0)
    temp1=where(ySTDEV eq 0,ct,comp=temp2)
    if ct ne 0 then $
        if ct eq ny then ySTDEV+=1 $ ; All zero
        else ySTDEV[temp1]=min(ySTDEV[temp2])
    endcase
3:     begin ; counts cut off at the lowest number
    ySTDEV=y
    mi=min(ySTDEV)
    if mi lt 0 then ySTDEV-=mi
    ySTDEV=sqrt(ySTDEV)
    temp1=where(ySTDEV eq 0,ct,comp=temp2)
    if ct ne 0 then $
        if ct eq ny then ySTDEV+=1 $ ; All zero
        else ySTDEV[temp1]=min(ySTDEV[temp2])
    endcase
else:return,make_array(ny,value=y[0]*0.+1.) ; no weighting
endcase

return,ySTDEV
end;function NLLSPoissonW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NLLS,xin, yin, A, fname,$
    Ainfo=Ainfo,ySTDEV=ySTDEVin,sigma=sigma,weight=weight,$
    ITMAX=itmax, TOL=tol,NUMDER=numder,$
    CHISQ=chisq, DOF=dof,  ITER=iter, Npegged=Npegged,$
    CError=CError,ErrorMsg=ErrorMsg,debug=debug, correlation=correlation,_REF_EXTRA=extra

;                Input                                Output
; Variables
; x*:             independent variable                no change
; y*:             dependent variable                    no change
; A*:             initial guess of parameters         fitted parameters
; fname*:        name procedure (calc pder and yfit)    no change
;        e.g. gfunc_test,X(independent variable),A(param),F(yfit),PDER(dyfit/dparam),struc=struc
;
; Keywords
; Ainfo:        array of structures                    no change
;         ainfo.fixed:        byte          Fixed=1b, free=0b
;         ainfo.nonls:        byte          This means: not refined by least squares. It is allowed for fname to
;                                           change these parameters to initialize them for the next iteration.
;                                           If there is no next iteration, this change will not be saved in the result.
;
;       Single value limits: L >= A <= H
;        ainfo.limited:      bytarr(2)     Use limit (Low/High)=1b, Don't use limit (Low/High)=0b
;        ainfo.limits:       fltarr(2)     Lower/Higher limit value
;
;       Linear combination limits: L >= A[0]+2*A[1] <= H
;        ainfo.limitedcombo: bytarr(2)     Use combo limit (Low/High)=1b, Don't use combo limit (Low/High)=0b
;        ainfo.comboID:      ulong         ID of linear combination  constraint (ID=0 -> not in any constraint)
;        ainfo.combocoeff:   float         Coefficient for linear combination
;        ainfo.limitscombo:  fltarr(2)     Lower/Higher limit value
;
;        ainfo.tied:         long          A[i]=A[tied]
;
; ySTDEV        y standard deviation                no change unless nothing was given
; sigma:        none                                standard deviations
; weight:        type of weighting                    none
; ITMAX            max number of iterations allowed    no change
; TOL            Convergence tolerance                no change
; NUMDER        set or not                            none
; CHISQ            none                                Chisq of fit
; DOF            none                                degrees of freedom
; ITER            none                                number of iter performed
; Npegged        none                                number of parameters pegged at there hard limits
; CERROR        none                                convertion error code
; ErrorMsg        none                                error message
; debug            set or not                            none
; correlation    set or not                            correlation info
; _REF_EXTRA    structure of gfunc input            no change
;
;(*=required)

CError=-1
ErrorMsg=''
chisq=0.
iter=0
Npegged=0

catch, catcherror
if catcherror NE 0 then begin
    catch, /cancel
    CError=-1
    ErrorMsg=!err_string
    return, -1
endif

; ----Check input----
nyorg=n_elements(yin)
if nyorg eq 0 or $ ;n_elements(x) ne ny or
    n_elements(A) eq 0 or n_elements(fname) eq 0 then begin
    ErrorMsg='Check input parameters.'
    return,-1
endif

; ----Skip missing data----
indvalid=where(finite(yin),nvalid)
if nvalid eq 0 then begin
    ErrorMsg='Invalid data (NaN or Inf).'
    return,-1
endif
y=yin[indvalid]
if n_elements(xin) mod nyorg eq 0 then begin
    s=size(xin,/dim)
    tmp=where(s eq nyorg,cttmp)
    case tmp of
    0: x=xin[indvalid,*,*,*,*,*,*,*]
    1: x=xin[*,indvalid,*,*,*,*,*,*]
    2: x=xin[*,*,indvalid,*,*,*,*,*]
    3: x=xin[*,*,*,indvalid,*,*,*,*]
    4: x=xin[*,*,*,*,indvalid,*,*,*]
    5: x=xin[*,*,*,*,*,indvalid,*,*]
    6: x=xin[*,*,*,*,*,*,indvalid,*]
    7: x=xin[*,*,*,*,*,*,*,indvalid]
    endcase
endif else x=xin
if n_elements(ySTDEVin) eq nyorg then ySTDEVin=ySTDEVin[indvalid]

; ----Check input----
nA = n_elements(A)
ny=n_elements(y)

if n_elements(Ainfo) ne nA then Ainfo=MakeAInfo(nA)
if not keyword_set(ySTDEVin) then begin
    ; Assume Poisson statistics
    if not keyword_set(weight) then weight=0
    ySTDEV=NLLSPoissonW(y,weight)
endif else ySTDEV=float(ySTDEVin) ; Because: -ySTDEVin gives problems for unsigned types

; ----Routine control flags----
bnumder=keyword_set(numder)
if n_elements(tol) EQ 0 then tol = 1.D-10 else tol=double(tol)      ;Convergence tolerance
gtol=tol
ftol=tol
xtol=tol
if n_elements(itmax) EQ 0 then itmax = 200 else itmax=fix(itmax)>0     ;Maximum # iterations

; ----Check tied parameters----
indtied=where(Ainfo.tied ne -1,ntied)
btied=ntied ne 0
if btied then begin
    Ainfo[indtied].fixed=1
    tied=Ainfo[indtied].tied
endif

; ----Check non-least squares parameters----
indnonls=where(Ainfo.nonls,nnonls)
bnonls=nnonls ne 0
if bnonls then begin
    Ainfo[indnonls].fixed=1
endif

; ----Check fixed parameters and initialize output----
indfree=where(~Ainfo.fixed,nfree)
dof = ny - nfree
sigma=A*0
blimits=0
if (dof le 0) then begin
    ErrorMsg='Not enough data points.'
    return,-1
endif
if nfree eq 0 then indfree=lindgen(nA)

; ----Initialize----
Areturn = A      ;; the set of parameters to be returned
Afree = Areturn[indfree]  ;; the set of free parameters
if bnonls then Anonls=Areturn[indnonls]
if btied then NLLS_tie,Areturn,indtied,tied
call_procedure, fname, x, Areturn, yfit,_EXTRA=extra
wResid=(y-yfit)/ySTDEV

; ----Deal with precision----
sz = size(wResid[0])
isdouble = (sz[sz[0]+1] EQ 5)
MACHEP0=floatpres(double=isdouble,DWARF=DWARF,GIANT=GIANT,rDWARF=rDWARF,rGIANT=rGIANT)
machvals={MACHEP0:MACHEP0,DWARF:DWARF,GIANT:GIANT,rDWARF:rDWARF,rGIANT:rGIANT}

szx = size(Afree)
tp = szx[szx[0]+1]
if isdouble AND tp NE 5 then begin
    Areturn=double(Areturn)
    Afree=double(Afree)
    if bnonls then Anonls=double(Anonls)
endif
if bnonls then Anonls_end=Anonls

; ----Initialize(continue...)----
wResidnorm = NLLS_norm(wResid,machvals)

; ----Define parameters for premature ending----
blimitscombo=0b

; ----No refined asked by user----
if nfree eq 0 or itmax eq 0 then begin
    CError=0
    goto,NLLS_TERMINATE
endif

; ----Check hard limits----
qllim = (Ainfo.limited)[0, indfree] eq 1
qulim = (Ainfo.limited)[1, indfree] eq 1
temp1 = where(qulim or qllim,ct)
blimits = ct ne 0
if blimits then begin
    llim  = (Ainfo.limits) [0, indfree]
    ulim  = (Ainfo.limits) [1, indfree]

    temp1 = where(qllim AND qulim AND llim ge ulim, ct)
    if ct GT 0 then begin
        ErrorMsg='Parameter limits are not consistent'
        return,-1
    endif
    temp1 = where((qllim AND Afree LT llim) OR $
        (qulim AND Afree GT ulim), ct)
    if ct GT 0 then begin
        temp1=strtrim(string(Afree[temp1],format='('+string(ct)+'(F,:,","))'),2)
        ErrorMsg='Parameters are not within limits: '+temp1
        return,-1
    endif
endif

; ----Check hard limits for linear combinations----
hlimits=histogram([Ainfo[indfree].comboID],binsize=1,rev=revlim)
nhlimits=n_elements(hlimits)-1
blimitscombo=nhlimits ne 0

if blimitscombo then begin
    qllimcombo=(Ainfo.limitedcombo)[0,indfree] eq 1 and Ainfo[indfree].comboID ne 0
    qulimcombo=(Ainfo.limitedcombo)[1,indfree] eq 1 and Ainfo[indfree].comboID ne 0
    combocoeff=Ainfo[indfree].combocoeff
    llimcombo=(Ainfo.limitscombo)[0,indfree]
    ulimcombo=(Ainfo.limitscombo)[1,indfree]

    ; Upper limits must be greater than lower limits
    temp1=where(qllimcombo and qulimcombo and llimcombo ge ulimcombo,ct)
    if ct GT 0 then begin
        ErrorMsg='Parameter limits are not consistent'
        return,-1
    endif

    for i=1,nhlimits do begin
        ni=revlim[i+1]-revlim[i]
        ; Missing linear constraint ID
        if ni eq 0 then begin
            ErrorMsg='Parameter limits are not consistent'
            return,-1
        endif

        temp1=revlim[revlim[i]:revlim[i+1]-1]

        ; Limits the same for this combo?
        if (total(llimcombo[temp1] eq llimcombo[temp1[0]],/pres) ne ni) or $
        (total(ulimcombo[temp1] eq ulimcombo[temp1[0]],/pres) ne ni) then begin
            ErrorMsg='Parameter limits are not consistent'
            return,-1
        endif

        temp2=(total(qllimcombo[temp1],/pres) eq ni)
        temp3=(total(qulimcombo[temp1],/pres) eq ni)
        ; All low and/or high constraints with the same limit?
        if temp2 eq 0 and temp3 eq 0 then begin
            ErrorMsg='Parameter limits are not consistent'
            return,-1
        endif

        ; Within constraints?
        temp4=total(Afree[temp1]*combocoeff[temp1],/pres)
        if (temp2 and temp4 lt llimcombo[temp1[0]]) or (temp3 and temp4 gt ulimcombo[temp1[0]]) then begin
            temp1=strtrim(string(Afree[temp1],format='('+string(ni)+'(F,:,","))'),2)
            ErrorMsg='Parameters are not within limits: '+temp1
            return,-1
        endif
    endfor

    ; Keep only needed info
    revlimall=revlim[nhlimits+2:*] ; all indices, sorted according to combo ID
    temp1=revlim[revlim[0:nhlimits]] ; first index for each bin (include combo ID=0)
    qllimcombo=qllimcombo[temp1]
    qulimcombo=qulimcombo[temp1]
    llimcombo=llimcombo[temp1]
    ulimcombo=ulimcombo[temp1]
endif

; ----Initialize Levelberg-Marquardt algo----
zero = Afree[0] * 0.
one  = zero + 1.
par = zero
QT_wResid = Afree * 0.
CError=0

; ----Outher loop----
iter = 1L
repeat begin
    ; Accept free parameters from previous loop
    Areturn[indfree]=Afree
    if bnonls then Areturn[indnonls]=Anonls ; Anonls contains the nonLS changed parameters
    if btied then NLLS_tie,Areturn,indtied,tied

    ; Jacobian of the weighted residuals: d(wResid)/dA = -1/ySTDEV*d(yfit)/dA
    if bnumder then begin
        JacResidT=make_array(ny,nfree,value=wResid[0]*0.)
        JacResidT=reform(JacResidT,ny,nfree,/overwrite)

        ; Small parameter increment
        eps=sqrt(MACHEP0)
        dA=eps*abs(Afree)
        temp1 = where(dA EQ 0, ct)
          if ct GT 0 then dA[temp1] = eps
          if blimits then begin
              temp1=where(qulim AND (Afree GT ulim-dA),ct)
              if ct GT 0 then dA[temp1] = -dA[temp1]
        endif

        ; One sided derivatives
        for i=0l,nfree-1 do begin
            P=Areturn
            P[indfree[i]]+=dA[i]
            if btied then NLLS_tie,P,indtied,tied
            call_procedure, fname, x, P, yfit1,_EXTRA=extra
            JacResidT[0,i] = (yfit-yfit1)/ySTDEV/dA[i]
        endfor
        ; Two sided derivatives
;        for i=0l,nfree-1 do begin
;            P=Areturn
;            P[indfree[i]]+=dA[i]
;            if btied then NLLS_tie,P,indtied,tied
;            call_procedure, fname, x, P, yfit1,_EXTRA=extra
;            P=Areturn
;            P[indfree[i]]-=dA[i]
;            if btied then NLLS_tie,P,indtied,tied
;            call_procedure, fname, x, P, yfit2,_EXTRA=extra
;            JacResidT[0,i] = (yfit2-yfit1)/ySTDEV/(2*dA[i])
;        endfor

        ; TEST: Compare Analytical and numerical derivative
        ; => Take dyfit/dA not dwResid/dA
        if keyword_set(debug) then begin
            P_=Areturn
            JacResidTnum=JacResidT*rebin(-ySTDEV,ny,nfree,/sample)
            JacResidT_=~Ainfo.fixed
            call_procedure, fname, x, P_, yfit_, JacResidT_,_EXTRA=extra
            if n_elements(JacResidT_) eq nA*ny and nfree ne nA then JacResidT_=JacResidT_[*,indfree]
            JacResidT_=reform(JacResidT_,ny,nfree,/overwrite)
            window
            yrange=[min(JacResidT_)<min(JacResidTnum),max(JacResidT_)>max(JacResidTnum)]
            for i=0l,nfree-1 do begin
                loadctrev,-39
                ;if i le 6 then xrange=[12.5,15] else xrange=[17,18]
                ;xrange=[40,60]
                plot,x,JacResidTnum[*,i],psym=4,/xs,xrange=xrange,title='Param'+string(i,format='(I0)')+': '+string(Afree[i]),/ys;,yrange=yrange
                oplot,x,JacResidT_[*,i],color=100
                wait,0.8
            endfor
            ;wdelete,!d.WINDOW
            goto,NLLS_TERMINATE
        endif
        
    endif else begin
        JacResidT=~Ainfo.fixed
        call_procedure, fname, x, Areturn, yfit, JacResidT,_EXTRA=extra
        ; User might calculate all partial derivatives: keep only those of the refined
        if n_elements(JacResidT) eq nA*ny and nfree ne nA then JacResidT=JacResidT[*,indfree]
        JacResidT=reform(JacResidT,ny,nfree,/overwrite)
        JacResidT/=rebin(-ySTDEV,ny,nfree,/sample)
    endelse

    ; Zero derivatives of parameters pegged at there limits
    if blimits then begin
        ; The direction of the steepest descent is -Grad(Chi-square)
        ; Grad(Chi-square) = 0.5*JT.wResid
        whlpeg = where(qllim AND (Afree EQ llim), nlpeg)
        whupeg = where(qulim AND (Afree EQ ulim), nupeg)

        ; See if any "pegged" values should keep their derivatives
        if (nlpeg GT 0) then begin
          ; Total derivative of sum wrt lower pegged parameters
          for i = 0L, nlpeg-1 do begin
              sum = total(wResid * JacResidT[*,whlpeg[i]])
              ; Chi-square gradient at parameter "whlpeg[i]" greater than zero
              ; so the parameter will become smaller and it is already at its 
              ; smallest limit => make 0
              if sum GT 0 then JacResidT[*,whlpeg[i]] = 0
          endfor
        endif
        if (nupeg GT 0) then begin
          ; Total derivative of sum wrt upper pegged parameters
          for i = 0L, nupeg-1 do begin
              sum = total(wResid * JacResidT[*,whupeg[i]])
              ; see above but now for highest limit
              if sum LT 0 then JacResidT[*,whupeg[i]] = 0
          endfor
        endif
    endif

    ; Do the same for linear combination constraints
    if blimitscombo then begin

        ; Get pegged values and zero derivatives if necessary (when going in the wrong direction)
        temp1=chunktotal((Afree*combocoeff)[revlimall],hlimits,/pres)

        wh = where(qllimcombo AND (temp1 eq llimcombo), nlpegcombo)
        if nlpegcombo ne 0 then begin
            whlpegcombo=-1L
            for k=0l,nlpegcombo-1 do begin
                i=wh[k]
                ni=revlim[i+1]-revlim[i]
                temp2=revlim[revlim[i]:revlim[i+1]-1]
                whlpegcombo=[whlpegcombo,temp2]
                for j=0l,ni-1 do begin
                    sum = total(wResid * JacResidT[*,temp2[j]])
                    if sum GT 0 then JacResidT[*,temp2[j]] = 0
                endfor
            endfor
            nlpegcombo=total(hlimits[wh],/pres)
            whlpegcombo=whlpegcombo[1:*]
        endif

        wh = where(qulimcombo AND (temp1 eq ulimcombo), nupegcombo)
        if nupegcombo ne 0 then begin
            whupegcombo=-1L
            for k=0l,nupegcombo-1 do begin
                i=wh[k]
                ni=revlim[i+1]-revlim[i]
                temp2=revlim[revlim[i]:revlim[i+1]-1]
                whupegcombo=[whupegcombo,temp2]
                for j=0l,ni-1 do begin
                    sum = total(wResid * JacResidT[*,temp2[j]])
                    if sum LT 0 then JacResidT[*,temp2[j]] = 0
                endfor
            endfor
            nupegcombo=total(hlimits[wh],/pres)
            whupegcombo=whupegcombo[1:*]
        endif
    endif

;    covarkeep=invert(JacResidT##transpose(JacResidT))

    ; Equation to solve: Gauss-Newton step x
    ; JacResid^T.JacResid.x = -JacResid^T.wResid
    ; JacResid = Q.R
    ; => R.x = -Q^T.wResid

    ; Compute the QR factorization of the jacobian
    ; JacResid (ny x nA) = Q(ny x ny) ## R(ny x na)
    ; JacResidT: nA x ny
    QRfact, JacResidT, ipvt, Rdiag, JacResidT_norm
    ; JacResidT: nA x ny => Q and the strict upper triangle of R (the upper triangle without the diagonal terms)
    ; Rdiag: nA => diagonal terms of R
    ; JacResidT_norm: nA => Jacobian column norms  (row norms of JacResidT before QRfact)
    ; ipvt: nA => decreasing order of JacResidT_norm

    ; Initialize step bound delta
    if (iter EQ 1) then begin
        diag = JacResidT_norm
        temp1 = where(diag EQ 0, ct)
        if ct GT 0 then diag[temp1] = one
        wa3 = diag * Afree
        Afree_norm=NLLS_norm(wa3,machvals)
        delta = 100.*Afree_norm
        if delta EQ zero then delta = zero + 100.
    endif

    ; Compute QT_wResid = the first n terms of Q^T ## wResid (which has m terms)
    ; and put Rdiag in JacResidT so JacResidT contains R
    temp1 = wResid
    for j=0L, nfree-1 do begin
        lj = ipvt[j]
        temp2 = JacResidT[j,lj]
        if temp2 NE 0 then begin
            fj = JacResidT[j:*,lj]
            wj = temp1[j:*]
            temp1[j] = wj - fj * total(fj*wj) / temp2
        endif
        JacResidT[j,lj] = Rdiag[j]
        QT_wResid[j] = temp1[j]
    endfor

    ; From this point on, only the square matrix, consisting of the
    ; upper triangle of R, is needed.
    JacResidT = JacResidT[0:nfree-1, 0:nfree-1]
    JacResidT = reform(JacResidT, nfree, nfree, /overwrite)
    JacResidT = JacResidT[*, ipvt]
    JacResidT = reform(JacResidT, nfree, nfree, /overwrite)

    ; Check for overflow
    temp1 = where(finite(JacResidT) EQ 0, ct)
    if ct gt 0 then begin
        ErrorMsg='Parameter or function values have become infinite.'
        CError=-1
        goto,NLLS_TERMINATE
    endif

    ; Compute the norm of the scaled gradient and check convergence
    gnorm = zero
    if wResidnorm NE 0 then begin
        for j=0L, nfree-1 do begin
          l = ipvt[j]
          if JacResidT_norm[l] NE 0 then begin
              sum = total(JacResidT[0:j,j]*QT_wResid[0:j])/wResidnorm
              gnorm >= abs(sum/JacResidT_norm[l])
          endif
        endfor
    endif
    if gnorm LE gtol then begin
        CError=4
        goto,NLLS_TERMINATE
    endif

    ; Rescale
    diag >= JacResidT_norm

    ; ----Inner loop----
    repeat begin
        ; Determine the levenberg-marquardt step:
        ; Remember that the outer loop calculates the QR-factorisation of JacResid
        ;     JacResid^T.JacResid.x = -JacResid^T.wResid
        ;     JacResid = R.Q
        ;     => R.x = -Q^T.wResid
        ;     
        ;     (Actually, we solve R.x = Q^T.wResid and x=-x afterwards)
        ;
        ; Now in the inner loop: additional LM-step condition
        ;     sqrt(par).D.x = 0
        ;
        ; So we need to solve A.x=B and C.x = 0
        ;
        ; NLLS_LMpar input
        ;     JacResidT: upper triangle of R
        ;     ipvt: pivoting array
        ;     diag: diagonal terms of D
        ;     QT_wResid: Q^T.wResid (first n terms)
        ;     delta: upperbound on dxnorm
        ;     par: initial estimate
        ; 
        ; NLLS_LMpar output
        ;     P^T(R^T.R+par.D.D).P=S^T.S
        ;     
        ;     JacResidT: lower triangle contains strict upper triangle of S
        ;     Rdiag: least squares solution x of the outer+inner system of equations
        ;     JacResidT_norm: diagonal terms of S (not used)
        ;     par: final estimate
        par = NLLS_LMpar(JacResidT, ipvt, diag, QT_wResid, delta, Rdiag, JacResidT_norm, machvals, par=par)
        
        ; Store the direction p and Afree+p. Calculate the norm of p
        Rdiag = -Rdiag

        if blimits or blimitscombo then begin
            ; Respect the limits.  If a step were to go out of bounds, then
            ; we should take a step in the same direction but shorter distance.
            ; The step should take us right to the limit in that case.
            alpha = one

            ; Do not allow any steps out of bounds
            if blimits then begin
                ; xnew = x + alpha * abs(Rdiagx) <= ulim and llim

                if nlpeg GT 0 then Rdiag[whlpeg] >= 0
                if nupeg GT 0 then Rdiag[whupeg] <= 0

                dRdiag = abs(Rdiag) GT MACHEP0

                wh = where(dRdiag AND qllim AND (Afree + Rdiag LT llim), ct)
                if ct GT 0 then alpha = min([alpha,(llim[wh]-Afree[wh])/Rdiag[wh]])

                wh = where(dRdiag AND qulim AND (Afree + Rdiag GT ulim), ct)
                if ct GT 0 then alpha = min([alpha,(ulim[wh]-Afree[wh])/Rdiag[wh]])
            endif
            if blimitscombo then begin
                ; a(x+alpha.|Rdiagx|)+b(y+alpha.|Rdiagy|)+... <= llim and ulim

                if nlpegcombo GT 0 then Rdiag[whlpegcombo] >= 0
                if nupegcombo GT 0 then Rdiag[whupegcombo] <= 0

                temp1=chunktotal(Afree[revlimall]*combocoeff[revlimall],hlimits,/pres)
                temp2=chunktotal(Rdiag[revlimall]*combocoeff[revlimall],hlimits,/pres)
                temp3=temp1+temp2

                dRdiag = abs(temp2) GT MACHEP0

                wh = where(dRdiag AND qllimcombo AND (temp3 LT llimcombo), ct)
                if ct GT 0 then alpha = min([alpha,(llimcombo[wh]-temp1[wh])/temp2[wh]])

                wh = where(dRdiag AND qulimcombo AND (temp3 GT ulimcombo), ct)
                if ct GT 0 then alpha = min([alpha,(ulimcombo[wh]-temp1[wh])/temp2[wh]])
            endif

            ; When alpha is too small, this might cause X-convergence and/or F-convergence while it's not converging
            ; Therefor, the fraction alpha of the lmstep taken must never be smaller than a certain number (take 1%)
            alpha>=0.01
            ; Scale the resulting vector
            Rdiag *= alpha
            JacResidT_norm = Afree + Rdiag

            ; Adjust the final output values.  If the step put us exactly
            ; on a boundary, make sure it is exact.
            if blimits then begin
                sgnu = (ulim GE 0)*2d - 1d
                sgnl = (llim GE 0)*2d - 1d

                wh = where(qulim AND JacResidT_norm GE (ulim*(1-sgnu*MACHEP0) - (ulim EQ 0)*MACHEP0), ct)
                if ct GT 0 then JacResidT_norm[wh] = ulim[wh]

                wh = where(qllim AND JacResidT_norm LE (llim*(1+sgnl*MACHEP0) + (llim EQ 0)*MACHEP0), ct)
                if ct GT 0 then JacResidT_norm[wh] = llim[wh]
            endif
            
            if blimitscombo then begin
                ; Change the parameter which needs the least change to set the linear ombination on its boundary
                sgnu = (ulimcombo GE 0)*2d - 1d
                sgnl = (llimcombo GE 0)*2d - 1d

                temp5=chunktotal(JacResidT_norm[revlimall]*combocoeff[revlimall],hlimits,/pres)

                wh = where(qllimcombo AND temp5 LE (llimcombo*(1+sgnl*MACHEP0) + (llimcombo EQ 0)*MACHEP0) and $
                                          temp5 ne llimcombo, ct)
                for i=0l,ct-1 do begin
                    j=wh[i]
                    nj=revlim[j+1]-revlim[j]
                    temp1=revlim[revlim[j]:revlim[j+1]-1]

                    if nj gt 1 then begin
                        temp2=rebin(JacResidT_norm[temp1]*combocoeff[temp1],nj,nj,/sample)
                        temp3=indgen(nj)
                        temp2[temp3,temp3]=0
                        temp2=(rebin([llimcombo[j]],nj,/sample)-total(temp2,1))/combocoeff[temp1]
                        temp4=min(abs(JacResidT_norm[temp1]-temp2),temp3)
                        JacResidT_norm[temp1[temp3]] = temp2[temp3]
                    endif else JacResidT_norm[temp1] = llimcombo[j]
                endfor

                wh = where(qulimcombo AND temp5 GE (ulimcombo*(1-sgnu*MACHEP0) - (ulimcombo EQ 0)*MACHEP0) and $
                                          temp5 ne ulimcombo, ct)
                for i=0l,ct-1 do begin
                    j=wh[i]
                    nj=revlim[j+1]-revlim[j]
                    temp1=revlim[revlim[j]:revlim[j+1]-1]

                    if nj gt 1 then begin
                        temp2=rebin(JacResidT_norm[temp1]*combocoeff[temp1],nj,nj,/sample)
                        temp3=indgen(nj)
                        temp2[temp3,temp3]=0
                        temp2=(rebin([ulimcombo[j]],nj,/sample)-total(temp2,1))/combocoeff[temp1]
                        temp4=min(abs(JacResidT_norm[temp1]-temp2),temp3)
                        JacResidT_norm[temp1[temp3]] = temp2[temp3]
                    endif else JacResidT_norm[temp1] = ulimcombo[j]
                endfor
            endif

        endif else begin
            ; No parameter limits, so just move to new position JacResidT_norm
            alpha = one
            JacResidT_norm = Afree + Rdiag
        endelse

        wa3 = diag * Rdiag
        pnorm = NLLS_norm(wa3,machvals)

        ; On the first iteration, adjust the initial step bound
        if iter EQ 1 then delta <= pnorm

        ; Evaluate the function at Afree+p and calculate its norm
        Areturn[indfree] = JacResidT_norm
        if bnonls then Areturn[indnonls]=Anonls ; Anonls contains the nonLS changed parameters
        if btied then NLLS_tie,Areturn,indtied,tied
        call_procedure, fname, x, Areturn, yfit1,_EXTRA=extra
        wResid1=(y-yfit1)/ySTDEV
        wResidnorm1=NLLS_norm(wResid1,machvals)

        ; Compute the scaled actual reduction
        actred = -one
        if 0.1D * wResidnorm1 LT wResidnorm then actred = - (wResidnorm1/wResidnorm)^2 + 1.

        ; Compute the scaled predicted reduction and the scaled directional
        ; derivative
        for j = 0L, nfree-1 do begin
            wa3[j] = 0
            wa3[0:j] += JacResidT[0:j,j]*Rdiag[ipvt[j]]
        endfor

        ; Remember, alpha is the fraction of the full LM step actually
        ; taken
        lmstep=alpha*par
        temp1 = NLLS_norm(alpha*wa3,machvals)/wResidnorm
        temp2 = (sqrt(lmstep)*pnorm)/wResidnorm
        half  = zero + 0.5
        prered = temp1*temp1 + (temp2*temp2)/half
        dirder = -(temp1*temp1 + temp2*temp2)

        ; Compute the ratio of the actual to the predicted reduction.
        ratio = zero
        tenth = zero + 0.1
        if prered NE 0 then ratio = actred/prered

        ; Update the step bound
        if ratio LE 0.25D then begin
            if actred GE 0 then temp1 = half $
            else temp1 = half*dirder/(dirder + half*actred)
            if ((0.1D*wResidnorm1) GE wResidnorm) OR (temp1 LT 0.1D) then temp1 = tenth
            delta = temp1*(delta < (pnorm/tenth))
            par /= temp1
        endif else begin
            if (par EQ 0) OR (ratio GE 0.75) then begin
                delta = pnorm/half
                par *= half
            endif
        endelse

        ; Test for successful iteration
        if ratio GE 0.0001 then begin
            ; Successful iteration.  Update Afree, wResid, and their norms
            Afree = JacResidT_norm
            JacResidT_norm *= diag

            yfit = yfit1
            wResid = wResid1
            Afree_norm = NLLS_norm(JacResidT_norm,machvals)
            wResidnorm = wResidnorm1
            
            ; Fixed parameters: accept for next iteration
            if bnonls then begin
                Anonls_end=Anonls ; Accept the nonLS changed parameters
                Anonls=Areturn[indnonls] ; next nonLS changed parameters
            endif
            iter++
        endif

        ; Tests for convergence
        if (abs(actred) LE ftol) AND (prered LE ftol) $
            AND  (0.5D * ratio LE 1) then CError = 1
        if delta LE xtol*Afree_norm then CError = 2
        if (CError EQ 1) AND (CError EQ 2) then CError = 3
        if CError NE 0 then goto, NLLS_TERMINATE

        ; Tests for termination and stringent tolerances
        if iter GE itmax then CError = 5
        if (abs(actred) LE MACHEP0) AND (prered LE MACHEP0) $
            AND (0.5*ratio LE 1) then CError = 6
        if delta LE MACHEP0*Afree_norm then CError = 7
        if gnorm LE MACHEP0 then CError = 8
        if lmstep lt 0 or ~finite(lmstep) or ~finite(ratio) then CError = 9
        if CError NE 0 then goto, NLLS_TERMINATE

        ; End of inner loop. Repeat if iteration unsuccessful
    endrep until ratio ge 0.0001

    ; Check for over/underflow
    wh = where(finite(Rdiag) EQ 0 OR finite(JacResidT_norm) EQ 0 OR finite(Afree) EQ 0, ct)
    if ct GT 0 OR finite(ratio) EQ 0 then begin
        ErrorMsg='Parameter or function values have become infinite.'
        CError=-1
    endif

endrep until CError ne 0

NLLS_TERMINATE:

; Finish last step
Areturn[indfree]=Afree
if bnonls then Areturn[indnonls]=Anonls_end
if btied then NLLS_tie,Areturn,indtied,tied

; TEST: visualize CHISQ space
;if keyword_set(debug) then begin
;    window
;    for k=26,26 do begin;nfree-1
;        j=indfree[k]
;        P=Areturn
;        if P[j] ne 0 then Pj=P[j]*makerange(0.1,1.9,0.01) $
;        else Pj=makerange(-1.5,1.5,0.01)
;        n=n_elements(Pj)
;        chi2=fltarr(n)
;        for i=0l,n-1 do begin
;            P=Areturn
;            P[j]=Pj[i]
;            if btied then NLLS_tie,P,indtied,tied
;            call_procedure, fname, x, P, yfit2,_EXTRA=extra
;            wResid2=(y-yfit2)/ySTDEV
;            chi2[i]=total(wResid2^2)
;        endfor
;        P=Areturn
;
;        plot,Pj,chi2,xstyle=1,ystyle=1,title='CHISQ space',xtitle='param'+stringr(j),ytitle='chisq'
;        plots,[P[j],P[j]],!y.crange,linestyle=2
;        wait,0.5
;    endfor
;endif

; TEST: visualize CHISQ space
;if keyword_set(debug) then begin
;    common chisqblock,Pi,Pj,Pm,chi2
;    i=47
;    j=48
;    m=49
;    P=Areturn
;    Pi=P[i]*makerange(-0.10,1.1,0.02)
;    Pj=P[j]*makerange(0.95,1.01,0.001)
;    Pm=P[m]*makerange(0.96,1.01,0.001)
;    nPi=n_elements(Pi)
;    nPj=n_elements(Pj)
;    nPm=n_elements(Pm)
;    chi2=fltarr(nPi,nPj,nPm)
;    s=[nPi,nPj,nPm]
;    ntot=nPi*nPj*nPm
;    T=systime(1)
;    for tind=0L,ntot-1 do begin
;        P=Areturn
;        tmp=array_indices(s,tind,/dim)
;        P[i]=Pi[tmp[0]]
;        P[j]=Pj[tmp[1]]
;        P[m]=Pm[tmp[2]]
;        if btied then NLLS_tie,P,indtied,tied
;        call_procedure, fname, x, P, yfit2,_EXTRA=extra
;        wResid2=(y-yfit2)/ySTDEV
;        chi2[tind]=total(wResid2^2)
;        temp1=tind*100./ntot
;        temp1=(systime(1)-T)*(100-temp1)/temp1
;        ConvertTime,temp1,timeaf
;        print,' Time left:',temp1,timeaf
;    endfor
;endif

; Number of parameters pegged at a hard limit value
if blimits then $
      wh = where((qulim AND (Afree EQ ulim)) OR $
                 (qllim AND (Afree EQ llim)), npegged)
if blimitscombo then begin
    temp1=chunktotal((Afree*combocoeff)[revlimall],hlimits,/pres)
    wh = where((qllimcombo and (temp1 eq llimcombo)) or $
               (qulimcombo and (temp1 eq ulimcombo)), temp2)
    if temp2 ne 0 then npegged+=total(hlimits[wh],/pres)
endif

; Chisq and A STDDEV
if CError ge 0 then begin
;    if n_elements(wResid1) then chisq=total((wResid > wResid1)^2) $
;    else chisq=wResidnorm^2
    chisq=wResidnorm^2

    if CError gt 0 then begin
        ; Covariance matrix
        covar = replicate(zero, nA, nA)
        cv=COvar(JacResidT[0:nfree-1,0:nfree-1],ipvt[0:nfree-1])
        cv = reform(cv, nfree, nfree, /overwrite)
        for i = 0L, nfree-1 do covar[indfree, indfree[i]] = cv[*,i]
;        covar=covarkeep
;        print,cerror

        ; Compute errors in parameters
        i = lindgen(nA)
        sigma = replicate(abs(covar[0])*0., nA)
        wh = where(covar[i,i] GE 0, ct)
        if ct GT 0 then sigma[wh] = sqrt(covar[wh, wh])
        
        if keyword_set(correlation) then begin
            ; Correlation matrix
            corr = covar * 0
            for i = 0, nA-1 do for j = 0, nA-1 do $
                corr[i,j] = covar[i,j]/sqrt(covar[i,i]*covar[j,j])
            ind=where(~finite(corr),ct)
            if ct ne 0 then corr[ind]=0
            corr[indgen(nA),indgen(nA)]=1

            correlation={corr:corr,Areturn:Areturn}
        endif
    endif
endif

; Error codes
case CError of
-1: ; Real error, message given for specific case
0: ErrorMsg='Nothing to refine.'
1: ErrorMsg='F-convergence: Both actual and predicted relative reductions in the sum of squares are at most FTOL.'
2: ErrorMsg='X-convergence: Relative error between two consecutive iterates is at most XTOL.'
3: ErrorMsg='F/X-convergence: Both actual and predicted relative reductions in the sum of squares are at most FTOL. Relative error between two consecutive iterates is at most XTOL.'
4: ErrorMsg='G-convergence: The cosine of the angle between wResid and any column of the jacobian is at most GTOL in absolute value.'
5: ErrorMsg='Stopped: Maximal number of iterations reached.'
6: ErrorMsg='Stopped: FTOL is too small. No further reduction in the sum of squares is possible.'
7: ErrorMsg='Stopped: XTOL is too small. No further improvement in the approximate solution Afree is possible.'
8: ErrorMsg='Stopped: GTOL is too small. wResid is orthogonal to the columns of the jacobian to machine precision.'
9: ErrorMsg='Error: Levenberg-marquardt step not valid.'
endcase

; Return A , yfit and ySTDEV
if nyorg ne ny then begin
    tmp=temporary(yfit)
    yfit=make_array(nyorg,value=tmp[0])
    yfit[*]=!Values.F_NAN
    yfit[indvalid]=temporary(tmp)
    
    ySTDEVin=yfit
    ySTDEVin[indvalid]=ySTDEV
endif else ySTDEVin=ySTDEV

A=Areturn
return,yfit
end;function NLLS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
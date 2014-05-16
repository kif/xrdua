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

function ExtractbitCode32,code,bits
; cfr MakebitCode32

nbits=n_elements(bits)
ncode=n_elements(code)
nvals=nbits*ncode

bits2=0
if nbits gt 1 then bits2=[bits2,bits[0:nbits-2]]
shft=total(bits2,/CUMULATIVE,/PRESERVE_TYPE)
band=2UL^(bits)-1
shft=rebin(shft,nvals,/sample)
band=rebin(band,nvals,/sample)
code2=reform(code#replicate(1,nbits),nvals)
vals=ishft(code2,-shft) and band

return,vals
end;function ExtractbitCode32
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakebitCode32,bits,vals

nbits=n_elements(bits)
nvals=n_elements(vals)
ncode=nvals/nbits

bits2=0
if nbits gt 1 then bits2=[bits2,bits[0:nbits-2]]
shft=total(bits2,/CUMULATIVE,/PRESERVE_TYPE)
shft=rebin(shft,nvals,/sample)
code2=reform(ishft(vals,shft),ncode,nbits)

code=ulonarr(ncode)
for i=0l,nbits-1 do code or= code2[*,i]

return,code
end;function MakebitCode32
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ExtractbitCode,ptr,codei,ibit=ibit
; cfr SetbitCode

bits=8
parts=4
i=codei/parts
j=codei mod parts

if n_elements(ibit) ne 0 then begin
    code=ishft((*ptr).code[i],-bits*j-ibit) mod 2
endif else begin
    extraxtbits=['FF'XUL,'FF00'XUL,'FF0000'XUL,'FF000000'XUL]
    code=ishft((*ptr).code[i] and extraxtbits[j],-bits*j)
endelse

return,code
end;function ExtractbitCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SetbitCode,ptr,codei,codeval,ibit=ibit,changed=changed
; code is a 32bit (unsigned long) array => n elements
; codei is a 8bit part of codeval => value from 0 to 8*n-1
; ibit is a 1bit part of codei => value 0,1,2,3,5,6,7

; when ibit: codeval is a 1bit number, to be set in codei 8bit part => value 0 or 1
; when no ibit: codeval is a 8bit number, to be set in codei 8bit part => value 0...31

bits=8 ; bits in a part
parts=4 ; parts in an element: 32/bits
i=codei/parts
j=codei mod parts

codeval=ulong(codeval)
if n_elements(ibit) ne 0 then begin
    newcode=(*ptr).code[i]
    if codeval ne 0 then newcode or= ishft(codeval<1UL,bits*j+ibit)  $
    else newcode and= not (ishft(1UL,bits*j+ibit))
endif else begin
    makezero=['FFFFFF00'XUL,'FFFF00FF'XUL,'FF00FFFF'XUL,'00FFFFFF'XUL]
    newcode=((*ptr).code[i] and makezero[j])
    if codeval ne 0 then newcode or= ishft(codeval<31UL,bits*j)
endelse

if n_elements(changed) eq 0 then changed=0b
changed or= newcode ne (*ptr).code[i]
(*ptr).code[i]=newcode
end;pro SetbitCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GlobalR,R,keep,xu,bpder=bpder,noglobal=noglobal
; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: Rjm(Rglobalj,xu)

bpder=keyword_set(bpder)

if keyword_set(noglobal) then begin
    if bpder then keep=R*0+1
    return,R
endif

c=!pi/180
cu=xu*c

if bpder then begin
    keep=[    [-(R[1]+2*R[2]/cu)/cu/xu],$ ; dR/du_H
            [xu*0+1],$                    ; dR/dr1
            [1./cu],$                    ; dR/dr2
            [1./(cu*cu)]]                ; dR/dr3
endif

return,R[0]+R[1]/cu+R[2]/cu/cu
end;function GlobalR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GlobalRpder,keep,off,pder,mult,toggle
; jm: jth phase, mth reflection
; return: dF/dRglobalj=dF/dRjm . dRjm/dRglobalj= mult*dRjm/dRglobalj

for i=0l,n_elements(keep)-2 do $
        if toggle[i] then pder[*,off+toggle[i,1]]+=mult*keep[i+1] ; skip dR/du_H
        
end;pro GlobalRpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro fitGlobalR,ptr,A,x,R

; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: fit Rjm(Rglobalj,xu)

x=180./!pi/x ; 1/2-theta (in rad^-1)
y=reform(R[0,*]) ;R

na=n_elements(A)
if n_elements(y) lt na then return

A=reform(POLY_FIT(x,y,2,yfit=yfit))

window,winid()
plot,x,y,xtitle='1/two-theta (in rad^-1)',ytitle='decay param',psym=2
oplot,x,yfit

if na gt 3 then begin
    y=reform(R[1,*]) ;R
    A=[A,reform(POLY_FIT(x,y,2,yfit=yfit))]

    window,winid()
    plot,x,y,xtitle='1/two-theta (in rad^-1)',ytitle='decay param',psym=2
    oplot,x,yfit
endif

end;pro fitGlobalR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GlobalAsym,A,keep,xu,bpder=bpder,noglobal=noglobal
; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: Ajm(Aglobalj,xu)

bpder=keyword_set(bpder)

if keyword_set(noglobal) then begin
    if bpder then keep=A*0+1
    return,A
endif

si=sin(xu/180.*!pi)
tmp1=sqrt(2)-1./si
tmp2=sqrt(2)-1./si^2
if bpder then begin
    keep=[    [!pi/180*cos(xu/180.*!pi)/si^2*(A[1]+2*A[2]/si)],$; dA/du_H
            [make_array(n_elements(xu),value=1)],$; dA/dA0
            [tmp1],$; dA/dA1
            [tmp2]]; dA/dA2
endif

return,A[0]+A[1]*tmp1+A[2]*tmp2

end;function GlobalAsym
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GlobalAsympder,keep,off,pder,mult,toggle,i0=i0,i1=i1
; jm: jth phase, mth reflection
; return: dF/dAsymglobalj=dF/dAsymjm . dAsymjm/dAsymglobalj= mult*dAsymjm/dAsymglobalj

if n_elements(i0) eq 0 or n_elements(i1) eq 0 then begin
    i0=0
    i1=(size(pder,/dim))[0]-1
endif

for i=0l,n_elements(keep)-2 do $
        if toggle[i] then pder[i0:i1,off+toggle[i,1]]+=mult*keep[i+1] ;i+1 to skip dA/du_H

end;pro GlobalAsympder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncGlobalAsym, X, A, F, pder, struc=struc
; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: Asymjm(Aglobalj,xu)
;          dAsymjm/dAsymglobalj

na=n_elements(A)
nx=n_elements(X)

sx=sin(X)
tmp1=sqrt(2)-1/sx
tmp2=sqrt(2)-1/sx^2.
F=A[0]+A[1]*tmp1+A[2]*tmp2

if N_PARAMS() ge 4 then begin
    pder=make_array(nx,na,type=size(A,/type))
    pder[*,0]=1
    pder[*,1]=tmp1
    pder[*,2]=tmp2
endif
end;pro gfuncGlobalAsym
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro fitGlobalAsym,ptr,A,x,As

; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: fit Ajm(Aglobalj,xu)

na=n_elements(A)
x*=!pi/180 ; 2-theta (in rad)
y=reform(As[0,*]) ; A

if n_elements(y) le na then return

yfit=NLLS(x,y, A,'gfuncGlobalAsym',errormsg=errormsg)

x*=180/!pi
ma=max(y,min=mi)
maf=max(yfit,min=mif)
ma=max([mi,ma,maf,mif],min=mi)
window,winid()
plot,x,y,xtitle=ChiXtitle(1),ytitle='A',psym=2,yrange=[mi,ma]+[-1,1]*(ma-mi)*0.05
oplot,x,yfit

end;pro fitGlobalAsym
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GlobalnmixSplit,nmix,keep,xu,bpder=bpder,noglobal=noglobal
; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: nmixjm(nmixglobalj,xu)

bpder=keyword_set(bpder)

if keyword_set(noglobal) then begin
    if bpder then begin
        keep=xu*0+1                ; dnL/du_H
        keep=[[keep],[keep]]    ; dnR/du_H
    endif
    
    return,nmix
endif

tmpL=nmix[1]/180.*!pi
tmpR=nmix[3]/180.*!pi
if bpder then begin
    keep=[    [replicate(tmpL,n_elements(xu))],$; dnL/du_H
            [replicate(tmpR,n_elements(xu))],$; dnR/du_H
            [xu*0+1],$                        ; dnL/dn0 = dnR/dn0
            [xu/180.*!pi]]                     ; dnL/dX = dnR/dX
endif

return,[nmix[0]+tmpL*xu,nmix[2]+tmpR*xu]
end;function GlobalnmixSplit
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GlobalnmixSplitpder,keep,off,pder,multL,multR,toggle,i0=i0,i1=i1
; jm: jth phase, mth reflection
; return: dF/dnmixglobalj=dF/dnmixjm . dnmixjm/dnmixglobalj= mult*dnmixjm/dnmixglobalj

if n_elements(i0) eq 0 or n_elements(i1) eq 0 then begin
    i0=0
    i1=(size(pder,/dim))[0]-1
endif

for i=0l,1 do $ 
    if toggle[i] then $
        ;dnL/dn0L, dnL/dXL 
        pder[i0:i1,off+toggle[i,1]]+=multL*keep[i+2] ;i+2 to skip dnL/du_H and dnR/du_H

for i=2,3 do $ 
    if toggle[i] then $
        ;dnR/dn0R, dnR/dXR 
        pder[i0:i1,off+toggle[i,1]]+=multR*keep[i] ;dnL/dn0 = dnR/dn0, dnL/dX = dnR/dX
    
end;pro GlobalnmixSplitpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Globalnmix,nmix,keep,xu,bpder=bpder,noglobal=noglobal
; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: nmixjm(nmixglobalj,xu)

bpder=keyword_set(bpder)

if keyword_set(noglobal) then begin
    if bpder then keep=nmix*0+1
    return,nmix
endif

tmp=nmix[1]/180.*!pi
if bpder then begin
    keep=[    [replicate(tmp,n_elements(xu))],$; dn/du_H
            [xu*0+1],$                        ; dn/dn0
            [xu/180.*!pi]]                     ; dn/dX
endif

return,nmix[0]+tmp*xu
end;function Globalnmix
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Globalnmixpder,keep,off,pder,mult,toggle,i0=i0,i1=i1
; jm: jth phase, mth reflection
; return: dF/dnmixglobalj=dF/dnmixjm . dnmixjm/dnmixglobalj= mult*dnmixjm/dnmixglobalj

if n_elements(i0) eq 0 or n_elements(i1) eq 0 then begin
    i0=0
    i1=(size(pder,/dim))[0]-1
endif

for i=0l,n_elements(keep)-2 do $
        if toggle[i] then pder[i0:i1,off+toggle[i,1]]+=mult*keep[i+1] ;i+1 to skip dn/du_H
        
end;pro Globalnmixpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro fitGlobalnmix,ptr,A,x,nmix

; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: fit nmixjm(nmixglobalj,xu)

na=n_elements(A)
x*=!pi/180. ; 2-theta (in rad)
y=reform(nmix[0,*]) ;nmix

A=linfit(x,y,yfit=yfit)

window,winid()
x*=180/!pi
plot,x,y,xtitle=ChiXtitle(1),ytitle='mixing param',psym=2
oplot,x,yfit

if na gt 2 then begin
    y=reform(nmix[1,*]) ;nmix
    A=[A,linfit(x,y,yfit=yfit)]

    window,winid()
    x*=180/!pi
    plot,x,y,xtitle=ChiXtitle(1),ytitle='mixing param',psym=2
    oplot,x,yfit
endif

end;pro fitGlobalnmix
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GlobalFWHMSplitPV,FWHM,typeinfo,keep,xu,bpder=bpder,noglobal=noglobal,MODSPV=MODSPV
; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: FWHMjm(FWHMglobalj,xu)

bpder=keyword_set(bpder)
if typeinfo.itable eq 0 then m=1 else m=typeinfo.ltable[2,typeinfo.itable]/typeinfo.ltable[2,0]

if noglobal then begin
    if bpder then begin
        keep=xu*0+m                ; dWL/du_H
        keep=[[keep],[keep]]    ; dWG/du_H
    endif
    
    if keyword_set(MODSPV) then return,FWHM*m $ ; WL and WG are given
    else return,[FWHM,FWHM]*m  ; WL=WG
endif

c=!pi/360.
cu=xu*c
sincu=sin(cu)
coscu=cos(cu)
tancu=tan(cu)

WL=sqrt(FWHM[0]*tancu^2+FWHM[1]*tancu+FWHM[2]+FWHM[3]/coscu^2)
if keyword_set(MODSPV) then begin
    WG=sqrt(FWHM[4]*tancu^2+FWHM[5]*tancu+FWHM[6]+FWHM[3]/coscu^2)
endif else begin
    WG=WL
endelse
FWHMret=[WL,WG]

if bpder then begin
    if keyword_set(MODSPV) then begin
        tmpL=(m/2.)/WL
        tmpG=(m/2.)/WG
        
        keep=[    [(2*c*(FWHM[0]+FWHM[3])*sincu/coscu^3+c*FWHM[1]/coscu^2)*tmpL],$ ; dWL/du_H
                [(2*c*(FWHM[4]+FWHM[3])*sincu/coscu^3+c*FWHM[5]/coscu^2)*tmpG],$ ; dWG/du_H
                [tmpL*tancu^2],$; dWL/dU
                [tmpL*tancu],$    ; dWL/dV
                [tmpL],$        ; dWL/dW
                [tmpL/coscu^2],$; dWL/dIG
                [tmpG*tancu^2],$; dWG/dU
                [tmpG*tancu],$    ; dWG/dV
                [tmpG],$        ; dWG/dW
                [tmpG/coscu^2]]    ; dWG/dIG
    endif else begin
        tmpL=(m/2.)/WL
        
        keep=[    [(2*c*(FWHM[0]+FWHM[3])*sincu/coscu^3+c*FWHM[1]/coscu^2)*tmpL],$ ; dWL/du_H = dWG/du_H
                [tmpL*tancu^2],$; dWL/dU = dWG/dU
                [tmpL*tancu],$    ; dWL/dV = dWG/dV
                [tmpL],$        ; dWL/dW = dWG/dW
                [tmpL/coscu^2]]; dWL/dIG = dWG/dIG
;        keep=keep[*,[0,0,1,2,3,0,1,2,3]]
    endelse
endif

return,FWHMret*m
end;function GlobalFWHMSplitPV
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GlobalFWHMmodSplitPVpder,keep,off,pder,multH1,multH2,toggle,i0=i0,i1=i1
; jm: jth phase, mth reflection
; return: dF/dFWHMglobalj=dF/dFWHMjm . dFWHMjm/dFWHMglobalj= mult*dFWHMjm/dFWHMglobalj

if n_elements(i0) eq 0 or n_elements(i1) eq 0 then begin
    i0=0
    i1=(size(pder,/dim))[0]-1
endif

for i=0l,2 do $ ; dWL/dU, dWL/dV, dWL/dW
        if toggle[i] then pder[i0:i1,off+toggle[i,1]]+=multH1*keep[i+2] ;i+2 to skip dWL/du_H and dWG/du_H
for i=4,6 do $ ; dWG/dU, dWG/dV, dWG/dW
        if toggle[i] then pder[i0:i1,off+toggle[i,1]]+=multH2*keep[i+2] ;i+2 to skip dWL/du_H and dWG/du_H

; dF/dIG = dF/dWL.dWL/dIG + dF/dWG.dWG/dIG
i=3
if toggle[i] then pder[i0:i1,off+toggle[i,1]]+=multH1*keep[i+2]+multH2*keep[i+6] ;i+2 to skip dWL/du_H and dWG/du_H

end;pro GlobalFWHMmodSplitPVpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GlobalFWHM,FWHM,typeinfo,keep,xu,bpder=bpder,noglobal=noglobal
; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: FWHMjm(FWHMglobalj,xu)

bpder=keyword_set(bpder)
if typeinfo.itable eq 0 then m=1 else m=typeinfo.ltable[2,typeinfo.itable]/typeinfo.ltable[2,0]

if noglobal then begin
    if bpder then keep=FWHM*0+m
    return,FWHM*m
endif

c=!pi/360.
cu=xu*c
sincu=sin(cu)
coscu=cos(cu)
tancu=tan(cu)

FWHMret=sqrt(FWHM[0]*tancu^2+FWHM[1]*tancu+FWHM[2]+FWHM[3]/coscu^2)

if bpder then begin
    tmp=(m/2.)/FWHMret
    
    keep=[    [(2*c*(FWHM[0]+FWHM[3])*sincu/coscu^3+c*FWHM[1]/coscu^2)*tmp],$ ; dW/du_H
            [tmp*tancu^2],$    ; dW/dU
            [tmp*tancu],$        ; dW/dV
            [tmp],$                ; dW/dW
            [tmp/coscu^2]]; dW/dIG
endif

return,FWHMret*m
end;function GlobalFWHM
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GlobalFWHMpder,keep,off,pder,mult,toggle,i0=i0,i1=i1
; jm: jth phase, mth reflection
; return: dF/dFWHMglobalj=dF/dFWHMjm . dFWHMjm/dFWHMglobalj= mult*dFWHMjm/dFWHMglobalj

if n_elements(i0) eq 0 or n_elements(i1) eq 0 then begin
    i0=0
    i1=(size(pder,/dim))[0]-1
endif

for i=0l,n_elements(keep)-2 do $
        if toggle[i] then pder[i0:i1,off+toggle[i,1]]+=mult*keep[i+1] ;i+1 to skip dW/du_H

end;pro GlobalFWHMpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GlobalFWHMntchz,FWHM,typeinfo,keep,xu,nmix,pderL=pderL,pderG=pderG,LG=LG,noglobal=noglobal,bpder=bpder
; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: FWHMjm(FWHMglobalj,xu)
;          nmixjm(FWHMglobalj,xu)

;FWHM=[FWHM_L,FWHM_G]

bpder=keyword_set(bpder)
if typeinfo.itable eq 0 then m=1 else m=typeinfo.ltable[2,typeinfo.itable]/typeinfo.ltable[2,0]

c1=2.69269
c2=2.42843
c3=4.47163
c4=0.07842
c5=1.36603
c6=-0.47719
c7=0.1116

if keyword_set(noglobal) then begin
    n=n_elements(FWHM)/2
    FWHML=FWHM[0:n-1]
    FWHMG=FWHM[n:*]
    FWHMret=FWHM*m
    
    if bpder then begin
        keep=xu*0+m                ; dWL/du_H
        keep=[[keep],[keep]]    ; dWG/du_H
    endif
endif else begin
    c=!pi/360.
    cu=xu*c
    sincu=sin(cu)
    coscu=cos(cu)
    tancu=tan(cu)

    FWHMG=sqrt(FWHM[0]*tancu*tancu+FWHM[1]*tancu+FWHM[2]+FWHM[3]/coscu^2)
    FWHML=FWHM[4]*tancu+FWHM[5]/coscu
    FWHMret=[FWHML,FWHMG]*m

    if bpder then begin
        keep=[[c*(FWHM[4]+FWHM[5]*sincu)/coscu^2],$ ; dWL/du_H
            [(2*c*(FWHM[0]+FWHM[3])*sincu/coscu^3+c*FWHM[1]/coscu^2)/(2*FWHMG)]] ; dWG/du_H
    endif
endelse

T=FWHMG^5.+c1*FWHMG^4.*FWHML+c2*FWHMG^3.*FWHML^2.+c3*FWHMG^2.*FWHML^3.+c4*FWHMG*FWHML^4.+FWHML^5.
if keyword_set(LG) then return,T^(0.2)
pderL=T^(-0.2)
T^=(-1.2)
T*=(-0.2)*FWHML

q=FWHML*pderL
nmix=c5*q+c6*q*q+c7*q*q*q

if ~bpder then return,FWHMret

pder=c5+2*c6*q+3*c7*q*q

pderL+=T*(c1*FWHMG^4.+2*c2*FWHMG^3.*FWHML+3*c3*FWHMG^2.*FWHML^2.+4*c4*FWHMG*FWHML^3.+5*FWHML^4.)
pderL*=pder
pderG=pder*T*(5*FWHMG^4.+4*c1*FWHMG^3.*FWHML+3*c2*FWHMG^2.*FWHML^2.+2*c3*FWHMG*FWHML^3.+c4*FWHML^4.)

if ~keyword_set(noglobal) then begin
    tmp=1./(2.*FWHMG)
    keep=[    [keep],$                ; dWL/du_H, dWG/du_H
            [tancu*tancu*tmp],$        ; dWG/dU
            [tancu*tmp],$            ; dWG/dV
            [tmp],$                    ; dWG/dW
            [tmp/coscu^2],$            ; dWG/dIG
            [tancu],$                ; dWL/dX
            [1./coscu]]                ; dWL/dY
    keep*=m
endif

return,FWHMret

end;function GlobalFWHMntchz
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GlobalFWHMtchzpder,keep,off,pder,multL,multG,toggle
; jm: jth phase, mth reflection
; return: dF/dFWHMglobalj

for i=0l,3 do $ ; dWG/dU, dWG/dV, dWG/dW, dWG/dIG
        if toggle[i] then pder[*,off+toggle[i,1]]+=multG*keep[i+2] ;i+2 to skip dWL/du_H and dWG/du_H
for i=4,5 do $ ; dWL/dX, dWL/dY
        if toggle[i] then pder[*,off+toggle[i,1]]+=multL*keep[i+2] ;i+2 to skip dWL/du_H and dWG/du_H

end;pro GlobalFWHMtchzpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncGlobalFWHM, X, A, F, pder, struc=struc

; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: FWHMjm(FWHMglobalj,xu)
;          dFWHMjm/dFWHMglobalj

na=n_elements(A)
nx=n_elements(X)

xu2=X/360.*!pi
tmp=tan(xu2)
if na gt 4 then begin
    nm=nx/2
    F=A[0]*tmp[0:nm-1]*tmp[0:nm-1]+A[1]*tmp[0:nm-1]+A[2]+A[3]/(cos(xu2[0:nm-1]))^2.
    if na eq 7 then F=[F,A[4]*tmp[nm:*]*tmp[nm:*]+A[5]*tmp[nm:*]+A[6]+A[3]/(cos(xu2[nm:*]))^2.] $
    else F=[F,A[4]*tmp[nm:*]+A[5]/cos(xu2[nm:*])]
endif else F=A[0]*tmp*tmp+A[1]*tmp+A[2]+A[3]/(cos(xu2))^2.

if N_PARAMS() ge 4 then begin
    pder=make_array(nx,na,type=size(A,/type))
    if na gt 4 then begin
        pder[0:nm-1,2]=1 ;W
        pder[0:nm-1,3]=1./(cos(xu2[0:nm-1]))^2.; IG
        pder[0:nm-1,1]=tmp[0:nm-1] ;V
        pder[0:nm-1,0]=tmp[0:nm-1]*tmp[0:nm-1] ;U
        if nA eq 7 then begin
            pder[nm:*,6]=1 ;W'
            pder[nm:*,5]=tmp[nm:*] ;V'
            pder[nm:*,4]=tmp[nm:*]*tmp[nm:*] ;U'
        endif else begin
            pder[nm:*,4]=tmp[nm:*] ;X
            pder[nm:*,5]=1./cos(xu2[nm:*]) ;Y
        endelse
    endif else begin
        pder[*,2]=1 ;W
        pder[*,3]=1./(cos(xu2))^2.; IG
        pder[*,1]=tmp ;V
        pder[*,0]=tmp*tmp ;U
    endelse
endif
end;pro gfuncGlobalFWHM
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro fitGlobalFWHM,ptr,A,x,FWHM

; xu: two-theta in degrees
; jm: jth phase, mth reflection
; return: fit FWHMjm(FWHMglobalj,xu)

na=n_elements(A)
;x: 2-theta (in deg)
y=reform(FWHM[0,*]) ; FWHM (in degrees)

if n_elements(y) le na then return

y*=y

if na gt 4 then begin ; i.e. 6 or 8
    x=[x,x]
    FWHM=reform(FWHM[1,*])
    if n_elements(A) eq 7 then y=[y,FWHM*FWHM] else y=[y,FWHM]
endif


yfit=NLLS(x,y, A,'gfuncGlobalFWHM',errormsg=errormsg)
A>=0
y=sqrt(y)
yfit=sqrt(yfit)

ma=max(y,min=mi)
maf=max(yfit,min=mif)
ma=max([mi,ma,maf,mif],min=mi)
window,winid()
plot,x,y,xtitle=ChiXtitle(1),ytitle='FWHM',psym=2,yrange=[mi,ma]+[-1,1]*(ma-mi)*0.05
oplot,x,yfit

end;pro fitGlobalFWHM
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GlobalPos,xu,typeinfo,keep,hkl,toggleu,bpder=bpder,noglobal=noglobal,xunoshift=tt

; xu: ...
; jm: jth phase, mth reflection
; return: ttjm(ttglobalj,xu)

; sin(theta_int) = sin(theta_init)/gain
; sin(theta_int) = lambda/(2*d)
; dist*tan(2theta) = (dist+delta_dist)*tan(2theta_int)

; zero=delta_dist/(dist*10^6)

bpder=keyword_set(bpder)

case typeinfo.type of
10: begin
    if keyword_set(noglobal) then begin
        if typeinfo.itable ne 0 then begin
            m=typeinfo.lambda/typeinfo.ltable[0,0]
            
            tmp=m*sin(!pi/360.*xu) ; sin(theta0)
            tt=360./!pi*asin(tmp) ; 2theta (in deg)
            if bpder then keep=m*cos(!pi/360.*xu)/sqrt(1-tmp*tmp)
            return,tt
        endif else begin
            if bpder then keep=xu*0+1
            tt=xu
            return,xu
        endelse
    endif
    
    if typeinfo.itable eq 0 then m=1 else m=typeinfo.lambda/typeinfo.ltable[0,0]
        
    ; xu: zero, gain, two-theta
    zero=xu[0]
    gain=xu[1]
    
    tmp=m*sin(!pi/360.*xu[2:*]) ; sin(theta0)
    tt=360./!pi*asin(tmp/gain) ; 2theta (in deg)
    
    if bpder then keep=tmp/sqrt(1.-(tmp/gain)^2.)/(-!pi/360.*gain^2.)
    ; keep rows (ncolumns == npeaks):
    ; du/dg
    nrow=1
    npeak=n_elements(xu[2:*])
        
    endcase
else:begin
    ; xu: zero,celparam, h,k,l
    zero=xu[0]
    npeak=(n_elements(xu)-7)/3
    hkl=transpose(reform(xu[7:*],npeak,3))
    if n_elements(toggleu) ne 0 then connect=toggleu[1:*,2]-1
    tt=reform(ttpdercell(hkl,xu[1:6],typeinfo.lambda,pder=keep,bpder=bpder,connect=connect))
    ; keep rows (ncolumns == nhkl):
    ; duhkl/da
    ; duhkl/db
    ; duhkl/dc
    ; duhkl/dalpha
    ; duhkl/dbeta
    ; duhkl/dgamma
    nrow=6

    endcase
endcase

case typeinfo.zerotype of
; TODO: Check McCuster guidelines for other zero-shifts
1:  begin ; Zero: delta sample-detector
    ttrad = tt/180.*!pi
    
    if zero eq 0 then begin
      a = 1
      ret = tt
      
      if bpder then begin
        keep=[[(90./!pi)*sin(2*ttrad)/(typeinfo.dist*1000.)],[keep]]
      endif
    endif else begin
      a = 1+zero/(typeinfo.dist*1000.)
      ret = 180./!pi*atan(a*tan(ttrad))
      
      ind = where(ret lt 0,ct)
      if ct ne 0 then ret[ind] += 180
      
      if bpder then begin
        tantt = tan(ttrad)
        b = (180./!pi)/((a*tantt)^2.+1)

        keep*=rebin(!pi/180.*a*b/(cos(ttrad))^2.,npeak,nrow,/sample)
        keep=[[b*tantt/(typeinfo.dist*1000.)],[keep]]
      endif
    endelse
    
    
    endcase
else:begin ; Zero: zero shift
    ret=tt+zero
    if bpder then keep=[[replicate(1,npeak)],[keep]]
    endelse
endcase

    ; keep rows (ncolumns == npeaks):
    ; du/dzero
    ; du/dg

    ; OR
    ; keep rows (ncolumns == nhkl):
    ; du/dzero
    ; duhkl/da
    ; duhkl/db
    ; duhkl/dc
    ; duhkl/dalpha
    ; duhkl/dbeta
    ; duhkl/dgamma

return,ret
end;function GlobalPos
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GlobalPospder,keep,off,pder,mult,toggle,i0=i0,i1=i1

; keep: dttj/dttglobalj
; mult: dF/dttj
; jm: jth phase, mth reflection
; return: dF/dttglobalj

if n_elements(i0) eq 0 or n_elements(i1) eq 0 then begin
    i0=0
    i1=(size(pder,/dim))[0]-1
endif

for i=0l,n_elements(keep)-1 do $
    if toggle[i] then pder[i0:i1,off+toggle[i,1]]+=mult*keep[i]

end;pro GlobalPospder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GlobalInt,I,typeinfo,keep,twotheta,hkl,bpder=bpder,noglobal=noglobal,bFcalc=bFcalc,Fcalc=Fcalc

; I: ...
; jm: jth phase, mth reflection
; return: Ijm(Iglobalj,xu)

if typeinfo.itable eq 0 then m=1 else m=typeinfo.ltable[1,typeinfo.itable]/typeinfo.ltable[1,0]

case typeinfo.type of
10:    begin
    ret=GetPDRelIntpder(I,twotheta,typeinfo,keep=keep,bpder=bpder,noglobal=noglobal)
    ret*=m
    if bpder then keep*=m
    if bFcalc then begin
        Fcalc=[[make_array(n_elements(twotheta),2,value=1)],[ret]];[1,1,I]
    endif

    ; keep:
    ; dIhkl/dS
    ; dIhkl/duhkl
    ; dIhkl/dIparam(=1)
    endcase
11:    begin
    ret=GetPawleyRelIntpder(I,twotheta,typeinfo,keep=keep,bpder=bpder,noglobal=noglobal,Fcalc=tmp)
    ret*=m
    if bpder then keep*=m
    if bFcalc then begin
        tmp[*,1]*=m
        Fcalc=[[temporary(tmp)],[ret]] ;[|F|,S.mH.LGP,I]
    endif
    
    ; keep:
    ; dIhkl/dS
    ; dIhkl/duhkl
    ; dIhkl/dIparam(=mhkl)
    endcase
12:    begin
    ret=GetRietveldRelIntpder(I,twotheta,hkl,typeinfo,keep=keep,bpder=bpder,Fcalc=tmp)
    ret*=m
    if bpder then keep*=m
    if bFcalc then begin
        tmp[*,1]*=m
        Fcalc=[[temporary(tmp)],[ret]] ;[|F|,S.mH.LGP,I]
    endif
    
    ; ret: I_hkl=S.m_hkl.|F_hkl|^2.LGP
    ; keep rows (ncolumns == nhkl):
    ; dIhkl/dS
    ; dIhkl/duhkl
    ; dIhkl/dxfit0
    ; dIhkl/dxfit1
    ; ...
    ; dIhkl/dyfit0
    ; dIhkl/dyfit1
    ; ...
    ; dIhkl/dzfit0
    ; dIhkl/dzfit1
    ; ...
    ; dIhkl/dSOF0
    ; dIhkl/dSOF1
    ; ...
    ; dIhkl/dBiso0
    ; dIhkl/dBiso1
    ; ...
    endcase
endcase

return,ret
end;function GlobalInt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GlobalIntpder,keep,off,pder,mult,toggle,typeinfo,i0=i0,i1=i1

; jm: jth phase, mth reflex
; return: dF/dIglobalj

if n_elements(i0) eq 0 or n_elements(i1) eq 0 then begin
    i0=0
    i1=(size(pder,/dim))[0]-1
endif

if toggle[0] then pder[i0:i1,off+toggle[0,1]]+=mult*keep[0] ; dI_H/dS

; keep[1]: dI_H/du_H
; keep[2]: dI_H/dI_Hparam

case typeinfo.type of
12: begin

    for i=1,n_elements(keep)-2 do $
        if toggle[i] then pder[i0:i1,off+toggle[i,1]]+=mult*keep[i+1] ;i+1 to skip dI_H/du_H

    endcase
else:
endcase

end;pro GlobalIntpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncOpvII, X, xu, I, FWHM, Rin, F, pder, off, bref, toggleu, toggleI, toggleW, toggleR, typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc

; Initiate
bpder=n_elements(pder) ne 0

; Convert global parameters to peakparameters
xu=GlobalPos(xu,typeinfo,keepxu,hkl,toggleu,bpder=bpder,noglobal=~bref[0,1],xunoshift=xunoshift)
I=GlobalInt(I,typeinfo,keepI,xunoshift,hkl,bpder=bpder,noglobal=~bref[1,1],bFcalc=bFcalc,Fcalc=Fcalc)
W=GlobalFWHM(FWHM,typeinfo,keepW,xunoshift,bpder=bpder,noglobal=~bref[2,1])
R=GlobalR(Rin,keepR,xunoshift,bpder=bpder,noglobal=~bref[3,1])

; Derived parameters
peaknr=n_elements(xu)
if ptr_valid(typeinfo.pPeakParam) then begin
    *typeinfo.pPeakParam=[reform(xu,1,peaknr),reform(I,1,peaknr),reform(W,1,peaknr),reform(R,1,peaknr)]
    return
endif

c=2.^(1./R)-1
gR=exp(lngamma(R)-lngamma(R-0.5)) ; == gamma(R)/gamma(R-0.5) (prevent Overflow)
dgR=digamma(R)-digamma(R-0.5)
A=2*gR*sqrt(c)/sqrt(!pi)

; Reliability factors
if bFcalc then begin
    Fcalc=[[temporary(Fcalc)],[make_array(peaknr,n_elements(X),type=size(Fcalc,/type))]]
    
    ; Fcalc
    ; cte=Iobs/Fobs^2
    ; Icalc
    ; Profiles
endif

for k=0l,peaknr-1 do begin
    ; Normalized Profile function (area=1)
    Barg=1+(4*c[k]/W[k]/W[k])*(X-xu[k])^(2.)
    B=Barg^(-R[k])
    multpderI=(A[k]/W[k])*B

    ; Profile function
    if keyword_set(SepF) then F=[[F],[I[k]*multpderI]] $
    else F+=I[k]*multpderI
    if bFcalc then Fcalc[k,3:*]=multpderI

    IF bpder THEN BEGIN
        ; dF/duH
        if bref[0] then multpderu=(A[k]*I[k]/W[k]^3*8*R[k]*c[k])*Barg^(-R[k]-1)*(X-xu[k])
        
        ; dF/dW
        if bref[2] or (bref[0] and bref[2,1]) then multpderW=(A[k]*I[k]/W[k]^4*8*R[k]*c[k])*Barg^(-R[k]-1)*(X-xu[k])^2. - (A[k]*I[k]/W[k]^2)*B
        
        ; dF/dR
        if bref[3] or (bref[0] and bref[3,1]) then begin
            if c[k] eq 0 then multpderR=0 else $
            multpderR=(I[k]*A[k]/W[k])*B*(dgR[k]-alog(Barg)+(4*2^(1./R[k])*alog(2)/(R[k]*W[k]^2))*(X-xu[k])^2./Barg)- $
                        (gR[k]*I[k]*2^(1./R[k])*alog(2)/(sqrt(!pi*c[k])*W[k]*R[k]^2))*B
        endif

        ; dF/duH += dF/dI.dI/duH
        if bref[0] and bref[1,1] then multpderu+=multpderI*keepI[k,1]
        ; dF/duH += dF/dW.dW/duH
        if bref[0] and bref[2,1] then multpderu+=multpderW*keepW[k,0]
        ; dF/duH += dF/dR.dR/duH
        if bref[0] and bref[3,1] then multpderu+=multpderR*keepR[k,0]

        ; pder to peak or global position parameters
        if bref[0] then begin
            if bref[0,1] then GlobalPospder,keepxu[k,*],off[0],pder,multpderu,toggleu $
            else pder[*,off[0]+k]+=multpderu*keepxu[k,0]; pder xu[k]
        endif
        
        ; pder to peak or global intensity parameters
        if bref[1] then begin
            if bref[1,1] then GlobalIntpder,keepI[k,*],off[1],pder,multpderI,toggleI,typeinfo $
            else pder[*,off[1]+k]+=multpderI*keepI[k,2] ; pder I[k]
        endif

        ; pder to peak or global FWHM parameters
        if bref[2] then begin
            if bref[2,1] then GlobalFWHMpder,keepW[k,*],off[2],pder,multpderW,toggleW $
            else pder[*,off[2]+k]+=multpderW*keepW[k,0] ; pder W[k]
        endif

        ; pder to peak or global decay parameters
        if bref[3] then begin
            if bref[3,1] then GlobalRpder,keepR[k,*],off[3],pder,multpderR,toggleR $
            else pder[*,off[3]+k]+=multpderR*keepR[k,0] ; pder R[k]
        endif

    endif
endfor

end;pro gfuncOpvII
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncOtchz, X, xu, I, FWHM, F, pder, off, bref, toggleu, toggleI, toggleW, typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc

; Initiate
bpder=n_elements(pder) ne 0

; Convert global parameters to peakparameters
xu=GlobalPos(xu,typeinfo,keepxu,hkl,toggleu,bpder=bpder,noglobal=~bref[0,1],xunoshift=xunoshift)
I=GlobalInt(I,typeinfo,keepI,xunoshift,hkl,bpder=bpder,noglobal=~bref[1,1],bFcalc=bFcalc,Fcalc=Fcalc)
;FWHM=[FWHM_L,FWHM_G]
W=GlobalFWHMntchz(FWHM,typeinfo,keepW,xunoshift,nmix,pderL=pderL,pderG=pderG,bpder=bpder,noglobal=~bref[2,1])
; When "not bref[2,1]" then W == FWHM (however nmix,pderL and pderG are calculated)

; Derived parameters
peaknr=n_elements(xu)
if ptr_valid(typeinfo.pPeakParam) then begin
    *typeinfo.pPeakParam=[reform(xu,1,peaknr),reform(I,1,peaknr),transpose(reform(W,peaknr,2)),reform(nmix,1,peaknr)]
    return
endif

WL=W[0:peaknr-1]
n=-1
cW2=4./(WL*WL)
dW=2/!pi/WL
sw=1./(2*SQRT(2*ALOG(2)))

WG=W[peaknr:2*peaknr-1]
sigmaG=WG*sw
sigmaG2=sigmaG^2.
c=1./(sqrt(2*!dpi)*sigmaG)

; Reliability factors
if bFcalc then begin
    Fcalc=[[temporary(Fcalc)],[make_array(peaknr,n_elements(X),type=size(Fcalc,/type))]]
    ; Fcalc
    ; cte=Iobs/Fobs^2
    ; Icalc
    ; Profiles
endif

for k=0l,peaknr-1 do begin
    ; Normalized Profile function (area=1)
    LorArg=(X-xu[k])^2*cW2[k]+1
    multpderI_L=dW[k]*LorArg^n
    multpderI_G=c[k]*exp(-(X-xu[k])^2/(2*sigmaG2[k]))
    multpderI=multpderI_G*(1-nmix[k]) + multpderI_L*nmix[k]

    ; Profile function
    if keyword_set(SepF) then F=[[F],[multpderI*I[k]]] $
    else F+= multpderI*I[k]
    if bFcalc then Fcalc[k,3:*]=multpderI

    IF bpder THEN BEGIN
        ; dF/duH
        if bref[0] or (bref[2] or (bref[0] and bref[2,1])) then begin
            multpderu_L=((-2)*I[k]*n*cw2[k]*dW[k])*(X-xu[k])*LorArg^(n-1)
            multpderu_G=(I[k]/sigmaG2[k])*multpderI_G*(X-xu[k])
            multpderu=(1-nmix[k])*multpderu_G+nmix[k]*multpderu_L
        endif
        
        ; dF/dW_L ; dF/dW_G
        if bref[2] or (bref[0] and bref[2,1]) then begin
            multpderW_L=multpderu_L*(X-xu[k])/WL[k]-multpderI_L*(I[k]/WL[k])
            multpderW_G=multpderI_G*((X-xu[k])^2/sigmaG2[k]-1)*(sw*I[k]/sigmaG[k])
            multpdern=(multpderI_L-multpderI_G)*I[k]
            
            multpderW_L=multpderW_L*nmix[k]+multpdern*pderL[k]
            multpderW_G=multpderW_G*(1-nmix[k])+multpdern*pderG[k]
        endif
        
        ; dF/duH += dF/dI.dI/duH
        if bref[0] and bref[1,1] then multpderu+=multpderI*keepI[k,1]
        ; dF/duH += dF/dW_L.dW_L/duH + dF/dW_G.dW_G/duH
        if bref[0] and bref[2,1] then multpderu+=multpderW_L*keepW[k,0]+multpderW_G*keepW[k,1]
        
        ; pder to peak or global position parameters
        if bref[0] then begin
            if bref[0,1] then GlobalPospder,keepxu[k,*],off[0],pder,multpderu,toggleu $
            else pder[*,off[0]+k]+=multpderu*keepxu[k,0]; pder xu[k]
        endif
        
        ; pder to peak or global intensity parameters
        if bref[1] then begin
            if bref[1,1] then GlobalIntpder,keepI[k,*],off[1],pder,multpderI,toggleI,typeinfo $
            else pder[*,off[1]+k]+=multpderI*keepI[k,2] ; pder I[k]
        endif

        ; pder to peak or global FWHM parameters
        if bref[2] then begin
            if bref[2,1] then GlobalFWHMtchzpder,keepW[k,*],off[2],pder,multpderW_L,multpderW_G,toggleW $
            else begin
                pder[*,off[2]+k]+=multpderW_L*keepW[k,0]
                pder[*,off[3]+k]+=multpderW_G*keepW[k,1]
            endelse; pder WL[k] and WG[k]
        endif
    endif
endfor

end;pro gfuncOtchz
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncOspvhalf,I,xu,H1,H2,nL,nR,A,F,pder,SWAP=SWAP,MODSPV=MODSPV,$
                                    keepxu,keepI,keepW,keepnmix,keepA,indk,$
                                    toggleu, toggleI, toggleW, togglen, toggleA,$
                                    bref,typeinfo,off,xu0,xu1
;     I: integrated peak area
;     xu: X - peak position u
;     H1: Lorentzian FWHM
;     H2: Gaussian FWHM  (=H1 when MODSPV not set)
;     nL: left mixing parameter
;     nR: right mixing parameter
;     A: asymmetry parameter
;     
;     Calculates 1 half of the split PV. The other half is obtained
;     by calling this function with /SWAP.

case n_elements(off) of
5:     begin
    ; off: u,I,H,[nL,nR],A
    ; SplitPV with global nL and nR
    ; OR
    ; off: u,I,[H1,H2],[nL,nR],A
    ; Mod. SplitPV with global H1 and H2 and HL and HR
    offset=off[[0,1,2,2,3,3,4]]
    endcase
6:     begin
    ; off: u,I,H,nL,nR,A
    ; SplitPV
    ; OR
    ; off: u,I,[H1,H2],nL,nR,A
    ; Mod. SplitPV with global H1 and H2
    ; OR
    ; off: u,I,H1,H2,[nL,nR],A
    ; Mod. SplitPV with global nL and nR
    if keyword_set(MODSPV) and bref[3,1] then offset=off[[0,1,2,3,4,4,5]] $
    else offset=off[[0,1,2,2,3,4,5]]
    endcase
7:     begin
    ; u,I,H1,H2,nL,nR,A
    ; Mod. SplitPV
    offset=off 
    endcase
endcase
; off: u,I,H1,H2,nL,nR,A

if keyword_set(SWAP) then begin
    A=1/A
    tmp=nL
    nL=nR
    nR=tmp
endif

c=sqrt(!pi*alog(2))
fL=nL+c*(1-nL)
fR=nR+c*(1-nR)
Kn=fL+A*fR
K=2*fL/Kn

cWH=2/(1+A)
W1=H1*cWH
W2=H2*cWH

Gnorm=2/W2*sqrt(alog(2)/!pi)*exp(-xu^2*4*alog(2)/W2^2)
LB=1+4*xu^2/W1^2
Lnorm=2/(!pi*W1)/LB

PVnorm=nR*Lnorm+(1-nR)*Gnorm
F=K*PVnorm ; *I must be done in caller!

IF n_elements(pder) eq 0 THEN return

; dSPV/du
if bref[0] then multpderu = K*I*8*xu*(nR/(LB*W1^2)*Lnorm+(1-nR)*alog(2)/W2^2*Gnorm)

; dSPV/dI
;multpderI = F

; dSPV/dH1; dSPV/dH2
if bref[2] or (bref[0] and bref[2,1]) then begin
    pderW1=K*nR*I/W1*(xu^2*8/(LB*W1^2)-1)*Lnorm
    multpderH1 = cWH*pderW1
    
    pderW2=K*(1-nR)*I/W2*(xu^2*8*alog(2)/W2^2-1)*Gnorm
    multpderH2 = cWH*pderW2
endif

; dSPV/dnL; dSPV/dnR
if bref[3] or (bref[0] and bref[3,1]) then begin
    multpdern_L = ((1-c)*(2-K)*I/Kn) * PVnorm
    multpdern_R = (K*I)*(Lnorm-Gnorm) - (A*K*(1-c)*I/Kn) * PVnorm
    
    if keyword_set(SWAP) then begin
        tmp=multpdern_L
        multpdern_L=temporary(multpdern_R)
        multpdern_R=temporary(tmp)
    endif
endif

; dSPV/dA
if bref[4] or (bref[0] and bref[4,1]) then begin
    if keyword_set(SWAP) then mA=A^2 else mA=-1
    if n_elements(pderW1) eq 0 then begin
        pderW1=K*nR*I/W1*(xu^2*8/(LB*W1^2)-1)*Lnorm
        pderW2=K*(1-nR)*I/W2*(xu^2*8*alog(2)/W2^2-1)*Gnorm
    endif
    multpderA = mA*2/(1+A)^2*(H1*pderW1+H2*pderW2) + (mA*fR*K*I/Kn) * PVnorm
endif

; dSPV/du += dSPV/dI.dI/du
if bref[0] and bref[1,1] then multpderu += F*keepI[indk,1]

; dSPV/du += dSPV/dH1.dH1/du + dSPV/dH2.dH2/du
if bref[0] and bref[2,1] then begin
    multpderu += multpderH1*keepW[indk,0]
    multpderu += multpderH2*keepW[indk,1]
endif

; dSPV/du += dSPV/dnL.dnL/du + dSPV/dnR.dnR/du
if bref[0] and bref[3,1] then begin
    multpderu += multpdern_L*keepnmix[indk,0]
    multpderu += multpdern_R*keepnmix[indk,1]
endif

; dSPV/du += dSPV/dA.dA/du
if bref[0] and bref[4,1] then multpderu += multpderA*keepA[indk,0]

; pder to peak or global position parameters
if bref[0] then begin
    if bref[0,1] then GlobalPospder,keepxu[indk,*],offset[0],pder,multpderu,toggleu,i0=xu0,i1=xu1 $
    else pder[xu0:xu1,offset[0]+indk]+=multpderu*keepxu[indk,0]; pder xu[k]
endif

; pder to peak or global intensity parameters
if bref[1] then begin
    if bref[1,1] then GlobalIntpder,keepI[indk,*],offset[1],pder,F,toggleI,typeinfo,i0=xu0,i1=xu1 $
    else pder[xu0:xu1,offset[1]+indk]+=F*keepI[indk,2] ; pder I[k]
endif

; pder to peak or global FWHM parameters
if bref[2] then begin
    if bref[2,1] then begin
        if keyword_set(MODSPV) then begin
            GlobalFWHMmodSplitPVpder,keepW[indk,*],offset[2],pder,multpderH1,multpderH2,toggleW,i0=xu0,i1=xu1
        endif else begin
            GlobalFWHMpder,keepW[indk,*],offset[2],pder,multpderH1,toggleW,i0=xu0,i1=xu1
            GlobalFWHMpder,keepW[indk,*],offset[3],pder,multpderH2,toggleW,i0=xu0,i1=xu1
        endelse
    endif else begin
        pder[xu0:xu1,offset[2]+indk]+=multpderH1*keepW[indk,0] ; pder H1[k]
        pder[xu0:xu1,offset[3]+indk]+=multpderH2*keepW[indk,1] ; pder H2[k]
    endelse
endif
        
; pder to peak or global mixing parameters
if bref[3] then begin
    if bref[3,1] then begin
        GlobalnmixSplitpder,keepnmix[indk,*],offset[4],pder,multpdern_L,multpdern_R,togglen,i0=xu0,i1=xu1
    endif else begin
        pder[xu0:xu1,offset[4]+indk]+=multpdern_L*keepnmix[indk,0] ; pder nmix_L[k]
        pder[xu0:xu1,offset[5]+indk]+=multpdern_R*keepnmix[indk,1] ; pder nmix_R[k]
    endelse
endif

; pder to peak or global asymmetry parameters
if bref[4] then begin
    if bref[4,1] then GlobalAsympder,keepA[indk,*],offset[6],pder,multpderA,toggleA,i0=xu0,i1=xu1 $
    else pder[xu0:xu1,offset[6]+indk]+=multpderA*keepA[indk,0] ; pder I[k]
endif

end;pro gfuncOspvhalf
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncOspv, X, xu, I, FWHM, nmixin, Ain, F, pder, off, bref, toggleu, toggleI, toggleW, togglen, toggleA, typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc, MODSPV=MODSPV

; Initiate
bpder=n_elements(pder) ne 0
nX=n_elements(X)

; Convert global parameters to peakparameters
xu=GlobalPos(xu,typeinfo,keepxu,hkl,toggleu,bpder=bpder,noglobal=~bref[0,1],xunoshift=xunoshift)
peaknr=n_elements(xu)

I=GlobalInt(I,typeinfo,keepI,xunoshift,hkl,bpder=bpder,noglobal=~bref[1,1],bFcalc=bFcalc,Fcalc=Fcalc)

W=GlobalFWHMSplitPV(FWHM,typeinfo,keepW,xunoshift,bpder=bpder,noglobal=~bref[2,1],MODSPV=MODSPV)
W=transpose(reform(W,peaknr,2))

nmix=GlobalnmixSplit(nmixin,keepnmix,xunoshift,bpder=bpder,noglobal=~bref[3,1])
nmix=transpose(reform(nmix,peaknr,2))

A=GlobalAsym(Ain,keepA,xunoshift,bpder=bpder,noglobal=~bref[4,1])

; Derived parameters
if ptr_valid(typeinfo.pPeakParam) then begin
    *typeinfo.pPeakParam=[reform(xu,1,peaknr),reform(I,1,peaknr),W,nmix,reform(A,1,peaknr)]
    return
endif

; Reliability factors
if bFcalc then begin
    Fcalc=[[temporary(Fcalc)],[make_array(peaknr,n_elements(X),type=size(Fcalc,/type))]]
    ; Fcalc
    ; cte=Iobs/Fobs^2
    ; Icalc
    ; Profiles
endif

if bpder then pderkhalf=make_array(1,type=size(pder,/type))

for k=0l,peaknr-1 do begin
    ; Left and right of peak position:
    Xshift=X-xu[k]
    indL=where(Xshift lt 0,ctL,comp=indR,ncomp=ctR)

    Fknorm=make_array(nX,type=size(F,/type))
    
    ; Right side of the split function:
    if ctR ne 0 then begin
        xu0=indR[0]
        xu1=indR[ctR-1]
        gfuncOspvhalf,I[k],Xshift[xu0:xu1],W[0,k],W[1,k],nmix[0,k],nmix[1,k],A[k],Fkhalf,pder,MODSPV=MODSPV,$
                        keepxu,keepI,keepW,keepnmix,keepA,k,$
                        toggleu, toggleI, toggleW, togglen, toggleA,$
                        bref,typeinfo,off,xu0,xu1
        Fknorm[xu0:xu1] = temporary(Fkhalf)
    endif
    
    ; Left side of the split function:
    if ctL ne 0 then begin
        xu0=indL[0]
        xu1=indL[ctL-1]
        gfuncOspvhalf,I[k],Xshift[xu0:xu1],W[0,k],W[1,k],nmix[0,k],nmix[1,k],A[k],Fkhalf,pder,/SWAP,MODSPV=MODSPV,$
                        keepxu,keepI,keepW,keepnmix,keepA,k,$
                        toggleu, toggleI, toggleW, togglen, toggleA,$
                        bref,typeinfo,off,xu0,xu1
        Fknorm[xu0:xu1] = temporary(Fkhalf)
    endif

    ; Profile function
    if keyword_set(SepF) then F=[[F],[Fknorm*I[k]]] $
    else F += Fknorm*I[k]
    if bFcalc then Fcalc[k,3:*]=Fknorm

endfor

end;pro gfuncOspv
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncOpv, X, xu, I, FWHM, nmixin, F, pder, off, bref, toggleu, toggleI, toggleW, togglen, typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc

; Initiate
bpder=n_elements(pder) ne 0

; Convert global parameters to peakparameters
xu=GlobalPos(xu,typeinfo,keepxu,hkl,toggleu,bpder=bpder,noglobal=~bref[0,1],xunoshift=xunoshift)
I=GlobalInt(I,typeinfo,keepI,xunoshift,hkl,bpder=bpder,noglobal=~bref[1,1],bFcalc=bFcalc,Fcalc=Fcalc)
W=GlobalFWHM(FWHM,typeinfo,keepW,xunoshift,bpder=bpder,noglobal=~bref[2,1])
nmix=Globalnmix(nmixin,keepnmix,xunoshift,bpder=bpder,noglobal=~bref[3,1])

; Derived parameters
peaknr=n_elements(xu)
if ptr_valid(typeinfo.pPeakParam) then begin
    *typeinfo.pPeakParam=[reform(xu,1,peaknr),reform(I,1,peaknr),reform(W,1,peaknr),reform(nmix,1,peaknr)]
    return
endif

n=-1
cW2=4./(W*W)
dW=2/!pi/W
sw=1./(2*SQRT(2*ALOG(2)))

sigma=W*sw
sigma2=sigma^2.
c=1./(sqrt(2*!dpi)*sigma)

; Reliability factors
if bFcalc then begin
    Fcalc=[[temporary(Fcalc)],[make_array(peaknr,n_elements(X),type=size(Fcalc,/type))]]
    ; Fcalc
    ; cte=Iobs/Fobs^2
    ; Icalc
    ; Profiles
endif

for k=0l,peaknr-1 do begin
    ; Normalized Profile function (area=1)
    LorArg=(X-xu[k])^2*cW2[k]+1
    multpderI_L=dW[k]*LorArg^n
    multpderI_G=c[k]*exp(-(X-xu[k])^2/(2*sigma2[k]))
    multpderI=multpderI_G*(1-nmix[k]) + multpderI_L*nmix[k]

    ; Profile function
    if keyword_set(SepF) then F=[[F],[multpderI*I[k]]] $
    else F+= multpderI*I[k]
    if bFcalc then Fcalc[k,3:*]=multpderI

    IF bpder THEN BEGIN
        ; dF/duH
        if bref[0] or (bref[2] or (bref[0] and bref[2,1])) then begin
            multpderu_L=((-2)*I[k]*n*cw2[k]*dW[k])*(X-xu[k])*LorArg^(n-1)
            multpderu_G=(I[k]/sigma2[k])*multpderI_G*(X-xu[k])
            multpderu=(1-nmix[k])*multpderu_G+nmix[k]*multpderu_L
        endif
        
        ; dF/dW
        if bref[2] or (bref[0] and bref[2,1]) then begin
            multpderW_L=multpderu_L*(X-xu[k])/W[k]-multpderI_L*(I[k]/W[k])
            multpderW_G=multpderI_G*((X-xu[k])^2/sigma2[k]-1)*(sw*I[k]/sigma[k])
            multpderW=(1-nmix[k])*multpderW_G+nmix[k]*multpderW_L
        endif
        
        ; dF/dn
        if bref[3] or (bref[0] and bref[3,1]) then multpdern=(multpderI_L-multpderI_G)*I[k]
        
        ; dF/duH += dF/dI.dI/duH
        if bref[0] and bref[1,1] then multpderu+=multpderI*keepI[k,1]
        ; dF/duH += dF/dW.dW/duH
        if bref[0] and bref[2,1] then multpderu+=multpderW*keepW[k,0]
        ; dF/duH += dF/dn.dn/duH
        if bref[0] and bref[3,1] then multpderu+=multpdern*keepnmix[k,0]
        
        ; pder to peak or global position parameters
        if bref[0] then begin
            if bref[0,1] then GlobalPospder,keepxu[k,*],off[0],pder,multpderu,toggleu $
            else pder[*,off[0]+k]+=multpderu*keepxu[k,0]; pder xu[k]
        endif

        ; pder to peak or global intensity parameters
        if bref[1] then begin
            if bref[1,1] then GlobalIntpder,keepI[k,*],off[1],pder,multpderI,toggleI,typeinfo $
            else pder[*,off[1]+k]+=multpderI*keepI[k,2] ; pder I[k]
        endif

        ; pder to peak or global FWHM parameters
        if bref[2] then begin
            if bref[2,1] then GlobalFWHMpder,keepW[k,*],off[2],pder,multpderW,toggleW $
            else pder[*,off[2]+k]+=multpderW*keepW[k,0] ; pder W[k]
        endif
        
        ; pder to peak or global mixing parameters
        if bref[3] then begin
            if bref[3,1] then Globalnmixpder,keepnmix[k,*],off[3],pder,multpdern,togglen $
            else pder[*,off[3]+k]+=multpdern*keepnmix[k,0] ; pder nmix[k]
        endif
        
    endif
endfor

end;pro gfuncOpv
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncOlor, X, xu, I, FWHM, F, pder, off, bref, type, toggleu, toggleI, toggleW, typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc

; Initiate
bpder=n_elements(pder) ne 0

; Convert global parameters to peakparameters
xu=GlobalPos(xu,typeinfo,keepxu,hkl,toggleu,bpder=bpder,noglobal=~bref[0,1],xunoshift=xunoshift)
I=GlobalInt(I,typeinfo,keepI,xunoshift,hkl,bpder=bpder,noglobal=~bref[1,1],bFcalc=bFcalc,Fcalc=Fcalc)
W=GlobalFWHM(FWHM,typeinfo,keepW,xunoshift,bpder=bpder,noglobal=~bref[2,1])

; Derived parameters
peaknr=n_elements(xu)
if ptr_valid(typeinfo.pPeakParam) then begin
    *typeinfo.pPeakParam=[reform(xu,1,peaknr),reform(I,1,peaknr),reform(W,1,peaknr)]
    return
endif

case type of
1:    begin
    n=-1
    c=4
    d=2/!pi
    endcase
2:    begin
    n=-2
    c=sqrt(2)-1
    d=4*sqrt(c)/!pi
    c*=4
    endcase
3:    begin
    n=-3./2
    c=2.^(2/3.)-1
    d=sqrt(c)
    c*=4
    endcase
endcase

cW2=c/(W*W)
dW=d/W

; Reliability factors
if bFcalc then begin
    Fcalc=[[temporary(Fcalc)],[make_array(peaknr,n_elements(X),type=size(Fcalc,/type))]]
    ; Fcalc
    ; cte=Iobs/Fobs^2
    ; Icalc
    ; Profiles
endif

; For all peaks
for k=0l,peaknr-1 do begin
    ; Normalized Profile function (area=1)
    LorArg=(X-xu[k])^2*cW2[k]+1
    multpderI=dW[k]*LorArg^n

    ; Profile function
    if keyword_set(SepF) then F=[[F],[I[k]*multpderI]] $
    else F+=I[k]*multpderI
    if bFcalc then Fcalc[k,3:*]=multpderI

    IF bpder THEN BEGIN
        ; dF/duH
        if bref[0] or (bref[2] or (bref[0] and bref[2,1])) then multpderu=((-2)*I[k]*n*cw2[k]*dW[k])*(X-xu[k])*LorArg^(n-1)
        ; dF/dW
        if bref[2] or (bref[0] and bref[2,1]) then multpderW=multpderu*(X-xu[k])/W[k]-multpderI*(I[k]/W[k])

        ; dF/duH += dF/dI.dI/duH
        if bref[0] and bref[1,1] then multpderu+=multpderI*keepI[k,1]
        ; dF/duH += dF/dW.dW/duH
        if bref[0] and bref[2,1] then multpderu+=multpderW*keepW[k,0]
        
        ; pder to peak or global position parameters
        if bref[0] then begin
            if bref[0,1] then GlobalPospder,keepxu[k,*],off[0],pder,multpderu,toggleu $
            else pder[*,off[0]+k]+=multpderu*keepxu[k,0]; pder xu[k]
        endif
        
        ; pder to peak or global intensity parameters
        if bref[1] then begin
            if bref[1,1] then GlobalIntpder,keepI[k,*],off[1],pder,multpderI,toggleI,typeinfo $
            else pder[*,off[1]+k]+=multpderI*keepI[k,2] ; pder I[k]
        endif

        ; pder to peak or global FWHM parameters
        if bref[2] then begin
            if bref[2,1] then GlobalFWHMpder,keepW[k,*],off[2],pder,multpderW,toggleW $
            else pder[*,off[2]+k]+=multpderW*keepW[k,0] ; pder W[k]
        endif
    endif
endfor

end;pro gfuncOlor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncOgauss, X, xu, I, FWHM, F, pder, off, bref, toggleu, toggleI, toggleW, typeinfo ,SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc

; Initiate
bpder=n_elements(pder) ne 0
;if bpder then stop

; Convert global parameters to peakparameters
xu=GlobalPos(xu,typeinfo,keepxu,hkl,toggleu,bpder=bpder,noglobal=~bref[0,1],xunoshift=xunoshift)
I=GlobalInt(I,typeinfo,keepI,xunoshift,hkl,bpder=bpder,noglobal=~bref[1,1],bFcalc=bFcalc,Fcalc=Fcalc)
W=GlobalFWHM(FWHM,typeinfo,keepW,xunoshift,bpder=bpder,noglobal=~bref[2,1])

; Derived parameters
peaknr=n_elements(xu)
if ptr_valid(typeinfo.pPeakParam) then begin
    *typeinfo.pPeakParam=[reform(xu,1,peaknr),reform(I,1,peaknr),reform(W,1,peaknr)]
    return
endif

sw=1./(2*SQRT(2*ALOG(2)))
sigma=W*sw
sigma2=sigma^2.
c=1/(sqrt(2*!dpi)*sigma)

; Reliability factors
if bFcalc then begin
    Fcalc=[[temporary(Fcalc)],[make_array(peaknr,n_elements(X),type=size(Fcalc,/type))]]
    ; Fcalc
    ; cte=Iobs/Fobs^2
    ; Icalc
    ; Profiles
endif

; For all peaks
for k=0l,peaknr-1 do begin
    ; Normalized Profile function (area=1)
    multpderI=c[k]*exp(-(X-xu[k])^2/(2*sigma2[k]))

    ; Profile function
    if keyword_set(SepF) then F=[[F],[I[k]*multpderI]] $
    else F+=I[k]*multpderI
    if bFcalc then Fcalc[k,3:*]=multpderI

    IF bpder THEN BEGIN
        ; dF/duH
        if bref[0] then multpderu=(I[k]/sigma2[k])*multpderI*(X-xu[k])
        ; dF/dW
        if bref[2] or (bref[0] and bref[2,1]) then multpderW=multpderI*((X-xu[k])^2/sigma2[k]-1)*(sw*I[k]/sigma[k])

        ; dF/duH += dF/dI.dI/duH
        if bref[0] and bref[1,1] then multpderu+=multpderI*keepI[k,1]
        ; dF/duH += dF/dW.dW/duH
        if bref[0] and bref[2,1] then multpderu+=multpderW*keepW[k,0]
        
        ; pder to peak or global position parameters
        if bref[0] then begin
            if bref[0,1] then GlobalPospder,keepxu[k,*],off[0],pder,multpderu,toggleu $
            else pder[*,off[0]+k]+=multpderu*keepxu[k,0]; pder xu[k]
        endif

        ; pder to peak or global intensity parameters
        if bref[1] then begin
            if bref[1,1] then GlobalIntpder,keepI[k,*],off[1],pder,multpderI,toggleI,typeinfo $
            else pder[*,off[1]+k]+=multpderI*keepI[k,2] ; pder I[k]
        endif

        ; pder to peak or global FWHM parameters
        if bref[2] then begin
            if bref[2,1] then GlobalFWHMpder,keepW[k,*],off[2],pder,multpderW,toggleW $
            else pder[*,off[2]+k]+=multpderW*keepW[k,0] ; pder W[k]
        endif
    endif
endfor

end;pro gfuncOgauss
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TypexExtractVar,ptr,j,peaknr,off,A,struc,_u,offu,toggleu

_u='dummy'

if (*ptr).codeproc[j,0] then begin
    bAref=(*ptr).codeproc[j,2]
    toggleoff=0

    global=(*ptr).codeproc[j,1] ne 0

    ; this part is different from similar codes below!!!!!
    ARRind=*(*ptr).paramind[j,global] ; column indices in (*p)
    nARRind=n_elements(ARRind) ; n columns
    
    offcol=0
    if global then begin
        if ~bAref then n=[0,nARRind] else $
        n=histogram((*ptr).togglefix[ARRind],min=0,max=1,REVERSE_INDICES=rtoggle) ;[nrefined,nfixed]
        if n[0] eq 0 then begin
            ARRind0=[-1]
            ARRind1=lindgen(nARRind)
        endif else begin
            ARRind0=rtoggle[rtoggle[0]:rtoggle[1]-1]
            if n[1] ne 0 then ARRind1=rtoggle[rtoggle[1]:rtoggle[2]-1] else ARRind1=[-1]
        endelse
    endif else begin
        nARRind*=peaknr
        if bAref then begin
            ARRind1=where((*ptr).peakhelpparam[ARRind],ct2,comp=ARRind0,ncomp=ct3) ;2>fixed, 3>refined
            offcol=lindgen(ct3)*peaknr
            n=[ct3,ct2]*peaknr ;[nrefined,nfixed]
            if ct2 ne 0 and ct3 ne 0 then add0=ARRind0[0] gt ARRind1[0] else add0=0b
            if n[0] eq 0 then ARRind0=[-1] else ARRind0=lindgen(n[0])
            if n[1] eq 0 then ARRind1=[-1] else ARRind1=lindgen(n[1])
            if add0 then ARRind0+=n[1] else ARRind1+=n[0]
        endif else begin
            n=[0,nARRind] ;[nrefined,nfixed]
            ARRind0=[-1]
            ARRind1=lindgen(n[1])
        endelse
    endelse

    _u=make_array(nARRind,type=size(A,/type))
    if global then toggleu=intarr(nARRind,3)
    
    if n[0] ne 0 then begin
        _u[ARRind0]=A[off[0]:off[0]+n[0]-1]
        if global then begin
            toggleu[ARRind0,0]=1 ; used gobal parameters: 1
            toggleu[ARRind0,1]=indgen(n[0]) ; offset for used global parameters
            ; toggleu[*,2] indices of connected, global parameters
            toggleoff+=n[0]
        endif
    endif
    if n[1] ne 0 then begin
        _u[ARRind1]=struc.Afix[off[1]:off[1]+n[1]-1]
    endif
    offu=off[0]+offcol; Offset in pder
    off+=n ; Offset in A and Afix
    
    ; Unit cell
    if (*ptr).type ne 10 and j eq 0 then begin
        _u[1:*]=celparamsetcon(_u[1:*],*(*ptr).pstrinfo,connect=connect)
        if n_elements(toggleu) ne 0 then toggleu[*,2]=[0,connect+1]
    endif
        
    ; Extra peak parameters from afix and/or A
    if global and (*ptr).codeproc[j,6] ne 0 then begin
        ARRind=*(*ptr).paramind[j,0] ; column indices in (*p)
        nARRind=n_elements(ARRind) ; n columns
        nARRind*=peaknr
        if bAref and (*ptr).codeproc[j,6] eq 2 then begin
            ARRind1=where((*ptr).peakhelpparam[ARRind],ct2,comp=ARRind0,ncomp=ct3) ;2>fixed, 3>refined
            n=[ct3,ct2]*peaknr ;[nrefined,nfixed]
            if ct2 ne 0 and ct3 ne 0 then add0=ARRind0[0] gt ARRind1[0] else add0=0b
            if n[0] eq 0 then ARRind0=[-1] else ARRind0=lindgen(n[0])
            if n[1] eq 0 then ARRind1=[-1] else ARRind1=lindgen(n[1])
            if add0 then ARRind0+=n[1] else ARRind1+=n[0]
        endif else begin
            n=[0,nARRind] ;[nrefined,nfixed]
            ARRind0=[-1]
            ARRind1=lindgen(n[1])
        endelse
        
        _uadd=make_array(nARRind,type=size(A,/type))
        if n[0] ne 0 then _uadd[ARRind0]=A[off[0]:off[0]+n[0]-1]
        if n[1] ne 0 then _uadd[ARRind1]=struc.Afix[off[1]:off[1]+n[1]-1]
        _u=[_u,_uadd]
            
        off+=n
    endif

    ; ASU (no info in Afix)
    if (*ptr).type eq 12 and j eq 1 then begin
        struc.typeinfo.pasuinfo=(*ptr).asu.I
    
        ; xyzfit
        ind=(*(*ptr).preverse_indices)[0:(*ptr).nasupos-1]
        xyz=transpose((*(*ptr).asu.I)[ind].wyck.xyz)
    
        xyz999=xyz ne 999
        ind2=where(xyz999,n)
        if n ne 0 then begin
            ind_=where(xyz999 and ~transpose((*(*ptr).asu.I)[ind].wyck.togglefix),ctref)
            xyz999=bytarr(n)
            if ctref ne 0 and bAref then begin
                xyz[ind_]=A[off[0]:off[0]+ctref-1]
                xyz999[ind_]=1
                off[0]+=ctref
            endif
            _u=[_u ,xyz[ind2]]
    
            ind=intarr(n)
            if ctref ne 0 then ind[ind_]=indgen(ctref)+toggleoff
            toggleu=[toggleu,[[xyz999],[ind],[ind]]]
            toggleoff+=ctref
        endif
    
        ; SOF and B
        other=[[(*(*ptr).asu.I).SOF],[(*(*ptr).asu.I).B]]
        n=2*(*ptr).nasuatoms
        togg=bytarr(n)
        ind_=where(~transpose((*(*ptr).asu.I).togglefix),ctref)
        if ctref ne 0 and bAref then begin
            other[ind_]=A[off[0]:off[0]+ctref-1]
            togg[ind_]=1
            off[0]+=ctref
        endif
        _u=[_u,reform(other,n)]
        ind=intarr(n)
        if ctref ne 0 then ind[ind_]=indgen(ctref)+toggleoff
        toggleu=[toggleu,[[togg],[ind],[ind]]]
    endif

endif

end;pro TypexExtractVar
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncTypex, X, A, F, struc, pder, ptr, off, SepF=SepF, finalize=finalize, initialize=initialize, Rfactorinfo=Rfactorinfo

peaknr=(*ptr).peaknr
if peaknr eq 0 then return
if (*ptr).type eq 12 then $
if (*ptr).nasupos eq 0 then return

ind=where((*ptr).codeproc[*,0])
bref=(*ptr).codeproc[ind,*]
bref=bref[*,[2,1]] ; refined?, global?

struc.typeinfo.type=(*ptr).type
if (*ptr).type eq 12 then struc.typeinfo.pstrinfo=(*ptr).pstrinfo

; ----Calculate reliability factors? (prepare for LeBail)----
bFcalc=keyword_set(finalize) and (*ptr).type ne 10 ; R-factors
bFcalc or= keyword_set(initialize) and bref[0,1] and bref[1,1] ; Initialize scaling factor
if (*ptr).type eq 11 then begin ; Le Bail fitting
    bLeBail=ExtractbitCode(ptr,2) eq 1
    bFcalc or= bLeBail and (n_elements(pder) eq 0) and n_elements(SepF) eq 0
endif

; ----Extract variables----
; u
TypexExtractVar,ptr,0,peaknr,off,A,struc,__u,offu,toggleu
; I
TypexExtractVar,ptr,1,peaknr,off,A,struc,__I,offI,toggleI
; W
TypexExtractVar,ptr,2,peaknr,off,A,struc,__W,offW,toggleW
; n
TypexExtractVar,ptr,3,peaknr,off,A,struc,__n,offn,togglen
; R
TypexExtractVar,ptr,4,peaknr,off,A,struc,__R,offR,toggleR
; A
TypexExtractVar,ptr,5,peaknr,off,A,struc,__A,offA,toggleA
; Add ASU to _I

; ----Choose type---
type=ExtractbitCode(ptr,0)

; ----Add reliability factor info----
if bFcalc then begin
    if n_elements(Rfactorinfo) eq 0 then begin
        jlast=-1
        Rfactorinfo={Profiles:ptr_new(),list:ptr_new({ptr:ptr,offA:offI[0],i0:0l,i1:-1l,bref:(*ptr).codeproc[1,2]})}
    endif else begin
        jlast=n_elements(*Rfactorinfo.list)-1
        i1=((*Rfactorinfo.list)[jlast]).i1
        *Rfactorinfo.list=[*Rfactorinfo.list,{ptr:ptr,offA:offI[0],i0:i1+1,i1:i1,bref:(*ptr).codeproc[1,2]}]
    endelse
endif

; ----Add contribution to the fit----
for i=0l,struc.typeinfo.maxtable do begin
    _u=__u
    _I=__I
    _W=__W
    _n=__n
    _R=__R
    _A=__A
    struc.typeinfo.itable=i
    struc.typeinfo.lambda=struc.typeinfo.ltable[0,i]
    case type of
    0:    gfuncOgauss, X, _u, _I, _W, F, pder, [offu,offI,offW], bref, toggleu, toggleI, toggleW, struc.typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc
    1:     gfuncOlor, X, _u, _I, _W, F, pder, [offu,offI,offW], bref, type, toggleu, toggleI, toggleW, struc.typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc
    2:     gfuncOlor, X, _u, _I, _W, F, pder, [offu,offI,offW], bref, type, toggleu, toggleI, toggleW, struc.typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc
    3:     gfuncOlor, X, _u, _I, _W, F, pder, [offu,offI,offW], bref, type, toggleu, toggleI, toggleW, struc.typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc
    4:    gfuncOpv, X, _u, _I, _W, _n, F, pder, [offu,offI,offW,offn], bref, toggleu, toggleI, toggleW, togglen, struc.typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc
    5:    gfuncOtchz, X, _u, _I, _W, F, pder, [offu,offI,offW], bref, toggleu, toggleI, toggleW, struc.typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc
    6:    gfuncOpvII, X, _u, _I, _W, _R, F, pder, [offu,offI,offW,offR], bref, toggleu, toggleI, toggleW, toggleR, struc.typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc
    7:  gfuncOspv, X, _u, _I, _W, _n, _A, F, pder, [offu,offI,offW,offn,offA], bref, toggleu, toggleI, toggleW, togglen, toggleA, struc.typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc
    8:  gfuncOspv, X, _u, _I, _W, _n, _A, F, pder, [offu,offI,offW,offn,offA], bref, toggleu, toggleI, toggleW, togglen, toggleA, struc.typeinfo, SepF=SepF, bFcalc=bFcalc, Fcalc=Fcalc, /MODSPV
    else:
    endcase
    
    if bFcalc then begin
        ; Fcalc: # colomns = # peaks
        ;  row0 = |FH|               (or 1 for PD)
        ;  row1 = S.mH.LGP.Ii/I0     (or 1 for PD)
        ;  row2 = IH                 (not PD: row2 = row0^2*row1)
        ;  row3->row3+npts = normalized peak profiles
        ;  
        ; Ii/I0: source emission line ratio
        ;
        ; window & plot,Fcalc[*,2],psym=1 & oplot,Fcalc[*,0]*Fcalc[*,0]*Fcalc[*,1],color=100
        ; window & plot,Fcalc[0,3:*] & for i=1,peaknr-1 do oplot,Fcalc[i,3:*]
        
        (*Rfactorinfo.list)[jlast+1].i1+=(size(Fcalc))[1]
        if ptr_valid(Rfactorinfo.Profiles) then *Rfactorinfo.Profiles=[*Rfactorinfo.Profiles,temporary(Fcalc)] $
        else Rfactorinfo.Profiles=ptr_new(temporary(Fcalc))
    endif
endfor


end;pro gfuncTypex
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;pro gfuncFitModel, X, Ain, F, struc, pder, SepF=SepF, finalize=finalize, initialize=initialize ;marquardtc
;pro gfuncFitModel, X, Ain, F, pder, _EXTRA=struc, SepF=SepF, finalize=finalize, initialize=initialize ;mpfit
pro gfuncFitModel, X, Ain, F, pder, struc=struc, SepF=SepF, finalize=finalize, initialize=initialize ;NLLS
; If SepF then the first row the background

; ----Initialize----
nx=n_elements(X)
F = make_array(nx,type=size(Ain,/type))
;bpder=N_PARAMS() GE 5 ;marquardtc
bpder=N_PARAMS() GE 4 ; mpfit,NLLS
IF bpder THEN BEGIN
    nA=n_elements(Ain)
    pder=make_array(nX,nA,type=size(Ain,/type))
ENDIF

A=Ain
off=intarr(2) ;A, Afix
ptrFit=struc.ptr

; ----Background----
Fbkg=F
bbkg=0b
if (*(*ptrFit).nparam)[0] ne 0 then begin
    bbkg=1b
    for k=0l,(*ptrFit).backparamIO[0] do Fbkg+=A[k]*(*(*ptrFit).POL[0])[*,k]
    IF bpder THEN $
        for k=0l,(*ptrFit).backparamIO[0] do pder[*,k]=(*(*ptrFit).POL[0])[*,k]
    off[0]+=(*(*ptrFit).nparam)[0]
    *(*ptrFit).yb=Fbkg
endif else if ~(*ptrFit).fitcontrol[4] then begin
    bbkg=1b
    Fbkg[*]=*(*ptrFit).yb
endif

if keyword_set(SepF) and bbkg then F=temporary(Fbkg)
offorg=off

; ----Highlight peak for plot----
ind=where(tag_names(struc) eq 'HIGHLIGHTPEAK',ct)
bthick=ct eq 1 and keyword_set(SepF)

; ----Loop over all models----
ptr=(*ptrFit).next
while PTR_VALID(ptr) do begin
    if (*ptr).include then begin
    
        ; highlight peak
        if bthick then $
            if ptr eq struc.highlightpeak.ptr then begin
                s=DimSize(F,2)
                struc.thicki=s[1]+struc.highlightpeak.i
            endif
            
        gfuncTypex, X, A, F, struc, pder, ptr, off, SepF=SepF, finalize=finalize, initialize=initialize, Rfactorinfo=Rfactorinfo
    endif
    ptr=(*ptr).next
endwhile

; ----Construct total model----
; SepF==1 : F=[[bkg],[peak1],[peak2],[...]]
; SepF==0 : F=peak1+peak2+...
; 
; peak1, peak2, ... must be multiplied by some f(X)
; and bkg must be added to F when SepF==0

case struc.typeinfo.LGPtype of
;3:    intmult=1/sin(X*!pi/360)
; It should be already removed in 2D
else:
endcase

if n_elements(intmult) ne 0 then begin
    if keyword_set(SepF) then begin
        s=DimSize(F,2)-[0,1]
        F[*,1:*]*=rebin(intmult,s,/sample)
    endif else begin
        F*=intmult
    endelse
    
    IF bpder THEN begin
        s=DimSize(pder,2)-[0,offorg[0]]
        pder[*,offorg[0]:*]*=rebin(intmult,s,/sample)
    endif
endif

; ----Reliability factors----
if n_elements(Rfactorinfo) ne 0 then begin
    jlast=n_elements(*Rfactorinfo.list)-1
    ratio=F ; model without background
    ind=where(abs(ratio) gt 1e-7,ny)
    
    if ny ne 0 then begin
        Fcalc=(*Rfactorinfo.Profiles)[*,0] ; |FH|
        Fobs=(*Rfactorinfo.Profiles)[*,1] ; S.mH.LGP         containes const=Iobs/Fobs^2
        Icalc=(*Rfactorinfo.Profiles)[*,2] ; IH
        Iobs=Icalc ; Icalc
        
        tmp=struc.y
        if bbkg then tmp-=Fbkg
        ratio=(temporary(tmp)>1)[ind]/ratio[ind]*(Binwidths(x))[ind]  ; y/ycalc*dx   (both without background)
        
        nref=(*Rfactorinfo.list)[jlast].i1+1
        ratio=rebin(reform(abs(ratio),1,ny),nref,ny,/sample)
        *Rfactorinfo.Profiles=((*Rfactorinfo.Profiles)[*,ind+3]) ; keep normalized profiles of all peaks
        ; i=760 & window & plot,(*Rfactorinfo.Profiles)[i,*] & print,total((*Rfactorinfo.Profiles)[i,*],2)*0.025
        *Rfactorinfo.Profiles*=temporary(ratio) ; multiply each normalized profile with y/ycalc*dx: normprofile*y/ycalc*dx
        
        nlines=struc.typeinfo.maxtable+1
        ; Loop over the phases
        for j=0l,jlast do begin
            i0=(*Rfactorinfo.list)[j].i0
            i1=(*Rfactorinfo.list)[j].i1
            nref=i1-i0+1 ; number of peaks of phase j
            nref/=nlines ; number of peaks of phase j / number of source emission lines
            i1=nref+i0-1
            ptr=(*Rfactorinfo.list)[j].ptr
            
            if keyword_set(initialize) then begin
                ; Fobs = Fcalc
                ; Snew = Sold * mean(Iobs/Icalc)
                if (*Rfactorinfo.list)[j].bref then begin
                    tmp=abs(initialize)
                    Ain[(*Rfactorinfo.list)[j].offA ]*=tmp*median(total((*Rfactorinfo.Profiles)[i0:i1,*],2))
                    initialize=-tmp
                endif
            endif else begin
                Iobs[i0:i1]*=total((*Rfactorinfo.Profiles)[i0:i1,*],2) ; Iobs=Icalc*normprofile*y/ycalc*dx
                ind2=where(Fobs[i0:i1] ne 0,ct)
                if ct ne 0 then begin
                    ind2+=i0
                    Fobs[ind2]=Iobs[ind2]/Fobs[ind2]
                endif
                
                if keyword_set(finalize) then begin
                    Fobs[i0:i1]=sqrt(Fobs[i0:i1]>0)
                    p=(*(*ptr).factors)[(*struc.ptr).icurrent]
                    rb_dev=total(Iobs[i0:i1])
                    rf_dev=total(Fobs[i0:i1])
                    (*p)[0]=0
                    (*p)[1]=0
                    if rb_dev ne 0 then (*p)[0]=100*total(abs(Iobs[i0:i1]-Icalc[i0:i1]))/rb_dev<100 ;Rb
                    if rf_dev ne 0 then (*p)[1]=100*total(abs(Fobs[i0:i1]-Fcalc[i0:i1]))/rf_dev<100 ;Rf
                endif else if ExtractbitCode(ptr,2) eq 1 then begin
                    ; LeBail
                    if (*Rfactorinfo.list)[j].bref then begin
                        tmp=(*Rfactorinfo.list)[j].offA 
                        ; +1: Skip scaling factor
                        Ain[tmp+1:tmp+nref]=Fobs[i0:i1]
                    endif
                endif
            endelse
        endfor
    endif
    
    (*Rfactorinfo.list)[0:jlast].ptr=ptr_new()
    heap_free,Rfactorinfo
endif

; ----Add background to fit if necessary----
if ~keyword_set(SepF) and bbkg then F+=temporary(Fbkg)
F=float(F) ; is already float ...

; ----For display----
if ((*ptrFit).fitcontrol[2] eq 2) and (not keyword_set(SepF)) and widget_info(struc.refreshev.top,/valid) and ~bpder then begin
    (*(*(*ptrFit).yfit)[(*ptrFit).icurrent])[(*ptrFit).backranind[0]:(*ptrFit).backranind[1]]=F
    if (*ptrFit).fitcontrol[4] then (*(*(*ptrFit).yfit)[(*ptrFit).icurrent])[(*ptrFit).backranind[0]:(*ptrFit).backranind[1]]+=(*(*ptrFit).yb)
    ;widget_control,struc.refreshev.top,send_event=struc.refreshev
    widget_control,struc.refreshev.top,get_uvalue=tmp
    RefreshDisplayCHI,tmp
    ;wait,0.5
endif

end;pro gfuncFitModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MakeLabelsFromTypex,ptr,Labels,std=std

std=keyword_set(std)
if std then stradd='_e.s.d.' else stradd=''
indi=where((*ptr).codeproc[*,0],ct)
peaknr=(*ptr).peaknr
if peaknr eq 0 then return
if (*ptr).type eq 12 then $
if (*ptr).nasupos eq 0 then return

bPhysicalPhaseInfo=checkPhysicalPhaseInfo(ptr)

for ii=0,ct-1 do begin
    bAref=(*ptr).codeproc[indi[ii],2]

    if bAref then begin
    
    global=(*ptr).codeproc[indi[ii],1] ne 0
    
    ; point to global or peakparam
    if global then begin
        p=(*ptr).globalparamIO.I
        p1=(*ptr).globallabels
    endif else begin
        p=(*ptr).peakparamIO.I
        p1=(*ptr).peaklabels
    endelse
    
    ; Add to A or/and Afix,constr
    ARRind=*(*ptr).paramind[indi[ii],global] ; column indices in (*p)
    nARRind=n_elements(ARRind) ; n columns
    if global then begin
        n=histogram((*ptr).togglefix[ARRind],min=0,max=1,REVERSE_INDICES=rtoggle) ;[nrefined,nfixed]
        if n[0] eq 0 then begin
            ARRind0=[-1]
        endif else begin
            ARRind0=ARRind[rtoggle[rtoggle[0]:rtoggle[1]-1]]
        endelse
    endif else begin
        ARRind1=where((*ptr).peakhelpparam[ARRind],ct2,comp=ARRind0,ncomp=ct3) ;2>fixed, 3>refined
        n=[ct3,ct2]*peaknr ;[nrefined,nfixed]
        if n[0] ne 0 then ARRind0=ARRind[ARRind0]
    endelse
    
    ; Peak indices
    if global then strId='' $
    else strId='_'+stringr(lindgen(1,(*ptr).peaknr))
    
    if n[0] ne 0 then begin ; only A, not Afix

        ; add to Labels
        tmp=string((*p)[ARRind0,*]) ; get columns
        for j=0l,n_elements(ARRind0)-1 do tmp[j,*]=(*ptr).name+strId+' '+p1[ARRind0[j]]+strId+stradd+' 1DFit '+tmp[j,*]
        if ~global then begin ; if not global: max sure that the column are not mixed
            tmp=transpose(tmp)
            tmp=reform(tmp,n[0],/overwrite)
        endif
    
        Labels=[Labels,tmp]
    endif
        
    ; Refined parameters but not in a least squares sense
    if global and (*ptr).codeproc[indi[ii],6] eq 2 then begin
        p=(*ptr).peakparamIO.I
        p1=(*ptr).peaklabels
            
        ARRind=*(*ptr).paramind[indi[ii],0] ; column indices in (*p)
        nARRind=n_elements(ARRind)
        ARRind1=where((*ptr).peakhelpparam[ARRind],ct2,comp=ARRind0,ncomp=ct3) ;2>fixed, 3>refined
        n=[ct3,ct2]*peaknr ;[nrefined,nfixed]
            
        if n[0] ne 0 then begin ; only A, not Afix
            ARRind0=ARRind[ARRind0]
                
            ; add to Labels
            tmp=string((*p)[ARRind0,*]) ; get columns
            for j=0l,n_elements(ARRind0)-1 do tmp[j,*]=(*ptr).name+strId+' '+p1[ARRind0[j]]+strId+stradd+' 1DFit '+tmp[j,*]
            ; make sure that the column are not mixed
            tmp=transpose(tmp)
            tmp=reform(tmp,n[0],/overwrite)
        
            Labels=[Labels,tmp]
        endif
    endif
    
    ; ASU labels
    if (*ptr).type eq 12 and (indi[ii] eq 1) then begin
        p=(*ptr).asu.I
        p1=(*ptr).asulabels
    
        ind=(*(*ptr).preverse_indices)[0:(*ptr).nasupos-1]
        togglefix=(*p)[ind].wyck.togglefix
            
        ; xyzfit
        asunames='pos'+stringr(indgen((*ptr).nasupos))+'_'+(*p)[ind].name
            
        xyz=(*p)[ind].wyck.xyz
        for i=0l,2 do begin
            ind=where(xyz[i,*] ne 999 and ~togglefix[i,*],ct)
            if ct ne 0 then begin
                indadd=asunames[ind];stringr(indgen(ct))
                tmp=(*ptr).name+' '+p1[((*ptr).indasused)[i]]+indadd+stradd+' 1DFit '+stringr(xyz[i,ind])
                Labels=[Labels,tmp]
            endif
        endfor
    
        ; SOF
        ind=where(~(*p).togglefix[0,*],ct)
        asunames='pos'+stringr(indgen((*ptr).nasuatoms))+'_'+(*p).name
        if ct ne 0 then begin
            indadd=asunames[ind];stringr(indgen(ct))
            tmp=(*ptr).name+' '+p1[((*ptr).indasused)[3]]+indadd+stradd+' 1DFit '+stringr((*p)[ind].SOF)
            Labels=[Labels,tmp]
        endif
    
        ; Biso
        ind=where(~(*p).togglefix[1,*],ct)
        if ct ne 0 then begin
            indadd=asunames[ind];stringr(indgen(ct))
            tmp=(*ptr).name+' '+p1[((*ptr).indasused)[4]]+indadd+stradd+' 1DFit '+stringr((*p)[ind].B)
            Labels=[Labels,tmp]
        endif
    endif
    
    ; Scaling factor derived parameters
    if ~std and (indi[ii] eq 1) and bPhysicalPhaseInfo then begin
        tmp=(*ptr).name+' '+strcompress((*ptr).propname,/remove_all)+' 1DFit -'
        Labels=[Labels,tmp]
    endif
        
    endif ; bAref
endfor

; Add statistics
if ptr_valid((*ptr).factors) and ~std then begin
    Labels=[Labels,(*ptr).name+[' Rb(%) 1DFit -',' Rf(%) 1DFit -']]
endif

end;pro MakeLabelsFromTypex
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MakeAFromTypex,ptr,A,nout,icurrent,constr,Afix,refinedA=refinedA,std=std,batchexport=batchexport

; In A and Afix: all positions next to eachother, all FWHM's next...
; refinedA=0: kset=0b    nkset2=1b    std=0b (get the initial A + Afix + constr)
; refinedA=1: kset=1b    nkset2=1b    std=0b (get the refined A + Afix; Add standard deviation)
; refinedA=2: kset=1b    nkset2=0b    std=.. (get the refined A; Add standard deviation)

batchexport=keyword_set(batchexport)
kset=keyword_set(refinedA)
std=keyword_set(std)
if kset then begin
    nkset2=refinedA eq 1
    if nkset2 then std=0b
endif else begin
    nkset2=1b
    std=0b
endelse
nout=0
bits=[1,1,1,1,16]
nbits=5

; Add peakparameters
ind=where((*ptr).codeproc[*,0],ct)
peaknr=(*ptr).peaknr
if peaknr eq 0 then return
if (*ptr).type eq 12 then $
if (*ptr).nasupos eq 0 then return

bPhysicalPhaseInfo=checkPhysicalPhaseInfo(ptr)

for i=0l,ct-1 do begin
    bAref=(*ptr).codeproc[ind[i],2]

    if nkset2 or bAref then begin

    global=(*ptr).codeproc[ind[i],1]
    
    ; Point to global or peakparam
    if global then begin
        if kset then p=std?(*(*ptr).globalparamIO.SD)[icurrent]:(*(*ptr).globalparamIO.O)[icurrent] $
        else p=(*ptr).globalparamIO.I
    endif else begin
        if kset then p=std?(*(*ptr).peakparamIO.SD)[icurrent]:(*(*ptr).peakparamIO.O)[icurrent] $
        else p=(*ptr).peakparamIO.I
    endelse
    
    ; Add to A or/and Afix,constr
    ARRind=*(*ptr).paramind[ind[i],global] ; column indices in (*p)
    nARRind=n_elements(ARRind) ; n columns
    if global then begin
        if ~bAref then n=[0,nARRind] else $
        n=histogram((*ptr).togglefix[ARRind],min=0,max=1,REVERSE_INDICES=rtoggle) ;[nrefined,nfixed]
        if n[0] eq 0 then begin
            ARRind0=[-1]
            ARRind1=ARRind
        endif else begin
            ARRind0=ARRind[rtoggle[rtoggle[0]:rtoggle[1]-1]]
            if n[1] ne 0 then ARRind1=ARRind[rtoggle[rtoggle[1]:rtoggle[2]-1]] else ARRind1=[-1]
        endelse
    endif else begin
        if bAref then begin
            ARRind1=where((*ptr).peakhelpparam[ARRind],ct2,comp=ARRind0,ncomp=ct3) ;2>fixed, 3>refined
            n=[ct3,ct2]*peaknr ;[nrefined,nfixed]
            if n[0] ne 0 then ARRind0=ARRind[ARRind0]
            if n[1] ne 0 then ARRind1=ARRind[ARRind1]
        endif else begin
            n=[0,nARRind*peaknr] ;[nrefined,nfixed]
            ARRind0=[-1]
            ARRind1=ARRind
        endelse
    endelse
    
    ; Add to A, constr
    if n[0] ne 0 then begin
        Aadd=reform(transpose((*p)[ARRind0,*]),n[0])
        
        ; Stabilize fit parameters
        stabmult=1
        if global and ind[i] eq 1 and ~batchexport then stabmult=10.^(*ptr).ScalingMultiplier ; Scaling factor
        Aadd*=stabmult

        ; Constraints
        if (~kset) then begin
            if global then pconstr=(*ptr).globalparamIO.constr $
            else pconstr=(*ptr).peakparamIO.constr
            constradd=MakeAInfo(n[0])
    
            bconstr=(*ptr).codeproc[ind[i],3:5]
            if total(bconstr,/pres) ne 0 then begin
                breject=replicate((*ptr).codeproc[ind[i],7],1,n[0])+1
    
                if bconstr[0] ne 0 then begin
                    tmp=(*pconstr)[ARRind0,0]*stabmult
                    constradd.limited[0,*]=breject
                        
                    if ~global then begin ; if not global: make sure that the columns are not mixed
                        tmp#=replicate(1b,peaknr)
                        tmp=transpose(tmp)
                    endif
    
                    constradd.limits[0,*]=reform(tmp,1,n[0])
                endif
                if bconstr[1] ne 0 then begin
                    tmp=(*pconstr)[ARRind0,1]*stabmult
                    constradd.limited[1,*]=breject
    
                    if ~global then begin ; if not global: make sure that the columns are not mixed
                        tmp#=replicate(1b,peaknr)
                        tmp=transpose(tmp)
                    endif
    
                    constradd.limits[1,*]=reform(tmp,1,n[0])
                endif
                if bconstr[2] ne 0 then begin
                    tmp0=(*pconstr)[ARRind0,2]*stabmult
                    tmp1=tmp0
                    constradd.limited=[breject,breject]
    
                    if ~global then begin ; if not global: make sure that the columns are not mixed
                        tmp0#=replicate(1b,peaknr)
                        tmp0=transpose(tmp0)
                        tmp1#=replicate(1b,peaknr)
                        tmp1=transpose(tmp1)
                    endif
    
                    constradd.limits[0,*]=reform(Aadd-tmp0,1,n[0])
                    constradd.limits[1,*]=reform(Aadd+tmp1,1,n[0])
                endif
    
            endif else begin
                ; No constraints given by user: default constraints for some parameters
                pres=floatpres()
    
                case ind[i] of
                ; Global and peak position >= 0, except for zero
                0:    begin
                    constradd.limited[0,*]=1
                    if global then constradd[0].limited[0]=0
                    endcase
    
                ; Global and peak I >= 0
                1:     begin
                    constradd.limited[0,*]=1
                    endcase
    
                ; FWHM > 0
                2:    begin
                    if global then begin
                        tmplt0=where(Arrind0 eq Arrind[1],cttmplt0,comp=tmpgt0,ncomp=cttmpgt0) ; V<0
                        if cttmplt0 ne 0 then constradd[tmplt0].limited[1]=1 ;<=0
                        if cttmpgt0 ne 0 then constradd[tmpgt0].limited[0]=1 ;>=0
                    endif else begin
                        constradd.limited[0,*]=1
                        constradd.limits[0,*]=pres
                    endelse
                    endcase
    
                ; 0 <= nmixing <= 1
                ; 0 = Gaussian
                ; 1 = Lorentzian
                3:    begin
                    if global then begin
                        constradd.limited=1  ;TODO: only on n0 (now also on X)
                        constradd.limits[1]=1
                    endif else begin
                        constradd.limited=1
                        constradd.limits[1,*]=1
                    endelse
                    endcase
    
                ; 1 <= R
                ; 1 = Cauchy
                ; 2 = mod1 Lorentzian
                ; inf = Gaussian
                4:    begin
                    if global then begin
                        constradd.limited[0,*]=1 ;>=0
                    endif else  begin
                        constradd.limited[0,*]=1
                        constradd.limits[0,*]=1
                    endelse
                    endcase
                
                ; Asym > 0
                5:    begin
                    if global then begin
                        constradd.limited[0,*]=1 ;>=0
                    endif else begin
                        constradd.limited[0,*]=1
                        constradd.limits[0,*]=pres
                    endelse
                    endcase
                else:
                endcase
            endelse

            constr=[constr,constradd]
        endif;~kset
    
        A=[A,Aadd]
        nout+=n[0]
    endif
    
    ; Add to Afix
    if (n[1] ne 0) and nkset2 then begin
        tmp=transpose((*p)[ARRind1,*])
        Afix=[Afix,reform(tmp,n[1])]
    endif
    
    ; Add peakparam used in global to afix and/or A
    if global and (*ptr).codeproc[ind[i],6] ne 0 then begin
        if kset then p=std?(*(*ptr).peakparamIO.SD)[icurrent]:(*(*ptr).peakparamIO.O)[icurrent] $
        else p=(*ptr).peakparamIO.I
            
        ARRind=*(*ptr).paramind[ind[i],0] ; column indices in (*p)
        nARRind=n_elements(ARRind) ; n columns
            
        if bAref and (*ptr).codeproc[ind[i],6] eq 2 then begin
            ARRind1=where((*ptr).peakhelpparam[ARRind],ct2,comp=ARRind0,ncomp=ct3) ;2>fixed, 3>refined
            n=[ct3,ct2]*peaknr ;[nrefined,nfixed]
            if n[0] ne 0 then ARRind0=ARRind[ARRind0]
            if n[1] ne 0 then ARRind1=ARRind[ARRind1]
        endif else begin
            n=[0,nARRind*peaknr] ;[nrefined,nfixed]
            ARRind0=[-1]
            ARRind1=ARRind
        endelse

        ; Add to A: these are refined parameters but not in the least squares sense
        if n[0] ne 0 then begin
            Aadd=(*p)[ARRind0,*]
            
            ; Stabilize fit parameters
            ; Scaling factor was stabilized above and LeBail:
            if global and ind[i] eq 1  and ~batchexport then Aadd/=10.^(*ptr).ScalingMultiplier ; Peak intensities

            ; Constraints
            if (~kset) then begin
                constradd=MakeAInfo(n[0])
                constradd.nonls=1
                
                ; Le Bail: fix scaling factor => it's faster not to fix it
                ;if (*ptr).type eq 11 then if ExtractbitCode(ptr,2) eq 1 then constr[n_elements(A)-1-n[0]].fixed=1
                constr=[constr,constradd]
            endif
            
            A=[A,reform(transpose(Aadd),n[0])]
            nout+=n[0]
        endif
            
        ; Add to Afix
        if (n[1] ne 0) and nkset2 then begin
            tmp=(*p)[ARRind1,*]
            
            ; Stabilize fit parameters
            ; Scaling factor was stabilized above and no LeBail:
            if global and ind[i] eq 1 and n[0] eq 0  and ~batchexport then tmp[n_elements(ARRind1)-1,*]/=10.^(*ptr).ScalingMultiplier ; Peak intensities
            
            Afix=[Afix,reform(transpose(tmp),n[1])]
        endif
    endif
    
    ; Add ASU parameters (no Afix needed because pointer to initial ASU is passed)
    if (*ptr).type eq 12 and (ind[i] eq 1) and bAref then begin
    
        if kset then p=std?(*(*ptr).ASU.SD)[icurrent]:(*(*ptr).ASU.O)[icurrent] $
        else p=(*ptr).ASU.I
    
        ; xyzfit
        ind2=(*(*ptr).preverse_indices)[0:(*ptr).nasupos-1]
        xyz=transpose((*p)[ind2].wyck.xyz)
        ind2=where(xyz ne 999 and ~transpose((*p)[ind2].wyck.togglefix),n)
        if n ne 0 then begin
            A=[A,reform(xyz[ind2])]
            if (~kset) then constr=[constr,MakeAInfo(n)]
            nout+=n
        endif
    
        ; SOF
        ind2=where(~(*p).togglefix[0,*],n)
        if n ne 0 then begin
            A=[A,reform(((*p).SOF)[ind2])]
            
            if (~kset) then begin
                    
                comboID=lonarr((*ptr).nasuatoms)
                comboULim=make_array(1,(*ptr).nasuatoms,value=1.)
                toggle=(*p).togglefix[0,*]
                SOF=toggle*((*p).SOF)
                ind3=(*(*ptr).preverse_indices)
                for j=0l,(*ptr).nasupos-1 do begin
                    ; If for this asu position more than 1 SOF refined: add comboID
                    ; Subtract the sum of the fixed SOFs from the upper limits for each asu position
                    i0=ind3[j]
                    i1=ind3[j+1]-1
                    if long(total(~toggle[i0:i1])) gt 1 then comboID[i0:i1]=max(comboID)+1
                    comboULim[i0:i1]-=total(SOF[i0:i1],/pres)
                endfor
                ind3=where(comboID ne 0,tmp)
                if tmp ne 0 then comboID[ind3]+=max(constr.comboID)
                comboID=comboID[ind2]
                comboULim=comboULim[ind2]
                    
                constradd=MakeAInfo(n)
                ; 0 <= SOF <= (1-fixed SOFs)
                constradd.limited=1
                for j=0l,n-1 do constradd[j].limits[1]=comboULim[j]
                    
                ; 0 <= SUM(SOF) <= (1-fixed SOFs)
                constradd.limitedcombo=1
                constradd.comboID=comboID
                constradd.combocoeff=1
                for j=0l,n-1 do constradd[j].limitscombo[1]=comboULim[j]
                    
                constr=[constr,constradd]
            endif
            nout+=n
        endif
    
        ; Biso
        ind2=where(~(*p).togglefix[1,*],n)
        if n ne 0 then begin
            A=[A,reform(((*p).B)[ind2])]
            
            if (~kset) then begin
                constradd=MakeAInfo(n)
                ; 0 <= Biso
                constradd.limited[0,*]=1
                constr=[constr,constradd]
            endif
            nout+=n
        endif

    endif ;Add ASU parameters
    
    ; Add properties derived from scaling factor
    if ~nkset2 and ~std and (ind[i] eq 1) and bAref and bPhysicalPhaseInfo then begin
        tmp=*(*(*ptr).propO)[icurrent]
        A=[A,tmp]
        nout+=n_elements(tmp)
    endif

    endif ;nkset2 or bAref
endfor

; Add statistics
if ~nkset2 and ptr_valid((*ptr).factors) and ~std then begin
    tmp=reform(*(*(*ptr).factors)[icurrent])
    A=[A,tmp]
    nout+=n_elements(tmp)
endif

end;pro MakeAFromTypex
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro typexincrease,ptr
if (*ptr).peaknr eq 0 then return
p=(*ptr).globalparamIO.O
(*p)=[(*p),ptr_new((*(*p)[0])*0)]
p=(*ptr).globalparamIO.SD
(*p)=[(*p),ptr_new((*(*p)[0])*0)]
p=(*ptr).peakparamIO.O
(*p)=[(*p),ptr_new((*(*p)[0])*0)]
p=(*ptr).peakparamIO.SD
(*p)=[(*p),ptr_new((*(*p)[0])*0)]
p=(*ptr).propO
(*p)=[(*p),ptr_new((*(*p)[0])*0)]

if (*ptr).type eq 12 then begin
    if (*ptr).nasupos ne 0 then begin
        p=(*ptr).asu.O
        (*p)=[(*p),ptr_new(*(*p)[0])]
        p=(*ptr).asu.SD
        (*p)=[(*p),ptr_new(*(*p)[0])]
    endif
endif

p=(*ptr).factors
if ptr_valid(p) then (*p)=[(*p),ptr_new((*(*p)[0])*0)]
end;pro typexincrease
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ptrfitincrease,ptrFit
p=(*ptrFit).factors
(*p)=[(*p),ptr_new((*(*p)[0])*0)]
end;pro ptrfitincrease
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MakeTypexFromA,ptr,A,Afix,sigma,off,indAReject,icurrent,typeinfo,refit,interactive=interactive

off0=off
ind=where((*ptr).codeproc[*,0],ct)
peaknr=(*ptr).peaknr
if peaknr eq 0 then return
if (*ptr).type eq 12 then $
if (*ptr).nasupos eq 0 then return

profiletype=ExtractbitCode(ptr,0)

ireject=-1
bglobalreject=0b

; Copy not-used params
(*(*(*ptr).globalparamIO.O)[icurrent])[*]=0
(*(*(*ptr).globalparamIO.SD)[icurrent])[*]=0
(*(*(*ptr).peakparamIO.O)[icurrent])=(*(*ptr).peakparamIO.I)
(*(*(*ptr).peakparamIO.SD)[icurrent])[*]=0

; Extract used params
for i=0l,ct-1 do begin

    ; Point to global or peakparam
    global=(*ptr).codeproc[ind[i],1]
    bAref=(*ptr).codeproc[ind[i],2]
    
    p3=(*(*ptr).peakparamIO.O)[icurrent]
    p4=(*(*ptr).peakparamIO.SD)[icurrent]
    if global then begin
        p0=(*(*ptr).globalparamIO.O)[icurrent]
        p1=(*(*ptr).globalparamIO.SD)[icurrent]
    endif else begin
        p0=p3
        p1=p4
    endelse
    
    ; Extract from A or/and Afix
    ARRind=*(*ptr).paramind[ind[i],global] ; column indices in (*p)
    nARRind=n_elements(ARRind) ; n columns
    if global then begin
        if ~bAref then n=[0,nARRind] else $
        n=histogram((*ptr).togglefix[ARRind],min=0,max=1,REVERSE_INDICES=rtoggle) ;[nrefined,nfixed]
        if n[0] eq 0 then begin
            ARRind0=[-1]
            ARRind1=ARRind
        endif else begin
            ARRind0=ARRind[rtoggle[rtoggle[0]:rtoggle[1]-1]]
            if n[1] ne 0 then ARRind1=ARRind[rtoggle[rtoggle[1]:rtoggle[2]-1]] else ARRind1=[-1]
        endelse
    endif else begin
        if bAref then begin
            ARRind1=where((*ptr).peakhelpparam[ARRind],ct2,comp=ARRind0,ncomp=ct3) ;2>fixed, 3>refined
            n=[ct3,ct2]*peaknr ;[nrefined,nfixed]
            if n[0] ne 0 then ARRind0=ARRind[ARRind0]
            if n[1] ne 0 then ARRind1=ARRind[ARRind1]
        endif else begin
            n=[0,nARRind*peaknr] ;[nrefined,nfixed]
            ARRind0=[-1]
            ARRind1=ARRind
        endelse
    endelse
    
    ; Extract from A
    if n[0] ne 0 then begin
        tmp0=A[off[0]:off[0]+n[0]-1]
        tmp1=sigma[off[0]:off[0]+n[0]-1]

        if ~global then begin
            tmp0=reform(tmp0,peaknr,n[0]/peaknr)
            tmp0=transpose(tmp0)
            tmp1=reform(tmp1,peaknr,n[0]/peaknr)
            tmp1=transpose(tmp1)
        endif
        
        ; Destabilize fit parameters
        if global and ind[i] eq 1 then begin
            tmp0/=10.^(*ptr).ScalingMultiplier ; Scaling factor
            tmp1/=10.^(*ptr).ScalingMultiplier
        endif
    
        (*p0)[ARRind0,*]=tmp0
        (*p1)[ARRind0,*]=tmp1
    
        tmp=where((indAReject ge off[0]) and (indAReject lt off[0]+n[0]),ctrej)
        if ctrej ne 0 then begin
            if global then bglobalreject=1b $
            else ireject=[ireject,(indAReject[tmp]-off[0]) mod n[0]] ; peaks to reject
        endif
    endif
    
    ; Extract from Afix
    if n[1] ne 0 then begin
        tmp0=Afix[off[1]:off[1]+n[1]-1]
        if ~global then begin
            tmp0=reform(tmp0,peaknr,n[1]/peaknr)
            tmp0=transpose(tmp0)
        endif
        (*p0)[ARRind1,*]=tmp0
        (*p1)[ARRind1,*]=0
    endif

    off+=n
        
    ; Connect celparam
    if (*ptr).type ne 10 and ind[i] eq 0 then begin
        tmp0=(*p0)[ARRind,*]
        (*p0)[ARRind,*]=[tmp0[0],celparamsetcon(tmp0[1:*],*(*ptr).pstrinfo)]
    endif
    
    ; Extract peakparam used in global from afix and/or A
    if global and (*ptr).codeproc[ind[i],6] ne 0 then begin
        ARRind=*(*ptr).paramind[ind[i],0] ; column indices in (*p)    
        nARRind=n_elements(ARRind) ; n columns

        if bAref and (*ptr).codeproc[ind[i],6] eq 2 then begin
            ARRind1=where((*ptr).peakhelpparam[ARRind],ct2,comp=ARRind0,ncomp=ct3) ;2>fixed, 3>refined
            n=[ct3,ct2]*peaknr ;[nrefined,nfixed]
            if n[0] ne 0 then ARRind0=ARRind[ARRind0]
            if n[1] ne 0 then ARRind1=ARRind[ARRind1]
        endif else begin
            n=[0,nARRind*peaknr] ;[nrefined,nfixed]
            ARRind0=[-1]
            ARRind1=ARRind
        endelse
            
        ; These are refined parameters but not in the least squares sense
        if n[0] ne 0 then begin
            tmp0=A[off[0]:off[0]+n[0]-1]
            tmp1=sigma[off[0]:off[0]+n[0]-1]
        
            tmp0=reform(tmp0,peaknr,n[0]/peaknr)
            tmp0=transpose(tmp0)
            tmp1=reform(tmp1,peaknr,n[0]/peaknr)
            tmp1=transpose(tmp1)
            
            ; Destabilize fit parameters
            if global and ind[i] eq 1 then begin
                ; Scaling factor was destabilized above and LeBail:
                tmp0*=10.^(*ptr).ScalingMultiplier ; Peak intensities
                tmp1*=10.^(*ptr).ScalingMultiplier
            endif

            (*p3)[ARRind0,*]=tmp0
            (*p4)[ARRind0,*]=tmp1
        endif
    
        if n[1] ne 0 then begin
            tmp0=Afix[off[1]:off[1]+n[1]-1]

            colnr=n[1]/peaknr
            tmp0=reform(tmp0,peaknr,colnr)
            tmp0=transpose(tmp0)
            
            ; Destabilize fit parameters
            if global and ind[i] eq 1 and n[0] eq 0 then begin
                ; Scaling factor was destabilized above and no LeBail:
                tmp0[colnr-1,*]*=10.^(*ptr).ScalingMultiplier ; Peak intensities
            endif

            (*p3)[ARRind1,*]=tmp0
            (*p4)[ARRind1,*]=0
        endif
            
        off+=n
    endif

    ; Extract ASU parameters
    if (*ptr).type eq 12 and ind[i] eq 1 then begin
        p=(*ptr).asu.I
        p0=(*(*ptr).asu.O)[icurrent]
        p1=(*(*ptr).asu.SD)[icurrent]
        *p0=*p
        *p1=*p
        
        (*p1).xyz=0
        (*p1).SOF=0
        (*p1).B=0
        (*p1).wyck.xyz*=(*p1).wyck.xyz eq 999
    
        if bAref then begin
        ; xyz
        ind_=(*(*ptr).preverse_indices)[0:(*ptr).nasupos-1]
        xyz=transpose((*p)[ind_].wyck.xyz)
        induse=where(xyz ne 999 and ~transpose((*p)[ind_].wyck.togglefix),ctuse)
        if ctuse ne 0 then begin
            xyzsd=xyz
            xyz[induse]=A[off[0]:off[0]+ctuse-1]
            xyzsd[induse]=sigma[off[0]:off[0]+ctuse-1]
            off[0]+=ctuse
            xyzsd=transpose(xyzsd)
            xyz=transpose(xyz)
    
            ind2=(*(*ptr).preverse_indices)
            for j=0l,(*ptr).nasupos-1 do begin
                RT=Symop64((*p)[ind_[j]].wyck.Rep)
                xyzj=reform(RT.rot##transpose(xyz[*,j]))+RT.trn
                ZOCompress0,xyzj
                xyzsdj=sqrt(reform((RT.rot*RT.rot)##transpose(xyzsd[*,j]*xyzsd[*,j])))
    
                nsub=ind2[j+1]-ind2[j]
                (*p0)[ind2[j]:ind2[j+1]-1].xyz=rebin(xyzj,3,nsub,/sample)
                (*p0)[ind2[j]:ind2[j+1]-1].wyck.xyz=rebin(xyz[*,j],3,nsub,/sample)
                (*p1)[ind2[j]:ind2[j+1]-1].xyz=rebin(xyzsdj,3,nsub,/sample)
                (*p1)[ind2[j]:ind2[j+1]-1].wyck.xyz=rebin(xyzsd[*,j],3,nsub,/sample)
            endfor
        endif
    
        ; SOF
        induse=where(~(*p).togglefix[0,*],ctuse)
        if ctuse ne 0 then begin
            (*p0)[induse].SOF=A[off[0]:off[0]+ctuse-1]
            (*p1)[induse].SOF=sigma[off[0]:off[0]+ctuse-1]
            off[0]+=ctuse
        endif
    
        ; Biso
        induse=where(~(*p).togglefix[1,*],ctuse)
        if ctuse ne 0 then begin
            (*p0)[induse].B=A[off[0]:off[0]+ctuse-1]
            (*p1)[induse].B=sigma[off[0]:off[0]+ctuse-1]
            off[0]+=ctuse
        endif
    
        endif ;bAref
    endif
endfor

; ----Reject peaks/phases and ask for refit----
; ireject: indices of the peaks that cause rejection
; bglobalreject: some global parameter caused rejection
bphasereject=0b
if (*ptr).codeproc[1,1] then begin
    ; Global peak intensity used
    if bglobalreject or n_elements(ireject) ne 1 then begin
        ; Reject phase
        ind=(*(*ptr).paramind[1,1]) ; Scaling factor
        
        ; Set refined scaling factor to zero
        p0=(*(*ptr).globalparamIO.O)[icurrent]
        p1=(*(*ptr).globalparamIO.SD)[icurrent]
        (*p0)[ind]=0
        (*p1)[ind]=0
        
        ; Add phase to refit info
        *refit.pkeepinit=[*refit.pkeepinit,ptr_new()]
        *refit.presetinclude=[*refit.presetinclude,ptr]
        
        ; Exclude phase and ask for refit
        (*ptr).include=0
        bphasereject=1b
        refit.brefit or=1b
        
        printw,interactive,'Rejecting '+(*ptr).name+'...'
    endif
endif else begin
    ; Single peak intensities used
    if n_elements(ireject) ne 1 then begin
        ; Reject single peaks
        ind=(*(*ptr).paramind[1,0]) ; Intensity parameters index
        ind=ind[n_elements(ind)-1] ; get the last one
        
        ; Set peak intensities to zero
        (*p3)[ind,ireject[1:*]]=0
        (*p4)[ind,ireject[1:*]]=0
        
        ; Add phase to refit info (if not already there)
        tmp=where(*refit.presetinclude eq ptr,ct)
        if ct eq 0 then begin
            *refit.pkeepinit=[*refit.pkeepinit,ptr_new(*(*ptr).peakparamIO.I)]
            *refit.presetinclude=[*refit.presetinclude,ptr]
        endif
            
        ; Set initial peaks to zero and exclude phase when all peaks are zero
        if ~array_equal((*(*ptr).peakparamIO.I)[ind,ireject[1:*]],0) then begin
            ; Not all "ireject"-peak intensities already 0: set to 0
            (*(*ptr).peakparamIO.I)[ind,ireject[1:*]]=0
            
            ; Exclude phase when all zero
            tmp=where((*(*ptr).peakparamIO.I)[ind,*] eq 0,ct)
            if ct eq (*ptr).peaknr then begin
                (*ptr).include=0
                bphasereject=1b
            endif
            
            ; Ask for refit
            refit.brefit or=1b
        endif
    endif
endelse

if bphasereject then begin
    ; When phase rejected, set R-factors to 0
    if ptr_valid((*ptr).factors) then begin
        p=(*(*ptr).factors)[icurrent]
        (*p)[*]=0
    endif
    
    ; Set refined parameters to zero:
    ;MakeTypexFromA,ptr,A*0,Afix,sigma*0,off0,[-1],icurrent,typeinfo,refit
endif

end;pro MakeTypexFromA
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FitInfoLabels,all=all
labels=['npts', 'nParamRef', 'nParamCon', 'nPegged', 'nReject', 'nParamFix','iter', 'DOF', 'Time(sec)', 'CHI2', 'CHI2red', 'Rp(%)', 'Rwp(%)', 'Re(%)', 'cRp(%)', 'cRwp(%)', 'cRe(%)', 'D', 'dw', 'Qd']
if keyword_set(all) then labels=[labels,'D+p.ln(p)','D+2p','Qd<d<4Qd']
return,labels
end;function FitInfoLabels
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FitModelLabels,ptrFit,Labels

Labels=''

ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if (*ptr).include then begin
        MakeLabelsFromTypex,ptr,Labels
        MakeLabelsFromTypex,ptr,Labels,/std
    endif
    ptr=(*ptr).next
endwhile

n=n_elements(Labels)
b=n ne 1
if b then Labels=[Labels[1:*],'- '+FitInfoLabels(/all)+' 1Dfit -']

return,b

end;function FitModelLabels
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InitFitModelShrink,ptr,icurrent,A,nparam

A=0.
nparam=0
if (*ptr).include then begin
    MakeAFromTypex,ptr,A,nparam1,icurrent,refinedA=2,/batchexport
    MakeAFromTypex,ptr,A,nparam,icurrent,refinedA=2,/std,/batchexport
    nparam+=nparam1
endif

b=nparam ne 0
if b then A=A[1:*]

return,b

end;function InitFitModelShrink
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InitFitModel,x,y,ptrModel,A,Afix,nAfix,constr,initinfo=initinfo,interactive=interactive,type=type

; Extract fit parameters from model
; Type can be:
;    0. start of fit (refinedA not set)
;    1. plotting initial model (refinedA not set)
;    2. plotting refined model (refinedA set)
if ~keyword_set(type) then type=0
case type of
0:     begin
    refinedA=0
    recall=1b
    endcase
1:     begin
    refinedA=0
    recall=1b
    endcase
2:     begin
    refinedA=1
    recall=0b
    endcase
endcase
nkset= refinedA eq 0

; ----Check and initialize----
ptrFit=(*ptrModel).ptrFit
nparamarr=total((*ptrFit).n)
nAfix=0
if nparamarr eq 0 then return,0
nparam=intarr(nparamarr+1)
iparam=0
nparam[iparam]=0
delvar2,A

; ----Background----
; recall -> recalculate background
code=ExtractbitCode(ptrFit,0)
case code of
0:    begin ;None -> set yb=0
    ptr_free,(*ptrFit).yb
    (*ptrFit).yb=ptr_new(0)
    endcase
1:    begin ;Calc+substract (strip) -> set yb=stripped y
    ptr_free,(*ptrFit).yb
    yb=Strip1D(y,(*ptrFit).backparamIO[2],(*ptrFit).backparamIO[1])
    if (*ptrFit).fitcontrol[4] then y-=yb
    (*ptrFit).yb=ptr_new(yb)
    endcase
2:    begin ;Calc+substract (ortpol) -> set yb=ortpol
    ptr_free,(*ptrFit).yb
    tmp=(*ptrFit).backparamIO[0]
    if tmp ge 0 then degree=tmp
    yb=ortpol(x,y,degree=degree,GR=GR,print=interactive)
    (*ptrFit).backparamIO[0]=GR
    if (*ptrFit).fitcontrol[4] then y-=yb
    (*ptrFit).yb=ptr_new(yb)
    endcase
3:    begin ;Fit (ortpol) -> set yb=0 and POL[0]=POL,POL[1]=COEFF (POL[1] is reset after fit)
    (*ptrFit).fitcontrol[4]=0 ;"backsub before fit" makes no sense
    ptr_free,(*ptrFit).yb
    (*ptrFit).yb=ptr_new(0)

    if recall then begin
        tmp=(*ptrFit).backparamIO[0]
        if tmp ge 0 then degree=tmp
        *(*ptrFit).yb=ortpol(x,y,COEFF=COEFF,POL=POL,degree=degree,GR=GR,print=interactive)
        ;if type eq 0 then begin
            (*ptrFit).backparamIO[0]=GR
            PTR_FREE,(*ptrFit).POL[0]
            (*ptrFit).POL[0]=PTR_NEW(POL)
            PTR_FREE,(*ptrFit).POL[1]
            (*ptrFit).POL[1]=PTR_NEW(COEFF)
        ;endif
    endif else begin
        GR=(*ptrFit).backparamIO[0]
        COEFF=*(*ptrFit).POL[1]
    endelse

    if nkset then constr=MakeAInfo(GR+1)

    nparam[iparam]=GR+1
    A=COEFF
    
    endcase
endcase
iparam++

; ----Loop over all models----
ptr=(*ptrFit).next
bA=n_elements(A) eq 0
if bA then begin
    A=0.
    if nkset then begin
        constr=MakeAInfo(1)
        constr.fixed=1b
    endif
endif
Afix=0.

while PTR_VALID(ptr) do begin
    if (*ptr).include then begin
        if keyword_set(initinfo) then begin
            ;tmp=InitScaling(ptr,initinfo,x,y-(*(*ptrFit).yb),/set,icurrent=(*ptrFit).icurrent)
            initinfo.finalize or= (*ptr).type ne 10
        endif
        MakeAFromTypex,ptr,A,n,(*ptrFit).icurrent,constr,Afix,refinedA=refinedA
        nparam[iparam]+=n
        
        ; Prepare output
        if (*ptrFit).icurrent ne 0 then typexincrease,ptr
    endif
    iparam++
    ptr=(*ptr).next
endwhile

bret=total(nparam[1:*]) ne 0
if bA and bret then begin
    A=A[1:*]
    if nkset then constr=constr[1:*]
endif
if ~bret then delvar2,A
nAfix=n_elements(Afix)-1
if nAfix ne 0 then Afix=Afix[1:*]

; ----Finish----
if nkset then begin
    PTR_FREE,(*ptrFit).nparam
    (*ptrFit).nparam=PTR_NEW(nparam)
endif

return,bret

end;function InitFitModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FinFitModel,ptrModel,x,yfit,A,Afix,sigma,factors,indAReject,typeinfo,refit,interactive=interactive

; Put A, yfit and factors in model

ptrFit=(*ptrModel).ptrFit
icurrent=(*ptrFit).icurrent

; ----Background----
if (*ptrFit).fitcontrol[4] then yfit+=(*(*ptrFit).yb)

; ----Check and initialize----
(*(*(*ptrFit).yfit)[icurrent])[(*ptrFit).backranind[0]:(*ptrFit).backranind[1]]=yfit
factors=transpose(factors)
if icurrent ne 0 then ptrfitincrease,ptrFit
(*(*(*ptrFit).factors)[icurrent])=factors

off=intarr(2) ; A,Afix

; ----Background----
if (*(*ptrFit).nparam)[0] ne 0 then begin
    off[0]+=(*(*ptrFit).nparam)[0]
    PTR_FREE,(*ptrFit).POL[1]
    (*ptrFit).POL[1]=PTR_NEW(A[0:off[0]-1])
endif

; ----Loop over all models----
ptr=(*ptrFit).next
while PTR_VALID(ptr) do begin
    if (*ptr).include then begin
        typeinfo.type=(*ptr).type
        MakeTypexFromA,ptr,A,Afix,sigma,off,indAReject,(*ptrFit).icurrent,typeinfo,refit,interactive=interactive
    endif
    ptr=(*ptr).next
endwhile

end;pro FinFitModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FitModelReset,refit
n=n_elements(*refit.presetinclude)-1
for i=1,n do begin
    ptr=(*refit.presetinclude)[i]
    (*ptr).include=1
    
    p=(*refit.pkeepinit)[i]
    if ptr_valid(p) then begin
        ptr_free,(*ptr).peakparamIO.I
        (*ptr).peakparamIO.I=p
    endif
endfor

*refit.presetinclude=0
*refit.pkeepinit=0
heap_free,refit
end;pro FitModelReset
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FitModelDerived,ptrModel,typeinfo,refit

ptrFit=(*ptrModel).ptrFit

; Reset refit info
FitModelReset,refit

; Quantitative info
SetPhysicalPhaseInfo,ptrFit,typeinfo.lambda,/O
end;pro FitModelDerived
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FitModel,x,y,yerror,ptrModel,typeinfo,refreshev,interactive=interactive,error=error,correlation=correlation,debug=debug

error=1b
T=systime(1)

; ----Refitting loop----
refit={brefit:0b,presetinclude:ptr_new(ptr_new()),pkeepinit:ptr_new(ptr_new())}
itertot=0L
initialize=1
ykeep=y

repeat begin
    ; ----Init----
    refit.brefit=0b
    factors=fltarr(23)
    indAReject=[-1]
    nconstrRej=0
    
    y=ykeep
    b=InitFitModel(x,y,ptrModel,A,Afix,nAfix,constr,initinfo=typeinfo,interactive=interactive)
    
    struc={typeinfo:typeinfo,$
        y:y,$
        ptr:(*ptrModel).ptrFit,$
        refreshev:refreshev}
    if nAfix ne 0 then struc=create_struct(struc,'Afix',Afix)
    
    ; ----Nothing to refine----
    if n_elements(A) eq 0 then begin
        if keyword_set(interactive) then temp=dialog_message('Nothing to refine.',/info)
        FitModelDerived,ptrModel,struc.typeinfo,refit
        error=0b
        return
    endif
    
    ; ----Some parameter initialization (e.g. scaling factor)----
;    gfuncFitModel, x, A, tmp,struc,initialize=initialize
;    gfuncFitModel, x, A, tmp,_extra=struc,initialize=initialize
    gfuncFitModel, x, A, tmp,struc=struc,initialize=initialize
    
    ; Force parameters within limits
    ind=where((constr.limited)[0,*],ct)
    if ct ne 0 then A[ind]>=(constr.limits)[0,ind]
    ind=where((constr.limited)[1,*],ct)
    if ct ne 0 then A[ind]<=(constr.limits)[1,ind]
    
    ; ----Fitting----
    ; Based on IDL's curvefit
;        temp1=where(ySTDEV eq 0,ct,comp=temp2,ncomp=ct2)
;        if ct ne 0 then $
;            if ct2 eq 0 then ySTDEV+=1 $ ; All zero
;            else ySTDEV[temp1]=min(ySTDEV[temp2])
;    yfit=MarquardtC(x,y, 1./ySTDEV^2., A, sigma,FUNCTION_NAME='gfuncFitModel',$
;        CHISQ=CHISQ,ITER=ITER,itmax=(*(*ptrModel).ptrFit).fitcontrol[0],nfree=dof,struc=struc,$
;        tol=(*(*ptrModel).ptrFit).fitcontrol[1],CError=CError,numder=(*(*ptrModel).ptrFit).fitcontrol[3])
;    ErrorMsg=''
;    chisq*=dof
;    npegged=0
    
    ; MPfit
;        temp1=where(ySTDEV eq 0,ct,comp=temp2,ncomp=ct2)
;        if ct ne 0 then $
;            if ct2 eq 0 then ySTDEV+=1 $ ; All zero
;            else ySTDEV[temp1]=min(ySTDEV[temp2])
;    yfit= mpcurvefit(x,y, 1./ySTDEV^2., A, sigma,FUNCTION_NAME='gfuncFitModel',CHISQ=CHISQ,ITER=ITER,$;,parinfo=constr,$
;                    itmax=(*(*ptrModel).ptrFit).fitcontrol[0],dof=dof,status=cerror,errmsg=errmsg,/nocatch,$
;                    tol=(*(*ptrModel).ptrFit).fitcontrol[1],FUNCTARGS=struc,/quiet,noderivative=(*(*ptrModel).ptrFit).fitcontrol[3])
;    ErrorMsg=errmsg
;    npegged=0

    if (*(*ptrModel).ptrFit).fitcontrol[5] eq 1 then begin
      ySTDEV = yerror
      ind = where(ySTDEV eq 0,ct)
      if ct ne 0 then ySTDEV[ind] = 1
    endif
    
    ; Based on MPfit
    numder=(*(*ptrModel).ptrFit).fitcontrol[3]
    if keyword_set(debug) then numder=1
    yfit=NLLS(x, y, A, 'gfuncFitModel',sigma=sigma,struc=struc,$
            Ainfo=constr,iter=iter,npegged=npegged,CError=CError,ErrorMsg=ErrorMsg,chisq=chisq,$
            dof=dof,numder=numder,$
            tol=(*(*ptrModel).ptrFit).fitcontrol[1],ySTDEV=ySTDEV,$
            itmax=(*(*ptrModel).ptrFit).fitcontrol[0],weight=(*(*ptrModel).ptrFit).fitcontrol[5],$
            correlation=correlation,debug=debug)
    itertot+=iter
    
    ; ----Calculate R-factors----
    if struc.typeinfo.finalize then begin
        P=A
        gfuncFitModel, x, P, tmp,struc=struc,/finalize
    endif
    
    ; ----Handle fit errors----
    if CError lt 0 then begin
        if keyword_set(interactive) then begin
            temp=dialog_message(ErrorMsg,/error)
            FitModelDerived,ptrModel,struc.typeinfo,refit
            return
        endif
    endif else if CError eq 0 then if keyword_set(interactive) then temp=dialog_message(ErrorMsg)
    printw,interactive,ErrorMsg
    
    if n_elements(dof) eq 0 then dof=0
    if n_elements(ySTDEV) eq 0 then ySTDEV=1
    if n_elements(sigma) eq 0 then sigma=a*0
    
    ;----Set param. rejection criteria----
    tmp=where(constr.limited[0,*] ne 0 or constr.limited[1,*] ne 0,nconstr)
    if nconstr ne 0 then $
        indAReject=where(((constr.limited[0,*] eq 2) and (A lt constr.limits[0,*])) or $
                         ((constr.limited[1,*] eq 2) and (A gt constr.limits[1,*])),nconstrRej)
    
    ; ----factors----
    yresid=y-yfit
    
    factors[0]=n_elements(yfit) ;n
    factors[1]=factors[0]-dof ;p
    factors[2]=nconstr
    factors[3]=npegged
    factors[4]=nconstrRej
    factors[5]=nAfix+n_elements(A)-factors[1] ; The A's that are fixed
    factors[6]=itertot
    factors[7]=dof ;n-p
    
    factors[8]=systime(1)-T
    
    ; RELIABILITY FACTORS:
    factors[9]=chisq ;chi2
    factors[10]=chisq/dof ;chi2red
    factors[11]= 100*total(abs(yresid))/total(y) ;Rp(%)
    factors[13]= 100*sqrt(dof/total((y/ySTDEV)^2));Re(%)
    factors[12]=factors[13]*sqrt(factors[10]) ;Rwp(%)
    
    ;if ~chisqr_test(chisq,dof) then print,"Data can't be described by the model."
    
    if (*(*ptrModel).ptrFit).fitcontrol[4] then begin
        factors[14]=factors[11]
        factors[16]=factors[13]
        factors[15]=factors[12]
    endif else begin
        tmp=abs(y-(*(*(*ptrModel).ptrFit).yb))
        factors[14]= 100*total(abs(yresid))/total(tmp);cRp(%)
        factors[16]= 100*sqrt(dof/total((tmp/ySTDEV)^2));cRe(%)
        factors[15]=factors[16]*sqrt(factors[10]);cRwp(%)
    endelse
    
    tmp=alog(y/yfit)
    ind=where(FINITE(tmp) ne 1,ct)
    if ct ne 0 then tmp[ind]=0 ; What should be done here?
    factors[17]= 2*total(y*tmp-yresid);D
    
    keep=yresid/ySTDEV^2.
    factors[18]=total((keep[1:*]-keep[0:factors[0]-2])^2) ;dw
    factors[18]/=total(keep^2);dw
    factors[19]= 2*((factors[0]-1.)/dof-3.0901/sqrt(factors[0]+2));Qd
    
    factors[20]=factors[17]+factors[1]*alog(factors[1]);D+p lnp
    factors[21]=factors[17]+2*factors[1];D+2p
    factors[22]= (factors[18] gt factors[19]) + 2*(factors[18] lt 4-factors[19]);Qd<d<4-Qd
    
    ; ----Finish----
    FinFitModel,ptrModel,x,yfit,A,Afix,sigma,factors,indAReject,struc.typeinfo,refit,interactive=interactive
    
    if itertot eq 1 and initialize lt 0 then begin
        ; Stopped after first step and initial scaling calculated: try initial scaling 10x lower
        ;initialize=0.1
        ;refit.brefit=1b
    endif
    ;initialize=0
endrep until ~refit.brefit

FitModelDerived,ptrModel,struc.typeinfo,refit
if (*(*ptrModel).ptrFit).fitcontrol[2] ne 0 then widget_control,refreshev.top,send_event=refreshev
error=0b

end;pro FitModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CreateFitModelInfo,list,noinitsource=noinitsource,peakparam=peakparam
ptrFit=(*list.ptrModel).ptrFit
if ~keyword_set(noinitsource) then (*(*ptrFit).sourceemission)[0]=list.lambda
if keyword_set(peakparam) then pPeakParam=ptr_new(/alloc) else pPeakParam=ptr_new()
return,{lambda:list.lambda,$
        ltable:(*(*ptrFit).sourceemission),$
        itable:0,$
        maxtable:n_elements((*(*ptrFit).sourceemission))/3-1,$
        pol_dmono:*(*ptrFit).misc[5],$
        pol_Plindeg:*(*ptrFit).misc[6],$
        pol_n:*(*ptrFit).misc[8],$
        pol_mmult:*(*ptrFit).misc[9],$
        pol_s:*(*ptrFit).misc[10],$
        pol_azimuth:*(*ptrFit).misc[11],$
        LGPtype:ExtractbitCode(ptrFit,1),$
        zerotype:ExtractbitCode(ptrFit,3),$
        type:0,$
        finalize:0b,$
        pPeakParam:pPeakParam,$
        dist:list.dist,anomal:*(*ptrFit).misc[4],$
        SmallDev:*(*ptrFit).misc[3],pstrinfo:PTR_NEW(),pasuinfo:PTR_NEW()}
end;function CreateFitModelInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FitModelCHI,top,correlation=correlation,debug=debug

widget_control,top,get_uvalue=list

; Prepare Fit
x=(*list.xvalt)[list.backranind[0]:list.backranind[1],1] ;two-theta
y=(*list.spet)[list.backranind[0]:list.backranind[1],1]
yerror=(*list.speerrort)[list.backranind[0]:list.backranind[1],1]

ptrFit=(*list.ptrModel).ptrFit
p=(*ptrFit).yfit
ptr_free,(*p)[(*ptrFit).icurrent]
(*p)[(*ptrFit).icurrent]=ptr_new((*list.spet)[*,list.xtype]*0)
(*ptrFit).backranind=list.backranind

FitModelStabilize,list,x,y

; Fit
;y=double(y)
FitModel,x,y,yerror,list.ptrModel,CreateFitModelInfo(list),$
    {ID:top,TOP:top,HANDLER:top,ACTION:'refresh'},interactive=list.OutID,correlation=correlation,debug=debug

; Finish fit
FitModelDestabilize,list

end;pro FitModelCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FitModelBatch,x,y,yerror,list,dummy=dummy,error=error

ptrFit=(*list.ptrModel).ptrFit
icurrent=(*ptrFit).icurrent
n=n_elements(y)
if icurrent eq 0 then begin
    (*(*ptrFit).yfit)[icurrent]=ptr_new(fltarr(n))
endif else begin
    p=(*ptrFit).yfit
    (*p)=[(*p),ptr_new(fltarr(n))]
endelse

if keyword_set(dummy) then begin
    if icurrent ne 0 then begin
        ptr=(*ptrFit).next
        while PTR_VALID(ptr) do begin
            if (*ptr).include then typexincrease,ptr
            ptr=(*ptr).next
        endwhile
        ptrfitincrease,ptrFit
    endif
endif else begin
    ; Prepare Fit
    (*ptrFit).backranind=[0,n-1]
    (*ptrFit).fitcontrol[2]=0
    if icurrent eq 0 then FitModelStabilize,list,x,y

    ; Fit
    FitModel,x,y,yerror,list.ptrModel,CreateFitModelInfo(list),{top:0L},error=error
;    if error then return
endelse

; Plot
if list.debug then begin
    if winidvalid(list.outid2) then wset,list.outid2 $
    else begin
        window,winid(),title='Batch fitting:'
        list.outid2=!d.window
    endelse
    if keyword_set(dummy) then begin
        erase
    endif else begin
        LoadctRev,-39,/silent
        if (*ptrFit).fitcontrol[4] then begin
            plot,x,y,xstyle=1,ystyle=1,psym=1,symsize=0.5
            oplot,x,(*(*(*ptrFit).yfit)[icurrent])-(*(*ptrFit).yb),color=10
        endif else begin
            plot,x,y,xstyle=1,ystyle=1,psym=1,symsize=0.5
            oplot,x,*(*(*ptrFit).yfit)[icurrent],color=10
        endelse
    endelse
endif

(*ptrFit).icurrent++
end;pro FitModelBatch
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SavefitmodelFitResultBegin_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_uvalue=uval

case uval.type of
0:    BEGIN
    widget_control,ev.top,get_uvalue=info
    info.sbutton[uval.i]=ev.select
    widget_control,ev.top,set_uvalue=info

    nSEsec=n_elements(info.SE)/3
    indselect=where(info.sbutton[0:nSEsec-1],ct)
    if ct eq 0 then begin
        info.sbutton[uval.i]=1b
        widget_control,ev.id,/set_button
        return
    endif
    
    *uval.ptr=(info.SE)[*,indselect]
    endcase
1:    BEGIN
    (*uval.ptr).include=ev.select
    endcase
3:    begin
    (*uval.ptr)=ev.select
    endcase
else: widget_control,ev.top,/destroy
endcase

end;pro SavefitmodelFitResultBegin_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SavefitmodelFitResultBegin,ptrFit,top,info,background

ioff=0L

; Select source emission lines
ptrlist=(*ptrFit).sourceemission
SE=*ptrlist
nlines=n_elements(SE)/3
sbutton=replicate(1b,nlines)
lbutton='Primary source emission line '+stringr(reform(SE[0,0]))
if nlines gt 1 then lbutton=[lbutton,'Secondary source emission line '+stringr(reform(SE[0,1:*]))]
ubutton=replicate({type:0,ptr:ptrlist,i:ioff++},nlines)
ubutton.i=lindgen(nlines)

; Select background
sbutton_=0
brackground_ptr=ptr_new(sbutton_)
sbutton=[sbutton,sbutton_]
lbutton=[lbutton,'Background']
ubutton=[ubutton,{type:3,ptr:brackground_ptr,i:ioff++}]

; Select group
ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if n_elements(include) eq 0 then include=(*ptr).include $
    else include=[include,(*ptr).include]
    
    sbutton=[sbutton,(*ptr).include]
    lbutton=[lbutton,(*ptr).name]
    ubutton=[ubutton,{type:1,ptr:ptr,i:ioff++}]

    ptr=(*ptr).next
endwhile

; User selection
ret=n_elements(include) ne 0
if ret then begin
    info={SE:SE,include:include,sbutton:sbutton}
    
    ; Make base
    base=widget_base(title='Select ',/column,uvalue=info,/modal,GROUP_LEADER=top)
    button=widget_button(base,value='OK',uvalue={type:2})
    baset=widget_base(base,/column,/nonexclusive)
    bbutton=lonarr(n_elements(sbutton))
    for i=0l,n_elements(sbutton)-1 do $
        bbutton[i]=widget_button(baset,value=lbutton[i],uvalue=ubutton[i])
    
    
    WIDGET_CONTROL,base, /REALIZE
    for i=0l,n_elements(bbutton)-1 do $
        widget_control,bbutton[i],set_button=sbutton[i]
    Xmanager,'SavefitmodelFitResultBegin',base, event_handler='SavefitmodelFitResultBegin_Event',$
        GROUP_LEADER=top
        
    background=*brackground_ptr
    ptr_free,brackground_ptr
endif

return,ret
end;function SavefitmodelFitResultBegin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SavefitmodelFitResultEnd,ptrFit,info
*(*ptrFit).sourceemission=info.SE

ptr=(*ptrFit).next
i=0
while ptr_valid(ptr) do begin
    (*ptr).include=info.include[i++]
    ptr=(*ptr).next
endwhile

end;pro SavefitmodelFitResultEnd
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakeYforFitModel,list,x,y
; x and y in two-theta
; return y in xtype

if list.xtype eq 1 then return,y

tmp=CHI_xval(x,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,1,spe=y,tspe=y2)
return,y2[*,list.xtype]
    
end;function MakeYforFitModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SavefitmodelFitResult,list,top,O=O
ptrFit=(*list.ptrModel).ptrFit
if ~SavefitmodelFitResultBegin(ptrFit,top,info,background) then return

CATCH, Error_status
IF Error_status NE 0 THEN begin
    SavefitmodelFitResultEnd,ptrFit,info
    return
endif

x=(*list.xvalt)[list.backranind[0]:list.backranind[1],1] ;two-theta
y=(*list.spet)[list.backranind[0]:list.backranind[1],1]
b=InitFitModel(x,y,list.ptrModel,A,Afix,nAfix,interactive=interactive,type=1+keyword_set(O))

if b then begin
    struc={ptr:ptrFit,typeinfo:CreateFitModelInfo(list,/noinitsource),y:y,refreshev:{top:0L}}
    if nAfix ne 0 then struc=create_struct(struc,'Afix',Afix)

    if ~background then begin
        ; Prevent background addition
        nbackparam=(*(*ptrFit).nparam)[0]
        fitcontrol=(*ptrFit).fitcontrol[4]
        (*ptrFit).fitcontrol[4]=1
        (*(*ptrFit).nparam)[0]=0
        if nbackparam ne 0 then A=A[nbackparam:*]
    endif
    
    gfuncFitModel, x, A, F, struc=struc
    
    if ~background then begin
        (*(*ptrFit).nparam)[0]=nbackparam
        (*ptrFit).fitcontrol[4]=fitcontrol
    endif else if (*ptrFit).fitcontrol[4] then F+=(*(*ptrFit).yb)

    y2=MakeYforFitModel(list,x,F)
;    y2/=100
    x2=(*list.xvalt)[list.backranind[0]:list.backranind[1],list.xtype]
    
    path=SelFile(list.path,list.file,ChiFormat(),'Save Fit as Profile...')
    if path eq '' then return
    
    temp=WriteChi(path,transpose([[x2],[y2]]),list.xtype,'Empty',list.OutID,list.ROIAzimuth)
endif

SavefitmodelFitResultEnd,ptrFit,info

end;pro SavefitmodelFitResult
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gfuncPeakparamBegin,ptrFit,ptrI,info

; Select first source emission line
SE=*(*ptrFit).sourceemission
*(*ptrFit).sourceemission=(*(*ptrFit).sourceemission)[*,0]

; Select group
ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if n_elements(include) eq 0 then include=(*ptr).include $
    else include=[include,(*ptr).include]
    
    (*ptr).include=ptrI eq ptr
    ptr=(*ptr).next
endwhile

; Keep info
ret=n_elements(include) ne 0
if ret then info={SE:SE,include:include}

return,ret
end;function gfuncPeakparamBegin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro gfuncPeakparamEnd,ptrFit,info
*(*ptrFit).sourceemission=info.SE

ptr=(*ptrFit).next
i=0
while ptr_valid(ptr) do begin
    (*ptr).include=info.include[i++]
    ptr=(*ptr).next
endwhile

end;pro gfuncPeakparamEnd
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gfuncPeakparam,list,ptr,O=O,d=d,tt=tt,x=x,y=y

ptrFit=(*list.ptrModel).ptrFit
if ~gfuncPeakparamBegin(ptrFit,ptr,info) then return,0

CATCH, Error_status
IF Error_status NE 0 THEN begin
    CATCH, /cancel
    if n_elements(struc) ne 0 then heap_free,struc.typeinfo.pPeakParam
    gfuncPeakparamEnd,ptrFit,info
    return,0
endif

if n_elements(x) eq 0 then xx=(*list.xvalt)[list.backranind[0]:list.backranind[1],1] else xx=x ;two-theta
if n_elements(y) eq 0 then yy=(*list.spet)[list.backranind[0]:list.backranind[1],1] else yy=y
b=InitFitModel(xx,yy,list.ptrModel,A,Afix,nAfix,interactive=interactive,type=1+keyword_set(O))

if b then begin
    struc={ptr:ptrFit,typeinfo:CreateFitModelInfo(list,/noinitsource,/peakparam),y:yy,refreshev:{top:0L}}
    if nAfix ne 0 then struc=create_struct(struc,'Afix',Afix)

    gfuncFitModel, xx, A, F, struc=struc
    
    ret=*struc.typeinfo.pPeakParam
    heap_free,struc.typeinfo.pPeakParam
endif

gfuncPeakparamEnd,ptrFit,info

if n_elements(ret) eq 0 then return,0

ttheta=ret[0,*]
ret=ret[1:*,*]

if keyword_set(d) then begin
    d=BraggTtoX(ttheta,struc.typeinfo.lambda,/onlyd,/angledeg)
    ret=[reform(d,1,n_elements(d)),ret]
endif
if keyword_set(tt) then begin
    ret=[reform(ttheta,1,n_elements(ttheta)),ret]
endif

return,ret
end;function gfuncPeakparam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ScalingFactorMultiplier,ptr,list,x,y,O=O,icurrent=icurrent

; Fitting the scaling factor is not always stable:
;    f(x) = S.mH.|FH|^2.LGP(x).G(x)
;            when mH.|FH|^2 very big, than S will be very small and v.v.
; This function returns a multiplier m so that
;    f(x) = (S.10^m).(mH/10^m).|FH|^2.LGP(x).G(x)
;            and (mH/10^m).|FH|^2.LGP(x).G(x) ~= 1

if n_elements(icurrent) eq 0 then icurrent=0
if keyword_set(O) then begin
    pg=(*(*ptr).globalparamIO.O)[icurrent]
endif else begin
    pg=(*ptr).globalparamIO.I
endelse

; Use global intensity?
if ~(*ptr).codeproc[1,1] then return,0

; Get peak intensities for scaling factor 1
S=(*pg)[*(*ptr).paramind[1,1]]
(*pg)[*(*ptr).paramind[1,1]]=1
ret=gfuncPeakparam(list,ptr,O=O,x=x,y=y)
(*pg)[*(*ptr).paramind[1,1]]=S

; Estimate multiplier
m=mean(abs(ret[0,*]))
if m eq 0 then return,0 else return,round(alog10(m)/2.)

end;function ScalingFactorMultiplier
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FitModelDestabilize,list,O=O,icurrent=icurrent

ptr=(*(*list.ptrModel).ptrFit).next
while PTR_VALID(ptr) do begin
    (*ptr).ScalingMultiplier=0
    ptr=(*ptr).next
endwhile

end;pro FitModelDestabilize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FitModelStabilize,list,x,y,O=O,icurrent=icurrent

ptr=(*(*list.ptrModel).ptrFit).next
while PTR_VALID(ptr) do begin
    (*ptr).ScalingMultiplier=0
    if (*ptr).include then (*ptr).ScalingMultiplier=ScalingFactorMultiplier(ptr,list,x,y,O=O,icurrent=icurrent)
    ptr=(*ptr).next
endwhile

end;pro FitModelStabilize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InitScaling,ptr,list,xtt,ytt,O=O,set=set,icurrent=icurrent

if keyword_set(set) and (*ptr).codeproc[1,1] ne 1 then return,0

if n_elements(icurrent) eq 0 then icurrent=0
if keyword_set(O) then begin
    p=(*(*ptr).peakparamIO.O)[icurrent]
    pg=(*(*ptr).globalparamIO.O)[icurrent]
endif else begin
    p=(*ptr).peakparamIO.I
    pg=(*ptr).globalparamIO.I
endelse

S=(*pg)[*(*ptr).paramind[1,1]]
(*pg)[*(*ptr).paramind[1,1]]=1
ret=gfuncPeakparam(list,ptr,O=O,/tt)
(*pg)[*(*ptr).paramind[1,1]]=S

; Estimate scaling
m=max(ret[1,*],indm)
mult=ret[2,*]*1.06447; FWHM*sqrt(!pi/alog(2))/2 (assume gaussian)
S=INTERPOL( ytt>0, xtt, ret[0,*])*mult/m
if n_elements(S) ne 1 then S=mean(S)+0.5*STDDEV(S)

if keyword_set(set) then (*pg)[*(*ptr).paramind[1,1]]=S
return,S

end;function InitScaling
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hkllabels,ptr,list,sd,O=O

icurrent=0
if keyword_set(O) then begin
    p=(*(*ptr).peakparamIO.O)[icurrent]
    pg=(*(*ptr).globalparamIO.O)[icurrent]
endif else begin
    p=(*ptr).peakparamIO.I
    pg=(*ptr).globalparamIO.I
endelse

labels=strcompress(string((*p)[*(*ptr).paramind[0,0],*],format='("(",I,I,I,")")'),/remove_all)

ret=gfuncPeakparam(list,ptr,O=O,/tt)
tt=ret[0,*]
;tt=ttfromhkl((*p)[*(*ptr).paramind[0,0],*],(*pg)[indCelParam(ptr)],list.lambda,$
;    list.dist,(*pg)[indZeroParam(ptr)],ExtractbitCode((*list.ptrModel).ptrFit,3))

n=n_elements(labels)
thick=make_array(n,value=1b)
if ptr eq list.highlightpeak.ptr then $
    if list.highlightpeak.i lt n and list.highlightpeak.i gt 0 then thick[list.highlightpeak.i]=2

for i=0l,n-1 do begin
    if labels[i] eq '' then continue
    ; Find neighbours within [-sd,+sd]
    ind=where(tt gt tt[i]-sd and tt lt tt[i]+sd,ct)
    if ct gt 1 then begin
        ; Add neighbours to current peaks
        ind=ind[where(ind ne i)]
        labels[i]+=strjoin(labels[ind],' ')
        labels[ind]=''
        thick[i]>=max(thick[ind])
    endif
endfor
ind=where(labels ne '')
labels=labels[ind]
tt=tt[ind]
thick=thick[ind]

return,{labels:labels,tt:tt,thick:thick}
end;function hkllabels
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitmodel2D,x,list,top,background

ptrFit=(*list.ptrModel).ptrFit
if ~SavefitmodelFitResultBegin(ptrFit,top,info,background) then return,x*0

CATCH, Error_status
IF Error_status NE 0 THEN begin
    SavefitmodelFitResultEnd,ptrFit,info
    return,x*0
endif

; No background
code=(*ptrFit).code
yb=*(*ptrFit).yb
nbackparam=(*(*ptrFit).nparam)[0]
SetbitCode,ptrFit,0,0
(*(*ptrFit).nparam)[0]=0
    
;x: two-theta
y=x*0
b=InitFitModel(x,y,list.ptrModel,A,Afix,nAfix,interactive=interactive,type=2)

if b then begin
    struc={ptr:ptrFit,typeinfo: CreateFitModelInfo(list),y:y,refreshev:{top:0L}}
    if nAfix ne 0 then struc=create_struct(struc,'Afix',Afix)

    gfuncFitModel, x, A, F, struc=struc
endif else F=y

SavefitmodelFitResultEnd,ptrFit,info

; Reset background
(*(*ptrFit).nparam)[0]=nbackparam
(*ptrFit).code=code
*(*ptrFit).yb=yb

return,F
end;function fitmodel2D
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro fitmodelplot,list,col,residual=residual

CATCH, Error_status
IF Error_status NE 0 THEN return

ptrFit=(*list.ptrModel).ptrFit

; ----Plot refined model----
p=(*(*ptrFit).yfit)[(*ptrFit).icurrent]
if ptr_valid(p) and (list.fitresolve eq 1) then begin
    np=n_elements(*p)-1
    b=0>list.backranind[0]<np
    e=0>list.backranind[1]<np
    if b eq e then return
    
    y=MakeYforFitModel(list,(*list.xvalt)[b:e,1],(*p)[b:e])
    powoplot,(*list.xvalt)[b:e,list.xtype],y,NOCLIP = 0,color=col[1],pow=list.ypower
endif

; ----Plot single peaks?----
if !d.name eq 'PS' then begin
    bsingle=dialog_message('Plot individual peaks?',/question) eq 'Yes'
endif else bsingle=1b

; ----Plot resolved peaks----
x_=(*list.xvalt)[list.backranind[0]:list.backranind[1],1] ;two-theta
y_=(*list.spet)[list.backranind[0]:list.backranind[1],1]
yback=(*list.yback)[list.backranind[0]:list.backranind[1]] ; initial background
x2=(*list.xvalt)[list.backranind[0]:list.backranind[1],list.xtype]

SepF=[1b,1b,0b]
residual=keyword_set(residual)
bool=[list.peaksflag and ~residual,(list.fitresolve eq 2) and ptr_valid(p) and ~residual,residual]
; initial, fitted, residual

for k=0l,2 do begin
if bool[k] ne 0 then begin
    ; Prepare Fit
    case k of
    0: type=1 ; initial parameters
    1: type=2 ; final parameters
    2: type=2 ; final parameters
    endcase
    
    x=x_
    y=y_
    b=InitFitModel(x,y,list.ptrModel,A,Afix,nAfix,interactive=interactive,type=type)
    
    if b then begin

        struc={ptr:ptrFit,typeinfo:CreateFitModelInfo(list),y:y,refreshev:{top:0L},highlightpeak:list.highlightpeak,thicki:-1l}
        if nAfix ne 0 then struc=create_struct(struc,'Afix',Afix)

        gfuncFitModel, x, A, F, struc=struc, SepF=SepF[k]

        ; Plot profile
        if SepF[k] then begin
            if (*ptrFit).fitcontrol[4] then $
                if type eq 1 then F[*,0]+=yback else F[*,0]+=(*(*ptrFit).yb)
            
            s=size(F)
            nrow=(s[0] eq 1)?1:s[2]
            F=reform(F,s[1],nrow,/overwrite)
    
            ; Plot background
            y2=MakeYforFitModel(list,x,F[*,0])
            powoplot,x2,y2,NOCLIP = 0,color=col[k]+1,thick=(struc.thicki eq 0)+1,linestyle=2,pow=list.ypower
                
            if ~bsingle then begin
                delvar2,coladd
                indoff=1
                
                ptr=(*ptrFit).next
                while ptr_valid(ptr) do begin
                    if (*ptr).include and (*ptr).peaknr ne 0 then begin
                        if n_elements(coladd) eq 0 then coladd=-col[k]+2 else coladd++
                        
                        y2=MakeYforFitModel(list,x,mtotal(F[*,indoff:indoff+(*ptr).peaknr-1],2)+F[*,0])
                        powoplot,x2,y2,NOCLIP = 0,color=col[k]+coladd,pow=list.ypower

                        indoff+=(*ptr).peaknr
                    endif
                    ptr=(*ptr).next
                endwhile
            
            endif else begin
    
                ; Plot all peaks
                for i=1,nrow-1 do begin
                    y2=MakeYforFitModel(list,x,F[*,i]+F[*,0])
                    powoplot,x2,y2,NOCLIP = 0,color=col[k],thick=(struc.thicki eq i)+1,linestyle=struc.thicki ne i,pow=list.ypower
                endfor

            endelse
            
            ; Plot total
            F=total(F,2)
            y2=MakeYforFitModel(list,x,F)
            powoplot,x2,y2,NOCLIP = 0,color=col[k],pow=list.ypower
            
        endif else begin
            if k eq 2 then begin
                ystdev=NLLSPoissonW(y,(*ptrFit).fitcontrol[5])
                F=(y-F)/ystdev

                y2=MakeYforFitModel(list,x,F)
                chisq=total(y2*y2) ; Chisq as a function of list.xtype
;                powplot,x2,y2,xtitle=list.xtitle,$
;                    ytitle='Weighted Residuals (a.u.) '+list.ypowerstr,psym=0,xrange=(*list.xvalt)[list.xrange,list.xtype],xstyle=1,ystyle=1,pow=list.ypower
;                powplots,[!X.CRange[0],!X.CRange[1]],[0,0],NOCLIP = 0,pow=list.ypower
;                ;powplots,[!X.CRange[0],!X.CRange[1]],[-3,-3],NOCLIP = 0,pow=list.ypower
;
;                powxyouts,0.05,0.05,'CHISQ: '+stringr(chisq),/normal,pow=list.ypower
                
                plot,x2,y2,xtitle=list.xtitle,$
                    ytitle='Weighted Residuals (a.u.) ',psym=0,xrange=(*list.xvalt)[list.xrange,list.xtype],xstyle=1,ystyle=1
                plots,[!X.CRange[0],!X.CRange[1]],[0,0],NOCLIP = 0

                xyouts,0.05,0.05,'CHISQ: '+stringr(chisq),/normal
            endif else begin
                if (*ptrFit).fitcontrol[4] then $
                    if type eq 1 then F[*,0]+=yback else F+=(*(*ptrFit).yb)
                y2=MakeYforFitModel(list,x,F)
                powoplot,x2,y2,NOCLIP = 0,color=col[k],pow=list.ypower
            endelse
        endelse

        ; Plot labels
        mi=list.ylog*(k ne 2)
        sd=abs((*list.xvalt)[1,1]-(*list.xvalt)[0,1])
        if list.sLabels and ~(k eq 2 and list.fitresolve ne 2) then begin
            delvar2,coladd
            ptr=(*ptrFit).next
            while ptr_valid(ptr) do begin
                if (*ptr).include and (*ptr).type ne 10 then begin
                    ret=hkllabels(ptr,list,sd,O = k ne 0)
                    
                    ind=value_locate(x,ret.tt)
                    retI=y2[ind];Already converted to xtype
                    
                    if list.xtype ne 1 then begin
                        tmp1=CHI_xval(ret.tt,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,1)
                        ret.tt=tmp1[*,list.xtype]
                    endif
                    
                    if ~bsingle then begin
                        if n_elements(coladd) eq 0 then coladd=-col[k]+1 else coladd++
                    endif else coladd=0

                    for i=0l,n_elements(retI)-1 do begin
                        powplots,[ret.tt[i],ret.tt[i]],[mi,retI[i]],NOCLIP = 0,color=col[k]+coladd,$
                        /data,thick=ret.thick[i],pow=list.ypower
                        powxyouts,ret.tt[i],retI[i],ret.labels[i],/data,color=col[k]+coladd,NOCLIP = 0,$
                        orientation=90,charthick=ret.thick[i],pow=list.ypower
                    endfor
       
                endif
                ptr=(*ptr).next
            endwhile
        endif
    endif
    delvar2,A
endif
endfor

end;pro fitmodelplot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

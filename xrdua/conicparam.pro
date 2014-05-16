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

pro ArcRecalcCorners,xp,yp,center,xs,ys,valid
; Handle invalid corners
; xp: xa0,xa1,xb0,xb1
; yp: ya0,ya1,yb0,yb1
; valid: va0,va1,vb0,vb1
; a: azi_a, b:azi_b
for i=0,3 do begin
    if valid[i] eq 0 then begin
        j=i-2*(i mod 2)+1 ; the valid one
        rico=(yp[j]-center[1])/(xp[j]-center[0])
        interc=center[1]-rico*center[0]

        xoff=[-interc/rico,xs,(ys-interc)/rico,0]
        yoff=[0,rico*xs+interc,ys,interc]
        
        ind=where(xoff le xs+1 and xoff ge 0 and yoff le ys+1 and yoff ge 0,ct)
        if ct eq 0 then begin
            xoff=[xp[j],xp[j]]
            yoff=[0.,ys]
        endif else begin
            xoff=xoff[ind]
            yoff=yoff[ind]
        endelse

        sx1=sign(xp[j]-center[0])
        sy1=sign(yp[j]-center[1])
        sx2=sign(xoff[0]-center[0])
        sy2=sign(yoff[0]-center[1])
        
        if sx1 eq sx2 and sy1 eq sy2 then ind=0 else ind=1
        xp[i]=xoff[ind]
        yp[i]=yoff[ind]
    endif
endfor
end;pro ArcRecalcCorners
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ConicTypeTT,a,tt,str=str

; 0: ellipse
; 1: hyperbola
; 2: parabola
; 3: point, line or double line
; 4: error

r=(sin(a))^2/(cos(tt))^2

type=4b
if r lt 1 then type=0b
if r gt 1 then type=1b
if r eq 1 then type=2b
if abs(a) eq !pi/2 or tt eq !pi/2. then type=3b

case type of
0:    str='                     '
1:    str='Warning: hyperbola!!!'
2:    str='Warning: parabola!!! '
3:    str='Warning: line!!!     '
4:    str='Warning: undefined!!!'
endcase

return,type
end;function ConicTypeTT
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function azMmtoPixel,phi,scrx,scry
if scrx eq scry then return,phi
return,atantan(phi,scrx/scry)
end;function azMmtoPixel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function azPixeltoMm,phi,scrx,scry
if scrx eq scry then return,phi
return,atantan(phi,scry/scrx)
end;function azPixeltoMm
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function azTiltToCone,azIN,dega
; in and out in mm
return,atantan(azIN,cos(dega))
end;function azTiltToCone
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function azConeToTilt,azIN,dega
; in and out in mm
return,atantan(azIN,cos(dega),/divide)
end;function azConeToTilt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AziMinInc,a,b,dist

; Since there is no such thing as "the minimum increment angle"
; (because this depends on two-theta), we'll take the angle along
; the major axis and at the shortest side cut off by the beam center.
; The true minima are lying on the left and on the right of this angle.
; How far to the left/right, depends on a and tt-difference.

; Start from b (angle from Y-axis)
degb=b-!dpi/2 ; 90 counter-clockwise (in b-angle frame)
if sin(a) lt 0 xor dist lt 0 then degb+=!dpi ; 180 clockwise

; map between [360,0[
add=degb-(degb mod (2*!dpi))
ind=where(degb le 0,ct)
if ct ne 0 then add[ind]-=2*!dpi ; 360 clockwise (see next line)
degb-=add ; add counter-clockwise

; map between [0,360[
degb=2*!dpi-degb ;(clockwise becomes counter-clockwise and v.v.)

return,degb ; in mm-frame
end;function AziMinInc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MinIncr,R
;return,acos(1-0.5/R/R)
return,0.01
end;function MinIncr
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function azConeRange,bazi,eazi
inc=MinIncr()
if eazi lt bazi then inc*=(-1)
delta=eazi-bazi
nphi=ceil(delta/inc)+1
if nphi eq 1 then return,bazi
inc=delta/(nphi-1.)
return,inc*lindgen(nphi)+bazi
end;function azConeRange
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AzimuthalShift,phipixphi0,scrx,scry,dega,degb
; in in pixel
; out in mm
; 
; delta = tan(phidphi0+degb).cos(dega)
return,azTiltToCone(azPixeltoMm(phipixphi0,scrx,scry)+degb,dega)
end;function AzimuthalShift
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function azTiltToConeShift,azIN,phipixphi0,scrx,scry,dega,degb
; in and out in mm
return,azTiltToCone(azIN+degb,dega)-AzimuthalShift(phipixphi0,scrx,scry,dega,degb)
end;function azTiltToConeShift
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function azConeShiftToTilt,azIN,phipixphi0,scrx,scry,dega,degb
; in and out in mm
return,azConeToTilt(azIN+AzimuthalShift(phipixphi0,scrx,scry,dega,degb),dega)-degb
end;function azConeShiftToTilt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function azDetRange,phi0,phi1,phipixphi0,scrx,scry,dega,degb
; phi0,phi1: azimuth in detector frame (mm)

; phidphi0 shift not needed
azCFrame=azTiltToCone(double([phi0,phi1])+degb,dega)
azCFrame=azConeRange(azCFrame[0],azCFrame[1]) ; without the shift
return,float(azConeToTilt(azCFrame,dega)-degb)

end;function azDetRange
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RealToCalibAngles,dega,beta,degb

; dega (alpha) : -90 -> 90
; beta (beta)  : -90 -> 90
; degb (-gamma): -90 -> 90

cosalpha2 = cos(dega)*cos(beta) ; 0 -> 1

if cosalpha2 eq 1 then begin
    dega=0
    beta=0
    delta=0
endif else begin
    
    tangamma2 = (tan(degb)*tan(dega)-sin(beta))/(tan(dega)+tan(degb)*sin(beta))
    tangammadelta = sin(beta)/tan(dega) ; shift in gamma due to beta alone
    
    gamma2=atan(tangamma2)  ; -90 -> 90
    ; Gamma shift due to beta can cause the conic section to look like -dega
    binvert = ~floatequal(float(degb-atan(tangammadelta)),float(gamma2))
    
    if binvert xor (dega lt 0) then alpha2=-acos(cosalpha2) else alpha2=acos(cosalpha2)  ; -90 -> 90

    asindelta = -sin(beta)/sin(alpha2)
    delta = asin(asindelta)  ; -90 -> 90
    
    if binvert then begin
        delta=!dpi-delta ; 90 -> 270
        if sin(delta) lt 0 then delta-=2*!dpi ; 90 -> 180 and -90 -> -180
    endif
    
    dega=alpha2
    beta=0
    degb=gamma2
endelse

return,delta
end;function RealToCalibAngles
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CalibToRealAngles,dega,degb,phidphi0

; dega (alpha) : -90 -> 90
; degb (-gamma): -90 -> 90
; phidphi0     : 0 -> 360

delta = azTiltToCone(double(phidphi0+degb),dega)

binvert = cos(delta) lt 0
if binvert then begin
    delta=!dpi-delta
    if delta gt !dpi/2 then delta-=2*!dpi
    if delta lt !dpi/2 then delta+=2*!dpi
endif
; delta: -90 -> 90

sinbeta = -sin(delta)*sin(dega)
beta = asin(sinbeta)  ; -90 -> 90

cosalpha = cos(dega)/cos(beta)
if binvert xor (dega lt 0) then alpha = -acos(cosalpha) else alpha = acos(cosalpha)

tanmingamma = (tan(degb)*tan(alpha)+sinbeta)/(tan(alpha)-tan(degb)*sinbeta)
mingamma = atan(tanmingamma)  ; -90 -> 90

dega=alpha
degb=mingamma

return,beta
end;function CalibToRealAngles
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Getphipixphi0,phipixphiX,phiX,dega,degb,scrx,scry
; tan(phiX+delta) = tan(phidphiX+degb).cos(alpha)
; delta = atan(tan(phidphiX+degb).cos(dega))-phiX
;       = atan(tan(phidphi0+degb).cos(dega))-0
;       
; <=> tan(phidphi0+degb).cos(dega) = tan(atan(tan(phidphiX+degb).cos(dega))-phiX)

return,azMmtoPixel(atantan(atantan(azPixeltoMm(double(phipixphiX),scrx,scry)+degb,cos(dega))-phiX,cos(dega),/divide)-degb,scrx,scry)

end;function Getphipixphi0
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TTPConvertionConstants,phid,phipixphi0,scrx,scry,dega,degb,dist,K1=K1,K2=K2,K3=K3,K4=K4,Kpol=Kpol,Katt=Katt

; Remark:
; if dega=0 then K1=1, K2=0, K3=1, K4=0

; phid: angle in detector frame (mm frame)
; phipixphi0: angle in the detector frame (pixel frame) which corresponds to azCFrame=0

; Angle in tilted frame
azTFrame=double(phid+degb)

; Angle in cone frame (without offset!)
azCFrame=azTiltToCone(azTFrame,dega)

; R = sqrt(xd^2 + yd^2)
;   = sqrt(xt^2 + yt^2)
;    = |xt/cos(azTFrame)|
;    = |xc(u)/cos(azTFrame)|
;    = ... ; filled in xc and u
T=cos(azTFrame)
N=cos(azCFrame)
ind=where(abs(N) lt 0.0001,ct)
if ct ne 0 then begin
    ; cfr. above:
    ; cos(azTFrame)/cos(azCFrame) = sin(azTFrame)*cos(a)/cos(azCFrame)
    T[ind]=sin(azTFrame[ind])*cos(dega)
    N[ind]=sin(azCFrame[ind])
endif
K1=float(T/N)
K2=float(-sin(azTFrame)*sin(dega))
K3=float(cos(dega))
K4=float(sin(azCFrame)*sin(dega))

; Offset due to calibration with 2 angles instead of 3
azCFrame-=AzimuthalShift(phipixphi0,scrx,scry,dega,degb)
Kpol=float(cos(2*azCFrame))
Katt=float(azCFrame)
end;pro TTPConvertionConstants
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LorentzFactor,tt,Linfo

case Linfo.Ltype of
1:    begin
    return,1/sin(tt/2)
    endcase
else: return,1. ; Handle in 1D
endcase

end;function LorentzFactor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PolarizationFactor,tt,Kpol,lambda,Pinfo

; Pinfo:
; P: linear degree of polarization synchrotron = (Iy-Ix)/(Ix+Iy)
; s: -1 = horizontally diffracting monochromator, +1 = vertically diffracting monochromator
; n: number of monochromator crystals
; mmult: multiplier for the monochromator crystal (mosaic: mmult=2, perfect: m=1)
; dmono: d-spacing of monochromator reflection

; phi=findgen(361)*!pi/180
; Kpol = cos(2*phi)
; tt=20*!pi/180
; Pinfo={Ptype:1,Plindeg:0.9,n:2,mmult:1.,dmono:3.1356,s:1}
; P=PolarizationFactor(tt,Kpol,0.4,Pinfo)
; window,0 & plot,phi*180/!pi,P,/xs

case Pinfo.Ptype of
1:    begin
    ; Plindeg between 0 and 1
    K = (1-Pinfo.Plindeg*Kpol)/2
    return, 1-K+K*(cos(tt))^2
    endcase
2:    begin
    m = Pinfo.mmult*Pinfo.n
    mono = (abs(cos(2*asin(lambda/(2*Pinfo.dmono)))))^m
    Plindeg = (Pinfo.s*(1+Pinfo.s*Pinfo.Plindeg)-Pinfo.s*mono*(1-Pinfo.s*Pinfo.Plindeg))/((1+Pinfo.s*Pinfo.Plindeg)+mono*(1-Pinfo.s*Pinfo.Plindeg))
    
    ; Plindeg between 0 and 1
    K = (1-Plindeg*Kpol)/2
    return, 1-K+K*(cos(tt))^2
    endcase
else: return,1. ; polarization independent of azimuth, handle in 1D
endcase

end;function PolarizationFactor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AttenuationFactor,tt,phi,Ainfo

; Ainfo:
;  flatnormphi: the azimuth of the surface normal in the cone frame 
;  flatnormpolar: the angle between beam and surface normal
;  flatD (mm): flat plate sample thickness
;  muL (cm^(-1)): linear attenuation coeff. of the mixture

; phi=findgen(361)*!pi/180
; tt=20*!pi/180
; Ainfo={atype:1,flatnormphi:!pi,flatnormpolar:85*!pi/180,flatD:10.,muL:1.}
; A=AttenuationFactor(tt,phi,Ainfo)
; window,0 & plot,phi*180/!pi,A,/xs

; Test rounding errors close to zero:
; Ainfo={atype:1,flatnormphi:0.,flatnormpolar:0*!pi/180,flatD:10.,muL:1.}
; tt=findgen(91)*!pi/180
; phi=20*!pi/180
; A=AttenuationFactor(tt,phi,Ainfo)
; window,0 & plot,tt*180/!pi,A,/xs,psym=1,xrange=[-5,10],/ys,yrange=[0.365,0.37]
; costt=cos(tt)
; y=costt/(costt-1)*(exp(-Ainfo.muL*Ainfo.flatD/10./costt)-exp(-Ainfo.muL*Ainfo.flatD/10.) )
; oplot,tt*180/!pi,y,color=100
; plots,tt*180/!pi,Ainfo.muL*Ainfo.flatD/10.*exp(-Ainfo.muL*Ainfo.flatD/10./costt)

case Ainfo.Atype of
1:    begin ; flat homogeneous sample
    ; reflection (cosb<0): xe=0
    ; transmission (cosb>0): xe=D
    
    if Ainfo.muL eq 0 or Ainfo.flatD eq 0 then return,1.
    
    costt=cos(tt)
    cosb=cos(phi-Ainfo.flatnormphi)*sin(Ainfo.flatnormpolar)*sin(tt)+cos(Ainfo.flatnormpolar)*costt
    cosa=cos(Ainfo.flatnormpolar)
    
    chi=Ainfo.muL*(1/cosa-1/cosb)

    A=(1-exp(-chi*Ainfo.flatD/10))/(1-cosa/cosb)
    ind=where(finite(A,/nan),ct)
    if ct ne 0 then A[ind]=Ainfo.muL*Ainfo.flatD/10/cosa ; cosa=cosb
    
    ind=where(cosb gt 0,ct) ; transmission
    if ct ne 0 then A[ind]*=exp(-Ainfo.muL*Ainfo.flatD/10/cosb[ind])
    
    ind=where(finite(A,/nan),ct)
    if ct ne 0 then A[ind]=0
    return,A ; between 0 and 1

    endcase
else: return,1. ; attenuation independent of azimuth
endcase

end;function AttenuationFactor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DefaultIntCor2Dinfo
return,{fftype:0,Pff:1e-4,$ ; no flat field correction and Pff Eunits/DU
        Ltype:0,$
        Ptype:0,Plindeg:0.9,n:2,mmult:1.,dmono:3.1356,s:1,$ ; vertical diffracting Si(111) DCM at a synchrotron with P linear = 0.9
        Atype:0,flatnormphi:0.,flatnormpolar:0.,flatD:0.1,muL:5.} ; flat sample perpendicular to the beam, 100um thick and muL=5 cm^(-1)
end;function DefaultIntCor2Dinfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LabelPff,Pff
case Pff of
1:      return,'Total power in Energy units: '
2:      return,'Power in Energy units / m^2: '
else: return,'Energy units/Detector units: '
endcase
end;function LabelPff
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IntCor2D,tt,dist,scrx,scry,IntCor2Dinfo,K3,K4,Kpol,Katt,lambda
; dist: m
; srcx,scry: m
; Eunits: energy units (e.g. Joule)
; DU: detector units

; Lorentz factor (at least part of it)
Lfactor = LorentzFactor(tt,IntCor2Dinfo)

; Polarization factor
Pfactor = PolarizationFactor(tt,Kpol,lambda,IntCor2Dinfo)

; Attenuation factor
Afactor = AttenuationFactor(tt,Katt,IntCor2Dinfo)

case IntCor2Dinfo.fftype of
1:    begin ; spherical (isotropic) flood field
    ; Pff: total power of the flood field (Eunits/s)
    Jcor = IntCor2Dinfo.Pff/(4*!pi)
    endcase
2:    begin ; plane (isotropic) flood field
    ; Pff: power of the flood field (Eunits/s/m^2)
    Jcor = IntCor2Dinfo.Pff*dist^2*K3^2/((K3*cos(tt)+K4*sin(tt))^3)
    endcase
else:begin ; no flood field
    ; Pff: Eunits/DU
    Jcor = IntCor2Dinfo.Pff*dist^2*K3^2/(scrx*scry*(K3*cos(tt)+K4*sin(tt))^3)
    endelse
endcase

return,Jcor/(Lfactor*Pfactor*Afactor)

end;function IntCor2D
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PixelNegR,x1,x2,a,b,dist,scrx,scry,rphi=rphi

if keyword_set(rphi) then begin
    ; x1: radius in detector frame
    ; x2: azimuth in detector frame
    xd=abs(x1)
    yd=xd*sin(x2)
    xd*=cos(x2)
endif else begin
    ; x1: x-coord in detector frame (-center) in pixels
    ; x2: y-coord ...
    xd=x1
    yd=x2
endelse

; Sign of R => Further from the center than the 90 degrees line: R < 0

if (b mod (!pi/2)) eq 0 and (b mod !pi) ne 0 then begin
    ; b = 90
    xd_90=dist*sin(b)/(sin(a)*scrx)
    return,(xd_90 gt 0 and xd gt xd_90) or $
                (xd_90 lt 0 and xd lt xd_90)
endif else begin
    yd_90=-tan(b)*scrx/scry*xd+dist/(sin(a)*cos(b)*scry)
    return,(yd_90 gt 0 and yd gt yd_90) or $
                (yd_90 lt 0 and yd lt yd_90)
endelse
end;function PixelNegR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggReform,A,nx,ny
n=n_elements(A)
if n eq 1 then return,make_array(nx,ny,value=A) else $
if n eq nx*ny then return,reform(A,nx,ny) else $
return,rebin(A,nx,ny,/sample)
end;function BraggReform
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggTtoX,tt,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,onlyd=onlyd,onlyp=onlyp,$
                    onlyq=onlyq,nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi

; scr: m/pixel
; D: m
; lambda: angstrom
; phi,dega,degb: rad
;
; phi: always in the (non-orthogonal) detector frame 
;        in pixel units when pixm=0 and convertion to P, otherwise in mm
;      it can be invalid when a hyperbola => validphi
;      
; P: (no need for convertion between the separate P types)
;     pixm=1, nonortp=1: (Rd(mm),phi(mm)) in detector
;     pixm=0, nonortp=1: (Rd(pixels),phi(pixels)) in detector
;     pixm=1, nonortp=0: (Rort(mm),phi(mm)) in virtual orthogonal detector
;     pixm=0, nonortp=0: (Rort(pixels),phi(pixels)) in virtual orthogonal detector
;     
; tt: in rad or degrees when angledeg
; 
; d: angstrom
; 
; Q: 1/nm
;
; The phi can be convert from mm to pixel and back when P is involved in the transformation:
;     P -> X : convert phi(pixels) to phi(mm) if not already in mm
;     X -> P : convert ? to pixm if not already pixm
;
; Order: P,T,D,Q

; Prepare return value
bool=[keyword_set(onlyp),keyword_set(onlyd),keyword_set(onlyq)]
tbool=total(bool,/pres)
if tbool eq 0 then begin
    bool=[1b,1b,1b]
    tbool=3
endif

ret=make_array(n_elements(tt)>n_elements(phi),tbool,value=tt[0]*0)
i=0

; ----Two-theta > Pix----
if keyword_set(angledeg) then mtt=180/!pi else mtt=1.
if bool[0] then begin
    ; Default azimuth in detector frame
    if n_elements(phi) eq 0 then phi=AziMinInc(dega,degb,dist)
    
    if keyword_set(pixm) then mPix=1000. else $
        if scrx ne scry then mPix=sqrt((cos(phi)/scrx)^2.+(sin(phi)/scry)^2.) else mPix=1./scrx
    
    if keyword_set(nonortp) then begin
        TTPConvertionConstants,phi,phipixphi0,scrx,scry,dega,degb,dist,K1=K1,K2=K2

        Pix=(dist*mPix)/(K1/tan(tt/mtt)-K2)
        validphi=Pix ge 0
        
        ret[*,i++]=temporary(Pix)
    endif else begin
        ret[*,i++]=(dist*mPix)*tan(tt/mtt)
    endelse
    
    if not keyword_set(pixm) then phi=azMmtoPixel(phi,scrx,scry)
endif

; ----Two-theta > D----
if bool[1] then ret[*,i++] = lambda/(2.*sin(tt/(mtt*2.)))

; ----Two-theta > Q----
if bool[2] then ret[*,i++] = (40*!pi/lambda)*sin(tt/(mtt*2.))

if n_elements(ret) eq 1 then return,ret[0] else return,ret
end;function BraggTtoX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggTtoX_I,tt,spe,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,onlyd=onlyd,onlyp=onlyp,$
                    onlyq=onlyq,nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi

; For a transformation H: x -> y
; I(x) = I(H^(-1)(y)).JH^(-1)(y)
;      = I(x).1/JH(x)
; For one dimension: JH(x)=dy/dx
; We return 1./(dy/dx)

; Prepare return value
bool=[keyword_set(onlyp),keyword_set(onlyd),keyword_set(onlyq)]
tbool=total(bool,/pres)
if tbool eq 0 then begin
    bool=[1b,1b,1b]
    tbool=3
endif

s=size(spe)
n=n_elements(tt)>n_elements(phi)
if s[0] eq 2 then nspe=s[2] else nspe=1

ret0=make_array(n,value=tt[0]*0)
ret1=BraggReform(spe,n,nspe)

case tbool of
1:ret={x1:ret0,y1:ret1}
2:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1}
3:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1,x3:ret0,y3:ret1}
endcase

i=0

; ----Two-theta > Pix----
if keyword_set(angledeg) then mtt=180/!pi else mtt=1.
if bool[0] then begin
    ; Default azimuth in detector frame
    if n_elements(phi) eq 0 then phi=AziMinInc(dega,degb,dist)
    
    if keyword_set(pixm) then mPix=1000. else $
        if scrx ne scry then mPix=sqrt((cos(phi)/scrx)^2.+(sin(phi)/scry)^2.) else mPix=1./scrx
    
    if keyword_set(nonortp) then begin
        TTPConvertionConstants,phi,phipixphi0,scrx,scry,dega,degb,dist,K1=K1,K2=K2
        ret0=K1/tan(tt/mtt)-K2
        ret.(i)=(dist*mPix)/ret0
        validphi=ret.(i) ge 0
        sintt=sin(tt/mtt)
        ret1=(mtt/(dist*mPix*K1))*(ret0*sintt)^2.
    endif else begin
        ret.(i)=(dist*mPix)*tan(tt/mtt)
        ret1=(mtt/(dist*mPix))*(cos(tt/mtt))^2.
    endelse
    
    if not keyword_set(pixm) then phi=azMmtoPixel(phi,scrx,scry)
    
    ret.(i+1)*=BraggReform(abs(ret1),n,nspe)
    i+=2
endif

; ----Two-theta > D----
if bool[1] then begin
    ret0=tt/(mtt*2.)
    ret.(i)=lambda/(2.*sin(ret0))
    ret.(i+1)*=BraggReform(abs((4.*mtt)/lambda*(sin(ret0))^2./cos(ret0)),n,nspe)
    i+=2
endif

; ----Two-theta > Q----
if bool[2] then begin
    ret0=tt/(mtt*2.)
    ret.(i)=(40*!pi/lambda)*sin(ret0)
    ret.(i+1)*=BraggReform(abs((mtt*lambda)/(20*!pi)/cos(ret0)),n,nspe)
    i+=2
endif

return,ret
end;function BraggTtoX_I
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggQtoX,q,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,onlyd=onlyd,onlyp=onlyp,$
                onlyt=onlyt,nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi

; Prepare return value
bool=[keyword_set(onlyp),keyword_set(onlyt),keyword_set(onlyd)]
tbool=total(bool,/pres)
if tbool eq 0 then begin
    bool=[1b,1b,1b]
    tbool=3
endif

ret=make_array(n_elements(q)>n_elements(phi),tbool,value=q[0]*0)
i=tbool-1

; ----Q > D----
if bool[2] then ret[*,i--] = (20*!pi)/q

if bool[0] or bool[1] then begin
    ; ----Q > Two-theta----
    if keyword_set(angledeg) then mtt=180/!pi else mtt=1.
    tt=2*mtt*asin((lambda/(40.*!pi))*q)
    if bool[1] then ret[*,i--] = tt
    
    ; ----Two-theta > Pix----
    if bool[0] then ret[*,i--]=BraggTtoX(tt,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyp,$
                                        nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi)
endif

if n_elements(ret) eq 1 then return,ret[0] else return,ret
end;function BraggQtoX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggQtoX_I,q,spe,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,onlyd=onlyd,onlyp=onlyp,$
                onlyt=onlyt,nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi

; Prepare return value
bool=[keyword_set(onlyp),keyword_set(onlyt),keyword_set(onlyd)]
tbool=total(bool,/pres)
if tbool eq 0 then begin
    bool=[1b,1b,1b]
    tbool=3
endif

s=size(spe)
n=n_elements(q)>n_elements(phi)
if s[0] eq 2 then nspe=s[2] else nspe=1

ret0=make_array(n,value=q[0]*0)
ret1=BraggReform(spe,n,nspe)

case tbool of
1:ret={x1:ret0,y1:ret1}
2:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1}
3:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1,x3:ret0,y3:ret1}
endcase

i=2*(tbool-1)

; ----Q > D----
if bool[2] then begin
    ret.(i)=(20*!pi)/q
    ret.(i+1)*=BraggReform(q*q/(20*!pi),n,nspe)
    i-=2
endif

if bool[1] or bool[0] then begin
    ; ----Q > Two-theta----
    if keyword_set(angledeg) then mtt=180/!pi else mtt=1
    ret0=2*mtt*asin((lambda/(40.*!pi))*q)
    ret1=BraggReform(sqrt(1-lambda*lambda/(1600*!pi*!pi)*q*q)/(mtt*lambda/(20*!pi)),n,nspe)
    if bool[1] then begin
        ret.(i)=ret0
        ret.(i+1)*=ret1
        i-=2
    endif
    
    ; ----Two-theta > Pix----
    if bool[0] then begin
        ret1=BraggTtoX_I(ret0,ret1,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyp,$
                        nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi)
        ret.(i)=ret1.x1
        ret.(i+1)*=ret1.y1
        i-=2
    endif
endif

return,ret
end;function BraggQtoX_I
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggDtoX,D,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,onlyt=onlyt,onlyp=onlyp,$
                onlyq=onlyq,nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi

; Prepare return value
bool=[keyword_set(onlyp),keyword_set(onlyt),keyword_set(onlyq)]
tbool=total(bool,/pres)
if tbool eq 0 then begin
    bool=[1b,1b,1b]
    tbool=3
endif

ret=make_array(n_elements(D)>n_elements(phi),tbool,value=D[0]*0)
i=tbool-1

; ----D > Q----
if bool[2] then ret[*,i--]=20*!pi/D

if bool[1] or bool[0] then begin
    ; ----D > Two-theta----
    if keyword_set(angledeg) then mtt=180/!pi else mtt=1
    tt=2*mtt*asin(lambda/2./D)
    if bool[1] then ret[*,i--]=tt

    ; ----Two-theta > Pix----
    if bool[0] then ret[*,i--]=BraggTtoX(tt,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyp,$
                                        nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi)
endif

if n_elements(ret) eq 1 then return,ret[0] else return,ret
end;function BraggDtoX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggDtoX_I,D,spe,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,onlyt=onlyt,onlyp=onlyp,$
                onlyq=onlyq,nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi

; Prepare return value
bool=[keyword_set(onlyp),keyword_set(onlyt),keyword_set(onlyq)]
tbool=total(bool,/pres)
if tbool eq 0 then begin
    bool=[1b,1b,1b]
    tbool=3
endif

s=size(spe)
n=n_elements(D)>n_elements(phi)
if s[0] eq 2 then nspe=s[2] else nspe=1

ret0=make_array(n,value=D[0]*0)
ret1=BraggReform(spe,n,nspe)

case tbool of
1:ret={x1:ret0,y1:ret1}
2:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1}
3:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1,x3:ret0,y3:ret1}
endcase

i=2*(tbool-1)

; ----D > Q----
if bool[2] then begin
    ret.(i)=20*!pi/D
    ret.(i+1)*=BraggReform(D*D/(20*!pi),n,nspe)
    i-=2
endif

if bool[1] or bool[0] then begin
    ; ----D > Two-theta----
    if keyword_set(angledeg) then mtt=180/!pi else mtt=1
    ret0=2*mtt*asin(lambda/2./D)
    ret1=BraggReform(sqrt(D^4.-D^2.*(lambda^2./4))/(lambda*mtt),n,nspe)
    if bool[1] then begin
        ret.(i)=ret0
        ret.(i+1)*=ret1
        i-=2
    endif
    
    ; ----Two-theta > Pix----
    if bool[0] then begin
        ret1=BraggTtoX_I(ret0,ret1,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyp,$
                        nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi)
        ret.(i)=ret1.x1
        ret.(i+1)*=ret1.y1
        i-=2
    endif
endif

return, ret
end;function BraggDtoX_I
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggPtoX,Pix,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,onlyd=onlyd,onlyt=onlyt,$
                nonortp=nonortp,onlyq=onlyq,pixm=pixm,angledeg=angledeg,validphi=validphi

; Prepare return value
bool=[keyword_set(onlyt),keyword_set(onlyd),keyword_set(onlyq)]
tbool=total(bool,/pres)
if tbool eq 0 then begin
    bool=[1b,1b,1b]
    tbool=3
endif

ret=make_array(n_elements(Pix)>n_elements(phi),tbool,value=Pix[0]*0)
i=0

; ----Pix > Two-theta----
if keyword_set(angledeg) then mtt=180/!pi else mtt=1

; Default azimuth in detector frame
if n_elements(phi) eq 0 then phi=AziMinInc(dega,degb,dist) $
else if not keyword_set(pixm) then phi=azPixeltoMm(phi,scrx,scry)
if keyword_set(pixm) then mPix=1000. else $
    if scrx ne scry then mPix=sqrt((cos(phi)/scrx)^2.+(sin(phi)/scry)^2.) else mPix=1./scrx

if keyword_set(nonortp) then begin
    TTPConvertionConstants,phi,phipixphi0,scrx,scry,dega,degb,dist,K1=K1,K2=K2
    tt=atan(K1/(K2+(dist*mPix)/Pix))
    validphi=Pix ge 0
endif else begin
    tt=atan(Pix/(dist*mPix))
endelse
tt+=!pi*(tt lt 0)
tt*=mtt
if bool[0] then ret[*,i++]=tt

; ----Two-theta > D/Q----
if bool[1] or bool[2] then begin
    add=bool[1]+bool[2]
    ret[*,i:i+add-1]=BraggTtoX(tt,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,onlyd=bool[1],$
                                onlyq=bool[2],nonortp=nonortp,pixm=pixm,angledeg=angledeg)
    i+=add
endif

if n_elements(ret) eq 1 then return,ret[0] else return,ret
end;function BraggPtoX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggPtoX_I,Pix,spe,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,onlyd=onlyd,onlyt=onlyt,$
                onlyq=onlyq,nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi

; Prepare return value
bool=[keyword_set(onlyt),keyword_set(onlyd),keyword_set(onlyq)]
tbool=total(bool,/pres)
if tbool eq 0 then begin
    bool=[1b,1b,1b]
    tbool=3
endif

s=size(spe)
n=n_elements(Pix)>n_elements(phi)
if s[0] eq 2 then nspe=s[2] else nspe=1

ret0=make_array(n,value=Pix[0]*0)
ret1=BraggReform(spe,n,nspe)

case tbool of
1:ret={x1:ret0,y1:ret1}
2:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1}
3:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1,x3:ret0,y3:ret1}
endcase

i=0

; ----Pix > Two-theta----
if keyword_set(angledeg) then mtt=180/!pi else mtt=1
; Default azimuth in detector frame
if n_elements(phi) eq 0 then phi=AziMinInc(dega,degb,dist) $
else if not keyword_set(pixm) then phi=azPixeltoMm(phi,scrx,scry)
if keyword_set(pixm) then mPix=1000. else $
    if scrx ne scry then mPix=sqrt((cos(phi)/scrx)^2.+(sin(phi)/scry)^2.) else mPix=1./scrx

if keyword_set(nonortp) then begin
    TTPConvertionConstants,phi,phipixphi0,scrx,scry,dega,degb,dist,K1=K1,K2=K2
    keep=(dist*mPix)/Pix+K2
    ret0=atan(K1/keep)
    ret1=(1./(mtt*dist*mPix*K1))*Pix^2.*(keep^2.+K1^2.)
    validphi=Pix ge 0
endif else begin
    keep=Pix/(dist*mPix)
    ret0=atan(keep)
    ret1=(dist*mPix/mtt)*(1+keep^2.)
endelse
ret0+=!pi*(ret0 lt 0)
ret0*=mtt
ret1=BraggReform(abs(ret1),n,nspe)

if bool[0] then begin
    ret.(i)=ret0
    ret.(i+1)*=ret1
    i+=2
endif

; ----Two-theta > D/Q----
if bool[1] or bool[2] then begin
    ret1=BraggTtoX_I(ret0,ret1,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,onlyd=bool[1],$
                    onlyq=bool[2],nonortp=nonortp,pixm=pixm,angledeg=angledeg)
    j=0
    for k=1,2 do $
    if bool[k] then begin
        ret.(i)=ret1.(j)
        ret.(i+1)*=ret1.(j+1)
        i+=2
        j+=2
    endif
endif

return,ret
end;function BraggPtoX_I
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggXYtoX,x,y,lambda,dist,scrx,scry,dega,degb,phipixphi0,center,onlyp=onlyp,onlyd=onlyd,$
    onlyt=onlyt,onlyq=onlyq,nonortp=nonortp,pixm=pixm,angledeg=angledeg,aziout=phi

; Prepare return value
bool=[keyword_set(onlyp),keyword_set(onlyt),keyword_set(onlyd),keyword_set(onlyq)]
tbool=total(bool,/pres)
if tbool eq 0 then begin
    bool=[1b,1b,1b,1b]
    tbool=4
endif

nx=n_elements(x)
ret=make_array(nx,tbool,value=x[0]*0.)
i=0

; ----X,Y > Pix----
xd=x-center[0]
yd=y-center[1]
if keyword_set(pixm) then begin
    xd*=scrx*1000
    yd*=scry*1000
endif
Pix=sqrt(xd*xd+yd*yd) ; nonort distance (mm or pixel frame), so should be positive
phi=atan(yd,xd) ; azimuth on detector (mm or pixel frame)

; ----Pix(nonort) > Pix----
if bool[0] then begin
    if ~keyword_set(nonortp) then begin
        ; From non-orthogonal to orthogonal (pixm: both in mm, ~pixm: both in pixel)
        phi_=phi
        if not keyword_set(pixm) then phi_=azPixeltoMm(phi_,scrx,scry)
        ; _phi in the mm-frame
        
        if keyword_set(pixm) then mPix=1000. else $
            if scrx ne scry then mPix=sqrt((cos(phi_)/scrx)^2.+(sin(phi_)/scry)^2.) else mPix=1./scrx
        TTPConvertionConstants,phi_,phipixphi0,scrx,scry,dega,degb,dist,K1=K1,K2=K2
        ret[*,i++]=K1/(K2/(dist*mPix)+1./Pix)
    endif else ret[*,i++]=Pix
endif
; phi in mm or pixels

; ----Pix(nonort) > tt/D/Q----
if bool[1] or bool[2] or bool[3] then begin
    add=bool[1]+bool[2]+bool[3]
    phi_=phi
    ret[*,i:i+add-1]=BraggPtoX(Pix,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi_,/nonortp,$
                onlyt=bool[1],onlyd=bool[2],onlyq=bool[3],pixm=pixm,angledeg=angledeg)
    if ~bool[0] then phi=phi_
    i+=add
endif

if n_elements(ret) eq 1 then return,ret[0] else return,ret

end;function BraggXYtoX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggXYtoX_I,x,y,spe,lambda,dist,scrx,scry,dega,degb,phipixphi0,center,IntCor2Dinfo,onlyp=onlyp,onlyd=onlyd,onlyt=onlyt,$
    onlyq=onlyq,nonortp=nonortp,pixm=pixm,angledeg=angledeg,aziout=phi

; Prepare return value
bool=[keyword_set(onlyp),keyword_set(onlyt),keyword_set(onlyd),keyword_set(onlyq)]
tbool=total(bool,/pres)
if tbool eq 0 then begin
    bool=[1b,1b,1b,1b]
    tbool=4
endif

; ----X,Y > Pix(nonort)----
xd=x-center[0]
yd=y-center[1]
if keyword_set(pixm) then begin
    xd*=scrx*1000
    yd*=scry*1000
endif
ret0=sqrt(xd*xd+yd*yd) ; nonort distance (mm or pixel frame), so should be positive
ret1=spe
phi=atan(yd,xd) ; azimuth on detector (mm or pixel frame)

; ----Pix(nonort) > tt(rad)----
; phi in mm or pixels
tt=BraggPtoX(ret0,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyt,/nonortp,pixm=pixm,angledeg=0)
tt=reform(tt,size(ret0,/dim),/overwrite)
; phi in mm
TTPConvertionConstants,phi,phipixphi0,scrx,scry,dega,degb,dist,K1=K1,K2=K2,K3=K3,K4=K4,Kpol=Kpol,Katt=Katt

ret1*=IntCor2D(tt,dist,scrx,scry,IntCor2Dinfo,K3,K4,Kpol,Katt,lambda)

case tbool of
1:ret={x1:ret0,y1:ret1}
2:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1}
3:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1,x3:ret0,y3:ret1}
4:ret={x1:ret0,y1:ret1,x2:ret0,y2:ret1,x3:ret0,y3:ret1,x4:ret0,y4:ret1}
endcase

i=0

; ----Two-theta(rad) > Pix/D/Q----
; phi (mm)
if bool[0] or bool[1] or bool[2] or bool[3] then begin
    if bool[0] or bool[2] or bool[3] then $
    ret1=BraggTtoX_I(reform(tt,n_elements(tt)),reform(ret1,n_elements(ret1)),lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,nonortp=nonortp,$
                onlyp=bool[0],onlyd=bool[2],onlyq=bool[3],pixm=pixm,angledeg=0)
                
    j=0
    for k=0,3 do $
    if bool[k] then begin
        if k eq 1 then begin
            ; ----tt(rad) > tt----
            if keyword_set(angledeg) then begin
                ret.(i)=tt*180/!pi
                ret.(i+1)*=!pi/180
            endif else begin
                ret.(i)=tt
            endelse
        endif else begin
            ret.(i)=ret1.(j)
            ret.(i+1)=ret1.(j+1)
            j+=2
        endelse
        i+=2
    endif
endif

return,ret
end;function BraggXYtoX_I
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggXtoXY,x,phi,lambda,dist,scrx,scry,dega,degb,phipixphi0,center,typein,$
            nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi

; phi is given in the detector frame

; ----X > P----
; Get non-orthogonal R (in pixel or mm, depending on pixm)
case typein of
0:     begin
    if not keyword_set(nonortp) then begin
        ; From orthogonal to non-orthogonal
        phi_=phi
        if n_elements(phi_) eq 0 then phi_=AziMinInc(dega,degb,dist) $
        else if not keyword_set(pixm) then phi_=azPixeltoMm(phi_,scrx,scry)
        if keyword_set(pixm) then mPix=1000. else $
            if scrx ne scry then mPix=sqrt((cos(phi_)/scrx)^2.+(sin(phi_)/scry)^2.) else mPix=1./scrx
        TTPConvertionConstants,phi_,phipixphi0,scrx,scry,dega,degb,dist,K1=K1,K2=K2
        R=1./(K1/x-K2/(dist*mPix))
    endif else R=x
    validphi=R ge 0
    endcase
1:  R=BraggTtoX(x,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyp,/nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi)
2:  R=BraggDtoX(x,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyp,/nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi)
3:  R=BraggQtoX(x,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyp,/nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi)
endcase

; ----[R,phi] > [x,y]----
; R and phi in mm or pixel
xd=R*cos(phi)
yd=R*sin(phi)
if keyword_set(pixm) then begin
    xd/=scrx*1000
    yd/=scry*1000
    phi=azMmtoPixel(phi,scrx,scry)
endif
xd+=center[0]
yd+=center[1]
; phi in pixel

return,[[xd],[yd]]
end;function BraggXtoXY
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BraggXtoXY_I,x,phi,spe,lambda,dist,scrx,scry,dega,degb,phipixphi0,center,typein,IntCor2Dinfo,$
            nonortp=nonortp,pixm=pixm,angledeg=angledeg,validphi=validphi

; phi is given in the detector frame (mm or pixel)

; ----X > Two-theta(rad)----
case typein of
0:     begin
    tt=BraggPtoX_I(x,spe,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyt,nonortp=nonortp,pixm=pixm,angledeg=0)
    ret1=tt.y1
    tt=tt.x1
    endcase
1:     begin
    s=size(spe)
    n=n_elements(x)>n_elements(phi)
    if s[0] eq 2 then nspe=s[2] else nspe=1
    ret1=BraggReform(spe,n,nspe)

    if keyword_set(angledeg) then begin
        tt=x*!pi/180
        ret1=BraggReform(spe*180/!pi,n,nspe)
    endif else begin
        tt=x
        ret1=BraggReform(spe,n,nspe)
    endelse
    endcase
2:     begin
    tt=BraggDtoX_I(x,spe,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyt,angledeg=0)
    ret1=tt.y1
    tt=tt.x1
    endcase
3:     begin
    tt=BraggQtoX_I(x,spe,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyt,angledeg=0)
    ret1=tt.y1
    tt=tt.x1
    endcase
endcase
; phi is given in the mm frame
; ret1 = dx/dtt

; ----Two-theta(rad) > X,Y----
TTPConvertionConstants,phi,phipixphi0,scrx,scry,dega,degb,dist,K3=K3,K4=K4,Kpol=Kpol,Katt=Katt
ret1/=IntCor2D(tt,dist,scrx,scry,IntCor2Dinfo,K3,K4,Kpol,Katt,lambda)

Rpix=BraggTtoX(tt,lambda,dist,scrx,scry,dega,degb,phipixphi0,phi=phi,/onlyP,pixm=0,/nonortp,angledeg=0,validphi=validphi)
; phi in pixels frame

return,{x11:Rpix*cos(phi)+center[0],$
        x12:Rpix*sin(phi)+center[1],$
        y1:ret1}

end;function BraggXtoXY_I
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro BraggXtoXY_I_optimized,x,phi,xmap,ymap,Jcor,Jac1D,$
                        lambda,dist,scrx,scry,dega,degb,phipixphi0,center,$
                        typein,IntCor2Dinfo,validphi=Rpix;,debug=debug

; Same as BraggXtoXY_I but optimized for the following conditions:
;     nonortp=0b
;     pixm=1b
;     angledeg=1b
;    x (row vector): "type" in (mm_ort, degrees, Angstrom or 1/nm)
;    phi (column vector): azimuth in radians (in mm frame)
;    
; Note that this function will be used for the inverse transformation (azint):
;   IX = Icor(X,Y) . Jcor / Jac1D
;     Jac1D = dx/dtt

; ----X > Two-theta(rad)----
case typein of
0:     begin
    ; tt = mtt.atan(x/(dist.1000))
    ; dtt/dx = mtt / (1+x^2/(dist.1000)^2) / (dist.1000)
    Jac1D=x/(dist*1000)
    tt=atan(Jac1D)
    Jac1D^=2
    Jac1D++
    Jac1D*=abs(dist*1000)
    tt+=!pi*(tt lt 0)
    endcase
1:     begin
    ; tt = x/mtt
    ; dtt/dx = 1/mtt
    tt=x*!pi/180
    Jac1D=180/!pi
    endcase
2:     begin
    ; tt = 2.mtt.asin(lambda/2/x)
    ; dtt/dx = -mtt.lambda/sqrt(x^4-lambda^2/4.x^2)
    tt=2*asin(lambda/2/x)
    Jac1D=sqrt(x^4-x^2*(lambda^2/4))/lambda
    endcase
3:     begin
    ; tt = 2.mtt.asin(lambda/40/pi.x)
    ; dtt/dx = mtt/sqrt(1-lambda^2/1600/pi^2.x^2).lambda/20/pi
    tt=2*asin((lambda/(40.*!pi))*x)
    Jac1D=sqrt(1-lambda^2/(1600*!pi^2)*x^2)/(lambda/(20*!pi))
    endcase
endcase

; ----Two-theta(rad) > X,Y----
TTPConvertionConstants,phi,phipixphi0,scrx,scry,dega,degb,dist,K1=K1,K2=K2,K3=K3,K4=K4,Kpol=Kpol,Katt=Katt

axbin=n_elements(x)
aybin=n_elements(phi)
xmap=rebin(tt,axbin,aybin,/sample)
ymap=rebin(phi,axbin,aybin,/sample)

; 2D intensity correction
Jcor=1.

; Polarization factor
if IntCor2Dinfo.Ltype ne 0 then Jcor /= LorentzFactor(xmap,IntCor2Dinfo)

; Polarization factor
if IntCor2Dinfo.Ptype ne 0 then Jcor /= PolarizationFactor(xmap,rebin(Kpol,axbin,aybin,/sample),lambda,IntCor2Dinfo)

; Attenuation factor
if IntCor2Dinfo.Atype ne 0 then Jcor /= AttenuationFactor(xmap,rebin(Katt,axbin,aybin,/sample),IntCor2Dinfo)

; Flood field factor
case IntCor2Dinfo.fftype of
1:    Jcor *= IntCor2Dinfo.Pff/(4*!pi)
2:    Jcor *= IntCor2Dinfo.Pff*dist^2*K3^2/((K3*cos(xmap)+rebin(K4,axbin,aybin,/sample)*sin(xmap))^3)
else:Jcor *= IntCor2Dinfo.Pff*dist^2*K3^2/(scrx*scry*(K3*cos(xmap)+rebin(K4,axbin,aybin,/sample)*sin(xmap))^3)
endcase

;debug={tt:tt,phi:phi,Jac1D:Jac1D,P:1/PolarizationFactor(xmap,rebin(Kpol,axbin,aybin,/sample),lambda,IntCor2Dinfo),$
;        L:1/LorentzFactor(xmap,IntCor2Dinfo),Jcor:IntCor2Dinfo.Pff*dist^2*K3^2/(scrx*scry*(K3*cos(xmap)+rebin(K4,axbin,aybin,/sample)*sin(xmap))^3)}

; Non-orthogonal radius in pixels
Rpix = rebin(K1,axbin,aybin,/sample)
Rpix /= tan(xmap)
Rpix -= rebin(K2,axbin,aybin,/sample)
if scrx ne scry then Rpix = dist*sqrt((cos(ymap)/scrx)^2+(sin(ymap)/scry)^2)/Rpix $
else Rpix = (dist/scrx)/Rpix

; Corresponding azimuth
if scrx ne scry then ymap=atantan(ymap,scrx/scry) ; mm to pixel

; (non-integer) pixel positions
xmap = Rpix*cos(ymap)+center[0]
ymap = Rpix*sin(temporary(ymap))+center[1]

; valid pixels
Rpix ge= 0

end;pro BraggXtoXY_I_optimized
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SemiAxisCoord,tt,a,b,dist,scrx,scry,center

;ytroot1=dist*tan(tt)*(sin(-a)*tan(tt)+cos(-a))/(1.-sin(-a)^2./cos(tt)^2.)
;ytroot2=dist*tan(tt)*(sin(-a)*tan(tt)-cos(-a))/(1.-sin(-a)^2./cos(tt)^2.)

ytroot1=dist/(sin(a)+cos(a)/tan(tt))
ytroot2=dist/(sin(a)-cos(a)/tan(tt))

H=(abs(ytroot1)+abs(ytroot2))/2.
c=(ytroot1>ytroot2)-H
h=sqrt(tan(tt)^2.*(c*sin(a)-dist)^2.-c*c*cos(a)^2.)
p1=[h,c]
p2=[0,ytroot1]
p3=[-h,c]
p4=[0,ytroot2]
rot=[[cos(b),sin(b)],[-sin(b),cos(b)]]
p1=rot##transpose(p1)
p2=rot##transpose(p2)
p3=rot##transpose(p3)
p4=rot##transpose(p4)
p=[p1,p3,p2,p4]
xx=p[*,0]
yy=p[*,1]

return,[[xx/scrx+center[0]],[yy/scry+center[1]]]
end;function SemiAxisCoord
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EllipsePlotCoord,tt,phi,a,b,phipixphi0,dist,scrx,scry,center,$
    forceopen=forceopen,reforminput=reforminput,gapshift=gapshift,nvalid=ctkeep

; phi is defined in the detector frame (mm)
phi_=phi
xy=BraggXtoXY(tt,phi_,0,dist,scrx,scry,a,b,phipixphi0,center,1,validphi=validphi)
if keyword_set(reforminput) then phi=phi_ ; mm -> pixel

; ----Keep only the correct side when a hyperbola----
indkeep=where(validphi,ctkeep)
if ctkeep ne 0 then begin
    xy=xy[indkeep,*]
    if keyword_set(reforminput) then phi=phi[indkeep]
endif else xy=fltarr(2)

; ----Force open when a hyperbola----
gapshift=0
if keyword_set(forceopen) then $
    if ConicTypeTT(a,tt) ne 0 and ctkeep gt 2 then begin
        ; Find the biggest gap in distance between points and break there
        m=max(  total((xy-shift(xy,-1,0))^2.,2),ind)
        gapshift=ctkeep-ind-1
        xy=shift(xy,gapshift,0)
        if keyword_set(reforminput) then phi=shift(phi,gapshift)
    endif
    
return,xy
end;function EllipsePlotCoord
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

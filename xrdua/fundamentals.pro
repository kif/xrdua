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

function speedoflight
return,299792458.;(m/sec)
end;function speedoflight
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function planckconstant
return,4.13566743E-15 ;(eV.sec)
end;function planckconstant
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function avogadroconstant,noexp=noexp
if keyword_set(noexp) then return,0.60221367
return,6.0221367E23 ; (#particles/mole) = # atoms in 12g 12C
end;function avogadroconstant
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function electronradius,noexp=noexp
if keyword_set(noexp) then ret=2.81794092 $
else ret=2.81794092E-15 ;(m)
return,ret
end;function electronradius
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function electronrestenergy
; E = mc^2 (m = rest mass of an electron)
return,0.510998910d6
end;function electronrestenergy
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LambdaEnergy,LorE,kev=kev

; E(keV) = h(eV sec) . c(m/sec)
;          ------------------
;             lambda(Angstrom).10^-10
ret=(planckconstant()*speedoflight()*1E10)/LorE
if keyword_set(kev) then ret/=1E3

return,ret

end;function LambdaEnergy
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LambdaFreq,LorF,n
; n = refractive index of the medium
; c (m/sec) = speed of light in vacuum
; v (m/sec) = speed of light in medium
; L0 (m) = wavelength of light in vacuum
; L (m) = wavelength of light in medium
; F (Hz=1/sec) = frequency of light
; 
;     n = c/v = L0/L
;     v = L.F
;     c = L0.F
; <=> F = c/L0 = c / (n.L)

if n_elements(n) eq 0 then n=1.0003 ;(air)

return,speedoflight()/(n*LorF)

end;function LambdaFreq
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DefaultIon,Z,Ion,repeats
if n_elements(Z) ne 1 then return

if Ion ne 0 or Z lt 0 or Z ge 104 then return

repeats=1

Ions=[0,1,0,$ ; H-He
    1,2,3,4,5,-2,-1,0,$ ; Li-Ne
    1,2,3,4,5,6,-1,0,$ ; Na-Ar
    1,2,3,4,5,6,2,2,2,2,2,2,3,2,5,6,-1,0,$; K-Kr
    1,2,3,4,5,6,4,4,3,2,1,2,3,4,5,6,7,0,$ ; Rb-Xe
    1,2,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,5,6,4,4,4,4,3,2,3,2,3,4,1,0,$; Cs-Rn
    1,2,3,4,5,6,5,4,3,3,3,3,3,3,3,2,3]    ; Fr-Lr

Ion=Ions[Z]
end;pro DefaultIon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ZElement,In,repeats=repeats,ion=ion

Mat1 = [ "",$
              "H" ,"He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na",$
              "Mg", "Al", "Si", "P", "S", "Cl", "Ar", "K", "Ca", "Sc", "Ti",$
              "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",$
              "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru",$
              "Rh", "Pd", "Ag", "Cd", "In", "Sn", "Sb", "Te", "I", "Xe", "Cs"]
Mat2 = [$
             "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb", "Dy",$
              "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir",$
              "Pt", "Au", "Hg", "Tl", "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra",$
              "Ac", "Th", "Pa", "U", "Np", "Pu", "Am", "Cm", "Bk", "Cf", "Es",$
              "Fm", "Md", "No", "Lr",""]
Mat = [temporary(Mat1),temporary(Mat2)]

nmax=n_elements(Mat)-1

if size(In,/type) eq 7 then begin
    repeats=1
    ion=0
    if In eq '' then return,-1
    
;    ; Check for single element without charge
;    Eli=strlowcase(In)
;    Eli=strupcase(strmid(Eli,0,1))+strmid(Eli,1)
;    Z=where(Mat eq Eli,ct)
;    if ct eq 1 then begin
;        In=Eli
;        return,Z
;    endif
    
    ; Extract elements and charge+repeat: e.g. "Zx+y" becomes "Z" and "x+y"
    poschars=strsplit(In,'0123456789+-',count=ct,length=lenchars)
    if ct eq 0 then return,-1
    posnumbers=poschars[0:ct-1]+lenchars
    if ct eq 1 then lennumbers=strlen(In)-posnumbers $
    else lennumbers=[poschars[1:*],strlen(In)]-posnumbers
    elements=strmid(In,poschars,lenchars)
    repion=strmid(In,posnumbers,lennumbers)
    
    ; Loop over all groups: get Z and check for composites
    Z=0
    El=''
    repion_new=''
    for i=0l,ct-1 do begin
        Eli=strlowcase(elements[i])
        Eli=strupcase(strmid(Eli,0,1))+strmid(Eli,1)
        Zi=where(Mat eq Eli,ct2)
        if ct2 ne 1 then begin
            ; Composite element: e.g. "OH"
            poschars=byte(elements[i])
            poschars=poschars ge 65 and poschars le 90
            poschars[0]=1
            poschars=where(poschars,ct2)
            if ct2 eq 1 then begin
                In=''
                repeats=1
                ion=0
                return,-1
            endif
            lenchars=[poschars[1:*],strlen(Eli)]-poschars
            Eli=strmid(elements[i],poschars,lenchars)
            Zi=intarr(ct2)
            for j=0l,ct2-1 do begin
                Zi[j]=where(Mat eq Eli[j],ct3)
                if ct3 ne 1 then begin
                    In=''
                    repeats=1
                    ion=0
                    return,-1
                endif
            endfor
        endif
        
        El=[El,Eli]
        Z=[Z,Zi]
        if ct2 gt 1 then repion_new=[repion_new,replicate('',ct2-1)]
        repion_new=[repion_new,repion[i]]
    endfor
    
    ct=n_elements(Z)-1
    Z=Z[1:*]
    El=El[1:*]
    repion=repion_new[1:*]
    repeats=intarr(ct)
    ion=repeats
    
    ; Loop over all groups and get charge/repeat: "x+y" becomes +x and y
    In=''
    for i=0l,ct-1 do begin
        tmp=strsplit(repion[i],'+-',count=ct2,/extract)
        case ct2 of
        0:    begin
            if (strpos(repion[i],'-'))[0] ne -1 then ion[i]=-1
            if (strpos(repion[i],'+'))[0] ne -1 then ion[i]=1
            repeats[i]=1
            endcase
        1:     begin
            if strlen(tmp) eq strlen(repion[i]) then begin
                repeats[i]=fix(tmp[0])
                ion[i]=0
            endif else begin
                p0=(strpos(repion[i],'+'))[0]
                m=1
                if p0 eq -1 then begin
                    p0=(strpos(repion[i],'-'))[0]
                    m=-1
                endif
                
                p1=strpos(repion[i],tmp[0])
                if p0 lt p1 then begin
                    repeats[i]=fix(tmp[0])
                    ion[i]=m
                endif else begin
                    repeats[i]=1
                    ion[i]=fix(tmp[0])*m
                endelse
                
            endelse
            endcase
        2:    begin
            repeats[i]=fix(tmp[1])
            ion[i]=fix(tmp[0])
            if (strpos(repion[i],'-'))[0] ne -1 then ion[i]=-ion[i]
            endcase
        else:begin
            In=''
            repeats=1
            ion=0
            return,-1
            endelse
        endcase
        
        In+=El[i]
        if ion[i] ne 0 then In+=string(abs(ion[i]),format='(I0)')+((ion[i] gt 0)?'+':'-')
        if repeats[i] gt 1 then In+=string(repeats[i],format='(I0)')
    endfor

    if In eq '' then return,-1 else return,Z
    
endif else begin
    out=Mat[0>In<nmax]
    
    for i=0l,n_elements(ion)-1 do $
        if ion[i] ne 0 then out[i]+=string(abs(ion[i]),format='(I0)')+((ion[i] gt 0)?'+':'-')
    
    for i=0l,n_elements(repeats)-1 do if repeats[i] gt 1 then out[i]+=string(repeats[i],format='(I0)')
    
    for i=1,n_elements(out)-1 do out[0]+=out[i]
    return,out[0]
endelse

end;function ZElement
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ICSDRadii,Symb;Element symbol, e.g. "Ca2+"

; Data adapted from:
; 
; R.W. Grosse-Kunstleve, N.K. Sauter, N.W. Moriarty and P.D. Adams
; The Computational Crystallography Toolbox: crystallographic 
; algorithms in a reusable software framework
; J. Appl. Cryst. (2002). 35, 126-136 

Elements0=[ $
    "","H","H1+","H1-","D","D1+","D1-","He","Li","Li1+", $
    "Be","Be2+","B","B1+","B2+","B3+","B2-","B3-","C","C1+", $
    "C2+","C3+","C4+","C1-","C2-","C4-","N","N1+","N2+","N3+", $
    "N4+","N5+","N1-","N2-","N3-","O","O1-","O2-","F","F7+", $
    "F1-","Ne","Na","Na1+","Mg","Mg2+","Al","Al3+","Si","Si2+", $
    "Si3+","Si4+","Si1-","Si4-","P","P1+","P2+","P3+","P4+","P5+", $
    "P1-","P2-","P3-","S","S1+","S2+","S3+","S4+","S5+","S6+", $
    "S1-","S2-","Cl","Cl1+","Cl3+","Cl5+","Cl7+","Cl1-","Ar","K", $
    "K1+","Ca","Ca1+","Ca2+","Sc","Sc1+","Sc2+","Sc3+","Ti","Ti2+"]
Elements1=[$
    "Ti3+","Ti4+","V","V1+","V2+","V3+","V4+","V5+","Cr","Cr1+", $
    "Cr2+","Cr3+","Cr4+","Cr5+","Cr6+","Mn","Mn1+","Mn2+","Mn3+","Mn4+", $
    "Mn5+","Mn6+","Mn7+","Mn1-","Fe","Fe1+","Fe2+","Fe3+","Fe4+","Fe6+", $
    "Co","Co1+","Co2+","Co3+","Co4+","Co1-","Ni","Ni1+","Ni2+","Ni3+", $
    "Ni4+","Cu","Cu1+","Cu2+","Cu3+","Zn","Zn2+","Ga","Ga1+","Ga2+", $
    "Ga3+","Ge","Ge2+","Ge3+","Ge4+","Ge4-","As","As2+","As3+","As4+", $
    "As5+","As1-","As2-","As3-","Se","Se1+","Se2+","Se4+","Se6+","Se1-", $
    "Se2-","Br","Br1+","Br3+","Br5+","Br7+","Br1-","Kr","Kr2+","Rb", $
    "Rb1+","Sr","Sr2+","Y","Y1+","Y2+","Y3+","Zr","Zr1+","Zr2+"]
Elements2=[$
    "Zr3+","Zr4+","Nb","Nb1+","Nb2+","Nb3+","Nb4+","Nb5+","Mo","Mo1+", $
    "Mo2+","Mo3+","Mo4+","Mo5+","Mo6+","Tc","Tc2+","Tc3+","Tc4+","Tc6+", $
    "Tc7+","Ru","Ru2+","Ru3+","Ru4+","Ru5+","Ru6+","Ru7+","Rh","Rh1+", $
    "Rh2+","Rh3+","Rh4+","Rh5+","Rh1-","Pd","Pd1+","Pd2+","Pd3+","Pd4+", $
    "Ag","Ag1+","Ag2+","Ag3+","Cd","Cd2+","In","In1+","In2+","In3+", $
    "Sn","Sn2+","Sn3+","Sn4+","Sb","Sb2+","Sb3+","Sb4+","Sb5+","Sb6+", $
    "Sb2-","Sb3-","Te","Te1+","Te2+","Te4+","Te6+","Te1-","Te2-","I", $
    "I1+","I3+","I5+","I7+","I1-","Xe","Xe2+","Xe4+","Xe6+","Xe8+", $
    "Cs","Cs1+","Ba","Ba2+","La","La1+","La2+","La3+","La4+","Ce"]
Elements3=[$
    "Ce2+","Ce3+","Ce4+","Pr","Pr2+","Pr3+","Pr4+","Nd","Nd2+","Nd3+", $
    "Nd4+","Pm","Pm3+","Sm","Sm2+","Sm3+","Eu","Eu2+","Eu3+","Eu4+", $
    "Gd","Gd1+","Gd2+","Gd3+","Gd4+","Tb","Tb1+","Tb2+","Tb3+","Tb4+", $
    "Dy","Dy2+","Dy3+","Ho","Ho2+","Ho3+","Er","Er1+","Er2+","Er3+", $
    "Tm","Tm2+","Tm3+","Yb","Yb2+","Yb3+","Lu","Lu2+","Lu3+","Hf", $
    "Hf2+","Hf3+","Hf4+","Ta","Ta1+","Ta2+","Ta3+","Ta4+","Ta5+","W", $
    "W2+","W3+","W4+","W5+","W6+","Re","Re1+","Re2+","Re3+","Re4+", $
    "Re5+","Re6+","Re7+","Os","Os1+","Os2+","Os4+","Os5+","Os6+","Os7+", $
    "Os8+","Ir","Ir1+","Ir2+","Ir3+","Ir4+","Ir5+","Pt","Pt2+","Pt3+"]
Elements4=[$
    "Pt4+","Pt5+","Pt6+","Au","Au1+","Au2+","Au3+","Au5+","Hg","Hg1+", $
    "Hg2+","Tl","Tl1+","Tl3+","Pb","Pb2+","Pb4+","Bi","Bi1+","Bi2+", $
    "Bi3+","Bi5+","Bi2-","Po","Po2+","Po4+","Po6+","At","At7+","Rn", $
    "Fr","Fr1+","Ra","Ra2+","Ac","Ac3+","Th","Th2+","Th3+","Th4+", $
    "Pa","Pa3+","Pa4+","Pa5+","U","U1+","U2+","U3+","U4+","U5+", $
    "U6+","Np","Np2+","Np3+","Np4+","Np5+","Np6+","Np7+","Pu","Pu2+", $
    "Pu3+","Pu4+","Pu5+","Pu6+","Am","Am2+","Am3+","Am4+","Am5+","Am6+", $
    "Cm","Cm3+","Cm4+","Bk","Bk3+","Bk4+","Cf","Cf3+","Es","Fm", $
    "Md","No","Lr"]
Elements=[ temporary(Elements0) ,$
temporary(Elements1) ,$
temporary(Elements2) ,$
temporary(Elements3) ,$
temporary(Elements4)]

Rads0=[ $
    0.0,0.78,-0.38,1.40,0.78,-0.24,1.40,1.00,1.56,0.59, $
    1.13,0.17,0.95,0.58,0.40,0.02,1.06,1.22,0.86,0.79, $
    0.60,0.55,-0.08,1.10,1.38,1.77,0.80,0.59,0.37,0.16, $
    0.15,-0.12,1.10,1.29,1.48,0.66,0.93,1.21,0.64,0.08, $
    1.15,1.00,1.91,0.97,1.60,0.49,1.43,0.39,1.34,1.25, $
    1.17,0.26,1.41,2.72,1.30,1.01,0.73,0.44,0.40,0.17, $
    1.59,1.20,2.17,1.04,1.26,0.87,0.62,0.37,0.34,0.12, $
    1.44,1.84,1.62,1.30,1.05,0.12,0.20,1.81,1.00,2.34, $
    1.33,1.97,1.70,0.99,1.64,1.36,1.09,0.73,1.45,0.86]
Rads1=[$
    0.67,0.53,1.35,1.02,0.79,0.64,0.59,0.36,1.27,1.07, $
    0.73,0.62,0.44,0.35,0.30,1.32,0.88,0.67,0.58,0.54, $
    0.55,0.27,0.26,1.06,1.27,0.84,0.61,0.49,0.54,0.30, $
    1.26,0.80,0.65,0.52,0.54,1.30,1.24,0.68,0.69,0.60, $
    0.56,1.28,0.46,0.62,0.60,1.39,0.60,1.40,1.14,0.88, $
    0.47,1.40,0.73,0.63,0.40,2.72,1.50,0.52,0.58,0.64, $
    0.33,1.59,1.85,2.11,1.60,1.39,1.08,0.50,0.29,1.77, $
    1.98,1.11,1.06,0.82,0.59,0.39,1.96,1.14,0.74,2.50, $
    1.47,2.15,1.12,1.80,1.11,1.30,0.89,1.60,1.42,1.21]
Rads2=[$
    0.89,0.72,1.48,1.00,0.71,0.70,0.69,0.32,1.40,1.00, $
    0.92,0.67,0.65,0.63,0.42,1.35,1.00,0.80,0.64,0.56, $
    0.98,1.32,0.90,0.68,0.62,0.52,0.37,0.54,1.34,0.82, $
    0.75,0.67,0.62,0.55,1.54,1.37,0.59,0.64,0.76,0.62, $
    1.44,0.67,0.89,0.65,1.57,0.84,1.66,1.35,1.08,0.79, $
    1.58,0.93,0.82,0.69,1.60,0.83,0.76,0.69,0.61,0.75, $
    2.16,2.44,1.70,1.45,1.20,0.52,0.56,1.95,2.21,1.95, $
    1.70,1.39,0.62,0.50,2.20,1.33,1.10,0.83,0.55,0.40, $
    2.71,1.67,2.24,1.34,1.87,1.40,1.27,1.06,1.01,1.82]
Rads3=[$
    1.30,1.03,0.80,1.83,1.00,1.01,0.78,1.82,1.30,1.00, $
    0.90,1.63,0.98,1.80,1.10,0.96,2.04,1.17,0.95,0.65, $
    1.80,0.91,0.94,0.94,1.00,1.78,1.50,1.22,0.92,0.76, $
    1.77,1.10,0.91,1.77,1.10,0.89,1.76,1.50,1.20,0.88, $
    1.75,1.16,0.87,1.94,0.90,0.86,1.73,1.20,0.85,1.59, $
    1.10,0.97,0.71,1.48,0.88,0.83,0.67,0.66,0.64,1.41, $
    0.80,0.75,0.65,0.66,0.41,1.46,1.23,1.00,0.77,0.63, $
    0.52,0.52,0.40,1.34,1.20,1.05,0.63,0.51,0.33,0.27, $
    0.20,1.36,1.37,1.00,0.73,0.63,0.68,1.39,0.80,0.73]
Rads4=[$
    0.63,0.58,0.50,1.44,1.37,1.11,0.70,0.70,1.62,0.97, $
    0.69,1.73,1.47,0.88,1.75,0.94,0.77,1.70,1.45,1.16, $
    0.96,0.74,1.70,1.70,1.40,1.10,0.67,1.53,0.62,1.53, $
    1.53,1.80,1.53,1.43,1.88,1.18,1.80,0.80,0.90,1.00, $
    1.61,1.13,0.98,0.89,1.55,1.40,1.30,1.06,0.97,0.76, $
    0.45,1.58,1.10,1.04,0.95,0.80,0.80,0.71,1.64,0.90, $
    1.00,0.80,0.70,0.60,1.73,1.20,1.01,0.92,0.69,0.50, $
    1.42,0.98,0.95,1.42,0.96,0.93,1.42,0.95,1.42,1.42, $
    1.42,1.42,1.42]
Rads=[ temporary(Rads0) ,$
temporary(Rads1) ,$
temporary(Rads2) ,$
temporary(Rads3) ,$
temporary(Rads4)]

ind=where(Elements eq Symb,ct)
if ct eq 0 then return,0 $
else return,Rads[ind[0]]

end;function ICSDRadii
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AtomColRad,Z,Ion,Repeats

; http://www.crystalmaker.com/support/tutorials/crystalmaker/AtomicRadii.html
;
; He, Ne, Ar, Kr, Rn, Np: Atomic-Ionic Radii
;    Slater JC (1964) Journal of Chemical Physics 39:3199-
; Rest: Shannon & Prewitt Ionic "Crystal" Radii
;    Shannon RD Prewitt CT (1969) Acta Crystallographica B25:925-946
;    Shannon RD (1976) Acta Crystallographica A23:751-761

dR=1.00
R1 = [ dR,$
              0.10 ,0.31, 0.9, 0.41, 0.25, 0.29, 0.30, 1.21, 1.19, 0.38, 1.16,$
              0.86, 0.53, 0.40, 0.31, 0.43, 1.67, 0.71, 1.52, 1.14, 0.89, 0.75,$
              0.68, 0.76, 0.81, 0.69, 0.54, 0.70, 0.71, 0.74, 0.76, 0.53, 0.72,$
              0.56, 1.82, 0.88, 1.66, 1.32, 1.04, 0.86, 0.78, 0.79, 0.79, 0.82,$
              0.81, 0.78, 1.29, 0.92, 0.94, 0.69, 0.90, 1.11, 2.06, 0.62, 1.81]
R2 = [$
              1.49, 1.36, 1.15, 1.32, 1.30, 1.28, 1.10, 1.31, 1.08, 1.18, 1.05,$
              1.04, 1.03, 1.02, 1.13, 1.00, 0.85, 0.78, 0.74, 0.77, 0.77, 0.77,$
              0.74, 1.51, 0.83, 1.03, 1.49, 1.17, 1.08, 0.76, 1.20, 1.94, 1.62,$
              1.26, 1.19, 1.09, 0.87, 1.75, 1.00, 1.12, 1.11, dR, dR, dR,$
              dR, dR, dR, dR,dR]
R=[R1,R2]

dcol=[255,19,147]
Col1 = [ dcol,$
              [255,188,190] ,[255,199,199], [76,136,173], [157,255,152], [133,137,173], [16,6,0], [128,176,255], [255,4,0], [8,119,8], [255,19,147], [201,202,4],$
              [251,201,6], [22,175,255], [0,0,254], [97,97,99], [255,251,2], [21,254,12], [255,19,147], [114,1,255], [95,158,199], [249,107,251], [124,204,255],$
              [229,27,7], [255,9,159], [253,1,150], [117,117,0], [0,2,175], [221,219,222], [0,0,251], [142,143,127], [165,223,120], [127,109,167], [73,41,80],$
              [156,237,18], [120,50,0], [255,19,147], [255,1,149], [0,255,39], [102,151,145], [0,255,1], [78,178,118], [180,134,171], [207,174,201], [209,182,171],$
              [204,208,171], [192,195,186], [183,186,193], [141,129,255], [215,125,189], [154,140,189], [217,133,86], [178,161,83], [85,0,132], [150,160,255], [0,4,253]]
Col2 = [$
             [188,104,117], [96,193,78], [211,251,4], [253,224,6], [253,140,8], [27,127,127], [250,8,120], [248,8,217], [193,3,255], [111,1,255], [48,6,253],$
              [7,67,251], [74,115,59], [40,202,251], [41,253,238], [32,255,181], [181,177,86], [184,156,83], [142,138,129], [178,178,137], [199,177,120], [201,206,114],$
              [194,193,198], [255,175,63], [212,183,201], [149,138,108], [83,83,91], [163,109,187], [127,127,127], [123,123,123], [255,19,147], [128,128,128], [109,171,91],$
              [104,157,113], [34,255,115], [48,249,57], [123,158,176], [255,19,147], [92,94,169], [126,126,126], [127,127,127], dcol, dcol, dcol,$
              dcol, dcol, dcol, dcol,dcol]
Col=[Col1,Col2]

n=n_elements(R)-1
doIon=n_elements(Ion) ne 0
doRepeats=n_elements(Repeats) ne 0

Colout=[0,0,0]
Rout=0.
for i=0l,n_elements(Z)-1 do begin
    RadICSD=abs(ICSDRadii(ZElement(Z[i],Ion=doIon?Ion[i]:0))) ; abs because some are negative

    off=3*(0>Z[i]<n)
    Colout or= Col[off:off+2]
    
    if RadICSD eq 0 then RadICSD=R[0>Z[i]<n]
    if doRepeats then RadICSD*=Repeats[i]
    Rout+=RadICSD
endfor

return,{Col:Colout,R:Rout}
end;function AtomColRad
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AtScatDis,Z,Lamb

;----------------------------------------------------------
; Complex atomic form factor:
;     f  = f0 + f' + fNT + i.f''
; 
; Symbols
;       Etot = total binding energy of atomic electrons (negative number)
;       mc^2 = electron rest mass energy
;       E = photon energy
;       Jensen = Z/2(E/mc^2)^2 = Stibius-Jensen term (Creagh, IT2004)
;       
; Relativistic corrections:
;         frel(CL) = 5/3.Etot/mc^2                    (Cromer Liberman 1970)
;         frel(3/5CL) = Etot/mc^2                     (Chantler 2000)           => this is currently accepted
;         frel(Sasaki) = frel(CL)                     (Sasaki 1989)
;         frel(H82) = - 2.19e-6*Z^3 - 1.03e-4*Z^2     (Henke 1982 w.o. Jensen)
;         frel(H93) = - (Z/82.5)^2.37                 (Henke 1993)
;         frel(IT) = Etot/mc^2                        (IT, Creagh 2004)         => should be the same as frel(3/5CL) but noticed some differences
;       
;----------------------------------------------------------
; FFAST (Chantler, 2000): f1, f2, frel(3/5CL), frel(H82) and fNT
; 
; f = f0 + f1 - Z + frel(3/5CL) + fNT + i.f2
;
; => f' = f1 - Z + frel(3/5CL)
;    f'' = f2
;
;----------------------------------------------------------
; Sasaki (1989): f'(Sasaki) and f''(Sasaki)
;  f  = f0 + f'(Sasaki) + i.f''(Sasaki)
;       f'(Sasaki) = f+ + frel(Sasaki)
;
; => f' = f'(Sasaki) - frel(Sasaki) + frel(3/5CL)
;       = f'(Sasaki) - 2/3.frel(3/5CL)
;    f''= f''(Sasaki)
;
;----------------------------------------------------------
; Henke (1993): f1(0) and f2(0)
;  f = f1 + i.f2
;      f1 = f1(0) - Z + f0
;      f1(0) = Z* + f+
;              Z* = Z + frel(Henke)
;              frel(H93) = Etot/mc^2 ~= - (Z/82.5)^2.37    (fit to data from Kissel & Pratt (1990))
;      f2 = f2(0)
;
; => f' = f1(0) - Z - frel(H93) + frel(3/5CL)
;       = f1(0) - Z + (Z/82.5)^2.37 + frel(3/5CL)
;    f''= f2(0)
;
;----------------------------------------------------------
; Henke (1982): not used in XRDUA
;     Z* = Z + frel(H82) - Jensen
;          frel(H82) = 5/3Etot/mc^2 ~= - 2.19e-6*Z^3 - 1.03e-4*Z^2   (fit to frel(CL))
;
;----------------------------------------------------------
; Cu(12):
;     frel(CL)    = -0.146
;     frel(3/5CL) = -0.0876     (=3/5CL)
;     frel(H82)   = -0.14003    (idem to formula)
;     frel(H93)   = -0.0839242  (only formula)
;     frel(IT)    = -0.0878     (different from frel(3/5CL))
;     
; Hg(80):
;     frel(CL)    = -1.743
;     frel(3/5CL) = -1.04580    (=3/5CL)
;     frel(H82)   = -1.78050    (DIFFERENT from formula !!??)
;     frel(H93)   = -0.929667   (only formula)
;     frel(IT)    = -1.046      (different from frel(3/5CL))
;
;----------------------------------------------------------
; Henke tables: experimental
; Sasaki tables: Cromer and Liberman calculations with precission near absorption edges
; FFAST tables (Chantler): ???

fan=fltarr(2)
Zstr=stringr(Z)
fcor=FFAST_fcor(Zstr)

bSas=SasakiUse(Z,Lamb) ; Use Sasaki when energy around an edge
if bSas then begin
    fan=Sasaki(Zstr,Lamb) ; binterpol not needed because checked already with SasakiUse
    fan[0]+=fcor[1]-2/3.*fcor[0]
endif else begin
    fan=Henke(Zstr,LambdaEnergy(Lamb),binterpol=binterpol)
    if binterpol then begin
        fan[0]+=fcor[1]-Z+(Z/82.5)^2.37+fcor[0]
    endif else begin
        fan=FFAST(Zstr,2,LambdaEnergy(Lamb),binterpol=binterpol)
        if n_elements(fan) ne 2 then fan=fltarr(2)
        if binterpol then fan[0]+=fcor[1]-Z+fcor[0] else begin
            fan[0]=fcor[0]+fcor[1]
            fan[1]=0
        endelse
    endelse
endelse

; [f'+fNT(e/atom),f''(e/atom)]
return,fan
end;function AtScatDis
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AtScat,Z,Ion,Repeats,d,anomal=anomal,Lamb=Lamb

;cut=total((d gt [1./12,0.25])*[2.,1.])

f0=0.
fan=fltarr(2)
for i=0l,n_elements(Z)-1 do begin
    Zsymb=ZElement(Z[i],Ion=Ion[i])
    
    ;case cut of
    ;3: f0+=CromerMannf0(Zsymb,d)*Repeats[i] ; 1/4 < d
    ;2: f0+=WaasKirff0(Zsymb,d)*Repeats[i]   ; 1/12 < d
    ;0: f0+=WaasKirff0(Zsymb,d)*Repeats[i] ; 0 < d ; Interpolation of Kissel data not implemented yet: use WaasKirff0
    ;endcase
    
    add=WaasKirff0(Zsymb,d)
    if add eq 0 then add=CromerMannf0(Zsymb,d)
    
    if add eq 0 then begin
        Zsymb=ZElement(Z[i])
        add=WaasKirff0(Zsymb,d)
        if add eq 0 then add=CromerMannf0(Zsymb,d)
    endif
    
    f0+=add*Repeats[i]
    
    ; fan=[f',f'']
    if keyword_set(anomal) then fan+=AtScatDis(Z[i],Lamb)*Repeats[i]
endfor

return,complex(f0+fan[0],fan[1])
; Re(f) = f0+f'
; Im(f) = f''
end;function AtScat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AtMassUse,Z,Ion,Repeats
; Ion is not used!!

ret=0.
for i=0l,n_elements(Z)-1 do ret+=AtMass(Z[i])*Repeats[i]

return,ret
end;function AtMassUse
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MassAtCoef,Z,Ion,Repeats,lambda,mmatom=mmatom
; Ion is not used!!

ret=0.
En=LambdaEnergy(lambda)
mmatom=0.
for i=0l,n_elements(Z)-1 do begin
    mmatomi=AtMass(Z[i])*Repeats[i]
    ret+=mmatomi*(FFAST(strtrim(Z[i],2),1,En))[0]
    mmatom+=mmatomi
endfor
return,ret/mmatom

; Part of the mass attenutation coefficient is the
; Mass photoelectric attenuation coefficient:
; mu(cm�/g)= Na(1/mol).�a(cm�) = Na(1/mol).2.re(m).lambda(Angstrom).10^(-6).f2
;           -----------------   --------------------------------------
;              MM(g/mol)            MM(g/mol)
; (cfr. http://physics.nist.gov/PhysRefData/FFast/Text2000/sec02.html)
;
;fan=AtScatDis(Z,lambda)
;mu=avogadroconstant()*2.*electronradius()*10^(-6.)*lambda*fan[1]/AtMass(Z)

end;function MassAtCoef
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CohCrossSection,Z,Ion,Repeats,lambda,mmatom=mmatom
; Ion is not used!!

ret=0.
En=LambdaEnergy(lambda)
mmatom=0.
for i=0l,n_elements(Z)-1 do begin
    mmatomi=AtMass(Z[i])*Repeats[i]
    ret+=mmatomi*XCOM(strtrim(Z[i],2),1,En)
    mmatom+=mmatomi
endfor
return,ret/mmatom

end;function CohCrossSection
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SourceEmission,context

; Data adapted from:
; 
; Coehlo, A. A. (2008). Topas-Academic, 
; http://members.optusnet.com.au/alancoelho/
; 
; Kraus, W.; Nolze, G..  Powder Cell - a program 
; for the representation and manipulation of crystal 
; structures and calculation of the resulting 
; x-ray powder patterns.
; Journal of Applied Crystallography  (1996),  29(3),  301-303

t = widget_button(context,value='Cr-Ka(2)',uvalue=[[2.28970,1,1],$
                        [2.293606,0.515,1]])
t = widget_button(context,value='Cr-Ka(7)',uvalue=[[2.2896977,0.378,0.6160974],$
                        [2.2936467,0.271,0.9436726],$
                        [2.2900254,0.132,0.7444350],$
                        [2.2905983,0.084,1.3279580],$
                        [2.2915136,0.073,2.1807273],$
                        [2.2943110,0.054,2.0124101],$
                        [2.2882482,0.009,0.8395685]])
t = widget_button(context,value='Cr-Kb(5)',uvalue=[[2.0848200,0.307,0.5959633],$
                        [2.0889262,0.236,5.6241545],$
                        [2.0850864,0.172,0.6662470],$
                        [2.0865602,0.148,2.3492088],$
                        [2.0855459,0.137,1.1822326]])
t = widget_button(context,value='Mn-Ka(2)',uvalue=[[2.101820,1,1],$
                        [2.10578,0.519,1]])
t = widget_button(context,value='Mn-Ka(7)',uvalue=[[2.1018365,0.353,0.6110764],$
                        [2.1058026,0.229,0.8444323],$
                        [2.1021879,0.141,0.7281904],$
                        [2.1062490,0.110,1.5085288],$
                        [2.1032712,0.079,1.6052406],$
                        [2.1026638,0.066,0.9496080],$
                        [2.1016355,0.005,0.3452010]])
t = widget_button(context,value='Mn-Kb(5)',uvalue=[[1.9101270,0.254,0.5385290],$
                        [1.9114758,0.234,2.7701238],$
                        [1.9140076,0.234,3.9061842],$
                        [1.9103713,0.164,0.5327797],$
                        [1.9107334,0.114,0.8274467]])
t = widget_button(context,value='Fe-Ka(2)',uvalue=[[1.936042,1,1],$
                        [1.939980,0.500,1]])
t = widget_button(context,value='Fe-Ka(7)',uvalue=[[1.9359990,0.278,0.4876162],$
                        [1.9399242,0.207,0.7548816],$
                        [1.9362569,0.182,0.5941855],$
                        [1.9370562,0.106,1.4626310],$
                        [1.9366253,0.094,0.8479062],$
                        [1.9405570,0.066,0.7104222],$
                        [1.9402020,0.065,1.3459384]])
t = widget_button(context,value='Fe-Kb(4)',uvalue=[[1.7594154,0.301,3.5378594],$
                        [1.7568450,0.279,0.7767031],$
                        [1.7565588,0.241,0.4902585],$
                        [1.7574577,0.179,1.5893663]])
t = widget_button(context,value='Co-Ka(2)',uvalue=[[1.788965,1,1],$
                        [1.792850,0.497,1]])
t = widget_button(context,value='Co-Ka(3)',uvalue=[[1.788960,0.5981,0.4700],$
                        [1.792850,0.3516,0.8052],$
                        [1.789879,0.0503,0.3357]])
t = widget_button(context,value='Co-Ka(7)',uvalue=[[1.7889847,0.378,0.4633522],$
                        [1.7927905,0.197,0.6237179],$
                        [1.7892524,0.144,0.6958819],$
                        [1.7896946,0.127,1.1767380],$
                        [1.7930637,0.095,0.7190761],$
                        [1.7888515,0.088,0.2085420],$
                        [1.7934738,0.050,1.1578452]])
t = widget_button(context,value='Co-Kb(6)',uvalue=[[1.6207938,0.449,0.6462327],$
                        [1.6211689,0.189,0.7588800],$
                        [1.6228580,0.153,2.0774644],$
                        [1.6216651,0.103,1.0372054],$
                        [1.6236359,0.082,2.8895524],$
                        [1.6198346,0.025,0.8020733]])
t = widget_button(context,value='Ni-Ka(2)',uvalue=[[1.657910,1,1],$
                        [1.661747,0.495,1]])
t = widget_button(context,value='Ni-Ka(5)',uvalue=[[1.6579244,0.487,0.4462793],$
                        [1.6617352,0.250,0.5955505],$
                        [1.6583129,0.171,1.0449117],$
                        [1.6620153,0.064,0.6770710],$
                        [1.6624264,0.028,0.9977195]])
t = widget_button(context,value='Ni-Kb(4)',uvalue=[[1.5001100,0.450,0.6824449],$
                        [1.5004731,0.258,0.7880970],$
                        [1.5016253,0.203,2.4915953],$
                        [1.4994406,0.089,0.9393379]])
t = widget_button(context,value='Cu-Ka(2)',uvalue=[[1.540598,0.66050,0.5],$
                        [1.544426,0.33950,0.5]])
t = widget_button(context,value='Cu-Ka(4)',uvalue=[[1.5405909,0.579,0.4374157],$
                        [1.5443990,0.236,0.5128764],$
                        [1.5446855,0.105,0.6872322],$
                        [1.5410639,0.080,0.6432140]])
t = widget_button(context,value='Cu-Ka(5)',uvalue=[[1.540596,0.5791,0.437],$
                        [1.54441,0.2417,0.52],$
                        [1.544721,0.0871,0.62],$
                        [1.541058,0.0762,0.6],$
                        [1.534753,0.0159,3.6854]])
t = widget_button(context,value='Cu-Kb(4)',uvalue=[[1.3922160,0.485,0.5502872],$
                        [1.3925949,0.248,0.5505868],$
                        [1.3917581,0.110,0.5546122],$
                        [1.3934905,0.100,1.2654733],$
                        [1.3913004,0.055,0.8290293]])
t = widget_button(context,value='Zr-Ka(2)',uvalue=[[0.78593,1,1],$
                        [0.79015,0.502,1]])
t = widget_button(context,value='Nb-Ka(2)',uvalue=[[0.74620,1,1],$
                        [0.75044,0.498,1]])
t = widget_button(context,value='Mo-Ka(2)',uvalue=[[0.709300,0.6533,0.2695],$
                        [0.713574,0.3467,0.2795]])
t = widget_button(context,value='Ag-Ka(2)',uvalue=[[0.5594075,1,1],$
                        [0.563798,0.499,1]])
t = widget_button(context,value='W-La(2)',uvalue=[[1.4767,1,1],$
                        [1.4877,0.115,1]])

end;pro SourceEmission
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


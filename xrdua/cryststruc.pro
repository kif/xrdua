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

function AtomPlotObject::Init, POS=pos, RADIUS=radius, DENSITY=density, $
                       TEX_COORDS=tex_coords, POBJ=pobj, _EXTRA=e

    ; Set properties of superclass
    IF (self->IDLgrModel::Init(_EXTRA=e) NE 1) THEN RETURN, 0
    
    ; Initialize properties
    self.pos = [0.0,0.0,0.0]
    self.radius = 1.0
    self.density = 1.0
    
    ; Set properties
    bcopy=0b
    if n_elements(e) ne 0 then begin
        t=tag_names(e)
        ind=where(t eq 'POLYGONS',ct)
        bcopy=(ct eq 1) and (n_elements(pobj) ne 0)
    endif
    
    IF (N_ELEMENTS(pos) EQ 3) THEN $
        self.pos = pos
    IF (N_ELEMENTS(radius) EQ 1) THEN $
        self.radius = radius
    IF (N_ELEMENTS(density) EQ 1) THEN $
        self.density = density
    IF (N_ELEMENTS(tex_coords) EQ 1) THEN $
        self.texture = tex_coords

    ; Set the polygon vertices and connectivity based on property settings.
    if bcopy then begin
        self.oPoly = OBJ_NEW('IDLgrPolygon', SHARE_DATA=pobj, /REJECT, _EXTRA=e)
        self->Add,self.oPoly
    endif else begin
        self.oPoly = OBJ_NEW('IDLgrPolygon', /REJECT, _EXTRA=e)
        self->Add,self.oPoly
        self->BuildPoly
    endelse
    
    RETURN, 1
end;function AtomPlotObject::Init
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AtomPlotObject::Cleanup

    ; Cleanup the polygon object used to represent the atom.
    OBJ_DESTROY, self.oPoly

    ; Cleanup the superclass.
    self->IDLgrModel::Cleanup
end;pro AtomPlotObject::Cleanup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AtomPlotObject::SetProperty, POS=pos, RADIUS=radius, DENSITY=density, $
                       TEX_COORDS=tex_coords, _EXTRA=e

    ; Pass along extraneous keywords to the superclass and/or to the
    ; polygon used to represent the atom.
    self->IDLgrModel::SetProperty, _EXTRA=e
    self.oPoly->SetProperty, _EXTRA=e
    
    IF (N_ELEMENTS(pos) EQ 3) THEN $
        self.pos = pos
    IF (N_ELEMENTS(radius) EQ 1) THEN $
        self.radius = radius
    IF (N_ELEMENTS(density) EQ 1) THEN $
        self.density = density
    IF (N_ELEMENTS(tex_coords) EQ 1) THEN $
        self.texture = tex_coords

    ; Rebuild the polygon according to keyword settings.
    self->BuildPoly
    
end;pro AtomPlotObject::SetProperty
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AtomPlotObject::GetProperty, POS=pos, RADIUS=radius, DENSITY=density,$
                         TEX_COORDS=tex_coords, POBJ=pobj, ALL=all, _REF_EXTRA=re
    
    ball= ARG_PRESENT(ALL)
    if ball then begin
        self->IDLgrModel::GetProperty, ALL=all, _EXTRA=re
        self.oPoly->GetProperty, ALL=all, _EXTRA=re
    endif else begin
        self->IDLgrModel::GetProperty, _EXTRA=re
        self.oPoly->GetProperty, _EXTRA=re
    endelse
    
    pos = self.pos
    radius = self.radius
    density = self.density
    pobj = self.oPoly
    tex_coords = self.texture
    if ball then all=create_struct('POS', pos, 'RADIUS', radius, 'DENSITY', density,$
                         'TEX_COORDS', tex_coords, all)
        
end;pro AtomPlotObject::GetProperty
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AtomPlotObject::Copy

    self->GetProperty, POBJ=pobj, ALL=e
    oCopy = OBJ_NEW('AtomPlotObject', POBJ=pobj ,_EXTRA=e)
    
    return,oCopy
end;function AtomPlotObject::Copy
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AtomPlotObject::Print
    PRINT, 'Position: ',self.pos
    PRINT, 'Radius: ',self.radius
    PRINT, 'Density: ',self.density
end;pro AtomPlotObject::Print
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AtomPlotObject::BuildPoly
    ; Build the atom.

    ; Number of rows and columns of vertices is based upon the density
    ; property.
    nrows = LONG(20.0*self.density) > 2
    ncols = LONG(20.0*self.density) > 2

    ; Create the vertex list and the connectivity array.
    nverts = nrows*ncols + 2
    nconn = (ncols*(nrows-1)*5) + (2*ncols*4)
    conn = LONARR(ncols*(nrows-1)*5 + 2*ncols*4)
    verts = FLTARR(3, nverts)
    IF (self.texture NE 0) THEN $
        tex = FLTARR(2,nverts)

    ; Fill in the vertices.
    i = 0L
    j = 0L
    k = 0L
    tzinc = (!PI)/FLOAT(nrows+1)
    tz = (!PI/2.0) - tzinc
    for k=0l,(nrows-1) DO BEGIN
        z = SIN(tz)*self.radius
        r = COS(tz)*self.radius
        t = 0
        tinc = (2.*!PI)/FLOAT(ncols-(self.texture NE 0))

        for j=0l,(ncols-1) DO BEGIN
            verts[0,i] = r*COS(t) + self.pos[0]
            verts[1,i] = r*SIN(t) + self.pos[1]
            verts[2,i] = z + self.pos[2]

            IF (self.texture NE 0) THEN BEGIN
                tex[0,i] = t/(2.*!PI)
                tex[1,i] = (tz+(!PI/2.0))/(!PI)
            endIF

            t += tinc
            i ++
        endFOR
        tz -= tzinc
    endFOR
    top = i
    verts[0,i] = self.pos[0]
    verts[1,i] = self.pos[1]
    verts[2,i] = self.radius*1.0 + self.pos[2]
    i++
    bot = i
    verts[0,i] = self.pos[0]
    verts[1,i] = self.pos[1]
    verts[2,i] = self.radius*(-1.0) + self.pos[2]

    IF (self.texture NE 0) THEN BEGIN
        tex[0,i] = 0.5
        tex[1,i] = 0.0
        tex[0,i-1] = 0.5
        tex[1,i-1] = 1.0
    endIF

    ; Fill in the connectivity array.
    i = 0L
    for k=0l,(nrows-2) DO BEGIN
        for j=0l,(ncols-1) DO BEGIN
            conn[i] = 4

            conn[i+4] = k*ncols + j

            w = k*ncols + j + 1L
            IF (j EQ (ncols-1)) THEN w = k*ncols
            conn[i+3] = w

            w = k*ncols + j + 1L + ncols
            IF (j EQ (ncols-1)) THEN w = k*ncols + ncols
            conn[i+2] = w

            conn[i+1] = k*ncols + j + ncols

            i += 5
            IF ((self.texture NE 0) AND (j EQ (ncols-1))) THEN $
                i -= 5
        endFOR
    endFOR
    for j=0l,(ncols-1) DO BEGIN
        conn[i] = 3
        conn[i+3] = top
        conn[i+2] = j+1L
        IF (j EQ (ncols-1)) THEN conn[i+2] = 0
        conn[i+1] = j
        i += 4
        IF ((self.texture NE 0) AND (j EQ (ncols-1))) THEN $
            i -= 4
    endFOR
    for j=0l,(ncols-1) DO BEGIN
        conn[i] = 3
        conn[i+3] = bot
        conn[i+2] = j+(nrows-1L)*ncols
        conn[i+1] = j+(nrows-1L)*ncols+1L
        IF (j EQ (ncols-1)) THEN conn[i+1] = (nrows-1L)*ncols
        i += 4
        IF ((self.texture NE 0) AND (j EQ (ncols-1))) THEN $
            i -= 4
    endFOR

    self.oPoly->SetProperty, DATA=verts, POLYGONS=conn

    IF (self.texture NE 0) THEN $
        self.oPoly->SetProperty, TEXTURE_COORD=tex
end;pro AtomPlotObject::BuildPoly
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AtomPlotObject__define
; From IDL's orb object
    struct = { AtomPlotObject, $
               INHERITS IDLgrModel, $
               pos: [0.0,0.0,0.0], $
               radius: 1.0, $
               density: 1.0, $
               texture: 0, $
               oPoly: OBJ_NEW() $
             }
end;pro AtomPlotObject__define
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function getptrFit,ptr
ptrFit=ptr
while (where(tag_names(*ptrFit) eq 'PREV'))[0] ne -1 do ptrFit=(*ptrFit).prev

return,ptrFit
end;function getptrFit
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indCelParam,ptr
return,(*(*ptr).paramind[0,1])[1:*]
end;function indCelParam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function indZeroParam,ptr
return,(*(*ptr).paramind[0,1])[0]
end;function indZeroParam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ZOCompress0,xyz
xyz mod=1
xyz += xyz lt 0
end;pro ZOCompress0
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ZOCompress1,xyz
ind=where(xyz lt 0 or xyz gt 1,ct)
if (ct ne 0) then xyz[ind]-=floor(xyz[ind])
end;pro ZOCompress1
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dfromhkl,hkl,celparam
n=n_elements(hkl)/3
M=MetricTensor(celparam,/reciprocal)
M=rebin(M[[0,1,2,4,5,8]],6,n,/sample)
M*=[hkl[0,*]*hkl[0,*],2*hkl[0,*]*hkl[1,*],2*hkl[0,*]*hkl[2,*],$
    hkl[1,*]*hkl[1,*],2*hkl[1,*]*hkl[2,*],hkl[2,*]*hkl[2,*]]
return,1./sqrt(total(M,1))
end;function dfromhkl
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ttfromhkl,hkl,celparam,lambda,dist,zero,zerotype
d=dfromhkl(hkl,celparam)
tt=BraggDtoX(d,lambda,/onlyt,/angledeg)

case zerotype of
0:    tt+=zero ; Zero: zero shift
1:    begin
    ; Zero: delta sample-detector
    if zero ne 0 then tt=180./!pi*atan((1+zero/(dist*1000.))*tan(tt/180.*!pi))
    endcase
else:
endcase

return,tt
end;function ttfromhkl
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ttpdercell,hkl,param,lambda,pder=pder,bpder=bpder,connect=connect

; tt = 2*asin(lambda/2*sqrt(S1))*180/pi
;
; S1 = X^2 + Y^2 + Z^2 + 2XY cos(gamma') + 2XZ cos(beta') + 2YZ cos(alpha')
;
; X = hbc sin(alpha)/V       Y = kac sin(beta)/V       Z = lab sin(gamma)/V
;   = h sin(alpha)/a sqrt(S2)  = k sin(beta)/b sqrt(S2)  = l sin(gamma)/c sqrt(S2)
;
; V = abc sqrt(S2)
; S2 = 1+2cos(alpha)cos(beta)cos(gamma)-cos(alpha)^2-cos(beta)^2-cos(gamma)^2
;
; cos(alpha') = (cos(beta)cos(gamma) - cos(alpha))/(sin(beta)sin(gamma))
; cos(beta') = (cos(alpha)cos(gamma) - cos(beta))/(sin(alpha)sin(gamma))
; cos(gamma') = (cos(alpha)cos(beta) - cos(gamma))/(sin(alpha)sin(beta))
;
; dtt/dp = M*0.5*dS1/Dp
;        = M1*dX/dp + M2*dY/dp + M3*dZ/dp + M4*dcos(alpha')/dp + M5*dcos(beta')/dp + M6*dcos(gamma')/dp
;
; M = lambda*180/pi / sqrt( S1-lambda^2*S1^2 / 4)
; M1 = M*(X+Ycos(gamma')+Zcos(beta'))
; M2 = M*(Y+Xcos(gamma')+Zcos(alpha'))
; M3 = M*(Z+Xcos(beta')+Ycos(alpha'))
; M4 = M*YZ
; M5 = M*XZ
; M6 = M*XY
;
; Depends on a: X
; Depends on b: Y
; Depends on c: Z
; Depends on alpha: all
; Depends on beta: all
; Depends on gamma: all

; ----unit cell parameters----
ct=!dtor
aa=param[0]
bb=param[1]
cc=param[2]
p1=param[3]*ct
p2=param[4]*ct
p3=param[5]*ct

; ----derived parameters(1)----
cp1=cos(p1)
cp2=cos(p2)
cp3=cos(p3)
sp1=sin(p1)
sp2=sin(p2)
sp3=sin(p3)
cq1=(cp2*cp3-cp1)/(sp2*sp3)
cq2=(cp3*cp1-cp2)/(sp3*sp1)
cq3=(cp1*cp2-cp3)/(sp1*sp2)

; ----derived parameters(2)----
S2=1.+2*cp1*cp2*cp3-cp1*cp1-cp2*cp2-cp3*cp3
SS2=sqrt(S2)

X=sp1/aa*hkl[0,*]/SS2
Y=sp2/bb*hkl[1,*]/SS2
Z=sp3/cc*hkl[2,*]/SS2

S1=X*X+Y*Y+Z*Z+2*X*Y*cq3+2*X*Z*cq2+2*Y*Z*cq1

; ----function----
F=360./!pi*asin(0.5*lambda*sqrt(S1))
if not keyword_set(bpder) then return,F
nhkl=n_elements(F)

; ----derived parameters(3)----
M=lambda*180./!pi/sqrt(S1-lambda*lambda/4.*S1*S1)
M1=M*(X+Y*cq3+Z*cq2)
M2=M*(Y+X*cq3+Z*cq1)
M3=M*(Z+X*cq2+Y*cq1)
M4=M*Y*Z
M5=M*X*Z
M6=M*X*Y
X*=M1 ; X*M1
Y*=M2 ; Y*M2
Z*=M3 ; Z*M3

; ----partial derivatives----
if n_elements(connect) eq 0 then connect=indgen(6)

; a,b,c,alpha,beta,gamma
pder=dblarr(nhkl,6)

; dtt_H/da
pder[*,connect[0]]-=X/aa
; dtt_H/db
pder[*,connect[1]]-=Y/bb
; dtt_H/dc
pder[*,connect[2]]-=Z/cc
; dtt_H/dalpha
pder[*,connect[3]]+=cp1/aa*hkl[0,*]/SS2*M1+(X+Y+Z)*(sp1*cp2*cp3-cp1*sp1)/S2
pder[*,connect[3]]+=M4*(sp1/sp2/sp3)
pder[*,connect[3]]-=M5*(cp3/sp3+(cp1*cp3-cp2)*cp1/sp1/sp1/sp3)
pder[*,connect[3]]-=M6*(cp2/sp2+(cp1*cp2-cp3)*cp1/sp1/sp1/sp2)
; dtt_H/dbeta
pder[*,connect[4]]+=cp2/bb*hkl[1,*]/SS2*M2+(X+Y+Z)*(cp1*sp2*cp3-cp2*sp2)/S2
pder[*,connect[4]]-=M4*(cp3/sp3+(cp2*cp3-cp1)*cp2/sp2/sp2/sp3)
pder[*,connect[4]]+=M5*(sp2/sp1/sp3)
pder[*,connect[4]]-=M6*(cp1/sp1+(cp2*cp1-cp3)*cp2/sp2/sp2/sp1)
; dtt_H/dgamma
pder[*,connect[5]]+=cp3/cc*hkl[2,*]/SS2*M3+(X+Y+Z)*(cp1*cp2*sp3-cp3*sp3)/S2
pder[*,connect[5]]-=M4*(cp2/sp2+(cp3*cp2-cp1)*cp3/sp3/sp3/sp2)
pder[*,connect[5]]-=M5*(cp1/sp1+(cp3*cp1-cp2)*cp3/sp3/sp3/sp1)
pder[*,connect[5]]+=M6*sp3/sp1/sp2

pder[*,3]*=ct
pder[*,4]*=ct
pder[*,5]*=ct

return,F
end;function ttpdercell
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GenHKL,param,dlim,sg
; ----hkl list----
; 1) Restrict the reciprocal lattice points hkl to
;     the once in a box around the sphere with the maximal reciprocal
;     lattice vector as radius i.e. radius=1/dlim[0]
;     => hmax,kmax,lmax
; 2) Restrict hkl to min and max sphere and exclude origin
;     => slim
; 3) Restrict by checking spacegroup systematic extinctions
;     => SysAbs
; 4) Restrict hkl by checking whether it is in the standard reciprocal space ASU
;     => ASUfn

paramrecip=DirectToReciprocalU(param)

tmp=ceil(1/(dlim[0]*paramrecip[0:2]))
;tmp=ceil(param[0:2]/dlim[0]) ; this might be enough?
hmax=tmp[0]
kmax=tmp[1]
lmax=tmp[2]
slim=1./dlim/dlim
hkl=[0,0,0]
d=[0.]
mult=[0]

M=MetricTensor(param,/reciprocal)
M=M[[0,1,2,4,5,8]]

lg=pgdata(sg.lghash,0)
if lg.pghash eq 0 then return, {hkl:0,d:0,mult:0,n:0}

;nn=0l
;mm=0l
for h=-hmax,hmax do begin ; 1)
    for k=-kmax,kmax do begin ; 1)
        for l=-lmax,lmax do begin ; 1)
            ; Sort the following 3 checks in accending order of processing time
             tmp=total([h*h,2*h*k,2*h*l,k*k,2*k*l,l*l]*M) ; H^(T).M^(-1).H
             if (tmp le slim[0]) and (tmp ge slim[1]) and (tmp ne 0) then begin; 2)
;                 nn++
                 if ASUfn(lg.ASU,h,k,l) then begin ; 4)
;                     mm++
                     if ~SysAbs(h,k,l,sg.allops,mult=multadd) then begin ; 3)
                         hkl=[[hkl],[h,k,l]]
                         d=[d,1./sqrt(tmp)]
                         mult=[mult,multadd]
                     endif
                 endif
             endif
         endfor
    endfor
endfor

if n_elements(hkl) eq 3 then return, {hkl:0,d:0,mult:0,n:0}
hkl=hkl[*,1:*]
d=d[1:*]
mult=mult[1:*]
n=n_elements(d)

;print,'HKL from 2theta limits:',nn
;print,'HKL after ASU:',mm
;print,'HKL after sysabs:',n
;print,'HKL with multiplicity:',total(mult)

sorti=reverse(sort(d))

return,{hkl:hkl[*,sorti],d:d[sorti],mult:mult[sorti],n:n}
end;function GenHKL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WyckLabel,Wyck
n=n_elements(Wyck)
ret=strarr(n)
for i=0l,n-1 do ret[i]=stringr(Wyck[i].M)+Wyck[i].L+' ('+StringRTop(Symop64(Wyck[i].Rep))+')'
return,ret
end;function WyckLabel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakeASUpos,natoms,name,Z,Ion,Repeats,xyz,SOF,B

return,{    natoms:long(natoms),$            ; number of atoms on this position in the ASU
            name:name,$                     ; name for this position in the asymmetric unit
            
            Z:ptr_new(fix(Z)),$                ; atomnumbers for this position in the asymmetric unit
            Ion:ptr_new(fix(Ion)),$            ; ion charges for this position in the asymmetric unit
            Repeats:ptr_new(fix(Repeats)),$    ; element repeat count
            
            Wyck:{L:'',M:0,Rep:CodeType64(0,0),xyz:float(xyz),togglefix:replicate(1b,3)},$
                                            ; Wyckoff position letter
                                            ; Wyckoff position multiplicity
                                            ; Wyckoff representative: e.g. (2x,0,0)
                                            ; Wyckoff xyzfit: e.g. (0.1,0,0)
                                            ; Wyckoff togglefix: fixed parameter
            xyz:float(xyz),$                 ; fractional coordinates for this position
            SOF:float(SOF),$                ; there multiplied substitution and occupation factor
            B:float(B),$                    ; there isotropic temperature factors
            togglefix:replicate(1b,2)}         ; fix SOF,B ?
end;function MakeASUpos
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddASUpos,ASU,name,Z,Ion,repeats,xyz,SOF,B

bfirst=1b
bundefined=n_elements(name) eq 0

; Are there already positions in ASU?
if n_elements(ASU) ne 0 then begin
    if bundefined then return
    if ASU.nasuatoms ne 0 then begin
        bfirst=0b
        nasupos=ASU.nasupos
        nasuatoms=ASU.nasuatoms

        if n_elements(xyz) eq 0 then xyz=ASU.asupos.xyz[ASU.nasuatoms-1] $ ; substitute last position
        else ZOCompress0,xyz

        ; New of substitution?
        ind=where(total(ASU.asupos.xyz eq rebin(xyz,3,ASU.nasuatoms),1) eq 3,ct)
        if ct ne 0 then begin ; substitution
            ASU.asupos[ind[0]].natoms++
            if bundefined then name=ASU.asupos[ind[0]].name

            asupos=[ASU.asupos[0:ind[ct-1]],MakeASUpos(-1L,name,Z,Ion,repeats,xyz,SOF,B)]
            if ind[ct-1] ne ASU.nasuatoms-1 then $
                asupos=[asupos,ASU.asupos[ind[ct-1]+1:*]]
        endif else begin ; new position
            nasupos++
            asupos=[ASU.asupos,MakeASUpos(1L,name,Z,Ion,repeats,xyz,SOF,B)]
        endelse
        nasuatoms++

    endif
endif

; First position or create empty structure
if bfirst then begin
    if bundefined then begin
        nasupos=0L
        nasuatoms=0L
        asupos=0L
    endif else begin
        nasupos=1L
        nasuatoms=1L
        ZOCompress0,xyz
        asupos=MakeASUpos(1L,name,Z,Ion,repeats,xyz,SOF,B)
    endelse
endif

; Calculate reverence indices
if nasuatoms ne 0 then begin
    reverse_indices=where(asupos.natoms ne -1,ct)
    reverse_indices=[reverse_indices,nasuatoms]
endif else reverse_indices=0L

; Make new ASU
ASU={nasupos:nasupos,$         ; number of positions in the asymmetric unit
    nasuatoms:nasuatoms,$    ; number of atoms in the asymmetric unit
    asupos:asupos,$            ; array of structures describing ASU positions
    reverse_indices:reverse_indices}
                            ; atoms on position i (i goes from 0 to nasupos-1)
                            ;     print,asu.asupos[asu.reverse_indices[i]:asu.reverse_indices[i+1]-1]
                            ; all ASU fractional coordinates
                            ;    print,asu.asupos[asu.reverse_indices[0:asu.nasupos-1]].xyz
                            ; atom j sits on position i
                            ;    i=value_locate(asu.reverse_indices,j)
end;pro AddASUpos
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ResetWyckASU,sg,asu,SmallDev

if asu.nasupos eq 0 then return

asupos=asu.asupos[asu.reverse_indices[0:asu.nasupos-1]]
SetASUwyckoff,sg,asupos,SmallDev

for i=0l,asu.nasupos-1 do $
    for j=asu.reverse_indices[i],asu.reverse_indices[i+1]-1 do asu.asupos[j].wyck=asupos[i].wyck

end;pro ResetWyckASU
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ResetWyckpASU,sg,ptr,SmallDev

if (*ptr).nasupos eq 0 then return
p=(*ptr).asu.I
prev=(*ptr).preverse_indices

asupos=(*p)[(*prev)[0:(*ptr).nasupos-1]]
SetASUwyckoff,sg,asupos,SmallDev

for i=0l,(*ptr).nasupos-1 do $
    for j=(*prev)[i],(*prev)[i+1]-1 do (*p)[j].wyck=asupos[i].wyck

end;pro ResetWyckpASU
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro KeepUniXYZ,xyz,row,SmallDev

; ----Keep only the unique positions----
xyznew=xyz[*,0]
rownew=1
n=n_elements(xyznew)

if row ne 1 then $
    for j=1,row-1 do begin
        temp=total(abs(xyznew - rebin(xyz[*,j],n,n_elements(xyznew)/n,/sample)) le SmallDev,1)
        ind=where(temp eq 3,ct)
        if ct eq 0 then begin
            xyznew=[[xyznew],[xyz[*,j]]]
            rownew++
        endif
    endfor

xyz=xyznew
row=rownew

end;function KeepUniXYZ
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CompressXYZ,xyz,row,SmallDev

; ----Convert to coordinates between 0 and 1----
ZOCompress0,xyz

; ----Keep only the unique positions----
KeepUniXYZ,xyz,row,SmallDev

end;pro CompressXYZ
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CoordPolyhedra,xyzi,labels,labeli,orbit=orbit,polyh=polyh,triangle=triangle

; xyzi: unit cell positions
; labels: orbit to which each position belongs to
; labeli: get the coordination of an atom in orbit "labeli"

; Duplicate unit cell in 3 directions
xyz=xyzi
xyz=[[xyz],[xyz+rebin([1,0,0],3,n)],[xyz+rebin([-1,0,0],3,n)]]
xyz=[[xyz],[xyz+rebin([0,1,0],3,n)],[xyz+rebin([0,-1,0],3,n)]]
xyz=[[xyz],[xyz+rebin([0,0,1],3,n)],[xyz+rebin([0,0,-1],3,n)]]

; Get the Delaunay triangles and the connectivity list
QHULL, xyz, triangle, /DELAUNAY, CONNECTIVITY = LIST

; Get the indices of the adjacent atoms
ind=where(labels eq labeli)
tmp=max(LIST[ind+1]-LIST[ind],i)
i=ind[i[0]]
indAdj = LIST[LIST[i] : LIST[i+1]-1]

; Adjacent atoms must belong to the same orbit
orbit=labels[indAdj[0]]
ind=where(labels[indAdj] ne orbit,ct)
if ct ne 0 then return,0

; Create the coordination polyhedron from the adjacent atoms
polyh = xyz[*,indAdj]
QHULL, polyh, triangle, /DELAUNAY

return,1
end;function CoordPolyhedra
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AppendXYZ,xyz,mult,xyz2,mult2,SmallDev,expanded

mult2=0
if n_elements(expanded) eq 0 then return

case expanded.plottype of
2:    begin
    xyz2=xyz
    
    ind=where(xyz2[0,*] ge -SmallDev and xyz2[0,*] le SmallDev,ct)
    if ct ne 0 then xyz2=[[xyz2],[xyz2[*,ind]+rebin([1,0,0],3,ct)]]
    mult2+=ct
    
    ind=where(xyz2[1,*] ge -SmallDev and xyz2[1,*] le SmallDev,ct)
    if ct ne 0 then xyz2=[[xyz2],[xyz2[*,ind]+rebin([0,1,0],3,ct)]]
    mult2+=ct
    
    ind=where(xyz2[2,*] ge -SmallDev and xyz2[2,*] le SmallDev,ct)
    if ct ne 0 then xyz2=[[xyz2],[xyz2[*,ind]+rebin([0,0,1],3,ct)]]
    mult2+=ct
    
    if mult2 eq 0 then xyz2=0L else xyz2=xyz2[*,mult:*]
    
    endcase
else: return
endcase

end;pro AppendXYZ
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddUnitCellpos,UC,SmallDev,xyzin,expanded=expanded

; xyz: all symmetry operations applied to one ASU position

bfirst=1b
bundefined=n_elements(xyzin) eq 0

if ~bundefined then begin
    xyz=xyzin
    nxyz=n_elements(xyz)/3.

    AppendXYZ,xyz,nxyz,xyz2,nxyz2,SmallDev,expanded

    t=replicate(-1,1,nxyz)
    t[0]=nxyz
    if nxyz2 ne 0 then begin
        t2=replicate(-2,1,nxyz2)
    endif
endif

; Are there already positions in UC?
if n_elements(UC) ne 0 then begin
    if bundefined then return
    if UC.nucpos ne 0 then begin

        UCxyz=[[UC.ucpos],[t,xyz]]
        if nxyz2 ne 0 then UCxyz=[[UCxyz],[t2,xyz2]]
        nucpos=UC.nucpos+1
        nucatoms=UC.nucatoms+nxyz+nxyz2
        reverse_indices2=UC.reverse_indices2

        bfirst=0b
    endif
endif

; First position or create empty structure
if bfirst then begin
    if bundefined then begin
        UCxyz=0L
        nucpos=0L
        nucatoms=0L
        reverse_indices2=0L
    endif else begin
        UCxyz=[t,xyz]
        if nxyz2 ne 0 then UCxyz=[[UCxyz],[t2,xyz2]]
        nucpos=1L
        nucatoms=nxyz+nxyz2
        reverse_indices2=0L
    endelse
endif

; Calculate reverence indices
if nucpos ne 0 then begin
    tmp=where(UCxyz[0,*] ge -1,nucatoms2)

    reverse_indices=where(UCxyz[0,*] gt -1,ct)
    reverse_indices=[reverse_indices,nucatoms]
    reverse_indices2=[reverse_indices2,nucatoms-nxyz2]
endif else begin
    nucatoms2=0L
    reverse_indices=0L
endelse

; Make new Unit Cell
UC={nucpos:nucpos,$         ; number of positions in the unit cell
    nucatoms:nucatoms,$        ; number of atoms in the expanded unit cell
    nucatoms2:nucatoms2,$    ; number of atoms in the unit cell
    ucpos:UCxyz,$            ; UC positions (columns):
    reverse_indices2:reverse_indices2,$
    reverse_indices:reverse_indices}
                            ; atoms on position i (i goes from 0 to nucpos-1)
                            ;     print,uc.ucpos[1:*,uc.reverse_indices[i]:uc.reverse_indices2[i+1]-1]
                            ; atoms on position i expanded (i goes from 0 to nucpos-1)
                            ;     print,uc.ucpos[1:*,uc.reverse_indices[i]:uc.reverse_indices[i+1]-1]
                            ; atoms on position i only expanded (i goes from 0 to nucpos-1)
                            ;     print,uc.ucpos[1:*,uc.reverse_indices2[i+1]:uc.reverse_indices[i+1]-1]
                            ;   => first check whether uc.reverse_indices2[i+1] ne uc.reverse_indices[i+1]

end;pro AddUnitCellpos
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ASUtoUnitCell,ptr,O=O,expanded=expanded

if keyword_set(O) then p=(*(*ptr).asu.O)[0] else p=(*ptr).asu.I
psg=(*ptr).pstrinfo
pri=(*ptr).preverse_indices

; ----Get precision----
tmp=getptrFit(ptr)
SmallDev=*((*tmp).misc[3])

; ----Get wyckoff operators----
str=sgdata((*psg).sghash,0)
wyck=call_function('wyckoff'+stringr(str.set),mult=mult)
tags=strlowcase(tag_names(wyck))

for j=0l,(*ptr).nasupos-1 do begin
    indj=(*pri)[j]
    ind=(where(tags eq (*p)[indj].wyck.L))[0]

    xyz=fltarr(3,mult[ind])
    for jj=0,mult[ind]-1 do begin
        RTjj=Symop64((wyck.(ind))[jj])
        xyz[*,jj]=RTjj.Rot##(*p)[indj].wyck.xyz+RTjj.Trn
    endfor
    ZOCompress0,xyz
    
    AddUnitCellpos,UC,SmallDev,xyz,expanded=expanded
endfor
AddUnitCellpos,UC,SmallDev,expanded=expanded

return,UC

end;function ASUtoUnitCell
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CopyASUstruc,ASU,ptr,recalc=recalc
if (n_elements(ASU) eq 0) or (*ptr).type ne 12 then return

if keyword_set(recalc) then begin ; only ASU.asupos
    (*ptr).nasuatoms=n_elements(ASU.asupos)
    (*ptr).nasupos=0
    ptr_free,(*ptr).preverse_indices
    if (*ptr).nasuatoms ne 0 then begin
        reverse_indices=where(ASU.asupos.natoms ne -1,nasupos)
        (*ptr).nasupos=nasupos
        reverse_indices=[reverse_indices,(*ptr).nasuatoms]
        (*ptr).preverse_indices=ptr_new(reverse_indices)
    endif
endif else begin
    (*ptr).nasuatoms=ASU.nasuatoms
    (*ptr).nasupos=ASU.nasupos
    ptr_free,(*ptr).preverse_indices
    if ASU.nasuatoms ne 0 then (*ptr).preverse_indices=ptr_new(ASU.reverse_indices)
endelse

heap_free,(*ptr).asu
if (*ptr).nasuatoms ne 0 then begin
    (*ptr).asu.I=ptr_new(ASU.asupos)
    (*ptr).asu.O=ptr_new([ptr_new(ASU.asupos)])

    tmp=ASU.asupos
    tmp.xyz=0
    tmp.SOF=0
    tmp.B=0
    (*ptr).asu.SD=ptr_new([ptr_new(temporary(tmp))])
endif

end;pro CopyASUstruc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CopystrucASU,ptr,O=O
if (*ptr).type ne 12 then return,{nasupos:0L,nasuatoms:0L,asupos:0L,reverse_indices:0L}

nasupos=(*ptr).nasupos
nasuatoms=(*ptr).nasuatoms
if ptr_valid((*ptr).preverse_indices) then reverse_indices=*(*ptr).preverse_indices $
else reverse_indices=0L
if keyword_set(O) then p=(*(*ptr).asu.O)[0] else p=(*ptr).asu.I
if ptr_valid(p) then asupos=*p else asupos=0L

ASU={nasupos:nasupos,$
    nasuatoms:nasuatoms,$
    asupos:asupos,$
    reverse_indices:reverse_indices}

return,ASU

end;function CopystrucASU
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CopystrucASU2,array,pasuinfo

; Array contains: xfit, yfit, zfit, SOF, B
off=0
asupos=(*pasuinfo)
nasuatoms=n_elements(asupos)
reverse_indices=where(asupos.natoms ne -1,nasupos)
reverse_indices=[reverse_indices,nasuatoms]

; xyzfit
ind=reverse_indices[0:nasupos-1]
induse=where(transpose(asupos[ind].wyck.xyz) ne 999,ctuse)
if ctuse ne 0 then begin
    xyz=replicate(999.,nasupos,3)
    xyz[induse]=array[off:off+ctuse-1]
    off+=ctuse
    xyz=transpose(xyz)

    ; xyz
    for i=0l,nasupos-1 do begin
        RT=Symop64(asupos[reverse_indices[i]].wyck.Rep)
        xyzi=reform(RT.rot##transpose(xyz[*,i]))+RT.trn
        ZOCompress0,xyzi

        nsub=reverse_indices[i+1]-reverse_indices[i]
        asupos[reverse_indices[i]:reverse_indices[i+1]-1].xyz=rebin(xyzi,3,nsub,/sample)
        asupos[reverse_indices[i]:reverse_indices[i+1]-1].wyck.xyz=rebin(xyz[*,i],3,nsub,/sample)
    endfor
endif

asupos.SOF=array[off:off+nasuatoms-1]
off+=nasuatoms
asupos.B=array[off:off+nasuatoms-1]
off+=nasuatoms

ASU={nasupos:nasupos,$
    nasuatoms:nasuatoms,$
    asupos:asupos,$
    reverse_indices:reverse_indices}

return,ASU
end;function CopystrucASU2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro printfasu,lun,ptr
if ptr_valid(ptr) then begin
    printf,lun,(*ptr).nasupos
    if (*ptr).nasupos ne 0L then begin
        printf,lun,(*ptr).nasuatoms
        p=(*ptr).asu.I
        printf,lun,(*p).natoms
        printf,lun,(*p).name,format='(A)'
        
        for i=0l,(*ptr).nasuatoms-1 do begin
            printf,lun,n_elements(*(*p)[i].Z)
            printf,lun,*(*p)[i].Z
            printf,lun,*(*p)[i].Ion
            printf,lun,*(*p)[i].Repeats
        endfor
        
        printf,lun,(*p).xyz
        printf,lun,(*p).SOF
        printf,lun,(*p).B
        printf,lun,(*p).Wyck.L,format='(A)'
        printf,lun,(*p).Wyck.M
        printf,lun,(*p).Wyck.Rep.R
        printf,lun,(*p).Wyck.Rep.T
        printf,lun,(*p).Wyck.xyz
        printf,lun,(*p).Wyck.togglefix
        printf,lun,(*p).togglefix
        printf,lun,*(*ptr).preverse_indices
    endif
endif else printf,lun,0L
end;pro printfasu
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro readfasu,lun,ptr
nasupos=0L
readf,lun,nasupos
if nasupos eq 0 then begin
    CopyASUstruc,{nasupos:0L,nasuatoms:0L},ptr
    return
endif

nasuatoms=0L
readf,lun,nasuatoms

struc=MakeASUpos(0L,'',0,0,0,fltarr(3),0.,0.)
for i=1,nasuatoms-1 do struc=[struc,MakeASUpos(0L,'',0,0,0,fltarr(3),0.,0.)]

tmp=lonarr(nasuatoms)
readf,lun,tmp
struc.natoms=tmp
tmp=strarr(nasuatoms)
readf,lun,tmp,format='(A)'
struc.name=tmp

for i=0l,nasuatoms-1 do begin
    n=0L
    readf,lun,n
    if n ne 0 then begin
        tmp=intarr(n)
        readf,lun,tmp
        *struc[i].Z=tmp
        readf,lun,tmp
        *struc[i].Ion=tmp
        readf,lun,tmp
        *struc[i].Repeats=tmp
    endif
endfor

tmp=fltarr(3,nasuatoms)
readf,lun,tmp
struc.xyz=tmp
tmp=fltarr(nasuatoms)
readf,lun,tmp
struc.SOF=tmp
readf,lun,tmp
struc.B=tmp

tmp=strarr(nasuatoms)
readf,lun,tmp,format='(A)'
struc.Wyck.L=tmp

tmp=intarr(nasuatoms)
readf,lun,tmp
struc.Wyck.M=tmp

tmp=lon64arr(nasuatoms)
readf,lun,tmp
struc.Wyck.Rep.R=tmp
readf,lun,tmp
struc.Wyck.Rep.T=tmp

tmp=fltarr(3,nasuatoms)
readf,lun,tmp
struc.Wyck.xyz=tmp

tmp=bytarr(3,nasuatoms)
readf,lun,tmp
struc.Wyck.togglefix=tmp

tmp=bytarr(2,nasuatoms)
readf,lun,tmp
struc.togglefix=tmp

reverse_indices=lonarr(nasupos+1)
readf,lun,reverse_indices

ASU={nasupos:nasupos,$
    nasuatoms:nasuatoms,$
    asupos:struc,$
    reverse_indices:reverse_indices}
CopyASUstruc,ASU,ptr

end;pro readfasu
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PlotASUAddLine,coord,oGroup,list
coord = list.NonOrtM#coord
oLine = OBJ_NEW('IDLgrPolyline' ,coord,XCOORD_CONV=NormCoord((list.list2).xr,list.PosX),$
                   YCOORD_CONV=NormCoord((list.list2).yr,list.PosY),$
                   ZCOORD_CONV=NormCoord((list.list2).zr,list.PosZ),thick=1)
oGroup->Add, oLine
end;pro PlotASUAddLine
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PlotASU, in

key=in.key
if key eq 0 then begin ; Change Z-range
    list=in.list

    ; Reset axis
    zr=(list.list2).zr
    list.oZAxis->SetProperty,range=zr,ZCoord_Conv=NormCoord(zr,list.PosZ)

endif else begin ; Add Objects
    list=in.list
    oGroup=list.oGroup

    ; ----Get data----
    plottype=*(list.list2).plottype
    PlotLabels=*(list.list2).PlotLabels
    ptr=(list.list2).ptr
    pri=(*ptr).preverse_indices
    if (list.list2).O then pasupos=(*(*ptr).asu.O)[0] else pasupos=(*ptr).asu.I
    UC=ASUtoUnitCell(ptr,O=(list.list2).O,expanded={plottype:plottype})
    celparam=(*(*ptr).globalparamIO.I)[indCelParam(ptr)]
    a=celparam[0]
    b=celparam[1]
    c=celparam[2]
    rmult=(*(list.list2).radmult)/max(celparam[0:2])

    ; ----Update list----
    list.NonOrtM=NonOrt2Ort(celparam[3],celparam[4],celparam[5])
    list.list2.xr=[0,a]
    list.oXAxis->SetProperty,range=(list.list2).xr,XCoord_Conv=NormCoord((list.list2).xr,list.PosX)
    list.list2.yr=[0,b]
    list.oYAxis->SetProperty,range=(list.list2).yr,YCoord_Conv=NormCoord((list.list2).yr,list.PosY)
    list.list2.zr=[0,c]
    list.oZAxis->SetProperty,range=(list.list2).zr,ZCoord_Conv=NormCoord((list.list2).zr,list.PosZ)

    ; ----Add unit cell----
    for ii=(list.list2).celx[0],(list.list2).celx[1] do $
        for jj=(list.list2).cely[0],(list.list2).cely[1] do $
            for kk=(list.list2).celz[0],(list.list2).celz[1] do begin
            
            ; ----Add box----
            
            ; Add faces separate
            x=reform(rebin([0,a],2,4,/sample),1,8)+ii*a
            y=reform(rebin([0,b],4,2,/sample),1,8)+jj*b
            z=reform(rebin([0,c],8,/sample),1,8)+kk*c
            verts = [x,y,z]
            faces=[[3,2,0,1,3],$ ; bottom (in default orientation)
                    [2,6,4,0,2],$ ; left
                    [2,3,7,6,2],$ ; back
                    [7,3,1,5,7],$ ; right
                    [6,7,5,4,6],$ ; top
                    [4,5,1,0,4]] ; front
            for facei=0,5 do begin
                coord = verts[*,faces[*,facei]]
                PlotASUAddLine,coord,oGroup,list
            endfor
            
;            ; Add centerings
;            x=reform([a/2, a/2,a/2, a/2,a/2, 0  ,  a],1,7)+ii*a
;            y=reform([b/2, b/2,b/2,   0,  b, b/2,b/2],1,7)+jj*b
;            z=reform([c/2,   0,  c, c/2,c/2, c/2,c/2],1,7)+kk*c
;            vertscen = [x,y,z]
;            
;            ; body-centered
;            for segi=0,7 do begin
;                coord = [[verts[*,segi]],[vertscen[*,0]]]
;                PlotASUAddLine,coord,oGroup,list
;            endfor
;            
;            ; base-centered C
;            for segi=0,7 do begin
;                coord = [[verts[*,segi]],[vertscen[*,1+(segi gt 3)]]]
;                PlotASUAddLine,coord,oGroup,list
;            endfor
;
;            ; base-centered B
;            for segi=0,7 do begin
;                coord = [[verts[*,segi]],[vertscen[*,3+(segi eq 2 or segi eq 3 or segi eq 6 or segi eq 7)]]]
;                PlotASUAddLine,coord,oGroup,list
;            endfor
;            
;            ; base-centered A
;            for segi=0,7 do begin
;                coord = [[verts[*,segi]],[vertscen[*,5+(segi mod 2)]]]
;                PlotASUAddLine,coord,oGroup,list
;            endfor
            

;            ;Add filled faces
;            if ii eq 0 and jj eq 0 and kk eq 0 then begin
;                
;                verts2=list.NonOrtM#verts
;                polygons = [4,6,7,5,4,$ ; top
;                            4,3,2,0,1,$ ; bottom
;                            4,7,3,1,5,$ ; right
;                            4,2,6,4,0,$ ; left
;                            4,2,3,7,6,$ ; back
;                            4,4,5,1,0]    ; front
;                oCel = OBJ_NEW('IDLgrPolygon' ,verts2,polygons=polygons,XCOORD_CONV=NormCoord((list.list2).xr,list.PosX),$
;                           YCOORD_CONV=NormCoord((list.list2).yr,list.PosY),$
;                           ZCOORD_CONV=NormCoord((list.list2).zr,list.PosZ),$
;                           thick=1,COLOR=[0b,255b,0b],ALPHA_CHANNEL=0.7,SHININESS=128,REJECT=0)
;                oGroup->Add, oCel
;            endif
            
            ; ----Add atoms and coordination polyhedra----
            ntype=(*ptr).nasupos
            for i=0l,ntype-1 do begin; Loop over all ASU positions
        
                result=AtomColRad(*(*pasupos)[(*pri)[i]].Z,*(*pasupos)[(*pri)[i]].Ion,*(*pasupos)[(*pri)[i]].Repeats)
                radiusi=result.R*rmult
                colori=result.Col
                labeli=(*pasupos)[(*pri)[i]].name
        
                ; ----Make atom of type i----
                oAtom = OBJ_NEW('AtomPlotObject', COLOR=colori, $
                    DENSITY=0.9, RADIUS=radiusi, STYLE=2, SHADING=1)
        
                ; ----Get all equivalent positions----
                case plottype of
                0: xyz=UC.ucpos[1:*,UC.reverse_indices[i]] ; ASU
                1: xyz=UC.ucpos[1:*,UC.reverse_indices[i]:UC.reverse_indices2[i+1]-1] ; unit cell
                2: xyz=UC.ucpos[1:*,UC.reverse_indices[i]:UC.reverse_indices[i+1]-1] ; expanded unit cell
                endcase
        
                n=n_elements(xyz)/3
                
                xyz+=rebin([ii,jj,kk],3,n) ; shift atoms for unit cell [ii,jj,kk]
        
                for j=0l,n-1 do begin; Loop over all positions generated for this ASU position
                    ; ----Make a new instance of the atom----
                    oModelAtom=oAtom->Copy()
                    ; ----Add a label to oModelAtom for the first atom in this ASU position----
                    if (PlotLabels eq 1) and (j eq 0) then begin
                        oModelAtom->SetProperty, UVALUE=labeli
                        oLabel=OBJ_NEW('IDLgrText',labeli,FONT=list.oLabelFont1,$
                                LOCATIONS=[list.PosX[1]*radiusi,0,0])
                        oModelAtom->Add, oLabel
                    endif
                    ; ----Add oModelAtom to oGroup and put on correct position----
                    PosAt=xyz[*,j]*[list.PosX[1],list.PosY[1],list.PosZ[1]]
                    PosAt=list.NonOrtM#PosAt
                    oModelAtom->Translate, PosAt[0], PosAt[1], PosAt[2]
                    oGroup->Add, oModelAtom
                endfor
        
                ; oAtom was not added to something, it was just used to generate some data
                obj_destroy,oAtom
            endfor    
    endfor

endelse

end;pro PlotASU
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PlotUC,ptr,base,update=update,O=O

if (*ptr).type ne 12 then return

celparam=(*(*ptr).globalparamIO.I)[indCelParam(ptr)]

if keyword_set(update) then begin
    TopBase=widget_info(base,find_by_uname='fplotTopBase')
    ID=widget_info(base,find_by_uname='fplotrebuildwidget')
    if widget_info(TopBase,/valid_id) and widget_info(ID,/valid_id) then begin
        WIDGET_CONTROL, TopBase, GET_UVALUE=list, /NO_COPY
        list.list2.xr[1]=celparam[0]
        list.list2.yr[1]=celparam[1]
        list.list2.zr[1]=celparam[2]
        WIDGET_CONTROL, TopBase, SET_UVALUE=list, /NO_COPY

        TLB=widgetTLB(base)
        widget_control,ID,send_event={ID:ID,TOP:TLB,HANDLER:TopBase}
    endif
endif else begin
    p=getptrFit(ptr)

    nunits=(*(*p).misc[7])[1,*]-(*(*p).misc[7])[0,*]
    
    nunitoff=abs(max((*(*p).misc[7]),/abs))
    zoom=0.5*nunitoff
    zclip=nunitoff+3
    

    list={    xr:[0,celparam[0]],$
        yr:[0,celparam[1]],$
        zr:[0,celparam[2]],$
        proplot:'PlotASU',$
        O:keyword_set(O),$
        plottype:(*p).misc[1],$    ; 0: assymmetric unit, 1: unit cel, 2: expanded unit cel
        PlotLabels:(*p).misc[2],$ ; 0 or 1
        celx:(*(*p).misc[7])[0:1,0],$
        cely:(*(*p).misc[7])[0:1,1],$
        celz:(*(*p).misc[7])[0:1,2],$
        radmult:(*p).misc[0],$
        ptr:ptr}

    flplot_obj,list,xtitle='X('+msymbols('angstroms')+')',ytitle='Y('+msymbols('angstroms')+')',ztitle='Z('+msymbols('angstroms')+')',/xtrue,/ytrue,/ztrue,$
            ayz=celparam[3],bxz=celparam[4],cxy=celparam[5],parentbase=base,xdim=400,ydim=400,zoom=zoom,zclip=zclip
endelse

end;pro PlotUC
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RefinedToInitialASU,ptr,inverse=inverse
update=0b
if (*ptr).type eq 12 then $
if (*ptr).nasuatoms ne 0 then begin
    icurrent=0
    if keyword_set(inverse) then begin
        pI=(*ptr).asu.I
        pO=(*(*ptr).asu.O)[icurrent]
        pSD=(*(*ptr).asu.SD)[icurrent]
    endif else begin
        pO=(*ptr).asu.I
        pI=(*(*ptr).asu.O)[icurrent]
    endelse

    (*pO)=(*pI)
    if n_elements(pSD) ne 0 then begin
        (*pSD)=(*pI)
        (*pSD).xyz=0
        (*pSD).SOF=0
        (*pSD).B=0
    endif
    update=1b
endif
return,update
end;function RefinedToInitialASU
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ConvertASUtotable,asupos

nx=9
ny=n_elements(asupos)
value=strarr(nx,ny)

value[0,*]=stringr(asupos.natoms)
value[1,*]=asupos.name
for i=0l,ny-1 do value[2,i]=ZElement(*asupos[i].Z,Ion=*asupos[i].Ion,Repeats=*asupos[i].Repeats)
value[3,*]=WyckLabel(asupos.Wyck)
value[4:6,*]=stringr(asupos.xyz)
value[7,*]=stringr(asupos.SOF)
value[8,*]=stringr(asupos.B)

return,value
end;function ConvertASUtotable
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditASUposXYZ,k,ipos,ifield,value,p,preverse_indices,psg,SmallDev,resval=resval

ind=where((*preverse_indices)-ipos le 0,ct)
ind=ind[ct-1]

; The same XYZ for substituted atoms
b=(*preverse_indices)[ind]
e=(*preverse_indices)[ind+1]-1
n=e-b+1
ind=b+indgen(n)

; Keep in [0,1[
ZOCompress0,value
(*p)[ind].xyz[k]=value[0]

; Recalculate wyckoff
asupos=(*p)[b]
SetASUwyckoff,*psg,asupos,SmallDev
for i=b,e do (*p)[i].wyck=asupos.wyck

resval=[[replicate(value,n)],[replicate(ifield,n)],[ind]]
end;pro EditASUposXYZ
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditASUpos,p,preverse_indices,psg,ifield,ipos,value,SmallDev,resval=resval,obsolete=obsolete

; ifield: table column index
; ipis: table row index

obsolete=1b

case ifield of
0:     begin ; not editable
    resval=(*p)[ipos].natoms
    resval=[[resval],[ifield],[ipos]]
    obsolete and= 0b
    endcase
1:     (*p)[ipos].name=value
2:     begin
    Z=ZElement(value,ion=ion,repeats=repeats)
    *(*p)[ipos].Z=Z
    *(*p)[ipos].Ion=ion
    *(*p)[ipos].Repeats=repeats
    endcase
3:     begin ; not editable
    resval=WyckLabel((*p)[ipos].Wyck)
    resval=[[resval],[ifield],[ipos]]
    obsolete=1b
       endcase
4:     EditASUposXYZ,0,ipos,ifield,rationalStringToFloat(value),p,preverse_indices,psg,SmallDev,resval=resval
5:     EditASUposXYZ,1,ipos,ifield,rationalStringToFloat(value),p,preverse_indices,psg,SmallDev,resval=resval
6:     EditASUposXYZ,2,ipos,ifield,rationalStringToFloat(value),p,preverse_indices,psg,SmallDev,resval=resval
7:     (*p)[ipos].SOF=float(value)
8:     (*p)[ipos].B=float(value)
endcase

end;pro EditASUpos
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddASUposI,ptr,ind

if ~ptr_valid((*ptr).pstrinfo) then return

p=(*ptr).asu.I
if n_elements(ind) ne 0 then begin
    ; Reset .natoms field
    i=ind[0]
    j=where((*(*ptr).preverse_indices)-i le 0,ct)
    j=j[ct-1]

    (*p)[(*(*ptr).preverse_indices)[j]].natoms++
    add=(*p)[(*(*ptr).preverse_indices)[j]]
    add.natoms=-1

    ; Insert atom
    asupos=[(*p)[0:i],add]
    if i ne (*ptr).nasuatoms-1 then asupos=[asupos,(*p)[i+1:*]]
    heap_copy,asupos
endif else begin
    asupos=MakeASUpos(1L,'',1,1,1,fltarr(3),1.,1.)

    ; ----Get precision----
    tmp=getptrFit(ptr)
    SmallDev=*((*tmp).misc[3])

    SetASUwyckoff,*(*ptr).pstrinfo,asupos,SmallDev
    if (*ptr).nasupos ne 0 then begin
        tmp=(*p)
        heap_copy,tmp
        asupos=[tmp,asupos]
    endif
endelse

CopyASUstruc,{asupos:asupos},ptr,/recalc
end;pro AddASUposI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeleteASUposI,ptr,ind

; Delete All?
n=n_elements(ind)
if n ge (*ptr).nasuatoms then begin
    CopyASUstruc,{nasupos:0L,nasuatoms:0L},ptr
    return
endif

; Reset .natoms fields
icurrent=0
p=(*ptr).asu.I
for i=0l,n-1 do begin
    j=where((*(*ptr).preverse_indices)-ind[i] le 0,ct)
    j=j[ct-1]

    a=(*(*ptr).preverse_indices)[j]
    b=(*(*ptr).preverse_indices)[j+1]-1
    (*p)[a].natoms--

    ; Handle deleting first atom of group
    if (a eq ind[i]) and (a ne b) then (*p)[a+1].natoms=(*p)[a].natoms
endfor

; Delete atoms
b=replicate(1b,(*ptr).nasuatoms)
b[ind]=0b

keep=(*p)[where(b)]
heap_copy,keep
CopyASUstruc,{asupos:keep},ptr,/recalc

end;pro DeleteASUposI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StructureFactorpder,ASU,HKL,psg,SmallDev,u,anomal=anomal,Lamb=Lamb,$
    pder=pder,bpder=bpder

bpder=keyword_set(bpder)
; HKL={hkl:hkl,d:d,mult:mult,n:n} ; cfr. above
; ASU=[nasuatoms,...] ; cfr. above
; Return: I_hkl = mult.|F_hkl|^2.
; pder: dI_hkl/u_hkl,xfit,yfit,zfit,SOF,B

; Get Unit cell positions xyz and info for derivative to xfit,yfit,zfit
str=sgdata((*psg).sghash,0)
wyck=call_function('wyckoff'+stringr(str.set),mult=mult)
tags=strlowcase(tag_names(wyck))

ind=ASU.reverse_indices[0:ASU.nasupos-1]
nxyz=total(ASU.asupos[ind].wyck.M)
xyz=fltarr(3,nxyz)
if bpder then Tjmpder=fltarr(3,3,nxyz)
tmp=total(ASU.asupos[ind].wyck.M,/cumul,/pres)
if ASU.nasupos gt 1 then m0=[0,tmp[0:ASU.nasupos-2]] else m0=[0]
m1=tmp-1

for j=0l,ASU.nasupos-1 do begin
    indj=ASU.reverse_indices[j]
    ind=(where(tags eq ASU.asupos[indj].wyck.L))[0]

    for jj=0,mult[ind]-1 do begin
        RTjj=Symop64((wyck.(ind))[jj])
        off=m0[j]+jj
        xyz[*,off]=RTjj.Rot##ASU.asupos[indj].wyck.xyz+RTjj.Trn
        if bpder then Tjmpder[*,*,off]=RTjj.Rot
    endfor
endfor
ZOCompress0,xyz

; Structure Factor
F_HKL=DCOMPLEXARR(HKL.n)

if bpder then begin
    ct=!pi/360.

    ; Partial derivatives
    ; dI_hkl/dp:
    ; dI_H0/du0, dI_H1/du1, ... dI_Hn/dun (u=2theta(deg))
    ; dI_H0/dx0, dI_H1/dx0, ... dI_Hn/dx0
    ; ... other x's
    ; dI_H0/dy0, dI_H1/dy0, ... dI_Hn/dy0
    ; ... other y's
    ; dI_H0/dz0, dI_H1/dz0, ... dI_Hn/dz0
    ; ... other z's
    ; dI_H0/dSOF0, dI_H1/dSOF0, ... dI_Hn/dSOF0
    ; ... other SOF's
    ; dI_H0/dB0, dI_H1/dB0, ... dI_Hn/dB0
    ; ... other B's

    ; Row indices in pder
    tmp=ASU.reverse_indices[0:ASU.nasupos-1]
    xyzref=long(ASU.asupos[tmp].wyck.xyz ne 999)
    if ASU.nasupos gt 1 then tmp=total(xyzref,2,/pres) else tmp=xyzref
    npderrow=[1L,tmp,ASU.nasuatoms,ASU.nasuatoms]

    pder=dblarr(HKL.n,total(npderrow,/pres))

    ; Indices u
    outindu=[0]
    off=npderrow[0]
    ; Indices xfit
    indpderx=replicate(-1,ASU.nasupos)
    if npderrow[1] ne 0 then begin
        pderx_A=dblarr(npderrow[1])
        pderx_B=dblarr(npderrow[1])
        outindx=indgen(npderrow[1])
        indpderx[where(xyzref[0,*])]=outindx
        outindx+=off
    endif
    off+=npderrow[1]
    ; Indices yfit
    indpdery=replicate(-1,ASU.nasupos)
    if npderrow[2] ne 0 then begin
        pdery_A=dblarr(npderrow[2])
        pdery_B=dblarr(npderrow[2])
        outindy=indgen(npderrow[2])
        indpdery[where(xyzref[1,*])]=outindy
        outindy+=off
    endif
    off+=npderrow[2]
    ; Indices zfit
    indpderz=replicate(-1,ASU.nasupos)
    if npderrow[3] ne 0 then begin
        pderz_A=dblarr(npderrow[3])
        pderz_B=dblarr(npderrow[3])
        outindz=indgen(npderrow[3])
        indpderz[where(xyzref[2,*])]=outindz
        outindz+=off
    endif
    off+=npderrow[3]
    ; Indices SOF
    if npderrow[4] ne 0 then begin
        pderSOF_A=dblarr(npderrow[4])
        pderSOF_B=dblarr(npderrow[4])
        indpderSOF=ASU.reverse_indices[0:ASU.nasupos-1]
        outindSOF=indgen(npderrow[4])+off
    endif else indpderSOF=replicate(-1,ASU.nasupos)
    off+=npderrow[4]
    ; Indices BISO
    if npderrow[5] ne 0 then begin
        pderBISO_A=dblarr(npderrow[5])
        pderBISO_B=dblarr(npderrow[5])
        indpderBISO=ASU.reverse_indices[0:ASU.nasupos-1]
        outindBISO=indgen(npderrow[5])+off
    endif else indpderBISO=replicate(-1,ASU.nasupos)

endif

for i=0l,HKL.n-1 do begin ; loop over all peaks
    hi=HKL.hkl[0,i]
    ki=HKL.hkl[1,i]
    li=HKL.hkl[2,i]
    di=HKL.d[i]

    A_H=0.
    B_H=0.

    if bpder then begin
        pderui_A=0d
        pderui_B=0d
    endif

    for j=0l,ASU.nasupos-1 do begin ; loop over all positions j in the asu
        ; Number of generated positions from ASU position j
        nj=m1[j]-m0[j]+1

        ; Number of substitutions on position j
        k0=ASU.reverse_indices[j]
        k1=ASU.reverse_indices[j+1]-1
        nsubj=k1-k0+1

        ; Get all (=nj) xyz for position j
        Tjm=2*!dpi*(hi*xyz[0,m0[j]:m1[j]]+ki*xyz[1,m0[j]:m1[j]]+li*xyz[2,m0[j]:m1[j]])
        cosTjm=rebin(cos(Tjm),nsubj,nj,/sample)
        sinTjm=rebin(sin(Tjm),nsubj,nj,/sample)

        ; Get all (=nsubj) atomic formfactors for position j
        f=dcomplexarr(nsubj)
        for k=k0,k1 do $
            f[k-k0]=AtScat(*ASU.asupos[k].Z,*ASU.asupos[k].Ion,*ASU.asupos[k].Repeats,di,anomal=anomal,Lamb=Lamb)
        ; TODO: Save CromerMann/WaasKirff coeff for subsequent calls to StructureFactorpder

        f_A=rebin(double(f),nsubj,nj,/sample)
        f_B=rebin(imaginary(f),nsubj,nj,/sample)

        if bpder then begin
            ; Get all Debye-Waller factors * SOF for position j
            DW=rebin([exp(-ASU.asupos[k0:k1].B/4./di/di)],nsubj,nj,/sample)
            Biso=rebin([ASU.asupos[k0:k1].B],nsubj,nj,/sample)
            SDW=rebin([ASU.asupos[k0:k1].SOF],nsubj,nj,/sample)*DW

            ; Get elements of the structure factor due to position j
            argA=f_A*cosTjm-f_B*sinTjm
            argB=f_A*sinTjm+f_B*cosTjm
            Aj=SDW*argA
            Bj=SDW*argB

            ; Get all j and m dependent forfactors
            ; xfitj
            if indpderx[j] ne -1 then begin
                tmp= rebin(reform($
                    hi*(Tjmpder[0,0,m0[j]:m1[j]])+$
                    ki*(Tjmpder[0,1,m0[j]:m1[j]])+$
                    li*(Tjmpder[0,2,m0[j]:m1[j]]),1,nj),nsubj,nj,/sample)

                pderx_A[indpderx[j]]=-total(Bj*tmp)
                pderx_B[indpderx[j]]=total(Aj*tmp)
            endif
            ; yfitj
            if indpdery[j] ne -1 then begin
                tmp= rebin(reform($
                    hi*(Tjmpder[1,0,m0[j]:m1[j]])+$
                    ki*(Tjmpder[1,1,m0[j]:m1[j]])+$
                    li*(Tjmpder[1,2,m0[j]:m1[j]]),1,nj),nsubj,nj,/sample)

                pdery_A[indpdery[j]]=-total(Bj*tmp)
                pdery_B[indpdery[j]]=total(Aj*tmp)
            endif
            ; zfitj
            if indpderz[j] ne -1 then begin
                tmp= rebin(reform($
                    hi*(Tjmpder[2,0,m0[j]:m1[j]])+$
                    ki*(Tjmpder[2,1,m0[j]:m1[j]])+$
                    li*(Tjmpder[2,2,m0[j]:m1[j]]),1,nj),nsubj,nj,/sample)

                pderz_A[indpderz[j]]=-total(Bj*tmp)
                pderz_B[indpderz[j]]=total(Aj*tmp)
            endif
            ; uH
            pderui_A+=total(Aj*Biso)
            pderui_B+=total(Bj*Biso)
            ; SOF
            if indpderSOF[j] ne -1 then begin
                pderSOF_A[indpderSOF[j]:indpderSOF[j]+nsubj-1]=(nj eq 1)?(DW*argA):total(DW*argA,2)
                pderSOF_B[indpderSOF[j]:indpderSOF[j]+nsubj-1]=(nj eq 1)?(DW*argB):total(DW*argB,2)
            endif
            ; Biso
            if indpderBISO[j] ne -1 then begin
                pderBISO_A[indpderBISO[j]:indpderBISO[j]+nsubj-1]=(nj eq 1)?Aj:total(Aj,2)
                pderBISO_B[indpderBISO[j]:indpderBISO[j]+nsubj-1]=(nj eq 1)?Bj:total(Bj,2)
            endif

        endif else begin
            ; Get all Debye-Waller factors * SOF for position j
            SDW=rebin([ASU.asupos[k0:k1].SOF*exp(-ASU.asupos[k0:k1].B/4./di/di)],nsubj,nj,/sample)

            ; Get elements of the structure factor due to position j
            Aj=SDW*(f_A*cosTjm-f_B*sinTjm)
            Bj=SDW*(f_A*sinTjm+f_B*cosTjm)
        endelse

        ; Sum over three dimensions:
        A_H+=total(Aj,/pres)
        B_H+=total(Bj,/pres)

    endfor

    F_HKL[i]=DCOMPLEX(A_H,B_H)

    if bpder then begin
        A_H*=2*HKL.mult[i]
        B_H*=2*HKL.mult[i]

        pder[i,outindu]=(-ct/Lamb^2.*sin(2*ct*u[i]))*(A_H*pderui_A+B_H*pderui_B)
        if npderrow[1] ne 0 then pder[i,outindx]=(A_H*pderx_A+B_H*pderx_B)*2*!pi
        if npderrow[2] ne 0 then pder[i,outindy]=(A_H*pdery_A+B_H*pdery_B)*2*!pi
        if npderrow[3] ne 0 then pder[i,outindz]=(A_H*pderz_A+B_H*pderz_B)*2*!pi
        if npderrow[4] ne 0 then pder[i,outindSOF]=A_H*pderSOF_A+B_H*pderSOF_B
        if npderrow[5] ne 0 then pder[i,outindBISO]=-(sin(ct*u[i])/Lamb)^2.*(A_H*pderBISO_A+B_H*pderBISO_B)
    endif

endfor

return,[[abs(F_HKL)],[reform(hkl.mult)]]
end;function StructureFactorpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Afactorpder,tti,info,pder=pder,bpder=bpder
; Not used yet!

tt=tti/180*!pi
case info.LGPtype of
1:    F=tt*0+0.5 ; Bragg-Brentano with infinitely thick sample
2:    F=0.5*(1-exp(-2*info.muL*info.D/sin(tt/2.))) ; Bragg-Brentano
3:  F=cos(tt)/(cos(tt)-1)*(exp(-info.muL*info.D/cos(tt))-exp(-info.muL*info.D)) ; Transmission perpendicular beam
else: F=tt*0+1 ; attenuation correction is 1
endcase

if n_elements(bpder) eq 0 then bpder=0
if bpder then begin
    case info.LGPtype of
    2: pder=-0.5*info.muL*info.D*cos(tt/2)*exp(-2*info.muL*info.D/sin(tt/2))/(sin(tt/2))^2
    3: pder=-sin(tt)*(exp(-info.muL*info.D/cos(tt))-exp(-info.muL*info.D))/(cos(x)-1) $
            +sin(tt)*cos(tt)*(exp(-info.muL*info.D/cos(tt))-exp(-info.muL*info.D))/(cos(x)-1)^2 $
            -info.muL*info.D*sin(tt)*exp(-info.muL*info.D/cos(tt))/(cos(tt)*(cos(x)-1))
    else: pder=tt*0
    endcase
endif

return,F

end;function Afactorpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LGPfactorpder,tti,info,pder=pder,bpder=bpder

tt=tti*!pi/180
case info.LGPtype of
1:    begin
    K = (1-info.pol_Plindeg*cos(info.pol_azimuth*!pi/90))/2

    Pol = 1-K+K*(cos(tt))^2
    F = Pol/(sin(tt/2.)*sin(tt))
    endcase
2:    begin
    m = info.pol_mmult*info.pol_n
    mono = (abs(cos(2*asin(info.lambda/(2*info.pol_dmono)))))^m
    Plindeg = (info.pol_s*(1+info.pol_s*info.pol_Plindeg)-info.pol_s*mono*(1-info.pol_s*info.pol_Plindeg))/((1+info.pol_s*info.pol_Plindeg)+mono*(1-info.pol_s*info.pol_Plindeg))
    K = (1-Plindeg*cos(info.pol_azimuth*!pi/90))/2
    
    Pol = 1-K+K*(cos(tt))^2
    F = Pol/(sin(tt/2.)*sin(tt))
    endcase
3:    begin
    m = info.pol_mmult*info.pol_n
    mono = (abs(cos(2*asin(info.lambda/(2*info.pol_dmono)))))^m
    K = (1-info.pol_Plindeg*info.pol_s)/2

    Pol = 1-K+K*mono*cos(tt)^2.
    F = Pol/(sin(tt/2.)*sin(tt))
    endcase
4:    begin
    F=1/sin(tt)
    endcase
else: F=tt*0+1
endcase

if n_elements(bpder) eq 0 then bpder=0
if bpder then begin
    case info.LGPtype of
    1:    begin
        t=tt/2.
        pder=-!pi/180.*(  2*K*cos(tt)/sin(t)+Pol*cos(tt)/((sin(tt))^2*sin(t))+Pol*cos(t)/((sin(t))^2*sin(tt))  )
        endcase
    2:    begin
        t=tt/2.
        pder=-!pi/180.*(  2*K*cos(tt)/sin(t)+Pol*cos(tt)/((sin(tt))^2*sin(t))+Pol*cos(t)/((sin(t))^2*sin(tt))  )
        endcase
    3:    begin
        t=tt/2.
        pder=-!pi/180.*(  2*K*mono*cos(tt)/sin(t)+Pol*cos(tt)/((sin(tt))^2*sin(t))+Pol*cos(t)/((sin(t))^2*sin(tt))  )
        endcase
    4:    begin
        pder=-!pi/180.*cos(tt)/(sin(tt))^2
        endcase
    else: pder=tt*0
    endcase
endif

return,F
end;function LGPfactorpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetPDRelIntpder,I,twotheta,typeinfo,keep=keep,bpder=bpder,noglobal=noglobal

; I: Scaling, m.|Fhkl|^2
if keyword_set(noglobal) then begin
    S=1.
    param=I
endif else begin
    S=I[0]
    param=I[1:*]
endelse

; S.m_H.|F_H|^2.LGP
LGP=LGPfactorpder(twotheta,typeinfo,pder=keep,bpder=bpder)
ret=param*LGP*S

; I_H=S.m_H.|F_H|^2.LGP
if keyword_set(bpder) then begin
    keep=[    [param*LGP],$    ; dI_H/dS
            [keep*param*S],$; dI_H/duH
            [S*LGP]]         ; dI_H/d(m_H.|F_H|^2)
endif

return,ret
end;function GetPDRelIntpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetPawleypeakparam,ptr,lambda,icurrent,O=O

if keyword_set(O) then begin
    pgi=(*(*ptr).peakparamIO.O)[icurrent]
    pg=(*(*ptr).globalparamIO.O)[icurrent]
endif else begin
    pgi=(*ptr).peakparamIO.I
    pg=(*ptr).globalparamIO.I
endelse

ind=*(*ptr).paramind[1,0]
; I=m_H.|F_H|^2
d=dfromhkl((*pgi)[*(*ptr).paramind[0,0],*],(*pg)[indCelParam(ptr)])
tt=BraggDtoX(d,lambda,/onlyt,/angledeg)
ret=(*pgi)[ind[0],*]*(*pgi)[ind[1],*]

n=(*ptr).peaknr
return,[reform(tt,1,n),ret]
end;function GetPawleypeakparam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetPawleyRelIntpder,I,twotheta,typeinfo,Fcalc=Fcalc,keep=keep,bpder=bpder,noglobal=noglobal

; I: Scaling, |Fhkl|^2, mhkl
if keyword_set(noglobal) then begin
    n=n_elements(I)/2
    S=1.
    m=I[0:n-1]
    F=I[n:2*n-1]
endif else begin
    n=(n_elements(I)-1)/2
    S=I[0]
    m=I[1:n]
    F=I[n+1:2*n]
endelse

; m_H.|F_H|^2.LGP
LGP=LGPfactorpder(twotheta,typeinfo,pder=keep,bpder=bpder)
ret=S*m*F*LGP

; Fcalc
Fcalc=[[sqrt(abs(F))],$; |F_H|
       [S*m*LGP]]; S.m_H.LGP
       
; I_H=S.m_H.|F_H|^2.LGP
if keyword_set(bpder) then begin
    keep=[    [m*F*LGP],$        ; dI_H/dS
            [S*m*F*keep],$    ; dI_H/duH
            [S*m*LGP]]         ; dI_H/d|F_H|^2
endif

return,ret
end;function GetPawleyRelIntpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetRietveldpeakparam,ptr,lambda,icurrent,O=O,PD=PD
if keyword_set(O) then begin
    p=(*(*ptr).asu.O)[icurrent]
    pgi=(*(*ptr).peakparamIO.O)[icurrent]
    pg=(*(*ptr).globalparamIO.O)[icurrent]
endif else begin
    p=(*ptr).asu.I
    pgi=(*ptr).peakparamIO.I
    pg=(*ptr).globalparamIO.I
endelse
psg=(*ptr).pstrinfo

; HKL data
n=(*ptr).peaknr
hkl=(*pgi)[*(*ptr).paramind[0,0],*]
mult=(*pgi)[*(*ptr).paramind[1,0],*]
d=dfromhkl(hkl,(*pg)[indCelParam(ptr)])
HKL={hkl:hkl,d:d,mult:mult,n:n}
; ASU data
ASU=CopystrucASU(ptr,O=O)
if ASU.nasupos eq 0 then return,[0,1]
; Other data
tmp=getptrFit(ptr)
SmallDev=*((*tmp).misc[3])
anomal=*((*tmp).misc[4])

tt=BraggDtoX(d,lambda,/onlyt,/angledeg)
ret=StructureFactorpder(ASU,HKL,(*ptr).pstrinfo,SmallDev,anomal=anomal,Lamb=lambda)

ret=transpose(ret)
if keyword_set(PD) then return,[reform(tt,1,n),ret[0,*]*ret[0,*]*ret[1,*]] $ ; tt,|F_H|^2.m_H
else return,[ret[1,*],ret[0,*]*ret[0,*]] ; m_H, |F_H|^2

end;function GetRietveldpeakparam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetRietveldRelIntpder,I,twotheta,hkl,typeinfo,Fcalc=Fcalc,keep=keep,bpder=bpder

; I: Scaling, m, xfit, yfit, zfit, SOF, B
; HKLdata
d=BraggTtoX(twotheta,typeinfo.lambda,/onlyd,/angledeg)
n=n_elements(d)
HKL={hkl:hkl,d:d,mult:I[1:n],n:n}
; ASU data
ASU=CopystrucASU2(I[n+1:*],typeinfo.pasuinfo)

; S.m_H.|F_H|^2.LGP
LGP=LGPfactorpder(twotheta,typeinfo,pder=keep2,bpder=bpder)
mFH2=StructureFactorpder(ASU,HKL,typeinfo.pstrinfo,typeinfo.SmallDev,twotheta,anomal=typeinfo.anomal,$
    Lamb=typeinfo.lambda,pder=keep,bpder=bpder)

Fcalc=mFH2; [|F_H|,m_H]
Fcalc[*,1]*=I[0]*LGP ; m_H.S.LGP
mFH2=mFH2[*,0]*mFH2[*,0]*mFH2[*,1];m_H.|F_H|^2
ret=I[0]*mFH2*LGP

; I_H=S.m_H.|F_H|^2.LGP
if keyword_set(bpder) then begin
    
    keep2=[    [mFH2*LGP],$ ; dI_H/dS
            [I[0]*(LGP*keep[*,0]+mFH2*keep2)]] ;dI_H/duH

    s=DimSize(keep,2)
    s[1]--
    if s[1] gt 0 then begin
        keep2=[    [keep2],$
                [keep[*,1:s[1]]*rebin(I[0]*LGP,s[0],s[1],/sample)]] ; dI_H/d...(atoms)        
    endif
    keep=temporary(keep2)
endif

return,ret
end;function GetRietveldRelIntpder
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calcElasticScatteringFraction,ptr,energy,thickness,icurrent,O=O
; Analog to UnitCellFundamentals12
; lambda: keV
; thickness: mm

if keyword_set(O) then begin
    p=(*(*ptr).asu.O)[icurrent]
    psd=(*(*ptr).asu.SD)[icurrent]
    pg=(*(*ptr).globalparamIO.O)[icurrent]
    pgsd=(*(*ptr).globalparamIO.SD)[icurrent]
    bI=0b
endif else begin
    p=(*ptr).asu.I
    psd=(*(*ptr).asu.SD)[icurrent]
    pg=(*ptr).globalparamIO.I
    pgsd=(*(*ptr).globalparamIO.SD)[icurrent]
    bI=1b
endelse

; Calculate energy dependent parameters
lambda=LambdaEnergy(energy,/keV)
n=n_elements(lambda)
UCMMass=0.
UCMu=fltarr(n)
UCCoh=UCMu
for j=0l,n-1 do begin
    UCMuj=0.
    UCCohj=0.
    
    ; Loop over positions in asym. unit
    for i=0l,(*ptr).nasupos-1 do begin
        k0=(*(*ptr).preverse_indices)[i]
        kn=(*(*ptr).preverse_indices)[i+1]-1
    
        ; Loop over substituted atoms on this ASU position
        m=((*p)[k0]).wyck.m ; multiplicity of position
        for k=k0,kn do begin
            mu=MassAtCoef(*(*p)[k].Z,*(*p)[k].Ion,*(*p)[k].Repeats,lambda[j],mmatom=mmatom)
            coh=CohCrossSection(*(*p)[k].Z,*(*p)[k].Ion,*(*p)[k].Repeats,lambda[j]);,mmatom=mmatom
    
            tmp=(*p)[k].SOF*mmatom*m
            if j eq 0 then UCMMass+=tmp
            UCMuj+=tmp*mu
            UCCohj+=tmp*coh
        endfor
    endfor
    
    UCMu[j]=UCMuj
    UCCoh[j]=UCCohj
endfor
UCMu/=UCMMass
UCCoh/=UCMMass

celparam=(*pg)[indCelParam(ptr)]
celparam_sd=double((*pgsd)[indCelParam(ptr)])
if bI then celparam_sd*=0
UCV=UCVolume(celparam,pder=UCV_sd,connect=celparamconnect(*(*ptr).pstrinfo))

ct=avogadroconstant(/noexp)
UCrho=UCMMass / ct / UCV

; Calculate the fraction of I0 that is elastically scattered
; (assuming all scattering is in the forward direction) with
; attenuation taking into account

nx = n_elements(energy)
ny = n_elements(thickness)
UCrho_r = rebin([UCrho],nx,ny,/sample)
UCCoh_r = rebin([UCCoh],nx,ny,/sample)
UCMu_r = rebin([UCMu],nx,ny,/sample)
thickness_r = rebin(reform([thickness],1,ny)/10.,nx,ny,/sample) ; mm -> cm

return, thickness_r*UCrho_r*UCCoh_r*exp(-UCMu_r*thickness_r*UCrho_r) ; number between 0 and 1
end;function calcElasticScatteringFraction
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_ElasticScatteringFractionWidget,ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
if WIDGET_INFO(list.top,/VALID_ID) then WIDGET_CONTROL, list.top, sensitive=1
end;pro CleanUp_ElasticScatteringFractionWidget
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ElasticScatteringFractionWidget_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
0:begin ; base
  endcase
1:begin ; button
  widget_control,ev.id,get_uvalue=uval

  case uval of
  'energy': begin
            ID=widget_info(ev.top,FIND_BY_UNAME='O')
            O = widget_info(ID,/BUTTON_SET)
            
            ID=widget_info(ev.top,FIND_BY_UNAME='emin')
            widget_control,ID,get_value=emin
            emin = double(emin)
            ID=widget_info(ev.top,FIND_BY_UNAME='emax')
            widget_control,ID,get_value=emax
            emax = double(emax)
            ID=widget_info(ev.top,FIND_BY_UNAME='einc')
            widget_control,ID,get_value=einc
            einc = double(einc)
            ID=widget_info(ev.top,FIND_BY_UNAME='thickness')
            widget_control,ID,get_value=thickness
            thickness = double(thickness)
            
            n = round((emax-emin)/einc+1)
            energy = emin[0] + einc[0]*lindgen(n)
            
            widget_control,ev.top,get_uvalue=list
            fraction=calcElasticScatteringFraction(list.ptr,energy,thickness,0,O=O)
            
            p=plot(energy,fraction*100,title='Interaction with '+stringr(thickness)+' mm '+(*list.ptr).name,$
                xtitle='Energy (keV)',$
                ytitle='Forward elastic scattering (% of I0)')

            m=max(fraction,indm)
            yrange=[min(fraction),m]*100
            line=POLYLINE(energy[[indm,indm]],yrange,/data,/current,color=!color.red,thick=2)
            txt=TEXT(energy[indm],(yrange[1]-yrange[0])/2.,stringr(energy[indm]),/data,/current,color=!color.red)
            
            ID=widget_info(ev.top,FIND_BY_UNAME='energy')
            widget_control,ID,set_value=stringr(energy[indm])
            (*list.p).energy = energy[indm]
            endcase
  'thickness':begin
            ID=widget_info(ev.top,FIND_BY_UNAME='O')
            O = widget_info(ID,/BUTTON_SET)
            
            ID=widget_info(ev.top,FIND_BY_UNAME='tmin')
            widget_control,ID,get_value=tmin
            tmin = double(tmin)
            ID=widget_info(ev.top,FIND_BY_UNAME='tmax')
            widget_control,ID,get_value=tmax
            tmax = double(tmax)
            ID=widget_info(ev.top,FIND_BY_UNAME='tinc')
            widget_control,ID,get_value=tinc
            tinc = double(tinc)
            ID=widget_info(ev.top,FIND_BY_UNAME='energy')
            widget_control,ID,get_value=energy
            energy = double(energy)
          
            n = round((tmax-tmin)/tinc+1)
            thickness = tmin[0] + tinc[0]*lindgen(n)
            
            widget_control,ev.top,get_uvalue=list
            fraction=calcElasticScatteringFraction(list.ptr,energy,thickness,0,O=O)
            
            p=plot(thickness,fraction*100,title='Interaction with '+(*list.ptr).name+' at '+strtrim(energy,2)+' keV.',$
                xtitle='Thickness (mm)',$
                ytitle='Forward elastic scattering (% of I0)')

            m=max(fraction,indm)
            yrange=[min(fraction),m]*100
            line=POLYLINE(thickness[[indm,indm]],yrange,/data,/current,color=!color.red,thick=2)
            txt=TEXT(thickness[indm],(yrange[1]-yrange[0])/2.,stringr(thickness[indm]),/data,/current,color=!color.red)
            
            ID=widget_info(ev.top,FIND_BY_UNAME='thickness')
            widget_control,ID,set_value=stringr(thickness[indm])
            (*list.p).thickness = thickness[indm]
            endcase
  'both':   begin
            ID=widget_info(ev.top,FIND_BY_UNAME='emin')
            widget_control,ID,get_value=emin
            emin = double(emin)
            ID=widget_info(ev.top,FIND_BY_UNAME='emax')
            widget_control,ID,get_value=emax
            emax = double(emax)
            ID=widget_info(ev.top,FIND_BY_UNAME='einc')
            widget_control,ID,get_value=einc
            einc = double(einc)
            
            ID=widget_info(ev.top,FIND_BY_UNAME='tmin')
            widget_control,ID,get_value=tmin
            tmin = double(tmin)
            ID=widget_info(ev.top,FIND_BY_UNAME='tmax')
            widget_control,ID,get_value=tmax
            tmax = double(tmax)
            ID=widget_info(ev.top,FIND_BY_UNAME='tinc')
            widget_control,ID,get_value=tinc
            tinc = double(tinc)
            
            nx = round((emax-emin)/einc+1)
            energy = emin[0] + einc[0]*lindgen(nx)
            
            ny = round((tmax-tmin)/tinc+1)
            thickness = tmin[0] + tinc[0]*lindgen(ny)
            
            widget_control,ev.top,get_uvalue=list
            fraction=calcElasticScatteringFraction(list.ptr,energy,thickness,0,O=O)

            loadctrev,-3
            tvlct,v,/get
            p=surface(fraction*100,energy,thickness,title='Interaction with '+(*list.ptr).name,$
                xtitle='Energy (keV)',$
                ytitle='Thickness (mm)',$
                ztitle='Forward elastic scattering (% of I0)',$
                RGB_TABLE=v,VERT_COLORS=bytscl(fraction))
            p.Rotate,/reset
            
            m = max(fraction,indm)
            indx = indm mod nx
            indy = indm/nx

            mfrac = fraction[indx,indy]*100
            line=POLYLINE(energy[[indx,indx]],[min(thickness),max(thickness)],[mfrac,mfrac],/data,/current,color=!color.red)
            line=POLYLINE([min(energy),max(energy)],thickness[[indy,indy]],[mfrac,mfrac],/data,/current,color=!color.red)

            ID=widget_info(ev.top,FIND_BY_UNAME='energy')
            widget_control,ID,set_value=stringr(energy[indx])
            ID=widget_info(ev.top,FIND_BY_UNAME='thickness')
            widget_control,ID,set_value=stringr(thickness[indy])
            ID=widget_info(ev.top,FIND_BY_UNAME='fraction')
            widget_control,ID,set_value=stringr(mfrac)
            (*list.p).energy = energy[indx]
            (*list.p).thickness = thickness[indy]
            (*list.p).fraction = mfrac
 
            endcase
  else: return
  endcase
  endcase
3:begin ; text widget
  widget_control,ev.id,get_uvalue=uval
  widget_control,ev.id,get_value=val
  val=double(val)
  widget_control,ev.top,get_uvalue=list
  case uval of
  'energy': (*list.p).energy = val
  'thickness': (*list.p).thickness = val
  'emin': (*list.p).emin = val
  'emax': (*list.p).emax = val
  'einc': (*list.p).einc = val
  'tmin': (*list.p).tmin = val
  'tmax': (*list.p).tmax = val
  'tinc': (*list.p).tinc = val
  else:
  endcase
  endcase
endcase

end;pro ElasticScatteringFractionWidget_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ElasticScatteringFractionWidget,top,ptr,param

WIDGET_CONTROL, top, sensitive=0
GROUP_LEADER=top
modal=1

device,get_screen_size=screen
p = ptr_new(param)
base=widget_base(title='Elastic Scattering Fraction',/column,uvalue={top:top,ptr:ptr,p:p},$
                GROUP_LEADER=GROUP_LEADER,modal=modal,xoffset=screen[0]/2.,yoffset=screen[1]/2.)

xs=150

baset = widget_base(base,/row)
  w = widget_label(baset, value=' ', xsize=xs)
  w = widget_label(baset, value='begin', xsize=xs)
  w = widget_label(baset, value='end', xsize=xs)
  w = widget_label(baset, value='incr', xsize=xs)
  
baset = widget_base(base,/row)
  w = widget_label(baset, value='Energy (keV):', xsize=xs)
  w = widget_text(baset,/editable,value=stringr(param.emin),uname='emin',uvalue='emin',scr_xsize=xs)
  w = widget_text(baset,/editable,value=stringr(param.emax),uname='emax',uvalue='emax',scr_xsize=xs)
  w = widget_text(baset,/editable,value=stringr(param.einc),uname='einc',uvalue='einc',scr_xsize=xs)
  b = widget_button(baset,value='plot',uvalue='energy')

baset = widget_base(base,/row)
  w = widget_label(baset, value='Thickness (mm):', xsize=xs)
  w = widget_text(baset,/editable,value=stringr(param.tmin),uname='tmin',uvalue='tmin',scr_xsize=xs)
  w = widget_text(baset,/editable,value=stringr(param.tmax),uname='tmax',uvalue='tmax',scr_xsize=xs)
  w = widget_text(baset,/editable,value=stringr(param.tinc),uname='tinc',uvalue='tinc',scr_xsize=xs)
  b = widget_button(baset,value='plot',uvalue='thickness')

baset = widget_base(base,/row)
  w = widget_label(baset, value=' ', xsize=xs)
  w = widget_label(baset, value='keV', xsize=xs)
  w = widget_label(baset, value='mm', xsize=xs)
  w = widget_label(baset, value='fraction(%)', xsize=xs)

baset = widget_base(base,/row)
  w = widget_label(baset, value=' ', xsize=xs)
  w = widget_text(baset,/editable,value=stringr(param.energy),uname='energy',uvalue='energy',scr_xsize=xs)
  w = widget_text(baset,/editable,value=stringr(param.thickness),uname='thickness',uvalue='thickness',scr_xsize=xs)
  w = widget_text(baset,value=stringr(param.fraction),uname='fraction',uvalue='fraction',scr_xsize=xs)
  b = widget_button(baset,value='calc',uvalue='both')
 
baset = widget_base(base,/row,/NONEXCLUSIVE)
  b = widget_button(baset,value='Use refined parameters',uname='O',uvalue='O')
  
widget_control,base,/realize
Xmanager,'ElasticScatteringFractionWidget',base, event_handler='ElasticScatteringFractionWidget_Event',$
    cleanup='CleanUp_ElasticScatteringFractionWidget',GROUP_LEADER=GROUP_LEADER

param = *p
ptr_free,p
end;pro ElasticScatteringFractionWidget
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotElasticScatteringFraction,ptr,lambda,thickness,nextra,icurrent,O=O
; Analog to UnitCellFundamentals12

if keyword_set(O) then begin
    p=(*(*ptr).asu.O)[icurrent]
    psd=(*(*ptr).asu.SD)[icurrent]
    pg=(*(*ptr).globalparamIO.O)[icurrent]
    pgsd=(*(*ptr).globalparamIO.SD)[icurrent]
    bI=0b
endif else begin
    p=(*ptr).asu.I
    psd=(*(*ptr).asu.SD)[icurrent]
    pg=(*ptr).globalparamIO.I
    pgsd=(*(*ptr).globalparamIO.SD)[icurrent]
    bI=1b
endelse

n=n_elements(lambda)
UCMMass=0.
UCMu=fltarr(n)
UCCoh=UCMu
for j=0l,n-1 do begin
    UCMuj=0.
    UCCohj=0.
    
    ; Loop over positions in asym. unit
    for i=0l,(*ptr).nasupos-1 do begin
        k0=(*(*ptr).preverse_indices)[i]
        kn=(*(*ptr).preverse_indices)[i+1]-1
    
        ; Loop over substituted atoms on this ASU position
        m=((*p)[k0]).wyck.m ; multiplicity of position
        for k=k0,kn do begin
            mu=MassAtCoef(*(*p)[k].Z,*(*p)[k].Ion,*(*p)[k].Repeats,lambda[j],mmatom=mmatom)
            coh=CohCrossSection(*(*p)[k].Z,*(*p)[k].Ion,*(*p)[k].Repeats,lambda[j]);,mmatom=mmatom
    
            tmp=(*p)[k].SOF*mmatom*m
            if j eq 0 then UCMMass+=tmp
            UCMuj+=tmp*mu
            UCCohj+=tmp*coh
        endfor
    endfor
    
    UCMu[j]=UCMuj
    UCCoh[j]=UCCohj
endfor
UCMu/=UCMMass
UCCoh/=UCMMass

celparam=(*pg)[indCelParam(ptr)]
celparam_sd=double((*pgsd)[indCelParam(ptr)])
if bI then celparam_sd*=0
UCV=UCVolume(celparam,pder=UCV_sd,connect=celparamconnect(*(*ptr).pstrinfo))

ct=avogadroconstant(/noexp)
UCrho=UCMMass / ct / UCV

; Fraction (%) of elastic scattering detected, assuming all scattering is in the forward direction
energy=LambdaEnergy(lambda,/keV)
if n_elements(thickness) gt 1 then $
  frac=thickness*UCrho[0]*UCCoh[0]*exp(-UCMu[0]*thickness*UCrho[0])*100 $
else $
  frac=thickness[0]*UCrho*UCCoh*exp(-UCMu*thickness[0]*UCrho)*100

if n_elements(thickness) gt 1 then begin
  ; Print
  print,'Thickness(cm)     Forward elastic scattering (%)'
  print,transpose([[thickness[-nextra:*]],[frac[-nextra:*]]])
  frac=frac[0:(-nextra-1)]
  thickness=thickness[0:(-nextra-1)]
  
  ; Plot
  plot,thickness*10,frac,yrange=[1,2],/xs,/ylog,$
      xtitle='Thickness(mm)',ytitle='Forward elastic scattering (% of I0)',$
      title='Interaction with '+(*ptr).name+' at '+strtrim(energy,2)+' keV'
endif else begin
  ; Print
  print,'Energy(keV)     Forward elastic scattering (%)'
  print,transpose([[energy[-nextra:*]],[frac[-nextra:*]]])
  frac=frac[0:(-nextra-1)]
  energy=energy[0:(-nextra-1)]
  
  ; Plot
  plot,energy,frac,yrange=[1e-3,100],/xs,/ylog,$
      xtitle='Energy(keV)',ytitle='Forward elastic scattering (% of I0)',$
      title='Interaction with '+strtrim(round(thickness*1d4),2)+' micron '+(*ptr).name
endelse

end;pro plotElasticScatteringFraction
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotElasticScatteringFractionhelp,ptr,icurrent,O=O

if keyword_set(O) then return

; Wavelength
b=1d
e=100d
inc=0.1
n=round((e-b)/inc+1)
En=b+inc*lindgen(n)
En=[En,17.479,22.163] ; keV
nextra=2
n+=nextra
lambda=LambdaEnergy(En,/keV)

; Thickness
thickness=0.08 ; cm

; Plot
window,keyword_set(O)
plotElasticScatteringFraction,ptr,lambda,thickness,nextra,icurrent,O=O
end;pro plotElasticScatteringFractionhelp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotElasticScatteringFractionhelp2,ptr,icurrent,O=O

if keyword_set(O) then return

; Wavelength
En=21. ; keV
lambda=LambdaEnergy(En,/keV)

; Thickness (cm)
b=0.01
e=0.15
inc=0.001
n=round((e-b)/inc+1)
thickness=b+inc*lindgen(n) ; cm
thickness=[thickness,1e-3]
nextra=1
n+=nextra

; Plot
window,keyword_set(O)
plotElasticScatteringFraction,ptr,lambda,thickness,nextra,icurrent,O=O
end;pro plotElasticScatteringFractionhelp2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UnitCellFundamentals12,ptr,lambda,icurrent,O=O,pg=pg

; Attenuation corrected forward elastic scattered fraction
;plotElasticScatteringFractionhelp,ptr,icurrent,O=O ; as a function of energy
;plotElasticScatteringFractionhelp2,ptr,icurrent,O=O ; as a function of thickness

if (*ptr).nasupos eq 0 then return,{name:(*ptr).name,RietScl:0.,RietScl_sd:0d,UCMMass:0.,$
        UCMMass_sd:0d,UCMu:0.,UCMu_sd:0d,UCV:0.,UCV_sd:0d,UCrho:0.,UCrho_sd:0d,UCCoh:0d,UCCoh_sd:0d}

; Unit cell relative mass
;     UCMMass(g/mol) = Sum(asym unit) [mult * Sum(substit.)[SOF*MMatom(g/mol)] ]
;
; Unit cell absolute mass
;              UCMMass(g/mol)
;    Mass(g) = --------------
;                  Na(1/mol)
;
; Unit Cell mass attenuation coefficient
;                    Sum(asym unit) [mult * Sum(substit.)[SOF*MMatom(g/mol)*mu(cm^2/g)] ]
;      UCMu(cm^2/g) = ------------------------------------------------------------------
;                                             UCMMass(g/mol)
;
; Unit Cell volume (Angstrom^3)
;     UCV
;
; Unit cell density (g/cm^3)
;                        Mass(g)              MM(g/mol) . 10^24
;    UCrho(g/cm^3)= ----------------       =   -----------------
;                    V (Angstrom^3).10^(-24)  Na(1/mol) . V (Angstrom^3)

; UC data
if keyword_set(O) then begin
    p=(*(*ptr).asu.O)[icurrent]
    psd=(*(*ptr).asu.SD)[icurrent]
    pg=(*(*ptr).globalparamIO.O)[icurrent]
    pgsd=(*(*ptr).globalparamIO.SD)[icurrent]
    bI=0b
endif else begin
    p=(*ptr).asu.I
    psd=(*(*ptr).asu.SD)[icurrent]
    pg=(*ptr).globalparamIO.I
    pgsd=(*(*ptr).globalparamIO.SD)[icurrent]
    bI=1b
endelse

; Loop over positions in asym. unit
UCMMass=0.
UCMMass_sd=0d
UCMu=0.
UCMu_sd=0d
UCCoh=0.
UCCoh_sd=0d
for i=0l,(*ptr).nasupos-1 do begin
    k0=(*(*ptr).preverse_indices)[i]
    kn=(*(*ptr).preverse_indices)[i+1]-1

    ; Loop over substituted atoms on this ASU position
    m=((*p)[k0]).wyck.m ; multiplicity of position
    for k=k0,kn do begin
        mu=MassAtCoef(*(*p)[k].Z,*(*p)[k].Ion,*(*p)[k].Repeats,lambda,mmatom=mmatom)
        coh=CohCrossSection(*(*p)[k].Z,*(*p)[k].Ion,*(*p)[k].Repeats,lambda);,mmatom=mmatom

        tmp=(*p)[k].SOF*mmatom*m
        UCMMass+=tmp
        UCMu+=tmp*mu
        UCCoh+=tmp*coh
    
        tmp=((*psd)[k].SOF)^2.*mmatom*mmatom*m*m
        if bI then tmp*=0
        UCMMass_sd+=tmp
        UCMu_sd+=tmp*mu*mu
        UCCoh_sd+=tmp*coh*coh
    endfor
endfor
UCMMass_sd=sqrt(UCMMass_sd)
UCMu/=UCMMass
UCMu_sd=sqrt(UCMu_sd+ (UCMu*UCMMass_sd)^2.)/UCMMass
UCCoh/=UCMMass
UCCoh_sd=sqrt(UCCoh_sd+ (UCCoh*UCMMass_sd)^2.)/UCMMass

celparam=(*pg)[indCelParam(ptr)]
celparam_sd=double((*pgsd)[indCelParam(ptr)])
if bI then celparam_sd*=0
UCV=UCVolume(celparam,pder=UCV_sd,connect=celparamconnect(*(*ptr).pstrinfo))
UCV_sd=sqrt(total(UCV_sd*UCV_sd*celparam_sd*celparam_sd))

;ct=avogadroconstant()
;Mass=UCMMass / ct
;Mass_sd=UCMMass_sd / ct

ct=avogadroconstant(/noexp)
UCrho=UCMMass / ct / UCV
UCrho_sd=sqrt(UCMMass_sd*UCMMass_sd+(UCMMass/UCV*UCV_sd)^2.)/(ct*UCV)

S=(*pg)[*(*ptr).paramind[1,1],*]
S_sd=double((*pgsd)[*(*ptr).paramind[1,1],*])
if bI then S_sd*=0

;Y=ScatYield(UCMu,UCcoh,UCrho,100)
;print,Y

return,{name:(*ptr).name,RietScl:S,RietScl_sd:S_sd,UCMMass:UCMMass,UCMMass_sd:UCMMass_sd,$
    UCMu:UCMu,UCMu_sd:UCMu_sd,UCV:UCV,UCV_sd:UCV_sd,UCrho:UCrho,UCrho_sd:UCrho_sd,UCCoh:UCCoh,UCCoh_sd:UCCoh_sd}

end;function UnitCellFundamentals12
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UnitCellFundamentals11,ptr,icurrent,O=O,pg=pg

if keyword_set(O) then begin
    pg=(*(*ptr).globalparamIO.O)[icurrent]
    pgsd=(*(*ptr).globalparamIO.SD)[icurrent]
    bI=0b
endif else begin
    pg=(*ptr).globalparamIO.I
    pgsd=(*(*ptr).globalparamIO.SD)[icurrent]
    bI=1b
endelse

UCMMass=(*ptr).propfix[0]
UCMMass_sd=double((*ptr).propfix[1])
UCMu=(*ptr).propfix[2]
UCMu_sd=double((*ptr).propfix[3])
UCCoh=(*ptr).propfix[4]
UCCoh_sd=double((*ptr).propfix[5])

celparam=(*pg)[indCelParam(ptr)]
celparam_sd=double((*pgsd)[indCelParam(ptr)])
if bI then celparam_sd*=0
UCV=UCVolume(celparam,pder=UCV_sd,connect=celparamconnect(*(*ptr).pstrinfo))
UCV_sd=sqrt(total(UCV_sd*UCV_sd*celparam_sd*celparam_sd))

ct=avogadroconstant(/noexp)
UCrho=UCMMass / ct / UCV
UCrho_sd=sqrt(UCMMass_sd*UCMMass_sd+(UCMMass/UCV*UCV_sd)^2.)/(ct*UCV)

S=(*pg)[*(*ptr).paramind[1,1],*]
S_sd=double((*pgsd)[*(*ptr).paramind[1,1],*])
if bI then S_sd*=0

return,{name:(*ptr).name,RietScl:S,RietScl_sd:S_sd,UCMMass:UCMMass,UCMMass_sd:UCMMass_sd,$
    UCMu:UCMu,UCMu_sd:UCMu_sd,UCV:UCV,UCV_sd:UCV_sd,UCrho:UCrho,UCrho_sd:UCrho_sd,UCCoh:UCCoh,UCCoh_sd:UCCoh_sd}
end;function UnitCellFundamentals11
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UnitCellFundamentals10,ptr,icurrent,O=O,pg=pg

if keyword_set(O) then begin
    pg=(*(*ptr).globalparamIO.O)[icurrent]
    pgsd=(*(*ptr).globalparamIO.SD)[icurrent]
    bI=0b
endif else begin
    pg=(*ptr).globalparamIO.I
    pgsd=(*(*ptr).globalparamIO.SD)[icurrent]
    bI=1b
endelse

UCMMass=(*ptr).propfix[0]
UCMMass_sd=double((*ptr).propfix[1])
UCV=(*ptr).propfix[2]
UCV_sd=double((*ptr).propfix[3])
UCMu=(*ptr).propfix[4]
UCMu_sd=double((*ptr).propfix[5])
UCCoh=(*ptr).propfix[6]
UCCoh_sd=double((*ptr).propfix[7])

if (*ptr).codeproc[0,1] ne 0 then begin
    g=(*pg)[*(*ptr).paramind[0,1],*]
    g_sd=double((*pgsd)[*(*ptr).paramind[0,1],*])
    if bI then g_sd*=0
    g=g[1]
    g_sd=g_sd[1]
    UCV_sd=sqrt(g^6.*UCV_sd^2.+4*g^4.*UCV^2.*g_sd^2.)
    UCV*=g^3.
endif

ct=avogadroconstant(/noexp)
UCrho=UCMMass / ct / UCV
UCrho_sd=sqrt(UCMMass_sd*UCMMass_sd+(UCMMass/UCV*UCV_sd)^2.)/(ct*UCV)

S=(*pg)[*(*ptr).paramind[1,1],*]
S_sd=double((*pgsd)[*(*ptr).paramind[1,1],*])
if bI then S_sd*=0

return,{name:(*ptr).name,RietScl:S,RietScl_sd:S_sd,UCMMass:UCMMass,UCMMass_sd:UCMMass_sd,$
    UCMu:UCMu,UCMu_sd:UCMu_sd,UCV:UCV,UCV_sd:UCV_sd,UCrho:UCrho,UCrho_sd:UCrho_sd,UCCoh:UCCoh,UCCoh_sd:UCCoh_sd}
end;function UnitCellFundamentals10
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UnitCellFundamentals,ptr,lambda,icurrent,O=O,pg=pg
case (*ptr).type of
10: return,UnitCellFundamentals10(ptr,icurrent,O=O,pg=pg)
11: return,UnitCellFundamentals11(ptr,icurrent,O=O,pg=pg)
12: return,UnitCellFundamentals12(ptr,lambda,icurrent,O=O,pg=pg)
endcase
end;function UnitCellFundamentals
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ScalingFromWeight,X,ptrFit,K,lambda,O=O
; X: weight fractions
; ptrs: corresponding structure
;
; Intensity of a single Bragg peak
;    I(jH) = Sj . K(jH)
;    Sj = K . Xj / (Vj^2 . rho(j) . u)
;    where u = X1.u1+X2.u2+...+Xj.uj+...+Xn.un
;    
;    j: phase (for example calcite or iron or ...)
;    H: Miller index of the Bragg peak
;    I(jH): intensity of Bragg peak H from phase j
;    Sj: Rietveld scaling factor of phase j
;    K(jH): contains structure factor, LP factor etc.
;    K: constant depending on instrument
;    Vj: unit cell volume of phase j
;    rho(j): density of phase j
;    u: mass attenuation coefficient of the mixture of n phases
;    Xj: weight percentage of phase j in a mixture of n phases
;    uj: mass attenuation coefficient of phase j

; Sum(X)=1
X/=total(X)

; Get current V^2.rho.mu
V2rho=X*0.
mu=V2rho
istop=n_elements(X)
i=0l
ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if i eq istop then break
    if ~(*ptr).include then begin
        ptr=(*ptr).next
        continue
    endif
    
    info=UnitCellFundamentals(ptr,lambda,(*ptrFit).icurrent,O=O)
    V2rho[i]=info.UCV*info.UCV*info.UCrho
    mu[i]=info.UCMu
    
    ptr=(*ptr).next
    i++
endwhile

; Calculate scaling factors 
S=K[0]*X/(V2rho*total(X*mu))

; Change scaling factors
istop=n_elements(X)
i=0l
ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if i eq istop then break
    if ~(*ptr).include then begin
        ptr=(*ptr).next
        continue
    endif
    
    if keyword_set(O) then begin
        pg=(*(*ptr).globalparamIO.O)[(*ptrFit).icurrent]
    endif else begin
        pg=(*ptr).globalparamIO.I
    endelse

    (*pg)[*(*ptr).paramind[1,1],*]=S[i++]
    
    ptr=(*ptr).next
endwhile

end;pro ScalingFromWeight
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ExperimentalKfactor,ptrFit,lambda,O=O
; cfr. ScalingFromWeight
; K = Sj . Vj^2 . rho(j) . u / Xj = Iin . re^2 . lamda^3 . Sbeam / 8pi

ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if ~(*ptr).include then begin
        ptr=(*ptr).next
        continue
    endif
    
    ; Weight fraction of this phase
    if keyword_set(O) then X_i=(*(*(*ptr).propO)[(*ptrFit).icurrent])[0] $
    else X_i=(*(*ptr).propI)[0]
    
    ; Unit cell volume, mass att. coeff. and density of the phase
    info=UnitCellFundamentals(ptr,lambda,(*ptrFit).icurrent,O=O)
    mu_i=info.UCMu
    
    ; Scaling factor of this phase
    if keyword_set(O) then begin
        pg=(*(*ptr).globalparamIO.O)[(*ptrFit).icurrent]
    endif else begin
        pg=(*ptr).globalparamIO.I
    endelse
    SV2rho_i=info.UCV*info.UCV*info.UCrho*(*pg)[*(*ptr).paramind[1,1],*]
    
    ; Append
    if n_elements(X) eq 0 then begin
        X=X_i
        SV2rho=SV2rho_i
        mu=mu_i
    endif else begin
        X=[X,X_i]
        SV2rho=[SV2rho,SV2rho_i]
        mu=[mu,mu_i]
    endelse
    
    ptr=(*ptr).next
endwhile

return,SV2rho*total(X*mu)/X
end;function ExperimentalKfactor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkPhysicalPhaseInfo,ptr
b=(*ptr).include

case (*ptr).type of
10:    b and= total((*ptr).propfix[[0,2,4]] ne 0,/pres) eq 3 $
            and (*ptr).codeproc[1,1] ne 0
11: b and= total((*ptr).propfix[[0,2]] ne 0,/pres) eq 2 $
            and (*ptr).codeproc[1,1] ne 0
12:
endcase
return,b
end;function checkPhysicalPhaseInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SetPhysicalPhaseInfo,ptrFit,lambda,O=O
n=0
ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if checkPhysicalPhaseInfo(ptr) then begin
        tmp=UnitCellFundamentals(ptr,lambda,(*ptrFit).icurrent,O=O)
        if n_elements(info)    eq 0 then begin
            info=tmp
            ptrkeep=ptr
        endif else begin
            info=[info,tmp]
            ptrkeep=[ptrkeep,ptr]
        endelse
    endif
    ptr=(*ptr).next
endwhile

n=n_elements(info)
if n ne 0 then begin
    S=reform(info.RietScl)>0
    S_sd=reform(info.RietScl_sd)
    
    ; Weight fractions
    SMV=info.UCMMass*info.UCV*S
    SMV_sd=sqrt((info.UCMMass_sd*info.UCV*S)^2.+(info.UCMMass*info.UCV_sd*S)^2.+(info.UCMMass*info.UCV*S_sd)^2.)
    m=total(SMV)/100
    m_sd=sqrt(total(SMV_sd*SMV_sd))/100
    if m ne 0 then begin
        SMV/=m
        SMV_sd=sqrt(SMV_sd*SMV_sd+(SMV*m_sd)^2.)/m
    endif else begin
        SMV[*]=0
        SMV_sd[*]=0
    endelse
    
    ; Volume fractions
    SVV=info.UCV*info.UCV*S
    SVV_sd=info.UCV*sqrt((2*S*info.UCV_sd)^2.+(info.UCV*S_sd)^2.)
    m=total(SVV)/100
    m_sd=sqrt(total(SVV_sd*SVV_sd))/100
    if m ne 0 then begin
        SVV/=m
        SVV_sd=sqrt(SVV_sd*SVV_sd+(SVV*m_sd)^2.)/m
    endif else begin
        SVV[*]=0
        SVV_sd[*]=0
    endelse
    
    ; Mole fractions
    ; ni/n = wi/MMi/sum(wj/MMj) = Si.Vi/sum(Sj.Vj)
    SV=info.UCV*S
    SV_sd=sqrt((S*info.UCV_sd)^2.+(info.UCV*S_sd)^2.)
    m=total(SV)/100
    m_sd=sqrt(total(SV_sd*SV_sd))/100
    if m ne 0 then begin
        SV/=m
        SV_sd=sqrt(SV_sd*SV_sd+(SV*m_sd)^2.)/m
    endif else begin
        SV[*]=0
        SV_sd[*]=0
    endelse
    
    ; Construct table of fundamental
    data=[[SMV],[SMV_sd],[SVV],[SVV_sd],[SV],[SV_sd],[info.UCMMass],[info.UCMMass_sd],$
        [info.UCV],[info.UCV_sd],[info.UCrho],[info.UCrho_sd],[info.UCMu],[info.UCMu_sd],[info.UCCoh],[info.UCCoh_sd]]
    
    if keyword_set(O) then $
        for i=0l,n-1 do *(*(*ptrkeep[i]).propO)[(*ptrFit).icurrent]=reform(float(data[i,*])) $
    else $
        for i=0l,n-1 do *(*ptrkeep[i]).propI=reform(float(data[i,*]))
    
endif

end;pro SetPhysicalPhaseInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetPhysicalPhaseInfo,ptrFit,n,O=O
n=0
ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if checkPhysicalPhaseInfo(ptr) then begin
        
        if n_elements(data)    eq 0 then begin
            if keyword_set(O) then data=*(*(*ptr).propO )[(*ptrFit).icurrent] $
            else data=*(*ptr).propI

            colname=(*ptr).name
            rowname=(*ptr).propname
        endif else begin
            if keyword_set(O) then data=[[data],[*(*(*ptr).propO )[(*ptrFit).icurrent]]] $
            else data=[[data],[*(*ptr).propI]]
            
            colname=[colname,(*ptr).name]
        endelse
    endif
    ptr=(*ptr).next
endwhile

n=n_elements(colname)
if n ne 0 then begin
    return,{data:transpose(double(data)),$
        colname:colname,$
        rowname:rowname}
endif else return,0b
end;function GetPhysicalPhaseInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetNumberPhysicalPhaseInfo,ptrFit
n=0
nmult=0
ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if checkPhysicalPhaseInfo(ptr) then begin
        if n eq 0 then nmult=n_elements((*ptr).propname)
        n++
    endif
    ptr=(*ptr).next
endwhile

return,n*nmult
end;function GetNumberPhysicalPhaseInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CalcStoichiometry,ptr,O=O
str=''
if (*ptr).type ne 12 then return,str
if (*ptr).nasuatoms eq 0 then return,str

icurrent=0
if keyword_set(O) then p=(*(*ptr).asu.O)[icurrent] $
else p=(*ptr).asu.I

Z=''
factors=0.
charge=0
UCMMass=0.
w=0.
for i=0l,(*ptr).nasuatoms-1 do begin
    Z_=*(*p)[i].Z
    ion=*(*p)[i].Ion
    repeats=*(*p)[i].Repeats
    for j=0l,n_elements(Z_)-1 do begin
        str=ZElement(Z_[j],ion=ion[j])
        f=repeats[j]*(*p)[i].wyck.m*(*p)[i].SOF
        
        m=AtMass(Z_[j])*f
        UCMMass+=m
    
        ind=where(Z eq str,ct)
        if ct eq 0 then begin
            Z=[Z,str]
            factors=[factors,f]
            w=[w,m]
        endif else begin
            factors[ind[0]]+=f
            w[ind[0]]+=m
        endelse
        
        charge+=ion[j]*f
    endfor
endfor

nZ=n_elements(Z)-1
if nZ eq 0 then return,''
Z=Z[1:*]
factors=factors[1:*]
w=w[1:*]
w/=total(w)
w*=100

str=['Charge: '+string(charge,format='(f0.2)'),$
    'Stoichiometry: ',$
    '   '+Z+': '+string(factors,format='(f0.2)'),$
    'Weight%: ',$
    '   '+Z+': '+string(w,format='(f0.1)')+' (%)']

return,str
end;function CalcStoichiometry
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
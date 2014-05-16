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

pro updatefilecounter,ev,list
ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
nfiles=widget_info(ID,/list_number)
ID=widget_info(ev.top,FIND_BY_UNAME='FilesText')
widget_control,ID,set_value=string([nfiles,list.nfiles,list.nrtif],format='(I0," files (progress: ",I0,"/",I0,")")')
end;pro updatefilecounter
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CalcNrGroupsFAST,list

ptrCur=(*(*list.ptrModel).ptrSum).next
if ~ptr_valid(ptrCur) then return,0
Profilex=(*ptrCur).sum
Profilex=(*Profilex)[*,0]

tmp=BraggDtoX((*ptrCur).mask[3:4,0],list.lambda,/onlyt,/angledeg)
t=min(abs(Profilex-tmp[0]),x1);convert from two-theta to index in Profilex
t=min(abs(Profilex-tmp[1]),x2)

t=round(0>[x1<x2,x1>x2]<(n_elements(Profilex)-1))

return,t[1]-t[0]+1
end;function CalcNrGroupsFAST
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CalcNrGroups,ptrModel

pfit=(*ptrModel).ptrFit
nrgroups=TotalTypenParam(pfit,/batch)
nrgroups+=total((*(*ptrModel).ptrSum).n)

return,fix(nrgroups)
end;function CalcNrGroups
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Fillblock,block,off,ptrFit

s=DimSize(block,3)
dim0=s[0]
dim1=s[1]
dim2=s[2]

nfilesdone=(*ptrFit).icurrent<(dim1*dim2)
if nfilesdone le 0 then return

; Add refined parameters
ptr=(*ptrFit).next
off0=off
while ptr_valid(ptr) do begin

    for i=0l,nfilesdone-1 do begin ; loop over files
        b=InitFitModelShrink(ptr,i,A,nparam)
        if b then begin
            i0=i mod dim1
            i1=i/dim1
            block[off:off+nparam-1,i0,i1]=A
        endif
    endfor

    off+=nparam
    ptr=(*ptr).next
endwhile

; Add fit parameters
; Make sure FitModelLabels and TotalTypenParam are changed accordingly
if off ne off0 then begin
    for i=0l,nfilesdone-1 do begin ; loop over files
        A=(*(*(*ptrFit).factors)[i])
        offadd=n_elements(A)-1
        i0=i mod dim1
        i1=i/dim1
        block[off:off+offadd,i0,i1]=A
    endfor
endif

end;pro Fillblock
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MakeDatablock,list,block,dim0,dim1,dim2,transform=transform,FASTMODE=FASTMODE,statgroupi=statgroupi
; Datablock:
;    dim0: number of groups
;    dim1: first map/scan dimension
;    dim2: second map/scan dimension

; Default orientation
dim0=CalcNrGroups(list.ptrModel)
case list.SCANDIM of
1:    begin
    dim1=list.N
    dim2=1
    endcase
2:    begin
    dim1=list.CR[0]
    dim2=list.CR[1]
    endcase
endcase
nfiles=dim1*dim2 ; same as list.nrtif

; make datablock
if keyword_set(FASTMODE) then begin
    ; datablock not [dim0,dim1,dim2] but [dim0,dim1*dim2]
    ptrCur=(*(*list.ptrModel).ptrSum).next
    if list.AzIntInfo.calcerror then block=[[(*(*ptrCur).sum)],[(*(*ptrCur).sum2)[*,1:*]]] $ ; [[x],[spe1],[spe2],...,[errorspe1],[errorspe2]]
    else block=(*(*ptrCur).sum)
    
    dim0=(size(block,/dim))[0]
    
    nfilesadd=(list.nFiles-nfiles)>0
    if nfilesadd ne 0 then $
        block=[[block],[fltarr(dim0,nfilesadd)]]

    transform=0
endif else begin
    block=fltarr(dim0,dim1,dim2)
    
    ; Fill datablock with ROI's
    ns=total((*(*list.ptrModel).ptrSum).n)
    if ns ne 0 then begin
        ptrCur=(*(*list.ptrModel).ptrSum).next
        keep=fltarr(dim1,dim2) ; use this to handle size issues
        nfiles <= n_elements(*(*ptrCur).sum)
        for i=0l,ns-1 do begin
            keep[0]=(*(*ptrCur).sum)[0:nfiles-1]
            block[i,*,*]=keep
            if list.NORMROIs then begin
                ma=max(block[i,*,*],min=mi)
                block[i,*,*]=(block[i,*,*]-mi)/(ma-mi)
            endif
            ptrCur=(*ptrCur).next
        endfor
    endif
    
    ; Fill datablock with fit results
    nf=total((*(*list.ptrModel).ptrFit).n)
    if nf ne 0 then Fillblock,block,ns,(*list.ptrModel).ptrFit
    
    if keyword_set(statgroupi) then begin
        if list.nfiles eq 0 then begin
            statgroupi.may=0
            statgroupi.miy=0
        endif else begin
            may=max((block[0>list.groupi<(dim0-1),*,*])[0:(nfiles<list.nfiles)-1],min=miy)
            statgroupi.may=may
            statgroupi.miy=miy
        endelse
    endif
endelse

if keyword_set(transform) then ReformXDIDatablock,block,[0b,list.TRAN],[0b,list.REV1],[0b,list.REV2],[0b,list.REV3],0b,dim0,dim1,dim2

end;pro MakeDatablock
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Save1Dxdiheader,lun,list

printf,lun,'$HEADING:'
case list.unit1 of
0: begin
    unit=msymbols('micro')+'m'
    mul=1000.
    endcase
1: begin
    mul=1
    unit='deg'
    endcase
endcase

fi='I'+stringr(max(strlen(stringr([list.N,list.ct])))+1)
if list.TRAN then begin
    printf,lun,format='("dscan samv ",F10.3,F10.3,'+fi+','+fi+')',$
       list.coord1,list.N-1,list.ct
    printf,lun,format='("size: ", '+fi+'," pixels")',list.N
    dif=((list.coord1[0,0]>list.coord1[1,0])-(list.coord1[0,0]<list.coord1[1,0]))*mul
    printf,lun,format='("vertical: ",F10.3," '+unit+' in steps of ",F10.3," '+unit+'")',dif,dif/(list.N-1)
endif else begin
    printf,lun,format='("dscan samh ",F10.3,F10.3,'+fi+','+fi+')',$
       list.coord1,list.N-1,list.ct
    printf,lun,format='("size: ", '+fi+'," pixels")',list.N
    dif=((list.coord1[0,0]>list.coord1[1,0])-(list.coord1[0,0]<list.coord1[1,0]))*mul
    printf,lun,format='("horizontal: ",F10.3," '+unit+' in steps of ",F10.3," '+unit+'")',dif,dif/(list.N-1)
endelse

end;pro Save1Dxdiheader
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Save2Dxdiheader,lun,list

printf,lun,'$HEADING:'
unit=strarr(2)
mul=fltarr(2)
for i=0l,1 do begin
    case list.unit[i] of
    0: begin
        unit[i]=msymbols('micro')+'m'
        mul[i]=1000.
        endcase
    1: begin
        mul[i]=1
        unit[i]='deg'
        endcase
    endcase
endfor
if list.REV3 then mesh='meshc' else mesh='mesh'

fi='I'+stringr(max(strlen(stringr([list.CR,list.ct])))+1)
if list.TRAN then begin
    printf,lun,format='("'+mesh+' samv ",F10.3,F10.3,'+fi+'," samh ",F10.3,F10.3,'+fi+','+fi+')',$
       list.coord[*,0],list.CR[0]-1,list.coord[*,1],list.CR[1]-1,list.ct
    printf,lun,format='("size: ", '+fi+', " x ", '+fi+', " pixels")',reverse(list.CR)
    dif=((list.coord[0,0]>list.coord[1,0])-(list.coord[0,0]<list.coord[1,0]))*mul[0]
    printf,lun,format='("vertical: ",F10.3," '+unit[0]+' in steps of ",F10.3," '+unit[0]+'")',dif,dif/(list.CR[0]-1)
    dif=((list.coord[0,1]>list.coord[1,1])-(list.coord[0,1]<list.coord[1,1]))*mul[1]
    printf,lun,format='("horizontal: ",F10.3," '+unit[1]+' in steps of ",F10.3," '+unit[1]+'")',dif,dif/(list.CR[1]-1)
endif else begin
    printf,lun,format='("'+mesh+' samh ",F10.3,F10.3,'+fi+'," samv ",F10.3,F10.3,'+fi+','+fi+')',$
       list.coord[*,0],list.CR[0]-1,list.coord[*,1],list.CR[1]-1,list.ct
    printf,lun,format='("size: ", '+fi+', " x ", '+fi+', " pixels")',list.CR
    dif=((list.coord[0,0]>list.coord[1,0])-(list.coord[0,0]<list.coord[1,0]))*mul[0]
    printf,lun,format='("horizontal: ",F10.3," '+unit[0]+' in steps of ",F10.3," '+unit[0]+'")',dif,dif/(list.CR[0]-1)
    dif=((list.coord[0,1]>list.coord[1,1])-(list.coord[0,1]<list.coord[1,1]))*mul[1]
    printf,lun,format='("vertical: ",F10.3," '+unit[1]+' in steps of ",F10.3," '+unit[1]+'")',dif,dif/(list.CR[1]-1)
endelse

end;pro Save2Dxdiheader
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Makefiles,list
if list.nrtif eq 0 then return,''
files=replicate('Dummy',list.nrtif)

if ptr_valid(list.PTR_files) then begin
    nFiles=n_elements(*list.PTR_files)<list.nrtif
    files[0]=(*list.PTR_files)[0:nFiles-1]
endif

return,files
end;function Makefiles
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SaveXDI,list,xdifile,nodefaultorient=nodefaultorient

if openw_safe(lun, xdifile) then return,0b
MakeDatablock,list,block,dim0,dim1,dim2,transform=~keyword_set(nodefaultorient)

case list.SCANDIM of
1:    Save1Dxdiheader,lun,list
2:    Save2Dxdiheader,lun,list
endcase

SaveModelResult,lun,list

printf,lun,'$IMS:'
printf,lun,3
printf,lun,dim0,dim1,dim2
writeu,lun,block

printf,lun,''
printf,lun,'$FILES:'
; Default orientation: TRAN=0,REV1=0,...
printf,lun,Makefiles(list),format='(A)'

free_lun, lun

return,1b
end;function SaveXDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CornersAzimuth,aziarr,azi0,azi1

; aziarr in ]-180,180]
; azi0 map to ]-180,180]
; azi1 = azi0 + |xazi| (left turn)
;       = azi0 - |xazi| (right turn)

; test angles:
;azi0=[45,45,45,45,45,45]
;azi1=[45,135,225,315,370,405]
;azi0=[azi0,45,45,45,45,45,45,45]
;azi1=[azi1,45,10,-45,-135,-225,-300,-315]
;azi0=[azi0,-135,-135,-135,-135,-135,-135]
;azi1=[azi1,-135,-45,45,135,190,225]
;azi0=[azi0,-135,-135,-135,-135,-135,-135,-135]
;azi1=[azi1,-135,-170,-255,-315,-405,-480,-495]

; convert right turning to left turning
t0=azi0<azi1
t1=azi0>azi1

; return when the whole angular range is covered
if (t1-t0) ge 2*!dpi then return,1b

; map azi0 in ]-180,180] and change azi1 accordingly
if (azi0 lt -!dpi) or (azi0 gt !dpi) then begin
    diff = azi1-azi0
    a0=azi0 mod (2*!dpi) ; map in ]-360,360[
    if a0 le -!dpi then a0+=2*!dpi ; map in ]-180,360[
    if a0 gt !dpi then a0-=2*!dpi ; map in ]-180,180]
    a1=a0+diff
    t0=a0<a1
    t1=a0>a1
    bleft=t0 eq a0
endif else bleft=t0 eq azi0

; everything between t0 and t1 set to1
if bleft then begin ;left turning
    b = aziarr ge t0
    if t1 lt !dpi then b and= aziarr le t1 $
    else b or= aziarr le (t1-2*!dpi)
endif else begin ;right turning
    b = aziarr le t1
    if t0 gt -!dpi then b and= aziarr ge t0 $
    else b or= aziarr ge (t0+2*!dpi)
endelse

return,b
end;function CornersAzimuth
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CornersMask, mask, list, imgmult=imgmult
bFast=keyword_set(imgmult)

case mask[0] of
1:  begin
    out=mask[3:6]
    if bFast then begin
       col=mask[4]-mask[3]+1
       row=mask[6]-mask[5]+1
       imgmult=lindgen(col*row)
    endif
    return,out
    endcase
2:  begin
    tt=BraggDtoX(mask[[3,5]],list.lambda,/onlyt)
    azi0=mask[4]
    azi1=mask[6]

    ; Calc phi range
    phi=azDetRange(azi0,azi1,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

    ; two circles
    xEllipse=EllipsePlotCoord(tt[0],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center)
    xli=min(xEllipse[*,0],max=xla)
    yli=min(xEllipse[*,1],max=yla)
    xEllipse=EllipsePlotCoord(tt[1],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center)
    xli<=min(xEllipse[*,0],max=tmp)
    xla>=tmp
    yli<=min(xEllipse[*,1],max=tmp)
    yla>=tmp

    endcase
3:  begin
    res=BraggDtoX(mask[[3,4]],list.lambda,/onlyt)
    tt=[res[0]<res[1],res[0]>res[1]]
    if tt[0] eq 0 then tt=tt[1]
    azi0=0
    azi1=2*!pi

    ; Calc phi range
    phi=azDetRange(azi0,azi1,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

    ; two circles
    xEllipse=EllipsePlotCoord(tt[0],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center)
    xli=min(xEllipse[*,0],max=xla)
    yli=min(xEllipse[*,1],max=yla)
    if n_elements(tt) ne 1 then begin
        xEllipse=EllipsePlotCoord(tt[1],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center)
        xli<=min(xEllipse[*,0],max=tmp)
        xla>=tmp
        yli<=min(xEllipse[*,1],max=tmp)
        yla>=tmp
    endif else begin
        xli<=0
        xla>=0
        yli<=0
        yla>=0
    endelse

    endcase
endcase

xli=floor(xli)
yli=floor(yli)
xla=ceil(xla)
yla=ceil(yla)
out=[xli,xla,yli,yla]

if bFast then begin
    col=(xla-xli+1)
    row=(yla-yli+1)

    ; ----make x,y matrices----
    xrange=indgen(col)+xli
    yrange=indgen(row)+yli
    xarr=(xrange)#replicate(1,row)
    yarr=replicate(1,col)#(yrange)

    ; ----make tt,azi matrices----
    ttarr=BraggXYtoX(xarr,yarr,0.,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/onlyt,aziout=aziarr)
    
    ; ----indices to be used----
    b=ttarr ge (tt[0]<tt[1])
    b AND= ttarr le (tt[0]>tt[1])
    b AND= CornersAzimuth(azMmtoPixel(aziarr,list.scrx,list.scry),azi0,azi1)
    ind=where(b,ct)
    if (ct ne 0) then imgmult=ind else imgmult=-1L
endif

return,out
end;function CornersMask
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ROIlabels,list,labels

ret=0b
labels=''

ptrCur=(*(*list.ptrModel).ptrSum).next
while PTR_VALID(ptrCur) do begin
    case (*ptrCur).type of
    0:    begin
        for i=0l,(*ptrCur).masknr-1 do begin
            case (*ptrCur).mask[0,i] of
            1:    begin
                type='square'
                XY=(*ptrCur).mask[3:*,i]
                res=BraggXYtoX(XY[[0,0,1,1]],XY[[2,3,2,3]],list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/onlyt,/onlyd)
                tthetaB=min(res[*,0],max=tthetaE)
                dspacB=min(res[*,1],max=dspacE)
                endcase
            2:    begin
                type='arc   '
                dspacB=(*ptrCur).mask[3,i]<(*ptrCur).mask[5,i]
                dspacE=(*ptrCur).mask[3,i]>(*ptrCur).mask[5,i]
                res=BraggDtoX([dspacB,dspacE],list.lambda,/onlyt)
                tthetaB=res[1]
                tthetaE=res[0]
                endcase
            3:    begin
                type='ellipse'
                dspacB=(*ptrCur).mask[3,i]<(*ptrCur).mask[4,i]
                dspacE=(*ptrCur).mask[3,i]>(*ptrCur).mask[4,i]
                res=BraggDtoX([dspacB,dspacE],list.lambda,/onlyt)
                tthetaB=res[1]
                tthetaE=res[0]
                endcase
            endcase
            tthetaA=(tthetaB+tthetaE)/2
            dspacA=(dspacB+dspacE)/2

            str=string(format='(I5,I5," '+type+' ",F10.4,F10.4,F10.4,F10.4,F10.4,F10.4)',$
                (*ptrCur).group,(*ptrCur).mask[1,i],dspacB,dspacE,dspacA,tthetaB*180/!dpi,tthetaE*180/!dpi,tthetaA*180/!dpi)
            labels=[labels,str]
        endfor
        endcase

    1:    begin
        for i=0l,(*ptrCur).masknr-1 do begin
            dspacB=(*ptrCur).mask[3,i]<(*ptrCur).mask[4,i]
            dspacE=(*ptrCur).mask[3,i]>(*ptrCur).mask[4,i]
            res=BraggDtoX([dspacB,dspacE],list.lambda,/onlyt)
            tthetaB=res[1]
            tthetaE=res[0]
            tthetaA=(tthetaB+tthetaE)/2
            dspacA=(dspacB+dspacE)/2
            str=string(format='(I5,I5," 1D ",F10.4,F10.4,F10.4,F10.4,F10.4,F10.4)',$
                (*ptrCur).group,(*ptrCur).mask[1,i],dspacB,dspacE,dspacA,tthetaB*180/!dpi,tthetaE*180/!dpi,tthetaA*180/!dpi)
            labels=[labels,str]
        endfor
        endcase

    2:    begin
        str=string(format='(I5," -  restimage - - - - - -")',(*ptrCur).group)
        labels=[labels,str]
        endcase
    3:    begin
        str=string(format='(I5," ",A,"  ",A," - - - - - -")',(*ptrCur).group,(*ptrCur).name,(*ptrCur).name)
        labels=[labels,str]
        endcase
    endcase
    ptrCur=(*ptrCur).next
endwhile

if n_elements(labels) gt 1 then begin
    labels=labels[1:*]
    ret=1b
endif

return, ret
end;function ROIlabels
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function getcurrentlabel,list

case list.FASTMODE of
1: return,(list.groupi le 0)?'Selected ROI: sum':string(list.groupi,format='("Selected ROI: PC",I0)')
2: return,''
else:
endcase

i=list.groupi

b=ROIlabels(list,labels)
if b ne 0 then begin
    n=n_elements(labels)
    n=long((strsplit(labels[n-1],' ',/extract))[0])+1
    if i lt n then begin
        ind=where(strpos(strtrim(labels,2),stringr(i)) eq 0,ct)
        if ct ne 0 then ret=((ct gt 1)?'(Group)':'')+labels[ind[0]]
    endif
endif else n=0

if n_elements(ret) eq 0 then begin
    i-=n
    b=FitModelLabels((*list.ptrModel).ptrFit,labels)
    if b ne 0 then begin
        n=n_elements(labels)
        if i lt n then ret=labels[i]
    endif else n=0
endif

if n_elements(ret) eq 0 then $
    ret=string(i,format='("Group ",I0)') $
else begin
    ret=strsplit(ret,' ',/extract,count=ct)
    ct=string(ct,format='(I0)')
    ret=string(ret,format='('+ct+'(A," | "))')
endelse

return,ret
end;function getcurrentlabel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SaveModelResult,lun,list

printf,lun,'$ORIENTATION:'
printf,lun,[list.TRAN,list.REV1,list.REV2,list.REV3]
printf,lun,list.Hsign,list.Vsign
printf,lun,list.ScanDim

printf,lun,'$PHYSMASK:'
printf,lun,(*list.ptrModel).n
printf,lun,'Qualitative analysis:'
printf,lun,(*(*list.ptrModel).ptrSum).n
printf,lun,'Group  Area  Type  Dmin('+msymbols('angstroms')+')  Dmax('+msymbols('angstroms')+')  Davg('+msymbols('angstroms')+')  ttmin('+msymbols('degrees')+')  ttmax('+msymbols('degrees')+')  ttavg('+msymbols('degrees')+')'
b=ROIlabels(list,labels)
if b then $
    printf,lun,labels,format='(A)'
printf,lun,'Quantitative analysis:'
printf,lun,(*(*list.ptrModel).ptrFit).n
printf,lun,TotalTypenParam((*list.ptrModel).ptrFit,/batch)
printf,lun,'Model  Symbol  Type  Init'
b=FitModelLabels((*list.ptrModel).ptrFit,labels)
if b then $
    printf,lun,labels,format='(A)'

end;pro SaveModelResult
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FASTMODE2GetSize,list,xsize,ysize
ptrCur=(*(*list.ptrModel).ptrSum).next
if not ptr_valid(ptrCur) then return,0b
    
s=size(*(*ptrCur).sum)
case list.SCANDIM of
1:    begin
    if (list.N lt 1) then begin
       printw,list.prid,'Unlikely dimensions of linescan!'
       return,0B
    endif
    list.xsx=list.N
    list.ysy=1
    endcase
2:    begin
    if (list.CR[0] lt 1) or (list.CR[1] lt 1) then begin
       printw,list.prid,'Unlikely dimensions of map!'
       return,0B
    endif
    list.xsx=list.CR[0]
    list.ysy=list.CR[1]
    endcase
endcase

; xsx;ysy: original size (HxV)
if list.TRAN then begin
    xsize=list.ysy*s[1]
    ysize=list.xsx*s[2]
endif else begin
    xsize=list.xsx*s[1]
    ysize=list.ysy*s[2]
endelse

return,1b
end;function FASTMODE2GetSize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PrepOut,list,ev
T=systime(1)

widget_control,ev.top,update=0

if list.FASTMODE eq 1 then list.groupi=0>list.groupi<(CalcNrGroupsFAST(list)-1) $
else list.groupi=0>list.groupi<(CalcNrGroups(list.ptrModel)-1)
list.plottitle=getcurrentlabel(list)

if widget_info(list.outdraw,/valid) then $
    widget_control,list.outdraw,/destroy

if widget_info(list.outdrawsum,/valid) then begin
    widget_control,list.outdrawsum,/destroy
    list.outdrawsum=0
    if winidvalid(list.pixindexsum) then wdelete,list.pixindexsum
endif

if list.FASTMODE eq 2 then begin
    if ~FASTMODE2GetSize(list,xsize,ysize) then begin
        widget_control,ev.top,update=1
        return,0b
    endif
    
    xsize*=list.nmag
    ysize*=list.nmag
    
    uname='DrawBase'
    ID=widget_info(ev.top,FIND_BY_UNAME=uname)

    draw=widget_draw(ID,xsize=xsize,ysize=ysize,/scroll,$
        x_scroll_size=list.xsm<(xsize+1),y_scroll_size=list.ysm<(ysize+1),$
        Event_Pro='FILESELECT_PROCESS_EVENTS',$
        KEYBOARD_EVENTS=2,/BUTTON_EVENTS,/MOTION_EVENTS)
    widget_control,draw,get_value=drawid
    list.outdrawid=drawid
    list.outdraw=draw
    list.title='Grid: Processing...'
    widget_control,ev.top,TLB_SET_TITLE=list.title+(list.char)[list.chari]
    printw,list.prid,'Make Plotting window: '+string(systime(1)-T)+' sec'

    widget_control,ev.top,update=1
    return,1b
endif

case list.SCANDIM of
1:    begin
    if (list.N lt 1) then begin
       printw,list.prid,'Unlikely dimensions of linescan!'
       widget_control,ev.top,update=1
       return,0B
    endif
    if list.PLOT then begin
        if list.FASTMODE eq 1 then uname='DrawBase2' else uname='DrawBase'
        ID=widget_info(ev.top,FIND_BY_UNAME=uname)
        draw=widget_draw(ID,xsize=list.xsm,ysize=list.ysm,$
           Event_Pro='FILESELECT_PROCESS_EVENTS',$
           KEYBOARD_EVENTS=2,/BUTTON_EVENTS,/MOTION_EVENTS)
        widget_control,draw,get_value=drawid
        list.outdrawid=drawid
        list.outdraw=draw
        list.title='Line-scan: Processing...'
        widget_control,ev.top,TLB_SET_TITLE=list.title+(list.char)[list.chari]
        printw,list.prid,'Make Plotting window: '+string(systime(1)-T)+' sec'
    endif
    endcase
2:    begin
    if (list.CR[0] lt 1) or (list.CR[1] lt 1) then begin
       printw,list.prid,'Unlikely dimensions of map!'
       widget_control,ev.top,update=1
       return,0B
    endif
    if list.PLOT then begin
       if list.FASTMODE eq 1 then uname='DrawBase2' else uname='DrawBase'
       ID=widget_info(ev.top,FIND_BY_UNAME=uname)
       if list.TRAN then CR=reverse(list.CR) else CR=list.CR
       
       xsize=list.nmag*CR[0]
       ysize=list.nmag*CR[1]
       ID2=widget_info(ev.top,FIND_BY_UNAME='plottitle')
       if widget_info(ID2,/valid) then widget_control,ID2,/destroy
       label=widget_label(ID,value=list.plottitle,uname='plottitle')
       draw=widget_draw(ID,xsize=xsize,ysize=ysize,/scroll,$
         x_scroll_size=list.xsm<(xsize+1),y_scroll_size=list.ysm<(ysize+1),$
         Event_Pro='FILESELECT_PROCESS_EVENTS',$
         KEYBOARD_EVENTS=2,/BUTTON_EVENTS,/MOTION_EVENTS)
       widget_control,draw,get_value=drawid
       list.outdrawid=drawid
       list.outdraw=draw
       list.title='Map: Processing...'
       widget_control,ev.top,TLB_SET_TITLE=list.title+(list.char)[list.chari]
       printw,list.prid,'Make Plotting window: '+string(systime(1)-T)+' sec'
    endif
    endcase
endcase

if list.PLOT then $
if list.FASTMODE eq 1 then begin
    ID=widget_info(ev.top,FIND_BY_UNAME='DrawBase')
    draw=widget_draw(ID,xsize=list.xsm,ysize=list.ysm,/BUTTON_EVENTS,/MOTION_EVENTS,$
        KEYBOARD_EVENTS=2,Event_Pro='ROIFASTMODE_PROCESS_EVENTS')
    widget_control,draw,get_value=drawid
    list.outdrawidsum=drawid
    list.outdrawsum=draw

    Window, /Free, xsize=list.xsm, ysize=list.ysm, /Pixmap
    list.pixindexsum = !D.Window
endif

widget_control,ev.top,update=1
return,1B
end;function PrepOut
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro resetmag_xrdbp,list
xsm=list.xsm/2
ysm=list.ysm/2
if list.SCANDIM eq 2 then begin
    list.nmag=((xsm/float(list.CR[~ list.TRAN]))<(ysm/float(list.CR[list.TRAN])))>1
endif
if list.FASTMODE eq 2 then begin
    if FASTMODE2GetSize(list,xsize,ysize) then begin
        list.nmag=((xsm/float(xsize))<(ysm/float(ysize)))>1
    endif
endif
if list.nmag eq 0 then list.nmag=1
end;pro resetmag_xrdbp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BatchInit,list,ev

; ----Destroy existing dynamic list:
; When pressed 'Go' => i.e. new scan----
ptr_free,list.PTR_files
list.PTR_files=PTR_NEW()
list.nFiles=0
list.AzIntParam[*]=0
DestroyModel,list.ptrModel
heap_free,list.Subfiles
heap_free,list.WarpStruc
list.Subfiles.nleft=0
list.Subfiles.n=0
list.CorMethod[*]=0
list.masknr=0
list.masknr1=0
list.backmod=-1
case list.SCANDIM of
1:    list.nrtif=list.N
2:    list.nrtif=list.CR[0]*list.CR[1]
endcase
ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
widget_control,ID,set_value=''
updatefilecounter,ev,list
printw,list.prid,'Reduced data dynamic list: Deleted'

; ----Load Mask----
Time=systime(1)
printw,list.prid,'Load Mask....'
LoadMask,list,0,list.prid,flagi=['$PDDPath:','$TILTP:','$LABELS:'],$
    error=error,/otD
if error then begin
    printw,list.prid,'Provide correct Mask File!'
    return,0B
endif

printw,list.prid,'Computation time:'+string(systime(1)-Time)+' sec'
Time=systime(1)

; ----Prepare spatial distortion correction----
if list.CorMethod[(*(list.CorSeq)).spatial] ne 0 then begin
    PTR_FREE,list.SpaDisStruc
    if (list.nrTieRing ne 0) or ((*list.nrTie)[0] eq 0) then begin

        ; Calc TieO point if necessary
        if (list.TieGrid eq 0) then begin
            ptr_free,list.TieO
            list.TieO=CalcTieO(list.nrTieRing,list.nrTie,list.TieI,list.center,$
                            list.lambda,list.dist,list.scrx,list.scry,list.TieRingD,list.a,list.b)
        endif

        ; Do the rest on first file because we need the size
    endif
    printw,list.prid,'Output Tiepoints for spatial distortion: '+string(systime(1)-Time)+' sec'
    Time=systime(1)
endif

; ----Load background image if necessary, For certain cases: SMALL2D=0----
bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
case bkcode[0] of
0:    begin
    printw,list.prid,'No 2D Background subtraction'
    list.SMALL2D=0
    endcase
1:    begin
    ; Don't strip ROI's if other corrections
    ; has to be performed afterwards
    i=(where(list.SeqTags eq 'back'))[0]
    j=n_elements(list.SeqTags)-1
    bool=0B
    while (j gt i) do begin
        bool or= (list.CorMethod[(*(list.CorSeq)).(j)]<1)
        j=j-1
    endwhile
    if bool then begin
        list.SMALL2D=0
        ID=widget_info(ev.top,FIND_BY_UNAME='Smallbutton2D')
        widget_control,ID,set_button=list.SMALL2D
    endif
    printw,list.prid,'2D Background stripping'
    endcase
2:    begin
    PTR_FREE,list.tiffb
    tiffb=ReadCCD(list.darkpath,list.darkfile,info=list.readinfo)
    if ptr_valid(tiffb) eq 0 then begin
       printw,list.prid,'Dark image not found:'+list.darkpath+list.darkfile
       return,0B
    endif
    if list.bcko eq 0 then list.tiffb=tiffb $
    else list.tiffb=mSavgol(data=tiffb,r=list.bcko,nr=list.bckw/2,nl=list.bckw/2,order=[0,0])
    (*list.tiffb)*=list.bckm
    list.SMALL2D=0
    ID=widget_info(ev.top,FIND_BY_UNAME='Smallbutton2D')
    widget_control,ID,set_button=list.SMALL2D
    printw,list.prid,'Dark image loaded: '+string(systime(1)-Time)+' sec'
    Time=systime(1)
    endcase
endcase
if bkcode[1] or bkcode[2] then printw,list.prid,'Thresholding in 2D Background subtraction'

; ----Load flat field image if nescessary----
if list.CorMethod[(*(list.CorSeq)).ff] ne 0 then begin
    ff=ReadCCD(list.ffpath,list.fffile,info=list.readinfo)
    if ptr_valid(ff) eq 0 then begin
        printw,list.prid,'Flat Field image not found: '+list.ffpath+list.fffile
        return,0B
    endif
    
    (*ff)=mean(*ff)/(*ff)
    PTR_FREE,list.ff
    list.ff=ff
    
    printw,list.prid,'Flat Field image loaded: '+string(systime(1)-Time)+' sec'
    Time=systime(1)
endif

; ----Load mask image if nescessary----
if list.CorMethod[(*(list.CorSeq)).maskoff] ne 0 and list.maskofffile ne '' then begin
    maskoff=ReadCCD(list.maskoffpath,list.maskofffile,info=list.readinfo)
    if ~ptr_valid(maskoff) and ~list.maskoffdyn then begin
        printw,list.prid,'Mask image not found: '+list.maskoffpath+list.maskofffile
        return,0B
    endif
    
    compareself,maskoff,list.maskoffthres
    ptr_free,list.maskoff
    list.maskoff=maskoff
        
    printw,list.prid,'Mask image loaded: '+string(systime(1)-Time)+' sec'
    Time=systime(1)
endif

; ----Init Sum model----
printw,list.prid,'Building the analyzing model:....'
case list.FASTMODE of
0:  sum=fltarr(list.nrtif)
1:    begin
    DestroyModel,list.ptrModel

    d0=list.AzIntParam[0]
    d1=list.AzIntParam[1]
    if d0 eq d1 and list.masknr ne 0 then begin
        ind=where((*list.mask)[0,*] ne 1,ct)
        if ct eq 0 then begin
            d0=0
            d1=0
        endif else begin
            ind=ind[0]
            if (*list.mask)[0,ind] eq 2 then begin
                d0=(*list.mask)[3,ind]<(*list.mask)[5,ind]
                d1=(*list.mask)[3,ind]>(*list.mask)[5,ind]
                a0=(*list.mask)[4,ind]<(*list.mask)[6,ind]
                a1=(*list.mask)[4,ind]>(*list.mask)[6,ind]
            endif else begin
                d0=(*list.mask)[3,ind]<(*list.mask)[4,ind]
                d1=(*list.mask)[3,ind]>(*list.mask)[4,ind]
                a0=0
                a1=2*!pi
            endelse
            if a0 eq a1 then begin
                a0=0
                a1=2*!pi
            endif

            tmp=BraggDtoX([d0,d1],list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/onlyt)
            nd=ceil(abs(tmp[1]-tmp[0])/list.AzIntInfo.BWT)+1

            list.AzIntParam=[d0,d1,nd,a0,a1]
        endelse
    endif
    if d0 eq d1 then begin
        printw,list.prid,'Mask file must contain azimuthal integration info!'
        return,0B
    endif

    list.masknr=0
    ptr_free,list.mask
    list.masknr1=1
    ptr_free,list.mask1
    tmp=(d0+d1)/2.
    list.mask1=ptr_new([1,0,0,tmp,tmp])

    sum=0
    endcase
2:    begin
    ; Model
    DestroyModel,list.ptrModel

    ; 2D mask
    if list.masknr eq 0 then begin
        printw,list.prid,'No 2D area found'
        return,0B
    endif
    ind=where( (*list.mask)[0,*] eq 1,ct)
    if ct eq 0 then begin
        printw,list.prid,'No 2D area found'
        return,0B
    endif
    list.masknr=1
    (*list.mask)=(*list.mask)[*,ind[0]]
    (*list.mask)[1]=0

    xs=(*list.mask)[4]-(*list.mask)[3]+1
    ys=(*list.mask)[6]-(*list.mask)[5]+1

    ; image size memory issue
    ; TODO: fix memory issue
    s=(xs>ys)<round(4000/sqrt(100))
    sum=fltarr(s,s,list.nrtif)

    ; 1D mask
    list.masknr1=0
    ptr_free,list.mask1
    endcase
endcase

ptrCur=(*list.ptrModel).ptrSum
ptrN=list.ptrModel
ptrNS=ptrCur

; ----Add model type0----
if (list.masknr ne 0) then begin
    imgmult=1b
    if list.masknr gt 1 then linktemp=reform((*list.mask)[1,*]) else linktemp=(*list.mask)[1,0:list.masknr-1]
    nrgroups=fix(max(linktemp))+1 ; sum type2
    linkind=lindgen(nrgroups)

    for i=0l,nrgroups-1 do begin
        ind=where(linktemp eq linkind[i],ct)
        if ct ne 0 then begin
            masknr=ct

            ; ----Add model entree----
            ; 0: 2D ROI summation
            ; => cryst. info: none
            ; => 1 group = a number of 2D ROIs
            ; => fit parameters: none
            ; => growth: no
            tmp=(*list.mask)[*,ind]
            tmp[1,*]=ind ; exchange groupnumber for masknumber

            (*ptrCur).next=PTR_NEW($
            {type:0,$                ; Model type
            group:linkind[i],$        ; Group number for type 0
            masknr:masknr,$            ; number of 2D ROI's
            mask:tmp,$                ; ROI specifications => second column changed from groupnumber to area number
            cor:fltarr(4,masknr),$    ; Corners of the ROI's
            imgmult:PTRARR(masknr),$; indices of the array defined by cor that are included in the summation
            code:replicate(list.SMALL2D,masknr),$    ; calc+sub background for each ROI?
            sum:ptr_new(sum),$        ; group sum for each file
            next:PTR_NEW()})        ; points to next model

            ptrCur=(*ptrCur).next

            ; ----Fill entree structure----
            for j=0l,masknr-1 do begin
                Cor=CornersMask((*list.mask)[*,ind[j]],list,imgmult=imgmult)
                (*ptrCur).Cor[*,j]=Cor
                (*ptrCur).imgmult[j]=PTR_NEW(imgmult)
            endfor

            (*ptrNS).n[0]++
            (*ptrN).n++
        endif else begin
           printw,list.prid,'LINK GROUPNUMBER NOT PRESENT!!!!'
           return,0B
        endelse
    endfor

    ; ----Add model type2----
    (*ptrCur).next=PTR_NEW($
        {type:2,$
        group:nrgroups,$
        code:list.SMALL2D,$
        sum:ptr_new(sum),$
        next:PTR_NEW()})
    (*ptrNS).n[2]++
    (*ptrN).n++
    ptrCur=(*ptrCur).next
endif

; ----Add model type1----
if (list.masknr1 ne 0) then begin
    ; AzIntParam: d1 d2 noutbins azbegin azend
    res=BraggDtoX(list.AzIntParam[0:1],list.lambda,/onlyt)
    b=res[0]<res[1]
    e=res[0]>res[1]
    Profilex=b+indgen(list.AzIntParam[2])*(e-b)/(list.AzIntParam[2]-1) ; in two-theta (radians)

    if list.FASTMODE eq 1 then sum=Profilex*180/!pi
    if list.masknr1 gt 1 then linktemp=reform((*list.mask1)[1,0:list.masknr1-1]) $
        else linktemp=(*list.mask1)[1,0:list.masknr1-1]
    nrgroups=fix(max(linktemp))+1
    linkind=lindgen(nrgroups)

    for i=0l,nrgroups-1 do begin
        ind=where(linktemp eq linkind[i],ct)
        if ct ne 0 then begin
            masknr=ct

            ; ----Add model entree----
            ; 1: 1D ROI summation
            ; => cryst. info: none
            ; => 1 group = a number of 1D ROIs
            ; => fit parameters: none
            ; => growth: no
            tmp=(*list.mask1)[*,ind]
            tmp[1,*]=ind ; exchange groupnumber for masknumber

            (*ptrCur).next=PTR_NEW($
            {type:1,$                ; Model type
            group:linkind[i],$        ; Group number for type1
            masknr:masknr,$            ; number of 1D ROI's
            mask:tmp,$                ; ROI specifications => first column changed from groupnumber to area number
            cor:fltarr(2,masknr),$    ; Corners of the ROI's
            code:replicate(list.SMALL1D,masknr),$    ; calc+sub background for each ROI?
            sum:ptr_new(sum),$        ; group sum for each file
            sum2:ptr_new(list.AzIntInfo.calcerror?sum:0),$        
            next:PTR_NEW()})        ; points to next model

            ; ----Fill entree structure----
            tmp[3:4,*]=BraggDtoX(tmp[3:4,*],list.lambda,/onlyt)
            for j=0l,masknr-1 do begin
                t=min(abs(Profilex-tmp[3,j]),x1);convert from two-theta to index in Profilex
                t=min(abs(Profilex-tmp[4,j]),x2)
                (*(*ptrCur).next).cor[*,j]=[x1<x2,x1>x2]
            endfor

            ptrCur=(*ptrCur).next
            (*ptrNS).n[1]++
            (*ptrN).n++
        endif else begin
           printw,list.prid,'LINK GROUPNUMBER NOT PRESENT!!!!'
           return,0B
        endelse
    endfor

    ; ----Add model type2----
    ; 2: 2D/1D ROI summation of whole pattern/profile minus the other ROI's (-> rest)
    ; => cryst. info: none
    ; => 1 group = the rest image
    ; => fit parameters: none
    ; => growth: no
    if list.FASTMODE eq 0 then begin
        (*ptrCur).next=PTR_NEW($
            {type:2,$                ; Model type
            group:nrgroups,$        ; Group number for type2
            code:list.SMALL1D,$        ; calc+sub background for each ROI?
            sum:ptr_new(sum),$        ; group sum for each file
            next:PTR_NEW()})        ; points to next model
        (*ptrNS).n[2]++
        (*ptrN).n++
    endif

    ptrCur=(*ptrCur).next
endif

; ----Edit model type10, type 11 and type12----
; allocation already done

; ----Proceed?----
if ((*ptrN).n eq 0) then begin
    printw,list.prid,'Mask file must contain 2D or 1D processing info!'
    return,0B
endif

printw,list.prid,'Computation time:'+string(systime(1)-Time)+' sec'
Time=systime(1)

; ----Make output window----
resetmag_xrdbp,list
result=PrepOut(list,ev)
Time=systime(1)

; ----Number of files to process----
printw,list.prid,'Directory: '+list.path
if n_elements(*list.filter) eq 1 then $
printw,list.prid,'Files: '+(*list.filter)
printw,list.prid,'Number of patterns left:'+string(list.nrtif-list.nFiles)

return,1B
end;function BatchInit
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Process2DGroups,Pattern,Patternvalid,list

ptrCur=(*(*list.ptrModel).ptrSum).next

ptrFirst=PTR_NEW()
while PTR_VALID(ptrCur) do begin
    if (*ptrCur).type ne 0 then ptrCur=(*ptrCur).next else begin
        ptrFirst=ptrCur
        ptrCur=PTR_NEW()
    endelse
endwhile

if PTR_VALID(ptrFirst) then ptrCur=ptrFirst else return,1b

if list.FASTMODE eq 2 then begin
    i=0
    x0=0>(*ptrCur).cor[0,i]<((*list.tiffs)[1]-1)
    x1=0>(*ptrCur).cor[1,i]<((*list.tiffs)[1]-1)
    y0=0>(*ptrCur).cor[2,i]<((*list.tiffs)[2]-1)
    y1=0>(*ptrCur).cor[3,i]<((*list.tiffs)[2]-1)
    img=(*Pattern)[x0:x1,y0:y1]
    if ptr_valid(Patternvalid) then img*=(*Patternvalid)[x0:x1,y0:y1]
    s=size(*(*ptrCur).sum)
    img=congrid(img,s[1],s[2])
    (*(*ptrCur).sum)[*,*,list.nFiles]=img
    return,1b
endif

TotPat=*Pattern
if ptr_valid(Patternvalid) then begin
    if list.NORMROIs then nvalidtotal=total(*Patternvalid)
    TotPat*=(*Patternvalid)
endif else if list.NORMROIs then nvalidtotal=n_elements(*Pattern)

repeat begin

    ; Loop over all ROIs in this group
    for i=0l,(*ptrCur).masknr-1 do begin

        ; Select area of Pattern
        c0=(*ptrCur).cor[0,i]
        c1=(*ptrCur).cor[2,i]
        s0=(*ptrCur).cor[1,i]-c0+1
        s1=(*ptrCur).cor[3,i]-c1+1
        img=extrac(*Pattern,c0,c1,s0,s1)
        if ptr_valid(Patternvalid) then begin
            imgvalid=extrac(*Patternvalid,c0,c1,s0,s1)
            img*=imgvalid
        endif
        
        ; Calc background from this area
       if (*ptrCur).code[i] then begin
         bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
         BackSub,img,Strip2D(img,list.wb,list.niter),bkcode[0]
         if bkcode[1] then begin
             ind=where(img le list.thres[0],ct)
             if ct ne 0 then img[ind]=0
         endif
         if bkcode[2] then begin
             ind=where(img ge list.thres[1],ct)
             if ct ne 0 then img[ind]=0
         endif
       endif

       ; Add to groupsum
       A=total(img[*((*ptrCur).imgmult[i])])
       if list.debug then begin
              if list.outid2 ne -1 then wset,list.outid2 $
              else begin
                    window,winid(),title='ROI: masked off image'
                    list.outid2=!d.window
              endelse
              imgtmp=make_array(size=size(img),value=0b)
              imgtmp[*((*ptrCur).imgmult[i])]=1
              imgtmp*=img
              if ptr_valid(Patternvalid) then imgtmp*=imgvalid
              tvscl,congrid(imgtmp<mean(imgtmp[*((*ptrCur).imgmult[i])]),!d.x_size<!d.y_size,!d.x_size<!d.y_size)
       endif
       if list.NORMROIs then begin
           if ptr_valid(Patternvalid) then A/=total(imgvalid[*((*ptrCur).imgmult[i])])$
           else A/=n_elements(*((*ptrCur).imgmult[i]))
       endif
       (*(*ptrCur).sum)[list.nFiles]+=A
       
       ; Rest
       ind=*((*ptrCur).imgmult[i])
       s=size(img,/dim)
       indx=ind mod s[1] + c0
       indy=ind/s[1] + c1
       TotPat[indx,indy]=0
    endfor

    ; Stop when the next is NULL or type1
    ptrCur=(*ptrCur).next
    bool=PTR_VALID(ptrCur)
    if bool then bool=(*ptrCur).type eq 0
endrep until (not bool)

if PTR_VALID(ptrCur) then begin ; Points to type2
    (*(*ptrCur).sum)[list.nFiles]=(total(TotPat)>0)
    if list.NORMROIs then $
        (*(*ptrCur).sum)[list.nFiles]/=nvalidtotal
endif

return,1b
end;function Process2DGroups
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro    Process1DROI,list,i,Profilex,Profiley,ptrCur,x,y,yb,cb,ce

; Recalc ROI in bins
tmp=BraggDtoX((*ptrCur).mask[3:4,i],list.lambda,/onlyt,/angledeg)
t=min(abs(Profilex-tmp[0]),x1);convert from two-theta to index in Profilex
t=min(abs(Profilex-tmp[1]),x2)
(*ptrCur).cor[*,i]=[x1<x2,x1>x2]

; Select area of Pattern
n=n_elements(Profiley)-1
cb=0>(*ptrCur).cor[0,i]<n
ce=0>(*ptrCur).cor[1,i]<n
y=Profiley[cb:ce]
x=Profilex[cb:ce]

; Calc background for this area
if (*ptrCur).code[i] then begin
    case list.backmod of
    -1:
    2:    begin
         yb=Strip1D(y,(*(*list.ptrModel).ptrFit).backparamIO[2],(*(*list.ptrModel).ptrFit).backparamIO[1])
         BackSub,y,yb,list.backmod,offset=offset
         if n_elements(offset) ne 0 then yb=(yb+offset)>0
         endcase
    0:    begin
         yb=ortpol(x,y,degree=(*(*list.ptrModel).ptrFit).backparamIO[0],GR=GR)
         BackSub,y,yb,list.backmod,offset=offset
         if n_elements(offset) ne 0 then yb=(yb+offset)>0
         endcase
    1:    begin
         yb=ortpol(x,y,degree=(*(*list.ptrModel).ptrFit).backparamIO[0],GR=GR)
         BackSub,y,yb,list.backmod,offset=offset
         if n_elements(offset) ne 0 then yb=(yb+offset)>0
         endcase
    endcase
endif
end;pro    Process1DROI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Process1DGroupsS,Profile,list,dummy=dummy

Profilex=reform(Profile[0,*])
Profiley=reform(Profile[1,*])
if list.AzIntInfo.calcerror then Profileyerror=reform(Profile[2,*])
ptrCur=(*(*list.ptrModel).ptrSum).next

ptrFirst=PTR_NEW()
while PTR_VALID(ptrCur) do begin
    if (*ptrCur).type ne 1 then ptrCur=(*ptrCur).next else begin
        ptrFirst=ptrCur
        ptrCur=PTR_NEW()
    endelse
endwhile

if PTR_VALID(ptrFirst) then ptrCur=ptrFirst else return,1b

if list.FASTMODE eq 1 then begin
    s=size(*(*ptrCur).sum,/dimensions)
    nxnew=n_elements(Profiley)
    if s[0] ne nxnew then begin
        ; When using .chi files: ignore x from mask file
        if n_elements(s) eq 1 then begin ; first file: use X
            list.AzIntParam[2]=n_elements(Profilex)
            *(*ptrCur).sum=Profilex
            s=size(*(*ptrCur).sum,/dimensions)
            if list.AzIntInfo.calcerror then *(*ptrCur).sum2=Profilex
        endif else begin ; not first file: treat as dummy
            Profiley=fltarr(s[0])
            printw,list.prid,'Use dummy: not the same size...'
        endelse
    endif
    if keyword_set(dummy) then begin
        Profiley=fltarr(s[0])
        if list.AzIntInfo.calcerror then Profileyerror=Profiley
    endif
    *(*ptrCur).sum=[[*(*ptrCur).sum],[Profiley]]
    if list.AzIntInfo.calcerror then *(*ptrCur).sum2=[[*(*ptrCur).sum2],[Profileyerror]]
    return,1b
endif

if list.debug then begin
    if list.outid2 ne -1 then wset,list.outid2 $
    else begin
        window,winid(),title='ROI: 1D'
        list.outid2=!d.window
    endelse
    plot,Profilex,Profiley,xstyle=1,ystyle=1
    ROIy=!Y.CRange
endif

TotPro=Profiley
repeat begin

    ; Loop over all ROIs in this group
    for i=0l,(*ptrCur).masknr-1 do begin

       Process1DROI,list,i,Profilex,Profiley,ptrCur,x,y,yb,cb,ce

       ; Add to groupsum
       A=total(y)
       if list.debug then begin
              n=n_elements(x)
                 if n_elements(yb) ne 0 then oplot,x,yb,thick=2
              PlotS, [x[0],x[0]], ROIy, /Data,linestyle=2,NOCLIP = 0
              PlotS, [x[n-1],x[n-1]], ROIy, /Data,linestyle=2,NOCLIP = 0
       endif
       if list.NORMROIs then A/=n_elements(y)
       (*(*ptrCur).sum)[list.nFiles]+=A
       
       ; Rest
       TotPro[cb:ce]=0
    endfor

    ; Stop when the next is NULL or type2
    ptrCur=(*ptrCur).next
    bool=PTR_VALID(ptrCur)
    if bool then bool=(*ptrCur).type eq 1
endrep until (not bool)
TotPro=total(TotPro)

if PTR_VALID(ptrCur) then begin ; Points to type2
    (*(*ptrCur).sum)[list.nFiles]=(TotPro>0)
    if list.NORMROIs then $
        (*(*ptrCur).sum)[list.nFiles]/=n_elements(Profiley)
endif

return,1b

end;function Process1DGroupsS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ProcessFile,list,ev

; ----Open tif file----
T=systime(1)

bdoTD=1b
bFreePattern=1b
if list.Subfiles.nleft ne 0 then begin
    ; File already in memory
    
    ; Read fixed block
    if list.Subfiles.fileoffset ne 0 then begin
        ind=indgen(list.Subfiles.fileoffset)
        case list.Subfiles.fileind of
        0: xval=reform((*list.Subfiles.ptr)[ind,*,*])
        1: xval=reform((*list.Subfiles.ptr)[*,ind,*])
        2: xval=(*list.Subfiles.ptr)[*,*,ind]
        else: return,0
        endcase
    endif

    ; Read current spectrum from block1
    ind=list.Subfiles.n-list.Subfiles.nleft+list.Subfiles.fileoffset
    case list.Subfiles.fileind of
        0:     spe=reform((*list.Subfiles.ptr)[ind,*,*])
        1:     spe=reform((*list.Subfiles.ptr)[*,ind,*])
        2:     spe=(*list.Subfiles.ptr)[*,*,ind]
        else: return,0
    endcase
    
    ; Read current spectrum from block2
    if list.Subfiles.filemult ge 2 then begin
        ind=2*list.Subfiles.n-list.Subfiles.nleft+list.Subfiles.fileoffset
        case list.Subfiles.fileind of
            0:     speerror=reform((*list.Subfiles.ptr)[ind,*,*])
            1:     speerror=reform((*list.Subfiles.ptr)[*,ind,*])
            2:     speerror=(*list.Subfiles.ptr)[*,*,ind]
            else: return,0
        endcase
    endif

    case list.Subfiles.ptrtype of
    1:    begin
        bdoTD=0b
        bFreePattern=0b
        xtype=list.Subfiles.xtype
        if n_elements(speerror) eq 0 then speerror=sqrt(spe)
        endcase
    2:    Pattern=ptr_new(temporary(spe))
    endcase

    list.Subfiles.nleft--
    if list.Subfiles.nleft eq 0 then heap_free,list.Subfiles
    error=0b
endif else begin

    filename=(*list.PTR_files)[list.nfiles]
    printw,list.prid,'Reading: '+filename+'...'
    f=FileFormat(0,filename)

    case f[0] of
    2:    begin
        Pattern=ReadCCD(filename,info=list.readinfo)
        error=ptr_valid(Pattern) eq 0
        endcase
    else:begin
        bdoTD=0b
        bFreePattern=0b
        error=CutPath(filename,path=path,file=file,ext=ext)
        error=~ReadScan(path,file+ext,ext,list.prid,xval=xval,spe=spe,stdspe=speerror,type=xtype)

        if f[0] eq 3 and ~error then begin
            heap_free,list.Subfiles

            s=size(spe,/dimensions)
            if n_elements(s) eq 1 then list.Subfiles.n=1 $
            else list.Subfiles.n=s[1]

            list.Subfiles.nleft=list.Subfiles.n-1
            if list.Subfiles.nleft ne 0 then begin
                list.Subfiles.fileind=1
                list.Subfiles.ptrtype=1
                list.Subfiles.fileoffset=1
                list.Subfiles.filemult=2
                list.Subfiles.xtype=xtype
                list.Subfiles.ptr=ptr_new([[xval],$ ; fixed block
                                            [spe],$ ; block1 (spectra)
                                            [speerror]]) ; block2 (spectra stdev's)

                ; Expand file list
                if list.nfiles lt (n_elements(*list.PTR_files)-1) then $
                    *list.PTR_files=[(*list.PTR_files)[0:list.nfiles],$
                                replicate(filename,list.Subfiles.nleft),$
                                (*list.PTR_files)[list.nfiles+1:*]] $
                else $
                    *list.PTR_files=[(*list.PTR_files)[0:list.nfiles],$
                                replicate(filename,list.Subfiles.nleft)]
                                
                                
                ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
                if ptr_valid(list.PTR_files) then widget_control,ID,set_value=*list.PTR_files
                widget_control,ID,SET_LIST_TOP=n_elements(*list.PTR_files)-1
                updatefilecounter,ev,list
        
            endif

            spe=spe[*,0]
            speerror=speerror[*,0]
        endif

        endelse
    endcase
endelse

nTypesS=(*(*list.ptrModel).ptrSum).n
tnTypesF=total((*(*list.ptrModel).ptrFit).n)

if error then begin

    if (list.filtermethod eq 2) or (list.nocheck eq 1) then begin
        ; list-of-files contains a dummy file: all areas are zero
        printw,list.prid,'Dummy file: take 0 for all areas'

        if (nTypesS[1] ne 0) and list.FASTMODE eq 1 then $
            temp=Process1DGroupsS(intarr(2,2),list,/dummy)

        if tnTypesF ne 0 then FitModelBatch,[0,1],[0,0],[0,0],list,/dummy
    endif else begin
        ; delete filename from list
        *list.PTR_files=shrinkarray(*list.PTR_files,list.nFiles,count=nfiles)
        if nfiles eq 0 then begin
            ptr_free,list.PTR_files
            list.PTR_files=PTR_NEW()
            widget_control,ev.top,set_uvalue=list
        endif

        ; ----Update list----
        ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
        if ptr_valid(list.PTR_files) then widget_control,ID,set_value=*list.PTR_files
        widget_control,ID,SET_LIST_TOP=nfiles-1

        ; ----Update file-counter----
        updatefilecounter,ev,list
        
        return,0
    endelse
endif else begin

if (list.update eq 1) and widget_info(list.top,/valid_id) and bdoTD then begin
    MainUpdateCManual,filename,list.top,list.ev,fileptr=Pattern
endif

printw,list.prid,'Reading time: '+string(systime(1)-T)+' sec'
T=systime(1)

; ----Apply corrections (+ calc 2D background)----
if bdoTD then begin
ptr_free,list.tiffs
list.tiffs = PTR_NEW(size(*Pattern))
Patternvalid=ptr_new()

for i=0l,n_elements(list.CorMethod)-1 do begin
    index=(*(list.CorSeq)).(i)
    if (list.CorMethod[index] ne 0) then $
    case list.SeqTags[i] of
    'zinger':    begin;zinger removal
        RemoveZingers,Pattern,Patternvalid,width=list.ZingerWidth,threshold=list.ZingerThreshold,count=count
        printw,list.prid,'Zingers found: '+stringr(count)+' in '+string(systime(1)-T)+' sec'
        T=systime(1)
        endcase
    'satur':    begin;saturation detection
        RemoveSatur,Pattern,Patternvalid,list.SaturMax,list.SaturAspect,list.SaturNMax,count=count
        printw,list.prid,'Saturated pixels found: '+stringr(count)+' in '+string(systime(1)-T)+' sec'
        T=systime(1)
        endcase
    'back':    begin;background correction
        bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
        bthres=1b
        case bkcode[0] of
        0:
        1:    begin
            ; Calc background of whole image
            ; When doing this, make sure that background calculations
            ; aren't done on the ROI's
            if not list.SMALL2D then begin
                p=Strip2D(Pattern,list.wb,list.niter)
                BackSub,Pattern,p,bkcode[0]
                ptr_free,p
                printw,list.prid,'Background stripping: '+string(systime(1)-T)+' sec'
                T=systime(1)
            endif else begin
                bthres=0b
                printw,list.prid,'Background stripping for each ROI.'
            endelse
            endcase
        2:    begin
            BackSub,Pattern,list.tiffb,bkcode[0]
            printw,list.prid,'Background subtraction: '+string(systime(1)-T)+' sec'
            T=systime(1)
            endcase
        endcase
        if not (bkcode[0] eq 1 and list.SMALL2D) then begin
        if bkcode[1] and bthres then begin
              ind=where((*Pattern) le list.thres[0],ct)
              if ct ne 0 then (*Pattern)[ind]=0
              printw,list.prid,'Lower thresholding.'
        endif
        if bkcode[2] and bthres then begin
              ind=where((*Pattern) ge list.thres[1],ct)
              if ct ne 0 then (*Pattern)[ind]=0
              printw,list.prid,'Higher thresholding.'
        endif
        endif

        endcase
    'spatial':    begin;spatial distortion
    
        if not PTR_VALID(list.SpaDisStruc) then begin
            if list.Tiegrid eq 2 then begin
                if not PTR_VALID(list.SpaDisStruc) then begin
                    struc=ReadFit2DSpline(list.TiePath+list.TieFile,error=error)
                    if error then begin
                        printw,list.prid,'Error on loading '+list.TiePath+list.TieFile
                    endif else begin
                        list.scrx=struc.scrx
                        list.scry=struc.scry
                        list.SpaDisStruc=PTR_NEW(struc)
                    endelse
                endif
                if PTR_VALID(list.SpaDisStruc) then $
                    list.SpaDisStruc=WarpSpaDistPrepSpline(list.tiffs,(*list.SpaDisStruc).tck1,(*list.SpaDisStruc).tck2,vflip=list.vFlipSpline)
            endif else begin
                if (list.nrTieRing eq 0) and ((*list.nrTie)[0] eq 0) then begin
                    printw,list.prid,'No tie points defined'
                endif else begin
                    list.SpaDisStruc=WarpSpaDistPrep(list.tiffs,list.TieI,list.TieO,fast=list.bFastSpline)
                endelse
            endelse
        endif
                
        if not PTR_VALID(list.SpaDisStruc) then begin
            printw,list.prid,'No Spa Dist information!!!!!!!!'
        endif else begin
            if list.Tiegrid eq 2 then WarpSpaDistRunSpline,Pattern,list.SpaDisStruc,tiffvalid=Patternvalid,fast=list.bFastSpline $
            else WarpSpaDistRun,Pattern,list.SpaDisStruc,tiffvalid=Patternvalid

            printw,list.prid,'Spa Dist Correction: '+string(systime(1)-T)+' sec'
        endelse
        T=systime(1)

        endcase
    'ff':    begin;ff
        type=size(*Pattern,/type)
        (*Pattern)*=(*list.ff)
        ;p=TypeConvert(Pattern,type)
        printw,list.prid,'Flat Field correction: '+string(systime(1)-T)+' sec'
        T=systime(1)
        endcase
    'maskoff':    begin;maskff
        
        if list.maskoffdyn then begin
            if not ptr_valid(Patternvalid) then Patternvalid = ptr_new(compare(Pattern,list.maskoffdynthres)) $
            else (*Patternvalid) and= compare(Pattern,list.maskoffdynthres)
        endif
        
        if ptr_valid(list.maskoff) then begin
            if not ptr_valid(Patternvalid) then Patternvalid = ptr_new(*list.maskoff) $
            else (*Patternvalid) and= (*list.maskoff)
        endif
        
        printw,list.prid,'Masking off: '+string(systime(1)-T)+' sec'
        T=systime(1)
        endcase
    endcase
endfor
endif

; Debug: check corrections
;pathtmp='C:\projects\idl\xrdua\xrdua\woutdenolf\test\correction\'
;PatternComp=ReadCCDtiff(pathtmp+'compareCCD.tif')
;PatternvalidComp=ReadCCDtiff(pathtmp+'compareValid.tif')
;help,where((*Patternvalid) ne PatternvalidComp)
;help,where((*Pattern) ne PatternComp)
printw,list.prid,'Reduction:....'

; ----Process All 2D groups----
if (nTypesS[0] ne 0) and bdoTD then temp=Process2DGroups(Pattern,Patternvalid,list)

; ----Azimuthal integration ----
if (nTypesS[1]+tnTypesF ne 0) then begin
    if bdoTD then begin
        if ptr_valid(list.WarpStruc) then begin
            Profile=AzIntWarpBatch(list.WarpStruc,Pattern,Patternvalid)
        endif else begin
            batchreturnp=ptr_new()
            res=BraggPtoX(list.AzIntInfo.MinPix,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
            rmin=[res[1],res[0],list.AzIntInfo.MinPix,res[2]]
            dr=[list.AzIntInfo.BWD,list.AzIntInfo.BWT,list.AzIntInfo.BWP,list.AzIntInfo.BWQ]
            Profile=AzIntWarp([2,0,0,list.AzIntParam[0],list.AzIntParam[3],list.AzIntParam[1],list.AzIntParam[4]],$
                {tiff:Pattern,tiffs:list.tiffs,dist:list.dist,scrx:list.scrx,scry:list.scry,center:list.center,a:list.a,b:list.b,phipixphi0:list.phipixphi0,lambda:list.lambda},$
                list.IntCor2Dinfo,1,rmin,dr,list.AzIntParam[2],list.AzIntInfo.AzWidth,0,[0B,list.AzIntInfo.bool],imagevalid=Patternvalid,batchreturnp=batchreturnp,median=list.AzIntInfo.median,$
                percentile=list.AzIntInfo.percentile,updatenorm=list.AzIntInfo.updatenorm,calcerror=list.AzIntInfo.calcerror)
            heap_free,list.WarpStruc
            list.WarpStruc=batchreturnp
        endelse
        
        if ~list.AzIntInfo.calcerror then Profile=[Profile,sqrt(Profile[1,*])]
        
        ; Debug: check azimuthal integration
;        tmp=read_chi('C:\projects\idl\xrdua\xrdua\woutdenolf\test\correction\grid1.chi')
;        help,where(Profile ne tmp)
        
    endif else begin
        if xtype ne 1 then begin
            xvalt=CHI_xval(xval,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,xtype,spe=spe,tspe=spet,stdspe=speerror,tstdspe=speerrort)
            xval=xvalt[*,1]
            spe=spet[*,1]
            speerror=speerrort[*,1]
        endif

        temp=value_locater(xval,list.backranindtt)>0
        temp=temp[sort(temp)]
        if temp[0] ne temp[1] then begin
            xval=xval[temp[0]:temp[1]] ;two-theta
            spe=spe[temp[0]:temp[1]]
            speerror=speerror[temp[0]:temp[1]]
        endif
        Profile=transpose([[xval],[spe],[speerror]])
    endelse

    if (list.update eq 1) and widget_info(list.top,/valid_id) then begin
        temp=min(Profile[1,*],y1)
        temp=max(Profile[1,*],y2)
        xtype=1

        fileptr={xval:reform(Profile[0,*]),spe:reform(Profile[1,*]),speerror:reform(Profile[2,*]),ROIAzimuth:list.AzIntParam[3:4],ytitle:'Intensity (a.u.)',$
            xtitle:ChiXtitle(xtype),title:(*list.PTR_files)[list.nfiles],xtype:xtype,$
            xrange:[0,n_elements(xval)-1],yrange:[y1[0],y2[0]]}

        MainUpdateCManual,'c:\dummy.chi',list.top,list.ev,fileptr=fileptr
    endif
endif

; ----Process All 1D groups of type 1 and 2----
if (nTypesS[1] ne 0) then begin
    ProfileS=Profile
    if not list.SMALL1D then begin
        y=reform(ProfileS[1,*])
        x=reform(ProfileS[0,*])
        case list.backmod of
               -1:
               2:    begin
                     yb=Strip1D(y,(*(*list.ptrModel).ptrFit).backparamIO[2],(*(*list.ptrModel).ptrFit).backparamIO[1])
                     BackSub,y,yb,list.backmod,offset=offset
                     if n_elements(offset) ne 0 then yb=(yb+offset)>0
                     endcase
               0:    begin
                     yb=ortpol(x,y,degree=(*(*list.ptrModel).ptrFit).backparamIO[0],GR=GR)
                     BackSub,y,yb,list.backmod,offset=offset
                     if n_elements(offset) ne 0 then yb=(yb+offset)>0
                     endcase
               1:    begin
                     yb=ortpol(x,y,degree=(*(*list.ptrModel).ptrFit).backparamIO[0],GR=GR)
                     BackSub,y,yb,list.backmod,offset=offset
                     if n_elements(offset) ne 0 then yb=(yb+offset)>0
                     endcase
        endcase
        ProfileS[1,*]=y
    endif

    temp=Process1DGroupsS(ProfileS,list)
endif

; ----Process All 1D groups of other types----
if tnTypesF ne 0 then begin
    FitModelBatch,reform(Profile[0,*]),reform(Profile[1,*]),reform(Profile[2,*]),list,error=error
    if error then begin
        printw,list.prid,'Batch fitting error!'
        ;return,0b
    endif
endif

printw,list.prid,'Computation time:'+string(systime(1)-T)+' sec'
T=systime(1)
endelse

; ----Save and plot reduced data----
list.nFiles++
RefreshDisplayBatch,list
printw,list.prid,'Save and plot results:'+string(systime(1)-T)+' sec'

if bFreePattern then ptr_free,Pattern
return,1B

end;function ProcessFile
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro bbbplot,list

LoadctRev,0,/silent

ptrCur=(*(*list.ptrModel).ptrSum).next
psum=(*ptrCur).sum

s=dimsize(*psum,3)
nrgroups=s[2]

; image positions
nx=list.xsx
ny=list.ysy
ind=indgen(nx,ny)
if list.REV3 and (list.ysy gt 1) then begin
    indodd=2*lindgen(list.ysy/2)+1
    ind[*,indodd]=reverse(ind[*,indodd],1)
endif
if list.TRAN then begin
    ind=transpose(ind)
    ny=list.xsx
    nx=list.ysy
endif
if list.REV1 and nx gt 1 then ind=reverse(ind,1)
if list.REV2 and ny gt 1 then ind=reverse(ind,2)

; plot
sx=s[0]*list.nmag
sy=s[1]*list.nmag
n=nrgroups<(nx*ny)

mi=min(*psum)
ma=median((*psum))+2*stdev((*psum))
ma-=(ma-mi)*list.nmult

for i=0l,n-1 do begin
    j=ind[i]
    ;if j lt n then tv,congrid(bytscl((*psum)[*,*,j],min=mi,max=ma),sx,sy),i
    if j lt n then begin
;        mi=min((*psum)[*,*,j])
;        ma=median((*psum)[*,*,j])+2*stdev((*psum)[*,*,j])
        tv,congrid(bytscl((*psum)[*,*,j],min=mi,max=ma),sx,sy),i
    endif
endfor

end;pro bbbplot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro bbplot,list,fastinfo

LoadctRev,-1,/silent

plot,fastinfo.x,sqrt(fastinfo.y),xtitle=ChiXtitle(1),ytitle='sqrt(I)',xrange=list.xrange,yrange=sqrt(list.yrange),/xs,/ys

tt=braggDtoX((*fastinfo.ptrCur).mask[3:4],/onlyt,list.lambda,/angledeg)
ROIy=!Y.CRange
PlotS, [tt[0],tt[0]], ROIy, /Data,linestyle=2,NOCLIP = 0
PlotS, [tt[1],tt[1]], ROIy, /Data,linestyle=2,NOCLIP = 0

ind1=where(fastinfo.x ge (tt[0]<tt[1]),ct)
if ct eq 0 then ind1=0 else ind1=ind1[0]
ind2=where(fastinfo.x le (tt[0]>tt[1]),ct)
if ct eq 0 then ind2=n_elements(fastinfo.x)-1 else ind2=ind2[ct-1]

n=ind2-ind1+1
if n gt 1 then begin
    pxval=fastinfo.x[ind1:ind2]
    pyval=sqrt(fastinfo.y[ind1:ind2])
    b=0
    if ~list.FORCE then begin
        pxval=[pxval[0],pxval,pxval[n-1]]
        tmp=min(!y.crange)
        pyval=[tmp,pyval,tmp]
    endif
    polyfill,pxval,pyval,/line_fill,orientation=45,col=16*10,/data
endif

end;pro bbplot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bplotgetxrange,list,nfiles,unit=unit
xrange=list.coord1
unit='file index'
if xrange[0] eq xrange[1] then xrange=[0,nfiles-1.] else $
case list.unit1 of
0: begin
    unit=msymbols('micro')+'m'
    xrange*=1000
    endcase
1: begin
    unit='deg'
    endcase
endcase
if list.TRAN then sign=list.Vsign else sign=list.Hsign
if sign gt 0 then xrange=xrange[sort(xrange)] else xrange=xrange[reverse(sort(xrange))]
return,xrange
end;function bplotgetxrange
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro bplot,list

statgroupi={may:0.,miy:0.}
MakeDatablock,list,block,nrgroups,dim1,dim2,/transform,statgroupi=statgroupi
groupi=0>list.groupi<(nrgroups-1)

case list.SCANDIM of
1:    begin
    LoadctRev,-1,/silent

    ; ----Display----
    nfiles=dim1*dim2
    xrange=bplotgetxrange(list,nfiles,unit=unit)
    x=MakeRange(xrange[0],xrange[1],(xrange[1]-xrange[0])/(nfiles-1.))
    if n_elements(x) eq 1 then return

    yrange=[statgroupi.miy,statgroupi.may]
    y=reform(block[groupi,*,*])
    
    if list.TRAN then begin
        plot,y,x,yrange=reverse(xrange),xrange=yrange,$
            xstyle=1,ystyle=1,ytitle='Spatial distance ('+unit+')',$
            xtitle='ROI value',title=list.plottitle
    endif else begin
        plot,x,y,xrange=xrange,yrange=yrange,$
             xstyle=1,ystyle=1,xtitle='Spatial distance ('+unit+')',$
            ytitle='ROI value',title=list.plottitle
    endelse

    endcase
2:    begin
    LoadctRev,-3,/silent
    yrange=[statgroupi.miy,statgroupi.may]
    range = yrange[1]-yrange[0]
    tv,rebin((!d.TABLE_SIZE-1)*((((reform(block[groupi,*,*])-yrange[0])>0)/(range*list.nmult))<1),$
                dim1*list.nmag,dim2*list.nmag,/sample),/order
    endcase
endcase

end;pro bplot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshDisplayBatch,list

if not list.PLOT or list.nFiles eq 0 then return

T=systime(1)

FastModelBegin,list,fastinfo

if n_elements(fastinfo) ne 0 then begin
    Wset,list.outdrawidsum
    bbplot,list,fastinfo
    GetSysVar,list.sysvar

    wset, list.pixindexsum
    h=!D.x_size
    v=!D.y_size
    Device, Copy = [0,0, h,v, 0,0,list.outdrawidsum]
endif

Wset,list.outdrawid
if list.FASTMODE eq 2 then bbbplot,list $
else bplot,list
GetSysVar,list.sysvar2

FastModelEnd,list,fastinfo

printw,list.prid,'Plotting: '+string(systime(1)-T)+' sec'

end;pro RefreshDisplayBatch
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MainUpdateDo,list,ev

if ~widget_info(list.top,/valid) then return

; Update 1D
nTypesS=(*(*list.ptrModel).ptrSum).n
if list.FASTMODE eq 1 and nTypesS[1] ne 0 then begin
    ptrCur=(*(*list.ptrModel).ptrSum).next

    ptrFirst=PTR_NEW()
    while PTR_VALID(ptrCur) do begin
        if (*ptrCur).type ne 1 then ptrCur=(*ptrCur).next else begin
            ptrFirst=ptrCur
            ptrCur=PTR_NEW()
        endelse
    endwhile
    if PTR_VALID(ptrFirst) then ptrCur=ptrFirst else return

    tmp=CutPath((*list.PTR_files)[ev.index],file=title)
    s=DimSize(*(*ptrCur).sum,2)
    iextract=ev.index+1
    if iextract lt s[1] then begin
        temp=min((*(*ptrCur).sum)[*,iextract],y1)
        temp=max((*(*ptrCur).sum)[*,iextract],y2)
        xtype=1
        
        if list.AzIntInfo.calcerror and (where(tag_names(*ptrCur) eq 'SUM2'))[0] ne -1 then speerror=(*(*ptrCur).sum2)[*,iextract] else speerror=sqrt((*(*ptrCur).sum)[*,iextract])
        
        fileptr={xval:(*(*ptrCur).sum)[*,0],spe:(*(*ptrCur).sum)[*,iextract],speerror:speerror,ROIAzimuth:list.AzIntParam[3:4],ytitle:'Intensity (a.u.)',$
            xtitle:ChiXtitle(xtype),title:title,xtype:xtype,$
            xrange:[0,n_elements(xval)-1],yrange:[y1[0],y2[0]]}
    
        MainUpdateCManual,'c:\dummy.chi',list.top,list.ev,fileptr=fileptr
    endif
endif

; Update 2D
MainUpdate,list,ev,/tdrequired

end;pro MainUpdateDo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro exportfastprofile,list

if list.update eq 0 then return

FastModelBegin,list,fastinfo
FastModelEnd,list,fastinfo
if n_elements(fastinfo) eq 0 then return

xtype=1
fileptr={xval:fastinfo.x,spe:fastinfo.y,speerror:sqrt(fastinfo.y),ROIAzimuth:list.AzIntParam[3:4],ytitle:'Intensity (a.u.)',$
            xtitle:ChiXtitle(xtype),title:'Explorative '+(list.fastmodesuperimposed?"superimposed":"average"),xtype:xtype,$
            xrange:lonarr(2),yrange:lonarr(2)}
    
MainUpdateCManual,'c:\dummy.chi',list.top,list.ev,fileptr=fileptr
end;pro exportfastprofile
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ROIThresholdPattern,list
if list.FASTMODE ne 1 then return

; Get ROI map
FastModelBegin,list,fastinfo
if n_elements(fastinfo) eq 0 then return
statgroupi={may:0.,miy:0.}
MakeDatablock,list,block,nrgroups,dim1,dim2,/transform,statgroupi=statgroupi
FastModelEnd,list,fastinfo

; Get threshold
thres='> 50%'
thres=PromptNumber(thres,list.top,'Pixel intensity (value or %):')
p=strpos(thres,'%')
if p ne -1 then begin
    thres=strmid(thres,0,p)
    ind=where(byte(thres) ge 48 and byte(thres) le 57,ct)
    if ct eq 0 then return
    thres=strmid(thres,0,ind[0])+stringr(statgroupi.miy+(statgroupi.may-statgroupi.miy)*float(strmid(thres,ind[0]))/100)
endif

; Plot threshold
ind=compare(block,thres,thresvalue)
Wset,list.outdrawid

case list.SCANDIM of
1:  begin
    SetSysVar,list.sysvar2
    if list.TRAN then plots,[thresvalue,thresvalue],!y.crange,/data $
    else plots,!x.crange,[thresvalue,thresvalue],/data
    endcase
2:  begin
    tvscl,rebin(reform(ind[0,*,*]),dim1*list.nmag,dim2*list.nmag,/sample),/order
    endcase
endcase

; File indices
ind=where(ind,n)
if n eq 0 then begin
  printw,list.prid,'No pixels below this threshold'
  return
endif

case list.SCANDIM of
1:  begin
    if (list.TRAN and list.REV2) or (~list.TRAN and list.REV1) then ind = list.N-1-ind
    endcase
2:  begin
    if list.TRAN then begin
        ncol=list.CR[1]
        nrow=list.CR[0]
    endif else begin
        ncol=list.CR[0]
        nrow=list.CR[1]
    endelse
    ind = ConvertFileNumber(ind mod ncol,ind/ncol,ncol,nrow,list.TRAN,list.REV1,list.REV2,list.REV3)
    endcase
endcase

; Indices that are already processed
ptrCur=(*(*list.ptrModel).ptrSum).next
pblock=(*ptrCur).sum
s=size(*pblock)
if s[0] eq 2 then ny=(s[2]-1)<list.nrtif $
else ny=1

ind2=where(ind lt ny,n)
if n eq 0 then begin
  printw,list.prid,'No pixels that are already processed.'
  return
endif
ind=ind[ind2]

; Combine profiles 
if n eq 1 then Profiley=(*pblock)[*,ind[0]] $
else begin
  if list.fastmodesuperimposed then Profiley=max((*pblock)[*,ind+1],dim=2) $
  else Profiley=total((*pblock)[*,ind+1],2)/n
endelse

xtype=1
fileptr={xval:(*pblock)[*,0],spe:Profiley,speerror:sqrt(Profiley),ROIAzimuth:list.AzIntParam[3:4],ytitle:'Intensity (a.u.)',$
            xtitle:ChiXtitle(xtype),title:'Explorative ROI '+(list.fastmodesuperimposed?"superimposed":"average"),xtype:xtype,$
            xrange:lonarr(2),yrange:lonarr(2)}

MainUpdateCManual,'c:\dummy.chi',list.top,list.ev,fileptr=fileptr

end;pro ROIThresholdPattern
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FILESELECT_PROCESS_EVENTS,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,/INPUT_FOCUS
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL', 'EXPOSE', 'KEY', 'CKEY', 'WHEEL' ]
thisEvent = possibleEventTypes[ev.type]

IF thisEvent eq 'MOTION' THEN RETURN

Widget_Control, ev.top, get_UValue=list
if list.update eq 0 then return

if thisEvent eq 'KEY' then begin
  IF ~ev.press then return
  if ev.ch eq 13 then ROIThresholdPattern,list
  return
endif

if thisEvent ne 'DOWN' THEN return

CATCH, Error_status
IF Error_status NE 0 THEN begin
    printw,list.prid,"Error in selecting the file..."
    return
endif

case list.SCANDIM of
1:    begin
    WSet, list.outdrawid
    SetSysVar,list.sysvar2
    coord = Convert_Coord(ev.x, ev.y, /Device, /To_Data)
    if list.TRAN then filei=coord[1] $
    else filei=coord[0]
    
    xrange=bplotgetxrange(list,list.N)
    if (list.TRAN and list.REV2) or (~list.TRAN and list.REV1) then xrange=reverse(xrange)
    bin=(xrange[1]-xrange[0])/(list.N-1)
    if bin ne 0 then filei=(filei-xrange[0])/bin
    filei=0>round(filei)<(list.N-1)

    endcase
2:    begin
    if list.TRAN then begin
        ncol=list.CR[1]
        nrow=list.CR[0]
    endif else begin
        ncol=list.CR[0]
        nrow=list.CR[1]
    endelse
    if list.FASTMODE eq 2 then begin
        ptrCur=(*(*list.ptrModel).ptrSum).next
        if not ptr_valid(ptrCur) then return
        s=size(*(*ptrCur).sum)
        evx=ev.x/s[1]
        evy=ev.y/s[2]
    endif else begin
        evx=ev.x
        evy=ev.y
    endelse
    nr=CurPos(evx,evy,ncol,nrow,1./list.nmag,1./list.nmag,XX=XX,YY=YY)
    filei=ConvertFileNumber(XX,YY,ncol,nrow,list.TRAN,list.REV1,list.REV2,list.REV3)
    endcase
endcase

ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
widget_control,ID,SET_LIST_SELECT=filei

MainUpdateDo,list,create_struct(ev,'index',filei)

end;pro FILESELECT_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ROIFASTMODE_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,/INPUT_FOCUS
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL', 'EXPOSE', 'KEY', 'CKEY', 'WHEEL' ]
thisEvent = possibleEventTypes[ev.type]

IF thisEvent eq 'MOTION' THEN RETURN

if thisEvent eq 'CKEY' THEN begin
    IF ~ev.press then return
    ROIFASTMODE_DRAWBOX,ev
    return
endif

widget_control,ev.top,get_uvalue=list
ptrCur=(*(*list.ptrModel).ptrSum).next
if ~PTR_VALID(ptrCur) then return

if thisEvent eq 'KEY' THEN begin
  IF ~ev.press then return
  if ev.ch eq 13 then exportfastprofile,list
  return
endif

IF thisEvent NE 'DOWN' THEN RETURN
possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
bPress = possibleButtons[ev.press]

list.xs=ev.x
list.ys=ev.y

if bpress eq 'LEFT' then begin
  WSet, list.outdrawidsum
  
  SetSysVar,list.sysvar
  coords = Convert_Coord(ev.x, ev.y, /Device, /To_Data)
  x1 = !X.CRange[0] > coords[0] < !X.CRange[1]
  ROIy=!Y.CRange
  ROIy=ROIy[sort(ROIy)]
  y1 = ROIy[0] > coords[1] < ROIy[1]
  PlotS, [x1, x1], ROIy
  
  h=!D.x_size
  v=!D.y_size
  WSet, list.pixindexsum
  Device, Copy = [0,0, h,v, 0,0,list.outdrawidsum]
  
  Widget_Control, ev.id, Event_Pro='ROIFASTMODE_DRAWBOX', draw_motion=1 
endif else $
  Widget_Control, ev.id, Event_Pro='ROIFASTMODE_ZOOMBOX', draw_motion=1
  
Widget_Control, ev.top, Set_UValue=list

end;pro ROIFASTMODE_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ROIFASTMODE_ZOOMBOX,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]

WSet, list.outdrawidsum
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.pixindexsum]
x = [list.xs, ev.x]
y = [list.ys, ev.y]

GetCRange,ROIx=ROIx,ROIy=ROIy

IF thisEvent EQ 'UP' THEN BEGIN
    Widget_Control, ev.id,Event_Pro='ROIFASTMODE_PROCESS_EVENTS'

    IF list.xs GT ev.x THEN x = [ev.x, list.xs]
    IF list.ys GT ev.y THEN y = [ev.y, list.ys]
    SetSysVar,list.sysvar
    coords = Convert_Coord(x, y, /Device, /To_Data)
    x1 = ROIx[0] > coords[0,0] < ROIx[1]
    x2 = ROIx[0] > coords[0,1] < ROIx[1]
    y1 = ROIy[0] > coords[1,0] < ROIy[1]
    y2 = ROIy[0] > coords[1,1] < ROIy[1]
    y1 ^= 2
    y2 ^= 2

    list.xrange=[x1,x2]
    list.yrange=[y1,y2]
    
    RefreshDisplayBatch,list
    Widget_Control, ev.top, Set_UValue=list
    widget_control,ev.id,/INPUT_FOCUS
    widget_control,ev.id,/clear_events
    return
endif

SetSysVar,list.sysvar
coords = Convert_Coord(x, y, /Device, /To_Data)
x1 = ROIx[0] > coords[0,0] < ROIx[1]
x2 = ROIx[0] > coords[0,1] < ROIx[1]
y1 = ROIy[0] > coords[1,0] < ROIy[1]
y2 = ROIy[0] > coords[1,1] < ROIy[1]
WSet, list.outdrawidsum
PlotS, [x1, x1], [y1,y2],Linestyle=1
PlotS, [x2, x2], [y1,y2],Linestyle=1
PlotS, [x1,x2],[y1, y1],Linestyle=1
PlotS, [x1,x2],[y2, y2],Linestyle=1
Widget_Control, ev.top, Set_UValue=list
end;pro ROIFASTMODE_ZOOMBOX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ROIFASTMODE_DRAWBOX,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL', 'EXPOSE', 'KEY', 'CKEY', 'WHEEL' ]
thisEvent = possibleEventTypes[ev.type]

WSet, list.outdrawidsum
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.pixindexsum]

x = [list.xs, ev.x]
y = [list.ys, ev.y]
SetSysVar,list.sysvar
coords = Convert_Coord(x, y, /Device, /To_Data)
x1 = !X.CRange[0] > coords(0,0) < !X.CRange[1]
x2 = !X.CRange[0] > coords(0,1) < !X.CRange[1]

IF thisEvent EQ 'UP' or thisEvent eq 'CKEY' THEN BEGIN
    Widget_Control, ev.id, Event_Pro='ROIFASTMODE_PROCESS_EVENTS'

    ptrCur=(*(*list.ptrModel).ptrSum).next
    x=(*((*ptrCur).sum))[*,0]

    if thisEvent eq 'CKEY' then begin
        x1=(*ptrCur).cor[0]
        x2=(*ptrCur).cor[1]
        case ev.key of
        5:    begin
            x1--
            x2--
            endcase
        6:    begin
            x1++
            x2++
            endcase
        7:    begin
            x1+=10
            x2+=10
            endcase
        8:    begin
            x1-=10
            x2-=10
            endcase
        else: return
        endcase
    endif else begin
        tt=[x1<x2,x1>x2]
        t=min(abs(x-tt[0]),x1)
        t=min(abs(x-tt[1]),x2)
    endelse

    (*ptrCur).cor=[x1<x2,x1>x2]
    tt=x[[x1<x2,x1>x2]]/180.*!pi
    (*ptrCur).mask[3:4]=reverse(braggTtoX(tt,/onlyd,list.lambda))
    list.groupi=0
    list.plottitle=getcurrentlabel(list)
    ID2=widget_info(ev.top,FIND_BY_UNAME='plottitle')
    if widget_info(ID2,/valid) then widget_control,ID2,set_value=list.plottitle

    RefreshDisplayBatch,list
    Widget_Control, ev.top, Set_UValue=list
    widget_control,ev.id,/INPUT_FOCUS
    widget_control,ev.id,/clear_events

    if list.update then begin
        ; Update 2D
        widget_control,list.top,get_uvalue=toplist
        toplist.ttmellipse=mean(tt)
        widget_control,list.top,set_uvalue=toplist
        RefreshDisplay,list.top,toplist
        TTWindowUpdate,toplist.ttmellipse,toplist,list.top
        
        ; Update 1D
        if ptr_valid(toplist.CHIchild) and toplist.CrossUpdate then begin
            val=toplist.ttmellipse*180./!pi
            for i=0l,n_elements(*toplist.CHIchild)-1 do $
                SetAnySelect,val,1,(*toplist.CHIchild)[i]
        endif
    endif

    return
endif

list.xd=ev.x
list.yd=ev.y

ROIy=!Y.CRange
ROIy=ROIy[sort(ROIy)]
y2 = ROIy[0] > coords(1,1) < ROIy[1]
WSet, list.outdrawidsum
Arrow, x1, y2, x2, y2, /Data, /Solid, HSize=12
PlotS, [x2, x2], ROIy

Widget_Control, ev.top, Set_UValue=list
end;pro ROIFASTMODE_DRAWBOX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EndBatchSequenceSave ,list, manual=manual

if list.nfiles eq 0 then begin
    printw,list.prid,'Nothing processed yet.'
    return,''
endif

T=systime(1)

ext='.xdi'
xdipath=list.xdipath+list.xdifile+ext
Result = file_search(xdipath,count=count)
if (count ne 0) and (list.AUTO eq 0) then begin
       Result = DIALOG_MESSAGE(xdipath+' is already present: overwrite it?',/Question)
       if Result eq 'No' then begin
            xdipath=''
         xdipath=DIALOG_PICKFILE(path=list.xdipath,file=list.xdifile+ext,filter='*'+ext,title='Save as...')
       endif
endif
if xdipath eq '' then begin
       printw,list.prid,'Reduced data dynamic list: Not Saved'
endif else begin
        if list.FASTMODE eq 2 then begin
            error=CutPath(xdipath,path=path,file=file)
            jpgpath = path+file+'.bmp'
            Wset,list.outdrawid
            write_bmp,jpgpath,tvrd(true=1),/rgb
            printw,list.prid,'Reduced data dynamic list: Saved'
               printw,list.prid,jpgpath+' :'+string(systime(1)-T)+' sec'
            return,''
        endif

        if list.FASTMODE eq 1 then begin
            error=CutPath(xdipath,path=path,file=file)
            
            MakeDatablock,list,block,/FASTMODE
            if not keyword_set(manual) then $
                result=WriteChi(path+file+'.tiff',block,1,'...',$
                                list.prid,list.AzIntParam[3:4],berror=list.AzIntInfo.calcerror)

            FastModelBegin,list,fastinfo
        endif

           result=SaveXDI(list,xdipath)

        if list.FASTMODE eq 1 then FastModelEnd,list,fastinfo

           if result then begin
               printw,list.prid,'Reduced data dynamic list: Saved'
               printw,list.prid,xdipath+' :'+string(systime(1)-T)+' sec'
           endif else begin
               printw,list.prid,'Reduced data dynamic list: Not Saved'
               xdipath = ''
           endelse
endelse

return,    xdipath
end;function EndBatchSequenceSave
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EndBatchSequence ,list,ev,prompt=prompt

; ----Save data to *.xdi file----
if keyword_set(prompt) then Resultp = DIALOG_MESSAGE('Save reduced data dynamic list?'+$
          ' When pressing Go this list will be flushed!',/Question,/DEFAULT_NO,/cancel)$
else Resultp='Yes'

if Resultp eq 'Cancel' then return, 0B

printw,list.prid,'----Terminate----'
T=systime(1)
if Resultp eq 'Yes' then begin
    xdipath=EndBatchSequenceSave(list)
    T=systime(1)
endif else begin
    printw,list.prid,'Reduced data dynamic list: Not Saved'
    xdipath = ''
endelse

; ----Save Display to .jpg----
if list.PLOT and (xdipath ne '') then begin
error=CutPath(xdipath,path=jpgpath,file=file)
jpgpath = jpgpath+file+'.bmp'
;write_jpeg,jpgpath,tvrd(true=3),true=3,QUALITY=100
write_bmp,jpgpath,tvrd(true=1),/rgb
printw,list.prid,jpgpath+' saved:'+string(systime(1)-T)+' sec'
T=systime(1)
endif

; ----Update Window----
ID=widget_info(ev.top,FIND_BY_UNAME='Go')
widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='Exit')
widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='Pause')
widget_control,ID,sensitive=0
widget_control,ID,set_value='Pause'
ID=widget_info(ev.top,FIND_BY_UNAME='Abort')
widget_control,ID,sensitive=0
XRD_BP_PROCESSSENSITIVE,ev.top,sensitive=1
;ID=widget_info(ev.top,FIND_BY_UNAME='Maskbase')
;widget_control,ID,sensitive=1
widget_control,ev.id,sensitive=0
list.timer=0B
list.pause=0B
list.title='Batch Processing: '
widget_control,ev.top,TLB_SET_TITLE=list.title+(list.char)[list.chari]

; ----Finish----
ETIME=systime(1)-list.TotTime[0]
ConvertTime,ETIME,timeaf
printw,list.prid,'Total time:'+string(ETIME)+timeaf
printw,list.prid,'-----FINISHED-----'
return,1B

end;function EndBatchSequence
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FastModelBegin,list,fastinfo

if list.FASTMODE ne 1 then return

; Point to Fast model (type1)
ptrCur=(*(*list.ptrModel).ptrSum).next
if ~PTR_VALID(ptrCur) then return

; Slow but better version:
; 
;    if list.FORCE then begin
;        keepforce=[(*ptrCur).code[0],list.backmod,(*(*list.ptrModel).ptrFit).backparamIO[0]]
;        cb=0>(*ptrCur).cor[0,0]<(list.AzIntParam[2]-1)
;        ce=0>(*ptrCur).cor[1,0]<(list.AzIntParam[2]-1)
;        if cb ne ce then begin
;            ; Ortpol with degree 1
;            (*ptrCur).code[0]=1
;            list.backmod=0
;            (*(*list.ptrModel).ptrFit).backparamIO[0]=1
;        endif else begin
;            (*ptrCur).code[0]=0
;            list.backmod=-1
;            (*(*list.ptrModel).ptrFit).backparamIO[0]=-1
;        endelse
;    endif
;    
;    pblock=(*ptrCur).sum
;    Profilex=(*pblock)[*,0]
;    Profiley=(*pblock)[*,1:*]
;    s=size(Profiley)
;    if s[0] eq 2 then begin
;        ny=s[2]
;        Profiley=total(Profiley,2)/ny
;    endif else ny=1
;    
;    ; Modify pblock to normal mode scan
;    ny<=list.nrtif
;    for j=1l,ny do begin
;        Process1DROI,list,0,Profilex,(*pblock)[*,j],ptrCur,x,y,yb
;        
;        ; Add to matrix for PCA
;        if j eq 1 then keepy=y else keepy=[[keepy],[y]]
;    endfor


; Mean pattern: sqrt(I) vs. two-theta
pblock=(*ptrCur).sum
Profilex=(*pblock)[*,0]
s=size(*pblock)
if s[0] eq 2 then begin
    ny=(s[2]-1)<list.nrtif
    if ny eq 1 then Profiley=(*pblock)[*,1] $
    else if list.fastmodesuperimposed then Profiley=max((*pblock)[*,1:ny],dim=2) $
    else Profiley=total((*pblock)[*,1:ny],2)/ny
endif else begin
    ny=1
    Profiley=(*pblock)[*,1]
endelse

; ROI in bins
tmp=BraggDtoX((*ptrCur).mask[3:4,0],list.lambda,/onlyt,/angledeg)
t=min(abs(Profilex-tmp[0]),x1);convert from two-theta to index in Profilex
t=min(abs(Profilex-tmp[1]),x2)
(*ptrCur).cor[*,0]=[x1<x2,x1>x2]

; Select ROI
n=n_elements(Profilex)-1
cb=0>(*ptrCur).cor[0,0]<n
ce=0>(*ptrCur).cor[1,0]<n
nx=ce-cb+1
keepy=(*pblock)[cb:ce,1:ny]

; Subtract background
if list.FORCE and nx gt 2 then begin
    ; Simple linear
    y1=rebin(keepy[0,*],nx,ny,/sample)
    dy=rebin(keepy[ce-cb,*],nx,ny,/sample)-y1
    keepy-=rebin(lindgen(nx),nx,ny,/sample)*temporary(dy)/(ce-cb)+temporary(y1)
    keepy>=0
endif

if list.groupi eq 0 then begin
    sum=total(keepy,1)
endif else begin
    list.groupi<=nx
    sum=mPCA(keepy,error=error,/standardize)
    if error then sum=total(keepy,1) $
    else sum=sum[list.groupi-1,*]
endelse
sum=reform(sum)

pkeep=(*ptrCur).sum
(*ptrCur).sum=ptr_new(sum)

fastinfo={pkeep:pkeep,ptrCur:ptrCur,x:Profilex,y:Profiley}

; Force background
;    if list.FORCE then begin
;        (*ptrCur).code[0]=keepforce[0]
;        list.backmod=keepforce[1]
;        (*(*list.ptrModel).ptrFit).backparamIO[0]=keepforce[2]
;    endif

end;pro FastModelBegin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FastModelEnd,list,fastinfo
if n_elements(fastinfo) eq 0 then return

ptr_free,(*fastinfo.ptrCur).sum
(*fastinfo.ptrCur).sum=fastinfo.pkeep

end;pro FastModelEnd
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PauseBatchExe,list,ev

list.TotTime[1]=systime(1)

; ----Update Window----
XRD_BP_PROCESSSENSITIVE,ev.top,sensitive=1,/pause
widget_control,ev.id,set_value='Resume'
list.timer=0B
list.Pause=1B

return,1B
end;function PauseBatchExe
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ResumeBatchExe ,list,ev

list.TotTime[0]=list.TotTime[0]+systime(1)-list.TotTime[1]
list.TotTime[1]=0

; ----Check if changes are OK for processed data so far----
case list.SCANDIM of
1:    NNew=list.N
2:    NNew=list.CR[0]*list.CR[1]
endcase
if list.nrtif ne NNew then begin
    list.nrtif=NNew
    ; TODO: resize everything
    if list.nFiles gt list.nrtif then begin
        printw,list.prid,'Dimensions are too large!'
        return,0B
    endif
endif

; ----Update Window----
XRD_BP_PROCESSSENSITIVE,ev.top,sensitive=0,/pause
widget_control,ev.id,set_value='Pause'
list.timer=1B
list.Pause=0B
ID=widget_info(ev.top,FIND_BY_UNAME='TimerLabel')
; ----Make output window----
result=PrepOut(list,ev)
; ----Send timer events----
list.ACQSTAMP=systime(1)
;WIDGET_CONTROL, ev.top, TIMER=0 ; Hit files timer
WIDGET_CONTROL, ID, TIMER=0 ; Hit process timer

return,1B
end;function ResumeBatchExe
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro UpdateCText,list,ev,val
values=['C1Text0','C1Text1','CText0','CText1','CText2','CText3']
ind=where(values eq val,count)
A=bytarr(3)
if count eq 0 then A[*]=1B else begin
    case ind[0] of
    0:A[0]=1B
    1:A[0]=1B
    2:A[1]=1B
    3:A[1]=1B
    4:A[2]=1B
    5:A[2]=1B
    endcase
endelse

if A[0] then begin
    if list.TRAN then begin
       BOOL=list.REV2
       sign=list.Vsign
    endif else begin
       BOOL=list.REV1
       sign=list.Hsign
    endelse
    if sign eq -1 then BOOL=BOOL ne 1
    if BOOL then list.coord1=list.coord1[reverse(sort(list.coord1))]$; Large to small
    else list.coord1=list.coord1[sort(list.coord1)]; Small To Large

    ID=widget_info(ev.top,FIND_BY_UNAME='C1Text0')
          widget_control,ID,set_value=stringr(list.coord1[0])
    ID=widget_info(ev.top,FIND_BY_UNAME='C1Text1')
          widget_control,ID,set_value=stringr(list.coord1[1])
endif

if A[1] then begin
    if list.TRAN then begin
       BOOL=list.REV2
       sign=list.Vsign
    endif else begin
       BOOL=list.REV1
       sign=list.Hsign
    endelse
    if sign eq -1 then BOOL=BOOL ne 1
    vector=list.coord[*,0]
    if BOOL then list.coord[*,0]=vector[reverse(sort(vector))]$; Large to small
    else list.coord[*,0]=vector[sort(vector)]; Small To Large
    ID=widget_info(ev.top,FIND_BY_UNAME='CText0')
          widget_control,ID,set_value=stringr(list.coord[0])
    ID=widget_info(ev.top,FIND_BY_UNAME='CText1')
          widget_control,ID,set_value=stringr(list.coord[1])
endif

if A[2] then begin
    if list.TRAN then begin
       BOOL=list.REV1
       sign=list.Hsign
    endif else begin
       BOOL=list.REV2
       sign=list.Vsign
    endelse
    if sign eq -1 then BOOL=BOOL ne 1
    vector=list.coord[*,1]
    if BOOL then list.coord[*,1]=vector[reverse(sort(vector))]$; Large to small
    else list.coord[*,1]=vector[sort(vector)]; Small To Large
    ID=widget_info(ev.top,FIND_BY_UNAME='CText2')
          widget_control,ID,set_value=stringr(list.coord[2])
    ID=widget_info(ev.top,FIND_BY_UNAME='CText3')
          widget_control,ID,set_value=stringr(list.coord[3])
endif

end;pro UpdateCText
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SSession_XRD_BP,ev,error=error,auto=auto

widget_control,/hourglass
widget_control,ev.top,get_uvalue=list
error=1b

if list.MSESSION eq 1 then begin
    datpath=list.xdipath+list.xdifile+'.dat'

    Result = file_search(datpath,count=count)
    if (count ne 0) and (list.AUTO eq 0) and (~keyword_set(auto)) then begin
           Result = DIALOG_MESSAGE(datpath+' is already present: overwrite it?',/Question)
           if Result eq 'No' then begin
                error=CutPath(datpath,path=path,file=file,ext=ext)
                datpath=''
             datpath=DIALOG_PICKFILE(path=path,file=file+ext,filter='*.dat',title='Save as...')
             if datpath eq '' then return
           endif
    endif
endif else datpath='idlsave_XRD_BP.dat'

CATCH, Error_status
IF Error_status NE 0 THEN begin
    printw,list.prid,datpath+' NOT saved.'
    return
endif

if FILE_TEST(datpath) then begin
    datpath2=datpath+'.tmp'
    SAVE,list,FILENAME=datpath2
    file_delete,datpath,/quiet
    file_move,datpath2,datpath
endif else SAVE,list,FILENAME=datpath

error=0b

printw,list.prid,datpath+' saved.'
end;pro SSession_XRD_BP
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_XRD_BP, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
if winidvalid(list.pixindexsum) then wdelete,list.pixindexsum

ptr_free,list.PTR_files
heap_free,list

;result=UpdateINI('$CrossUpdate:',{n:list.update})
end;pro CleanUp_XRD_BP
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Refresh_XRD_BP,ev

widget_control,ev.top,get_uvalue=list

; Set it in the Paused state when saved during processing state!
if (list.timer eq 1) and (list.pause eq 0) then begin
    list.timer=0B
    list.pause=1B
endif
widget_control,ev.top,TLB_SET_TITLE=list.title+(list.char)[list.chari]
ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
    if ptr_valid(list.PTR_files) then begin
        widget_control,ID,set_value=*list.PTR_files
        n=n_elements(*list.PTR_files)
    endif else n=0l
ID=widget_info(ev.top,FIND_BY_UNAME='FilesText')
    widget_control,ID,set_value=string(n,format='(I0," files")')
XRD_BP_PROCESSSENSITIVE,ev.top,sensitive=0
XRD_BP_PROCESSSENSITIVE,ev.top,sensitive=list.timer ne 1,pause=list.pause
ID=widget_info(ev.top,FIND_BY_UNAME='FilterText')
if list.filtermethod eq 2 then begin
    widget_control,ID,sensitive=0
endif else begin
    widget_control,ID,set_value=(*list.filter),sensitive=1
endelse
ID=widget_info(ev.top,FIND_BY_UNAME='PathText')
    widget_control,ID,set_value=list.path
;ID=widget_info(ev.top,FIND_BY_UNAME='Maskbase')
;    widget_control,ID,sensitive=(list.timer XOR list.pause) ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='MaskText')
    widget_control,ID,set_value=list.mskpath+list.mskfile
ID=widget_info(ev.top,FIND_BY_UNAME='OutFileText')
    widget_control,ID,set_value=list.xdifile
ID=widget_info(ev.top,FIND_BY_UNAME='OutDirText')
    widget_control,ID,set_value=list.xdipath
ID=widget_info(ev.top,FIND_BY_UNAME='SecSlider')
    widget_control,ID,set_value=list.sec
ID=widget_info(ev.top,FIND_BY_UNAME='Go')
    widget_control,ID,sensitive=(list.timer XOR list.pause) ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='Exit')
    widget_control,ID,sensitive=(list.timer XOR list.pause) ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='Abort')
    widget_control,ID,sensitive=list.pause
ID=widget_info(ev.top,FIND_BY_UNAME='Pause')
    widget_control,ID,sensitive=list.timer or list.pause
ID=widget_info(ev.top,FIND_BY_UNAME='Pause')
    if (list.timer eq 0) and (list.pause eq 1) then pval='Resume'$
       else pval='Pause'
    widget_control,ID,set_value=pval

ID=widget_info(ev.top,FIND_BY_UNAME='NBase')
    widget_control,ID,sensitive=list.SCANDIM eq 1
ID=widget_info(ev.top,FIND_BY_UNAME='ColBase')
    widget_control,ID,sensitive=list.SCANDIM eq 2
ID=widget_info(ev.top,FIND_BY_UNAME='RowBase')
    widget_control,ID,sensitive=list.SCANDIM eq 2
ID=widget_info(ev.top,FIND_BY_UNAME='REV1')
    widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
ID=widget_info(ev.top,FIND_BY_UNAME='REV2')
    widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
ID=widget_info(ev.top,FIND_BY_UNAME='REV3')
    widget_control,ID,sensitive=(list.SCANDIM eq 2)
ID=widget_info(ev.top,FIND_BY_UNAME='HORIZ')
    widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
ID=widget_info(ev.top,FIND_BY_UNAME='VERT')
    widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
ID=widget_info(ev.top,FIND_BY_UNAME='NText')
    widget_control,ID,set_value=stringr(list.N-1)
ID=widget_info(ev.top,FIND_BY_UNAME='ColText')
    widget_control,ID,set_value=stringr(list.CR[0]-1)
ID=widget_info(ev.top,FIND_BY_UNAME='RowText')
    widget_control,ID,set_value=stringr(list.CR[1]-1)
ID=widget_info(ev.top,FIND_BY_UNAME='CT')
    widget_control,ID,set_value=stringr(list.ct)
if list.TRAN then ar=['MotV','MotH']$
          else ar=['MotH','MotV']
    ID=widget_info(ev.top,FIND_BY_UNAME='NLabel')
       widget_control,ID,set_value=ar[0]
    ID=widget_info(ev.top,FIND_BY_UNAME='ColLabel')
       widget_control,ID,set_value=ar[0]
    ID=widget_info(ev.top,FIND_BY_UNAME='RowLabel')
       widget_control,ID,set_value=ar[1]
    UpdateCText,list,ev,'dummy'
ID=widget_info(ev.top,FIND_BY_UNAME='Unit1')
    widget_control,ID,SET_DROPLIST_SELECT=list.unit1
ID=widget_info(ev.top,FIND_BY_UNAME='Unit2')
    widget_control,ID,SET_DROPLIST_SELECT=list.unit[0]
ID=widget_info(ev.top,FIND_BY_UNAME='Unit3')
    widget_control,ID,SET_DROPLIST_SELECT=list.unit[1]
ID=widget_info(ev.top,FIND_BY_UNAME='First vertical')
    widget_control,ID,set_button=list.TRAN
ID=widget_info(ev.top,FIND_BY_UNAME='First horizontal')
    widget_control,ID,set_button=list.TRAN ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='Right-Left')
    widget_control,ID,set_button=list.REV1
ID=widget_info(ev.top,FIND_BY_UNAME='Left-Right')
    widget_control,ID,set_button=list.REV1 ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='Bottom-Top')
    widget_control,ID,set_button=list.REV2
ID=widget_info(ev.top,FIND_BY_UNAME='Top-Bottom')
    widget_control,ID,set_button=list.REV2 ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='Jump')
    widget_control,ID,set_button=list.REV3 ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='Continue')
    widget_control,ID,set_button=list.REV3
ID=widget_info(ev.top,FIND_BY_UNAME='L--R')
    if list.Hsign eq 1 then bval='L-->R' else bval='L<--R'
    widget_control,ID,set_value=bval
ID=widget_info(ev.top,FIND_BY_UNAME='T--B')
    if list.Vsign eq 1 then bval='T-->B' else bval='T<--B'
    widget_control,ID,set_value=bval
ID=widget_info(ev.top,FIND_BY_UNAME='Update Main Window')
    widget_control,ID,set_button=list.update
ID=widget_info(ev.top,FIND_BY_UNAME='checkbutton')
    widget_control,ID,set_button=list.nocheck
ID=widget_info(ev.top,FIND_BY_UNAME='Autobutton')
    widget_control,ID,set_button=list.AUTO
ID=widget_info(ev.top,FIND_BY_UNAME='Outbutton')
    widget_control,ID,set_button=list.SR eq 1
ID=widget_info(ev.top,FIND_BY_UNAME='Outbutton2')
    widget_control,ID,set_button=list.SR eq 2
;ID=widget_info(ev.top,FIND_BY_UNAME='MSESSIONbutton')
;    widget_control,ID,set_button=list.MSESSION
ID=widget_info(ev.top,FIND_BY_UNAME='Plotbutton')
    widget_control,ID,set_button=list.PLOT
ID=widget_info(ev.top,FIND_BY_UNAME='Debugbutton')
    widget_control,ID,set_button=list.DEBUG
ID=widget_info(ev.top,FIND_BY_UNAME='Forcebutton')
    widget_control,ID,set_button=list.FORCE
ID=widget_info(ev.top,FIND_BY_UNAME='Supbutton')
    widget_control,ID,set_button=list.fastmodesuperimposed
ID=widget_info(ev.top,FIND_BY_UNAME='Smallbutton2D')
    widget_control,ID,set_button=list.SMALL2D
ID=widget_info(ev.top,FIND_BY_UNAME='Smallbutton1D')
    widget_control,ID,set_button=list.SMALL1D
ID=widget_info(ev.top,FIND_BY_UNAME='NORMROIsbutton')
    widget_control,ID,set_button=list.NORMROIs
ID=widget_info(ev.top,FIND_BY_UNAME='Scanbutton')
    widget_control,ID,set_button=list.SCANDIM eq 1
ID=widget_info(ev.top,FIND_BY_UNAME='Mapbutton')
    widget_control,ID,set_button=list.SCANDIM eq 2
case list.fastmode of
0: ID=widget_info(ev.top,FIND_BY_UNAME='fullmode')
1: ID=widget_info(ev.top,FIND_BY_UNAME='fastmode')
2: ID=widget_info(ev.top,FIND_BY_UNAME='gridmode')
endcase
widget_control,ID,set_button=1

; ----Update plot window----
result=PrepOut(list,ev)
RefreshDisplayBatch,list
widget_control,ev.top,set_uvalue=list
end;pro Refresh_XRD_BP
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RSession_XRD_BP,ev

widget_control,/hourglass
widget_control,ev.top,get_uvalue=list2
; ---- Select session----
if list2.MSESSION eq 1 then begin
    datpath=DIALOG_PICKFILE(path=list2.xdipath,file=list2.xdifile+'.dat',$
         filter='*.dat',title='Load session....')
endif else datpath='idlsave_XRD_BP.dat'
if datpath eq '' then return

; ----Restore list----
RESTORE,FILENAME=datpath

; ----Copy fields except widget IDs----
RestoreStruc,list2,list,exclude=['outdraw','outdrawid','outdrawsum','outdrawidsum','pixindexsum','top','outid1','outid2','outid3','ev','prid']

; ----Convert from old to new----
if ptr_valid(list2.PTR_list) then begin
    ; Old version: 5.2.17.3 or earlier
    ptr_free,list2.PTR_files
    list2.PTR_files=ptr_new(PrintLL(list2.PTR_list))
    DestroyLL,list2.PTR_list,-1
endif

list2.nmult = 1 ; Old version: before 6.3.1.1

; ----Save structure----
widget_control,ev.top,set_uvalue=list2

; ----Refresh display----
Refresh_XRD_BP,ev

if list2.PLOT then printw,list2.prid,datpath+' restored.'

end;pro RSession_XRD_BP
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRD_BP_PROCESSSENSITIVE,top,pause=pause,sensitive=sensitive
if keyword_set(pause) then unames=['Initbase'] $
else unames=['Initbase','modebase']

for i=0l,n_elements(unames)-1 do begin
    ID=widget_info(top,FIND_BY_UNAME=unames[i])
    widget_control,ID,sensitive=sensitive
endfor
end;pro XRD_BP_PROCESSSENSITIVE
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SelectGroupMainBP,list
if ~list.update then return
widget_control,list.top,get_uvalue=toplist
toplist.maskP=list.groupi
widget_control,list.top,set_uvalue=toplist
RefreshDisplay,list.top,toplist
end;pro SelectGroupMainBP
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRD_BP_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list

case widget_info(ev.id,/type) of
4:  begin;draw
    end
0:  begin;Base (File Search Timer Event)

    ; Stop timer
    if list.timer eq 0 then return
    
    ; Next timer event when still processing
    if ptr_valid(list.PTR_files) then $
        if list.nfiles lt n_elements(*list.PTR_files) then begin
            ; Still processing
            WIDGET_CONTROL, ev.id, TIMER=list.sec ; generate files timer
            return
        endif
        
    widget_control,/hourglass
    
    ; Number of files detected
    if ptr_valid(list.PTR_files) then nfilesold=n_elements(*list.PTR_files) $
    else nfilesold=0l
    
    ; Detect new files and generate files timer
    update=list.update
    list.update=0b; Update on file processing, not on loading
    result=DetectFiles2(list,ev,/nodel)
    list.update=update

    ; Number of new files detected
    if ptr_valid(list.PTR_files) then nfilesnew=n_elements(*list.PTR_files) $
    else nfilesnew=0
    
    ; Hit process timer when new files
    if nfilesnew ne nfilesold then begin
        ID=widget_info(ev.top,FIND_BY_UNAME='TimerLabel')
        WIDGET_CONTROL, ID, TIMER=0
    endif

    endcase
1:  begin;Button event
    widget_control,ev.id,get_value=val
    case val of
    'Full processing mode':begin
          list.FASTMODE=~ev.select
          endcase
    'Explorative mode':begin
          list.FASTMODE=ev.select
          endcase
    'Grid mode':begin
          if ev.select then list.FASTMODE=2
          endcase
    'Save Session':begin
          SSession_XRD_BP,ev
          endcase
    'Restore Session':begin
          RSession_XRD_BP,ev
          widget_control,ev.top,get_uvalue=list
          return
          endcase
    'L-->R':   begin
          list.Hsign=-1
          widget_control,ev.id,set_value='L<--R'
          UpdateCText,list,ev,val
          RefreshDisplayBatch,list
          endcase
    'L<--R':   begin
          list.Hsign=1
          widget_control,ev.id,set_value='L-->R'
          UpdateCText,list,ev,val
          RefreshDisplayBatch,list
          endcase
    'T-->B':   begin
          list.Vsign=-1
          widget_control,ev.id,set_value='T<--B'
          UpdateCText,list,ev,val
          RefreshDisplayBatch,list
          endcase
    'T<--B':   begin
          list.Vsign=1
          widget_control,ev.id,set_value='T-->B'
          UpdateCText,list,ev,val
          RefreshDisplayBatch,list
          endcase
    'Pause':   begin
          result=PauseBatchExe(list,ev)
          endcase
    'Resume':  begin
          result=ResumeBatchExe(list,ev)
          endcase
    'Abort':   begin
          result=EndBatchSequence(list,ev,/prompt)
          endcase
    'Save current result':begin
            xdipath=EndBatchSequenceSave(list)
            endcase
    'Save current result (no tiff)':begin
            xdipath=EndBatchSequenceSave(list,/manual)
            endcase
    '+':    begin
            list.nmag++
            result=PrepOut(list,ev)
            RefreshDisplayBatch,list
            endcase
    '-':    begin
            list.nmag=(list.nmag-1)>1
            result=PrepOut(list,ev)
            RefreshDisplayBatch,list
            endcase
    '*':    begin
            list.nmult=(list.nmult+0.05)<1
            printw,list.prid,'Intensity scaling at '+stringr(list.nmult*100)+'%'
            result=PrepOut(list,ev)
            RefreshDisplayBatch,list
            endcase
    '/':    begin
            list.nmult=(list.nmult-0.05)>0
            printw,list.prid,'Intensity scaling at '+stringr(list.nmult*100)+'%'
            result=PrepOut(list,ev)
            RefreshDisplayBatch,list
            endcase
    '<<':    begin
            list.groupi--
            result=PrepOut(list,ev)
            RefreshDisplayBatch,list
            SelectGroupMainBP,list
            endcase
    '>>':    begin
            list.groupi++
            result=PrepOut(list,ev)
            RefreshDisplayBatch,list
            SelectGroupMainBP,list
            endcase
    'Go':   begin
          list.prid=widget_info(ev.top,FIND_BY_UNAME='Status')
          list.TotTime[0]=systime(1)
          printw,list.prid,'----Initialisation----'
          error=BatchInit(list,ev)
          if not error then return
          printw,list.prid,'-----Initialised-----'

          ; ----Prepare Display----
          ID=widget_info(ev.top,FIND_BY_UNAME='Abort')
          widget_control,ID,sensitive=1
          ID=widget_info(ev.top,FIND_BY_UNAME='Exit')
          widget_control,ID,sensitive=0
          ID=widget_info(ev.top,FIND_BY_UNAME='Pause')
          widget_control,ID,sensitive=1
          XRD_BP_PROCESSSENSITIVE,ev.top,sensitive=0
          list.timer=1
          widget_control,ev.id,sensitive=0

          ; ----Send timer events----
          list.ACQSTAMP=systime(1)
          WIDGET_CONTROL, ev.top, TIMER=0 ; Hit files timer (list.sec will be set after first search)
          endcase
    'Exit':    begin
          widget_control,ev.top,/destroy
          return
          endcase
    'File sorting':begin
          MakeBOptionsWindow,ev
          return
          endcase
    'Normalized ROI`s':list.NORMROIs=ev.select
    '2D background on ROI`s':list.SMALL2D=ev.select
    '1D background on ROI`s':list.SMALL1D=ev.select
    'Monitor processing':   list.DEBUG=ev.select
    'Explorative background': begin
                    list.FORCE=ev.select
                    RefreshDisplayBatch,list
          endcase
    'Explorative superimposed profile': begin
          list.fastmodesuperimposed=ev.select
          RefreshDisplayBatch,list
          endcase
    'View progress':begin
;            if ev.select then begin
                ID=widget_info(ev.top,FIND_BY_UNAME='Status')
                list.prid=ID
;            endif else list.prid=0
            list.PLOT=ev.select
            result=PrepOut(list,ev)
            RefreshDisplayBatch,list
          endcase
    'Autosave Session/Point':begin
            list.SR=ev.select
            ID=widget_info(ev.top,FIND_BY_UNAME='Outbutton2')
            widget_control,ID,set_button=0
          endcase
    'Autosave Session/Line':begin
            list.SR=ev.select*2
            ID=widget_info(ev.top,FIND_BY_UNAME='Outbutton')
            widget_control,ID,set_button=0
          endcase
    'Update Main Window':list.update=ev.select
    'Search now':begin
            WIDGET_CONTROL, ev.top, TIMER=0 ; Hit files timer
            endcase
    "Don't check files":list.nocheck=ev.select
    'Overwrite existing outputfiles':list.AUTO=ev.select
    'Multiple sessions':list.MSESSION=ev.select
    'First horizontal':   begin
          list.TRAN=ev.select ne 1
          if list.TRAN then ar=['MotV','MotH']$
          else ar=['MotH','MotV']
          ID=widget_info(ev.top,FIND_BY_UNAME='NLabel')
          widget_control,ID,set_value=ar[0]
          ID=widget_info(ev.top,FIND_BY_UNAME='ColLabel')
          widget_control,ID,set_value=ar[0]
          ID=widget_info(ev.top,FIND_BY_UNAME='RowLabel')
          widget_control,ID,set_value=ar[1]
          ID=widget_info(ev.top,FIND_BY_UNAME='REV1')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
          ID=widget_info(ev.top,FIND_BY_UNAME='REV2')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
          ID=widget_info(ev.top,FIND_BY_UNAME='HORIZ')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
          ID=widget_info(ev.top,FIND_BY_UNAME='VERT')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
          UpdateCText,list,ev,val
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          endcase
    'First vertical':   begin
          list.TRAN=ev.select
          if list.TRAN then ar=['MotV','MotH']$
          else ar=['MotH','MotV']
          ID=widget_info(ev.top,FIND_BY_UNAME='ColLabel')
          widget_control,ID,set_value=ar[0]
          ID=widget_info(ev.top,FIND_BY_UNAME='RowLabel')
          widget_control,ID,set_value=ar[1]
          ID=widget_info(ev.top,FIND_BY_UNAME='NLabel')
          widget_control,ID,set_value=ar[0]
          ID=widget_info(ev.top,FIND_BY_UNAME='REV1')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
          ID=widget_info(ev.top,FIND_BY_UNAME='REV2')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
          ID=widget_info(ev.top,FIND_BY_UNAME='HORIZ')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
          ID=widget_info(ev.top,FIND_BY_UNAME='VERT')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
          UpdateCText,list,ev,val
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          endcase
    'Left-Right':   begin
          list.REV1=ev.select ne 1
          UpdateCText,list,ev,val
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          endcase
    'Right-Left':   begin
          list.REV1=ev.select
          UpdateCText,list,ev,val
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          endcase
    'Top-Bottom':   begin
          list.REV2=ev.select ne 1
          UpdateCText,list,ev,val
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          endcase
    'Bottom-Top':   begin
          list.REV2=ev.select
          UpdateCText,list,ev,val
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          endcase
    'Jump':begin
            list.REV3=ev.select ne 1
            result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          endcase
    'Continue':begin
            list.REV3=ev.select
            result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          endcase
    'Line':       begin
          if ev.select then list.SCANDIM=1 else list.SCANDIM=2
          ID=widget_info(ev.top,FIND_BY_UNAME='NBase')
          widget_control,ID,sensitive=list.SCANDIM eq 1
          ID=widget_info(ev.top,FIND_BY_UNAME='ColBase')
          widget_control,ID,sensitive=list.SCANDIM eq 2
          ID=widget_info(ev.top,FIND_BY_UNAME='RowBase')
          widget_control,ID,sensitive=list.SCANDIM eq 2
          ID=widget_info(ev.top,FIND_BY_UNAME='REV1')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
          ID=widget_info(ev.top,FIND_BY_UNAME='REV2')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
          ID=widget_info(ev.top,FIND_BY_UNAME='REV3')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)
          ID=widget_info(ev.top,FIND_BY_UNAME='HORIZ')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
          ID=widget_info(ev.top,FIND_BY_UNAME='VERT')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          updatefilecounter,ev,list
          endcase
    'Map':     begin
          if ev.select then list.SCANDIM=2 else list.SCANDIM=1
          ID=widget_info(ev.top,FIND_BY_UNAME='NBase')
          widget_control,ID,sensitive=list.SCANDIM eq 1
          ID=widget_info(ev.top,FIND_BY_UNAME='ColBase')
          widget_control,ID,sensitive=list.SCANDIM eq 2
          ID=widget_info(ev.top,FIND_BY_UNAME='RowBase')
          widget_control,ID,sensitive=list.SCANDIM eq 2
          ID=widget_info(ev.top,FIND_BY_UNAME='REV1')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
          ID=widget_info(ev.top,FIND_BY_UNAME='REV2')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
          ID=widget_info(ev.top,FIND_BY_UNAME='REV3')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)
          ID=widget_info(ev.top,FIND_BY_UNAME='HORIZ')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
          ID=widget_info(ev.top,FIND_BY_UNAME='VERT')
          widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          updatefilecounter,ev,list
          endcase
    '...':     begin
          widget_control,ev.id,get_uvalue=uval
          case uval of
          'OutButton':begin
                   path=''
                   path=DIALOG_PICKFILE(path=list.xdipath,title='Select output directory...',/directory)
                   if path eq '' then return
                   list.xdipath=path
                   ID=widget_info(ev.top,FIND_BY_UNAME='OutDirText')
                   widget_control,ID,set_value=list.xdipath
                   endcase
          'MaskButton':begin
                   path=''
                   path=DIALOG_PICKFILE(path=list.mskpath,file=list.mskfile,filter='*.msk',title='Select .msk File...')
                   if path eq '' then return
                   error=CutPath(path,path=mskpath,file=file,ext=ext)
                   list.mskpath=mskpath
                   list.mskfile=file+ext
                   ID=widget_info(ev.top,FIND_BY_UNAME='MaskText')
                   widget_control,ID,set_value=list.mskpath+list.mskfile
                   list.xdipath=list.mskpath
                   list.xdifile=file
                   ID=widget_info(ev.top,FIND_BY_UNAME='OutFileText')
                   widget_control,ID,set_value=list.xdifile
                   ID=widget_info(ev.top,FIND_BY_UNAME='OutDirText')
                   widget_control,ID,set_value=list.xdipath
                   endcase
          'PathButton':begin
                      path=''
                   path=DIALOG_PICKFILE(path=list.path,file='',filter='*.*',title='Select path or .lof file...')
                   if path eq '' then return
                   error=CutPath(path,path=pathtemp,ext=ext)
                   list.path=pathtemp
                   ID1=widget_info(ev.top,FIND_BY_UNAME='PathText')
                   ID2=widget_info(ev.top,FIND_BY_UNAME='FilterText')
                   PTR_FREE, list.filter
                   if ext eq '.lof' then begin
                        widget_control,ID1,set_value=path
                        error=ReadLOF(path,lof,list.prid)
                        if error eq 0 then begin
                            printw,list.prid,path+' Reading Error!!!'
                            widget_control,ID2,sensitive=1
                            widget_control,ID2,get_value=val
                            list.filter=PTR_NEW(val[0])
                            list.filtermethod=1
                        endif else begin
                            list.filter=PTR_NEW(lof)
                            widget_control,ID2,sensitive=0
                            list.filtermethod=2
                        endelse
                   endif else begin
                        widget_control,ID1,set_value=list.path
                        widget_control,ID2,sensitive=1
                        widget_control,ID2,get_value=val
                        list.filter=PTR_NEW(val[0])
                        
                        if ptr_valid(list.PTR_files) then begin
                            *list.PTR_files=list.path+FILE_BASENAME(*list.PTR_files)
                            ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
                            widget_control,ID,set_value=*list.PTR_files
                        endif
                   endelse

                   endcase
          endcase
          endcase
    endcase
    endcase
6:  begin;list
    MainUpdateDo,list,ev
    endcase
3:  begin
    widget_control,ev.id,get_uvalue=uval
    widget_control,ev.id,get_value=val
    val=val[0]
    case uval of
    'NText':   begin
          list.N=long(val)+1>1
          case list.SCANDIM of
            1:    list.nrtif=list.N
            2:    list.nrtif=list.CR[0]*list.CR[1]
          endcase
          resetmag_xrdbp,list
          widget_control,ev.top,set_uvalue=list
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          updatefilecounter,ev,list
          endcase
    'ColText': begin
          list.CR[0]=long(val)+1>1
          case list.SCANDIM of
            1:    list.nrtif=list.N
            2:    list.nrtif=list.CR[0]*list.CR[1]
          endcase
          resetmag_xrdbp,list
          widget_control,ev.top,set_uvalue=list
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          updatefilecounter,ev,list
          endcase
    'RowText': begin
          list.CR[1]=long(val)+1>1
          case list.SCANDIM of
            1:    list.nrtif=list.N
            2:    list.nrtif=list.CR[0]*list.CR[1]
          endcase
          resetmag_xrdbp,list
          widget_control,ev.top,set_uvalue=list
          result=PrepOut(list,ev)
          RefreshDisplayBatch,list
          updatefilecounter,ev,list
          endcase
    'CT':   begin
          list.ct=float(val)
          endcase
    'PathText':    begin
          path=val[0]
          error=CutPath(path,path=pathtemp,ext=ext)
          list.path=pathtemp
          ID1=widget_info(ev.top,FIND_BY_UNAME='PathText')
          ID2=widget_info(ev.top,FIND_BY_UNAME='FilterText')
          PTR_FREE, list.filter
          if ext eq '.lof' then begin
                    widget_control,ID1,set_value=path
                    error=ReadLOF(path,lof,list.prid)
                    if error eq 0 then begin
                        printw,list.prid,path+' Reading Error!!!'
                        widget_control,ID2,sensitive=1
                        widget_control,ID2,get_value=val
                        list.filter=PTR_NEW(val[0])
                        list.filtermethod=1
                    endif else begin
                        list.filter=PTR_NEW(lof)
                        widget_control,ID2,sensitive=0
                        list.filtermethod=2
                    endelse
          endif else begin
                    widget_control,ID1,set_value=list.path
                    widget_control,ID2,sensitive=1
                    widget_control,ID2,get_value=val
                    list.filter=PTR_NEW(val[0])
                   
                    if ptr_valid(list.PTR_files) then begin
                        *list.PTR_files=list.path+FILE_BASENAME(*list.PTR_files)
                        ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
                        widget_control,ID,set_value=*list.PTR_files
                    endif
          endelse

          endcase
    'FilterText':begin
          pos=strpos(val[0],',')
          if pos ne -1 then list.filtermethod=3 else list.filtermethod=1
          (*list.filter)=val[0]
          endcase
    'MaskText':    begin
          if val eq '' then return
          error=CutPath(val,path=mskpath,file=file,ext=ext)
          list.mskpath=mskpath
          list.mskfile=file+ext
          
          list.xdipath=list.mskpath
          list.xdifile=file
          ID=widget_info(ev.top,FIND_BY_UNAME='OutFileText')
          widget_control,ID,set_value=list.xdifile
          ID=widget_info(ev.top,FIND_BY_UNAME='OutDirText')
          widget_control,ID,set_value=list.xdipath
          endcase
    'OutFileText':    list.xdifile=val
    'OutDirText': list.xdipath=val
    'C1Text0': begin
          list.coord1[0]=float(val)
          UpdateCText,list,ev,uval
          RefreshDisplayBatch,list
          endcase
    'C1Text1': begin
          list.coord1[1]=float(val)
          UpdateCText,list,ev,uval
          RefreshDisplayBatch,list
          endcase
    'CText0':  begin
          list.coord[0]=float(val)
          UpdateCText,list,ev,uval
          RefreshDisplayBatch,list
          endcase
    'CText1':  begin
          list.coord[1]=float(val)
          UpdateCText,list,ev,uval
          RefreshDisplayBatch,list
          endcase
    'CText2':  begin
          list.coord[2]=float(val)
          UpdateCText,list,ev,uval
          RefreshDisplayBatch,list
          endcase
    'CText3':  begin
          list.coord[3]=float(val)
          UpdateCText,list,ev,uval
          RefreshDisplayBatch,list
          endcase
    endcase
    endcase
2:  begin; Slider
    widget_control,ev.id,get_value=val
    list.sec=val ; this will be set as timeout when the current running timer times out
    ; Make the current timer time out by hitting it
    WIDGET_CONTROL, ev.top, TIMER=0
    endcase
5:  begin; Label (Process Timer Event)

    ; ----Stop timer?----
    if list.timer eq 0 then return

    ; ----Return when no files present or waiting for file----
    if ptr_valid(list.PTR_files) then bsearch=list.nfiles ge n_elements(*list.PTR_files) $
    else bsearch=1b
    if bsearch then begin
       ; Stop process timer and hit the files timer
       WIDGET_CONTROL, ev.top, TIMER=0
       return
    endif
    T=systime(1)

    ; ----Acquisition time----
    ACQ=T-list.ACQSTAMP
    printw,list.prid,'Detect file: '+string(ACQ)+' sec'

    ; ----Process----
    error=ProcessFile(list,ev)
    if error eq 0 then return ;File couldn't be loaded

    list.ACQSTAMP=systime(1)
    PROC=(list.ACQSTAMP-T)+ACQ
    printw,list.prid,'Total file time:'+string(PROC)+' sec'

    ; ----Time left----
    printw,list.prid,'Number of patterns left:'+string(list.nrtif-list.nFiles)
    ETIME=(list.nrtif-list.nFiles)*PROC
    ConvertTime,ETIME,timeaf
    printw,list.prid,'Estimated time left:'+string(ETIME)+timeaf
    printw,list.prid,'-----------------------'

    ; ----Update file counter----
    updatefilecounter,ev,list

    ; ----End?----
    if (list.nrtif eq list.nFiles) then result=EndBatchSequence(list,ev)
    
    ; ----Hit process timer----
    WIDGET_CONTROL, ev.id, TIMER=0

    ; ----Save Session?----
    if list.SR ne 0 then begin
        if list.SR eq 2 then begin
            case list.SCANDIM of
            1:    n=list.N
            2:    n=list.CR[0]>list.CR[1]
            endcase
            bSR=(list.nFiles mod n) eq 0
        endif else bSR=1b
        if bSR then begin
               widget_control,ev.top,set_uvalue=list
               SSession_XRD_BP,ev,/auto
               return
           endif
    endif else begin
        if (list.nFiles eq list.nrtif) and (list.FASTMODE ne 0) then begin
            widget_control,ev.top,set_uvalue=list
               SSession_XRD_BP,ev
               return
        endif
    endelse
    

    endcase
8:    begin
    widget_control,ev.id,get_uvalue=uval
    case uval of
    'Unit1': list.unit1=ev.index
    'Unit2': list.unit[0]=ev.index
    'Unit3': list.unit[1]=ev.index
    endcase
    RefreshDisplayBatch,list
    endcase
else:return
endcase
widget_control,ev.top,set_uvalue=list
end;pro XRD_BP_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRD_BP,ev=ev
if keyword_set(ev) then begin
    widget_control,ev.top,get_uvalue=list
    path=list.path
    mskpath=list.mskpath
    mskfile=list.mskfile
    tmp=CutPath(mskfile,file=xdifile)
    Platform=list.Platform
    top=ev.top
    sclmax=list.sclmax
    Format=list.Format
    filter='*'+Format
    TiePath=list.TiePath
    TieFile=list.TieFile
    
    AzIntInfo=list.AzIntInfo
    readinfo=list.readinfo
endif else begin
    top=0L
    ev={WIDGET_BUTTON,ID:-1L,TOP:top,HANDLER:top,SELECT:0L}
    path=''
    mskpath=''
    mskfile=''
    xdifile=''
    TiePath=''
    TieFile=''
    ControlDevice,Platform,Visual

    result=ReadINI('$Scale:')
    clipslider=[0.,(result.(0))[2],(result.(1))[1],(result.(0))[0]]
    sclmax=2.^(clipslider[3]*clipslider[2])-1

    Format='.tif'
    filter='*'+Format
    
    readinfo=ReadINI('$ReadCCDInfo:')
    
    result=ReadINI('$AzInt:')
    AzIntInfo={oD:result.(0)[0],soD:result.(0)[1],Fill:result.(0)[2],type:result.(0)[3],$
            bool:result.(0)[4],SaveRingPercentage:result.(0)[5],Batch:result.(0)[6],$
            tiff:result.(0)[7],origformat:result.(0)[8],$
            MinPix:result.(1)[0],BWP:result.(1)[1],BWT:result.(1)[2],BWD:result.(1)[3],$
            BWQ:result.(1)[4],AzWidth:result.(1)[5],$
            outbin:result.(2)[0],median:result.(3),$
            stD:result.(4),nsectors:result.(5),$
            aoD:result.(6)[0],loD:result.(6)[1],percentile:result.(7),$
            updatenorm:result.(8)[0],calcerror:result.(8)[1]}
endelse

CorSeq=PTR_NEW({zinger:1B,satur:2B,ff:3B,back:4B,spatial:5B,maskoff:0B})
SeqTags=STRLOWCASE(tag_names(*CorSeq))

;result=ReadINI('$CrossUpdate:')
;CrossUpdate=result.(0)
CrossUpdate=0b

; ----Setup Environment list----
list={listtype:40,$

    ; obsolete file handling (too slow)
    PTR_list:PTR_NEW(),$    ; Linked list with processed file info
    llist:{str:'',next:PTR_NEW()},$ ; Format of linked list
    ; new file handling
    PTR_files:PTR_NEW(),$    ; files to process
    
    ptrModel:CreateModel(),$    ; analyzing model
    
    Subfiles:{n:0l,nleft:0l,fileind:0l,ptrtype:0,fileoffset:0l,filemult:0l,xtype:-1,ptr:PTR_NEW()},$
                            ; One file contains the info of several files:
                            ; n: number of subfiles
                            ; nleft: number of subfiles unprocessed
                            ; fileind: file index in datablock
                            ; ptrtype: 1(1D), 2(2D)
                            ; fileoffset: fileind offset
                            ; filemult: fileind multiplier
                            ; xtype: X type for 1D
                            ; ptr: point to read subfiles
    
    tiffb:PTR_NEW(),$       ; Dark image when necessary
    ff:PTR_NEW(),$            ; inverse + scaled flat field image
    tiffs:PTR_NEW(),$        ; Size info on current processed image
    maskoff:PTR_NEW(),$        ; Mask image
    readinfo:readinfo,$     ; Normalize EDF on reading
    
    CorMethod:lonarr(6),$      ; zinger removal: 0=no, 1=yes
                            ; saturation removal: 0=no, 1=yes
                               ; background correction: 0=no, 1=stripping, 2=dark image
                               ; spatial distortion corr.: 0=no, 1=yes
                               ; flat field correction: 0=no, 1=yes
                               ; mask off: 0=no, 1=yes
    CorSeq:CorSeq,$         ; Positions of the fields in the structure gives sequence,
                               ; the value in the field gives the position in CorMethode/Done
    SeqTags:SeqTags,$        ; Tag names of CorSeq
    
    ; -Flags
    update:CrossUpdate,$       ; Update Main applet
    SR:0B,$                   ; Save processed data so far
    PLOT:1B,$                  ; Show image or linescan-plot when processing
    SCANDIM:1,$                ; 1=> Scan 2=> Map
    SMALL2D:0B,$              ; set if small ROI's (i.e. no big circles)
    SMALL1D:0B,$              ; set if small ROI's (i.e. no big circles)
    NORMROIs:0B,$            ; normalize ROI's
    DEBUG:0B,$              ; see fitted data
    timer:0B,$              ; Timer running or not
    pause:0B,$              ; To make the distinction between timer=0 'Not started yet or aborded'
                             ; and timer=0 'Paused'
    TRAN:0B,$                  ; 1=> First samv and then samh
    REV1:0B,$                  ; 1=> RL and not LR
    REV2:0B,$                  ; 1=> BU and not UB
    REV3:0B,$                ; 1=> stages don't jump after finishing a row or column
    sortmethod:1B,$            ; For sorting incoming files
    Hsign:1,$                  ; 1 => LR positive movement, -1 => RL positive movement
    Vsign:1,$                  ; 1 => UB positive movement, -1 => BU positive movement
    AUTO:0B,$                  ; Don't prompt when .xdi overwritten
    MSESSION:1B,$            ; Save multiple sessions
    FASTMODE:0B,$            ; No full processing mode
    FORCE:1B,$                ; Use background in explorative mode? (if yes, override with ortpol degree 1)
    ; -Indices
    prid:0L,$                  ; Output for status messages
    outdraw:0L,$             ; display window
    outdrawid:0L,$           ; Index of display window
    outdrawsum:0L,$         ; display window (plot sumspectrum)
    outdrawidsum:0L,$       ; Index of display window (plot sumspectrum)
    pixindexsum:0L,$        ; Index of pixmap of outdrawsum
    top:top,$                  ; ID of Main applet
    outid1:-1L,$             ; Window for plotting things when debug=1 "Warping(1D) or Peaksearch(2D)"
    outid2:-1L,$             ; Window for plotting things when debug=1 "1D Fit result"
    outid3:-1L,$             ; Window for plotting things when debug=1 "2D Fit result"
    xs:0L, $                ; X static corner of the zoom box.
    ys:0L, $                ; Y static corner of the zoom box.
    xd:0L, $                ; X dynamic corner of the zoom box.
    yd:0L,$                 ; Y dynamic corner of the zoom box.
    sysvar:PointSysVar(),$    ; System variable for coord. convertion
    sysvar2:PointSysVar(),$    ; System variable for coord. convertion
    ; -Initialisation data
    filter:PTR_NEW(filter),$; Filter for processed files
    filtermethod:1,$           ; 1: '*.*'
                             ; 2: list of files
                             ; 3: '*.*,26,10' (i.e. skip first 26, read only 10 files)
    nocheck:1b,$            ; Don't check files before adding to batch que
    N:1L,$                    ; Number of points in linescan
    CR:[1L,1L],$             ; Dimension of mapping
                             ; E.g. samh ... ... 20 samv ... ... 31
                             ; or samv ... ... 20 samh ... ... 31
                             ; CR=[21,32]
    coord:fltarr(2,2),$        ; Scan coordinates in mm or degrees for map
                             ; E.g. samh -0.38 -0.43 ... samv -0.665 -0.615 ...
                             ; or samv -0.38 -0.43 ... samh -0.665 -0.615 ...
                             ; -0.380 -0.430  (first motor)
                             ; -0.665 -0.615  (second motor)
    coord1:fltarr(2),$       ; Scan coordinates in mm or degrees for linescan
    unit:intarr(2),$        ; map units
    unit1:0,$                ; scan units
    ct:10,$                   ; Collection time in seconds
    nFiles:0L,$              ; Number of files processed
    nrtif:0L,$                ; Number of files to process
    nmag:1,$                  ; Magnification of display for mapping
    nmult:0.,$                ; Scaling of display for mapping (% of maximum)
    xsx:0,$                    ; maps in x-direction for displaying
    ysy:0,$                    ; maps in y-direction for displaying
    xsm:600,$                  ; Max xsize of drawwidget
    ysm:284,$                  ; Max ysize of drawwidget
    ACQSTAMP:0D,$              ; For estimating acquisition time
    TotTime:dblarr(2),$        ; Total time | paused/resume time lost
    ; -Begin data from .msk file
    mskpath:mskpath,$         ; Path of .msk file
    mskfile:mskfile,$         ; Name of .msk file
    lambda:0.,$                ; X-ray lambda(Angstroms)
    scrx:0.,$                  ; Pixel x-size (m, input in um)
    scry:0.,$                  ; Pixel y-size (m, input in um)
    dist:0.,$                  ; Distance sample-detector (m, input in cm)
    wb:3,$                   ; Strip-width for background calc
    niter:20,$              ; Number of iterations for background calc
    thres:[0.,0.],$            ; Cutoff value
    bcko:0,$                ; Order of dark-image smoothing
    bckw:5,$                ; Width of dark-image smoothing
    bckm:1.,$                ; Background multiplication
    v:3,$                      ; side windows width for Top-hat (peaksearch)
    w:4,$                     ; center window width for Top-hat (peaksearch)
    r:3,$                     ; criterium for peaksearch
    fit:0,$                  ; fit with: real(0),suggestive(1),no(2) constraints
    masknr:0,$                ; number of areas in mask
    mask:PTR_NEW(),$              ; first column: maskpart type (1:square, 2:arc, 3: circle)
                             ; second column: index of first part of a linked group of parts
                             ; third column: 1=set, 0=not set
                             ; rest: parameters of maskpart
                             ;     square: x1,x2,y1,y2 => LLcorner(1) and URcorner(2) coord
                             ;     arc:r1,theta1,r2,theta2
                             ;     circle:r1,r2,0,0
    masknr1:0,$                ; number of areas in 1D mask
    mask1:PTR_NEW(),$        ; first column: maskpart type (1:ROI)
                            ; second column: index of first part of a linked group of parts
                               ; third column: 1=set, 0=not set
                               ; rest: parameters of maskpart
                               ;     left,right coord in two-theta
    center:fltarr(2)-1,$       ; coordinates (2D pixel index= 2D coord in tiffmatrix) of center
    SpaDisStruc:PTR_NEW(),$    ; Spatial distortion correction information
    WarpStruc:PTR_NEW(),$    ; Azimuthal integration information
    TiePath:TiePath,$        ; spline file
    TieFile:TieFile,$        ; spline file
    nrTieRing:0,$           ; Number of Tiepoint rings
    nrTie:PTR_NEW([0]),$    ; Number of Tie points in each ring
    TieRingD:PTR_NEW(),$       ; D-spacing for each Tie-ring
    TieI:PTR_NEW(),$           ; Tie point coordinates of inputimage
    TieO:PTR_NEW(),$        ; Tie point coordinates of outputimage
    TieGrid:1B,$            ; Tie points on grid or circles
    bFastSpline:1b,$        ; Perform spline spatial distortion in fast mode
    vFlipSpline:1b,$        ; Flip image before spatial distortion (when using fit2d spline file)
    a:0.,$                     ; a Tilt angle: RotXc(a)
    b:0.,$                     ; b Rotation angle: RotZt(b)
    phipixphi0:0.,$         ; azimuthal shift due to 2-angle calibration
    IntCor2Dinfo:DefaultIntCor2Dinfo(),$ ; 2D intensity correction parameters
    xdipath:mskpath,$        ; output directory
    xdifile:xdifile,$        ; output filename
    darkpath:mskpath,$      ; Path of dark image
    darkfile:'',$            ; Dark image file
    ffpath:mskpath,$          ; Path of flat field image
    fffile:'',$                ; Flat field image file
    maskoffpath:mskpath,$   ; Path of mask image
    maskofffile:'',$        ; Mask image file
    maskoffthres:'> 0',$    ; Mask image threshold
    maskoffdynthres:'> 0',$    ; Mask image dynamic threshold (on the image itself)
    maskoffdyn:0b,$            ; do dynamic mask off
    vin:2,$                  ; side windows width for 1D Top-hat (peaksearch)
    win:3,$                  ; center window width for 1D Top-hat (peaksearch)
    r1:0.5,$                   ; criterium for 1D peaksearch
    filorder:0B,$             ; peaksearch before (0) or after subtraction of background
    fitmethod:5,$            ; 1D fit method
    backmod:-1,$             ; 1D background method
    AzIntParam:fltarr(5),$    ; Parameters for Azimuthal integration: d1, d2, noutbins, azbegin, azend
    backranindtt:fltarr(2),$; ROI in two-theta
    AzIntInfo:AzIntInfo,$    ; Azimuthal integration info
    ; -Misc
    fastmodesuperimposed:0b,$
    xrange:dblarr(2),$
    yrange:dblarr(2),$
    plottitle:'',$            ; Title of 1D plot
    groupi:0L,$                ; Group to plot
    profileptr:ptr_new(/alloc),$ ; 1D profiles
    path:path,$              ; Path of files to process
    Platform:Platform,$     ; Platform currently on
    ev:ev,$                   ; Button in main window that generated the event
    SaturMax:0.,$            ; Maximum Pattern value
    SaturAspect:2.,$        ; Saturation aspect ratio
    SaturNMax:0l,$            ; Saturation maximal area
    ZingerWidth:3,$            ; Smoothing width in zinger removal
    ZingerThreshold:1.2,$    ; Ratio threshold in zinger removal
    sclmax:sclmax,$            ; For display of 2D fitting (peaksearch)
    sec:-1,$                   ; Timer event speed: -1 not searching for files, 0 as quick as possible,...
    sortseparator:'_',$     ; For sorting incoming files
    chari:0,$                  ; Process character index
    char:['  \  ','  |  ','  /  '],$; Process chars
    title:'Batch Processing: '}; Applet base-title

; ----Main base----
TopBase=widget_base(MBAR=bar,/row,title=list.title+(list.char)[list.chari],xoffset=0,yoffset=0)

FilesBase=widget_base(TopBase,/column)
    ModeBase=widget_base(FilesBase,/row,/exclusive,uname='modebase')
        buttonfull=widget_button(ModeBase,value='Full processing mode',uname='fullmode')
        buttonfast=widget_button(ModeBase,value='Explorative mode',uname='fastmode')
        buttongrid=widget_button(ModeBase,value='Grid mode',uname='gridmode')

    dlxs=50
    dlys=20
    droplist=widget_list(FilesBase,uname='FilesList',ysize=dlys,xsize=dlxs)
    text=widget_text(FilesBase,value='0 files',uname='FilesText')
    StatusText=widget_text(FilesBase,uname='Status',uvalue=[0,500,5],ysize=dlys,xsize=dlxs,/scroll);,/wrap
    list.prid=StatusText

WrapBase=widget_base(TopBase,/column)
IABase=widget_base(Wrapbase,/row)
WrapInitBase=widget_base(IABase,/row,uname='Initbase')

WrapInitBase=widget_tab(WrapInitBase)
InitBase=widget_base(WrapInitBase,/column,title='Files')
    xs=85

    base=widget_base(InitBase,/row)
    label=widget_label(base,value='*.*(,nskip,nmax)',xsize=xs)
    text=widget_text(base,/editable,value=*(list.filter),uvalue='FilterText',uname='FilterText')

    base=widget_base(InitBase,/row)
    label=widget_label(base,value='Dir or list of files:',xsize=xs)
    text=widget_text(base,/editable,value=path,uvalue='PathText',uname='PathText')
    button=widget_button(base,value='...',uvalue='PathButton')

    button=widget_button(InitBase,value='File sorting')

    base=widget_base(InitBase,/row,uname='Maskbase')
    label=widget_label(base,value='Mask File:',xsize=xs,uname='TimerLabel')
    text=widget_text(base,/editable,value=list.mskpath+list.mskfile,uvalue='MaskText',uname='MaskText')
    button=widget_button(base,value='...',uvalue='MaskButton')
    
    base=widget_base(InitBase,/row)
    label=widget_label(base,value='Output dir:',xsize=xs)
    text=widget_text(base,/editable,value=list.xdipath,uvalue='OutDirText',uname='OutDirText')
    button=widget_button(base,value='...',uvalue='OutButton')

    base=widget_base(InitBase,/row)
    label=widget_label(base,value='Output filename:',xsize=xs)
    text=widget_text(base,/editable,value=list.xdifile,uvalue='OutFileText',uname='OutFileText')

    base=widget_base(InitBase,/column,/nonexclusive)
       Outbutton=widget_button(base,value='Autosave Session/Point',uname='Outbutton')
       Outbutton2=widget_button(base,value='Autosave Session/Line',uname='Outbutton2')
;       MSESSIONbutton=widget_button(base,value='Multiple sessions',uname='MSESSIONbutton')
       Autobutton=widget_button(base,value='Overwrite existing outputfiles',uname='Autobutton')
       checkbutton=widget_button(base,value="Don't check files",uname='checkbutton')

InitBase=widget_base(WrapInitBase,/column,title='Scan Dimensions')
    base=widget_base(InitBase,/row,/exclusive)
    Scanbutton=widget_button(base,value='Line',uname='Scanbutton')
    Mapbutton=widget_button(base,value='Map',uname='Mapbutton')

    xss=5
    xsss=50
    label=widget_label(Initbase,value='Direction     start     stop     steps     unit')
    base=widget_base(InitBase,/row,uname='NBase',sensitive=list.SCANDIM eq 1)
    label=widget_label(base,value='MotH',xsize=xsss,uname='NLabel')
    text=widget_text(base,/editable,value=stringr(list.coord1[0]),$
       uvalue='C1Text0',uname='C1Text0',xsize=xss)
    text=widget_text(base,/editable,value=stringr(list.coord1[1]),$
       uvalue='C1Text1',uname='C1Text1',xsize=xss)
    text=widget_text(base,/editable,value=stringr(list.N-1),$
       uvalue='NText',uname='NText',xsize=xss)
    drop=widget_droplist(base,value=['mm','deg'],uvalue='Unit1',uname='Unit1')

    base=widget_base(InitBase,/row,uname='ColBase',sensitive=list.SCANDIM eq 2)
    label=widget_label(base,value='MotH',xsize=xsss,uname='ColLabel')
    text=widget_text(base,/editable,value=stringr(list.coord[0]),$
       uvalue='CText0',uname='CText0',xsize=xss)
    text=widget_text(base,/editable,value=stringr(list.coord[1]),$
       uvalue='CText1',uname='CText1',xsize=xss)
    text=widget_text(base,/editable,value=stringr(list.CR[0]-1),$
       uvalue='ColText',uname='ColText',xsize=xss)
    drop=widget_droplist(base,value=['mm','deg'],uvalue='Unit2',uname='Unit2')

    base=widget_base(InitBase,/row,uname='RowBase',sensitive=list.SCANDIM eq 2)
    label=widget_label(base,value='MotV',xsize=xsss,uname='RowLabel')
    text=widget_text(base,/editable,value=stringr(list.coord[2]),$
       uvalue='CText2',uname='CText2',xsize=xss)
    text=widget_text(base,/editable,value=stringr(list.coord[3]),$
       uvalue='CText3',uname='CText3',xsize=xss)
    text=widget_text(base,/editable,value=stringr(list.CR[1]-1),$
       uvalue='RowText',uname='RowText',xsize=xss)
    drop=widget_droplist(base,value=['mm','deg'],uvalue='Unit3',uname='Unit3')

    base=widget_base(InitBase,/row)
    label=widget_label(base,value='Exposure time (sec):')
    text=widget_text(base,/editable,value=stringr(list.ct),$
       uvalue='CT',uname='CT',xsize=xss)

    base=widget_base(InitBase,/row,/exclusive)
    Movehbutton=widget_button(base,value='First horizontal',uname='First horizontal')
    Movevbutton=widget_button(base,value='First vertical',uname='First vertical')

    base=widget_base(InitBase,/row,/exclusive,uname='REV1',sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1))
    LRbutton=widget_button(base,value='Left-Right',uname='Left-Right')
    RLbutton=widget_button(base,value='Right-Left',uname='Right-Left')

    base=widget_base(InitBase,/row,/exclusive,uname='REV2',sensitive=(list.SCANDIM eq 2)or(list.TRAN))
    UBbutton=widget_button(base,value='Top-Bottom',uname='Top-Bottom')
    BUbutton=widget_button(base,value='Bottom-Top',uname='Bottom-Top')

    base=widget_base(InitBase,/row,/exclusive,uname='REV3',sensitive=(list.SCANDIM eq 2))
    Jbutton=widget_button(base,value='Jump',uname='Jump')
    Cbutton=widget_button(base,value='Continue',uname='Continue')

    base=widget_base(InitBase,/row)
    baset=widget_base(base,/row,uname='HORIZ',sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1))
    label=widget_label(baset,value='Hori+:')
    button=widget_button(baset,value='L-->R',uname='L--R')
    baset=widget_base(base,/row,uname='VERT',sensitive=(list.SCANDIM eq 2)or(list.TRAN))
    label=widget_label(baset,value='Vert+:')
    button=widget_button(baset,value='T-->B',uname='T--B')

InitBase=widget_base(WrapInitBase,/column,title='Process options')
    base=widget_base(InitBase,/column,/nonexclusive)
    Smallbutton2D=widget_button(base,value='2D background on ROI`s',uname='Smallbutton2D')
    Smallbutton1D=widget_button(base,value='1D background on ROI`s',uname='Smallbutton1D')
    NORMROIsbutton=widget_button(base,value='Normalized ROI`s',uname='NORMROIsbutton')
    Plotbutton=widget_button(base,value='View progress',uname='Plotbutton')
    Debugbutton=widget_button(base,value='Monitor processing',uname='Debugbutton')
    Forcebutton=widget_button(base,value='Explorative background',uname='Forcebutton')
    Supbutton=widget_button(base,value='Explorative superimposed profile',uname='Supbutton')

DrawBase=widget_base(WrapBase,/column,uname='DrawBase')

WrapActiveBase=widget_base(IABase,/column,uname='DrawBase2')
WrapActiveBase=widget_base(WrapActiveBase,/row)
ActiveBase=widget_base(WrapActiveBase,/column)
    button=widget_button(ActiveBase,value='+',uname='+')
    button=widget_button(ActiveBase,value='-',uname='-')
    button=widget_button(ActiveBase,value='<<',uname='<<')
ActiveBase=widget_base(WrapActiveBase,/column)
    button=widget_button(ActiveBase,value='*',uname='*')
    button=widget_button(ActiveBase,value='/',uname='/')
    button=widget_button(ActiveBase,value='>>',uname='>>')
ActiveBase=widget_base(WrapActiveBase,/column)
    button=widget_button(ActiveBase,value='Go',uname='Go')
    button=widget_button(ActiveBase,value='Abort',uname='Abort',sensitive=0)
    button=widget_button(ActiveBase,value='Pause',uname='Pause',sensitive=0)
ActiveBase=widget_base(WrapActiveBase,/column)
    base=widget_base(ActiveBase,/row,/nonexclusive)
       Updatebutton=widget_button(base,value='Update Main Window',uname='Update Main Window')
    label=widget_label(ActiveBase,value='File Search Speed (sec.):')
    slider=widget_slider(ActiveBase,value=list.sec,minimum=-1,maximum=300,uname='SecSlider')
    button=widget_button(ActiveBase,value='Search now')

menu1=WIDGET_BUTTON(bar, VALUE='File', /MENU)
       button = WIDGET_BUTTON(menu1, VALUE='Save Session')
       button = widget_button(menu1, value='Restore Session')
       button = widget_button(menu1, value='Save current result')
       button = widget_button(menu1, value='Save current result (no tiff)')
       button=widget_button(menu1,value='Exit',uname='Exit',/separator)

WIDGET_CONTROL, TopBase, /REALIZE,set_uvalue=list
widget_control,Updatebutton,set_button=list.update
widget_control,checkbutton,set_button=list.nocheck
widget_control,Autobutton,set_button=list.AUTO
widget_control,Outbutton,set_button=list.SR eq 1
widget_control,Outbutton2,set_button=list.SR eq 2
widget_control,Plotbutton,set_button=list.PLOT
widget_control,Debugbutton,set_button=list.DEBUG
widget_control,Scanbutton,set_button=list.SCANDIM eq 1
widget_control,Mapbutton,set_button=list.SCANDIM eq 2
widget_control,Smallbutton2D,set_button=list.SMALL2D
widget_control,Smallbutton1D,set_button=list.SMALL1D
widget_control,NORMROIsbutton,set_button=list.NORMROIs
widget_control,Movehbutton,set_button=1
widget_control,LRbutton,set_button=1
widget_control,UBbutton,set_button=1
widget_control,Jbutton,set_button=1
;widget_control,MSESSIONbutton,set_button=list.MSESSION
widget_control,Forcebutton,set_button=list.FORCE
widget_control,Supbutton,set_button=list.fastmodesuperimposed
case list.fastmode of
0: widget_control,buttonfull,set_button=1
1: widget_control,buttonfast,set_button=1
2: widget_control,buttongrid,set_button=1
endcase

;widget_control,/DELAY_DESTROY
Xmanager,'XRD_BP',TopBase,/no_block,cleanup='CleanUp_XRD_BP'
end;pro XRD_BP
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
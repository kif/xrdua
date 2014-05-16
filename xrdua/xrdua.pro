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

function PointInsidePolygon,x,y,nodes
nnodes=n_elements(nodes)/2
if nnodes le 1 then return,0b

object = Obj_New('IDLanROI', nodes[0,1:*], nodes[1,1:*])
b=object->ContainsPoints(x, y)
Obj_Destroy, object
return,b
   
end;function PointInsidePolygon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PixelsInPolygon,sx,sy,px,py
return,POLYFILLV( round(px), round(py), sx, sy )

; Possible pixels
bx=floor(min(px,max=ex))>0
by=floor(min(py,max=ey))>0
ex=ceil(ex)<(sx-1)
ey=ceil(ey)<(sy-1)
nx=ex-bx+1
ny=ey-by+1
x=rebin(bx+lindgen(nx),nx,ny,/sample)
y=rebin(by+lindgen(1,ny),nx,ny,/sample)
end;function PixelsInPolygon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DrawPolygon,list

nnodes=n_elements(*list.maskpolygon)/2
if nnodes le 1 then return
x=(*list.maskpolygon)[0,1:*]
y=(*list.maskpolygon)[1,1:*]
nnodes--

case nnodes of
1:    begin
    xx=RTStat(x,list.sfh)
    yy=RTStat(y,list.sfv)
    WSet, list.drawindex
    Plots,[xx-list.CW1,xx+list.CW1],[yy,yy], COLOR=150,/device
    Plots,[xx,xx],[yy-list.CW1,yy+list.CW1], COLOR=150,/device
    if widget_info(list.drawdyn,/VALID_ID) then begin
        xx=RTDyn(x,list.sf,list.hrange[0])
        yy=RTDyn(y,list.sf,list.vrange[0])
        WSet, list.drawdynindex
        Plots,[xx-list.CW2,xx+list.CW2],[yy,yy], COLOR=150,/device
        Plots,[xx,xx],[yy-list.CW2,yy+list.CW2], COLOR=150,/device
    endif
    endcase
2:    begin
    WSet, list.drawindex
    plots, RTStat(x,list.sfh),RTStat(y,list.sfv), COLOR=150,/device
    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        plots, RTDyn(x,list.sf,list.hrange[0]),RTDyn(y,list.sf,list.vrange[0]), COLOR=150,/device
    endif
    endcase
else:begin
    WSet, list.drawindex
    POLYFILL, RTStat(x,list.sfh),RTStat(y,list.sfv), COLOR=150,/device,/fill
    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        POLYFILL, RTDyn(x,list.sf,list.hrange[0]),RTDyn(y,list.sf,list.vrange[0]), COLOR=150,/device,/fill
    endif
    endcase
endcase

end;pro DrawPolygon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddPolygon_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

bdown=ev.type eq 0
bkey=ev.type eq 5 or ev.type eq 6

if bdown then begin
    bleft=ev.press eq 1
    bright=ev.press eq 4
    if ~bleft and ~bright then return
    
    ; Whipe out previously drawn polygon
    widget_control,ev.top,get_uvalue=list
    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    endif

    if bleft then begin ; Add node
        ; calcultate image position of cursor
        x = ev.x
        y = ev.y
        if ev.id eq list.draw then begin
            x=0>StatTR(x,list.sfh)<((*list.tiffs)[1]-1)
            y=0>StatTR(y,list.sfv)<((*list.tiffs)[2]-1)
        endif else begin
            x=0>DynTR(x,list.sf,list.hrange[0])<((*list.tiffs)[1]-1)
            y=0>DynTR(y,list.sf,list.vrange[0])<((*list.tiffs)[2]-1)
        endelse
;        if ~PointInsidePolygon(x,y,*list.maskpolygon) then $
            *list.maskpolygon=[[*list.maskpolygon],[x,y]]
    endif else begin; Remove node
        widget_control,ev.top,get_uvalue=list
        nnodes=n_elements(*list.maskpolygon)/2
        if nnodes le 1 then return
        *list.maskpolygon=(*list.maskpolygon)[*,0:nnodes-2]
    endelse
    
    ; Redraw polygon
    DrawPolygon,list
    
endif else if bkey then begin
    if ~(ev.press and (ev.ch eq 13 or ev.ch eq 27)) then return

    ; Whipe out drawn polygon
    widget_control,ev.top,get_uvalue=list
    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    endif

    case ev.ch of
    13: begin ; ENTER => accept
        if not ptr_valid(list.tiffvalid) then list.tiffvalid=ptr_new(make_array((*list.tiffs)[1:2],value=1b))
        
        (*list.tiffvalid)[PixelsInPolygon((*list.tiffs)[1],(*list.tiffs)[2],(*list.maskpolygon)[0,1:*],(*list.maskpolygon)[1,1:*])]=0
        RefreshDisplay,ev.top,list
        endcase
    27: begin ; ESC => cancel
        endcase
    endcase
    *list.maskpolygon=fltarr(2)
endif

end;pro AddPolygon_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ConvertCCDimages,ev
widget_control,ev.top,get_uvalue=list
files=Select_Files(ev.top,list,outpath=pathred,outfile=file,sortmethod=list.sortmethod,separator=list.sortseparator)
if files[0] eq '' then return

path=DIALOG_PICKFILE(path=pathred,title='Save in directory...',/directory)
if path EQ '' then return

outformat=FileFormat(2)
outformat=outformat[promptnumber(outformat,ev.top,'Select output format:',/choice)]

if dialog_message('Invert intensity?',/question) eq 'Yes' then begin
    nbit=long64(promptnumber(stringr(list.clipslider[1]*list.clipslider[2]),ev.top,'Number of bits:'))>0
    nbit=2^nbit-1
endif else nbit=0

progressBar = Obj_New("PROGRESSBAR",title="Converting CCD images...")
progressBar -> Start

n=n_elements(files)
for i=0l,n-1 do begin

    IF progressBar -> CheckCancel() THEN BEGIN
        progressBar -> Destroy
        printw,list.outid,'Process CANCELLED!'
        RETURN
    ENDIF
    
    data=ReadCCD(files[i],/static,error=error,nfiles=nfiles,info=list.readinfo)
    if error then continue
    if nbit gt 0 then data=nbit-data
    
    tmp=CutPath(files[i],file=file)
    if nfiles eq 1 then begin
        error=writeccd(data,path,file+outformat)
    endif else begin
        format='("_",I0'+stringr(strlen(stringr(nfiles-1)))+')'
        
        j=0
        filej=file+string(j,format=format)+outformat
        error=writeccd(data,path,filej)

        for j=1l,nfiles-1 do begin
            data=ReadCCD(files[i],/static,error=error,xdinr=j,info=list.readinfo)
            if error then continue
            if nbit gt 0 then data=nbit-data
            
            filej=file+string(j,format=format)+outformat
            error=writeccd(data,path,filej)
        endfor
    endelse
    
    progressBar -> Update, (i+1.)/n*100
endfor

progressBar -> Destroy

end;pro ConvertCCDimages
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RestoreCorrectionev,ev
widget_control,ev.top,get_uvalue=list
widget_control,ev.id,get_value=val

RestoreCorrection,list

RefreshDisplay, ev.top,list
mode,ev,val
end;pro RestoreCorrectionev
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WMean,x,w

; cfr. Wincross for variance estimator
n=n_elements(x)
m=total(x)/n
vari=total((x-m)^2.)/(n-1)
wtot=total(w)
wm=total(w*x)/wtot
varwm=vari*total(w*w)/(wtot*wtot)

return,[wm,varwm]
end;function WMean
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function OverlapnProfiles,Profiles,xEx=xEx,yEx=yEx

; Profiles is an array of pointers

; ---- Make x extended ----
n=n_elements(Profiles)
for i=0l,n-1 do begin
    incc=((*Profiles[i])[0,1]-(*Profiles[i])[0,0])
    bc=min((*Profiles[i])[0,*],max=ec)
    if n_elements(inc) ne 0 then begin
        b=b<bc
        e=e>ec
        inc=inc<incc
    endif else begin
        b=bc
        e=ec
        inc=incc
    endelse
endfor
nEx=ceil((e-b)/inc)+1

; ---- Make y extented ----
yEx=fltarr(nEx)
for i=0l,n-1 do begin
    tmp=EqualBin((*Profiles[i])[0,*], (*Profiles[i])[1,*], b, inc, nEx, X2=xEx)
    yEx=yEx>(tmp/max(tmp)) ;scale between 0-1 and superimpose
endfor

return,WMean(xEx,yEx)
end;function OverlapnProfiles
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Overlap2Profiles,x1,y1,x2,y2

x1=reform(x1)
x2=reform(x2)

; ---- Make x extented ----
inc1=x1[1]-x1[0]
inc2=x2[1]-x2[0]
inc=inc1<inc2
n1=n_elements(y1)
n2=n_elements(y2)
b1=min(x1,max=e1)
b2=min(x2,max=e2)
if b1 gt e2 or b2 gt e1 then return,[0,0]
b=b1<b2
e=e1>e2
nEx=ceil((e-b)/inc)+1

; ---- Make y extended ----
yEx1=EqualBin(x1, y1, b, inc, nEx)
yEx2=EqualBin(x2, y2, b, inc, nEx)

; ---- Overlap ----
barr=byte((yEx1*yEx2)<1b)
Overlap1by2=total(barr*yEx1)/total(yEx1)
Overlap2by1=total(barr*yEx2)/total(yEx2)

; ---- out ----
return,[Overlap1by2,Overlap2by1]
end;function Overlap2Profiles
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CreateOverlapTable,top
widget_control,top,get_uvalue=list

if list.masknr eq 0 then return

IndMask=where(((*list.smask) eq 1) and $
    (((*list.mask)[0,*] eq 2) or $
    ((*list.mask)[0,*] eq 3)),nProf)

if nProf le 1 then return

str=PromptNumber(stringr(list.OverLapCrit),top,'Enter overlap threshold (between 0 and 1):')
val=0.
reads,str,val
list.OverLapCrit=0>val<1

Window, winid()

; ----prepare azint----
res=BraggPtoX(list.AzIntInfo.MinPix,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
rmin=[res[1],res[0],list.AzIntInfo.MinPix,res[2]]
dr=[list.AzIntInfo.BWD,list.AzIntInfo.BWT,list.AzIntInfo.BWP,list.AzIntInfo.BWQ]

; ----azint----
Profiles=PTRARR(nProf)
for i=0l,nProf-1 do begin
    ; ---- Integrate ----
    d=AzIntWarp((*list.mask)[*,IndMask[i]],list,list.IntCor2Dinfo,list.AzIntInfo.type,rmin,dr,list.AzIntInfo.outbin,list.AzIntInfo.AzWidth,0,$
                [list.AzIntInfo.Fill,list.AzIntInfo.bool],imagevalid=list.tiffvalid,median=list.AzIntInfo.median,percentile=list.AzIntInfo.percentile,$
                updatenorm=list.AzIntInfo.updatenorm,calcerror=list.AzIntInfo.calcerror)

    if i eq nProf-1 then str='overlap table is beeing generated...' else str=stringr(nProf-i-1)+' ahead...'
    plot,d[0,*],d[1,*],xtitle=ChiXtitle(list.AzIntInfo.type),ytitle='Intensity (a.u.)',psym=-2,$
    xstyle=1,ystyle=1,title='Integrated area'+stringr(IndMask[i])+', '+str
    ; ---- Store ----
    Profiles[i]=PTR_NEW(d)
endfor

; ----Calc overlap table----
tmp=lindgen(nProf)
one=replicate(1b,nProf)
col=(tmp#one)
row=(one#tmp)
CombInd=where(col gt row)
index1=col[CombInd]
index2=row[CombInd]

Table=replicate(1.,nProf,nProf)
Table2=Table
; loop over Table indices:
; e.g. Table -> 4x4
;        index1=1,2,3,2,3,3
;        index2=0,0,0,1,1,2
;        => right-top half of the Table
for i=0l,n_elements(index1)-1 do begin
    tmp=Overlap2Profiles((*Profiles[index1[i]])[0,*],(*Profiles[index1[i]])[1,*],$
    (*Profiles[index2[i]])[0,*],(*Profiles[index2[i]])[1,*])
    Table[index2[i],index1[i]]=tmp[0] ; store 1by2 in bottom-left half
    Table[index1[i],index2[i]]=tmp[1] ; store 2by1 in top-right half
    tmp=total(tmp)/2
    Table2[index2[i],index1[i]]=tmp
    Table2[index1[i],index2[i]]=tmp
endfor

; ---- Group ----
index12b=Table2[CombInd] ge list.OverLapCrit
; e.g. Table -> 4x4
;        index1=1,2,3,2,3,3
;        index2=0,0,0,1,1,2
;      index12b=1,0,0,1,0,0
AvoidArray=[-1]
TableMean=fltarr(3)

; Loop over all areas
for i=0l,nProf-1 do begin

    ; last area
    if i eq nProf-1 then begin
        tmp=where(AvoidArray eq i,ct)
        if ct eq 0 then begin ; not already linked in previous loops
            LinkAreas,list.mask,list.masknr,IndMask[i]
            TableMean=[[TableMean],[IndMask[i],WMean((*Profiles[i])[0,*],(*Profiles[i])[1,*])]]
        endif
    endif else begin

        ; Search for overlap with group i
        ind=where(index2 eq i,ct)
        ind1=where(index12b[ind] eq 1,ct)+ind[0]

        if ct eq 0 then begin ; no overlap=> group on his own
            tmp=where(AvoidArray eq i,ct)
            if ct eq 0 then begin ; not already linked in previous loops
                LinkAreas,list.mask,list.masknr,IndMask[i]
                TableMean=[[TableMean],[IndMask[i],WMean((*Profiles[i])[0,*],(*Profiles[i])[1,*])]]
                AvoidArray=[AvoidArray,i]
            endif
        endif else begin ; overlap
            groupind=[i,index1[ind1]]

            ; Exclude the ones that are already linked in previous loops
            ind=-1
            for j=0l,n_elements(groupind)-1 do begin
                tmp=where(AvoidArray eq groupind[j],ct)
                if ct eq 0 then ind=[ind,groupind[j]]
            endfor

            if n_elements(ind) ne 1 then begin  ; at least one area not already linked in previous loops
                ind=ind[1:*]
                LinkAreas,list.mask,list.masknr,IndMask[ind]
                TableMean=[[TableMean],[IndMask[ind[0]],OverlapnProfiles(Profiles[ind])]]
                AvoidArray=[AvoidArray,ind]
            endif
        endelse

    endelse
endfor

nGroups=(size(TableMean,/dimensions))[1]-1
TableMean[0,*]=(*list.mask)[1,TableMean[0,*]]
TableMean[2,*]=sqrt(TableMean[2,*])

maskgroup=long(reform((*list.mask)[1,IndMask]))
tmp=(long(reform(TableMean[0,*])))[1:*]
for i=0l,nProf-1 do begin
    ind=where(tmp eq maskgroup[i])
    maskgroup[i]=ind
endfor

; ---- Output ----
tmp=CutPath(list.mskfile,file=file)
path=SelFile(list.mskpath,file+'.txt','*.*','Save overlap table...')

if path ne '' then $
if ~openw_safe(lun,path) then begin
printf,lun,'$nProf:'
printf,lun,nProf
printf,lun,'$XrowByYcol:'
printf,lun,Table,format='('+stringr(nProf)+'f15.7)'
printf,lun,'$XY:'
printf,lun,Table2,format='('+stringr(nProf)+'f15.7)'
printf,lun,'$Profiles:'
for i=0l,nProf-1 do begin
    n=n_elements(*Profiles[i])/2
    printf,lun,n
    printf,lun,transpose(*Profiles[i]),format='('+stringr(n)+'f15.7)'
endfor
printf,lun,'$Grouping:'
printf,lun,list.OverLapCrit
printf,lun,IndMask
printf,lun,maskgroup
printf,lun,'$',ChiXtitle(list.AzIntInfo.type),' (groupnr, d, SD):'
printf,lun,nGroups
tmp=strlen(stringr(fix(list.masknr)))+1
printf,lun,tmp
printf,lun,TableMean[*,1:*],format='(I'+stringr(tmp)+',2f15.7)'
free_lun,lun
endif

PTR_FREE,Profiles
wdelete,!D.window
RefreshDisplay,top,list

end;pro CreateOverlapTable
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_AreaRange,ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=listt
widget_control,listt.top,sensitive=1
widget_control,listt.top,get_uvalue=list

AreaRangeStruc={i:listt.unit,b:[listt.azi,listt.auto],str:listt.str}
list.AreaRangeStruc.i=listt.unit
list.AreaRangeStruc.b=[listt.azi,listt.auto]
list.AreaRangeStruc.str=listt.str

widget_control,listt.top,set_uvalue=list
end;pro CleanUp_AreaRange
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AreaRange_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list2
widget_control,list2.top,get_uvalue=list
widget_control,ev.id,get_value=val

case widget_info(ev.id,/type) of
1 : begin ;button event
    case val of
    'OK' :    begin
            ; ---- radius range ----
            d=[list2.str[0]<list2.str[1],list2.str[0]>list2.str[1]]
            case list2.unit of
            0:    d=BraggPtoX(d/(list.scrx*1000),list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/onlyd); mm -> angstrom
            1:    
            2:    d=BraggTtoX(d,list.lambda,/onlyd,/angledeg) ; degrees -> angstrom
            3:    d=BraggQtoX(d,list.lambda,/onlyd) ; 1/nm -> angstrom
            endcase

            ; ---- Mask ----
            if list2.azi eq 1 then begin
                ; ---- azimuthal range ----
                azi1=list2.str[2]
                azi2=list2.str[3]
                azi=[azi1,azi2]*!dpi/180

                out=[2,0,list2.auto,d[0],azi[0],d[1],azi[1]]

            endif else begin
                out=[3,0,list2.auto,d[0],d[1],0,0]
            endelse

            ; ---- save and exit ----
            if list.masknr eq 0 then begin ; New
                out[1]=0
                list.mask=ptr_new(out)
                list.smask=ptr_new([1b])
                list.masknr++
            endif else begin
                ind=where((*list.mask)[2,*] eq 0,count)

                if count eq 0 then begin ; Add
                    L=max((*list.mask)[1,*])+1
                    out[1]=L
                    (*list.mask)=[[(*list.mask)],[out]]
                    (*list.smask)=[(*list.smask),1b]
                    list.masknr++
                endif else begin ; Replace last
                    out[1]=(*list.mask)[1,ind[0]]
                    (*list.mask)[*,ind[0]]=out
                endelse
            endelse

              RefreshDisplay,list2.top,list

            WIDGET_CONTROL, ev.top, /DESTROY
            endcase
    else :    begin
            widget_control,ev.id,get_uvalue=uval
            case uval of
            'mm':         list2.unit=0
            msymbols('angstroms'): list2.unit=1
            'degrees':    list2.unit=2
            '1/nm':        list2.unit=3
            'auto':        list2.auto=ev.select
            'arc':        begin
                        list2.azi=1
                        ID=widget_info(ev.top,FIND_BY_UNAME='azibase')
                        widget_control,ID,sensitive=1
                        endcase
            'Ellipse':    begin
                        list2.azi=0
                        ID=widget_info(ev.top,FIND_BY_UNAME='azibase')
                        widget_control,ID,sensitive=0
                        endcase
            endcase
            widget_control,ev.top,set_uvalue=list2
            endelse
    endcase
    endcase
3 : begin ;text event
    widget_control,ev.id,get_uvalue=i
    f=0.
    reads,val,f
    list2.str[i]=f
    widget_control,ev.top,set_uvalue=list2
    endcase
endcase

end;pro AreaRange_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AreaRangeWindow,top

widget_control,top,get_uvalue=list
widget_control,top,sensitive=0

listt={top:top,unit:list.AreaRangeStruc.i,azi:(list.AreaRangeStruc.b)[0],auto:(list.AreaRangeStruc.b)[1],str:list.AreaRangeStruc.str}
topbase=widget_base(/column,title='Add area:',uvalue=listt)
xs=100

tbase=widget_base(topbase,/row,/exclusive)
    buttonp0=widget_button(tbase,value='mm',uvalue='mm',xsize=xs)
    buttonp1=widget_button(tbase,value=msymbols('angstroms'),uvalue=msymbols('angstroms'),xsize=xs)
    buttonp2=widget_button(tbase,value='degrees',uvalue='degrees',xsize=xs)
    buttonp3=widget_button(tbase,value='1/nm',uvalue='1/nm',xsize=xs)

tbase=widget_base(topbase,/row)
    label=widget_label(tbase,value='Begin value:',xsize=xs)
    text=widget_text(tbase,value=stringr((listt.str)[0]),scr_xsize=xs,/editable,uvalue=0)

tbase=widget_base(topbase,/row)
    label=widget_label(tbase,value='End value:',xsize=xs)
    text=widget_text(tbase,value=stringr((listt.str)[1]),scr_xsize=xs,/editable,uvalue=1)

tbase=widget_base(topbase,/row,/exclusive)
    buttonc0=widget_button(tbase,value='Ellipse',uvalue='Ellipse',xsize=xs)
    buttonc1=widget_button(tbase,value='arc',uvalue='arc',xsize=xs)

tbase=widget_base(topbase,/row,uname='azibase',sensitive=listt.azi)
    label=widget_label(tbase,value='Azimuth('+msymbols('degrees')+'):',xsize=xs)
    text=widget_text(tbase,value=stringr((listt.str)[2]),scr_xsize=xs,/editable,uvalue=2)
    text=widget_text(tbase,value=stringr((listt.str)[3]),scr_xsize=xs,/editable,uvalue=3)

tbase=widget_base(topbase,/row,/nonexclusive)
    buttona=widget_button(tbase,value='Auto set area',uvalue='auto',xsize=xs)

button=widget_button(topbase,value='OK')

WIDGET_CONTROL, topbase, /REALIZE
case listt.unit of
0:    widget_control,buttonp0,set_button=1
1:    widget_control,buttonp1,set_button=1
2:    widget_control,buttonp2,set_button=1
3:    widget_control,buttonp3,set_button=1
endcase
case listt.azi of
0: widget_control,buttonc0,set_button=1
1: widget_control,buttonc1,set_button=1
endcase
widget_control,buttona,set_button=listt.auto

Xmanager,'AreaRangeWindow',topbase, event_handler='AreaRange_event',$
cleanup='CleanUp_AreaRange',GROUP_LEADER=top,/no_block
end;pro AreaRangeWindow
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_PointOptions,ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
widget_control,list.top,sensitive=1
end;pro CleanUp_PointOptions
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PointOptions_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list2
widget_control,list2.top,get_uvalue=list

case widget_info(ev.id,/type) of
1 : begin ;button event
    widget_control,ev.id,get_value=val
    case val of
    'OK' :    begin
            WIDGET_CONTROL, ev.top, /DESTROY
            return
            endcase
    else :    begin
            widget_control,ev.id,get_uvalue=uval
            case uval of
            'buttonmass':    list.DelMass=ev.select
            'buttonmslider': list.cslideruseM=ev.select
            'UpMain':list.CrossUpdate=ev.select
            endcase
            endelse
    endcase
    endcase
3 : begin ;text event
    widget_control,ev.id,get_value=val,get_uvalue=uval
    in=0.

    if uval ne 12 then reads,val,in
    case uval of
    0:    list.letter=in
    1:    list.CW1=in
    2:    list.CW2=in
    3:    begin
        list.DelWidth=in
        ID=widget_info(ev.top,FIND_BY_UNAME='DelArea')
        WIDGET_CONTROL, ID, /DESTROY
        ID=widget_info(ev.top,FIND_BY_UNAME='drawbase')
        w=2*list.DelWidth+1
        base=widget_draw(ID,xsize=w,ysize=w,uname='DelArea')
        endcase
    4:    begin
        if (in ne 0) then begin
            list.cslider=list.cslider*list.cslider[2]/in
            list.cslider[2]=in
            ID=widget_info(list2.top,FIND_BY_UNAME='circle slider')
            WIDGET_CONTROL, ID,SET_SLIDER_MIN=list.cslider[0]
            WIDGET_CONTROL, ID,SET_SLIDER_MAX=list.cslider[1]
            WIDGET_CONTROL, ID,SET_VALUE=list.cslider[3]
            ID=widget_info(list2.top,FIND_BY_UNAME='csliderlabel')
            WIDGET_CONTROL, ID,SET_VALUE='Azimuth: '+stringr(list.cslider[3]*list.cslider[2])+msymbols('degrees')
        endif
        endcase
    5:    begin
        if (in ne 0) then begin
            list.clipslider*=list.clipslider[2]/in
            list.clipslider[3:4]=list.clipslider[0]>list.clipslider[3:4]<list.clipslider[1]
            list.clipslider[2]=in
            list.sclmax=2.^(list.clipslider[3:4]*list.clipslider[2])-1
            RefreshDisplay,list2.top,list
            RefreshScalingSlider,list2.top
            return
        endif
        endcase
    6:    begin
        if (in ne 0) then begin
            list.clipslider[1]=in/list.clipslider[2]
            list.clipslider[3:4]=list.clipslider[0]>list.clipslider[3:4]<list.clipslider[1]
            list.sclmax=2.^(list.clipslider[3:4]*list.clipslider[2])-1
            RefreshDisplay,list2.top,list
            RefreshScalingSlider,list2.top
            return
        endif
        endcase
    7:    list.dplotthick=in
    8:    list.RadAxMin=in
    9:    list.RadAxMax=in
    10:    list.RadAxAzimuth=in
    11:    list.RadAxBins=in>1
    12:    list.RadAxFormat=val
    endcase
    endcase
8:    begin
    widget_control,ev.id,get_uvalue=uval
    case uval of
    1:    list.RadAxType=ev.index
    endcase
    endcase
else: return
endcase

RefreshDisplay,list2.top,list

end;pro PointOptions_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PointOptions,top

widget_control,top,get_uvalue=list
widget_control,top,sensitive=0

topbase=widget_base(/column,title='Interface Options',uvalue={top:top})
xs=200

tab=widget_tab(topbase)
tab1=widget_base(tab,/column,title='Misc')

tbase=widget_base(tab1,/row)
    label=widget_label(tbase,value='Size of labels',xsize=xs)
    text=widget_text(tbase,value=stringr(list.letter),scr_xsize=xs,/editable,uvalue=0)

tbase=widget_base(tab1,/row)
    label=widget_label(tbase,value='Dynamic window line thickness',xsize=xs)
    text=widget_text(tbase,value=stringr(list.dplotthick),scr_xsize=xs,/editable,uvalue=7)
    
tbase=widget_base(tab1,/row)
    label=widget_label(tbase,value='Size points in static window',xsize=xs)
    text=widget_text(tbase,value=stringr(list.CW1),scr_xsize=xs,/editable,uvalue=1)

tbase=widget_base(tab1,/row)
    label=widget_label(tbase,value='Size points in dynamic window (pixels)',xsize=xs)
    text=widget_text(tbase,value=stringr(list.CW2),scr_xsize=xs,/editable,uvalue=2)

tbase=widget_base(tab1,/row,uname='drawbase')
    label=widget_label(tbase,value='Point selection width (pixels)',xsize=xs)
    text=widget_text(tbase,value=stringr(list.DelWidth),scr_xsize=xs,/editable,uvalue=3)
    base=widget_draw(tbase,xsize=list.DelWidth,ysize=list.DelWidth,uname='DelArea')

tbase=widget_base(tab1,/row)
    label=widget_label(tbase,value='Circle slider increment ('+msymbols('degrees')+')',xsize=xs)
    text=widget_text(tbase,value=stringr(list.cslider[2]),scr_xsize=xs,/editable,uvalue=4)

tbase=widget_base(tab1,/row)
    label=widget_label(tbase,value='Image clipping (min/max)= 2^(clipslidervalue*increment)-1')
tbase=widget_base(tab1,/row)
    label=widget_label(tbase,value='Image clip slider increment',xsize=xs)
    text=widget_text(tbase,value=stringr(list.clipslider[2]),scr_xsize=xs,/editable,uvalue=5)

tbase=widget_base(tab1,/row)
    label=widget_label(tbase,value='Image depth for clipping (bits)',xsize=xs)
    text=widget_text(tbase,value=stringr(list.clipslider[1]*list.clipslider[2]),scr_xsize=xs,/editable,uvalue=6)

tbase=widget_base(tab1,/column,/nonexclusive)
    buttonmslider=widget_button(tbase,value='Use Debye Marker for circle slider',uvalue='buttonmslider')
    buttonmass=widget_button(tbase,value='Delete multiple points without prompt',uvalue='buttonmass')
    buttonupmain=widget_button(tbase,value='Update 1D window',uvalue='UpMain')
    
tab2=widget_base(tab,/column,title='Radial Axis')
    
tbase=widget_base(tab2,/row)
    label=widget_label(tbase,value='Radial axis type',xsize=xs)
    RadAxType=widget_droplist(tbase,value=['Radial (mm)','2-theta (deg)','D-spacing ('+msymbols('angstroms')+')','Q (1/nm)'],uvalue=1,scr_xsize=xs)

tbase=widget_base(tab2,/row)
    label=widget_label(tbase,value='Radial axis minimum',xsize=xs)
    text=widget_text(tbase,value=stringr(list.RadAxMin),scr_xsize=xs,/editable,uvalue=8)
    
tbase=widget_base(tab2,/row)
    label=widget_label(tbase,value='Radial axis maximum',xsize=xs)
    text=widget_text(tbase,value=stringr(list.RadAxMax),scr_xsize=xs,/editable,uvalue=9)

tbase=widget_base(tab2,/row)
    label=widget_label(tbase,value='Radial axis bins',xsize=xs)
    text=widget_text(tbase,value=stringr(list.RadAxBins),scr_xsize=xs,/editable,uvalue=11)
    
tbase=widget_base(tab2,/row)
    label=widget_label(tbase,value='Radial axis azimuth',xsize=xs)
    text=widget_text(tbase,value=stringr(list.RadAxAzimuth),scr_xsize=xs,/editable,uvalue=10)
    
tbase=widget_base(tab2,/row)
    label=widget_label(tbase,value='Radial axis format',xsize=xs)
    text=widget_text(tbase,value=list.RadAxFormat,scr_xsize=xs,/editable,uvalue=12)

button=widget_button(topbase,value='OK')

WIDGET_CONTROL, topbase, /REALIZE
widget_control,buttonmass,set_button=list.DelMass
widget_control,buttonmslider,set_button=list.cslideruseM
widget_control,buttonupmain,set_button=list.CrossUpdate
widget_control,RadAxType,SET_DROPLIST_SELECT=list.RadAxType

Xmanager,'PointOptions',topbase, event_handler='PointOptions_event',$
cleanup='CleanUp_PointOptions',GROUP_LEADER=top,/no_block

end;pro PointOptions
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CMethodAllow,done,SeqTags,CorSeq,CorMethod,tag
; bool: changing method allowed?

; k is the index of the method to change
ind=where(SeqTags eq tag)
k=ind[0]

; j is the index of the last operation performed
j=n_elements(done)-1
while ((done[(*CorSeq).(j)] eq 0) and (j ne 0)) do j=j-1

; If k is less then j then we can't do change
if (k lt j) then return,0B

return,(done[(*CorSeq).(k)] eq 0)

end;function CMethodAllow
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CorrAllow,done,SeqTags,CorSeq,CorMethod,tag,Rest=Rest
; Rest: restoration allowed?
; bool: correction allowed?

; k is the index of the operation to (un)do
ind=where(SeqTags eq tag)
k=ind[0]

; j is the index of the last operation performed
j=n_elements(done)-1
while ((done[(*CorSeq).(j)] eq 0) and (j ne 0)) do j=j-1

; If k is less then j then we can't do anything (do nor undo)
if (k lt j) then begin
    Rest=0B
    return,0B
endif

; Check whether all operations before are performed or not required
bool=1B
if (k gt 0) then $
for i=0l,k-1 do begin
    bool=bool and (done[(*CorSeq).(i)] or ((CorMethod[(*CorSeq).(i)]<1) eq 0))
endfor
Rest=done[(*CorSeq).(k)]
bool=bool and (Rest eq 0)
return,bool

end;function CorrAllow
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_Corrections, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
widget_control,list.top,sensitive=1
mode,{top:list.top},list.mode,/refresh
end;pro CleanUp_Corrections
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CorrOptPaint,ev,list,n

; Everything before j: sensitive=0
j=n-1
while ((list.CorDone[(*(list.CorSeq)).(j)] eq 0) and (j ne 0)) do j=j-1

; Set sensitive and button select
bool=1
for i=0l,n-1 do begin
      ID=widget_info(ev.top,FIND_BY_UNAME=stringr(i))
       done=list.CorDone[(*(list.CorSeq)).(i)]
       widget_control,ID,sensitive=bool and (i ge j),set_button=done
    ID=widget_info(ev.top,FIND_BY_UNAME='d'+stringr(i))
    widget_control,ID,sensitive=(done eq 0) and (i ge j)
    ID=widget_info(ev.top,FIND_BY_UNAME='e'+stringr(i))
    widget_control,ID,sensitive=(done eq 0) and (i ge j)
    bool=bool and (done or ((list.CorMethod[(*(list.CorSeq)).(i)]<1) eq 0))
endfor

end;pro CorrOptPaint
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AutoCorrect,list
printw,list.OutID,'AutoCorrect...'
ind=where(list.CorDone eq 1,ct)
nlast=n_elements(list.CorMethod)-1
nfirst=0
if ct gt 0 then nfirst=ind[ct-1]+1
if nfirst ge nlast then begin
    printw,list.OutID,'No corrections left'
endif else begin
    for i=nfirst,nlast do begin
        index=(*(list.CorSeq)).(i)
        if (list.CorMethod[index] ne 0) then begin
        case list.SeqTags[i] of
        'zinger':    begin ;Zinger
                T=systime(1)
                tmp=list.tiffvalid
                  RemoveZingers,list.tiff,tmp,width=list.ZingerWidth,threshold=list.ZingerThreshold,count=count
                  list.tiffvalid=tmp
                printw,list.OutID,'Zingers found: '+stringr(count)+' in '+string(systime(1)-T)+' sec'
                endcase
        'satur':    begin ;Saturation
                T=systime(1)
                tmp=list.tiffvalid
                RemoveSatur,list.tiff,tmp,list.SaturMax,list.SaturAspect,list.SaturNMax,count=count
                list.tiffvalid=tmp
                printw,list.OutID,'Saturated pixels found: '+stringr(count)+' in '+string(systime(1)-T)+' sec'
                endcase
        'back':    begin ;Background
                T=systime(1)
                bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
                case bkcode[0] of
                0:    begin
                    endcase
                1:    begin
                    ptr_free,list.tiffb
                    list.tiffb=Strip2D(list.tiff,list.wb,list.niter)
                      BackSub,list.tiff,list.tiffb,bkcode[0]
                    endcase
                2:    begin
                    if ~ptr_valid(tiffb) then begin
                        tiffb=ReadCCD(list.darkpath,list.darkfile,info=list.readinfo)
                        if ptr_valid(tiffb) eq 0 then begin
                           printw,list.OutID,'Dark image not found:'+list.darkpath+list.darkfile
                           return
                        endif
                        ptr_free,list.tiffb
                        if list.bcko eq 0 then list.tiffb=tiffb $
                          else begin
                              list.tiffb=mSavgol(data=tiffb,r=list.bcko,nr=list.bckw/2,nl=list.bckw/2,order=[0,0])
                              ptr_free,tiffb
                          endelse
                          (*list.tiffb)*=list.bckm
                      endif
                      BackSub,list.tiff,list.tiffb,bkcode[0]
                    endcase
                endcase
                if bkcode[1] then begin
                      ind=where((*list.tiff) le list.thres[0],ct)
                      if ct ne 0 then (*list.tiff)[ind]=0
                  endif
                  if bkcode[2] then begin
                      ind=where((*list.tiff) ge list.thres[1],ct)
                      if ct ne 0 then (*list.tiff)[ind]=0
                  endif
                  printw,list.OutID,'Background subtraction: time (sec):'+string(systime(1)-T)
                endcase
        'spatial':    begin ;Spatial distortion
                T=systime(1)
                
                if not PTR_VALID(list.SpaDisStruc) then begin
                    if list.Tiegrid eq 2 then begin
                        if not PTR_VALID(list.tck1) then begin
                            struc=ReadFit2DSpline(list.TiePath+list.TieFile,error=error)
                            if error then begin
                                printw,list.OutID,'Error on loading '+list.TiePath+list.TieFile
                            endif else begin
                                list.GridDim=struc.GridDim
                                list.GridOffset=struc.GridOffset
                                list.GridOrigin=struc.GridOrigin
                                list.GridTilt=struc.GridTilt
                                list.GridSpReal=struc.GridSpReal
                                list.GridSpPix=struc.GridSpPix
                                list.scrx=struc.scrx
                                list.scry=struc.scry
                                PTR_FREE,list.tck1
                                PTR_FREE,list.tck2
                                list.tck1=PTR_NEW(struc.tck1)
                                list.tck2=PTR_NEW(struc.tck2)
                            endelse
                          endif
                          if PTR_VALID(list.tck1) then $
                              list.SpaDisStruc=WarpSpaDistPrepSpline(list.tiffs,*list.tck1,*list.tck2,vflip=list.vFlipSpline)
                    endif else begin
                        if (list.nrTieRing eq 0) and ((*list.nrTie)[0] eq 0) then begin
                            printw,list.OutID,'No tie points defined'
                        endif else begin
                            if list.Tiegrid then TieO=list.TieO $
                            else TieO=CalcTieO(list.nrTieRing,list.nrTie,list.TieI,list.center,list.lambda,list.dist,list.scrx,list.scry,list.TieRingD,list.a,list.b)
                            list.SpaDisStruc=WarpSpaDistPrep(list.tiffs,list.TieI,TieO,fast=list.bFastSpline)
                            if list.Tiegrid eq 0 then ptr_free,TieO
                        endelse
                    endelse
                endif
                
                if not PTR_VALID(list.SpaDisStruc) then return
                tmp=list.tiffvalid
                if list.Tiegrid eq 2 then WarpSpaDistRunSpline,list.tiff,list.SpaDisStruc,tiffvalid=tmp,fast=list.bFastSpline $
                else WarpSpaDistRun,list.tiff,list.SpaDisStruc,tiffvalid=tmp
                list.tiffvalid=tmp
                
                t2=systime(1)
                printw,list.OutID,'Spatial distortion: time (sec):'+string(t2-T)
                endcase
        'ff':    begin ;Flat field
                T=systime(1)
                
                if ~ptr_valid(list.ff) then begin
                    ff=ReadCCD(list.ffpath,list.fffile,info=list.readinfo)
                    if ptr_valid(ff) eq 0 then begin
                        printw,list.OutID,'Flat Field image not found:'+list.ffpath+list.fffile
                        return
                    endif
                    (*ff)=mean(*ff)/(*ff)
                      ptr_free,list.ff
                      list.ff=ff
                  endif
                  
                  type=size(*list.tiff,/type)
                (*list.tiff)*=(*list.ff)
                ;p=TypeConvert(list.tiff,type)
                t2=systime(1)
                printw,list.OutID,'Flat Field: time (sec):'+string(t2-T)
                endcase
        'maskoff':begin ;Mask off
                T=systime(1)
                
                if ~ptr_valid(list.maskoff) then begin
                    maskoff=ReadCCD(list.maskoffpath,list.maskofffile,info=list.readinfo)
                    if ptr_valid(maskoff) eq 0 then begin
                        if ~list.maskoffdyn then begin
                            printw,list.OutID,'Mask image not found:'+list.maskoffpath+list.maskofffile
                            return
                        endif
                    endif else begin
                        compareself,maskoff,list.maskoffthres
                        ptr_free,list.maskoff
                        list.maskoff=maskoff
                    endelse
                endif
                
                if list.maskoffdyn then begin
                    if not ptr_valid(list.tiffvalid) then list.tiffvalid = ptr_new(compare(list.tiff,list.maskoffdynthres)) $
                    else (*list.tiffvalid) and= compare(list.tiff,list.maskoffdynthres)
                endif
                
                if ptr_valid(list.maskoff) then begin
                    if not ptr_valid(list.tiffvalid) then list.tiffvalid = ptr_new(*list.maskoff) $
                    else (*list.tiffvalid) and= (*list.maskoff)
                endif
        
                t2=systime(1)
                printw,list.OutID,'Mask off: time (sec):'+string(t2-T)
                endcase
            endcase
            list.CorDone[index]=1
            endif
    endfor
endelse
end;pro AutoCorrect
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReadWithAutoStart,list

; Make sure the current corrections are made for all files
list.CorMethod*=list.CorDone
ind=where((list.CorDone eq 1) and (list.CorMethod eq 0),ct) ; exception: background -> if done, method=1 or 2
if ct ne 0 then list.CorMethod[ind]=1

list.tiffs=ptr_new(*list.tiffs)
list.tiff=ptr_new(*list.tiff)
if ptr_valid(list.tiffb) then list.tiffb=ptr_new(*list.tiffb)
if ptr_valid(list.ff) then list.ff=ptr_new(*list.ff)
if ptr_valid(list.maskoff) then list.maskoff=ptr_new(*list.maskoff)
if ptr_valid(list.tiffvalid) then list.tiffvalid=ptr_new(*list.tiffvalid)
end;pro ReadWithAutoStart
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReadWithAutoStop,list
ptr_free,list.tiffs
ptr_free,list.tiff
ptr_free,list.tiffb
ptr_free,list.ff
ptr_free,list.maskoff
ptr_free,list.tiffvalid
end;pro ReadWithAutoStop
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadWithAutoCorrect,pathccd,list,path=path,file=file
tiff=ReadCCD(pathccd,info=list.readinfo)
if ptr_valid(tiff) then begin
    tiffs=size(*tiff)

    ; THIS MESSES WITH LIST, MAKE SURE IT'S NOT SAVED!!!!!!

    res=CutPath(pathccd,path=path,file=file,ext=ext)
    list.path=path
    list.file=file+ext
    ptr_free,list.tiff
    list.tiff=tiff
    *list.tiffs=tiffs

    list.CorDone=0b

    AutoCorrect,list
    return,1b
endif else return,0b

end;function ReadWithAutoCorrect
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Corrections_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
1 : begin ;button event

    WIDGET_CONTROL, ev.id, GET_VALUE=val
    case val of
    'OK':     begin
            WIDGET_CONTROL, ev.top, /DESTROY
            return
            endcase
    'Auto correct now':begin
            widget_control, ev.top, get_uvalue=list2
            widget_control,/HOURGLASS
            widget_control, list2.top, get_uvalue=list
            AutoCorrect,list
            RefreshDisplay,list2.top,list
            CorrOptPaint,ev,list,list2.n
            return
            endcase
    else:    begin
            WIDGET_CONTROL, ev.id, GET_UVALUE=uval
            widget_control, ev.top, get_uvalue=list2
            endcase
    endcase


    case uval of
    -1:    begin
        list2.from=long((str_sep(val,' '))[1])
        widget_control, ev.top, set_uvalue=list2
        endcase
    -2:    begin
        list2.to=long((str_sep(val,' '))[1])
        widget_control, ev.top, set_uvalue=list2
        endcase
    -3:    begin
        widget_control, list2.top, get_uvalue=list
        list.autocor=ev.select
        widget_control, list2.top, set_uvalue=list
        endcase
    else:begin
        widget_control, list2.top, get_uvalue=list
        widget_control,/HOURGLASS

        if ev.select then begin ; Correct image (Do)
            case list.SeqTags[uval] of
            'zinger': begin ;Zinger
                T=systime(1)
                tmp=list.tiffvalid
                  RemoveZingers,list.tiff,tmp,width=list.ZingerWidth,threshold=list.ZingerThreshold,count=count
                  list.tiffvalid=tmp
                printw,list.OutID,'Zingers found: '+stringr(count)+' in '+string(systime(1)-T)+' sec'
                endcase
            'satur': begin ;Saturation
                T=systime(1)
                tmp=list.tiffvalid
                RemoveSatur,list.tiff,tmp,list.SaturMax,list.SaturAspect,list.SaturNMax,count=count
                list.tiffvalid=tmp
                printw,list.OutID,'Saturated pixels found: '+stringr(count)+' in '+string(systime(1)-T)+' sec'
                endcase
            'back':    begin ;Background
                T=systime(1)
                bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
                case bkcode[0] of
                0:    begin
                    ;printw,list.OutID,'No subtraction method selected'
                    ;widget_control,ev.id,set_button=0
                    ;return
                    endcase
                1:    begin
                    ptr_free,list.tiffb
                    list.tiffb=Strip2D(list.tiff,list.wb,list.niter)
                      BackSub,list.tiff,list.tiffb,bkcode[0]
                    endcase
                2:    begin
                    if ~ptr_valid(list.tiffb) then begin
                        tiffb=ReadCCD(list.darkpath,list.darkfile,info=list.readinfo)
                        if ptr_valid(tiffb) eq 0 then begin
                           printw,list.OutID,'Dark image not found:'+list.darkpath+list.darkfile
                           widget_control,ev.id,set_button=0
                           return
                        endif
                        ptr_free,list.tiffb
                        if list.bcko eq 0 then list.tiffb=tiffb $
                          else begin
                              list.tiffb=mSavgol(data=tiffb,r=list.bcko,nr=list.bckw/2,nl=list.bckw/2,order=[0,0])
                              ptr_free,tiffb
                          endelse
                          (*list.tiffb)*=list.bckm
                      endif
                      BackSub,list.tiff,list.tiffb,bkcode[0]
                    endcase
                endcase
                if bkcode[1] then begin
                      ind=where((*list.tiff) le list.thres[0],ct)
                      if ct ne 0 then (*list.tiff)[ind]=0
                  endif
                  if bkcode[2] then begin
                      ind=where((*list.tiff) ge list.thres[1],ct)
                      if ct ne 0 then (*list.tiff)[ind]=0
                  endif
                  printw,list.OutID,'Background subtraction: time (sec):'+string(systime(1)-T)
                endcase
            'spatial':    begin ;Spatial distortion
                T=systime(1)
                
                if not PTR_VALID(list.SpaDisStruc) then begin
                    if list.Tiegrid eq 2 then begin
                        if not PTR_VALID(list.tck1) then begin
                            struc=ReadFit2DSpline(list.TiePath+list.TieFile,error=error)
                            if error then begin
                                printw,list.OutID,'Error on loading '+list.TiePath+list.TieFile
                            endif else begin
                                list.GridDim=struc.GridDim
                                list.GridOffset=struc.GridOffset
                                list.GridOrigin=struc.GridOrigin
                                list.GridTilt=struc.GridTilt
                                list.GridSpReal=struc.GridSpReal
                                list.GridSpPix=struc.GridSpPix
                                list.scrx=struc.scrx
                                list.scry=struc.scry
                                PTR_FREE,list.tck1
                                PTR_FREE,list.tck2
                                list.tck1=PTR_NEW(struc.tck1)
                                list.tck2=PTR_NEW(struc.tck2)
                            endelse
                          endif
                          if PTR_VALID(list.tck1) then $
                              list.SpaDisStruc=WarpSpaDistPrepSpline(list.tiffs,*list.tck1,*list.tck2,vflip=list.vFlipSpline)
                    endif else begin
                        if (list.nrTieRing eq 0) and ((*list.nrTie)[0] eq 0) then begin
                            printw,list.OutID,'No tie points defined'
                        endif else begin
                            if list.Tiegrid then TieO=list.TieO $
                            else TieO=CalcTieO(list.nrTieRing,list.nrTie,list.TieI,list.center,list.lambda,list.dist,list.scrx,list.scry,list.TieRingD,list.a,list.b)
                            list.SpaDisStruc=WarpSpaDistPrep(list.tiffs,list.TieI,TieO,fast=list.bFastSpline)
                            if list.Tiegrid eq 0 then ptr_free,TieO
                        endelse
                    endelse
                endif
                
                if not PTR_VALID(list.SpaDisStruc) then begin
                    widget_control,ev.id,set_button=0
                    return
                endif
                tmp=list.tiffvalid
                if list.Tiegrid eq 2 then WarpSpaDistRunSpline,list.tiff,list.SpaDisStruc,tiffvalid=tmp,fast=list.bFastSpline $
                else WarpSpaDistRun,list.tiff,list.SpaDisStruc,tiffvalid=tmp
                list.tiffvalid=tmp

                t2=systime(1)
                printw,list.OutID,'Spatial distortion: time (sec):'+string(t2-T)
                endcase
            'ff':    begin ;Flat field
                T=systime(1)
                
                if ~ptr_valid(list.ff) then begin
                    ff=ReadCCD(list.ffpath,list.fffile,info=list.readinfo)
                    if ptr_valid(ff) eq 0 then begin
                        printw,list.OutID,'Flat Field image not found:'+list.ffpath+list.fffile
                        widget_control,ev.id,set_button=0
                        return
                    endif
                    (*ff)=mean(*ff)/(*ff)
                      ptr_free,list.ff
                      list.ff=ff
                  endif
                  
                type=size(*list.tiff,/type)
                (*list.tiff)*=(*list.ff)
                ;p=TypeConvert(list.tiff,type)
                t2=systime(1)
                printw,list.OutID,'Flat Field: time (sec):'+string(t2-T)
                endcase
            'maskoff':begin ;Mask off
                T=systime(1)
                if ~ptr_valid(list.maskoff) then begin
                    maskoff=ReadCCD(list.maskoffpath,list.maskofffile,info=list.readinfo)
                    if ptr_valid(maskoff) eq 0 then begin
                        if ~list.maskoffdyn then begin
                            printw,list.OutID,'Mask image not found:'+list.maskoffpath+list.maskofffile
                            widget_control,ev.id,set_button=0
                            return
                        endif
                    endif else begin
                        compareself,maskoff,list.maskoffthres
                        ptr_free,list.maskoff
                        list.maskoff=maskoff
                    endelse
                endif
    
                if list.maskoffdyn then begin
                    if not ptr_valid(list.tiffvalid) then list.tiffvalid = ptr_new(compare(list.tiff,list.maskoffdynthres)) $
                    else (*list.tiffvalid) and= compare(list.tiff,list.maskoffdynthres)
                endif
                
                if ptr_valid(list.maskoff) then begin
                    if not ptr_valid(list.tiffvalid) then list.tiffvalid = ptr_new(*list.maskoff) $
                    else (*list.tiffvalid) and= (*list.maskoff)
                endif
                
                t2=systime(1)
                printw,list.OutID,'Mask off: time (sec):'+string(t2-T)
                endcase
            endcase
            list.CorDone[((*(list.CorSeq))).(uval)]=ev.select
            CorrOptPaint,ev,list,list2.n
        endif else begin ;Restore image (Undo)
            result=DIALOG_MESSAGE('This will undo all corrections made, proceed?',/question)
            case result of
            'Yes':    begin
                    RestoreCorrection,list
                    CorrOptPaint,ev,list,list2.n
                    endcase
            'No':    begin
                    widget_control,ev.id,set_button=1
                    endcase
            endcase
        endelse
        RefreshDisplay,list2.top,list
        endcase
    endcase
    endcase

8:    begin ;droplist
    WIDGET_CONTROL, ev.id, GET_UVALUE=uval
    widget_control, ev.top, get_uvalue=list2
    widget_control, list2.top, get_uvalue=list
    
    case strmid(uval,0,1) of
    'd':    begin
            uval=long(strmid(uval,1))
            
            if list.SeqTags[uval] eq 'back' then begin
                bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).(uval)])
                bkcode[0]=ev.index mod 3
                if ev.index eq 0 then bkcode[1:2]=0
                if ev.index eq 3 then bkcode[1:2]=1
                list.CorMethod[(*(list.CorSeq)).back]=Set2DbkCode(bkcode)
            endif else list.CorMethod[((*(list.CorSeq))).(uval)]=ev.index
            
            widget_control, list2.top, set_uvalue=list
            CorrOptPaint,ev,list,list2.n
            endcase
    'e':    begin
            widget_control, ev.top, get_uvalue=list2
            widget_control, list2.top, get_uvalue=list
            uval=long(strmid(uval,1))
            if ev.index eq uval then return

            ; Rearrange array
            indorder=indgen(list2.n)
            
            if ev.index gt uval then begin
                indorder[uval:ev.index]++
                indorder[ev.index]=uval
            endif else begin
                indorder[ev.index:uval]--
                indorder[ev.index]=uval
            endelse

            ; Make new CorSeq with SeqTags as fields
            SeqTags=list.SeqTags[indorder]
            print,list.SeqTags
            print,SeqTags
            CorSeq=CREATE_STRUCT(SeqTags[0], (*list.CorSeq).(indorder[0]))
            for i=1,list2.n-1 do CorSeq=CREATE_STRUCT(CorSeq,SeqTags[i], (*list.CorSeq).(indorder[i]))
            PTR_FREE, list.CorSeq
            list.CorSeq=PTR_NEW(CorSeq)
            list.SeqTags=SeqTags

            ; Check for 'new order' vs. 'corrections applied' conflicts
            ; change CorMethod to 'No' if necessary
            bool=0B
            for i=(list2.n-1),0,-1 do begin
                bool or= list.CorDone[CorSeq.(i)]
                if bool and (list.CorDone[CorSeq.(i)] eq 0) then $
                    list.CorMethod[CorSeq.(i)]=0
            endfor

            RefreshDisplay,list2.top,list
            WIDGET_CONTROL, ev.top, /DESTROY
            CorrectOptions,list2.top
            endcase
    endcase

    endcase
endcase

end; pro Corrections_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CorrectOptions,top

widget_control,top,get_uvalue=list
widget_control,top,sensitive=0
n=n_tags((*(list.CorSeq)))
names=list.SeqTags

topbase=widget_base(/column,title='Corrections',uvalue={top:top,n:n,from:-1,to:-1,mode:list.mode})

xs=160
baset=widget_base(topbase,/row)
    
    label=widget_label(baset,value='(Un)do  ')
    label=widget_label(baset,value='Correction type',xsize=xs)
    label=widget_label(baset,value='Method options  ')
    label=widget_label(baset,value='Order',/ALIGN_RIGHT)

; Everything before j: sensitive=0
j=n-1
while ((list.CorDone[(*(list.CorSeq)).(j)] eq 0) and (j ne 0)) do j=j-1

; Loop over corrections in order performindex
lonbut=lonarr(n)
donear=bytarr(n)
londrop=lonbut
londrop2=lonbut
val2=stringr(indgen(n))
dropar=donear
bool=1
for i=0l,n-1 do begin
    basew=widget_base(topbase,/row,uname='wrap'+stringr(i))
    
    
    baset=widget_base(basew,/column,/nonexclusive)
          lonbut[i]=widget_button(baset,value=' ',uvalue=i,uname=stringr(i),sensitive=bool and (i ge j))
          donear[i]=list.CorDone[(*(list.CorSeq)).(i)]
          
    case names[i] of
    'zinger':begin
             button=widget_label(basew,value='Zingers',xsize=xs)
             val=['Not required','Required']
             endcase
    'satur':     begin
             button=widget_label(basew,value='Saturation',xsize=xs)
             val=['Not required','Required']
             endcase
    'back':     begin
             button=widget_label(basew,value='Background',xsize=xs)
             val=['Not required','Stripping','Dark image','thresholding']
             endcase
    'spatial':     begin
             button=widget_label(basew,value='Spatial distortion',xsize=xs)
             val=['Not required','Required']
             endcase
    'ff':     begin
             button=widget_label(basew,value='Flat field',xsize=xs)
             val=['Not required','Required']
             endcase
   'maskoff':begin
             button=widget_label(basew,value='Mask off',xsize=xs)
             val=['Not required','Required']
             endcase
    else:
    endcase
    londrop[i]=widget_droplist(basew,value=val,sensitive=(donear[i] eq 0) and (i ge j),uname='d'+stringr(i),uvalue='d'+stringr(i))
    londrop2[i]=widget_droplist(basew,value=val2,sensitive=(donear[i] eq 0) and (i ge j),uname='e'+stringr(i),uvalue='e'+stringr(i))
    
    dropar[i]=list.CorMethod[(*(list.CorSeq)).(i)]
    if names[i] eq 'back' and dropar[i] eq 0 then begin
        bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).(i)])
        if bkcode[1] or bkcode[2] then dropar[i]=3
    endif

    bool=bool and (donear[i] or ((dropar[i]<1) eq 0))
endfor

baset=widget_base(topbase,/column,/nonexclusive)
    buttona=widget_button(baset,value='Auto correct pattern on load (performs only required ones)',uvalue=-3)
button=widget_button(topbase,value='Auto correct now')
button=widget_button(topbase,value='OK')

WIDGET_CONTROL, topbase, /REALIZE
for i=0l,n-1 do begin
      widget_control,lonbut[i],set_button=donear[i]
       widget_control,londrop[i],SET_DROPLIST_SELECT=dropar[i]
       widget_control,londrop2[i],SET_DROPLIST_SELECT=i
endfor
widget_control,buttona,set_button=list.autocor
Xmanager,'CorrectOptions',topbase, event_handler='Corrections_event',$
cleanup='CleanUp_Corrections',GROUP_LEADER=top,/no_block

end; pro CorrectOptions
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_ReadOptions, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
widget_control,list.top,sensitive=1
end;pro CleanUp_ReadOptions
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReadOptions_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_uvalue=uval
widget_control,ev.id,get_value=val
widget_control,ev.top,get_uvalue=list2
widget_control,list2.top,get_uvalue=list

case uval of
'showheader': list.readinfo.showheader=ev.select
'bnorm':    begin
            list.readinfo.bnorm=ev.select
            ID=widget_info(ev.top,find_by_uname='normbase')
            widget_control,ID,sensitive=ev.select
            endcase
'OK':         begin
            WIDGET_CONTROL, ev.top, /DESTROY
            return
            endcase
'int':         list.readinfo.int=val
'gain':     list.readinfo.gain=val
'offset':     list.readinfo.offset=val
'xs':         list.readinfo.binxs=fix(val)
'ys':         list.readinfo.binys=fix(val)
'endian':    list.readinfo.binlendian=ev.select
'type':        list.readinfo.bintype=([1,2,12,3,13,14,15,4,5])[ev.index]
'xpad3':    list.readinfo.binnxpad3=fix(val)
endcase

widget_control,list2.top,set_uvalue=list
end;pro ReadOptions_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReadOptions,top
widget_control,top,get_uvalue=list
widget_control,top,sensitive=0

topbase=widget_base(/row,title='Read options',uvalue={top:top})

baset=widget_base(topbase,/column,/nonexclusive)
buttonshowheader=widget_button(baset,value='Show DM3 header on reading',uvalue='showheader')
buttonbnorm=widget_button(baset,value='Normalize on reading',uvalue='bnorm')

normbase=widget_base(topbase,/column,uname='normbase',/frame,sensitive=list.readinfo.bnorm)
label=widget_label(normbase,value='normalization: image/(gain*intensity+offset)')
baset=widget_base(normbase,/row)
label=widget_label(baset,value='Intensity label:')
label=widget_text(baset,value=list.readinfo.int,uvalue='int',/editable)
baset=widget_base(normbase,/row)
label=widget_label(baset,value='Gain label:     ')
label=widget_text(baset,value=list.readinfo.gain,uvalue='gain',/editable)
baset=widget_base(normbase,/row)
label=widget_label(baset,value='Offset label:   ')
label=widget_text(baset,value=list.readinfo.offset,uvalue='offset',/editable)

binbase=widget_base(topbase,/column,uname='binbase',/frame)
label=widget_label(binbase,value='Binary format (.raw)')
baset=widget_base(binbase,/row)
label=widget_label(baset,value='XSize:')
label=widget_text(baset,value=stringr(list.readinfo.binxs),uvalue='xs',/editable)
baset=widget_base(binbase,/row)
label=widget_label(baset,value='YSize:')
label=widget_text(baset,value=stringr(list.readinfo.binys),uvalue='ys',/editable)
baset=widget_base(binbase,/row)
label=widget_label(baset,value='Type:')
dltype=widget_droplist(baset,value=['8bit','16bit singed','16bit unsinged','32bit singed','32bit unsinged',$
                    '64bit signed','64bit unsigned','32bit floating point','64bit floating point'],uvalue='type')
baset=widget_base(binbase,/row,/nonexclusive)
buttonlendian=widget_button(baset,value='Little Endian',uvalue='endian')
;baset=widget_base(binbase,/row)
;label=widget_label(baset,value='XPAD3 blocks:')
;label=widget_text(baset,value=stringr(list.readinfo.binnxpad3),uvalue='xpad3',/editable)

WIDGET_CONTROL, topbase, /REALIZE

widget_control,buttonshowheader,set_button=list.readinfo.showheader
widget_control,buttonbnorm,set_button=list.readinfo.bnorm
widget_control,buttonlendian,set_button=list.readinfo.binlendian

ind=where([1,2,12,3,13,14,15,4,5] eq list.readinfo.bintype)
widget_control,dltype,set_droplist_select=ind[0]

Xmanager,'ReadOptions',topbase, event_handler='ReadOptions_event',$
cleanup='CleanUp_ReadOptions',GROUP_LEADER=top,/no_block
end;pro ReadOptions
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro BatchProc2DTiff,ev,type

;type=0 -- all
;type=1 -- average
;type=2 -- superimpose
;type=3 -- summate

; Select files
widget_control,ev.top,get_uvalue=list
path=Select_Files(ev.top,list,outpath=pathred,outfile=file,sortmethod=list.sortmethod,separator=list.sortseparator)
if path[0] eq '' then return
path2=SelFile(pathred,file+list.Format,'*.*','Save Image....')
if path2 EQ '' then return
error=CutPath(path2,path=path3,file=file3,ext=ext)

n=n_elements(path)
Result = DIALOG_MESSAGE(stringr(n)+' files: proceed?',/question)
if result eq 'No' then return
widget_control,/hourglass

; Prevent messing with list
ReadWithAutoStart,list

; Read first file
if not ReadWithAutoCorrect(path[0],list) then begin
    ReadWithAutoStop,list
       printw,list.OutID,'First file can not be read.'
       return
endif
printw,list.OutID,path[0]

; Initialize output frames
tiff=ptrarr(3)
btiff=[type eq 1,type eq 2,type eq 3] or type eq 0
if btiff[0] then tiff[0]=ptr_new(double(*list.tiff)/n)
if btiff[1] then tiff[1]=ptr_new(*list.tiff)
if btiff[2] then tiff[2]=ptr_new(double(*list.tiff))

; Read other files
if n gt 0 then begin
    progressBar = Obj_New("PROGRESSBAR",title="Processing files...")
    progressBar -> Start

    for i=1ll,n-1 do begin
    
        IF progressBar -> CheckCancel() THEN BEGIN
            progressBar -> Destroy
            printw,list.outid,'Process CANCELLED!'
            RETURN
        ENDIF
    
        if ReadWithAutoCorrect(path[i],list) then begin
            ; add
            if btiff[0] then (*tiff[0])+=(*list.tiff)/double(n)
            if btiff[1] then (*tiff[1])>=(*list.tiff)
            if btiff[2] then (*tiff[2])+=(*list.tiff)
        endif
        printw,list.OutID,path[i]
        
        progressBar -> Update, (i+1.)/n*100
    endfor
    
    progressBar -> Destroy
endif

; Prevent messing with list
ReadWithAutoStop,list

; Save result
bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])

struc={center:list.center,dist:list.dist,file:list.file,scrx:list.scrx,$
    lambda:list.lambda,darkfile:(bkcode[0] eq 2)?list.darkfile:''}

if btiff[0] then begin
    ; Save average
    path4=path3+'avg'+file3+ext
    if WriteCCD(tiff[0],path4,struc=struc) then printw,list.OutID,'Pattern saved: '+path4 $
    else printw,list.OutID,'Pattern not saved: '+path4
endif

if btiff[1] then begin
    ; Save sup
    path4=path3+'sup'+file3+ext
    if WriteCCD(tiff[1],path4,struc=struc) then printw,list.OutID,'Pattern saved: '+path4 $
    else printw,list.OutID,'Pattern not saved: '+path4
endif

if btiff[2] then begin
    ; Save sum
    path4=path3+'sum'+file3+ext
    if WriteCCD(tiff[2],path4,struc=struc) then printw,list.OutID,'Pattern saved: '+path4 $
    else printw,list.OutID,'Pattern not saved: '+path4
endif

heap_free,tiff

end;pro BatchProc2DTiff
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DoublePresent,C,coord=coord,count=count
dind=-1
count=-1
coord=[-1,-1]

for i=0l,n_elements((*C))/2-1 do begin
    ind=where(((*C)[0,*] eq (*C)[0,i]) and ((*C)[1,*] eq (*C)[1,i]),ct)
    if ct gt 1 then begin
       ind2=where(ind lt i,ct)
       if ct eq 0 then begin
         dind=[dind,ind]
         coord=[[coord],[(*C)[*,i]]]
         count=[count,n_elements(ind)]
       endif
    endif
endfor

if n_elements(count) ne 1 then begin
    dind=dind[1:*]
    coord=coord[*,1:*]
    count=count[1:*]
endif else count=0

return,dind
end;function DoublePresent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GridNodeCoord,x,y,GridDim,GridSpReal,GridSpPix,GridOffset,GridOrigin,GridTilt
; convert node indices to TieO

; Convert grid node indices grid node positions in pixels
x*=GridSpPix[0]
y*=GridSpPix[1]
x+=GridOffset[0]
y+=GridOffset[1]

; Rotate over GridTilt
if GridTilt ne 0 then begin
    k=GridSpReal/GridSpPix
    if k[0] eq k[1] then k=[1,1] else k/=mean(k)
    cosa=cos(GridTilt*!dpi/180)
    sina=sin(GridTilt*!dpi/180)
    tmp=temporary(x)
    x=k[0]*cosa*tmp-k[1]*sina*y
    y=k[0]*sina*temporary(tmp)+k[1]*cosa*y
endif

; Shift the origin
x+=GridOrigin[0]
y+=GridOrigin[1]

end;pro GridNodeCoord
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GridNodeCoordInv,x,y,GridDim,GridSpReal,GridSpPix,GridOffset,GridOrigin,GridTilt
; convert TieO to node indices

; Shift the origin back
x-=GridOrigin[0]
y-=GridOrigin[1]

; Rotate over -GridTilt
if GridTilt ne 0 then begin
    k=GridSpReal/GridSpPix
    if k[0] eq k[1] then k=[1,1] else k/=mean(k)
    cosa=cos(GridTilt*!dpi/180)
    sina=sin(GridTilt*!dpi/180)
    tmp=temporary(x)
    x=cosa/k[0]*tmp+sina/k[0]*y
    y=-sina/k[1]*temporary(tmp)+cosa/k[1]*y
endif

; Pixel positions to grid node indices
x-=GridOffset[0]
y-=GridOffset[1]
x/=GridSpPix[0]
y/=GridSpPix[1]
x=round(x)
y=round(y)

end;pro GridNodeCoordInv
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GridCoordKeep,XY,nGridDim,nGridOffset,nGridOrigin,nGridTilt,nGridSpReal,nGridSpPix,oGridDim,oGridOffset,oGridOrigin,oGridTilt,oGridSpReal,oGridSpPix
; XY: TieO of the previous grid => TieO of the new grid
n=n_elements(*XY)/2
XYGrid=(*XY)*0

; Grid node indices of TieO
indx=(*XY)[0,*]
indy=(*XY)[1,*]
GridNodeCoordInv,indx,indy,oGridDim,oGridSpReal,oGridSpPix,oGridOffset,oGridOrigin,oGridTilt

; TieO of the new grid
xn=rebin(lindgen(nGridDim[0]),nGridDim,/sample)
yn=rebin(lindgen(1,nGridDim[1]),nGridDim,/sample)
GridNodeCoord,xn,yn,nGridDim,nGridSpReal,nGridSpPix,nGridOffset,nGridOrigin,nGridTilt

; Select the corresponding TieO
indx<=(nGridDim[0]-1)
indy<=(nGridDim[1]-1)
XYGrid[0,*]=xn[indx,indy]
XYGrid[1,*]=yn[indx,indy]

return,XYGrid

end;function GridCoordKeep
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GridCoord,XY,GridDim,GridOffset,GridOrigin,GridTilt,GridSpReal,GridSpPix
; XY=TieI => calculate corresponding TieO

n=n_elements(XY)/2
XYGrid=XY*0

; Calculate TieO
x=rebin(lindgen(GridDim[0]),GridDim,/sample)
y=rebin(lindgen(1,GridDim[1]),GridDim,/sample)
GridNodeCoord,x,y,GridDim,GridSpReal,GridSpPix,GridOffset,GridOrigin,GridTilt

; Find closest TieO for each TieI
for i=0l,n-1 do begin
    m=min((XY[0,i]-x)^2.+(XY[1,i]-y)^2.,ind)
    XYGrid[0,i]=x[ind]
    XYGrid[1,i]=y[ind]
endfor

return,XYGrid

end;function GridCoord
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AutoGridCoord,GridDim,GridOffset,GridOrigin,GridTilt,GridSpReal,GridSpPix

n=GridDim[0]*GridDim[1]

x=rebin(lindgen(GridDim[0]),GridDim,/sample)
y=rebin(lindgen(1,GridDim[1]),GridDim,/sample)
GridNodeCoord,x,y,GridDim,GridSpReal,GridSpPix,GridOffset,GridOrigin,GridTilt

return,[reform(x,1,n),reform(y,1,n)]

end;function AutoGridCoord
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AutoSortTie,p,GridDim

n=GridDim[0]*GridDim[1]
ind=reform(sort((*p)[1,*]),GridDim) ; group rows together
indkeep=0L
for j=0l,GridDim[1]-1 do begin
    ind2=ind[*,j]
    indkeep=[indkeep,ind2[sort((*p)[0,ind2])]] ; sort within row
endfor
*p=(*p)[*,indkeep[1:*]]

end;pro AutoSortTie
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MP_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]

IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent EQ 'DOWN' THEN BEGIN

    possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
    bPress = possibleButtons[ev.press]
    widget_control,ev.id,get_uvalue=list2
    if bPress eq 'RIGHT' then begin
       list2.point=list2.point+1 ; Do not change destination TieO
    endif else begin
       ; ----Check for zoom----
       if ev.id eq list.drawdyn then begin
       point=DynTR([ev.x,ev.y],list.sf,[list.hrange[0],list.vrange[0]])
       endif else begin
       point=StatTR([ev.x,ev.y],[list.sfh,list.sfv])
       endelse

       ;----Calc Output Tie Point----
       peakso=GridCoord(point,list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix)

       ;----Temp overcome rounding problems----
       ind=where(round((*list.TieO)[0,*]*1000) eq round(peakso[0]*1000) and $
                    round((*list.TieO)[1,*]*1000) eq round(peakso[1]*1000), ct)
       if ct ne 0 then peakso=(*list.TieO)[*,ind[0]]

       ;----Edit Output Point----
       (*list.TieO)[*,list2.dind[list2.point]]=peakso
       list2.point=list2.point+1
    endelse
    widget_control,list.draw,set_uvalue=list2
    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    booldrawdyn=widget_info(list.drawdyn,/VALID_ID)
    if booldrawdyn then begin
        widget_control,list.drawdyn,set_uvalue=list2
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    endif

    ; ----Still points left?----
    if list2.point lt n_elements(list2.dind) then begin
       list.xs=list2.dTieI[0,list2.point]
       list.ys=list2.dTieI[1,list2.point]
        RefreshDisplay, ev.top,list
       return
    endif
    widget_control,list.draw,draw_button=0,draw_motion=0
    if booldrawdyn then widget_control,list.drawdyn,draw_button=0,draw_motion=0

    ; ----Finish----
    ID=widget_info(list2.apptop,FIND_BY_UNAME='TieO')
    dind=DoublePresent(list.TieO,coord=coord)
    widget_control,ID,set_value=coord

    RefreshDisplay, ev.top,list
    widget_control,list2.apptop,sensitive=1
    widget_control,ev.top,sensitive=0
    return
endif

if ev.id eq list.drawdyn then begin
    pointx=[list.xs,DynTR(ev.x,list.sf,list.hrange[0])]
    pointy=[list.ys,DynTR(ev.y,list.sf,list.vrange[0])]
endif else begin
    pointx=[list.xs,StatTR(ev.x,list.sfh)]
    pointy=[list.ys,StatTR(ev.y,list.sfv)]
endelse

if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    x=RTDyn(pointx,list.sf,list.hrange[0])
    y=RTDyn(pointy,list.sf,list.vrange[0])
    arrow,x[0],y[0],x[1],y[1]
endif

WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
x=RTStat(pointx,list.sfh)
y=RTStat(pointy,list.sfv)
arrow,x[0],y[0],x[1],y[1]
end;pro MP_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MC_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN
possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
bPress = possibleButtons[ev.press]
if bPress eq 'RIGHT' then begin
    widget_control,ev.top,get_uvalue=list,sensitive=0
    Widget_Control, ev.id ,draw_button=0,get_uvalue=top
    widget_control,top,sensitive=1
    return
endif
widget_control,ev.top,get_uvalue=list
list.xs=ev.x
Widget_Control, ev.id, Event_Pro='MC_DRAWSEL', $
    draw_motion=1
Widget_Control, ev.top, Set_UValue=list
end;pro MC_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MC_DRAWSEL ,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF

IF thisEvent EQ 'DOWN' THEN BEGIN

    list.ys=list.xs
    list.yd=ev.x
    list.xs=(list.xs+ev.x)/2
    y=[ev.y,ev.y]
    x=[list.xs,ev.x]

    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
        PlotS, [x[0],x[0]],[0,list.dynv-1],/device
        PlotS, [x[1],x[1]],[0,list.dynv-1],/device
        arrow,x[0],y[0],x[1],y[1]
    endif

    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    x=DynTR(x,list.sf,list.hrange[0])
    y=DynTR(y,list.sf,list.vrange[0])
    x=RTStat(x,list.sfh)
    y=RTStat(y,list.sfv)
    vrange=RTStat(list.vrange,list.sfv)
    PlotS, [x[0],x[0]],vrange,/device
    PlotS, [x[1],x[1]],vrange,/device
    arrow,x[0],y[0],x[1],y[1]

    Widget_Control, ev.id, Event_Pro='MC_DRAWARROW'
    Widget_Control, ev.top, Set_UValue=list
    return
endif

y = [ev.y, ev.y]
x = [list.xs, ev.x]

if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    PlotS, [x[0],x[0]],[0,list.dynv-1],/device
    PlotS, [x[1],x[1]],[0,list.dynv-1],/device
    arrow,x[0],y[0],x[1],y[1]
endif

WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
x=DynTR(x,list.sf,list.hrange[0])
y=DynTR(y,list.sf,list.vrange[0])
x=RTStat(x,list.sfh)
y=RTStat(y,list.sfv)
vrange=RTStat(list.vrange,list.sfv)
PlotS, [x[0],x[0]],vrange,/device
PlotS, [x[1],x[1]],vrange,/device
arrow,x[0],y[0],x[1],y[1]

end;pro MC_DRAWSEL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MC_DRAWARROW,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent EQ 'DOWN' THEN BEGIN
    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    endif
    Widget_Control, ev.id, Event_Pro='MC_PROCESS_EVENTS'

    ; Process x coordinates
    xnew = ev.x
    xr = [list.ys,list.yd]

    xnew=DynTR(xnew,list.sf,list.hrange[0])
    xr=DynTR(xr,list.sf,list.hrange[0])

    ; Search for TieI within xr:
    ind=where(((*list.TieI)[0,*] gt (xr[0]<xr[1])) and ((*list.TieI)[0,*] lt (xr[0]>xr[1]))$
    and ((*list.TieI)[1,*] gt list.vrange[0]) and ((*list.TieI)[1,*] lt list.vrange[1]),ct)
    if ct eq 0 then return

    keep=(*list.TieI)[0,ind]
    (*list.TieI)[0,ind]=xnew
    (*list.TieO)[*,ind]=GridCoord((*list.TieI)[*,ind],$
           list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix)
    (*list.TieI)[0,ind]=keep
    
    RefreshDisplay, ev.top,list
    return
endif

y = [ev.y, ev.y]
x = [list.xs, ev.x]

if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    PlotS, [x[0],x[0]],[0,list.dynv-1],/device
    PlotS, [x[1],x[1]],[0,list.dynv-1],/device
    arrow,x[0],y[0],x[1],y[1]
endif

WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
x=DynTR(x,list.sf,list.hrange[0])
y=DynTR(y,list.sf,list.vrange[0])
x=RTStat(x,list.sfh)
y=RTStat(y,list.sfv)
vrange=RTStat(list.vrange,list.sfv)
PlotS, [x[0],x[0]],vrange,/device
PlotS, [x[1],x[1]],vrange,/device
arrow,x[0],y[0],x[1],y[1]

end;pro MC_DRAWARROW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MR_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN
possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
bPress = possibleButtons[ev.press]
if bPress eq 'RIGHT' then begin
    widget_control,ev.top,get_uvalue=list
    Widget_Control, ev.id ,draw_button=0
    widget_control,ev.id,get_uvalue=top
    widget_control,ev.top,sensitive=0
    widget_control,top,sensitive=1
    return
endif
widget_control,ev.top,get_uvalue=list
list.ys=ev.y
Widget_Control, ev.id, Event_Pro='MR_DRAWSEL', $
    draw_motion=1
Widget_Control, ev.top, Set_UValue=list
end;pro MR_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MR_DRAWSEL ,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent EQ 'DOWN' THEN BEGIN

    list.xs=list.ys
    list.xd=ev.y
    list.ys=(list.ys+ev.y)/2
    x=[ev.x,ev.x]
    y=[list.ys,ev.y]

    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
        PlotS,[0,list.dynh-1], [y[0],y[0]],/device
        PlotS,[0,list.dynh-1], [y[1],y[1]],/device
        arrow,x[0],y[0],x[1],y[1]
    endif

    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    x=DynTR(x,list.sf,list.hrange[0])
    y=DynTR(y,list.sf,list.vrange[0])
    x=RTStat(x,list.sfh)
    y=RTStat(y,list.sfv)
    hrange=RTStat(list.hrange,list.sfh)
    PlotS, hrange,[y[0],y[0]],/device
    PlotS, hrange,[y[1],y[1]],/device
    arrow,x[0],y[0],x[1],y[1]

    Widget_Control, ev.id, Event_Pro='MR_DRAWARROW'
    Widget_Control, ev.top, Set_UValue=list
    return
endif


x = [ev.x, ev.x]
y = [list.ys, ev.y]

if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    PlotS,[0,list.dynh-1], [y[0],y[0]],/device
    PlotS,[0,list.dynh-1], [y[1],y[1]],/device
    arrow,x[0],y[0],x[1],y[1]
endif

WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
x=DynTR(x,list.sf,list.hrange[0])
y=DynTR(y,list.sf,list.vrange[0])
x=RTStat(x,list.sfh)
y=RTStat(y,list.sfv)
hrange=RTStat(list.hrange,list.sfh)
PlotS, hrange,[y[0],y[0]],/device
PlotS, hrange,[y[1],y[1]],/device
arrow,x[0],y[0],x[1],y[1]

end;pro MR_DRAWSEL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MR_DRAWARROW,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent EQ 'DOWN' THEN BEGIN
    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    endif
    Widget_Control, ev.id, Event_Pro='MR_PROCESS_EVENTS'

    ; Process y coordinates
    ynew = ev.y
    yr = [list.xs,list.xd]

    ynew=DynTR(ynew,list.sf,list.vrange[0])
    yr=DynTR(yr,list.sf,list.vrange[0])

    ; Search for TieI within yr:
    ind=where(((*list.TieI)[1,*] gt (yr[0]<yr[1])) and ((*list.TieI)[1,*] lt (yr[0]>yr[1]))$
    and ((*list.TieI)[0,*] gt list.hrange[0]) and ((*list.TieI)[0,*] lt list.hrange[1]) ,ct)
    if ct eq 0 then return

    keep=(*list.TieI)[1,ind]
    (*list.TieI)[1,ind]=ynew
    (*list.TieO)[*,ind]=GridCoord((*list.TieI)[*,ind],$
           list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix)
    (*list.TieI)[1,ind]=keep
    
    RefreshDisplay, ev.top,list
    return
endif

x = [ev.x, ev.x]
y = [list.ys, ev.y]

if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    PlotS,[0,list.dynh-1], [y[0],y[0]],/device
    PlotS,[0,list.dynh-1], [y[1],y[1]],/device
    arrow,x[0],y[0],x[1],y[1]
endif

WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
x=DynTR(x,list.sf,list.hrange[0])
y=DynTR(y,list.sf,list.vrange[0])
x=RTStat(x,list.sfh)
y=RTStat(y,list.sfv)
hrange=RTStat(list.hrange,list.sfh)
PlotS, hrange,[y[0],y[0]],/device
PlotS, hrange,[y[1],y[1]],/device
arrow,x[0],y[0],x[1],y[1]

end;pro MR_DRAWARROW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ConnectGrid,XYi,GridDim,GridOffset,GridOrigin,GridTilt,GridSpReal,GridSpPix,$
          count=count,ngaps=ngaps,gapcoord=gapcoord

; returns "full connected" output grid couples
; gapcoord: "not fully connected" output grid couples

; ----Grid node indices----
x=(*XYi)[0,*]
y=(*XYi)[1,*]
GridNodeCoordInv,x,y,GridDim,GridSpReal,GridSpPix,GridOffset,GridOrigin,GridTilt
Edind=x+y*GridDim[0]
n=n_elements(Edind)

; ----Grid nodes that are not connected----
temp=bytarr(GridDim[0],GridDim[1])
temp[x,y]=1B
indhole=where(temp eq 0,nhole)

; ----Generate all lines from indb to inde----
indr=indgen(GridDim[0])
indc=indgen(GridDim[1])*GridDim[0]
indb=-1
inde=-1
; Connect between columns
for i=0l,GridDim[0]-2 do begin
    indb=[indb,indc+i]
    inde=[inde,indc+i+1]
endfor
; Connect between rows
for i=0l,GridDim[1]-2 do begin
    indb=[indb,indr+i*GridDim[0]]
    inde=[inde,indr+(i+1)*GridDim[0]]
endfor
count=n_elements(indb)-1
indb=indb[1:count]
inde=inde[1:count]

; ----Test whether all lines are there----
;mag=10
;xb=indb mod GridDim[0]
;yb=indb / GridDim[0]
;xe=inde mod GridDim[0]
;ye=inde / GridDim[0]
;window,0,xsize=mag*GridDim[0],ysize=mag*GridDim[1]
;for i=0l, count-1 do plots,mag*[xb[i],xe[i]],mag*[yb[i],ye[i]],/device
;xb=indhole mod GridDim[0]
;yb=indhole / GridDim[0]
;for i=0l,nhole-1 do plots,xb[i]*mag,yb[i]*mag,psym=2,/device,/continu,symsize=0.5

; ----Search holes----
indh=-1
if nhole ne 0 then begin
    ; indh: all the lines that contain an endpoint which is a "hole"
    for i=0l,nhole-1 do begin
        indh=[indh,where((indb eq indhole[i]) or $
             (inde eq indhole[i]))]
    endfor
    temp=bytarr(count)
    if n_elements(indh) ne 1 then temp[indh[1:*]]=1
    indh=where(temp eq 1,ngaps)

    ; Calculate hole coordinates
    xhb=indb[indh] mod GridDim[0]
    yhb=indb[indh]/GridDim[0]
    xhe=inde[indh] mod GridDim[0]
    yhe=inde[indh]/GridDim[0]
    
    ; Test holes lines
    ;erase
    ;for i=0l,ngaps-1 do plots,[xhb[i],xhe[i]]*mag,[yhb[i],yhe[i]]*mag,/device
    
    GridNodeCoord,xhb,yhb,GridDim,GridSpReal,GridSpPix,GridOffset,GridOrigin,GridTilt
    GridNodeCoord,xhe,yhe,GridDim,GridSpReal,GridSpPix,GridOffset,GridOrigin,GridTilt
    gapcoord=[[xhb],[xhe],[yhb],[yhe]]

    ; Kick out holes
    indb=Shrinkarray(indb,indh,count=count)
    inde=Shrinkarray(inde,indh)
endif else ngaps=0
if count eq 0 then return,[-1,-1]

; Points will not be calculated because we need
; the indices of TieO to do perhaps something
; with the corresponding TieI
; Loop over all present grid points Edind
; to get the wanted indices
Oindb=[0,indb] ; this is done to prevent 2 if statements in for loop
Oinde=[0,inde]
for i=0l,n-1 do begin
    ind=where(indb eq Edind[i])
    Oindb[ind+1]=i
    ind=where(inde eq Edind[i])
    Oinde[ind+1]=i
endfor
Oindb=Oindb[1:count]
Oinde=Oinde[1:count]


return,[[Oindb],[Oinde]]
end;function ConnectGrid
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RecalcPixelSize,list

k=list.GridSpReal/list.GridSpPix ; mm/pixel

if list.GridTilt eq 0 then begin
    list.scrx=k[0]/1000
    list.scry=k[1]/1000
endif else begin
    k=mean(k)/1000
    list.scrx=k
    list.scry=k
endelse

end;pro RecalcPixelSize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ChangeGridParam,newvalue,type,list,fixcon=fixcon

; ----Primary grid properties:----
; GridDim:        Number of grid nodes horizontal and vertical
; GridSpReal:    Spacing between grid points in mm
; GridOffset:    Border offset of the straight grid in the image coordinate system
; GridOrigin:    Origin of the straight or tilted grid in the image coordinate system
; GridTilt:        Tilt of the grid in the image coordinate system
; 
; ----Secondary grid properties:----
; GridSpPix:    Spacing between grid nodes in pixels
; scrx,scry:    Pixel size in m of the straight/tilted grid
; 
; ----Relations:----
; GridOffset = (tiffs-1 - GridSpPix*(GridDim-1))/2
; GridSpPix = GridSpReal/(1000*[scrx,scry])
;
; ----Node indices vs. coordinates in the image coordinate system----
; cfr. GridNodeCoord

if n_elements(fixcon) ne 0 then $
if keyword_set(fixcon) then begin
    ; Prepare for new TieO calculation
    oGridDim=list.GridDim
    oGridOffset=list.GridOffset
    oGridOrigin=list.GridOrigin
    oGridTilt=list.GridTilt
    oGridSpPix=list.GridSpPix
    oGridSpReal=list.GridSpReal
endif

case type of
0:    begin;GridSpReal
    list.GridSpReal=newvalue
    endcase
1:    begin;GridDim
    list.GridDim=newvalue
    tiffs=*list.tiffs
    list.GridSpPix=(tiffs[1:2]-1-2.*list.GridOffset)/(list.GridDim-1.)
    endcase
2:    begin;GridOrigin    
    list.GridOrigin=newvalue
    endcase
3:    begin;GridOffset
    list.GridOffset=newvalue
    tiffs=*list.tiffs
    list.GridSpPix=(tiffs[1:2]-1-2.*list.GridOffset)/(list.GridDim-1.)
    endcase
4:    begin; scrx,scry
    list.scrx=newvalue[0]/10^6. ; um->m
    list.scry=newvalue[1]/10^6.
    tiffs=*list.tiffs
    list.GridSpPix=list.GridSpReal*1000/newvalue ; mm->um
    list.GridOffset=(tiffs[1:2]-1-list.GridSpPix*(list.GridDim-1.))/2.
    endcase
5:    begin;GridTilt
    list.GridTilt=newvalue
    endcase
else:
endcase

RecalcPixelSize,list

if n_elements(fixcon) ne 0 then begin
    if keyword_set(fixcon) then begin
        (*list.TieO)=GridCoordKeep(list.TieO,$
        list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix,$
        oGridDim,oGridOffset,oGridOrigin,oGridTilt,oGridSpReal,oGridSpPix)
    endif else begin
        (*list.TieO)=GridCoord(*list.TieI,$
        list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix)
    endelse
endif

end;pro ChangeGridParam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_RealGrid,ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
if WIDGET_INFO(list.top,/VALID_ID) eq 1 then $
WIDGET_CONTROL, list.top, sensitive=1
end;pro CleanUp_RealGrid
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RealGrid_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=tlist
top=tlist.top
widget_control,top,get_uvalue=list
case widget_info(ev.id,/type) of
1 : begin ;button event
    widget_control,ev.id,get_value=val
    case val of
    'AutoReconnect':begin
           if list.GridDim[0]*list.GridDim[1] ne n_elements(*list.TieO)/2. then return
           (*list.TieO)=AutoGridCoord(list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix)
           AutoSortTie,list.TieI,list.GridDim
           RefreshDisplay,top,list
           
           ID=widget_info(ev.top,FIND_BY_UNAME='TieO')
           widget_control,ID,set_value=[-1,-1]
           return
           endcase
    'Find Double Output':begin
           ID=widget_info(ev.top,FIND_BY_UNAME='TieO')
           dind=DoublePresent(list.TieO,coord=coord)
           widget_control,ID,set_value=coord
           endcase
    'Remove Overlapping Input':begin
           removedoubleTieI,ev,list
           ID=widget_info(ev.top,FIND_BY_UNAME='nGP')
           widget_control,ID,set_value=stringr((*list.nrTie)[0])
           ID=widget_info(ev.top,FIND_BY_UNAME='TieO')
           dind=DoublePresent(list.TieO,coord=coord)
           widget_control,ID,set_value=coord
           RefreshDisplay,top,list
           endcase
    'Move Row':begin
            if widget_info(list.drawdyn,/VALID_ID) then begin
                widget_control,top,sensitive=1
                widget_control,ev.top,sensitive=0
               
                   Widget_Control, list.drawdyn,Event_Pro='MR_PROCESS_EVENTS',draw_button=1
                   Widget_Control, list.drawdyn, set_uvalue=ev.top
               
                   Widget_Control, list.draw,draw_button=0,draw_motion=0
                return
            endif
            endcase
    'Move Column':begin
            if widget_info(list.drawdyn,/VALID_ID) then begin
               widget_control,top,sensitive=1
               widget_control,ev.top,sensitive=0
               Widget_Control, list.drawdyn, Event_Pro='MC_PROCESS_EVENTS',draw_button=1
               Widget_Control, list.drawdyn, set_uvalue=ev.top
               Widget_Control, list.draw,draw_button=0,draw_motion=0
               return
            endif
           endcase
    'View Input Points':begin
           list.ShowTie[0]=ev.select
           RefreshDisplay,top,list
           return
           endcase
    'View Output Points':begin
           list.ShowTie[1]=ev.select
           RefreshDisplay,top,list
           endcase
    'View I-O Connects':begin
           list.ShowTie[2]=ev.select
           RefreshDisplay,top,list
           return
           endcase
    'View O-O Connects':begin
           list.ShowTie[3]=ev.select
           RefreshDisplay,top,list
           return
           endcase
    'Reconnect':begin
           ID=widget_info(ev.top,FIND_BY_UNAME='TieO')
           widget_control,ID,set_uvalue=ev.select
           endcase
    'Delete':begin
           ID=widget_info(ev.top,FIND_BY_UNAME='TieO')
           widget_control,ID,set_uvalue=ev.select ne 1
           endcase
    'Fix connections when calculating output nodes':begin
           tlist.fixcon=ev.select
           widget_control,ev.top,set_uvalue=tlist
           return
           endcase
    'Calculate Output Nodes':begin
           ChangeGridParam,dummy,-1,list,fixcon=tlist.fixcon
           RefreshDisplay,top,list
           return
           endcase
    endcase
    endcase
3:  begin ;text event
    
    widget_control,ev.id,get_uvalue=uval
    widget_control,ev.id,get_value=val
    
    case uval of
    'GridSpReal':begin
            ChangeGridParam,float(strsplit(strcompress(val,/remove_all),',',/extract)),0,list,fixcon=tlist.fixcon
            ID=widget_info(ev.top,FIND_BY_UNAME='PixSize')
            widget_control,ID,set_value=stringr(list.scrx*10^6.)+','+stringr(list.scry*10^6.)
            endcase
    'GridDim':begin
            ChangeGridParam,float(strsplit(strcompress(val,/remove_all),',',/extract)),1,list,fixcon=tlist.fixcon
            ID=widget_info(ev.top,FIND_BY_UNAME='PixSize')
            widget_control,ID,set_value=stringr(list.scrx*10^6.)+','+stringr(list.scry*10^6.)
            ID=widget_info(ev.top,FIND_BY_UNAME='nGPmax')
            widget_control,ID,set_value=stringr(list.GridDim[0]*list.GridDim[1])
            endcase
    'GridOrigin':begin
            ChangeGridParam,float(strsplit(strcompress(val,/remove_all),',',/extract)),2,list,fixcon=tlist.fixcon
            endcase
    'GridOffset':begin
            ChangeGridParam,float(strsplit(strcompress(val,/remove_all),',',/extract)),3,list,fixcon=tlist.fixcon
            ID=widget_info(ev.top,FIND_BY_UNAME='PixSize')
            widget_control,ID,set_value=stringr(list.scrx*10^6.)+','+stringr(list.scry*10^6.)
            endcase
     'PixSize':begin
             ChangeGridParam,float(strsplit(strcompress(val,/remove_all),',',/extract)),4,list,fixcon=tlist.fixcon
             ID=widget_info(ev.top,FIND_BY_UNAME='GridOffset')
            widget_control,ID,set_value=stringr(list.GridOffset[0])+','+stringr(list.GridOffset[1])
            ID=widget_info(ev.top,FIND_BY_UNAME='PixSize')
            widget_control,ID,set_value=stringr(list.scrx*10^6.)+','+stringr(list.scry*10^6.)
             endcase
     'GridTilt':begin
             ChangeGridParam,float(val),5,list,fixcon=tlist.fixcon
             ID=widget_info(ev.top,FIND_BY_UNAME='PixSize')
            widget_control,ID,set_value=stringr(list.scrx*10^6.)+','+stringr(list.scry*10^6.)
             endcase
    endcase
    RefreshDisplay,top,list
    return
    endcase
9:  begin ;tabel event
    if (ev.type eq 4) then begin
       if ev.sel_top ne -1 then begin
         widget_control,ev.id,get_uvalue=uval

         ; ----Indices of Tie points----
         dind=DoublePresent(list.TieO,count=count)
         if (count[0] eq 0) or (ev.sel_top ge n_elements(count)) then return
           if ev.sel_top eq 0 then i=0 else i=total(count[0:ev.sel_top-1])
           dind=dind[i:i+count[ev.sel_top]-1]
           dTieI=(*list.TieI)[*,dind]

         ; ----Reconnect or delete TieI----
         if uval then begin

          list.xs=dTieI[0,0]
          list.ys=dTieI[1,0]
          widget_control,top,sensitive=1
          widget_control,ev.top,sensitive=0
          if widget_info(list.drawdyn,/VALID_ID) then begin
              Widget_Control, list.drawdyn,Event_Pro='MP_PROCESS_EVENTS',$
                     draw_button=1,draw_motion=1
              Widget_Control, list.drawdyn, set_uvalue={apptop:ev.top,dTieI:dTieI,dind:dind,point:0}
              Widget_Control, list.draw,draw_button=0,draw_motion=0
          endif
                  Widget_Control, list.draw,Event_Pro='MP_PROCESS_EVENTS',$
                 draw_button=1,draw_motion=1
                Widget_Control, list.draw, set_uvalue={apptop:ev.top,dTieI:dTieI,dind:dind,point:0}
           endif else begin

             str=stringr(dTieI)
             str='['+str[0,*]+';'+str[1,*]+']'
             sel=ChooseFromList(str,title='Select Tie Points to delete',$
                       nsel=nsel,top=ev.top)
             if nsel eq 0 then return
             ind=dind[sel]

             ;----Delete Point----
             (*list.TieI)=ShrinkArray(*list.TieI,ind,/row)
             (*list.TieO)=ShrinkArray(*list.TieO,ind,/row)

             ;----Find out in which ring----
             s=0l
             i=0l
             repeat begin
                 s+=(*list.nrTie)[i]
                i++
             endrep until (ind[0] lt s)
             Ring=i-1

             (*list.nrTie)[Ring]-=nsel
             if (*list.nrTie)[Ring] eq 0 then begin
                 (*list.nrTie)=ShrinkArray((*list.nrTie),Ring)
                 list.nrTieRing=(list.nrTieRing-1)>0
                 (*list.TieRingD)=ShrinkArray((*list.TieRingD),Ring)
             endif

             ID=widget_info(ev.top,FIND_BY_UNAME='nGP')
             widget_control,ID,set_value=stringr((*list.nrTie)[0])
             ID=widget_info(ev.top,FIND_BY_UNAME='TieO')
             dind=DoublePresent(list.TieO,coord=coord)
             widget_control,ID,set_value=coord
             
             RefreshDisplay,top,list
             return
             endelse
        endif
    endif
    endcase
else:
endcase
widget_control,top,set_uvalue=list

end;pro RealGrid_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotcomp, in

key=in.key
if key eq 0 then begin ; Change Z-range
    list=in.list
    oGroup=list.oGroup
    zr=in.zr

    oChildArr = oGroup->IDL_Container::Get(ISA='IDLgrSurface',/all, COUNT=nKids)
    oChildArr[0]->Setproperty, ZCOORD_CONV=NormCoord(zr,list.PosZ)
    oChildArr = oGroup->IDL_Container::Get(ISA='IDLgrPolyline',/all, COUNT=nKids)
    oChildArr[0]->Setproperty, ZCOORD_CONV=NormCoord(zr,list.PosZ)

endif else begin ; Add Objects
    list=in.list

    ;----plot background----
    oSurface = OBJ_NEW('IDLgrSurface', (list.list2).imgfit,(list.list2).x,(list.list2).y , STYLE=2, SHADING=0, $
               COLOR=[60,60,255], BOTTOM=[64,192,128], $
               XCOORD_CONV=NormCoord((list.list2).xr,list.PosX),$
               YCOORD_CONV=NormCoord((list.list2).yr,list.PosY),$
               ZCOORD_CONV=NormCoord((list.list2).zr,list.PosZ))
    list.oGroup->Add, oSurface


    ;----plot image----
    x=(list.list2).x
    y=(list.list2).y
    axbin=n_elements(x)
    aybin=n_elements(y)
    n=axbin*aybin
    surf=(list.list2).surf
    xmap=(x)#replicate(1,aybin)
    ymap=replicate(1,axbin)#(y)
    xmap=reform(xmap,n)
    ymap=reform(ymap,n)
    oSymbol = OBJ_NEW('IDLgrSymbol',2)
    oPLine = OBJ_NEW('IDLgrPolyline' ,xmap,ymap,reform(surf,n),XCOORD_CONV=NormCoord((list.list2).xr,list.PosX),$
              YCOORD_CONV=NormCoord((list.list2).yr,list.PosY),$
              ZCOORD_CONV=NormCoord((list.list2).zr,list.PosZ),$
              LINESTYLE=6,symbol=oSymbol)
    list.oGroup->Add, oPLine
    list.oHolder->Add, oSymbol


endelse

end;pro plotcomp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro backplot, in

key=in.key
if key eq 0 then begin ; Change Z-range
    list=in.list
    oGroup=list.oGroup
    zr=in.zr

    oChildArr = oGroup->IDL_Container::Get(ISA='IDLgrSurface',/all, COUNT=nKids)
    oChildArr[0]->Setproperty, ZCOORD_CONV=NormCoord(zr,list.PosZ)
    oChildArr[1]->Setproperty, ZCOORD_CONV=NormCoord(zr,list.PosZ)

endif else begin ; Add Objects
    list=in.list
    oGroup=list.oGroup
    ;

    ;----plot background----
    oSurface = OBJ_NEW('IDLgrSurface', (list.list2).imgfit,(list.list2).x,(list.list2).y , STYLE=2, SHADING=0, $
               COLOR=[60,60,255], BOTTOM=[64,192,128], $
               XCOORD_CONV=NormCoord((list.list2).xr,list.PosX),$
               YCOORD_CONV=NormCoord((list.list2).yr,list.PosY),$
               ZCOORD_CONV=NormCoord((list.list2).zr,list.PosZ))
    oGroup->Add, oSurface


    ;----plot image----

    oSurface = OBJ_NEW('IDLgrSurface',(list.list2).surf,(list.list2).x,(list.list2).y ,XCOORD_CONV=NormCoord((list.list2).xr,list.PosX),$
              YCOORD_CONV=NormCoord((list.list2).yr,list.PosY),$
              ZCOORD_CONV=NormCoord((list.list2).zr,list.PosZ),$
              STYLE=2, SHADING=0,COLOR=[60,60,60], BOTTOM=[100,100,100],ALPHA_CHANNEL=0.5)
              ;,STYLE=0,THICK=2)
    oGroup->Add, oSurface
    
;    ,$
endelse

end;pro backplot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotsurf, in

key=in.key
if key eq 0 then begin ; Change Z-range
    list=in.list
    oGroup=list.oGroup
    zr=in.zr

    oChildArr = oGroup->IDL_Container::Get(ISA='IDLgrSurface',/all, COUNT=nKids)
    oChildArr[0]->Setproperty, ZCOORD_CONV=NormCoord(zr,list.PosZ)

endif else begin ; Add Objects
    list=in.list
    oGroup=list.oGroup
    ;
    xx=indgen((list.list2).xr[1]-(list.list2).xr[0]+1)+(list.list2).xr[0]
    yy=indgen((list.list2).yr[1]-(list.list2).yr[0]+1)+(list.list2).yr[0]

    ; ----Create the surface----
    oSurface = OBJ_NEW('IDLgrSurface', (list.list2).tiff,xx,yy , STYLE=2, SHADING=0, $
               COLOR=[60,60,255], BOTTOM=[64,192,128], $
               XCOORD_CONV=NormCoord((list.list2).xr,list.PosX),$
               YCOORD_CONV=NormCoord((list.list2).yr,list.PosY),$
               ZCOORD_CONV=NormCoord((list.list2).zr,list.PosZ))
    oGroup->Add, oSurface
endelse

end;pro plotsurf
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro red_res,chisqr,areas2,positions2,sx2,n2,sigma
; save as:      sigma of area must be finite
;          area > sigma
;          area > 0
n=n_elements(areas2)
if n eq 1 and areas2[0] eq 0 then return

; finite?
mask=finite(sigma)
ind=where(mask ne 1,count)
if count ne 0 then begin
ind=ind / 4
areas2=ShrinkArray(areas2,ind)
;positions2=ShrinkArray(positions2,ind)
;sx2=ShrinkArray(sx2,ind)
;n2=ShrinkArray(n2,ind)
sigma=ShrinkArray(sigma,ind,/row)
n=n_elements(areas2)
if n eq 1 and areas2[0] eq 0 then return
endif

; areas > sigma and > 0?
ind=where((areas2 gt reform(sigma[0,*])) and $
         areas2 gt 0,count)
if count eq 0 then begin
    areas2=[0.]
    return
endif
areas2=areas2[ind]
;positions2=positions2[ind]
;sx2=sx2[ind]
;n2=n2[ind]
sigma=sigma[*,ind]

if chisqr gt 1000 then begin
areas2=[0.]
return
endif

end;pro red_res
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_AutoLoad, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
ptr_free,list.PTR_files
end;pro CleanUp_AutoLoad
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AutoLoad_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list

case widget_info(ev.id,/type) of
0:  begin;Timer event (base widget)
    ; ----Stop timer?----
    if list.binit then begin
        list.binit=0b
        widget_control,ev.top,set_uvalue=list
        result=DetectFiles2(list,ev)
    endif else begin
        if list.timer eq 0 then return
        ;result=DetectFiles(list,ev)
        result=DetectFiles2(list,ev)
    endelse
    endcase
1:  begin;Button event
    widget_control,ev.id,get_value=val
    case val of
    'Update Main Window':begin
                 list.update=ev.select
                 widget_control,ev.top,set_uvalue=list
                 endcase
    'Run timer':   begin
              list.timer=ev.select
              widget_control,ev.top,set_uvalue=list
              if list.timer then WIDGET_CONTROL, ev.top, TIMER=list.sec
              endcase
    'Stop':    begin
          widget_control,ev.top,/destroy
          endcase
    endcase
    endcase
3:    begin; text
    widget_control,ev.id,get_value=val
    d=0.
    reads,val,d
    list.sec=(-1)>d
    widget_control,ev.top,set_uvalue=list
    ; Make the current timer time out
    WIDGET_CONTROL, ev.top, TIMER=0
    endcase
6:  begin;list
    MainUpdate,list,ev
    widget_control,ev.id,/input_focus
    endcase
endcase
end;pro AutoLoad_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_SelMask, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
if WIDGET_INFO(list.ID1,/VALID_ID) ne 0 then widget_control,listtemp.ID1,sensitive=1
end;pro CleanUp_SelMask
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SelMask_event , ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control, ev.top, get_uvalue=listtemp
widget_control,ev.id,get_value=val
case val of
'Select All':begin
         (*listtemp.select)[*]=ev.select
         nbutton=n_elements(*listtemp.select)
         uname=stringr(indgen(nbutton))
         for i=0l,nbutton-1 do begin
          ID=widget_info(ev.top,FIND_BY_UNAME=uname[i])
          widget_control,ID,set_button=ev.select
         endfor
         widget_control,ev.top,set_uvalue=listtemp
         endcase
    'OK' : begin
           WIDGET_CONTROL, ev.top, /DESTROY
           endcase
  else:     begin
         i=1L
         widget_control,ev.id,get_uvalue=uval
         (*listtemp.select)[uval]=ev.select
         widget_control,ev.top,set_uvalue=listtemp
         endcase
endcase
end;pro SelMask_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function eigenns,A,EIGENVECTORS=evecs,DOUBLE=DOUBLE,RESIDUAL=RESIDUAL

hes = ELMHES(A)    ;Reduce to upper Hessenberg format.
evals = HQR(hes)    ;Compute the eigenvalues.

; Note that complex eigenvalues of an n x n real,
; nonsymmetric array always occur in complex conjugate pairs.

;Compute the eigenvectors and the residual for each eigenvalue/eigenvector pair.
evecs = EIGENVEC(A, evals, DOUBLE=DOUBLE, RESIDUAL=RESIDUAL); each row is an eigenvector
evecs=float(evecs); real part

return, float(evals); real part
end;function eigenns
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ConicSectionLoci,par,azistep=azistep
; Cfr. ConicTEllipse

; AA x^2 + BB x y + CC y^2 + DD x + EE y + FF = 0
; Auu u^2 + Avv v^2 + Au u + Av v + A0 = 0 (rotated)
; Auu s^2 + Avv t^2 + A1 = 0 (centered)

if (par[2] eq par[0]) then thetarad=0. $
else thetarad = 0.5*atan(par[1]/(par[2] - par[0]))

cost = cos(thetarad)
sint = sin(thetarad)
sin_squared = sint*sint
cos_squared = cost*cost
cos_sin = sint * cost

A0 = par[5]
Au =   par[3] * cost - par[4] * sint
Av =   par[3] * sint + par[4] * cost
Auu = par[0] * cos_squared + par[2] * sin_squared - par[1] * cos_sin
Avv = par[0] * sin_squared + par[2] * cos_squared + par[1] * cos_sin
A1 = A0 - Au^2/(4*Auu) - Av^2/(4*Avv)

uCentre = - Au/(2*Auu)
vCentre = - Av/(2*Avv)

b=0
e=2*!dpi
if ~keyword_set(azistep) then azistep=!pi/180
azi=azistep*findgen(ceil(2*!dpi/azistep)+1)
cosazi = cos(azi)
sinazi = sin(azi)
R = sqrt(-A1/(Auu*cosazi^2+Avv*sinazi^2))

return, [[cost,sint],[-sint,cost]]##[[R*cosazi + uCentre],[R*sinazi + vCentre]]

end;function ConicSectionLoci
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ConicTEllipse,par,Xin,Yin,parrot=parrot

; Cte.[AA x^2 + BB x y + CC y^2 + DD x + EE y + FF] = 0

; Rotation angle
if (par[2] eq par[0]) then thetarad=0. $ ; this should be 180 or -180, but because this means a circle => take 0
else thetarad = 0.5*atan(par[1]/(par[2] - par[0])) ; [-45,45]

repeat begin
    cost = cos(thetarad)
    sint = sin(thetarad)
    sin_squared = sint*sint
    cos_squared = cost*cost
    cos_sin = sint * cost
    
    ; **** Rotated XY frame to UV frame with thetarad: ****
    ; u = x * cost - y * sint
    ; v = x * sint + y * cost
    ; Auu u^2 + Avv v^2 + Au u + Av v + A0 = 0
    A0 = par[5]
    Au =   par[3] * cost - par[4] * sint
    Av =   par[3] * sint + par[4] * cost
    Auu = par[0] * cos_squared + par[2] * sin_squared - par[1] * cos_sin
    Avv = par[0] * sin_squared + par[2] * cos_squared + par[1] * cos_sin
    ;Auv = (par[0]-par[2]) * 2*cos_sin + par[1]*(cos_squared-sin_squared)
    
    ; Ellipse center
    uCentre = - Au/(2*Auu)
    vCentre = - Av/(2*Avv)
    wCentre = Auu*uCentre*uCentre + Avv*vCentre*vCentre - A0
    
    ; Geometric radii
    Rx = abs(wCentre/Auu)
    Ry = abs(wCentre/Avv)

    bstop=Ry ge Rx
    if ~bstop then $
        if (thetarad lt 0) then thetarad+=!dpi/2 else thetarad-=!dpi/2
endrep until bstop

parrot=[Auu,0,Avv,Au,Av,A0]

; Rotate vector [uCentre,vCentre] with b to [xCentre,yCentre]
xCentre = uCentre * cost + vCentre * sint
yCentre = -uCentre * sint + vCentre * cost

; Geometric radii
Rx = sqrt(Rx)
Ry = sqrt(Ry)

;xt=Xin*cost-Yin*sint
;yt=Xin*sint+Yin*cost
;loadctrev,-39
;window,1
;m=abs(max(xt,/abs))>abs(max(yt,/abs))
;plot,[-1.2,1.2]*m,[-1.2,1.2]*m,/xs,/ys,/nodata,/iso
;plots,xt,yt
;plots,Xin,Yin,psym=3

return, [xCentre, yCentre, Rx, Ry, thetarad]; xcenel,ycenel,h,H,Beta
end;function ConicTEllipse
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitellipse,Xin,Yin,parraw=par,parrot=parrot

; **** normalize data ****
; A xn^2 + B xn yn + C yn^2 + D xn + E yn + F = 0
; coeff=> A
mx = mean(Xin)
my = mean(Yin)
sx = (max(Xin)-min(Xin))/2
sy = (max(Yin)-min(Yin))/2

x = (Xin-mx)/sx
y = (Yin-my)/sy

; force into rows
x=reform(x)
y=reform(y)

; Build transposed design matrix
DT = [ [x*x],[x*y],[y*y],[x],[y],[replicate(1,n_elements(Xin))] ]

; Build scatter matrix
S = DT ## transpose(DT)

; Build 6x6 constraint matrix
C=intarr(6,6)
C[2,0] = -2
C[1,1] = 1
C[0,2] = -2

; **** Solve eigensystem: S.A=L.C.A ****
; Break into blocks
tmpA = S[0:2,0:2]
tmpB = S[3:5,0:2]
tmpC = S[3:5,3:5]
tmpD = C[0:2,0:2]
tmpE = invert(tmpC) ## transpose(tmpB)
tmp=invert(tmpD) ## (tmpA - tmpB##tmpE)
eval_x=eigenns(tmp,EIGENVECTORS=evec_x)

; Find the positive (as det(tmpD) < 0) eigenvalue
I = where(((eval_x lt 1e-8) and finite(eval_x)) eq 1)
I>=0

; Extract eigenvector corresponding to negative eigenvalue
A = transpose(evec_x[*,I[0]])

; Recover the bottom half...
evec_y = (-tmpE) ## A
A = [[A], [evec_y]]

; **** unnormalize ****
; replace in conic equation xn by (x-mx)/sx and yn by (y-my)/sy
; AA x^2 + BB x y + CC y^2 + DD x + EE y + FF = 0
; coeff=> param
AA=A[0]*sy*sy
BB=A[1]*sx*sy
CC=A[2]*sx*sx
DD=-2*AA*mx - BB*my + A[3]*sx*sy*sy
EE=-BB*mx - 2*CC*my + A[4]*sx*sx*sy
FF=AA*mx*mx + BB*mx*my + CC*my*my - A[3]*sx*sy*sy*mx - A[4]*sx*sx*sy*my + A[5]*sx*sx*sy*sy
; par=param*sx*sx*sy*sy
par = [AA,BB,CC,DD,EE,FF]

; **** convert to elliptical coordinates ****
return,ConicTEllipse(par,Xin,Yin,parrot=parrot)

end;function fitellipse
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RefTilt,M,nj,x,select,tt,scrx,scry,NEGATIVE=NEGATIVE

; ----Find a,b,xc,yc,dist for each ring----
indselect=where(select,ct)
a=dblarr(ct,2)
b=a
xc=a
yc=a
dist=a

ind=0l
off=0l
for j=0l,M-1 do begin; loop over all rings
    if select[j] eq 1 then begin; ring selected?

        ; Fit elliptical param.
        params=fitellipse(scrx*x[ind:ind+nj[j]-1,0],scry*x[ind:ind+nj[j]-1,1]);[xCentre, yCentre, Rx, Ry, thetarad]
        if params[2] gt params[3] then begin
            ; This should never happen...
            Hmajor = params[2]
            Hminor = params[3]
        endif else begin
            Hmajor = params[3]
            Hminor = params[2]
        endelse

        ; alpha
        a2=abs(asin(sqrt(1-Hminor^2/Hmajor^2)*cos(tt[j])))

        ; tz
        di = Hmajor*abs((cos(a2))^2-(sin(a2))^2*(tan(tt[j]))^2)/(cos(a2)*tan(tt[j]))

        ; c
        ytroot=di*tan(tt[j])/(sin(a2)*abs(tan(tt[j]))+cos(a2))
        c=Hmajor-ytroot

        ; Save params
        c_=c
        if a2*di lt 0 then c_=-c_
        xc[off,0]=params[0]+c_*sin(params[4])
        yc[off,0]=params[1]+c_*cos(params[4])
        a[off,0]=a2
        b[off,0]=params[4]
        dist[off,0]=di
        
        ; Alternative
        a2=-a2
        c_=c
        if a2*di lt 0 then c_=-c_
        
        xc[off,1]=params[0]+c_*sin(params[4])
        yc[off,1]=params[1]+c_*cos(params[4])
        a[off,1]=a2
        b[off,1]=params[4]
        dist[off,1]=di

        off++
    endif
    ind+=nj[j]
endfor

if ct eq 1 then begin
    if keyword_set(NEGATIVE) then i=1 else i=0
endif else begin
    i=0
    xs1=stdev(xc[*,0])
    xs2=stdev(xc[*,1])
    ys1=stdev(yc[*,0])
    ys2=stdev(yc[*,1])
    
    xs=xs1>xs2
    ys=ys1>ys2
    
    if xs gt ys then i=xs1 gt xs2 $
    else i=ys1 gt ys2
endelse

return,{a:total(a[*,i])/ct,b:total(b[*,i])/ct,xc:total(xc[*,i])/(ct*scrx),yc:total(yc[*,i])/(ct*scry),dist:total(dist[*,i])/ct}
end;function RefTilt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SortPoints,points,nj,center
M=n_elements(nj)
ind=0
for j=0l,M-1 do begin
    x=points[0,ind:ind+nj[j]-1]
    y=points[1,ind:ind+nj[j]-1]
    azi=atan(y-center[1],x-center[0])
    indazi=where(azi lt 0,count)
    if count ne 0 then azi[indazi]=azi[indazi]+2.*!dpi
    indsort=sort(azi)
    points[0,ind:ind+nj[j]-1]=x[indsort]
    points[1,ind:ind+nj[j]-1]=y[indsort]
    ind=ind+nj
endfor
end;function SortPoints
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FitTilt,tlist,evtop
; ----Pushed on 'Go'----
widget_control,tlist.top,get_uvalue=list

; ----refbool gives refinement flags----
;refbutton[0]: Refine tilt angle?
;refbutton[1]: Refine rotation angle?
;refbutton[2]: Refine center?
;refbutton[3]: Refine Distance?
;refbutton[4]: Refine two-theta?
;refbutton[5]: Refine Wavelength?
;refbutton[6]: Refine D-spacings?
;refbutton[7]: Use intensities as weights?
;refbutton[8]: Reject outliers?
;refbutton[9]: Fit with constraints?
;refbutton[10]: Plot?
if total(tlist.refbool[0:6]) eq 0 then begin
    printw,list.OutID,'Nothing to refine'
    return
endif

tmp=where((*tlist.TiltRingD)[0:list.nrTiltRing-1] eq 0)
if (tmp[0] ne -1) then begin
    printw,list.OutID,'D-spacing 0 is not valid'
    return
endif

; ----Set fit mode----
if tlist.refbool[4] then fitmode=4 else $
fitmode=tlist.refbool[5]*2+tlist.refbool[6]

; ----Start rejection-loop----
; ----Take rings----
refit:M=list.nrTiltRing
nj=(*list.nrTilt)[0:M-1]
ntotal=total(nj)
TiltRingD=(*tlist.TiltRingD)[0:M-1]

; ----Remove non-selected rings----
indsel=where(tlist.Nrings eq 1,count) ;indices of rings to keep
if count eq 0 then return
genind=-1
for i=0l,count-1 do begin
    if indsel[i] eq 0 then offset=0 else offset=total(nj[0:indsel[i]-1])
    genind=[genind,offset+indgen(nj[indsel[i]])]
endfor
genind=genind[1:*] ;indices of points to keep (i.e. points on selected rings)
nj=nj[indsel]
TiltRingD=TiltRingD[indsel]
x=(*list.TiltI)[*,genind]
ntotal=total(nj)
M=count

; ----Sort points in rings by azimuth----
ind=0l
for i=0l,M-1 do begin
    points=x[*,ind:ind+nj[i]-1]
    SortPoints,points,nj[i],tlist.center
    x[*,ind:ind+nj[i]-1]=points
    ind=ind+nj[i]
endfor
(*list.TiltI)[*,genind]=x
x=transpose(float(x)) ; in pixels

; ----Initialize first five parameters----
; Remark: the sign of the tilt and rotation angle
; used in the fitfunction are the opposite of the real ones
A=-1
Afix=-1
constr=[-1,-1]
if tlist.refbool[0] then begin
    A=[A,+tlist.a]
    constr=[[constr],[-!dpi/2,!dpi/2]]
endif else Afix=[Afix,+tlist.a]
if tlist.refbool[1] then begin
    A=[A,+tlist.b]
    constr=[[constr],[-!dpi/2,!dpi/2]]
endif else Afix=[Afix,+tlist.b]
if tlist.refbool[2] then begin
     A=[A,tlist.center*[list.scrx,list.scry]]
    s=*list.tiffs
    constr=[[constr],[[0,s[1]-1]*list.scrx],[[0,s[2]-1]*list.scry]]
endif else Afix=[Afix,tlist.center*[list.scrx,list.scry]]
di=tlist.dist
if tlist.refbool[3] then begin
    A=[A,di]
    if di lt 0 then constr=[[constr],[2*di,0]] $
    else constr=[[constr],[0,2*di]]
endif else Afix=[Afix,di]


; ----Values of fit function in x----
case fitmode of

0:    begin;d and lambda fixed (i.e. two-theta fixed)
    ; ----y----
    y=fltarr(ntotal)
    tt=BraggDtoX(TiltRingD[0:M-1],tlist.lambda,/onlyt)
    ind=0l
    for j=0l,M-1 do begin
        y[ind:ind+nj[j]-1]=tt[j]
        ind+=nj[j]
    endfor
    endcase

1:    begin;refine d-spacings, lambda fixed
    ; ----y----
    y=replicate(tlist.lambda/2.,ntotal)
    ; ----param----
    A=[A,TiltRingD[0:M-1]]
    constr=[[constr],[replicate(1,1,M),reform(2*TiltRingD[0:M-1],1,M)]]
    Afix=[Afix,tlist.lambda]
    endcase

2:    begin;refine wavelength, d fixed
    ; ----y----
    y=fltarr(ntotal)
    ind=0l
    for j=0l,M-1 do begin
        y[ind:ind+nj[j]-1]=2*TiltRingD[j]
        ind=ind+nj[j]
    endfor
    ; ----param----
    A=[A,tlist.lambda]
    constr=[[constr],[0,2*tlist.lambda]]
    Afix=[Afix,TiltRingD[0:M-1]]
    endcase

3:    begin;refine d and lambda
    ; ----y----
    y=replicate(0.5,ntotal)
    ; ----param----
    A=[A,tlist.lambda]
    constr=[[constr],[tlist.lambda,tlist.lambda]]
    A=[A,TiltRingD[0:M-1]]
    constr=[[constr],[replicate(1,1,M),reform(2*TiltRingD[0:M-1],1,M)]]
    endcase

4:    begin;refine two-theta
    ; ----y----
    y=fltarr(ntotal)
    ; ----param----
    tt=fltarr(M)
    ttexp=BraggXYtoX(x[*,0],x[*,1],tlist.lambda,tlist.dist,tlist.scrx,tlist.scry,tlist.a,tlist.b,tlist.phipixphi0,tlist.center,/onlyt)
    ind=0l
    for j=0l,M-1 do begin
       if tlist.refbool[2] then $
         tt[j]=total(ttexp[ind:ind+nj[j]-1])/nj[j] $
       else tt[j]=min(ttexp[ind:ind+nj[j]-1])
       ind+=nj[j]
    endfor
    A=[A,tt]
    constr=[[constr],[replicate(0.,1,M),replicate(!dpi/2.,1,M)]]
    endcase
endcase

; ----Prepare fit parameters-----
A=A[1:*]
if n_elements(Afix) ne 1 then Afix=Afix[1:*]
constr=constr[*,1:*]
struc={nj:nj,M:M,Afix:Afix,refbool:tlist.refbool}
func='gfuncTilt'+stringr(fitmode)
A0=A
nparam=n_elements(A)
if tlist.refbool[7] then weights=(*list.tiff)[x[*,0],x[*,1]] else weights=fltarr(ntotal)+1

; ----Constraints (physical boundaries)----
; a,b => [-pi/2,pi/2]
; center => in image
; tt => [-pi/2,pi/2]
; di => [0,2*di]
; lambda => [0,2*lambda]
; d => [0,2*d]
if tlist.refbool[9] then begin
    Ainfo=MakeAinfo(n_elements(A))
    Ainfo.limited=1
    Ainfo.limits=constr
endif else constr[*]=0

; ----Fit----
widget_control, /hourglass
xfit=x
xfit[*,0]*=list.scrx
xfit[*,1]*=list.scry

A=double(A)
yfit=NLLS(xfit,y, weights=weights, A, sigma=sigma,func,$
    CHISQ=CHISQ,ITER=ITER,nfree=nfree,struc=struc,$
                 errormsg=errormsg,CError=CError,Ainfo=Ainfo);,/numder,/debug)

if CError lt 0 then begin
    temp=dialog_message(ErrorMsg,/error)
    return
endif

Areject=A

; ----Construct fit output----
ID=widget_info(evtop,FIND_BY_UNAME='wrapbase')
if tlist.outid ne 0 then widget_control,tlist.outid,/destroy
rowlabels=['a(tilt'+msymbols('degrees')+')','b(rotation'+msymbols('degrees')+')','Xcenter (pixels)','Ycenter (pixels)',$
'Distance (cm)','Lambda ('+msymbols('angstroms')+')']
if tlist.refbool[4] then rowlabels=[rowlabels,'Two-Theta ('+msymbols('degrees')+')'] else rowlabels=[rowlabels,'D-spacing ('+msymbols('angstroms')+')']
collabels=['initial','lowest','highest','refined','SD for algorithm']
idlabel=widget_info(evtop,FIND_BY_UNAME='labelout')
tablevalue=fltarr(5,6+M)
col=replicate(255b,3,5,6+M)


; ----Display and save results----
if tlist.refbool[7] then pref='Int. weights, Iter:' else pref='No weights, Iter:'
widget_control,idlabel,set_value=pref+stringr(iter)+', Reduced chisq:'+stringr(chisq)+', datapoints:'+stringr(long(ntotal))
offset=0
if tlist.refbool[0] then begin
    tlist.a=+A[offset]
    tablevalue[*,0]=([A0[offset],constr[*,offset],A[offset],sigma[offset]])*(180./!dpi)
    offset++
endif else begin
    tablevalue[0,0]=+tlist.a*180./!dpi
    col[*,*,0]=200
endelse

if tlist.refbool[1] then begin
    tlist.b=+A[offset] mod 360
    tablevalue[*,1]=([A0[offset],constr[*,offset],A[offset],sigma[offset]])*(180./!dpi)
    offset++
endif else begin
    tablevalue[0,1]=+tlist.b*180./!dpi
    col[*,*,1]=200
endelse

if tlist.refbool[2] then begin
    tlist.center=A[offset:offset+1]/[list.scrx,list.scry]
    tablevalue[*,2]=[A0[offset],constr[*,offset],A[offset],sigma[offset]]/list.scrx
    offset++
    tablevalue[*,3]=[A0[offset],constr[*,offset],A[offset],sigma[offset]]/list.scry
    offset++
endif else begin
    tablevalue[0,2:3]=transpose(tlist.center)
    col[*,*,2:3]=200
endelse

if tlist.refbool[3] then begin
    tlist.dist=A[offset]
    tablevalue[*,4]=[A0[offset],constr[*,offset],A[offset],sigma[offset]]*100
    offset++
endif else begin
    tablevalue[0,4]=tlist.dist*100
    col[*,*,4]=200
endelse

if tlist.refbool[5] then begin
    tlist.lambda=A[offset]
    tablevalue[*,5]=[A0[offset],constr[*,offset],A[offset],sigma[offset]]
    offset++
endif else begin
    tablevalue[0,5]=tlist.lambda
    col[*,*,5]=200
endelse

if tlist.refbool[4] then begin
    tablevalue[*,6:M+5]=transpose([[A0[offset:offset+M-1]],[transpose(constr[*,offset:offset+M-1])],$
                    [A[offset:offset+M-1]],[sigma[offset:offset+M-1]]])*(180./!dpi)
    TiltRingD=BraggTtoX(A[offset:offset+M-1],tlist.lambda,/onlyd) ; this is just for output
endif else begin
    if tlist.refbool[6] then begin
        (*tlist.TiltRingD)[indsel]=A[offset:offset+M-1]
        tablevalue[*,6:M+5]=transpose([[A0[offset:offset+M-1]],[transpose(constr[*,offset:offset+M-1])],$
                        [A[offset:offset+M-1]],[sigma[offset:offset+M-1]]])
    endif else begin
        tablevalue[0,6:M+5]=transpose(TiltRingD[0:M-1])
        col[*,*,6:M+5]=200
    endelse
endelse

if ~tlist.refbool[9] then col[*,1:2,*]=200

table=widget_table(ID,COLUMN_LABELS=collabels,ROW_LABELS=rowlabels,BACKGROUND_COLOR=col,$
    xsize=5,ysize=M+6,/RESIZEABLE_COLUMNS,/RESIZEABLE_ROWS,uname='table',value=tablevalue)
tlist.outid=widget_info(evtop,FIND_BY_UNAME='table')

; ----Rejection and fit quality evaluation----
; Fit quality: compare r(azimuth) of
;            1. Fit: two-theta circles --(rotate and tilt parameter eq.)--> A
;            2. Data: (x,y) --(in polar coord.)--> B
; Residual=B-A
; Remark: Residual=f-cte

; Rejection: (x,y) --(gfuncTilt)--> f (d-spacing or two-theta)
; f must in <f>msymbols('plus_minus')k.SD(f)
; f is a constant when good fit (except for noise)!
; In the actual plot: f - <f>

LoadctRev,0,/silent,rev=-1

; ----Reject outliers and re-fit----
if (tlist.refbool[8] eq 1) then begin
    ; ----Find outliers method----
    struc.refbool[9]=0 ;Areject are physical parameters
    if (fitmode eq 0) or (fitmode eq 4) then begin
        gfuncTilt0,xfit, Areject,f,pder,struc=struc
        ; f = 2.theta(rad)
        f*=180./!dpi
        ytitle='f - <f>('+msymbols('degrees')+')'
        legend=[msymbols('plus_minus')+'(k.SD(f))','(f - <f>)','f = Calculated two-theta ']
    endif else begin
        gfuncTilt2,xfit, Areject,f,struc=struc
        ; f = 2.d
        f/=2.
        ytitle='f - <f> ('+msymbols('angstroms')+')'
        legend=[msymbols('plus_minus')+'(k.SD(f))','(f - <f>)','f = Calculated d-spacing']
    endelse
    f0=f*0
    sdf0=f0
    ind=0l
    for j=0l,M-1 do begin
           f0[ind:ind+nj[j]-1]=total(f[ind:ind+nj[j]-1])/nj[j] ; mean
           sdf0[ind:ind+nj[j]-1]=sqrt(total((f[ind:ind+nj[j]-1]-f0[ind])^2.)/(nj[j]-1)) ; standard deviation
        ind+=nj[j]
    endfor
    sdf0=tlist.k*sdf0 ; k.standard deviation
    f1=sdf0
    f2=-sdf0
    f-=f0
    indrej=where((f gt f1) or (f lt f2),countoutl)

    ; ----Visualisation of rejecting----
    if tlist.refbool[10] then begin
        DEVICE, WINDOW_STATE=wstate
        if tlist.Index[0] ne -1 then $
            if wstate[tlist.Index[0]] ne 0 then wdelete,tlist.Index[0]
        tlist.Index[0]=winid()
        window,tlist.Index[0]

        ymax=max([f1,f])
        ymin=min([f2,f])
        plot,f1,title='Outliers',xtitle='DataPoints',ytitle=ytitle,$
        yrange=[ymin,ymax],linestyle=2,xstyle=1
        oplot,f2,linestyle=2
        oplot,f
        LoadctRev,3,/silent
        if countoutl ne 0 then begin
            tvlct,RoldCT,GoldCT,BoldCT,/get
            plots,indrej,f[indrej],psym=2,color=16.*7
            tvlct,RoldCT,GoldCT,BoldCT
        endif
        plots,indrej,f[indrej],psym=2,color=16.*7
        
;        p=plot(f1,title='Outliers',xtitle='DataPoints',ytitle=ytitle,$
;            yrange=[ymin,ymax],linestyle=2,xstyle=1,xrange=[0,2000])
;        p1=plot(f2,linestyle=2,/overplot)
;        p2=plot(f,/overplot)
;        p3=plot(indrej,f[indrej],symbol=2,color='ff0000'xl,/overplot,linestyle=6)
    
        xyouts,0.2,0.02,'____ '+legend[1],/normal
        xyouts,0.35,0.02,'- - - '+legend[0],/normal
        xyouts,0.01,0.02,legend[2],/normal
    endif

    ; ----Delete outliers----
    if (countoutl gt ntotal-5*M) then $
    result=DIALOG_MESSAGE('To many outliers. Operation will stop and keep points.',/information) $
    else if (countoutl ne 0) then begin
        Result = DIALOG_MESSAGE(stringr(countoutl)+' outliers found. Reject and refit?',/question)
        if Result eq 'Yes' then begin
            sum=0L
            for j=0l,M-1 do begin
               ind=where((indrej ge sum) and (indrej lt sum+nj[j]),count)
               sum+=nj[j]
               nj[j]-=count
            endfor
            ind=where(nj eq 0,count)
            if (count ne 0) then begin
                result=DIALOG_MESSAGE('Rings are deleted. Operation will stop and keep points.',/information)
            endif else begin
                (*list.nrTilt)[indsel]=nj
                (*list.TiltI)=ShrinkArray((*list.TiltI),genind[indrej],/ROW)
                goto,refit
            endelse
        endif
    endif; else No outliers found
endif

; ----Visualisation of fit----
if tlist.refbool[10] then begin
    ; ----Calculate Radius vs azimuth for each ring: observed and fitted----
    case fitmode of
        0:ttdisplay=tt
        1:ttdisplay=BraggDtoX(A[offset:offset+M-1],tlist.lambda,/onlyt)
        2:ttdisplay=BraggDtoX(TiltRingD[0:M-1],tlist.lambda,/onlyt)
        3:ttdisplay=BraggDtoX(A[offset:offset+M-1],tlist.lambda,/onlyt)
        4:ttdisplay=A[offset:offset+M-1]
    endcase

    phiDet=azDetRange(0,2*!pi,tlist.phipixphi0,tlist.scrx,tlist.scry,tlist.a,tlist.b)
    nphi=n_elements(phiDet)
    phiexp=fltarr(ntotal)
    rexp=phiexp
    phifit=fltarr(M*nphi)
    rfit=phifit
    ind1=0L
    ind2=0L
    for j=0l,M-1 do begin
        ; For this ring: calculate the coordinates of the conic section
        ; with the angles, lambda, apex,...
        xEllipse=EllipsePlotCoord(ttdisplay[j],phiDet,tlist.a,tlist.b,tlist.phipixphi0,tlist.dist,tlist.scrx,tlist.scry,tlist.center,nvalid=nvalid)
        if nvalid ne 0 then begin
            xout=xEllipse[*,0]-tlist.center[0]
            yout=xEllipse[*,1]-tlist.center[1]
    
            ; From the coordinates of the conic section we calculate
            ; their radius and azimuth
            r=sqrt(xout*xout+yout*yout)
            phi=atan(yout,xout)
            indtemp=where(phi lt 0,count)
            if count ne 0 then phi[indtemp]+=2*!dpi
            indtemp=sort(phi)
            phi=phi[indtemp]
            r=r[indtemp]
            nphi2=n_elements(phi)
            phifit[ind1:ind1+nphi2-1]=phi
            rfit[ind1:ind1+nphi2-1]=r
    
            ; Observed radius and azimuth of this ring
            xout=x[ind2:ind2+nj[j]-1,0]-tlist.center[0]
            yout=x[ind2:ind2+nj[j]-1,1]-tlist.center[1]
            r=sqrt(xout*xout+yout*yout)
            phi=atan(yout,xout)
            indtemp=where(phi lt 0,count)
            if count ne 0 then phi[indtemp]+=2*!dpi
            indtemp=sort(phi)
            phi=phi[indtemp]
            r=r[indtemp]
            phiexp[ind2:ind2+nj[j]-1]=phi
            rexp[ind2:ind2+nj[j]-1]=r
        endif

        ind1+=nphi
        ind2+=nj[j]
    endfor

    ; Get standard window size
    i=winid()
    window,i,/pixmap
    xs=!D.x_size
    ys=!D.y_size
    wdelete,i
    !P.MULTI = [0, 1, M]

    ; Make draw widget
    if widget_info(tlist.Index[1],/VALID_ID) eq 1 then WIDGET_CONTROL, tlist.Index[1], /destroy

    tlist.Index[1]=widget_base(title='Fit result')
    draw = WIDGET_DRAW(tlist.Index[1], XSIZE=xs, YSIZE=M*ys,/scroll,$
            x_scroll_size=xs, y_scroll_size=ys)
    WIDGET_CONTROL, tlist.Index[1], /REALIZE

    ; ----Plot Radius vs azimuth: observed and fitted----
    ind1=0L
    ind2=0L

    
    for j=0l,M-1 do begin
        plot,phiexp[ind2:ind2+nj[j]-1],rexp[ind2:ind2+nj[j]-1],$
        title='Ring '+stringr(j),xtitle='Azimuthal (rad)',$
        ytitle='Radius (pixels)',psym=2,ystyle=1,xstyle=1,charsize=M<2
        tvlct,RoldCT,GoldCT,BoldCT,/get
        LoadctRev,3,/silent
        oplot,phifit[ind1:ind1+nphi-1],rfit[ind1:ind1+nphi-1],color=16.*10,thick=2
        tvlct,RoldCT,GoldCT,BoldCT

;        if j eq 8 then begin
;             p=plot(phiexp[ind2:ind2+nj[j]-1],rexp[ind2:ind2+nj[j]-1],$
;                title='Ring '+stringr(j),xtitle='Azimuthal (rad)',$
;                ytitle='Radius (pixels)',symbol=2,ystyle=1,xstyle=1,linestyle=6)
;            p1=plot(phifit[ind1:ind1+nphi-1],rfit[ind1:ind1+nphi-1],color='ff0000'xl,thick=2,/overplot)
;        endif
        
        ind1=ind1+nphi
        ind2=ind2+nj[j]
    endfor
    
    ;WRITE_BMP, 'C:\Misc\LecPresPub\EXRS2006\Calibration.bmp', tvrd(true=1),/RGB

    ; Make draw widget
    if widget_info(tlist.Index[2],/VALID_ID) eq 1 then WIDGET_CONTROL, tlist.Index[2], /destroy
    tlist.Index[2]=widget_base(title='Residual')
    draw = WIDGET_DRAW(tlist.Index[2], XSIZE=xs, YSIZE=M*ys,/scroll,$
            x_scroll_size=xs, y_scroll_size=ys)
    WIDGET_CONTROL, tlist.Index[2], /REALIZE

    ; ----Calculate and plot Residual----
    ind1=0L
    ind2=0L
    zero=0.00001
    ResVar=fltarr(M)
    for j=0l,M-1 do begin
        ; Interpolate rfit for phiexp datapoints
        x=phifit[ind1:ind1+nphi-1]
        y=rfit[ind1:ind1+nphi-1]
        ; Problems when phifit datapoint are to close (almost zero)
        ; leave them out for interpolation
        ind=where(abs(x-shift(x,1)) gt zero,ct)
        if ct ne 0 then begin
            x=x[ind]
            y=y[ind]
            y2=spl_init(x,y)
            rfit2=SPL_INTERP(x,y, y2, phiexp[ind2:ind2+nj[j]-1])
            resid=(rexp[ind2:ind2+nj[j]-1]-rfit2);/rexp[ind2:ind2+nj[j]-1]
    
            ; SD on residual in d-spacing
            Dtmp1=BraggPtoX(rexp[ind2:ind2+nj[j]-1],tlist.lambda,tlist.dist,tlist.scrx,tlist.scry,tlist.a,tlist.b,tlist.phipixphi0,phi=phiexp[ind2:ind2+nj[j]-1],/onlyd,/nonortp)
            Dtmp2=BraggPtoX(rfit2,tlist.lambda,tlist.dist,tlist.scrx,tlist.scry,tlist.a,tlist.b,tlist.phipixphi0,phi=phiexp[ind2:ind2+nj[j]-1],/onlyd,/nonortp)
            ResVar[j]=stddev(Dtmp1-Dtmp2)
    
            tmp=abs(max(resid,/abs))
            plot,phiexp[ind2:ind2+nj[j]-1],resid,title='Ring '+stringr(j),xtitle='Azimuthal (rad)',$
            ytitle='Difference between obs. and fit (pixels)',ystyle=1,xstyle=1,charsize=M<2,yrange=[-tmp,tmp]*1.1
        endif
        
        ind1+=nphi
        ind2+=nj[j]
    endfor
    !P.MULTI = 0

    ; ----Plot standard deviation on residuals----
    DEVICE, WINDOW_STATE=wstate
    if tlist.Index[3] ne -1 then $
        if wstate[tlist.Index[3]] ne 0 then wdelete,tlist.Index[3]
    tlist.Index[3]=winid()
    window,tlist.Index[3]

    x=TiltRingD[0:M-1]
    indsort=sort(TiltRingD[0:M-1])
    plot,x[indsort],ResVar[indsort],xtitle=ChiXtitle(0),psym=-2,$
    ytitle='Standard deviation on d-spacing residual ('+msymbols('angstroms')+')',title='Standard deviation on d-spacing residual for each ring ('+msymbols('angstroms')+')'

;    p=plot(x[indsort],ResVar[indsort],xtitle=ChiXtitle(0),symbol=2,$
;    ytitle='Standard deviation on d-spacing residual ('+msymbols('angstroms')+')',title='Standard deviation on d-spacing residual for each ring ('+msymbols('angstroms')+')')

endif

; ----Finish----
RefreshDisplay,tlist.top,list
end; pro FitTilt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_tilt, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
WIDGET_CONTROL, list.top, sensitive=1
if list.Index[0] ne -1 then wdelete,list.Index[0]
if WIDGET_INFO(list.Index[1],/VALID_ID) ne 0 then widget_control,list.Index[1],/destroy
if WIDGET_INFO(list.Index[2],/VALID_ID) ne 0 then widget_control,list.Index[2],/destroy
if list.Index[3] ne -1 then wdelete,list.Index[3]
end;pro CleanUp_tilt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro tilt_event, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

evorg=widget_info(ev.id,/type)
if evorg eq 10 then return

widget_control,ev.id,get_value=val
widget_control,ev.top,get_uvalue=tlist

case evorg of
3:    begin;text
    widget_control,ev.id,get_uvalue=uval
    fval=0.
    reads,val,fval
    case uval of
    'a':tlist.a=val/180.*!dpi
    'b':tlist.b=val/180.*!dpi
    'xc':tlist.center[0]=val
    'yc':tlist.center[1]=val
    'ddist':tlist.dist=val*0.01
    else:tlist.k=val
    endcase
    widget_control,ev.top,set_uvalue=tlist
    return
    endcase

1:    begin;button

    case val of
    'Default nr of SD`s':begin
             ID=widget_info(ev.top,FIND_BY_UNAME='k')
             widget_control,ID,set_value='4'
             tlist.k=4
             endcase
    'Refine tilt angle':begin
             tlist.refbool[0]=ev.select
             endcase
    'Refine rotation angle':begin
             tlist.refbool[1]=ev.select
             endcase
    'Use intensities as weights':begin
             tlist.refbool[7]=ev.select
             endcase
    'Refine center':begin
             tlist.refbool[2]=ev.select
             endcase
    'Refine Distance':begin
             tlist.refbool[3]=ev.select
             endcase
    'Refine D-spacings':begin
             tlist.refbool[6]=ev.select
             endcase
    'Refine Wavelength':begin
             tlist.refbool[5]=ev.select
             endcase
    'Refine two-theta':begin
             tlist.refbool[4]=ev.select
             tlist.refbool[[5,6]]=ev.select eq 0
             widget_control,ev.top,set_uvalue=tlist
             ID=widget_info(ev.top,FIND_BY_UNAME='Ref2')
             widget_control,ID,sensitive=ev.select eq 0,$
                 set_button=tlist.refbool[5]
             ID=widget_info(ev.top,FIND_BY_UNAME='Ref3')
             widget_control,ID,sensitive=(ev.select eq 0) and tlist.ref3,$
                 set_button=tlist.refbool[6]
             endcase
    'Reject outliers':begin
             tlist.refbool[8]=ev.select
             endcase
    'Fit with constraints':begin
             tlist.refbool[9]=ev.select
             endcase
    'Visualize fit and outliers':begin
             tlist.refbool[10]=ev.select
             endcase
    'Initialize parameters':begin
             indsel=where(tlist.Nrings eq 1,count)
             if count eq 0 then return
             widget_control,tlist.top,get_uvalue=list

             M=list.nrTiltRing
             nj=(*list.nrTilt)[0:M-1]
             ntotal=total(nj)
             x=transpose((*list.TiltI)[*,0:ntotal-1])
             tt=BraggDtoX((*list.TiltRingD)[0:M-1],tlist.lambda,/onlyt)
             
             ID=widget_info(ev.top,FIND_BY_UNAME='a<0')
             NEGATIVE=widget_info(ID,/button_set)

             list2=RefTilt(M,nj,x,tlist.nRings,tt,tlist.scrx,tlist.scry,NEGATIVE=NEGATIVE)
             tlist.a=list2.a
             tlist.b=list2.b
             tlist.center[0]=list2.xc
             tlist.center[1]=list2.yc
             tlist.dist=list2.dist

             ID=widget_info(ev.top,FIND_BY_UNAME='a')
             widget_control,ID,set_value=stringr(tlist.a*180./!dpi)
             ID=widget_info(ev.top,FIND_BY_UNAME='b')
             widget_control,ID,set_value=stringr(tlist.b*180./!dpi)
             ID=widget_info(ev.top,FIND_BY_UNAME='xc')
             widget_control,ID,set_value=stringr(tlist.center[0])
             ID=widget_info(ev.top,FIND_BY_UNAME='yc')
             widget_control,ID,set_value=stringr(tlist.center[1])
             ID=widget_info(ev.top,FIND_BY_UNAME='ddist')
             widget_control,ID,set_value=stringr(tlist.dist*100)
             endcase
    'Default tilt':begin
             widget_control,tlist.top,get_uvalue=list
             tlist.a=list.a
             ID=widget_info(ev.top,FIND_BY_UNAME='a')
             widget_control,ID,set_value=stringr(tlist.a*180./!dpi)
             endcase
    'Default rotation':begin
             widget_control,tlist.top,get_uvalue=list
             tlist.b=list.b
             ID=widget_info(ev.top,FIND_BY_UNAME='b')
             widget_control,ID,set_value=stringr(tlist.b*180./!dpi)
             endcase
    'Default Xcen':begin
             widget_control,tlist.top,get_uvalue=list
             tlist.center[0]=list.center[0]
             ID=widget_info(ev.top,FIND_BY_UNAME='xc')
             widget_control,ID,set_value=stringr(tlist.center[0])
             endcase
    'Default Ycen':begin
             widget_control,tlist.top,get_uvalue=list
             tlist.center[1]=list.center[1]
             ID=widget_info(ev.top,FIND_BY_UNAME='yc')
             widget_control,ID,set_value=stringr(tlist.center[1])
             endcase
    'Default distance':begin
             widget_control,tlist.top,get_uvalue=list
             tlist.dist=list.dist
             ID=widget_info(ev.top,FIND_BY_UNAME='ddist')
             widget_control,ID,set_value=stringr(tlist.dist*100)
             endcase
      'Calibrate':  begin
             FitTilt,tlist,ev.top
             ID=widget_info(ev.top,FIND_BY_UNAME='a')
             widget_control,ID,set_value=stringr(tlist.a*180./!dpi)
             ID=widget_info(ev.top,FIND_BY_UNAME='b')
             widget_control,ID,set_value=stringr(tlist.b*180./!dpi)
             ID=widget_info(ev.top,FIND_BY_UNAME='xc')
             widget_control,ID,set_value=stringr(tlist.center[0])
             ID=widget_info(ev.top,FIND_BY_UNAME='yc')
             widget_control,ID,set_value=stringr(tlist.center[1])
             endcase
      'OK':  begin
               widget_control,tlist.top,get_uvalue=list
             list.a=tlist.a
             list.b=tlist.b
             list.center=tlist.center
             list.lambda=tlist.lambda
             list.dist=tlist.dist
             ; ----Finish----
             RefreshDisplay,tlist.top,list
             widget_control,ev.top,/destroy
             return
             endcase
    'Cancel':begin
             widget_control,ev.top,/destroy
             return
             endcase
    'Select All':begin
             tlist.Nrings[*]=ev.select
             nbutton=n_elements(tlist.Nrings)
             uname=stringr(indgen(nbutton))
             for i=0l,nbutton-1 do begin
              ID=widget_info(ev.top,FIND_BY_UNAME=uname[i])
              widget_control,ID,set_button=ev.select
             endfor
             indsel=where(tlist.Nrings eq 1,count)
             ID=widget_info(ev.top,FIND_BY_UNAME='Tiltsign')
             widget_control,ID,sensitive=count eq 1
             endcase
       'Tilt > 0': return
       'Tilt < 0': return
       'Distance > 0': return
       'Distance < 0': return
      else:  begin
             i=1L
             reads,val,i
             tlist.Nrings[i]=ev.select
             
             indsel=where(tlist.Nrings eq 1,count)
             ID=widget_info(ev.top,FIND_BY_UNAME='Tiltsign')
             widget_control,ID,sensitive=count eq 1
             endcase
    endcase
    widget_control,ev.top,set_uvalue=tlist
    endcase

endcase

end;pro tilt_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_del_ring, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
widget_control,list.ID1,sensitive=1
end;pro CleanUp_del_ring
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro del_ring_event , ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control, ev.top, get_uvalue=listtemp
list=listtemp.(0)
ID1=listtemp.(1)
select=listtemp.(2)
call=listtemp.(3)
widget_control,ev.id,get_value=val
case val of
'Select All':begin
         listtemp.select[*]=ev.select
         nbutton=n_elements(listtemp.select)
         uname=stringr(indgen(nbutton))
         for i=0l,nbutton-1 do begin
          ID=widget_info(ev.top,FIND_BY_UNAME=uname[i])
          widget_control,ID,set_button=ev.select
         endfor
         widget_control,ev.top,set_uvalue=listtemp
         endcase
    'OK' : begin
         ind=where(select eq 0,count)
         case count of
          0 : begin
              if call eq 'Tilt' then begin
                 list.nrTiltRing=0
                 (*list.nrTilt)=[0]
                 ptr_free,list.TiltI
                 ptr_free,list.TiltRingD
              endif else begin
                 list.nrTieRing=0
                 (*list.nrTie)=[0]
                 ptr_free,list.TieI
                 ptr_free,list.TieO
                 ptr_free,list.TieRingD
              endelse
              endcase
              else:begin
              if Call eq 'Tilt' then begin
                 genind=-1
                 for i=0l,count-1 do begin
                   if ind[i] eq 0 then offset=0 else offset=total((*list.nrTilt)[0:ind[i]-1])
                   if (*list.nrTilt)[ind[i]] ne 0 then genind=[genind,offset+indgen((*list.nrTilt)[ind[i]])]
                 endfor
                 if n_elements(genind) ne 1 then begin
                   genind=genind[1:*]
                   (*list.nrTilt)=(*list.nrTilt)[ind]
                   (*list.TiltRingD)=(*list.TiltRingD)[ind]
                   (*list.TiltI)=(*list.TiltI)[*,genind]
                   list.nrTiltRing=n_elements(*list.nrTilt)-1
                 endif else begin
                     list.nrTiltRing=0
                     (*list.nrTilt)[0]=0
                     ptr_free,list.TiltI
                     ptr_free,list.TiltRingD
                 endelse
              endif else begin
                 genind=-1
                 for i=0l,count-1 do begin
                   if ind[i] eq 0 then offset=0 else offset=total((*list.nrTie)[0:ind[i]-1])
                   if (*list.nrTie)[ind[i]] ne 0 then genind=[genind,offset+indgen((*list.nrTie)[ind[i]])]
                 endfor
                 if n_elements(genind) ne 1 then begin
                   genind=genind[1:*]
                   (*list.nrTie)=(*list.nrTie)[ind]
                   (*list.TieRingD)=(*list.TieRingD)[ind]
                   (*list.TieI)=(*list.TieI)[*,genind]
                   if ptr_valid(list.TieO) then (*list.TieO)=(*list.TieO)[*,genind]
                   list.nrTieRing=n_elements(*list.nrTie)-1
                 endif else begin
                     list.nrTieRing=0
                     (*list.nrTie)=[0]
                     ptr_free,list.TieI
                     ptr_free,list.TieO
                     ptr_free,list.TieRingD
                 endelse
              endelse
              endelse
           endcase
         RefreshDisplay,ID1,list
           WIDGET_CONTROL, ev.top, /DESTROY
           endcase
  else:     begin
         i=1L
         reads,strmid(val,4,strlen(val)-4),i
         listtemp.select[i]=ev.select
         widget_control,ev.top,set_uvalue=listtemp
         endcase
endcase
end;pro del_ring_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DTIE_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN
widget_control,ev.top,get_uvalue=list

if (*list.nrTie)[0] eq 0 then return

; ----Check for zoom----
if ev.id eq list.drawdyn then begin
point=DynTR([ev.x,ev.y],list.sf,[list.hrange[0],list.vrange[0]])
endif else begin
point=StatTR([ev.x,ev.y],[list.sfh,list.sfv])
endelse

;----Find Point selected----
xrange=[point[0]-list.DelWidth,point[0]+list.DelWidth]
yrange=[point[1]-list.DelWidth,point[1]+list.DelWidth]

ind=where( ((*list.TieI)[0,*] ge xrange[0]) and ((*list.TieI)[0,*] le xrange[1]) and $
       ((*list.TieI)[1,*] ge yrange[0]) and ((*list.TieI)[1,*] le yrange[1]) ,count)
if count eq 0 then return
ind=ind[0]

;----Find out in which ring----
s=0l
i=0l
repeat begin
    s+=(*list.nrTie)[i]
    i++
endrep until (ind[0] lt s); Take first ring
Ring=i-1

;----Prompt for d-spacing input----
dstr=PromptNumber(stringr((*list.TieRingD)[ring]),ev.top,'Enter d-spacing ('+msymbols('angstroms')+'):')
d=0.
reads,dstr,d
(*list.TieRingD)[ring]=d
widget_control,ev.top,set_uvalue=list

end;pro DTIE_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro removetiepoints,list,indP
;----Delete Point----
(*list.TieI)=ShrinkArray((*list.TieI),indP,/row)
(*list.TieO)=ShrinkArray((*list.TieO),indP,/row)

;----Find out in which rings----
Ringct=lonarr(list.nrTieRing+1)
ind=0
for i=0l,list.nrTieRing do begin
    indt=where((indP ge ind) and (indP lt ind+(*list.nrTie)[i]),ct)
    Ringct[i]=ct
    ind+=(*list.nrTie)[i]
endfor
(*list.nrTie)-=Ringct

; ----Detect rings with no points----
ind0=where((*list.nrTie) eq 0,ct)
; last ring can have 0 points
if (ct eq 1) and (ind0[0] eq list.nrTieRing) then ct=0
if (ct ne 0) then begin
    if (ind0[ct-1] eq list.nrTieRing) then ind0=ind0[0:ct-2]
    (*list.nrTie)=ShrinkArray((*list.nrTie),ind0)
    list.nrTieRing=(list.nrTieRing-ct)>0
    (*list.TieRingD)=ShrinkArray((*list.TieRingD),ind0)
endif
end;pro removetiepoints
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro removedoubleTieI,ev,list

x=reform((*list.TieI)[0,*])
y=reform((*list.TieI)[1,*])

dd=PromptNumber('5',ev.top,'Delete input tie points within distance (pixels):')
dd=ulong(dd)
if dd lt 0 then return

; ----Find tie points that are "close" together----
n=n_elements(x)

i=rebin(lindgen(n),n,n,/sample)
j=rebin(lindgen(1,n),n,n,/sample)
d=((*list.TieI)[0,i]-(*list.TieI)[0,j])^2.+((*list.TieI)[1,i]-(*list.TieI)[1,j])^2.
d[(n+1)*lindgen(n)]=dd^2+0.1
d lt= dd^2
d=reform(d,n,n,/overwrite)

indP=-1l
for i=0l,n-1 do begin
    ind=where(d[i,*],ct)
    if ct ne 0 then begin
        d[ind,*]=0
        indP=[indP,ind]
    endif
endfor

if n_elements(indP) eq 1 then return
indP=indP[1:*]

; ----Remove tie points----
removetiepoints,list,indP
end;pro removedoubleTieI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DELTIE_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]

case thisEvent of
'UP' :  begin
       widget_control,ev.id,draw_motion=0
       return
       endcase
'DOWN' :begin
       widget_control,ev.id,draw_motion=1
       endcase
'SCROLL' : begin
        XRDUA_event, ev
        return
        endcase
else:   begin
;       possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
;        bPress = possibleButtons[ev.press]
;        if bPress eq 'NONE' then begin
;           widget_control,ev.id,draw_motion=0
;            return
;        endif
       endelse
endcase

widget_control,ev.top,get_uvalue=list
if (*list.nrTie)[0] eq 0 then return

if ev.id eq list.drawdyn then begin
point=DynTR([ev.x,ev.y],list.sf,[list.hrange[0],list.vrange[0]])
endif else begin
point=StatTR([ev.x,ev.y],[list.sfh,list.sfv])
endelse

; ----Find Tie Point selected----
xrange=[point[0]-list.DelWidth,point[0]+list.DelWidth]
yrange=[point[1]-list.DelWidth,point[1]+list.DelWidth]
indP=where( ((*list.TieI)[0,*] ge xrange[0]) and ((*list.TieI)[0,*] le xrange[1]) and $
       ((*list.TieI)[1,*] ge yrange[0]) and ((*list.TieI)[1,*] le yrange[1]) ,count)
if count eq 0 then return
if count ne 1 then begin
    str=stringr((*list.TieI)[*,indP])
    str='['+str[0,*]+';'+str[1,*]+']'
    if (list.DelMass eq 0) then begin
        widget_control,ev.id,draw_motion=0
        sel=ChooseFromList(str,title='Select Tie Points to delete',$
              nsel=nsel,top=ev.top)
        if nsel eq 0 then return
        indP=indP[sel]
    endif
endif

; ----Remove tie points----
removetiepoints,list,indP

RefreshDisplay,ev.top,list
end;pro DELTIE_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TIE_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN
widget_control,ev.top,get_uvalue=list

if ev.id eq list.drawdyn then begin
point=DynTR([ev.x,ev.y],list.sf,[list.hrange[0],list.vrange[0]])
endif else begin
point=StatTR([ev.x,ev.y],[list.sfh,list.sfv])
endelse

if (*list.nrTie)[0] eq 0 then begin
    list.TieI=ptr_new(point)
    if list.TieGrid then point=GridCoord(point,list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix)
    list.TieO=ptr_new(point)
    list.TieRingD=ptr_new([0.])
endif else begin
    (*list.TieI)=[[(*list.TieI)],[point]]
    if list.TieGrid then point=GridCoord(point,list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix)
    (*list.TieO)=[[(*list.TieO)],[point]]
    if (*list.nrTie)[list.nrTieRing] eq 0 then (*list.TieRingD)=[(*list.TieRingD),0.]
endelse
(*list.nrTie)[list.nrTieRing]++

RefreshDisplay,ev.top,list
end;pro TIE_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DTILT_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN
widget_control,ev.top,get_uvalue=list

if (*list.nrTilt)[0] eq 0 then return

if ev.id eq list.drawdyn then begin
point=DynTR([ev.x,ev.y],list.sf,[list.hrange[0],list.vrange[0]])
endif else begin
point=StatTR([ev.x,ev.y],[list.sfh,list.sfv])
endelse

;----Find Point selected----
xrange=[point[0]-list.DelWidth,point[0]+list.DelWidth]
yrange=[point[1]-list.DelWidth,point[1]+list.DelWidth]
ind=where( ((*list.TiltI)[0,*] ge xrange[0]) and ((*list.TiltI)[0,*] le xrange[1]) and $
       ((*list.TiltI)[1,*] ge yrange[0]) and ((*list.TiltI)[1,*] le yrange[1]) ,count)
if count eq 0 then return
ind=ind[0]

;----Find out in which ring----
s=0
i=0
repeat begin
s=s+(*list.nrTilt)[i]
i=i+1
endrep until (ind[0] lt s); take first ring
Ring=i-1

;----Prompt for d-spacing input----
dstr=PromptNumber(stringr((*list.TiltRingD)[ring]),ev.top,'Enter d-spacing ('+msymbols('angstroms')+'):')
d=0.
reads,dstr,d
(*list.TiltRingD)[ring]=d
widget_control,ev.top,set_uvalue=list

end;pro DTILT_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DELTILT_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]

case thisEvent of
'UP' :  begin
           widget_control,ev.id,draw_motion=0
           return
           endcase
'DOWN' :begin
        widget_control,ev.id,draw_motion=1
        endcase
'SCROLL' : begin
        XRDUA_event, ev
        return
        endcase
else:   begin
;       possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
;        bPress = possibleButtons[ev.press]
;        if bPress eq 'NONE' then begin
;           widget_control,ev.id,draw_motion=0
;            return
;        endif
       endelse
endcase

widget_control,ev.top,get_uvalue=list
if (*list.nrTilt)[0] eq 0 then return

if ev.id eq list.drawdyn then begin
point=DynTR([ev.x,ev.y],list.sf,[list.hrange[0],list.vrange[0]])
endif else begin
point=StatTR([ev.x,ev.y],[list.sfh,list.sfv])
endelse

;----Find Point selected----
xrange=[point[0]-list.DelWidth,point[0]+list.DelWidth]
yrange=[point[1]-list.DelWidth,point[1]+list.DelWidth]
indP=where( ((*list.TiltI)[0,*] ge xrange[0]) and ((*list.TiltI)[0,*] le xrange[1]) and $
       ((*list.TiltI)[1,*] ge yrange[0]) and ((*list.TiltI)[1,*] le yrange[1]) ,count)
if count eq 0 then return
if count ne 1 then begin
    str=stringr(Point[*,indP])
    str='['+str[0,*]+';'+str[1,*]+']'
    if (list.DelMass eq 0) then begin
        widget_control,ev.id,draw_motion=0
        sel=ChooseFromList(str,title='Select Tilt Points to delete',$
           nsel=nsel,top=ev.top)
        if nsel eq 0 then return
        indP=indP[sel]
    endif
endif

;----Delete Point----
(*list.TiltI)=ShrinkArray((*list.TiltI),indP,/row)

;----Find out in which rings----
Ringct=lonarr(list.nrTiltRing+1)
ind=0
for i=0l,list.nrTiltRing do begin
    indt=where((indP ge ind) and (indP lt ind+(*list.nrTilt)[i]),ct)
    Ringct[i]=ct
    ind=ind+((*list.nrTilt))[i]
endfor
(*list.nrTilt)-=Ringct

; ----Detect rings with no points----
ind0=where((*list.nrTilt) eq 0,ct)
; last ring can have 0 points
if (ct eq 1) and (ind0[0] eq list.nrTiltRing) then ct=0
if (ct ne 0) then begin
    if (ind0[ct-1] eq list.nrTieRing) then ind0=ind0[0:ct-2]
    (*list.nrTilt)=ShrinkArray((*list.nrTilt),ind0)
    list.nrTiltRing=(list.nrTiltRing-ct)>0
    (*list.TiltRingD)=ShrinkArray((*list.TiltRingD),ind0)
endif

RefreshDisplay, ev.top,list
end;pro DELTILT_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TILT_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN
widget_control,ev.top,get_uvalue=list

if ev.id eq list.drawdyn then begin
point=DynTR([ev.x,ev.y],list.sf,[list.hrange[0],list.vrange[0]])
endif else begin
point=StatTR([ev.x,ev.y],[list.sfh,list.sfv])
endelse

if (*list.nrTilt)[0] eq 0 then begin
    list.TiltI=ptr_new(point)
    list.TiltRingD=ptr_new([0.])
endif else begin
    (*list.TiltI)=[[(*list.TiltI)],[point]]
    if (*list.nrTilt)[list.nrTiltRing] eq 0 then (*list.TiltRingD)=[(*list.TiltRingD),0.]
endelse
(*list.nrTilt)[list.nrTiltRing]++

RefreshDisplay, ev.top,list
end;pro TILT_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fittdGrid,img,v,w,r,sqr,imgvalid,reduce,OutID,drawi=drawi,mag=mag,$
    sclmax=sclmax,debug=debug,search=search

; ----Estimate peak parameters----
; A,mx,my,sx,sy,sxy=0,n
peaks=imgpeaks(img,v,w,drawi,mag,r,imgvalid,sclmax=sclmax)
s=size(peaks)
if s[1] ne 7 then begin
    printw,OutID,'No peaks found'
    flplotid=0L
    return,0
endif
if s[0] eq 1 then apeaks=1 else apeaks=s[2]


if keyword_set(search) then begin
    peaks[1,*]=peaks[1,*]+sqr[3]
    peaks[2,*]=peaks[2,*]+sqr[5]
    return,peaks
endif

; ----Fit img by fitting each peak seperatly----
; No saturation removal implemented for this
c=5 ; Define number of std to make sub-image
xl=sqr[4]-sqr[3]
yl=sqr[6]-sqr[5]
simg=size(img)
for i=0l,apeaks-1 do begin
    AA=peaks[0:5,i]
    ;xb=round(AA[1]-c*AA[3])>0
    ;yb=round(AA[2]-c*AA[4])>0
    ;xe=round(AA[1]+c*AA[3])<xl
    ;ye=round(AA[2]+c*AA[4])<yl
    xb=round(AA[1]-reduce/2.)>0
    yb=round(AA[2]-reduce/2.)>0
    xe=round(AA[1]+reduce/2.)<xl
    ye=round(AA[2]+reduce/2.)<yl
    x=indgen(xe-xb+1)+xb
    y=indgen(ye-yb+1)+yb
    axbin=n_elements(x)
    aybin=n_elements(y)
    n=axbin*aybin
    surf=reform(img[xb:xe,yb:ye],n)
    weight=1./(surf>1)
    xx={x:x,y:y}
    imgfit = Marquardt(xx,surf, weight, [AA,0], SIGMA,FUNCTION_NAME='gfunc2dg',$
        CHISQ=CHISQ,ITER=ITER,itmax=100,nfree=nfree)
    peaks[0:5,i]=AA[0:5]
    ; ----compare fit with data----
    if keyword_set(debug) then begin
        struct={surf:reform(surf,axbin,aybin),imgfit:reform(imgfit,axbin,aybin),$
           x:x,y:y,xr:[min(x),max(x)],yr:[min(y),max(y)],$
           zr:[min(img),max(img)],proplot:'plotcomp'}
        flplot_obj,struct,topid=flplotid,/xtrue,/ytrue
        widget_control,flplotid,/destroy
    endif
endfor

peaks[1,*]=peaks[1,*]+sqr[3]
peaks[2,*]=peaks[2,*]+sqr[5]
return,peaks

end;function fittdGrid
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_GridSearch,ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
ptr_free,list.GP
end;pro CleanUp_GridSearch
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GridSearch_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list2
if widget_info(ev.id,/type) EQ 3 THEN BEGIN ; event van text_widget
    widget_control,ev.id,get_uvalue=uvalue
    widget_control,ev.id,get_value=b
    ID=list2.ID1
    list=list2.list
    case uvalue of
    'reduce':begin
           reduce=0.
           reads,b,reduce
           list2.reduce=reduce
           widget_control,ev.top,set_uvalue=list2
           endcase
    'v':     begin
          vin=0
             reads,b,vin
             list.v=vin
             list2.list=list
            widget_control,ID,set_uvalue=list
            widget_control,ev.top,set_uvalue=list2
            endcase
    'w':        begin
          win=0
             reads,b,win
             list.w=win
            list2.list=list
            widget_control,ID,set_uvalue=list
            widget_control,ev.top,set_uvalue=list2
            endcase
    'r':        begin
          rin=0
             reads,b,rin
             list.r=rin
            list2.list=list
            widget_control,ID,set_uvalue=list
            widget_control,ev.top,set_uvalue=list2
            endcase
    endcase
endif ELSE  begin
    widget_control,ev.id,get_value=val
    ID=list2.ID1
    list=list2.list
    case val of
    'Exit' :   begin
              widget_control,ev.top,/destroy
              return
              endcase
    'Add':  begin
           if list2.nGP ne 0 then begin
               n=list2.nGP
               offset=(*list.nrTie)[0]
               peakso=GridCoord(*list2.GP,list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix)
               if (*list.nrTie)[0] eq 0 then begin
                       list.TieI=ptr_new(*list2.GP)
                       list.TieO=ptr_new(peakso)
                       (*list.nrTie)[0]=n
               endif else begin
                   (*list.TieI)=[[(*list.TieI)],[*list2.GP]]
                   (*list.TieO)=[[(*list.TieO)],[peakso]]
                      (*list.nrTie)[0]+=n
               endelse

               ; Filter out close pixels
               halfshift=0.5
               TieImin=(*list.TieI)-halfshift
               TieIMax=(*list.TieI)+halfshift
               nrTie=(*list.nrTie)[0]
               i=0
               while i lt nrTie-1 do begin
                    ind=where( ((*list.TieI)[0,*] ge TieImin[0,i]) and ((*list.TieI)[0,*] le TieIMax[0,i]) and $
                     ((*list.TieI)[1,*] ge TieImin[1,i]) and ((*list.TieI)[1,*] le TieIMax[1,i]),ct )
                    if ct gt 1 then begin
                        this=where(ind ne i,ct)
                        (*list.TieI)=ShrinkArray((*list.TieI),ind[this],/row,count=nrTie)
                        (*list.TieO)=ShrinkArray((*list.TieO),ind[this],/row)
                        TieImin=ShrinkArray(TieImin,ind[this],/row)
                        TieImax=ShrinkArray(TieImax,ind[this],/row)
                        ;ind2=[ind2,ind[this]]
                    endif
                    i=i+1
               endwhile

               (*list.nrTie)[0]=nrTie

               RefreshDisplay,ID,list
            endif
           endcase
    'Search':  begin
          tvlct,RoldCT,GoldCT,BoldCT,/get
          LoadctRev,3,/silent
          v=list.v
          w=list.w
          sqr=list2.sqr
          drawi=list2.drawi
          mag=list2.mag
          r=list.r

          img=(*list.tiff)[sqr[3]:sqr[4],sqr[5]:sqr[6]]
          if ptr_valid(list.tiffvalid) then imgvalid=(*list.tiffvalid)[sqr[3]:sqr[4],sqr[5]:sqr[6]]

          printw,list.OutID,'Finding spots...'
          widget_control, /hourglass
          T1=systime(1)
          peaks=fittdGrid(img,v,w,r,sqr,imgvalid,list2.reduce,list.OutID,$
            drawi=drawi,mag=mag,sclmax=list.sclmax)
          T2=systime(1)

          s=size(peaks)
          if s[1] ne 7 then apeaks=0 else $
          if s[0] eq 1 then apeaks=1 else apeaks=s[2]

          if apeaks ne 0 then begin
              ; ----Reduce peaks by user supplied criterium----
              ; Reduce on basis of width
              if list2.reduce ne 0 then begin
                 c=2*sqrt(2*alog(2))
                 ind=where(c*peaks[3,*] lt list2.reduce and c*peaks[4,*] lt list2.reduce,count)
                 if count eq 0 then begin
                     apeaks=0
                 endif else begin
                     peaks=peaks[*,ind]
                     apeaks=count
                 endelse
              endif
              if apeaks ne 0 then begin
                  list2.nGP=apeaks
                  ptr_free,list2.GP
                  list2.GP=ptr_new(peaks[1:2,*])
                  list2.list=list
                  WIDGET_CONTROL, ev.top,set_uvalue=list2
              endif
          endif

         ; ----Plot results----
         if apeaks ne 0 then begin
            s=size(img)
             xx=reform(peaks[1,*])-sqr[3]
             yy=reform(peaks[2,*])-sqr[5]
             outimg=img*0.
             outimg[xx,yy]=255.
             LoadctRev,-3,/silent
             wset,drawi
             if nsatur ne 0 then img[saturarray]=0
             tv,rebin(outimg,s[1]*mag,s[2]*mag)+ $
               rebin(WeakScl(img,list.sclmax,list.displaypower),s[1]*mag,s[2]*mag,/sample)
         endif else begin
            s=size(img)
             LoadctRev,-3,/silent
             wset,drawi
             if nsatur ne 0 then img[saturarray]=0
             tv,rebin(WeakScl(img,list.sclmax,list.displaypower),s[1]*mag,s[2]*mag,/sample)
         endelse

          printw,list.OutID,'Number of peaks detected: '+string(apeaks)
          printw,list.OutID,'Peaksearch time in sec: '+string(T2-T1)
          tvlct,RoldCT,GoldCT,BoldCT
          endcase
    endcase
endelse
end;pro GridSearch_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro surffit_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list2
if widget_info(ev.id,/type) EQ 3 THEN BEGIN ; event van text_widget
    widget_control,ev.id,get_uvalue=uvalue
    widget_control,ev.id,get_value=b
    ID=list2.ID1
    list=list2.list
    case uvalue of
    'v':     begin
          vin=0
             reads,b,vin
             list.v=vin
             list2.list=list
            widget_control,ID,set_uvalue=list
            widget_control,ev.top,set_uvalue=list2
            endcase
    'w':        begin
          win=0
             reads,b,win
             list.w=win
            list2.list=list
            widget_control,ID,set_uvalue=list
            widget_control,ev.top,set_uvalue=list2
            endcase
    'r':        begin
          rin=0
             reads,b,rin
             list.r=rin
            list2.list=list
            widget_control,ID,set_uvalue=list
            widget_control,ev.top,set_uvalue=list2
            endcase
    endcase
endif ELSE  begin
    widget_control,ev.id,get_value=val
    ID1=list2.ID1
    list=list2.list
    case val of
    'Constraints':begin
             list.fit=0
            list2.list=list
            widget_control,ID1,set_uvalue=list
            widget_control,ev.top,set_uvalue=list2
            endcase
    'Free':   begin
             list.fit=2
            list2.list=list
            widget_control,ID1,set_uvalue=list
            widget_control,ev.top,set_uvalue=list2
            endcase
    'Exit' :   begin
          widget_control,ev.top,/destroy
          return
          endcase
    'Search':  begin
          widget_control, /hourglass
          tvlct,RoldCT,GoldCT,BoldCT,/get
          LoadctRev,3,/silent
          v=list.v
          w=list.w
          sqr=list2.sqr
          drawi=list2.drawi
          mag=list2.mag
          r=list.r

          img=(*list.tiff)[sqr[3]:sqr[4],sqr[5]:sqr[6]]
          if ptr_valid(list.tiffvalid) then imgvalid=(*list.tiffvalid)[sqr[3]:sqr[4],sqr[5]:sqr[6]]

          T1=systime(1)
          peaks=imgpeaks(img,v,w,drawi,mag,r,imgvalid,sclmax=list.sclmax)
          T2=systime(1)
          s=size(peaks)
          if s[1] ne 7 then apeaks=0 else $
          if s[0] eq 1 then apeaks=1 else apeaks=s[2]
          printw,list.OutID,'Number of peaks detected: '+string(apeaks)
          printw,list.OutID,'Peaksearch time in sec: '+string(T2-T1)
          tvlct,RoldCT,GoldCT,BoldCT
          endcase
    'Calc':       begin
          widget_control, /hourglass
          tvlct,RoldCT,GoldCT,BoldCT,/get
          LoadctRev,3,/silent

          sqr=list2.sqr
          drawi=list2.drawi
          mag=list2.mag

          img=float((*list.tiff)[sqr[3]:sqr[4],sqr[5]:sqr[6]])
          if ptr_valid(list.tiffvalid) then imgvalid=(*list.tiffvalid)[sqr[3]:sqr[4],sqr[5]:sqr[6]]

          peaks=fittd(img,list.v,list.w,list.r,0>fix(list.fit)<2,sqr,$
              imgvalid,sclmax=list.sclmax,drawi=drawi,mag=mag,$
              file=list.file,path=list.path,printres=list.OutID,/debug)
          tvlct,RoldCT,GoldCT,BoldCT
          endcase
    endcase
endelse
end;pro surffit_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_surf_int, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
WIDGET_CONTROL, list.(1), sensitive=1
end;pro CleanUp_surf_int
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro surf_int_event, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_value=val
widget_control,ev.top,get_uvalue=tlist
widget_control,ev.top,/destroy
list=tlist.(0)
ID1=tlist.(1)
process=tlist.(2)

if val eq 'Exit' then begin
return
endif

v=list.v
w=list.w
r=list.r
fit=list.fit
sqr=(*list.mask)[*,val]
mag=6
xs=sqr[4]-sqr[3]
ys=sqr[6]-sqr[5]
s=float(xs>ys)
mag=mag-s/200
mag=fix(mag)>1
xs=xs+1
ys=ys+1

case process of
0:  begin
    list2={list:list,ID1:ID1,sqr:sqr,drawi:0L,mag:mag}
    base=widget_base(/column,title='Surface Fitting',xoffset=200,yoffset=200,uvalue=list2)
    base1=widget_base(base,/row)
    basep=widget_base(base1,/column)
       label=widget_label(basep,value='PeakSearch Parameters:')
       basepp=widget_base(basep,/row)
       base2=widget_base(basepp,/column)
         wlabel=widget_label(base2,value='w = FWHM peak in pixels:',ysize=25)
         vlabel=widget_label(base2,value='v in [FWHM/3,FWHM/2]:',ysize=25)
         rlabel=widget_label(base2,value='r = threshold r.stdev:',ysize=25)
       base3=widget_base(basepp,/column)
         wtext=widget_text(base3,value=stringr(w),/editable,uvalue='w')
         vtext=widget_text(base3,value=stringr(v),/editable,uvalue='v')
         rtext=widget_text(base3,value=stringr(r),/editable,uvalue='r')
    basep=widget_base(base1,/column)
       label=widget_label(basep,value='Parameter Boundaries:')
       base4=widget_base(basep,/column,/exclusive)
         button=lonarr(3)
         button[0]=widget_button(base4,value='Constraints')
         button[2]=widget_button(base4,value='Free')
    Sbutton=widget_button(base,value='Search')
    calcbutton=widget_button(base,value='Calc')
    releasebutton=widget_button(base,value='Exit')

    wset, list.drawindex
    draw = WIDGET_DRAW(base,xsize=mag*xs, ysize=mag*ys, /scroll,$
         x_scroll_size=list.stath, y_scroll_size=list.statv,uname='draw')

    WIDGET_CONTROL, base, /REALIZE
    widget_control,button[list.fit],set_button=1
    Xmanager,'surf_int_event',base, event_handler='surffit_event',GROUP_LEADER=ID1
    widget_control,draw,get_value=drawi
    list2.drawi=drawi
    widget_control, base, set_uvalue=list2
    endcase
1:  begin
    widget_control, /hourglass
    tvlct,RoldCT,GoldCT,BoldCT,/get
    LoadctRev,3,/silent

    wset, list.drawindex
    base=widget_base(title='Center Fit')
    draw = WIDGET_DRAW(base, XSIZE=xs*mag, YSIZE=ys*mag,/scroll,$
        x_scroll_size=list.stath, y_scroll_size=list.statv)
    WIDGET_CONTROL, base, /REALIZE
    widget_control,draw,get_value=drawi

    ; ----Select ROI and background substr----
    img=float((*list.tiff)[sqr[3]:sqr[4],sqr[5]:sqr[6]])
    if ptr_valid(list.tiffvalid) then tiffvalid=(*list.tiffvalid)[sqr[3]:sqr[4],sqr[5]:sqr[6]]

    peaks=fittd(img,list.v,list.w,list.r,0>fix(list.fit)<2,sqr,$
       tiffvalid,sclmax=list.sclmax,drawi=drawi,mag=mag,$
       printres=list.OutID,/debug)
    tvlct,RoldCT,GoldCT,BoldCT,/get
    s=size(peaks)
    case s[0] of
    0:    printw,list.OutID,'No peaks found'
    1:    begin
        list.center=peaks[1:2]
        str=stringr(list.center)
        ID=widget_info(ID1,FIND_BY_UNAME='center')
        widget_control,ID,set_value=str[0]+','+str[1]
        RefreshDisplay,ID1,list
        endcase
    else: begin
        str=stringr(peaks[1:2,*])
        str='['+str[0,*]+';'+str[1,*]+']'
        sel=ChooseFromList(str,title='Select center',$
                       nsel=nsel,top=ev.top)
        if nsel eq 1 then begin
            list.center=peaks[1:2]
            str=stringr(list.center)
            ID=widget_info(ID1,FIND_BY_UNAME='center')
            widget_control,ID,set_value=str[0]+','+str[1]
            RefreshDisplay,ID1,list
        endif else begin
            printw,list.OutID,'Select 1 center'
        endelse
        endcase
    endcase
    endcase
2:  begin
    widget_control,ID1,sensitive=ID1
    list2={list:list,ID1:ID1,sqr:sqr,drawi:0L,mag:mag,GP:PTR_NEW(),nGP:0L,reduce:10.}
    base=widget_base(/column,title='Find Grid Points',xoffset=200,yoffset=200,uvalue=list2)
    base1=widget_base(base,/row)
    basep=widget_base(base1,/column)
       label=widget_label(basep,value='GridSearch Parameters:')
       basepp=widget_base(basep,/row)
       base2=widget_base(basepp,/column)
         vlabel=widget_label(base2,value='v:',ysize=25)
         wlabel=widget_label(base2,value='w:',ysize=25)
         rlabel=widget_label(base2,value='r:',ysize=25)
         redlabel=widget_label(base2,value='maximum peak width (pixels):',ysize=25)
       base3=widget_base(basepp,/column)
         vtext=widget_text(base3,value=stringr(v),/editable,uvalue='v')
         wtext=widget_text(base3,value=stringr(w),/editable,uvalue='w')
         rtext=widget_text(base3,value=stringr(r),/editable,uvalue='r')
         reducetext=widget_text(base3,value=stringr(list2.reduce),/editable,uvalue='reduce')

    Sbutton=widget_button(base,value='Search')
    Add=widget_button(base,value='Add')
    releasebutton=widget_button(base,value='Exit')

    wset, list.drawindex
    draw = WIDGET_DRAW(base,xsize=mag*xs, ysize=mag*ys, /scroll,$
         x_scroll_size=list.stath, y_scroll_size=list.statv,uname='draw')

    WIDGET_CONTROL, base, /REALIZE
    Xmanager,'surf_int_event',base, event_handler='GridSearch_event',GROUP_LEADER=ID1,cleanup='CleanUp_GridSearch'
    widget_control,draw,get_value=drawi
    list2.drawi=drawi
    widget_control, base, set_uvalue=list2
    endcase
endcase

end;pro surf_int_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro surf_int,list,process,ID1,title

if list.masknr eq 0 then return
WIDGET_CONTROL, ID1, sensitive=0
base=widget_base(/column,title=title,xoffset=200,yoffset=200,uvalue={list:list,ID1:ID1,process:process})

; ----Mask buttons----
nr=10
nbase=ceil(float(list.masknr)/nr)
base1=lonarr(nbase)
for i=0l,nbase-1 do base1[i]=widget_base(base,/row)
button = lonarr(list.masknr)
j=0
for i=0l,list.masknr-1 do begin
    button[i]=widget_button(base1[j],value=stringr(i))
    if (*list.mask)[0,i] ne 1 then widget_control,button[i],sensitive=0
    if ((i+1) mod nr) eq 0 then j=j+1
endfor

releasebutton=widget_button(base,value='Exit',uname='Exit',xsize=250)
WIDGET_CONTROL, base, /REALIZE
Xmanager,'surf_int',base, event_handler='surf_int_event',cleanup='CleanUp_surf_int',GROUP_LEADER=ID1
end;pro surf_int
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_rectangle_sum, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=top
WIDGET_CONTROL, top, sensitive=1
end;pro CleanUp_rectangle_sum
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro rectangle_sum_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_value=val
case val of
'Horizontal summation': return
'Exit':    begin
        widget_control,ev.top,/destroy
        return
        endcase
else:
endcase

ID=widget_info(ev.top,find_by_uname='Horizontal summation')
bhor=widget_info(ID,/button_set)
widget_control,ev.top,get_uvalue=top
widget_control,top,get_uvalue=list

sqr=(*list.mask)[*,val]
img=float((*list.tiff)[sqr[3]:sqr[4],sqr[5]:sqr[6]])
s = *list.tiffs
if ptr_valid(list.tiffvalid) then img*=(*list.tiffvalid)[sqr[3]:sqr[4],sqr[5]:sqr[6]]

if bhor then begin
    img=reform(total(img,1))
    x=lindgen(n_elements(img))*list.scry*1000
endif else begin
    img=total(img,2)
    x=lindgen(n_elements(img))*list.scrx*1000
endelse

window
plot,x,img,xtitle='Distance (mm)',ytitle='Summed intensity',/xs

path=selfile(list.chipath,list.file+'.txt','*.txt','Save profile...')
if path eq '' then return
tmp=CutPath(path,path=chipath,file=chifile,ext=ext)
list.chipath=chipath
list.chifile=chifile+ext
widget_control,top,set_uvalue=list

if openw_safe(lun,path) then return
printf,lun,transpose([[x],[img]])
free_lun,lun

printw,list.outid,'File saved: '+path

end;pro rectangle_sum_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro rectangle_sum,list,top,title

if list.masknr eq 0 then return
WIDGET_CONTROL, top, sensitive=0
base=widget_base(/column,title=title,xoffset=200,yoffset=200,uvalue=top)

; ----Mask buttons----
nr=10
nbase=ceil(float(list.masknr)/nr)
base1=lonarr(nbase)
for i=0l,nbase-1 do base1[i]=widget_base(base,/row)
button = lonarr(list.masknr)
j=0
for i=0l,list.masknr-1 do begin
    button[i]=widget_button(base1[j],value=stringr(i))
    if (*list.mask)[0,i] ne 1 then widget_control,button[i],sensitive=0
    if ((i+1) mod nr) eq 0 then j=j+1
endfor

base1=widget_base(base,/nonexclusive)
direction=widget_button(base1,value='Horizontal summation',uname='Horizontal summation',xsize=250)
releasebutton=widget_button(base,value='Exit',uname='Exit',xsize=250)
WIDGET_CONTROL, base, /REALIZE
Xmanager,'rectangle_sum',base, event_handler='rectangle_sum_event',cleanup='CleanUp_rectangle_sum',GROUP_LEADER=top
end;pro rectangle_sum
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bin_findpeaks,image,rimage,imgvalid,xoff,yoff,count=ct

debug=0b

s=size(image,/dimensions)
sr=size(rimage,/dimensions)
bcheckvalid=n_elements(imgvalid) ne 0 

; Get blob indices:
b = LABEL_REGION(image)

; Get population and members of each blob:
h = HISTOGRAM(b, REVERSE_INDICES=r)
hval=image[r[r[indgen(n_elements(h))]]]

; Each region
ind=where(h gt 1 and hval eq 1,ct)

if ct ne 0 then begin
    m=fltarr(2,ct)
    mkeep=make_array(ct,value=1b)
    if debug then all=fltarr(9,ct)
    for j=0l, ct-1 DO BEGIN
        i=ind[j]
       ;Find subscripts of members of region i.
       p = r[r[i]:r[i+1]-1]
    
       ; Pixels of region i
       x=(p mod s[0])-1
       y=(p/s[0])-1
       sx=stdev(x,mx)>1
       sy=stdev(y,my)>1
       area=rimage[mx,my]*2.*!dpi*sx*sy
       AA=[area,mx,my,sx,sy,0.,0]
    
       ; Fit part of real image
       mix=min(x,max=max)
       miy=min(y,max=may)
       dx=(max-mix)>10
       dy=(may-miy)>10
       
       xb=(mix-dx)>0
       xe=(max+dx)<(sr[0]-1)
       yb=(miy-dy)>0
       ye=(may+dx)<(sr[1]-1)
    
       x=indgen(xe-xb+1)+xb
       y=indgen(ye-yb+1)+yb
       axbin=n_elements(x)
       aybin=n_elements(y)
       n=axbin*aybin
       surf=reform(rimage[x[0]:x[axbin-1],y[0]:y[aybin-1]],n)
       AA[6]=min(surf)
       
       if debug then begin
           xmap=x#replicate(1,aybin)
           ymap=replicate(1,axbin)#y
           xmap=reform(xmap,axbin*aybin,/overwrite)
           ymap=reform(ymap,axbin*aybin,/overwrite)
           xx=[[temporary(xmap)],[temporary(ymap)]]
           tmp1=where(surf gt (max(surf)-40),tmp2)
           surf2=float(surf)
           if tmp2 ne 0 then surf2[tmp1]=!values.f_nan
           imgfit = NLLS(xx,surf2,AA,'gfunc2dgbis',weight=1,sigma=sigma,iter=iter,chisq=chisq,dof=dof)
           gfunc2dgbis,xx,AA,imgfit
           AA[1:2]+=[xoff,yoff]
           
           ; Show fit:
           all[*,j]=[AA,iter,chisq/dof]
           if j eq 0 then begin
               x+=xoff
               y+=yoff
               xr=[min(x),max(x)]
            yr=[min(y),max(y)]
            list={surf:reform(surf,axbin,aybin),imgfit:reform(imgfit,axbin,aybin),x:x,y:y,xr:xr,yr:yr,zr:[min(surf),max(surf)],proplot:'backplot'}
           
               LoadctRev,1,/silent
            flplot_obj,list,/xtrue,/ytrue,handle=handle
    ;        wait,0.5
    ;           WIDGET_CONTROL, handle,/destroy
           endif
       endif else begin

           xx={x:x,y:y}
           weight=1./(surf>1)
           if bcheckvalid then weight*=reform(imgvalid[x[0]:x[axbin-1],y[0]:y[aybin-1]],n)
           
           imgfit = marquardt(xx,surf, weight, AA, SIGMA,FUNCTION_NAME='gfunc2dg',$
                CHISQ=CHISQ,ITER=ITER,itmax=100)

            if bcheckvalid then mkeep[j]=imgvalid[round(AA[1]),round(AA[2])]

;            ;plot
;           if iter lt 0 then begin
;                erase
;               loadctrev,0,rev=-1
;               tmp=reform(surf,axbin,aybin);*imgvalid[x[0]:x[axbin-1],y[0]:y[aybin-1]]
;               mag=5
;               tvscl,rebin(tmp,mag*axbin,mag*aybin,/sample)
;               
;               loadctrev,3
;               xy=RTDyn(AA[1:2],1./mag,[x[0],y[0]])
;               plots,xy[[0,0]],[0,mag*aybin],/device,color=150
;               plots,[0,mag*axbin],xy[[1,1]],/device,color=150
;
;               x+=xoff
;               y+=yoff
;               xr=[min(x),max(x)]
;               yr=[min(y),max(y)]
;               list={surf:reform(surf,axbin,aybin),imgfit:reform(imgfit,axbin,aybin),x:x,y:y,xr:xr,yr:yr,zr:[min(surf),max(surf)],proplot:'backplot'}
;               
;               LoadctRev,1,/silent
;               flplot_obj,list,/xtrue,/ytrue,handle=handle
;               WIDGET_CONTROL, handle,/destroy
;                stop
;            endif
            
            AA[1:2]+=[xoff,yoff]
    
       endelse
       
       m[*,j]=AA[1:2]
    ENDFOR
    if debug then begin
        file=selfile('C:\tmp\evert\','pts.csv','*.csv','')
        mwrite_csv,file,all,colnames=['Volume','mx','my','sx','sy','sxy','const','Niter','Chisq/DOF'],sep=';',dec=','
    endif
    
endif else m=[-1,-1]

ind=where(mkeep,ct)
if ct eq 0 then return,[-1,-1] else return,m[*,ind]
end;function bin_findpeaks
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro easy_grid,ev

widget_control,ev.top,get_uvalue=list
ptr_free,list.SpaDisStruc

; Select ROI
if list.masknr ne 0 then begin
    ind=where((*list.mask)[0,*] eq 1,ct)
    if ct eq 0 then sqr=[0,(*list.tiffs)[1]-1,0,(*list.tiffs)[2]-1] $
    else sqr=(*list.mask)[3:6,ind[0]]
endif else sqr=[0,(*list.tiffs)[1]-1,0,(*list.tiffs)[2]-1]
xoff=sqr[0]
xend=sqr[1]
yoff=sqr[2]
yend=sqr[3]
image=(*list.tiff)[xoff:xend,yoff:yend]
if ptr_valid(list.tiffvalid) then imgvalid=(*list.tiffvalid)[xoff:xend,yoff:yend]
pixmax=[xend-xoff,yend-yoff]
tiffs=pixmax+1

; Number of subimages in x and y
nmax=(*list.tiffs)[1:2]/float(list.GridDim)*5
if list.GridDim[0] eq 0 then nmax[0]=10
if list.GridDim[1] eq 0 then nmax[1]=10
n=floor(tiffs/nmax)>1

; Subimage size
npixblock=ceil(tiffs/float(n))

; Loop over subimages
bimage=bytarr(tiffs[0],tiffs[1])
m=n[0]*n[1]
for i=0l,m-1 do begin
    ; Subimage boundries
      x=i mod n[0]
      x=npixblock[0]*[x,x+1]
      x[1]--
      x<=pixmax[0]

      y=i/n[0]
      y=npixblock[1]*[y,y+1]
      y[1]--
      y<=pixmax[1]

    ; Binarize subimage
    if x[1]-x[0] lt 5 or y[1]-y[0] lt 5 then img=image[x[0]:x[1],y[0]:y[1]] $
    else img=smooth(image[x[0]:x[1],y[0]:y[1]],5,/EDGE_TRUNCATE,missing=0)

    bimage[x[0]:x[1],y[0]:y[1]]=img gt (max(img)+median(img))*0.5
endfor

; find peaks
peaks=bin_findpeaks(bimage,image,imgvalid,xoff,yoff,count=n)
if n eq 0 then return

; Add peaks
ind=where(peaks[0,*] ge 0 and peaks[0,*] lt (*list.tiffs)[1] and peaks[1,*] ge 0 and peaks[1,*] lt (*list.tiffs)[2],n)
if n eq 0 then return
peaks=peaks[*,ind]

peakso=GridCoord(peaks,list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix)
if (*list.nrTie)[0] eq 0 then begin
    list.TieI=ptr_new(peaks)
    list.TieO=ptr_new(peakso)
    (*list.nrTie)[0]=n
endif else begin
    (*list.TieI)=[[(*list.TieI)],[peaks]]
    (*list.TieO)=[[(*list.TieO)],[peakso]]
    (*list.nrTie)[0]+=n
endelse

; Filter out close pixels
halfshift=1
TieImin=(*list.TieI)-halfshift
TieIMax=(*list.TieI)+halfshift
nrTie=(*list.nrTie)[0]
i=0
while i lt nrTie-1 do begin
    ind=where( ((*list.TieI)[0,*] ge TieImin[0,i]) and ((*list.TieI)[0,*] le TieIMax[0,i]) and $
               ((*list.TieI)[1,*] ge TieImin[1,i]) and ((*list.TieI)[1,*] le TieIMax[1,i]),ct )
    if ct gt 1 then begin
          this=where(ind ne i,ct)
          (*list.TieI)=ShrinkArray((*list.TieI),ind[this],/row,count=nrTie)
          (*list.TieO)=ShrinkArray((*list.TieO),ind[this],/row)
          TieImin=ShrinkArray(TieImin,ind[this],/row)
          TieImax=ShrinkArray(TieImax,ind[this],/row)
          ;ind2=[ind2,ind[this]]
    endif
    i++
endwhile
(*list.nrTie)[0]=nrTie

RefreshDisplay, ev.top,list
end;pro easy_grid
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeletePreviousAutoSet,list,call
ptr_free,list.smask
ptr_free,list.mask
list.masknr=0

if call eq 'Tilt' then begin
    list.nrTiltRing=0
    (*list.nrTilt)=[0]
    ptr_free,list.TiltI
    ptr_free,list.TiltRingD
endif else begin
    list.nrTieRing=0
    (*list.nrTie)=[0]
    ptr_free,list.TieI
    ptr_free,list.TieO
    ptr_free,list.TieRingD
endelse
end;pro DeletePreviousAutoSet
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AutoSetRing,mask,list,nSections

tstart=systime(1)

; ----(tt,azimuth)space----
res=BraggPtoX(list.AzIntInfo.MinPix,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
rmin=[res[1],res[0],list.AzIntInfo.MinPix,res[2]]
; User defined dr:
dr=[list.AzIntInfo.BWD,list.AzIntInfo.BWT,list.AzIntInfo.BWP,list.AzIntInfo.BWQ]/10
; Adapt dr to the current ring:
;dr=fltarr(4)

tmp=AzIntWarp(mask,list,list.IntCor2Dinfo,1,rmin,dr,list.AzIntInfo.outbin,list.AzIntInfo.AzWidth,nSections,$
            [list.AzIntInfo.Fill,list.AzIntInfo.bool],imagevalid=list.tiffvalid,median=list.AzIntInfo.median,$
            percentile=list.AzIntInfo.percentile,printres=list.OutID,$
            updatenorm=list.AzIntInfo.updatenorm,calcerror=list.AzIntInfo.calcerror,$
            xrange=tt,yrange=phi,oimage=oimage)

; ----Get [tt,azimuth] centers----
axbin=(size(oimage,/dim))[0]
aybin=(size(oimage,/dim))[1]
X=indgen(axbin)
sigma=2.
areamult=sqrt(2.*!pi)*sigma
ttmax=max(tt,min=ttmin)

;window ;PLOT

n_nonzero=axbin*0.9
tt_i=0.
phi_i=0.
for i=0l,aybin-1 do begin
    ind1=where(oimage[*,i] ne 0,count1)

    if count1 ge n_nonzero then begin
        ; Find center of reflection for this azimuth
        ma=max(oimage[*,i],cen,min=mi)
        cen=tt[cen[0]]
        tmp=where(oimage[*,i] ge (ma+mi)/2,ct)
        
        if ct gt 0 then begin
            ma2=max(tt[tmp],min=mi2)
            FWHM=ma2-mi2
            sigma=0.424661*FWHM
            height=ma
            area=height*sqrt(2.*!pi)*sigma

            ; Use Gaussian
            keep = [height,cen,sigma,mi,0]
            
            ;plot,tt,oimage[*,i],psym=1;PLOT
            ;gfunc_gauss,tt,keep,result;PLOT
            ;oplot,tt,result,color=100;PLOT
            
            ; Use IDL's gaussfit
            ;Result = GAUSSFIT( tt, oimage[*,i], keep ,nterms=5, ESTIMATES=keep,MEASURE_ERRORS=1/sqrt(oimage[*,i]))
            
            ; Use curvefit
            Result = CURVEFIT(tt,oimage[*,i],1/oimage[*,i],keep,FUNCTION_NAME = 'gfunc_gauss',iter=iter,status=status,chisq=chisq)
            
            ; Use NLLS (too slow)
            ;keep = [area,cen,sigma,mi,0]
            ;Result = NLLS(tt, oimage[*,i], keep, 'gfunc_gauss',CError=CError,weight=1)
            ;status = CError le 0
            
            ; Use Pearson VII (takes to long, constraint R<=1)
;            keep = [area,cen,FWHM,1,mi,0]
;            window
;            plot,tt,oimage[*,i],psym=3
;            gfunc_PearsonVII,tt,keep,result
;            oplot,tt,result,color=100
;            ;Result = CURVEFIT(tt,oimage[*,i],1/oimage[*,i],keep,FUNCTION_NAME = 'gfunc_PearsonVII')
;            ;Result = NLLS(tt, oimage[*,i], keep, 'gfunc_PearsonVII',CError=CError,weight=1)
            
            ;oplot,tt,result;PLOT
            ;wait,0.1;PLOT
            
            cen=keep[1]
            if cen le ttmax and cen ge ttmin and keep[0] gt 0 and keep[1] gt 0 and keep[2] gt 0 and status eq 0 and iter ne 0 then begin
                tt_i=[tt_i,cen]
                phi_i=[phi_i,phi[i]]
            endif
        endif
    endif
endfor

nSections=n_elements(tt_i)-1
if nSections eq 0 then out=0 $
else begin
    tt_i=tt_i[1:nSections]
    phi_i=phi_i[1:nSections]
endelse

; ----Convert [tt,azimuth] to [x,y]----
out=BraggXtoXY(tt_i,phi_i,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,1,/angledeg)

printw,list.OutID,'AutoSetRing processing time:'+string(systime(1)-tstart)
return,out
end;function AutoSetRing
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotoimage,list,PS=PS
LoadctRev,0,/silent,rev=list.revcol
soimage=size(list.z,/dim)
nmag=[!d.x_size,!d.y_size]*0.8/soimage
imgtvscl,list.z,0.05,0.1,nmag=nmag,xtitle=list.str,ytitle='Azimuth(degrees)',xx=list.x,yy=list.y,sclmax=list.sclmax,displaypower=list.displaypower
end;pro plotoimage
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_az_int, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=tlist
WIDGET_CONTROL, tlist.top, sensitive=1,get_uvalue=list,/no_copy
list.AzIntInfo=tlist.AzIntInfo
list.IntCor2Dinfo=tlist.IntCor2Dinfo
WIDGET_CONTROL, tlist.top,set_uvalue=list,/no_copy
end;pro CleanUp_az_int
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro az_int_event, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

evorg=widget_info(ev.id,/type)
if evorg eq 10 then return

; Main + local structure
widget_control,ev.id,get_value=val
widget_control,ev.top,get_uvalue=tlist
widget_control,tlist.top,get_uvalue=list

; Text widget
if evorg eq 3 then begin
    fval=0.
    reads,val,fval
    widget_control,ev.id,get_uvalue=uval
    case uval of
    'outbins':    tlist.AzIntInfo.outbin=fval
    'MaxD':     begin
                tlist.MaxD=fval
                res=BraggDtoX(tlist.MaxD,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
                tlist.Mintt=res[1]
                tlist.AzIntInfo.MinPix=res[0]
                tlist.MinQ=res[2]
                ID=widget_info(ev.top,FIND_BY_UNAME='Mintt')
                 widget_control,ID,set_value=stringr(tlist.Mintt*180/!dpi)
                 ID=widget_info(ev.top,FIND_BY_UNAME='MinPix')
                 widget_control,ID,set_value=stringr(tlist.AzIntInfo.MinPix)
                 ID=widget_info(ev.top,FIND_BY_UNAME='MinQ')
                 widget_control,ID,set_value=stringr(tlist.MinQ)
                endcase
    'Mintt':     begin
                tlist.Mintt=fval/180.*!dpi
                res=BraggTtoX(tlist.Mintt,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
                tlist.MaxD=res[1]
                tlist.AzIntInfo.MinPix=res[0]
                tlist.MinQ=res[2]
                ID=widget_info(ev.top,FIND_BY_UNAME='MaxD')
                 widget_control,ID,set_value=stringr(tlist.MaxD)
                 ID=widget_info(ev.top,FIND_BY_UNAME='MinPix')
                 widget_control,ID,set_value=stringr(tlist.AzIntInfo.MinPix)
                 ID=widget_info(ev.top,FIND_BY_UNAME='MinQ')
                 widget_control,ID,set_value=stringr(tlist.MinQ)
                endcase
    'MinPix':     begin
                tlist.AzIntInfo.MinPix=fval
                res=BraggPtoX(tlist.AzIntInfo.MinPix,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
                tlist.Mintt=res[0]
                 tlist.MaxD=res[1]
                 tlist.MinQ=res[2]
                 ID=widget_info(ev.top,FIND_BY_UNAME='MaxD')
                 widget_control,ID,set_value=stringr(tlist.MaxD)
                 ID=widget_info(ev.top,FIND_BY_UNAME='Mintt')
                 widget_control,ID,set_value=stringr(tlist.Mintt*180/!dpi)
                 ID=widget_info(ev.top,FIND_BY_UNAME='MinQ')
                 widget_control,ID,set_value=stringr(tlist.MinQ)
                endcase
    'MinQ':     begin
                tlist.MinQ=fval
                res=BraggQtoX(tlist.MinQ,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
                tlist.AzIntInfo.MinPix=res[0]
                tlist.Mintt=res[1]
                 tlist.MaxD=res[2]
                 ID=widget_info(ev.top,FIND_BY_UNAME='MaxD')
                 widget_control,ID,set_value=stringr(tlist.MaxD)
                 ID=widget_info(ev.top,FIND_BY_UNAME='Mintt')
                 widget_control,ID,set_value=stringr(tlist.Mintt*180/!dpi)
                 ID=widget_info(ev.top,FIND_BY_UNAME='MinPix')
                 widget_control,ID,set_value=stringr(tlist.AzIntInfo.MinPix)
                endcase
    'BWD':         begin
                tlist.AzIntInfo.BWD=fval
                endcase
    'BWT':         begin
                tlist.AzIntInfo.BWT=fval/180*!dpi
                endcase
    'BWP':         begin
                tlist.AzIntInfo.BWP=fval
                endcase
    'BWQ':         begin
                tlist.AzIntInfo.BWQ=fval
                endcase
    'ChangeNameText':begin
                tlist.AzIntInfo.AzWidth=fval/180*!dpi
                endcase
    'ChangeSectors': tlist.AzIntInfo.nsectors=round(fval)>0
    'IntCor2Dinfo.Plindeg': tlist.IntCor2Dinfo.Plindeg=fval
    'IntCor2Dinfo.n': tlist.IntCor2Dinfo.n=fval
    'IntCor2Dinfo.mmult': tlist.IntCor2Dinfo.mmult=fval
    'IntCor2Dinfo.dmono': tlist.IntCor2Dinfo.dmono=fval
    'IntCor2Dinfo.flatnormphi': tlist.IntCor2Dinfo.flatnormphi=fval/180*!dpi
    'IntCor2Dinfo.flatnormpolar': tlist.IntCor2Dinfo.flatnormpolar=fval/180*!dpi
    'IntCor2Dinfo.flatD': tlist.IntCor2Dinfo.flatD=fval
    'IntCor2Dinfo.muL': tlist.IntCor2Dinfo.muL=fval
    'percentile': tlist.AzIntInfo.percentile=fval
    endcase
    widget_control,ev.top,set_uvalue=tlist
    return
endif

; droplis widget
if evorg eq 8 then begin
    widget_control,ev.id,get_uvalue=uval
    case uval of
    'IntCor2Dinfo.Ptype':     begin
                            tlist.IntCor2Dinfo.Ptype=ev.index
                            tlist.IntCor2Dinfo.Ltype=tlist.IntCor2Dinfo.Ptype ne 0
                            ID=widget_info(ev.top,find_by_uname='basePnone')
                            widget_control,ID,sensitive=tlist.IntCor2Dinfo.Ptype ne 0
                            ID=widget_info(ev.top,find_by_uname='basePmono')
                            widget_control,ID,sensitive=tlist.IntCor2Dinfo.Ptype eq 2
                            endcase
    'IntCor2Dinfo.s': tlist.IntCor2Dinfo.s=ev.index?-1:1
    'IntCor2Dinfo.Atype':     begin
                            tlist.IntCor2Dinfo.Atype=ev.index
                            ID=widget_info(ev.top,find_by_uname='baseFlat')
                            widget_control,ID,sensitive=tlist.IntCor2Dinfo.Atype eq 1
                            endcase
    endcase
    widget_control,ev.top,set_uvalue=tlist
    return
endif
    
; Button widget
case val of
    'Correct for continuous part of the Lorentz factor (1/sin(theta))': begin
         tlist.IntCor2Dinfo.Ltype=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
    'Exit' :begin
         widget_control,ev.top,/destroy
         return
         endcase
'D-spacing ('+msymbols('angstroms')+')':begin
         tlist.AzIntInfo.type=0
         widget_control,ev.top,set_uvalue=tlist
         ID=widget_info(ev.top,FIND_BY_UNAME='BWP')
         widget_control,ID,sensitive=0
         ID=widget_info(ev.top,FIND_BY_UNAME='BWT')
         widget_control,ID,sensitive=0
         ID=widget_info(ev.top,FIND_BY_UNAME='BWD')
         widget_control,ID,sensitive=1
         ID=widget_info(ev.top,FIND_BY_UNAME='BWQ')
         widget_control,ID,sensitive=0
         return
         endcase
'2-theta (deg)':begin
         tlist.AzIntInfo.type=1
         widget_control,ev.top,set_uvalue=tlist
         ID=widget_info(ev.top,FIND_BY_UNAME='BWP')
         widget_control,ID,sensitive=0
         ID=widget_info(ev.top,FIND_BY_UNAME='BWT')
         widget_control,ID,sensitive=1
         ID=widget_info(ev.top,FIND_BY_UNAME='BWD')
         widget_control,ID,sensitive=0
         ID=widget_info(ev.top,FIND_BY_UNAME='BWQ')
         widget_control,ID,sensitive=0
         return
         endcase
'Radial (mm)':begin
         tlist.AzIntInfo.type=2
         widget_control,ev.top,set_uvalue=tlist
         ID=widget_info(ev.top,FIND_BY_UNAME='BWP')
         widget_control,ID,sensitive=1
         ID=widget_info(ev.top,FIND_BY_UNAME='BWT')
         widget_control,ID,sensitive=0
         ID=widget_info(ev.top,FIND_BY_UNAME='BWD')
         widget_control,ID,sensitive=0
         ID=widget_info(ev.top,FIND_BY_UNAME='BWQ')
         widget_control,ID,sensitive=0
         return
         endcase
'Q (1/nm)':begin
         tlist.AzIntInfo.type=3
         widget_control,ev.top,set_uvalue=tlist
         ID=widget_info(ev.top,FIND_BY_UNAME='BWP')
         widget_control,ID,sensitive=0
         ID=widget_info(ev.top,FIND_BY_UNAME='BWT')
         widget_control,ID,sensitive=0
         ID=widget_info(ev.top,FIND_BY_UNAME='BWD')
         widget_control,ID,sensitive=0
         ID=widget_info(ev.top,FIND_BY_UNAME='BWQ')
         widget_control,ID,sensitive=1
         return
         endcase
'Nearest Neighbor':begin
         tlist.AzIntInfo.bool=0
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Bilinear': begin
         tlist.AzIntInfo.bool=1
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Cubic' :   begin
         tlist.AzIntInfo.bool=2
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Default OutputBins':begin
         result=ReadINI('$AzInt:')
         tlist.AzIntInfo.outbin=result.(2)
         ID=widget_info(ev.top,FIND_BY_UNAME='outbins')
         widget_control,ID,set_value=stringr(tlist.AzIntInfo.outbin)
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Default Min Inner Radius':begin
         result=ReadINI('$AzInt:')
         tlist.AzIntInfo.MinPix=result.(1)[0]
         res=BraggPtoX(tlist.AzIntInfo.MinPix,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
         tlist.Mintt=res[0]
         tlist.MaxD=res[1]
         tlist.MinQ=res[2]
         ID=widget_info(ev.top,FIND_BY_UNAME='MaxD')
         widget_control,ID,set_value=stringr(tlist.MaxD)
         ID=widget_info(ev.top,FIND_BY_UNAME='MinQ')
         widget_control,ID,set_value=stringr(tlist.MinQ)
         ID=widget_info(ev.top,FIND_BY_UNAME='Mintt')
         widget_control,ID,set_value=stringr(tlist.Mintt*180/!dpi)
         ID=widget_info(ev.top,FIND_BY_UNAME='MinPix')
         widget_control,ID,set_value=stringr(tlist.AzIntInfo.MinPix)
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Default Binwidth':begin
         result=ReadINI('$AzInt:')
         tlist.AzIntInfo.BWP=result.(1)[1]
         tlist.AzIntInfo.BWT=result.(1)[2]
         tlist.AzIntInfo.BWD=result.(1)[3]
         tlist.AzIntInfo.BWQ=result.(1)[4]
         ID=widget_info(ev.top,FIND_BY_UNAME='BWD')
         widget_control,ID,set_value=stringr(tlist.AzIntInfo.BWD)
         ID=widget_info(ev.top,FIND_BY_UNAME='BWT')
         widget_control,ID,set_value=stringr(tlist.AzIntInfo.BWT*180/!dpi)
         ID=widget_info(ev.top,FIND_BY_UNAME='BWP')
         widget_control,ID,set_value=stringr(tlist.AzIntInfo.BWP)
         ID=widget_info(ev.top,FIND_BY_UNAME='BWQ')
         widget_control,ID,set_value=stringr(tlist.AzIntInfo.BWQ)
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Default Azimuthal binwidth':begin
         result=ReadINI('$AzInt:')
         tlist.AzIntInfo.AzWidth=result.(1)[5]
         textin=stringr(tlist.AzIntInfo.AzWidth*180/!dpi)
         ID=widget_info(ev.top,FIND_BY_UNAME='ChangeNameText')
         widget_control,ID,set_value=textin
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Default Sectors':begin
         result=ReadINI('$AzInt:')
         tlist.AzIntInfo.nsectors=result.(5)
         textin=stringr(tlist.AzIntInfo.nsectors)
         ID=widget_info(ev.top,FIND_BY_UNAME='ChangeSectors')
         widget_control,ID,set_value=textin
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Show 1D scan':begin
         tlist.AzIntInfo.oD=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Open in 1D Editor':begin
         tlist.AzIntInfo.loD=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Automatic integration on load':begin
         tlist.AzIntInfo.aoD=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Save 1D scan':begin
         tlist.AzIntInfo.soD=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Save Azimuth vs. Radial image':begin
         tlist.AzIntInfo.stD=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Fill pattern':  begin
         tlist.AzIntInfo.Fill=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Use median instead of average':  begin
         tlist.AzIntInfo.median=ev.select
         widget_control,ev.top,set_uvalue=tlist
         ID=widget_info(ev.top,FIND_BY_UNAME='percentile')
         widget_control,ID,sensitive=ev.select
         return
         endcase
'Error propagation':  begin
         tlist.AzIntInfo.calcerror=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Changing bad pixels':  begin
         tlist.AzIntInfo.updatenorm=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Batch Integration':begin
         tlist.AzIntInfo.Batch=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
'Save with ring percentage':begin
         tlist.AzIntInfo.SaveRingPercentage=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
         endcase
  'Save original CHI':begin
           tlist.AzIntInfo.origformat=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
           endcase
  'Save in tiff file':begin
           tlist.AzIntInfo.tiff=ev.select
         widget_control,ev.top,set_uvalue=tlist
         return
           endcase
  else:  begin
           widget_control, /hourglass
         endcase
endcase

; ----Prepare Integrate----
s=*list.tiffs
rmin=[tlist.MaxD,tlist.Mintt,tlist.AzIntInfo.MinPix,tlist.MinQ]
dr=[tlist.AzIntInfo.BWD,tlist.AzIntInfo.BWT,tlist.AzIntInfo.BWP,tlist.AzIntInfo.BWQ]
;dr=0.5 pixels  ;oversampling of two
;=> twice as many output bins of half the radial width

mask=(*list.mask)[*,val]

if tlist.AzIntInfo.oD or tlist.AzIntInfo.stD then begin
    oDindex=winid()
    Window,oDindex
endif

if tlist.AzIntInfo.soD or tlist.AzIntInfo.stD or tlist.AzIntInfo.loD then begin
    path=CutPath(list.file,file=file)

    if tlist.AzIntInfo.soD then begin
        filter=ChiFormat(/both)
        path=SelFile(list.chipath,file+strmid(filter[0],1),filter,'Save 1D scan....')
        if path eq '' then return
        path=CutPath(path,path=chipath,file=file,ext=ext)
        if ext eq '' then ext='.chi'
        list.chipath=chipath
        list.chifile=file+ext
    endif
    
    if tlist.AzIntInfo.stD then begin
        path=SelFile(list.chipath,file+'.tiff','*.tiff','Save 2D scan....')
        if path eq '' then return
        path=CutPath(path,path=oimagepath,file=file,ext=ext2)
        if ext2 eq '' then ext2='.tiff'
    endif
    
    widget_control,tlist.top,set_uvalue=list
    case mask[0] of
    2:    begin
        ROIAzimuth=mask[[4,6]]
        endcase
    3:    begin
        ROIAzimuth=[0,2*!pi]
        endcase
    endcase
    str=ChiXtitle(tlist.AzIntInfo.type)
endif

if tlist.AzIntInfo.Batch then begin
    pathlist=Select_Files(ev.top,list)
    if pathlist[0] eq '' then return

    n=n_elements(pathlist)
    Result = DIALOG_MESSAGE(stringr(n)+' files for azimuthal integration: proceed?',/question)
    if result eq 'No' then return
    pathlist=['this',pathlist]
    n++
endif else begin
    pathlist='this'
    n=1
endelse

; ----Perform Integrate----
if tlist.AzIntInfo.Batch then begin
    progressBar = Obj_New("PROGRESSBAR",title="Azimuthally integrating files...")
    progressBar -> Start
    ReadWithAutoStart,list
endif

printw,list.OutID,'Azimuthal integration:'
i0=tlist.AzIntInfo.Batch ; don't process 'this' when in batch mode
batchreturnp=ptr_new()
for i=i0, n-1 do begin

    if tlist.AzIntInfo.Batch then $
    IF progressBar -> CheckCancel() THEN BEGIN
        progressBar -> Destroy
        ReadWithAutoStop,list
        printw,list.outid,'Process CANCELLED!'
        RETURN
    ENDIF

    ; Load CCD image + correct
    if i ne 0 then ok=ReadWithAutoCorrect(pathlist[i],list,file=file) else ok=1b

    ; integrate
    if ok then begin
        
        if ptr_valid(batchreturnp) then begin
            d=AzIntWarpBatch(batchreturnp,list.tiff,list.tiffvalid,oimage=oimage)
            bring=0b
        endif else begin
            d=AzIntWarp(mask,list,tlist.IntCor2Dinfo,tlist.AzIntInfo.type,rmin,dr,tlist.AzIntInfo.outbin,tlist.AzIntInfo.AzWidth,0,$
                [tlist.AzIntInfo.Fill,tlist.AzIntInfo.bool],imagevalid=list.tiffvalid,median=tlist.AzIntInfo.median,percentile=tlist.AzIntInfo.percentile,$
                printres=list.OutID,$
                listsaverp={wb:list.wb,niter:list.niter,thres:list.thres,done:list.CorDone[((*(list.CorSeq))).back]},saverp=tlist.AzIntInfo.SaveRingPercentage,$
                batchreturnp=batchreturnp,$
                splitnsectors=tlist.AzIntInfo.nsectors,$
                updatenorm=tlist.AzIntInfo.updatenorm,calcerror=tlist.AzIntInfo.calcerror,$
                oimage=oimage,xrange=xrange,yrange=yrange)
            bring=tlist.AzIntInfo.SaveRingPercentage
        endelse
        
        ; Plot 1D scan
        if tlist.AzIntInfo.oD then begin
            wset,oDindex
            loadctrev,-39
            plot,d[0,*],d[1,*],xtitle=str,ytitle='Intensity (a.u.)',psym=-2,xstyle=1,ystyle=1
        endif
        
        ; Save 1D scan
        if tlist.AzIntInfo.soD then begin
            temp=WriteChi(chipath+file+ext,d,tlist.AzIntInfo.type,list.path+list.file,list.OutID,ROIAzimuth,$
                            origformat=tlist.AzIntInfo.origformat,nsectors=tlist.AzIntInfo.nsectors,$
                            bring=bring,berror=tlist.AzIntInfo.calcerror,dist=list.dist)
        endif
        
        if tlist.AzIntInfo.loD then begin
            temp=min(d[1,*],y1)
            temp=max(d[1,*],y2)
            if tlist.AzIntInfo.calcerror then speerror=reform(d[2,*]) else speerror=reform(sqrt(d[1,*]))
            fileptr={xval:reform(d[0,*]),spe:reform(d[1,*]),speerror:speerror,ROIAzimuth:ROIAzimuth,ytitle:'Intensity (a.u.)',$
                xtitle:ChiXtitle(tlist.AzIntInfo.type),title:file,xtype:tlist.AzIntInfo.type,$
                xrange:[0,n_elements(d)/2-1],yrange:[y1[0],y2[0]]}
    
            MainUpdateCManual,'c:\dummy.chi',tlist.top,ev,fileptr=fileptr
        endif
        
        ; Save 2D image
        if tlist.AzIntInfo.stD then begin
            result=WriteChi(oimagepath[0]+file+ext2,[[reform(xrange)],[oimage]],tlist.AzIntInfo.type,'...',list.OutID,[yrange[0],yrange[n_elements(yrange)-1]])
            
;            if i ne i0 then plotoimagestruc.z=oimage $
;            else plotoimagestruc={z:oimage,str:str,x:reform(xrange),y:reform(yrange)*180/!pi,sclmax:list.sclmax,displaypower:list.displaypower,revcol:list.revcol}
;            saveimage,oimagepath[0],file+ext2,'plotoimage',plotoimagestruc,$
;                /noprompt,/forcereplot,outid=list.OutID
        endif
        
    endif else printw,list.outid,'Something wrong with '+pathlist[i]
    
    loadctrev,0,rev=-1
    if tlist.AzIntInfo.Batch then progressBar -> Update, (i+1.)/n*100
endfor

if tlist.AzIntInfo.Batch then begin
    progressBar -> Destroy
    ReadWithAutoStop,list
endif

if size(batchreturnp,/type) eq 10 then begin
    widget_control,tlist.top,get_uvalue=list
    heap_free,list.batchreturnp
    list.batchreturnp=batchreturnp
    widget_control,tlist.top,set_uvalue=list
endif

end;pro az_int_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro az_int,ev

widget_control,ev.top,get_uvalue=list

WIDGET_CONTROL, ev.top, sensitive=0
s=*list.tiffs

res=BraggPtoX(list.AzIntInfo.MinPix,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
tlist={top:ev.top,Mintt:res[0],MaxD:res[1],MinQ:res[2],AzIntInfo:list.AzIntInfo,IntCor2Dinfo:list.IntCor2Dinfo}

base=widget_base(/column,title='Azimuthal integration',xoffset=200,yoffset=200,uvalue=tlist)

; ----Mask buttons----
if list.masknr ne 0 then begin
    label=widget_label(base,value='Select area to integrate:',/ALIGN_LEFT)
    nr=10
    nbase=ceil(float(list.masknr)/nr)
    base1=lonarr(nbase)
    for i=0l,nbase-1 do base1[i]=widget_base(base,/row)
    button = lonarr(list.masknr)
    j=0
    for i=0l,list.masknr-1 do begin
        button[i]=widget_button(base1[j],value=stringr(i))
        if (*list.mask)[0,i] eq 1 then widget_control,button[i],sensitive=0
        if ((i+1) mod nr) eq 0 then j=j+1
    endfor
endif

; ----Type Buttons----
tab=widget_tab(base)
basea=widget_base(tab,/column,title='Basic')
basee=widget_base(tab,/column,title='Speed')
baseb=widget_base(tab,/column,title='Advanced')
basec=widget_base(tab,/column,title='Polarization')
based=widget_base(tab,/column,title='Attenuation')

; ----Basic----
base2=widget_base(basea,/row)
base22=widget_base(base2,/column)
base222=widget_base(base22,/column,/exclusive)
button30=widget_button(base222,value='D-spacing ('+msymbols('angstroms')+')')
button31=widget_button(base222,value='2-theta (deg)')
button32=widget_button(base222,value='Radial (mm)')
button33=widget_button(base222,value='Q (1/nm)')

base5=widget_base(base22,/column,/nonexclusive)
button1=widget_button(base5,value='Show 1D scan')
button8=widget_button(base5,value='Open in 1D Editor')
button2=widget_button(base5,value='Save 1D scan')
button6b=widget_button(base5,value='Batch Integration')

base2b=widget_base(base2,/column)
button=widget_button(base2b,value='Default Binwidth')
baset=widget_base(base2b,/row)
    text=widget_text(baset,value=stringr(tlist.AzIntInfo.BWD),/editable,uname='BWD',uvalue='BWD',sensitive=tlist.AzIntInfo.type eq 0)
    label=widget_label(baset,value=msymbols('angstroms'))
baset=widget_base(base2b,/row)
    text=widget_text(baset,value=stringr(tlist.AzIntInfo.BWT*180/!dpi),/editable,uname='BWT',uvalue='BWT',sensitive=tlist.AzIntInfo.type eq 1)
    label=widget_label(baset,value='deg')
baset=widget_base(base2b,/row)
    text=widget_text(baset,value=stringr(tlist.AzIntInfo.BWP),/editable,uname='BWP',uvalue='BWP',sensitive=tlist.AzIntInfo.type eq 2)
    label=widget_label(baset,value='flat pixels')
baset=widget_base(base2b,/row)
    text=widget_text(baset,value=stringr(tlist.AzIntInfo.BWQ),/editable,uname='BWQ',uvalue='BWQ',sensitive=tlist.AzIntInfo.type eq 3)
    label=widget_label(baset,value='1/nm')
    
base2b=widget_base(base2b,/column)
button=widget_button(base2b,value='Default OutputBins')
text=widget_text(base2b,value=stringr(tlist.AzIntInfo.outbin),/editable,uname='outbins',uvalue='outbins')

; ----Speed----
base4=widget_base(basee,/row)
basett=widget_base(base4,/row)
    button=widget_button(basett,value='Default Azimuthal binwidth')
    baset=widget_base(basett,/row)
        text=widget_text(baset,value=stringr(tlist.AzIntInfo.AzWidth*180/!dpi),/editable,uname='ChangeNameText',uvalue='ChangeNameText')
        label=widget_label(baset,value='deg',uname='ChangeNameLabel')
baseoptions=widget_base(basee,/row)
base5=widget_base(baseoptions,/column,/nonexclusive)
button6fb=widget_button(base5,value='Error propagation')
button6fc=widget_button(base5,value='Changing bad pixels')

; ----Advanced----
base3=widget_base(baseb,/row)

base2a=widget_base(base3,/column)
button=widget_button(base2a,value='Default Min Inner Radius')
baset=widget_base(base2a,/row)
    text=widget_text(baset,value=stringr(tlist.MaxD),/editable,uname='MaxD',uvalue='MaxD')
    label=widget_label(baset,value=msymbols('angstroms'))
baset=widget_base(base2a,/row)
    text=widget_text(baset,value=stringr(tlist.Mintt*180/!dpi),/editable,uname='Mintt',uvalue='Mintt')
    label=widget_label(baset,value='deg')
baset=widget_base(base2a,/row)
    text=widget_text(baset,value=stringr(tlist.AzIntInfo.MinPix),/editable,uname='MinPix',uvalue='MinPix')
    label=widget_label(baset,value='flat pixels')
baset=widget_base(base2a,/row)
    text=widget_text(baset,value=stringr(tlist.MinQ),/editable,uname='MinQ',uvalue='MinQ')
    label=widget_label(baset,value='1/nm')
    
basett=widget_base(base3,/column)
    button=widget_button(basett,value='Default Sectors')
    baset=widget_base(basett,/row)
        text=widget_text(baset,value=stringr(tlist.AzIntInfo.nsectors),/editable,uname='ChangeSectors',uvalue='ChangeSectors')
    
baseoptions=widget_base(baseb,/row)
base5=widget_base(baseoptions,/column,/nonexclusive)
button6fa=widget_button(base5,value='Use median instead of average')
button6=widget_button(base5,value='Fill pattern')
button6c=widget_button(base5,value='Save with ring percentage')
button6d=widget_button(base5,value='Save original CHI')
button2b=widget_button(base5,value='Save Azimuth vs. Radial image')
;button6e=widget_button(base5,value='Save in tiff file')
button7=widget_button(base5,value='Automatic integration on load')

base4bb=widget_base(baseoptions,/column)
base4b=widget_base(base4bb,/column)
    label=widget_label(base4b,value='Percentile (%):')
    textpercentile=widget_text(base4b,value=stringr(tlist.AzIntInfo.percentile),/editable,uname='percentile',uvalue='percentile')
base4b=widget_base(base4bb,/column,/exclusive)
    button40=widget_button(base4b,value='Nearest Neighbor')
    button41=widget_button(base4b,value='Bilinear')
    button42=widget_button(base4b,value='Cubic')


; ----Corrections----
xs=200

baset=widget_base(basec,/row)
    label=widget_label(baset,value='Type of polarization correction',xsize=xs)
    dropPtype=widget_droplist(baset,value=['None','Without optics','With monochromator'],uvalue='IntCor2Dinfo.Ptype')

basect=widget_base(basec,/column,sensitive=tlist.IntCor2Dinfo.Ptype ne 0,uname='basePnone')

baset=widget_base(basect,/row)
    label=widget_label(baset,value='Linear degree of polarization (horizontal)',xsize=xs)
    text=widget_text(baset,value=stringr(tlist.IntCor2Dinfo.Plindeg),/editable,uvalue='IntCor2Dinfo.Plindeg')
    
basect=widget_base(basec,/column,sensitive=tlist.IntCor2Dinfo.Ptype eq 2,uname='basePmono')

baset=widget_base(basect,/row)
    label=widget_label(baset,value='Number of monochromator crystals',xsize=xs)
    text=widget_text(baset,value=stringr(tlist.IntCor2Dinfo.n),/editable,uvalue='IntCor2Dinfo.n')
    
baset=widget_base(basect,/row)
    label=widget_label(baset,value='Mosaic (2) or perfect (1) crystals',xsize=xs)
    text=widget_text(baset,value=stringr(tlist.IntCor2Dinfo.mmult),/editable,uvalue='IntCor2Dinfo.mmult')

baset=widget_base(basect,/row)
    label=widget_label(baset,value='Monochromator crystal d-spacing ('+msymbols('angstroms')+')',xsize=xs)
    text=widget_text(baset,value=stringr(tlist.IntCor2Dinfo.dmono),/editable,uvalue='IntCor2Dinfo.dmono')
    
baset=widget_base(basect,/row)
    label=widget_label(baset,value='Monochromator orientation',xsize=xs)
    drops=widget_droplist(baset,value=['Vertical','Horizontal'],uvalue='IntCor2Dinfo.s')

;baset=widget_base(basec,/row,/nonexclusive)
;    button9=widget_button(baset,value='Correct for continuous part of the Lorentz factor (1/sin(theta))')
    
baset=widget_base(based,/row)
    label=widget_label(baset,value='Type of attenuation correction',xsize=xs)
    dropAtype=widget_droplist(baset,value=['None','Flat plate'],uvalue='IntCor2Dinfo.Atype')
    
basedt=widget_base(based,/column,sensitive=tlist.IntCor2Dinfo.Atype eq 1,uname='baseFlat')

baset=widget_base(basedt,/row)
    label=widget_label(baset,value='Surface normal polar angle (deg)',xsize=xs)
    text=widget_text(baset,value=stringr(tlist.IntCor2Dinfo.flatnormpolar*180/!pi),/editable,uvalue='IntCor2Dinfo.flatnormpolar')

baset=widget_base(basedt,/row)
    label=widget_label(baset,value='Surface normal azimuthal angle (deg)',xsize=xs)
    text=widget_text(baset,value=stringr(tlist.IntCor2Dinfo.flatnormphi*180/!pi),/editable,uvalue='IntCor2Dinfo.flatnormphi')

baset=widget_base(basedt,/row)
    label=widget_label(baset,value='Thickness (mm)',xsize=xs)
    text=widget_text(baset,value=stringr(tlist.IntCor2Dinfo.flatD),/editable,uvalue='IntCor2Dinfo.flatD')

baset=widget_base(basedt,/row)
    label=widget_label(baset,value='Linear attenuation coefficient (cm^(-1))',xsize=xs)
    text=widget_text(baset,value=stringr(tlist.IntCor2Dinfo.muL),/editable,uvalue='IntCor2Dinfo.muL')

releasebutton=widget_button(base,value='Exit',uname='Exit')


WIDGET_CONTROL, base, /REALIZE
widget_control,button1,set_button=tlist.AzIntInfo.oD
widget_control,button2,set_button=tlist.AzIntInfo.soD
widget_control,button2b,set_button=tlist.AzIntInfo.stD
widget_control,button7,set_button=tlist.AzIntInfo.aoD
widget_control,button8,set_button=tlist.AzIntInfo.loD
case tlist.AzIntInfo.type of
0:widget_control,button30,set_button=1
1:widget_control,button31,set_button=1
2:widget_control,button32,set_button=1
3:widget_control,button33,set_button=1
endcase
case tlist.AzIntInfo.bool of
0:widget_control,button40,set_button=1
1:widget_control,button41,set_button=1
2:widget_control,button42,set_button=1
endcase
widget_control,button6,set_button=tlist.AzIntInfo.Fill
widget_control,button6fa,set_button=tlist.AzIntInfo.median
widget_control,button6fb,set_button=tlist.AzIntInfo.calcerror
widget_control,button6fc,set_button=tlist.AzIntInfo.updatenorm
widget_control,textpercentile,sensitive=tlist.AzIntInfo.median
widget_control,button6b,set_button=tlist.AzIntInfo.Batch
widget_control,button6c,set_button=tlist.AzIntInfo.SaveRingPercentage
widget_control,button6d,set_button=tlist.AzIntInfo.origformat
;widget_control,button9,set_button=tlist.IntCor2Dinfo.Ltype
;widget_control,button6e,set_button=tlist.AzIntInfo.tiff
widget_control,dropPtype,SET_DROPLIST_SELECT=list.IntCor2Dinfo.Ptype
widget_control,drops,SET_DROPLIST_SELECT=list.IntCor2Dinfo.s eq -1
widget_control,dropAtype,SET_DROPLIST_SELECT=list.IntCor2Dinfo.Atype

Xmanager,'XRDUA_event',base, event_handler='az_int_event',$
    cleanup='CleanUp_az_int',GROUP_LEADER=ev.top

end;pro az_int
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_auto_set, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
WIDGET_CONTROL, list.top, sensitive=1
end;pro CleanUp_auto_set
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddTi,TiI,nSections,Call,list,top,d
if nSections eq 0 then return
TiI=transpose(TiI)
if Call eq 'Tilt' then begin
    if ptr_valid(list.TiltI) then (*list.TiltI)=[[(*list.TiltI)],[TiI]] else list.TiltI=ptr_new(TiI)
    if ptr_valid(list.TiltRingD) then (*list.TiltRingD)=[(*list.TiltRingD),0.] else list.TiltRingD=ptr_new([0.])
    (*list.nrTilt)=[(*list.nrTilt),0l]
    if (*list.nrTilt)[list.nrTiltRing] eq 0 then begin
       (*list.nrTilt)[list.nrTiltRing]=nSections
       if n_elements(d) eq 1 then (*list.TiltRingD)[list.nrTiltRing]=d
    endif else begin
       (*list.nrTilt)[list.nrTiltRing+1]=nSections
       if n_elements(d) eq 1 then (*list.TiltRingD)[list.nrTiltRing+1]=d
    endelse
    list.nrTiltRing++
endif else begin
    if ptr_valid(list.TieI) then (*list.TieI)=[[(*list.TieI)],[TiI]] else list.TieI=ptr_new(TiI)
    if ptr_valid(list.TieO) then (*list.TieO)=[[(*list.TieO)],[TiI]] else list.TieO=ptr_new(TiI)
    if ptr_valid(list.TieRingD) then (*list.TieRingD)=[(*list.TieRingD),0.] else list.TieRingD=ptr_new([0.])
    (*list.nrTie)=[(*list.nrTie),0l]
    if (*list.nrTie)[list.nrTieRing] eq 0 then begin
       (*list.nrTie)[list.nrTieRing]=nSections
       if n_elements(d) eq 1 then (*list.TieRingD)[list.nrTieRing]=d
    endif else begin
       (*list.nrTie)[list.nrTieRing+1]=nSections
       if n_elements(d) eq 1 then (*list.TieRingD)[list.nrTieRing+1]=d
    endelse
    list.nrTieRing++
endelse
RefreshDisplay,top,list
end;pro AddTi
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro auto_set_event, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_value=val
if val[0] eq 'Exit' then begin
    widget_control,ev.top,/destroy
    return
endif

widget_control,ev.top,get_uvalue=tlist
nSections=tlist.nSections
widget_control,tlist.top,get_uvalue=list

if widget_info(ev.id,/type) EQ 3 THEN BEGIN
    reads,val,nSections
    tlist.(1)=nSections
    widget_control,ev.top,set_uvalue=tlist
    return
endif

mask=(*list.mask)[*,val]
TiI=AutoSetRing(mask,list,nSections)
AddTi,TiI,nSections,tlist.Call,list,tlist.top
end;pro auto_set_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro auto_set,ev,Call
widget_control,ev.top,get_uvalue=list
if (list.masknr eq 0) then return
WIDGET_CONTROL, ev.top, sensitive=0
base=widget_base(/row,title='Select Ring',xoffset=200,yoffset=200,$
                    uvalue={top:ev.top,nSections:300l,Call:Call})
button = lonarr(list.masknr)
for i=0l,list.masknr-1 do begin
    button[i]=widget_button(base,value=stringr(i))
    if ((*list.mask)[0,i] ne 2) and ((*list.mask)[0,i] ne 3) then widget_control,button[i],sensitive=0
endfor
label=widget_label(base,value='Sections:')
text=widget_text(base,value='300',/editable)
releasebutton=widget_button(base,value='Exit',uname='Exit')
WIDGET_CONTROL, base, /REALIZE
Xmanager,'auto_set',base, event_handler='auto_set_event',$
          cleanup='CleanUp_auto_set',GROUP_LEADER=ev.top
end;pro auto_set
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetFullArc,list

s=(*list.tiffs)[1:2]

centeroffimage=list.center[0] lt 0 or list.center[0] gt s[0]-1 or $
   list.center[1] lt 0 or list.center[1] gt s[1]-1

; d-spacings
x=[lindgen(s[0]),make_array(s[1],value=s[0]-1),lindgen(s[0]),lonarr(s[1])]
y=[lonarr(s[0]),lindgen(s[1]),make_array(s[0],value=s[1]-1),lindgen(s[1])]
d=BraggXYtoX(x,y,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/onlyd,aziout=aziout)

; NaN around
res=BraggPtoX(list.AzIntInfo.MinPix,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
ind=where(d gt res[1],ct)
if ct ne 0 then d[ind]=!values.F_NAN

; Azimuth
mi=min(aziout,max=ma)
if ma-mi gt 0.995*!pi then begin
    if centeroffimage then begin
        ind=where(aziout lt 0,ct)
        if ct ne 0 then aziout[ind]+=2*!pi
        mi=min(aziout,max=ma)
        aziout=[mi,ma]
    endif else delvar2,aziout
endif else aziout=[mi,ma]

mi=min(d,max=ma,/nan)
poff=list.AzIntInfo.MinPix
if ~centeroffimage then $
    ma=BraggXYtoX(list.center[0]+poff,list.center[1]+poff,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/onlyd)

if n_elements(aziout) eq 0 then return,[3,0,0,ma,mi,0,0] else return,[2,0,0,ma,aziout[0],mi,aziout[1]]
end;function GetFullArc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DisablePDFrings,index,list,aziout=aziout

; ----Get visible d-spacings of phase "index"----
if index ne 0 then i0=total((*list.PDFn)[0:index-1]) else i0=0
i1=total((*list.PDFn)[0:index])-1
ind=where((*list.PDFs)[i0:i1],ct)+i0
if ct eq 0 then return

; ----Disable d-spacings that fall of the image----
s=(*list.tiffs)[1:2]

centeroffimage=list.center[0] lt 0 or list.center[0] gt s[0]-1 or $
   list.center[1] lt 0 or list.center[1] gt s[1]-1

x=[lindgen(s[0]),make_array(s[1],value=s[0]-1),lindgen(s[0]),lonarr(s[1])]
y=[lonarr(s[0]),lindgen(s[1]),make_array(s[0],value=s[1]-1),lindgen(s[1])]
d=BraggXYtoX(x,y,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/onlyd,aziout=aziout)
mi=min(aziout,max=ma)
if ma-mi gt 0.995*!pi then begin
    if centeroffimage then begin
        ind=where(aziout lt 0,ct)
        if ct ne 0 then aziout[ind]+=2*!pi
        mi=min(aziout,max=ma)
        aziout=[mi,ma]
    endif else aziout=[0,2*!pi]
endif else aziout=[mi,ma]

(*list.PDFs)[ind] and= (*list.PDFd)[ind] ge min(d)
if centeroffimage then (*list.PDFs)[ind] and= (*list.PDFd)[ind] le max(d)

; ----Take only the unique rings----
ind=where((*list.PDFs)[i0:i1],ct)+i0
if ct eq 0 then return

ind=ind[uniq((*list.PDFd)[ind],sort((*list.PDFd)[ind]))]
(*list.PDFs)[*]=0
(*list.PDFs)[ind]=1

end;pro DisablePDFrings
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CalibFromPDFHandler,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        if n_elements(progressBar) ne 0 then progressBar -> Destroy
        return
    endif
endif

case widget_info(ev.id,/type) of
8:    begin
        widget_control,ev.top,get_uvalue=tlist
        widget_control,tlist.top,get_uvalue=list
        
        DeletePreviousAutoSet,list,tlist.Call
        widget_control,tlist.top,set_uvalue=list
        
        DisablePDFrings,ev.index,list,aziout=aziout

        if ev.index ne 0 then i0=total((*list.PDFn)[0:ev.index-1]) else i0=0
        i1=total((*list.PDFn)[0:ev.index])-1

        ind=where((*list.PDFs)[i0:i1],ct)+i0
        if ct ne 0 then begin
            width=PromptNumber('1',ev.top,'Ring width in 2-theta degrees:')
            width=float(width)*!pi/360
            
            azi0=PromptNumber(stringr(aziout[0]*180/!pi),ev.top,'Begin azimuth (degrees):')
            azi0=float(azi0)/180.*!pi
            azi1=PromptNumber(stringr(aziout[1]*180/!pi),ev.top,'End azimuth (degrees):')
            azi1=float(azi1)/180.*!pi
        endif else widget_control,tlist.top,set_uvalue=list
        
        
        progressBar = Obj_New("PROGRESSBAR",title="Detecting Debye rings...")
        progressBar -> Start
        for i=0l,ct-1 do begin
            IF progressBar -> CheckCancel() THEN break
    
            d=(*list.PDFd)[ind[i]]
            tt=BraggDToX(d,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b)
            
            ;tmp=azMmtoPixel([azi0<azi1,azi0>azi1],list.scrx,list.scry)
            ;Larc=tt[0]*(tmp[1]-tmp[0])
            ;nSections=ceil(Larc/3)+1
            nSections=0l
            
            tt=tt[1]
            tt=[tt-width,tt+width]
            dspac=BraggTToX(tt,list.lambda,/onlyd)
            ;mask=[3,0,1,dspac[0],dspac[1],0,0]
            mask=[2,0,0,dspac[0],azi0,dspac[1],azi1]
            
            TiI=AutoSetRing(mask,list,nSections)
            if nSections ne 0 then begin
                if list.masknr eq 0 then begin
                    list.mask=ptr_new(mask)
                    list.smask=ptr_new([1b])
                endif else begin
                    (*list.mask)[2,0:list.masknr-1]=1
                    L=max((*list.mask)[1,*])+1
                    mask[1]=L
                    (*list.mask)=[[(*list.mask)],[mask]]
                    (*list.smask)=[(*list.smask),1b]
                endelse
                list.masknr++
                widget_control,tlist.top,set_uvalue=list
            endif
            AddTi,TiI,nSections,tlist.Call,list,tlist.top,d
            
            progressBar -> Update, (i+1.)/ct*100
        endfor
        progressBar -> Destroy
        
    endcase
endcase

WIDGET_CONTROL, ev.top, /DESTROY
end; pro CalibFromPDFHandler
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CalibFromPDF,ev,Call

widget_control,ev.top,get_uvalue=list
if (list.nPDF eq 0) then return

WIDGET_CONTROL, ev.top, sensitive=0
base=widget_base(/column,title='Choose PDF:',uvalue={top:ev.top,Call:Call})
droplist=widget_droplist(base,value=(*list.PDFname))

WIDGET_CONTROL, base, /REALIZE

Xmanager,'CalibFromPDF',base, event_handler='CalibFromPDFHandler',cleanup='CleanUp_auto_set',GROUP_LEADER=ev.top
end;pro CalibFromPDF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro OBJECT_PROCESS_EVENTS, ev, EventPro, i

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN
widget_control,ev.top,get_uvalue=list

list.xs=ev.x
list.ys=ev.y

Widget_Control, ev.id, Event_Pro=EventPro, $
    draw_motion=1
Widget_Control, ev.top, Set_UValue=list
LoadctRev,i,/silent
end;pro OBJECT_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ARC_PROCESS_EVENTS, ev
OBJECT_PROCESS_EVENTS, ev, 'ARC_DRAW', 2
end;pro ARC_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CIRCLE_PROCESS_EVENTS, ev
OBJECT_PROCESS_EVENTS, ev, 'CIRCLE_DRAW', 2
end;pro CIRCLE_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SQUARE_PROCESS_EVENTS, ev
OBJECT_PROCESS_EVENTS, ev, 'SQUARE_DRAW', 2
end;pro SQUARE_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ELLIPSEM_PROCESS_EVENTS, ev
OBJECT_PROCESS_EVENTS, ev, 'ELLIPSEM_DRAW', 3
end;pro ELLIPSEM_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ARC_DRAW ,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]

IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF

; calcultate image position of cursor
x = [list.xs, ev.x]
y = [list.ys, ev.y]
if ev.id eq list.draw then begin
    x=0>StatTR(x,list.sfh)<((*list.tiffs)[1]-1)
    y=0>StatTR(y,list.sfv)<((*list.tiffs)[2]-1)
endif else begin
    x=0>DynTR(x,list.sf,list.hrange[0])<((*list.tiffs)[1]-1)
    y=0>DynTR(y,list.sf,list.vrange[0])<((*list.tiffs)[2]-1)
endelse

; calculate azimuths
xm=x-list.center[0]
ym=y-list.center[1]
azi0=atan(ym[0],xm[0])
azi1=atan(ym[1],xm[1])

; change the calculated azimuths when passing through 180deg
if (xm[1] lt 0) and $
(list.thetadyn lt 0) and (azi1 gt 0) then $
    list.dangle-=2*!dpi
if (xm[1] lt 0) and $
(list.thetadyn gt 0) and (azi1 lt 0) then $
    list.dangle+=2*!dpi
list.thetadyn=azi1
azi1+=list.dangle

; Convert axi from pixel to mm frame
tmp=azPixeltoMm([azi0,azi1],list.scrx,list.scry)
azi0=tmp[0]
azi1=tmp[1]

; Calc two-theta
tt=BraggXYtoX(x,y,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/onlyt)

IF thisEvent EQ 'DOWN' THEN BEGIN
    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
        Widget_Control, list.drawdyn, draw_motion=0, $
           Event_Pro='ARC_PROCESS_EVENTS'
    endif
    Widget_Control, list.draw, draw_motion=0, $
       Event_Pro='ARC_PROCESS_EVENTS'

    ; calc d-spacings
    list.thetadyn=0
    list.dangle=0
    dspac=BraggTtoX(tt,list.lambda,/onlyd)

    ; Set mask entry
    if list.masknr eq 0 then begin ; New
        list.mask=ptr_new([2,0,0,dspac[0],azi0,dspac[1],azi1])
        list.smask=ptr_new([1b])
        list.masknr++
    endif else begin
        ind=where((*list.mask)[2,*] eq 0,count)
        if count eq 0 then begin ; Add
            L=max((*list.mask)[1,*])+1
            (*list.mask)=[[(*list.mask)],[2,L,0,dspac[0],azi0,dspac[1],azi1]]
            (*list.smask)=[(*list.smask),1b]
            list.masknr++
        endif else begin ; Replace last
            L=(*list.mask)[1,ind[0]]
            (*list.mask)[*,ind[0]]=[2,L,0,dspac[0],azi0,dspac[1],azi1]
        endelse
    endelse

    RefreshDisplay, ev.top,list
    return
endif

; whipe out previously drawn arc
WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
endif

; Calc corner coordinates
xya=BraggXtoXY(reverse(tt),[azi0,azi1],list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,1,validphi=valid)
xp=[x[0],xya[0,0],x[1],xya[1,0]] ; list.xs, ev.x(Azi0), ev.x, list.xs(Azi1)
yp=[y[0],xya[0,1],y[1],xya[1,1]]
ArcRecalcCorners,xp,yp,list.center,(*list.tiffs)[1]-1,(*list.tiffs)[2]-1,[1,valid[0],1,valid[1]]

; Calc phi range
phi=azDetRange(azi0,azi1,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

; plot 2 arcs
col=160
for i=0l,1 do begin
    xEllipse=EllipsePlotCoord(tt[i],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,gapshift=gapshift,nvalid=nvalid)
    if nvalid ne 0 then begin
        WSet, list.drawindex
        x2=RTStat(xEllipse[*,0],list.sfh)
        y2=RTStat(xEllipse[*,1],list.sfv)
        
        if gapshift eq 0 then $
            PlotS,transpose([[x2],[y2]]),/device,color=col $
        else begin
            PlotS,transpose([[x2[gapshift:*]],[y2[gapshift:*]]]),/device,color=col
            PlotS,transpose([[x2[0:gapshift-1]],[y2[0:gapshift-1]]]),/device,color=col
        endelse
    
        if widget_info(list.drawdyn,/VALID_ID) then begin
            WSet, list.drawdynindex
            x2=RTDyn(xEllipse[*,0],list.sf,list.hrange[0])
            y2=RTDyn(xEllipse[*,1],list.sf,list.vrange[0])
            PlotS,transpose([[x2],[y2]]),/device,color=col
        endif
    endif
endfor

; plot lines
WSet, list.drawindex
x2=RTStat(xp,list.sfh)
y2=RTStat(yp,list.sfv)
PlotS,transpose([[x2[0:1]],[y2[0:1]]]),/device,color=col
PlotS,transpose([[x2[2:3]],[y2[2:3]]]),/device,color=col
if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    x2=RTDyn(xp,list.sf,list.hrange[0])
    y2=RTDyn(yp,list.sf,list.vrange[0])
    PlotS,transpose([[x2[0:1]],[y2[0:1]]]),/device,color=col
    PlotS,transpose([[x2[2:3]],[y2[2:3]]]),/device,color=col
endif

Widget_Control, ev.top, Set_UValue=list
end;pro ARC_DRAW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CIRCLE_DRAW ,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]

IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF

; calcultate image position of cursor
if list.circle2 then begin
    x = [list.xs, ev.x]
    y = [list.ys, ev.y]
endif else begin
    x = ev.x
    y = ev.y
endelse
if ev.id eq list.draw then begin
    x=0>StatTR(x,list.sfh)<((*list.tiffs)[1]-1)
    y=0>StatTR(y,list.sfv)<((*list.tiffs)[2]-1)
endif else begin
    x=0>DynTR(x,list.sf,list.hrange[0])<((*list.tiffs)[1]-1)
    y=0>DynTR(y,list.sf,list.vrange[0])<((*list.tiffs)[2]-1)
endelse

; Calc two-theta and corner coordinates
tt=BraggXYtoX(x,y,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/onlyt)

IF thisEvent EQ 'DOWN' THEN BEGIN
    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
        Widget_Control, list.drawdyn, draw_motion=0, $
           Event_Pro='CIRCLE_PROCESS_EVENTS'
    endif
    Widget_Control, list.draw, draw_motion=0, $
       Event_Pro='CIRCLE_PROCESS_EVENTS'

    if list.circle2 eq 0 then tt=[0,tt]
    dspac=BraggTtoX(tt,list.lambda,/onlyd)

    ; Set mask entry
    if list.masknr eq 0 then begin ; New
        list.mask=ptr_new([3,0,0,dspac[0],dspac[1],0,0])
        list.smask=ptr_new([1b])
        list.masknr++
    endif else begin
        ind=where((*list.mask)[2,*] eq 0,count)
        if count eq 0 then begin ; Add
            L=max((*list.mask)[1,*])+1
            (*list.mask)=[[(*list.mask)],[3,L,0,dspac[0],dspac[1],0,0]]
            (*list.smask)=[(*list.smask),1b]
            list.masknr++
        endif else begin ; Replace last
            L=(*list.mask)[1,ind[0]]
            (*list.mask)[*,ind[0]]=[3,L,0,dspac[0],dspac[1],0,0]
        endelse
    endelse

    RefreshDisplay, ev.top,list
    return
endif

; whipe out previously drawn arc
WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
endif

; Calc phi range
phi=azDetRange(0,2*!pi,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

; plot 2 circles
col=160
for i=0l,list.circle2 do begin
    xEllipse=EllipsePlotCoord(tt[i],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
    if nvalid ne 0 then begin
        WSet, list.drawindex
        x2=RTStat(xEllipse[*,0],list.sfh)
        y2=RTStat(xEllipse[*,1],list.sfv)
        PlotS,transpose([[x2],[y2]]),/device,color=col
    
        if widget_info(list.drawdyn,/VALID_ID) then begin
            WSet, list.drawdynindex
            x2=RTDyn(xEllipse[*,0],list.sf,list.hrange[0])
            y2=RTDyn(xEllipse[*,1],list.sf,list.vrange[0])
            PlotS,transpose([[x2],[y2]]),/device,color=col
        endif
    endif
endfor

Widget_Control, ev.top, Set_UValue=list
end;pro CIRCLE_DRAW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SQUARE_DRAW ,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]

IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF

; calcultate image position
x=[list.xs,ev.x]
y=[list.ys,ev.y]
if ev.id eq list.draw then begin
    x=0>StatTR(x,list.sfh)<((*list.tiffs)[1]-1)
    y=0>StatTR(y,list.sfv)<((*list.tiffs)[2]-1)
endif else begin
    x=0>DynTR(x,list.sf,list.hrange[0])<((*list.tiffs)[1]-1)
    y=0>DynTR(y,list.sf,list.vrange[0])<((*list.tiffs)[2]-1)
endelse

IF thisEvent EQ 'DOWN' THEN BEGIN
    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
        Widget_Control, list.drawdyn, draw_motion=0, $
           Event_Pro='SQUARE_PROCESS_EVENTS'
    endif
    Widget_Control, list.draw, draw_motion=0, $
       Event_Pro='SQUARE_PROCESS_EVENTS'

    ; Set mask
    x=round(x)
    y=round(y)
    x=x[sort(x)]
    y=y[sort(y)]
    if list.masknr eq 0 then begin ; New
        list.mask=ptr_new([1,0,0,x,y])
        list.smask=ptr_new([1b])
        list.masknr++
    endif else begin
        ind=where((*list.mask)[2,*] eq 0,count)
        if count eq 0 then begin ; Add
            L=max((*list.mask)[1,*])+1
            (*list.mask)=[[(*list.mask)],[1,L,0,x,y]]
            (*list.smask)=[(*list.smask),1b]
            list.masknr++
        endif else begin ; Replace last
            L=(*list.mask)[1,ind[0]]
            (*list.mask)[*,ind[0]]=[1,L,0,x,y]
        endelse
    endelse

    RefreshDisplay,ev.top,list
    LoadctRev,0,/silent,rev=list.revcol
    return
endif

; whipe out previously drawn square
WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
endif

WSet, list.drawindex
xx=RTStat(x,list.sfh)
yy=RTStat(y,list.sfv)
col=160
PlotS, [xx[0], xx[0]], yy,/device,color=col
PlotS, [xx[1], xx[1]], yy,/device,color=col
PlotS, xx,[yy[0], yy[0]],/device,color=col
PlotS, xx,[yy[1], yy[1]],/device,color=col
if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    xx=RTDyn(x,list.sf,list.hrange[0])
    yy=RTDyn(y,list.sf,list.vrange[0])
    PlotS, [xx[0], xx[0]], yy,/device,color=col
    PlotS, [xx[1], xx[1]], yy,/device,color=col
    PlotS, xx,[yy[0], yy[0]],/device,color=col
    PlotS, xx,[yy[1], yy[1]],/device,color=col
endif

Widget_Control, ev.top, Set_UValue=list
end;pro SQUARE_DRAW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TTWindowUpdate,tt,list,top
ID=widget_info(top,FIND_BY_UNAME='MEllipsed')
if ~widget_info(ID,/valid) then return
res=BraggTtoX(tt,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/pixm)
widget_control,ID,set_value=stringr(res[1])
ID=widget_info(top,FIND_BY_UNAME='MEllipseq')
widget_control,ID,set_value=stringr(res[2])
ID=widget_info(top,FIND_BY_UNAME='MEllipsefr')
widget_control,ID,set_value=stringr(res[0])
ID=widget_info(top,FIND_BY_UNAME='MEllipset')
widget_control,ID,set_value=stringr(tt*180./!dpi)
end;pro TTWindowUpdate
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ELLIPSEM_DRAW ,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]

IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF

; calcultate image position of cursor
if ev.id eq list.draw then begin
    x=0>StatTR(ev.x,list.sfh)<((*list.tiffs)[1]-1)
    y=0>StatTR(ev.y,list.sfv)<((*list.tiffs)[2]-1)
endif else begin
    x=0>DynTR(ev.x,list.sf,list.hrange[0])<((*list.tiffs)[1]-1)
    y=0>DynTR(ev.y,list.sf,list.vrange[0])<((*list.tiffs)[2]-1)
endelse

ret=BraggXYtoX(x,y,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/pixm)
tt=ret[1]
azi=atan(y-list.center[1],x-list.center[0])
Rnonort=BraggTtoX(tt,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,/pixm,/onlyp,/nonortp,phi=azPixeltoMm(azi,list.scrx,list.scry))
azi*=180/!pi
if azi lt 0 then azi+=360

ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsed')
widget_control,ID,set_value=stringr(ret[2])
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipseq')
widget_control,ID,set_value=stringr(ret[3])
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsep')
widget_control,ID,set_value=stringr(x)+','+stringr(y)
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsefr')
widget_control,ID,set_value=stringr(ret[0])
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipser')
widget_control,ID,set_value=stringr(Rnonort)
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipset')
widget_control,ID,set_value=stringr(tt*180./!dpi)
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsea')
widget_control,ID,set_value=stringr(azi)
ID=widget_info(ev.top,FIND_BY_UNAME='warning')
tmp=ConicTypeTT(list.a,tt,str=str)
widget_control,ID,set_value=str

ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsei')
int=(*list.tiff)[round(x),round(y)]
if (size(int,/type)) eq 1 then int=fix(int); fix for if tiff is in bytes
widget_control,ID,set_value=stringr(int)

if ptr_valid(list.CHIchild) and list.CrossUpdate then begin
    val=tt*180./!pi
    for i=0l,n_elements(*list.CHIchild)-1 do $
        SetAnySelect,val,1,(*list.CHIchild)[i]
endif

IF thisEvent EQ 'DOWN' THEN BEGIN
    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    if widget_info(list.drawdyn,/VALID_ID) then begin
        WSet, list.drawdynindex
        Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
        Widget_Control, list.drawdyn, draw_motion=0, $
           Event_Pro='ELLIPSEM_PROCESS_EVENTS'
    endif
    Widget_Control, list.draw, draw_motion=0, $
       Event_Pro='ELLIPSEM_PROCESS_EVENTS'

    list.ttmellipse=tt
    list.smellipse=list.smellipse or 1
    ID=widget_info(ev.top,FIND_BY_UNAME='Show Debye Marker')
    widget_control,ID,set_value='# Show Debye Marker'

    RefreshDisplay, ev.top,list
    return
endif

; whipe out previously drawn ellipse
WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
endif

; Calc coordinates
phi=azDetRange(0,2*!pi,list.phipixphi0,list.scrx,list.scry,list.a,list.b)
xEllipse=EllipsePlotCoord(tt,phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
if nvalid eq 0 then return

; Draw
WSet, list.drawindex
col=160
x2=RTStat(xEllipse[*,0],list.sfh)
y2=RTStat(xEllipse[*,1],list.sfv)
PlotS,transpose([[x2],[y2]]),/device,color=col

if (list.smellipse and 2) eq 2 then begin
    xAxis=SemiAxisCoord(tt,list.a,list.b,list.dist,list.scrx,list.scry,list.center)
    xx2=RTStat(xAxis[*,0],list.sfh)
    yy2=RTStat(xAxis[*,1],list.sfv)
    Plots,transpose([[xx2[0:1]],[yy2[0:1]]]),/device,color=col,linestyle=2
    Plots,transpose([[xx2[2:3]],[yy2[2:3]]]),/device,color=col
endif

if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    x2=RTDyn(xEllipse[*,0],list.sf,list.hrange[0])
    y2=RTDyn(xEllipse[*,1],list.sf,list.vrange[0])
    PlotS,transpose([[x2],[y2]]),/device,color=col

    if (list.smellipse and 2) eq 2 then begin
        xx2=RTDyn(xAxis[*,0],list.sf,list.hrange[0])
        yy2=RTDyn(xAxis[*,1],list.sf,list.vrange[0])
        Plots,transpose([[xx2[0:1]],[yy2[0:1]]]),/device,color=col,linestyle=2
        Plots,transpose([[xx2[2:3]],[yy2[2:3]]]),/device,color=col
    endif
endif

Widget_Control, ev.top, Set_UValue=list
end;pro ELLIPSEM_DRAW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AZIOFF_DRAW_MOVE,ev,list,x,y,tt

LoadctRev,39,/silent

; Calcultate image position of cursor
if ev.id eq list.draw then begin
    x=0>StatTR(ev.x,list.sfh)<((*list.tiffs)[1]-1)
    y=0>StatTR(ev.y,list.sfv)<((*list.tiffs)[2]-1)
endif else begin
    x=0>DynTR(ev.x,list.sf,list.hrange[0])<((*list.tiffs)[1]-1)
    y=0>DynTR(ev.y,list.sf,list.vrange[0])<((*list.tiffs)[2]-1)
endelse

; Calculate 2theta (radians)
ret=BraggXYtoX(x,y,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/pixm)
tt=ret[1]

; whipe out previously drawn ellipse
WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
endif

; Calc coordinates
phi=azDetRange(0,2*!pi,list.phipixphi0,list.scrx,list.scry,list.a,list.b)
xEllipse=EllipsePlotCoord(tt,phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
if nvalid eq 0 then return

; Draw
WSet, list.drawindex
col=250
col2=150
col3=200
x2=RTStat(xEllipse[*,0],list.sfh)
y2=RTStat(xEllipse[*,1],list.sfv)
PlotS,transpose([[x2],[y2]]),/device,color=col

xx2=RTStat([0,(*list.tiffs)[1]-1,list.center[0],list.center[0]],list.sfh)
yy2=RTStat([list.center[1],list.center[1],0,(*list.tiffs)[2]-1],list.sfv)
arrow,xx2[0],yy2[0],xx2[1],yy2[1],color=col2
arrow,xx2[2],yy2[2],xx2[3],yy2[3],color=col3

xx2=RTStat([list.center[0],x],list.sfh)
yy2=RTStat([list.center[1],y],list.sfv)
Plots,transpose([[xx2[0:1]],[yy2[0:1]]]),/device,color=col

if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    x2=RTDyn(xEllipse[*,0],list.sf,list.hrange[0])
    y2=RTDyn(xEllipse[*,1],list.sf,list.vrange[0])
    PlotS,transpose([[x2],[y2]]),/device,color=col

    xx2=RTDyn([0,(*list.tiffs)[1]-1,list.center[0],list.center[0]],list.sf,list.hrange[0])
    yy2=RTDyn([list.center[1],list.center[1],0,(*list.tiffs)[2]-1],list.sf,list.vrange[0])
    arrow,xx2[0],yy2[0],xx2[1],yy2[1],color=col2
    arrow,xx2[2],yy2[2],xx2[3],yy2[3],color=col3
    
    xx2=RTDyn([list.center[0],x],list.sf,list.hrange[0])
    yy2=RTDyn([list.center[1],y],list.sf,list.vrange[0])
    Plots,transpose([[xx2[0:1]],[yy2[0:1]]]),/device,color=col
endif

end;pro AZIOFF_DRAW_MOVE
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AZIOFF_DRAW_MOVE2,tt,list
; Convert from cone frame to detector frame (mm)
phid=azConeShiftToTilt([!dpi,0,3*!dpi/2,!dpi/2],list.phipixphi0,list.scrx,list.scry,list.a,list.b)
xy=BraggXtoXY(replicate(tt,4),phid,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,1)

; Draw
WSet, list.drawindex
col2=150
col3=200

xx2=RTStat(xy[*,0],list.sfh)
yy2=RTStat(xy[*,1],list.sfv)
arrow,xx2[0],yy2[0],xx2[1],yy2[1],color=col2
arrow,xx2[2],yy2[2],xx2[3],yy2[3],color=col3

if widget_info(list.drawdyn,/VALID_ID) then begin
    WSet, list.drawdynindex
    
    xx2=RTDyn(xy[*,0],list.sf,list.hrange[0])
    yy2=RTDyn(xy[*,1],list.sf,list.vrange[0])
    arrow,xx2[0],yy2[0],xx2[1],yy2[1],color=col2
    arrow,xx2[2],yy2[2],xx2[3],yy2[3],color=col3
endif
end;pro AZIOFF_DRAW_MOVE2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AZIOFF_DRAW ,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]

IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF

; Get pixel positions and draw ellipse
AZIOFF_DRAW_MOVE,ev,list,x,y,tt

; Azimuth in pixel or mm
phipix=atan(y-list.center[1],x-list.center[0])
phid=azPixeltoMm(phipix,list.scrx,list.scry)

; Calculate phipixphi0
IF thisEvent EQ 'DOWN' THEN BEGIN
    ID=widget_info(ev.top,FIND_BY_UNAME='phi0')
    widget_control,ID,get_value=phi0
    list.phipixphi0=Getphipixphi0(phipix,phi0*!dpi/180,list.a,list.b,list.scrx,list.scry)
    Widget_Control, ev.top, Set_UValue=list
endif

; Plot axis of the orthogonal frame
AZIOFF_DRAW_MOVE2,tt,list

; Get the azimuth in the orthogonal frame (a.k.a. cone frame, sample frame, flat frame)
phi=azTiltToConeShift(phid,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

; Output
phipix*=180/!pi
if phipix lt 0 then phipix+=360
phipix mod= 360

phid*=180/!pi
if phid lt 0 then phid+=360
phid mod= 360

phi*=180/!pi
if phi lt 0 then phi+=360
phi mod= 360

ID=widget_info(ev.top,FIND_BY_UNAME='phipix')
widget_control,ID,set_value=stringr(phipix)
ID=widget_info(ev.top,FIND_BY_UNAME='phid')
widget_control,ID,set_value=stringr(phid)
ID=widget_info(ev.top,FIND_BY_UNAME='phi')
widget_control,ID,set_value=stringr(phi)
ID=widget_info(ev.top,FIND_BY_UNAME='phipixphi0')
widget_control,ID,set_value=stringr(list.phipixphi0*180/!pi)

end;pro AZIOFF_DRAW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro usemouse_events, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN
widget_control,ev.top,get_uvalue=list

if ev.id eq list.drawdyn then begin
list.center=DynTR([ev.x,ev.y],list.sf,[list.hrange[0],list.vrange[0]])
endif else begin
list.center=StatTR([ev.x,ev.y],[list.sfh,list.sfv])
endelse

str=stringr(list.center)
ID=widget_info(ev.top,FIND_BY_UNAME='center')
widget_control,ID,set_value=str[0]+','+str[1]
RefreshDisplay, ev.top,list
end;pro usemouse_events
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro circle_events, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN
widget_control,ev.top,get_uvalue=list

if ev.id eq list.drawdyn then begin
xy=DynTR([ev.x,ev.y],list.sf,[list.hrange[0],list.vrange[0]])
endif else begin
xy=StatTR([ev.x,ev.y],[list.sfh,list.sfv])
endelse

if list.bP eq 1 then begin
list.P1=xy
str=stringr(list.P1)
ID=widget_info(ev.top,FIND_BY_UNAME='P1')
widget_control,ID,set_value=str[0]+','+str[1]
endif
if list.bP eq 2 then begin
list.P2=xy
str=stringr(list.P2)
ID=widget_info(ev.top,FIND_BY_UNAME='P2')
widget_control,ID,set_value=str[0]+','+str[1]
endif
if list.bP eq 3 then begin
list.P3=xy
str=stringr(list.P3)
ID=widget_info(ev.top,FIND_BY_UNAME='P3')
widget_control,ID,set_value=str[0]+','+str[1]
endif

RefreshDisplay, ev.top,list
end;pro circle_events
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function calc_middle,P1,P2,P3
;line P1P2: vgl: y=(x-x1).(y2-y1)/(x2-x1)+y1=(x-x1).r12+y1
;      middle of the two points: [(x1+x2)/2,(y1+y2)/2]
;      rico bisector: -1/r12
;      vgl bisector: y=-(x-xm)/r12+ym with (xm,ym) middle of the two points
;do this for the other 2 combinations
;eliminate bisectors with rico 0 or infinite
;choose two biscetor equations and solve to x and y to calculate the center

;make floating
P1=float(P1)
P2=float(P2)
P3=float(P3)

;take two rico's a and b how are not 0 of infinite
r12=(P2[1]-P1[1])/(P2[0]-P1[0])
r23=(P3[1]-P2[1])/(P3[0]-P2[0])
r13=(P3[1]-P1[1])/(P3[0]-P1[0])
r=[r12,r23,r13]
ind1=finite(r)
ind2=r ne 0
ind3=ind1*ind2
ind=where(ind3 eq 1,count)
if count lt 2 then return,[0,0]
a=r[ind[0]]
b=r[ind[1]]

;take the two corresponding middles
m12=[(P1[0]+P2[0])/2,(P1[1]+P2[1])/2]
m23=[(P3[0]+P2[0])/2,(P3[1]+P2[1])/2]
m13=[(P1[0]+P3[0])/2,(P1[1]+P3[1])/2]
m={m12:m12,m23:m23,m13:m13}
ma=m.(ind[0])
mb=m.(ind[1])

;result when solving two bisector equations
c1=ma[0]/a+ma[1]
c2=mb[0]/b+mb[1]
x=(c2-c1)/(1/b-1/a)
y=-x/a+c1
return,[x,y]
end;function calc_middle
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro tplot, list
tvlct,RoldCT,GoldCT,BoldCT,/get
LoadctRev,0,/silent,rev=list.revcol
xorcol = byte(-list.revcol)

;----Display image----
image=congrid(*list.tiff, list.stath, list.statv)
tv, WeakScl(image,list.sclmax,list.displaypower)

;----Display point to determine middle----
if list.cross eq 1 then begin
point=list.P1
if point[0] ne -1 and point[1] ne -1 then begin
x=RTstat(point[0],list.sfh)
y=RTstat(point[1],list.sfv)
Plots,[x-list.CW1,x+list.CW1],y,/device
Plots,x,[y-list.CW1,y+list.CW1],/device
endif
point=list.P2
if point[0] ne -1 and point[1] ne -1 then begin
x=RTstat(point[0],list.sfh)
y=RTstat(point[1],list.sfv)
Plots,[x-list.CW1,x+list.CW1],y,/device
Plots,x,[y-list.CW1,y+list.CW1],/device
endif
point=list.P3
if point[0] ne -1 and point[1] ne -1 then begin
x=RTstat(point[0],list.sfh)
y=RTstat(point[1],list.sfv)
Plots,[x-list.CW1,x+list.CW1],y,/device
Plots,x,[y-list.CW1,y+list.CW1],/device
endif
point=list.center
if point[0] ne -1 and point[1] ne -1 then begin
x=RTstat(point[0],list.sfh)
y=RTstat(point[1],list.sfv)
Plots,[x-list.CW1,x+list.CW1],y,/device
Plots,x,[y-list.CW1,y+list.CW1],/device
endif
endif

;----Display bad pixels----
if ptr_valid(list.tiffvalid) and (list.MarkBadPixels eq 1) then begin
    bad=where((*list.tiffvalid) eq 0,ct)
    if ct ne 0 then begin
;        DEVICE, GET_GRAPHICS = oldg, SET_GRAPHICS=6          ;Xor drawing mode
        LoadctRev,3,/silent
        x=RTStatExt(bad mod (*list.tiffs)[1],list.sfh,/x)
        y=RTStatExt(bad / (*list.tiffs)[1],list.sfv)
        plots,x,y,/device,psym=3,color=150
;        DEVICE, SET_GRAPHICS=oldg          ;Copy mode
    endif
endif

;----Display mask----
if (list.masknr ge 1) then begin
LoadctRev,2,/silent
center=list.center
color=10*5

for i=0l,list.masknr-1 do begin
if (*list.smask)[i] eq 1 then begin
    if ((*list.mask)[1,i] eq list.maskP) then color=10*2 else color=10*16
    arr=reform((*list.mask)[*,i])
    str= stringr(fix(i))+';'+stringr(fix(arr[1]))

    case arr[0] of
    1: begin
       x=RTstat(arr[3:4],list.sfh)
       y=RTstat(arr[5:6],list.sfv)
       PlotS, [x[0], x[0]], y,/device,color=color
       PlotS, [x[1], x[1]], y,/device,color=color
       PlotS, x,[y[0], y[0]],/device,color=color
       PlotS, x,[y[1], y[1]],/device,color=color
       ;plot numbers
       if arr[2] eq 1 then begin
       xlabel=x[1]
       ylabel=y[0]
       xyouts,xlabel,ylabel,str,/device,charsize=list.letter,color=color
       endif
       endcase
    2: begin
        tt=BraggDtoX(arr[[3,5]],list.lambda,/onlyt)
        azi0=arr[4]
        azi1=arr[6]

        xya=BraggXtoXY([tt,tt],[azi0,azi0,azi1,azi1],list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,1,validphi=valid)
        xp=xya[*,0]
        yp=xya[*,1]
        ArcRecalcCorners,xp,yp,list.center,(*list.tiffs)[1]-1,(*list.tiffs)[2]-1,valid
        
        ; Calc phi range
        phi=azDetRange(azi0,azi1,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

        ; plot 2 arcs
        for j=0l,1 do begin
            xEllipse=EllipsePlotCoord(tt[j],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,gapshift=gapshift,nvalid=nvalid)
            if nvalid ne 0 then begin
                x2=RTStat(xEllipse[*,0],list.sfh)
                y2=RTStat(xEllipse[*,1],list.sfv)
                
                if gapshift eq 0 then $
                    PlotS,transpose([[x2],[y2]]),/device,color=color $
                else begin
                    PlotS,transpose([[x2[gapshift:*]],[y2[gapshift:*]]]),/device,color=color
                    PlotS,transpose([[x2[0:gapshift-1]],[y2[0:gapshift-1]]]),/device,color=color
                endelse
            endif
        endfor

        ; plot lines
        x2=RTStat(xp,list.sfh)
        y2=RTStat(yp,list.sfv)
        PlotS,transpose([[x2[0:1]],[y2[0:1]]]),/device,color=color
        PlotS,transpose([[x2[2:3]],[y2[2:3]]]),/device,color=color

       ;plot numbers
       if arr[2] eq 1 then xyouts,x2[0],y2[0],str,/device,charsize=list.letter,color=color
       endcase
    3: begin
           tt=BraggDtoX(arr[[3,4]],list.lambda,/onlyt)
           if tt[0] eq 0 then tt=tt[1]

        ; Calc phi range
        phi=azDetRange(0,2*!pi,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

        ; plot 2 lines
        for j=0l,n_elements(tt)-1 do begin
            xEllipse=EllipsePlotCoord(tt[j],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
            if nvalid ne 0 then begin
                x2=RTStat(xEllipse[*,0],list.sfh)
                y2=RTStat(xEllipse[*,1],list.sfv)
                PlotS,transpose([[x2],[y2]]),/device,color=color
            endif
        endfor

       ;plot numbers
       if arr[2] eq 1 then begin
       xya=BraggXtoXY(tt,tt*0,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,1)
       x2=RTstat(xya[0,0],list.sfh)
       y2=RTstat(xya[0,1],list.sfv)
       xyouts,x2,y2,str,/device,charsize=list.letter,color=color
       endif
       endcase
    endcase
endif
endfor
endif

;----Display mark ellipse----
if list.smellipse ne 0 then begin
    LoadctRev,3,/silent
    tt=list.ttmellipse
    phi=azDetRange(0,2*!pi,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

    if (list.smellipse and 1) eq 1 then begin
        xEllipse=EllipsePlotCoord(tt,phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
        if nvalid ne 0 then begin
            x2=RTStat(xEllipse[*,0],list.sfh)
            y2=RTStat(xEllipse[*,1],list.sfv)
            PlotS,transpose([[x2],[y2]]),/device,color=16*10
        endif
    endif

    if (list.smellipse and 2) eq 2 then begin
        xAxis=SemiAxisCoord(tt,list.a,list.b,list.dist,list.scrx,list.scry,list.center)
        xx2=RTStat(xAxis[*,0],list.sfh)
        yy2=RTStat(xAxis[*,1],list.sfv)
        Plots,transpose([[xx2[0:1]],[yy2[0:1]]]),/device,color=16*10,linestyle=2
        Plots,transpose([[xx2[2:3]],[yy2[2:3]]]),/device,color=16*10
    endif
endif

;----Display Tie points----
if total(list.ShowTie) ne 0 and (*list.nrTie)[0] ne 0 then begin
LoadctRev,1,/silent
i=0l
; Input Tie points
if list.ShowTie[0] then begin
if list.nrTieRing ne 0 then $
    for i=0l,total((*list.nrTie)[0:list.nrTieRing-1])-1 do begin
       x=RTstat((*list.TieI)[0,i],list.sfh)
       y=RTstat((*list.TieI)[1,i],list.sfv)
       Plots,[x-list.CW1,x+list.CW1],y,/device,color=7*16
       Plots,x,[y-list.CW1,y+list.CW1],/device,color=7*16
    endfor
for j=0l+i,total((*list.nrTie)[list.nrTieRing])+i-1 do begin
    x=RTstat((*list.TieI)[0,j],list.sfh)
    y=RTstat((*list.TieI)[1,j],list.sfv)
    Plots,[x-list.CW1,x+list.CW1],y,/device,color=10*16
    Plots,x,[y-list.CW1,y+list.CW1],/device,color=10*16
endfor
endif
if list.TieGrid then begin
    LoadctRev,3,/silent
    ; Output Tie points
    if list.ShowTie[1] then begin
        if list.ShowTie[3] then begin
            ind=ConnectGrid(list.TieO,$
                list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix,$
                count=ct,ngaps=ngaps,gapcoord=gapcoord)
            x=RTstat((*list.TieO)[0,*],list.sfh)
            y=RTstat((*list.TieO)[1,*],list.sfv)
            for i=0l,ct-1 do Plots,x[[ind[i,0],ind[i,1]]],$
                  y[[ind[i,0],ind[i,1]]],/device,color=10*16
            if ngaps ne 0 then begin
            gapcoord[*,0:1]=RTstat(gapcoord[*,0:1],list.sfh)
            gapcoord[*,2:3]=RTstat(gapcoord[*,2:3],list.sfv)
            for i=0l,ngaps-1 do Plots,gapcoord[i,0:1],$
                  gapcoord[i,2:3],/device,color=10*16,linestyle=1
            endif
        endif else begin
        for j=0l,total((*list.nrTie)[0])-1 do begin
            x=RTstat((*list.TieO)[0,j],list.sfh)
            y=RTstat((*list.TieO)[1,j],list.sfv)
            Plots,[x-list.CW1,x+list.CW1],y,/device,color=10*16
            Plots,x,[y-list.CW1,y+list.CW1],/device,color=10*16
        endfor
        endelse
    endif
    ; Connect I-O
    if list.ShowTie[2] then $
     for j=0l,total((*list.nrTie)[0])-1 do begin
        point=[[(*list.TieI)[*,j]], [(*list.TieO)[*,j]]]
        x=RTstat(point[0,*],list.sfh)
        y=RTstat(point[1,*],list.sfv)
        Plots,x,y,/device,color=10*16
    endfor
endif
endif

;----Display Tilt points----
if list.ShowTilt eq 1 and (*list.nrTilt)[0] ne 0 then begin
LoadctRev,3,/silent
i=0l
if list.nrTiltRing ne 0 then $
    for i=0l,total((*list.nrTilt)[0:list.nrTiltRing-1])-1 do begin
       x=RTstat((*list.TiltI)[0,i],list.sfh)
       y=RTstat((*list.TiltI)[1,i],list.sfv)
       Plots,[x-list.CW1,x+list.CW1],y,/device,color=7*16
       Plots,x,[y-list.CW1,y+list.CW1],/device,color=7*16
    endfor
for j=0l+i,total((*list.nrTilt)[list.nrTiltRing])-1+i do begin
    x=RTstat((*list.TiltI)[0,j],list.sfh)
    y=RTstat((*list.TiltI)[1,j],list.sfv)
    Plots,[x-list.CW1,x+list.CW1],y,/device,color=10*16
    Plots,x,[y-list.CW1,y+list.CW1],/device,color=10*16
endfor
endif

;----Display PDF data----
if (list.PDFnTotal ne 0) then begin
    LoadctRev,ntableadd(1),/silent
    off=0
    for i=0l,list.nPDF-1 do begin
        ind=where((*list.PDFs)[off:off+(*list.PDFn)[i]-1] ne 0,ct)+off
        if ct ne 0 then begin
            tt=BraggDtoX((*list.PDFd)[ind],list.lambda,/onlyt)
            phi=azDetRange(0,2*!pi,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

            ;plot circles
            if ((*list.PDFni)[i] ne 0) then PDFi=ScaleArr((*list.PDFi)[ind],1,3) $
            else PDFi=replicate(1,ct)

            ; Labels
            hkl=fix((*list.PDFhkl)[*,ind])
            mihkl=min(hkl)
            n=strlen(stringr(max(hkl)))>(strlen(stringr(mihkl)))
            label=string(hkl,format='(3(I'+stringr(n)+'))')
            code=(*list.PDFs)[ind]

            col=i+1
            PDFangle=45*((i-(i/8)*8)+1)/180.*!dpi
            for j=0l,ct-1 do begin
               if (code[j] and 1) eq 1 then begin
                       xEllipse=EllipsePlotCoord(tt[j],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
                       if nvalid ne 0 then begin
                        x2=RTStat(xEllipse[*,0],list.sfh)
                        y2=RTStat(xEllipse[*,1],list.sfv)
                           PlotS,transpose([[x2],[y2]]),/device,color=col,thick=PDFi[j]
                       endif
               endif
               if (code[j] and 2) eq 2 then begin
                       xEllipse=EllipsePlotCoord(tt[j],PDFangle,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
                       if nvalid ne 0 then begin
                           x2=RTstat([xEllipse[*,0],xEllipse[*,0]+0.5*list.CW1],list.sfh)
                           y2=RTstat([xEllipse[*,1],xEllipse[*,1]+1.5*list.CW1],list.sfv)
                           Plots,x2,y2,/device,color=col
                           xyouts,x2[1],y2[1],label[j],color=col,/device,charsize=list.letter
                       endif
               endif
            endfor
        endif
        off+=(*list.PDFn)[i]
    endfor
endif

;----Display dynamic window middle marker----
if list.bPdyn ne 0 then begin
LoadctRev,4,/silent
point=list.Pdyn
if point[0] ne -1 and point[1] ne -1 then begin
    x=RTstat(point[0],list.sfh)
    y=RTstat(point[1],list.sfv)
    hrange=RTstat(list.hrange,list.sfh)
    vrange=RTstat(list.vrange,list.sfv)
    col=15*16
    Plots,hrange,y,/device,color=col
    Plots,x,vrange,/device,color=col
    Plots,hrange,[vrange[0],vrange[0]],/device,color=col
    Plots,hrange,[vrange[1],vrange[1]],/device,color=col
    Plots,[hrange[0],hrange[0]],vrange,/device,color=col
    Plots,[hrange[1],hrange[1]],vrange,/device,color=col
endif
endif

; ----Display Axis----
if list.RadAxShow then begin
    LoadctRev,3,/silent
    col=16*10
    xx=list.RadAxMin+(list.RadAxMax-list.RadAxMin)/list.RadAxBins*lindgen(list.RadAxBins+1)
    phi=make_array(list.RadAxBins+1,value=list.RadAxAzimuth*!pi/180)
    out=BraggXtoXY(xx,phi,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,list.RadAxType,/angledeg,/pixm,/nonort)
    x=RTstat(out[*,0],list.sfh)
    y=RTstat(out[*,1],list.sfv)
    Plots,x[[0,list.RadAxBins]],y[[0,list.RadAxBins]],/device,color=col
    
    angle=(90-list.RadAxAzimuth)*!pi/180
    addx=list.CW2*cos(angle)
    addy=list.CW2*sin(angle)
    for i=0,list.RadAxBins do begin
        Plots,[x[i],x[i]-addx],[y[i],y[i]+addy],/device,color=col
        xyouts,x[i]+2.5*addx,y[i]-2.5*addy,string(xx[i],format=list.RadAxFormat),orientation=list.RadAxAzimuth,/device,color=col
    endfor
endif

tvlct,RoldCT,GoldCT,BoldCT
end; pro tplot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro dplot, list, dummy=dummy
tvlct,RoldCT,GoldCT,BoldCT,/get
LoadctRev,0,/silent,rev=list.revcol
xorcol = byte(-list.revcol)

;----Display image----
image=congrid((*list.tiff)[list.hrange[0]:list.hrange[1],list.vrange[0]:list.vrange[1]], list.dynh, list.dynv)
; tv and tvscl puts value with matrixindex [x,y]
;             at point [x,y] of display (lower-left = [0,0])
tv, WeakScl(image,list.sclmax,list.displaypower)

;----Display point to determine middle----
if list.cross eq 1 then begin
point=list.P1
if point[0] ne -1 and point[1] ne -1 then begin
x=RTDyn(point[0],list.sf,list.hrange[0])
y=RTDyn(point[1],list.sf,list.vrange[0])
Plots,[x-list.CW2,x+list.CW2],y,/device,thick=list.dplotthick
Plots,x,[y-list.CW2,y+list.CW2],/device,thick=list.dplotthick
endif
point=list.P2
if point[0] ne -1 and point[1] ne -1 then begin
x=RTDyn(point[0],list.sf,list.hrange[0])
y=RTDyn(point[1],list.sf,list.vrange[0])
Plots,[x-list.CW2,x+list.CW2],y,/device,thick=list.dplotthick
Plots,x,[y-list.CW2,y+list.CW2],/device,thick=list.dplotthick
endif
point=list.P3
if point[0] ne -1 and point[1] ne -1 then begin
x=RTDyn(point[0],list.sf,list.hrange[0])
y=RTDyn(point[1],list.sf,list.vrange[0])
Plots,[x-list.CW2,x+list.CW2],y,/device,thick=list.dplotthick
Plots,x,[y-list.CW2,y+list.CW2],/device,thick=list.dplotthick
endif
point=list.center
if point[0] ne -1 and point[1] ne -1 then begin
x=RTDyn(point[0],list.sf,list.hrange[0])
y=RTDyn(point[1],list.sf,list.vrange[0])
Plots,[x-list.CW2,x+list.CW2],y,/device,thick=list.dplotthick
Plots,x,[y-list.CW2,y+list.CW2],/device,thick=list.dplotthick
endif
endif

;----Display bad pixels----
if ptr_valid(list.tiffvalid) and (list.MarkBadPixels eq 1) then begin
    bad=where((*list.tiffvalid) eq 0,ct)
    if ct ne 0 then begin
;        DEVICE, GET_GRAPHICS = oldg, SET_GRAPHICS=6          ;Xor drawing mode
        LoadctRev,3,/silent
        x=bad mod (*list.tiffs)[1]
        y=bad / (*list.tiffs)[1]
        ind=where(x ge list.hrange[0] and x le list.hrange[1] and y ge list.vrange[0] and y le list.vrange[1],ct)
        if ct ne 0 then begin
            x=RTDynExt(x[ind],list.sf,list.hrange[0],/x)
            y=RTDynExt(y[ind],list.sf,list.vrange[0])
            plots,x,y,/device,psym=3,color=150
        endif
;        DEVICE, SET_GRAPHICS=oldg          ;Copy mode
    endif
endif

;----Display mask----
if (list.masknr ge 1) then begin
LoadctRev,2,/silent
center=list.center

for i=0l,list.masknr-1 do begin
if (*list.smask)[i] eq 1 then begin
    if ((*list.mask)[1,i] eq list.maskP) then color=10*2 else color=10*16
    arr=reform((*list.mask)[*,i])
    str=stringr(fix(i))+';'+stringr(fix(arr[1]))

    case arr[0] of
    1: begin
       x=RTDyn(arr[3:4],list.sf,list.hrange[0])
       y=RTDyn(arr[5:6],list.sf,list.vrange[0])
       PlotS, [x[0], x[0]], y,/device,color=color,thick=list.dplotthick
       PlotS, [x[1], x[1]], y,/device,color=color,thick=list.dplotthick
       PlotS, x,[y[0], y[0]],/device,color=color,thick=list.dplotthick
       PlotS, x,[y[1], y[1]],/device,color=color,thick=list.dplotthick
       ;plot numbers
       if arr[2] eq 1 then begin
       xlabel=x[1]
       ylabel=y[0]
       xyouts,xlabel,ylabel,str,/device,charsize=list.letter,color=color,charthick=list.dplotthick
       endif
       endcase
    2: begin
           tt=BraggDtoX(arr[[3,5]],list.lambda,/onlyt)
        azi0=arr[4]
        azi1=arr[6]

        xya=BraggXtoXY([tt,tt],[azi0,azi0,azi1,azi1],list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,1,validphi=valid)
        xp=xya[*,0]
        yp=xya[*,1]
        ArcRecalcCorners,xp,yp,list.center,(*list.tiffs)[1]-1,(*list.tiffs)[2]-1,valid
        
        ; Calc phi range
        phi=azDetRange(azi0,azi1,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

        ; plot 2 arcs
        for j=0l,1 do begin
            xEllipse=EllipsePlotCoord(tt[j],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,gapshift=gapshift,nvalid=nvalid)
            if nvalid ne 0 then begin
                x2=RTDyn(xEllipse[*,0],list.sf,list.hrange[0])
                y2=RTDyn(xEllipse[*,1],list.sf,list.vrange[0])
                
                if gapshift eq 0 then $
                    PlotS,transpose([[x2],[y2]]),/device,color=color,thick=list.dplotthick $
                else begin
                    PlotS,transpose([[x2[gapshift:*]],[y2[gapshift:*]]]),/device,color=color,thick=list.dplotthick
                    PlotS,transpose([[x2[0:gapshift-1]],[y2[0:gapshift-1]]]),/device,color=color,thick=list.dplotthick
                endelse
            endif
        endfor

        ; plot lines
        x2=RTDyn(xp,list.sf,list.hrange[0])
        y2=RTDyn(yp,list.sf,list.vrange[0])
        PlotS,transpose([[x2[0:1]],[y2[0:1]]]),/device,color=color,thick=list.dplotthick
        PlotS,transpose([[x2[2:3]],[y2[2:3]]]),/device,color=color,thick=list.dplotthick

       ;plot numbers
       if arr[2] eq 1 then xyouts,x2[0],y2[0],str,/device,charsize=list.letter,color=color,charthick=list.dplotthick
       endcase
    3: begin
           tt=BraggDtoX(arr[[3,4]],list.lambda,/onlyt)
           if tt[0] eq 0 then tt=tt[1]

        ; Calc phi range
        phi=azDetRange(0,2*!pi,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

        ; plot 2 lines
        for j=0l,n_elements(tt)-1 do begin
            xEllipse=EllipsePlotCoord(tt[j],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
            if nvalid ne 0 then begin
                x2=RTDyn(xEllipse[*,0],list.sf,list.hrange[0])
                y2=RTDyn(xEllipse[*,1],list.sf,list.vrange[0])
                PlotS,transpose([[x2],[y2]]),/device,color=color,thick=list.dplotthick
            endif
        endfor

       ;plot numbers
       if arr[2] eq 1 then begin
       xya=BraggXtoXY(tt,tt*0,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,1)
       x2=RTDyn(xya[0,0],list.sf,list.hrange[0])
       y2=RTDyn(xya[0,1],list.sf,list.vrange[0])
       xyouts,x2,y2,str,/device,charsize=list.letter,color=color,charthick=list.dplotthick
       endif
       endcase
    endcase

endif
endfor
endif

;----Display mark ellipse----
if list.smellipse ne 0 then begin
    LoadctRev,3,/silent
    tt=list.ttmellipse
    phi=azDetRange(0,2*!pi,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

    if (list.smellipse and 1) eq 1 then begin
        xEllipse=EllipsePlotCoord(tt,phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
        if nvalid ne 0 then begin
            x2=RTDyn(xEllipse[*,0],list.sf,list.hrange[0])
            y2=RTDyn(xEllipse[*,1],list.sf,list.vrange[0])
            PlotS,transpose([[x2],[y2]]),/device,color=10*16,thick=list.dplotthick
        endif
    endif

    if (list.smellipse and 2) eq 2 then begin
        xAxis=SemiAxisCoord(tt,list.a,list.b,list.dist,list.scrx,list.scry,list.center)
        xx2=RTDyn(xAxis[*,0],list.sf,list.hrange[0])
        yy2=RTDyn(xAxis[*,1],list.sf,list.vrange[0])
        Plots,transpose([[xx2[0:1]],[yy2[0:1]]]),/device,color=16*10,linestyle=2,thick=list.dplotthick
        Plots,transpose([[xx2[2:3]],[yy2[2:3]]]),/device,color=16*10,thick=list.dplotthick
    endif
endif

;----Display Tie points----
if total(list.ShowTie) ne 0 and (*list.nrTie)[0] ne 0 then begin
LoadctRev,1,/silent
i=0l
; Input Tie points
if list.ShowTie[0] then begin
if list.nrTieRing ne 0 then $
    for i=0l,total((*list.nrTie)[0:list.nrTieRing-1])-1 do begin
       x=RTDyn((*list.TieI)[0,i],list.sf,list.hrange[0])
       y=RTDyn((*list.TieI)[1,i],list.sf,list.vrange[0])
       Plots,[x-list.CW2,x+list.CW2],y,/device,color=7*16,thick=list.dplotthick
       Plots,x,[y-list.CW2,y+list.CW2],/device,color=7*16,thick=list.dplotthick
    endfor
for j=0+i,total((*list.nrTie)[list.nrTieRing])+i-1 do begin
    x=RTDyn((*list.TieI)[0,j],list.sf,list.hrange[0])
    y=RTDyn((*list.TieI)[1,j],list.sf,list.vrange[0])
    Plots,[x-list.CW2,x+list.CW2],y,/device,color=10*16,thick=list.dplotthick
    Plots,x,[y-list.CW2,y+list.CW2],/device,color=10*16,thick=list.dplotthick
endfor
endif
if list.TieGrid then begin
    LoadctRev,3,/silent
    ; Output Tie points
    if list.ShowTie[1] then begin
        if list.ShowTie[3] then begin
               ind=ConnectGrid(list.TieO,$
                list.GridDim,list.GridOffset,list.GridOrigin,list.GridTilt,list.GridSpReal,list.GridSpPix,$
                count=ct,ngaps=ngaps,gapcoord=gapcoord)
            x=RTDyn((*list.TieO)[0,*],list.sf,list.hrange[0])
            y=RTDyn((*list.TieO)[1,*],list.sf,list.vrange[0])
            for i=0l,ct-1 do Plots,x[[ind[i,0],ind[i,1]]],$
                  y[[ind[i,0],ind[i,1]]],/device,color=10*16,thick=list.dplotthick
            if ngaps ne 0 then begin
            gapcoord[*,0:1]=RTDyn(gapcoord[*,0:1],list.sf,list.hrange[0])
            gapcoord[*,2:3]=RTDyn(gapcoord[*,2:3],list.sf,list.vrange[0])
            for i=0l,ngaps-1 do Plots,gapcoord[i,0:1],$
                  gapcoord[i,2:3],/device,color=10*16,linestyle=1,thick=list.dplotthick
            endif
        endif else begin
        for j=0l,total((*list.nrTie)[0])-1 do begin
            x=RTDyn((*list.TieO)[0,j],list.sf,list.hrange[0])
            y=RTDyn((*list.TieO)[1,j],list.sf,list.vrange[0])
            Plots,[x-list.CW2,x+list.CW2],y,/device,color=10*16,thick=list.dplotthick
            Plots,x,[y-list.CW2,y+list.CW2],/device,color=10*16,thick=list.dplotthick
        endfor
        endelse
    endif
    ; Connect I-O
    if list.ShowTie[2] then $
     for j=0l,total((*list.nrTie)[0])-1 do begin
        point=[[(*list.TieI)[*,j]], [(*list.TieO)[*,j]]]
        x=RTDyn(point[0,*],list.sf,list.hrange[0])
        y=RTDyn(point[1,*],list.sf,list.vrange[0])
        Plots,x,y,/device,color=10*16,thick=list.dplotthick
    endfor
endif
endif

;----Display Tilt points----
if list.ShowTilt eq 1 and (*list.nrTilt)[0] ne 0 then begin
LoadctRev,3,/silent
i=0l
if list.nrTiltRing ne 0 then $
    for i=0l,total((*list.nrTilt)[0:list.nrTiltRing-1])-1 do begin
       x=RTDyn((*list.TiltI)[0,i],list.sf,list.hrange[0])
       y=RTDyn((*list.TiltI)[1,i],list.sf,list.vrange[0])
       Plots,[x-list.CW2,x+list.CW2],y,/device,color=7*16,thick=list.dplotthick
       Plots,x,[y-list.CW2,y+list.CW2],/device,color=7*16,thick=list.dplotthick
    endfor
for j=0l+i,total((*list.nrTilt)[list.nrTiltRing])-1+i do begin
    x=RTDyn((*list.TiltI)[0,j],list.sf,list.hrange[0])
    y=RTDyn((*list.TiltI)[1,j],list.sf,list.vrange[0])
    Plots,[x-list.CW2,x+list.CW2],y,/device,color=10*16,thick=list.dplotthick
    Plots,x,[y-list.CW2,y+list.CW2],/device,color=10*16,thick=list.dplotthick
endfor
endif

;----Display PDF data----
if (list.PDFnTotal ne 0) then begin
    LoadctRev,ntableadd(1),/silent
    off=0
    for i=0l,list.nPDF-1 do begin
        ind=where((*list.PDFs)[off:off+(*list.PDFn)[i]-1] ne 0,ct)+off
        if ct ne 0 then begin
            tt=BraggDtoX((*list.PDFd)[ind],list.lambda,/onlyt)
            phi=azDetRange(0,!pi*2,list.phipixphi0,list.scrx,list.scry,list.a,list.b)

            ;plot circles
            if ((*list.PDFni)[i] ne 0) then PDFi=ScaleArr((*list.PDFi)[ind],1,3) $
            else PDFi=replicate(1,ct)

            ; Labels
            hkl=fix((*list.PDFhkl)[*,ind])
            mihkl=min(hkl)
            n=strlen(stringr(max(hkl)))>(strlen(stringr(mihkl)))
            label=string(hkl,format='(3(I'+stringr(n)+'))')
            code=(*list.PDFs)[ind]

            col=i+1
            PDFangle=45*((i-(i/8)*8)+1)/180.*!dpi
            ;PDFangle=[45,45,45]/180.*!dpi
            for j=0l,ct-1 do begin
               if (code[j] and 1) eq 1 then begin
                       xEllipse=EllipsePlotCoord(tt[j],phi,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
                       if nvalid ne 0 then begin
                        x2=RTDyn(xEllipse[*,0],list.sf,list.hrange[0])
                        y2=RTDyn(xEllipse[*,1],list.sf,list.vrange[0])
                           PlotS,transpose([[x2],[y2]]),/device,color=col,thick=list.dplotthick*PDFi[j]
                       endif
               endif
               if (code[j] and 2) eq 2 then begin
                       xEllipse=EllipsePlotCoord(tt[j],PDFangle,list.a,list.b,list.phipixphi0,list.dist,list.scrx,list.scry,list.center,/forceopen,nvalid=nvalid)
                       if nvalid ne 0 then begin
                           x2=RTDyn([xEllipse[*,0],xEllipse[*,0]+0.5*list.CW1],list.sf,list.hrange[0])
                           y2=RTDyn([xEllipse[*,1],xEllipse[*,1]+1.5*list.CW1],list.sf,list.vrange[0])
                           Plots,x2,y2,/device,color=col,thick=list.dplotthick
                           xyouts,x2[1],y2[1],label[j],color=col,/device,charsize=list.letter,charthick=list.dplotthick;,ORIENTATION=45
                       endif
               endif
            endfor
        endif
        off+=(*list.PDFn)[i]
    endfor
endif

; ----Display Axis----
if list.RadAxShow then begin
    LoadctRev,3,/silent
    col=16*10
    xx=list.RadAxMin+(list.RadAxMax-list.RadAxMin)/list.RadAxBins*lindgen(list.RadAxBins+1)
    phi=make_array(list.RadAxBins+1,value=list.RadAxAzimuth*!pi/180)
    out=BraggXtoXY(xx,phi,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,list.RadAxType,/angledeg,/pixm,/nonort)
    x=RTDyn(out[*,0],list.sf,list.hrange[0])
    y=RTDyn(out[*,1],list.sf,list.vrange[0])
    Plots,x[[0,list.RadAxBins]],y[[0,list.RadAxBins]],/device,color=col
    
    angle=(90-list.RadAxAzimuth)*!pi/180
    addx=list.CW2*cos(angle)
    addy=list.CW2*sin(angle)
    for i=0,list.RadAxBins do begin
        Plots,[x[i],x[i]-addx],[y[i],y[i]+addy],/device,color=col
        xyouts,x[i]+2.5*addx,y[i]-2.5*addy,string(xx[i],format=list.RadAxFormat),orientation=list.RadAxAzimuth,/device,color=col
    endfor
endif

tvlct,RoldCT,GoldCT,BoldCT
end; pro dplot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshScalingSlider,top
ID=widget_info(top,find_by_uname='drawslider')
widget_control,top,get_uvalue=list
widget_control,ID,get_uvalue=info,get_value=setid

wset,setid
LoadctRev,0,/silent,rev=list.revcol
erase,100

; Legend in full bitrange
b=2ll^(list.clipslider[0]*list.clipslider[2])-1
e=2ll^(list.clipslider[1]*list.clipslider[2])-1
norm=e-b
n=info.xs
inc=norm/(n-1)
img=b+inc*l64indgen(n)
img=rebin(img,info.xs,info.ys,/sample)

tv, WeakScl(img,list.sclmax,list.displaypower),info.xoff,info.yoff,/normal

; Lines
LoadctRev,39,/silent

mi=(2^(list.clipslider[3]*list.clipslider[2])-2^(list.clipslider[0]*list.clipslider[2]))/norm
ma=(2^(list.clipslider[4]*list.clipslider[2])-2^(list.clipslider[0]*list.clipslider[2]))/norm
m=(mi+ma)/2

plots,[mi,mi],[0,1],/normal,thick=1.5,col=250
plots,[m,m],[0,1],/normal,thick=1.5,col=200
plots,[ma,ma],[0,1],/normal,thick=1.5,col=50

; Refresh slider
ID=widget_info(top,FIND_BY_UNAME='clip slider min')
WIDGET_CONTROL, ID,SET_SLIDER_MIN=list.clipslider[0]
WIDGET_CONTROL, ID,SET_SLIDER_MAX=list.clipslider[1]
WIDGET_CONTROL, ID,SET_VALUE=list.clipslider[3]
ID=widget_info(top,FIND_BY_UNAME='clip slider max')
WIDGET_CONTROL, ID,SET_SLIDER_MIN=list.clipslider[0]
WIDGET_CONTROL, ID,SET_SLIDER_MAX=list.clipslider[1]
WIDGET_CONTROL, ID,SET_VALUE=list.clipslider[4]
ID=widget_info(top,FIND_BY_UNAME='clip slider power')
WIDGET_CONTROL, ID,SET_VALUE=round(list.displaypower*100)
end;pro RefreshScalingSlider
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SL_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent NE 'DOWN' THEN RETURN

widget_control,ev.top,get_uvalue=list
widget_control,ev.id,get_uvalue=info

norm=2^(list.clipslider[1]*list.clipslider[2])-2^(list.clipslider[0]*list.clipslider[2])
mi=(2^(list.clipslider[3]*list.clipslider[2])-2^(list.clipslider[0]*list.clipslider[2]))/norm
ma=(2^(list.clipslider[4]*list.clipslider[2])-2^(list.clipslider[0]*list.clipslider[2]))/norm
m=(mi+ma)/2

tmp=min((convert_coord(ev.x,ev.y,/device,/to_normal))[0]-[mi,m,ma],type,/abs)
info.type=type
Widget_Control, ev.id, draw_motion=1, Event_Pro='SL_DRAG_EVENTS',set_uvalue=info
end;pro SL_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SL_DRAG_EVENTS ,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]

widget_control,ev.top,get_uvalue=list
widget_control,ev.id,get_uvalue=info


case info.type of
0:    begin
    p=(convert_coord(ev.x,ev.y,/device,/to_normal))[0]*(2^(list.clipslider[1]*list.clipslider[2])-2^(list.clipslider[0]*list.clipslider[2]))
    p=alog(p+2^(list.clipslider[0]*list.clipslider[2]))/(alog(2)*list.clipslider[2])
    list.clipslider[3]=list.clipslider[0]>round(p)<list.clipslider[1]
    list.clipslider[4]>=list.clipslider[3]
    endcase
1:    begin
    norm=2^(list.clipslider[1]*list.clipslider[2])-2^(list.clipslider[0]*list.clipslider[2])
    mi=(2^(list.clipslider[3]*list.clipslider[2])-2^(list.clipslider[0]*list.clipslider[2]))/norm
    ma=(2^(list.clipslider[4]*list.clipslider[2])-2^(list.clipslider[0]*list.clipslider[2]))/norm
    p=(convert_coord(ev.x,ev.y,/device,/to_normal))[0]-(mi+ma)/2
    mi+=p
    ma+=p
    p=mi*norm
    p=alog(p+2^(list.clipslider[0]*list.clipslider[2]))/(alog(2)*list.clipslider[2])
    list.clipslider[3]=list.clipslider[0]>round(p)<list.clipslider[1]
    p=ma*norm
    p=alog(p+2^(list.clipslider[0]*list.clipslider[2]))/(alog(2)*list.clipslider[2])
    list.clipslider[4]=list.clipslider[0]>round(p)<list.clipslider[1]
    list.clipslider[3]<=list.clipslider[4]
    endcase
2:    begin
    p=(convert_coord(ev.x,ev.y,/device,/to_normal))[0]*(2^(list.clipslider[1]*list.clipslider[2])-2^(list.clipslider[0]*list.clipslider[2]))
    p=alog(p+2^(list.clipslider[0]*list.clipslider[2]))/(alog(2)*list.clipslider[2])
    list.clipslider[4]=list.clipslider[0]>round(p)<list.clipslider[1]
    list.clipslider[3]<=list.clipslider[4]
    endcase
endcase
list.sclmax=2.^(list.clipslider[3:4]*list.clipslider[2])-1

widget_control,ev.top,set_uvalue=list
RefreshDisplay, ev.top,list
RefreshScalingSlider,ev.top

IF thisEvent EQ 'UP' THEN BEGIN
    info.type=-1
    Widget_Control, ev.id, draw_motion=0, Event_Pro='SL_PROCESS_EVENTS',set_uvalue=info
endif

end;pro SL_DRAG_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshCircleSlider,top,list,update=update

; Middle of the dynamic screen in image coordinates
widget_control,list.drawdyn,get_draw_view=XY
XY=DynTR((XY>0)+[(list.stath<list.dynh)/2.,(list.statv<list.dynv)/2.],list.sf,[list.hrange[0],list.vrange[0]])

; azimuth for XY
list.Pdyn=XY
XY=XY-list.center
azi=atan(XY[1],XY[0])
if azi lt 0 then azi+=2*!dpi
list.cslider[3]=azi*180/!dpi/list.cslider[2]

; update slider
ID=widget_info(top,FIND_BY_UNAME='circle slider')
WIDGET_CONTROL, ID,SET_VALUE=list.cslider[3]
ID=widget_info(top,FIND_BY_UNAME='csliderlabel')
WIDGET_CONTROL, ID,SET_VALUE='Azimuth: '+stringr(list.cslider[3]*list.cslider[2])+msymbols('degrees')

; dynamic marker plot
if list.bPdyn ne 0 then begin
    Wset,list.drawindex
    tplot,list
    Wset,list.pixindex1
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.drawindex]
    update=1B
endif else update=0B

end;pro RefreshCircleSlider
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshDisplay,top,list

if ptr_valid(list.tiff) eq 0 then return

update=0B
if widget_info(list.drawdyn,/VALID_ID) then begin
    Wset,list.drawdynindex
    dplot,list
    wset,list.pixindex2
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.drawdynindex]
    RefreshCircleSlider,top,list,update=update
endif

if update ne 1 then begin
    Wset,list.drawindex
    tplot,list
    Wset,list.pixindex1
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.drawindex]
endif

if total((list.CorDone eq 0)*list.CorMethod) ne 0 then addstr=' >> NOT CORRECTED!' else addstr=''
widget_control,top,TLB_SET_TITLE='XRDUA: '+list.file+((list.xdinr gt 0)?('_'+stringr(list.xdinr)):'')+addstr

Widget_Control, top, Set_UValue=list
end;pro RefreshDisplay
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MIDDLE_PROCESS_EVENTS, ev
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN
widget_control,ev.top,get_uvalue=list
list.middle=[ev.x,ev.y]
Widget_Control, ev.top, Set_UValue=list
wset,list.drawindex
end;pro MIDDLE_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro mode, ev, val, refresh=refresh
widget_control,ev.top,update=0
widget_control,ev.top,get_uvalue=list
strmat=['Beam position','Mask Setup','Background','Zoom','Spatial Distortion','Calibration','Azimuthal offset','Debye Marker','Zingers','Saturation','Flat Field','Mask off']
ind=where(strmat ne val,count)
if count ne n_elements(strmat) then begin
    list.mode=val
    widget_control,ev.top,set_uvalue=list
endif
if count ne 0 then strmat=strmat[ind]
widget_control,list.draw,draw_button=0,draw_motion=0
booldrawdyn=widget_info(list.drawdyn,/VALID_ID)
if booldrawdyn then widget_control,list.drawdyn,draw_button=0,draw_motion=0

case val of
    'Zoom': begin
             widget_control,list.draw,Event_Pro='ZOOM_PROCESS_EVENTS',draw_button=1

               text=strarr(2)
               text[0]='X:'+stringr(list.hrange[0])+' - '+stringr(list.hrange[1])
               text[1]='Y:'+stringr(list.vrange[0])+' - '+stringr(list.vrange[1])

            ID3=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID3 ne 0 then widget_control,ID3,/destroy
            ID3=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID3,/row,uname='Mbase')
            textw=widget_text(base,xsize=10,ysize=2,uname='coord',value=text)
            endcase
    'Beam position':begin

            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/row,uname='Mbase')
            base1=widget_base(base,/column)
            button0=widget_button(base1,value='Fit Center',uname='Fit Center')
            button1=widget_button(base1,value='Use Mouse',uname='Use Mouse')
            button2=widget_button(base1,value='3 Points',uname='3 Points')
            center=list.center
            str=stringr(center)
            text1=widget_text(base1,xsize=20,/editable,value=str[0]+','+str[1],uvalue='center',uname='center')
            base2=widget_base(base,/column)
            P1=list.P1
            str=stringr(P1)
            text2=widget_text(base2,xsize=10,value=str[0]+','+str[1],uname='P1')
            P2=list.P2
            str=stringr(P2)
             text3=widget_text(base2,xsize=10,value=str[0]+','+str[1],uname='P2')
            P3=list.P3
            str=stringr(P3)
            text4=widget_text(base2,xsize=10,value=str[0]+','+str[1],uname='P3')
            button3=widget_button(base2,value='Calc')
            base3=widget_base(base,/column)
            button4=widget_button(base3,value='modify...',uvalue='bP1',uname='bP1')
            button5=widget_button(base3,value='modify...',uvalue='bP2',uname='bP2')
            button6=widget_button(base3,value='modify...',uvalue='bP3',uname='bP3')
            endcase
    'Use Mouse':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='usemouse_events',$
                draw_button=1
            widget_control,list.draw,Event_Pro='usemouse_events',$
                draw_button=1
            strmat=['3 Points','Fit Center']
            endcase
    '3 Points':   begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='circle_events',$
              draw_button=1
            widget_control,list.draw,Event_Pro='circle_events',$
              draw_button=1
            strmat=['Use Mouse','Fit Center']
            endcase
    'Fit Center': begin
            surf_int,list,1,ev.top,'Select area of Beam-Center'
            strmat=['Use Mouse','3 Points']
            refresh=1b
            endcase
    'Mask Setup':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy

            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/row,uname='Mbase',uvalue=fltarr(5))
            base1=widget_base(base,/column)
            label1=widget_label(base1,value='Define area:')
            button1=widget_button(base1,value='Rectangle',uname='Rectangle')
            
            button2=widget_button(base1,value='Ellipse',uname='Ellipse')
            button22=widget_button(base1,value='Ellipse (0,0)',uname='Ellipse (0,0)')
            button3=widget_button(base1,value='Arc',uname='Arc')
            
            base2=widget_base(base,/column)
            label2=widget_label(base2,value=' ')
            button3=widget_button(base2,value='Whole Image',uname='Whole Image')
            button3=widget_button(base2,value='Full Arc',uname='Full Arc')
            button=widget_button(base2,value='Range',uname='Range')
            
            base2=widget_base(base,/column)
            label2=widget_label(base2,value='Perform:')
            button4=widget_button(base2,value='Set',uname='Set')
            button6=widget_button(base2,value='Edit Mask',uname='Edit Mask')
            endcase
    'Arc':  begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='ARC_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='ARC_PROCESS_EVENTS',$
            draw_button=1
            strmat=['Rectangle','Ellipse','Ellipse (0,0)']
            endcase
    'Ellipse':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='CIRCLE_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='CIRCLE_PROCESS_EVENTS',$
            draw_button=1
            list.circle2=1
            widget_control,ev.top,set_uvalue=list
            strmat=['Rectangle','Arc','Ellipse (0,0)']
            endcase
    'Ellipse (0,0)':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='CIRCLE_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='CIRCLE_PROCESS_EVENTS',$
            draw_button=1
            list.circle2=0
            widget_control,ev.top,set_uvalue=list
            strmat=['Ellipse','Arc','Rectangle']
            endcase
    'Rectangle':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='SQUARE_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='SQUARE_PROCESS_EVENTS',$
            draw_button=1
            strmat=['Ellipse','Arc','Ellipse (0,0)']
            endcase
    'Background':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/column,uname='Mbase')

            baset=widget_base(base,/row)
            bool=1b;CMethodAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'back')
            base3=widget_base(baset,/row,/exclusive,sensitive=bool,uname='exclbase')
                button1=widget_button(base3,value='Background Stripping',uvalue=1)
                button2=widget_button(base3,value='Dark Image subtraction',uvalue=2)
                button0=widget_button(base3,value='Already Subtracted',uvalue=0)
            base3=widget_base(baset,/row,/nonexclusive,sensitive=bool,uname='nonexclbase')
                button3=widget_button(base3,value='Lower threshold',uvalue=3)
                button4=widget_button(base3,value='Higher threshold',uvalue=3)

            base4=widget_base(base,/row)
            base1=widget_base(base4,/column,uname='BackBase1')
                button=widget_button(base1,value='Strip')
                label=widget_label(base1,value='Width of filter:')
                text=widget_text(base1,value=stringr(list.wb),$
                        uname='wstrip',uvalue='wstrip',/editable)
                label=widget_label(base1,value='Number of iterations:')
                text=widget_text(base1,value=stringr(list.niter),$
                        uname='istrip',uvalue='istrip',/editable)

            base2=widget_base(base4,/column,uname='BackBase2')
              button=widget_button(base2,value='Load Dark Image')
              label=widget_label(base2,value='Smoothing order:')
                text=widget_text(base2,value=stringr(list.bcko),$
                        uname='bcko',uvalue='bcko',/editable)
                label=widget_label(base2,value='Smoothing width:')
                text=widget_text(base2,value=stringr(list.bckw),$
                        uname='bckw',uvalue='bckw',/editable)
                label=widget_label(base2,value='Multiplication factor:')
                text=widget_text(base2,value=stringr(list.bckm),$
                        uname='bckm',uvalue='bckm',/editable)
            base5=widget_base(base4,/column)

              bool=CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'back',Rest=bool2)
              button=widget_button(base5,value='Subtract backgr.',uname='Subtract backgr.',sensitive=bool)
              button=widget_button(base5,value='Reset backgr.',uname='Reset backgr.',sensitive=bool2)
              button=widget_button(base5,value='Show backgr.')
              button=widget_button(base5,value='Save backgr.')
            base6=widget_base(base4,/row)
              label=widget_label(base6,value='Lower:')
              text1=widget_text(base6,value=stringr(list.thres[0]),$
                  uname='tstrip1',uvalue='tstrip1',/editable)
              label=widget_label(base6,value='Higher:')
              text2=widget_text(base6,value=stringr(list.thres[1]),$
                  uname='tstrip2',uvalue='tstrip2',/editable)

            bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])

            case bkcode[0] of
            0:   begin
            widget_control,button0,set_button=1
            widget_control,base1,sensitive=0
            widget_control,base2,sensitive=0
            endcase
            1:   begin
            widget_control,button1,set_button=1
            widget_control,base2,sensitive=0
            endcase
            2:   begin
            widget_control,button2,set_button=1
            widget_control,base1,sensitive=0
            endcase
            endcase

            if bkcode[1] eq 1 then widget_control,button3,set_button=1
            if bkcode[2] eq 1 then widget_control,button4,set_button=1
            widget_control,text1,sensitive=bkcode[1]
            widget_control,text2,sensitive=bkcode[2]

            endcase
    'Subtract backgr.':begin
            strmat=['Reset backgr.']
            endcase
    'Reset backgr.':begin
            if CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'back') then strmat=['Subtract backgr.']
            ID=widget_info(ev.top,FIND_BY_UNAME='exclbase')
            widget_control,ID,sensitive=1
            ID=widget_info(ev.top,FIND_BY_UNAME='nonexclbase')
            widget_control,ID,sensitive=1
            endcase
    'Debye Marker':begin
            widget_control,list.draw,Event_Pro='ELLIPSEM_PROCESS_EVENTS',$
            draw_button=1
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='ELLIPSEM_PROCESS_EVENTS',$
            draw_button=1
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/row,uname='Mbase')

            base3=widget_base(base,/row)
            xs=130

            base5=widget_base(base3,/column)
            res=BraggTtoX(list.ttmellipse,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/pixm)
            baset=widget_base(base5,/row)
            label=widget_label(baset,value='d-Spacing ('+msymbols('angstroms')+'):',xsize=xs)
            text=widget_text(baset,/editable,value=stringr(res[1]),$
              uname='MEllipsed',uvalue='MEllipsed')
            
            baset=widget_base(base5,/row)
            label=widget_label(baset,value='Q-space (1/nm):',xsize=xs)
            text=widget_text(baset,/editable,value=stringr(res[2]),$
              uname='MEllipseq',uvalue='MEllipseq')

            baset=widget_base(base5,/row)
            label=widget_label(baset,value='2-theta ('+msymbols('degrees')+'):',xsize=xs)
            text=widget_text(baset,/editable,value=stringr(list.ttmellipse*180./!dpi),$
              uname='MEllipset',uvalue='MEllipset')

            baset=widget_base(base5,/row)
            label=widget_label(baset,value='Flat Rad. dist. (mm):',xsize=xs)
            text=widget_text(baset,/editable,value=stringr(res[0]),$
              uname='MEllipsefr',uvalue='MEllipsefr')

            tmp=ConicTypeTT(list.a,list.ttmellipse,str=vstr)
            label=widget_label(base5,value=vstr,uname='warning',xsize=xs)
            
            base5=widget_base(base3,/column)
            
            baset=widget_base(base5,/row)
            label=widget_label(baset,value='Position (pixels):',xsize=xs)
            text=widget_text(baset,uname='MEllipsep',uvalue='MEllipsep')

            baset=widget_base(base5,/row)
            label=widget_label(baset,value='Rad. dist. (mm):',xsize=xs)
            text=widget_text(baset,uname='MEllipser',uvalue='MEllipser')

            baset=widget_base(base5,/row)
            label=widget_label(baset,value='Intensity ():',xsize=xs)
            text=widget_text(baset,uname='MEllipsei',uvalue='MEllipsei')

            baset=widget_base(base5,/row)
            label=widget_label(baset,value='Azimuth ('+msymbols('degrees')+'):',xsize=xs)
            text=widget_text(baset,uname='MEllipsea',uvalue='MEllipsea')
            
            endcase
    'Spatial Distortion':begin
            if (list.TieGrid ne 2) then begin
                if booldrawdyn then widget_control,list.drawdyn,Event_Pro='SQUARE_PROCESS_EVENTS',$
                draw_button=1
                widget_control,list.draw,Event_Pro='SQUARE_PROCESS_EVENTS',$
                draw_button=1
            endif
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/column,uname='Mbase')
            baset=widget_base(base,/column,/nonexclusive)
                button=widget_button(baset,value='Fast mode')
                widget_control,button,set_button=list.BFastSpline

            baset=widget_base(base,/row,/exclusive)
            buttong1=widget_button(baset,value='With Grid')
            buttong2=widget_button(baset,value='Without Grid')
            buttong3=widget_button(baset,value='Spline file')
            case list.TieGrid of
            1:widget_control,buttong1,set_button=1
            0:widget_control,buttong2,set_button=1
            2:widget_control,buttong3,set_button=1
            endcase


            baset=widget_base(base,/row)
            base3=widget_base(baset,/column,uname='WithGridBase',sensitive=list.TieGrid eq 1)
            button=widget_button(base3,value='Find Tie Points')
            button=widget_button(base3,value='Grid Properties')

            base1=widget_base(baset,/column,uname='WithoutGridBase',sensitive=list.TieGrid eq 0)
            button7=widget_button(base1,value='Auto Set ',uname='Auto Set ')
            button7=widget_button(base1,value='Auto Set with PDF ',uname='Auto Set with PDF ')
            button6= widget_button(base1,value='Delete Ring',uname='Delete Ring')
            button6= widget_button(base1,value='Edit d-spacing ',uname='Edit d-spacing ')
            button1=widget_button(base1,value='Close Ring')

            base2=widget_base(baset,/column,uname='SPbase',sensitive=list.TieGrid eq 2)
            button=widget_button(base2,value='Load spline file',uname='Load spline file')
            basett=widget_base(base2,/column,/nonexclusive)
                button=widget_button(basett,value='Flip vertical')
                widget_control,button,set_button=list.vFlipSpline


            base2=widget_base(baset,/column,uname='GridBase')
            button6= widget_button(base2,value='Select area',uname='Select area',sensitive=(list.TieGrid ne 2))
            widget_control,button6,set_button=list.TieGrid ne 2
            button3=widget_button(base2,value='Add Tie Point',uname='Add Tie Point',sensitive=(list.TieGrid ne 2))
            button2= widget_button(base2,value='Delete Tie Point',uname='Delete Tie Point',sensitive=(list.TieGrid ne 2))

            button=widget_button(base2,value='Delete All Tie Points',uname='Delete All Tie Points',sensitive=(list.TieGrid ne 2))
            bool=CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'spatial',Rest=bool2)
            button4=widget_button(base2,value='Correct Spatial',uname='Correct Spatial',sensitive=bool)
            button5=widget_button(base2,value='Restore Spatial',uname='Restore Spatial',sensitive=bool2)

            endcase
    'Calibration': begin
;            widget_control,list.draw,Event_Pro='TILT_PROCESS_EVENTS',$
;            draw_button=1
;            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='TILT_PROCESS_EVENTS',$
;            draw_button=1

            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/row,uname='Mbase')
            base1=widget_base(base,/column)
            button8=widget_button(base1,value='Auto Set with PDF',uname='Auto Set with PDF')
            button8=widget_button(base1,value='Auto Set',uname='Auto Set')
            button3=widget_button(base1,value='Set Point',uname='Set Point')
            button1=widget_button(base1,value='Close TiltRing',uname='Close TiltRing')
            button6= widget_button(base1,value='Edit d-spacing',uname='Edit d-spacing')
 
            base2=widget_base(base,/column)
            button2= widget_button(base2,value='Delete Point',uname='Delete Point')
            button7= widget_button(base2,value='Delete TiltRing',uname='Delete TiltRing')
            button1=widget_button(base2,value='Delete small',uname='Delete small')
            button9= widget_button(base2,value='Calibrate',uname='Calibrate')

;            base3=widget_base(base,/column,/exclusive,sensitive=bool)
;            buttonNN=widget_button(base3,value='Nearest Neighbor')
;            buttonB=widget_button(base3,value='Bilinear')
;            buttonC=widget_button(base3,value='Cubic')
;            case list.mipol of
;            1: widget_control,buttonNN,set_button=1
;            2: widget_control,buttonB,set_button=1
;            3: widget_control,buttonC,set_button=1
;            endcase
            endcase
    'Delete Tie Point':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='DELTIE_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='DELTIE_PROCESS_EVENTS',$
            draw_button=1
            strmat=['Add Tie Point','Select area']
            if not list.TieGrid then strmat=[strmat,'Edit d-spacing ']
            endcase
    'Add Tie Point':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='TIE_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='TIE_PROCESS_EVENTS',$
            draw_button=1
            strmat=['Delete Tie Point','Select area']
            if not list.TieGrid then strmat=[strmat,'Edit d-spacing ']
            endcase
    'Select area':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='SQUARE_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='SQUARE_PROCESS_EVENTS',$
            draw_button=1
            strmat=['Delete Tie Point','Add Tie Point']
            if not list.TieGrid then strmat=[strmat,'Edit d-spacing ']
            endcase
    'Delete Point':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='DELTILT_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='DELTILT_PROCESS_EVENTS',$
            draw_button=1
            strmat=['Set Point','Edit d-spacing']
            endcase
    'Set Point':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='TILT_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='TILT_PROCESS_EVENTS',$
            draw_button=1
            strmat=['Delete Point','Edit d-spacing']
            endcase
    'Edit d-spacing':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='DTILT_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='DTILT_PROCESS_EVENTS',$
            draw_button=1
            strmat=['Delete Point','Set Point']
            endcase
    'Azimuthal offset':begin
        
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/row,uname='Mbase')
            
            xs=100
            
            base3=widget_base(base,/column)
            
            baset=widget_base(base3,/row)
            label=widget_label(baset,value='',xsize=xs)
            label=widget_label(baset,value='Azimuth in degrees: ',xsize=xs)

            baset=widget_base(base3,/row)
            label=widget_label(baset,value='Image frame:',xsize=xs)
            text=widget_text(baset,uname='phipix',uvalue='phipix')
            
            baset=widget_base(base3,/row)
            label=widget_label(baset,value='Detector frame:',xsize=xs)
            text=widget_text(baset,uname='phid',uvalue='phid')
            
            baset=widget_base(base3,/row)
            label=widget_label(baset,value='Orthogonal frame:',xsize=xs)
            text=widget_text(baset,uname='phi',uvalue='phi')
            
            baset=widget_base(base3,/row)
            label=widget_label(baset,value='Azimuthal offset:',xsize=xs)
            text=widget_text(baset,value=stringr(list.phipixphi0*180/!pi),uname='phipixphi0',uvalue='phipixphi0',/editable)
            
            base3=widget_base(base,/column)
            
            baset=widget_base(base3,/row)
            label=widget_label(baset,value='',xsize=xs)
            label=widget_label(baset,value='Azimuth in degrees: ',xsize=xs)

            baset=widget_base(base3,/row)
            label=widget_label(baset,value='Reference:',xsize=xs)
            text=widget_text(baset,value='0',uname='phi0',uvalue='phi0',/editable)

            button2= widget_button(base3,value='Azimuthal default',uname='Azimuthal default')
        
        
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='AZIOFF_DRAW',$
            draw_button=1,draw_motion=1
            widget_control,list.draw,Event_Pro='AZIOFF_DRAW',$
            draw_button=1,draw_motion=1
            endcase
    'Edit d-spacing ':begin
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='DTIE_PROCESS_EVENTS',$
            draw_button=1
            widget_control,list.draw,Event_Pro='DTIE_PROCESS_EVENTS',$
            draw_button=1
            strmat=['Delete Tie Point','Add Tie Point','Select area']
            endcase
'Correct Spatial':begin
            strmat=['Restore Spatial']
            endcase
'Restore Spatial':begin
            if CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'spatial') then strmat=['Correct Spatial']
            endcase
'Zingers':    begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy

            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/row,uname='Mbase',uvalue=fltarr(5))
            base1=widget_base(base,/column)
            
            bool=CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'zinger',Rest=bool2)
            
            baset=widget_base(base1,/row)
            label=widget_label(baset,value='Zinger width:')
            text=widget_text(baset,value=stringr(list.ZingerWidth),uvalue='zingerwidth',/editable)
            baset=widget_base(base1,/row)
            label=widget_label(baset,value='Image/Smoothed threshold:')
            text=widget_text(baset,value=stringr(list.ZingerThreshold),uvalue='zingerthres',/editable)
            
            button=widget_button(base1,value='Remove Zingers',uname='Remove Zingers',sensitive=bool)
            button=widget_button(base1,value='Restore Zingers',uname='Restore Zingers',sensitive=bool2)
            endcase
'Remove Zingers': strmat=['Restore Zingers']
'Restore Zingers': if CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'zinger') then strmat=['Remove Zingers']
'Saturation':    begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy

            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/row,uname='Mbase',uvalue=fltarr(5))
            base1=widget_base(base,/column)

            bool=CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'satur',Rest=bool2)
            
            baset=widget_base(base1,/row)
                label=widget_label(baset,value='Higher Threshold:')
                text=widget_text(baset,value=stringr(list.SaturMax),uname='saturmax',uvalue='saturmax',/editable)
            baset=widget_base(base1,/row)
                label=widget_label(baset,value='Minimal aspect ratio:')
                text=widget_text(baset,value=stringr(list.SaturAspect),uname='saturaspect',uvalue='saturaspect',/editable)
            baset=widget_base(base1,/row)
                label=widget_label(baset,value='Minimal area:')
                text=widget_text(baset,value=stringr(list.SaturNMax),uname='saturnmax',uvalue='saturnmax',/editable)
            button=widget_button(base1,value='Remove Saturation',uname='Remove Saturation',sensitive=bool)
            button=widget_button(base1,value='Restore Saturation',uname='Restore Saturation',sensitive=bool2)
            endcase
'Remove Saturation': strmat=['Restore Saturation']
'Restore Saturation': if CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'satur') then strmat=['Remove Saturation']
'Flat Field':begin
             ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy

            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/row,uname='Mbase',uvalue=fltarr(5))
            base1=widget_base(base,/column)

            bool=CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'ff',Rest=bool2)
            button=widget_button(base1,value='Load Flat-Field Image',sensitive=~ptr_valid(list.ff),uname='Load Flat-Field Image')
            button=widget_button(base1,value='Reset Flat-Field Image',sensitive=ptr_valid(list.ff),uname='Reset Flat-Field Image')
            
            base1=widget_base(base,/column)
            button=widget_button(base1,value='Apply Flat-Field Correction',uname='Apply Flat-Field Correction',sensitive=bool)
            button=widget_button(base1,value='Restore Flat-Field',uname='Restore Flat-Field',sensitive=bool2)
            
            base1=widget_base(base,/column)
            
            xs=150
            baset=widget_base(base1,/row)
            label=widget_label(baset,value='Type of flood field: ',xsize=xs)
            drop=widget_droplist(baset,value=['None','Spherical','Plane'],uname='fftype',uvalue='fftype')
            widget_control,drop,SET_DROPLIST_SELECT=list.IntCor2Dinfo.fftype
           
            baset=widget_base(base1,/row)
            label=widget_label(baset,value=LabelPff(list.IntCor2Dinfo.fftype),uname='labelPff',xsize=xs)
            text=widget_text(baset,value=stringr(list.IntCor2Dinfo.Pff),uname='Pff',uvalue='Pff',/editable)
            
             endcase
'Apply Flat-Field Correction': strmat=['Restore Flat-Field']
'Restore Flat-Field': if CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'ff') then strmat=['Apply Flat-Field Correction']
'Load Flat-Field Image':strmat=['Reset Flat-Field Image']
'Reset Flat-Field Image':strmat=['Load Flat-Field Image']
'Mask off':begin
             ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy

            ID2=widget_info(ev.top,FIND_BY_UNAME='dynbase')
            base=widget_base(ID2,/row,uname='Mbase',uvalue=fltarr(5))
            base0=widget_base(base,/row)
            
            base2=widget_base(base0,/column)
            baset=widget_base(base2,/column,/nonexclusive)
            button=widget_button(baset,value='Valid pixels:')
            widget_control,button,set_button=list.maskoffdyn
            text=widget_text(base2,value=list.maskoffdynthres,uname='maskoffdynthres',uvalue='maskoffdynthres',editable=list.maskoffdyn)
            
            base1=widget_base(base0,/column)
            bool=CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'maskoff',Rest=bool2)
            button=widget_button(base1,value='Load Mask Image',uname='Load Mask Image',sensitive=~ptr_valid(list.maskoff))
            button=widget_button(base1,value='Reset Mask Image',uname='Reset Mask Image',sensitive=ptr_valid(list.maskoff))
            
            base3=widget_base(base0,/column)
            button=widget_button(base3,value='Apply Mask',uname='Apply Mask',sensitive=bool)
            button=widget_button(base3,value='Restore Mask',uname='Restore Mask',sensitive=bool2)
            
            DrawPolygon,list
            if booldrawdyn then widget_control,list.drawdyn,Event_Pro='AddPolygon_EVENTS', draw_button=1,draw_keyboard_events=1,/INPUT_FOCUS
            widget_control,list.draw,Event_Pro='AddPolygon_EVENTS', draw_button=1,draw_keyboard_events=1,/INPUT_FOCUS
             endcase
'Apply Mask': strmat=['Restore Mask']
'Restore Mask': if CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'maskoff') then strmat=['Apply Mask']
'Load Mask Image': strmat=['Reset Mask Image']
'Reset Mask Image': strmat=['Load Mask Image']
    else:
endcase
if not keyword_set(refresh) then begin
    setsens,ev.top,strmat,sen=sen
    widget_control,ev.id,sensitive=0
endif
widget_control,ev.top,update=1
end ;pro mode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_Expd, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=IDtop
widget_control,IDtop,sensitive=1
end;pro CleanUp_Expd
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Expd_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_uvalue=uval
widget_control,ev.id,get_value=val
widget_control,ev.top,get_uvalue=top
widget_control,top,get_uvalue=list
case uval of
    'OK' :  begin
            ID=widget_info(ev.top,find_by_uname='a')
            widget_control,ID,get_value=val
            list.a=float(val)*!pi/180.
            ID=widget_info(ev.top,find_by_uname='b')
            widget_control,ID,get_value=val
            list.b=float(val)*!pi/180.
            ID=widget_info(ev.top,find_by_uname='scrx')
            widget_control,ID,get_value=val
            list.scrx=float(val)*10^(-6.)
            ID=widget_info(ev.top,find_by_uname='scry')
            widget_control,ID,get_value=val
            list.scry=float(val)*10^(-6.)
            ID=widget_info(ev.top,find_by_uname='dist')
            widget_control,ID,get_value=val
            list.dist=float(val)*0.01
            ID=widget_info(ev.top,find_by_uname='xcen')
            widget_control,ID,get_value=val
            list.center[0]=float(val)
            ID=widget_info(ev.top,find_by_uname='ycen')
            widget_control,ID,get_value=val
            list.center[1]=float(val)
            ID=widget_info(ev.top,find_by_uname='xlambda')
            widget_control,ID,get_value=val
            list.lambda=float(val)
            RefreshDisplay,top,list
            WIDGET_CONTROL, ev.top, /DESTROY
            return
            endcase
    'a':   list.a=float(val)*!dpi/180. ;from deg to rad
    'b':   list.b=float(val)*!dpi/180. ;from deg to rad
    'xlambda':begin
         ID=widget_info(ev.top,FIND_BY_UNAME='xen')
         list.lambda=float(val)
         widget_control,ID,set_value=stringr(list.cconv/list.lambda)
         endcase
    'xen': begin
         ID=widget_info(ev.top,FIND_BY_UNAME='xlambda')
         list.lambda=list.cconv/float(val)
         widget_control,ID,set_value=stringr(list.lambda)
         endcase
    'scrx': list.scrx=float(val)*10^(-6.) ; from um to m
    'scry': list.scry=float(val)*10^(-6.) ; from um to m
    'dist': list.dist=float(val)*0.01 ; from cm to m
    'xcen': list.center[0]=float(val)
    'ycen': list.center[1]=float(val)
endcase
RefreshDisplay,top,list
end; pro Expd_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ZOOM_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent NE 'DOWN' THEN RETURN

possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
bPress = possibleButtons[ev.press]

widget_control,ev.top,get_uvalue=list
if (bPress eq 'RIGHT') then begin
    if (list.zoom eq 0) then begin ;move scrollbar
        x=StatTR(ev.x,list.sfh)
        y=StatTR(ev.y,list.sfv)
        xy=[RTDyn(x,list.sf,list.hrange[0]),RTDyn(y,list.sf,list.vrange[0])]
        WIDGET_CONTROL, list.drawdyn, SET_DRAW_VIEW=[(xy[0]-list.stath/2.)>0,(xy[1]-list.statv/2.)>0]
        RefreshCircleSlider,ev.top,list
        Widget_Control, ev.top, Set_UValue=list
    endif else begin ;move zoom area but never resize it
        m=(list.hrange[1]-list.hrange[0])/2.
        list.hrange=StatTR(ev.x,list.sfh)+[-m,m]
        tmp=-list.hrange[0]
        if (tmp gt 0) then list.hrange=list.hrange+tmp
        tmp=list.hrange[1]-(*list.tiffs)[1]+1
        if (tmp gt 0) then list.hrange=list.hrange-tmp

        m=(list.vrange[1]-list.vrange[0])/2.
        list.vrange=(StatTR(ev.y,list.sfv)+[-m,m])
        tmp=-list.vrange[0]
        if (tmp gt 0) then list.vrange=list.vrange+tmp
        tmp=list.vrange[1]-(*list.tiffs)[2]+1
        if (tmp gt 0) then list.vrange=list.vrange-tmp

        xs=(list.hrange[1]-list.hrange[0]+1)
        ys=(list.vrange[1]-list.vrange[0]+1)
        list.dynh=round(xs/list.sf)
        list.dynv=round(ys/list.sf)
        ZOOM_CHANGE,ev,list,xs,ys,/USLIDER
    endelse
endif else begin ; zoom or unzoom
    list.xs=ev.x
    list.ys=ev.y
    Widget_Control, ev.id, Event_Pro='ZOOM_DRAWBOX', $
        draw_motion=1
    Widget_Control, ev.top, Set_UValue=list
endelse

end;pro ZOOM_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ZOOM_DRAWBOX ,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]

IF thisEvent EQ 'SCROLL' THEN BEGIN
    XRDUA_event, ev
    return
ENDIF
IF thisEvent EQ 'UP' THEN BEGIN
    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    Widget_Control, ev.id, draw_motion=0, $
       Event_Pro='ZOOM_PROCESS_EVENTS'
    x = [list.xs, ev.x]
    y = [list.ys, ev.y]

    if (x[0] eq x[1]) and (y[0] eq y[1]) then begin ;unzoom and move to click coord
        if list.sf lt 1 then list.sf=1.
        xs=(*list.tiffs)[1]
        ys=(*list.tiffs)[2]
        list.hrange=[0,xs-1]
        list.vrange=[0,ys-1]
        list.zoom=0
        list.dynh=round((*list.tiffs)[1]/list.sf)
        list.dynv=round((*list.tiffs)[2]/list.sf)

        xcc=StatTR(ev.x,list.sfh)
        ycc=StatTR(ev.y,list.sfv)
        xy=[RTDyn(xcc,list.sf,list.hrange[0]),RTDyn(ycc,list.sf,list.vrange[0])]
    endif else begin ; change zoom area and center
        IF list.xs GT ev.x THEN x = [ev.x, list.xs]
        IF list.ys GT ev.y THEN y = [ev.y, list.ys]
        
        list.hrange=0>StatTR(x,list.sfh)<((*list.tiffs)[1]-1)
        list.vrange=0>StatTR(y,list.sfv)<((*list.tiffs)[2]-1)
;        list.hrange=[1660,2047l]
;        list.vrange=[1660,2047l]
        
        list.zoom=1
        xs=(list.hrange[1]-list.hrange[0]+1)
        ys=(list.vrange[1]-list.vrange[0]+1)
        list.dynh=round(xs/list.sf)
        list.dynv=round(ys/list.sf)
        xy=0
    endelse

    ZOOM_CHANGE,ev,list,xs,ys,xy=xy,/USLIDER
    return
endif

WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
x = [list.xs, ev.x]
y = [list.ys, ev.y]
WSet, list.drawindex
x=0>x<(list.stath-1)
y=0>y<(list.statv-1)
;device,get_graphics=oldg,set_graphics=6
PlotS, [x[0], x[0]], y,Linestyle=2,/device
PlotS, [x[1], x[1]], y,Linestyle=2,/device
PlotS, x,[y[0], y[0]],Linestyle=2,/device
PlotS, x,[y[1], y[1]],Linestyle=2,/device
;device,set_graphics=oldg
Widget_Control, ev.top, Set_UValue=list
end;pro ZOOM_DRAWBOX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ZOOM_CHANGE,ev,list,xs,ys,xy=xy,USLIDER=USLIDER
drawdyn=list.drawdyn
drawdynindex=list.drawdynindex
pixindex2=list.pixindex2
dynh=list.dynh
dynv=list.dynv
sf=list.sf
tmp=CreateDyn(drawdyn,drawdynindex,pixindex2,dynh,dynv,list.stath,list.statv,xs,ys,$
    sf,ev.id,widget_info(ev.top,FIND_BY_UNAME='window scale'))
list.drawdyn=drawdyn
list.drawdynindex=drawdynindex
list.pixindex2=pixindex2
list.dynh=dynh
list.dynv=dynv
list.sf=sf

; Set view area to middle
if keyword_set(xy) then WIDGET_CONTROL, drawdyn, SET_DRAW_VIEW=[(xy[0]-list.stath/2.)>0,(xy[1]-list.statv/2.)>0] $
else WIDGET_CONTROL, drawdyn, SET_DRAW_VIEW=[(dynh-list.stath)/2.,(dynv-list.statv)/2.]

if widget_info(list.drawdyn,/VALID_ID) then begin
    Wset,list.drawdynindex
    dplot,list
    wset,list.pixindex2
    Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.drawdynindex]
    if keyword_set(USLIDER) then RefreshCircleSlider,ev.top,list
endif

ID1=widget_info(ev.top,FIND_BY_UNAME='coord')
if (widget_info(ID1,/VALID_ID) eq 1) then begin
    text=strarr(2)
    text[0]='X:'+stringr(list.hrange[0])+' - '+stringr(list.hrange[1])
    text[1]='Y:'+stringr(list.vrange[0])+' - '+stringr(list.vrange[1])
    widget_control,ID1,set_value=text
endif

Widget_Control, ev.top, Set_UValue=list
end;pro ZOOM_CHANGE
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUpInfo,ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=top
widget_control,top,sensitive=1
end;pro CleanUpInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro info_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_uvalue=v
case v of
0:xdisplayfile,'history'
1:OpenURL,'xrdua.html','http://xrdua.ua.ac.be'
2:xdisplayfile,'COPYING'
else:
endcase
;WIDGET_CONTROL, ev.top, /DESTROY
end ;pro info_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRDUASaveINI,list
result=UpdateINI('$XRD path:',{path:list.path,file:list.file})
result=UpdateINI('$XRD mskpath:',{path:list.mskpath,file:list.mskfile})
result=UpdateINI('$XRD pddpath:',{path:list.pddpath,file:list.pddfile})
result=UpdateINI('$XRD_CHI path:',{path:list.chipath,file:list.chifile})
result=UpdateINI('$XRD_XDI path:',{path:list.xdipath,file:list.xdifile})
result=UpdateINI('$Show:',{n:[list.smellipse,list.cross,list.ShowTie,list.ShowTilt,list.MarkBadPixels]})
result=UpdateINI('$Scale:',{n1:long([list.clipslider[3],list.clipslider[4],list.revcol,list.clipslider[1]]),n2:[1./list.sf,list.clipslider[2],list.displaypower]})
result=UpdateINI('$Pointer:',{n:[list.DelMass,list.DelWidth]})
result=UpdateINI('$CW:',{n:[list.CW1,list.CW2]})
result=UpdateINI('$Misc:',{n:[list.letter,list.OverlapCrit,list.mipol]})
result=UpdateINI('$CSlider:',{n:[list.cslideruseM,list.bPdyn]})
result=UpdateINI('$AutoCor:',{n:list.autocor})
result=UpdateINI('$XRD tiepath:',{path:list.TiePath,file:list.TieFile})
result=UpdateINI('$CrossUpdate:',{n:list.CrossUpdate})
result=UpdateINI('$AzInt:',{n1:[list.AzIntInfo.oD,list.AzIntInfo.soD,list.AzIntInfo.Fill,list.AzIntInfo.type,$
    list.AzIntInfo.bool,list.AzIntInfo.SaveRingPercentage,list.AzIntInfo.Batch,list.AzIntInfo.tiff,list.AzIntInfo.origformat],$
    n2:[list.AzIntInfo.MinPix, list.AzIntInfo.BWP,list.AzIntInfo.BWT,list.AzIntInfo.BWD,list.AzIntInfo.BWQ,$
    list.AzIntInfo.AzWidth],n3:list.AzIntInfo.outbin,n4:list.AzIntInfo.median,$
    n5:list.AzIntInfo.stD,n6:list.AzIntInfo.nsectors,n7:[list.AzIntInfo.aoD,list.AzIntInfo.loD],n8:list.AzIntInfo.percentile,$
    n9:[list.AzIntInfo.updatenorm,list.AzIntInfo.calcerror]})
result=UpdateINI('$ReadCCDInfo:',list.readinfo)
result=UpdateINI('$DplotThick:',{n:list.dplotthick})
result=UpdateINI('$RadAx:',{n1:[list.RadAxType,list.RadAxShow,list.RadAxMin,list.RadAxMax,list.RadAxAzimuth,list.RadAxBins],n2:list.RadAxFormat})

end;pro XRDUASaveINI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_XRDUA, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
wdelete,list.pixindex1
set_plot, list.Platform
;device,list.Visual
heap_free,list

if (list.listtype ne 0) then begin
    wdelete,list.pixindex2
    XRDUASaveINI,list
endif

cd,XRDUAworkingDIR()
error=UnSubscribeToOutID()

print,'=========================================STOP XRDUA'
end;pro CleanUp_XRDUA
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRDUA_event, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
0:    begin
    ; This is in developement
;    TopBase=widget_info(ev.top,find_by_uname='XRDUAtop')
;    widget_control,ev.top,get_uvalue=list
;    TLBsize=list.TLBsize
;    propagateresize,TopBase,list.TLBsize,'draw',list.statv/float(list.stath)
;    ID=widget_info(ev.top,find_by_uname='draw')
;    if widget_info(ID,/valid) then begin
;        geom=widget_info(long(ID),/GEOMETRY)
;        help,geom,/struc
;        ;list.stath=!d.x_size
;        ;list.statv=!d.y_size
;        if ptr_valid(list.tiffs) then begin
;            list.sfh=(*list.tiffs)[1]/float(list.stath)
;            list.sfv=(*list.tiffs)[2]/float(list.statv)
;        endif
;    endif
;    list.TLBsize=TLBsize
;    widget_control,ev.top,set_uvalue=list
    endcase
1 : begin ;button event
    WIDGET_CONTROL, ev.id, GET_VALUE=val
    case val of
    'Read 2D Pattern':begin
          OpenFile,ev,1,/info
          endcase
    'Write 2D Pattern':begin
          widget_control,ev.top,get_uvalue=list
          if n_elements(list) eq 0 then return

          path=SelFile(list.path,list.file,'*'+FileFormat(2),'Save CCD...')
          if path EQ '' then return

          ; Write file
          bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
          struc={center:list.center,dist:list.dist,file:list.file,scrx:list.scrx,$
          lambda:list.lambda,darkfile:(bkcode[0] eq 2)?list.darkfile:''}
          if WriteCCD(list.tiff,path,struc=struc) then printw,list.OutID,'Pattern saved: '+path
          endcase
    'Auto Load':begin
          WIDGET_CONTROL, ev.top, get_uvalue=list
          if n_elements(list) eq 0 then return
          ID1=ev.top
          filter='*'+list.Format
          list2={top:ID1,update:1B,timer:0B,PTR_files:PTR_NEW(),path:list.path,sec:10.,filter:filter,$
              ev:ev,chari:0,char:['  \  ','  |  ','  /  '],title:'Find new '+list.Format+' files',$
              sortseparator:list.sortseparator,sortmethod:list.sortmethod,$
              filtermethod:-1,binit:1b}
          base=widget_base(/column,title=list2.title+(list2.char)[list2.chari],xoffset=200,yoffset=200,$
          uvalue=list2)

          droplist=widget_list(base,uname='FilesList',ysize=20)
          text=widget_text(base,value='0 files',uname='FilesText')

          base1=widget_base(base,/row,/nonexclusive)
              button1=widget_button(base1,value='Update Main Window')
              button2=widget_button(base1,value='Run timer')

          base2=widget_base(base,/row)
              lbl=widget_label(base2,value='Refresh rate(sec):')
              txt=widget_text(base2,value=stringr(list2.sec),/editable)
          releasebutton=widget_button(base,value='Stop',uname='StopButton')

          WIDGET_CONTROL, base, /REALIZE
          widget_control,button1,set_button=list2.update
          widget_control,button2,set_button=list2.timer
          WIDGET_CONTROL, base, TIMER=0
          Xmanager,'XRDUA_event',base, event_handler='AutoLoad_event',$
            /no_block,cleanup='CleanUp_AutoLoad',GROUP_LEADER=ID1

          endcase
    'Exit':  begin
          WIDGET_CONTROL, ev.top, /DESTROY
          endcase
    'Online Support':begin
          OpenURL,'xrdua.html','http://xrdua.ua.ac.be'
          endcase
       'Online Update':begin
             widget_control,ev.top,get_uvalue=list
             ;if list.Platform eq 'X' then add='L' else add='W'
             add='B'
             URL='http://xrdua.ua.ac.be/index.php?update='+add+'.'+Version()
             OpenURL,'update.html',URL
             endcase
    'About XRDUA':begin
          widget_control,ev.top,sensitive=0
          widget_control,ev.top,get_uvalue=list

          base=widget_base(title='About XRDUA v'+Version()+' ('+((!version.memory_bits eq 64)?'64':'32')+'bit)',/column,uvalue=ev.top)
          text=['Copyright (C) 2005-'+(strsplit(SYSTIME(),' ',/extract))[4]+ ' Wout De Nolf',$
            '',$
            'University of Antwerp - Dep. Chemistry',$
            'Universiteitsplein 1',$
            'B-2610 Wilrijk BELGIUM',$
            'E-mail: woutdenolf@users.sf.net',$
            '',$
            'XRDUA is free software; you can redistribute it and/or',$
            'modify it under the terms of the GNU General Public License',$
            'as published by the Free Software Foundation; either version 3',$
            'of the License, or (at your option) any later version.',$
            '',$
            'XRDUA is distributed in the hope that it will be useful,',$
            'but WITHOUT ANY WARRANTY; without even the implied warranty of',$
            'MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the',$
            'GNU General Public License for more details.',$
            '',$
            'You should have received a copy of the GNU General Public License',$
            'along with XRDUA.  If not, see <http://www.gnu.org/licenses/>.']

          bdraw=QUERY_BMP(ProgramRootDir()+'ua.bmp')
          if bdraw then begin
              logo=READ_IMAGE(ProgramRootDir()+'ua.bmp')
              s=size(logo)
              draw=widget_draw(base,xsize=s[2],ysize=s[3],uvalue=-1)
              ;button3=widget_button(base2,value=ProgramRootDir()+'ua.bmp',/bitmap,uvalue=-1)
          endif

          reftext=['Please refer to: DOI 10.1002/sia.3125',$
                      '    W. De Nolf and K. Janssens. Micro X-ray diffraction and fluorescence ',$
                      '    tomography for the study of multilayered automotive paints.',$
                      '    Surf. Interface Anal., 2010, 42, 411-418']
          ref=widget_text(base,value=reftext,ysize=n_elements(reftext),xsize=70)
          
          text=widget_text(base,ysize=n_elements(text),xsize=70,value=text)
          button=widget_button(base,value='History...',uvalue=0)
          button=widget_button(base,value='GNU General Public License ...',uvalue=2)
          button=widget_button(base,value='Home page...',uvalue=1)

          widget_control,base,/realize
          if bdraw then tvscl,logo,true=1
          Xmanager,'XRDUA_event',base, event_handler='info_event',GROUP_LEADER=ev.top,cleanup='CleanUpInfo'
          endcase
    'Manual'  :    begin
          CATCH, Error_status ; Only in IDL6.0
          IF Error_status EQ 0 THEN ONLINE_HELP, BOOK='XRDUAManual.pdf'
          endcase
     'Configuration File':begin
          inifile=ReadINI('',/EDIT)
          xdisplayfile,inifile, /EDITABLE
          endcase
 'Save Configuration':begin
           widget_control,ev.top,get_uvalue=list
           XRDUASaveINI,list
           endcase
'Register DLMs': RegisterDLMs,ev.top
 'Zoom'  :begin
          mode,ev,val
          endcase
 'Zingers'  :begin
          mode,ev,val
          endcase
 'Saturation'  :begin
          mode,ev,val
          endcase
 'Experimental Geometry':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.top,sensitive=0
          base=widget_base(title=val,uvalue=ev.top,/row)
          base1=widget_base(base,/column)
          xs=150
              baset=widget_base(base1,/row)
              label=widget_label(baset,value='X-ray energy (keV)  :        ',xsize=xs)
              text=widget_text(baset,/editable,value= $
              stringr(list.cconv/list.lambda),uvalue='xen',uname='xen')

              baset=widget_base(base1,/row)
              label=widget_label(baset,value='X-ray wavelength ('+msymbols('angstroms')+')  :      ',xsize=xs)
              text=widget_text(baset,/editable,value= $
              stringr(list.lambda),uvalue='xlambda',uname='xlambda')

              baset=widget_base(base1,/row)
              label=widget_label(baset,value='Pixel X-size ('+msymbols('micro')+'m) :          ',xsize=xs)
              text=widget_text(baset,/editable,value= $
              stringr(list.scrx*10^6.),uvalue='scrx',uname='scrx')

              baset=widget_base(base1,/row)
              label=widget_label(baset,value='Pixel Y-size ('+msymbols('micro')+'m) :          ',xsize=xs)
              text=widget_text(baset,/editable,value= $
              stringr(list.scry*10^6.),uvalue='scry',uname='scry')

              baset=widget_base(base1,/row)
              label=widget_label(baset,value='Distance Sample-Detector(cm):',xsize=xs)
              text=widget_text(baset,/editable,value= $
              stringr(list.dist*10^2.),uvalue='dist',uname='dist')

              baset=widget_base(base1,/row)
              label=widget_label(baset,value='CenterX(pixels):',xsize=xs)
              text=widget_text(baset,value=stringr(list.center[0]),/editable,uvalue='xcen',uname='xcen')

              baset=widget_base(base1,/row)
              label=widget_label(baset,value='CenterY(pixels):',xsize=xs)
              text=widget_text(baset,value=stringr(list.center[1]),/editable,uvalue='ycen',uname='ycen')

              baset=widget_base(base1,/row)
              label=widget_label(baset,value='Tilt a('+msymbols('degrees')+'):',xsize=xs)
              text=widget_text(baset,/editable,value= $
              stringr(list.a/!dpi*180.),uvalue='a',uname='a')

              baset=widget_base(base1,/row)
              label=widget_label(baset,value='Rot b('+msymbols('degrees')+'):',xsize=xs)
              text=widget_text(baset,/editable,value=$
              stringr(list.b/!dpi*180.),uvalue='b',uname='b')

          button=widget_button(base,value='OK',uvalue='OK')
          WIDGET_CONTROL, base, /REALIZE
          Xmanager,'XRDUA_event',base, event_handler='Expd_event',cleanup='CleanUp_Expd',GROUP_LEADER=ev.top
          endcase
    'Beam position':  begin
          mode,ev,val
          endcase
    'Use Mouse':begin
          mode,ev,val
          endcase
    '3 Points':   begin
          mode,ev,val
          endcase
    'Fit Center':       begin
          mode,ev,val
          endcase
    'Surface': begin
          widget_control,ev.top,get_uvalue=list
          if list.zoom then begin
            xr=list.hrange
            yr=list.vrange
            img=(*list.tiff)[xr[0]:xr[1],yr[0]:yr[1]]
          endif else begin
              xr=[0,round((*list.tiffs)[1])-1]
            yr=[0,round((*list.tiffs)[2])-1]
            img=*list.tiff
          endelse

          list={tiff:img,zr:[min(img),max(img)],$
              xr:xr,yr:yr,proplot:'plotsurf'}
          flplot_obj,list,GROUP_LEADER=ev.top,xtitle='X(Column)',ytitle='Y(Row)',ztitle='Counts',$
          /xtrue,/ytrue
          endcase
    'Calc':       begin
          widget_control,ev.top,get_uvalue=list
          center=calc_middle(list.P1,list.P2,list.P3)
          list.center=center
          str=stringr(center)
          ID=widget_info(ev.top,FIND_BY_UNAME='center')
          widget_control,ID,set_value=str[0]+','+str[1]
          RefreshDisplay, ev.top,list
          endcase
    'Strip':   begin
          widget_control,/hourglass
          widget_control,ev.top,get_uvalue=list
          T1=systime(1)
          ;----Smooth image----
          imgb=Strip2D(list.tiff,list.wb,list.niter)
          T3=systime(1)
          printw,list.OutID,string(list.wb,format='("Width of filter:",I)')
          printw,list.OutID,string(list.niter,format='("Iterations performed:",I)')
          printw,list.OutID,'Stripping: '+string(T3-T1)+' sec'
          ptr_free,list.tiffb
          list.tiffb=imgb
          widget_control,ev.top,set_uvalue=list
          endcase
    'Subtract backgr.':begin
          widget_control,ev.top,get_uvalue=list

          T=systime(1)
          bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])

          done=0b
          if ptr_valid(list.tiffb) then begin
              BackSub,list.tiff,list.tiffb,bkcode[0]
              done=1b
          endif
          if bkcode[1] then begin
              ind=where((*list.tiff) le list.thres[0],ct)
              if ct ne 0 then (*list.tiff)[ind]=0
              done=1b
          endif
          if bkcode[2] then begin
              ind=where((*list.tiff) ge list.thres[1],ct)
              if ct ne 0 then (*list.tiff)[ind]=0
              done=1b
          endif
          if not done then return
          printw,list.OutID,'Background subtraction: time (sec):'+string(systime(1)-T)

          list.CorDone[((*(list.CorSeq))).back]=1
          widget_control,ev.top,set_uvalue=list
          mode,ev,val
          widget_control,ev.top,get_uvalue=list
          RefreshDisplay,ev.top,list
          endcase
    'Reset backgr.': begin
          RestoreCorrectionev,ev
          endcase
    'Show backgr.':     begin
          widget_control,ev.top,get_uvalue=list
          if not ptr_valid(list.tiffb) then return
          x1=list.hrange[0]
          x2=list.hrange[1]
          y1=list.vrange[0]
          y2=list.vrange[1]
          img=(*list.tiff)[x1:x2,y1:y2]
          imgb=(*list.tiffb)[x1:x2,y1:y2]
          s=size(img)
          x=indgen(s[1])+x1
          y=indgen(s[2])+y1
          xr=[x1,x2]
          yr=[y1,y2]
          list={surf:img,imgfit:imgb,x:x,y:y,xr:xr,yr:yr,zr:[min(img),max(img)],proplot:'backplot'}
          LoadctRev,1,/silent
          flplot_obj,list,/xtrue,/ytrue
          endcase
    'Save backgr.': begin
          widget_control,ev.top,get_uvalue=list
          if not ptr_valid(list.tiffb) then return

          path=SelFile(list.path,list.file,'*'+FileFormat(2),'Save Background....')
          if path EQ '' then return

          ; Write file
          bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
          struc={center:list.center,dist:list.dist,file:list.file,scrx:list.scrx,$
          lambda:list.lambda,darkfile:(bkcode[0] eq 2)?list.darkfile:''}
          error=WriteCCD(list.tiffb,path,struc=struc)
          printw,list.OutID,'Background saved: '+path
          endcase
    'modify...':     begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_uvalue=uvalue
          ID1=widget_info(ev.top,FIND_BY_UNAME='bP1')
          ID2=widget_info(ev.top,FIND_BY_UNAME='bP2')
          ID3=widget_info(ev.top,FIND_BY_UNAME='bP3')
          case uvalue of
          'bP1':begin
              widget_control,ID1,sensitive=0
              widget_control,ID2,sensitive=1
              widget_control,ID3,sensitive=1
              list.bP=1
              endcase
          'bP2':begin
              widget_control,ID1,sensitive=1
              widget_control,ID2,sensitive=0
              widget_control,ID3,sensitive=1
              list.bP=2
              endcase
          'bP3':begin
              widget_control,ID1,sensitive=1
              widget_control,ID2,sensitive=1
              widget_control,ID3,sensitive=0
              list.bP=3
              endcase
          endcase
          widget_control,ev.top,set_uvalue=list
          endcase
 '  Show Points':begin
          widget_control,ev.top,get_uvalue=list
          list.cross=1
          widget_control,ev.id,set_value='# Show Points'
          RefreshDisplay, ev.top,list
          endcase
'  Show Debye Marker':begin
          widget_control,ev.top,get_uvalue=list
          list.smellipse=list.smellipse or 1
          widget_control,ev.id,set_value='# Show Debye Marker'
          RefreshDisplay, ev.top,list
          endcase
'  Show Debye Marker Axis':begin
          widget_control,ev.top,get_uvalue=list
          list.smellipse=list.smellipse or 2
          widget_control,ev.id,set_value='# Show Debye Marker Axis'
          RefreshDisplay, ev.top,list
          endcase
    'View Mask':begin
          EditMaskWindow,ev,0,0
          endcase
    'Edit Mask':begin
          EditMaskWindow,ev,1,0
          endcase
'# Show Points':begin
          widget_control,ev.top,get_uvalue=list
          list.cross=0
          widget_control,ev.id,set_value='  Show Points'
          RefreshDisplay, ev.top,list
          endcase
'# Show Debye Marker':begin
          widget_control,ev.top,get_uvalue=list
          list.smellipse=list.smellipse and 2
          widget_control,ev.id,set_value='  Show Debye Marker'
          RefreshDisplay, ev.top,list
          endcase
'# Show Debye Marker Axis':begin
          widget_control,ev.top,get_uvalue=list
          list.smellipse=list.smellipse and 1
          widget_control,ev.id,set_value='  Show Debye Marker Axis'
          RefreshDisplay, ev.top,list
          endcase
'# Show SpaDist':   begin
          widget_control,ev.top,get_uvalue=list
          list.ShowTie=bytarr(4)
          widget_control,ev.id,set_value='  Show SpaDist'
          RefreshDisplay, ev.top,list
          endcase
'  Show SpaDist':   begin
          widget_control,ev.top,get_uvalue=list
          list.ShowTie=[1b,0b,0b,0b]
          widget_control,ev.id,set_value='# Show SpaDist'
          RefreshDisplay, ev.top,list
          endcase
'# Show Calib':  begin
          widget_control,ev.top,get_uvalue=list
          list.ShowTilt=0
          widget_control,ev.id,set_value='  Show Calib'
          RefreshDisplay, ev.top,list
          endcase
'  Show Calib':  begin
          widget_control,ev.top,get_uvalue=list
          list.ShowTilt=1
          widget_control,ev.id,set_value='# Show Calib'
          RefreshDisplay, ev.top,list
          endcase
'# Show bad pixels':begin
          widget_control,ev.top,get_uvalue=list
          list.MarkBadPixels=0
          widget_control,ev.id,set_value='  Show bad pixels'
          RefreshDisplay, ev.top,list
          endcase
'  Show bad pixels':begin
          widget_control,ev.top,get_uvalue=list
          list.MarkBadPixels=1
          widget_control,ev.id,set_value='# Show bad pixels'
          RefreshDisplay, ev.top,list
          endcase
'# Show DynMark':begin
          widget_control,ev.top,get_uvalue=list
          list.bPdyn=0
          widget_control,ev.id,set_value='  Show DynMark'
          RefreshDisplay, ev.top,list
          endcase
'  Show DynMark':begin
          widget_control,ev.top,get_uvalue=list
          list.bPdyn=1
          widget_control,ev.id,set_value='# Show DynMark'
          RefreshDisplay, ev.top,list
          endcase
'# Show Radial Axis':begin
          widget_control,ev.top,get_uvalue=list
          list.RadAxShow=0
          widget_control,ev.id,set_value='  Show Radial Axis'
          RefreshDisplay, ev.top,list
          endcase
'  Show Radial Axis':begin
          widget_control,ev.top,get_uvalue=list
          list.RadAxShow=1
          widget_control,ev.id,set_value='# Show Radial Axis'
          RefreshDisplay, ev.top,list
          endcase
'View PDF':   begin
          EditPDFWindow,ev,'RefreshDisplay',/topid
          endcase
    'Mask Setup':begin
          mode,ev,val
          endcase
    'Arc':     begin
          mode,ev,val
          endcase
    'Ellipse':  begin
          mode,ev,val
          endcase
    'Ellipse (0,0)':begin
          mode,ev,val
          endcase
    'Rectangle':  begin
          mode,ev,val
          endcase
    'Range':begin
          AreaRangeWindow,ev.top
          endcase
    'Full Arc':begin
          Widget_Control, ev.top, get_UValue=list
          
          add=GetFullArc(list)
          
          s=(*list.tiffs)
          if list.masknr eq 0 then begin ; New
            list.mask=ptr_new(add)
            list.smask=ptr_new([1b])
            list.masknr++
          endif else begin
            ind=where((*list.mask)[2,*] eq 0,count)

            if count eq 0 then begin ; Add
                L=max((*list.mask)[1,*])+1
                add[1]=L
                (*list.mask)=[[(*list.mask)],[add]]
                (*list.smask)=[(*list.smask),1b]
                list.masknr++
            endif else begin ; Replace last
                L=(*list.mask)[1,ind[0]]
                add[1]=L
                (*list.mask)[*,ind[0]]=add
            endelse
          endelse

          RefreshDisplay, ev.top,list
          endcase
    'Whole Image':begin
          Widget_Control, ev.top, get_UValue=list
          s=(*list.tiffs)
          if list.masknr eq 0 then begin ; New
            list.mask=ptr_new([1,0,0,0,s[1]-1,0,s[2]-1])
            list.smask=ptr_new([1b])
            list.masknr++
          endif else begin
            ind=where((*list.mask)[2,*] eq 0,count)

            if count eq 0 then begin ; Add
                L=max((*list.mask)[1,*])+1
                (*list.mask)=[[(*list.mask)],[1,L,0,0,s[1]-1,0,s[2]-1]]
                (*list.smask)=[(*list.smask),1b]
                list.masknr++
            endif else begin ; Replace last
                L=(*list.mask)[1,ind[0]]
                (*list.mask)[*,ind[0]]=[1,L,0,0,s[1]-1,0,s[2]-1]
            endelse
          endelse

          RefreshDisplay, ev.top,list
          endcase
    'Set':     begin
          widget_control,ev.top,get_uvalue=list
          if list.masknr ne 0 then begin
          (*list.mask)[2,0:list.masknr-1]=1
          RefreshDisplay, ev.top,list
          endif
          endcase
    'Save Mask':begin
          widget_control,ev.top,get_uvalue=list
          if n_elements(list) eq 0 then return
          SaveMask,list,1
          widget_control,ev.top,set_uvalue=list
          endcase
    'Load Mask':begin
          widget_control,ev.top,get_uvalue=list
          if n_elements(list) eq 0 then return
          LoadMask,list,1,list.OutID
          RefreshDisplay, ev.top,list
          mode,ev,list.mode,/refresh
          endcase
    'Add Mask':begin
          widget_control,ev.top,get_uvalue=list
          if n_elements(list) eq 0 then return
          LoadMask,list,1,list.OutID,/flags,flagi=['$CENTER:','$SOURCE ('+msymbols('angstroms')+'):','$PIXSIZE (m):',$
              '$DET-SAM (m):','$BACK:','$SEARCH:','$TILTA:','$AZIOFF:','$IntCor2Dinfo:','$CORRECTIONS:','$FF:','$MASKOFF:',$
              '$SATURATION:','$ZINGERS:','$AzInt:','$ReadCCDInfo:']
          RefreshDisplay, ev.top,list
          endcase
    'Azimuthal Integration': az_int,ev
    'Automatic Integration': az_int_auto,ev
    'Automatic Correction':    begin
                            widget_control,ev.top,get_uvalue=list
                            AutoCorrect,list
                            RefreshDisplay,ev.top,list
                            mode,{top:ev.top},list.mode,/refresh
                            endcase
    'Auto Set with PDF': CalibFromPDF,ev,'Tilt'
    'Auto Set with PDF ': CalibFromPDF,ev,'Tie'
    'Auto Set': auto_set,ev,'Tilt'
    'Auto Set ': auto_set,ev,'Tie'
'Surface Integration': begin
          widget_control,ev.top,get_uvalue=list
          if (list.masknr eq 0) then return
          surf_int,list,0,ev.top,'Select area to integrate'
          endcase
'Rectangle summation': begin
          widget_control,ev.top,get_uvalue=list
          if (list.masknr eq 0) then return
          rectangle_sum,list,ev.top,'Select area to summate'
          endcase
    'Background':begin
          mode,ev,val
          endcase
    'Debye Marker':    begin
          mode,ev,val
          endcase
    'Spatial Distortion':begin
          mode,ev,val
          endcase
    'Horizontal streaking':begin
          widget_control,ev.top,get_uvalue=list
          list.SaturHor=ev.select
          widget_control,ev.top,set_uvalue=list
          endcase
    'Vertical streaking':begin
          widget_control,ev.top,get_uvalue=list
          list.SaturHor=~ev.select
          widget_control,ev.top,set_uvalue=list
          endcase
    'Remove Saturation':begin
          widget_control,ev.top,get_uvalue=list
;          if CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'satur') eq 0 then begin
;              printw,list.OutID,"Can't perfom Saturation removal"
;              return
;          endif

          widget_control, /hourglass
          T=systime(1)
          tmp=list.tiffvalid
          RemoveSatur,list.tiff,tmp,list.SaturMax,list.SaturAspect,list.SaturNMax,count=count
          list.tiffvalid=tmp
          printw,list.OutID,'Saturated pixels found: '+stringr(count)+' in '+string(systime(1)-T)+' sec'

          list.CorDone[((*(list.CorSeq))).satur]=1
          RefreshDisplay, ev.top,list
          mode,ev,val
          endcase
'Restore Saturation':begin
          RestoreCorrectionev,ev
          endcase
'Restore Zingers':begin
          RestoreCorrectionev,ev
          endcase
'Remove Zingers':begin
          widget_control,ev.top,get_uvalue=list
;          if CorrAllow(list.CorDone,list.SeqTags,list.CorSeq,list.CorMethod,'zinger') eq 0 then begin
;              printw,list.OutID,"Can't perfom zinger removal"
;              return
;          endif

          widget_control, /hourglass
          T=systime(1)
          tmp=list.tiffvalid
          RemoveZingers,list.tiff,tmp,width=list.ZingerWidth,threshold=list.ZingerThreshold,count=count
          list.tiffvalid=tmp
          printw,list.OutID,'Zingers found: '+stringr(count)+' in '+string(systime(1)-T)+' sec'

          list.CorDone[((*(list.CorSeq))).zinger]=1
          RefreshDisplay, ev.top,list
          mode,ev,val
          endcase
'Next file':begin
          widget_control,ev.top,get_uvalue=list
          if (Systime(1)-list.QueTimeStamp) le list.QueTres then return
          path = file_search(list.path+'*'+list.Format)
          nfiles = n_elements(path)
          path = SortListFiles(path,list.sortmethod,list.sortseparator)
          ind = where(path eq list.path+list.file,count)
          if (count eq 0) or (ind[0] eq nfiles-1) then ind=0 else ind=ind+1
          path=path[ind]
          path=path[0]
          error=CutPath(path,file=file,ext=ext)
          list.file=file+ext
          widget_control,ev.top,set_uvalue=list
          OpenFile,ev,0,/info
          endcase
'Previous file':begin
          widget_control,ev.top,get_uvalue=list
          if (Systime(1)-list.QueTimeStamp) le list.QueTres then return
          path = file_search(list.path+'*'+list.Format)
          nfiles = n_elements(path)
          path = SortListFiles(path,list.sortmethod,list.sortseparator)
          ind = where(path eq list.path+list.file,count)
          if (count eq 0) or (ind[0] eq 0) then ind=nfiles-1 else ind=ind-1
          path=path[ind]
          path=path[0]
          error=CutPath(path,file=file,ext=ext)
          list.file=file+ext
          widget_control,ev.top,set_uvalue=list
          OpenFile,ev,0,/info
          endcase
'Close Ring':   begin
          widget_control,ev.top,get_uvalue=list
          if (*list.nrTie)[list.nrTieRing] eq 0 then return
          ; ----Point to next ring----
          list.nrTieRing++
          (*list.nrTie)=[(*list.nrTie),0]
          printw,list.OutID,'Pointer to ring:'+string(list.nrTieRing+1)
          RefreshDisplay, ev.top,list
          endcase
'Delete small':begin
          widget_control,ev.top,get_uvalue=list
          if list.nrTiltRing eq 0 then return
          perc=float(promptnumber('1',ev.top,'Delete the smallest tie points (threshold = x% of maximum):'))/100.
          if perc le 0 or perc ge 100 then return
          
          ; ----Get ring info----
          x=(*list.TiltI)[0,*]
          y=(*list.TiltI)[1,*]
          z=(*list.tiff)[x,y]
          n=list.nrTiltRing
          ni=*list.nrTilt
          dspac=*list.TiltRingD
          
          ; ----Delete all rings----
          list.nrTiltRing=0
          (*list.nrTilt)=[0]
          ptr_free,list.TiltI
          ptr_free,list.TiltRingD
          
          ; ----Remove small tie points----
          ma=max(z)
          b=bytarr(n_elements(x))
          off=0ul
          for i=0l,n-1 do begin
              add=ni[i]
              ;ma=max(z[off:off+add-1])
              b[off:off+add-1]=z[off:off+add-1] gt perc*ma
              ni[i]=total(b[off:off+add-1])
              off+=add
          endfor
          
          ; ----Keep non empty rings----
          ind=where(ni ne 0,n)
          if n ne 0 then begin
              last=n_elements(ni)-1
              if ni[last] eq 0 then begin
                  ind=[ind,last]
                  n++
              endif
              list.nrTiltRing=n-1
              (*list.nrTilt)=ni[ind]
              list.TiltRingD=ptr_new(dspac[ind])
              ind=where(b,ct)
              ind=reform(ind,1,ct,/overwrite)
              list.TiltI=ptr_new([x[ind],y[ind]])
          endif
          
          RefreshDisplay, ev.top,list
          
          endcase
'Close TiltRing':begin
          widget_control,ev.top,get_uvalue=list
          if (*list.nrTilt)[list.nrTiltRing] eq 0 then return
          ; ----Sort points by azimuth----
          if list.nrTiltRing ne 0 then n=total((*list.nrTilt)[0:list.nrTiltRing-1]) $
          else n=0
          nj=(*list.nrTilt)[list.nrTiltRing]
          points=(*list.TiltI)[*,n:n+nj-1]
          SortPoints,points,nj,list.center
          (*list.TiltI)[*,n:n+nj-1]=points
          ; ----Point to next ring----
          list.nrTiltRing++
          (*list.nrTilt)=[(*list.nrTilt),0]
          printw,list.OutID,'Pointer to ring:'+string(list.nrTiltRing+1)
          RefreshDisplay, ev.top,list
          endcase
'Delete Tie Point':begin
          mode,ev,val
          endcase
'Add Tie Point':begin
          mode,ev,val
          endcase
'Delete Point':begin
          mode,ev,val
          endcase
'Select area':begin
          mode,ev,val
          endcase
'Delete TiltRing':begin
          widget_control,ev.top,get_uvalue=list
          if (list.nrTiltRing eq 0) then return
          ID1=ev.top
          WIDGET_CONTROL, ID1, sensitive=0

          base=widget_base(/column,title='Select ring to delete',xoffset=200,yoffset=200,$
          uvalue={list:list,ID1:ID1,select:bytarr(list.nrTiltRing+1),Call:'Tilt'})
          ; ----Ring buttons----
          nr=5
          nbase=ceil(float(list.nrTiltRing+1)/nr)
          base6=lonarr(nbase)
          for i=0l,nbase-1 do base6[i]=widget_base(base,/row,/nonexclusive)
          rbutton = lonarr(list.nrTiltRing+2)
          j=0
          for i=0l,list.nrTiltRing do begin
              val=stringr(i)
              rbutton[i]=widget_button(base6[j],value='ring '+val,uname=val)
              if ((i+1) mod nr) eq 0 then j=j+1
          endfor
          rbutton[i]=widget_button(base6[nbase-1],value='Select All')
          releasebutton=widget_button(base,value='OK',uname='OK')

          WIDGET_CONTROL, base, /REALIZE
          for i=0l,list.nrTiltRing+1 do widget_control,rbutton[i],set_button=0
          Xmanager,'XRDUA_event',base, event_handler='del_ring_event',$
          cleanup='CleanUp_del_ring',GROUP_LEADER=ev.top
          endcase
'Delete Ring':begin
          widget_control,ev.top,get_uvalue=list
          if (list.nrTieRing eq 0) then return
          ID1=ev.top
          WIDGET_CONTROL, ID1, sensitive=0
          base=widget_base(/column,title='Select ring to delete',xoffset=200,yoffset=200,$
          uvalue={list:list,ID1:ID1,select:bytarr(list.nrTieRing+1),Call:'Tie'})
          ; ----Ring buttons----
          nr=5
          nbase=ceil(float(list.nrTieRing+1)/nr)
          base6=lonarr(nbase)
          for i=0l,nbase-1 do base6[i]=widget_base(base,/row,/nonexclusive)
          rbutton = lonarr(list.nrTieRing+2)
          j=0
          for i=0l,list.nrTieRing do begin
              val=stringr(i)
              rbutton[i]=widget_button(base6[j],value='ring '+val,uname=val)
              if ((i+1) mod nr) eq 0 then j=j+1
          endfor
          rbutton[i]=widget_button(base6[nbase-1],value='Select All')
          releasebutton=widget_button(base,value='OK',uname='OK')

          WIDGET_CONTROL, base, /REALIZE
          for i=0l,list.nrTieRing+1 do widget_control,rbutton[i],set_button=0
          Xmanager,'XRDUA_event',base, event_handler='del_ring_event',$
          cleanup='CleanUp_del_ring',GROUP_LEADER=ev.top
          endcase
'Set Point':    begin
          mode,ev,val
          endcase
'Azimuthal offset':    begin
          mode,ev,val
          endcase
'Load spline file':begin
          widget_control,ev.top,get_uvalue=list
          stopl=0

          path=''
          repeat begin
              path=DIALOG_PICKFILE(path=list.TiePath,file=list.TieFile,filter='*.*',title='Load Spline File...')
              if path EQ '' then return
              error=CutPath(path,path=TiePath,file=TieFile,ext=ext)
              list.TiePath=TiePath
              list.TieFile=TieFile+ext
            struc=ReadFit2DSpline(path,error=error)
            if error then begin
                printw,list.OutID,'Error on loading '+path
            endif else begin
                list.GridDim=struc.GridDim
                list.GridOffset=struc.GridOffset
                list.GridOrigin=struc.GridOrigin
                list.GridTilt=struc.GridTilt
                list.GridSpReal=struc.GridSpReal
                list.GridSpPix=struc.GridSpPix
                list.scrx=struc.scrx
                list.scry=struc.scry
                PTR_FREE,list.tck1
                PTR_FREE,list.tck2
                list.tck1=PTR_NEW(struc.tck1)
                list.tck2=PTR_NEW(struc.tck2)
                ptr_free,list.SpaDisStruc
                stopl=1
            endelse
          endrep until stopl

          widget_control,ev.top,set_uvalue=list
          endcase
'Correct Spatial':begin
          widget_control,ev.top,get_uvalue=list
                  T=systime(1)
                  
                  ptr_free,list.SpaDisStruc

                if not PTR_VALID(list.SpaDisStruc) then begin
                    if list.Tiegrid eq 2 then begin
                        if not PTR_VALID(list.tck1) then begin
                            struc=ReadFit2DSpline(list.TiePath+list.TieFile,error=error)
                            if error then begin
                                printw,list.OutID,'Error on loading '+list.TiePath+list.TieFile
                            endif else begin
                                list.GridDim=struc.GridDim
                                list.GridOffset=struc.GridOffset
                                list.GridOrigin=struc.GridOrigin
                                list.GridTilt=struc.GridTilt
                                list.GridSpReal=struc.GridSpReal
                                list.GridSpPix=struc.GridSpPix
                                list.scrx=struc.scrx
                                list.scry=struc.scry
                                PTR_FREE,list.tck1
                                PTR_FREE,list.tck2
                                list.tck1=PTR_NEW(struc.tck1)
                                list.tck2=PTR_NEW(struc.tck2)
                            endelse
                          endif
                          if PTR_VALID(list.tck1) then $
                              list.SpaDisStruc=WarpSpaDistPrepSpline(list.tiffs,*list.tck1,*list.tck2,vflip=list.vFlipSpline)
                    endif else begin
                        if (list.nrTieRing eq 0) and ((*list.nrTie)[0] eq 0) then begin
                            printw,list.OutID,'No tie points defined'
                        endif else begin
                            if list.Tiegrid then TieO=list.TieO $
                            else TieO=CalcTieO(list.nrTieRing,list.nrTie,list.TieI,list.center,list.lambda,list.dist,list.scrx,list.scry,list.TieRingD,list.a,list.b)
                            list.SpaDisStruc=WarpSpaDistPrep(list.tiffs,list.TieI,TieO,fast=list.bFastSpline)
                            if list.Tiegrid eq 0 then ptr_free,TieO
                        endelse
                    endelse
                endif
                
                if not PTR_VALID(list.SpaDisStruc) then return
                tmp=list.tiffvalid
                if list.Tiegrid eq 2 then WarpSpaDistRunSpline,list.tiff,list.SpaDisStruc,tiffvalid=tmp,fast=list.bFastSpline $
                else WarpSpaDistRun,list.tiff,list.SpaDisStruc,tiffvalid=tmp
                list.tiffvalid=tmp

                t2=systime(1)
                printw,list.OutID,'Spatial distortion: time (sec):'+string(t2-T)

          list.CorDone[((*(list.CorSeq))).spatial]=1
          RefreshDisplay, ev.top,list
          mode,ev,val
          endcase
'Calibrate':   begin
          widget_control,ev.top,get_uvalue=list
          if (list.nrTiltRing eq 0) then return
          ID1=ev.top
          WIDGET_CONTROL, ID1, sensitive=0

          nref=11
          refbool=bytarr(nref)
          ind=where((*list.TiltRingD)[0:list.nrTiltRing-1] ne 0,count)
          refbool[[0,1,2,3,5,8,11]]=1
          if (count eq 0) then begin
              ref3=0
              refbool[6]=1
          endif else begin
              ref3=1
          endelse


          tlist={a:list.a,b:list.b,phipixphi0:list.phipixphi0,center:list.center,lambda:list.lambda,dist:list.dist,TiltRingD:list.TiltRingD,scrx:list.scrx,$
                  scry:list.scry,top:ID1,k:4.,method:1,ref3:ref3,$
                  Nrings:(bytarr(list.nrTiltRing)+1),refbool:refbool,outid:0L,Index:replicate(-1L,4)}
          topbase=widget_base(/row,title='Calibration Options',xoffset=200,yoffset=200,uvalue=tlist)

          base=widget_base(topbase,/column)

          ; ----Ring buttons----
          label=widget_label(base,value='Debye rings:',/align_left)
          nr=5
          nbase=ceil(float(list.nrTiltRing+1)/nr)
          base6=lonarr(nbase)
          for i=0l,nbase-1 do base6[i]=widget_base(base,/row,/nonexclusive)
          rbutton = lonarr(list.nrTiltRing+1)
          j=0
          for i=0l,list.nrTiltRing-1 do begin
              val=stringr(i)
              rbutton[i]=widget_button(base6[j],value=val,uname=val,uvalue=0B)
              if ((i+1) mod nr) eq 0 then j=j+1
          endfor
          rbutton[i]=widget_button(base6[nbase-1],value='Select All')

          tab=widget_tab(base)
          Estbase=widget_base(tab,/column,title='Initial estimation')
          Wrapbase=widget_base(tab,/row,title='Advanced')
              base1=widget_base(Wrapbase,/column,uname='base1')
              base11=widget_base(base1,/row)
              base110=widget_base(base11,/column,/nonexclusive)
              base111=widget_base(base11,/column,/nonexclusive)

                refbutton=lonarr(nref)
                refbutton[0]=widget_button(base110,value='Refine tilt angle')
                refbutton[1]=widget_button(base110,value='Refine rotation angle')
              refbutton[2]=widget_button(base110,value='Refine center')
              refbutton[3]=widget_button(base110,value='Refine Distance')
              refbutton[4]=widget_button(base110,value='Refine two-theta',uname='Ref1')
              refbutton[5]=widget_button(base110,value='Refine Wavelength',$
                              uname='Ref2',sensitive=refbool[4] eq 0)
              refbutton[6]=widget_button(base110,value='Refine D-spacings',$
                              uname='Ref3',sensitive=(refbool[4] eq 0) and ref3)
              refbutton[7]=widget_button(base111,value='Use intensities as weights')
              refbutton[8]=widget_button(base111,value='Reject outliers')
              refbutton[9]=widget_button(base111,value='Fit with constraints')
              refbutton[10]=widget_button(base111,value='Visualize fit and outliers')

          base2=widget_base(Estbase,/row)
              text=widget_text(base2,/editable,value=stringr(tlist.k),uname='k',uvalue='k')
              button=widget_button(base2,value='Default nr of SD`s')

          base3=widget_base( Estbase,/row)
              text=widget_text(base3,/editable,value=stringr(tlist.a*180./!dpi),$
                        uname='a',uvalue='a')
              button=widget_button(base3,value='Default tilt')

          base4=widget_base( Estbase,/row)
              text=widget_text(base4,/editable,value=stringr(tlist.b*180./!dpi),$
                        uname='b',uvalue='b')
              button=widget_button(base4,value='Default rotation')

          base6=widget_base( Estbase,/row)
              text=widget_text(base6,/editable,value=stringr(tlist.center[0]),$
                        uname='xc',uvalue='xc')
              button=widget_button(base6,value='Default Xcen')
          base7=widget_base( Estbase,/row)
              text=widget_text(base7,/editable,value=stringr(tlist.center[1]),$
                        uname='yc',uvalue='yc')
              button=widget_button(base7,value='Default Ycen')

          base8=widget_base( Estbase,/row)
              text=widget_text(base8,/editable,value=stringr(tlist.dist*100),$
                        uname='ddist',uvalue='ddist')
              button=widget_button(base8,value='Default distance')
              
;          base11=widget_base( Estbase,/row)
            tmp=where(tlist.Nrings eq 1,count)
           base11a=widget_base( Estbase,/row,/exclusive,sensitive=count eq 1,uname='Tiltsign')
           button=widget_button(base11a,value='Tilt > 0',uname='a>0')
           widget_control,button,/set_button
           button=widget_button(base11a,value='Tilt < 0',uname='a<0')
;            base11a=widget_base( base11,/column,/exclusive)
;           button=widget_button(base11a,value='Distance > 0',uname='dist>0')
;           widget_control,button,/set_button
;           button=widget_button(base11a,value='Distance < 0',uname='dist<0')

          base10=widget_base( Estbase,/column)
          button=widget_button(base10,value='Initialize parameters')
          
          base9=widget_base(base,/row)
          gobutton=widget_button(base9,value='Calibrate')
          releasebutton=widget_button(base9,value='OK')
          releasebutton=widget_button(base9,value='Cancel')

          baset=widget_base(topbase,/column,uname='wrapbase')
          label=widget_label(baset,value='',uname='labelout')

          WIDGET_CONTROL, topbase, /REALIZE

          for i=0l,list.nrTiltRing do widget_control,rbutton[i],set_button=1
          for i=0l,nref-1 do widget_control,refbutton[i],set_button=refbool[i]

          Xmanager,'XRDUA_event',topbase, event_handler='tilt_event',$
          cleanup='CleanUp_tilt',GROUP_LEADER=ev.top,/NO_BLOCK

          endcase
'Restore Spatial':begin
          RestoreCorrectionev,ev
             endcase
 'Load PDF':    begin
           widget_control,ev.top,get_uvalue=list
          if n_elements(list) eq 0 then return
           LoadPDD,ev,1
           widget_control,ev.top,get_uvalue=list
          RefreshDisplay, ev.top,list
           endcase
'Load Multiple PDF':begin
              widget_control,ev.top,get_uvalue=list
          if n_elements(list) eq 0 then return
           LoadPDD,ev,2
           widget_control,ev.top,get_uvalue=list
          RefreshDisplay, ev.top,list
              endcase
'Calibration':     begin
          mode,ev,val
          endcase
'Edit d-spacing':begin
          mode,ev,val
          endcase
'Edit d-spacing ':begin
          mode,ev,val
          endcase
'Save Image':    begin
          WIDGET_CONTROL, ev.top, get_uvalue=list
          if n_elements(list) eq 0 then return
          if widget_info(list.drawdyn,/VALID_ID) then Wset,list.drawdynindex $
          else wset,list.drawindex
          error=CutPath(list.file,file=file)
          SaveImage,list.path,file+'.jpg','dplot',list,OutID=list.OutID
          endcase
'Clear Display':begin
          WIDGET_CONTROL, ev.top, get_uvalue=list
          ID=widget_info(ev.top,FIND_BY_UNAME='Show Points')
          widget_control,ID,set_value='  Show Points'
          ID=widget_info(ev.top,FIND_BY_UNAME='Show Debye Marker')
          widget_control,ID,set_value='  Show Debye Marker'
          ID=widget_info(ev.top,FIND_BY_UNAME='Show Debye Marker Axis')
          widget_control,ID,set_value='  Show Debye Marker Axis'
          ID=widget_info(ev.top,FIND_BY_UNAME='Show Tie')
          widget_control,ID,set_value='  Show SpaDist'
          ID=widget_info(ev.top,FIND_BY_UNAME='Show Tilt')
          widget_control,ID,set_value='  Show Calib'
             ID=widget_info(ev.top,FIND_BY_UNAME='Show DynMark')
          widget_control,ID,set_value='  Show DynMark'
          ID=widget_info(ev.top,FIND_BY_UNAME='Show Radial Axis')
          widget_control,ID,set_value='  Show Radial Axis'
          ID=widget_info(ev.top,FIND_BY_UNAME='Show bad pixels')
          widget_control,ID,set_value='  Show bad pixels'
          list.cross=0
          list.smellipse=0
          list.ShowTie=bytarr(4)
          list.ShowTilt=0
          if ptr_valid(list.PDFs) then (*list.PDFs)[*]=0
          list.bPdyn=0
          list.RadAxShow=0b
          if ptr_valid(list.smask) then (*list.smask)[*]=0
          list.MarkBadPixels=0
          RefreshDisplay, ev.top,list
          endcase
'File Browsing':begin
          MakeBOptionsWindow,ev
          endcase
'Auto scale':  begin
          WIDGET_CONTROL, ev.top, get_uvalue=list
          if ~ptr_valid(list.tiff) then return
          
          ; Select image for optimal clipping
          img=(*list.tiff)[list.hrange[0]:list.hrange[1],list.vrange[0]:list.vrange[1]]

          ; Get the highest upper bound that clips 99% of the pixels
          nclip=total(histogram(img,rev=rev),/cumul)
          s=size(img)
          ind=where(nclip gt ((s[1]*s[2])*0.99),ct)
          if ct eq 0 then return
          ma=(img[rev[rev[ind[0]]]])

          ; Get the highest lower bound that doesn't clip anything
          mi=min(img)

          ; Set lower and upper bound sliders
          list.clipslider[3:4]=list.clipslider[0]>(alog([mi,ma]+1)/(alog(2)*list.clipslider[2]))<list.clipslider[1]
          list.clipslider[3]>=0
          list.clipslider[3]<=list.clipslider[4]
          
          ; New scaling and update
          list.sclmax=2.^(list.clipslider[3:4]*list.clipslider[2])-1
          RefreshDisplay, ev.top,list
          RefreshScalingSlider,ev.top
          endcase
'Reverse Colors':begin
          WIDGET_CONTROL, ev.top, get_uvalue=list
          list.revcol=list.revcol*(-1)
          RefreshDisplay, ev.top,list
          RefreshScalingSlider,ev.top
          endcase
'Lower threshold':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_uvalue=uval
          bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
          bkcode[1]=ev.select
          list.CorMethod[(*(list.CorSeq)).back]=Set2DbkCode(bkcode)
          widget_control,ev.top,set_uvalue=list
          ID=widget_info(ev.top,FIND_BY_UNAME='tstrip1')
          widget_control,ID,sensitive=ev.select
          RefreshDisplay, ev.top,list
          endcase
'Higher threshold':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_uvalue=uval
          bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
          bkcode[2]=ev.select
          list.CorMethod[(*(list.CorSeq)).back]=Set2DbkCode(bkcode)
          widget_control,ev.top,set_uvalue=list
          ID=widget_info(ev.top,FIND_BY_UNAME='tstrip2')
          widget_control,ID,sensitive=ev.select
          RefreshDisplay, ev.top,list
          endcase
'Background Stripping':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_uvalue=uval
          bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
          bkcode[0]=uval
          list.CorMethod[(*(list.CorSeq)).back]=Set2DbkCode(bkcode)
          widget_control,ev.top,set_uvalue=list
          ID=widget_info(ev.top,FIND_BY_UNAME='BackBase1')
          widget_control,ID,sensitive=1
          ID=widget_info(ev.top,FIND_BY_UNAME='BackBase2')
          widget_control,ID,sensitive=0
          RefreshDisplay, ev.top,list
          endcase
'Dark Image subtraction':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_uvalue=uval
          bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
          bkcode[0]=uval
          list.CorMethod[(*(list.CorSeq)).back]=Set2DbkCode(bkcode)
          widget_control,ev.top,set_uvalue=list
          ID=widget_info(ev.top,FIND_BY_UNAME='BackBase1')
          widget_control,ID,sensitive=0
          ID=widget_info(ev.top,FIND_BY_UNAME='BackBase2')
          widget_control,ID,sensitive=1
          RefreshDisplay, ev.top,list
          endcase
'Already Subtracted':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_uvalue=uval
          bkcode=Get2DbkCode(list.CorMethod[(*(list.CorSeq)).back])
          bkcode[0]=uval
          list.CorMethod[(*(list.CorSeq)).back]=Set2DbkCode(bkcode)
          ptr_free,list.tiffb
          widget_control,ev.top,set_uvalue=list
          ID=widget_info(ev.top,FIND_BY_UNAME='BackBase1')
          widget_control,ID,sensitive=0
          ID=widget_info(ev.top,FIND_BY_UNAME='BackBase2')
          widget_control,ID,sensitive=0
          RefreshDisplay, ev.top,list
          endcase
'Load Dark Image':begin
          widget_control,ev.top,get_uvalue=list
          stopl=0

          path=''
          repeat begin
          path=CAPDIALOG_PICKFILE(path=list.darkpath,file=list.darkfile,filter='*'+list.Format,title='Load Dark Image...')
          if path EQ '' then return
          error=CutPath(path,path=darkpath,file=file,ext=ext)
          list.darkpath=darkpath
          list.darkfile=file+ext

          tiffb=ReadCCD(list.darkpath,list.darkfile,info=list.readinfo)
          if not ptr_valid(tiffb) then begin
              printw,list.OutID,'Dark image not loaded: '+list.darkpath+list.darkfile
              return
          endif
          s1=size(*tiffb)
          s2=size(*list.tiff)
          if (s1[1] ne s2[1]) or (s1[2] ne s2[2]) then printw,list.OutID,'Dark image has wrong dimensions' else stopl=1
          endrep until stopl

          ; smooth
          ptr_free,list.tiffb
          if list.bcko eq 0 then list.tiffb=tiffb $
          else begin
              list.tiffb=mSavgol(data=tiffb,r=list.bcko,nr=list.bckw/2,nl=list.bckw/2,order=[0,0])
              ptr_free,tiffb
          endelse
          *list.tiffb*=list.bckm

          widget_control,ev.top,set_uvalue=list
          endcase
'Superimpose Patterns':begin
          BatchProc2DTiff,ev,2
          endcase
'Average Patterns':begin
          BatchProc2DTiff,ev,1
          endcase
'Summate Patterns':begin
          BatchProc2DTiff,ev,3
          endcase
'Do all 3':begin
          BatchProc2DTiff,ev,0
          endcase
'Rebin Image':  begin
            n=PromptNumber('2',ev.top,'Binning (n or 1/n):')
            tmp=strsplit(n,'/',/extract,count=ct)
            case ct of
            0: return
            1: n=0>fix(n)
            2: n=-(0>fix(tmp[1]))
            else:return
            endcase

            if n eq 0 then return
            widget_control,ev.top,get_uvalue=list
            OpenFile,ev,0,fileptr=list.tiff,rebin=n
          endcase
'Make MPEG':    begin
          if xregistered("XInterAnimate") then begin
              ret=dialog_message("Only one animation at a time is allowed.",/info)
              return
          endif

          widget_control,ev.top,get_uvalue=list

          ; ----Select files to load and to save----
          path=Select_Files(ev.top,list,outpath=pathred,outfile=file,sortmethod=list.sortmethod,separator=list.sortseparator)
          if path[0] eq '' then return
          n=n_elements(path)
          Result = DIALOG_MESSAGE(stringr(n)+' files: proceed?',/question)
          if result eq 'No' then return
          ;pathmpeg=SelFile(pathred,file+'.mpeg','*.mpeg','Save MPEG....')
          ;if pathmpeg EQ '' then return
          widget_control,/hourglass

          ; ----Dimensions and scaling----
          hrange=list.hrange
          vrange=list.vrange
          if not list.zoom then begin
              dyn=0B ; Use static draw window ROI
              wset,list.drawindex
              dim=[list.stath,list.statv]
          endif else begin
              dyn=1B ; Use dynamic draw window ROI
              dim=[list.dynh,list.dynv]
          endelse

          ; ----Make sequence----
          tvlct,RoldCT,GoldCT,BoldCT,/get
          LoadctRev,0,/silent,rev=list.revcol
          XINTERANIMATE, SET=[dim,n], /SHOWLOAD;,MPEG_FILENAME=pathmpeg
          for i=0l,n-1 DO begin
              tiff=ReadCCD(path[i],format=list.Format,info=list.readinfo)
              if ptr_valid(tiff) then begin
                  if dyn then $
                  (*tiff)=congrid((*tiff)[list.hrange[0]:list.hrange[1],list.vrange[0]:list.vrange[1]], list.dynh, list.dynv) $
                  else (*tiff)=congrid(*tiff, list.stath, list.statv)
                  XINTERANIMATE, FRAME = i, IMAGE = WeakScl(*tiff,list.sclmax,list.displaypower)
              endif
              ptr_free,tiff
              printw,list.OutID,path[i]
          endfor
          XINTERANIMATE,1;,/MPEG_CLOSE
          tvlct,RoldCT,GoldCT,BoldCT
          endcase
'Experimental Setup':begin
          widget_control,ev.top,get_uvalue=list
          ExpSetup,ev=ev
          endcase
'Batch Processing':begin
          widget_control,ev.top,get_uvalue=list
          XRD_BP,ev=ev
          endcase
'Close All Windows':    begin
          DEVICE, WINDOW_STATE=wstate
          wind = where(wstate eq 1,nopen)
          if nopen ne 0 then $
             for i=0l,nopen-1 do begin
                if wind[i] lt 32 then wdelete,wind[i]
             endfor
          endcase
'Edit XDI files':begin
          widget_control,ev.top,get_uvalue=list
          XRD_XDI,ev=ev
          endcase
'Edit 1D profile':begin
          widget_control,ev.top,get_uvalue=list
          XRD_CHI,ev=ev,retID=retID
          widget_control,ev.top,get_uvalue=list
          if ptr_valid(list.CHIchild) then (*list.CHIchild)=[(*list.CHIchild),retID] $
          else begin
              list.CHIchild=ptr_new(retID)
              widget_control,ev.top,set_uvalue=list
          endelse
          endcase
'Compare 1D profiles':begin
          widget_control,ev.top,get_uvalue=list
          oplotchi,ev=ev
          endcase
'Flat Field':   begin
          mode,ev,val
          endcase
'Mask off':   begin
          mode,ev,val
          endcase
'Load Flat-Field Image':begin
          widget_control,ev.top,get_uvalue=list
          stopl=0

          path=''
          repeat begin
          path=CAPDIALOG_PICKFILE(path=list.ffpath,file=list.fffile,filter='*'+list.Format,title='Load Flat-Field Image...')
          if path EQ '' then return
          error=CutPath(path,path=ffpath,file=file,ext=ext)
          list.ffpath=ffpath
          list.fffile=file+ext

          ff=ReadCCD(list.ffpath,list.fffile,info=list.readinfo)
          if ptr_valid(ff) eq 0 then return
          s1=size(*ff)
          s2=size(*list.tiff)
          if (s1[1] ne s2[1]) or (s1[2] ne s2[2]) then printw,list.OutID,'Flat-Field image has wrong dimensions' else stopl=1
          endrep until stopl

          (*ff)=mean(*ff)/(*ff)
          ptr_free,list.ff
          list.ff=ff

          widget_control,ev.top,set_uvalue=list
          mode,ev,val
          endcase
'Reset Flat-Field Image':begin
          widget_control,ev.top,get_uvalue=list
          ptr_free,list.ff
          widget_control,ev.top,set_uvalue=list
          mode,ev,val
          endcase
'Reset Mask Image':begin
          widget_control,ev.top,get_uvalue=list
          ptr_free,list.maskoff
          list.maskofffile=''
          widget_control,ev.top,set_uvalue=list
          mode,ev,val
          endcase
'Valid pixels:':begin
            widget_control,ev.top,get_uvalue=list
            list.maskoffdyn=ev.select
            ID=widget_info(ev.top,FIND_BY_UNAME='maskoffdynthres')
            widget_control,ID,editable=list.maskoffdyn
            widget_control,ev.top,set_uvalue=list
          endcase
'Load Mask Image':begin
          widget_control,ev.top,get_uvalue=list
          stopl=0

          path=''
          filter='*'+FileFormat(2)
          sortformat,filter,list.Format
          repeat begin
          path=CAPDIALOG_PICKFILE(path=list.maskoffpath,file=list.maskofffile,filter=filter,title='Load Mask Image...')
          if path EQ '' then return
          error=CutPath(path,path=maskoffpath,file=file,ext=ext)
          list.maskoffpath=maskoffpath
          list.maskofffile=file+ext

          maskoff=ReadCCD(list.maskoffpath,list.maskofffile,info=list.readinfo)
          if ptr_valid(maskoff) eq 0 then return
          s1=size(*maskoff)
          s2=size(*list.tiff)
          if (s1[1] ne s2[1]) or (s1[2] ne s2[2]) then printw,list.OutID,'Mask image has wrong dimensions' else stopl=1
          endrep until stopl

          list.maskoffthres='> '+stringr(long(min(*maskoff)))
          thres=PromptNumber(list.maskoffthres,ev.top,'Mask image valid pixels:')
          list.maskoffthres=thres
          compareself,maskoff,list.maskoffthres
          ptr_free,list.maskoff
          list.maskoff=maskoff

          widget_control,ev.top,set_uvalue=list
          mode,ev,val
          endcase
'Save Valid Pixels':begin
          widget_control,ev.top,get_uvalue=list
          if ptr_valid(list.tiffvalid) then begin
            struc={center:list.center,dist:list.dist,file:list.maskofffile,scrx:list.scrx,$
                lambda:list.lambda,darkfile:''}
            path=SelFile(list.maskoffpath,'mask.tif','*'+FileFormat(2),'Save CCD...')
            if path EQ '' then return
            error=WriteCCD(list.tiffvalid,path,struc=struc)
            if ~error then printw,list.OutID,'Error in saving '+path
          endif
          endcase
'Apply Flat-Field Correction':begin
          T=systime(1)
          widget_control,ev.top,get_uvalue=list
          if ptr_valid(list.ff) eq 0 then return
          type=size(*list.tiff,/type)
          (*list.tiff)*=(*list.ff)
          ;p=TypeConvert(list.tiff,type)
          list.CorDone[((*(list.CorSeq))).ff]=1
          widget_control,ev.top,set_uvalue=list
          RefreshDisplay, ev.top,list
          mode,ev,val
          t2=systime(1)
          printw,list.OutID,'Flat Field: time (sec):'+string(t2-T)
          endcase
'Apply Mask':begin
          T=systime(1)
          widget_control,ev.top,get_uvalue=list
          
          if list.maskoffdyn then begin
                if not ptr_valid(list.tiffvalid) then list.tiffvalid = ptr_new(compare(list.tiff,list.maskoffdynthres)) $
                else (*list.tiffvalid) and= compare(list.tiff,list.maskoffdynthres)
            endif
                
            if ptr_valid(list.maskoff) then begin
                if not ptr_valid(list.tiffvalid) then list.tiffvalid = ptr_new(*list.maskoff) $
                else (*list.tiffvalid) and= (*list.maskoff)
            endif
            
            if ~ptr_valid(list.tiffvalid) then return
            
          list.CorDone[((*(list.CorSeq))).maskoff]=1
          widget_control,ev.top,set_uvalue=list
          RefreshDisplay, ev.top,list
          mode,ev,val
          t2=systime(1)
          printw,list.OutID,'Mask off: time (sec):'+string(t2-T)
          endcase
'Restore Flat-Field':begin
          RestoreCorrectionev,ev
          endcase
'Restore Mask':begin
          RestoreCorrectionev,ev
          endcase
'Fast mode':begin
          widget_control,ev.top,get_uvalue=list
          list.bFastSpline=ev.select
          ptr_free,list.SpaDisStruc
          widget_control,ev.top,set_uvalue=list
          endcase
'Flip vertical':begin
          widget_control,ev.top,get_uvalue=list
          list.vFlipSpline=ev.select
          widget_control,ev.top,set_uvalue=list
          endcase
'Delete All Tie Points':begin
          widget_control,ev.top,get_uvalue=list
          list.nrTieRing=0
          (*list.nrTie)=[0]
          ptr_free,list.TieI
          ptr_free,list.TieO
          ptr_free,list.TieRingD
          RefreshDisplay, ev.top,list
          endcase
'Spline file':begin
          widget_control,ev.top,get_uvalue=list
          if widget_info(list.drawdyn,/VALID_ID) then widget_control,list.drawdyn,$
            draw_button=0
          ID1=widget_info(ev.top,FIND_BY_UNAME='WithGridBase')
          ID2=widget_info(ev.top,FIND_BY_UNAME='WithoutGridBase')
          ID3=widget_info(ev.top,FIND_BY_UNAME='SPbase')
          ID4=widget_info(ev.top,FIND_BY_UNAME='Add Tie Point')
          ID5=widget_info(ev.top,FIND_BY_UNAME='Delete Tie Point')
          ID6=widget_info(ev.top,FIND_BY_UNAME='Select area')
          ID7=widget_info(ev.top,FIND_BY_UNAME='Delete All Tie Points')
          list.TieGrid=2
          widget_control,ID1,sensitive=0
          widget_control,ID2,sensitive=0
          widget_control,ID3,sensitive=1
          widget_control,ID4,sensitive=0
          widget_control,ID5,sensitive=0
          widget_control,ID6,sensitive=0
          widget_control,ID7,sensitive=0
          list.nrTieRing=0
          (*list.nrTie)=[0]
          ptr_free,list.TieI
          ptr_free,list.TieO
          ptr_free,list.TieRingD
          ptr_free,list.SpaDisStruc
          RefreshDisplay, ev.top,list
          endcase
'With Grid':    begin
          widget_control,ev.top,get_uvalue=list
          if widget_info(list.drawdyn,/VALID_ID) then widget_control,list.drawdyn,$
            Event_Pro='TIE_PROCESS_EVENTS',draw_button=1
          ID1=widget_info(ev.top,FIND_BY_UNAME='WithGridBase')
          ID2=widget_info(ev.top,FIND_BY_UNAME='WithoutGridBase')
          ID3=widget_info(ev.top,FIND_BY_UNAME='SPbase')
          ID4=widget_info(ev.top,FIND_BY_UNAME='Add Tie Point')
          ID5=widget_info(ev.top,FIND_BY_UNAME='Delete Tie Point')
          ID6=widget_info(ev.top,FIND_BY_UNAME='Select area')
          ID7=widget_info(ev.top,FIND_BY_UNAME='Delete All Tie Points')
          list.TieGrid=ev.select
          widget_control,ID1,sensitive=ev.select
          widget_control,ID2,sensitive=ev.select eq 0
          widget_control,ID3,sensitive=0
          widget_control,ID4,sensitive=1
          widget_control,ID5,sensitive=1
          widget_control,ID6,sensitive=1
          widget_control,ID7,sensitive=1
          list.nrTieRing=0
          (*list.nrTie)=[0]
          ptr_free,list.TieI
          ptr_free,list.TieO
          ptr_free,list.TieRingD
          ptr_free,list.SpaDisStruc
          RefreshDisplay, ev.top,list
          endcase
'Without Grid': begin
          widget_control,ev.top,get_uvalue=list
          if widget_info(list.drawdyn,/VALID_ID) then widget_control,list.drawdyn,$
            Event_Pro='TIE_PROCESS_EVENTS',draw_button=1
          ID1=widget_info(ev.top,FIND_BY_UNAME='WithGridBase')
          ID2=widget_info(ev.top,FIND_BY_UNAME='WithoutGridBase')
          ID3=widget_info(ev.top,FIND_BY_UNAME='SPbase')
          ID4=widget_info(ev.top,FIND_BY_UNAME='Add Tie Point')
          ID5=widget_info(ev.top,FIND_BY_UNAME='Delete Tie Point')
          ID6=widget_info(ev.top,FIND_BY_UNAME='Select area')
          ID7=widget_info(ev.top,FIND_BY_UNAME='Delete All Tie Points')
          list.TieGrid=ev.select eq 0
          widget_control,ID1,sensitive=ev.select eq 0
          widget_control,ID2,sensitive=ev.select
          widget_control,ID3,sensitive=0
          widget_control,ID4,sensitive=1
          widget_control,ID5,sensitive=1
          widget_control,ID6,sensitive=1
          widget_control,ID7,sensitive=1
          list.nrTieRing=0
          (*list.nrTie)=[0]
          ptr_free,list.TieI
          ptr_free,list.TieO
          ptr_free,list.TieRingD
          ptr_free,list.SpaDisStruc
          RefreshDisplay, ev.top,list
          endcase
'Find Tie Points':begin
          widget_control,/hourglass
          easy_grid,ev
;          widget_control,ev.top,get_uvalue=list
;          if (list.masknr eq 0) then return
;          surf_int,list,2,ev.top,'Select area of Grid Points'
          endcase
'Grid Properties':begin
          widget_control,ev.top,get_uvalue=list
          if (*list.nrTie)[0] eq 0 then return
          ptr_free,list.SpaDisStruc
          widget_control,ev.top,sensitive=0
          fixcon=1B

          base=widget_base(/row,title=val,uvalue={top:ev.top,fixcon:fixcon})
          baset=widget_base(base,/column)
          xs=150

          baserow=widget_base(baset,/row)
          label=widget_label(baserow,value='Input Tie points:',xsize=xs)
          text=widget_text(baserow,value=stringr((*list.nrTie)[0]),editable=0,uvalue='nGP',scr_xsize=xs,uname='nGP')

          baserow=widget_base(baset,/row)
          label=widget_label(baserow,value='Output Grid nodes:',xsize=xs)
          text=widget_text(baserow,value=stringr(list.GridDim[0]*list.GridDim[1]),$
                  editable=0,uname='nGPmax',scr_xsize=xs)

          baserow=widget_base(baset,/row)
          label=widget_label(baserow,value='Grid Dimension:',xsize=xs)
          text=widget_text(baserow,value=stringr(list.GridDim[0])+',' $
            +stringr(list.GridDim[1]),editable=1,uvalue='GridDim',scr_xsize=xs)

          baserow=widget_base(baset,/row)
          label=widget_label(baserow,value='Grid spacing (mm):',xsize=xs)
          text=widget_text(baserow,value=stringr(list.GridSpReal[0])+',' $
            +stringr(list.GridSpReal[1]),editable=1,uvalue='GridSpReal',scr_xsize=xs)

          baserow=widget_base(baset,/row)
          label=widget_label(baserow,value='Grid offset (pixels):',xsize=xs)
          text=widget_text(baserow,value=stringr(list.GridOffset[0])+',' $
            +stringr(list.GridOffset[1]),editable=1,uvalue='GridOffset',uname='GridOffset',scr_xsize=xs)
            
          baserow=widget_base(baset,/row)
          label=widget_label(baserow,value='Grid origin (pixels):',xsize=xs)
          text=widget_text(baserow,value=stringr(list.GridOrigin[0])+',' $
            +stringr(list.GridOrigin[1]),editable=1,uvalue='GridOrigin',uname='GridOrigin',scr_xsize=xs)
            
          baserow=widget_base(baset,/row)
          label=widget_label(baserow,value='Grid tilt (degrees):',xsize=xs)
          text=widget_text(baserow,value=stringr(list.GridTilt) $
            ,editable=1,uvalue='GridTilt',uname='GridTilt',scr_xsize=xs)
            
          baserow=widget_base(baset,/row)
          label=widget_label(baserow,value='Pixel size ('+msymbols('micro')+'m):',xsize=xs)
          text=widget_text(baserow,value=stringr(list.scrx*10^6.)+',' $
            +stringr(list.scry*10^6.),editable=1,uvalue='PixSize',uname='PixSize',scr_xsize=xs)

          ;button=widget_button(baset,value='Sort') Not useful anymore
          button=widget_button(baset,value='Calculate Output Nodes')
          button=widget_button(baset,value='AutoReconnect')

          baset=widget_base(baset,/column,/nonexclusive)
          bfix=widget_button(baset,value='Fix connections when calculating output nodes')
          bhold=lonarr(4)
          bhold[0]=widget_button(baset,value='View Input Points')
          bhold[1]=widget_button(baset,value='View Output Points')
          bhold[2]=widget_button(baset,value='View I-O Connects')
          bhold[3]=widget_button(baset,value='View O-O Connects')

          baset=widget_base(base,/column)
          
          button=widget_button(baset,value='Find Double Output')
          button=widget_button(baset,value='Remove Overlapping Input')
          button=widget_button(baset,value='Move Column')
          button=widget_button(baset,value='Move Row')

          dind=DoublePresent(list.TieO,coord=coord)
          label=widget_label(baset,value='Double Output')
          table=widget_table(baset,all_events=1,$
          X_SCROLL_SIZE=2,Y_SCROLL_SIZE=10,uname='TieO',uvalue=1,$
          COLUMN_LABELS=['Output X-coord','Output Y-coord'],value=coord)

          baset=widget_base(baset,/row,/exclusive)
          button1=widget_button(baset,value='Reconnect')
          button2=widget_button(baset,value='Delete')

          WIDGET_CONTROL, base, /REALIZE
          for i=0l,3 do widget_control,bhold[i],set_button=list.ShowTie[i]
          widget_control,bfix,set_button=fixcon
          widget_control,button1,set_button=1
          Xmanager,'XRDUA_event',base, event_handler='RealGrid_event',$
          cleanup='CleanUp_RealGrid',GROUP_LEADER=ev.top,/no_block
          endcase
'Refresh Display':begin
          widget_control,ev.top,get_uvalue=list
          RefreshDisplay, ev.top,list
          endcase
'Corrections':begin
          CorrectOptions,ev.top
          endcase
'Read Options': ReadOptions,ev.top
'Interface Options':begin
          PointOptions,ev.top
          endcase
'Area overlap statistics':begin
          CreateOverlapTable,ev.top
          endcase
'View overlap statistics':begin
          EditOverlap
          endcase
'Convert CCD format': ConvertCCDimages,ev
'Prepare for tomorecon: SkyScan':begin
          widget_control,ev.top,get_uvalue=list
          SkyScan_NRecon,ev.top,list,$
            sortmethod=list.sortmethod,separator=list.sortseparator
          endcase
'Nearest Neighbor':begin
          widget_control,ev.top,get_uvalue=list
          list.mipol=1
          widget_control,ev.top,set_uvalue=list
          endcase
'Bilinear':begin
          widget_control,ev.top,get_uvalue=list
          list.mipol=2
          widget_control,ev.top,set_uvalue=list
          endcase
'Cubic':  begin
          widget_control,ev.top,get_uvalue=list
          list.mipol=3
          widget_control,ev.top,set_uvalue=list
          endcase
'Azimuthal default':begin
          widget_control,ev.top,get_uvalue=list
          list.phipixphi0=azMmtoPixel(-list.b,list.scrx,list.scry) ; causes the third angle to be zero
          ID=widget_info(ev.top,FIND_BY_UNAME='phipixphi0')
          widget_control,ID,set_value=stringr(list.phipixphi0*180/!dpi)
          widget_control,ev.top,set_uvalue=list
          endcase
    endcase
    endcase
2:    begin;slider event
    widget_control,ev.id,get_uvalue=uvalue
    WIDGET_CONTROL, ev.id, GET_VALUE=val
    WIDGET_CONTROL, ev.top, get_uvalue=list
    if n_elements(list) eq 0 then return
    case uvalue of
    'clip slider min':  begin
            list.clipslider[3]=val
            list.clipslider[4]>=val
            list.sclmax=2.^(list.clipslider[3:4]*list.clipslider[2])-1
            RefreshDisplay, ev.top,list
            RefreshScalingSlider,ev.top
            endcase
    'clip slider max':  begin
            list.clipslider[3]<=val
            list.clipslider[4]=val
            list.sclmax=2.^(list.clipslider[3:4]*list.clipslider[2])-1
            RefreshDisplay, ev.top,list
            RefreshScalingSlider,ev.top
            endcase
    'clip slider power':  begin
            minmax=widget_info(ev.id,/SLIDER_MIN_MAX)
            list.displaypower=float(val)/minmax[1]
            RefreshDisplay, ev.top,list
            RefreshScalingSlider,ev.top
            endcase
    'circle slider':begin
            list.cslider[3]=val
            ID=widget_info(ev.top,FIND_BY_UNAME='csliderlabel')
            azi=list.cslider[3]*list.cslider[2]
            WIDGET_CONTROL, ID,SET_VALUE='Azimuth: '+stringr(azi)+msymbols('degrees')
            azi*=!pi/180.
            
            ; Size of the dynamic viewport
            wsizex=list.stath<list.dynh
            wsizey=list.statv<list.dynv
            
            ; Viewport center
            viewcenter=[wsizex,wsizey]/2.
            
            ; Center of dynamic viewport in image coordinates
            widget_control,list.drawdyn,get_draw_view=viewoffset ; offset off the viewport with respect to the dynamic window
            viewoffset>=0
            XY=viewcenter+viewoffset ; center of the viewport in the dynamic window
            XY=DynTR(XY,list.sf,[list.hrange[0],list.vrange[0]]) ; center of the dynamic viewport in image coordinates

            ; Get radius and azimuth of the new viewport center
            validphi=1b
            phi=azPixeltoMm(azi,list.scrx,list.scry)
            if list.cslideruseM then begin
                R=BraggTtoX(list.ttmellipse,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,phi=phi,/onlyp,/nonortp,validphi=validphi)
            endif else begin
                tt=BraggXYtoX(XY[0],XY[1],0,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/onlyt)
                ; tt should be the same when "zoom" is not touched (except for rounding errors)
                R=BraggTtoX(tt,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,phi=phi,/onlyp,/nonortp,validphi=validphi)
            endelse
            if finite(R) eq 0 or ~validphi then return

            ; New viewport center in image coordinates
            XY=abs(R)*[cos(azi),sin(azi)]+list.center
            
            ; New viewport center in dynamic window coordinates
            XYnew=RTDyn(XY,list.sf,[list.hrange[0],list.vrange[0]])
            
            ; Within scroll range?
            scrollmaxx=(list.dynh-wsizex)>0 ; maximum x-offset
            scrollmaxy=(list.dynv-wsizey)>0 ; maximum y-offset
            viewoffsetnew=XYnew-viewcenter
            bscroll=floor(viewoffsetnew) gt 0 and ceil(viewoffsetnew) lt [scrollmaxx,scrollmaxy]
            bscroll=(bscroll[0] and bscroll[1])
            if list.zoom eq 1 then bscroll=0
            
            ; Scroll or move
            if bscroll then begin
                list.Pdyn=XY
                  WIDGET_CONTROL, list.drawdyn, SET_DRAW_VIEW=viewoffsetnew
            endif else begin
                ; New zoom area
                m=(list.hrange[1]-list.hrange[0])/2
                list.hrange=XY[0]+[-m,m]
                if (list.hrange[0] lt 0) then list.hrange-=list.hrange[0]
                if (list.hrange[1] ge (*list.tiffs)[1]) then list.hrange-=list.hrange[1]-(*list.tiffs)[1]+1
                
                m=(list.vrange[1]-list.vrange[0])/2
                list.vrange=XY[1]+[-m,m]
                if (list.vrange[0] lt 0) then list.vrange-=list.vrange[0]
                if (list.vrange[1] ge (*list.tiffs)[2]) then list.vrange-=list.vrange[1]-(*list.tiffs)[2]+1
                
                if (list.hrange[0] lt 0) or ((list.hrange[1] ge (*list.tiffs)[1])) or $
                   (list.vrange[0] lt 0) or ((list.vrange[1] ge (*list.tiffs)[2])) then return
                    
                xs=(list.hrange[1]-list.hrange[0]+1)
                ys=(list.vrange[1]-list.vrange[0]+1)
                list.dynh=round(xs/list.sf)
                list.dynv=round(ys/list.sf)
                ZOOM_CHANGE,ev,list,xs,ys
                
                ; Set scroll in the middle
                widget_control,list.drawdyn,SET_DRAW_VIEW=[scrollmaxx,scrollmaxy]/2
                
                ; Save new viewport center in image coordinates
                widget_control,list.drawdyn,get_draw_view=viewoffset
                list.Pdyn=DynTR((viewoffset>0)+viewcenter,list.sf,[list.hrange[0],list.vrange[0]])                    
            endelse
            WIDGET_CONTROL, ev.top, set_uvalue=list

            ; dynamic marker plot
            if list.bPdyn ne 0 then begin
                Wset,list.drawindex
                tplot,list
                Wset,list.pixindex1
                Device, Copy = [0,0, list.stath,list.statv, 0,0,list.drawindex]
            endif

            endcase
    endcase
    endcase
3:     BEGIN ; text_widget
    widget_control,ev.id,get_uvalue=uvalue
    case uvalue of
    'center':  begin
          widget_control,ev.top,get_uvalue=list
          center=list.center
          widget_control,ev.id,get_value=str
          str=strcompress(str,/remove_all)
          str=str[0]
          starr=strsplit(str,',',/extract,count=count)
          if count ne 2 then return
          list.center=float(starr)
          RefreshDisplay, ev.top,list
          endcase
    'wstrip':  begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          wb=0
          reads,str,wb
          list.wb=wb
          widget_control,ev.top,set_uvalue=list
          endcase
    'istrip':  begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          niter=0
          reads,str,niter
          list.niter=niter
          widget_control,ev.top,set_uvalue=list
          endcase
    'tstrip1':  begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          thres=0.
          reads,str,thres
          list.thres[0]=thres
          widget_control,ev.top,set_uvalue=list
          endcase
   'tstrip2':  begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          thres=0.
          reads,str,thres
          list.thres[1]=thres
          widget_control,ev.top,set_uvalue=list
          endcase
    'bcko':  begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          bcko=0
          reads,str,bcko
          list.bcko=bcko
          widget_control,ev.top,set_uvalue=list
          endcase
    'bckw':  begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          bckw=0
          reads,str,bckw
          list.bckw=bckw
          widget_control,ev.top,set_uvalue=list
          endcase
    'bckm':  begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          bckm=0.
          reads,str,bckm
          list.bckm=bckm
          widget_control,ev.top,set_uvalue=list
          endcase
    'MEllipsed':begin
          widget_control,ev.id,get_value=dspac
          MEllipseSet,float(dspac),0,ev.top
          endcase
    'MEllipset':begin
          widget_control,ev.id,get_value=ttheta
          MEllipseSet,float(ttheta),1,ev.top
          endcase
   'MEllipsefr':begin
          widget_control,ev.id,get_value=R
          MEllipseSet,float(R),2,ev.top
          endcase
    'MEllipseq':begin
          widget_control,ev.id,get_value=Q
          MEllipseSet,float(Q),3,ev.top
          endcase
    'zingerthres':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          ZingerThreshold=0.
          reads,str,ZingerThreshold
          list.ZingerThreshold=ZingerThreshold
          widget_control,ev.top,set_uvalue=list
          endcase
    'zingerwidth':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          ZingerWidth=0
          reads,str,ZingerWidth
          list.ZingerWidth=ZingerWidth
          widget_control,ev.top,set_uvalue=list
          endcase
    'saturmax':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          saturmax=0.
          reads,str,saturmax
          list.SaturMax=saturmax
          widget_control,ev.top,set_uvalue=list
          endcase
    'saturnmax':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          saturnmax=0l
          reads,str,saturnmax
          list.SaturNMax=saturnmax
          widget_control,ev.top,set_uvalue=list
          endcase
    'maskoffdynthres':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          list.maskoffdynthres=str
          widget_control,ev.top,set_uvalue=list
          endcase
    'saturaspect':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=str
          saturaspect=0.
          reads,str,saturaspect
          list.SaturAspect=saturaspect
          widget_control,ev.top,set_uvalue=list
          endcase
    'window scale':begin
          widget_control,ev.top,get_uvalue=list

          ; Middle of the dynamic screen in image coordinates
          widget_control,list.drawdyn,get_draw_view=XY
          XY=DynTR((XY>0)+[(list.stath<list.dynh)/2.,(list.statv<list.dynv)/2.],list.sf,[list.hrange[0],list.vrange[0]])

          widget_control,ev.id,get_value=val
          sf=float(val)
          if sf eq 0 then begin
              sf=1.
              widget_control,ev.id,set_value=stringr(sf)
          endif
          list.sf=1./sf

          if list.zoom then begin
             xs=list.hrange[1]-list.hrange[0]+1
             ys=list.vrange[1]-list.vrange[0]+1
          endif else begin
             xs=(*list.tiffs)[1]
             ys=(*list.tiffs)[2]
          endelse
          list.dynh=round(xs/list.sf)
          list.dynv=round(ys/list.sf)

          drawdyn=list.drawdyn
          drawdynindex=list.drawdynindex
          pixindex2=list.pixindex2
          dynh=list.dynh
          dynv=list.dynv
          sf=list.sf
          tmp=CreateDyn(drawdyn,drawdynindex,pixindex2,dynh,dynv,list.stath,list.statv,xs,ys,$
            sf,widget_info(ev.top,FIND_BY_UNAME='drawbasedyn'),ev.id)
          list.drawdyn=drawdyn
          list.drawdynindex=drawdynindex
          list.pixindex2=pixindex2
          list.dynh=dynh
          list.dynv=dynv
          list.sf=sf

          x=RTDyn(XY[0],list.sf,list.hrange[0])
          y=RTDyn(XY[1],list.sf,list.vrange[0])
          WIDGET_CONTROL, list.drawdyn, SET_DRAW_VIEW=[(x-list.stath/2.)>0,(y-list.statv/2.)>0]

          RefreshDisplay, ev.top,list
          endcase
'phipixphi0':begin
         widget_control,ev.top,get_uvalue=list
         widget_control,ev.id,get_value=str
         list.phipixphi0=float(str)*!pi/180
          RefreshDisplay, ev.top,list
          endcase
'Pff':    begin
        widget_control,ev.top,get_uvalue=list
        widget_control,ev.id,get_value=str
        list.IntCor2Dinfo.Pff=float(str)
        widget_control,ev.top,set_uvalue=list
        endcase
    else:
    endcase
    endcase
4:    begin; draw widget
    widget_control,ev.top,get_uvalue=list
    RefreshCircleSlider,ev.top,list
    widget_control,ev.top,set_uvalue=list
    endcase
8:    begin ; droplist
    widget_control,ev.id,get_uvalue=uvalue
    case uvalue of
    'fftype':    begin
                widget_control,ev.top,get_uvalue=list
                list.IntCor2Dinfo.fftype=ev.index
                widget_control,ev.top,set_uvalue=list
                ID=widget_info(ev.top,FIND_BY_UNAME='labelPff')
                widget_control,ID,set_value=LabelPff(list.IntCor2Dinfo.fftype)
                endcase
    else:
    endcase
    endcase
else:
endcase

end;pro XRDUA_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRDUA

print,'=========================================START XRDUA'

ControlDevice,Platform,Visual

result=ReadINI('$Scale:')
cslider=[0.,360,1,0]
clipslider=[0.,(result.(0))[3],(result.(1))[1],(result.(0))[0],(result.(0))[1]]
revcol=(result.(0))[2]
sclmax=2.^(clipslider[3:4]*clipslider[2])-1
sf=1./((result.(1))[0]<1)
displaypower=(result.(1))[2]

result2=ReadINI('$StatSize:')
stath=result2.(0)
statv=stath

base = WIDGET_BASE(MBAR=bar,/row,title='XRDUA:',/TLB_SIZE_EVENTS,uname='XRDUAtop')
    base1=widget_base(base,/column,sensitive=0,uname='base1')
       label3= widget_label(base1,value='IMAGE SCALING:')
       baset=widget_base(base1)
       slider1=widget_slider(base1,minimum=clipslider[0],maximum=clipslider[1],uname='clip slider min',uvalue='clip slider min',value=clipslider[0]>clipslider[3]<clipslider[1])
       slider1=widget_slider(base1,minimum=clipslider[0],maximum=clipslider[1],uname='clip slider max',uvalue='clip slider max',value=clipslider[0]>clipslider[4]<clipslider[1])
       slider1=widget_slider(base1,minimum=0,maximum=100,uname='clip slider power',uvalue='clip slider power',value=round(displaypower*100))
       button34=widget_button(base1,value='Auto scale')
       button34=widget_button(base1,value='Reverse Colors')
       label4=widget_label(base1,value='WINDOW SCALING:')
       text1=widget_text(base1,value=stringr((result.(1))[0]<1),uname='window scale',uvalue='window scale',/editable)
       label2= widget_label(base1,value='BROWSE:')
       button21= widget_button(base1,value='Next file')
       button22= widget_button(base1,value='Previous file')
       label1= widget_label(base1,value='MODE:')
       button5= widget_button(base1,value='Zoom',uname='Zoom')
       button7= widget_button(base1,value='Beam position',uname='Beam position')
       button11= widget_button(base1,value='Mask Setup',uname='Mask Setup')
       button23= widget_button(base1,value='Spatial Distortion',uname='Spatial Distortion')
       button27= widget_button(base1,value='Calibration',uname='Calibration')
       button27= widget_button(base1,value='Azimuthal offset',uname='Azimuthal offset')
       button30= widget_button(base1,value='Debye Marker',uname='Debye Marker')
       button17= widget_button(base1,value='Background',uname='Background')
       button27= widget_button(base1,value='Flat Field',uname='Flat Field')
       button17= widget_button(base1,value='Saturation',uname='Saturation')
       button55= widget_button(base1,value='Zingers',uname='Zingers')
       button27= widget_button(base1,value='Mask off',uname='Mask off')
       
       geom = Widget_Info(button27, /Geometry)
       drawslider = WIDGET_DRAW(baset,uname='drawslider',button_events=1,event_pro='SL_PROCESS_EVENTS',xsize=geom.scr_xsize,ysize=1.5*geom.scr_ysize,uvalue={xoff:0.,yoff:1./6,xs:geom.scr_xsize,ys:geom.scr_ysize,type:-1})

    base2=widget_base(base,/column)
       base3=widget_base(base2,/row,uname='drawbase')
            base7=widget_base(base3,/column,uname='drawbasestat')
         draw = WIDGET_DRAW(base7, XSIZE=stath, YSIZE=statv,uname='draw')
         ; make dummy window1
         Window, /Free, xsize=stath, ysize=statv, /Pixmap
         pixindex1 = !D.Window
         base8=widget_base(base3,/column,uname='drawbasedyn')
       base6=widget_base(base2,/row, xsize=stath)
         label5=widget_label(base6,value='Azimuth: 0'+msymbols('degrees'),uname='csliderlabel', xsize=stath)
       base5=widget_base(base2,/row,uname='csliderbase',sensitive=0)
           slider2=widget_slider(base5,minimum=0,maximum=360,uvalue='circle slider',$
           uname='circle slider',value=0,xsize=stath,/DRAG,SCROLL=1,/SUPPRESS_VALUE)
       base4=widget_base(base2,/column,uname='dynbase')

    menu1 = WIDGET_BUTTON(bar, VALUE='File', /MENU)
       button1 = WIDGET_BUTTON(menu1, VALUE='Read 2D Pattern',uname='Read 2D Pattern')
       button50= widget_button(menu1,value='Read Options')
       button43 = WIDGET_BUTTON(menu1, VALUE='Write 2D Pattern',uname='Write 2D Pattern',sensitive=0)
       button36= widget_button(menu1,value='Auto Load',uname='Auto Load',sensitive=0)
       button12 = WIDGET_BUTTON(menu1, VALUE='Load Mask',uname='Load Mask',/separator,sensitive=0)
       button12 = WIDGET_BUTTON(menu1, VALUE='Add Mask',uname='Add Mask',sensitive=0)
       button13 = WIDGET_BUTTON(menu1, VALUE='Save Mask',uname='Save Mask',sensitive=0)
       button13 = WIDGET_BUTTON(menu1, VALUE='Save Valid Pixels',uname='Save Valid Pixels',sensitive=0)
       button26 = WIDGET_BUTTON(menu1, VALUE='Load PDF',uname='Load PDF',sensitive=0)
       button26 = WIDGET_BUTTON(menu1, VALUE='Load Multiple PDF',uname='Load Multiple PDF',sensitive=0)
       button29 = WIDGET_BUTTON(menu1, VALUE='Save Image',uname='Save Image',/separator,sensitive=0)
       button2 = WIDGET_BUTTON(menu1, VALUE='Exit',uname='Exit',/separator)

    menu2 = WIDGET_BUTTON(bar, VALUE='Options', /MENU,sensitive=0,uname='menu2')
       button6 = WIDGET_BUTTON(menu2, VALUE='Experimental Geometry')
       button33= widget_button(menu2,value='File Browsing')
       button49= widget_button(menu2,value='Corrections')
       button50= widget_button(menu2,value='Interface Options')
    menu5 = WIDGET_BUTTON(bar, VALUE='Perform', /MENU,uname='menu5')
       button15= widget_button(menu5,value='Azimuthal Integration',uname='menu50a',sensitive=0)
       button16= widget_button(menu5,value='Surface Integration',uname='menu51',sensitive=0)
       button15= widget_button(menu5,value='Rebin Image',uname='menu57',sensitive=0)
       
       button15= widget_button(menu5,value='Automatic Correction',uname='menu50c',sensitive=0,ACCELERATOR = "F6",/separator)
       button15= widget_button(menu5,value='Automatic Integration',uname='menu50b',sensitive=0,ACCELERATOR = "F7")
       
       button35= widget_button(menu5,value='Superimpose Patterns',/separator,uname='menu53',sensitive=0)
       button40= widget_button(menu5,value='Average Patterns',uname='menu54',sensitive=0)
       button41= widget_button(menu5,value='Summate Patterns',uname='menu55',sensitive=0)
       button41= widget_button(menu5,value='Do all 3',uname='menu58',sensitive=0)
       
       button37= widget_button(menu5,value='Experimental Setup',/separator)
       button38= widget_button(menu5,value='Batch Processing',ACCELERATOR = "F4")
       button42= widget_button(menu5,value='Edit XDI files',ACCELERATOR = "F3")
       button43= widget_button(menu5,value='Edit 1D profile',ACCELERATOR = "F1")
       button57= widget_button(menu5,value='Compare 1D profiles',ACCELERATOR = "F2")
       
    menu6 = WIDGET_BUTTON(bar, VALUE='Custom', /MENU,uname='menu6')
        button53= widget_button(menu6,value='Area overlap statistics',uname='menu52',sensitive=0)
        button54= widget_button(menu6,value='View overlap statistics')
        button45= widget_button(menu6,value='Make MPEG',uname='menu56',sensitive=0,/separator)
        button37= widget_button(menu6,value='Prepare for tomorecon: SkyScan')
        button16= widget_button(menu6,value='Rectangle summation',uname='menu59',sensitive=0)
        button16= widget_button(menu6,value='Convert CCD format')
        

    menu4 = WIDGET_BUTTON(bar, VALUE='Display', /MENU,sensitive=0,uname='menu4')
       button8= widget_button(menu4,value='Surface',uname='Surface')
       button9= widget_button(menu4,value='  Show Points',/separator,uname='Show Points' )
       button31= widget_button(menu4,value='  Show Debye Marker',uname='Show Debye Marker')
       button56= widget_button(menu4,value='  Show Debye Marker Axis',uname='Show Debye Marker Axis')
       button24= widget_button(menu4,value='  Show SpaDist',uname='Show Tie')
       button28= widget_button(menu4,value='  Show Calib',uname='Show Tilt')
       button46= widget_button(menu4,value='  Show bad pixels',uname='Show bad pixels')
       button52= widget_button(menu4,value='  Show DynMark',uname='Show DynMark')
       button52= widget_button(menu4,value='  Show Radial Axis',uname='Show Radial Axis')
       button25= widget_button(menu4,value='View PDF')
       button10= widget_button(menu4,value='View Mask')
       button32= widget_button(menu4,value='Clear Display',/separator)
       button48= widget_button(menu4,value='Refresh Display',ACCELERATOR = "F5")

    menu6 = WIDGET_BUTTON(bar, VALUE='System', /MENU,uname='menu6')
       button39= widget_button(menu6,value='Close All Windows')
       button14 = widget_button(menu6, value='Configuration File',uname='Configuration File')
       button14 = widget_button(menu6, value='Save Configuration',uname='Save Configuration')
       button14 = widget_button(menu6, value='Register DLMs',uname='Register DLMs')

    menu3=WIDGET_BUTTON(bar, VALUE='Help', /MENU)
       button3 = WIDGET_BUTTON(menu3, VALUE='Manual',uname='Manual')
       button3 = WIDGET_BUTTON(menu3, VALUE='Online Support',uname='Online Support')
       button3 = WIDGET_BUTTON(menu3, VALUE='Online Update',uname='Online Update')
       button4 = widget_button(menu3, value='About XRDUA',uname='About XRDUA')

; ----Arrange OutID----
OutID=SubscribeToOutID()

; ----Make datastructure----
args = COMMAND_LINE_ARGS(COUNT=argsct)
bopenfile=argsct ne 0
if bopenfile then begin
    error=CutPath(args[0],path=path,file=file,ext=ext)
    file+=ext
    Format=ext
endif else begin
    result=ReadINI('$XRD path:')
    path=result.(0)
    file=result.(1)
    PPos = RSTRPOS(file, '.')
    Format=strmid(file,PPos)
endelse

result=ReadINI('$XRD mskpath:')
mskpath=result.(0)
mskfile=result.(1)

result=ReadINI('$XRD pddpath:')
pddpath=result.(0)
pddfile=result.(1)
darkpath = path
darkfile = ''
ffpath = path
fffile = ''
maskoffpath = path
maskofffile = ''

result=ReadINI('$XRD_CHI path:')
chipath=result.(0)
chifile=result.(1)

result=ReadINI('$XRD_XDI path:')
xdipath=result.(0)
xdifile=result.(1)

result=(ReadINI('$Show:')).(0)
smellipse=result[0]
cross=result[1]
ShowTie=result[2:5]
ShowTilt=result[6]
MarkBadPixels=result[7]

result=ReadINI('$XRD tiepath:')
TiePath=result.(0)
TieFile=result.(1)

result=(ReadINI('$Pointer:')).(0)
DelMass=byte(result[0])
DelWidth=result[1]

result=(ReadINI('$CW:')).(0)
CW1=result[0]
CW2=result[1]

result=(ReadINI('$CSlider:')).(0)
cslideruseM=result[0]
bPdyn=result[1]

result=(ReadINI('$Misc:')).(0)
letter=result[0]
OverLapCrit=result[1]
mipol=fix(result[2])

result=ReadINI('$CrossUpdate:')
CrossUpdate=result.(0)

CorSeq=PTR_NEW({zinger:1B,satur:2B,back:4B,ff:3B,spatial:5B,maskoff:0B})
SeqTags=STRLOWCASE(tag_names(*CorSeq))
autocor=(ReadINI('$AutoCor:')).(0)

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

readinfo=ReadINI('$ReadCCDInfo:')

dplotthick=(ReadINI('$DplotThick:')).(0)

TLBsize=widgettreesizeinit(base,'draw')

RadAx=ReadINI('$RadAx:')

compilers=AllCompilerConfigurations()

list={listtype:0    ,$                ; listtype: 0=initial, 1=normal

    tiff:PTR_NEW(),$                   ; 2D tiff image
    tiffb:PTR_NEW(),$               ; Background
    ff:PTR_NEW(),$                    ; inverse Flat-field image scaled so mean=1
    maskoff:PTR_NEW(),$                ; binarized mask off image
    maskpolygon:PTR_NEW(fltarr(2)),$ ; mask polygon in progress            
    tiffvalid:PTR_NEW(),$            ; Mask image(spa.dist.)
    tiffs:PTR_NEW(),$                  ; Info on size and format of image
    revcol:revcol,$                    ; 1 => don't reverse colors, -1 => reverse
    readinfo:readinfo,$             ; Normalize EDF on reading

    CorMethod:lonarr(6),$              ; zinger removal: 0=no, 1=yes
                                    ; saturation removal: 0=no, 1=yes
                                       ; background correction: 0=no, 1=stripping, 2=dark image
                                       ; spatial distortion corr.: 0=no, 1=yes
                                       ; flat field correction: 0=no, 1=yes
                                       ; mask off: 0=no, 1=yes
    CorDone:bytarr(6),$               ; Corrections are performed
    CorSeq:CorSeq,$                 ; Positions of the fields in the structure gives sequence,
                                       ; the value in the field gives the position in CorMethod/Done
    SeqTags:SeqTags,$                ; Tag names of CorSeq
    autocor:autocor,$                ; Perform all required corrections on load

    AzIntInfo:AzIntInfo,$            ; Azimuthal integration info
    batchreturnp:ptr_new(),$        ; Last azimuthal integration information
    wb:3,$                           ; Strip-width for background calc
    niter:20,$                      ; Number of iterations for background calc
    bcko:0,$                        ; Order of dark-image smoothing
    bckw:5,$                        ; Width of dark-image smoothing
    bckm:1.,$                        ; Multiplication factor
    thres:[0.,0.],$                    ; Lower/Upper Background thresholding
    maskoffthres:'> 0',$            ; Mask image threshold
    maskoffdynthres:'> 0',$        ; Mask image dynamic threshold (on the image itself)
    maskoffdyn:0b,$                    ; do dynamic mask off
    
    path:path,$                     ; Directory of tiff image
    file:file,$                     ; Name of tiff image
    xdinr:0l,$                        ; Image index in a multi page image
    chipath:chipath,$                ; Path of .chi
    chifile:chifile,$                ; 1D profile name
    xdipath:xdipath,$                ; xdi path
    xdifile:xdifile,$                ; xdi file
    pddpath:pddpath,$               ; Path of .pdd file
    pddfile:pddfile,$               ; Powder diffraction markers from .pdd file
    mskpath:mskpath,$               ; Path of .msk file
    mskfile:mskfile,$               ; Display, fit and other options from .msk file
    darkpath:darkpath,$             ; Path of dark image
    darkfile:darkfile,$             ; Dark image file
    ffpath:ffpath,$                    ; Flat-field path
    fffile:fffile,$                    ; Flat-field image file
    maskoffpath:maskoffpath,$        ; Mask-off path
    maskofffile:maskofffile,$        ; Mask image file
    Format:Format,$                 ; Format of CCD patterns

    Platform:Platform,$                ; Platform
    Visual:Visual,$                    ; Name of the visual class
    mode:'',$                        ; Dynamic window mode

    sortseparator:'_',$                ; for sorting files
    sortmethod:1b,$                 ; 1 => sort strings, 0 sort strings with offset

    lambda:0.,$                     ; X-ray lambda(msymbols('angstroms'))
    cconv:double(4.13566743*299792458.*10.^(-8)),$; E(keV)=12.3984/lambda(msymbols('angstroms'))
    scrx:0.,$                         ; Pixel x-size (m, input in um)
    scry:0.,$                         ; Pixel y-size (m, input in um)
    dist:0.,$                         ; Distance sample-detector (m, input in cm)

    thetadyn:0,$                      ; theta of dynamic cursor when defining arc mask
    dangle:0.,$                     ; for passing 180 when defining arc mask
    circle2:1b,$                    ; 2 circles or 1

    P1:[-1L,-1L],$                  ; point to determine center circle
    P2:[-1L,-1L],$                  ; point to determine center circle
    P3:[-1L,-1L],$                  ; point to determine center circle
    bP:1,$                             ; which point to change P1,P2 or P3
    CW1:CW1,$                       ; Size of Points in stat.
    CW2:CW2,$                         ; Size of Points in dyn.
    letter:letter,$                 ; Charsize of labels
    center:[-1.,-1.],$              ; coordinates (2D pixel index= 2D coord in tiffmatrix) of center
    cross:cross,$                      ; show point and center
    DelMass:DelMass,$                ; Prompt when deleting multiple points or not
    DelWidth:DelWidth,$                ; Half width of point detection
    Pdyn:[-1L,-1L],$                ; point on stat. to indicate center of dyn.
    bPdyn:bPdyn,$                    ; show Pdyn?

    clipslider:clipslider,$            ; min, max, increment, val_slidermin, val_slidermax
    sclmax:sclmax,$                 ; minimum/maximum value for grayscaling: 2^(min_or_max*clipsliderincr)-1
    displaypower:displaypower,$        ; image^power before appying sclmax
    cslider:cslider,$                ; min, max, increment and value for circle slider
    cslideruseM:cslideruseM,$        ; use Mellipse for cslider
    
    RadAxType:byte((RadAx.(0))[0]),$; Radial axis type
    RadAxShow:byte((RadAx.(0))[1]),$; Display radial axis?
    RadAxMin:(RadAx.(0))[2],$        ; Radial axis minimum
    RadAxMax:(RadAx.(0))[3],$        ; Radial axis maximum
    RadAxAzimuth:(RadAx.(0))[4],$    ; Radial axis azimuth
    RadAxBins:long((RadAx.(0))[5]),$; Radial axis bins
    RadAxFormat:(RadAx.(1))[0],$    ; Radial axis format    

    smask:PTR_NEW(),$                  ; show mask areas
    masknr:0,$                         ; number of areas in mask
    mask:PTR_NEW(),$                   ; first column: maskpart type (1:square, 2:arc, 3: circle)
                                     ; second column: index of first part of a linked group of parts
                                       ; third column: 1=set, 0=not set
                                       ; rest: parameters of maskpart
                                       ;     square: x1,x2,y1,y2 => LLcorner(1) and URcorner(2) coord
                                       ;     arc:tt1,azi1,tt2,azi2 (azi defined in cone frame!!!)
                                       ;     circle:tt1,tt2,0,0
    MaskP:-1L,$                        ; Highlight this mask group
    OverLapCrit:OverLapCrit,$        ; threshold for grouping areas by there overlap
    AreaRangeStruc:{i:0,b:[0b,1b],str:fltarr(4)},$; Remember area range definition settings

    xs:0L, $                         ; X static corner of the zoom box.
    ys:0L, $                         ; Y static corner of the zoom box.
    xd:0L, $                         ; X dynamic corner of the zoom box.
    yd:0L,$                          ; Y dynamic corner of the zoom box.
    sfh:0.,$                          ; scale factor for overviewwindow horizontal
    sfv:0.,$                          ; scale factor for overviewwindow vertical
    sf:sf,$                           ; scale factor for dynamic window
    hrange:lonarr(2),$              ; horizontal view-range in image coord. if the same => not zoomed
    vrange:lonarr(2),$              ; vertical view-range in image coord.
    zoom:0b,$                          ; zoomed or not
    stath:stath,$                      ; horizontal stat. window size
    statv:statv,$                      ; vertical stat. window size
    dynh:0L,$                          ; horizontal dyn. window size
    dynv:0L,$                          ; vertical dyn. window size
    dplotthick:dplotthick,$            ; line thickness in dynamic window

    draw:draw,$                     ; ID of overview drawwidget
    drawindex:0L,$                   ; index of overview drawwidget
    drawdyn:0L,$                       ; ID of dynamic scrollwindow
    drawdynindex:0L,$                ; index of dynamic scrollwindow
    pixindex2:0L,$                   ; index of pixwindow2 (contains dynamic window)
    pixindex1:pixindex1,$           ; index of pixwindow1 (contains overview window)
    TLBsize:TLBsize,$                ; Resize info

    smellipse:smellipse,$           ; first bit: show ellipse, second: show axis
    ttmellipse:0.,$                 ; two-theta of mark ellipse (in radians)

    SaturMax:0.,$                    ; Maximum Pattern value
    SaturAspect:2.,$                ; Saturation aspect ratio
    SaturNMax:0l,$                    ; Saturation maximal area
    ZingerWidth:3,$                    ; Smoothing width in zinger removal
    ZingerThreshold:1.2,$            ; Ratio threshold in zinger removal
    MarkBadPixels:MarkBadPixels,$     ; Show bad pixels on display

    nrTieRing:0,$                   ; Number of Tiepoint rings
    nrTie:PTR_NEW([0l]),$              ; Number of Tie points in each ring
    TieRingD:PTR_NEW(),$            ; d-spacings for closed rings
    TieI:PTR_NEW(),$                 ; Tie point coordinates of inputimage
    TieO:PTR_NEW(),$                 ; Tie point coordinates for output image
    ShowTie:ShowTie,$               ; Show Tiepoints or not: [TieI,TieO,ConnectIO,ConnectOO]
    TiePath:Tiepath,$                ; Fit2D spline file
    TieFile:TieFile,$                ; Fit2D spline file
    tck1:PTR_NEW(),$                ; X-distortion
    tck2:PTR_NEW(),$                ; Y-distortion
    TieGrid:1b,$                       ; Tie points on grid(1) or circles(0) or fit2d spline(2)
    bFastSpline:1b,$                ; Perform spline spatial distortion in fast mode
    vFlipSpline:0b,$                ; Flip image before spatial distortion (when using fit2d spline file)

    GridDim:[51,51],$               ; Number of grid points horizontal and vertical
    GridSpReal:[70/31.,70/31.],$     ; Spacing between grid points in mm
    GridSpPix:[0.,0.],$               ; Spacing between grid points in pixels
                                    ; will be calculated as soon as we have the image dimensions
    GridOffset:[0.,0.],$             ; Border offset of the grid
    GridOrigin:[0.,0.],$             ; Origin of the grid in the image coordinate system
    GridTilt:0.,$                    ; Tilt of the grid in the image coordinate system
    
    SpaDisStruc:PTR_NEW(),$            ; Spatial distortion correction information
    WarpStruc:PTR_NEW(),$            ; Azimuthal integration information

    a:0.,$                          ; a Tilt angle: RotXc(a) => this is a tilt in it's own frame
    b:0.,$                          ; b Rotation angle: RotZt(b)=> this is a rotation in it's own frame
    phipixphi0:0.,$                 ; azimuthal shift due to 2-angle calibration
    IntCor2Dinfo:DefaultIntCor2Dinfo(),$; 2D intensity correction parameters
    nrTiltRing:0,$                     ; Number of rings used
    nrTilt:PTR_NEW([0l]),$             ; Number of points in each ring
    TiltI:PTR_NEW(),$                  ; point coordinates of inputimage
    ShowTilt:ShowTilt,$             ; Show points or not
    TiltRingD:PTR_NEW(),$           ; d-spacings for closed rings
    mipol:mipol,$                    ; 1: NN, 2: bilinear, 3: cubic

    nPDF:0,$                        ; number of phases
    PDFn:PTR_NEW(),$                 ; number of d-spacings for each phase
    PDFnTotal:0,$                    ; total number of d-spacings
    PDFd:PTR_NEW(),$                   ; d-spacings in phases
    PDFni:PTR_NEW(),$                  ; number of intensities in phases
    PDFi:PTR_NEW(),$                   ; intensities in phases
    PDFs:PTR_NEW(),$                   ; first bit: show lines or not, second bit: show labels or not
    PDFname:PTR_NEW(),$               ; names of loaded phases
    PDFhkl:PTR_NEW(),$                ; miller indices
    PDFSearch:{type:0,PDFSubstr:'',PDFdi:0.,PDFde:0.},$    ; PDF subsearch

    v:3,$                             ; side windows width for Top-hat (peaksearch)
    w:4,$                              ; center window width for Top-hat (peaksearch)
    r:3,$                              ; criterium for peaksearch
    fit:0,$                          ; fit with: real(0),suggestive(1),no(2) constraints

    OutID:OutID,$                    ; Output log ID
    CHIchild:PTR_NEW(),$            ; CHI Child ID's
    CrossUpdate:CrossUpdate,$        ; Update CHIChild ?

    compilers:compilers,$            ; compiler settings for DLMS
    
    QueTres:0.1,$                      ; Que-treshold for rejecting events (in seconds)
                                       ; In this time the overflow que events must be rejected
    QueTimeStamp:Systime(1)}        ; browser Que-flush criterium


WIDGET_CONTROL, base, /REALIZE,set_uvalue=list
;widget_control,/DELAY_DESTROY
Xmanager,'XRDUA',base,/NO_BLOCK,cleanup='CleanUp_XRDUA'

RefreshScalingSlider,base
if bopenfile then OpenFile,{top:base,id:base},0,/info
end ;XRDUA
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

pro RefreshSetupWindow,ev,list
res=BraggTtoX(list.ttmellipse,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/pixm)
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsed')
widget_control,ID,set_value=stringr(res[1])
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsefr')
widget_control,ID,set_value=stringr(res[0])
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipseq')
widget_control,ID,set_value=stringr(res[2])
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipset')
widget_control,ID,set_value=stringr(list.ttmellipse*180./!dpi)
ID=widget_info(ev.top,FIND_BY_UNAME='warning')
tmp=ConicTypeTT(list.a,list.ttmellipse,str=str)
widget_control,ID,set_value=str

ID1=widget_info(ev.top,FIND_BY_UNAME='dist')
widget_control,ID1,set_value=stringr(list.dist*100)

ID1=widget_info(ev.top,FIND_BY_UNAME='ttrange')
widget_control,ID1,set_value=$
              stringr(list.ttrange[0]*180/!dpi)+' - '+$
              stringr(list.ttrange[1]*180/!dpi)
res=BraggTtoX(list.ttrange,list.lambda,/onlyd)
ID1=widget_info(ev.top,FIND_BY_UNAME='drange')
widget_control,ID1,set_value= $
              stringr(res[0])+' - '+$
              stringr(res[1])

ID1=widget_info(ev.top,FIND_BY_UNAME='xen')
          widget_control,ID1,set_value=stringr(list.cconv/list.lambda)
          ID1=widget_info(ev.top,FIND_BY_UNAME='xlambda')
          widget_control,ID1,set_value=stringr(list.lambda)

end;pro RefreshSetupWindow
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ELLIPSEM2_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]
IF thisEvent NE 'DOWN' THEN RETURN
widget_control,ev.top,get_uvalue=list
Widget_Control, ev.id, Event_Pro='ELLIPSEM2_DRAW', $
    draw_motion=1
Widget_Control, ev.top, Set_UValue=list
LoadctRev,3,/silent
end;pro ELLIPSEM2_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ELLIPSEM2_DRAW ,ev

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

; calcultate image position of cursor
x=0>StatTR(ev.x,list.sfh)<(list.nx-1)
y=0>StatTR(ev.y,list.sfv)<(list.ny-1)
ret=BraggXYtoX(x,y,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/pixm)
tt=ret[1]

azi=atan(y-list.center[1],x-list.center[0])
azi*=180/!pi
if azi lt 0 then azi+=360

Rflat=BraggTtoX(tt,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,/pixm,/onlyp)

ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsed')
widget_control,ID,set_value=stringr(ret[2])
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipseq')
widget_control,ID,set_value=stringr(ret[3])
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsep')
widget_control,ID,set_value=stringr(x)+','+stringr(y)
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsefr')
widget_control,ID,set_value=stringr(Rflat)
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipser')
widget_control,ID,set_value=stringr(ret[0])
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipset')
widget_control,ID,set_value=stringr(tt*180./!dpi)
ID=widget_info(ev.top,FIND_BY_UNAME='MEllipsea')
widget_control,ID,set_value=stringr(azi)
ID=widget_info(ev.top,FIND_BY_UNAME='warning')
tmp=ConicTypeTT(list.a,tt,str=str)
widget_control,ID,set_value=str

IF thisEvent EQ 'DOWN' THEN BEGIN
    WSet, list.drawindex
    Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]
    Widget_Control, list.draw, draw_motion=0, $
       Event_Pro='ELLIPSEM2_PROCESS_EVENTS'

    list.ttmellipse=tt
    list.smellipse=list.smellipse or 1

    RefreshDisplaySetup,list
    Widget_Control, ev.top, Set_UValue=list
    return
endif

; whipe out previously drawn ellipse
WSet, list.drawindex
Device, Copy = [0,0, list.stath,list.statv, 0,0,list.pixindex1]

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

Widget_Control, ev.top, Set_UValue=list
end;pro ELLIPSEM2_DRAW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro splot, list
LoadctRev,0,/silent,rev=list.revcol

;----Clean window----
h=!d.x_size
v=!d.y_size
tv,bytarr(h,v)

;----Display point to determine middle----
point=list.center
if point[0] ne -1 and point[1] ne -1 then begin
x=RTstat(point[0],list.sfh)
y=RTstat(point[1],list.sfv)
Plots,[x-list.CW2,x+list.CW2],y,/device
Plots,x,[y-list.CW2,y+list.CW2],/device
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
                           x2=RTstat([xEllipse[*,0],xEllipse[*,0]+1.5*list.CW1],list.sfh)
                           y2=RTstat([xEllipse[*,1],xEllipse[*,1]+0.5*list.CW1],list.sfv)
                           Plots,x2,y2,/device,color=col
                           xyouts,x2[1],y2[1],label[j],color=col,/device,charsize=list.letter
                       endif
               endif
            endfor
        endif
        off+=(*list.PDFn)[i]
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

end; pro splot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshDisplaySetup,list
Wset,list.drawindex
splot,list
Wset,list.pixindex1
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.drawindex]
end;pro RefreshDisplaySetup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_Exps, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
heap_free,list
wdelete,list.pixindex1

error=UnSubscribeToOutID()
end;pro CleanUp_Exps
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Exps_event,ev

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
    endcase
0:  begin
    widget_control,ev.id,get_uvalue=list2
    list.draw=list2.draw
    list.pixindex1=list2.pixindex1
    widget_control,list2.draw,get_value=drawindex
    list.drawindex=drawindex
    widget_control,ev.id,/destroy
    RefreshDisplaySetup,list
    endcase
1:  begin
    widget_control,ev.id,get_value=val
    case val of
'Exit'     :begin
          WIDGET_CONTROL, ev.top, /DESTROY
          return
          endcase
'Save Image':   begin
          Wset,list.drawindex
          error=CutPath(list.pddfile,file=file)
          path=list.pddpath
          SaveImage,path,file+'.jpg',OutID=list.OutID
          if path eq '' then return
          error=CutPath(path,path=path,file=file)
          path=path+file+'.txt'

          if openw_safe(lun, path, width=500) then begin
              printw,list.OutID,path+' not saved'
              return
          endif
          printf,lun,'$D:'
          printf,lun,list.nPDF
          off=0
          for i=0l,list.nPDF-1 do begin
              ind=where(list.PDFs[off:off+list.PDFn[i]-1],ct)
              printf,lun,list.PDFname[i]
            printf,lun,ct
              if ct ne 0 then printf,lun,list.PDFd[ind+off]
            off=off+list.PDFn[i]
          endfor
          printf,lun,'$SOURCE ('+msymbols('angstroms')+'):'
          printf,lun,list.lambda
          printf,lun,'$PIXSIZE (m):'
          printf,lun,list.scrx,list.scry
          printf,lun,'$NPIXELS:'
          printf,lun,list.nx,list.ny
          printf,lun,'$NPIXELSDISPLAY:'
          wset,list.drawindex
          h=!d.x_size
          v=!d.y_size
          printf,lun,h,v
          printf,lun,'$DET-SAM (m):'
          printf,lun,list.dist
          printf,lun,'$ANGLES ('+msymbols('degrees')+'):'
          printf,lun,[list.a,list.b]*180./!dpi
          free_lun,lun
          printw,list.OutID,path+' saved'
          endcase
'Load PDF':begin
          LoadPDD,ev,1
          widget_control,ev.top,get_uvalue=list
          RefreshDisplaySetup,list
          widget_control,ev.top,TLB_SET_TITLE=list.title
          endcase
'Range (E fixed)--> Best fitting distance':begin
          ;list.dist=
          RefreshDisplaySetup,list
          RefreshSetupWindow,ev,list
          endcase
'Distance and E--> Range':begin
          xmap=indgen(list.nx)#replicate(1,list.ny)
          ymap=replicate(1,list.nx)#indgen(list.ny)
          rmap1=(xmap-((list.nx-1)/2.))^2.+(ymap-((list.ny-1)/2.))^2.
          rmap2=(xmap-list.center[0])^2.+(ymap-list.center[1])^2.
          ind=where((rmap1 le ((list.nx-1)/2.)^2.) and (rmap2 ge (list.nx-1)*0.5))
          tt=BraggXYtoX(xmap,ymap,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/onlyt)
          m=min(tt[ind],max=ma)
          list.ttrange=[m,ma]
          RefreshSetupWindow,ev,list
          endcase
'Range (distance fixed)--> Best fitting E':begin
          ;list.lambda=
          RefreshDisplaySetup,list
          RefreshSetupWindow,ev,list
          endcase
'Show Ellipse':begin
          list.smellipse=ev.select
          if ev.select eq 1 then widget_control,list.draw,Event_Pro='ELLIPSEM2_PROCESS_EVENTS',draw_button=1 $
          else widget_control,list.draw,draw_button=0
          RefreshDisplaySetup,list
          endcase
'View PDF':begin
           EditPDFWindow,ev,'RefreshDisplaySetup'
           return
           endcase
    endcase
    endcase
3:  begin
    widget_control,ev.id,get_uvalue=uval
    case uval of
    'xlambda':begin
         widget_control,ev.id,get_value=xlambda
         xlambdat=0.
         reads,xlambda,xlambdat
         list.lambda=xlambdat
         RefreshDisplaySetup,list
         RefreshSetupWindow,ev,list
         endcase
    'xen': begin
         ID1=widget_info(ev.top,FIND_BY_UNAME='xlambda')
         widget_control,ev.id,get_value=xen
         xent=0.
         reads,xen,xent
         xlambdat=list.cconv/xent
         list.lambda=xlambdat
         RefreshDisplaySetup,list
         RefreshSetupWindow,ev,list
         endcase
    'scrx':    begin
         widget_control,ev.id,get_value=in
         int=0.
         reads,in,int
         list.scrx=int*10^(-6.) ; from um to m
         RefreshDisplaySetup,list
         RefreshSetupWindow,ev,list
         endcase
    'scry':    begin
         widget_control,ev.id,get_value=in
         int=0.
         reads,in,int
         list.scry=int*10^(-6.) ; from um to m
         endcase
    'nx':  begin
         widget_control,ev.id,get_value=in
         int=0.
         reads,in,int
         list.nx=int
         list.sfh=list.nx/float(list.stath)
         RefreshDisplaySetup,list
         endcase
    'ny':  begin
         widget_control,ev.id,get_value=in
         int=0.
         reads,in,int
         list.ny=int
         list.sfv=list.ny/float(list.statv)
         RefreshDisplaySetup,list
         endcase
    'dist':    begin
         widget_control,ev.id,get_value=in
         int=0.
         reads,in,int
         list.dist=int*10^(-2.) ; from cm to m
         RefreshDisplaySetup,list
         RefreshSetupWindow,ev,list
         endcase
    'drange':begin
         drange=fltarr(2)
         widget_control,ev.id,get_value=str
         str=strcompress(str,/remove_all)
         str=str[0]
         starr=str_sep(str,'-')
         reads,starr,drange
         list.ttrange=BraggDtoX(drange,list.lambda,/onlyt)
         RefreshSetupWindow,ev,list
         endcase
    'ttrange':begin
         ttrange=list.ttrange*0.
         widget_control,ev.id,get_value=str
         str=strcompress(str,/remove_all)
         str=str[0]
         starr=str_sep(str,'-')
         reads,starr,ttrange
         list.ttrange=ttrange/180*!dpi
         RefreshSetupWindow,ev,list
         endcase
   'MEllipsed':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=dspac
          res=BraggDtoX(float(dspac),list.lambda,/onlyt)
          list.ttmellipse=res
          RefreshDisplaySetup,list
          RefreshSetupWindow,ev,list
          endcase
    'MEllipseq':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=q
          res=BraggQtoX(float(q),list.lambda,/onlyt)
          list.ttmellipse=res
          RefreshDisplaySetup,list
          RefreshSetupWindow,ev,list
          endcase
    'MEllipset':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=ttheta
          list.ttmellipse=float(ttheta)/180.*!dpi
          RefreshDisplaySetup,list
          RefreshSetupWindow,ev,list
          endcase
   'MEllipsefr':begin
          widget_control,ev.top,get_uvalue=list
          widget_control,ev.id,get_value=R
          R=float(R)
          res=BraggPtoX(R,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/onlyt,/pixm)
          list.ttmellipse=res
          RefreshDisplaySetup,list
          RefreshSetupWindow,ev,list
          endcase
    'a':    begin
         widget_control,ev.id,get_value=in
         int=0.
         reads,in,int
         list.a=int/180*!dpi
         RefreshDisplaySetup,list
         RefreshSetupWindow,ev,list
         endcase
     'b':    begin
         widget_control,ev.id,get_value=in
         int=0.
         reads,in,int
         list.b=int/180*!dpi
         RefreshDisplaySetup,list
         RefreshSetupWindow,ev,list
         endcase
     'xcen':    begin
         widget_control,ev.id,get_value=in
         int=0.
         reads,in,int
         list.center[0]=int
         RefreshDisplaySetup,list
         RefreshSetupWindow,ev,list
         endcase
     'ycen':    begin
         widget_control,ev.id,get_value=in
         int=0.
         reads,in,int
         list.center[1]=int
         RefreshDisplaySetup,list
         RefreshSetupWindow,ev,list
         endcase
    endcase
    endcase
endcase
widget_control,ev.top,set_uvalue=list
end; pro Exps_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ExpSetup,ev=ev

if keyword_set(ev) then begin
widget_control,ev.top,get_uvalue=list
lambda=list.lambda
cconv=list.cconv
scrx=list.scrx
scry=list.scry
drange=fltarr(2)
ttrange=fltarr(2)
dist=list.dist
sfv=list.sfv
sfh=list.sfh
nPDF=list.nPDF
PDFn=list.PDFn
PDFd=list.PDFd
PDFhkl=list.PDFhkl
PDFni=list.PDFni
PDFi=list.PDFi
PDFname=list.PDFname
PDFnTotal=list.PDFnTotal
PDFs=list.PDFs
title='Experimental Setup: '
pddpath=list.pddpath
pddfile=list.pddfile
Platform=list.Platform
CW2=list.CW2
center=list.center
revcol=-1
smellipse=0B
ttmellipse=0.
ys=list.ys
xs=list.xs
a=list.a
b=list.b
phipixphi0=list.phipixphi0
h=list.stath
v=list.statv
if ptr_valid(list.tiffs) then begin
    s=*list.tiffs
    nx=s[1]
    ny=s[2]
endif else begin
    nx=h
    ny=v
endelse
endif else begin
ControlDevice,Platform,Visual,Msize
lambda=0D
cconv=double(4.13566743*299792458.*10.^(-8))
scrx=0D
scry=0D
drange=fltarr(2)
ttrange=fltarr(2)
dist=0D
result2=ReadINI('$StatSize:')
h=result2.(0)
v=h
nx=h
ny=v
sfh=nx/float(h)
sfv=ny/float(v)

title='Experimental Setup: '
pddpath=''
pddfile=''
Platform=Platform
center=[nx/2,ny/2]
revcol=-1
smellipse=0B
ttmellipse=0.
ys=0.
xs=0.
CW2=10
a=0.
b=0.
phipixphi0=0.
endelse

OutID=SubscribeToOutID()

list2={listtype:50,$
       lambda:lambda,$
       cconv:cconv,$
       scrx:scrx,$
       scry:scry,$
       xs:xs,$
       ys:ys,$
       draw:0L,$
       drawindex:0L,$
       pixindex1:0L,$
       title:title,$
       ttrange:ttrange,$
       smellipse:smellipse,$
       ttmellipse:ttmellipse,$
       dist:dist,$
       a:a,$
       b:b,$
       phipixphi0:phipixphi0,$
       nx:nx,$
       ny:ny,$
       center:center,$
       revcol:revcol,$
       stath:h,$
       statv:v,$
       sfv:sfv,$
       sfh:sfh,$
       pddpath:pddpath,$
       pddfile:pddfile,$
       Platform:Platform,$
       OutID:OutID,$
       CW2:CW2,$
       nPDF:0,$                        ; number of phases
    PDFn:PTR_NEW(),$                 ; number of d-spacings for each phase
    PDFnTotal:0,$                    ; total number of d-spacings
    PDFd:PTR_NEW(),$                   ; d-spacings in phases
    PDFni:PTR_NEW(),$                  ; number of intensities in phases
    PDFi:PTR_NEW(),$                   ; intensities in phases
    PDFs:PTR_NEW(),$                   ; first bit: show lines or not, second bit: show labels or not
    PDFname:PTR_NEW(),$               ; names of loaded phases
    PDFhkl:PTR_NEW(),$                ; miller indices
    PDFSearch:{type:0,PDFSubstr:'',PDFdi:0.,PDFde:0.}}    ; PDF subsearch


base=widget_base(title=title,uvalue=list2,/row)
    base2=widget_base(base,/row)

       base1=widget_base(base2,/column)
         xs=170

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='X-ray energy (keV)  :        ',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(cconv/lambda),uvalue='xen',uname='xen')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='X-ray wavelength ('+msymbols('angstroms')+')  :      ',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(lambda),uvalue='xlambda',uname='xlambda')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='Pixel X-size ('+msymbols('micro')+'m) :          ',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(scrx*10^6.),uvalue='scrx',uname='scrx')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='Pixel Y-size ('+msymbols('micro')+'m) :          ',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(scry*10^6.),uvalue='scry',uname='scry')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='Distance Sample-Detector(cm):',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(dist*10^2.),uvalue='dist',uname='dist')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='CenterX(pixels):',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(center[0]),uvalue='xcen',uname='xcen')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='CenterY(pixels):',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(center[1]),uvalue='ycen',uname='ycen')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='Tilt angle('+msymbols('degrees')+'):',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(a*180./!dpi),uvalue='a',uname='a')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='Rotation angle('+msymbols('degrees')+'):',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(b*180./!dpi),uvalue='b',uname='b')

    button= widget_button(base1,value='Load PDF')
    button= widget_button(base1,value='Save Image')
    button= widget_button(base1,value='Range (E fixed)--> Best fitting distance')
    button= widget_button(base1,value='Range (distance fixed)--> Best fitting E')
    button= widget_button(base1,value='Distance and E--> Range')
    button= widget_button(base1,value='Exit')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='Pixels in X-direction :    ',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(nx),uvalue='nx',uname='nx')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='Pixels in Y-direction :    ',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(ny),uvalue='ny',uname='ny')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='D-spacing range ('+msymbols('angstroms')+'):',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(lambda/2./sin(ttrange[1]/2.))+' - '+$
         stringr(lambda/2./sin(ttrange[0]/2.)),uvalue='drange',uname='drange')

         baset=widget_base(base1,/row)
         label=widget_label(baset,value='Two-Theta range ('+msymbols('degrees')+'):',xsize=xs)
         text=widget_text(baset,/editable,value= $
         stringr(ttrange[0]*180/!dpi)+' - '+$
         stringr(ttrange[1]*180/!dpi),uvalue='ttrange',uname='ttrange')


    base3=widget_base(base2,/column,uname='drawbase')
    draw = WIDGET_DRAW(base3, XSIZE=h, YSIZE=v,uname='draw',button_events=1)
    ; make dummy window1
    Window, /Free, xsize=h, ysize=v, /Pixmap
    pixindex1 = !D.Window
    basedummy=widget_base(base,/row,uvalue={draw:draw,pixindex1:pixindex1})

    base7=widget_base(base3,/column,/nonexclusive)
         button=widget_button(base7,value='Show Ellipse')

    button=widget_button(base3,value='View PDF')

    base3=widget_base(base3,/row)
    base5=widget_base(base3,/column)
            xs=130

            res=BraggTtoX(list2.ttmellipse,list2.lambda,list2.dist,list2.scrx,list2.scry,list2.a,list2.b,list2.phipixphi0,/pixm)

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
            text=widget_text(baset,/editable,value=stringr(list2.ttmellipse*180./!dpi),$
              uname='MEllipset',uvalue='MEllipset')

            baset=widget_base(base5,/row)
            label=widget_label(baset,value='Flat Rad. dist. (mm):',xsize=xs)
            text=widget_text(baset,/editable,value=stringr(res[0]),$
              uname='MEllipsefr',uvalue='MEllipsefr')
              
            tmp=ConicTypeTT(list2.a,list2.ttmellipse,str=vstr)
            label=widget_label(base5,value=vstr,uname='warning',xsize=xs)

base5=widget_base(base3,/column)
            
            baset=widget_base(base5,/row)
            label=widget_label(baset,value='Position (display pixels):',xsize=xs)
            text=widget_text(baset,uname='MEllipsep',uvalue='MEllipsep')

            baset=widget_base(base5,/row)
            label=widget_label(baset,value='Rad. dist. (mm):',xsize=xs)
            text=widget_text(baset,uname='MEllipser',uvalue='MEllipser')

            baset=widget_base(base5,/row)
            label=widget_label(baset,value='Azimuth ('+msymbols('degrees')+'):',xsize=xs)
            text=widget_text(baset,uname='MEllipsea',uvalue='MEllipsea')
            
widget_control,draw,get_value=drawindex
WIDGET_CONTROL, base, /REALIZE
WIDGET_CONTROL, basedummy, TIMER=0
;widget_control,/DELAY_DESTROY
Xmanager,'ExpSetup',base, event_handler='Exps_event',/no_block,cleanup='CleanUp_Exps'
end;pro ExpSetup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
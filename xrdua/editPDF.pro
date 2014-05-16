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

pro ApplyEditPDF,list

widget_control,list.top,get_uvalue=listmain
ind=where(list.del ne 1,ct)
if ct eq 0 then begin ; Nothing left
    listmain.nPDF=0
    ptr_free,listmain.PDFn
    listmain.PDFnTotal=0
    ptr_free,listmain.PDFd
    ptr_free,listmain.PDFni
    ptr_free,listmain.PDFi
    ptr_free,listmain.PDFhkl
    ptr_free,listmain.PDFname
    ptr_free,listmain.PDFs
    if list.od then begin
        ptr_free,listmain.PDFx
        ptr_free,listmain.PDFi2
        ptr_free,listmain.PDFScale
        ptr_free,listmain.PDFOffset
    endif
endif else begin ; Delete some Phases or Keep everything
      listmain.nPDF=ct
    (*listmain.PDFn)=list.PDFn[ind]
    (*listmain.PDFni)=list.PDFni[ind]
    listmain.PDFnTotal=total(list.PDFn[ind])
    (*listmain.PDFname)=list.names[ind]

    off=0
    bfirst=1b
    for i=0l,list.nPDF-1 do begin
        if list.del[i] eq 0 then begin
            i0=off
            i1=off+list.PDFn[i]-1

            if bfirst then begin
                (*listmain.PDFd)=list.PDFd[i0:i1]
                (*listmain.PDFi)=list.PDFi[i0:i1]
                (*listmain.PDFhkl)=list.PDFhkl[*,i0:i1]
                (*listmain.PDFs)=list.select[i0:i1]
                if list.od then begin
                    (*listmain.PDFx)=list.PDFx[i0:i1,*]
                    (*listmain.PDFi2)=list.PDFi2[i0:i1]
                    (*listmain.PDFScale)=list.PDFScale[i]
                    (*listmain.PDFOffset)=list.PDFOffset[i]
                endif
                bfirst=0b
            endif else begin
                (*listmain.PDFd)=[(*listmain.PDFd),list.PDFd[i0:i1]]
                (*listmain.PDFi)=[(*listmain.PDFi),list.PDFi[i0:i1]]
                (*listmain.PDFhkl)=[[(*listmain.PDFhkl)],[list.PDFhkl[*,i0:i1]]]
                (*listmain.PDFs)=[(*listmain.PDFs),list.select[i0:i1]]
                if list.od then begin
                    (*listmain.PDFx)=[(*listmain.PDFx),list.PDFx[i0:i1,*]]
                    (*listmain.PDFi2)=[(*listmain.PDFi2),list.PDFi2[i0:i1]]
                    (*listmain.PDFScale)=[(*listmain.PDFScale),list.PDFScale[i]]
                    (*listmain.PDFOffset)=[(*listmain.PDFOffset),list.PDFOffset[i]]
                endif
            endelse

        endif
        off+=list.PDFn[i]
    endfor

endelse

if list.topid then CALL_PROCEDURE, list.proc,list.top,listmain $
else begin
       CALL_PROCEDURE, list.proc,listmain
       widget_control,list.top,set_uvalue=listmain
endelse
end;pro ApplyEditPDF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_EditPDF, ID

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
end;pro CleanUp_EditPDF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EditPDFWindowP,ev,list

if size(ev,/type) eq 8 then begin
    TLB=ev.top
    ID=widget_info(TLB,FIND_BY_UNAME='basegroups')
    WIDGET_CONTROL, ID, /DESTROY
    wrapbase=widget_info(TLB,FIND_BY_UNAME='wrapbase')
endif else begin
    wrapbase=ev
endelse

result=ReadINI('$EditMask:')

if n_elements(TLB) ne 0 then widget_control,TLB,update=0
basegroups=widget_base(wrapbase,/row,uname='basegroups',sensitive=list.Del[list.CurPhase] eq 0,$
    x_scroll_size=result.(1),/scroll)

; Make different columns, depending on how many elements you want in a column
nr=20
nbase=ceil(float(list.PDFn[list.curPhase])/nr)
base6=lonarr(nbase)
for i=0l,nbase-1 do base6[i]=widget_base(basegroups,/column,/nonexclusive)

PDFi=fix(list.PDFi[list.i0:list.i1])
ind=list.sortond?(sort(list.PDFd[list.i0:list.i1])):(sort(PDFi))
ind=reverse(ind)

; Loop over all groups and place them in the right row
j=0
code=list.labelview?2:1
for i=0l,list.PDFn[list.curPhase]-1 do begin
    hkl=fix(list.PDFhkl[*,ind[i]+list.i0])
    str=stringr(list.PDFd[ind[i]+list.i0])+'  '+$
           stringr(PDFi[ind[i]])+'  '+'('+stringr(hkl[0])+' '+$
           stringr(hkl[1])+' '+stringr(hkl[2])+')'
   button=widget_button(base6[j],value=str,$
       uvalue=ind[i],uname=stringr(ind[i]))
   widget_control,button,set_button=list.select[ind[i]+list.i0] and code
   if ((i+1) mod nr) eq 0 then j=j+1
endfor

if n_elements(TLB) ne 0 then widget_control,TLB,update=1
return,basegroups
end;function EditPDFWindowP
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshPDFWindowP,ev,list

PDFi=fix(list.PDFi[list.i0:list.i1])
ind=list.sortond?(sort(list.PDFd[list.i0:list.i1])):(sort(PDFi))
ind=reverse(ind)

code=list.labelview?2:1
for i=0l,list.PDFn[list.curPhase]-1 do begin
    ID=widget_info(ev.top,FIND_BY_UNAME=stringr(ind[i]))
       widget_control,ID,set_button=list.select[ind[i]+list.i0] and code
endfor

end;pro RefreshPDFWindowP
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ChangeButtonState,ev,list

ID=widget_info(ev.top,FIND_BY_UNAME='View All')
    tmp=where((list.select[list.i0:list.i1] and 1) eq 1,ct)
    widget_control,ID,set_button=ct gt 0
ID=widget_info(ev.top,FIND_BY_UNAME='Use Intensities')
       if list.PDFni[list.curPhase] ne 0 then widget_control,ID,set_button=1 $
    else widget_control,ID,set_button=0
ID=widget_info(ev.top,FIND_BY_UNAME='View Labels')
    tmp=where((list.select[list.i0:list.i1] and 2) eq 2,ct)
    widget_control,ID,set_button=ct gt 0
ID=widget_info(ev.top,FIND_BY_UNAME='Delete')
       widget_control,ID,set_button=list.Del[list.curPhase]
ID=widget_info(ev.top,FIND_BY_UNAME='% Cutofftxt')
     widget_control,ID,set_value='0'
ID=widget_info(ev.top,FIND_BY_UNAME='Change Label View')
    widget_control,ID,set_button=list.labelview
ID=widget_info(ev.top,FIND_BY_UNAME='sclslider')
if widget_info(ID,/VALID_ID) then widget_control,ID,set_value=list.PDFScale[list.curPhase]*100

end;pro ChangeButtonState
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditPDFHandler,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control, ev.top, get_uvalue=list

case widget_info(ev.id,/type) of
1 : begin ;button event
    widget_control,ev.id,get_value=val
    case val of
      'Change Label View':begin
              list.labelview=ev.select
              widget_control,ev.top,set_uvalue=list
            RefreshPDFWindowP,ev,list
            ChangeButtonState,ev,list
               endcase
      'Sort on d-spacing':begin
              list.sortond=ev.select
              widget_control,ev.top,set_uvalue=list

              ; Display new buttons
            basegroups=EditPDFWindowP(ev,list)
              endcase
      'View Labels':begin
               if ev.select then list.select[list.i0:list.i1]=list.select[list.i0:list.i1] or 2 $
               else list.select[list.i0:list.i1]=list.select[list.i0:list.i1] and 253
               if list.labelview then begin
                 for i=0l,list.PDFn[list.curPhase]-1 do begin
                      ID=widget_info(ev.top,FIND_BY_UNAME=stringr(i))
                      widget_control,ID,set_button=ev.select
                 endfor
             endif
             widget_control,ev.top,set_uvalue=list
               endcase
      'Use Intensities':begin
             if (ev.select eq 0) then list.PDFni[list.curPhase]=0 $
             else begin
                  ind=where(list.PDFi[list.i0:list.i1] eq 0,ct)
                  if ct eq 0 then list.PDFni[list.curPhase]=list.PDFn[list.curPhase] $
                  else list.PDFni[list.curPhase]=-1
             endelse
             widget_control,ev.top,set_uvalue=list
             endcase
      'View All':begin
             if ev.select then list.select[list.i0:list.i1]=list.select[list.i0:list.i1] or 1 $
               else list.select[list.i0:list.i1]=list.select[list.i0:list.i1] and 254
               if list.labelview eq 0 then begin
                 for i=0l,list.PDFn[list.curPhase]-1 do begin
                      ID=widget_info(ev.top,FIND_BY_UNAME=stringr(i))
                      widget_control,ID,set_button=ev.select
                 endfor
             endif
             widget_control,ev.top,set_uvalue=list
             endcase
      'Delete':begin
               list.Del[list.curPhase]=ev.select
               ID=widget_info(ev.top,FIND_BY_UNAME='basegroups')
               widget_control,ID,sensitive=ev.select eq 0
               widget_control,ev.top,set_uvalue=list
               endcase
      'Apply':  begin
             ApplyEditPDF,list
             endcase
       'OK': begin
             ApplyEditPDF,list
             WIDGET_CONTROL, ev.top, /DESTROY
             endcase
      'Cancel': WIDGET_CONTROL, ev.top, /DESTROY
      else:     begin
                widget_control,ev.id,get_uvalue=i
                if ev.select then list.select[list.i0+i]=list.select[list.i0+i] or (list.labelview?2:1) $
                   else list.select[list.i0+i]=list.select[list.i0+i] and (list.labelview?253:254)
                widget_control,ev.top,set_uvalue=list
             endcase
    endcase
    endcase
2:    begin
    widget_control,ev.id,get_value=val
    list.PDFScale[list.curPhase]=val/100.
    widget_control,ev.top,set_uvalue=list
    widget_control,list.top,get_uvalue=listmain
    (*listmain.PDFScale)[list.curPhase]=list.PDFScale[list.curPhase]
    if list.topid then CALL_PROCEDURE, list.proc,list.top,listmain $
    else CALL_PROCEDURE, list.proc,listmain
    endcase
3:    begin
    widget_control,ev.id,get_value=val
    percent=0.
    reads,val,percent
    percent=0>percent<100
    ;scale to range 0-100
    PDFi=list.PDFi[list.i0:list.i1]
    PDFi=fix(PDFi)
    ind=where(PDFi ge percent,ct)+list.i0
    list.select[list.i0:list.i1]=list.select[list.i0:list.i1] and (list.labelview?253:254)
    if ct ne 0 then list.select[ind]=list.select[ind] or (list.labelview?2:1)

    c=list.labelview+1
    for i=0l,list.PDFn[list.curPhase]-1 do begin
         ID=widget_info(ev.top,FIND_BY_UNAME=stringr(i))
         widget_control,ID,set_button=(list.select[i+list.i0] and c) eq c
    endfor
    widget_control,ev.top,set_uvalue=list
    endcase
8:    begin
    ; Point to new Phase
    if ev.index eq list.curPhase then return
    list.curPhase=ev.index
    if list.curPhase ne 0 then list.i0=total(list.PDFn[0:list.curPhase-1]) else list.i0=0
    list.i1=total(list.PDFn[0:list.curPhase])-1
    widget_control,ev.top,set_uvalue=list

    ; Display new buttons
    basegroups=EditPDFWindowP(ev,list)

    ; Change button states
    ChangeButtonState,ev,list

    endcase
endcase

end;pro EditPDFHandler
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditPDFWindow,ev,proc,topid=topid,od=od

widget_control,ev.top,get_uvalue=list
if (list.nPDF eq 0) then return

WIDGET_CONTROL, ev.top, sensitive=0

if keyword_set(od) then begin
    PDFx=od.PDFx
    PDFi2=od.PDFi2
    PDFScale=od.PDFScale
    PDFOffset=od.PDFOffset
    od=1b
endif else begin
    od=0b
    PDFx=0
    PDFi2=0
    PDFScale=0
    PDFOffset=0.
endelse

list={top:ev.top,$
       proc:proc,$
       topid:keyword_set(topid),$
       od:od,$
       curPhase:0,$
       labelview:0b,$
       sortond:0b,$
       i0:0,$
       i1:0,$
       PDFx:PDFx,$
       PDFi2:PDFi2,$
       PDFScale:PDFScale,$
       PDFOffset:PDFOffset,$
       nPDF:list.nPDF,$
       del:bytarr(list.nPDF),$
       names:(*list.PDFname),$
       PDFhkl:(*list.PDFhkl),$
       PDFni:(*list.PDFni),$
       PDFn:(*list.PDFn),$
       select:(*list.PDFs),$
       PDFi:(*list.PDFi),$
       PDFd:(*list.PDFd)}
if list.curPhase ne 0 then list.i0=total(list.PDFn[0:list.curPhase-1])
list.i1=total(list.PDFn[0:list.curPhase])-1

base=widget_base(/column,title='Edit PDF',uvalue=list)

basebuttons=widget_base(base,/row)
    baset=widget_base(basebuttons,/row,/nonexclusive)
       vbutton=widget_button(baset,value='View All',uname='View All')
       ubutton=widget_button(baset,value='Use Intensities',uname='Use Intensities')
       xbutton=widget_button(baset,value='View Labels',uname='View Labels')
       wbutton=widget_button(baset,value='Delete',uname='Delete')
       button=widget_button(baset,value='Change Label View',uname='Change Label View')
       button=widget_button(baset,value='Sort on d-spacing',uname='Sort on d-spacing')
    label=widget_label(basebuttons,value='% Cutoff:',uname='% Cutoff')
    text=widget_text(basebuttons,value='0',/editable,uname='% Cutofftxt')
    button=widget_button(basebuttons,value='OK',uname='OK')
    button=widget_button(basebuttons,value='Apply',uname='Apply')
    button=widget_button(basebuttons,value='Cancel',uname='Cancel')

wrapbase=widget_base(base,uname='wrapbase',/column)
baset=widget_base(wrapbase,/row)
    droplist=widget_droplist(baset,value=list.names)
    if keyword_set(od) then begin
        baset=widget_base(wrapbase,/column)
            label=widget_label(baset,value='Scale maximum peak:')
            slider=widget_slider(baset,minimum=-1,maximum=100,/SUPPRESS_VALUE,uname='sclslider',value=list.PDFScale[list.curPhase]*100)
    endif
basegroups=EditPDFWindowP(wrapbase,list)

WIDGET_CONTROL, base, /REALIZE

tmp=where((list.select[list.i0:list.i1] and 1) eq 1,ct)
widget_control,vbutton,set_button=ct gt 0
tmp=where((list.select[list.i0:list.i1] and 2) eq 2,ct)
widget_control,xbutton,set_button=ct gt 0
if list.PDFni[list.curPhase] ne 0 then widget_control,ubutton,set_button=1 $
else widget_control,ubutton,set_button=0
widget_control,droplist,SET_DROPLIST_SELECT=list.curPhase
widget_control,wbutton,set_button=list.Del[list.curPhase]

Xmanager,'EditPDFWindow',base, event_handler='EditPDFHandler',$
    cleanup='CleanUp_EditPDF',GROUP_LEADER=ev.top

end;pro EditPDFWindow
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
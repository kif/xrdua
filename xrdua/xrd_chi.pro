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

pro CleanUp_ProfileRange,ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list2
widget_control,list2.top,sensitive=1
end;pro CleanUp_ProfileRange
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ProfileRange_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list2
widget_control,ev.id,get_value=val

case widget_info(ev.id,/type) of
1 : begin ;button event
    case val of
    'Cancel' :    begin
            ptr_free,list2.ptr
            WIDGET_CONTROL, ev.top, /DESTROY
            endcase
    'OK' :    WIDGET_CONTROL, ev.top, /DESTROY
    else :    begin
            widget_control,ev.id,get_uvalue=uval
            case uval of
            'default':  begin
                        (*list2.ptr)[5]=ev.select
                        ID=widget_info(ev.top,find_by_uname='senbase')
                        widget_control,ID,sensitive=~ev.select
                        endcase
            'mm':         (*list2.ptr)[0]=2
            msymbols('angstroms'):        (*list2.ptr)[0]=0
            'degrees':    (*list2.ptr)[0]=1
            '1/nm':         (*list2.ptr)[0]=3
            endcase
            endelse
    endcase
    endcase
3 : begin ;text event
    widget_control,ev.id,get_uvalue=i
    (*list2.ptr)[i]=float(val)
    case i of
    3:    begin
        (*list2.ptr)[2]=((*list2.ptr)[4]-(*list2.ptr)[1])/((*list2.ptr)[3]-1)
        ID=widget_info(ev.top,find_by_uname='inc')
        widget_control,ID,set_value=stringr((*list2.ptr)[2])
        endcase
    else:begin
        (*list2.ptr)[3]=((*list2.ptr)[4]-(*list2.ptr)[1])/(*list2.ptr)[2]+1
        ID=widget_info(ev.top,find_by_uname='bins')
        widget_control,ID,set_value=stringr(ulong((*list2.ptr)[3]))
        endelse
    endcase

    endcase
endcase

end;pro ProfileRange_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ProfileRangeWindow,top,xtypeconv,x0,xi,xn,xdefault,error=error

error=0b
widget_control,top,sensitive=0

ptr=ptr_new([1,5,0.05,501,30,1])
list2={top:top,ptr:ptr}
topbase=widget_base(/column,title='Profile range:',uvalue=list2)
xs=100

tbase=widget_base(topbase,/row,/nonexclusive)
    buttonp=widget_button(tbase,value='default',uvalue='default',xsize=xs)

tbase=widget_base(topbase,/row,/exclusive)
    buttonp0=widget_button(tbase,value=msymbols('angstroms'),uvalue=msymbols('angstroms'),xsize=xs)
    buttonp1=widget_button(tbase,value='degrees',uvalue='degrees',xsize=xs)
    buttonp2=widget_button(tbase,value='mm',uvalue='mm',xsize=xs)
    buttonp3=widget_button(tbase,value='1/nm',uvalue='1/nm',xsize=xs)

senbase=widget_base(topbase,uname='senbase',sensitive=0,/column)
tbase=widget_base(senbase,/row)
    label=widget_label(tbase,value='Minimum:',xsize=xs)
    text=widget_text(tbase,value=stringr((*ptr)[1]),scr_xsize=xs,/editable,uvalue=1)

tbase=widget_base(senbase,/row)
    label=widget_label(tbase,value='Increment:',xsize=xs)
    text=widget_text(tbase,value=stringr((*ptr)[2]),scr_xsize=xs,/editable,uvalue=2,uname='inc')

tbase=widget_base(senbase,/row)
    label=widget_label(tbase,value='Number of bins:',xsize=xs)
    text=widget_text(tbase,value=stringr(ulong((*ptr)[3])),scr_xsize=xs,/editable,uvalue=3,uname='bins')

tbase=widget_base(senbase,/row)
    label=widget_label(tbase,value='Maximum:',xsize=xs)
    text=widget_text(tbase,value=stringr((*ptr)[4]),scr_xsize=xs,/editable,uvalue=4)

button=widget_button(topbase,value='OK')
button=widget_button(topbase,value='Cancel')

WIDGET_CONTROL, topbase, /REALIZE
widget_control,buttonp,set_button=1
case (*ptr)[0] of
0:    widget_control,buttonp0,set_button=1
1:    widget_control,buttonp1,set_button=1
2:    widget_control,buttonp2,set_button=1
3:    widget_control,buttonp3,set_button=1
endcase

Xmanager,'ProfileRangeWindow',topbase, event_handler='ProfileRange_event',$
cleanup='CleanUp_ProfileRange',GROUP_LEADER=top

if not ptr_valid(ptr) then error=1b else begin
    xtypeconv=fix((*ptr)[0])
    x0=(*ptr)[1]
    xi=(*ptr)[2]
    xn=(*ptr)[3]
    xdefault=(*ptr)[5]
    ptr_free,ptr
endelse
end;pro ProfileRangeWindow
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MASK1D_PROCESS_EVENTS, ev

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
list.xs=ev.x
list.ys=ev.y

WSet, list.drawindex
SetSysVar,list.sysvar

coords = Convert_Coord(ev.x, ev.y, /Device, /To_Data)
GetCRange,ROIx=ROIx,ROIy=ROIy
x1 = ROIx[0] > coords[0] < ROIx[1]
y1 = ROIy[0] > coords[1] < ROIy[1]
PlotS, [x1, x1], ROIy

h=!D.x_size
v=!D.y_size
WSet, list.pixindex
Device, Copy = [0,0, h,v, 0,0,list.drawindex]

Widget_Control, ev.id, Event_Pro='MASK1D_DRAW', $
    draw_motion=1
Widget_Control, ev.top, Set_UValue=list
end;pro MASK1D_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MASK1D_DRAW ,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL']
thisEvent = possibleEventTypes[ev.type]

WSet, list.drawindex
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.pixindex]

x = [list.xs, ev.x]
y = [list.ys, ev.y]
SetSysVar,list.sysvar
coords = Convert_Coord(x, y, /Device, /To_Data)
GetCRange,ROIx=ROIx,ROIy=ROIy
x1 = ROIx[0] > coords[0,0] < ROIx[1]
x2 = ROIx[0] > coords[0,1] < ROIx[1]

IF thisEvent EQ 'UP' THEN BEGIN
    Widget_Control, ev.id, draw_motion=0, $
       Event_Pro='MASK1D_PROCESS_EVENTS'

    if x1 ne x2 then begin
        temp=value_locater((*list.xvalt)[*,list.xtype],[x1,x2])>0
        x1=temp[0]
        x2=temp[1]
    endif else begin
          x1=0
          x2=list.nxval-1
    endelse

    x=(*list.xvalt)[[x1,x2],0]

    if list.masknr eq 0 then begin ; New
        list.mask=ptr_new([1,0,0,x[sort(x)]])
        list.smask=ptr_new([1b])
        list.masknr++
    endif else begin
        ind=where((*list.mask)[2,*] eq 0,count)
        if count eq 0 then begin ; Add
            L=max((*list.mask)[1,*])+1
            (*list.mask)=[[(*list.mask)],[1,L,0,x[sort(x)]]]
            (*list.smask)=[(*list.smask),1b]
            list.masknr++
        endif else begin ; Replace last
            L=(*list.mask)[1,ind[0]]
            (*list.mask)[*,ind[0]]=[1,L,0,x[sort(x)]]
        endelse
    endelse

    RefreshDisplayCHI,list
    Widget_Control, ev.top, Set_UValue=list
    return
endif

list.xd=ev.x
list.yd=ev.y

y2 = ROIy[0] > coords[1,1] < ROIy[1]
WSet, list.drawindex
Arrow, x1, y2, x2, y2, /Data, /Solid, HSize=12
PlotS, [x2, x2], ROIy

Widget_Control, ev.top, Set_UValue=list
end;pro MASK1D_DRAW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MakePeaksout,list,base4
if (list.apeaksout ne 0) then begin
    table=widget_table(base4,all_events=1,$
              X_SCROLL_SIZE=1,Y_SCROLL_SIZE=10,uname='peaksout',uvalue='peaksout',$
             COLUMN_LABELS=['m('+list.xunits[list.xtypeorig]+')','a','s','n'],$
              value=list.peakout[0:3,0:list.apeaksout-1])
    widget_control,table,SET_TABLE_SELECT=[ 0, list.fitselect, 2, list.fitselect ]
    widget_control,table,SET_TABLE_VIEW=[ 0, list.fitselect]
endif
end;pro MakePeaksout
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MakePeaksin,list,base4
if (list.apeaks ne 0) then begin
    if ptr_valid(list.speringpr) then begin
        mat=transpose([[list.peakxen[0:list.apeaks-1]],$
        [list.peakarea[0:list.apeaks-1]],[list.peaksigma[0:list.apeaks-1]],$
        [(*list.speringpr)[list.peakxkan[0:list.apeaks-1]]]])
        lab=['m('+list.xunits[list.xtypeorig]+')','a','s','ring%']
    endif else begin
        mat=transpose([[list.peakxen[0:list.apeaks-1]],$
        [list.peakarea[0:list.apeaks-1]],[list.peaksigma[0:list.apeaks-1]]])
        lab=['m('+list.xunits[list.xtypeorig]+')','a','s']
    endelse
    table=widget_table(base4,all_events=1,$
           X_SCROLL_SIZE=1,Y_SCROLL_SIZE=10,uname='peaksin',uvalue='peaksin',$
          COLUMN_LABELS=lab,$
           value=mat,/editable)
       widget_control,table,SET_TABLE_SELECT=[ 0, list.peakselect, 2, list.peakselect ]
    widget_control,table,SET_TABLE_VIEW=[ 0, list.peakselect]
endif
end;pro MakePeaksin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshPDD_CHI,list

if list.PDFnTotal eq 0 then return

; ----Convert to all possible units----
(*list.PDFx)=CHI_xval((*list.PDFd),list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,0)

; ----Find position intensities----
x=(*list.xvalt)[*,list.xtype]
y=(*list.spet)[*,list.xtype]
ind=sort(x)
x=x[ind]
y=y[ind]
Y2 = SPL_INIT(x, y)
hights = SPL_INTERP(x, y, Y2, (*list.PDFx)[*,list.xtype])

(*list.PDFi2)=hights

end;RefreshPDD_CHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LoadPDD_CHI,ev,prompt,search=search
LoadPDD,ev,prompt,searh=search
widget_control,ev.top,get_uvalue=list
if list.PDFnTotal ne 0 then begin
    RefreshPDD_CHI,list
    widget_control,ev.top,set_uvalue=list
       RefreshDisplayCHI,list
       RefreshWindowCHI,ev
endif
end;pro LoadPDD_CHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_ParamW, ID

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
end;pro CleanUp_ParamW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ParamW_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if widget_info(ev.id,/type) eq 10 then return

widget_control,ev.id,get_uvalue=uval
widget_control,ev.top,get_uvalue=top
widget_control,top,get_uvalue=list
widget_control,ev.id,get_value=val
case uval of
    'OK' :  begin
;            ID=widget_info(ev.top,find_by_uname='a')
;            widget_control,ID,get_value=val
;            list.a=float(val)*!pi/180.
;            ID=widget_info(ev.top,find_by_uname='b')
;            widget_control,ID,get_value=val
;            list.b=float(val)*!pi/180.
            ID=widget_info(ev.top,find_by_uname='scrx')
            widget_control,ID,get_value=val
            list.scrx=float(val)*10^(-6.)
            ID=widget_info(ev.top,find_by_uname='scry')
            widget_control,ID,get_value=val
            list.scry=float(val)*10^(-6.)
            ID=widget_info(ev.top,find_by_uname='dist')
            widget_control,ID,get_value=val
            list.dist=float(val)*0.01
            ID=widget_info(ev.top,find_by_uname='xlambda')
            widget_control,ID,get_value=val
            list.lambda=float(val)
            ID=widget_info(ev.top,find_by_uname='vin')
            widget_control,ID,get_value=val
            list.vin=fix(val)
            ID=widget_info(ev.top,find_by_uname='win')
            widget_control,ID,get_value=val
            list.win=fix(val)
            ID=widget_info(ev.top,find_by_uname='r1')
            widget_control,ID,get_value=val
            list.r1=float(val)
            ID=widget_info(ev.top,find_by_uname='letter')
            widget_control,ID,get_value=val
            list.letter=float(val)
            ;ID=widget_info(ev.top,find_by_uname='kDisp')
            ;widget_control,ID,get_value=val
            ;list.kDisp=fix(val)
            ID=widget_info(ev.top,find_by_uname='DelWidth0')
            widget_control,ID,get_value=val
            list.DelWidth[0]=fix(val)
            ID=widget_info(ev.top,find_by_uname='DelWidth1')
            widget_control,ID,get_value=val
            list.DelWidth[1]=fix(val)
            
            widget_control,top,set_uvalue=list
            RefreshDataCHI,{top:top}
            WIDGET_CONTROL, ev.top, /DESTROY
            return
            endcase
    'vin':list.vin=fix(val)
    'win':list.win=fix(val)
    'r1':list.r1=float(val)
    'letter':list.letter=float(val)
    'kDisp':list.kDisp=fix(val)
    'DelWidth0':begin
         list.DelWidth[0]=fix(val)
         ID=widget_info(ev.top,FIND_BY_UNAME='DelArea')
         WIDGET_CONTROL, ID, /DESTROY
         ID=widget_info(ev.top,FIND_BY_UNAME='drawbase')
         w=2*list.DelWidth+1
         base=widget_draw(ID,xsize=w[0],ysize=w[1],uname='DelArea')
         endcase
    'DelWidth1':begin
         list.DelWidth[1]=fix(val)
         ID=widget_info(ev.top,FIND_BY_UNAME='DelArea')
         WIDGET_CONTROL, ID, /DESTROY
         ID=widget_info(ev.top,FIND_BY_UNAME='drawbase')
         w=2*list.DelWidth+1
         base=widget_draw(ID,xsize=w[0],ysize=w[1],uname='DelArea')
         endcase
    'DelMass':list.DelMass=ev.select
    'UpMain':list.CrossUpdate=ev.select
   'ShortLegend':list.blegendshort=ev.select
   'PDFpsym': list.bPDFpsym=ev.select
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
    'a':   list.a=float(val)*!dpi/180. ;from deg to rad
    'b':   list.b=float(val)*!dpi/180. ;from deg to rad
    'power':begin
            list.ypower=float(val)
            if list.ypower eq 1 or list.ypower eq 0 then begin
                list.ypowerstr=''
            endif else begin
                list.ypowerstr='^'+val
                list.ylog=0
                id=widget_info(top,find_by_uname='Y')
                widget_control,id,set_value='Y: linear'
            endelse
            endcase
endcase
widget_control,top,set_uvalue=list
RefreshDataCHI,{top:top}
end; pro ParamW_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ParamWindow,ev
widget_control,ev.top,get_uvalue=list
widget_control,ev.top,sensitive=0
base=widget_base(title=val,uvalue=ev.top,/row)

wtab=widget_tab(base)

base1=widget_base(wtab,/column,title='Experimental geometry')
xs=150
    baset=widget_base(base1,/row)
    label=widget_label(baset,value='X-ray energy (keV)  :        ',xsize=xs)
    text=widget_text(baset,/editable,value= $
    stringr(list.cconv/list.lambda),uvalue='xen',uname='xen',scr_xsize=xs)

    baset=widget_base(base1,/row)
    label=widget_label(baset,value='X-ray wavelength ('+msymbols('angstroms')+')  :      ',xsize=xs)
    text=widget_text(baset,/editable,value= $
    stringr(list.lambda),uvalue='xlambda',uname='xlambda',scr_xsize=xs)

    baset=widget_base(base1,/row)
    label=widget_label(baset,value='Pixel X-size ('+msymbols('micro')+'m) :          ',xsize=xs)
    text=widget_text(baset,/editable,value= $
          stringr(list.scrx*10^6.),uvalue='scrx',uname='scrx',scr_xsize=xs)

    baset=widget_base(base1,/row)
    label=widget_label(baset,value='Pixel Y-size ('+msymbols('micro')+'m) :          ',xsize=xs)
    text=widget_text(baset,/editable,value= $
          stringr(list.scry*10^6.),uvalue='scry',uname='scry',scr_xsize=xs)
    
    baset=widget_base(base1,/row)
    label=widget_label(baset,value='Distance Sample-Detector(cm):',xsize=xs)
    text=widget_text(baset,/editable,value= $
          stringr(list.dist*10^2.),uvalue='dist',uname='dist',scr_xsize=xs)
          
    baset=widget_base(base1,/row)
    label=widget_label(baset,value='Tilt a('+msymbols('degrees')+'):',xsize=xs)
    text=widget_text(baset,/editable,value= $
    stringr(list.a/!dpi*180.),uvalue='a',uname='a',scr_xsize=xs)

    baset=widget_base(base1,/row)
    label=widget_label(baset,value='Rot b('+msymbols('degrees')+'):',xsize=xs)
    text=widget_text(baset,/editable,value=$
    stringr(list.b/!dpi*180.),uvalue='b',uname='b',scr_xsize=xs)

xs=200
base1=widget_base(wtab,/column,title='Others')
    baset=widget_base(base1,/row)
    label=widget_label(baset,value='Peak search w (= FWHM in bins) :      ',xsize=xs)
    text=widget_text(baset,/editable,value= $
    stringr(list.win),uvalue='win',uname='win',scr_xsize=xs)
    
    baset=widget_base(base1,/row)
    label=widget_label(baset,value='Peak search v (in [FWHM/3,FWHM/2]) :      ',xsize=xs)
    text=widget_text(baset,/editable,value= $
    stringr(list.vin),uvalue='vin',uname='vin',scr_xsize=xs)

    baset=widget_base(base1,/row)
    label=widget_label(baset,value='Peak search r (threshold r.stdev) :      ',xsize=xs)
    text=widget_text(baset,/editable,value= $
    stringr(list.r1),uvalue='r1',uname='r1',scr_xsize=xs)

    baset=widget_base(base1,/row)
    label=widget_label(baset,value='Label size  :      ',xsize=xs)
    text=widget_text(baset,/editable,value= $
    stringr(list.letter),uvalue='letter',uname='letter',scr_xsize=xs)

    baset=widget_base(base1,/row)
    label=widget_label(baset,value='Y-axis power  :    ',xsize=xs)
    text=widget_text(baset,/editable,value= $
    strmid(list.ypowerstr,1),uvalue='power',uname='power',scr_xsize=xs)
    
    ;baset=widget_base(base1,/row)
    ;label=widget_label(baset,value='Peak display width  :      ',xsize=xs)
    ;text=widget_text(baset,/editable,value= $
    ;stringr(list.kDisp),uvalue='kDisp',uname='kDisp',scr_xsize=xs)

    baset=widget_base(base1,/row,uname='drawbase')
    label=widget_label(baset,value='Detect labels width  :      ',xsize=xs)
    text=widget_text(baset,/editable,value= $
    stringr(list.DelWidth[0]),uvalue='DelWidth0',uname='DelWidth0',scr_xsize=xs/2)
    text=widget_text(baset,/editable,value= $
    stringr(list.DelWidth[1]),uvalue='DelWidth1',uname='DelWidth1',scr_xsize=xs/2)
         w=2*list.DelWidth+1
         dr=widget_draw(baset,xsize=w[0],ysize=w[1],uname='DelArea')

    baset=widget_base(base1,/column,/nonexclusive)
    buttonmass=widget_button(baset,value='Delete multiple labels without prompt',uvalue='DelMass')
    buttonupmain=widget_button(baset,value='Update Main window',uvalue='UpMain')
    buttonlegend=widget_button(baset,value='Short Legend',uvalue='ShortLegend')
    buttonPDFpsym=widget_button(baset,value='PDF symbols',uvalue='PDFpsym')

button=widget_button(base,value='OK',uvalue='OK')
WIDGET_CONTROL, base, /REALIZE

widget_control,buttonmass,set_button=list.DelMass
widget_control,buttonupmain,set_button=list.CrossUpdate
widget_control,buttonlegend,set_button=list.blegendshort
widget_control,buttonPDFpsym,set_button=list.bPDFpsym
Xmanager,'ParamWindow',base, event_handler='ParamW_event',cleanup='CleanUp_ParamW',GROUP_LEADER=ev.top
end;pro ParamWindow
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro peaks,list,out=out,apeaks=apeaks,x=x,y=y,offx=offx,offy=offy

if not(keyword_set(x)) then x=(*list.xvalt)[*,list.xtypeorig]
if not(keyword_set(y)) then y=(*list.spet)[*,list.xtypeorig]
if not(keyword_set(offx)) then offx=0
if not(keyword_set(offy)) then offy=0
;A,m,s,n
out=estspe(y,x,list.r1,list.vin,list.win,/details)

s=size(out)
if s[1] ne 8 then begin
    printw,list.OutID,'No peaks found'
    apeaks=0
    return
endif
if s[0] eq 1 then apeaks=1 else apeaks=s[2]

list.apeaks=apeaks
list.peakarea=list.peakarea*0
list.peakxen=list.peakxen*0
list.peaksigma=list.peaksigma*0
list.peakxkan=list.peakxkan*0
list.peakbegin=list.peakbegin*0
list.peakend=list.peakend*0
list.peaky=list.peaky*0

list.peakarea[0:list.apeaks-1]=out[0,*]
list.peakxen[0:list.apeaks-1]=out[1,*]
list.peaksigma[0:list.apeaks-1]=out[2,*]
list.peakxkan[0:list.apeaks-1]=out[4,*]+offx
list.peakbegin[0:list.apeaks-1]=out[5,*]+offx
list.peakend[0:list.apeaks-1]=out[6,*]+offx
list.peaky[0:list.apeaks-1]=out[7,*]

out=Shrinkarray(out,[4,5,6,7])
end;pro peaks
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FitScan,list

;----Clean output----
list.peakout=list.peakout*0
list.peakouts=list.peakouts*0
list.apeaksout=0

;----Set ROI----
x=(*list.xvalt)[list.backranind[0]:list.backranind[1],list.xtypeorig]
y=(*list.spet)[list.backranind[0]:list.backranind[1],list.xtypeorig]

out=fitod(x,y,list,interactive=list.OutID)

end;pro FitScan
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DRAG_DROP_EVENTS , ev

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
IF thisEvent EQ 'DOWN' THEN BEGIN
    WSet, list.drawindex
    Device, Copy = [0,0, !D.x_size,!D.y_size, 0,0,list.pixindex]
    widget_control,ev.id,draw_motion=0, $
              Event_Pro='VERPM_PROCESS_EVENTS'

    SetSysVar,list.sysvar
    coords = Convert_Coord(ev.x, ev.y, /Device, /To_Data)
    GetCRange,ROIx=ROIx,ROIy=ROIy
    x = ROIx[0] > coords[0] < ROIx[1]
    y = ROIy[0] > coords[1] < ROIy[1]
    ysave = y
    if list.ypower ne 0 and list.ypower ne 1 then ysave^=1/list.ypower

    WSet, list.drawindex
    xyouts,x,y,list.dragdrop,color=[255, 255, 255],charsize=list.letter
    WSet, list.pixindex
    xyouts,x,y,list.dragdrop,color=[255, 255, 255],charsize=list.letter
    if list.alabels eq 0 then begin
        ptr_free,list.labels
        ptr_free,list.xlabels
        ptr_free,list.ylabels
        list.labels=ptr_new(list.dragdrop)
        list.ylabels=ptr_new(ysave)
        list.xlabels=ptr_new(CHI_xval(x,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtype))
    endif else begin
        (*list.labels)=[(*list.labels),list.dragdrop]
        (*list.ylabels)=[(*list.ylabels),ysave]
        (*list.xlabels)=[(*list.xlabels),CHI_xval(x,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtype)]
    endelse
    list.alabels++
    list.dragdrop=''
    Widget_Control, ev.top, Set_UValue=list
endif
WSet, list.drawindex
Device, Copy = [0,0, !D.x_size,!D.y_size, 0,0,list.pixindex]

SetSysVar,list.sysvar
coords = Convert_Coord(ev.x, ev.y, /Device, /To_Data)
GetCRange,ROIx=ROIx,ROIy=ROIy
x = ROIx[0] > coords[0] < ROIx[1]
y = ROIy[0] > coords[1] < ROIy[1]
WSet, list.drawindex
xyouts,x,y,list.dragdrop,color=[255, 255, 255],charsize=list.letter
end; pro DRAG_DROP_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DELM_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]

case thisEvent of
'UP' :  begin
       widget_control,ev.id,draw_motion=0
       return
       endcase
'DOWN' :begin
       widget_control,ev.id,draw_motion=1
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
if list.alabels eq 0 then return

; ----Find in label list which label is clicked----
x1=ev.x-list.DelWidth[0]
x2=ev.x+list.DelWidth[0]
y1=ev.y-list.DelWidth[1]
y2=ev.y+list.DelWidth[1]
SetSysVar,list.sysvar
ylabels=(*list.ylabels)[0:list.alabels-1]
if list.ypower ne 0 and list.ypower ne 1 then ylabels^=list.ypower
coords = Convert_Coord((*list.xlabels)[0:list.alabels-1,list.xtype], ylabels, /Data, /To_Device)
indP=where( (coords[0,*] gt x1 and coords[0,*] lt x2 and $
             coords[1,*] gt y1 and coords[1,*] lt y2),count)

if count eq 0 then return
labels=(*list.labels)
xlabels=(*list.xlabels)[*,list.xtype]
ylabels=(*list.ylabels)
if count ne 1 then begin
    str=[reform(labels,1,list.alabels),stringr([reform(xlabels,1,list.alabels),$
                reform(ylabels,1,list.alabels)])]
    str=str[0,*]+':'+'['+str[1,*]+';'+str[2,*]+']'
    if (list.DelMass eq 0) then begin
        widget_control,ev.id,draw_motion=0
        sel=ChooseFromList(str,title='Select labels to delete',$
              nsel=nsel,top=ev.top)
        if nsel eq 0 then return
        indP=indP[sel]
    endif
endif
xlabels=(*list.xlabels)

;----Delete Point----
labels=ShrinkArray(labels,indP)
xlabels=ShrinkArray(xlabels,indP)
ylabels=ShrinkArray(ylabels,indP)
list.alabels=list.alabels-n_elements(indP)
ptr_free,list.labels
ptr_free,list.xlabels
ptr_free,list.ylabels
if list.alabels ne 0 then begin
    list.labels=ptr_new(labels)
    list.xlabels=ptr_new(xlabels)
    list.ylabels=ptr_new(ylabels)
endif
Widget_Control, ev.top, Set_UValue=list
RefreshDisplayCHI,list
end ;pro DELM_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro VERPM_PROCESS_EVENTS,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]

case thisEvent of
'UP' :  begin
       widget_control,ev.id,draw_motion=0
       return
       endcase
'DOWN' :begin
       widget_control,ev.id,draw_motion=1
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
if list.alabels eq 0 then return

; ----Find in label list which label is clicked----
x1=ev.x-list.DelWidth[0]
x2=ev.x+list.DelWidth[0]
y1=ev.y-list.DelWidth[1]
y2=ev.y+list.DelWidth[1]
SetSysVar,list.sysvar
ylabels=(*list.ylabels)[0:list.alabels-1]
if list.ypower ne 0 and list.ypower ne 1 then ylabels^=list.ypower
coords = Convert_Coord((*list.xlabels)[0:list.alabels-1,list.xtype], ylabels, /Data, /To_Device)
indP=where( (coords[0,*] gt x1 and coords[0,*] lt x2 and $
             coords[1,*] gt y1 and coords[1,*] lt y2),count)

if count eq 0 then return
labels=(*list.labels)
xlabels=(*list.xlabels)[*,list.xtype]
ylabels=(*list.ylabels)
if count ne 1 then begin
    str=[reform(labels,1,list.alabels),stringr([reform(xlabels,1,list.alabels),$
                reform(ylabels,1,list.alabels)])]
    str=str[0,*]+':'+'['+str[1,*]+';'+str[2,*]+']'
    if (list.DelMass eq 0) then begin
        widget_control,ev.id,draw_motion=0
        sel=ChooseFromList(str,title='Select labels to delete',$
              nsel=nsel,top=ev.top)
        if nsel eq 0 then return
        indP=indP[sel]
    endif
endif
xlabels=(*list.xlabels)

;----Delete Point----
list.dragdrop=labels[indP]
labels=ShrinkArray(labels,indP)
xlabels=ShrinkArray(xlabels,indP)
ylabels=ShrinkArray(ylabels,indP)
list.alabels=list.alabels-n_elements(indP)
ptr_free,list.labels
ptr_free,list.xlabels
ptr_free,list.ylabels
if list.alabels ne 0 then begin
    list.labels=ptr_new(labels)
    list.xlabels=ptr_new(xlabels)
    list.ylabels=ptr_new(ylabels)
endif

Widget_Control, ev.top, Set_UValue=list
RefreshDisplayCHI,list
widget_control,list.draw,Event_Pro='DRAG_DROP_EVENTS',draw_motion=1
end ;pro VERPM_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro modeCHI, ev,val
;widget_control,ev.top,update=0
widget_control,ev.top,get_uvalue=list
strmat=['Data range','Zoom','Mask Setup','Select Peak','Select Any','Delete label','Move label','Refine','Background','PDF markers','Edit Peaks','Move Legend']
ind=where(strmat ne val,count)
if count ne 0 then strmat=strmat[ind]
widget_control,list.draw,draw_button=0,draw_motion=0
bool=0B
bsetsen=1b

case val of
'Move Legend':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
             widget_control,list.draw,Event_Pro='CHANGE_LEGENDCHI_EVENT',draw_button=1,draw_motion=1
            endcase
    'Data range':     begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
             widget_control,list.draw,Event_Pro='ROI_PROCESS_EVENTS_CHI',draw_button=1
             endcase
    'Zoom': begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
             widget_control,list.draw,Event_Pro='ZOOM_PROCESS_EVENTS_CHI',draw_button=1
               endcase
'Mask Setup':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
             widget_control,list.draw,Event_Pro='MASK1D_PROCESS_EVENTS',draw_button=1

            ID2=widget_info(ev.top,FIND_BY_UNAME='base5')
            base=widget_base(ID2,/row,uname='Mbase',uvalue=fltarr(5))
            base1=widget_base(base,/column)
            label1=widget_label(base1,value='Define area:')
            button4=widget_button(base1,value='Set',uname='Set')
            button6=widget_button(base1,value='Edit Mask',uname='Edit Mask')

            endcase
'Select Peak':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
             widget_control,list.draw,Event_Pro='Select_Peak_Event',draw_button=1
               endcase
'Select Any':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
             widget_control,list.draw,Event_Pro='Select_Any_Event',draw_button=1
               endcase
'Delete label':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
             widget_control,list.draw,Event_Pro='DELM_PROCESS_EVENTS',draw_button=1
               endcase
'Move label':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
             widget_control,list.draw,Event_Pro='VERPM_PROCESS_EVENTS',draw_button=1
               endcase
'Refine':    begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
            ID2=widget_info(ev.top,FIND_BY_UNAME='base5')
            base=widget_base(ID2,/column,uname='Mbase')
               baset=widget_base(base,/column,/exclusive)
               b=lonarr(6)
            b[0]=widget_button(baset,value='MLR(G) substr' ,uname='Fit0',xsize=br1)
            b[1]=widget_button(baset,value='NLLS(G) ortpolc' ,uname='Fit1',xsize=br1)
            b[2]=widget_button(baset,value='NLLS(G) substr' ,uname='Fit2',xsize=br1)
            b[3]=widget_button(baset,value='NLLS(PVn) substr' ,uname='Fit3',xsize=br1)
            b[5]=widget_button(baset,value='NLLS(PVnarr) substr' ,uname='Fit5')
            if list.filorder then button12=widget_button(base,value='Peaksearch after substr',uname='filorder') $
            else button12=widget_button(base,value='Peaksearch before substr',uname='filorder')
            button12=widget_button(base,value='Fit Spectrum',uname='Fit Spectrum')
            button12=widget_button(base,value='Reset Fit',uname='Reset Fit')
            label6=widget_label(base,value='',uname='Chisq',xsize=br1)
            label6=widget_label(base,value='',uname='Iter',xsize=br1)
            baset=widget_base(base,/column,uname='FitSearch')
                button=widget_button(baset,value='>',uname='>',uvalue=0)
                button=widget_button(baset,value='<',uname='<',uvalue=0)
                 table=widget_table(baset,ROW_LABELS=['m','a','s','n'],EDITABLE=0,$
                          xsize=2,ysize=4,X_SCROLL_SIZE=1,Y_SCROLL_SIZE=4,uname='table',$
                          COLUMN_LABELS=['Val('+list.xunits[list.xtypeorig]+')','SD'])
            widget_control,b[list.fitmethod],set_button=1
            bool=1B
            endcase
'Background':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
            ID2=widget_info(ev.top,FIND_BY_UNAME='base5')
            base4=widget_base(ID2,/column,uname='Mbase')

            tab=widget_tab(base4,uvalue='1Dback')
            baset=widget_base(tab,title='none',/column)
            baset=widget_base(tab,title='ortpol',/column)
                 label5=widget_label(baset,value='Degree Orthogonal Pol:')
                 text2= widget_text(baset,/EDITABLE,value=stringr((*(*list.ptrModel).ptrFit).backparamIO[0]),uname='degree',uvalue='degree',xsize=br)
             baset=widget_base(tab,title='strip',/column)
                 label4=widget_label(baset,value='Peak Stripping:')
                 button34=widget_button(baset,value='Strip',uname='Strip',xsize=br)
                 label5=widget_label(baset,value='Width:')
                 text4= widget_text(baset,/EDITABLE,value=stringr((*(*list.ptrModel).ptrFit).backparamIO[2]),uname='Width',uvalue='Width',xsize=br)
                 label5=widget_label(baset,value='Niter:')
                 text5= widget_text(baset,/EDITABLE,value=stringr((*(*list.ptrModel).ptrFit).backparamIO[1]),uname='Niter',uvalue='Niter',xsize=br)
            tabc=list.backmod
             if tabc le 0 then tabc++
             widget_control,tab,set_tab_current=tabc
             button32=widget_button(base4,value='Subtract',uname='Subtract',xsize=br1,sensitive=list.backcor eq 0)
             button31=widget_button(base4,value='Restore',uname='Restore',xsize=br1,sensitive=list.backcor)
            endcase
'PDF markers':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
            ID2=widget_info(ev.top,FIND_BY_UNAME='base5')
            base=widget_base(ID2,/column,uname='Mbase')

            label=widget_text(base,value=(ptr_valid(list.PDFname)?(*list.PDFname):''),uname='PDFname',YSIZE=list.nPDF<20)
            label=widget_button(base,value='PDF >>')
            label=widget_button(base,value='PDF <<')
            tbase=widget_base(base,/nonexclusive)
                button=widget_button(tbase,value='Keep Last Added',uname='Keep PDF')
                widget_control,button,set_button=~list.ReplaceLastPDF
            tab=widget_tab(base,uvalue='PDFsearch')
                tbase=widget_base(tab,title='Search string',/column)
                    label=widget_text(tbase,value=list.PDFSearch.PDFSubstr,uvalue='PDFSubstr',uname='PDFSubstr',/editable)
                tbase=widget_base(tab,title='Search d-spacing',/column)
                    label=widget_label(tbase,value='From')
                    label=widget_text(tbase,value=string(list.PDFSearch.PDFdi),uvalue='PDFdi',uname='PDFdi',/editable)
                    label=widget_label(tbase,value='To')
                    label=widget_text(tbase,value=string(list.PDFSearch.PDFde),uvalue='PDFde',uname='PDFde',/editable)
            endcase
'Edit Peaks':begin
            ID2=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            if ID2 ne 0 then widget_control,ID2,/destroy
            ID2=widget_info(ev.top,FIND_BY_UNAME='base5')
            base=widget_base(ID2,/column,uname='Mbase')

            base5=widget_base(base,/column,/exclusive)
                button=widget_button(base5,value='Initial')
                widget_control,button,set_button=1
                button=widget_button(base5,value='Refined')
            base4=widget_base(base,/column,uname='RepBaseI')
            button=widget_button(base4,value='Recall',uname='Recall')
            button=widget_button(base4,value='Add iPeak',uname='Add iPeak')
            button=widget_button(base4,value='Delete iPeak',uname='Delete iPeak')
            button=widget_button(base4,value='Select iPeak',uname='Select iPeak')
            button=widget_button(base4,value='Delete iSelect')
            button=widget_button(base4,value='Peak-10',uname='Peak-10')
            button=widget_button(base4,value='Peak-1',uname='Peak-1')
            button=widget_button(base4,value='Peak+1',uname='Peak+1')
            button=widget_button(base4,value='Peak+10',uname='Peak+10')
            MakePeaksin,list,base4
            endcase
'Initial':    begin
            if ev.select eq 0 then return
            ID2=widget_info(ev.top,FIND_BY_UNAME='RepBaseR')
            if widget_info(ID2,/VALID_ID) then widget_control,ID2,/destroy
            base=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            base4=widget_base(base,/column,uname='RepBaseI')
            button=widget_button(base4,value='Recall',uname='Recall')
               button=widget_button(base4,value='Add iPeak',uname='Add iPeak')
            button=widget_button(base4,value='Delete iPeak',uname='Delete iPeak')
            button=widget_button(base4,value='Select iPeak',uname='Select iPeak')
            button=widget_button(base4,value='Delete iSelect')
            button=widget_button(base4,value='Peak-10',uname='Peak-10')
            button=widget_button(base4,value='Peak-1',uname='Peak-1')
            button=widget_button(base4,value='Peak+1',uname='Peak+1')
            button=widget_button(base4,value='Peak+10',uname='Peak+10')
            MakePeaksin,list,base4
            return
            endcase
'Refined':    begin
            if ev.select eq 0 then return
            ID2=widget_info(ev.top,FIND_BY_UNAME='RepBaseI')
            widget_control,ID2,/destroy
            base=widget_info(ev.top,FIND_BY_UNAME='Mbase')
            base4=widget_base(base,/column,uname='RepBaseR')
            button=widget_button(base4,value='Delete rSelect')
            button=widget_button(base4,value='>',uname='>',uvalue=1)
            button=widget_button(base4,value='<',uname='<',uvalue=1)
            MakePeaksout,list,base4
            return
            endcase
'Delete iPeak':begin
            strmat=['Add iPeak','Select iPeak']
            widget_control,list.draw,Event_Pro='Delete_Peak_Event',draw_button=1
            endcase
'Delete iSelect':begin
            if list.apeaks eq 0 then return
            arr=[[list.peaky[0:list.apeaks-1]],$
                [list.peakxen[0:list.apeaks-1]],$
                [list.peakxkan[0:list.apeaks-1]],$
                [list.peaksigma[0:list.apeaks-1]],$
                [list.peakarea[0:list.apeaks-1]],$
                [list.peakbegin[0:list.apeaks-1]],$
                [list.peakend[0:list.apeaks-1]]]
            arr2=ShrinkArray(arr,list.peakselect)

            list.apeaks=list.apeaks-1
            peaky=list.peaky*0
            list.peaky=peaky
            list.peakxen=peaky
            list.peakxkan=peaky
            list.peaksigma=peaky
            list.peakarea=peaky
            list.peakbegin=peaky
            list.peakend=peaky

            if list.apeaks ne 0 then begin
            list.peaky=arr2[*,0]
            list.peakxen=arr2[*,1]
            list.peakxkan=arr2[*,2]
            list.peaksigma=arr2[*,3]
            list.peakarea=arr2[*,4]
            list.peakbegin=arr2[*,5]
            list.peakend=arr2[*,6]
            endif

            list.peakselect=0>(list.peakselect<(list.apeaks-1))
            RefreshDisplayCHI,list
            widget_control,ev.top,Set_UValue=list
            return
            endcase
'Add iPeak':begin
            strmat=['Delete iPeak','Select iPeak']
            widget_control,list.draw,Event_Pro='Add_Peak_Event',draw_button=1
            endcase
'Select iPeak':begin
            strmat=['Delete iPeak','Add iPeak']
             widget_control,list.draw,Event_Pro='Select_Peak_Event',draw_button=1
               endcase
'Delete rSelect':begin
            if (list.apeaksout ne 0) then begin
            arr=[list.peakout[*,0:list.apeaksout-1],list.peakouts[*,0:list.apeaksout-1]]
            arr2=Shrinkarray(arr,list.fitselect,/ROW)
            list.apeaksout=list.apeaksout-1
            list.peakout=list.peakout*0
            list.peakouts=list.peakout
            if (list.apeaksout ne 0) then begin
                list.peakout=arr2[0:3,*]
                list.peakouts=arr2[4:7,*]
            endif
            list.fitselect=0>(list.fitselect<(list.apeaksout-1))
            RefreshDisplayCHI,list
            widget_control,ev.top,Set_UValue=list
            endif
            return
            endcase
    else:    bsetsen=0b
endcase
widget_control,ev.top,Set_UValue=list
setsens,ev.top,strmat
if bsetsen then widget_control,ev.id,sensitive=0
;widget_control,ev.top,update=1
if bool then RefreshWindowCHI,ev
end ;pro modeCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CHANGE_LEGENDCHI_EVENT,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
if list.PDFnTotal eq 0 then return
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]

WSet, list.drawindex
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.pixindex]

IF thisEvent EQ 'UP' THEN BEGIN
    temp=convert_coord(ev.x,ev.y,/device,/to_normal)
    list.xpos=temp[0]
    list.ypos=temp[1]
    RefreshDisplayCHI,list
    modeCHI,ev,'dummy'
    Widget_Control, ev.top, Set_UValue=list
    return
ENDIF

WSet, list.drawindex
xpos=ev.x
ypos=ev.y
incpos=list.incpos*list.letter*!d.y_size
for i=0l,list.nPDF-1 do begin
    xyouts,xpos,i*incpos+ypos,'_____ '+(*list.PDFname)[i],/device,col=i+1,charsize=list.letter
endfor

end;pro CHANGE_LEGENDCHI_EVENT
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotCHI, list,PS=PS,CGM=CGM
; TODO: add sqrt-scale for Y

if keyword_set(PS) then LoadctRev,39,/silent else LoadctRev,-39,/silent
; TODO: replace all '(o)plot' without color with color=255 and remove this line

; 39:
; white:         0
; dark purple:   1
; light purple:  2
; dark blue a:   3
; dark blue b    4
; medium blue    5
; dark cyan      6
; light cyan     7
; green-blue     8
; medium green   9
; light green   10
; yellow-green  11
; yellow        12
; beige         13
; orange        14
; red           15
;col0=16*0;background
;col1=16*0;spec
col2=16*5;yback
col3=16*1;fit
col4=16*7;peaks
;col5=16*0;roi
col6=16*12;PDD
col7=16*13;mask
col8=16*7;mask

; ----YRange----
ny=n_elements((*list.spet)[*,list.xtype])
yrange=0>list.yrange<(ny-1)
yrange=(*list.spet)[list.yrange,list.xtype]
if yrange[0] eq yrange[1] then begin
    xrange=0>list.xrange<(ny-1)
    if xrange[0] eq xrange[1] then xrange=[0,ny-1] $
    else xrange=xrange[sort(xrange)]
    
    tmp1=max((*list.spet)[xrange[0]:xrange[1],list.xtype],min=tmp2)
    yrange=[tmp2,tmp1]
    yrange[1]+=0.05*(yrange[1]-yrange[0])
    yrange[0]-=0.01*(yrange[1]-yrange[0])
    yrange[0]<=0
endif else begin
    d=yrange[1]-yrange[0]
    yrange[0]-=0.01*d
    yrange[1]+=0.05*d
endelse

; ----plot spectrum----
if keyword_set(CGM) then position=[0.18,0.1,0.95,0.65]
if (list.ylog eq 1) then $
    if yrange[0] ne yrange[1] then begin
        ind=where((*list.spet)[*,list.xtype] ne 0,ct)
        if ct ne 0 then m=min((*list.spet)[ind,list.xtype]) else m=1
        yrange[0]>=m
    endif

if list.byerror then yerror=((*list.speerrort)[*,list.xtype])

powplot,(*list.xvalt)[*,list.xtype],((*list.spet)[*,list.xtype]),XRange=(*list.xvalt)[list.xrange,list.xtype], YRange=yrange,ylog=list.ylog,$
       xtitle=list.xtitle,ytitle=list.ytitle+' '+list.ypowerstr,xstyle=1,ystyle=1,psym=list.psym,position=position,charsize=list.letter,xlog=list.xlog,$
       SYMSIZE=0.5*list.letter,pow=list.ypower,yerror=yerror

GetCRange,ROIx=ROIx,ROIy=ROIy

; ----plot background----
if ptr_valid(list.yback) then begin
temp=where((*list.yback) ne 0,ct)
if (ct ne 0) and (list.sBackground) then begin
    powoplot,(*list.xvalt)[list.backranind[0]:list.backranind[1],list.xtype],$
       (*list.yback)[list.backranind[0]:list.backranind[1]],NOCLIP = 0,linestyle=2,color=col2,pow=list.ypower
endif
endif

; ----plot model----
FitModelplot,list,[col4,col3,col3]

; ----plot yfit----
if ptr_valid(list.yfit) then begin
temp=where((*list.yfit) ne 0,ct)
if (ct ne 0) and (list.fitresolve ne 0) then begin
    if (list.fitresolve eq 2)  and (list.xtype eq list.xtypeorig) then begin
       x=(*list.xvalt)[list.backranind[0]:list.backranind[1],list.xtype]
       yb=(*list.yback)[list.backranind[0]:list.backranind[1]]
       m=list.peakout[0,0:list.apeaksout-1]
       s=list.peakout[2,0:list.apeaksout-1]
       a=list.peakout[1,0:list.apeaksout-1]
       n=list.peakout[3,0:list.apeaksout-1]
       for i=0l,list.apeaksout-1 do begin
         if i eq list.fitselect then thick=3 else thick=1
         powoplot,x,(yb+PVconstr(x,a[i],m[i],s[i],n[i]))>1,NOCLIP = 0,color=col3,thick=thick,pow=list.ypower
       endfor
    endif else begin
       powoplot,(*list.xvalt)[list.backranind[0]:list.backranind[1],list.xtype],$
         (*list.yfit)[list.backranind[0]:list.backranind[1]]>1,NOCLIP = 0,color=col3,pow=list.ypower
    endelse
endif
endif

; ----plot Peaks----
if (list.peaksflag ne 0) and (list.apeaks ne 0) then begin
    for i=0l,list.apeaks-1 do begin
       if (i eq list.peakselect) then begin
               if (list.kDisp ne 0) then begin
                   m=list.peakxen[i]
                s=list.peaksigma[i]
                   a=list.peakarea[i]
                   temp=min(abs((*list.xvalt)[*,list.xtype]-(m-list.kDisp*s)),index0)
                   temp=min(abs((*list.xvalt)[*,list.xtype]-(m+list.kDisp*s)),index1)
                   nx=list.nxval-1
                   index0=0>index0<nx
                   index1=0>index1<nx
                   if index0 ne index1 then begin
                       index=[index0<index1,index0>index1]
                    x=(*list.xvalt)[index[0]:index[1],list.xtype]
                    yb=(*list.yback)[index[0]:index[1]]
                    powoplot,x,(yb+PVconstr(x,a,m,s,0))>1,NOCLIP = 0,color=col4,pow=list.ypower
                endif
            endif
            line=0
       endif else line=1

       powplots,(*list.xvalt)[[list.peakxkan[i],list.peakxkan[i]],list.xtype],[list.ylog,list.peaky[i]],$
       linestyle=line,NOCLIP = 0,color=col4,pow=list.ypower

       x1=(*list.xvalt)[list.peakbegin[list.peakselect],list.xtype]
       y1=(*list.spet)[list.peakbegin[list.peakselect],list.xtype]
       x2=(*list.xvalt)[list.peakend[list.peakselect],list.xtype]
       y2=(*list.spet)[list.peakend[list.peakselect],list.xtype]
       powoplot,[x1,x2],[y1,y2],psym=2,NOCLIP = 0,color=col4,pow=list.ypower

       if list.sringpr and ptr_valid(list.speringpr) then xyouts,(*list.xvalt)[list.peakxkan[i],list.xtype],$
            list.peaky[i],stringr(fix((*list.speringpr)[list.peakxkan[i]])),charsize=list.letter,NOCLIP = 0
    endfor
endif

; ----plot marker----
if (list.peakselectflag ne 0) then begin
    powplots,list.anyselect[[list.xtype,list.xtype]],[list.ylog,yrange[1]],$
       linestyle=0,NOCLIP = 0,color=col3,pow=list.ypower
    if list.sringpr and ptr_valid(list.speringpr) then powxyouts,list.anyselect[list.xtype],$
       (*list.spet)[list.anyselect[4]<(ny-1),list.xtype],stringr(fix((*list.speringpr)[list.anyselect[4]])),charsize=list.letter,NOCLIP = 0,pow=list.ypower
endif

; ----plot ROI-lines----
;if (list.backranind[0] ne list.backranind[1]) then begin
if list.plotROI then begin
    ROIx2=(*list.xvalt)[list.backranind,list.xtype]
    PlotS, [ROIx2[0], ROIx2[0]], ROIy,linestyle=2,NOCLIP = 0;,pow=list.ypower
    PlotS, [ROIx2[1], ROIx2[1]], ROIy,linestyle=2,NOCLIP = 0;,pow=list.ypower
endif
;endif

; ----plot labels----
if (list.alabels ne 0) and (list.sLabels ne 0) then for i=0l,list.alabels-1 do powxyouts,(*list.xlabels)[i,list.xtype],$
       (*list.ylabels)[i],(*list.labels)[i],charsize=list.letter,NOCLIP = 0,pow=list.ypower;,color=col2

;----Display PDF data----
if (list.PDFnTotal ne 0) then begin
    LoadctRev,ntableadd(1),/silent
    off=0
    mi=min((*list.xvalt)[*,list.xtype])
    ma=max((*list.xvalt)[*,list.xtype])
    labellinex=!d.y_size*list.letter*0.015
    labelliney=!d.y_size*list.letter*0.010
    xpos=list.xpos
    ypos=list.ypos
    incpos=list.incpos*list.letter
    psymoffset=0.02*(ROIy[1]-ROIy[0])

    for i=0l,list.nPDF-1 do begin
        ; Find PDF markers within the data range
        x=(*list.PDFx)[off:off+(*list.PDFn)[i]-1,list.xtype]
        code=(*list.PDFs)[off:off+(*list.PDFn)[i]-1]
        ind=where((x ge mi)and(x le ma)and(code ne 0),ct)
        col=i+1
        psym=(i+4)<7

        if ct ne 0 then begin
            code=code[ind]
            x=x[ind]
            ind+=off

            ; Use relative intensities or not
            bpow=list.ypower ne 0 and list.ypower ne 1
            PDFOffset=ROIy[0]
            if ((*list.PDFni)[i] ne 0) then begin
                scl=(*list.PDFScale)[i]
                PDFi=(*list.PDFi)[ind]
                PDFi/=max(PDFi)
                if bpow then PDFi^=list.ypower
                if (scl ge 0) and (scl le 1) then begin
                    PDFi*=(ROIy[1]-PDFOffset)*scl
                    PDFi+=PDFOffset
                endif else begin

                    PDFi2=(*list.PDFi2)[ind]
                    if bpow then PDFi2^=list.ypower
                    PDFi2-=PDFOffset  ; Pattern values at Bragg lines - offset
                    m=max(PDFi,indm)    ; most intens line
                    m=min(PDFi2[indm])    ; m = Pattern value of the most intense line
                        
                    PDFi*=m
                    PDFi+=PDFOffset
                endelse
            endif else begin
                PDFi=(*list.PDFi2)[ind]
                if bpow then PDFi^=list.ypower
            endelse
            
            ; Labels
            hkl=fix((*list.PDFhkl)[*,ind])
            mihkl=min(hkl)
            n=strlen(stringr(max(hkl)))>(strlen(stringr(mihkl)))
            label=string(hkl,format='(3(I'+stringr(n)+'))')

            for j=0l,ct-1 do begin
                if (code[j] and 1) eq 1 then begin
                    PlotS,[x[j],x[j]],[PDFOffset,PDFi[j]],color=col,NOCLIP = 0,/data;,pow=list.ypower
                    if list.bPDFpsym then PlotS,x[j],PDFi[j]+psymoffset,psym=psym,NOCLIP = 0, color=col,/data,symsize=0.7;,pow=list.ypower
                endif
                if (code[j] and 2) eq 2 then begin
                    cout=CONVERT_COORD(x[j],PDFi[j],/DATA,/TO_DEVICE)
                    PlotS,[cout[0],cout[0]+labellinex],[cout[1],cout[1]+labelliney],color=col,/device,NOCLIP = 0
                    xyouts,cout[0]+labellinex,cout[1]+labelliney,label[j],color=col,/device,charsize=list.letter,NOCLIP = 0 ;,ORIENTATION =90
                endif
            endfor
            if list.showlegend then begin
                if list.blegendshort then begin
                    p=strpos((*list.PDFname)[i],';')
                    if p ne -1 then begin
                        p2=strpos((*list.PDFname)[i],';',p+1)
                        if p2 ne -1 then p = p2
                        tmp=strmid((*list.PDFname)[i],0,p)
                    endif else tmp=(*list.PDFname)[i]
                endif else tmp=(*list.PDFname)[i]
                if list.bPDFpsym then begin
                    xyouts,xpos,ypos,'!3'+tmp,/normal,col=col,charsize=list.letter;,pow=list.ypower
                    plots,xpos*0.95,ypos*1.01,/normal,col=col,psym=psym,symsize=0.7;,pow=list.ypower
                endif else $
                    xyouts,xpos,ypos,'!3_____ '+tmp,/normal,col=col,charsize=list.letter;,pow=list.ypower
                
                ypos+=incpos
            endif
        endif

        off+=(*list.PDFn)[i]
    endfor
endif

;----Display mask----
if (list.masknr ge 1) then begin

GetCRange,ROIy=ROIy
b=ROIy[0]+0.1*(ROIy[1]-ROIy[0])
e=ROIy[0]+0.9*(ROIy[1]-ROIy[0])
y=b+indgen(list.masknr)*(e-b+1.)/list.masknr

for i=0l,list.masknr-1 do begin
if (*list.smask)[i] eq 1 then begin
    if ((*list.mask)[1,i] eq list.maskP) then color=col8 else color=col7
    arr=reform((*list.mask)[*,i])
    str=stringr(fix(i))+';'+stringr(fix(arr[1]))

    t=min(abs((*list.xvalt)[*,0]-arr[3]),x1)
    t=min(abs((*list.xvalt)[*,0]-arr[4]),x2)
    x=(*list.xvalt)[[x1,x2],list.xtype]
    x=x[sort(x)]

    Arrow, x[0], y[i], x[1], y[i], /Data, /Solid, HSize=12,color=color,NOCLIP = 0;,pow=list.ypower
    PlotS, [x[0], x[0]], ROIy, /Data,linestyle=2,NOCLIP = 0,color=color;,pow=list.ypower
    PlotS, [x[1], x[1]], ROIy, /Data,linestyle=2,NOCLIP = 0,color=color;,pow=list.ypower

    if arr[2] eq 1 then xyouts,x[0],y[i],str,/data,charsize=list.letter,color=color,NOCLIP = 0;,pow=list.ypower

endif
endfor
endif

end;pro plotCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotres,list
fitmodelplot,list,[0,0,16],/residual
end;plotres
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro peaksearch,ev,extra
widget_control,ev.top,get_uvalue=list
if list.apeaks EQ 0 then return

index=list.peakselect+extra
if index lt 0 then index = list.apeaks+index
if index gt (list.apeaks-1) then index=index-list.apeaks
index=0>index<(list.apeaks-1)
list.peakselect=index

widget_control,ev.top,set_uvalue=list
RefreshDisplayCHI,list

ID1=widget_info(ev.top,FIND_BY_UNAME='x-value')
xen=CHI_xval(list.peakxen[list.peakselect],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtypeorig)
str1=strtrim(string(format='(F15.3)',xen[list.xtype]),2)+list.xunits[list.xtype]
str2=strtrim(string(format='(F15.3)',list.peaky[list.peakselect]),2)
if ptr_valid(list.speringpr) then str3=strtrim(string(format='(F15.3,"%")',(*list.speringpr)[list.peakxkan[list.peakselect]]),2) else str3=''
widget_control,ID1,set_value=[str1,str2,str3]

ID1=widget_info(ev.top,FIND_BY_UNAME='peaksin')
widget_control,ID1,SET_TABLE_SELECT=[ 0, list.peakselect, 2, list.peakselect ]
widget_control,ID1,SET_TABLE_VIEW=[ 0, list.peakselect]

end; pro peaksearch
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro fpeaksearch,ev,extra,type
widget_control,ev.top,get_uvalue=list
if list.apeaksout EQ 0 then return

index=list.fitselect+extra
if index lt 0 then index = list.apeaksout+index
if index gt (list.apeaksout-1) then index=index-list.apeaksout
index=0>index<(list.apeaksout-1)
list.fitselect=index

widget_control,ev.top,set_uvalue=list
RefreshDisplayCHI,list

case type of
0:    begin
    ID1=widget_info(ev.top,FIND_BY_UNAME='table')
    widget_control,ID1,set_value=transpose([[list.peakout[*,list.fitselect]],$
       [list.peakouts[*,list.fitselect]]]),USE_TABLE_SELECT=[0,0,1,2]
    endcase
1:    begin
    ID1=widget_info(ev.top,FIND_BY_UNAME='peaksout')
    widget_control,ID1,SET_TABLE_SELECT=[ 0, list.fitselect, 2, list.fitselect ]
    widget_control,ID1,SET_TABLE_VIEW=[ 0, list.fitselect]
    endcase
endcase

end; pro fpeaksearch
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Add_Peak_Event, ev

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

list.xs=ev.x
list.ys=ev.y

SetSysVar,list.sysvar
coords = Convert_Coord(ev.x,ev.y, /Device, /To_Data)
temp=min(abs((*list.xvalt)[*,list.xtype]-coords[0]),index)

WSet, list.drawindex

x=(*list.xvalt)[index,list.xtype]
y=(*list.spet)[index,list.xtype]
plots,[x,x],[list.ylog,y],$
       linestyle=line,NOCLIP = 0
oplot,[x],[y],psym=2,NOCLIP = 0

h=!D.x_size
v=!D.y_size
WSet, list.pixindex
Device, Copy = [0,0, h,v, 0,0,list.drawindex]

Widget_Control, ev.id, Event_Pro='Add_Move_Event', $
    draw_motion=1
widget_control,ev.top,set_uvalue=list
end;pro Add_Peak_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Add_Move_Event,ev

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

WSet, list.drawindex
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.pixindex]
WSet, list.drawindex
SetSysVar,list.sysvar
coords = Convert_Coord([list.xs,ev.x],[list.ys,ev.y], /Device, /To_Data)
temp=min(abs((*list.xvalt)[*,list.xtype]-coords[0,0]),index0)
temp=min(abs((*list.xvalt)[*,list.xtype]-coords[0,1]),index1)

IF thisEvent EQ 'UP' THEN BEGIN
    Widget_Control, ev.id, draw_motion=0, $
       Event_Pro='Add_Peak_Event'
    x=(*list.xvalt)[*,list.xtypeorig]
    y=(*list.spet)[*,list.xtypeorig]
    peakbegin=index0<index1
    peakend=index0>index1
    peakval=estpeak(x,y,peakbegin,peakend,debug=0,drawi=0)

    list.peakbegin[list.apeaks]=peakbegin
    list.peakend[list.apeaks]=peakend
    list.peakxen[list.apeaks]=peakval[0]
    list.peakxkan[list.apeaks]=peakval[1]
    list.peakarea[list.apeaks]=peakval[2]
    list.peaksigma[list.apeaks]=peakval[3]
    list.peaky[list.apeaks]=peakval[4]
    list.apeaks=list.apeaks+1

    RefreshDisplayCHI,list
    Widget_Control, ev.top, Set_UValue=list
    RefreshWindowCHI,ev
    return
endif

WSet, list.drawindex
x=(*list.xvalt)[index1,list.xtype]
y=(*list.spet)[index1,list.xtype]
plots,[x,x],[list.ylog,y],linestyle=line,NOCLIP = 0
oplot,[x],[y],psym=2,NOCLIP = 0

end;pro Add_Move_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Delete_Peak_Event, ev

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
if list.apeaks eq 0 then return

SetSysVar,list.sysvar
coords = Convert_Coord(ev.x, ev.y, /Device, /To_Data)
x=CHI_xval(coords[0],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtype)
temp=min(abs(list.peakxen-x[list.xtypeorig]),index)

arr=[[list.peaky[0:list.apeaks-1]],$
    [list.peakxen[0:list.apeaks-1]],$
    [list.peakxkan[0:list.apeaks-1]],$
    [list.peaksigma[0:list.apeaks-1]],$
    [list.peakarea[0:list.apeaks-1]],$
    [list.peakbegin[0:list.apeaks-1]],$
    [list.peakend[0:list.apeaks-1]]]
arr2=ShrinkArray(arr,index)

list.apeaks=list.apeaks-n_elements(index)
peaky=list.peaky*0
list.peaky=peaky
list.peakxen=peaky
list.peakxkan=peaky
list.peaksigma=peaky
list.peakarea=peaky
list.peakbegin=peaky
list.peakend=peaky

if list.apeaks ne 0 then begin
list.peaky=arr2[*,0]
list.peakxen=arr2[*,1]
list.peakxkan=arr2[*,2]
list.peaksigma=arr2[*,3]
list.peakarea=arr2[*,4]
list.peakbegin=arr2[*,5]
list.peakend=arr2[*,6]
endif

list.peakselect=0>(list.peakselect<(list.apeaks-1))

widget_control,ev.top,set_uvalue=list
RefreshDisplayCHI,list
RefreshWindowCHI,ev
end;pro Delete_Peak_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro UpdatePeakSelectText,ev,list
ID1=widget_info(ev.top,FIND_BY_UNAME='x-value2')
str1=strtrim(string(format='(F15.3)',list.anyselect[list.xtype]),2)+' '+list.xunits[list.xtype]
str2=strtrim(string(format='(F15.3)',(*list.spet)[0>list.anyselect[4]<(list.nxval-1),list.xtype]),2)+' a.u.'
str3=msymbols('plus_minus',/widget)+strtrim(string(format='(F15.3)',(*list.speerrort)[0>list.anyselect[4]<(list.nxval-1),list.xtype]),2)+' a.u.'
if ptr_valid(list.speringpr) then str4=strtrim(string(format='(F15.3,"%")',(*list.speringpr)[list.anyselect[4]]),2) else str4='... %'
widget_control,ID1,set_value=[str1,str2,str3,str4]
end;pro UpdatePeakSelectText
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Select_Any_Event, ev

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
wset,list.drawindex
SetSysVar,list.sysvar
coords = Convert_Coord(ev.x, ev.y, /Device, /To_Data)
x=CHI_xval(coords[0],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtype)
temp=min(abs((*list.xvalt)[*,list.xtypeorig]-x[list.xtypeorig]),index)
list.anyselect=[[x],[index]]

widget_control,ev.top,set_uvalue=list
RefreshDisplayCHI,list
UpdatePeakSelectText,ev,list

if widget_info(list.main,/valid) and list.CrossUpdate then begin
    SRDirectGraph,struct
    MEllipseSet,list.anyselect[list.xtype],list.xtype,list.main,setself=ev.top
    SRDirectGraph,struct
endif

end;pro Select_Any_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Select_Peak_Event, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

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
if list.apeaks eq 0 then return
SetSysVar,list.sysvar
coords = Convert_Coord(ev.x, ev.y, /Device, /To_Data)
x=CHI_xval(coords[0],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtype)
temp=min(abs(list.peakxen-x[list.xtypeorig]),index)
list.peakselect=index

widget_control,ev.top,set_uvalue=list
RefreshDisplayCHI,list

ID1=widget_info(ev.top,FIND_BY_UNAME='x-value')
xen=CHI_xval(list.peakxen[list.peakselect],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtypeorig)
str1=strtrim(string(format='(F15.3)',xen[list.xtype]),2)+list.xunits[list.xtype]
str2=strtrim(string(format='(F15.3)',list.peaky[list.peakselect]),2)
if ptr_valid(list.speringpr) then str3=strtrim(string(format='(F15.3,"%")',(*list.speringpr)[list.peakxkan[list.peakselect]]),2) else str3=''
widget_control,ID1,set_value=[str1,str2,str3]
end;pro Select_Peak_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ROI_PROCESS_EVENTS_CHI, ev

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
list.xs=ev.x
list.ys=ev.y

WSet, list.drawindex
SetSysVar,list.sysvar
coords = Convert_Coord(ev.x, ev.y, /Device, /To_Data)
GetCRange,ROIx=ROIx,ROIy=ROIy
x1 = ROIx[0] > coords[0] < ROIx[1]
y1 = ROIy[0] > coords[1] < ROIy[1]
PlotS, [x1, x1], ROIy

h=!D.x_size
v=!D.y_size
WSet, list.pixindex
Device, Copy = [0,0, h,v, 0,0,list.drawindex]

Widget_Control, ev.id, Event_Pro='ROI_DRAWBOX', $
    draw_motion=1
Widget_Control, ev.top, Set_UValue=list
end;pro ROI_PROCESS_EVENTS_CHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ROI_DRAWBOX,ev

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

WSet, list.drawindex
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.pixindex]

x = [list.xs, ev.x]
y = [list.ys, ev.y]
SetSysVar,list.sysvar
coords = Convert_Coord(x, y, /Device, /To_Data)
GetCRange,ROIx=ROIx,ROIy=ROIy
x1 = ROIx[0] > coords[0,0] < ROIx[1]
x2 = ROIx[0] > coords[0,1] < ROIx[1]

IF thisEvent EQ 'UP' THEN BEGIN
    Widget_Control, ev.id, draw_motion=0, $
       Event_Pro='ROI_PROCESS_EVENTS_CHI'

    if x1 ne x2 then begin
        temp=value_locater((*list.xvalt)[*,list.xtype],[x1,x2])>0
        list.backranind=temp[sort(temp)]
        temp=(*list.xvalt)[temp,0]
        list.backranindD=temp[sort(temp)]
    endif else begin
          list.backranind=[0,list.nxval-1]
          temp=(*list.xvalt)[list.backranind,0]
          list.backranindD=temp[sort(temp)]
    endelse

    RefreshDisplayCHI,list
    Widget_Control, ev.top, Set_UValue=list
    return
endif


list.xd=ev.x
list.yd=ev.y

y2 = ROIy[0] > coords[1,1] < ROIy[1]
WSet, list.drawindex
Arrow, x1, y2, x2, y2, /Data, /Solid, HSize=12
PlotS, [x2, x2], ROIy

Widget_Control, ev.top, Set_UValue=list
end;pro ROI_DRAWBOX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ZOOM_PROCESS_EVENTS_CHI, ev

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
list.xs=ev.x
list.ys=ev.y
Widget_Control, ev.id, Event_Pro='ZOOM_DRAWBOX_CHI', $
    draw_motion=1
Widget_Control, ev.top, Set_UValue=list
end;pro ZOOM_PROCESS_EVENTS_CHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ZOOM_DRAWBOX_CHI ,ev

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

WSet, list.drawindex
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.pixindex]
x = [list.xs, ev.x]
y = [list.ys, ev.y]

GetCRange,ROIx=ROIx,ROIy=ROIy

IF thisEvent EQ 'UP' THEN BEGIN
    Widget_Control, ev.id, draw_motion=0, $
       Event_Pro='ZOOM_PROCESS_EVENTS_CHI'

    IF list.xs GT ev.x THEN x = [ev.x, list.xs]
    IF list.ys GT ev.y THEN y = [ev.y, list.ys]
    SetSysVar,list.sysvar
    coords = Convert_Coord(x, y, /Device, /To_Data)
    x1 = ROIx[0] > coords[0,0] < ROIx[1]
    x2 = ROIx[0] > coords[0,1] < ROIx[1]
    y1 = ROIy[0] > coords[1,0] < ROIy[1]
    y2 = ROIy[0] > coords[1,1] < ROIy[1]
    if list.ypower ne 0 and list.ypower ne 1 then begin
        y1^=1/list.ypower
        y2^=1/list.ypower
    endif
    temp=min(abs((*list.xvalt)[*,list.xtype]-x1),ix1)
    temp=min(abs((*list.xvalt)[*,list.xtype]-x2),ix2)
    list.xrange = [ix1,ix2]
    temp=min(abs((*list.spet)[*,list.xtype]-y1),iy1)
    temp=min(abs((*list.spet)[*,list.xtype]-y2),iy2)
    list.yrange=[iy1,iy2]
    RefreshDisplayCHI,list

    Widget_Control, ev.top, Set_UValue=list
    return
endif

SetSysVar,list.sysvar
coords = Convert_Coord(x, y, /Device, /To_Data)
x1 = ROIx[0] > coords[0,0] < ROIx[1]
x2 = ROIx[0] > coords[0,1] < ROIx[1]
y1 = ROIy[0] > coords[1,0] < ROIy[1]
y2 = ROIy[0] > coords[1,1] < ROIy[1]
WSet, list.drawindex
PlotS, [x1, x1], [y1,y2],Linestyle=1
PlotS, [x2, x2], [y1,y2],Linestyle=1
PlotS, [x1,x2],[y1, y1],Linestyle=1
PlotS, [x1,x2],[y2, y2],Linestyle=1
Widget_Control, ev.top, Set_UValue=list
end;pro ZOOM_DRAWBOX_CHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_XRD_CHI, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list

if n_elements(list) ne 0 then begin
    wdelete,list.pixindex
    result=UpdateINI('$XRD_CHI path:',{path:list.path,file:list.file})
    result=UpdateINI('$XRD mskpath:',{path:list.mskpath,file:list.mskfile})
    result=UpdateINI('$XRD pddpath:',{path:list.pddpath,file:list.pddfile})
    result=UpdateINI('$CrossUpdate:',{n:list.CrossUpdate})

    if widget_info(list.main,/VALID_ID) then begin
        widget_control,list.main,get_uvalue=list2
        ind=where((*list2.CHIchild) eq ID)
        (*list2.CHIchild)=ShrinkArray((*list2.CHIchild),ind,count=ct)
        if ct eq 0 then begin
            ptr_free,list2.CHIchild
            widget_control,list.main,set_uvalue=list2
        endif
        list2.chipath=list.path
        list2.chifile=list.file
    endif
    ;set_plot, list.Platform
    ;device,list.Visual
    heap_free,list
endif

error=UnSubscribeToOutID()
print,'=========================================STOP 1D Profile Editor'
end;pro CleanUp_XRD_CHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AutoCorrect1D,list
if list.backmod ne 0 then begin
    x1=list.backranind[0]
    x2=list.backranind[1]
    
    case list.backmod of
    1:    BEGIN
        degree=(*(*list.ptrModel).ptrFit).backparamIO[0]
        if degree ge 0 then y=ortpol((*list.xvalt)[x1:x2,list.xtype],(*list.spet)[x1:x2,list.xtype],degree=degree,print=list.OutID) $
        else y=ortpol((*list.xvalt)[x1:x2,list.xtype],(*list.spet)[x1:x2,list.xtype],GR=degree,print=list.OutID)
        (*(*list.ptrModel).ptrFit).backparamIO[0]=degree    
        endcase
    2:    begin
        y=Strip1D((*list.spet)[x1:x2,list.xtype],(*(*list.ptrModel).ptrFit).backparamIO[2],(*(*list.ptrModel).ptrFit).backparamIO[1])
        endcase
    else: return
    endcase
    
    if not ptr_valid(list.yback) then list.yback=ptr_new((*list.spet)[*,list.xtype]*0)
    (*list.yback)[x1:x2]=y
    
    if list.backcor ne 0 then begin
        spe=(*list.spet)[*,list.xtype]
        BackSub,spe,*list.yback,list.backmod,offset=offset
        (*list.spet)[*,list.xtype]=spe
    endif
endif
end;pro AutoCorrect1D
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro NewUnits,ev
widget_control,ev.top,get_uvalue=list

; Determine type
list.xtype++
list.xtype mod= n_elements(list.xmodes)

; Set X-axis: values, title
ind=sort((*list.xvalt)[list.xrange,list.xtype])
list.xrange=list.xrange[ind]
list.xtitle=ChiXtitle(list.xtype)
widget_control,ev.id,set_value='X: '+list.xmodes[list.xtype]

; Corrections
(*list.spet)=*list.speorg
(*list.speerrort)=(*list.speerrororg)
AutoCorrect1D,list

RefreshPDD_CHI,list
widget_control,ev.top,set_uvalue=list
;RefreshUnitsCHI,ev
RefreshDisplayCHI,list
RefreshWindowCHI,ev
end;pro NewUnits
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshUnitsCHI,ev
widget_control,ev.top,get_uvalue=list

RefreshPDD_CHI,list
RefreshDisplayCHI,list
RefreshWindowCHI,ev
end;pro RefreshUnitsCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshDataCHI,ev
widget_control,ev.top,get_uvalue=list
; ----x values----
xvalt=CHI_xval((*list.xvalt)[*,list.xtypeorig],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtypeorig,$
                spe=(*list.speorg)[*,list.xtypeorig],tspe=spet,$
                stdspe=(*list.speerrororg)[*,list.xtypeorig],tstdspe=speerrort)
(*list.xvalt)=temporary(xvalt)
(*list.speorg)=temporary(spet)
(*list.spet)=*list.speorg
(*list.speerrororg)=temporary(speerrort)
(*list.speerrort)=*list.speerrororg

temp=value_locater((*list.xvalt)[*,0],list.backranindD)>0
list.backranind=temp[sort(temp)]

; ----Labels x-values----
if list.alabels ne 0 then $
    (*list.xlabels)=CHI_xval((*list.xlabels)[*,list.xtypeorig],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtypeorig)

; ----Corrections----
AutoCorrect1D,list
RefreshPDD_CHI,list

widget_control,ev.top,set_uvalue=list
;RefreshUnitsCHI,ev
RefreshDisplayCHI,list
RefreshWindowCHI,ev
end;pro RefreshDataCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshWindowCHI,ev
widget_control,ev.top,update=0
widget_control,ev.top,get_uvalue=list
ID=widget_info(ev.top,FIND_BY_UNAME='PDFname')
if widget_info(ID,/VALID_ID) then widget_control,ID,set_value=(ptr_valid(list.PDFname)?(*list.PDFname):''),ysize=(list.nPDF)<20

;ID1=widget_info(ev.top,FIND_BY_UNAME='x-value')
;xen=CHI_xval(list.peakxen[list.peakselect],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtypeorig)
;str1=strtrim(string(format='(F15.3)',xen[list.xtype]),2)+list.xunits[list.xtype]
;str2=strtrim(string(format='(F15.3)',list.peaky[list.peakselect]),2)
;if ptr_valid(list.speringpr) then str3=strtrim(string(format='(F15.3,"%")',(*list.speringpr)[list.peakxkan[list.peakselect]]),2) else str3=''
;widget_control,ID1,set_value=[str1,str2,str3]

if ptr_valid(list.spet) then UpdatePeakSelectText,ev,list

ID=widget_info(ev.top,FIND_BY_UNAME='degree')
if widget_info(ID,/VALID_ID) then widget_control,ID,set_value=stringr((*(*list.ptrModel).ptrFit).backparamIO[0])

ID=widget_info(ev.top,FIND_BY_UNAME='Width')
if widget_info(ID,/VALID_ID) then widget_control,ID,set_value=stringr((*(*list.ptrModel).ptrFit).backparamIO[2])

ID=widget_info(ev.top,FIND_BY_UNAME='Niter')
if widget_info(ID,/VALID_ID) then widget_control,ID,set_value=stringr((*(*list.ptrModel).ptrFit).backparamIO[1])

ID=widget_info(ev.top,FIND_BY_UNAME='Y')
if list.ylog then widget_control,ID,set_value='Y: logarithmic' $
else widget_control,ID,set_value='Y: linear'

ID=widget_info(ev.top,FIND_BY_UNAME='X')
if list.xlog then widget_control,ID,set_value='X: logarithmic' $
else widget_control,ID,set_value='X: linear'

ID=widget_info(ev.top,FIND_BY_UNAME='peaks')
if list.peaksflag then widget_control,ID,set_value='Hide peaks' $
else widget_control,ID,set_value='Show peaks'

ID=widget_info(ev.top,FIND_BY_UNAME='peakselect')
if list.peakselectflag then widget_control,ID,set_value='Hide peakselect' $
else widget_control,ID,set_value='Show peakselect'

ID=widget_info(ev.top,FIND_BY_UNAME='view background')
if list.sBackground then widget_control,ID,set_value='Hide background' $
else widget_control,ID,set_value='Show background'

ID=widget_info(ev.top,FIND_BY_UNAME='ringpr')
if list.sringpr then widget_control,ID,set_value='Hide ring% values' $
else widget_control,ID,set_value='Show ring% values'

ID=widget_info(ev.top,FIND_BY_UNAME='Labels')
if list.sLabels then widget_control,ID,set_value='Hide Labels' $
else widget_control,ID,set_value='Show Labels'

ID=widget_info(ev.top,FIND_BY_UNAME='Legend')
if list.showlegend then widget_control,ID,set_value='Hide Legend' $
else widget_control,ID,set_value='Show Legend'

ID=widget_info(ev.top,FIND_BY_UNAME='sROI')
if list.plotROI then widget_control,ID,set_value='Hide Data range' $
else widget_control,ID,set_value='Show Data range'

ID=widget_info(ev.top,FIND_BY_UNAME='Xunit')
widget_control,ID,set_value='X: '+list.xmodes[list.xtype]

ID=widget_info(ev.top,FIND_BY_UNAME='PlotFit')
case list.fitresolve of
    0:widget_control,ID,set_value='Fit hide'
    1:widget_control,ID,set_value='Fit total'
    2:widget_control,ID,set_value='Fit resolved'
endcase

base4=widget_info(ev.top,FIND_BY_UNAME='RepBaseI')
if widget_info(base4, /VALID_ID) then begin
    ID=widget_info(ev.top,FIND_BY_UNAME='peaksin')
    if widget_info(ID, /VALID_ID) then widget_control,ID,/destroy
    MakePeaksin,list,base4
endif

base4=widget_info(ev.top,FIND_BY_UNAME='RepBaseR')
if widget_info(base4, /VALID_ID) then begin
    ID=widget_info(ev.top,FIND_BY_UNAME='peaksout')
    if widget_info(ID, /VALID_ID) then widget_control,ID,/destroy
    MakePeaksout,list,base4
endif

ID=widget_info(ev.top,FIND_BY_UNAME='Chisq')
if widget_info(ID, /VALID_ID) then begin
    widget_control,ID,set_value='Chisq Red: '+stringr(list.chisq)

    ID=widget_info(ev.top,FIND_BY_UNAME='Iter')
    widget_control,ID,set_value='Iter: '+stringr(list.iter)

    if (list.apeaksout ne 0) then begin
        ID=widget_info(ev.top,FIND_BY_UNAME='table')
        widget_control,ID,set_value=transpose([[list.peakout[*,list.fitselect]],$
                     [list.peakouts[*,list.fitselect]]]),USE_TABLE_SELECT=[0,0,1,2]
    endif else begin
        ID=widget_info(ev.top,FIND_BY_UNAME='table')
        widget_control,ID,set_value=fltarr(1,4),USE_TABLE_SELECT=[0,0,1,2]
    endelse
endif

ID=widget_info(ev.top,FIND_BY_UNAME='Subtract')
if widget_info(ID, /VALID_ID) then widget_control,ID,sensitive=~list.backcor
ID=widget_info(ev.top,FIND_BY_UNAME='Restore')
if widget_info(ID, /VALID_ID) then widget_control,ID,sensitive=list.backcor

widget_control,ev.top,update=1
end;pro RefreshWindowCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshDisplayCHI,list
wset, list.drawindex
plotCHI,list
GetSysVar,list.sysvar
wset, list.pixIndex
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.drawindex]
if list.resindex ne -1 then begin
    wset,list.resindex
    plotres,list
endif
end;pro RefreshDisplayCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro BatchConvertCHI,top,list

; Files
list2={path:list.path,file:list.file,Format:fileformat(1),OutID:list.OutID}
path=Select_Files(top,list2,sortmethod=list.sortmethod,separator=list.sortseparator)
if path[0] eq '' then return
nfiles=n_elements(path)
Result = DIALOG_MESSAGE(stringr(nfiles)+' files to convert: proceed?',/question)
if result eq 'No' then return

; New X-axis
ProfileRangeWindow,top,xtypeconv,x0,xi,xn,xdefault,error=error
if error then return

; Convert
brev=xtypeconv eq 0

for i=0l,nfiles-1 do begin
    printw,list.OutID,path[i]
    tmp=CutPath(path[i],path=p,file=f,ext=ext)
    error=ReadScan(p,f+ext,ext,list.OutID,xval=xval,spe=spe,stdspe=speerror,type=xtype,ROIAzimuth=ROIAzimuth,asciitemplate=list.asciitemplate)
    if error ne 0 then begin
        xvalt=CHI_xval(xval,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,xtype,spe=spe,tspe=spet,stdspe=speerror,tstdspe=speerrort)
        xval=xvalt[*,xtypeconv]
        spe=spet[*,xtypeconv]
        speerror=speerrort[*,xtypeconv]
        if brev then begin
            xval=reverse(xval)
            spe=reverse(spe)
            speerror=reverse(speerror)
        endif

        if xdefault eq 1 then begin
            xn=n_elements(xval)
            x0=min(xval,max=xi)
            xi=(xi-x0)/(xn-1.)
        endif
        xconv = x0+xi*lindgen(xn)
        yconv = minterpol( spe, xval, xconv ,/SPLINE)
        stdyconv = minterpol( speerror, xval, xconv ,/SPLINE)

        temp=WriteChi(p+'_'+f+'.chi',transpose([[xconv],[yconv],[stdyconv]]),xtypeconv,path[i],list.OutID,ROIAzimuth,berror=1)
    endif
endfor

end;pro BatchConvertCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SaveSingleChi,list
path=SelFile(list.path,list.file,ChiFormat(),'Save Profile...')
if path eq '' then return

if WriteChi(path,transpose([[(*list.xvalt)[*,list.xtype]],[(*list.spet)[*,list.xtype]],[(*list.speerrort)[*,list.xtype]]]),list.xtype,'Empty',list.OutID,list.ROIAzimuth,berror=1b) then begin
    error=CutPath(path,path=path,file=file,ext=Format)
    list.path=path
endif

end;pro SaveSingleChi
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro OpenCHI,ev,prompt,fileptr=fileptr,xdinr=xdinr,dummy=dummy
widget_control,ev.top,get_uValue=list

; ----Find file to read----
if keyword_set(dummy) then begin
    path=''
    file='dummy.chi'
    Format='.chi'
endif else if prompt then begin
    path=''
    filter=fileformat(1)
    sortformat,filter,list.Format
    path=CAPDIALOG_PICKFILE(path=list.path,file=list.file,filter='*'+filter)
    if path EQ '' then return
    error=CutPath(path,path=path,file=file,ext=Format)
    file=file+Format
endif else begin
    path=list.path
    file=list.file
    Format=list.Format
endelse

; ----Read file----
if keyword_set(dummy) then begin
    error = ReadScanDummy(ev.top,xval=xval,spe=spe,stdspe=speerror,title=title,$
                          YTITLE=YTITLE,xrange=xrange,yrange=yrange,type=xtype,ROIAzimuth=ROIAzimuth)
endif else if n_elements(fileptr) ne 0 then begin
    xval=fileptr.xval
    spe=fileptr.spe
    speerror=fileptr.speerror
    title=fileptr.title
    xtitle=fileptr.xtitle
    ytitle=fileptr.ytitle
    xrange=fileptr.xrange
    yrange=fileptr.yrange
    xtype=fileptr.xtype
    ROIAzimuth=fileptr.ROIAzimuth
    file=fileptr.title+'.chi'
endif else begin
    error=ReadScan(path,file,Format,list.OutID,xval=xval,spe=spe,stdspe=speerror,title=title,$
        YTITLE=YTITLE,xrange=xrange,yrange=yrange,type=xtype,ROIAzimuth=ROIAzimuth,ringpr=ringpr,xdinr=xdinr,$
        asciitemplate=list.asciitemplate)
    if error eq 0 then return
endelse

; ----x values----
if (size(spe))[0] gt 1 then begin
    spe=total(spe,2)
    speerror=sqrt(total(speerror^2,2))
endif
xvalt=CHI_xval(xval,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,xtype,spe=spe,tspe=spet,stdspe=speerror,tstdspe=speerrort)

; ----Edit data structure----
list.path=path
list.file=file
list.Format=Format
ptr_free,list.xvalt
list.xvalt=ptr_new(temporary(xvalt))
ptr_free,list.spet
list.spet=ptr_new(temporary(spet))
ptr_free,list.speorg
list.speorg=ptr_new(*list.spet)
ptr_free,list.speerrort
list.speerrort=ptr_new(temporary(speerrort))
ptr_free,list.speerrororg
list.speerrororg=ptr_new(*list.speerrort)
ptr_free,list.speringpr
if n_elements(ringpr) ne 0 then list.speringpr=ptr_new(ringpr)
ptr_free,list.yback
ptr_free,list.yfit
;ptrFit=(*list.ptrModel).ptrFit
;ptr_free,(*(*ptrFit).yfit)[(*ptrFit).icurrent]
list.title=title
list.ytitle=ytitle
if list.listtype eq 10 then begin
    list.xrange=xrange
    list.yrange=yrange
endif else begin
    ; Full scale
    if prompt then list.yrange=0
    ;list.yrange=0
endelse
;list.xtype=xtype
list.xtypeorig=xtype
list.xtitle=ChiXtitle(list.xtype)
list.ROIAzimuth=ROIAzimuth
list.nxval=n_elements(xval)

widget_control,ev.top,TLB_SET_TITLE='1D Profile Editor: '+list.title
if list.listtype eq 10 then begin
    list.backranind=[0,n_elements(spe)-1]
    temp=(*list.xvalt)[[0,list.nxval-1],0]
    list.backranindD=temp[sort(temp)]

    list.listtype=11
    widget_control,list.draw,get_value=drawindex
    list.drawindex=drawindex
    ; ----Refresh these only in the procedure----
    ID3=widget_info(ev.top,FIND_BY_UNAME='base1')
    ID16=widget_info(ev.top,FIND_BY_UNAME='Save...')
    ID17=widget_info(ev.top,FIND_BY_UNAME='Load...')
    ID19=widget_info(ev.top,FIND_BY_UNAME='base2')
    ID22=widget_info(ev.top,FIND_BY_UNAME='MOptions')
    ID23=widget_info(ev.top,FIND_BY_UNAME='MDisplay')
    ID24=widget_info(ev.top,FIND_BY_UNAME='Perform')
    widget_control,ID3,sensitive=1
    widget_control,ID16,sensitive=1
    widget_control,ID17,sensitive=1
    widget_control,ID19,sensitive=1
    widget_control,ID22,sensitive=1
    widget_control,ID23,sensitive=1
    widget_control,ID24,sensitive=1
    
    RefreshPDD_CHI,list

    ; ----Load Mask and Refresh----
    if list.main ne 0 then error=0b else LoadMask,list,0,list.OutID,error=error,/oD
    if error then begin
       path=file_search(list.mskpath+'*.msk',count=count)
       if count ne 0 then begin
         path=path[0]
         error=CutPath(path,path=mskpath,file=file,ext=ext)
         list.mskpath=mskpath
         list.mskfile=file+ext
         LoadMask,list,0,list.OutID,/oD
         widget_control,ev.top,set_uvalue=list
         RefreshDataCHI,ev
       endif else begin
            widget_control,ev.top,set_uvalue=list
         RefreshDisplayCHI,list
         RefreshWindowCHI,ev
       endelse
    endif else begin
        widget_control,ev.top,set_uvalue=list
        RefreshDataCHI,ev
    endelse
endif else begin
    temp=value_locater((*list.xvalt)[*,0],list.backranindD)>0
    list.backranind=temp[sort(temp)]

    AutoCorrect1D,list
    RefreshPDD_CHI,list
    
    widget_control,ev.top,set_uvalue=list
    RefreshDisplayCHI,list
    RefreshWindowCHI,ev
endelse

end;pro OpenCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRD_CHI_event, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
0:    begin
    case ev.action of
    'refresh':    begin
                widget_control,ev.top,get_uvalue=list
                RefreshDisplayCHI,list
                endcase
    'refreshdata':RefreshDataCHI,ev
    endcase
    endcase
3:    BEGIN ; event van text_widget
    widget_control,ev.top,get_uvalue=list
    widget_control,ev.id,get_uvalue=uvalue
    case uvalue of
    'degree':  begin
              widget_control,ev.id,get_value=degree
            (*(*list.ptrModel).ptrFit).backparamIO[0]=long(degree)
            AutoCorrect1D,list
            widget_control,ev.id,set_value=strtrim(string((*(*list.ptrModel).ptrFit).backparamIO[0]),2)
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list
            endcase
    'Niter':begin
              widget_control,ev.id,get_value=niter
              niter=long(niter)
            (*(*list.ptrModel).ptrFit).backparamIO[1]=niter
            AutoCorrect1D,list
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list
            endcase
    'Width':begin
              widget_control,ev.id,get_value=width
              width=long(width)
            (*(*list.ptrModel).ptrFit).backparamIO[2]=width
            AutoCorrect1D,list
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list
            endcase
    'ulab': begin
            WIDGET_CONTROL, ev.top, get_uvalue=list
              widget_control,ev.id,get_value=lab
              if list.alabels eq 0 then begin
                ptr_free,list.labels
                ptr_free,list.xlabels
                ptr_free,list.ylabels
                list.labels=ptr_new(lab)
                if list.ylog eq 1 then list.ylabels=ptr_new((*list.spet)[(list.yrange[0]+list.yrange[1])/2,list.xtype]>1) $
                  else list.ylabels=ptr_new((*list.spet)[(list.yrange[0]+list.yrange[1])/2,list.xtype])
                list.xlabels=ptr_new((*list.xvalt)[(list.xrange[0]+list.xrange[1])/2,*])
            endif else begin
                (*list.labels)=[(*list.labels),lab]
                if list.ylog eq 1 then (*list.ylabels)=[(*list.ylabels),(*list.spet)[(list.yrange[0]+list.yrange[1])/2,list.xtype]>1] $
                  else (*list.ylabels)=[(*list.ylabels),(*list.spet)[(list.yrange[0]+list.yrange[1])/2,list.xtype] ]
                (*list.xlabels)=[(*list.xlabels),(*list.xvalt)[(list.xrange[0]+list.xrange[1])/2,*]]
            endelse
            list.alabels++

              widget_control,ev.top,set_uvalue=list
              RefreshDisplayCHI,list
            endcase
 'PDFSubstr':begin
             WIDGET_CONTROL, ev.top, get_uvalue=list,/no_copy
              widget_control,ev.id,get_value=lab
              list.PDFSearch.PDFSubstr=lab
              widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
 'PDFdi':    begin
             WIDGET_CONTROL, ev.top, get_uvalue=list,/no_copy
              widget_control,ev.id,get_value=lab1
              lab1=float(lab1)
              lab2=list.PDFSearch.PDFde
              list.PDFSearch.PDFdi=lab1<lab2
              list.PDFSearch.PDFde=lab1>lab2
              widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
 'PDFde':    begin
             WIDGET_CONTROL, ev.top, get_uvalue=list,/no_copy
              widget_control,ev.id,get_value=lab1
              lab1=float(lab1)
              lab2=list.PDFSearch.PDFdi
              list.PDFSearch.PDFdi=lab1<lab2
              list.PDFSearch.PDFde=lab1>lab2
              widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
    endcase
    endcase
1:     begin ;event from button
    WIDGET_CONTROL, ev.id, GET_VALUE=val
    case val OF
    'Open dummy':begin
            OpenCHI,ev,1,/dummy
            endcase
    'Open...':begin
            OpenCHI,ev,1
            endcase
     'Exit':begin
            widget_control,ev.top,/destroy
            endcase
'Move Legend':begin
            modeCHI,ev,val
            endcase
'Data range' :     begin
            modeCHI,ev,val
            endcase
'Zoom'  :    begin
            modeCHI,ev,val
            endcase
'Y: logarithmic' :begin
            widget_control,ev.top,get_uvalue=list
            list.ylog=0
            widget_control,ev.top,set_uvalue=list
            RefreshDisplayCHI,list
            widget_control,ev.id,set_value='Y: linear'
            endcase
'Y: linear': begin
            widget_control,ev.top,get_uvalue=list
            list.ylog=1
            list.ypower=0
            list.ypowerstr=''
            widget_control,ev.top,set_uvalue=list
            RefreshDisplayCHI,list
            widget_control,ev.id,set_value='Y: logarithmic'
            endcase
'X: logarithmic' :begin
            widget_control,ev.top,get_uvalue=list
            list.xlog=0
            widget_control,ev.top,set_uvalue=list
            RefreshDisplayCHI,list
            widget_control,ev.id,set_value='X: linear'
            endcase
'X: linear': begin
            widget_control,ev.top,get_uvalue=list
            list.xlog=1
            widget_control,ev.top,set_uvalue=list
            RefreshDisplayCHI,list
            widget_control,ev.id,set_value='X: logarithmic'
            endcase
'>':          begin
            widget_control,ev.id,get_uvalue=uval
            fpeaksearch,ev,1,uval
            endcase
'<':          begin
            widget_control,ev.id,get_uvalue=uval
            fpeaksearch,ev,-1,uval
               endcase
'Peak+1':     begin
            peaksearch,ev,1
            endcase
'Peak-1':     begin
            peaksearch,ev,-1
            endcase
'Peak+10':    begin
            peaksearch,ev,10
            endcase
'Peak-10':    begin
            peaksearch,ev,-10
            endcase
'Select Peak':begin
            modeCHI,ev,val
            endcase
'Select Any':begin
            modeCHI,ev,val
            endcase
'Refine':    begin
            modeCHI,ev,val
            endcase
'Background':begin
            modeCHI,ev,val
            endcase
'Show peaks':begin
            widget_control,ev.top,get_uvalue=list
            list.peaksflag=1
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list
            widget_control,ev.id,set_value='Hide peaks'
            endcase
'Show peakselect':begin
            widget_control,ev.top,get_uvalue=list
            list.peakselectflag=1
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list
            widget_control,ev.id,set_value='Hide peakselect'
;            ID5=widget_info(ev.top,find_by_uname='x-value')
;            str1=strtrim(string(format='(F15.3)',list.peakxen[list.peakselect]),2)+list.xunits[list.xtype]
;            str2=strtrim(string(format='(F15.3)',list.peaky[list.peakselect]),2)
;            if ptr_valid(list.speringpr) then str3=strtrim(string(format='(F15.3,"%")',(*list.speringpr)[list.peakxkan[list.peakselect]]),2) else str3=''
;            widget_control,ID5,set_value=[str1,str2,str3]
            UpdatePeakSelectText,ev,list
            endcase
'Recall'  : begin
            widget_control,ev.top,get_uvalue=list
            x=(*list.xvalt)[list.backranind[0]:list.backranind[1],list.xtypeorig]
            y=(*list.spet)[list.backranind[0]:list.backranind[1],list.xtypeorig]
            peaks, list,x=x,y=y,offx=list.backranind[0],offy=list.backranind[0]
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list
            RefreshWindowCHI,ev
            endcase
'Hide peaks': begin
            widget_control,ev.top,get_uvalue=list
            list.peaksflag=0
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list
            widget_control,ev.id,set_value='Show peaks'
            endcase
'Hide peakselect': begin
            widget_control,ev.top,get_uvalue=list
            list.peakselectflag=0
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list
            widget_control,ev.id,set_value='Show peakselect'
            endcase
'MLR(G) substr':begin;Multiple Linear Regression
            if ev.select eq 0 then return
            widget_control,ev.top,get_uvalue=list
            list.fitmethod=0
            widget_control,ev.top,set_uvalue=list
            endcase
'NLLS(G) ortpolc':begin
               if ev.select eq 0 then return
               widget_control,ev.top,get_uvalue=list
               list.fitmethod=1
               widget_control,ev.top,set_uvalue=list
               endcase
'NLLS(G) substr':begin
            if ev.select eq 0 then return
             widget_control,ev.top,get_uvalue=list
              list.fitmethod=2
             widget_control,ev.top,set_uvalue=list
                endcase
'NLLS(PVn) substr':begin
              if ev.select eq 0 then return
             widget_control,ev.top,get_uvalue=list
               list.fitmethod=3
               widget_control,ev.top,set_uvalue=list
               endcase
'NLLS(PVnarr) substr':begin
             if ev.select eq 0 then return
              widget_control,ev.top,get_uvalue=list
               list.fitmethod=5
              widget_control,ev.top,set_uvalue=list
              endcase
'Fit Spectrum':   begin
              widget_control,/hourglass
              widget_control,ev.top,get_uvalue=list
              FitScan,list
              RefreshDisplayCHI,list
             widget_control,ev.top,set_uvalue=list
             RefreshWindowCHI,ev
               endcase
'Reset Fit':begin
            WIDGET_CONTROL, ev.top, get_uvalue=list
            ptr_free,list.yfit
            list.apeaksout=0
            list.chisq=0
            list.iter=0
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list
            RefreshWindowCHI,ev
            endcase
'Save as Image...':begin
            WIDGET_CONTROL, ev.top, get_uvalue=list
            if n_elements(list) EQ 0 then return
            wset,list.drawindex
            res=CutPath(list.file,file=name)
            SaveImage,list.path,name+'.jpg','plotCHI',list,OutID=list.OutID
            endcase
'Subtract':begin
            WIDGET_CONTROL, ev.top, get_uvalue=list
            if not ptr_valid(list.yback) then return
            spe=(*list.spet)[*,list.xtype]
               BackSub,spe,*list.yback,list.backmod,offset=offset
               (*list.spet)[*,list.xtype]=spe
             list.backcor=1
             RefreshPDD_CHI,list
             WIDGET_CONTROL, ev.top, set_uvalue=list
             RefreshDisplayCHI,list
             ID=widget_info(ev.top,FIND_BY_UNAME='Restore')
              widget_control,ID,sensitive=1
               widget_control,ev.ID,sensitive=0
              endcase
'Restore':     begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
               (*list.spet)=(*list.speorg)
               (*list.speerrort)=(*list.speerrororg)
             list.backcor=0
             RefreshPDD_CHI,list
            WIDGET_CONTROL, ev.top, set_uvalue=list
               RefreshDisplayCHI,list
              ID=widget_info(ev.top,FIND_BY_UNAME='Subtract')
             widget_control,ID,sensitive=1
               widget_control,ev.ID,sensitive=0
              endcase
'Strip':    begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
               AutoCorrect1D,list
             RefreshDisplayCHI,list
              widget_control,ev.top,set_uvalue=list
             endcase
;'Save PKS': begin
;              WIDGET_CONTROL, ev.top, get_uvalue=list
;             if n_elements(list) EQ 0 then return
;              error=SavePKS(list,/oD)
;              endcase
'Save Mask...':begin
               WIDGET_CONTROL, ev.top, get_uvalue=list
            if n_elements(list) EQ 0 then return

              if n_tags(ev) ne 5 then SaveMask,list,1,/oD $
             else SaveMask,list,0,/oD

            widget_control,ev.top,set_uvalue=list
               endcase
'Load Mask':begin
             WIDGET_CONTROL, ev.top, get_uvalue=list
              LoadMask,list,1,list.OutID,error=error,/oD
              widget_control,ev.top,set_uvalue=list
              if widget_info(list.NewWindowID[0],/valid) then begin
                  widget_control,list.NewWindowID[0],/destroy
                  EditModel,ev.top
              endif
              RefreshDataCHI,ev
              endcase
'Add Mask':    begin
            WIDGET_CONTROL, ev.top, get_uvalue=list
              LoadMask,list,1,list.OutID,error=error,/oD,/flags,flagi=[$
                  '$FIT_ROI:','$BACKoD:','$FIL:',$
                  '$SOURCE ('+msymbols('angstroms')+'):','$PIXSIZE (m):','$DET-SAM (m):']
              widget_control,ev.top,set_uvalue=list
              if widget_info(list.NewWindowID[0],/valid) then begin
                  widget_control,list.NewWindowID[0],/destroy
                  EditModel,ev.top
              endif
              RefreshDataCHI,ev
            endcase
'Residual >>': begin
               WIDGET_CONTROL, ev.top, get_uvalue=list
               if n_elements(list) EQ 0 then return
               WIDGET_CONTROL, ev.id, get_uvalue=base
               draw = WIDGET_DRAW(base, XSIZE=list.xsize, YSIZE=list.ysize/2.5,uname='drawres')
              widget_control,draw,get_value=resindex
               list.resindex=resindex
               wset,list.resindex
              plotres,list
              WIDGET_CONTROL, ev.id, set_value='Residual <<'
              WIDGET_CONTROL, ev.top, set_uvalue=list
            endcase
'Residual <<': begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              if n_elements(list) EQ 0 then return
              ID=widget_info(ev.top,FIND_BY_UNAME='drawres')
              WIDGET_CONTROL,ID, /DESTROY
              list.resindex=-1L
            WIDGET_CONTROL, ev.id, set_value='Residual >>'
              WIDGET_CONTROL, ev.top, set_uvalue=list
              endcase
'Add Label':begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              if n_elements(list) EQ 0 then return
              x=list.anyselect[list.xtype]
            y=(*list.spet)[list.anyselect[4],list.xtype]
            
              str=strtrim(string(format='(F15.3)',x),2)+list.shortxunits[list.xtype]
             ;str=strtrim(string(format='(F15.3)',y),2)
              if list.alabels eq 0 then begin
                ptr_free,list.labels
                ptr_free,list.xlabels
                ptr_free,list.ylabels
                list.labels=ptr_new(str)
                list.ylabels=ptr_new(y)
                list.xlabels=ptr_new(CHI_xval(x,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtype))
            endif else begin
                (*list.labels)=[(*list.labels),str]
                (*list.ylabels)=[(*list.ylabels),y]
                (*list.xlabels)=[(*list.xlabels),CHI_xval(x,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtype)]
            endelse
            list.alabels++
              widget_control,ev.top,set_uvalue=list
              RefreshDisplayCHI,list
              endcase
'Delete label':begin
              widget_control,ev.top,get_uvalue=list
             if list.alabels eq 0 then begin
                  printw,list.OutID, 'No labels present'
                  return
               endif
             modeCHI,ev,val
              endcase
'Move label':begin
              widget_control,ev.top,get_uvalue=list
               if list.alabels eq 0 then begin
                  printw,list.OutID, 'No labels present'
                  return
               endif
               modeCHI,ev,val
              endcase
'Peaksearch before substr':begin
              widget_control,ev.top,get_uvalue=list
             list.filorder=1
              widget_control,ev.top,set_uvalue=list
              widget_control,ev.id,set_value='Peaksearch after substr'
              endcase
'Peaksearch after substr':begin
               widget_control,ev.top,get_uvalue=list
               list.filorder=0
               widget_control,ev.top,set_uvalue=list
              widget_control,ev.id,set_value='Peaksearch before substr'
               endcase
'Fit total':begin
              widget_control,ev.top,get_uvalue=list
               list.fitresolve=2
              widget_control,ev.top,set_uvalue=list
              RefreshDisplayCHI,list
               widget_control,ev.id,set_value='Fit resolved'
               endcase
'Fit resolved':begin
               widget_control,ev.top,get_uvalue=list
              list.fitresolve=0
            widget_control,ev.top,set_uvalue=list
               RefreshDisplayCHI,list
               widget_control,ev.id,set_value='Fit hide'
               endcase
'Fit hide':    begin
               widget_control,ev.top,get_uvalue=list
              list.fitresolve=1
            widget_control,ev.top,set_uvalue=list
               RefreshDisplayCHI,list
               widget_control,ev.id,set_value='Fit total'
               endcase
'Spectrum Line':begin
              widget_control,ev.top,get_uvalue=list
               list.psym=1
             widget_control,ev.top,set_uvalue=list
              RefreshDisplayCHI,list
               widget_control,ev.id,set_value='Spectrum Points'
              endcase
'Spectrum Points':  begin
              widget_control,ev.top,get_uvalue=list
               list.psym=0
               widget_control,ev.top,set_uvalue=list
               RefreshDisplayCHI,list
              widget_control,ev.id,set_value='Spectrum Line'
              endcase
'Hide error bars':begin
              widget_control,ev.top,get_uvalue=list
               list.byerror=1
             widget_control,ev.top,set_uvalue=list
              RefreshDisplayCHI,list
               widget_control,ev.id,set_value='Show error bars'
              endcase
'Show error bars':  begin
              widget_control,ev.top,get_uvalue=list
               list.byerror=0
               widget_control,ev.top,set_uvalue=list
               RefreshDisplayCHI,list
              widget_control,ev.id,set_value='Hide error bars'
              endcase
'Next file':begin
            widget_control,ev.top,get_uvalue=list
              pathi = file_search(list.path+'*'+list.Format)
              nfiles = n_elements(pathi)
              path = SortListFiles(pathi,list.sortmethod,list.sortseparator)
              ind = where(path eq list.path+list.file,count)
              if (count eq 0) or (ind[0] eq nfiles-1) then ind=0 else ind=ind+1
              path=path[ind]
              path=path[0]
              error=CutPath(path,file=file,ext=ext)
              list.file=file+ext
              widget_control,ev.top,set_uvalue=list
              OpenCHI,ev,0
              endcase
'Previous file':begin
              widget_control,ev.top,get_uvalue=list
              pathi = file_search(list.path+'*'+list.Format)
              nfiles = n_elements(pathi)
              path = SortListFiles(pathi,list.sortmethod,list.sortseparator)
              ind = where(path eq list.path+list.file,count)
              if (count eq 0) or (ind[0] eq 0) then ind=nfiles-1 else ind=ind-1
              path=path[ind]
              path=path[0]
              error=CutPath(path,file=file,ext=ext)
              list.file=file+ext
              widget_control,ev.top,set_uvalue=list
              OpenCHI,ev,0
              endcase
'SortOptions':begin
              MakeBOptionsWindow,ev
              endcase
'Parameters':begin
              ParamWindow,ev
              endcase
'PDF markers':begin
              modeCHI,ev,val
              endcase
'Load PDF':begin
              LoadPDD_CHI,ev,1
              endcase
'Load Multiple PDF':begin
              LoadPDD_CHI,ev,2
              endcase
'Keep Last Added':    begin
            widget_control,ev.top,get_uvalue=list
            list.ReplaceLastPDF=~ev.select
            widget_control,ev.top,set_uvalue=list
            endcase
'PDF >>':    begin
            widget_control,ev.top,get_uvalue=list
            b=FindNextPDD(list)
            if b then widget_control,ev.top,set_uvalue=list

              ID=widget_info(ev.top,FIND_BY_UNAME='Keep PDF')
            widget_control,ID,set_button=0
              LoadPDD_CHI,ev,0,search=1
              endcase
'PDF <<':    begin
            widget_control,ev.top,get_uvalue=list
            b=FindNextPDD(list,/rev)
            if b then widget_control,ev.top,set_uvalue=list

              ID=widget_info(ev.top,FIND_BY_UNAME='Keep PDF')
            widget_control,ID,set_button=0
              LoadPDD_CHI,ev,0,search=2
              endcase
'Toggle PDF':begin
            widget_control,ev.top,get_uvalue=list
            if list.PDFnTotal ne 0 then begin
                  if ptr_valid(list.PDFs) then (*list.PDFs)[*]=list.PDFsToggle
                  list.PDFsToggle++
                  list.PDFsToggle mod= 4
                  widget_control,ev.top,set_uvalue=list
                  RefreshDisplayCHI,list
              endif
              endcase
'View PDF':   begin
            widget_control,ev.top,get_uvalue=list
            if list.PDFnTotal ne 0 then $
              EditPDFWindow,ev,'RefreshDisplayCHI',od={PDFx:(*list.PDFx),PDFi2:(*list.PDFi2),PDFScale:(*list.PDFScale),PDFOffset:(*list.PDFOffset)}
              endcase
'X: D-spacing':begin
            NewUnits,ev
            endcase
'X: 2-theta':begin
            NewUnits,ev
            endcase
'X: Radial Distance':begin
            NewUnits,ev
            endcase
'X: Q-space':begin
            NewUnits,ev
            endcase
'Show Legend':    begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.showlegend=1
            widget_control,ev.id,set_value='Hide Legend'
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
'Hide Legend':    begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.showlegend=0
            widget_control,ev.id,set_value='Show Legend'
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
'Show Data range':    begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.plotROI=1
            widget_control,ev.id,set_value='Hide Data range'
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
'Hide Data range':    begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.plotROI=0
            widget_control,ev.id,set_value='Show Data range'
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
'Show Labels':    begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.sLabels=1
            widget_control,ev.id,set_value='Hide Labels'
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
'Hide Labels':    begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.sLabels=0
            widget_control,ev.id,set_value='Show Labels'
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
'Show ring% values':begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.sringpr=1
            widget_control,ev.id,set_value='Hide ring% values'
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
'Hide ring% values':begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.sringpr=0
            widget_control,ev.id,set_value='Show ring% values'
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
'Show background':    begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.sBackground=1
            widget_control,ev.id,set_value='Hide background'
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
'Hide background':    begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.sBackground=0
            widget_control,ev.id,set_value='Show background'
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
'Refresh Display':    begin
            widget_control,ev.top,get_uvalue=list
            RefreshDisplayCHI,list
            endcase
'Clear Display':    begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.fitresolve=0
            if ptr_valid(list.PDFs) then (*list.PDFs)[*]=0
            list.sLabels=0
            list.peaksflag=0
            list.peakselectflag=0
            list.sBackground=0
            list.showlegend=0
            list.plotROI=0
            if ptr_valid(list.smask) then (*list.smask)[*]=0
            RefreshDisplayCHI,list
            widget_control,ev.top,set_uvalue=list,/no_copy
            RefreshWindowCHI,ev
            endcase
'Edit Peaks':begin
            modeCHI,ev,val
            endcase
'Initial':    begin
            modeCHI,ev,val
            endcase
'Refined':    begin
            modeCHI,ev,val
            endcase
'Delete iPeak':    begin
            modeCHI,ev,val
            endcase
'Delete iSelect':begin
            modeCHI,ev,val
            RefreshWindowCHI,ev
            endcase
'Add iPeak':begin
            modeCHI,ev,val
            endcase
'Select iPeak':begin
            modeCHI,ev,val
            endcase
'Delete rSelect':begin
            modeCHI,ev,val
            RefreshWindowCHI,ev
            endcase
'View Mask':begin
              EditMaskWindow,ev,0,1
              endcase
'Edit Mask':begin
              EditMaskWindow,ev,1,1
              endcase
'Mask Setup':begin
            modeCHI,ev,val
            endcase
'Set':        begin
              widget_control,ev.top,get_uvalue=list
              if list.masknr ne 0 then begin
                  (*list.mask)[2,*]=1
                  widget_control,ev.top,set_uvalue=list
                  RefreshDisplayCHI,list
              endif
              endcase
'Model':    begin
            widget_control,ev.top,get_uvalue=list
            if widget_info(list.NewWindowID[0],/valid) then widget_control,list.NewWindowID[0],/show $
            else EditModel,ev.top
            endcase
'Convert profile...':begin
            widget_control,ev.top,get_uvalue=list
            BatchConvertCHI,ev.top,list
            endcase
'Save profile...':begin
            widget_control,ev.top,get_uvalue=list
            SaveSingleChi,list
            widget_control,ev.top,set_uvalue=list
            endcase
'Save fitted profile...': begin
            widget_control,ev.top,get_uvalue=list
            SavefitmodelFitResult,list,ev.top,/O
            endcase
'Generate 2D pattern from Fit': begin
            widget_control,ev.top,get_uvalue=list
            if ~widget_info(list.main,/valid) then return
            
            printw,list.OutID,'Start 2D pattern generation ...'
            
            widget_control,list.main,get_uvalue=list2,/hourglass
            if ptr_valid(list2.tiffs) then begin
                imgsize=(*list2.tiffs)[1:2]
                center=list2.center
                a=list2.a
                b=list2.b
                phipixphi0=list2.phipixphi0
                IntCor2Dinfo=list2.IntCor2Dinfo
            endif else begin
                imgsize=[2048,2048]
                center=imgsize/2.
                a=list.a
                b=list.b
                phipixphi0=list.phipixphi0
                IntCor2Dinfo=list.IntCor2Dinfo
            endelse
            
            res=BraggXYtoX_I(rebin(lindgen(imgsize[0]),imgsize,/sample),rebin(lindgen(1,imgsize[1]),imgsize,/sample),1.,$
                            list.lambda,list.dist,list.scrx,list.scry,a,b,phipixphi0,center,IntCor2Dinfo,/onlyt,/angledeg,/pixm)
            
            img=reform( fitmodel2D(res.x1,list,ev.top,background)/res.(1),imgsize)
            ind=where(~finite(img),ct)
            if ct ne 0 then img[ind]=0

            if keyword_set(background) then begin
                p=(*(*list.ptrModel).ptrFit).yb
                if ptr_valid(p) then begin
                    img+=inverseAzIntWarp(1,(*list.xvalt)[list.backranind[0]:list.backranind[1],1],*p,imgsize,$
                                list.lambda,list.dist,list.scrx,list.scry,a,b,phipixphi0,center,IntCor2Dinfo)
                endif
            endif
            
            printw,list.OutID,'Min/Max = '+stringr(min(img))+"  / "+stringr(max(img))
            
            ; For testing
            ;img[*]=100
            
;            x={x:lindgen(imgsize[0]),y:lindgen(imgsize[1])}
;            p=[[1024,1700],$
;                [1259,1660],$
;                [1480,1529],$
;                [1656,1275],$
;                [1697,927],$
;                [1615,689],$
;                [1413,476],$
;                [1222,370],$
;                [935,349],$
;                [673,439],$
;                [480,615],$
;                [370,845],$
;                [353,1123],$
;                [427,1349],$
;                [583,1541]]
;            i=15
;            gfunc2dg,x,[1000.^2,p[0,i],p[1,i],20,20,1,100],f
;            img=reform(f,imgsize)
            
            
            ;AddNoise,img
            
            ;if ptr_valid(list2.tiffs) then $
            ;    img=TypeConvert(img,(*list2.tiffs)[3])
            
            MainUpdateCManual,list.path+list.file+'.tif',list.main,list.ev,fileptr=img
            endcase
'Generate 2D pattern from Data':begin
            widget_control,ev.top,get_uvalue=list
            if ~widget_info(list.main,/valid) then return
            
            printw,list.OutID,'Start 2D pattern generation ...'
            
            widget_control,list.main,get_uvalue=list2,/hourglass
            if ptr_valid(list2.tiffs) then begin
                imgsize=(*list2.tiffs)[1:2]
                center=list2.center
                a=list2.a
                b=list2.b
                phipixphi0=list2.phipixphi0
                IntCor2Dinfo=list2.IntCor2Dinfo
            endif else begin
                imgsize=[2048,2048]
                center=imgsize/2.
                a=list.a
                b=list.b
                phipixphi0=list.phipixphi0
                IntCor2Dinfo=list.IntCor2Dinfo
            endelse
            
            img=inverseAzIntWarp(list.xtypeorig,(*list.xvalt)[*,list.xtypeorig],(*list.spet)[*,list.xtypeorig],imgsize,$
                                list.lambda,list.dist,list.scrx,list.scry,a,b,phipixphi0,center,IntCor2Dinfo)

            printw,list.OutID,'Min/Max = '+stringr(min(img))+"  / "+stringr(max(img))
            ;if ptr_valid(list2.tiffs) then $
            ;    img=TypeConvert(img,(*list2.tiffs)[3])
                
            MainUpdateCManual,list.path+list.file+'.tif',list.main,list.ev,fileptr=img

            endcase
'Generate 2D pattern from cte':begin
            widget_control,ev.top,get_uvalue=list
            if ~widget_info(list.main,/valid) then return
            
            printw,list.OutID,'Start 2D pattern generation ...'
            
;            path=SelFile(list.path,list.file+'.tif','*'+fileformat(2),'Save CCD...')
;            if path EQ '' then return
            
            widget_control,list.main,get_uvalue=list2,/hourglass
            if ptr_valid(list2.tiffs) then begin
                imgsize=(*list2.tiffs)[1:2]
                center=list2.center
                a=list2.a
                b=list2.b
                phipixphi0=list2.phipixphi0
                IntCor2Dinfo=list2.IntCor2Dinfo
            endif else begin
                imgsize=[2048,2048]
                center=imgsize/2.
                a=list.a
                b=list.b
                phipixphi0=list.phipixphi0
                IntCor2Dinfo=list.IntCor2Dinfo
            endelse
            
            img=inverseAzIntWarp(list.xtype,(*list.xvalt)[*,list.xtype],1.,imgsize,$
                                list.lambda,list.dist,list.scrx,list.scry,a,b,phipixphi0,center,IntCor2Dinfo)

            printw,list.OutID,'Min/Max = '+stringr(min(img))+"  / "+stringr(max(img))
            ;if ptr_valid(list2.tiffs) then $
            ;    img=TypeConvert(img,(*list2.tiffs)[3])
                
;            struc={center:center,dist:list.dist,file:list.file,scrx:list.scrx,$
;                lambda:list.lambda,darkfile:''}
;            if WriteCCD(img,path,struc=struc) then printw,list.OutID,'Pattern saved: '+path

            MainUpdateCManual,list.path+list.file+'.tif',list.main,list.ev,fileptr=img

            endcase
              else :  return
    endcase
    endcase
9:    begin
    widget_control,ev.id,get_uvalue=uval
    case uval of
    'peaksin':    begin
                case (ev.type>2) of
                2:    begin ;Change cel content
                    selectold=widget_info(ev.id,/TABLE_SELECT)
                    widget_control,ev.id,SET_TABLE_SELECT=[ ev.x, ev.y, ev.x, ev.y ]
                    widget_control,ev.id,get_value=val,/USE_TABLE_SELECT
                    widget_control,ev.top,get_uvalue=list
                    case ev.x of
                    0:    list.peakxen[ev.y]=val
                    1:    list.peakarea[ev.y]=val
                    2:    list.peaksigma[ev.y]=val
                    else: return
                    endcase
                    widget_control,ev.top,set_uvalue=list
                    RefreshDisplayCHI,list
                    widget_control,ev.id,SET_TABLE_SELECT=selectold
                    endcase
                4:    begin ;Select
                       if ev.sel_top ne -1 then begin
                           widget_control,ev.top,get_uvalue=list
                        list.peakselect=ev.sel_top
                        widget_control,ev.top,set_uvalue=list
                        RefreshDisplayCHI,list

                        ID1=widget_info(ev.top,FIND_BY_UNAME='x-value')
                        xen=CHI_xval(list.peakxen[list.peakselect],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtypeorig)
                        str1=strtrim(string(format='(F15.3)',xen[list.xtype]),2)+list.xunits[list.xtype]
                        str2=strtrim(string(format='(F15.3)',list.peaky[list.peakselect]),2)
                        if ptr_valid(list.speringpr) then str3=strtrim(string(format='(F15.3,"%")',(*list.speringpr)[list.peakxkan[list.peakselect]]),2) else str3=''
                        widget_control,ID1,set_value=[str1,str2,str3]
                    endif
                    endcase
                else:
                endcase
                endcase
    'peaksout':    begin
                case [ev.type] of
                4:    begin ;Select
                       if ev.sel_top ne -1 then begin
                           widget_control,ev.top,get_uvalue=list
                        list.fitselect=ev.sel_top
                        widget_control,ev.top,set_uvalue=list
                        RefreshDisplayCHI,list
                    endif
                    endcase
                else:
                endcase
                endcase
    endcase
    endcase
10:    begin
    widget_control,ev.id,get_uvalue=uval
    case uval of
    '1Dback':begin
            WIDGET_CONTROL, ev.top, get_uvalue=list,/no_copy
            if ev.tab eq 0 then list.backmod=-1 else list.backmod=ev.tab
             WIDGET_CONTROL, ev.top, set_uvalue=list,/no_copy
            endcase
    'PDFsearch':begin
            WIDGET_CONTROL, ev.top, get_uvalue=list,/no_copy
            list.PDFSearch.type=ev.tab
             WIDGET_CONTROL, ev.top, set_uvalue=list,/no_copy
            endcase
    endcase
    endcase
else:
endcase
end ;pro XRD_CHI_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRD_CHI,ev=ev,retID=Topbase
print,'=========================================START 1D Profile Editor'

if keyword_set(ev) then begin
    widget_control,ev.top,get_uvalue=list
    path=list.chipath
    file=list.chifile
    mskpath=list.mskpath
    mskfile=list.mskfile
    pddpath=list.pddpath
    pddfile=list.pddfile
    Platform=list.Platform
    main=ev.top
    lambda=list.lambda
    cconv=list.cconv
    scrx=list.scrx
    scry=list.scry
    dist=list.dist
    a=list.a
    b=list.b
    phipixphi0=list.phipixphi0
    IntCor2Dinfo=list.IntCor2Dinfo
    sortseparator=list.sortseparator
    sortmethod=list.sortmethod
endif else begin
    ControlDevice,Platform,Visual,Msize
    result=ReadINI('$XRD_CHI path:')
    path=result.(0)
    file=result.(1)
    result=ReadINI('$XRD mskpath:')
    mskpath=result.(0)
    mskfile=result.(1)
    result=ReadINI('$XRD pddpath:')
    pddpath=result.(0)
    pddfile=result.(1)
    main=0L
    lambda=0.
    cconv=double(4.13566743*299792458.*10.^(-8))
    scrx=0.
    scry=0.
    dist=0.
    a=0.
    b=0.
    phipixphi0=0.
    IntCor2Dinfo=DefaultIntCor2Dinfo()
    sortseparator='_'
    sortmethod=1B
    ev={top:0l,id:0l}
endelse

result=ReadINI('$CrossUpdate:')
CrossUpdate=result.(0)

error=CutPath(file,ext=ext)
Format=ext

br=15
br1=25

result=ReadINI('$ChiSize:')
xsize=(result.(0))[0]
ysize=(result.(0))[1]
Topbase = WIDGET_BASE(MBAR=bar,/row,title='1D Profile Editor:')
    base1=widget_base(Topbase,/column,sensitive=0,uname='base1')
       label1= widget_label(base1,value='BROWSE:')
       button21= widget_button(base1,value='Next file')
       button22= widget_button(base1,value='Previous file')

       label1= widget_label(base1,value='MODE:')
       button1= widget_button(base1,value='Data range',uname='Data range',xsize=br)
       button2= widget_button(base1,value='Zoom',uname='Zoom',xsize=br)
       button11= widget_button(base1,value='Mask Setup',uname='Mask Setup',xsize=br)
       ;button33=widget_button(base1,value='Edit Peaks',uname='Edit Peaks',xsize=br)
       ;button37=widget_button(base1,value='Refine',uname='Refine',xsize=br)
       button37=widget_button(base1,value='Background',uname='Background',xsize=br)
       button37=widget_button(base1,value='PDF markers',uname='PDF markers',xsize=br)
       button37=widget_button(base1,value='Model',uname='Model',xsize=br)
       button37=widget_button(base1,value='Move Legend',uname='Move Legend',xsize=br)

       label3= widget_label(base1,value='LABELS:')
       base3=widget_base(base1,/column)
             button30=widget_button(base3,value='Select Any',uname='Select Any',xsize=br)
          text1=widget_text(base3,uname='x-value2',xsize=br,ysize=4)
             ;button30=widget_button(base3,value='Select Peak',uname='Select Peak',xsize=br)
          ;text1=widget_text(base3,uname='x-value',xsize=br,ysize=3)
          text8=widget_text(base3,/editable,xsize=10,uvalue='ulab')
          button43=widget_button(base3,value='Add Label',uname='Add Label',xsize=br)
          button44=widget_button(base3,value='Delete label',uname='Delete label',xsize=br)
             button37=widget_button(base3,value='Move label',uname='Move label',xsize=br)
          

    base2=widget_base(Topbase,/column,sensitive=0,uname='base2')
       draw = WIDGET_DRAW(base2, XSIZE=xsize, YSIZE=ysize,uname='draw')
       Window, /Free, xsize=xsize, ysize=ysize, /Pixmap
       pixindex = !D.Window
       button42=widget_button(base2,value='Residual >>' ,uname='Expand',xsize=br1,uvalue=base2)
    base5=widget_base(Topbase,/column,uname='base5')

OutID=SubscribeToOutID()

tempP=(ReadINI('$#Peaks:')).(0)
far500=fltarr(tempP)

ptrModel=CreateModel()

list={listtype:10,$            ; list type
    xvalt:PTR_NEW(),$        ; d-values, two-theta and radial distance
    nxval:0L,$                ; number of bins
    spet:PTR_NEW(),$        ; spectrum data
    yfit:PTR_NEW(),$        ; fitted spectrum
    speorg:PTR_NEW(),$      ; spectrum data original
    speerrort:PTR_NEW(),$   ; spectrum standard deviation
    speerrororg:PTR_NEW(),$   ; spectrum standard deviation    
    speringpr:PTR_NEW(),$    ; ring percentage for each point in spe
    ypower:0.,$                ; intensity^power for plotting
    ypowerstr:'',$            ; power string
    Format:Format,$         ; contains extension of file to read
    Platform:Platform,$     ; OS Platform specifications
    xtype:1,$                ; 0=> d, 1=> two-theta, 2=> radial distance, 3=> Q-space
    xtypeorig:1,$            ; Original x-type
    main:main,$                ; Main ID
    ev:ev,$                    ; Button that generated this window
    CrossUpdate:CrossUpdate,$; Update Main?
    ;
    masknr:0,$                ; number of areas in mask
    mask:PTR_NEW(),$        ; first column: maskpart type (1:ROI)
                             ; second column: index of first part of a linked group of parts
                               ; third column: 1=set, 0=not set
                               ; rest: parameters of maskpart
                               ;     left,right coord in Angstroms
    smask:PTR_NEW(),$        ; show mask areas
    MaskP:-1L,$                ; Highlight this mask group
    ;
    path:path,$             ; path of d-file(vb.D:\wout\files2\)
    file:file,$             ; name of d-file(vb.ca24.d)
    mskfile:mskfile,$       ; name of .in file (vb.ca24.in)
    mskpath:mskpath,$       ; path where all files are saved and loaded(vb.D:\wout\files2\)
    pddpath:pddpath,$       ; Path of .pdd file
    pddfile:pddfile,$       ; Powder diffraction markers from .pdd file
    ;
    sortseparator:sortseparator,$ ; for sorting files
    sortmethod:sortmethod,$ ; 1 => sort strings, 0 sort strings with offset
    ;
    asciitemplate:{    rowskip:2l,$
                    delimiter:string(32b)+string(9b),$
                    ind:[1,2l]},$
    ;
    lambda:lambda,$         ; X-ray lambda(Angstroms)
    cconv:cconv,$            ; E(keV)=12.3984/lambda(Angstroms)
    scrx:scrx,$                ; Pixel x-size (m, input in um)
    scry:scry,$                ; Pixel y-size (m, input in um)
    dist:dist,$                ; Distance sample-detector (m, input in cm)
    a:a,$                    ; a Tilt angle: RotXc(a)
    b:b,$                    ; b Rotation angle: RotZt(b)
    phipixphi0:phipixphi0,$ ; azimuthal shift due to 2-angle calibration
    IntCor2Dinfo:IntCor2Dinfo,$; 2D intensity correction parameters
    ;
    labels:ptr_new(),$      ; names of added labels
    xlabels:ptr_new(),$        ; x-values of labels in different units
    ylabels:ptr_new(),$     ; y-values of added labels ()
    alabels:0,$             ; number of added labels
    DelMass:0b,$            ; Prompt when deleting multiple labels or not
    DelWidth:[15,5],$        ; Half width of label detection
    sLabels:1b,$            ; Show/hide labels
    xpos:0.5,$                ; Legend x-pos
    ypos:0.5,$                ; Legend y-pos
    incpos:0.03,$            ; Legend increment position (normal)
    showlegend:1b,$            ; Show legend
    plotROI:1b,$            ; Show ROI lines
    ;
    PDFn:ptr_new(),$        ; number of d-spacings in PDF file
    PDFd:ptr_new(),$        ; d-spacings in PDF file
    PDFx:ptr_new(),$        ; d, tt and radial distance
    PDFni:ptr_new(),$       ; number of intensities in PDF file
    PDFi:ptr_new(),$        ; intensities in PDF file scaled from 0-100
    PDFi2:ptr_new(),$        ; position intensities
    PDFhkl:ptr_new(),$        ; miller indices
    PDFname:ptr_new(),$     ; name of loaded compount
    PDFs:ptr_new(),$        ; first bit: show lines or not, second bit: show labels or not
    PDFsToggle:0,$            ; cfr PDFs
    nPDF:0,$                ; Number of PDF files loaded
    PDFnTotal:0,$            ; Total number of reflections loaded
    PDFScale:ptr_new(),$    ; Use spe values when not in [0,1]
    PDFOffset:ptr_new(),$    ; Start of PDF lines
    ReplaceLastPDF:1b,$        ; Replace last PDF or Replace all?
    PDFSearch:{type:0,PDFSubstr:'',PDFdi:0.,PDFde:0.},$    ; PDF subsearch
    BEXpos:-1L,$            ; Current BEX cards when browsing PDFs
    BEXeof:1b,$                ; 1: BEXpos last card
    ;
    psym:0,$                 ; Psym for spectrum plot
    byerror:0b,$            ; Yerror bars
    title:'',$               ; title of window
    xtitle:'',$             ; title of x-axis
    xmodes:ChiXtitle(mode=1),$;Possible x-titles
    xunits:ChiXtitle(mode=2),$    ;Possible x units
    shortxunits:ChiXtitle(mode=3),$    ;Possible x units
    ytitle:'',$             ; title of y-axis
    xrange:lonarr(2), $     ; The current X range of the plot.
    yrange:lonarr(2), $     ; The current Y range of the plot.
    xs:0L, $                ; X static corner of the zoom box.
    ys:0L, $                ; Y static corner of the zoom box.
    xd:0L, $                ; X dynamic corner of the zoom box.
    yd:0L,$                 ; Y dynamic corner of the zoom box.
    sysvar:PointSysVar(),$    ; System variable for coord. convertion
    xsize:xsize,$            ; draw widget xsize
    ysize:ysize,$            ; draw widget ysize
    draw:draw,$             ; ID van drawwidget
    drawindex:-1L,$           ; index van drawwidget
    pixindex:pixindex,$     ; index van pixwindow
    resindex:-1L,$             ; index van 2nd drawwidget
    OutID:OutID,$            ; output log
    NewWindowID:lonarr(1),$    ; id's of windows that might be open: [EditModel]
    ;
    peakout:fltarr(4,tempP),$ ; fit results for output: positions, area, sigma, n
    peakouts:fltarr(4,tempP),$; SD on param
    apeaksout:0,$           ; Number of peaks in fitting model
    chisq:0.,$              ; chisquare van fit
    iter:0L,$                ; Number of iterations for fit
    nfree:0,$                ; n_elements(y in ROI)-n_elements(fitparameters)
    fitmethod:5,$           ; fitmethod
    fitresolve:1,$             ; 0: plot fit, 1: plot seperate peaks
    ;
    tree:spacegrouptree(),$    ; Spacegroup tree
    ptrModel:ptrModel,$        ; fitmodel
    scatyieldstruct:{emin:1d,emax:100d,einc:0.1d,tmin:0.001d,tmax:1d,tinc:0.001d,energy:20d,thickness:1d,fraction:0d},$
    ;peaky:far500,$          ; y-values of found peaks in spectrum
    ;peakxen:far500,$           ; x-values of found peaks in spectrum
    ;peakxkan:intarr(tempP),$; x-values of found peaks in spectrum
    ;peaksigma:far500,$      ; sigma of found peaks in spectrum
    ;peakarea:far500,$       ; area of found peaks in spectrum
    ;peakbegin:intarr(tempP),$ ; indeces in list.en/spe of begin of found peaks in spectrum
    ;peakend:intarr(tempP),$ ; indeces in list.en/spe of end of found peaks in spectrum
    apeaks:0,$              ; number of found peaks in spectrum
    kDisp:0,$                ; k times sigma for display of peaks
    ;
    peakselect:0,$             ; index in list.peakx/... of selected peak
    anyselect:fltarr(5),$   ; index in list.xval/... of selected peak
    fitselect:0,$           ; index in list.peakout of selected peak (fitted)
    peaksflag:0b,$           ; 0->don't show peaks
                             ; 1->show peaks
    peakselectflag:1b,$        ; Show peak selector
    vin:10,$                   ; param of top-hat filter
    win:30,$                   ; param of top-hat filter
    r1:0.5,$                   ; param of top-hat filtering
    filorder:0b,$             ; find peaks before (0) or after (1) subtraction
    dragdrop:'',$              ; name of the label we're moving
    highlightpeak:{ptr:ptr_new(),i:0l},$; highlighted peak
    ;
    bPDFpsym:0b,$            ; PDF with psym
    blegendshort:1b,$        ; Short version of legend
    letter:1.,$                ; size of labels
    backcor:0B,$              ; Background subtracted (1B) or not (0B)
    yback:ptr_new(),$       ; backgroundvalue ()
                             ; this is the same outside backranind as list.spe
    backmod:-1,$             ; 0=>ortpol subtracted (Gaussian)
                             ; 1=>ortpol and fit its coeff. (Gaussian)
                             ; 2=>stripping(always substracted) (PV)
                             ; -1=>don't substract
    sBackground:1b,$        ; show/hide background
    sringpr:1b,$            ; show/hide ring percentages
    backranind:lonarr(2),$     ; indices in list.spe/yback in with background is calculated
    backranindD:fltarr(2),$    ; d-spacing for backranind
    ROIAzimuth:fltarr(2),$    ; azimuth over which the 1D pattern was integrated
    xlog:0b,$                ; x-axis log?
    ylog:0b}                ; 0->y-as linear / 1->y-as logarit.

    menu1 = WIDGET_BUTTON(bar, VALUE='File', /MENU)
       buttonopen = WIDGET_BUTTON(menu1, VALUE='Open...',uname='Open...')
       button = WIDGET_BUTTON(menu1, VALUE='Open dummy',uname='Open dummy')
       buttun14= widget_button(menu1,value='Save...',uname='Save...',/menu,sensitive=0)
            button39=widget_button(buttun14,value='Save profile...')
            button39=widget_button(buttun14,value='Save fitted profile...')
         button39=widget_button(buttun14,value='Save Mask...',uname='Save Mask...')
         ;button41=widget_button(buttun14,value='Save PKS',uname='Save PKS')
         button15=widget_button(buttun14,value='Save as Image...',uname='Save as Image...')
         button17=widget_button(buttun14,value='Convert profile...',uname='Convert profile...')
         
       button18= widget_button(menu1,value='Load...',uname='Load...',/menu,sensitive=0)
         button40=widget_button(button18,value='Load Mask',uname='Load Mask')
         button40=widget_button(button18,value='Add Mask',uname='Add Mask')
         button40=widget_button(button18,value='Load PDF',uname='Load PDF')
         button40=widget_button(button18,value='Load Multiple PDF',uname='Load Multiple PDF')
       button19 = WIDGET_BUTTON(menu1, VALUE='Exit',uname='Exit')

    menu2=WIDGET_BUTTON(bar, VALUE='Options', /MENU,uname='MOptions',sensitive=0)
        button3= widget_button(menu2,value='X: D-spacing ('+msymbols('angstroms')+')',uname='Xunit')
        button3= widget_button(menu2,value='X: linear',uname='X')
        button3= widget_button(menu2,value='Y: linear',uname='Y')
        button3= widget_button(menu2,value='SortOptions')
        button3= widget_button(menu2,value='Parameters')

    menu2=WIDGET_BUTTON(bar, VALUE='Perform', /MENU,uname='Perform',sensitive=0)
        button17=widget_button(menu2,value='Generate 2D pattern from Data',uname='Generate 2D pattern from Data')
        button17=widget_button(menu2,value='Generate 2D pattern from Fit',uname='Generate 2D pattern from Fit')
        button17=widget_button(menu2,value='Generate 2D pattern from cte',uname='Generate 2D pattern from cte')
         

    menu2=WIDGET_BUTTON(bar, VALUE='Display', /MENU,uname='MDisplay',sensitive=0)
       button12=widget_button(menu2,value='Fit hide',uname='PlotFit')
       button12=widget_button(menu2,value='Spectrum Line',uname='Spectrum Line')
       button12=widget_button(menu2,value='Hide error bars',uname='Hide error bars')
       button12=widget_button(menu2,value='View PDF',uname='View PDF')
       button10= widget_button(menu2,value='View Mask')
       button12=widget_button(menu2,value='Show Data range',uname='sROI')
       button12=widget_button(menu2,value='Show Legend',uname='Legend')
       button12=widget_button(menu2,value='Show Labels',uname='Labels')
       button12=widget_button(menu2,value='Show peaks',uname='peaks')
       button12=widget_button(menu2,value='Show peakselect',uname='peakselect')
       button12=widget_button(menu2,value='Show background',uname='view background')
       button12=widget_button(menu2,value='Show ring% values',uname='ringpr')
       button12=widget_button(menu2,value='Toggle PDF',ACCELERATOR = "F1")
       button12=widget_button(menu2,value='Refresh Display',ACCELERATOR = "F5")
       button12=widget_button(menu2,value='Clear Display')


WIDGET_CONTROL,Topbase, /REALIZE,set_Uvalue=list
Xmanager,'XRD_CHI',Topbase,/NO_BLOCK,cleanup='CleanUp_XRD_CHI'
widget_control,buttonopen,send_event={ID:buttonopen,TOP:Topbase,HANDLER:Topbase,SELECT:long(1)}
end ;pro XRD_CHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
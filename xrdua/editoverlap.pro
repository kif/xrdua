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

pro IMG_PROCESS_EVENTS, ev

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
tmp=DynTR([ev.x,ev.y],1./list.nmag,0); col,row
tmp[1]=list.nProf-1-tmp[1]
tmp=long(tmp)

ID=widget_info(ev.top,FIND_BY_UNAME='label3')
widget_control,ID,set_value='Last motion position: COL'+stringr(tmp[0])+'-ROW'+stringr(tmp[1])

IF thisEvent NE 'DOWN' THEN RETURN

list.x=tmp[0]
list.y=tmp[1]

widget_control,list.TableID[0],SET_TABLE_SELECT=[list.x,list.y,list.x,list.y]
widget_control,list.TableID[0],set_table_view=[list.x,list.y]
widget_control,list.TableID[1],SET_TABLE_SELECT=[list.x,list.y,list.x,list.y]
widget_control,list.TableID[1],set_table_view=[list.x,list.y]
WIDGET_CONTROL, ev.top, SET_UVALUE = list
RefreshEditOverlap,ev.top

end;pro IMG_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nFracDigits ,str,bE=bE,Form=Form

nAfter=0
posE=strpos(str , 'E')
posdec=strpos(str , '.')
nlen=strlen(str)
bE=posE eq -1
if bE then begin ; No scientific notation
    Form='F'
    if (posdec ne -1) then nAfter=nAfter+nlen-1-posdec $
    else begin
        i=nlen-1
        pos=strpos(str,'0',i)
        while i eq pos do begin
            nAfter=nAfter-1
            i=i-1
            pos=strpos(str,'0',i)
        endwhile
    endelse
    if (posdec eq nlen-1) then str=strmid(str,0,nlen-1); remove point
endif else begin ; Scientific notation
    Form='E'
    E=long(strmid(str,posE+1))
    nAfter=nAfter-E
    if (posdec ne -1) then nAfter=nAfter+posE-1-posdec
endelse

return,nAfter

end;function nFracDigits
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FormatnFracDigits ,xin,nAfter,Form,bE,strrep

if nAfter gt 0 then begin
    if bE eq 0 then nAfter=strlen(stringr(long(xin)))-1+nAfter
    if nAfter lt 0 then out=strrep else $
    out=string(xin,format='('+Form+'255.'+stringr(nAfter)+')')
endif else begin
    if nAfter ne 0 then begin
        exp10=10.^nAfter
        tmp=round(xin*exp10)/exp10
        nAfter=strlen(stringr(long(xin)))-1+nAfter
        if nAfter lt 0 then out=strrep else $
        out=string(tmp,format='('+Form+'255.'+stringr(nAfter)+')')
    endif else begin
        out=string(xin,format='('+Form+'255.0)')
        if bE then out=strmid(out,0,strlen(out)-1); remove point
    endelse
endelse

return,out
end;function FormatnFracDigits
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RoundWithSD ,xin,sxin

n=n_elements(xin)

; 1. Round standard deviation to one significant digit
;
; # significant bits = start from the first non-zero digit and end at the last digit
; Exception: e.g. 100
;    This has at least 1 significant bit, but could have 2 or 3
;    For this, use scientific notation: 1E+2 has 1 significant digit
;    Since 100E+0 has three significant digits, we might say that 100 has also three,
;    but this is just a personal convention.

SXstr=string(sxin,format='(G255.1)')
Xstr=SXstr

; 2. Round x so that it has the same number of digits
;    after the decimal point as the standard deviation
for i=0l,n-1 do begin

    ; How many digits after point?
    str=SXstr[i]
    nAfter=nFracDigits (str,Form=Form,bE=bE)
    SXstr[i]=str

    ; Format x
    Xstr[i]= FormatnFracDigits (xin[i],nAfter,Form,bE,SXstr[i])

endfor

Xstr=strcompress(Xstr,/remove_all)
SXstr=strcompress(SXstr,/remove_all)
return,{X:Xstr,SX:SXstr}
end;function RoundWithSD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PlotGroupProfile,ev

WIDGET_CONTROL, ev.top, GET_UVALUE = list
if widget_info(list.TableID[0],/VALID_ID) eq 0 then return

g=list.maskgroup[list.x]
ind=where(list.maskgroup eq g,ct)

if ct eq 1 then begin
    title='GROUP'+stringr(long(list.TableMean[0,g]))+' contains area '+stringr(list.IndMask[ind[0]])
    plot,(*list.Profiles[ind[0]])[*,0],$
        (*list.Profiles[ind[0]])[*,1],title=title
endif else begin
    title='GROUP'+stringr(long(list.TableMean[0,g]))+' contains areas: '
    for i=0l,ct-1 do begin
        xmi=min((*list.Profiles[ind[i]])[*,0],max=xma)
        ymi=min((*list.Profiles[ind[i]])[*,1],max=yma)
        if n_elements(b) ne 0 then begin
            bx=bx<xmi
            ex=ex>xma
            by=by<ymi
            ey=ey>yma
        endif else begin
            bx=xmi
            ex=xma
            by=ymi
            ey=yma
        endelse
        title=title+' ,'+stringr(list.IndMask[ind[i]])
    endfor

    i=0
    plot,(*list.Profiles[ind[i]])[*,0],$
        (*list.Profiles[ind[i]])[*,1],xrange=[bx,ex],yrange=[by,ey],title=title
    for i=1,ct-1 do begin
        oplot,(*list.Profiles[ind[i]])[*,0],$
        (*list.Profiles[ind[i]])[*,1]
    endfor
    bot=list.TableMean[1,g]-3*list.TableMean[2,g]
    up=list.TableMean[1,g]+3*list.TableMean[2,g]
    plots,[bot,bot],[!Y.crange[0],!Y.crange[1]]
    plots,[up,up],[!Y.crange[0],!Y.crange[1]]
endelse

end;pro PlotGroupProfile
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PlotProfiles,list,top
    tvlct,RoldCT,GoldCT,BoldCT,/get

    ; ---- Extract 2 profiles ----
    ProfRow=*list.Profiles[list.y]
    ProfCol=*list.Profiles[list.x]
    x1=ProfRow[*,0]
    n1=n_elements(x1)
    x2=ProfCol[*,0]
    n2=n_elements(x2)

    ; ---- Normalize ----
    xmi1=min(x1,max=xma1)
    xmi2=min(x2,max=xma2)
    ymi1=min(ProfRow[*,1],max=yma1)
    ymi2=min(ProfCol[*,1],max=yma2)
    xmi=xmi1<xmi2
    xma=xma1>xma2
    y1=ProfRow[*,1]/yma1
    y2=ProfCol[*,1]/yma2

    ; ---- Group info ----
    g1=list.maskgroup[list.y]
    g2=list.maskgroup[list.x]
    out=RoundWithSD([list.TableMean[1,g1],list.TableMean[1,g2]],[list.TableMean[2,g1],list.TableMean[2,g2]])
    m=out.(0)
    s=out.(1)

    ID=widget_info(top,FIND_BY_UNAME='label1')
    widget_control,ID,set_value='OVERLAP ROW'+stringr(list.y)+' (AREA'+stringr(list.IndMask[list.y])+') from GROUP'+stringr(long(list.TableMean[0,g1]))+': '+m[0]+msymbols('plus_minus')+s[0]+msymbols('angstroms')
    ID=widget_info(top,FIND_BY_UNAME='label2')
    widget_control,ID,set_value='BY COLUMN'+stringr(list.x)+' (AREA'+stringr(list.IndMask[list.x])+') from GROUP'+stringr(long(list.TableMean[0,g2]))+': '+m[1]+msymbols('plus_minus')+s[1]+msymbols('angstroms')

    ; ---- Plot ----
    LoadctRev,-1,/silent
    col=10*16
    plot,x1,y1,xrange=[xmi,xma],yrange=[-0.1,1.1],xstyle=1,ystyle=1,$
    title='Overlap two profiles (boarders:d'+msymbols('plus_minus')+'3.SD)',xtitle=ChiXtitle(0),ytitle='Intensity (a.u.)',col=col
    pxval=[x1[0],x1,x1[n1-1]]
    pyval=[-0.1 ,y1,-0.1 ]
    polyfill,pxval,pyval,/line_fill,col=col,/data,LINESTYLE=1,SPACING=0.05

    XYOUTS, 10,15, 'OVERLAP ROW'+stringr(list.y),/device,col=col

    LoadctRev,3,/silent
    col=10*12
    oplot,x2,y2,col=col
    pxval=[x2[0],x2,x2[n2-1]]
    pyval=[-0.1 ,y2,-0.1 ]
    polyfill,pxval,pyval,/line_fill,orientation=45,col=col,/data

    XYOUTS, 10,5, 'BY COLUMN'+stringr(list.x),/device,col=col

    LoadctRev,0,/silent,rev=-1
    bot=list.TableMean[1,g2]-3*list.TableMean[2,g2]
    up=list.TableMean[1,g2]+3*list.TableMean[2,g2]
    plots,[bot,bot],[!Y.crange[0],!Y.crange[1]]
    plots,[up,up],[!Y.crange[0],!Y.crange[1]]

    tvlct,RoldCT,GoldCT,BoldCT

end;pro PlotProfiles,list
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshEditOverlap,top

Widget_Control, top, Get_UValue=list

if widget_info(list.draw,/VALID_ID) then begin
    widget_control,list.draw,get_value=drawindex
    Wset,drawindex
    PlotProfiles,list,top

endif

end;pro RefreshEditOverlap
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION EditOverlapField, base, val, NAME = NAME, $
                        X_SCROLL_SIZE = X_SCROLL_SIZE, $
                        Y_SCROLL_SIZE = Y_SCROLL_SIZE

; ---- Is val a 2D array? ----
varsize = size(val)                 ;determine the size and
vardims = N_ELEMENTS(varsize) - 2           ;type of the variable
numelements = varsize[vardims + 1]
ID=0

IF (numelements le 1) or (varsize[0] NE 1 AND varsize[0] NE 2) THEN return,0L ; No array

abase = WIDGET_BASE(base, /FRAME, /COLUMN, XPAD = 8, YPAD = 8)

; ---- Title of array ----
IF(numelements GT 1) THEN BEGIN
  suffix = " ["
  FOR j = 1, varsize[0] DO BEGIN            ;
    suffix = suffix + strtrim(varsize[j], 2)        ;
    IF j NE varsize[0] THEN suffix = suffix + ", "
  ENDFOR
  suffix = suffix+"]"
ENDIF ELSE suffix = ""

IF(KEYWORD_SET(NAME)) THEN $
      lbl = WIDGET_LABEL(abase, VALUE = NAME + suffix) $
ELSE lbl = WIDGET_LABEL(abase, value = suffix)

; ---- Make Table ----
IF(N_ELEMENTS(X_SCROLL_SIZE) EQ 0) THEN $
  XSCROLL_SIZE = 4 ELSE XSCROLL_SIZE = X_SCROLL_SIZE
IF(N_ELEMENTS(Y_SCROLL_SIZE) EQ 0) THEN $
  YSCROLL_SIZE = 4 ELSE YSCROLL_SIZE = Y_SCROLL_SIZE
table = WIDGET_TABLE(abase, value = val, /ALL_EVENTS,$
                             X_SCROLL_SIZE = XSCROLL_SIZE, $
                             Y_SCROLL_SIZE = YSCROLL_SIZE)

RETURN, table

END ;FUNCTION EditOverlapField
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_EditOverlap, ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
PTR_FREE,list.Profiles
end;pro CleanUp_EditOverlap
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PRO EditOverlap_event, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
0 : begin ;buttongroup event

    CASE ev.value OF
        0: BEGIN                                            ;the user chose the
            WIDGET_CONTROL, ev.top, /DESTROY     ;return the initial
          ENDCASE                       ;variable
        1: BEGIN                        ;read another variable
            WIDGET_CONTROL, ev.top, GET_UVALUE = list,/NO_COPY

            path=list.path
            file=list.file
            if ReadOverlapFile(path,file,result=result) then begin
                ; CLean up
                if widget_info(list.ID,/VALID_ID) then $
                    WIDGET_CONTROL, list.ID,/destroy
                PTR_FREE,list.Profiles

                ; New table
                ID = WIDGET_BASE(list.dynID, /COLUMN)
                ID1 = EditOverlapField(ID, result.XY,NAME='Overlap XY')
                ID2 = EditOverlapField(ID, result.XbY,NAME='Overlap XRowByYCol')
                ID3 = widget_draw(ID,xsize=result.nProf*list.nmag,$
                    ysize=result.nProf*list.nmag,button=1,motion=1,$
                    Event_Pro='IMG_PROCESS_EVENTS')

                list={ID:ID,path:path,file:file,Profiles:result.Profiles,x:0L,y:0L,$
                draw:list.draw,dynID:list.dynID,TableID:[ID1,ID2,ID3],nmag:list.nmag,nProf:result.nProf, $
                nGroups:result.nGroups,maskgroup:result.maskgroup,TableMean:result.TableMean,IndMask:result.IndMask}

                widget_control,list.TableID[2],get_value=drawindex
                Wset,drawindex
                img=rebin(result.XY,list.nProf*list.nmag,list.nProf*list.nmag,/sample)
                LoadctRev,0,/silent,rev=-1
                tvscl,img,/order

            endif
            WIDGET_CONTROL, ev.top, SET_UVALUE = list
        ENDCASE
        2:    begin
            WIDGET_CONTROL, ev.top, GET_UVALUE = list
            if widget_info(list.TableID[0],/VALID_ID) eq 0 then return
            if widget_info(list.draw,/VALID_ID) eq 0 then return
            widget_control,list.draw,get_value=drawindex
            wset,drawindex
              error=CutPath(list.file,file=file)
              path=SelFile(list.path,file+'.ps','*.ps','Save Display To....')
              if path eq '' then return

            old_dev = !d.name
            old_fancy = !fancy
            set_plot, 'PS'
            !fancy = 0
            device, file=path, /landscape,/COLOR
            PlotProfiles,list
            device,/close
            !fancy = old_fancy
            set_plot, old_dev

            endcase
        3:    begin
            WIDGET_CONTROL, ev.top, GET_UVALUE = list
            if widget_info(list.TableID[0],/VALID_ID) eq 0 then return
            tmp=list.y
            list.y=list.x
            list.x=tmp
            widget_control,list.TableID[0],SET_TABLE_SELECT=[list.x,list.y,list.x,list.y]
            widget_control,list.TableID[0],set_table_view=[list.x,list.y]
            widget_control,list.TableID[1],SET_TABLE_SELECT=[list.x,list.y,list.x,list.y]
            widget_control,list.TableID[1],set_table_view=[list.x,list.y]
            WIDGET_CONTROL, ev.top, SET_UVALUE = list
            RefreshEditOverlap,ev.top
            endcase
        4:    begin
            WIDGET_CONTROL, ev.top, GET_UVALUE = list
            if widget_info(list.TableID[0],/VALID_ID) eq 0 then return

            str=PromptNumber(stringr(list.y)+','+stringr(list.x),ev.top,'Enter table coord [row,col]:')
            starr=str_sep(str,',')
            val=lonarr(2)
            reads,str,val
            n=list.nProf-1
            list.y=0>val[0]<n
            list.x=0>val[1]<n

            widget_control,list.TableID[0],SET_TABLE_SELECT=[list.x,list.y,list.x,list.y]
            widget_control,list.TableID[0],set_table_view=[list.x,list.y]
            widget_control,list.TableID[1],SET_TABLE_SELECT=[list.x,list.y,list.x,list.y]
            widget_control,list.TableID[1],set_table_view=[list.x,list.y]
            WIDGET_CONTROL, ev.top, SET_UVALUE = list
            RefreshEditOverlap,ev.top
            endcase
        5:    begin
            window,winid()
            PlotGroupProfile,ev
            endcase
        6:    begin
            WIDGET_CONTROL, ev.top, GET_UVALUE = list
            if widget_info(list.TableID[0],/VALID_ID) eq 0 then return
            if widget_info(list.draw,/VALID_ID) eq 0 then return
            widget_control,list.draw,get_value=drawindex
            wset,drawindex
              error=CutPath(list.file,file=file)
              path=SelFile(list.path,file+'.ps','*.ps','Save Group Profile To....')
              if path eq '' then return

            old_dev = !d.name
            old_fancy = !fancy
            set_plot, 'PS'
            !fancy = 0
            device, file=path, /landscape,/COLOR
            PlotGroupProfile,ev
            device,/close
            !fancy = old_fancy
            set_plot, old_dev
            endcase
      ELSE:

    ENDCASE
    endcase
9:    begin
    if ev.type eq 4 then begin
        if (ev.SEL_LEFT ne -1 and ev.SEL_TOP ne -1) then begin
            WIDGET_CONTROL, ev.top, GET_UVALUE = list
            case ev.id of
            list.TableID[0]:begin
                            ID=list.TableID[1]
                            endcase
            list.TableID[1]:begin
                            ID=list.TableID[0]
                            endcase
            endcase
            widget_control,ID,SET_TABLE_SELECT=[ev.SEL_LEFT, ev.SEL_TOP, ev.SEL_RIGHT, ev.SEL_BOTTOM]
            viewoff=widget_info(ev.id,/TABLE_VIEW)
            widget_control,ID,set_table_view=viewoff
            list.x=ev.SEL_LEFT
            list.y=ev.SEL_TOP
            WIDGET_CONTROL, ev.top, SET_UVALUE = list
            RefreshEditOverlap,ev.top
        endif
    endif
    endcase
endcase

END;PRO EditOverlap_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PRO EditOverlap

EditOverlapbase = WIDGET_BASE(TITLE = "EditOverlap", /COLUMN)    ;create the main base
menu = Cw_Bgroup(EditOverlapbase, ['Exit', 'Read', 'Save PS','Swap Col<>Row','GoTo','GroupProfile','GroupProfilePS'], /ROW, UVALUE="THEBUTTON")

ControlDevice,Platform,Visual,WD
result=ReadINI('$XRD mskpath:')
path=result.(0)
tmp=CutPath(result.(1),file=file)
file=file+'.txt'

basewrap=WIDGET_BASE(EditOverlapbase, /ROW)
dynID = WIDGET_BASE(basewrap, /COLUMN)
window,/Free,/pixmap
xs=!d.x_size-40
ys=!d.y_size-40
wdelete,!d.window
draw=widget_draw(basewrap,xsize=xs,ysize=ys)
label=widget_label(dynID,value='',uname='label1',xsize=xs/2)
label=widget_label(dynID,value='',uname='label2',xsize=xs/2)
label=widget_label(dynID,value='',uname='label3',xsize=xs/2)

bool=ReadOverlapFile(path,file,result=result)
nmag=3
if bool then begin
    ID = WIDGET_BASE(dynID, /COLUMN)
    ID1 = EditOverlapField(ID, result.XY,NAME='Overlap XY')
    ID2 = EditOverlapField(ID, result.XbY,NAME='Overlap XRowByYCol')
    ID3 = widget_draw(ID,xsize=result.nProf*nmag,ysize=result.nProf*nmag,$
        button=1,motion=1,Event_Pro='IMG_PROCESS_EVENTS')
    list={ID:ID,path:path,file:file,Profiles:result.Profiles,x:0L,y:0L,draw:draw,dynID:dynID,TableID:[ID1,ID2,ID3],$
    nGroups:result.nGroups,maskgroup:result.maskgroup,TableMean:result.TableMean,nmag:nmag,nProf:result.nProf,IndMask:result.IndMask}

endif else list={ID:0L,path:path,file:file,Profiles:PTR_NEW(),x:0L,y:0L,draw:draw,dynID:dynID,TableID:0L,nmag:nmag}

WIDGET_CONTROL, EditOverlapbase, SET_UVALUE=list
WIDGET_CONTROL, EditOverlapbase, /REALIZE
if bool then begin
    widget_control,list.TableID[2],get_value=drawindex
    Wset,drawindex
    img=rebin(result.XY,list.nProf*list.nmag,list.nProf*list.nmag,/sample)
    LoadctRev,0,/silent,rev=-1
    tvscl,img,/order

    RefreshEditOverlap,EditOverlapbase
endif
XManager, "EditOverlap", EditOverlapbase,/NO_BLOCK,cleanup='CleanUp_EditOverlap'

END ;PRO editoverlap
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

pro printwerrorstate,id,inc=inc

printw,id,!ERROR_STATE.MSG_PREFIX+'IDL error state:',inc=inc

Catch, /Cancel
help, /last_message, output=errors
for i=0,n_elements(errors)-1 do $
    printw,id,errors[i],inc=inc

;printw,id,!ERROR_STATE.MSG_PREFIX+!ERROR_STATE.name,inc=inc
;
;str=!ERROR_STATE.MSG_PREFIX+strsplit(!ERROR_STATE.msg,string([10b,13b]),/extract,count=n)
;for i=0,n-1 do $
;    printw,id,str[i],inc=inc
;
;str=!ERROR_STATE.MSG_PREFIX+strsplit(!ERROR_STATE.sys_msg,string([10b,13b]),/extract,count=n)
;for i=0,n-1 do $
;    printw,id,str[i],inc=inc
    
end;pro printw,id
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro printw,id,str,inc=inc
if n_elements(id) eq 0 then return

if not widget_info(id,/valid_id) then begin
    print,keyword_set(inc)?string(replicate(byte(' '),3*inc))+str:str
    return
endif
widget_control,id,get_uvalue=n,/NO_COPY
nextra=n_elements(str)
n[0]+=nextra

if n[0] gt n[1] then begin
    ;widget_control,id,get_value=buffer
    ;n2=n_elements(str)<(n_elements(buffer)-1)
    ;widget_control,id,set_value=buffer[n2:*]
    n[0]=nextra
    append=0
endif else append=1

widget_control,id,set_uvalue=n
widget_control,id,set_value=keyword_set(inc)?string(replicate(byte(' '),3*inc))+string(str):string(str),append=append,/NO_COPY
if !version.os_family eq 'Windows' then widget_control,id,SET_TEXT_TOP_LINE=(n[0]>0);-n[2]
end;pro printw
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro makeoutlogcleanup,ID
end;makeoutlogcleanup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro makeoutlogevent,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

TopBase=widget_info(ev.top,FIND_BY_UNAME='TopBase')
; This is not the top-level base!!!!!
; Top-level base: ev.top is resized, not TopBase

case widget_info(ev.id,/type) of
0:    begin
    IF TAG_NAMES(ev, /STRUCTURE_NAME) EQ 'WIDGET_KILL_REQUEST' THEN begin
        return
    endif else begin
        WIDGET_CONTROL, TopBase, GET_UVALUE=list, /NO_COPY
        TLBsize=list.TLBsize
        propagateresize,TopBase,TLBsize,'Log'
        list.TLBsize=TLBsize
        WIDGET_CONTROL, TopBase, SET_UVALUE=list, /NO_COPY
    endelse
    endcase
1:    begin
    widget_control,ev.id,get_value=val
    case val of
    'clear':begin
            ; Clear log
            widget_control,ev.id,get_uvalue=id
            widget_control,id,get_uvalue=n
            n[0]=1
            widget_control,id,set_uvalue=n
            widget_control,id,set_value=''
            widget_control,id,SET_TEXT_TOP_LINE=n[0]
            endcase
    'save':    begin
            ; Save log
            widget_control,ev.id,get_uvalue=id
            widget_control,id,get_uvalue=n
            if n[0] ne 0 then begin
                file=''
                file = DIALOG_PICKFILE(/WRITE, FILTER = ['*.txt','*.*'],/OVERWRITE_PROMPT)
                if file ne '' then begin
                    widget_control,id,get_value=buffer
                    if openw_safe(lun,file) then return
                    for i=0l,n[0]-1 do printf,lun,buffer[i]
                    free_lun,lun
                endif
            endif
            endcase
    endcase
    endcase
endcase

end;makeoutlogevent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function makeoutlog,ParentID=ParentID,standalone=standalone

; Maximum base tree (* is optional):
; ----------------------------------
; TLB*(ev.top, produces the resize event)
;  |
;  |_\  ParentID*(empty child base)
;    /        |
;              |_\  TopBase(made in this application, contains all widgets of makeoutlog)
;                /

; If ParentID is set, then group_leader will be the same.

; ----Screen Dimensions----
TLBRestSize=lonarr(2)
TopBaseRestSize=[15,65]
if not keyword_set(ParentID) then ParentID=-1L

if widget_info(ParentID,/VALID_ID) then begin

    TLB_ID=widget_info(ParentID,/PARENT)
    temp=TLB_ID
    while (temp ne 0) do begin
        TLB_ID=temp
        temp=widget_info(TLB_ID,/PARENT)
    endwhile
    if (TLB_ID eq 0) then return,-1L ; Parent base must be a child base.
    group_leader=TLB_ID
    
    temp=widgetsize(ParentID)-[15,65]
    xdim=temp[0]
    ydim=temp[1]
    base=widget_base(ParentID,/column,uname='TopBase')
endif else begin
    xdim = 500
    ydim = 500
    base=widget_base(/column,title='Output Log',/TLB_SIZE_EVENTS,$
                    uname='TopBase',TLB_KILL_REQUEST_EVENTS=~keyword_set(standalone),$
                    TLB_FRAME_ATTR=(~keyword_set(standalone))*8)
endelse

StatusText=widget_text(base,uname='Log',uvalue=[0,500,10],/scroll,/wrap,scr_xsize=xdim,scr_ysize=ydim)
button=widget_button(base,value='clear',uvalue=StatusText)
button=widget_button(base,value='save',uvalue=StatusText)

TLBsize=widgettreesizeinit(base,'Log')

WIDGET_CONTROL, base, /REALIZE

randomid=randlong()
list = {ParentID:ParentID,$    ; Parent ID (-1 if no parent)
        TLBsize:TLBsize,$
        randomid:randomid,$
        delete:0b}

WIDGET_CONTROL, base, SET_UVALUE=list, /NO_COPY
Xmanager,'makeoutlog',base,/NO_BLOCK,cleanup='makeoutlogcleanup',event_handler='makeoutlogevent',GROUP_LEADER=group_leader

return,[StatusText,base,randomid]
end;function makeoutlog
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
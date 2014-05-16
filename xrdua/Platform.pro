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

function Version
return,'6.4.3.1'
; In principle, in subsequent releases, the major number is increased when there are significant 
; jumps in functionality, the minor number is incremented when only minor features or 
; significant fixes have been added, and the revision number is incremented when minor bugs are fixed.

; DON'T CHANGE THIS FUNCTION UNLESS YOU KNOW WHAT YOU'RE DOING
end; function Version
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CompareVersion,version1,version2
; Return:
; equal = 0
; greater than = 1
; less than = -1
if version1 eq version2 then return,0
tmp1=long(strsplit(version1,'.',/extract))
tmp2=long(strsplit(version2,'.',/extract))
for i=0,3 do begin
    if tmp1[i] gt tmp2[i] then return,1
    if tmp1[i] lt tmp2[i] then return,-1
endfor
return,0
end;function CompareVersion
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function winidvalid,id
DEVICE, WINDOW_STATE=wstate
return,wstate[id>0]
end;function winidvalid
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function winid
DEVICE, WINDOW_STATE=wstate
ind=where(wstate eq 0,ct)
if ct eq 0 then begin
    window,/free
    return,!d.window
endif else begin
    return,ind[0]
endelse
end;function winid
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LoadctX,list
device,DECOMPOSE=0
case list.ctl.type of
0:    LoadctRev,list.ctl.nr,/silent,rev=-1
1:    TVLCT, (*list.ctl.R), (*list.ctl.G), (*list.ctl.B)
endcase
end;LoadctX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ntableadd,nadd
return,n_tags(!color)+nadd-1
end;function ntableadd
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LoadctRev, table_number,silent=silent,rev=rev,gamma=gamma,low=low,high=high
if not keyword_set(rev) then rev=-1

if table_number eq 0 then sign=rev else begin
    sign=table_number/abs(table_number)
    table_number=abs(table_number)
endelse

tablemax=n_tags(!color)-1
if table_number gt tablemax then begin
    case table_number-tablemax of
    1:    begin
        ; white(background),...
        r=[255,0,0,225,0,255,128,123,35]
        g=[255,195,0,0,225,0,128,25,123]
        b=[255,0,255,225,225,0,128,13,200]
        endcase
    2:    begin
        r=congrid([0,255,255,255],256,/interp)
        g=congrid([0,0,255,255],256,/interp)
        b=congrid([0,0,0,255],256,/interp)
        endcase
    else: return
    endcase

endif else begin
    filename = filepath('colors1.tbl', subdir=['resource', 'colors'])
    if openr_safe(lun, filename, /block) then return

;    if not keyword_set(silent) then begin
;        ntables = 0b
;        readu, lun, ntables
;        names = bytarr(32, ntables)
;        point_lun, lun, ntables * 768L + 1  ;Read table names
;        readu, lun, names
;        names = strtrim(names, 2)
;
;        if sign eq -1 then str='-' else str='+'
;        message,'Loading table ' +str+ names[table_number],/INFO
;    endif

    aa=assoc(lun, bytarr(256),1)    ;Read 256 long ints (8-bit)
    r = aa[table_number*3]
    g = aa[table_number*3+1]
    b = aa[table_number*3+2]

    free_lun,lun
endelse

if sign eq -1 then begin
    r=reverse(r,/overwrite)
    g=reverse(g,/overwrite)
    b=reverse(b,/overwrite)
endif


if n_elements(low) eq 0 then low=0
if n_elements(high) eq 0 then high=1
if n_elements(gamma) eq 0 then gamma=1
if low ne 0 or high ne 1 or gamma ne 1 then begin
    x = findgen(256)/255
    x = (0>(x-low))/high<1
    x ^= gamma
    x = 0>fix(x*255)<255
    
    r=r[x]
    g=g[x]
    b=b[x]
endif

tvlct,r,g,b
end;pro LoadctRev
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function enlarge_dyn,xsize,ysize,pixmap,ID,draw,drawindex,xss,yss

ON_ERROR, 3 ; skip line when error

xs0=xsize
ys0=ysize

catch, errstat
if (errstat ne 0) then begin
    xsize=xss
    ysize=yss
endif

DEVICE, WINDOW_STATE=wstate
if wstate[pixmap] then wdelete, pixmap
Window, /Free, xsize=xsize, ysize=ysize, /Pixmap
pixmap=!d.window
widget_control, draw, draw_xsize=xsize, draw_ysize=ysize
; => when error, user is supplied with messagebox. Press 'OK' and the program proceeds

; Check whether enlarging failed:
if not widget_info(draw,/VALID_ID) then begin
    ; window is deleted after user pressed 'OK'
    return,create_dyn(xsize,ysize,pixmap,ID,draw,drawindex,xss,yss)
endif else widget_control,draw,get_value=drawindex

return,xsize*ysize

end;function enlarge_dyn
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function create_dyn,xsize,ysize,pixmap,ID,draw,drawindex,xss,yss

ON_ERROR, 3 ; skip line when error

catch, errstat
if (errstat ne 0) then begin
    xsize=xss
    ysize=yss
endif

; Destroy
if widget_info(draw,/valid) then WIDGET_CONTROL, draw, /DESTROY
DEVICE, WINDOW_STATE=wstate
if wstate[pixmap] then wdelete, pixmap

; Make new
Window, /Free, xsize=xsize, ysize=ysize, /Pixmap
pixmap=!d.window
draw = WIDGET_DRAW(ID,xsize=xsize, ysize=ysize, /scroll,$
       x_scroll_size=xss, y_scroll_size=yss,uname='drawdyn',uvalue=0B,$
       /VIEWPORT_EVENTS)
; => when error, user is supplied with messagebox. Press 'OK' and the program proceeds

; Handle made or not made
if not widget_info(draw,/VALID_ID) then begin
    if (xsize ge xss) and (ysize ge yss) then a[0]=1; create error
    draw=0L ; dynamic window smaller than static window => no window
    drawindex=0L
    wdelete,pixmap
    pixmap=0L
endif else widget_control,draw,get_value=drawindex

return,xsize*ysize

end;function create_dyn
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CreateDyn,draw,drawindex,pixmap,xsize,ysize,xss,yss,xs,ys,sf,ID,IDtext,rebuild=rebuild

; Prepare: delete pixmap and check for draw widget
bEnlarge=widget_info(draw,/VALID_ID) and ~keyword_set(rebuild)

if bEnlarge then Msize=enlarge_dyn(xsize,ysize,pixmap,ID,draw,drawindex,xss,yss) $
else Msize=create_dyn(xsize,ysize,pixmap,ID,draw,drawindex,xss,yss)

; If to large, then adapt scale factor
sfnew=sqrt(float(xs)*ys/Msize)
if sfnew gt sf then sf=sfnew
widget_control,IDtext,set_value=stringr(1./sf)

return,1B
end;function CreateDyn
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function test_pixmap_size

s0=2000L
s=s0

catch, errstat
if (errstat ne 0) then s=s-1

window, /pixmap, /free, xsize=s, ysize=s
wdelete, !d.window

return,s*s

end;function test_pixmap_size
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SetXRDUAEnv,Platform,WD

WD=XRDUAworkingDIR()
return

if Platform eq 'WIN' then begin
    if WD eq '' then begin
        WD=''
        while WD eq '' do begin
            path=''
            path=DIALOG_PICKFILE(path=XRDUAworkingDIR(),/DIRECTORY,title='Select XRDUA installation directory...')
            B=CutPath(path,path=WD)
        endwhile
    endif

    AddRegistry,'xrdua.reg','HKEY_CURRENT_USER\Environment','xrduadir',WD
    r=dialog_message('Added %xrduaini% to the registry. To make it permanent, edit the "Environment Variables" is the "Control Pannel" and press "OK" without changing anything or reboot your PC.',/information)
endif else begin
    if WD eq '' then begin
        WD=''
        while WD eq '' do begin
            path=''
            path=DIALOG_PICKFILE(path=XRDUAworkingDIR(),/DIRECTORY,title='Select XRDUA installation directory...')
            B=CutPath(path,path=WD)
        endwhile
    endif

    shell=GETENV('SHELL')
    shell=strmid(shell,rstrpos(shell,Path_Sep())+1)
    case strlowcase(shell) of
    'csh':    begin ;C shell
            file='~/.cshrc'
            command='setenv xrduadir '
            endcase
    'tcsh':    begin ;TC shell
            file='~/.cshrc' ; '.tcshrc'
            command='export xrduadir='
            endcase
    'sh':    begin ;Bourne shell
            file='~/.profile'
            command='export xrduadir='
            endcase
    'bash':    begin ;Bash
            file='~/.bashrc'
            command='export xrduadir='
            endcase
    'ksh':    begin ;Korn shell
            file='~/.kshrc'
            command='export xrduadir='
            endcase
    else:     begin
            r=dialog_message('Unkown shell. Add %xrduaini% manually to environment.',/information)
            return
            endcase
    endcase
    file=file_search(file,count=bfound)
    if bfound eq 0 then begin
        r=dialog_message('Unkown shell. Add %xrduaini% manually to environment.',/information)
        return
    endif
    file=file[0]

    if openu_safe(lun,file) then return
    line=''
    len=strlen(command)
    pos0=-1
    while (strmid(line,0,len) ne command) and (not eof(lun)) do begin
        point_lun,-lun,pos0
        readf,lun,line
    endwhile
    if (strmid(line,0,len) eq command) then begin
        header=''
        while (not eof(lun)) do begin
            readf,lun,line
            header=[header,line]
        endwhile
        point_lun,lun,pos0
        TRUNCATE_LUN, lun
        printf,lun,command+WD
        for i=1,n_elements(header)-1 do printf,lun,header[i]
    endif else begin
        free_lun,lun
        if openw_safe(lun,file,/append) then return
        printf,lun,command+WD
    endelse
    free_lun,lun

    r=dialog_message('Added %xrduaini% to shell script '+file+'. Exit this shell and start a new one to make it permanent.',/information)
endelse
end;pro SetXRDUAEnv
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ControlDevice,Platform,Visual,SHORT=SHORT

ON_ERROR, 3 ; skip line when error

; Save device name
PossiblePlatforms=['MacOS', 'unix', 'vms', 'Windows']
DevicePlatforms=['MAC','X','X','WIN']
Platform=!version.os_family
ind=where(strlowcase(PossiblePlatforms) eq strlowcase(Platform))
if (ind[0] eq 0) or (ind[0] eq 2) then $
    r=dialog_message('This program was not tested on this platform!',/information)
Platform=DevicePlatforms[ind[0]]
if keyword_set(SHORT) then return

; Forse IDL to:
; - use Indexed Color Model (i.e. handle TrueColor displays)
;         Decomposed 24bit colors (decompose=1): color='00FFFF'x
;         Undecomposed 8bit colors (decompose=0 or 256bit screen): color=123b
; - perform backing store by IDL
DEVICE, DECOMPOSE = 0, retain=2
!EXCEPT=0

; When X platform: use TrueColor to prevent flashing
;if Platform eq 'X' then device, TRUE_COLOR=16
Device, Get_Visual_Name=Visual

; Get maximum pixels in a pixmap
;Msize=test_pixmap_size()

; Set working directory
WD=GETENV('xrduadir')
if WD eq '' then SetXRDUAEnv,Platform,WD
char=Path_Sep()
if strmid(WD,strlen(WD)-1,1) ne char then WD+=char
print,'Working directory: '+WD
cd,WD

end;pro ControlDevice
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

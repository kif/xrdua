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

pro fplot_obj_proplotdefault, in

key=in.key
if key eq 0 then begin ; Change Z-range
    list=in.list
    oGroup=list.oGroup
    zr=in.zr

    oChildArr = oGroup->IDL_Container::Get(ISA='IDLgrSurface',/all, COUNT=nKids)
    oChildArr[0]->Setproperty, ZCOORD_CONV=NormCoord(zr)

endif else begin ; Add Objects
    list=in.list
    oGroup=list.oGroup
    ; ----Create the surface----
    oSurface = OBJ_NEW('IDLgrSurface', (list.list2).z , STYLE=2, SHADING=0, $
               COLOR=[60,60,255], BOTTOM=[64,192,128], $
               XCOORD_CONV=NormCoord((list.list2).xr),$
               YCOORD_CONV=NormCoord((list.list2).yr),$
               ZCOORD_CONV=NormCoord((list.list2).zr))
    oGroup->Add, oSurface
endelse

end;pro fplot_obj_proplotdefault
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NormCoord,range,position

IF (N_Elements(position) EQ 0) THEN position = [0.0, 1.0] ELSE $
   position=Float(position)
range = Float(range)

d=range[1]-range[0]
if d eq 0 then return,[0.,1]
scale = [((position[0]*range[1])-(position[1]*range[0])) / d, (position[1]-position[0])/d]

RETURN, scale

end;function NormCoord
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NonOrt2Ort,a,b,c

; See UnitCelltoRep

a*=!pi/180.
b*=!pi/180.
c*=!pi/180.
cosa=cos(a)
sina=sin(a)
cosb=cos(b)
sinb=sin(b)
cosc=cos(c)
sinc=sin(c)
V=sqrt(1.-cosa*cosa-cosb*cosb-cosc*cosc+2*cosa*cosb*cosc)

; Transpose M to use it like this:
; xyz##=M where xyz has three columns
M=[[1,0,0],$
   [cosc,sinc,0],$
   [cosb,(cosa-cosb*cosc)/sinc,V/sinc]]

return,M
end;function NonOrt2Ort
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rotx,a
; Change of basis matrix
return,[[1,0,0,0],$
    [0,cos(a),sin(a),0],$
    [0,-sin(a),cos(a),0],$
    [0,0,0,1]]
end;function Rotx
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Roty,b
; Change of basis matrix
return,[[cos(b),0,-sin(b),0],$
    [0,1,0,0],$
    [sin(b),0,cos(b),0],$
    [0,0,0,1]]
end;function Roty
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rotz,c
; Change of basis matrix
return,[[cos(c),sin(c),0,0],$
    [-sin(c),cos(c),0,0],$
    [0,0,1,0],$
    [0,0,0,1]]
end;function Rotz
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DefTrans

;oGroup->Rotate, [1,0,0], -90
;oGroup->Rotate, [0,1,0], 30
;oGroup->Rotate, [1,0,0], 30

M=Rotx(75.*!dpi/180.)
M##=Rotz(30.*!dpi/180.)

M=Rotx(70.*!dpi/180.)
M##=Rotz(20.*!dpi/180.)


;M=Rotz(-atan(1,2))
;M=Rotx(0)

return,M

end;function DefTrans
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_flplot_obj,id

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list

if n_elements(list) eq 0 then return
; Destroy the objects.
if OBJ_VALID(list.oHolder) then OBJ_DESTROY, list.oHolder
if OBJ_VALID(list.oWindow) then OBJ_DESTROY, list.oWindow
end;pro CleanUp_flplot_obj
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro flplot_obj_refresh,list,zr,rebuild=rebuild

key=keyword_set(rebuild)

; Delete all user objects
if key then begin
    list.oZAxis->GetProperty,range=zr
    oChildArr = list.oGroup->IDL_Container::Get(/all, COUNT=nKids)
    for i=list.noGroupFix,nkids-1 do begin
        list.oGroup->IDL_Container::Remove,oChildArr[i]
        OBJ_DESTROY,oChildArr[i]
    endfor
endif

; Axis ranges
xr=list.list2.xr
yr=list.list2.yr
if n_elements(zr) eq 0 then zr=list.list2.zr

; Recall plot
tmp={list:list,zr:zr,key:key}
call_procedure,(list.list2).proplot,tmp
if key then begin
    tmp.key=0
    call_procedure,(list.list2).proplot,tmp
endif
list=tmp.list

; Reset axis and draw
list.oXAxis->SetProperty,range=xr,XCoord_Conv=NormCoord(xr,list.PosX)
list.oYAxis->SetProperty,range=yr,YCoord_Conv=NormCoord(yr,list.PosY)
list.oZAxis->SetProperty,range=zr,ZCoord_Conv=NormCoord(zr,list.PosZ)
list.oWindow->Draw, list.oView
end;pro flplot_obj_refresh
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro flplot_obj_event8,ev,TopBase
WIDGET_CONTROL, TopBase, GET_UVALUE=list, /NO_COPY
if (ev.index eq 0) then begin
    list.oTrack->Reset,list.tbCenter,list.tbRadius,CONSTRAIN=0,AXIS=0
    list.bDataConstr=0b
endif else if ((ev.index le 3) and (ev.index ge 1)) then begin
    list.oTrack->Reset,list.tbCenter,list.tbRadius,CONSTRAIN=1,AXIS=ev.index-1
    list.bDataConstr=0b
endif else begin
    list.oTrack->Reset,list.tbCenter,list.tbRadius,CONSTRAIN=1,AXIS=ev.index-4
    list.bDataConstr=1b
endelse
WIDGET_CONTROL, TopBase, SET_UVALUE=list, /NO_COPY
end;pro flplot_obj_event8
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro flplot_obj_event4,ev,TopBase
    WIDGET_CONTROL, TopBase, GET_UVALUE=list, /NO_COPY

    ; ----Expose event----
    IF (ev.type EQ 4) THEN BEGIN
       list.oWindow->Draw, list.oView
       WIDGET_CONTROL, TopBase, SET_UVALUE=list, /NO_COPY
       RETURN
    ENDIF

     possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
    bType = possibleEventTypes[ev.type]
    possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
    bPress = possibleButtons[ev.press]

    ; ----Init change----
    list.oGroup->GetProperty, TRANSFORM=t
    IF (bType EQ 'DOWN') THEN begin
        list.bState=bPress
        WIDGET_CONTROL, list.draw, /DRAW_MOTION
        IF (list.bState eq 'RIGHT') THEN BEGIN
            list.xs=ev.x
            list.ys=ev.y
        endif
        if list.bDataConstr then list.starttransform=t
    endif

    ; ----Handle rotation----
    IF list.bDataConstr THEN starttransform=list.starttransform
    ; Update will calculate the transformation matrix to transform
    ; 'DOWN'-coordinate to current coordinate on the trackball
    ;
    ; When constrainted to an axis rotation, the two coordinates are
    ; first projected to the plain perpendicular to it, e.g.
    ;     about ScreenX: project on plane perp to [1,0,0] => don't set START_TRANSFORM
    ;     about DataX: project on plane perp to [1,0,0]#t => set START_TRANSFORM
    bHaveTransform = list.oTrack->Update( ev, START_TRANSFORM=starttransform, TRANSFORM=qmat )
    IF (bHaveTransform NE 0) THEN BEGIN
        list.oGroup->SetProperty, TRANSFORM=t#qmat
        list.oWindow->Draw, list.oView
    ENDIF

    ; ----Handle translation----
    if (bType EQ 'MOTION') THEN BEGIN
        IF (list.bState eq 'RIGHT') THEN BEGIN; translation
            list.oWindow->GetProperty, SCREEN_DIMENSIONS=s

            dx = 6*(list.xs - ev.x) / Float(s[0])
             dy = 6*(list.ys - ev.y) / Float(s[1])
             list.xs=ev.x
            list.ys=ev.y

            list.oView->GetProperty, VIEWPLANE_RECT=rect
            list.oView->SetProperty, VIEWPLANE_RECT=rect+[dx,dy,0,0]

            list.oWindow->Draw, list.oView
        ENDIF
    ENDIF

    ; ----Finish----
    if (bType eq 'UP') then WIDGET_CONTROL, list.draw, DRAW_MOTION=0
    WIDGET_CONTROL, TopBase, SET_UVALUE=list, /NO_COPY
end;pro flplot_obj_event4
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro flplot_obj_event1,ev,TopBase
     widget_control,ev.id,get_value=val
     WIDGET_CONTROL, TopBase, GET_UVALUE=list, /NO_COPY
     case val of
         'OK':     WIDGET_CONTROL, TopBase, /DESTROY
        '+Z':    begin
                list.oZAxis->GetProperty, Range=zr
                zr[1]-=(zr[1]-zr[0])*0.2
                ;list.oZAxis->SetProperty,range=zr,ZCoord_Conv=NormCoord(zr,list.PosZ)
                flplot_obj_refresh,list,zr
                endcase
        '-Z':    begin
                list.oZAxis->GetProperty, Range=zr
                zr[1]+=(zr[1]-zr[0])*0.2
                ;list.oZAxis->SetProperty,range=zr,ZCoord_Conv=NormCoord(zr,list.PosZ)
                flplot_obj_refresh,list,zr
                endcase
        'Reset':begin
                list.oView->SetProperty, VIEWPLANE_RECT=list.myview
                list.oGroup -> SetProperty, TRANSFORM =DefTrans()
                zr=(list.list2).zr
                ;list.oZAxis->SetProperty,range=zr,ZCoord_Conv=NormCoord(zr,list.PosZ)
                flplot_obj_refresh,list,zr
                endcase
        'XY':begin
                M=identity(4)
                list.oGroup -> SetProperty, TRANSFORM =M
                zr=(list.list2).zr
                ;list.oZAxis->SetProperty,range=zr,ZCoord_Conv=NormCoord(zr,list.PosZ)
                flplot_obj_refresh,list,zr
                endcase
        'YZ':begin
                M=Rotx(90.*!dpi/180.)
                M##=Rotz(90.*!dpi/180.)
                list.oGroup -> SetProperty, TRANSFORM =M
                zr=(list.list2).zr
                ;list.oZAxis->SetProperty,range=zr,ZCoord_Conv=NormCoord(zr,list.PosZ)
                flplot_obj_refresh,list,zr
                endcase
        'XZ':begin
                M=Rotx(90.*!dpi/180.)
                list.oGroup -> SetProperty, TRANSFORM =M
                zr=(list.list2).zr
                ;list.oZAxis->SetProperty,range=zr,ZCoord_Conv=NormCoord(zr,list.PosZ)
                flplot_obj_refresh,list,zr
                endcase
        'Axis on/off':begin
                list.bAxisHide=~list.bAxisHide
                oChildArr = list.oGroup->IDL_Container::Get(ISA='IDLgrAxis',/all, COUNT=nKids)
                for i=0l,n_elements(oChildArr)-1 do oChildArr[i]->SetProperty,hide=list.bAxisHide
                flplot_obj_refresh,list,zr
                endcase
        'Rebuild':begin
                flplot_obj_refresh,list,/rebuild
                endcase
        'Save Image...':begin
                path=''
                path=DIALOG_PICKFILE(path=list.path,file='cap.bmp',$
                     filter=['*.png','*.tif;*.tiff','*.jpg;*.jpeg','*.bmp','*.eps'],title='Save image ....',/OVERWRITE_PROMPT,/WRITE)
                if path ne '' then begin
                    p=CutPath(path,ext=ext)
                    ext=strlowcase(ext)
                    case ext of
                    '.eps':    begin
                            list.oWindow->GetProperty, units=units
                            list.oWindow->SetProperty, units=1
                            list.oWindow->GetProperty, dimensions=xy
                            list.oWindow->SetProperty, units=units

                            clipboard = Obj_New("IDLgrClipboard", Dimensions=xy, Units=1, $
                                        Resolution=[2.54/300., 2.54/300.]) ;
                            clipboard->Draw, list.oView, Filename=path, /PostScript, /Vector
                            obj_destroy  ,clipboard
                            endcase
                    else:    begin
                            flplot_obj_refresh,list,/rebuild
                            list.oWindow->GetProperty,image_data=image_data
                            case ext of
                            '.jpg': WRITE_JPEG, path, image_data,QUALITY=100,true=1
                            '.jpeg': WRITE_JPEG, path, image_data,QUALITY=100,true=1
                            '.tif': WRITE_TIFF,path,image_data,ORIENTATION=0
                            '.tiff': WRITE_TIFF,path,image_data,ORIENTATION=0
                            '.png': WRITE_PNG,path,image_data
                            '.bmp': WRITE_BMP, path, image_data,/RGB
                            else:
                            endcase
                            
                            endelse
                    endcase

                endif
                endcase
        else:     begin; Zoom
        
                ID=widget_info(ev.top,find_by_uname='zoom')
                case widget_info(ID,/DROPLIST_SELECT) of
                0: add=0.05
                1: add=0.1
                2: add=1
                endcase
                if val eq '+' then add=-add
                add=[-add,-add,2*add,2*add]

                list.oView->GetProperty, VIEWPLANE_RECT=rect
                list.oView->SetProperty, VIEWPLANE_RECT=rect+add
                list.oWindow->Draw, list.oView
                endelse
     endcase
     WIDGET_CONTROL, TopBase, SET_UVALUE=list, /NO_COPY
end;pro flplot_obj_event1
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro flplot_obj_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

TopBase=find_by_uname(ev.top,'fplotTopBase')
if ~widget_info(TopBase[0],/VALID_ID) then return
nTopBase=n_elements(TopBase)
; This is not the top-level base!!!!!
; Top-level base: ev.top is resized, not TopBase

case widget_info(ev.id,/type) of
 0 :begin
     for i=0l,nTopBase-1 do begin
         WIDGET_CONTROL, TopBase[i], GET_UVALUE=list, /NO_COPY
         TLBsize=list.TLBsize
         propagateresize,TopBase[i],TLBsize,'fplotdrawwidget',list.yxratio
        list.TLBsize=TLBsize
         WIDGET_CONTROL, TopBase[i], SET_UVALUE=list, /NO_COPY
    endfor
     endcase
 4 :begin
     ;ind=wheresibling(TopBase,ev.id,count=ct)
     ;if ct eq 0 then return
     ;TopBase=TopBase[ind[0]]
    for i=0l,nTopBase-1 do flplot_obj_event4,ev,TopBase[i]
     endcase
 1 :begin
     widget_control,ev.id,get_value=val
     if val eq 'Save Image...' then begin
         ind=wheresibling(TopBase,ev.id,count=ct)
         if ct eq 0 then return
         TopBase=TopBase[ind[0]]
         flplot_obj_event1,ev,TopBase
     endif else $
        for i=0l,nTopBase-1 do flplot_obj_event1,ev,TopBase[i]
     endcase
 8: begin
     widget_control,ev.id,get_uvalue=uval
     case uval of
     'constraint':    begin
                 ;ind=wheresibling(TopBase,ev.id,count=ct)
                 ;if ct eq 0 then return
                 ;TopBase=TopBase[ind[0]]
                for i=0l,nTopBase-1 do flplot_obj_event8,ev,TopBase[i]
                 endcase
    'zoom':        begin
                endcase
    endcase
    endcase
endcase
end; pro flplot_obj_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro flplot_obj,list2,topid=topid,GROUP_LEADER=group_leader,xtitle=xtitle,$
    ytitle=ytitle,ztitle=ztitle,xtrue=xtrue,ytrue=ytrue,ztrue=ztrue,$
    ayz=ayz,bxz=bxz,cxy=cxy,parentbase=parentbase,title=title,xdim=xdim,ydim=ydim,handle=base,zoom=zoom,zclip=zclip

; Required input:
; ---------------
; list2 should have following fields:
;    'zr','xr','yr','proplot'

; Optional input:
; ---------------
; topid: output ID of top base
; group_leader: top ID that called this routine
; xtitle,ytitle,ztitle,title
; xtrue,ytrue,ztrue: set these => non-isotropic axis
;        (e.g. xr=[0,5] and yr=[0,1] then X-axis 5 times longer then Y-axis)
; ayz,bxz,cxy: angles between axis (if not set: orthogonal)
; parentbase: EMPTY child of top level base
; xdim,ydim: draw_widget size

; Object Hierarchy:
; ----------------
; oWindow
; |
; |_____oTrack
;    |__oView__oTop
;                 |__oLight1
;                 |__oLight2
;                 |__oGroup
;                       |__oXAxis__oXtitle
;                         |                  |__otitlefont
;                         |                  |__oticfont
;                      |
;                         |__oYAxis__oYtitle
;                         |                  |__otitlefont
;                         |                  |__oticfont
;                      |
;                         |__oZAxis__oZtitle
;                         |                  |__otitlefont
;                         |                  |__oticfont
;                      |
;                         |__...(user defined proplot)
;
; oHolder keeps all objects that must be destroyed explicitly:
;    oTrack, oView, OX/Y/Ztitle, otitlefont, oticfont


; ----Extract 3 ranges from input----
names = TAG_NAMES(list2)

ind=where(names eq 'XR')
if ind[0] eq -1 then list2=create_struct(list2,'xr',[0L,1])
if list2.xr[0] eq list2.xr[1] then list2.xr[1]=list2.xr[0]+1
xr=list2.xr
ind=where(names eq 'YR')
if ind[0] eq -1 then list2=create_struct(list2,'yr',[0L,1])
if list2.yr[0] eq list2.yr[1] then list2.yr[1]=list2.yr[0]+1
yr=list2.yr
ind=where(names eq 'ZR')
if ind[0] eq -1 then list2=create_struct(list2,'zr',[0L,1])
if list2.zr[0] eq list2.zr[1] then list2.zr[1]=list2.zr[0]+1
zr=list2.zr
ind=where(names eq 'PROPLOT')
if ind[0] eq -1 then begin
    list2=create_struct(list2,'proplot','fplot_obj_proplotdefault')
endif

xs=xr[1]-xr[0]
ys=yr[1]-yr[0]
zs=zr[1]-zr[0]

; ----Create Applet----
if not keyword_set(title) then title='3D rep.'
if n_elements(parentbase) eq 0 then parentbase=0L
if widget_info(parentbase,/VALID_ID) then $
    base=widget_base(parentbase,/column,title=title,uname='fplotTopBase') $
else $
    base=widget_base(/column,title=title,/TLB_SIZE_EVENTS,uname='fplotTopBase')

if not keyword_set(xdim) then xdim=500
if not keyword_set(ydim) then ydim=500
draw=widget_draw(base,xsize=xdim,ysize=ydim,GRAPHICS_LEVEL=2,$
    /BUTTON_EVENTS, /EXPOSE_EVENTS,uname='fplotdrawwidget')
base2=widget_base(base,/row)
button1=widget_button(base2,value='+')
button2=widget_button(base2,value='-')
button1=widget_button(base2,value='+Z')
button2=widget_button(base2,value='-Z')
button3=widget_button(base2,value='XY')
button3=widget_button(base2,value='YZ')
button3=widget_button(base2,value='XZ')
button3=widget_button(base2,value='Reset')
button3=widget_button(base2,value='Rebuild',uname='fplotrebuildwidget')
button3=widget_button(base2,value='Axis on/off')
base2=widget_base(base,/row)
button4=widget_droplist(base2,value=['Unconstrained','RotX Screen','RotY Screen',$
                                'RotZ Screen','RotX Data','RotY Data','RotZ Data'],uvalue='constraint')
button6=widget_droplist(base2,value=['Slow zoom','Medium Zoom','Fast Zoom'],uvalue='zoom',uname='zoom')
button5=widget_button(base2,value='Save Image...')
;button=widget_button(base2,value='OK')
WIDGET_CONTROL, base, /REALIZE

topid=base
TLBsize=widgettreesizeinit(base,'fplotdrawwidget')

;widget_control,button,/INPUT_FOCUS
widget_control,draw,get_value=oWindow

; ----Create view----
off=0.5
if keyword_set(zoom) then off+=zoom
myview=[-off,-off,1+2*off,1+2*off]

if keyword_set(ZCLIP) then ZCLIP=abs(ZCLIP)*[1.,-1.] else ZCLIP=[3.,-3.]
oView = OBJ_NEW('IDLgrView',PROJECTION=1,$
                    ;COLOR=[255,0,0],$ ; background color
                    EYE=ZCLIP[0]+1,$ ; coord as for zclip (doesn't change anything for PROJ=1)
                            ; EYE must be greater than ZCLIP[0]
                    ZCLIP=ZCLIP,$; normal. coord. relative to
                                      ; surface plot range (x, y or z)
                    VIEWPLANE_RECT=myview ,$; normal. coord. relative to
                                                ; surface plot x and y range
                    LOCATION=[0,0],dimensions=[1,1],$; Occupy whole window if units=3
                    units=3); 3: normalised relative to the graphics destination's rect

; ----Create model----
oTop = OBJ_NEW('IDLgrModel')
oGroup = OBJ_NEW('IDLgrModel')
oTop->Add, oGroup

; ----Create axis----
otitlefont= OBJ_NEW('IDLgrFont','times',size=11-ZCLIP[0])
oticfont = OBJ_NEW('IDLgrFont','times',size=10-ZCLIP[0])
if not keyword_set(xtitle) then xtitle='X'
if not keyword_set(ytitle) then ytitle='Y'
if not keyword_set(ztitle) then ztitle='Z'
oXtitle=OBJ_NEW('IDLgrText',xtitle,FONT=otitlefont, Recompute_Dimensions=2)
oYtitle=OBJ_NEW('IDLgrText',ytitle,FONT=otitlefont, Recompute_Dimensions=2)
oZtitle=OBJ_NEW('IDLgrText',ztitle,FONT=otitlefont, Recompute_Dimensions=2)

RelLength=[xs,ys,zs]
TrueBool=[keyword_set(xtrue),keyword_set(ytrue),keyword_set(ztrue)]
ind=where(TrueBool eq 1,ct)
PosX=[0,1.]
PosY=[0,1.]
PosZ=[0,1.]
if ct ge 2 then begin
    Pos=[[PosX],[PosY],[PosZ]]
    RelLength=RelLength[ind]
    RelLength=RelLength/float(max(RelLength))
    Pos[1,ind]=RelLength
    PosX=Pos[*,0]
    PosY=Pos[*,1]
    PosZ=Pos[*,2]
endif
;ticklen=0.001*[yr[1] - yr[0],xr[1] - xr[0],xr[1] - xr[0]]
ticklen=bytarr(3)
oXAxis = OBJ_NEW('IDLgrAxis', 0, RANGE=xr,TICKLEN=ticklen[0],$
title=oXtitle,/exact,TEXTPOS=0, XCOORD_CONV=NormCoord(xr,PosX))
oYAxis = OBJ_NEW('IDLgrAxis', 1, RANGE=yr,TICKLEN=ticklen[1],$
title=oYtitle,/exact,TEXTPOS=0, YCOORD_CONV=NormCoord(yr,PosY))
oZAxis = OBJ_NEW('IDLgrAxis', 2, RANGE=zr,TICKLEN=ticklen[2],$
title=oZtitle,/exact,TEXTPOS=0, ZCOORD_CONV=NormCoord(zr,PosZ))

oXAxis->GetProperty, TICKTEXT = xtick_text
oYAxis->GetProperty, TICKTEXT = ytick_text
oZAxis->GetProperty, TICKTEXT = ztick_text
xtick_text->SetProperty, font=oticfont, Recompute_Dimensions=2
ytick_text->SetProperty, font=oticfont, Recompute_Dimensions=2
ztick_text->SetProperty, font=oticfont, Recompute_Dimensions=2

oGroup->Add, oXAxis
oGroup->Add, oYAxis
oGroup->Add, oZAxis

; ----Create some lights----
oLight1 = OBJ_NEW('IDLgrLight', LOCATION=replicate(ZCLIP[0],3), TYPE=1)
oTop->Add, oLight1
oLight2 = OBJ_NEW('IDLgrLight', TYPE=0, INTENSITY=0.5)
oTop->Add, oLight2

; ----Place the model in the view----
oView->Add, oTop

; ----Rotate to standard view for first draw----
oGroup -> SetProperty, TRANSFORM =DefTrans()
if not keyword_set(ayz) then ayz=90
if not keyword_set(bxz) then bxz=90
if not keyword_set(cxy) then cxy=90
if (ayz eq 90) and (bxz eq 90) and $
 (cxy eq 90) then NonOrtM=identity(3) else NonOrtM=NonOrt2Ort(ayz,bxz,cxy)

; ----Create a trackball----
tbCenter=[xdim/2, ydim/2.]
tbRadius=xdim/2.
oTrack = OBJ_NEW('trackball',tbCenter,tbRadius)

; ----Create some fonts----
oLabelFont1= OBJ_NEW('IDLgrFont','times',size=5)

; ----Create a holder object for easy destruction----
oHolder = OBJ_NEW('IDL_Container')
oHolder->Add, oView
oHolder->Add, oTrack
oHolder->Add, oXtitle
oHolder->Add, oYtitle
oHolder->Add, oZtitle
oHolder->Add, otitlefont
oHolder->Add, oticfont
oHolder->Add, oLabelFont1

; ----Save state----
tmp = oGroup->IDL_Container::Get(/all, COUNT=noGroupFix)

CD, CURRENT=path
list = {list2:list2,$
        xs:0L,$
        ys:0L,$
        myview:myview,$
        draw: draw,$
        path:path,$
        bAxisHide:0b,$

        oHolder: oHolder,$ ; Pointer to objects
        oXAxis:oXAxis,$
        oYAxis:oYAxis,$
        oZAxis:oZAxis,$
        oTrack:oTrack,$
        oWindow: oWindow,$
        oView: oView,$
        oGroup: oGroup,$
        noGroupFix:noGroupFix,$ ; the n first kids are non-user
        oZtitle:oZtitle,$
        oLabelFont1:oLabelFont1,$

        TLBsize:TLBsize,$
        yxratio:ydim/float(xdim),$

        tbCenter:tbCenter,$
        tbRadius:tbRadius,$
        bDataConstr:0b,$
        startTransform:identity(4),$

        NonOrtM:NonOrtM,$
        PosX:PosX,$
        PosY:PosY,$
        PosZ:PosZ,$
        bstate:''}

; ----Add user defined objects----
tmp={list:list,key:1}
call_procedure,list2.proplot,tmp
list=tmp.list

WIDGET_CONTROL, base, SET_UVALUE=list, /NO_COPY

oWindow->draw, oView
Xmanager,'flplot_obj',base,/NO_BLOCK,cleanup='CleanUp_flplot_obj',event_handler='flplot_obj_event',GROUP_LEADER=group_leader

end; pro flplot_obj
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



;pro Example_Event,ev
;
;case widget_info(ev.id,/type) of
; 1 :begin
;     ID=widget_info(ev.top,FIND_BY_UNAME='fplotdrawwidget')
;     WIDGET_CONTROL, ID,xsize=600
;     endcase
; 0:    begin
;     if (widget_info(ev.id,/type)eq 0) then flplot_obj_event,ev
;     endcase
;endcase
;
;end; pro Example_Event
;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
;
;pro Example,list
;zData = BESELJ(SHIFT(DIST(40),20,20)/2,0)
;
;child=1
;
;if (child eq 1) then begin
;    base=widget_base(/column,title='3D rep.',/TLB_SIZE_EVENTS)
;    button=widget_button(base,value='test')
;    base3=widget_base(base,/row)
;    ;base3=base
;        button=widget_button(base3,value='test')
;        base2=widget_base(base3,/column)
;        button=widget_button(base3,value='test')
;        button=widget_button(base3,value='test')
;    button=widget_button(base,value='test')
;    button=widget_button(base,value='test')
;    WIDGET_CONTROL, base, /REALIZE
;    Xmanager,'Example',base,/NO_BLOCK,event_handler='Example_event'
;endif
;
;list={z:zData,zr:[min(zData),max(zData)],xr:[0,39],yr:[0,39]}
;
;if (child eq 1) then flplot_obj,list,parentbase=base2 $
;else flplot_obj,list
;
;end;pro Example
;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
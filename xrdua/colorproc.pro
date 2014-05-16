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

pro ConvertBlackToWhite,r,g,b,val
; This does not give the required result
;COLOR_CONVERT,r,g,b,h,s,v,/RGB_HSV
;COLOR_CONVERT,h,v,s,r,g,b,/HSV_RGB

COLOR_CONVERT,r,g,b,h,s,v,/RGB_HSV
t=250
ind=where(v le val,ct)
if ct ne 0 then begin
    r[ind]=255
    g[ind]=255
    b[ind]=255
endif
end;pro ConvertBlackToWhite
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function gencolor,i,n
hue=360./n*i
sat=1.
val=1.
color_convert,hue,sat,val,R,G,B,/HSV_RGB

if G gt 250 and B gt 250 and R lt 5 then begin
    ; Yellow -> Red
    G=255-G
    B=255-B
    R=255
endif

return,ishft(long(B),16) or ishft(long(G),8) or long(R)
end;function gencolor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro getnamedcolors,names,colors

names=[     'ALICE_BLUE',         'ANTIQUE_WHITE',                  'AQUA',            'AQUAMARINE',                 'AZURE',$
                 'BEIGE',                'BISQUE',                 'BLACK',       'BLANCHED_ALMOND',                  'BLUE',$
           'BLUE_VIOLET',                 'BROWN',             'BURLYWOOD',            'CADET_BLUE',            'CHARTREUSE',$
             'CHOCOLATE',                 'CORAL',            'CORNFLOWER',              'CORNSILK',               'CRIMSON',$
                  'CYAN',             'DARK_BLUE',             'DARK_CYAN',        'DARK_GOLDENROD',             'DARK_GRAY',$
            'DARK_GREEN',             'DARK_GREY',            'DARK_KHAKI',          'DARK_MAGENTA',      'DARK_OLIVE_GREEN',$
           'DARK_ORANGE',           'DARK_ORCHID',              'DARK_RED',           'DARK_SALMON',        'DARK_SEA_GREEN',$
       'DARK_SLATE_BLUE',       'DARK_SLATE_GRAY',       'DARK_SLATE_GREY',        'DARK_TURQUOISE',           'DARK_VIOLET',$
             'DEEP_PINK',         'DEEP_SKY_BLUE',              'DIM_GRAY',              'DIM_GREY',           'DODGER_BLUE',$
             'FIREBRICK',          'FLORAL_WHITE',          'FOREST_GREEN',               'FUCHSIA',             'GAINSBORO',$
           'GHOST_WHITE',                  'GOLD',             'GOLDENROD',                  'GRAY',                 'GREEN',$
          'GREEN_YELLOW',                  'GREY',              'HONEYDEW',              'HOT_PINK',            'INDIAN_RED',$
                'INDIGO',                 'IVORY',                 'KHAKI',              'LAVENDER',        'LAVENDER_BLUSH',$
            'LAWN_GREEN',         'LEMON_CHIFFON',            'LIGHT_BLUE',           'LIGHT_CORAL',            'LIGHT_CYAN',$
       'LIGHT_GOLDENROD',           'LIGHT_GREEN',            'LIGHT_GRAY',            'LIGHT_GREY',            'LIGHT_PINK',$
          'LIGHT_SALMON',       'LIGHT_SEA_GREEN',        'LIGHT_SKY_BLUE',      'LIGHT_SLATE_GRAY',      'LIGHT_SLATE_GREY',$
      'LIGHT_STEEL_BLUE',          'LIGHT_YELLOW',                  'LIME',            'LIME_GREEN',                 'LINEN',$
               'MAGENTA',                'MAROON',     'MEDIUM_AQUAMARINE',           'MEDIUM_BLUE',         'MEDIUM_ORCHID',$
         'MEDIUM_PURPLE',      'MEDIUM_SEA_GREEN',     'MEDIUM_SLATE_BLUE',   'MEDIUM_SPRING_GREEN',      'MEDIUM_TURQUOISE',$
     'MEDIUM_VIOLET_RED',         'MIDNIGHT_BLUE',            'MINT_CREAM',            'MISTY_ROSE',              'MOCCASIN',$
          'NAVAJO_WHITE',                  'NAVY',              'OLD_LACE',                 'OLIVE',            'OLIVE_DRAB',$
                'ORANGE',            'ORANGE_RED',                'ORCHID',        'PALE_GOLDENROD',            'PALE_GREEN',$
        'PALE_TURQUOISE',       'PALE_VIOLET_RED',           'PAPAYA_WHIP',            'PEACH_PUFF',                  'PERU',$
                  'PINK',                  'PLUM',           'POWDER_BLUE',                'PURPLE',                   'RED',$
            'ROSY_BROWN',            'ROYAL_BLUE',          'SADDLE_BROWN',                'SALMON',           'SANDY_BROWN',$
             'SEA_GREEN',              'SEASHELL',                'SIENNA',                'SILVER',              'SKY_BLUE',$
            'SLATE_BLUE',            'SLATE_GRAY',            'SLATE_GREY',                  'SNOW',          'SPRING_GREEN',$
            'STEEL_BLUE',                   'TAN',                  'TEAL',               'THISTLE',                'TOMATO',$
             'TURQUOISE',                'VIOLET',                 'WHEAT',                 'WHITE',           'WHITE_SMOKE',$
                'YELLOW',          'YELLOW_GREEN']


colors=[[240b,248b,255b],[250b,235b,215b],[  0b,255b,255b],[127b,255b,212b],[240b,255b,255b],$
        [245b,245b,220b],[255b,228b,196b],[  0b,  0b,  0b],[255b,235b,205b],[  0b,  0b,255b],$
        [138b, 43b,226b],[165b, 42b, 42b],[222b,184b,135b],[ 95b,158b,160b],[127b,255b,  0b],$
        [210b,105b, 30b],[255b,127b, 80b],[100b,149b,237b],[255b,248b,220b],[220b, 20b, 60b],$
        [  0b,255b,255b],[  0b,  0b,139b],[  0b,139b,139b],[184b,134b, 11b],[169b,169b,169b],$
        [  0b,100b,  0b],[169b,169b,169b],[189b,183b,107b],[139b,  0b,139b],[ 85b,107b, 47b],$
        [255b,140b,  0b],[153b, 50b,204b],[139b,  0b,  0b],[233b,150b,122b],[143b,188b,143b],$
        [ 72b, 61b,139b],[ 47b, 79b, 79b],[ 47b, 79b, 79b],[  0b,206b,209b],[148b,  0b,211b],$
        [255b, 20b,147b],[  0b,191b,255b],[105b,105b,105b],[105b,105b,105b],[ 30b,144b,255b],$
        [178b, 34b, 34b],[255b,250b,240b],[ 34b,139b, 34b],[255b,  0b,255b],[220b,220b,220b],$
        [248b,248b,255b],[255b,215b,  0b],[218b,165b, 32b],[127b,127b,127b],[  0b,127b,  0b],$
        [173b,255b, 47b],[127b,127b,127b],[240b,255b,240b],[255b,105b,180b],[205b, 92b, 92b],$
        [ 75b,  0b,130b],[255b,255b,240b],[240b,230b,140b],[230b,230b,250b],[255b,240b,245b],$
        [124b,252b,  0b],[255b,250b,205b],[173b,216b,230b],[240b,128b,128b],[224b,255b,255b],$
        [250b,250b,210b],[144b,238b,144b],[211b,211b,211b],[211b,211b,211b],[255b,182b,193b],$
        [255b,160b,122b],[ 32b,178b,170b],[135b,206b,250b],[119b,136b,153b],[119b,136b,153b],$
        [176b,196b,222b],[255b,255b,224b],[  0b,255b,  0b],[ 50b,205b, 50b],[250b,240b,230b],$
        [255b,  0b,255b],[127b,  0b,  0b],[102b,205b,170b],[  0b,  0b,205b],[186b, 85b,211b],$
        [147b,112b,219b],[ 60b,179b,113b],[123b,104b,238b],[  0b,250b,154b],[ 72b,209b,204b],$
        [199b, 21b,133b],[ 25b, 25b,112b],[245b,255b,250b],[255b,228b,225b],[255b,228b,181b],$
        [255b,222b,173b],[  0b,  0b,128b],[253b,245b,230b],[128b,128b,  0b],[107b,142b, 35b],$
        [255b,165b,  0b],[255b, 69b,  0b],[218b,112b,214b],[238b,232b,170b],[152b,251b,152b],$
        [175b,238b,238b],[219b,112b,147b],[255b,239b,213b],[255b,218b,185b],[205b,133b, 63b],$
        [255b,192b,203b],[221b,160b,221b],[176b,224b,230b],[127b,  0b,127b],[255b,  0b,  0b],$
        [188b,143b,143b],[ 65b,105b,225b],[139b, 69b, 19b],[250b,128b,114b],[244b,164b, 96b],$
        [ 46b,139b, 87b],[255b,245b,238b],[160b, 82b, 45b],[192b,192b,192b],[135b,206b,235b],$
        [106b, 90b,205b],[112b,128b,144b],[112b,128b,144b],[255b,250b,250b],[  0b,255b,127b],$
        [ 70b,130b,180b],[210b,180b,140b],[  0b,128b,128b],[216b,191b,216b],[255b, 99b, 71b],$
        [ 64b,224b,208b],[238b,130b,238b],[245b,222b,179b],[255b,255b,255b],[245b,245b,245b],$
        [255b,255b,  0b],[154b,205b, 50b]]

end;pro getnamedcolors
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function namedcolors,name
CATCH, Error_status
IF Error_status NE 0 THEN return,bytarr(3)
getnamedcolors,names,colors
if size(name,/type) ne 7 then return,colors[*,name]
ind=where(names eq strupcase(name))
return,colors[*,ind[0]]
end;function namedcolors
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotcolorpickerSV,hin,xsize,ysize,_extra=ex
h=make_array(xsize,ysize,value=hin)
s=rebin(findgen(xsize)/(xsize-1),xsize,ysize,/sample)
v=rebin(findgen(1,ysize)/(ysize-1),xsize,ysize,/sample)
COLOR_CONVERT,h,s,v,r,g,b,/HSV_RGB
tv,[[[r]],[[g]],[[b]]],true=3,_extra=ex
end;pro plotcolorpickerSV
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotcolorpickerH,xsize,ysize,_extra=ex
h=rebin(findgen(1,ysize)/(ysize-1)*360,xsize,ysize,/sample)
s=make_array(xsize,ysize,value=1)
v=make_array(xsize,ysize,value=1)
COLOR_CONVERT,h,s,v,r,g,b,/HSV_RGB
tv,[[[r]],[[g]],[[b]]],true=3,_extra=ex
end;pro plotcolorpickerH
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_ColorPicker,ID

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
end;pro CleanUp_ColorPicker
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshColorPicker,ev,init=init
binit=keyword_set(init)
if binit then begin
    list={h:init.h,s:init.s,v:init.v,pColor:init.pColor}
endif else begin
    widget_control,ev.top,get_uvalue=list
endelse

; Show color
if binit then begin
    ID=init.draw
endif else begin
    ID=widget_info(ev.top,FIND_BY_UNAME='draw')
endelse
widget_control,ID,get_value=drawid
wset,drawid
img=bytarr(!d.x_size,!d.y_size,3)
img[*,*,0]=(*list.pColor)[0]
img[*,*,1]=(*list.pColor)[1]
img[*,*,2]=(*list.pColor)[2]
tv,img,true=3

if binit then begin
    ID=init.drawSV
endif else begin
    ID=widget_info(ev.top,FIND_BY_UNAME='drawSV')
endelse
widget_control,ID,get_value=drawid
wset,drawid
plotcolorpickerSV,list.h,!d.x_size,!d.y_size

plots,[0,!d.x_size-1],list.v*(!d.y_size-1)*[1,1],/device,color=[0,0,0]
plots,list.s*(!d.x_size-1)*[1,1],[0,!d.y_size-1],/device,color=[0,0,0]

if binit then begin
    ID=init.drawH
endif else begin
    ID=widget_info(ev.top,FIND_BY_UNAME='drawH')
endelse
widget_control,ID,get_value=drawid
wset,drawid
plotcolorpickerH,!d.x_size,!d.y_size

plots,[0,!d.x_size-1],list.h/360*(!d.y_size-1)*[1,1],/device,color=[0,0,0]

end;pro RefreshColorPicker
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ColorPickerCall,list
if list.refreshname ne '' then call_procedure,list.refreshname,list.refreshvar1,colorpicker=list
end;pro ColorPickerCall
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ColorPicker_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list

type=widget_info(ev.id,/type)
case type of
1:    begin
    widget_control,ev.id,get_value=val
    if val eq 'Cancel' then begin
        *list.pColor=0
        ColorPickerCall,list
    endif
    widget_control,ev.top,/destroy
    return
    endcase
6:    begin
    *list.pColor=namedcolors(ev.index)
    ColorPickerCall,list
    endcase
4:    begin
    widget_control,ev.id,get_uvalue=uval
    case uval of
    'drawH':    begin
                ID=widget_info(ev.top,FIND_BY_UNAME='drawH')
                widget_control,ID,get_value=drawid
                wset,drawid
                list.h=ev.y/(!d.y_size-1.)*360
                widget_control,ev.top,set_uvalue=list
                endcase
    'drawSV':    begin
                ID=widget_info(ev.top,FIND_BY_UNAME='drawSV')
                widget_control,ID,get_value=drawid
                wset,drawid
                h=list.h
                s=ev.x/(!d.x_size-1.)
                v=ev.y/(!d.y_size-1.)
                list.s=s
                list.v=v
                widget_control,ev.top,set_uvalue=list
                COLOR_CONVERT,h,s,v,r,g,b,/HSV_RGB
                *list.pColor=[r,g,b]
                ColorPickerCall,list
                endcase
    endcase
    endcase
else:
endcase

RefreshColorPicker,ev
end;pro ColorPicker_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ColorPicker,top,refreshname=refreshname,refreshvar1=refreshvar1,colorinit=colorinit

if n_elements(colorinit) ne 3 then colorinit=[255b,255b,255b]

COLOR_CONVERT,colorinit[0],colorinit[1],colorinit[2],h,s,v,/RGB_HSV

pColor=PTR_NEW(colorinit)
if WIDGET_INFO(top,/VALID_ID) eq 1 then begin
    WIDGET_CONTROL, top, sensitive=0
    GROUP_LEADER=top
    modal=1
endif else begin
    modal=keyword_set(modal)
endelse

if not keyword_set(refreshname) then refreshname=''
if not keyword_set(refreshname) then refreshvar1=0
list={top:top,pColor:pColor,h:h,s:s,v:v,refreshname:refreshname,refreshvar1:refreshvar1}

device,get_screen_size=screen
base=widget_base(title='Color Picker',/row,uvalue=list,$
                GROUP_LEADER=GROUP_LEADER,modal=modal,xoffset=screen[0]/2.,yoffset=screen[1]/2.)

tab=widget_tab(base)

base2=widget_base(tab,/row,title='Picker')
drawSV=widget_draw(base2,xsize=200,ysize=200,uname='drawSV',uvalue='drawSV',/BUTTON_EVENTS)
drawH=widget_draw(base2,xsize=50,ysize=200,uname='drawH',uvalue='drawH',/BUTTON_EVENTS)

base1=widget_base(tab,/row,title='Named')
getnamedcolors,names
drop=widget_list(base1,value=names,uname='colors',xsize=50,ysize=10)

baset=widget_base(base,/column)
draw=widget_draw(baset,xsize=100,ysize=50,uname='draw')

b=widget_button(baset,value='Cancel')
b=widget_button(baset,value='OK')

widget_control,base,/realize

RefreshColorPicker,dummy,init={draw:draw,drawH:drawH,drawSV:drawSV,h:list.h,s:list.s,v:list.v,pColor:list.pColor}

widget_control,b,/INPUT_FOCUS
Xmanager,'ColorPicker',base, event_handler='ColorPicker_Event',$
    cleanup='CleanUp_ColorPicker',GROUP_LEADER=GROUP_LEADER

out=(*pColor)
PTR_FREE,pColor
return,out
end;function ColorPicker
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GradientFillMissing,R,G,B,rgbback,niter

filter=make_array(3,3,value=1/8.)
filter[1,1]=0

ind=where(R eq rgbback[0] and G eq rgbback[1] and B eq rgbback[2],ct)
if ct eq 0 then return

for i=0,niter-1 do begin
    ;if ct ne 0 then begin
        R[ind]=round((convol(R,filter,/edge_zero))[ind])<255
        G[ind]=round((convol(G,filter,/edge_zero))[ind])<255
        B[ind]=round((convol(B,filter,/edge_zero))[ind])<255
    ;endif
endfor
end;pro GradientFillMissing
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AlphaBlendingWithGradient,R,G,B,A,intensity_,gradient,rgbback=rgbback,normal=normal

; Background
if n_elements(R) eq 0 or n_elements(G) eq 0 or n_elements(B) eq 0 or n_elements(A) eq 0 then begin
    ; Background
    s=dimsize(intensity_,2)
    rgbback=gradient.rgbback
    R = make_array(s,value=rgbback[0])
    G = make_array(s,value=rgbback[1])
    B = make_array(s,value=rgbback[2])
    A = gradient.opaqback ; opaque
endif

; Scale to [0,1]
intensity=intensity_
intensity-=min(intensity)
intensity/=max(intensity)

; Find clipped pixels
indlow=where(intensity lt gradient.cliplow,ncliplow)
indhigh=where(intensity gt gradient.cliphigh,ncliphigh)

; Clip
intensity-=gradient.cliplow
intensity/=gradient.cliphigh-gradient.cliplow
intensity>=0
intensity<=1

; Gamma
gamma=alog(0.5)/alog(gradient.gammap)
gammaopaq=alog(0.5)/alog(gradient.gammapopaq)

; Linear interpolation of RGB colors with gamma correction
R_ = gradient.rgblow[0] + intensity^gamma*(long(gradient.rgbhigh[0])-gradient.rgblow[0])
G_ = gradient.rgblow[1] + intensity^gamma*(long(gradient.rgbhigh[1])-gradient.rgblow[1])
B_ = gradient.rgblow[2] + intensity^gamma*(long(gradient.rgbhigh[2])-gradient.rgblow[2])

; Linear interpolation of opacity with gamma correction
A_ = gradient.opaqlow + intensity^gammaopaq*(gradient.opaqhigh-gradient.opaqlow)

; Opacity of clipped pixels
if ncliplow ne 0 then A_[indlow]=gradient.opaqcliplow
if ncliphigh ne 0 then A_[indhigh]=gradient.opaqcliphigh

; Blend using alpha channels
; https://dvcs.w3.org/hg/FXTF/rawfile/default/compositing/index.html#blendingnormal
; 
; Formula:
;    Anew.Rnew = A_.(1-A).R_ + A_.A.f(R,R_) + A.(1-A_).R

Anew = 1-(1-A)*(1-A_)
bapply = 1b
mode=1
mode=gradient.mode
if keyword_set(normal) then mode=0
case mode of
1:    begin ; addition
    bapply = 0b
    R = ((A_*R_ + A*R)/Anew)<255
    G = ((A_*G_ + A*G)/Anew)<255
    B = ((A_*B_ + A*B)/Anew)<255
    endcase
2:    begin; difference
    ; f(R,R_) = abs(R-R_)
    fR=abs(R-R_)
    fG=abs(G-G_)
    fB=abs(B-B_)
    endcase
3:    begin ; test
    bapply = 1b
    fR=R>R_
    fG=G>G_
    fB=B>B_
    endcase
else:begin; normal
    ; f(R,R_) = R_
    bapply = 0b
    R = (A_*R_ + A*(1-A_)*R)/Anew
    G = (A_*G_ + A*(1-A_)*G)/Anew
    B = (A_*B_ + A*(1-A_)*B)/Anew
    endelse
endcase

if keyword_set(bapply) then begin
    R = (A_*(1-A)*R_ + A_*A*fR + A*(1-A_)*R)/Anew
    G = (A_*(1-A)*G_ + A_*A*fG + A*(1-A_)*G)/Anew
    B = (A_*(1-A)*B_ + A_*A*fB + A*(1-A_)*B)/Anew
endif

A = temporary(Anew)

end;pro AlphaBlendingWithGradient
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function blendmodes
return,['Normal','Addition','Difference']
end;function blendmodes
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DefaultGradient,type
; All floating point values between 0 and 1
; All byte values between 0 and 255

if type eq -1 then col=[255b,255b,255b] else col=namedcolors(type)

return,{rgblow:[0b,0b,0b],rgbhigh:col,rgbback:[0b,0b,0b],$
    opaqlow:1.,opaqhigh:1.,opaqback:1.,$
    gammap:0.5,gammapopaq:0.5,$
    cliplow:0.01,cliphigh:1.,opaqcliplow:0.,opaqcliphigh:1.,$
    buse:0b,index:-1l,mode:1}
end;function DefaultGradient
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Defaultgradients
ret=ptrarr(3,/alloc)
*ret[0]=DefaultGradient(119)
*ret[1]=DefaultGradient(82)
*ret[2]=DefaultGradient(9)

return,ptr_new(ret)
end;function Defaultgradients
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ResizeGradients,pgradients,n,pgradient
n>=1
n=n[0]
ngrad=n_elements(*pgradients)
if ngrad eq n then return
if ngrad lt n then begin
    add=ptrarr(n-ngrad,/alloc)
    for i=0l,n-ngrad-1 do *add[i]=DefaultGradient(-1)
    *pgradients=[(*pgradients),add]
endif
if ngrad gt n then begin
    heap_free,(*pgradients)[n:*]
    (*pgradients)=(*pgradients)[0:n-1]
    if not ptr_valid(*pgradient) then *pgradient=(*pgradients)[0]
endif
end;pro ResizeGradients,pgradients
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MoveGradientLayer,list,move

ind=where(*list.pgradients eq *list.pgradient,ct)
if ct ne 1 then return
ind=ind[0]
n=n_elements(*list.pgradients)

ind=MoveElementInArray(n,ind,ind+move)

*list.pgradients=(*list.pgradients)[ind]

end;pro MoveGradientLayer
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SaveGradients,list
; Select file
path=list.path
tmp=CutPath(list.file,file=file)
path=DIALOG_PICKFILE(path=path,file=file+'.bld',filter='*.bld',title='Save Blend Settings ...',/OVERWRITE_PROMPT,/WRITE)
if path eq '' then return

arr1=*list.pgradients
arr2=list.gradientfill
save,arr1,arr2,filename=path
end;pro SaveGradients
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RestoreGradients,list
; Select file
path=list.path
tmp=CutPath(list.file,file=file)
path=DIALOG_PICKFILE(path=path,file=file+'.bld',filter='*.bld',title='Load Blend Settings ...')
if path eq '' then return

ind=where(*list.pgradients eq *list.pgradient,ct)

restore,filename=path

ind=where(tag_names(*(arr1[0])) eq 'MODE',ct)
if ct eq 0 then $
  for i=0,n_elements(arr1)-1 do $
    *(arr1[i])=create_struct(*(arr1[i]),{mode:0})

heap_free,*list.pgradients
*list.pgradients=arr1

ind>=0
ind<=n_elements(*list.pgradients)-1
*list.pgradient=(*list.pgradients)[ind[0]]

list.gradientfill=arr2
end;pro RestoreGradients
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SetGradientTableValues,pgradient
values=[transpose((**pgradient).rgblow),transpose((**pgradient).rgbhigh),transpose((**pgradient).rgbback)]
values=[[values],[[(**pgradient).opaqlow,(**pgradient).opaqhigh,(**pgradient).opaqback]*100]]
values=[[values],[[(**pgradient).cliplow,(**pgradient).cliphigh,!values.F_NAN]*100]]
values=[[values],[[(**pgradient).opaqcliplow,(**pgradient).opaqcliphigh,!values.F_NAN]*100]]
return,values
end;function SetGradientTableValues
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GetGradientTableValues,ev,pgradient
widget_control,ev.id,get_value=values

(**pgradient).rgblow=0>round(transpose(values[0,0:2]))<255
(**pgradient).rgbhigh=0>round(transpose(values[1,0:2]))<255
(**pgradient).rgbback=0>round(transpose(values[2,0:2]))<255

(**pgradient).opaqlow=0>values[0,3]/100.<1
(**pgradient).opaqhigh=0>values[1,3]/100.<1
(**pgradient).opaqback=0>values[2,3]/100.<1

(**pgradient).cliplow=0>values[0,4]/100.<1
(**pgradient).cliphigh=0>values[1,4]/100.<1

(**pgradient).opaqcliplow=0>values[0,5]/100.<1
(**pgradient).opaqcliphigh=0>values[1,5]/100.<1

end;pro GetGradientTableValues
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SelectClipMarker,ev,xs,ys,list

fracclipmarker=[0.025,4/30.]
nx=round(fracclipmarker[0]*xs)
ny=round(fracclipmarker[1]*ys)

; Low intensity marker slected?
x=[-1,1]*nx/2.+(**list.pgradient).cliplow*xs
y=[0,ny-1]
if ev.x ge x[0] and ev.x le x[1] and ev.y ge y[0] and ev.y le y[1] then return,1 

; High intensity marker selected?
x=[-1,1]*nx/2.+(**list.pgradient).cliphigh*xs
y=[0,ny-1]
if ev.x ge x[0] and ev.x le x[1] and ev.y ge y[0] and ev.y le y[1] then return,2

; Gamma marker selected?
x0=(**list.pgradient).cliplow*xs
x1=(**list.pgradient).cliphigh*xs
x=x0+(**list.pgradient).gammap*(x1-x0)+[-5,5]
y=ny+[-5,5]
if ev.x ge x[0] and ev.x le x[1] and ev.y ge y[0] and ev.y le y[1] then return,3

; Low opacity marker selected?
x2=x0+(**list.pgradient).opaqlow*(x1-x0)
x=[-1,1]*nx/2.+x2
y=ys-1-[ny-1,0]
if ev.x ge x[0] and ev.x le x[1] and ev.y ge y[0] and ev.y le y[1] then return,4

; High opacity marker selected?
x3=x0+(**list.pgradient).opaqhigh*(x1-x0)
x=[-1,1]*nx/2.+x3
y=ys-1-[ny-1,0]
if ev.x ge x[0] and ev.x le x[1] and ev.y ge y[0] and ev.y le y[1] then return,5

; Opacity Gamma marker selected?
x=x2+(**list.pgradient).gammapopaq*(x3-x2)+[-5,5]
y=ys-1-ny+[-5,5]
if ev.x ge x[0] and ev.x le x[1] and ev.y ge y[0] and ev.y le y[1] then return,6

return,0
end;function SelectClipMarker
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GradientClipMarker,xs,ys,x,y,boarder=boarder
fracclipmarker=[0.025,4/30.]
nx=round(fracclipmarker[0]*xs)
nx+=1-(nx mod 2) ; make odd
ny=round(fracclipmarker[1]*ys)

image=fltarr(nx,ceil(1.5*ny))

if keyword_set(boarder) then begin
    image[0:nx-1,0]=1
    image[0:nx-1,ny-1]=1
    image[0,0:ny-1]=1
    image[nx-1,0:ny-1]=1
    
    x=[0,(nx-1)/2.]
    y=[ny,1.5*ny]
    m=(y[1]-y[0])/(x[1]-x[0])
    b=ny
    x=lindgen(x[1]-x[0])+x[0]
    y=m*x+b
    image[x,y]=1
    
    x=[(nx-1)/2.,nx-1]
    y=[1.5*ny,ny]
    m=(y[1]-y[0])/(x[1]-x[0])
    b=2*ny
    x=lindgen(x[1]-x[0]+1)+x[0]
    y=m*x+b
    image[x,y]=1
    
;    window & tvscl,rebin(image,nx*10,ceil(1.5*ny)*10,/sample)
endif else begin
    image[0:nx-1,0:ny-1]=1
endelse

ind=where(image)
x=ind mod nx
y=ind/nx
x-=nx/2

end;pro GradientClipMarker
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PlotGradientImage,pgradient

; Fractions of the x and y size
xs=!d.x_size
ys=!d.y_size

; Color band
fracband=0.6
ysband=fracband*ys
image=[[fltarr(xs,floor((ys-ysband)/2.))],[rebin(findgen(xs)/(xs-1),xs,ysband,/sample)],[fltarr(xs,ceil((ys-ysband)/2.))]]
AlphaBlendingWithGradient,R,G,B,A,image,**pgradient,/normal
tmp=size(image,/dim)
xs=tmp[0]
ys=tmp[1]

; Clip markers
gradientmarker=DefaultGradient(-1)

for i=0,1 do begin
    boarder=i eq 1
    GradientClipMarker,xs,ys,x,y,boarder=boarder
    
    xx=x+(**pgradient).cliplow*xs
    ind=where(xx ge 0 and x lt xs,ct)
    if ct ne 0 then begin
        image=fltarr(xs,ys)
        image[xx[ind],y[ind]]=1
        gradientmarker.rgbhigh=boarder?(255b-(**pgradient).rgbback):(**pgradient).rgblow
        gradientmarker.opaqhigh=1
        AlphaBlendingWithGradient,R,G,B,A,image,gradientmarker,/normal
    endif
    
    xx=x+(**pgradient).cliphigh*xs
    yy=y
    ind=where(xx ge 0 and x lt xs,ct)
    if ct ne 0 then begin
        image=fltarr(xs,ys)
        image[xx[ind],yy[ind]]=1
        gradientmarker.rgbhigh=boarder?(255b-(**pgradient).rgbback):(**pgradient).rgbhigh
        gradientmarker.opaqhigh=1
        AlphaBlendingWithGradient,R,G,B,A,image,gradientmarker,/normal
    endif
    
    x0=(**pgradient).cliplow*xs
    x1=(**pgradient).cliphigh*xs
    xx=x+x0+(**pgradient).opaqlow*(x1-x0)
    ind=where(xx ge 0 and x lt xs,ct)
    if ct ne 0 then begin
        image=fltarr(xs,ys)
        image[xx[ind],ys-1-y[ind]]=1
        gradientmarker.rgbhigh=boarder?(255b-(**pgradient).rgbback):(**pgradient).rgblow
        gradientmarker.opaqhigh=boarder?100:(**pgradient).opaqlow
        AlphaBlendingWithGradient,R,G,B,A,image,gradientmarker,/normal
    endif
    
    x0=(**pgradient).cliplow*xs
    x1=(**pgradient).cliphigh*xs
    xx=x+x0+(**pgradient).opaqhigh*(x1-x0)
    ind=where(xx ge 0 and x lt xs,ct)
    if ct ne 0 then begin
        image=fltarr(xs,ys)
        image[xx[ind],ys-1-y[ind]]=1
        gradientmarker.rgbhigh=boarder?(255b-(**pgradient).rgbback):(**pgradient).rgbhigh
        gradientmarker.opaqhigh=boarder?100:(**pgradient).opaqhigh
        AlphaBlendingWithGradient,R,G,B,A,image,gradientmarker,/normal
    endif
endfor

tv,[[[r]],[[g]],[[b]]],true=3

; Plot gamma
device,decompos=1
x0=(**pgradient).cliplow*xs
x1=(**pgradient).cliphigh*xs
x=x0+(**pgradient).gammap*(x1-x0)
y=4/30.*ys
plots,x,y,psym=4,/device,color=255b-(**pgradient).rgbback

; Plot opacity gamma
x2=x0+(**pgradient).opaqlow*(x1-x0)
x3=x0+(**pgradient).opaqhigh*(x1-x0)
x=x2+(**pgradient).gammapopaq*(x3-x2)
y=ys-1-4/30.*ys
plots,x,y,psym=4,/device,color=255b-(**pgradient).rgbback
device,decompos=0
end;pro PlotGradientImage
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GradientWidgetCall,list
if list.gradientrefreshname ne '' then call_procedure,list.gradientrefreshname,list,/blendupdate
end;pro GradientWidgetCall
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshGradientImage,ev,colorpicker=colorpicker,save=save

widget_control,ev.top,get_uvalue=list

if keyword_set(colorpicker) then begin
    if list.gradientclipev eq 1 then (**list.pgradient).rgblow=*colorpicker.pcolor $
    else (**list.pgradient).rgbhigh=*colorpicker.pcolor
endif

; Plot gradient
Wset,list.drawindexgradient
PlotGradientImage,list.pgradient

; Save as image
if keyword_set(save) then begin
    path=list.path
    tmp=CutPath(list.file,file=file)
    path=DIALOG_PICKFILE(path=path,file=file+'.bmp',filter='*.bmp',title='Save Gradient ...',/OVERWRITE_PROMPT,/WRITE)
    if path eq '' then return

    img=tvrd(true=1)
    write_bmp,path,img,/rgb
endif

; Copy to pixmap
Wset,list.pixindexgradient
Device, Copy = [0,0, !d.x_size, !d.y_size, 0,0,list.drawindexgradient]

; Refresh table
ID=widget_info(ev.top,FIND_BY_UNAME='gradienttable')
widget_control,ID,set_value=SetGradientTableValues(list.pgradient)

; Refresh main window
GradientWidgetCall,list

end;pro RefreshGradientImage
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GradientWidgetEventDraw,ev

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
possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
thisEvent = possibleButtons[ev.press]

widget_control,ev.top,get_uvalue=list
WSet, list.drawindexgradient
list.gradientclipev=SelectClipMarker(ev,!d.x_size,!d.y_size,list)
widget_control,ev.top,set_uvalue=list
if list.gradientclipev eq 0 then return

IF thisEvent EQ 'RIGHT' THEN BEGIN
    if list.gradientclipev eq 1 or list.gradientclipev eq 2 then begin
        if list.gradientclipev eq 1 then colorkeep=(**list.pgradient).rgblow $
        else colorkeep=(**list.pgradient).rgbhigh
        
        color=ColorPicker(ev.top,refreshname='RefreshGradientImage',refreshvar1=ev,colorinit=colorkeep)
        if n_elements(color) ne 3 then color=colorkeep
        
        if list.gradientclipev eq 1 then (**list.pgradient).rgblow=color $
        else (**list.pgradient).rgbhigh=color
        
        RefreshGradientImage,ev
    endif
ENDIF

IF thisEvent EQ 'LEFT' THEN BEGIN
    Widget_Control, ev.id, Event_Pro='GradientWidgetEventMotion',draw_motion=1
ENDIF
end;pro GradientWidgetEventDraw
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GradientWidgetEventMotion,ev
if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
WSet, list.drawindexgradient

case list.gradientclipev of
1:    begin
    clipnew=0>(ev.x/(!d.x_size-1.))<1
    (**list.pgradient).cliplow=clipnew;<(**list.pgradient).cliphigh
    endcase
2:    begin
    clipnew=0>(ev.x/(!d.x_size-1.))<1
    (**list.pgradient).cliphigh=clipnew;>(**list.pgradient).cliplow
    endcase
3:    begin
    x0=(**list.pgradient).cliplow*!d.x_size
    x1=(**list.pgradient).cliphigh*!d.x_size
    gammanew=0>((ev.x-x0)>0)/(x1-x0)<0.999
    (**list.pgradient).gammap=gammanew
    endcase
4:    begin
    x0=(**list.pgradient).cliplow*!d.x_size
    x1=(**list.pgradient).cliphigh*!d.x_size
    clipnew=0>(((ev.x-x0)>0)/(x1-x0))<1
    (**list.pgradient).opaqlow=clipnew
    endcase
5:    begin
    x0=(**list.pgradient).cliplow*!d.x_size
    x1=(**list.pgradient).cliphigh*!d.x_size
    clipnew=0>((ev.x-x0)/(x1-x0))<1
    (**list.pgradient).opaqhigh=clipnew
    endcase
6:    begin
    x0=(**list.pgradient).cliplow*!d.x_size
    x1=(**list.pgradient).cliphigh*!d.x_size
    x2=x0+(**list.pgradient).opaqlow*(x1-x0)
    x3=x0+(**list.pgradient).opaqhigh*(x1-x0)
    gammanew=0>(ev.x-x2)/(x3-x2)<0.999
    (**list.pgradient).gammapopaq=gammanew
    endcase
endcase

; Flip?
if (**list.pgradient).cliplow gt (**list.pgradient).cliphigh then begin
      tmp=(**list.pgradient).cliplow
       (**list.pgradient).cliplow=(**list.pgradient).cliphigh 
       (**list.pgradient).cliphigh =tmp
    
       tmp=(**list.pgradient).opaqlow
       (**list.pgradient).opaqlow=(**list.pgradient).opaqhigh
       (**list.pgradient).opaqhigh=tmp
    
       tmp=(**list.pgradient).rgblow
       (**list.pgradient).rgblow=(**list.pgradient).rgbhigh 
       (**list.pgradient).rgbhigh =tmp
       
       (**list.pgradient).gammap=1-(**list.pgradient).gammap
       (**list.pgradient).gammapopaq=1-(**list.pgradient).gammapopaq
       
       if list.gradientclipev eq 2 then list.gradientclipev=1 else list.gradientclipev=2
       widget_control,ev.top,set_uvalue=list
endif

RefreshGradientImage,ev

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ,'EXPOSE']
thisEvent = possibleEventTypes[ev.type]
IF thisEvent EQ 'UP' THEN BEGIN
    Widget_Control, ev.id, draw_motion=0,Event_Pro='GradientWidgetEventDraw'
    return
ENDIF
       
end;pro GradientWidgetEventMotion
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GradientWidgetEventTable,ev
widget_control,ev.top,get_uvalue=list
GetGradientTableValues,ev,list.pgradient
RefreshGradientImage,ev
end;pro GradientWidgetEventTable
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GradientWidgetEventList,ev
widget_control,ev.top,get_uvalue=list
buse=(**list.pgradient).buse
index=(**list.pgradient).index
**list.pgradient=Defaultgradient(ev.index)
(**list.pgradient).buse=buse
(**list.pgradient).index=index
RefreshGradientImage,ev
end;pro GradientWidgetEventList
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GradientWidget,base,xsize,ysize,list

draw=widget_draw(base,xsize=xsize,ysize=ysize,uname='drawgradient',/BUTTON_EVENTS,EVENT_PRO='GradientWidgetEventDraw')
Window,xsize=xsize,ysize=ysize, /Free, /Pixmap
list.pixindexgradient=!D.Window

tbl=widget_table(base,value=SetGradientTableValues(list.pgradient),$
                column_labels=['Low','High','Back'],$
                row_labels=['R','G','B','Opacity (%)','Clipping (%)','Clipped Opacity (%)'],$
                /editable,X_SCROLL_SIZE=3,Y_SCROLL_SIZE=6,uname='gradienttable',uvalue='gradienttable',EVENT_PRO='GradientWidgetEventTable')

;getnamedcolors,names ; Cfr. Defaultgradients
;wlist=widget_list(base,value=names,ysize=6,EVENT_PRO='GradientWidgetEventList')
end;pro GradientWidget
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
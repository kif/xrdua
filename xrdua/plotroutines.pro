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

pro powplot,p1,p2,pow=pow,ylog=ylog_,yrange=yrange_,max_value=max_value_,min_value=min_value_,$
                data=data,clip=clip_,yerror=yerror,_REF_EXTRA=extra

bploterror=keyword_set(yerror)

if (keyword_set(pow)?(pow eq 1):1b) then begin
    if n_params() eq 1 then begin
        if bploterror then ploterror,p1,yerror,ylog=ylog_,yrange=yrange_,max_value=max_value_,min_value=min_value_,$
                data=data,clip=clip_,_EXTRA=extra $
        else plot,p1,ylog=ylog_,yrange=yrange_,max_value=max_value_,min_value=min_value_,$
                data=data,clip=clip_,_EXTRA=extra 
    endif else begin
        if bploterror then ploterror,p1,p2,yerror,ylog=ylog_,yrange=yrange_,max_value=max_value_,min_value=min_value_,$
                data=data,clip=clip_,_EXTRA=extra $
        else plot,p1,p2,ylog=ylog_,yrange=yrange_,max_value=max_value_,min_value=min_value_,$
                data=data,clip=clip_,_EXTRA=extra
    endelse
    return
endif

bpositive=abs(pow) lt 1

; skip ylog

; Power keywords:

if keyword_set(yrange_) then begin
    yrange=yrange_
    if bpositive then yrange>=0
    yrange^=pow
endif

if keyword_set(max_value_) then begin
    max_value=max_value_
    if bpositive then max_value>=0
    max_value^=pow
endif

if keyword_set(min_value_) then begin
    min_value=min_value_
    if bpositive then min_value>=0
    min_value^=pow
endif

if keyword_set(data) then begin
    if keyword_set(clip_) then begin
        clip=clip_
        if bpositive then clip[[1,3]]>=0
        clip[[1,3]]^=pow
    endif
endif

if n_params() eq 1 then begin
    if bploterror then ploterror,p1^pow,yerror*pow*p1^(pow-1),yrange=yrange,max_value=max_value,min_value=min_value,data=data,clip=clip,_EXTRA=extra $
    else plot,p1^pow,yrange=yrange,max_value=max_value,min_value=min_value,data=data,clip=clip,_EXTRA=extra
endif else begin
    if bploterror then ploterror,p1,p2^pow,yerror*pow*p2^(pow-1),yrange=yrange,max_value=max_value,min_value=min_value,data=data,clip=clip,_EXTRA=extra $
    else plot,p1,p2^pow,yrange=yrange,max_value=max_value,min_value=min_value,data=data,clip=clip,_EXTRA=extra
endelse

end;pro powplot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro powoplot,p1,p2,pow=pow,max_value=max_value_,min_value=min_value_,$
                clip=clip_,_REF_EXTRA=extra

if (keyword_set(pow)?(pow eq 1):1b) then begin
    if n_params() eq 1 then oplot,p1,max_value=max_value_,min_value=min_value_,$
                clip=clip_,_EXTRA=extra $
    else oplot,p1,p2,max_value=max_value_,min_value=min_value_,$
                clip=clip_,_EXTRA=extra
    return
endif

bpositive=abs(pow) lt 1

; Power keywords:

if keyword_set(max_value_) then begin
    max_value=max_value_
    if bpositive then max_value>=0
    max_value^=pow
endif

if keyword_set(min_value_) then begin
    min_value=min_value_
    if bpositive then min_value>=0
    min_value^=pow
endif

if keyword_set(clip_) then begin
    clip=clip_
    if bpositive then clip[[1,3]]>=0
    clip[[1,3]]^=pow
endif


if n_params() eq 1 then oplot,p1^pow,max_value=max_value,min_value=min_value,clip=clip,_EXTRA=extra $
else oplot,p1,p2^pow,max_value=max_value,min_value=min_value,clip=clip,_EXTRA=extra
end;pro powoplot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro powplots,p1,p2,pow=pow,data=data,normal=normal,device=device,clip=clip_,_REF_EXTRA=extra

if (keyword_set(pow)?(pow eq 1):1b) or keyword_set(normal) or keyword_set(device) then begin
    if n_params() eq 1 then plots,p1,data=data,normal=normal,device=device,clip=clip_,_EXTRA=extra $
    else plots,p1,p2,data=data,normal=normal,device=device,clip=clip_,_EXTRA=extra
    return
endif

bpositive=abs(pow) lt 1

; Power keywords:
if keyword_set(clip_) then begin
    clip=clip_
    if bpositive then clip[[1,3]]>=0
    clip[[1,3]]^=pow
endif

if n_params() eq 1 then begin
    p=p1
    p[1,*]^=pow
    plots,p,data=data,normal=normal,device=device,clip=clip,_EXTRA=extra
endif else plots,p1,p2^pow,data=data,normal=normal,device=device,clip=clip,_EXTRA=extra
end;pro powplots
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro powxyouts,p1,p2,str,pow=pow,data=data,normal=normal,device=device,clip=clip_,_REF_EXTRA=extra

if (keyword_set(pow)?(pow eq 1):1b) or keyword_set(normal) or keyword_set(device) then begin
    xyouts,p1,p2,str,data=data,normal=normal,device=device,clip=clip_,_EXTRA=extra
    return
endif

bpositive=abs(pow) lt 1

; Power keywords:
if keyword_set(clip_) then begin
    clip=clip_
    if bpositive then clip[[1,3]]>=0
    clip[[1,3]]^=pow
endif

xyouts,p1,p2^pow,str,data=data,normal=normal,device=device,clip=clip,_EXTRA=extra
end;pro powxyouts
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetBackgroundContrastColor
; Postscript background = white
if !d.name eq 'PS' or !D.NAME EQ 'CGM' then begin
    ; background = white
    TVLCT, R, G, B , /GET
    col=ishft(long(B[[0,255]]),16) or ishft(long(G[[0,255]]),8) or long(R[[0,255]])
    if col[0] eq 16777215 then return,255b else return,0b
endif else begin
    ; !p.BACKGROUND = colorindex
    if !p.BACKGROUND eq 0 then return,255b else return,0b
endelse
end;function GetBackgroundContrastColor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro imgtvscl,image,px,py,xx=x,yy=y,nmag=nmag,$
    xtitle=xtitle,ytitle=ytitle,title=title,order=order,noerase=noerase,$
    xtickinterval=xtickinterval,ytickinterval=ytickinterval,$
    sclmax=sclmax,displaypower=displaypower,$
    xtickname=xtickname,ytickname=ytickname,_REF_EXTRA = _extra

if not keyword_set(nmag) then nmag=1.
if n_elements(nmag) eq 1 then nmag=[nmag,nmag]
order=keyword_set(order)
if n_elements(displaypower) ne 1 then displaypower=1

; Display image
s=DimSize(image,2)
nx=s[0]*nmag[0]
ny=s[1]*nmag[1]
    
color=GetBackgroundContrastColor()
if not keyword_set(noerase) then erase
IF !D.NAME EQ 'PS' or !D.NAME EQ 'CGM' THEN begin
    xper=(1-2.*px)*!d.x_vsize
    yper=(1-2.*py)*!d.y_vsize
    m=xper/nx<yper/ny
    nx*=m
    ny*=m

    if keyword_set(sclmax) then tv, WeakScl(image,sclmax,displaypower),px*!d.x_vsize,py*!d.y_vsize,xsize=nx,ysize=ny,order=order, _EXTRA=_extra $
    else tvscl,image,px*!d.x_vsize,py*!d.y_vsize,xsize=nx,ysize=ny,order=order, _EXTRA=_extra
    if !D.NAME EQ 'PS' then thick=3.5277 ; for a line thickness of 1pt
endif else begin

    if nmag[0] eq long(nmag[0]) and nmag[1] eq long(nmag[1]) then begin
        nx = long(nx)
        ny = long(ny)
        if keyword_set(sclmax) then tv, rebin(WeakScl(image,sclmax,displaypower),nx,ny,/sample),px,py,/norm,order=order, _EXTRA=_extra $
        else tvscl,rebin(image,nx,ny,/sample),px,py,/norm,order=order, _EXTRA=_extra
    endif else begin
        if keyword_set(sclmax) then tv, congrid(WeakScl(image,sclmax,displaypower),nx,ny),px,py,/norm,order=order, _EXTRA=_extra $
        else tvscl,congrid(image,nx,ny),px,py,/norm,order=order, _EXTRA=_extra
    endelse
    
endelse

; Axis values
if n_elements(x) eq s[0] then x=[2*x[0]-x[1],x]+(x[1]-x[0])/2. else $ ; assume equaly spaced
if n_elements(x) ne s[0]+1 then x=lindgen(s[0]+1)-0.5
if n_elements(y) eq s[1] then y=[2*y[0]-y[1],y]+(y[1]-y[0])/2. else $ ; assume equaly spaced
if n_elements(y) ne s[1]+1 then y=lindgen(s[1]+1)-0.5

; Axis positions
loc=[px,py, px+float(nx)/!d.x_vsize, py+float(ny)/!d.y_vsize]

; Y ticks
yrange0=min(y,max=yrange1)
if keyword_set(ytickinterval) then begin
    n=(yrange1-yrange0)/ytickinterval+1
    bname=n_elements(ytickname) eq n
    
    ytickinterval0=ytickinterval
    i=1
    while n gt 40 do begin
        ytickinterval=ytickinterval0*i
        n=(yrange1-yrange0)/ytickinterval+1
        i++
    endwhile
    
    if bname then ytickname=ytickname[indgen(n)*i]
endif
if order then yrange=[yrange1,yrange0] $
else yrange=[yrange0,yrange1]

; X ticks
if keyword_set(xtickinterval) then begin
    xrange0=min(x,max=xrange1)
    bname=n_elements(xtickname) eq n
        
    n=(xrange1-xrange0)/xtickinterval+1
    xtickinterval0=xtickinterval
    i=1
    while n gt 40 do begin
        xtickinterval=xtickinterval0*i
        n=(xrange1-xrange0)/xtickinterval+1
        i++
    endwhile
    
    if bname then ytickname=ytickname[indgen(n)*i]
endif

; Maximum 60 thicknames
if n_elements(xtickname) gt 60 then xtickname=''
if n_elements(ytickname) gt 60 then ytickname=''

; Plot axis
PLOT, x, y, /NODATA, /NOERASE, $
      POSITION=loc,xtickname=xtickname,ytickname=ytickname, $
      XSTYLE=1,YSTYLE=1,color=color,yrange=yrange,$
      title=title,xtitle=xtitle,ytitle=ytitle,/norm,thick=thick,$
      xthick=thick,ythick=thick,xtickinterval=xtickinterval,$
      ytickinterval=ytickinterval, _EXTRA=_extra
axis,loc[0],loc[3],xaxis=-1,/norm,/xs,xthick=thick,color=color,$
                    xtickinterval=xtickinterval,ytickinterval=ytickinterval,$
                    xtickname=xtickname,ytickname=ytickname,_EXTRA=_extra
end;pro imgtvscl
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TVcolorbar,loc,range=range,bar=bar
; loc = xstart,ystart,xstop,ystop ; normal coordinates
; loc = [0.15, 0.80, 0.95, 0.90] ; horizontal
; loc=[0.8,0.15,0.9,0.95] ; vertical

xsize = (loc[2] - loc[0]) * !D.X_VSIZE
ysize = (loc[3] - loc[1]) * !D.Y_VSIZE 
xstart = loc[0] * !D.X_VSIZE
ystart = loc[1] * !D.Y_VSIZE

hor=xsize gt ysize
if hor then begin
    if not keyword_set(bar) then bar = rebin(BINDGEN(256),256,2,/sample)
    xticklen=!P.TICKLEN*5
    yticklen=!P.TICKLEN/100000.
    ytickinterval=1
    ytickname=[' ',' ']
    xrange=range
endif else begin
    if not keyword_set(bar) then bar = rebin(BINDGEN(1,256),2,256,/sample)
    xticklen=!P.TICKLEN/100000.
    yticklen=!P.TICKLEN*5
    xtickname=[' ',' ']
    xtickinterval=1
    yrange=range
endelse

color=GetBackgroundContrastColor()
erase
IF !D.NAME EQ 'PS' or !D.NAME EQ 'CGM' THEN begin
    TV, bar, xstart, ystart, XSIZE=xsize, YSIZE=ysize
    if !D.NAME EQ 'PS' then thick=3.5277 ; for a line thickness of 1pt
endif else begin
    TV, CONGRID(bar, xsize+1, ysize+1), xstart, ystart
endelse

plot,[0,1],position=loc,/normal,/nodata,/noerase,xticklen=xticklen,yticklen=yticklen,$
    xtickinterval=xtickinterval,ytickinterval=ytickinterval,xtickname=xtickname,$
    ytickname=ytickname,color=color,xrange=xrange,yrange=yrange,/xs,/ys,thick=thick,xthick=thick,ythick=thick

end;pro TVcolorbar
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
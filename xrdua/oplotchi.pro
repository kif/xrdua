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

function ScaleYCHI,ptr,list,calcy=calcy

y=(*ptr).y*(*ptr).mult+(*ptr).off
if list.ypower ne 0 and list.ypower ne 1 then y^=list.ypower

if list.normalizey and not keyword_set(calcy) then begin

    ; Calc min/max from zoomed area
    bnormal=(*ptr).x[1,list.xtype] gt (*ptr).x[0,list.xtype]
    ind1=where(((*ptr).x[*,list.xtype]-list.xrange[0,list.xtype]) ge 0,ct1)
    ind2=where(((*ptr).x[*,list.xtype]-list.xrange[1,list.xtype]) ge 0,ct2)
    if bnormal then begin
        ind1=ind1[0]
        ind2=ind2[0]
        if ind1 eq -1 then ind1=n_elements(y)-1
        if ind2 eq -1 then ind2=n_elements(y)-1
    endif else begin
        ind1=ind1[ct1-1>0]
        ind2=ind2[ct2-1>0]
        if ind1 eq -1 then ind1=0
        if ind2 eq -1 then ind2=0
    endelse

    mi=min(y[0>(ind1<ind2):(ind1>ind2)<(n_elements(y)-1)],max=ma)

    ; Scale
    y-=mi
    y*=(list.yrange[1]-list.yrange[0])/(ma-mi)
    y+=list.yrange[0]
endif

return,y
end;function ScaleYCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro oplotchi_sum,list
if list.ndata le 1 then return

; ---- Make x extended ----
ptr=list.data
while PTR_VALID(ptr) do begin
    if ~((*ptr).hide) then begin
        x=(*ptr).x[*,list.xtype]
        xr=list.xrange[*,list.xtype]
        ind=where(x ge (xr[0]<xr[1]) and x le (xr[0]>xr[1]),ct)
        if ct ge 2 then begin
            x=x[ind]
            incc=x[1]-x[0]
            bc=min(x,max=ec)
            if n_elements(inc) ne 0 then begin
                b=b<bc
                e=e>ec
                inc=inc<incc
            endif else begin
                b=bc
                e=ec
                inc=incc
            endelse
        endif
    endif
    ptr=(*ptr).next
endwhile
if n_elements(inc) eq 0 then return
nEx=ceil((e-b)/inc)+1

; ---- Make y extented ----
yEx=fltarr(nEx)
ptr=list.data
while PTR_VALID(ptr) do begin
    if ~((*ptr).hide) then begin
        yEx >= EqualBin((*ptr).x[*,list.xtype], ScaleYCHI(ptr,list), b, inc, nEx, X2=xEx)
    endif
    ptr=(*ptr).next
endwhile

; ----Save sum----
window
plot,xEx,yEx,xs=1,ys=1
;temp=WriteChi(list.path+'sum_'+list.file,transpose([[xEx],[yEx]]),list.xtype,'oplotchisum',list.OutID,[0,2*!pi])

end;pro oplotchi_sum
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshXVal,list

first=1b
ptr=list.data
while ptr_valid(ptr) do begin
    (*ptr).x=CHI_xval((*ptr).x[*,(*ptr).xtype],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,(*ptr).xtype)

    if ~(*ptr).hide then begin
        xrange=transpose([[min((*ptr).x,DIMENSION=1)],[max((*ptr).x,DIMENSION=1)]])
        if first then begin
            first=0b
            list.xrange=xrange
        endif else list.xrange=[list.xrange[0,*]<xrange[0,*],list.xrange[1,*]>xrange[1,*]]
    endif

    ptr=(*ptr).next
endwhile

end;pro RefreshXVal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakePlotType1Img,list,error=error,xrange=xrange,yrange=yrange,xnew=xconv

; xnew is for surface plot
; xrange and yrange are for plottype1 (-> shifts so that the axis values are in the middle of the pixels)

error=1b
ptr=list.data
if not PTR_VALID(ptr) then return,0

inc=(list.xrange[1,list.xtype]-list.xrange[0,list.xtype])/(list.nbins-1)
xconv = list.xrange[0,list.xtype]+inc*lindgen(list.nbins)
xrange=list.xrange[*,list.xtype]+[-0.5*inc,0.5*inc]
yrange=[0,list.ndata]-0.5


bnnormal=(*ptr).x[0,list.xtype] gt (*ptr).x[1,list.xtype]
if bnnormal then begin
    x=reverse((*ptr).x[*,list.xtype])
    y=reverse(ScaleYCHI(ptr,list))
endif else begin
    x=(*ptr).x[*,list.xtype]
    y=ScaleYCHI(ptr,list)
endelse
img = minterpol( y, x, xconv ,/SPLINE)

while PTR_VALID((*ptr).next) do begin
    ptr=(*ptr).next

    bnnormal=(*ptr).x[0,list.xtype] gt (*ptr).x[1,list.xtype]
    if bnnormal then begin
        x=reverse((*ptr).x[*,list.xtype])
        y=reverse(ScaleYCHI(ptr,list))
    endif else begin
        x=(*ptr).x[*,list.xtype]
        y=ScaleYCHI(ptr,list)
    endelse
    yconv = minterpol( y, x, xconv ,/SPLINE)

    img=[[img],[yconv]]
endwhile

error=0b
return,img
end;function MakePlotType1Img
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotchis3d, in

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
    xx=(list.list2).xx
    yy=indgen((list.list2).yr[1]-(list.list2).yr[0]+1)+(list.list2).yr[0]

    ; ----Create the surface----
    oSurface = OBJ_NEW('IDLgrSurface', (list.list2).tiff,xx,yy , STYLE=3, $
               COLOR=[60,60,255], BOTTOM=[64,192,128], $
               XCOORD_CONV=NormCoord((list.list2).xr,list.PosX),$
               YCOORD_CONV=NormCoord((list.list2).yr,list.PosY),$
               ZCOORD_CONV=NormCoord((list.list2).zr,list.PosZ))
    oGroup->Add, oSurface
endelse

end;pro plotchis3d
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshPDD_OPLOTCHI,list

if list.PDFnTotal eq 0 then return

; ----Convert to all possible units----
(*list.PDFx)=CHI_xval((*list.PDFd),list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,0)

; ----Find position intensities----
ptr=list.data
if ptr_valid(ptr) then begin
    (*list.PDFi2)[*]=0
    (*list.PDFOffset)[*]=0
endif

case list.plottype of
1:    begin
;    ; Loop over all patterns
;    while PTR_VALID(ptr) do begin
;        temp=(*list.PDFx)[*,list.xtype]
;        mi=min((*ptr).x[*,list.xtype],max=ma)
;        ind=where((temp ge mi) and (temp le ma),ct)
;        if ct ne 0 then begin
;            ptry=ScaleYCHI(ptr,list)
;            hights = minterpol( ptry, (*ptr).x[*,list.xtype], temp[ind] ,/SPLINE)
;            (*list.PDFi2)[ind]>=hights
;        endif
;        ptr=(*ptr).next
;    endwhile
    (*list.PDFi2)[*]=list.ndata-0.5
    (*list.PDFOffset)[*]=-0.5
    endcase
0:    begin
    ; Loop over all patterns
    while PTR_VALID(ptr) do begin
        ; Get connected PDFs
        PDFj=strsplit((*ptr).PDFconnect,' ,;',/extract,count=n)
        
        if n ne 0 then begin
            mi=min((*ptr).x[*,list.xtype],max=ma)
            ptry=ScaleYCHI(ptr,list)
        endif
        
        ; Loop over connected PDFs
        for i=0l,n-1 do begin
            PDFi=long(PDFj[i])
            if PDFi ge list.nPDF then continue
            if PDFi eq 0 then i0=0l else i0=total((*list.PDFn)[0:PDFi-1],/int)
            i1=i0-1+(*list.PDFn)[PDFi]
            
            ; Interpolate reflections of this PDF in the x range
            temp=(*list.PDFx)[i0:i1,list.xtype]
            ind=where((temp ge mi) and (temp le ma),ct)
            if ct ne 0 then begin
                hights = minterpol( ptry, (*ptr).x[*,list.xtype], temp[ind] ,/SPLINE)
                (*list.PDFi2)[ind+i0]=hights
                
                (*list.PDFOffset)[PDFi]=min(ptry)
            endif
        endfor
    
        ptr=(*ptr).next
    endwhile
    endcase
endcase

end;pro RefreshPDD_OPLOTCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro plotchis,list,CGM=CGM,PS=PS

if keyword_set(CGM) then position=[0.20,0.15,0.90,0.60]

case list.plottype of
1:    begin
    ; ---make image----
    img=MakePlotType1Img(list,error=error,xrange=xrange,yrange=yrange)
    if error then return

    ; ----calc coord----
    ; plot corner coordinates (normalized)
    if keyword_set(CGM) then POS=position else $
    POS=[list.imgoffpr,list.imgoffpr, 1-list.imgoffpr/2., 1-list.imgoffpr/2.]

    ; new img x and y size to fit in these corners
    xsize=fix((POS[2]-POS[0])*list.xsize)
    ysize=fix((POS[3]-POS[1])*list.ysize)

    ; ----Show image----
    LoadctRev,list.ctl
    erase
    color=GetBackgroundContrastColor()
    
    imgbound=([(list.scaleimg-100)>0,list.scaleimg]<100)/100.
    m=list.yrange[0]+(list.yrange[1]-list.yrange[0])*imgbound

    tvscl,congrid(m[0]>img<m[1],xsize,ysize),POS[0],POS[1],xsize=POS[2]-POS[0],ysize=POS[3]-POS[1],/normal

    if list.ytitle eq '' then ytitle='' else ytitle='scannr'
    plot,xrange,yrange,xstyle=1,ystyle=1,xtitle=list.xtitle[list.xtype],ytitle=ytitle,$
    title=list.title,xrange=xrange,yrange=yrange,/nodata,/noerase,$
    POSITION=POS,color=color,charsize=list.char

    ROIy=!Y.CRange
    endcase
0:    begin
    LoadctRev,-39,/silent
    erase
    color=GetBackgroundContrastColor()

    ; ----Data----
    i=0
    xpos=list.xpos*!d.x_size
    ypos=list.ypos*!d.y_size
    incpos=list.incpos*list.char*!d.y_size
    GetCRange,ROIy=ROIy

    yrange=list.yrange
;    yrange+=(yrange[1]-yrange[0])*[0,0.08]
    plot,list.xrange[*,list.xtype],yrange,xstyle=1,ystyle=1,ylog=list.ylog,xlog=list.xlog,xtitle=list.xtitle[list.xtype],ytitle=list.ytitle+' '+list.ypowerstr,$
    title=list.title,xrange=list.xrange[*,list.xtype],yrange=yrange,/nodata,position=position,charsize=list.char,color=color
    if not ptr_valid(list.data) then return
    ptr=list.data
    while PTR_VALID(ptr) do begin
        if ~((*ptr).hide) then begin
            oplot,(*ptr).x[*,list.xtype],ScaleYCHI(ptr,list),col=(*ptr).col
            ;x=(*ptr).x[*,list.xtype]
            ;y=ScaleYCHI(ptr,list)
            ;test,5*x[20:57],y[20:57],list.outid
            ;test,15*x[50:70],y[50:70],list.outid
            ;return
            if list.slegend then begin
                xyouts,xpos,ypos,'_____ '+(*ptr).label,/device,col=(*ptr).col,charsize=list.char
                ypos+=incpos
            endif
        endif
        i++
        ptr=(*ptr).next
    endwhile
    endcase
endcase

    ;----Display PDF data----
    if (list.PDFnTotal ne 0) then begin
        LoadctRev,ntableadd(1),/silent
        off=0
        mi=list.xrange[0,list.xtype]
        ma=list.xrange[1,list.xtype]
        labellinex=!d.y_size*list.char*0.015
        labelliney=!d.y_size*list.char*0.010
        xpos=list.xpos+0.2
        ypos=list.ypos
        incpos=list.incpos*list.char
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
                PDFOffset=(*list.PDFOffset)[i]
                if ((*list.PDFni)[i] ne 0) then begin
                    scl=(*list.PDFScale)[i]
                    PDFi=(*list.PDFi)[ind]
                    PDFi/=max(PDFi)
                    if (scl ge 0) and (scl le 1) then begin
                        PDFi*=(ROIy[1]-PDFOffset)*scl
                        PDFi+=PDFOffset
                    endif else begin
                        PDFi2=(*list.PDFi2)[ind]-PDFOffset  ; Pattern values at Bragg lines - offset
                        m=max(PDFi,indm)    ; most intens line
                        m=min(PDFi2[indm])    ; m = Pattern value of the most intense line
                        
                        PDFi*=m
                        PDFi+=PDFOffset
                    endelse
                endif else PDFi=(*list.PDFi2)[ind]
                
                ; Labels
                hkl=fix((*list.PDFhkl)[*,ind])
                mihkl=min(hkl)
                n=strlen(stringr(max(hkl)))>(strlen(stringr(mihkl)))
                label=string(hkl,format='(3(I'+stringr(n)+'))')

                for j=0l,ct-1 do begin
                    if (code[j] and 1) eq 1 then begin
                        PlotS,[x[j],x[j]],[PDFOffset,PDFi[j]],color=col,NOCLIP = 0,/data
                        if list.bPDFpsym then PlotS,x[j],PDFi[j]+psymoffset,psym=psym,NOCLIP = 0, color=col,/data,symsize=0.7
                    endif
                    if (code[j] and 2) eq 2 then begin
                        cout=CONVERT_COORD(x[j],PDFi[j],/DATA,/TO_DEVICE)
                        PlotS,[cout[0],cout[0]+labellinex],[cout[1],cout[1]+labelliney],color=col,/device,NOCLIP = 0
                        xyouts,cout[0]+labellinex,cout[1]+labelliney,label[j],color=col,/device,charsize=list.char,NOCLIP = 0
                    endif
                endfor
                
                ;if list.showlegend then begin
                    ;if list.blegendshort then begin
                        p=strpos((*list.PDFname)[i],';')
                        if p ne -1 then begin
                            p2=strpos((*list.PDFname)[i],';',p+1)
                            if p2 ne -1 then p = p2
                            tmp=strmid((*list.PDFname)[i],0,p)
                        endif else tmp=(*list.PDFname)[i]
                    ;endif else tmp=(*list.PDFname)[i]
                    if list.bPDFpsym then begin
                        xyouts,xpos,ypos,'!3'+tmp,/normal,col=col,charsize=list.char
                        plots,xpos*0.95,ypos*1.01,/normal,col=col,psym=psym,symsize=0.7
                    endif else $
                        xyouts,xpos,ypos,'!3_____ '+tmp,/normal,col=col,charsize=list.char
                    
                    ypos+=incpos
                ;endif
            endif

            off+=(*list.PDFn)[i]
        endfor
    endif

end;pro plotchis
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshDisplayCHIS,list
wset, list.drawindex
plotchis,list
GetSysVar,list.sysvar
wset, list.pixIndex
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.drawindex]
end;pro RefreshDisplayCHIS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshDataOPLOTCHI,ev
widget_control,ev.top,get_uvalue=list

; ----PDF markers----
RefreshPDD_OPLOTCHI,list

widget_control,ev.top,set_uvalue=list
RefreshDisplayCHIS,list
end;pro RefreshDataOPLOTCHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditOplotchiParam_cleanup,ID

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
end;pro EditOplotchiParam_cleanup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditOplotchiParam_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_uvalue=uval
widget_control,ev.top,get_uvalue=ID
widget_control,ID,get_uvalue=list
widget_control,ev.id,get_value=val
dval=0D
case uval of
    'OK' :  begin
         WIDGET_CONTROL, ev.top, /DESTROY
         return
         endcase
    'xlambda':begin
         ID1=widget_info(ev.top,FIND_BY_UNAME='xen')
         reads,val,dval
         list.lambda=dval
         widget_control,ID1,set_value=stringr(list.cconv/dval)
         RefreshXVal,list
         endcase
    'xen': begin
         ID1=widget_info(ev.top,FIND_BY_UNAME='xlambda')
         reads,val,dval
         list.lambda=list.cconv/dval
         widget_control,ID1,set_value=stringr(list.lambda)
         RefreshXVal,list
         endcase
    'scrx':begin
         reads,val,dval
         list.scrx=dval*10^(-6.) ; from um to m
         RefreshXVal,list
         endcase
  'scry':begin
         reads,val,dval
         list.scry=dval*10^(-6.) ; from um to m
         RefreshXVal,list
         endcase
  'dist':begin
         reads,val,dval
         list.dist=dval*10^(-2.) ; from cm to m
         RefreshXVal,list
         endcase
    'char':   begin
         reads,val,dval
         list.char=dval
         endcase
    'inc':   begin
         reads,val,dval
         list.incpos=dval
         endcase
    'yrange0':begin
        reads,val,dval
        list.yrange[0]=dval
        endcase
    'yrange1':begin
        reads,val,dval
        list.yrange[1]=dval
        endcase
    'power':begin
            list.ypower=float(val)
            if list.ypower eq 1 or list.ypower eq 0 then begin
                list.ypowerstr=''
            endif else begin
                list.ypowerstr='^'+val
                list.ylog=0
                id2=widget_info(ID,find_by_uname='Y')
                widget_control,id2,set_value='Y: linear'
            endelse
            endcase
    'title':list.title=val
    'xtitle':list.xtitle[list.xtype]=val
    'ytitle':list.ytitle=val
    'Show Cross':list.scross=ev.select
    'Normalize Y':list.normalizey=ev.select
    'Show Legend':list.slegend=ev.select
    'boarder':list.imgoffpr=0>float(val)/100.<1
     'ctl':list.ctl=fix(val)
     'PDFpsym': list.bPDFpsym=ev.select
endcase
printw,list.OutID,'Edit '+uval
widget_control,ID,set_uvalue=list
RefreshDataOPLOTCHI,{top:ID}
end;pro EditOplotchiParam_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditOplotchiParam,ev
widget_control,ev.top,get_uvalue=list
widget_control,ev.top,sensitive=0
base=widget_base(title='oplotchi parameters',uvalue=ev.top,/row)
base1=widget_base(base,/column)
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
label=widget_label(baset,value='Charsize  :        ',xsize=xs)
text=widget_text(baset,/editable,value= $
stringr(list.char),uvalue='char',uname='char')

baset=widget_base(base1,/row)
label=widget_label(baset,value='Legend spacing (*char)  :        ',xsize=xs)
text=widget_text(baset,/editable,value= $
stringr(list.incpos),uvalue='inc',uname='inc')

baset=widget_base(base1,/row)
label=widget_label(baset,value='Title  :        ',xsize=xs)
text=widget_text(baset,/editable,value= $
list.title,uvalue='title',uname='title')

baset=widget_base(base1,/row)
label=widget_label(baset,value='XTitle  :        ',xsize=xs)
text=widget_text(baset,/editable,value= $
list.xtitle[list.xtype],uvalue='xtitle',uname='xtitle')

baset=widget_base(base1,/row)
label=widget_label(baset,value='YTitle  :        ',xsize=xs)
text=widget_text(baset,/editable,value= $
list.ytitle,uvalue='ytitle',uname='ytitle')

baset=widget_base(base1,/row)
label=widget_label(baset,value='Y-range  :        ',xsize=xs)
text=widget_text(baset,/editable,value= $
string(list.yrange[0]),uvalue='yrange0',uname='yrange0')
text=widget_text(baset,/editable,value= $
string(list.yrange[1]),uvalue='yrange1',uname='yrange1')

baset=widget_base(base1,/row)
label=widget_label(baset,value='Y-axis power  :    ',xsize=xs)
text=widget_text(baset,/editable,value= strmid(list.ypowerstr,1),uvalue='power',uname='power',scr_xsize=xs)

baset=widget_base(base1,/column,/nonexclusive)
button=widget_button(baset,value='Show Cross',uvalue='Show Cross')
widget_control,button,set_button=list.scross
button=widget_button(baset,value='Normalize Y',uvalue='Normalize Y')
widget_control,button,set_button=list.normalizey
button=widget_button(baset,value='Show Legend',uvalue='Show Legend')
widget_control,button,set_button=list.slegend
button=widget_button(baset,value='PDF symbols',uvalue='PDFpsym')
widget_control,button,set_button=list.bPDFpsym

baset=widget_base(base1,/row)
label=widget_label(baset,value='Boarder offset% :     ',xsize=xs)
text=widget_text(baset,value=stringr(list.imgoffpr*100),uvalue='boarder',/editable)

baset=widget_base(base1,/row)
label=widget_label(baset,value='Color table     :          ',xsize=xs)
text=widget_text(baset,/editable,value= $
stringr(list.ctl),uvalue='ctl',uname='ctl')

button=widget_button(base,value='OK',uvalue='OK')
WIDGET_CONTROL, base, /REALIZE
Xmanager,'EditOplotchiParam',base, event_handler='EditOplotchiParam_event',$
cleanup='EditOplotchiParam_cleanup',GROUP_LEADER=ev.top
end;pro EditOplotchiParam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ConstructDataList,ev

widget_control,ev.top,get_uValue=list
listtype=list.listtype
if list.listtype eq 30 then begin
    list.listtype=31
    widget_control,list.draw,get_value=drawindex
    list.drawindex=drawindex

;    widget_control,ev.top,set_uValue=list
;    ID=widget_info(ev.top,FIND_BY_UNAME='Edit')
;    widget_control,ID,sensitive=1
;    ID=widget_info(ev.top,FIND_BY_UNAME='Perform')
;    widget_control,ID,sensitive=1
endif

pathlist=Select_Files(ev.top,list,sortmethod=list.sortmethod,$
    separator=list.sortseparator,filter='*'+FileFormat(3))
if pathlist[0] eq '' then return
n=n_elements(pathlist)
Result = DIALOG_MESSAGE(stringr(n)+' files to add: proceed?',/question)
if result eq 'No' then return

; ----Insert pathlist----
ptr=list.data
if ptr_valid(ptr) then $
    while PTR_VALID((*ptr).next) do ptr=(*ptr).next

for i=0l,n-1 do begin
    res=CutPath(pathlist[i],path=path,file=label,ext=ext)
    list.Format=ext
    file=label+ext
    error=ReadScan(path,file,list.Format,list.OutID,xval=xval,spe=spe,type=xtype,asciitemplate=list.asciitemplate)

    if (error eq 0) then begin
        printw,list.OutID,file+' not loaded.'
    endif else begin
        first=not ptr_valid(list.data)

        xvalt=CHI_xval(xval,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,xtype)
        xrange=transpose([[min(xvalt,DIMENSION=1)],[max(xvalt,DIMENSION=1)]])
        yrange=[min(spe),max(spe)]

        s=size(spe)
        nbins=s[1]
        nadd=s[(s[0] eq 1)?0:2]

        if first then begin
            list.xrange=xrange
            list.yrange=yrange
            list.nbins=nbins
        endif else begin
            list.xrange=[list.xrange[0,*]<xrange[0,*],list.xrange[1,*]>xrange[1,*]]
            list.yrange=[list.yrange[0]<yrange[0],list.yrange[1]>yrange[1]]
            list.nbins>=nbins
        endelse

        col=0
        maxcol=n_elements(list.colors)-1
        if not first then begin
            col=(where(list.colors eq (*ptr).col))[0]
            if col ne maxcol then col++ else col=0
        endif
        
        Ix='I'+stringr(strlen(stringr(nadd-1)))
        format='(A,"[",'+Ix+',"]")'

        label2=label
        for j=0l,nadd-1 do begin
            label2=string(label,j,format)

            GenericListInsertAfter,ptr,PTR_NEW({x:xvalt,$
                y:spe[*,j],xtype:xtype,off:0.,mult:1.,col:list.colors[col],$
                label:label2,path:path,file:file,hide:0b,PDFconnect:'',$
                next:PTR_NEW(),prev:PTR_NEW()})
            list.ndata++
            if first then list.data=ptr
            first=0b
            if col ne maxcol then col++ else col=0
        endfor

        printw,list.OutID,file+' added.'

    endelse
endfor

; Last loaded file
list.path=path
list.file=file
widget_control,ev.top,TLB_SET_TITLE=list.bartitle+list.file

; ----Load Mask and Refresh----
if listtype eq 30 then begin
    if list.main ne 0 then error=0b else LoadMask,list,0,list.OutID,error=error,/oD,flagi=['$IntCor2Dinfo:','$FIT_ROI:','$BACKoD:','$FIL:','$LABELS:','$MASK1:','$FITMODEL:']
    if error then begin
           path=file_search(list.mskpath+'*.msk',count=count)
           if count ne 0 then begin
             path=path[0]
             error=CutPath(path,path=mskpath,file=file,ext=ext)
             list.mskpath=mskpath
             list.mskfile=file+ext
             LoadMask,list,0,list.OutID,/oD,flagi=['$IntCor2Dinfo:','$FIT_ROI:','$BACKoD:','$FIL:','$LABELS:','$MASK1:','$FITMODEL:']

             RefreshXVal,list
             widget_control,ev.top,set_uvalue=list
             RefreshDataOPLOTCHI,ev
           endif else begin
                widget_control,ev.top,set_uvalue=list
             RefreshDisplayCHIS,list
           endelse
    endif else begin
        RefreshXVal,list
        widget_control,ev.top,set_uvalue=list
        RefreshDisplayCHIS,list
    endelse
endif else begin
    widget_control,ev.top,set_uvalue=list
    RefreshDisplayCHIS,list
endelse
end;pro ConstructDataList
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditChi_cleanup,ID

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
end;pro EditChi_cleanup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditChiPart,wrapbase,list
result=ReadINI('$EditMask:')

basegroups=widget_base(wrapbase,/column,uname='basegroups')
if list.ndata eq 0 then return

; Make different rows, depending on how many (nr) groups you want on a row
nr=result.(0)
nbase=ceil(float(list.ndata)/nr)
base6=lonarr(nbase)
for i=0l,nbase-1 do base6[i]=widget_base(basegroups,/row)

; Loop over all groups and place them in the right row
j=0
ptr=list.data
for i=0l,list.ndata-1 do begin
    si=stringr(i)

    ; Base for the group
    bgroup=widget_base(base6[j],/column,/frame)

    ; Text boxes
    bgroupt=widget_base(bgroup,/column)
        text=widget_text(bgroupt,value=(*ptr).label,/editable,uvalue=[i,0])
        text=widget_text(bgroupt,value=stringr((*ptr).off),/editable,uvalue=[i,1],uname='T'+si)
        text=widget_text(bgroupt,value=stringr((*ptr).mult),/editable,uvalue=[i,2])
        text=widget_text(bgroupt,value=stringr((*ptr).col),/editable,uvalue=[i,3])

    ; Buttons
    bgroupne=widget_base(bgroup,/column,/NONEXCLUSIVE)
        button=widget_button(bgroupne,value='Hide',uvalue=i,uname='H'+si)
        widget_control,button,set_button=(*ptr).hide
        button=widget_button(bgroupne,value='Delete',uvalue=[i,bgroupt],uname='D'+si)

    ; PDF Droplist
    pdfnames='Add PDF...'
    if list.nPDF ne 0 then pdfnames=[pdfnames,*list.PDFname]
    bgroupp=widget_base(bgroup,/column)
        drop=widget_droplist(bgroupp,value=pdfnames,uvalue=i)
        text=widget_text(bgroupp,value=(*ptr).PDFconnect,/editable,uvalue=[i,4],uname='P'+si)
        
    if ((i+1) mod nr) eq 0 then j=j+1
    ptr=(*ptr).next
endfor

end;pro EditChiPart
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro oplotchi_hide,list,show=show,prev=prev
if list.plottype eq 1 then return

bhide=~keyword_set(show)
bprev=keyword_set(prev)

ptr=list.data
ptr2=ptr_new()
bool=0b
blast=0b
bfirst=1b
for i=0l,list.ndata-1 do begin
    bool or= (*ptr).hide
    if ~(*ptr).hide and bfirst then begin
        bfirst=0b
        ptr2=bprev?(*ptr).prev:(*ptr).next
        if ~ptr_valid(ptr2) then $
            if bprev then blast=1b else ptr2=list.data
    endif
    (*ptr).hide=bhide
    if blast then ptr2=ptr
    ptr=(*ptr).next
endfor

if bhide and ptr_valid(ptr2) then (*ptr2).hide=0b
end;pro oplotchi_hide
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditChi_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
1:    begin; button
    widget_control,ev.id,get_value=val
    case val of
    'Exit':    begin
            widget_control,ev.top,get_uvalue=list2
            ind=where(list2.markdel,ct)
            if ct ne 0 then begin
                widget_control,list2.top,get_uvalue=list,/no_copy
                ptr=list.data
                GenericListDeleteI,ptr,ind
                list.data=ptr
                list.ndata-=ct
                RefreshPDD_OPLOTCHI,list
                RefreshDisplayCHIS,list
                widget_control,list2.top,set_uvalue=list,/no_copy
            endif
            widget_control,ev.top,/destroy
            endcase
    'Offsets++':begin
            widget_control,ev.top,get_uvalue=list2
            widget_control,list2.top,get_uvalue=list
            ptr=list.data
            d=(list.yrange[1]-list.yrange[0])/100.
            for i=0l,list.ndata-1 do begin
                (*ptr).off+=i*d
                ID=widget_info(ev.top,FIND_BY_UNAME='T'+stringr(i))
                widget_control,ID,set_value=stringr((*ptr).off)
                ptr=(*ptr).next
            endfor
            RefreshDisplayCHIS,list
            endcase
    'Offsets--':begin
            widget_control,ev.top,get_uvalue=list2
            widget_control,list2.top,get_uvalue=list
            ptr=list.data
            d=list.yrange[1]-list.yrange[0]
            for i=0l,list.ndata-1 do begin
                (*ptr).off-=i*d
                ID=widget_info(ev.top,FIND_BY_UNAME='T'+stringr(i))
                widget_control,ID,set_value=stringr((*ptr).off)
                ptr=(*ptr).next
            endfor
            RefreshDisplayCHIS,list
            endcase
    'Hide': begin
            widget_control,ev.top,get_uvalue=list2
            widget_control,list2.top,get_uvalue=list
            widget_control,ev.id,get_uvalue=uval
            ptr=list.data
            for j=0l,uval-1 do ptr=(*ptr).next
            (*ptr).hide=ev.select
            RefreshDisplayCHIS,list
            endcase

    'Delete':begin
            widget_control,ev.id,get_uvalue=uval
            widget_control,uval[1],sensitive= ~ev.select
            widget_control,ev.top,get_uvalue=list2,/no_copy
            list2.markdel[uval[0]]=ev.select
            widget_control,ev.top,set_uvalue=list2,/no_copy
            endcase

    'Hide All':begin
            widget_control,ev.top,get_uvalue=list2
            widget_control,list2.top,get_uvalue=list
            ptr=list.data
            for i=0l,list.ndata-1 do begin
                ID=widget_info(ev.top,FIND_BY_UNAME='H'+stringr(i))
                widget_control,ID,set_button=ev.select
                (*ptr).hide=ev.select
                ptr=(*ptr).next
            endfor
            RefreshDisplayCHIS,list
            endcase

    'Delete All':begin
            widget_control,ev.top,get_uvalue=list2,/no_copy
            widget_control,list2.top,get_uvalue=list
            ptr=list.data
            for i=0l,list.ndata-1 do begin
                ID=widget_info(ev.top,FIND_BY_UNAME='D'+stringr(i))
                widget_control,ID,set_button=ev.select
                widget_control,ID,get_uvalue=uval
                widget_control,uval[1],sensitive= ~ev.select
                list2.markdel[uval[0]]=ev.select
                ptr=(*ptr).next
            endfor
            widget_control,ev.top,set_uvalue=list2,/no_copy
            RefreshDisplayCHIS,list
            endcase

    endcase
    endcase
3:    begin; text
    widget_control,ev.top,get_uvalue=list2
    widget_control,list2.top,get_uvalue=list
    widget_control,ev.id,get_uvalue=uval
    widget_control,ev.id,get_value=valstr

    i=uval[0]
    ptr=list.data
    for j=0l,i-1 do ptr=(*ptr).next
    case uval[1] of
    0:    begin
        (*ptr).label=valstr
        endcase
    1:    begin
        val=0.
        reads,valstr,val
        (*ptr).off=val
        ;ChangeXRange,0,0,list
        widget_control,list2.top,set_uvalue=list
        endcase
    2:    begin
        val=0.
        reads,valstr,val
        (*ptr).mult=val
        ;ChangeXRange,0,0,list
        widget_control,list2.top,set_uvalue=list
        endcase
    3:    begin
        val=0b
        reads,valstr,val
        (*ptr).col=val
        endcase
    4:    begin
        (*ptr).PDFconnect=valstr
        endcase
    endcase
    RefreshPDD_OPLOTCHI,list
    RefreshDisplayCHIS,list
    endcase
8:    begin
    if ev.index eq 0 then return
    PDFi=ev.index-1
    widget_control,ev.id,get_uvalue=i
    
    widget_control,ev.top,get_uvalue=list2
    widget_control,list2.top,get_uvalue=list
    
    ; Remove PDFi from all patterns j
    ptr=list.data
    j=0l
    while ptr_valid(ptr) do begin
        PDFj=strsplit((*ptr).PDFconnect,' ,;',/extract,count=ct)
        if ct ne 0 then begin
            PDFj=long(PDFj)
            ind=where(PDFj ne PDFi,ct)
            if ct eq 0 then (*ptr).PDFconnect='' else (*ptr).PDFconnect=strjoin(stringr(PDFj[ind]),';')
            
            ID=widget_info(ev.top,FIND_BY_UNAME='P'+stringr(j))
            widget_control,ID,set_value=(*ptr).PDFconnect
        endif
        
        if j eq i then ptrkeep=ptr
        
        ptr=(*ptr).next
        j++
    endwhile 
    
    ; Add PDFi to pattern i
    tmp=strsplit((*ptrkeep).PDFconnect,' ,;',/extract,count=ct)
    if ct eq 0 then (*ptrkeep).PDFconnect=stringr(PDFi) $
    else (*ptrkeep).PDFconnect=strjoin([tmp,stringr(PDFi)],';')

    
    ID=widget_info(ev.top,FIND_BY_UNAME='P'+stringr(i))
    widget_control,ID,set_value=(*ptrkeep).PDFconnect
    
    RefreshPDD_OPLOTCHI,list
    RefreshDisplayCHIS,list
    endcase
endcase
end;pro EditChi_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditChi,ev

widget_control,ev.top,get_uvalue=list
if not PTR_VALID(list.data) then return
WIDGET_CONTROL, ev.top, sensitive=0

result=ReadINI('$EditMask:')

base=widget_base(/column,title='Edit Chi',xoffset=100,yoffset=300,$
uvalue={top:ev.top,markdel:bytarr(list.ndata)},/scroll,y_scroll_size=result.(1),x_scroll_size=result.(1))

basebuttons=widget_base(base,/row)
    button=widget_button(basebuttons,value='Exit',uname='Exit')
    button=widget_button(basebuttons,value='Offsets++')
    button=widget_button(basebuttons,value='Offsets--')
    basebuttons2=widget_base(basebuttons,/row,/nonexclusive)
        button=widget_button(basebuttons2,value='Hide All')
        button=widget_button(basebuttons2,value='Delete All')

wrapbase=widget_base(base,uname='wrapbase')
EditChiPart,wrapbase,list

WIDGET_CONTROL, base, /REALIZE
Xmanager,'EditChi',base, event_handler='EditChi_event',$
    cleanup='EditChi_cleanup',GROUP_LEADER=ev.top

end;pro EditChi
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ModeOChi,ev,val

widget_control,ev.top,get_uvalue=list

modes=['Edit xrange','Legend Position','PlotType1','Restore']
ind=where(modes eq val,ct)
if ct eq 0 then begin
    widget_control,list.draw,draw_button=0,draw_motion=0
    return
endif

ID=widget_info(ev.top,FIND_BY_UNAME='basedyn')

if ind[0] le 1 then begin
    if list.plottype eq 1 then return
    widget_control,list.draw,draw_button=0,draw_motion=0
    case val of
        'Edit xrange': begin
                widget_control,list.draw,Event_Pro='CHANGE_XRANGE_EVENT',draw_button=1
                widget_control,ID,set_uvalue=val
                 endcase
        'Legend Position':begin
                widget_control,list.draw,Event_Pro='CHANGE_LEGEND_EVENT',draw_button=1,draw_motion=1
                widget_control,ID,set_uvalue=''
                endcase
    endcase
endif else begin
    widget_control,list.draw,draw_button=0,draw_motion=0
    tmp=widget_info(ev.top,FIND_BY_UNAME='Mbase')
    if tmp ne 0 then widget_control,tmp,/destroy
    base=widget_base(ID,/column,uname='Mbase')

    case val of
    'PlotType1':begin
               text=widget_text(base,ysize=7,uname='plotstyle1info')

                baset=widget_base(base,/row)
                    label=widget_label(baset,value='Image Scale: ',xsize=xs)
                    slider=widget_slider(baset,minimum=0,maximum=200,uname='scale',uvalue='scale',value=list.scaleimg)

                widget_control,list.draw,Event_Pro='PLOTTYPE1_EVENT',draw_button=1,draw_motion=1
                endcase
    'Restore':begin
                baset=widget_base(base,/row)
                    b=widget_button(baset,value='Show All',xsize=xs)
                    b=widget_button(baset,value='Show Next',xsize=xs)
                    b=widget_button(baset,value='Show Prev',xsize=xs)

                widget_control,ID,get_uvalue=val2
                ModeOChi,ev,val2
                endcase
    endcase
endelse

end;pro ModeOChi
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CHANGE_XRANGE_EVENT,ev

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

Widget_Control, ev.id, Event_Pro='CHANGE_XRANGE_DRAW', $
    draw_motion=1
Widget_Control, ev.top, Set_UValue=list
end;pro CHANGE_XRANGE_EVENT
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ChangeXRange,x1,x2,list

ptr=list.data
if not PTR_VALID(ptr) then return


if x1 eq x2 then begin
    first=1b
    while ptr_valid(ptr) do begin

        if ~(*ptr).hide then begin
            xrange=transpose([[min((*ptr).x,DIMENSION=1)],[max((*ptr).x,DIMENSION=1)]])

            if list.ylog eq 1 then begin
                tmp=ScaleYCHI(ptr,list,/calcy)
                ind=where(tmp ne 0,ct)
                if ct ne 0 then mi=min(tmp[ind],max=ma) else mi=min(tmp,max=ma)
            endif else mi=min(ScaleYCHI(ptr,list,/calcy),max=ma)

            if first then begin
                first=0b
                list.xrange=xrange
                list.yrange=[mi,ma]
            endif else begin
                list.xrange=[list.xrange[0,*]<xrange[0,*],list.xrange[1,*]>xrange[1,*]]
                list.yrange=[mi<list.yrange[0],ma>list.yrange[1]]
            endelse
        endif

        ptr=(*ptr).next
    endwhile
endif else begin
    xr=[x1<x2,x1>x2]

    while PTR_VALID(ptr) do begin
        if ~(*ptr).hide then begin
            bnormal=(*ptr).x[1,list.xtype] gt (*ptr).x[0,list.xtype]
            ind1=where(((*ptr).x[*,list.xtype]-xr[0]) ge 0,ct1)
            ind2=where(((*ptr).x[*,list.xtype]-xr[1]) ge 0,ct2)
            if bnormal then begin
                ind1=ind1[0]
                ind2=ind2[0]
                if ind1 eq -1 then ind1=n_elements((*ptr).x[*,list.xtype])-1
                if ind2 eq -1 then ind2=n_elements((*ptr).x[*,list.xtype])-1
            endif else begin
                ind1=ind1[ct1-1>0]
                ind2=ind2[ct2-1>0]
                if ind1 eq -1 then ind1=0
                if ind2 eq -1 then ind2=0
            endelse
            ma=n_elements((*ptr).x[*,list.xtype])-1
            ind1=0>ind1<ma
            ind2=0>ind2<ma
            if ind1 ne ind2 then begin
                y=ScaleYCHI(ptr,list,/calcy)
                y=y[ind1<ind2:ind1>ind2]
                
                if list.ylog eq 1 then begin
                    ind=where(y ne 0,ct)
                    if ct ne 0 then mi=min(y[ind],max=ma) else mi=min(y,max=ma)
                endif else mi=min(y,max=ma)
            
                if n_elements(yr) eq 0 then yr=[mi,ma] else yr=[mi<yr[0],ma>yr[1]]
            endif
        endif
        ptr=(*ptr).next
    endwhile

    list.xrange=CHI_xval(xr,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtype)
    if list.xrange[0] gt list.xrange[1] then list.xrange[[0,1]]=list.xrange[[1,0]]
    if n_elements(yr) ne 0 then list.yrange=yr
endelse
end;pro ChangeXRange
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CHANGE_XRANGE_DRAW,ev

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
       Event_Pro='CHANGE_XRANGE_EVENT'

    ChangeXRange,x1,x2,list
    widget_control,ev.top,set_uvalue=list
    RefreshDisplayCHIS,list
    return
endif

y2 = ROIy[0] > coords[1,1] < ROIy[1]
WSet, list.drawindex
Arrow, x1, y2, x2, y2, /Data, /Solid, HSize=12
PlotS, [x2, x2], ROIy
end;pro CHANGE_XRANGE_DRAW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PLOTTYPE1_EVENT,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

ID=widget_info(ev.top,FIND_BY_UNAME='plotstyle1info')
if not widget_info(ID,/VALID_ID) then return

widget_control,ev.top,get_uvalue=list

; plot cross
WSet, list.drawindex
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.pixindex]

SetSysVar,list.sysvar
xy=CONVERT_COORD(ev.x,ev.y,/device,/to_data)
if list.scross then begin
    PlotS, [xy[0],xy[0]],!y.crange,/data,noclip=0
    PlotS, !x.crange,[xy[1],xy[1]],/data,noclip=0
endif

; point to correct file
y=0>fix(xy[1]+0.5)<list.ndata

ptr=list.data
if not ptr_valid(ptr) then return
for i=0l,y-1 do if ptr_valid((*ptr).next) then ptr=(*ptr).next

; recalc x and y values
inc=(list.xrange[1,list.xtype]-list.xrange[0,list.xtype])/(list.nbins-1)
xconv = list.xrange[0,list.xtype]+inc*lindgen(list.nbins)
yconv= ScaleYCHI(ptr,list)
yconv = minterpol( yconv, (*ptr).x[*,list.xtype], xconv ,/SPLINE)

; find nearest x-index
tmp=abs(xconv-xy[0])
tmp=min(tmp,x)
x=x[0]

; x,y
y=string(yconv[x])
x=CHI_xval(xconv[x],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtype)
x=string(x)+list.shortxunits

; display info
widget_control,ID,set_value=[(*ptr).file,(*ptr).label,'x = '+x[0],'    '+x[1],'    '+x[2],'    '+x[3], 'I = '+y ]

end;pro PLOTTYPE1_EVENT
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CHANGE_LEGEND_EVENT,ev

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

IF thisEvent EQ 'UP' THEN BEGIN
    list.xpos=ev.x/float(!d.x_size)
    list.ypos=ev.y/float(!d.y_size)
    widget_control,ev.top,set_uvalue=list
    RefreshDisplayCHIS,list
    ModeOChi,ev,'dummy'
    return
ENDIF

WSet, list.drawindex
ptr=list.data
xpos=ev.x
ypos=ev.y
incpos=list.incpos*list.char*!d.y_size
i=0
while PTR_VALID(ptr) do begin
    xyouts,xpos,i*incpos+ypos,'_____ '+(*ptr).label,/device,col=(*ptr).col,charsize=list.char
    ptr=(*ptr).next
    i=i+1
endwhile

end;pro CHANGE_LEGEND_EVENT
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro oplotchi_cleanup,ID

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
    result=UpdateINI('$XRD pddpath:',{path:list.pddpath,file:list.pddfile})
    result=UpdateINI('$XRD mskpath:',{path:list.mskpath,file:list.mskfile})

    if widget_info(list.main,/VALID_ID) then begin
        widget_control,list.main,get_uvalue=list2
        list2.chipath=list.path
        list2.chifile=list.file
        widget_control,list.main,set_uvalue=list2
    endif
    
    heap_free,list
endif

error=UnSubscribeToOutID()

;help,/heap
print,'=========================================STOP 1D Profile Compare'
end;pro oplotchi_cleanup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro oplotchi_savesession,ev
WIDGET_CONTROL, ev.top, get_uvalue=list
datpath=selfile(list.path,'oplotchi.dat','*.dat','Save session...')
if datpath eq '' then return
SAVE,list,FILENAME=datpath
end;pro oplotchi_savesession
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro oplotchi_restoresession,ev

widget_control,ev.top,get_uvalue=list2

datpath=DIALOG_PICKFILE(path=list2.path,file='oplotchi.dat',$
         filter='*.dat',title='Load session....')
if datpath eq '' then return

; Restore list
RESTORE,FILENAME=datpath

; Copy fields except widget IDs
RestoreStruc,list2,list,exclude=['draw','OutID','main']

; Set display as restored
widget_control,list2.draw,get_value=drawindex,xsize=list2.xsize,ysize=list2.ysize
list2.drawindex=drawindex
wdelete,list2.pixindex
Window, /Free, xsize=list2.xsize, ysize=list2.ysize, /Pixmap
list2.pixindex = !D.Window

widget_control,ev.top,set_uvalue=list2

RefreshPDD_OPLOTCHI,list2
RefreshDisplayCHIS,list2
end;pro oplotchi_restoresession
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro oplotchi_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
8:    begin
    widget_control,ev.id,get_uvalue=uval
    case uval of
        'plotstyle':begin
            WIDGET_CONTROL, ev.top, get_uvalue=list
            list.plottype=ev.index
            widget_control,ev.top,set_uvalue=list
            ModeOChi,ev,((list.plottype eq 1)?'PlotType1':'Restore')
            RefreshDataOPLOTCHI,ev
            endcase
    endcase
    endcase
2:    begin
    widget_control,ev.id,get_uvalue=uval
    widget_control,ev.id,get_value=val
    case uval of
    'scale':begin
            WIDGET_CONTROL, ev.top, get_uvalue=list
            list.scaleimg=val
            WIDGET_CONTROL, ev.top, set_uvalue=list
            RefreshDataOPLOTCHI,ev
            endcase
    endcase
    endcase
1:    begin; button
    widget_control,ev.id,get_value=val
    case val of
    'Load chi': ConstructDataList,ev
    'Exit':    begin
            widget_control,ev.top,/destroy
            endcase
    'Save Image':begin
            WIDGET_CONTROL, ev.top, get_uvalue=list
            if n_elements(list) EQ 0 then return
            wset,list.drawindex
            res=CutPath(list.file,file=name)
            SaveImage,list.path,name+'.jpg','plotchis',list,OutID=list.OutID
            endcase
    'Load PDF':begin
            LoadPDD,ev,1
            widget_control,ev.top,get_uvalue=list
            RefreshPDD_OPLOTCHI,list
            RefreshDisplayCHIS,list
            endcase
  'Load Multiple PDF':begin
            LoadPDD,ev,2
            widget_control,ev.top,get_uvalue=list
            RefreshPDD_OPLOTCHI,list
            RefreshDisplayCHIS,list
            endcase
    'Resize':begin
            widget_control,ev.top,get_uvalue=list
            
            txt=string(list.xsize,list.ysize,format='("[",I0,",",I0,"]")')
            txt=promptnumber(txt,ev.top,'Dimensions:')
            tmp=strsplit(txt,'[], \t',/extract,count=ct)
            if ct ne 2 then return
            if tmp[0] le 0 or tmp[1] le 0 then return
            list.xsize=long(tmp[0])
            list.ysize=long(tmp[1])
            
            widget_control,list.draw,xsize=list.xsize,ysize=list.ysize
            wdelete,list.pixindex
            Window, /Free, xsize=list.xsize, ysize=list.ysize, /Pixmap
            list.pixindex = !D.Window
            
            widget_control,ev.top,set_uvalue=list
            
            RefreshPDD_OPLOTCHI,list
            RefreshDisplayCHIS,list
            endcase
    'View PDF':   begin
            WIDGET_CONTROL, ev.top, get_uvalue=list
            if list.PDFnTotal ne 0 then $
              EditPDFWindow,ev,'RefreshDisplayCHIS',od={PDFx:(*list.PDFx),PDFi2:(*list.PDFi2),PDFScale:(*list.PDFScale),PDFOffset:(*list.PDFOffset)}
              endcase
    'Load Mask':begin
             WIDGET_CONTROL, ev.top, get_uvalue=list
              LoadMask,list,1,list.OutID,error=error,/oD,flagi=['$IntCor2Dinfo:','$FIT_ROI:','$BACKoD:','$FIL:','$LABELS:','$MASK1:','$FITMODEL:']

              RefreshXVal,list
              widget_control,ev.top,set_uvalue=list
              RefreshDataOPLOTCHI,ev
              endcase
    'Y: logarithmic' :begin
            widget_control,ev.top,get_uvalue=list
            list.ylog=0
            widget_control,ev.top,set_uvalue=list
            RefreshDisplayCHIS,list
            widget_control,ev.id,set_value='Y: linear'
            endcase
    'Y: linear': begin
            widget_control,ev.top,get_uvalue=list
            list.ylog=1
            list.ypower=0
            list.ypowerstr=''
            widget_control,ev.top,set_uvalue=list
            RefreshDisplayCHIS,list
            widget_control,ev.id,set_value='Y: logarithmic'
            endcase
   'X: logarithmic' :begin
            widget_control,ev.top,get_uvalue=list
            list.xlog=0
            widget_control,ev.top,set_uvalue=list
            RefreshDisplayCHIS,list
            widget_control,ev.id,set_value='X: linear'
            endcase
    'X: linear': begin
            widget_control,ev.top,get_uvalue=list
            list.xlog=1
            widget_control,ev.top,set_uvalue=list
            RefreshDisplayCHIS,list
            widget_control,ev.id,set_value='X: logarithmic'
            endcase
    'Sum displayed': begin
            widget_control,ev.top,get_uvalue=list
            oplotchi_sum,list
            endcase
    'Surface plot': begin
            widget_control,ev.top,get_uvalue=list
            img=MakePlotType1Img(list,error=error,xnew=xnew)
            if error then return
            s=size(img)
            xr=list.xrange[*,list.xtype]
            yr=[0,list.ndata-1]
            if list.ndata le 1 then return
            zr=[min(img),max(img)]

            list2={tiff:img,zr:zr,xr:xr,yr:yr,xx:xnew,proplot:'plotchis3d'}
              flplot_obj,list2,GROUP_LEADER=ev.top,xtitle=list.xtitle[list.xtype],ytitle='Y(profilenr)',ztitle='Instensity'
            endcase
'Edit xrange': ModeOChi,ev,val
'Edit Parameters': EditOplotchiParam,ev
'Legend Position': ModeOChi,ev,val
'Edit Chi': EditChi,ev
'File Browsing':begin
              MakeBOptionsWindow,ev
              endcase
'Show All':    begin
            widget_control,ev.top,get_uvalue=list
            oplotchi_hide,list,/show
            RefreshDisplayCHIS,list
            endcase
'Show Next':begin
            widget_control,ev.top,get_uvalue=list
            oplotchi_hide,list
            RefreshDisplayCHIS,list
            endcase
'Show Prev':begin
            widget_control,ev.top,get_uvalue=list
            oplotchi_hide,list,/prev
            RefreshDisplayCHIS,list
            endcase
'Save session': oplotchi_savesession,ev
'Restore session': oplotchi_restoresession,ev
 else:        begin
             widget_control,ev.top,get_uvalue=list
            list.xtype++
            list.xtype mod= n_elements(list.xmodes)
            widget_control,ev.top,set_uvalue=list
            RefreshDisplayCHIS,list
            widget_control,ev.id,set_value='X: '+list.xmodes[list.xtype]
             endelse
    endcase
    endcase
endcase
end;pro oplotchi_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro oplotchi,ev=ev
print,'=========================================START 1D Profile Compare'

if keyword_set(ev) then begin
    widget_control,ev.top,get_uvalue=list
    path=list.chipath
    Platform=list.Platform
    file=list.chifile
    sortseparator=list.sortseparator
    sortmethod=list.sortmethod
    pddpath=list.pddpath
    pddfile=list.pddfile
    lambda=list.lambda
    cconv=list.cconv
    scrx=list.scrx
    scry=list.scry
    dist=list.dist
    a=list.a
    b=list.b
    phipixphi0=list.phipixphi0
    mskpath=list.mskpath
    mskfile=list.mskfile
    main=ev.top
endif else begin
    ControlDevice,Platform,Visual,Msize
    result=ReadINI('$XRD mskpath:')
    mskpath=result.(0)
    mskfile=result.(1)
    result=ReadINI('$XRD_CHI path:')
    path=result.(0)
    file=result.(1)
    sortseparator='_'
    sortmethod=1b
    result=ReadINI('$XRD pddpath:')
    pddpath=result.(0)
    pddfile=result.(1)
    lambda=0.
    cconv=double(4.13566743*299792458.*10.^(-8))
    scrx=0.
    scry=0.
    dist=0.
    a=0.
    b=0.
    phipixphi0=0.
    main=0L
endelse
error=CutPath(file,ext=ext)
Format=ext

result=ReadINI('$ChiSize:')
xsize=(result.(0))[0]
ysize=(result.(0))[1]
bartitle='1D Profile Compare: '
Topbase = WIDGET_BASE(MBAR=bar,/column,title=bartitle)

    baset=widget_base(Topbase,/column,uname='base2')
       draw = WIDGET_DRAW(baset, XSIZE=xsize, YSIZE=ysize,uname='draw')
       Window, /Free, xsize=xsize, ysize=ysize, /Pixmap
       pixindex = !D.Window
    basedyn=widget_base(Topbase,/column,uname='basedyn',uvalue='')
        baset=widget_base(basedyn,/row)
            label=widget_label(baset,value='PlotStyle  :     ',xsize=xs)
            text=widget_droplist(baset,value= ['Normal','Image']$
            ,uvalue='plotstyle',uname='plotstyle')
            plottype=0
            widget_control,text,SET_DROPLIST_SELECT=plottype
        baset=widget_base(basedyn,/column,uname='Mbase')
        baset=widget_base(baset,/row)
            b=widget_button(baset,value='Show All',xsize=xs)
            b=widget_button(baset,value='Show Next',xsize=xs)
            b=widget_button(baset,value='Show Prev',xsize=xs)

OutID=SubscribeToOutID()

list={listtype:30,$
        draw:draw,$
        drawindex:0L,$
        xsize:xsize,$
        ysize:ysize,$
        pixindex:pixindex,$
        path:path,$
        file:file,$
        pddpath:pddpath,$
        pddfile:pddfile,$
        mskfile:mskfile,$       ; name of .in file (vb.ca24.in)
        mskpath:mskpath,$       ; path where all files are saved and loaded(vb.D:\wout\files2\)
        Platform:Platform,$
        Format:Format,$
        OutID:OutID,$
        main:main,$

        sortseparator:sortseparator,$
        sortmethod:sortmethod,$
        lambda:lambda,$             ; X-ray lambda(Angstroms)
        cconv:cconv,$; E(keV)=12.3984/lambda(Angstroms)
        scrx:scrx,$                  ; Pixel x-size (m, input in um)
        scry:scry,$                  ; Pixel y-size (m, input in um)
        dist:dist,$                  ; Distance sample-detector (m, input in cm)
        a:a,$
        b:b,$
        phipixphi0:phipixphi0,$
        
        asciitemplate:{    rowskip:2l,$
                    delimiter:string(32b)+string(9b),$
                    ind:[1,2l]},$

        xs:0,$
        ys:0,$
        sysvar:PointSysVar(),$    ; System variable for coord. convertion
        colors:fix(16*[15.99,14,12,10,8,2]),$
        plottype:plottype,$        ; 0: normal plot, 1: tvscl
        scaleimg:100,$            ; Scale plottype1
        scross:1b,$                ; show cross on plottype1
        slegend:1b,$            ; show legend on plottype0
        bPDFpsym:0b,$            ; PDF with psym
        imgoffpr:0.08,$            ; Offset for plottype1
        normalizey:0b,$            ; Normalize profiles so there Y-range is the same
        ypower:0.,$                ; intensity^power for plotting
        ypowerstr:'',$            ; power string

        data:PTR_NEW(),$
        ndata:0,$
        nbins:0,$                ; maximum number of bins

        PDFn:ptr_new(),$        ; number of d-spacings in PDF file
        PDFd:ptr_new(),$        ; d-spacings in PDF file
        PDFx:ptr_new(),$        ; d, tt and radial distance
        PDFni:ptr_new(),$       ; number of intensities in PDF file
        PDFi:ptr_new(),$        ; intensities in PDF file scaled from 0-100
        PDFi2:ptr_new(),$        ; position intensities
        PDFhkl:ptr_new(),$        ; miller indices
        PDFname:ptr_new(),$     ; name of loaded compount
        PDFs:ptr_new(),$        ; first bit: show lines or not, second bit: show labels or not
        nPDF:0,$                ; Number of PDF files loaded
        PDFnTotal:0,$            ; Total number of reflections loaded
        PDFScale:ptr_new(),$    ; Use spe values when not in [0,1]
        PDFOffset:ptr_new(),$    ; Start of PDF lines
        PDFSearch:{type:0,PDFSubstr:'',PDFdi:0.,PDFde:0.},$    ; PDF subsearch

        ctl:-3,$    ; Color table
        xpos:0.8,$
        ypos:0.5,$
        incpos:0.03,$
        char:1.,$
        xlog:0b,$
        ylog:0b,$
        xrange:fltarr(2,4),$
        yrange:fltarr(2),$
        xtype:1,$
        xmodes:ChiXtitle(mode=1),$
        xunits:ChiXtitle(mode=2),$    ;Possible x units
        shortxunits:ChiXtitle(mode=3),$    ;Possible x units
        xtitle:ChiXtitle(),$;Possible x-titles
        ytitle:'Intensity (a.u.)',$
        bartitle:bartitle,$
        title:'Powder patterns'}

menu1 = WIDGET_BUTTON(bar, VALUE='File', /MENU)
       buttonopen = WIDGET_BUTTON(menu1, VALUE='Load chi',uname='Load chi')
       button = WIDGET_BUTTON(menu1, VALUE='Load Mask',uname='Load Mask')
       button = WIDGET_BUTTON(menu1, VALUE='Load PDF',uname='Load PDF')
       button = WIDGET_BUTTON(menu1, VALUE='Load Multiple PDF',uname='Load Multiple PDF')
       button = WIDGET_BUTTON(menu1, VALUE='Save session',uname='Save session',/separator)
       button = WIDGET_BUTTON(menu1, VALUE='Restore session',uname='Restore session')
       button = WIDGET_BUTTON(menu1, VALUE='Save Image',uname='Save Image',/separator)
       button = WIDGET_BUTTON(menu1, VALUE='Exit',uname='Exit',/separator)

menu2=WIDGET_BUTTON(bar, VALUE='Edit', /MENU,uname='Edit')
        
        button= widget_button(menu2,value='X: '+list.xmodes[list.xtype],uname='Xax')
        button= widget_button(menu2,value='Y: linear',uname='Y')
        button= widget_button(menu2,value='X: linear',uname='X')
        button= widget_button(menu2,value='Edit xrange',uname='Xrange',/separator)
        button= widget_button(menu2,value='Edit Parameters',uname='Z')
        button= widget_button(menu2,value='Legend Position',uname='Legend')
        button= widget_button(menu2,value='Edit Chi',uname='Edit Chi')
        button= widget_button(menu2,value='File Browsing')
        button= widget_button(menu2,value='View PDF',uname='View PDF')
        button= widget_button(menu2,value='Resize',uname='Resize')
        
menu3=WIDGET_BUTTON(bar, VALUE='Perform', /MENU,uname='Perform')
        button= widget_button(menu3,value='Surface plot',uname='Surface plot')
        button= widget_button(menu3,value='Sum displayed',uname='Sum displayed')

WIDGET_CONTROL,Topbase, /REALIZE,set_uvalue=list
widget_control,buttonopen,send_event={ID:buttonopen,TOP:Topbase,HANDLER:Topbase,SELECT:long(1)}
Xmanager,'oplotchi',Topbase,/NO_BLOCK,cleanup='oplotchi_cleanup'
end;pro oplotchi
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
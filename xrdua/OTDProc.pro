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

pro BackSub,image,back,Method,offset=offset
bs=size(image,/type) eq 10
if bs then begin
    s=size(*image)
    if ~array_equal(size(*image,/dim),size(*back,/dim)) then return
endif else begin
    s=size(image)
    if ~array_equal(size(image,/dim),size(back,/dim)) then return
endelse


if s[0] eq 2 then begin
    case Method of
    0: return ;No background
    else:begin ;1 Strip, 2 Dark image
       ; Scale
       ;itemp=image-back
       ;minim=min(itemp,indi)
       ;if minim ne 0 then begin
       ; indi=indi[0]
       ; scl=image[indi]/back[indi]
       ; image=image-scl*back
       ;endif else image=itemp
       ; Shift
       ;offset=min(image-back)
       ;image=image-((back+offset)>0)
       ; Cutoff
       if bs then (*image)=((*image)-(*back))>0 $
       else image=(image-back)>0
       endcase
    endcase
endif else begin
    case Method of
    -1:        return ; No background
    else:    begin ; 0 Ortpol, 2 Strip (1: Ortpol fit)
           ; Scale
           ;itemp=image-back
           ;minim=min(itemp,indi)
           ;if minim ne 0 then begin
           ; indi=indi[0]
           ; scl=image[indi]/back[indi]
           ; image=image-scl*back
           ;endif else image=itemp
           ; Shift
           ;offset=min(image-back)
           ;image=image-((back+offset)>0)
           ; Cutoff
              if bs then (*image)=((*image)-(*back))>0 $
              else image=(image-back)>0
           endcase
    endcase
endelse
end ;pro BackSub
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function minterpol,V, X, U, outer=outer,_EXTRA=e

Y = INTERPOL( V, X, U ,_EXTRA=e)

mi=min(X,max=ma)

if not keyword_set(outer) then begin
    ind=where(U ge mi and U le ma,ct)
    if ct eq 0 then outer=0 else outer=min(Y[ind])
endif

ind=where(U lt mi,ct)
if ct ne 0 then Y[ind]=outer
ind=where(U gt ma,ct)
if ct ne 0 then Y[ind]=outer

return,Y
end;function minterpol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SetAnySelect,val,type,top
widget_control,top,get_uvalue=list
if not ptr_valid(list.spet) then return

x=CHI_xval(val,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,type)
temp=min(abs((*list.xvalt)[*,list.xtypeorig]-x[list.xtypeorig]),index)
list.anyselect=[[x],[index]]

widget_control,top,set_uvalue=list
RefreshDisplayCHI,list

ID1=widget_info(top,FIND_BY_UNAME='x-value2')
str1=strtrim(string(format='(F15.3)',list.anyselect[list.xtype]),2)+list.xunits[list.xtype]
str2=strtrim(string(format='(F15.3)',(*list.spet)[list.anyselect[4],list.xtype]),2)
widget_control,ID1,set_value=[str1,str2]

end;pro SetAnySelect
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MEllipseSet,val,type,top,setself=setself
widget_control,top,get_uvalue=list
case type of
0:        begin
        res=BraggDtoX(val,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/pixm)
        list.ttmellipse=res[1]
        ID=widget_info(top,FIND_BY_UNAME='MEllipset')
        if widget_info(ID,/VALID_ID) then begin
            widget_control,ID,set_value=stringr(list.ttmellipse*180./!dpi)
            ID=widget_info(top,FIND_BY_UNAME='MEllipsefr')
            widget_control,ID,set_value=stringr(res[0])
            ID=widget_info(top,FIND_BY_UNAME='MEllipseq')
            widget_control,ID,set_value=stringr(res[2])
            if keyword_set(setself) then begin
                ID=widget_info(top,FIND_BY_UNAME='MEllipsed')
                widget_control,ID,set_value=stringr(val)
            endif
        endif
        endcase
1:        begin
        list.ttmellipse=val/180.*!dpi
        ID=widget_info(top,FIND_BY_UNAME='MEllipsed')
        if widget_info(ID,/VALID_ID) then begin
          res=BraggTtoX(list.ttmellipse,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/pixm)
          widget_control,ID,set_value=stringr(res[1])
          ID=widget_info(top,FIND_BY_UNAME='MEllipsefr')
          widget_control,ID,set_value=stringr(res[0])
          ID=widget_info(top,FIND_BY_UNAME='MEllipseq')
           widget_control,ID,set_value=stringr(res[2])
          if keyword_set(setself) then begin
              ID=widget_info(top,FIND_BY_UNAME='MEllipset')
              widget_control,ID,set_value=stringr(val)
          endif
        endif
        endcase
2:        begin
        res=BraggPtoX(val,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/pixm)
        list.ttmellipse=res[0]
        ID=widget_info(top,FIND_BY_UNAME='MEllipsed')
        if widget_info(ID,/VALID_ID) then begin
          widget_control,ID,set_value=stringr(res[1])
          ID=widget_info(top,FIND_BY_UNAME='MEllipset')
          widget_control,ID,set_value=stringr(list.ttmellipse*180/!dpi)
          ID=widget_info(top,FIND_BY_UNAME='MEllipseq')
          widget_control,ID,set_value=stringr(res[2])
          if keyword_set(setself) then begin
              ID=widget_info(top,FIND_BY_UNAME='MEllipsefr')
              widget_control,ID,set_value=stringr(val)
          endif
        endif
        endcase
3:        begin
        res=BraggQtoX(val,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/pixm)
        list.ttmellipse=res[1]
        ID=widget_info(top,FIND_BY_UNAME='MEllipset')
        if widget_info(ID,/VALID_ID) then begin
            widget_control,ID,set_value=stringr(list.ttmellipse*180./!dpi)
            ID=widget_info(top,FIND_BY_UNAME='MEllipsefr')
            widget_control,ID,set_value=stringr(res[0])
            ID=widget_info(top,FIND_BY_UNAME='MEllipsed')
            widget_control,ID,set_value=stringr(res[2])
            if keyword_set(setself) then begin
                ID=widget_info(top,FIND_BY_UNAME='MEllipseq')
                widget_control,ID,set_value=stringr(val)
            endif
        endif
        endcase
endcase
list.smellipse=1b
ID=widget_info(top,FIND_BY_UNAME='warning')
if widget_info(ID,/VALID_ID) then begin
    tmp=ConicTypeTT(list.a,list.ttmellipse,str=str)
    widget_control,ID,set_value=str
endif
RefreshDisplay, top,list

if ptr_valid(list.CHIchild) and list.CrossUpdate then begin
    if not keyword_set(setself) then setself=-1L
    for i=0l,n_elements(*list.CHIchild)-1 do begin
        if (*list.CHIchild)[i] ne setself then SetAnySelect,val,type,(*list.CHIchild)[i]
    endfor
endif

end;pro MEllipseSet
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro UnlinkAreas,mask,masknr,ind
count=n_elements(ind)

; Extract list of linknr for each mask
if (*mask)[2,masknr-1] eq 1 then offset=masknr else offset=masknr-1
linknr=(*mask)[1,0:offset-1]
if offset ne 1 then linknr=transpose(linknr)

; Assign new link numbers
m=max(linknr)
newlinknr=indgen(count)+m+1
linknr[ind]=newlinknr
; Loop over all link number and cope
; with blank missing numbers
for i=0l,m+count do begin
    ind2=where(linknr eq i,count2)
    if count2 eq 0 then begin
        ind2=where(linknr gt i,count2)
        if count2 ne 0 then begin
            linknr[ind2]=linknr[ind2]-1
            i=i-1
        endif
    endif
endfor
(*mask)[1,0:offset-1]=linknr

end;pro UnlinkAreas
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LinkAreas,mask,masknr,ind

; Extract list of linknr for each mask
if (*mask)[2,masknr-1] eq 1 then offset=masknr else offset=masknr-1
linknr=(*mask)[1,0:offset-1]
if offset ne 1 then linknr=transpose(linknr)

; Assign new link number
newlinknr=max(linknr)+1
linknr[ind]=newlinknr

; Loop over all link number and cope
; with blank missing numbers
for i=0l,newlinknr do begin
    ind=where(linknr eq i,count)
    if count eq 0 then begin
        ind=where(linknr gt i,count)
        if count ne 0 then begin
            linknr[ind]=linknr[ind]-1
            i=i-1
        endif
    endif
endfor
(*mask)[1,0:offset-1]=linknr

end;pro LinkAreas
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_EditMask, ID

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
widget_control,list.top,sensitive=1
end;pro CleanUp_EditMask
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditMaskWindowP,ev,list

if size(ev,/type) eq 8 then begin
    TLB=ev.top
    ID=widget_info(TLB,FIND_BY_UNAME='basegroups')
    WIDGET_CONTROL, ID, /DESTROY
    wrapbase=widget_info(TLB,FIND_BY_UNAME='wrapbase')
endif else begin
    wrapbase=ev
endelse

result=ReadINI('$EditMask:')

basegroups=widget_base(wrapbase,/column,uname='basegroups',/scroll,$
    y_scroll_size=result.(1)/2,x_scroll_size=result.(1))
    if list.masknr eq 0 then return
    if n_elements(TLB) ne 0 then widget_control,TLB,update=0
    ; Make different rows, depending on how many (nr) groups you want on a row
    nr=result.(0)
    nbase=ceil(float(list.masknr)/nr)
    base6=lonarr(nbase)
    for i=0l,nbase-1 do base6[i]=widget_base(basegroups,/row)


    ; Loop over all groups and place them in the right row
    rbutton = lonarr(list.masknr)
    j=0
    for i=0l,list.ngr-1 do begin
       val=stringr((*list.linkdiff)[i])
       ; Base for the group
       bgroups=widget_base(base6[j],/column)
         ind=where((*list.mask)[1,*] eq (*list.linkdiff)[i],ct)

         ; Group Base
         baset=widget_base(bgroups,/nonexclusive,/row,uvalue=ind)
          button=widget_button(baset,value='group '+val,$
              uvalue=(*list.linkdiff)[i],/ALIGN_CENTER,sensitive=type,uname='1g'+val)
          button=widget_button(baset,value='',uvalue=(*list.linkdiff)[i],uname='2g'+val)
         ; Area bases

         for k=0l,ct-1 do begin
          val=stringr(ind[k])
          baset=widget_base(bgroups,/nonexclusive,/row)
              button=widget_button(baset,value='area '+val,$
                 uvalue=ind[k],sensitive=type,uname='1a'+val)
              rbutton[ind[k]]=widget_button(baset,value='',uvalue=ind[k],uname='2a'+val)
         endfor
       if ((i+1) mod nr) eq 0 then j=j+1
    endfor
for i=0l,list.masknr-1 do widget_control,rbutton[i],set_button=(*list.smask)[i]

if n_elements(TLB) ne 0 then widget_control,TLB,update=1
end;pro EditMaskWindowP
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ExtractLinkNr,list

link=reform((*list.mask)[1,*])
linkdiff=[0]
for i=0l,list.masknr-1 do begin
     ind=where(linkdiff eq link[i],count)
     if count eq 0 then linkdiff=[linkdiff,link[i]]
endfor
linkdiff=long(linkdiff)
linkdiff=linkdiff[sort(linkdiff)]
list.ngr=n_elements(linkdiff)
(*list.linkdiff)=temporary(linkdiff)

end;pro ExtractLinkNr
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditMaskHandler,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control, ev.top, get_uvalue=list
widget_control,ev.id,get_value=val

case val of
 'Sort':begin

        ind=where((*list.select) eq 1,count)
         case count of
          0 : return
              else: begin

              ; Calc mean d-spacings
              D=fltarr(count)
              case list.refreshcall of
              0:    begin
                      for i=0l,count-1 do begin
                            case (*list.mask)[0,ind[i]] of
                            1:    begin
                                XY=(*list.mask)[3:*,ind[i]]
                                res=BraggXYtoX(XY[[0,0,1,1]],XY[[2,3,2,3]],list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,/onlyd)
                                dspacB=min(res,max=dspacE)
                                D[i]=(dspacB+dspacE)/2.
                                endcase
                            2:    begin
                                D[i]=((*list.mask)[3,ind[i]]+(*list.mask)[5,ind[i]])/2
                                endcase
                            3:    begin
                                D[i]=((*list.mask)[3,ind[i]]+(*list.mask)[4,ind[i]])/2
                                endcase
                            endcase
                      endfor
                     endcase
              1:    begin
                      D=total((*list.mask)[3:4,ind],1)/2.
                      endcase
              endcase

              ; Sort
              tmp1=(*list.mask)[*,ind[reverse(sort(D))]]
              tmp2=ShrinkArray((*list.mask),ind,/ROW)
              if n_elements(tmp2) le 1 then (*list.mask)=tmp1 else (*list.mask)=[[tmp1],[tmp2]]
              UnlinkAreas,list.mask,list.masknr,indgen(count)

              endelse
         endcase

        ExtractLinkNr,list
        (*list.select)[*]=0
        widget_control, ev.top, set_uvalue=list
        EditMaskWindowP,ev,list
         end
 'Link':    begin
       ind=where((*list.select) eq 1,count)
         case count of
          0 : return
              else: begin
              LinkAreas,list.mask,list.masknr,ind
              endelse
         endcase

         ExtractLinkNr,list
        (*list.select)[*]=0
        widget_control, ev.top, set_uvalue=list

         ; Make new window
         EditMaskWindowP,ev,list
       endcase
 'UnLink':  begin
       ind=where((*list.select) eq 1,count)
         case count of
          0 : return
             else: begin
              UnlinkAreas,list.mask,list.masknr,ind
              endelse
           endcase

         ExtractLinkNr,list
        (*list.select)[*]=0
        widget_control, ev.top, set_uvalue=list

         ; Make new window
         EditMaskWindowP,ev,list
       endcase
 'Del':     begin
         ; Make new nask, delete the selected ones
         ind=where((*list.select) eq 0,count)
         case count of
          0  : begin
                  ptr_free,list.linkdiff
                  ptr_free,list.select
                  ptr_free,list.smask
                  ptr_free,list.mask
                  list.masknr=0
                  list.ngr=0

                ; Make new window
                widget_control, ev.top, set_uvalue=list
                EditMaskWindowP,ev,list
                return
               endcase
             1  : begin
                 (*list.mask)=(*list.mask)[*,ind]
                 (*list.smask)=(*list.smask)[ind]
                 list.masknr=count
                 (*list.mask)[1,0]=0
               endcase
             else: begin
                 (*list.mask)=(*list.mask)[*,ind]
                 (*list.smask)=(*list.smask)[ind]
                 list.masknr=count
                 link=(*list.mask)[1,*]
                 indtemp=where((link-(shift(link,1))) eq 0,count2)
                 ind=link*0
                 if count2 ne 0 then ind[indtemp]=1
                 ind[0]=0
                 link=0
                 for i=0l,count-1 do begin
                   if ind[i] eq 1 then (*list.mask)[1,i]=(*list.mask)[1,i-1] $
                   else begin
                   (*list.mask)[1,i]=link
                   link=link+1
                   endelse
                 endfor
               endelse
           endcase

        ExtractLinkNr,list
        (*list.select)=bytarr(list.masknr)
        widget_control, ev.top, set_uvalue=list

        ; Make new window
        EditMaskWindowP,ev,list
       endcase
 'Select All':begin
         (*list.select)[*]=ev.select
         for i=0l,list.ngr-1 do begin
          ID=widget_info(ev.top,FIND_BY_UNAME='1g'+stringr(i))
          widget_control,ID,set_button=ev.select
         endfor
         for i=0l,list.masknr-1 do begin
          ID=widget_info(ev.top,FIND_BY_UNAME='1a'+stringr(i))
          widget_control,ID,set_button=ev.select
         endfor
         widget_control,ev.top,set_uvalue=list
         endcase
  'View All':begin
         (*list.smask)[*]=ev.select
         for i=0l,list.ngr-1 do begin
          ID=widget_info(ev.top,FIND_BY_UNAME='2g'+stringr((*list.linkdiff)[i]))
          widget_control,ID,set_button=ev.select
         endfor
         for i=0l,list.masknr-1 do begin
          ID=widget_info(ev.top,FIND_BY_UNAME='2a'+stringr(i))
          widget_control,ID,set_button=ev.select
         endfor
         widget_control,ev.top,set_uvalue=list
         endcase
  'Apply':  begin
         widget_control,list.top,get_uvalue=listmain
         ptr_free,listmain.smask
         if ptr_valid(list.smask) then listmain.smask=PTR_NEW(*list.smask)
         listmain.masknr=list.masknr
         ptr_free,listmain.mask
         if ptr_valid(list.mask) then listmain.mask=PTR_NEW(*list.mask)
         case list.refreshcall of
         0:    RefreshDisplay,list.top,listmain
         1: begin
             widget_control,list.top,set_uvalue=listmain
             RefreshDisplayCHI,listmain
             endcase
         endcase
         endcase
   'OK':    begin
         widget_control,list.top,get_uvalue=listmain
         ptr_free,listmain.smask
         if ptr_valid(list.smask) then listmain.smask=PTR_NEW(*list.smask)
         listmain.masknr=list.masknr
         ptr_free,listmain.mask
         if ptr_valid(list.mask) then listmain.mask=PTR_NEW(*list.mask)
         case list.refreshcall of
         0:    RefreshDisplay,list.top,listmain
         1: begin
             widget_control,list.top,set_uvalue=listmain
             RefreshDisplayCHI,listmain
             endcase
         endcase
           WIDGET_CONTROL, ev.top, /DESTROY
         endcase
  'Cancel': WIDGET_CONTROL, ev.top, /DESTROY
  else:     begin
         ID=widget_info(ev.id,/parent)
         widget_control,ID,get_uvalue=uval
         nuval=n_elements(uval)

         if nuval eq 0 then begin ; Area selection
              widget_control,ev.id,get_uvalue=i
              if val eq '' then (*list.smask)[i]=ev.select else $
              (*list.select)[i]=ev.select
           endif else begin ; Whole group selection
              if val eq '' then begin
                pref='2a'
                (*list.smask)[uval]=ev.select
              endif else begin
                pref='1a'
                (*list.select)[uval]=ev.select
              endelse
              for i=0l,nuval-1 do begin
                ID=widget_info(ev.top,FIND_BY_UNAME=pref+stringr(uval[i]))
                widget_control,ID,set_button=ev.select
              endfor
           endelse
           widget_control,ev.top,set_uvalue=list
         endcase
endcase

end;pro EditMaskHandler
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditMaskWindow,ev,type,refreshcall

; type=1 => editing
; type=0 => viewing

widget_control,ev.top,get_uvalue=list
if (list.masknr eq 0) then return
mask=(*list.mask)
smask=(*list.smask)

; Extract different link numbers
link=reform(mask[1,*])
linkdiff=[0]
for i=0l,list.masknr-1 do begin
    ind=where(linkdiff eq link[i],count)
    if count eq 0 then linkdiff=[linkdiff,link[i]]
endfor
linkdiff=long(linkdiff)
linkdiff=linkdiff[sort(linkdiff)]
ngr=n_elements(linkdiff)

WIDGET_CONTROL, ev.top, sensitive=0

case refreshcall of
0:    begin
    center=list.center
    a=list.a
    b=list.b
    endcase
1:    begin
    center=0b
    a=0b
    b=0b
    endcase
endcase

list={top:ev.top,$
       masknr:list.masknr,$
       ngr:ngr,$
       select:PTR_NEW(bytarr(list.masknr)),$
       smask:PTR_NEW(smask),$
       linkdiff:PTR_NEW(linkdiff),$
       refreshcall:refreshcall,$
       mask:PTR_NEW(mask),lambda:list.lambda,dist:list.dist,scrx:list.scrx,scry:list.scry,center:center,a:a,b:b}
base=widget_base(/column,title='Edit Mask',uvalue=list)

basebuttons=widget_base(base,/row)
    button=widget_button(basebuttons,value='Link',uname='Link',sensitive=type)
    button=widget_button(basebuttons,value='UnLink',uname='UnLink',sensitive=type)
    button=widget_button(basebuttons,value='Del',uname='Del',sensitive=type)
    button=widget_button(basebuttons,value='Sort',uname='Sort',sensitive=type)
    baset=widget_base(basebuttons,/row,/nonexclusive)
       sbutton=widget_button(baset,value='Select All',uname='Select All',sensitive=type)
       vbutton=widget_button(baset,value='View All',uname='View All')
    button=widget_button(basebuttons,value='OK',uname='OK')
    button=widget_button(basebuttons,value='Apply',uname='Apply')
    button=widget_button(basebuttons,value='Cancel',uname='Cancel',sensitive=type)

wrapbase=widget_base(base,uname='wrapbase')
EditMaskWindowP,wrapbase,list

WIDGET_CONTROL, base, /REALIZE
widget_control,sbutton,set_button=0
widget_control,vbutton,set_button=0

Xmanager,'EditMaskWindow',base, event_handler='EditMaskHandler',$
    cleanup='CleanUp_EditMask',GROUP_LEADER=ev.top

end;pro EditMaskWindow
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mSavgol,data=data,r=r,nr=nr,nl=nl,order=order,coeff=coeff,fvar=fvar,sigma2=sigma2

; r: degree of the interpolating polynomial
; order: make a kernel to calculate the $order$th derivative of a dataset (order=0 => smoothing) (two values for 2D datasets)
; nr,nl: number of points right and left to the central kernel point
; Data: pointer to 1D or 2D array to be filtered (if data is specified, mSavgol will return the filtered data, else it will return the kernel)
; coeff: calculated kernel
; fvar: when set, variance on filtered data is calculated
; sigma2: input = variance on data (set to data if not specified); output = variance on filtered data

; E.g.: kernel for smoothing: kernel=mSavGol(order=0,nr=4,nl=4,r=2)
;        use kernel: smdata=convol(data,kernel,/edge_trunc)
;
;        2D analog: kernel=mSavGol(order=[0,0],nr=4,nl=4,r=2)
;                   smdata=convol(data,kernel,/edge_trunc)
;
; Remark:    boxcar <> r=0
;            cubic <> r=3
;            quintic <> r=5
;
; mSavGol(data=data,order=[0,0],nr=1,nl=1,r=0) is the same as smooth(data,3)

;----Evaluate input----
DataIsPtr=0b
nData=n_elements(data)
if nData ne 0 then begin
    if size(data,/TYPE) eq 10 then begin
        sdata=size(*data,/STRUCTURE)
        DataIsPtr=1b
    endif else sdata=size(data,/STRUCTURE)
endif else sdata=size(data,/STRUCTURE)

if not keyword_set(r) then r=2
if not keyword_set(nr) then nr=2
if not keyword_set(nl) then nl=2
if not keyword_set(order) then order=intarr(sdata.N_DIMENSIONS>1)

sorder=n_elements(order)

nr=nr>1
nl=nl>1
r=r>2
order= 0 > order < r
dim=sorder
dim=1>dim<2
w=nr+nl+1

;----Set up Design Matrix----
case dim of
1:  begin
    order=order[0]
    ncol=r+1
    nrow=w
    s=indgen(nrow)-nl
    one=fltarr(ncol)+1
    A=one#s
    for j=0l,r do A[j,*]=A[j,*]^j
    endcase
2:  begin
    ncol=(r+1)*(r+2)/2
    nrow=w*w
    s=indgen(w)-nl
    one=intarr(w)+1
    xs=reform(s#one,w*w)
    ys=reform(one#s,w*w)
    ntr=0
    A=fltarr(ncol,nrow)
    searchmat=intarr(r+1,r+1)
    for j=0l,r do begin
        for i=0l,j do begin
        A[ntr,*] = xs^(j-i)*ys^i
        searchmat[i,j-i]=ntr
        ntr=ntr+1
        endfor
    endfor
    endcase
endcase

;----Calculate Coefficients----
AT=transpose(A)
C=invert(AT##A)##AT
; each row of C corresponds with an a
; for 1 dimensions:
; a = (a0 a1 a2 a3 ...)T
; corresponding filter  (f f' f'' ...)/order!
; order 0:
; for 2 dimensions:
; a = (a00 a10 a01 a20 a11 a02 a30 a21 a12 a03 ...)T
; corresponding filter  (f df/dx df/dx df/dy d2f/dx2 df/dxdy d2f/dy2 ...)/(order[0]!*order[1]!)
; order[0]: order of derivative to x
; order[1]: order of derivative to y
; eg: order=[1,0] <> df/dx <> a10 <> searchmat[0,1]

;----Extract desired coeff and multiplied by appropriate constant----
case dim of
1:  begin
    ofact=1.0d
    for k = order+0.0d, 1d, -1d do ofact = ofact * k
    fil=C[*,order]*ofact
    endcase
2:  begin
    ofact=dblarr(dim)+1.
    for l=0l,dim-1 do $
       for k = order[l]+0.0d, 1d, -1d do ofact[l] = ofact[l] * k
    ind=searchmat[order[1],order[0]]
    fil=reform(C[*,ind],w,w)*ofact[0]*ofact[1]
    endcase
endcase

;----Output----
coeff=fil
if nData ne 0 then begin
    if sdata.N_DIMENSIONS ne sorder then message,'Dimensions of Data and Order don�t correspond'

    if sdata.type ne 4 and sdata.type ne 5 then begin
        data2=TypeConvert(convol(float((DataIsPtr?(*data):data)),fil,/edge_trunc),sdata.type)
    endif else data2=convol((DataIsPtr?(*data):data),fil,/edge_trunc)
    if keyword_set(fvar) then begin
       if not keyword_set(sigma2) then sigma2=(DataIsPtr?(*data):data)
       ssigma2=size(sigma2,/STRUCTURE)
       if ssigma2.type ne 4 and ssigma2.type ne 5 then begin
               fvar=TypeConvert(convol(sigma2,fil^2.,/edge_trunc),ssigma2.type)
       endif else fvar=convol(sigma2,fil^2.,/edge_trunc)
    endif

    if DataIsPtr then begin
        p=PTR_NEW(0)
        *p=temporary(data2)
        return,p
    endif else return,data2
endif else return,fil

end;function mSavgol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
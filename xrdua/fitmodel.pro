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

function ExpData,list,xd=xd,xtt=xtt,y=y,dlim=dlim,lambda=lambda,$
                LGPtype=LGPtype,dist=dist,zerotype=zerotype,promptdlim=promptdlim
                

if size(list,/type) ne 8 then begin
if widget_info(list,/valid_id) then begin
    widget_control,list,get_uvalue=top
    widget_control,top,get_uvalue=list
endif
endif

result={expdata:1b}
if keyword_set(dlim) then begin
    dlim=(*list.xvalt)[[list.backranind[0],list.backranind[1]],0]
    if keyword_set(promptdlim) then begin
        ttmax=(*list.xvalt)[list.backranind[1],1]
        ttmax=promptnumber(stringr(ttmax),list.main,'Maximal 2theta (deg):')
        dlim[1]=BraggTtoX(float(ttmax),list.lambda,/onlyd,/angledeg)
    endif 
    result=create_struct(result,'dlim',dlim[sort(dlim)])
endif

breverse=0b
i1=list.backranind[0]
i2=list.backranind[1]
if keyword_set(xd) then begin
    if (*list.xvalt)[i1,0] gt (*list.xvalt)[i2,0] then breverse=1b
    tmp=(*list.xvalt)[i1:i2,0]
    if breverse then tmp=reverse(tmp)
    result=create_struct(result,'xd',tmp)
    if keyword_set(y) then begin
        tmp=(*list.spet)[i1:i2,0]
        if breverse then tmp=reverse(tmp)
        result=create_struct(result,'yd',tmp)
    endif
endif
if keyword_set(xtt) then begin
    if (*list.xvalt)[i1,1] gt (*list.xvalt)[i2,1] then breverse=1b
    tmp=(*list.xvalt)[i1:i2,1]
    if breverse then tmp=reverse(tmp)
    result=create_struct(result,'xtt',tmp)
    if keyword_set(y) then begin
        tmp=(*list.spet)[i1:i2,1]
        if breverse then tmp=reverse(tmp)
        result=create_struct(result,'ytt',tmp)
    endif
endif

if keyword_set(LGPtype) then begin
    result=create_struct(result,'LGPtype',ExtractbitCode((*list.ptrModel).ptrFit,1))
    result=create_struct(result,'pol_dmono',*((*(*list.ptrModel).ptrFit).misc[5]))
    result=create_struct(result,'pol_Plindeg',*((*(*list.ptrModel).ptrFit).misc[6]))
    result=create_struct(result,'pol_n',*((*(*list.ptrModel).ptrFit).misc[8]))
    result=create_struct(result,'pol_mmult',*((*(*list.ptrModel).ptrFit).misc[9]))
    result=create_struct(result,'pol_s',*((*(*list.ptrModel).ptrFit).misc[10]))
    result=create_struct(result,'pol_azimuth',*((*(*list.ptrModel).ptrFit).misc[11]))
    lambda=1b
endif

if keyword_set(lambda) then $
    result=create_struct(result,'lambda',list.lambda)
    
if keyword_set(dist) then $
    result=create_struct(result,'dist',list.dist)
    
if keyword_set(zerotype) then $
    result=create_struct(result,'zerotype',ExtractbitCode((*list.ptrModel).ptrFit,3))
    
return,result
end;function ExpData
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_SetScalingToWeight,ID
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
end;pro CleanUp_SetScalingToWeight
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SetScalingToWeight_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
1 : begin ; button event
    widget_control,ev.id,get_value=val
    case val of
    'Calculate':begin
            ID=widget_info(ev.top,find_by_uname='K')
            widget_control,ID,get_value=K
            ID=widget_info(ev.top,find_by_uname='X')
            widget_control,ID,get_value=X
            
            widget_control,ev.top,get_uvalue=list2
            widget_control,list2.top,get_uvalue=top
            widget_control,top,get_uvalue=list
            
            ScalingFromWeight,reform(X),(*list.ptrModel).ptrFit,float(K),list.lambda,O=list2.O
            endcase
    'Cancel':
    else: return
    endcase
    endcase
3:    begin ; text event
    return
    endcase
9:    begin ; tabel event
    widget_control,ev.id,get_value=X
    ID=widget_info(ev.top,find_by_uname='total')
    widget_control,ID,set_value=stringr(total(X))
    return
    endcase
else: return
endcase

widget_control,ev.top,/destroy

end;pro SetScalingToWeight_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SetScalingToWeight,top,O=O

if WIDGET_INFO(top,/VALID_ID) eq 1 then begin
    WIDGET_CONTROL, top, sensitive=0
    GROUP_LEADER=top
    modal=1
endif else begin
    modal=0 ; no modal when no group leader
endelse

if ~keyword_set(O) then O=0
device,get_screen_size=screen
base=widget_base(title='Initialize scaling factors',/column,uvalue={top:top,O:O},GROUP_LEADER=GROUP_LEADER,$
            modal=modal,xoffset=screen[0]/2.,yoffset=screen[1]/2.)

widget_control,top,get_uvalue=top
widget_control,top,get_uvalue=list

; Get weight fractions and names
ptrFit=(*list.ptrModel).ptrFit
ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if ~(*ptr).include then begin
        ptr=(*ptr).next
        continue
    endif
    
    ; Extract weight fraction (first element in propI or propO)
    if keyword_set(O) then tmp=(*(*(*ptr).propO)[(*ptrFit).icurrent])[0] $
    else tmp=(*(*ptr).propI)[0]
    
    ; Append name and weight fraction
    if n_elements(names) eq 0 then begin
        names=(*ptr).name
        X=tmp
    endif else begin
        names=[names,(*ptr).name]
        X=[X,tmp]
    endelse
    
    ptr=(*ptr).next
endwhile

if n_elements(X) le 1 then return

table=widget_table(base,column_labels=['Weight%'],row_labels=names,$
                                value=transpose(X),uname='X',XSIZE=1,YSIZE=n_elements(X),/editable)
baset=widget_base(base,/row)
    label=widget_label(baset,value='Total:',/ALIGN_LEFT)
    text=widget_text(baset,value=stringr(total(X)),uname='total')

K=median(ExperimentalKfactor(ptrFit,list.lambda,O=O))
if K eq 0 then K=1.
baset=widget_base(base,/row)
    label=widget_label(baset,value='Experimental multiplier:',/ALIGN_LEFT)
    text=widget_text(baset,value=stringr(K),uname='K',/editable)

baset=widget_base(base,/row)
    button=widget_button(baset,value='Calculate')
    button=widget_button(baset,value='Cancel')

widget_control,base,/realize
Xmanager,'SetScalingToWeight',base, event_handler='SetScalingToWeight_Event',$
    cleanup='CleanUp_SetScalingToWeight',GROUP_LEADER=GROUP_LEADER

end;pro SetScalingToWeight
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TestTypeObsolete,ptr,dlim
obsolete=0b
case (*ptr).type of
10:
else:begin
    celparam=(*(*ptr).globalparamIO.I)[indCelParam(ptr)]
    hkldata=GenHKL(celparam,dlim,*(*ptr).pstrinfo)
    obsolete=hkldata.n ne (*ptr).peaknr
    endcase
endcase
return,obsolete
end;function TestTypeObsolete
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ProfileQuantile,profiletype,param,p

; xcut = x cutoff value
;      = xpos + xoff
; P = probability
; p = area percentage of PDF to take into account
; K1 = P(X <= xcut) = (p+1)/2
; K2 = P(X > xcut) = (1-p)/2

; probability density function (PDF): PDF(x)
; cumulative distribution function (CDF): CDF(xcut) = P(X <= xcut) = K1
; quantile function (CDF^-1): CDF^-1(K1) = xcut

; IDL: *_PDF functions are "probability distribution functions" which are CDFs !
;      *_CVF functions are quantile functions


case profiletype of
0: begin
    ; Gaussian
    ;     CDF^-1(K1) = IMSL_ERF(p,/inverse)*sqrt(2)*sigma     (for xpos=0)
    ;               = GAUSS_CVF(K2)*sigma                    (GAUSS_CVF is for xpos=0 and sigma=1)
    ; param=FWHM
    xoff=GAUSS_CVF((1-(0>p<1))/2.)*param/(2*sqrt(2*alog(2)))
    endcase
1:
2:
3:
4:
5:
6:
endcase

if n_elements(xoff) eq 0 then xoff=0.

return,xoff

end;function ProfileQuantile
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TestSaveTypeX,list,ptr,type=type

peaknr=(*ptr).peaknr
if peaknr eq 0 then return
if ~keyword_set(type) then type=0

; ----File to save results----
path_pos = rstrpos(list.file,'.')  ; find first dot
if(path_pos eq -1) then path_pos = strlen(list.file)     ; no extention -> append

af=(type eq 2)?'*.csv':'*.csv;*.chi'
name = strmid(list.file,0,path_pos)+'_param.csv'
path=SelFile(list.path,name,af,'Save Results To....')
if path EQ '' then return
error=CutPath(path,path=pathtemp,ext=ext)
list.path=pathtemp

; ----Extract results----
matrix=gfuncPeakparam(list,ptr,/O,/tt,d=type eq 2)

; Save result
case type of
0:    begin; Save 2-theta and FWHM
    ; Replace area by FWHM
    matrix=matrix[[0,2],*]
    colnames=['2-theta','FWHM']
    endcase
1:    begin; Save 2-theta and Area/sqrt(background_area)
    ptrFit=(*list.ptrModel).ptrFit
    if ~ptr_valid((*ptrFit).yb) then begin
        printw,list.OutID,'No background in model.'
        return
    endif
    x=(*list.xvalt)[list.backranind[0]:list.backranind[1],1] ;two-theta
    y=(*list.spet)[list.backranind[0]:list.backranind[1],1]
    yb=(*(*ptrFit).yb)
    if n_elements(x) ne n_elements(yb) then begin
        printw,list.OutID,"Background doesn't have the right size."
        return
    endif
    
    ; Find cutoff values for background integration
    p=float(promptnumber('95',list.main,'Peak area to take into account (%):'))/100.
    if p eq 0 then return
    off=ProfileQuantile(ExtractbitCode(ptr,0),matrix[2,*],p)
    
    matrix=matrix[0:1,*]
    i0=value_locate(x,matrix[0,*]-off)
    i1=value_locate(x,matrix[0,*]+off)

    n=n_elements(i0)
    backint=fltarr(n)
;    window
    for i=0l,n-1 do begin
        j0=i0[i]
        j1=i1[i]
        if j1 gt j0 then begin
            backint[i]=int_tabulated(x[j0:j1],yb[j0:j1])
;            plot,x,y,xrange=x[[j0,j1]]*[0.95,1.05]
;            oplot,x,yb,color=100
;            plots,x[[j0,j0]],!y.crange
;            plots,x[[j1,j1]],!y.crange
;            pxval=x[[j0,j0,j1,j1]]
;            pyval=[0,yb[[j0,j1]],0]
;            polyfill,pxval,pyval,/line_fill,orientation=45,col=100,/data
            
            ; Should not be y but the fitted single peak
;            pxval=[x[j0],x[j0:j1],x[j1]]
;            pyval=[yb[j0],y[j0:j1],yb[j1]]
;            polyfill,pxval,pyval,/line_fill,orientation=45,col=10,/data
;            wait,0.1
        endif
    endfor
    
    ; Modify intensity
    matrix[1,*]*=p/sqrt(backint)
    colnames=['2-theta','Area/sqrt(BackgroundArea)']
    endcase
2:    begin; save all peak param
    ind=*((*ptr).paramind[6,0])
    colnames=(*ptr).peaklabels[ind]
    colnames=[colnames[0],'d-spacings',colnames[1:*]]
    endcase
endcase

; Sort using tt
ind=sort(matrix[0,*])
matrix=matrix[*,ind]

case ext of
'.csv': mWRITE_CSV,path,matrix,colnames=colnames
'.chi': temp=WriteChi(path,matrix,1,list.path,list.OutID,/orig)
endcase
printw,list.OutID,'Saved '+path
end;pro TestSaveTypeX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SavePKSTypeX,list,ptr

peaknr=(*ptr).peaknr
if peaknr eq 0 then return

path_pos = rstrpos(list.file,'.')  ; find first dot
if(path_pos eq -1) then path_pos = strlen(list.file)     ; no extention -> append

name = strmid(list.file,0,path_pos)+'.pdd'
path=SelFile(list.path,name,['*.pdd','*.txt'],'Save Results To....')
if path EQ '' then return
error=CutPath(path,path=pathtemp,ext=ext)
list.path=pathtemp

pcsiwin=ext eq '.txt'

; Open
if openw_safe(lun, path, width=pcsiwin?120:0) then begin
    printw,list.OutID,path+' not saved'
    return
endif

; Peak position and intensity
; TODO: change LGP with lambda=CuKa??
matrix=gfuncPeakparam(list,ptr,/O,/d)
matrix=matrix[0:1,*]

; HKL
bhkl=~pcsiwin and (*ptr).type ne 10
if bhkl then begin
    icurrent=0
    p=(*(*ptr).peakparamIO.O)[icurrent]
    matrix=[matrix,(*p)[*(*ptr).paramind[0,0],*]]
endif

; Only finite
ind=where(finite(matrix[0,*]) and finite(matrix[1,*]) and (matrix[0,*] ge 0) and (matrix[1,*] ge 0),ct)
if ct ne 0 then begin
    matrix=matrix[*,ind]
    peaknr=ct
endif

; Only 40 peaks (Pcsiwin)
if pcsiwin then nmax=40 else nmax=peaknr
if peaknr gt nmax then begin
    ind=reverse(sort(matrix[1,*]))
    ind=ind[0:nmax-1]
    matrix=matrix[*,ind]
    peaknr=40
endif

; Scale intensities
indf=where(finite(matrix[1,*]) eq 1)
m=abs(max(matrix[1,indf]))
matrix[1,*]=(matrix[1,*]/m*999)>0

; Sort d-spacings
matrix=matrix[*,reverse(sort(matrix[0,*]))]

if pcsiwin then begin
    for i=0l,peaknr-1 do printf,lun,matrix[0,i],matrix[1,i],1,format='(f15.4,I4,I2)'
endif else begin
    printf,lun,'$COMP:'
    printf,lun,(*ptr).name
    printf,lun,'$D:'
    printf,lun,peaknr
    printf,lun,matrix[0,*]
    printf,lun,'$I:'
    printf,lun,matrix[1,*]
    if bhkl then begin
        printf,lun,'$hkl:'
        printf,lun,round(matrix[2:4,*])
    endif
endelse

free_lun, lun
printw,list.OutID,path+' saved'

end;pro SavePKSTypeX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SendMainRefresh,top,data=data
widget_control,top,get_uvalue=top
widget_control,top,send_event={ID:top,TOP:top,HANDLER:top,ACTION:'refresh'+(keyword_set(data)?'data':'')}
end;pro SendMainRefresh
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ModelType10SetCode,ptr,ptrFit,codei,codeval,ibit=ibit,refresh=refresh,changed=changed

; ----Change code----
if not keyword_set(refresh) then begin
    case codei of
    0:    breturn=(codeval ge 7) and (ExtractbitCode(ptr,1) ge 1) ; asymmtype
    1:    breturn=(codeval ge 1) and (ExtractbitCode(ptr,0) ge 7) ; peaktype
    else: breturn=0b
    endcase
    if breturn then begin
        changed=1b ; for updating screen
        return
    endif

    if codei ge 4 and codeval eq 1 then $
        case ibit of
        2: begin
            SetbitCode,ptr,codei,0,ibit=4,changed=changed
            endcase
        3: begin
            SetbitCode,ptr,codei,0,ibit=4,changed=changed
            endcase
        4: begin
            SetbitCode,ptr,codei,0,ibit=2,changed=changed
            SetbitCode,ptr,codei,0,ibit=3,changed=changed
            endcase
        else:
        endcase

    SetbitCode,ptr,codei,codeval,ibit=ibit,changed=changed
endif else changed=0b

; ----Derived parameters----
PeakType=ExtractbitCode(ptr,0)
AsymType=ExtractbitCode(ptr,1)
ZeroType=ExtractbitCode(ptrFit,3)

; --paramind--

; paramind: indices used for u,I,W,n,R,A,All
; row0: in ARR1
; row1: in ARR2
ptr_free,(*ptr).paramind

(*ptr).paramind[0,0]=ptr_new([0])
(*ptr).paramind[1,0]=ptr_new([1])

(*ptr).paramind[0,1]=ptr_new([ZeroType,2])
(*ptr).paramind[1,1]=ptr_new([3])

codeproc0=[1b,1b,0b,0b,0b,0b]
codeproc6=codeproc0

case PeakType of
4:    begin
    codeproc0[[2,3]]=1b
    (*ptr).paramind[2,0]=ptr_new([2])
    (*ptr).paramind[3,0]=ptr_new([4])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,4])

    (*ptr).paramind[2,1]=ptr_new([4,5,6,7])
    (*ptr).paramind[3,1]=ptr_new([13,14])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,13,14])
    endcase
5:    begin
    codeproc0[[2]]=1b
    (*ptr).paramind[2,0]=ptr_new([2,3])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3])

    (*ptr).paramind[2,1]=ptr_new([4,5,6,7,8,9])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9])
    endcase
6:    begin
    codeproc0[[2,4]]=1b
    (*ptr).paramind[2,0]=ptr_new([2])
    (*ptr).paramind[4,0]=ptr_new([6])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,6])

    (*ptr).paramind[2,1]=ptr_new([4,5,6,7])
    (*ptr).paramind[4,1]=ptr_new([17,18,19])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,17,18,19])
    endcase
7:    begin
    codeproc0[[2,3,5]]=1b
    (*ptr).paramind[2,0]=ptr_new([2])
    (*ptr).paramind[3,0]=ptr_new([4,5])
    (*ptr).paramind[5,0]=ptr_new([8])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,4,5,8])

    (*ptr).paramind[2,1]=ptr_new([4,5,6,7])
    (*ptr).paramind[3,1]=ptr_new([13,14,15,16])
    (*ptr).paramind[5,1]=ptr_new([24,25,26])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,13,14,15,16,24,25,26])
    endcase
8:    begin
    codeproc0[[2,3,5]]=1b
    (*ptr).paramind[2,0]=ptr_new([2,3])
    (*ptr).paramind[3,0]=ptr_new([4,5])
    (*ptr).paramind[5,0]=ptr_new([8])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,5,8])

    (*ptr).paramind[2,1]=ptr_new([4,5,6,7,10,11,12])
    (*ptr).paramind[3,1]=ptr_new([13,14,15,16])
    (*ptr).paramind[5,1]=ptr_new([24,25,26])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,10,11,12,13,14,15,16,24,25,26])
    endcase
9:    begin
    codeproc0[[2,4,5]]=1b
    (*ptr).paramind[2,0]=ptr_new([2])
    (*ptr).paramind[4,0]=ptr_new([6,7])
    (*ptr).paramind[5,0]=ptr_new([8])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,6,7,8])

    (*ptr).paramind[2,1]=ptr_new([4,5,6,7])
    (*ptr).paramind[4,1]=ptr_new([20,21,22,23])
    (*ptr).paramind[5,1]=ptr_new([24,25,26])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,20,21,22,23,24,25,26])
    endcase
else: begin
    codeproc0[[2]]=1b
    (*ptr).paramind[2,0]=ptr_new([2])
    (*ptr).paramind[6,0]=ptr_new([0,1,2])

    (*ptr).paramind[2,1]=ptr_new([4,5,6,7])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7])
    endelse
endcase

case AsymType of
0:
1:     begin
    codeproc0[[5]]=1b
    if ptr_valid((*ptr).paramind[5,1]) then (*(*ptr).paramind[5,1])=[(*(*ptr).paramind[5,1]),27] $
    else (*ptr).paramind[5,1]=ptr_new([27])
    (*(*ptr).paramind[6,1])=[(*(*ptr).paramind[6,1]),27]
    endcase
2:     begin
    codeproc0[[5]]=1b
    if ptr_valid((*ptr).paramind[5,1]) then (*(*ptr).paramind[5,1])=[(*(*ptr).paramind[5,1]),28,29,30,31] $
    else (*ptr).paramind[5,1]=ptr_new([28,29,30,31])
    (*(*ptr).paramind[6,1])=[(*(*ptr).paramind[6,1]),28,29,30,31]
    endcase
endcase

; --codeproc--
code=ExtractbitCode(ptr,[4,5,6,7,8,9]) ; u, I, W, n, R, A
code1=((code and 1) eq 1)
code2=((code and 2) eq 2)
code4=((code and 4) eq 4)
code8=((code and 8) eq 8)
code16=((code and 16) eq 16)
(*ptr).codeproc=[[codeproc0],$ ; used?
    [codeproc0 and code1],$ ; global?
    [codeproc0 and (code2 or code4 or code8 or code16)],$ ; refined?
    [codeproc0 and code4],$ ; refine with constraints: a0>=CL
    [codeproc0 and code8],$ ; refine with constraints: a0<=CH
    [codeproc0 and code16],$ ; refine with constraints: a0+/-c
    [codeproc6],$ ; if global, function of peakparameters?
    [codeproc0 and code2 and (code4 or code8 or code16)]] ; use constraints as a rejection criterium

; --reset toggle--
if changed then begin
    if codei ge 4 then $
    if ibit eq 0 then begin ; bglobal changed
;        p=(*ptr).paramind[codei-4,1]
;        if ptr_valid(p) then (*ptr).togglefix[*p]=0
    endif
endif

end;pro ModelType10SetCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ModelType11SetCode,ptr,ptrFit,codei,codeval,ibit=ibit,refresh=refresh,changed=changed

; ----Change code----
if not keyword_set(refresh) then begin
    breturn=0b
    case codei of
    0:    breturn=(codeval ge 7) and (ExtractbitCode(ptr,1) ge 1) ; peaktype
    1:    breturn=(codeval ge 1) and (ExtractbitCode(ptr,0) ge 7) ; asymmtype
    2:    if (codeval ge 1) then SetbitCode,ptr,5,1,ibit=0,changed=changed ; LeBail
    4:    breturn=(codeval eq 0) and (ibit eq 0) ; u
    5:    breturn= (ExtractbitCode(ptr,2) eq 1) and (codeval eq 0) and (ibit eq 0) ; I  
    else:
    endcase
    if breturn then begin
        changed=1b ; for updating screen
        return
    endif

    if codei ge 4 and codeval eq 1 then $
        case ibit of
        2: begin
            SetbitCode,ptr,codei,0,ibit=4,changed=changed
            endcase
        3: begin
            SetbitCode,ptr,codei,0,ibit=4,changed=changed
            endcase
        4: begin
            SetbitCode,ptr,codei,0,ibit=2,changed=changed
            SetbitCode,ptr,codei,0,ibit=3,changed=changed
            endcase
        else:
        endcase

    SetbitCode,ptr,codei,codeval,ibit=ibit,changed=changed
endif else changed=0b

; ----Derived parameters----
PeakType=ExtractbitCode(ptr,0)
AsymType=ExtractbitCode(ptr,1)
LeBail=ExtractbitCode(ptr,2)
ZeroType=ExtractbitCode(ptrFit,3)

; --paramind--

; paramind: indices used for u,I,W,n,R,A,All
; row0: in ARR1
; row1: in ARR2
ptr_free,(*ptr).paramind

(*ptr).paramind[0,0]=ptr_new([0,1,2])
(*ptr).paramind[1,0]=ptr_new([3,4])

(*ptr).paramind[0,1]=ptr_new([ZeroType,2,3,4,5,6,7])
(*ptr).paramind[1,1]=ptr_new([8])

codeproc0=[1b,1b,0b,0b,0b,0b]
codeproc6=codeproc0
if LeBail then codeproc6[1]=2

case PeakType of
4:    begin
    codeproc0[[2,3]]=1b
    (*ptr).paramind[2,0]=ptr_new([5])
    (*ptr).paramind[3,0]=ptr_new([7])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,5,7])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12])
    (*ptr).paramind[3,1]=ptr_new([18,19])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,18,19])
    endcase
5:    begin
    codeproc0[[2]]=1b
    (*ptr).paramind[2,0]=ptr_new([5,6])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,5,6])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12,13,14])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,13,14])
    endcase
6:    begin
    codeproc0[[2,4]]=1b
    (*ptr).paramind[2,0]=ptr_new([5])
    (*ptr).paramind[4,0]=ptr_new([9])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,5,9])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12])
    (*ptr).paramind[4,1]=ptr_new([22,23,24])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,22,23,24])
    endcase
7:    begin
    codeproc0[[2,3,5]]=1b
    (*ptr).paramind[2,0]=ptr_new([5])
    (*ptr).paramind[3,0]=ptr_new([7,8])
    (*ptr).paramind[5,0]=ptr_new([11])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,5,7,8,11])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12])
    (*ptr).paramind[3,1]=ptr_new([18,19,20,21])
    (*ptr).paramind[5,1]=ptr_new([29,30,31])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,18,19,20,21,29,30,31])
    endcase
8:    begin
    codeproc0[[2,3,5]]=1b
    (*ptr).paramind[2,0]=ptr_new([5,6])
    (*ptr).paramind[3,0]=ptr_new([7,8])
    (*ptr).paramind[5,0]=ptr_new([11])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,5,6,7,8,11])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12,15,16,17])
    (*ptr).paramind[3,1]=ptr_new([18,19,20,21])
    (*ptr).paramind[5,1]=ptr_new([29,30,31])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,15,16,17,18,19,20,21,29,30,31])
    endcase
9:    begin
    codeproc0[[2,4,5]]=1b
    (*ptr).paramind[2,0]=ptr_new([5])
    (*ptr).paramind[4,0]=ptr_new([9,10])
    (*ptr).paramind[5,0]=ptr_new([11])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,5,9,10,11])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12])
    (*ptr).paramind[4,1]=ptr_new([25,26,27,28])
    (*ptr).paramind[5,1]=ptr_new([29,30,31])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,25,26,27,28,29,30,31])
    endcase
else: begin
    codeproc0[[2]]=1b
    (*ptr).paramind[2,0]=ptr_new([5])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,5])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12])
    endelse
endcase

case AsymType of
0:
1:     begin
    codeproc0[[5]]=1b
    if ptr_valid((*ptr).paramind[5,1]) then (*(*ptr).paramind[5,1])=[(*(*ptr).paramind[5,1]),32] $
    else (*ptr).paramind[5,1]=ptr_new([32])
    (*(*ptr).paramind[6,1])=[(*(*ptr).paramind[6,1]),32]
    endcase
2:     begin
    codeproc0[[5]]=1b
    if ptr_valid((*ptr).paramind[5,1]) then (*(*ptr).paramind[5,1])=[(*(*ptr).paramind[5,1]),33,34,35,36] $
    else (*ptr).paramind[5,1]=ptr_new([33,34,35,36])
    (*(*ptr).paramind[6,1])=[(*(*ptr).paramind[6,1]),33,34,35,36]
    endcase
endcase

; --codeproc--
code=ExtractbitCode(ptr,[4,5,6,7,8,9]) ; u, I, W, n, R, A
code1=((code and 1) eq 1)
code2=((code and 2) eq 2)
code4=((code and 4) eq 4)
code8=((code and 8) eq 8)
code16=((code and 16) eq 16)
(*ptr).codeproc=[[codeproc0],$ ; used?
    [codeproc0 and code1],$ ; global?
    [codeproc0 and (code2 or code4 or code8 or code16)],$ ; refined?
    [codeproc0 and code4],$ ; refine with constraints: a0>=CL
    [codeproc0 and code8],$ ; refine with constraints: a0<=CH
    [codeproc0 and code16],$ ; refine with constraints: a0+/-c
    [codeproc6],$ ; if global, function of peakparameters?
    [codeproc0 and code2 and (code4 or code8 or code16)]] ; use constraints as a rejection criterium

; --reset toggle--
if changed then begin
    if codei ge 4 then $
    if ibit eq 0 then begin ; bglobal changed
;        p=(*ptr).paramind[codei-4,1]
;        if ptr_valid(p) then (*ptr).togglefix[*p]=0
    endif
    if ptr_valid((*ptr).pstrinfo) then $
        (*ptr).togglefix[indCelParam(ptr)]or=celparamfix(*(*ptr).pstrinfo)
endif

end;pro ModelType11SetCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ModelType12SetCode,ptr,ptrFit,codei,codeval,ibit=ibit,refresh=refresh,changed=changed

; ----Change code----
if not keyword_set(refresh) then begin
    case codei of
    0:    breturn=(codeval ge 7) and (ExtractbitCode(ptr,1) ge 1) ; peaktype
    1:    breturn=(codeval ge 1) and (ExtractbitCode(ptr,0) ge 7) ; asymmtype
    4:    breturn=(codeval eq 0) and (ibit eq 0) ; u
    5:    breturn=(codeval eq 0) and (ibit eq 0) ; I
    else: breturn=0b
    endcase
    if breturn then begin
        changed=1b ; for updating screen
        return
    endif

    if codei ge 4 and codeval eq 1 then $
        case ibit of
        2: begin
            SetbitCode,ptr,codei,0,ibit=4,changed=changed
            endcase
        3: begin
            SetbitCode,ptr,codei,0,ibit=4,changed=changed
            endcase
        4: begin
            SetbitCode,ptr,codei,0,ibit=2,changed=changed
            SetbitCode,ptr,codei,0,ibit=3,changed=changed
            endcase
        else:
        endcase

    SetbitCode,ptr,codei,codeval,ibit=ibit,changed=changed
endif else changed=0b

; ----Derived parameters----
PeakType=ExtractbitCode(ptr,0)
AsymType=ExtractbitCode(ptr,1)
ZeroType=ExtractbitCode(ptrFit,3)

; --paramind--

; paramind: indices used for u,I,W,n,R,A,All
; row0: in ARR1
; row1: in ARR2
ptr_free,(*ptr).paramind

(*ptr).paramind[0,0]=ptr_new([0,1,2])
(*ptr).paramind[1,0]=ptr_new([3])

(*ptr).paramind[0,1]=ptr_new([ZeroType,2,3,4,5,6,7])
(*ptr).paramind[1,1]=ptr_new([8])

codeproc0=[1b,1b,0b,0b,0b,0b]
codeproc6=codeproc0

case PeakType of
4:    begin ; Pseudo-voigt
    codeproc0[[2,3]]=1b
    (*ptr).paramind[2,0]=ptr_new([4])
    (*ptr).paramind[3,0]=ptr_new([6])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,6])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12])
    (*ptr).paramind[3,1]=ptr_new([18,19])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,18,19])
    endcase
5:    begin ; TCHZ
    codeproc0[[2]]=1b
    (*ptr).paramind[2,0]=ptr_new([4,5])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,5])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12,13,14])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,13,14])
    endcase
6:    begin ; Pearson VII
    codeproc0[[2,4]]=1b
    (*ptr).paramind[2,0]=ptr_new([4])
    (*ptr).paramind[4,0]=ptr_new([8])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,8])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12])
    (*ptr).paramind[4,1]=ptr_new([22,23,24])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,22,23,24])
    endcase
7:    begin ; Split pseudo-voigt
    codeproc0[[2,3,5]]=1b
    (*ptr).paramind[2,0]=ptr_new([4])
    (*ptr).paramind[3,0]=ptr_new([6,7])
    (*ptr).paramind[5,0]=ptr_new([10])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,6,7,10])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12])
    (*ptr).paramind[3,1]=ptr_new([18,19,20,21])
    (*ptr).paramind[5,1]=ptr_new([29,30,31])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,18,19,20,21,29,30,31])
    endcase
8:    begin ; Mod Split Pseudo-Voigt
    codeproc0[[2,3,5]]=1b
    (*ptr).paramind[2,0]=ptr_new([4,5])
    (*ptr).paramind[3,0]=ptr_new([6,7])
    (*ptr).paramind[5,0]=ptr_new([10])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,5,6,7,10])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12,15,16,17])
    (*ptr).paramind[3,1]=ptr_new([18,19,20,21])
    (*ptr).paramind[5,1]=ptr_new([29,30,31])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,15,16,17,18,19,20,21,29,30,31])
    endcase
9:    begin ; Split Pearson VII
    codeproc0[[2,4,5]]=1b
    (*ptr).paramind[2,0]=ptr_new([4])
    (*ptr).paramind[4,0]=ptr_new([8,9])
    (*ptr).paramind[5,0]=ptr_new([10])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4,8,9,10])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12])
    (*ptr).paramind[4,1]=ptr_new([25,26,27,28])
    (*ptr).paramind[5,1]=ptr_new([29,30,31])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12,25,26,27,28,29,30,31])
    endcase
else: begin
    codeproc0[[2]]=1b
    (*ptr).paramind[2,0]=ptr_new([4])
    (*ptr).paramind[6,0]=ptr_new([0,1,2,3,4])

    (*ptr).paramind[2,1]=ptr_new([9,10,11,12])
    (*ptr).paramind[6,1]=ptr_new([ZeroType,2,3,4,5,6,7,8,9,10,11,12])
    endelse
endcase

case AsymType of
0:
1:     begin
    codeproc0[[5]]=1b
    if ptr_valid((*ptr).paramind[5,1]) then (*(*ptr).paramind[5,1])=[(*(*ptr).paramind[5,1]),32] $
    else (*ptr).paramind[5,1]=ptr_new([32])
    (*(*ptr).paramind[6,1])=[(*(*ptr).paramind[6,1]),32]
    endcase
2:     begin
    codeproc0[[5]]=1b
    if ptr_valid((*ptr).paramind[5,1]) then (*(*ptr).paramind[5,1])=[(*(*ptr).paramind[5,1]),33,34,35,36] $
    else (*ptr).paramind[5,1]=ptr_new([33,34,35,36])
    (*(*ptr).paramind[6,1])=[(*(*ptr).paramind[6,1]),33,34,35,36]
    endcase
endcase

; --codeproc--
code=ExtractbitCode(ptr,[4,5,6,7,8,9]) ; u, I, W, n, R, A
code1=((code and 1) eq 1)
code2=((code and 2) eq 2)
code4=((code and 4) eq 4)
code8=((code and 8) eq 8)
code16=((code and 16) eq 16)
(*ptr).codeproc=[[codeproc0],$ ; used?
    [codeproc0 and code1],$ ; global?
    [codeproc0 and (code2 or code4 or code8 or code16)],$ ; refined?
    [codeproc0 and code4],$ ; refine with constraints: a0>=CL
    [codeproc0 and code8],$ ; refine with constraints: a0<=CH
    [codeproc0 and code16],$ ; refine with constraints: a0+/-c
    [codeproc6],$ ; if global, function of peakparameters?
    [codeproc0 and code2 and (code4 or code8 or code16)]] ; use constraints as a rejection criterium

; --reset toggle--
if changed then begin
    if codei ge 4 then $
    if ibit eq 0 then begin ; bglobal changed
;        p=(*ptr).paramind[codei-4,1]
;        if ptr_valid(p) then (*ptr).togglefix[*p]=0
    endif
    if ptr_valid((*ptr).pstrinfo) then $
        (*ptr).togglefix[indCelParam(ptr)]or=celparamfix(*(*ptr).pstrinfo)
endif

end;pro ModelType12SetCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ModelTypeSetCode,ptr,codei,codeval,ibit=ibit,refresh=refresh,changed=changed
changed=0b

ptrFit=getptrFit(ptr)

case (*ptr).type of
10: ModelType10SetCode,ptr,ptrFit,codei,codeval,ibit=ibit,refresh=refresh,changed=changed
11: ModelType11SetCode,ptr,ptrFit,codei,codeval,ibit=ibit,refresh=refresh,changed=changed
12: ModelType12SetCode,ptr,ptrFit,codei,codeval,ibit=ibit,refresh=refresh,changed=changed
endcase
end;pro ModelTypeSetCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ModelTypeSetCodeAll,ptrFit
ptr=ptrFit
while ptr_valid((*ptr).next) do begin
    ptr=(*ptr).next
    ModelTypeSetCode,ptr,/refresh
endwhile
end;pro ModelTypeSetCodeAll
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MakeColorsTypeASU,ptr,colfor=colfor,colback=colback,constraints=constraints

COLuse=[255b,255b,255b] ; in use
COLnuse=[200b,200b,200b] ; not in use
COLfixed=[255b,0b,0b] ; fixed
COLrefined=[0b,0b,0b] ; not-fixed

nx=(*ptr).nasumax
ny=(*ptr).nasuatoms
colback=bytarr(3,nx,ny)
colback[0,*,*]=COLnuse[0]
colback[1,*,*]=COLnuse[1]
colback[2,*,*]=COLnuse[2]
colfor=bytarr(3,nx,ny)
colfor[0,*,*]=COLrefined[0]
colfor[1,*,*]=COLrefined[1]
colfor[2,*,*]=COLrefined[2]

; Used?
indasused=(*ptr).indasused
colback[0,indasused,*]=COLuse[0]
colback[1,indasused,*]=COLuse[1]
colback[2,indasused,*]=COLuse[2]

; Fixed?
if (*ptr).codeproc[1,2] eq 0 then begin
    colfor[0,indasused,*]=COLfixed[0]
    colfor[1,indasused,*]=COLfixed[1]
    colfor[2,indasused,*]=COLfixed[2]
endif else begin
    indasudep=(*ptr).indasudep

    ; Substition atoms
    indfixed=where((*(*ptr).asu.I).natoms eq -1,ct)
    if ct ne 0 then begin
        for i=0l,n_elements(indasudep)-1 do begin
            colfor[0,indasudep[i],indfixed]=COLfixed[0]
            colfor[1,indasudep[i],indfixed]=COLfixed[1]
            colfor[2,indasudep[i],indfixed]=COLfixed[2]
        endfor
    endif

    for j=0l,(*ptr).nasuatoms-1 do begin
        ; xyzfit
        if (*(*ptr).asu.I)[j].natoms ne 0 then begin
            RT=Symop64((*(*ptr).asu.I)[j].wyck.Rep)

            ; Same variables used
            Rotb=RT.Rot ne 0
            tmp=[array_equal(Rotb[*,0],Rotb[*,1]),$
                array_equal(Rotb[*,0],Rotb[*,2]),$
                array_equal(Rotb[*,1],Rotb[*,2])]
            dtmp=total(Rotb,1,/pres) eq 1
            if tmp[0] and dtmp[0] then begin
                ind=indasudep[0]+1
                colfor[0,ind,j]=COLfixed[0]
                colfor[1,ind,j]=COLfixed[1]
                colfor[2,ind,j]=COLfixed[2]
            endif
            if (tmp[1] and dtmp[1]) or (tmp[2] and dtmp[2]) then begin
                ind=indasudep[0]+2
                colfor[0,ind,j]=COLfixed[0]
                colfor[1,ind,j]=COLfixed[1]
                colfor[2,ind,j]=COLfixed[2]
            endif

            ; Fixed?
            ind=where((*(*ptr).asu.I)[j].wyck.togglefix,ct)
            if ct ne 0 then RT.Rot[ind,*]=0
            ind=where(round(total(abs(RT.rot),1)) eq 0,ct)
            if ct ne 0 then begin
                ind+=indasudep[0]
                colfor[0,ind,j]=COLfixed[0]
                colfor[1,ind,j]=COLfixed[1]
                colfor[2,ind,j]=COLfixed[2]
            endif
        endif

        ; SOF,Biso
        ind=where((*(*ptr).asu.I)[j].togglefix,ct)
        if ct ne 0 then begin
            ind+=indasused[3]
            colfor[0,ind,j]=COLfixed[0]
            colfor[1,ind,j]=COLfixed[1]
            colfor[2,ind,j]=COLfixed[2]
        endif
    endfor
endelse

end;pro MakeColorsTypeASU
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MakeColorsType,ptr,colfor=colfor,colback=colback,constraints=constraints

bconset=keyword_set(constraints)
nbconset=not bconset

COLuse=[255b,255b,255b] ; in use
COLnuse=[200b,200b,200b] ; not in use
COLfixed=[255b,0b,0b] ; fixed
COLrefined=[0b,0b,0b] ; not-fixed

n=(*ptr).nglobalmax+(*ptr).npeakmax
colback=bytarr(3,n,bconset?3:1)
colback[0,*,*]=COLnuse[0]
colback[1,*,*]=COLnuse[1]
colback[2,*,*]=COLnuse[2]
colfor=bytarr(3,n)
colfor[0,*]=COLrefined[0]
colfor[1,*]=COLrefined[1]
colfor[2,*]=COLrefined[2]

for i=0l,5 do begin
    if (*ptr).codeproc[i,0] then begin
    
        bglobal=(*ptr).codeproc[i,1]
        ind_=*(*ptr).paramind[i,bglobal]
    
        bfix = ~(*ptr).codeproc[i,2]
        if bglobal then begin
            toggle=(*ptr).togglefix[ind_]
            ind=ind_+(*ptr).npeakmax
                
            if bfix then indfix=ind else begin ; refined
                ; fix the ones that are toggled
                tmp=where(toggle,ct)
                bfix=ct ne 0
                if bfix then indfix=ind[tmp]
            endelse
            
            if (*ptr).codeproc[i,6] eq 2 and ~bconset then begin
                ind_=*(*ptr).paramind[i,0]
                ind2=where(~(*ptr).peakhelpparam[ind_],ct2)
                if ct2 ne 0 then begin
                    ; Given these peak parameters the same properties as global
                    ind=[ind,ind_[ind2]]
                    if bfix then indfix=[indfix,ind_[ind2]]
                endif
            endif
        endif else begin
            ind2=where(~(*ptr).peakhelpparam[ind_],ct2)
            if ct2 eq 0 then continue
            ind=ind_[ind2]
            indfix=ind
        endelse
    
        if bconset then begin ; constraints are displayed
    
            ; Low
            if (*ptr).codeproc[i,3] then begin
                colback[0,ind,0]=COLuse[0]
                colback[1,ind,0]=COLuse[1]
                colback[2,ind,0]=COLuse[2]
            endif
            ; High
            if (*ptr).codeproc[i,4] then begin
                colback[0,ind,1]=COLuse[0]
                colback[1,ind,1]=COLuse[1]
                colback[2,ind,1]=COLuse[2]
            endif
            ; PM
            if (*ptr).codeproc[i,5] then begin
                colback[0,ind,2]=COLuse[0]
                colback[1,ind,2]=COLuse[1]
                colback[2,ind,2]=COLuse[2]
            endif
    
        endif else begin ; parameters are displayed
            if bfix then begin; fixed?
                colfor[0,indfix]=COLfixed[0]
                colfor[1,indfix]=COLfixed[1]
                colfor[2,indfix]=COLfixed[2]
            endif
            colback[0,ind]=COLuse[0]
            colback[1,ind]=COLuse[1]
            colback[2,ind]=COLuse[2]
        endelse
    endif
endfor

end;pro MakeColorsType
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddModelTypeX,psource,pdest,ind

n=n_elements(ind)

icurrent=0
p=(*psource).peakparamIO
arrayI=(*p.I)[*,ind]
arrayO=(*(*p.O)[icurrent])[*,ind]
arraySD=(*(*p.SD)[icurrent])[*,ind]

if (*pdest).peaknr eq 0 then begin
    (*pdest).peakparamIO.I=PTR_NEW(arrayI)
    (*pdest).peakparamIO.constr=PTR_NEW(fltarr((*pdest).npeakmax,3))

    (*(*pdest).peakparamIO.O)[icurrent]=PTR_NEW(arrayO)
    (*(*pdest).peakparamIO.SD)[icurrent]=PTR_NEW(arraySD)

    arr=fltarr((*pdest).nglobalmax)

    (*pdest).globalparamIO.I=PTR_NEW(arr)
    (*pdest).globalparamIO.constr=PTR_NEW(fltarr((*pdest).nglobalmax,3))
    (*(*pdest).globalparamIO.O)[icurrent]=PTR_NEW(arr)
    (*(*pdest).globalparamIO.SD)[icurrent]=PTR_NEW(arr)
endif else begin
    p=(*pdest).peakparamIO
    *p.I=[[*p.I],[arrayI]]
    *(*(p.O))[icurrent]=[[*(*(p.O))[icurrent]],[arrayO]]
    *(*(p.SD))[icurrent]=[[*(*(p.SD))[icurrent]],[arraySD]]
endelse
(*pdest).peaknr+=n

end;pro AddModelTypeX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MovePeaksFromGroup,ptr,ind,top
if (*ptr).peaknr eq 0 then begin
    ind=-1
    return
endif

; Get all groups from the same category
ptrfit=getptrfit(ptr)
p=(*ptrfit).next
while ptr_valid(p) do begin
    if (*p).type eq (*ptr).type then $
    if n_elements(pkeep) eq 0 then pkeep=p else pkeep=[pkeep,p]
    p=(*p).next
endwhile

n=n_elements(pkeep)
if n le 1 then begin
    ind=-1
    return
endif

; User selects group where the peak(s) should be moved to
groups=strarr(n)
for i=0l,n-1 do groups[i]=(*pkeep[i]).name
i=promptnumber(groups,top,'Select destination group:',/choice)
pdest=pkeep[i]

if pdest eq ptr then begin
    ind=-1
    return
endif

; Add peaks to destination
AddModelTypeX,ptr,pdest,ind
end;pro MovePeaksFromGroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FilterDPeaksFromGroup,ptr,ind

; Select peaks with the same position
if (*ptr).type ne 10 then return
p=(*ptr).peakparamIO.I

ind=[-1]

n=(*ptr).peaknr
arr=indgen(n)
pind=*(*ptr).paramind[0,0]
tn=n_elements(pind)
for i=0l,n-1 do begin
    tmp=replicate(1b,n)
    for j=0l,tn-1 do tmp and= ((*p)[pind[j],*] eq (*p)[pind[j],i])
    tmp=where( tmp and (arr gt i),ct)
    if ct ne 0 then ind=[ind,tmp]
endfor

if n_elements(ind) gt 1 then ind=ind[1:*]

end; pro FilterDPeaksFromGroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FilterUPeaksFromGroupHandler,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

bdestroy=0b
case widget_info(ev.id,/type) of
    0:    bdestroy=1b
    1:    begin
        widget_control,ev.id,get_value=val
        if val ne 'Apply' then begin
            widget_control,ev.id,get_uvalue=uval
            ID=widget_info(ev.top,FIND_BY_UNAME=uval)
            widget_control,ID,sensitive=ev.select
        endif else bdestroy=1b
        endcase
    else:
endcase

if bdestroy then begin
    ID=widget_info(ev.top,FIND_BY_UNAME='buttonL')
    blower=widget_info(ID,/BUTTON_SET)
    ID=widget_info(ev.top,FIND_BY_UNAME='buttonU')
    bupper=widget_info(ID,/BUTTON_SET)
    if blower or bupper then begin
        ID=widget_info(ev.top,FIND_BY_UNAME='drop')
        i=widget_info(ID,/DROPLIST_SELECT)
        if i ne 0 then begin
            widget_control,ev.top,get_uvalue=list

            ind=(*((*list.ptr).paramind[6,0]))[i-1]

            arr=(*((*list.ptr).peakparamIO.I))[ind,*]
            b=bytarr((*list.ptr).peaknr)

            if blower then begin
                ID=widget_info(ev.top,FIND_BY_UNAME='textL')
                widget_control,ID,get_value=vlower
                vlower=float(vlower[0])
                ind=where(arr lt vlower,ct)
                if ct ne 0 then b[ind]=1
            endif

            if bupper then begin
                ID=widget_info(ev.top,FIND_BY_UNAME='textU')
                widget_control,ID,get_value=vupper
                vupper=float(vupper[0])
                ind=where(arr gt vupper,ct)
                if ct ne 0 then b[ind]=1
            endif

            ind=where(b,ct)
            if ct ne 0 then *(list.pind)=ind
        endif
    endif
    WIDGET_CONTROL, ev.top, /DESTROY
endif

end; pro FilterUPeaksFromGroupHandler
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FilterUPeaksFromGroup,ptr,ind,top

if (*ptr).peaknr eq 0 then return
pind=ptr_new([-1])

device,get_screen_size=screen
base=widget_base(/column,title='Set Conditions:',uvalue={ptr:ptr,top:top,pind:pind}, $
                xoffset=screen[0]/2.,yoffset=screen[1]/2.,/TLB_KILL_REQUEST_EVENTS)

ind=*((*ptr).paramind[6,0])

droplist=widget_droplist(base,value=['',(*ptr).peaklabels[ind]],uname='drop')
base1=widget_base(base,/row)
baset=widget_base(base1,/column,/nonexclusive)
    buttonL=widget_button(baset,value='Lower limit',uvalue='textL',uname='buttonL')
    buttonU=widget_button(baset,value='Upper limit',uvalue='textU',uname='buttonU')
baset=widget_base(base1,/column)
    textL=widget_text(baset,uname='textL',/editable,sensitive=1)
    textU=widget_text(baset,uname='textU',/editable,sensitive=1)
button=widget_button(baset,value='Apply')

WIDGET_CONTROL, base, /REALIZE
widget_control,buttonL,set_button=1
widget_control,buttonU,set_button=1

Xmanager,'FilterUPeaksFromGroup',base, event_handler='FilterUPeaksFromGroupHandler',GROUP_LEADER=top

ind=*pind
ptr_free,pind
if ind[0] eq -1 then delvar2,ind

end;pro FilterUPeaksFromGroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FilterVPeaksFromGroup,ptr,ind,top

widget_control,top,get_uvalue=top2
widget_control,top2,get_uvalue=list

ret=gfuncPeakparam(list,ptr)
I=ret[0,*]
I-=min(I)
if array_equal(I,0) then return
I*=100/max(I)

window,0
h=histogram(I)
plot,h,xrange=[0,100],yrange=[0,max(h)*1.1],psym=10,title='Peak intensity distribution'
tmp=I[reverse(sort(I))]
for j=0l,n_elements(I)-1 do xyouts,0.5+0.1*(j/10),0.8-0.02*(j mod 10),tmp[j],/normal
str=PromptNumber(stringr(1),top,'Enter peak intensity threshold (%):')
t=0>float(str)<100
if !d.window eq 0 then wdelete,!d.window

ind=where(I lt t)
end;pro FilterVPeaksFromGroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FilterPeaksFromGroup,ptr,ind

if (*ptr).type ne 10 then return
inddel=[-1]
p=(*ptr).peakparamIO.I

; filter out peaks with non-finite params
tmp=where(~finite(*p),ct)
if ct ne 0 then begin
    tmp/=(*ptr).npeakmax ; column indices
    inddel=[inddel,tmp]
endif

; filter out peaks with negative intensity or FWHM
pind=[*(*ptr).paramind[1,0],*(*ptr).paramind[2,0]]
tmp=where(total(((*p)[pind,*]) le 0,1) ne 0,ct)
if ct ne 0 then inddel=[inddel,tmp]

; finish
if n_elements(inddel) gt 1 then begin
    inddel=where(histogram(inddel[1:*],min=0,binsize=1) ne 0) ; remove all duplicates

    if n_elements(ind) ne 0 then begin
        ind=where(histogram([inddel,ind],min=0,binsize=1) eq 2,ct)  ; keep only when in both
        if ct eq 0 then delvar2,ind
    endif else ind=inddel

endif else delvar2,ind

end; pro FilterPeaksFromGroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeletePeaksFromGroup,ptr,ind,top=top,all=all,filter=filter,duplicates=duplicates,user=user,dlow=dlow,move=move

if (*ptr).peaknr eq 0 then return

if keyword_set(filter) then FilterPeaksFromGroup,ptr,ind else $
if keyword_set(duplicates) then FilterDPeaksFromGroup,ptr,ind else $
if keyword_set(user) then FilterUPeaksFromGroup,ptr,ind,top else $
if keyword_set(dlow) then FilterVPeaksFromGroup,ptr,ind,top else $
if keyword_set(move) then MovePeaksFromGroup,ptr,ind,top

n=n_elements(ind)
; ind not valid: n=0
if n ne 0 then if (where((ind lt 0) or (ind ge (*ptr).peaknr)))[0] ne -1 then n=0
if keyword_set(all) then n=(*ptr).peaknr
if n eq 0 then return

if n eq (*ptr).peaknr then begin ; delete all
    heap_free,(*ptr).peakparamIO
    heap_free,(*ptr).globalparamIO
    (*ptr).peakparamIO.O=ptr_new(PTRARR(1))
    (*ptr).peakparamIO.SD=ptr_new(PTRARR(1))
    (*ptr).globalparamIO.O=ptr_new(PTRARR(1))
    (*ptr).globalparamIO.SD=ptr_new(PTRARR(1))
    (*ptr).peaknr=0
    
    (*(*ptr).propI)[*]=0
    (*(*(*ptr).propO)[0])[*]=0
endif else begin
    nnew=(*ptr).peaknr-n
    icurrent=0
    p=(*ptr).peakparamIO.I
    (*p)=ShrinkArray(*p,ind,/ROW)
    p=(*(*ptr).peakparamIO.O)[icurrent]
    (*p)=ShrinkArray(*p,ind,/ROW)
    p=(*(*ptr).peakparamIO.SD)[icurrent]
    (*p)=ShrinkArray(*p,ind,/ROW)
    (*ptr).peaknr=nnew
endelse

end;pro DeletePeaksFromGroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FillModelType12,ptr,celparam,hkldata,FWHM

DeletePeaksFromGroup,ptr,/all

icurrent=0

if hkldata.n le 0 then return

arr=replicate(0.5,2,hkldata.n)
arr2=replicate(1,3,hkldata.n)
array=[hkldata.hkl,transpose(hkldata.mult),FWHM,FWHM,arr,arr2]
arr=fltarr((*ptr).npeakmax,hkldata.n)

(*ptr).peakparamIO.I=PTR_NEW(array)
(*ptr).peakparamIO.constr=PTR_NEW(fltarr((*ptr).npeakmax,3))

(*(*ptr).peakparamIO.O)[icurrent]=PTR_NEW(arr)
(*(*ptr).peakparamIO.SD)[icurrent]=PTR_NEW(arr)

arr=fltarr((*ptr).nglobalmax)
(*ptr).globalparamIO.I=PTR_NEW(arr)
(*ptr).globalparamIO.constr=PTR_NEW(fltarr((*ptr).nglobalmax,3))
(*(*ptr).globalparamIO.O)[icurrent]=PTR_NEW(arr)
(*(*ptr).globalparamIO.SD)[icurrent]=PTR_NEW(arr)

(*(*ptr).globalparamIO.I)[indCelParam(ptr)]=celparam
ind=[11,13,14,17,18,20,22,29]
;print,(*ptr).globallabels[ind]
(*(*ptr).globalparamIO.I)[ind]=[0.01,0.1,0.01,0.01,0.5,0.5,1,1]

(*ptr).peaknr=hkldata.n

end;pro FillModelType12
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FillModelType11,ptr,celparam,hkldata,area,FWHM

DeletePeaksFromGroup,ptr,/all

if hkldata.n le 0 then return

arr=replicate(0.5,2,hkldata.n)
arr2=replicate(1,3,hkldata.n)
array=[hkldata.hkl,transpose(hkldata.mult),area,FWHM,FWHM,arr,arr2]
arr=fltarr((*ptr).npeakmax,hkldata.n)

icurrent=0

(*ptr).peakparamIO.I=PTR_NEW(array)
(*ptr).peakparamIO.constr=PTR_NEW(fltarr((*ptr).npeakmax,3))

(*(*ptr).peakparamIO.O)[icurrent]=PTR_NEW(arr)
(*(*ptr).peakparamIO.SD)[icurrent]=PTR_NEW(arr)

arr=fltarr((*ptr).nglobalmax)
(*ptr).globalparamIO.I=PTR_NEW(arr)
(*ptr).globalparamIO.constr=PTR_NEW(fltarr((*ptr).nglobalmax,3))
(*(*ptr).globalparamIO.O)[icurrent]=PTR_NEW(arr)
(*(*ptr).globalparamIO.SD)[icurrent]=PTR_NEW(arr)

(*(*ptr).globalparamIO.I)[indCelParam(ptr)]=celparam
ind=[11,13,14,17,18,20,22,29]
;print,(*ptr).globallabels[ind]
(*(*ptr).globalparamIO.I)[ind]=[0.01,0.1,0.01,0.01,0.5,0.5,1,1]

(*ptr).peaknr=hkldata.n

end;pro FillModelType11
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FillModelType10,ptr,pos,area,FWHM

n=n_elements(pos)
arr=replicate(0.5,2,n)
arr2=replicate(1,3,n)
array=[pos,area,FWHM,FWHM,arr,arr2]
arr=fltarr((*ptr).npeakmax,n)

icurrent=0
if (*ptr).peaknr eq 0 then begin
    (*ptr).peakparamIO.I=PTR_NEW(array)
    (*ptr).peakparamIO.constr=PTR_NEW(fltarr((*ptr).npeakmax,3))

    (*(*ptr).peakparamIO.O)[icurrent]=PTR_NEW(arr)
    (*(*ptr).peakparamIO.SD)[icurrent]=PTR_NEW(arr)

    arr=fltarr((*ptr).nglobalmax)
    (*ptr).globalparamIO.I=PTR_NEW(arr)
    (*ptr).globalparamIO.constr=PTR_NEW(fltarr((*ptr).nglobalmax,3))
    (*(*ptr).globalparamIO.O)[icurrent]=PTR_NEW(arr)
    (*(*ptr).globalparamIO.SD)[icurrent]=PTR_NEW(arr)
    
    ind=[2,6,8,9,12,13,15,17,24]
;    print,(*ptr).globallabels[ind]
    (*(*ptr).globalparamIO.I)[ind]=[1,0.01,0.1,0.01,0.01,0.5,0.5,1,1]
    
endif else begin
    p=(*ptr).peakparamIO
    *p.I=[[*p.I],[array]]
    *(*(p.O))[icurrent]=[[*(*(p.O))[icurrent]],[arr]]
    *(*(p.SD))[icurrent]=[[*(*(p.SD))[icurrent]],[arr]]
endelse
(*ptr).peaknr+=n

end;pro FillModelType10
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ChangeSymmetry,ptr,sg,celparam,celdefadd,SmallDev=SmallDev

; no cell parameters given
if n_elements(celparam) eq 0 or (celdefadd ne 0) then begin
    celparamdef=celparamdefault(sg,celdefadd)
    if ptr_valid((*ptr).globalparamIO.I) then begin ; transform the old cel to the new one
        celparam=(*(*ptr).globalparamIO.I)[indCelParam(ptr)]
        ind=where(celparamfix(sg),ct)
        if ct ne 0 then celparam[ind]=celparamdef[ind]
    endif else begin ; there is no old cel
        celparam=float(celparamdef)
    endelse
endif

; connected cell parameters
celparam=celparamsetcon(celparam,sg)
bresetwyck=(*ptr).type eq 12 and ptr_valid((*ptr).pstrinfo)

; Recalculate ASU wyckoff
if bresetwyck then begin
    bresetwyck=(*(*ptr).pstrinfo).sghash ne sg.sghash and n_elements(SmallDev) ne 0
    if bresetwyck then ResetWyckpASU,sg,ptr,SmallDev
endif

ptr_free,(*ptr).pstrinfo
(*ptr).pstrinfo=ptr_new(sg)

return,1b
end;function ChangeSymmetry
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditTypeStructure,ev,ptr,sg,celparam,info,SmallDev=SmallDev,PDF=PDF,reset=reset
case (*ptr).type of
10:
else:    begin

    breset=keyword_set(reset)
    celdefadd=-1
    repeat begin
        celdefadd++

        ; change space group
        if breset then heap_free,(*ptr).globalparamIO.I
        tmp=ChangeSymmetry(ptr,sg,celparam,celdefadd,SmallDev=SmallDev)

        ; generate peaks
        hkldata=GenHKL(celparam,info.dlim,sg)

        ; reset: try again until (cell will be made larger)
        bHKL=(hkldata.n ne 0) or (~breset)
        
        if bHKL then if info.dlim[0] eq info.dlim[1] then begin
            tmp=dialog_message('No Bragg peaks within the current data range. Make it larger!',/error)
            bHKL=1
        endif
    endrep until bHKL

    ; peak intensities and FWHM
    if (hkldata.n ne 0) then begin
        ; Assume peak widths of 0.1
        FWHM=replicate(0.1,1,hkldata.n)

        ; peak intensities from PDF or data
        if ((*ptr).type eq 11) then begin

            if keyword_set(PDF) then begin
                ; Relative peak intensities from PDF-database
                ; Other peaks are set to 0
                ind=where(tag_names(PDF) eq 'PDFN',ct)
                if ct ne 0 then $
                if (PDF.PDFn ne 0) and (PDF.PDFni ne 0) then begin

                    area=fltarr(hkldata.n)
                    boolmatch=0b
                    boola=replicate(1b,hkldata.n)
                    boolb=replicate(1b,PDF.PDFn)
                    
                    ; Remove LGP factor from PDF intensities
                    PDFi=PDF.PDFi
                    tmp=float(promptnumber(stringr(PDF.lambda),ev.top,'Remove LGP-factor with lambda (0 when not applied)='))
                    if tmp ne 0 then begin
                        infoPDF={LGPtype:1,lambda:tmp,pol_Plindeg:0.5,pol_azimuth:45.}
                        PDFtt=BraggDtoX(PDF.PDFd,infoPDF.lambda,/onlyt,/angledeg)
                        PDFi/=LGPfactorpder(PDFtt,infoPDF)
                    endif

                    ; Find corresponding hkl
                    intersec=SetIntersection(hkldata.hkl,PDF.hkl,/nd,/secondind,count=ct)
                    if ct ne 0 then begin
                        area[intersec.a]=PDFi[intersec.b]
                        boola[intersec.a]=0
                        boolb[intersec.b]=0
                        boolmatch=1b
                    endif

                    ; Find equivalent hkl for others
                    inda=where(boola,cta)
                    indb=where(boolb,ctb)
                    if cta ne 0 and ctb ne 0 then begin
                        intersec=HKLequivalent(hkldata.hkl[*,inda],PDF.hkl[*,indb],sg.allops,count=ct)

                        if ct ne 0 then begin
                            area[inda[intersec.a]]=PDFi[indb[intersec.b]]
                            boola[inda[intersec.a]]=0
                            boolb[indb[intersec.b]]=0
                            boolmatch=1b
                        endif
                    endif
                    
                    ; Find equivalent d-spacings
                    inda=where(boola,cta)
                    indb=where(boolb,ctb)
                    if cta ne 0 and ctb ne 0 then begin
                        intersec=SetIntersection(hkldata.d[inda],PDF.PDFd[indb],/indices,/secondind,/uniq,count=ct,nafter=3)
        
                        if ct ne 0 then begin
                            area[inda[intersec.a]]=PDFi[indb[intersec.b]]
                            boolmatch=1b
                        endif
                    endif

                    ; Remove multiplicity from PDF intensities
                    if boolmatch then area/=hkldata.mult $
                    else delvar2,area

                endif
            endif
            
            if n_elements(area) eq 0 then begin
                ; Area from experimental pattern
                
                ; experimental peak height -> area, remove:
                ;    Area/Height ratio: Area = Height * [FWHM*sqrt(!pi/alog(2))/2]
                ;    multiplicity
                ;    LGP factor
                
                hkldatatt=BraggDtoX(hkldata.d,info.lambda,/onlyt,/angledeg)
                ym=hkldata.mult*LGPfactorpder(hkldatatt,info)*(20/sqrt(!pi/alog(2)))
                
                ; intensities from data
                ind=sort(info.xtt)
                x=info.xtt[ind]
                y=info.ytt[ind]
                Y2 = SPL_INIT(x,y)
                area = SPL_INTERP(x,y, Y2, hkldatatt)/ym
                
                if ~keyword_set(PDF) and (*ptr).peaknr ne 0 then begin
                    ; get intensities from previous structure (matching hkl)
                    ind=*((*ptr).paramind[0,0])
                    hkl=(*((*ptr).peakparamIO.I))[ind,*]
                    ind=*((*ptr).paramind[1,0])
                    I=(*((*ptr).peakparamIO.I))[ind[1],*]
    
                    intersec=SetIntersection(hkldata.hkl,fix(hkl),/nd,/secondind,/uniq,count=ct)
                    if ct ne 0 then area[intersec.a]=I[intersec.b]
                endif
            endif

            area=reform(area,1,hkldata.n)
        endif ; type11
    endif else begin; hkldata.n ne 0
      area=1.
      FWHM=0.1
    endelse
  
    if ((*ptr).type eq 11) then FillModelType11,ptr,celparam,hkldata,area,FWHM $
    else FillModelType12,ptr,celparam,hkldata,FWHM

    (*ptr).togglefix[indCelParam(ptr)] or= celparamfix(*(*ptr).pstrinfo)

    endcase
endcase
end;pro EditTypeStructure
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReCalculateHKL,ev,ptr
result=ExpData(ev.top,/dlim,/xd,/xtt,/y,/lambda,/LGPtype,/promptdlim)
;if ~TestTypeObsolete(ptr,result.dlim) then return

sg=*(*ptr).pstrinfo
EditTypeStructure,ev,ptr,sg,(*(*ptr).globalparamIO.I)[indCelParam(ptr)],result
BranchInfo,ev.top
end;pro ReCalculateHKL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PeakAutomatic,ptr,top
widget_control,top,get_uvalue=list

result=ExpData(list,/xtt,/y,/lambda,/LGPtype)
if n_elements(result.ytt) lt list.vin then return
out=estspe(result.ytt,result.xtt,list.r1,list.vin,list.win); area, position , sigma, mixing

s=size(out)
if s[1] ne 4 then begin
    printw,list.OutID,'No peaks found'
    return
endif

FillModelType10,ptr,out[1,*],out[0,*]/LGPfactorpder(out[1,*],result),out[2,*]*(2*SQRT(2*ALOG(2)))
end;pro PeakAutomatic
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PeakFromStruc,ptr,ev
widget_control,ev.top,get_uvalue=top
widget_control,top,get_uvalue=list
pddpath=list.pddpath
SmallDev=*(*(*list.ptrmodel).ptrfit).misc[3]
struc=LoadStructure(pddpath,'',nopos=(*ptr).type ne 12,SmallDev,error=error)
if error then return
list.pddpath=pddpath
widget_control,top,set_uvalue=list

; Change ASU
CopyASUstruc,struc.ASU,ptr

; Change Structure (SG and/or celparam)
result=ExpData(list,/xd,/xtt,/y,/dlim,/lambda,/LGPtype)

EditTypeStructure,ev,ptr,struc.SG,struc.celparam,result,PDF=struc

BranchInfo,ev.top
end;pro PeakFromStruc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro StrucFromPeak,ptr,ev

if ~ptr_valid((*ptr).pstrinfo) then return
widget_control,ev.top,get_uvalue=top
widget_control,top,get_uvalue=list
pddpath=list.pddpath

O=DIALOG_MESSAGE('Use refined values?',/question) eq 'Yes'

; List of peaks
matrix=gfuncPeakparam(list,ptr,O=O,/d)
ind=where(finite(matrix[0,*]) and finite(matrix[1,*]) and (matrix[0,*] ge 0) and (matrix[1,*] ge 0),ct)
if ct ne 0 then begin
    peaks=matrix[0:1,*]
    peaknr=ct
    
    ; Scale intensities
    indf=where(finite(peaks[1,*]) eq 1)
    m=abs(max(peaks[1,indf]))
    peaks[1,*]=(peaks[1,*]/m*999)>0
    
    ; Sort d-spacings
    peaks=peaks[*,reverse(sort(peaks[0,*]))]
endif else peaknr=0

if ~SaveStructure(pddpath,'',ptr,O=O,peaks=peaks,peaknr=peaknr) then return
list.pddpath=pddpath
widget_control,top,set_uvalue=list

end;pro StrucFromPeak
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PeakFromPDFHandler,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
    8:    begin
        widget_control,ev.top,get_uvalue=list2

        if ev.index ne 0 then i0=total((*list2.list.PDFn)[0:ev.index-1]) else i0=0
        i1=total((*list2.list.PDFn)[0:ev.index])-1

        PDFd=(*list2.list.PDFd)[i0:i1]
        xrange=(*list2.list.xvalt)[[list2.list.backranind[0],list2.list.backranind[1]],0]
        ind=where((PDFd ge (xrange[0]<xrange[1])) and (PDFd le xrange[0]>xrange[1]),n)
        if n ne 0 then begin
            ; Get peak positions
            PDFd=PDFd[ind]
            ind+=i0

            ; Assume peaks widths of 0.1
            FWHM=replicate(0.1,n)
            
            ; Peak relative intensities(areas) from PDF: PDFi
            ; Peaks heights from data: PDF2
            if ((*list2.list.PDFni)[ev.index] ne 0) then begin ; use PDFi
                PDFi=(*list2.list.PDFi)[ind]
                
                tmp=float(promptnumber('0',ev.top,'Remove LGP-factor with lambda (0 when not applied)='))
                if tmp ne 0 then begin
                    info={LGPtype:1,lambda:tmp,pol_Plindeg:0.5,pol_azimuth:45.}
                
                    tt=BraggDtoX(PDFd,info.lambda,/onlyt,/angledeg)
                    ym=LGPfactorpder(tt,info)
                    PDFi/=ym
                endif
                
            endif else begin ; use PDFi2
                ; Area = Height * [FWHM*sqrt(!pi/alog(2))/2]
                
                ; Convert height to area with removal of LGP
                info=ExpData(list2.list,/lambda,/LGPtype)
                tt=BraggDtoX(PDFd,info.lambda,/onlyt,/angledeg)
                ym=LGPfactorpder(tt,info)*(20/sqrt(!pi/alog(2)))
            
                PDFi=(*list2.list.PDFi2)[ind]/ym
            endelse
            
            ; Peak position in 2theta (deg)
            ind=reverse(sort(PDFd))
            PDFd=BraggDtoX(PDFd[ind],list2.list.lambda,/onlyt,/angledeg)
            PDFi=PDFi[ind]

            FillModelType10,list2.ptr,reform(PDFd,1,n),reform(PDFi,1,n),reform(FWHM,1,n)

            BranchInfo,list2.top
        endif
        endcase
endcase

WIDGET_CONTROL, ev.top, /DESTROY
end; pro PeakFromPDFHandler
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PeakFromPDF,ptr,ev

widget_control,ev.top,get_uvalue=top
widget_control,top,get_uvalue=list
if (list.nPDF eq 0) then return

device,get_screen_size=screen
base=widget_base(/column,title='Choose PDF:',uvalue={list:list,ptr:ptr,top:ev.top},$
                xoffset=screen[0]/2.,yoffset=screen[1]/2.)
droplist=widget_droplist(base,value=(*list.PDFname))

WIDGET_CONTROL, base, /REALIZE

Xmanager,'PeakFromPDF',base, event_handler='PeakFromPDFHandler',GROUP_LEADER=ev.top
end;pro PeakFromPDF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PeakManual_Event, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

; Check button event and return to fitmodel window when needed
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]
IF thisEvent NE 'DOWN' THEN RETURN
possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
bPress = possibleButtons[ev.press]

widget_control,ev.id,get_uvalue=top
ID=widget_info(top,FIND_BY_UNAME='topbranch')
widget_control,ID,get_uvalue=list2

if bPress eq 'RIGHT' or (~ptr_valid(list2.ptr)) then begin
    Widget_Control, ev.id ,draw_button=0,get_uvalue=top
    widget_control,top,sensitive=1
    return
endif

; Save first boarder
widget_control,ev.top,get_uvalue=list
list.xs=ev.x
list.ys=ev.y
widget_control,ev.top,set_uvalue=list

SetSysVar,list.sysvar
coords = Convert_Coord(ev.x,ev.y, /Device, /To_Data)
temp=min(abs((*list.xvalt)[*,list.xtype]-coords[0]),index)

; Plot first boarder
WSet, list.drawindex

x=(*list.xvalt)[index,list.xtype]
y=(*list.spet)[index,list.xtype]

if list.ypower ne 0 and list.ypower ne 1 then y^=list.ypower

plots,[x,x],[list.ylog,y],$
       linestyle=line,NOCLIP = 0
oplot,[x],[y],psym=2,NOCLIP = 0

h=!D.x_size
v=!D.y_size
WSet, list.pixindex
Device, Copy = [0,0, h,v, 0,0,list.drawindex]

; Go to next event handler
Widget_Control, ev.id, Event_Pro='PeakManualMove_Event', $
    draw_motion=1
end;pro PeakManual_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PeakManualMove_Event,ev

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
       Event_Pro='PeakManual_Event'
    x=(*list.xvalt)[*,1] ; two-theta
    y=(*list.spet)[*,1]
    if list.xtypeorig eq 0 then begin
        ny=n_elements(y)
        index0=ny-1-index0
        index1=ny-1-index1
    endif

    peakbegin=index0<index1
    peakend=index0>index1
    if x[0] gt x[1] then begin
        x=reverse(x)
        y=reverse(y)
    endif
    peakval=estpeak(x,y,peakbegin,peakend,debug=0,drawi=0); pos, index, area, sigma, height

    widget_control,ev.id,get_uvalue=top
    ID=widget_info(top,FIND_BY_UNAME='topbranch')
    widget_control,ID,get_uvalue=list2
    
    info=ExpData(list,/lambda,/LGPtype)

    FillModelType10,list2.ptr,peakval[0],peakval[2]/LGPfactorpder(peakval[0],info),peakval[3]*(2*SQRT(2*ALOG(2)))
    BranchInfo,top
    RefreshDisplayCHI,list
    return
endif

WSet, list.drawindex

x=(*list.xvalt)[index1,list.xtype]
y=(*list.spet)[index1,list.xtype]
if list.ypower ne 0 and list.ypower ne 1 then y^=list.ypower
plots,[x,x],[list.ylog,y],linestyle=line,NOCLIP = 0
oplot,[x],[y],psym=2,NOCLIP = 0

end;pro PeakManualMove_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_PromptModelEntreeName,ID

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
end;CleanUp_PromptModelEntreeName
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PromptModelEntreeName_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
1 : begin ;button event
    widget_control,ev.id,get_value=val
    case val of
    'Default':begin
            ID=widget_info(ev.top,FIND_BY_UNAME='text')
            widget_control,ID,set_value=''
            return
            endcase
    else:
    endcase
    endcase
else:
endcase

; Read from text widget
ID=widget_info(ev.top,FIND_BY_UNAME='text')
widget_control,ID,get_value=str
; Save dynamically
widget_control,ev.top,get_uvalue=list
(*list.filter)=str[0]
; Return
widget_control,ev.top,/destroy
end;pro PromptModelEntreeName_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PromptModelEntreeName,in,top,label,xsize=xsize

;----Prompt for filter----
filter=PTR_NEW(in)
if WIDGET_INFO(top,/VALID_ID) eq 1 then begin
    WIDGET_CONTROL, top, sensitive=0
    GROUP_LEADER=top
    modal=1
endif else begin
    modal=0 ; no modal when no group leader
endelse

device,get_screen_size=screen
base=widget_base(title='',/row,uvalue={top:top,filter:filter},GROUP_LEADER=GROUP_LEADER,$
            modal=modal,xoffset=screen[0]/2.,yoffset=screen[1]/2.)
label=widget_label(base,value=label)
text=widget_text(base,/editable,value=(*filter),uname='text')
b=widget_button(base,value='Default')
b=widget_button(base,value='OK')
widget_control,base,/realize
widget_control,b,/INPUT_FOCUS
Xmanager,'PromptModelEntreeName',base, event_handler='PromptModelEntreeName_Event',$
    cleanup='CleanUp_PromptModelEntreeName',GROUP_LEADER=GROUP_LEADER

out=(*filter)
PTR_FREE,filter
return,out
end;function PromptModelEntreeName
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ChangeModelEntreeName,ptr,id,top
str=PromptModelEntreeName((*ptr).name,top,'Enter new name:')
if str eq '' then begin
    case (*ptr).type of
    10:    str='PD'
    11:    str='Pawley'
    12:    str='Rietveld'
    endcase
    str=str+stringr((*ptr).group)
endif

(*ptr).name=str
widget_control,id,set_value=str

end;pro ChangeModelEntreeName
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetTableID,ev

ID=widget_info(ev.top,FIND_BY_UNAME='groupinfochildtab')
if not widget_info(ID,/VALID_ID) then return,0
tabsel=widget_info(ID,/TAB_CURRENT)

ID=widget_info(ID,/child)
for i=0l,tabsel-1 do ID = Widget_Info(ID, /Sibling)
if not widget_info(ID,/VALID_ID) then return,0

return,widget_info(ID,FIND_BY_UNAME='Peak parameters:')

end;function GetTableID,ev
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditTypeTableASU,type,x,y,ptr,value,resval=resval,obsolete=obsolete

icurrent=0

case type of
4:p=(*ptr).asu.I
5:return
6:p=(*(*ptr).asu.O)[icurrent]
7:return
endcase

; ----Get precision----
tmp=getptrFit(ptr)
SmallDev=*((*tmp).misc[3])

EditASUpos,p,(*ptr).preverse_indices,(*ptr).pstrinfo,x,y,value,SmallDev,resval=resval,obsolete=obsolete

if type eq 6 then begin
    pSB=(*(*ptr).asu.SD)[icurrent]
    (*pSB)=(*p)
    (*pSB).xyz=0
    (*pSB).SOF=0
    (*pSB).B=0
endif

end;pro EditTypeTableASU
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditTypeTable,tabletype,type,x,y,ptr,value,resval=resval,obsolete=obsolete,dlim=dlim

icurrent=0
obsolete=0b

if (tabletype eq 0) then begin
    case type of
    0:p=(*ptr).peakparamIO.I
    1:p=(*ptr).peakparamIO.constr
    2:p=(*(*ptr).peakparamIO.O)[icurrent]
    3:p=(*(*ptr).peakparamIO.SD)[icurrent]
    endcase

    i=(*(*ptr).paramind[6,0])[x[0]:x[1]]
endif else begin
    case type of
    0:p=(*ptr).globalparamIO.I
    1:p=(*ptr).globalparamIO.constr
    2:p=(*(*ptr).globalparamIO.O)[icurrent]
    3:p=(*(*ptr).globalparamIO.SD)[icurrent]
    endcase

    i=(*(*ptr).paramind[tabletype-1,1])[x[0]:x[1]]

    ; check for cell parameter connections
    if ((*ptr).type ne 10) then begin
        indCel=indCelParam(ptr)
        nindCel=n_elements(indCel)
        imin=indCel[0]
        imax=indCel[nindCel-1]
        ind=where(i ge imin and i le imax,ct)
        if ct ne 0 then begin
            indu=(*(*ptr).paramind[0,1])
            nindu=n_elements(indu)
            add=replicate(0,nindu-nindCel)
            
            ; cell parameter input
            valueorg=(*p)[indu,y[0]]
            resval=valueorg
            resval[x[0],0]=value

            ; change only values that are allowed to change
            tmp=[add,celparamfix(*(*ptr).pstrinfo)] ; add elements for "zero"
            ind=where(tmp eq 1,ct)
            if ct ne 0 then resval[ind]=valueorg[ind]

            ; connect values
            connect=celparamconnect(*(*ptr).pstrinfo)
            connect=[add,connect+1]; add elements for "zero"
            ind=where(histogram(connect,REVERSE_INDICES=R) gt 1,ct)
            for j=0l,ct-1 do begin ; change connected values + return all changed values
                tmp=R[R[ind[j]] : R[ind[j]+1]-1]
                resval[tmp]=resval[tmp[0]]
            endfor

            (*p)[indu,y[0]]=resval
            resval=[[resval],[indgen(nindu)],[replicate(y[0],nindu)]]
            if type eq 0 then obsolete=TestTypeObsolete(ptr,dlim)
            
            ; Calculate other parameters based on the unit cell parameter
            ;CapobiancoAdaptASU,ptr,type
            
            return
        endif
    endif
endelse

(*p)[i,y[0]:y[1]]=value
;if n_elements(dlim) ne 0 then obsolete=TestTypeObsolete(ptr,dlim)

end;pro EditTypeTable
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditTypeTableev,ev,tabletype,type,x,y,ptr,array,no_obsolete=no_obsolete

if type le 3 then begin
    result=ExpData(ev.top,/dlim,/xd,/xtt,/y,/lambda,/LGPtype)
    EditTypeTable,tabletype,type,x,y,ptr,array[x[0]:x[1],y[0]:y[1]],resval=resval,dlim=result.dlim,obsolete=obsolete
    if n_elements(resval) ne 0 then begin
        array[resval[*,1],resval[*,2]]=stringr(resval[*,0])
        widget_control,ev.id,set_value=array
    endif
    if obsolete and ~keyword_set(no_obsolete) then begin
        if dialog_message('Keep the current Bragg peaks?',/question) eq 'No' then begin
            sg=*(*ptr).pstrinfo
            EditTypeStructure,ev,ptr,sg,(*(*ptr).globalparamIO.I)[indCelParam(ptr)],result
            BranchInfo,ev.top
        endif
    endif
    return
endif

if type le 7 then begin
    nx=x[1]-x[0]+1
    ny=y[1]-y[0]+1
    xx=x[0]+indgen(nx)
    yy=y[0]+indgen(ny)
    obsolete=0b
    for i=0l,nx-1 do for j=0l,ny-1 do begin
        EditTypeTableASU,type,xx[i],yy[j],ptr,array[xx[i],yy[j]],resval=resval,obsolete=tmp
        obsolete or=tmp
        if n_elements(resval) ne 0 then begin
            array[resval[*,1],resval[*,2]]=stringr(resval[*,0])
            widget_control,ev.id,set_value=array
            delvar2,resval
        endif
    endfor
    if obsolete then BranchInfo,ev.top,/update
    return
endif

if type eq 8 then begin
    ny=y[1]-y[0]+1
    if ny eq 1 then if y[0] mod 2 eq 0 then begin
        if array[y[0]] eq 0 then array[y+1]=0 $
        else array[y[0]+1]=defaulterror(array[y[0]])
        widget_control,ev.id,set_value=array
    endif
    
    (*ptr).propfix=array
endif

end;pro EditTypeTableev
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro typeTable,base,ptr,refined=refined,constraints=constraints,sd=sd,update=update

if (*ptr).peaknr eq 0 then return
if keyword_set(sd) then type=3 else $
if keyword_set(refined) then type=2 else $
if keyword_set(constraints) then type=1 else type=0
noupdate=not keyword_set(update)

if noupdate then begin
    ; ----Delete all children----
    DeleteWBaseChilds,base

    ; ----Make table labels----
    label=widget_label(base,value='Cell color: Gray = not used, Red = fixed',/ALIGN_LEFT)
    label=widget_label(base,value=' ',/ALIGN_LEFT)
endif

; ----Make/edit tables----
MakeColorsType,ptr,colback=colback,colfor=colfor,constraints=constraints

editable=type ne 3
;ALL_EVENTS=type ne 1
ALL_EVENTS=0
icurrent=0

    ind=*((*ptr).paramind[6,0])
    n=n_elements(ind)
    case type of
    0:     value=(*((*ptr).peakparamIO.I))[ind,*]
    1:     begin
        value=[(*((*ptr).peakparamIO.constr))[ind,*]]
        ROW_LABELS=['>= C0','<= C1',msymbols('plus_minus')+'C']
        endcase
    2:     value=(*((*(*ptr).peakparamIO.O)[icurrent]))[ind,*]
    3:    value=(*((*(*ptr).peakparamIO.SD)[icurrent]))[ind,*]
    endcase
    BACKGROUND_COLOR=colback[*,ind,*]
    FOREGROUND_COLOR=colfor[*,ind,*]
    if noupdate then begin
        label=widget_label(base,/ALIGN_LEFT,value='Peak parameters:')
        table = WIDGET_TABLE(base, COLUMN_LABELS  = (*ptr).peaklabels[ind], EDITABLE=editable,CONTEXT_EVENTS=editable,$
         value=value,uname='Peak parameters:',uvalue='Peak parameters:',X_SCROLL_SIZE =n,Y_SCROLL_SIZE=3>((*ptr).peaknr<6),$
         BACKGROUND_COLOR=BACKGROUND_COLOR,FOREGROUND_COLOR=FOREGROUND_COLOR,ROW_LABELS=ROW_LABELS,ALL_EVENTS=ALL_EVENTS)
         widget_control,table,set_table_select=[-1,-1,-1,-1]
    endif else begin
        ID=widget_info(base,FIND_BY_UNAME='Peak parameters:')
        if widget_info(ID,/VALID_ID) then widget_control,ID,set_value=value,$
         BACKGROUND_COLOR=BACKGROUND_COLOR,FOREGROUND_COLOR=FOREGROUND_COLOR
    endelse


str=['Global position parameters:','Global intensity parameters:','Global FWHM parameters:','Global mixing parameters:','Global decay parameters:','Global asymmetry parameters:']

for i=0l,5 do begin
    if (*ptr).codeproc[i,0] then begin
        ind=(*((*ptr).paramind[i,1]))
        n=n_elements(ind)
        case type of
        0:     value=(*((*ptr).globalparamIO.I))[ind]
        1:     begin
            value=[(*((*ptr).globalparamIO.constr))[ind,*]]
            ROW_LABELS=['>= C0','<= C1',msymbols('plus_minus')+'C']
            endcase
        2:     value=(*((*(*ptr).globalparamIO.O)[icurrent]))[ind]
        3:    value=(*((*(*ptr).globalparamIO.SD)[icurrent]))[ind]
        endcase
        BACKGROUND_COLOR=colback[*,ind+(*ptr).npeakmax,*]
        FOREGROUND_COLOR=colfor[*,ind+(*ptr).npeakmax,*]
        if noupdate then begin
            label=widget_label(base,/ALIGN_LEFT,value=str[i])
            table = WIDGET_TABLE(base,ROW_LABELS=ROW_LABELS, COLUMN_LABELS  = (*ptr).globallabels[ind],$
             value=value,uname=str[i],uvalue=str[i],X_SCROLL_SIZE =n,$
             BACKGROUND_COLOR=BACKGROUND_COLOR,FOREGROUND_COLOR=FOREGROUND_COLOR, EDITABLE=editable,CONTEXT_EVENTS=editable)
            widget_control,table,set_table_select=[-1,-1,-1,-1]
        endif else begin
            ID=widget_info(base,FIND_BY_UNAME=str[i])
            if widget_info(ID,/VALID_ID) then widget_control,ID,set_value=value,$
             BACKGROUND_COLOR=BACKGROUND_COLOR,FOREGROUND_COLOR=FOREGROUND_COLOR
        endelse

    endif
endfor

end;pro typeTable
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro typeTableasu,base,ptr,refined=refined,constraints=constraints,sd=sd,update=update

if (*ptr).peaknr eq 0 then return
if keyword_set(sd) then type=3 else $
if keyword_set(refined) then type=2 else $
if keyword_set(constraints) then type=1 else type=0
noupdate=not keyword_set(update)

if noupdate then begin
    ; ----Delete all children----
    DeleteWBaseChilds,base
endif

; ----Make/edit tables----
if (*ptr).nasupos eq 0 then return

editable=type ne 3
;ALL_EVENTS=type ne 1
ALL_EVENTS=0
icurrent=0

n=(*ptr).nasumax
case type of
0:     begin
    value=ConvertASUtotable(*((*ptr).asu.I))
    endcase
1:     begin
    if noupdate then begin
        label=widget_label(base,value='Following physical constraints are always used:',/ALIGN_LEFT)
        label=widget_label(base,value='1. 0 <= xyz < 1',/ALIGN_LEFT)
        label=widget_label(base,value='2. 0 <= SOF <= 1',/ALIGN_LEFT)
        label=widget_label(base,value='3. som(SOF) for each ASU position <= 1',/ALIGN_LEFT)
        label=widget_label(base,value='4. 0 <= Biso',/ALIGN_LEFT)
    endif
    return
    endcase
2:     begin
    value=ConvertASUtotable(*((*(*ptr).asu.O)[icurrent]))
    endcase
3:    begin
    value=ConvertASUtotable(*((*(*ptr).asu.SD)[icurrent]))
    endcase
endcase

MakeColorsTypeASU,ptr,colback=BACKGROUND_COLOR,colfor=FOREGROUND_COLOR,constraints=constraints

case type of
0:stoich=CalcStoichiometry(ptr)
2:stoich=CalcStoichiometry(ptr,/O)
else:
endcase

if noupdate then begin
    labels=(*ptr).asulabels
    
    ; ----Make table----
    label=widget_label(base,value='Cell color: Gray = no parameter, Red = fixed',/ALIGN_LEFT)
    label=widget_label(base,value=' ',/ALIGN_LEFT)
    label=widget_label(base,/ALIGN_LEFT,value='ASU parameters:')
    
    table = WIDGET_TABLE(base, COLUMN_LABELS  = labels, EDITABLE=editable,CONTEXT_EVENTS=editable,$
     value=value,uname='ASU parameters:',uvalue='ASU parameters:',X_SCROLL_SIZE =n,Y_SCROLL_SIZE=((*ptr).nasuatoms<6),$
     BACKGROUND_COLOR=BACKGROUND_COLOR,FOREGROUND_COLOR=FOREGROUND_COLOR,ROW_LABELS=ROW_LABELS,ALL_EVENTS=ALL_EVENTS)
     widget_control,table,set_table_select=[-1,-1,-1,-1]

    ; ----Cell info----
    base=widget_base(base,/row)
    if n_elements(stoich) ne 0 then begin
        baset=widget_base(base,/column)
        label=widget_label(baset,value='Unit cell properties: ',/ALIGN_LEFT)
        label=widget_label(baset,value=' ',/ALIGN_LEFT)
        label=widget_text(baset,value=stoich,/ALIGN_LEFT,uname='stoichiometry',ysize=n_elements(stoich)<15,/scroll)
        label=widget_label(baset,value=' ',/ALIGN_LEFT)
    endif
    
    
endif else begin
    ID=widget_info(base,FIND_BY_UNAME='stoichiometry')
    if widget_info(ID,/VALID_ID) then widget_control,ID,set_value=stoich,ysize=n_elements(stoich)<15
    
    ID=widget_info(base,FIND_BY_UNAME='ASU parameters:')
    if widget_info(ID,/VALID_ID) then widget_control,ID,set_value=value,$
     BACKGROUND_COLOR=BACKGROUND_COLOR,FOREGROUND_COLOR=FOREGROUND_COLOR
endelse

; ----Plot----
case type of
0: PlotUC,ptr,base,update=update
2: PlotUC,ptr,base,update=update,/O
else:
endcase

end;pro typeTableasu
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ExcludedBMP

;i=READ_IMAGE('C:\work\C\excluded24.bmp') ;=> bytarr(3,16,16)
;i2=TRANSPOSE(i,[1,2,0]) ;=> bytarr(16,16,3)
;print,i2,format='(15(I3,"b,"),"$")' ;=> image
;print,where(i2 eq 0),format='(15(I3,","),"$")' ;=> special case

image=replicate(255b,16,16,3)
ind=[277,278,279,280,281,282,292,293,294,295,296,297,298,299,307,$
308,309,310,311,312,313,314,315,316,322,323,324,325,326,327,$
328,329,330,331,332,333,337,338,339,340,341,342,343,344,345,$
346,347,348,349,350,353,354,355,356,357,358,359,360,361,362,$
363,364,365,366,369,382,385,398,401,402,403,404,405,406,407,$
408,409,410,411,412,413,414,417,418,419,420,421,422,423,424,$
425,426,427,428,429,430,434,435,436,437,438,439,440,441,442,$
443,444,445,451,452,453,454,455,456,457,458,459,460,468,469,$
470,471,472,473,474,475,485,486,487,488,489,490,533,534,535,$
536,537,538,548,549,550,551,552,553,554,555,563,564,565,566,$
567,568,569,570,571,572,578,579,580,581,582,583,584,585,586,$
587,588,589,593,594,595,596,597,598,599,600,601,602,603,604,$
605,606,609,610,611,612,613,614,615,616,617,618,619,620,621,$
622,625,638,641,654,657,658,659,660,661,662,663,664,665,666,$
667,668,669,670,673,674,675,676,677,678,679,680,681,682,683,$
684,685,686,690,691,692,693,694,695,696,697,698,699,700,701,$
707,708,709,710,711,712,713,714,715,716,724,725,726,727,728,$
729,730,731,741,742,743,744,745,746]
image[ind]=0b

image=reform(image,16,16,3,/OVERWRITE)
return,image
end;function ExcludedBMP
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshMainTab,top

ID=widget_info(top,FIND_BY_UNAME='topbranch')
widget_control,ID,get_uvalue=list2

ID=widget_info(top,FIND_BY_UNAME='dropprofiletype')
widget_control,ID,SET_DROPLIST_SELECT=ExtractbitCode(list2.ptr,0)

;ID=widget_info(top,FIND_BY_UNAME='dropasymtype')
;widget_control,ID,SET_DROPLIST_SELECT=ExtractbitCode(list2.ptr,1)

ID=widget_info(top,FIND_BY_UNAME='dropItype')
if widget_info(ID,/valid) then widget_control,ID,SET_DROPLIST_SELECT=ExtractbitCode(list2.ptr,2)

ID=widget_info(top,FIND_BY_UNAME='dropfittype')
sel=widget_info(ID,/DROPLIST_SELECT)
widget_control,ID,get_uvalue=struct

code=ExtractbitCode(list2.ptr,(struct.i)[sel])
base=widget_info(id,/parent)
ID=widget_info(base,FIND_BY_UNAME='Global')
widget_control,ID,set_button=(code and 1) eq 1
ID=widget_info(base,FIND_BY_UNAME='Refined')
widget_control,ID,set_button=(code and 2) eq 2
ID=widget_info(base,FIND_BY_UNAME='Refined (constraints >= Clow)')
widget_control,ID,set_button=(code and 4) eq 4
ID=widget_info(base,FIND_BY_UNAME='Refined (constraints <= Chigh)')
widget_control,ID,set_button=(code and 8) eq 8
ID=widget_info(base,FIND_BY_UNAME='Refined (constraints +/-C)')
widget_control,ID,set_button=(code and 16) eq 16
end;pro RefreshMainTab
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MainCheckEvent,ev
widget_control,ev.id,get_uvalue=ibit

ID=widget_info(ev.top,FIND_BY_UNAME='dropfittype')
sel=widget_info(ID,/DROPLIST_SELECT)
widget_control,ID,get_uvalue=struct
ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
widget_control,ID,get_uvalue=list2

ModelTypeSetCode,list2.ptr,(struct.i)[sel],ev.select,ibit=ibit,changed=changed
if changed then BranchInfo,ev.top,/update

end;pro MainCheckEvent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MainDropEvent,ev,i
ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
widget_control,ID,get_uvalue=list2

ModelTypeSetCode,list2.ptr,i,ev.index,changed=changed
if changed then BranchInfo,ev.top
end;pro MainDropEvent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro InitGlobal,ev,type

; Get model info
widget_control,widget_info(ev.id,/parent),get_uvalue=struct
widget_control,widget_info(struct.ID,/parent),get_uvalue=uval2
bout =uval2 eq 'RefinedTab'

ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
widget_control,ID,get_uvalue=list2
list=ExpData(ev.top,/lambda,/xd,/y,/xtt,/LGPtype,/dist,/zerotype)

ptr=list2.ptr

widget_control,struct.ID,get_value=array

ind=*(*ptr).paramind[type,1]
icurrent=0
if bout then begin
    p=(*(*ptr).peakparamIO.O)[icurrent]
    pg=(*(*ptr).globalparamIO.O)[icurrent]
endif else begin
    p=(*ptr).peakparamIO.I
    pg=(*ptr).globalparamIO.I
endelse

; Initial peak positions and refined FWHM, nmix,...
case type of
0:    array=[0.,1.]
1:    begin
    widget_control,ev.top,get_uvalue=top
    widget_control,top,get_uvalue=list3
    array=[InitScaling(ptr,list3,list.xtt,list.ytt,O=bout)]
    endcase
else: begin
    case (*ptr).type of
    10:    tt=reform((*p)[*(*ptr).paramind[0,0],*])
    else: tt=ttfromhkl((*p)[*(*ptr).paramind[0,0],*],(*pg)[indCelParam(ptr)],list.lambda,$
                        list.dist,(*pg)[indZeroParam(ptr)],list.zerotype)
    endcase
    tmp=(*p)[*(*ptr).paramind[type,0],*]
    case type of
    2:    fitGlobalFWHM,p,array,tt,tmp
    3:    fitGlobalnmix,p,array,tt,tmp
    4:    fitGlobalR,p,array,tt,tmp
    5:    fitGlobalAsym,p,array,tt,tmp
    endcase
    endcase
endcase

; update table and model
(*pg)[ind]=array
widget_control,struct.ID,set_value=array
end;pro InitGlobal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Type12nASUParam,ptr,stop=stop

stop=0b
if (*ptr).type ne 12 then return,0
if (*ptr).nasuatoms eq 0 then begin
    stop=1b
    return,0
endif

if (*ptr).codeproc[1,2] then begin
    tmp=where((*(*ptr).asu.I).togglefix eq 0,n)
    ind=(*(*ptr).preverse_indices)[0:(*ptr).nasupos-1]
    n+=total((*(*ptr).asu.I)[ind].wyck.xyz ne 999 and ~(*(*ptr).asu.I)[ind].wyck.togglefix,/pres)
endif else n=0

return,n
end;function Type12nASUParam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro toggleASU,ptr,sel
;sel: [ left, top, right,  bottom ]
; TODO: don't toggle the same coordinate twice, e.g. (x,2x,0)

sel[[0,2]]-=(*ptr).indasused[0]
if sel[0] lt 0 and sel[2] lt 0 then return
sel>=0
p=(*ptr).ASU.I

ncolxyz=(sel[2]<2) - sel[0] +1
if ncolxyz gt 0 then begin
    ind=indgen(3,(*ptr).nasuatoms)
    ind=ind[sel[0]:(sel[2]<2),sel[1]:sel[3]]
    n=n_elements(ind)
    
    btoggle=bytarr(3,(*ptr).nasuatoms)
    
    for j=0l,n-1 do begin
        colj=ind[j] mod 3
        rowj=ind[j]/3
        jj=value_locate(*(*ptr).preverse_indices,rowj)
        b=(*(*ptr).preverse_indices)[jj]
        e=(*(*ptr).preverse_indices)[jj+1]-1

        RT=Symop64((*p)[b].wyck.Rep)
        ind_=where((RT.Rot)[*,colj] ne 0,ct)
        if ct ne 0 then for k=b,e do begin
            ind2=where(btoggle[ind_,k] eq 0,ct2)
            if ct2 ne 0 then begin
                ind_=ind_[ind2]
                (*p)[k].wyck.togglefix[ind_]= ~(*p)[k].wyck.togglefix[ind_]
                btoggle[ind_,k]=1
            endif
        endfor
    endfor
endif

; SOF
i=3
if sel[0] le i and sel[2] ge i then begin
    n=sel[3]-sel[1]+1
    for j=sel[1],sel[1]+n-1 do (*p)[j].togglefix[0] =~(*p)[j].togglefix[0]
endif

; Biso
i=4
if sel[0] le i and sel[2] ge i then begin
    n=sel[3]-sel[1]+1
    for j=sel[1],sel[1]+n-1 do (*p)[j].togglefix[1] =~(*p)[j].togglefix[1]
endif

; Reset toggle:
n=Type12nASUParam(ptr)
ind=*(*ptr).paramind[1,1]
tmp=where((*ptr).togglefix[ind] eq 0,ct)
n+=ct
if n eq 0 then begin
    (*ptr).togglefix[ind]=0
    (*p).togglefix=0
    (*p).wyck.togglefix=0
    
    ModelTypeSetCode,ptr,5,0,ibit=1
endif

end;pro toggleASU
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro updatespacegroup,ev
ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
if not widget_info(ID,/VALID_ID) then return
widget_control,ID,get_uvalue=list2
widget_control,ev.top,get_uvalue=top
widget_control,top,get_uvalue=list

sg=SpaceGroupStructure_event(ev,list.tree)

error=sg.sghash eq 0
if error then begin
    bdrop=widget_info(ev.id,/type) eq 8
    BranchInfo,ev.top,/update,nodropupdate=bdrop
endif else begin
    result=ExpData(list,/xd,/xtt,/y,/dlim,/lambda,/LGPtype)
    ind=where(tag_names(ev) eq 'ACTION',ct)
    SmallDev=*(*(*list.ptrmodel).ptrfit).misc[3]
    EditTypeStructure,ev,list2.ptr,sg,celparam,result,reset=ct ne 0,SmallDev=SmallDev
    BranchInfo,ev.top
endelse

end;pro updatespacegroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RefinedToInitial,ptr,inverse=inverse
update=0b
if (*(ptr)).peaknr ne 0 then begin
    icurrent=0
    if keyword_set(inverse) then begin
        *((*(*ptr).peakparamIO.O)[icurrent])=*((*ptr).peakparamIO.I)
        *((*(*ptr).globalparamIO.O)[icurrent])=*((*ptr).globalparamIO.I)
    endif else begin
        *((*ptr).peakparamIO.I)=*((*(*ptr).peakparamIO.O)[icurrent])
        *((*ptr).globalparamIO.I)=*((*(*ptr).globalparamIO.O)[icurrent])
    endelse
    
    update=1b
endif
update or= RefinedToInitialASU(ptr,inverse=inverse)
return,update
end;function RefinedToInitial        
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefinedToInitialType,ptr,type,ev,inverse=inverse

if size(type,/type) eq 10 then begin
    update=RefinedToInitial(type,inverse=inverse)
    if update then BranchInfo,ev.top,/update
endif else begin
    while PTR_VALID((*ptr).next) do begin
        ptr=(*ptr).next
        if (*ptr).type eq type and (*ptr).include then update=RefinedToInitial(ptr,inverse=inverse)
    endwhile
endelse

end;pro RefinedToInitialType
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro BranchFileEvent,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
0:    begin ;base
;    event=TAG_NAMES(ev, /STRUCTURE_NAME)
;    if event eq 'WIDGET_CONTEXT' then begin
;        ; Show context menu for this group
;        ID=WIDGET_INFO(ev.top, FIND_BY_UNAME = 'ContextPeak')
;        WIDGET_DISPLAYCONTEXTMENU, ev.ID, ev.X, ev.Y, ID
;        return
;    endif
    endcase
1:    begin ;button
    widget_control,ev.id,get_value=val
    case val of
    'Highlight':begin
                ID=GetTableID(ev)
                if not widget_info(ID,/VALID_ID) then return
                sel=widget_info(ID,/TABLE_SELECT)
                n=sel[3]-sel[1]+1
                ;if n eq 1 then return
                i=sel[1];+indgen(n)
                
                ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                widget_control,ID,get_uvalue=list2
                
                widget_control,ev.top,get_uvalue=top
                widget_control,top,get_uvalue=list
                list.highlightpeak.ptr=list2.ptr
                list.highlightpeak.i=i
                widget_control,top,set_uvalue=list
                
                BranchInfo,ev.top,/update
                endcase
    'Delete':    begin
                ID=GetTableID(ev)
                if not widget_info(ID,/VALID_ID) then return
                sel=widget_info(ID,/TABLE_SELECT)
                n=sel[3]-sel[1]+1
                if n eq 0 then return
                ind=indgen(n)+sel[1]
                ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                widget_control,ID,get_uvalue=list2
                DeletePeaksFromGroup,list2.ptr,ind
                BranchInfo,ev.top
                endcase
    'Move':        begin
                ID=GetTableID(ev)
                if not widget_info(ID,/VALID_ID) then return
                sel=widget_info(ID,/TABLE_SELECT)
                n=sel[3]-sel[1]+1
                if n eq 0 then return
                ind=indgen(n)+sel[1]
                ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                widget_control,ID,get_uvalue=list2
                DeletePeaksFromGroup,list2.ptr,ind,/move,top=ev.top
                BranchInfo,ev.top
                endcase
    'Add ASU atom':begin
                widget_control,widget_info(ev.id,/parent),get_uvalue=struct
                ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                widget_control,ID,get_uvalue=list2

                ; get selection info
                sel=widget_info(struct.ID,/TABLE_SELECT)
                ind=where(sel eq -1,ct) ;[ left, top, right,  bottom ]
                if (ct ne 0) then return
                ny=sel[3]-sel[1]+1
                if (ny lt 1) then return

                AddASUposI,list2.ptr,sel[1]
                BranchInfo,ev.top
                endcase
    'Delete ASU atom':begin
                widget_control,widget_info(ev.id,/parent),get_uvalue=struct
                ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                widget_control,ID,get_uvalue=list2

                ; get selection info
                sel=widget_info(struct.ID,/TABLE_SELECT)
                ind=where(sel eq -1,ct) ;[ left, top, right,  bottom ]
                if (ct ne 0) then return
                ny=sel[3]-sel[1]+1
                if (ny lt 1) then return
                ind=sel[1]+indgen(ny)

                DeleteASUposI,list2.ptr,ind
                BranchInfo,ev.top
                endcase
    'Calc':        begin
                widget_control,widget_info(ev.id,/parent),get_uvalue=struct
                   InitGlobal,ev,struct.i-1
                   endcase
    'Toggle fix':begin
                widget_control,widget_info(ev.id,/parent),get_uvalue=struct
                ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                widget_control,ID,get_uvalue=list2

                ; get selection info
                sel=widget_info(struct.ID,/TABLE_SELECT)
                ind=where(sel eq -1,ct) ;[ left, top, right,  bottom ]
                if (ct ne 0) then return

                if struct.i eq 7 then begin
                    toggleASU,list2.ptr,sel
                endif else begin
                    nx=sel[2]-sel[0]+1

                    ind=*((*list2.ptr).paramind[struct.i-1,1])
                    subind=ind[sel[0]:sel[2]]

                    (*list2.ptr).togglefix[subind] = ~(*list2.ptr).togglefix[subind]
                    if (*list2.ptr).type ne 10 then begin
                        if ptr_valid((*list2.ptr).pstrinfo) then $
                            (*list2.ptr).togglefix[indCelParam(list2.ptr)]or=celparamfix(*(*list2.ptr).pstrinfo)
                    endif

                    ; Reset toggle?
                    if (*list2.ptr).type eq 12 and struct.i eq 2 then begin
                        p=(*list2.ptr).ASU.I
                        ; Rietveld scaling factor
                        if ptr_valid(p) then begin
                            tmp=(*list2.ptr).togglefix[ind]
                            add=(*p).togglefix
                            tmp=[tmp,reform(add,n_elements(add))]
                            tmp2=where((*p).wyck.xyz ne 999,ct)
                            if ct ne 0 then begin
                                add=((*p).wyck.togglefix)[tmp2]
                                tmp=[tmp,reform(add,n_elements(add))]
                            endif
                            tmp=where(tmp eq 0,ct)
                        endif
                    endif else tmp=where((*list2.ptr).togglefix[ind] eq 0,ct)
                    if ct eq 0 then begin
                        (*list2.ptr).togglefix[ind]=0
                        if n_elements(p) ne 0 then begin
                            (*p).togglefix=0
                            (*p).wyck.togglefix=0
                        endif
                        
                        ModelTypeSetCode,list2.ptr,struct.i+3,0,ibit=1
                    endif
                endelse

                BranchInfo,ev.top,/update
                endcase
    'Init<>Ref':begin
                ; Get model info
                widget_control,widget_info(ev.id,/parent),get_uvalue=struct

                widget_control,widget_info(struct.ID,/parent),get_uvalue=uval2
                case uval2 of
                    'InitialTab':     begin
                                    type=0
                                    uname_get='RefinedTab'
                                    endcase
                    'RefinedTab':     begin
                                    type=2
                                    uname_get='InitialTab'
                                    endcase
                    'InitialTabASU':begin
                                    type=4
                                    uname_get='RefinedTabASU'
                                    endcase
                    'RefinedTabASU':begin
                                    type=6
                                    uname_get='InitialTabASU'
                                    endcase
                    else: return
                endcase
                ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                widget_control,ID,get_uvalue=list2

                ; get selection info
                sel=widget_info(struct.ID,/TABLE_SELECT)
                ind=where(sel eq -1,ct) ;[ left, top, right,  bottom ]
                if (ct ne 0) then return
                nx=sel[2]-sel[0]+1
                ny=sel[3]-sel[1]+1
                ;if ((nx ne 1) and (ny ne 1)) then return
                widget_control,struct.ID,get_value=array

                ; get selection info
                ID=widget_info(ev.top,FIND_BY_UNAME=uname_get)
                widget_control,struct.ID,get_uvalue=TableID
                ID=widget_info(ID,FIND_BY_UNAME=TableID)
                widget_control,ID,get_value=array2

                array[sel[0]:sel[2],sel[1]:sel[3]]=array2[sel[0]:sel[2],sel[1]:sel[3]]

                ; update table and model
                widget_control,struct.ID,set_value=array
                EditTypeTableev,struct,struct.i,type,sel[[0,2]],sel[[1,3]],list2.ptr,array,/no_obsolete
                endcase
    'Duplicate':        begin
                ; Get model info
                widget_control,widget_info(ev.id,/parent),get_uvalue=struct
                widget_control,widget_info(struct.ID,/parent),get_uvalue=uval2
                case uval2 of
                    'InitialTab': type=0
                    'ConstraintsTab':type=1
                    'RefinedTab': type=2
                    'SDTab': type=3
                endcase
                ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                widget_control,ID,get_uvalue=list2

                ; get selection info
                sel=widget_info(struct.ID,/TABLE_SELECT)
                ind=where(sel eq -1,ct) ;[ left, top, right,  bottom ]
                if (ct ne 0) then return
                nx=sel[2]-sel[0]+1
                ny=sel[3]-sel[1]+1
                if ((nx ne 1) and (ny ne 1)) then return

                ; make array
                widget_control,struct.ID,get_value=array
                s=size(array)
                nxmax=s[1]
                case s[0] of
                    1:    nymax=1
                    2:    nymax=s[2]
                endcase
                arraysel=array[sel[0]:sel[2],sel[1]:sel[3]]

                if (struct.i eq 0) and (type ne 1) then begin ; copy vertical
                    if ny ne 1 then return
                    array[sel[0]:sel[2],*]=arraysel#replicate(1b,nymax)
                    sel[1]=0
                    sel[3]=nymax-1
                endif else begin ; copy horizontal
                    if nx ne 1 then return
                    array[*,sel[1]:sel[3]]=replicate(1b,nxmax)#arraysel
                    sel[0]=0
                    sel[2]=nxmax-1
                endelse

                ; update table and model
                widget_control,struct.ID,set_value=array
                EditTypeTableev,struct,struct.i,type,sel[[0,2]],sel[[1,3]],list2.ptr,array
                endcase
    'Global':    begin
                MainCheckEvent,ev
                endcase
    'Refined':    begin
                MainCheckEvent,ev
                endcase
    'Refined (constraints +/-C)':begin
                MainCheckEvent,ev
                endcase
    'Refined (constraints >= Clow)':begin
                MainCheckEvent,ev
                endcase
    'Refined (constraints <= Chigh)':begin
                MainCheckEvent,ev
                endcase
    endcase

    endcase
3:    begin ;text
    widget_control,ev.id,get_uvalue=uval
    case size(uval,/type) of
    10:    begin
        widget_control,ev.id,get_value=val
        *uval=float(val[0])
        BranchInfo,ev.top,/update
        endcase
    else : updatespacegroup,ev
    endcase
    endcase
8:    begin ;droplist
    widget_control,ev.id,get_uvalue=uval
    if size(uval,/type) ne 8 then uval={type:'unknown'}
    case uval.type of
    'dropprofiletype':MainDropEvent,ev,uval.i
    'dropasymtype':MainDropEvent,ev,uval.i
    'dropItype':MainDropEvent,ev,uval.i
    'dropfittype':RefreshMainTab,ev.top
    'ContextPeak': begin
                    widget_control,ev.id,SET_DROPLIST_SELECT=0
                    widget_control,ev.id,get_value=val

                        case val[ev.index] of
                        'Perform...':
                        '':
                        'Add Peaks Manually':begin
                            widget_control,ev.top,get_uvalue=top
                            modeCHI,{top:top},'dummy'
                            widget_control,top,get_uvalue=list,sensitive=1
                            Widget_Control, list.draw, Event_Pro='PeakManual_Event',draw_button=1
                               Widget_Control, list.draw, set_uvalue=ev.top
                            endcase
                        'Add Peaks Automatically':begin
                            widget_control,ev.top,get_uvalue=top
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            PeakAutomatic,list2.ptr,top
                            DeletePeaksFromGroup,list2.ptr,/filter
                            BranchInfo,ev.top
                            endcase
                        'Load From PDF':    begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            PeakFromPDF,list2.ptr,ev
                            endcase
                        'Delete Selected Peaks':begin
                            ID=GetTableID(ev)
                            if not widget_info(ID,/VALID_ID) then return
                            sel=widget_info(ID,/TABLE_SELECT)
                            n=sel[3]-sel[1]+1
                            if n eq 0 then return
                            ind=indgen(n)+sel[1]
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            DeletePeaksFromGroup,list2.ptr,ind
                            BranchInfo,ev.top
                            endcase
                        'Delete All Peaks':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            DeletePeaksFromGroup,list2.ptr,/all
                            BranchInfo,ev.top
                            endcase
                        'Delete Invalid Peaks...':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            DeletePeaksFromGroup,list2.ptr,/filter
                            BranchInfo,ev.top
                            endcase
                        'Delete Peaks With...':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            DeletePeaksFromGroup,list2.ptr,top=ev.top,/user
                            BranchInfo,ev.top
                            endcase
                        'Delete Small Peaks...':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            DeletePeaksFromGroup,list2.ptr,top=ev.top,/dlow
                            BranchInfo,ev.top
                            endcase
                        'Delete duplicates':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            DeletePeaksFromGroup,list2.ptr,/duplicates
                            BranchInfo,ev.top
                            endcase
                        'Sort peaks':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            
                            ptr=list2.ptr
                            p=(*ptr).peakparamIO.I

                            ind=*(*ptr).paramind[0,0]
                            ind=sort((*p)[ind,*])
                            
                            (*p)=(*p)[*,ind]
                            p=(*ptr).peakparamIO.O
                            icurrent=0
                            (*((*p)[icurrent]))=(*((*p)[icurrent]))[*,ind]
                            p=(*ptr).peakparamIO.SD
                            (*((*p)[icurrent]))=(*((*p)[icurrent]))[*,ind]
                            
                            BranchInfo,ev.top,/update
                            endcase
                        'Isodeform=1':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            if not widget_info(ID,/VALID_ID) then return
                            widget_control,ID,get_uvalue=list2
                            widget_control,ev.top,get_uvalue=top
                            widget_control,top,get_uvalue=list
                            
                            ptr=list2.ptr
                            p=(*ptr).peakparamIO.I
                            pg=(*ptr).globalparamIO.I

                            ind=*(*ptr).paramind[0,0]
                            inds=*(*ptr).paramind[0,1]
                            
                            iso=(*pg)[inds[1],*]
                            if iso ne 1 then begin
                                zero=(*pg)[inds[0],*]
                                (*pg)[inds[0],*]=0
                                ret=gfuncPeakparam(list,ptr,/tt)
                                (*p)[ind,*]=ret[0,*]
                                (*pg)[inds[0],*]=zero
                                (*pg)[inds[1]]=1

                                BranchInfo,ev.top,/update
                            endif
                            
                            endcase
                        'Shift=0':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            if not widget_info(ID,/VALID_ID) then return
                            widget_control,ID,get_uvalue=list2
                            widget_control,ev.top,get_uvalue=top
                            widget_control,top,get_uvalue=list
                            
                            ptr=list2.ptr
                            p=(*ptr).peakparamIO.I
                            pg=(*ptr).globalparamIO.I

                            ind=*(*ptr).paramind[0,0]
                            inds=*(*ptr).paramind[0,1]
                            
                            zero=(*pg)[inds[0],*]
                            if zero ne 0 then begin
                                iso=(*pg)[inds[1],*]
                                (*pg)[inds[1],*]=1
                                ret=gfuncPeakparam(list,ptr,/tt)
                                (*p)[ind,*]=ret[0,*]
                                (*pg)[inds[1],*]=iso
                                (*pg)[inds[0]]=0

                                BranchInfo,ev.top,/update
                            endif
                            
                            endcase
                        'Scaling=1':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            if not widget_info(ID,/VALID_ID) then return
                            widget_control,ID,get_uvalue=list2
                            
                            ptr=list2.ptr
                            p=(*ptr).peakparamIO.I
                            pg=(*ptr).globalparamIO.I

                            ind=*(*ptr).paramind[1,0]
                            inds=*(*ptr).paramind[1,1]

                            S=(*pg)[inds,*]
                            if S ne 0 then begin
                                (*p)[ind,*]*=S[0]
                                (*pg)[inds]=1
    
                                BranchInfo,ev.top,/update
                            endif
                            endcase
                        'Save Peaks as PDF...':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            if not widget_info(ID,/VALID_ID) then return
                            widget_control,ID,get_uvalue=list2
                            widget_control,ev.top,get_uvalue=top
                            widget_control,top,get_uvalue=list

                            SavePKSTypeX,list,list2.ptr
                            endcase
                        'Save Peak param...':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            if not widget_info(ID,/VALID_ID) then return
                            widget_control,ID,get_uvalue=list2
                            widget_control,ev.top,get_uvalue=top
                            widget_control,top,get_uvalue=list

                            TestSaveTypeX,list,list2.ptr,type=2
                            endcase
                        'Save FWHM...':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            if not widget_info(ID,/VALID_ID) then return
                            widget_control,ID,get_uvalue=list2
                            widget_control,ev.top,get_uvalue=top
                            widget_control,top,get_uvalue=list

                            TestSaveTypeX,list,list2.ptr,type=0
                            endcase
                        'Save SNR...':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            if not widget_info(ID,/VALID_ID) then return
                            widget_control,ID,get_uvalue=list2
                            widget_control,ev.top,get_uvalue=top
                            widget_control,top,get_uvalue=list

                            TestSaveTypeX,list,list2.ptr,type=1
                            endcase
                        'Load Structure File':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            PeakFromStruc,list2.ptr,ev
                            endcase
                        'Save Structure File':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            StrucFromPeak,list2.ptr,ev
                            endcase
                        'Reset Unit Cell':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='hall')
                            if widget_info(ID,/valid) then begin
                                handler=widget_info(ev.top,FIND_BY_UNAME='groupinfochild')
                                widget_control,ID,send_event={ID:ID,TOP:ev.top,HANDLER:handler,action:'reset'}
                                return
                            endif
                            endcase
                        'Recalculate HKL':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            if (*(list2.ptr)).peaknr eq 0 then return
                            ReCalculateHKL,ev,list2.ptr
                            endcase
                        'Add ASU position':begin
                            ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
                            widget_control,ID,get_uvalue=list2
                            AddASUposI,list2.ptr
                            BranchInfo,ev.top
                            endcase
                        endcase
                        endcase
    else:     updatespacegroup,ev
    endcase
    endcase
9:    begin ;table
    widget_control,ev.id,get_uvalue=uval
    str=['Peak parameters:','Global position parameters:','Global intensity parameters:','Global FWHM parameters:','Global mixing parameters:',$
    'Global decay parameters:','Global asymmetry parameters:','ASU parameters:']
    i=(where(str eq uval))[0]

    widget_control,widget_info(ev.id,/parent),get_uvalue=uval2
    case uval2 of
        'InitialTab': type=0
        'ConstraintsTab':type=1
        'RefinedTab': type=2
        'SDTab': type=3
        'InitialTabASU': type=4
        'ConstraintsTabASU':type=5
        'RefinedTabASU': type=6
        'SDTabASU': type=7
        'Info': type=8
    endcase

    ID=widget_info(ev.top,FIND_BY_UNAME='topbranch')
    widget_control,ID,get_uvalue=list2

    event=TAG_NAMES(ev, /STRUCTURE_NAME)
    if event eq 'WIDGET_CONTEXT' then begin
        ; Show context menu for this group
        if type le 3 then begin
            if type eq 1 then begin
                if ((*list2.ptr).type eq 10) or (i ne 1) then uname='ContextTable'
            endif else begin

                    if i eq 0 then uname='ContextTable2' else begin ; peakparam
                        if ((*list2.ptr).type ne 10) and (i eq 1) then uname='ContextTable3' $ ; celparam
                        else uname='ContextTable1' ; other param
                    endelse
            endelse
        endif else begin ; ASU tabs
            uname='ContextTable4'
        endelse

        if n_elements(uname) ne 0 then begin
            ID=WIDGET_INFO(ev.top, FIND_BY_UNAME = uname)
            widget_control,ID,set_uvalue={ID:ev.id,top:ev.top,i:i}
            WIDGET_DISPLAYCONTEXTMENU, ev.ID, ev.X, ev.Y, ID
        endif
        return
    endif

    case ev.type of
    0:    begin; edit field
        widget_control,ev.id,get_value=array
        EditTypeTableev,ev,i,type,[ev.x,ev.x],[ev.y,ev.y],list2.ptr,array
        endcase
    else:
    endcase

    endcase
10:    begin ;tab
    return
    endcase
endcase

SendMainRefresh,ev.top
end;pro BranchFileEvent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TypenParam,ptr

n=Type12nASUParam(ptr,stop=stop)
if stop then return,0

bglobal=(*ptr).codeproc[*,1] ; global?
ind=where((*ptr).codeproc[*,2],ct) ; refined?
peaknr=(*ptr).peaknr

if ct eq 0 or peaknr eq 0 then return,n

for i=0l,ct-1 do begin
    j=ind[i]
    
    tmp=*(*ptr).paramind[j,bglobal[j]]
    if bglobal[j] then begin
        tmp=where((*ptr).togglefix[tmp] eq 0,ct)
        n+=ct
        if (*ptr).codeproc[j,6] eq 2 then begin
            tmp=*(*ptr).paramind[j,0]
            ind2=where(~(*ptr).peakhelpparam[tmp],ct2)
            n+=ct2*peaknr
        endif 
    endif else begin
        ind2=where(~(*ptr).peakhelpparam[tmp],ct2)
        n+=ct2*peaknr
    endelse
endfor

return,n

end;function TypenParam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TotalTypenParam,ptrFit,batch=batch

mult_stdev=keyword_set(batch)+1

n=0
ptr=(*ptrFit).next
while ptr_valid(ptr) do begin
    if (*ptr).include then begin
        n+=mult_stdev*TypenParam(ptr)
        if ptr_valid((*ptr).factors) and batch then n+=n_elements(*(*(*ptr).factors)[0])
    endif
    ptr=(*ptr).next
endwhile

if batch then begin
    n+=GetNumberPhysicalPhaseInfo(ptrFit) ; Rietveld Scaling derived
    if n ne 0 then n+=n_elements(*(*(*ptrFit).factors)[0]) ;fit info
endif
return,n
end;function TotalTypenParam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro BranchFileInfo,ptr,evtop,update=update,nodropupdate=nodropupdate

if keyword_set(update) then begin
    base=widget_info(evtop,FIND_BY_UNAME='groupinfochild')
    ID=widget_info(base,FIND_BY_UNAME='Group name:')
    widget_control,ID,set_value='Group name: '+(*ptr).name
    ID=widget_info(base,FIND_BY_UNAME='Group number:')
    widget_control,ID,set_value='Group number: '+string((*ptr).group)
    ID=widget_info(base,FIND_BY_UNAME='Included:')
    widget_control,ID,set_value='Included: '+(((*ptr).include)?'yes':'no ')
    ID=widget_info(base,FIND_BY_UNAME='Number of fit parameters:')
    widget_control,ID,set_value='Number of fit parameters:'+string(TypenParam(ptr))
    if (*ptr).type ne 10 then basestruc=widget_info(base,FIND_BY_UNAME='SymmTab')
    baseinit=widget_info(base,FIND_BY_UNAME='InitialTab')
    basecon=widget_info(base,FIND_BY_UNAME='ConstraintsTab')
    baseref=widget_info(base,FIND_BY_UNAME='RefinedTab')
    basesd=widget_info(base,FIND_BY_UNAME='SDTab')
    baseinitasu=widget_info(base,FIND_BY_UNAME='InitialTabASU')
    baseconasu=widget_info(base,FIND_BY_UNAME='ConstraintsTabASU')
    baserefasu=widget_info(base,FIND_BY_UNAME='RefinedTabASU')
    basesdasu=widget_info(base,FIND_BY_UNAME='SDTabASU')
    
    if ptr_valid((*ptr).factors) then begin
        ID=widget_info(base,FIND_BY_UNAME='rfactors')
        widget_control,ID,set_value=*(*(*ptr).factors)[0]
        
        ID=widget_info(base,FIND_BY_UNAME='physinfo')
        if widget_info(ID,/valid) then $
            widget_control,ID,set_value=transpose((*ptr).propfix)
    endif
                
endif else begin

    case (*ptr).type of
    ; Add 'test'
    10:    begin
        dropperform=['Perform...','Add Peaks Manually','Add Peaks Automatically','Load From PDF','',$
                'Delete All Peaks','Delete Invalid Peaks...','Delete Peaks With...', 'Delete Small Peaks...','Delete duplicates','','Sort peaks','Scaling=1','Isodeform=1','Shift=0','','Save Peak param...','Save FWHM...','Save SNR...','Save Peaks as PDF...']
        endcase
    11: begin
        dropperform=['Perform...','Load Structure File', 'Reset Unit Cell','Recalculate HKL','Delete Peaks With...', 'Delete Small Peaks...','','Scaling=1','','Save Peak param...','Save FWHM...','Save SNR...','Save Peaks as PDF...']
        endcase
    12: begin
        dropperform=['Perform...','Load Structure File','Save Structure File','Add ASU position','', 'Delete Peaks With...', 'Delete Small Peaks...','','Reset Unit Cell','Recalculate HKL','','Save Peak param...','Save FWHM...','Save SNR...','Save Peaks as PDF...']
        endcase
    endcase

    ID=widget_info(evtop,FIND_BY_UNAME='groupinfochild')
    tabsel=0
    dropfitsel=0
    if ID ne 0 then begin
        ID2=widget_info(ID,FIND_BY_UNAME='groupinfochildtab')
        if widget_info(ID2,/VALID_ID) then tabsel=widget_info(ID2,/TAB_CURRENT)
        ID2=widget_info(ID,FIND_BY_UNAME='dropfittype')
        if widget_info(ID2,/VALID_ID) then dropfitsel=widget_info(ID2,/DROPLIST_SELECT)
        widget_control,ID,/destroy
    endif

    base=widget_base(widget_info(evtop,FIND_BY_UNAME='groupinfoparent'),uname='groupinfochild',$
        /column,EVENT_PRO='BranchFileEvent')

    baset=widget_base(base,/frame,/column)
        label=widget_label(baset,value='Group name: '+(*ptr).name,/ALIGN_LEFT,uname='Group name:')
        label=widget_label(baset,value='Group number: '+string((*ptr).group),/ALIGN_LEFT,uname='Group number:')
        label=widget_label(baset,value='Included: '+(((*ptr).include)?'yes':'no '),/ALIGN_LEFT,uname='Included:')
        label=widget_label(baset,value='Number of fit parameters:'+string(TypenParam(ptr)),/ALIGN_LEFT,uname='Number of fit parameters:')

        drop=widget_droplist(baset,value=dropperform,uvalue={type:'ContextPeak'},uname='ContextPeak')

    baset = WIDGET_BASE(base,  /CONTEXT_MENU, UNAME='ContextTable')
        button=widget_button(baset,value='Duplicate')
    baset = WIDGET_BASE(base,  /CONTEXT_MENU, UNAME='ContextTable1')
        button=widget_button(baset,value='Duplicate')
        button=widget_button(baset,value='Calc')
        button=widget_button(baset,value='Toggle fix')
        button=widget_button(baset,value='Init<>Ref')
    baset = WIDGET_BASE(base,  /CONTEXT_MENU, UNAME='ContextTable2')
        button=widget_button(baset,value='Highlight')
        button=widget_button(baset,value='Duplicate')
        button=widget_button(baset,value='Delete')
        button=widget_button(baset,value='Move')
        button=widget_button(baset,value='Init<>Ref')
    baset = WIDGET_BASE(base,  /CONTEXT_MENU, UNAME='ContextTable3')
        button=widget_button(baset,value='Duplicate')
        button=widget_button(baset,value='Toggle fix')
        button=widget_button(baset,value='Init<>Ref')
    baset = WIDGET_BASE(base,  /CONTEXT_MENU, UNAME='ContextTable4')
        button=widget_button(baset,value='Add ASU atom')
        button=widget_button(baset,value='Delete ASU atom')
        button=widget_button(baset,value='Toggle fix')
        button=widget_button(baset,value='Init<>Ref')

    tab=widget_tab(base,uname='groupinfochildtab')
        baset=widget_base(tab,title='Main',/column)
            ; ----Make table labels----
            label=widget_label(baset,value='Peak profile type:',/ALIGN_LEFT)
            ; To be added
            ; 'Split Pearson VII'
            drop=widget_droplist(baset,value=['Gaussian','Lorentzian','Mod 1 Lorentzian',$
                    'Mod 2 Lorentzian','Pseudo-Voigt','TCHZ Pseudo-Voigt','Pearson VII','Split Pseudo-Voigt','Mod. Split Pseudo-Voigt'],uvalue={type:'dropprofiletype',i:0},uname='dropprofiletype')
                    
            ; To be added
            ;label=widget_label(baset,value='Peak asymmetry correction:',/ALIGN_LEFT)
            ; 'Rietveld','Berar-Baldinozzi'
            ;drop=widget_droplist(baset,value=['None'],uvalue={type:'dropasymtype',i:1},uname='dropasymtype')
            
            if (*ptr).type eq 11 then begin
                label=widget_label(baset,value='Intensity type:',/ALIGN_LEFT)
                drop=widget_droplist(baset,value=['Pawley','LeBail'],uvalue={type:'dropItype',i:2},uname='dropItype')
            endif
            
            label=widget_label(baset,value='Fit parameters:',/ALIGN_LEFT)
            baset2=widget_base(baset,/row)
            drop=widget_droplist(baset2,value=['Position','Area','FWHM','Mixing','Decay','Assym'],uvalue={type:'dropfittype',i:[4,5,6,7,8,9]},uname='dropfittype')
            widget_control,drop,SET_DROPLIST_SELECT=dropfitsel
            basett=widget_base(baset2,/column,/NONEXCLUSIVE)
                button=widget_button(basett,value='Global',uvalue=0,uname='Global')
                button=widget_button(basett,value='Refined',uvalue=1,uname='Refined')
                button=widget_button(basett,value='Refined (constraints >= Clow)',uvalue=2,uname='Refined (constraints >= Clow)')
                button=widget_button(basett,value='Refined (constraints <= Chigh)',uvalue=3,uname='Refined (constraints <= Chigh)')
                button=widget_button(basett,value='Refined (constraints +/-C)',uvalue=4,uname='Refined (constraints +/-C)')
                
            label=widget_label(baset,value=' ',/ALIGN_LEFT)
            label=widget_label(baset,value='Refined: uses some default constraints',/ALIGN_LEFT)
            label=widget_label(baset,value='Constraints +/-C: stay around initial values',/ALIGN_LEFT)
            label=widget_label(baset,value='Constraints >=Clow <=CHigh: limit parameters',/ALIGN_LEFT)
            label=widget_label(baset,value='Refined + refined with constraints (single intensities): reject peaks when outside limits',/ALIGN_LEFT)
            label=widget_label(baset,value='Refined + refined with constraints (scaling factor): reject phase and fit again',/ALIGN_LEFT)
            label=widget_label(baset,value=' ',/ALIGN_LEFT)
        
        baset=widget_base(tab,title='Info',/column,uvalue='Info')
            if (*ptr).indpropfix[0] ne -1 then begin
                table=widget_table(baset,column_labels=['Physical info'],row_labels=(*ptr).propname[(*ptr).indpropfix],$
                                value=transpose((*ptr).propfix),uname='physinfo',uvalue='physinfo',XSIZE=1,YSIZE=n_elements((*ptr).propfix),/editable)
            endif
            
            if ptr_valid((*ptr).factors) then $
                table=widget_table(baset,column_labels=['R-factors'],row_labels=['Rb(%)','Rf(%)'],value=*(*(*ptr).factors)[0],uname='rfactors',uvalue='rfactors',$
                                    XSIZE=1,YSIZE=2)

        if (*ptr).type ne 10 then basestruc=widget_base(tab,title='Symmetry',/column,uname='SymmTab',uvalue='SymmTab')
        baseinit=widget_base(tab,title='Initial',/column,uname='InitialTab',uvalue='InitialTab')
        basecon=widget_base(tab,title='Constraints',/column,uname='ConstraintsTab',uvalue='ConstraintsTab')
        baseref=widget_base(tab,title='Refined',/column,uname='RefinedTab',uvalue='RefinedTab')
        basesd=widget_base(tab,title='Standard Dev.',/column,uname='SDTab',uvalue='SDTab')
        if (*ptr).type eq 12 then baseinitasu=widget_base(tab,title='ASU(init)',/column,uname='InitialTabASU',uvalue='InitialTabASU')
        if (*ptr).type eq 12 then baseconasu=widget_base(tab,title='ASU(constr)',/column,uname='ConstraintsTabASU',uvalue='ConstraintsTabASU')
        if (*ptr).type eq 12 then baserefasu=widget_base(tab,title='ASU(ref)',/column,uname='RefinedTabASU',uvalue='RefinedTabASU')
        if (*ptr).type eq 12 then basesdasu=widget_base(tab,title='ASU(SD)',/column,uname='SDTabASU',uvalue='SDTabASU')
    widget_control,tab,SET_TAB_CURRENT=tabsel
endelse
RefreshMainTab,evtop

if (*ptr).type ne 10 then begin
    widget_control,evtop,get_uvalue=top
    widget_control,top,get_uvalue=list
    typeTablestruc,basestruc,list.tree,(*ptr).pstrinfo,update=update,nodropupdate=nodropupdate
endif
typeTable,baseinit,ptr,update=update
typeTable,basecon,ptr,/constraints,update=update
typeTable,baseref,ptr,/refined,update=update
typeTable,basesd,ptr,/sd,update=update
if (*ptr).type eq 12 then begin
    typeTableasu,baseinitasu,ptr,update=update
    typeTableasu,baseconasu,ptr,/constraints,update=update
    typeTableasu,baserefasu,ptr,/refined,update=update
    typeTableasu,basesdasu,ptr,/sd,update=update
endif

end;pro BranchFileInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro StartRefinement,ev
widget_control,ev.top,get_uvalue=top
ID=WIDGET_INFO(ev.top, FIND_BY_UNAME = 'correlation')
b=WIDGET_INFO(ID,/valid)
correlation=b

ID2=WIDGET_INFO(ev.top, FIND_BY_UNAME = 'debug')
debug=widget_info(ID2,/BUTTON_SET)

FitModelCHI,top,correlation=correlation,debug=debug
if b then widget_control,ID,set_uvalue=correlation,/no_copy
end;pro StartRefinement
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro BranchFolderEvent,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
1 : begin ;button event
    widget_control,ev.id,get_value=val
    case val of
    'Start refinement':begin
                    widget_control,/hourglass
                    StartRefinement,ev
                    BranchInfo,ev.top,/update
                    endcase
    'Scaling factors from weights ...':begin
                    widget_control,ev.id,get_uvalue=O
                    SetScalingToWeight,ev.top,O=O
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    SetPhysicalPhaseInfo,(*list.ptrModel).ptrFit,list.lambda,O=O
                    BranchInfo,ev.top,/update
                    SendMainRefresh,ev.top
                    endcase
    'Plot parameter correlation':begin
    
                    widget_control,ev.id,get_uvalue=info
                    if size(info,/type) ne 8 then return
                    nA=n_elements(info.Areturn)
                    
                    ; very weak: 0%-10%
                    ; weak: 10%-40%
                    ; medium: 40%-70%
                    ; strong: 70%-100%
                    img=info.corr ge 0.10
                    ind=where(info.corr ge 0.40,ct)
                    if ct ne 0 then img[ind]++
                    ind=where(info.corr ge 0.70,ct)
                    if ct ne 0 then img[ind]++
                    
                    loadctrev,-3
                    wset,0
                    erase
                    nmag=floor(0.9*!d.Y_VSIZE/(nA>4))>1
                    ;TVcolorbar,[0.8,0.15,0.9,0.95],range=[0,3]
                    imgtvscl,indgen(1,4),0.7,0.05,nmag=nmag,xtickinterval=1,ytickinterval=1,xtickname=' ',ytickname=['<10%','10%-40%','40%-70%','>70%']
                    imgtvscl,img,0.05,0.05,nmag=nmag,xtickinterval=1,ytickinterval=1,/order,/noerase

                    ; Strongly correlated parameters
                    yoff=0.9
                    for k=3,2,-1 do begin
                        corrb=ulong(img eq k)
                        corrb[indgen(nA),indgen(nA)]=1
                        for i=0l, nA-1 do begin
                            ind0=where(corrb[i,i:*],ct)
                            if ct gt 1 and corrb[i,i] then begin
                                x0=replicate(i,ct-1)
                                y0=ind0[1:*]+i
                                ct--
                                
                                for j=0l,ct-1 do begin
                                    a=[x0[j],y0[j]]
                                    b=info.Areturn[a]
                                    str=[stringr(fix(a)),stringr(b)]
                                    str=str[0]+"(="+str[2]+") <-> "+str[1]+"(="+str[3]+")"
                                    j0=a[0]
                                    j1=a[1]
                                    plots,[j0,j0,0],[0,j1,j1],color=100,/data
                                    plots,j0,j1,/data,color=0,psym=-4
                                    ;plots,-[0.5,1.5,1.5,0.5],[j0,j0,j1,j1],color=100,/data
                                    xyouts,0.7,yoff,str,/normal
                                    yoff-=0.03
                                endfor
                            endif
                        endfor
                        yoff-=0.03
                    endfor
                    
                    endcase
    'Numerical Derivatives': begin
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    (*(*list.ptrModel).ptrFit).fitcontrol[3]=ev.select
                    endcase
    'Subtract background before fitting': begin
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    (*(*list.ptrModel).ptrFit).fitcontrol[4]=ev.select
                    endcase
    'Compare numerical and calculated derivatives': begin
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    (*(*list.ptrModel).ptrFit).debug=ev.select
                    endcase
    'ASU':            begin
                    if ev.select eq 0 then return
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    *((*(*list.ptrModel).ptrFit).misc[1])=0
                    endcase
    'Unit Cell':    begin
                    if ev.select eq 0 then return
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    *((*(*list.ptrModel).ptrFit).misc[1])=1
                    endcase
    'Expanded Unit Cell':begin
                    if ev.select eq 0 then return
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    *((*(*list.ptrModel).ptrFit).misc[1])=2
                    endcase
    'Show ASU labels':begin
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    *((*(*list.ptrModel).ptrFit).misc[2])=ev.select
                    endcase
    'Anomalous dispersion':begin
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    *((*(*list.ptrModel).ptrFit).misc[4])=ev.select
                    endcase
    'Zero shift = sample detector distance shift':begin
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    (*(*list.ptrModel).ptrFit).fitcontrol[6]=ev.select
                    endcase
    'Zero shift = 2theta shift':begin
                    widget_control,ev.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    (*(*list.ptrModel).ptrFit).fitcontrol[6]=~ev.select
                    endcase
    'Add line':        begin
                    widget_control,widget_info(ev.id,/parent),get_uvalue=struc
                    
                    sel=widget_info(struc.ID,/TABLE_SELECT)
                    i=sel[1]
                
                    *struc.p=[[*struc.p],[(*struc.p)[*,i]]]
                    SendMainRefresh,struc.top
                    BranchInfo,struc.top,/update
                    endcase
    'Delete line':    begin
                    widget_control,widget_info(ev.id,/parent),get_uvalue=struc
                    sel=widget_info(struc.ID,/TABLE_SELECT)
                    n=sel[3]-sel[1]+1
                    i=sel[1]+indgen(n)
                    np=(Dimsize(*struc.p,2))[1]
                    if np le 1 or n ge np then return
                    tmp=ShrinkArray(struc.p,i,/row)
                    SendMainRefresh,struc.top
                    BranchInfo,struc.top,/update
                    endcase
    'Graphite(002)':begin
                    widget_control,widget_info(ev.id,/parent),get_uvalue=struc
                    *struc.p=3.3480
                    BranchInfo,struc.top,/update
                    endcase
    'Silicon(111)': begin
                    widget_control,widget_info(ev.id,/parent),get_uvalue=struc
                    *struc.p=3.1356
                    BranchInfo,struc.top,/update
                    endcase
    'non-polarized beam':begin
                    widget_control,widget_info(ev.id,/parent),get_uvalue=struc
                    *struc.p=0.
                    BranchInfo,struc.top,/update
                    endcase
    'synchrotron radiation':begin
                    widget_control,widget_info(ev.id,/parent),get_uvalue=struc
                    *struc.p=0.9
                    BranchInfo,struc.top,/update
                    endcase
    else:            begin
                    widget_control,widget_info(ev.id,/parent),get_uvalue=struc
                    widget_control,ev.id,get_uvalue=uval
                    
                    widget_control,struc.top,get_uvalue=top
                    widget_control,top,get_uvalue=list
                    list.lambda=uval[0]
                    widget_control,top,set_uvalue=list
                    RefreshDataCHI,{top:top}
                    widget_control,top,get_uvalue=list
                    list.backranind=[0,list.nxval-1]
                      temp=(*list.xvalt)[list.backranind,0]
                      list.backranindD=temp[sort(temp)]
                      widget_control,top,set_uvalue=list
                    
                    *struc.p=uval
                    SendMainRefresh,struc.top,/data
                    BranchInfo,struc.top,/update
                    endelse
    endcase
    endcase
8:    begin ;droplist
    widget_control,ev.id,get_uvalue=uval
    widget_control,ev.top,get_uvalue=top
    widget_control,top,get_uvalue=list
    case uval of
    'dropback':    SetbitCode,(*list.ptrModel).ptrFit,0,ev.index
    'dropLP':    begin
                SetbitCode,(*list.ptrModel).ptrFit,1,ev.index
                BranchInfo,ev.top,/update
                endcase
    'dropZS':    begin
                SetbitCode,(*list.ptrModel).ptrFit,3,ev.index
                ModelTypeSetCodeAll,(*list.ptrModel).ptrFit
                endcase
    'droprefresh':(*(*list.ptrModel).ptrFit).fitcontrol[2]=ev.index
    'dropweight':(*(*list.ptrModel).ptrFit).fitcontrol[5]=ev.index
    'dropOR': *((*(*list.ptrModel).ptrFit).misc[10])=ev.index?-1:1
    endcase
    SendMainRefresh,ev.top
    endcase
3:    begin
    widget_control,ev.top,get_uvalue=top
    widget_control,top,get_uvalue=list
    
    ; Context event
    widget_control,ev.id,get_uvalue=uval
    event=TAG_NAMES(ev, /STRUCTURE_NAME)
    if event eq 'WIDGET_CONTEXT' then begin
        case uval of
        'Monochromator':     begin
                                            uname='ContextMono'
                                            p=((*(*list.ptrModel).ptrFit).misc[5])
                                            endcase
        'Linear degree of polarization (horizontal)': begin
                                            uname='ContextPol'
                                            p=((*(*list.ptrModel).ptrFit).misc[6])
                                            endcase
        else:return
        endcase
        ID=WIDGET_INFO(ev.top, FIND_BY_UNAME = uname)
        widget_control,ID,set_uvalue={ID:ev.id,top:ev.top,p:p}
        WIDGET_DISPLAYCONTEXTMENU, ev.ID, ev.X, ev.Y, ID
        return
    endif

    ; Edit event
    widget_control,ev.id,get_value=val
    val=val[0]
    case uval of
    'maxiter':(*(*list.ptrModel).ptrFit).fitcontrol[0]=round(val)
    'convtol':(*(*list.ptrModel).ptrFit).fitcontrol[1]=float(val)
    'AtomRadiusMult':*((*(*list.ptrModel).ptrFit).misc[0])=float(val)
    'AtomPrecision':*((*(*list.ptrModel).ptrFit).misc[3])=float(val)
    'Monochromator':*((*(*list.ptrModel).ptrFit).misc[5])=float(val)
    'Linear degree of polarization (horizontal)':*((*(*list.ptrModel).ptrFit).misc[6])=float(val)
    'Number of monochromator crystals': *((*(*list.ptrModel).ptrFit).misc[8])=round(float(val))
    'Mosaic (2) or perfect (1) crystals': *((*(*list.ptrModel).ptrFit).misc[9])=float(val)
    'Scattering azimuth': *((*(*list.ptrModel).ptrFit).misc[11])=float(val)
    endcase
    endcase
9:    begin
    widget_control,ev.id,get_uvalue=uval
    case uval of
    'nCellTable':begin
                widget_control,ev.top,get_uvalue=top
                widget_control,top,get_uvalue=list
                widget_control,ev.id,get_value=val
                ; sort
                val=[val[0,*]<val[1,*],val[0,*]>val[1,*]]
                *((*(*list.ptrModel).ptrFit).misc[7])=val
                widget_control,ev.id,set_value=val
                
                endcase
    'elines':    begin
                widget_control,ev.top,get_uvalue=top
                widget_control,top,get_uvalue=list
                
                event=TAG_NAMES(ev, /STRUCTURE_NAME)
                if event eq 'WIDGET_CONTEXT' then begin
                    uname='ContextELines'
                    ID=WIDGET_INFO(ev.top, FIND_BY_UNAME = uname)
                    widget_control,ID,set_uvalue={ID:ev.id,top:ev.top,p:(*(*list.ptrModel).ptrFit).sourceemission}
                    WIDGET_DISPLAYCONTEXTMENU, ev.ID, ev.X, ev.Y, ID
                    return
                endif
                
                widget_control,ev.id,get_value=val
                *(*(*list.ptrModel).ptrFit).sourceemission=val
                
                if ev.x eq 0 and ev.y eq 0 then begin
                    list.lambda=(*(*(*list.ptrModel).ptrFit).sourceemission)[0]
                    widget_control,top,set_uvalue=list
                    RefreshDataCHI,{top:top}
                    widget_control,top,get_uvalue=list
                    list.backranind=[0,list.nxval-1]
                      temp=(*list.xvalt)[list.backranind,0]
                      list.backranindD=temp[sort(temp)]
                      widget_control,top,set_uvalue=list
                endif
                SendMainRefresh,ev.top
                endcase
    else:
    endcase
    endcase
else:
endcase
end;pro BranchFolderEvent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro BranchFolderInfo,val,evtop,update=update

update=keyword_set(update)
widget_control,evtop,get_uvalue=top
widget_control,top,get_uvalue=list

if not update then begin
    ID=widget_info(evtop,FIND_BY_UNAME='groupinfochild')
    if ID ne 0 then widget_control,ID,/destroy
    base=widget_base(widget_info(evtop,FIND_BY_UNAME='groupinfoparent'),uname='groupinfochild',/column,EVENT_PRO='BranchFolderEvent')
endif else base=widget_info(evtop,FIND_BY_UNAME='groupinfochild')

case val of
    'Model':    begin
                icurrent=(*(*list.ptrModel).ptrFit).icurrent
                pfactors=(*(*(*list.ptrModel).ptrFit).factors)[icurrent]
                if update then begin
                    dropback=widget_info(base,FIND_BY_UNAME='dropback')
                    dropLP=widget_info(base,FIND_BY_UNAME='dropLP')
                    dropZS=widget_info(base,FIND_BY_UNAME='dropZS')
                    dropOR=widget_info(base,FIND_BY_UNAME='dropOR')
                    ID=widget_info(base,FIND_BY_UNAME='relfactors1')
                    widget_control,ID,set_value=fix((*pfactors)[0,0:7])
                    ID=widget_info(base,FIND_BY_UNAME='relfactors2')
                    widget_control,ID,set_value=(*pfactors)[0,8:16]
                    ID=widget_info(base,FIND_BY_UNAME='relfactors3')
                    widget_control,ID,set_value=(*pfactors)[0,17:19]
                    text0=widget_info(base,FIND_BY_UNAME='chi2')
                    text1=widget_info(base,FIND_BY_UNAME='Schwarz')
                    text2=widget_info(base,FIND_BY_UNAME='Akaike')
                    text3=widget_info(base,FIND_BY_UNAME='Durbin')
                    text01=widget_info(base,FIND_BY_UNAME='maxiter')
                    text02=widget_info(base,FIND_BY_UNAME='convtol')
                    drop03=widget_info(base,FIND_BY_UNAME='droprefresh')
                    drop04=widget_info(base,FIND_BY_UNAME='dropweight')
                    button01=widget_info(base,FIND_BY_UNAME='numder')
                    button07=widget_info(base,FIND_BY_UNAME='chisqbk')
                    button08=widget_info(base,FIND_BY_UNAME='debug')
                    text03=widget_info(base,FIND_BY_UNAME='AtomRadiusMult')
                    button02=widget_info(base,FIND_BY_UNAME='ASU')
                    button03=widget_info(base,FIND_BY_UNAME='Unit Cell')
                    button04=widget_info(base,FIND_BY_UNAME='Expanded Unit Cell')
                    button05=widget_info(base,FIND_BY_UNAME='Show ASU labels')
                    text04=widget_info(base,FIND_BY_UNAME='AtomPrecision')
                    table01=widget_info(base,FIND_BY_UNAME='nCellTable')
                    text05=widget_info(base,FIND_BY_UNAME='Monochromator')
                    text06=widget_info(base,FIND_BY_UNAME='Linear degree of polarization (horizontal)')
                    text07=widget_info(base,FIND_BY_UNAME='Number of monochromator crystals')
                    text08=widget_info(base,FIND_BY_UNAME='Mosaic (2) or perfect (1) crystals')
                    text09=widget_info(base,FIND_BY_UNAME='Scattering azimuth')
                    button06=widget_info(base,FIND_BY_UNAME='Anomalous dispersion')
                    tab4a=widget_info(base,FIND_BY_UNAME='tab4a')
                    DeleteWBaseChilds,tab4a
                    tab4b=widget_info(base,FIND_BY_UNAME='tab4b')
                    DeleteWBaseChilds,tab4b
                    tab0=widget_info(base,FIND_BY_UNAME='tab0')
                    tmp=widget_info(tab0,FIND_BY_UNAME='elines')
                    widget_control,tmp,/destroy
                    tmp=widget_info(tab0,FIND_BY_UNAME='ContextELines')
                    widget_control,tmp,/destroy
                endif else begin

                    button=widget_button(base,value='Start refinement')

                    labels=FitInfoLabels()
                    tab=widget_tab(base)

                    tab0main=widget_base(tab,title='General profile parameters',/column)
                    tab0main=widget_tab(tab0main)
                        tab0=widget_base(tab0main,title='Main',/column,uname='tab0')
                            label=widget_label(tab0,value='Background type:',/ALIGN_LEFT)
                            dropback=widget_droplist(tab0,value=['None','Calc+substract (strip)','Calc+substract (ortpol)','Fit (ortpol)'],uvalue='dropback',uname='dropback')
                            
                            label=widget_label(tab0,value='Zero Shift type:',/ALIGN_LEFT)
                            dropZS=widget_droplist(tab0,value=['2theta shift','sample-detector shift'],uvalue='dropZS',uname='dropZS')
                            
                            baset=widget_base(tab0,/row,/NONEXCLUSIVE )
                                button06=widget_button(baset,value='Anomalous dispersion',uname='Anomalous dispersion')
                            label=widget_label(tab0,value='Source emission lines (first = max):',/ALIGN_LEFT)
                        
                        tab0b=widget_base(tab0main,title='Lorentz-polarization',/column)
                        label=widget_label(tab0b,value='Lorentz-polarization correction:',/ALIGN_LEFT)

                        dropLP=widget_droplist(tab0b,value=['None','Not removed in 2D (optics)','Not removed in 2D (pre-diff. monochromator)','Not removed in 2D (post-diff. monochromator)','Removed in 2D'],uvalue='dropLP',uname='dropLP')

                        label=widget_label(tab0b,value='Linear degree of polarization (horizontal):',/ALIGN_LEFT)
                        text06=widget_text(tab0b,uname='Linear degree of polarization (horizontal)',uvalue='Linear degree of polarization (horizontal)',/editable,/context_events)
                        context = WIDGET_BASE(tab0b,  /CONTEXT_MENU, UNAME='ContextPol',Y_SCROLL_SIZE=2)
                            t = widget_button(context,value='non-polarized beam')
                            t = widget_button(context,value='synchrotron radiation')
                            
                        label=widget_label(tab0b,value='Monochromator d-spacing ('+msymbols('angstroms')+'): ',/ALIGN_LEFT)
                        text05=widget_text(tab0b,uname='Monochromator',uvalue='Monochromator',/editable,/context_events)
                        context = WIDGET_BASE(tab0b,  /CONTEXT_MENU, UNAME='ContextMono',Y_SCROLL_SIZE=1)
                            t = widget_button(context,value='Silicon(111)')
                            t = widget_button(context,value='Graphite(002)')
                
                        label=widget_label(tab0b,value='Number of monochromator crystals:',/ALIGN_LEFT)
                        text07=widget_text(tab0b,uname='Number of monochromator crystals',uvalue='Number of monochromator crystals',/editable)
                        
                        label=widget_label(tab0b,value='Mosaic (2) or perfect (1) crystals:',/ALIGN_LEFT)
                        text08=widget_text(tab0b,uname='Mosaic (2) or perfect (1) crystals',uvalue='Mosaic (2) or perfect (1) crystals',/editable)
                        
                        label=widget_label(tab0b,value='Scattering azimuth (degrees):',/ALIGN_LEFT)
                        text09=widget_text(tab0b,uname='Scattering azimuth',uvalue='Scattering azimuth',/editable)
                        
                        label=widget_label(tab0b,value='Monochromator orientation:',/ALIGN_LEFT)
                        dropOR=widget_droplist(tab0b,value=['Vertical','Horizontal'],uvalue='dropOR',uname='dropOR')
    
                        

                    tab1=widget_base(tab,title='General fit results',/column)
                        baset=widget_base(tab1,/row)
                        table1=widget_table(baset,value=fix((*pfactors)[0,0:7]),row_labels=labels[0:7],COLUMN_LABELS=['Number of'],uname='relfactors1')
                        table2=widget_table(baset,value=(*pfactors)[0,8:16],row_labels=labels[8:16],COLUMN_LABELS=['NLLS factors'],uname='relfactors2')
                        table3=widget_table(baset,value=(*pfactors)[0,17:19],row_labels=labels[17:19],COLUMN_LABELS=['Stat. factors'],uname='relfactors3')

                        baset=widget_base(tab1,/row)
                            label=widget_label(baset,value='Chi2 p-value: ',/ALIGN_LEFT)
                            text0=widget_text(baset,value='',uname='chi2')
                        baset=widget_base(tab1,/row)
                            label=widget_label(baset,value='Schwarz criterion(model selection): D+p ln(p) ',/ALIGN_LEFT)
                            text1=widget_text(baset,value='',uname='Schwarz')
                        baset=widget_base(tab1,/row)
                            label=widget_label(baset,value='Akaike criterion(model selection): D+2p ',/ALIGN_LEFT)
                            text2=widget_text(baset,value='',uname='Akaike')
                        baset=widget_base(tab1,/row)
                            label=widget_label(baset,value='Durbin-Watson: Qd< dw <4-Qd: ',/ALIGN_LEFT)
                            text3=widget_text(baset,value='',xsize=40,uname='Durbin',ysize=3)
                            
                        button=widget_button(tab1,value='Plot parameter correlation',uname='correlation')
                    
                    tab4=widget_base(tab,title='Quantitative',/column,uname='tab4')
                    tab4t=widget_tab(tab4)
                    tab4b=widget_base(tab4t,title='Refined',/column,uname='tab4b')
                    tab4a=widget_base(tab4t,title='Initial',/column,uname='tab4a')
                    
                    tab2=widget_base(tab,title='Advanced fit parameters',/column)
                        baset=widget_base(tab2,/row)
                            label=widget_label(baset,value='Maximum number of iterations: ',/ALIGN_LEFT)
                            text01=widget_text(baset,uname='maxiter',uvalue='maxiter',/editable)
                        baset=widget_base(tab2,/row)
                            label=widget_label(baset,value='Convergence tolerance: ',/ALIGN_LEFT)
                            text02=widget_text(baset,uname='convtol',uvalue='convtol',/editable)
                        baset=widget_base(tab2,/row)
                            label=widget_label(baset,value='Refresh display: ',/ALIGN_LEFT)
                            drop03=widget_droplist(baset,value=['never','once','each iteration'],uvalue='droprefresh',uname='droprefresh')
                        baset=widget_base(tab2,/row)
                            label=widget_label(baset,value='NLLS Weights: ',/ALIGN_LEFT)
                            drop04=widget_droplist(baset,value=['none','from file','sqrt(y>0)'],uvalue='dropweight',uname='dropweight')
                        baset=widget_base(tab2,/column,/NONEXCLUSIVE )
                            button07=widget_button(baset,value='Subtract background before fitting',uname='chisqbk')
                            button01=widget_button(baset,value='Numerical Derivatives',uname='numder')
                            button08=widget_button(baset,value='Compare numerical and calculated derivatives',uname='debug')

                    tab3=widget_base(tab,title='Unit cell parameters',/column)
                    baset=widget_base(tab3,/row)
                            label=widget_label(baset,value='Atomic radius multiplier: ',/ALIGN_LEFT)
                            text03=widget_text(baset,uname='AtomRadiusMult',uvalue='AtomRadiusMult',/editable)
                    baset=widget_base(tab3,/row,/EXCLUSIVE )
                            button02=widget_button(baset,value='ASU',uname='ASU')
                            button03=widget_button(baset,value='Unit Cell',uname='Unit Cell')
                            button04=widget_button(baset,value='Expanded Unit Cell',uname='Expanded Unit Cell')
                    baset=widget_base(tab3,/row,/NONEXCLUSIVE )
                            button05=widget_button(baset,value='Show ASU labels',uname='Show ASU labels')
                    baset=widget_base(tab3,/row)
                            label=widget_label(baset,value='Atomic position precision for Wyckoff determination: ',/ALIGN_LEFT)
                            text04=widget_text(baset,uname='AtomPrecision',uvalue='AtomPrecision',/editable)
                    baset=widget_base(tab3,/row)
                            label=widget_label(baset,value='Number of unit cells: ',/ALIGN_LEFT)
                            table01=widget_table(baset,value=(*((*(*list.ptrModel).ptrFit).misc[7])),row_labels=['a','b','c'],COLUMN_LABELS=['Index0','Index1'],uname='nCellTable',uvalue='nCellTable',/editable)
                endelse
                widget_control,dropback,SET_DROPLIST_SELECT=ExtractbitCode((*list.ptrModel).ptrFit,0)
                LPtmp=ExtractbitCode((*list.ptrModel).ptrFit,1)
                widget_control,dropLP,SET_DROPLIST_SELECT=LPtmp
                widget_control,text05,sensitive=LPtmp eq 3 or LPtmp eq 2
                widget_control,text06,sensitive=LPtmp eq 1 or LPtmp eq 2 or LPtmp eq 3
                widget_control,text07,sensitive=LPtmp eq 3 or LPtmp eq 2
                widget_control,text08,sensitive=LPtmp eq 3 or LPtmp eq 2
                widget_control,text09,sensitive=LPtmp eq 1 or LPtmp eq 2

                widget_control,dropZS,SET_DROPLIST_SELECT=ExtractbitCode((*list.ptrModel).ptrFit,3)
                widget_control,dropOR,SET_DROPLIST_SELECT= fix(*((*(*list.ptrModel).ptrFit).misc[10])) eq -1,sensitive=LPtmp eq 3 or LPtmp eq 2
                
                tmp=chisqr_test((*pfactors)[9],(*pfactors)[7],pout=pvalue)
                widget_control,text0,set_value= stringr(pvalue)
                widget_control,text1,set_value= stringr((*pfactors)[20])
                widget_control,text2,set_value= stringr((*pfactors)[21])

                case (*pfactors)[22] of
                1:    tmp=['Negative serial correlation: ','successive values of the residuals','tend to have opposite sign.']
                2:    tmp=['Positive serial correlation: ','successive values of the residuals','tend to have the same sign.']
                3:    tmp=['No serial correlation.']
                else:tmp=''
                endcase

                widget_control,text3,set_value=tmp
                widget_control,text01,set_value= stringr(fix((*(*list.ptrModel).ptrFit).fitcontrol[0]))
                widget_control,text02,set_value= stringr((*(*list.ptrModel).ptrFit).fitcontrol[1])
                widget_control,drop03,SET_DROPLIST_SELECT=fix((*(*list.ptrModel).ptrFit).fitcontrol[2])
                widget_control,drop04,SET_DROPLIST_SELECT=fix((*(*list.ptrModel).ptrFit).fitcontrol[5])
                widget_control,button01,set_button=fix((*(*list.ptrModel).ptrFit).fitcontrol[3])
                widget_control,button07,set_button=fix((*(*list.ptrModel).ptrFit).fitcontrol[4])
                widget_control,button08,set_button=(*(*list.ptrModel).ptrFit).debug
                widget_control,text03,set_value= stringr(*((*(*list.ptrModel).ptrFit).misc[0]))
                widget_control,button02,set_button=fix(*((*(*list.ptrModel).ptrFit).misc[1])) eq 0
                widget_control,button03,set_button=fix(*((*(*list.ptrModel).ptrFit).misc[1])) eq 1
                widget_control,button04,set_button=fix(*((*(*list.ptrModel).ptrFit).misc[1])) eq 2
                widget_control,button05,set_button=fix(*((*(*list.ptrModel).ptrFit).misc[2]))
                widget_control,text04,set_value= stringr(*((*(*list.ptrModel).ptrFit).misc[3]))
                widget_control,button06,set_button=fix(*((*(*list.ptrModel).ptrFit).misc[4]))
                widget_control,text05,set_value= stringr(*((*(*list.ptrModel).ptrFit).misc[5]))
                widget_control,text06,set_value= stringr(*((*(*list.ptrModel).ptrFit).misc[6]))
                widget_control,text07,set_value= stringr(*((*(*list.ptrModel).ptrFit).misc[8]))
                widget_control,text08,set_value= stringr(*((*(*list.ptrModel).ptrFit).misc[9]))
                widget_control,text09,set_value= stringr(*((*(*list.ptrModel).ptrFit).misc[11]))
                widget_control,table01,set_value= *((*(*list.ptrModel).ptrFit).misc[7])

                SetPhysicalPhaseInfo,(*list.ptrModel).ptrFit,list.lambda
                SetPhysicalPhaseInfo,(*list.ptrModel).ptrFit,list.lambda,/O
                info=GetPhysicalPhaseInfo((*list.ptrModel).ptrFit,nphase)
                if nphase ne 0 then begin
                    table=widget_table(tab4a,value=info.data,row_labels=[info.rowname],column_labels=[info.colname],X_SCROLL_SIZE=3,Y_SCROLL_SIZE=16)
                    button=widget_button(tab4a,value='Scaling factors from weights ...',uvalue=0)
                endif
                info=GetPhysicalPhaseInfo((*list.ptrModel).ptrFit,nphase,/O)
                if nphase ne 0 then begin
                    table=widget_table(tab4b,value=info.data,row_labels=[info.rowname],column_labels=[info.colname],X_SCROLL_SIZE=3,Y_SCROLL_SIZE=16)
                    button=widget_button(tab4b,value='Scaling factors from weights ...',uvalue=1)
                endif
                
                (*(*(*list.ptrModel).ptrFit).sourceemission)[0]=list.lambda
                table=widget_table(tab0,value=*(*(*list.ptrModel).ptrFit).sourceemission,$
                    column_labels=['Wavelength('+msymbols('angstroms')+')','Intensity','FWHM'],$
                    X_SCROLL_SIZE=3,Y_SCROLL_SIZE=4,/editable,uvalue='elines',uname='elines',/context_events)
                context = WIDGET_BASE(tab0,  /CONTEXT_MENU, UNAME='ContextELines',Y_SCROLL_SIZE=2)
                t = widget_button(context,value='Add line')
                t = widget_button(context,value='Delete line')
                SourceEmission,context

                return
                endcase
    'PD':        begin
                if update then return
                str='The items here are groups of peaks with position defined in 2-theta (no crystallographic info required).'
                str2='Right-click the folder to add groups.'
                endcase
    'Pawley':    begin
                if update then return
                str='The items here are phases with known spacegroup.'
                str2='Right-click the folder to add phases.'
                endcase
    'Rietveld':    begin
                if update then return
                str='The items here are phases with known spacegroup and asymmetric unit.'
                str2='Right-click the folder to add groups of phases.'
                endcase
    else: return
endcase
label=widget_label(base,value=str,/ALIGN_LEFT)
label=widget_label(base,value=str2,/ALIGN_LEFT)
end;pro BranchFolderInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro BranchInfo,evtop,update=update,nodropupdate=nodropupdate

widget_control,evtop,update=keyword_set(update)

ID=widget_info(evtop,FIND_BY_UNAME='topbranch')
widget_control,ID,get_uvalue=list2 ; first field: pointer to group, second field: ID of selected branch item

if PTR_VALID(list2.ptr) then BranchFileInfo,list2.ptr,evtop,update=update,nodropupdate=nodropupdate $ ; file
else begin ; folder
    widget_control,list2.ID,get_value=val
    BranchFolderInfo,val,evtop,update=update
endelse

widget_control,evtop,update=1
end;pro BranchInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro IncludeModelEntree,ptr,id

if (*ptr).include then begin
    (*ptr).include=0b
    widget_control,id,SET_TREE_BITMAP=ExcludedBMP()
endif else begin
    (*ptr).include=1b
    widget_control,id,SET_TREE_BITMAP=0
endelse

end;pro IncludeModelEntree
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ExcludeAllButThis,ptr,id,ptrfit

; Structure
ptrcur=ptrfit
while PTR_VALID((*ptrcur).next) do begin
    ptrcur=(*ptrcur).next
    (*ptrcur).include=ptrcur eq ptr
endwhile

; Display
node=widget_info(widget_info(widget_info(id,/parent),/parent),/child)
repeat begin
    childID=widget_info(node,/child)
    while childID ne 0 do begin
        widget_control,childID,SET_TREE_BITMAP=(childID eq id)?0:ExcludedBMP()
        childID = Widget_Info(childID, /Sibling)
    endwhile
    node = Widget_Info(node, /Sibling)
endrep until (node eq 0)
    
end;pro ExcludeAllButThis
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro IncludeModelType,ptr,id,type,bool

childID=widget_info(id,/child)
if childID eq 0 then return

repeat begin
    if bool then widget_control,childID,SET_TREE_BITMAP=0 $
    else widget_control,childID,SET_TREE_BITMAP=ExcludedBMP()
    childID = Widget_Info(childID, /Sibling)
endrep until (childID eq 0)

while PTR_VALID((*ptr).next) do begin
    ptr=(*ptr).next
    if (*ptr).type eq type then (*ptr).include=bool
endwhile

end;pro IncludeModelType
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CreateModel

; Main structure:
return,PTR_NEW(    {n:0,$                                        ; number of groups
                    ptrSum:PTR_NEW(    {n:intarr(4),$            ; number of groups of type0, type1, type2 and type3
                                    next:PTR_NEW()}),$        ; pointer to first group of type0, type1, type2 or type3
                     ptrFit:PTR_NEW(    {n:intarr(3),$            ; number of groups of type10, type11 and type12
                                     nextra:0,$                ; total number of extra group representations (normally a group has 1 representation, but e.g. type12 has #=phasenr)
                                    code:16777473UL,$        ; bit code (32bit unsigned): |(4)|(4)|(4)|(4)|"Zero shift type"(4)|"<deleted>"(4)|"Lorentz-polarisation factor type"(4)|"calc background + substract for this region? or include in fit"(4)|
                                                            ; Lorentz-polarization correction:
                                                            ;    0: None (default)
                                                            ;    1: Bragg-Brentano
                                                            ;    2: Bragg-Brentano with graphite analyzer
                                                            ; Background type:
                                                            ;    0: None (default)
                                                            ;    1: Calc+substract (strip)
                                                            ;    2: Calc+substract (ortpol)
                                                            ;    3: Fit (ortpol)
                                                            ; Zero shift type:
                                                            ;    0: 2theta shift
                                                            ;    1: sample-detector shift (default)
                                    sourceemission:ptr_new([-1,1.,1.]),$; Additional source emission lines
                                    icurrent:0l,$            ; Current file number
                                    nparam:PTR_NEW(),$        ; number of background parameters and number of other param
                                    POL:PTRARR(2),$            ; 0: each row is a polynomial, 1: fitted ortpol coefficients
                                    backparamIO:[-1,5000l,1],$ ; background parameters    (these are not fit parameters!): ortpol degree, strip iter, strip width
                                    yfit:PTR_NEW(PTRARR(1)),$; fitted pattern
                                    yb:PTR_NEW(),$            ; background (when fixed)
                                    factors:PTR_NEW([ptr_new(fltarr(1,23))]),$    ; npts, nparam, nparamcon, npegged, nreject, nparamfixed, iter, nfree(npts-nparam), time, CHI2, CHI2red, Rp, Rwp, Re, cRp, cRwp, cRe, Rb, Rf, D, dw, Qd, D+p.ln(p), D+2p, Qd<d<4Qd
                                     fitcontrol:[50,1.e-10,1,0,1,1],$ ; maximum iterations, tollerance, refresh(0:no, 1:only once, 2:each iteration), numerical derivatives?, subtract background before fitting?, weighting scheme?
                                    misc:[ptr_new(0.2),$    ; atom radius multiplier
                                        ptr_new(2),$         ; plot ASU/UC/Expanded
                                        ptr_new(0b),$        ; show ASU labels
                                        ptr_new(0.00001),$    ; fractional coord. precision
                                        ptr_new(1b),$        ; anomalous dispersion
                                        ptr_new(3.1356),$    ; monochromator d-spacing ; Graphite = 3.355, Si(111) = 3.1356
                                        ptr_new(0.9),$        ; linear degree of polarization (horizontal)
                                        ptr_new(intarr(2,3)),$    ; number of unit cells displayed
                                        ptr_new(2),$        ; number of monochromator crystals
                                        ptr_new(1.),$        ; mosaic = 2, perfect = 1
                                        ptr_new(1),$        ; vertical = 1, horizontal = -1
                                        ptr_new(45.)],$        ; scattering azimuth in degrees
                                    debug:0b,$                ; Debug flag
                                    backranind:lonarr(2),$    ; part of yfit used
                                     next:PTR_NEW()})})        ; pointer to first group of type10, type11 and type12

end;function CreateModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DestroyModel,ptrModel,full=full

if keyword_set(full) then heap_free,ptrModel $
else begin
    heap_free,(*ptrModel).ptrSum
    heap_free,(*ptrModel).ptrFit
    p=CreateModel()
    (*ptrModel)=*p
endelse

end;pro DestroyModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SaveFitTypex,lun,ptr

printf,lun,(*ptr).group
printf,lun,(*ptr).name
printf,lun,(*ptr).include
printf,lun,(*ptr).peaknr
printf,lun,(*ptr).code

if (*ptr).peaknr ne 0 then begin
    printf,lun,float(*(*ptr).peakparamIO.I),format=format
    printf,lun,float(*(*ptr).peakparamIO.constr),format=format
    printf,lun,float(*(*ptr).globalparamIO.I),format=format
    printf,lun,float(*(*ptr).globalparamIO.constr),format=format
endif

printf,lun,(*ptr).togglefix

if (*ptr).indpropfix[0] ne -1 then $
    printf,lun,(*ptr).propfix

case (*ptr).type of
10:
11: if ptr_valid((*ptr).pstrinfo) then $
    printf,lun,(*(*ptr).pstrinfo).sghash else printf,lun,0L
12: begin
    if ptr_valid((*ptr).pstrinfo) then $
    printf,lun,(*(*ptr).pstrinfo).sghash else printf,lun,0L
    printfasu,lun,ptr
    endcase
endcase

end; pro SaveFitTypex
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LoadFitTypex,lun,ptr,type,versionmask

case type of
10: CreateModeltype10,0,p
11: CreateModeltype11,0,p
12: CreateModeltype12,0,p
endcase
InsertModelEntree,ptr,p

(*p).type=type
line=(*p).group
readf,lun,line
(*p).group=line
line=''
readf,lun,line
(*p).name=line
line=(*p).include
readf,lun,line
(*p).include=line
line=(*p).peaknr
readf,lun,line
(*p).peaknr=line
line=(*p).code
readf,lun,line
(*p).code=line

icurrent=0

if (*p).peaknr ne 0 then begin
    format='('+stringr((*p).npeakmax)+'F13)'
    line=fltarr((*p).npeakmax,(*p).peaknr)
    readf,lun,line,format=format
    (*p).peakparamIO.I=ptr_new(line)
    line*=0
    (*(*p).peakparamIO.O)[icurrent]=ptr_new(line)
    (*(*p).peakparamIO.SD)[icurrent]=ptr_new(line)
    line=fltarr((*p).npeakmax,3)
    readf,lun,line,format=format
    (*p).peakparamIO.constr=ptr_new(line)

    format='('+stringr((*p).nglobalmax)+'F13)'
    line=fltarr((*p).nglobalmax)
    readf,lun,line,format=format
    (*p).globalparamIO.I=ptr_new(line)
    line*=0
    (*(*p).globalparamIO.O)[icurrent]=ptr_new(line)
    (*(*p).globalparamIO.SD)[icurrent]=ptr_new(line)
    line=fltarr((*p).nglobalmax,3)
    readf,lun,line,format=format
    (*p).globalparamIO.constr=ptr_new(line)
endif

line=bytarr((*p).nglobalmax)
readf,lun,line
(*p).togglefix=line

if (*p).indpropfix[0] ne -1 then begin
    if CompareVersion(versionmask,'6.1.1.1') lt 0 then begin
        line=(*p).propfix
        nmax=n_elements(line)-3
        line=line[0:nmax]
        readf,lun,line
        (*p).propfix[0:nmax]=line
    endif else begin
        line=(*p).propfix
        readf,lun,line
        (*p).propfix=line
    endelse
endif

case type of
10:
11:    begin
    line=0LL
    readf,lun,line; Read as 64bit integer because Linux and Windows handle 32bit integer overflow different
    line=long(line)
    sg=SpaceGroup(line,6,error=error)
    if error ne 1 then (*p).pstrinfo=ptr_new(sg)
    endcase
12:    begin
    line=0LL
    readf,lun,line
    line=long(line)
    sg=SpaceGroup(line,6,error=error)
    if error ne 1 then (*p).pstrinfo=ptr_new(sg)
    readfasu,lun,p
    endcase
endcase

ModelTypeSetCode,p,/refresh
ptr=p
end;pro LoadFitTypex
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SaveFitModel,lun,ptrModel

; Remark: type0,type1 and type2 are never saved since they are derived from "$MASK:" and "$MASK1:"

; ----Save ptrFit----
ptrFit=(*ptrModel).ptrFit
printf,lun,(*ptrFit).n
printf,lun,(*ptrFit).nextra
printf,lun,(*ptrFit).code
printf,lun,(*ptrFit).fitcontrol

n=n_elements(*(*ptrFit).sourceemission)/3
printf,lun,n
printf,lun,*(*ptrFit).sourceemission

n=n_elements((*ptrFit).misc)
printf,lun,n
for i=0l,n-1 do printf,lun,*(*ptrFit).misc[i]

; ----Save types----
ptrCur=(*ptrFit).next

WHILE PTR_VALID(ptrCur) DO BEGIN
    printf,lun,(*ptrCur).type
    SaveFitTypex,lun,ptrCur
    ptrCur=(*ptrCur).next
ENDWHILE

end;pro SaveFitModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LoadFitModel,lun,ptrModel,versionmask,add=add

; ----Read ptrFit----
if add then begin
    ptrFit=(*ptrModel).ptrFit
    line=(*ptrFit).n
    readf,lun,line
    (*ptrFit).n+=line
    ntot=total(line)
    (*ptrModel).n+=ntot
    line=(*ptrFit).nextra
    readf,lun,line
    (*ptrFit).nextra+=line

    line=(*ptrFit).code
    readf,lun,line
    line=(*ptrFit).fitcontrol
    readf,lun,line
    
    line=0L
    readf,lun,line
    line=fltarr(3,line)
    readf,lun,line
    
    line=0L
    readf,lun,line
    for i=0l,line-1 do begin
        line=*(*ptrFit).misc[i]
        readf,lun,line
    endfor
    
    ptr=ptrFit
    while ptr_valid((*ptr).next) do ptr=(*ptr).next
endif else begin
    DestroyModel,ptrModel
    ptrFit=(*ptrModel).ptrFit

    line=(*ptrFit).n
    readf,lun,line
    (*ptrFit).n=line
    ntot=total(line)
    (*ptrModel).n=ntot
    line=(*ptrFit).nextra
    readf,lun,line
    (*ptrFit).nextra=line
    line=(*ptrFit).code
    readf,lun,line
    (*ptrFit).code=line
    line=(*ptrFit).fitcontrol
    readf,lun,line
    (*ptrFit).fitcontrol=line
    
    line=0L
    readf,lun,line
    line=fltarr(3,line)
    readf,lun,line
    *(*ptrFit).sourceemission=line
    
    line=0L
    readf,lun,line
    for i=0l,line-1 do begin
        line=*(*ptrFit).misc[i]
        readf,lun,line
        *(*ptrFit).misc[i]=line
    endfor

    ptr=ptrFit
endelse

; ----Read types----
for i=1,ntot do begin
    type=0
    readf,lun,type
    LoadFitTypex,lun,ptr,type,versionmask
endfor

end;pro LoadFitModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeleteModelEntree,ptr,id,ptrModel
ptrkeep=(*ptr).prev

; cut
(*(*ptr).prev).next=(*ptr).next
if PTR_VALID((*ptr).next) then (*(*ptr).next).prev=(*ptr).prev
(*ptr).prev=PTR_NEW()
(*ptr).next=PTR_NEW()

; free memory
type=(*ptr).type
i=type-10
heap_free,ptr

; set new group numbers
(*ptrModel).n--
ptrCur=(*ptrModel).ptrFit
(*ptrCur).n[i]--
ptrCur=(*ptrCur).next
i=0
while PTR_VALID(ptrCur) do begin
    if (*ptrCur).type eq type then (*ptrCur).group=i++
    ptrCur=(*ptrCur).next
endwhile

widget_control,id,/destroy
ptr=ptrkeep
end;pro DeleteModelEntree
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeleteModelType,ptrModel,id,type

childID=widget_info(id,/child)
if childID eq 0 then return

ptr=(*ptrModel).ptrFit
while PTR_VALID((*ptr).next) do begin
    ptr=(*ptr).next
    if (*ptr).type eq type then begin
        nextchildID = Widget_Info(childID, /Sibling)
        DeleteModelEntree,ptr,childID,ptrModel
        childID=nextchildID
    endif
endwhile

end;pro DeleteModelType
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeleteModelTypeSelect,ptrModel,id,type

childID=widget_info(id,/child)
if childID eq 0 then return

ptr=(*ptrModel).ptrFit
while PTR_VALID((*ptr).next) do begin
    ptr=(*ptr).next
    if (*ptr).type eq type then begin
        nextchildID = Widget_Info(childID, /Sibling)
        res=dialog_message('Delete '+(*ptr).name,/question,/default_no)
        if res eq 'Yes' then DeleteModelEntree,ptr,childID,ptrModel
        childID=nextchildID
    endif
endwhile

end;pro DeleteModelTypeSelect
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeleteAllButThis,ptr,id,ptrModel
type=(*ptr).type
childID=widget_info(widget_info(id,/parent),/child)

ptrcur=(*ptrModel).ptrFit
while PTR_VALID((*ptrcur).next) do begin
    ptrcur=(*ptrcur).next
    if (*ptrcur).type eq type then begin
        nextchildID = Widget_Info(childID, /Sibling)
        if ptrcur ne ptr and childID ne id then DeleteModelEntree,ptrcur,childID,ptrModel
        childID=nextchildID
    endif
endwhile
end;pro DeleteAllButThis
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro InsertModelEntree,ptrBefore,ptrInsertObj
(*ptrInsertObj).next=(*ptrBefore).next
(*ptrBefore).next=ptrInsertObj
(*ptrInsertObj).prev=ptrBefore
if PTR_VALID((*ptrInsertObj).next) then (*(*ptrInsertObj).next).prev=ptrInsertObj
end;pro InsertModelEntree
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DuplicateModel,ptr,n,p

pnext=(*ptr).next
pprev=(*ptr).prev
(*ptr).next=ptr_new()
(*ptr).prev=ptr_new()

p=ptr
Heap_Copy,p

(*ptr).next=pnext
(*ptr).prev=pprev

case (*p).type of
10:name='PD'+stringr(n)
11:name='Pawley'+stringr(n)
12:name='Rietveld'+stringr(n)
endcase
(*p).name=name

end;pro DuplicateModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ConvertModeltypex_10,psource,pdest,lambda,O=O
if (*psource).peaknr eq 0 then return

; peakparam
icurrent=0
case (*psource).type of
11:    begin
    gen=GetPawleypeakparam(psource,lambda,icurrent,O=O); tt,I
    i0=5
    i1=11
    endcase
12:    begin
    gen=GetRietveldpeakparam(psource,lambda,icurrent,O=O,/PD); tt,I
    i0=4
    i1=10
    endcase
endcase

data=fltarr((*pdest).npeakmax,(*pdest).peaknr)
data[0:1,*]=gen

p=(*psource).peakparamIO.I
data[2:8,*]=(*p)[i0:i1,*]
(*pdest).peakparamIO.I=ptr_new(data)
    
p=(*(*psource).peakparamIO.O)[icurrent]
data[2:8,*]=(*p)[i0:i1,*]
(*(*pdest).peakparamIO.O)[icurrent]=ptr_new(data)
    
p=(*(*psource).peakparamIO.SD)[icurrent]
data[0:1,*]=0
data[2:8,*]=(*p)[i0:i1,*]
(*(*pdest).peakparamIO.SD)[icurrent]=ptr_new(data)
    
data=fltarr((*pdest).npeakmax,3)
p=(*psource).peakparamIO.constr
data[2:8,*]=(*p)[i0:i1,*]
(*pdest).peakparamIO.constr=ptr_new(data)

; global param
data=fltarr((*pdest).nglobalmax)
data[2]=1

p=(*psource).globalparamIO.I
data[0:1]=(*p)[0:1]
data[3:31]=(*p)[8:36]
(*pdest).globalparamIO.I=ptr_new(data)

p=(*(*psource).globalparamIO.O)[icurrent]
data[0:1]=(*p)[0:1]
data[3:31]=(*p)[8:36]
(*(*pdest).globalparamIO.O)[icurrent]=ptr_new(data)

p=(*(*psource).globalparamIO.SD)[icurrent]
data[0:1]=(*p)[0:1]
data[2]=0
data[3:31]=(*p)[8:36]
(*(*pdest).globalparamIO.SD)[icurrent]=ptr_new(data)

data=fltarr((*pdest).nglobalmax,3)
p=(*psource).globalparamIO.constr
data[0:1,*]=(*p)[0:1,*]
data[3:31,*]=(*p)[8:36,*]
(*pdest).globalparamIO.constr=ptr_new(data)

(*pdest).togglefix[0:1]=(*psource).togglefix[0:1]
(*pdest).togglefix[3:31]=(*psource).togglefix[8:36]
    
end;pro ConvertModeltypex_10
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ConvertModeltype12_11,psource,pdest,lambda,O=O

if (*psource).peaknr ne 0 then begin
    ; peakparam
    icurrent=0
    gen=GetRietveldpeakparam(psource,lambda,icurrent,O=O); m,Fh^2

    data=fltarr((*pdest).npeakmax,(*pdest).peaknr)
    data[3:4,*]=gen
    
    p=(*psource).peakparamIO.I
    data[0:2,*]=(*p)[0:2,*]
    data[5:11,*]=(*p)[4:10,*]
    (*pdest).peakparamIO.I=ptr_new(data)
    
    p=(*(*psource).peakparamIO.O)[icurrent]
    data[0:2,*]=(*p)[0:2,*]
    data[5:11,*]=(*p)[4:10,*]
    (*(*pdest).peakparamIO.O)[icurrent]=ptr_new(data)
    
    p=(*(*psource).peakparamIO.SD)[icurrent]
    data[0:2,*]=(*p)[0:2,*]
    data[3:4,*]=0
    data[5:11,*]=(*p)[4:10,*]
    (*(*pdest).peakparamIO.SD)[icurrent]=ptr_new(data)
    
    data=fltarr((*pdest).npeakmax,3)
    p=(*psource).peakparamIO.constr
    data[0:2,*]=(*p)[0:2,*]
    data[5:11,*]=(*p)[4:10,*]
    (*pdest).peakparamIO.constr=ptr_new(data)

    ; global param
    p=(*psource).globalparamIO.I
    (*pdest).globalparamIO.I=ptr_new(*p)
    p=(*(*psource).globalparamIO.O)[icurrent]
    (*(*pdest).globalparamIO.O)[icurrent]=ptr_new(*p)
    p=(*(*psource).globalparamIO.SD)[icurrent]
    (*(*pdest).globalparamIO.SD)[icurrent]=ptr_new(*p)
    p=(*psource).globalparamIO.constr
    (*pdest).globalparamIO.constr=ptr_new(*p)
endif

(*pdest).togglefix=(*psource).togglefix
(*pdest).pstrinfo=ptr_new(*(*psource).pstrinfo)

end;pro ConvertModeltype12_11
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ConvertModeltype,psource,pdest,lambda,O=O

; Copy some general fields
icurrent=0

(*pdest).name=(*psource).name
(*pdest).include=(*psource).include
(*pdest).peaknr=(*psource).peaknr
(*pdest).code=(*psource).code
(*(*pdest).propI)=(*(*psource).propI)
(*(*(*pdest).propO)[icurrent])=(*(*(*psource).propO)[icurrent])

if (*pdest).indpropfix[0] ne -1 then $
if keyword_set(O) then (*pdest).propfix=(*(*(*pdest).propO)[icurrent])[(*pdest).indpropfix] $
else (*pdest).propfix=(*(*pdest).propI)[(*pdest).indpropfix]

; Copy type specific fields
case (*psource).type of
10: 
11:    case (*pdest).type of
    10: ConvertModeltypex_10,psource,pdest,lambda,O=O
    else:
    endcase
12: case (*pdest).type of
    10: ConvertModeltypex_10,psource,pdest,lambda,O=O
    11: ConvertModeltype12_11,psource,pdest,lambda,O=O
    else:
    endcase
endcase

end;pro ConvertModeltype
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CreateModeltype12,n,p

toggledefault=bytarr(37)
toggledefault[[2,3,4,5,6,7,9,10,12,15,16,19,21,23,24,30,31]]=1b

p=PTR_NEW({type:12,$                ; Model type
group:n,$                            ; Phase number for type12
name:'Rietveld'+stringr(n),$        ; Phase name
include:1b,$                        ; Include this phase?
peaknr:0,$                            ; Number of peaks
factors:PTR_NEW([ptr_new(fltarr(1,2))]),$; Statistical parameters: Rb(%), Rf(%)
code:[0ul,50529027ul,771ul],$        ; bit code (32bit unsigned)=>
                                    ; code[0]: profile   |(8)|(8)|"type of peak assym"(8)|"peak profile type"(8)|
                                    ; code[1:2]: fit     |(3)|"ref+H"(1)|"ref+L"(1)|"ref+PM"(1)|"ref?"(1)|"global?"(1)|
                                    ;    => |n(8)|W(8)|I(8)|u(8)|      |(8)|(8)|A(8)|R(8)|
npeakmax:11,$                        ; number of peak variables in memory (cfr. peakparamIO) for each peak: hkl(3),I(0),W(2),n(2),R(2),A(1) => ARR1
peakhelpparam:[1b,1b,1b,0b,0b,0b,0b,0b,0b,0b,0b],$     ; when 1 => this parameter is a help parameter and is never refined
nglobalmax:37,$                        ; number of global variables in memory (cfr. globalparamIO): zero(2),unit cell(6),scaling(1),W(4,2,3),n(2,2),R(3,2,2),A(3,1,4) => ARR2
nasumax:9,$                            ; number of asu variables in memory (cfr. asu) for each atom in ASU: name(1), Z(1), 'Ion'(1), xyz(3),SOF(1),B(1)
peaklabels:['h','k','l','m','FWHM_L','FWHM_G','n_l','n_r','R_l','R_r','A'],$
globallabels:['zero('+msymbols('degrees')+')','ddist(mm)','a','b','c','alpha','beta','gamma','scaling','U','V','W','IG','X','Y','UG','VG','WG','n0_l','X_l','n0_r','X_r',$
'r1','r2','r3','r1','r2','r3','r4','a1','a2','a3','P1','P1','P2','P3','P4'],$
asulabels:['Natoms','Name','Z','Wyck','x','y','z','SOF','Biso'],$
togglefix:toggledefault,$            ; e.g. UL->1: UL is fixed, even if the code says FWHM is refined (and global)
ScalingMultiplier:0l,$                ; Scaling factor multiplier (10^ScalingMultiplier) to make fit stable
codeproc:intarr(6,8),$                ; processed code for u,I,W,n,R,A
                                    ; row0: used?
                                    ; row1: global?
                                    ; row2: refined?
                                    ; row3: refined+L?
                                    ; row4: refined+H?
                                    ; row5: refined+PM?
                                    ; row6: global is (1) function of peakparameters (2) function of refined peakparameters (not in a least squares way)
                                    ; row7: refined+c(rejection)?
paramind:PTRARR(7,2),$                ; indices used for u,I,W,n,R,A,All
                                    ; row0: in ARR1
                                    ; row1: in ARR2
indasused:[4,5,6,7,8],$                ; indices of asulabels used as parameters
indasudep:[4,5,6],$                    ; indices of asulabels that might depend on other atoms in ASU
indpropfix:[-1],$                    ; fixed derived properties
propI:PTR_NEW(fltarr(16)),$            ; derived properties (initial)
propO:PTR_NEW([PTR_NEW(fltarr(16))]),$ ; derived properties (refined)
propname:['Weight frac.(w%)','E.s.d.(w%)','Volume frac.(v%)','E.s.d.(v%)','Mole frac.(n/n%)','E.s.d.(n/n%)','Rel. Mass(g/mol)','E.s.d.(g/mol)',$
            'Unit Cell Volume('+msymbols('angstroms')+msymbols('^3')+')','E.s.d.('+msymbols('angstroms')+msymbols('^3')+')',$
            'Density(g/cm'+msymbols('^3')+')','E.s.d.(g/cm'+msymbols('^3')+')','Mass Att.Coeff.(cm'+msymbols('^2')+'/g)','E.s.d.(cm'+msymbols('^2')+'/g)',$
            'Coh.Cross.(cm'+msymbols('^2')+'/g)','E.s.d.(cm'+msymbols('^2')+'/g)'],$
peakparamIO:{I:PTR_NEW(), constr:PTR_NEW(),$
        O:PTR_NEW(PTRARR(1)),SD:PTR_NEW(PTRARR(1))},$
                                    ; initial + refined/SD peak parameters
                                    ; dimensions: peaknr x npeakmax
globalparamIO:{I:PTR_NEW(), constr:PTR_NEW(),$
        O:PTR_NEW(PTRARR(1)),SD:PTR_NEW(PTRARR(1))},$
                                    ; initial + refined/SD global peak parameters
                                    ; dimensions: nglobalmax
nasupos:0L,$                        ; number of positions in the ASU
nasuatoms:0L,$                        ; number of atoms in the ASU
preverse_indices:ptr_new(),$        ; reference indices for ASU (dimension: nasupos)
asu:{I:PTR_NEW(), constr:PTR_NEW(),$        ; asymmetric unit
        O:PTR_NEW(PTRARR(1)),SD:PTR_NEW(PTRARR(1))},$    ; dimensions: nasuatoms x {...}
prev:PTR_NEW(),$                    ; points to previous model
next:PTR_NEW(),$                    ; points to next model
pstrinfo:PTR_NEW()})                ; SG:SG> {sghash,allops,...}
                                    ; celparam:celparam> A,A,A,deg,deg,deg

end;pro CreateModeltype12
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CreateModeltype11,n,p

toggledefault=bytarr(37)
toggledefault[[2,3,4,5,6,7,9,10,12,15,16,19,21,23,24,30,31]]=1b

p=PTR_NEW({type:11,$                ; Model type
group:n,$                            ; Phase number for type11
name:'Pawley'+stringr(n),$            ; Phase name
include:1b,$                        ; Include this phase?
peaknr:0,$                            ; Number of peaks
factors:PTR_NEW([ptr_new(fltarr(1,2))]),$; Statistical parameters: Rb(%), Rf(%)
code:[0ul,50529027ul,771ul],$        ; bit code (32bit unsigned)=>
                                    ; code[0]: profile   |(8)|"LeBail"(8)|"type of peak assym"(8)|"peak profile type"(8)|
                                    ; code[1:2]: fit     |(3)|"ref+H"(1)|"ref+L"(1)|"ref+PM"(1)|"ref?"(1)|"global?"(1)|
                                    ;    => |n(8)|W(8)|I(8)|u(8)|      |(8)|(8)|A(8)|R(8)|
npeakmax:12,$                        ; number of peak variables in memory (cfr. peakparamIO) for each peak: hkl(3),mI(2),W(2),n(2),R(2),A(1) => ARR1
peakhelpparam:[1b,1b,1b,1b,0b,0b,0b,0b,0b,0b,0b,0b],$ ; when 1 => this parameter is a help parameter and is never refined
nglobalmax:37,$                        ; number of global variables in memory (cfr. globalparamIO): zero(2),unit cell(6),scaling(1),W(4,2,3),n(2,2),R(3,2,2),A(3,1,4) => ARR2
peaklabels:['h','k','l','m','|Fhkl|^2','FWHM_L','FWHM_G','n_l','n_r','R_l','R_r','A'],$
globallabels:['zero('+msymbols('degrees')+')','ddist(mm)','a','b','c','alpha','beta','gamma','scaling','U','V','W','IG','X','Y','UG','VG','WG','n0_l','X_l','n0_r','X_r',$
'r1','r2','r3','r1','r2','r3','r4','a1','a2','a3','P1','P1','P2','P3','P4'],$
togglefix:toggledefault,$            ; e.g. UL->1: UL is fixed, even if the code says FWHM is refined (and global)
ScalingMultiplier:0l,$                ; Scaling factor multiplier (10^ScalingMultiplier) to make fit stable
codeproc:intarr(6,8),$                ; processed code for u,I,W,n,R,A
                                    ; row0: used?
                                    ; row1: global?
                                    ; row2: refined?
                                    ; row3: refined+L?
                                    ; row4: refined+H?
                                    ; row5: refined+PM?
                                    ; row6: global is (1) function of peakparameters (2) function of refined peakparameters (not in a least squares way)
                                    ; row7: refined+c(rejection)?
paramind:PTRARR(7,2),$                ; indices used for u,I,W,n,R,A,All
                                    ; row0: in ARR1
                                    ; row1: in ARR2
indpropfix:[6,7,12,13,14,15],$        ; fixed derived properties in prop
propfix:fltarr(6),$                    ; fixed derived properties
propI:PTR_NEW(fltarr(16)),$            ; derived properties (initial)
propO:PTR_NEW([PTR_NEW(fltarr(16))]),$            ; derived properties (refined)
propname:['Weight frac.(w%)','E.s.d.(w%)','Volume frac.(v%)','E.s.d.(v%)','Mole frac.(n/n%)','E.s.d.(n/n%)','Rel. Mass(g/mol)','E.s.d.(g/mol)',$
            'Unit Cell Volume('+msymbols('angstroms')+msymbols('^3')+')','E.s.d.('+msymbols('angstroms')+msymbols('^3')+')',$
            'Density(g/cm'+msymbols('^3')+')','E.s.d.(g/cm'+msymbols('^3')+')','Mass Att.Coeff.(cm'+msymbols('^2')+'/g)','E.s.d.(cm'+msymbols('^2')+'/g)',$
            'Coh.Cross.(cm'+msymbols('^2')+'/g)','E.s.d.(cm'+msymbols('^2')+'/g)'],$
peakparamIO:{I:PTR_NEW(), constr:PTR_NEW(),$
        O:PTR_NEW(PTRARR(1)),SD:PTR_NEW(PTRARR(1))},$
                                    ; initial + refined/SD peak parameters
                                    ; dimensions: peaknr x npeakmax
globalparamIO:{I:PTR_NEW(), constr:PTR_NEW(),$
        O:PTR_NEW(PTRARR(1)),SD:PTR_NEW(PTRARR(1))},$
                                    ; initial + refined/SD global peak parameters
                                    ; dimensions: nglobalmax
prev:PTR_NEW(),$                    ; points to previous model
next:PTR_NEW(),$                    ; points to next model
pstrinfo:PTR_NEW()})                ; SG:SG> {sghash,allops,...}
                                    ; celparam:celparam> A,A,A,deg,deg,deg

end;pro CreateModeltype11
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CreateModeltype10,n,p

toggledefault=bytarr(32)
toggledefault[[2,4,5,7,10,11,14,16,18,19,25,26]]=1b

p=PTR_NEW({type:10,$                ; Model type
group:n,$                            ; Group number for type10
name:'PD'+stringr(n),$                ; Group name
include:1b,$                        ; Include this group?
peaknr:0,$                            ; Number of peaks
factors:PTR_NEW(),$                    ; Statistical parameters
code:[0ul,33686018ul,514ul],$        ; bit code (32bit unsigned)=>
                                    ; code[0]: profile   |(8)|(8)|"type of peak assym"(8)|"peak profile type"(8)|
                                    ; code[1:2]: fit     |(3)|"ref+H"(1)|"ref+L"(1)|"ref+PM"(1)|"ref?"(1)|"global?"(1)|
                                    ;    => |n(8)|W(8)|I(8)|u(8)|      |(8)|(8)|A(8)|R(8)|
npeakmax:9,$                        ; number of peak variables in memory (cfr. peakparamIO) for each peak: u(1),I(1),W(2),n(2),R(2),A(1) => ARR1
peakhelpparam:bytarr(9),$             ; when 1 => this parameter is a help parameter and is never refined
nglobalmax:32,$                        ; number of global variables in memory (cfr. globalparamIO): zero shift(2),isotropic deformation(1),scaling,W(4,2,3),n(2,2),R(3,2,2),A(3,1,4) => ARR2
peaklabels:['2-theta','I','FWHM_L','FWHM_G','n_l','n_r','R_l','R_r','A'],$
globallabels:['zero('+msymbols('degrees')+')','ddist(mm)','isodeform','scaling','U','V','W','IG','X','Y','UG','VG','WG','n0_l','X_l','n0_r','X_r',$
'r1','r2','r3','r1','r2','r3','r4','a1','a2','a3','P1','P1','P2','P3','P4'],$
togglefix:toggledefault,$            ; e.g. UL->1: UL is fixed, even if the code says FWHM is refined (and global)
ScalingMultiplier:0l,$                ; Scaling factor multiplier (10^ScalingMultiplier) to make fit stable
codeproc:intarr(6,8),$                ; processed code for u,I,W,n,R,A
                                    ; row0: used?
                                    ; row1: global?
                                    ; row2: refined?
                                    ; row3: refined+L?
                                    ; row4: refined+H?
                                    ; row5: refined+PM?
                                    ; row6: global is (1) function of peakparameters (2) function of refined peakparameters (not in a least squares way)
                                    ; row7: refined+c(rejection)?
paramind:PTRARR(7,2),$                ; indices used for u,I,W,n,R,A,All
                                    ; row0: in ARR1
                                    ; row1: in ARR2
indpropfix:[6,7,8,9,12,13,14,15],$    ; fixed derived properties in prop
propfix:fltarr(8),$                    ; fixed derived properties
propI:PTR_NEW(fltarr(16)),$            ; derived properties (initial)
propO:PTR_NEW([PTR_NEW(fltarr(16))]),$            ; derived properties (refined)
propname:['Weight frac.(w%)','E.s.d.(w%)','Volume frac.(v%)','E.s.d.(v%)','Mole frac.(n/n%)','E.s.d.(n/n%)','Rel. Mass(g/mol)','E.s.d.(g/mol)',$
            'Unit Cell Volume('+msymbols('angstroms')+msymbols('^3')+')','E.s.d.('+msymbols('angstroms')+msymbols('^3')+')',$
            'Density(g/cm'+msymbols('^3')+')','E.s.d.(g/cm'+msymbols('^3')+')','Mass Att.Coeff.(cm'+msymbols('^2')+'/g)','E.s.d.(cm'+msymbols('^2')+'/g)',$
            'Coh.Cross.(cm'+msymbols('^2')+'/g)','E.s.d.(cm'+msymbols('^2')+'/g)'],$
peakparamIO:{I:PTR_NEW(), constr:PTR_NEW(),$
        O:PTR_NEW(PTRARR(1)),SD:PTR_NEW(PTRARR(1))},$
                                    ; initial + refined/SD peak parameters
                                    ; dimensions: peaknr x npeakmax
                                    ; O and SD are arrays of pointers (batch processing)
globalparamIO:{I:PTR_NEW(), constr:PTR_NEW(),$
        O:PTR_NEW(PTRARR(1)),SD:PTR_NEW(PTRARR(1))},$
                                    ; initial + refined/SD global peak parameters
                                    ; dimensions: nglobalmax
prev:PTR_NEW(),$                    ; points to previous model
next:PTR_NEW()})                    ; points to next model

end;pro CreateModeltype10
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MoveModelGroup,ev
; ev.id: reference node (destination)
; wDraggedNode: dragged node (source)
wParent = WIDGET_INFO( ev.id, /PARENT )
wDraggedNode = WIDGET_INFO( ev.drag_id, /TREE_SELECT )

; Make sure source and destination are in the same folder
if WIDGET_INFO( wDraggedNode, /PARENT ) ne wParent then return

; Get indices of source and reference node
index = WIDGET_INFO( ev.id, /TREE_INDEX )
IF ( ev.position EQ 4 ) THEN index++
curindex=WIDGET_INFO( wDraggedNode, /TREE_INDEX )
if curindex lt index then index--
if index eq curindex then return
widget_control,wDraggedNode,SET_TREE_INDEX=index

; Move model
widget_control,wDraggedNode,get_uvalue=ptrThis
widget_control,ev.id,get_uvalue=ptrRef
if ptr_valid((*ptrThis).next) then (*(*ptrThis).next).prev=(*ptrThis).prev
if ptr_valid((*ptrThis).prev) then (*(*ptrThis).prev).next=(*ptrThis).next
IF ( ev.position EQ 4 ) then GenericListInsertAfter,ptrRef,ptrThis $
else GenericListInsertBefore,ptrRef,ptrThis

end;pro MoveModelGroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddModeltype,type,ptrModel,parentbranch,pduplicate=pduplicate,convert=convert

if PTR_VALID(ptrModel) eq 0 then return

; ----Browse for insertion position----
ptrCur=(*ptrmodel).ptrFit
n=(*ptrCur).n[type-10]
ptrNext=(*ptrCur).next
while PTR_VALID(ptrNext) do begin
    if (*ptrNext).type gt type then ptrNext=0 $
    else begin
        ptrCur=ptrNext
        ptrNext=(*ptrCur).next
    endelse
endwhile

; ----Insert----
if keyword_set(pduplicate) then begin
    DuplicateModel,pduplicate,n,p
endif else begin
    case type of
        10: CreateModeltype10,n,p
        11: CreateModeltype11,n,p
        12: CreateModeltype12,n,p
    endcase

    if keyword_set(convert) then ConvertModeltype,convert.ptr,p,convert.lambda,O=convert.O
endelse
InsertModelEntree,ptrCur,p
ModelTypeSetCode,p,/refresh

; ----Add count----
(*ptrModel).n++
ptr=(*ptrmodel).ptrFit
(*ptr).n[type-10]++

branch = WIDGET_TREE(parentbranch,value=(*p).name,uvalue=p,/drop_events,/draggable)
widget_control,branch,/SET_TREE_VISIBLE

if (*p).include eq 0 then $
    widget_control,branch,SET_TREE_BITMAP=ExcludedBMP()

end;pro AddModeltype
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_EditModel,ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=top
if WIDGET_INFO(top,/VALID_ID) then begin
    widget_control,top,get_uvalue=list
    list.NewWindowID[0]=0
    widget_control,top,set_uvalue=list
endif
end;CleanUp_EditModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditModel_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
0 : begin ;base
    ; resizing event
    flplot_obj_event,ev
    endcase
1 : begin ;button event
    widget_control,ev.id,get_value=val
    widget_control,ev.id,get_uvalue=uval
    widget_control,ev.top,get_uvalue=top
    widget_control,top,get_uvalue=list
    ID=widget_info(ev.id,/parent)
    widget_control,ID,get_uvalue=list2

    case val of
    'Add':    AddModeltype,list2.ptr,list.ptrModel,list2.ID
    'Duplicate': AddModeltype,(*list2.ptr).type,list.ptrModel,widget_info(list2.ID,/parent),pduplicate=list2.ptr
    'Copy to PD':     begin
                    O=DIALOG_MESSAGE('Use refined parameters to calculate properties?',/question,/DEFAULT_NO) eq 'Yes'
                    keep=(*list2.ptr).include
                    (*list2.ptr).include=1
                    SetPhysicalPhaseInfo,(*list.ptrModel).ptrFit,list.lambda,O=O
                    (*list2.ptr).include=keep
                    
                    branchID=widget_info(widgetparent(list2.ID,nup=2),find_by_uname='PD')
                    AddModeltype,10,list.ptrModel,branchID,convert={ptr:list2.ptr,O:O,lambda:list.lambda}
                    endcase
    'Scattering yield': begin
                    param = list.scatyieldstruct
                    ElasticScatteringFractionWidget,top,list2.ptr,param
                    list.scatyieldstruct = param
                    widget_control,top,set_uvalue=list
                    endcase
    'Copy to Pawley': begin
                    O=DIALOG_MESSAGE('Use refined parameters to calculate properties?',/question,/DEFAULT_NO) eq 'Yes'
                    keep=(*list2.ptr).include
                    (*list2.ptr).include=1
                    SetPhysicalPhaseInfo,(*list.ptrModel).ptrFit,list.lambda,O=O
                    (*list2.ptr).include=keep
                    
                    branchID=widget_info(widgetparent(list2.ID,nup=2),find_by_uname='Pawley')
                    AddModeltype,11,list.ptrModel,branchID,convert={ptr:list2.ptr,O:O,lambda:list.lambda}
                    endcase
    'Include All': IncludeModelType,(*list.ptrModel).ptrFit,list2.ID,list2.ptr,1b
    'Exclude All': IncludeModelType,(*list.ptrModel).ptrFit,list2.ID,list2.ptr,0b
    'Refined = Initial': RefinedToInitialType,(*list.ptrModel).ptrFit,list2.ptr,ev,/inverse
    'Initial = Refined': RefinedToInitialType,(*list.ptrModel).ptrFit,list2.ptr,ev
    'Delete All':    begin
                    ;ID=list2.ID
                    DeleteModelType,list.ptrModel,list2.ID,list2.ptr
                    ;widget_control,ID,SET_TREE_EXPANDED=0
                    ;widget_control,ID,send_event={ID:ID,TOP:ev.top,HANDLER:ev.handler, TYPE:1, EXPAND:0L}
                    endcase
    'Delete':    begin
                ID=widget_info(list2.ID,/parent)
                DeleteModelEntree,list2.ptr,list2.ID,list.ptrModel
                ;widget_control,ID,SET_TREE_EXPANDED=0
                widget_control,ID,send_event={ID:ID,TOP:ev.top,HANDLER:ev.handler, TYPE:1, EXPAND:0L}
                endcase
    'Delete ...':begin
                DeleteModelTypeSelect,list.ptrModel,list2.ID,list2.ptr
                endcase
    'Change Name':begin
            ChangeModelEntreeName,list2.ptr,list2.ID,ev.top
            BranchInfo,ev.top
            endcase
    'Include/Exclude':begin
            IncludeModelEntree,list2.ptr,list2.ID
            BranchInfo,ev.top,/update
            endcase
    'Exclude all but this':begin
            ExcludeAllButThis,list2.ptr,list2.ID,(*list.ptrModel).ptrFit
            BranchInfo,ev.top,/update
            endcase
    'Delete all but this':begin
            DeleteAllButThis,list2.ptr,list2.ID,list.ptrModel
            BranchInfo,ev.top,/update
            endcase
    else:
    endcase
    SendMainRefresh,ev.top
    endcase
11:    begin ;tree
    widget_control,ev.id,get_uvalue=uval
    event=TAG_NAMES(ev, /STRUCTURE_NAME)

    case event of
    'WIDGET_CONTEXT':begin
                    ; Show context menu for this group or type or main
                    if PTR_VALID(uval.ptr) then ID=WIDGET_INFO(ev.top, FIND_BY_UNAME = 'ContextGroup'+stringr((*uval.ptr).type))$
                    else begin
                        if uval.ptr eq -1 then return
                        ID=WIDGET_INFO(ev.top, FIND_BY_UNAME = 'ContextType')
                    endelse
                    widget_control,ID,set_uvalue=uval
                    WIDGET_DISPLAYCONTEXTMENU, ev.ID, ev.X, ev.Y, ID
                    endcase
    'WIDGET_DROP':    MoveModelGroup,ev
    else:             begin
                    ; Show info for this group or type or main
                    topbranch=widget_info(ev.id,/TREE_ROOT)
                    widget_control,topbranch,get_uvalue=list2
                    if list2.ID eq ev.id then return
                    widget_control,topbranch,set_uvalue={ptr:uval,ID:ev.id}
                    BranchInfo,ev.top
                    endelse
    endcase

    endcase
endcase
end;pro EditModel_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddModelEntree,ptri,parentbranch,type

ptr=ptri
while PTR_VALID(ptr) do begin
    if (*ptr).type eq type then begin
        branch = WIDGET_TREE(parentbranch,value=(*ptr).name,uvalue=ptr,/drop_events,/draggable)
        if (*ptr).include then widget_control,branch,SET_TREE_BITMAP=0 $
        else widget_control,branch,SET_TREE_BITMAP=ExcludedBMP()
    endif
    ptr=(*ptr).next
endwhile

end;pro AddModelEntree
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EditModel,top

base=widget_base(/row,title='Pattern Modelling', uvalue=top,/TLB_SIZE_EVENTS)
widget_control,top,get_uvalue=list
list.NewWindowID[0]=base
widget_control,top,set_uvalue=list
ptrfirst=(*(*list.ptrModel).ptrFit).next

topbranch = WIDGET_TREE(base,/CONTEXT_EVENTS,uvalue={ptr:0L,ID:0L},uname='topbranch')
Mbranch = WIDGET_TREE(topbranch,value='Model', /FOLDER, uvalue=-1L)
branch = WIDGET_TREE(Mbranch,value='PD', /FOLDER, uvalue=10L, uname='PD')
    AddModelEntree,ptrfirst,branch,10
branch = WIDGET_TREE(Mbranch,value='Pawley', /FOLDER, uvalue=11L, uname='Pawley')
    AddModelEntree,ptrfirst,branch,11
branch = WIDGET_TREE(Mbranch,value='Rietveld', /FOLDER, uvalue=12L, uname='Rietveld')
    AddModelEntree,ptrfirst,branch,12

context = WIDGET_BASE(base,  /CONTEXT_MENU, UNAME='ContextType')
t = widget_button(context,value='Add')
t = widget_button(context,value='Delete All')
t = widget_button(context,value='Delete ...')
t = widget_button(context,value='Include All')
t = widget_button(context,value='Exclude All')
t = widget_button(context,value='Refined = Initial')
t = widget_button(context,value='Initial = Refined')
context = WIDGET_BASE(base,  /CONTEXT_MENU, UNAME='ContextGroup10')
t = widget_button(context,value='Include/Exclude')
t = widget_button(context,value='Exclude all but this')
t = widget_button(context,value='Change Name')
t = widget_button(context,value='Duplicate')
t = widget_button(context,value='Delete')
t = widget_button(context,value='Delete all but this')
t = widget_button(context,value='Refined = Initial')
t = widget_button(context,value='Initial = Refined')
context = WIDGET_BASE(base,  /CONTEXT_MENU, UNAME='ContextGroup11')
t = widget_button(context,value='Include/Exclude')
t = widget_button(context,value='Exclude all but this')
t = widget_button(context,value='Change Name')
t = widget_button(context,value='Delete')
t = widget_button(context,value='Delete all but this')
t = widget_button(context,value='Duplicate')
t = widget_button(context,value='Copy to PD')
t = widget_button(context,value='Refined = Initial')
t = widget_button(context,value='Initial = Refined')
context = WIDGET_BASE(base,  /CONTEXT_MENU, UNAME='ContextGroup12')
t = widget_button(context,value='Include/Exclude')
t = widget_button(context,value='Exclude all but this')
t = widget_button(context,value='Change Name')
t = widget_button(context,value='Delete')
t = widget_button(context,value='Delete all but this')
t = widget_button(context,value='Duplicate')
t = widget_button(context,value='Copy to Pawley')
t = widget_button(context,value='Copy to PD')
t = widget_button(context,value='Scattering yield')
t = widget_button(context,value='Refined = Initial')
t = widget_button(context,value='Initial = Refined')


info = widget_base(base,uname='groupinfoparent')

WIDGET_CONTROL, base, /REALIZE
widget_control,branch,/SET_TREE_VISIBLE

Xmanager,'EditModel',base,/NO_BLOCK,cleanup='CleanUp_EditModel',GROUP_LEADER=top
end;pro EditModel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
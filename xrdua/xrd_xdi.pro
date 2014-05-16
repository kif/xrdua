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

function getxdicoord,list,coord

case list.ScanDim of
1:    begin
    ; coord must be in default orientation: moveh,LR,LR(+)
    ; convert to current orientation: TRAN=0,REV1,Hsign or TRAN=1,REV2,Vsign
    ; = [move0,move1]
    coord=coord[*,0]
    if list.TRAN then brev=list.Vsign lt 0 xor list.REV2 else brev=list.Hsign lt 0 xor list.REV1
    if brev then coord=coord[reverse(sort(coord))] else coord=coord[sort(coord)]
    endcase
2:    begin
    ; coord must be in default orientation: moveh,LR,TB,LR(+),TB(+)
    ; convert to current orientation: REV1,REV2,Hsign,Vsign
    ; = [[moveh0,moveh1],[movev0,movev1]]
    if list.Hsign lt 0 xor list.REV1 then coord[*,0]=coord[reverse(sort(coord[*,0])),0] else coord[*,0]=coord[sort(coord[*,0]),0]
    if list.Vsign lt 0 xor list.REV2 then coord[*,1]=coord[reverse(sort(coord[*,1])),1] else coord[*,1]=coord[sort(coord[*,1]),1]
    endcase
endcase

return,coord
end;function getxdicoord
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CorrectICR,ICR,OCR,DT,plot=plot,outid=outid

; Deadtime regions:
; color regions(only for display): 0--red--thresLOWdt--green--thresHIGHdt--blue--100
; linear region for DT as a function of ICR: 0--thresLIN
DT<=100
thresLOWdt=30
thresHIGHdt=90
thresLIN=50

; ----Get linear relation ICR as a function of DT ----
h=histogram(DT,bin=1,rev=rev)
indh=where(h ne 0,cth)

; Loop over all DT values
for j=0l,cth-1 do begin
    i=indh[j]
    k=rev[rev[i]:rev[i+1]-1]
    DTi=DT[k[0]]
    if DTi le thresLIN then begin
        ; This deadtime is in the linear region
        if n_elements(mean0) eq 0 then mean0=mean(ICR[k])

        tmp1=where(ICR[k] gt mean0,tmp2, comp=tmp3, ncomp=tmp4)
        if tmp2 ne 0 then if n_elements(indLINFIT) eq 0 then indLINFIT=k[tmp1] else indLINFIT=[indLINFIT,k[tmp1]]
        if tmp4 ne 0 then if n_elements(indLOWER) eq 0 then indLOWER=k[tmp3] else indLOWER=[indLOWER,k[tmp3]]
    endif
endfor
if n_elements(indLINFIT) eq 0 then return

; Linear fit to DT(ICR) in the linear region
res = LADFIT(DT[indLINFIT],ICR[indLINFIT])

if keyword_set(plot) then begin
    pkeep=!P
    old_dev = !d.name
    old_fancy = !fancy
    
    a=PSWINDOW()
    set_plot, 'PS'
    !fancy = 0
    device, file=plot,/COLOR,/ENCAPSULATED, /isolatin1
    loadct,40,/silent

    col0=16*15;'FF'XL
    col1=16*9;'FF00'XL
    col2=16*4;'FF0000'XL
    ;!P.BACKGROUND=0
    ;!P.COLOR='FFFFFF'XL
    !P.MULTI = [0, 2, 2]
    
    indLOW=where(DT le thresLOWdt,ct0,comp=indMID,ncomp=ct1)
    indHIGH=where(DT ge thresHIGHdt,ct2)

    plot,ICR,OCR,/nodata,xtitle='ICR(cps)',ytitle='OCR(cps)',title='OCR vs. ICR';,/iso
    if ct0 ne 0 then oplot,ICR[indLOW],OCR[indLOW],psym=2,col=col0
    if ct1 ne 0 then oplot,ICR[indMID],OCR[indMID],psym=2,col=col1
    if ct2 ne 0 then oplot,ICR[indHIGH],OCR[indHIGH],psym=2,col=col2

    plot,DT,ICR,/nodata,xtitle='DT(%)',ytitle='ICR(cps)',title='ICR vs. DT',yrange=[0,max(ICR)*1.5];,/iso
    if ct0 ne 0 then oplot,DT[indLOW],ICR[indLOW],psym=2,col=col0
    if ct1 ne 0 then oplot,DT[indMID],ICR[indMID],psym=2,col=col1
    if ct2 ne 0 then oplot,DT[indHIGH],ICR[indHIGH],psym=2,col=col2
    if n_elements(indLOWER) ne 0 then oplot,DT[indLOWER],ICR[indLOWER],psym=2

    p=[min(DT),max(DT)]
    oplot,p,res[0]+res[1]*p
endif

; ----Correct for non-linearity----
for j=0l,cth-1 do begin
    i=indh[j]
    k=rev[rev[i]:rev[i+1]-1]
    DTi=DT[k[0]]
    if DTi le thresLIN then begin
        ;tmp1=where(ICR[k] gt mean(ICR[k]),tmp2)
        ;if tmp2 gt 0 then ICR[k[tmp1]]=(res[0]+res[1]*DTi)
    endif else begin
        ; Make linear in the non-linear region
        ICR[k]=(res[0]+res[1]*DTi)
    endelse
endfor

if keyword_set(plot) then begin
    plot,ICR,OCR,/nodata,xtitle='ICR(cps)',ytitle='OCR(cps)',title='OCR vs. corrected ICR';,/iso
    if ct0 ne 0 then oplot,ICR[indLOW],OCR[indLOW],psym=2,col=col0
    if ct1 ne 0 then oplot,ICR[indMID],OCR[indMID],psym=2,col=col1
    if ct2 ne 0 then oplot,ICR[indHIGH],OCR[indHIGH],psym=2,col=col2

    plot,DT,OCR,/nodata,xtitle='DT(%)',ytitle='OCR(cps)',title='OCR vs. DT';,/iso
    if ct0 ne 0 then oplot,DT[indLOW],OCR[indLOW],psym=2,col=col0
    if ct1 ne 0 then oplot,DT[indMID],OCR[indMID],psym=2,col=col1
    if ct2 ne 0 then oplot,DT[indHIGH],OCR[indHIGH],psym=2,col=col2

    device,/close
    !fancy = old_fancy
    set_plot, old_dev
    !P=pkeep
    printw,outid,'Saved '+plot
endif
end;pro CorrectICR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ParseXiast,path,outid=outid,prompt=prompt

CATCH, Error_status
IF Error_status NE 0 THEN begin
    printw,outid,!error_state.MSG
    printw,outid,'ParseXIA failed.'
    return,0b
endif

; Select file(s) of the scan
if keyword_set(prompt) then $
    sfile=dialog_pickfile(path=path,filter='*.edf',/multi) $
else sfile=file_search(path+path_sep()+'*.edf')
if sfile[0] eq '' then return,0
sfile=sfile[sort(sFile)]

for i=0l,n_elements(sfile)-1 do begin
    file=sfile[i]
    
    tmp=CutPath(file,path=path,file=file)
    
    ; Parse filename
    ; radix_xiast_0009_0000_0004: scannumber=9, rownumber=4
    sp=strsplit(file,'_',/ext,count=ct)
    if ct lt 5 then begin
        printw,outid,'ParseXIA failed (unknown filename format).'
        return,0b
    endif
    if ct gt 5 then begin
        ex=ct-5
        for j=1,ex do sp[0]+='_'+sp[j]
        sp=sp[[0,indgen(4)+ex+1]]
    endif
    sp[1]='xiast'
    mask=path+sp[0]+'_'+sp[1]+'_'+sp[2]+'_'+sp[3]+'_*.edf'
    file=path+sp[0]+'_'+sp[1]+'_'+sp[2]+'_'+sp[3]+'_'+sp[4]+'.edf'
    
    ; Get number of detectors and number of columns in the scan
    if OPENR_safe(lun, file) then begin
        printw,outid,"ParseXIA failed (can't open edf file)."
        return,0b
    endif
    line='dummy'
    while strmid(line,strlen(line)-1,1) ne '}' do begin
        readf,lun,line
        if strmid(line,0,5) eq 'Dim_1' then begin
            p1=strpos(line,'=', /REVERSE_SEARCH)+1
            p2=strpos(line,';', /REVERSE_SEARCH)
            ndet=fix(strmid(line,p1,p2-p1))/6
        endif
        if strmid(line,0,5) eq 'Dim_2' then begin
            p1=strpos(line,'=', /REVERSE_SEARCH)+1
            p2=strpos(line,';', /REVERSE_SEARCH)
            ncols=fix(strmid(line,p1,p2-p1))
        endif
    endwhile
    free_lun,lun
    if n_elements(ndet) eq 0 or n_elements(ncols) eq 0 then return,0 
    if i eq 0 then begin
        ndetkeep=ndet
        ncolskeep=ncols
    endif else begin
        if ndet ne ndetkeep or ncols ne ncolskeep then return,0
    endelse
    
    ; Get all xiast files of this scan
    tmp=file_search(mask,count=ct)
    if ct ne 0 then $
        if i eq 0 then begin
            files=tmp
            nrows=ct
        endif else begin
            files=[files,tmp]
            nrows+=ct
        endelse

endfor

; Read ICR,OCR,DT
ICR=fltarr(ndet,ncols,nrows)
OCR=ICR
DT=ICR
EVT=ICR
LiveT=ICR
off=lindgen(ndet)*6
for j=0l,nrows-1 do begin
    xiast=readedf(files[j])
    EVT[*,*,j]=xiast[1+off,*]
    ICR[*,*,j]=xiast[2+off,*]
    OCR[*,*,j]=xiast[3+off,*]
    LiveT[*,*,j]=xiast[4+off,*]
    DT[*,*,j]=xiast[5+off,*]
endfor
ind=where(ICR eq 0,ct)
if ct ne 0 then begin
    OCR[ind]=1
    ICR[ind]=1
endif

; Output filenames
pathout=path+'COR'+path_sep()
file_mkdir,pathout

for k=0l,ndet-1 do begin
    base=pathout+sp[0]+'_'+'xia'+string(k,format='(I02)')+'_'+sp[2]+'_'
    xiapath1=base+'cor.xdi'
    xiapath2=base+'corcor.xdi'
    xiapath3=base+'corcor_info.eps'
    xiapath4=base+'stat.xdi'

    ICR_=ICR[k,*,*]
    OCR_=OCR[k,*,*]
    DT_=DT[k,*,*]
    
    ndim=(nrows gt 1) + 1
    
    Type3XDI,xiapath1,OCR_/ICR_,'OCR/ICR',ndim,outid
    Type3XDI,xiapath4,[EVT[k,*,*],ICR_,OCR_,LiveT[k,*,*],DT_],$
    ['Events','ICR(cts/sec)','OCR(cts/sec)','LiveTime(ms)','DT(%)'],ndim,outid
    
    CorrectICR,ICR_,OCR_,DT_,plot=xiapath3,outid=outid
    
    ind2=where(ICR_ eq 0,ct2)
    if ct2 ne 0 then begin
        OCR_[ind2]=1
        ICR_[ind2]=1
    endif
    Type3XDI,xiapath2,OCR_/ICR_,'OCR/ICRcorrected',ndim,outid
endfor

return,1b
end;function ParseXiast
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDIdrawfree,ID
widget_control,ID,get_uvalue=uval
if n_elements(uval) ne 0 then heap_free,uval
end;pro XDIdrawfree,ID
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GetXDIdatadim,list,nmaps,ncol,nrow,default=default
nmaps=list.NMaps
if keyword_set(default) and list.TRAN then BEGIN
    ncol=list.MapRow
    nrow=list.MapCol
endif else begin
    nrow=list.MapRow
    ncol=list.MapCol
endelse
end;pro GetXDIdatadim
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetGroupLabel,list,i
labeli=0>list.labeli<(n_elements(((*list.pmstruc).(i+1)).data)-1)
return,stringr(((*list.pmstruc).(i+1)).data[labeli])
end;function GetGroupLabel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RGBtriangle,ev
; Make RGB Chromaticity Triangle

widget_control,ev.top,get_uvalue=list
file=SelFile(list.path,'RGBtriangle.tif','*.tif','Save RGB Chromaticity Triangle...')
if file eq '' then return
xsize=long(promptnumber(stringr(list.xsize),ev.top,'Image size for RGB triangle'))
ysize=xsize*sqrt(3)/2

x=rebin(lindgen(xsize)/(xsize-1.),xsize,ysize,/sample)
y=rebin(lindgen(1,ysize)/(xsize-1.),xsize,ysize,/sample)

A=2*y/sqrt(3)   ; top
B=1-y/sqrt(3)-x ; left
C=x-y/sqrt(3)   ; right

background=255
ind=where(B lt 0 or C lt 0)
A[ind]=background
B[ind]=background
C[ind]=background

;Maximum possible total intensity
M=A>B>C
A/=M
B/=M
C/=M

img=bytscl(transpose([[[A]],[[B]],[[C]]],[2,0,1]))

write_tiff_safe,file,img,ORIENTATION=3
end;pro RGBtriangle
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function alpha_blending,list,indices=indices,nextra=nextra,ncombi=ncombi,error=error

; Get indices
error=1b
indices=-1L
n=n_elements(*list.pgradients)
iuse=-1L

for i=0l,n-1 do begin
    if ~(*(*list.pgradients)[i]).buse then continue
    index=(*(*list.pgradients)[i]).index
    if index lt 0 or index ge list.NMaps then continue
    indices=[indices,index]
    iuse=[iuse,i]
endfor

ncombi=n_elements(indices)-1
if ncombi le 0 then return,0
indices=indices[1:*]
iuse=iuse[1:*]

; Blend
for i=0l,ncombi-1 do begin
    image=reform((*list.data)[indices[i],*,*],list.MapCol,list.MapRow)
    AlphaBlendingWithGradient,R,G,B,A,image,(*(*list.pgradients)[iuse[i]]),rgbback=rgbback
endfor
GradientFillMissing,R,G,B,rgbback,list.gradientfill
;ConvertBlackToWhite,r,g,b,0.1

; Return
R=long(0>round(R)<255)
G=long(0>round(G)<255)
B=long(0>round(B)<255)
img=R or ishft(G,8) or ishft(B,16)

nextra=1
error=0b
return,reform(img,1,list.MapCol,list.MapRow,/overwrite)
end;function alpha_blending
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro autofindpixels,image,plist,preplacep,cutoff,direction,operator,top,lofind=lofind,plistall=plistall,preplacepall=preplacepall,replacevalkeep=replacevalkeep

ball=keyword_set(plistall) and keyword_set(preplacepall)

if n_elements(*plist) ne 0 then image[(*plist)[0,*]]=image[(*plist)[1,*]]
case operator of
'==':ind=where(image eq cutoff,ct)
'!=':ind=where(image ne cutoff,ct)
'<':ind=where(image lt cutoff,ct)
'>':ind=where(image gt cutoff,ct)
'>=':ind=where(image ge cutoff,ct)
else: ind=where(image le cutoff,ct);'<='
endcase

if ct ne 0 and direction eq -1 then begin
    if n_elements(replacevalkeep) eq 0 then begin
        val=float(mean(image[ind]))
        tmp=PromptNumber(stringr(val),top,'Enter replacing value:')
        reads,tmp,val
        replacevalkeep=val
    endif else val=replacevalkeep
    tmp=[reform(ind,1,ct),replicate(val,1,ct)]
    
    if n_elements(*preplacep) eq 0 then *preplacep=tmp $
    else *preplacep=[[*preplacep],[tmp]]
    
    if ball then for i=0l,n_elements(*preplacepall)-1 do begin
        p=(*preplacepall)[i]
        if p eq preplacep then continue
        if n_elements(*p) eq 0 then *p=tmp $
        else *p=[[*p],[tmp]]
    endfor

    return
endif

if ct ne 0 and direction eq -2 then begin
    ; Make LOF
    lofind=temporary(ind)
    return
endif

if (ct ne 0) and (ct ne n_elements(image)) then begin
    list=lonarr(2,ct)
    list[0,*]=ind
    s=size(image)
    x=ind mod s[1]
    y=ind/s[1]
    xorig=x
    yorig=y

    ind2=lindgen(ct)
    case direction of
    0:add0=-1
    1:add0=1
    2:add0=-1
    3:add0=1
    endcase
    add=replicate(add0,ct)
    ctold=ct

    while ct ne 0 do begin
        case direction of
        0:    begin
            y[ind2]+=add[ind2]
            ind=where(y[ind2] ge s[2],ct)
            if ct ne 0 then begin
                y[ind2[ind]]+=2
                add[ind2[ind]]=1
            endif
            endcase
        1:    begin
            y[ind2]+=add[ind2]
            ind=where(y[ind2] lt 0,ct)
            if ct ne 0 then begin
                y[ind2[ind]]-=2
                add[ind2[ind]]=-1
            endif
            endcase
        2:    begin
            x[ind2]+=add[ind2]
            ind=where(x[ind2] ge s[1],ct)
            if ct ne 0 then begin
                x[ind2[ind]]+=2
                add[ind2[ind]]=1
            endif
            endcase
        3:    begin
            x[ind2]+=add[ind2]
            ind=where(x[ind2] lt 0,ct)
            if ct ne 0 then begin
                x[ind2[ind]]-=2
                add[ind2[ind]]=-1
            endif
            endcase
        endcase

        case operator of
        '==':ind2=where(image[x,y] eq cutoff[x,y],ct)
        '!=':ind2=where(image[x,y] ne cutoff[x,y],ct)
        '<':ind2=where(image[x,y] lt cutoff[x,y],ct)
        '>':ind2=where(image[x,y] gt cutoff[x,y],ct)
        '>=':ind2=where(image[x,y] ge cutoff[x,y],ct)
        else: ind2=where(image[x,y] le cutoff[x,y],ct);'<='
        endcase

        if ct eq ctold then ct=0 else ctold=ct
    endwhile

    list[1,*]=x+y*s[1]
    if ptr_valid(plist) then (*plist)=[[(*plist)],[list]] else plist=ptr_new(list)
    
    if n_elements(*plist) eq 0 then *plist=list $
    else *plist=[[*plist],[list]]
    
    if ball then for i=0l,n_elements(*plistall)-1 do begin
        p=(*plistall)[i]
        if p eq plist then continue
        if n_elements(*p) eq 0 then *p=list $
        else *p=[[*p],[list]]
    endfor
    
endif
end;pro autofindpixels
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckSubCommand,str
char=byte(str)
return, total(((char ge 48) and (char le 57)) or (char eq 64)) eq n_elements(char)
end;function CheckSubCommand
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ParseSubCommand,subcmdsplit,add,ct,badd

badd=1b
add=''
case ct of
0:    badd=0b
1:     begin
    bool = subcmdsplit[0] eq '@'
    badd = ~bool
    bool or= CheckSubCommand(subcmdsplit[0])
    add=subcmdsplit[0]
    ct=0
    endcase
2:     begin
    bool = (subcmdsplit[0] eq 'not')
    bool and= CheckSubCommand(subcmdsplit[1])
    add=subcmdsplit[0]+' '+subcmdsplit[1]
    ct=0
    endcase
else:begin

    if (subcmdsplit[0] eq 'not') then begin
        bool = CheckSubCommand(subcmdsplit[1])

        add=subcmdsplit[0]+' '+subcmdsplit[1]

        subcmdsplit=['@',subcmdsplit[2:*]]
        ct-=1
    endif else if (subcmdsplit[2] eq 'not') then begin
        bool = ct gt 3
        if bool then begin
            bool and= CheckSubCommand(subcmdsplit[3])

            add=subcmdsplit[2]+' '+subcmdsplit[3]

            if ct gt 4 then subcmdsplit=[subcmdsplit[0:1],'@',subcmdsplit[4:*]] $
            else subcmdsplit=[subcmdsplit[0:1],'@']
            ct-=1
        endif
    endif else begin
        bool = (subcmdsplit[1] eq 'and') or (subcmdsplit[1] eq 'or') or (subcmdsplit[1] eq 'xor')

        bool and= CheckSubCommand(subcmdsplit[0])
        bool and= CheckSubCommand(subcmdsplit[2])

        add=subcmdsplit[0]+' '+subcmdsplit[1]+' '+subcmdsplit[2]

        ct-=3
        if ct gt 0 then begin
            subcmdsplit=['@',subcmdsplit[3:*]]
            ct++
        endif

    endelse

    endcase
endcase

badd and= bool
return,bool
end; function ParseSubCommand
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ParseComboCommand,str,commands=subcmds,indices=indices,error=error

; Definitions:
; D1. Separator (Sc): <space> character
; D2. Grouper (Gc): Gc1='(', Gc2=')' characters
; D3. Block: a piece of the command, separated from the rest by and Sc or Gc
; D4. Number (Nr): block with only characters 0-9
; D5. Operator (Op): block with only characters a-z
; D6. @: special character used in parsing
; D7. Command (Cmd): a sequence of all previous defined elements
; D8. SubCommand (SubCmd): X AND Y, X OR Y, X XOR Y, NOT Y  (X and Y = Nr or @)
;
; Command Syntax rules:
; S1. There must be something else then only Sc or Gc
; S2. The last Block of a Command must be a Number
; S3. The first Block of a Command must be a Number or the 'not'-Operator
; S4. Blocks can contain only Numbers or Operators
; S5. No Blocks containing a Number, may be positioned next to eachother
; S6. No @ signs allowed in a Command
; S7. Total number of Gc1 must equal that of Gc2 in a Command
; S8. No more Gc2's then Gc1's at any given position in a Command
;
; Command Parse rules:
; P1. Search in Command for the first Gc2 and look back for the corresponding Gc1
; P2. Extract what is between and replace it with @. This is the rest.
; P3. Split the extracted piece with Sc in units.
; P4. Number of units 0: don't do anything
;                        1: unit must be equal to Nr or '@', if Nr: save as subcommand
;                      2: first unit must be 'not' and the second must be Nr or @, if so: save as subcommand
;                      3 or more: - If the first unit is 'not', the second must be Nr or @,
;                                    if so: save subcommand and replace units with '@'
;                                 - If the third unit is 'not', the forth must be Nr or @,
;                                    if so: save subcommand and replace units with '@'
;                                 - Else: the second subcommand must be 'and' or 'or' and the first and
;                                    thrid must be Nr or @, if so: save subcommand and replace units with '@'
; P5. Repeat P1-P4 on the Rest until no more Gc's are present
; P6. P3-P4 on the Rest

error=''

; Preprocess string
str=strtrim(strcompress(strlowcase(str)),2)

; Split string
stra=strsplit(str,' ()',/extract,COUNT=nstra)
bool=nstra ne 0

if ~bool then begin
    error='Syntax error: not known command.'
    return,0
endif
;------------------------Syntax OK: S1

; Number of chars
char=byte(stra)
sum=total(char ne 0,1) ; number of chars on each row
tmp=total((char ge 48) and (char le 57),1) ; number of 0-9 chars on each row

; S2: end with a number
bool = (tmp[nstra-1] ne 0)

; S3: begin with a number unless the first is 'not'
bool and= (tmp[0] ne 0) or (stra[0] eq 'not')

tmp1=tmp eq sum ; rows with only numbers
tmp2=tmp eq 0 ; rows with only operators

; S4: no number/operator mix
bool and= total(tmp1 or tmp2) eq nstra

; S5: no numbers next to eachother
tmp1=where(tmp1,ct)
if ct ne 0 then indices=stra[tmp1]
if (ct gt 1) then begin
    bool and= total((tmp1 eq shift(tmp1+1,1))[1:*]) eq 0
    bool and= total((tmp1 eq shift(tmp1-1,-1))[1:*]) eq 0
endif

; S6: no @ signs
bool and= strpos(str,'@') eq -1

if ~bool then begin
    error='Syntax error: check order of numbers and operators and spaces.'
    return,0
endif
;------------------------Syntax OK: S1-S6

; check ( and ) sequence
char=byte(str)
n=n_elements(char)
bopen=fix(char eq 40)
bclose=fix(char eq 41)
cum=total(bopen-bclose,/cumulative)

; S7: same number of ( and )
bool and= cum[n-1] eq 0

; S8: never more ) then ( at any given position in the string
tmp=where(cum lt 0,ct)
bool and= ct eq 0

if ~bool then begin
    error="Syntax error: check parentheses."
    return,0
endif
;------------------------Syntax OK: S7-S8

; Split command in subcommands with following syntax:
;     X and Y
;     X or Y
;     not Y
;    Y
; where X and Y are indices or @

subcmds=''
rest=str

; P1: find Gc2
pc=strpos(rest,')')
bparse=(pc ne -1) and bool

while bparse do begin ; extract all pieces between parentheses
    ; P1: find corresponding Gc1
    po=strpos(strmid(rest,0,pc),'(',/REVERSE_SEARCH)
    n=pc-po+1

    if n gt 2 then begin
        ; P2: extract piece between parentheses
        subcmd=strmid(rest,po+1,n-2)
        subcmd=strtrim(strcompress(subcmd),2)

        ; P3: split with Sc
        subcmdsplit=strsplit(subcmd,' ',/extract,count=ct)
        if ct eq 0 then bool and= 0b

        ; P4: parse units into Subcommands
        while bool and (ct gt 0) do begin
            bool and= ParseSubCommand(subcmdsplit,add,ct,badd)
            if badd then subcmds=[subcmds,add]
        endwhile
    endif

    ; P2: replace piece and parentheses with @
    repl='@'
    if n gt 1 then repl+=string(replicate(32b,n-1))
    strput,rest,repl,po

    ; P5: repeat
    pc=strpos(rest,')')
    bparse=(pc ne -1) and bool
endwhile

; P6: parse the rest
if bool then begin
    subcmdsplit=strsplit(rest,' ',/extract,count=ct)
    if ct eq 0 then bool and= 0b

    while bool and (ct gt 0) do begin
        bool and= ParseSubCommand(subcmdsplit,add,ct,badd)
        if badd then subcmds=[subcmds,add]
    endwhile

    if n_elements(subcmds) gt 1 then subcmds=subcmds[1:*]
endif

if ~bool then begin
    error="Syntax error: allowed operations are 'X and Y', 'X or Y', 'X xor Y' and 'not X'."
    return,0
endif
;------------------------Command split in subcommands

return,bool
end;function ParseComboCommand
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ParseIndices,indices,more=more

; Chop in pieces
indices=strcompress(indices,/remove_all)
indices=STRSPLIT(indices,',',/EXTRACT)
tmp=where(indices ne '',ct)
if ct ne 0 then indices=indices[tmp]

; Make integers + expand ranges
pos=strpos(indices,'-')
ind=where(pos ne -1,ct,COMPLEMENT=indc,NCOMPLEMENT=ctc)
if ct ne 0 then begin
    keep=indices[ind]
    if ctc ne 0 then indices=fix(indices[indc]) else delvar2,indices
    for i=0l,ct-1 do begin
        tmp=fix(STRSPLIT(keep[i],'-',/EXTRACT,COUNT=ct2))
        if ct2 eq 2 then begin
            b=tmp[0]<tmp[1]
            e=tmp[0]>tmp[1]
            tmp=b+indgen(e-b+1)
            if n_elements(indices) eq 0 then indices=tmp else indices=[indices,tmp]
        endif
    endfor
endif else begin
    indices=long(indices)
endelse

;indices=indices[sort(indices)]
return,(size(indices,/type) eq 7) or (keyword_set(more)?(n_elements(indices) eq 1):0b)
end;function ParseIndices
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ImageOverlap,list,cmds,indices,thres,lthen,RGB,nextra=nextra,ncombi=ncombi,error=error

; Returns TrueColor image (true=3)
error=1b

; No combi==2 allowed
nindices=n_elements(indices)
for i=0l,nindices-1 do $
    if (*list.pmstruc).(indices[i]+1).combi eq 2 then begin
        printw,list.OutID,'Image Combination error: no overlapped images allowed.'
        return,0
    endif

; Make colors (must be as long as indices and thres)
if n_elements(thres) lt nindices then begin
    printw,list.OutID,'Image Overlap error: not enough thresholds specified.'
    return,0
endif

r=reform(RGB[0,*])
g=reform(RGB[1,*])
b=reform(RGB[2,*])

if n_elements(r) lt nindices then begin
    printw,list.OutID,'Image Overlap error: not enough colors specified.'
    return,0
endif

; Interpret commands
p=ptrarr(2)
indp=0
pX=ptr_new()
pY=ptr_new()

for i=0l,n_elements(cmds)-1 do begin
    ; cmd of form: ((X) op) Y
    cmd=strsplit(cmds[i],' ',count=ct,/extract)
    case ct of
    1:    begin; Y
        strX=''
        strY=cmd[0]
        cmd=0
        endcase
    2:    begin; NOT Y
        strX=''
        strY=cmd[1]
        cmd=3
        endcase
    3:    begin; X AND/OR/XOR Y
        strX=cmd[0]
        strY=cmd[2]
        case cmd[1] of
        'and':    cmd=1
        'or':    cmd=2
        'xor':    cmd=4
        endcase
        endcase
    endcase
    ; cmd==0 -> <nothing>
    ; cmd==1 -> AND
    ; cmd==2 -> OR
    ; cmd==3 -> NOT
    ; cmd==4 -> XOR
    iX=(strX eq '@')
    iY=(strY eq '@')

    ; Image Y
    if iY eq 0 then begin
        ind=long(strY)
        j=(where(indices eq ind,ct))[0]
        if (j eq -1) or (ind ge list.NMaps) then begin
            heap_free,p
            return,0
        endif
        if lthen[j] then pY=ptr_new(reform((*list.dataorg)[ind,*,*],list.MapCol,list.MapRow) lt thres[j]) $
        else pY=ptr_new(reform((*list.dataorg)[ind,*,*],list.MapCol,list.MapRow) gt thres[j])

        if cmd eq 3 then (*pY) = ~ (*pY)
        (*pY)=[ [[*pY]], [[(*pY)*r[j]]], [[(*pY)*g[j]]], [[(*pY)*b[j]]], [[*pY]] ]

        (*pY)=long((*pY))
    endif else begin
        indp--
        pY=p[indp]
    endelse

    ; Image X
    if strX ne '' then $
    if iX eq 0 then begin
        ind=long(strX)
        j=(where(indices eq ind,ct))[0]
        if (j eq -1) or (ind ge list.NMaps) then begin
            heap_free,p
            return,0
        endif
        if lthen[j] then pX=ptr_new(reform((*list.dataorg)[ind,*,*],list.MapCol,list.MapRow) lt thres[j]) $
        else pX=ptr_new(reform((*list.dataorg)[ind,*,*],list.MapCol,list.MapRow) gt thres[j])

        (*pX)=[ [[*pX]], [[(*pX)*r[j]]], [[(*pX)*g[j]]], [[(*pX)*b[j]]], [[*pX]] ]
        (*pX)=long((*pX))
    endif else begin
        indp--
        pX=p[indp]
    endelse

    ; Evaluate: blend colors + keep binary image
    ; pY and pX: binary image, red, green and blue channel image, color multiplicity
    case cmd of
    0:    ;nothing to do
    1:    begin
        (*pY)[*,*,0] and= (*pX)[*,*,0]
        (*pY)[*,*,1:4] += (*pX)[*,*,1:4]
        endcase
    2:    begin
        (*pY)[*,*,0] or= (*pX)[*,*,0]
        (*pY)[*,*,1:4] += (*pX)[*,*,1:4]
        endcase
    3:    ;already done
    4:    begin
        (*pY)[*,*,0] xor= (*pX)[*,*,0]
        (*pY)[*,*,1:4] += (*pX)[*,*,1:4]
        endcase
    endcase
    ptr_free,pX

    if pY ne p[indp] then begin
        ptr_free,p[indp]
        p[indp]=pY
    endif
    indp++
endfor

img=*p[0]
heap_free,p

; Collapse into 2D lonarr
frac=img[*,*,0]*float(img[*,*,4])
ind=where(frac eq 0,ctwhite)
if ctwhite ne 0 then frac[ind]=1

img[*,*,1]/=frac
img[*,*,2]/=frac
img[*,*,3]/=frac

if ctwhite ne 0 then begin
    s=size(img)
    indcol=ind mod s[1]
    indrow=ind/s[1]
    ind3=replicate(3,ctwhite)

    img[indcol,indrow,ind3--]=255
    img[indcol,indrow,ind3--]=255
    img[indcol,indrow,ind3--]=255
endif

img=img[*,*,1] or ishft(img[*,*,2],8) or ishft(img[*,*,3],16)
; Extract
; img=[[[img]],[[ishft(img,-8)]],[[ishft(img,-16)]]] and 255b

ncombi=nindices
nextra=1
error=0b
return,reform(img,1,list.MapCol,list.MapRow)

end;function ImageOverlap
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ImageCombine,list,indices,nextra=nextra,ncombi=ncombi,error=error,pairs=pairs

error=1b

; No combi==2 allowed
nindices=n_elements(indices)
for i=0l,nindices-1 do $
    if (*list.pmstruc).(indices[i]+1).combi eq 2 then begin
        printw,list.OutID,'Image Combination error: no overlapped images allowed.'
        return,0
    endif

if keyword_set(pairs) then begin
    ; Sum or superimpose pairs, normalized or not
    ncombi=2
    nextra=nindices/2
    if nextra le 0 then return,0
    
    ; normalize
    n=replicate(1.,2,nextra)
    if list.CombineNorm then begin
        for i=0l,nextra-1 do begin
            tmp1=max((*list.dataorg)[indices[0,i],*,*])
            tmp2=max((*list.dataorg)[indices[1,i],*,*])
            n[*,i]=(tmp1>tmp2)/[tmp1,tmp2]
        endfor
    endif
    
    ; operation
    block=make_array(nextra,list.MapCol,list.MapRow,type=size(*list.dataorg,/type))
    case list.CombineType of
      0: for i=0l,nextra-1 do block[i,*,*]=(*list.dataorg)[indices[0,i],*,*]*n[0,i]+(*list.dataorg)[indices[1,i],*,*]*n[1,i]
      1: for i=0l,nextra-1 do block[i,*,*]=(*list.dataorg)[indices[0,i],*,*]*n[0,i]>(*list.dataorg)[indices[1,i],*,*]*n[1,i]
      2: for i=0l,nextra-1 do block[i,*,*]=(*list.dataorg)[indices[0,i],*,*]*n[0,i]/((*list.dataorg)[indices[1,i],*,*]*n[1,i])
      3: for i=0l,nextra-1 do block[i,*,*]=(*list.dataorg)[indices[0,i],*,*]*n[0,i]-(*list.dataorg)[indices[1,i],*,*]*n[1,i]
      4: for i=0l,nextra-1 do block[i,*,*]=(*list.dataorg)[indices[0,i],*,*]*n[0,i]*(*list.dataorg)[indices[1,i],*,*]*n[1,i]
    endcase
    
endif else begin
    ; normalize
    n=replicate(1.,nindices)
    if list.CombineNorm then begin
        for i=0l,nindices-1 do begin
            tmp=max((*list.dataorg)[indices[i],*,*])
            n[i]=tmp
        endfor
        n=max(n)/n
    endif
    
    if list.CombineType eq 2 or list.CombineType eq 3 then begin
        ; Combine with last one
        ncombi=2
        nextra=nindices-1
        if nextra le 0 then return,0
    
        block=(*list.dataorg)[indices[nindices-1],*,*]*n[nindices-1]
        if list.CombineType eq 2 then $
            for i=0l,nindices-2 do $
                block=[block,((*list.dataorg)[indices[i],*,*]*n[i])/block[0,*,*]] $
        else $
            for i=0l,nindices-2 do $
                block=[block,((*list.dataorg)[indices[i],*,*]*n[i])-block[0,*,*]]
        block= block[1:*,*,*]
        
        if nextra gt 1 then $
            indices=[reform(indices[0:nextra-1],1,nextra),replicate(indices[nextra],1,nextra)]
    
    endif else begin
        ; Combine all
        ncombi=nindices
        nextra=1
    
        block=(*list.dataorg)[indices[0],*,*]*n[0]
        case list.CombineType of
          0: for i=1,nindices-1 do block+=(*list.dataorg)[indices[i],*,*]*n[i]
          1: for i=1,nindices-1 do block>=(*list.dataorg)[indices[i],*,*]*n[i]
          4: for i=1,nindices-1 do block*=(*list.dataorg)[indices[i],*,*]*n[i]
        endcase
   
    endelse
endelse

error=0b
return,block
end; function ImageCombine
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro injcorrect,image,ic,n,zero,meanpix=meanpix

;ic: cut before ic's
;n: cut n before ic's

meanpix=keyword_set(meanpix)
npix=n_elements(image)
npixmax=npix
lastpixmax=npixmax-1

s=size(image)
if s[0] eq 2 then image=reform(image,npixmax) $
else meanpix=0b
bresize=0

; some restrictions
ic >= 0
ic <= npixmax
tmp=ic + n
tmp=where((tmp lt 0) or (tmp gt npixmax),ct)
if ct ne 0 then return

nn=n_elements(n)-1
for i=0l,nn do begin
    if n[i] ne 0 then begin
        lastpix=npix-1

        ; Cut pixels
        if n[i] lt 0 then begin
            image[ic[i]+n[i]:lastpix+n[i]]=image[ic[i]:lastpix]

        ; Insert pixels
        endif else begin
            overflow=image[npix-n[i]:lastpix]
            image[ic[i]+n[i]+1:lastpix]=image[ic[i]+1:lastpix-n[i]]
            meanpix=0
            if meanpix then begin
                ind=ic[i]+1+indgen(n[i])
                rep=replicate(image[ic[i]],n[i])
                ind2=ind-s[0]
                ind3=ind+s[0]
                ind4=where(ind2 gt 0 and ind3 le lastpix,ct)
                if ct ne 0 then rep[ind4]=(image[ind2[ind4]]+image[ind3[ind4]])/2.

                image[ind]=rep
            endif else image[ic[i]+1:ic[i]+n[i]]=image[ic[i]]
            image=[image,overflow]
            bresize=1b
        endelse

        if i ne nn then ic[i+1:*]+=n[i]
        npix+=n[i]
    endif
endfor

if bresize then image=image[0:npixmax-1]
if npix lt npixmax then image[npix:*]=zero
if s[0] eq 2 then image=reform(image,s[1],s[2])
end;pro injcorrect
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SinAngles,list,bsinovert=bsinovert,nviews=nviews,coordtrans=coordtrans,transunit=transunit

irot=where(list.unit eq 'deg',ct)
if ct ne 1 then begin
    nviews=list.MapRow
    bsinovert=1b
    irot=1
endif else begin
    irot=irot[0]
    bsinovert=irot eq 1
    nviews=bsinovert?list.MapRow:list.MapCol
endelse

coord=getxdicoord(list,list.coord)

; Translation range
itrans = 1-irot
coordtrans=coord[*,itrans]
transunit=list.unit[itrans]

; Angle range
c1=coord[0,irot]
c2=coord[1,irot]
cc1=(c1<c2)*!pi/180
cc2=(c1>c2)*!pi/180
add=-list.addangle
if cc1 eq cc2 then begin
    ; default: from 0 to 180-ainc degrees
    ;ainc=!pi/nviews
    ; default: from 0 to 180 degrees
    ainc=!pi/(nviews-1)
endif else begin
    ; e.g. nviews=90 and cc1=0 then the range will be from 0deg to 178deg
    add+=cc1
    ainc=(cc2-cc1)/(nviews-1)
endelse

if bsinovert then s=list.Vsign else s=list.Hsign
angles=add+ainc*indgen(nviews)
if s lt 0 then return,reverse(angles) else return,angles

end;function SinAngles
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function projcenter,list,angles,Nrho,drho

; Coordinates of a point with polar coordinates [B,C] in the tomogram:
;   x = xcen + B.cos(angles+C) + R(angles).cos(angles+C)
;   y = ycen + B.sin(angles+C) + R(angles).sin(angles+C)
;   
;   Radial distance of ellipse:
;   R(X) = D.E/sqrt(D^2.sin(X-F)^2+E^2.cos(X-F)^2)
;
; For the center of rotation B=0 and C=anything
;   x = xcen + R(angles).cos(angles+C)
;   y = ycen + R(angles).sin(angles+C)
;
; IDL's Radon projects the tomogram on the X-axis when angle = 0rad
;   projcen = xcen
;   projcen2 = xcen + R(angles).cos(angles+C)

if list.Corb[1] then projcen=list.Cor[1] else projcen=(Nrho-1)/2.
return,drho*projcen
end;function projcenter
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro tomoinfo,list,angles=angles,Nx=Nx,Ny=Ny,xmin=xmin,ymin=ymin,Nrho=Nrho,drho=drho,$
    nviews=nviews,projcen=projcen,bsinovert=bsinovert,$
    coordtrans=coordtrans,transunit=transunit
; Nx, Ny: tomogram size
; xmin,ymin: coordinates of pixel [0,0]
; Nrho: number of ray sums
; drho: ratio between "spacing between the ray sums" and "tomogram pixel size"
; nviews: number of projection

angles=SinAngles(list,bsinovert=bsinovert,nviews=nviews,coordtrans=coordtrans,transunit=transunit)

; Sinogram dimensions and spacing
;drho = list.tomodrho ; Gives artifacts in the ART algorithm due to ???: therefore divide tomogram by tomodrho instead
drho = 1
Nrho = bsinovert?list.MapCol:list.MapRow

; Tomogram dimensions (spacing = 1)
n = drho*(Nrho-1)+1
Nx=round(n*list.mXSizeTomo)
Ny=round(n*list.mYSizeTomo)

; Center of projection in sinogram
projcen=projcenter(list,angles,Nrho,drho)

; Coordinates of first pixel in tomogram (so that the center of rotation has coordinates [0,0] in the tomogram)
xmin=-projcen[0]-(((Nx-1)-drho*(Nrho-1))>0)/2.-list.mXShiftTomo
ymin=-(Ny-1)/2.-list.mYShiftTomo

; Tomogram coordinate range in transunits
mi=coordtrans[0]<coordtrans[1]
ma=coordtrans[0]>coordtrans[1]
dtrans=(ma-mi)/(Nrho-1.)/drho
bx=xmin
by=ymin
ex=bx+(Nx-1)
ey=by+(Ny-1)
coordtrans=[[bx,ex],[by,ey]]*dtrans

end;pro tomoinfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TomoSize,list,X=X,Y=Y
tomoinfo,list,Nx=Nx,Ny=Ny,bsinovert=bsinovert

if keyword_set(Y) then return,1>round(Ny/list.sfv)
return,1>round(Nx/list.sfh)
end;function TomoSize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RhoFromTomogram,evx,evy,list,nviews=nviews,bsinovert=bsinovert

tomoinfo,list,angles=angles,xmin=xmin,ymin=ymin,Ny=Ny,projcen=projcen,nviews=nviews,bsinovert=bsinovert,drho=drho

x=StatTR(evx,list.sfh)+xmin
y=Ny-1-StatTR(evy,list.sfv)+ymin  ; Flip Y-axis because Y points down (i.e. tvscl,/order)
rho=x*cos(angles)+y*sin(angles)+projcen

return,rho/drho ; Sinogram(rho(mm),angles(rad))
end;function RhoFromTomogram
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ddistCheck,ev,pixelsize=pixelsize
; TODO: update tomo recon stuff

;pixelsize= ...(mm)

; Initialize
widget_control,ev.top,get_uvalue=list
if not keyword_set(pixelsize) then pixelsize=1. ; (mm)
pixelsize=pixelsize[0]

label='scaling'
for index_scaling=list.nselect+1,list.NMaps do $
    if string(((*list.pmstruc).(index_scaling)).data[1]) eq label then break
if index_scaling gt list.NMaps then return

label=string(((*list.pmstruc).(list.nselect+1)).data[1])+'_e.s.d.'
for index_error=list.nselect+1,list.NMaps do $
    if string(((*list.pmstruc).(index_error)).data[1]) eq label then break
if index_error gt list.NMaps then index_error=-1

; Coordinates of a point with polar coordinates [B,C] in the tomogram:
;   x = xcen + B.cos(angles+C)
;   y = ycen + B.sin(angles+C)
;   Wobble is included in the sinogram but shifting rows

; Mouse at point P in tomogram (X pointing right, Y pointing down)
; Get coordinates of P at angles[0]
tomoinfo,list,angles=angles,xmin=xmin,ymin=ymin,Ny=Ny,projcen=projcen
xx=StatTR(ev.x,list.sfh)+xmin
yy=Ny-1-StatTR(ev.y,list.sfh)+ymin  ; tvscl,/order (flip Y axis from device to data coord. system)

; Get coordinates of P for all angles (in pixels)
x=xx*cos(angles)+yy*sin(angles)-xmin
y=-xx*sin(angles)+yy*cos(angles)-ymin

; Observed sample-detector distance
case list.detectorpos of
0:     begin
    ; Detector on the top
    ddist_calc=y
    endcase
1:    begin
    ; Detector on the bottom
    ddist_calc=-y
    endcase
2:    begin
    ; Detector on the left
    ddist_calc=x
    endcase
3:    begin
    ; Detector on the right
    ddist_calc=-x
    endcase
endcase

; Fitted sample-detector distance
weights=x
XX=0>[[floor(weights)],[ceil(weights)]]<(list.MapCol-1)
weights-=XX[*,0]
weights>=0
weights=[[weights],[1-weights]]
YY=indgen(list.MapRow)
YY=[[YY],[YY]]
nr=ConvertFileNumber(XX,YY,list.MapCol,list.MapRow,list.TRAN,list.REV1,list.REV2,list.REV3)

ddist_obs=(*list.dataorg)[list.nselect,*,*]
ddist_obs=(ddist_obs[nr[*,0]]+ddist_obs[nr[*,1]])/2.

if index_error ne -1 then begin
    ddist_obsesd=(*list.dataorg)[index_error,*,*] ; Standard deviation
    ddist_obsesd=sqrt(ddist_obsesd[nr[*,0]]^2.+ddist_obsesd[nr[*,1]]^2.)/2.
endif else ddist_obsesd=ddist_obs*0

; Convert from pixels to mm and get the calculated sample-detector distance offset
ddist_calc*=pixelsize
ddist_calc+=ddist_obs[0]-ddist_calc[0]

; Plot fitted sample-detector distance
mi=min(ddist_calc)<min(ddist_obs)
ma=max(ddist_calc)>max(ddist_obs)
loadctrev,39
plot,angles*180/!pi,ddist_calc,yrange=[mi,ma],/xs,/ys,xtitle='Rotation angle('+msymbols('degrees')+')',ytitle=msymbols('Delta')+'sample detector distance(mm)'
oploterror,angles*180/!pi,ddist_obs,ddist_obsesd*0,ddist_obsesd,color=100,errcolor=100

if !d.name ne 'PS' then begin
    ; Plot sinus
    wset,list.drawindex
    Device, Copy = [0,0, list.XSize,list.YSize, 0,0,list.pixindex]
    plots,RTStat(x,list.sfh),RTStat(reverse(indgen(list.MapRow)),list.sfv),/device
    
;    ; Plot [x,y]
;    Window, /Free, xsize=TomoSize(list,/X), ysize=TomoSize(list,/Y), /Pixmap
;    pixindex=!d.window
;    Device, Copy = [0,0, TomoSize(list,/X),TomoSize(list,/Y), 0,0,list.TomoReconIndex]
;    wset,list.TomoReconIndex
;    xx=RTStat(x,list.sfh)
;    yy=RTStat(list.MapCol-1-y,list.sfv)
;    plots,xx,yy,/device
;    wait,0.5
;    Device, Copy = [0,0, TomoSize(list,/X),TomoSize(list,/Y), 0,0,pixindex]
;    wdelete,pixindex
endif
end;pro ddistCheck
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ddistCheck_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
CursorPropTomo,widget_info(ev.top,FIND_BY_UNAME='position'),ev.x,ev.y,list,XX=XX

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]

wset,0

j=where(list.unit eq 'mm',ct)
if ct ne 1 then return
n=(j eq 0)?list.MapCol:list.MapRow
coord=getxdicoord(list,list.coord)
pixelsize=((coord[0,j]>coord[1,j])-(coord[0,j]<coord[1,j]))/(n-1.)
IF thisEvent eq 'DOWN' THEN begin
    list.detectorpos=PromptNumber(['Top','Bottom','Left','Right'],ev.top,'Detector in tomogram:',/choice)
    widget_control,ev.top,set_uvalue=list
    pixelsize=float(PromptNumber(stringr(pixelsize),ev.top,'Pixel size in mm:'))
    SaveImage,list.path,list.file+'_ddist.eps','ddistCheck',ev,pixelsize=pixelsize
endif else begin
    ddistCheck,ev,pixelsize=pixelsize ;(mm)
endelse

end;pro ddistCheck_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DRAWSINUS,rho,nviews,bsinovert,list

DEVICE, GET_GRAPHICS = oldg, SET_GRAPHICS=6          ;Xor drawing mode
col = 0;*16
;while col lt !d.table_size do col = col + col
;if (col ge !d.table_size) then col = !d.table_size - 1.

if bsinovert then plots,RTStat(rho,list.sfh),RTStat(reverse(indgen(nviews)),list.sfv),COLOR=col,/device,THICK=1 $
else plots,RTStat(indgen(nviews),list.sfh),list.YSize-1-RTStat(rho,list.sfv),COLOR=col,/device,THICK=1

DEVICE, SET_GRAPHICS=oldg          ;Copy mode

end;pro DRAWSINUS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DRAWSINUS_EVENTS,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
CursorPropTomo,widget_info(ev.top,FIND_BY_UNAME='position'),ev.x,ev.y,list,XX=XX

possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
bPress = possibleButtons[ev.press]
if bPress eq 'RIGHT' then begin
    modeXDI, ev, 'Dummy'
    return
endif

widget_control,ev.top,get_uvalue=list
rho=RhoFromTomogram(ev.x,ev.y,list,nviews=nviews,bsinovert=bsinovert)

wset,list.drawindex
Device, Copy = [0,0, list.XSize,list.YSize, 0,0,list.pixindex]
DRAWSINUS,rho,nviews,bsinovert,list
end;pro DRAWSINUS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReconRamlak, x, pointSpacing
    zero = where(x EQ 0.0, count)
    q = x
    if (count NE 0) then q[zero] = .01
    y = -SIN(!pi*x/2)^2 / (!pi^2 * q^2 * pointSpacing)
    if (count NE 0) then y[zero] = 1./(4.*pointSpacing)
    RETURN, y
end;function ReconRamlak
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReconFilter,type,kernelSize,NProj,NAngles,fourier=fourier,npad=npad

fourier=type eq 5
if ~fourier then begin
    ; Filter in direct space
    kernelSize=0>kernelSize<((NProj-2)/2)
    
    x = FINDGEN(2*kernelSize+1)-kernelSize ;x coordinates
    pointSpacing=1. ;point spacing
    npad=1
endif else npad=3

case type of
1:    begin ;RAMLAK
     y=ReconRamlak(x,pointSpacing)
    endcase
2:    begin ;Shepp_Logan
    pointSpacing = !pi^2 * pointSpacing * (1.-4.*x^2)
    zeros = where(abs(pointSpacing) LE 1.0e-6, count)
    if (count ne 0) then pointSpacing[zeros] = .001
    y= 2./pointSpacing
    endcase
3:    begin ;Lp_Cosine
    y= 0.5 * (ReconRamlak(x-1, pointSpacing) + ReconRamlak(x+1, pointSpacing))
    endcase
4:    begin ;Gen_Hamming
    alpha = 0.5
    y= alpha * ReconRamlak(x, pointSpacing) + ((1.-alpha)/2) * $
        (ReconRamlak(x-1, pointSpacing) + ReconRamlak(x+1, pointSpacing))
    endcase
5:    begin ; FTRam (Ramp filter: |freq| in fourier domain)
    N=NProj*npad
    drho=1
    nfreq=N/2+1 ; # positive frequencies
    freq=findgen(nfreq)/(N*drho)
    nfreq_m=nfreq-1-(~(N mod 2)) ; # negative frequencies
    freq=[freq,reverse(freq[1:nfreq_m])]
    return,COMPLEX(freq,0.0)
    endcase
else:begin ;None
    return,[1.0]
    endelse
endcase

;y/=(TOTAL(y) * NAngles * SQRT(NProj)) ; Normalize it
return,y
end;function ReconFilter
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PixelJuggle,image,list,zero,meanpix=meanpix,imageindex=imageindex
if n_elements(imageindex) ne 1 then imageindex=0

; ----Corrections on current orientation----
if ptr_valid(list.rowshift) then begin
    s=size(*list.rowshift,/dim)
    simg=size(image,/dim)
    for i=0l,s[0]-1 do begin
        j=(*list.rowshift)[i,0]
        nj=(*list.rowshift)[i,1]
        if j lt 0 then begin
            ; column shift
            j=-j
            if j lt simg[0] and abs(nj) lt simg[1] then begin
                b=image[j,0]
                e=image[j,simg[1]-1]
                f=shift(image[j,*],nj)
                if nj gt 0 then f[0:nj-1]=b
                if nj lt 0 then f[simg[1]+nj:simg[1]-1]=e
                image[j,*]=f
            endif
        endif else begin
            ; row shift
            if j lt simg[1] and abs(nj) lt simg[0] then begin
                b=image[0,j]
                e=image[simg[0]-1,j]
                f=shift(image[*,j],nj)
                if nj gt 0 then f[0:nj-1]=b
                if nj lt 0 then f[simg[0]+nj:simg[0]-1]=e
                image[*,j]=f
            endif
        endelse
    endfor
endif

; ----Corrections on default orientation----
image=reform(image,1,list.MapCol,list.MapRow,/overwrite)
ReformXDIdatablock,image,[list.TRAN,0b],[list.REV1,0b],[list.REV2,0b],[list.REV3,0b],0b,s0,s1,s2
image=reform(image,s1,s2,/overwrite)
    
if ptr_valid(list.cutpos) then $
    injcorrect,image,*list.cutpos,*list.cutn,zero,meanpix=meanpix

if ptr_valid(list.pautofill) then begin
    p=(*list.pautofill)[imageindex]
    if n_elements(*p) ne 0 then image[(*p)[0,*]]=image[(*p)[1,*]]
endif

if ptr_valid(list.preplacep) then begin
    p=(*list.preplacep)[imageindex]
    if n_elements(*p) ne 0 and size(zero,/type) ne 7 then image[(*p)[0,*]]=(*p)[1,*]
endif

image=reform(image,s0,s1,s2,/overwrite)
ReformXDIdatablock,image,[0b,list.TRAN],[0b,list.REV1],[0b,list.REV2],[0b,list.REV3],0b
image=reform(image,list.MapCol,list.MapRow,/overwrite)

end;pro PixelJuggle
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TomoConvergence,diff,threshold,reldiff=reldiff
n=n_elements(diff)

if n lt 2 then return,0b
reldiff=abs(diff[n-2]-diff[n-1])/diff[n-2]*100

return,reldiff lt threshold

end;function TomoConvergence
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ApplyTomogramConstraints,g,constraints,ind=ind
if constraints.nfix ne 0 then g[constraints.indfix]=g[constraints.gfix]
if keyword_set(ind) then begin
    if constraints.blow then g[ind]>=constraints.glow
    if constraints.bhigh then g[ind]<=constraints.ghigh
endif else begin
    if constraints.blow then g>=constraints.glow
    if constraints.bhigh then g<=constraints.ghigh
endelse
end;pro ApplyTomogramConstraints
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FilterSinogram,sinogram,type,kernelSize,Nrho,nviews

; Sinogram vertical
s=size(sinogram,/dim)
btranspose=s[1] ne nviews
if btranspose then sinogram=transpose(sinogram)

; Filter
ker=ReconFilter(type,kernelSize,Nrho,nviews,fourier=fourier,npad=npad)

; Zero-padding
if npad ne 1 then sinogram=extrac(sinogram,0,0,Nrho*npad,nviews)
    
; Filter sinogram
if fourier then begin
    for i=0l, nviews-1 do sinogram[*,i]=fconvol(sinogram[*,i],ker)
endif else begin
    sinogram=convol(sinogram,ker,/EDGE_ZERO)
endelse

; Remove padding
if npad ne 1 then sinogram=sinogram[0:Nrho-1,0:nviews-1]

; Reset sinogram
if btranspose then sinogram=transpose(sinogram)
end;pro FilterSinogram
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SHEPPLOGANPHANTOM,N,M

b=-1
e=1
inc=(e-b)/(N-1.)
x=rebin(b+inc*lindgen(N),N,M,/sample)
inc=(e-b)/(M-1.)
y=rebin(b+inc*lindgen(1,M),N,M,/sample)

f=fltarr(N,M)

;i0=2
;i1=-0.98
;i2=-0.02
;i3=0.01

i0=2
i1=-1
i2=-2./3
i3=1./3

ind=where(((x/0.69)^2+(y/0.92)^2) le 1,ct)
if ct ne 0 then f[ind]=i0

ind=where(((x/0.6624)^2+((y+0.0184)/0.874)^2) le 1,ct)
if ct ne 0 then f[ind]+=i1

c1=(x-0.22)*cos(0.4*!pi)+y*sin(0.4*!pi)
c2=-(x-0.22)*sin(0.4*!pi)+y*cos(0.4*!pi)

ind=where(((c1/0.31)^2+(c2/0.11)^2) le 1,ct)
if ct ne 0 then f[ind]+=i2

c1=(x+0.22)*cos(0.6*!pi)+y*sin(0.6*!pi)
c2=-(x+0.22)*sin(0.6*!pi)+y*cos(0.6*!pi)

ind=where(((c1/0.41)^2+(c2/0.16)^2) le 1,ct)
if ct ne 0 then f[ind]+=i2

ind=where(((x/0.21)^2+((y-0.35)/0.25)^2) le 1,ct)
if ct ne 0 then f[ind]+=i3

ind=where(((x/0.046)^2+((y-0.1)/0.046)^2) le 1,ct)
if ct ne 0 then f[ind]+=i3

ind=where(((x/0.046)^2+((y+0.1)/0.046)^2) le 1,ct)
if ct ne 0 then f[ind]+=i3

ind=where((((x+0.08)/0.046)^2+((y+0.605)/0.023)^2) le 1,ct)
if ct ne 0 then f[ind]+=i3

ind=where(((x/0.023)^2+((y+0.605)/0.023)^2) le 1,ct)
if ct ne 0 then f[ind]+=i3

ind=where((((x-0.06)/0.023)^2+((y+0.605)/0.046)^2) le 1,ct)
if ct ne 0 then f[ind]+=i3

return,f
end;function SHEPPLOGANPHANTOM
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SHEPPLOGANPHANTOMSINOGRAM,N,M,anglesrange,nangles,xc=xc,yc=yc
tomogram=SHEPPLOGANPHANTOM(N,M)
if n_elements(xc) eq 0 then xc=(N-1)/2.
if n_elements(yc) eq 0 then yc=(M-1)/2.

nrho=N>M
rmin=-(nrho-1)/2.

b=0.
e=anglesrange
angles=b+(e-b)/(nangles-1)*lindgen(nangles)
angles*=!pi/180

return,transpose(RADON(tomogram,THETA=angles,xmin=-xc,ymin=-yc,drho=1,NRHO=nrho,RMIN=rmin,/LINEAR))
end;function SHEPPLOGANPHANTOMSINOGRAM
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FBPrecon,sinogram,list,tomogram=tomogram,constraints=constraints

; Some constants
tomoinfo,list,angles=angles,Nx=Nx,Ny=Ny,xmin=xmin,ymin=ymin,Nrho=Nrho,drho=drho,nviews=nviews,projcen=projcen,bsinovert=bsinovert

; Filter sinogram
FilterSinogram,sinogram,list.filtertype,list.filtersize,Nrho,nviews

; Inverse Radon transform
rho = drho*indgen(Nrho)-projcen
tomogram = radon(transpose(sinogram),/BACKPROJECT,rho=rho,theta=angles,/LINEAR,nx=Nx,ny=Ny,xmin=xmin,ymin=ymin)/drho

if keyword_set(constraints) then $
    ApplyTomogramConstraints,tomogram,constraints

end;pro FBPrecon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MLEMrecon,sinogram,list,tomogram=tomogram,convergence=diff,relerrorphantom=relerrorphantom,phantom=phantom,constraints=constraints

; Some constants
if keyword_set(phantom) then begin
    tomogramtheory=SHEPPLOGANPHANTOM(250,250)
    tomogramtheorytotal=total(tomogramtheory)/100.
endif
tomoinfo,list,angles=angles,Nx=Nx,Ny=Ny,xmin=xmin,ymin=ymin,Nrho=Nrho,drho=drho,nviews=nviews,projcen=projcen

; Sinogram in correct orientation: columns are the views
sinogram=transpose(sinogram)
ind=where(~finite(sinogram),ct)
if ct ne 0 then sinogram[ind]=0

; Detector alignement
rho=drho*indgen(Nrho)-projcen

; Initial tomogram
tomogram=replicate(1.,Nx,Ny)

; Tomogram normalisation factor
;tomonorm=replicate(1.,Nrho,nviews)
;FilterSinogram,tomonorm,list.filtertype,list.filtersize,Nrho,nviews
;tomonorm=radon(transpose(tomonorm),rho=rho,theta=angles,nx=Nx,ny=Ny,xmin=xmin,ymin=ymin,/BACKPROJECT,/LINEAR);/drho
tomonorm=radon(replicate(1.,nviews,Nrho),rho=rho,theta=angles,nx=Nx,ny=Ny,xmin=xmin,ymin=ymin,/BACKPROJECT,/LINEAR);/drho
;ind1=where(tomonorm eq 0,ct1,comp=ind2,ncomp=ct2)
;if ct1 ne 0 and ct2 ne 0 then tomonorm[ind1]=min(tomonorm[ind2],/nan)

bCheck=~list.bTomoIter
if bCheck then niter=list.TomoIterMax else niter=list.TomoIter
bShow=list.bshowconvergence

T0=systime(1)
for k=0l,niter-1 do begin
    if bCheck or bShow then tomogramprevious=tomogram
    
    ; Forward projection of tomogram >> Estimated sinogram
    estsinogram=radon(tomogram,theta=angles,xmin=xmin,ymin=ymin,RMIN=-projcen,drho=drho,NRHO=Nrho,/linear)
    
    ; Handle zero's in estimated sinogram
    if k lt list.TomoIter/2 then m=0 $
    else m=list.estcutoff*mean(estsinogram,/nan)
    ind1=where(estsinogram le m,ct1,comp=ind2,ncomp=ct2)
    if ct1 ne 0 and ct2 ne 0 then estsinogram[ind1]=max(estsinogram[ind2])

    ; Backprojection of normalised experimental sinogram >> new tomogram >> normalised tomogram
;    tmp=sinogram/estsinogram
;    FilterSinogram,tmp,list.filtertype,list.filtersize,Nrho,nviews
;    tomogram*=radon(tmp,rho=rho,theta=angles,nx=Nx,ny=Ny,xmin=xmin,ymin=ymin,/BACKPROJECT,/LINEAR)/tomonorm
    tomogram*=radon(sinogram/estsinogram,rho=rho,theta=angles,nx=Nx,ny=Ny,xmin=xmin,ymin=ymin,/BACKPROJECT,/LINEAR)/tomonorm;/drho
    
    if keyword_set(constraints) then $
        ApplyTomogramConstraints,tomogram,constraints
    
    ; Check convergence
    if bCheck or bShow then begin
        tmp=total(abs(tomogramprevious-tomogram))
        if keyword_set(phantom) then tmp2=total(abs(tomogramtheory-tomogram))/tomogramtheorytotal
        if k eq 0 then begin
            diff=tmp
            if keyword_set(phantom) then relerrorphantom=tmp
        endif else begin
            diff=[diff,tmp]
            if keyword_set(phantom) then relerrorphantom=[relerrorphantom,tmp2]
            if bCheck then if TomoConvergence(diff,list.TomoConvergence,reldiff=reldiff) then break
        endelse
    endif
endfor

printw,list.OutID,'MLEM(Radon) #iterations: '+stringr(k)
printw,list.OutID,'MLEM(Radon) time(sec)/iteration: '+stringr((systime(1)-T0)/k)
if n_elements(reldiff) ne 0 then printw,list.OutID,'Relative difference (%): '+stringr(reldiff)

end;pro MLEMrecon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SIRTrecon,sinogram,list,tomogram=tomogram,convergence=diff,relerrorphantom=relerrorphantom,phantom=phantom,constraints=constraints

; Some constants
if keyword_set(phantom) then begin
    tomogramtheory=SHEPPLOGANPHANTOM(250,250)
    tomogramtheorytotal=total(tomogramtheory)/100.
endif
tomoinfo,list,angles=angles,Nx=Nx,Ny=Ny,xmin=xmin,ymin=ymin,Nrho=Nrho,drho=drho,nviews=nviews,projcen=projcen

; Sinogram in correct orientation: columns are the views
sinogram=transpose(sinogram)
ind=where(~finite(sinogram),ct)
if ct ne 0 then sinogram[ind]=0

; Detector alignement
rho=drho*indgen(Nrho)-projcen

; Initial tomogram
tomogram=fltarr(Nx,Ny)

; Tomogram normalisation factor
tomonorm=radon(make_array(Nx,Ny,value=1.),theta=angles,xmin=xmin,ymin=ymin,RMIN=-projcen,drho=drho,NRHO=Nrho,/linear)
;FilterSinogram,tomonorm,list.filtertype,list.filtersize,Nrho,nviews
tomonorm=radon(tomonorm,rho=rho,theta=angles,nx=Nx,ny=Ny,xmin=xmin,ymin=ymin,/BACKPROJECT,/LINEAR);/drho
;ind1=where(tomonorm eq 0,ct1,comp=ind2,ncomp=ct2)
;if ct1 ne 0 and ct2 ne 0 then tomonorm[ind1]=min(tomonorm[ind2],/nan)

bCheck=~list.bTomoIter
if bCheck then niter=list.TomoIterMax else niter=list.TomoIter
bShow=list.bshowconvergence

T0=systime(1)
for k=0l,niter-1 do begin
    if bCheck or bShow then tomogramprevious=tomogram
    
    ; Forward projection of tomogram >> Estimated sinogram
    estsinogram=radon(tomogram,theta=angles,xmin=xmin,ymin=ymin,RMIN=-projcen,drho=drho,NRHO=Nrho,/linear)

    ; Backprojection of normalised experimental sinogram >> new tomogram >> normalised tomogram
;    tmp=sinogram-estsinogram
;    FilterSinogram,tmp,list.filtertype,list.filtersize,Nrho,nviews
;    tomogram+=radon(tmp,rho=rho,theta=angles,nx=Nx,ny=Ny,xmin=xmin,ymin=ymin,/BACKPROJECT,/LINEAR)/tomonorm
    tomogram+=radon(sinogram-estsinogram,rho=rho,theta=angles,nx=Nx,ny=Ny,xmin=xmin,ymin=ymin,/BACKPROJECT,/LINEAR)/tomonorm;/drho
    
    if keyword_set(constraints) then $
        ApplyTomogramConstraints,tomogram,constraints
    
    ; Check convergence
    if bCheck or bShow then begin
        tmp=total(abs(tomogramprevious-tomogram))
        if keyword_set(phantom) then tmp2=total(abs(tomogramtheory-tomogram))/tomogramtheorytotal
        if k eq 0 then begin
            diff=tmp
            if keyword_set(phantom) then relerrorphantom=tmp
        endif else begin
            diff=[diff,tmp]
            if keyword_set(phantom) then relerrorphantom=[relerrorphantom,tmp2]
            if bCheck then if TomoConvergence(diff,list.TomoConvergence,reldiff=reldiff) then break
        endelse
    endif
endfor

printw,list.OutID,'SIRT(Radon) #iterations: '+stringr(k)
printw,list.OutID,'SIRT(Radon) time(sec)/iteration: '+stringr((systime(1)-T0)/k)
if n_elements(reldiff) ne 0 then printw,list.OutID,'Relative difference (%): '+stringr(reldiff)

end;pro SIRTrecon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ARTRays,angles,Ny,xc,yc,projcen,drho

; Grid: Nx x Ny pixels with origin at [xc,yc]
;  Rotate line going through [r0,0] and [r0,Ny-1]
;  over all angles and get linear parameters
; 
; Line through [r0,0] and [r0,Ny-1] is the projection
; on the X-axis (ray sum along the Y-axis)
;
; r0: x-coordinate of the first ray

x=[-projcen,-projcen]
y=[0,Ny-1]-yc
ca=cos(angles)
sa=-sin(angles)
binc=-drho/sa
x1=x[0]*ca+y[0]*sa+xc
x2=x[1]*ca+y[1]*sa+xc
y1=-x[0]*sa+y[0]*ca+yc
y2=-x[1]*sa+y[1]*ca+yc

m=(y2-y1)/(x2-x1)
b=y1-m*x1
ind=where(round(x1) eq round(x2) or abs(x1-x2) lt 1e-2,ct) ; vertical ray stays within a pixel
if ct ne 0 then begin
    binc[ind]=drho*(2*(m[ind] ge 0)-1)
    m[ind]=!values.F_INFINITY
    b[ind]=x1[ind]
endif

ind=where(sa eq 0,ct) ; m=0
if ct ne 0 then binc[ind]=drho

return,{m:m,b:b,binc:binc}
end;function ARTRays
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ARTweightsArea,m,b,binc,Nx,Ny

; Ray with left and right boundaries
if binc lt 0 then begin
    b1=b+binc/2.
    b2=b-binc/2.
endif else begin
    b1=b-binc/2.
    b2=b+binc/2.
endelse

; Handle special cases
if m eq 0 then begin
    y1=round(b1)
    if b1 eq -0.5 then y1=0
    y2=round(b2)
    my=y2-0.5
    
    by1=y1 ge 0 and y1 lt Ny
    by2=y2 ge 0 and y2 lt Ny
    by=by1+2*by2
    case by of
    0: return,{n:0}
    1: return,{n:Nx,ind:lindgen(Nx)+Nx*y1,d:replicate(my-b1,Nx)}
    2: return,{n:Nx,ind:lindgen(Nx)+Nx*y2,d:replicate(b2-my,Nx)}
    3: begin
        if y1 eq y2 then return,{n:Nx,ind:lindgen(Nx)+Nx*y1,d:replicate(1,Nx)} $
        else return,{n:2*Nx,ind:[lindgen(Nx)+Nx*y1,lindgen(Nx)+Nx*y2],$
                d:[replicate(my-b1,Nx),replicate(b2-my,Nx)]}
       endcase
    endcase
endif
if ~finite(m) then begin
    ; (b has to be set to x in ARTRays!!!!)
    x1=round(b1)
    if b1 eq -0.5 then x1=0
    x2=round(b2)
    mx=x2-0.5
    
    bx1=x1 ge 0 and x1 lt Nx
    bx2=x2 ge 0 and x2 lt Nx
    bx=bx1+2*bx2
    case bx of
    0: return,{n:0}
    1: return,{n:Ny,ind:x1+Nx*lindgen(Ny),d:replicate(mx-b1,Ny)}
    2: return,{n:Ny,ind:x2+Nx*lindgen(Ny),d:replicate(b2-mx,Ny)}
    3: begin
        if x1 eq x2 then return,{n:Ny,ind:x1+Nx*lindgen(Ny),d:replicate(1,Ny)} $
        else return,{n:2*Ny,ind:[x1+Nx*lindgen(Ny),x2+Nx*lindgen(Ny)],$
                d:[replicate(mx-b1,Ny),replicate(b2-mx,Ny)]}
       endcase
    endcase
endif

; All intersection points between grid and ray (2 boundaries): 4*(N+1) nodes
; + Add all grid nodes between the two boundaries

; intersections
tx=lindgen(Nx+1)-0.5
ty=lindgen(Ny+1)-0.5
x=[tx,(ty-b1)/m,tx,(ty-b2)/m]
y1=m*tx+b1
y2=m*tx+b2
y=[y1,ty,y2,ty]

; add grid nodes
y1=round(y1)
y2=round(y2)
nadd=y2-y1
if total(nadd) ne 0 then begin
    xadd=chunkindex(nadd)
    yadd=chunklindgen(nadd)
    x=[x,xadd-0.5]
    y=[y,y1[xadd]+yadd+0.5]
endif

; Triangles and area of triangles
TRIANGULATE, x, y, Triangles

x=x[Triangles]
y=y[Triangles]
;xkeep=x
;ykeep=y

small=-0.6
bigx=Nx-0.4
bigy=Ny-0.4
t=x ge small and x le bigx and y ge small and y le bigy
t=t[0,*] and t[1,*] and t[2,*]

; area
d=abs(x[0,*]*(y[1,*]-y[2,*])+x[1,*]*(y[2,*]-y[0,*])+x[2,*]*(y[0,*]-y[1,*]))/2.

; pixel coord.
x=round(total(x,1)/3.)
y=round(total(y,1)/3.)

t=where(x ge 0 and x lt Nx and y ge 0 and y lt Ny and d gt 0 and t,Mnodes)
if Mnodes eq 0 then return,{n:0}

;for i=0l,Mnodes-1 do $
;    plots,xkeep[[0,1,2,0],t[i]],ykeep[[0,1,2,0],t[i]],psym=-2,/data

x=x[t]+Nx*y[t]
d=d[t]
Drizzle,x,d
return,{n:n_elements(x),ind:x,d:d}
end;function ARTweightsArea
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ARTweightsLine,m,b,Nx,Ny

; Handle special cases
if m eq 0 then begin
    y=round(b)
    if y lt 0 or y ge Ny then return,{n:0}
    return,{n:Nx,ind:lindgen(Nx)+Nx*y,d:replicate(1,Nx)}
endif
if ~finite(m) then begin
    x=round(b) ; (b has to be set to x in ARTRays!!!!)
    if x lt 0 or x ge Nx then return,{n:0}
    return,{n:Ny,ind:x+Nx*lindgen(Ny),d:replicate(1,Ny)}
endif

; All intersection points between grid and ray(line): (Nx+Ny+2) nodes
tx=lindgen(Nx+1)-0.5
ty=lindgen(Ny+1)-0.5
x=[tx,(ty-b)/m]
y=[m*tx+b,ty]
t=sort(x)
x=x[t]
y=y[t]

; Distances between nodes: (Nx+Ny+1) dnodes
Mnodes=Nx+Ny+1
dx=x[1:Mnodes]-x[0:Mnodes-1]
dy=y[1:Mnodes]-y[0:Mnodes-1]
d=sqrt(dy*dy+dx*dx)

; Node pixels
x=round(x[0:Mnodes-1]+dx/2.)
y=round(y[0:Mnodes-1]+dy/2.)

; Weights
t=where(x ge 0 and x lt Nx and y ge 0 and y lt Ny and d gt 0,Mnodes)
if Mnodes eq 0 then return,{n:0}
return,{n:Mnodes,ind:x[t]+Nx*y[t],d:d[t]}

; For duplicate ind:
x=x[t]+Nx*y[t]
d=d[t]
Drizzle,x,d
return,{n:n_elements(x),ind:x,d:d}
end;function ARTweightsLine
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ARTrecon,sinogram,list,tomogram=tomogram,convergence=diff,relerrorphantom=relerrorphantom,phantom=phantom,$
    SIRT=SIRT,SART=SART,MLEM=MLEM,OSEM=OSEM,constraints=constraints

; Some constants
if keyword_set(phantom) then begin
    tomogramtheory=SHEPPLOGANPHANTOM(250,250)
    tomogramtheorytotal=total(tomogramtheory)/100.
endif
tomoinfo,list,angles=angles,Nx=Nx,Ny=Ny,xmin=xmin,ymin=ymin,Nrho=Nrho,drho=drho,nviews=NAngles,projcen=projcen
;rho=drho*indgen(Nrho)-projcen

; center of rotation in tomogram pixel coordinates
xc=-xmin
yc=-ymin

; unknowns of the linear system of equations
tomogram=fltarr(Nx*Ny)
if keyword_set(MLEM) or keyword_set(OSEM) then tomogram++
bChunk=keyword_set(SIRT) or keyword_set(SART) or keyword_set(MLEM) or keyword_set(OSEM)
if bChunk then begin
    RayInd=lindgen(Nrho,NAngles) ; [inner loop, outer loop]
    deltakeepNUM=tomogram
    deltakeepDENOM=tomogram
    if keyword_set(OSEM) then begin
        ; Number of projections in each subset
        nproj=(NAngles/list.OSEMsubsets)>1
        
        ; Number of subsets
        nsubsets=NAngles/nproj
        
        ; Projections that are not used
        nleft=NAngles-nsubsets*nproj
        bleft=nleft ne 0

        ; Group subsets
        ind=lindgen(nproj)*nsubsets
        n=nproj*Nrho
        keep=make_array(n+bleft*Nrho,nsubsets,value=-1L)
        for i=0,nsubsets-1 do $
            keep[0:n-1,i]=reform(RayInd[*,ind+i],n)

        ; Add unused projections
        if bleft then $
            keep[n:*,0:nleft-1]=RayInd[*,NAngles-nleft:NAngles-1]
        
        RayInd=temporary(keep)
    endif
endif else RayInd=lindgen(Nrho,NAngles) ; ART: a randomized 1D array doesn't converge as fast

; Parallel set of rays for each angles:
; m: rico
; b: cutoff
; binc: cutoff increment of the parallel set
RayStruc=ARTRays(angles,Ny,xc,yc,projcen,drho)

; Weights for each ray
pW=list.ARTweights
brecal=list.recalcweights
if ~brecal then begin
    s=dimsize(*pW,2)
    brecal=s[0] ne Nrho or s[1] ne NAngles
    if brecal then heap_free,*pW
endif

if brecal then begin
    ; For debugging
;    Window, /Free, xsize=TomoSize(list,/X), ysize=TomoSize(list,/Y), /Pixmap
;    pixindex=!d.window
;    Device, Copy = [0,0, TomoSize(list,/X),TomoSize(list,/Y), 0,0,list.TomoReconIndex]
;    wset,list.TomoReconIndex

    *pW=ptrarr(Nrho,NAngles,/allocate)
    for j=0l,NAngles-1  do begin
        ; Ray characteristics
        m=RayStruc.m[j]
        b=RayStruc.b[j]
        binc=RayStruc.binc[j]
                
        for i=0l,Nrho-1 do begin
            ; Calculate intersections of ray [i,j] with a tomogram
            if list.arthigh then w=ARTweightsArea(m,b,binc,Nx,Ny) $
            else w=ARTweightsLine(m,b,Nx,Ny)

            *(*pW)[i,j]=w
            b+=binc
        endfor
        
        ; For debugging
;        Device, Copy = [0,0, TomoSize(list,/X),TomoSize(list,/Y), 0,0,pixindex]
;        xx=RTStat((*(*pW)[0,j]).ind mod Nx,list.sfh)
;        yy=RTStat(Ny-1-(*(*pW)[0,j]).ind/Nx,list.sfv)
;        plots,xx,yy,/device
;        xx=RTStat((*(*pW)[-1,j]).ind mod Nx,list.sfh)
;        yy=RTStat(Ny-1-(*(*pW)[-1,j]).ind/Nx,list.sfv)
;        plots,xx,yy,/device
;        wait,0.02
    endfor
    
    ; For debugging
;    wdelete,pixindex
endif

; Iterative adaptation of the tomogram
bCheck=~list.bTomoIter
if bCheck then niter=list.TomoIterMax else niter=list.TomoIter
bShow=list.bshowconvergence

T0=systime(1)
sRayInd=dimsize(RayInd,2)
for k=0l,niter-1 do begin
    if bCheck or bShow then tomogramprevious=tomogram
    
    ; Randomize chunks
    if sRayInd[1] gt 1 then $
        RayInd=RayInd[*,sort(randomu(seed,sRayInd[1]))]
    
    ; Outer loop
    for jRay=0l,sRayInd[1]-1  do begin
        ; Inner loop
        for iRay=0l,sRayInd[0]-1 do begin
            ; Index of ray in sinogram
            if RayInd[iRay,jRay] eq -1 then continue
            i=RayInd[iRay,jRay] mod Nrho ; rho
            j=RayInd[iRay,jRay]/Nrho ; angle
            
            ; Intersection weights of this ray with the tomogram (only the non-zero weights)
            w=*(*pW)[i,j]
            
            if w.n ne 0 then begin
                p=sinogram[i,j]
                q=total(w.d*tomogram[w.ind])
                if bChunk then begin
                    if keyword_set(MLEM) or keyword_set(OSEM) then begin
                        if k lt list.TomoIter/2 then thres=0 $
                        else thres=list.estcutoff*mean(q,/nan)
                        ind1=where(q le thres,ct1,comp=ind2,ncomp=ct2)
                        if ct1 ne 0 and ct2 ne 0 then q[ind1]=max(q[ind2])
    
                        deltakeepNUM[w.ind]+=p/q*w.d
                        deltakeepDENOM[w.ind]+=w.d
                    endif else begin
                        if list.altavg then begin
                            deltakeepNUM[w.ind]+=((p-q)/total(w.d))*w.d
                            deltakeepDENOM[w.ind]+=w.d
                        endif else begin
                            deltakeepNUM[w.ind]+=((p-q)/total(w.d*w.d))*w.d
                            deltakeepDENOM[w.ind]+=w.d
                        endelse
                    endelse
                endif else begin
                    tomogram[w.ind]+=((p-q)/total(w.d*w.d))*w.d
                    if keyword_set(constraints) then $
                        ApplyTomogramConstraints,tomogram,constraints,ind=w.ind
                endelse
            endif
        endfor
        
        ; Adapt tomogram after inner loop
        if keyword_set(SART) then begin
            ind=where(deltakeepDENOM eq 0,ct)
            if ct ne 0 then deltakeepDENOM[ind]=1
            tomogram+=deltakeepNUM/deltakeepDENOM
            deltakeepNUM=fltarr(Nx*Ny)
            deltakeepDENOM=deltakeepNUM

            if keyword_set(constraints) then $
                ApplyTomogramConstraints,tomogram,constraints
        endif
        
        if keyword_set(OSEM) then begin
            ind=where(deltakeepDENOM eq 0,ct)
            if ct ne 0 then deltakeepDENOM[ind]=1
            tomogram*=deltakeepNUM/deltakeepDENOM
            deltakeepNUM=fltarr(Nx*Ny)
            deltakeepDENOM=deltakeepNUM

            if keyword_set(constraints) then $
                ApplyTomogramConstraints,tomogram,constraints
        endif
    endfor

    ; Adapt tomogram after outer loop
    if keyword_set(SIRT) then begin
        ind=where(deltakeepDENOM eq 0,ct)
        if ct ne 0 then deltakeepDENOM[ind]=1
        tomogram+=deltakeepNUM/deltakeepDENOM
        deltakeepNUM=fltarr(Nx*Ny)
        deltakeepDENOM=deltakeepNUM

        if keyword_set(constraints) then $
            ApplyTomogramConstraints,tomogram,constraints
    endif
    
    if keyword_set(MLEM) then begin
        ind=where(deltakeepDENOM eq 0,ct)
        if ct ne 0 then deltakeepDENOM[ind]=1
        tomogram*=deltakeepNUM/deltakeepDENOM
        deltakeepNUM=fltarr(Nx*Ny)
        deltakeepDENOM=deltakeepNUM

        if keyword_set(constraints) then $
            ApplyTomogramConstraints,tomogram,constraints
    endif
    
    ; Check convergence
    if bCheck or bShow then begin
        tmp=total(abs(tomogramprevious-tomogram))
        if keyword_set(phantom) then tmp2=total(abs(tomogramtheory-tomogram))/tomogramtheorytotal
        if k eq 0 then begin
            diff=tmp
            if keyword_set(phantom) then relerrorphantom=tmp
        endif else begin
            diff=[diff,tmp]
            if keyword_set(phantom) then relerrorphantom=[relerrorphantom,tmp2]
            if bCheck then if TomoConvergence(diff,list.TomoConvergence,reldiff=reldiff) then break
        endelse
    endif
endfor

; Some statistics
name='ART'
if keyword_set(SIRT) then name='SIRT'
if keyword_set(SART) then name='SART'
if keyword_set(MLEM) then name='MLEM'
if keyword_set(OSEM) then name='OSEM'

printw,list.OutID,name+' #iterations: '+stringr(k)
printw,list.OutID,name+' time(sec)/iteration: '+stringr((systime(1)-T0)/k)
if n_elements(reldiff) ne 0 then printw,list.OutID,'Relative difference (%): '+stringr(reldiff)

;window,0
;for j=0l,NAngles-1 do begin
;    ; Plot ray tomogram pixel weights
;    weights=fltarr(Nx,Ny)
;        
;    w=*pW[0,j]
;    if w.n ne 0 then weights[w.ind]=w.d
;    w=*pW[Nrho-1,j]
;    if w.n ne 0 then weights[w.ind]=w.d
;        
;    for i=0l,Nrho-1 do begin
;        w=*pW[i,j]
;        ww=weights
;        if w.n ne 0 then ww[w.ind]=w.d
;        tvscl,ww
;    endfor
;endfor

tomogram=reform(tomogram,Nx,Ny,/overwrite)
end;pro ARTrecon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TOMOrecon,sinogram,list,tomogram=tomogram,average=average,weightsptr=weightsptr

T0=systime(1)

if ~list.updatesinogram then keep=sinogram

; Reconstruction info
tomoinfo,list,angles=angles,Nx=Nx,Ny=Ny,xmin=xmin,ymin=ymin,Nrho=Nrho,drho=drho,nviews=nviews,projcen=projcen,bsinovert=bsinovert

; Weighted reconstruction?
if keyword_set(average) then begin
    ; Weights
    average=1b
    if keyword_set(weightsptr) then begin
        s=dimsize(*weightsptr,2)
        average=Nx eq s[0] and Ny eq s[1]
        if average then weights=*weightsptr
    endif else average=0b
endif else average=0b

; sinogram orientation: 
;    horizontal: translation
;    vertical: rotation
if ~bsinovert then sinogram=transpose(sinogram)
phantom=0
if phantom then sinogram=SHEPPLOGANPHANTOMSINOGRAM(250,250,180,181)

; Apply weights to sinogram
if average then begin
    ; Possitive weights
    weights>=0
    
    if list.tomoequalweights then begin
        ; Binarize
        mi=min(weights,max=ma)
        weights gt= mi+(ma-mi)*list.tomoavgcutoff/100.
    endif else begin
        ; Remove small weights
        mi=min(weights,max=ma)
        ind=where(weights lt mi+(ma-mi)*list.tomoavgcutoff/100.,ct)
        if ct ne 0 then weights[ind]=0
    endelse
    ;common cbin,weights
    
    ; Apply weights
    sinogram*=transpose(radon(float(weights),theta=angles,xmin=xmin,ymin=ymin,RMIN=-projcen,drho=drho,NRHO=Nrho,/linear))
    if ~bsinovert then sinogram=transpose(sinogram)
endif

; Constraints
glow=0
blow=list.tomopositive
nfix=0l
indfix=-1l
gfix=0

if list.tomoedge then begin
    tmp=bytarr(Nx,Ny)
    mi=min(tmp,/nan)
    offb=0>((list.tomoedgecutoff-1)<(Nx-1))
    offe=0>(Nx-list.tomoedgecutoff)
    tmp[0:offb,*]=1b
    tmp[offe:*,*]=1b
    offb=0>((list.tomoedgecutoff-1)<(Ny-1))
    offe=0>(Ny-list.tomoedgecutoff)
    tmp[*,0:offb]=1b
    tmp[*,offe:*]=1b
    ind=where(tmp,ct)
    nfix+=ct
    indfix=[indfix,ind]
endif

if list.tomocircle then begin
    ; Center of rotation
    x0=-xmin
    y0=-ymin
    
    ; Radius of clipping circle
    R2cut=list.tomocircleR
    if R2cut le 0 then R2cut=(Nx<Ny)/2.
    R2cut^=2.
    
    ; Clipping
    ind=where((rebin(lindgen(Nx)-x0,Nx,Ny,/sample)^2.+$
      rebin(lindgen(1,Ny)-y0,Nx,Ny,/sample)^2.) gt R2cut,ct)

    if ct ne 0 then begin
        nfix+=ct
        indfix=[indfix,ind]
    endif
endif

if nfix ne 0 then indfix=indfix[1:*]
if ~list.tomopostcorrect then constraints={nfix:nfix,indfix:indfix,gfix:gfix,blow:blow,glow:glow,bhigh:0b}

; Reconstruct
case list.ReconType of
0:    FBPrecon,sinogram,list,tomogram=tomogram,constraints=constraints
1:    MLEMrecon,sinogram,list,tomogram=tomogram,convergence=convergence,relerrorphantom=relerrorphantom,phantom=phantom,constraints=constraints
2:    SIRTrecon,sinogram,list,tomogram=tomogram,convergence=convergence,relerrorphantom=relerrorphantom,phantom=phantom,constraints=constraints
3:    ARTrecon,sinogram,list,tomogram=tomogram,convergence=convergence,relerrorphantom=relerrorphantom,phantom=phantom,constraints=constraints,/MLEM
4:    ARTrecon,sinogram,list,tomogram=tomogram,convergence=convergence,relerrorphantom=relerrorphantom,phantom=phantom,constraints=constraints,/OSEM
5:    ARTrecon,sinogram,list,tomogram=tomogram,convergence=convergence,relerrorphantom=relerrorphantom,phantom=phantom,constraints=constraints
6:    ARTrecon,sinogram,list,tomogram=tomogram,convergence=convergence,relerrorphantom=relerrorphantom,phantom=phantom,constraints=constraints,/SART
7:    ARTrecon,sinogram,list,tomogram=tomogram,convergence=convergence,relerrorphantom=relerrorphantom,phantom=phantom,constraints=constraints,/SIRT
endcase

; drho=1 because otherwise the ART algorithm introduces artifacts (because of ???): instead, divide by drho here
tomogram/=list.tomodrho

; Apply weights to tomogram
if average then begin
    ; Correct tomogram
    indzero=where(weights eq 0,nzero,comp=indnotzero,ncomp=nindnotzero)
    if nzero ne 0 then tomogram[indzero]=0
    if nindnotzero ne 0 then tomogram[indnotzero]/=weights[indnotzero]
endif

; Post reconstruction correction
if list.tomopostcorrect then begin
    if nfix ne 0 then tomogram[indfix]=gfix
    if blow then tomogram>=glow
endif

; Marker tomogram for debugging
;    tomogram[0,0]=max(tomogram)
;    tomogram[0,-1]=max(tomogram)
;    tomogram[-1,0]=max(tomogram)
;    tomogram[-1,-1]=max(tomogram)

; New sinogram
if list.updatesinogram then sinogram=transpose(radon(tomogram,theta=angles,xmin=xmin,ymin=ymin,RMIN=-projcen,drho=drho,NRHO=Nrho,/LINEAR)) $
else sinogram=temporary(keep)

printw,list.OutID,'Total reconstruction time (sec): '+stringr(systime(1)-T0)

; Debug
if n_elements(convergence) gt 1 and list.bshowconvergence then begin
    convergence[0]=!values.F_NAN
    
    old=!d.WINDOW
    window,0
    ntmp=n_elements(convergence)
    plot,lindgen(ntmp)+1,convergence,psym=-1,/ylog,xtitle='Iteration',ytitle='Adjustment difference',title='Convergence'
    
    if phantom then begin
        window,1
        plot,lindgen(ntmp)+1,relerrorphantom,psym=-1,/ylog,xtitle='Iteration',ytitle='Relative difference with theory',title='Convergence'
        
        window,2
        real=SHEPPLOGANPHANTOM(250,250)
        plot,real[*,150],yrange=[min(real[*,150])<min(tomogram[*,150]),max(real[*,150])>max(tomogram[*,150])]
        oplot,tomogram[*,150],color=200

;        case list.ReconType of
;        1:    name='MLEMRadon'
;        2:    name='SIRTRadon'
;        3:    name='MLEM'
;        4:    name='OSEM'
;        5:    name='ART'
;        6:    name='SART'
;        7:    name='SIRT'
;        endcase
;        file='C:\Theorie\Thesis\Manuscript\Main\Figures\Tomophantom\iterconv\convergence\'+name+'.sav'
;        save,convergence,relerrorphantom,filename=file
    endif
    
    wset,old
endif

end;pro TOMOrecon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Makewplot3DImage,list,index,tomogram=tomogram,PATH_XY=PATH_XY,btomoedit=btomoedit

bTomo=widget_info(list.TomoReconID,/VALID_ID) eq 1 and not keyword_set(btomoedit)

; Prepare image/sinogram
image=reform((*list.data)[index[0],*,*],list.MapCol,list.MapRow)

;return,image 

if bTomo then $
    TOMOrecon,image,list,tomogram=tomogram,average=list.tomoavg,weight=list.tomoweights

; Try all
;if bTomo then begin
;    for i=1,7 do begin
;        image=reform((*list.data)[index[0],*,*],list.MapCol,list.MapRow)
;        list.ReconType=i
;        TOMOrecon,image,list,tomogram=tomogram,average=list.tomoavg,weight=list.tomoweights
;    endfor
;endif

return,image
end;function Makewplot3DImage
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TomogramWeightFix,ev
widget_control,ev.top,get_uvalue=list

sinogram=reform((*list.data)[list.NSelect,*,*],list.MapCol,list.MapRow)

TOMOrecon,sinogram,list,tomogram=tomogram,average=0

*list.tomoweights=tomogram

WorkingDataCorrect,list
WIDGET_CONTROL, ev.top, set_uvalue=list
RefreshDisplayXDI,list
RefreshWindowXDI,ev
end;pro TomogramWeightFix
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GetGroupindices,NMapsi,top,indices=NSelect

prom='0-'+stringr(NMapsi-1)
str=PromptNumber(prom,top,'Select groups:')
str=str_sep(strtrim(strcompress(str),2),',')
NMaps=n_elements(str)
ind=where(strlen(str) eq 0,ct)
if ct ne 0 then return,0

if NMaps ne 0 then begin
    NSelect=-1
    for i=0l,NMaps-1 do begin
        stri=str_sep(str[i],'-')
        case n_elements(stri) of
        1: NSelect=[NSelect,long(str[i])>0]
        2:     begin
            b=(long(stri[0])<long(stri[1]))>0
            e=(long(stri[0])>long(stri[1]))>0
            ni=e-b+1
            if ni gt 0 then NSelect=[NSelect,b+lindgen(ni)]
            endcase
        else:
        endcase
    endfor
    NSelect=NSelect[1:*]
    NMaps=n_elements(NSelect)
endif

return,NMaps
end;function GetGroupindices
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro savecolorbar,list,dummy=dummy

imgbound=([(list.scaleimg-100)>0,list.scaleimg]<100)/100.

bar = INDGEN(1,256)
if list.scaleimg gt 100 then begin
    rico=100./(200.-list.scaleimg)
    intercept=(list.scaleimg-100.)/(list.scaleimg-200.)*255
endif else begin
    rico=100./list.scaleimg
    intercept=0.
endelse

bar = byte(0>(rico*bar+intercept)<255)
bar=[bar,bar]

ma=max((*list.data)[list.nselect,*,*],min=mi)

TVcolorbar,[0.8,0.15,0.9,0.95],range=[mi,ma],bar=bar
end;pro savecolorbar
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro wplot3D_tomosep,list,dummy=dummy

ID=widget_info(list.top,FIND_BY_UNAME='Binned pixels')
widget_control,ID,get_uvalue=sample
LoadctX,list
imgbound=([(list.scaleimg-100)>0,list.scaleimg]<100)/100.

ma=max((*list.ptrtomocurrent),min=mi)
m=mi+(ma-mi)*imgbound
tvscl,congrid(m[0]>(*list.ptrtomocurrent)<m[1],TomoSize(list,/X),TomoSize(list,/Y),interp=~sample),/order

ID=widget_info(list.top,FIND_BY_UNAME='Save with labels')
widget_control,ID,get_uvalue=bLabel
if bLabel then begin
    labeli=0>list.labeli<(n_elements(((*list.pmstruc).(list.NSelect+1)).data)-1)
    labeli=string(((*list.pmstruc).(list.NSelect+1)).data[labeli])
    
    if bLabel then xyouts,list.legendpos[0],list.legendpos[1],labeli,/normal,CHARSIZE=list.charpr*TomoSize(list,/X),charthick=1.5,col=list.labelcol
endif

end;pro wplot3D_tomosep
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro wplot3D,list,tomosep=tomosep,blendupdate=blendupdate

ID=widget_info(list.top,FIND_BY_UNAME='Binned pixels')
widget_control,ID,get_uvalue=sample

brgb=0b
if keyword_set(blendupdate) then begin
    img=reform(alpha_blending(list,error=error))
    if error then return
    brgb=1b
endif else $
if (*list.pmstruc).(list.NSelect+1).combi eq 2 then begin
    img=reform(long((*list.data)[list.NSelect,*,*]),list.MapCol,list.MapRow)
    brgb=1b
endif

if brgb then begin
;    device,GET_DECOMPOSED=cState
;    device,DECOMPOSE = 1
    if sample then begin
        img=congrid(img,list.XSize,list.YSize,interp=~sample)
        img=[[[img]],[[ishft(img,-8)]],[[ishft(img,-16)]]] and 255b
        tv,img,/order,true=3
    endif else begin
        img=[[[img]],[[ishft(img,-8)]],[[ishft(img,-16)]]] and 255b
        tvscl,congrid(img,list.XSize,list.YSize,3,interp=~sample),/order,true=3
    endelse
;    device,DECOMPOSE = cState
    return
endif

LoadctX,list

ID=widget_info(list.top,FIND_BY_UNAME='Suspend update')
widget_control,ID,get_uvalue=btomoedit
image=Makewplot3DImage(list,list.NSelect,tomogram=tomogram,btomoedit=btomoedit)
imgbound=([(list.scaleimg-100)>0,list.scaleimg]<100)/100.
bTomo= n_elements(tomogram) ne 0

; ----Plot Image (/order to surpress flipping)----
ID=widget_info(list.top,FIND_BY_UNAME='Log Scale')
widget_control,ID,get_uvalue=log
if log eq 1 then image=alog10(image>0)
ma=max(image,min=mi)
m=mi+(ma-mi)*imgbound

tvscl,congrid(m[0]>image<m[1],list.XSize,list.YSize,interp=~sample),/order

lthick=1.5

ID=widget_info(list.top,FIND_BY_UNAME='Add contour')
widget_control,ID,get_uvalue=bcontour
if bcontour and m[0] ne m[1] then begin
    position=[RTStat([0,list.MapCol-1],list.sfh),RTStat([0,list.MapRow-1],list.sfv)]
    position=position[[0,2,1,3]]
    CONTOUR,reverse(image,2),xstyle=1,ystyle=1,levels=m,pos=position,/device,/noerase,xthick = lthick,ythick = lthick, thick=lthick,$
            color=(list.scaleimg lt 100)?list.boardercol:(not list.boardercol);,/PATH_DATA_COORDS,PATH_XY=PATH_XY,PATH_INFO=PATH_INFO

    ;for i=0l,n_tags(path_info)-1 do begin
    ;    if path_info[i].n gt 50 then begin
    ;        i0=path_info[i].offset
    ;        i1=path_info[i].offset+path_info[i].n-1
    ;        ;plots,PATH_XY[*,i0:i1]
    ;        for j=i0+1,i1,20 do begin
    ;;            xy=PATH_XY[*,j]
    ;;            rico=(PATH_XY[1,j-1]-PATH_XY[1,j+1])/(PATH_XY[0,j-1]-PATH_XY[0,j+1])
    ;;            rico=-1/rico
    ;;            b=xy[1]-rico*xy[0]
    ;;            off=10/(rico*rico+1)
    ;;            x1=xy[0]+off
    ;;            x2=xy[0]-off
    ;;            y1=rico*x1+b
    ;;            y2=rico*x2+b
    ;;            plots,[x1,x2],[y1,y2]
    ;            d=total((PATH_XY[*,i0:i1]-rebin(PATH_XY[*,j],2,path_info[i].n,/sample))^2.,1)
    ;            plot,max(d)-d
    ;            tmp=min(d,ind)
    ;            ;plots,[[PATH_XY[*,j]],[PATH_XY[*,i0+ind[0]]]],linestyle=0
    ;            ;AddGaussFit,list.data,list.NSelect,x1,x2,y1,y2
    ;        endfor
    ;    endif
    ;endfor
endif

ID=widget_info(list.top,FIND_BY_UNAME='Save with labels')
widget_control,ID,get_uvalue=bLabel
if bLabel then begin
    labeli=0>list.labeli<(n_elements(((*list.pmstruc).(list.NSelect+1)).data)-1)
    ;labeli=((*list.pmstruc).(list.NSelect+1)).data[labeli]
    ;labeli=string(float(labeli),format='(f4.2)')
    labeli=string(((*list.pmstruc).(list.NSelect+1)).data[labeli])
    xyouts,list.legendpos[0],list.legendpos[1],labeli,/normal,CHARSIZE=list.charpr*list.XSize,charthick=1.5,col=list.labelcol
endif

; ----Plot ROI's----
if list.pverts ne -1 then begin
ID=widget_info(list.top,FIND_BY_UNAME='Fill')
widget_control,ID,get_uvalue=fill
ID=widget_info(list.top,FIND_BY_UNAME='Continuous')
widget_control,ID,get_uvalue=cont
ind=0
if list.nrROI ne 0 then $
    for i=0l,list.nrROI-1 do begin
       ind2=ind+(*list.nrPoints)[i]-1
       x=(*list.xverts)[ind:ind2]
       y=(*list.yverts)[ind:ind2]
;       print,x
;       print,y
       if cont eq 1 then begin
           x=[x,(*list.xverts)[ind]]
           y=[y,(*list.yverts)[ind]]
       endif
       DRAW_POL, RTStat(x,list.sfh),RTStat(y,list.sfv),(*list.highlight)[i], FILL = fill
       ind=ind2+1
    endfor
if (list.pverts ge ind) then DRAW_POL, RTStat((*list.xverts)[ind:list.pverts],list.sfh),RTStat((*list.yverts)[ind:list.pverts],list.sfv),0
endif

; ----Plot reconstruct from sinogram----
if bTomo then begin
    *list.ptrtomocurrent=tomogram
    if ~keyword_set(tomosep) then begin
        ma=max(tomogram,min=mi)
        m=mi+(ma-mi)*imgbound
        wset,list.TomoReconIndex
        tvscl,congrid(m[0]>tomogram<m[1],TomoSize(list,/X),TomoSize(list,/Y),interp=~sample),/order
        if bLabel then xyouts,list.legendpos[0],list.legendpos[1],labeli,/normal,CHARSIZE=list.charpr*TomoSize(list,/X),charthick=1.5,col=list.labelcol
    endif
endif

end;pro wplot3D
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wplot2Dgetxrange,list,n
xrange=list.coord[*,0]
if xrange[0] eq xrange[1] then xrange=[0,n-1.] else $
if list.unit[0] ne 'deg' then xrange*=1000
xrange=getxdicoord(list,xrange)
return,xrange
end;function wplot2Dgetxrange
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ConvertXDIunit,unit
if unit eq 'mm' then return,msymbols('micro')+'m' else return,unit
end;function ConvertXDIunit
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro wplot2D,list,PS=PS,CGM=CGM

; ----Plot Profile----
if keyword_set(PS) then LoadctRev,39,/silent else LoadctRev,0,/silent,rev=-1
if keyword_set(CGM) then position=[0.1,0.1,0.95,0.58]

y=reform((*list.data)[list.NSelect,*,*])
n=n_elements(y)

xrange=wplot2Dgetxrange(list,n)
x=MakeRange(xrange[0],xrange[1],(xrange[1]-xrange[0])/(n-1))

labeli=0>list.labeli<(n_elements(((*list.pmstruc).(list.NSelect+1)).data)-1)
title=string(((*list.pmstruc).(list.NSelect+1)).data[labeli])

ID=widget_info(list.top,FIND_BY_UNAME='Log Scale')
widget_control,ID,get_uvalue=log

ind=where(~finite(y),ct)
if ct ne 0 then y[ind]=0
y^=list.scaleimg/100.

CHARSIZE=list.charpr*100

if list.TRAN then plot,y,x,xlog=log,yrange=reverse(xrange),xstyle=1,ystyle=1,ytitle='Spatial distance ('+ConvertXDIunit(list.unit[0])+')',$
         xtitle='I',linestyle=0,title=title,CHARSIZE=CHARSIZE,position=position $
else plot,x,y,ylog=log,xrange=xrange,xstyle=1,ystyle=1,xtitle='Spatial distance ('+ConvertXDIunit(list.unit[0])+')',$
         ytitle='I',linestyle=0,title=title,CHARSIZE=CHARSIZE,position=position

end;pro wplot2D
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro wplot2D1,list,NSelect=NSelect,PS=PS,CGM=CGM

; ----Plot Profile----
if keyword_set(PS) then LoadctRev,39,/silent else LoadctRev,0,/silent,rev=-1
if keyword_set(CGM) then position=[0.1,0.1,0.60,0.95]

TDblock=(*list.data)^(list.scaleimg/100.)
sTDblock=size(TDblock)
nrtif=list.MapCol*list.MapRow
;acol=list.Nmaps
acol=n_elements(NSelect)

xrange=wplot2Dgetxrange(list,nrtif)
x=MakeRange(xrange[0],xrange[1],(xrange[1]-xrange[0])/(nrtif-1.))

; Make a vector of 16 points, A[i] = 2pi/16:
A = FINDGEN(17) * (!PI*2/16.)
; Define the symbol to be a circle with 16 points,
; and set the filled flag:
USERSYM, 0.8*COS(A), 0.8*SIN(A), /FILL

if list.markers eq '' then nmarkers=0 else begin
    markers=float(str_sep(list.markers,','))
    nmarkers=n_elements(markers)
endelse

CHARSIZE=list.charpr*100
SYMSIZE=CHARSIZE/2.

if list.TRAN then begin
   for k=0l,acol-1 do begin
        jsel=NSelect[k]
     may=max(TDblock[jsel,*,*],min=miy)
     if miy eq may then y=replicate(0.5,nrtif) else y=(TDblock[jsel,*,*]-miy)/(may-miy)
     y+=2*k

     if k eq 0 then plot,y,x,yrange=reverse(xrange),xrange=[0,2*acol],$
     xstyle=1,ystyle=1,ytitle='Spatial distance ('+ConvertXDIunit(list.unit[0])+')',CHARSIZE=CHARSIZE,$
     xtitle='Normalised within each group',psym=-8,position=position,SYMSIZE=SYMSIZE $
     else oplot,y,x,psym=-8,SYMSIZE=SYMSIZE

     plots,[2*k,2*k],!Y.crange,NOCLIP = 0
     ;plots,[2*k+2,2*k+2],!Y.crange,NOCLIP = 0

     labeli=0>list.labeli<(n_elements(((*list.pmstruc).(jsel+1)).data)-1)
     xyouts,2*k+1.25,0.9*(!Y.crange[1]-!Y.crange[0])+!Y.crange[0],string(((*list.pmstruc).(jsel+1)).data[labeli]),$
         CHARSIZE=CHARSIZE,orientation=-90

     for m=0,nmarkers-1 do plots,!X.crange,[markers[m],markers[m]],NOCLIP = 0,linestyle=2
   endfor
endif else begin
   for k=0l,acol-1 do begin
     jsel=NSelect[k]
     may=max(TDblock[jsel,*,*],min=miy)
     if miy eq may then y=replicate(0.5,nrtif) else y=(TDblock[jsel,*,*]-miy)/(may-miy)
     y+=2*k

     if k eq 0 then plot,x,y,xrange=xrange,yrange=[0,2*acol],$
     xstyle=1,ystyle=1,xtitle='Spatial distance ('+ConvertXDIunit(list.unit[0])+')',CHARSIZE=CHARSIZE,$
     ytitle='Normalised within each group',psym=-8,position=position,SYMSIZE=SYMSIZE $
     else oplot,x,y,psym=-8,SYMSIZE=SYMSIZE

     plots,!X.crange,[2*k,2*k],NOCLIP = 0
     ;plots,!X.crange,[2*k+2,2*k+2],NOCLIP = 0

     labeli=0>list.labeli<(n_elements(((*list.pmstruc).(jsel+1)).data)-1)
     xyouts,0.9*(!X.crange[1]-!X.crange[0])+!X.crange[0],2*k+1.25,string(((*list.pmstruc).(jsel+1)).data[labeli]),$
         CHARSIZE=CHARSIZE

     for m=0,nmarkers-1 do plots,[markers[m],markers[m]],!Y.crange,NOCLIP = 0,linestyle=2
   endfor
endelse

end;pro wplot2D1
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshDisplayXDI,list,blendupdate=blendupdate

if not ptr_valid(list.data) then return

; Plot
wset,list.drawindex
case list.ScanDim of
1: wplot2D,list
2: wplot3D,list,blendupdate=blendupdate
else:return
endcase
GetSysVar,list.sysvar

; Statistics
image=(*list.data)[list.NSelect,*,*]

ID=widget_info(list.top,FIND_BY_UNAME='Stat_Mean')
widget_control,ID,set_value='0'+stringr(177b)+'0'
ind=where(image ne 0,ct)
if ct gt 1 then begin
    a=mean(image[ind])
    b=stddev(image[ind])
    CATCH, Error_status
    if Error_status eq 0 then widget_control,ID,set_value=rounderror(a,b,format=1)
endif
mi=min(image,max=ma)
ID=widget_info(list.top,FIND_BY_UNAME='Stat_Min')
widget_control,ID,set_value=stringr(mi)
ID=widget_info(list.top,FIND_BY_UNAME='Stat_Max')
widget_control,ID,set_value=stringr(ma)
ID=widget_info(list.top,FIND_BY_UNAME='Stat_FWHM')
widget_control,ID,set_value=''

wset,list.pixindex
Device, Copy = [0,0, list.XSize,list.YSize, 0,0,list.drawindex]

end;pro RefreshDisplayXDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AnimateTomo,ev
widget_control,ev.top,get_uvalue=list
if ~widget_info(list.TomoReconID,/valid_id) then return

temp=CutPath(list.file,file=file)
path2=SelFile(list.path,file+'.tif',SaveImageFilter(/copy),'Save Animation series....')
bsave = path2 ne ''
if bsave then temp=CutPath(path2,path=path,file=file,ext=ext)

; Tomogram coordinate system: X points to the right, Y points down
; Object rotation: start=Angle offset -> end=Angle offset - 180
;     - visible: counter clockwise
;    - flip Y axis: clockwise
;    - Negate rotation: counter clockwise
; So in the tomogram coordinate system, a positive rotation causes the object to rotate counter clockwise

tomoinfo,list,angles=angles,nviews=n,xmin=xmin,ymin=ymin,bsinovert=bsinovert

keep=list.addangle

nmag=3
indy=rebin(indgen(nmag),nmag,nmag,/sample)-(nmag-1)/2.
indx=transpose(indy)

indx+=RTStat(-xmin,list.sfh)
indy+=RTStat(-ymin,list.sfv)

progressBar = Obj_New("PROGRESSBAR",title="Rotating ...")
progressBar -> Start

format='(I0'+stringr(strlen(stringr(n-1)))+')'
if bsinovert then offmax=list.YSize-1 else offmax=list.XSize-1
for i=0l,n-1 do begin
    
    IF progressBar -> CheckCancel() THEN BEGIN
        printw,list.outid,'Process CANCELLED!'
        break
    ENDIF
    
    ; Tomogram
    list.addangle=angles[i]
    widget_control,ev.top,set_uvalue=list
    RefreshDisplayXDI,list

    ; Original sinogram
    keep2=list.TomoReconID
    list.TomoReconID=0
    RefreshDisplayXDI,list
    list.TomoReconID=keep2

    ; Sinogram progress
    wset,list.drawindex
    if bsinovert then begin
        off=RTStat(i+1,list.sfv)
        if off lt offmax then begin
            img=intarr(list.XSize,offmax-off+1)
            tvscl,img,/order
        endif
    endif else begin
        off=RTStat(i+1,list.sfh)
        if off lt offmax then begin
            tmp=offmax-off+1
            img=intarr(tmp,list.YSize)
            tvscl,img,list.XSize-tmp,0,/device
        endif
    endelse
    
    ; Center of rotation
    wset,list.TomoReconIndex
    plots,indx,indy,/device

    ; Save image
    if bsave then begin
        SaveAllXDI,ev,current=path+file+string(i,format=format)+ext
    endif
    
    progressBar -> Update, (i+1.)/n*100
endfor

progressBar -> Destroy

list.addangle=keep
widget_control,ev.top,set_uvalue=list
RefreshDisplayXDI,list

end;pro AnimateTomo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshDisplayXDI_LoadctX,DATA=DATA

keep=!d.window
TVLCT, R, G, B, /GET
(*DATA.ctl.R)=R
(*DATA.ctl.G)=G
(*DATA.ctl.B)=B
RefreshDisplayXDI,DATA
wset,keep
end;pro RefreshDisplayXDI_LoadctX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeleteXDIGroupsHelper,list,val

; Cut groups
list.dataorg=ShrinkArray(list.dataorg,val,count=NMAP)
list.data=ShrinkArray(list.data,val,count=NMAP)

; Rearrange group descriptions
pmstruc=create_struct('t',-1)
pmtags=-1
for i=0l,list.NMaps-1 do begin
    ind=where(val eq i,ct)
    if ct eq 0 then begin
        j=i+1
        pmtags=[pmtags,(*list.pmtags)[j]]
        pmstruc=create_struct(pmstruc,'t'+stringr((*list.pmtags)[j]),(*list.pmstruc).(j))
    endif
endfor

(*list.pmstruc)=pmstruc
(*list.pmtags)=pmtags

list.NMaps=NMAP

end;pro DeleteXDIGroupsHelper
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeleteXDIGroups,ev

; Get group indices
widget_control,ev.id,get_value=val
val=strcompress(val,/remove_all)
val=STRSPLIT(val,',',/EXTRACT)
tmp=where(val ne '',ct)
if ct le 0 then return
val=fix(val[tmp])
widget_control,ev.top,get_uvalue=list
if not ptr_valid(list.data) then return
tmp=where((val lt 0) or (val gt list.NMaps-1),ct)
if ct ne 0 then return

; Delete
DeleteXDIGroupsHelper,list,val

; Finish
list.NSelect<=(list.NMaps-1)
ID=widget_info(ev.top,FIND_BY_UNAME='NRslider')
widget_control,ID,set_slider_max=list.NMaps-1,set_value=list.NSelect

widget_control,ev.top,set_uvalue=list
RefreshDisplayXDI,list
RefreshWindowXDI,ev
end; pro DeleteXDIGroups
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PlotCorrelateXDIGroups,list,PS=PS,CGM=CGM

if keyword_set(PS) then LoadctRev,39 else LoadctRev,-39

title='Pearson/Spearman R: '+string(list.r.r_pearson[0],list.r.r_spearman[0],format='(F0.2,"/",F0.2)')
plot,list.x,list.y,/nodata,xtitle=list.xtitle,ytitle=list.ytitle,/iso,title=title
oplot,list.x,list.y,psym=1,color=keyword_set(PS)?250:5b

max=max(list.x,min=mix)
x=[mix,max]
oplot,x,list.res[0]+list.res[1]*x,thick=2

;may=max(list.y,min=miy)
;y=[miy,may]
;plots,x,y,linestyle=2

end;pro PlotCorrelateXDIGroups
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CorrelateXDIGroups,ev

widget_control,ev.top,get_uvalue=list
ID=widget_info(ev.top,find_by_uname='CorGroup1')
widget_control,ID,get_value=gr1
ID=widget_info(ev.top,find_by_uname='CorGroup2')
widget_control,ID,get_value=gr2

ID=widget_info(ev.top,find_by_uname='CorGroupRtol')
widget_control,ID,get_value=Rtol
ID=widget_info(ev.top,find_by_uname='CorGroupshifttol')
widget_control,ID,get_value=shifttol
ID=widget_info(ev.top,find_by_uname='CorGroupnpix')
widget_control,ID,get_value=npix
Rtol=float(Rtol)
shifttol=float(shifttol)
npix=fix(npix)

gr1=0>fix(gr1)<(list.NMaps-1)
gr2=0>fix(gr2)<(list.NMaps-1)
ind=[gr1,gr2]

x=reform((*list.data)[gr1,*,*],list.MapCol,list.MapRow)
y=reform((*list.data)[gr2,*,*],list.MapCol,list.MapRow)

ind=where(~FINITE(x),ct)
if ct ne 0 then x[ind]=0
ind=where(~FINITE(y),ct)
if ct ne 0 then y[ind]=0

LoadctX,list
ret=image_equal(x,y,outid=list.outid,Rtol=Rtol,npix=npix,shifttol=shifttol,$
    xsize=list.XSize,ysize=list.YSize,quick=dialog_message('Quick correlation?',/question) eq 'Yes')

; Scatter plot + linfit
x-=min(x)
y-=min(y)
x*=100./max(x)
y*=100./max(y)
res=LADFIT(x,y)
L={x:x,y:y,res:res,xtitle:GetGroupLabel(list,gr1),ytitle:GetGroupLabel(list,gr2),r:ret}
window,1
PlotCorrelateXDIGroups,L

; Save result
;if ret.bpass then $
SaveImage,list.path,list.file+'.eps','PlotCorrelateXDIGroups',L,OutID=list.OutID

end;pro CorrelateXDIGroups
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddXDIGroups,ev

widget_control,ev.id,get_uvalue=uval
bdeleteoriginal = 0b
bresetselect=0b
case uval of
'Create Overlap':begin
    widget_control,ev.top,get_uvalue=list
    if list.ScanDim ne 2 then return
    if not ptr_valid(list.data) then return

    ; Get thresholds and group indices
    ID=widget_info(ev.top,find_by_uname='OverlapGroups')
    widget_control,ID,get_value=str
    b=ParseComboCommand(str,commands=cmds,indices=indices,error=error)
    if ~b then begin
        printw,list.OutID,error
        return
    endif

    ID=widget_info(ev.top,find_by_uname='Thresholds')
    widget_control,ID,get_value=thres
    thres=strtrim(strcompress(thres,/remove_all),2)
    thres=strsplit(thres,',',/extract,count=ct)
    tmp=strmid(thres,0,1)
    lthen=tmp eq '<'
    gthen=tmp eq '>'
    lgthen=lthen or gthen
    for i=0l,ct-1 do if lgthen[i] then thres[i]=strmid(thres[i],1)
    thres=float(thres)

    ID=widget_info(ev.top,find_by_uname='RGB')
    widget_control,ID,get_value=RGB
    RGB=strtrim(strcompress(RGB,/remove_all),2)
    RGB=strsplit(RGB,'][,',/extract,count=ct)
    RGB=float(RGB)
    n=n_elements(RGB)
    if (n mod 3) ne 0 then begin
        printw,list.OutID,'Color syntax error.'
        return
    endif
    RGB=reform(RGB,3,n/3)

    ; Make overlap
    block=ImageOverlap(list,cmds,indices,thres,lthen,RGB,nextra=nextra,ncombi=ncombi,error=error)
    if error then return

    combi=2b
    bresetselect=1b
    endcase
'AddGroups':begin
    ; Get group indices
    widget_control,ev.id,get_value=indices

    if ParseIndices(indices,/more) then return

    widget_control,ev.top,get_uvalue=list
    if not ptr_valid(list.data) then return
    tmp=where((indices lt 0) or (indices gt list.NMaps-1),ct)
    if ct ne 0 then return

    ; Make combination(s)
    block=ImageCombine(list,indices,nextra=nextra,ncombi=ncombi,error=error)
    if error then return

    combi=1b
    bresetselect=1b
    endcase
'AddGroupsXDI':begin
    WIDGET_CONTROL, ev.top, get_uvalue=list
    
    ; Append groups
    if list.bopenformat then format=list.openformat
    if list.NMaps eq 0 then return
    n = list.NMaps
    OpenXDI,ev,1,/add,format=format
    WIDGET_CONTROL, ev.top, get_uvalue=list
    if list.NMaps ne 2*n then begin
        if list.NMaps ne n then $
            DeleteXDIGroupsHelper,list,n+indgen(list.NMaps-n)
        return
    endif
    
    ; Make combination(s)
    indices=indgen(1,n)
    indices=[indices,indices+n]
    block=ImageCombine(list,indices,nextra=nextra,ncombi=ncombi,error=error,/pairs)
    if error then begin
        DeleteXDIGroupsHelper,list,n+indgen(n)
        return
    endif
    
    ; Delete all groups
    bdeleteoriginal=1b
    inddelete=indgen(list.NMaps)

    combi=1b
    endcase
'Create Blend Group':begin
    widget_control,ev.top,get_uvalue=list
    block=alpha_blending(list,indices=indices,nextra=nextra,ncombi=ncombi,error=error)
    if error then return

    combi=2b
    bresetselect=1b
    endcase
endcase

; Add new group(s)
(*list.dataorg)=[(*list.dataorg),block]
(*list.data)=(*list.dataorg)
list.NMaps+=nextra

; Add new group description(s)
nrow=nextra
ncol=ncombi

for k=0l,nrow-1 do begin
    pmadd=strarr(list.npmfields)
    npmadd=0
    for i=0l,ncol-1 do begin
        j=indices[i,k]+1

        npmadd+=(*list.pmstruc).(j).n
        tmp=(*list.pmstruc).(j).data

        s=size(tmp)
        if s[1] ne list.npmfields then tmp=[tmp,replicate('-',list.npmfields-s[1],((s[0] le 1)?1:s[2]))]
        pmadd=[[pmadd],[tmp]]
    endfor

    tmp=max(*list.pmtags)+1
    (*list.pmtags)=[(*list.pmtags),tmp]
    (*list.pmstruc)=create_struct((*list.pmstruc),'t'+stringr(tmp),$
        {n:npmadd,data:pmadd[*,1:*],labeli:0b,combi:combi}) ; Only labeli=0 supports grouping when saved/reloaded
endfor

; Delete added groups
if bdeleteoriginal then DeleteXDIGroupsHelper,list,inddelete

; Finish
if bresetselect then list.NSelect=list.NMaps-1
ID=widget_info(ev.top,FIND_BY_UNAME='NRslider')
widget_control,ID,set_value=list.NSelect,set_slider_max=list.NMaps-1

WorkingDataCorrect,list
widget_control,ev.top,set_uvalue=list
RefreshDisplayXDI,list
RefreshWindowXDI,ev
end; pro AddXDIGroups
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_XDIparam,ID

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
end;pro CleanUp_XDIparam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDIparam_event,ev

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
    'nhorizontal':   begin
         reads,val,dval
         list.nhorizontal=1>dval
         endcase
    'charpr':   begin
         reads,val,dval
         list.charpr=dval/100
         endcase
    'boardercol':begin
         reads,val,dval
         list.boardercol=byte(dval)
         endcase
    'labelcol':begin
         reads,val,dval
         list.labelcol=fix(dval)
         endcase
    'labeli':begin
         reads,val,dval
         list.labeli=0>dval<(list.npmfields-1)
         endcase
    'markers':begin
         list.markers=val
         endcase
    'ofor':list.openformat=val
    'buttonFOR':begin
         list.bopenformat=ev.select
         tmp=widget_info(ev.top,FIND_BY_UNAME='ofor')
         widget_control,tmp,sensitive=ev.select
         endcase
    'ctl':begin
         list.ctl.nr=fix(val)
         endcase
    'buttonCT':begin
         list.ctl.type=~ev.select
         ctlID=widget_info(ev.top,FIND_BY_UNAME='ctl')
         widget_control,ctlID,sensitive=list.ctl.type eq 0
         if list.ctl.type eq 1 then begin
             heap_free,list.ctl
             TVLCT, R, G, B, /GET
             list.ctl.R=ptr_new(R)
             list.ctl.G=ptr_new(G)
             list.ctl.B=ptr_new(B)
             widget_control,ID,set_uvalue=list
             XLOADCT, GROUP=ev.top, UPDATECALLBACK='RefreshDisplayXDI_LoadctX', UPDATECBDATA=list,/modal
         endif

         endcase
endcase
printw,list.OutID,'Edit '+uval
widget_control,ID,set_uvalue=list
RefreshDisplayXDI,list
end;pro XDIparam_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDIparam,ev

widget_control,ev.top,get_uvalue=list
widget_control,ev.top,sensitive=0
base=widget_base(title='XDI parameters',uvalue=ev.top,/row)
base1=widget_base(base,/column)
xs=150

baset=widget_base(base1,/row)
label=widget_label(baset,value='Charsize (%xsize)  :        ',xsize=xs)
text=widget_text(baset,/editable,value= $
stringr(list.charpr*100),uvalue='charpr',uname='charpr')

baset=widget_base(base1,/row)
label=widget_label(baset,value='Label-index  :      ',xsize=xs)
text=widget_text(baset,/editable,value= $
stringr(list.labeli),uvalue='labeli',uname='labeli')

baset=widget_base(base1,/row)
label=widget_label(baset,value='Output nx :          ',xsize=xs)
text=widget_text(baset,/editable,value= $
stringr(list.nhorizontal),uvalue='nhorizontal',uname='nhorizontal')

baset=widget_base(base1,/row)
label=widget_label(baset,value='1D plot markers :          ',xsize=xs)
text=widget_text(baset,/editable,value= $
stringr(list.markers),uvalue='markers',uname='markers')

base2=widget_base(base1,/row)
baset=widget_base(base2,/row,/nonexclusive,xsize=xs)
buttonFR=widget_button(baset,value='Open formated: ',uvalue='buttonFOR')
baset=widget_base(base2,/row)
text=widget_text(baset,/editable,value= $
stringr(list.openformat),uvalue='ofor',uname='ofor',sensitive=list.bopenformat)

base2=widget_base(base1,/row)
baset=widget_base(base2,/row,/nonexclusive,xsize=xs)
buttonCT=widget_button(baset,value='Manual color table: ',uvalue='buttonCT')
baset=widget_base(base2,/row)
text=widget_text(baset,/editable,value= $
stringr(list.ctl.nr),uvalue='ctl',uname='ctl',sensitive=list.ctl.type eq 0)

baset=widget_base(base1,/row)
label=widget_label(baset,value='Color index for boarders  :        ',xsize=xs)
text=widget_text(baset,/editable,value= $
stringr(fix(list.boardercol)),uvalue='boardercol',uname='boardercol')

baset=widget_base(base1,/row)
label=widget_label(baset,value='Color index for labels  :        ',xsize=xs)
text=widget_text(baset,/editable,value= $
stringr(fix(list.labelcol)),uvalue='labelcol',uname='labelcol')

button=widget_button(base,value='OK',uvalue='OK')
WIDGET_CONTROL, base, /REALIZE
if list.ctl.type eq 0 then widget_control,buttonCT,set_button=1
if list.bopenformat then widget_control,buttonFR,set_button=1

Xmanager,'XDIparam',base, event_handler='XDIparam_event',cleanup='CleanUp_XDIparam',GROUP_LEADER=ev.top
end;pro XDIparam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro synchrotron_run,X,A,F,PDER,struc=struc

; A: decay | intercept | rico | intercept
; Y: intercept * exp(-decay*x)
;    or
;    rico*x + intercept

n = n_elements(A)
nx = n_elements(X)
F = fltarr(nx)

bpder=N_PARAMS(0) gt 3
if bpder then pder=fltarr(nx,n)

off1=struc.ndecay
for i=0l,struc.ndecay-1 do begin
    b=struc.ind_decay[0,i]
    e=struc.ind_decay[1,i]
    ind=b+indgen(e-b+1)
    F[ind]=A[i+off1]*exp(-A[i]*x[ind])
    
    if bpder then begin
        pder[ind,i]=-x*F[ind]
        pder[ind,i+off1]=exp(-A[i]*x[ind])
    endif
endfor

off0=2*struc.ndecay
off1=off0+struc.nlin
for i=0l,struc.nlin-1 do begin
    b=struc.ind_lin[0,i]
    e=struc.ind_lin[1,i]
    ind=b+indgen(e-b+1)
    F[ind]=A[i+off0]*x[ind]+A[i+off1]
    if bpder then begin
        pder[ind,i+off0]=x[ind]
        pder[ind,i+off1]=1
    endif
endfor

end;pro synchrotron_run
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function synchrotron_run_joints,y,k
n=n_elements(y)

; Find points where derivative high
tmp=deriv(y)
mi=min(tmp,max=ma,/nan)
if n_elements(k) ne 1 then k=3
ind=where(tmp gt k*stdev(tmp),ct)

if ct eq 0 then joints=[[0ul,0],[n-1,n-1]] else begin

    ; Separate different injection points
    if ct eq 1 then diff=[0l] else begin
        diff=((ind-shift(ind,1))<2)-1
        diff=total([0,diff[1:*]],/cum)
    endelse
    
    h=histogram(diff,rev=rev)
    njoints=n_elements(h)
    joints=ulonarr(2,njoints)
    for i=0l,njoints-1 do begin
    
        ; Get end of the previous run and begin of the next
        ni=rev[i+1]-rev[i]
        indi=ind[rev[rev[i]:rev[i+1]-1]]
        indi0=indi[0]-reverse(indgen(ni>2)+1)
        indi1=indi[ni-1]+indgen(ni>2)+1
        m0=mean(y[indi0])+[-1,1]*stddev(y[indi0])*4
        m1=mean(y[indi1])+[-1,1]*stddev(y[indi1])*4
        indi=[indi0,indi,indi1]
                
        tmp=where(y[indi] ge m0[0] and y[indi] le m0[1],ct)
        imin=tmp[ct-1]
        tmp=where(y[indi] ge m1[0] and y[indi] le m1[1],ct)
        imax=tmp[0]

        joints[*,i]=[indi[imin],indi[imax]]
    endfor
    
    if joints[njoints-1] ne n-1 then joints=[[joints],[n-1,n-1]]
    if joints[0] ne 0 then joints=[[0,0],[joints]]
endelse

return,joints
end;function synchrotron_run_joints
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro synchrotron_run_init,x,y,A,struc,k

binit=n_elements(struc) eq 0

if binit then begin
    ; Find joints
    joints=synchrotron_run_joints(y,k)
    njoints=n_elements(joints)/2
    
    ; Initial parameters of decay
    ndecay=njoints-1
    ind_decay=reform(joints[1:2*njoints-2],2,ndecay)
    
    init_decayct=alog(y[ind_decay[0,*]]/y[ind_decay[1,*]])/(x[ind_decay[1,*]]-x[ind_decay[0,*]])
    init_N0=y[ind_decay[0,*]]*exp(init_decayct*x[ind_decay[0,*]])
    init_N0+=y[ind_decay[1,*]]*exp(init_decayct*x[ind_decay[1,*]])
    init_N0/=2
    
    A=[reform(init_decayct),reform(init_N0)]

    ; Initial parameters of linear joints
    ind=where(joints[1,*]-joints[0,*] gt 1,nlin)
    if nlin ne 0 then begin
        ind_lin=joints[*,ind]
        init_rico=(y[ind_lin[1,*]]-y[ind_lin[0,*]])/(x[ind_lin[1,*]]-x[ind_lin[0,*]])
        init_intercept=y[ind_lin[0,*]]-init_rico*x[ind_lin[0,*]]
        init_intercept+=y[ind_lin[1,*]]-init_rico*x[ind_lin[1,*]]
        init_intercept/=2
        
        A=[A,reform(init_rico),reform(init_intercept)]
        ind_lin[0,*]++
        ind_lin[1,*]--
    endif else ind_lin=[-1,-1]
    
    AInfo=MakeAInfo(n_elements(A))
;    off0=2*struc.ndecay
;    off1=off0+2*struc.nlin-1
;    AInfo[off0:off1].fixed=1
        
    struc={ind_decay:ind_decay,ndecay:ndecay,ind_lin:ind_lin,nlin:nlin,AInfo:AInfo}
endif else begin
    if struc.nlin ne 0 then begin
        synchrotron_run,x,A,y,struc=struc
        
        ind_lin=struc.ind_lin
        ind_lin[0,*]--
        ind_lin[1,*]++
        
        init_rico=(y[ind_lin[1,*]]-y[ind_lin[0,*]])/(x[ind_lin[1,*]]-x[ind_lin[0,*]])
        init_intercept=y[ind_lin[0,*]]-init_rico*x[ind_lin[0,*]]
        init_intercept+=y[ind_lin[1,*]]-init_rico*x[ind_lin[1,*]]
        init_intercept/=2
        
        off0=2*struc.ndecay
        off1=off0+struc.nlin
        off2=off1+struc.nlin
        A[off0:off1-1]=init_rico
        A[off1:off2-1]=init_intercept
    endif
endelse

end;pro synchrotron_run_init
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function synchrotron_run_model,current,struc,A,k=k
; current: synchrotron current
s=size(current,/dim)
n=n_elements(current)
y=reform(float(current),n)
x=findgen(n)

if n_elements(struc) eq 0 then begin
    synchrotron_run_init,x,y,A,struc,k
    F=NLLS(x,y,A,'synchrotron_run',weights=y*0+1,struc=struc,AInfo=struc.AInfo)
endif
synchrotron_run_init,x,y,A,struc
synchrotron_run,x,A,F,struc=struc
F=reform(F,s,/overwrite)

return,F
end;function synchrotron_run_model
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro WriteDC,imsfile,datain,OutID

data=datain
; ----Make IMS----
s=size(reform(data))
wDimensions=s[0]
if wDimensions eq 1 then data=transpose(data)
if wDimensions eq 2 then data=reverse(data,2)
writeims, wDimensions, data, 'DC', imsfile
printw,OutID,'Saved '+imsfile
end;pro WriteDC
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadDC,path,ext,top,OutID,error=error,iCounter=iCounter

DC=0.
error=0b

if openr_safe(lun, path)then begin
    printw,OutID,path+' not found'
    error=1B
    return,DC
endif

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    error=1B
ENDIF ELSE BEGIN
    case strlowcase(ext) of
    '.spe':    begin
            block=['$DORIS_CURRENT:','$MEAS_TIM0:']
            if n_elements(iCounter) eq 0 then iCounter=PromptNumber(block,top,'Select counter:',/choice)
            block=block[iCounter]
            endcase
    '.fio': block='%c'
    '.mca':    block='#L'
    else:    block='FORMAT' ; SMART frame
    endcase
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if (strmid(line,0,len) eq block) then begin
        case strlowcase(ext) of
        '.spe':    begin
                case iCounter of
                0:    begin
                    in=fltarr(3)
                    readf,lun,in
                    DC=in[0]
                    endcase
                1:     begin
                    ; Read live time
                    Dc=0.
                    readf,lun,DC
                    endcase
                endcase
                endcase
        '.fio': begin
                block='  Gain'
                n=strlen(block)
                keep=''
                line=''
                while strmid(line,0,n) ne block and not eof(lun) do begin
                    keep+=line
                    readf,lun,line
                endwhile
                line=keep

                tmp=strsplit(line,' ,',/extract,count=ct)
                if n_elements(iCounter) eq 0 then begin
                    iCounter=PromptNumber(tmp[2*indgen(ct/2)],top,'Select counter:',/choice)
                    iCounter=2*iCounter+1
                endif
                DC=float(tmp[iCounter])

                endcase
        '.mca':    begin
                line2=''
                readf,lun,line2
                
                tmp=strsplit(line2,' ',/extract,count=ct)
                if n_elements(iCounter) eq 0 then begin
                    tmp2=strsplit(line,' ',/extract)
                    iCounter=PromptNumber(tmp2[1:ct-1],top,'Select counter:',/choice) ; Removes #L
                endif
                DC=float(tmp[iCounter])
                
                endcase
        else:    begin
                ; Get live time
                line=STRCOMPRESS(line)
                pos = STRPOS(line, 'ELAPSDA')
                DC=1.
                if pos eq -1 then error=1b else begin
                    DC=0.
                    reads,STRMID(line, pos+8,72),DC
                endelse
                endelse
        endcase

    endif else error=1B

ENDELSE

free_lun,lun

return,DC
end;function ReadDC
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PlotDCCor,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_uvalue=list2
widget_control,ev.top,get_uvalue=list

if n_elements((*list.filter)) eq 0 then return
i=*list.filter

case list2.ScanDim of
1:    begin
    window,0
    plot,list2.cor[i,*],/xs,/ys
    endcase
2:    begin
    s=Dimsize(list2.cor,3)
    xs=s[1]*3
    ys=s[2]*3
    window,0,xsize=xs,ysize=ys
    tvscl,rebin(list2.cor[i,*,*],1,xs,ys,/sample),/order
    endcase
endcase

end;pro PlotDCCor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DCNormal_cleanup,ID

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
    heap_free,list
endif

end;pro DCNormal_cleanup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DCNormal_refresh,ev,listmain,bsave=bsave
; sets uvalue of ev.top

widget_control,listmain.main,get_uvalue=list

; ----Flux----
if listmain.busemodel then cor_apply=*listmain.cormodel $
else cor_apply=(*listmain.cor)/(*listmain.cormodelorg)*(*listmain.cormodel)
ind=where(~finite(cor_apply),ct)
if ct ne 0 then cor_apply[ind]=1

; ----Plot flux----
if n_elements(ev) ne 0 then begin
    Wset,listmain.drawindex
    
    xrange=listmain.xrange
    xrange=[xrange[0]<xrange[1],xrange[0]>xrange[1]]
    if xrange[0] eq xrange[1] then xrange=[0,n_elements(cor_apply)-1] else psym=-1
    
    LoadctRev,-39,/silent
    mi=min((*listmain.cor)[xrange[0]:xrange[1]],max=ma)
    yrange=[mi,ma]
    mi=min((*listmain.cormodel)[xrange[0]:xrange[1]],max=ma)
    yrange=[mi<yrange[0],ma>yrange[1]]
    mi=min(cor_apply[xrange[0]:xrange[1]],max=ma)
    yrange=[mi<yrange[0],ma>yrange[1]]
    
    plot,*listmain.cor,ytitle='Flux',xtitle='Filenr',xrange=xrange,yrange=yrange,/xs,psym=psym
    oplot,*listmain.cormodel,color=200
    oplot,cor_apply,color=10
    
    yoff=0.8
    xoff=0.5
    d=0.05
    dd=0.03
    plots,[xoff,xoff+d],[yoff,yoff],/normal
    xyouts,xoff+d*1.2,yoff,'Flux',/normal
    yoff+=dd
    plots,[xoff,xoff+d],[yoff,yoff],/normal,color=200
    xyouts,xoff+d*1.2,yoff,'Modelled Flux',/normal,color=200
    yoff+=dd
    plots,[xoff,xoff+d],[yoff,yoff],/normal,color=10
    xyouts,xoff+d*1.2,yoff,'Used Flux',/normal,color=10
    
    GetSysVar,listmain.sysvar
    widget_control,ev.top,set_uvalue=listmain
    Wset,listmain.pixindex
    Device, Copy = [0,0, listmain.xsize,listmain.ysize, 0,0,listmain.drawindex]
endif

; ----Save flux----
if keyword_set(bsave) then begin
    path=SelFile(list.path,'flux.xdi','*.xdi','Save Flux....')
    
    coord=getxdicoord(list,list.coord)
    unit=list.unit eq 'deg'
    if list.ScanDim eq 2 and list.TRAN then begin
        coord=reverse(coord,2)
        unit=reverse(unit)
    endif
    
       if path ne '' then Type3XDI,path,cor_apply,['Flux'],list.ScanDim,list.OutID,$
           TRAN=list.TRAN,REV1=list.REV1,REV2=list.REV2,REV3=list.REV3,ct=list.expo,$
           ounit=unit[0],unit=unit,coord=coord,ocoord=coord,$
           Hsign=list.Hsign,Vsign=list.Vsign
       ;WriteDC,path,cor_apply,list.OutID
endif

; ---Normalize flux----
cor_apply/=listmain.normalize
ind=where(~finite(cor_apply),ct)
if ct ne 0 then cor_apply[ind]=1

; ----Get correction----
cor_apply=1./cor_apply
ind=where(~finite(cor_apply),ct)
if ct ne 0 then cor_apply[ind]=1

; ----Order like xdi orientation----
ReformXDIDatablock,cor_apply,[0b,list.TRAN],[0b,list.REV1],[0b,list.REV2],[0b,list.REV3],0b

; ---Show correction----
window,0,xsize=list.XSize,ysize=list.YSize,title='Correction (data is multiplied by this)'
case list.ScanDim of
2:    begin
    LoadctX,list
    tvscl,congrid(reform(cor_apply[0,*,*],list.MapCol,list.MapRow),list.XSize,list.YSize),/ORDER
    endcase
1:    begin
    LoadctRev,0,/silent,rev=-1
    n=n_elements(cor_apply)

    x=indgen(n)
    xrange=[0,n-1]
    if list.TRAN then BOOL=list.REV2 else BOOL=list.REV1
    if BOOL then begin
        x=reverse(x)
        xrange=reverse(xrange)
    endif

    if list.TRAN then plot,cor_apply,x,yrange=reverse(xrange),xstyle=1,ystyle=1,ytitle='Spatial distance (arbitrarily)',xcharsize=.7,$
             xtitle='I',ycharsize=.7,linestyle=0 $
    else plot,x,cor_apply,xrange=xrange,xstyle=1,ystyle=1,xtitle='Spatial distance (arbitrarily)',xcharsize=.7,$
             ytitle='I',ycharsize=.7,linestyle=0

    endcase
endcase

; ----Apply correction----
ptr=list.dataorg
*ptr=*listmain.dataorg
for i=0l,list.NMaps-1 do (*ptr)[i,*,*]*=cor_apply
*list.data=*ptr

WorkingDataCorrect,list
RefreshDisplayXDI,list

end;pro DCNormal_refresh
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DCNormal_DRAW,ev

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
       Event_Pro='DCNormal_EVENT'

    list.xrange=[x1,x2]
    DCNormal_refresh,ev,list
    return
endif

y2 = ROIy[0] > coords[1,1] < ROIy[1]
WSet, list.drawindex
Arrow, x1, y2, x2, y2, /Data, /Solid, HSize=12
PlotS, [x2, x2], ROIy
end;pro DCNormal_DRAW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DCNormal_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=listmain
case widget_info(ev.id,/type) of
0:    begin
    IF TAG_NAMES(ev, /STRUCTURE_NAME) EQ 'WIDGET_KILL_REQUEST' THEN bcancel=1b else begin
        if listmain.listtype eq 60 then begin
            listmain.listtype=61
            widget_control,listmain.draw,get_value=drawindex,draw_button=1
            listmain.drawindex=drawindex
        endif
    endelse
    endcase
4:    begin
    possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
    thisEvent = possibleEventTypes[ev.type]
    IF thisEvent NE 'DOWN' THEN RETURN
    
    listmain.xs=ev.x
    listmain.ys=ev.y
    
    WSet, listmain.drawindex
    SetSysVar,listmain.sysvar
    coords = Convert_Coord(ev.x, ev.y, /Device, /To_Data)
    GetCRange,ROIx=ROIx,ROIy=ROIy
    x1 = ROIx[0] > coords[0] < ROIx[1]
    y1 = ROIy[0] > coords[1] < ROIy[1]
    PlotS, [x1, x1], ROIy
    
    h=!D.x_size
    v=!D.y_size
    WSet, listmain.pixindex
    Device, Copy = [0,0, h,v, 0,0,listmain.drawindex]
    
    Widget_Control, ev.id, Event_Pro='DCNormal_DRAW', draw_motion=1
    widget_control,ev.top,set_uvalue=listmain
    return
    endcase
1:    begin
    widget_control,ev.id,get_value=val
    case val of
    'OK':     begin
            bcancel=0b
            bsave=1b
            endcase
    'Cancel': bcancel=1b
    'Use modelled values':listmain.busemodel=ev.select
    'Save Flux ...':bsave=1b
    'Default':    begin
                cormodel=synchrotron_run_model((*listmain.cor),struc,A,k=promptnumber('3.',ev.top,'Jump search: k times the stdev'))
                *listmain.cormodel=cormodel
                *listmain.cormodelorg=cormodel
                *listmain.struc=struc
                *listmain.A=A
                ID=widget_info(ev.top,find_by_uname='table')
                widget_control,ID,set_value=transpose(reform(A[0:2*struc.ndecay-1],struc.ndecay,2))
                endcase
    'Add':     begin
            widget_control,ev.id,get_uvalue=uval
            sel=widget_info(uval.id,/TABLE_SELECT);[ left, top, right,  bottom ]
            ; Not implemented yet

            endcase
    'Delete':     begin
                widget_control,ev.id,get_uvalue=uval
                sel=widget_info(uval.id,/TABLE_SELECT);[ left, top, right,  bottom ]
                ; Not implemented yet
                
                endcase
    endcase
    endcase
3:    begin
    widget_control,ev.id,get_uvalue=uval
    widget_control,ev.id,get_value=val
    case uval of
    'normalize': listmain.normalize=float(val)
    endcase
    endcase
9:    begin

    event=TAG_NAMES(ev, /STRUCTURE_NAME)
    if event eq 'WIDGET_CONTEXT' then begin
        ID=WIDGET_INFO(ev.top, FIND_BY_UNAME = 'ContextTable')
        widget_control,ID,set_uvalue={ID:ev.id,top:ev.top}
        WIDGET_DISPLAYCONTEXTMENU, ev.ID, ev.X, ev.Y, ID
        return
    endif

    widget_control,ev.id,get_value=Ain
    A=*listmain.A
    A[0:2*(*listmain.struc).ndecay-1]=transpose(Ain)
    
    cormodel=synchrotron_run_model((*listmain.cor),*listmain.struc,A)
    *listmain.A=A
    *listmain.cormodel=cormodel
    endcase
else: 
endcase

DCNormal_refresh,ev,listmain,bsave=bsave
if n_elements(bcancel) eq 0 then return

widget_control,ev.top,get_uvalue=listmain

; ----Keep normalization?----
widget_control,listmain.main,get_uvalue=list
if bcancel then begin
    printw,list.OutID,'Normalization canceled.'
    *list.dataorg=*listmain.dataorg
    *list.data=*listmain.dataorg
endif else begin
    printw,list.OutID,'Normalization done.'
    list.normflag=1b
    WIDGET_CONTROL, listmain.main, set_uvalue=list
endelse

WorkingDataCorrect,list
RefreshDisplayXDI,list

widget_control,ev.top,/destroy
end;pro DCNormal_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DCNormal_interactive,cor,list

mean=mean(cor)

binteractive=dialog_message('Model with synchrotron runs?',/question,/default_no) eq 'Yes'
if binteractive then begin
    cormodel=synchrotron_run_model(cor,struc,A)
    
    result=ReadINI('$ChiSize:')
    xsize=(result.(0))[0]
    ysize=(result.(0))[1]
    Topbase = WIDGET_BASE(/column,title='Normalization to flux (divide by the flux)',/TLB_KILL_REQUEST_EVENTS)
    baset=widget_base(Topbase,/column)
            draw = WIDGET_DRAW(baset, XSIZE=xsize, YSIZE=ysize)
            Window, /Free, xsize=xsize, ysize=ysize, /Pixmap
            pixindex = !D.Window
    basedynn=widget_base(Topbase,/row)
    baset=widget_base(basedynn,/column)
            b=widget_table(baset,value=transpose(reform(A[0:2*struc.ndecay-1],struc.ndecay,2)),$
            column_labels=['Decay','Intercept'],uname='table',/editable,xsize=2,CONTEXT_EVENTS=0)
            b=widget_button(baset,value='Default')
            
    basedyn=widget_base(basedynn,/column)
            baset=widget_base(basedyn,/row,/nonexcl)
            b=widget_button(baset,value='Use modelled values')
            baset=widget_base(basedyn,/row)
            b=widget_label(baset,value='Flux normalized to:')
            b=widget_text(baset,value=string(mean,format='(f0.0)'),uvalue='normalize',/editable)
            baset=widget_base(basedyn,/row)
            b=widget_button(baset,value='OK')
            b=widget_button(baset,value='Cancel')
            b=widget_button(baset,value='Save Flux ...')
    
    baset = WIDGET_BASE(Topbase,  /CONTEXT_MENU, UNAME='ContextTable')
        button=widget_button(baset,value='Add')
        button=widget_button(baset,value='Delete')
    
    listmain={listtype:60,$
            cor:ptr_new(cor),$
            cormodelorg:ptr_new(cormodel),$
            cormodel:ptr_new(cormodel),$
            struc:ptr_new(struc),$
            A:ptr_new(A),$
            normalize:mean,$
            busemodel:0b,$
            dataorg:ptr_new(*list.dataorg),$
            pixindex:pixindex,$
            xrange:[0L,0],$
            draw:draw,$
            drawindex:0L,$
            xsize:xsize,$
            ysize:ysize,$
            xs:0,$
            ys:0,$
            sysvar:PointSysVar(),$    ; System variable for coord. convertion
            main:list.top}
    
    WIDGET_CONTROL,Topbase, /REALIZE,set_uvalue=listmain
    widget_control,Topbase,send_event={ID:Topbase,TOP:Topbase,HANDLER:Topbase,SELECT:long(1)}
    Xmanager,'DCNormal_interactive',Topbase,/NO_BLOCK,event='DCNormal_event',cleanup='DCNormal_cleanup'
endif else begin
    p=ptr_new(cor)
    listmain={listtype:60,$
            cor:p,$
            cormodelorg:p,$
            cormodel:p,$
            struc:ptr_new(),$
            A:ptr_new(),$
            normalize:mean,$
            busemodel:1b,$
            dataorg:list.dataorg,$
            pixindex:0L,$
            xrange:[0L,0],$
            draw:0L,$
            drawindex:0L,$
            xsize:0L,$
            ysize:0L,$
            xs:0,$
            ys:0,$
            sysvar:ptr_new(),$
            main:list.top}
    DCNormal_refresh,ev,listmain,bsave=1
    ptr_free,p
endelse

end;pro DCNormal_interactive
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DCNormal,ev,type

widget_control,ev.top,get_uvalue=list
if list.normflag then begin
    printw,list.OutID,'Already normalized!'
    return
endif

bgroups=1b
case type of
0:     begin
    ext='.SPE'
    bgroups=0b
    endcase
1:     ext='.ims'
2:     begin
    ext='.fio'
    bgroups=0b
    endcase
3:     ext='.spec'
4:     ext='.tiff'
5:    ext='.xdi'
6:    begin
    ext='.*' ;FRAME header normalization
    bgroups=0b
    endcase
7:     begin
    ext='.MCA'
    bgroups=0b
    endcase
8:    ext='.tiff'
9:    ext='.edf'
10: ext='.edf'
endcase

if bgroups then begin

    ; ----Read DC file----
    t=CutPath(list.file,file=file)
    imsfile=''
    imsfile=DIALOG_PICKFILE(path=list.path,file=file+ext,filter='*'+ext,title='Read file...')
    if imsfile EQ '' then return

    orientation=bytarr(4)
    case type of
    1:    cor=readims(imsfile, roi_names, ndimensions)
    3:    cor=readspec(imsfile,ev.top, roi_names, ndimensions,orientation)
    4:    begin
        cor=read_tiff(imsfile)
        s=Dimsize(cor,2)
        cor=reform(cor,1,s[0],s[1],/overwrite)
        roi_names=['cor']
        endcase
    5:    begin
        tmp=CutPath(imsfile,path=path,file=file,ext=ext)
        if list.bopenformat then format=list.openformat
        struct=ReadXDI(path,file+ext,list.OutID,error=error,format=format)
        cor=struct.data
        orientation=struct.orientation
        n=n_tags(struct.pmstruc)-1
        roi_names=strarr(n)
        for i=1,n do roi_names[i-1]=struct.pmstruc.(i).data[0>list.labeli<(n_elements(struct.pmstruc.(i).data)-1)]
        endcase
    8:    begin
        cor=read_tiff(imsfile)
        s=Dimsize(cor,2)
        
        t=CutPath(imsfile,path=path,file=file)
        imsfile=DIALOG_PICKFILE(path=path,file=file+ext,filter='*'+ext,title='Read ICR...')
        if imsfile EQ '' then return
        ICR=read_tiff(imsfile)
        if ~array_equal(Dimsize(ICR,2),s) then return
        
        t=CutPath(imsfile,path=path,file=file)
        imsfile=DIALOG_PICKFILE(path=path,file=file+ext,filter='*'+ext,title='Read OCR...')
        if imsfile EQ '' then return
        OCR=read_tiff(imsfile)
        if ~array_equal(Dimsize(OCR,2),s) then return
        
        OCR/=float(ICR)
        ind=where(~finite(OCR) or OCR eq 0,ct)
        if ct ne 0 then OCR[ind]=1
        
        cor*=OCR
        cor=reform(cor,1,s[0],s[1],/overwrite)
        roi_names=['cor']
        endcase
    9:    begin
        cor=readedf(imsfile)
        s=Dimsize(cor,2)
        cor=reform(cor,1,s[0],s[1],/overwrite)
        roi_names=['cor']
        endcase
    endcase
    
    ; ----Set to default orientation and check dimensions----
    ReformXDIDatablock,cor,[orientation[0],0b],[orientation[1],0b],[orientation[2],0b],[orientation[3],0b],0b
    GetXDIdatadim,list,nmaps,ncol,nrow,/default
    s=DimSize(cor,3)
    bool=(s[1] eq ncol) and (s[2] eq nrow)
    if ~bool then begin
        bool=(s[2] eq ncol) and (s[1] eq nrow)
        if bool then begin
            ; Transposed
            s[1]=ncol
            s[2]=nrow
            cor=transpose(cor,[0,2,1])
        endif else begin
            bool=ncol*nrow eq s[1]*s[2]
            if bool then begin
                ; Reformed
                s[1]=ncol
                s[2]=nrow
                cor=reform(cor,s[0],s[1],s[2],/overwrite)
            endif
        endelse
    endif

    if ~bool then begin
        printw,list.OutID,'Wrong dimensions in '+imsfile
        return
    endif

    ; ----Select flux counter----
    i=PromptNumber(roi_names,ev.top,'Select Flux counter:',/choice,username='Plot',userpro='PlotDCCor',userdata={cor:cor,ScanDim:list.SCanDim})
    if i eq -1 then return
    cor=cor[i,*,*]

endif else begin

    ; ----Prompt for SPE files----
    Result = DIALOG_MESSAGE('Use filenames from the xdi file?',/question)
    widget_control,/hourglass
    t=CutPath(list.file,file=file)
    if (result eq 'Yes') and PTR_VALID(list.filenames) then begin
        tmp=CutPath(*list.filenames,path=path2,file=path,ext=ext2)
        if type eq 6 then ext=ext2
        n=n_elements(path)

        path3=DIALOG_PICKFILE(path=path2[0],title='Select Directory of files...',/DIRECTORY)
        if path3 eq '' then return
        if (path3 ne path2[0]) then path2=replicate(path3,n)
        path=temporary(path2)+path+ext2
    endif else begin
        path=Select_Files(ev.top,{path:list.path,file:file+ext,Format:ext,Platform:0},$
            sortmethod=list.sortmethod,separator=list.sortseparator)
        if path[0] eq '' then begin
            printw,list.OutID,'No files found.'
            return
        endif
        n=n_elements(path)
    endelse
    
    ; ----Make cor in default orientation----
    GetXDIdatadim,list,nmaps,ncol,nrow,/default
    cor=fltarr(1,ncol,nrow)

    if n ne n_elements(cor) then begin
          printw,list.OutID,'Wrong number of files.'
           n<=n_elements(cor)
    endif

    ; ----Read *.SPE----
    progressBar = Obj_New("PROGRESSBAR",title="Reading files...")
    progressBar -> Start
    for i=0ul,n-1 do begin
        IF progressBar -> CheckCancel() THEN BEGIN
            progressBar -> Destroy
            printw,list.outid,'Process CANCELLED!'
            RETURN
        ENDIF
    
        cor[i]=ReadDC(path[i],ext,ev.top,list.OutID,error=error,iCounter=iCounter)
        if error then begin
              printw,list.OutID,'Error in reading '+path[i]
              return
        endif else printw,list.OutID,path[i]
        
        progressBar -> Update, (i+1.)/n*100
    endfor
    progressBar -> Destroy
    
endelse

DCNormal_interactive,cor,list

end;pro DCNormal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDI2ATHENA,path,file,edffile,OutID,format=format

; ----Initialize----
struct=ReadXDI(path,file,OutID,error=error,format=format)
if error then return

; ----Make DAT files----
s=dimsize(struct.data,3)
nEl=s[0]
nEnergy=s[1]
nrep=s[2]

; Header
header=['ENERGY','SUM REP','SUM REP/I0',reform([replicate('REP',1,nrep),replicate('I0',1,nrep)],2*nrep)]+'('+stringr(lindgen(2*nrep+3)+1)+')'

; Get energy and I0 from spec file
specfile=dialog_pickfile(path=path,filter=['*.spec','*.dat'],title='Select spec file...')
if specfile[0] eq '' then return

; Read spec
n=GetGroupindices(nrep,0l,indices=ind)
if n ne nrep then begin
    printw,OutID,"Number of scans should be "+stringr(nrep)
    return
endif
ind=ind[sort(ind)]

for i=0l,n-1 do begin
    spec=readspec(specfile, 0l, roi_names, ndimensions,orientation,scannumber=ind[i],lun=lun)
    if ndimensions ne 1 then begin
        free_lun,lun
        printw,OutID,"Spec file scan doesn't have the correct dimensions"
        return
    endif

    if i eq 0 then begin
        iE=PromptNumber(roi_names,0l,'Select Energy:',/choice)
        iI0=PromptNumber(roi_names,0l,'Select I0 counter:',/choice)
        energy=spec[iE,*]
        if n_elements(energy) ne nEnergy then begin
            free_lun,lun
            printw,OutID,"Number of energy points should be "+stringr(nEnergy)
            return
        endif
        I0=spec[iI0,*]
    endif else begin
        I0=[I0,spec[iI0,*]]
    endelse
endfor
if n_elements(lun) ne 0 then free_lun,lun

; Write .dat
ind=lindgen(1,nrep)
ind=reform([ind,ind+nrep],2*nrep)
for i=0l,nEl-1 do begin
    wIrname=strjoin(strsplit(struct.pmstruc.(i+1).data[1],'/\',/extract),'_')
    edffilei=edffile+stringr(i)+'_'+wIrname+'.dat'
    
    tmp=transpose(reform(struct.data[i,*,*],nEnergy,nrep))
    tot=reform(total(tmp,1),1,nEnergy)
    totcor=reform(total(tmp/I0,1),1,nEnergy)
    tmp=[tmp,I0]
    
    ok=WriteAthenaDat(edffilei,header,[energy,tot,totcor,tmp[ind,*]])
    if ok then printw,OutID,'Saved '+edffilei $
    else printw,OutID,'Not saved '+edffilei
endfor
end;pro XDI2ATHENA
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDI2IMS,path,file,imsfile,OutID,format=format

; ----Initialize----
struct=ReadXDI(path,file,OutID,error=error,format=format)
if error then return

; ----Make IMS----
s=dimsize(struct.data,3)
nEl=s[0]
wIrnames=strarr(nEl)
for i=1,nEl do wIrnames[i-1]=struct.pmstruc.(i).data[1]
writeims, struct.ScanDim, struct.data, wIrnames, imsfile
printw,OutID,'Saved '+imsfile

end;pro XDI2IMS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDI2TIFF,path,file,edffile,OutID,format=format

; ----Initialize----
struct=ReadXDI(path,file,OutID,error=error,format=format)
if error then return

; ----Make EDF files----
s=dimsize(struct.data,3)
nEl=s[0]

for i=0l,nEl-1 do begin
    wIrname=struct.pmstruc.(i+1).data[0]+'_'+struct.pmstruc.(i+1).data[1]
    edffilei=edffile+stringr(i)+'_'+wIrname+'.tiff'
    
    write_tiff_safe,edffilei,reform(struct.data[i,*,*],s[1],s[2]),error=error
    if ~error then printw,OutID,'Saved '+edffilei
endfor
end;pro XDI2TIFF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDI2EDF,path,file,edffile,OutID,format=format

; ----Initialize----
struct=ReadXDI(path,file,OutID,error=error,format=format)
if error then return

; ----Make EDF files----
s=dimsize(struct.data,3)
nEl=s[0]

for i=0l,nEl-1 do begin
    wIrname=struct.pmstruc.(i+1).data[0]+'_'+struct.pmstruc.(i+1).data[1]
    edffilei=edffile+stringr(i)+'_'+wIrname+'.tiff'
    if WriteEDF(edffilei,reform(struct.data[i,*,*],s[1],s[2])) then $
        printw,OutID,'Saved '+edffilei
endfor
end;pro XDI2EDF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Type3XDI,xdifile,block,roi_names,ndimensions,OutID,$
    TRAN=TRAN,REV1=REV1,REV2=REV2,REV3=REV3,ct=ct,ounit=ounit,unit=unit,$
    coord=coord,ocoord=ocoord,Hsign=Hsign,Vsign=Vsign,nodefaultorient=nodefaultorient
; block in default orientation unless nodefaultorient is set

; ----Initialize----
sblock=Dimsize(block,3)
block=reform(block,sblock,/overwrite)
if n_elements(TRAN) eq 0 then TRAN=0
if n_elements(REV1) eq 0 then REV1=0
if n_elements(REV2) eq 0 then REV2=0
if n_elements(REV3) eq 0 then REV3=0
if n_elements(Hsign) eq 0 then Hsign=1
if n_elements(Vsign) eq 0 then Vsign=1
if n_elements(ct) eq 0 then ct=0
if n_elements(ounit) eq 0 then ounit=0
if n_elements(unit) eq 0 then unit=intarr(2)
if n_elements(ocoord) eq 0 then ocoord=fltarr(2)
if n_elements(coord) eq 0 then coord=fltarr(2,2)

; ----Make structure----
list={ ptrModel:CreateModel(),$
        TRAN:TRAN,$
        REV2:REV2,$
        REV1:REV1,$
        REV3:REV3,$
        Hsign:Hsign,$
        Vsign:Vsign,$
        scandim:ndimensions,$
        nFiles:0L,$
        nrtif:0L,$
        CR:lonarr(2),$
        N:0L,$
        ct:ct,$
        unit1:ounit,$
        unit:unit,$
        FASTMODE:0b,$
        NORMROIs:0b,$
        PTR_files:PTR_NEW(),$
        coord:coord,$
        coord1:ocoord}

ptrSum=(*list.ptrModel).ptrSum
ptrCur=ptrSum

; ----Make xdi----
nmaps=sblock[0]
ncol=sblock[1]
nrow=sblock[2]

(*ptrSum).n=[0,0,0,nmaps]
(*list.ptrModel).n=nmaps
list.nfiles=ncol*nrow
list.nrtif=list.nfiles
list.N=list.nfiles
list.CR=[ncol,nrow]

if n_elements(roi_names) ne nmaps then roi_names='ROI'+string(indgen(nmaps),format='(I0)')
for i=0l,nmaps-1 do begin
    ; type 3: General purpose (e.g. XRF)
    (*ptrCur).next=PTR_NEW({type:3,$
                        group:i,$
                        name:roi_names[i],$
                        sum:ptr_new(reform(block[i,*,*],list.nfiles)),$
                        next:PTR_NEW()})
    ptrCur=(*ptrCur).next
endfor
if SaveXDI(list,xdifile,nodefaultorient=nodefaultorient) eq 1 then printw,OutID,'Saved '+xdifile

; ----Cleanup----
DestroyModel,list.ptrModel,/full

end;pro Type3XDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EDF2XDI,xdifile,edffile,OutID,doreform=doreform

block=readedf(edffile)

s=dimsize(block,2)
n=s[0]*s[1]
if n_elements(doreform) eq 0 then begin
    seprows=0b
    if s[0] ne 1 and s[1] ne 1 then $
        seprows=DIALOG_MESSAGE('Separate rows?',/question,/default_no) eq 'Yes'
endif else seprows=doreform.seprows
if seprows then s=[s,1] else s=[1,s]

if n_elements(doreform) eq 0 then $
    breset=DIALOG_MESSAGE('Reset dimensions?',/question,/default_no) eq 'Yes' $
else breset=doreform.breset
if breset then begin
    title='First dimension (total='+stringr(n)+'):'
    repeat s[1]=ulong(promptnumber(stringr(s[1]),0,title[0])) until n mod s[1] eq 0
    s[2]=n/s[1]
endif

if n_elements(doreform) eq 0 then $
    doreform={seprows:seprows,breset:breset}

block=reform(block,s,/overwrite)
if s[1] eq 1 or s[2] eq 1 then ndimensions=1 else ndimensions=2

temp=CutPath(edffile,file=file)
roi_names=file
if s[0] gt 1 then roi_names+=stringr(indgen(n))

Type3XDI,xdifile,block,roi_names,ndimensions,OutID

end;pro EDF2XDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PyMCADAT2XDI,xdifile,imsfile,OutID

block=readPyMCADAT(imsfile, roi_names, ndimensions)
block>=0
if ndimensions eq 0 then return
if ndimensions eq 2 then begin
    s=dimsize(block,3)
    if s[1] eq 1 then begin
        ndimensions=1
        block=reform(block,s[0],s[2],/overwrite)
    endif
    if s[2] eq 1 then begin
        ndimensions=1
        block=reform(block,s[0],s[1],/overwrite)
    endif
endif
Type3XDI,xdifile,block,roi_names,ndimensions,OutID

end;pro PyMCADAT2XDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro IMS2XDI,xdifile,imsfile,OutID

block=readims(imsfile, roi_names, ndimensions)
block>=0 ; Do this or not?
if ndimensions eq 0 then return
;TODO convert 2D line to 1D
Type3XDI,xdifile,block,roi_names,ndimensions,OutID

end;pro IMS2XDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MICROXASFASTXRF2XDI,xdifile,imsfile,OutID

block=readMICROXASFASTXRF(imsfile, roi_names, ndimensions, scandim)
if ndimensions eq 0 then return

ct=max(scandim[3,*])

coord1=scandim[0:1,0]
Hsign=2*(coord1[1] gt coord1[0])-1

if ndimensions gt 1 then begin
    coord=scandim[0:1,*]
    Vsign=2*(coord[1,1] gt coord[0,1])-1
endif

Type3XDI,xdifile,block,roi_names,ndimensions,OutID,ct=ct,coord=coord,ocoord=ocoord,Hsign=Hsign,Vsign=Vsign

end;pro MICROXASFASTXRF2XDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TIFF2XDI,xdifile,files,OutID

n=n_elements(files)
block=read_tiff(files[0],INTERLEAVE=2)

s=size(block)
case s[0] of
0: return
1:
2:
3:    block=block[*,*,0]/3. + block[*,*,1]/3. + block[*,*,2]/3.
else: return
endcase

s2=dimsize(block,2)
ncol=s2[0]
nrow=s2[1]
block=reform(block,1,ncol,nrow,/overwrite)

for i=1,n-1 do begin
    add=read_tiff(files[i],INTERLEAVE=2)
    if s[0] eq 3 then add=0.3*add[*,*,0] + 0.59*add[*,*,1] + 0.11*add[*,*,2]
    s3=dimsize(add,2)
    if ~array_equal(s2,s3) then stop
    if array_equal(s2,s3) then $
        block=[block,reform(add,1,ncol,nrow)]
endfor

tmp=CutPath(files,file=files)
Type3XDI,xdifile,block,files,2,OutID
end;pro TIFF2XDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Type3CHI,chifile,block,list,x=x

tmp=cutpath(chifile,path=path,file=file,ext=ext)
chifile=path+file
xdifile=list.path+list.file

n=list.Nmaps-1
for i=0l,n do begin
    labeli=0>list.labeli<(n_elements(((*list.pmstruc).(i+1)).data)-1)
    file=chifile+stringr(((*list.pmstruc).(i+1)).data[labeli])+'_'+stringr(i)
    tmp=file_search(file+'*'+ext,count=ct)
    if ct ne 0 then file+=stringr(ct)+ext else file+=ext
    y=reform((block)[i,*,*])
    ny=n_elements(y)
    if n_elements(x) ne ny then x=indgen(ny)

    temp=WriteChi(file,transpose([[x],[y]]),1,xdifile,list.OutID,[0,2*!pi])
endfor
end;pro Type3CHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDI2CHI,ev

widget_control,ev.top,get_uvalue=list
if (list.ScanDim ne 1) then begin
    Result = DIALOG_MESSAGE('XDI has wrong dimensions.',/information)
    return
endif

; ----Select file and open----
chifile=SelFile(list.path,list.file,ChiFormat(),'Save To....')
if chifile eq '' then return

n=list.MapCol*list.MapRow
xrange=list.coord[*,0]
if xrange[0] eq xrange[1] then xrange=[0,n-1.] else $
if list.unit[0] ne 'deg' then xrange*=1000
xrange=getxdicoord(list,xrange)
x=MakeRange(xrange[0],xrange[1],(xrange[1]-xrange[0])/(n-1))

Type3CHI,chifile,*list.data,list,x=x

end;pro XDI2CHI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Type3TIFF,file,blocki,outid,x=x
block=blocki
s=dimsize(block,3)
n=s[1]*s[2]
block=reform(block,s[0],n,/overwrite)
block=[[(n_elements(x) ne n)?indgen(n):x],[transpose(block)]]

result=WriteChi(file,block,1,'...',outid,[0,2*!pi])
end;pro Type3TIFF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDI2CHIb,ev

widget_control,ev.top,get_uvalue=list
if (list.ScanDim ne 1) then begin
    Result = DIALOG_MESSAGE('XDI has wrong dimensions.',/information)
    return
endif

; ----Select file and open----
chifile=SelFile(list.path,list.file,'*.tiff','Save To....')
if chifile eq '' then return

n=list.MapCol*list.MapRow
xrange=list.coord[*,0]
if xrange[0] eq xrange[1] then xrange=[0,n-1.] else $
if list.unit[0] ne 'deg' then xrange*=1000
xrange=getxdicoord(list,xrange)
x=MakeRange(xrange[0],xrange[1],(xrange[1]-xrange[0])/(n-1))

Type3TIFF,chifile,*list.data,list.outid,x=x

end;pro XDI2CHIb
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDIconvert,ev,type
case type of
0:    begin
    extin='.xdi'
    extout='.ims'
    endcase
1:    begin
    extin='.ims'
    extout='.xdi'
    endcase
2:    begin
    extin='.edf'
    extout='.xdi'
    endcase
3:    begin
    extin='.tiff'
    extout='.xdi'
    endcase
4:    begin
    extin='.xdi'
    extout='.edf'
    endcase
5:    begin
    extin='.dat'
    extout='.xdi'
    endcase
6:    begin
    extin='.xdi'
    extout='.tiff'
    endcase
7:    begin
    extin='.xdi'
    extout='.dat'
    endcase
8:    begin
    extin='.xml'
    extout='.xdi'
    endcase
endcase
WIDGET_CONTROL, ev.top, get_uvalue=list

t=CutPath(list.file,file=file)
filein=''
filein=MDIALOG_PICKFILE(path=list.path,file=file+extin,filter='*'+extin,$
                            title='Select files to convert...',/MULTIPLE_FILES)
if filein[0] EQ '' then return

for i=0l,n_elements(filein)-1 do begin
    tmp=CutPath(filein[i],path=pathout,file=fileout)
    fileout2=pathout+fileout+extout
    
    list.path=pathout
    list.file=fileout+'.xdi'
    WIDGET_CONTROL, ev.top, set_uvalue=list
    
    case type of
    0:  begin
        tmp=Cutpath(filein[i],path=path,file=file,ext=ext)
        if list.bopenformat then format=list.openformat
        XDI2IMS,path,file+ext,fileout2,list.OutID,format=format
        endcase
    1:  IMS2XDI,fileout2,filein[i],list.OutID
    2:  EDF2XDI,fileout2,filein[i],list.OutID,doreform=doreform
    3:    begin
        TIFF2XDI,fileout2,filein,list.OutID
        return ; All in 1 xdi file
        endcase
    4:  begin
        tmp=Cutpath(filein[i],path=path,file=file,ext=ext)
        if list.bopenformat then format=list.openformat
        XDI2EDF,path,file+ext,pathout+fileout,list.OutID,format=format
        endcase
    5:    PyMCADAT2XDI,fileout2,filein[i],list.OutID
    6:    begin
        tmp=Cutpath(filein[i],path=path,file=file,ext=ext)
        XDI2TIFF,path,file+ext,pathout+fileout,list.OutID
        endcase
    7:    begin
        tmp=Cutpath(filein[i],path=path,file=file,ext=ext)
        XDI2ATHENA,path,file+ext,pathout+fileout,list.OutID
        endcase
    8:    begin
        MICROXASFASTXRF2XDI,fileout2,filein[i],list.OutID
        endcase
    endcase
endfor
end;pro XDIconvert
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SinZingers,image,width=width,threshold=threshold

; Default kernel width for smoothing
if (keyword_set(width) eq 0) then width=9
width=fix(width)

; Default threshold for detecting zingers
if (keyword_set(threshold) eq 0) then threshold=1.2

; Compute ratio of unsmoothed image to smoothed image
smoothed=mSavGol(data=image,order=[0,0],nr=width/2,nl=width/2,r=3)
;ratio = float(image)/smoothed
ratio=(smoothed-image)/smoothed

; Find indices of pixels which are above threshold
;inger = where(ratio gt threshold, count)
zinger=where((ratio lt 1./threshold) or (ratio gt threshold),count)

; Replace pixels with smoothed value
if (count gt 0) then image[zinger] = smoothed[zinger]
return

; Replace pixels with zingers by the average of the pixels 2 to the left and
; 2 to the right.  We don't use nearest neighbor, since zingers sometimes hit 2
; adjacent pixels.
if (count gt 0) then image[zinger] = $
        (image[zinger-2] + image[zinger+2])/2.
return

; Replace pixels with zingers by the spline interpolated values
; using the none zinger pixels
if (count gt 0) then begin
s=size(image,/dimensions)
ncol=s[0]
nrow=s[1]
czinger=zinger mod ncol
rzinger=zinger / ncol
ind=indgen(ncol)

for i=0l,nrow-1 do begin
    indz=where(rzinger eq i,ct)
    if ct gt 0 then begin
        indz=czinger[indz]
        indn=ShrinkArray(ind,indz,count=ct)
        case ct of
        0:    image[*,i]=0
        1:    image[*,i]=image[indn,i]
        2:    image[*,i]=(image[indn[0],i]+image[indn[0],i])/2
        else: image[indz,i]=SPLINE(indn,image[indn,i],indz)
        endcase
    endif
endfor
endif

end;pro SinZingers
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro BacklashCorrect,image,shft

; Check dimensions
nshft=n_elements(shft)
s=size(image,/dimensions)
if nshft eq 1 then begin
    if shft eq 0 then return
endif else if nshft ne s[1] then return

ind=where(~finite(image),ct)
if ct ne 0 then image[ind]=0

if nshft eq 1 then begin
    ; Shift all even rows
    ind=indgen(s[0])
    indS=indgen(s[0])+shft
    for i=0l,s[1]-1,2 do image[*,i]=SPLINE(ind,image[*,i],indS)
    ;for i=0l,s[1]-1,2 do image[*,i]=interpol(image[*,i],ind,indS,/quad)
    ;for i=0l,s[1]-1,2 do image[*,i]=image[indS,i]
endif else begin
    for i=0l,s[1]-1 do begin
        ind=indgen(s[0])
        indS=indgen(s[0])+shft[i]
        tmp=SPLINE(ind,image[*,i],indS)

        shfti=ceil(shft[i])
        if shfti lt 0 then tmp[0:-shfti-1]=image[0,i] else $
        if shfti gt 0 then tmp[s[0]-shfti:*]=image[s[0]-1,i]
        
        image[*,i]=tmp
    endfor
endelse

end;pro BacklashCorrect
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SinNormal,image

s=size(image,/dimensions)
ncol=s[0]
nviews=s[1]

totar=total(image,1)
m=total(totar)/nviews
totar=reform(m/totar,1,nviews)
norm=replicate(1b,ncol)#totar
image=image*norm

end;pro SinNormal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_FindCOPManual,ID
end;CleanUp_FindCOPManual
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FindCOPManualPlot,list
wset,list.drawindex
LoadctRev,-39

plot,list.proj0

; rotate projection at 180 around the center
ind=2*(*list.pcen)-indgen(n_elements(list.proj180))
proj180=interpolate(list.proj180,ind,missing=0)

oplot,proj180*max(list.proj0)/max(proj180),col=100

plots,[*list.pcen,*list.pcen],!Y.crange
plots,[list.cen,list.cen],!Y.crange,linestyle=2

end;pro FindCOPManualPlot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FindCOPManual_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

type=widget_info(ev.id,/type)
widget_control,ev.top,get_uvalue=list

case type of
1:    widget_control,ev.top,/destroy
3:    begin
    widget_control,ev.id,get_value=cen
    *list.pcen=float(cen[0])
    FindCOPManualPlot,list
    
;    p=plot(list.proj0)
;    ind=2*(*list.pcen)-indgen(n_elements(list.proj180))
;    proj180=interpolate(list.proj180,ind,missing=0)
;    p2=plot(proj180*max(list.proj0)/max(proj180),/overplot,color=[255b,0,0])
    
    endcase
4:    begin
    cen=CONVERT_COORD(ev.x,ev.y,/device,/to_data)
    cen=cen[0]
    *list.pcen=cen
    FindCOPManualPlot,list
    ID=widget_info(ev.top,find_by_uname='cen')
    widget_control,ID,set_value=stringr(cen)
    endcase
endcase

end;pro FindCOPManual_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FindCOPManual,top,list
angles=SinAngles(list,nviews=nviews,bsinovert=bsinovert)*180/!pi

diff=(angles - (angles[0]+180)) mod 360
ind=where( abs(diff) lt 1e-3,ct)
if ct ne 1 then begin
    ind=(sort(abs(diff)))[0:1]
    if ~((diff[ind[0]] lt 0 and diff[ind[1]] gt 0) or (diff[ind[0]] gt 0 and diff[ind[1]] lt 0)) then return
    weight=abs(diff[ind])
    weight/=total(weight)
endif

image=reform((*list.dataorg)[list.NSelect,*,*],list.MapCol,list.MapRow)
if ~bsinovert then image=transpose(image)
s=size(image,/dimensions)
cen=(s[0]-1)/2.

pcen=ptr_new(list.Cor[1])
base=widget_base(title='',/column,GROUP_LEADER=top,/modal)

window,/pixmap
xsize=!d.x_size
ysize=!d.y_size
wdelete
draw=widget_draw(base,xsize=xsize,ysize=ysize,/BUTTON_EVENTS)

baset=widget_base(base,/row)
    text=widget_text(baset,value=stringr(*pcen),/editable,uname='cen')
    b=widget_button(baset,value='OK')

widget_control,base,/realize
widget_control,b,/INPUT_FOCUS
widget_control,draw,get_value=drawindex

case n_elements(ind) of 
1: proj180=image[*,ind[0]]
2: proj180=weight[0]*image[*,ind[0]]+weight[1]*image[*,ind[1]]
endcase

list2={proj0:image[*,0],proj180:proj180,pcen:pcen,cen:cen,drawindex:drawindex}
widget_control,base,set_uvalue=list2

FindCOPManualPlot,list2

Xmanager,'FindCOPManual',base, event_handler='FindCOPManual_Event',$
    cleanup='CleanUp_FindCOPManual',GROUP_LEADER=top

list.Cor[1]=*pcen
list.Corb[1]=1
ptr_free,pcen
end;pro FindCOPManual
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro WorkingDataCorrect,list,noTomo=noTomo

if list.ScanDim ne 2 then return

; Perform corrections:
;    1. Pixel juggling
;    2. Zinger removal
;    3. Backlash correction
;    4. Normalization

width=list.Cor[2]
threshold=list.Cor[3]
bTomo=widget_info(list.TomoReconID,/VALID_ID) eq 1

for i=0l,list.NMaps-1 do begin
    image=(*list.dataorg)[i,*,*]
    image=reform(image,list.MapCol,list.MapRow,/overwrite)

    ; Pixel juggling
    PixelJuggle,image,list,0,/meanpix,imageindex=i
    
    ; Manual backlash
    if list.ManualBacklash[0] ne 0 then begin
        BacklashCorrect,image,list.ManualBacklash[0] ; horizontal backlash
    endif
    if list.ManualBacklash[1] ne 0 then begin
        image=transpose(image)
        BacklashCorrect,image,list.ManualBacklash[1] ; vertical backlash
        image=transpose(image)
    endif
    
    if bTomo and ~keyword_set(noTomo) then begin
        ; Zinger removal
        if list.Corb[3] then SinZingers,image,width=width,threshold=threshold

        ; Correct sinogram for backlash
        if list.Corb[0] then BacklashCorrect,image,list.Cor[0]
        
        ; Elliptical wobble correction
        if list.Corb[4] then begin
            BacklashCorrect,image,wobbleshift(list)
        endif
        
        ; Manual wobble correction
        if list.Corb[5] then begin
            ;shft=[21.7940,27.1951]-list.Cor[1]
            ;shft=[make_array(61,value=shft[0]),make_array(60,value=shft[1])]
            ;BacklashCorrect,image,shft
            BacklashCorrect,image,*list.WobblePtr
        endif
        
        ; Normalize sinogram
        if list.Corb[2] then SinNormal,image
    endif

    image=reform(image,1,list.MapCol,list.MapRow,/overwrite)
    (*list.data)[i,*,*]=image
endfor

if list.Corb[3] then begin
    list.Cor[2]=width
    list.Cor[3]=threshold

    ID=widget_info(list.top,FIND_BY_UNAME='Zwidth')
        widget_control,ID,set_value=stringr(list.Cor[2])
    ID=widget_info(list.top,FIND_BY_UNAME='Zthreshold')
        widget_control,ID,set_value=stringr(list.Cor[3])
endif

end;pro WorkingDataCorrect
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ModelSinogram,list,rho=rho,nviews=nviews,bsinovert=bsinovert

; No backlash and center correction or wobble correction
bfit=list.Corb[1] or list.Corb[4] or list.Corb[5]
if (list.Corb[0] eq 0) and (bfit eq 0) then return

; Use correct sinogram, but without typical sinogram corrections
WorkingDataCorrect,list,/noTomo
image=reform((*list.data)[list.NSelect,*,*],list.MapCol,list.MapRow)
angles=SinAngles(list,nviews=nviews,bsinovert=bsinovert)
if ~bsinovert then image=transpose(image)

; Zinger removal
if list.Corb[3] then SinZingers,image,width=width,threshold=threshold

; Define some constants
s=size(image,/dimensions)
nrow=s[1]
ncol=s[0]
;if (ncol mod 2) eq 0 then begin ; make number of columns uneven
;    ncol=ncol-1
;    image=image[0:ncol-1,*]
;endif

indeven=2*indgen(nrow/2 + (nrow mod 2))
indodd=2*indgen(nrow/2)+1
weights=replicate(1b,nrow)

; Calculate center of gravity for each row
w = indgen(ncol)+1
CG = reform(image ## w)/total(image,1)

; Reject outliers
ind=where(abs(CG-smooth(CG,2,/nan)) gt mean(abs(CG[1:*]-CG[0:nrow-2])) or ~finite(CG),ct)
if ct ne 0 then CG[ind]=!Values.F_NAN

; Initial sinus model:
; 
;  Tomogram coordinate system x-y
;  Sinogram coordinate system X-Y   (i.e. angles-rho)
; 
;  Point at the center of rotation
;   x = Ax
;   y = Ay
; 
;  Rotation of a point with polar coord. [B,C] around center of rotation
;   x = Ax + B.cos(X+C)
;   y = Ay + B.sin(X+C)
; 
;  Rotation of a point around a wobbling center of rotation
;   x = Ax + B.cos(X+C) + R(X).cos(X+C)
;   y = Ay + B.sin(X+C) + R(X).sin(X+C)
;   
;   Radial distance of ellipse:
;   R(X) = D.E/sqrt(D^2.sin(X-F)^2+E^2.cos(X-F)^2)
; 
;  Projection on the X axis:
;  Y=A+(B+R(X)).cos(X+C)
;    A: center of projection
;    B: amplitude
;    C: phase
;    D: first axis of wobble-ellipse <-- list.Corb[4]
;    E: second axis of wobble-ellipse <-- list.Corb[4]
;    F: wobble-ellipse tilt with respect to the tomogram coord. system <-- list.Corb[4]
;        F=0 => ellipse alligned to tomogram coordinate system with D along x
;    
; TODO: check whether wobble is valid
bfix=[0b,0b,0b,~list.Corb[4],~list.Corb[4],~list.Corb[4]]
A=[(ncol+1)/2.,(max(CG)-min(CG))/2,0.]
if list.Corb[4] then begin
    A=[A,list.Cor[4:6]*0+1]
    struc={bfix:bfix}
endif else begin
    struc={bfix:bfix,Afix:[0.,0,0]}
endelse

; Model even and odd separate and shift even with the difference
if list.Corb[0] then begin
    ; Model sinus to CG's for odd rows
    CGo = nlls(angles[indodd], CG[indodd], A,'gfuncSin', struc=struc,ystdev=weights[indodd],CHISQ=CHISQ,SIGMA=SIGMA,iter=iter,itmax=50)
    Ceno = A[0]-1

    ; Plot results
    window,0
    plot,angles[indodd]*180/!pi, CG[indodd],psym=2,xtitle='Rotation angle (deg)',ytitle='Center-of-gravity (odd rows)',/xs
    oplot,angles[indodd]*180/!pi,CGo
    printw,list.OutID,'Red. Chisq + niter: '+string(CHISQ)+string(ITER)
    printw,list.OutID,'Center of projection for odd: '+string(Ceno)+string(SIGMA[0])

    ; Model sinus to CG's for even rows
    CGe = nlls(angles[indeven], CG[indeven], A,'gfuncSin', struc=struc,ystdev=weights[indeven],CHISQ=CHISQ,SIGMA=SIGMA,iter=iter,itmax=50)
    Cene = A[0]-1

    ; Plot results
    window,1
    plot,angles[indeven]*180/!pi, CG[indeven],psym=2,xtitle='Rotation angle (deg)',ytitle='Center-of-gravity (even rows)',/xs
    oplot,angles[indeven]*180/!pi,CGe
    printw,list.OutID,'Red. Chisq + niter: '+string(CHISQ)+string(ITER)
    printw,list.OutID,'Center of projection for even: '+string(Cene)+string(SIGMA[0])

    ; Calculate and perform backlash correction
    Backlashshift = (Cene - Ceno)
    BacklashCorrect,image,Backlashshift

    ; Calculate new centers of gravity for calculation of the center of projection
    if bfit then begin
        ; Normalize sinogram
        if list.Corb[2] then SinNormal,image

        CG = reform(image ## w)/total(image,1)
    endif

    list.Cor[0]=Backlashshift
endif else list.Cor[0]=0.

; Model sinus to CG's for all rows to get the projection center
if bfit then begin

    CGn = nlls(angles, CG, A,'gfuncSin', struc=struc,ystdev=weights,CHISQ=CHISQ,SIGMA=SIGMA,iter=iter,errormsg=errormsg)
    ;A[3:5]=[50,40,20/180.*!pi]
    ;gfuncSin,angles,A,CGn, struc=struc
    
    list.Cor[1]=A[0]-1 ; Center of projection
    if list.Corb[4] then list.Cor[4:6]=A[3:5] ; Elliptical wobble
    if list.Corb[5] then begin ; Row-shift wobble
        shft=CG-CGn
        ind=where(~finite(shft),ct)
        if ct ne 0 then shft[ind]=0
        *list.WobblePtr=shft
        CG=CGn
    endif
    rho=CGn-1

    ; Plot results
    window,2
    plot,angles*180/!pi, CG,psym=2,xtitle='Rotation angle (deg)',ytitle='Center-of-gravity (all rows)',/xs
    oplot,angles*180/!pi,CGn,color=100
    printw,list.OutID,'Red. Chisq + niter: '+string(CHISQ)+string(ITER)
    printw,list.OutID,'Center of projection: '+string(list.Cor[1])+string(SIGMA[0])

endif else list.Cor[1]=0

end;pro ModelSinogram
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro printfphysmask,lun,list,NSelect

printf,lun,'$ORIENTATION:'
printf,lun,[list.TRAN,list.REV1,list.REV2,list.REV3]
printf,lun,list.Hsign,list.Vsign
printf,lun,list.ScanDim

NMaps=n_elements(NSelect)

printf,lun,'$PHYSMASK:'
printf,lun,NMaps
printf,lun,'Qualitative analysis:'

ntype3=0l
datatype3=['Group ',(*list.pmfields)[*,0]]
ind=where(datatype3 ne '-',ct3)
if ct3 ne 0 then datatype3=datatype3[ind]
ntype3tot=0l
indtype3=-1l

ntype10=0l
datatype10=(*list.pmfields)[*,1]
ind=where(datatype10 ne '-',ct10)
if ct10 ne 0 then datatype10=datatype10[ind]
ntype10tot=0l
indtype10=-1l

datacombi=datatype3
ncombi=0l
indcombi=-1l

for j=0l,NMaps-1 do begin
    i=NSelect[j]+1
    case (*list.pmstruc).(i).combi of
    0:    begin ; groups as loaded
        case (*list.pmstruc).(i).labeli of
        0:    begin
            datatype3=[[datatype3],[replicate(string(ntype3,format='(I5)'),1,(*list.pmstruc).(i).n),(*list.pmstruc).(i).data[0:ct3-2,*]]]  ; ct3-2 because of 'Group'
            ntype3++
            ntype3tot+=(*list.pmstruc).(i).n
            indtype3=[indtype3,j]
            endcase
        1:    begin
            datatype10=[[datatype10],[(*list.pmstruc).(i).data[0:ct10-1,*]]]
            ntype10++
            ntype10tot+=(*list.pmstruc).(i).n
            indtype10=[indtype10,j]
            endcase
        endcase
        endcase
    1:    begin ; groups as combined
        datacombi=[[datacombi],[replicate(string(ncombi,format='(I5)'),1,(*list.pmstruc).(i).n),(*list.pmstruc).(i).data[0:ct3-2,*]]]  ; ct3-2 because of 'Group'
        ncombi++
        indcombi=[indcombi,j]
        endcase
    2:    begin    ; groups as overlapped
        datacombi=[[datacombi],[replicate(string(ncombi,format='(I5)')+'*',1,(*list.pmstruc).(i).n),(*list.pmstruc).(i).data[0:ct3-2,*]]]  ; ct3-2 because of 'Group'
        ncombi++
        indcombi=[indcombi,j]
        endcase
    endcase
endfor

indtype=[indtype3,indtype10,indcombi]
indtype=[0,indtype[where(indtype ne -1)]+1]

printf,lun,[0,0,0,ntype3]
printf,lun,datatype3[*,0]
if ntype3tot ne 0 then printf,lun,datatype3[*,1:*]
printf,lun,'Quantitative analysis:'
printf,lun,[ntype10,0,0]
printf,lun,ntype10tot
printf,lun,datatype10[*,0]
if ntype10tot ne 0 then printf,lun,datatype10[*,1:*]
if ncombi ne 0 then begin
    printf,lun,'Combined:'
    printf,lun,datacombi[*,1:*]
endif
printf,lun,'$Reorder:'
printf,lun,indtype

end;pro printfphysmask
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro UpdateXDITD,lun,list

printf,lun,'$HEADING:'
coord=getxdicoord(list,(list.coord)[*,0])
unit=list.unit[0]
case unit of
'mm': begin
    mul=1000.
    unit=msymbols('micro')+'m'
    endcase
'deg': begin
    mul=1
    endcase
endcase

fi='I'+stringr(max(strlen(stringr([list.MapCol>list.MapRow,list.expo])))+1)
if list.TRAN then begin
    printf,lun,format='("dscan samv ",F10.3,F10.3,'+fi+','+fi+')',$
       coord,list.MapRow-1,list.expo
    printf,lun,format='("size: ", '+fi+'," pixels")',list.MapRow
    dif=((coord[0]>coord[1])-(coord[0]<coord[1]))*mul
    printf,lun,format='("vertical: ",F10.3," '+unit+' in steps of ",F10.3," '+unit+'")',dif,dif/(list.MapRow-1)
endif else begin
    printf,lun,format='("dscan samh ",F10.3,F10.3,'+fi+','+fi+')',$
       coord,list.MapCol-1,list.expo
    printf,lun,format='("size: ", '+fi+'," pixels")',list.MapCol
    dif=((coord[0]>coord[1])-(coord[0]<coord[1]))*mul
    printf,lun,format='("horizontal: ",F10.3," '+unit+' in steps of ",F10.3," '+unit+'")',dif,dif/(list.MapCol-1)
endelse

end;pro UpdateXDITD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro UpdateXDIDD,lun,list

printf,lun,'$HEADING:'
if list.REV3 then mesh='meshc' else mesh='mesh'
coord=getxdicoord(list,list.coord)
unit=list.unit
if list.TRAN then begin
    coord=reverse(coord,2)
    unit=reverse(unit)
endif

mul=fltarr(2)
for i=0l,1 do begin
    case unit[i] of
    '': begin
        unit[i]='mm'
        mul[i]=1
        endcase
    msymbols('micro')+'m': begin
        mul[i]=1
        endcase
    'mm': begin
        mul[i]=1000.
        unit[i]=msymbols('micro')+'m'
        endcase
    'deg': begin
        mul[i]=1
        endcase
    endcase
endfor

fi='I'+stringr(max(strlen(stringr([list.MapCol,list.MapRow,list.expo])))+1)
if list.TRAN then begin
    printf,lun,format='("'+mesh+' samv ",F10.3,F10.3,'+fi+'," samh ",F10.3,F10.3,'+fi+','+fi+')',$
       coord[*,0],list.MapRow-1,coord[*,1],list.MapCol-1,list.expo
    printf,lun,format='("size: ", '+fi+', " x ", '+fi+', " pixels")',list.MapCol,list.MapRow
    dif=((coord[0,0]>coord[1,0])-(coord[0,0]<coord[1,0]))*mul[0]
    printf,lun,format='("vertical: ",F10.3," '+unit[0]+' in steps of ",F10.3," '+unit[0]+'")',dif,dif/(list.MapRow-1)
    dif=((coord[0,1]>coord[1,1])-(coord[0,1]<coord[1,1]))*mul[1]
    printf,lun,format='("horizontal: ",F10.3," '+unit[1]+' in steps of ",F10.3," '+unit[1]+'")',dif,dif/(list.MapCol-1)
endif else begin
    printf,lun,format='("'+mesh+' samh ",F10.3,F10.3,'+fi+'," samv ",F10.3,F10.3,'+fi+','+fi+')',$
       coord[*,0],list.MapCol-1,coord[*,1],list.MapRow-1,list.expo
    printf,lun,format='("size: ", '+fi+', " x ", '+fi+', " pixels")',list.MapCol,list.MapRow
    dif=((coord[0,0]>coord[1,0])-(coord[0,0]<coord[1,0]))*mul[0]
    printf,lun,format='("horizontal: ",F10.3," '+unit[0]+' in steps of ",F10.3," '+unit[0]+'")',dif,dif/(list.MapCol-1)
    dif=((coord[0,1]>coord[1,1])-(coord[0,1]<coord[1,1]))*mul[1]
    printf,lun,format='("vertical: ",F10.3," '+unit[1]+' in steps of ",F10.3," '+unit[1]+'")',dif,dif/(list.MapRow-1)
endelse

end;pro UpdateXDIDD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro UpdateXDI,ev,tomo=tomo
;TODO: use SaveXDI from xrd_bp
widget_control,ev.top,get_uvalue=list

; ----Check for tomography data----
if keyword_set(tomo) then $
if (list.ScanDim ne 2) or (widget_info(list.TomoReconID,/VALID_ID) ne 1) then begin
    Result = DIALOG_MESSAGE('No tomography data found.',/information)
    return
endif

; ----Select file and open----
xdipath=SelFile(list.path,list.file,'*.xdi','Save To....')
if xdipath eq '' then return

; ----Select groups----
prom='0-'+stringr(list.NMaps-1)
NSelect=PromptNumber(prom,ev.top,'Select groups:')
if ParseIndices(NSelect) then return
NMaps=n_elements(NSelect)

; ----Tomography data----
if keyword_set(tomo) then begin
    ; Data block
    tomoinfo,list,Nx=Nx,Ny=Ny,coordtrans=coordtrans,transunit=transunit
    
    block=fltarr(NMaps,Nx,Ny)
    for j=0l,NMaps-1 do begin
        i=NSelect[j]
        tmp=Makewplot3DImage(list,i,tomogram=tomogram)
        block[j,*,*]=temporary(tomogram)
    endfor
    list.MapCol=Nx
    list.MapRow=Ny
    
    ; Change coordinates
    list.coord=coordtrans
    list.unit=transunit
endif

; ----Update----
if openw_safe(lun,xdipath) then return
case list.ScanDim of
1:    UpdateXDITD,lun,list
2:    UpdateXDIDD,lun,list
endcase

printfphysmask,lun,list,NSelect

printf,lun,'$IMS:'
printf,lun,3
printf,lun,NMaps,list.MapCol,list.MapRow
if keyword_set(tomo) then begin
    writeu,lun,block
endif else begin
    Result = DIALOG_MESSAGE('Save with corrections?',/question)
    if Result eq 'Yes' then pdata=list.data else pdata=list.dataorg
    writeu,lun,(*pdata)[NSelect,*,*]
endelse

if PTR_VALID(list.filenames) and ~keyword_set(tomo) then begin
    printf,lun,''
    printf,lun,'$FILES:'
    ; Default orientation
    block=*list.filenames
    ReformXDIDatablock,block,[list.TRAN,0b],[list.REV1,0b],[list.REV2,0b],[list.REV3,0b],0b
    printf,lun,block,format='(A)'
endif

free_lun,lun

end;pro UpdateXDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro xdi_setdim,list,limit=limit
XSize0=list.XSize
YSize0=list.YSize

if keyword_set(limit) then begin
    ; Make sure it fits on screen, keeps the aspect ratio and 
    ; uses an integer magnification if possible
    device,get_screen_size=screen
    screen=round(screen*0.8)
    
    list.XSize=1>round(list.MapCol/list.sfh)<screen[0]
    list.YSize=1>round(list.MapRow/list.sfv)<screen[1]

    mag=float(list.XSize)/list.MapCol < float(list.YSize)/list.MapRow
    if mag gt 1 then mag=floor(mag)
    list.sfh=1./mag
    list.sfv=1./mag
endif

list.XSize=1>round(list.MapCol/list.sfh)
list.YSize=1>round(list.MapRow/list.sfv)
list.sfh=list.MapCol/float(list.XSize)
list.sfv=list.MapRow/float(list.YSize)

end;pro xdi_setdim
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro xdi_changesize,list
; Don't forget to save list afterwards

wdelete,list.pixindex
if widget_info(list.TomoReconID,/VALID_ID) eq 1 then widget_control,list.TomoReconID,/destroy

case list.ScanDim of
1:    begin
    ; No tomo
    list.TomoReconID=0L
    list.TomoReconIndex=0L
    
    ; Pixmap
    Window, /Free, /Pixmap
    list.pixindex = !D.Window
    list.XSize=!d.x_size
    list.YSize=!d.y_size
    
    ; Draw widget
    widget_control,list.draw,/destroy
    ID=widget_info(list.top,FIND_BY_UNAME='drawadd')
    device,get_screen_size=screen
    list.draw = WIDGET_DRAW(ID,uname='draw',uvalue=0,xsize=list.XSize,ysize=list.YSize)
    widget_control,list.draw,get_value=drawindex,draw_motion_event=1,event_pro='DEFAULT_PROCESS_EVENTS'
    list.drawindex=drawindex
    endcase
2:    begin
    device,get_screen_size=screen
    screen[0]*=0.3
    screen[1]*=0.6
    
    ; Dimension from multipliers
    xdi_setdim,list
    
    ; Pixmap
    Window, /Free, xsize=list.XSize, ysize=list.YSize, /Pixmap
    list.pixindex = !D.Window
    
    ; Draw widget
    widget_control,list.draw,/destroy
    IDparent=widget_info(list.top,FIND_BY_UNAME='drawadd')
    list.draw = WIDGET_DRAW(IDparent,uname='draw',uvalue=0,xsize=list.XSize,ysize=list.YSize,$
                            /scroll,x_scroll_size=screen[0]<list.XSize, y_scroll_size=screen[1]<list.YSize)
    widget_control,list.draw,get_value=drawindex,draw_motion_event=1,event_pro='DEFAULT_PROCESS_EVENTS'
    list.drawindex=drawindex
    
    ID=widget_info(list.top,FIND_BY_UNAME='Tomo Reconstruct')
    widget_control,ID,get_uvalue=uval
    if uval ne 0 then begin
        xsize=TomoSize(list,/X)
        ysize=TomoSize(list,/Y)
          list.TomoReconID = WIDGET_DRAW(IDparent, xsize=xsize, ysize=ysize,/scroll,x_scroll_size=screen[0]<xsize, y_scroll_size=screen[1]<ysize)
          widget_control,list.TomoReconID,get_value=TomoReconIndex,draw_motion_event=1,event_pro='DEFAULTTOMO_PROCESS_EVENTS'
          list.TomoReconIndex=TomoReconIndex
    endif
        
    endcase
endcase

end;pro xdi_changesize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ChangeXDIcoord_Event,ev
WIDGET_CONTROL, ev.top, get_uvalue=top
WIDGET_CONTROL, top, get_uvalue=list

widget_control,ev.id,get_uvalue=uval
case uval of
'c00':    begin
        widget_control,ev.id,get_value=val
        list.coord[0,0]=float(val)
        endcase
'c01':    begin
        widget_control,ev.id,get_value=val
        list.coord[0,1]=float(val)
        endcase
'c10':    begin
        widget_control,ev.id,get_value=val
        list.coord[1,0]=float(val)
        endcase
'c11':    begin
        widget_control,ev.id,get_value=val
        list.coord[1,1]=float(val)
        endcase
'Unit0':begin
        list.unit[0]=(['mm','deg'])[ev.index]
        endcase
'Unit1':begin
        list.unit[1]=(['mm','deg'])[ev.index]
        endcase
'Swap':    begin
        list.coord=list.coord[*,[1,0]]
        list.unit=list.unit[[1,0]]
        endcase
endcase

list.coord[*,0]=list.coord[sort(list.coord[*,0]),0]
list.coord[*,1]=list.coord[sort(list.coord[*,1]),1]

id=widget_info(ev.top,find_by_uname='c00')
widget_control,id,set_value=stringr(list.coord[0,0])
id=widget_info(ev.top,find_by_uname='c10')
widget_control,id,set_value=stringr(list.coord[1,0])
id=widget_info(ev.top,find_by_uname='c01')
widget_control,id,set_value=stringr(list.coord[0,1])
id=widget_info(ev.top,find_by_uname='c11')
widget_control,id,set_value=stringr(list.coord[1,1])
id=widget_info(ev.top,find_by_uname='Unit0')
widget_control,id,set_droplist_select=list.unit[0] eq 'deg'
id=widget_info(ev.top,find_by_uname='Unit1')
widget_control,id,set_droplist_select=list.unit[1] eq 'deg'

xdi_changesize,list
WIDGET_CONTROL, top, set_uvalue=list
RefreshDisplayXDI,list
end;pro ChangeXDIcoord_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ChangeXDIcoord,ev
WIDGET_CONTROL, ev.top, get_uvalue=list

base=widget_base(title='Change coordinates',/column,uvalue=ev.top,GROUP_LEADER=ev.top)

xs=15
base1=widget_base(base,/row)
button=widget_label(base1,value='Horizontal motor:',xsize=xs*5)
text=widget_text(base1,value=stringr(list.coord[0,0]),xsize=xs,/editable,uvalue='c00',uname='c00')
text=widget_text(base1,value=stringr(list.coord[1,0]),xsize=xs,/editable,uvalue='c10',uname='c10')
drop1=widget_droplist(base1,value=['mm','deg'],uvalue='Unit0',uname='Unit0')

base1=widget_base(base,/row,sensitive=list.ScanDim eq 2)
button=widget_label(base1,value='Vertical motor:',xsize=xs*5)
text=widget_text(base1,value=stringr(list.coord[0,1]),xsize=xs,/editable,uvalue='c01',uname='c01')
text=widget_text(base1,value=stringr(list.coord[1,1]),xsize=xs,/editable,uvalue='c11',uname='c11')
drop2=widget_droplist(base1,value=['mm','deg'],uvalue='Unit1',uname='Unit1')

button=widget_button(base,value='Swap',uvalue='Swap')

widget_control,base,/realize
widget_control,drop1,SET_DROPLIST_SELECT=list.unit[0] eq 'deg'
widget_control,drop2,SET_DROPLIST_SELECT=list.unit[1] eq 'deg'

Xmanager,'ChangeXDIcoord',base, event_handler='ChangeXDIcoord_Event',GROUP_LEADER=ev.top
end;pro ChangeXDIcoord
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EnlargeDataBlock_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

type=widget_info(ev.id,/type)
case type of
4:
else:    begin
    WIDGET_CONTROL, ev.top, get_uvalue=list
    
    ; Clear display
    WSet, list.drawindex
    h=!D.x_size
    v=!D.y_size
    Device, Copy = [0,0, h,v, 0,0,list.pixindex]
        
    ; Get new size
    dim=[list.NMaps,list.MapCol,list.MapRow]
    off=[0,0,0]

    id=find_by_uname(ev.top,'col')
    deltadim1=ev.id eq id
    if widget_info(id,/valid) then begin
        widget_control,id,get_value=val
        dim[1]=long(val)
    endif
    id=find_by_uname(ev.top,'row')
    if widget_info(id,/valid) then begin
        widget_control,id,get_value=val
        dim[2]=long(val)
    endif
    id=find_by_uname(ev.top,'coloff')
    if widget_info(id,/valid) then begin
        widget_control,id,get_value=val
        off[1]=long(val)
    endif
    id=find_by_uname(ev.top,'rowoff')
    if widget_info(id,/valid) then begin
        widget_control,id,get_value=val
        off[2]=long(val)
    endif
    
    id=find_by_uname(ev.top,'type')
    bresize=0
    if widget_info(id,/valid) then $
        widget_control,id,get_uvalue=bresize
    
    ; check new size
    r=0
    if bresize then begin
        n=list.MapCol*list.MapRow
        if deltadim1 then dim[2]=n/dim[1] else dim[1]=n/dim[2]
        r=dim[1]*dim[2] ne list.MapCol*list.MapRow
    endif
    r or= dim[1] le 0 or dim[2] le 0
    if r then begin
        dim[1] = list.MapCol
        dim[2] = list.MapRow
        r=dialog_message("Can't resize to these dimensions.")
    endif
    id=find_by_uname(ev.top,'col')
    if widget_info(id,/valid) then widget_control,id,set_value=stringr(dim[1])
    id=find_by_uname(ev.top,'row')
    if widget_info(id,/valid) then widget_control,id,set_value=stringr(dim[2])
    
    case type of
    1:    begin
        modeXDI, {id:0l,top:ev.top}, 'Dummy'

        ; Coordinates
        s1=list.MapCol
        s2=list.MapRow
        i1=0
        i2=1
        if list.TRAN then begin
            i1=1
            i2=0
        endif
        if s1 ne 1 then begin
            inc=(list.coord[1,i1]-list.coord[0,i1])/(s1-1)
            list.coord[*,i1]=list.coord[0,i1]+inc*[off[1],off[1]+dim[1]-1]
        endif
        if s2 ne 1 then begin
            inc=(list.coord[1,i2]-list.coord[0,i2])/(s2-1)
            list.coord[*,i2]=list.coord[0,i2]+inc*[off[2],off[2]+dim[2]-1]
        endif
        
        ; Resize
        list.MapCol=dim[1]
        list.MapRow=dim[2]
        if dim[1] eq 1 or dim[2] eq 1 then begin
            ScanDim=(DIALOG_MESSAGE("Two-dimensional map?",/question) eq "Yes")+1
            if ScanDim eq 1 then list.TRAN=list.MapCol eq 1
        endif else ScanDim=2
        ;if ScanDim ne list.ScanDim then begin
        ;TODO: do something with list.coord
        ;endif
        list.ScanDim=ScanDim
        xdi_changesize,list
        WIDGET_CONTROL, list.top, set_uvalue=list
        
        ; Resize data and filenames
        case bresize of
        0:    begin
            (*list.dataorg)=extrac(*list.dataorg,off[0],off[1],off[2],dim[0],dim[1],dim[2])
        
            if PTR_VALID(list.filenames) then $
                (*list.filenames)=extrac(*list.filenames,0,off[1],off[2],1,dim[1],dim[2])
            endcase
        1:    begin
            (*list.dataorg)=reform(*list.dataorg,dim[0],dim[1],dim[2],/overwrite)
        
            if PTR_VALID(list.filenames) then $
                (*list.filenames)=reform(*list.filenames,1,dim[1],dim[2],/overwrite)
            endcase
        endcase
        (*list.data)=(*list.dataorg)
        
        ; Finish
        WorkingDataCorrect,list
        WIDGET_CONTROL, list.top, set_uvalue=list
        RefreshDisplayXDI,list
        if bresize eq 1 then RefreshWindowXDI,{top:list.top}
        
        ; Destroy dialog
        widget_control,ev.top,/destroy
        endcase
    3:    begin
        if bresize ne 0 then return
        
        WIDGET_CONTROL, ev.top, get_uvalue=list
        
        ev=[off,off+dim-1]
        case list.ScanDim of
        1:    begin
            if list.TRAN then begin
                tmp=(!X.crange[1]-!X.crange[0])/2.
                PlotS, [tmp,tmp],ev[[2,5]],Linestyle=0,/data
            endif else begin
                tmp=(!Y.crange[1]-!Y.crange[0])/2.
                PlotS, ev[[1,4]],[tmp,tmp],Linestyle=0,/data
            endelse
            endcase
        2:    begin
            ev[[1,4]]=RTStat(ev[[1,4]],list.sfh)
            ev[[2,5]]=!d.Y_SIZE-1-RTStat(ev[[2,5]],list.sfv)
            PlotS, ev[[1,1]],ev[[2,5]],Linestyle=0,/device
            PlotS, ev[[4,4]],ev[[2,5]],Linestyle=0,/device
            PlotS, ev[[1,4]],ev[[2,2]],Linestyle=0,/device
            PlotS, ev[[1,4]],ev[[5,5]],Linestyle=0,/device
            endcase
        endcase
        
        endcase
    endcase
    
    endelse
endcase
        
end;pro EnlargeDataBlock_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EnlargeDataBlock,ev,type

modeXDI, {id:0l,top:ev.top}, 'Dummy'
WIDGET_CONTROL, ev.top, get_uvalue=list

; Change dimensions
base=widget_base(title='',/column,uvalue=list,GROUP_LEADER=ev.top)

baset=widget_base(base,/row)
    b=widget_label(baset,value='Column offset: ')
    b=widget_text(baset,value=stringr(0),uname='coloff',editable=type eq 0)
    b=widget_label(baset,value='Columns: ')
    b1=widget_text(baset,value=stringr(list.MapCol),uname='col',/editable)

baset=widget_base(base,/row)
    b=widget_label(baset,value='Row offset: ')
    b=widget_text(baset,value=stringr(0),uname='rowoff',editable=type eq 0)
    b=widget_label(baset,value='Rows: ')
    b=widget_text(baset,value=stringr(list.MapRow),uname='row',/editable)

b=widget_button(base,value='OK',uname='type',uvalue=type)
widget_control,base,/realize
widget_control,b,/INPUT_FOCUS

XDIdrawfree,list.draw
widget_control,list.draw,draw_motion_event=1,$
              draw_button_event=1,event_pro='LOF_PROCESS_EVENTS',set_uvalue={type:3,top:base}

widget_control,b1,send_event={ID:b1,TOP:base,HANDLER:base}
Xmanager,'EnlargeDataBlock',base, event_handler='EnlargeDataBlock_Event',GROUP_LEADER=ev.top

end;pro EnlargeDataBlock
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ConvertDataBlock,ev,chcode
WIDGET_CONTROL, ev.top, get_uvalue=list

; chcode:
;    first 5 bits => which variable changed (bit0->deltaREV2,bit1->deltaREV1,bit2->deltaTRAN,bit3->deltaREV3,bit4->SWAP)
;    bit6 => new value

new=(chcode and 32) eq 32
btran=((chcode and 4) eq 4)
b1=((chcode and 2) eq 2)
b2=((chcode and 1) eq 1)
b3=((chcode and 8) eq 8)
b4=((chcode and 16) eq 16)
TRAN=[list.TRAN,btran?new:list.TRAN]
REV1=[list.REV1,b1?new:list.REV1]
REV2=[list.REV2,b2?new:list.REV2]
REV3=[list.REV3,b3?new:list.REV3]
SWAP=b4?new:0

block=(*list.dataorg)
ReformXDIDatablock,block,TRAN,REV1,REV2,REV3,SWAP,s0,s1,s2
(*list.data)=block
(*list.dataorg)=block
list.MapCol=s1
list.MapRow=s2

if ptr_valid(list.filenames) then begin
    block=*list.filenames
    ReformXDIDatablock,block,TRAN,REV1,REV2,REV3,SWAP
    *list.filenames=block
endif

; Change others
case list.ScanDim of
1:
2:    begin
    if btran or SWAP then begin
        xdi_changesize,list
        if list.pverts ne -1 then begin
            t=(*list.xverts)
            (*list.xverts)=(*list.yverts)
            (*list.yverts)=t
        endif
    endif
    endcase
endcase

list.TRAN=TRAN[1]
list.REV1=REV1[1]
list.REV2=REV2[1]
list.REV3=REV3[1]
list.SWAP=~list.SWAP

WorkingDataCorrect,list
WIDGET_CONTROL, ev.top, set_uvalue=list
RefreshDisplayXDI,list

end;pro ConvertDataBlock
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CursorPropTomo,ID,X,Y,list,XX=XX,YY=YY,nr=nr

tomoinfo,list,Nx=Nx,Ny=Ny,coordtrans=coordtrans,transunit=transunit
nr=CurPos(X,Y,Nx,Ny,list.sfh,list.sfh,XX=XX,YY=YY)

CATCH, Error_status
IF Error_status NE 0 THEN return

pX=coordtrans[0,0]+(coordtrans[1,0]-coordtrans[0,0])/(Nx-1)*XX
pY=coordtrans[0,1]+(coordtrans[1,1]-coordtrans[0,1])/(Ny-1)*YY
pX=stringr(pX)+transunit
pY=stringr(pY)+transunit

if widget_info(ID,/valid) then begin
    text=strarr(3)
    text[0]='POSITION=['+stringr(XX)+','+stringr(YY)+']'+', COORD.=['+pX+','+pY+']'
    text[1]='INT.='+stringr((*list.ptrtomocurrent)[XX,YY])
    widget_control,ID,set_value=text
endif
end;pro CursorPropTomo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CursorProp,ID,X,Y,list,XX=XX,YY=YY,nr=nr,_nrX=nrX,_nrY=nrY

text=strarr(3)
case list.ScanDim of
2:    begin
    nr=CurPos(X,Y,list.MapCol,list.MapRow,list.sfh,list.sfv,XX=XX,YY=YY)
    nr=ConvertFileNumber(XX,YY,list.MapCol,list.MapRow,list.TRAN,list.REV1,list.REV2,list.REV3,XX=nrX,YY=nrY)
    coord=getxdicoord(list,list.coord)
    unit=list.unit
    if list.TRAN then begin
        coord=reverse(coord,2)
        unit=reverse(unit)
    endif
    
    ; First motor
    pX=coord[0,0]+(coord[1,0]-coord[0,0])/((list.TRAN?list.MapRow:list.MapCol)-1)*nrX
    
    ; Second motor
    pY=coord[0,1]+(coord[1,1]-coord[0,1])/((list.TRAN?list.MapCol:list.MapRow)-1)*nrY
    
    pX=stringr(pX)+unit[0]
    pY=stringr(pY)+unit[1]
        
    if list.TRAN then coordstr=pY+','+pX $
    else coordstr=pX+','+pY
    
    if widget_info(ID,/valid) then begin
        text[0]='POSITION=['+stringr(XX)+','+stringr(YY)+'], COORD.=['+coordstr+']'
        text[1]='FILEPOS.['+stringr(nrX)+','+stringr(nrY)+'], FILENR.='+stringr(nr)
        text[2]=['INT.='+stringr((*list.data)[list.NSelect,XX,YY])]
    endif
    
    endcase
1:    begin
    wset,list.drawindex
    SetSysVar,list.sysvar
    coord=CONVERT_COORD(X,Y,/device,/to_data)
    coordi=coord
    n=list.MapCol*list.MapRow
    
    xrange=wplot2Dgetxrange(list,n)
    xrangei=xrange
    
    if (list.TRAN and list.REV2) or (~list.TRAN and list.REV1) then xrange=reverse(xrange)
    bin=(xrange[1]-xrange[0])/(n-1.)
    if bin ne 0 then coord=(coord-xrange[0])/bin
    
    bin=(xrangei[1]-xrangei[0])/(n-1.)
    if bin ne 0 then coordi=(coordi-xrangei[0])/bin
    
    if list.TRAN then begin
        n=list.MapRow
        nr=0>round(coordi[1,*])<(n-1)
        YY=nr
        XX=lonarr(n_elements(nr))
        
        nr=0>round(coord[1,*])<(n-1)
        nrX=0
        nrY=nr
    endif else begin
        n=list.MapCol
        nr=0>round(coordi[0,*])<(n-1)
        XX=nr
        YY=lonarr(n_elements(nr))
        
        nr=0>round(coord[0,*])<(n-1)
        nrY=0
        nrX=nr
    endelse
    
    coord=getxdicoord(list,list.coord)
    p=coord[0,0]+(coord[1,0]-coord[0,0])/(n-1)*nr
    p=stringr(p)+list.unit[0]
    
    if widget_info(ID,/valid) then begin
        text[0]='FILENR.='+stringr(nr)+', COORD.='+p
        text[1]=['INT.='+stringr((*list.data)[list.NSelect,XX,YY])]
    endif
    endcase
endcase

XX=XX[sort(XX)]
YY=YY[sort(YY)]
nr=nr[sort(nr)]
if widget_info(ID,/valid) then widget_control,ID,set_value=text
end;pro CursorProp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SelectGroupMain,ev

ID=widget_info(ev.top,FIND_BY_UNAME='Update Main')
widget_control,ID,get_uvalue=bUpdate
if (bUpdate eq 0) then return

widget_control,ev.top,get_uvalue=list
if (list.main eq 0) then return
widget_control,list.main,get_uvalue=toplist
toplist.maskP=list.NSelect
widget_control,list.main,set_uvalue=toplist
RefreshDisplay,list.main,toplist

end;pro SelectGroupMain
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DDXDIplot, in

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
    image=(list.list2).data

    s=size(image)
    ncol=s[1]
    nrow=s[2]

    ;plotting
    ;    --------> X
    ; Y � [...]
    ;  \/
    xx=findgen(ncol)
    yy=findgen(nrow)


    ; ----Create the surface----
    oSurface = OBJ_NEW('IDLgrSurface', image,xx,yy , STYLE=(list.list2).style, SHADING=0, $
               COLOR=[60,60,255], BOTTOM=[64,192,128], $
               XCOORD_CONV=NormCoord((list.list2).xr,list.PosX),$
               YCOORD_CONV=NormCoord((list.list2).yr,list.PosY),$
               ZCOORD_CONV=NormCoord((list.list2).zr,list.PosZ))
    oGroup->Add, oSurface

    ; ----Display image in corner----
    ;oImage = OBJ_NEW('IDLgrImage', image,/order, INTERLEAVE=0)
     ;oGroup->Add, oImage

endelse

end;pro DDXDIplot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro INSERT_PROCESS_EVENTS,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list

if ev.type eq 5 or ev.type eq 6 then begin
    thisEvent='MOTION'
    bPress='NONE'
    if ev.press eq 1 then begin
        widget_control,ev.id,get_uvalue=list2
        bprint=0ull
        bmult=2ull
        borientation=4ull
        bexpand=8ull
        bcur=16ull
        case ev.type of
        5:    begin
            case ev.ch of
            42:    begin
                list2.mult+=list2.multperc
                bprint=bmult
                endcase
            47:    begin
                list2.mult-=list2.multperc
                bprint=bmult
                endcase
            43: begin
                list2.multperc*=10
                bprint=bmult
                endcase
            45: begin
                list2.multperc*=0.1
                bprint=bmult
                endcase
            104: begin
                if list2.orientation[0] ne 0 then begin
                    list2.orientation[0]=0
                    tmp=list2.MapRow
                    list2.MapRow=list2.MapCol
                    list2.MapCol=tmp
                    ;tmp=list2.curx
                    ;list2.curx=list2.cury
                    ;list2.cury=tmp
                    tmp=list2.hexpand
                    list2.hexpand=list2.vexpand
                    list2.vexpand=tmp
                    bprint=borientation or bexpand
                endif else bprint=borientation
                endcase
            118: begin
                if list2.orientation[0] ne 1 then begin
                    list2.orientation[0]=1
                    tmp=list2.MapRow
                    list2.MapRow=list2.MapCol
                    list2.MapCol=tmp
                    ;tmp=list2.curx
                    ;list2.curx=list2.cury
                    ;list2.cury=tmp
                    tmp=list2.hexpand
                    list2.hexpand=list2.vexpand
                    list2.vexpand=tmp
                    bprint=borientation or bexpand
                endif else bprint=borientation
                endcase
            108: begin
                if list2.orientation[1] ne 0 then begin
                    list2.orientation[1]=0
                    ;list2.curx=list2.MapCol-1-list2.curx
                    bprint=borientation
                endif else bprint=borientation
                endcase
            114: begin
                if list2.orientation[1] ne 1 then begin
                    list2.orientation[1]=1
                    ;list2.curx=list2.MapCol-1-list2.curx
                    bprint=borientation
                endif else bprint=borientation
                endcase
            116: begin
                if list2.orientation[2] ne 0 then begin
                    list2.orientation[2]=0
                    ;list2.cury=list2.MapRow-1-list2.cury
                    bprint=borientation
                endif else bprint=borientation
                endcase
            98: begin
                if list2.orientation[2] ne 1 then begin
                    list2.orientation[2]=1
                    ;list2.cury=list2.MapRow-1-list2.cury
                    bprint=borientation
                endif else bprint=borientation
                endcase
            106: begin
                list2.orientation[3]=0
                bprint=borientation
                endcase
            99: begin
                list2.orientation[3]=1
                bprint=borientation
                endcase
            54:    begin
                list2.hexpand+=1+round(0.1*list2.MapCol*(list2.expandfast ne 0))
                bprint=bexpand
                endcase
            52:    begin
                list2.hexpand-=1+round(0.1*list2.MapCol*(list2.expandfast ne 0))
                list2.hexpand>=1-list2.MapCol
                bprint=bexpand
                endcase
            50:    begin
                list2.vexpand+=1+round(0.1*list2.MapRow*(list2.expandfast ne 0))
                bprint=bexpand
                endcase
            56:    begin
                list2.vexpand-=1+round(0.1*list2.MapRow*(list2.expandfast ne 0))
                list2.vexpand>=1-list2.MapRow
                bprint=bexpand
                endcase
            53:    begin
                list2.expandfast=~list2.expandfast
                bprint=bexpand
                endcase
            27:    begin
                thisEvent='DOWN'
                bPress='RIGHT'
                endcase
            13:    begin
                thisEvent='DOWN'
                bPress='LEFT'
                endcase
            110:begin
                list.NSelect=0>(list.NSelect+1)<(list.NMaps-1)
                ID=widget_info(ev.top,FIND_BY_UNAME='NRslider')
                widget_control,ID,set_value=list.NSelect
                widget_control,ev.top,set_uvalue=list
                SelectGroupMain,ev
                RefreshWindowXDI,ev
                endcase
            112:begin
                list.NSelect=0>(list.NSelect-1)<(list.NMaps-1)
                ID=widget_info(ev.top,FIND_BY_UNAME='NRslider')
                widget_control,ID,set_value=list.NSelect
                widget_control,ev.top,set_uvalue=list
                SelectGroupMain,ev
                RefreshWindowXDI,ev
                endcase
            else: return
            endcase
            endcase
        6:    begin
            case ev.key of
            5:    list2.curx+=1+round(0.1*list2.MapCol*(ev.modifiers ne 0))
            6:    list2.curx-=1+round(0.1*list2.MapCol*(ev.modifiers ne 0))
            7:    list2.cury+=1+round(0.1*list2.MapRow*(ev.modifiers ne 0))
            8:    list2.cury-=1+round(0.1*list2.MapRow*(ev.modifiers ne 0))
            else: return
            endcase
            bprint=bcur
            endcase
        else: return
        endcase
        
        if bprint ne 0 then printw,list.OutID,''
        if (bprint and bmult) eq bmult then printw,list.OutID,'Insert image multiplier: '+string(list2.mult)
        if (bprint and bmult) eq bmult then printw,list.OutID,'Insert image multiplier-step: '+string(list2.multperc)
        if (bprint and borientation) eq borientation then printw,list.OutID,'Insert image orientation: '+(list2.orientation[0]?'Movev; ':'Moveh; ')+$
                (list2.orientation[1]?'RL; ':'LR; ')+$
                (list2.orientation[2]?'BT; ':'TB; ')+$
                (list2.orientation[3]?'Continue; ':'Jump; ')
        if (bprint and bexpand) eq bexpand then printw,list.OutID,'Insert image size expansion (hor x ver):'+string(list2.hexpand)+string(list2.vexpand)
        if (bprint and bcur) eq bcur then printw,list.OutID,'Insert cursor offset: '+string(list2.curx)+string(list2.cury)
        
        widget_control,ev.id,set_uvalue=list2
    endif
endif else begin
    widget_control,list.draw,/INPUT_FOCUS
    possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
    thisEvent = possibleEventTypes[ev.type]
    possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
    bPress = possibleButtons[ev.press]
endelse

; ----Display current cursor position----
X=ev.x
Y=ev.y
if (X lt 0) or (Y lt 0) or $
    X gt (list.XSize-1) or Y gt (list.YSize-1) then return
ID=widget_info(ev.top,FIND_BY_UNAME='position')
CursorProp,ID,X,Y,list,XX=XX,YY=YY

; ----Insert data----
(*list.data)=(*list.dataorg)
widget_control,ev.id,get_uvalue=list2
XX=XX[0]-list2.curx
YY=YY[0]-list2.cury

MapCol=list2.MapCol
MapRow=list2.MapRow
bcongrid=list2.hexpand ne 0 or list2.vexpand ne 0
if bcongrid then begin
    MapCol+=list2.hexpand
    MapRow+=list2.vexpand
endif

nx=(XX lt 0)?((MapCol+XX)<(list.MapCol)):((list.MapCol-XX)<(MapCol))
ny=(YY lt 0)?((MapRow+YY)<(list.MapRow)):((list.MapRow-YY)<(MapRow))
if (nx gt 0) and (ny gt 0) then begin
    bx=XX>0
    ex=(XX>0)+nx-1
    by=YY>0
    ey=(YY>0)+ny-1
    bx2=-(XX<0)
    ex2=-(XX<0)+nx-1
    by2=-(YY<0)
    ey2=-(YY<0)+ny-1

    data=*list2.data
    ReformXDIdatablock,data,[list.TRAN,list2.orientation[0]],[list.REV1,list2.orientation[1]],$
                            [list.REV2,list2.orientation[2]],[list.REV3,list2.orientation[3]],0b
    if bcongrid then begin
        data2=make_array(list.NMaps,MapCol,MapRow,value=data[0]*0)
        for i=0l,list.NMaps-1 do data2[i,*,*]=congrid(reform(data[i,*,*]),MapCol,MapRow,/interp)
        data=temporary(data2)
    endif
    
    (*list.data)[*,bx:ex,by:ey]=data[*,bx2:ex2,by2:ey2]*list2.mult
    if list2.stickyx then begin
        for i=0l,bx-1 do (*list.data)[*,i,by:ey]=data[*,0,by2:ey2]*list2.mult
        for i=ex+1,list.MapCol-1 do (*list.data)[*,i,by:ey]=data[*,MapCol-1,by2:ey2]*list2.mult
    endif
    if list2.stickyy then begin
        for i=0l,by-1 do (*list.data)[*,bx:ex,i]=data[*,bx2:ex2,0]*list2.mult
        for i=ey+1,list.MapRow-1 do (*list.data)[*,bx:ex,i]=data[*,bx2:ex2,MapRow-1]*list2.mult
    endif
    IF (thisEvent EQ 'DOWN') and bPress eq 'LEFT' and list2.bfiles and ~bcongrid THEN begin
        filenames=*list2.filenames
        ReformXDIdatablock,filenames,[list.TRAN,list2.orientation[0]],[list.REV1,list2.orientation[1]],$
                            [list.REV2,list2.orientation[2]],[list.REV3,list2.orientation[3]],0b
                            
        (*list.filenames)[*,bx:ex,by:ey]=filenames[*,bx2:ex2,by2:ey2]
        if list2.stickyx then begin
            for i=0l,bx-1 do (*list.filenames)[*,i,by:ey]=filenames[*,0,by2:ey2]
            for i=ex+1,list.MapCol-1 do (*list.filenames)[*,i,by:ey]=filenames[*,MapCol-1,by2:ey2]
        endif
        if list2.stickyy then begin
            for i=0l,by-1 do (*list.filenames)[*,bx:ex,i]=filenames[*,bx2:ex2,0]
            for i=ey+1,list.MapRow-1 do (*list.filenames)[*,bx:ex,i]=filenames[*,bx2:ex2,MapRow-1]
        endif
    endif
endif

RefreshDisplayXDI,list

IF thisEvent EQ 'DOWN' THEN BEGIN
    ; ----Save changes----
    if bPress eq 'LEFT' then (*list.dataorg)=(*list.data) $
    else (*list.data)=(*list.dataorg)

    ; ---Post processes---
    heap_free,list2
    WorkingDataCorrect,list
    RefreshDisplayXDI,list
    WIDGET_CONTROL, ev.top, set_uvalue=list
    modeXDI, ev, 'Dummy'
    return
endif
end;pro INSERT_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DEFAULT_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
if size(list,/type) eq 0 then return
X=ev.x
Y=ev.y
if (X lt 0) or (Y lt 0) or $
    X gt (list.XSize-1) or Y gt (list.YSize-1) then return
CursorProp,widget_info(ev.top,FIND_BY_UNAME='position'),X,Y,list
end;pro DEFAULT_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DEFAULTTOMO_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
if size(list,/type) eq 0 then return
X=ev.x
Y=ev.y
if (X lt 0) or (Y lt 0) or $
    X gt (TomoSize(list,/X)-1) or Y gt (TomoSize(list,/Y)-1) then return

CursorPropTomo,widget_info(ev.top,FIND_BY_UNAME='position'),X,Y,list
end;pro DEFAULTTOMO_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XDILEGEND_EVENTS,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.top,get_uvalue=list
if size(list,/type) eq 0 then return
X=ev.x
Y=ev.y
if (X lt 0) or (Y lt 0) or $
    X gt (list.XSize-1) or Y gt (list.YSize-1) then return
CursorProp,widget_info(ev.top,FIND_BY_UNAME='position'),X,Y,list

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]
if thisEvent ne 'DOWN' then return

list.legendpos=[X,Y]/float([list.XSize,list.YSize])
RefreshDisplayXDI,list
widget_control,ev.top,set_uvalue=list

end;pro XDILEGEND_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro UPDATE_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]

; ----Display current cursor position----
widget_control,ev.top,get_uvalue=list
if PTR_VALID(list.filenames) eq 0 then return
X=ev.x
Y=ev.y
if (X lt 0) or (Y lt 0) or $
    X gt (list.XSize-1) or Y gt (list.YSize-1) then return
CursorProp,widget_info(ev.top,FIND_BY_UNAME='position'),X,Y,list,XX=XX,YY=YY,nr=nr,_nrX=nrX,_nrY=nrY

IF thisEvent NE 'DOWN' THEN RETURN
error = WIDGET_INFO(list.main,/VALID_ID)
if error eq 0 then  return

filenames=reform(*list.filenames)
PixelJuggle,filenames,list,'c:\dummy.tif'

path=filenames[XX,YY]
if path eq '' then return

; Check for EDF files which contain a row
tmp=cutpath(path,ext=ext)
if ext eq '.edf' then $
    if EDFformat(path) eq 3 then nr=nrX

MainUpdateCore,list,list.main,path,nr=nr

widget_control,list.main,get_uvalue=toplist

end;pro UPDATE_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function blockprofile,p,x1i,x2i,y1i,y2i,block,blockx=blockx,dx=dx,dy=dy,round=round

; Spacing between profile points = 1 pixel

s=size(*p,/dimensions)
if n_elements(s) eq 2 then s=[s,1]
if n_elements(s) ne 3 then return,1b

x1=x1i[0]
x2=x2i[0]
y1=y1i[0]
y2=y2i[0]
if keyword_set(round) then begin
    x1=round(x1)
    x2=round(x2)
    y1=round(y1)
    y2=round(y2)
endif

D=makerange(0,sqrt((x2-x1)^2. + (y2-y1)^2.),1)
nPoints = n_elements(D)
if nPoints le 1 then return,0b

T = y2 - y1
N = x2 - x1
if T eq 0 then begin
    ; Horizontal
    if x1 gt x2 then D=-D
    xloc=x1+D
    yloc=replicate(y1,nPoints)
endif else if N eq 0 then begin
    ; Vertical
    if y1 gt y2 then D=-D
    xloc=replicate(x1,nPoints)
    yloc=y1+D
endif else begin
    ; Random
    if x1 gt x2 then D=-D
    m=T/float(N)
    D/=sqrt(1+m*m)
    xloc=x1+D
    yloc=y1+m*D
endelse

block=fltarr(s[0],nPoints)
for i=0l,s[0]-1 do block[i,*]= Interpolate(reform((*p)[i,*,*]), xloc, yloc)

if arg_present(blockx) then begin
    if not keyword_set(dx) then dx=1.
    if not keyword_set(dy) then dy=1.
    blockx=sqrt(((xloc-xloc[0])*dx)^2.+((yloc-yloc[0])*dy)^2.)
endif

return,0b
end;function blockprofile
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function closeregion,list
if list.nrROI ne 0 then $
              ind=total((*list.nrPoints)) else ind=0

nr=list.pverts-ind+1 ; number of points 1 more then number of vertices

case nr of
0:    return,0b
1:    begin
    list.pverts--
    if list.pverts eq -1 then begin
        ptr_free,list.xverts
        ptr_free,list.yverts
    endif else begin
        (*list.xverts)=(*list.xverts)[0:list.pverts]
        (*list.yverts)=(*list.yverts)[0:list.pverts]
    endelse

    return,1b
    endcase
else:    begin

    if list.nrROI eq 0 then begin
        ptr_free,list.nrPoints
        list.nrPoints=ptr_new(nr)
        ptr_free,list.highlight
        list.highlight=ptr_new(0b)
        ptr_free,list.names
        list.names=ptr_new('NoName')
    endif else begin
        (*list.nrPoints)=[(*list.nrPoints),nr]
        (*list.highlight)=[(*list.highlight),0b]
        (*list.names)=[(*list.names),'NoName']
    endelse
    list.nrROI++

    WIDGET_CONTROL, list.top, set_uvalue=list
    ID=widget_info(list.top,FIND_BY_UNAME='Highslider')
    widget_control,ID,set_slider_max=list.nrROI+1

    return,1b
    endelse
endcase

end;function closeregion
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddVert,list,X,Y
if list.pverts eq -1 then begin
    ptr_free,list.xverts
    ptr_free,list.yverts
    list.xverts=ptr_new(float(X))
    list.yverts=ptr_new(float(Y))
    n=n_elements(X)
endif else begin
    ind=where((X ne (*list.xverts)[list.pverts]) and (Y ne (*list.yverts)[list.pverts]),n)
    if n eq 0 then return
    (*list.xverts)=[(*list.xverts),X[ind]]
    (*list.yverts)=[(*list.yverts),Y[ind]]
endelse
list.pverts+=n
end;pro AddVert
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro printArray_meanstdev,arri,label,mult,outid
n=n_elements(arri)/2

x=arri[0,*]*mult
sx=arri[1,*]*mult

mx=total(x)/n
smx=sqrt(total(sx*sx))/n
        
diff=x-mx
FWHM_stdev=sqrt(total(diff*diff)/(n-1))
FWHM_stdev_s=sqrt(total(sx*sx*diff*diff)+smx*smx*(total(diff))^2.) / ((n-1)*FWHM_stdev)

printw,outid,'number of datapoints: '+stringr(n)
printw,outid,'Max '+label+': '+stringr(max(x))
printw,outid,'Min '+label+': '+stringr(min(x))

printw,outid,'<'+label+'>: '+rounderror(mx,smx)+' '+msymbols('micro')+'m'
printw,outid,'std'+label+': '+rounderror(FWHM_stdev,FWHM_stdev_s)+' '+msymbols('micro')+'m'

printw,outid,''
end;pro printArray_meanstdev
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro wplot2DgetFWHM,ev

widget_control,ev.top,get_uvalue=list
if list.ScanDim ne 1 then return

; Get X-Y
y=reform((*list.data)[list.NSelect,*,*])
n=n_elements(y)

xrange=wplot2Dgetxrange(list,n)
x=MakeRange(xrange[0],xrange[1],(xrange[1]-xrange[0])/(n-1))

; Fit step or gaussian
my=mean(y)
b1=y[0] lt my and y[n-1] gt my
b2=y[0] gt my and y[n-1] lt my
if b1 or b2 then begin
    maxy=max(y)
    if y[0] gt my and y[n-1] lt my then maxy=-maxy
    s=(x[n-1]-x[0])/8.
    A=[maxy,x[n/4.],s,min(y)]
    f=NLLS(x,y,A,'gfunc_stepfunc',sigma=sigma,chisq=chisq)
    FWHM=A[2]*2.*SQRT(2*ALOG(2))
    FWHMs=sigma[2]*2.*SQRT(2*ALOG(2))
endif else begin
    f=gaussfit(x,y,A,chisq=chisq,sigma=sigma,nterms=5)
    FWHM=A[2]*2.*SQRT(2*ALOG(2))
    FWHMs=sigma[2]*2.*SQRT(2*ALOG(2))
endelse

; Show result
RefreshDisplayXDI,list
wset,list.drawindex
oplot,x,f,color=100

ID=widget_info(list.top,FIND_BY_UNAME='Stat_FWHM')
widget_control,ID,set_value=rounderror(FWHM,FWHMs)+' '+list.unit[0]

end;pro wplot2DgetFWHM
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AddGaussFit,arr,pdata,idata,XX1,XX2,YY1,YY2,dx=dx,dy=dy,title=title,outid=outid,type=type

if blockprofile(pdata,XX1,XX2,YY1,YY2,block,blockx=x,dx=dx,dy=dy,/round) then return,0b
y=reform(block[idata,*])

nx=n_elements(x)
x2=makerange(x[0],x[nx-1],(x[nx-1]-x[0])/(10.*nx))

window,0

if not keyword_set(type) then type=0
case type of
0:    begin
    ; Gaussian FWHM
    f=gaussfit(x,y,A,chisq=chisq,sigma=sigma,nterms=5)
    
    gfunc_gauss,x2,A,f
    
    plot,x,y,psym=4,title=title,/xs
    oplot,x2,f
    plots,A[1]-A[2]*SQRT(2*ALOG(2)),!y.crange,color=100
    plots,A[1]+A[2]*SQRT(2*ALOG(2)),!y.crange,color=100
    
    ; reject
    if total(A[0:2] lt 0) ne 0 then begin
        A[0]=0
    endif
    
    i=2
    label='FWHM'
    mult=2.*SQRT(2*ALOG(2));sigma -> FWHM
    endcase
1:    begin
    ; Step function FWHM's
    maxy=max(y)
    s=(x[nx-1]-x[0])/8.
    A=[maxy,x[nx/4.],s,maxy,x[3*nx/4.],s,min(y)]
    f=NLLS(x,y,A,'gfunc_bellfunc',sigma=sigma,chisq=chisq)
    gfunc_bellfunc,x2,A,f
    
    plot,x,y,psym=4,title=title,/xs
    oplot,x2,f
    
    gfunc_gauss,x2,[(sqrt(2*!pi)*A[2])*(max(y)-A[6]),A[1],A[2],A[6]],tmp,struc=1
    oplot,x2,tmp,color=100
    gfunc_gauss,x2,[(sqrt(2*!pi)*A[5])*(max(y)-A[6]),A[4],A[5],A[6]],tmp,struc=1
    oplot,x2,tmp,color=100
    plots,A[[1,1]],!y.crange,color=100
    plots,A[[4,4]],!y.crange,color=100
    
    i=[[2],[5]]
    sigma[i]=0.1
    ;A[5]=A[2]
    
    ; reject
    b=A[2]<A[5]
    e=A[2]>A[5]
    if total(A[0:5] lt 0) ne 0 then begin ;abs((e-b)/A[2]) gt 0.5
        A[0]=0
    endif
    
    label='sigma'
    mult=1.
    endcase
2:    begin
    ; Wall thickness
    maxy=max(y)
    s=(x[nx-1]-x[0])/8.
    A=[maxy,x[nx/4.],s,maxy,x[3*nx/4.],s,min(y)]
    f=NLLS(x,y,A,'gfunc_bellfunc',sigma=sigma,chisq=chisq)
    gfunc_bellfunc,x2,A,f
    
    plot,x,y,psym=4,title=title,/xs
    oplot,x2,f
    
    gfunc_gauss,x2,[(sqrt(2*!pi)*A[2])*(max(y)-A[6]),A[1],A[2],A[6]],tmp,struc=1
    oplot,x2,tmp,color=100
    gfunc_gauss,x2,[(sqrt(2*!pi)*A[5])*(max(y)-A[6]),A[4],A[5],A[6]],tmp,struc=1
    oplot,x2,tmp,color=100
;    B=A
;    B[3]=0
;    gfunc_bellfunc,x2,B,tmp,struc=1
;    oplot,x2,tmp,color=100
;    B=A
;    B[0]=0
;    B[3]*=-1
;    B[5]*=-1
;    gfunc_bellfunc,x2,B,tmp,struc=1
;    oplot,x2,tmp,color=200
    
    plots,A[[1,1]],!y.crange,color=100
    plots,A[[4,4]],!y.crange,color=100
    
    i=4
    A[i]-=A[1]
    sigma[i]=0.1
    
    ; reject
    b=A[2]<A[5]
    e=A[2]>A[5]
    if total(A[0:5] lt 0) ne 0 then begin ;abs((e-b)/A[2]) gt 0.5 ; A[i] lt b
        A[0]=0
    endif
    
    label='Thickness(inflection points)'
    mult=1.
    endcase
endcase

; Add point
if n_elements(arr) eq 0 then arr=fltarr(2)
add=[A[i],sigma[i]]
arr=[[arr],[add]]

; Print new point(s)
printw,outid,'Added point:'
printArray_meanstdev,add,label,mult,outid

; Keep only valids
ind=where(arr[0,*] gt 0,ct)
if ct eq 0 then return,0b
arr2=arr[*,ind]

; Reject outliers
if ct ge 2 then begin
    k=10.
    limit=k*stdev(arr2[0,*],m)
    ind=where(abs(arr2[0,*]-m) lt limit,ct)
    if ct eq 0 then return,0b
    window,1
    plot,arr2[0,*],yrange=m+[-limit,limit]
    plots,!x.crange,m+[limit,limit],color=100
    plots,!x.crange,m-[limit,limit],color=100
    oplot,ind,arr2[0,ind],psym=1,color=100
    
    arr2=arr2[*,ind]
endif

; Print
printArray_meanstdev,arr2,label,mult,outid

return,1b
end;function AddGaussFit
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro updatexdicoord,list,dx=dx,dy=dy
coord=getxdicoord(list,list.coord)
dpix=abs(coord[1,*]-coord[0,*])*1000.
ncol=list.MapCol-1
nrow=list.MapRow-1
dpix[0]/=ncol
dpix[1]/=nrow
dx=float(promptnumber(stringr(dpix[0]),list.top,'X pixel size in micrometer:'))
dy=float(promptnumber(stringr(dpix[1]),list.top,'Y pixel size in micrometer:'))
coord[1,*]=coord[0,*]+[dx*ncol,dy*nrow]/1000.
coord[*,0]=coord[sort(coord[*,0]),0]
coord[*,1]=coord[sort(coord[*,1]),1]
list.coord=coord
end;pro updatexdicoord
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GaussFitROIs,list,type=type

if list.nrROI eq 0 then return

updatexdicoord,list,dx=dx,dy=dy

i0=0
for i=0l,list.nrROI-1 do begin
    i1=i0+(*list.nrPoints)[i]-1
    
    if (*list.nrPoints)[i] eq 2 then begin
        tmp=AddGaussFit(arr,list.data,list.NSelect,(*list.xverts)[i0],(*list.xverts)[i1],$
            (*list.yverts)[i0],(*list.yverts)[i1],title='ROI index'+string(i+1),dx=dx,dy=dy,outid=list.outid,type=type)
        ;wait,1
    endif
    i0=i1+1
endfor

end;pro GaussFitROIs
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro EXTRACT1D_PROCESS_EVENTS,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]
possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
thisEvent = thisEvent+possibleButtons[ev.press]

; ----Display current cursor position----
widget_control,ev.top,get_uvalue=list
if list.ScanDim ne 2 then return
X=ev.x
Y=ev.y
if (X lt 0) or (Y lt 0) or $
    X gt (list.XSize-1) or Y gt (list.YSize-1) then return
IDpos=widget_info(ev.top,FIND_BY_UNAME='position')
CursorProp,IDpos,X,Y,list

case thisEvent of
'DOWNLEFT': begin
            list.xs=X
            list.ys=Y
            Widget_Control, ev.top, Set_UValue=list
            endcase
'MOTIONNONE': begin
            WSet, list.drawindex
            h=!D.x_size
            v=!D.y_size
            Device, Copy = [0,0, h,v, 0,0,list.pixindex]
            x = [list.xs, X]
            y = [list.ys, Y]
            
            case ev.modifiers of
            1: begin;shift
                x=[0,h-1]
                y[0]=y[1]
               endcase
            2: begin;ctrl
                y=[0,v-1]
                x[0]=x[1]
               endcase
            else:
            endcase
            
            x=0>x<(h-1)
            y=0>y<(v-1)
            PlotS, x, y,Linestyle=0,/device
            endcase
'DOWNRIGHT': begin
            widget_control,ev.id,get_uvalue=type
            WSet, list.drawindex
            h=!D.x_size
            v=!D.y_size
            X=[X,list.xs]
            Y=[Y,list.ys]
            
            case ev.modifiers of
            1: begin;shift
                x=[0,h-1]
                y[1]=y[0]
               endcase
            2: begin;ctrl
                y=[0,v-1]
                x[1]=x[0]
               endcase
            else:
            endcase
            
            ; ----Add to vertices----
            if type eq 2 then begin
                AddVert,list,StatTR(X,list.sfh),StatTR(list.YSize-Y,list.sfv)
                if closeregion(list) then RefreshDisplayXDI,list
                return
            endif

            ; ----Select file----
            xdifile=SelFile(list.path,list.file,['*.xdi','*.chi','*.tiff'],'Save To....')
            if xdifile eq '' then return
            tmp=CutPath(xdifile,ext=ext)

            ; ----Create data----
            Result = DIALOG_MESSAGE('Save with corrections?',/question)
            if Result eq 'Yes' then p=list.data else p=list.dataorg
            
            tmp=CurPos(X,Y,list.MapCol,list.MapRow,list.sfh,list.sfv,XX=XX,YY=YY)
            updatexdicoord,list,dx=dx,dy=dy
            widget_control,ev.top,set_uvalue=list
            if blockprofile(p,XX[0],XX[1],YY[0],YY[1],block,dx=dx,dy=dy,blockx=x,/round) then return

            ; ----Save XDI----
            case strlowcase(ext) of
            '.chi':    begin
                    Type3CHI,xdifile,block,list,x=x
                    endcase
            '.xdi':    begin
                    names=strarr(list.NMaps)
                    for i=0l,list.NMaps-1 do begin
                        labeli=0>list.labeli<(n_elements(((*list.pmstruc).(i+1)).data)-1)
                        names[i]=stringr(((*list.pmstruc).(i+1)).data[labeli])
                    endfor
                    Type3XDI,xdifile,block,names,1,list.OutID,ounit=0,ocoord=x[[0,n_elements(x)-1]]/1000.
                    endcase
            else:    begin
                    Type3TIFF,xdifile,block,list.outid,x=x
                    endelse
            endcase
            
            endcase
else:
endcase

end;pro EXTRACT1D_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LOF_PROCESS_EVENTS, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]

; ----Display current cursor position----
widget_control,ev.top,get_uvalue=list
X=ev.x
Y=ev.y
if (X lt 0) or (Y lt 0) or $
    X gt (list.XSize-1) or Y gt (list.YSize-1) then return
CursorProp,widget_info(ev.top,FIND_BY_UNAME='position'),X,Y,list

IF thisEvent NE 'DOWN' THEN RETURN
list.xs=ev.x
list.ys=ev.y
Widget_Control, ev.id, Event_Pro='LOF_DRAWBOX', $
    draw_motion=1
Widget_Control, ev.top, Set_UValue=list
end;pro LOF_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LOF_DRAWBOX ,ev

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
    h=!D.x_size
    v=!D.y_size
    Device, Copy = [0,0, h,v, 0,0,list.pixindex]
    Widget_Control, ev.id, Event_Pro='LOF_PROCESS_EVENTS'
    x = [list.xs, ev.x]
    y = [list.ys, ev.y]
    IF list.xs GT ev.x THEN x = [ev.x, list.xs]
    IF list.ys LT ev.y THEN y = [ev.y, list.ys]
    x=0>x<(h-1)
    y=0>y<(v-1)

    ; upper-left and lower-right corner:
    ; 0-----> X
    ; |
    ; |     *------
    ; |     |      |
    ; V  Y  |      |
    ;       |      |
    ;        ------*

    ; ----Make file indices (current orientation)----
    CursorProp,0L,x,y,list,XX=XX,YY=YY
    nx=XX[1]-XX[0]+1
    ny=YY[1]-YY[0]+1
    if nx le 0 and ny le 0 then return
    xind=rebin(lindgen(nx)+XX[0],nx,ny,/sample)
    yind=rebin(lindgen(1,ny)+YY[0],nx,ny,/sample)
    ind=yind*list.MapCol+xind

    widget_control,ev.id,get_uvalue=list2
    case list2.type of
    1:    begin

        header=strarr(4)
        header[0]='$# columns:'
        header[1]=stringr(nx)
        header[2]='$# rows:'
        header[3]=stringr(ny)

        error=CutPath(list.file,file=file)
        filenames=reform(*list.filenames)
        tmp=where(filenames eq '',ct)
        if ct ne 0 then filenames[tmp]='c:\dummy.tif'
        Result = DIALOG_MESSAGE('Save with corrections?',/question)
        if Result eq 'Yes' then PixelJuggle,filenames,list,'c:\dummy.tif'

        error=WriteLOF(list.path,file+'.lof',ind,ev.top,filenames,list.OutID,header=header)
        if error eq 0 then printw,list.OutID,'Couldn`t save .lof file!!!'

        Widget_Control, ev.top, Set_UValue=list
        return
        endcase
    2:    begin
        tmp=(*list.data)[list.NSelect,*,*]
        val=float(mean(tmp[ind]))
        tmp=PromptNumber(stringr(val),ev.top,'Enter replacing value:')
        reads,tmp,val
        
        ; Convert ind to default orientation
        ind=ConvertFileNumber(xind,yind,list.MapCol,list.MapRow,list.TRAN,list.REV1,list.REV2,list.REV3)
        nind=n_elements(ind)
        
        tmp=[reform(ind,1,nind),replicate(val,1,nind)]
        
        if ~ptr_valid(list.preplacep) then list.preplacep=ptr_new(ptrarr(list.NMaps,/alloc)) else $
        if n_elements(*list.preplacep) ne list.NMaps then begin
            heap_free,*list.preplacep
            *list.preplacep=ptrarr(list.NMaps,/alloc)
        endif
        
        for i=0l,list.NMaps-1 do begin
            p=(*list.preplacep)[i]
            if n_elements(*p) eq 0 then *p=tmp $
            else *p=[[*p],[tmp]]
        endfor

        WorkingDataCorrect,list
        widget_control,ev.top,set_uvalue=list
        RefreshDisplayXDI,list
        RefreshWindowXDI,ev
        return
        endcase
    3:    begin
        id=find_by_uname(list2.top,'col')
        if widget_info(id,/valid) then begin
            widget_control,id,set_value=stringr(nx)
        endif
        id=find_by_uname(list2.top,'row')
        if widget_info(id,/valid) then begin
            widget_control,id,set_value=stringr(ny)
        endif
        id=find_by_uname(list2.top,'coloff')
        if widget_info(id,/valid) then begin
            widget_control,id,set_value=stringr(XX[0])
        endif
        id=find_by_uname(list2.top,'rowoff')
        if widget_info(id,/valid) then begin
            widget_control,id,set_value=stringr(YY[0])
        endif
        endcase
    else: return
    endcase
endif

WSet, list.drawindex
h=!D.x_size
v=!D.y_size
Device, Copy = [0,0, h,v, 0,0,list.pixindex]
x = [list.xs, ev.x]
y = [list.ys, ev.y]
x=0>x<(h-1)
y=0>y<(v-1)
;DEVICE, GET_GRAPHICS = oldg, SET_GRAPHICS=6          ;Xor drawing mode
PlotS, [x[0], x[0]], y,Linestyle=0,/device
PlotS, [x[1], x[1]], y,Linestyle=0,/device
PlotS, x,[y[0], y[0]],Linestyle=0,/device
PlotS, x,[y[1], y[1]],Linestyle=0,/device
;DEVICE, SET_GRAPHICS=oldg          ;Copy mode
Widget_Control, ev.top, Set_UValue=list

; ----Display current cursor position----
X=ev.x
Y=ev.y
if (X lt 0) or (Y lt 0) or $
    X gt (list.XSize-1) or Y gt (list.YSize-1) then return
CursorProp,widget_info(ev.top,FIND_BY_UNAME='position'),X,Y,list
end;pro LOF_DRAWBOX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LOF_From_PixelInt,ev

Widget_Control, ev.top, Get_UValue=list

img=(*list.data)[list.NSelect,*,*]

;thres='> '+stringr(min(img))
thres='> 50%'
thres=PromptNumber(thres,ev.top,'Pixel intensity (value or %) in LOF:')

p=strpos(thres,'%')
if p ne -1 then begin
    thres=strmid(thres,0,p)
    ind=where(byte(thres) ge 48 and byte(thres) le 57,ct)
    if ct eq 0 then return
    mi=min(img,max=ma,/nan)
    thres=strmid(thres,0,ind[0])+stringr(mi+(ma-mi)*float(strmid(thres,ind[0]))/100)
endif

ind=compare(img,thres)
wset,list.drawindex
tvscl,congrid(reform(ind),list.XSize,list.YSize,interp=0),/order
ind=where(ind,n)

header=strarr(4)
header[0]='$# columns:'
header[1]=stringr(n)
header[2]='$# rows:'
header[3]='1'

error=CutPath(list.file,file=file)
filenames=reform(*list.filenames)
tmp=where(filenames eq '',ct)
if ct ne 0 then filenames[tmp]='c:\dummy.tif'
Result = DIALOG_MESSAGE('Save with corrections?',/question)
if Result eq 'Yes' then PixelJuggle,filenames,list,'c:\dummy.tif'

error=WriteLOF(list.path,file+'.lof',ind,ev.top,filenames,list.OutID,header=header)
if error eq 0 then printw,list.OutID,'Couldn`t save .lof file!!!'

Widget_Control, ev.top, Set_UValue=list
end;pro LOF_From_PixelInt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DRAW_POL, x,y,highlight, FILL = fill

; Draw the outline (or polygon if FILL is set)
; Use the XOR drawing mode.

DEVICE, GET_GRAPHICS = oldg, SET_GRAPHICS=6          ;Xor drawing mode
col = 0;*16
;while col lt !d.table_size do col = col + col
;if (col ge !d.table_size) then col = !d.table_size - 1.

if highlight eq 1 then thick = 4. else thick = 1.5
y2=!d.Y_SIZE-1-y
plots,x,y2,COLOR=col,/device,THICK=thick
IF KEYWORD_SET(FILL) then POLYFILL, x,y2, COLOR=col,/device

DEVICE, SET_GRAPHICS=oldg          ;Copy mode
end;pro DRAW_POL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ROI_PROCESS_EVENTS,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]
possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
thisEvent = thisEvent+possibleButtons[ev.press]

widget_control,ev.top,get_uvalue=list
X=ev.x
Y=ev.y
if (X lt 0) or (Y lt 0) or $
    X gt (list.XSize-1) or Y gt (list.YSize-1) then return
if (list.MState eq 'DOWNLEFT') and (thisEvent ne 'UPNONE') then $
 thisEvent='DOWNLEFT' else list.MState=''

case thisEvent of
'DOWNLEFT': begin
          list.MState=thisEvent
          ; ----Add point to vertices----
          AddVert,list,StatTR(X,list.sfh),StatTR(list.YSize-Y,list.sfv)
          ; ----Clear ROI in process----
          WSet, list.drawindex
          Device, Copy = [0,0, list.XSize,list.YSize, 0,0,list.pixindex]
          ; ----Display ROI in process----
          if list.nrROI ne 0 then $
          ind=total((*list.nrPoints)) else $
          ind=0
          DRAW_POL, RTStat((*list.xverts)[ind:list.pverts],list.sfh),RTStat((*list.yverts)[ind:list.pverts],list.sfv),0
         endcase
'DOWNRIGHT':begin
         ; ----Remove last point----
         if list.nrROI ne 0 then $
          ind=total((*list.nrPoints)) else $
          ind=0
         if list.pverts ge ind then begin
             if list.pverts ne 0 then begin
                  (*list.xverts)=(*list.xverts)[0:list.pverts-1]
                  (*list.yverts)=(*list.yverts)[0:list.pverts-1]
              endif else begin
                  ptr_free,list.xverts
                  ptr_free,list.yverts
              endelse
          list.pverts--
          ; ----Clear ROI in process----
          WSet, list.drawindex
              Device, Copy = [0,0, list.XSize,list.YSize, 0,0,list.pixindex]
              ; ----Display ROI in process----
              if list.pverts ge ind then $
              DRAW_POL, RTStat((*list.xverts)[ind:list.pverts],list.sfh),RTStat((*list.yverts)[ind:list.pverts],list.sfv),0
         endif
         endcase
else: a=0B
endcase
widget_control,ev.top,set_uvalue=list

; ----Display current cursor position----
CursorProp,widget_info(ev.top,FIND_BY_UNAME='position'),X,Y,list

end;pro ROI_PROCESS_EVENTS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RemakePatternBlock,spe,list,nx
GetXDIdatadim,list,nmaps,ncol,nrow,/default
if ncol*nrow ne n_elements(spe)/nx then return

spe=reform(spe,nx,ncol,nrow,/overwrite)
ReformXDIDatablock,spe,[0b,list.TRAN],[0b,list.REV1],[0b,list.REV2],[0b,list.REV3],0b
end;pro RemakePatternBlock
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro VoxelPattern,j,list
; THIS ROUTINE IS NOT CORRECT ANYMORE

if PTR_VALID(list.filenames) eq 0 then return
if WIDGET_INFO(list.main,/VALID_ID) eq 0 then return

if FileFormat(0,(*list.filenames)[0]) ne 3 then return

; Read data and reform
tmp=CutPath((*list.filenames)[0],path=path,file=file,ext=ext)
error=ReadScan(path,file+ext,ext,list.OutID,xval=xval,spe=spe)
RemakePatternBlock,spe,list,n_elements(xval)

; Make voxel-pattern
nx=n_elements(xval)
spevoxel=make_array(nx,value=spe[0])

time=systime(1)
progprev=0
pos=0
window
for i=0l,nx-1 do begin
    image=spe[*,*,i]

    TOMOrecon,image,list,tomogram=tomogram;,average=list.tomoavg
    spevoxel[i]=tomogram[j]
    
    ; Progress
    prog=fix(i*100./nx)
    if prog ne progprev then begin
        tvscl,image,pos,/order
        tvscl,tomogram,pos+1,/order
        pos=(pos+2) mod 6
        progprev=prog
        timeleft=(systime(1)-time)*(100.-prog)/prog
        ConvertTime,timeleft,timeaf
        str='Progress: '+stringr(prog)+'% (left: '+stringr(timeleft)+timeaf
            
        if prog eq 1 then $
            if DIALOG_MESSAGE([str,'Proceed?'],/question) eq 'No' then return
        printw,list.OutID,str
    endif
endfor

plot,xval,spevoxel
end;pro VoxelPattern
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro UPDATE_PROCESS_EVENTS2, ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

; Update position info
widget_control,ev.top,get_uvalue=list
CursorPropTomo,widget_info(ev.top,FIND_BY_UNAME='position'),ev.x,ev.y,list,XX=XX,YY=YY

; Draw sinus
rho=RhoFromTomogram(ev.x,ev.y,list,nviews=nviews,bsinovert=bsinovert)
wset,list.drawindex
Device, Copy = [0,0, list.XSize,list.YSize, 0,0,list.pixindex]
DRAWSINUS,rho,nviews,bsinovert,list

; Check event
possibleEventTypes = [ 'DOWN', 'UP', 'MOTION', 'SCROLL' ]
thisEvent = possibleEventTypes[ev.type]
IF thisEvent NE 'DOWN' THEN RETURN

possibleButtons = ['NONE', 'LEFT', 'MIDDLE', 'NONE', 'RIGHT']
bPress = possibleButtons[ev.press]
if bPress eq 'RIGHT' then begin
    modeXDI, ev, 'Dummy'
    return
endif

; Make pattern
;print,xx,yy
;j=CurPos(ev.x,ev.y,list.MapCol,list.MapRow,list.sfh,list.sfv,XX=XX,YY=YY)
;print,xx,yy
;VoxelPattern,j,list
;return

if PTR_VALID(list.filenames) eq 0 then return
if WIDGET_INFO(list.main,/VALID_ID) eq 0 then return

weights=rho
XX=0>[[floor(weights)],[ceil(weights)]]<(list.MapCol-1)
weights-=XX[*,0]
weights>=0
weights=[[weights],[1-weights]]
YY=indgen(list.MapRow)
YY=[[YY],[YY]]
nr=ConvertFileNumber(XX,YY,list.MapCol,list.MapRow,list.TRAN,list.REV1,list.REV2,list.REV3)

filenames=reform(*list.filenames)
PixelJuggle,filenames,list,'c:\dummy.tif'
path=filenames[nr]
if path[0] eq '' then return
MainUpdateCoreCalc,list,list.main,path,nr=nr,weights=weights
widget_control,list.main,get_uvalue=toplist

end;pro UPDATE_PROCESS_EVENTS2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro OpenXDI,ev,prompt,addxdi=addxdi,format=format,path=path

widget_control,ev.top,get_uValue=list

if keyword_set(addxdi) and (list.listtype eq 20) then addxdi=0

; ----Choose file, with or without dialog----
callagain=0b
if prompt then begin
    path=''
    path=DIALOG_PICKFILE(path=list.path,file=list.file,filter='*.xdi',title='Read XDI...',/multi)
    if path[0] EQ '' then return

    if keyword_set(format) then begin
        format=PromptNumber(format,ev.top,'Number of characters for each value:')
        if format eq '' then delvar2,format else list.openformat=format
    endif

    callagain=n_elements(path) gt 1
    if callagain then begin
        pathagain=path[1:*]
        path=path[0]
    endif
endif
error=CutPath(path,path=path,file=file,ext=ext)
file+=ext

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,list.OutID,'Error in reading '+path+file
    return
ENDIF

; ----Read file----
struct=ReadXDI(path,file,list.OutID,error=error,format=format)
if error then begin
    temp=dialog_message('Did you try Options > parameters > open formatted?',/info)
    return
endif

; ----Edit list----
if list.listtype eq 20 then begin
    list.top=ev.top
    
    ID=widget_info(ev.top,FIND_BY_UNAME='drawgradient')
    widget_control,ID,get_value=drawindexgradient
    list.drawindexgradient=drawindexgradient
endif

if keyword_set(addxdi) then begin
    if list.ScanDim ne struct.ScanDim then begin
        printw,list.OutID,"XDI hasn't got the same dimensions."
        return
    endif
    
    data=struct.data
    ReformXDIDatablock,data,[struct.orientation[0],list.TRAN],[struct.orientation[1],list.REV1],[struct.orientation[2],list.REV2],[struct.orientation[3],list.REV3],0b

    ; Add data
    s=size(data,/dim)
    (*list.dataorg)=[(*list.dataorg),temporary(data)]
    (*list.data)=(*list.dataorg)
    list.NMaps+=s[0]

    ; Add info
    pmtags=lindgen(s[0])+max(*list.pmtags)+1
    (*list.pmtags)=[(*list.pmtags),pmtags]
    for i=1,s[0] do $
        (*list.pmstruc)=create_struct((*list.pmstruc),'t'+stringr(pmtags[i-1]),struct.pmstruc.(i))

endif else begin
    list.path=path
    list.file=file
    ptr_free,list.data
    ptr_free,list.dataorg
    list.data=ptr_new(struct.data)
    list.dataorg=ptr_new(struct.data)
    ptr_free,list.filenames
    list.filenames=((size(struct.filenames,/type) eq 7)?PTR_NEW(struct.filenames):PTR_NEW())
    ptr_free,list.pmstruc
    list.pmstruc=ptr_new(struct.pmstruc)
    ptr_free,list.pmtags
    list.pmtags=ptr_new(struct.pmtags)
    ptr_free,list.pmfields
    list.pmfields=ptr_new(struct.pmfields)
    list.npmfields=struct.npmfields
    list.ScanDim=struct.Scandim
    header=struct.header

    ; ----Prepare Windows and extract data ordening----
    s=dimsize(*list.data,3)
    list.NMaps=s[0]
    list.MapCol=s[1]
    list.MapRow=s[2]
    list.TRAN=struct.orientation[0]
    list.REV1=struct.orientation[1]
    list.REV2=struct.orientation[2]
    list.REV3=struct.orientation[3]
    list.Hsign=struct.sign[0]
    list.Vsign=struct.sign[1]
    
    if widget_info(list.TomoReconID,/VALID_ID) eq 1 then widget_control,list.TomoReconID,/destroy
    list.TomoReconID=0L
    list.TomoReconIndex=0L
    
    if list.TRAN then sign=list.Vsign else sign=list.Hsign
    header_scan=str_sep(strtrim(strcompress(header[0]),2),' ')
    case list.ScanDim of
    2:    begin
        if list.listtype eq 20 then begin
            device,get_screen_size=screen
            screen[0]*=0.3
            screen[1]*=0.6
            list.sfh>=2*float((list.MapCol>list.MapRow))/(screen[0]<screen[1])
            list.sfv>=list.sfh
        endif
        
        tmpp=widget_info(ev.top,find_by_uname='Magx')
        widget_control,tmpp,set_value=stringr(1./list.sfh)
        tmpp=widget_info(ev.top,find_by_uname='Magy')
        widget_control,tmpp,set_value=stringr(1./list.sfv)

        ; coord and unit in default orientation: moveh,LR,UB,LR(+),UB(+)
        tmp2=str_sep(strtrim(strcompress(header[2]),2),' ')
        tmp3=str_sep(strtrim(strcompress(header[3]),2),' ')
        list.unit=[tmp2[2],tmp3[2]]
        list.coord=float([[header_scan[2:3]],[header_scan[6:7]]])
        if list.TRAN then begin
            list.coord=reverse(list.coord,2)
            list.unit=reverse(list.unit)
        endif
        list.coord[*,0]=list.coord[sort(list.coord[*,0]),0]
        list.coord[*,1]=list.coord[sort(list.coord[*,1]),1]

        list.expo=header_scan[9]
        endcase
    1:    begin
        ; coord and unit in default orientation: moveh,LR,LR(+)
        tmp2=str_sep(strtrim(strcompress(header[2]),2),' ')
        list.unit=[tmp2[2],'']
        list.coord=float(header_scan[2:3])
        list.coord[*,0]=list.coord[sort(list.coord[*,0]),0]

        list.expo=header_scan[5]
        endcase
    endcase
    
    ind=where(list.unit eq msymbols('micro')+'m',ct)
    if ct ne 0 then list.unit[ind]='mm'
    
    list.normflag=0b
    xdi_changesize,list
endelse

; ----Set slider for selecting propper map----
list.NSelect<=list.NMaps-1
ID=widget_info(ev.top,FIND_BY_UNAME='NRslider')
widget_control,ID,set_slider_max=list.NMaps-1,set_value=list.NSelect

; ----Set slider for selecting propper region----
if list.nrROI ne 0 then begin
ID=widget_info(ev.top,FIND_BY_UNAME='Highslider')
widget_control,ID,set_value=0,set_slider_max=list.nrROI+1
ID=widget_info(ev.top,FIND_BY_UNAME='Names')
widget_control,ID,set_value=''
endif

; ----Finish----
list.listtype=21

WorkingDataCorrect,list
widget_control,ev.top,Set_UValue=list
if callagain then begin
    for i=0l,n_elements(pathagain)-1 do begin
        OpenXDI,ev,0,/add,format=format,path=pathagain[i]
    endfor
endif
widget_control,ev.top,TLB_SET_TITLE='XDI Editor: '+list.file
RefreshDisplayXDI,list
RefreshWindowXDI,ev
RefreshGradientImage,ev

end;pro OpenXDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SaveAllXDI,ev,current=current,noprompt=noprompt

WIDGET_CONTROL, ev.top, get_uvalue=list

ID=widget_info(list.top,FIND_BY_UNAME='Save Groups Separate')
widget_control,ID,get_uvalue=SaveSep
ID=widget_info(list.top,FIND_BY_UNAME='Save without boarders')
widget_control,ID,get_uvalue=NoBoard
if SaveSep or keyword_set(current) then NoBoard=1b
ID=widget_info(list.top,FIND_BY_UNAME='Save 1D data in 1 plot')
widget_control,ID,get_uvalue=OnePlot

; ----Prompt for save sequence----
path = list.path
file = strmid(list.file,0,STRLEN(list.file)-4)
forcepath=0b

bcurrent=keyword_set(current)
breplot=~bcurrent or keyword_set(noprompt)
if bcurrent then begin
    NMaps=1
    NSelect=list.NSelect
    if size(current,/type) eq 7 then begin
        forcepath=1b
        t=CutPath(current,path=path,file=file,ext=ext)
    endif
endif else begin
    if keyword_set(noprompt) then begin
        NMaps=list.NMaps
        NSelect=lindgen(NMaps)
    endif else NMaps=GetGroupindices(list.NMaps,ev.top,indices=NSelect)
endelse
if (NMaps eq 0) or (NMaps gt list.NMaps) or (max(NSelect) gt list.NMaps-1) or (min(NSelect) lt 0) then begin
    printw,list.OutID,'Wrong save index sequence.'
    return
endif

; ----Save all 1D in 1 plot----
if (SaveSep eq 0) and (list.ScanDim eq 1) and (OnePlot eq 1) then begin

    window,xsize=600,ysize=800;,/pixmap
    wplot2D1,list,NSelect=NSelect

    SaveImage,path,file+'.tif','wplot2D1',list,OutID=list.OutID,NSelect=NSelect,noprompt=noprompt

    wdelete,!d.window
    return
endif

; ----Select output----
if ~breplot then SaveSep=0
bTomo=widget_info(list.TomoReconID,/VALID_ID) eq 1
if bTomo then begin
    ID=widget_info(list.top,FIND_BY_UNAME='Save Tomo Separate')
    widget_control,ID,get_uvalue=SaveTomoSep
    NewDevice=SaveSep and SaveTomoSep
endif else begin
    SaveTomoSep=0
    NewDevice=SaveSep
endelse

if ~forcepath then begin
    if keyword_set(noprompt) then begin
        ext='.tif'
    endif else begin
        path=SelFile(path,file+'.tif',SaveImageFilter(copy=~NewDevice),'Save Image....')
        if path eq '' then return
        t=CutPath(path,path=path,file=file,ext=ext)
    endelse
endif

if bTomo then fileTomo=file+'Recon'

; ----Select output device----
if breplot then window,xsize=list.XSize,ysize=list.YSize,/pixmap,/free $
else wset,list.drawindex
index=!d.window

if bTomo then begin
    if breplot then window,xsize=TomoSize(list,/X), ysize=TomoSize(list,/Y),/pixmap,/free $
    else wset,list.TomoReconIndex
    list.TomoReconIndex=!d.window
endif

; ----Prepare 1D or 2D offsets----
if SaveSep eq 0 then begin
    xs=list.nhorizontal<NMaps
    ys=ceil(float(NMaps)/xs)
    nempty=xs*ys-NMaps
    xsize=list.XSize
    ysize=list.YSize
    xsize2=TomoSize(list,/X)
    ysize2=TomoSize(list,/Y)
    if NoBoard eq 0 then begin
        xsize+=2
        ysize+=2
        xsize2+=2
        ysize2+=2
        Boff=1
    endif else Boff=0
    boardercol=list.boardercol

    if bTomo then begin
        if SaveTomoSep eq 1 then begin
            imgS2=bytarr(3,xsize2*xs,ysize2*ys)+boardercol
            xoff2=((xsize2*lindgen(xs))#replicate(1B,ys))+Boff
            yoff2=(replicate(1B,xs)#(xsize2*reverse(lindgen(ys))))+Boff
        endif else begin
            xsize>=xsize2
            ysize+=ysize2-(NoBoard eq 0) ; -1 causes the boarder to be thin between tomo and sinogram
        endelse
    endif
    imgS=bytarr(3,xsize*xs,ysize*ys)+boardercol
    xoff=((xsize*lindgen(xs))#replicate(1B,ys))+Boff
    yoff=(replicate(1B,xs)#(ysize*reverse(lindgen(ys))))+Boff
    if NoBoard eq 0 then begin
        xsize-=2
        ysize-=2
        xsize2-=2
        ysize2-=2
    endif
endif else begin
    nempty=0
    if SaveTomoSep eq 0 then begin
        xsize=list.XSize>TomoSize(list,/X)
    endif
endelse

; ----HTML wrapper----
bhtmlwrapper=0b
subdir=''
if SaveSep eq 1 and NMaps gt 1 then begin
    if keyword_set(noprompt) then bhtmlwrapper=1b $
    else bhtmlwrapper=dialog_message('Produce HTML wrapper?',/question,/default_no) eq 'Yes'
    if bhtmlwrapper then begin
        if keyword_set(noprompt) then htmlfile=path+file+'.html' $
        else htmlfile=SelFile(path,file+'.html','*.html','Save HTML wrapper....')
        xs=list.nhorizontal
        ys=ceil(float(NMaps)/xs)
        ymag=2+(bTomo and SaveTomoSep)
        table=strarr(xs,ymag*ys)
        tableindex=0L
    endif
    subdir=file+path_sep()
    file_mkdir,path+subdir
endif

; ----Save 1D or 2D----
nameadd=MakeNumber(NSelect)
for i=0l,NMaps-1 do begin
    list.NSelect=NSelect[i]
    labeli=0>list.labeli<(n_elements(((*list.pmstruc).(list.NSelect+1)).data)-1)
    nr=nameadd[i]+'_'+string(((*list.pmstruc).(list.NSelect+1)).data[list.labeli])
    fileOut=subdir+file+'_'+nr+ext
    
    ; Plot and save
    if NewDevice then begin
        
        wset,index
        path2=path
        case list.ScanDim of
        1:     SaveImage,path2,fileOut,'wplot2D',list,OutID=list.OutID,/forcereplot,/noprompt
        2:    begin
            SaveImage,path2,fileOut,'wplot3D',list,OutID=list.OutID,/tomosep,/forcereplot,/noprompt
            if bTomo then begin
                path2=path
                fileTomoOut=subdir+fileTomo+'_'+nr+ext
                wset,list.TomoReconIndex
                SaveImage,path2,fileTomoOut,'wplot3D_tomosep',list,OutID=list.OutID,/forcereplot,/noprompt
            endif
            endcase
        endcase
        
    endif else begin

        if breplot then begin
            wset,index
            case list.ScanDim of
                1: wplot2D,list
                2: wplot3D,list
            endcase
        endif
        
        wset,index
        img=tvrd(true=1)
        if bTomo then begin
            wset,list.TomoReconIndex
            img2=tvrd(true=1)
            if SaveTomoSep eq 1 then begin
                imgS2[*,xoff2[i]:xoff2[i]+xsize2-1,yoff2[i]:yoff2[i]+ysize2-1]=img2
            endif else begin
                if NoBoard eq 0 then img2=transpose([[transpose(img2,[1,2,0])],[bytarr(xsize2,1,3)]],[2,0,1]) ; Separation boarder
                
                s=dimsize(img,3)
                if s[1] ne xsize then img=congrid(img,s[0],xsize,s[2])
                s=dimsize(img2,3)
                if s[1] ne xsize then img2=congrid(img,s[0],xsize,s[2])
        
                img=transpose([[transpose(img2,[1,2,0])],[transpose(img,[1,2,0])]],[2,0,1])
            endelse
        endif
        if SaveSep eq 1 then begin
            temp=path+fileOut
            SaveRGBImage,temp,img,OutID=list.OutID
        endif else begin
            imgS[*,xoff[i]:xoff[i]+xsize-1,yoff[i]:yoff[i]+ysize-1]=img
        endelse
    endelse
    
    if bhtmlwrapper then begin
        xi=tableindex mod list.nhorizontal
        yi=tableindex/list.nhorizontal*ymag
        tmp=string(((*list.pmstruc).(list.NSelect+1)).data[list.labeli,*])
        for k=1,n_elements(tmp)-1 do tmp[0]+='; '+tmp[k]
        table[xi,yi]=tmp[0]
        table[xi,yi+1]='<img src=".'+path_sep()+fileOut+'" alt="'+path_sep()+fileOut+'">'
        if n_elements(fileTomoOut) ne 0 then table[xi,yi+2]='<img src=".'+path_sep()+fileTomoOut+'" alt="'+path_sep()+fileTomoOut+'">'

        tableindex++
    endif
endfor

for i=0l,nempty-1 do begin
    j=i+NMaps
    if (SaveSep ne 1) then begin
        if bTomo and (SaveTomoSep eq 1) then imgS2[*,xoff2[j]:xoff2[j]+xsize2-1,yoff2[j]:yoff2[j]+ysize2-1]=255b
        imgS[*,xoff[j]:xoff[j]+xsize-1,yoff[j]:yoff[j]+ysize-1]=255b
    endif
endfor

if breplot then begin
    wdelete,index
    if bTomo then wdelete,list.TomoReconIndex
endif

if SaveSep eq 0 then begin
    if bTomo and (SaveTomoSep eq 1) then begin
        temp=path+fileTomo+ext
        SaveRGBImage,temp,imgS2,OutID=list.OutID
    endif
    temp=path+file+ext
    SaveRGBImage,temp,imgS,OutID=list.OutID
endif

if bhtmlwrapper then begin
    if ~openw_safe(lun,htmlfile) then begin
        printf,lun,'<html>'
        printf,lun,'<head>'
        printf,lun,'<title>'
        printf,lun,file
        printf,lun,'</title>'
        printf,lun,'</head>'
        printf,lun,'<body>'
        printf,lun,'Results from '+file+'.xdi'+'<br>'
        printf,lun,'<table border="1">'
        n=string(list.nhorizontal,format='(I0)')
        table=string(table,format='('+n+'("<td>",A,"</td>"))')
        printf,lun,table,format='("<tr>",A,"</tr>")'
        printf,lun,'</table>'
        printf,lun,'</body>'
        printf,lun,'</html>'
        free_lun,lun
        printw,list.OutID,'Saved '+htmlfile
    endif else printw,list.OutID,'NOT saved '+htmlfile
endif
end;pro SaveAllXDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RefreshWindowXDI,ev
widget_control,ev.top,get_uvalue=list

DDbool=list.ScanDim eq 2

; ----Set the correct sensitivity----
ID=widget_info(ev.top,FIND_BY_UNAME='magnification')
    widget_control,ID,sensitive=DDbool
ID=widget_info(ev.top,FIND_BY_UNAME='Show 3D')
    widget_control,ID,sensitive=DDbool
ID=widget_info(ev.top,FIND_BY_UNAME='3DROI')
    widget_control,ID,sensitive=DDbool
ID=widget_info(ev.top,FIND_BY_UNAME='Browse')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='Open...')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='Save...')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='Save As...')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='Save image...')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='Load...')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='Normalization')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='Options')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='basemode')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='BaseShow')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='BlendBase')
    widget_control,ID,sensitive=list.ScanDim eq 2
ID=widget_info(ev.top,FIND_BY_UNAME='InitBase')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='Backlash base')
    widget_control,ID,sensitive=1
ID=widget_info(ev.top,FIND_BY_UNAME='REV1')
       widget_control,ID,sensitive=DDbool or(list.TRAN ne 1)
ID=widget_info(ev.top,FIND_BY_UNAME='REV2')
       widget_control,ID,sensitive=DDbool or(list.TRAN)
ID=widget_info(ev.top,FIND_BY_UNAME='REV3')
       widget_control,ID,sensitive=DDbool
ID=widget_info(ev.top,FIND_BY_UNAME='First vertical')
    widget_control,ID,set_button=list.TRAN
ID=widget_info(ev.top,FIND_BY_UNAME='First horizontal')
    widget_control,ID,set_button=list.TRAN ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='Right-Left')
    widget_control,ID,set_button=list.REV1
ID=widget_info(ev.top,FIND_BY_UNAME='Left-Right')
    widget_control,ID,set_button=list.REV1 ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='Bottom-Top')
    widget_control,ID,set_button=list.REV2
ID=widget_info(ev.top,FIND_BY_UNAME='Top-Bottom')
    widget_control,ID,set_button=list.REV2 ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='Jump')
    widget_control,ID,set_button=list.REV3 ne 1
ID=widget_info(ev.top,FIND_BY_UNAME='Continue')
    widget_control,ID,set_button=list.REV3
ID=widget_info(ev.top,FIND_BY_UNAME='L--R')
    if list.Hsign eq 1 then bval='L-->R' else bval='L<--R'
    widget_control,ID,set_value=bval
ID=widget_info(ev.top,FIND_BY_UNAME='T--B')
    if list.Vsign eq 1 then bval='T-->B' else bval='T<--B'
    widget_control,ID,set_value=bval
ID=widget_info(ev.top,FIND_BY_UNAME='HORIZ')
    widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN ne 1)
ID=widget_info(ev.top,FIND_BY_UNAME='VERT')
    widget_control,ID,sensitive=(list.SCANDIM eq 2)or(list.TRAN)
;ID=widget_info(ev.top,FIND_BY_UNAME='Swap Dim')
;    widget_control,ID,set_button=list.SWAP
ID=widget_info(ev.top,FIND_BY_UNAME='Swap Dim')
       widget_control,ID,sensitive=DDbool
ID=widget_info(ev.top,FIND_BY_UNAME='Tomo base')
    widget_control,ID,sensitive=DDbool
ID=widget_info(ev.top,FIND_BY_UNAME='CutBase')
    widget_control,ID,sensitive=DDbool
ID=widget_info(ev.top,FIND_BY_UNAME='AddBase')
    widget_control,ID,sensitive=1
ID=widget_info(list.top,FIND_BY_UNAME='addangle')
    widget_control,ID,set_value=stringr(list.addangle*180/!pi)
ID=widget_info(list.top,FIND_BY_UNAME='Backlash cor.')
    widget_control,ID,set_button=list.Corb[0]
ID=widget_info(list.top,FIND_BY_UNAME='Proj. center cor')
    widget_control,ID,set_button=list.Corb[1]
ID=widget_info(list.top,FIND_BY_UNAME='Zinger removal')
    widget_control,ID,set_button=list.Corb[3]
ID=widget_info(list.top,FIND_BY_UNAME='Normalize')
    widget_control,ID,set_button=list.Corb[2]
ID=widget_info(list.top,FIND_BY_UNAME='Wobble')
    widget_control,ID,set_button=list.Corb[4]
ID=widget_info(list.top,FIND_BY_UNAME='WobbleMan')
    widget_control,ID,set_button=list.Corb[5]
ID=widget_info(list.top,FIND_BY_UNAME='Backlash')
    widget_control,ID,set_value=stringr(list.Cor[0]),sensitive=list.Corb[0]
ID=widget_info(list.top,FIND_BY_UNAME='Projection center')
    widget_control,ID,set_value=stringr(list.Cor[1]),sensitive=list.Corb[1]
ID=widget_info(list.top,FIND_BY_UNAME='Zwidth')
    widget_control,ID,set_value=stringr(list.Cor[2]),sensitive=list.Corb[3]
ID=widget_info(list.top,FIND_BY_UNAME='Zthreshold')
    widget_control,ID,set_value=stringr(list.Cor[3]),sensitive=list.Corb[3]
;ID=widget_info(list.top,FIND_BY_UNAME='Wa')
;    widget_control,ID,set_value=stringr(list.Cor[4]),sensitive=list.Corb[4]
;ID=widget_info(list.top,FIND_BY_UNAME='Wb')
;    widget_control,ID,set_value=stringr(list.Cor[5]),sensitive=list.Corb[4]
;ID=widget_info(list.top,FIND_BY_UNAME='Wangle')
;    widget_control,ID,set_value=stringr(list.Cor[6]*180/!pi),sensitive=list.Corb[4]
bTomo=widget_info(list.TomoReconID,/VALID_ID) eq 1
ID=widget_info(list.top,FIND_BY_UNAME='SinCorBase')
    widget_control,ID,sensitive=bTomo
ID=widget_info(list.top,FIND_BY_UNAME='CombineNorm')
    widget_control,ID,set_button=list.CombineNorm
ID=widget_info(list.top,FIND_BY_UNAME='CombineType')
    widget_control,ID,SET_DROPLIST_SELECT=list.CombineType
ID=widget_info(list.top,FIND_BY_UNAME='CCLEAR')
    widget_control,ID,sensitive=ptr_valid(list.pautofill)
ID=widget_info(list.top,FIND_BY_UNAME='CCLEAR2')
    widget_control,ID,sensitive=ptr_valid(list.preplacep)
ID=widget_info(list.top,FIND_BY_UNAME='FBPi')
    widget_control,ID,sensitive=list.ReconType eq 0
ID=widget_info(list.top,FIND_BY_UNAME='MLEMi')
    widget_control,ID,sensitive=list.ReconType eq 1 or list.ReconType eq 3 or list.ReconType eq 4
ID=widget_info(list.top,FIND_BY_UNAME='ARTi')
    widget_control,ID,sensitive=list.ReconType ne 0
ID=widget_info(list.top,FIND_BY_UNAME='FBP')
    widget_control,ID,set_button=list.ReconType eq 0
ID=widget_info(list.top,FIND_BY_UNAME='MLEM(radon)')
    widget_control,ID,set_button=list.ReconType eq 1
ID=widget_info(list.top,FIND_BY_UNAME='ART')
    widget_control,ID,set_button=list.ReconType eq 5
ID=widget_info(list.top,FIND_BY_UNAME='SIRT')
    widget_control,ID,set_button=list.ReconType eq 7
ID=widget_info(list.top,FIND_BY_UNAME='SART')
    widget_control,ID,set_button=list.ReconType eq 6
ID=widget_info(list.top,FIND_BY_UNAME='SIRT(radon)')
    widget_control,ID,set_button=list.ReconType eq 2
ID=widget_info(list.top,FIND_BY_UNAME='MLEM')
    widget_control,ID,set_button=list.ReconType eq 3
ID=widget_info(list.top,FIND_BY_UNAME='OSEM')
    widget_control,ID,set_button=list.ReconType eq 4
ID=widget_info(list.top,FIND_BY_UNAME='Heavy weights')
    widget_control,ID,set_button=list.arthigh,sensitive=list.ReconType ge 3
ID=widget_info(list.top,FIND_BY_UNAME='L1-norm')
    widget_control,ID,set_button=list.altavg,sensitive=list.ReconType ge 5
ID=widget_info(list.top,FIND_BY_UNAME='Fixed iterations')
    widget_control,ID,set_button=list.bTomoIter
ID=widget_info(list.top,FIND_BY_UNAME='OSEMsubsetsbase')
    widget_control,ID,sensitive=list.ReconType eq 4
ID=widget_info(list.top,FIND_BY_UNAME='Recalculate weights')
    widget_control,ID,set_button=list.recalcweights,sensitive=list.ReconType ge 3
ID=widget_info(list.top,FIND_BY_UNAME='Show Convergence')
    widget_control,ID,set_button=list.bshowconvergence,sensitive=list.ReconType ne 0
ID=widget_info(list.top,FIND_BY_UNAME='tomoavg')
    widget_control,ID,set_button=list.tomoavg
ID=widget_info(list.top,FIND_BY_UNAME='tomoavgcutoff')
    widget_control,ID,editable=list.tomoavg
ID=widget_info(list.top,FIND_BY_UNAME='tomoedge')
    widget_control,ID,set_button=list.tomoedge
ID=widget_info(list.top,FIND_BY_UNAME='tomocircle')
    widget_control,ID,set_button=list.tomocircle
ID=widget_info(list.top,FIND_BY_UNAME='tomopositive')
    widget_control,ID,set_button=list.tomopositive
ID=widget_info(list.top,FIND_BY_UNAME='tomopostcorrect')
    widget_control,ID,set_button=list.tomopostcorrect
ID=widget_info(list.top,FIND_BY_UNAME='Equal weights')
    widget_control,ID,set_button=list.tomoequalweights
ID=widget_info(list.top,FIND_BY_UNAME='tomoedgecutoff')
    widget_control,ID,editable=list.tomoedge
ID=widget_info(list.top,FIND_BY_UNAME='tomocircleR')
    widget_control,ID,editable=list.tomocircle
ID=widget_info(list.top,FIND_BY_UNAME='filterlist')
    widget_control,ID,SET_DROPLIST_SELECT=list.filtertype
ID=widget_info(list.top,FIND_BY_UNAME='filtersize')
    widget_control,ID,set_value=stringr(list.filtersize)
ID=widget_info(list.top,FIND_BY_UNAME='TomoIter')
    widget_control,ID,set_value=stringr(list.bTomoIter?list.TomoIter:list.TomoIterMax)
ID=widget_info(list.top,FIND_BY_UNAME='Convergence')
    widget_control,ID,set_value=stringr(list.TomoConvergence)
ID=widget_info(list.top,FIND_BY_UNAME='ConvergenceBase')
    widget_control,ID,sensitive=~list.bTomoIter
ID=widget_info(list.top,FIND_BY_UNAME='OSEMsubsets')
    widget_control,ID,set_value=stringr(list.OSEMsubsets)
ID=widget_info(list.top,FIND_BY_UNAME='estcutoff')
    widget_control,ID,set_value=stringr(list.estcutoff)
ID=widget_info(list.top,FIND_BY_UNAME='tomoavgcutoff')
    widget_control,ID,set_value=stringr(list.tomoavgcutoff)
ID=widget_info(list.top,FIND_BY_UNAME='tomoedgecutoff')
    widget_control,ID,set_value=stringr(list.tomoedgecutoff)
ID=widget_info(list.top,FIND_BY_UNAME='tomocircleR')
    widget_control,ID,set_value=stringr(list.tomocircleR)
ID=widget_info(list.top,FIND_BY_UNAME='tomoXmult')
    widget_control,ID,set_value=stringr(list.mXSizeTomo)
ID=widget_info(list.top,FIND_BY_UNAME='tomoYmult')
    widget_control,ID,set_value=stringr(list.mYSizeTomo)
ID=widget_info(list.top,FIND_BY_UNAME='tomoXshift')
    widget_control,ID,set_value=stringr(list.mXShiftTomo)
ID=widget_info(list.top,FIND_BY_UNAME='tomoYshift')
    widget_control,ID,set_value=stringr(list.mYShiftTomo)
ID=widget_info(list.top,FIND_BY_UNAME='tomodrho')
    widget_control,ID,set_value=stringr(list.tomodrho)

ID=widget_info(list.top,FIND_BY_UNAME='Enable')
    widget_control,ID,set_button=(**list.pgradient).buse
ID=widget_info(list.top,FIND_BY_UNAME='BlendIndex')
    widget_control,ID,set_value=stringr((**list.pgradient).index)
ID=widget_info(list.top,FIND_BY_UNAME='GradientFill')
    widget_control,ID,set_value=list.gradientfill
ID=widget_info(list.top,FIND_BY_UNAME='ngradients')
    ngradients=n_elements(*list.pgradients)
    widget_control,ID,set_value=stringr(ngradients)
ID=widget_info(list.top,FIND_BY_UNAME='ngradientslist')
    ind=where(*list.pgradients eq *list.pgradient,ct)
    names=string(lindgen(ngradients)+1,format='(I30)')
    for i=0l,ngradients-1 do begin
        group=(*(*list.pgradients)[i]).index
        if group ge list.NMaps or group lt 0 then continue
        names[i]=(*list.pmstruc).(group+1).data[0]+' '+(*list.pmstruc).(group+1).data[1]
    endfor
    widget_control,ID,set_value=names,SET_DROPLIST_SELECT=ind[0]
ID=widget_info(list.top,FIND_BY_UNAME='modes')
    widget_control,ID,set_value=string(blendmodes(),format='(A25)'),SET_DROPLIST_SELECT=(**list.pgradient).mode

; ----Print correct group's mask numbers----
ID=widget_info(ev.top,FIND_BY_UNAME='MaskTable')
i=list.NSelect+1
labeli=((*list.pmstruc).(i)).labeli
value=((*list.pmstruc).(i)).data
if widget_info(ID,/VALID_ID) eq 0 then begin
    ID=widget_info(ev.top,FIND_BY_UNAME='Mask')

    table = WIDGET_TABLE(ID, value = value,$
                            COLUMN_LABELS=(*list.pmfields)[*,labeli],$
                             X_SCROLL_SIZE = (list.npmfields)<4, $
                             Y_SCROLL_SIZE = ((*list.pmstruc).(i)).n<4,uname='MaskTable',uvalue='MaskTable',/editable)
endif else begin
    widget_control,ID,set_value=value,TABLE_YSIZE=((*list.pmstruc).(i)).n,COLUMN_LABELS=(*list.pmfields)[*,labeli]
endelse

end;pro RefreshWindowXDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro modeXDI, ev, val
widget_control,ev.top,get_uvalue=list
strmat=['Draw Region','Draw 1D Region','Make LOF','Pattern','Insert XDI','Draw Sinus','Check ddist','Extract 1D','Select Pixels','Move Legend']
ind=where(strmat ne val,count)
if count ne 0 then strmat=strmat[ind]

widget_control,list.draw,draw_motion_events=0,DRAW_KEYBOARD_EVENTS=0;,draw_button=0 ; setting draw_button=0 causes problems
if widget_info(list.TomoReconID,/VALID_ID) eq 1 then widget_control,list.TomoReconID,draw_motion_event=1,event_pro='DEFAULTTOMO_PROCESS_EVENTS'

case val of
   'Dummy':begin
             widget_control,list.draw,draw_motion_event=1,$
              event_pro='DEFAULT_PROCESS_EVENTS'
             setsens,ev.top,strmat,sen=sen
             return
             endcase
   'Draw Region':    begin
             if list.ScanDim ne 2 then begin
                 widget_control,list.draw,draw_motion_event=1,$
              event_pro='DEFAULT_PROCESS_EVENTS'
                 setsens,ev.top,strmat,sen=sen
                 return
             endif
          widget_control,list.draw,draw_motion_event=1,$
              draw_button_event=1,event_pro='ROI_PROCESS_EVENTS'
          RefreshDisplayXDI,list
          endcase
    'Select Pixels':    begin
          XDIdrawfree,list.draw
          widget_control,list.draw,draw_motion_event=1,$
              draw_button_event=1,event_pro='LOF_PROCESS_EVENTS',set_uvalue={type:2}
          endcase
    'Make LOF':    begin
          if PTR_VALID(list.filenames) eq 0 then return
          XDIdrawfree,list.draw
          widget_control,list.draw,draw_motion_event=1,$
              draw_button_event=1,event_pro='LOF_PROCESS_EVENTS',set_uvalue={type:1}
          endcase
    'Pattern': begin
          if PTR_VALID(list.filenames) eq 0 then return
          widget_control,list.draw,draw_motion_event=1,$
              draw_button_event=1,event_pro='UPDATE_PROCESS_EVENTS'
          if widget_info(list.TomoReconID,/VALID_ID) eq 1 then $
              widget_control,list.TomoReconID,draw_motion_event=1,$
                  draw_button_event=1,event_pro='UPDATE_PROCESS_EVENTS2'
          endcase
    'Insert XDI': begin
          ; ----Read xdi----
          path2=''
          path2=DIALOG_PICKFILE(path=list.path,file=list.file,filter='*.xdi',title='Read XDI...')
          if path2 EQ '' then return
          error=CutPath(path2,path=path,file=file,ext=ext)
          file=file+ext
          if list.bopenformat then format=list.openformat
          struct=ReadXDI(path,file,list.OutID,error=error,format=format)
          if error then return

          ; ----Make same orientation----
          data=struct.data
          o=struct.orientation
          TRAN=o[0]
          REV1=o[1]
          REV2=o[2]
          REV3=o[3]
          filenames=struct.filenames
          ReformXDIdatablock,data,[TRAN,list.TRAN],[REV1,list.REV1],[REV2,list.REV2],[REV3,list.REV3],0b
          s=size(data,/dim)
          Nmaps=s[0]
          MapCol=s[1]
          MapRow=s[2]
          if Nmaps ne list.Nmaps then begin
                printw,list.OutID,'XDI has wrong dimensions.'
                return
          endif
          
          bfiles=size(filenames,/type) eq 7 and ptr_valid(list.filenames)
          if bfiles then $
              ReformXDIdatablock,filenames,[TRAN,list.TRAN],[REV1,list.REV1],[REV2,list.REV2],[REV3,list.REV3],0b
    
          ; ----Make data structure and start mouse event handler----
          ID=widget_info(ev.top,FIND_BY_UNAME='Sticky xboarders for insert XDI')
          widget_control,ID,get_uvalue=stickyx
          ID=widget_info(ev.top,FIND_BY_UNAME='Sticky yboarders for insert XDI')
          widget_control,ID,get_uvalue=stickyy
          
          mult=1.
          multperc=0.1
          hexpand=0l
          vexpand=0l
          orientation=[list.TRAN,list.REV1,list.REV2,list.REV3]
          
          XDIdrawfree,list.draw
          widget_control,list.draw,set_uvalue={data:ptr_new(data),$
          bfiles:bfiles,$
          filenames:ptr_new(filenames),$
          stickyx:stickyx,$
          stickyy:stickyy,$
          ScanDim:struct.Scandim,$
          mult:mult,$
          multperc:multperc,$
          curx:MapCol/2,$
          cury:MapRow/2,$
          orientation:orientation,$
          hexpand:hexpand,$
          vexpand:vexpand,$
          expandfast:0b,$
          MapCol:MapCol,$
          MapRow:MapRow}

          (*list.data)=(*list.dataorg)
          printw,list.OutID,''
          printw,list.OutID,['Insert keyboard events:','Arrow keys (w/wo SHIFT) to move cursor offset',$
          'Change multiplier with * and / keys','Change multiplier step with + and - keys',$
          'Orientation h=Moveh,v=Movev,l=LR,r=RL,t=TB,b=BT,j=JUMP,c=CONTINU',$
          'Shrink/expand size with 4,6,2,8 (5 for faster/slower expanding)',$
          'Browse groups: n(next), p(previous)']
          printw,list.OutID,'Insert image multiplier: '+string(mult)
          printw,list.OutID,'Insert image multiplier-step: '+string(multperc)
          printw,list.OutID,'Insert image size expansion (hor x ver):'+string(hexpand)+string(vexpand)
          printw,list.OutID,'Insert image orientation: '+(orientation[0]?'Movev; ':'Moveh; ')+$
                (orientation[1]?'RL; ':'LR; ')+$
                (orientation[2]?'BT; ':'TB; ')+$
                (orientation[3]?'Continue; ':'Jump; ')
          printw,list.OutID,'Insert cursor offset: '+string(MapCol/2)+string(MapRow/2)
          widget_control,list.draw,draw_motion_event=1,$
              draw_button_event=1,DRAW_KEYBOARD_EVENTS=1,event_pro='INSERT_PROCESS_EVENTS',/INPUT_FOCUS
          
          RefreshDisplayXDI,list
          endcase
    'Draw Sinus':begin
          if widget_info(list.TomoReconID,/VALID_ID) ne 1 then return
          widget_control,list.TomoReconID,draw_motion_event=1,$
              draw_button_event=1,event_pro='DRAWSINUS_EVENTS'
          widget_control,list.draw,draw_motion_event=1,$
              event_pro='DEFAULT_PROCESS_EVENTS'
          endcase
    'Move Legend':begin
          widget_control,list.draw,draw_motion_event=1,$
              draw_button_event=1,event_pro='XDILEGEND_EVENTS'
          endcase
   'Check ddist':begin
          if widget_info(list.TomoReconID,/VALID_ID) ne 1 then return
          widget_control,list.TomoReconID,draw_motion_event=1,$
              draw_button_event=1,event_pro='ddistCheck_event'
          widget_control,list.draw,draw_motion_event=1,$
              event_pro='DEFAULT_PROCESS_EVENTS'
          endcase
   'Extract 1D': begin
             if list.ScanDim ne 2 then return
             XDIdrawfree,list.draw
          widget_control,list.draw,draw_motion_event=1,set_uvalue=1,$
              draw_button_event=1,DRAW_KEYBOARD_EVENTS=1,event_pro='EXTRACT1D_PROCESS_EVENTS'
          printw,list.OutID,['Extract keyboard events:','hold SHIFT-key: horizontal line','hold CTRL-key: vertical line','']
          endcase
   'Draw 1D Region': begin
             if list.ScanDim ne 2 then return
             XDIdrawfree,list.draw
          widget_control,list.draw,draw_motion_event=1,set_uvalue=2,$
              draw_button_event=1,event_pro='EXTRACT1D_PROCESS_EVENTS'
          endcase
    else: return
endcase
setsens,ev.top,strmat,sen=sen
if widget_info(ev.id,/valid) then widget_control,ev.id,sensitive=0
end ;pro modeXDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRD_XDI_CleanUp,ID

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
    wdelete,list.pixindexgradient
    result=UpdateINI('$XRD_XDI path:',{path:list.path,file:list.file})

    if widget_info(list.main,/valid_id) ne 0 then begin
        widget_control,list.main,get_uvalue=toplist
        toplist.maskP=-1
        toplist.xdipath=list.path
        toplist.xdifile=list.file
        widget_control,list.main,set_uvalue=toplist
        RefreshDisplay,list.main,toplist
    endif

    HEAP_FREE,list
endif
WIDGET_CONTROL, ID, /DESTROY
ControlDevice,Platform,/SHORT
set_plot, Platform

error=UnSubscribeToOutID()

print,'=========================================STOP XDI Editor'

end;XRD_XDI_CleanUp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CutInsertEvent,ev
widget_control,ev.id,get_uvalue=uval
widget_control,ev.top,get_uvalue=list
case uval of
'CutPos':begin
        widget_control,ev.id,get_value=cutpos
        ID=widget_info(ev.top,find_by_uname='Cutn')
        widget_control,ID,get_value=cutn
        endcase
'Cutn':    begin
        widget_control,ev.id,get_value=cutn
        ID=widget_info(ev.top,find_by_uname='CutPos')
        widget_control,ID,get_value=cutpos
        endcase
endcase

cutpos=strcompress(strsplit(cutpos, ',',/EXTRACT),/remove_all)
cutn=strcompress(strsplit(cutn, ',',/EXTRACT),/remove_all)
ind=where(cutpos eq '',ct1)
ind=where(cutn eq '',ct2)

ptr_free,list.cutpos
ptr_free,list.cutn

if ct1+ct2 eq 0 then begin
    n=n_elements(cutpos)<n_elements(cutn)
    cutpos=cutpos[0:n-1]
    cutn=cutn[0:n-1]
    cutpos=long(cutpos)
    cutn=long(cutn)
    ind=sort(cutpos)

    cutpos=cutpos[ind]
    cutn=cutn[ind]

    list.cutpos=ptr_new(cutpos)
    list.cutn=ptr_new(cutn)
    WorkingDataCorrect,list
endif
widget_control,ev.top,set_uvalue=list
RefreshDisplayXDI,list
end;pro CutInsertEvent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ShiftRowEvent,ev
widget_control,ev.id,get_uvalue=uval
widget_control,ev.top,get_uvalue=list
case uval of
'ShiftRownr':begin
        widget_control,ev.id,get_value=ShiftRownr
        ID=widget_info(ev.top,find_by_uname='ShiftRown')
        widget_control,ID,get_value=ShiftRown
        endcase
'ShiftRown':    begin
        widget_control,ev.id,get_value=ShiftRown
        ID=widget_info(ev.top,find_by_uname='ShiftRownr')
        widget_control,ID,get_value=ShiftRownr
        endcase
endcase

ShiftRownr=strcompress(strsplit(ShiftRownr, ',',/EXTRACT),/remove_all)
ShiftRown=strcompress(strsplit(ShiftRown, ',',/EXTRACT),/remove_all)
ind=where(ShiftRownr eq '',ct1)
ind=where(ShiftRown eq '',ct2)

ptr_free,list.rowshift

if ct1+ct2 eq 0 then begin
    n=n_elements(ShiftRownr)<n_elements(ShiftRown)
    ShiftRownr=ShiftRownr[0:n-1]
    ShiftRown=ShiftRown[0:n-1]
    
    a=0L
    b=0L
    for i=0l,n-1 do begin
        ind=strcompress(strsplit(ShiftRownr[i], ':',/EXTRACT,count=ct),/remove_all)
        if ct eq 2 then begin
            ind=long(ind)
            s=sign(ind[0])
            ind=abs(ind)
            ind=ind[sort(ind)]
            
            nind=ind[1]-ind[0]+1
            
            a=[a,s*(ind[0]+lindgen(nind))]
            b=[b,replicate(long(ShiftRown[i]),nind)]
        endif else begin
            a=[a,long(ind)]
            b=[b,long(ShiftRown[i])]
        endelse
    endfor

    list.rowshift=ptr_new([[a[1:*]],[b[1:*]]])
    WorkingDataCorrect,list
endif
widget_control,ev.top,set_uvalue=list
RefreshDisplayXDI,list
end;pro ShiftRowEvent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro bCorEvent,ev
widget_control,ev.id,get_value=val
WIDGET_CONTROL, ev.top, get_uvalue=list
case val of
'WobbleMan':    begin
            list.Corb[5]=ev.select
            if ev.select then ModelSinogram,list,rho=rho,nviews=nviews,bsinovert=bsinovert $
            else delvar2,*list.WobblePtr
            endcase
'Wobble':    begin
            list.Corb[4]=ev.select
            endcase
'Normalize':list.Corb[2]=ev.select
'Backlash cor.':list.Corb[0]=ev.select
'Proj. center cor':list.Corb[1]=ev.select
'Zinger removal':list.Corb[3]=ev.select
'Fit':      ModelSinogram,list,rho=rho,nviews=nviews,bsinovert=bsinovert
'COPManual': FindCOPManual,ev.top,list
'FBP':    begin
        if ~ev.select then return
        list.ReconType=0
        endcase
'MLEM(radon)':    begin
        if ~ev.select then return
        list.ReconType=1
        endcase
'SIRT(radon)':    begin
        if ~ev.select then return
        list.ReconType=2
        endcase
'MLEM':    begin
        if ~ev.select then return
        list.ReconType=3
        endcase
'OSEM':    begin
        if ~ev.select then return
        list.ReconType=4
        endcase
'ART':    begin
        if ~ev.select then return
        list.ReconType=5
        endcase
'SART':    begin
        if ~ev.select then return
        list.ReconType=6
        endcase
'SIRT':    begin
        if ~ev.select then return
        list.ReconType=7
        endcase
'Area intersect':begin
        keep=list.arthigh
        list.arthigh=ev.select
        if list.arthigh then begin
            no=dialog_message('This is just for testing and might take very long. Continue?',/question) eq 'No'
            if no then begin
                list.arthigh=0b
                widget_control,ev.id,set_button=0
            endif
        endif
        if keep ne list.arthigh then list.recalcweights=1b
        endcase
'Recalculate weights':begin
        list.recalcweights=ev.select
        endcase
'Show Convergence':begin
        list.bshowconvergence=ev.select
        endcase
'L1-norm':begin
        list.altavg=ev.select
        list.recalcweights=1b
        endcase
'Fixed iterations': begin
        list.bTomoIter=ev.select
        endcase
'Ray average (no sum)': list.tomoavg=ev.select
'Remove edge': list.tomoedge=ev.select
'Clip circle': list.tomocircle=ev.select
'Positive': list.tomopositive=ev.select
'Constraints after reconstruction': list.tomopostcorrect=ev.select
'Binarize weights': list.tomoequalweights=ev.select
endcase

WorkingDataCorrect,list
WIDGET_CONTROL, ev.top, set_uvalue=list
RefreshDisplayXDI,list
if n_elements(rho) ne 0 then begin
    wset,list.drawindex
    DRAWSINUS,rho,nviews,bsinovert,list
endif
RefreshWindowXDI,ev
end;pro bCorEvent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CorEvent,ev
widget_control,ev.id,get_uvalue=uval
widget_control,ev.top,get_uvalue=list
widget_control,ev.id,get_value=str
case uval of
'Backlash':list.Cor[0]=float(str)
'Hbacklash':list.ManualBacklash[0]=float(str)
'Vbacklash':list.ManualBacklash[1]=float(str)
'TomoIter':    begin
            if list.bTomoIter then list.TomoIter=0>fix(str) else list.TomoIterMax=0>fix(str)
            endcase
'Convergence':list.TomoConvergence=0>float(str)<100
'OSEMsubsets':list.OSEMsubsets=1>fix(str)
'estcutoff':list.estcutoff=float(str)
'tomoavgcutoff':list.tomoavgcutoff=float(str)
'tomoedgecutoff':list.tomoedgecutoff=fix(str)>1
'tomocircleR':list.tomocircleR=float(str)
'tomoXmult':begin
            list.mXSizeTomo=float(str)
            xdi_changesize,list
            endcase
'tomoYmult':begin
            list.mYSizeTomo=float(str)
            xdi_changesize,list
            endcase
'tomoXshift':list.mXShiftTomo=float(str)
'tomoYshift':list.mYShiftTomo=float(str)
'tomodrho':begin
          list.tomodrho=float(str)
          xdi_changesize,list
          endcase
'filtersize':begin
        val=0>round(float(str))<((list.MapCol-2)/2)
        widget_control,ev.id,set_value=stringr(val)
        list.filtersize=val
        endcase
'addangle':begin
        val=0.
        reads,str,val
        widget_control,ev.id,set_value=stringr(val)
        list.addangle=val*!pi/180
        endcase
'Projection center':list.Cor[1]=float(str)
'Zwidth':list.Cor[2]=float(str)
'Wa':    list.Cor[4]=float(str)
'Wb':    list.Cor[5]=float(str)
'Wangle':list.Cor[6]=float(str)/180.*!pi
'Zthreshold':list.Cor[3]=float(str)
endcase

WorkingDataCorrect,list
widget_control,ev.top,set_uvalue=list
RefreshDisplayXDI,list
RefreshWindowXDI,ev
end;pro CorEvent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AutoFillEvent,ev
widget_control,ev.id,get_value=val
widget_control,ev.top,get_uvalue=list
if list.ScanDim ne 2 then return

ID=widget_info(ev.top,find_by_uname='Cutoff')
widget_control,ID,get_value=cutoff

; Get operator
cutoff=byte(cutoff)
ind1=where(((cutoff ge 60)and(cutoff le 62))or(cutoff eq 33) ,ct1,COMPLEMENT=ind2,NCOMPLEMENT=ct2)
if ct1 eq 0 then operator='<=' else operator=string(cutoff[ind1])

; Groups to compare
indgroup=list.NSelect
bprocess=0
if ct2 eq 0 then cutoff=0 else begin
    cutoff=string(cutoff[ind2])
    tmp=strpos(cutoff,'g')
    cutoff=strsplit(cutoff,'*g',/extract,count=ct)
    case ct of
    1:     begin
        cutoff=float(cutoff[0])
        if tmp ne -1 then begin
            ; Compare groups 0-1 and 2-3 and ...
            indgroup=lindgen(list.NMaps)
            cutofforg=cutoff
            
            bprocess=1
        endif
        endcase
    2:     begin
        ; Compare one group to another and apply the changes on all groups
        ind=long(cutoff[1])
        if ind ge list.NMAPS then return
        cutoff=float(cutoff[0])
        
        cutoff*=(*list.data)[ind,*,*]
        ReformXDIdatablock,cutoff,[list.TRAN,0b],[list.REV1,0b],[list.REV2,0b],[list.REV3,0b],0b,s0,s1,s2
        cutoff=reform(cutoff,s1,s2,/overwrite)
        
        endcase
    else:return
    endcase
endelse

case val of
'DOWN':direction=0
'UP':direction=1
'LEFT':direction=2
'RIGHT':direction=3
'REPLACE':direction=-1
'LOF':direction=-2
endcase

if ~ptr_valid(list.pautofill) then list.pautofill=ptr_new(ptrarr(list.NMaps,/alloc)) else $
if n_elements(*list.pautofill) ne list.NMaps then begin
    heap_free,*list.pautofill
    *list.pautofill=ptrarr(list.NMaps,/alloc)
endif

if ~ptr_valid(list.preplacep) then list.preplacep=ptr_new(ptrarr(list.NMaps,/alloc)) else $
if n_elements(*list.preplacep) ne list.NMaps then begin
    heap_free,*list.preplacep
    *list.preplacep=ptrarr(list.NMaps,/alloc)
endif

n=n_elements(indgroup)
if n gt 1 then n-=n mod 2

for i=0l,n-1 do begin
    if bprocess eq 1 and ((i mod 2) eq 1) then continue
    
    ; Default orientation
    image=(*list.data)[indgroup[i],*,*]
    delvar2,s0
    delvar2,s1
    delvar2,s2
    ReformXDIdatablock,image,[list.TRAN,0b],[list.REV1,0b],[list.REV2,0b],[list.REV3,0b],0b,s0,s1,s2
    image=reform(image,s1,s2,/overwrite)
    
    ; Group as cutoff
    if bprocess eq 1 then begin
        cutoff=cutofforg*(*list.data)[indgroup[i]+1,*,*]
        delvar2,s0
        delvar2,s1
        delvar2,s2
        ReformXDIdatablock,cutoff,[list.TRAN,0b],[list.REV1,0b],[list.REV2,0b],[list.REV3,0b],0b,s0,s1,s2
        cutoff=reform(cutoff,s1,s2,/overwrite)
    endif
    
    plist=(*list.pautofill)[indgroup[i]]
    preplacep=(*list.preplacep)[indgroup[i]]
    if bprocess eq 0 then begin
        plistall=list.pautofill
        preplacepall=list.preplacep
    endif
    
    autofindpixels,image,plist,preplacep,cutoff,direction,operator,ev.top,lofind=lofind,plistall=plistall,preplacepall=preplacepall,replacevalkeep=replacevalkeep
endfor

if direction eq -2  and n_elements(lofind) gt 0 then begin
    header=strarr(4)
    header[0]='$# columns:'
    header[1]=stringr(n_elements(lofind))
    header[2]='$# rows:'
    header[3]='1'

    error=CutPath(list.file,file=file)
    filenames=reform(*list.filenames)
    tmp=where(filenames eq '',ct)
    if ct ne 0 then filenames[tmp]='c:\dummy.tif'
    Result = DIALOG_MESSAGE('Save with corrections?',/question)
    if Result eq 'Yes' then PixelJuggle,filenames,list,'c:\dummy.tif'

    error=WriteLOF(list.path,file+'.lof',lofind,ev.top,filenames,list.OutID,header=header)
    if error eq 0 then printw,list.OutID,'Couldn`t save .lof file!!!'
endif else begin
    WorkingDataCorrect,list
    widget_control,ev.top,set_uvalue=list
    RefreshDisplayXDI,list
    RefreshWindowXDI,ev
endelse
end;pro AutoFillEvent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRD_XDI_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
2:     begin; ----Event from widget_slider----
    WIDGET_CONTROL, ev.id, GET_VALUE=val
    WIDGET_CONTROL, ev.top, get_uvalue=list
    WIDGET_CONTROL, ev.id, get_uvalue=uval
    
    case uval of
    'scale':begin
           list.scaleimg=val
           RefreshDisplayXDI,list
           WIDGET_CONTROL, ev.top, set_uvalue=list
           endcase
    'NRslider':begin
           list.NSelect=0>val<(list.NMaps-1)
           RefreshDisplayXDI,list
           WIDGET_CONTROL, ev.top, set_uvalue=list
           SelectGroupMain,ev
           RefreshWindowXDI,ev
           endcase
    'Highslider':begin
           if list.nrROI eq 0 then return
           ID=widget_info(ev.top,FIND_BY_UNAME='Names')
           (*list.highlight)[*]=0
           if val eq list.nrROI+1 then begin
             (*list.highlight)[*]=1
             widget_control,ID,set_value='Selected all'
           endif else $
               if val ne 0 then begin
                 (*list.highlight)[val-1]=1B
                 widget_control,ID,set_value=(*list.names)[val-1]
               endif else widget_control,ID,set_value='No selection'
    
           RefreshDisplayXDI,list
           WIDGET_CONTROL, ev.top, set_uvalue=list
           endcase
    'GradientFill':begin
            widget_control,ev.top,get_uvalue=list
            list.gradientfill=val
            widget_control,ev.top,set_uvalue=list
            RefreshGradientImage,ev
            RefreshWindowXDI,ev
            endcase
    endcase
    endcase

1:    begin; ----Event from widget_button----
    WIDGET_CONTROL, ev.id, GET_VALUE=val
    case val of
        'Open XDI':begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              if list.bopenformat then format=list.openformat
              OpenXDI,ev,1,format=format
              endcase
        'Default image output...': begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              if list.bopenformat then format=list.openformat
              t=CutPath(list.file,file=file)
              files=MDIALOG_PICKFILE(path=list.path,file=file+'.xdi',filter='*.xdi',$
                            title='Select files to process...',/MULTIPLE_FILES)
              if files[0] eq '' then return
              for i=0l,n_elements(files)-1 do begin
                  OpenXDI,ev,0,format=format,path=files[i]
                  SaveAllXDI,ev,/noprompt
              endfor
              endcase
        'Add XDI':begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              if list.bopenformat then format=list.openformat
              OpenXDI,ev,1,/add,format=format
              endcase
        'Show 3D': begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              case list.ScanDim of
              2:  begin
                  data=reform(((*list.data))[list.NSelect,*,*],list.MapCol,list.MapRow)
                  list2={data:data,xr:[0,list.MapCol-1],yr:[0,list.MapRow-1],zr:[min(data),max(data)],$
                     proplot:'DDXDIplot',xtitle:'X(column)',ytitle:'Y(row)',ztitle:'Integrated Intensity',style:2}
                  flplot_obj,list2,GROUP_LEADER=ev.top,/xtrue,/ytrue
                  endcase
              1:  begin
                    data=transpose(reform(*list.data,list.MapCol,list.MapRow))
                    getXDIdatadim,list,nmaps,ncol,nrow
                    list2={data:data,xr:[0,ncol*nrow-1],yr:[0,nmaps-1],zr:[min(data),max(data)],$
                     proplot:'DDXDIplot',xtitle:'X(files)',ytitle:'Y(group)',ztitle:'Integrated Intensity',style:3}
                  flplot_obj,list2,GROUP_LEADER=ev.top
                    endcase
              endcase
              endcase
        'Save RGB triangle': RGBtriangle,ev
        'Exit':       begin
              widget_control,ev.top,/destroy
              endcase
        'Fill':       begin
              widget_control,ev.id,set_uvalue=ev.select
              WIDGET_CONTROL, ev.top, get_uvalue=list
              RefreshDisplayXDI,list
              endcase
        'Continuous':       begin
              widget_control,ev.id,set_uvalue=ev.select
              WIDGET_CONTROL, ev.top, get_uvalue=list
              RefreshDisplayXDI,list
              endcase
        'Log Scale':begin
              widget_control,ev.id,set_uvalue=ev.select
              WIDGET_CONTROL, ev.top, get_uvalue=list
              RefreshDisplayXDI,list
              endcase
           'Binned pixels':begin
              widget_control,ev.id,set_uvalue=ev.select
              WIDGET_CONTROL, ev.top, get_uvalue=list
              RefreshDisplayXDI,list
              endcase
        'Add contour':begin
              widget_control,ev.id,set_uvalue=ev.select
              WIDGET_CONTROL, ev.top, get_uvalue=list
              RefreshDisplayXDI,list
              endcase
        'Save Groups Separate':begin
              widget_control,ev.id,set_uvalue=ev.select
              return
              endcase
        'Save Tomo Separate':begin
              widget_control,ev.id,set_uvalue=ev.select
              return
              endcase
        'Save without boarders':begin
              widget_control,ev.id,set_uvalue=ev.select
              return
              endcase
        'Save 1D data in 1 plot':begin
              widget_control,ev.id,set_uvalue=ev.select
              return
              endcase
        'Save with labels':begin
              widget_control,ev.id,set_uvalue=ev.select
              widget_control,ev.top,get_uvalue=list
              RefreshDisplayXDI,list
              return
              endcase
        'Update Main':begin
              widget_control,ev.id,set_uvalue=ev.select
              return
              endcase
           'Sticky xboarders for insert XDI':begin
              widget_control,ev.id,set_uvalue=ev.select
              return
              endcase
           'Sticky yboarders for insert XDI':begin
              widget_control,ev.id,set_uvalue=ev.select
              return
              endcase
           'Dynamic reconstruct (slow)':begin
                 WIDGET_CONTROL, ev.top, get_uvalue=list
                 list.dynrecon=ev.select
              WIDGET_CONTROL, ev.top, set_uvalue=list
                 endcase
        'Update sinogram':begin
                 WIDGET_CONTROL, ev.top, get_uvalue=list
                 list.updatesinogram=ev.select
              WIDGET_CONTROL, ev.top, set_uvalue=list
              RefreshDisplayXDI,list
                 endcase
           'Tomo Reconstruct':begin
                 widget_control,ev.id,set_uvalue=ev.select
              WIDGET_CONTROL, ev.top, get_uvalue=list
              xdi_changesize,list
              WIDGET_CONTROL, ev.top, set_uvalue=list ;???
              WorkingDataCorrect,list
              widget_control,ev.top,set_uvalue=list
              RefreshDisplayXDI,list
              ID=widget_info(list.top,FIND_BY_UNAME='SinCorBase')
                widget_control,ID,sensitive=ev.select
                 endcase
        'Suspend update':begin
              widget_control,ev.id,set_uvalue=ev.select
              WIDGET_CONTROL, ev.top, get_uvalue=list
              RefreshDisplayXDI,list
              endcase
        'Close Region':begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              if closeregion(list) then RefreshDisplayXDI,list
              endcase
        'Fit Gaussians':begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              GaussFitROIs,list,type=0
              WIDGET_CONTROL, ev.top, set_uvalue=list
              endcase
        'Fit Steps':begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              ;SaveImage,'c:\','test.eps','GaussFitROIs',list,type=1,OutID=list.OutID
              
              GaussFitROIs,list,type=1
              WIDGET_CONTROL, ev.top, set_uvalue=list
              endcase
        'Wall thickness':begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              GaussFitROIs,list,type=2
              WIDGET_CONTROL, ev.top, set_uvalue=list
              endcase
        'Save RGN':    begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              if list.pverts eq -1 then return
              file = strmid(list.file,0,STRLEN(list.file)-4)
              path=SelFile(list.path,file+'.rgn','*.rgn','Save Regions To....')
              if path eq '' then return
              if openw_safe(lun, path) then begin
                  printw,list.OutID,path+' not saved'
                  return
              endif
    
              printf,lun,list.pverts
              printf,lun,list.nrRoi
              if list.nrRoi ne 0 then printf,lun,*list.nrPoints
              for i=0l,(list.nrRoi-1)>0 do printf,lun,(*list.names)[i]
              if ptr_valid(list.xverts) then printf,lun,(*list.xverts)
              if ptr_valid(list.yverts) then printf,lun,(*list.yverts)
    
              free_lun,lun
              printw,list.OutID,'Regions saved to: '+path
              endcase
        'Load RGN':    begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              file = strmid(list.file,0,STRLEN(list.file)-4)
              path=''
              path=DIALOG_PICKFILE(path=list.path,file=file+'.rgn',filter='*.rgn',title='Load Regions...')
              if openr_safe(lun, path) then begin
                  printw,list.OutID,'.rgn file not found'
                  return
              endif
              ; ----Clear memory----
              ptr_free,list.nrPoints
              ptr_free,list.xverts
              ptr_free,list.yverts
              ptr_free,list.highlight
              ptr_free,list.names
              ; ----Read data----
              in=0L
              readf,lun,in
              list.pverts=in
    
              in=0
              readf,lun,in
              list.nrRoi=in
    
              if list.nrRoi ne 0 then begin
                  in=lonarr(list.nrRoi)
                  readf,lun,in
                  list.nrPoints=ptr_new(in)
    
                  in=''
                  names=strarr(list.nrRoi)
                  for i=0l,(list.nrRoi-1) do begin
                  readf,lun,in
                  names[i]=in
                  endfor
                  list.names=ptr_new(names)
    
                  list.highlight=ptr_new(bytarr(list.nrROI))
              endif
    
              if list.pverts ne -1 then begin
                  in=fltarr(list.pverts+1)
                  readf,lun,in
                  list.xverts=ptr_new(in)
    
                  readf,lun,in
                  list.yverts=ptr_new(in)
              endif
    
              ; ----Finish----
              free_lun,lun
              printw,list.OutID,'Mask loaded from: '+path
              ID=widget_info(ev.top,FIND_BY_UNAME='Highslider')
              widget_control,ID,set_value=0,set_slider_max=list.nrROI+1
              WIDGET_CONTROL, ev.top, set_uvalue=list
              RefreshDisplayXDI,list
              endcase
        'Save Group':    begin
              SaveAllXDI,ev,/current
              endcase
        'Save All':begin
              SaveAllXDI,ev
              endcase
        'Save colorbar':begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              if list.ScanDim eq 2 then begin
                  window,/free;,/pixmap
                  erase
                  LoadctX,list
                  keep=!d.window
                  savecolorbar,list
                  saveimage,list.path,'colorbar.eps','savecolorbar',list
                  wdelete,keep
              endif
              endcase
        'Draw Region':    begin
              modeXDI, ev, val
              endcase
        'Draw 1D Region': modeXDI, ev, val
        'Make LOF':    begin
              modeXDI, ev, val
              endcase
        'Select Pixels': modeXDI, ev, val
        'Pattern': begin
              modeXDI, ev, val
              endcase
        'Delete Select':begin
              WIDGET_CONTROL, ev.top, get_uvalue=list
              if list.nrROI eq 0 then return
              ind=where((*list.highlight) eq 1,count)
              ; count can be 0, 1 or everything
              if count eq 0 then return
              if count ne 1 then begin
                  ; ----Clear memory----
                  list.nrROI=0
                  ptr_free,list.nrPoints
                  list.pverts=-1
                  ptr_free,list.xverts
                  ptr_free,list.yverts
                  ptr_free,list.highlight
                  ptr_free,list.names
                  ID=widget_info(ev.top,FIND_BY_UNAME='Highslider')
                  widget_control,ID,set_value=0,set_slider_max=1
                  ID=widget_info(ev.top,FIND_BY_UNAME='Names')
                  widget_control,ID,set_value=''
                  WIDGET_CONTROL, ev.top, set_uvalue=list
                  RefreshDisplayXDI,list
                  return
              endif
              ; ----Collect region data----
              delnr=(*list.nrPoints)[ind[0]]
              if ind[0] ne 0 then indstart=total((*list.nrPoints)[0:ind[0]-1]) else indstart=0
              ind2=lindgen(delnr)+indstart
              ; ----Reduce Arrays----
              list.nrROI-=count
              list.pverts-=delnr
              (*list.nrPoints)=ShrinkArray(*list.nrPoints,ind)
              (*list.names)=ShrinkArray(*list.names,ind)
              (*list.xverts)=ShrinkArray(*list.xverts,ind2)
              (*list.yverts)=ShrinkArray(*list.yverts,ind2)
              (*list.highlight)=ShrinkArray(*list.highlight,ind)
    
              ID=widget_info(ev.top,FIND_BY_UNAME='Highslider')
              widget_control,ID,set_value=0,set_slider_max=list.nrROI+1
              ID=widget_info(ev.top,FIND_BY_UNAME='Names')
              widget_control,ID,set_value=''
              WIDGET_CONTROL, ev.top, set_uvalue=list
              RefreshDisplayXDI,list
              endcase
    'First horizontal':  begin
              ; Movev will create event
              endcase
    'First vertical':  begin
              ConvertDataBlock,ev,32*ev.select or 4
              RefreshWindowXDI,ev
              endcase
    'Left-Right':     begin
              ; RL will create event
              endcase
    'Right-Left':     begin
              ConvertDataBlock,ev,32*ev.select or 2
              endcase
    'Top-Bottom':     begin
              ; BU will create event
              endcase
    'Bottom-Top':     begin
              ConvertDataBlock,ev,32*ev.select or 1
              endcase
    'Jump':      begin
              ; CONTINU will create event
              endcase
    'Continue':begin
              ConvertDataBlock,ev,32*ev.select or 8
              endcase
    'L-->R':       begin
                WIDGET_CONTROL, ev.top, get_uvalue=list
                list.Hsign=-1
                WIDGET_CONTROL, ev.top, set_uvalue=list
                widget_control,ev.id,set_value='L<--R'
                RefreshDisplayXDI,list
                endcase
    'L<--R':       begin
                WIDGET_CONTROL, ev.top, get_uvalue=list
                list.Hsign=1
                WIDGET_CONTROL, ev.top, set_uvalue=list
                widget_control,ev.id,set_value='L-->R'
                RefreshDisplayXDI,list
                endcase
    'T-->B':       begin
                WIDGET_CONTROL, ev.top, get_uvalue=list    
                list.Vsign=-1
                WIDGET_CONTROL, ev.top, set_uvalue=list
                widget_control,ev.id,set_value='T<--B'
                RefreshDisplayXDI,list
                endcase
    'T<--B':       begin
                WIDGET_CONTROL, ev.top, get_uvalue=list    
                list.Vsign=1
                WIDGET_CONTROL, ev.top, set_uvalue=list
                widget_control,ev.id,set_value='T-->B'
                RefreshDisplayXDI,list
                endcase
    'Resize':begin
              EnlargeDataBlock,ev,0
                endcase
    'Reform':begin
              EnlargeDataBlock,ev,1
                endcase
    'Swap Dim':begin
              ConvertDataBlock,ev,48
              endcase
    'Make LOF from Int': LOF_From_PixelInt,ev
    'Change coord.': ChangeXDIcoord,ev
    'Save XDI': UpdateXDI,ev
    'Save Tomo': UpdateXDI,ev,/tomo
    'Save CHI':XDI2CHI,ev
    'Save CHI/TIFF':XDI2CHIb,ev
    'Wobble': bCorEvent,ev
    'WobbleMan': bCorEvent,ev
    'Normalize':bCorEvent,ev
    'Backlash cor.':bCorEvent,ev
    'Proj. center cor':bCorEvent,ev
    'Zinger removal':bCorEvent,ev
    'Fit':      bCorEvent,ev
    'COPManual': bCorEvent,ev
    'XDI to IMS': XDIconvert,ev,0
    'IMS to XDI': XDIconvert,ev,1
    'EDF to XDI': XDIconvert,ev,2
    'TIFF to XDI': XDIconvert,ev,3
    'XDI to EDF': XDIconvert,ev,4
    'XDI to TIFF': XDIconvert,ev,6
    'XDI to ATHENA': XDIconvert,ev,7
    'PyMCA DAT to XDI': XDIconvert,ev,5
    'MicroXAS FastXRF to XDI': XDIconvert,ev,8
    'LiveTime from CCD images':begin
              DCNormal,ev,6
              endcase
    'Flux from XDI':begin
              DCNormal,ev,5
              endcase
    'Flux from TIFF':begin
              DCNormal,ev,4
              endcase
    'Flux from SPEC':begin
              DCNormal,ev,3
              endcase
    'Flux from IMS':begin
              DCNormal,ev,1
              endcase
    'Flux or LiveTime from SPE':begin
              DCNormal,ev,0
              endcase
    'Flux from FIO':begin
              DCNormal,ev,2
              endcase
    'Flux from MCA':begin
              DCNormal,ev,7
              endcase
    'Flux and DT from TIFF':begin
              DCNormal,ev,8
              endcase
    'Flux from EDF':begin
              DCNormal,ev,9
              endcase
    'Insert XDI':begin
              modeXDI,ev,val
              endcase
    'Draw Sinus':begin
              modeXDI,ev,val
              endcase
    'Move Legend':begin
              modeXDI,ev,val
              endcase
    'Check ddist': modeXDI,ev,val
    'Extract 1D':begin
              modeXDI,ev,val
              endcase
    'SortOptions':MakeBOptionsWindow,ev
    'XIAParse':begin
                widget_control,ev.top,get_uvalue=list
                 ret=ParseXiast(list.path,outid=list.OutID,/prompt)
                 endcase
    'Prepare for tomorecon: SkyScan':begin
                widget_control,ev.top,get_uvalue=list
                if ~ptr_valid(list.data) then return
                sino=(*list.data)[list.NSelect,*,*]
                sino*=0.9*(2.^16-1)/max(sino)
                sino=uint(sino)
                sino=rebin(sino,10,list.MapCol,list.MapRow,/sample)
                sino=transpose(sino,[2,1,0])
                
                list2={path:list.path,file:list.file,Format:'',Platform:0,OutID:list.OutID}
                tmp=cutpath(list.file,file=filename)
                SkyScan_NRecon,ev.top,list2,data=sino,filename=filename,sortmethod=list.sortmethod,separator=list.sortseparator
                endcase
    'Parameters': XDIparam,ev
    'Combine with other XDI...':begin
              widget_control,ev.top,get_uvalue=list
              AddXDIGroups,ev
              endcase
    'normalized':begin
              widget_control,ev.top,get_uvalue=list,/no_copy
              list.CombineNorm = ev.select
              widget_control,ev.top,set_uvalue=list,/no_copy
              endcase
    'Create Overlap': AddXDIGroups,ev
    'Correlate': CorrelateXDIGroups,ev
    'Create Blend Group': AddXDIGroups,ev
    'Reverse RGB channels':
    'Choose...':begin
            color=ColorPicker(ev.top)
            if n_elements(color) eq 3 then begin
                widget_control,ev.id,get_uvalue=uval
                case uval of
                'Start': ID=widget_info(ev.top,find_by_uname='BlendColorsStart')
                'Stop': ID=widget_info(ev.top,find_by_uname='BlendColorsStop')
                'Back': ID=widget_info(ev.top,find_by_uname='BlendColorsBack')
                endcase
                if uval eq 'Back' then colors='' else widget_control,ID,get_value=colors
                colors+=string(color,format='("[",I3,",",I3,",",I3,"]")')
                widget_control,ID,set_value=colors
            endif
            endcase
    'CLEAR ':begin
            widget_control,ev.top,get_uvalue=list
            heap_free,list.preplacep
            WorkingDataCorrect,list
            RefreshDisplayXDI,list
            RefreshWindowXDI,ev
            endcase
    'CLEAR':begin
            widget_control,ev.top,get_uvalue=list
            heap_free,list.pautofill
            WorkingDataCorrect,list
            RefreshDisplayXDI,list
            RefreshWindowXDI,ev
            endcase
    'DOWN': AutoFillEvent,ev
    'UP':AutoFillEvent,ev
    'LEFT':AutoFillEvent,ev
    'RIGHT':AutoFillEvent,ev
    'REPLACE':AutoFillEvent,ev
    'LOF':AutoFillEvent,ev
    'Load Filenames':begin
            widget_control,ev.top,get_uvalue=list
            ; ----Read xdi----
            path2=''
            path2=DIALOG_PICKFILE(path=list.path,file=list.file,filter=['*.xdi','*.tiff'],title='Read XDI or batch tiff...')
            if path2 EQ '' then return
            error=CutPath(path2,path=path,file=file,ext=ext)
            file+=ext
    
            blist=1b
            case strlowcase(ext) of
            '.xdi':    begin
                    if list.bopenformat then format=list.openformat
                    struc=ReadXDI(path,file,list.OutID,error=error,format=format)
                    if error then return
                    filenames=struc.filenames
                    orientation=struc.orientation
                    blist=0b
                    endcase
            '.tiff':begin
                    error=ReadScan(path,file,ext,list.OutID)
                    if ~error then return
                    filenames=replicate(path2,list.MapCol*list.MapRow)
                    orientation=bytarr(4)
                    endcase
            else:     begin
                    filenames=Select_Files(ev.top,{path:path,file:file,Format:ext,Platform:0},$
                            outpath=path,outfile=file,sortmethod=list.sortmethod,separator=list.sortseparator)
                    orientation=bytarr(4)
                    endelse
            endcase
            
            nfiles=n_elements(filenames)
            npixels=list.MapCol*list.MapRow
            if nfiles ne npixels then begin
                if nfiles eq list.MapRow then begin
                    filenames=reform(rebins(reform(filenames,1,list.MapRow),list.MapCol,list.MapRow),npixels)
                    nfiles=npixels
                endif
            endif
    
            if nfiles eq npixels then begin
                if blist then GetXDIdatadim,list,s0,s1,s2,/default
                ReformXDIDatablock,filenames,[orientation[0],list.TRAN],[orientation[1],list.REV1],$
                                            [orientation[2],list.REV2],[orientation[3],list.REV3],0b,s0,s1,s2
                heap_free,list.filenames
                list.filenames=ptr_new(temporary(filenames))
                ind=where(*list.filenames eq '',ct)
                if ct ne 0 and (*list.filenames)[0] ne '' then (*list.filenames)[ind]=(*list.filenames)[0]
                widget_control,ev.top,set_uvalue=list
                printw,list.OutID,'Filenames loaded.'
            endif else printw,list.OutID,'Not the same dimension.'
            endcase
     'FBP':    bCorEvent,ev
     'MLEM':bCorEvent,ev
     'OSEM':bCorEvent,ev
     'MLEM(radon)':bCorEvent,ev
     'ART':bCorEvent,ev
     'SIRT(radon)':bCorEvent,ev
     'SIRT':bCorEvent,ev
     'SART':bCorEvent,ev
     'Area intersect':bCorEvent,ev
     'L1-norm':bCorEvent,ev
     'Fixed iterations':bCorEvent,ev
     'Recalculate weights':bCorEvent,ev
     'Show Convergence':bCorEvent,ev
     'Ray average (no sum)':bCorEvent,ev
     'Remove edge':bCorEvent,ev
     'Clip circle':bCorEvent,ev
     'Positive':bCorEvent,ev
     'Constraints after reconstruction':bCorEvent,ev
     'Binarize weights': bCorEvent,ev
     'Animate Tomo': AnimateTomo,ev
     'FWHM:': wplot2DgetFWHM,ev
     'Current tomogram as weights': TomogramWeightFix,ev
     'Enable':    begin
                 widget_control,ev.top,get_uvalue=list
                 (**list.pgradient).buse=ev.select
                 RefreshGradientImage,ev
                 endcase
     'Bring Forward':begin
                 widget_control,ev.top,get_uvalue=list
                 MoveGradientLayer,list,1
                 RefreshGradientImage,ev
                 RefreshWindowXDI,ev
                 endcase
     'Send Backward':begin
                 widget_control,ev.top,get_uvalue=list
                 MoveGradientLayer,list,-1
                 RefreshGradientImage,ev
                 RefreshWindowXDI,ev
                 endcase
     'Reset Gamma':begin
                 widget_control,ev.top,get_uvalue=list
                 (**list.pgradient).gammap=0.5
                 (**list.pgradient).gammapopaq=0.5
                 RefreshGradientImage,ev
                 RefreshWindowXDI,ev
                 endcase
     'Reset Opacity':begin
                 widget_control,ev.top,get_uvalue=list
                 (**list.pgradient).opaqlow=1
                 (**list.pgradient).opaqhigh=1
                 (**list.pgradient).opaqback=1
                 (**list.pgradient).opaqcliplow=0
                 (**list.pgradient).opaqcliphigh=1
                 RefreshGradientImage,ev
                 RefreshWindowXDI,ev
                 endcase
     'Reset Clipping':begin
                 widget_control,ev.top,get_uvalue=list
                 (**list.pgradient).cliplow=0.01
                 (**list.pgradient).cliphigh=1
                 RefreshGradientImage,ev
                 RefreshWindowXDI,ev
                 endcase
     'Save Blend':begin
                 widget_control,ev.top,get_uvalue=list
                SaveGradients,list
                 endcase
    'Restore Blend':begin
                widget_control,ev.top,get_uvalue=list
                RestoreGradients,list
                widget_control,ev.top,set_uvalue=list
                RefreshGradientImage,ev
                 RefreshWindowXDI,ev
                 endcase
     'Save Gradient image':begin
                 RefreshGradientImage,ev,/save
                 endcase
    endcase
    endcase

3:    begin; ----Event from widget_text----
    widget_control,ev.id,get_uvalue=uvalue
    case uvalue of
      'Magx': begin
           widget_control,ev.top,get_uvalue=list
           widget_control,ev.id,get_value=str
           list.sfh=1./float(str)
           xdi_changesize,list
    
           widget_control,ev.top,set_uvalue=list
           RefreshDisplayXDI,list
           endcase
      'Magy': begin
           widget_control,ev.top,get_uvalue=list
           widget_control,ev.id,get_value=str
           list.sfv=1./float(str)
           xdi_changesize,list
    
           widget_control,ev.top,set_uvalue=list
           RefreshDisplayXDI,list
           endcase
    'Names': begin
           widget_control,ev.top,get_uvalue=list
           if list.nrRoi eq 0 then return
           ID=widget_info(ev.top,FIND_BY_UNAME='Highslider')
           widget_control,ID,get_value=select
           if (select eq 0) or (select eq list.nrROI+1) then $
             widget_control,ev.id,set_value='' else begin
             widget_control,ev.id,get_value=str
             (*list.names)[select-1]=str
             widget_control,ev.top,set_uvalue=list
             endelse
           endcase
    'Backlash':CorEvent,ev
    'Hbacklash':CorEvent,ev
    'Vbacklash':CorEvent,ev
    'TomoIter':CorEvent,ev
    'Convergence':CorEvent,ev
    'OSEMsubsets':CorEvent,ev
    'estcutoff':CorEvent,ev
    'tomoavgcutoff':CorEvent,ev
    'tomoedgecutoff':CorEvent,ev
    'tomocircleR':CorEvent,ev
    'tomoXmult':CorEvent,ev
    'tomoYmult':CorEvent,ev
    'tomoXshift':CorEvent,ev
    'tomoYshift':CorEvent,ev
    'tomodrho':CorEvent,ev
    'filtersize':CorEvent,ev
    'addangle': CorEvent,ev
    'Projection center':CorEvent,ev
    'Zwidth':CorEvent,ev
    'Wa':    CorEvent,ev
    'Wb':    CorEvent,ev
    'Wangle':CorEvent,ev
    'Zthreshold':CorEvent,ev
    'CutPos':CutInsertEvent,ev
    'Cutn':    CutInsertEvent,ev
    'ShiftRownr':ShiftRowEvent,ev
    'ShiftRown':ShiftRowEvent,ev
    'Cutoff':
    'AddGroups':AddXDIGroups,ev
    'DelGroups':DeleteXDIGroups,ev
    'BlendIndex':    begin
                    widget_control,ev.top,get_uvalue=list
                    widget_control,ev.id,get_value=str
                    (**list.pgradient).index=long(str)
                    RefreshGradientImage,ev
                    RefreshWindowXDI,ev
                    endcase
    'ngradients':    begin
                    widget_control,ev.top,get_uvalue=list
                    widget_control,ev.id,get_value=str
                    n=long(str)
                    ResizeGradients,list.pgradients,n,list.pgradient
                    widget_control,ev.id,set_value=string(n)
                    RefreshGradientImage,ev
                    RefreshWindowXDI,ev
                    endcase
    else:
    endcase
    endcase
8:    begin; ----Event from widget_droplist----
    widget_control,ev.id,get_uvalue=uval
    case uval of
    'filterlist':begin
            widget_control,ev.top,get_uvalue=list
            list.filtertype=ev.index
            widget_control,ev.top,set_uvalue=list
               RefreshDisplayXDI,list
            endcase
    'ngradients':begin
            widget_control,ev.top,get_uvalue=list
            *list.pgradient=(*list.pgradients)[ev.index]
            RefreshGradientImage,ev
            RefreshWindowXDI,ev
            endcase
    'modes':begin
            widget_control,ev.top,get_uvalue=list
            (**list.pgradient).mode=ev.index
            RefreshGradientImage,ev
            RefreshWindowXDI,ev
            endcase
    'CombineType':begin
            widget_control,ev.top,get_uvalue=list,/no_copy
            list.CombineType = ev.index
            widget_control,ev.top,set_uvalue=list,/no_copy
            endcase
    endcase
    endcase
9:    begin; ----Event from widget_table----
    widget_control,ev.id,get_uvalue=uval
    case uval of
    'MaskTable':begin
            widget_control,ev.id,get_value=value
            widget_control,ev.top,get_uvalue=list
            i=list.NSelect+1
            (*list.pmstruc).(i).data=value
            endcase
    endcase
    endcase
else:
endcase

end; pro XRD_XDI_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro XRD_XDI,ev=ev

print,'=========================================START XDI Editor'
ControlDevice

if keyword_set(ev) then begin
    widget_control,ev.top,get_uvalue=list
    xdipath=list.xdipath
    xdifile=list.xdifile
    main=ev.top
    ev=ev
    sortseparator=list.sortseparator
    sortmethod=list.sortmethod
    readinfo=list.readinfo
endif else begin
    result=ReadINI('$XRD_XDI path:')
    xdipath=result.(0)
    xdifile=result.(1)
    main=0L
    ev=-1L
    sortseparator='_'
    sortmethod=1B
    
    readinfo=ReadINI('$ReadCCDInfo:')
endelse

OutID=SubscribeToOutID()

pgradients=DefaultGradients()
pgradient=ptr_new((*pgradients)[0])

list={listtype:20,$
    data:PTR_NEW(),$         ; 2D or 3D databox, containing XRD mappings
    dataorg:PTR_NEW(),$        ; original data (NULL when 2D databox)
    ptrtomocurrent:PTR_NEW(/alloc),$ ; pointer to the last reconstructed tomogram
    tomoweights:PTR_NEW(/alloc),$ ; pointer to tomogram weights

; dataorg corrections: DC normalization, reorientation
; additional corrections for data: backlash cor., zinger removal, normalization
; additional corrections for displayed data: filter

    normflag:0b,$            ; 1b => dataorg normalized
    filenames:PTR_NEW(),$    ; 1D or 2D databox, containing CCD image paths

    pmstruc:PTR_NEW(),$        ; physical mask data -> {...,ti:{n: number of areas in this group,data:...,labeli: row index in pmfields }}
    pmtags:PTR_NEW(),$        ; physical mask data tags (1 for each group) -> [...,'ti',...]
    pmfields:PTR_NEW(),$    ; physical mask datafields    -> [[strarr(npmfields)],[strarr(npmfields)]]
    npmfields:0,$            ; number of physical mask datafields

    path:xdipath,$          ; Directory of xdi file
    file:xdifile,$          ; Name of xdi file

    ScanDim:0,$                ; Scan or Map
    NMaps:0L,$                 ; Number of mappings
    MapCol:0L,$             ; Number of mappixels horizontal
    MapRow:0L,$             ; Number of mappixels vertical
    TRAN:0b,$                ; Dataset is transposed (default: first horiz. then vert.)
    REV1:0b,$                ; Dataset columns are reversed (default: LR)
    REV2:0b,$                ; Dataset rows are reversed (default: UB)
    REV3:0b,$                ; Dataset odd rows are reversed (default: JUMP)
    SWAP:0b,$                ; Dimensions swapped w.r.t. original?
    Hsign:1,$                ; LR positive, RL negative
    Vsign:1,$                ; UB positive, BU negative
    XSize:0L,$               ; Draw window size horizontal
    YSize:0L,$               ; Draw window size vertical
    mXSizeTomo:1.,$            ; X-size of tomogram: XSize*mXSizeTomo
    mYSizeTomo:1.,$            ; Y-size of tomogram: XSize*mYSizeTomo
    mXShiftTomo:0.,$        ; Shift the center of rotation
    mYShiftTomo:0.,$        ; Shift the center of rotation
    tomodrho:1.,$           ; sinogram translation step size divided by the tomogram pixel size
    
    NSelect:0,$              ; Selected mapping to display

    ReconType:0b,$            ; FBP:0, MLEM:1, ...
    arthigh:0b,$            ; Use area weights instead of distance weights
    altavg:1b,$                ; L1-norm for SIRT/SART
    recalcweights:1b,$        ; Recalculate tomography weights
    bshowconvergence:0,$    ; Show tomographic reconstruction convergence
    ARTweights:ptr_new(/alloc),$; Tomography weights
    Corb:bytarr(6),$        ; Do [backlash correction?, projection center correction?,Normalize so that each view in the
                            ;     sinogram has the same number of counts as the average total counts for a view, zinger removal, wobble correction, new wobble correction]
    Cor:fltarr(7),$            ; [Shift the even rows with this amount to correct for Backlash , Projection center, Zinger filter width,
                            ; Zinger threshold, Wobble-ellipse major axis, Wobble-ellipse minor axis, Wobble-ellipse angle]
    ManualBacklash:fltarr(2),$; Horizontal/vertical backlash correction (prior to sinogram backlash)
    WobblePtr:ptr_new(/alloc),$    ; Pointer to wobble row-shifts
    dynrecon:0b,$            ; Reconstruct image?
    filtersize:500,$        ; Sinogram FBP filter half width
    filtertype:1,$            ; Sinogram FBP filter type
    TomoIter:1,$            ; Tomography iterations
    TomoIterMax:100,$        ; Tomography iterations (maximal)
    bTomoIter:1b,$            ; Use fixed number of iterations
    TomoConvergence:10.,$    ; Convergence criterion for tomographic reconstruction
    OSEMsubsets:1,$            ; OSEM subsets
    estcutoff:0.,$            ; MLEM estimated sinogram cutoff
    addangle:0.,$            ; First angle in sinogram
    tomoavg:0b,$            ; Ray sums are no sums but averages
    tomoavgcutoff:80.,$        ; Normalization tomogram cutoff
    tomoedge:0b,$            ; Remove edge after reconstruction
    tomoedgecutoff:1ul,$    ; Edge width to remove
    tomocircle:0b,$            ; Circle clipping
    tomocircleR:0.,$        ; Radius of circle clip
    tomopositive:1b,$        ; Positive tomogram
    tomopostcorrect:1b,$    ; Tomogram constraints after reconstruction
    tomoequalweights:0b,$    ; Equal tomogram weights
    
    updatesinogram:1b,$        ; Update sinogram
    detectorpos:0,$            ; Position of the detector in the tomogram

    CombineNorm:0b,$
    CombineType:0,$     

    expo:0.,$                ; Exposure time
    coord:fltarr(2,2),$        ; scan coordinates for default orientation: moveh,LR,UB,LR(+),UB(+)
    unit:strarr(2),$        ; scan coordinate units

    nrROI:0,$               ; Number of closed ROI's
    nrPoints:PTR_NEW(),$     ; Number of Points for each closed ROI
    pverts:-1L,$               ; Points to last added point
    xverts:PTR_NEW(),$       ; X-coord of vertices
    yverts:PTR_NEW(),$       ; Y-coord of vertices
    highlight:PTR_NEW(),$    ; 1 => highlight region
    names:PTR_NEW(),$       ; Names of regions
    
    pgradient:pgradient,$    ; Current color blend for 1 map
    pgradients:pgradients,$ ; Color blends for multiple maps
    gradientfill:0l,$        ; Fill background pixels (iterations)
    pixindexgradient:0l,$   ; 
    drawindexgradient:0l,$  ; 
    gradientclipev:0,$      ; 
    gradientrefreshname:'RefreshDisplayXDI',$ ; 
  
    scaleimg:100,$             ; grayscaling: 100 -> [min,max], 50 -> [min,min+((max-min)*50%)], 150 -> [min+((max-min)*50%),min+((max-min)*100%)],
                            ;  general -> min+(max-min)*[(scaleimg-100)>0,scaleimg]<100
    ctl:{type:0,nr:-3,R:ptr_new(),G:ptr_new(),B:ptr_new()},$ ; Color table
    sfh:0.2,$                ; sfh=MapCol/XSize
    sfv:0.2,$                ; sfv=MapRow/YSize
    MState:'',$             ; Mouse state
    OutID:OutID,$            ; Output log
    charpr:0.01,$            ; charpr*list.XSize => charsize
    labeli:0,$                ; output label index
    nhorizontal:4,$            ; Save all nx
    markers:'',$            ; Positions of 1D plot markers
    labelcol:255,$            ; Color table index for label colors
    boardercol:0b,$            ; Color table index for boarder colors
    bopenformat:0b,$        ; Open formatted
    openformat:'unformatted',$    ; Open formatted

    sortseparator:sortseparator,$ ; for sorting files
    sortmethod:sortmethod,$ ; 1 => sort strings, 0 sort strings with offset

    cutpos:ptr_new(),$        ; Position before which to cut
    cutn:ptr_new(),$        ; cut n pixels
    rowshift:ptr_new(),$    ; Rows/columns to shift
    pautofill:ptr_new(),$    ; Autofill redirections
    preplacep:ptr_new(),$    ; Pixel replace

    xs:0L,$                  ; Static corner of LOF selection
    ys:0L,$                  ; Static corner of LOF selection
    legendpos:fltarr(2),$    ; Position of legend in normal coord.
    sysvar:PointSysVar(),$    ; System variable for coord. convertion

    readinfo:readinfo,$     ; Normalize EDF on reading

    main:main,$             ; ID of Main Program
    ev:ev,$                  ; Event structure that triggered the start of this application
    top:0L,$                 ; ID of Top Level Base
    drawindex:0L,$            ; index of draw widget
    draw:0L,$                 ; ID of overview drawwidget
    TomoReconID:0L,$        ; ID of reconstruct draw widget
    TomoReconIndex:0L,$        ; ID of reconstruct draw widget
    pixindex:0L}             ; index of pixwindow

base = WIDGET_BASE(MBAR=bar,/row,title='XDI Editor:')

    base2=widget_base(base,/column)
        baset=widget_base(base2,/row,uname='drawadd')
        
       draw = WIDGET_DRAW(baset,uname='draw',uvalue=0)
         Window, /Free, /Pixmap
         pixindex = !D.Window
         list.draw=draw
         list.pixindex=pixindex
       text = WIDGET_TEXT(base2, YSIZE=3, XSIZE=60,uname='position')
       base1c=widget_base(base2,/column,/frame,sensitive=0,uname='Browse')
               label = widget_label(base1c,value='Display Group:')
               slider=widget_slider(base1c,minimum=0,maximum=1,uname='NRslider',uvalue='NRslider',value=0)
               label = widget_label(base1c,value='Group information:')
               tmp = WIDGET_base(base1c,uname='Mask')

    tab=widget_tab(base,MULTILINE =5)
        base1a=widget_base(tab,/column,title='Image Scaling')
               label = widget_label(base1a,value='Intensity Scaling:')
            slider=widget_slider(base1a,minimum=0,maximum=200,uname='scale',uvalue='scale',value=list.scaleimg)

            base1a=widget_base(base1a,/column,sensitive=0,uname='magnification')
               label=widget_label(base1a,value='Magnification:')
               text = widget_text(base1a,value=stringr(1/list.sfh),uvalue='Magx',uname='Magx',/editable)
               text = widget_text(base1a,value=stringr(1/list.sfv),uvalue='Magy',uname='Magy',/editable)

        base1b=widget_base(tab,/column,sensitive=0,uname='3DROI',title='Region Properties')
               button=widget_button(base1b,value='Close Region')
               button=widget_button(base1b,value='Fit Gaussians')
               button=widget_button(base1b,value='Fit Steps')
               button=widget_button(base1b,value='Wall thickness')
               button=widget_button(base1b,value='Delete Select')
               text = WIDGET_TEXT(base1b, uname='Names',uvalue='Names',/editable)
               label = widget_label(base1b,value='Select Closed Region:')
               slider=widget_slider(base1b,minimum=0,maximum=1,uname='Highslider',uvalue='Highslider',value=0)
               baset=widget_base(base1b,/column,/nonexclusive)
                 button = WIDGET_BUTTON(baset, VALUE='Fill',uvalue=0,uname='Fill')
                 button = WIDGET_BUTTON(baset, VALUE='Continuous',uvalue=0,uname='Continuous')

        base3a=widget_base(tab,/column,title='Options')
               base4=widget_base(base3a,/column,/nonexclusive)
                 buttonSGsep = WIDGET_BUTTON(base4, VALUE='Save Groups Separate',uvalue=1,uname='Save Groups Separate')
                 buttonSTsep = WIDGET_BUTTON(base4, VALUE='Save Tomo Separate',uvalue=1,uname='Save Tomo Separate')
                 button = WIDGET_BUTTON(base4, VALUE='Save without boarders',uvalue=0,uname='Save without boarders')
                 button = WIDGET_BUTTON(base4, VALUE='Save 1D data in 1 plot',uvalue=0,uname='Save 1D data in 1 plot')
                 button = WIDGET_BUTTON(base4, VALUE='Save with labels',uvalue=0,uname='Save with labels')
                 button = WIDGET_BUTTON(base4, VALUE='Update Main',uvalue=0,uname='Update Main')
                 button = WIDGET_BUTTON(base4, VALUE='Sticky xboarders for insert XDI',uvalue=0,uname='Sticky xboarders for insert XDI')
                 button = WIDGET_BUTTON(base4, VALUE='Sticky yboarders for insert XDI',uvalue=0,uname='Sticky yboarders for insert XDI')
                 ;baset=widget_base(base1a,/column,/nonexclusive)
                 button = WIDGET_BUTTON(base4, VALUE='Log Scale',uvalue=0,uname='Log Scale')
                 buttonSP = WIDGET_BUTTON(base4, VALUE='Binned pixels',uvalue=1,uname='Binned pixels')
                 button= WIDGET_BUTTON(base4, VALUE='Add contour',uvalue=0,uname='Add contour')    

        base3b=widget_base(tab,/column,uname='basemode',sensitive=1,title='Mode')
               baset=widget_base(base3b,/column)
                 button = WIDGET_BUTTON(baset, VALUE='Draw Region',uname='Draw Region')
                 button = WIDGET_BUTTON(baset, VALUE='Draw 1D Region',uname='Draw 1D Region')
                 button = WIDGET_BUTTON(baset, VALUE='Make LOF',uname='Make LOF')
                 button = WIDGET_BUTTON(baset, VALUE='Make LOF from Int',uname='Make LOF from Int')
                 button = WIDGET_BUTTON(baset, VALUE='Pattern',uname='Pattern')
                 button = widget_button(baset,value='Draw Sinus',uname='Draw Sinus')
                 button = WIDGET_BUTTON(baset, VALUE='Check ddist',uname='Check ddist')
                 button = widget_button(baset,value='Extract 1D',uname='Extract 1D')
                 button = WIDGET_BUTTON(baset, VALUE='Resize',uname='Resize')
                 button = WIDGET_BUTTON(baset, VALUE='Reform',uname='Reform')
                 button = WIDGET_BUTTON(baset, VALUE='Move Legend',uname='Move Legend')

           base3c=widget_base(tab,/column,sensitive=0,uname='BaseShow',title='Show')
               button=widget_button(base3c,value='Show 3D',uname='Show 3D')

        base3d=widget_base(tab,/column,uname='CutBase',title='Pixel juggling',sensitive=0)
            base3dt=widget_base(base3d,/frame,/column)
            tmp = WIDGET_LABEL(base3dt, VALUE='Cut/Insert pixels: ')
            baset=widget_base(base3dt,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Positions: ')
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='CutPos',uname='CutPos',/editable)
            baset=widget_base(base3dt,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Number of pixels (- to cut): ')
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='Cutn',uname='Cutn',/editable)
            base3dt=widget_base(base3d,/frame,/column)
            tmp = WIDGET_LABEL(base3dt, VALUE='Shift rows/columns: ')
            baset=widget_base(base3dt,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Row number (- for column): ')
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='ShiftRownr',uname='ShiftRownr',/editable)
            baset=widget_base(base3dt,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Number of pixels to shift: ')
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='ShiftRown',uname='ShiftRown',/editable)
            base3dt=widget_base(base3d,/frame,/column)
            tmp = WIDGET_LABEL(base3dt, VALUE='Fill pixels: ')
            baset=widget_base(base3dt,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Pixel value (==,!=,<,>,<=,>=): ')
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='Cutoff',uname='Cutoff',/editable)
            baset=widget_base(base3dt,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Example (value or x*group): ')
                tmp = WIDGET_TEXT(baset, VALUE='>=10.23 or <=2*g3')
            baset=widget_base(base3dt,/row)
                tmp = WIDGET_BUTTON(baset, VALUE='DOWN',uvalue='CDOWN',uname='CDOWN')
                tmp = WIDGET_BUTTON(baset, VALUE='UP',uvalue='CUP',uname='CUP')
                tmp = WIDGET_BUTTON(baset, VALUE='LEFT',uvalue='CLEFT',uname='CLEFT')
                tmp = WIDGET_BUTTON(baset, VALUE='RIGHT',uvalue='CRIGHT',uname='CRIGHT')
                tmp = WIDGET_BUTTON(baset, VALUE='REPLACE',uvalue='CREPLACE',uname='CREPLACE')
                tmp = WIDGET_BUTTON(baset, VALUE='LOF')
                tmp = WIDGET_BUTTON(baset, VALUE='CLEAR',uvalue='CCLEAR',uname='CCLEAR')
            base3dt=widget_base(base3d,/frame,/column)
            tmp = WIDGET_LABEL(base3dt, VALUE='Replace pixels: ')
            baset=widget_base(base3dt,/row)
                tmp = WIDGET_BUTTON(baset, VALUE='Select Pixels',uvalue='Select Pixels',uname='Select Pixels')
                tmp = WIDGET_BUTTON(baset, VALUE='CLEAR ',uvalue='CCLEAR2',uname='CCLEAR2')


        base3e=widget_base(tab,/column,uname='AddBase',title='Combine Groups')
            baset=widget_base(base3e,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Groups to combine: ')
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='AddGroups',uname='AddGroups',/editable)
            baset=widget_base(base3e,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Example: ')
                tmp = WIDGET_TEXT(baset, VALUE='0,3,6')
            tmp = WIDGET_BUTTON(base3e, VALUE='Combine with other XDI...',uvalue='AddGroupsXDI')
            tmp = WIDGET_LABEL(base3e, VALUE='The default operation is summation. Following options are available:')
            baset=widget_base(base3e,/column,/nonexclusive)
                tmp=widget_button(baset,VALUE='normalized',uname='CombineNorm')
            tmp = widget_droplist(base3e,value=['Add','Maximum','Divide (by last)','Subtract (by last)','Multiply'],uname='CombineType',uvalue='CombineType')

        base3f=widget_base(tab,/column,uname='OverlapBase',title='Overlap Groups')
            baset=widget_base(base3f,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Groups to overlap: ')
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='OverlapGroups',uname='OverlapGroups',/editable)
            baset=widget_base(base3f,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Example: ')
                tmp = WIDGET_TEXT(baset, VALUE='0 or (5 and 6) and not 7')
            baset=widget_base(base3f,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Thresholds: ')
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='Thresholds',uname='Thresholds',/editable)
            baset=widget_base(base3f,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Example: ')
                tmp = WIDGET_TEXT(baset, VALUE='<1000,>5.3,0,20')
            baset=widget_base(base3f,/row)
                tmp = WIDGET_LABEL(baset, VALUE='RGB Colors: ')
                tmp = WIDGET_TEXT(baset, VALUE='[255,0,0][0,255,0][0,0,255]',uvalue='RGB',uname='RGB',/editable)
            button=widget_button(base3f, value='Create Overlap', uvalue='Create Overlap')

        base3h=widget_base(tab,/column,uname='BlendBase',title='Group blending',sensitive=0)
            baset=widget_base(base3h,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Number of layers:',xs=150)
                ngradients=n_elements(*list.pgradients)
                tmp=widget_text(baset,value=stringr(ngradients),/editable,uvalue='ngradients',uname='ngradients')
            baset=widget_base(base3h,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Layer to adapt:',xs=150)
                tmp=widget_droplist(baset,value=string(lindgen(ngradients)+1,format='(I30)'),uvalue='ngradients',uname='ngradientslist')
                baset=widget_base(baset,/row,/nonexclusive)
                    tmp=widget_button(baset,value='Enable',uname='Enable')
            baset=widget_base(base3h,/row)
            tmp = WIDGET_LABEL(baset, VALUE='Layer mode:',xs=150)
                tmp=widget_droplist(baset,value=string(blendmodes(),format='(A25)'),uvalue='modes',uname='modes')
            baset=widget_base(base3h,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Group index:',xs=150)
                tmp = WIDGET_TEXT(baset, VALUE=stringr((**list.pgradient).index),uname='BlendIndex',uvalue='BlendIndex',/editable)
            baset=widget_base(base3h,/row)
                tmp = widget_button(baset,value='Bring Forward')
                tmp = widget_button(baset,value='Send Backward')
            baset=widget_base(base3h,/row)
                tmp = widget_button(baset,value='Reset Gamma')
                tmp = widget_button(baset,value='Reset Clipping')
                tmp = widget_button(baset,value='Reset Opacity')

            GradientWidget,base3h,450,80,list
            
            baset=widget_base(base3h,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Fill iterations:',xs=150)
                tmp = widget_slider(baset,minimum=0,maximum=50,uname='GradientFill',uvalue='GradientFill',value=list.gradientfill)
            baset=widget_base(base3h,/row)
                tmp = widget_button(baset,value='Create Blend Group',uvalue='Create Blend Group')
                tmp = widget_button(baset,value='Save Blend')
                tmp = widget_button(baset,value='Restore Blend')
                tmp = widget_button(baset,value='Save Gradient image')

        base3g=widget_base(tab,/column,uname='AddBase',title='Delete Groups')
            baset=widget_base(base3g,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Groups to delete: ')
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='DelGroups',uname='DelGroups',/editable)
            baset=widget_base(base3g,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Example: ')
                tmp = WIDGET_TEXT(baset, VALUE='0,3,6')

        base3g=widget_base(tab,/column,uname='CorBase',title='Correlate Groups')

            baset=widget_base(base3g,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Subimage cross-correlation tollerance: ')
                tmp = WIDGET_TEXT(baset, VALUE='0.9',uvalue='CorGroupRtol',uname='CorGroupRtol',/editable)
            
            baset=widget_base(base3g,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Subimage size: ')
                tmp = WIDGET_TEXT(baset, VALUE='0',uvalue='CorGroupnpix',uname='CorGroupnpix',/editable)
                
            baset=widget_base(base3g,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Subimage shift tollerance: ')
                tmp = WIDGET_TEXT(baset, VALUE='1',uvalue='CorGroupshifttol',uname='CorGroupshifttol',/editable)

            baset=widget_base(base3g,/row)
                tmp = WIDGET_LABEL(baset, VALUE='Groups to correlate: ')
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='CorGroup1',uname='CorGroup1',/editable)
                tmp = WIDGET_TEXT(baset, VALUE='',uvalue='CorGroup2',uname='CorGroup2',/editable)

            button=widget_button(base3g, value='Correlate', uvalue='Correlate')

        InitBase=widget_base(tab,/column,sensitive=0,uname='InitBase',title='Orientation')
             baset=widget_base(InitBase,/row,/exclusive)
            Movehbutton=widget_button(baset,value='First horizontal',uname='First horizontal')
            Movevbutton=widget_button(baset,value='First vertical',uname='First vertical')

            baset=widget_base(InitBase,/row,/exclusive,uname='REV1')
            LRbutton=widget_button(baset,value='Left-Right',uname='Left-Right')
            RLbutton=widget_button(baset,value='Right-Left',uname='Right-Left')

            baset=widget_base(InitBase,/row,/exclusive,uname='REV2')
            UBbutton=widget_button(baset,value='Top-Bottom',uname='Top-Bottom')
            BUbutton=widget_button(baset,value='Bottom-Top',uname='Bottom-Top')

            baset=widget_base(InitBase,/row,/exclusive,uname='REV3')
            Jbutton=widget_button(baset,value='Jump',uname='Jump')
            Cbutton=widget_button(baset,value='Continue',uname='Continue')
            
            basee=widget_base(InitBase,/row)
            baset=widget_base(basee,/row,uname='HORIZ')
            label=widget_label(baset,value='Hori+:')
            button=widget_button(baset,value='L-->R',uname='L--R')
            baset=widget_base(basee,/row,uname='VERT')
            label=widget_label(baset,value='Vert+:')
            button=widget_button(baset,value='T-->B',uname='T--B')

            baset=widget_base(InitBase,/row,uname='Swap Dim base')
            button = WIDGET_BUTTON(baset, VALUE='Swap Dim',uname='Swap Dim')
            button = WIDGET_BUTTON(baset, VALUE='Change coord.',uname='Change coord.')

        base4a=widget_base(tab,/column,uname='Statistics',title='Statistics')
            baset=widget_base(base4a,/row)
            label = widget_label(baset,value='Mean +/- stddev:',xsize=100)
            text=widget_text(baset,uname='Stat_Mean')
            baset=widget_base(base4a,/row)
            label = widget_label(baset,value='Min:',xsize=100)
            text=widget_text(baset,uname='Stat_Min')
            baset=widget_base(base4a,/row)
            label = widget_label(baset,value='Max:',xsize=100)
            text=widget_text(baset,uname='Stat_Max')
            baset=widget_base(base4a,/row)
            label = widget_button(baset,value='FWHM:',xsize=100)
            text=widget_text(baset,uname='Stat_FWHM')
            
        base4a=widget_base(tab,/column,uname='Backlash base',sensitive=0,title='Backlash')
            baset=widget_base(base4a,/row)
            label = widget_label(baset,value='Horizontal backlash:',xsize=100)
            text=widget_text(baset,uname='Hbacklash',uvalue='Hbacklash',value='0.0000',/editable)
            baset=widget_base(base4a,/row)
            label = widget_label(baset,value='Vertical backlash:',xsize=100)
            text=widget_text(baset,uname='Vbacklash',uvalue='Vbacklash',value='0.0000',/editable)
            
        base4a=widget_base(tab,/column,uname='Tomo base',sensitive=0,title='Tomography')
               baset=widget_base(base4a,/column,/nonexclusive)
               button = WIDGET_BUTTON(baset, VALUE='Tomo Reconstruct',uvalue=0,uname='Tomo Reconstruct')
               button = WIDGET_BUTTON(baset, VALUE='Suspend update',uvalue=0,uname='Suspend update')
               button = WIDGET_BUTTON(baset, VALUE='Update sinogram',uvalue=0,uname='Update sinogram')
               if list.updatesinogram then widget_control,button,set_button=1
               ;button = WIDGET_BUTTON(baset, VALUE='Dynamic reconstruct (slow)',uname='Dynamic reconstruct')


            subtab=widget_tab(base4a,MULTILINE =4,uname='SinCorBase')
            
            baset = widget_base(subtab,/column,title='Sinogram corrections')
                baset3 = widget_base(baset,/row)

                baset1 = widget_base(baset3,/column,/nonexclusive)
                    button = WIDGET_BUTTON(baset1, VALUE='Backlash cor.',uvalue=0,uname='Backlash cor.')
                       button = WIDGET_BUTTON(baset1, VALUE='Proj. center cor',uvalue=0,uname='Proj. center cor')
                       button = WIDGET_BUTTON(baset1, VALUE='Zinger removal',uname='Zinger removal')
                       button = WIDGET_BUTTON(baset1, VALUE='Normalize',uname='Normalize')
                       button = WIDGET_BUTTON(baset1, VALUE='Wobble',uname='Wobble')
                       button = WIDGET_BUTTON(baset1, VALUE='WobbleMan',uname='WobbleMan')
                baset2 = widget_base(baset3,/column)
                    text = WIDGET_TEXT(baset2,uvalue='Backlash',uname='Backlash',/editable)
                    text = WIDGET_TEXT(baset2,uvalue='Projection center',uname='Projection center',/editable)
                    text = WIDGET_TEXT(baset2,uvalue='Zwidth',uname='Zwidth',/editable)
                    text = WIDGET_TEXT(baset2,uvalue='Zthreshold',uname='Zthreshold',/editable)
                    ;text = WIDGET_TEXT(baset2,uvalue='Wa',uname='Wa',/editable)
                    ;text = WIDGET_TEXT(baset2,uvalue='Wb',uname='Wb',/editable)
                    ;text = WIDGET_TEXT(baset2,uvalue='Wangle',uname='Wangle',/editable)
                baset4 = widget_base(baset2,/row)
                    button = WIDGET_BUTTON(baset4, VALUE='Fit',uname='Fit')
                    button = WIDGET_BUTTON(baset4, VALUE='COPManual',uname='COPManual')
                    
            base4b = widget_base(subtab,/column,title='Tomogram corrections')
                baset = widget_base(base4b,/column,/nonexclusive)
                    button=widget_button(baset,value='Remove edge',uvalue='tomoedge',uname='tomoedge')
                    button=widget_button(baset,value='Clip circle',uvalue='tomocircle',uname='tomocircle')
                    button=widget_button(baset,value='Positive',uvalue='tomopositive',uname='tomopositive')
                    button=widget_button(baset,value='Constraints after reconstruction',uvalue='tomopostcorrect',uname='tomopostcorrect')

                baset1 = widget_base(base4b,/column)
                    baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='Edge(pixels):')
                        text = WIDGET_TEXT(baset,uname='tomoedgecutoff',uvalue='tomoedgecutoff',editable=list.tomoedge)
                    baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='Clip circle radius(pixels):')
                        text = WIDGET_TEXT(baset,uname='tomocircleR',uvalue='tomocircleR',editable=list.tomocircle)
                    baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='X enlarge:')
                        text = WIDGET_TEXT(baset,uname='tomoXmult',uvalue='tomoXmult',/editable)
                    baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='Y enlarge:')
                        text = WIDGET_TEXT(baset,uname='tomoYmult',uvalue='tomoYmult',/editable)
                    baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='X shift:')
                        text = WIDGET_TEXT(baset,uname='tomoXshift',uvalue='tomoXshift',/editable)
                    baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='Y shift:')
                        text = WIDGET_TEXT(baset,uname='tomoYshift',uvalue='tomoYshift',/editable)
                     baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='Pixel size (units):')
                        text = WIDGET_TEXT(baset,uname='tomodrho',uvalue='tomodrho',/editable)
                        
            base4b = widget_base(subtab,/row,title='Reconstruction')
                baset = widget_base(base4b,/column,/exclusive)
                    button=widget_button(baset,value='FBP',uvalue='FBP',uname='FBP')
                    button=widget_button(baset,value='MLEM(radon)',uvalue='MLEM(radon)',uname='MLEM(radon)')
                    button=widget_button(baset,value='ART',uvalue='ART',uname='ART')
                    button=widget_button(baset,value='SIRT',uvalue='SIRT',uname='SIRT')
                    button=widget_button(baset,value='SART',uvalue='SART',uname='SART')
                    button=widget_button(baset,value='SIRT(radon)',uvalue='SIRT(radon)',uname='SIRT(radon)')
                    button=widget_button(baset,value='MLEM',uvalue='MLEM',uname='MLEM')
                    button=widget_button(baset,value='OSEM',uvalue='OSEM',uname='OSEM')
            base4b = widget_base(base4b,/column)
    
                baset1 = widget_base(base4b,/column,uname='FBPi',sensitive=list.ReconType eq 0)
                       baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='Filter type:')
                        drop = WIDGET_DROPLIST(baset, VALUE=['None', 'RamLak', 'Shepp_Logan', $
                        'LP_Cosine', 'Gen_Hamming', 'Ramp'],uname='filterlist',uvalue='filterlist')
                     baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='Filter size:')
                        text = WIDGET_TEXT(baset,uname='filtersize',uvalue='filtersize',/editable)
                baset1 = widget_base(base4b,/column,uname='MLEMi',sensitive=list.ReconType eq 1 or list.ReconType eq 3 or list.ReconType eq 4)
                    baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='EstSin cutoff:')
                        text = WIDGET_TEXT(baset,uname='estcutoff',uvalue='estcutoff',/editable)
                baset1 = widget_base(base4b,/column,uname='ARTi',sensitive=list.ReconType ne 0)
                    baset = widget_base(baset1,/column)
                        basett = widget_base(baset,/column,/nonexclusive)
                            button=widget_button(basett,value='Fixed iterations',uvalue='Fixed iterations',uname='Fixed iterations')
                    baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='Iterations:')
                        text = WIDGET_TEXT(baset,uname='TomoIter',uvalue='TomoIter',/editable)
                    baset = widget_base(baset1,/row,uname='ConvergenceBase',sensitive=~list.bTomoIter)
                        label = WIDGET_LABEL(baset, VALUE='Convergence (%):')
                        text = WIDGET_TEXT(baset,uname='Convergence',uvalue='Convergence',/editable)
                        
                    baset = widget_base(baset1,/row,uname='OSEMsubsetsbase',sensitive=list.ReconType eq 4)
                        label = WIDGET_LABEL(baset, VALUE='Subsets:')
                        text = WIDGET_TEXT(baset,uname='OSEMsubsets',uvalue='OSEMsubsets',/editable)
                baset1 = widget_base(base4b,/column)
                    baset = widget_base(baset1,/column,/nonexclusive)
                        button=widget_button(baset,value='L1-norm',uvalue='L1-norm',uname='L1-norm',sensitive=list.ReconType ge 5)
                        button=widget_button(baset,value='Area intersect',uvalue='Heavy weights',uname='Heavy weights',sensitive=list.ReconType ge 3)
                        button=widget_button(baset,value='Recalculate weights',uvalue='Recalculate weights',uname='Recalculate weights',sensitive=list.ReconType ge 3)
                        button=widget_button(baset,value='Show Convergence',uvalue='Show Convergence',uname='Show Convergence',sensitive=list.ReconType ne 0)
                        

            base4b = widget_base(subtab,/column,title='Weights')
                baset = widget_base(base4b,/column,/nonexclusive)
                    button=widget_button(baset,value='Ray average (no sum)',uvalue='tomoavg',uname='tomoavg')
                    button=widget_button(baset,value='Binarize weights',uvalue='Equal weights',uname='Equal weights')
                baset1 = widget_base(base4b,/column)
                    baset = widget_base(baset1,/row)
                        label = WIDGET_LABEL(baset, VALUE='Weights cutoff(%):')
                        text = WIDGET_TEXT(baset,uname='tomoavgcutoff',uvalue='tomoavgcutoff',editable=list.tomoavg)
                    button = widget_button(baset1,value='Current tomogram as weights')
                    
            baset = widget_base(base4a,/row)
                label = WIDGET_LABEL(baset, VALUE='Angle offset(deg):')
                text = WIDGET_TEXT(baset,uname='addangle',uvalue='addangle',/editable)
                button = widget_button(baset,value='Animate Tomo')


    menu1 = WIDGET_BUTTON(bar, VALUE='File', /MENU)

       buttonopen = WIDGET_BUTTON(menu1, VALUE='Open XDI')
       menu2 = WIDGET_BUTTON(menu1, VALUE='Open...',/menu,UNAME='Open...',sensitive=0)
           button1 = WIDGET_BUTTON(menu2, VALUE='Add XDI')
        button1 = widget_button(menu2,value='Insert XDI',uname='Insert XDI')
       
      menu2 = WIDGET_BUTTON(menu1, VALUE='Save As...',/menu,UNAME='Save As...',sensitive=0)
        button = WIDGET_BUTTON(menu2, VALUE='Save XDI')
           button = WIDGET_BUTTON(menu2, VALUE='Save CHI')
           button = WIDGET_BUTTON(menu2, VALUE='Save CHI/TIFF')
       
       menu2 = WIDGET_BUTTON(menu1, VALUE='Save...',/menu,UNAME='Save...',sensitive=0)
           button = WIDGET_BUTTON(menu2, VALUE='Save Tomo')
           button4 = WIDGET_BUTTON(menu2, VALUE='Save RGN',uname='Save RGN')
           
       menu2 = WIDGET_BUTTON(menu1, VALUE='Save image...',/menu,UNAME='Save image...',sensitive=0)
           button9 = WIDGET_BUTTON(menu2, VALUE='Save Group',uname='Save Group')
           button9 = WIDGET_BUTTON(menu2, VALUE='Save All',uname='Save All')
           button9 = WIDGET_BUTTON(menu2, VALUE='Save colorbar',uname='Save colorbar')
           button9 = WIDGET_BUTTON(menu2, VALUE='Save RGB triangle',uname='Save RGB triangle')
           
       menu3 = WIDGET_BUTTON(menu1, VALUE='Load...',/menu,UNAME='Load...',sensitive=0)
           button3 = WIDGET_BUTTON(menu3, VALUE='Load RGN',uname='Load RGN')
           button3 = WIDGET_BUTTON(menu3, VALUE='Load Filenames',uname='Load Filenames')
       menu4 = WIDGET_BUTTON(menu1, VALUE='Convert...',/menu,UNAME='Convert...')
           button3 = WIDGET_BUTTON(menu4, VALUE='XDI to IMS')
           button3 = WIDGET_BUTTON(menu4, VALUE='XDI to EDF')
           button3 = WIDGET_BUTTON(menu4, VALUE='XDI to TIFF')
           button3 = WIDGET_BUTTON(menu4, VALUE='XDI to ATHENA')
           button3 = WIDGET_BUTTON(menu4, VALUE='IMS to XDI')
           button3 = WIDGET_BUTTON(menu4, VALUE='EDF to XDI')
           button3 = WIDGET_BUTTON(menu4, VALUE='TIFF to XDI')
           button3 = WIDGET_BUTTON(menu4, VALUE='PyMCA DAT to XDI')
           button3 = WIDGET_BUTTON(menu4, VALUE='MicroXAS FastXRF to XDI')
       button2 = WIDGET_BUTTON(menu1, VALUE='Exit',uname='Exit',/separator)
   menu1 = WIDGET_BUTTON(bar, VALUE='Normalization', /MENU,uname='Normalization',sensitive=0)
        button = WIDGET_BUTTON(menu1, VALUE='Flux or LiveTime from SPE')
        button = WIDGET_BUTTON(menu1, VALUE='Flux from FIO')
        button = WIDGET_BUTTON(menu1, VALUE='Flux from MCA')
        button = WIDGET_BUTTON(menu1, VALUE='Flux from IMS')
        button = WIDGET_BUTTON(menu1, VALUE='Flux from SPEC')
        button = WIDGET_BUTTON(menu1, VALUE='Flux from TIFF')
        button = WIDGET_BUTTON(menu1, VALUE='Flux from EDF')
        button = WIDGET_BUTTON(menu1, VALUE='Flux from XDI')
        button = WIDGET_BUTTON(menu1, VALUE='Flux and DT from TIFF')
        button = WIDGET_BUTTON(menu1, VALUE='LiveTime from CCD images')
   
   menu1 = WIDGET_BUTTON(bar, VALUE='Processing',/MENU,uname='Processing',sensitive=1)
    button = widget_button(menu1,value='Default image output...')
       button = widget_button(menu1,value='XIAParse')
       button = widget_button(menu1,value='Prepare for tomorecon: SkyScan')
   
   menu1 = WIDGET_BUTTON(bar, VALUE='Options',/MENU,uname='Options',sensitive=1)
     button = WIDGET_BUTTON(menu1, VALUE='Parameters')
     button = widget_button(menu1,value='SortOptions')
     
WIDGET_CONTROL, base, /REALIZE,set_uvalue=list
widget_control,buttonSP,set_button=1
widget_control,buttonSGsep,set_button=1
widget_control,buttonSTsep,set_button=1

Xmanager,'XRD_XDI',base,/NO_BLOCK,CLEANUP='XRD_XDI_CleanUp'
widget_control,buttonopen,send_event={ID:buttonopen,TOP:base,HANDLER:base,SELECT:long(1)}
end;pro XRD_XDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
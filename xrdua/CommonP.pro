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

function XRDUAresourceDIR, AppDirname, AppDesc, AppReadmeText_
; Directory that is specific to and writable by the user running XRDUA

AuthorDirname = 'xrdua'
AuthorDesc = 'XRDUA.'
AuthorReadmeText = ['This is the user resource directory for XRDUA.', $
        '', $
        'It is safe to remove this directory, as it', $
        'will be recreated on demand. Note that all', $
        'settings will revert to their default settings.']
AuthorReadmeVersion = 1

AppReadmeText = [$
        AppReadmeText_ , $
        '', $
        'It is safe to remove this directory, as it', $
        'will be recreated on demand. Note that all', $
        'settings will revert to their default settings.']
AppReadmeVersion = 1

return, APP_USER_DIR(AuthorDirname, AuthorDesc, AppDirname, AppDesc, $
                  AppReadmeText, AppReadmeVersion, $
                  AUTHOR_README_TEXT=AuthorReadmeText, AUTHOR_README_VERSION=AuthorReadmeVersion)+path_sep()

end;function XRDUAresourceDIR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XRDUAconfigDIR
return, XRDUAresourceDIR('config','XRDUA configurations.','This is the user configuration directory for XRDUA.')
end;function XRDUAconfigDIR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XRDUAdlmDIR
return, XRDUAresourceDIR('DLM','DLM builds.','This is the user DLM build output directory for XRDUA.')
end;function XRDUAdlmDIR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XRDUAworkingDIR
return, XRDUAresourceDIR('WD','Working directory.','This is the XRDUA working directory.')
end;function XRDUAworkingDIR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function XCOMInterpol,x,y,ind1,ind2,xout,type
; xout is located between x[ind1] and x[ind2]
; 
; type:
; 1. Coherent scattering => log-log cubic-spline
; 2. Incoherent scattering => log-log cubic-spline
; 3. Photoelectric Absorption
; 4. Pair Production in nuclear field
; 5. Pair Production in electron field

x1=alog(x[ind1])
x2=alog(x[ind2])
y1=alog(y[ind1])
y2=alog(y[ind2])
return,exp((y2-y1)/(x2-x1)*(alog(xout)-x1)+y1)

type=keyword_set(type)
if type le 0 or type ge 5 then type=0

; http://physics.nist.gov/PhysRefData/Xcom/Text/chap3.html
case type of
0:    begin ; log-log linear interpolation
    x1=alog(x[ind1])
    x2=alog(x[ind2])
    y1=alog(y[ind1])
    y2=alog(y[ind2])
    return,exp((y2-y1)/(x2-x1)*(alog(xout)-x1)+y1)
    endcase
1:    begin ; log-log cubic-spline
    return,exp(spline(alog(x),alog(y),alog(xout)))
    endcase
2:    begin ; log-log cubic-spline
    return,exp(spline(alog(x),alog(y),alog(xout)))
    endcase
3:    begin ; log-log cubic-spline
    ; problem with edges
    return,exp(spline(alog(x),alog(y),alog(xout)))
    endcase
4:    begin ; log-log cubic-spline
    return,exp(spline(alog(x),alog(y*(1-x/1022d6)^3),alog(xout)))
    endcase
5:    begin ; log-log cubic-spline
    return,exp(spline(alog(x),alog(y*(1-x/2044d6)^3),alog(xout)))
    endcase
endcase

end;function XCOMInterpol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckDLMloaded,dlmname
CATCH, Error_status
IF Error_status ne 0 THEN return,0
return,call_function(dlmname+'_check') eq 0
end;function CheckDLMloaded
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isLittleEndian
a=uint(1)
swap_endian_inplace,a,/SWAP_IF_LITTLE_ENDIAN
return, a eq 256
end;function isLittleEndian
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function checkrelease
stack=SCOPE_TRACEBACK()
if n_elements(stack) lt 2 then return,0b
return,stack[1] eq 'IDLRTMAIN'
end;function checkrelease
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro printLastError
Catch, /Cancel
help, /last_message, output=errors
print, errors,format='(A)'
heap_gc
end;pro printLastError
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro OpWithDuplicate,ind,dest1,addvalue1,op1,dest2,addvalue2,op2,dest3,addvalue3,op3

bsecond=n_params() gt 4
bthird=n_params() gt 7

; Find repeats
h1=histogram(ind,REVERSE_INDICES=r1)

; Group by repeat count
h2=histogram(h1,min=1,omax=omax,REVERSE_INDICES=r2)

; No repeats
if r2[0] ne r2[1] then begin
    ; Indices in h1 with repeat count=1
    inddest=r2[r2[0] : r2[1]-1]
    
    ; Indices in dest with repeat count=1
    ;     r1[ r1[inddest[k]] : r1[inddest[k]+1]-1 ]
    ;    = r1[ r1[inddest[k]] : r1[inddest[k]] ]
    ;    = r1[ r1[inddest[k]] ]
    indval=r1[r1[inddest]]
    inddest=ind[indval]
    
    ; Add values to destination
    case op1 of
    '+':    dest1[inddest]+=addvalue1[indval]
    'or':    dest1[inddest]or=addvalue1[indval]
    'and':    dest1[inddest]and=addvalue1[indval]
    'mean': dest1[inddest]+=addvalue1[indval]
    endcase
    if bsecond then $
    case op2 of
    '+':    dest2[inddest]+=addvalue2[indval]
    'or':    dest2[inddest]or=addvalue2[indval]
    'and':    dest2[inddest]and=addvalue2[indval]
    '+mean':dest2[inddest]+=addvalue2[indval]
    endcase
    if bthird then $
    case op3 of
    '+':    dest3[inddest]+=addvalue3[indval]
    'or':    dest3[inddest]or=addvalue3[indval]
    'and':    dest3[inddest]and=addvalue3[indval]
    '+mean':dest3[inddest]+=addvalue3[indval]
    endcase
endif

; Handle repeat count >1
for j=1,omax-1 do begin
    if r2[j] ne r2[j+1] then begin
        ; Indices in h1 with repeat count=nrep
        inddest=r2[r2[j] : r2[j+1]-1]
        
        ; Indices in dest with repeat count=nrep
        ;     r1[ r1[inddest[k]] : r1[inddest[k]+1]-1 ]
        ;    = r1[ r1[inddest[k]] : r1[inddest[k]]+nrep-1 ]
        ;    = r1[ r1[inddest[k]] + lindgen(nrep) ]
        nrep=j+1
        indval=r1[rebin(r1[inddest],h2[j],nrep,/SAMPLE)+ $
                    rebin(lindgen(1,nrep),h2[j],nrep,/SAMPLE)]
        inddest=ind[indval[*,0]]
        
        ; Add values to destination
        case op1 of
        '+':    dest1[inddest]+=total(addvalue1[indval],2)
        'or':    dest1[inddest]or=total(addvalue1[indval],2) ne 0
        'and':    dest1[inddest]and=total(addvalue1[indval],2) eq nrep
        'mean':    dest1[inddest]+=total(addvalue1[indval],2)/nrep
        endcase
        if bsecond then $
        case op2 of
        '+':    dest2[inddest]+=total(addvalue2[indval],2)
        'or':    dest2[inddest]or=total(addvalue2[indval],2) ne 0
        'and':    dest2[inddest]and=total(addvalue2[indval],2) eq nrep
        'mean':    dest2[inddest]+=total(addvalue2[indval],2)/nrep
        endcase
        if bthird then $
        case op3 of
        '+':    dest3[inddest]+=total(addvalue3[indval],2)
        'or':    dest3[inddest]or=total(addvalue3[indval],2) ne 0
        'and':    dest3[inddest]and=total(addvalue3[indval],2) eq nrep
        'mean':    dest3[inddest]+=total(addvalue3[indval],2)/nrep
        endcase
    endif
endfor

end;pro OpWithDuplicate
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckFile,path

;res=file_search(path,/TEST_READ)
;return,res[0] ne ''

;if openr_safe(lun, path) then return,0b
;free_lun,lun
;return,1B

return,FILE_TEST(path,/read)
end;function CheckFile
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RestoreStruc,list,listsource,exclude=exclude

CATCH, Error_status
IF Error_status ne 0 THEN return
            
names=tag_names(list)
namessource=tag_names(listsource)
if n_elements(exclude) eq 0 then exclude=''
exclude=strupcase(exclude)

for destfield=0l,n_tags(list)-1 do begin

    ; Skip field?
    ind=where(exclude eq names[destfield],ct)
    if ct ne 0 then continue
    
    ; Copy field?
    ind=where(namessource eq names[destfield],ct)
    if ct eq 1 then begin
        sourcefield=ind[0]
        
        ; Free original field
        s=size(listsource.(sourcefield),/type)
        if s eq 8 then begin
            ; Field is a structure
            list_=list.(destfield)
            listsource_=listsource.(sourcefield)
            RestoreStruc,list_,listsource_,exclude=exclude
            list.(destfield)=list_
            listsource.(sourcefield)=listsource_
        endif else begin
            ; Free pointers and object references
            if s eq 10 or s eq 11 then heap_free,list.(destfield)
            
            ; Copy field
            list.(destfield)=listsource.(sourcefield)
            
            ; "Neutralize" source field
            case s of
            10:    listsource.(sourcefield)[*]=ptr_new()
            11: listsource.(sourcefield)[*]=obj_new()
            else:
            endcase
        endelse
    endif
endfor

; Free unused fields
heap_free,listsource
end;pro RestoreStruc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function capitalization,strarr
b=byte(strarr)
ind=where(b ge 97 and b le 122,ct)
if ct eq 0 then return,strarr
b[ind]-=32
return,string(b)
end;function capitalization
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MDIALOG_PICKFILE,_REF_EXTRA=extra
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    res=dialog_message(!ERROR_STATE.MSG,/info)
    return,''
ENDIF

return,DIALOG_PICKFILE(_EXTRA=extra)
end;function MDIALOG_PICKFILE
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CAPDIALOG_PICKFILE,_EXTRA=extra
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    res=dialog_message(!ERROR_STATE.MSG,/info)
    return,''
ENDIF

if n_elements(extra) ne 0 and !version.os_family ne 'Windows' then begin
    ind=where(tag_names(extra) eq 'FILTER',ct)
    filter=capitalization(extra.filter)
    if ct eq 1 then begin
        ind=where(filter ne extra.filter,ct)
        if ct ne 0 then extra.filter[ind]+=';'+filter[ind]
    endif
endif

return,DIALOG_PICKFILE(_EXTRA=extra)
end;function CAPDIALOG_PICKFILE
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function strucbytesize,struc
help,struc,/struc,output=res
res=res[0]
p0=strpos(res,'data length=')+12
p1=strpos(res,',',p0)
return,ulong64(strmid(res,p0,p1-p0))
end;function strucbytesize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function arraybytesize,a,float=float,signed=signed

float=0b
signed=1b
case size(a,/type) of
1:     n=1
2:     n=2
3:     n=4
4:     begin
    n=4
    float=1b
    endcase
5:     begin
    n=8
    float=1b
    endcase
6:     begin
    n=8
    float=1b
    endcase
7:     n=strlen(a[0])
9:     begin
    n=16
    float=1b
    endcase
12: begin
    n=2
    signed=0b
    endcase
13: begin
    n=4
    signed=0b
    endcase
14: n=8
15: begin
    n=8
    signed=0b
    endcase
else: n=0
endcase

return,n_elements(a)*ulong64(n)
end;function arraybytesize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro typetobyte,arr
; BIG ENDIAN
case size(arr,/type) of
1:    ; already byte
2:  arr=byte(arr,0,2,n_elements(arr))
3:  arr=byte(arr,0,4,n_elements(arr))
4:  arr=byte(arr,0,4,n_elements(arr))
5:    arr=byte(arr,0,8,n_elements(arr))
12: arr=byte(arr,0,2,n_elements(arr))
13: arr=byte(arr,0,4,n_elements(arr))
14: arr=byte(arr,0,8,n_elements(arr))
15: arr=byte(arr,0,8,n_elements(arr))
else: stop,'Type conversion not supported.'
endcase
end;pro typetobyte
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro bytetotype,arr,type
; BIG ENDIAN
case type of
1:    ; already byte
2:  arr=fix(arr,0,n_elements(arr)/2)
3:  arr=long(arr,0,n_elements(arr)/4)
4:  arr=float(arr,0,n_elements(arr)/4)
5:    arr=double(arr,0,n_elements(arr)/8)
12: arr=uint(arr,0,n_elements(arr)/2)
13: arr=ulong(arr,0,n_elements(arr)/4)
14: arr=long64(arr,0,n_elements(arr)/8)
15: arr=ulong64(arr,0,n_elements(arr)/8)
else: stop,'Type conversion not supported.'
endcase
end;pro bytetotype
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function basetodecimal,fraction,base,status=status,msg=msg
; str contains a fraction (number with base and point) 
; fraction = mantissa.base^exp
; 
; example:
; binary number (base=2) with point=4
; 1101.101     = 1.2^3 + 1.2^2 + 0.2^1 + 1.2^0 + 1.2^(-1) + 0.2^(-2) + 1.2^(-3)
;            = 1.2^m + 1.2^(m-1) + ...                             + 1.2^(m-n+1)
;            = [1.2^(n-1) + 1.2^(n-2) + ...                             + 1.2^0]*2^(m-n+1)

status=0

; Get point position, n digits, m leading digit exponent
str=fraction
point=strpos(str,'.')
if point eq -1 then begin
    n=strlen(str)
    point=n
endif else begin
    str= strsplit(str,'.',count=ct,/extract)
    case ct of 
    1: str=str[0]
    2: str=str[0]+str[1]
    else: return,0
    endcase
    n=strlen(str)
endelse
m=point-1

; Base in the correct format
n_=alog(2)/alog(base)*[64,32,16,8]
if n gt n_[0] then begin
    base_=ulong64(base)^ul64indgen(n)
    status or=1
endif else $
if n gt n_[1] then base_=ulong64(base)^ul64indgen(n) else $
if n gt n_[2] then base_=ulong(base)^ulindgen(n) else $
if n gt n_[3] then base_=uint(base)^uindgen(n) else base_=byte(base)^bindgen(n)

; Digits
str=strlowcase(str)
digits=byte(str)-(byte('0'))[0]

off=(byte('a')-byte('0'))[0]
ind=where(digits ge off and digits lt off+base-10,ct)
if ct ne 0 then digits[ind]+=10-off
    
ind=where(digits ge base,ct)
if ct ne 0 then begin
    digits[ind]=0
    status or= 2
endif

; Mantissa:
mantissa=total(base_*reverse(digits),/pres)
exp=m-n+1

case status of
0: msg='no error'
1: msg='too many digits'
2: msg='some digits set to 0'
3: msg='too many digits and some digits set to 0'
endcase

return,{mantissa:mantissa,exp:exp}
end;function basetodecimal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function floattobinary,f
; Binary == string with 1's and 0's
if size(f,/type) eq 5 then nbytes=8 else nbytes=4
return,string(reverse(byte(f,0,nbytes)),format='('+string(nbytes,format='(I0)')+'b08)')
end;function floattobinary,f
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function decomposefloat,fin

; USAGE: decompose floating-point to sign, mantissa, base and exponent
; 
;    f=125848.516416
;    a=decomposefloat(f)
;    if a.s*a.mantissa*a.base^a.exp eq f then print,'DECOMPOSITION: OK'
;    
;    f=5546.6468d
;    a=decomposefloat(f)
;    if a.s*a.mantissa*a.base^a.exp eq f then print,'DECOMPOSITION: OK'
;
;
; Floating point:
;                    Sign        Exponent    Mantissa    Bias
; Half Precision    1 [15]        5 [14-10]    10 [09-00]    15
; Single Precision    1 [31]        8 [30-23]    23 [22-00]    127
; Double Precision    1 [63]       11 [62-52]    52 [51-00]    1023

; f = sign.fraction.base^(exponent)
;   = sign.mantissa.base^(exponent+exponent_fraction)
;
;     - sign: 1 => negative number; 0 => positive number
;    - base(radix): 2
;    - fraction: a fractional number with a specific base and nfraction digits: e.g. 10.23 (base=10, nfraction=4)
;            fraction = mantissa.base^exponent_fraction
;    - mantissa: the integer part of the fraction (all digits without the radix point)
;            normalized representation: take radix point just after the first non-zero digit
;                                                => exponent_fraction = -nfraction+1
;            normalized representation (radix-2): first non-zero digit is always 1 so we will omit this
;                                                => exponent_fraction = -nfraction
;    - exponent: subtract bias from the exponent stored in the machine to see the real exponent
;    - significand: other name for mantissa
; f(binary) = sign | exponent + bias | mantissa without leading 1

; MACHAR fields
;    - ibeta: base (e.g. 2)
;    - it: number of bits for the mantissa (including the omitted bit: nfraction=it-1 )
;    - iexp: number if bits for the exponent
;    - machep: the smallest n for which (1+base^n) ne 1 (machep=-nfraction)
;        - eps: base^machep
;    - negep: the smallest n for which (1-base^n) ne 1
;        - epsneg: base^negep => this is bigger than eps, so use eps
;    - minexp: minimum exponent (take normalized representation into account)
;        - xmin: base^minexp
;    - maxexp: maximal exponent
;        - xmax: (1-epsneg).base^maxexp

; Rounding errors when storing floating point numbers:
;    freal -> fmach = sign.mantissa.base^(exponent+machep) = sign.mantissa.eps.base^exponent
;    1 ulp (unit in last place)     = this is the distance between a stored floating point f1
;                                  number and its closest neighbour f2:
;                                              - exponent1 = exponent2 = exp
;                                              - mantissa1 = mantissa2 + 1
;                                = abs(f1 - f2)
;                                = eps.(mantissa1.base^exponent1-mantissa2.base^exponent2)
;                                = eps.base^(exp)
;    Absolute error when storing: abs(freal-fmach)
;                        - truncation: abs(freal-fmach) <= eps.base^(exp)
;                        - rounding: abs(freal-fmach) <= 0.5.eps.base^(exp)
;                        abs(freal-fmach) <= c.eps.base^(exp) = c ulp
;    Relative error when storing: abs(freal-fmach)/abs(freal)
;                        abserr  <= c.eps.base^(exp)/abs(freal)
;                                 = [c.eps.base^(exp)]/[(m.eps+/-c.eps).base^(exp)]
;                                 = c.eps/[m.eps+/-c.eps]
;                                 = c.eps/[1.xxxxxx(y+/-c)]
;                                         for truncation and radix-2, these are the extremes:
;                                             1.1111111 + 1 =10.0000000  >1
;                                             1.1111111 - 1 = 1.1111110  >1
;                                             1.0000000 + 1 = 1.0000001  >1
;                                             1.0000000 - 1 = 0.1111111  <1  problem?
;                                <= c.eps
;
; Two numbers are equal within the sorting error if
;    abs(f1-f2)/(abs(f1)>abs(f2)) <= eps   (divide by the biggest number to )

f=fin
type=size(f,/type)
if type eq 6 then begin
    re=float(f)
    im=imaginary(f)
    return,[[decomposefloat(re)],[decomposefloat(im)]]
endif else if type eq 9 then begin
    re=double(f)
    im=imaginary(f)
    return,[[decomposefloat(re)],[decomposefloat(im)]]
endif else if type ne 4 and type ne 5 then $
    if type eq 14 or type eq 15 then begin
        f=double(f)
        type=5
    endif else begin
        f=float(f)
        type=4
    endelse

double=type eq 5
res=machar(double=double)

machine=floattobinary(f)
base=res.ibeta
exp=basetodecimal(strmid(machine,1,res.iexp),base)
exp=exp.mantissa-res.maxexp+1 ; bias
fraction=basetodecimal('1.'+strmid(machine,res.iexp+1,res.it-1),base)
exp+=fraction.exp
mantissa=fraction.mantissa
s=(strmid(machine,0,1) eq '1')?-1:1

; f=s*mantissa*base^exp
return,{s:s,mantissa:mantissa,base:double?double(base):float(base),exp:exp}
end;function decomposefloat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function binarytofloat,integer,precision=precision
; Binary == integer that contains the 1's and 0's

; USAGE: convert binary to floating-point
; 
;    integer='3555'x
;    f=binarytofloat(integer,precision=0)
;    print,f,' (should be 0.33325...)'
;    
;    integer='3eaaaaab'xl
;    f=binarytofloat(integer,precision=1)
;    print,f,' (should be 0.33333...)'
;    
;    integer='3fd5555555555555'xll
;    f=binarytofloat(integer,precision=2)
;    print,f,' (should be 0.33333...)'


if n_elements(precision) eq 0 then precision=1
case precision of
1:     begin
    res=machar() ;single
    nbit=32
    double=0b
    endcase
2:     begin
    res=machar(/double) ;double
    nbit=64
    double=1b
    endcase
3:     begin
    res={ibeta:2l,it:11,iexp:5l,maxexp:16l}; half
    nbit=16
    double=0b
    endcase
else: return,0
endcase

machine=string(integer,format='(b0'+string(nbit,format='(I0)')+')')

base=res.ibeta
exp=basetodecimal(strmid(machine,1,res.iexp),base)
exp=exp.mantissa-res.maxexp+1 ; bias
fraction=basetodecimal('1.'+strmid(machine,res.iexp+1,res.it-1),base)
exp+=fraction.exp
mantissa=fraction.mantissa
s=(strmid(machine,0,1) eq '1')?-1:1

return,s*mantissa*(double?double(base):float(base))^exp
end;function binarytofloat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function floatpres,f,double=double,DWARF=DWARF,GIANT=GIANT,rDWARF=rDWARF,rGIANT=rGIANT

if n_elements(f) ne 0 then begin
    double=size(f,/type)
    double=double eq 5 or double eq 9
endif

res = machar(double=double)
DWARF=res.xmin
GIANT=res.xmax
rDWARF=sqrt(DWARF*1.5) * 10
rGIANT=sqrt(GIANT) * 0.1

return, res.eps
end;function floatpres
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function floatequal,x,y,epsm=epsm,array=array

; check type
type=size(x,/type)
float = type eq 4 or type eq 6
double = type eq 5 or type eq 9
complex = type eq 6 or type eq 9
type=size(y,/type)
float or= type eq 4 or type eq 6
double and= type eq 5 or type eq 9
complex or= type eq 6 or type eq 9

; handle complex
if complex then begin
    re_x=double?double(x):float(x)
    im_x=imaginary(x)
    re_y=double?double(y):float(y)
    im_y=imaginary(y)
    return,floatequal(re_x,re_y,epsm=epsm) and floatequal(im_x,im_y,epsm=epsm)
endif else begin
    if ~float and ~double then begin
        if keyword_set(array) then return,x eq y $
        else return,array_equal(x,y)
    endif
endelse

; threshold
if not keyword_set(epsm) then epsm=1
threshold=(machar(double=double)).eps*epsm
inf=double?double('inf'):float('inf')
zero=double?0d:0.

; relative error
relerr=abs(x-y)/(abs(x)>abs(y))

; both inf: equal
binf = finite(x,/inf,sign=1) and finite(y,/inf,sign=1)
binf or= finite(x,/inf,sign=-1) and finite(y,/inf,sign=-1)
ind=where(binf,ct)
if ct ne 0 then relerr[ind]=zero

; both zero: equal
bzero=(abs(x) le threshold) and (abs(y) le threshold)
ind=where(bzero,ct)
if ct ne 0 then relerr[ind]=zero

; check nan: not equal
ind=where(finite(relerr,/nan),ct)
if ct ne 0 then relerr[ind]=inf

if keyword_set(array) then return,relerr le threshold $
else return,total(relerr gt threshold,/int) eq 0
end;function floatequal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function floatequalangle,x,y,epsm=epsm

; Map x between 0 and 360
type=size(x,/type)
floatx = type eq 4 or type eq 6
if floatx then pi=!pi else pi=!dpi
x2 = x mod (2*pi)
ind=where(x2 lt 0,ct)
x2[ind] += 2*pi

; Map y between 0 and 360
type=size(y,/type)
floaty = type eq 4 or type eq 6
if floaty then pi=!pi else pi=!dpi
y2 = y mod (2*pi)
ind=where(y2 lt 0,ct)
y2[ind] += 2*pi

; Difference in one direction
diff = (x2 > y2) - (x2 < y2)

; Difference in the other direction
if floatx or floaty then pi=!pi else pi=!dpi
diff <= 2*pi - diff
if floatx or floaty then diff=float(diff)

return,floatequal(diff,0d,epsm=epsm)
end;function floatequalangle
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TypeConvert,data,type
if size(data,/type) eq 10 then begin
    case type of
        1:    (*data)=byte(0>round(*data)<255)
        2:    (*data)=fix(-32768>round(*data)<32767)
        3:    (*data)=long(-2147483648>round(*data)<2147483647)
        4:    (*data)=float(*data)
        5:    (*data)=double(*data)
        6:    (*data)=complex(*data)
        7:    (*data)=string(*data)
        9:    (*data)=dcomplex(*data)
        12:    (*data)=uint(0>round(*data)<65535)
        13:    (*data)=ulong(0>round(*data,/l64)<4294967295)
        14:    (*data)=long64(-9223372036854775808>round(*data,/l64)<9223372036854775807)
        15:    (*data)=ulong64(*data)
        else:
    endcase
    return,data
endif else begin
    case type of
        1:    return,byte(0>round(data)<255)
        2:    return,fix(-32768>round(data)<32767)
        3:    return,long(-2147483648>round(data)<2147483647)
        4:    return,float(data)
        5:    return,double(data)
        6:    return,complex(data)
        7:    return,string(data)
        9:    return,dcomplex(data)
        12:    return,uint(0>round(data)<65535)
        13:    return,ulong(0>round(data,/l64)<4294967295)
        14:    return,long64(-9223372036854775808>round(data,/l64)<9223372036854775807)
        15:    return,ulong64(data)
        else: return,data
    endcase
endelse
end;function TypeConvert
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SetType,A,type
CATCH, Error_status
IF Error_status NE 0 THEN return

case type of
1: A=byte(A)
2: A=fix(A)
3: A=long(A)
4: A=float(A)
5: A=double(A)
6: A=complex(A)
7: A=string(A)
9: A=dcomplex(A)
12: A=uint(A)
13: A=ulong(A)
14: A=long64(A)
15: A=ulong64(A)
else:
endcase
end;pro SetType
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Heap_Copy,p
; Make sure p is already a copy
; No objects should be involved

s=size(p,/type)
case s of
8:    begin ;structure
    for i=0l,n_tags(p)-1 do begin
        ; Prevent copying large array's of data
        s=size(p.(i),/type)
        if s eq 8 or s eq 10 then begin
            tmp=p.(i)
            Heap_Copy,tmp
            p.(i)=temporary(tmp)
        endif
    endfor
    endcase
10:    begin ;pointer
    ind=where(ptr_valid(p),ct)
    for i=0l,ct-1 do begin
        j=ind[i]
        p[j]=ptr_new(*p[j])
        Heap_Copy,*p[j]
    endfor
    endcase
else:
endcase
end;pro Heap_Copy
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function openw_safe,lun,path,_REF_EXTRA=extra
ON_IOERROR, NULL

; return 1 when error and 0 when ok
openw,lun,path,/get_lun,error=err,_EXTRA=extra
;if err eq 0 then begin
;    info=fstat(lun)
;    print,lun,' :',info.name
;endif
return,err ne 0
end;function openw_safe
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function openr_safe,lun,path,_REF_EXTRA=extra
ON_IOERROR, NULL

openr,lun,path,/get_lun,error=err,_EXTRA=extra
;if err eq 0 then begin
;    info=fstat(lun)
;    print,lun,' :',info.name
;endif
;message,!ERR_STRING,/info
return,err ne 0
end;function openr_safe
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function openu_safe,lun,path,_REF_EXTRA=extra
ON_IOERROR, NULL

openu,lun,path,/get_lun,error=err,_EXTRA=extra
;if err eq 0 then begin
;    info=fstat(lun)
;    print,lun,' :',info.name
;endif
return,err ne 0
end;function openu_safe
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro WRITE_TIFF_SAFE, file, image, error=error, _REF_EXTRA=extra

error=1

; Set type-keywords
s=size(image,/type)
case s of
1:  begin
    type='BYTE'
    endcase
2:  begin
    type='16bit INT'
    short=1
    signed=1
    endcase
3:  begin
    type='32bit INT'
    long=1
    signed=1
    endcase
4:  begin
    type='32bit FLOAT'
    float=1
    endcase
5:  begin
    type='64bit FLOAT'
    double=1
    endcase
6:  begin
    type='32bit COMPLEX'
    complex=1
    endcase
9:  begin
    type='64bit COMPLEX'
    dcomplex=1
    endcase
12: begin
    type='16bit UINT'
    short=1
    endcase
13: begin
    type='32bit UINT'
    long=1
    endcase
14: begin
    type='64bit INT'
    L64=1
    signed=1
    endcase
15: begin
    type='64bit UINT'
    L64=1
    endcase
else: begin
    type='unknown'
    return
    endelse
endcase

ct=0
CATCH, Error_status
IF Error_status NE 0 THEN begin
    if ct eq 1 then return
    ct++
    WRITE_TIFF,file,reform(Image,n_elements(image),1),short=short,long=long,float=float,double=double,$
        complex=complex,dcomplex=dcomplex,L64=L64,signed=signed,_EXTRA=extra
endif

WRITE_TIFF,file,image,short=short,long=long,float=float,double=double,$
    complex=complex,dcomplex=dcomplex,L64=L64,signed=signed,_EXTRA=extra
    
error=0
end;pro WRITE_TIFF_SAFE
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Binwidths,x
n=n_elements(x)
if n le 1 then return,0.
dx=x[1:n-1]-x[0:n-2]
dx=[dx[0],dx,dx[n-2]]
return,(dx[0:n-1]+dx[1:n])/2.
end;function Binwidths
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ConvertTime,ETIME,timeaf
timeaf=' sec'
if ETIME gt 60 then begin
   ETIME/=60.
   timeaf=' min'
   if ETIME gt 60 then begin
     ETIME/=60.
     timeaf=' hrs'
     if ETIME gt 24 then begin
         ETIME/=24.
        timeaf=' days'
        if ETIME gt 365 then begin
             ETIME/=365.
            timeaf=' years'
         endif
     endif
   endif
endif
end;pro ConvertTime
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ConvertSize,s,unit,bits=bits,HD=HD
; s: number of bits or bytes
; HD: harddrive size convention

if keyword_set(bits) then begin
    unit=' bits'
    if s gt 8 then begin
        s/=8.
        unit=' bytes'
    endif else return
endif else unit=' bytes'

if keyword_set(HD) then begin
    i=1000
    f=1000.
endif else begin
    i=1024
    f=1024.
endelse

if s ge i then begin
    s/=f
    unit=' kB'
    if s ge i then begin
        s/=f
        unit=' MB'
        if s ge i then begin
            s/=f
            unit=' GB'
            if s ge i then begin
                s/=f
                unit=' TB'
                if s ge i then begin
                    s/=f
                    unit=' PB'
                    if s ge i then begin
                        s/=f
                        unit=' EB'
                        if s ge i then begin
                            s/=f
                            unit=' ZB'
                            if s ge i then begin
                                s/=f
                                unit=' YB'
                            endif
                        endif
                    endif
                endif
            endif
        endif
    endif
endif

end;pro ConvertSize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function convertCRange,range,log=log
if ~keyword_set(log) then log=0
case log of
2:range=sqrt(range)
else:
endcase

return,range
end;function convertCRange
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GetCRange,ROIx=ROIx,ROIy=ROIy,ROIz=ROIz
if arg_present(ROIx) then begin
    ROIx=!X.CRange
    if !X.type eq 1 then ROIx=10.^ROIx
    ROIx=ROIx[sort(ROIx)]
endif
if arg_present(ROIy) then begin
    ROIy=!Y.CRange
    if !Y.type eq 1 then ROIy=10.^ROIy
    ROIy=ROIy[sort(ROIy)]
endif
if arg_present(ROIz) then begin
    ROIz=!Z.CRange
    if !Z.type eq 1 then ROIz=10.^ROIz
    ROIz=ROIz[sort(ROIz)]
endif
end;pro GetCRange
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakeRange, b,e,inc,incper=incper,rangeext=rangeext
; b: begin
; e: end
; inc: increment
; incper: increment is incper% of initial range
; rangeext: expand range with rangeext%

if b eq e then return,b

if n_elements(incper) ne 0 then inc=(e-b)*incper/100.
if n_elements(rangeext) ne 0 then begin
    rangeext=0>rangeext
    e=e+(e-b)*rangeext/100.
    b=b-(e-b)*rangeext/100.
endif

abin=ceil((e-b)/float(inc))+1
inc2=(e-b)/float((abin-1)>1)

return,inc2*lindgen(abin)+b
end; function MakeRange
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function msymbols,name,widget=widget,wfont=widgetfont

; A. Direct graphics fonts
;         1. Normal Device (!p.font=-1): Hershey fonts
;         2. PS Device with /isolatin (!p.font=0): Postscript fonts (i.e. a hardware or device font)
;     Embedded formatting is supported (introduced by the exclamation mark, !)
; 
; B. Widget fonts
;         1. set font keyword in the widget creation: the font name is different for each platform!!
;             b = Widget_Button(tlb, Value=msymbols('twotheta',/widget,wfont=wfont),font=wfont)
;     Embedded formatting is not supported
; 
; The default font !3 looks the same for all fonts

bPS=!d.name eq 'PS' and !p.font eq 0
bwidget=keyword_set(widget) ; Font 

case name of
'angstroms':     return,string(197b)
'degrees':         return,string(176b)
'plus_minus':     return,string(177b)
'micro':        return,string(181b)
'^1':             return,string(185b)
'^2':             return,string(178b)
'^3':             return,string(179b)
'twotheta':        begin
                if bwidget then begin
                    widgetfont='symbol'
                    return,'2q'
                endif else return,bPS?'!32!9q!X':'!32!4h!X'
                endcase
'Delta':        begin
                if bwidget then begin
                    widgetfont='symbol'
                    return,'D'
                endif else return,bPS?'!9D!X':'!7D!X'
                endcase
'lambda':        begin
                if bwidget then begin
                    widgetfont='symbol'
                    return,'lambda'
                endif else return,bPS?'!9l!X':'!7k!X'
                endcase
'nu':            begin
                if bwidget then begin
                    widgetfont='symbol'
                    return,'v'
                endif else return,bPS?'!9n!X':'!4m!X'
                endcase
else:            return,''
endcase

end;function msymbols
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DimSize,img,dim
if n_elements(dim) eq 0 then return,size(img,/dim)
sout=replicate(1l,dim)
s=size(img)
if s[0] eq 0 then sout[*]=1 else sout[0]=s[1:s[0]<dim]>1
return,sout
end;function DimSize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DimSizeFull,img,dim
sout=replicate(1l,dim)
s=size(img)
if s[0] eq 0 then begin
    sout[*]=1
    return,[dim,sout,s[2:*]]
endif else begin
    e=s[0]<dim
    sout[0]=s[1:e]>1
    return,[dim,sout,s[e+1:*]]
endelse
end;function DimSizeFull
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rebins,array,d1,d2,d3,d4,d5,d6,d7,d8
nd=N_PARAMS()-1
if nd le 1 then return,array
s=DimSize(array,nd)

case nd of
1: return,array[rebin(lindgen(s),d1,/sample)]
2: return,array[rebin(lindgen(s),d1,d2,/sample)]
3: return,array[rebin(lindgen(s),d1,d2,d3,/sample)]
4: return,array[rebin(lindgen(s),d1,d2,d3,d4,/sample)]
5: return,array[rebin(lindgen(s),d1,d2,d3,d4,d5,/sample)]
6: return,array[rebin(lindgen(s),d1,d2,d3,d4,d5,d6,/sample)]
7: return,array[rebin(lindgen(s),d1,d2,d3,d4,d5,d6,d7,/sample)]
8: return,array[rebin(lindgen(s),d1,d2,d3,d4,d5,d6,d7,d8,/sample)]
endcase

end;function rebins
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;function dimEXTRAC, Array, P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15
;ndim=(n_params()-1)/2
;s=size(Array)
;if s[0] lt ndim then Array=reform(Array,Dimsize(Array,ndim),/overwrite)
;return,EXTRAC(Array, P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15)
;end;function dimEXTRAC
;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function subarray,array,dims

; more a less the same as extrac

; Check new size
ndims=n_elements(dims)
if ndims eq 0 or ndims gt 8 then return,size(dummy)

; Check old size
if n_elements(array) eq 0 then array=0
s=size([array],/dim)
if n_elements(s) lt ndims then s=[s,lonarr(ndims-n_elements(s))]
dims>=1

; Shrink array
switch ndims of
8:    if dims[7] lt s[7] then array=array[*,*,*,*,*,*,*,0:dims[7]-1]
7:    if dims[6] lt s[6] then array=array[*,*,*,*,*,*,0:dims[6]-1,*]
6:    if dims[5] lt s[5] then array=array[*,*,*,*,*,0:dims[5]-1,*,*]
5:    if dims[4] lt s[4] then array=array[*,*,*,*,0:dims[4]-1,*,*,*]
4:    if dims[3] lt s[3] then array=array[*,*,*,0:dims[3]-1,*,*,*,*]
3:    if dims[2] lt s[2] then array=array[*,*,0:dims[2]-1,*,*,*,*,*]
2:    if dims[1] lt s[1] then array=array[*,0:dims[1]-1,*,*,*,*,*,*]
1:    if dims[0] lt s[0] then array=array[0:dims[0]-1,*,*,*,*,*,*,*]
endswitch

; Expand array
temp=temporary(array)
array=make_array(dimension=dims,type=size(temp,/type))

case ndims of
1:    array[0]=temporary(temp)
2:    array[0,0]=temporary(temp)
3:    array[0,0,0]=temporary(temp)
4:    array[0,0,0,0]=temporary(temp)
5:    array[0,0,0,0,0]=temporary(temp)
6:    array[0,0,0,0,0,0]=temporary(temp)
7:    array[0,0,0,0,0,0,0]=temporary(temp)
8:    array[0,0,0,0,0,0,0,0]=temporary(temp)
endcase

; Return size
return,size(array)
end;function subarray
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function value_locater,arr,val

n=n_elements(val)
ret=lonarr(n)
for i=0l,n-1 do begin
    temp=min(arr-val[i],p,/absolute)
    ret[i]=p
endfor

return,ret

;n=n_elements(arr)
;
;if arr[0] eq arr[n-1] and n_elements(val) eq 2 then return,[0,n-1] else begin
;    if arr[1] lt arr[0] then $
;        return,(n-1)-(value_locate(reverse(arr),val)>0) $
;    else return,value_locate(arr,val)
;endelse
end;function value_locater
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function chunklindgen,n
; Chunk lindgen:
;  input:    [1, 2, 1, 4]
;  output:   [0, 0, 1, 0, 0, 1, 2, 3]

nn=n_elements(n)
if nn eq 0 then return,-1l

nres=total(n,/int)
if nres eq 0 then return,-1l

res=lonarr(nres)
indexarr=findgen(max(n))
i2=0l
FOR i=0l,n_elements(n)-1 DO BEGIN
   IF n[i] GE 1 THEN BEGIN
      res[i2:i2+n[i]-1]=indexarr[0:n[i]-1]
      i2+=n[i]
   ENDIF
ENDFOR

return,res
end;function chunklindgen
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function chunkindex,n
; Chunk indexing:
;  input:  [1, 2, 1, 4]
;  output: [0, 1, 1, 2, 3, 3, 3, 3]

nn=n_elements(n)
if nn eq 0 then return,-1

CATCH, Error_status
IF Error_status NE 0 THEN return,-1 ; this happens all n are zero

if nn eq 1 then return,lonarr(n)

h=histogram(total(n>0,/CUMULATIVE,/int)-1,/BINSIZE,MIN=0,REVERSE_INDICES=ri)
return,ri[0:n_elements(h)-1]-ri[0]
end;function chunkindex
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function chunktotal,x,n,_REF_EXTRA=extra
; Chunk total:
;  input:  x=[1, 2, 1, 4]
;          n=[2,2]
;  output: [3,5]
;  
; (Uses chunk indexing)

nn=n_elements(n)
if nn eq 0 then return,0

CATCH, Error_status
IF Error_status NE 0 THEN return,lonarr(nn) ; this happens all n are zero

nx=n_elements(x)
if nx eq 0 then begin
    tmp=where(n ne 0,nn)
    if nn eq 0 then return,0l else return,lonarr(nn)
endif

if nn eq 1 then return,total(x[0:(n<nx)-1],_EXTRA=extra)

h=histogram(total(n>0,/CUMULATIVE,/int)-1,/BINSIZE,MIN=0,REVERSE_INDICES=ri)

sum=replicate(x[0]*0,nx,nn)
nh=n_elements(h)
sum[lindgen(nh),ri[0:nh-1]-ri[0]]=1
sum*=rebin(x,nx,nn,/sample)

return,total(sum,1,_EXTRA=extra)
end;function chunktotal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Drizzle,ind,val
; input: ind=[1,0,0,1]
;        val=[2,3,1,5]
; output: ind=[0,1]
;         val=[4,7]

indout=ind[0]
valout=val[0]

h1=histogram(ind,reverse_indices=ri1,OMIN=om)
h2=histogram(h1,reverse_indices=ri2,MIN=1)

if h2[0] ne 0 then begin 
    ;single values w/o duplication
    vec_inds=ri2[ri2[0]:ri2[1]-1] 
    indout=[indout,om+vec_inds]
    valout=[valout,val[ri1[ri1[vec_inds]]]]
endif
for j=1,n_elements(h2)-1 do begin 
    if h2[j] eq 0 then continue ;none with that many duplicates
    vec_inds=ri2[ri2[j]:ri2[j+1]-1] ;indices into h1
    indout=[indout,om+vec_inds]
    vec_inds=rebin(ri1[vec_inds],h2[j],j+1,/SAMPLE)+ $
             rebin(transpose(lindgen(j+1)),h2[j],j+1,/SAMPLE)
    valout=[valout,total(val[ri1[vec_inds]],2)]
endfor

ind=indout[1:*]
val=valout[1:*]
end;pro Drizzle
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SetIntersection, a, b, count=count, indices=indices, nd=nd, secondind=secondind, uniq=uniq, nafter=nafter
; Return: elements in A that are also in B
; /indices: return indices in A instead of the elements itself
; /nd: each row is a multi-dimensional point. "indices" is set automatically.
; /secondind: calculate indices in b and return structure
; /uniq: give only the first occuring index when double elements are present. Works only with /indices.
; nafter: number of positions after the comma that have to be taken into account. Not give: consider integers.

count=0
; Only need intersection of ranges
bnd=keyword_set(nd)
if bnd then begin
    ; Check dimensions
    sa=size(a)
    sb=size(b)
    if (sa[1] ne sb[1]) $
        or (sa[1] gt 8) or (sb[1] gt 8) then return,-1

    dimm=2
    dimm=2
    if sa[0] ne 2 then begin
        sa=[sa,1]
        a=reform(a,sa[1],1)
    endif
    if sb[0] ne 2 then begin
        sb=[sb,1]
        b=reform(b,sb[1],1)
    endif
endif
mina = Min(a, Max=maxa, DIMENSION=dimm)
minb = Min(b, Max=maxb, DIMENSION=dimm)
minab = mina > minb
maxab = maxa < maxb

; If either set is empty, or their ranges don't intersect:
IF ~array_equal(((maxa LT minab) AND (minb GT maxab)) OR ((maxb LT minab) AND (mina GT maxab)),0) THEN RETURN, -1

; Intersection of histograms
if keyword_set(nafter) then begin
    ; never define BINSIZE together with NBINS
    _BINSIZE=10.^(-fix(nafter))
    
    ; enlarge range to avoid rounding errors
    minab-=_BINSIZE
    maxab+=_BINSIZE
    
    NBINS=round((maxab-minab)/_BINSIZE+1)
endif else begin
    BINSIZE=1
    _BINSIZE=1
endelse
if bnd then begin
    ; Make histograms
    ha=hist_nd(a, Min=minab, Max=maxab,BINSIZE,NBINS=NBINS, REVERSE_INDICES=ra)
    if keyword_set(secondind) then hb=hist_nd(b, Min=minab, Max=maxab,BINSIZE,NBINS=NBINS, REVERSE_INDICES=rb) $
    else hb=hist_nd(b, Min=minab, Max=maxab,BINSIZE,NBINS=NBINS)

    indices=1
endif else begin
    ; Make histograms
    ha=Histogram(a, Min=minab, Max=maxab,BINSIZE=BINSIZE,NBINS=NBINS, REVERSE_INDICES=ra)
    if keyword_set(secondind) then hb=Histogram(b, Min=minab, Max=maxab,BINSIZE=BINSIZE,NBINS=NBINS, REVERSE_INDICES=rb) $
    else hb=Histogram(b, Min=minab, Max=maxab,BINSIZE=BINSIZE,NBINS=NBINS)

endelse

r = Where((ha NE 0) AND (hb NE 0), count)
IF count EQ 0 THEN RETURN, -1

if keyword_set(indices) then begin
    uniq=keyword_set(uniq)
    if keyword_set(secondind) then begin
        if uniq then begin
            r={a:ra[ra[r]],b:rb[rb[r]]}
        endif else begin
            buffera=-1
            bufferb=-1
            for i=0l,count-1 do begin
                buffera=[buffera,ra[ra[r[i]]:ra[r[i]+1]-1]]
                bufferb=[bufferb,rb[rb[r[i]]:rb[r[i]+1]-1]]
            endfor
            r={a:buffera[1:*],b:bufferb[1:*]}
        endelse
    endif else begin
        if uniq then begin
            r=ra[ra[r]]
        endif else begin
            buffer=-1
            for i=0l,count-1 do buffer=[buffer,ra[ra[r[i]] : ra[r[i]+1]-1]]
            r=buffer[1:*]
        endelse
    endelse
endif else begin
    r *= _BINSIZE
    r += minab
endelse

return,r

end;function SetIntersection
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SetSysVar,p
!X=(*p).x
!Y=(*p).y
!P=(*p).p
end;pro SetSysVar
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GetSysVar,p
(*p).x=!X
(*p).y=!Y
(*p).p=!P
end;pro GetSysVar
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PointSysVar
return,ptr_new({x:!X,y:!Y,p:!P})
end;function PointSysVar
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro delvar2, varname
IF (N_Elements(varname) NE 0) THEN tempvar = Temporary(varname)
end;pro delvar2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro delstruct,keywords,struc
if n_elements(struc) ne 0 then begin
    tags=tag_names(struc)
    bdel=bytarr(n_tags(struc))
    for i=0l,n_elements(keywords)-1 do $
        bdel or= tags eq strupcase(keywords[i])
    ind=where(~bdel,ct)
    if ct eq 0 then delvar2,struc $
    else begin
        j=ind[0]
        struc2=create_struct(tags[j],struc.(j))
        for i=1,ct-1 do begin
            j=ind[i]
            struc2=create_struct(struc2,tags[j],struc.(j))
        endfor
        struc=struc2
    endelse
endif
end;pro delstruct
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro addstruct,tags,val,struc
for i=0l,n_elements(tags)-1 do begin
    if (N_ELEMENTS(struc) eq 0) then begin
        struc = CREATE_STRUCT(tags[i], val[i])
        continue
    endif
    if (where(tag_names(struc) eq strupcase(tags[i])))[0] ne -1 then continue
    struc = CREATE_STRUCT(struc, tags[i], val[i])
endfor
end;pro addstruct
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SRDirectGraph,struct
if n_tags(struct) eq 0 then begin
    struct={xs:!x.s,ys:!y.s}
endif else begin
    !x.s=struct.xs
    !y.s=struct.ys
endelse
end;pro SRDirectGraph
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro markdir,dir
j=strlen(dir)-1

for i=0l,n_elements(dir)-1 do begin
    j=strlen(dir[i])-1
    last=strmid(dir[i],j,1)
    while (last eq '\' or last eq '/') and j gt 0 do begin
        j--
        last=strmid(dir,j,1)
    endwhile
    
    dir[i]=strmid(dir[i],0,j+1)+path_sep()
endfor

end;pro markdir
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CutPath,pathii,path=path,file=file,ext=ext

;e.g. "D:\data\esrf0304\diff_scan53\hold\23-6_51814.mccd"
; path='D:\data\esrf0304\diff_scan53\hold\'
; file='23-6_51814'
; ext='.mccd'

; Extention must be the same for all paths
pathi = pathii
PPos = strpos(pathi[0], '.',/reverse_search)
ext = (PPos ne -1)?strmid(pathi[0],PPos):''

path = FILE_DIRNAME(pathi, /MARK_DIRECTORY)
file = FILE_BASENAME(pathi , ext)

n=strlen(pathi)-1
for i=0l,n_elements(pathi)-1 do $
if strmid(pathi[i],n[i]) eq path_sep() then begin
    path[i]=pathi
    file[i]=''
endif

return, 1B
end;function CutPath
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddRegistry,file,path,var,data
SPAWN,'REG ADD '+path+' /v '+var+' /d '+data,EXIT_STATUS=S
if S ne 0 then begin ; windows 2000 doesn't have reg.exe
    if openw_safe(lun,file) then return
    printf,lun,'Windows Registry Editor Version 5.00'
    printf,lun,'['+path+']'
    printf,lun,'"'+var+'"="'+data+'"'
    free_lun,lun
    SPAWN, 'regedit/s '+file
endif
end;pro AddRegistry
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro OpenURL,file,URL
if openw_safe(lun,file) then return
printf,lun,'<HTML><HEAD><META HTTP-EQUIV="refresh" CONTENT="0; URL='+URL+'"></HEAD>If this page does not refresh, click <A HREF="'+URL+'">here</A>.</HTML>'
free_lun,lun
file=file_search(file,/full)
online_help,book=file
end;pro OpenURL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_PromptNumber,ID

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
end;pro CleanUp_PromptNumber
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro PromptNumber_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

type=widget_info(ev.id,/type)
ID=widget_info(ev.top,FIND_BY_UNAME='text')
if (type eq 1) and (widget_info(ID,/valid_id)) then type=3

widget_control,ev.top,get_uvalue=list
case type of
1:    begin
    endcase
8:    begin
    ; Save dynamically
    (*list.filter)=ev.index
    return
    endcase
3:    begin
    ; Read from text widget
    widget_control,ID,get_value=str
    if strlen(str[0]) eq 0 then return
    ; Save dynamically
    (*list.filter)=str[0]
    endcase
endcase

; Return
widget_control,ev.top,/destroy

end;pro PromptNumber_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PromptNumber,in,top,label,xsize=xsize,choice=choice,username=username,userpro=userpro,userdata=userdata,modal=modal

;----Prompt for filter----
choice=keyword_set(choice)
filter=PTR_NEW(in)
if WIDGET_INFO(top,/VALID_ID) eq 1 then begin
    WIDGET_CONTROL, top, sensitive=0
    GROUP_LEADER=top
    modal=1
endif else begin
    modal=keyword_set(modal)
endelse

device,get_screen_size=screen
base=widget_base(title='',/row,uvalue={top:top,filter:filter,choice:choice},$
                GROUP_LEADER=GROUP_LEADER,modal=modal,xoffset=screen[0]/2.,yoffset=screen[1]/2.)
label=widget_label(base,value=label)
if choice then begin
    drop=widget_droplist(base,value=(*filter),uname='drop')
    *filter=0
endif else begin
    ;s=widget_info(base,STRING_SIZE=(*filter))
    ;text=widget_text(base,/editable,value=(*filter),uname='text',xsize=s[0]*1.1)
    if keyword_set(xsize) then text=widget_text(base,/editable,value=(*filter),uname='text',xsize=xsize<50) $
    else text=widget_text(base,/editable,value=(*filter),uname='text')
endelse

buser=keyword_set(userdata) and keyword_set(userpro) and keyword_set(username)
if buser then b2=widget_button(base,value=username,uvalue=userdata,EVENT_PRO=userpro)

b=widget_button(base,value='OK')

widget_control,base,/realize
widget_control,b,/INPUT_FOCUS
Xmanager,'PromptNumber',base, event_handler='PromptNumber_Event',$
    cleanup='CleanUp_PromptNumber',GROUP_LEADER=GROUP_LEADER

out=(*filter)
PTR_FREE,filter
return,out
end;function PromptNumber
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WeakScl,img,sclmax,power

ma=max(img,min=mi)

if power ne 1 then begin
    mi=sclmax[0]>mi
    ret=(img-mi)>0
    ret<=(ma<sclmax[1])-mi
    return,bytscl(ret^power)

    ;if mi lt 0 then return,bytscl((img-mi)^power,min=sclmax[0]>0,max=((ma-mi)^power)<sclmax[1]) $
    ;else return,bytscl(img^power,min=sclmax[0]>(mi^power),max=(ma^power)<sclmax[1])
endif else begin
    return,bytscl(img,min=sclmax[0]>mi,max=ma<sclmax[1])
endelse

end;function WeakScl
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function stringr,in
return,strtrim(in,2)
;return,strcompress(string(in),/remove_all)
end;function stringr
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakeNumber,nr,prompt=prompt

prompt=keyword_set(prompt)

; Fill with zero's or not?
if prompt then Result = DIALOG_MESSAGE('Fill with zero`s or not?',/question) else result='Yes'
if result eq 'Yes' then begin
    m=strlen(string(max(nr),format='(I0)'))
    if prompt then begin
        repeat begin
            Result = DIALOG_MESSAGE('Expand number to '+stringr(m)+' digits?',/question) 
            m2=long(PromptNumber(stringr(m),0L,'Enter number of digits to use:'))
        endrep until m2 ge m
        m=m2
    endif
    
    nr=string(nr,format='(I0'+string(m,format='(I0)')+')')
endif

return,nr
end;function MakeNumber
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ScaleArr,arrin,b,e

mi=min(arrin)
arr=arrin-mi
m=max(arr)
if (m ne 0) then out=arr/float(m) else out=arr

; arr is here scale between 0 and e-b and shifted with b
out=b+out*(e-b)

return,out
end;function ScaleArr
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wherer,arin,operator,elin,n,count

; inline:
; ind=where(round(arin*1000) eq round(elin*1000),count)
mult=10.^(fix(n)>0)
ar=round(mult*arin)
el=round(mult*elin)

case STRUPCASE(operator) of
'LT':ind=where(ar lt el,count)
'GT':ind=where(ar gt el,count)
'LE':ind=where(ar le el,count)
'GE':ind=where(ar ge el,count)
'EQ':ind=where(ar eq el,count)
'NE':ind=where(ar ne el,count)
endcase

return,ind

end;function wherer
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro parseoperator,operator,comp
operator=strcompress(operator,/remove_all)
if n_elements(comp) ne 0 then return

p=strpos(operator,'<=')
if p ne -1 then begin
    comp=float(strmid(operator,p+2))
    operator='LE'
    return
endif

p=strpos(operator,'>=')
if p ne -1 then begin
    comp=float(strmid(operator,p+2))
    operator='GE'
    return
endif

p=strpos(operator,'<')
if p ne -1 then begin
    comp=float(strmid(operator,p+1))
    operator='LT'
    return
endif

p=strpos(operator,'>')
if p ne -1 then begin
    comp=float(strmid(operator,p+1))
    operator='GT'
    return
endif

p=strpos(operator,'!=')
if p ne -1 then begin
    comp=float(strmid(operator,p+2))
    operator='NE'
    return
endif

p=strpos(operator,'==')
if p ne -1 then begin
    comp=float(strmid(operator,p+2))
    operator='EQ'
    return
endif

end;pro parseoperator
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compare,arr,operator,comp
parseoperator,operator,comp
if size(arr,/type) eq 10 then begin
case STRUPCASE(operator) of
'LT':return,*arr lt comp
'GT':return,*arr gt comp
'LE':return,*arr le comp
'GE':return,*arr ge comp
'EQ':return,*arr eq comp
'NE':return,*arr ne comp
else: return,make_array(size=size(*arr),/byte)
endcase
endif else begin
case STRUPCASE(operator) of
'LT':return,arr lt comp
'GT':return,arr gt comp
'LE':return,arr le comp
'GE':return,arr ge comp
'EQ':return,arr eq comp
'NE':return,arr ne comp
else: return,make_array(size=size(arr),/byte)
endcase
endelse
end;function compare
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro compareself,arr,operator,comp
parseoperator,operator,comp
if size(arr,/type) eq 10 then begin
case STRUPCASE(operator) of
'LT':*arr lt= comp
'GT':*arr gt= comp
'LE':*arr le= comp
'GE':*arr ge= comp
'EQ':*arr eq= comp
'NE':*arr ne= comp
else: 
endcase
endif else begin
case STRUPCASE(operator) of
'LT':arr lt= comp
'GT':arr gt= comp
'LE':arr le= comp
'GE':arr ge= comp
'EQ':arr eq= comp
'NE':arr ne= comp
else: 
endcase
endelse
end;pro compareself
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro setsens,top,strmat,sen=sen
n=n_elements(strmat)
if n eq 0 then return
if not keyword_set(sen) then sen=bytarr(n)+1
for i=0l,n-1 do begin
    ID1=widget_info(top,FIND_BY_UNAME=strmat[i])
    if ID1 ne 0 then widget_control,ID1,sensitive=1
endfor
end;pro setsens
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UniqValues,array
return,array[uniq(array,sort(array))]
end;function UniqValues
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ShrinkArray, array,ind_,ROW=ROW,count=ct
ind=ind_

; Init
null=size(array,/type)
bp=null eq 10

if bp then begin
    null=size(*array,/type)
    s=size(*array)
endif else begin
    s=size(array)
endelse

if null eq 7 then null='' else null=0

if (s[0] eq 0) or (s[0] gt 3) then begin
    ct=0
    return, (bp?array:null) ; scalar or higher then dim three
endif
if not keyword_set(ROW) then INDS=0 else INDS=ROW<(s[0]-1)

; Remove indices out of range
n=s[INDS+1]
ind2=where(ind ge 0 and ind lt n,ct)
if ct eq 0 then return,array
ind=ind[ind2]

; Indices to keep
b=bytarr(n)
b[ind]=1
index=where(~b,ct)
if ct eq 0 then begin
    if bp then begin
        ptr_free,array
        return,array
    endif else return,null
endif

; Return
if bp then begin
    case INDS of
    0: (*array)=(*array)[index,*,*] ;s[0]=1, 2 or 3
    1: (*array)=(*array)[*,index,*] ;s[0]=2 or 3
    2: (*array)=(*array)[*,*,index] ; s[0]=3
    endcase
    return,array
endif else begin
    case INDS of
    0: return,array[index,*,*] ;s[0]=1, 2 or 3
    1: return,array[*,index,*] ;s[0]=2 or 3
    2: return,array[*,*,index] ; s[0]=3
    endcase
endelse

end;function ShrinkArray
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MoveElementInArray,nel,indold,indnew
indold=0>indold<(nel-1)
indnew=0>indnew<(nel-1)

ind=[-1,lindgen(nel),-1]
ind=[ind[0:indold],ind[indold+2:*]]
ind=[ind[0:indnew],indold,ind[indnew+1:*]]
ind=ind[1:nel]

return,ind
end;function MoveElementInArray
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_ChooseFromList,ID
if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
WIDGET_CONTROL, list.top, sensitive=1
end;CleanUp_ChooseFromList
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ChooseFromList_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_value=val
if val eq 'OK' then begin
    widget_control,ev.top,/destroy
    return
endif

widget_control,ev.id,get_uvalue=uval
widget_control,ev.top,get_uvalue=list
(*list.bool)[uval]=ev.select

end;pro ChooseFromList_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ChooseFromList,str,title=title,nsel=nsel,top=top,select_one=select_one

n=n_elements(str)
if not keyword_set(title) then title=''

;----Prompt for filter----
bool=PTR_NEW(bytarr(n))
if WIDGET_INFO(top,/VALID_ID) eq 1 then WIDGET_CONTROL, top, sensitive=0
base=widget_base(title=title,/column,uvalue={top:top,bool:bool})
if keyword_set(select_one) then baset=widget_base(base,/column,/exclusive) $
else baset=widget_base(base,/column,/nonexclusive)
for i=0l,n-1 do button=widget_button(baset,value=str[i],uvalue=i)
b=widget_button(base,value='OK')
widget_control,base,/realize
widget_control,b,/INPUT_FOCUS
Xmanager,'ChooseFromList',base, event_handler='ChooseFromList_Event',$
    cleanup='CleanUp_ChooseFromList',GROUP_LEADER=top

sel=where(*bool eq 1,nsel)
PTR_FREE,bool

return,sel
end;function ChooseFromList
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro strreplace,str,char1,char2
;+
; :Description:
;    Replace characters in a string.
;-
str=string(255b)+str+string(255b)
for i=0l,n_elements(str)-1 do begin
    tmp=STRJOIN(STRSPLIT(str[i],char1, /EXTRACT), char2)
    str[i]=strmid(tmp,1,strlen(tmp)-2)
endfor
end;pro strreplace
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ResetDim_cleanup,ID
if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
if WIDGET_INFO(list.top,/VALID_ID) then WIDGET_CONTROL, list.top, sensitive=1
end;pro ResetDim_cleanup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ResetDim_event,ev
type=widget_info(ev.id,/type)
case type of
1:     widget_control,ev.top,/destroy
3:    begin
    widget_control,ev.top,get_uvalue=list
    widget_control,ev.id,get_value=dimi,get_uvalue=i
    dimi=ulong(dimi)
    i=long(i)
    
    dim=ulong64(*list.pdim)
    ndim=n_elements(dim)
    
    if i eq ndim-1 then j=i-1 else j=i+1
    n=dim[i]*dim[j]
    if n mod dimi eq 0 and dimi ne 0 then begin
        dim[i]=dimi
        dim[j]=n/dimi
    endif
    
    *list.pdim=dim
    for i=0l,ndim-1 do begin
        uval=stringr(i)
        ID=widget_info(ev.top,find_by_uname=uval)
        widget_control,ID,set_value=stringr(dim[i])
    endfor

    endcase
endcase

end;pro ResetDim_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ResetDim,dim,top,label
pdim=PTR_NEW(dim)
if n_elements(top) eq 0 then top=0L
if WIDGET_INFO(top,/VALID_ID) eq 1 then begin
    WIDGET_CONTROL, top, sensitive=0
    GROUP_LEADER=top
    modal=1
endif

device,get_screen_size=screen
base=widget_base(title='Edit dimensions',/column,uvalue={top:top,pdim:pdim},$
                GROUP_LEADER=GROUP_LEADER,modal=modal,xoffset=screen[0]/2.,yoffset=screen[1]/2.)

label=widget_label(base,value=label)
for i=0l,n_elements(dim)-1 do begin
    uval=stringr(i)
    text=widget_text(base,/editable,value=stringr(dim[i]),uname=uval,uvalue=uval)
endfor

b=widget_button(base,value='OK')
widget_control,base,/realize
widget_control,b,/INPUT_FOCUS
Xmanager,'ResetDim',CLEANUP='ResetDim_cleanup',base,GROUP_LEADER=GROUP_LEADER

out=*pdim
PTR_FREE,pdim
return,out
end;function ResetDim
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IndexStruct,struc,ind

if n_elements(ind) eq 0 then return,struc

n=n_tags(struc)
names=tag_names(struc)
namesadd=lonarr(n)

j=ind[0]
strucnew=create_struct(names[j],struc.(j))

for i=1l,n_elements(ind)-1 do begin
    j=ind[i]
    tmp=where(tag_names(strucnew) eq names[j],ct)
    if ct ne 0 then namesadd[j]++
    strucnew=create_struct(strucnew,names[j]+((namesadd[j] eq 0)?'':stringr(namesadd[j])),struc.(j))
endfor

return,strucnew
end;function IndexStruct
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
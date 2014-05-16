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

function ReadCCDInfo
return,{showheader:0b,bnorm:0b,int:'',gain:'',offset:'',binxs:2048,binys:2048,binlendian:1b,bintype:12,binnxpad3:0}
end;function ReadCCDInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro mWRITE_CSV,file,data,colnames=colnames,rownames=rownames,sep=sep,dec=dec
;+
; :Description:
;    WRITE_CSV with different seperator.
;-

if ~keyword_set(sep) then sep=','
if ~keyword_set(dec) then dec='.'
if sep eq '\t' then sep=string(9b)

s=size(data)
case s[0] of
1:    begin
        ncol=s[1]
        nrow=1
        endcase
2:    begin
        ncol=s[1]
        nrow=s[2]
        endcase
else: return
endcase

if openw_safe(lun,file) then return

; Column names
if keyword_set(colnames) then begin
    str=string(colnames,format='('+string(ncol,format='(I0)')+'(A,"'+sep+'"))')
    if keyword_set(rownames) then str=' '+sep+str
    if dec ne '.' then strreplace,str,'.',','
    printf,lun,str,format='(A)'
endif

; Data
str=string(data,format='('+string(ncol,format='(I0)')+'(F,"'+sep+'"))')
if dec ne '.' then strreplace,str,'.',dec

; Row names
if keyword_set(rownames) then str=rownames+sep+str

; Output
printf,lun,str,format='(A)'

free_lun,lun
end;pro mWRITE_CSV
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_tiff_tagvalue,type,nbytes=nbytes

case type of
1:     begin
    value=0b
    nbytes=1
    endcase
2:  begin
    value=0b
    nbytes=1
    endcase
3:  begin
    value=0u
    nbytes=2
    endcase
4:  begin
    value=0ul
    nbytes=4
    endcase
5:  begin
    value=[0ul,0ul]
    nbytes=8
    endcase
6:  begin
    value=0b
    nbytes=1
    endcase
7:  begin
    value=0b
    nbytes=1
    endcase
8:  begin
    value=0
    nbytes=2
    endcase
9:  begin
    value=0l
    nbytes=4
    endcase
10:  begin
    value=[0l,0l]
    nbytes=8
    endcase
11:  begin
    value=0.
    nbytes=4
    endcase
12:  begin
    value=0d
    nbytes=8
    endcase
endcase

return,value
end;function read_tiff_tagvalue
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro encodeLZW_initializestringtable,stringtable
stringtable=list()
stringtable.add,indgen(258),/extract
end;pro encodeLZW_initializestringtable
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro encodeLZW_writecode,stream,value,format
; format='(b09)', format='(b010)', format='(b011)' or format='(b012)'
stream+=strjoin(string(value,format=format))
end;pro encodeLZW_writecode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro encodeLZW_binarystream,streamin
; Convert string '10110000110110...' to binary

zero=(byte('0'))[0]

stream=streamin
n=strlen(stream)
nzeros=(8-(n mod 8)) mod 8
if nzeros ne 0 then stream+=string(make_array(nzeros,value=zero))
n+=nzeros
nbytes=n/8

stream=reform(byte(stream),8,nbytes)-zero
stream*=rebin(reverse(2b^bindgen(8)),8,nbytes,/sample)
stream=total(stream,1,/pres)
streamin=stream
end;pro encodeLZW_binarystream
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro encodeLZW_codesfromstream,streamin
; This is for debugging
zero=(byte('0'))[0]

stream=streamin
n=strlen(stream)
nzeros=(9-(n mod 9)) mod 9
if nzeros ne 0 then stream+=string(make_array(nzeros,value=zero))
n+=nzeros
nbytes=n/9

stream=reform(byte(stream),9,nbytes)-zero
stream*=rebin(reverse(2^indgen(9)),9,nbytes,/sample)
stream=total(stream,1,/pres)
streamin=stream
end;pro encodeLZW_codesfromstream
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro encodeLZW_tiff,strip
; Example:
; strip=byte([7,7,7,8,8,7,7,6,6]) & encodeLZW_tiff,strip
; Codes should be:
;  256   7 258   8   8 258   6   6 257

bitstream=''
;codestream=!null
encodeLZW_initializestringtable,stringtable
ClearCode=256
EndOfInformation=257
format='(b09)'
encodeLZW_writecode,bitstream,ClearCode,format
;codestream=[codestream,ClearCode]

Omega=!null
foreach K, strip, index do begin
    Omeganew=[Omega,K] ; new string
    tmp=stringtable.where(Omeganew,count=ct) ; New string already in table?
    if ct eq 0 then begin
        ; NO
        ; 
        ; Add index of the previous string to the bitstream:
        encodeLZW_writecode,bitstream,stringtable.where(Omega),format
;        codestream=[codestream,stringtable.where(Omega)]

        ; Add the new string to the table:
        stringtable.add,Omeganew
        
        ; Start new string
        Omega=K
        
        ; Increment code length in bits (not sure when exactly)
        n=n_elements(stringtable)
        if n eq 4095 then begin
            encodeLZW_writecode,bitstream,ClearCode,format
;            codestream=[codestream,ClearCode]
            encodeLZW_initializestringtable,stringtable
            format='(b09)'
            
;            print,n_elements(codestream)
;            print,codestream[-5:*],format='(6I5)'
        endif else if (n mod 512) eq 0 then begin
            if ispoweroftwo(n) then $
            case n of
            512: format='(b010)'
            1024:format='(b011)'
            2048:format='(b012)'
            endcase
            
;            print,n_elements(codestream)
;            print,codestream[-5:*],format='(6I5)'
        endif

    endif else Omega=Omeganew ; YES
endforeach

encodeLZW_writecode,bitstream,[stringtable.where(Omega),EndOfInformation],format
;codestream=[codestream,stringtable.where(Omega),EndOfInformation]

encodeLZW_binarystream,bitstream
;encodeLZW_codesfromstream,bitstream ; for debugging
strip=temporary(bitstream)

;print,n_elements(codestream)
;print,codestream[-5:*],format='(6I5)'

end;pro encodeLZW_tiff
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function decodeLZW_getnextcode,strip,byteoffset,bitoffset,ncodebits
; ncodebits between 9 and 12

nbits=bitoffset+ncodebits
nbytes=nbits/8+((nbits mod 8) ne 0)
byteoffsetnew=byteoffset+nbytes-1

bits=byte(strmid(string(strip[byteoffset:byteoffsetnew],format='('+stringr(nbytes)+'b08)'),bitoffset,ncodebits))-(byte('0'))[0]

bitoffset=nbits mod 8
byteoffset=byteoffsetnew + (bitoffset eq 0)

return,total(reverse(2^indgen(ncodebits))*bits,/pres)
end;function decodeLZW_getnextcode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro decodeLZW_tiff,strip
; http://www.fileformat.info/format/tiff/corion-lzw.htm
; Tested with LZW compressed images of IrfanView and Photoshop

bytestream=!null
;codestream=!null
ClearCode=256
EndOfInformation=257

byteoffset=0l ; in strip
bitoffset=0b ; in current byte (0,1,...,7)
ncodebits=9 ; number of bits in each code

code=decodeLZW_getnextcode(strip,byteoffset,bitoffset,ncodebits)
;codestream=[codestream,code]
while code ne EndOfInformation do begin
    if code eq ClearCode then begin
        encodeLZW_initializestringtable,stringtable
        ncodebits=9 ; number of bits in each code
        code=decodeLZW_getnextcode(strip,byteoffset,bitoffset,ncodebits)
;        codestream=[codestream,code]
        if code eq EndOfInformation then break
        bytestream=[bytestream,stringtable[code]]
        oldcode=code
    endif else begin
        if code lt n_elements(stringtable) then begin
            ; Code already in table
            bytestream=[bytestream,stringtable[code]]
            stringtable.add,[stringtable[oldcode],(stringtable[code])[0]]
            oldcode=code
        endif else begin
            ; Code not in table yet
            outstring=stringtable[oldcode]
            outstring=[outstring,outstring[0]]
            bytestream=[bytestream,outstring]
            stringtable.add,outstring
            oldcode=code
        endelse
    endelse
    
    ; Increment code length in bits (not sure when exactly)
    n=n_elements(stringtable)+1
    if (n mod 512) eq 0 then $
        if ispoweroftwo(n) then ncodebits++
    
    ; Next code
    code=decodeLZW_getnextcode(strip,byteoffset,bitoffset,ncodebits)
;    codestream=[codestream,code]
endwhile

strip=temporary(bytestream)

end;pro decodeLZW_tiff
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro reform_tiff_format,strip,bitspersample,bigendian,float=float,LZW=LZW

if keyword_set(LZW) then decodeLZW_tiff,strip

case bitspersample[0] of
8:    
12: begin
    bswap=(~isLittleEndian()) ne bigendian

    ; prepare bytes for conversion to 16bit
    nbytes=n_elements(strip)
    nval=nbytes*8/bitspersample[0]
    nadd=(6-(nbytes mod 6)) mod 6
    if nadd ne 0 then strip=[strip,bytarr(nadd)]
    nbytes+=nadd

    ;i=(where(strip ne 0))[0] & i=i/6*6 & print,strip[i:i+6*2-1],format='(6Z02)'

    ; convert to 16bit array
    strip=fix(strip,0,3,nbytes/6) ; BIG ENDIAN
    if ~bigendian then SWAP_ENDIAN_INPLACE, strip

    ; convert to 12bit array
    ; Each tripple of 16bit numbers: [4 4 4][4 \ 4 4][4 4 \ 4][4 4 4]
    strip=[ishft(strip[0,*],-4),$
        ishft(strip[0,*] and '000F'x,8) or ishft(strip[1,*] and fix('FF00'x),-8),$
        ishft(strip[1,*] and '00FF'x,4) or ishft(strip[2,*] and fix('F000'x),-12),$
        strip[2,*] and '0FFF'x]
    strip=strip[0:nval-1]
    
    ;print,strip[round((i)/1.5):round((i)/1.5)+5],format='(6Z04)'
    
    ;if ~bswap then SWAP_ENDIAN_INPLACE, strip ; ?????? Photonic Science (MicroXAS)
    ;strip=ishft(strip,-4)
    
    ;if ~bswap then begin
    ;    strip=ishft(strip and '000F'x,8) or (strip and '00F0'x) or ishft(strip and '0F00'x,-8)
    ;endif
    
    SWAP_ENDIAN_INPLACE, strip, SWAP_IF_BIG_ENDIAN=~bigendian, SWAP_IF_LITTLE_ENDIAN=bigendian
    
    endcase
16:    begin
    bytetotype,strip,12
    SWAP_ENDIAN_INPLACE, strip, SWAP_IF_BIG_ENDIAN=~bigendian, SWAP_IF_LITTLE_ENDIAN=bigendian
    endcase
32:    begin
    bytetotype,strip,keyword_set(float)?4:13
    SWAP_ENDIAN_INPLACE, strip, SWAP_IF_BIG_ENDIAN=~bigendian, SWAP_IF_LITTLE_ENDIAN=bigendian
    endcase
64:    begin
    bytetotype,strip,keyword_set(float)?5:15
    SWAP_ENDIAN_INPLACE, strip, SWAP_IF_BIG_ENDIAN=~bigendian, SWAP_IF_LITTLE_ENDIAN=bigendian
    endcase
endcase
end;pro reform_tiff_format
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_tiff_structure,file,error=error,format=format

;  TIFF format: http://www.fileformat.info/format/tiff/
;                  http://en.wikipedia.org/wiki/Tagged_Image_File_Format
;
;  Internal TIFF structure fields:
;        IFH  .byteorder
;             .tiffversion
;             .idflength (always 0)
;             .idfoffset
;        IFD0 .idflength
;             .tag0--------------->.code
;             .tag1                .type
;             ...                  .length
;             .idfoffset           .value (this is an fileptr offset when heap is specified)
;             .image              (.heap)
;        IFD1...

bigendian=0b
error=1b

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    print,'Error in reading ',file
    if n_elements(lun) ne 0 then free_lun,lun
    message,!ERR_STRING,/info
    return,0
ENDIF

; Open file with the correct byte order
openr,lun,file,/get_lun
byteorder=0
readu,lun,byteorder
bigendian=byteorder eq '4D4D'x ; 4949h (II: Intel) or 4D4Dh (MM: Motorola)
close,lun
openr,lun,file, SWAP_IF_BIG_ENDIAN=~bigendian, SWAP_IF_LITTLE_ENDIAN=bigendian

; Read Image File Header (IFH)
byteorder=0
readu,lun,byteorder
tiffversion=0
readu,lun,tiffversion
ifdoffset=0L ;Offset of the first Image File
readu,lun,ifdoffset
IFH={byteorder:byteorder,tiffversion:tiffversion,ifdlength:0,ifdoffset:ifdoffset}
TIFF={IFH:IFH}

; Read Image File Directory (IFD)
ImageStrucIDn=0
code=0u
type=0u
length=0L
while ifdoffset ne 0 do begin
    ; Read Image File header
    point_lun,lun,ifdoffset
    ifdlength=0 ; number of tags
    readu,lun,ifdlength
    IFD={ifdlength:ifdlength}

    ; Read tags + extended tag data
    float=0b
    for i=0l,ifdlength-1 do begin
        tagnamei='TAG'+strtrim(i,2)

        ; Read part of tag
        readu,lun,code
        readu,lun,type
        readu,lun,length

        ; Last 32-bit parameter: Value or pointer
        value=read_tiff_tagvalue(type,nbytes=nbytes)
        if (length)*nbytes gt 4 then begin
            ; Data pointed by pointer
            heap=rebin([value],length*n_elements(value),/sample)

            ; Read pointer
            value=0L
            readu,lun,value

            ; Read data pointed by pointer
            point_lun,-lun,pos
            point_lun,lun,value
            readu,lun,heap
            point_lun,lun,pos

            ; Keep information for reading the image
            if code eq 273 then StripOffsets=heap
            if code eq 279 then StripByteCounts=heap
            if code eq 256 then imagewidth=heap
            if code eq 257 then imagelength=heap
            if code eq 258 then bitspersample=heap
            if code eq 339 then float=value[0] eq 3
            if code eq 259 then LZW=heap eq 5

            ; Add tag
            IFD=create_struct(IFD,tagnamei,{code:code,type:type,length:length,value:value,heap:heap})
        endif else begin
            ; Read value
            if nbytes lt 4 then value=replicate(value,4/nbytes)
            readu,lun,value
    
            ; Keep information for reading the image
            if code eq 273 then StripOffsets=value
            if code eq 279 then StripByteCounts=value
            if code eq 256 then imagewidth=value
            if code eq 257 then imagelength=value
            if code eq 258 then bitspersample=value
            if code eq 339 then float=value[0] eq 3
            if code eq 259 then LZW=value[0] eq 5
    
            ; Add tag
            IFD=create_struct(IFD,tagnamei,{code:code,type:type,length:length,value:value})
        endelse
    endfor
    
    ; Read offset to next IFD
    readu,lun,ifdoffset
    IFD=create_struct(IFD,'ifdoffset',ifdoffset)
    
    ; Read image data
    point_lun,lun,StripOffsets[0]
    image=bytarr(StripByteCounts[0])
    readu,lun,image
    reform_tiff_format,image,bitspersample,bigendian,float=float,LZW=LZW
    
    for i=1,n_elements(StripOffsets)-1 do begin
        point_lun,lun,StripOffsets[i]
        strip=bytarr(StripByteCounts[i])
        readu,lun,strip
        reform_tiff_format,strip,bitspersample,bigendian,float=float,LZW=LZW
        image=[image,temporary(strip)]
    endfor
        
    ; Format image
    if keyword_set(format) then image=reform(image,imagewidth[0],imagelength[0],/overwrite)
    
    IFD=create_struct(IFD,'image',temporary(image))
    
    ; Add IFD to TIFF structure
    TIFF=create_struct(TIFF,'IFD'+strtrim(ImageStrucIDn,2),IFD)
    ImageStrucIDn++
endwhile

size=(fstat(lun)).size
free_lun,lun

error=0b;strucbytesize(TIFF) ne size
return,TIFF

end;function read_tiff_structure
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro sortformat,filter,ext
if ext eq '' then ind=where(filter eq '*',ct,comp=comp) $
else ind=where(strpos(strlowcase(filter), strlowcase(ext)) ne -1,ct,comp=comp)
if ct eq 1 then filter=[filter[ind],filter[comp]]
end;pro sortformat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ChiXtitle,ind,mode=mode,count=count

if n_elements(mode) eq 0 then mode=0
case mode of
1:    ret=['D-spacing','2-theta','Radial Distance','Q-space'] ; Mode
2:    ret=['Angstrom','Degrees','mm','1/nm'] ; Units
3:    ret=[msymbols('angstroms'),msymbols('degrees'),'mm','1/nm'] ; Short Units
4:    ret=['D-spacing (Angstroms)','2-Theta Angle (Degrees)','Radial Distance (mm)','Q (Inverse Nanometres)'] ; Fit2D
else: ret=['D-spacing','2-theta','Radial Distance','Q-space']+' ('+['Angstrom','Degrees','mm','1/nm']+')' ; XTitle
endcase

if size(ind,/type) eq 7 then return,(where(ret eq ind,count))[0] $
else if n_elements(ind) ne 0 then return,ret[ind] else return,ret
        
end;function ChiXtitle
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function specblockextract,lun,block,nlines,out=out

; Find block:
; return: 1 if found, 0 if not found
; out: whole block of nlines. If nlines not given: nlines=1 and ID is cut

out=''
nblock=strlen(block)
line=''
while strmid(line,0,nblock) ne block and not eof(lun) do readf,lun,line

b=strmid(line,0,nblock) eq block
if ~b then print,'No block:"',block,'"' else begin
    if n_elements(nlines) eq 0 then out=strmid(line,nblock) $
    else begin
        out=line
        for i=1,nlines-1 do begin
            readf,lun,line
            out=[out,line]
        endfor
    endelse
endelse
return,b

end;function specblockextract
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function readspec,specfile, top, roi_names, ndimensions,orientation,scannumber=scannumber,lun=lun

ndimensions=0
out=0
orientation=bytarr(4)

keepopen=arg_present(lun)
if n_elements(lun) eq 0 then $
if OPENR_safe(lun, specfile)then begin
    print,'Open error: ',specfile
    return,0
endif

if n_elements(scannumber) eq 0 then scannumber=PromptNumber('1',top,'Scan number:')
scannumber=strtrim(scannumber,2)
if specblockextract(lun,'#S '+scannumber,out=out) then begin

    out=strsplit(out,' ',/extract,count=ct)
    nmca=1

    case out[0] of
    ;ascan  samv -5.54705 -5.44705  10 10
    'ascan': begin
        tran=out[1] eq 'samv'
        if tran then begin
            nx=1
            revx=0b
            ny=fix(out[4])+1
            revy=float(out[2]) gt float(out[3])
        endif else begin
            ny=1
            revy=0b
            nx=fix(out[4])+1
            revx=float(out[2]) gt float(out[3])
        endelse
        ndimensions=1
        endcase
    ;mesh  samh -0.292 -0.282 10  samv 9.738 9.748 10  100
    'mesh': begin
        nx=fix(out[4])+1
        ny=fix(out[8])+1
        revx=float(out[2]) gt float(out[3])
        revy=float(out[6]) gt float(out[7])
        tran=out[1] eq 'samv'
        ndimensions=2
        endcase
    ;zapimage  sampy 0 298 1490 200 rots 93 138 30 0
    'zapimage': begin
        nx=fix(out[4])+1
        ny=fix(out[9])+1
        revx=float(out[2]) gt float(out[3])
        revy=float(out[7]) gt float(out[8])
        tran=0b
        ndimensions=2
        endcase
    ;zapline mono 41.7367 39.0997 660 100
    'zapline':begin
        tran=out[1] eq 'samv'
        if tran then begin
            nx=1
            revx=0b
            ny=fix(out[4]) ; !!!!!
            revy=float(out[2]) gt float(out[3])
        endif else begin
            ny=1
            revy=0b
            nx=fix(out[4]) ; !!!!!
            if out[1] eq 'mono' then revx=0b $
            else revx=float(out[2]) gt float(out[3])
        endelse
        ndimensions=1
        endcase
    endcase
    nmca=nx*ny

    ; MCA format
    if specblockextract(lun,'#L ',out=out) then begin
        roi_names=strsplit(out,' ',/extract,count=nroi_names)
        out=fltarr(nroi_names,nmca)

        line=''
        i=0
        repeat begin
            readf,lun,line
            tmp=strsplit(line,' ',/extract,count=ct)
            if ct eq nroi_names then out[*,i++]=float(tmp)
        endrep until (strmid(line,0,1) eq '#') or eof(lun) or (i eq nmca)


        out=reform(out,nroi_names,nx,ny,/overwrite)
        if revx then out=reverse(out,2)
        if revy then out=reverse(out,3)
        if tran then out=transpose(out,[0,2,1])
        orientation[0:2]=[tran,revx,revy]
    endif

endif

if ~keepopen then free_lun,lun
return,out
end;pro readspec
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TIFFformat,file
tmp=QUERY_TIFF ( file, GEOTIFF=struc)
if tmp and size(struc,/type) eq 8 then if struc.PCSCITATIONGEOKEY eq 'XRDUA_X.1D=2D' then return,3
return,2
end;function TIFFformat    
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RAWformat,file

if OPENR_safe(lun, file) then return,2

line=bytarr(3)
readu,lun,line
if string(line) eq 'RAW' then ret=1 else ret=2

free_lun,lun
return,ret
end;function RAWformat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EDFformat,file

if OPENR_safe(lun, file) then return,2
line='dummy'
while strmid(line,strlen(line)-1,1) ne '}' and not eof(lun) do begin
    readf,lun,line
    if strmid(line,0,4) eq 'xcfg' then begin
        free_lun,lun
        return,3
    endif
endwhile

free_lun,lun
return,2
end;function EDFformat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SPEformat,file
if OPENR_safe(lun, file) then return,2

line=''
readf,lun,line
if line eq '$SPEC_ID:' then ret=1 else ret=2

free_lun,lun
return,ret
end;function SPEformat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ChiFormat,multi=multi,both=both
if keyword_set(both) then return,['*.chi','*.x_y','*.q','*.asc','*.tiff','*.esg']
if keyword_set(multi) then return,['*.tiff','*.esg'] $
else return,['*.chi','*.x_y','*.q','*.asc']
end;function ChiFormat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FileFormat,type,checkfiles

; No extention: .mccd
; Unknown extention: SMART format

; 1D
ext1=['.chi','.x_y','.asc','.d','.q','.plt','.drn','.raw','.uxd','.dat','.prn','.mca','.fio','.spe','.nor']
; 2D
ext2=['.tif','.tiff','.mccd','.edf','.spe','.raw','.dm3','.im','.sif','.jpg','.jpeg','.bmp']
; X.1D = 2D
ext3=['.tif','.tiff','.edf']

n=n_elements(checkfiles)
if n ne 0 then begin
    f=bytarr(n)
    for i=0l,n-1 do begin
        error=CutPath(checkfiles[i],ext=ext)
        ext=strlowcase(ext)
        
        ind1=where(ext1 eq ext,ct1)
        ind2=where(ext2 eq ext,ct2)
        ind3=where(ext3 eq ext,ct3)
        ct=ct1+ct2+ct3
        ind=[ind1,ind2,ind3]
        
        case ct of
        0:    f[i]=2 ; .mccd or .SMART
        1:  f[i]=where(ind ne -1)+1 ; unambiguous format
        else: begin ;ambiguous format
            ext=strlowcase(ext)
            case ext of
            '.raw':    f[i]=RAWformat(checkfiles[i])
            '.edf':    f[i]=EDFformat(checkfiles[i])
            '.spe': f[i]=SPEformat(checkfiles[i])
            '.tif': f[i]=TIFFformat(checkfiles[i])
            '.tiff': f[i]=TIFFformat(checkfiles[i])
            endcase
            endelse
        endcase
    endfor
endif else begin
    case type of
    1: f=ext1 ; 1D (single)
    2: f=ext2 ; 2D
    3: f=[ext1,ext3] ; 1D (single or multiple)
    else:f=''
    endcase
endelse

return,f
end;function FileFormat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SortListFiles,path,sortmethod,separator

case sortmethod of
0:  begin
    n=n_elements(path)
    pos=reform(strpos(path,separator,/REVERSE_SEARCH),1,n)
    ppos=reform(strpos(path,'.',/REVERSE_SEARCH),1,n)
    ind=where(ppos eq -1,ct)
    if ct ne 0 then ppos[ind]=reform(strlen(path[ind]),1,ct)
    ppos-=pos+1

    s=sort(string(FORMAT='(I)',strmid(path,pos+strlen(separator),ppos)))
    out=path[s]
    endcase
1:  out=path[sort(path)]
;2:    begin
;    pos=transpose(strpos(path,separator,/REVERSE_SEARCH))
;    ppos=transpose(strpos(path,'.',/REVERSE_SEARCH))
;    ind=where(ppos eq -1,ct)
;    if ct ne 0 then ppos[ind]=transpose(strlen(path[ind]))
;    ppos-=pos+1
;    
;    formatbase='(A'+string(max(pos),format='(I0)')+')'
;    s=sort(string(strmid(path,0,pos),FORMAT=formatbase)+string(FORMAT='(I)',strmid(path,pos+strlen(separator),ppos)))
;    out=path[s]
;    endcase
2:    begin
    work=file_basename(path)
    nr=make_array(n_elements(path),value='0')
    
    p=stregex(work,'[0-9][0-9]*',length=len)
    ind=where(p ne -1,ct)
    
    while ct ne 0 do begin
        p=reform(p[ind],1,ct,/overwrite)
        len=reform(len[ind],1,ct,/overwrite)
        add=strmid(work[ind],p,len)
        nr[ind]+=string(add,format='(I0'+string(max(strlen(add)),format='(I0)')+')')
        work[ind]=strmid(work[ind],p+len)
        
        p=stregex(work,'[0-9][0-9]*',length=len)
        ind=where(p ne -1,ct)
    endwhile
    
    out=path[sort(nr)]
    
    endcase
endcase

return, out
end;function SortListFiles
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadFit2DSpline,name,error=error

error=1

CATCH, Error_status
IF Error_status NE 0 THEN return,0

if openr_safe(lun,name) then return,0
line = ''
format='(5E14.7)'
readf,lun,line;SPATIAL DISTORTION SPLINE INTERPOLATION COEFFICIENTS
readf,lun,line;
readf,lun,line;  VALID REGION
region=fltarr(4)
readf,lun,region,format=format
readf,lun,line;
readf,lun,line;  GRID SPACING, X-PIXEL SIZE, Y-PIXEL SIZE
grid=fltarr(3)
readf,lun,grid,format=format
readf,lun,line;
readf,lun,line;  X-DISTORTION
n=intarr(2)
readf,lun,n
nx1=n[0]
ny1=n[1]
tx1=fltarr(nx1)
ty1=fltarr(ny1)
c1=fltarr((ny1-4),(nx1-4))
readf,lun,tx1,format=format
readf,lun,ty1,format=format
readf,lun,c1,format=format
readf,lun,line;
readf,lun,line;  Y-DISTORTION
n=intarr(2)
readf,lun,n
nx2=n[0]
ny2=n[1]
tx2=fltarr(nx2)
ty2=fltarr(ny2)
c2=fltarr((ny2-4),(nx2-4))
readf,lun,tx2,format=format
readf,lun,ty2,format=format
readf,lun,c2,format=format

free_lun, lun

GridSpReal=[grid[0],grid[0]]/1000.;um to mm
scrx=grid[1]/1000000.;um to m
scry=grid[2]/1000000.;um to m
GridSpPix=GridSpReal/[scrx,scry]/1000.
GridOffset=[0.,0] ; choose offset to be zero
GridOrigin=[0.,0] ; choose origin to be zero
GridTilt=0. ; choose tilt to be zero
GridDim =(region[2:3]-1-2.*GridOffset)/GridSpPix+1

tck1={uvec:tx1,vvec:ty1,CP:transpose(c1),p:3,q:3,h:nx1-1,k:ny1-1,n:nx1-5,m:ny1-5}; knot locations and coefficients
tck2={uvec:tx2,vvec:ty2,CP:transpose(c2),p:3,q:3,h:nx2-1,k:ny2-1,n:nx2-5,m:ny2-5}

error=0
return,{tck1:tck1,tck2:tck2,GridDim:GridDim,GridOffset:GridOffset,GridOrigin:GridOrigin,GridTilt:GridTilt,$
    GridSpReal:GridSpReal,GridSpPix:GridSpPix,scrx:scrx,scry:scry}
end;function ReadFit2DSpline
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function readPyMCADAT, imsfile, roi_names, ndimensions

ndimensions = 0
rscan_data = -1

if OPENR_safe(lun,imsfile) then return,rscan_data

; Get ROI names
line=''
readf,lun,line
roi_names=strsplit(line,' ',/ext,count=nROI)

; Get number of points in the scan
point_lun,-lun,p
n=0ul
while not eof(lun) do begin
    readf,lun,line
    if strlen(line) ne 0 then n++
endwhile

; Read scan data
rscan_data = fltarr(nROI,n)
point_lun,lun,p
readf,lun,rscan_data

free_lun,lun

; Reform scan data
if roi_names[0] ne 'row' or roi_names[1] ne 'column' then return,-1

ma=max(rscan_data[0,*],min=mi)
nrow=ma-mi+1
ma=max(rscan_data[1,*],min=mi)
ncol=ma-mi+1

rscan_data=rscan_data[2:*,*]
roi_names=roi_names[2:*]
nROI-=2

rscan_data=reform(rscan_data,nROI,ncol,nrow,/overwrite)
ndimensions=2

return,rscan_data
end;function readPyMCADAT
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PRO writeims, wDimensions, wImageData, wIrnames, imsfile

if openw_safe(lun,imsfile) then return

wDimensions=fix(wDimensions)
n=wDimensions+1
datasize = SIZE(wImageData,/dim)
b=n_elements(datasize) ne n
if b then begin
    tmp=lonarr(n)+1
    
    tmp2=datasize
    tmp2[0]=2
    ind=where(tmp2 gt 1,ct)
    ct<=n
    tmp[0:ct-1]=datasize[ind[0:ct-1]]
    
    datasize=temporary(tmp)
    wImageData=reform(wImageData,datasize,/overwrite)
endif

CASE wDimensions OF
  1:     wImageData=transpose(wImageData)
  2:     begin
          wImageData=transpose(wImageData,[1,2,0])
          wImageData=reverse(wImageData,2)
          endcase
  3:     wImageData=transpose(wImageData,[1,2,3,0])
  ELSE:
ENDCASE


PRINTF,lun,wDimensions
datasize = DIMSIZE(wImageData,n)
PRINTF,lun,datasize
PRINTF,lun,wImageData
PRINTF,lun,wIrnames,FORMAT='(A)'

free_lun,lun
END;PRO writeims
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION readims, imsfile, roi_names, ndimensions, noreform=noreform

ndimensions = 0
rscan_data = -1
    
    ; read existing ims file and return
    if OPENR_safe(lun,imsfile) then return,rscan_data
    READF,lun,ndimensions
    roi_names = 0
    d1=0
    d2=0
    d3=0
    d4=0
    CASE ndimensions OF
      1: BEGIN
         READF,lun,d1,d2
         rscan_data = FLTARR(d1,d2)
         roi_names = strarr(d2)
         READF,lun,rscan_data
         READF,lun,roi_names,FORMAT='(A10)'
         rscan_data=reform(rscan_data,d1,d2,/overwrite)
         rscan_data=transpose(rscan_data)
         END
      2: BEGIN
         READF,lun,d1,d2,d3
         rscan_data = FLTARR(d1,d2,d3)
         roi_names = strarr(d3)
         READF,lun,rscan_data
         READF,lun,roi_names,FORMAT='(A10)'
         rscan_data=reform(rscan_data,d1,d2,d3,/overwrite)
         rscan_data=reverse(rscan_data,2) ; reverse Y-axis
         rscan_data=reform(rscan_data,d1,d2,d3,/overwrite)
         rscan_data=transpose(rscan_data,[2,0,1])
         END
      3: BEGIN
         READF,lun,d1,d2,d3,d4
         rscan_data = FLTARR(d1,d2,d3,d4)
         roi_names = strarr(d4)
         READF,lun,rscan_data
         READF,lun,roi_names,FORMAT='(A10)'
         rscan_data=reform(rscan_data,d1,d2,d3,d4,/overwrite)
         rscan_data=transpose(rscan_data,[3,0,1,2])
         END
      ELSE:
      ENDCASE
    free_lun,lun
    ndimensions = LONG(ndimensions)

RETURN,rscan_data
END;FUNCTION readims
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function readMICROXASFASTXRF, imsfile, roi_names, ndimensions, scandim

; Open XML file
oDoc = OBJ_NEW( 'IDLffXMLDOMDocument', FILENAME=imsfile,/quiet )
if not obj_valid(oDoc) then begin
    ndimensions = 0
    return,-1
endif

; Structure:
;...
;<scan>
;    ...
;    <dimension zigZag="true" ...> ; CHILD OF "scan"
;        <positioner id="posX" ...> ; CHILD OF "dimension"
;            <start>0.0</start> ; CHILD OF "positioner"
;            <end>360.0</end> ; SIBLING OF "start"
;            <stepSize>2.0</stepSize>
;            <integrationTime>0.2</integrationTime>
;        </positioner>
;        <detector scaler="0" id="ROI0" arraySize="-1"/>
;        <detector scaler="1" id="ROI1" arraySize="-1"/>
;        ...
;    </dimension>
;    <dimension zigZag="false" ...>
;        <positioner settlingTime="0.0" readback="X05LA-ES2-MA1:TRX.RBV" id="posY" channel="X05LA-ES2-MA1:TRX">
;            <start>-1.3</start>
;            <end>-0.7</end>
;            <stepSize>0.01</stepSize>
;        </positioner>
;    </dimension>
;</scan>
;
; The first dimension is the inner loop
; The second dimension is the outer loop

; Node walker
oPointer=oDoc->CreateTreeWalker(OBJ_NEW())
oNode = oPointer->FirstChild()
WHILE OBJ_VALID( oNode ) DO BEGIN
    
    
;    print,oNode->GetNodeName()
;    
;    type=oNode->GetNodeType()
;    case type of
;    3:    begin ; IDLffXMLDOMText
;        print,oNode->GetNodeValue()
;        endcase
;    1:    begin ; IDLffXMLDOMElement
;        oAttributes=oNode->GetAttributes()
;        nattibutes=oAttributes->GetLength()
;        for i=0l,nattibutes-1 do begin
;            oAttibute=oAttributes->Item(i) ;IDLffXMLDOMAttr
;            
;            print,oAttibute->GetName(),'=',oAttibute->GetValue()
;        endfor
;        endcase
;    else:
;    endcase
;    
;    oNode=oPointer->NextNode()
;    continue
    

    if strlowcase(oNode->GetNodeName()) eq 'scan' then begin
        ; <scan>
        oNode = oPointer->FirstChild()
        WHILE OBJ_VALID( oNode ) DO BEGIN
        
            ; <dimension>
            if strlowcase(oNode->GetNodeName()) eq 'dimension' then begin
                ;oAttributes=oNode->GetAttributes()
                ;oAttribute=oAttributes->GetNamedItem('zigZag')
            
                oNode=oPointer->FirstChild()
                WHILE OBJ_VALID( oNode ) DO BEGIN
                    case strlowcase(oNode->GetNodeName()) of
                    'detector':    begin
                                ; Store detector ID
                                oAttributes=oNode->GetAttributes()
                                oAttribute=oAttributes->GetNamedItem('id')
                                if OBJ_VALID( oAttribute ) then $
                                if n_elements(roi_names) eq 0 then roi_names=oAttribute->GetValue() else roi_names=[roi_names,oAttribute->GetValue()]
                                endcase
                    'positioner':begin
                                ; Store motor ID and properties
                                oAttributes=oNode->GetAttributes()
                                oAttribute=oAttributes->GetNamedItem('id')
                                if OBJ_VALID( oAttribute ) then $
                                if n_elements(motor_names) eq 0 then motor_names=oAttribute->GetValue() else motor_names=[motor_names,oAttribute->GetValue()]
                
                                if n_elements(scandim) eq 0 then begin
                                    scandim=fltarr(4)
                                    iscandim=0
                                endif else begin
                                    scandim=[[scandim],[fltarr(4)]]
                                    iscandim++
                                endelse
                                
                                oNode=oPointer->FirstChild()
                                WHILE OBJ_VALID( oNode ) DO BEGIN
                                    case strlowcase(oNode->GetNodeName()) of
                                    'start':     begin
                                                oNode=oPointer->FirstChild()
                                                scandim[0,iscandim]=float(oNode->GetNodeValue())
                                                oNode=oPointer->ParentNode() 
                                                endcase
                                    'end':         begin
                                                oNode=oPointer->FirstChild()
                                                scandim[1,iscandim]=float(oNode->GetNodeValue())
                                                oNode=oPointer->ParentNode() 
                                                endcase
                                    'stepsize': begin
                                                oNode=oPointer->FirstChild()
                                                scandim[2,iscandim]=float(oNode->GetNodeValue())
                                                oNode=oPointer->ParentNode() 
                                                endcase
                                    'integrationtime': begin
                                                oNode=oPointer->FirstChild()
                                                scandim[3,iscandim]=float(oNode->GetNodeValue())
                                                oNode=oPointer->ParentNode() 
                                                endcase
                                    else:
                                    endcase
                                    
                                    oNode=oPointer->NextSibling()
                                ENDWHILE
                                oNode=oPointer->ParentNode() 
                                endcase
                    else:
                    endcase
                    
                    oNode=oPointer->NextSibling()
                ENDWHILE
                
                oNode=oPointer->ParentNode() 
                
            endif
            ; </dimension>
            
            oNode=oPointer->NextSibling()
        ENDWHILE
        ; </scan>

        break
    endif else oNode=oPointer->NextSibling()
ENDWHILE

OBJ_DESTROY, oDoc

if n_elements(roi_names) eq 0 or n_elements(motor_names) eq 0 then begin
    ndimensions = 0
    return,-1
endif

; Find corresponding text file
datafile=file_search(file_dirname(imsfile),file_basename(imsfile,'.xml')+'*.txt',count=n)
if n eq 0 then begin
    ndimensions = 0
    return,-1
endif
if OPENR_safe(lun,datafile[0]) then begin
    ndimensions = 0
    return,-1
endif

line='#'
bstart=1b
p=0l
while strmid(line,0,1) eq '#' do begin
    point_lun,-lun,p
    readf,lun,line
    line=strlowcase(line)
    if strpos(line,'dimension') ne -1 then begin ; motor identifier (1=inner motor, 2=outer motor, ...)
        tmp=strsplit(line,':',/extract,count=ct)
        if ct ne 2 then continue
        if bstart then motorcol=long(tmp[1])-1 else motorcol=[motorcol,long(tmp[1])-1]
        
        point_lun,-lun,p
        readf,lun,line
        line=strlowcase(line)
        if strpos(line,'size') ne -1 then begin ; number of datapoints for this motor
            tmp=strsplit(line,':',/extract,count=ct)
            if ct ne 2 then continue
            if bstart then begin
                dim=long(tmp[1])
                bstart=0b
            endif else dim=[dim,long(tmp[1])]
        endif
    endif
endwhile
point_lun,lun,p

; motor_names,scandim = inner loop, outer loop
; motorcol = [1,0] => outer loop in the first column, inner loop in the second column
roi_names=[motor_names[motorcol],roi_names]
dim=dim[motorcol]
ndimensions=n_elements(dim)

; Read data
nROI=n_elements(roi_names)
data=fltarr(nROI,product(dim,/int))
readf,lun,data
data=reform(data,[nROI,dim],/overwrite)

; Adapt scandim to real motor positions
for i=0,ndimensions-1 do begin

    ; Get motor positions from encoders
    bskip=0b
    case motorcol[i] of
    0:    begin
        startstop=reform(data[i,[0,dim[motorcol[i]]-1],0])
        endcase
    1:    begin
        startstop=reform(data[i,0,[0,dim[motorcol[i]]-1]])
        endcase
    else: bskip=1b
    endcase
    if bskip then continue
    
    ; Replace user supplied values by encoder values?
    diff=(scandim[1,motorcol[i]]>scandim[0,motorcol[i]])-(scandim[1,motorcol[i]]<scandim[0,motorcol[i]])
    dimorg=round(diff/abs(scandim[2,motorcol[i]])+1)
    if dim[motorcol[i]] ne dimorg then begin
        scandim[0:1,motorcol[i]]=startstop
        diff=(scandim[1,motorcol[i]]>scandim[0,motorcol[i]])-(scandim[1,motorcol[i]]<scandim[0,motorcol[i]])
        scandim[2,motorcol[i]]=diff/(dim[motorcol[i]]-1.)
    endif
    
endfor

free_lun,lun
RETURN,data
end;function readMICROXASFASTXRF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WriteAthenaDat,file,header,data
openw,lun,file,/get_lun,error=err
if err ne 0 then return,0
n=strtrim(n_elements(header),2)
printf,lun,header,format='('+n+'(A25))'
printf,lun,double(data),format='('+n+'E)'
free_lun,lun
return,1b
end;function WriteAthenaDat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro write_esg,file,d,ROIAzimuth,nsectors,dist=dist,error=error

error=1b
if openw_safe(lun,file) then return

if n_elements(dist) eq 0 then dist='?' else dist=stringr(dist*100)

b=(ROIAzimuth[1]<ROIAzimuth[0])*180/!pi
e=(ROIAzimuth[1]>ROIAzimuth[0])*180/!pi
inc=(e-b)/nsectors
ind=b+inc*lindgen(nsectors)
for i=0l,nsectors-1 do begin
    printf,lun,'_pd_block_id noTitle|#',stringr(i)
    printf,lun,''
    if i eq 0 then begin
        printf,lun,'_diffrn_detector Image Plate'
        printf,lun,'_diffrn_detector_type ?'
        printf,lun,'_pd_meas_step_count_time ?'
        printf,lun,'_diffrn_measurement_method ?'
        printf,lun,'_diffrn_measurement_distance_unit cm'
        printf,lun,'_pd_instr_dist_spec/detc ',dist
        printf,lun,'_diffrn_radiation_wavelength ?'
        printf,lun,'_diffrn_source_target ?'
        printf,lun,'_diffrn_source_power ?'
        printf,lun,'_diffrn_source_current ?'
        printf,lun,'_pd_meas_angle_omega 0.0'
        printf,lun,'_pd_meas_angle_chi 0.0'
        printf,lun,'_pd_meas_angle_phi 0.0'
        printf,lun,'_riet_par_spec_displac_x 0'
        printf,lun,'_riet_par_spec_displac_y 0'
        printf,lun,'_riet_par_spec_displac_z 0'
        printf,lun,'_riet_meas_datafile_calibrated false'
    endif
    printf,lun,'_pd_meas_angle_eta ',ind[i]
    printf,lun,''
    printf,lun,'loop_'
    printf,lun,'_pd_proc_2theta_corrected'
    printf,lun,'_pd_calc_intensity_total'
    printf,lun,d[[0,i+1],*],format='(f0.14," ",f0.05)'
    printf,lun,''
endfor

free_lun, lun
error=0b

end;pro write_esg
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro WriteChi_,path,d,xtype,tiffpath,ROIAzimuth,origformat=origformat,error=error,bring=bring,berror=berror

error=1b
if openw_safe(lun,path) then return

str=ChiXtitle(xtype,mode=4)
strind=strpos(str,' ')
s=size(d)
printf,lun,tiffpath+': '+strmid(str,0,strind)+' Scan'
printf,lun,str
printf,lun,' '
printf,lun,s[2]
printf,lun,format='(2E16.7)',d[0:1,*]
if not keyword_set(origformat) then begin
    printf,lun,'$'
    printf,lun,ROIAzimuth
    off=2
    if keyword_set(berror) then begin
        printf,lun,'$YError:'
        printf,lun,format='(E16.7)',d[off,*]
        off++
    endif
    if keyword_set(bring) then begin
        printf,lun,'$Ring%:'
        printf,lun,format='(E16.7)',d[off,*]
    endif
endif
free_lun,lun
error=0b

end;pro WriteChi_
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WriteChi,path_,d,xtype,tiffpath,OutID,ROIAzimuth,origformat=origformat,nsectors=nsectors,bring=bring,berror=berror,_extra=_extra

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
;message,!ERR_STRING,/info
return,0b
ENDIF

; Path
path=path_
if path eq '' then return,0B
temp=cutpath(path,path=p,file=f,ext=ext)
sext=strlowcase(ext)

; Sectors
if n_elements(nsectors) eq 0 then nsectors=1
nsectors>=1
if nsectors eq 1 then begin
    sec=''
    secind=ROIAzimuth
endif else begin
    if n_elements(ROIAzimuth) eq 0 then begin
        b=0.
        e=360.
    endif else begin
        b=(ROIAzimuth[1]<ROIAzimuth[0])*180/!pi
        e=(ROIAzimuth[1]>ROIAzimuth[0])*180/!pi
    endelse
    inc=(e-b)/nsectors
    secind=b+inc*lindgen(nsectors+1)
    sec=string(secind,format='("_",I04)')
    sec=string(lindgen(nsectors+1),format='("_azi",I0'+stringr(strlen(stringr(nsectors)))+')')+sec
    secind*=!pi/180
endelse

; d: X, [Y] x nsectors, [Yerror] x nsectors, Ring% x nsectors
; => normally in column format, except for tiff
bring=keyword_set(bring)
berror=keyword_set(berror)
typed=berror+2*bring

case sext of
'.chi'    :begin
        for i=0,nsectors-1 do begin
            path=p+f+sec[i]+ext
            
            case typed of
            0: tmp=d[[0,i+1],*]
            1: tmp=d[[0,i+1,i+1+nsectors],*]
            2: tmp=d[[0,i+1,i+1+nsectors],*]
            3: tmp=d[[0,i+1,i+1+nsectors,i+1+2*nsectors],*]
            endcase

            WriteChi_,path,tmp,xtype,tiffpath,secind[i:i+1],origformat=origformat,error=error,bring=bring,berror=berror
            if error then return, 0b
        endfor
        endcase
'.x_y'    :begin
        for i=0,nsectors-1 do begin
            path=p+f+sec[i]+ext
            if openw_safe(lun,path) then return,0
            printf,lun,format='(2E16.7)',d[[0,i+1],*]
            free_lun,lun
        endfor
        endcase
'.q'    :begin
        for i=0,nsectors-1 do begin
            path=p+f+sec[i]+ext
            if openw_safe(lun,path) then return,0
            printf,lun,format='(2E16.7)',d[[0,i+1],*]
            free_lun,lun
        endfor
        endcase
'.asc'    :begin
        for i=0,nsectors-1 do begin
            path=p+f+sec[i]+ext
            if openw_safe(lun,path) then return,0
            printf,lun,format='(2E16.7)',d[[0,i+1],*]
            free_lun,lun
        endfor
        endcase
'.mca':    begin
        for i=0,nsectors-1 do begin
            path=p+f+sec[i]+ext
            write_mca,path,reform(d[i+1,*]),path,error=error
            if error then return, 0b
        endfor
        endcase
'.spe':    begin
        for i=0,nsectors-1 do begin
            path=p+f+sec[i]+ext
            write_spe,path,reform(d[i+1,*]),path,error=error
            if error then return, 0b
        endfor
        endcase
'.fio':    begin
        for i=0,nsectors-1 do begin
            path=p+f+sec[i]+ext
            write_fio,path,d[i+1,*],path,error=error
            if error then return, 0b
        endfor
        endcase
'.esg':    begin
        write_esg,path,d,ROIAzimuth,nsectors,error=error,_extra=_extra
        if error then return, 0b
        endcase
else:    begin
        if sext eq '.tif' or sext eq '.tiff' then begin
            str=ChiXtitle(xtype,mode=4)
            strind=strpos(str,' ')
            GEOTIFF={GTCITATIONGEOKEY:tiffpath+': '+strmid(str,0,strind)+' Scan',$
                    GEOGCITATIONGEOKEY:str,$
                    GTMODELTYPEGEOKEY:berror,$
                    GTRASTERTYPEGEOKEY:bring,$
                    GEOGLINEARUNITSIZEGEOKEY:ROIAzimuth[0],$
                    GEOGANGULARUNITSIZEGEOKEY:ROIAzimuth[1],$
                    PCSCITATIONGEOKEY:'XRDUA_X.1D=2D'}

            ; d should be in row format
            write_tiff_safe,path,d,GEOTIFF=GEOTIFF,error=error
            if error then return,0b
        endif else begin
            printw,OutID,path+' unknown extention'
            return, 0B
        endelse
        endelse
endcase

printw,OutID,path+' saved'
return, 1B
end;function WriteChi
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadLOF,path,lof,prid

CATCH, Error_status
IF Error_status NE 0 THEN return,0

if OPENR_safe(lun, path) then begin
    printw,prid,path+' NOT loaded'
    return,0B
endif

line=''

block='$# columns:'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
if not eof(lun) then begin
    ct=0L
    readf,lun,ct
    printw,prid,'Number of columns: '+stringr(ct)
endif

point_lun,lun,0
block='$# rows:'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
if not eof(lun) then begin
    ct=0L
    readf,lun,ct
    printw,prid,'Number of rows: '+stringr(ct)
endif

point_lun,lun,0
block='$Files:'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
if not eof(lun) then begin
    ct=0L
    readf,lun,ct
    lof=strarr(ct)
    for i=0l,ct-1 do begin
       readf,lun,line
       lof[i]=line
    endfor
endif

FREE_LUN, lun

printw,prid,path+' loaded'
return,1B
end; function ReadLOF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReformXDIDatablock,block,TRAN,REV1,REV2,REV3,SWAP,s0,s1,s2
; block: datablock
; TRAN,REV1,...: TRAN[0] -> TRAN[1]; REV[0] -> REV[1]; ...
; s0,s1,s2: optional input, output dimensions

; Block size
s=DimSize(block,3)
dim0=s[0]
dim1=s[1]
dim2=s[2]

if n_elements(s0) ne 0 and n_elements(s1) ne 0 and n_elements(s2) ne 0 then begin
    if dim1*dim2 eq 1 and s1*s2 ne 1 then begin
        ; Filenames array
        dim0=1
        dim1=s1
        dim2=s2
    endif else begin 
        dim0=s0
        dim1=s1
        dim2=s2
    endelse
endif

block=reform(block,dim0,dim1,dim2,/overwrite)

; Convert to default orientation
if REV2[0] and dim2 gt 1 then begin
    block=reverse(block,3)
    block=reform(block,dim0,dim1,dim2,/overwrite)
endif
if REV1[0] and dim1 gt 1 then begin
    block=reverse(block,2)
    block=reform(block,dim0,dim1,dim2,/overwrite)
endif
if TRAN[0] then begin
    block=transpose(block,[0,2,1])
    t=dim1
    dim1=dim2
    dim2=t
    block=reform(block,dim0,dim1,dim2,/overwrite)
endif
if REV3[0] and dim2 gt 1 then begin
    indodd=2*lindgen(dim2/2)+1
    block[*,*,indodd]=reverse(block[*,*,indodd],2)
endif

; SWAP dim1 and dim2
if SWAP then begin
    t=dim1
    dim1=dim2
    dim2=t
    block=reform(block,dim0,dim1,dim2,/overwrite)
endif

; Convert to new orientation
if REV3[1] and dim2 gt 1 then begin
    indodd=2*lindgen(dim2/2)+1
    block[*,*,indodd]=reverse(block[*,*,indodd],2)
endif
if TRAN[1] then begin
    block=transpose(block,[0,2,1])
    t=dim1
    dim1=dim2
    dim2=t
    block=reform(block,dim0,dim1,dim2,/overwrite)
endif
if REV1[1] and dim1 gt 1 then begin
    block=reverse(block,2)
    block=reform(block,dim0,dim1,dim2,/overwrite)
endif
if REV2[1] and dim2 gt 1 then begin
    block=reverse(block,3)
    block=reform(block,dim0,dim1,dim2,/overwrite)
endif

s0=dim0
s1=dim1
s2=dim2

end;pro ReformXDIDatablock
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadXDI,path,file,OutID,error=error,format=format

nError=0
error=0B

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,OutID,path+file+' not found.'
    error=1B
    return,0
ENDIF

dim=0
if openr_safe(lun, path+file) then begin
    printw,OutID,path+file+' not found.'
    error=1B
    return,0
endif
line=''

; ----Read Header block----
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,OutID,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
    block='$HEADING:'
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
        readf,lun,line
        header=line
        while strmid(line,0,1) ne '$' do begin
            readf,lun,line
            header=[header,line]
        endwhile
        line=''
    endif
ENDELSE

; ----Read PhysMask block----
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,OutID,block+' corrupted or older version'
    ; Don't issue error
    ;nError=nError+1
    ;error=1B
ENDIF ELSE BEGIN
    point_lun,lun,0
    block='$PHYSMASK:'
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
        nr=0
        readf,lun,nr

        readf,lun,line ;"Qualitative analysis:"

        narr=lonarr(4)
        readf,lun,narr
        narr=total(narr)

        readf,lun,line
        pmfields=strsplit(line,' ',/extract)
        pmfields=pmfields[1:*]
        npmfields=n_elements(pmfields) ; all fields except group number

        pmstruc=create_struct('t',-1)
        pmtags=-1
        val=-1

        ; read qualitative
        if narr ne 0 then begin
            pmgroup=-1
            pmdata=strarr(npmfields)
            add=0

            block2='Quantitative'
            len=strlen(block2)
            readf,lun,line
            while (strmid(line,0,len) ne block2) do begin
                tmp=strsplit(line,' ',/extract,count=ntmp)
                
                pmdata=[[pmdata],[tmp[1:npmfields]]] ; add group info

                tmp=fix(tmp[0]) ; group number
                if tmp lt val then add+=val+1
                val=tmp

                pmgroup=[pmgroup,val+add] ; add group number

                readf,lun,line
            endwhile
            val+=add

            for i=1,n_elements(pmgroup)-1 do begin ; loop over all groups
                ind=where(pmtags eq pmgroup[i],ct)
                if ct eq 0 then begin ; group number not in pmtags
                    ind=where(pmgroup eq pmgroup[i],ct) ; find all masks with the same group number
                    struc={n:ct,data:[pmdata[*,ind]],labeli:0b,combi:0b} ; add group to pmstruc
                    pmstruc=create_struct(pmstruc,'t'+stringr(pmgroup[i]),struc)
                    pmtags=[pmtags,pmgroup[i]] ; add group to pmtags
                endif
            endfor

        endif else readf,lun,line ;"Quantitative analysis:"

        narr=lonarr(3)
        readf,lun,narr
        narr=0
        readf,lun,narr

        readf,lun,line
        pmfields2=strsplit(line,' ',/extract)
        npmfields2=n_elements(pmfields2)

        ; read quantitative
        for i=0l,narr-1 do begin
            readf,lun,line
            tmp=strsplit(line,' ',/extract,count=ntmp)
            if ntmp gt npmfields2 then begin
                tmp=[strjoin(tmp[0:ntmp-npmfields2],'_'),tmp[ntmp-npmfields2+1:*]]
            endif
            val++ ; group number
            struc={n:1,data:tmp,labeli:1,combi:0b} ; add group to pmstruc
            pmstruc=create_struct(pmstruc,'t'+stringr(val),struc)
            pmtags=[pmtags,val] ; add group to pmtags
        endfor

        if npmfields gt npmfields2 then pmfields2=[pmfields2,replicate('-',npmfields-npmfields2)]
        if npmfields2 gt npmfields then pmfields=[pmfields,replicate('-',npmfields2-npmfields)]
        pmfields=[[pmfields],[pmfields2]]
        npmfields>=npmfields2

        ; additional?
        readf,lun,line
        if line eq 'Combined:' then begin
            pmgroup=-1
            pmdata=strarr(npmfields)
            combi=0b
            add=0
            val++ ; 

            block2='$'
            len=strlen(block2)
            readf,lun,line
            while (strmid(line,0,len) ne block2) do begin
                tmp=strsplit(line,' ',/extract,count=ntmp)
                pmdata=[[pmdata],[tmp[1:npmfields]]] ; add group info

                p=strpos(tmp[0],'*')
                if p ne -1 then begin
                    combi=[combi,2]
                    tmp=fix(strmid(tmp[0],0,p)) ; group number
                endif else begin
                    combi=[combi,1]
                    tmp=fix(tmp[0]) ; group number
                endelse
                if tmp lt val then add+=val+1
                val=tmp

                pmgroup=[pmgroup,val+add] ; add group number

                readf,lun,line
            endwhile

            for i=1,n_elements(pmgroup)-1 do begin ; loop over all groups
                ind=where(pmtags eq pmgroup[i],ct)
                if ct eq 0 then begin ; group number not in pmtags
                    ind=where(pmgroup eq pmgroup[i],ct) ; find all masks with the same group number
                    struc={n:ct,data:[pmdata[*,ind]],labeli:0b,combi:combi[i]} ; add group to pmstruc
                    pmstruc=create_struct(pmstruc,'t'+stringr(pmgroup[i]),struc)
                    pmtags=[pmtags,pmgroup[i]] ; add group to pmtags
                endif
            endfor
        endif
        
        ; Reorder?
        if line eq '$Reorder:' then begin
            in=lonarr(n_elements(pmtags))
            readf,lun,in
            in=sort(in)
            pmstruc=IndexStruct(pmstruc,in)
            pmtags=pmtags[in]
        endif

    endif
ENDELSE

; ----Read IMS block----
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,OutID,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
    point_lun,lun,0
    block='$IMS:'
    len=strlen(block)
    line=''
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line

    ; Datablock size (3 dimenions)
    readf,lun,dim
    if dim ne 2 and dim ne 3 then begin
        printw,OutID,'XDI: wrong dimensions!'
        nError=nError+1
        error=1B
    endif
    s=ulonarr(dim)
    readf,lun,s
    if dim eq 2 then s=[s,1]
    
    ; Read 3D datablock
    ; s[0]: number of groups
    ; s[1]: nCol
    ; s[2]: nRow
    data=fltarr(s[0],s[1],s[2])
    data=reform(data,s,/overwrite)
    if keyword_set(format) then begin
        
        ; This is for old XDI formats
        if strpos(strlowcase(format),'unformatted') ne -1 then readu,lun,data else begin
            tmp=strsplit(format,'FDE()',/extract)
            if tmp[0] eq '?' then begin
                format_in='('+stringr(s[0])+strmid(format,strlen(tmp[0]))+')'
                readf,lun,data,format=format_in
            endif else begin
                if n_elements(tmp) ne 2 or s[0] lt tmp[0] then createerror[0]=0
                format_in='('+format+')'
                row=fltarr(s[0])
                datai=0l
                line=''
                for j=1,s[2] do begin
                    for i=1,s[1] do begin
                        readf,lun,row,format=format_in
                        data[datai:datai+s[0]-1]=row
                        datai+=s[0]
                    endfor
                    ; New group
                    if j ne s[2] and tmp[0] ne s[0] then readf,lun,line
                endfor
            endelse
        endelse
        
        ; 2D block: ngroups,nfiles
        ; 3D block: ncol,nrow,ngroups
        if dim eq 3 then begin
            ind=[2,0,1]
            data=transpose(data,ind)
            s=s[ind]
            data=reform(data,s,/overwrite)
        endif
    endif else begin
        readu,lun,data
    endelse
ENDELSE

; ----Read Files block----
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,OutID,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
    point_lun,lun,0
    block='$FILES:'
    len=strlen(block)
    line=''
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
        ; ----Make array to save filenames----
        n=s[1]*s[2]
        filenames=strarr(n)

        ; ----Read all strings and put them in the right place----
        readf,lun,filenames,format='(A)'
        
        ; Get first and check for exist
        tmp = file_search(filenames[0],count=ct)
        bpath=0b
        if ct eq 0 and widget_info(OutID,/valid) then begin
            t=CutPath(filenames[0],file=filef,ext=ext)
            pathf=path
            tmp = file_search(pathf+filef+ext,count=ct)
            bpath=ct ne 0
            if ~bpath then begin
                pathf=DIALOG_PICKFILE(filter='*'+ext,title='Select path of datafiles...',/DIRECTORY)
                bpath=pathf ne ''
                if bpath then begin
                    tmp = file_search(pathf+filef+ext,count=ct)
                    bpath=ct ne 0
                    if bpath then filenames[0]=tmp[0]
                endif
            endif else filenames[0]=tmp[0]
        endif

        if bpath then $
        for k=1l,n-1 do begin ;Loop over all points in the map
               t=CutPath(filenames[k],file=filef,ext=ext)
               filenames[k]=pathf+filef+ext
        endfor
    endif else filenames=-1
ENDELSE

; ----Read orientation----
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,OutID,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
    point_lun,lun,0
    block='$ORIENTATION:'
    orientation=bytarr(4)
    sign=[1,1]
    scandim=0
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
        readf,lun,orientation
        readf,lun,sign
        readf,lun,scandim
    endif else begin
        ; This is for old xdi files
        header_scan=strsplit(header[0],' ',/extract)
        orientation[0]=header_scan[1] eq 'samv'
        case dim of
        3:    begin
            scandim=2
            if header_scan[0] eq 'meshc' then orientation[3]=1b
            if orientation[0] then begin
                orientation[1]=float(header_scan[6]) gt float(header_scan[7])
                orientation[2]=float(header_scan[2]) gt float(header_scan[3])
            endif else begin
                orientation[1]=float(header_scan[2]) gt float(header_scan[3])
                orientation[2]=float(header_scan[6]) gt float(header_scan[7])
            endelse
            endcase
        2:    begin
            scandim=1
            if orientation[0] then $
                orientation[2]=float(header_scan[2]) gt float(header_scan[3]) $
            else $
                orientation[1]=float(header_scan[2]) gt float(header_scan[3])
            endcase
        endcase
    endelse
    
    ; Convert default to current orientation:
    if n_elements(filenames) ne 0 then $
    if (size(filenames[0],/type)) eq 7 then begin
        ncol=s[1]
        nrow=s[2]
        if orientation[0] then begin
            ncol=s[2]
            nrow=s[1]
        endif
        ReformXDIDatablock,filenames,[0b,orientation[0]],[0b,orientation[1]],[0b,orientation[2]],[0b,orientation[3]],0b,1,ncol,nrow
    endif
ENDELSE

free_lun,lun
printw,OutID,'XDI loaded from: '+path+file+' with '+stringr(nError)+' errors'

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,OutID,!error_state.msg
    nError=nError+1
    error=1B
    return,0
ENDIF ELSE BEGIN
    return,{data:temporary(data),filenames:temporary(filenames),header:header,$
    pmstruc:pmstruc,pmtags:pmtags,pmfields:pmfields,npmfields:npmfields,scandim:scandim,$
    orientation:orientation,sign:sign}
ENDELSE

end;function ReadXDI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Select_Files,top,list,$
    sortmethod=sortmethod,separator=separator,filter=filter,$
    outpath=outpath,outfile=outfile,description=description

if keyword_set(description) then title='Select '+description+' manually?' $
else title='Select files manually?'
Result = DIALOG_MESSAGE(title,/question)

path=''
if result eq 'Yes' then begin
    if not keyword_set(filter) then filter=['*'+list.Format,'*.lof','*.xdi']
    path=CAPDIALOG_PICKFILE(path=list.path,file=list.file,filter=filter,$
       /MULTIPLE_FILES,title='Select files...')
    if path[0] EQ '' then return,''
    if (n_elements(path) eq 1) then begin
              error=CutPath(path[0],path=path2,file=outfile,ext=ext)
              ext=strlowcase(ext)
              case ext of
               '.lof':    begin
                    error=ReadLOF(path,lof,list.OutID)
                    if error eq 0 then return,'' else path=lof
                    endcase
            '.xdi':    begin
                    filenames=ReadXDI(path2,outfile+ext,list.OutID,error=error)
                    if error then return,''
                    path=filenames.filenames
                    endcase
            else:
            endcase
    endif
    path=path[sort(path)]
    error=CutPath(path[0],path=outpath,file=outfile)
endif else begin
    outpath=DIALOG_PICKFILE(path=list.path,title='Select Directory...',/DIRECTORY)
    if outpath[0] EQ '' then return,''
    outpath=outpath[0]
    ;----Prompt for filter----
    nFormat=n_elements(list.Format)
    if nFormat eq 1 then begin
        filter=PromptNumber(outpath+'*'+list.Format,top,'Enter filter string:')
        path=file_search(filter)
    endif else begin
        filter=outpath+';*'+list.Format[0]
        for i=1,nFormat-1 do filter+=';*'+list.Format[i]
        filter=PromptNumber(filter,top,'Enter filter string:')
        tmp=strsplit(filter,';',/extract,count=ct)
        if ct lt 2 then return,''
        outpath=tmp[0]
        filter=tmp[1:*]
        path=file_search(outpath+filter)
    endelse
    path=path[sort(path)]
    error=CutPath(path[0],path=outpath,file=outfile)
endelse

if keyword_set(separator) then path=SortListFiles(path,sortmethod,separator)
return,path
end;function Select_Files
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadOverlapFile,path,file,result=result

CATCH, Error_status
IF Error_status NE 0 THEN return,0b

pathin=''
pathin=DIALOG_PICKFILE(path=path,file=file,$
         filter=['*.txt','*.*'],title='Read overlap table...')

if pathin eq '' then return,0B
tmp=CutPath(pathin,path=path,file=file,ext=ext)
file=file+ext

if openr_safe(lun,pathin) then return,0b

block='$nProf:'
len=strlen(block)
line=''
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
nProf=0L
if not eof(lun) then begin
    readf,lun,nProf
endif

if nProf eq 0 then return,0B
Profiles=PTRARR(nProf)
TableXbY=fltarr(nProf,nProf)
TableXY=fltarr(nProf,nProf)

point_lun,lun,0
block='$XrowByYcol:'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
if not eof(lun) then readf,lun,TableXbY,format='('+stringr(nProf)+'f15.7)'

point_lun,lun,0
block='$XY:'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
if not eof(lun) then readf,lun,TableXY,format='('+stringr(nProf)+'f15.7)'

point_lun,lun,0
block='$Profiles:'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
if not eof(lun) then begin
    for i=0l,nProf-1 do begin
        n=0L
        readf,lun,n
        d=fltarr(n,2)
        readf,lun,d,format='('+stringr(n)+'f15.7)'
        Profiles[i]=PTR_NEW(d)
    endfor
endif

point_lun,lun,0
block='$Grouping:'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
if not eof(lun) then begin
    overlap=0.
    readf,lun,overlap
    IndMask=lonarr(nProf)
    readf,lun,IndMask
    maskgroup=lonarr(nProf)
    readf,lun,maskgroup
endif

point_lun,lun,0
block='(groupnr, d, SD):'
while (strpos(line,block) eq -1) and (not eof(lun)) do readf,lun,line
if not eof(lun) then begin
    nGroups=0L
    readf,lun,nGroups
    TableMean=fltarr(3,nGroups)
    tmp=0L
    readf,lun,tmp
    readf,lun,TableMean,format='(I'+stringr(tmp)+',2f15.7)'
endif

result={XY:TableXY,XbY:TableXbY,Profiles:Profiles,$
nGroups:nGroups,maskgroup:maskgroup,TableMean:TableMean,nProf:nProf,IndMask:IndMask}
free_lun,lun

return,1B
end;function ReadOverlapFile
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SelFile,path,file,filter,title
;repeat begin
;    patht=DIALOG_PICKFILE(path=path,file=file,$
;         filter=filter,title=title)
;    if patht EQ '' then return,patht
;    Result = file_search(patht,count=ct)
;    if ct ne 0 then begin
;       Result = DIALOG_MESSAGE('Overwrite file?',/question)
;       if result eq 'Yes' then ct=0
;    endif
;endrep until (ct eq 0)

patht=''
patht=CAPDIALOG_PICKFILE(path=path,file=file,$
         filter=filter,title=title,/OVERWRITE_PROMPT,/WRITE)
return,patht
end;function SelFile
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION PSWINDOW_ASPECT, aspectRatio, MARGIN=margin, WindowAspect=wAspectRatio

ON_ERROR, 1

   ; Check for aspect ratio parameter and possibilities.

IF N_PARAMS() EQ 0 THEN aspectRatio = 1.0

IF aspectRatio EQ 0 THEN BEGIN
   MESSAGE, 'Aspect Ratio of 0. Changing to 1...', /Informational
   aspectRatio = 1.0
ENDIF

s = SIZE(aspectRatio)
IF s(s(0)+1) NE 4 THEN $
   MESSAGE, 'Aspect Ratio is not a FLOAT. Take care...', /Informational

   ; Check for margins.

IF N_ELEMENTS(margin) EQ 0 THEN margin = 0.15

   ; Error checking.

IF margin LT 0 OR margin GE 0.5 THEN $
   MESSAGE, 'The MARGIN keyword value must be between 0.0 and 0.5.'

   ; Calculate the aspect ratio of the current window.

IF N_Elements(wAspectRatio) EQ 0 THEN wAspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE

   ; Calculate normalized positions in window.

IF (aspectRatio LE wAspectRatio) THEN BEGIN
   xstart = margin
   ystart = 0.5 - (0.5 - margin) * (aspectRatio / wAspectRatio)
   xend = 1.0 - margin
   yend = 0.5 + (0.5 - margin) * (aspectRatio / wAspectRatio)
ENDIF ELSE BEGIN
   xstart = 0.5 - (0.5 - margin) * (wAspectRatio / aspectRatio)
   ystart = margin
   xend = 0.5 + (0.5 - margin) * (wAspectRatio / aspectRatio)
   yend = 1.0 - margin
ENDELSE

position = [xstart, ystart, xend, yend]

RETURN, position
END;FUNCTION PSWINDOW_ASPECT
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FUNCTION PSWINDOW, LANDSCAPE=landscape, CM=cm, MARGIN=margin, $
   PageSize=pagesize, Printer=printer, Fudge=fudge, XFudge=xfudge, YFudge=yfudge, $
   DefaultPageSize=defaultpagesize, aspectpage=aspectpage

   ; Set up default values.

landscape = Keyword_Set(landscape)
cm = Keyword_Set(cm)
printer = Keyword_Set(printer)
inches = 1

   ; Set up printer fudge factors, if necessary.

IF N_Elements(fudge) NE 0 THEN BEGIN
   xfudge = fudge
   yfudge = fudge
ENDIF
IF N_Elements(xfudge) EQ 0 THEN xfudge = 0.0
IF N_Elements(yfudge) EQ 0 THEN yfudge = 0.0

   ; Get the page size.

ptype=size(pagesize,/type)
if ptype eq 0 or ptype eq 7 then begin
    IF N_Elements(pagesize) EQ 0 THEN pagesize = 'DEFAULT' $
       ELSE pagesize = StrUpCase(pagesize)
    CASE pagesize OF
       'LETTER': BEGIN
          shortside = 8.5
          longside = 11.0
          ENDCASE
       'LEGAL': BEGIN
          shortside = 8.5
          longside = 14.0
          ENDCASE
        'A4': BEGIN
          shortside = 8.27
          longside = 11.7
          ENDCASE
        ELSE: BEGIN
            if n_elements(defaultpagesize) eq 2 then begin
                shortside = float(defaultpagesize[0]<defaultpagesize[1])
                longside = float(defaultpagesize[1]>defaultpagesize[0])
            endif else begin
                shortside = 8.5
                longside = 11.0
            endelse
          ENDCASE
    ENDCASE
endif else begin
    shortside = float(pagesize[0]<pagesize[1])
    longside = float(pagesize[1]>pagesize[0])
    if KEYWORD_SET(cm) then inches = 0
    cm=0
endelse

   ; Need measurements in centimeters?

IF KEYWORD_SET(cm) THEN BEGIN
      shortside = shortside * 2.54
      longside = longside * 2.54
      inches = 0
ENDIF

   ; Determine the margin of the window on the page.

IF N_ELEMENTS(margin) EQ 0 THEN margin=0.05

   ; Get the aspect ratio of the current display window. Aspect ratio
   ; is ratio of xsize/ysize.

aspectRatio = FLOAT(!D.Y_VSIZE) / !D.X_VSIZE

   ; Get the aspect ratio of the page.

IF Keyword_Set(landscape) THEN pAspectRatio = shortside / longside $
   ELSE pAspectRatio = longside / shortside
   
   ; Aspect ratio of the page
IF Keyword_Set(aspectpage) THEN aspectRatio = pAspectRatio

   ; Get the position on the page for this window.

pos = PSWindow_Aspect(aspectRatio, Margin=margin, WindowAspect=pAspectRatio)

   ; Convert normalized position coordinates to size units.

IF KEYWORD_SET(landscape) THEN BEGIN
   IF printer THEN BEGIN
      xsize = (pos[2] - pos[0]) * longside
      ysize = (pos[3] - pos[1]) * shortside
      yoffset = pos[1] * shortside - yfudge
      xoffset = pos[0] * longside - xfudge
      landscape = 1
      portrait = 0
   ENDIF ELSE BEGIN
      xsize = (pos[2] - pos[0]) * longside
      ysize = (pos[3] - pos[1]) * shortside
      xoffset = pos[1] * shortside
      yoffset = longside - (pos[0] * longside)
      landscape = 1
      portrait = 0
   ENDELSE
ENDIF ELSE BEGIN
   xsize = (pos[2] - pos[0]) * shortside
   ysize = (pos[3] - pos[1]) * longside
   IF printer THEN BEGIN
      xoffset = pos[0] * shortside - xfudge
      yoffset = pos[1] * longside - yfudge
   ENDIF ELSE BEGIN
      xoffset = pos[0] * shortside
      yoffset = pos[1] * longside
   ENDELSE
   landscape = 0
   portrait = 1
ENDELSE

   ; Return the proper DEVICE data structure.

RETURN, {PSWINDOW_STRUCT, XSIZE:xsize, YSIZE:ysize, $
   XOFFSET:xoffset, YOFFSET:yoffset, INCHES:inches, $
   PORTRAIT:portrait, LANDSCAPE:landscape}

END;FUNCTION PSWINDOW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro insert_cgm_bin, infile, outfile,xsize=xsize,ysize=ysize

; scale xsize and ysize to maximum
xsize=!d.x_size
ysize=!d.y_size ; 0.7*!d.x_size
mult=32767./(xsize>ysize)
xsize=round(xsize*mult)
ysize=round(ysize*mult)

command=uint([8392,0,0,0>xsize<32767,0>ysize<32767]) ; VDCEXT xll yll xur yur (number from 0 to 32767)
insert_length = 10  ; in bytes
insert_string=bytarr(insert_length)
insert_string[[0,2,4,6,8]]=ishft(command and 'ff00'XU,-8)
insert_string[[1,3,5,7,9]]=command and '00ff'XU

if openu_safe(luni, infile) then return
if openw_safe(luno, outfile) then begin
    free_lun, luni
    return
endif

; read first 48 bytes

nbytes = 48L
buf1 = bytarr(48)
readu, luni, buf1
writeu, luno, buf1

; insert sequence

writeu, luno, insert_string
nbytes += insert_length

; copy rest of the file

buffer_length = 1 ; in bytes
buffer = bytarr(buffer_length)
while not eof(luni) do begin
  readu, luni, buffer
  writeu, luno, buffer
  nbytes += buffer_length
endwhile

free_lun, luni
free_lun, luno

end;pro insert_cgm_bin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SaveImageFilter,copy=copy
if keyword_set(copy) then return,['*.jpg;*.jpeg', '*.tif;*.tiff', '*.png', '*.bmp','*.gif','*.ppm','*.srf'] $
else return,['*.jpg;*.jpeg','*.ps','*.eps', '*.tif;*.tiff', '*.png', '*.bmp','*.gif','*.ppm','*.srf','*.cgm']
end;function SaveImageFilter
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SaveRGBImage,path,img,OutID=OutID

; Expect true=1

CATCH, Error_status
IF Error_status NE 0 THEN begin
    if keyword_set(OutID) then printw,OutID,path+' NOT saved.'
    return
endif

res=CutPath(path,ext=ext)
ext=strupcase(strmid(ext,1))

case ext of
'JPG':WRITE_JPEG, path, img,QUALITY=100,true=1
'JPEG':WRITE_JPEG, path, img,QUALITY=100,true=1
'TIF':write_tiff_safe,path,reverse(img,3)
'TIFF':write_tiff_safe,path,reverse(img,3)
'PNG':WRITE_PNG,path,img
'BMP':WRITE_BMP, path, img,/RGB
'GIF':WRITE_GIF,path,img ; NO COLORS
'PPM':WRITE_PPM, path,img
'SRF':WRITE_SRF, path,img
else:return
endcase

printw,OutID,path+' saved.'
end;pro SaveRGBImage
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SaveImage,path,file,procname,list,$
    noprompt=noprompt,forcereplot=forcereplot,decomposed=decomposed,$
    object=object,OutID=OutID,_EXTRA = _extra

; necessary: path, file
; necessary for PS,EPS and CGM: procname,list

if keyword_set(noprompt) then path+=file $
else path=SelFile(path,file,['*.jpg;*.jpeg','*.ps','*.eps', '*.tif;*.tiff', '*.png', '*.bmp','*.gif','*.ppm','*.srf','*.cgm'] ,'Save image ....')
if path EQ '' then return

CATCH, Error_status
IF Error_status NE 0 THEN begin
    if keyword_set(OutID) then printw,OutID,path+' NOT saved.'
    return
endif

res=CutPath(path,ext=ext)
ext=strupcase(strmid(ext,1))

delstruct,['PS','CGM'],_extra
if keyword_set(forcereplot) and ext ne 'CGM' and ext ne 'PS' and ext ne 'EPS' then begin
    if obj_valid(object) then CALL_METHOD, procname,object,list,_extra $
    else CALL_PROCEDURE, procname,list,_EXTRA = _extra
endif

case ext of
'JPG':WRITE_JPEG, path, tvrd(true=1),QUALITY=100,true=1
'JPEG':WRITE_JPEG, path, tvrd(true=1),QUALITY=100,true=1
'TIF':write_tiff_safe,path,tvrd(true=1,/order)
'TIFF':write_tiff_safe,path,tvrd(true=1,/order)
'PNG':WRITE_PNG,path,tvrd(true=1)
'BMP':WRITE_BMP, path, tvrd(true=1),/RGB
'GIF':WRITE_GIF,path,tvrd() ; NO COLORS
'PPM':WRITE_PPM, path,tvrd(true=1)
'SRF':WRITE_SRF, path
'CGM':    begin
        tmp=cutpath(path,path=dir)
        tmpcgmfile = dir+'tmp.cmg'

        ; Prepare device
        old_dev = !d.name
        ;old_fancy = !fancy
        set_plot, 'CGM'
        ;!fancy = 0
        device, /binary, file=tmpcgmfile;, colors=255, set_character_size=[1000,1000]
        ; Plot
        addstruct,'CGM',1,_extra
        if obj_valid(object) then CALL_METHOD, procname,object,list,_extra $
        else CALL_PROCEDURE, procname,list,_EXTRA = _extra
        xsize=!d.x_size
        ysize=!d.y_size

        ; Reset device
        device,/close
        ;!fancy = old_fancy
        set_plot, old_dev

        ; Finish
        insert_cgm_bin,tmpcgmfile, path,xsize=xsize,ysize=ysize
        FILE_DELETE,tmpcgmfile,/quiet

        endcase
else:     begin
        case ext of
        'PS':ENCAPSULATED=0b
        'EPS':ENCAPSULATED=1b
        else: return
        endcase

        ; Prepare device
        old_dev = !d.name
        old_fancy = !fancy
        a=PSWINDOW(_EXTRA = _extra)
        set_plot, 'PS'
        !fancy = 0
        device, file=path,ENCAPSULATED=ENCAPSULATED,/isolatin1,$
            /COLOR,BITS_PER_PIXEL=8,decomposed=decomposed,_EXTRA=a

        ; Plot
        addstruct,'PS',1,_extra
        if obj_valid(object) then CALL_METHOD, procname,object,list,_extra $
        else CALL_PROCEDURE, procname,list,_EXTRA = _extra

        ; Reset device
        device,/close
        !fancy = old_fancy
        set_plot, old_dev

        endcase
endcase

printw,OutID,path+' saved.'
end;pro SaveImage
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadEDF,path2,norm=norm

if OPENR_safe(lun, path2) then return,0

; Read Header
line='dummy'
bnorm=keyword_set(norm)
if bnorm then bnorm=norm.bnorm
while strmid(line,strlen(line)-1,1) ne '}' do begin
    readf,lun,line
    
    if strmid(line,0,8) eq 'DataType' then begin
        p1=strpos(line,'=', /REVERSE_SEARCH)+1
        p2=strpos(line,';', /REVERSE_SEARCH)
        DataType=strlowcase(strtrim(strmid(line,p1,p2-p1),2))
    endif
    if strmid(line,0,9) eq 'ByteOrder' then begin
        LEndian = strpos(line,'LowByteFirst') ne -1
    endif
    if strmid(line,0,5) eq 'Dim_1' then begin
        p1=strpos(line,'=', /REVERSE_SEARCH)+1
        p2=strpos(line,';', /REVERSE_SEARCH)
        xs=fix(strmid(line,p1,p2-p1))
    endif
    if strmid(line,0,5) eq 'Dim_2' then begin
        p1=strpos(line,'=', /REVERSE_SEARCH)+1
        p2=strpos(line,';', /REVERSE_SEARCH)
        ys=fix(strmid(line,p1,p2-p1))
    endif
    if bnorm then begin
        if strmid(line,0,strlen(norm.int)) eq norm.int then begin
            p1=strpos(line,'=', /REVERSE_SEARCH)+1
            p2=strpos(line,';', /REVERSE_SEARCH)
            norm_int=float(strmid(line,p1,p2-p1))
        endif
        if strmid(line,0,strlen(norm.gain)) eq norm.gain then begin
            p1=strpos(line,'=', /REVERSE_SEARCH)+1
            p2=strpos(line,';', /REVERSE_SEARCH)
            norm_gain=float(strmid(line,p1,p2-p1))
        endif
        if strmid(line,0,strlen(norm.offset)) eq norm.offset then begin
            p1=strpos(line,'=', /REVERSE_SEARCH)+1
            p2=strpos(line,';', /REVERSE_SEARCH)
            norm_offset=float(strmid(line,p1,p2-p1))
        endif
    endif
endwhile

if n_elements(xs) eq 0 or n_elements(ys) eq 0 $
or n_elements(DataType) eq 0 or n_elements(LEndian) eq 0 then begin
    FREE_LUN, lun
    return,0
endif 

; Data type
stat=FSTAT(lun)
nbytes=(stat.SIZE-stat.CUR_PTR)/(float(xs)*ys)
if nbytes mod 1 ne 0 then begin
    FREE_LUN, lun
    return,0
endif 

case DataType of
'signedshort'    : signed=1b
'unsignedshort'  : signed=0b
'signedinteger'  : signed=1b
'unsignedinteger': signed=0b
'signedlong'     : signed=1b
'unsignedlong'   : signed=0b
else:
endcase

if n_elements(signed) ne 0 then $
case nbytes of
1: DataType=signed?'signedbyte':'unsignedbyte'
2: DataType=signed?'signedshort':'unsignedshort'
4: DataType=signed?'signedinteger':'unsignedinteger'
8: DataType=signed?'signedlong':'unsignedlong'
else:
endcase

case DataType of
'signedbyte'     :  tiff =BYTARR(xs,ys)
'unsignedbyte'   :  tiff =BYTARR(xs,ys)
'signedshort'    :  tiff =INTARR(xs,ys)
'unsignedshort'  :  tiff =UINTARR(xs,ys)
'signedinteger'  :  tiff =LONARR(xs,ys)
'unsignedinteger':  tiff =ULONARR(xs,ys)
'signedlong'     :  tiff =LON64ARR(xs,ys)
'unsignedlong'   :  tiff =ULON64ARR(xs,ys)
'floatvalue'     :  tiff =FLTARR(xs,ys)
'float'          :  tiff =FLTARR(xs,ys)
'doublevalue'    :  tiff =DBLARR(xs,ys)
'double'         :  tiff =DBLARR(xs,ys)
else             :  return,0
endcase

; Read image
READU, lun, tiff
SWAP_ENDIAN_INPLACE, tiff,SWAP_IF_BIG_ENDIAN =LEndian,SWAP_IF_LITTLE_ENDIAN = LEndian eq 0
FREE_LUN, lun

; Normalize
if bnorm then begin
    mult=1.
    add=0.
    if n_elements(norm_int) ne 0 then mult*=norm_int
    if n_elements(norm_gain) ne 0 then mult*=norm_gain
    if n_elements(norm_offset) ne 0 then add+=norm_offset
    tiff/=mult+add
endif

return,tiff
end;function ReadEDF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WriteEDF,path,pdata

type=size(pdata,/type)
bptr=type eq 10
if bptr then begin
    type=size(*pdata,/type)
    dim=size(*pdata,/dim)
endif else dim=size(pdata,/dim)
if n_elements(dim) ne 2 then begin
    if n_elements(dim) eq 1 then dim=[dim,1] $
    else return,0b
endif

case type of
1:    begin
    DataType='UnsignedByte'
    nbytes=1
    endcase
2:    begin
    DataType='SignedShort'
    nbytes=2
    endcase
12:    begin
    DataType='UnsignedShort'
    nbytes=2
    endcase
3:    begin
    DataType='SignedInteger'
    nbytes=4
    endcase
13:    begin
    DataType='UnsignedInteger'
    nbytes=4
    endcase
14:    begin
    DataType='SignedLong'
    nbytes=8
    endcase
15:    begin
    DataType='UnsignedLong'
    nbytes=8
    endcase
4:    begin
    DataType='FloatValue'
    nbytes=4
    endcase
5:    begin
    DataType='DoubleValue'
    nbytes=8
    endcase
else: return,0b
endcase

tmp=CutPath(path,path=dir,ext=suffix,file=prefix)
nr=0
tmp=strsplit(prefix,'_',/extract,count=ct)
if ct gt 1 then begin
    for i=1,ct-2 do tmp[0]+='_'+tmp[i]
    prefix=tmp[0]+'_'
    nr=long(tmp[ct-1])
endif

if openw_safe(lun, path) then return,0b

LF=string(byte(10))
header='{'+LF
header+=string('HeaderID','= ','EH:000001:000000:000000',' ;',format='(A-15,A,A,A)')+LF
header+=string('Image','= ','1',' ;',format='(A-15,A,A,A)')+LF
header+=string('ByteOrder','= ','LowByteFirst',' ;',format='(A-15,A,A,A)')+LF
header+=string('DataType','= ',DataType,' ;',format='(A-15,A,A,A)')+LF
header+=string('Dim_1','= ',dim[0],' ;',format='(A-15,A,I0,A)')+LF
header+=string('Dim_2','= ',dim[1],' ;',format='(A-15,A,I0,A)')+LF
header+=string('Size','= ',dim[0]*dim[1]*nbytes,' ;',format='(A-15,A,I0,A)')+LF
header+=string('count_time','= ','Na',' ;',format='(A-15,A,A,A)')+LF
header+=string('point_no','= ',0,' ;',format='(A-15,A,I0,A)')+LF
header+=string('preset','= ','Na',' ;',format='(A-15,A,A,A)')+LF
header+=string('col_end','= ',dim[0]-1,' ;',format='(A-15,A,I0,A)')+LF
header+=string('col_beg','= ',0,' ;',format='(A-15,A,I0,A)')+LF
header+=string('row_end','= ',dim[1]-1,' ;',format='(A-15,A,I0,A)')+LF
header+=string('row_beg','= ',0,' ;',format='(A-15,A,I0,A)')+LF
header+=string('col_bin','= ',1,' ;',format='(A-15,A,I0,A)')+LF
header+=string('row_bin','= ',1,' ;',format='(A-15,A,I0,A)')+LF
header+=string('time','= ',systime(),' ;',format='(A-15,A,A,A)')+LF
header+=string('time_of_day','= ',systime(/seconds),' ;',format='(A-15,A,D0.0,A)')+LF
header+=string('dir','= ',dir,' ;',format='(A-15,A,A,A)')+LF
header+=string('suffix','= ',suffix,' ;',format='(A-15,A,A,A)')+LF
header+=string('prefix','= ',prefix,' ;',format='(A-15,A,A,A)')+LF
header+=string('run','= ',nr,' ;',format='(A-15,A,I0,A)')+LF
header+=string('title','= ','XRDUA version '+version(),' ;',format='(A-15,A,A,A)')+LF
addend='}'+LF
nspaces=strlen(header)+strlen(addend)
nspaces=ceil(nspaces/512.)*512-nspaces
if nspaces gt 0 then header+=string(replicate(32b,nspaces))
header+=addend

writeu,lun,header
if bptr then $
    writeu,lun,SWAP_ENDIAN(*pdata, /SWAP_IF_BIG_ENDIAN) $
else $
    writeu,lun,SWAP_ENDIAN(pdata, /SWAP_IF_BIG_ENDIAN)
free_lun,lun

return,1b
end;function WriteEDF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_d, name,header=header

CATCH, Error_status
IF Error_status NE 0 THEN return,0

 line = ''
 header=line
 if openr_safe(lun,name) then return,0
 i=0
 repeat begin
    readf,lun,line
    header=[header,line]
    i=i+1
 endrep until ((strpos(line,'Intensity') GE 0)) OR eof(lun)

 if eof(lun) then return, 0
 j=0
 readf,lun,line
 reads,line,j
 data=fltarr(2,j)

 readf,lun,data

 free_lun, lun
 return, data
end;function read_d
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_spe,path,error=error

error=1B

CATCH, Error_status
IF Error_status NE 0 THEN return,0

if openr_safe(lun, path) then return,0
line = ''
while(strmid(line,0,6) ne '$DATA:') and (not eof(lun)) do readf, lun, line
if eof(lun) then return,0
n0=0
n1 = 0
readf, lun, n0, n1
n_channels=n1
spe = fltarr(n_channels)
readf, lun, spe
free_lun, lun

error=0b
return,[findgen(1,n_channels),transpose(spe)]
end;function read_spe
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro write_spe,file,spe,fileorg,error=error
error=1b
if openw_safe(lunj,file) then return
printf,lunj,'$SPEC_ID:'
printf,lunj,fileorg
printf,lunj,'$DATA:'
printf,lunj,'1  '+strtrim(n_elements(spe),2)
printf,lunj,spe,format='(10(I0," "))'
free_lun,lunj
error=0b
end;pro write_spe
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_fio,path,error=error,nspec=nspec,nchan=nchan,header=header,nochan=nochan,counters=counters,hcounters=hcounters

error=1B
CATCH, Error_status
IF Error_status NE 0 THEN return,0

if openr_safe(lun,path) then return,0

; Find comment block
line=''
block='! Comments'
nblock=strlen(block)
while (strmid(line,0,nblock) ne block) or eof(lun) do readf,lun,line
if eof(lun) then return,0

; Find data block
block='! Data'
nblock=strlen(block)
buffer=''
while (strmid(line,0,nblock) ne block) or eof(lun) do begin
    readf,lun,line
    buffer=[buffer,line]
endwhile
if eof(lun) then return,0

; Find data size
readf,lun,line ; " %d"
block=' Col'
nblock=strlen(block)
nspec=0
readf,lun,line
header=''
while (strmid(line,0,nblock) eq block) and ~eof(lun) do begin
    header=[header,(strsplit(line,' ',/extract))[2]]
    nspec++
    readf,lun,line
endwhile
if nspec eq 0 then return,0
header=header[1:*]

; Read data
arr=fltarr(nspec)
reads,line,arr
spe=arr
while not eof(lun) do begin
    readf,lun,arr
    spe=[[spe],[arr]]
endwhile

; Finish
free_lun,lun

; Counters
C2=0.
C4=0.
IDORIS=0.
ICR=fltarr(nspec)
OCR=ICR
LIVETIME=ICR
mcanames=stregex(strlowcase(header),'mca[0-9]+',/extract)
buffer=strlowcase(buffer)
hcounters=['C2','C4','IDORIS',mcanames+'_ICR',mcanames+'_OCR',mcanames+'_LT']
for i=0l,n_elements(buffer)-1 do begin
    line=buffer[i]
    p=strpos(line,'c2 ')
    if p ne -1 then C2=float(stregex(strmid(line,p+2),' -*[0-9.]+',/extract))
    p=strpos(line,'c4 ')
    if p ne -1 then C4=float(stregex(strmid(line,p+2),' -*[0-9.]+',/extract))
    p=strpos(line,'doris ')
    if p ne -1 then IDORIS=float(stregex(strmid(line,p+5),' -*[0-9.]+',/extract))
    
    p=strpos(line,'livetime ')
    if p ne -1 then begin
        ind=where(mcanames eq stregex(line,'mca[0-9]+',/extract),ct)
        if ct eq 1 then LIVETIME[ind]=float(stregex(strmid(line,p+8),' -*[0-9.]+',/extract))
    endif
    p=strpos(line,'icr ')
    if p ne -1 then begin
        ind=where(mcanames eq stregex(line,'mca[0-9]+',/extract),ct)
        if ct eq 1 then ICR[ind]=float(stregex(strmid(line,p+3),' -*[0-9.]+',/extract))
    endif
    p=strpos(line,'ocr ')
    if p ne -1 then begin
        ind=where(mcanames eq stregex(line,'mca[0-9]+',/extract),ct)
        if ct eq 1 then OCR[ind]=float(stregex(strmid(line,p+3),' -*[0-9.]+',/extract))
    endif
endfor

; Return
counters=[C2,C4,IDORIS,ICR,OCR,LIVETIME]
nchan=n_elements(spe)/nspec
error=0b
if keyword_set(nochan) then return, spe $
else return,[findgen(1,nchan),spe]
end;function read_fio
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro write_fio,file,spe,fileorg,error=error

error=1b
if openw_safe(lun,file) then return

printf,lun,'!'
printf,lun,'! Comments '
printf,lun,'%c '
printf,lun,' ',fileorg
printf,lun,'! '
printf,lun,'! Data '
printf,lun,'%d '

s=size(spe)
tmp=indgen(s[1])
printf,lun,tmp,tmp,format='(" Col ",I0," ROI",I0,"  FLOAT ")'
printf,lun,spe

free_lun, lun
error=0b
end;pro write_fio
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_mca,path,error=error

error=1b
CATCH, Error_status
IF Error_status NE 0 THEN return,0

if openr_safe(lun, path) then return,0
line = ''

while strmid(line,0,6) ne '#@MCA ' and not eof(lun) do readf,lun,line
nline=fix((strsplit(line,'#@MCA% ',/extract))[0])
if nline eq 0 then return,0

point_lun,lun,0
while strmid(line,0,7) ne '#@CHANN' and not eof(lun) do readf,lun,line
nch=long((strsplit(line,' ',/extract))[1])
if nch eq 0 then return,0

point_lun,lun,0
while strmid(line,0,2) ne '@A' do readf,lun,line

i = 0
imax = nch-1
d = lonarr(1,nch)
d16 = lonarr(nline)
while i lt nch do begin
    if i eq 0 then line = strmid(line,2,strlen(line)-2) $
    else readf, lun,line
    
    reads, line, d16
    d[i:i+nline-1] = d16
    i += nline
    
    if eof(lun) then break
endwhile 

free_lun, lun
error=0b
return,[findgen(1,nch),d[0,0:nch-1]]
end;function read_mca
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro write_mca,file,spe,fileorg,error=error
error=1b
if openw_safe(lun,file) then return

nspe=n_elements(spe)
nadd=2^ceil(alog(nspe)/alog(2.))-nspe
if nadd gt 0 then begin
    spe=[spe,make_array(nadd,value=spe[0]*0)]
    nspe+=nadd
endif

format='(A,"'+string(10b)+'")'

buffer=''
buffer+=string('#F '+fileorg,format=format)
buffer+=string('#S 1 Single spectrum',format=format)
buffer+=string('#@MCA %16C',format=format)
buffer+=string('#@CHANN '+strtrim(nspe,2)+' 0 '+strtrim(nspe-1,2)+' 1',format=format)

t=strcompress(string(spe,format='(16I)'))
buffer+='@A'+strjoin(t,'\'+string(10b))+string(10b)
writeu,lun,buffer

free_lun, lun
error=0b
end;pro write_mca
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_chi, name,XTITLE=XTITLE,header=header,error=error,ROIAzimuth=ROIAzimuth,berror=berror,bring=bring

error=1B
berror=0b
bring=0b

CATCH, Error_status
IF Error_status NE 0 THEN return,0

 line = ''
 header=line
 if openr_safe(lun,name) THEN return,0

 readf,lun,line
 header=[header,line]
 readf,lun,line
 header=[header,line]
 XTITLE=line
 readf,lun,line
 header=[header,line]
 readf,lun,line

 j=0l
 reads,line,j
 data=fltarr(2,j)

 readf,lun,data

 ROIAzimuth=[0,2*!pi]
 while (strmid(line,0,1) ne '$') and (not eof(lun)) do readf,lun,line
 if (not eof(lun)) then readf,lun,ROIAzimuth
 
 point_lun,lun,0
 while (strmid(line,0,8) ne '$YError:') and (not eof(lun)) do readf,lun,line
 if (not eof(lun)) then begin
     tmp=fltarr(1,j)
     readf,lun,tmp
     data=[data,tmp]
     berror=1b
 endif
 
 point_lun,lun,0
 while (strmid(line,0,7) ne '$Ring%:') and (not eof(lun)) do readf,lun,line
 if (not eof(lun)) then begin
     tmp=fltarr(1,j)
     readf,lun,tmp
     data=[data,tmp]
     bring=1b
 endif

 free_lun, lun
 error=0B
 return, data
end;function read_chi
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_plt, name,header=header,error=error,ROIAzimuth=ROIAzimuth

error=1b

CATCH, Error_status
IF Error_status NE 0 THEN return,0

header=''
if openr_safe(lun,name) THEN BEGIN
;message,!ERR_STRING,/info
error=1B
return,0
ENDIF

line='!@'
while (strmid(line,0,2) eq '!@') and (not eof(lun)) do begin
    readf,lun,line
    header=[header,line]
    pos=strpos(line,'Gamma:')
    if pos ne -1 then begin
        str=strmid(line,pos+6)
        str=str_sep(str,'to')
        ROIAzimuth=[0,2*!pi]
        reads,str,ROIAzimuth
        ROIAzimuth=ROIAzimuth/180*!pi
    endif
endwhile

header=header[0:n_elements(header)-2]

data=fltarr(2)
j=data
while not eof(lun) do begin
    readf,lun,j
    data=[[data],[j]]
endwhile
free_lun, lun

if n_elements(data) eq 2 then return,0

data=data[*,1:*]
error=0B
return, data
end;function read_plt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mread_ascii, name,header=header,error=error,asciitemplate=asciitemplate

error=1b

CATCH, Error_status
IF Error_status NE 0 THEN return,0

header=''
line=''
if openr_safe(lun,name) then return,0

if n_elements(asciitemplate) eq 0 then begin
    data=fltarr(2)
    line=data
    readf,lun,data
    CATCH, Error_status
    IF Error_status eq 0 THEN $
    while not eof(lun) do begin
        readf,lun,line
        data=[[data],[line]]
    endwhile
endif else begin
    if asciitemplate.rowskip ne 0 then begin
        header=strarr(asciitemplate.rowskip)
        readf,lun,header,format='(A)'
    endif
    
    line=''
    readf,lun,line
    data=strsplit(line,asciitemplate.delimiter,/extract)
    CATCH, Error_status
    IF Error_status eq 0 THEN $
    while not eof(lun) do begin
        readf,lun,line
        data=[[data],[strsplit(line,asciitemplate.delimiter,/extract)]]
    endwhile
    
    data=data[asciitemplate.ind,*]
    strreplace,data,',','.'
    data=float(data)
    
endelse


free_lun,lun

error=0B
return, data
end;function mread_ascii
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_drn, name,header=header,error=error

error=1b

CATCH, Error_status
IF Error_status NE 0 THEN return,0

header=''
if openr_safe(lun,name) then return,0

data=fltarr(2)
j=data
line=''
while not eof(lun) do begin
    readf,lun,line
    line=byte(line)
    ind=where(line eq 44,ct)
    if ct ne 0 then line[ind]=46
    line=string(line)
    reads,line,j
    data=[[data],[j]]
endwhile
free_lun, lun

if n_elements(data) eq 2 then return,0

data=data[*,1:*]
error=0B
return, data
end;function read_drn
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_x_y, name,header=header,error=error

error=1b

CATCH, Error_status
IF Error_status NE 0 THEN return,0

header=''
if openr_safe(lun,name) then return,0

data=fltarr(2)
line=data
while not eof(lun) do begin
    readf,lun,line
    data=[[data],[line]]
endwhile

free_lun,lun

error=0B
return, data[*,1:*]
end;function read_x_y
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_q, name,header=header,error=error

error=1b

CATCH, Error_status
IF Error_status NE 0 THEN return,0

header=''
if openr_safe(lun,name) then return,0

data=fltarr(2)
line=data
while not eof(lun) do begin
    readf,lun,line
    data=[[data],[line]]
endwhile

free_lun,lun

error=0B
return, data[*,1:*]
end;function read_q
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_raw101,lun,header

; ----Read header (712 bytes, 4 bytes already read)----
in=bytarr(4)
readu,lun,in

in=0ul
readu,lun,in
case in of
1:header={file_status:'done'}
2:header={file_status:'active'}
3:header={file_status:'aborted'}
4:header={file_status:'interrupted'}
else:header={file_status:'unknown'}
endcase

range_cnt=0ul
readu,lun,range_cnt

in=bytarr(10)
readu,lun,in
header=create_struct(header,'MEASURE_DATE',string(in))

in=bytarr(10)
readu,lun,in
header=create_struct(header,'MEASURE_TIME',string(in))

in=bytarr(72)
readu,lun,in
header=create_struct(header,'USER',string(in))

in=bytarr(218)
readu,lun,in
header=create_struct(header,'SITE',string(in))

in=bytarr(60)
readu,lun,in
header=create_struct(header,'SAMPLE_ID',string(in))

in=bytarr(160)
readu,lun,in
header=create_struct(header,'COMMENT',string(in))

in=bytarr(62)
readu,lun,in

in=bytarr(4)
readu,lun,in
header=create_struct(header,'ANODE_MATERIAL',string(in))

in=bytarr(4)
readu,lun,in

in=0d
readu,lun,in
header=create_struct(header,'ALPHA_AVERAGE',in)

in=0d
readu,lun,in
header=create_struct(header,'ALPHA1',in)

in=0d
readu,lun,in
header=create_struct(header,'ALPHA2',in)

in=0d
readu,lun,in
header=create_struct(header,'BETA',in)

in=0d
readu,lun,in
header=create_struct(header,'ALPHA_RATIO',in)

in=bytarr(8)
readu,lun,in

in=0.
readu,lun,in
header=create_struct(header,'measurement time',in)

in=bytarr(44)
readu,lun,in

; ----Range header----
range_cnt<=1 ; Only read the first
for i=0l,range_cnt-1 do begin
    headerlen=0ul
    steps=0ul
    start_theta=0d
    start_2theta=0d
    readu,lun,headerlen
    readu,lun,steps
    readu,lun,start_theta
    readu,lun,start_2theta
    
    in=bytarr(76)
    readu,lun,in
    
    in=0.
    readu,lun,in
    header=create_struct(header,'HIGH_VOLTAGE',in)
    in=0.
    readu,lun,in
    header=create_struct(header,'AMPLIFIER_GAIN',in)
    in=0.
    readu,lun,in
    header=create_struct(header,'DISCRIMINATOR_1_LOWER_LEVEL',in)
    
    in=bytarr(64)
    readu,lun,in
    
    step_size=0d
    readu,lun,step_size
    
    in=bytarr(8)
    readu,lun,in
    
    in=0.
    readu,lun,in
    header=create_struct(header,'TIME_PER_STEP',in)
    
    in=bytarr(12)
    readu,lun,in

    in=0.
    readu,lun,in
    header=create_struct(header,'ROTATION_SPEED_rpm',in)
    
    in=bytarr(12)
    readu,lun,in

    in=0ul
    readu,lun,in
    header=create_struct(header,'GENERATOR_VOLTAGE',in)
    in=0ul
    readu,lun,in
    header=create_struct(header,'GENERATOR_CURRENT',in)
    
    in=bytarr(8)
    readu,lun,in
    
    in=0d
    readu,lun,in
    header=create_struct(header,'USED_LAMBDA',in)
    
    in=bytarr(8)
    readu,lun,in
    
    supplementary_headers_size=0ul
    readu,lun,supplementary_headers_size
    
    in=bytarr(44)
    readu,lun,in
    
    point_lun,-lun,p
    if p-712 ne headerlen then continue
    
    if supplementary_headers_size ne 0 then begin
        in=bytarr(supplementary_headers_size)
        readu,lun,in
    endif
    
    data=fltarr(steps)
    readu,lun,data
    data=[start_2theta+step_size*lindgen(1,steps),reform(data,1,steps,/overwrite)]
endfor

if n_elements(data) eq 0 then return,0
return,data
end;function read_raw101
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_raw, name,header=header,error=error
; Siemens/Bruker RAW

error=1b
header=''

CATCH, Error_status
IF Error_status NE 0 THEN begin
    if n_elements(lun) ne 0 then free_lun,lun
    return,0
endif

if openr_safe(lun,name,/SWAP_IF_BIG_ENDIAN) then return,0

v=bytarr(4)
readu,lun,v
v=string(v)
case v of
'RAW ':
'RAW1': data=read_raw101(lun,header)
'RAW2':
else:
endcase

free_lun,lun
if n_elements(data) le 0 then return,0
error=0B
return, data
end;function read_raw
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_dat, name,header=header,error=error

error=1b

CATCH, Error_status
IF Error_status NE 0 THEN return,0

header=''
if openr_safe(lun,name) then return,0

x=fltarr(3)
readf,lun,x
x=makerange(x[0],x[2],x[1])

y=fltarr(n_elements(x))
readf,lun,y

free_lun,lun

error=0B
return, transpose([[x],[y]])
end;function read_dat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function read_nor, name,header=header,error=error

error=1b

CATCH, Error_status
IF Error_status NE 0 THEN return,0

header=''
if openr_safe(lun,name) then return,0

line='#'
while strmid(line,0,1) eq '#' and not eof(lun) do readf,lun,line
if eof(lun) then begin
    free_lun,lun
    return,0
endif

x=fltarr(4)
reads,line,x
ret=x
while not eof(lun) do begin
    readf,lun,x
    ret=[[ret],[x]]
endwhile

free_lun,lun

error=0B
return, ret[0:1,*]
end;function read_nor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadScanDummy,top,xval=xval,spe=spe,stdspe=speerror,title=title,$
  YTITLE=YTITLE,xrange=xrange,yrange=yrange,type=type,ROIAzimuth=ROIAzimuth

tt0 = 0>float(PromptNumber('0',top,'2-theta start (degrees:)'))<180
tt1 = 0>float(PromptNumber('180',top,'2-theta end (degrees:)'))<180
tti = 0>float(PromptNumber('0.001',top,'2-theta increment (degrees:)'))<180

tmp = tt0
tt0 <= tt1
tt1 >= tmp

n = ceil((tt1-tt0)/tti)+1

ROIAzimuth = [0,2*!pi]
type = 1
yrange = [0,1]
xrange = [0,n-1]
title = 'dummy'
YTITLE='Intensity (a.u.)'

xval = tt0 + tti*findgen(n)
spe = fltarr(n)
speerror = spe

return,1b
end;function ReadScanDummy
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadScan,path,file,Format,OutID,xval=xval,spe=spe,stdspe=speerror,title=title,$
    YTITLE=YTITLE,xrange=xrange,yrange=yrange,header=header,type=type,ROIAzimuth=ROIAzimuth,$
    ringpr=ringpr,xdinr=xdinr,asciitemplate=asciitemplate

CATCH, Error_status
IF Error_status NE 0 THEN return,0

ROIAzimuth=[0,2*!pi]
berror=0b
bring=0b
sFormat=strlowcase(Format)
case sFormat of
'.d' :     begin
        m=read_d(path+file,header=header)
        type=0
        endcase
'.chi': begin
        m=read_chi(path+file,XTITLE=XTITLE,header=header,error=error,ROIAzimuth=ROIAzimuth,berror=berror,bring=bring)
        if error then begin
            printw,OutID,'Can`t open '+path+file
            return,0
        endif
        type=ChiXtitle(XTITLE,mode=4,count=ct)
        if (ct eq 0) then printw,OutID,'Unkown data type in '+path+file
        endcase
'.plt': begin
        m=read_plt(path+file,header=header,error=error,ROIAzimuth=ROIAzimuth)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.drn': begin
        m=read_drn(path+file,header=header,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.raw': begin
        m=read_raw(path+file,header=header,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.x_y':    begin
        m=read_x_y(path+file,header=header,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.q':    begin
        m=read_q(path+file,header=header,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=3
        endcase
'.asc':    begin
        m=read_x_y(path+file,header=header,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.dat':    begin
        m=read_dat(path+file,header=header,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.nor':    begin
        m=read_nor(path+file,header=header,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.prn':    begin
        m=mread_ascii(path+file,header=header,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.uxd':    begin
        m=read_x_y(path+file,header=header,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.mca':    begin
        m=read_mca(path+file,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.spe':    begin
        m=read_spe(path+file,error=error)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif
        type=1
        endcase
'.fio':    begin
        m=read_fio(path+file,error=error,nspec=nspec)
        if error then begin
            printw,OutID,'Can`t open '+path+file+ ' or no data in file.'
            return,0
        endif else begin
            if nspec ne 1 then $
                m=[m[0,*],transpose(total(m[1:*,*],1))]
        endelse
        type=1
        endcase
else:    begin
        if sFormat eq '.tiff' then sFormat='.tif'
        case sFormat of
        '.tif':    begin
                m = read_tiff(path+file,geotiff=geotiff)
                if n_elements(m) eq 1 then return,0B

                tmp=where(tag_names(geotiff) eq 'GTMODELTYPEGEOKEY',tmp2)
                if tmp2 eq 0 then begin
                    berror=0b
                    bring=0b
                endif else begin
                    berror=geotiff.GTMODELTYPEGEOKEY
                    bring=geotiff.GTRASTERTYPEGEOKEY
                endelse
                nblocks=1+berror+bring
                
                s=dimsize(m,2)
                blocklen=(s[1]-1)/nblocks
                ind1=indgen(nblocks)*blocklen+1
                ind2=ind1+blocklen-1

                if n_elements(xdinr) then begin
                    title=file+'_'+stringr(xdinr)
                    
                    spe=m[*,ind1[0]+xdinr]
                    
                    case nblocks of
                    1:    begin
                        speerror=sqrt(spe)
                        endcase
                    2:    begin
                        if berror then speerror=m[*,ind1[1]+xdinr]
                        if bring then ringpr=m[*,ind1[1]+xdinr]
                        endcase
                    3:    begin
                        speerror=m[*,ind1[1]+xdinr]
                        ringpr=m[*,ind1[2]+xdinr]
                        endcase
                    endcase
                    
                endif else begin
                    title = file

                    spe=m[*,ind1[0]:ind2[0]]
                    
                    case nblocks of
                    1:    begin
                        speerror=sqrt(spe)
                        endcase
                    2:    begin
                        if berror then speerror=m[*,ind1[1]:ind2[1]]
                        if bring then ringpr=m[*,ind1[1]:ind2[1]]
                        endcase
                    3:    begin
                        speerror=m[*,ind1[1]:ind2[1]]
                        ringpr=m[*,ind1[2]:ind2[2]]
                        endcase
                    endcase
                endelse

                type=ChiXtitle(geotiff.GEOGCITATIONGEOKEY,mode=4,count=ct)
                ROIAzimuth=[geotiff.GEOGLINEARUNITSIZEGEOKEY,$
                            geotiff.GEOGANGULARUNITSIZEGEOKEY]
    
                YTITLE='Intensity (a.u.)'
                xval=reform(m[*,0])
                xrange=[0,n_elements(xval)-1]
                s=size(spe)
                temp=min(min(spe,dimension=s[0]),y1)
                temp=max(max(spe,dimension=s[0]),y2)
                yrange=[y1[0],y2[0]]
                
                return,1B
                endcase
        '.edf':    begin
                m = readedf(path+file)
                if n_elements(m) eq 1 then return,0B
                s=DimSize(m,2)
                xval=lindgen(s[0])
                if n_elements(xdinr) then begin
                    spe=reform(m[*,1+xdinr])
                    title=file+'_'+stringr(xdinr)
                endif else begin
                    spe=reform(m[*,1:*])
                    title = file+'_sum'
                endelse
    
                type=1
                YTITLE='Intensity (a.u.)'
                xrange=[0,n_elements(xval)-1]
                s=size(spe)
                temp=min(min(spe,dimension=s[0]),y1)
                temp=max(max(spe,dimension=s[0]),y2)
                yrange=[y1[0],y2[0]]
                
                speerror=sqrt(spe)
                return,1B
                
                endcase
        else:    begin
                m=mread_ascii(path+file,header=header,error=error,asciitemplate=asciitemplate)
                if error then begin
                    printw,OutID,'No data in file or no valid file extension.'
                    return,0B
                endif
                type=1
                endelse
        endcase
        endcase
endcase

if n_elements(m) eq 1 then return,0B

s=size(m)
xval=reform(m[0,*])
ind=sort(xval)
xval=xval[ind]
spe=reform(m[1,ind])

if berror then speerror=reform(m[2,ind]) else speerror=sqrt(spe) ; Poisson statistics
if bring then ringpr=reform(m[2+berror,ind])

title = file

YTITLE='Intensity (a.u.)'
xrange=[0,n_elements(xval)-1]
temp=min(spe,y1)
temp=max(spe,y2)
yrange=[y1[0],y2[0]]

return,1B
end;fucntion ReadScan
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CHI_xval,xval,lambda,scrx,scry,a,b,dist,xtype,spe=spe,tspe=spet,stdspe=speerror,tstdspe=speerrort

n=n_elements(xval)
xvalt=make_array(n,4,value=xval[0]*0)
xvalt[*,xtype]=xval
berror=n_elements(speerror) eq n
if keyword_set(spe) then begin
    spet=make_array(n,4,value=spe[0]*0)
    spet[*,xtype]=spe
    if berror then begin
        speerrort=spet
        speerrort[*,xtype]=speerror
    endif
    case xtype of
        0:    begin
            res=BraggDtoX_I(xval,berror?[[spe],[speerror]]:spe,lambda,dist,scrx,scry,a,b,/pixm,/angledeg)
            xvalt[*,1]=res.x2
            xvalt[*,2]=res.x1
            xvalt[*,3]=res.x3
            spet[*,1]=res.y2[*,0]
            spet[*,2]=res.y1[*,0]
            spet[*,3]=res.y3[*,0]
            if berror then begin
                speerrort[*,1]=res.y2[*,1]
                speerrort[*,2]=res.y1[*,1]
                speerrort[*,3]=res.y3[*,1]
            endif
            endcase
        1:    begin
            res=BraggTtoX_I(xval,berror?[[spe],[speerror]]:spe,lambda,dist,scrx,scry,a,b,/pixm,/angledeg)
            xvalt[*,0]=res.x2
            xvalt[*,2]=res.x1
            xvalt[*,3]=res.x3
            spet[*,0]=res.y2[*,0]
            spet[*,2]=res.y1[*,0]
            spet[*,3]=res.y3[*,0]
            if berror then begin
                speerrort[*,0]=res.y2[*,1]
                speerrort[*,2]=res.y1[*,1]
                speerrort[*,3]=res.y3[*,1]
            endif
            endcase
        2:    begin
            res=BraggPtoX_I(xval,berror?[[spe],[speerror]]:spe,lambda,dist,scrx,scry,a,b,/pixm,/angledeg)
            xvalt[*,0]=res.x2
            xvalt[*,1]=res.x1
            xvalt[*,3]=res.x3
            spet[*,0]=res.y2[*,0]
            spet[*,1]=res.y1[*,0]
            spet[*,3]=res.y3[*,0]
            if berror then begin
                speerrort[*,0]=res.y2[*,1]
                speerrort[*,1]=res.y1[*,1]
                speerrort[*,3]=res.y3[*,1]
            endif
            endcase
        3:    begin
            res=BraggQtoX_I(xval,berror?[[spe],[speerror]]:spe,lambda,dist,scrx,scry,a,b,/pixm,/angledeg)
            xvalt[*,0]=res.x3
            xvalt[*,1]=res.x2
            xvalt[*,2]=res.x1
            spet[*,0]=res.y3[*,0]
            spet[*,1]=res.y2[*,0]
            spet[*,2]=res.y1[*,0]
            if berror then begin
                speerrort[*,0]=res.y3[*,1]
                speerrort[*,1]=res.y2[*,1]
                speerrort[*,2]=res.y1[*,1]
            endif
            endcase
    endcase
endif else begin
    case xtype of
        0:    begin
            res=BraggDtoX(xval,lambda,dist,scrx,scry,a,b,/pixm,/angledeg)
            xvalt[*,1]=res[*,1]
            xvalt[*,2]=res[*,0]
            xvalt[*,3]=res[*,2]
            endcase
        1:    begin
            res=BraggTtoX(xval,lambda,dist,scrx,scry,a,b,/pixm,/angledeg)
            xvalt[*,0]=res[*,1]
            xvalt[*,2]=res[*,0]
            xvalt[*,3]=res[*,2]
            endcase
        2:    begin
            res=BraggPtoX(xval,lambda,dist,scrx,scry,a,b,/pixm,/angledeg)
            xvalt[*,0]=res[*,1]
            xvalt[*,1]=res[*,0]
            xvalt[*,3]=res[*,2]
            endcase
        3:    begin
            res=BraggQtoX(xval,lambda,dist,scrx,scry,a,b,/pixm,/angledeg)
            xvalt[*,0]=res[*,2]
            xvalt[*,1]=res[*,1]
            xvalt[*,2]=res[*,0]
            endcase
    endcase
endelse
return,xvalt
end;function CHI_xval
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SavePKS,list,oD=oD
if keyword_set(oD) then begin
    if list.apeaksout eq 0 then begin
       printw,list.OutID, 'There are no fitresults'
       return,0B
    endif
    path_pos = rstrpos(list.file,'.')  ; find first dot
    if(path_pos eq -1) then path_pos = strlen(list.file)     ; no extention -> append

    ; ----Save in format for PDF-2 database?----
    Result = DIALOG_MESSAGE('Save for Pcsiwin?',/question)
    if Result eq 'No' then begin
        name = strmid(list.file,0,path_pos)+'.pks'
        path=SelFile(list.path,name,'*.pks','Save Results To....')
        if path EQ '' then return,0B
        error=CutPath(path,path=pathtemp)
        list.path=pathtemp

        ; read header
        error=ReadScan(list.path,list.file,list.Format,list.OutID,header=header)
        if n_elements(header) ne 1 then header=header[1:*]

        ;writing (all fitted peaks)
        if openw_safe(lun, path, width=120) then begin
           printw,list.OutID, 'Not saved: '+path
           return,0B
        endif

        for i=0l,n_elements(header)-1 do printf,lun,header[i]

        printf,lun,list.apeaksout
        printf,lun,'$Xval|Area|Sx|n:'
        matrix=list.peakout[*,0:list.apeaksout-1]
        ind=reverse(sort(matrix[0,*]))
        for i=0l,list.apeaks-1 do printf,lun,matrix[*,ind[i]]
        printf,lun,'$SD:'
        matrix=list.peakouts[*,0:list.apeaksout-1]
        for i=0l,list.apeaks-1 do printf,lun,matrix[*,ind[i]]
        printf,lun,'$CHI_SQR:'
        printf,lun,list.chisq

    endif else begin
        name = strmid(list.file,0,path_pos)+'.txt'
        path=SelFile(list.path,name,'*.txt','Save Results To....')
        if path EQ '' then return,0B
        error=CutPath(path,path=pathtemp)
        list.path=pathtemp

        ;writing (all fitted peaks)
        if openw_safe(lun, path, width=120) then begin
           printw,list.OutID, 'Not saved: '+path
           return,0B
        endif

        matrix=list.peakout[*,0:list.apeaksout-1]

        ; save d-spacing
        xen=CHI_xval(reform(matrix[0,*]),list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,list.xtypeorig)
        xen=xen[*,0]
        matrix[0,*]=xen

        ind=reverse(sort(matrix[0,*]))
        indf=where(finite(matrix[1,*]) eq 1)
        m=abs(max(matrix[1,indf]))
        matrix[1,*]=(matrix[1,*]/m*999)>0
        matrix=matrix[*,ind]
        for i=0l,list.apeaksout-1 do printf,lun,matrix[0,i],long(matrix[1,i]),1
    endelse

    free_lun, lun
    printw,list.OutID,path+' saved'

endif else begin
    path=SelFile(list.pathpks,list.filepks,'*.pks','Save Results To....')
    if path EQ '' then return,0B
    if openw_safe(lun, path, width=120) then  begin
       printw,list.OutID, 'Not saved: '+path
       return,0B
    endif
    printf,lun,list.pathorig+': 2D Fitresults'
    printf,lun,'$ROI: '
    printf,lun,list.sqr[3:*]
    printf,lun,list.apeaks
    printf,lun,'$volume, mx, my, sx, sy, sxy, n:'
    ;printf,lun,format='(7F15.4)',list.AA
    for i=0l,list.apeaks-1 do printf,lun,list.AA[*,i]
    printf,lun,'$SD:'
    ;printf,lun,format='(7F15.4)',list.sigma
    for i=0l,list.apeaks-1 do printf,lun,list.sigma[*,i]
    free_lun,lun
    printw,list.OutID,path+' saved'
endelse
return,1B
end;function SavePKS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WriteLOF,path,file,nr,top,filenames,OutID,header=header

savefile=SelFile(path,file,'*.lof','Save list of files to...')
if savefile eq '' then return,0B

; ----Generate new names or take the old ones----
if n_elements(filenames) eq 1 then Result='Yes' else $
Result = DIALOG_MESSAGE('Generate new file names (current ones otherwise)?',/question,/default_no)

if Result eq 'Yes' then begin
    ; ----Make BaseString----
    UPos = RSTRPOS(file, '_')
    if UPos eq -1 then UPos = RSTRPOS(file, '.')
    file=strmid(file,0,UPos)
    basestring=path+file+'*.tif'

    ;----Modify basestring----
    basestring=PromptNumber(basestring,top,'Enter filter string:')

    ; ----Some options----
    ; Check for number insert point
    SPos = RSTRPOS(basestring, '*')
    if SPos eq -1 then begin
       printw,OutID,'Do not delete the "*"'
       return,0
    endif
    ; Offset?
    Result = DIALOG_MESSAGE('Offset for file numbering: 0?',/question)
    if result eq 'No' then begin
       offset=0L
       filter=PromptNumber('0',top,'Enter numberoffset:')
       reads,filter,offset
       nr=nr+offset
    endif
    ; Make numbers in string format
    nr=MakeNumber(nr,/prompt)
    ; Make file-strings
    files=strmid(basestring,0,SPos)+nr+strmid(basestring,SPos+1)
endif else files=filenames[nr]

; ----Save Files----
if openw_safe(lun,savefile) then begin
    printw,OutID,savefile+' not saved'
    return,0B
endif

n=n_elements(files)
if keyword_set(header) then $
    for i=0l,n_elements(header)-1 do printf,lun,header[i]
printf,lun,'$Files:'
printf,lun,n
printf,lun,files,format='(A)'
free_lun,lun
printw,OutID,savefile+' saved'
return,1B

end; function WriteLOF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CleanUp_BrowserOptions, ID

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
end;pro CleanUp_BrowserOptions
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro BrowserOptions_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_value=val
val=val[0]
widget_control,ev.top,get_uvalue=tlist
case val of
'Sort by string (file_10.tif,file_2.tif,...)':   begin
              widget_control,ev.id,get_uvalue=uval
              tlist.method=uval
              widget_control,ev.top,set_uvalue=tlist
              ID=widget_info(ev.top,FIND_BY_UNAME='text')
              widget_control,ID,sensitive=0
              endcase
'Sort with sep. (file_2.tif,file_10.tif,...)': begin
              widget_control,ev.id,get_uvalue=uval
              tlist.method=uval
              widget_control,ev.top,set_uvalue=tlist
              ID=widget_info(ev.top,FIND_BY_UNAME='text')
              widget_control,ID,sensitive=1
              endcase
'Sort by number (file1_10.tif,file2_2.tif,file2_10.tif,...)': begin
              widget_control,ev.id,get_uvalue=uval
              tlist.method=uval
              widget_control,ev.top,set_uvalue=tlist
              ID=widget_info(ev.top,FIND_BY_UNAME='text')
              widget_control,ID,sensitive=0
              endcase
'OK':          begin
              widget_control,tlist.top,get_uvalue=list
              list.sortseparator=tlist.separator
              list.sortmethod=tlist.method
              widget_control,tlist.top,set_uvalue=list
              widget_control,ev.top,/destroy
              endcase
else:          begin
              tlist.separator=val
              widget_control,ev.top,set_uvalue=tlist
              endcase
endcase

end;pro BrowserOptions_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MakeBOptionsWindow,ev
WIDGET_CONTROL, ev.top, get_uvalue=list
ID1=ev.top
WIDGET_CONTROL, ID1, sensitive=0
tlist={top:ID1,method:list.sortmethod,separator:list.sortseparator}
base=widget_base(/column,title='File Browsing Options',xoffset=200,yoffset=200,$
uvalue=tlist,xsize=300)
base4c=widget_base(base,/column,/exclusive)
    button4=widget_button(base4c,value='Sort by string (file_10.tif,file_2.tif,...)',uvalue=1)
    button5=widget_button(base4c,value='Sort with sep. (file_2.tif,file_10.tif,...)',uvalue=0)
    button6=widget_button(base4c,value='Sort by number (file1_10.tif,file2_2.tif,file2_10.tif,...)',uvalue=2)
label=widget_label(base,value='Separator for sort by number:')
text=widget_text(base,value=list.sortseparator,sensitive=tlist.method eq 0,uname='text',/editable)
OKbutton=widget_button(base,value='OK')
WIDGET_CONTROL, base, /REALIZE
case tlist.method of
0:widget_control,button5,set_button=1
1:widget_control,button4,set_button=1
2:widget_control,button6,set_button=1
endcase

widget_control,OKbutton,/INPUT_FOCUS
Xmanager,'MakeBOptionsWindow',base, event_handler='BrowserOptions_event',$
    cleanup='CleanUp_BrowserOptions',GROUP_LEADER=ev.top
end;pro MakeBOptionsWindow
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadINI,blocktype,blocks=blocks,EDIT=EDIT

inifile=XRDUAconfigDIR()+'XRDUA.ini'

Result = file_search(inifile,count=ct)
if keyword_set(EDIT) then return,inifile
blocks=['$XRD path:',$
'$XRD mskpath:',$
'$XRD pddpath:',$
'$XRD_CHI path:',$
'$XRD_XDI path:',$
'$#Peaks:',$
'$EditMask:',$
'$Scale:',$
'$Show:',$
'$Pointer:',$
'$CW:',$
'$Misc:',$
'$CSlider:',$
'$AutoCor:',$
'$AzInt:',$
'$StatSize:',$
'$OutID:',$
'$ChiSize:',$
'$XRD tiepath:',$
'$CrossUpdate:',$
'$ReadCCDInfo:',$
'$DplotThick:',$
'$RadAx:']
type=(where(blocks eq blocktype))[0]
    
; Defaults
Tags='f'+stringr(indgen(n_elements(blocks)))
CD, CURRENT=dir
datadef=CREATE_STRUCT(Tags,{p:dir,f:'*.tif'},$ ;'$XRD path:'
        {p:dir,f:'*.msk'},$ ;'$XRD mskpath:'
        {p:dir,f:'*.pdd'},$ ;'$XRD pddpath:'
        {p:dir,f:'*.chi'},$ ;'$XRD_CHI path:'
        {p:dir,f:'*.xdi'},$ ;'$XRD_XDI path:'
        {n:500},$ ;'$#Peaks:'
        {n:8,y:800},$ ;'$EditMask:'
        {n1:[0L,80L,1,160L],n2:[1.,0.1,1.]},$ ;'$Scale:'
        {n:[1B,1B,1B,0B,0B,0B,1B,1B]},$ ;'$Show:'
        {n:[0B,5]},$ ;'$Pointer:'
        {n:[5,5]},$ ;'$CW:'
        {n:[1.,0.9,2]},$ ;'$Misc:'
        {n:[1B,0B]},$ ;'$CSlider:'
        {n:0B},$ ;'$AutoCor:'
        {n1:[1,0,0,1,1,0,0,0,0],n2:[10.0,1,!pi/9000,0.01,0.01,0],n3:0L,n4:0b,n5:0b,n6:0,n7:[0b,0b],n8:50.,n9:[1b,1b]},$ ;'$AzInt:'
        {n:500},$ ;'$StatSize:'
        {n1:1L,n2:[-1L,-1L,-2147093248,0L]},$ ;'$OutID:'
        {n:[900,500]},$ ;'$ChiSize:'
        {p:dir,f:'*.spline'},$ ;'$XRD tiepath:'
        {n:1b},$ ;'$CrossUpdate:'
        ReadCCDInfo(),$;'$ReadCCDInfo:'
        {n:1},$;'$DplotThick:'
        {n1:[1,0,0,0,0,10.],n2:'(F0.1)'});'$RadAx:'

; ----Make .ini file when not present----
if ct eq 0 then begin

    ; Make Setup file
    if openw_safe(lun,inifile,width=500) then begin
        print,'INI file not created'
        return,''
    endif
    for i=0l,n_elements(blocks)-1 do begin
       printf,lun,blocks[i]
       n=n_tags(datadef.(i))
       for j=0l,n-1 do printf,lun,(datadef.(i)).(j)
    endfor
    free_lun,lun

    ; Return path and file
    print,'INI file created'
    if type lt 0 then return,inifile $
    else return,datadef.(type)
endif

if type lt 0 then return,inifile

; ----Read ini file when present----
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    print,'INI file corrupted (corrupted field data): return default values.'
    return,datadef.(type)
ENDIF

if openr_safe(lun,inifile) then begin
    print,'INI file not read: return default values.'
    return,datadef.(type)
endif
line=''
len=strlen(blocks[type])
readf,lun,line
while (strmid(line,0,len) ne blocks[type]) and (not eof(lun)) do readf,lun,line
if not eof(lun) then readf,lun,line else begin
    free_lun,lun
    print,'INI file corrupted (missing field): return default values.'
    return,datadef.(type)
endelse

if (strmid(line,0,1) ne '$') then begin
    data=line
    while (strmid(line,0,1) ne '$') and (not eof(lun)) do begin
       readf,lun,line
       data=[data,line]
    endwhile
endif else begin
  free_lun,lun
  print,'INI file corrupted (missing field data): return default values.'
  return,datadef.(type)
endelse

; ----Post-process data block if strings has to be converted----
vlong=0L
case type of
0:    odata={p:data[0],f:data[1]} ;'$XRD path:'
1:    odata={p:data[0],f:data[1]} ;'$XRD mskpath:'
2:    odata={p:data[0],f:data[1]} ;'$XRD pddpath:'
3:    odata={p:data[0],f:data[1]} ;'$XRD_CHI path:'
4:    odata={p:data[0],f:data[1]} ;'$XRD_XDI path:'
5: begin ;'$#Peaks:'
    reads,data[0],vlong
    odata={n:vlong}
    endcase
6:    begin ;'$EditMask:'
    reads,data[0],vlong
    n=vlong
    reads,data[1],vlong
    y=vlong
    odata={n:n,y:y}
    endcase
7:    begin ;'$Scale:'
    tmp1=lonarr(4)
    reads,data[0],tmp1
    tmp2=fltarr(3)
    reads,data[1],tmp2
    odata={n1:tmp1,n2:tmp2}
    endcase
8:    begin ;'$Show:'
    tmp=bytarr(8)
    reads,data,tmp
    odata={n:tmp}
    endcase
9:    begin ;'$Pointer:'
    tmp=intarr(2)
    reads,data,tmp
    odata={n:tmp}
    endcase
10:    begin ;'$CW:'
    tmp=intarr(2)
    reads,data,tmp
    odata={n:tmp}
    endcase
11:    begin ;'$Misc:'
    tmp=fltarr(3)
    reads,data,tmp
    odata={n:tmp}
    endcase
12:    begin ;'$CSlider:'
    tmp=bytarr(2)
    reads,data,tmp
    odata={n:tmp}
    endcase
13:    begin ;'$AutoCor:'
    tmp=0B
    reads,data,tmp
    odata={n:tmp}
    endcase
14:    begin ;'$AzInt:'
    tmp1=intarr(9)
    reads,data[0],tmp1
    tmp2=fltarr(6)
    reads,data[1],tmp2
    tmp3=0L
    reads,data[2],tmp3
    tmp4=0b
    reads,data[3],tmp4
    tmp5=0b
    reads,data[4],tmp5
    tmp6=0
    reads,data[5],tmp6
    tmp7=bytarr(2)
    reads,data[6],tmp7
    tmp8=0.
    reads,data[7],tmp8
    tmp9=bytarr(2)
    reads,data[8],tmp9
    odata={n1:tmp1,n2:tmp2,n3:tmp3,n4:tmp4,n5:tmp5,n6:tmp6,n7:tmp7,n8:tmp8,n9:tmp9}
    endcase
15: begin ;'$StatSize:'
    reads,data[0],vlong
    odata={n:vlong}
    endcase
16: begin ;'$OutID:'
    tmp=0L
    reads,data[0],tmp
    tmp2=lonarr(4)
    if tmp ne 0 then begin
        reads,data[1],tmp2
        tmp3=lonarr(4)
        for l=1,tmp-1 do begin
            reads,data[l+1],tmp3
            tmp2=[[tmp2],[tmp3]]
        endfor
    endif
    odata={n1:tmp,n2:tmp2}
    endcase
17:    begin ;'$ChiSize:'
    tmp=intarr(2)
    reads,data[0],tmp
    odata={n:tmp}
    endcase
18:    odata={p:data[0],f:data[1]} ;'$XRD tiepath:'
19:    begin ;'$CrossUpdate:'
    tmp=0B
    reads,data,tmp
    odata={n:tmp}
    endcase
20:    begin;'$ReadCCDInfo:'
    odata=datadef.(type)
    for i=0l,9 do odata.(i)=data[i]
    endcase
21:    begin;'$DplotThick:'
    tmp=0
    reads,data,tmp
    odata={n:tmp}
    endcase
22:    begin;'$RadAx:'
    tmp1=fltarr(6)
    reads,data[0],tmp1
    tmp2=''
    reads,data[1],tmp2
    odata={n1:tmp1,n2:tmp2}
    endcase
endcase

free_lun,lun
return,odata
end;function ReadINI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UpdateINI,blocktype,data

; ----Create file if necessary----
inifile=ReadINI('',blocks=blocks)
type=(where(blocks eq blocktype))[0]

Result = file_search(inifile,count=ct)

SaveBlock=''
line=''
if openr_safe(lun,inifile) then return,0b
readf,lun,line

; ----Store previous blocks----
len=strlen(blocks[type])
while (strmid(line,0,len) ne blocks[type]) and (not eof(lun)) do begin
    SaveBlock=[SaveBlock,line]
    readf,lun,line
endwhile

; ----Skip block----
if (strmid(line,0,len) eq blocks[type]) then begin ; line='$...:'
    if (not eof(lun)) then begin
        readf,lun,line ; skip line='$...:'
        while (strmid(line,0,1) ne '$') and (not eof(lun)) do readf,lun,line ; skip field data

        ; ----Store next blocks----
        if (strmid(line,0,1) eq '$') then begin ; line='$...:'
            SaveBlock=[SaveBlock,line]
            while not eof(lun) do begin
                readf,lun,line
                SaveBlock=[SaveBlock,line]
            endwhile
        endif
    endif
endif else SaveBlock=[SaveBlock,line] ; end of file and no blocks left

free_lun,lun

; ----Make Setup file----
if openw_safe(lun,inifile,width=500) then return,0b
for i=1,n_elements(SaveBlock)-1 do printf,lun,SaveBlock[i]
printf,lun,blocks[type]
n=n_tags(data)
for i=0l,n-1 do printf,lun,data.(i)
free_lun,lun

return,1b
end;function UpdateINI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SubscribeToOutID

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    tmp=UpdateINI('$OutID:',{n1:0,n2:lonarr(4)})
    return,0
ENDIF

result=ReadINI('$OutID:')
n=result.(0)
ID=result.(1)

ind=where(widget_info(ID[0,*],/VALID_ID),ct)
if ct ne 0 then begin
    widget_control,ID[1,ind[0]],get_uvalue=list
    ind=where(ID[2,*] eq list.randomid,ct)

    OutID=ID[0,ind[0]]
    ID[3,ind[0]]++
endif else begin
    result=makeoutlog()
    OutID=result[0]
    if ID[3,0] eq 0 then begin
        ID=[result,1]
        n=1
    endif else begin
        ID=[[ID],[result,1]]
        n++
    endelse
endelse

tmp=UpdateINI('$OutID:',{n1:n,n2:ID})
return,OutID
end;function SubscribeToOutID
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UnSubscribeToOutID

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    tmp=UpdateINI('$OutID:',{n1:0,n2:lonarr(4)})
    return,0
ENDIF

result=ReadINI('$OutID:')
n=result.(0)
ID=result.(1)

ind=where(widget_info(ID[0,*],/VALID_ID),ct)
if ct eq 0 then return,1b

widget_control,ID[1,ind[0]],get_uvalue=list
ind=where(ID[2,*] eq list.randomid,ct,complement=indcp)
if ct ne 1 then return,1b
ind=ind[0]

if ID[3,ind] le 1 then begin
    widget_control,ID[1,ind],/destroy
    if n gt 1 then ID=ID[*,indcp] $
    else ID[*]=0
    n--
endif else ID[3,ind]--

tmp=UpdateINI('$OutID:',{n1:n,n2:ID})
return,0b
end;function UnSubscribeToOutID
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadPDDefault,lun,struc,PDFname

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    return,0b
ENDIF

PDFn=0
PDFd=0
PDFi=-1

line=''
while (not eof(lun)) do begin
    readf,lun,line
    val=float((str_sep(strtrim(strcompress(line),2),' ')))
    if val[0] ne 0 then begin
        PDFd=[PDFd,val[0]]
        PDFi=[PDFi,val[1]]
        PDFn++
    endif
endwhile
if PDFn eq 0 then return,0b

PDFd=PDFd[1:*]
PDFi=PDFi[1:*]
struc={type:0,PDFname:PDFname,PDFd:PDFd,PDFn:PDFn,PDFi:(PDFi*100./max(PDFi)),$
    PDFni:PDFn,hkl:intarr(3,PDFn)}
return,1b

end;function ReadPDDefault
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadPDDefault2,lun,struc,PDFname

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    return,0b
ENDIF

; Colomn numbers of d, I, h, k, l
line=''
readf,lun,line
ind=long(strsplit(line,' \t',/extract,count=ct))

PDFn=0
PDFd=0
PDFi=-1
hkl=intarr(3)

while (not eof(lun)) do begin
    readf,lun,line
    
    p=strpos(line,',')
    while p ne -1 do begin
        strput,line,'.',p
        p=strpos(line,',',p)
    endwhile
    
    val=float(strsplit(line,' \t',/extract,count=ct2))
    if array_equal(val,fltarr(ct2>1)) then return,0b
    
    if ct ge 1 then PDFd=[PDFd,val[ind[0]]] else PDFd=[PDFd,0]
    if ct ge 2 then PDFi=[PDFi,val[ind[1]]] else PDFi=[PDFi,0]
    hkladd=intarr(3)
    if ct ge 3 then hkladd[0]=val[ind[2]]
    if ct ge 4 then hkladd[1]=val[ind[3]]
    if ct ge 5 then hkladd[2]=val[ind[4]]
    hkl=[[hkl],[hkladd]]
    PDFn++

endwhile
if PDFn eq 0 then return,0b

PDFd=PDFd[1:*]
PDFi=PDFi[1:*]
hkl=hkl[*,1:*]
struc={type:0,PDFname:PDFname,PDFd:PDFd,PDFn:PDFn,PDFi:(PDFi*100./max(PDFi)),$
    PDFni:PDFn,hkl:hkl}
return,1b

end;function ReadPDDefault2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadDS,lun,struc,PDFname

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    return,0b
ENDIF

PDFn=0l
val=0.
while not eof(lun) do begin
    readf,lun,val
    if PDFn eq 0 then PDFd=val else PDFd=[PDFd,val]
    PDFn++
endwhile
if PDFn eq 0 then return,0b

struc={type:0,PDFname:PDFname,PDFd:PDFd,PDFn:PDFn,PDFi:PDFd*0,$
    PDFni:0,hkl:intarr(3,PDFn)}
return,1b
end;function ReadDS
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadAID2,lun,struc

; 82 byte lines:
;    0-70:
;    71: P
;    72-77: ID
;    78-79: XY (X:crystal system, Y:line code)
;    80-81: char13, char10
;  ____________________________________
; |               | |      |  |        |
; |               |P|..ID..|XY|'13''10'|
; |_______________|_|______|__|________|


; line codes Y:
; Y=3: spacegroup | IT number
; Y=6: compound name
; Y=7: chemical formula
; Y=E: unit cell
; Y=F: geometry

point_lun,lun,0
line=bytarr(82)
readu,lun,line
X=line[78]

; ----Name----
; name = 'name; formula; ID'
PDFname='Unknown'
Y=(byte('6'))[0]
repeat readu,lun,line until (((line[78] eq X)and(line[79] eq Y))or(eof(lun)))
if not eof(lun) then begin
    name=strtrim(strcompress(string(line[0:67])),2)
    Y2=(byte('7'))[0]

    readu,lun,line
    while (((line[78] eq X)and((line[79] eq Y)or(line[79] eq Y2)))and(not eof(lun))) do begin
        name=[name,strtrim(strcompress(string(line[0:67])),2)]
        readu,lun,line
    endwhile
    ;name=shift(name,1)
    PDFname=name[0]
    for i=1,n_elements(name)-1 do PDFname+='; '+name[i]
    PDFname+='; ('+strtrim(strcompress(string(line[72:77])),2)+')'
endif

; ----Wavelength----
lambda=0.
Y=(byte('F'))[0]
point_lun,lun,0
repeat readu,lun,line until (((line[78] eq X)and(line[79] eq Y))or(eof(lun)))
if not eof(lun) then begin
    tmp=strsplit(string(line[0:67]),' ',count=ct,/extract)
    f0=stregex(tmp[0],'[0-9.]+',/ext)
    if ct gt 1 then f1=stregex(tmp[1],'[0-9.]+',/ext) else f1=''
    
    if f0 ne '' then lambda=float(f0) else $
    if f1 ne '' then lambda=float(f1)
endif

; ----Reflections----
Y=(byte('I'))[0]
point_lun,lun,0
while (line[79] ne Y) and (not eof(lun)) do readu,lun,line

fl=0.
i=0
h=0
k=0
l=0
PDFd=fl
PDFi=i
hkl=[h,k,l]
c=''

repeat begin
    if eof(lun) then return,0b

    ; 0-19: ref1
    ;    |___________\  0-6: d-spacing
    ;               /  7-9: int
    ; 20-22: space     10:  * or L (i.e. less then)
    ; 23-42: ref2      11-13: h
    ; 43-45: space     14-16: k
    ; 46-65: ref3      17-19: l
    ; 66-70: number
    reads,string(line[0:19]),format='(f7.5,I3,A1,I3,I3,I3)',fl,i,c,h,k,l
    PDFd=[PDFd,fl]
    PDFi=[PDFi,i]
    hkl=[[hkl],[h,k,l]]
    reads,string(line[23:42]),format='(f7.5,I3,A1,I3,I3,I3)',fl,i,c,h,k,l
    PDFd=[PDFd,fl]
    PDFi=[PDFi,i]
    hkl=[[hkl],[h,k,l]]
    reads,string(line[46:65]),format='(f7.5,I3,A1,I3,I3,I3)',fl,i,c,h,k,l
    PDFd=[PDFd,fl]
    PDFi=[PDFi,i]
    hkl=[[hkl],[h,k,l]]

    readu,lun,line
endrep until (line[79] ne Y)

; extract valid entrees
ind=where(PDFd ne 0,PDFn)
if PDFn eq 0 then return,0b

struc={PDFname:PDFname,PDFd:PDFd[ind],PDFn:PDFn,PDFi:(PDFi[ind]*100./max(PDFi[ind])),$
    PDFni:PDFn,hkl:hkl[*,ind],lambda:lambda}

return,1b

end;function ReadAID2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadAID,file

if openr_safe(lun,file) then return,0b

; ----Spacegroup and unit cell----
error=1b
line=bytarr(82)
readu,lun,line
linecel=line
point_lun,lun,0
ITN=0
Y=(byte('3'))[0]
repeat readu,lun,line until (((line[79] eq Y))or(eof(lun)))
if not eof(lun) then begin
    tmp=strsplit(string(line),' ',/extract,count=ct)
    if ct ge 2 then ITN=fix(tmp[1])
endif
if ITN ne 0 then begin
    SG=spacegroup(ITN,2,error=error)
    if ~error then begin
        ; Unit cell parameters
        celparam=fltarr(6)
        reads,strmid(linecel,0,50),celparam,format='(f9,f9,f9,f8,f8,f7)'
        if total(celparam eq 0,/pres) ne 6 then begin
            celparamdef=celparamdefault(SG)
            ifix=celparamfix(SG)
            ind=where(ifix,ct)
            if ct ne 0 then celparam[ind]=celparamdef[ind]
            celparam=celparamsetcon(celparam,SG)
        endif else error=1b
    endif
endif
if error then begin
    free_lun,lun
    return,0b
endif

; ----Atomic positions----
; no available in this format
AddASUpos,ASU

; ----PDF----
bPDF=ReadAID2(lun,struc)

; ----Finish----
free_lun,lun

if bPDF then return,create_struct(struc,'SG',SG,'celparam',celparam,'ASU',ASU) $
else return,create_struct('SG',SG,'celparam',celparam,'ASU',ASU)

end;function ReadAID
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadPDD,lun,struc,OutID

line=''
block='$COMP'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
PDFname='Unknown'
if not eof(lun) then begin
    readf,lun,line
    PDFname=line
    printw,OutID,block+' loaded'
endif

block='$PDFNR'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
if not eof(lun) then begin
    readf,lun,line
    PDFname+='; ('+line+')'
    printw,OutID,block+' loaded'
endif

point_lun,lun,0
block='$D'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
if not eof(lun) then begin
    PDFn=0
    readf,lun,PDFn
    PDFd=fltarr(PDFn)
    readf,lun,PDFd
    printw,OutID,block+' loaded'
endif else return,0b

point_lun,lun,0
block='$I'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
PDFi=fltarr(PDFn)
PDFni=0
if not eof(lun) then begin
    PDFni=PDFn
    readf,lun,PDFi
    PDFi=PDFi/max(PDFi)*100
    printw,OutID,block+' loaded'
endif

point_lun,lun,0
block='$hkl'
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
hkl=intarr(3,PDFn)
if not eof(lun) then begin
    readf,lun,hkl
    printw,OutID,block+' loaded'
endif

struc={type:0,PDFname:PDFname,PDFd:PDFd,PDFn:PDFn,PDFi:PDFi,PDFni:PDFni,hkl:hkl}
return,1b

end;function ReadPDD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadBEX,lun,struc,OutID,rever=rev,revskip=revskip,nskip=nskip,revoff=revoff

; rever: get card before lun pointer
; revskip: set here to number of cards in the BEX file
; nskip: number of cards to skip
; revoff: skip cards in reverse reading

if keyword_set(rev) then begin ; point lun to previous 'TI' line
    if n_elements(revoff) eq 0 then revoff=1

    ; count number of cards
    if eof(lun) then begin
        point_lun,lun,0
        revskip=0

        line=''
        while not eof(lun) do begin
            readf,lun,line
            if strmid(line,0,2) eq 'TI' then revskip++
        endwhile
    endif

    ; point to previous card
    point_lun,-lun,pos
    if pos ne 0 then begin
        barr=bytarr(512)
        for i=0l,revoff do begin
            p=-1
            while (p eq -1) and (pos gt 0) do begin
                pos-=512
                if pos lt 0 then begin
                    point_lun,lun,0
                    barr=barr[0:512+pos-1]
                    readu,lun,barr
                    pos=0
                endif else begin
                    point_lun,lun,pos
                    readu,lun,barr
                endelse

                p=STRPOS( string(barr), 'TI', /REVERSE_SEARCH)
            endwhile
        endfor
        point_lun,lun,(pos+p)>0
    endif

    nskip=0
endif

; First line
if n_elements(nskip) eq 0 then nskip=0 else nskip>=0
for i=0l,nskip do begin
    line=''
    while (strmid(line,0,2) ne 'TI') and (not eof(lun)) do readf,lun,line
    if eof(lun) then return,0
endfor

temp=str_sep(STRTRIM(strcompress(line),2),' ')
if n_elements(temp) gt 3 then PDFname=STRJOIN(temp[3:n_elements(temp)-1],' ')+'; ' else PDFname='' ; Name
PDFname+=temp[2] ; Formula
PDFname+='; ('+temp[1]+')' ; PDF number

; Second line
readf,lun,line
if strmid(line,0,2) ne 'WL' then begin
    printw,OutID,'Unknown BEX format.'
    return,0
endif

; Data block
temp=fltarr(6)
readf,lun,line
reads,line,temp
PDFd=temp[2]
PDFn=1
PDFi=temp[4]
PDFni=1
hkl=intarr(3)

readf,lun,line
while strlen(line) ne 0 do begin
    reads,line,temp
    PDFd=[PDFd,temp[2]]
    PDFn++
    PDFi=[PDFi,temp[4]]
    PDFni++
    hkl=[[hkl],[intarr(3)]]
    if eof(lun) then line='' else readf,lun,line
endwhile

printw,OutID,'...'+PDFname
struc={type:0,PDFname:PDFname,PDFd:PDFd,PDFn:PDFn,PDFi:PDFi,PDFni:PDFni,hkl:hkl}
return,1b
end;function ReadBEX
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadCif2,lun,struc

; ----Name----
; name = 'name; formula'
point_lun,lun,0
PDFname=''
blocks='_chemical_name'+['_mineral','_common','_systematic','']
i=0
repeat begin
    line=''
    point_lun,lun,0
    block=blocks[i]
    nblock=strlen(block)
    while strmid(line,0,nblock) ne block and not eof(lun) do readf,lun,line
    if not eof(lun) then begin
        line=strsplit(line,"'",/extract,count=ct)
        PDFname=line[ct eq 2]
    endif else if i eq 3 then PDFname=' '
    i++
endrep until PDFname ne ''

point_lun,lun,0
line=''
block='_chemical_formula_sum'
nblock=strlen(block)
while strmid(line,0,nblock) ne block and not eof(lun) do readf,lun,line
if not eof(lun) then begin
    line=strsplit(line,"'",/extract,count=ct)
    PDFname+='; '+line[ct eq 2]
endif

; ----Wavelength----
lambda=0.
point_lun,lun,0
line=''
block='wavelength'
while strpos(line,block) eq -1 and not eof(lun) do readf,lun,line
if not eof(lun) then begin
    line=strsplit(line," ",/extract,count=ct)
    if ct eq 2 then lambda=float(line[1])
endif

; ----Reflections----
;loop_
;_pd_peak_d_spacing
;_pd_peak_calc_intensity_net
;     5.415646        108.19
point_lun,lun,0
line=''
block='_pd_peak_d_spacing'
nblock=strlen(block)
while strmid(line,0,nblock) ne block and not eof(lun) do readf,lun,line
while strmid(line,0,1) eq '_' and not eof(lun) do readf,lun,line
if eof(lun) then return,0b

data=fltarr(2)
repeat begin
    tmp=strsplit(line," ",/extract,count=ct)
    if ct ge 2 then data=[[data],[float(tmp[0:1])]]
    
    bool=ct lt 2 or eof(lun)
    if not eof(lun) then readf,lun,line
endrep until bool
if n_elements(data) eq 2 then return,0b

PDFd=reform(data[0,1:*])
n=n_elements(PDFd)
PDFi=reform(data[1,1:*])
struc={type:0,PDFname:PDFname,PDFd:PDFd,PDFn:n,PDFi:PDFi/max(PDFi)*100,PDFni:n,hkl:intarr(3,n),lambda:lambda}

return,1b

end;function ReadCif2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadCif,file,SmallDev,nopos=nopos

if openr_safe(lun,file) then return,0b
line=''

; ----Space group----
;_symmetry_space_group_name_Hall      'R 3 -2"'
;_symmetry_space_group_name_H-M 'F m -3 m'
;_space_group_IT_number 225
error=1b
block="_symmetry_space_group_name_Hall"
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do begin
    readf,lun,line
    line=strtrim(line,2)
endwhile
if not eof(lun) then begin
    tmp=strsplit(line,block,/extract,/regex)
    tmp=strtrim(tmp[0],2)
    tmp=strmid(tmp,1,strlen(tmp)-2)
    SG=spacegroup(tmp,0,error=error,Pconv=Pconv)
endif
if error then begin

    point_lun,lun,0
    line=''
    block="_symmetry_space_group_name_H-M"
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do begin
        readf,lun,line
        line=strtrim(line,2)
    endwhile
    if not eof(lun) then begin
        tmp=strsplit(line,block,/extract,/regex)
        tmp=strtrim(tmp[0],2)
        tmp=strmid(tmp,1,strlen(tmp)-2)
        SG=spacegroup(tmp,1,error=error,Pconv=Pconv)
    endif
    
    if error then begin
        point_lun,lun,0
        line=''
        block="_space_group_IT_number"
        len=strlen(block)
        while (strmid(line,0,len) ne block) and (not eof(lun)) do begin
            readf,lun,line
            line=strtrim(line,2)
        endwhile
        if not eof(lun) then begin
            tmp=strsplit(line," ",/extract,count=ct)
            SG=spacegroup(fix(tmp[ct-1]),2,error=error,Pconv=Pconv)
        endif
    endif
endif

; Get space group from symmetry operations:
bsymm=0b
ops=''
point_lun,lun,0
line=''
block0='_symmetry_equiv_pos_site_id'
len0=strlen(block0)
bskipcolumn=0b
block="_symmetry_equiv_pos_as_xyz"
len=strlen(block)
while (strmid(line,0,len) ne block) and (not eof(lun)) do begin
    readf,lun,line
    line=strtrim(line,2)
    
    bskipcolumn or= strmid(line,0,len0) eq block0
endwhile
bsymm=not eof(lun)

if eof(lun) then begin
    point_lun,lun,0
    line=''
    block0='_space_group_symop_id'
    len0=strlen(block0)
    bskipcolumn=0b
    
    block="_space_group_symop_operation_xyz"
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do begin
        readf,lun,line
        line=strtrim(line,2)
        bskipcolumn or= strmid(line,0,len0) eq block0
    endwhile
    bsymm=not eof(lun)
endif

if bsymm then begin
    error or= DIALOG_MESSAGE('Use symmetry operations instead of S.G. symbol?',/question) eq 'Yes'
    
    if error then begin
        repeat begin
            readf,lun,line
            line=strtrim(line,2)
            
            if bskipcolumn then begin
                p=strpos(line,' ')
                if p ne -1 then line=strmid(line,p+1)
            endif
                
            p1=strpos(line,',')
            p2=strpos(line,',',/reverse_search)
            bstop=p1 eq -1 or p2 eq -1 or p1 eq p2
            if ~bstop then ops+=(strlen(ops) eq 0)?line:(";"+line)
        endrep until bstop
        
        SG=spacegroup(ops,3,error=error,Pconv=Pconv)
    endif
endif

if error then begin
    free_lun,lun
    return,0b
endif

; ----Cell param----
;_cell_length_a                       12.6963024(3)
;_cell_length_b                       12.6963024(3)
;_cell_length_c                       7.9210911(7)
;_cell_angle_alpha                    90.00000
;_cell_angle_beta                     90.00000
;_cell_angle_gamma                    120.0(3)
celparam=fltarr(6)
point_lun,lun,0
line=''
block='_cell_length'
blockcommon='_cell_'
len=strlen(block)
lencommon=strlen(blockcommon)
while (strmid(line,0,len) ne block) and (not eof(lun)) do begin
    readf,lun,line
    line=strtrim(line,2)
endwhile
if not eof(lun) then begin
    param=blockcommon+['length_a','length_b','length_c','angle_alpha','angle_beta','angle_gamma']
    repeat begin
        tmp=strsplit(line," ",/extract,count=ct)
        j=where(param eq tmp[0],ct2)
        if ct2 eq 1 then begin
            tmp=strsplit(tmp[ct-1],"()",/extract)
            celparam[j]=float(tmp[0])
        endif
        readf,lun,line
    endrep until (strmid(line,0,lencommon) ne blockcommon)
endif else begin
    free_lun,lun
    return,0b
endelse
TransformCelparam,celparam,Pconv
Pconv=invert(Pconv) ; for atomic coordinates

; ----Atomic positions----
;loop_
;    _atom_site_label
;    _atom_site_fract_x
;    _atom_site_fract_y
;    _atom_site_fract_z
;    _atom_site_U_iso_or_equiv
;    _atom_site_occupancy
;    _atom_site_adp_type              # Not in version 2.0.1
;    _atom_site_type_symbol
; al6    0.33760     0.29120     0.34450     0.01064     1.00000    Uiso al
; cr1    0.00000     0.00000     0.17950     0.00215     0.40000    Uiso cr
if ~keyword_set(nopos) then begin
    point_lun,lun,0
    
    i_name=-1
    while i_name eq -1 and not eof(lun) do begin
        line=''
        block='_atom_site_'
        len=strlen(block)
        while (strmid(line,0,len) ne block) and (not eof(lun)) do begin
            readf,lun,line
            line=strtrim(line,2)
        endwhile
        if not eof(lun) then begin
            ; Default table indices
            i_name=-1
            i_ZIon=-1
            i_x=-1
            i_y=-1
            i_z=-1
            i_SOF=-1
            i_Biso=-1
            name=''
            ZIon=[1,1]
            xyz=fltarr(1,4)
            xyz[3]=1
            SOF=1.
            Biso=1.
            isomult=1. ; Biso=Uiso*8*!pi*!pi
            
            ; CIF table indices
            ncol=0
            repeat begin
                case (strsplit(line,block,/extract,/regex))[0] of
                'label':        i_name=ncol
                'type_symbol':    i_ZIon=ncol
                'fract_x':        i_x=ncol
                'fract_y':        i_y=ncol
                'fract_z':        i_z=ncol
                'occupancy':    i_SOF=ncol
                'B_iso_or_equiv': i_Biso=ncol
                'Uiso_or_equiv': begin
                                i_Biso=ncol
                                isomult=8*!pi*!pi
                                endcase
                'U_iso_or_equiv': begin
                                i_Biso=ncol
                                isomult=8*!pi*!pi
                                endcase
                else:
                endcase
                ncol++
                readf,lun,line
                line=strtrim(line,2)
            endrep until strmid(line,0,len) ne block
            line=strsplit(line,' ',/extract,count=ct)
            if i_ZIon eq -1 then i_ZIon=i_name
        endif
    endwhile
    
    if not eof(lun) then begin
        ; Read atom table
        while ct eq ncol do begin
            if i_name ne -1 then name=line[i_name]
            if i_ZIon ne -1 then Z=ZElement(line[i_ZIon],Ion=Ion,repeats=repeats)
            if i_x ne -1 then xyz[0]=float(line[i_x])
            if i_y ne -1 then xyz[1]=float(line[i_y])
            if i_z ne -1 then xyz[2]=float(line[i_z])
            if i_SOF ne -1 then SOF=float(line[i_SOF])
            if i_Biso ne -1 then Biso=float(line[i_Biso])*isomult
            
            DefaultIon,Z,Ion,repeats
            
            xyz=Pconv##xyz
            AddASUpos,ASU,name,Z,Ion,repeats,reform(xyz[0:2]),SOF,Biso
            
            line=''
            if eof(lun) then ct=ncol+1 else begin
                readf,lun,line
                line=strsplit(line,' ',/extract,count=ct)
            endelse
        endwhile
        
        ; Anisotropic factors
        point_lun,lun,0
        line=''
        block='_atom_site_aniso_'
        len=strlen(block)
        while (strmid(line,0,len) ne block) and (not eof(lun)) do begin
            readf,lun,line
            line=strtrim(line,2)
        endwhile
        if not eof(lun) then begin
            ; Default table indices
            i_name=-1
            i_B11=-1
            i_B22=-1
            i_B33=-1
            i_B12=-1
            i_B13=-1
            i_B23=-1
            isomult=1. ; B=U*8*!pi*!pi
            
            ; CIF table indices
            ncol=0
            repeat begin
                case (strsplit(line,block,/extract,/regex))[0] of
                'label':i_name=ncol
                'B_11': i_B11=ncol
                'B_22': i_B22=ncol
                'B_33': i_B33=ncol
                'B_12': i_B12=ncol
                'B_13': i_B13=ncol
                'B_23': i_B23=ncol
                'U_11': begin
                        i_B11=ncol
                        isomult=8*!pi*!pi
                        endcase
                'U_22': begin
                        i_B22=ncol
                        isomult=8*!pi*!pi
                        endcase
                'U_33': begin
                        i_B33=ncol
                        isomult=8*!pi*!pi
                        endcase
                'U_12': begin
                        i_B12=ncol
                        isomult=8*!pi*!pi
                        endcase
                'U_13': begin
                        i_B13=ncol
                        isomult=8*!pi*!pi
                        endcase
                'U_23': begin
                        i_B23=ncol
                        isomult=8*!pi*!pi
                        endcase
                else:
                endcase
                ncol++
                readf,lun,line
                line=strtrim(line,2)
            endrep until strmid(line,0,len) ne block
            line=strsplit(line,' ',/extract,count=ct)
            
            ; Read atom table
            if i_B11 ne -1 and i_B22 ne -1 and i_B33 ne -1 and $
                   i_B12 ne -1 and i_B13 ne -1 and i_B23 ne -1 and i_name ne -1 then begin
                   
                celparamr = DirectToReciprocalU(celparam)
                celparamr[3:5] = 0
;                tmp = celparam & tmp[3:5] = 0 
                ; The above should be done according to
                ;     http://www.iucr.org/__data/iucr/cifdic_html/1/cif_core.dic/Iatom_site_B_iso_or_equiv.html
                ; but not according to 
                ;    http://scripts.iucr.org/cgi-bin/paper?es0238
                isomult *= MetricTensor(celparam)*MetricTensor(celparamr)
                while ct eq ncol do begin
                    name=line[i_name]
                    Baniso=float(line[[i_B11,i_B22,i_B33,i_B12,i_B13,i_B23]])
                    Biso=total(isomult*Baniso[[[0,3,4],[3,1,5],[4,5,2]]])/3
                    
                    ; Edit ASU
                    tmp=where(ASU.asupos.name eq name)
                    if tmp[0] ne -1 then ASU.asupos[tmp[0]].B=Biso
                    
                    line=''
                    if eof(lun) then ct=ncol+1 else begin
                        readf,lun,line
                        line=strsplit(line,' ',/extract,count=ct)
                    endelse
                endwhile
            endif
        endif
        
    endif else begin
        free_lun,lun
        return,0b
    endelse
endif
AddASUpos,ASU

; ----PDF----
bPDF=ReadCif2(lun,struc)

; ----Finish----
free_lun,lun
ResetWyckASU,SG,ASU,SmallDev

if bPDF then return,create_struct(struc,'ASU',ASU,'SG',SG,'celparam',celparam) $
else return,create_struct('ASU',ASU,$ ; asymmetric unit
                'SG',SG,$    ; space group information
                'celparam',celparam); 6 celparameters: A,A,A,deg,deg,deg

end;function ReadCif
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StorePDD,list,struc,badd,PDFSearch

case PDFSearch.type of
0:    begin
    tmp=strpos(strlowcase(struc.PDFname),PDFSearch.PDFSubstr)
    if tmp[0] eq -1 and PDFSearch.PDFSubstr ne '' then return,0b
    endcase
1:    begin
    tmp=where((struc.PDFd ge PDFSearch.PDFdi) and (struc.PDFd le PDFSearch.PDFde),ct)
    if ct eq 0 then return,0b
    endcase
endcase

; sort
ind=reverse(sort(struc.PDFd))
struc.PDFd=(struc.PDFd)[ind]
struc.PDFi=(struc.PDFi)[ind]
struc.hkl=(struc.hkl)[*,ind]

; Inherit "Use intensities" setting from the majority of the loaded PDFs
;if ptr_valid(list.PDFni) then begin
;    temp=where((*list.PDFni) eq 0,ct1)
;    temp=where((*list.PDFni) ne 0,ct2)
;    if ct1 ge ct2 then struc.PDFni=0
;endif

if (list.listtype eq 11) or (list.listtype eq 31) then begin
    if (list.nPDF gt 0) then PDFScale=(*list.PDFScale)[list.nPDF-1] $
    else PDFScale=-0.01
    PDFOffset=0.
endif

if (list.listtype eq 11) or (list.listtype eq 31) then nmodes=n_elements(list.xmodes)
if (badd ne 0) and (list.nPDF gt 0) then begin ; Add new phase or replace last one
    if (badd eq 2) then begin
        ; Replace last one
        if (list.nPDF eq 1) then begin
            (*list.PDFname)=struc.PDFname
            (*list.PDFn)=struc.PDFn
            list.PDFnTotal=struc.PDFn
            (*list.PDFd)=struc.PDFd
            (*list.PDFhkl)=struc.hkl
            (*list.PDFni)=struc.PDFni
            (*list.PDFi)=struc.PDFi
            (*list.PDFs)=replicate(1b,struc.PDFn)

            if (list.listtype eq 11) or (list.listtype eq 31) then begin
                (*list.PDFx)=fltarr(struc.PDFn,nmodes)
                (*list.PDFi2)=struc.PDFi
                (*list.PDFScale)=PDFScale
                (*list.PDFOffset)=PDFOffset
            endif
        endif else begin
            tmp=where((*list.PDFname)[0:list.nPDF-2] eq struc.PDFname,ct)
            if ct ne 0 then begin
                printw,list.OutID,struc.PDFname+' is already loaded'
                return,0b
            endif
        
            list.PDFnTotal-=(*list.PDFn)[list.nPDF-1]
            (*list.PDFname)=[(*list.PDFname)[0:list.nPDF-2],struc.PDFname]
            (*list.PDFn)=[(*list.PDFn)[0:list.nPDF-2],struc.PDFn]
            (*list.PDFd)=[(*list.PDFd)[0:list.PDFnTotal-1],struc.PDFd]
            (*list.PDFhkl)=[[(*list.PDFhkl)[*,0:list.PDFnTotal-1]],[struc.hkl]]
            (*list.PDFni)=[(*list.PDFni)[0:list.nPDF-2],struc.PDFni]
            (*list.PDFi)=[(*list.PDFi)[0:list.PDFnTotal-1],struc.PDFi]
            (*list.PDFs)=[(*list.PDFs)[0:list.PDFnTotal-1],replicate(1b,struc.PDFn)]

            if (list.listtype eq 11) or (list.listtype eq 31) then begin
                (*list.PDFx)=[(*list.PDFx)[0:list.PDFnTotal-1,*],fltarr(struc.PDFn,nmodes)]
                (*list.PDFi2)=[(*list.PDFi2)[0:list.PDFnTotal-1],struc.PDFi]
                (*list.PDFScale)=[(*list.PDFScale)[0:list.nPDF-2],PDFScale]
                (*list.PDFOffset)=[(*list.PDFOffset)[0:list.nPDF-2],PDFOffset]
            endif
            list.PDFnTotal+=struc.PDFn
        endelse
    endif else begin
        ; Add new
        tmp=where((*list.PDFname)[0:list.nPDF-1] eq struc.PDFname,ct)
        if ct ne 0 then begin
            printw,list.OutID,struc.PDFname+' is already loaded'
            return,0b
        endif
    
        (*list.PDFname)=[(*list.PDFname),struc.PDFname]
        (*list.PDFn)=[(*list.PDFn),struc.PDFn]
        list.PDFnTotal+=struc.PDFn
        (*list.PDFd)=[(*list.PDFd),struc.PDFd]
        (*list.PDFhkl)=[[(*list.PDFhkl)],[struc.hkl]]
        (*list.PDFni)=[(*list.PDFni),struc.PDFni]
        (*list.PDFi)=[(*list.PDFi),struc.PDFi]
        list.nPDF++
        (*list.PDFs)=[(*list.PDFs),replicate(1b,struc.PDFn)]

        if (list.listtype eq 11) or (list.listtype eq 31) then begin
            (*list.PDFx)=[(*list.PDFx),fltarr(struc.PDFn,nmodes)]
            (*list.PDFi2)=[(*list.PDFi2),struc.PDFi]
            (*list.PDFScale)=[(*list.PDFScale),PDFScale]
            (*list.PDFOffset)=[(*list.PDFOffset),PDFOffset]
        endif
    endelse
endif else begin ; Delete all phases
    ptr_free,list.PDFname
    ptr_free,list.PDFn
    ptr_free,list.PDFd
    ptr_free,list.PDFhkl
    ptr_free,list.PDFni
    ptr_free,list.PDFi
    ptr_free,list.PDFs

    list.PDFnTotal=struc.PDFn
    list.nPDF=1
    list.PDFname=ptr_new(struc.PDFname)
    list.PDFn=ptr_new(struc.PDFn)
    list.PDFd=ptr_new(struc.PDFd)
    list.PDFhkl=ptr_new(struc.hkl)
    list.PDFni=ptr_new(struc.PDFni)
    list.PDFi=ptr_new(struc.PDFi)
    list.PDFs=ptr_new(replicate(1b,struc.PDFn))

    if (list.listtype eq 11) or (list.listtype eq 31) then begin
        ptr_free,list.PDFx
        ptr_free,list.PDFi2
        ptr_free,list.PDFScale
        ptr_free,list.PDFOffset
        list.PDFx=ptr_new(fltarr(struc.PDFn,nmodes))
        list.PDFi2=ptr_new(struc.PDFi)
        list.PDFScale=ptr_new([PDFScale])
        list.PDFOffset=ptr_new([PDFOffset])
    endif
endelse

return,1b
end;function StorePDD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FindNextPDD,list,rev=rev,i0=i0

krev=keyword_set(rev)
if krev then bnext=(list.BEXpos le 0) else bnext=list.BEXeof

if bnext then begin
    pathi=file_search(list.pddpath+['*.pdd','*.aid','*.ds','*.bex','*.txt','*.dat','*.cif'],count=nfiles,/FOLD_CASE)
    if nfiles eq 0 then return,0
    path = SortListFiles(pathi,0,list.sortseparator)
    ind = where(path eq list.pddpath+list.pddfile,count)
    ind=ind[0]
    b=n_elements(i0) eq 0
    if b then i0=ind

    if krev then begin
        if ind[0] le 0 then ind=nfiles-1 else ind--
    endif else ind++

    ind mod= nfiles
    if ~b then i0=ind
    path=path[ind]
    error=CutPath(path,file=file,ext=ext)
    list.pddfile=file+ext

    list.BEXpos=-1
    list.BEXeof=1b
endif else begin
    i0=-1
endelse

return,bnext

end;function FindNextPDD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LoadPDD,ev,prompt,searh=search

widget_control,ev.top,get_uvalue=list
if n_elements(list) eq 0 then return
bsubfilereading=(list.listtype eq 10) or (list.listtype eq 11)

; ----Different loading options----
; search: 1 or 2 => jump to next file when no substring comp.
PDFSearch=list.PDFSearch
case prompt of ; PDF >> or <<
0:    begin
    ; files
    nfiles=1
    path=list.pddpath+list.pddfile

    ; subsearch
    bsearch=keyword_set(search)
    if bsearch then searchdir=search-1 else searchdir=0 ; 0: forward, 1: backwards
    PDFSearch.PDFSubstr=strlowcase(PDFSearch.PDFSubstr)

    ; add/replace
    badd=list.ReplaceLastPDF+1b
    list.ReplaceLastPDF=1
    endcase
1:    begin ; load PDF
    ; files
    path=''
    tmp=CutPath(list.pddfile,ext=ext)
    filter=['*.pdd','*.aid','*.ds','*.bex','*.txt','*.dat','*.cif']
    sortformat,filter,ext
    path=CAPDIALOG_PICKFILE(path=list.pddpath,file=list.pddfile,filter=filter,title='Load PDF data...')
    if path[0] eq '' then return
    nfiles=1

    ; subsearch
    PDFSearch.PDFSubstr=''
    PDFSearch.type=0
    bsearch=0b
    searchdir=0

    ; add/replace
    Result = DIALOG_MESSAGE('Delete current PDF data?',/question,/DEFAULT_NO)
    badd= result eq 'No'

    ; BEX sub
    if bsubfilereading then begin
        list.BEXpos=-1
        list.BEXeof=1b
    endif
    endcase
2:    begin ; load multiple PDF
    tmp=CutPath(list.pddfile,ext=ext)
    filter=['.pdd','.aid','.ds','.bex','.txt','.dat','.cif']
    sortformat,filter,ext
    
    ; files
    list2={path:list.pddpath,file:list.pddfile,Platform:list.Platform,Format:filter,OutID:list.OutID}
    path=Select_Files(ev.top,list2,sortmethod=list.sortmethod,separator=list.sortseparator)
    if path[0] eq '' then return
    nfiles=n_elements(path)
    Result = DIALOG_MESSAGE(stringr(nfiles)+' files to add: proceed?',/question)
    if result eq 'No' then return

    ; subsearch
    substr=list.PDFSearch.PDFSubstr
    if substr eq '' then substr=' '
    substr=PromptNumber(substr,ev.top,'Add only PDF containing:')
    substr=strlowcase(substr)
    if substr eq ' ' then substr=''
    PDFSearch.PDFSubstr=substr
    PDFSearch.type=0
    bsearch=0b
    searchdir=0

    ; add/replace
    badd= 1b

    ; BEX sub
    if bsubfilereading then begin
        list.BEXpos=-1
        list.BEXeof=1b
    endif
    endcase
endcase
; badd:
; 0: Delete all phases
; 1: Add new phase
; 2: Replace last phase

; ----Read File(s)----
repeat begin ; LOOP: search all files until substring is found
    for i=0l,nfiles-1 do begin ; LOOP: all files to read

        ; ----Open file----
        if openr_safe(lun, path[i]) then begin
            printw,list.OutID,path[i] + ' not found'
            return
        endif
        error=CutPath(path[i],path=pddpath,file=file,ext=ext)
        list.pddpath=pddpath
        list.pddfile=file+ext
        ext=strlowcase(ext)

        CATCH, Error_status
        IF Error_status eq 0 THEN BEGIN
              ; ----Read and store----
              ; Make sure to reset BEX flags
              case ext of
              '.pdd': begin
                      result=ReadPDD(lun,struc,list.OutID)
                      if result ne 0 then bstored=StorePDD(list,struc,badd,PDFSearch) else bstored=0b
                      endcase
              '.aid': begin
                      result=ReadAID2(lun,struc)
                      if result ne 0 then bstored=StorePDD(list,struc,badd,PDFSearch) else bstored=0b
                      endcase
              '.cif':    begin
                      result=ReadCif2(lun,struc)
                      if result ne 0 then bstored=StorePDD(list,struc,badd,PDFSearch) else bstored=0b
                      endcase
              '.txt': begin
                      result=ReadPDDefault(lun,struc,list.pddfile)
                      if result ne 0 then bstored=StorePDD(list,struc,badd,PDFSearch) else bstored=0b
                      endcase
               '.dat': begin
                      result=ReadPDDefault2(lun,struc,list.pddfile)
                      if result ne 0 then bstored=StorePDD(list,struc,badd,PDFSearch) else bstored=0b
                      endcase
              '.ds':    begin
                      result=ReadDS(lun,struc,list.pddfile)
                      if result ne 0 then bstored=StorePDD(list,struc,badd,PDFSearch) else bstored=0b
                      endcase
              '.bex': begin
                      ; BEX is an multi-card format, therefore do a local substring search if asked for
                      if prompt eq 0 then begin ; read 1 card in BEX
      
                          ; read and check 1 card
                          if searchdir then begin
                              BEXpos0=list.BEXpos eq -1
      
                              if BEXpos0 then begin
                                  tmp = FSTAT(lun)
                                  POINT_LUN, lun, tmp.size
      
                                  ; skip and read
                                  result=ReadBEX(lun,struc,list.OutID,/rever,revoff=0,revskip=revskip)
      
                                  list.BEXpos=revskip-1
                                  list.BEXeof=1b
                              endif else begin
      
                                  ; skip and read
                                  list.BEXpos--
                                  result=ReadBEX(lun,struc,list.OutID,nskip=list.BEXpos)
                                  list.BEXeof=0b
                              endelse
      
                          endif else begin
      
                              ; skip and read
                              list.BEXpos++
                              result=ReadBEX(lun,struc,list.OutID,nskip=list.BEXpos)
                              list.BEXeof=eof(lun)
                          endelse
      
                          ; store result
                          if result ne 0 then bstored=StorePDD(list,struc,badd,PDFSearch) else bstored=0b
      
                          ; search next cards
                          if bsearch then begin
                              bloop = (~bstored) and (searchdir?(list.BEXpos gt 0):(~list.BEXeof))
      
                              while bloop do begin
                                  result=ReadBEX(lun,struc,list.OutID,rever=searchdir)
                                  list.BEXpos+=(searchdir?(-1):1)
                                  list.BEXeof=eof(lun)
                                  if result ne 0 then bstored=StorePDD(list,struc,badd,PDFSearch) else bstored=0b
                                  bloop = (~bstored) and (searchdir?(list.BEXpos gt 0):(~list.BEXeof))
                              endwhile
                          endif
      
                      endif else begin ; read all cards in BEX
                          while not eof(lun) do begin
                              result=ReadBEX(lun,struc,list.OutID,rever=searchdir)
                              if result ne 0 then bstored=StorePDD(list,struc,badd,PDFSearch) else bstored=0b
                          endwhile
                      endelse
                      endcase
              else:
              endcase
        ENDIF ELSE BEGIN
          result=0b
          bstored=0b
        ENDELSE
    
        ; ----Close file----
        free_lun,lun
        printw,list.OutID,'Tried loading: '+path[i]
    endfor

    ; ----Jump to next file when substring not found----
    if bsearch then begin
        if (not bstored) then begin
            b=FindNextPDD(list,rev=searchdir,i0=i0)
            if n_elements(isearch0) eq 0 then isearch0=i0 else begin
                if isearch0 eq i0 then begin
                    bstored=1b
                    case PDFSearch.type of
                    0:    printw,list.OutID,'Nothing matched to "'+PDFSearch.PDFSubstr+'"'
                    1:    printw,list.OutID,'Nothing matched to ['+string(PDFSearch.PDFdi)+' -'+string(PDFSearch.PDFde)+'] '+msymbols('angstroms')+'"'
                    endcase
                endif
            endelse
            path=[list.pddpath+list.pddfile]
        endif
    endif else bstored=1b
endrep until bstored

widget_control,ev.top,set_uvalue=list

end;pro LoadPDD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SaveCif,file,ptr,O=O,peaks=peaks,peaknr=peaknr

if openw_safe(lun,file) then return,0b

CALDAT, SYSTIME(/JULIAN), Month , Day , Year 
vandaag=string(Day,Month,Year,format='(I02,"/",I02,"/",I04)')
vandaag2=string(Day,Month,Year,format='(I02,"_",I02,"_",I04)')

printf,lun,"##############################################################################"
printf,lun,"###                 Entry exported on "+vandaag+" from XRDUA                ###"
printf,lun,"###                          http://xrdua.ua.ac.be/                        ###"
printf,lun,"##############################################################################"

; Data block
printf,lun,"data_xrdua_"+vandaag2

; Name
printf,lun,"_chemical_name_common  '"+(*ptr).name+"'"

; Symmetry
sg=sgdata((*(*ptr).pstrinfo).sghash,0)
lg=pgdata((*(*ptr).pstrinfo).lghash,0)

printf,lun,"_symmetry_space_group_name_Hall '"+sg.hall+"'"
printf,lun,"_symmetry_space_group_name_H-M '"+sg.hm+"'"
printf,lun,"_space_group_IT_number ",sg.num
printf,lun,"_space_group_crystal_system '"+strlowcase(lg.crystsys)+"'"
printf,lun,"loop_"
printf,lun,"_symmetry_equiv_pos_as_xyz"
printf,lun,"  '"+strcompress(StringRTopcode((*(*ptr).pstrinfo).allops),/remove_all)+"'",format='(A)'

; Unit cell
if keyword_set(O) then p=(*(*ptr).globalparamIO.O)[0] $
else p=(*ptr).globalparamIO.I

if ~ptr_valid(p) then begin
    free_lun,lun
    return,0b
endif
celparam=(*p)[indCelParam(ptr)]
printf,lun,"_cell_length_a       ",celparam[0]
printf,lun,"_cell_length_b       ",celparam[1]
printf,lun,"_cell_length_c       ",celparam[2]
printf,lun,"_cell_angle_alpha    ",celparam[3]
printf,lun,"_cell_angle_beta     ",celparam[4]
printf,lun,"_cell_angle_gamma    ",celparam[5]
printf,lun,"_cell_volume         ",UCVolume(celparam)

; ASU
ASU=CopystrucASU(ptr,O=O)
if ASU.nasupos gt 0 then begin
    printf,lun,"loop_"
    printf,lun,"_atom_site_label"
    printf,lun,"_atom_site_fract_x"
    printf,lun,"_atom_site_fract_y"
    printf,lun,"_atom_site_fract_z"
    printf,lun,"_atom_site_occupancy"
    printf,lun,"_atom_site_B_iso_or_equiv"
endif
for i=0l,ASU.nasupos-1 do begin
    k0=ASU.reverse_indices[i]
    k1=ASU.reverse_indices[i+1]-1
    for k=k0,k1 do begin
        str=ZElement(*ASU.asupos[k].Z,Ion=*ASU.asupos[k].Ion,Repeats=*ASU.asupos[k].Repeats)
        printf,lun,str,ASU.asupos[k].xyz,ASU.asupos[k].SOF,ASU.asupos[k].B,format='(A,3F,F,F)'
    endfor
endfor

; Peaks
if peaknr gt 0 then begin
    printf,lun,"loop_"
    printf,lun,"_pd_peak_d_spacing"
    printf,lun,"_pd_peak_intensity"
    printf,lun,peaks
endif

free_lun,lun
return,1b

end;function SaveCif
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SaveCel,file,ptr,O=O

if openw_safe(lun,file) then return,0b

if keyword_set(O) then p=(*(*ptr).globalparamIO.O)[0] $
else p=(*ptr).globalparamIO.I
if ~ptr_valid(p) then begin
    free_lun,lun
    return,0b
endif
celparam=(*p)[indCelParam(ptr)]
printf,lun,celparam,format='("CELL",6(f10.4))'

ASU=CopystrucASU(ptr,O=O)
for i=0l,ASU.nasupos-1 do begin
    k0=ASU.reverse_indices[i]
    k1=ASU.reverse_indices[i+1]-1
    k=k0
    
    label=strmid(ASU.asupos[k].name,0,4)

    str=ZElement(*ASU.asupos[k].Z,Ion=*ASU.asupos[k].Ion,Repeats=*ASU.asupos[k].Repeats)
    n=12-strlen(str)-strlen(label)
    if n gt 0 then str=string(replicate(32b,n))+str
    str=label+str
    printf,lun,str,ASU.asupos[k].xyz,ASU.asupos[k].SOF,$
        ASU.asupos[k].B,format='(A12,5(f10.4))'

    dummy=strarr(3)
    for k=k0+1,k1 do begin
        label=strmid(ASU.asupos[k].name,0,4)
        str=ZElement(*ASU.asupos[k].Z,Ion=*ASU.asupos[k].Ion,Repeats=*ASU.asupos[k].Repeats)
        n=12-strlen(str)-strlen(label)
        if n gt 0 then str=string(replicate(32b,n))+str
        str=label+str

        str=ZElement(*ASU.asupos[k].Z,Ion=*ASU.asupos[k].Ion,Repeats=*ASU.asupos[k].Repeats)
        printf,lun,str,dummy,ASU.asupos[k].SOF,$
            ASU.asupos[k].B,format='(A12,3(A10),2(f10.4))'
    endfor

endfor


tmp=sgdata((*(*ptr).pstrinfo).sghash,0)
SG=stringr(tmp.num)
if tmp.ext ne '' then begin
    case tmp.ext of
    'R':SG+=' 1'
    'H':SG+=' 2'
    else:SG+=' '+tmp.ext
    endcase
endif

printf,lun,SG,format='("RGNR ",A)'
free_lun,lun
return,1b
end;function SaveCel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadCel,file,SmallDev,nopos=nopos,error=error

; .cel file example:
;CELL   10.0000   10.0000   10.0000    90.000    90.000    90.000
;Na       Na+    0.0000    0.0000    0.0000    0.9000     0.3
;          K+                                  0.1000
;Cl       Cl-    0.5000    0.5000    0.5000
;RGNR   1
;
; returns:
;     SG: {...}
;     celparam: [a,b,c,aa,bb,cc]
;    ASU: {...}

error=1b
if openr_safe(lun,file) then return,0b
line=''
readf,lun,line

; ----Read Cell param----
SearchStr='CELL'
len=strlen(SearchStr)
while (STRUPCASE(strmid(line,0,len)) ne SearchStr) and (not eof(lun)) do readf,lun,line
if eof(lun) then return,0b
line1 = strsplit(line,' ',/extract,count=ct)
if ct lt 7 then return,0b
celparam=fltarr(6)
reads,line1[1:6],celparam

; ----Read atoms----

; first line might be "natoms 2" e.g.
readf,lun,line
SearchStr='NATOM'
len=strlen(SearchStr)
if STRUPCASE(strmid(line,0,len)) eq SearchStr then readf,lun,line

; stop reading when RGNR
SearchStr='RGNR'
len=strlen(SearchStr)
while ((STRUPCASE(strmid(line,0,len)) ne SearchStr) and (not eof(lun))) do begin
    if n_elements(buffer) eq 0 then buffer=line else buffer=[buffer,line]
    readf,lun,line
endwhile

; ----Space group----
line1 = strsplit(line,' ',/extract,count=ct)
ITN=0
reads,line1[1],ITN
SG=spacegroup(ITN,2,error=error,Pconv=Pconv)
free_lun,lun

; ----Convert----
TransformCelparam,celparam,Pconv
Pconv=invert(Pconv) ; for atomic coordinates

; ----Atomic positions----
if ~keyword_set(nopos) then begin
    csep=' '+string(9b)
    for i=0l,n_elements(buffer)-1 do begin
        line1 = strsplit(buffer[i],csep,/extract,count=n)
        if n ge 5 then begin ; new position
            name=line1[0]
            Z=ZElement(line1[1],Ion=Ion,repeats=repeats)
            xyz=rationalStringToFloat(line1[2:4])
            xyz=(Pconv##reform([xyz,1],1,4))[0:2]
            if n gt 5 then SOF=float(line1[5]) else SOF=1.
            if n gt 6 then B=float(line1[6]) else B=0.
        endif else begin ; substition
            Z=ZElement(line1[0],Ion=Ion,repeats=repeats)
            if n gt 1 then SOF=float(line1[1]) else SOF=1.
            if n gt 2 then B=float(line1[2]) else B=0.
        endelse

        AddASUpos,ASU,name,Z,Ion,repeats,xyz,SOF,B
    endfor
    AddASUpos,ASU
endif 

; ----Finish----
if ~error then ResetWyckASU,SG,asu,SmallDev
return,create_struct('ASU',ASU,$ ; asymmetric unit
                'SG',SG,$    ; space group information
                'celparam',celparam); 6 celparameters: A,A,A,deg,deg,deg
end;function ReadCel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LoadStructure,path,file,SmallDev,nopos=nopos,error=error

error=1b
file=CAPDIALOG_PICKFILE(path=path,file=file,filter=['*.cel','*.aid','*.cif'],title='Load Structure File')
if file eq '' then return,0B

tmp=CutPath(file,path=path,ext=ext)
case strlowcase(ext) of
'.cel':    ret=ReadCel(file,SmallDev,nopos=nopos)
'.cif':    ret=ReadCif(file,SmallDev,nopos=nopos)
'.aid':    ret=ReadAID(file)
else:     return,0B
endcase

if size(ret,/type) ne 8 then return,0B

error=0b
return,ret
end;function LoadStructure
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SaveStructure,path,file,ptr,O=O,peaks=peaks,peaknr=peaknr

error=1b
file=CAPDIALOG_PICKFILE(path=path,file=file,filter=['*.cel','*.cif'],title='Save Structure File')
if file eq '' then return,error

tmp=CutPath(file,path=path,ext=ext)
case strlowcase(ext) of
'.cel':    error=SaveCel(file,ptr,O=O)
'.cif':    error=SaveCif(file,ptr,O=O,peaks=peaks,peaknr=peaknr)
else:
endcase

return,error
end;function SaveStructure
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadDM3SimpleTypes,lun,encoding,iencoding,dummy=dummy
; field 1 = type

if keyword_set(dummy) then begin
    iencoding++
    return,0
endif

case encoding[iencoding++] of
2: data=0    ;SHORT
3: data=0L   ;LONG
4: data=0u   ;USHORT
5: data=0ul  ;ULONG
6: data=0.   ;FLOAT
7: data=0d   ;DOUBLE
8: data=0b   ;BOOLEAN
9: data=0b   ;CHAR
10: data=0b  ;OCTET
endcase

readu,lun,data

return,data
end;function ReadDM3SimpleTypes
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadDM3String,lun,encoding,iencoding,dummy=dummy
; field 1 = number of 2byte Unicode characters (UTF-16)

if keyword_set(dummy) then begin
    iencoding++
    return,''
endif

n=encoding[iencoding++]
if n eq 0 then return,''

str=bytarr(2*n)
readu,lun,str
return,string(str)

; Do this later, because endian swap:
;Z=hexadecimal
;!Z= Hershey-font format command for 2-byte Unicode char
;string(uintarr(n), FORMAT='('+stringr(n)+'("!Z(",Z04,")"))')

end;function ReadDM3String
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FORWARD_FUNCTION ReadDM3Element

function ReadDM3Struct,lun,encoding,iencoding,dummy=dummy
; field 1 = struct_namelength
; field 2 = number of fields (n)
;        (field_namelength, fieldtype) x n
;        where each of the n fieldtypes is one of the EncodedTypes.

if keyword_set(dummy) then begin
    iencoding++
    n=encoding[iencoding++]
    for i=0l,n-1 do begin
        iencoding++
        temp=ReadDM3Element(lun,encoding,iencoding,dummy=dummy)
    endfor
    return,{tag0:0}
endif

temp=encoding[iencoding++]
if temp ne 0 then begin
    strucname=string(bytarr(temp))
    readu,lun,strucname
endif else strucname=''
n=encoding[iencoding++]

for i=0l,n-1 do begin
    temp=encoding[iencoding++]
    if temp ne 0 then begin
        tagname=string(bytarr(temp))
        readu,lun,tagname
        tagname=strcompress(tagname,/remove_all)
        tagname=IDL_VALIDNAME(tagname,/convert_all)
    endif else tagname='tag'+stringr(i)

    datai=ReadDM3Element(lun,encoding,iencoding,dummy=dummy)

    if i eq 0 then data=create_struct(tagname,datai) $
    else begin
        ind=where(tag_names(data) eq strupcase(tagname),ct)
        if ct ne 0 then tagname+=stringr(i)
        data=create_struct(data,tagname,datai)
    endelse
endfor

return,data
end;function ReadDM3Struct
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadDM3Array,lun,encoding,iencoding,dummy=dummy
; field 1 = array_type, which is one of the EncodedTypes.
;    ...
; field last = array_length.

; First a dummy read (lun pointer doesn't change) to get the array length
iencoding_=iencoding
temp=ReadDM3Element(lun,encoding,iencoding,/dummy)
n=encoding[iencoding++]
if n eq 0 or keyword_set(dummy) then return,[0]
iencoding=iencoding_

; Read first element
point_lun,-lun,pos
data=ReadDM3Element(lun,encoding,iencoding)
n=encoding[iencoding++]
if n eq 1 then return,data

; Add other elements
if encoding[iencoding_] eq 15 or encoding[iencoding_] eq 20 then begin
    for i=1,n-1 do begin
        iencoding=iencoding_
        data=[data,ReadDM3Element(lun,encoding,iencoding)]
    endfor
endif else begin
    ; Read array of simple type at once
    point_lun,lun,pos
    data=make_array(n,value=data)
    readu,lun,data
endelse

return,data
end;function ReadDM3Array
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadDM3Element,lun,encoding,iencoding,dummy=dummy

case encoding[iencoding] of
15: begin;STRUCT
    iencoding++
    return,ReadDM3Struct(lun,encoding,iencoding,dummy=dummy)
    endcase
18: begin;STRING
    iencoding++
    return,ReadDM3String(lun,encoding,iencoding,dummy=dummy)
    endcase
20: begin ;ARRAY
    iencoding++
    return,ReadDM3Array(lun,encoding,iencoding,dummy=dummy)
    endcase
else:return,ReadDM3SimpleTypes(lun,encoding,iencoding,dummy=dummy) ;Simple type
endcase

end;function ReadDM3Element
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadDM3TagType,lun
;1.  4 bytes, equals to "%%%%".
;2.  4 bytes, n = length of definiton of encoded type:
;    for a simple type this will = 1,
;    for a string this will = 2,
;    an array of a simple type will = 3,
;    structs have 1+2*f where f=number of fields in struct
;
;NB arrays of structs and arrays of arrays etc. will require additional
;fields for the complete definition (see Further details: TagType 3.
;EncodedType)
;
;3.  4 x n bytes: where n is number of encoded types defined above;
;    each 4 byte value helps defines the Encoded Type, either as an ID
;    to be looked up in the table below or (for complex types) indicating
;    field numbers etc as described below
;4.  data, its size depends on encoded types.

N=0L
str=bytarr(4)

readu,lun,str
if string(str) ne '%%%%' then return,'error'

readu,lun,N
encoding=lonarr(N)
readu,lun,encoding
iencoding=0

data=ReadDM3Element(lun,encoding,iencoding)

return,data
end;function ReadDM3TagType
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FORWARD_FUNCTION ReadDM3TagGroup

function ReadDM3TagEntry,lun,tagname
;1.  1 byte, specifying whether it is data (21) or tag group (20)
;2.  2 bytes, length of tag's label
;3.  n bytes, tag's label as a string
;4.  tag instance: if tag is a group, it's a TagGroup instance.
;    if tag is a data tag, it's a TagType instance.
I8=0b
I16=0
readu,lun,I8
bdata=I8 eq 21b
readu,lun,I16
if I16 ne 0 then begin
    tagname=bytarr(I16)
    readu,lun,tagname
    tagname=string(tagname)
    tagname=strcompress(tagname,/remove_all)
    tagname=IDL_VALIDNAME(tagname,/convert_all)
endif

if bdata then data=ReadDM3TagType(lun) $
else data=ReadDM3TagGroup(lun)

return,data

end;function ReadDM3TagEntry
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadDM3TagGroup,lun
;1.  1 byte, is this group sorted?
;2.  1 byte, is this group open?
;3.  4 bytes, number of tags in this group
;4.  all tag entries. each tag entry is a TagEntry instance.

I8=0b
N=0L

readu,lun,I8
readu,lun,I8
readu,lun,N

if N eq 0 then return,'empty'
defaulttagnames='tag'+stringr(indgen(N))

tagname=defaulttagnames[0]
temp=ReadDM3TagEntry(lun,tagname)
data=create_struct(tagname,temp)

for i=1,N-1 do begin
    tagname=defaulttagnames[i]
    temp=ReadDM3TagEntry(lun,tagname)
    ind=where(tag_names(data) eq strupcase(tagname),ct)
    if ct ne 0 then tagname+=stringr(i)

    data=create_struct(data,tagname,temp)
endfor

return,data
end;function ReadDM3TagGroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DM3taglabel,data,topbranch
if size(data,/type) ne 8 then begin
    help,data,output=output
    Mbranch = WIDGET_TREE(topbranch,value=strcompress(output[0]))
    return
endif

tags=tag_names(data)
n=n_tags(data)
for i=0l,n-1 do begin
    Mbranch = WIDGET_TREE(topbranch,value=tags[i],/folder)
    DM3taglabel,data.(i),Mbranch
endfor
end;pro DM3taglabel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DM3INFO_EVENT,ev
end;DM3INFO_EVENT
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadDM3,path2,info=info
;http://rsb.info.nih.gov/ij/plugins/DM3Format.gj.html
;
;1.  4 bytes, image file version, should be 3
;2.  4 bytes, number of bytes in the file
;3.  4 bytes, byte-ordering of the tag data (0 = bigendian, 1 = littleendian)
;4.  a single TagGroup instance including all the tags and data

I32=0L

; All data "big-endian" byte order, except tag data -> LEndian
if OPENR_safe(lun, path2, /SWAP_IF_LITTLE_ENDIAN) then return,0
readu,lun,I32
if I32 ne 3 then return,0
readu,lun,I32
readu,lun,I32
LEndian=I32 ne 0
data=ReadDM3TagGroup(lun)
FREE_LUN, lun

if LEndian then SWAP_ENDIAN_INPLACE, data

; Extract data: take image with biggest size
n=n_tags(data.imagelist)
s=lonarr(2,n)
for i=0l,n-1 do begin
    temp=data.imagelist.(i).imagedata.dimensions
    s[*,i]=[temp.(0),temp.(1)]
endfor
m=max(s[0,*]*s[1,*],i)

if keyword_set(info) then begin
    if info.showheader then begin
        base=widget_base(/row,title=path2)
        topbranch = WIDGET_TREE(base,xsize=300,ysize=300)
        DM3taglabel,data,topbranch
        WIDGET_CONTROL, base, /REALIZE
        Xmanager,'DM3info',base,/NO_BLOCK
    endif
endif

return,reform(data.imagelist.(i).imagedata.data,s[0,i],s[1,i])
end;function ReadDM3
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadRAW,path2,info=info
if OPENR_safe(lun, path2) then return,0

tiff=make_array(info.binxs,info.binys,type=info.bintype)

; Read image
READU, lun, tiff
SWAP_ENDIAN_INPLACE, tiff,SWAP_IF_BIG_ENDIAN =info.binlendian,SWAP_IF_LITTLE_ENDIAN = info.binlendian eq 0
FREE_LUN, lun

if info.binnxpad3 ne 0 then begin
    if (info.binxs mod info.binnxpad3) eq 0 then begin
        ; Indices of boundary columns of the different blocks
        nbounds=info.binnxpad3-1
        wblock=info.binxs / info.binnxpad3
        boundind=wblock*lindgen(nbounds)+wblock
        boundind=[boundind,boundind-1]
        ind=[lindgen(info.binxs),boundind]
        ind=ind[sort(ind)]
        
        ; Value/2 for boundary columns
        tiff[boundind,*]/=2
        
        ; Duplicate the boundary columns
        tiff=temporary(tiff[ind,*])
        
    endif
endif

return,tiff
end;function ReadRAW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadSIFblock,lun
line=''
for i=0l,25 do readf,lun,line

; Image area
readf,lun,line
tmp=strsplit(line,' ',/extract,count=n)
imageArea=(long(tmp[n-8:n-1]))[[0,3,5,2,1,4]]

; Frame area
readf,lun,line
tmp=strsplit(line,' ',/extract,count=n)
frameArea=(long(tmp[1:*]))[[0,3,2,1,5,4]]
frameBins=frameArea[4:5]
frameArea=frameArea[0:3]

; Image
tiff=fltarr(frameArea[2:3]-frameArea[0:1]+1)
readf,lun,line
readu,lun,tiff

return,tiff
end;function ReadSIFblock
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadSIF,path2,xdinr=xdinr,nfiles=nfiles
;Andor Technology Multi-Channel File
;http://www.mathworks.com/matlabcentral/fileexchange/11224

if OPENR_safe(lun, path2) then return,0

line=''
readf,lun,line
readf,lun,line

if n_elements(xdinr) eq 0 then xdinr=0

nfiles=0l
repeat begin
    img=ReadSIFblock(lun)
    if nfiles eq xdinr then tiff=temporary(img)
    nfiles++
    readf,lun,line
    readf,lun,line
endrep until fix(line) eq 0

FREE_LUN, lun
return,tiff
end;function ReadSIF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RigakuRaxisHeader,hversion
; http://usaxs.xor.aps.anl.gov/staff/ilavsky/nika.html

if n_params() eq 0 then hversion=2
case hversion of
0:     begin ; Nika old header
    XrayFocus=bytarr(10)
    Reserve2=bytarr(58)
    cmnt=''
    smpl=''
    init=''
    ipus=''
    dexp=''
    expn=''
    posx=''
    posy=''
    xray=''
    res58=bytarr(100)
    res59=bytarr(100)
    res60=bytarr(56)
    endcase
1:    begin ; Nika new header
    XrayFocus=bytarr(12)
    Reserve2=bytarr(56)
    cmnt=bytarr(20)
    smpl=bytarr(20)
    init=0l
    ipus=0l
    dexp=0l
    expn=0l
    posx=lonarr(20)
    posy=lonarr(20)
    xray=0
    res58=bytarr(68)
    res59=''
    res60=''
    endcase
2:    begin ; modifications to 1
    XrayFocus=bytarr(12)
    Reserve2=bytarr(56)
    cmnt=bytarr(20)
    smpl=bytarr(20)
    init=0l
    ipus=0l
    dexp=0l
    expn=0l
    posx=lonarr(20)
    posy=lonarr(20)
    xray=0
    res58=bytarr(68)
    res59=''
    res60=''
    endcase
endcase

return,{DeviceName:bytarr(10),$
    Version:bytarr(10),$
    CrystalName:bytarr(20),$
    CrystalSystem:bytarr(12),$
    LatticeA:0.,$
    LatticeB:0.,$
    LatticeC:0.,$
    LatticeAlpha:0.,$
    LatticeBeta:0.,$
    LatticeGamma:0.,$
    SpaceGroup:bytarr(12),$
    MosaicSpread:0.,$
    Memo:bytarr(80),$
    Reserve1:bytarr(84),$

    Date_:bytarr(12),$
    MeasurePerson:bytarr(20),$
    Xraytarget:bytarr(4),$
    Wavelength:0.,$
    Monochromator:bytarr(20),$
    MonocromatorDq:0.,$
    Collimator:bytarr(20),$
    Filter:bytarr(4),$
    CameraLength_mm:0.,$
    XrayTubeVoltage:0.,$
    XrayTubeCurrent:0.,$
    XrayFocus:XrayFocus,$
    XrayOptics:bytarr(80),$
    CameraShape:0l,$
    WeissenbergOscillation:0.,$
    Reserve2:Reserve2,$

    MountAxis:bytarr(4),$
    BeamAxis:bytarr(4),$
    something7:0.,$
    StartSomething :0.,$
    EndSomething:0.,$
    TimesOfOscillation:0l,$
    ExposureTime:0.,$
    DirectBeamPositionX:0.,$
    DirectBeamPositionY:0.,$
    Something8 :0.,$
    Something9:0.,$
    Something10:0.,$
    Something11:0.,$
    Reserve3:bytarr(100),$
    Reserve3a:bytarr(100),$
    Reserve3b:bytarr(4),$

    xDirectionPixNumber:0l,$
    yDirectionPixNumber:0l,$
    xDirectionPixelSizeMM:0.,$
    yDirectionPixelSizeMM:0.,$
    RecordLengthByte:0l,$
    NumberOfRecord:0l,$
    ReadStartLine:0l,$
    IPNumber:0l,$
    OutPutRatioHighLow:0.,$
    FadingTime1:0.,$
    FadingTime2:0.,$
    HostComputerClass:bytarr(10),$
    IPClass:bytarr(10),$
    DataDirectionHorizontal:0l,$
    DataDirectionVertical:0l,$
    DataDirectionFrontBack:0l,$

    shft:0.,$;         /* pixel shift, R-AXIS V */
    ineo:0.,$;         /* intensity ratio E/O R-AXIS V */
    majc:0l,$;         /* magic number to indicate next values are legit */
    naxs:0l,$;         /* Number of goniometer axes */
    gvec1:fltarr(5),$;   /* Goniometer axis vectors */
    gvec2:fltarr(5),$;   /* Goniometer axis vectors */
    gvec3:fltarr(5),$;   /* Goniometer axis vectors */
    gst:fltarr(5),$;       /* Start angles for each of 5 axes */
    gend:fltarr(5),$;      /* End angles for each of 5 axes */
    goff:fltarr(5),$;      /* Offset values for each of 5 axes */ 
    saxs:0l,$;         /* Which axis is the scan axis? */
    gnom:bytarr(40),$;     /* Names of the axes (space or comma separated?) */

; Most of below is program dependent.  Different programs use
; this part of the header for different things.  So it is essentially 
; a big "common block" area for dumping transient information.

   file:bytarr(16),$;     /* */
   cmnt:cmnt,$;     /* */
   smpl:smpl,$;     /* */
   iext:0l,$;         /* */
   reso:0l,$;         /* */
   save_:0l,$;         /* */
   dint:0l,$;         /* */
   byte:0l,$;         /* */
   init:init,$;         /* */
   ipus:ipus,$;         /* */
   dexp:dexp,$;         /* */
   expn:expn,$;         /* */
   posx:posx,$;     /* */
   posy:posy,$;     /* */
   xray:xray,$;         /* */
   res51:bytarr(100),$;    /* reserved space for future use */
   res52:bytarr(100),$;    /* reserved space for future use */
   res53:bytarr(100),$;    /* reserved space for future use */
   res54:bytarr(100),$;    /* reserved space for future use */
   res55:bytarr(100),$;    /* reserved space for future use */
   res56:bytarr(100),$;    /* reserved space for future use */
   res57:bytarr(100),$;    /* reserved space for future use */
   res58:res58,$;    /* reserved space for future use */
   res59:res59,$;    /* reserved space for future use */
   res60:res60};    /* reserved space for future use */
end;function RigakuRaxisHeader
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadRigakuRaxis,path2
if OPENR_safe(lun, path2,/SWAP_IF_LITTLE_ENDIAN) then return,0

; Get header
h=RigakuRaxisHeader()
readu,lun,h

; Get image format
nx=h.xDirectionPixNumber
ny=h.yDirectionPixNumber
nbytes=h.RecordLengthByte/h.NumberOfRecord

case nbytes of
1: tiff=bytarr(nx,ny)
2: tiff=intarr(nx,ny)
4: tiff=lonarr(nx,ny)
endcase

; Read image
offset=h.RecordLengthByte
if offset lt 250 then offset=2*nx
point_lun,lun,offset
readu,lun,tiff

free_lun,lun

return,tiff
end;function ReadRigakuRaxis
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ParseSPEHeader,new_head
print,'Header size: ',strucbytesize(new_head)
print,'----Dimensions----'
print,'XY-dimension: ', new_head.x_chip, new_head.y_chip
print,'XY-virtual:   ', new_head.x_vir, new_head.y_vir
print,'XY-actual:    ', new_head.x_actual, new_head.y_actual
print,'XY-prescan:   ', new_head.x_pre, new_head.y_pre
print,'XY-postscan:  ', new_head.x_post, new_head.y_post
print,'NFILES:       ', new_head.nfiles

print,'----Time----'
day   = float(new_head.date1)
monthstring = new_head.date2
monthstrarr = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
montharr = where(monthstrarr EQ monthstring)
month = montharr[0] + 1
year  = float(new_head.date3)
DATE_OBS =  string(day,  format='(I2)') + '/' + $
            string(month,format='(I2)') + '/' + $
            string(year, format='(I4)')
print,'DATE-OBS:     ', DATE_OBS
print,'Exposure, s: ', new_head.EXP_TIME
LT_START = new_head.ltime1 + ':' +  new_head.ltime2 + ':' +  new_head.ltime3
UT_START = new_head.utime1 + ':' +  new_head.utime2 + ':' +  new_head.utime3
print,'LT-START:     ', LT_START;(HH:MM:SS)  Local time?
print,'UT-START:     ', UT_START;(HH:MM:SS)  Universal Time
print,'TIMEUNITS:   ', new_head.TIMUNITS

print,'----Properties----'
CASE new_head.data_type OF
  0: print, 'DATATYPE: FLOATING POINT'
  1: print, 'DATATYPE: LONG INTEGER'
  2: print, 'DATATYPE: INTEGER'
  3: print, 'DATATYPE: UNSIGNED INTEGER'
ELSE: print,'DATATYPE: UNKNOWN'
ENDCASE
print,'Back. sub.? ', new_head.back_sub
print,'GAIN:       ', new_head.gain
print,'CCD-TEMP:    ', new_head.temper

print,'----Comments----'
print,'Comments_1: ', new_head.comment1
print,'Comments_2: ', new_head.comment2
print,'Comments_3: ', new_head.comment3
print,'Comments_4: ', new_head.comment4
print,'Comments_5: ', new_head.comment5
end;pro ParseSPEHeader
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadSPE,path,xdinr=xdinr,nfiles=nfiles,getnfiles=getnfiles
if OPENR_safe(lun, path) then return,0

; ----Read Header----
empty_comment = ''
atmp = ' '
a01 = string(atmp, format='(A1)')
a02 = string(atmp, format='(A2)')
a03 = string(atmp, format='(A3)')
a04 = string(atmp, format='(A4)')
a06 = string(atmp, format='(A6)')
a07 = string(atmp, format='(A7)')
a10 = string(atmp, format='(A10)')
a20 = string(atmp, format='(A20)')
a40 = string(atmp, format='(A40)')
a50 = string(atmp, format='(A50)')
a54 = string(atmp, format='(A54)')
a56 = string(atmp, format='(A56)')
a80 = string(atmp, format='(A80)')

new_head = {skip1:a06, x_chip:1, skip2:0, EXP_TIME:0.0, x_vir:1,y_vir:1,y_chip:1, $
            date1:a02, date2:a03, date3:a04, date4:a01, skip3:a06, temper:0.0, skip4:0,$
            x_actual:0, skip5:a54, x_pre:0, x_post:0, y_pre:0, y_post:0, skip6:0, $
            data_type:0, skip7:a40, back_sub:0,  skip8:a20, $
            ltime1:a02,ltime2:a02,ltime3:a02,ltime4:a01, $
            utime1:a02,utime2:a02,utime3:a02,utime4:a01, $
            TIMUNITS:0,$
            skip9:a10, gain:0, comment1:a80, comment2:a80, $
            comment3:a80, comment4:a80, comment5:a80, $
            skip10:a56, y_actual:0, $
            skip11:replicate(' ',788), nfiles:1,$
            skip12:replicate(' ',2652)}

readu,lun, new_head
;ParseSPEHeader,new_head

; ----Find out the dimensions of the data----
nx = new_head.x_actual
ny = new_head.y_actual

CASE new_head.data_type OF
    0:     begin
        ima = fltarr(nx,ny)
        nbytes=4l
        endcase
    1:     begin
        ima = lonarr(nx,ny)
        nbytes=4l
        endcase
    2:     begin
        ima = intarr(nx,ny)
        nbytes=2l
        endcase
    3:     begin
        ima = uintarr(nx,ny)
        nbytes=2l
        endcase
ELSE:     begin
        free_lun,lun
        return,0
        endelse
ENDCASE

nbytes*=long(nx)*ny
point_lun,-lun,pos
stot=(fstat(lun)).size
sdata=stot-pos
nfiles=sdata/nbytes
if nfiles ne new_head.nfiles then begin
    free_lun,lun
    return,0
endif

; ----Point to image----
if n_elements(xdinr) ne 0 then begin
    xdinr>=0
    xdinr<=nfiles-1
    if xdinr ne 0 then begin
        posnew=pos+xdinr*nbytes
        point_lun,lun,posnew
    endif
endif

; ----Only looking for the number of files?----
if keyword_set(getnfiles) then begin
    FREE_LUN, lun
    return,nfiles
endif

; ----Get data----
readu, lun, ima
FREE_LUN, lun

return,ima
end;function ReadSPE
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadCCDtiff,path
CATCH, Error_status
IF Error_status eq 0 THEN begin
    tiff = read_tiff(path,orientation=orientation)
    ; type returned:
    ;    byte(8bit unsigned)
    ;    int(16bit signed)
    ;    int + /unsigned(16bit unsigned, stored in 32bit signed format)
    ;    long(32bit signed)
    ;    float
    
    
;    i=(where(tiff[0:7781] ne 0))[0]
;    print,tiff[i:i+10],format='(11Z04)'
;    stop
    
    
    o=[0,7,2,5,0,1,6,3,4]
    tiff=rotate(temporary(tiff),o[orientation])
endif else begin
    CATCH, Error_status
    IF Error_status ne 0 THEN return,0
    tiff = read_tiff_structure(path,error=error,/format)
    if error then tiff=0 else tiff=tiff.ifd0.image
endelse

return,tiff
end;function ReadCCDtiff
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadCCDjpg,path
CATCH, Error_status
IF Error_status ne 0 THEN return,0
;READ_JPEG,path,tiff,tbl,COLORS=!D.N_COLORS-1,/TWO_PASS_QUANTIZE
READ_JPEG,path,tiff,true=3

s=size(tiff,/dim)
if n_elements(s) le 2 then return,tiff
if n_elements(s) ne 3 then return,0
tiff=0.299 * tiff[*,*,0] + 0.587 * tiff[*,*,1] + 0.114 * tiff[*,*,2]

return,tiff
end;function ReadCCDjpg
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadCCDbmp,path
CATCH, Error_status
IF Error_status ne 0 THEN return,0
tiff=READ_BMP(path)

s=size(tiff,/dim)
if n_elements(s) le 2 then return,tiff
if n_elements(s) ne 3 then return,0

ind=where(s eq 3,ct)
if ct ne 1 then return,0
case ind[0] of
0: tiff=reform(0.299 * tiff[0,*,*] + 0.587 * tiff[1,*,*] + 0.114 * tiff[2,*,*])
1: tiff=reform(0.299 * tiff[*,0,*] + 0.587 * tiff[*,1,*] + 0.114 * tiff[*,2,*])
2: tiff=0.299 * tiff[*,*,0] + 0.587 * tiff[*,*,1] + 0.114 * tiff[*,*,2]
endcase
    
return,tiff
end;function ReadCCDbmp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadMarCCD,file

; http://www.sb.fsu.edu/~xray/Manuals/marCCD165header.html
; Mar ccd frames (always 16bit unsigned)
;|-- 1024 bytes TIFF HEADER -------------|
;|-- 3072 byte frame_header structure ---|
;|-- nfast*nslow*depth byte image -------|
; frame_header structure:
;File/header format parameters (256 bytes)
;Data statistics (128)
;More statistics (256)
;Goniostat parameters (128 bytes)
;Detector parameters (128 bytes)
;X-ray source and optics parameters (128 bytes)
;File parameters (1024 bytes)
;Dataset parameters (512 bytes)
;Empty (512 bytes)

;case !VERSION.OS of
;    'linux' : OPENR, lun, path2, /GET_LUN
;    'Win32' : OPENR, lun, path2, /GET_LUN, BUFSIZE=0;?
;   else : OPENR, lun, path2, /GET_LUN, /SWAP_ENDIAN;?
;endcase

if OPENR_safe(lun, file) then return,0

;POINT_LUN, lun, 1024
;tmp=bytarr(128)
;readu,lun,tmp
;tmp=bytarr(128)
;readu,lun,tmp
;tmp=bytarr(64)
;readu,lun,tmp
;tmp=bytarr(32)
;readu,lun,tmp
;tmp=bytarr(32)
;readu,lun,tmp
;tmp=bytarr(32)
;readu,lun,tmp

POINT_LUN, lun, 4096;Skip header
    
Result = QUERY_TIFF (file, Info)
if Result then begin
    xs=(Info.dimensions)[0]
    ys=(Info.dimensions)[1]
endif else begin
    Info=fstat(lun)
    nPixels=(Info.size-4096)/2
    xs=ulong(sqrt(nPixels))
    ys=xs
endelse
    
tiff =UINTARR(xs,ys)
READU, lun, tiff
FREE_LUN, lun

return,tiff

end;function ReadMarCCD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function decodeccp4_extractbits,raw,bitoff,nbits

byte1=bitoff/8
bit1=bitoff mod 8
    
bitoff+=nbits
bit2=bitoff-1
byte2=bit2/8
bit2 mod= 8

nbytes=byte2-byte1+1

return,byte(strmid(string(raw[byte2:byte1:-1],format='('+strtrim(nbytes,2)+'b08)'),7-bit2,nbits))-(byte('0'))[0]
end;function decodeccp4_extractbits
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro decodeccp4, raw, nx, ny

nraw      = n_elements(raw)
bitoff    = 0ULL
max       = 8ULL*nraw-5 ; there should be at least 6 bits left to continue

npixels   = long(nx) * ny; number of pixels
img       = intarr(npixels); uncompressed image
pixel     = 0L ; current image index

bitdecode = [0L, 4L, 5L, 6L, 7L, 8L, 16L, 32L] ; array translating the number of bits per error

bitconv3  = [4b,2b,1b]
bitconv32 = 2l^reverse(lindgen(32))

while bitoff lt max do begin
    ; Extract number of errors in chunk
    nerrors = ishft(1L,total(decodeccp4_extractbits(raw,bitoff,3)*bitconv3,/pres)) ; maximal 128
    
    ; Extract number of bits per error
    nerrorbits = bitdecode[total(decodeccp4_extractbits(raw,bitoff,3)*bitconv3,/pres)] ; maximal 32
    
    ; Extract errors
    if nerrorbits eq 0 then errors=lonarr(nerrors) else begin
        errors=decodeccp4_extractbits(raw,bitoff,nerrors*nerrorbits) ; maximal 4096 bits
        errors=reform(errors,nerrorbits,nerrors,/overwrite)
        ind=where(errors[0,*],ct)
        errors=total(errors*rebin(bitconv32[32-nerrorbits:*],nerrorbits,nerrors),1,/pres) ; maximal 128
        if ct ne 0 then errors[ind] or= not (ishft(1L,nerrorbits)-1L)
        errors=reverse(errors)
    endelse
    
    ; Fill the image with these errors
    for i=0l,nerrors-1 do begin
        pixeli=pixel+i
        
        ;if errors[i] ne 0 then print,'error:',errors[i]
        
        if (pixeli gt nx) then begin
            img[pixeli] = (errors[i] + $
                          (2L + img[pixeli-1] + img[pixeli-nx+1] + $
                           img[pixeli-nx] + img[pixeli-nx-1]) / 4L) ; what about overflow?
        endif else if (pixeli ne 0) then begin
            img[pixeli] = img[pixeli-1] + long(errors[i]) ; what about overflow?
        endif else begin
            img[pixeli++] = errors[i]
        endelse
    endfor
    
    pixel+=nerrors

endwhile

raw=temporary(uint(img))

end;pro decodeccp4
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro decodeccp4org, raw, nx, ny
; Adapted from Mark Rivers, 2006
; 
; IT-G: 5.6.3. Compression schemes
;       5.6.3.2. CCP4-style compression
;       
; Pixel "error" is the different between the pixel value and the average of 4 neighbours:
;     E[x,y] = I[x,y] - (I[x-1,y-1] + I[x-1,y] + I[x,y-1] + I[x-1,y+1])/4


npixels   = long(nx) * ny; number of pixels
img       = intarr(npixels); uncompressed image

spill     = 0L ; receives the raw bytes
spillbits = 0L ; number of bits in the spill

window    = 0L ; receives the spill bits
valids    = 0L ; number of bits in the window

next      = 0L ; current raw index
pixel     = 0L ; current image index
error     = 0L ; current error

bitnum    = 0L ; number of bits per error
pixnum    = 0L ; number of errors in a chunk
bitdecode = [0L, 4L, 5L, 6L, 7L, 8L, 16L, 32L] ; array translating the number of bits per error

bitmask   = ['00000000'xl, '00000001'xl, '00000003'xl, '00000007'xl, $
             '0000000F'xl, '0000001F'xl, '0000003F'xl, '0000007F'xl, $
             '000000FF'xl, '000001FF'xl, '000003FF'xl, '000007FF'xl, $
             '00000FFF'xl, '00001FFF'xl, '00003FFF'xl, '00007FFF'xl, $
             '0000FFFF'xl, '0001FFFF'xl, '0003FFFF'xl, '0007FFFF'xl, $
             '000FFFFF'xl, '001FFFFF'xl, '003FFFFF'xl, '007FFFFF'xl, $
             '00FFFFFF'xl, '01FFFFFF'xl, '03FFFFFF'xl, '07FFFFFF'xl, $
             '0FFFFFFF'xl, '1FFFFFFF'xl, '3FFFFFFF'xl, '7FFFFFFF'xl, $
             'FFFFFFFF'xl]

   while (pixel lt npixels) do begin
      if (valids lt 6) then begin
         if (spillbits gt 0)  then begin
            window or= ishft(spill, valids) ; merge the spill bits into the window, not overwriting the valid bits in window
            valids += spillbits ; valid bits in window increases with the spill bits
            spillbits = 0L ; all spill bits are used
         endif else begin
            spill = long(raw[next++]) ; get the next raw byte
            spillbits = 8L ; spill contains therefore 8 bits
         endelse
      endif else begin
           ; Bits 0-2 = number of errors in this chunk
         pixnum = ishft(1, (window and bitmask[3])) ; 2^n (maximal 2^7)
         window = ishft(window, -3)
         
         ; Bits 3-5 = number of bits per error
         bitnum = bitdecode[window and bitmask[3]]
         window = ishft(window, -3)
         valids -= 6L ; we used 6 bits of the window
         
         ; Loop over the pixels in a chunk
         while ((pixnum gt 0) and (pixel lt npixels)) do begin
            if (valids lt bitnum) then begin
                ; The window contains less bits than there are bits in one error
                
               if (spillbits gt 0) then begin
                     ; Some spill bits are not merged in the window: perform the merging
                     
                  window or= ishft(spill, valids) ; merge the spill bits into the window, not overwriting the valid bits in window
                  if ((32L - valids) gt spillbits) then begin
                       ; No window overflow
                     valids += spillbits ; valid bits in window increases with the spill bits
                     spillbits = 0L ; all spill bits are used
                  endif else begin
                       ; Overflow
                     usedbits = 32L - valids ; spill bits merged in window
                     spill = ishft(spill, -usedbits) ; remove those bits
                     spillbits -= usedbits ; adapt the counter accordingly
                     valids = 32L ; window is full
                  endelse
               endif else begin
                     ; All spill bits have been merged in the window: get new bits
                     
                  spill = long(raw[next++]) ; get the next raw byte
                  spillbits = 8L ; the spill contains therefore 8 bits
               endelse
            endif else begin
               ; The window contains enough bits to contain an error
               ; 
               ; Get the error from the window:
               if (bitnum eq 0) then begin
                  error = 0L ; the only number that does not require any bits to be set to 1
               endif else begin
                  error = window and bitmask[bitnum] ; get error from window
                  valids -= bitnum ; remove window bits
                  window = ishft(window, -bitnum) ; remove window bits
                  if ((error and (ishft(1, (bitnum - 1)))) ne 0) then $ ; first bit of the error is 1
                     error or= not bitmask[bitnum] ; fill all unused bits with '1'
               endelse

               ; Calculate "img" at index "pixel", based on previous pixels (if any)
               if (pixel gt nx) then begin
                  img[pixel] = (error + $
                               (2L + img[pixel-1] + img[pixel-nx+1] + $
                                img[pixel-nx] + img[pixel-nx-1]) / 4L) ; what about overflow?
                  ++pixel
               endif else if (pixel ne 0) then begin
                  img[pixel] = img[pixel-1] + long(error) ; what about overflow?
                  ++pixel
               endif else begin
                  img[pixel++] = error
               endelse
               
               ; One pixel done
               --pixnum
            endelse
         endwhile
      endelse
   endwhile

raw=temporary(uint(img))

end;pro decodeccp4org
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadMarIP345,file
; Adapted from Mark Rivers, 2006

if OPENR_safe(lun, file) then return,0

; ---- Read 4096 bytes header ----
head = lonarr(16)
readu, lun, head
if head[0] eq 1234 then swap=0b else swap=1b
if swap then head = swap_endian(head)

text=bytarr(4032)
readu, lun, text
text=strtrim(string(text),2)

header = {nx:              head[1], $
          ny:              head[1], $
          nhigh:           head[2], $
          format:          head[3], $
          collection_mode: head[4], $
          npixels:         head[5], $
          x_pixel_size:    head[6]/1000., $
          y_pixel_size:    head[7]/1000., $
          wavelength:      head[8]/1.e6, $
          distance:        head[9]/1000., $
          phi_start:       head[10]/1000., $
          phi_end:         head[11]/1000., $
          omega_start:     head[12]/1000., $
          omega_end:       head[13]/1000., $
          chi:             head[14]/1000., $
          two_theta:       head[15]/1000., $
          text_lines:      text }

; ---- Read the high-intensity pixel table ----
if (header.nhigh gt 0) then begin
    table_entries = long(header.nhigh/8.0 + .9) * 8
    records = lonarr(2, table_entries)
    readu, lun, records
    high_offsets = reform(records[0, 0:header.nhigh-1])
    high_values  = reform(records[1, 0:header.nhigh-1])
endif

; ---- Read the CBF header ----
;  \nCCP4 packed image, X: %04d, Y: %04d\n
temp = bytarr(37)
readu, lun, temp
if (string(temp[1:4]) ne "CCP4") then begin
    print, "Error, did not find CCP4 header where expected"
    return,0
endif

; ---- Read data ----
if CheckDLMloaded('ReadMarIP345') then begin
    FREE_LUN, lun
    data = uintarr(header.nx*header.ny, /nozero)
    status = XRDUADLM_ReadMarIP345(file, data)
endif else begin
    fs = fstat(lun)
    nread = fs.size - fs.cur_ptr
    data = bytarr(nread)
    readu, lun, data
    FREE_LUN, lun
    
    ; ---- Unpack ----
    decodeccp4org, data, header.nx, header.ny 
    ;decodeccp4, data, header.nx, header.ny ; slower
endelse

; ---- Convert 16bit unit to 32bit long when high pixels ----
if (header.nhigh gt 0) then begin
    data = long(data)
    data[high_offsets]=high_values
endif

; ---- Reform array to 2-D inplace ----
data = reform(data, header.nx, header.ny, /overwrite)
return,data

end;function ReadMarIP345
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadSMART,file

if OPENR_safe(lun, file) then return,0

line=''
readf,lun,line ;This is the header
line=STRCOMPRESS(line)

; *Real Time/ Live Time
bnorm=keyword_set(info)
if bnorm then bnorm=info.bnorm
if bnorm then begin
    pos = STRPOS(line, 'ELAPSDR')
    realtime=0d
    reads,STRMID(line, pos+8,72),realtime
    pos = STRPOS(line, 'ELAPSDA')
    livetime=0d
    reads,STRMID(line, pos+8,72),livetime
endif
    
; Find size of frame
pos = STRPOS(line, 'NROWS')
ys=0L
reads,STRMID(line, pos+8,72),ys
pos = STRPOS(line, 'NCOLS')
xs=0L
reads,STRMID(line, pos+8,72),xs
    
; Read frame with correct number of bytes per pixel (1, 2 or 4)
pos = STRPOS(line, 'NPIXELB')
nbp=0L
reads,STRMID(line, pos+8,72),nbp
case nbp of
1:  begin
     tiff =BYTARR(xs,ys)
     FullScale=255B
    endcase
2:  begin
      tiff =UINTARR(xs,ys)
      FullScale=65535U
    endcase
4:  begin
      tiff =ULONARR(xs,ys)
      FullScale=2UL^32-1
    endcase
endcase
    
; Read Frame
pos = STRPOS(line, 'HDRBLKS')
hs=0L
reads,STRMID(line, pos+8,72),hs
pos = STRPOS(line, 'NEXP')
if pos ne (-1) then begin
    pos2= STRPOS(line, 'CCDPARM')
    orientation=intarr(5)
    reads,STRMID(line, pos+6,pos2-pos-6),orientation
    orientation=orientation[3]
endif else orientation=8
POINT_LUN, lun,hs*512
READU, lun, tiff
    
; Read Overflow table
pos = STRPOS(line, 'NOVERFL')
nOver=0L
reads,STRMID(line, pos+8,72),nOver
if (nOver ne 0) and (not eof(lun)) then begin
    OTable=bytarr(16)
    val=0.
    ind=0L
    line=''
    repeat begin
        readu,lun,OTable
        v=float(string(OTable[0:8]))
        val=[val,float(v)]
        in=long(string(OTable[9:15]))
        ind=[ind,in]
    endrep until eof(lun) or (v eq 0 and in eq 0)
    val=val[1:*]
    ind=ind[1:*]
    nOver=n_elements(val)
    while (val[nOver-1] eq 0) and (ind[nOver-1] eq 0) do begin
        val=val[0:nOver-2]
        ind=ind[0:nOver-2]
        nOver=nOver-1
    endwhile
    
    ; Edit tiff with OverFlow table
    indFullScale=where(tiff ge FullScale,ctFullScale)
    if (ctFullScale ne 0) then begin
        valFullScale=fltarr(ctFullScale)
        for i=0l,ctFullScale-1 do begin
            temp=(where(ind eq indFullScale[i]))[0]
            valFullScale[i]=val[temp]
        endfor
        m=max(valFullScale)
        nbpnew=0
        repeat begin
            nbpnew=nbpnew+1
            two=m^(1./(8.*nbpnew))
        endrep until (two le 2)
    
        if nbpnew gt nbp then begin
            case nbpnew of
            1: begin
               tiff=byte(tiff)
               endcase
            2: begin
               tiff=uint(tiff)
               endcase
            else:begin
               tiff=ulong(tiff)
               endcase
            endcase
        endif
        tiff[indFullScale]=valFullScale
    endif
endif
    
FREE_LUN, lun

if bnorm then tiff*=realtime/livetime
    
o=[0,6,3,1,4,2,7,5,0]
return,rotate(temporary(tiff),o[orientation])
            
end;function ReadSMART
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReadCCD,path,file,format=format,info=info,static=static,error=error,xdinr=xdinr,nfiles=nfiles

; Reading a tiff: /unsigned keyword: 16bit
; Writing a tiff: /short keyword:16bit

error=1b
nfiles=1
if N_PARAMS() eq 2 then path2=path+file else path2=path
if not keyword_set(format) then begin
    PPos = RSTRPOS(path2, '.')
    if PPos eq -1 then format='.mccd' else format=strmid(path2,PPos)
endif
bstatic=keyword_set(static)

;CATCH, Error_status
;IF Error_status NE 0 THEN BEGIN
;    ;message,!ERR_STRING,/info
;    if n_elements(lun) ne 0 then free_lun,lun
;    
;    if bstatic then return,0 else return,PTR_NEW()
;ENDIF

sFormat=strlowcase(format)
if sFormat eq '.tiff' then sFormat='.tif'
if sFormat eq '.jpeg' then sFormat='.jpg'
if strmid(sFormat,0,4) eq '.mar' then sFormat='.mar345'
if strmid(sFormat,0,4) eq '.pck' then sFormat='.mar345'

bmultipage=0b
case sFormat of
'.tif':     begin
            tiff=ReadCCDtiff(path2)
            endcase
'.edf':        begin
            tiff=ReadEDF(path2,norm=info)
            endcase
'.raw':        begin
            tiff=ReadRAW(path2,info=info)
            endcase
'.spe':        begin
            tiff=ReadSPE(path2,xdinr=xdinr,nfiles=nfiles)
            bmultipage=nfiles gt 1
            endcase
'.dm3':        begin
            tiff=ReadDM3(path2,info=info)
            endcase
'.im':        begin
            tiff=ReadRigakuRaxis(path2)
            endcase
'.sif':        begin
            tiff=ReadSIF(path2,xdinr=xdinr,nfiles=nfiles)
            bmultipage=nfiles gt 1
            endcase
'.mccd':    begin
            tiff=ReadMarCCD(path2)
;            tiff=read_tiff(path2)
            endcase
'.jpg':        begin
            tiff=ReadCCDjpg(path2)
            endcase
'.bmp':        begin
            tiff=ReadCCDbmp(path2)
            endcase
'.mar345':    begin
            tiff=ReadMarIP345(path2)
            endcase
else:       begin; SMART frames '.xrd';'.001';...
            tiff=ReadSMART(path2)
            endcase
endcase

;tiff = not tiff
if bmultipage eq 0 then xdinr=-1l

if bstatic then begin
    if n_elements(tiff) le 4 then return,0
    error=0b
    return,tiff
endif else begin
    if n_elements(tiff) le 4 then return,ptr_new()
    ptiff=PTR_NEW(/allocate_heap)
    *ptiff=temporary(tiff)
    error=0b
    return,ptiff
endelse

; ----Orientation----

; A. Three Operations:
;tiff=reverse(tiff,1) ; => -X,+Y
;tiff=reverse(tiff,2) ;    => +X,-Y
;tiff=transpose(tiff) ;    => +Y,+X
; combinations: e.g. -Y,+X => tiff=transpose(tiff) and tiff=reverse(tiff,1)
;                or reverse(tiff,2) and tiff=transpose(tiff)
;                we choose the one with reverse first!

; B. Different formats:
; .tif file orientation values:
;0     Column 0 represents the left-hand side, and row 0 represents the bottom (same as 4)
;    => +X,-Y (default: convert to this) ->rotate(...,0)
;1     Column 0 represents the left-hand side, and row 0 represents the top.
;    => +X,+Y ->rotate(...,7)
;2     Column 0 represents the right-hand side, and row 0 represents the top.
;    => -X,+Y ->rotate(...,2)
;3     Column 0 represents the right-hand side, and row 0 represents the bottom.
;    => -X,-Y ->rotate(...,5)
;4     Column 0 represents the left-hand side, and row 0 represents the bottom (same as 0)
;    => +X,-Y ->rotate(...,0)
;5     Column 0 represents the top, and row 0 represents the left-hand side.
;    => +Y,+X ->rotate(...,3)
;6     Column 0 represents the top, and row 0 represents the right-hand side.
;    => +Y,-X ->rotate(...,4)
;7     Column 0 represents the bottom, and row 0 represents the right-hand side.
;    => -Y,-X ->rotate(...,1)
;8     Column 0 represents the bottom, and row 0 represents the left-hand side.
;    => -Y,+X ->rotate(...,6)

; SMART file orientation values:
;1    +X,+Y ->rotate(...,6)
;2    +X,-Y ->rotate(...,3)
;3    -X,+Y ->rotate(...,1)
;4    -X,-Y ->rotate(...,4)
;5    +Y,+X ->rotate(...,2)
;6    +Y,-X ->rotate(...,7)
;7    -Y,+X ->rotate(...,5)
;8    -Y,-X    (default: convert to this) ->rotate(...,0)

; IDL rotation function:
;0    +X,+Y
;1    -Y,+X    REV2+TRANS
;2    -X,-Y    REV1+REV2
;3    +Y,-X    REV1+TRANS
;4    +Y,+X    TRANS
;5    -X,+Y    REV1
;6    -Y,-X    REV1+REV2+TRANS
;7    +X,-Y    REV2

; C. XRDUA:
; We want an array with indices like this
; looking from behind the detector:
;  _______
; |0 1 2 3
; |4 5 5 7
; |....
;
; IDL flips this when plotting (bottom-up) so that the cursor
; position (x,y) is also the tiff-matrix col and row index.
; Remark: since in general we don't know the connection between the
; software determined orientation number and the physical orientation
; (the detector may be turned), we will choose the "default" marked
; number as "backview-orientation".

end;function ReadCCD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro rebinCCD,list,tiff,b

inverse=b lt 0
b=abs(fix(b))

s=*list.tiffs
if s[0] ne 2 then return

printw,list.OutID,'Original image: '
printw,list.OutID,'X size: '+string(s[1])
printw,list.OutID,'X size: '+string(s[2])
printw,list.OutID,'X Pixel size: '+string(list.scrx)
printw,list.OutID,'Y Pixel size: '+string(list.scry)
printw,list.OutID,'X center: '+string(list.center[0])
printw,list.OutID,'Y center: '+string(list.center[1])

if inverse then begin
    s1=s[1]*b
    s2=s[2]*b
    
    pkeep=tiff
    tiff=ptr_new(rebin(*tiff,s1,s2,/sample))
    ptr_free,pkeep
    
    list.scrx/=b
    list.scry/=b
    list.center=(list.center+0.5)*b-0.5
    list.GridOffset=(list.GridOffset+0.5)*b-0.5
    list.GridOrigin=(list.GridOrigin+0.5)*b-0.5
    list.GridSpPix=(list.GridSpPix+0.5)*b-0.5
    
    if ptr_valid(list.TieI) then (*list.TieI)=((*list.TieI)+0.5)*b-0.5
    if ptr_valid(list.TieO) then (*list.TieO)=((*list.TieO)+0.5)*b-0.5
    if ptr_valid(list.TiltI) then (*list.TiltI)=((*list.TiltI)+0.5)*b-0.5

endif else begin
    ;if (s[1] mod b) ne 0 or (s[2] mod b) ne 0 then return
    s1=s[1]/b
    s2=s[2]/b
    
    t1=rebin(b*indgen(s1)+b-1,s1,s2,/sample)
    t2=rebin(transpose(b*indgen(s2)+b-1),s1,s2,/sample)
    
    kernel=replicate(1,b,b)
    c=convol(float(*tiff),kernel)
    ptr_free,tiff
    tiff=ptr_new(c[t1,t2])
    
    list.scrx*=b
    list.scry*=b
    list.center=(list.center+0.5)/b-0.5
    list.GridOffset=(list.GridOffset+0.5)/b-0.5
    list.GridOrigin=(list.GridOrigin+0.5)/b-0.5
    list.GridSpPix=(list.GridSpPix+0.5)/b-0.5
    
    if ptr_valid(list.TieI) then (*list.TieI)=((*list.TieI)+0.5)/b-0.5
    if ptr_valid(list.TieO) then (*list.TieO)=((*list.TieO)+0.5)/b-0.5
    if ptr_valid(list.TiltI) then (*list.TiltI)=((*list.TiltI)+0.5)/b-0.5
endelse

printw,list.OutID,'Rebined image: '
printw,list.OutID,'X size: '+string(s1)
printw,list.OutID,'X size: '+string(s2)
printw,list.OutID,'X Pixel size: '+string(list.scrx)
printw,list.OutID,'Y Pixel size: '+string(list.scry)
printw,list.OutID,'X center: '+string(list.center[0])
printw,list.OutID,'Y center: '+string(list.center[1])

printw,list.OutID,'Calibration and spatial distortion tie points recalculated'
end;pro rebinCCD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WriteMarCCD,file,tiff
if openw_safe(lun, file) then return,0

type=size(tiff,/type)
bptr=type eq 10

offset=bytarr(4096)
writeu,lun,offset ;Make header
writeu,lun,uint((bptr?*tiff:tiff))
FREE_LUN, lun
return,1
end;function WriteMarCCD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WriteSMART,file,tiff,struc,path
if openw_safe(lun, file) then return,0

type=size(tiff,/type)
bptr=type eq 10

; ----Write header----
; 96 lines (15*512 bytes)
IntFormat='(I10)'
IntFormat2='(I9)'
IntFormat3='(I7)'
IntFormat4='(I12)'
AddInt=replicate(' ',62)
AddInt2=replicate(' ',52)
AddInt2b=replicate(' ',48)
AddInt3=replicate(' ',42)
AddInt3b=replicate(' ',36)
AddInt5=replicate(' ',22)
AddInt5b=replicate(' ',12)
AddInt6=AddInt5b
FltFormat='(f14.7)'
FltFormat2='(f10.6)'
AddFlt=replicate(' ',58)
AddFlt2=replicate(' ',44)
AddFlt3=replicate(' ',30)
AddFlt4=replicate(' ',16)
AddFlt5=replicate(' ',2)
AddUnknown=replicate(' ',65)
Add72=replicate(' ',72)
AddFill=[26B,4B,replicate(46B,78)]

; Save as 8-bit image with overflow table if necessary
; When to mush overflow: save as 16-bit or 32-bit
FullScale=255B
indOver=where((bptr?*tiff:tiff) gt FullScale,nOverfl)
bpp=1
if (nOverfl gt 4096) then begin
    FullScale=65535U
    indOver=where((bptr?*tiff:tiff) gt FullScale,nOverfl)
    bpp=2
    if (nOverfl gt 4096) then begin
        FullScale=2UL^32-1
        indOver=where((bptr?*tiff:tiff) gt FullScale,nOverfl)
        bpp=4
    endif
endif
case bpp of
    1:tiffo=byte((bptr?*tiff:tiff)<FullScale)
    2:tiffo=uint((bptr?*tiff:tiff)<FullScale)
    4:tiffo=ulong((bptr?*tiff:tiff)<FullScale)
endcase

; Current time
CALDAT, SYSTIME(/julian), Month, Day, Year, Hour, Minute, Second
created=string(Month,format='(I2,"/")')+string(Day,format='(I2,"/")')+strmid(string(Year,format='(I4,"   ")'),2,5)+$
    string(Hour,format='(I2,":")')+string(Minute,format='(I2,":")')+string(Second,format='(I2)')

; Dim
s=size((bptr?*tiff:tiff),/dimensions)
ncol=s[0]
nrow=s[1]
mx=ncol-1
my=nrow-1
mix=0
miy=0
;Min X, Min (X+Y), Min Y, Max (X-Y), Max X, Max (X+Y), Max Y, Max (Y-X)
octmask=[mix,mix+miy,miy,mx-miy,mx,mx+my,my,my-mix]

writeu,lun,'FORMAT :',string(86,FORMAT = IntFormat),AddInt
writeu,lun,'VERSION:',string(9,FORMAT = IntFormat),AddInt
writeu,lun,'HDRBLKS:',string(15,FORMAT = IntFormat),AddInt
writeu,lun,'TYPE   :','Add frame',replicate(' ',63)
writeu,lun,'SITE   :','UNKNOWN',AddUnknown
writeu,lun,'MODEL  :','P4',replicate(' ',70)
writeu,lun,'USER   :','UNKNOWN',AddUnknown
writeu,lun,'SAMPLE :','UNKNOWN',AddUnknown
writeu,lun,'SETNAME:','UNKNOWN',AddUnknown
writeu,lun,'RUN    :',string(0,FORMAT = IntFormat),AddInt
writeu,lun,'SAMPNUM:',string(1,FORMAT = IntFormat),AddInt
if strlen(path) gt 72 then writeu,lun,'TITLE  :',strmid(path,0,69)+'...' $
else writeu,lun,'TITLE  :',path,replicate(' ',72-strlen(path))
writeu,lun,'TITLE  :','?',replicate(' ',71)
writeu,lun,'TITLE  :',Add72
writeu,lun,'TITLE  :',Add72
writeu,lun,'TITLE  :',Add72
writeu,lun,'TITLE  :',Add72
writeu,lun,'TITLE  :',Add72
writeu,lun,'TITLE  :',Add72
writeu,lun,'NCOUNTS:',string(total((bptr?*tiff:tiff)),FORMAT = IntFormat),AddInt
writeu,lun,'NOVERFL:',string(ceil(nOverfl/512.),FORMAT = IntFormat),AddInt
writeu,lun,'MINIMUM:',string(min((bptr?*tiff:tiff)),FORMAT = IntFormat),AddInt
writeu,lun,'MAXIMUM:',string(max((bptr?*tiff:tiff),indmax),FORMAT = IntFormat),AddInt
writeu,lun,'NONTIME:',string(0,FORMAT = IntFormat),AddInt
writeu,lun,'NLATE  :',string(0,FORMAT = IntFormat),AddInt
writeu,lun,'FILENAM:',struc.file,replicate(' ',72-strlen(struc.file))
writeu,lun,'CREATED:',created,replicate(' ',72-strlen(created))
writeu,lun,'CUMULAT:',string(0.,FORMAT = FltFormat),AddFlt
writeu,lun,'ELAPSDR:',string(0.,FORMAT = FltFormat),AddFlt
writeu,lun,'ELAPSDA:',string(0.,FORMAT = FltFormat),AddFlt
writeu,lun,'OSCILLA:',string(0,FORMAT = IntFormat),AddInt
writeu,lun,'NSTEPS :',string(1,FORMAT = IntFormat),AddInt
writeu,lun,'RANGE  :',string(0.25,FORMAT = FltFormat),AddFlt
writeu,lun,'START  :',string(0.,FORMAT = FltFormat),AddFlt
writeu,lun,'INCREME:',string(0.25,FORMAT = FltFormat),AddFlt
writeu,lun,'NUMBER :',string(0,FORMAT = IntFormat),AddInt
writeu,lun,'NFRAMES:',string(1,FORMAT = IntFormat),AddInt
writeu,lun,'ANGLES :',string([0.,0,0,0],FORMAT = FltFormat),AddFlt4
writeu,lun,'NOVER64:',string(intarr(3),FORMAT = IntFormat4),AddInt3b
writeu,lun,'NPIXELB:',string(bpp,FORMAT = IntFormat),AddInt
writeu,lun,'NROWS  :',string(nrow,FORMAT = IntFormat),AddInt
writeu,lun,'NCOLS  :',string(ncol,FORMAT = IntFormat),AddInt
writeu,lun,'WORDORD:',string(0,FORMAT = IntFormat),AddInt
writeu,lun,'LONGORD:',string(0,FORMAT = IntFormat),AddInt
writeu,lun,'TARGET :','UNKNOWN',AddUnknown
writeu,lun,'SOURCEK:',string(0.,FORMAT = FltFormat),AddFlt
writeu,lun,'SOURCEM:',string(0.,FORMAT = FltFormat),AddFlt
writeu,lun,'FILTER :','UNKNOWN',AddUnknown
writeu,lun,'CELL   :',string([1.,1,1,90,90],FORMAT = FltFormat),AddFlt5
writeu,lun,'CELL   :',string(90.,FORMAT = FltFormat),AddFlt
writeu,lun,'MATRIX :',string([1.,0,0,0,1],FORMAT = FltFormat),AddFlt5
writeu,lun,'MATRIX :',string([0.,0,0,1],FORMAT = FltFormat),AddFlt4
writeu,lun,'LOWTEMP:',string(intarr(3),FORMAT = IntFormat4),AddInt3b
writeu,lun,'ZOOM   :',string([0.5,0.5,1],FORMAT = FltFormat),AddFlt3
writeu,lun,'CENTER :',string(struc.center,FORMAT = FltFormat),AddFlt2
writeu,lun,'DISTANC:',string(float(100*struc.dist),FORMAT = FltFormat),AddFlt
writeu,lun,'TRAILER:',string(-1,FORMAT = IntFormat),AddInt
writeu,lun,'COMPRES:','NONE',replicate(' ',68)
writeu,lun,'LINEAR :',string([1.,0],FORMAT = FltFormat),AddFlt2
writeu,lun,'PHD    :',string([0.,0],FORMAT = FltFormat),AddFlt2
writeu,lun,'PREAMP :',string(1.,FORMAT = FltFormat),AddFlt
writeu,lun,'CORRECT:','LINEAR',replicate(' ',66)
writeu,lun,'WARPFIL:','LINEAR',replicate(' ',66)
writeu,lun,'WAVELEN:',string(float([struc.lambda,struc.lambda,struc.lambda]),FORMAT = FltFormat),AddFlt3
writeu,lun,'MAXXY  :',string(float([indmax mod ncol,indmax/ncol]),FORMAT = FltFormat),AddFlt2
writeu,lun,'AXIS   :',string(2,FORMAT = IntFormat),AddInt
writeu,lun,'ENDING :',string([0.,0.25,0,0],FORMAT = FltFormat),AddFlt4
writeu,lun,'DETPAR :',string([0.,0,0,0,0],FORMAT = FltFormat),AddFlt5
writeu,lun,'DETPAR :',string(0.,FORMAT = FltFormat),AddFlt
writeu,lun,'LUT    :','$DEFAULT',replicate(' ',64)
writeu,lun,'DISPLIM:',string([0.,0],FORMAT = FltFormat),AddFlt2
writeu,lun,'PROGRAM:','xrdua',replicate(' ',67)
writeu,lun,'ROTATE :',string(0,FORMAT = IntFormat),AddInt
writeu,lun,'BITMASK:','$NULL',replicate(' ',67)
writeu,lun,'OCTMASK:',string(octmask[0:5],FORMAT = IntFormat4)
writeu,lun,'OCTMASK:',string(octmask[6:7],FORMAT = IntFormat4),AddInt2b
writeu,lun,'ESDCELL:',string([0.02,0.02,0.02,0.02,0.02],FORMAT = FltFormat),AddFlt5
writeu,lun,'ESDCELL:',string(0.02,FORMAT = FltFormat),AddFlt
writeu,lun,'DETTYPE:','CCD-UNKNOWN',string(float([1./struc.scrx/100.,0.0,0,0,0]),FORMAT = FltFormat2),replicate(' ',11);pixels/cm, cm to grid(distance between front face and phosphor),...
writeu,lun,'NEXP   :',string([1,0,0,8,0],FORMAT = IntFormat4),AddInt5b ; 1=single/2=correlated sum, per-exposure bias level, baseline offset, orientation (1-8), 1=overscan/0=no overscan data
writeu,lun,'CCDPARM:',string(fltarr(5),FORMAT = FltFormat),AddFlt5 ;readnoise, e/ADU, e/photon, bias, full scale
writeu,lun,'CHEM   :','?',replicate(' ',71)
writeu,lun,'MORPH  :','?',replicate(' ',71)
writeu,lun,'CCOLOR :','?',replicate(' ',71)
writeu,lun,'CSIZE  :','?|?|?|?|?',replicate(' ',63)
writeu,lun,'DNSMET :','?',replicate(' ',71)
writeu,lun,'DARK   :',struc.darkfile,replicate(' ',72-strlen(struc.darkfile))
writeu,lun,'AUTORNG:',string(fltarr(5),FORMAT = FltFormat),AddFlt5 ; gain, high-speed time, scale, offset, full scale
writeu,lun,'ZEROADJ:',string(fltarr(4),FORMAT = FltFormat),AddFlt4
writeu,lun,'XTRANS :',string(fltarr(3),FORMAT = FltFormat),AddFlt3
writeu,lun,'HKL&XY :',string(fltarr(5),FORMAT = FltFormat),AddFlt5
writeu,lun,'AXES2  :',string(fltarr(4),FORMAT = FltFormat),AddFlt4
writeu,lun,'ENDING2:',string(fltarr(4),FORMAT = FltFormat),AddFlt4
writeu,lun,AddFill,AddFill,AddFill

; ----write data----
writeu,lun,tiffo

; ----write overflow table----
if nOverFl ne 0 then begin
    for i=0l,nOverFl-1 do begin
        writeu,lun,string((bptr?*tiff:tiff)[indOver[i]],FORMAT = IntFormat2),string(indOver[i],FORMAT = IntFormat3)
    endfor
    padd=-nOverFl*16
    while (padd lt 0) do padd+=512
    if padd ne 0 then writeu,lun,bytarr(padd)
endif

FREE_LUN, lun
return,1
end;function WriteSMART
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WriteCCD,tiff,path,file,struc=struc,format=format

type=size(tiff,/type)
bptr=type eq 10

if N_PARAMS() eq 3 then path2=path+file else path2=path
if not keyword_set(format) then begin
    PPos = RSTRPOS(path2, '.')
    if PPos eq -1 then format='.mccd' else format=strmid(path2,PPos)
endif

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    ;message,!ERR_STRING,/info
    if n_elements(lun) ne 0 then free_lun,lun
    
    return,0
ENDIF

tmp=strlowcase(format)
if tmp eq '.tiff' then tmp='.tif'
if tmp eq '.jpeg' then tmp='.jpg'
case tmp of
'.tif':    begin
        write_tiff_safe,path2,(bptr?*tiff:tiff),COMPRESSION=1,ORIENTATION=0
        endcase
'.mccd':begin
           if WriteMarCCD(path2,tiff) eq 0 then return,0b
        endcase
'.edf':    begin
        if WriteEDF(path2,tiff) eq 0 then return,0b
        endcase
else:   begin
        if tmp eq '.raw' or tmp eq '.dm3' or tmp eq '.jpg' or tmp eq '.bmp' or tmp eq '.spe' or tmp eq '.im' or tmp eq '.sif' then begin
            path2+='.tif'    
            write_tiff_safe,path2,(bptr?*tiff:tiff),COMPRESSION=0,ORIENTATION=0
        endif else begin
            if WriteSMART(path2,tiff,struc,path) eq 0 then return,0b
        endelse
        endcase
endcase

return,1B

end;function WriteCCD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CopyFields,path,fields,OutID
if openr_safe(lun, path) then return,''
n=n_elements(fields)

line=''
buffer=''
for i=0l,n-1 do begin
    block=fields[i]
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin

       repeat begin
         buffer=[buffer,line]
         readf,lun,line
       endrep until (strmid(line,0,1) eq '$') or eof(lun)
       if (strmid(line,0,1) ne '$') then buffer=[buffer,line]

       printw,OutID,block+' Copied'
    endif
endfor

free_lun,lun
return,buffer[1<(n_elements(buffer)-1):*]
end;function CopyFields
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SaveMask,list,prompt,oD=oD

; ----File for output----
if prompt then begin
    path=SelFile(list.mskpath,list.mskfile,'*.msk','Save Mask To....')
endif else path=list.mskpath+list.mskfile

if path eq '' then return

error=CutPath(path,path=mskpath,file=file,ext=ext)
list.mskpath=mskpath
list.mskfile=file+ext
result=file_search(path,count=ct)
if ct eq 0 then exist=0B else exist=1B

; ----Define fields----
; 1D and 2D fields
fields0=['$Version:','$SOURCE ('+msymbols('angstroms')+'):','$PIXSIZE (m):','$DET-SAM (m):','$TILTA:','$AZIOFF:','$IntCor2Dinfo:']
; 1D fields
fields1=['$FIT_ROI:','$BACKoD:','$FIL:','$LABELS:','$MASK1:','$FITMODEL:']
; 2D fields
fields2=['$CENTER:','$MASK:','$BACK:','$SEARCH:','$SPADIS:','$TILTP:','$CORRECTIONS:','$FF:','$MASKOFF:',$
         '$SATURATION:','$ZINGERS:','$AzInt:','$ReadCCDInfo:']

; ----Copy fields----
if exist then begin
   if keyword_set(oD) then buffer=CopyFields(path,fields2,list.OutID) $
   else buffer=CopyFields(path,fields1,list.OutID)
endif

; ----Open and print copied fields----
if openw_safe(lun, path, width=1000) then begin
    printw,list.OutID,'Mask NOT saved to: '+path
    return
endif
if exist then begin
    if (buffer[0] ne '') or (n_elements(buffer) ne 1) then $
    for i=0l,n_elements(buffer)-1 do $
        printf,lun,buffer[i]
endif

; ----Print common fields----
printf,lun,'$Version:'
printf,lun,version()
printf,lun,'$SOURCE ('+msymbols('angstroms')+'):'
printf,lun,list.lambda,format='(E)'
printf,lun,'$PIXSIZE (m):'
printf,lun,list.scrx,list.scry,format='(2E)'
printf,lun,'$DET-SAM (m):'
printf,lun,list.dist,format='(E)'
if list.dist gt 100 then begin
    writeu,lun,list.dist
    printf,lun,''
endif
printf,lun,'$TILTA:'
printf,lun,[list.a,list.b],format='(2E)'
printf,lun,'$AZIOFF:'
printf,lun,[list.phipixphi0],format='(E)'
printf,lun,'$IntCor2Dinfo:'
printf,lun,[list.IntCor2Dinfo.fftype,list.IntCor2Dinfo.Atype,list.IntCor2Dinfo.Ptype,list.IntCor2Dinfo.Ltype,list.IntCor2Dinfo.n,list.IntCor2Dinfo.s]
printf,lun,[list.IntCor2Dinfo.Pff,list.IntCor2Dinfo.Plindeg,list.IntCor2Dinfo.mmult,list.IntCor2Dinfo.dmono,$
            list.IntCor2Dinfo.flatnormphi,list.IntCor2Dinfo.flatnormpolar,list.IntCor2Dinfo.flatD,list.IntCor2Dinfo.muL],format='(8E)'
                
if keyword_set(oD) then begin
    printf,lun,'$FIT_ROI:'
    printf,lun,list.backranind
    printf,lun,list.backranindD,list.nxval
    printf,lun,list.ROIAzimuth
    printf,lun,'$BACKoD:'
    printf,lun,[list.backmod,(*(*list.ptrModel).ptrFit).backparamIO]
    printf,lun,'$FIL:'
    printf,lun,list.vin,list.win,list.r1,list.fitmethod,list.filorder
    printf,lun,'$LABELS:'
    printf,lun,list.alabels
    if (list.alabels ne 0) then begin
        printf,lun,list.xtypeorig
        for i=0l,list.alabels-1 do printf,lun,(*list.labels)[i],(*list.xlabels)[i,list.xtypeorig],(*list.ylabels)[i]
    endif
    printf,lun,'$MASK1:'
    printf,lun,list.masknr
    if list.masknr ne 0 then begin
       m=stringr(strlen(stringr(fix(max((*list.mask)[1,*]))))+1)
       printf,lun,format='(I2,I'+m+',I2,F10.4,F10.4)',(*list.mask)
    endif
    printf,lun,'$FITMODEL:'
    SaveFitModel,lun,list.ptrModel
endif else begin
    printf,lun,'$CENTER:'
    printf,lun,list.center,format='(E,E)'
    ind=where(abs(list.center) gt 1000000,ct)
    if ct ne 0 then begin
        writeu,lun,list.center
        printf,lun,''
    endif
    printf,lun,'$MASK:'
    printf,lun,list.masknr
    if list.masknr ne 0 then begin
       m=stringr(strlen(stringr(fix(max((*list.mask)[1,*]))))+1)
       printf,lun,format='(I2,I'+m+',I2,F10.4,F10.4,F10.4,F10.4)',(*list.mask)
    endif
    printf,lun,'$BACK:'
    printf,lun,list.wb,list.niter,list.bcko,list.bckw,list.bckm,list.thres
    printf,lun,list.darkpath
    printf,lun,list.darkfile
    printf,lun,'$SEARCH:'
    printf,lun,list.v,list.w,list.r,list.fit
    printf,lun,'$SPADIS:'
    printf,lun,list.TieGrid
    case list.TieGrid of
    0:    begin
        printf,lun,list.nrTieRing
        printf,lun,(*list.nrTie)
        if ptr_valid(list.TieI) then printf,lun,(*list.TieI)
        if ptr_valid(list.TieRingD) then printf,lun,(*list.TieRingD)
        endcase
    1:    begin
        printf,lun,list.nrTieRing
        printf,lun,(*list.nrTie)
        if ptr_valid(list.TieI) then printf,lun,(*list.TieI)
        if ptr_valid(list.TieO) then printf,lun,(*list.TieO)
        printf,lun,list.GridDim
        printf,lun,list.GridSpReal
        printf,lun,list.GridOffset
        printf,lun,list.GridOrigin
        printf,lun,list.GridTilt
        endcase
       2:    begin
           printf,lun,list.TiePath
           printf,lun,list.TieFile
           endcase
       endcase
       printf,lun,list.bFastSpline,list.vFlipSpline
    printf,lun,'$TILTP:'
    printf,lun,list.nrTiltRing
    printf,lun,(*list.nrTilt)
    if ptr_valid(list.TiltRingD) then printf,lun,(*list.TiltRingD)
    if ptr_valid(list.TiltI) then printf,lun,(*list.TiltI)
    printf,lun,'$CORRECTIONS:'
    n=n_elements(list.CorDone)
    printf,lun,n
    for i=0l,n-1 do begin
        printf,lun,list.SeqTags[i]
        printf,lun,(*list.CorSeq).(i)
    endfor
    printf,lun,list.CorMethod
    printf,lun,'$FF:'
    printf,lun,list.ffpath
    printf,lun,list.fffile
    printf,lun,'$MASKOFF:'
    printf,lun,list.maskoffpath
    printf,lun,list.maskofffile
    printf,lun,list.maskoffthres
    printf,lun,list.maskoffdyn
    printf,lun,list.maskoffdynthres
    printf,lun,'$SATURATION:'
    printf,lun,list.SaturMax
    printf,lun,list.SaturAspect
    printf,lun,list.SaturNMax
    printf,lun,'$ZINGERS:'
    printf,lun,list.ZingerWidth
    printf,lun,list.ZingerThreshold
    printf,lun,'$AzInt:'
    printf,lun,[list.AzIntInfo.oD,list.AzIntInfo.soD,list.AzIntInfo.Fill,list.AzIntInfo.type,$
        list.AzIntInfo.bool,list.AzIntInfo.SaveRingPercentage,list.AzIntInfo.Batch,list.AzIntInfo.tiff,list.AzIntInfo.origformat]
    printf,lun,[list.AzIntInfo.MinPix, list.AzIntInfo.BWP,list.AzIntInfo.BWT,list.AzIntInfo.BWD,list.AzIntInfo.BWQ,$
        list.AzIntInfo.AzWidth],format='(6E)'
    printf,lun,list.AzIntInfo.outbin
    printf,lun,list.AzIntInfo.median
    printf,lun,list.AzIntInfo.stD
    printf,lun,list.AzIntInfo.nsectors
    printf,lun,[list.AzIntInfo.aoD,list.AzIntInfo.loD]
    printf,lun,list.AzIntInfo.percentile
    printf,lun,list.AzIntInfo.updatenorm,list.AzIntInfo.calcerror
    printf,lun,'$ReadCCDInfo:'
    printf,lun,list.readinfo.showheader
    printf,lun,list.readinfo.bnorm
    if strlen(list.readinfo.int) eq 0 then printf,lun,' ' $
    else printf,lun,list.readinfo.int
    if strlen(list.readinfo.gain) eq 0 then printf,lun,' ' $
    else printf,lun,list.readinfo.gain
    if strlen(list.readinfo.offset) eq 0 then printf,lun,' ' $
    else printf,lun,list.readinfo.offset
    printf,lun,list.readinfo.binxs
    printf,lun,list.readinfo.binys
    printf,lun,list.readinfo.binlendian
    printf,lun,list.readinfo.bintype
    printf,lun,list.readinfo.binnxpad3
endelse

free_lun,lun
printw,list.OutID,'Mask saved to: '+path
end;pro SaveMask
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AdaptDirectory,dir,refdir,file=file
; For example the directory of a background image is
;     dir='C:\data\results\sample1\results\'
; 
; The directory of the mask file is
;     refdir='/home/data/user/data/results/calib/results/'
;
; Return '/home/data/user/data/results/sample1/results/'

; Don't change dir when it exists
if ~keyword_set(file) then file=''
if file_test(dir+file) then return,dir

; Prepare search/match
dirsep=strsplit(dir,'/\',/extract,count=n)
refsep=strsplit(refdir,'/\',/extract,count=nref)
sep=(total(byte(refdir) eq byte('/')) gt total(byte(refdir) eq byte('\')))?'/':'\'
if strmid(refdir,0,1) eq '/' then init='/' else init=''

; Find the correct subdirectory of refdir
bfound=0b
for i=n-1,0l,-1l do begin
    ind=where(refsep eq dirsep[i],ct)
    for j=ct-1,0l,-1l do begin
        iref=ind[j]
        retdir=refsep[0:iref]
        if i lt n-1 then retdir=[retdir,dirsep[i+1:*]]
        retdir=init+strjoin(retdir,sep)+sep

        bfound=file_test(retdir+file)
        if bfound then break
    endfor
    if bfound then break
endfor

; Subdirectory not found
if ~bfound then begin
    if file eq '' then return,dir
    ; Maybe refdir itself?
    if file_test(refdir+file) then return,refdir
    return,dir
endif

; Subdirectory found
return,retdir
end;function AdaptDirectory
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro LoadMask,list,prompt,id,flags=flags,top=top,$
          flagi=flagi,error=error,oD=oD,otD=otD

flags=keyword_set(flags)
error=0B
; flags=0 -> 'Load data'
; flags=1 -> 'Add data': Pop-up Window
; flagi -> 'Fields not to use'
; flagi not set -> 'Use all fields'

; ----Select File----
path=''
if prompt then $
path=DIALOG_PICKFILE(path=list.mskpath,file=list.mskfile,filter='*.msk',title='Load Mask...') $
else path=list.mskpath+list.mskfile
if openr_safe(lun, path) then begin
    printw,id,path+' not found'
    error=1B
    return
endif
errorr=CutPath(path,path=mskpath,file=file,ext=ext)
list.mskpath=mskpath
list.mskfile=file+ext

; ----Select which things should be loaded----
; 1D and 2D fields
fields0=['$Version:','$SOURCE ('+msymbols('angstroms')+'):','$PIXSIZE (m):','$DET-SAM (m):','$TILTA:','$AZIOFF:','$IntCor2Dinfo:']
n0=n_elements(fields0)
otDfields=indgen(n0)
; 1D fields
fields1=['$FIT_ROI:','$BACKoD:','$FIL:','$LABELS:','$MASK1:','$FITMODEL:']
n1=n_elements(fields1)
oDfields=indgen(n1)+n0
; 2D fields
fields2=['$CENTER:','$MASK:','$BACK:','$SEARCH:',$
        '$SPADIS:','$TILTP:','$CORRECTIONS:','$FF:','$MASKOFF:','$SATURATION:','$ZINGERS:','$AzInt:','$ReadCCDInfo:']
n2=n_elements(fields2)
tDfields=indgen(n2)+n0+n1


fields=[fields0,fields1,fields2]
nfields=n_elements(fields)
; Create heap variable with pointer to it
; The heap variable is global, the pointer local
; Load Everything
PTR_flagi = PTR_NEW(bytarr(nfields)+1)
; Don't load fields supplied in flagi to 0
if keyword_set(flagi) then begin
    for i=0l,n_elements(flagi)-1 do begin
       ind=where(fields eq flagi[i],count)
       if count ne 0 then (*PTR_flagi)[ind[0]]=0
    endfor
endif

; Set fields to zero
if not keyword_set(otD) then $
if keyword_set(oD) then (*PTR_flagi)[tDfields]=0 $
else (*PTR_flagi)[oDfields]=0

; ----Pop-up Window for changing selection of fields----
if flags then begin
    if keyword_set(top) then WIDGET_CONTROL, top, sensitive=0 else top=0L

    base=widget_base(/column,title='Select fields to load',xoffset=200,yoffset=200,$
    uvalue={ID1:top,select:PTR_flagi})
    ; PTR_flagi is local in scope so we make an other pointer to the same heap variable
    nr=5
    nbase=ceil(float(nfields)/nr)
    base6=lonarr(nbase)
    for i=0l,nbase-1 do base6[i]=widget_base(base,/row,/nonexclusive)
    rbutton = lonarr(nfields+1)
    j=0
    for i=0l,nfields-1 do begin
       val=stringr(i)
       rbutton[i]=widget_button(base6[j],value=fields[i],uname=val,uvalue=i)
       if ((i+1) mod nr) eq 0 then j=j+1
    endfor
    rbutton[i]=widget_button(base6[nbase-1],value='Select All')
    releasebutton=widget_button(base,value='OK',uname='OK')

    WIDGET_CONTROL, base, /REALIZE
    for i=0l,nfields-1 do widget_control,rbutton[i],set_button=(*PTR_flagi)[i]
    if not keyword_set(otD) then $
    if keyword_set(oD) then $
       for i=0l,n2-1 do  widget_control,rbutton[tDfields[i]],sensitive=0 $
    else for i=0l,n1-1 do  widget_control,rbutton[oDfields[i]],sensitive=0
    Xmanager,'LoadMask',base, event_handler='SelMask_event',$
        cleanup='CleanUp_SelMask',GROUP_LEADER=top

    ; Set fields to zero
    if not keyword_set(otD) then $
    if keyword_set(oD) then (*PTR_flagi)[tDfields]=0 $
    else (*PTR_flagi)[oDfields]=0
endif

flagi=*PTR_flagi
PTR_FREE,PTR_flagi; destroy the heap variable
PTR_flagi=0; destroy the pointer


; ----Start reading file----
nError=0

versionmask='0.0.0.0'
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$Version:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,line
       versionmask=line
       if line ne version() then begin
           tmp1=long(strsplit(line,'.',/extract))
           tmp2=long(strsplit(version(),'.',/extract))
           ind=where(tmp1 ne tmp2)
           printw,id,'Made with '+((tmp1[ind[0]] gt tmp2[ind[0]])?'newer':'previous')+' version: '+line
       endif
    endif else printw,id,'Made with previous version.'
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$CENTER:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    center=list.center
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,center
       list.center=center
       
       ind=where(abs(center) gt 1000000,ct)
       if ct ne 0 then begin
           readu,lun,center
           list.center=center
       endif
       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$MASK:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    masknr=0
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,masknr
       if flags then begin ;Add data
         if masknr ne 0 then begin
          mask=fltarr(7,masknr)
          readf,lun,mask
          if list.masknr ne 0 then (*list.mask)[2,*]=1
          if list.masknr eq 0 then L=0 else L=max((*list.mask)[1,0:(list.masknr-1)])+1
          mask[1,*]=mask[1,*]+L
          if list.masknr eq 0 then begin
              list.mask=PTR_NEW(mask)
              list.smask=PTR_NEW(replicate(1b,masknr))
          endif else begin
              (*list.mask)=[[(*list.mask)],[mask]]
              (*list.smask)=[(*list.smask),replicate(1b,masknr)]
          endelse
          list.masknr+=masknr
         endif
       endif else begin ;Replace data
            ptr_free,list.mask
            bsmask=(flagi[where(fields eq '$MASK1:')])[0] eq 0
            if bsmask then ptr_free,list.smask
            list.masknr=masknr
         if masknr ne 0 then begin
          mask=fltarr(7,masknr)
          readf,lun,mask
          list.mask=PTR_NEW(mask)
          if bsmask then list.smask=PTR_NEW(replicate(1b,masknr))
         endif
       endelse
       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$SOURCE ('+msymbols('angstroms')+'):'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    lambda=0.
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,lambda
       list.lambda=lambda
       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$PIXSIZE (m):'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    scr=fltarr(2)
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,scr
       list.scrx=scr[0]
       list.scry=scr[1]
       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$DET-SAM (m):'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    dist=list.dist
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,dist
       list.dist=dist
       
       if abs(dist) gt 100 then begin
           readu,lun,dist
           list.dist=dist
       endif 
       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$FITMODEL:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       LoadFitModel,lun,list.ptrModel,versionmask,add=flags
       printw,id,block+' loaded'
    endif else begin
        *(*(*list.ptrmodel).ptrFit).sourceemission=float([list.lambda,1,1])
    endelse
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$BACK:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
         in=fltarr(7)
         readf,lun,in
         list.wb=in[0]
         list.niter=in[1]
         list.bcko=in[2]
         list.bckw=in[3]
         list.bckm=in[4]
         list.thres=in[5:6]
         in=''
         readf,lun,in
         list.darkpath=in
         readf,lun,in
         list.darkfile=in
         list.darkpath=AdaptDirectory(list.darkpath,mskpath,file=list.darkfile)
         ptr_free,list.tiffb
         printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$SEARCH:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    in=intarr(4)
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,in
       list.v=in[0]
       list.w=in[1]
       list.r=in[2]
       list.fit=in[3]
       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$SPADIS:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    in=0
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,in
       TieGrid=in

       if TieGrid eq 2 then begin
               list.TieGrid=TieGrid
               readf,lun,line
               list.TiePath=line
               readf,lun,line
               list.TieFile=line
            list.TiePath=AdaptDirectory(list.TiePath,mskpath,file=list.TieFile)
       endif else begin
           readf,lun,in
           nrTieRing=in
           nrTie=lonarr(nrTieRing+1)
           readf,lun,nrTie

             if flags then begin ;Add data

              if list.TieGrid ne TieGrid then begin
                  printw,id,'Not the same Tie settings! (TieGrid= 0 or 1)'
              endif else begin
                    if nrTie[0] ne 0 then begin
                  ; Point to the next ring if one not finished
                  if (*list.nrTie)[list.nrTieRing] ne 0 then begin
                      list.nrTieRing++
                      (*list.nrTie)=[(*list.nrTie),nrTie]
                  endif else begin
                      if list.nrTieRing eq 0 then (*list.nrTie)=nrTie else $
                      (*list.nrTie)=[(*list.nrTie)[0:list.nrTieRing-1],nrTie]
                  endelse

                  n=total(nrTie)
                  in=fltarr(2,n)
                  readf,lun,in
                  if ptr_valid(list.TieI) then (*list.TieI)=[[(*list.TieI)],[in]] else list.TieI=ptr_new(in)

                  if TieGrid then begin
                      readf,lun,in
                      if ptr_valid(list.TieO) then (*list.TieO)=[[(*list.TieO)],[in]] else list.TieO=ptr_new(in)
                      in=lonarr(2)
                      readf,lun,in
                      list.GridDim=in
                      in=fltarr(2)
                      readf,lun,in
                      list.GridSpReal=in
                      readf,lun,in
                      list.GridOffset=in
                      if CompareVersion(versionmask,'5.3.1.1') ne -1 then begin
                          readf,lun,in
                        list.GridOrigin=in
                        in=0.
                        readf,lun,in
                        list.GridTilt=in
                      endif else begin
                          list.GridOrigin=0
                          list.GridTilt=0
                      endelse

                       tiffs=*list.tiffs
                       list.GridSpPix=[(tiffs[1]-1-2.*list.GridOffset[0])/(list.GridDim[0]-1.),$
                      (tiffs[2]-1-2.*list.GridOffset[1])/(list.GridDim[1]-1.)]
                      ; The following must not be done here!
                      ;list.scrx=list.GridSpReal[0]/list.GridSpPix[0]/1000
                      ;list.scry=list.GridSpReal[1]/list.GridSpPix[1]/1000
                  endif else begin
                      in=fltarr(nrTieRing+(nrTie[nrTieRing]<1))
                      readf,lun,in
                      if ptr_valid(list.TieRingD) then (*list.TieRingD)=[[(*list.TieRingD)],in] else list.TieRingD=ptr_new(in)
                  endelse

                  list.nrTieRing+=nrTieRing
                  endif
              endelse
             endif else begin ;Load data, not add
              list.TieGrid=TieGrid
              (*list.nrTie)=nrTie
              if nrTie[0] ne 0 then begin
              n=total(nrTie)
              in=fltarr(2,n)
              readf,lun,in
              ptr_free,list.TieI
              list.TieI=ptr_new(in)

              if list.TieGrid then begin
                  readf,lun,in
                  ptr_free,list.TieO
                    list.TieO=ptr_new(in)
                    if (list.listtype ne 40) and (list.listtype ne 41) then begin
                      in=lonarr(2)
                      readf,lun,in
                      list.GridDim=in
                      in=fltarr(2)
                      readf,lun,in
                      list.GridSpReal=in
                      readf,lun,in
                      list.GridOffset=in
                      if CompareVersion(versionmask,'5.3.1.1') ne -1 then begin
                          readf,lun,in
                        list.GridOrigin=in
                        in=0.
                        readf,lun,in
                        list.GridTilt=in
                      endif else begin
                          list.GridOrigin=0
                          list.GridTilt=0
                      endelse

                         tiffs=*list.tiffs
                       list.GridSpPix=[(tiffs[1]-1-2.*list.GridOffset[0])/(list.GridDim[0]-1.),$
                      (tiffs[2]-1-2.*list.GridOffset[1])/(list.GridDim[1]-1.)]
                      ; The following must not be done here!
                      ;list.scrx=list.GridSpReal[0]/list.GridSpPix[0]/1000
                      ;list.scry=list.GridSpReal[1]/list.GridSpPix[1]/1000
                  endif else begin
                      in=lonarr(2)
                      readf,lun,in
                      in=fltarr(2)
                      readf,lun,in
                      readf,lun,in
                  endelse
              endif else begin
                  in=fltarr(nrTieRing+(nrTie[nrTieRing]<1))
                  readf,lun,in
                  ptr_free,list.TieRingD
                  list.TieRingD=ptr_new(in)
              endelse

              endif
              list.nrTieRing=nrTieRing
             endelse
       endelse

        CATCH, Error_status ;Because of old msk versions
        IF Error_status eq 0 THEN BEGIN
           line=[1b,1b]
           readf,lun,line
           list.bFastSpline=line[0]
           list.vFlipSpline=line[1]
        ENDIF

       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$TILTA:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    in=fltarr(2)
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,in
       list.a=in[0]
       list.b=in[1]
       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$AZIOFF:'
if (flagi[where(fields eq block)])[0] then begin
    list.phipixphi0=azMmtoPixel(-list.b,list.scrx,list.scry)
    point_lun,lun,0
    line=''
    in=0.
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,in
       list.phipixphi0=in
       printw,id,block+' loaded'
    endif else list.phipixphi0=azMmtoPixel(-list.b,list.scrx,list.scry)
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$IntCor2Dinfo:'
if (flagi[where(fields eq block)])[0] then begin
    list.IntCor2Dinfo=DefaultIntCor2Dinfo()
    list.IntCor2Dinfo.Pff=list.scrx*list.scrx*1d6
    
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
        in=intarr(6)
        readf,lun,in
        list.IntCor2Dinfo.fftype=in[0]
        list.IntCor2Dinfo.Atype=in[1]
        list.IntCor2Dinfo.Ptype=in[2]
        list.IntCor2Dinfo.Ltype=in[3]
        list.IntCor2Dinfo.n=in[4]
        list.IntCor2Dinfo.s=in[5]
        
        in=fltarr(8)
        readf,lun,in
        list.IntCor2Dinfo.Pff=in[0]
        list.IntCor2Dinfo.Plindeg=in[1]
        list.IntCor2Dinfo.mmult=in[2]
        list.IntCor2Dinfo.dmono=in[3]
        list.IntCor2Dinfo.flatnormphi=in[4]
        list.IntCor2Dinfo.flatnormpolar=in[5]
        list.IntCor2Dinfo.flatD=in[6]
        list.IntCor2Dinfo.muL=in[7]

       printw,id,block+' loaded'
    endif else list.IntCor2Dinfo=DefaultIntCor2Dinfo()
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$TILTP:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       in=0
       readf,lun,in
       nrTiltRing=in

         nrTilt=lonarr(nrTiltRing+1)
         readf,lun,nrTilt
         if flags then begin ;Add
          if nrTilt[0] ne 0 then begin
          ; Point to the next ring if one not finished
          if (*list.nrTilt)[list.nrTiltRing] ne 0 then begin
                      list.nrTiltRing++
                      (*list.nrTilt)=[(*list.nrTilt),nrTilt]
          endif else begin
                      if list.nrTiltRing eq 0 then (*list.nrTilt)=nrTilt else $
                      (*list.nrTilt)=[(*list.nrTilt)[0:list.nrTiltRing-1],nrTilt]
          endelse

          n=total(nrTilt)
          in=fltarr(nrTiltRing+(nrTilt[nrTiltRing]<1))
          readf,lun,in
          if ptr_valid(list.TiltRingD) then (*list.TiltRingD)=[[(*list.TiltRingD)],in] else list.TiltRingD=ptr_new(in)
          in=fltarr(2,n)
          readf,lun,in
          if ptr_valid(list.TiltI) then (*list.TiltI)=[[(*list.TiltI)],[in]] else list.TiltI=ptr_new(in)
          list.nrTiltRing+=nrTiltRing
          endif
         endif else begin
          (*list.nrTilt)=nrTilt
          if nrTilt[0] ne 0 then begin
          n=total(nrTilt)
          in=fltarr(nrTiltRing+(nrTilt[nrTiltRing]<1))
          readf,lun,in
          ptr_free,list.TiltRingD
          list.TiltRingD=ptr_new(in)
          in=fltarr(2,n)
          readf,lun,in
          ptr_free,list.TiltI
          list.TiltI=ptr_new(in)
          endif
          list.nrTiltRing=nrTiltRing
         endelse

       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$CORRECTIONS:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
        names=list.SeqTags
        n=n_elements(names)
        newSeqTags=strarr(n)
        newSeqVals=bindgen(n)
        newCorMethod=lonarr(n)
        namesused=bytarr(n)
        
        ; Read (possibly incomplete corrections list)
        n_in=0L
           readf,lun,n_in
           string_in=''
           byte_in=0b
           for i=0l,n_in-1 do begin
               readf,lun,string_in
               readf,lun,byte_in
               newSeqTags[i]=string_in
            newSeqVals[i]=byte_in
            namesused[(where(names eq string_in))[0]]=1b
           endfor
           
           in=lonarr(n_in)
          readf,lun,in
          newCorMethod[0]=in
          
          ; Add missing correction names
        ind=where(~namesused,ct)
        for i=0l,ct-1 do newSeqTags[n_in+i]=names[ind[i]]
           
        ; Save
        list.SeqTags=newSeqTags
        CorSeq=CREATE_STRUCT(list.SeqTags[0], newSeqVals[0])
        for i=1,n-1 do CorSeq=CREATE_STRUCT(CorSeq,list.SeqTags[i], newSeqVals[i])

          PTR_FREE, list.CorSeq
        list.CorSeq=PTR_NEW(CorSeq)
        list.CorMethod=newCorMethod
       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$FIT_ROI:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       ROI1=intarr(2)
       readf, lun, ROI1 ; Obsolete
       ROI2=fltarr(5)
       readf, lun, ROI2

       t=where(TAG_NAMES(list) eq 'AZINTPARAM',ct)
       if ct ne 0 then begin ;xrd_bp
           list.AzIntParam=ROI2
           tmp=CHI_xval(list.AzIntParam[0:1],list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,0)
        list.backranindtt=tmp[*,1]
       endif else begin
           list.backranindD=ROI2[0:1]
           temp=value_locater((*list.xvalt)[*,0],list.backranindD)>0
        list.backranind=temp[sort(temp)]
       endelse

       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$BACKoD:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
        vi=fltarr(4)
        readf,lun,vi
        list.backmod=vi[0]
        (*(*list.ptrModel).ptrFit).backparamIO=round(vi[1:*])

           printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$FIL:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    ct=fltarr(5)
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf, lun, ct
       list.vin=ct[0]
       list.win=ct[1]
       list.r1=ct[2]
       list.fitmethod=ct[3]
       list.filorder=ct[4]
       printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$LABELS:' ;These are always replaced!!!!!!!!!!!!!
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
        alabels=0L
        readf,lun,alabels

        if list.alabels ne 0 then begin
            ptr_free,list.labels
            ptr_free,list.xlabels
            ptr_free,list.ylabels

            list.labels=ptr_new(strarr(alabels))
            list.xlabels=ptr_new(fltarr(alabels,3))
            list.ylabels=ptr_new(fltarr(alabels))
            list.alabels=alabels

            type=0B
            readf,lun,type

            for k=0l,alabels-1 do begin
                   line = ''
                 readf, lun, line
                 line=STRCOMPRESS(line) ; reduce white spaces to 1 white space
                 line=strtrim(line,2) ; delete spaces at the beginning and end
                 line=STR_SEP(line,' ') ; chop line
                 temp=0.
                   (*list.labels)[k]=line[0]
                   reads,line[1],temp
                   (*list.xlabels)[k,*]=CHI_xval(temp,list.lambda,list.scrx,list.scry,list.a,list.b,list.dist,type)
                   reads,line[2],temp
                   (*list.ylabels)[k]=temp
           endfor
        endif
           printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$MASK1:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    masknr=0
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
       readf,lun,masknr
       ;Replace data
       if masknr ne 0 then begin
          mask=fltarr(5,masknr)
          readf,lun,mask
       endif
       printw,id,block+' loaded'
    endif

    if (flagi[where(fields eq '$MASK:')])[0] then begin
        list.masknr1=masknr
        ptr_free,list.mask1
        if masknr ne 0 then list.mask1=PTR_NEW(mask)
    endif else begin
        list.masknr=masknr
        ptr_free,list.mask
        ptr_free,list.smask
        if masknr ne 0 then begin
            list.mask=PTR_NEW(mask)
            list.smask=PTR_NEW(replicate(1b,masknr))
        endif
    endelse
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$FF:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
         in=''
         readf,lun,in
         list.ffpath=in
         readf,lun,in
         list.fffile=in
         ptr_free,list.ff
         list.ffpath=AdaptDirectory(list.ffpath,mskpath,file=list.fffile)
         printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$MASKOFF:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
         in=''
         readf,lun,in
         list.maskoffpath=in
         readf,lun,in
         list.maskofffile=in
         list.maskoffpath=AdaptDirectory(list.maskoffpath,mskpath,file=list.maskofffile)
         readf,lun,in
         list.maskoffthres=in
         in=0b
         readf,lun,in
         list.maskoffdyn=in
         in=''
         readf,lun,in
         list.maskoffdynthres=in
         ptr_free,list.maskoff
         printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$SATURATION:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
         in=0.
         readf,lun,in
         list.SaturMax=in
         in=0.
         readf,lun,in
         list.SaturAspect=in
         in=0l
         readf,lun,in
         list.SaturNMax=in
         printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$ZINGERS:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
         in=0
         readf,lun,in
         list.ZingerWidth=in
         in=0.
         readf,lun,in
         list.ZingerThreshold=in
         printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$AzInt:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
         in=intarr(9)
         readf,lun,in
         list.AzIntInfo.oD=in[0]
         list.AzIntInfo.soD=in[1]
         list.AzIntInfo.Fill=in[2]
         list.AzIntInfo.type=in[3]
         list.AzIntInfo.bool=in[4]
         list.AzIntInfo.SaveRingPercentage=in[5]
         list.AzIntInfo.Batch=in[6]
         list.AzIntInfo.tiff=in[7]
         list.AzIntInfo.origformat=in[8]
         
         in=fltarr(6)
         readf,lun,in
         list.AzIntInfo.MinPix=in[0]
         list.AzIntInfo.BWP=in[1]
         list.AzIntInfo.BWT=in[2]
         list.AzIntInfo.BWD=in[3]
         list.AzIntInfo.BWQ=in[4]
         list.AzIntInfo.AzWidth=in[5]

         in=0L
         readf,lun,in
         list.AzIntInfo.outbin=in
         
         in=0b
         readf,lun,in
         list.AzIntInfo.median=in
         
         in=0b
         readf,lun,in
         list.AzIntInfo.stD=in
         
         in=0
         readf,lun,in
         list.AzIntInfo.nsectors=in
         
         if CompareVersion(versionmask,'5.3.2.1') ne -1 then begin
             in=bytarr(2)
             readf,lun,in
             list.AzIntInfo.aoD=in[0]
             list.AzIntInfo.loD=in[1]
         endif else begin
              list.AzIntInfo.aoD=0b
             list.AzIntInfo.loD=0b
         endelse
         
         if CompareVersion(versionmask,'6.0.0.1') ne -1 then begin
            in=0.
            readf,lun,in
            list.AzIntInfo.percentile=in
         endif else begin
             list.AzIntInfo.percentile=50
         endelse
         
         if CompareVersion(versionmask,'6.1.2.3') ne -1 then begin
            in=bytarr(2)
            readf,lun,in
            list.AzIntInfo.updatenorm=in[0]
            list.AzIntInfo.calcerror=in[1]
         endif else begin
             list.AzIntInfo.updatenorm=0
            list.AzIntInfo.calcerror=0
         endelse
         
         printw,id,block+' loaded'
    endif
endif
ENDELSE

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    printw,id,block+' corrupted or older version'
    nError=nError+1
    error=1B
ENDIF ELSE BEGIN
block='$ReadCCDInfo:'
if (flagi[where(fields eq block)])[0] then begin
    point_lun,lun,0
    line=''
    len=strlen(block)
    while (strmid(line,0,len) ne block) and (not eof(lun)) do readf,lun,line
    if not eof(lun) then begin
         in=0b
         readf,lun,in
         list.readinfo.showheader=in
         readf,lun,in
         list.readinfo.bnorm=in
         in=''
         readf,lun,in
         list.readinfo.int=in
         readf,lun,in
         list.readinfo.gain=in
         readf,lun,in
         list.readinfo.offset=in
         in=0
         readf,lun,in
         list.readinfo.binxs=in
         readf,lun,in
         list.readinfo.binys=in
         in=0b
         readf,lun,in
         list.readinfo.binlendian=in
         in=0
         readf,lun,in
         list.readinfo.bintype=in
         readf,lun,in
         list.readinfo.binnxpad3=in
         printw,id,block+' loaded'
    endif
endif
ENDELSE

free_lun,lun
printw,id,'Mask loaded from: '+path+' with '+stringr(nError)+' errors'
end;pro LoadMask
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RestoreCorrection,list,tiff=tiff,xdinr=xdinr
ptr_free,list.tiff
list.tiff=ReadCCD(list.path,list.file,info=list.readinfo,xdinr=xdinr)
if n_elements(tiff) ne 0 then tiff=list.tiff
list.CorDone[*]=0

ptr_free,list.tiffvalid
ptr_free,list.tiffb
end;pro RestoreCorrection
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro OpenFile,ev,prompt,fileptr=fileptr,rebin=rebin,info=info,xdinr=xdinr

;look for previous existing data (open new file)
widget_control,ev.top,get_uValue=list
if keyword_set(info) then begin
    info=list.readinfo
    if prompt eq 0 then info.showheader=0
endif

; ----Choose file----
if prompt then begin
    path=''
    filter='*'+FileFormat(2)
    sortformat,filter,list.Format

    path=CAPDIALOG_PICKFILE(path=list.path,file=list.file,title='Read diffraction pattern...',filter=filter)
    if path EQ '' then return
    error=CutPath(path,path=path,file=file,ext=ext)
    Format=ext
    file=file+ext
    if ptr_valid(fileptr) then ptr_free,fileptr
endif else begin
    path=list.path
    file=list.file
    Format=list.Format
endelse

; ----Read file----
if ptr_valid(fileptr) then begin
    tiff=fileptr
    if ptr_valid(tiff) eq 0 then return
    if keyword_set(rebin) then begin
        RestoreCorrection,list,tiff=tiff,xdinr=xdinr
        rebinCCD,list,tiff,rebin
    endif
endif else tiff=ReadCCD(path,file,info=info,xdinr=xdinr)

if ptr_valid(tiff) eq 0 then return

; ----File statistics----
tiffs=DimSizeFull(*tiff,2)
;list.xdinr=xdinr

; ----Make dynamic and static window----
result2=ReadINI('$StatSize:')
stat=result2.(0)
xscroll=stat;<tiffs[1]
yscroll=stat;<tiffs[2]

ratio=tiffs[1]/float(tiffs[2])
rebuild=round(list.stath/ratio) ne list.statv ;Ratio changed
if rebuild then begin
    ; rebuilt because we can't change dynamic window scroll size

    if (ratio gt 1) then begin; H size  gt V size
        list.stath=stat
        list.statv=round(stat/ratio)
    endif else begin
        list.stath=round(stat*ratio)
        list.statv=stat
    endelse
    widget_control,list.draw,xsize=list.stath,ysize=list.statv
endif

list.hrange=0>list.hrange<(tiffs[1]-1)
list.vrange=0>list.vrange<(tiffs[2]-1)
list.sfh=tiffs[1]/float(list.stath)
list.sfv=tiffs[2]/float(list.statv)

if list.zoom then begin
    xs=list.hrange[1]-list.hrange[0]+1
    ys=list.vrange[1]-list.vrange[0]+1
endif else begin
    xs=tiffs[1]
    ys=tiffs[2]
    list.hrange=[0,tiffs[1]-1]
    list.vrange=[0,tiffs[2]-1]
endelse
list.dynh=round(xs/list.sf)
list.dynv=round(ys/list.sf)

drawdyn=list.drawdyn
drawdynindex=list.drawdynindex
pixindex2=list.pixindex2
dynh=list.dynh
dynv=list.dynv
sf=list.sf
tmp=CreateDyn(drawdyn,drawdynindex,pixindex2,dynh,dynv,xscroll,yscroll,xs,ys,$
    sf,widget_info(ev.top,FIND_BY_UNAME='drawbasedyn'),widget_info(ev.top,FIND_BY_UNAME='window scale'),rebuild=rebuild)
list.drawdyn=drawdyn
list.drawdynindex=drawdynindex
list.pixindex2=pixindex2
list.dynh=dynh
list.dynv=dynv
list.sf=sf

; ----Update----
listtype=list.listtype
case listtype of
0:  begin
    list.listtype=1
    list.tiffs=PTR_NEW(tiffs)

    ; ----Calculate image dimension dependent variables----
    list.GridSpPix=[(tiffs[1]-1-2.*list.GridOffset[0])/(list.GridDim[0]-1.),$
           (tiffs[2]-1-2.*list.GridOffset[1])/(list.GridDim[1]-1.)]
    list.scrx=list.GridSpReal[0]/list.GridSpPix[0]/1000.
    list.scry=list.GridSpReal[1]/list.GridSpPix[1]/1000.

    ; ----Update Window----
    ID1=widget_info(ev.top,FIND_BY_UNAME='base1')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu2')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu4')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu50a')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu50b')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu50c')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu51')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu52')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu53')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu54')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu55')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu56')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu57')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu58')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='menu59')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='Write 2D Pattern')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='Save Valid Pixels')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='Auto Load')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='Load Mask')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='Add Mask')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='Save Mask')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='Load PDF')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='Load Multiple PDF')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='Save Image')
    widget_control,ID1,sensitive=1
    ID1=widget_info(ev.top,FIND_BY_UNAME='csliderbase')
    widget_control,ID1,sensitive=1

    if (list.smellipse and 1) eq 1 then begin
        ID1=widget_info(ev.top,FIND_BY_UNAME='Show Debye Marker')
        widget_control,ID1,set_value='# Show Debye Marker'
    endif
    if list.cross eq 1 then begin
        ID1=widget_info(ev.top,FIND_BY_UNAME='Show Points')
        widget_control,ID1,set_value='# Show Points'
    endif
    if total(list.ShowTie) ge 1 then begin
        ID1=widget_info(ev.top,FIND_BY_UNAME='Show Tie')
        widget_control,ID1,set_value='# Show SpaDist'
    endif
    if list.ShowTilt eq 1 then begin
           ID1=widget_info(ev.top,FIND_BY_UNAME='Show Tilt')
           widget_control,ID1,set_value='# Show Calib'
    endif
    if list.MarkBadPixels eq 1 then begin
           ID1=widget_info(ev.top,FIND_BY_UNAME='Show bad pixels')
           widget_control,ID1,set_value='# Show bad pixels'
    endif
    if list.bPdyn eq 1 then begin
        ID1=widget_info(ev.top,FIND_BY_UNAME='Show DynMark')
        widget_control,ID1,set_value='# Show DynMark'
    endif

    WIDGET_CONTROL, list.drawdyn, SET_DRAW_VIEW=[(list.dynh-list.stath)/2.,(list.dynv-list.statv)/2.]
    widget_control,list.draw,get_value=drawindex
    list.drawindex=drawindex

    endcase
1:    begin
    ptr_free,list.tiff
    ptr_free,list.tiffvalid
    *list.tiffs=tiffs
    ptr_free,list.tiffb
    ptr_free,list.ff
    endcase
endcase
list.tiff=tiff
list.path=path
list.file=file
list.Format=Format

list.thetadyn=0
list.dangle=0.
list.bP=1
list.xs=0
list.ys=0
list.xd=0
list.yd=0
list.CorDone[*]=0
list.QueTimeStamp=Systime(1)

; ----Finalize----
if list.autocor then AutoCorrect,list

RefreshDisplay, ev.top,list
if listtype eq 1 then mode,ev,list.mode,/refresh
if list.AzIntInfo.aoD then az_int_auto,ev

; Flush Event-queue
widget_control,ev.id,/clear_events
end; pro OpenFile
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UpdateFiles,files,list,delflag=delflag,delind=delind,nodel=nodel,nocheck=nocheck

delflag=0B
delind=-1l

; No files in list from previous search
if ~ptr_valid(list.PTR_files) then begin
    ; Check files
    if ~keyword_set(nocheck) then begin
        ind=where(CheckFile(files),ct)
        if ct eq 0 then return,0b
        files=files[ind]
    endif

    list.PTR_files=ptr_new(temporary(files))
    return,1b
endif

; Check which files are new and which of the old list can be kept
bnew=bytarr(n_elements(files))
bkeep=bytarr(n_elements(*list.PTR_files))
for i=0l,n_elements(files)-1 do begin
    ind=where(*list.PTR_files eq files[i],ct)
    if ct eq 0 then bnew[i]=1b $
    else bkeep[ind]=1b
endfor

; Delete files for list
if not keyword_set(nodel) then begin
    ind=where(bkeep,nkeep,comp=delind,ncomp=ndel)
    if ndel ne 0 then delflag=1
    if nkeep eq 0 then ptr_free,list.PTR_files $
    else *list.PTR_files=(*list.PTR_files)[ind]
endif

; New files?
ind=where(bnew,ct)
if ct eq 0 then return,0
files=files[ind]

; Check files
if ~keyword_set(nocheck) then begin
    ind=where(CheckFile(files),ct)
    if ct eq 0 then return,0b
    files=files[ind]
endif

; Add files
if ptr_valid(list.PTR_files) then $
    *list.PTR_files=[*list.PTR_files,temporary(files)] $
else list.PTR_files=ptr_new(temporary(files))

return,1b
end;function UpdateFiles
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DetectFiles2,list,ev,notimer=notimer,delind=delind,nodel=nodel
; ----Visualize timer----
if not keyword_set(notimer) then begin
    if list.chari eq 2 then list.chari=0 else list.chari++
    widget_control,ev.top,set_uvalue=list
    widget_control,ev.top,TLB_SET_TITLE=list.title+(list.char)[list.chari]
endif

; ----Find new entrees in directory----
nocheck=0B
bsort=1B
case list.filtermethod of
    -3:        begin
            filesi=notimer
            count=n_elements(filesi)
            nocheck=1B
            bsort=0B
            endcase
    -2:        begin
            count=0
            filesi=''
            for i=0l,n_elements(list.filter)-1 do begin
                filesi=[filesi,file_search(list.path+list.filter[i],count=ct)]
                count=count+ct
            endfor
            endcase
    -1:        filesi=file_search(list.path+list.filter,count=count)
    1:          begin
            filesi=file_search(list.path+(*list.filter),count=count)
            nocheck=list.nocheck
            endcase
    2:         begin
            ; Accept all files, even if they don't exist:
            filesi=(*list.filter)
            count=n_elements(filesi)
            nocheck=list.nocheck
            bsort=0B
               ;count=0
               ;filesi=''
               ;lof=(*list.filter)
               ;for i=0l,n_elements(lof)-1 do begin
            ;     file=file_search(lof[i],count=count2)
            ;     if count2 eq 1 then begin
            ;          filesi=[filesi,file]
            ;          count=count+1
            ;     endif
               ;endfor
               ;if count ne 0 then filesi=filesi[1:*]
            endcase
    3:        begin
            skipinfo=str_sep((*list.filter),',')
            filesi=file_search(list.path+skipinfo[0],count=count)
            nskip=long(skipinfo[1])
            if count gt nskip then begin
                filesi = SortListFiles(filesi,list.sortmethod,list.sortseparator)
                bsort=0b

                if n_elements(skipinfo) ge 3 then begin
                    nmax=long(skipinfo[2])
                    filesi=filesi[nskip:(count-1)<(nskip+nmax-1)]
                endif else filesi=filesi[nskip:*]
                count=n_elements(filesi)

            endif else count=0
            nocheck=list.nocheck
            endcase
endcase

tmp=where(TAG_NAMES(list) eq 'NFILES',uflag)

if count eq 0 then begin ;No files in directory
    ret=1
    if PTR_VALID(list.PTR_files) and ~keyword_set(nodel) then begin ;There were files previously
       ptr_free,list.PTR_files

       if uflag ne 0 then begin
         list.nfiles=0
         widget_control,ev.top,set_uvalue=list
       endif

       ; ----Update droplist----
       ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
       widget_control,ID,set_value=''
       ; ----Update file-counter----
       ID=widget_info(ev.top,FIND_BY_UNAME='FilesText')
       widget_control,ID,set_value='0 files'
       ret=2
    endif

       if not keyword_set(notimer) then  WIDGET_CONTROL, ev.id, TIMER=list.sec
    return,ret
endif

; Sort files
if bsort then files = SortListFiles(filesi,list.sortmethod,list.sortseparator) $
else files = temporary(filesi)

; Update the files list
bnew=UpdateFiles(files,list,delflag=delflag,delind=delind,nodel=nodel,nocheck=nocheck)
widget_control,ev.top,set_uvalue=list
n=n_elements(*list.PTR_files)

if bnew or delflag then begin
    ; Update Widgets
    ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
    widget_control,ID,set_value=ptr_valid(list.PTR_files)?(*list.PTR_files):''
    widget_control,ID,SET_LIST_TOP=n-1
    ID=widget_info(ev.top,FIND_BY_UNAME='FilesText')
    widget_control,ID,set_value=string(n,format='(I0," files")')
    
    ; Update main window
    if (list.update eq 1) and widget_info(list.top,/valid_id) then $
        MainUpdateCManual,(*list.PTR_files)[n-1],list.top,list.ev
endif

; Generate new file search event
if not keyword_set(notimer) then WIDGET_CONTROL, ev.id, TIMER=list.sec

if ~bnew and delflag and ptr_valid(list.PTR_files) then return,2 ; Files deleted but still some left
return,1 ; New files (and possibily some delete)
end;function DetectFiles2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DetectFiles,list,ev,notimer=notimer,delind=delind

; ----Visualize timer----
if not keyword_set(notimer) then begin
    if list.chari eq 2 then list.chari=0 else list.chari++
    widget_control,ev.top,set_uvalue=list
    widget_control,ev.top,TLB_SET_TITLE=list.title+(list.char)[list.chari]
endif

; ----Find new entrees in directory----
nocheck=0B
bsort=1B
case list.filtermethod of
    -3:        begin
            filesi=notimer
            count=n_elements(filesi)
            nocheck=1B
            bsort=0B
            endcase
    -2:        begin
            count=0
            filesi=''
            for i=0l,n_elements(list.filter)-1 do begin
                filesi=[filesi,file_search(list.path+list.filter[i],count=ct)]
                count=count+ct
            endfor
            endcase
    -1:        filesi=file_search(list.path+list.filter,count=count)
    1:          begin
            filesi=file_search(list.path+(*list.filter),count=count)
            nocheck=list.nocheck
            endcase
    2:         begin
            ; Accept all files, even if they don't exist:
            filesi=(*list.filter)
            count=n_elements(filesi)
            nocheck=list.nocheck
            bsort=0B
               ;count=0
               ;filesi=''
               ;lof=(*list.filter)
               ;for i=0l,n_elements(lof)-1 do begin
            ;     file=file_search(lof[i],count=count2)
            ;     if count2 eq 1 then begin
            ;          filesi=[filesi,file]
            ;          count=count+1
            ;     endif
               ;endfor
               ;if count ne 0 then filesi=filesi[1:*]
            endcase
    3:        begin
            skipinfo=str_sep((*list.filter),',')
            filesi=file_search(list.path+skipinfo[0],count=count)
            nskip=long(skipinfo[1])
            if count gt nskip then begin
                filesi = SortListFiles(filesi,list.sortmethod,list.sortseparator)
                bsort=0b

                if n_elements(skipinfo) ge 3 then begin
                    nmax=long(skipinfo[2])
                    filesi=filesi[nskip:(count-1)<(nskip+nmax-1)]
                endif else filesi=filesi[nskip:*]
                count=n_elements(filesi)

            endif else count=0
            nocheck=list.nocheck
            endcase
endcase

tmp=where(TAG_NAMES(list) eq 'NFILES',uflag)

if count eq 0 then begin ;No files in directory
    ret=1
    if PTR_VALID(list.PTR_list) then begin ;There were files previously
       DestroyLL,list.PTR_list,-1

       if uflag ne 0 then begin
         list.PTR_cur=list.PTR_list
         widget_control,ev.top,set_uvalue=list
       endif

       ; ----Update droplist----
       ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
       widget_control,ID,set_value=''
       ; ----Update file-counter----
       ID=widget_info(ev.top,FIND_BY_UNAME='FilesText')
       widget_control,ID,set_value='0 files'
       ret=2
    endif

    WIDGET_CONTROL, ev.id, TIMER=list.sec
    return,ret
endif
if bsort then files = SortListFiles(filesi,list.sortmethod,list.sortseparator) $
else files = filesi

new_entree=UpdateLL(files,list,delflag=delflag,delind=delind,nocheck=nocheck)
if uflag ne 0 then begin
    ; Point to first file
    list.PTR_cur=list.PTR_list
    ; Point to where we were (or NULL when to few files)
    for i=1,list.nFiles do if ptr_valid(list.PTR_cur) then $
        list.PTR_cur=(*list.PTR_cur).next
endif

widget_control,ev.top,set_uvalue=list

path=new_entree[0]
if (path eq '') then begin ;No New files
    ret=1
    if delflag and (files[0] ne '') then begin ;Files deleted but still some left
       ; ----Update droplist----
       ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
       widget_control,ID,set_value=files
       widget_control,ID,SET_LIST_TOP=n_elements(files)-1
       ; ----Update file-counter----
       ID=widget_info(ev.top,FIND_BY_UNAME='FilesText')
       widget_control,ID,set_value=string(n_elements(files),format='(I0," files")')
       ret=2
    endif
    WIDGET_CONTROL, ev.id, TIMER=list.sec
    return,ret
endif

; ----Update list----
ID=widget_info(ev.top,FIND_BY_UNAME='FilesList')
if files[0] ne '' then files=[files,new_entree] else files=new_entree
widget_control,ID,set_value=files
widget_control,ID,SET_LIST_TOP=n_elements(files)-1

; ----Update file-counter----
ID=widget_info(ev.top,FIND_BY_UNAME='FilesText')
widget_control,ID,set_value=string(n_elements(files),format='(I0," files")')

; ----Update main window----
if (list.update eq 1) and widget_info(list.top,/valid_id) then $
    MainUpdateCManual,path,list.top,list.ev

; ----Generate new event----
if not keyword_set(notimer) then WIDGET_CONTROL, ev.id, TIMER=list.sec
return,1
end;function DetectFiles
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MainUpdateCore,list,main,path,nr=nr,tdrequired=tdrequired

if ~widget_info(main,/valid) then return
widget_control,main,get_uvalue=toplist
error=CutPath(path,path=path2,file=temp,ext=ext)

case FileFormat(0,path) of
2:    begin
    toplist.path=path2
    toplist.file=temp+ext
    toplist.Format=ext
    widget_control,main,set_uvalue=toplist

    OpenFile,list.ev,0,/info,xdinr=nr
    endcase
else:begin
    if ptr_valid(toplist.CHIchild) and not keyword_set(tdrequired) then begin
        ID=(*toplist.CHIchild)[0]
        widget_control,ID,get_uvalue=toplist
        toplist.path=path2
        toplist.file=temp+ext
        toplist.Format=ext
        widget_control,ID,set_uvalue=toplist

        OpenCHI,{ID:ID,top:ID},0,xdinr=nr
    endif
    endcase
endcase

end;pro MainUpdateCore
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MainUpdate,list,ev,tdrequired=tdrequired
; Update main window
if (list.update eq 1) then begin
    
    if ~ptr_valid(list.PTR_files) then return
    path=(*list.PTR_files)[ev.index]
    
    error = WIDGET_INFO(list.top,/VALID_ID)
    if error eq 0 then return
    if ~widget_info(list.top,/valid) then return
    MainUpdateCore,list,list.top,path,nr=ev.index,tdrequired=tdrequired
endif
end;pro MainUpdate
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MainUpdateCManual,path,main,ev,fileptr=fileptr

if ~widget_info(main,/valid) then return
widget_control,main,get_uvalue=toplist
error=CutPath(path,path=path2,file=temp,ext=ext)

case FileFormat(0,path) of
2:    begin
    toplist.path=path2
    toplist.file=temp+ext
    toplist.Format=ext
    widget_control,main,set_uvalue=toplist

    ; Copy file because the main window might autocorrect it
    if arg_present(fileptr) then fileptrnew=ptr_new((size(fileptr,/type) eq 10)?(*fileptr):fileptr)
    OpenFile,ev,0,fileptr=fileptrnew
    endcase
else:begin
    if ptr_valid(toplist.CHIchild) then begin
        ID=(*toplist.CHIchild)[0]
        OpenCHI,{ID:ID,top:ID},0,fileptr=fileptr
    endif
    endcase
endcase

end;pro MainUpdateCManual
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MainUpdateCoreCalc,list,main,path,nr=nr,weights=weights

if ~widget_info(main,/valid) then return
widget_control,main,get_uvalue=toplist

n=n_elements(weights)-1
bbreak=keyword_set(nr)
for i=0l,n do begin

    error=CutPath(path[i],path=path2,file=temp,ext=ext)
    file=temp+ext
    
    case FileFormat(0,path[i]) of
    2:    begin
        tiff=ReadCCD(path2,file,info=list.readinfo);,xdinr=i)
        if ~ptr_valid(tiff) then begin
            bbreak=1b 
            printw,list.OutID,'Error reading: '+path2+file
        endif else begin
            bbreak=0b 
            printw,list.OutID,path2+file
            if i eq 0 then begin
                fileptr=tiff
                *fileptr*=weights[i]/float(n)
            endif else begin
                (*fileptr)+=(*tiff)*weights[i]/float(n)
                ptr_free,tiff
            endelse
        endelse
        endcase
    else:begin
        if ~ptr_valid(toplist.CHIchild) then return

        error=ReadScan(path2,file,ext,list.OutID,spe=spe,stdspe=speerror,xval=xval,title=title,$
        YTITLE=YTITLE,xrange=xrange,yrange=yrange,type=xtype,ROIAzimuth=ROIAzimuth)

        if bbreak then begin
            s=size(spe)
            coeff=rebin(reform(weights,1,n+1),s[1],n+1,/sample)
            spe=total(coeff*spe[*,nr]/float(n),2,/pres)
            speerror=sqrt(total(coeff^2*(spe[*,nr])^2,2,/pres))/float(n)
        endif else begin
            spe*=weights[i]/float(n)
            speerror*=weights[i]/float(n)
        endelse
        
        if i eq 0 then begin
            fileptr={xval:xval,$
            spe:spe,$
            speerror:speerror,$
            title:title+'_tomoavg',$
            xtitle:ChiXtitle(xtype)  ,$
            ytitle:ytitle,$
            xrange:xrange,$
            yrange:yrange,$
            xtype:xtype,$
            ROIAzimuth:ROIAzimuth}
        endif else begin
            fileptr.spe+=spe
            fileptr.speerror=sqrt((fileptr.speerror)^2+speerror^2)
        endelse
        
        endelse
    endcase

    if bbreak then break
endfor

MainUpdateCManual,path[0],main,list.ev,fileptr=fileptr
end;pro MainUpdateCoreCalc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro az_int_auto,ev
FORWARD_FUNCTION AzIntWarpBatch

widget_control,ev.top,get_uvalue=list
printw,list.OutID,'Azimuthal integration:'
if ~ptr_valid(list.batchreturnp) then begin
    printw,list.OutID,'One manual azimuthal integration needs to be performed in order for automatic integration to work.'
    return
endif else printw,list.OutID,'Uses settings of previous integration.'

tstart=systime(1)
d=AzIntWarpBatch(list.batchreturnp,list.tiff,list.tiffvalid)
printw,list.OutID,'AzIntWarp processing time:'+string(systime(1)-tstart)

ROIAzimuth=[0,2*!pi]
temp=min(d[1,*],y1)
temp=max(d[1,*],y2)

if list.AzIntInfo.calcerror then speerror=reform(d[2,*]) else speerror=reform(sqrt(d[1,*]))
    
fileptr={xval:reform(d[0,*]),spe:reform(d[1,*]),speerror:speerror,ROIAzimuth:(*list.batchreturnp).ROIAzimuth,ytitle:'Intensity (a.u.)',$
            xtitle:ChiXtitle(list.AzIntInfo.type),title:file_basename(list.file,list.Format),xtype:list.AzIntInfo.type,$
            xrange:[0,n_elements(d)/2-1],yrange:[y1[0],y2[0]]}
    
MainUpdateCManual,'c:\dummy.chi',ev.top,ev,fileptr=fileptr
            
end;pro az_int_auto
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
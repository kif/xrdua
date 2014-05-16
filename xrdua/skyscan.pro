pro SkyScanLogFile,file,nfiles,rows,cols,depth,dir,prefix,pixelsize,range

if openw_safe(lun,file) then return

; Skyscan:
;  range=360 degrees
;  step=0.4 degrees/step
;  => nfiles = 900 (= range/step)
; XRDUA:
;  range=360 degrees
;  step=0.4 degrees/step
;  => nfiles = 901 (= range/step + 1)
; XRDUA -> Skyscan:
;  range=360.4
;  step=0.4 degrees/step
;  => nfiles = 901 (= range/step)

range=float(range)
range+=range/(nfiles-1)
pixelsize=float(pixelsize)

; Instrument information
printf,lun,'[System]'
printf,lun,'Scanner=Skyscan1172'
printf,lun,'Instrument S/N=001'
printf,lun,'Hardware version=A'
printf,lun,'Software=Version 1. 5 (build 6)'
printf,lun,'Home directory=C:\dummy'
printf,lun,'Source Type=Hamamatsu 100/250'
printf,lun,'Camera=Hamamatsu 10Mp camera'
printf,lun,'Camera Pixel Size (um)=   11.73' ; Camera pixel size
printf,lun,'CameraXYRatio='+string(pixelsize,format='(f0)') ; Camera pixel aspect ratio
printf,lun,'Incl.in lifting (um/mm)=0.0000' ; Compensation for object horizontal shift

; Measurement information
printf,lun,'[Acquisition]'
printf,lun,'Data directory='+dir
printf,lun,'Filename Prefix='+prefix
printf,lun,'Number of Files= '+string(nfiles,format='(I0)'); Number of CCD images
printf,lun,'Source Voltage (kV)= 100'
printf,lun,'Source Current (uA)= 100'
printf,lun,'Number of Rows= '+string(rows,format='(I0)'); Number of pixels (vertical) for CCD image
printf,lun,'Number of Columns= '+string(cols,format='(I0)'); Number of pixels (horizontal) for CCD image
printf,lun,'Image Pixel Size (um)=    '+string(pixelsize,format='(f0)'); um/pixel
printf,lun,'Object to Source (mm)=100000' ; Focal length
printf,lun,'Camera to Source (mm)=100001' ; Object is at the focal point
printf,lun,'Vertical Object Position (mm)=0.000' ;?
printf,lun,'Optical Axis (line)= '+string(round(rows/2.),format='(I0)') ; Optical axis corresponds to the line in the image with highest resolution,
                                                        ; this is the place where the X-ray beam hits the camera with exactly 90° angle.
printf,lun,'Filter=No'
printf,lun,'Image Format=TIFF'
printf,lun,'Depth (bits)='+string(depth,format='(I0)');CCD image pixel depth
printf,lun,'Screen LUT=0'
printf,lun,'Exposure (ms)=  1000'
printf,lun,'Rotation Step (deg)='+string(range/nfiles,format='(f0)')
printf,lun,'Frame Averaging=ON (4)'
printf,lun,'Random Movement=OFF (10)'
printf,lun,'Use 360 Rotation=NO'
printf,lun,'Geometrical Correction=ON'
printf,lun,'Camera Offset=OFF'
printf,lun,'Median Filtering=ON'
printf,lun,'Flat Field Correction=ON'
printf,lun,'Rotation Direction=CC'
printf,lun,'Scanning Trajectory=ROUND'
printf,lun,'Type Of Motion=STEP AND SHOOT'
printf,lun,'Study Date and Time=Feb 27, 2008  10:33:40'
printf,lun,'Scan duration=01:00:00'

; Reconstruction information
printf,lun,'[Reconstruction]'
printf,lun,'Reconstruction Program=NRecon'
printf,lun,'Program Version=Version: 1.5.1.3'
printf,lun,'Program Home Directory=C:\dummy'
printf,lun,'Dataset Origin=Skyscan1172'
printf,lun,'Dataset Prefix='+prefix
printf,lun,'Dataset Directory='+dir
printf,lun,'Time and Date=Feb 27, 2008  11:19:53'
printf,lun,'First Section=0'
printf,lun,'Last Section='+string(rows-1,format='(I0)')
printf,lun,'Reconstruction duration per slice (seconds)=1.000000'
printf,lun,'Postalignment=0.00'
printf,lun,'Section to Section Step=1'
printf,lun,'Sections Count='+string(rows,format='(I0)')
printf,lun,'Result File Type=BMP'
printf,lun,'Result File Header Length (bytes)=1134'
printf,lun,'Result Image Width (pixels)='+string(cols,format='(I0)')
printf,lun,'Result Image Height (pixels)='+string(cols,format='(I0)')
printf,lun,'Pixel Size (um)='+string(pixelsize,format='(f0)')
printf,lun,'Reconstruction Angular Range (deg)='+string(range,format='(f0)')
printf,lun,'Use 180+=OFF'
printf,lun,'Angular Step (deg)='+string(range/nfiles,format='(f0)')
printf,lun,'Smoothing=2'
printf,lun,'Ring Artifact Correction=7'
printf,lun,'Draw Scales=ON'
printf,lun,'Object Bigger than FOV=OFF'
printf,lun,'Reconstruction from ROI=OFF'
;printf,lun,'ROI Top (pixels)=1865'
;printf,lun,'ROI Bottom (pixels)=0'
;printf,lun,'ROI Left (pixels)=0'
;printf,lun,'ROI Right (pixels)=1812'
;printf,lun,'ROI reference length=2000'
printf,lun,'Undersampling factor=1'
printf,lun,'Threshold for defect pixel mask (%)=0'
printf,lun,'Beam Hardening Correction (%)=0'
printf,lun,'CS Static Rotation (deg)=0.0'
printf,lun,'Mininum for CS to Image Conversion=0.0000'
printf,lun,'Maximum for CS to Image Conversion=0.5000'
printf,lun,'HU Calibration=OFF'
printf,lun,'BMP LUT=0'
printf,lun,'Cone-beam Angle Horiz.(deg)=0.000000'
printf,lun,'Cone-beam Angle Vert.(deg)=0.000000'

free_lun,lun
end;pro SkyScanLogFile
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SkyScanTIFF_write,file,TIFF,error=error,offset=offset,insertdata=insertdata

error=1b

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
  print,'Error in writing ',file
  message,!ERR_STRING,/info
  return
ENDIF else begin
  bigendian=TIFF.ifh.byteorder eq '4D4D'x

  ; File with the required size
  openw,lun,file,/get_lun
  writeu,lun,bytarr(strucbytesize(TIFF))
  close,lun
  openu,lun,file, SWAP_IF_BIG_ENDIAN=~bigendian, SWAP_IF_LITTLE_ENDIAN=bigendian

  writeu,lun,TIFF.IFH.byteorder
  writeu,lun,TIFF.IFH.tiffversion
  ifdoffset=TIFF.IFH.ifdoffset
  writeu,lun,ifdoffset

  ImageStrucIDn=1
  while ifdoffset ne 0 do begin
    point_lun,lun,ifdoffset
    n=TIFF.(ImageStrucIDn).ifdlength
    writeu,lun,n

    for i=1,n do begin
      ; Write extended tag data when value is a pointer
      if n_tags(TIFF.(ImageStrucIDn).(i)) gt 4 then begin
        writeu,lun,TIFF.(ImageStrucIDn).(i).code
        writeu,lun,TIFF.(ImageStrucIDn).(i).type
        writeu,lun,TIFF.(ImageStrucIDn).(i).length
        writeu,lun,TIFF.(ImageStrucIDn).(i).value

        point_lun,-lun,pos
        point_lun,lun,TIFF.(ImageStrucIDn).(i).value
        writeu,lun,TIFF.(ImageStrucIDn).(i).heap
        point_lun,lun,pos

        if TIFF.(ImageStrucIDn).(i).code eq 273 then StripOffsets=TIFF.(ImageStrucIDn).(i).heap
        if TIFF.(ImageStrucIDn).(i).code eq 279 then StripByteCounts=TIFF.(ImageStrucIDn).(i).heap
      endif else begin
        writeu,lun,TIFF.(ImageStrucIDn).(i)

        if TIFF.(ImageStrucIDn).(i).code eq 273 then StripOffsets=TIFF.(ImageStrucIDn).(i).value
        if TIFF.(ImageStrucIDn).(i).code eq 279 then StripByteCounts=TIFF.(ImageStrucIDn).(i).value
      endelse
    endfor

    ifdoffset=TIFF.(ImageStrucIDn).ifdoffset
    writeu,lun,ifdoffset

    ; Write image data
    n=0L
    for i=0l,n_elements(StripOffsets)-1 do begin
      point_lun,lun,StripOffsets[i]
      writeu,lun,TIFF.(ImageStrucIDn).image[n:n+StripByteCounts[i]-1]
      n+=StripByteCounts[i]
    endfor

    ImageStrucIDn++
  endwhile

  if n_elements(offset) ne 0 and n_elements(insertdata) ne 0 then begin
    point_lun,lun,offset
    writeu,lun,insertdata
  endif
  free_lun,lun
endelse

error=0b
end;pro SkyScanTIFF_write
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SkyScanTIFF_insert,TIFF,offset,nbytes

for i=0l,n_tags(TIFF)-1 do begin
  if TIFF.(i).ifdoffset ge offset then TIFF.(i).ifdoffset+=offset
  n=TIFF.(i).ifdlength
  for j=1,n do begin
    if n_tags(TIFF.(i).(j)) gt 4 then $
      if TIFF.(i).(j).value ge offset then TIFF.(i).(j).value+=offset
    if TIFF.(i).(j).code eq 273 then begin
      if n_tags(TIFF.(i).(j)) gt 4 then begin
        ind=where(TIFF.(i).(j).heap ge offset,ct)
        if ct ne 0 then TIFF.(i).(j).heap[ind]+=offset
      endif else begin
        ind=where(TIFF.(i).(j).value ge offset,ct)
        if ct ne 0 then TIFF.(i).(j).value[ind]+=offset
      endelse
    endif
  endfor
endfor

end;pro SkyScanTIFF_insert
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SkyScanTIFF_patch,file

; Read Tiff
TIFF=read_tiff_structure(file,error=error)
if error then return

; Patch tiff for Skyscan
; Offset 8: the 8byte zero terminating string "SkyScan" is inserted
offset=8
insertdata=[83b, 107b, 121b,  83b,  99b,  97b, 110b, 0b]
SkyScanTIFF_insert,TIFF,offset,n_elements(insertdata)

; Write Tiff
SkyScanTIFF_write,file,TIFF,error=error,offset=offset,insertdata=insertdata

end;pro SkyScanTIFF_patch
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SkyScan_SelectFiles,top,list,$
    sortmethod=sortmethod,separator=separator,filter=filter,$
    outpath=outpath,description=description,files=files
    
files=Select_Files(top,list,$
    sortmethod=sortmethod,separator=separator,filter=filter,$
    outpath=outpath,description=description)
if n_elements(files) eq 0 then return,ptr_new()

; Read first file
return,ReadCCD(files[0])
end;function SkyScan_SelectFiles
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SkyScan_BockInterpol,block,iknown,i
; i and iknown between [0,n-1]
; block: k x l x (m <= n)
; 
; Interpolate between the different block[*,*,j]
; to get block[*,*,i]

if n_elements(iknown) eq 1 then return,block

j=value_locate(iknown,i)>0
if iknown[j] eq i then return,block[*,*,j] $
else begin
    x1=iknown[j]
    x2=iknown[j+1]
    
    ; Interpolate (nearest neighbour)
    if abs(i-x2) lt abs(i-x1) then j++
    return,block[*,*,j]
    
    ; Interpolate (linear)
    ;return,(block[*,*,j+1]-block[*,*,j])*((i-x1)/(x2-x1))+block[*,*,j]
endelse
end;function SkyScan_BockInterpol    
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SkyScan_NRecon,top,list,$
        data=data,dark=dark,flat=flat,filename=filename,$
        sortmethod=sortmethod,separator=separator

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
  printw,list.OutID,'Error in preparing Skyscan reconstruct.'
  message,!ERR_STRING,/info
  return
ENDIF

list2={path:list.path,file:list.file,Format:list.Format,Platform:list.Platform,OutID:list.OutID}

bmemory=n_elements(data) ne 0
if bmemory then begin
    ; Check datasize
    ; first dim: number of projections
    ; second dim: projection col size
    ; second dim: projection row size
    sdata=DimSize(data,3)
    nproj=sdata[0]
    sdata=sdata[1:2]
    bytesize=arraybytesize(data[0],FLOAT=FLOAT,SIGNED=SIGNED)
    
    ; dark images
    bdark=n_elements(dark) ne 0
    
    ; flat field
    bflat=n_elements(flat) ne 0
    
endif else begin
    
    ; Select projections
    files=Select_Files(top,list2,sortmethod=sortmethod,$
        separator=separator,description='projections')
    if files[0] eq '' then return
    nproj=n_elements(files)
    
    ; Check datasize
    data=ReadCCD(files[0],/static,error=error)
    if error then return
    bytesize=arraybytesize(data[0],FLOAT=FLOAT,SIGNED=SIGNED)
    sdata=DimSize(data,2)
    
    ; Read dark images
    printw,list2.OutID,'Reading dark images...'
    tmp=CutPath(files[0],path=path,file=file,ext=ext)
    list2.path=path
    list2.file=file
    list2.Format=ext
    files2=Select_Files(top,list2,sortmethod=sortmethod,$
        separator=separator,description='dark images')
    bdark=files2[0] ne ''
    if bdark then begin
        dark=ReadCCD(files2[0],/static)
        bdark=array_equal(sdata,DimSize(dark,2))
        if bdark then begin
            for i=1,n_elements(files2)-1 do begin
                add=ReadCCD(files2[i],/static)
                if array_equal(sdata,DimSize(add,2)) then dark=[[[dark]],[[add]]]
            endfor
        endif
    endif
    
    ; Read flat field images
    printw,list2.OutID,'Reading flat fields...'
    files2=Select_Files(top,list2,sortmethod=sortmethod,$
        separator=separator,description='flat fields')
    bflat=files2[0] ne ''
    if bflat then begin
        flat=ReadCCD(files2[0],/static)
        bflat=array_equal(sdata,DimSize(flat,2))
        if bflat then begin
            for i=1,n_elements(files2)-1 do begin
                add=ReadCCD(files2[i],/static)
                if array_equal(sdata,DimSize(add,2)) then flat=[[[flat]],[[add]]]
            endfor
        endif
    endif
endelse

; Output data format
if bflat then begin
    bytesize>=4
    FLOAT=1
endif
if bdark and ~SIGNED then begin
    SIGNED=1b
    if bytesize eq 4 then bytesize=8
    if bytesize eq 2 then bytesize=4
endif

; SkyScan format restrictions
FLOAT=0b
bytesize<=2

; Prepare for dark and flat field correction
if bdark then begin
    ndark=(DimSize(dark,3))[2] 
    idark=findgen(ndark)*(nproj-1.)/(ndark-1.)
endif else ndark=0
if bflat then begin
    nflat=(DimSize(flat,3))[2] 
    iflat=findgen(nflat)*(nproj-1.)/(nflat-1.)
endif else nflat=0

; Split in chunks
chunksizemin=bytesize*sdata[0]
filesize=chunksizemin*sdata[1]
datasize=filesize*ulong64(nproj)
nbin=round(2.*datasize/1024ULL^3)>1
repeat begin
    ConvertSize,filesize,unitfile
    totalsize=datasize
    ConvertSize,totalsize,unittotal
    chunksize=(datasize/float(nbin)) > chunksizemin
    ConvertSize,chunksize,unitchunk
    
    msg='Process '+string(nproj,format='(I)')+' files of '+string(filesize,format='(f5.1)')+unitfile+' (total='+$
    string(totalsize,format='(f5.1)')+unittotal+', split in '+string(nbin,format='(I)')+' times '+string(chunksize,format='(f5.1)')+unitchunk+')'
    msg=[msg,'Dark field correction: '+((bdark)?'YES':'NO')+' (interpolate between '+strtrim(ndark)+' dark images)']
    msg=[msg,'Flat field correction: '+((bflat)?'YES':'NO')+' (interpolate between '+strtrim(nflat)+' flat fields)']
    
    answer=DIALOG_MESSAGE(strcompress(msg), /QUESTION, /CANCEL)
    if answer eq 'Cancel' then return
    if answer eq 'No' then begin
        nbin=fix(PromptNumber(string(nbin,format='(I0)'),top,'Number of chunks:'))
        nbin>=1
    endif
endrep until answer eq 'Yes'

ind=sdata[1]/float(nbin)*indgen(nbin+1)
ind0=ceil(ind[0:nbin-1])<(sdata[1]-1)
ind1=floor(ind[1:*]-1)<(sdata[1]-1)

; Outpath file names
;     files: base+4-digit number+'.tif'
;    logfile: base+'.log'
pathout=DIALOG_PICKFILE(path=list2.path,title='Save in directory...',/directory)
if pathout EQ '' then return
pathout+=path_sep()+'vol'+MakeNumber(indgen(nbin)+1)+path_sep()
file_mkdir,pathout
if bmemory then prefix=filename+'_' else prefix='f_'
filesout=string(lindgen(nproj),format='(I04)')+'.tif'

; Process
progressBar = Obj_New("PROGRESSBAR",title="Converting images...")
progressBar -> Start

for i=0L,nproj-1 do begin

    IF progressBar -> CheckCancel() THEN BEGIN
        progressBar -> Destroy
        printw,list2.OutID,'Process CANCELLED!'
        RETURN
    ENDIF

    if bmemory then datai=reform(data[i,*,*],sdata) $
    else datai=ReadCCD(files[i],/static)
    
    if array_equal(DimSize(datai,2),sdata) then begin
        ; Dark and flatfield correction: (I-DA)/(FF-DA)
        if bdark then begin
            DA=SkyScan_BockInterpol(dark,idark,i)
            datai-=DA
        endif
        if bflat then begin
            FF=float(SkyScan_BockInterpol(flat,iflat,i))
            if bdark then FF-=DA
            FF=mean(FF)/FF
            datai*=FF
        endif
        if bdark or bflat then begin
            ; Remove nan and inf
            tmp=where(~finite(datai),tmp2)
            if tmp2 ne 0 then datai[tmp]=0b
            
            ; Output format
            case bytesize of
            1: datai=byte(datai)
            2: datai=SIGNED?fix(datai):uint(datai)
            4: datai=FLOAT?float(datai):(SIGNED?long(datai):ulong(datai))
            8: datai=FLOAT?double(datai):(SIGNED?long64(datai):ulong64(datai))
            endcase
        endif

        ; Write chunks
        for j=0l,nbin-1 do begin
            out=pathout[j]+prefix+filesout[i]
              WRITE_TIFF_SAFE,out,datai[*,ind0[j]:ind1[j]],error=error
              if error then printw,list2.OutID,'File not written: '+out ;else patchtiff,out
          endfor
    endif else printw,list2.OutID,'File not read: '+files[i]
    
    ; Display progress
    progressBar -> Update, (i+1.)/nproj*100
endfor

progressBar -> Destroy

; Write logfiles
printw,list2.OutID,'Writing Skyscan log files...'
logfile=pathout+prefix+'.log'
nrows=ind1-ind0+1
pixelsize=1.0
range=180.
for j=0l,nbin-1 do $
      SkyScanLogFile,logfile[j],nproj,nrows[j],sdata[0],bytesize*8,pathout[j],prefix,pixelsize,range

printw,list2.OutID,'Done.'

list.path=list2.path
widget_control,top,set_uvalue=list
end;pro SkyScan_NRecon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

function LoadDLMs,dlms,outid

CATCH, Error_status
IF Error_status NE 0 THEN begin
    printw,outid,'Failed to load DLM: '
    printwerrorstate,outid
    return,1
endif

cd,c=path
path+=path_sep()+'DLM'+path_sep()

dest=!dlm_path+path_sep()
for i=0,n_elements(dlms)-1 do begin
    files=file_search(path+dlms[i],dlms[i]+'*.{so,dll,dlm}',count=n)
    if n eq 0 then continue
    
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printw,outid,'Failed to load DLM: '
        printwerrorstate,outid
        return,1
    endif

    file_copy,files,dest,/overwrite
    DLM_REGISTER,path+dlms[i]+path_sep()+dlms[i]+'.dlm'
    DLM_LOAD,dlms[i]
endfor

return,0

end;function LoadDLMs
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CreateDLLMainCode,path,name,functions,messages,includes,defines,outid

cd,path,c=c

file=name+'.c'
openw,lun,file,/get_lun,error=error
if error ne 0 then begin
    printw,outid,'Could not open '+file
    printwerrorstate,outid
    cd,c
    return,1
endif

printf,lun,'#include <stdio.h>'
printf,lun,'#include "idl_export.h"'

printf,lun,''
printf,lun,''
printf,lun,''
printf,lun,'/* Custom includes */'
for i=0,n_elements(includes)-1 do $
    printf,lun,'#include ',includes[i]

printf,lun,''
printf,lun,''
printf,lun,''
printf,lun,'/* Mark functions as exports */'
printf,lun,'#if defined(_WIN32)'
printf,lun,'#define MY_EXPORT __declspec(dllexport)'
printf,lun,'#elif defined(_WIN64)'
printf,lun,'#define MY_EXPORT __declspec(dllexport)'
printf,lun,'#else'
printf,lun,'#define MY_EXPORT'
printf,lun,'#endif'

printf,lun,''
printf,lun,''
printf,lun,''
printf,lun,'/* Custom definitions */'

for i=0,n_elements(defines)-1 do $
    printf,lun,'#define ',defines[i]
    
printf,lun,''
printf,lun,''
printf,lun,''
printf,lun,'/* prototype for IDL_Load */'
printf,lun,'MY_EXPORT int IDL_Load( void );'
printf,lun,''
printf,lun,'/*'
printf,lun,' * Define message codes and their corresponding printf(3) format'
printf,lun,' * strings. Note that message codes start at zero and each one is'
printf,lun,' * one less that the previous one. Codes must be monotonic and'
printf,lun,' * contiguous.'
printf,lun,' */'

if n_elements(messages) eq 0 then printf,lun,'/* No messages */' else begin
    printf,lun,'static IDL_MSG_DEF msg_arr[] ='
    printf,lun,'{'
endelse
for i=0,n_elements(messages)-1 do begin
    msg=messages[i]
    printf,lun,'#define ',msg.name,'                       ',strtrim(i eq 0?0:-i,2)
    printf,lun,'  {  "',msg.name,'",   "%N',msg.str,'" },'
endfor
if n_elements(messages) ne 0 then printf,lun,'};'

printf,lun,''
printf,lun,''
printf,lun,''
printf,lun,'/*'
printf,lun,' * The load function fills in this message block handle with the'
printf,lun,' * opaque handle to the message block used for this module. The other'
printf,lun,' * routines can then use it to throw errors from this block.'
printf,lun,' */'
printf,lun,'static IDL_MSG_BLOCK msg_block;'

format='("'+string(9b)+'",A)'
for i=0,n_elements(functions)-1 do begin
    func=functions[i]
    printf,lun,''
    printf,lun,''
    printf,lun,''
    printf,lun,'/* Implementation of the ',func.name,' IDL procedure */'
    printf,lun,''
    printf,lun,'MY_EXPORT ',((func.type eq 1)?'void':'IDL_VPTR'),' ',func.name,'(int argc, IDL_VPTR *argv)'
    printf,lun,'{'
    printf,lun,*func.code,format=format
    printf,lun,'}'
endfor

printf,lun,''
printf,lun,''
printf,lun,''
printf,lun,'MY_EXPORT int IDL_Load(void)'
printf,lun,'{'
printf,lun,'  /*'
printf,lun,'   * These tables contain information on the functions and procedures'
printf,lun,'   * that make up the ',name,' DLM. The information contained in these'
printf,lun,'   * tables must be identical to that contained in ',name,'.dlm'
printf,lun,'   */'

indpro=where(functions.type,npro,comp=indfun,ncomp=nfun)

if nfun ne 0 then printf,lun,'  static IDL_SYSFUN_DEF2 function_addr[] = {'
for i=0,nfun-1 do begin
    func=functions[indfun[i]]
    printf,lun,'    { ',func.name,', "',strupcase(func.name),'", ',func.minargs,', ',func.maxargs,', 0, 0},'
endfor
if nfun ne 0 then printf,lun,'  };'

if npro ne 0 then printf,lun,'  static IDL_SYSFUN_DEF2 procedure_addr[] = {'
for i=0,npro-1 do begin
    func=functions[indpro[i]]
    printf,lun,'    { (IDL_FUN_RET) ',func.name,', "',strupcase(func.name),'", ',func.minargs,', ',func.maxargs,', 0, 0},'
endfor
if npro ne 0 then printf,lun,'  };'

printf,lun,''
printf,lun,'  /*'
printf,lun,'   * Create a message block to hold our messages. Save its handle where'
printf,lun,'   * the other routines can access it.'
printf,lun,'   */'
printf,lun,'  if (!(msg_block = IDL_MessageDefineBlock("'+name+'",'
printf,lun,'                       IDL_CARRAY_ELTS(msg_arr), msg_arr)))'
printf,lun,'    return IDL_FALSE;'
printf,lun,''
printf,lun,'  /*'
printf,lun,'   * Register our routine. The routines must be specified exactly the same'
printf,lun,'   * as in '+name+'.dlm.'
printf,lun,'   */'

if nfun eq 0 then $
    printf,lun,'  return IDL_SysRtnAdd(procedure_addr, FALSE, IDL_CARRAY_ELTS(procedure_addr));' $
else if npro eq 0 then $
    printf,lun,'  return IDL_SysRtnAdd(function_addr, TRUE, IDL_CARRAY_ELTS(function_addr));' $
else begin
    printf,lun,'  return IDL_SysRtnAdd(function_addr, TRUE, IDL_CARRAY_ELTS(function_addr))'
    printf,lun,'    && IDL_SysRtnAdd(procedure_addr, FALSE, IDL_CARRAY_ELTS(procedure_addr));'
endelse

printf,lun,'}'

free_lun,lun

cd,c
return,0
end;function CreateDLLMainCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CreateDLM,path,name,description,functions,outid

cd,path,c=c

file=name+'.dlm'
openw,lun,file,/get_lun,error=error
if error ne 0 then begin
    printw,outid,'Could not open '+file
    printwerrorstate,outid
    cd,c
    return,1
endif

printf,lun,"MODULE "+name
printf,lun,"DESCRIPTION "+description
printf,lun,"VERSION "+version()
printf,lun,"SOURCE XRDUA"
printf,lun,string(SYSTIME(/JULIAN),FORMAT='(C("BUILD_DATE ",CMOA," ",CDI," ",CYI))')
for i=0,n_elements(functions)-1 do begin
    func=functions[i]
    printf,lun,(func.type eq 1?"PROCEDURE ":"FUNCTION "),strupcase(func.name)," ",func.minargs," ",func.maxargs;," ",func.options
endfor
free_lun,lun

cd,c
return,0
end;function CreateDLM
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CreateDLMfunction,name,minargs,maxargs,code,procedure=procedure
return,{name:name,type:keyword_set(procedure),minargs:strtrim(minargs,2),maxargs:strtrim(maxargs,2),code:ptr_new(code)}
end;function CreateDLMfunction
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CreateDLMmessage,name,str
return,{name:name,str:str}
end;function CreateDLMmessage
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InitMakeDllEnv,path,ccinfo

ccinfo->GetProperty,SDKbasedir=SDKbasedir,$
    SDKsubdirpath=SDKsubdirpath,$
    SDKsubdirinclude=SDKsubdirinclude,$
    SDKsubdirlib_=SDKsubdirlib_,$
    SDKsubdirlibpath=SDKsubdirlibpath,$
    basedir=basedir,$
    subdirpath=subdirpath,$
    subdirinclude=subdirinclude,$
    subdirlib_=subdirlib_,$
    subdirlibpath=subdirlibpath

dir_path=''
if strtrim(SDKsubdirpath,2) ne '' then dir_path+=strjoin(SDKbasedir+strsplit(SDKsubdirpath,';',/extract),';')
if strtrim(subdirpath,2) ne '' then dir_path+=(dir_path eq ''?'':';')+strjoin(basedir+strsplit(subdirpath,';',/extract),';')

dir_include=''    
if strtrim(SDKsubdirinclude,2) ne '' then dir_include+=strjoin(SDKbasedir+strsplit(SDKsubdirinclude,';',/extract),';')
if strtrim(subdirinclude,2) ne '' then dir_include+=(dir_include eq ''?'':';')+strjoin(basedir+strsplit(subdirinclude,';',/extract),';')

dir_lib=''
if strtrim(SDKsubdirlib_,2) ne '' then dir_lib+=strjoin(SDKbasedir+strsplit(SDKsubdirlib_,';',/extract),';')
if strtrim(subdirlib_,2) ne '' then dir_lib+=(dir_lib eq ''?'':';')+strjoin(basedir+strsplit(subdirlib_,';',/extract),';')

dir_libpath=''
if strtrim(SDKsubdirlibpath,2) ne '' then dir_libpath+=strjoin(SDKbasedir+strsplit(SDKsubdirlibpath,';',/extract),';')
if strtrim(subdirlibpath,2) ne '' then dir_libpath+=(dir_libpath eq ''?'':';')+strjoin(basedir+strsplit(subdirlibpath,';',/extract),';')

envold_PATH = GETENV( 'PATH' )
envnew_PATH = 'PATH='+dir_path+';'+envold_PATH
SETENV, 'PATH='+envnew_PATH

envold_INCLUDE = GETENV( 'INCLUDE' )
envnew_INCLUDE = 'INCLUDE='+dir_include+';'+envold_INCLUDE
SETENV, 'INCLUDE='+envnew_INCLUDE

envold_LIB = GETENV( 'LIB' )
envnew_LIB = 'LIB='+dir_lib+';'+envold_LIB
SETENV, 'LIB='+envnew_LIB

envold_LIBPATH = GETENV( 'LIBPATH' )
envnew_LIBPATH = 'LIBPATH='+dir_lib+';'+envold_LIBPATH
SETENV, 'LIBPATH='+envnew_LIBPATH

cd,path,c=c

return,{c:c,envold_PATH:envold_PATH,envold_INCLUDE:envold_INCLUDE,envold_LIB:envold_LIB,envold_LIBPATH:envold_LIBPATH}
end;function InitMakeDllEnv
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FinishMakeDllEnv,info
cd,info.c
SETENV, 'PATH='+info.envold_PATH
SETENV, 'INCLUDE='+info.envold_INCLUDE
SETENV, 'LIB='+info.envold_LIB
SETENV, 'LIBPATH='+info.envold_LIBPATH
end;pro FinishMakeDllEnv
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CreateDLML,name,ccinfo,description,functions,messages,includes,defines,outid,debug=debug

; Copy files to user directory
cd,c=sourcepath
sourcepath+=path_sep()+'DLM'+path_sep()+name+path_sep()

path = XRDUAdlmDIR()
file_copy,sourcepath,path,/recursive,/force,/overwrite

; Create DLM
path += name+path_sep()
error=CreateDLM(path,name,description,functions,outid)
if error ne 0 then begin
    heap_free,functions
    return,1
endif

error=CreateDLLMainCode(path,name,functions,messages,includes,defines,outid)
if error ne 0 then begin
    heap_free,functions
    return,1
endif

info=InitMakeDllEnv(path,ccinfo)

sourcefiles=name
files=file_search(path,'*.c',count=n)
if n ne 0 then sourcefiles=FILE_BASENAME(files,'.c')

ccinfo->GetProperty,cflags=cflags,lflags=lflags,cc=cc,ld=ld

pathkeep = path
make_dll, sourcefiles, functions.name, $
              input_directory=path, $
              output_directory=path, $
              compile_directory=path,$
              dll_path=path,$
              cc=cc,$
              ld=ld,$
              EXTRA_CFLAGS=cflags,$
              EXTRA_LFLAGS=lflags,$
              /platform_extension,nocleanup=keyword_set(debug)
path = pathkeep 

FinishMakeDllEnv,info

; Distribute build
files=file_search(path,'*.{so,dll,dlm}',count=n)
if n ne 2 then return,1
file_copy,files,!dlm_path,/force,/overwrite

return,0

end;function CreateDLML
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DLMDefaultDefines
return,['MY_NOERROR IDL_GettmpLong(M_TM_NOERROR)',$
         'MY_ERROR IDL_GettmpLong(M_TM_ERROR)']
end;function DLMDefaultDefines
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DLMDefaultMessages,name
return,[CreateDLMmessage('M_TM_NOERROR','No error occured in '+name+'.'),$
        CreateDLMmessage('M_TM_ERROR','Error occured in '+name+'.'),$
        CreateDLMmessage('M_TM_BYTESNOTCORRECT','Array should be 16 bit.')]
end;function DLMDefaultMessages
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DLMDefaultFunctions,name
return,[CreateDLMfunction(name+'_check','0','0','return MY_NOERROR;')]
end;function DLMDefaultFunctions
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CreateReadMarIP345,ccinfo,outid,debug=debug

name='ReadMarIP345'
description='Read packed Marresearch images. Based on procedures by Mark Rivers (March, 2006).'

includes=['"mar3xx_pck.h"']
defines=[DLMDefaultDefines(),'MAXBUFF 1024']
messages=DLMDefaultMessages(name)

tab=string(9b)
codeReadMarIP345=[$
    'char buffer [MAXBUFF];',$
    'INT16 *data;',$
    'char *filename;',$
    'FILE *input;',$
    '',$
    '// First parameter: string',$
    'filename =  IDL_VarGetString(argv[0]);',$
    '',$
    '// Second parameter: 16bit array',$
    'IDL_ENSURE_ARRAY(argv[1]);',$
    'if (argv[1]->value.arr->elt_len != 2){',$
    tab+'IDL_MessageFromBlock(msg_block, M_TM_BYTESNOTCORRECT, IDL_MSG_RET);',$
    tab+'return MY_ERROR;',$
    '}',$
    'data = (INT16 *) argv[1]->value.arr->data;',$
    '',$
    '// Open file',$
    'input = fopen(filename, "rb");',$
    'if (input == NULL) {',$
    tab+'sprintf_s(buffer,MAXBUFF,"Error opening input file %s\n", filename);',$
    tab+'IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_RET, buffer);',$
    tab+'return MY_ERROR;',$
    '}',$
    '',$
    '// Get data',$
    'get_pck(input, data);',$
    '',$
    '// Close file',$
    'fclose(input);',$
    '',$
    'return MY_NOERROR;'$
            ]

    
functions=[ DLMDefaultFunctions(name),$
            CreateDLMfunction('XRDUADLM_ReadMarIP345','2','2',codeReadMarIP345)]

return,CreateDLML('ReadMarIP345',ccinfo,description,functions,messages,includes,defines,outid,debug=debug)

end;function CreateReadMarIP345
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function compilerconfig::Init, _REF_EXTRA = _extra

; Initialize the superclass.
IF (self->IDLitComponent::Init() ne 1) THEN $
   RETURN, 0

; Register properties.
;
; * Only registered properties will show up in the property sheet.
; * <identifier> must match self.<identifier>.

self->RegisterProperty, 'CC', /STRING, $
   NAME = 'Compiler', DESCRIPTION = 'Compiler with flags.'
self->RegisterProperty, 'LD', /STRING, $
   NAME = 'Linker', DESCRIPTION = 'Linker with flags.'
self->RegisterProperty, 'CFLAGS', /STRING, $
   NAME = 'Extra compiler flags', DESCRIPTION = 'Extra compiler flags.'
self->RegisterProperty, 'LFLAGS', /STRING, $
   NAME = 'Extra linker flags', DESCRIPTION = 'Extra linker flags.'
   
self->RegisterProperty, 'SDKBASEDIR', /STRING, $
   NAME = 'SDK directory', DESCRIPTION = 'Location of the Software Development Kit.'
self->RegisterProperty, 'SDKSUBDIRPATH', /STRING, $
   NAME = 'SDK to %PATH%', DESCRIPTION = 'The SDK subdirectories added to environement variable %PATH%.'
self->RegisterProperty, 'SDKSUBDIRINCLUDE', /STRING, $
   NAME = 'SDK to %INCLUDE%', DESCRIPTION = 'The SDK subdirectories added to environement variable %INCLUDE%.'
self->RegisterProperty, 'SDKSUBDIRLIB_', /STRING, $
   NAME = 'SDK to %LIB%', DESCRIPTION = 'The SDK subdirectories added to environement variable %LIB%.'
self->RegisterProperty, 'SDKSUBDIRLIBPATH', /STRING, $
   NAME = 'SDK to %LIBPATH%', DESCRIPTION = 'The SDK subdirectories added to environement variable %LIBPATH%.'
self->RegisterProperty, 'BASEDIR', /STRING, $
   NAME = 'IDE', DESCRIPTION = 'Location of the Integrated Development Environment.'
self->RegisterProperty, 'SUBDIRPATH', /STRING, $
   NAME = 'IDE to %PATH%', DESCRIPTION = 'The IDE subdirectories added to environement variable %PATH%.'
self->RegisterProperty, 'SUBDIRINCLUDE', /STRING, $
   NAME = 'IDE to %INCLUDE%', DESCRIPTION = 'The IDE subdirectories added to environement variable %INCLUDE%.'
self->RegisterProperty, 'SUBDIRLIB_', /STRING, $
   NAME = 'IDE to %LIB%', DESCRIPTION = 'The IDE subdirectories added to environement variable %LIB%.'
self->RegisterProperty, 'SUBDIRLIBPATH', /STRING, $
   NAME = 'IDE to %LIBPATH%', DESCRIPTION = 'The IDE subdirectories added to environement variable %LIBPATH%.'

   
; Set any property values.
self->SetProperty, _EXTRA = _extra
self->SetDefaults

RETURN, 1

end;function compilerconfig::Init
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro compilerconfig::Cleanup
self->IDLitComponent::Cleanup
end;pro compilerconfig::Cleanup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro compilerconfig::GetProperty, $
    OS=OS,nbit=nbit,$
    SDKbasedir=SDKbasedir,$
    SDKsubdirpath=SDKsubdirpath,$
    SDKsubdirinclude=SDKsubdirinclude,$
    SDKsubdirlib_=SDKsubdirlib_,$
    SDKsubdirlibpath=SDKsubdirlibpath,$
    basedir=basedir,$
    subdirpath=subdirpath,$
    subdirinclude=subdirinclude,$
    subdirlib_=subdirlib_,$
    subdirlibpath=subdirlibpath, $
    cflags=cflags,$
    lflags=lflags,$
    cc=cc,$
    ld=ld,$
   _REF_EXTRA = _extra

IF (arg_present(OS)) THEN OS = self.OS
IF (arg_present(nbit)) THEN nbit = self.nbit
IF (arg_present(SDKbasedir)) THEN SDKbasedir = self.SDKbasedir
IF (arg_present(SDKsubdirpath)) THEN SDKsubdirpath = self.SDKsubdirpath
IF (arg_present(SDKsubdirinclude)) THEN SDKsubdirinclude = self.SDKsubdirinclude
IF (arg_present(SDKsubdirlib_)) THEN SDKsubdirlib_ = self.SDKsubdirlib_
IF (arg_present(SDKsubdirlibpath)) THEN SDKsubdirlibpath = self.SDKsubdirlibpath
IF (arg_present(basedir)) THEN basedir = self.basedir
IF (arg_present(subdirpath)) THEN subdirpath = self.subdirpath
IF (arg_present(subdirinclude)) THEN subdirinclude = self.subdirinclude
IF (arg_present(subdirlib_)) THEN subdirlib_ = self.subdirlib_
IF (arg_present(subdirlibpath)) THEN subdirlibpath = self.subdirlibpath
IF (arg_present(cflags)) THEN cflags = self.cflags
IF (arg_present(lflags)) THEN lflags = self.lflags
IF (arg_present(cc)) THEN cc = self.cc
IF (arg_present(ld)) THEN ld = self.ld

IF (n_elements(_extra) gt 0) THEN $
   self->IDLitComponent::GetProperty, _EXTRA = _extra

end;pro compilerconfig::GetProperty
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro compilerconfig::SetProperty, $
    OS=OS,nbit=nbit,$
    SDKbasedir=SDKbasedir,$
    SDKsubdirpath=SDKsubdirpath,$
    SDKsubdirinclude=SDKsubdirinclude,$
    SDKsubdirlib_=SDKsubdirlib_,$
    SDKsubdirlibpath=SDKsubdirlibpath,$
    basedir=basedir,$
    subdirpath=subdirpath,$
    subdirinclude=subdirinclude,$
    subdirlib_=subdirlib_,$
    subdirlibpath=subdirlibpath, $
    cflags=cflags,$
    lflags=lflags,$
    cc=cc,$
    ld=ld,$
   _REF_EXTRA = _extra

IF (n_elements(OS) ne 0) THEN self.OS = OS
IF (n_elements(nbit) ne 0) THEN self.nbit = nbit
IF (n_elements(SDKbasedir) ne 0) THEN self.SDKbasedir = SDKbasedir
IF (n_elements(SDKsubdirpath) ne 0) THEN self.SDKsubdirpath = SDKsubdirpath
IF (n_elements(SDKsubdirinclude) ne 0) THEN self.SDKsubdirinclude = SDKsubdirinclude
IF (n_elements(SDKsubdirlib_) ne 0) THEN self.SDKsubdirlib_ = SDKsubdirlib_
IF (n_elements(SDKsubdirlibpath) ne 0) THEN self.SDKsubdirlibpath = SDKsubdirlibpath
IF (n_elements(basedir) ne 0) THEN self.basedir = basedir
IF (n_elements(subdirpath) ne 0) THEN self.subdirpath = subdirpath
IF (n_elements(subdirinclude) ne 0) THEN self.subdirinclude = subdirinclude
IF (n_elements(subdirlib_) ne 0) THEN self.subdirlib_ = subdirlib_
IF (n_elements(subdirlibpath) ne 0) THEN self.subdirlibpath = subdirlibpath
IF (n_elements(cflags) ne 0) THEN self.cflags = cflags
IF (n_elements(lflags) ne 0) THEN self.lflags = lflags
IF (n_elements(cc) ne 0) THEN self.cc = cc
IF (n_elements(ld) ne 0) THEN self.ld = ld

self->IDLitComponent::SetProperty, _EXTRA = _extra

end;pro compilerconfig::SetProperty
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro compilerconfig::SetDefaults
self.default_name = self.name
self.default_OS = self.OS
self.default_nbit = self.nbit
self.default_SDKbasedir = self.SDKbasedir
self.default_SDKsubdirpath = self.SDKsubdirpath
self.default_SDKsubdirinclude = self.SDKsubdirinclude
self.default_SDKsubdirlib_ = self.SDKsubdirlib_
self.default_SDKsubdirlibpath = self.SDKsubdirlibpath
self.default_basedir = self.basedir
self.default_subdirpath = self.subdirpath
self.default_subdirinclude = self.subdirinclude
self.default_subdirlib_ = self.subdirlib_
self.default_subdirlibpath = self.subdirlibpath
self.default_cflags = self.cflags
self.default_lflags = self.lflags
self.default_cc = self.cc
self.default_ld = self.ld
end;pro compilerconfig::SetDefaults
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro compilerconfig::UseDefaults
self.name = self.default_name
self.OS = self.default_OS
self.nbit = self.default_nbit
self.SDKbasedir = self.default_SDKbasedir
self.SDKsubdirpath = self.default_SDKsubdirpath
self.SDKsubdirinclude = self.default_SDKsubdirinclude
self.SDKsubdirlib_ = self.default_SDKsubdirlib_
self.SDKsubdirlibpath = self.default_SDKsubdirlibpath
self.basedir = self.default_basedir
self.subdirpath = self.default_subdirpath
self.subdirinclude = self.default_subdirinclude
self.subdirlib_ = self.default_subdirlib_
self.subdirlibpath = self.default_subdirlibpath
self.cflags = self.default_cflags
self.lflags = self.default_lflags
self.cc = self.default_cc
self.ld = self.default_ld
end;pro compilerconfig::USEDefaults
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro compilerconfig::Export,path,file

file=SelFile(path,file,'*.cfg','Export compiler settings')
if openw_safe(lun,file) then return

printf,lun,'NAME='+self.name
printf,lun,'OS='+self.OS
printf,lun,'NBIT='+strtrim(self.nbit,2)
printf,lun,'SDKBASEDIR='+self.SDKbasedir
printf,lun,'SDKSUBDIRPATH='+self.SDKsubdirpath
printf,lun,'SDKSUBDIRINCLUDE='+self.SDKsubdirinclude
printf,lun,'SDKSUBDIRLIB_='+self.SDKsubdirlib_
printf,lun,'SDKSUBDIRLIBPATH='+self.SDKsubdirlibpath
printf,lun,'BASEDIR='+self.basedir
printf,lun,'SUBDIRPATH='+self.subdirpath
printf,lun,'SUBDIRINCLUDE='+self.subdirinclude
printf,lun,'SUBDIRLIB_='+self.subdirlib_
printf,lun,'SUBDIRLIBPATH='+self.subdirlibpath
printf,lun,'CFLAGS='+self.cflags
printf,lun,'CFLAGS='+self.lflags
printf,lun,'CC='+self.cc
printf,lun,'LD='+self.ld

free_lun,lun

tmp=cutpath(file,path=path,file=file,ext=ext)
file+=ext

end;pro compilerconfig::Export
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro compilerconfig::Import,path,file

file=MDIALOG_PICKFILE(path=path,file=file,filter='*.cfg',title='Export compiler settings')
if openr_safe(lun,file) then return

ex={dummy:0}
line=''
while not eof(lun) do begin
    readf,lun,line
    tmp=strsplit(line,'=',/extract,count=n)
    if n le 1 then continue
    if n gt 2 then tmp[1]=strjoin(tmp[1:*],'=')
    prop=tmp[0]
    val=tmp[1]
    case prop of
    'nbit':    val=fix(val)
    else:
    else:
    endcase
    ex=create_struct(ex,prop,val)
endwhile

self->SetProperty , _EXTRA=ex
free_lun,lun

tmp=cutpath(file,path=path,file=file,ext=ext)
file+=ext

end;pro compilerconfig::Import
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro compilerconfig__define
    class={compilerconfig,INHERITS IDLitComponent,$
        OS:'',nbit:0,SDKbasedir:'',SDKsubdirpath:'',SDKsubdirinclude:'',SDKsubdirlib_:'',$
        SDKsubdirlibpath:'',basedir:'',subdirpath:'',subdirinclude:'',subdirlib_:'',$
        subdirlibpath:'',cflags:'',lflags:'',cc:'',ld:'',$
        ; Defaults
        default_name:'',default_OS:'',default_nbit:0,default_SDKbasedir:'',default_SDKsubdirpath:'',default_SDKsubdirinclude:'',default_SDKsubdirlib_:'',$
        default_SDKsubdirlibpath:'',default_basedir:'',default_subdirpath:'',default_subdirinclude:'',default_subdirlib_:'',$
        default_subdirlibpath:'',default_cflags:'',default_lflags:'',default_cc:'',default_ld:''}
end;pro compilerconfig__define
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CompilerConfigurations,name,OS,nbit,SDKbasedir,SDKsubdirpath,SDKsubdirinclude,$
        SDKsubdirlib_,SDKsubdirlibpath,basedir,subdirpath,subdirinclude,subdirlib_,subdirlibpath,$
        cflags,lflags,cc,ld
; OS: linux, darwin, sunos, Win32

return,OBJ_NEW('compilerconfig',name=name,description='Compiler & linker configuration',OS=OS,nbit=nbit,$
    SDKbasedir=SDKbasedir,$
    SDKsubdirpath=SDKsubdirpath,$
    SDKsubdirinclude=SDKsubdirinclude,$
    SDKsubdirlib_=SDKsubdirlib_,$
    SDKsubdirlibpath=SDKsubdirlibpath,$
    basedir=basedir,$
    subdirpath=subdirpath,$
    subdirinclude=subdirinclude,$
    subdirlib_=subdirlib_,$
    subdirlibpath=subdirlibpath,$
    cflags=cflags,$
    lflags=lflags,$
    cc=cc,$
    ld=ld)

end;function CompilerConfigurations
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CompilerConfigurationCustom
name='Custom'
OS=!version.OS
nbit=!version.memory_bits
SDKbasedir=''
SDKsubdirpath=''
SDKsubdirinclude=''
SDKsubdirlib_=''
SDKsubdirlibpath=''
basedir=''
subdirpath=''
subdirinclude=''
subdirlib_=''
subdirlibpath=''
cflags=''
lflags=''
cc=!make_dll.cc
ld=!make_dll.ld
return,CompilerConfigurations(name,OS,nbit,SDKbasedir,SDKsubdirpath,SDKsubdirinclude,$
        SDKsubdirlib_,SDKsubdirlibpath,basedir,subdirpath,subdirinclude,subdirlib_,subdirlibpath,$
        cflags,lflags,cc,ld)
end;function CompilerConfigurationCustom
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CompilerConfigurationVCExpress2008
name='Visual C++ Express 2008'
OS='Win32'
nbit=32
SDKbasedir='C:\Program Files\Microsoft SDKs\Windows\v6.0A\'
SDKsubdirpath='bin\'
SDKsubdirinclude='include\'
SDKsubdirlib_='lib\'
SDKsubdirlibpath=''
basedir='C:\Program Files (x86)\Microsoft Visual Studio 9.0\'
subdirpath='Common7\IDE\;VC\BIN\'
subdirinclude='VC\INCLUDE\'
subdirlib_='VC\LIB\'
subdirlibpath='VC\LIB\'
cflags=''
lflags='user32.lib'
cc=!make_dll.cc
ld=!make_dll.ld
return,CompilerConfigurations(name,OS,nbit,SDKbasedir,SDKsubdirpath,SDKsubdirinclude,$
        SDKsubdirlib_,SDKsubdirlibpath,basedir,subdirpath,subdirinclude,subdirlib_,subdirlibpath,$
        cflags,lflags,cc,ld)
end;function CompilerConfigurationVCExpress2008
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CompilerConfigurationVCExpress2010
name='Visual C++ Express 2010'
OS='Win32'
nbit=32
SDKbasedir='C:\Program Files\Microsoft SDKs\Windows\v6.0A\'
SDKsubdirpath='bin\'
SDKsubdirinclude='include\'
SDKsubdirlib_='lib\'
SDKsubdirlibpath=''
basedir='C:\Program Files (x86)\Microsoft Visual Studio 10.0\'
subdirpath='Common7\IDE\;VC\BIN\'
subdirinclude='VC\INCLUDE\'
subdirlib_='VC\LIB\'
subdirlibpath='VC\LIB\'
cflags=''
lflags='user32.lib'
cc=!make_dll.cc
ld=!make_dll.ld
return,CompilerConfigurations(name,OS,nbit,SDKbasedir,SDKsubdirpath,SDKsubdirinclude,$
        SDKsubdirlib_,SDKsubdirlibpath,basedir,subdirpath,subdirinclude,subdirlib_,subdirlibpath,$
        cflags,lflags,cc,ld)
end;function CompilerConfigurationVCExpress2008
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CompilerConfigurationVCExpress2008_64
name='Visual C++ Express 2008 (64 bit)'
OS='Win32'
nbit=64
SDKbasedir='C:\Program Files\Microsoft SDKs\Windows\v6.0A\'
SDKsubdirpath='bin\'
SDKsubdirinclude='include\'
SDKsubdirlib_='lib\'
SDKsubdirlibpath=''
basedir='C:\Program Files (x86)\Microsoft Visual Studio 9.0\'
subdirpath='Common7\IDE\;VC\BIN\'
subdirinclude='VC\INCLUDE\'
subdirlib_='VC\LIB\'
subdirlibpath='VC\LIB\'
cflags=''
lflags=''
cc=!make_dll.cc
ld=!make_dll.ld
return,CompilerConfigurations(name,OS,nbit,SDKbasedir,SDKsubdirpath,SDKsubdirinclude,$
        SDKsubdirlib_,SDKsubdirlibpath,basedir,subdirpath,subdirinclude,subdirlib_,subdirlibpath,$
        cflags,lflags,cc,ld)
end;function CompilerConfigurationVCExpress2008_64
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CompilerConfigurationGCC
name='GCC'
OS='linux'
nbit=32
cflags=''
lflags=''
cc=!make_dll.cc
ld=!make_dll.ld
return,CompilerConfigurations(name,OS,nbit,SDKbasedir,SDKsubdirpath,SDKsubdirinclude,$
        SDKsubdirlib_,SDKsubdirlibpath,basedir,subdirpath,subdirinclude,subdirlib_,subdirlibpath,$
        cflags,lflags,cc,ld)
end;function CompilerConfigurationGCC
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CompilerConfigurationGCC64
name='GCC'
OS='linux'
nbit=64
cflags=''
lflags=''
cc=!make_dll.cc
ld=!make_dll.ld
return,CompilerConfigurations(name,OS,nbit,SDKbasedir,SDKsubdirpath,SDKsubdirinclude,$
        SDKsubdirlib_,SDKsubdirlibpath,basedir,subdirpath,subdirinclude,subdirlib_,subdirlibpath,$
        cflags,lflags,cc,ld)
end;function CompilerConfigurationGCC64
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AllCompilerConfigurations
all=[CompilerConfigurationCustom(),$
    CompilerConfigurationVCExpress2008(),$
    CompilerConfigurationVCExpress2008_64(),$
    CompilerConfigurationVCExpress2010(),$
    CompilerConfigurationGCC(),$
    CompilerConfigurationGCC64()]

; keep only those for this platform
n=n_elements(all)
b=bytarr(n)
for i=0,n-1 do begin
    all[i]->GetProperty,OS=OS,nbit=nbit
    b[i]=OS eq !version.OS and nbit eq !version.memory_bits
endfor

ind=where(b, ct)
return,all[ind]
end;function AllCompilerConfigurations
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RegisterDLMs_Cleanup,ID
;if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
;widget_control,ID,get_uvalue=list
end;pro RegisterDLMs_Cleanup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RegisterDLMs_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
1:    begin
    widget_control,ev.id,get_value=val
    case val of
    'Exit':    begin
            WIDGET_CONTROL, ev.top, /DESTROY
            endcase
    'Keep intermediate files': begin
            widget_control,ev.top,get_uvalue=list
            list.debug=ev.select
            widget_control,ev.top,set_uvalue=list
            endcase
    'Import config':begin
            widget_control,ev.top,get_uvalue=list
            path=list.cfgpath
            file=list.cfgfile
            list.compilers[list.icompiler]->Import,path,file
            list.cfgpath=path
            list.cfgfile=file
            widget_control,ev.top,set_uvalue=list
            ID=widget_info(ev.top,find_by_uname='config')
            widget_control,ID,/REFRESH_PROPERTY
            endcase
    'Export config':begin
            widget_control,ev.top,get_uvalue=list
            path=list.cfgpath
            file=list.cfgfile
            list.compilers[list.icompiler]->Export,path,file
            list.cfgpath=path
            list.cfgfile=file
            widget_control,ev.top,set_uvalue=list
            endcase
    'Reset':begin
            widget_control,ev.top,get_uvalue=list
            list.compilers[list.icompiler]->UseDefaults
            widget_control,ev.top,set_uvalue=list
            ID=widget_info(ev.top,find_by_uname='config')
            widget_control,ID,/REFRESH_PROPERTY
            endcase
    else:    begin
            widget_control,ev.top,get_uvalue=list
            i=where(list.dlms eq val)

            ; Create DLL and DLM
            error=call_function(('Create'+list.dlms[i])[0],list.compilers[list.icompiler],list.outid,debug=list.debug)
            if error ne 0 then begin
                widget_control,ev.id,set_button=0
                return
            endif
            
            ; Register and Load DLM
            error=LoadDLMs(list.dlms[i],list.outid)
            if error ne 0 then begin
                widget_control,ev.id,set_button=0
                return
            endif

            ; Check
            widget_control,ev.id,set_button=CheckDLMloaded(list.dlms[i])
            endcase
    endcase
    endcase
13:    begin
    if ev.type eq 0 then begin
        value = WIDGET_INFO(ev.ID, COMPONENT = ev.component, $
                            PROPERTY_VALUE = ev.identifier)
        ev.component->SetPropertyByIdentifier, ev.identifier, value
    endif
    endcase
8:    begin
    widget_control,ev.top,get_uvalue=list
    list.icompiler=ev.index
    widget_control,ev.top,set_uvalue=list
    
    ID=widget_info(ev.top,find_by_uname='config')
    widget_control,ID,set_value=list.compilers[list.icompiler]
    endcase
else:
endcase

end;pro RegisterDLMs_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RegisterDLMs,top

dlms='ReadMarIP345'

widget_control,top,get_uvalue=list
list2={outid:list.outid,dlms:dlms,debug:0b,compilers:list.compilers,icompiler:0,cfgpath:'',cfgfile:''}

GROUP_LEADER=top
base=widget_base(title='Register DLMs',GROUP_LEADER=GROUP_LEADER,/column,mbar=bar)

menu1 = WIDGET_BUTTON(bar, VALUE='File', /MENU)
button = WIDGET_BUTTON(menu1, VALUE='Import config')
button = WIDGET_BUTTON(menu1, VALUE='Export config')
button = WIDGET_BUTTON(menu1, VALUE='Exit')

menu2 = WIDGET_BUTTON(bar, VALUE='Help', /MENU)

tab=widget_tab(base)
base1=widget_base(tab,/column,title='Create and register DLMs')

text=widget_label(base1,value='Restart after registring new DLMs.',/align_left)
text=widget_label(base1,value='Checkboxes indicate which DLMs are loaded.',/align_left)

baset=widget_base(base1,/nonexclusive,/column)
bdlms=lonarr(n_elements(dlms))
for i=0,n_elements(dlms)-1 do $
    bdlms[i]=widget_button(baset,value=dlms[i])

dbb=widget_button(baset,value='Keep intermediate files')

base2=widget_base(tab,/column,title='Compiler config')

n=n_elements(list2.compilers)
names=strarr(n)
for i=0,n-1 do begin
    list2.compilers[i]->GetProperty,name=name
    names[i]=name
endfor

baset=widget_base(base2,/row)
    dropcc=widget_droplist(baset,value=names,title='Preconfigured settings:')
    b=widget_button(baset,value='Reset')
    
device,get_screen_size=screen
prop=widget_propertysheet(base2,value=list2.compilers[list2.icompiler],uname='config',scr_xsize=screen[0]/2)


WIDGET_CONTROL, base, /REALIZE,set_uvalue=list2

for i=0,n_elements(dlms)-1 do $
    widget_control,bdlms[i],set_button=CheckDLMloaded(dlms[i])

widget_control,dbb,set_button=list2.debug

;widget_control,/DELAY_DESTROY
Xmanager,'RegisterDLMs',base,/NO_BLOCK,GROUP_LEADER=GROUP_LEADER,CLEANUP='RegisterDLMs_Cleanup'
end;pro RegisterDLMs
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
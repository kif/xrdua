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

function PrintLL,currentin

str=''
current=currentin
WHILE PTR_VALID(current) do begin
    str=[str,(*current).str]
    current=(*current).next
endwhile

if n_elements(str) gt 1 then str=str[1:*]

return,str
end;function PrintLL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DestroyLLEntreeI,current,index,llistAl,ptrbefore,before

topptr=PTR_NEW() ;This will get a pointer to the first not-deleted entree
ptrprev=PTR_NEW() ;This is a pointer to the previous entree in the while loop (next field must be updated when deleting)
before=0l ;Counts all deleted entrees before the first non-processed entree (ptrbefore)
baddbefore=1 ;Increment to add to before

index=index[sort(index)]
imax=index[n_elements(index)-1]
i=0l
WHILE PTR_VALID(current) and (i le imax) do begin
    bCur=current eq ptrbefore
    if bCur then baddbefore=0l

       ind=where(index eq i,count)
       if count ne 0 then begin ; Delete this entree
           if bCur then ptrbefore=(*current).next

           DestroyLLEntree,current,llistAl ; Delete entree and point to next entree
           if PTR_VALID(ptrprev) then (*ptrprev).next=current else topptr=current
           before=before+baddbefore
       endif else begin ; Keep this entree
           if PTR_VALID(topptr) eq 0 then topptr=current ; Save first non-deleted entree
           ptrprev=current ; Save current for next loop
           current=(*current).next ; Point to next entree
       endelse
       i++
endwhile

current=topptr

end;pro DestroyLLEntreeI
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FindLL,currentin,index
current=currentin
if index[0] lt 0 then begin; Take last entree
    WHILE PTR_VALID(current) do current=(*current).next
    out=(*current).str
endif else begin
    out=strarr(n_elements(index))
    i=0l
    WHILE PTR_VALID(current) do begin
       ind=where(index eq i,count)
       if count ne 0 then out[ind[0]]=(*current).str
       current=(*current).next
       i++
    endwhile
endelse
return,out
end;function FindLL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReplaceStrLL,currentin,ncut,strreplace,files=files
current=currentin
files=''
WHILE PTR_VALID(current) do begin
    (*current).str=strreplace+strmid((*current).str,ncut)
    files=[files,(*current).str]
    current=(*current).next
endwhile
if n_elements(files) gt 1 then files=files[1:*]
end;function ReplaceStrLL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeleteLLEntree,ptr,ind
prev=ptr_new()
current=ptr
ind=ind[sort(ind)]

j=0l
for i=0l,n_elements(ind)-1 do begin
    while j ne ind[i] do begin
        if ptr_valid(current) then begin
            prev=current
            current=(*current).next
        endif
        j++
    endwhile
    if not ptr_valid(current) then return

    ; cut chain
    next=(*current).next
    (*current).next=ptr_new()
    ; glue chain
    if ptr_valid(prev) then (*prev).next=next $
    else ptr=next
    ; destroy cut piece
    heap_free,current
    ; point to next
    current=next
    j++
endfor

end;pro DeleteLLEntree
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DestroyLLEntree,current,llistAl
if llistAl[0] ne -1 then $
    for i=0l,n_elements(llistAl)-1 do begin
       PTR_FREE,(*current).(llistAl[i])
    endfor
next=(*current).next
PTR_FREE,current
current=next
end;pro DestroyLLEntree
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DestroyLL,current,llistAl
;WHILE PTR_VALID(current) do DestroyLLEntree,current,llistAl

ptr=current
while ptr_valid(ptr) do BEGIN
    del=ptr
    ptr=(*ptr).next
    ptr_free,del
endwhile

; FATAL error when linked list is too long
;heap_free,current
end;pro DestroyLL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UpdateLL,files,list,delflag=delflag,delind=delind,nocheck=nocheck

delflag=0B

; ----Delete entrees----
i=0l
current=list.PTR_list
WHILE PTR_VALID(current) do begin
    ind=where(files eq (*current).str,count)
    if count eq 0 then begin ;Entree deleted from directory
       DestroyLLEntree,current,[-1];list.llistAl
       if PTR_VALID(last) then (*last).next=current $; Relink broken chain
       else list.PTR_list=current ; First entree deleted
       if delflag eq 0 then delind=i else begin
           delind=[delind,i]
        delflag=1B
       endelse
    endif else begin ;Entree already in directory before
       ind=ind[0]
       if n_elements(entree) eq 0 then entree=files[ind] $
       else entree=[entree,files[ind]]
       files=Shrinkarray(files,ind)
       last=current
       current=(*current).next
       if files[0] eq '' then begin
               WHILE PTR_VALID(current) do begin
                   last=current
                   current=(*current).next
               endwhile
       endif
    endelse
    i++
endwhile
if PTR_VALID(last) then current=last ;Pointer to last Entree

new_entree=files
if n_elements(entree) ne 0 then files=entree else files=''
if new_entree[0] eq '' then return,new_entree

; ----Add entrees----
n=n_elements(new_entree)
NotFinInd=0l
nocheck=keyword_set(nocheck)
for i=0l,n-1 do begin
    ; ----Remove entree if it is still been written----
    if nocheck then result=1B else result=CheckFile(new_entree[i])
    if result then begin
    if PTR_VALID(current) then begin ;current points to last entree (i.e. Not NULL)
       (*current).next = PTR_NEW(list.llist)
       current=(*current).next
       (*current).str=new_entree[i]
    endif else begin ;current points to NULL because no entrees present
              ; i.e. all deleted or no entrees in the beginning
       list.PTR_list= PTR_NEW(list.llist)
       (*list.PTR_list).str=new_entree[i]
       current=list.PTR_list
    endelse
    endif else NotFinInd=[NotFinInd,i]
endfor

if n_elements(NotFinInd) ne 1 then new_entree=Shrinkarray(new_entree,NotFinInd)

return,new_entree
end;function UpdateLL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GenericListInsertAfter,ptr,ptrnew

; ptr will point to the new entree

if PTR_VALID(ptrnew) then begin
    if PTR_VALID(ptr) then begin
        (*ptrnew).next=(*ptr).next
        if PTR_VALID((*ptr).next) then (*(*ptr).next).prev=ptrnew
        (*ptr).next=ptrnew
        (*ptrnew).prev=ptr
    endif
    ptr=ptrnew
endif

end;pro GenericListInsertAfter
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GenericListInsertBefore,ptr,ptrnew

; ptr will point to the new entree

if PTR_VALID(ptrnew) then begin
    if PTR_VALID(ptr) then begin
        (*ptrnew).next=ptr
        (*ptrnew).prev=(*ptr).prev
        (*ptr).prev=ptrnew
        if ptr_valid((*ptrnew).prev) then (*(*ptrnew).prev).next=ptrnew
    endif
    ptr=ptrnew
endif

end;pro GenericListInsertBefore
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GenericListDelete,ptr

; ptr will point to the previous entree or the next if there is no previous
if PTR_VALID(ptr) then begin
    ptrkeep=PTR_NEW()
    if PTR_VALID((*ptr).next) then begin
        ptrkeep=(*ptr).next
        (*(*ptr).next).prev=(*ptr).prev
    endif
    if PTR_VALID((*ptr).prev) then begin
        ptrkeep=(*ptr).prev
        (*(*ptr).prev).next=(*ptr).next
    endif
    (*ptr).prev=PTR_NEW()
    (*ptr).next=PTR_NEW()
    heap_free,ptr
    ptr=ptrkeep
endif

end;pro GenericListDelete
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GenericListDeleteI,ptr,ind

; ptr will point to the first entree when returned
if not ptr_valid(ptr) then return

ind=ind[sort(ind)]
j=0l
for i=0l,n_elements(ind)-1 do begin
    while j ne ind[i] do begin
        ptr=(*ptr).next
        j++
    endwhile

    ptrkeep=(*ptr).next
    GenericListDelete,ptr
    if ptr eq ptrkeep then j++
endfor

if ptr_valid(ptr) then $
while ptr_valid((*ptr).prev) do ptr=(*ptr).prev

end;pro GenericListDelete,ptr
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro GenericListDestroy,ptr
heap_free,ptr
end;pro GenericListDestroy
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GenericListSearch,ptrin,fieldname,fieldvalue

if not ptr_valid(ptrin) then return,ptr_new()
ffieldname=STRUPCASE(fieldname)
ptr=ptrin
bool=0b

repeat begin
    fields=TAG_NAMES(*ptr)
    ind=where(fields eq ffieldname,ct)
    if ct eq 1 then bool=(*ptr).(ind[0]) eq fieldvalue
    if not bool then ptr=(*ptr).next
endrep until (bool or (not PTR_VALID(ptr)))

return,ptr
end;function GenericListSearch
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
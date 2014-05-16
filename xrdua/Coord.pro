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

function RTDynExt,xy,sf,offset,x=x
invsf=1/sf
xy2=(xy-offset)*invsf
if invsf gt 1 then begin
    n=floor(invsf)
    nn=n*n
    ind=indgen(n)
    if keyword_set(x) then begin
        ind=reform(ind#replicate(1b,n),nn)
    endif else begin
        ind=rebin(ind,nn,/sample)
    endelse
    nxy2=n_elements(xy2)
    xy3=replicate(xy2[0],nxy2*nn)
    for i=0L,nxy2-1 do xy3[i*nn]=xy2[i]+ind
    return,xy3
endif else return,xy2
end;function RTDynExt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RTDyn,xy,sf,offset
; offset in real
xy2=(xy-offset+0.5)/sf-0.5
return,xy2
end;function RTDyn
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DynTR,xy,sf,offset
; offset in real
xy2=(xy+0.5)*sf-0.5+offset
return,xy2
end;function DynTR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RTStatExt,xy,sf,x=x
invsf=1/sf
xy2=xy*invsf
if invsf gt 1 then begin
    n=floor(invsf)
    nn=n*n
    ind=indgen(n)
    if keyword_set(x) then begin
        ind=reform(ind#replicate(1b,n),nn)
    endif else begin
        ind=rebin(ind,nn,/sample)
    endelse
    nxy2=n_elements(xy2)
    xy3=replicate(xy2[0],nxy2*nn)
    for i=0L,nxy2-1 do xy3[i*nn]=xy2[i]+ind
    return,xy3
endif else return,xy2
end;function RTStatExt
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RTStat,xy,sf
xy2=(xy+0.5)/sf-0.5
return,xy2
end;function RTStat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StatTR,xy,sf
xy2=(xy+0.5)*sf-0.5
return,xy2
end;function StatTR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ConvertFileNumber,X,Y,nCol,nRow,TRAN,REV1,REV2,REV3,XX=XX,YY=YY
; From an image nColxnRow the X and Y coordinates are given
; Return the corresponding coordinate in the default orientation

; Datablock was converted as follows: REV3 x TRAN x REV1 x REV2
; Remark: REV1 and REV2 are commutative

XX=X
YY=Y

if REV2 then YY=nRow-1-YY
if REV1 then XX=nCol-1-XX
if TRAN then begin
    t=XX
    XX=YY
    YY=t
    t=nCol
    nCol=nRow
    nRow=t
endif
if REV3 then begin
    ind=where((YY mod 2) eq 1)
    if ind[0] ne -1 then XX[ind]=nCol-1-XX[ind]
endif

return,XX+YY*nCol

end;function ConvertFileNumber
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CurPos,X,Y,MapCol,MapRow,sfh,sfv,XX=XX,YY=YY

XX=round(StatTR(X,sfh))
YY=(MapRow-1)-round(StatTR(Y,sfv))
return,YY*MapCol+XX

end;function CurPos
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
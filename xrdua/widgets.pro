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

function widgetparent,id,nup=nup
if keyword_set(nup) eq 0 then nup=1 else nup=long(nup)>0

idp=id
for i=1,nup do $
    if widget_info(idp,/valid) then $
        idp=widget_info(idp,/parent)

return,idp
end;function widgetparent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro DeleteWBaseChilds,base
; ----Delete all children----
ChildID=widget_info(base,/child)
while ChildID ne 0 do begin
    childID2 = Widget_Info(childID, /Sibling)
    widget_control,ChildID,/destroy
    childID = childID2
endwhile
end;pro DeleteWBaseChilds
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wheresibling,TopSearchIDs,ChildID,count=count

; get widget tree
tree=ChildID
pID=widget_info(ChildID,/PARENT)
while pID ne 0 do begin
    tree=[tree,pID]
    pID=widget_info(pID,/PARENT)
endwhile

ind=SetIntersection([TopSearchIDs],tree,count=count,/indices)

return,ind
end;function wheresibling
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function find_by_uname,TopSearchID,uname

; find first
ID=widget_info(TopSearchID,FIND_BY_UNAME=uname)
if ~widget_info(ID,/VALID_ID) then return,ID
uname2=uname+'2'
widget_control,ID,SET_UNAME=uname2

; find others
ID2=widget_info(TopSearchID,FIND_BY_UNAME=uname)
while widget_info(ID2,/VALID_ID) do begin
    ID=[ID,ID2]
    widget_control,ID2,SET_UNAME=uname2
    ID2=widget_info(TopSearchID,FIND_BY_UNAME=uname)
endwhile

; Reset unames
for i=0l,n_elements(ID)-1 do widget_control,ID[i],SET_UNAME=uname

return,ID
end;function find_by_uname
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function widgetsize,ID
geom=widget_info(long(ID),/GEOMETRY)
return,[geom.SCR_XSIZE + (2* geom.MARGIN),geom.SCR_YSIZE + (2* geom.MARGIN)]
end;function widgetsize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro widgetsetsize,ID,xs,ys
geom=widget_info(long(ID),/GEOMETRY)

CATCH, Error_status
IF Error_status ne 0 THEN return

widget_control,ID,scr_xsize=xs-2* geom.MARGIN,scr_ysize=ys-2* geom.MARGIN
end;pro widgetsetsize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function widgetTLB,ID

if n_elements(ID) eq 0 then return,0L
wID=ID
pID=widget_info(wID,/PARENT)
while pID ne 0 do begin
    wID=pID
    pID=widget_info(pID,/PARENT)
endwhile

return,wID
end;function widgetTLB
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function widgettreesize,TopSearchID,uname

wID=widget_info(TopSearchID,find_by_uname=uname)
S=[wID,widgetsize(wID)]
pID=widget_info(wID,/PARENT)
while pID ne 0 do begin
    S=[[S],[pID,widgetsize(pID)]]
    pID=widget_info(pID,/PARENT)
endwhile

return,S
end;function widgettreesize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function widgettreesizeinit,TopSearchID,uname
TLBsize=widgettreesize(TopSearchID,uname)
s=size(TLBsize)
if s[0] eq 2 then TLBsize=TLBsize[1:*,s[2]-1] $
else TLBsize=TLBsize[1:*]

return,TLBsize
end;function widgettreesizeinit
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro propagateresize,TopSearchID,TLBsize,uname,yxratio

; Get size of widget and all his parents
widgettree=widgettreesize(TopSearchID,uname)
n=n_elements(widgettree)/3

; How much did the top-level base change in size?
dx=widgettree[1,n-1]-TLBsize[0]
dy=widgettree[2,n-1]-TLBsize[1]

; Adapt resizing to conserve xy ratio
if n_elements(yxratio) eq 1 then $
    if abs(dx) gt abs(dy) then dy=dx*yxratio else dx=dy/yxratio
    
; New top-level size
TLBsize+=[dx,dy]
TLBsize>=0

; Resize all none-top-level widgets
widgettree[1,0:n-2]+=dx
widgettree[2,0:n-2]+=dy
for i=0l,n-2 do $
    widgetsetsize,widgettree[0,i],widgettree[1,i],widgettree[2,i]

; Resize top-level widget again
widgetsetsize,widgettree[0,n-1],TLBsize[0],TLBsize[1]

end;pro propagateresize
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro baseresizeevent,ev,topuname,drawuname
TopBase=find_by_uname(ev.top,topuname)
if ~widget_info(TopBase[0],/VALID_ID) then return
nTopBase=n_elements(TopBase)
for i=0l,nTopBase-1 do begin
     WIDGET_CONTROL, TopBase[i], GET_UVALUE=list, /NO_COPY
     TLBsize=list.TLBsize
     propagateresize,TopBase[i],TLBsize,drawuname,list.yxratio
     list.TLBsize=TLBsize
     WIDGET_CONTROL, TopBase[i], SET_UVALUE=list, /NO_COPY
endfor
end;pro baseresizeevent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
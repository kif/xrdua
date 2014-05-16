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

function Savgol1D_r2,y,w
;second degree pol
nr=fix(w)/2
nl=fix(w)/2
;----Evaluate input----
r=2
order=0
sy=size(y)
nr=nr>1
nl=nl>1
w=nr+nl+1
;----Set up Design Matrix----
ncol=r+1
nrow=w
s=indgen(nrow)-nl
one=fltarr(ncol)+1
A=one#s
for j=0l,r do A[j,*]^=j
;----Calculate Coefficients----
AT=transpose(A)
C=invert(AT##A)##AT
; each row of C corresponds with an a
; a = (a0 a1 a2 a3 ...)T
; corresponding filter  (f f' f'' ...)/order!
;----Extract desired coeff and multiplied by appropriate constant----
fil=C[*,order]
;----Output----
data=convol(y,fil,/edge_trunc)
return,data
end;Savgol1D_r2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Strip1D,y,w,niter,nlast
if n_elements(y) le w then return,y*0
yt=y>0
;----default values----
if n_elements(w) eq 0 then w=1
if n_elements(niter) eq 0 then niter=500
;if n_elements(nlast) eq 0 then nlast=8
;----transform (in this case square root)----
;yy=alog10(sqrt(yt))
yy=sqrt(yt)
;filter with width=w
yt=Savgol1D_r2(yy,w)
;yt=smooth(yy,w)
yb=yt
;----stripping----
;determine window
;ind0=w
;ind1=n_elements(y)-1-w
;wr=1.
;for n=1,niter do begin
;if n gt (niter-nlast) then wr=wr/sqrt(2)
;w=round(ind0*wr)
;for i=ind0,ind1 do yb[i]=((yt[i-w]+yt[i+w])/2)<yt[i]
;yt=yb
;endfor
h=fltarr(w)
h[0]=0.5
h[w-1]=0.5
for n=1l,niter do begin
    yb=convol(yb,h,/edge_trunc)
    yb<=yt
    yt=yb
endfor
;----inverse transform (in this case square)----
return,yb*yb
;return,(10^yb)^2
end;function Strip1D
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TopHatFilter1D, y, v, w
; w = FWHM
; v in [FWHM/3,FWHM/2]
; r in [2,4]
a=-v-w/2
b=-a
k=indgen(2*b+1)-b
h=dblarr(2*b+1)
ind1=fix((where(k eq (-w/2-1)))[0])
ind2=fix((where(k eq w/2))[0])
h[0:ind1]=-1./(2*v)
h[ind1+1:ind2]=1./w
h[ind2+1:*]=-1./(2*v)
h=transpose(h)
yfilter=y*0.
yfiltervar=yfilter
for i=0l,n_elements(y)-1 do begin
    yfilter[i]=h#y[i+k]
    yfiltervar[i]=(h^2)#y[i+k]
endfor
return, {yfilter:yfilter,yfiltervar:yfiltervar}
end; function TopHatFilter1D
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function estpeak,xval,spe,peakbegin,peakend,debug=debug,drawi=drawi

if peakbegin eq peakend then begin
    peakbegin-=5
    peakend+=5
    n=n_elements(xval)
    peakbegin=0>peakbegin<(n-1)
    peakend=0>peakend<(n-1)
endif

x=xval[peakbegin:peakend]
y=spe[peakbegin:peakend]
w = peakend-peakbegin+1

; ----First estimation: numerical integration----
; window= curvefit-window
;    -> estimation of peakhight and position by maximum y value
;    -> Trapezoidal Rule + linear background: area and sigma
temp=w/2
peaky=y[temp]
peakxkan=temp+peakbegin
peakxen=xval[peakxkan]
h=x-shift(x,1)
h=h[1:*]
f=y+shift(y,1)
f=f[1:*]
;peakarea=0.5*(total(h*f)-(xval[peakend]-xval[peakbegin])*(2*min([spe[peakend],spe[peakbegin]])))
peakarea=(0.5*(total(h*f)-(xval[peakend]-xval[peakbegin])*(spe[peakend]+spe[peakbegin])))>0
b1=(spe[peakend]-spe[peakbegin])/(xval[peakend]-xval[peakbegin])
b0=spe[peakend]-b1*xval[peakend]
hight=peaky-(b1*x[temp]+b0)
if hight le 0 then hight=peaky
sx=peakarea/(sqrt(2*!pi)*hight)

; ----Refine----
; window= view-window
;    ->estimate area, position and sx again
A=[peakarea,peakxen,sx,b0,b1]
Akeep=A
yfit=NLLS(x,y, A, 'gfunc2',CError=CError)
temp=abs((A-Akeep)/Akeep)
if CError gt 0 and temp[1] lt 0.5 then begin
    ; No error and peakshift not to high
    peakarea=A[0]
    peakxen=A[1]
    sx=A[2]
    b0=A[3]
    b1=A[4]
    hulpf=min(abs(xval-peakxen),peakxkan)
endif

if keyword_set(debug) then begin
    LoadctRev,0,/silent,rev=-1
    wset,drawi
    plot,x,y,psym=2,xstyle=1,ystyle=1
    gfunc2, x, Akeep, f
    oplot,x,f,linestyle=2
    gfunc2, x, A, f
    oplot,x,f
    wait,0.5
endif

return,[peakxen,peakxkan,peakarea,sx,peaky]
end; function estpeak
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function estspe,spe,xval,r,vin,win,debug=debug,drawi=drawi,details=details
vinn=vin
winn=win
rn=r

;----Calc sec deriv----
;fvar1=1
;fvar2=1
;spet=mSavgol(data=spe,nr=winn/2,nl=winn/2,order=2,fvar=fvar1,sigma2=sigma2)+$
;mSavgol(data=spe,nr=winn/2,nl=winn/2,order=2,fvar=fvar2,sigma2=sigma2)
;spet=-spet
;fvar=(fvar1+fvar2)
;spet=fvar>spet
fvar=1
spet=mSavgol(data=spe,nr=winn/2,nl=winn/2,order=2,fvar=fvar,sigma2=sigma2)
spet=-spet
spet=fvar>spet
;spet=spe

; ----Top-hat filter to spectrum----
;first param:v (width in channels of sidelobs)
;second param:w (width in channels of middlelob)
yf=TopHatFilter1D(spet,vinn,winn)
yfilter=yf.(0)
yfiltervar=yf.(1)
yfiltertemp=yfilter
ind=where(yfilter le rn*sqrt(abs(yfiltervar)),count)
if count ne 0 then yfilter[ind]=0.

; ----Loop over datapoints----
temp=0
peaky=0D
peakxen=0.
peakxkan=0L
sx=0.
peakarea=0.
ipeak=0L
epeak=0L
apeaks=0L
for i=1,n_elements(xval)-1 do begin
    if (yfilter[i-1] eq 0) and (yfilter[i] ne 0) then temp=i else $
    if (yfilter[i-1] ne 0) and (yfilter[i] eq 0) and (i-1 ne temp) then begin
       ; temp->index of begincut peak
       ; (i-1)->index of endcut peak

       ; ----Determine peakborders:----
       j=temp
       while (yfiltertemp[j] gt yfiltertemp[(j-1)>0]) do begin
               j=j-1
               if j eq 0 then break
       endwhile
       peakbegin=j
       j=i-1
       while yfiltertemp[j] gt yfiltertemp[j+1] do begin
           j=j+1
           if j eq (n_elements(yfilter)-1) then break
       endwhile
      peakend=j

       ; ----Prepare and save data----
       if apeaks ne 0 then begin
         peaky=[peaky,0]
         peakxen=[peakxen,0]
         peakxkan=[peakxkan,0]
         sx=[sx,0]
         peakarea=[peakarea,0]
         ipeak=[ipeak,0]
         epeak=[epeak,0]
       endif
       ipeak[apeaks]=peakbegin
       epeak[apeaks]=peakend

       ; ----Estimate parameters----
       peakval=estpeak(xval,spe,peakbegin,peakend,debug=debug,drawi=drawi)
       peakxen[apeaks]=peakval[0]
       peakxkan[apeaks]=peakval[1]
       peakarea[apeaks]=peakval[2]
       sx[apeaks]=peakval[3]
       peaky[apeaks]=peakval[4]

       apeaks=apeaks+1
    endif
endfor

if apeaks ne 0 then begin
if keyword_set(details) then out=fltarr(8,apeaks)$
    else out=fltarr(4,apeaks)
out[0,*]=peakarea
out[1,*]=peakxen
out[2,*]=sx
out[3,*]=0.5
if keyword_set(details) then begin
    out[4,*]=peakxkan
    out[5,*]=ipeak
    out[6,*]=epeak
    out[7,*]=peaky
endif
return,out
endif else return,0

end; pro estspe
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ortpol, xi, yi, weights=wi, r=r, COEFF=coeff, degree=degree, SIGMA=sigma,$
 POL=pol, GR=gr, maxdegree=maxdegree, print=print
;input: xi en yi
;   no weights value => weights=1/yi
;   no r value (used in location of background datapoints) => 2
;   no degree value => degree=4 and counting, autonomous determined degree=GR
;output: backgroupsmodel
;   coeff (with SD=SIGMA) and pol


;examine input
if n_elements(maxdegree) eq 0 then maxdegree=30
if n_elements(degree) eq 0 then begin
    degree=4
    superv=0
endif else begin
    if degree lt 0 then begin
        degree=4
        superv=0
    endif else begin
        superv=1
    endelse
endelse
nxi = n_elements(xi)
xxi = double(xi)
yyi = double(yi)
xmed = 0.
avx = total(xxi)/nxi
xmed = median(xxi)
xxi -= avx; shift x to the middle
if n_elements(r) eq 0 then r=0.2


itt=0
if(n_elements(weights) le 0) then wi = 1./(yi>1.)
wwi = double(wi)
nbpnts=nxi
;loop for succesively improve coeff by adjusting the weights
;until the stop changing
repeat begin
    ;variables
    P=dblarr(nxi,degree+1) ;each row is a pol:P0,P1,...
    a=dblarr(degree+1)
    b=a
    g=a
    C=a
    sigma=a
    ;j=0 -> P0,C0
    P[*,0]=1.
    g[0]=total(wwi*P[*,0]*P[*,0])
    a[0]=(total(wwi*xxi*P[*,0]*P[*,0]))/g[0]
    b[0]=0.
    C[0]=(total(wwi*yyi*P[*,0]))/g[0]
    ;j=1 -> P1,C1
    if degree ge 1 then begin
        P[*,1]=(xxi-a[0])*P[*,0]
        g[1]=total(wwi*P[*,1]*P[*,1])
        a[1]=(total(wwi*xxi*P[*,1]*P[*,1]))/g[1]
        b[1]=g[1]/g[0];
        C[1]=(total(wwi*yyi*P[*,1]))/g[1]
    endif
    ;j -> Pj,Cj
    for j=2,degree do begin
        P[*,j]=(xxi-a[j-1])*P[*,j-1]-b[j-1]*P[*,j-2]
        g[j]=total(wwi*P[*,j]*P[*,j])
        a[j]=(total(wwi*xxi*P[*,j]*P[*,j]))/g[j]
        b[j]=g[j]/g[j-1]
        C[j]=(total(wwi*yyi*P[*,j]))/g[j]
    endfor
    sigma=(g*0.+1.)/sqrt(abs(g))
    yfit=C##P
    chisq=total(wwi*(yyi-yfit)^2)
    ;reduced chisq(divided by degree of freedom)=rchisq => 1 for perfect fit
    if nbpnts gt (degree+1) then rchisq=chisq/double(nbpnts-degree-1) else rchisq=1.
    s=sqrt(rchisq)

    ;if superv=0: add polynomials until the coeff of the added polynomial statisticaly doesn't differ from 0
    while (abs(C[degree]) ge abs(2*s*sigma[degree])) and (degree le maxdegree) and (superv eq 0) do begin
     degree++
     GR=degree
     hoopP=P
     hoopg=g
     hoopa=a
     hoopb=b
     hoopC=C
     hoops=sigma
     P=dblarr(nxi,degree+1) ;each row is a pol:P0,P1,...
     a=dblarr(degree+1)
     b=a
     g=a
     C=a
     sigma=a
     P[*,0:degree-1]=hoopP
     g[0:degree-1]=hoopg
     a[0:degree-1]=hoopa
     b[0:degree-1]=hoopb
     C[0:degree-1]=hoopC
     sigma[0:degree-1]=hoops
     P[*,degree]=(xxi-a[degree-1])*P[*,degree-1]-b[degree-1]*P[*,degree-2]
     g[degree]=total(wwi*P[*,degree]*P[*,degree])
     a[degree]=(total(wwi*xxi*P[*,degree]*P[*,degree]))/g[degree]
     b[degree]=g[degree]/g[degree-1]
     C[degree]=(total(wwi*yyi*P[*,degree]))/g[degree]
     sigma[degree]=1/sqrt(abs(g[degree]))
     yfit=C##P
     chisq=total(wwi*(yyi-yfit)^2)
     if nbpnts gt (degree+1) then rchisq=chisq/double(nbpnts-degree-1) else rchisq=1.
     s=sqrt(rchisq)
    endwhile
    superv=1
    coeff=C
    pol=P
    ;we performed an other iteration
    itt++

    ;test whether to perform another weight adjustment (iteration) or not
    ;stop if difference between old and new coeff is smaller the SD on new coeff * sqrt(reduced chisq)
    crit=1; =stop
    if itt ge 2 then begin
     for j=0l,degree do begin
         if abs(coeff[j]-coeffold[j]) gt (s*sigma[j]) then begin
            crit=0; =perform an other iteration
            goto ,jump
         endif
     endfor
    endif else crit=0
    jump:crit=crit
    ;stop when max iterations
    if itt ge 20 then crit=1

    ;weights for next itertion
    ;wwi[i]=1/yfit[i] if datapoint belong to background
    ;wwi[i]=1/(yfit[i]-yyi[i])^2 of 0 if datapoint DOESN'T belong to background
    nbpnts=0
    for i=0l,nxi-1 do begin
     if yyi[i] le yfit[i]+r*sqrt(abs(yfit[i])) then begin
         wwi[i]=1/(yfit[i]>1.)
         nbpnts=nbpnts+1
     endif else wwi[i]=0;1/(yfit[i]-yyi[i])^2 ;wwi[i]=0
    endfor

    ;save coeff form determination of end of loop
    coeffold=coeff
endrep until crit ; repeat until crit=1

if keyword_set(print) then begin
    printw,print, 'Number of iterations for polynomial of degree '+strtrim(string(degree),2)+':'+string(itt)+' (max=20)'
endif
GR=degree
return, yfit
end; function ortpol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fitod,xval,spe,list,$
    chisqr=chisqr,debug=debug,drawi=drawi,interactive=interactive

; ----Possible modes----
;    Interactive or in batch mode: keyword_set(interactive)
;    Estimate peaks before or after background subtraction: list.filorder
;    Background: substract nothing (-1), substract ortpol (0), fit ortpol (1), substract strip (2)
;    Fit:     MLR, G, substr
;            NLLS, G, fit ortpol
;            NLLS, G, substr
;            NLLS, PV(n), substr
;            NLLS, PV(narr),substr

x=xval
y=spe

; ----Background already substracted?----
if not keyword_set(interactive) then backcor=0 else begin
    backcor=list.backcor
    ; If background is already substracted, then we
    ; can't perform some of the fitting properties
    bool1=(list.fitmethod eq 1) and (backcor eq 1)
    bool2=(list.filorder eq 0) and (backcor eq 1)
    if bool1 or bool2 then begin
       printw,interactive,'Restore spectrum first'
       return,0
    endif
endelse

; ----Fitmethod and backmod conflict?----
if (list.fitmethod eq 1) then list.backmod=1
if (list.backmod eq 1) and (list.fitmethod ne 1) then list.backmod=0

; ----Estimate peaks in ROI before subtraction----
; A,m,s,n
if list.filorder eq 0 then $
if keyword_set(interactive) then begin
    if list.apeaks eq 0 then begin
       printw,interactive,'No peaks present: search for them'
       peaks,list,out=peaks,apeaks=apeaks,x=x,y=y,$
         offx=list.backranind[0],offy=list.backranind[0]
    endif else begin
       ;----Parameters of detected peaks in ROI:----
       ;temp=list.peakxkan[0:list.apeaks-1]-list.backranind[0]
       ;hulpf=min(abs(temp),index1)
       ;if temp[index1] lt 0 then index1=index1+1
       ;temp=list.peakxkan[0:list.apeaks-1]-list.backranind[1]
       ;hulpf=min(abs(list.peakxkan-list.backranind[1]),index2)
       ;if temp[index2] gt 0 then index2=index2-1
       sx=list.peaksigma[0:list.apeaks-1]
       positions=list.peakxen[0:list.apeaks-1]
       areas=list.peakarea[0:list.apeaks-1]
       apeaks=n_elements(positions)
       peaks=transpose([[areas],[positions],[sx],[intarr(apeaks)+0.5]])
    endelse
    if apeaks eq 0 then begin
        printw,interactive,'No peaks present in ROI'
        return,0
    endif
endif else begin
    peaks=estspe(y,x,list.r,list.vin,list.w,debug=debug,drawi=drawi)
    s=size(peaks)
    if s[1] ne 4 then return,0 ;No peaks found
    if s[0] eq 1 then apeaks=1 else apeaks=s[2]
endelse

; ----Calculate 1D background and substract----
if backcor eq 0 then begin
    case list.backmod of
       -1:   yb=y*0
       2:    begin
             yb=Strip1D(y,(*(*list.ptrModel).ptrFit).backparamIO[2],(*(*list.ptrModel).ptrFit).backparamIO[1])
             BackSub,y,yb,list.backmod,offset=offset
             if n_elements(offset) ne 0 then yb=(yb+offset)>0
             endcase
       0:    begin
             yb=ortpol(x,y,COEFF=COEFF,POL=POL,SIGMA=SIGMA,degree=(*(*list.ptrModel).ptrFit).backparamIO[0],GR=GR,print=interactive)
             BackSub,y,yb,list.backmod,offset=offset
             if n_elements(offset) ne 0 then yb=(yb+offset)>0
             endcase
       1:    begin
             yb=ortpol(x,y,COEFF=COEFF,POL=POL,SIGMA=SIGMA,degree=(*(*list.ptrModel).ptrFit).backparamIO[0],GR=GR,print=interactive)
             endcase
             ;POL: each row is a polynomial
    endcase
endif else yb=y*0
weights=1./(y>1.)

; ----Estimate peaks in ROI after subtraction----
; A,m,s,n
if list.filorder eq 1 then $
if keyword_set(interactive) then begin
    printw,interactive,'Peaksearch after background calc.'
    peaks,list,out=peaks,apeaks=apeaks,x=x,y=y,$
       offx=list.backranind[0],offy=list.backranind[0]
    if apeaks eq 0 then begin
       printw,interactive,'No peaks found'
       return,0
    endif
endif else begin
    peaks=estspe(y,x,list.r,list.vin,list.w,debug=debug,drawi=drawi)
    s=size(peaks)
    if s[1] ne 4 then return,0
    if s[0] eq 1 then apeaks=1 else apeaks=s[2]
endelse

; ----Status description----
; y: spectrum to fit (background substracted if necessary)
; yb: background (not used further when substracted)
; peaks: estimated (before or after subtraction) peak parameters (PV, i.e. 4 per peak)

; ----Choose between different fit procedures----
case list.fitmethod of
0:  begin
    ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ;multiple linear regression fit:
    ;-------------------------------
    ;   background  -> substracted
    ;   parameters  -> area of each detected Gaussian peak
    ;   output      -> new areas and backgroundcoeff (nothing saved!!!!)
    ;   constraints -> none
    ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    areas=peaks[0,*]
    positions=peaks[1,*]
    sx=peaks[2,*]

    XX=dblarr(apeaks,n_elements(x))
    for k=0l,apeaks-1 do XX[k,*]=1/(sqrt(2*!dpi)*sx[k])*exp(-(X-positions[k])^2/(2*sx[k]^2))
    result=regress(XX,Y,weights,yfit,const,sigma,ftest,r,rmul,chisq)
    positions2=positions
    areas2=result[0:apeaks-1]
    sx2=sx
    n2=0
    spositions2=0
    sareas2=sigma[0:apeaks-1]
    ssx2=0
    sn2=0
    ITER=0

    if keyword_set(interactive) then begin
       ; Save fit
       if not ptr_valid(list.yfit) then list.yfit=ptr_new((*list.spet)[*,list.xtype]*0)
       if not ptr_valid(list.yback) then list.yback=ptr_new((*list.spet)[*,list.xtype]*0)
       (*list.yfit)[list.backranind[0]:list.backranind[1]]=yfit+yb
       (*list.yback)[list.backranind[0]:list.backranind[1]]=yb+const
       ; Print peakdata
       printw,interactive,'Old areas:|New areas:|Stdev:'
       m=transpose([[reform(areas)],[areas2],[sareas2]])
       printw,interactive,m
    endif

    endcase
1:  begin
    ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ;nonlinear least-squares fit1:
    ;---------------------------
    ;   background -> ortpol
    ;   parameters  -> area,FWHM en position of each detected Gaussian peak + coeff of ortpol
    ;   output      -> new areas,FWHMs, positions en backgroundcoeff (nothing saved!!!!)
    ;   constraints -> none
    ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    areas=peaks[0,*]
    positions=peaks[1,*]
    sx=peaks[2,*]
    A=dblarr(apeaks*3+GR+1)
    indar=indgen(apeaks)*3
    indpos=indar+1
    indsx=indar+2
    indcoeff=indgen(GR+1)+3*apeaks
    A[indar]=areas
    A[indpos]=positions
    A[indsx]=sx
    A[indcoeff]=COEFF
    listnew={POL:POL,GR:GR}

    ; View init. guess
    if keyword_set(interactive) then begin
       loadct,39,/silent
       gfunc1, x, a, yfit,listnew;->initial fit
       wset,list.drawindex
       oplot,x,yfit,color=15*16
    endif

    yfit = marquardt(x,y, weights, A, SIGMA,FUNCTION_NAME='gfunc1',CHISQ=CHISQ,$
       ITER=ITER,itmax=100,struc=listnew,nfree=nfree)

    areas2=A[indar]
    positions2=A[indpos]
    sx2=A[indsx]
    n2=0
    sareas2=sigma[indar]
    spositions2=sigma[indpos]
    ssx2=sigma[indsx]
    sn2=0
    coeff2=A[indcoeff]

    if keyword_set(interactive) then begin
       printw,interactive,'Number of iterations performed by Non-Linear Least Squares fit:'+string(iter)
       if not ptr_valid(list.yfit) then list.yfit=ptr_new((*list.spet)[*,list.xtype]*0)
       if not ptr_valid(list.yback) then list.yback=ptr_new((*list.spet)[*,list.xtype]*0)
       (*list.yfit)[list.backranind[0]:list.backranind[1]]=yfit
       (*list.yback)[list.backranind[0]:list.backranind[1]]=coeff2##POL
       ; Print peakdata
       printw,interactive,'Old positions:|New positions:|Stdev:'
       m=transpose([[reform(positions)],[positions2],[spositions2]])
       printw,interactive,m
       printw,interactive,'Old areas:|New areas:|Stdev:'
       m=transpose([[reform(areas)],[areas2],[sareas2]])
       printw,interactive,m
       printw,interactive,'Old sx:|New sx:|Stdev:'
       m=transpose([[reform(sx)],[sx2],[ssx2]])
       printw,interactive,m
       printw,interactive,'Old coeff:|New coeff:|Stdev:'
       m=transpose([[reform(coeff)],[coeff2],[sigma[indcoeff]]])
       printw,interactive,m
    endif

    endcase
2:  begin
    ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ;nonlinear least-squares fit2:

    ;---------------------------
    ;   background  -> substracted
    ;   parameters  -> area,FWHM and position of each detected Gaussian peak
    ;   output      -> new areas,FWHMs, positions (nothing saved!!!!)
    ;   constraints -> suggestive
    ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    areas=peaks[0,*]
    positions=peaks[1,*]
    sx=peaks[2,*]
    A=dblarr(apeaks*3)
    indar=indgen(apeaks)*3
    indpos=indar+1
    indsx=indar+2
    A[indar]=areas
    A[indpos]=positions
    A[indsx]=sx
    constr=dblarr(3*apeaks)+1
    constr[indar]=100
    constr[indpos]=0.2
    constr[indsx]=10

    ; View init. guess
    if keyword_set(interactive) then begin
       loadct,39,/silent
       gfunc3, x, a, yfit;->initial fit
       wset,list.drawindex
       oplot,x,yfit+yb,color=15*16
    endif

    ;yfit = marquardtc(x,y, weights, A, SIGMA,FUNCTION_NAME='gfunc3',$
    ;   CHISQ=CHISQ,ITER=ITER,itmax=100,nfree=nfree,constr=constr)
    yfit = marquardt(x,y, weights, A, SIGMA,FUNCTION_NAME='gfunc3',$
       CHISQ=CHISQ,ITER=ITER,itmax=100,nfree=nfree)

    areas2=A[indar]
    positions2=A[indpos]
    sx2=A[indsx]
    n2=0
    sareas2=sigma[indar]
    spositions2=sigma[indpos]
    ssx2=sigma[indsx]
    sn2=0

    if keyword_set(interactive) then begin
       printw,interactive,'Number of iterations performed by Non-Linear Least Squares fit:'+string(iter)
       if not ptr_valid(list.yfit) then list.yfit=ptr_new((*list.spet)[*,list.xtype]*0)
       if not ptr_valid(list.yback) then list.yback=ptr_new((*list.spet)[*,list.xtype]*0)
       (*list.yfit)[list.backranind[0]:list.backranind[1]]=yfit+yb
       (*list.yback)[list.backranind[0]:list.backranind[1]]=yb
       ; Print peakdata
       printw,interactive,'Old positions:|New positions:|Stdev:'
       m=transpose([[reform(positions)],[positions2],[spositions2]])
       printw,interactive,m
       printw,interactive,'Old areas:|New areas:|Stdev:'
       m=transpose([[reform(areas)],[areas2],[sareas2]])
       printw,interactive,m
       printw,interactive,'Old sx:|New sx:|Stdev:'
       m=transpose([[reform(sx)],[sx2],[ssx2]])
       printw,interactive,m
    endif

    endcase
3:  begin
    ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ;nonlinear least-squares fit3:

    ;---------------------------
    ;   background  -> substracted
    ;   parameters  -> area,FWHM and position of each detected PV peak + n
    ;   output      -> new areas,FWHMs, positions (nothing saved!!!!)
    ;   constraints -> suggestive
    ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    areas=peaks[0,*]
    positions=peaks[1,*]
    sx=peaks[2,*]
    n=peaks[3,0]
    A=dblarr(apeaks*3)
    indar=indgen(apeaks)*3
    indpos=indar+1
    indsx=indar+2
    A[indar]=areas
    A[indpos]=positions
    A[indsx]=sx
    constra=dblarr(4*apeaks)+1
    constra[indar]=100
    constra[indpos]=0.2
    constra[indsx]=10
    constra[3*apeaks:*]=0.5
    A=[A,n]

    ; View init. guess
    if keyword_set(interactive) then begin
       loadct,39,/silent
       gfunc4, x, a, yfit;->initial fit
       wset,list.drawindex
       oplot,x,yfit+yb,color=15*16
    endif

    yfit = marquardtc(x,y, weights, A, SIGMA,FUNCTION_NAME='gfunc4',CHISQ=CHISQ,$
       ITER=ITER,itmax=100,nfree=nfree,constr=constra)

    areas2=A[indar]
    positions2=A[indpos]
    sx2=A[indsx]
    n2=A[n_elements(A)-1]
    sareas2=sigma[indar]
    spositions2=sigma[indpos]
    ssx2=sigma[indsx]
    sn2=sigma[n_elements(sigma)-1]

    if keyword_set(interactive) then begin
       printw,interactive,'Number of iterations performed by Non-Linear Least Squares fit:'+string(iter)
       if not ptr_valid(list.yfit) then list.yfit=ptr_new((*list.spet)[*,list.xtype]*0)
       if not ptr_valid(list.yback) then list.yback=ptr_new((*list.spet)[*,list.xtype]*0)
       (*list.yfit)[list.backranind[0]:list.backranind[1]]=yfit+yb
       (*list.yback)[list.backranind[0]:list.backranind[1]]=yb
       ; Print peakdata
       printw,interactive,'Old positions:|New positions:|Stdev:'
       m=transpose([[reform(positions)],[positions2],[spositions2]])
       printw,interactive,m
       printw,interactive,'Old areas:|New areas:|Stdev:'
       m=transpose([[reform(areas)],[areas2],[sareas2]])
       printw,interactive,m
       printw,interactive,'Old sx:|New sx:|Stdev:'
       m=transpose([[reform(sx)],[sx2],[ssx2]])
       printw,interactive,m
       printw,interactive,'% Lorentz old:|% Lorentz new|Stdev:'
       printw,interactive,100*[n,n2,sn2]
    endif
    endcase
5:  begin
    ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    ;nonlinear least-squares fit5:

    ;---------------------------
    ;   background  -> substracted
    ;   parameters  -> area,sx,position and n of each detected PV peak
    ;   output      -> new areas,sxs, positions (nothing saved!!!!)
    ;   constraints -> real
    ;@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    constrar=fltarr(8,apeaks)
    constrar[0,*]=peaks[0,*]
    constrar[1,*]=total(spe)
    constrar[2,*]=peaks[1,*]
    constrar[3,*]=2
    constrar[4,*]=peaks[2,*]
    constrar[5,*]=peaks[2,*]
    constrar[6,*]=0.5
    constrar[7,*]=0.5

    ; physical meaningful parameters AA
    AA=peaks
    ; calculate fit parameters A
    A=AA;=> physical meaningful parameters
    ;ar0,Car,m0,Cm,F0,Cf,n0,Cn
    j=bindgen(4)
    for i=0l,n_elements(j)-1 do begin
       A[j[i],*]=tan((!dpi/2)*(A[j[i],*]-constrar[2*j[i],*])/constrar[2*j[i]+1,*])
    endfor
    A=reform(A,4*apeaks);=> fit parameters
    struc={constrar:constrar,j:j}

    ; View init. guess
    if keyword_set(interactive) then begin
       loadct,39,/silent
       gfunc6c, x, a, yfit,struc;->initi�le fit
       wset,list.drawindex
       oplot,x,yfit+yb,color=15*16
    endif

    yfit = marquardt(x,y, weights, A, SIGMA,FUNCTION_NAME='gfunc6c',CHISQ=CHISQ,$
                 ITER=ITER,itmax=100,nfree=nfree,struc=struc)

    A=reform(A,4,apeaks)
    AA=A
    sigma=reform(sigma,4,apeaks)
    sigma2=sigma^2.
    for i=0l,n_elements(j)-1 do begin
       AA[j[i],*]=constrar[2*j[i],*]+2.*constrar[2*j[i]+1,*]/!dpi*atan(A[j[i],*])
       sigma2[j[i],*]=sigma2[j[i],*]*(2/!dpi*constrar[2*j[i]+1,*]/(1+A[j[i],*]^2))^2
    endfor
    sigma=sqrt(sigma2)

    areas2=reform(AA[0,*])
    positions2=reform(AA[1,*])
    sx2=reform(AA[2,*])
    n2=reform(AA[3,*])
    sareas2=reform(sigma[0,*])
    spositions2=reform(sigma[1,*])
    ssx2=reform(sigma[2,*])
    sn2=reform(sigma[3,*])

    ; ----Reduce fitting results----
    ;red_res,chisq,areas2,positions2,sx2,n2,sigma

    if keyword_set(interactive) then begin
        if not ptr_valid(list.yfit) then list.yfit=ptr_new((*list.spet)[*,list.xtype]*0)
       if not ptr_valid(list.yback) then list.yback=ptr_new((*list.spet)[*,list.xtype]*0)
       (*list.yfit)[list.backranind[0]:list.backranind[1]]=yfit+yb
       (*list.yback)[list.backranind[0]:list.backranind[1]]=yb
       printw,interactive,'Number of iterations performed by Non-Linear Least Squares fit:'+string(iter)
       ; Print peakdata
       printw,interactive,'Old positions:|New positions:|Stdev:'
       m=transpose([[reform(peaks[1,*])],[positions2],[spositions2]])
       printw,interactive,m
       printw,interactive,'Old areas:|New areas:|Stdev:'
       m=transpose([[reform(peaks[0,*])],[areas2],[sareas2]])
       printw,interactive,m
       printw,interactive,'Old sx:|New sx:|Stdev:'
       m=transpose([[reform(peaks[2,*])],[sx2],[ssx2]])
       printw,interactive,m
       printw,interactive,'% Lorentz old:|% Lorentz new:|Stdev:'
       m=transpose([[reform(peaks[3,*])],[n2],[sn2]])
       printw,interactive,100*m
    endif
    endcase
endcase

; ----Save Data in interactive mode----
if keyword_set(interactive) then begin
    list.peakout=list.peakout*0
    list.peakouts=list.peakouts*0
    ; Save peakdata
    list.peakout[0,0:apeaks-1]=positions2
    list.peakout[1,0:apeaks-1]=areas2
    list.peakout[2,0:apeaks-1]=sx2
    list.peakout[3,0:apeaks-1]=n2
    list.peakouts[0,0:apeaks-1]=spositions2
    list.peakouts[1,0:apeaks-1]=sareas2
    list.peakouts[2,0:apeaks-1]=ssx2
    list.peakouts[3,0:apeaks-1]=sn2
    ; Save fit
    if not ptr_valid(list.resid) then list.resid=ptr_new((*list.spet)[*,list.xtype]*0)
    (*list.resid)[list.backranind[0]:list.backranind[1]]=((*list.spet)[list.backranind[0]:list.backranind[1],list.xtype]- $
    (*list.yfit)[list.backranind[0]:list.backranind[1]])/sqrt((*list.spet)[list.backranind[0]:list.backranind[1],list.xtype])
    list.chisq=chisq
    list.iter=ITER
    list.apeaksout=apeaks
endif

; ----Visualize----
if keyword_set(debug) then begin
    LoadctRev,0,/silent,rev=-1
    wset,drawi
    plot,x,y+yb,linestyle=1,title="Chisq red:"+string(chisq)
    oplot,x,yfit+yb
    oplot,x,yb,linestyle=2,COLOR=16*7
    for i=0l,apeaks-1 do $
    PlotS, [positions2[i],positions2[i]], $
       [!Y.CRange[0], !Y.CRange[1]],linestyle=2,NOCLIP = 0,color=16*5
    wait,0.5
endif

chisqr=chisq
out=areas2
return,out
end; function fitod
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
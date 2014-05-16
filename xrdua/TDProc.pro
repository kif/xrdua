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

function image_equal,img1,img2,npix=npix,shifttol=shifttol,Rtol=Rtol,outid=outid,xsize=xsize,ysize=ysize,quick=quick

; Image offsets or scales don't matter
; npix: subimage pixels for cross-correlation
; shifttol: subimage shift tollerance
; Rtol: cross-correlation tollerance
;
; nmag: for display
; outid: output widget

s1=dimsize(img1,2)
s2=dimsize(img2,2)
msizex=s1[0]<s2[0]
msizey=s2[1]<s2[1]
msize=msizex<msizey

; ----------
; Equal size
; ----------
img1_resized=img1
img2_resized=img2
if ~array_equal(s1,s2) then begin
    s1_resized=s1
    s2_resized=s2
    if s1[0] gt s2[0] then s2_resized[0]=s1[0] $
    else s1_resized[0]=s2[0]
    if s1[1] gt s2[1] then s2_resized[1]=s1[1] $
    else s1_resized[1]=s2[1]
    if ~array_equal(s1,s1_resized) then img1_resized=congrid(img1,s1_resized[0],s1_resized[1],/interp)
    if ~array_equal(s2,s2_resized) then img2_resized=congrid(img2,s2_resized[0],s2_resized[1],/interp)
endif
n_resized=n_elements(img1_resized)
mi=min(img1_resized,max=ma)
img1_resized_scl=bytscl(img1_resized,min=mi,max=ma)
mi=min(img2_resized,max=ma)
img2_resized_scl=bytscl(img2_resized,min=mi,max=ma)

; ----------------------------------------------
; Pearson product-moment correlation coefficient
; ----------------------------------------------
r_pearson = pearson_correlate(img1,img2)

;----------------------------------------
; Spearman's rank correlation coefficient
;----------------------------------------
r_spearman = R_CORRELATE(img1_resized,img2_resized)

;----------------------------------------
; Kendalls's (tau) rank correlation
;----------------------------------------
if n_resized le 32767 and ~keyword_set(quick) then $
    r_kendall = R_CORRELATE(img1_resized,img2_resized,/kendall) $
else r_kendall=fltarr(2)

;----------------------------------------
; RMS difference
;----------------------------------------
diff=img1_resized_scl-img2_resized_scl
RMS=sqrt(total(diff^2.)/n_resized)

;----------------------------------------
; Normalized difference
;----------------------------------------
; 1: equal
; 0: complement
SIMDIFF=1-total(abs(diff))/(255.*n_resized)

;-----------------------------
; Normalized Cross-correlation
;-----------------------------
r_cross=c_correlate(img1[0:msizex-1,0:msizey-1],img2[0:msizex-1,0:msizey-1],0)
r_cross=[r_cross,0] ; no p-value calculation

;----------
; Subimages
;----------
if not keyword_set(npix) then npix=fix(msize*0.4)>10 else npix=npix[0] ; 40% of the size
npix=npix<[msizex,msizey]
if not keyword_set(nmag) then nmag=1
if n_elements(shifttol) eq 0 then shifttol=(msize*0.01)>1 else shifttol=shifttol[0] ; 1% of the size
if not keyword_set(Rtol) then Rtol=0.9 else Rtol=Rtol[0]

; Subimages in img2
nsub=s2/npix
nx=nsub[0]
ny=nsub[1]
x0=npix[0]*indgen(nx)
if nx gt 1 then x1=[x0[1:*],s2[0]]-1 else x1=s2[0]-1
y0=npix[1]*indgen(ny)
if ny gt 1 then y1=[y0[1:*],s2[1]]-1 else y1=s2[1]-1

; img2 subimages in img1
xoff=lonarr(nsub)
yoff=xoff
xyccor=fltarr(nsub)

; Cross-correlate subimages of img2 with img1
for i=0l,nx-1 do $
    for j=0l,ny-1 do begin
        ; Extract subimage
        sub=img2[x0[i]:x1[i],y0[j]:y1[j]]

        ; Number of sub-shifts in img1
        nsub=n_elements(sub)
        ssub=dimsize(sub,2)-1
        noffx=s1[0]-ssub[0]
        noffy=s1[1]-ssub[1]
        ccor=fltarr(noffx,noffy)
        
        ; Correlate sub with img1
        for k=0l,noffx-1 do $
            for l=0l,noffy-1 do $
                ccor[k,l]=c_correlate(img1[k:k+ssub[0],l:l+ssub[1]],sub,0)
                ;ccor[k,l]=mutual_information(img1[k:k+ssub[0],l:l+ssub[1]],sub,/norm,/scale)
                ;ccor[k,l]=total((img1[k:k+ssub[0],l:l+ssub[1]]-sub)^2.)/nsub

        ; Sub image offset and cross-correlation
        mccor=max(ccor,moff)
        xoff[i,j]=moff mod noffx
        yoff[i,j]=moff/noffx
        xyccor[i,j]=mccor
        
        if keyword_set(outid) then $
            printw,outid,'Progress: '+stringr((i*ny+j+1.)/(nx*ny)*100)+'%'
    endfor

; Check whether img2 and img1 are equal
bool = xyccor ge Rtol
if nx eq 1 and ny eq 1 then tmp=0. else $
tmp=rebin(total(xoff,2)/ny,nx,ny)-xoff
mashiftx=max(tmp,abs,min=mishiftx)
bool and= tmp le shifttol
tmp=rebin(reform(total(yoff,1),1,ny)/nx,nx,ny)-yoff
mashifty=max(tmp,abs,min=mishifty)
bool and= tmp le shifttol
bsame=total(~bool,/pres) eq 0

if keyword_set(outid) then begin
    img2recon=img1*0
    for i=0l,nx-1 do $
        for j=0l,ny-1 do $
            if bool[i,j] then img2recon[xoff[i,j],yoff[i,j]]=img2[x0[i]:x1[i],y0[j]:y1[j]]

    window,title='Images are '+(bsame?'':'not ')+'the same.'
    
    tvscl,congrid(img1,xsize,ysize),0,/order
    tvscl,congrid(img2,xsize,ysize),1,/order
    tvscl,congrid(img2recon,xsize,ysize),2,/order
    
    s=dimsize(diff,2)
    tv,congrid(abs(diff),xsize,ysize),3,/order

    xyouts,0.1,0.1,'img1',/normal,color=100
    xyouts,0.3,0.1,'img2',/normal,color=100
    xyouts,0.5,0.1,'reconstructed img2',/normal,color=100
    xyouts,0.7,0.1,'img1 - img2',/normal,color=100
    
    printw,outid,'Cross-correlation tollerance: '+string(Rtol)
    printw,outid,'Subimage shift tollerance: '+string(shifttol)
    printw,outid,'Subimage size (pixels): '+string(npix)
    mi=min(xyccor,max=ma)
    printw,outid,'Cross-correlation of subimages between: '+string(mi,ma,format='(f0.2," - ",f0.2)')
    printw,outid,'X-shift of subimages: '+string(abs(mishiftx),abs(mashiftx),format='(f0.2," - ",f0.2)')
    printw,outid,'Y-shift of subimages: '+string(abs(mishifty),abs(mashifty),format='(f0.2," - ",f0.2)')
    printw,outid,'Total cross-correlation (p-value):'+string(r_cross,format='(f," (",f0,")")')
    printw,outid,'Pearson R (p-value):'+string(r_pearson,format='(f," (",f0,")")')
    printw,outid,'Spearman R (p-value):'+string(r_spearman,format='(f," (",f0,")")')
    printw,outid,'Kendall tau (p-value):'+string(r_kendall,format='(f," (",f0,")")')
    printw,outid,'RMS difference:'+string(RMS)
    printw,outid,'Sim measure (1=equal,0=complement):'+string(SIMDIFF)
    ;isurface,xyccor
endif

return,{bpass:bsame,r_pearson:r_pearson,r_spearman:r_spearman,r_kendall:r_kendall,$
        r_cross:r_cross,Rtol:Rtol,shifttol:shifttol,npix:npix,RMS:RMS,SIMDIFF:SIMDIFF}
end;function image_equal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Set2DbkCode,arr
; arr=[type,lower thresholding?,higher thresholding?]

b1=arr[0]
b2=ishft(long(arr[1] or 2*arr[2]),8)
return,b1 or b2

end;pro Set2DbkCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Get2DbkCode,in32

; first byte: 0(nothing),1(strip),2(image)
; second byte: 0(nothing),1(lower thresholding),2(higher thresholding),3(both)

zero='000F'XU
b1=in32 and zero
b2=ishft(in32,-8) and zero

return,[b1,(b2 and 1) eq 1,(b2 and 2) eq 2]

end;function Get2DbkCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakeStrip2DFilter,w
h=fltarr(w,w)

; w=3,5,7,...

; Intersects
R=w/2
y=indgen(R)+0.5
x=sqrt(R^2-y^2)
azi=atan(y,x)
azi=[azi,!pi/2-azi]
azi=[azi,!pi-azi]
azi=[azi,2*!pi-azi]
azi=azi[sort(azi)]

; Intersect weights
nazi=n_elements(azi)
azidiff=[azi[1:*]-azi[0:nazi-2],azi[0]*2]
weights=R*azidiff/(2*!pi*R)

; Intersect locations
aziavg=[azi[0:nazi-2]+azidiff[0:nazi-2]/2,0]
xpix=round(R*cos(aziavg))+R
ypix=round(R*sin(aziavg))+R

h[xpix,ypix]=weights
return,h

; Old filter:
;pval=1./(4*(w-1))
;h[0,*]=pval
;h[w-1,*]=pval
;h[*,0]=pval
;h[*,w-1]=pval
;return,h
end;function MakeStrip2DFilter
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Strip2D, img, w, niter

niter=0>niter
w=3>w
if (w mod 2) eq 0 then w=w+1

;transform (in this case square root)
bptr=size(img,/type) eq 10
if bptr then imgg=sqrt((*img)>0) else imgg=sqrt(img>0)

;filter with width=w
;imgg=mSavgol(data=imgg,nr=w/2,nl=w/2,r=2,order=[0,0])
imgg=smooth(imgg,w)
imgb=imgg

;stripping
h=MakeStrip2DFilter(w)
for n=1,niter do begin
    imgb=convol(imgb,h,/edge_trunc)
    imgb=imgb<imgg
    imgg=imgb
endfor

;inverse transform (in this case square)
imgb*=imgb

if bptr then return,ptr_new(temporary(imgb)) else return,imgb

return,imgb
end;function Strip2D
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TopHatFilter2D, img, v, w,fvar=fvar,sigma2=sigma2
; w = FWHM
; v in [FWHM/3,FWHM/2]
; r in [2,4]
s=size(img)

stop
a=-v-w/2
b=-a
abin=fix(2*b+1)
mid=abin/2
k=indgen(abin)-b
h=fltarr(abin,abin)
ind1=fix((where(k eq (-w/2-1)))[0])
ind2=fix((where(k eq w/2))[0])
; cross:
h[ind1+1:ind2,ind1+1:ind2]=-1./(2*v)
h[mid,0:ind1]=-1./(2*v)
h[mid,ind1+1:ind2]=1./w
h[mid,ind2+1:*]=-1./(2*v)
h[0:ind1,mid]=-1./(2*v)
h[ind1+1:ind2,mid]=1./w
h[ind2+1:*,mid]=-1./(2*v)
; square:
;h=h-1./(2*v)
;h[ind1+1:ind2,ind1+1:ind2]=1./w
;window,winid()
;shade_surf,h

imgf=convol(img,h,/edge_trunc)
if keyword_set(fvar) then begin
    if not keyword_set(sigma2) then sigma2=img
    fvar=convol(sigma2,h^2.,/edge_trunc)
endif

return, imgf
end; function TopHatFilter2D
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function estims,yb,xb,ybf,mu

CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
;message,!ERR_STRING,/info
return,0
ENDIF

; ----Find boarders of peak----
yb=float(yb)
; start->index of begincut peak
; stopp->index of endcut peak
; determine peakborders:
j=mu
while ybf[j] gt ybf[(j-1)>0] do begin
j=j-1
if j eq 0 then goto,jumpend1
endwhile
jumpend1:start=j

j=mu
while ybf[j] gt ybf[j+1] do begin
j=j+1
if j eq (n_elements(yb)-1) then goto,jumpend2
endwhile
jumpend2:stopp=j

; ----Make different ranges for viewing and fitting---
; x,y   => curvefit-window: [start:stopp]       (width w)
; xw,yw => view-window:     [start-ww:stopp+ww]    (ww=w/2 widen on 2 sides)
; xp,yp => polyfit-window:  [start+www:stopp-www]   (www=w/4 narrow on 2 sides)
w = stopp-start+1
ww = w/2
www = w/4
x=xb[start:stopp];for curvefit
y=yb[start:stopp];for curvefit
board1=start-ww
board2=stopp+ww
xw = xb[0>board1:board2<n_elements(xb)-1];for viewwindow
yw = yb[0>board1:board2<n_elements(xb)-1];for viewwindow
xp=xb[start+www:stopp-www];for polfit
yp=yb[start+www:stopp-www];for polfit

; calculate linear background for x,y:
b1=(yb[stopp]-yb[start])/(xb[stopp]-xb[start])
b0=yb[stopp]-b1*xb[stopp]
; view lin background
;plot,x,y
;oplot,x,b0+b1*x

; ----Model peak with second degree polynomial -> estimation of sigma and area----
weights=1/(yp^2)
xxp=xp-total(xp)/n_elements(xp)
aa=polyfitw(xxp,alog(yp+1),weights,2)
sx=sqrt(-1./(2.*aa[2]))
area=sqrt(-!dpi/aa[2])*exp(aa[0]-aa[1]^2/(4*aa[2]))
; Fit broader window when no solution (curvefit window)
if (finite(sx) eq 0) or (finite(area) eq 0) then begin
    xp=x
    yp=y
    weights=1/(yp^2)
    xxp=xp-total(xp)/n_elements(xp)
    aa=polyfitw(xxp,alog(yp+1),weights,2)
    sx=sqrt(-1./(2.*aa[2]))
    area=sqrt(-!dpi/aa[2])*exp(aa[0]-aa[1]^2/(4*aa[2]))
endif
; Fit broader window when no solution (view window)
if (finite(sx) eq 0) or (finite(area) eq 0) then begin
    xp=xw
    yp=yw
    weights=1/(yp^2)
    xxp=xp-total(xp)/n_elements(xp)
    aa=polyfitw(xxp,alog(yp+1),weights,2)
    sx=sqrt(-1./(2.*aa[2]))
    area=sqrt(-!dpi/aa[2])*exp(aa[0]-aa[1]^2/(4*aa[2]))
endif
; Take arbit. values when still no solution
if finite(sx) eq 0 then sx=1
if finite(area) eq 0 then area=(sqrt(2*!dpi))*10.
;Result = CHECK_MATH( )

; ----Refine peak ->estimate area, position and Sigma again for x,y----
; shift until middle = 0
xx=x-total(x)/n_elements(x)
b0=b0+b1*total(x)/n_elements(x)
A=[area,0,sx,b0,b1]
weights=1/y
if n_elements(y)-n_elements(A) gt 0 then yfit=Marquardt(xx,y, weights, A,$
        SIGMA,FUNCTION_NAME='gfunc2',CHISQ=CHISQ,ITER=ITER,itmax=500)

;difference between refinement and polfit bigger then 95% -> position: x-value of ymax
;                                  -> area and sigma: as determined with polynom. fit
diff=0.95
b0=A[3]
b1=A[4]
b0=b0-b1*total(x)/n_elements(x)
x1=A[1]+total(x)/n_elements(x);mean of gauss
if (abs(A[2]-sx)/sx lt diff) and $
   (abs(A[0]-area)/area lt diff) then begin
       sx=A[2]
       area=A[0]
endif

; ----view peaks + fit: (add breakpoint to see this)----
;window,0
       ; watch peak in widen window
       ;plot,xw,yw,psym=2
       ; fit with polynomial
       ;oplot,xp,exp(aa[0]+aa[1]*xxp+aa[2]*xxp^2),linestyle=2
       ; fit with marquardt
       ;oplot,x,yfit
       ; fitted gaussian (plotted without background)
       ;oplot,x,A[0]/(sqrt(2*!dpi)*A[2])*exp(-(x-x1)^2/(2*A[2]^2)),psym=-4
       ; fitted gaussian (plotted with background)=fit with marquardt
       ;oplot,x,A[0]/(sqrt(2*!dpi)*A[2])*exp(-(x-x1)^2/(2*A[2]^2))+b0+b1*x,psym=-4
       ; fitted gaussian (middle of peak by polynom if not converge)
       ;oplot,x,A[0]/(sqrt(2*!dpi)*A[2])*exp(-(x-list.peakxen[list.apeaks])^2/(2*A[2]^2))
       ; begin of curvefit (not the case when curvefit performed in widen window,cfr supra)
       ;PlotS, [xb[start], xb[start]], [!Y.CRange[0], !Y.CRange[1]],linestyle=2,NOCLIP = 0
       ; end of curvefit
        ;PlotS, [xb[stopp], xb[stopp]], [!Y.CRange[0], !Y.CRange[1]],linestyle=2,NOCLIP = 0

return,sx
end; function estims
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function imgpeaks,img,v,w,drawi,mag,r,imgvalid,sclmax=sclmax,nr=nr,nl=nl

s=size(img)
img=float(img)
if n_elements(imgvalid) eq 0 then imgvalid=1b

; ----Calc sec deriv----
if not keyword_set(nr) then nr=w;/2
if not keyword_set(nl) then nl=w;/2
fvar1=1
fvar2=1
imgg=mSavgol(data=img,nr=w/2,nl=w/2,order=[2,0],fvar=fvar1,sigma2=sigma2)+$
mSavgol(data=img,nr=w/2,nl=w/2,order=[0,2],fvar=fvar2,sigma2=sigma2)
imgg=-imgg
fvar=(fvar1+fvar2)
imgg=fvar>imgg

; ----Filter image with Top-hat filter----
imgf=TopHatFilter2D(imgg,v,w,fvar=fvar)
imgft=imgf
ind=where(imgf le r*sqrt(abs(fvar)),count); r between 2 and 4
if count ne 0 then imgf[ind]=0.

; ----Determine peaks x and y coord.----
s=size(imgf)
dxR = shift(imgf,1,0)
dxL = shift(imgf,-1,0)
dyD = shift(imgf,0,1)
dyU = shift(imgf,0,-1)
m = (imgf gt dxR) and (imgf gt dxL) and (imgf gt dyD) and (imgf gt dyU)
m[0,*]=0
m[*,0]=0
m[s[1]-1,*]=0
m[*,s[2]-1]=0
ind = where(m eq 1 and imgvalid, count)

; ----Return when no peaks----
if count ne 0 then begin
    mx=ind mod s[1]
    my=ind / s[1]
endif else begin
    if keyword_set(sclmax) then begin
       LoadctRev,-3,/silent
       wset,drawi
       tv,rebin(WeakScl(img*imgvalid,sclmax,1),s[1]*mag,s[2]*mag,/sample)
    endif
    return,0
endelse

; ----Refine parameters----
sx=fltarr(count)
sy=sx
A=sx
for i=0L,count-1 do begin
    ; sx and mx
    yb=img[*,my[i]]
    ybf=imgft[*,my[i]]
    xb=indgen(s[1])
    sx[i]=estims(yb,xb,ybf,mx[i])
    ; sy and my
    yb=reform(img[mx[i],*])
    ybf=reform(imgft[mx[i],*])
    xb=indgen(s[2])
    sy[i]=estims(yb,xb,ybf,my[i])
    A[i]=img[mx[i],my[i]]*2.*!dpi*sx[i]*sy[i]
endfor

; ----Plot results of peak detection----
xx=reform(mx)
yy=reform(my)
outimg=img*0.
outimg[xx,yy]=255.
if keyword_set(sclmax) then begin
   LoadctRev,-3,/silent
   wset,drawi
   tv,rebin(outimg,s[1]*mag,s[2]*mag)+ $
     rebin(WeakScl(img*imgvalid,sclmax,1),s[1]*mag,s[2]*mag,/sample)
endif

; ----Return results: A,mx,my,sx,sy,sxy=0,n=0.5----
output=fltarr(7,count)
output[0,*]=A
output[1,*]=mx
output[2,*]=my
output[3,*]=sx
output[4,*]=sy
output[6,*]=0.5
return,output
end;function imgpeaks
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fittd,img,v,w,r,fit,sqr,imgvalid,drawi=drawi,mag=mag,$
    sclmax=sclmax,chisq=chisq,file=file,path=path,printres=printres,$
    debug=debug,flplotid=flplotid

T1=systime(1)
; ----Estimate peak parameters----
; A,mx,my,sx,sy,sxy=0,n
peaks=imgpeaks(img,v,w,drawi,mag,r,imgvalid,sclmax=sclmax)
s=size(peaks)
if s[1] ne 7 then begin
    if keyword_set(printres) then printw,printres,'No peaks found'
    flplotid=0L
    return,0
endif
if s[0] eq 1 then apeaks=1 else apeaks=s[2]

;----Fit img with initial parameters----
peaks[1,*]=peaks[1,*]+sqr[3]
peaks[2,*]=peaks[2,*]+sqr[5]
axbin=sqr[4]-sqr[3]+1
aybin=sqr[6]-sqr[5]+1
n=axbin*aybin
x=indgen(axbin)+sqr[3]
y=indgen(aybin)+sqr[5]
surf=reform(img,n)
weight=1./(surf>1)
xx={x:x,y:y}

T2=systime(1)

;----Fitting----
case fit of
0: begin
    ; ----fit with real constraints----
    strp='Fit with constraints'
    ;ar0,Car,mx0,Cmx,my0,Cmy,sx0,Csx,sy0,Csy,sxy0,Csxy,n0,Cn
    constrar=fltarr(14,apeaks)
    ; large constraints
    ;constrar[0,*]=total(surf)
    ;constrar[1,*]=total(surf)
    ;constrar[2,*]=(sqr[4]-sqr[3])/2.+sqr[3]
    ;constrar[3,*]=(sqr[4]-sqr[3])/2.
    ;constrar[4,*]=(sqr[6]-sqr[5])/2.+sqr[5]
    ;constrar[5,*]=(sqr[6]-sqr[5])/2.
    ;onstrar[6,*]=(sqr[4]-sqr[3])/4.
    ;constrar[7,*]=(sqr[4]-sqr[3])/4.
    ;constrar[8,*]=(sqr[6]-sqr[5])/4.
    ;constrar[9,*]=(sqr[6]-sqr[5])/4.
    ;maxsmult=((sqr[4]-sqr[3])/2.)*((sqr[6]-sqr[5])/2.)
    ;constrar[10,*]=0.
    ;constrar[11,*]=maxsmult
    ;constrar[12,*]=0.5
    ;constrar[13,*]=0.5
    ; small constraints
    constrar[0,*]=total(surf)
    constrar[1,*]=total(surf)
    constrar[2,*]=peaks[1,*]
    constrar[3,*]=5
    constrar[4,*]=peaks[2,*]
    constrar[5,*]=5
    constrar[6,*]=peaks[3,*]
    constrar[7,*]=2
    constrar[8,*]=peaks[4,*]
    constrar[9,*]=2
    maxsmult=(peaks[3,*]+2)*(peaks[4,*]+2)
    constrar[10,*]=0.
    constrar[11,*]=maxsmult
    constrar[12,*]=0.5
    constrar[13,*]=0.5

    printw,printres,'Constraints:'
    printw,printres,constrar

    ; physical meaningful parameters AA
    AA=peaks
    ; calculate fit parameters A
    A=AA;=> physical meaningful parameters
    ;ar0,Car,mx0,Cmx,my0,Cmy,sx0,Csx,sy0,Csy,sxy0,Csxy,n0,Cn
    j=bindgen(7)
    j=[0,1,2,3,4,5,6]
    for i=0l,n_elements(j)-1 do begin
    A[j[i],*]=tan((!dpi/2)*(A[j[i],*]-constrar[2*j[i],*])/constrar[2*j[i]+1,*])
    endfor
    A=reform(A,7*apeaks);=> fit parameters
    struc={constrar:constrar,j:j}

    imgfit = Marquardt(xx,surf, weight, A, SIGMA,FUNCTION_NAME='gfunc2dPVc',$
    CHISQ=CHISQ,ITER=ITER,itmax=100,nfree=nfree,struc=struc)
    A=reform(A,7,apeaks)
    AA=A
    sigma=reform(sigma,7,apeaks)
    for i=0l,n_elements(j)-1 do begin
       sigma[j[i],*]=sigma[j[i],*]*abs(2./!dpi*constrar[2*j[i]+1,*]/(1+A[j[i],*]^2))
       AA[j[i],*]=constrar[2*j[i],*]+2.*constrar[2*j[i]+1,*]/!dpi*atan(A[j[i],*])
    endfor
    endcase
1:begin
    ; ----fit with suggestive constraints (i.e. restraints)----
    strp='Fit with restraints'
    AA=reform(peaks,7*apeaks)

    constr=peaks*0.
    ;Car,Cmx,Cmy,Csx,Csy,Csxy,Cn
    ; large constraints
    constr[0,*]=armin(abs(2*total(surf)-peaks[0,*]),abs(0-peaks[0,*]))
    constr[1,*]=armin(abs(sqr[4]-peaks[1,*]),abs(sqr[3]-peaks[1,*]))
    constr[2,*]=armin(abs(sqr[6]-peaks[2,*]),abs(sqr[5]-peaks[2,*]))
    constr[3,*]=armin(abs((sqr[4]-sqr[3])/2.-peaks[3,*]),abs(0-peaks[3,*]))
    constr[4,*]=armin(abs((sqr[6]-sqr[5])/2.-peaks[4,*]),abs(0-peaks[4,*]))
    constr[5,*]=armin(abs((peaks[3,*]-constr[3,*])*(peaks[4,*]-constr[4,*])-peaks[5,*]),$
    abs((peaks[3,*]+constr[3,*])*(peaks[4,*]+constr[4,*])-peaks[5,*]))
    constr[6,*]=0.5
    ; small constraints
    constr[0,*]=armin(abs(2*total(surf)-peaks[0,*]),abs(0-peaks[0,*]))
    constr[1,*]=5
    constr[2,*]=5
    constr[3,*]=2
    constr[4,*]=2
    maxsmult=(peaks[3,*]+2)*(peaks[4,*]+2)
    constr[5,*]=maxsmult
    constr[6,*]=0.5

    constr=reform(constr,7*apeaks)
    printw,printres,'Restraints Constraints:'
    printw,printres,constr
    imgfit = MarquardtC(xx,surf, weight, AA, SIGMA,FUNCTION_NAME='gfunc2dPV',$
    CHISQ=CHISQ,ITER=ITER,itmax=100,nfree=nfree,constr=constr)
    AA=reform(AA,7,apeaks)
    sigma=reform(sigma,7,apeaks)
    endcase
2:begin
    ; ----fit without constraints----
    strp='Fit without constraints or restraints'
    AA=reform(peaks,7*apeaks)
    imgfit = Marquardt(xx,surf, weight, AA, SIGMA,FUNCTION_NAME='gfunc2dPV',$
    CHISQ=CHISQ,ITER=ITER,itmax=100,nfree=nfree,struc=struc)
    AA=reform(AA,7,apeaks)
    sigma=reform(sigma,7,apeaks)
    endcase
endcase
T3=systime(1)

;----Output result----
if keyword_set(printres) then begin
    printw,printres,'Number of peaks detected: '+string(apeaks)
    printw,printres,'volume, mx, my, sx, sy, sxy, n'
    printw,printres,string(peaks,format='("Initial  ",7F15.4)')
    printw,printres,string(AA,format='("New      ",7F15.4)')
    printw,printres,string(sigma,format='("Stdev    ",7F15.4)')
    printw,printres,'Number of iterations:'+string(iter)
    printw,printres,'Initiation time in sec: '+string(T2-T1)
    printw,printres,'Fitting time in sec/peak:'+string((T3-T2)/apeaks)
    printw,printres,strp
    printw,printres,'Chisq Red:'+string(chisq)
endif

; ----save Fit----
if keyword_set(file) and keyword_set(path) then begin
    error=CutPath(file,file=filepks)
    pathpks = path
    filepks = filepks+'.pks'
    listout={pathpks:pathpks,$
         filepks:filepks,$
         OutID:printres,$
         pathorig:path+file,$
         sqr:sqr,$
         apeaks:apeaks,$
         sigma:sigma,$
         AA:AA}
    error=SavePKS(listout)
endif

; ----compare fit with data----
if keyword_set(debug) then begin
    struct={surf:reform(surf,axbin,aybin),imgfit:reform(imgfit,axbin,aybin),$
       x:x,y:y,xr:[min(x),max(x)],yr:[min(y),max(y)],$
       zr:[min(img),max(img)],proplot:'plotcomp'}
    flplot_obj,struct,topid=flplotid,/xtrue,/ytrue
endif

return,AA
end;function fittd
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function cc,y1,y2,y3,y4,d
return,y2+d*((-y1+y3)+d*((2*y1-2*y2+y3-y4)+d*(-y1+y2-y3+y4)))
end;function cc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bi,y1,y2,d
return,y1*(1-d)+y2*d
end;function bi
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mInterpolate, image,xmap,ymap,axbin,aybin,method=method,printres=printres

if not keyword_set(method) then method=1 else $
method=1>fix(method)<3

; Assign gray level to non-integer coord (xmap,ymap)
case method of
1:  begin
    ;---- Nearest neighbor resampling (zero order interpolation) ----
    if keyword_set(printres) then printw,printres,'Interpolation: Nearest Neighbor'
    oimage=(*image)[round(xmap),round(ymap)]
    endcase
2:  begin
    ;---- Bilinear ----
    if keyword_set(printres) then printw,printres,'Interpolation: Bilinear'
    oimage=Interpolate((size(*image,/type) eq 5)?*image:float(*image),xmap,ymap,missing=0)
;    Kf=floor(xmap)
;    Kc=Kf+1
;    Lf=floor(ymap)
;    Lc=Lf+1
;    c=xmap-Kf
;    d=ymap-Lf
;    t1=bi((*image)[Kf,Lf],(*image)[Kc,Lf],c)
;    t2=bi((*image)[Kf,Lc],(*image)[Kc,Lc],c)
;    oimage=bi(t1,t2,d)
    endcase
3:  begin
    ;---- Cubic ----
    ; approximates the theoretically optimum sinc interpolation function
    if keyword_set(printres) then printw,printres,'Interpolation: Cubic'
    oimage=Interpolate((size(*image,/type) eq 5)?*image:float(*image),xmap,ymap,/cubic,missing=0)

;    Kf=floor(xmap)
;    Kc=Kf+1
;    Kff=Kf-1
;    Kcc=Kc+1
;    Lf=floor(ymap)
;    Lc=Lf+1
;    Lff=Lf-1
;    Lcc=Lc+1
;    c=xmap-Kf
;    d=ymap-Lf
;    t1=cc((*image)[Kff,Lff],(*image)[Kff,Lf],(*image)[Kff,Lc],(*image)[Kff,Lcc],c)
;    t2=cc((*image)[Kf,Lff],(*image)[Kf,Lf],(*image)[Kf,Lc],(*image)[Kf,Lcc],c)
;    t3=cc((*image)[Kc,Lff],(*image)[Kc,Lf],(*image)[Kc,Lc],(*image)[Kc,Lcc],c)
;    t4=cc((*image)[Kcc,Lff],(*image)[Kcc,Lf],(*image)[Kcc,Lc],(*image)[Kcc,Lcc],c)
;    oimage=cc(t1,t2,t3,t4,d)
    endcase
endcase

return,reform(oimage,axbin,aybin)
end;function mInterpolate
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AzIntWarpBatch,pstruc,pimage,pimagevalid,oimage=oimage

; Initialize
nsectors=(*pstruc).nsectors
axbin=(*pstruc).axbin
aybin=(*pstruc).aybin
percentile=(*pstruc).percentile
method=(*pstruc).method
median=(*pstruc).median
calcerror=(*pstruc).calcerror

out=fltarr(1+nsectors+nsectors*calcerror,axbin)
out[0,*]=(*pstruc).xrange

; Jcor
bnormchanged=0b
if n_elements(pimagevalid) ne 0 and (*pstruc).updatenorm then $
if ptr_valid(pimagevalid) then begin
    norm = (*pstruc).normfixed and (mInterpolate(pimagevalid,(*pstruc).xmap,(*pstruc).ymap,axbin,aybin,method=method) eq 1)
    
    ; 2D intensity correction
    (*pstruc).Jcorn=(*pstruc).Jcor*norm
    
    if median then begin
        ; Pixels that are zeroed out, will not affect total()
        ; but they do affect median(). Therefore make them NaN.
        ind=where(norm eq 0,ct)
        if ct ne 0 then ((*pstruc).Jcorn)[ind]=!values.F_NAN
    endif

    bnormchanged=1b
endif

; Interpolate diffraction pattern
oimage=mInterpolate(pimage,(*pstruc).xmap,(*pstruc).ymap,axbin,aybin,method=method)

; Applied 2D correction
oimage*=(*pstruc).Jcorn

; Azimuthal mean/median + 1D correction
if nsectors gt 1 then begin
    ; n + 1D correction
    if bnormchanged then (*pstruc).Jac1Dn=(*pstruc).Jac1D*(total(reform(norm,axbin,aybin/nsectors,nsectors),2)>1)
    
    ; Error for azimuthal average
    if calcerror then out[nsectors+1:2*nsectors,*]=transpose(sqrt(total(reform(oimage*(*pstruc).Jcorn,axbin,aybin/nsectors,nsectors),2,/nan))/(*pstruc).Jac1Dn)
    
    ; Azimuthal integration
    if keyword_set(median) then begin
        out[1:nsectors,*]=transpose(percentile(reform(oimage,axbin,aybin/nsectors,nsectors),pct=percentile,dim=2,/even)/(*pstruc).Jac1D)
        if calcerror then out[nsectors+1:2*nsectors,*]*=sqrt(!pi/2)
    endif else begin
        out[1:nsectors,*]=transpose(total(reform(oimage,axbin,aybin/nsectors,nsectors),2)/(*pstruc).Jac1Dn)
    endelse

endif else begin
    ; n + 1D correction
    if bnormchanged then (*pstruc).Jac1Dn=(*pstruc).Jac1D*(total(temporary(norm),2)>1)
    
    ; Error for azimuthal average
    if calcerror then out[2,*]=sqrt(total(oimage*(*pstruc).Jcorn,2,/nan))/(*pstruc).Jac1Dn
    
    ; Azimuthal integration
    if keyword_set(median) then begin
        out[1,*]=percentile(oimage,pct=percentile,dim=2,/even)/(*pstruc).Jac1D
        if calcerror then out[2,*]*=sqrt(!pi/2)
    endif else begin
        out[1,*]=total(oimage,2)/(*pstruc).Jac1Dn
    endelse
endelse

return,out
end;function AzIntWarpBatch
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AzIntWarp,mask,list,IntCor2Dinfo,type,rmin_,dr_,outbins,AzWidth,azbins,flags,printres=printres,$
    saverp=saverp,listsaverp=saverplist,imagevalid=imagevalid,batchreturnp=batchreturnp,$
    splitnsectors=nsectors,updatenorm=updatenorm,calcerror=calcerror,$
    oimage=oimage,xrange=xrange,yrange=phi,median=median,percentile=percentile

; type: output x-axis format (d(Angstrom),tt(deg),R(ort,mm),Q(1/nm))
; dr: output x-axis binwidth (A,rad,(ort,pix),1/nm)
; rmin: this is to handle ellipses with two-theta=0 or 180 (A,rad,(ort,pix),1/nm)
; outbins: if this is 0, the number of bins is choosen according to dr
; AzWidth: azimuthal binwidth (radians)
; azbins: if this is 0, the number of bins is choosen according to AzWidth
; flags[0]:fill used pixels on image
; flags[1]:
;     0. NN (+/- equivalent to Bresenham's algorithm)
;     1. Bilinear interpol.
;     2. Cubic interpol. (equivalent to ideal inc interpol.)
; printres: ID to print process info

tstart=systime(1)
saverp=keyword_set(saverp)
updatenorm=keyword_set(updatenorm)
calcerror=keyword_set(calcerror)
median=keyword_set(median)
s=*list.tiffs
if ~keyword_set(nsectors) then nsectors=1
if ~keyword_set(percentile) then percentile=50
nsectors=round(nsectors>1)

;Convert (ort,pix) to (ort,mm): depends on azimuth, so take the smallest
rmin=rmin_
dr=dr_
rmin[2]*=1000*(list.scrx<list.scry)
dr[2]*=1000*(list.scrx<list.scry)
; Convert radians to degrees
rmin[1]*=180/!pi
dr[1]*=180/!pi

; ----1) Prepare d-spacing and azimuthal limits----
; x: integration limits (d-spacing)
; azibegin,aziend: azimuthal integration limits
case mask[0] of
2: begin
    x=mask[[3,5]]
    azi=mask[[4,6]]
    azibegin=azi[0]<azi[1]
    aziend=azi[0]>azi[1]
    endcase
3: begin
    x=mask[[3,4]]
    azibegin=0
    aziend=2*!dpi
    endcase
endcase
if rmin[0] ne 0 then $
if list.dist lt 0 then x>=rmin[0]$
else x<=rmin[0] ; Since x is in Angstroms, use rmin[0] (which is a maximum for d-spacing)

; ----2) Profile x-axis----
; xrange, axbin, incx => 1D-scan properties (equally spaced)

; convert d-spacing limits to type
x_d=x
case type of
    0:     begin
        ; ----fixed d-spacing increment (Angstroms)----
           ; x = D-spacing (Angstroms)
           endcase
    1:     begin
           ; ----fixed 2-theta increment (degrees)----
           ; x = 2-Theta (degrees)
        x=BraggDtoX(x,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/onlyt,/angledeg)
           endcase
    2:     begin
           ; ----fixed radial increment (mm)----
           ; x = orthogonal radial distance (mm)
        x=BraggDtoX(x,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/onlyp,/pixm)
           endcase
    3:    begin
        ; ----scattering vector (1/nm)----
           ; x = Q-space (1/nm)
           x=BraggDtoX(x,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,/onlyq)
           endcase
endcase

; equally spaced "type" range
xbegin=x[0]<x[1]
xend=x[0]>x[1]
if (outbins ne 0) then axbin=long(outbins) else $
if (dr[type] eq 0) then begin
    xclose=abs(BraggDtoX(x_d,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,/onlyp,/nonortp))
    axbin=ceil(abs(xclose[1]-xclose[0]))+1
endif else axbin=ceil((xend-xbegin)/dr[type])+1

incx=(xend-xbegin)/(axbin-1)
xrange=incx*lindgen(axbin)+xbegin
if keyword_set(printres) then begin
    printw,printres,'Output range: '+string(xbegin)+string(xend)
    printw,printres,'Output spacing: '+string(incx)
    printw,printres,'Number of output bins:'+string(axbin)
endif

; ----3) Azimuthal range----
phipixel = [azibegin,aziend]
phimm = azPixeltoMm(phipixel,list.scrx,list.scry)
azidiff=phimm[1]-phimm[0]

if (azbins gt 0) then aybin=azbins else $
if AzWidth eq 0 then begin
    ; Estimate arc length of the maximal conic section in pixels
    Rarc =max(abs(BraggDtoX(x_d[0],list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,/onlyp,/nonortp,phi=findgen(360)*!pi/180)))
    Rarc>=max(abs(BraggDtoX(x_d[1],list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,/onlyp,/nonortp,phi=findgen(360)*!pi/180)))
    Larc=Rarc*(phipixel[1]-phipixel[0]) ;TODO: arc length for ellipse, not circle
    
    ; Choose bins so that the maximal conic section has 1-pixel bins
;    aybin=ceil(Larc/sqrt(2))+1
    aybin=ceil(Larc*0.8)
endif else aybin=ceil(azidiff/AzWidth)+1 ;round instead of ceil to avoid 360?

; Even times the number of sectors
if nsectors gt 1 then aybin=ceil(aybin/float(nsectors))*nsectors

incazi=azidiff/(aybin-1)
phi=incazi*indgen(1,aybin)+azibegin ; in mm-frame
if keyword_set(printres) then begin
    printw,printres,'Azimuthal binwidth ('+msymbols('degrees')+'):'+string(incazi*180/!dpi)
    printw,printres,'Number of azimuthal bins:'+string(aybin)
endif

; ----4) convert from (xrange,phi) to (xpix,ypix)----
; We have at this point
;    xrange: "type" in (Angstrom, degrees, mm_ort or 1/nm)
;    phi: azimuth in radians (in mm frame)
;
; Optimized code
BraggXtoXY_I_optimized,xrange,phi,xmap,ymap,Jcor,Jac1D,$
                        list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,$
                        type,IntCor2Dinfo,validphi=norm;,debug=debug

; Non optimized code for debugging
;xmap_=rebin(xrange,axbin,aybin,/sample)
;ymap_=rebin(phi,axbin,aybin,/sample)
;ret=BraggXtoXY_I(xmap_,ymap_,1.,list.lambda,list.dist,list.scrx,list.scry,list.a,list.b,list.phipixphi0,list.center,type,IntCor2Dinfo,/angledeg,/pixm,nonortp=0,validphi=norm)
;xmap_=reform(ret.x11,axbin,aybin)
;ymap_=reform(ret.x12,axbin,aybin)
;Jcor_=rebin([Jac1D],axbin,aybin,/sample)/reform(ret.y1,axbin,aybin)
;
;if ~floatequal(xmap,xmap_)then stop,max(xmap-xmap_,/abs)
;if ~floatequal(ymap,ymap_)then stop,max(ymap-ymap_,/abs)
;if ~floatequal(Jcor,Jcor_)then stop,max(Jcor-Jcor_,/abs)

; ----5) find (non-integer) values that fall out off 
;        input image or that fall in a masked off area
;        or are invalid----
; Norm indicates the valid pixels
norm and= (xmap ge 0) and (xmap le (s[1]-1))
norm and= (ymap ge 0) and (ymap le (s[2]-1))
fJcor=finite(Jcor)
norm and= fJcor
xmap=0>xmap<(s[1]-1)
ymap=0>ymap<(s[2]-1)
ind=where(~temporary(fJcor),ct)
if ct ne 0 then Jcor[ind]=0

n=axbin*aybin
xmap=reform(xmap,n,/overwrite)
ymap=reform(ymap,n,/overwrite)

if n_elements(batchreturnp) ne 0 and updatenorm then normfixed=norm else normfixed=0 ; non-variable part of norm

if n_elements(imagevalid) ne 0 then $
if ptr_valid(imagevalid) then $
    norm and= mInterpolate(imagevalid,xmap,ymap,axbin,aybin,method=flags[1]+1) eq 1 ; variable part of norm

; ----6) assign gray level to non-integer coord (xmap,ymap)----
; ie: non-integer coord. in 'image' but integer coord. in 'oimage'
oimage=mInterpolate(list.tiff,xmap,ymap,axbin,aybin,method=flags[1]+1,printres=printres)
;oimage=reform(Interpolate(*list.tiff,xmap,ymap,/cubic),axbin,aybin)

; ----7) display XY----
if flags[0] eq 1 then begin
    ; dynamic view of azimuthal integration
    wset,list.drawdynindex
    ind=where(norm ne 0,ct)
    if ct ne 0 then begin
        xpix=RTDynExt(xmap[ind],list.sf,list.hrange[0],/x)
        ypix=RTDynExt(ymap[ind],list.sf,list.vrange[0])
        Plots,xpix,ypix,/device,psym=3
        ;Device, Copy = [0,0, list.dynh,list.dynv, 0,0,list.pixindex2]
    endif
endif

; -----8) prepare output----
out=fltarr(1+nsectors*(1+calcerror+saverp),axbin)
out[0,*]=xrange

; -----9) ring percentage----
if saverp then begin
    ; background level
    if saverplist.done then begin
        ; count pixels above background level
        tmp=oimage*norm gt saverplist.thres[0]
    endif else begin
        oimageb=Strip2D(oimage,saverplist.wb,saverplist.niter)

        ; count pixels above background level
        if saverplist.thres[0] ge 0 then begin
            tmp=(oimage-oimageb)*norm gt saverplist.thres[0]
        endif else begin
               oimageb+=sqrt(oimageb); + k*Standard Deviation
               tmp=(oimage gt oimageb)*norm
           endelse
    endelse
    
    ; sectors?
    if nsectors gt 1 then begin
        tmp=reform(tmp,axbin,aybin/nsectors,nsectors)
        out[(calcerror+1)*nsectors+1:(calcerror+2)*nsectors,*]=total(tmp,2)/(total(reform(norm,axbin,aybin/nsectors,nsectors),2)>1)*100.
    endif else begin
        ; count for each column and normalize
        out[(calcerror+1)*nsectors+1,*]=total(tmp,2)/(total(norm,2)>1)*100.
    endelse

endif

; ----Debug output---
;AzintExampleOutput,create_struct(debug,'norm',norm,'oimage',oimage),list
;return,0

; ----10) convert I(x,y) to I(x)----
; 2D intensity correction
Jcorn=Jcor*norm
if ~updatenorm then Jcor=0

if median then begin
    ; Pixels that are zeroed out, will not affect total()
    ; but they do affect median(). Therefore make them NaN.
    ind=where(norm eq 0,ct)
    if ct ne 0 then Jcorn[ind]=!values.F_NAN
endif

; TODO: check whether Jcor and Jac1D are always possitive
oimage*=Jcorn ; I(tt,phi) in (rad,rad)

; Azimuthal averaging:
; yi = ai/ni . SUM_j (bj.Iij)
; sigma(yi) = ai/ni . sqrt(SUM_j (bj^2.Iij))
; 
;   bj -> Jcor
;   ai -> 1/Jac1D
;   ni -> norm
;
; Azimuthal median:
; yi = ai . bj.Iij   where j indicates the median
; sigma(yi) = sqrt(!pi/2)*sigma(yi previous)

; Azimuthal mean/median + 1D correction
if nsectors gt 1 then begin
    ; n + 1D correction
    Jac1D=rebin([Jac1D],axbin,nsectors,/sample)
    Jac1Dn=Jac1D*(total(reform(norm,axbin,aybin/nsectors,nsectors),2)>1)
    
    ; Error for azimuthal average
    if calcerror then out[nsectors+1:2*nsectors,*]=transpose(sqrt(total(reform(oimage*Jcorn,axbin,aybin/nsectors,nsectors),2,/nan))/Jac1Dn)
    
    ; Azimuthal integration
    if median then begin
        out[1:nsectors,*]=transpose(percentile(reform(oimage,axbin,aybin/nsectors,nsectors),pct=percentile,dim=2,/even)/Jac1D)
        if calcerror then out[nsectors+1:2*nsectors,*]*=sqrt(!pi/2)
    endif else begin
        out[1:nsectors,*]=transpose(total(reform(oimage,axbin,aybin/nsectors,nsectors),2)/Jac1Dn)
    endelse

endif else begin
    ; n + 1D correction
    Jac1Dn=Jac1D*(total(temporary(norm),2)>1)
    
    ; Error for azimuthal average
    if calcerror then out[2,*]=sqrt(total(oimage*Jcorn,2,/nan))/Jac1Dn
    
    ; Azimuthal integration
    if median then begin
        out[1,*]=percentile(oimage,pct=percentile,dim=2,/even)/Jac1D
        if calcerror then out[2,*]*=sqrt(!pi/2)
    endif else begin
        out[1,*]=total(oimage,2)/Jac1Dn
    endelse
endelse

if n_elements(batchreturnp) ne 0 then begin
    batchreturnp=ptr_new({axbin:axbin,aybin:aybin,method:flags[1]+1,xmap:temporary(xmap),ymap:temporary(ymap),$
                normfixed:temporary(normfixed),Jcor:temporary(Jcor),Jcorn:temporary(Jcorn),Jac1D:temporary(Jac1D),Jac1Dn:temporary(Jac1Dn),$
                xrange:xrange,median:median,percentile:percentile,nsectors:nsectors,$
                updatenorm:updatenorm,calcerror:calcerror,$
                ROIAzimuth:[azibegin,aziend]})
endif

if keyword_set(printres) then printw,printres,'AzIntWarp processing time:'+string(systime(1)-tstart)
return,out
end;function AzIntWarp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inverseAzIntWarp,type,xval,spe,imgsize,lambda,dist,scrx,scry,dega,degb,phipixphi0,center,IntCor2Dinfo

case type of
0: onlyd=1
1: onlyt=1
2: onlyp=1
3: onlyq=1
endcase

res=BraggXYtoX_I(rebin(lindgen(imgsize[0]),imgsize,/sample),rebin(lindgen(1,imgsize[1]),imgsize,/sample),1.,$
            lambda,dist,scrx,scry,dega,degb,phipixphi0,center,IntCor2Dinfo,onlyp=onlyp,onlyd=onlyd,onlyt=onlyt,$
            onlyq=onlyq,/pixm,/angledeg);,aziout=phi

if n_elements(spe) eq 1 then img=spe[0]/res.(1) $
else img=reform(minterpol( spe, xval, res.(0), /SPLINE),imgsize)/res.(1)

ind=where(~finite(img),ct)
if ct ne 0 then img[ind]=0

if n_elements(img) eq 1 then img=make_array(imgsize,value=img)

return,img
end;function inverseAzIntWarp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function EqualBin, x1, I1, x20, x2i, x2n, X2=x2

; Rebin array x1 so the total integrand stays the same
;
; x1 equally or unequally spaced with I1
; x2 equally spaced (x20: start, x2i: incr, x2n: #elements) with I2
;
; I1 and I2 might be called "the mean density" for each interval

;eg: print,EqualBin([0,1,3,6],[1,2,1,3.],0.5,1,4)
;   must give [1.5,2,1.5,1]

x2=x20+x2i*lindgen(1,x2n+1)
I2=fltarr(x2n)

n=n_elements(x1)
x1=reform(x1)
I1=reform(I1)

A=rebin(x1,n,x2n+1,/sample)
X=rebin(x2,n,x2n+1,/sample)
B=A-X ; n col, x2n+1 row

for i=0l,x2n-1 do begin ; loop over rows in pairs
    ind1=(where(B[*,i] ge 0,ct1))[0]
    ind2=(where(B[*,i+1] ge 0,ct2))[0]
    if ind1 eq ind2 then begin
        if ((ind1 ne 0)and (ct1 ne 0)) then I2[i]=I2[i]+I1[ind1-1]*x2i ; all from 1 bin
    endif else begin
        if ((ind1 ne 0) and (ct1 ne 0)) then I2[i]=I2[i]+I1[ind1-1]*B[ind1,i]
        if (ct2 ne 0) then I2[i]=I2[i]+I1[ind2-1]*(-B[ind1,i+1])$
        else I2[i]=I2[i]+I1[n-1]*(-B[ind1,i+1])
    endelse
endfor

x2=reform(x2[0,0:x2n-1])
return,I2/x2i
end;function EqualBin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RemoveSatur,image,imgvalid,threshold,aspect_threshold,npix_threshold,count=count

dims = Size(*image, /Dimensions)
count=0

; Label all blobs with value lower than threshold
lab = Label_Region((*image) gt threshold)
h = histogram(lab,omin=omin)
if n_elements(h) eq 1 then return
h = h[1: *] ; removes background pixel count

; Limit area
w = Where(h gt npix_threshold, n_areas)+1

; Limit aspect ratio
aspect_threshold = [aspect_threshold,1./aspect_threshold]

; Find saturated blobs with streaks
count=0l
For i = 0, n_areas - 1 Do Begin
    ; Blob pixels
    wa = Where(lab EQ w[i],ct)
    wx = wa Mod dims[0]
    wy = wa / dims[0]
    
    ; Check aspect ratio
    max=max(wx,min=mix)
    may=max(wy,min=miy)
    nx=max - mix + 1.0
    ny=may - miy + 1.0
    ratio=ny/nx
    if ratio gt aspect_threshold[0] or ratio lt aspect_threshold[1] then begin
        if not ptr_valid(imgvalid) then imgvalid=ptr_new(make_array(size(*image,/dim),value=1b))
        (*imgvalid)[wa]=0b
        count+=ct
    endif
EndFor

end;pro RemoveSatur
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RemoveZingers,image,imgvalid,width=width,threshold=threshold,count=count

; Default kernel width for smoothing
if (keyword_set(width) eq 0) then width=3
width=fix(width)

; Default threshold for detecting zingers
if (keyword_set(threshold) eq 0) then threshold=1.2

; Compute ratio of unsmoothed image to smoothed image
;smoothed=mSavGol(data=image,order=[0,0],nr=width/2,nl=width/2,r=3)
smoothed=smooth(*image,width)
ratio = (*image)/smoothed

; Find indices of pixels which are above threshold
zinger = where(ratio gt threshold, count)

; Replace pixels with zingers by the smoothed
if (count gt 0) then begin
    if not ptr_valid(imgvalid) then imgvalid=ptr_new(make_array(size(*image,/dim),value=1b))
    (*imgvalid)[zinger] = 0
;    (*image)[zinger] = smoothed[zinger]
endif

end;function RemoveZingers
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CalcTieO,nRings,nrTie,TieI,center,lambda,dist,scrx,scry,dspac,dega,degb
TieO=PTR_NEW((*TieI)*0.)
ind=0

for i=0l,nRings-1 do begin
    xi=(*TieI)[0,0+ind:(*nrTie)[i]+ind-1]-center[0]
    yi=(*TieI)[1,0+ind:(*nrTie)[i]+ind-1]-center[1]
    phi=atan(yi,xi)

    if (*dspac)[i] ne 0 then begin
        ri=BraggDtoX((*dspac)[i],lambda,dist,scrx,scry,/onlyp,dega,degb,phi=azPixeltoMm(phi,scrx,scry))
    endif else begin
        r=sqrt(xi^2.+yi^2.)
        ;rgem=total(r)/(*nrTie)[i]
        ri=min(r)
    endelse
    t=atan(yi,xi)
    xo=ri*cos(phi)+center[0]
    yo=ri*sin(phi)+center[1]

    (*TieO)[0,0+ind:(*nrTie)[i]+ind-1]=xo
    (*TieO)[1,0+ind:(*nrTie)[i]+ind-1]=yo
    ind=ind+(*nrTie)[i]
endfor

return,TieO
end;function CalcTieO
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bsplineBF, n, p, m, u, uvec

; uvec: clamped knot vector in domain [0,1]
; m+1: number of knots

; u: points (within [uvec[0],uvec[m]]) for which we want to know the BF

; p: degree
; n+1: number of control points and data points (n+1=m-p)

Nbf=fltarr(m+1); N0,p(u) ... Nm,p(u)

; Calculate the first order basis functions: N0,0(u) ... Nm,0(u)
if u lt uvec[0] then u=uvec[p] $
else if u gt uvec[m] then u=uvec[n+1]
if u eq uvec[0] then begin
    Nbf[0] = 1.0
    return,Nbf[0:n]
endif else if u eq uvec[m] then begin
    Nbf[n] = 1.0
    return,Nbf[0:n]
endif
k=(where((uvec-u) gt 0))[0]-1 ; u in knot span [uvec[k],uvec[k+1][
Nbf[k] = 1.0

; Calculate the higher order basis functions: N0,d(u) ... Nm,d(u)
for d=1,p do begin
    ; The only terms that might be non-zero for degree d are Nk-d,...,Nk
    for i = ((k-d)>0), k do begin
        if (Nbf[i] ne 0) then c1=(u-uvec[i])/(uvec[i+d]-uvec[i])*Nbf[i] $
        else c1=0
        if (Nbf[i+1] ne 0) then c2=(uvec[i+d+1]-u)/(uvec[i+d+1]-uvec[i+1])* Nbf[i+1] $
        else c2=0
        Nbf[i]=c1+c2
    endfor
endfor

return,Nbf[0:n]
end;function bsplineBF
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function bsplineint2Dcp,uvec,vvec,p,q,n,m,h,k,seval,teval,CPT

; Nbfu: (nseval)x(n+1)
; Nbfv: (nteval)x(m+1)
; CP: (n+1)x(m+1)
; => D=Nbfu##CP##NbfvT   (nseval)x(nteval)
; => DT=Nbfv##CPT##NbfuT (nteval)x(nseval)

; Evaluate b-spline for seval,teval
nDs=n_elements(seval)
NbfuT=fltarr(nDs,n+1) ; (n+1)x(...)
for i=0l,nDs-1 do NbfuT[i,*]=bsplineBF(n, p, h, seval[i], uvec)
nDt=n_elements(teval)
Nbfv=fltarr(m+1,nDt) ; (...)x(m+1)
for i=0l,nDt-1 do Nbfv[*,i]=bsplineBF(m, q, k, teval[i], vvec)

return,Nbfv##CPT##NbfuT
end;function bsplineint2Dcp
;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro resample_array, F, psrc, pdst, len, ipo, off

; len: length of column or row
; ipo: pixel increment in pscr or pdst
; off: offset when accessing psrc or pdst

; ----precompute input index for each output pixel----
inpos=fltarr(len+2) ; forward mapping function

; linear interpolation:
; p(x)=f(x0)+(x-x0)(f(x1)-f(x0))/(x1-x0)
;    x: F[ui] -> output pixel index for ui
;    f(x): ui -> input pixel index
; p(xi)=ui+(xi-F[ui])(ui+1-ui)/(F[ui+1]-F[ui])=ui+(xi-F[ui])/(F[ui+1]-F[ui])
;    xi: output pixel indices
;    inpos: input pixel index for xi
ui=0
for xi=0,len-1 do begin
    while (ui lt (len-1)) and (F[(ui+1)<(len-1)] lt xi) do ui++
    if (ui lt (len-1)) then inpos[xi] = ui + (xi - F[ui]) / (F[ui+1] - F[ui]) $
    else inpos[xi] = ui + 1.0
endfor
inpos[len] = len

ui=0 ; input pixel index
xi=0 ; output pixel index
acc=0. ; accumulator
intensity=0.
INSFAC=inpos[1] ; inverse scale factor
INSEG=1.0 ; entire input pixel is available
OUTSEG=inpos[1] ; # input pixels that map onto 1 output pixel

uimax=(len-1)*ipo
while (xi lt len) do begin
    ; use linear interpolation for reconstruction
    if ui le uimax then begin
        intensity = INSEG * (*psrc)[off+ui]
        if ui lt uimax then intensity+= (1-INSEG) * (*psrc)[off+ui+ipo]
    endif else intensity = 0

    ; INSEG < OUTSEG: input pixel is entirely consumed before output pixel
    if(INSEG lt OUTSEG) then begin
      acc += (intensity * INSEG);  /* accumulation of weighted contrib */
      OUTSEG -= INSEG;             /* INSEG portion has been filled */
      INSEG = 1.0;                 /* new input pixel will be available */
      ui += ipo;                   /* index to next input pixel */

    ; INSEG >= OUTSEG: input pixel is not consumed before output pixel */
    endif else begin
      acc += (intensity * OUTSEG); /* accumulate weighted contrib */

      if (INSFAC eq 0) then INSFAC = 1 ;It seems that when INSFAC is zero, so is acc

      (*pdst)[off+xi*ipo] = acc/INSFAC;    /* init output w/ normalized accumulator */
      acc = 0.0;                   /* reset accumulator for next output pixel */
      INSEG -= OUTSEG;             /* OUTSEG portion of input has been used */
      xi++;                        /* index to next output pixel */
      INSFAC = inpos[xi+1] - inpos[xi]; /* init spatially varying INSFAC */
      OUTSEG = INSFAC;             /* init spatially varying SIZFAC */
    endelse
endwhile

end;pro resample_array
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro WarpSpaDistRunSpline,img,struc,tiffvalid=imgvalid,fast=fast
; Forward warping: (*struc).x, (*struc).y are coordinates in the output image
; where the center of the corresponding input pixels are mapped.

imgs=size(*img)

; Remark: Some strange things happen when img is UINT and not e.g. float
if keyword_set(fast) then begin
    tmp=MAKE_ARRAY(imgs[1:2],type=size(*img,/type))
    
    roundx=round((*struc).x)
    roundy=round((*struc).y)
    indkeep=where(roundx ge 0 and roundx lt imgs[1] and roundy ge 0 and roundy lt imgs[2])
    inddest=roundx[indkeep] + imgs[1] * roundy[indkeep]
    OpWithDuplicate,inddest,tmp,(*img)[indkeep],'mean' ; missing pixels!!!
    *img = temporary(tmp)
    
    ; Bad pixels
    if ptr_valid(imgvalid) then begin
        tmp=bytarr(imgs[1:2])
        OpWithDuplicate,inddest,tmp,(*imgvalid)[indkeep],'or'
        *imgvalid = temporary(tmp)
    endif else begin
        imgvalid=ptr_new(bytarr(imgs[1:2]))
        (*imgvalid)[inddest]=1b
    endelse
    
endif else begin

; ---------Old algo (takes too long)----------
;    interimg=ptr_new(*img)
;
;    ; loop over rows
;    for i=0l,imgs[2]-1 do $
;        resample_array, (*struc).x[*,i], img, interimg, imgs[1], 1, i*imgs[1]
;    ; loop over columns
;    for i=0l,imgs[1]-1 do $
;        resample_array, (*struc).y[i,*], interimg, img, imgs[2], imgs[1], i
;
;    ptr_free,interimg
;    return
; ---------Old algo----------

    ; Intermediate image: add boarder of two pixels
    imgsinter=imgs+4
    type=size(*img,/type)
    interimg=MAKE_ARRAY(imgsinter[1:2],type=type)
    if ptr_valid(imgvalid) then interimgvalid=MAKE_ARRAY(imgsinter[1:2],value=1b)
    bconvert=type ne 4 or type ne 5
    if bconvert then fracsum=float(interimg) $
    else fracsum=interimg

    ; Divide pixels over 9 neighbours (at least 5 will get 0% of pixel)
    ; first pixel of the 4 receiving output pixels (the ouput pixel that contains
    ; the left bottom corner of the mapped input pixel)
    npix=imgs[1]*imgs[2]
    v=[reform((*struc).x,1,npix),reform((*struc).y,1,npix)]
    inddest=round(v-0.5)
    ; Overlap between mapped input pixel and first receiving output pixel
    dxy=inddest-v+1
    ; area of 4 rectangles (total area=1)
    frac=[dxy[0,*]*dxy[1,*], (1.-dxy[0,*])*dxy[1,*], $
        dxy[0,*]*(1.-dxy[1,*]), (1.-dxy[0,*])*(1.-dxy[1,*])]

    ; For each input pixel: 4 output pixels + values
    ; Pixels that fall off: clip (double boarder will catch them)
    indx=rebin([2,3,2,3],4,npix,/sample) ; [0,1,0,1]+2 for boarder
    indx+=rebin(inddest[0,*],4,npix)
    indx>=0
    indx<=imgsinter[1]-1
    
    indy=rebin([2,2,3,3],4,npix,/sample) ; [0,0,1,1]+2 for boarder
    indy+=rebin(inddest[1,*],4,npix)
    indy>=0
    indy<=imgsinter[2]-1
    
    inddest=temporary(indx) + imgsinter[1] * temporary(indy)
    add_pix=rebin(reform((*img),1,npix),4,npix,/sample)*frac
    
    ; Divide input pixels over output pixels
    if ptr_valid(imgvalid) then begin
        tmp=rebin(reform((*imgvalid),1,npix),4,npix,/sample)
        OpWithDuplicate,inddest,interimg,temporary(add_pix),'+',fracsum,temporary(frac),'+',interimgvalid,temporary(tmp),'and'
    endif else OpWithDuplicate,inddest,interimg,temporary(add_pix),'+',fracsum,temporary(frac),'+'

    ; Cut boarders
    interimg=interimg[2:1+imgs[1],2:1+imgs[2]]
    if ptr_valid(imgvalid) then interimgvalid=interimgvalid[2:1+imgs[1],2:1+imgs[2]]
    fracsum=fracsum[2:1+imgs[1],2:1+imgs[2]]
    
    ;    ; interpolation of the pixels that didn't receive anything
;    if ncomp ne 0 then begin
;        (*struc).x=comp/imgs[1]
;        (*struc).y=comp mod imgs[1]
;
;        ; This is slow and something is wrong
;;        TRIANGULATE, (*struc).x, (*struc).y, tr, bounds
;;        gs = [1,1]    ;Grid spacing: for pixels -> 1 in each direction
;;        b = [0,0, imgs[1]-1, imgs[2]-1] ;Bounds: 0,0,maxx,maxy
;;        interimg=TRIGRID((*struc).x, (*struc).y,interimg[(*struc).x, (*struc).y],tr, gs, b, /QUINT,EXTRA = bounds)
;
;        ; This is faster and because averaging is used, less accurate
;        h=[[0.125,0.125,0.125],[0.125,0,0.125],[0.125,0.125,0.125]]
;        interimg[(*struc).x,(*struc).y]=(convol(float(interimg),h,/edge_trunc))[(*struc).x,(*struc).y]
;    endif
;    
    ; Bad pixels
    validpixels=fracsum ne 0
    if ptr_valid(imgvalid) then begin
        interimgvalid and= validpixels
        *imgvalid=temporary(interimgvalid)
    endif else imgvalid=ptr_new(validpixels)

    ; Normalize by the total covered fraction (missing pixels=0)
    interimg/=fracsum
    ind=where(validpixels,comp=comp,ncomp=ncomp)
    if ncomp ne 0 then interimg[comp]=0
    ;if bconvert then interimg=TypeConvert(temporary(interimg),type)
    
    (*img)=temporary(interimg)
endelse

end;pro WarpSpaDistRunSpline
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WarpSpaDistPrepSpline,s,tck1,tck2,vflip=vflip

seval=indgen((*s)[1])
teval=indgen((*s)[2])

; Forward X and Y distortions
x=bsplineint2Dcp(tck1.uvec,tck1.vvec,tck1.p,tck1.q,tck1.n,tck1.m,tck1.h,tck1.k,seval,teval,tck1.CP)
y=bsplineint2Dcp(tck2.uvec,tck2.vvec,tck2.p,tck2.q,tck2.n,tck2.m,tck2.h,tck2.k,seval,teval,tck2.CP)

if keyword_set(vflip) then begin
    x=reverse(x,2)
    y=-reverse(y,2)
endif

; XY-destination
x+=rebin(seval,(*s)[1],(*s)[2],/sample)
y+=rebin(reform(teval,1,(*s)[2]),(*s)[1],(*s)[2],/sample)

p=PTR_NEW({x:x,y:y})
return,p
end;function WarpSpaDistPrepSpline
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function WarpSpaDistPrep,s,TieI,TieO,fast=fast
gs = [1,1]                ;Grid spacing: for pixels -> 1 in each direction

if keyword_set(fast) then begin
    TRIANGULATE, (*TieO)[0,*], (*TieO)[1,*], tr, bounds
    b = [0,0, (*s)[1]-1, (*s)[2]-1] ;Bounds: 0,0,maxx,maxy

    x=TRIGRID((*TieO)[0,*],(*TieO)[1,*],(*TieI)[0,*],tr, gs, b, /QUINT, missing=!values.F_NAN)
    y=TRIGRID((*TieO)[0,*],(*TieO)[1,*],(*TieI)[1,*],tr, gs, b, /QUINT, missing=!values.F_NAN)
endif else begin
    x=GRID_TPS((*TieO)[0,*],(*TieO)[1,*],(*TieI)[0,*], START=[0,0],DELTA=gs,NGRID=(*s)[1:2])
    if n_elements(x) eq 1 then return,ptr_new()
    y=GRID_TPS((*TieO)[0,*],(*TieO)[1,*],(*TieI)[1,*], START=[0,0],DELTA=gs,NGRID=(*s)[1:2])
    if n_elements(y) eq 1 then return,ptr_new()
    
    ; Remark: this gives the same result
    ; dx=GRID_TPS((*TieO)[0,*],(*TieO)[1,*],(*TieI)[0,*]-(*TieO)[0,*], START=[0,0],DELTA=gs,NGRID=(*s)[1:2])
    ; x=rebin(lindgen((*s)[1]),(*s)[1:2],/sample)+dx
    ; dy=GRID_TPS((*TieO)[0,*],(*TieO)[1,*],(*TieI)[1,*]-(*TieO)[1,*], START=[0,0],DELTA=gs,NGRID=(*s)[1:2])
    ; y=rebin(lindgen(1,(*s)[2]),(*s)[1:2],/sample)+dy
endelse

; Trigrid:
;    EXTRA = bounds to extrapolate for points outside the triangles
;    MISSING = -1 all points outside the triangles will get value -1 (remark: -1 is outside the image)
; Interpolate (cfr. WarpSpaDistRun):
;    MISSING = 0 for all points outside the image

; Trigrid+Interpolate = WARP_TRI

; Visualize
;hrange=[1660,2047l]
;vrange=[1660,2047l]
;xs=(hrange[1]-hrange[0]+1)
;ys=(vrange[1]-vrange[0]+1)
;        
;dx=GRID_TPS((*TieO)[0,*],(*TieO)[1,*],(*TieI)[0,*]-(*TieO)[0,*], START=[0,0],DELTA=gs,NGRID=(*s)[1:2])
;dy=GRID_TPS((*TieO)[0,*],(*TieO)[1,*],(*TieI)[1,*]-(*TieO)[1,*], START=[0,0],DELTA=gs,NGRID=(*s)[1:2])
;
;ind=where((*TieO)[0,*] ge hrange[0] and (*TieO)[0,*] le hrange[1] and $
;          (*TieO)[1,*] ge vrange[0] and (*TieO)[1,*] le vrange[1])
;xx=reform((*TieO)[0,ind])
;yy=reform((*TieO)[1,ind])
;zzx=reform((*TieI)[0,ind]-(*TieO)[0,ind])
;zzy=reform((*TieI)[1,ind]-(*TieO)[1,ind])
;
;m=min(dx[hrange[0]:hrange[1],vrange[0]:vrange[1]])
;plt=surface(dx[hrange[0]:hrange[1],vrange[0]:vrange[1]],findgen(xs)+hrange[0],findgen(ys)+vrange[0],$
;            xthick=2,ythick=2,zthick=2,xtickname=replicate('',20),ytickname=replicate('',20),$
;            ztickname=replicate('',20),xstyle=1,ystyle=1,color='bebebe'xl,$
;            xTICKINTERVAL=xs/5,yTICKINTERVAL=ys/5,zTICKINTERVAL=10)
;for i=0,n_elements(xx)-1 do $
;    line=POLYLINE(xx[[i,i]],yy[[i,i]],[m,zzx[i]],/data,target=plt,color='ff0000'xl,thick=5)
;
;m=min(dy[hrange[0]:hrange[1],vrange[0]:vrange[1]])
;plt=surface(dy[hrange[0]:hrange[1],vrange[0]:vrange[1]],findgen(xs)+hrange[0],findgen(ys)+vrange[0],$
;            xthick=2,ythick=2,zthick=2,xtickname=replicate('',20),ytickname=replicate('',20),$
;            ztickname=replicate('',20),xstyle=1,ystyle=1,color='bebebe'xl,$
;            xTICKINTERVAL=xs/5,yTICKINTERVAL=ys/5,zTICKINTERVAL=10)
;for i=0,n_elements(xx)-1 do $
;    line=POLYLINE(xx[[i,i]],yy[[i,i]],[m,zzy[i]],/data,target=plt,color='ff0000'xl,thick=5)

p=PTR_NEW({x:x,y:y})
return,p
end;function WarpSpaDistPrep
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro WarpSpaDistRun,tiff,struc,tiffvalid=tiffvalid

; Prepare warp
badpixels=n_elements(tiffvalid) eq 0
if badpixels then begin
    ; tiff == tiffvalid
    (*tiff)=float(*tiff)
    
    missing=0
endif else begin 
    if ptr_valid(tiffvalid) then $
        WarpSpaDistRun,tiffvalid,struc
        
    type=size(*tiff,/type)
    bconvert=type ne 4 or type ne 5
    if bconvert then *tiff=float(*tiff)
    
    missing=!values.F_NAN
endelse

; Warp
*tiff=INTERPOLATE(*tiff,(*struc).x,(*struc).y, MISSING=missing,cubic=-0.5)

; Handle missing data
if badpixels then begin
    ; tiff == tiffvalid
;    (*tiff) eq= 1 ; too stringent
    (*tiff) gt= 0.9
endif else begin
    valid=finite(*tiff)
    ind=where(~valid,ct)
    if ct ne 0 then begin
        (*tiff)[ind]=0
        if not ptr_valid(tiffvalid) then tiffvalid=ptr_new(temporary(valid)) $
        else (*tiffvalid) and= temporary(valid)
    endif
    ;if bconvert then *tiff=TypeConvert(temporary(*tiff),type)
endelse

end;pro WarpSpaDistRun
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
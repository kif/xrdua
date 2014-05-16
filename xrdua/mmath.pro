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

function CheckOrthogonal,M
return,abs(max(transpose(M)-invert(M,/double),/abs)) lt 1e-5 
end;function CheckOrthogonal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RotDecompose,Min,inversion=inversion

if ~array_equal(size(Min,/dim),[3,3]) then stop,'RotDecompose - Matrix must be 3x3.'
if ~CheckOrthogonal(Min) then stop,'RotDecompose - Matrix must be orthogonal in order to decompose it in rotations.'

M=Min
inv=-identity(3)

for inversion=0,1 do begin
    if inversion then M=Min##inv
    for i=0,7 do begin
        ; asin -> [-90,90]
        flipa=(i and 1) ne 0
        flipb=(i and 2) ne 0
        flipg=(i and 4) ne 0
    
        sinb=double(M[0,2])
        b=asin(-sinb)
        if flipb then b=!dpi-b
        cosb=cos(b)
        sina=M[1,2]/cosb
        ;cosa=M[2,2]/cosb
        sing=M[0,1]/cosb
        ;cosg=M[0,0]/cosb
        
        a1=asin(sina)
        if flipa then a1=!dpi-a1
        ;a2=acos(cosa)
        g1=asin(sing)
        if flipg then g1=!dpi-g1
        ;g2=acos(cosg)
        
        if abs(max((rotz(-g1)##roty(-b)##rotx(-a1))[0:2,0:2]-M,/abs)) lt 1e-5 then break
    endfor
    
    if i lt 8 then break
endfor

if i eq 8 then stop,'RotDecompose - decomposition failed.'

return,[a1,b,g1] ; angles for change of basis matrices
end;function RotDecompose
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ispoweroftwo,x
 return, ((x ne 0) and ((x and (not x + 1)) eq x))
end;function ispoweroftwo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mtotal,array,dim,_extra=ex
if n_elements(dim) eq 1 then begin
    s=size(array,/dim)
    if dim gt n_elements(s) or dim le 0 then return,array else begin
        if s[dim-1] eq 1 then return,array else begin
            ret=total(array,dim,_extra=ex)
            s[dim-1]=1
            ret=reform(ret,/overwrite)
            return,ret
        endelse
    endelse
endif else return,total(array,_extra=ex)
end;function mtotal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function percentile,data,pct=pct,dimension=dimension,even=even,indices=indices
; data=findgen(4,4) & data[1,1]=!values.F_NAN & data[1,2]=-!values.F_INFINITY
; print,median(data,dim=1),percentile(data,pct=50,dim=1)
; print,median(data,dim=2),percentile(data,pct=50,dim=2)

; Median
if ~keyword_set(pct) then pct=50
pct=0>pct<100
if pct eq 50 and ~arg_present(indices) then return,median(data,dimension=dimension,even=even)

; Percentile of all data
if ~keyword_set(dimension) then begin
    n=total(~finite(data,/nan),/int)
    indices=ceil(pct/100.*(n-1))
    indices=(sort(data))[indices]
    return,data[indices]
endif

s=dimsize(data,3)
if dimension gt n_elements(s) then return,data

s2=shrinkarray(s,dimension-1)
;ret=make_array(s2,value=data[0])

n=total(~finite(data,/nan),dimension,/int)
indices=ceil(pct/100.*(n-1))

; Works only for 2 dimensions for now
case dimension of
1:    begin
    for j=0l,s[1]-1 do $
        for k=0l,s[2]-1 do begin
            indices[j,k]=(sort(data[*,j,k]))[indices[j,k]]
            ;ret[j,k]=data[indices[j,k],j,k]
        endfor
    ret=data[indices,rebin(lindgen(s[1]),s[1],s[2],/sample),rebin(lindgen(1,s[2]),s[1],s[2],/sample)]
    endcase
2:    begin
    for i=0l,s[0]-1 do $
        for k=0l,s[2]-1 do begin
            indices[i,k]=(sort(data[i,*,k]))[indices[i,k]]
            ;ret[i,k]=data[i,indices[i,k],k]
        endfor
    ret=data[rebin(lindgen(s[0]),s[0],s[2],/sample),indices,rebin(lindgen(1,s[2]),s[0],s[2],/sample)]
    endcase
3:    begin
    for i=0l,s[0]-1 do $
        for j=0l,s[1]-1 do begin
            indices[i,j]=(sort(data[i,j,*]))[indices[i,j]]
            ;ret[i,j]=data[i,j,indices[i,j]]
        endfor
    ret=data[rebin(lindgen(s[0]),s[0],s[1],/sample),rebin(lindgen(1,s[1]),s[0],s[1],/sample),indices]
    endcase
endcase

return,ret
end;function percentile,data
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SphericalToCartesian,v_
; Each row is a vector

v=double(v_)
sinpolar=sin(v[2,*])
return,[v[0,*]*cos(v[1,*])*sinpolar,v[0,*]*sin(v[1,*])*sinpolar,v[0,*]*cos(v[2,*])]

end;function SphericalToCartesian
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CartesianToSpherical,v_
; Each row is a vector

v=double(v_)
s=dimsize(v)
if n_elements(s) ne 2 then R=sqrt(total(v*v,1)) $
else R=reform(sqrt(total(v*v,1)),1,s[1])

azimuth=atan(v[1,*],v[0,*]) ; between -180 and 180
ind=where(azimuth lt 0,ct)
if ct ne 0 then azimuth[ind]+=2*!dpi ; between 0 and 360

polar=acos(v[2,*]/R) ; between 0 and 180

return,[R,azimuth,polar]
end;function CartesianToSpherical
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InnerProduct,v,w,M=M,ROW=ROW
bROW=keyword_set(ROW)

sv=dimsize(v,2)
sw=dimsize(w,2)
if ~array_equal(sv,sw) then stop,'Vectors not compatible.'
n=sv[1-bROW]
nvec=sv[bROW]

if n_elements(M) eq 0 then begin
    ;M=identity(n,type=size(v,/type))
    if bROW then ret=v##transpose(w) $
    else ret=transpose(v)##w
endif else begin
    sM=dimsize(M,2)
    if sM[0] ne n or sM[1] ne n then stop,'Metric tensor not compatible.'
    
    if bROW then ret=v##M##transpose(w) $
    else ret=transpose(v)##M##w
endelse

return,ret[lindgen(nvec) * (nvec + 1L)]

end;function InnerProduct
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function VectorNorm,v,M=M,ROW=ROW
return,sqrt(InnerProduct(v,v,M=M,ROW=ROW))
end;function VectorNorm
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CrossProduct,v,w,C=C,ROW=ROW
; ROW does not apply to C
bROW=keyword_set(ROW)

sv=dimsize(v,2)
sw=dimsize(w,2)
if ~array_equal(sv,sw) then stop,'Vectors not compatible.'
n=sv[1-bROW]
nvec=sv[bROW]
if n ne 3 then stop,'Cross product only implemented for 3 dimensions.'

if n_elements(C) eq 0 then begin
    C=identity(n,type=size(v,/type))
endif else begin
    sM=dimsize(C,2)
    if sM[0] ne n or sM[1] ne n then stop,'Change of basis not compatible.'
endelse

if bROW then coeff=[v[1,*]*w[2,*]-v[2,*]*w[1,*],v[2,*]*w[0,*]-v[0,*]*w[2,*],v[0,*]*w[1,*]-v[1,*]*w[0,*]] $
else coeff=transpose([[v[*,1]*w[*,2]-v[*,2]*w[*,1]],[v[*,2]*w[*,0]-v[*,0]*w[*,2]],[v[*,0]*w[*,1]-v[*,1]*w[*,0]]])
coeff=reform(coeff,n*nvec,/overwrite)
coeff=rebin(coeff,n*nvec,n,/sample)

coeff=coeff*C[lindgen(n*nvec) mod n,*]

ret=[[chunktotal(coeff[*,0],replicate(n,nvec))],$
    [chunktotal(coeff[*,1],replicate(n,nvec))],$
    [chunktotal(coeff[*,2],replicate(n,nvec))]]
if bROW then ret=transpose(ret)

return,ret
end;function CrossProduct
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TripleProduct,u,v,w,C=C,ROW=ROW
; Transform coordinates to the canonical basis
if n_elements(C) ne 0 then begin
    bROW=keyword_set(ROW)

    sv=dimsize(v,2)
    sw=dimsize(w,2)
    if ~array_equal(sv,sw) then stop,'Vectors not compatible.'
    n=sv[1-bROW]
    nvec=sv[bROW]
    if n ne 3 then stop,'Cross product only implemented for 3 dimensions.'
    
    sM=dimsize(C,2)
    if sM[0] ne n or sM[1] ne n then stop,'Change of basis not compatible.'
    
    if bROW then begin
        return,innerproduct(matrix_multiply(u,C,/atranspose),CrossProduct(matrix_multiply(v,C,/atranspose),matrix_multiply(w,C,/atranspose)))
    endif else begin
        return,innerproduct(C##u,CrossProduct(C##v,C##w))
    endelse
    
endif else begin
    return,innerproduct(u,CrossProduct(v,w,ROW=ROW),ROW=ROW)
endelse

end;function TripleProduct
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SpannedVolume,u,v,w,C=C,ROW=ROW
; ROW does not apply to C
bROW=keyword_set(ROW)

sv=dimsize(v,2)
sw=dimsize(w,2)
if ~array_equal(sv,sw) then stop,'Vectors not compatible.'
n=sv[1-bROW]
nvec=sv[bROW]
if n ne 3 then stop,'Cross product only implemented for 3 dimensions.'

if n_elements(C) eq 0 then begin
    if bROW then Cuvw=[[[transpose(u)]],$
                        [[transpose(v)]],$
                        [[transpose(w)]]]$
    else Cuvw=[[[u]],[[v]],[[w]]]
endif else begin
    sM=dimsize(C,2)
    if sM[0] ne n or sM[1] ne n then stop,'Metric tensor not compatible.'
    
    if bROW then Cuvw=[[[matrix_multiply(u,C,/atranspose)]],$
                        [[matrix_multiply(v,C,/atranspose)]],$
                        [[matrix_multiply(w,C,/atranspose)]]]$
    else Cuvw=[[[C##u]],[[C##v]],[[C##w]]]
    
endelse

Cnew=transpose(Cuvw,[2,1,0])

ret=make_array(nvec,type=size(v,/type))
for i=0l,nvec-1 do ret[i]=sqrt(abs(determ(transpose(Cnew[*,*,i])##Cnew[*,*,i])))

return,ret    
end;function SpannedVolume
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TransformMetricTensor,M,C
; C: columns are the coordinates of the new basis vectors with respect to the old ones
if ~array_equal(dimsize(C,2),dimsize(M,2)) then stop,'Wrong dimensions.'

return,transpose(C)##M##C
end;function TransformMetricTensor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Combinations,k,n,error=error
; All possibilities (order does not matter) to select k elements from a set of n elements.
; http://compprog.wordpress.com/2007/10/17/generating-combinations-1/

error=1b
subsets=lindgen(k)
subset=subsets

if k gt n then return,subsets

repeat begin
    i = k - 1l
    subset[i]++
    while ((i ge 0)?(subset[i] ge (n - k + 1 + i)):0b) do subset[--i]++

    ; Combination (n-k, n-k+1, ..., n) reached
    ; No more combinations can be generated
    bstop=subset[0] gt (n - k)
    
    ; subset now looks like (..., x, n, n, n, ..., n).
    ; Turn it into (..., x, x + 1, x + 2, ...)
    if ~bstop then begin
        for i=i+1,k-1 do subset[i] = subset[i - 1] + 1
        subsets=[[subsets],[subset]]
    endif
endrep until bstop

s=dimsize(subsets,2)
if round(factorial(n)/(factorial(k)*factorial(n-k))) ne s[1] then return,subsets

error=0b
return,subsets
end;function Combinations
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Permutations,n,error=error
; http://compprog.wordpress.com/2007/10/05/generating-permutations-1/

error=1b
tupples=lindgen(n)+1
tupple=tupples

repeat begin
    i = n - 1l
    tupple[i]++
    while ((i ge 0)?(tupple[i] gt n):0b) do begin
        tupple[i] = 1
        i--
        if (i ge 0) then tupple[i]++
    endwhile
    
    if i ge 0 then $
        if n_elements(uniq(tupple,sort(tupple))) eq n then tupples=[[tupples],[tupple]]
endrep until i lt 0

s=dimsize(tupples,2)
if round(factorial(n)) ne s[1] then return,tupples

error=0b
return,tupples

end;function Permutations
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function kPermutations,k,n,error=error
; All possibilities (order does matter) to select k elements from a set of n elements.

subsets=Combinations(k,n,error=error)
if error then return,subsets
tupples=Permutations(k,error=error)-1
if error then return,subsets

s1=dimsize(subsets,2)
s2=dimsize(tupples,2)
nsubsets=s1[1]
ntupples=s2[1]

ct=nsubsets*ntupples
subsets=rebin(temporary(subsets),k,ct,/sample)
for i=0,nsubsets-1 do subsets[*,i*ntupples:(i+1)*ntupples-1]=(subsets[*,i*ntupples])[tupples]

return,subsets
end;function kPermutations
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mPCA,A,loadings=loadings,contribution=contribution,standardize=standardize,error=error

; INPUT and OUTPUT:
;
; A (input, n x m):         each column is a variable, each row is an element
; PC (returned):             each column is a principal component
; contribution (output):     contribution to the variance for each PC
; standardize (input):         use standard scaling conventions
; loadings (output, m x m):    each row i is a scaled eigenvector which coordinates
;                           give the contribution of each variable to the PCi
;                           each column j gives the contribution of variable j
;                           to all PCs
; 
; The m-dimensional row vectors of PC are the m-dimensional row vectors
; of A with respect to a translated and rotated coordinate system. The 
; transformation if chosen so that the X-axis lies along the greatest variance
; is the scatter plot, the Y-axis along the greatest perpendicular, etc.
; 
; 
; BACKGROUND:
; 
; 
; Eucledian vector space V over a field R:
;     basis B: {b1,b2,b3}
;     vector v: v=x.b1+y.b2+z.b3
;     representation of v w.r.t. B: vB=[x,y,z] or [v]B=[x,y,z]^T
;     standard representation F of V with respect to B: F:V->R^3:v->[v]B
; 
; Change of basis: linear transformation
;     second basis C: {c1,c2,c3}
;     The matrix of representations of the basis vectors bi w.r.t. C: MB_C=[[b1]C,[b2]C,[b3]C]
;     This transformation matrix can be used to get a vector as a
;     function of the new basis C : [v]C = MB_C.[v]B = MC_B^(-1).[v]B
;     
; Change of basis: orthogonal transformation
;    MC_B^(-1) = MC_B^T   => MB_C = MC_B^T
; 
; Singular value decomposition (cfr. SolveLinSystem): A = U.W.V^T
;   A (n x m): n experiments (rows), m variables (columns)
;   U (n x m): orthogonal matrix
;   W (m x m): diagonal matrix with singular values on the diagonal
;   V (m x m): orthogonal matrix
;   
;         A^T.A.V = V.W.W
;         R.V = V.W.W
;     <=> (R-W[i,i]^2.I).V=0
;     <=> V[i,*] are the eigenvectors of R=A^T.A with eigenvalues W[i,i]^2
; 
; Empirical covariance matrix of A: C = (Acen^T.Acen)/(n-1) = R/(n-1)
;     Acen (n x m): A with mean of each column zeroed (Acen[i,*]=A[i,*]-mean(A[i,*]))
;     C (m x m): covariance matrix
;
; Principal Component Analysis is an orthogonal transformation that transforms the 
; data to a new coordinate system such that the greatest variance by any projection of the 
; data comes to lie on the first coordinate (called the first principal component), the
; greatest perpendicular variance on the second coordinate, and so on. 
; 
; First transform the data so that the mean is at the origin
;    A (n x m): each of the n rows is an m-dimensional vector
;    Acen (n x m): sum of all m-dimensional vectors is the zero-vector
;    
; Secondly, rotate the coordinate system so that the greatest variance lies along the X-axis,
; the greatest perpendicular variance along the Y-axis, etc. The basis vectors of the new coordinate system
; with respect to the old are given by the eigenvectors V of the covariance matrix of A. The
; orthogonal transformation is given by:
; 
;        PC^T = V^T . Acen^T
;    <=> PC = Acen.V
;    
;    Acen (n x m): each row is a vector with respect to the old basis (which has already an origin shift included)
;    PC (n x m): each row is a vector with respect to the new basis
;    V (m x m): each column is a new basis vector with respect to the old system and is an eigenvector of R
;
; Other definitions are:
;     loadings = V.W (i.e. scaled eigenvectors)
;     component scores = Acen.V.W (i.e. PCs calculated from scaled eigenvectors)


error=1b
CATCH, Error_status
if Error_status ne 0 then return,0

dim = SIZE(A, /DIM)

; Centered A:
tmp = MOMENT(A, DIMENSION=2, Maxmoment=2)
Acen = A-rebin(tmp[*,0],dim,/sample) ; subtract mean from each column (i.e. each variable)

; Standardize A:
if keyword_set(standardize) then begin
    ; Use this when variables don't have the some units
    ind=where(tmp[*,1] eq 0,ct)
    if ct ne 0 then tmp[ind,1]=1
    Acen/=rebin(sqrt(tmp[*,1]),dim,/sample) ; divide each column by its standard deviation
endif

; Covariance matrix (which is symmetric)
R=MATRIX_MULTIPLY(Acen, Acen, /BTRANSPOSE)/(dim[1]-1) ; Used C instead of R, why?
; MATRIX_MULTIPLY(A, B, /BTRANSPOSE) = A # B^T = B^T ## A

; Eigenvalues and eigenvectors of R (sorted by descending eigenvalue)
eigenvalues = EIGENQL(R, EIGENVECTORS = eigenvectors) ; each row is an eigenvector

; Scaled eigenvectors or loadings = V.W
;     For each eigenvector v: v.v = 1
;     We can scale v so that v.v = eigenvalue
;    These scaled eigenvectors are called loadings.
W=identity(sqrt(eigenvalues),value=eigenvalues[0]*0)
eigenvectors=MATRIX_MULTIPLY(eigenvectors, W)

; Principal components: PC = Acen.V   (eigenvectors = V^T)
PC = MATRIX_MULTIPLY(eigenvectors, Acen, /ATRANSPOSE) ; each column is a principal component
; MATRIX_MULTIPLY(A, B, /ATRANSPOSE) = A^T # B = B ## A^T

if ARG_PRESENT(loadings) then loadings=temporary(eigenvectors) ; each row is an eigenvector

; Contribution of each PC to the total variance:
if ARG_PRESENT(contribution) then $
    contribution=eigenvalues/total(eigenvalues)

error=0b
return,PC

end;function mPCA
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function chisqr_test,chi2,dof,pvalue_cutoff=pvalue,pout=p
; NULL HYPOTHESIS H0:      DISTRIBUTION FITS THE DATA
; ALTERNATE HYPOTHESIS HA: DISTRIBUTION DOES NOT FIT THE DATA

if ~keyword_set(pvalue) then pvalue=0.95
p=CHISQR_PDF(chi2,dof) ; probability that chi2 is less than the chi2 given as input
return, p le pvalue ; 1=accept H0, 0=reject H0
end;function chisqr_test
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function entropy,p
;calculates the entropy of a distribution p
q=p/total(p)
ind=where(q ne 0,ct)
if ct eq 0 then return,0.
q=q[ind]
base=2
return,-total(q*alog(q)/alog(base))
end;function entropy
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function joint_entropy,x
;calculates the joint entropy of an m-D array x(n,m) with quantization q
s=size(x,/dim)
return,entropy(histogram(long(total(x*rebin((max(x)+1)^lindgen(1,s[1]),s[0],s[1],/sample),2))))
end;function joint_entropy
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mutual_information,x,y,scale=scale,normalized=normalized
if keyword_set(scale) then begin
    mix=min(x,max=max)
    miy=min(y,max=may)
    x1=bytscl(x,min=mix,max=max)
    y1=bytscl(y,min=miy,max=may)
endif else begin
    x1=x
    y1=y
endelse
n=n_elements(x1)

hx=entropy(histogram(x1))
hy=entropy(histogram(y1))
;hxy=joint_entropy([[reform(x1,n)],[reform(y1,n)]])
hxy=entropy(hist_nd([reform(x1,1,n),reform(y1,1,n)],1))

if keyword_set(normalized) then return,(hx+hy)/hxy $
else return,hx+hy-hxy
end;function mutual_information
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pearson_correlate,x,y,interval=interval
; Pearson product-moment correlation coefficient
; (http://davidmlane.com/hyperstat/B62223.html)
r=correlate(x,y) ; Pearson correlation coefficient (no normal sampling distribution)

; "Fisher's z transformation": convert Pearson's r in a normally distributed zvalue
;     zvalue = 0.5*alog((1+r)/(1-r))
; with a standard deviation
;     sigmaz = 1/sqrt(n_elements(x)-3.0)
; Remark: this only works when x and y are coming from a normal distribution
zvalue = alog((1.+r)/(1.-r))/2.            ; E(z)
nx=n_elements(x)
sigmaz = 1./sqrt(nx-3.0)        ; E(sqrt(VAR(zmean)))

; Confidence interval for r
;    zvalue +/- cutoff.sigmaz  <-- 99% of the time, this interval will contain the real z
P = 0.99 ; 99% confidence interval
cutoff = GAUSS_CVF((1-P)/2.)
interval=zvalue+[-cutoff,cutoff]*sigmaz
; convert interval back to pearson correlation
interval=exp(2.*interval)
interval=(interval-1)/(interval+1)

; Hypothesis testing:
; 1. Suppose we measure property x of a population with unknown distribution function P(x).
; 2. The null hypothesis (H0) makes an assumption on the distribution function, e.g. a normal distribution
;    with mean=0 and variance=1.
; 3. The alternative hypothesis (H1) states that the measured x doesn't come from the distribution of proposed by H0.
; 4. Compute the p-value of the measured x for P(x): probability that measure a value which is "more extreme" than x
; 5. Significance level alpha: the smallest acceptable probability for a value to be concidered "not extreme"
; 6. There are two possible outcomes:
;         a. p-value < alpha: x is concidered too extreme, H0 is rejected and H1 is accepted
;         b. p-value >= alpha: not conclusion possible, except for the fact that we can't reject H0

; Hypothesis testing for Pearson correlation:
; tests whether r is significantly different from zero (i.e. H0 rejected)
; 
;r0=0
;zvalue0=alog((1.+r0)/(1.-r0))/2.=0
;xx=(zvalue-zvalue0)/sigmaz
xx=zvalue/sigmaz

; H0: x has normal distribution with mean=0 and variance=1, or in other words r=r0=0
; H1: r!=0
pvalue1=1-gauss_pdf(abs(xx))

; Same test using the student-t distribution
pvalue2=1-T_PDF(abs(r)*sqrt(nx-2)/sqrt(1-r*r),nx-2)

return,[r,2*pvalue1,2*pvalue2] ; 2x for the two-tailed probability
end;function pearson_correlate
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function testnormal,a,KS=KS,KUIPER=KUIPER
; NOT IMPLEMENTED!!!: See IMSL

; Default: Pearson's chi-square test
; Other:     KS = Kolmogorov-Smirnov test
;            KUIPER = Kuiper test

type=keyword_set(KUIPER)?2:(keyword_set(KS)?1:0)

;case type of
;0:    begin
;    df = n_elements(b)-1
;    z = total((a - b)^2.0 / b)
;    pvalue = 1 - chisqr_pdf(z, df)
;    return, [z, prob] 
;    endcase
    
return,0b
end;function testnormal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function hashCRC32,data

; Hash function type: cyclic redundancy check (CRC), namely CRC-32 (IEEE 802.3) MPEG-2 specifications
; polynomial:
; x^32 + x^26 + x^23 + x^22 + x^16 + x^12 + x^11 + x^10 + x^8 + x^7 + x^5 + x^4 + x^2 + x + 1

polynomial = '04c11db7'XUL
remainder  = 'ffffffff'XUL
for word=0,n_elements(data)-1 do begin
    datum = ulong(data[word])
    remainder xor= datum
    for bit=0,31 do begin
      if ( remainder and '80000000'XUL ) ne 0 then $
        remainder = ishft(remainder , 1) xor polynomial $
      else remainder = ishft(remainder , 1)
    endfor
endfor
return,remainder
end;function hashCRC32
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sign,x
; Keep type
; Nan and -Nan => sign=1
y=x
y[*]=1-2*(x lt 0)
return,y
end;function sign
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function atantan,xin,ain,divide=divide
  ; Evaluate y = atan(a.tan(x))
  ; or       y = atan(tan(x)/a)
  ; 
  ; Problem: atan maps everything to [-90,90]
  ; 
  ; a>=0: |y-x| < 90
  ; a<0: transfer sign of a to x

  x = xin*sign(ain)
  a = abs(ain)
  
  if keyword_set(divide) then y=atan(tan(x)/a) $
  else y=atan(a*tan(x))
  
  ; The difference cannot be greater or equal to 90 degrees
  diff = y-x
  ind = where(abs(diff) ge !dpi/2,ct)
  if ct ne 0 then y[ind] += round(-diff[ind]/!dpi)*!dpi

  return,y
end;function atantan
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

;function atantan,xin,ain,divide=divide
;; Conserve quadrant of y in
;; y=atan(a.tan(x)) or y=atan(tan(x)/a)
;
;; Transfer sign of a to x
;x=xin*sign(ain)
;a=abs(ain)
;
;; Check input
;if a eq 1 then return,x
;
;double=size(x,/type) eq 5
;pi=double?!dpi:!pi
;
;; y in [-90,90]
;if keyword_set(divide) then y=atan(tan(x)/a) $
;else y=atan(a*tan(x))
;
;; Quadrant encoding: ...,-1=[-180,-90[,0=[-90,0[,1=[0,90],2=]90,180],...
;; Source quadrant:
;qy=(y gt 0)
;; Destination quadrant: 
;qx=(x gt 0) + fix(x*2/pi)
;
;; Map source to destination quadrant
;y+=pi/2*(qx-qy)
;
;; Round off errors
;cutoff=abs(atan(tan(acos(a))))
;cutoff+=2*floatpres(double=double)
;diff=abs(y-x)
;ind=where(diff ge cutoff,ct)
;if ct ne 0 then y[ind]=x[ind]
;
;return,y
;end;function atantan
;;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Equal0ModZ,x,prec
; x =? 0 (mod Z)
tmp = abs(x) mod 1 ; from ]-inf,+inf[ to [0,1[
tmp mod= 1-prec ; from [0,1[ to [0,1-prec[
tmp le= prec
return,tmp
end;function Equal0ModZ
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function seed
return,nothing;randomu(SYSTIME(1,/SECONDS))
end;function seed
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function randxy,x,y,n,long=long
; n random (integer) number between x and y
if keyword_set(long) then return,round(randomu(seed(),n,/uniform)*(y-x+1)+x-0.5) $
else return,randomu(seed(),n,/uniform)*(y-x)+x
end;function randxy
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function randlong
; random, positive integer number
return,randomu(seed(),/long)
end;function randlong
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddNoise,y
for i=0l,n_elements(y)-1 do $
    y[i]=RANDOMN(seed(),1,Poisson=y[i])
end;pro AddNoise
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function digamma,xin
;http://www.codecogs.com/d-ox/maths/special/gamma/psi.php

n=n_elements(xin)
x=reform(xin,1,n)
B=[-12.0,120.0,-252.0,240.0]
nB=n_elements(B)
pow=rebin(-2-2*findgen(nB),nB,n)
B = rebin(B,nB,n)
;  =  -2n/B_{2n} , even Bernoulli numbers

; Make x greater than 6 and keep digamma correction
; digamma(z+1) = digamma(z) + 1/z
; e.g. z=4.5 => digamma(5.5) = digamma(6.5) - 1/4.5 - 1/5.5
a=6
val=where(x le a,ct)
if ct ne 0 then begin
    nx = ceil(a-x)
    nx += nx eq 0
    nxmax = max(nx)
    if nxmax le 0 or n le 0 then stop
    val=rebin(indgen(nxmax),nxmax,n)

    bool=val lt rebin(nx,nxmax,n)

    val += rebin(x,nxmax,n)
    val = -total(bool/val,1)
    x += nx
endif else val=0.

; Get digamma of x+...(>6) without the Bernoulli sum
val += alog(x)-0.5/x

; Add Bernoulli sum
x=rebin(x,nB,n)
val += total(x^pow/B,1)

return, (n eq 1)?val[0]:val
end;function digamma
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function dergamma,x
return,gamma(x)*digamma(x)
end;function dergamma,x
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fround,in,n
; n: number of digits after the decimal point
nrnd=10.^n
return,round(in*nrnd,/L64)/nrnd
end;function fround
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ffloor,in,n
; n: number of digits after the decimal point
nrnd=10.^n
return,floor(in*nrnd,/L64)/nrnd
end;function ffloor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fceil,in,n
; n: number of digits after the decimal point
nrnd=10.^n
return,ceil(in*nrnd,/L64)/nrnd
end;function fceil
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sigdigits,value,sd,round=round
; sd: position of the least significant digit
;         value = ABC.DEF
;        A => sd=2, B => sd=1, C => sd=0, D => sd=-1, E = sd=-1, F = sd=-1
; value: truncate (or round) value at the least significant digit
; return: significant digits

bround=keyword_set(round)
if sd lt 0 then begin
    ; ABC.DEF -> truncate/round at ABC.D or ABC.DE or ABD.DEF or ...
    exp1=10ull^(-sd+bround) ; 1 digit too much when rounding
    exp2=10d^sd
endif else begin
    ; ABC.DEF -> truncate/round at ABC or AB0 or A00
    exp1=10d^(-sd+bround) ; 1 digit too much when rounding
    exp2=10ull^sd
endelse 

; Get significant digits (1 digit too much when rounding)
digits=string(abs(value)*exp1,format='(I0)')

; When rounding, remove the last digit and adapt the previous digit
if bround then begin
    n=strlen(digits)
    if fix(strmid(digits,n-1,1)) ge 5 then begin
        ; round
        if n eq 1 then digits='1' $
        else begin
            digits=strmid(digits,0,--n) ; remove last digit
            
            ; Add 1 to last digit (might need propagation)
            for i=n-1,0,-1 do begin
                rdigit=fix(strmid(digits,i,1))+1
                if rdigit eq 10 then insdigit='0' $
                else insdigit=string(rdigit,format='(I0)')
                strput,digits,insdigit,i
                if rdigit ne 10 then break
            endfor
            if i eq -1 then digits='1'+digits
            
        endelse
    endif else begin
        ; truncate
        if n eq 1 then digits='0' $
        else digits=strmid(digits,0,n-1) ; remove last digit
    endelse
endif

; Get value from significant digits + appropriate power of 10
value=sign(value)*long64(digits)*exp2

return,digits
end;function sigdigits
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rounderror,value,error,format=format
; Round "value" according to its error "error"
; Examples:
;    IDL> value=0.268995 & error=0.00196546 & print,rounderror(value,error)
;    0.2690(19)
;    IDL> value=0.268995 & error=0.00196546 & print,rounderror(value,error,format=1)
;    0.2690 ± 0.0019

; Absolute value of error
error=abs(error)

; Truncate the error according to the 10% rule:
;        (error - truncated_error)/error < 10%
;        <=> (ZX.YYYYY - ZX) < ZX.YYYYY/10    ; Z={1,...,9}, X={0,...,9}, Y={0,...,9}
;        <=> 0.YYYYY < Z.XYYYYY
;    The last expression is always true because 0.YYYYY < 1 and Z.XYYYYY >= 1,
;    so one can forget about the Y's and just check whether to use Z or ZX for
;    the truncated error. So the truncated error has 1 or 2 significant figures.
;
; Significant digit: ABC.DEF
;    A => 2, B => 1, C => 0, D => -1, E = -1, F = -1
;
if error ne 0 then begin
    sdfirst_error = floor(alog10(error)) ; first significant digit of the error
    
    ; Truncate error at the first significant digit
    errortrunc1=error
    errordigits=sigdigits(errortrunc1,sdfirst_error)

    if (error-errortrunc1)/error ge 0.10 then begin
        ; Keep two significant digits
        sdleast = sdfirst_error-1 ; last significant digit of the error and the value
        errordigits=sigdigits(error,sdleast)
    endif else begin
        ; Keep one significant digit
        sdleast = sdfirst_error ; last significant digit of the error and the value
        error = errortrunc1
    endelse
endif else begin
    sdfirst_error = 0 ; first significant digit of the error
    sdleast = 0 ; last significant digit of the error and the value
    errordigits = '0'
endelse

; Round the value to the least significant digit
valuedigits=sigdigits(value,sdleast,/round)

; Format output
if n_elements(format) eq 0 then format=0
case format of
1:    begin
    if sdleast eq 0 then output=valuedigits+' '+stringr(177b)+' '+errordigits $
    else $
    if sdleast lt -strlen(valuedigits)-4 or sdleast gt 0 then begin
        ; Use scientific notation, add a decimal point and add to exponent
        dec=1
        
        sdleast_=sdleast+strlen(valuedigits)-dec
        if sdleast_ eq 0 then exp='' else exp='E'+string(sdleast_,format='(I0)')
        valuedigits=strmid(valuedigits,0,dec)+'.'+strmid(valuedigits,dec)+exp
        
        sdleast_=sdleast+strlen(errordigits)-dec
        if sdleast_ eq 0 then exp='' else exp='E'+string(sdleast_,format='(I0)')
        errordigits=strmid(errordigits,0,dec)+'.'+strmid(errordigits,dec)+exp

        output=valuedigits+' '+stringr(177b)+' '+errordigits
    endif else begin
        ; Add insignificant zero's so that the negative exponent vanishes (maximal 5 insignificant zeros)
        dec=strlen(valuedigits)+sdleast
        if dec le 0 then begin
            valuedigits=string(replicate(byte('0'),-dec+1))+valuedigits
            dec=1
        endif
        valuedigits=strmid(valuedigits,0,dec)+'.'+strmid(valuedigits,dec)
        
        dec=strlen(errordigits)+sdleast
        if dec le 0 then begin
            errordigits=string(replicate(byte('0'),-dec+1))+errordigits
            dec=1
        endif
        errordigits=strmid(errordigits,0,dec)+'.'+strmid(errordigits,dec)
        
        output=valuedigits+' '+stringr(177b)+' '+errordigits
    endelse
    
    endcase
else: begin
    ; Default
    
    if sdleast eq 0 then output=valuedigits+'('+errordigits+')' $
    else $
    if sdleast lt -strlen(valuedigits)-4 or sdleast gt 0 then begin
        ; Use scientific notation
        if strlen(valuedigits) gt strlen(errordigits) then begin
            ; Add a decimal point and add to exponent
            dec=1
            sdleast+=strlen(valuedigits)-dec
            if sdleast eq 0 then exp='' else exp='E'+string(sdleast,format='(I0)')
            output=strmid(valuedigits,0,dec)+'.'+strmid(valuedigits,dec)+'('+errordigits+')'+exp
        endif else begin
            if sdleast eq 0 then exp='' else exp='E'+string(sdleast,format='(I0)')
            output=valuedigits+'('+errordigits+')'+exp
        endelse
    endif else begin
        ; Add insignificant zero's so that the negative exponent vanishes (maximal 5 insignificant zeros)
        dec=strlen(valuedigits)+sdleast
        if dec le 0 then begin
            valuedigits=string(replicate(byte('0'),-dec+1))+valuedigits
            dec=1
        endif
        output=strmid(valuedigits,0,dec)+'.'+strmid(valuedigits,dec)+'('+errordigits+')'
    endelse
    
    endelse
endcase
if sign(value) eq -1 then output='-'+output

return,output
end;function rounderror
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function defaulterror,val

if val eq 0 then return,val*0

str=string(val)

expon=0

p=strpos(str,'e',0)
if p ne -1 then begin
    expon+=fix(strmid(str,p+1))
    str=strmid(str,0,p-1)
endif

p=strpos(str,'.',0)
if p ne -1 then str=strmid(str,p+1)

ind=where(byte(str) ne 48b,ct)
if ct eq 0 then return,val*0
expon-=ind[ct-1]+1

return,10.^expon
end;function defaulterror
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function fconvol,array,kernel
return,fft(fft(array,-1)*kernel,1)
end;function fconvol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BlockDiag,AA,sind=sind,onlycolumns=onlycolumns,recursivecall=recursivecall

; AA: mxn array
; onlycolumns: swap only columns
; sind: new indices
; recursivecall: internal use only

; Swapping columns?
brows=keyword_set(recursivecall)
if brows then begin
    A=transpose(AA)
    sind=transpose(sind)
endif else begin
    A=AA
    sind=lindgen(size(AA,/dim))
endelse

; Sort columns
s=size(A,/dim)
indswap=indgen(s[0])
coli=indswap
for j=0l,s[1]-1 do begin
    ind=where(A[indswap,j] ne 0,ct,comp=comp,ncomp=ncomp)
    nindswap=n_elements(indswap)
    if ct ne 0 and ct ne nindswap then begin
        ind=indswap[ind]
        comp=indswap[comp]
        coli[s[0]-nindswap:*]=[ind,comp]
        indswap=comp
    endif
endfor
A=A[coli,*]
sind=sind[coli,*]

; Sort rows
if brows then begin
    A=transpose(A)
    sind=transpose(sind)
endif else if not keyword_set(onlycolumns) then A=BlockDiag(A,sind=sind,/recursivecall)

return,A
end;function BlockDiag
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UniqueRows,T,bunique=bunique

; T: row vectors
complex=size(T,/type)
complex=complex eq 6 or complex eq 9
s=Dimsize(T,2)
if s[0] eq 1 or s[1] eq 1 then begin
    bunique=replicate(1b,s[1]>1)
    return,T
endif

; Eliminate vectors which are a multiple of one of the others:
; scalar product equals to the product of the norms
if complex then begin
    T_=conj(T)
    epsm=2
endif else T_=T
norm=sqrt(total(abs(T_*T),1,/pres))

bunique=bytarr(s[1])
ind=where(norm ne 0,ct)
if ct eq 0 then return,T[*,0]*0
i0=ind[0]
Tout=T[*,i0]
normout=norm[i0]
bunique[i0]=1
n=1

for i=0l,s[1]-1 do begin
    bool=norm[i] ne 0
    for j=0l,n-1 do $
        bool and= ~floatequal(abs(total(T_[*,i]*Tout[*,j],/pres)),normout[j]*norm[i],epsm=epsm)

    if bool then begin
        Tout=[[Tout],[T[*,i]]]
        normout=[normout,norm[i]]
        bunique[i]=1b
        n++
    endif
endfor

return,Tout
end;function UniqueRows
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RealUniqueRows,T,bunique=bunique

; T: row vectors
complex=size(T,/type)
complex=complex eq 6 or complex eq 9
s=Dimsize(T,2)
if s[0] eq 1 or s[1] eq 1 then begin
    bunique=replicate(1b,s[1]>1)
    return,T
endif

; Eliminate vectors which are a multiple of one of the others:
; scalar product equals to the product of the norms
if complex then begin
    T_=conj(T)
    epsm=2
endif else T_=T

bunique=bytarr(s[1])
i0=0
Tout=T[*,i0]
bunique[i0]=1
n=1

for i=0l,s[1]-1 do begin
    bool=1b
    for j=0l,n-1 do $
        bool and= ~floatequal(T_[*,i],Tout[*,j],epsm=epsm)

    if bool then begin
        Tout=[[Tout],[T[*,i]]]
        bunique[i]=1b
        n++
    endif
endfor

return,Tout
end;function RealUniqueRows
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UniqueCols,T,bunique=bunique

; T: column vectors
complex=size(T,/type)
complex=complex eq 6 or complex eq 9
s=Dimsize(T,2)
if s[0] eq 1 or s[1] eq 1 then begin
    bunique=replicate(1b,s[0]>1)
    return,T
endif

; Eliminate vectors which are a multiple of one of the others:
; scalar product equals to the product of the norms
if complex then begin
    T_=conj(T)
    epsm=2
endif else T_=T
norm=sqrt(total(abs(T_*T),2,/pres))

bunique=bytarr(s[0])
ind=where(norm ne 0,ct)
if ct eq 0 then return,T[0,*]*0
i0=ind[0]
Tout=T[i0,*]
normout=norm[i0]
bunique[i0]=1
n=1

for i=0l,s[0]-1 do begin
    bool=norm[i] ne 0
    for j=0l,n-1 do $
        bool and= ~floatequal(abs(total(T_[i,*]*Tout[j,*],/pres)),normout[j]*norm[i],epsm=epsm)

    if bool then begin
        Tout=[Tout,T[i,*]]
        normout=[normout,norm[i]]
        bunique[i]=1b
        n++
    endif
endfor

return,Tout
end;function UniqueCols
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RealUniqueCols,T,bunique=bunique

; T: column vectors
complex=size(T,/type)
complex=complex eq 6 or complex eq 9
s=Dimsize(T,2)
if s[0] eq 1 or s[1] eq 1 then begin
    bunique=replicate(1b,s[0]>1)
    return,T
endif

; Eliminate vectors which are a multiple of one of the others:
; scalar product equals to the product of the norms
if complex then begin
    T_=conj(T)
    epsm=2
endif else T_=T

bunique=bytarr(s[0])
i0=0
Tout=T[i0,*]
bunique[i0]=1
n=1

for i=0l,s[0]-1 do begin
    bool=1b
    for j=0l,n-1 do $
        bool and= ~floatequal(T_[i,*],Tout[j,*],epsm=epsm)

    if bool then begin
        Tout=[Tout,T[i,*]]
        bunique[i]=1b
        n++
    endif
endfor

return,Tout
end;function RealUniqueCols
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function singular,A,Double = Double
CATCH, Error_status
if Error_status ne 0 then return,1b
tmp = LA_INVERT(A, Status=Status, Double = Double)
return,Status gt 0
end;function singular
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro QRfact,AT, ipvt, rdiag, Acol_norm,returnR=returnR,nopivot=nopivot

; QR decomposition:
;   A ## P = Q ## R
;     A(m x n) m >= n
;     Q(m x m) orthogonal (Q^T.Q=I)
;     R(m x n) upper triangular (the strict lower trapezoidal part = 0)
;     P(n x n) permutation matrix (also orthogonal)
;
; Input:
;   AT(n x m): transposed of A
;
; Output:
;   AT(n x m): the strict lower triangular part of AT contains the strict
;      upper triangular part of R, and the upper trapezoidal
;      part of AT contains the reflectors that generate Q (cfr. comments at the end)
;   ipvt(1 x n); column sort indices for A
;   rdiag(1 x n): diagonal elements of R
;   Acol_norm(1 x n): column norms of A

sAT = size(AT)
mcol = sAT[1]
nrow = sAT[2]

; Conver type to float or double
type = size(AT,/type)
if type eq 5 or type eq 9 or type eq 14 or type eq 15 then type=5 else type=4
SetType,AT,type
double=type EQ 5
fpres = floatpres(double=double)

; Compute the initial column norms of A and initialize arrays
Acol_norm=sqrt(total(AT^2,1))
rdiag = Acol_norm
wa = rdiag
ipvt = lindgen(nrow)
bpivot=not keyword_set(nopivot)

; Reduce AT to R with householder transformations
minmn = nrow<mcol
for j = 0L, minmn-1 do begin

    ; Bring the column of largest norm in A into the pivot position.
    ; Exchange via pivot only, the real exchange happens later.
    if bpivot then begin
        rmax = max(rdiag[j:*],kmax)
        if kmax ne 0 then begin
            kmax +=j
            temp = ipvt[j]
            ipvt[j] = ipvt[kmax]
            ipvt[kmax] = temp
            rdiag[kmax] = rdiag[j]
            wa[kmax] = wa[j]
        endif
    endif

    ; Compute the householder transformation to reduce the jth
    ; column of A to a multiple of the jth unit vector
    lj     = ipvt[j]
    ajj    = AT[j:*,lj]
    ajnorm = sqrt(total(ajj^2))
    if ajnorm eq 0 then rdiag[j]=0 else begin

        if AT[j,lj] LT 0 then ajnorm = -ajnorm
        ajj/=ajnorm
        ajj[0]++
        AT[j,lj] = ajj

        ; Apply the transformation to the remaining columns of A
        ; and update the norms
        for k=j+1, nrow-1 do begin ; loop over columns of A
            lk = ipvt[k]
            ajk = AT[j:*,lk]
            if AT[j,lj] NE 0 then $
                AT[j,lk] = ajk - ajj * total(ajk*ajj)/AT[j,lj]

            if bpivot and rdiag[k] NE 0 then begin
                temp = AT[j,lk]/rdiag[k]
                rdiag[k] *= sqrt((1.-temp^2) > 0)
                temp = rdiag[k]/wa[k]
                if 0.05D*temp*temp LE fpres then begin
                    rdiag[k] = sqrt(total(AT[j+1:*,lk]^2))
                    wa[k] = rdiag[k]
                endif
            endif
        endfor

        rdiag[j] = -ajnorm
    endelse
endfor

if keyword_set(returnR) then begin
    ; R-matrix
    rmat = make_array(nrow,mcol,type=type)
    for j = 1, minmn-1 do rmat[j,0:j-1] = AT[0:j-1,ipvt[j]]
    idiag = lindgen(minmn)
    rmat[idiag, idiag] = rdiag

    AT=temporary(rmat)
endif

;; Test: qmat ## rmat must give A ## P i.e. transpose(AT[*,ipvt])    
;; R-matrix
;rmat = make_array(nrow,mcol,type=type)
;for j = 1, minmn-1 do rmat[j,0:j-1] = AT[0:j-1,ipvt[j]]
;idiag = lindgen(minmn)
;rmat[idiag, idiag] = rdiag
;
;; Q-matrix
;ident = identity(mcol,double=double,float=~double)
;qmat = ident
;for i = 0L, nrow-1 do begin
;    v = AT[*,ipvt[i]]
;    if i GT 0 then v[0:i-1] = 0
;    qmat ##= (ident - 2*(v # v)/total(v * v))
;endfor
;    
;; Q.R == transpose(AT[*,ipvt])
;;Attention: AT changed!!!
;print,qmat##rmat

end;pro QRfact
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro QRsolve, r, ipvt, diag, qtb, x, sdiag
; Get X (least squares solution) for which:
;     A(m x n) ## X(n x 1) = B(m x 1) and D(n x n) ## X(n x 1) = 0 (D diagonal matrix)
;
; If X = P ## Z then
;    A ## P ## Z = B and D ## P ## Z = 0
;    Q ## R ## Z = B and D ## P ## Z = 0 ; QR fact.
;    R ## Z = Q^(T) ## B and D ## P ## Z = 0 ; Q orthogonal
;
; P^T(A^T.A+D.D).P=S^T.S where S an upper triangular matrix
;
; Input:
;     r(n x n): the full lower triangle must contain the full upper
;               triangle of the matrix R.
;     ipvt(1 x n): from QRfact
;     diag(1 x n): diagonal elements of matrix D
;     qtb(1 x n): must contain the first n elements of the vector Q^T ## B
;
; Output:
;     r(n x n): on output the full upper triangle is unaltered, and the
;         strict lower triangle contains the strict upper triangle
;         (transposed) of the upper triangular matrix S.
;     x(1 x n): contains the least squares solution.
;     sdiag(1 x n): contains the diagonal elements of the upper triangular matrix S.

sz = size(r)
mcol = sz[1]
nrow = sz[2]
delm = lindgen(nrow) * (mcol+1) ;; Diagonal elements of r

; copy r and (q transpose)*b to preserve input and initialize s.
; in particular, save the diagonal elements of r in x.
for j = 0L, nrow-1 do $
    r[j:nrow-1,j] = r[j,j:nrow-1]
x = r[delm]
wa = qtb
; Below may look strange, but it's so we can keep the right precision
zero = qtb[0]*0.
half = zero + 0.5
quart = zero + 0.25

; Eliminate the diagonal matrix d using a givens rotation
for j = 0L, nrow-1 do begin
    l = ipvt[j]
    if diag[l] ne 0 then begin
        sdiag[j:*] = 0
        sdiag[j] = diag[l]

        ; The transformations to eliminate the row of d modify only a
        ; single element of (q transpose)*b beyond the first nrow, which
        ; is initially zero.

        qtbpj = zero
        for k = j, nrow-1 do begin
            if sdiag[k] ne 0 then begin
                if abs(r[k,k]) LT abs(sdiag[k]) then begin
                      cotan  = r[k,k]/sdiag[k]
                      sine   = half/sqrt(quart + quart*cotan*cotan)
                      cosine = sine*cotan
                endif else begin
                      tang   = sdiag[k]/r[k,k]
                      cosine = half/sqrt(quart + quart*tang*tang)
                      sine   = cosine*tang
                endelse

                ; Compute the modified diagonal element of r and the
                ; modified element of ((q transpose)*b,0).
                r[k,k] = cosine*r[k,k] + sine*sdiag[k]
                temp = cosine*wa[k] + sine*qtbpj
                qtbpj = -sine*wa[k] + cosine*qtbpj
                wa[k] = temp

                ; Accumulate the transformation in the row of s
                if nrow GT k+1 then begin
                    temp = cosine*r[k+1:nrow-1,k] + sine*sdiag[k+1:nrow-1]
                    sdiag[k+1:nrow-1] = -sine*r[k+1:nrow-1,k] + cosine*sdiag[k+1:nrow-1]
                    r[k+1:nrow-1,k] = temp
                endif
            endif
        endfor

    endif
    sdiag[j] = r[j,j]
    r[j,j] = x[j]
endfor

; Solve the triangular system for z.  If the system is singular
; then obtain a least squares solution
nsing = nrow
wh = where(sdiag EQ 0, ct)
if ct GT 0 then begin
    nsing = wh[0]
    wa[nsing:*] = 0
endif

if nsing GE 1 then begin
    wa[nsing-1] = wa[nsing-1]/sdiag[nsing-1] ; Degenerate case
    for j=nsing-2,0,-1 do begin
        sum = total(r[j+1:nsing-1,j]*wa[j+1:nsing-1])
        wa[j] = (wa[j]-sum)/sdiag[j]
    endfor
endif

; Permute the components of z back to components of x
x[ipvt] = wa

end;pro QRsolve
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReducedRowEchelonFrom,M
s=DimSize(M,2)>1
mcol=s[0]
nrow=s[1]

; Three operations that don't change the system
; 1. The order of two rows can be interchanged.
; 2. A row can be multiplied by a nonzero constant.
; 3. A row can be replaced by the sum of that row and 
;    a nonzero multiple of any other row.

tolerance=floatpres(M)*(mcol>nrow)*MAX(TOTAL(ABS(M), 1, /pres))
j=0

bcol=bytarr(nrow)
for i=0l,mcol-1 do begin
    ; Find the pivot row
    ma=max(M[i,j:*],pivot,/abs)
    pivot += j
      
    if abs(ma) le tolerance then begin
        ;Skip column i, making sure the approximately zero terms are
        ;actually zero.
          M[i,j:*] = 0
    endif else begin
        ;Swap current row and pivot row (1)
        M[i:*,[pivot, j]] = M[i:*,[j, pivot]]

        ;Normalize pivot row (2)
        M[i:*,j] /= M[i,j]
    
        ;Eliminate the current column (3)
        bcol[*]=1b
        bcol[j]=0b
        ridx = where(bcol,nj)
        mi=mcol-i
        M[i:*,ridx] -= rebin(M[i,ridx],mi,nj,/s) * rebin(M[i:*,j],mi,nj,/s)

        ;Check if done
        if ++j eq nrow then break
    endelse

endfor

ind=where(M lt tolerance,ct)
if ct ne 0 then M[ind]=0

end;pro ReducedRowEchelonFrom
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RowEchelonFrom,M
s=DimSize(M,2)>1
mcol=s[0]
nrow=s[1]

; Three operations that don't change the system
; 1. The order of two rows can be interchanged.
; 2. A row can be multiplied by a nonzero constant.
; 3. A row can be replaced by the sum of that row and 
;    a nonzero multiple of any other row.

tolerance=floatpres(M)*(mcol>nrow)*MAX(TOTAL(ABS(M), 1, /pres))
j=0

for i=0l,mcol-1 do begin
    ; Find the pivot row
    ma=max(M[i,j:*],pivot,/abs)
    pivot += j
      
    if abs(ma) le tolerance then begin
        ;Skip column i, making sure the approximately zero terms are
        ;actually zero.
          M[i,j:*] = 0
    endif else begin
        ;Swap current row and pivot row (1)
        M[i:*,[pivot, j]] = M[i:*,[j, pivot]]

        ;Normalize pivot row (2)
        M[i:*,j] /= M[i,j]
    
        ;Eliminate the current column (3)
        nj=nrow-1-j
        if nj ne 0 then begin
            ridx=j+1+lindgen(nj)
            mi=mcol-i
            M[i:*,ridx] -= rebin(M[i,ridx],mi,nj,/s) * rebin(M[i:*,j],mi,nj,/s)
        endif
    
        ;Check if done
        if ++j eq nrow then break
    endelse
endfor

ind=where(M lt tolerance,ct)
if ct ne 0 then M[ind]=0

end;pro RowEchelonFrom
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ElimSolve,AA,BB

; This doesn't work (choose 1 for variable that might be 0...)
stop

A=AA
B=BB

s=DimSize(A,2)>1
mcol=s[0]
nrow=s[1]

; Default solution: choose 1 for all variables
X=make_array(mcol,value=A[0]*0+1)
bfree=bytarr(mcol)
    
repeat begin
    ; Variables on each row
    h=histogram([total(A ne 0,1,/int)],min=0,max=mcol,binsize=1,rev=rev)
        
    ; Row: a0.x0+a1.x1+... = b
    nvar=1
    while h[nvar] eq 0 and nvar lt mcol do nvar++

    if h[nvar] ne 0 then begin
        ; Take one of the rows with nvar variables:
        ;  i between 0 and h[nvar]-1
        i=0
        row=rev[rev[nvar]+i]
        col=where(A[*,row] ne 0,ct)
        
        ; Take the variable with the highest coefficient
        tmp1=max(A[col,row],tmp2,/abs)
        col=col[tmp2]
        
        ; If this is the only variable in the row: solve
        ; If there are otherwise, set this one to 1.
        if nvar eq 1 then r=B[row]/A[col,row] else begin
            r=1
            bfree[col]=1b
        endelse
        X[col]=r
        
        ; Eliminate this variable from all equations
        row=where(A[col,*] ne 0,ct)
        for ii=0,ct-1 do begin
            ; e.g.:      a0.x + a1.y + a2.z = b
            ;         and  x = r
            ;              0.x + a1.y + a2.z = b - a0.r
            B[row[ii]]-= r*A[col,row[ii]]
            A[col,row[ii]] = 0
        endfor
    endif
endrep until array_equal(h[1:*],0)

; Did we find a solution?
error=total(B ne 0,/int) ne 0
if error then begin
    ind=where(bfree,ct)
    for i=0l,ct-1 do begin
        ; Do something
        ;ElimSolve(AA,BB,error=error)
    endfor
endif

return,X
end;function ElimSolve
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GJSolve,A,B,M=M,error=error
; Solve through Gauss–Jordan elimination
; ElimSolve doesn't work yet!!!!!

; Augmented matrix
s=DimSize(A,2)>1
mcol=s[0]
nrow=s[1]
nvar=mcol
M=[A,reform(B,1,nrow)]
mcol++

; Create row echelon form and solve
if (mcol-1) eq nrow and floatequal(B,0) then begin
    ; Homogeneous system
    if singular(A) then begin
        ; Infinite number of solutions: find a nonzero one
        RowEchelonFrom,M
        X=ElimSolve(M[0:mcol-2,*],M[mcol-1,*],error=error)
    endif else begin
        X=make_array(nvar,value=A[0]*0)
        error=0b
    endelse
endif else begin
    ReducedRowEchelonFrom,M

    n=nvar<nrow
    X=M[mcol-1,0:n-1]
    if nrow lt nvar then X=[[X],[replicate(1,nvar-nrow)]]
    
    ; Check error
    AA=M[0:mcol-2,*]
    ind=indgen(n)
    AA[ind,ind]--
    ind=where(AA gt tolerance,ct)
    error=ct ne 0
endelse

return,X

end;function GJSolve
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MatrixColumnSpace,A,rank
; Get orthonormal basis of column space
s=DimSize(A,2)>1
mcol=s[0]
nrow=s[1]

type=size(A,/type)
complex=type EQ 6 or type EQ 9
double=type EQ 5 or type EQ 9
pres=floatpres(double=double)*(mcol>nrow)*MAX(TOTAL(ABS(A), 1, /pres))

if mcol eq 1 then begin
    U=A
    ind=where(abs(U) le pres,ct,ncomp=rank)
    if ct ne 0 then U[ind]=0
    rank <= 1
    return,U
endif

if complex then LA_SVD, A, W, U, V $
else SVDC, A, W, U, V
indRange=where(abs(W) gt pres,rank)

if rank eq 0 then return,A[0,*]*0

U=U[indRange,*]
ind=where(abs(U) le pres,ct)
if ct ne 0 then U[ind]=0
            
return,U
end;function MatrixColumnSpace
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro mSVDC, Ain, W, U, V, NullSpace=NullSpace, nullity=nullity

; Make A floating point: single or double precision
A=Ain
type1=size(A,/type)
complex=type1 EQ 6 or type1 EQ 9
double=keyword_set(Double) or type1 EQ 5 or type1 EQ 9
if double then type=5 else type=4
if complex then begin
  if double then begin
    A = [[DOUBLE(A), -DOUBLE(IMAGINARY(A))], $
      [DOUBLE(IMAGINARY(A)), DOUBLE(A)]]
  endif else begin
    A = [[FLOAT(A), -FLOAT(IMAGINARY(A))], $
      [FLOAT(IMAGINARY(A)), FLOAT(A)]]
  endelse
endif else begin
  SetType,A,type
endelse

s=DimSize(A,2)>1
mcol=s[0]
nrow=s[1]

; Singular value decomposition
SVDC, A, W, U, V

; Null space
pres=floatpres(Double=Double)*(mcol>nrow)*MAX(TOTAL(ABS(A), 1, /pres))
indNull=where(abs(W) le pres,nullity);,comp=indRange,ncomp=rank
if nullity ne 0 then NullSpace=V[indNull,*] else NullSpace=make_array(1,mcol,type=type)  ; each column is a basis vector of NullSpace

; Complex Null space
if complex then begin
    n=mcol/2
    
    NullSpace=double ? DCOMPLEX(NullSpace[*,0:n-1], NullSpace[*,n:*]) : COMPLEX(NullSpace[*,0:n-1], NullSpace[*,n:*])
    NullSpace=MatrixColumnSpace(NullSpace,nullity); nullity is too high (why?)
endif

end;pro mSVDC
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MatrixNullSpace,A,nullity
; Get orthonormal basis of column space

mSVDC, A, NullSpace=NullSpace, nullity=nullity

return,NullSpace
end;function MatrixNullSpace
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mSVSOL, U, W, V, Bin, A, Trivial=Trivial, Residual=Residual

; Make B floating point: single or double precision
B=Bin
type1=size(A,/type)
type2=size(B,/type)
complex=type1 EQ 6 or type1 EQ 9 or type2 EQ 6 or type2 EQ 9
double=keyword_set(Double) or type1 EQ 5 or type1 EQ 9 or type2 EQ 5 or type2 EQ 9
if double then type=5 else type=4
if complex then begin
  if double then begin
    B = [DOUBLE(B), DOUBLE(IMAGINARY(B))]
  endif else begin
    B = [FLOAT(B), FLOAT(IMAGINARY(B))]
  endelse
endif else begin
  SetType,B,type
endelse

s=DimSize(A,2)>1
mcol=s[0]
nrow=s[1]

; Particular solution
X=SVSOL(U, W, V, reform(B))

; Trivial solution?
Trivial=floatequal(X,0.)

; Residuals
if ARG_PRESENT(Residual) then Residual=reform(A##reform(X,1,mcol)-B)

; Complex solution
if complex then begin
    n=mcol/2
    X=double ? DCOMPLEX(X[0:n-1], X[n:*]) : COMPLEX(X[0:n-1], X[n:*])
    
    if ARG_PRESENT(Residual) then begin
        n=nrow/2
        Residual=double ? DCOMPLEX(Residual[0:n-1], Residual[n:*]) : COMPLEX(Residual[0:n-1], Residual[n:*])
    endif
endif

X=reform(X,DimSize(Bin),/overwrite)

return,X
end;function mSVSOL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SolveLinSystem,Ain,Bin,NonTrivial=NonTrivial,NullSpace=NullSpace,nullity=nullity,Residual=Residual,Double=Double

; Find a solution for A(n x m).XX(m x 1) = B(n x 1)
;   => XX = X + Null(A)
;       X: particular solution
;       Kernel or Null space of A: Null(A) = {X: A.X = 0}
;       Nullity: dimension of the Null space
;   => There is no solution if B doesn't belong to the column space of A
;       Column space of A (aka range or image): space spanned by the columns of A
;       Rank: dimension of the Column space (= number of independent rows and columns)
;   Rank-nullity theorem: rank(A) + nullity(A) = m (=number of columns)
;   
; Eigenvalues (λ) and eigenvectors (v) of A when A(n x n):
;         A.v = λ.v
;     <=>    (A-λI).v=0
;     
;     - E_λ(A): eigenspace corresponding to λ, spanned by the corresponding eigenvectors
;     - dim(E_λ(A)) = m_λ : dimension of eigenspace corresponding to λ = multiplicity of λ
;     - E_0(A) = 0 <=> A.v=0: eigenspace corresponding to 0 is the kernel of A
;       E_1(A) = 0 <=> (A-I).v=0: eigenspace corresponding to 1 is the kernel of A-1
;       ...
;     
;   
; Singular value decomposition: A = U.W.V^T      (A^T = V.W.U^T)
;   A (n x m)
;   U (n x m): orthogonal matrix
;   W (m x m): diagonal matrix with singular values on the diagonal
;   V (m x m): orthogonal matrix
;       
; eigenvalues and eigenvectors:
;         R = A^T.A = V.W.W.V^T
;     <=>    R.V = V.W.W
;     <=> (R-W[i,i]^2.I).V=0
;     <=> V[i,*] are the eigenvectors of R with eigenvalues W[i,i]^2
;     
;         Q = A.A^T = U.W.W.U^T
;     <=>    Q.U = U.W.W
;     <=> (Q-W[i,i]^2.I).U=0
;     <=> U[i,*] are the eigenvectors of Q with eigenvalues W[i,i]^2
;     
; sum of its eigenvalues = trace of a matrix
; <=> sum(W[i,i]^2) = sum(R[i,i]) = sum(Q[i,i])
;   
; column and null space:
;       A.V = U.W
;   <=> A.V[i,*] = U[i,*].W[i,i]
;       U^T.A = W.V^T
;   <=> (U[i,*])^T.A = W[i,i].(V[i,*])^T
;         V[i,*] is the right singular vector of A for singular value W[i,i]
;         U[i,*] is the left singular vector of A for singular value W[i,i]
;       All V[i,*] corresponding to singular value W[i,i]==0: span Null space
;       All U[i,*] corresponding to singular value W[i,i]!=0: span Column space
;
;    V = [Ker(A),KerC(A)]
;        Ker(A): span null space (columns of V with W==0)
;        KerC(A): complement to the null space (columns of V with W!=0)
;    U = [LIm(A),Im(A)]
;        LIm(A): linear combinations of Im(A) (columns of U with W==0)
;        Im(A): span column space (columns of U with W!=0)
;
; A = U.W.V^T
;   = [LIm(A),Im(A)].W.[Ker(A),KerC(A)]^T
;   = [LIm(A),Im(A)].([0,KerC(A).Wnz])^T
;       Wnz: all W!=0
;       KerC(A).Wnz=[W0.V0,W1.V1,...]
;       Im(A)=[U0,U1,...]
;   = Im(A).[W0.V0,W1.V1,...]^T
;
; => We reduced the matrices: A(n x m) = U(n x m-k).W(m-k x m-k).V(m-k x m)
;    with m-k non-zero singular values
; => SVDC: k=0    (i.e. some singular values can be zero)
;    LA_SVD: k=(m-n)>0 (i.e. avoid obvious zero-valued singular values, but some can still be zero)

; Convert input to float or double
; For a complex system:
; Are | -Aim   Xre   Bre
; ---------- . --- = ---
; Aim | Are    Xim   Bim

A=Ain
B=Bin
type1=size(A,/type)
type2=size(B,/type)
complex=type1 EQ 6 or type1 EQ 9 or type2 EQ 6 or type2 EQ 9
double=keyword_set(Double) or type1 EQ 5 or type1 EQ 9 or type2 EQ 5 or type2 EQ 9
if double then type=5 else type=4
if complex then begin
  if double then begin
    A = [[DOUBLE(A), -DOUBLE(IMAGINARY(A))], $
      [DOUBLE(IMAGINARY(A)), DOUBLE(A)]]
    B = [DOUBLE(B), DOUBLE(IMAGINARY(B))]
  endif else begin
    A = [[FLOAT(A), -FLOAT(IMAGINARY(A))], $
      [FLOAT(IMAGINARY(A)), FLOAT(A)]]
    B = [FLOAT(B), FLOAT(IMAGINARY(B))]
  endelse
endif else begin
  SetType,A,type
  SetType,B,type
endelse

s=DimSize(A,2)>1
mcol=s[0]
nrow=s[1]

; Singular value decomposition
SVDC, A, W, U, V

; Null space
pres=floatpres(Double=Double)*(mcol>nrow)*MAX(TOTAL(ABS(A), 1, /pres))
indNull=where(abs(W) le pres,nullity);,comp=indRange,ncomp=rank
if nullity ne 0 then NullSpace=V[indNull,*] else NullSpace=make_array(1,mcol,type=type)  ; each column is a basis vector of NullSpace

; Particular solution
X=SVSOL(U, W, V, reform(B))

; Non trivial solution?
NonTrivial=~(floatequal(X,0.) and nullity eq 0)

; Residuals
if ARG_PRESENT(Residual) then Residual=reform(A##reform(X,1,mcol)-B)

; Complex solution
if complex then begin
    n=mcol/2
    X=double ? DCOMPLEX(X[0:n-1], X[n:*]) : COMPLEX(X[0:n-1], X[n:*])

    NullSpace=double ? DCOMPLEX(NullSpace[*,0:n-1], NullSpace[*,n:*]) : COMPLEX(NullSpace[*,0:n-1], NullSpace[*,n:*])
    NullSpace=MatrixColumnSpace(NullSpace,nullity); nullity is too high (why?)
    
    if ARG_PRESENT(Residual) then begin
        n=nrow/2
        Residual=double ? DCOMPLEX(Residual[0:n-1], Residual[n:*]) : COMPLEX(Residual[0:n-1], Residual[n:*])
    endif
endif

return,X
end;function SolveLinSystem
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ScaleRoundColumns,V,pres
s=DimSize(V,2)
for i=0l,s[0]-1 do begin
    ind=where(abs(V[i,*]) gt pres,ct)
    if ct eq 0 then continue
    tmp=V[i,*]/V[i,ind[0]]
    if total(abs(tmp-round(tmp)) gt pres,/int) eq 0 then V[i,*]=round(tmp)
endfor
end;pro ScaleRoundColumns
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function JordanDecompHelper,A,Eval,nEval,malg,Evec,Jordan,pres

; Returns:
;     -1 -> error
;     0  -> OK
;    1  -> can't find the Jordan canonical form

CATCH, Error_status
IF Error_status NE 0 THEN return,-1

; Some constants
s=DimSize(A,2)
n=s[0] ; A[nxn]
diag=lindgen(n)*(n+1) ; diagonal

type=size(Eval,/type)
complex=type EQ 6 or type EQ 9
Jordan=make_array(size=size(A))

; Find (generalized) eigenvectors
Joff=0
for i=0l,nEval-1 do begin
    ; Characteristic matrix: M=A-λ.I
    M=A
    M[diag]-=Eval[i]
    Mkeep=M
    
    ; Loop over k in (A-λ.I)^k
    dimeigenspace=0
    k=1
    repeat begin
        
        ; Decompose M=(A-λ.I)^k
        if complex then begin
            LA_SVD, M, W, U, V, status=stat
            if stat ne 0 then return,-1
        endif else SVDC, M, W, U, V

        ; Null space of (A-λI)^k, i.e. solve (A-λI)^k.v = 0
        tmp=abs(max(W,/abs))>1
        indNull=where(abs(W)/tmp le pres,nullity)
        if nullity ne 0 then begin
            V=V[indNull,*] ; null space
            
            if dimeigenspace eq 0 then begin
                ; Geometric multiplicity of λ greater than the algebraic (happens when 2 eigenvalues are actually the same)
                if nullity gt malg[i] then return,1
                
                ; Scale+round eigenvectors when possible
                ScaleRoundColumns,V,pres

                ; First members of the eigenspace
                eigenspace = V
                dimeigenspace = nullity

                ; Prepare Jordan chains
                nchains=dimeigenspace ; Number of chains
                chaindepthmax=malg[i]-dimeigenspace+1 ; maximum jordan chain length
                chains=intarr(chaindepthmax+1,dimeigenspace) ; Each row is a Jordan chain
                chains[0,*]=1 ; First column: chain count
                chains[1,*]=indgen(dimeigenspace) ; Second column: eigenspace column indices, pointing to the first vectors of each Jordan chain
        
            endif else begin
                ; V contains the whole eigenspace and possibly something more
                ; For example:
                ;     eigenspace=[v1,v2] and V=[2.v1+3.v2,v1+5.v2+v3,v1+2.v3]
                
                ; (A-λI)^(k-1).v=0            ; combination of previous eigenvectors
                ; (A-λI)^(k-1).v=vprev        ; combination of previous and new eigenvectors
                ; For example:
                ;    eigenspace=[v1,v2] and (A-λI)^(k-1).V=[0,(A-λI)^(k-1).v3,2.(A-λI)^(k-1).v3]
                ScaleRoundColumns,V,pres
                Vprev=Mprev##V
                ScaleRoundColumns,Vprev,pres
                
                ; Keep only unique columns
                ; For example:
                ;    eigenspace=[v1,v2] and (A-λI).V=[(A-λI).v3]
                Vprev=UniqueCols(Vprev,bunique=bunique)
                V=V[where(bunique),*]
                ind=where(total(abs(Vprev) gt pres,2,/int) ne 0,ct)
                beigenspace=bytarr(dimeigenspace)
                dimeigenspace0=dimeigenspace
                for j=0l,ct-1 do begin
                    ; Find previous vector in Jordan chain
                    ; For example: (A-λI).v3=3.v2
                    Vadd=Vprev[ind[j],*]
                    indlink=where(abs(reform(Vadd)##eigenspace) gt pres,ct)
                    indlink=indlink[0]
                    if ct eq 0 then continue ; no previous vector (can this happen?)
                    if indlink ge dimeigenspace0 then continue ; no previous vector (can this happen?)
                    if beigenspace[indlink] then continue ; already found a vector
                    if chains[0,indlink] eq chaindepthmax then continue ; Jordan chain is full (can this happen?)
                    
                    ; Add vector to eigenspace and update chains
                    beigenspace[indlink]=1
                    eigenspace=[eigenspace,V[ind[j],*]]
                    chains[0,indlink]++
                    chains[chains[0,indlink],indlink]=dimeigenspace++
                endfor
                
            endelse
        endif
        
        ; M=(A-λ.I)^(k+1)
        Mprev=M
        M ##= Mkeep
        k++

    endrep until dimeigenspace eq malg[i] or k gt chaindepthmax
    
    ; Edit Jordan diagonal
    ind=Joff+indgen(malg[i])
    Jordan[ind,ind]=Eval[i]
    
    ; Group the Jordan chains: largest chains first
    ind=bsort(chains[0,*],/rev)
    chains=chains[*,ind]
    
    ; Group the Jordan chains: get vector indices
    ind=lonarr(dimeigenspace)
    off=0
    for j=0l,nchains-1 do begin
        ; vectors in this chain, starting with the normal eigenvector
        indchain=chains[1:chains[0,j],j]
        ind[off:off+chains[0,j]-1]=indchain
        
        ; Use highest eigenvector and generate lower vectors in chain
        for k=chains[0,j]-1,1,-1 do $
            eigenspace[indchain[k-1],*]=Mkeep##eigenspace[indchain[k],*]
        
        ; edit super-diagonal: 1's
        if chains[0,j] gt 1 then begin
            tmp=indgen(chains[0,j]-1)+off+Joff
            Jordan[tmp+1,tmp]=1
        endif
        off+=chains[0,j]
    endfor
    eigenspace=eigenspace[ind,*]

    ; Add eigenspace for this eigenvalue to the list of generalized eigenvectors of A
    if i eq 0 then Evec=eigenspace $
    else Evec=[Evec,eigenspace]
    
    Joff+=malg[i]
endfor

; ----Jordan canonical form----
if complex then J=la_invert(Evec)##A##Evec $
else J=invert(Evec)##A##Evec

ind=where(abs(J-Jordan) gt pres*10,ct)
if ct ne 0 then return,1
        
return,0
end;function JordanDecompHelper
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function mHQR,A,pres=pres,error=error,Double=Double
error=1b

; ----Constants----
s=DimSize(A,2)
if s[0] ne s[1] then return,0
n=s[0] ; A[nxn]
diag=lindgen(n)*(n+1) ; diagonal

; ----Get eigenvalues of A----
type1=size(A,/type)
if type1 eq 6 or type1 eq 9 then begin
    Eval = LA_HQR(LA_ELMHES(A,Q,Double = Double,PERMUTE_RESULT=permute,SCALE_RESULT=scale),$
                            Double = Double,PERMUTE_RESULT=permute,status=status)
    error = status ne 0
    if error then return,0
endif else Eval = HQR(ELMHES(A,Double = Double),Double = Double)

; ----Check type of the eigenvalues----
type2=size(Eval,/type)
double=keyword_set(Double) or type1 EQ 5 or type1 EQ 9 or type2 EQ 5 or type2 EQ 9
pres0=floatpres(Double = Double)
complex=(type1 EQ 6 or type1 EQ 9 or total(abs(imaginary(Eval)) gt pres0,/int) ne 0)
if complex then begin
    if double then type=9 else type=6
endif else begin
    if double then type=5 else type=4
endelse
AA=A
SetType,AA,type
SetType,Eval,type
epsm=n*MAX(TOTAL(ABS(A), 1, /pres))
pres=pres0*epsm

; ----Round eigenvalues when possible----
for i=0l,n-1 do begin
    ; Characteristic matrix: M=A-λ.I
    M=AA
    if complex then rEval=complex(round(Eval[i]),round(imaginary(Eval[i]))) $
    else rEval=round(Eval[i])
    M[diag]-=rEval
    
    ; Decompose M
    if complex then begin
        LA_SVD, M, W, U, V, status=stat
        if stat ne 0 then continue
    endif else SVDC, M, W, U, V

    ; Null space of A-λI, i.e. solve (A-λI).v = 0
    tmp=abs(max(W,/abs))>1
    indNull=where(abs(W)/tmp le pres,nullity)
    if nullity ne 0 then Eval[i]=rEval
endfor

; ----Check type of the eigenvalues----
type2=size(Eval,/type)
double=keyword_set(Double) or type1 EQ 5 or type1 EQ 9 or type2 EQ 5 or type2 EQ 9
complex=(type1 EQ 6 or type1 EQ 9 or total(abs(imaginary(Eval)) gt pres0,/int) ne 0)
if complex then begin
    if double then type=9 else type=6
endif else begin
    if double then type=5 else type=4
endelse
SetType,A,type
SetType,Eval,type
epsm=n*MAX(TOTAL(ABS(A), 1, /pres))
pres=pres0*epsm

; ----Sanity check: trace(A) = sum(eigenvalues) ----
;error = abs(total(Eval,Double = Double) - total(A[diag],Double = Double)) gt pres ; needed?
;if error then return,0

error=0b
return,Eval
end;function mHQR
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function JordanDecomp, Ain, Double = Double, Jordan=Jordan, error=error, forcereal=forcereal
; A = V.J.V^(-1)
; A: square matrix
; J: Jordan normal form of A (eigenvalues on the diagonal)
; V: similarity transform

A=Ain
s=DimSize(A,2)
error=s[0] ne s[1]
if error then return,0
n=s[0]

CATCH, Error_status
error=Error_status NE 0
IF error THEN return,0 ; HQR failed

; Characteristic polynomial
; -------------------------
; An eigenvector v of a linear transformation A (nxn) is a vector that is only scaled (with eigenvalue λ): 
;         A.v = λ.v
;     <=>    (A-λI).v=0
; This has a solution when there is no inverse (A-λI)^(-1), otherwise v=0 is the only solution: 
;         det(A-λI)=0            (no inverse == singular == det(...)=0)
;     <=> P(λ) = 0            (characteristic polynomial of A: P(λ) = λ^n + b1.λ^(n-1) + ...)
;     <=> (λ-c1).(λ-c2)...=0     (decomposition of the polynomial with roots {c1,...cn})
;     
;     
;     Remark: P(λ) = a.λ^n + b.λ^(n-1) + ... + c is monic if a==1, so the characteristic polynomial is monic
;     
; Eigen decomposition of a square matrix A
; ----------------------------------------
;         A = V.D.V^(-1)
;             D: diagonal matrix with eigenvalues of A on the diagonal (D = λI)
;             V: matrix with eigenvectors of A (column vectors)
;     <=>    D = V^(-1).A.V
;     
;     
; Generalized eigenvectors
; ------------------------
; The algebraic multiplicity of an eigenvalue is the multiplicity of the root λ=ci of the characteristic polynomial
; The geometric multiplicity of an eigenvalue is the dimension of the corresponding eigenspace (i.e. the Null space
;     of the matrix (A-ci.I)), spanned by the corresponding eigenvectors. (geometric multiplicity = nullity of (A-ci.I))
; 1. mg = ma for λ=ci: m-D Null space with m eigenvectors (basis of the Null space)
; 2. mg < ma for λ=ci: mg-D Null space with mg eigenvectors (A is called defective)
;      In the later case we have to find generalized eigenvectors so that mg = ma.
;      Suppose we have a λ with ma=m and mg=1 with one eigenvector v_1 then we have m
;      generalized eigenvectors {v_1,...,v_m} that satisfy:
;         (A-λI).v_k   = v_(k-1)    (v_0=0, v_1=eigenvector)
;     <=>    (A-λI)^k.v_k = 0
; 
; Jordan Decomposition of a square matrix A
; -----------------------------------------
;         A = V.J.V^(-1)
;            J: Jordan Normal form of A (block diagonal matrix)
;             V: matrix with generalized eigenvectors of A (column vectors)
;                i.e. the similarity transform matrix
;     <=>    J = V^(-1).A.V
;
;
; Example: A = [[3,0,1,0,1,0],[1,3,1,-1,0,1],[-1,0,1,0,-1,0],[1,1,1,1,0,1],[0,0,0,0,2,0],[1,0,1,0,0,2]]
;    P(A) = (λ-2)^6 = 0
;    V = [v1,v2,v3,v1',v2',v1'']
;        (A-2I).v3  = v2        <=>    A.v3  = 2.v3 + v2
;        (A-2I).v2  = v1        <=> A.v2  = 2.v2 + v1
;        (A-2I).v1  = 0        <=> A.v1  = 2.v1
;        (A-2I).v2' = v1'    <=> A.v2' = 2.v2' + v1'
;        (A-2I).v1' = 0        <=> A.v1' = 2.v1'
;        (A-2I).v1''= 0        <=> A.v1''= 2.v1''
;        A.V = V.J
;        [A.v1,A.v2,A.v3,A.v1',A.v2',A.v1''] = [2.v1,2.v2 + v1,2.v3 + v2,2.v1',2.v2' + v1',2.v1'']
;    => J = 
;       2       1       0       0       0       0
;       0       2       1       0       0       0
;       0       0       2       0       0       0
;       0       0       0       2       1       0
;       0       0       0       0       2       0
;       0       0       0       0       0       2
;    => V (1's above diagonal) = [v1,v2,v3,v1',v2',v1'']
;       V (1's below diagonal) = [v3,v2,v1,v2',v1',v1'']
;    Jordan chains (cyclic subspaces):
;        v1  -> v2 -> v3
;        v1' -> v2'
;        v1''


; ----Get eigenvalues of A----
Eval=mHQR(A,pres=pres,error=error,Double=Double)

; ----Algebraic multiplicity of eigenvalues----
Eval=Eval[reverse(sort(Eval))] ; decreasing order
malg=replicate(1,n)

indUniq = where(abs(Eval - shift(Eval, -1)) gt pres, nEval)
if nEval eq 0 then begin
  nEval=1
  indUniq=[0l]
  malg=[n]
endif else begin
  indUniq=[0,indUniq[0:nEval-2]+1]
  malg = [indUniq[1:*],n]-indUniq
endelse

Eval=Eval[indUniq]

; ----Find (generalized) eigenvectors for each eigenvalue----
repeat begin
    ; ----Find (generalized) eigenvectors----
    status=JordanDecompHelper(A,Eval,nEval,malg,Evec,Jordan,pres)
    
    ; ----Do we have the Jordan normal form?----
    case status of
    -1:    begin
        ; Jordan normal form can't be found
        error=1b
        return,0
        endcase
    0:    error=0b ; Jordan normal form found
    1:    begin
        ; The geometric multiplicity exceeded the algebraic multiplicity for an eigenvalue or
        ; the Jordan normal form sanity check failed.
        ; This might mean some eigenvalues are actually the same. Merge the two closest eigenvalues
        ; and try again to find the Jordan normal form of A.
        error = nEval eq 1
        if error then return,0
        
        ; Merge two closest eigenvalues
        tmp=min(Eval[1:*]-Eval[0:nEval-2],/abs,i)
        
        Eval[i]=(malg[i]*Eval[i]+malg[i+1]*Eval[i+1])/(malg[i]+malg[i+1])
        malg[i]+=malg[i+1]
        b=bytarr(nEval)
        b[i+1]=1
        ind=where(~b,nEval)
        Eval=Eval[ind]
        malg=malg[ind]
        
        ; Change type from complex to real?
        type=size(Eval,/type)
        if type EQ 6 or type EQ 9 then begin
            pres0=floatpres(Double = Double)
            complex=total(abs(imaginary(Eval)) gt pres0,/int) ne 0 and not keyword_set(forcereal)
            if ~complex then begin
                if double then type=5 else type=4
                SetType,A,type
                SetType,Eval,type
            endif
        endif

        endcase
    endcase
endrep until status eq 0

; ----Return the (generalized) eigenvectors----
return,Evec
end;function JordanDecomp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ShurDecomp, A, Double = Double, Eval=Eval, Shur=S, error=error
; A = V.S.V^T
; A: square matrix
; S: Shur form of A (upper triangular, eigenvalues on the diagonal)
; V: unitary

; M is normal <=> M^T.M = M.M^T   (T: conjugate transpose)
; special cases:
;    M is unitary <=> V.V^T=I            (this is called orthogonal when not complex)
;    M is hermitian <=> V[i,j]=V[j,i]^T  (this is called symmetric when not complex)

s=DimSize(A,2)>1
error=s[0] ne s[1]
if error then return,0

; Hessenberg decomposition:
; A = Q.H.Q^T
;    H: upper Hessenberg form of A
;    Q: unitary
S = LA_ELMHES(A,V,Double = Double,PERMUTE_RESULT=permute,NORM_BALANCE=L1normA,SCALE_RESULT=scale)

; Shur decomposition:
Eval = LA_HQR(S,V,Double = Double,PERMUTE_RESULT=permute,status=status)
error=status ne 0

return,V
end;function ShurDecomp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function COvar, A, ipvt
; Covar(A) = invert(A^(T) ## A)
;          = P ## invert(R^(T) ## R) ## P^(T) ; QR fact.

sz = size(A)
if sz[0] NE 2 then return, -1L
n = sz[1]
if n NE sz[2] then return, -1L

zero = A[0] * 0.
one  = zero  + 1.
if n_elements(ipvt) EQ 0 then ipvt = lindgen(n)
r = A
r = reform(A, n, n, /overwrite)

; Form the inverse of r in the full upper triangle of r
l = -1L
tol = one*1.E-14
tolr = tol * abs(r[0,0])

k=0
while (k ne n)?(abs(r[k,k]) gt tolr):0b do begin
    r[k,k] = one/r[k,k]
    for j = 0L, k-1 do begin
        temp = r[k,k] * r[j,k]
        r[j,k] = zero
        r[0,k] = r[0:j,k] - temp*r[0:j,j]
    endfor
    l = k
    k++
endwhile


; Form the full upper triangle of the inverse of (r transpose)*r
; in the full upper triangle of r
if l GE 0 then $
    for k = 0L, l do begin
      for j = 0L, k-1 do begin
          temp = r[j,k]
          r[0,j] = r[0:j,j] + temp*r[0:j,k]
      endfor
      temp = r[k,k]
      r[0,k] = temp * r[0:k,k]
endfor

; Form the full lower triangle of the covariance matrix
; in the strict lower triangle of r and in wa
wa = replicate(r[0,0], n)
for j = 0L, n-1 do begin
      jj = ipvt[j]
      sing = j GT l
      for i = 0L, j do begin
          if sing then r[i,j] = zero
          ii = ipvt[i]
          if ii GT jj then r[ii,jj] = r[i,j]
          if ii LT jj then r[jj,ii] = r[i,j]
      endfor
      wa[jj] = r[j,j]
endfor

; Symmetrize the covariance matrix in r
for j = 0L, n-1 do begin
    r[0:j,j] = r[j,0:j]
    r[j,j] = wa[j]
endfor

return, r
end;function COvar
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
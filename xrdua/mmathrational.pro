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

function GCD_extended,a,b,x,y
; Check type to prevent endless loops
tmp=size(a,/type)
if tmp eq 4 or tmp eq 5 then a=round(a)
tmp=size(b,/type)
if tmp eq 4 or tmp eq 5 then b=round(b)

; Greatest common divisor of a and b: extended Euclidean algorithm
; which also calculates x and y for which a.x+b.y = GCD(a,b)
; A pair Bézout coefficients (x,y) is not unique (but the minimal in absolute value are found)
r0=abs(a)>abs(b)
r1=abs(a)<abs(b)

x=r0*0
y=x+1
lastx=y
lasty=x

ind=where(r1 ne 0,ct)
while ct ne 0 do begin
    r2= r0[ind] mod r1[ind]
    q = r0[ind]/r1[ind]
    
    r0[ind]=r1[ind]
    r1[ind]=r2
    
    u=x[ind]
    v=y[ind]
    x[ind]=lastx[ind]-q*x[ind]
    y[ind]=lasty[ind]-q*y[ind]
    lastx[ind]=u
    lasty[ind]=v
    
    ind=where(r1 ne 0,ct)
endwhile

; x belongs to a>b and y belongs to a<b
; but it should be
; x belongs to a and y belongs to b
ind=where(abs(a) lt abs(b),ct)
if ct ne 0 then begin
    ; swap x and y
    tmp=lastx[ind]
    lastx[ind]=lasty[ind]
    lasty[ind]=temporary(tmp)
endif 

; All was done for absolute values of a and b
; This has no effect on r0 (except that it is already positive)
; but it has an effect on a and b. Adding sign back to a and b
; leads to adding the same sign to coefficient x and y
x=lastx*(a/abs(a))
y=lasty*(b/abs(b))

; Handle special cases
ind=where(a eq 0 and b eq 0,ct)
if ct ne 0 then begin
    r0[ind]=1
    x[ind]=0
    y[ind]=0
endif

tmp=CHECK_MATH() ; zero's are handled, clear math errors
return,r0 ; always positive
end;function GCD_extended
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GCD,a,b
; Check type to prevent endless loops
tmp=size(a,/type)
if tmp eq 4 or tmp eq 5 then a=round(a)
tmp=size(b,/type)
if tmp eq 4 or tmp eq 5 then b=round(b)

; Greatest common divisor of a and b: Euclidean algorithm
r0=a>b
r1=a<b

ind=where(r1 ne 0,ct)
while ct ne 0 do begin
    r2= r0[ind] mod r1[ind]
    q = r0[ind]/r1[ind]
    
    r0[ind]=r1[ind]
    r1[ind]=r2
    
    ind=where(r1 ne 0,ct)
endwhile

; Handle special cases
ind=where(a eq 0 and b eq 0,ct)
if ct ne 0 then r0[ind]=1

tmp=CHECK_MATH() ; zero's are handled, clear math errors
return,abs(r0)
end;function GCD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GCDmore,arr
; Greatest common divisor of a,b,c,...
; = GCD(a,GCD(b,GCD(c,...)))

if n_elements(arr) eq 2 then return,GCD(arr[0],arr[1])
return,GCD(arr[0],GCDmore(arr[1:*]))
end;function GCDmore
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LCM,a,b
; Least common multiple of a and b
return,abs(a * b) / GCD(a, b)
end;function LCM
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LCMmore,arr
; Least common multiple of a,b,c,...
; = LCM(a,LCM(b,LCM(c,...)))

if n_elements(arr) eq 2 then return,LCM(arr[0],arr[1])
return,LCM(arr[0],LCMmore(arr[1:*]))
end;function LCMmore
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prime_sieve,n
; Sieve of Eratosthenes: find all primes less or equal to n

; Possible primes: 0,1,2,3,...,n
primes=l64indgen(n+1)

; Zero out all not primes
primes[1]=0
for i=2ll,n-1 do begin
    ; Skip zero'ed out prime
    if primes[i] eq 0 then continue
    
    ; Zero out all multiples of primes[i]
    ind=where((primes[i+1:*] mod primes[i]) eq 0,ct,/L64)
    if ct ne 0 then primes[ind+i+1]=0
endfor

; Keep all primes
return,primes[where(primes ne 0,/L64)]
end;function prime_sieve
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IntFact,m
; Integer factorization by trial division

n=abs(m)
if n eq 0 then return,0
if n eq 1 then return,1

primes = prime_sieve(long64(n^0.5) + 1)
prime_factors = 0

for i=0ll,n_elements(primes)-1 do begin
    p=primes[i]
    if p*p gt n then break
    while (n mod p) eq 0 do begin
        prime_factors=[prime_factors,p]
        n /= p
    endwhile
endfor

if n gt 1 then prime_factors=[prime_factors,n]

return, prime_factors[1:*]
end;function IntFact
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isrational,M
if size(M[0],/type) ne 8 then return,0b
if n_tags(M) ne 2 then return,0b
names=tag_names(M)
if names[0] ne 'NUM' and names[1] ne 'DENOM' then return,0b
return,1b
end;function isrational
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function make_ratarray,s,value=value,_EXTRA=re
zero=(make_array(1,_EXTRA=re))[0]
bscalar=n_elements(s) eq 1 and s[0] eq 0
if bscalar then s=1
ret=make_array(s,value={NUM:zero,DENOM:zero+1b})
if bscalar then ret=ret[0]
if n_elements(value) ne 0 then begin
    if isrational(value) then begin
        ret.num=value.num
        ret.denom=value.denom
    endif else ret.num=value
endif
return,ret
end;function make_ratarray
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function identityrat,N_,_EXTRA=ex
if n_elements(N_) ne 1 then begin
    N=n_elements(N_)
    val=N_
endif else begin
    N=N_
    val=1
endelse
Array=make_ratarray([N,N],_EXTRA=ex)
Array[LINDGEN(N) * (N+1)].num = val
return,Array
end;function identityrat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro rational,M,L64=L64
if ~isrational(M) then M=make_ratarray(size(M,/dim),value=M,type=keyword_set(L64)?14:size(M,/type))
end;pro rational
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro printrational,M
if n_tags(M) ne 2 then begin
    print,M
    return
endif

names=tag_names(M)
b=names eq 'NUM'
b or=names eq 'DENOM'
if b[0]+b[1] ne 2 then begin
    print,M
    return
endif

s=size(M,/dim)
out=strarr(s)
ind=where(M.denom ne 1,ct)
if ct ne 0 then out[ind]='/'+stringr(M[ind].denom)
out=stringr(M.num)+out
format='('+stringr(s[0] eq 0?1:s[0])+'A'+stringr(max(strlen(out))+2)+')'
print,out,format=format
end;pro printrational
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReduceNumDenom,M
; Reduce 0/x to 0/1
ind=where(M.NUM eq 0 and M.DENOM ne 0,ct)
if ct ne 0 then M[ind].DENOM=1

; Smallest possible numerator and denominator giving the same rational number
gcd=gcd(M.NUM,M.DENOM)
M.NUM/=gcd
M.DENOM/=gcd
end;pro ReduceNumDenom
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Rational2Float,val,double=double
if keyword_set(double) then return,val.num/double(val.denom)$
else return,val.num/float(val.denom)
end;function Rational2Float
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro Float2Rational,trn,base_,L64=L64,thres=thres
; Convert floating point to a rational with maximal denominator "base"

base=base_
if keyword_set(thres) then begin
    num = make_array(dimsize(trn,2),value=1,type=keyword_set(L64)?14:3)
    denom = num
    trnabs=abs(trn)
    
    ; Threshold = minimal distance between two rationals in [-1,1] with n=maximal denominator
    ;   threshold = (n-1)/n - (n-2)/(n-1)
    ;   threshold = 1/(n.(n-1))
    ; The number might be smaller due to roundoff:
    ;   tmp=721/360. & Float2Rational,tmp,360,/thres & printrational,tmp
    ;   tmp=721/360d & Float2Rational,tmp,360,/thres & printrational,tmp
    base>=2
    threshold=1d/(base*double(base)) ; make the threshold a bit smaller
    
    ; Handle zero's
    bzero = trnabs lt threshold
    indzero = where(bzero,nzero)
    if nzero ne 0 then num[indzero]=0

    ; Select rationals to aren't close enough
    diff = abs(num/double(denom))-trnabs
    if nzero ne 0 then diff[indzero]=0
    indb = where(abs(diff) ge threshold,nb)
    
    while nb ne 0 do begin
        ind=where(diff[indb] lt 0,ct,comp=comp,ncomp=ncomp)
        if ct ne 0 then num[indb[ind]]++
        if ncomp ne 0 then begin
            ind=indb[comp]
            denom[ind]++
            num[ind]=trnabs[ind]*denom[ind]
            
            ; Denominator is too high: add to the list of zero's
            bmax=denom eq (base+1) ; +1 because one might want to add numbers to the numerator
            indmax=where(bmax,nmax)
            if nmax ne 0 then begin
                bzero or= bmax
                indzero = where(bzero,nzero)
            endif
        endif
        
        ; Select rationals to aren't close enough
        diff = abs(num/double(denom))-trnabs
        if nzero ne 0 then diff[indzero]=0
        indb = where(abs(diff) ge threshold,nb)
    endwhile
    
    num *= long(sign(trn))
    
    trn=make_ratarray(dimsize(trn,2),type=keyword_set(L64)?14:3)
    trn.num=temporary(num)
    trn.denom=temporary(denom)
    
endif else begin 
    ; Denominator must divide base
    trn*=base
    trn=round(trn)
    rational,trn,L64=L64
    trn.denom=base
endelse

ReduceNumDenom,trn
end;pro Float2Rational
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RationalModZ,a,m
SignNum,a
ReduceNumDenom,a

if n_elements(m) eq 0 then m=1
mult=a.denom*m
a.num -= (a.num/mult)*mult
a.num += (a.num lt 0)*mult
end;pro RationalModZ
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SignNum,M
; Negative sign in numerator
ind=where(M.DENOM lt 0,ct)
if ct ne 0 then begin
    M[ind].NUM=-M[ind].NUM
    M[ind].DENOM=-M[ind].DENOM
endif
end;pro SignNum
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rationalsign,M
; Sign of a rational number
return,(~(M.num lt 0 xor M.denom lt 0))*2-1
end;function rationalsign
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratarray_equal,M1,M2
if ~array_equal(dimsize(M1),dimsize(M2)) then return,0b
if ~array_equal(rationalsign(M1),rationalsign(M2)) then return,0b
if ~array_equal(abs(M1.num),abs(M2.num)) then return,0b
if ~array_equal(abs(M1.denom),abs(M2.denom)) then return,0b
return,1b
end;function ratarray_equal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratarray_compare,M1_,M2_,operator
M1=M1_
rational,M1
M2=M2_
rational,M2

case STRUPCASE(operator) of
'LT':    begin
        M1=M1.num/double(M1.denom)
        M2=M2.num/double(M2.denom)
        return,M1 lt M2
        endcase
'GT':    begin
        M1=M1.num/double(M1.denom)
        M2=M2.num/double(M2.denom)
        return,M1 gt M2
        endcase
'LE':    begin
        M1=M1.num/double(M1.denom)
        M2=M2.num/double(M2.denom)
        return,M1 le M2
        endcase
'GE':    begin
        M1=M1.num/double(M1.denom)
        M2=M2.num/double(M2.denom)
        return,M1 ge M2
        endcase
'EQ': return, (rationalsign(M1) eq rationalsign(M2)) and (abs(M1.num) eq abs(M2.num)) and (abs(M1.denom) eq abs(M2.denom))
'NE': return, (rationalsign(M1) ne rationalsign(M2)) or (abs(M1.num) ne abs(M2.num)) or (abs(M1.denom) ne abs(M2.denom))
endcase

end;function ratarray_compare
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReduceNumDenomPrime,N,D,sign
; N and D are factorized integers, return smallest N/D

ind=where(N eq 0,ctN)
ind=where(D eq 0,ctD)
if ctN ne 0 and ctD eq 0 then begin ; N==0, D!=0
    N=0
    D=1
endif else $
if ctN eq 0 and ctD ne 0 then begin ; N!=0, D==0
    N=1
    D=0
endif else $
if ctN ne 0 and ctD ne 0 then begin ; N==0, D==0
    N=0
    D=0
endif else begin ; N!=0, D!=0
    ; Remove identical primes
    for k=0,n_elements(N)-1 do begin
        ind=where(D eq N[k],ct)
        if ct ne 0 then begin
            N[k]=0
            D[ind[0]]=0
        endif
    endfor
        
    ind=where(N ne 0,ct)
    if ct ne 0 then N=product(N[ind],/int) else N=1
    
    ind=where(D ne 0,ct)
    if ct ne 0 then D=product(D[ind],/int) else D=1
    
endelse

N*=sign
end;pro ReduceNumDenomPrime
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReduceNumDenomPrime2,NA,NB,D,signA,signB
; NA, NB and D are factorized integers, return smallest (NA+NB)/D

ind=where(NA eq 0,ctNA)
ind=where(NB eq 0,ctNB)
ind=where(D eq 0,ctD)
if (ctNA ne 0 and ctNB ne 0) and ctD eq 0 then begin; NA==0, NB==0, D!=0
    NA=0
    D=1
endif else $
if (ctNA eq 0 xor ctNB eq 0) and ctD ne 0 then begin; NA!=0, NB==0, D==0  or NA==0, NB!=0, D==0
    NA=1
    D=0
endif else $
if (ctNA eq 0 and ctNB eq 0) and ctD ne 0 then begin; NA!=0, NB!=0, D==0
    NA=product(NA,/int)+product(NB,/int)
    if NA ne 0 then NA=1
    D=0
endif else $
if (ctNA ne 0 and ctNB ne 0) and ctD ne 0 then begin; NA==0, NB==0, D==0
    NA=0
    D=0
endif else $
if (ctNA eq 0 and ctNB ne 0) and ctD eq 0 then begin; NA!=0, NB==0, D!=0
    ReduceNumDenomPrime,NA,D,signA
endif else $
if (ctNA ne 0 and ctNB eq 0) and ctD eq 0 then begin; NA==0, NB!=0, D!=0
    NA=NB
    ReduceNumDenomPrime,NA,D,signB
endif else begin; NA!=0, NB!=0, D!=0
    ; Remove identical primes
    for k=0,n_elements(NA)-1 do begin
        ind1=where(NB eq NA[k],ct1)
        ind2=where(D eq NA[k],ct2)
        if ct1 ne 0 and ct2 ne 0 then begin
            NA[k]=0
            NB[ind1[0]]=0
            D[ind2[0]]=0
        endif
    endfor
    
    ind=where(NA ne 0,ct)
    if ct ne 0 then NA=product(NA[ind],/int) else NA=1
    
    ind=where(NB ne 0,ct)
    if ct ne 0 then NB=product(NB[ind],/int) else NB=1
    
    ind=where(D ne 0,ct)
    if ct ne 0 then D=product(D[ind],/int) else D=1
    
    NA=signA*NA+signB*NB

    ; Running ReduceNumDenomPrime,NA,D is pointless, N/D reduction will be done afterwards
endelse

end;pro ReduceNumDenomPrime2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ratinv,M
; Inverse of rational numbers
tmp=M.NUM
M.NUM=M.DENOM
M.DENOM=temporary(tmp)
SignNum,M
end;pro ratinv
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ratneg,M
; Negative of rational numbers
M.NUM=-M.NUM
end;pro ratneg
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratmultfact,M1,M2
; Multiply rational numbers (sign in numerator, reduced N/D)

; Check dimensions
bscalar=0b
s1=DimSize(M1)
s2=DimSize(M2)
if n_elements(M1) eq 1 then begin
    M3=M2
    N1=IntFact(M1.num)
    D1=IntFact(M1.denom)
    sign1=rationalsign(M1)
    bscalar=1b
endif else if n_elements(M2) eq 1 then begin
    M3=M1
    N2=IntFact(M2.num)
    D2=IntFact(M2.denom)
    sign2=rationalsign(M2)
    bscalar=2b
endif else if array_equal(s1,s2) then M3=M1 else stop,'Wrong dimensions.'

; Use prime factorization to avoid (if possible) overflow
for i=0l,n_elements(M3)-1 do begin
    ; Prime factorization to numerator and denominator
    if bscalar ne 1 then begin
        N1=IntFact(M1[i].num)
        D1=IntFact(M1[i].denom)
        sign1=rationalsign(M1[i])
    endif
    if bscalar ne 2 then begin
        N2=IntFact(M2[i].num)
        D2=IntFact(M2[i].denom)
        sign2=rationalsign(M2[i])
    endif
    N3=[N1,N2]
    D3=[D1,D2]
    
    ; Remove identical primes in numerator and denominator
    ReduceNumDenomPrime,N3,D3,sign1*sign2
    M3[i].num=N3
    M3[i].denom=D3
endfor

;ReduceNumDenom,M3 ; Not needed, done by ReduceNumDenomPrime
return,M3
end;function ratmultfact
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratmultfast,M1,M2
; Multiply rational numbers (sign in numerator, reduced N/D)

; Check dimensions
bscalar=0b
s1=DimSize(M1)
s2=DimSize(M2)
if n_elements(M1) eq 1 then begin
    M3=M2
    R1=M1
    bscalar=1b
endif else if n_elements(M2) eq 1 then begin
    M3=M1
    R2=M2
    bscalar=2b
endif else if array_equal(s1,s2) then M3=M1 else stop,'Wrong dimensions.'

; Use prime factorization to avoid (if possible) overflow
for i=0l,n_elements(M3)-1 do begin
    ; Prime factorization to numerator and denominator
    if bscalar ne 1 then R1=M1[i]
    if bscalar ne 2 then R2=M2[i]
    M3[i].num=R1.num*R2.num
    M3[i].denom=R1.denom*R2.denom
endfor

ReduceNumDenom,M3
return,M3
end;function ratmultfast
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratmult,M1,M2,fact=fact
if keyword_set(fact) then return,ratmultfact(M1,M2) $
else return,ratmultfast(M1,M2)
end;function ratmult
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratdiv,M1,M2,fact=fact
; Divide rational numbers
M3=M2
ratinv,M3
return,ratmult(M1,M3,fact=fact)
end;function ratdiv
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratsumfact,M1,M2
; Add rational numbers (sign in numerator, reduced N/D)

; Check dimensions
bscalar=0b
s1=DimSize(M1)
s2=DimSize(M2)
if n_elements(M1) eq 1 then begin
    M3=M2
    N1=IntFact(M1.num)
    D1=IntFact(M1.denom)
    signA=rationalsign(M1)
    bscalar=1b
endif else if n_elements(M2) eq 1 then begin
    M3=M1
    N2=IntFact(M2.num)
    D2=IntFact(M2.denom)
    signB=rationalsign(M2)
    bscalar=2b
endif else if array_equal(s1,s2) then M3=M1 else stop,'Wrong dimensions.'

; Use prime factorization to avoid (if possible) overflow
for i=0l,n_elements(M3)-1 do begin
    ; Prime factorization to numerator and denominator
    if bscalar ne 1 then begin
        N1=IntFact(M1[i].num)
        D1=IntFact(M1[i].denom)
        signA=rationalsign(M1[i])
    endif
    if bscalar ne 2 then begin
        N2=IntFact(M2[i].num)
        D2=IntFact(M2[i].denom)
        signB=rationalsign(M2[i])
    endif
    N3A=[N1,D2]
    N3B=[N2,D1]
    D3=[D1,D2]
    
    ; Remove identical primes in numerator and denominator
    ReduceNumDenomPrime2,N3A,N3B,D3,signA,signB
    M3[i].num=N3A
    M3[i].denom=D3
endfor

ReduceNumDenom,M3 ; needed, not done by ReduceNumDenomPrime2
return,M3

end;function ratsumfact
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratsumfast,M1,M2
; Add rational numbers (sign in numerator, reduced N/D)

; Check dimensions
bscalar=0b
s1=DimSize(M1)
s2=DimSize(M2)
if n_elements(M1) eq 1 then begin
    M3=M2
    R1=M1
    bscalar=1b
endif else if n_elements(M2) eq 1 then begin
    M3=M1
    R2=M2
    bscalar=2b
endif else if array_equal(s1,s2) then M3=M1 else stop,'Wrong dimensions.'

; Use prime factorization to avoid (if possible) overflow
for i=0l,n_elements(M3)-1 do begin
    ; Prime factorization to numerator and denominator
    if bscalar ne 1 then R1=M1[i]
    if bscalar ne 2 then R2=M2[i]
    
    M3[i].num=R1.num*R2.denom+R2.num*R1.denom
    M3[i].denom=R1.denom*R2.denom
endfor

ReduceNumDenom,M3
return,M3

end;function ratsumfast
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratsum,M1,M2,fact=fact
if keyword_set(fact) then return,ratsumfact(M1,M2) $
else return,ratsumfast(M1,M2)
end;function ratsum
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rattotal,array,dim,_extra=ex
g=LCMmore(array.denom)
if n_elements(dim) ne 0 then t=total(array.num*(g/array.denom),dim,_extra=ex) $
else t=total(array.num*(g/array.denom),_extra=ex)
rational,t
t.denom=g
ReduceNumDenom,t
return,t
end;function rattotal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratsub,M1,M2,fact=fact
; Subtract rational numbers
M3=M2
ratneg,M3
return,ratsum(M1,M3,fact=fact)
end;function ratsub
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratrebin,M,snew
sold=DimSize(M,n_elements(snew))
return,M[rebin(lindgen(sold),snew,/sample)]
end;function ratrebin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratmattranspose,A,P
bP=n_elements(P) ne 0
if bP then B=transpose(A.num,P) else B=transpose(A.num)
rational,B
if bP then B.denom=transpose(A.denom,P) else B.denom=transpose(A.denom)
return,B
end;function ratmattranspose
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratmatinvert,A,singular=singular,fact=fact
; Calculate inverse through Gauss–Jordan elimination of [A,I]

singular=1b

; Dimensions
s=DimSize(A,2)
mcol=s[0]
nrow=s[1]

; Square?
if mcol ne nrow then return,0

; Covert to rational array if needed
M=A
rational,M

; Augmented matrix
M=[M,identityrat(mcol,type=size(M.num,/type))]
mcol+=mcol

; Gauss–Jordan elimination of M
j=0l
bcol=bytarr(nrow)
for i=0l,mcol-1 do begin
    ; Find the pivot row
    ma=max(M[i,j:*].NUM/double(M[i,j:*].DENOM),pivot,/abs)
    pivot += j
      
    if abs(ma) eq 0 then begin
        ;Skip column i, making sure the approximately zero terms are
        ;actually zero.
          M[i,j:*].NUM = 0
          M[i,j:*].DENOM = 1
    endif else begin
        ;Swap current row and pivot row (1)
        if j ne pivot then $
            M[i:*,[pivot, j]] = M[i:*,[j, pivot]]

        ;Normalize pivot row (2)
        M[i:*,j] = ratdiv(M[i:*,j],M[i,j],fact=fact)

        ;Eliminate the current column (3)
        bcol[*]=1b
        bcol[j]=0b
        ridx = where(bcol,nj)
        mi=mcol-i

        tmp1 = ratrebin([M[i,ridx]],[mi,nj])
        tmp2 = ratrebin([M[i:*,j]],[mi,nj])
        tmp1 = ratmult(tmp1,tmp2,fact=fact)
        M[i:*,ridx] = ratsub(M[i:*,ridx],tmp1,fact=fact)

        ;Check if done
        if ++j eq nrow then break
    endelse
    
endfor

; Matrix M=[I,A^(-1)], unless A is singular
ind=where(M[0:nrow-1,*].num ne 0,ct)
if ct eq nrow then $
    singular=~array_equal(ind,LINDGEN(nrow) * (nrow+1))

; Keep A^(-1)
M=M[nrow:*,*]

return,M
end;function ratmatinvert
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratmatdet,A
; Calculate the determinant of A using the row echelon form of a matrix

; Dimensions
s=DimSize(A,2)
mcol=s[0]
nrow=s[1]

; Square?
if mcol ne nrow then return,0

; Covert to rational array if needed
M=A
rational,M

; Row echelon form of M
det=M[0]
det.num=1
det.denom=1

j=0l
for i=0l,mcol-1 do begin
    ; Find the pivot row
    ma=max(M[i,j:*].NUM/double(M[i,j:*].DENOM),pivot,/abs)
    pivot += j
      
    if abs(ma) eq 0 then begin
        ;Skip column i, making sure the approximately zero terms are
        ;actually zero.
          M[i,j:*].NUM = 0
          M[i,j:*].DENOM = 1
    endif else begin
        ;Swap current row and pivot row (1)
        if j ne pivot then begin
            M[i:*,[pivot, j]] = M[i:*,[j, pivot]]
            ratneg,det
        endif

        ;Normalize pivot row (2)
        det = ratmult(det,M[i,j],fact=fact)
        M[i:*,j] = ratdiv(M[i:*,j],M[i,j],fact=fact)
    
        ;Eliminate the current column (3)
        nj=nrow-1-j
        if nj ne 0 then begin
            ridx=j+1+lindgen(nj)
            mi=mcol-i
            
            tmp1 = ratrebin([M[i,ridx]],[mi,nj])
            tmp2 = ratrebin([M[i:*,j]],[mi,nj])
            tmp1 = ratmult(tmp1,tmp2,fact=fact)
            M[i:*,ridx] = ratsub(M[i:*,ridx],tmp1,fact=fact)
        
;            for k=0,nj-1 do det = ratmult(det,M[i,ridx[k]],fact=fact)
        endif
    
        ;Check if done
        if ++j eq nrow then break
    endelse
endfor

for k=0,mcol-1 do det = ratmult(det,M[k,k],fact=fact)

return,det
end;function ratmatdet
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ratmatmult,M1,M2,fact=fact

s1=dimsize(M1)
s2=dimsize(M2)
if n_elements(s1) ne 2 or n_elements(s2) ne 2 then stop,"Array's must be two dimensional."
if s1[0] ne s2[1] then stop,"Array's have incompatible dimensions."

M3=make_ratarray([s2[0],s1[1]],type=size(M1[0].num*M2[0].num,/type))

for j=0l,s1[1]-1 do $
    for i=0l,s2[0]-1 do begin
        tmp=ratmult(M1[*,j],reform(M2[i,*]),fact=fact)
        for k=0l,s1[0]-1 do M3[i,j]=ratsum(M3[i,j],tmp[k],fact=fact)
    endfor

return,M3
end;function ratmatmult
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function equalmodn,a,b,n
x=a mod n
ind=where(x lt 0,ct)
if ct ne 0 then x[ind]+=n

y=b mod n
ind=where(y lt 0,ct)
if ct ne 0 then y[ind]+=n

return,array_equal(x,y)
end;function equalmodn
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SolveLinearCongruence,a_,b_,n_,bsol=bsol
; Solve a.x = b (mod n)
;       a.x = b + n.Z   where Z all integers
; a: integer
; b: integer
; n: positive integer > 0
; x: all solutions
;    x = x0 (mod n)
;    x = x0 + n.Z
; x0: particular solution (smallest in absolute value)

x0=0
bsol=0b

a=a_
tmp=size(a,/type)
if tmp eq 4 or tmp eq 5 then a=round(a)
b=b_
tmp=size(b,/type)
if tmp eq 4 or tmp eq 5 then b=round(b)
n=abs(n_)
tmp=size(n,/type)
if tmp eq 4 or tmp eq 5 then n=round(n)
n>=1

; Solvable?
d=GCD_extended(a,n,r,s)
if (b mod d) ne 0 then return,x0 ; no solution

; Find particular solution
x0=r*(b mod n)/d

bsol=1b
return,x0 mod n

end;function SolveLinearCongruence
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function prquot,a,b,rem=r
; Quotient q=a/b with positive remainder r so that q.b+r=a

q=a/b
r=a mod b

negrem=r lt 0
posb=b gt 0

r+=negrem*(2*posb-1)*b
q+=negrem*(2*(~posb)-1)

return,q
end;function prquot
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MinAbsGtZero,array,xy

ind=where(array ne 0,ct)
if ct eq 0 then begin
    if n_params() gt 1 then xy=[-1l,-1]
    return,0
endif

omin=min(array[ind],ind2,/abs)
if n_params() gt 1 then begin
    s=size(array,/dim)
    xy=ind[ind2]
    xy=[xy mod s[0],xy/s[0]]
endif

return,omin
end;function MinAbsGtZero
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function find_non_division,B, P, Q, R, row, col

s=size(B,/dimensions)
mcol=s[0]
nrow=s[1]

for irow=R,nrow-1 do $
    for icol=R,mcol-1 do $
        if (B[icol,irow] mod B[col,row]) ne 0 then begin
            B[col,*]+=B[icol,*]
            Q[col,*]+=Q[icol,*]
            return,1b
        endif


if row ne R then begin
    tmp=B[*,R]
    B[*,R]=B[*,row]
    B[*,row]=tmp

    tmp=P[*,R]
    P[*,R]=P[*,row]
    P[*,row]=tmp
endif


if col ne R then begin
    tmp=B[R,*]
    B[R,*]=B[col,*]
    B[col,*]=tmp

    tmp=Q[R,*]
    Q[R,*]=Q[col,*]
    Q[col,*]=tmp
endif

return,0b
end;function find_non_division
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SmithNormalForm,A,P=P,Q=Q

; Frank Lübeck: On the Computation of Elementary Divisors of Integer Matrices
;
; SVD for integer matrices:
; A: mxn integers
; Find matrices P(mxm) in GL(m,Z) and Q(nxn) in GL(n,Z) so the P.A.Q = B with B(mxn) diagonal
;
; Diagonal matrix B: 
;     - Diagonal elements B[i,i] are called elementary divisors
;     - B[i,i] is a positive integer
;     - The last elementary divisors may be zero: B[j:*,j:*]=0
;     - B[i,i] is always a divider of B[i+1,i+1] unless B[i+1] is zero (i.e. i+1 >= j)
;    - When A is invertible then the biggest elementary divisor is given by
;        Q=ratmatinvert(A)
;        Bmax=LCMmore(Q.DENOM)
;

s=size(A,/dimensions)
mcol=s[0]
nrow=s[1]

Q=Identity(mcol,type=size(A,/type))
P=Identity(nrow,type=size(A,/type))
B = long(A)

row=0 ; row index
col=0 ; column index
R = 0 ; submatrix offset

repeat begin

    ; Find first, smallest, non-zero, absolute value in submatrix
    if (R eq mcol) or (R eq nrow) then return,B
    if MinAbsGtZero(B[R:*,R:*],xy) eq 0 then return,B
    col=xy[0]+R
    row=xy[1]+R

    ; If this element is less then zero: multiply this row in B and P by (-1)
    if B[col,row] lt 0 then begin
        B[*,row]*=(-1)
        P[*,row]*=(-1)
    endif

    min_irow = row
    min_icol = col

    bstop_inner=0b
    repeat begin
        min = 0

        quot = prquot(B[col,R:*],B[col,row])
        for irow=R,nrow-1 do $
            if(irow ne row) then begin
                B[*,irow]-=quot[irow-R]*B[*,row]
                P[*,irow]-=quot[irow-R]*P[*,row]

                if (B[col,irow] ne 0) and (min eq 0 or abs(B[col,irow]) lt min) then begin
                    min = abs(B(col,irow))
                    min_irow = irow
                endif
            endif

        if min eq 0 then begin
            quot=prquot(B[R:*,row],B[col,row],rem=rem)
            for icol=R,mcol-1 do $
                if icol ne col then begin
                    B[icol,row] = rem[icol-R]
                    Q[icol,*]-=quot[icol-R]*Q[col,*]

                    if (B[icol,row] ne 0) and (min eq 0 or abs(B[icol,row]) lt min) then begin
                        min = abs(B[icol,row])
                        min_icol = icol
                    endif
                endif
        endif

        bstop_inner or= min eq 0
        if ~bstop_inner then begin
            row = min_irow
            col = min_icol
        endif

    endrep until bstop_inner

    R+= ~find_non_division(B, P, Q, R, row, col)

endrep until 0b

end;function SmithNormalForm
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro printpol,p
n=n_elements(p)

; absolute value of rational numbers
out=strarr(n)
ind=where(p.denom ne 1,ct)
if ct ne 0 then out[ind]='/'+stringr(abs(p[ind].denom))
out=stringr(abs(p.num))+out

; remove '1'
if n gt 1 then begin
    ind=where(out[1:*] eq '1',ct)
    if ct ne 0 then out[ind+1]=''
endif

; sign of rational numbers
sign=replicate(' + ',n)
ind=where(rationalsign(p) eq -1,ct)
if ct ne 0 then sign[ind]=' - '
sign[n-1]=''

; polynomial exponent
exp='x^'+stringr(lindgen(n))
exp[0]=''
if n gt 1 then exp[1]='x'

; dot
ind=where(exp ne '' and out ne '',ct)
if ct ne 0 then exp[ind]='.'+exp[ind]

; Print
out=sign+out+exp
ind=where(p.num ne 0 or ~finite(p.denom),ct)
if ct eq 0 then print,0 else print,reverse(out[ind]),format='('+stringr(ct)+'A)'
end;pro printpol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PolDegree,p
; p[0] + p[1].x + p[2].x^2 + ... p[r].x^r
; degree r is given by the index of the last non-zero coefficient
ind=where(p.num ne 0 or ~finite(p.denom),ct)
if ct eq 0 then return,0
return,ind[ct-1]
end;function PolDegree
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReducePol,p
p=p[0:PolDegree(p)]
end;pro ReducePol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckZeroPol,p
if PolDegree(p) gt 0 then return,0b
if p.num eq 0 and finite(p.denom) then return,1b
return,0b
end;function CheckZeroPol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function polsum,p1,p2

np1=n_elements(p1)
np2=n_elements(p2)
if np1 gt np2 then $
    p3=ratsum(p1,[p2,make_ratarray(np1-np2,type=size(p2[0].num,/type))])
if np1 lt np2 then $
    p3=ratsum(p2,[p1,make_ratarray(np2-np1,type=size(p1[0].num,/type))])
if np1 eq np2 then $
    p3=ratsum(p1,p2)

ReducePol,p3
return,p3
end;function polsum
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function polsub,p1,p2
p3=p2
ratneg,p3
return,polsum(p1,p3)
end;function polsub
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function polmult,p1,p2
; (p1[0] + p1[1].x + p1[2].x^2 + ... p1[r1].x^r1) * (p2[0] + p2[1].x + p2[2].x^2 + ... p2[r2].x^r2)

np1=n_elements(p1)
np2=n_elements(p2)

p3=make_ratarray(np1+np2-1,type=size(p1[0].num*p2[0].num,/type))
for i=0l,np1-1 do $
    for j=0l,np2-1 do $
        p3[i+j]=ratsum(p3[i+j],ratmult(p1[i],p2[j]))

return,p3
end;function polmult
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function poldiv,p1,p2,rem
; (p1[0] + p1[1].x + p1[2].x^2 + ... p1[r1].x^r1) / (p2[0] + p2[1].x + p2[2].x^2 + ... p2[r2].x^r2)
; 
; p1 = quotient*p2 + rem
; 
; Polynomial long division

quotient=p1[0]
quotient.num=0
quotient.denom=1
rem=p1

degrem=PolDegree(rem)
degp2=PolDegree(p2)
while degrem ge degp2 do begin
    ; Divide coefficients of the highest power
    deg=degrem-degp2
    a=ratdiv(rem[degrem],p2[degp2])
    addq=make_ratarray(deg+1,type=size(a.num,/type))
    addq[deg]=a
    
    ; Add to quotient
    quotient=polsum(quotient,addq)
    
    ; Subtract from remainder
    rem=polsub(rem,polmult(addq,p2))
    degrem=PolDegree(rem)
    
    if CheckZeroPol(rem) then break
endwhile

return,quotient
end;function poldiv
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GCDPol,p1,p2
; Greatest common divisor of two polynomials, Euclidean algorithm

if PolDegree(p2) gt PolDegree(p1) then begin
    pp0=p2
    pp1=p1
endif else begin
    pp0=p1
    pp1=p2
endelse

repeat begin

    q=poldiv(pp0,pp1,pp2)
    
    bstop=CheckZeroPol(pp2)
    if ~bstop then begin
        pp0=pp1
        pp1=pp2
    endif
    
endrep until bstop

; Make it monic to make it unique
r=PolDegree(pp1)
for i=0l,r do pp1[i]=ratdiv(pp1[i],pp1[r])

return,pp1
end;function GCDPol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PolFact,pol

; Gauss's lemma: pol = poldeg0 * polint
poldeg0=pol[0]
poldeg0.num=LCMmore(pol.denom)
poldeg0.denom=GCDmore(pol.num)

polint=polmult(pol,poldeg0)
ratinv,poldeg0

; Factorize integer polynomial polint
stop

return,{field0:poldeg0,field1:polint}
end;function PolFact
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ESTPermute,A,P,Pinv,permut,i,j
; Ozello: P i,j
; 
; Elementary similarity transformation of matrix A -> P^(-1).A.P
; Switch rowi and rowj, switch coli and colj

; Similarity transformation (i.e. A=Pinv##A##P)
A[*,[i,j]]=A[*,[j,i]]
A[[i,j],*]=A[[j,i],*]

; permut stores the permutations of the columns of P
if n_elements(permut) ne 0 then permut[[i,j]]=permut[[j,i]]

if n_elements(P) eq 0 then return

; Elementary matrix P and its inverse
Pinv[*,[i,j]]=Pinv[*,[j,i]]
P[[i,j],*]=P[[j,i],*]

end;pro ESTPermute
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ESTMultiply,A,P,Pinv,i,m
; Ozello: M i,a
; 
; Elementary similarity transformation of matrix A -> P^(-1).A.P
; rowi*=m, coli/=m

; Similarity transformation (i.e. A=Pinv##A##P)
A[*,i]=ratmult(A[*,i],m)
A[i,*]=ratdiv(A[i,*],m)

if n_elements(P) eq 0 then return

; Elementary matrix P and its inverse
Pinv[*,i]=ratmult(Pinv[*,i],m)
P[i,*]=ratdiv(P[i,*],m)

end;pro ESTMultiply
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ESTAddMultiple,A,P,Pinv,i,m,j
; Ozello: L i,a,j
; 
; Elementary similarity transformation of matrix A -> P^(-1).A.P
; rowi+=m*rowj, colj-=m*coli

; Similarity transformation (i.e. A=Pinv##A##P)
A[*,i]=ratsum(A[*,i],ratmult(A[*,j],m))
A[j,*]=ratsub(A[j,*],ratmult(A[i,*],m))

if n_elements(P) eq 0 then return

; Elementary matrix P and its inverse
Pinv[*,i]=ratsum(Pinv[*,i],ratmult(Pinv[*,j],m))
P[j,*]=ratsub(P[j,*],ratmult(P[i,*],m))

end;pro ESTAddMultiple
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ESTAddMultiple2,A,P,Pinv,i,m,j
; Ozello: C i,a,j
; 
; Elementary similarity transformation of matrix A -> P^(-1).A.P
; coli+=m*colj, rowj-=m*rowi

; Same as ESTAddMultiple,A,P,Pinv,j,-m,i

; Similarity transformation (i.e. A=Pinv##A##P)
A[*,j]=ratsub(A[*,j],ratmult(A[*,i],m))
A[i,*]=ratsum(A[i,*],ratmult(A[j,*],m))

if n_elements(P) eq 0 then return

; Elementary matrix P and its inverse
Pinv[*,j]=ratsub(Pinv[*,j],ratmult(Pinv[*,i],m))
P[i,*]=ratsum(P[i,*],ratmult(P[j,*],m))

end;pro ESTAddMultiple2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SearchNonZeroRow,M,j
; Ozello: chercheli
; 
; search in column j from row j+1 till the end
s=DimSize(M,2)
n=s[0]
k=j+1
while (k lt n)?(M[j,k].num eq 0):0b do k++
return,k
end;function SearchNonZeroRow
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SearchNonZeroColumn,M,i,j
; Ozello: chercheco
; 
; search in row i from column j till the end
s=DimSize(M,2)
n=s[0]
k=j
while (k lt n)?(M[k,i].num eq 0):0b do k++
return,k
end;function SearchNonZeroColumn
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MinimalPolbasevec,M,P,Pinv,permut,joff
; Ozello: polymine1
;
; The columns of P are the coordinates of the current basis vectors (e1,...,en)
; with respect to the canonical basis of K^n.
; 
; MinimalPolbasevec transforms matrix M (joff points to the first column/row of M)
; to matrix M' = P^(-1).M.P (d points to the last column/row of C):
; 
; M' = C B1    (suppose joff=0)
;      0 B2
;   
; C =  0 0 0 0 ... 0 -a[r]
;      1 0 0 0 ... 0 -a[r-1]
;      0 1 0 0 ... 0 -a[r-2]
;      0 0 1 0 ... 0 -a[r-3]
;      ...
;      0 0 0 0 ... 1 -a[1]
;      
;     p(x) = x^r + a[1].x^(r-1) + ... + a[r-1].x + a[r]
;
; The polynomial associated with companion matrix C has degree r=d+1.
; Due to the nature of the change of basis P, we can write that
;          f_1 = e_1
;          f_2 = M.f_1 = M.e_1
;          f_3 = M.f_2 = M^2.e_1
;          ...
;          f_r = M.f_(r-1) = M^(r-1).e_1
;          
;          M.f_r = - a[r].f_1 -  a[r-1].f_2 - ... - a[1].f_r
;  
; <=>    M^r.e_1 + a[1].M^(r-1).e_1 + ... + a[r-1].M.e_1 + a[r].e_1 = 0
; <=>    p(M)(e_1) = 0
;
; The polynomial associated to C is the minimal polynomial of e_1
; or in general e_(joff+1) with respect to M. Or in other words, p is
; the minimal polynomial of the joff'th column of P, i.e. basis vector e_(joff+1).

s=DimSize(M,2)
n=s[0]
d=joff

; Search first non-zero element in column d and row > d
k=SearchNonZeroRow(M,d)

; Loop until row-pointer k overflows
while k lt n do begin
    
    ; M[k,k] to M[d+1,d+1]
    if k ne d+1 then ESTPermute,M,P,Pinv,permut,k,d+1
    
    ; Set M[d,d+1]=1
    a=M[d,d+1]
    ratinv,a
    ESTMultiply,M,P,Pinv,d+1,a
    
    ; M[d,i]=0 except for M[d,d+1]
    for i=0l,n-1 do $
        if i ne d+1 then begin
            a=M[d,i]
            ratneg,a
            ESTAddMultiple,M,P,Pinv,i,a,d+1
        endif
    
    ; M has changed depending on column index d:
    ;  when d=0 =>  column d in M  becomes [0,1,0,0,...]
    ;  or   d=1 =>  column d in M  becomes [0,0,1,0,...]
    ;  or   ...
    
    ; Jump to the next column
    d++
    
    ; Search first non-zero element in column d and row > d
    k=SearchNonZeroRow(M,d)
endwhile

return,d

end;function MinimalPolbasevec
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CompanionPol,M,k,j
; Ozello: polyassocbloc
;
; Get associated (monic) polynomial from the Companion matrix C=M[k:j,k:j]
; 
; C=   0 0 0 0 ... 0 -a[r]
;      1 0 0 0 ... 0 -a[r-1]
;      0 1 0 0 ... 0 -a[r-2]
;      0 0 1 0 ... 0 -a[r-3]
;      ...
;      0 0 0 0 ... 1 -a[1]
;      
; p(x) = x^r + a[1].x^(r-1) + ... + a[r-1].x + a[r]

one=M[0]
one.num=1
one.denom=1

coeff=reform(M[j,k:j])
ratneg,coeff

return,[coeff,one] ; a[r],a[r-1],...,a[1],1
end;function CompanionPol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MinimalPolvec,M,P,Pinv,permut,v,joff
; Ozello: polymine
;
; Minimal polynomial of v=v[0].f_1 + ... + v[n-1].f_n (integer coeff)
; with respect to M, where f_i are the basis vectors of R^n
; given with respect to the canonical basis as columns of P.

; Suppose v[r-1] is the first non-zero coefficient of v.

; Move f_r to the column joff of P 
; (r=ind[0])
ind=where(v ne 0,ct)
if ind[0] ne 0 then ESTPermute,M,P,Pinv,permut,joff,ind[0]

; Multiply f_r with v[r-1]
a=v[ind[0]]
rational,a
ratinv,a
ESTMultiply,M,P,Pinv,joff,a

; The joff'th column of P contains v[r-1].f_r
; Now add all other non-zero v[k].f_k
for k=1,ct-1 do begin
    a=v[ind[k]]
    rational,a
    ESTAddMultiple2,M,P,Pinv,joff,a,ind[k]
endfor

; Minimal polynomial of v (now in the joff'th column of P) with respect to M
d=MinimalPolbasevec(M,P,Pinv,permut,joff)
return,CompanionPol(M,joff,d) ; p(M)(v)=0

end;function MinimalPolvec
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MinimalPolSetvec,B,P,Pinv,permutB,joff,tbldim,famgen
; Ozello: extgen
;
; Input: permutB = lindgen(n)

s=DimSize(B,2)
n=s[0]

; Minimal polynomial of base vector f_(joff+1), make new companion matrix
; (Base vector f_(k+1) is given by the kth column of P)
k=joff
d=MinimalPolbasevec(B,P,Pinv,permutB,k)
;printpol,companionpol(B,k,d)
k=d+1 ; column/row after new companion matrix
tbldim=k ; dimension of the subspace spanned by (f_(l0+1),f_(l1+1),...,f_(joff+1))
famgen=joff

while k lt n do begin
    ; Polynomial of base vector f_(k+1), make new companion matrix
    d=MinimalPolbasevec(B,P,Pinv,permutB,k)
    ;printpol,companionpol(B,k,d)
    k=d+1 ; column/row after new companion matrix
    tbldim=[tbldim,k] ; add dimension of the subspace spanned by (f_(l0+1),f_(l0+1),...,f_(k0+1),f_(k1+1),...)
endwhile

if n_elements(tbldim) gt 1 then $
    famgen=[famgen,permutB[tbldim[0:n_elements(tbldim)-2]]]

; The base vectors f_(famgen+1) span an subspace, invariant under B
; The dimension of the subspace spanned by (f_(l0+1),f_(l1+1),...,f_(k0+1)) is given by tbldim[0]
; The dimension of the subspace spanned by (f_(l0+1),f_(l1+1),...,f_(k0+1),f_(k1+1)) is given by tbldim[1]
; etc...

end;pro MinimalPolSetvec
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function canonicalBasisvector,n,k,type=type
; Canonical basis vector k
v=make_array(n,type=type)
v[k]=1
return,v
end;function canonicalBasisvector
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CheckMinimalPolVec,A,ppol,v_,basevec=basevec
; Check whether ppol is the minimal polynomial of v with respect to A.
; 
; ppol[0].M^r.v + ppol[1].M^(r-1).v + ... + ppol[r-1].v + ppol[r].v = 0
;
; This is not complete because it doesn't prove there exists not other polynomial
; with lower degree for which this goes.

s=DimSize(A,2)
n=s[0]
type=size(A[0].num,/type)

if keyword_set(basevec) then v=canonicalBasisvector(n,v_,type=type) else v=v_

rational,v
r=n_elements(ppol)-1
v=ratmult(v,ppol[r])
for i=r-1,0,-1 do v=ratmult(ratsum(v,reform(ratmatmult(A,reform(v,1,n)),n)),ppol[i])

ind=where(v.num ne 0,ct)
if ct ne 0 then stop,'Vector v does not have minimal polynomial ppol!'

end;pro CheckMinimalPolVec
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MinimalPol,A,joff,fam=fam,poly=poly,B=B,permutB=permutB
; Ozello: polymina
; 
; Minimal polynomial of matrix A[joff:*,joff:*]
;
; Input: permutB = lindgen(n)

s=DimSize(A,2)
n=s[0]
type=size(A[0].num,/type)

B=A
MinimalPolSetvec,B,P,Pinv,permutB,joff,tbldim,famgen ; P,Pinv undefined because not used
ppol = CompanionPol(B,joff,tbldim[0]-1) ; ppol = minimal polynomial of base vector f_(joff+1)

fam=famgen[0]
poly={field0:ppol}
indpolk=n_elements(famgen)-1

pol=ppol
dimsubspace=n

;CheckMinimalPolVec,A,ppol,fam,/basevec

while indpolk gt 0 and PolDegree(ppol) lt n and PolDegree(pol) lt dimsubspace do begin
    C=A
    
    k=famgen[indpolk] ; no permutations in C
    v=canonicalBasisvector(n,k,type=type)
    pol=MinimalPolvec(C,P,Pinv,permut,v,joff) ; pol = minimal polynomial of base vector f_(k+1)
    
    g=GCDPol(pol,ppol)
    if PolDegree(g) lt PolDegree(pol) then begin
        ; Some factors of pol are not in ppol
        ppol=polmult(ppol,poldiv(pol,g)) ; multiplication of the minimal polynomials of several base vectors
        fam=[fam,k]
        poly=create_struct(poly,'field'+stringr(n_tags(poly)),ppol)
    endif
    
    dimsubspace=tbldim[indpolk]
    indpolk--
endwhile

; the polynomial ppol is the LCM of basis vectors "fam"

return,ppol

end;function MinimalPol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecMinimalPol2,A,v,w,joff,pLCM,B=B,permutB=permutB
; Ozello: veminimax2
; 
; pLCM(A) is the minimal polynomial given by the least common multiple of
; the minimal polynomials of vector v and w
; 
;    pLCM(A) = LCM(p(A)(v),p(A)(w))
;
; Return a vector u=v+k.w so then p(A)(v+k.w)=pLCM(A)

k=0l
repeat begin
    k++
    B=A
    u=v+k*w
    
    pol=MinimalPolvec(B,P,Pinv,permutB,u,joff)
endrep until ratarray_equal(pol,pLCM)

return,u
end;function vecMinimalPol2
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecMinimalPol,A,fam,poly,imax,joff,B=B,permutB=permutB
; Ozello: veminimax
;
; fam: base vectors f_(fam+1)
; poly.(i) = LCM(minpol(f_(fam[0:i]+1)))
; 
; Find v = lincom(f_(fam[0:imax]+1)) so that minpol(v) = poly.(imax)

s=DimSize(A,2)
n=s[0]
type=size(A[0].num,/type)

; Return base vector fam[0]
if imax eq 0 then return,canonicalBasisvector(n,fam[0],type=type)

; Return v'=v+k.w where v=fam[0] and w=fam[1], so that pv'==poly.(1)
if imax eq 1 then return,vecMinimalPol2(A,canonicalBasisvector(n,fam[0],type=type),$
                                    canonicalBasisvector(n,fam[1],type=type),joff,poly.(1),B=B,permutB=permutB)

; Trivial linear combination of vectors fam[0:imax]
v=canonicalBasisvector(n,fam[0],type=type)
for i=1l,imax do v+=canonicalBasisvector(n,fam[i],type=type)

; Minimal polynomial of this vector
B=A
pv = MinimalPolvec(B,P,Pinv,permutB,v,joff)

; Minimal polynomial equal to poly.(imax)?
if ratarray_equal(pv,poly.(imax)) then return,v

; Linear combnation v of fam[0:imax-1] so that pv==poly.(imax-1)
v=vecMinimalPol(A,fam,poly,imax-1,joff,B=B,permutB=permutB)

; Return v'=v+k.w where w is base vector fam[imax], so that pv'==poly.(imax)
return,vecMinimalPol2(A,v,canonicalBasisvector(n,fam[imax],type=type),joff,poly.(imax),B=B,permutB=permutB)

end;function vecMinimalPol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro FormInterm,A,F,P,blockstart

s=DimSize(A,2)
n=s[0]

; Prepare permutation tracker
permut=lindgen(n)

d=0
blockstart=d
while d lt n do begin
    ; Column/row joff point to the beginning of the submatrix of the intermediate
    ; matrix F which doesn't have diagonal companion matrices yet
    joff=d
    
    ; Minimal polynomial of A[joff:*,joff:*]
    permutB=lindgen(n)
    pp=MinimalPol(F,joff,fam=fam,poly=poly,B=B,permutB=permutB)
    
    ; Set of base vectors f_(fam+1) invariant under A and LCM of their minimal polynomials
    ; gives the minimal polynomial pp of A. Now find a linear combination of the f_(fam+1)
    ; vectors which has minimal polynomial pp.
    imax=n_elements(fam)-1
    v=vecMinimalPol(F,fam,poly,imax,joff,B=B,permutB=permutB)
    
    ; Add new companion matrix dimension
    d+=PolDegree(pp) ; pp = poly.(imax)
    blockstart=[blockstart,d]
    
    ; Update similarity transform P (only integers)
    ; P[joff:d-1,permut] = [v, A.v, A^2.v, ..., A^(d-1).v]
    ;P[joff,permut]=transpose(v)
    ;for k=joff+1,d-1 do P[k,*] += A##P[k-1,*]
    P[joff,permut].num=transpose(v)
    P[joff,permut].denom=1
    for k=joff+1,d-1 do P[k,*] = ratsum(P[k,*],ratmatmult(A,P[k-1,*]))
    
    ; Update F with new companion matrix and new permutations
    F=temporary(B)
    permut=permut[permutB]
endwhile

end;pro FormInterm
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TryRightZeroNewCompanion,M,P,Pinv,i1,i2
; Ozello: zeroadroite
;
; i1 column/row just before the beginning of this companion matrix
; i2 last column/row of this companion matrix

s=DimSize(M,2)
n=s[0]

for i=i2,i1+2,-1 do $
    for j=i2+1,n-1 do begin
        a=M[j,i]
        ratneg,a
        ESTAddMultiple2,M,P,Pinv,j,a,i-1
    endfor

end;pro TryRightZeroNewCompanion
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro RightZeroFormInterm,M,blockstart
; Ozello: zeroadroites

k=n_elements(blockstart)
ind=blockstart-1 ; last row/column of the the different blocks
for h=k-2,1,-1 do $
    TryRightZeroNewCompanion,M,P,Pinv,ind[h-1],ind[h] ; set elements next to block h to zero
    
end;pro RightZeroFormInterm
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NewCompanion,M,P,Pinv,i1
; Ozello: unbloc
; 
; i1 point to last row and column of the last companion matrix C.
; For example suppose M is of the form (B is just any matrix)
;  
;  C1  0 0
;   0 C2 0
;   0  0 B
;
; then i1 points to the last row and column of C2


s=DimSize(M,2)
n=s[0]

repeat begin
    i2=MinimalPolbasevec(M,P,Pinv,permut,i1+1)
    
    ; Matrix M:
    ;  C1  0  0   0
    ;   0 C2  0   0
    ;   0  0 C3' B1
    ;   0  0  0  B2
    ; 
    ; i1 points to last row/column of C2
    ; i2 points to last row/column of C3'
    ;
    ; rows of B1 from i1+1 till i2
    ; columns of B1 from i2+1 till n-1
    
    TryRightZeroNewCompanion,M,P,Pinv,i1,i2
    
    ; Matrix M (D has only zero's, except for the first row):
    ;  C1  0  0   0
    ;   0 C2  0   0
    ;   0  0 C3'  D
    ;   0  0  0  B3'
    
    ; Search first non-zero element in row i1+1 and column > i2,
    ; i.e. first non-zero elements in first row of D
    j=SearchNonZeroColumn(M,i1+1,i2+1)
    
    ; Move this elements to the first element of C3'
    if j lt n then ESTPermute,M,P,Pinv,permut,i1+1,j
    
    ; Matrix M:
    ;  C1  0  0
    ;   0 C2  0
    ;   0  0  B'
    
endrep until j ge n ; D has only zero elements

; Matrix M:
;  C1  0  0  0
;   0 C2  0  0
;   0  0 C3  0
;   0  0  0 B3
    
return,i2 ; points to last row/column of C3 (new companion matrix)

end;function NewCompanion
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DetectCompanions,F,i0,i1
; F: Frobenius normal form

s=DimSize(F)
n=s[0]

; Loop over subdiagonal
dimC=1
Ci=0
for i=0,n-2 do begin
    idiag = i * (n+1)
    isubdiag = idiag+n
    
    bsubdiag0 = F[isubdiag].num eq 0
    bsubdiag1 = F[isubdiag].num eq 1 and F[isubdiag].denom eq 1
    
    if ~bsubdiag0 and ~bsubdiag1 then return,0b
    
    if bsubdiag0 then begin
        Ci++
        dimC=[dimC,1]
    endif else dimC[Ci]++
    
endfor
nC=n_elements(dimC)

; Begin and end index
if nC eq 1 then begin
    i0=0
    i1=dimC-1
endif else begin
    i0=[0,total(dimC[0:n_elements(dimC)-2],/cumul,/int)]
    i1=total(dimC,/cumul,/int)-1
endelse

return,dimC
end;function DetectCompanions
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function IsFrobenius,F,weak=weak

; Companion matrices
DimC=DetectCompanions(F,i0,i1)
nC=n_elements(dimC)

; Everything next to the companion matrices must be zero
F2=F
for i=0,nC-1 do F2[i0[i]:i1[i],i0[i]:i1[i]].num=0
if total(F2.num ne 0,/int) ne 0 then return,0b

; Check Companion matrices
for i=0,nC-1 do begin
    Ci=F[i0[i]:i1[i],i0[i]:i1[i]]
    
    n=dimC[i]
    case n of
    1: 
    2: if Ci[0].num ne 0 then return,0b
    else: begin
        if total(Ci[lindgen(n-1) * (n+1)+n].num ne 1,/int) ne 0 then return,0b
        Ci[lindgen(n-1) * (n+1)+n].num=0
        Ci[n-1,*].num=0
        if total(Ci.num ne 0,/int) ne 0 then return,0b
        endelse
    endcase
endfor

if keyword_set(weak) then return,1b

; Check polynomials
for i=0,nC-2 do begin
    pol=CompanionPol(F,i0[i],i1[i])
    polnext=CompanionPol(F,i0[i+1],i1[i+1])
    
    ; polnext|pol
    if Poldegree(polnext) gt Poldegree(pol) then return,0b
    
    q=poldiv(pol,polnext,rem)
    if n_elements(rem) ne 1 then return,0b
    if rem[0].num ne 0 then return,0b
endfor

return,1b

end;function IsFrobenius
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FrobeniusDecomp, A, P=P, invP=Pinv, weak=weak, error=error
; Source:
;    
; P. Ozello: Calcul exact des formes de Jordan et de Frobenius d'une matrice, second algorithm.
; Thèse de l'Universite Scientifique Technologique et Medicale de Grenoble, 1987
; => written in Matlisp for Common lisp
;
; --------------------------------------
; Endomorphism:
; 
; f is a linear map from vector space V over field K to itself (i.e. an endomorphism)
;     f: V->V
; Choose a basis (e_1,...,e_n) and this can be written in coordinate space
;     f: K^n->K^n: A.X=Y
;     
; FrobeniusDecomp transforms A like this
;         F = P^(-1).A.P
; which is the same as changing the basis with change of basis matrix P
; where then columns are the coordinates new basis vectors with respect to the old ones.
; The endomorphism for the new basis (f_1,...,f_n) can then be written as
;     f: K^n->K^n: F'.X' = Y'
;             <=>  A.P.X' = P.Y'
; Indeed X = P.X' for the coordinates of a vector under a change of basis.
;
; --------------------------------------
; Minimal polynomial:
; 
; Consider a polynomial p with coefficients in field K. Then there is an endomorphism p(f)
;     p(f): V->V
; Choose a basis (e_1,...,e_n) and this can be written in coordinate space
;     p(A): V->V: p(A)(X) = Y
; The polynomial p is a minimal polynomial of f (or A in coordinate space) iff
;     1. Im(p(A)) = 0
;     2. p is monic
;     3. if p'(A) = 0 then the degree of p' is less or equal to p
; For each f (or A in coordinate space) there is precisly one minimal polynomial.
; When p'(A) = 0 then p divides p'.
; 
; --------------------------------------
; Minimal polynomial of a vector relative to an endomorphism:
; 
; The polynomial p is the minimal polynomial of a vector, with respect to A iff
;     1. p(A)(X) = 0
;     2. p is monic
;     3. if p'(A)(X) = 0 then the degree of p' is less or equal to p
; 
; --------------------------------------
; 
; A = P.F.P^(-1)
; A: square matrix (nxn) with rational elements
; F: Frobenius normal form of A
; P: similarity transform
;
; F=diag(C1,C2,...,Cm)
; 
; Ci=  0 0 0 0 ... 0 -c[r]
;      1 0 0 0 ... 0 -c[r-1]
;      0 1 0 0 ... 0 -c[r-2]
;      0 0 1 0 ... 0 -c[r-3]
;      ...
;      0 0 0 0 ... 1 -c[1]
; pi(x) = x^r + c[1].x^(r-1) + ... + c[r-1].x + c[r]
;
; pi(x) is a monic polynomial of degree r.
; Ci is the companion matrix of pi(x).
; 
; 1. pi is the characteristic polynomial of Ci, i.e. pi(x)=|Ci-x|
; 2. pi is the minimal polynomial of Ci:
;         pi(Ci)(X) = Ci^r.X + c[1].Ci^(r-1).X + ... + c[r-1].Ci.X + c[r].X = 0
; 3. Polynomial pi+1 divides pi
; 4. p1 is the minimal polynomial of A
; 5. p1.p2...pm is the characteristic polynomial of A
;
; --------------------------------------
; Calculate characteristic polynomial
; wims.unice.fr/wims/en_tool~linear~matrix.html
;

; Check square matrix
s=DimSize(A,2)
error=s[0] ne s[1]
if error then return,0
n=s[0]

; Initialize frobenius normal form
F=A ; rational or integer
rational,F,/l64
    
; Convert to rational array and initialize similarity transform
delvar2,P & delvar2,Pinv
if keyword_set(weak) then begin
    ; First algorithm (p. 27): weak Frobenius normal form ("Polynomial pi+1 divides pi" is not true)

    ; Use augemented matrix or P ind Pinv separate
    ;P=identityrat(n,/l64)
    ;Pinv=identityrat(n,/l64)
    F=[[F],[identityrat(n,/l64)]]
    
    ; Generate companion matrix by companion matrix through similarity transforms
    ; However pi+1 divides pi is no longer true.
    i1=-1
    repeat begin
        ; i1 points to last row/column of the last companion matrix
        i2=NewCompanion(F,P,Pinv,i1)
        ; i1 points to last row/column of the new companion matrix
        i1=i2
    endrep until i2 eq n-1
    
    ; Extract the Frobenius normal form and the similarity transform
    ; from the augemented matrix. (Not needed when P and Pinv are calculated separately)
    P=F[*,n:*]
    F=F[*,0:n-1]
    Pinv=ratmatinvert(P)
    
    error=~IsFrobenius(F,/weak)
endif else begin
    ; Second algorithm (p. 57): Frobenius normal form

    AA=F
    ;F=[[F],[identityrat(n,/l64)]] ; for testing

    ; Intermediate form (companion matrices on the diagonal) and P
    P=make_ratarray(s,/l64)
    FormInterm,AA,F,P,blockstart
    
    ; Set the blocks next to the companion matrices to zero without changing the companion matrices
    F=[[F],[P]]
    ;if ~ratarray_equal(F[*,n:*],P) then stop  ; for testing
    RightZeroFormInterm,F,blockstart
    P=F[*,n:*]
    F=F[*,0:n-1]
    
    ; Get inverse of P (which is an integer array)
    Pinv=ratmatinvert(P)
    
    error=~IsFrobenius(F)
endelse

return,F
end;function FrobeniusDecomp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
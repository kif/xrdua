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

function EXTRAC, Array, P0,P1,P2,P3,P4,P5,P6,P7,P8,P9,P10,P11,P12,P13,P14,P15

on_error, 2            ; Return to caller on error

asize = SIZE(Array)
ndim = asize[0]
orig_dims = asize[1:ndim]
ndim_org = ndim

; XRDUA change:
ndiminput=(n_params()-1)/2
if ndim lt ndiminput then begin
    orig_dims=[orig_dims,replicate(1,ndiminput-ndim)]
    ndim=ndiminput
endif

; Is it an array?
if (ndim eq 0) then message, 'Target argument must be an array.'

; Is there an appropriate number of arguments present?
if (n_params() ne (ndim * 2 + 1)) then message, 'Wrong number of arguments.'

; Convert the arguments to a more convenient form.
args = lonarr(2 * ndim)
CASE (ndim) of
; 8: BEGIN & args[15] = P15 & args[14] = P14 & GOTO, do_seven & END
  7: do_seven: BEGIN & args[13] = P13 & args[12] = P12 & GOTO, do_six & END
  6: do_six: BEGIN & args[11] = P11 & args[10] = P10 & GOTO, do_five & END
  5: do_five: BEGIN & args[9] = P9 & args[8] = P8 & GOTO, do_four & END
  4: do_four: BEGIN & args[7] = P7 & args[6] = P6 & GOTO, do_three & END
  3: do_three: BEGIN & args[5] = P5 & args[4] = P4 & GOTO, do_two & END
  2: do_two: BEGIN & args[3] = P3 & args[2] = P2 & GOTO, do_one & END
  1: do_one: BEGIN & args[1] = P1 & args[0] = P0 & END
ENDCASE
srt = args[0:ndim-1]
dims = args[ndim:*]

; Determine if the subarray extends beyond the edges of the original.
; If not, a simple expression will do the job.
srt_over = where(srt lt 0, s_cnt)
dims_over = where((srt + dims) gt orig_dims, b_cnt)

if ((s_cnt eq 0) and (b_cnt eq 0)) then begin
  ; The extracted array does not go beyond the array boundaries. Use the
  ; normal expression to extract the result.
  bnd = dims + srt - 1
  case (ndim) of
    1: result =  Array[P0:bnd[0]]
    2: result =  Array[P0:bnd[0],P1:bnd[1]]
    3: result =  Array[P0:bnd[0],P1:bnd[1],P2:bnd[2]]
    4: result =  Array[P0:bnd[0],P1:bnd[1],P2:bnd[2],P3:bnd[3]]
    5: result =  Array[P0:bnd[0],P1:bnd[1],P2:bnd[2],P3:bnd[3],P4:bnd[4]]
    6: result =  Array[P0:bnd[0],P1:bnd[1],P2:bnd[2],P3:bnd[3],P4:bnd[4], $
               P5:bnd[5]]
    7: result =  Array[P0:bnd[0],P1:bnd[1],P2:bnd[2],P3:bnd[3],P4:bnd[4], $
               P5:bnd[5],P6:bnd[6]]
;   8: result =  Array[P0:bnd[0],P1:bnd[1],P2:bnd[2],P3:bnd[3],P4:bnd[4], $
;               P5:bnd[5],P6:bnd[6],P7:bnd[7]]
   endcase
  goto, done
endif

; If we get this far, the sub array extends beyond the source array
; dimensions. Get a zeroed array of the correct type and extract the
; non-zero part of the original into it.
result = make_array(type=asize[ndim_org + 1], dimension=dims)  ; Get a result array

; Determine the insertion point for the subarray.
isrt = lonarr(ndim)
if (s_cnt ne 0) then isrt[srt_over] = abs(srt[srt_over])


; If any of the starting points exceed the dimensions, then we're done.
dims_over = where(isrt ge dims, b_cnt)
srt_over = where(srt ge orig_dims, s_cnt)
if ((b_cnt ne 0) or (s_cnt ne 0)) then goto, done

; Determine the size of the subarray to be inserted. This is the
; lesser of the original size and the room for insertion in the target.
;
; dims - isrt is the availible room in the result array
; orig_dims - srt - is the largest possible subarray we can pull out of ARRAY

t1 = dims - isrt
bnd = (t1 < (orig_dims - srt)) < orig_dims    ; Minimum of the two sizes
srt = srt > 0                ; Clip starting point to non-negative
bnd = srt + bnd - 1            ; Calcualte the actual outer boundary

; Insert the subarray from ARRAY into RESULT
case (ndim) of
  1: result[isrt[0]] =  Array[srt[0]:bnd[0]]
  2: result[isrt[0],isrt[1]] =  Array[srt[0]:bnd[0],srt[1]:bnd[1]]
  3: result[isrt[0],isrt[1],isrt[2]] = Array[srt[0]:bnd[0],srt[1]:bnd[1],srt[2]:bnd[2]]
  4: result[isrt[0],isrt[1],isrt[2],isrt[3]] =  $
    Array[srt[0]:bnd[0],srt[1]:bnd[1],srt[2]:bnd[2],srt[3]:bnd[3]]
  5: result[isrt[0],isrt[1],isrt[2],isrt[3],isrt[4]] =  $
    Array[srt[0]:bnd[0],srt[1]:bnd[1],srt[2]:bnd[2],srt[3]:bnd[3], $
          srt[4]:bnd[4]]
  6: result[isrt[0],isrt[1],isrt[2],isrt[3],isrt[4],isrt[5]] = $
    Array[srt[0]:bnd[0],srt[1]:bnd[1],srt[2]:bnd[2],srt[3]:bnd[3], $
          srt[4]:bnd[4], srt[5]:bnd[5]]
  7: result[isrt[0],isrt[1],isrt[2],isrt[3],isrt[4],isrt[5],isrt[6]] = $
    Array[srt[0]:bnd[0],srt[1]:bnd[1],srt[2]:bnd[2],srt[3]:bnd[3], $
          srt[4]:bnd[4], srt[5]:bnd[5],srt[6]:bnd[6]]
; 8: result[isrt[0],isrt[1],isrt[2],isrt[3],isrt[4],isrt[5],isrt[6],isrt[7]]= $
;    Array[srt[0]:bnd[0],srt[1]:bnd[1],srt[2]:bnd[2],srt[3]:bnd[3], $
;          srt[4]:bnd[4], srt[5]:bnd[5],srt[6]:bnd[6],srt[7]:bnd[7]]
endcase


done:
  return, result
end;function EXTRAC
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PRO ARROW, x0, y0, x1, y1, HSIZE = hsize, COLOR = color, HTHICK = hthick, $
    THICK = thick, DATA = data, NORMALIZED = norm, $
    SOLID = solid, NOCLIP=NOCLIP, pow=pow

COMPILE_OPT idl2

ON_ERROR, 2
;  Set up keyword params

if keyword_set(pow) then begin
    y0^=pow
    y1^=pow
endif


if n_elements(thick) eq 0 then thick = 1.
if n_elements(hthick) eq 0 then hthick = thick

                ;Head size in device units
if n_elements(hsize) eq 0 then arrowsize = !d.x_size/100. * (hthick/2. > 1) $
    else arrowsize = float(hsize)
if n_elements(color) eq 0 then color = !P.color

mcost = -.866d        ;We use 30 degrees for head angle
sint = .500d
msint = - sint

for i = 0L, n_elements(x0)-1 do begin        ;Each vector
    if keyword_set(data) then $        ;Convert?
        p = convert_coord([x0[i],x1[i]],[y0[i],y1[i]], /data, /to_dev) $
    else if keyword_set(norm) then $
        p = convert_coord([x0[i],x1[i]],[y0[i],y1[i]], /norm, /to_dev) $
    else p = [[x0[i], y0[i]],[x1[i], y1[i]]]

    xp0 = p[0,0]
    xp1 = p[0,1]
    yp0 = p[1,0]
    yp1 = p[1,1]

    dx = xp1 - xp0
    dy = yp1 - yp0
    zz = sqrt(dx^2d + dy^2d)    ;Length

    if zz gt 0 then begin
        dx = dx/zz        ;Cos th
        dy = dy/zz        ;Sin th
    endif else begin
        dx = 1.
        dy = 0.
        zz = 1.
    endelse
    if arrowsize gt 0 then a = arrowsize $  ;a = length of head
    else a = -zz * arrowsize

    xxp0 = xp1 + a * (dx*mcost - dy * msint)
    yyp0 = yp1 + a * (dx*msint + dy * mcost)
    xxp1 = xp1 + a * (dx*mcost - dy * sint)
    yyp1 = yp1 + a * (dx*sint  + dy * mcost)

    if keyword_set(solid) then begin    ;Use polyfill?
      b = a * mcost*.9d    ;End of arrow shaft (Fudge to force join)
      plots, [xp0, xp1+b*dx], [yp0, yp1+b*dy], /DEVICE, $
        COLOR = color, THICK = thick, NOCLIP=NOCLIP
      polyfill, [xxp0, xxp1, xp1, xxp0], [yyp0, yyp1, yp1, yyp0], $
        /DEVICE, COLOR = color, NOCLIP=NOCLIP
    endif else begin
      plots, [xp0, xp1], [yp0, yp1], /DEVICE, COLOR = color, THICK = thick, NOCLIP=NOCLIP
      plots, [xxp0,xp1,xxp1],[yyp0,yp1,yyp1], /DEVICE, COLOR = color, $
            THICK = hthick, NOCLIP=NOCLIP
    endelse
    ENDFOR
end;pro ARROW
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Identity,N_,_EXTRA=re
if n_elements(N_) ne 1 then begin
    N=n_elements(N_)
    val=N_
endif else begin
    N=N_
    val=1
endelse
Array=MAKE_ARRAY(N,N,_EXTRA=re)
Array[LINDGEN(N) * (N+1)] = val
return,Array
end;function Identity
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
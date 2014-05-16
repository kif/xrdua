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

; 45 ASU functions

function ASU_111,h,k,l
    return, (l gt 0 or (l eq 0 and (h gt 0 or (h eq 0 and k ge 0))))
end
function ASU_112,h,k,l
    return, (l ge 0 and (h gt 0 or (h eq 0 and k ge 0)))
end
function ASU_121,h,k,l
    return, (k ge 0 and (l gt 0 or (l eq 0 and h ge 0)))
end
function ASU_211,h,k,l
    return, (h ge 0 and (k gt 0 or (k eq 0 and l ge 0)))
end
function ASU_21U,h,k,l
    return, (h+k ge 0 and (l gt 0 or (l eq 0 and h-k ge 0)))
end
function ASU_21V,h,k,l
    return, (l+h ge 0 and (k gt 0 or (k eq 0 and l-h ge 0)))
end
function ASU_21W,h,k,l
    return, (k+l ge 0 and (h gt 0 or (h eq 0 and k-l ge 0)))
end
function ASU_21X,h,k,l
    return, (h-k ge 0 and (l gt 0 or (l eq 0 and h+k ge 0)))
end
function ASU_21Y,h,k,l
    return, (l-h ge 0 and (k gt 0 or (k eq 0 and l+h ge 0)))
end
function ASU_21Z,h,k,l
    return, (k-l ge 0 and (h gt 0 or (h eq 0 and k+l ge 0)))
end
function ASU_222,h,k,l
    return, (h ge 0 and k ge 0 and l ge 0)
end
function ASU_22U,h,k,l
    return, (h le k and h ge -k and l ge 0)
end
function ASU_22V,h,k,l
    return, (l le h and l ge -h and k ge 0)
end
function ASU_22W,h,k,l
    return, (k le l and k ge -l and h ge 0)
end
function ASU_114,h,k,l
    return, (l ge 0 and ((h ge 0 and k gt 0) or (h eq 0 and k eq 0)))
end
function ASU_141,h,k,l
    return, (k ge 0 and ((l ge 0 and h gt 0) or (l eq 0 and h eq 0)))
end
function ASU_411,h,k,l
    return, (h ge 0 and ((k ge 0 and l gt 0) or (k eq 0 and l eq 0)))
end
function ASU_224,h,k,l
    return, (h ge k and k ge 0 and l ge 0)
end
function ASU_242,h,k,l
    return, (l ge h and h ge 0 and k ge 0)
end
function ASU_422,h,k,l
    return, (k ge l and l ge 0 and h ge 0)
end
function ASU_113,h,k,l
    return, (h ge 0 and k gt 0) or (h eq 0 and k eq 0 and l  ge  0)
end
function ASU_131,h,k,l
    return, (l ge 0 and h gt 0) or (l eq 0 and h eq 0 and k  ge  0)
end
function ASU_311,h,k,l
    return, (k ge 0 and l gt 0) or (k eq 0 and l eq 0 and h  ge  0)
end
function ASU_11T,h,k,l
    return, (h le 0 and k gt 0) or (h eq 0 and k eq 0 and l  ge  0)
end
function ASU_1T1,h,k,l
    return, (l le 0 and h gt 0) or (l eq 0 and h eq 0 and k  ge  0)
end
function ASU_T11,h,k,l
    return, (k le 0 and l gt 0) or (k eq 0 and l eq 0 and h  ge  0)
end
function ASU_31A,h,k,l
    return, (k-l ge 0 and l-h gt 0) or (h eq l and k eq l and h+k+l ge 0)
end
function ASU_31B,h,k,l
    return, (k-l ge 0 and l+h gt 0) or (-h eq l and k eq l and -h+k+l ge 0)
end
function ASU_31C,h,k,l
    return, (-k-l ge 0 and l-h gt 0) or (h eq l and -k eq l and h-k+l ge 0)
end
function ASU_31D,h,k,l
    return, (k+l ge 0 and -l-h gt 0) or (h eq -l and k eq -l and h+k-l ge 0)
end
function ASU_223,h,k,l
    return, (h ge k and k ge 0 and (k gt 0 or l ge 0))
end
function ASU_232,h,k,l
    return, (l ge h and h ge 0 and (h gt 0 or k ge 0))
end
function ASU_322,h,k,l
    return, (k ge l and l ge 0 and (l gt 0 or h ge 0))
end
function ASU_32A,h,k,l
    return, (h ge k and k+l ge h+h and (k+l gt h+h or h+k+l ge 0))
end
function ASU_32B,h,k,l
    return, (-h ge k and k+l ge -h-h and (k+l gt -h-h or -h+k+l ge 0))
end
function ASU_32C,h,k,l
    return, (h ge -k and -k+l ge h+h and (-k+l gt h+h or h-k+l ge 0))
end
function ASU_32D,h,k,l
    return, (h ge k and k-l ge h+h and (k-l gt h+h or h+k-l ge 0))
end
function ASU_32U,h,k,l
    return, (h ge k and k ge 0 and (h gt k or l ge 0))
end
function ASU_32V,h,k,l
    return, (k ge l and l ge 0 and (k gt l or h ge 0))
end
function ASU_32W,h,k,l
    return, (l ge h and h ge 0 and (l gt h or k ge 0))
end
function ASU_32X,h,k,l
    return, (-h ge k and k ge 0 and (-h gt k or l ge 0))
end
function ASU_32Y,h,k,l
    return, (-k ge l and l ge 0 and (-k gt l or h ge 0))
end
function ASU_32Z,h,k,l
    return, (-l ge h and h ge 0 and (-l gt h or k ge 0))
end
function ASU_M3B,h,k,l
    return, (h ge 0 and ((l ge h and k gt h) or (l eq h and k eq h)))
end
function ASU_M3M,h,k,l
    return, (k ge l and l ge h and h ge 0)
end

function ASUfn,name,h,k,l
; Reciprocal space asymmetric unit: determined by the Laue symmetry of the spacegroup.
; Test if hkl is in default reciprocal ASU.
; ABCD are 4 dirns for body diagonal, UVWXYZ are 6 dirns for face diagonal

return,call_function(name,h,k,l)
end;function ASUfn
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pgdata,in,type

;    => 230 spacegroups tabulated or reference settings
;        => 530 spacegroup settings
;    => 32 crystalographic pointgroups
;        => 85 pointgroup settings
;    => 11 laue groups (and ASU's)
;        => 51 laue group settings (with 45 unique ASU's)
;    => 7 crystal systems

; Table origin (clipper/clipper/core/spacegroup_data.cpp):
;    => original table contained 51 Laue group settings with 45 corresponding ASU's
;       (refered too as "pathological", meaning there are more then necessary)
;    => some aren't derived from the 530 spacegroup settings
;    => non-laue groups are added to the original laue group tables to get 85 point groups settings:
;             - Hash code generation:
; symbx="-P 1" & allops=ExpandOpsList(SGFromHall(symbx)) & allops=allops[sort(allops)] & print,baseconvert(hashCRC32(allops),10,16)
;             - General positions check with Bilbao Server:
; symbx="-P 1" & allops=ExpandOpsList(SGFromHall(symbx)) & allops=allops[sort(allops)] & for i=0l,n_elements(allops)-1 do print,StringRTop(Symop(allops[i]))
;    => comment format: Hall, ASU(, crystal system, pointgroups from this system with brackets around the laue group)

tmp={strdefPG,pghash:0L,hall:'',hm:'',crystsys:'',sfname:'',ASU:''}
pgdata = [$
  { strdefPG,'2df60a45'XL, "-P 1",                    '-1' ,    'Triclinic','Ci','ASU_111'},     $; Original Laue setting
  { strdefPG,'c704dd7b'XL, "P 1",                    '1' ,    'Triclinic','C1','ASU_111'},     $
  { strdefPG,'9b92779f'XL, "-P 2",                    '2/m',     'Monoclinic','C2h','ASU_112'},    $; Original Laue setting
  { strdefPG,'068dfeb3'XL, "-P 2y",                    '2/m',     'Monoclinic','C2h','ASU_121'},    $; Original Laue setting
  { strdefPG,'b08d46cf'XL, "-P 2x",                    '2/m',     'Monoclinic','C2h','ASU_211'},    $; Original Laue setting
  { strdefPG,'12339040'XL, "-P 2"+'"',                '2/m',     'Monoclinic','C2h','ASU_21U'},    $; Original Laue setting
  { strdefPG,'44aa9a14'XL, "-P 2y"+'"',                '2/m',     'Monoclinic','C2h','ASU_21V'},    $; Original Laue setting
  { strdefPG,'53e4b366'XL, "-P 2x"+'"',                '2/m',     'Monoclinic','C2h','ASU_21W'},    $; Original Laue setting
  { strdefPG,'4321a07d'XL, "-P 2'",                 '2/m',     'Monoclinic','C2h','ASU_21X'},    $; Original Laue setting
  { strdefPG,'1c1b5411'XL, "-P 2y'",                '2/m',     'Monoclinic','C2h','ASU_21Y'},    $; Original Laue setting
  { strdefPG,'e34a99ed'XL, "-P 2x'",                '2/m',     'Monoclinic','C2h','ASU_21Z'},    $; Original Laue setting
  { strdefPG,'90a34743'XL, "P 2",                    '2',     'Monoclinic','C2','ASU_112'},    $
  { strdefPG,'902cf668'XL, "P 2y",                    '2',     'Monoclinic','C2','ASU_121'},    $
  { strdefPG,'9b02bfbf'XL, "P 2x",                    '2',     'Monoclinic','C2','ASU_211'},    $
  { strdefPG,'939daf66'XL, "P -2",                    'm',     'Monoclinic','Cs','ASU_112'},    $
  { strdefPG,'93121e4d'XL, "P -2y",                    'm',     'Monoclinic','Cs','ASU_121'},    $
  { strdefPG,'6dd7504b'XL, "P -2x",                    'm',     'Monoclinic','Cs','ASU_211'},    $
  { strdefPG,'e7243bbc'XL, "-P 2 2",                'mmm',     'Orthorhombic','D2h','ASU_222'},    $; Original Laue setting
  { strdefPG,'885920bf'XL, "-P 2 2"+'"',            'mmm',     'Orthorhombic','D2h','ASU_22U'},    $; Original Laue setting
  { strdefPG,'e980f874'XL, "-P 2 2"+'"'+"(y,z,x)",    'mmm',     'Orthorhombic','D2h','ASU_22V'},    $; Original Laue setting
  { strdefPG,'05c7f86e'XL, "-P 2 2"+'"'+"(z,x,y)",    'mmm',     'Orthorhombic','D2h','ASU_22W'},    $; Original Laue setting
  { strdefPG,'a959fc0b'XL, "P 2 2",                    '222',     'Orthorhombic','D2','ASU_222'},    $
  { strdefPG,'715a8a38'XL, "P 2 -2",                'mm2',     'Orthorhombic','D2v','ASU_222'},    $
  { strdefPG,'34467527'XL, "P -2 2",                'mm2',     'Orthorhombic','D2v','ASU_222'},    $
  { strdefPG,'454292ce'XL, "P -2 -2",                'mm2',     'Orthorhombic','D2v','ASU_222'},    $
  { strdefPG,'fdd759b5'XL, "-P 3",                    '-3',     'Trigonal','C3i','ASU_113'},     $; Original Laue setting
  { strdefPG,'f2769c28'XL, "-P 3 (y,z,x)",            '-3',     'Trigonal','C3i','ASU_131'},     $; Original Laue setting
  { strdefPG,'cd4b8428'XL, "-P 3 (z,x,y)",            '-3',     'Trigonal','C3i','ASU_311'},     $; Original Laue setting
  { strdefPG,'07fa5ca1'XL, "-P 3 (-x,y,z)",            '-3',     'Trigonal','C3i','ASU_11T'},     $; Original Laue setting
  { strdefPG,'0f070468'XL, "-P 3 (y,z,-x)",            '-3',     'Trigonal','C3i','ASU_1T1'},     $; Original Laue setting
  { strdefPG,'9bd3dcf4'XL, "-P 3 (z,-x,y)",            '-3',     'Trigonal','C3i','ASU_T11'},     $; Original Laue setting
  { strdefPG,'d9a29bac'XL, "-P 3*",                    '-3',     'Trigonal','C3i','ASU_31A'},     $; Original Laue setting
  { strdefPG,'627a2f8c'XL, "-P 3* (-x,y,z)",        '-3',     'Trigonal','C3i','ASU_31B'},     $; Original Laue setting
  { strdefPG,'cfb71d71'XL, "-P 3* (x,-y,z)",        '-3',     'Trigonal','C3i','ASU_31C'},     $; Original Laue setting
  { strdefPG,'8a86426e'XL, "-P 3* (x,y,-z)",        '-3',     'Trigonal','C3i','ASU_31D'},     $; Original Laue setting
  { strdefPG,'329c980e'XL, "P 3",                    '3',     'Trigonal','C3','ASU_113'},     $
  { strdefPG,'d5a0aa2d'XL, "P 3*",                    '3',     'Trigonal','C3','ASU_31A'},     $
  { strdefPG,'f74c7f83'XL, "-P 3 2",                '-3m',     'Trigonal','D3d','ASU_223'},    $; Original Laue setting
  { strdefPG,'573b981c'XL, "-P 3 2 (y,z,x)",        '-3m',     'Trigonal','D3d','ASU_232'},    $; Original Laue setting
  { strdefPG,'1799544d'XL, "-P 3 2 (z,x,y)",        '-3m',     'Trigonal','D3d','ASU_322'},    $; Original Laue setting
  { strdefPG,'1c80e47a'XL, "-P 3* 2",                '-3m',     'Trigonal','D3d','ASU_32A'},    $; Original Laue setting
  { strdefPG,'ea7284da'XL, "-P 3* 2 (-x,y,z)",        '-3m',     'Trigonal','D3d','ASU_32B'},    $; Original Laue setting
  { strdefPG,'b193db73'XL, "-P 3* 2 (x,-y,z)",        '-3m',     'Trigonal','D3d','ASU_32C'},    $; Original Laue setting
  { strdefPG,'04fecdf9'XL, "-P 3* 2 (-x,-y,z)",        '-3m',     'Trigonal','D3d','ASU_32D'},    $; Original Laue setting
  { strdefPG,'fc3edafb'XL, "-P 3 2"+'"',            '-3m',     'Trigonal','D3d','ASU_32U'},    $; Original Laue setting
  { strdefPG,'d60d11a0'XL, "-P 3 2"+'"'+"(z,x,y)",    '-3m',     'Trigonal','D3d','ASU_32V'},    $; Original Laue setting
  { strdefPG,'f7d5112f'XL, "-P 3 2"+'"'+"(y,z,x)",    '-3m',     'Trigonal','D3d','ASU_32W'},    $; Original Laue setting
  { strdefPG,'fbb8f18b'XL, "-P 3 2"+'"'+"(-x,y,z)",'-3m',     'Trigonal','D3d','ASU_32X'},    $; Original Laue setting
  { strdefPG,'530fcba9'XL, "-P 3 2"+'"'+"(z,-x,y)",'-3m',     'Trigonal','D3d','ASU_32Y'},    $; Original Laue setting
  { strdefPG,'a3d49592'XL, "-P 3 2"+'"'+"(y,z,-x)",'-3m',     'Trigonal','D3d','ASU_32Z'},    $; Original Laue setting
  { strdefPG,'65b7a72b'XL, "P 3 2",                    '32',     'Trigonal','D3','ASU_223'},    $
  { strdefPG,'a20b8591'XL, "P 3* 2",                '32',     'Trigonal','D3','ASU_32A'},    $
  { strdefPG,'c1840a7a'XL, "P 3 2"+'"',                '32',     'Trigonal','D3','ASU_32U'},    $
  { strdefPG,'39859b12'XL, "P 3 -2",                '3m',     'Trigonal','C3v','ASU_223'},    $
  { strdefPG,'b951b4f7'XL, "P 3* -2",                '3m',     'Trigonal','C3v','ASU_32A'},    $
  { strdefPG,'9f4cffaa'XL, "P 3 -2"+'"',            '3m',     'Trigonal','C3v','ASU_32U'},    $
  { strdefPG,'7f473453'XL, "-P 4",                    '4/m',     'Tetragonal','C4h','ASU_114'},    $; Original Laue setting
  { strdefPG,'081d78e5'XL, "-P 4 (y,z,x)",            '4/m',     'Tetragonal','C4h','ASU_141'},    $; Original Laue setting
  { strdefPG,'6a9fa6e5'XL, "-P 4 (z,x,y)",             '4/m',     'Tetragonal','C4h','ASU_411'},    $; Original Laue setting
  { strdefPG,'e194247e'XL, "P 4",                    '4',     'Tetragonal','C4','ASU_114'},    $
  { strdefPG,'5b1380ec'XL, "P -4",                    '-4',     'Tetragonal','S4','ASU_114'},    $
  { strdefPG,'b8f2113e'XL, "-P 4 2",                '4/mmm','Tetragonal','D4h','ASU_224'},  $; Original Laue setting
  { strdefPG,'b9ce4369'XL, "-P 4 2 (y,z,x)",         '4/mmm','Tetragonal','D4h','ASU_242'},  $; Original Laue setting
  { strdefPG,'c39787b7'XL, "-P 4 2 (z,x,y)",        '4/mmm','Tetragonal','D4h','ASU_422'},  $; Original Laue setting
  { strdefPG,'781566a5'XL, "P 4 2",                    '422',     'Tetragonal','D4','ASU_224'},    $
  { strdefPG,'93373957'XL, "P 4 -2",                '4mm',     'Tetragonal','C4v','ASU_224'},    $
  { strdefPG,'ae367c74'XL, "P -4 2",                '-42m',    'Tetragonal','D2d','ASU_224'},    $
  { strdefPG,'d50d13d5'XL, "P -4 -2",                '-42m',    'Tetragonal','D2d','ASU_224'},    $
  { strdefPG,'32dacfb6'XL, "-P 6",                     '6/m',     'Hexagonal','C6h','ASU_114'},    $; Original Laue setting
  { strdefPG,'e7987b0c'XL, "-P 6 (y,z,x)",             '6/m',     'Hexagonal','C6h','ASU_141'},    $; Original Laue setting
  { strdefPG,'b5b69658'XL, "-P 6 (z,x,y)",             '6/m',     'Hexagonal','C6h','ASU_411'},    $; Original Laue setting
  { strdefPG,'a2ddaf47'XL, "P 6",                     '6',     'Hexagonal','C6','ASU_114'},    $
  { strdefPG,'35eb6ccb'XL, "P -6",                     '-6',     'Hexagonal','C3h','ASU_114'},    $
  { strdefPG,'f1fc7952'XL, "-P 6 2",                '6/mmm','Hexagonal','D6h','ASU_224'},  $; Original Laue setting
  { strdefPG,'386f0ab4'XL, "-P 6 2 (y,z,x)",        '6/mmm','Hexagonal','D6h','ASU_242'},  $; Original Laue setting
  { strdefPG,'cd531b66'XL, "-P 6 2 (z,x,y)",         '6/mmm','Hexagonal','D6h','ASU_422'},  $; Original Laue setting
  { strdefPG,'90230e5f'XL, "P 6 2",                    '622',    'Hexagonal','D6','ASU_224'},  $
  { strdefPG,'3b0a2d17'XL, "P 6 -2",                '6mmm',    'Hexagonal','C6v','ASU_224'},  $
  { strdefPG,'eeacb736'XL, "P -6 2",                '-62m',    'Hexagonal','D3h','ASU_224'},  $
  { strdefPG,'55b5be6a'XL, "P -6 -2",                '-62m',    'Hexagonal','D3h','ASU_224'},  $
  { strdefPG,'72e55913'XL, "-P 2 2 3",                 'm-3',     'Cubic','Th','ASU_M3B'},    $; Original Laue setting
  { strdefPG,'5843870d'XL, "P 2 2 3",                 '23',     'Cubic','T','ASU_M3B'},    $
  { strdefPG,'74c407d3'XL, "-P 4 2 3",                'm-3m', 'Cubic','Oh','ASU_M3M'},    $; Original Laue setting
  { strdefPG,'93a6edeb'XL, "P 4 2 3",                '432',     'Cubic','O','ASU_M3M'},    $
  { strdefPG,'6fc31c7b'XL, "P -4 2 3",                '-43m', 'Cubic','Td','ASU_M3M'}]

; Type:
;    0: Hash in, all out

case type of
0:    begin
    ind=where(pgdata[*].pghash eq In[0],ct)
    case ct of
    0:    return,tmp
    else:    return,pgdata[ind[0]]
    endcase
    endcase
1:    return,pgdata.hall
endcase

end;function pgdata
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SGSFSymbol,In
; http://cci.lbl.gov/sginfo/ (sginfo.h)

Mat=[   '',    'C1^1',    'Ci^1',    'C2^1',    'C2^2',    'C2^3',$
    'Cs^1',    'Cs^2',    'Cs^3',    'Cs^4',   'C2h^1',   'C2h^2',$
   'C2h^3',   'C2h^4',   'C2h^5',   'C2h^6',    'D2^1',    'D2^2',$
    'D2^3',    'D2^4',    'D2^5',    'D2^6',    'D2^7',    'D2^8',$
    'D2^9',   'C2v^1',   'C2v^2',   'C2v^3',   'C2v^4',   'C2v^5',$
   'C2v^6',   'C2v^7',   'C2v^8',   'C2v^9',  'C2v^10',  'C2v^11',$
  'C2v^12',  'C2v^13',  'C2v^14',  'C2v^15',  'C2v^16',  'C2v^17',$
  'C2v^18',  'C2v^19',  'C2v^20',  'C2v^21',  'C2v^22',   'D2h^1',$
   'D2h^2',   'D2h^3',   'D2h^4',   'D2h^5',   'D2h^6',   'D2h^7',$
   'D2h^8',   'D2h^9',  'D2h^10',  'D2h^11',  'D2h^12',  'D2h^13',$
  'D2h^14',  'D2h^15',  'D2h^16',  'D2h^17',  'D2h^18',  'D2h^19',$
  'D2h^20',  'D2h^21',  'D2h^22',  'D2h^23',  'D2h^24',  'D2h^25',$
  'D2h^26',  'D2h^27',  'D2h^28',    'C4^1',    'C4^2',    'C4^3',$
    'C4^4',    'C4^5',    'C4^6',    'S4^1',    'S4^2',   'C4h^1',$
   'C4h^2',   'C4h^3',   'C4h^4',   'C4h^5',   'C4h^6',    'D4^1',$
    'D4^2',    'D4^3',    'D4^4',    'D4^5',    'D4^6',    'D4^7',$
    'D4^8',    'D4^9',   'D4^10',   'C4v^1',   'C4v^2',   'C4v^3',$
   'C4v^4',   'C4v^5',   'C4v^6',   'C4v^7',   'C4v^8',   'C4v^9',$
  'C4v^10',  'C4v^11',  'C4v^12',   'D2d^1',   'D2d^2',   'D2d^3',$
   'D2d^4',   'D2d^5',   'D2d^6',   'D2d^7',   'D2d^8',   'D2d^9',$
  'D2d^10',  'D2d^11',  'D2d^12',   'D4h^1',   'D4h^2',   'D4h^3',$
   'D4h^4',   'D4h^5',   'D4h^6',   'D4h^7',   'D4h^8',   'D4h^9',$
  'D4h^10',  'D4h^11',  'D4h^12',  'D4h^13',  'D4h^14',  'D4h^15',$
  'D4h^16',  'D4h^17',  'D4h^18',  'D4h^19',  'D4h^20',    'C3^1',$
    'C3^2',    'C3^3',    'C3^4',   'C3i^1',   'C3i^2',    'D3^1',$
    'D3^2',    'D3^3',    'D3^4',    'D3^5',    'D3^6',    'D3^7',$
   'C3v^1',   'C3v^2',   'C3v^3',   'C3v^4',   'C3v^5',   'C3v^6',$
   'D3d^1',   'D3d^2',   'D3d^3',   'D3d^4',   'D3d^5',   'D3d^6',$
    'C6^1',    'C6^2',    'C6^3',    'C6^4',    'C6^5',    'C6^6',$
   'C3h^1',   'C6h^1',   'C6h^2',    'D6^1',    'D6^2',    'D6^3',$
    'D6^4',    'D6^5',    'D6^6',   'C6v^1',   'C6v^2',   'C6v^3',$
   'C6v^4',   'D3h^1',   'D3h^2',   'D3h^3',   'D3h^4',   'D6h^1',$
   'D6h^2',   'D6h^3',   'D6h^4',     'T^1',     'T^2',     'T^3',$
     'T^4',     'T^5',    'Th^1',    'Th^2',    'Th^3',    'Th^4',$
    'Th^5',    'Th^6',    'Th^7',     'O^1',     'O^2',     'O^3',$
     'O^4',     'O^5',     'O^6',     'O^7',     'O^8',    'Td^1',$
    'Td^2',    'Td^3',    'Td^4',    'Td^5',    'Td^6',    'Oh^1',$
    'Oh^2',    'Oh^3',    'Oh^4',    'Oh^5',    'Oh^6',    'Oh^7',$
    'Oh^8',    'Oh^9',   'Oh^10',        '']

nmax=230

if size(In,/type) eq 7 then begin
    In=strlowcase(In)
    In=strupcase(strmid(In,0,1))+strmid(In,1)

    ind=where(Mat eq In[0],ct)
    case ct of
    0:    return,0
    else:    return,ind[0]
    endcase

endif else return, Mat[0>In<(nmax+1)]

end;function SGSFSymbol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function shmsymbol, In
; http://icsd.ccp14.ac.uk/icsd/lazy_spcgrp.html

Mat=[  '','P 1',       'P -1',        'P 2',       'P 21',        'C 2',        'P M',$
        'P C',        'C M',        'C C',      'P 2/M',     'P 21/M',      'C 2/M',$
      'P 2/C',     'P 21/C',      'C 2/C',    'P 2 2 2',   'P 2 2 21',  'P 21 21 2',$
 'P 21 21 21',   'C 2 2 21',    'C 2 2 2',    'F 2 2 2',    'I 2 2 2', 'I 21 21 21',$
    'P M M 2',   'P M C 21',    'P C C 2',    'P M A 2',   'P C A 21',    'P N C 2',$
   'P M N 21',    'P B A 2',   'P N A 21',    'P N N 2',    'C M M 2',   'C M C 21',$
    'C C C 2',    'A M M 2',    'A B M 2',    'A M A 2',    'A B A 2',    'F M M 2',$
    'F D D 2',    'I M M 2',    'I B A 2',    'I M A 2',    'P M M M',    'P N N N',$
    'P C C M',    'P B A N',    'P M M A',    'P N N A',    'P M N A',    'P C C A',$
    'P B A M',    'P C C N',    'P B C M',    'P N N M',    'P M M N',    'P B C N',$
    'P B C A',    'P N M A',    'C M C M',    'C M C A',    'C M M M',    'C C C M',$
    'C M M A',    'C C C A',    'F M M M',    'F D D D',    'I M M M',    'I B A M',$
    'I B C A',    'I M M A',        'P 4',       'P 41',       'P 42',       'P 43',$
        'I 4',       'I 41',       'P -4',       'I -4',      'P 4/M',     'P 42/M',$
      'P 4/N',     'P 42/N',      'I 4/M',     'I 41/A',    'P 4 2 2',   'P 4 21 2',$
   'P 41 2 2',  'P 41 21 2',   'P 42 2 2',  'P 42 21 2',   'P 43 2 2',  'P 43 21 2',$
    'I 4 2 2',   'I 41 2 2',    'P 4 M M',    'P 4 B M',   'P 42 C M',   'P 42 N M',$
    'P 4 C C',    'P 4 N C',   'P 42 M C',   'P 42 B C',    'I 4 M M',    'I 4 C M',$
   'I 41 M D',   'I 41 C D',   'P -4 2 M',   'P -4 2 C',  'P -4 21 M',  'P -4 21 C',$
   'I -4 M 2',   'P -4 C 2',   'P -4 B 2',   'P -4 N 2',   'P -4 M 2',   'I -4 C 2',$
   'P -4 2 M',   'I -4 2 D',  'P 4/M M M',  'P 4/M C C',  'P 4/N B M',  'P 4/N N C',$
  'P 4/M B M',  'P 4/M N C',  'P 4/N M M',  'P 4/N C C', 'P 42/M M C', 'P 42/M C M',$
 'P 42/N B C', 'P 42/N N M', 'P 42/M B C', 'P 42/M N M', 'P 42/N M C', 'P 42/N C M',$
  'I 4/M M M',  'I 4/M C M', 'I 41/A M D', 'I 41/A C D',        'P 3',       'P 31',$
       'P 32',        'R 3',       'P -3',       'R -3',    'P 3 1 2',    'P 3 2 1',$
   'P 31 1 2',   'P 31 2 1',   'P 32 1 2',   'P 32 2 1',      'R 3 2',    'P 3 M 1',$
    'P 3 1 M',    'P 3 C 1',    'P 3 1 C',      'R 3 M',      'R 3 C',   'P -3 1 M',$
   'P -3 1 C',   'P -3 M 1',   'P -3 C 1',     'R -3 M',     'R -3 C',        'P 6',$
       'P 61',       'P 65',       'P 62',       'P 64',       'P 63',       'P -6',$
      'P 6/M',     'P 63/M',    'P 6 2 2',   'P 61 2 2',   'P 65 2 2',   'P 62 2 2',$
   'P 64 2 2',   'P 63 2 2',    'P 6 M M',    'P 6 C C',   'P 63 C M',   'P 63 M C',$
   'P -6 M 2',   'P -6 C 2',   'P -6 2 M',   'P -6 2 C',  'P 6/M M M',  'P 6/M C C',$
 'P 63/M C M', 'P 63/M M C',      'P 2 3',      'F 2 3',      'I 2 3',     'P 21 3',$
     'I 21 3',      'P M 3',      'P N 3',      'F M 3',      'F D 3',      'I M 3',$
      'P A 3',      'I A 3',    'P 4 3 2',   'P 42 3 2',    'F 4 3 2',   'F 41 3 2',$
    'I 4 3 2',   'P 43 3 2',   'P 41 3 2',   'I 41 3 2',   'P -4 3 M',   'F -4 3 M',$
   'I -4 3 M',   'P -4 3 N',   'F -4 3 C',   'I -4 3 D',    'P M 3 M',    'P N 3 N',$
    'P M 3 N',    'P N 3 M',    'F M 3 M',    'F M 3 C',    'F D 3 M',    'F D 3 C',$
    'I M 3 M',    'I A 3 D','']

nmax=230

if size(In,/type) eq 7 then begin
    In=strupcase(In)

    ind=where(Mat eq In[0],ct)
    case ct of
    0:    return,0
    else:    return,ind[0]
    endcase

endif else return, Mat[0>In<(nmax+1)]

end;function shmsymbol
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sgdata,in,type,pref_12=pref_12,pref_hr=pref_hr

; Data adapted from:
; 
; S.R. Hall
; Space-Group Notation with an Explicit Origin
; Acta Cryst. (1981). A37, 517-525
; (http://cci.lbl.gov/sginfo/hall_symbols.html)

;Monoclinic                  code  =  <unique axis><cell choice>
;     Unique axis choices            b  -b  c  -c  a  -a
;     Cell choices                   1  2  3
;Orthorhombic                code  =  <origin choice><setting>
;     Origin choices                 1  2
;     Setting choices                abc  ba-c  cab  -cba  bca  a-cb
;Tetragonal, Cubic           code  =  <origin choice>
;     Origin choices                 1  2
;Trigonal                    code  =  <cell choice>
;     Cell choices                   H (hex)   R (rhomb)

tmp={strdefSG,sghash:0L,hall:'',hm:'',ext:'',choice:'',num:0,set:0}
; hash, hall symbol, Full Hermann-Mauguin symbol, extra, choice, IT number, tabulated setting number

; Some entrees have the same hash (i.e. same hall and symmetry).
; One of them is taken as default and the others are converted:
; 'C 2 2 -1bc': 'C c c b' > 'C c c a'
; 'A 2 2 -1ac': 'A c a a' > 'A b a a'
; 'B 2 2 -1bc': 'B b a b' > 'B b c b'

sgdata0=[ $
    {strdefSG,'00000000'XL,                '',                '',         '',         '',   0,   0}, $
    {strdefSG,'c704dd7b'XL,             'P 1',             'P 1',         '',         '',   1,   1}, $
    {strdefSG,'2df60a45'XL,            '-P 1',            'P -1',         '',         '',   2,   2}, $
    {strdefSG,'902cf668'XL,            'P 2y',         'P 1 2 1',         '',        'b',   3,   3}, $
    {strdefSG,'90a34743'XL,             'P 2',         'P 1 1 2',         '',        'c',   3,   4}, $
    {strdefSG,'9b02bfbf'XL,            'P 2x',         'P 2 1 1',         '',        'a',   3,   5}, $
    {strdefSG,'da168154'XL,           'P 2yb',        'P 1 21 1',         '',        'b',   4,   6}, $
    {strdefSG,'7851c0d1'XL,            'P 2c',        'P 1 1 21',         '',        'c',   4,   7}, $
    {strdefSG,'ae0e24db'XL,           'P 2xa',        'P 21 1 1',         '',        'a',   4,   8}, $
    {strdefSG,'deb5173b'XL,            'C 2y',         'C 1 2 1',         '',       'b1',   5,   9}, $
    {strdefSG,'1e4dbca6'XL,            'A 2y',         'A 1 2 1',         '',       'b2',   5,  10}, $
    {strdefSG,'6042abb6'XL,            'I 2y',         'I 1 2 1',         '',       'b3',   5,  11}, $
    {strdefSG,'2a55a450'XL,             'A 2',         'A 1 1 2',         '',       'c1',   5,  12}, $
    {strdefSG,'458c0b8c'XL,             'B 2',         'B 1 1 2',         '',       'c2',   5,  13}, $
    {strdefSG,'545ab340'XL,             'I 2',         'I 1 1 2',         '',       'c3',   5,  14}, $
    {strdefSG,'c4ef9a2d'XL,            'B 2x',         'B 2 1 1',         '',       'a1',   5,  15}, $
    {strdefSG,'6bce9e6c'XL,            'C 2x',         'C 2 1 1',         '',       'a2',   5,  16}, $
    {strdefSG,'d53922e1'XL,            'I 2x',         'I 2 1 1',         '',       'a3',   5,  17}, $
    {strdefSG,'93121e4d'XL,           'P -2y',         'P 1 m 1',         '',        'b',   6,  18}, $
    {strdefSG,'939daf66'XL,            'P -2',         'P 1 1 m',         '',        'c',   6,  19}, $
    {strdefSG,'6dd7504b'XL,           'P -2x',         'P m 1 1',         '',        'a',   6,  20}, $
    {strdefSG,'7be099df'XL,          'P -2yc',         'P 1 c 1',         '',       'b1',   7,  21}, $
    {strdefSG,'4eec02bb'XL,         'P -2yac',         'P 1 n 1',         '',       'b2',   7,  22}, $
    {strdefSG,'a61e8529'XL,          'P -2ya',         'P 1 a 1',         '',       'b3',   7,  23}, $
    {strdefSG,'a6913402'XL,           'P -2a',         'P 1 1 a',         '',       'c1',   7,  24}, $
    {strdefSG,'ecab433e'XL,          'P -2ab',         'P 1 1 n',         '',       'c2',   7,  25}, $
    {strdefSG,'d9a7d85a'XL,           'P -2b',         'P 1 1 b',         '',       'c3',   7,  26}, $
    {strdefSG,'27ed2777'XL,          'P -2xb',         'P b 1 1',         '',       'a1',   7,  27}, $
    {strdefSG,'2ac91f43'XL,         'P -2xbc',         'P n 1 1',         '',       'a2',   7,  28}, $
    {strdefSG,'8525d7d9'XL,          'P -2xc',         'P c 1 1',         '',       'a3',   7,  29}, $
    {strdefSG,'77b286e1'XL,           'C -2y',         'C 1 m 1',         '',       'b1',   8,  30}, $
    {strdefSG,'b74a2d7c'XL,           'A -2y',         'A 1 m 1',         '',       'b2',   8,  31}, $
    {strdefSG,'c9453a6c'XL,           'I -2y',         'I 1 m 1',         '',       'b3',   8,  32}, $
    {strdefSG,'8352358a'XL,            'A -2',         'A 1 1 m',         '',       'c1',   8,  33}, $
    {strdefSG,'ec8b9a56'XL,            'B -2',         'B 1 1 m',         '',       'c2',   8,  34}, $
    {strdefSG,'fd5d229a'XL,            'I -2',         'I 1 1 m',         '',       'c3',   8,  35}, $
    {strdefSG,'ddc257bb'XL,           'B -2x',         'B m 1 1',         '',       'a1',   8,  36}, $
    {strdefSG,'72e353fa'XL,           'C -2x',         'C m 1 1',         '',       'a2',   8,  37}, $
    {strdefSG,'cc14ef77'XL,           'I -2x',         'I m 1 1',         '',       'a3',   8,  38}, $
    {strdefSG,'38ef08aa'XL,          'C -2yc',         'C 1 c 1',         '',       'b1',   9,  39}, $
    {strdefSG,'54faf6d0'XL,         'A -2yac',         'A 1 n 1',         '',       'b2',   9,  40}, $
    {strdefSG,'9d54258d'XL,          'I -2ya',         'I 1 a 1',         '',       'b3',   9,  41}, $
    {strdefSG,'e35b329d'XL,          'A -2ya',         'A 1 a 1',         '',      '-b1',   9,  42}, $
    {strdefSG,'6cfe174b'XL,         'C -2ybc',         'C 1 n 1',         '',      '-b2',   9,  43}, $
    {strdefSG,'2af5e1c0'XL,          'I -2yc',         'I 1 c 1',         '',      '-b3',   9,  44}, $
    {strdefSG,'d7432a6b'XL,           'A -2a',         'A 1 1 a',         '',       'c1',   9,  45}, $
    {strdefSG,'0f3b41fa'XL,          'B -2bc',         'B 1 1 n',         '',       'c2',   9,  46}, $
    {strdefSG,'4afce6d7'XL,           'I -2b',         'I 1 1 b',         '',       'c3',   9,  47}, $
    {strdefSG,'5b2a5e1b'XL,           'B -2b',         'B 1 1 b',         '',      '-c1',   9,  48}, $
    {strdefSG,'60e2ee26'XL,          'A -2ac',         'A 1 1 n',         '',      '-c2',   9,  49}, $
    {strdefSG,'a94c3d7b'XL,           'I -2a',         'I 1 1 a',         '',      '-c3',   9,  50}, $
    {strdefSG,'6a6393f6'XL,          'B -2xb',         'B b 1 1',         '',       'a1',   9,  51}, $
    {strdefSG,'69afc250'XL,         'C -2xbc',         'C n 1 1',         '',       'a2',   9,  52}, $
    {strdefSG,'2fa434db'XL,          'I -2xc',         'I c 1 1',         '',       'a3',   9,  53}, $
    {strdefSG,'3dbeddb1'XL,          'C -2xc',         'C c 1 1',         '',      '-a1',   9,  54}, $
    {strdefSG,'3e728c17'XL,         'B -2xbc',         'B n 1 1',         '',      '-a2',   9,  55}, $
    {strdefSG,'7bb52b3a'XL,          'I -2xb',         'I b 1 1',         '',      '-a3',   9,  56}, $
    {strdefSG,'068dfeb3'XL,           '-P 2y',       'P 1 2/m 1',         '',        'b',  10,  57}, $
    {strdefSG,'9b92779f'XL,            '-P 2',       'P 1 1 2/m',         '',        'c',  10,  58}, $
    {strdefSG,'b08d46cf'XL,           '-P 2x',       'P 2/m 1 1',         '',        'a',  10,  59}, $
    {strdefSG,'54fa8558'XL,          '-P 2yb',      'P 1 21/m 1',         '',        'b',  11,  60}, $
    {strdefSG,'31194672'XL,           '-P 2c',      'P 1 1 21/m',         '',        'c',  11,  61}, $
    {strdefSG,'ce8251df'XL,          '-P 2xa',      'P 21/m 1 1',         '',        'a',  11,  62}, $
    {strdefSG,'09efdd05'XL,           '-C 2y',       'C 1 2/m 1',         '',       'b1',  12,  63}, $
    {strdefSG,'8e118804'XL,           '-A 2y',       'A 1 2/m 1',         '',       'b2',  12,  64}, $
    {strdefSG,'5c5d4c9f'XL,           '-I 2y',       'I 1 2/m 1',         '',       'b3',  12,  65}, $
    {strdefSG,'5a48cc10'XL,            '-A 2',       'A 1 1 2/m',         '',       'c1',  12,  66}, $
    {strdefSG,'97e84b5c'XL,            '-B 2',       'B 1 1 2/m',         '',       'c2',  12,  67}, $
    {strdefSG,'8804088b'XL,            '-I 2',       'I 1 1 2/m',         '',       'c3',  12,  68}, $
    {strdefSG,'05a3ecad'XL,           '-B 2x',       'B 2/m 1 1',         '',       'a1',  12,  69}, $
    {strdefSG,'4ffd3ee0'XL,           '-C 2x',       'C 2/m 1 1',         '',       'a2',  12,  70}, $
    {strdefSG,'1a4faf7a'XL,           '-I 2x',       'I 2/m 1 1',         '',       'a3',  12,  71}, $
    {strdefSG,'ac06cf5e'XL,          '-P 2yc',       'P 1 2/c 1',         '',       'b1',  13,  72}, $
    {strdefSG,'f817d0bf'XL,         '-P 2yac',       'P 1 2/n 1',         '',       'b2',  13,  73}, $
    {strdefSG,'529ce152'XL,          '-P 2ya',       'P 1 2/a 1',         '',       'b3',  13,  74}, $
    {strdefSG,'cf83687e'XL,           '-P 2a',       'P 1 1 2/a',         '',       'c1',  13,  75}, $
    {strdefSG,'9df41395'XL,          '-P 2ab',       'P 1 1 2/n',         '',       'c2',  13,  76}, $
    {strdefSG,'c9e50c74'XL,           '-P 2b',       'P 1 1 2/b',         '',       'c3',  13,  77}, $
    {strdefSG,'4aada62e'XL,          '-P 2xb',       'P 2/b 1 1',         '',       'a1',  13,  78}, $
    {strdefSG,'f45a1aa3'XL,         '-P 2xbc',       'P 2/n 1 1',         '',       'a2',  13,  79}, $
    {strdefSG,'e58ca26f'XL,          '-P 2xc',       'P 2/c 1 1',         '',       'a3',  13,  80}, $
    {strdefSG,'41572403'XL,         '-P 2ybc',      'P 1 21/c 1',         '',       'b1',  14,  81}, $
    {strdefSG,'15463be2'XL,          '-P 2yn',      'P 1 21/n 1',         '',       'b2',  14,  82}, $
    {strdefSG,'00eb9ab9'XL,         '-P 2yab',      'P 1 21/a 1',         '',       'b3',  14,  83}, $
    {strdefSG,'65085993'XL,          '-P 2ac',      'P 1 1 21/a',         '',       'c1',  14,  84}, $
    {strdefSG,'8859b2ce'XL,           '-P 2n',      'P 1 1 21/n',         '',       'c2',  14,  85}, $
    {strdefSG,'dc48ad2f'XL,          '-P 2bc',      'P 1 1 21/b',         '',       'c3',  14,  86}, $
    {strdefSG,'34a2b13e'XL,         '-P 2xab',      'P 21/b 1 1',         '',       'a1',  14,  87}, $
    {strdefSG,'8a550db3'XL,          '-P 2xn',      'P 21/n 1 1',         '',       'a2',  14,  88}]
sgdata1=[$
    {strdefSG,'9b83b57f'XL,         '-P 2xac',      'P 21/c 1 1',         '',       'a3',  14,  89}, $
    {strdefSG,'673d531e'XL,          '-C 2yc',       'C 1 2/c 1',         '',       'b1',  15,  90}, $
    {strdefSG,'428e62a8'XL,         '-A 2yac',       'A 1 2/n 1',         '',       'b2',  15,  91}, $
    {strdefSG,'b7f841ea'XL,          '-I 2ya',       'I 1 2/a 1',         '',       'b3',  15,  92}, $
    {strdefSG,'65b48571'XL,          '-A 2ya',       'A 1 2/a 1',         '',      '-b1',  15,  93}, $
    {strdefSG,'8c985e6b'XL,         '-C 2ybc',       'C 1 2/n 1',         '',      '-b2',  15,  94}, $
    {strdefSG,'90c2a633'XL,          '-I 2yc',       'I 1 2/c 1',         '',      '-b3',  15,  95}, $
    {strdefSG,'b1edc165'XL,           '-A 2a',       'A 1 1 2/a',         '',       'c1',  15,  96}, $
    {strdefSG,'5b77a1f0'XL,          '-B 2bc',       'B 1 1 2/n',         '',       'c2',  15,  97}, $
    {strdefSG,'af3eef52'XL,           '-I 2b',       'I 1 1 2/b',         '',       'c3',  15,  98}, $
    {strdefSG,'b0d2ac85'XL,           '-B 2b',       'B 1 1 2/b',         '',      '-c1',  15,  99}, $
    {strdefSG,'96d726bc'XL,          '-A 2ac',       'A 1 1 2/n',         '',      '-c2',  15, 100}, $
    {strdefSG,'63a105fe'XL,           '-I 2a',       'I 1 1 2/a',         '',      '-c3',  15, 101}, $
    {strdefSG,'e4301506'XL,          '-B 2xb',       'B 2/b 1 1',         '',       'a1',  15, 102}, $
    {strdefSG,'aecdc717'XL,         '-C 2xbc',       'C 2/n 1 1',         '',       'a2',  15, 103}, $
    {strdefSG,'62f43297'XL,          '-I 2xc',       'I 2/c 1 1',         '',       'a3',  15, 104}, $
    {strdefSG,'37e5a351'XL,          '-C 2xc',       'C 2/c 1 1',         '',      '-a1',  15, 105}, $
    {strdefSG,'7d187140'XL,         '-B 2xbc',       'B 2/n 1 1',         '',      '-a2',  15, 106}, $
    {strdefSG,'fbdc56d1'XL,          '-I 2xb',       'I 2/b 1 1',         '',      '-a3',  15, 107}, $
    {strdefSG,'a959fc0b'XL,           'P 2 2',         'P 2 2 2',         '',         '',  16, 108}, $
    {strdefSG,'03d2cde6'XL,          'P 2c 2',        'P 2 2 21',         '',         '',  17, 109}, $
    {strdefSG,'d756eb1b'XL,         'P 2a 2a',        'P 21 2 2',         '',      'cab',  17, 110}, $
    {strdefSG,'010e6701'XL,          'P 2 2b',        'P 2 21 2',         '',      'bca',  17, 111}, $
    {strdefSG,'2b106ff0'XL,         'P 2 2ab',       'P 21 21 2',         '',         '',  18, 112}, $
    {strdefSG,'ee8326bb'XL,         'P 2bc 2',       'P 2 21 21',         '',      'cab',  18, 113}, $
    {strdefSG,'82570fbb'XL,       'P 2ac 2ac',       'P 21 2 21',         '',      'bca',  18, 114}, $
    {strdefSG,'8f7a6eec'XL,       'P 2ac 2ab',      'P 21 21 21',         '',         '',  19, 115}, $
    {strdefSG,'d3dceae0'XL,          'C 2c 2',        'C 2 2 21',         '',         '',  20, 116}, $
    {strdefSG,'a3d855bc'XL,         'A 2a 2a',        'A 21 2 2',         '',      'cab',  20, 117}, $
    {strdefSG,'31f9a8c4'XL,          'B 2 2b',        'B 2 21 2',         '',      'bca',  20, 118}, $
    {strdefSG,'bd0e64fb'XL,           'C 2 2',         'C 2 2 2',         '',         '',  21, 119}, $
    {strdefSG,'3af031fa'XL,           'A 2 2',         'A 2 2 2',         '',      'cab',  21, 120}, $
    {strdefSG,'f750b6b6'XL,           'B 2 2',         'B 2 2 2',         '',      'bca',  21, 121}, $
    {strdefSG,'8a25cd68'XL,           'F 2 2',         'F 2 2 2',         '',         '',  22, 122}, $
    {strdefSG,'e8bcf561'XL,           'I 2 2',         'I 2 2 2',         '',         '',  23, 123}, $
    {strdefSG,'7ba265f9'XL,         'I 2b 2c',      'I 21 21 21',         '',         '',  24, 124}, $
    {strdefSG,'715a8a38'XL,          'P 2 -2',         'P m m 2',         '',         '',  25, 125}, $
    {strdefSG,'34467527'XL,          'P -2 2',         'P 2 m m',         '',      'cab',  25, 126}, $
    {strdefSG,'454292ce'XL,         'P -2 -2',         'P m 2 m',         '',      'bca',  25, 127}, $
    {strdefSG,'dbd1bbd5'XL,         'P 2c -2',        'P m c 21',         '',         '',  26, 128}, $
    {strdefSG,'8ed05f75'XL,        'P 2c -2c',        'P c m 21',         '',     'ba-c',  26, 129}, $
    {strdefSG,'1e587dd6'XL,        'P -2a 2a',        'P 21 m a',         '',      'cab',  26, 130}, $
    {strdefSG,'4a496237'XL,         'P -2 2a',        'P 21 a m',         '',     '-cba',  26, 131}, $
    {strdefSG,'ed1509c4'XL,        'P -2 -2b',        'P b 21 m',         '',      'bca',  26, 132}, $
    {strdefSG,'1735e925'XL,        'P -2b -2',        'P m 21 b',         '',     'a-cb',  26, 133}, $
    {strdefSG,'245b6e98'XL,         'P 2 -2c',         'P c c 2',         '',         '',  27, 134}, $
    {strdefSG,'60576ac6'XL,         'P -2a 2',         'P 2 a a',         '',      'cab',  27, 135}, $
    {strdefSG,'bf62722f'XL,       'P -2b -2b',         'P b 2 b',         '',      'bca',  27, 136}, $
    {strdefSG,'0f559d28'XL,         'P 2 -2a',         'P m a 2',         '',         '',  28, 137}, $
    {strdefSG,'8b7a6ad9'XL,         'P 2 -2b',         'P b m 2',         '',     'ba-c',  28, 138}, $
    {strdefSG,'66310ecc'XL,         'P -2b 2',         'P 2 m b',         '',      'cab',  28, 139}, $
    {strdefSG,'9ecd44ca'XL,         'P -2c 2',         'P 2 c m',         '',     '-cba',  28, 140}, $
    {strdefSG,'1043766e'XL,       'P -2c -2c',         'P c 2 m',         '',      'bca',  28, 141}, $
    {strdefSG,'3b4d85de'XL,       'P -2a -2a',         'P m 2 a',         '',     'a-cb',  28, 142}, $
    {strdefSG,'f0df4865'XL,       'P 2c -2ac',        'P c a 21',         '',         '',  29, 143}, $
    {strdefSG,'c427e492'XL,        'P 2c -2b',        'P b c 21',         '',     'ba-c',  29, 144}, $
    {strdefSG,'183e19dc'XL,        'P -2b 2a',        'P 21 a b',         '',      'cab',  29, 145}, $
    {strdefSG,'b4d34c3b'XL,       'P -2ac 2a',        'P 21 c a',         '',     '-cba',  29, 146}, $
    {strdefSG,'a7e2b223'XL,      'P -2bc -2c',        'P c 21 b',         '',      'bca',  29, 147}, $
    {strdefSG,'931a1ed4'XL,      'P -2a -2ab',        'P b 21 a',         '',     'a-cb',  29, 148}, $
    {strdefSG,'358dd654'XL,        'P 2 -2bc',         'P n c 2',         '',         '',  30, 149}, $
    {strdefSG,'5a547988'XL,        'P 2 -2ac',         'P c n 2',         '',     'ba-c',  30, 150}, $
    {strdefSG,'cadc5b2b'XL,        'P -2ac 2',         'P 2 n a',         '',      'cab',  30, 151}, $
    {strdefSG,'3220112d'XL,        'P -2ab 2',         'P 2 a n',         '',     '-cba',  30, 152}, $
    {strdefSG,'c16d653f'XL,     'P -2ab -2ab',         'P b 2 n',         '',      'bca',  30, 153}, $
    {strdefSG,'0195cea2'XL,     'P -2bc -2bc',         'P n 2 b',         '',     'a-cb',  30, 154}, $
    {strdefSG,'8fc0a434'XL,        'P 2ac -2',        'P m n 21',         '',         '',  31, 155}, $
    {strdefSG,'72570ce4'XL,      'P 2bc -2bc',        'P n m 21',         '',     'ba-c',  31, 156}, $
    {strdefSG,'b60fe6dc'XL,      'P -2ab 2ab',        'P 21 m n',         '',      'cab',  31, 157}, $
    {strdefSG,'1f488697'XL,        'P -2 2ac',        'P 21 n m',         '',     '-cba',  31, 158}, $
    {strdefSG,'464f1412'XL,       'P -2 -2bc',        'P n 21 m',         '',      'bca',  31, 159}, $
    {strdefSG,'4324f6c4'XL,       'P -2ab -2',        'P m 21 n',         '',     'a-cb',  31, 160}, $
    {strdefSG,'f5757dc9'XL,        'P 2 -2ab',         'P b a 2',         '',         '',  32, 161}, $
    {strdefSG,'739caf97'XL,        'P -2bc 2',         'P 2 c b',         '',      'cab',  32, 162}, $
    {strdefSG,'6e4c617e'XL,     'P -2ac -2ac',         'P c 2 a',         '',      'bca',  32, 163}, $
    {strdefSG,'04df4f0f'XL,        'P 2c -2n',        'P n a 21',         '',         '',  33, 164}, $
    {strdefSG,'ba28f382'XL,       'P 2c -2ab',        'P b n 21',         '',     'ba-c',  33, 165}, $
    {strdefSG,'0d93b887'XL,       'P -2bc 2a',        'P 21 n b',         '',      'cab',  33, 166}, $
    {strdefSG,'5982a766'XL,        'P -2n 2a',        'P 21 c n',         '',     '-cba',  33, 167}, $
    {strdefSG,'d9eda533'XL,      'P -2n -2ac',        'P c 21 n',         '',      'bca',  33, 168}, $
    {strdefSG,'c83b1dff'XL,      'P -2ac -2n',        'P n 21 a',         '',     'a-cb',  33, 169}, $
    {strdefSG,'4b82c144'XL,         'P 2 -2n',         'P n n 2',         '',         '',  34, 170}, $
    {strdefSG,'278db076'XL,         'P -2n 2',         'P 2 n n',         '',      'cab',  34, 171}, $
    {strdefSG,'7f9ad9b2'XL,       'P -2n -2n',         'P n 2 n',         '',      'bca',  34, 172}, $
    {strdefSG,'93d6f19f'XL,          'C 2 -2',         'C m m 2',         '',         '',  35, 173}, $
    {strdefSG,'eea975ee'XL,          'A -2 2',         'A 2 m m',         '',      'cab',  35, 174}, $
    {strdefSG,'b633cb54'XL,         'B -2 -2',         'B m 2 m',         '',      'bca',  35, 175}, $
    {strdefSG,'fd047f84'XL,         'C 2c -2',        'C m c 21',         '',         '',  36, 176}, $
    {strdefSG,'851ce235'XL,        'C 2c -2c',        'C c m 21',         '',     'ba-c',  36, 177}]
sgdata2=[$
    {strdefSG,'9c241cdd'XL,        'A -2a 2a',        'A 21 m a',         '',      'cab',  36, 178}, $
    {strdefSG,'778111a8'XL,         'A -2 2a',        'A 21 a m',         '',     '-cba',  36, 179}, $
    {strdefSG,'709ad526'XL,        'B -2 -2b',        'B b 21 m',         '',      'bca',  36, 180}, $
    {strdefSG,'91092c8d'XL,        'B -2b -2',        'B m 21 b',         '',     'a-cb',  36, 181}, $
    {strdefSG,'ebce6c2e'XL,         'C 2 -2c',         'C c c 2',         '',         '',  37, 182}, $
    {strdefSG,'050c789b'XL,         'A -2a 2',         'A 2 a a',         '',      'cab',  37, 183}, $
    {strdefSG,'57a032ff'XL,       'B -2b -2b',         'B b 2 b',         '',      'bca',  37, 184}, $
    {strdefSG,'1428a49e'XL,          'A 2 -2',         'A m m 2',         '',         '',  38, 185}, $
    {strdefSG,'d98823d2'XL,          'B 2 -2',         'B m m 2',         '',     'ba-c',  38, 186}, $
    {strdefSG,'2309f2a2'XL,          'B -2 2',         'B 2 m m',         '',      'cab',  38, 187}, $
    {strdefSG,'695720ef'XL,          'C -2 2',         'C 2 m m',         '',     '-cba',  38, 188}, $
    {strdefSG,'fc6d1919'XL,         'C -2 -2',         'C m 2 m',         '',      'bca',  38, 189}, $
    {strdefSG,'7b934c18'XL,         'A -2 -2',         'A m 2 m',         '',     'a-cb',  38, 190}, $
    {strdefSG,'f5bb5d35'XL,         'A 2 -2c',         'A b m 2',         '',         '',  39, 191}, $
    {strdefSG,'40a04794'XL,         'B 2 -2c',         'B m a 2',         '',     'ba-c',  39, 192}, $
    {strdefSG,'c8acffd7'XL,         'B -2c 2',         'B 2 c m',         '',      'cab',  39, 193}, $
    {strdefSG,'82f22d9a'XL,         'C -2b 2',         'C 2 m b',         '',     '-cba',  39, 194}, $
    {strdefSG,'65457d5f'XL,       'C -2b -2b',         'C m 2 a',         '',      'bca',  39, 195}, $
    {strdefSG,'9a00b5b3'XL,       'A -2c -2c',         'A c 2 m',         '',     'a-cb',  39, 196}, $
    {strdefSG,'8d00c0d8'XL,         'A 2 -2a',         'A m a 2',         '',         '',  40, 197}, $
    {strdefSG,'381bda79'XL,         'B 2 -2b',         'B b m 2',         '',     'ba-c',  40, 198}, $
    {strdefSG,'0433157b'XL,         'B -2b 2',         'B 2 m b',         '',      'cab',  40, 199}, $
    {strdefSG,'0785aef4'XL,         'C -2c 2',         'C 2 c m',         '',     '-cba',  40, 200}, $
    {strdefSG,'847584a8'XL,       'C -2c -2c',         'C c 2 m',         '',      'bca',  40, 201}, $
    {strdefSG,'e2bb285e'XL,       'A -2a -2a',         'A m 2 a',         '',     'a-cb',  40, 202}, $
    {strdefSG,'6c933973'XL,        'A 2 -2ac',         'A b a 2',         '',         '',  41, 203}, $
    {strdefSG,'a133be3f'XL,        'B 2 -2bc',         'B b a 2',         '',     'ba-c',  41, 204}, $
    {strdefSG,'ef96180e'XL,        'B -2bc 2',         'B 2 c b',         '',      'cab',  41, 205}, $
    {strdefSG,'ec20a381'XL,        'C -2bc 2',         'C 2 c b',         '',     '-cba',  41, 206}, $
    {strdefSG,'1d5de0ee'XL,     'C -2bc -2bc',         'C c 2 a',         '',      'bca',  41, 207}, $
    {strdefSG,'0328d1f5'XL,     'A -2ac -2ac',         'A c 2 a',         '',     'a-cb',  41, 208}, $
    {strdefSG,'18172b53'XL,          'F 2 -2',         'F m m 2',         '',         '',  42, 209}, $
    {strdefSG,'ad1a17db'XL,          'F -2 2',         'F 2 m m',         '',      'cab',  42, 210}, $
    {strdefSG,'a88fc3e9'XL,         'F -2 -2',         'F m 2 m',         '',      'bca',  42, 211}, $
    {strdefSG,'921865a5'XL,         'F 2 -2d',         'F d d 2',         '',         '',  43, 212}, $
    {strdefSG,'1b59dbd3'XL,         'F -2d 2',         'F 2 d d',         '',      'cab',  43, 213}, $
    {strdefSG,'22808d1f'XL,       'F -2d -2d',         'F d 2 d',         '',      'bca',  43, 214}, $
    {strdefSG,'c6646005'XL,          'I 2 -2',         'I m m 2',         '',         '',  44, 215}, $
    {strdefSG,'3ce5b175'XL,          'I -2 2',         'I 2 m m',         '',      'cab',  44, 216}, $
    {strdefSG,'a9df8883'XL,         'I -2 -2',         'I m 2 m',         '',      'bca',  44, 217}, $
    {strdefSG,'bedffde8'XL,         'I 2 -2c',         'I b a 2',         '',         '',  45, 218}, $
    {strdefSG,'d740bc00'XL,         'I -2a 2',         'I 2 c b',         '',      'cab',  45, 219}, $
    {strdefSG,'484c7128'XL,       'I -2b -2b',         'I c 2 a',         '',      'bca',  45, 220}, $
    {strdefSG,'5f4c0443'XL,         'I 2 -2a',         'I m a 2',         '',         '',  46, 221}, $
    {strdefSG,'27f799ae'XL,         'I 2 -2b',         'I b m 2',         '',     'ba-c',  46, 222}, $
    {strdefSG,'1bdf56ac'XL,         'I -2b 2',         'I 2 m b',         '',      'cab',  46, 223}, $
    {strdefSG,'f07a5bd9'XL,         'I -2c 2',         'I 2 c m',         '',     '-cba',  46, 224}, $
    {strdefSG,'d164156e'XL,       'I -2c -2c',         'I c 2 m',         '',      'bca',  46, 225}, $
    {strdefSG,'30f7ecc5'XL,       'I -2a -2a',         'I m 2 a',         '',     'a-cb',  46, 226}, $
    {strdefSG,'e7243bbc'XL,          '-P 2 2',         'P m m m',         '',         '',  47, 227}, $
    {strdefSG,'51e85b6c'XL,       'P 2 2 -1n',         'P n n n',        '1',        '1',  48, 228}, $
    {strdefSG,'78ff5b81'XL,      '-P 2ab 2bc',         'P n n n',        '2',        '2',  48, 229}, $
    {strdefSG,'353736e5'XL,         '-P 2 2c',         'P c c m',         '',         '',  49, 230}, $
    {strdefSG,'0c8136c9'XL,         '-P 2a 2',         'P m a a',         '',      'cab',  49, 231}, $
    {strdefSG,'02b17898'XL,        '-P 2b 2b',         'P b m b',         '',      'bca',  49, 232}, $
    {strdefSG,'1392524e'XL,      'P 2 2 -1ab',         'P b a n',        '1',        '1',  50, 233}, $
    {strdefSG,'e91475ed'XL,       '-P 2ab 2b',         'P b a n',        '2',        '2',  50, 234}, $
    {strdefSG,'4e9de450'XL,      'P 2 2 -1bc',         'P n c b',        '1',     '1cab',  50, 235}, $
    {strdefSG,'935a56f4'XL,       '-P 2b 2bc',         'P n c b',        '2',     '2cab',  50, 236}, $
    {strdefSG,'38bbae1f'XL,      'P 2 2 -1ac',         'P c n a',        '1',     '1bca',  50, 237}, $
    {strdefSG,'de923b90'XL,        '-P 2a 2c',         'P c n a',        '2',     '2bca',  50, 238}, $
    {strdefSG,'6179e0c6'XL,        '-P 2a 2a',         'P m m a',         '',         '',  51, 239}, $
    {strdefSG,'2be88448'XL,         '-P 2b 2',         'P m m b',         '',     'ba-c',  51, 240}, $
    {strdefSG,'ce7dc76c'XL,         '-P 2 2b',         'P b m m',         '',      'cab',  51, 241}, $
    {strdefSG,'b013e0d3'XL,        '-P 2c 2c',         'P c m m',         '',     '-cba',  51, 242}, $
    {strdefSG,'6200ed8a'XL,         '-P 2c 2',         'P m c m',         '',      'bca',  51, 243}, $
    {strdefSG,'8adcedb3'XL,         '-P 2 2a',         'P m a m',         '',     'a-cb',  51, 244}, $
    {strdefSG,'c514e0e2'XL,       '-P 2a 2bc',         'P n n a',         '',         '',  52, 245}, $
    {strdefSG,'fea280fb'XL,        '-P 2b 2n',         'P n n b',         '',     'ba-c',  52, 246}, $
    {strdefSG,'38e06b40'XL,        '-P 2n 2b',         'P b n n',         '',      'cab',  52, 247}, $
    {strdefSG,'637980f3'XL,       '-P 2ab 2c',         'P c n n',         '',     '-cba',  52, 248}, $
    {strdefSG,'15078d8e'XL,       '-P 2ab 2n',         'P n c n',         '',      'bca',  52, 249}, $
    {strdefSG,'a90b452c'XL,       '-P 2n 2bc',         'P n a n',         '',     'a-cb',  52, 250}, $
    {strdefSG,'89a5e0ff'XL,        '-P 2ac 2',         'P m n a',         '',         '',  53, 251}, $
    {strdefSG,'42ae4859'XL,      '-P 2bc 2bc',         'P n m b',         '',     'ba-c',  53, 252}, $
    {strdefSG,'84eca3e2'XL,      '-P 2ab 2ab',         'P b m n',         '',      'cab',  53, 253}, $
    {strdefSG,'58cfe0ea'XL,        '-P 2 2ac',         'P c n m',         '',     '-cba',  53, 254}, $
    {strdefSG,'2eb1ed97'XL,        '-P 2 2bc',         'P n c m',         '',      'bca',  53, 255}, $
    {strdefSG,'c04d893d'XL,        '-P 2ab 2',         'P m a n',         '',     'a-cb',  53, 256}, $
    {strdefSG,'b36aed9f'XL,       '-P 2a 2ac',         'P c c a',         '',         '',  54, 257}, $
    {strdefSG,'88dc8d86'XL,        '-P 2b 2c',         'P c c b',         '',     'ba-c',  54, 258}, $
    {strdefSG,'25d8ca19'XL,        '-P 2a 2b',         'P b a a',         '',      'cab',  54, 259}, $
    {strdefSG,'5bb6eda6'XL,       '-P 2ac 2c',         'P c a a',         '',     '-cba',  54, 260}, $
    {strdefSG,'d3456635'XL,       '-P 2bc 2b',         'P b c b',         '',      'bca',  54, 261}, $
    {strdefSG,'6f49ae97'XL,       '-P 2b 2ab',         'P b a b',         '',     'a-cb',  54, 262}, $
    {strdefSG,'a3851163'XL,        '-P 2 2ab',         'P b a m',         '',         '',  55, 263}, $
    {strdefSG,'8b3b9e72'XL,        '-P 2bc 2',         'P m c b',         '',      'cab',  55, 264}, $
    {strdefSG,'364e3ba9'XL,      '-P 2ac 2ac',         'P c m a',         '',      'bca',  55, 265}, $
    {strdefSG,'0e8156fc'XL,      '-P 2ab 2ac',         'P c c n',         '',         '',  56, 266}]
sgdata3=[$
    {strdefSG,'31173243'XL,      '-P 2ac 2bc',         'P n a a',         '',      'cab',  56, 267}, $
    {strdefSG,'bebdb03a'XL,      '-P 2bc 2ab',         'P b n b',         '',      'bca',  56, 268}, $
    {strdefSG,'3a7e15cd'XL,        '-P 2c 2b',         'P b c m',         '',         '',  57, 269}, $
    {strdefSG,'ddeb36dc'XL,       '-P 2c 2ac',         'P c a m',         '',     'ba-c',  57, 270}, $
    {strdefSG,'e45d36f0'XL,       '-P 2ac 2a',         'P m c a',         '',      'cab',  57, 271}, $
    {strdefSG,'46105247'XL,        '-P 2b 2a',         'P m a b',         '',     '-cba',  57, 272}, $
    {strdefSG,'48201c16'XL,       '-P 2a 2ab',         'P b m a',         '',      'bca',  57, 273}, $
    {strdefSG,'280f97bc'XL,       '-P 2bc 2c',         'P c m b',         '',     'a-cb',  57, 274}, $
    {strdefSG,'43493b98'XL,         '-P 2 2n',         'P n n m',         '',         '',  58, 275}, $
    {strdefSG,'609e9307'XL,         '-P 2n 2',         'P m n n',         '',      'cab',  58, 276}, $
    {strdefSG,'c4f39323'XL,        '-P 2n 2n',         'P n m n',         '',      'bca',  58, 277}, $
    {strdefSG,'57337891'XL,    'P 2 2ab -1ab',         'P m m n',        '1',        '1',  59, 278}, $
    {strdefSG,'adb55f32'XL,       '-P 2ab 2a',         'P m m n',        '2',        '2',  59, 279}, $
    {strdefSG,'2282419e'XL,    'P 2bc 2 -1bc',         'P n m m',        '1',     '1cab',  59, 280}, $
    {strdefSG,'dab23f36'XL,       '-P 2c 2bc',         'P n m m',        '2',     '2cab',  59, 281}, $
    {strdefSG,'e9d1ae0a'XL,  'P 2ac 2ac -1ac',         'P m n m',        '1',     '1bca',  59, 282}, $
    {strdefSG,'0ff83b85'XL,        '-P 2c 2a',         'P m n m',        '2',     '2bca',  59, 283}, $
    {strdefSG,'5518bd4f'XL,       '-P 2n 2ab',         'P b c n',         '',         '',  60, 284}, $
    {strdefSG,'c3aa9ac9'XL,        '-P 2n 2c',         'P c a n',         '',     'ba-c',  60, 285}, $
    {strdefSG,'a8ec36ed'XL,        '-P 2a 2n',         'P n c a',         '',      'cab',  60, 286}, $
    {strdefSG,'2f569e56'XL,       '-P 2bc 2n',         'P n a b',         '',     '-cba',  60, 287}, $
    {strdefSG,'d1db18b8'XL,       '-P 2ac 2b',         'P b n a',         '',      'bca',  60, 288}, $
    {strdefSG,'e5245b89'XL,       '-P 2b 2ac',         'P c n b',         '',     'a-cb',  60, 289}, $
    {strdefSG,'bc23ceb7'XL,      '-P 2ac 2ab',         'P b c a',         '',         '',  61, 290}, $
    {strdefSG,'45f741b3'XL,      '-P 2bc 2ac',         'P c a b',         '',     'ba-c',  61, 291}, $
    {strdefSG,'5cefe44c'XL,       '-P 2ac 2n',         'P n m a',         '',         '',  62, 292}, $
    {strdefSG,'e6c3487d'XL,       '-P 2bc 2a',         'P m n b',         '',     'ba-c',  62, 293}, $
    {strdefSG,'5786c3c2'XL,       '-P 2c 2ab',         'P b n m',         '',      'cab',  62, 294}, $
    {strdefSG,'ae524cc6'XL,       '-P 2n 2ac',         'P c m n',         '',     '-cba',  62, 295}, $
    {strdefSG,'0d664508'XL,        '-P 2n 2a',         'P m c n',         '',      'bca',  62, 296}, $
    {strdefSG,'b74ae939'XL,        '-P 2c 2n',         'P n a m',         '',     'a-cb',  62, 297}, $
    {strdefSG,'c0dd6737'XL,         '-C 2c 2',         'C m c m',         '',         '',  63, 298}, $
    {strdefSG,'a4e78801'XL,        '-C 2c 2c',         'C c m m',         '',     'ba-c',  63, 299}, $
    {strdefSG,'637ac153'XL,        '-A 2a 2a',         'A m m a',         '',      'cab',  63, 300}, $
    {strdefSG,'7af93f9c'XL,         '-A 2 2a',         'A m a m',         '',     '-cba',  63, 301}, $
    {strdefSG,'5ea5920d'XL,         '-B 2 2b',         'B b m m',         '',      'bca',  63, 302}, $
    {strdefSG,'df056b84'XL,         '-B 2b 2',         'B m m b',         '',     'a-cb',  63, 303}, $
    {strdefSG,'d95e99f8'XL,        '-C 2bc 2',         'C m c a',         '',         '',  64, 304}, $
    {strdefSG,'014f7f4e'XL,      '-C 2bc 2bc',         'C c m b',         '',     'ba-c',  64, 305}, $
    {strdefSG,'e2da38da'XL,      '-A 2ac 2ac',         'A b m a',         '',      'cab',  64, 306}, $
    {strdefSG,'3e8a450b'XL,        '-A 2 2ac',         'A c a m',         '',     '-cba',  64, 307}, $
    {strdefSG,'e28e9b8d'XL,        '-B 2 2bc',         'B b c m',         '',      'bca',  64, 308}, $
    {strdefSG,'c686954b'XL,        '-B 2bc 2',         'B m a b',         '',     'a-cb',  64, 309}, $
    {strdefSG,'cc263714'XL,          '-C 2 2',         'C m m m',         '',         '',  65, 310}, $
    {strdefSG,'c6d2361c'XL,          '-A 2 2',         'A m m m',         '',      'cab',  65, 311}, $
    {strdefSG,'1ad6e89a'XL,          '-B 2 2',         'B m m m',         '',      'bca',  65, 312}, $
    {strdefSG,'a81cd822'XL,         '-C 2 2c',         'C c c m',         '',         '',  66, 313}, $
    {strdefSG,'df51c8d3'XL,         '-A 2a 2',         'A m a a',         '',      'cab',  66, 314}, $
    {strdefSG,'9b761113'XL,        '-B 2b 2b',         'B b m b',         '',      'bca',  66, 315}, $
    {strdefSG,'d5a5c9db'XL,         '-C 2b 2',         'C m m a',         '',         '',  67, 316}, $
    {strdefSG,'698ec05b'XL,        '-C 2b 2b',         'C m m b',         '',     'ba-c',  67, 317}, $
    {strdefSG,'4772cf95'XL,        '-A 2c 2c',         'A b m m',         '',      'cab',  67, 318}, $
    {strdefSG,'82a14c8b'XL,         '-A 2 2c',         'A c m m',         '',     '-cba',  67, 319}, $
    {strdefSG,'a6fde11a'XL,         '-B 2 2c',         'B m c m',         '',      'bca',  67, 320}, $
    {strdefSG,'03551655'XL,         '-B 2c 2',         'B m a m',         '',     'a-cb',  67, 321}, $
    {strdefSG,'d12ab103'XL,      'C 2 2 -1bc',         'C c c a',        '1',        '1',  68, 322}, $
    {strdefSG,'0db42f6d'XL,       '-C 2b 2bc',         'C c c a',        '2',        '2',  68, 323}, $
    {strdefSG,'d12ab103'XL,      'C 2 2 -1bc',         'C c c b',        '1',    '1ba-c',  68, 324}, $
    {strdefSG,'b19f26ed'XL,        '-C 2b 2c',         'C c c b',        '2',    '2ba-c',  68, 325}, $
    {strdefSG,'6a64cf1d'XL,      'A 2 2 -1ac',         'A b a a',        '1',     '1cab',  68, 326}, $
    {strdefSG,'9b22b244'XL,        '-A 2a 2c',         'A b a a',        '2',     '2cab',  68, 327}, $
    {strdefSG,'6a64cf1d'XL,      'A 2 2 -1ac',         'A c a a',        '1',    '1-cba',  68, 328}, $
    {strdefSG,'5ef1315a'XL,       '-A 2ac 2c',         'A c a a',        '2',    '2-cba',  68, 329}, $
    {strdefSG,'b660119b'XL,      'B 2 2 -1bc',         'B b c b',        '1',     '1bca',  68, 330}, $
    {strdefSG,'82f5efdc'XL,       '-B 2bc 2b',         'B b c b',        '2',     '2bca',  68, 331}, $
    {strdefSG,'b660119b'XL,      'B 2 2 -1bc',         'B b a b',        '1',    '1a-cb',  68, 332}, $
    {strdefSG,'275d1893'XL,       '-B 2b 2bc',         'B b a b',        '2',    '2a-cb',  68, 333}, $
    {strdefSG,'cd91580e'XL,          '-F 2 2',         'F m m m',         '',         '',  69, 334}, $
    {strdefSG,'d31c4512'XL,       'F 2 2 -1d',         'F d d d',        '1',        '1',  70, 335}, $
    {strdefSG,'c3fdc94f'XL,      '-F 2uv 2vw',         'F d d d',        '2',        '2',  70, 336}, $
    {strdefSG,'4a8ca5fc'XL,          '-I 2 2',         'I m m m',         '',         '',  71, 337}, $
    {strdefSG,'b2d4d6eb'XL,         '-I 2 2c',         'I b a m',         '',         '',  72, 338}, $
    {strdefSG,'530f5b33'XL,         '-I 2a 2',         'I m c b',         '',      'cab',  72, 339}, $
    {strdefSG,'cb2c5c75'XL,        '-I 2b 2b',         'I c m a',         '',      'bca',  72, 340}, $
    {strdefSG,'770755f5'XL,        '-I 2b 2c',         'I b c a',         '',         '',  73, 341}, $
    {strdefSG,'177c21a4'XL,        '-I 2a 2b',         'I c a b',         '',     'ba-c',  73, 342}, $
    {strdefSG,'8f5f26e2'XL,         '-I 2b 2',         'I m m a',         '',         '',  74, 343}, $
    {strdefSG,'ef2452b3'XL,        '-I 2a 2a',         'I m m b',         '',     'ba-c',  74, 344}, $
    {strdefSG,'6e84ab3a'XL,        '-I 2c 2c',         'I b m m',         '',      'cab',  74, 345}, $
    {strdefSG,'0effdf6b'XL,         '-I 2 2b',         'I c m m',         '',     '-cba',  74, 346}, $
    {strdefSG,'f6a7ac7c'XL,         '-I 2 2a',         'I m c m',         '',      'bca',  74, 347}, $
    {strdefSG,'96dcd82d'XL,         '-I 2c 2',         'I m a m',         '',     'a-cb',  74, 348}, $
    {strdefSG,'e194247e'XL,             'P 4',             'P 4',         '',         '',  75, 349}, $
    {strdefSG,'eb09f6d0'XL,            'P 4w',            'P 41',         '',         '',  76, 350}, $
    {strdefSG,'1e1ef133'XL,            'P 4c',            'P 42',         '',         '',  77, 351}, $
    {strdefSG,'21968a94'XL,           'P 4cw',            'P 43',         '',         '',  78, 352}, $
    {strdefSG,'8d2db265'XL,             'I 4',             'I 4',         '',         '',  79, 353}, $
    {strdefSG,'04c2f7b9'XL,           'I 4bw',            'I 41',         '',         '',  80, 354}, $
    {strdefSG,'5b1380ec'XL,            'P -4',            'P -4',         '',         '',  81, 355}]
sgdata4=[$
    {strdefSG,'429018fc'XL,            'I -4',            'I -4',         '',         '',  82, 356}, $
    {strdefSG,'7f473453'XL,            '-P 4',           'P 4/m',         '',         '',  83, 357}, $
    {strdefSG,'07df08e7'XL,           '-P 4c',          'P 42/m',         '',         '',  84, 358}, $
    {strdefSG,'8d9739ab'XL,      'P 4ab -1ab',           'P 4/n',        '1',        '1',  85, 359}, $
    {strdefSG,'28870b39'XL,           '-P 4a',           'P 4/n',        '2',        '2',  85, 360}, $
    {strdefSG,'da4091d2'XL,        'P 4n -1n',          'P 42/n',        '1',        '1',  86, 361}, $
    {strdefSG,'cfc232bb'XL,          '-P 4bc',          'P 42/n',        '2',        '2',  86, 362}, $
    {strdefSG,'e426791e'XL,            '-I 4',           'I 4/m',         '',         '',  87, 363}, $
    {strdefSG,'dd05fe66'XL,      'I 4bw -1bw',          'I 41/a',        '1',        '1',  88, 364}, $
    {strdefSG,'f5a09fa3'XL,          '-I 4ad',          'I 41/a',        '2',        '2',  88, 365}, $
    {strdefSG,'781566a5'XL,           'P 4 2',         'P 4 2 2',         '',         '',  89, 366}, $
    {strdefSG,'529dbafe'XL,       'P 4ab 2ab',        'P 42 1 2',         '',         '',  90, 367}, $
    {strdefSG,'1129373c'XL,         'P 4w 2c',        'P 41 2 2',         '',         '',  91, 368}, $
    {strdefSG,'7abbc4b0'XL,      'P 4abw 2nw',       'P 41 21 2',         '',         '',  92, 369}, $
    {strdefSG,'00a13651'XL,          'P 4c 2',        'P 42 2 2',         '',         '',  93, 370}, $
    {strdefSG,'b11a4eb4'XL,         'P 4n 2n',       'P 42 21 2',         '',         '',  94, 371}, $
    {strdefSG,'c507cd54'XL,        'P 4cw 2c',        'P 43 2 2',         '',         '',  95, 372}, $
    {strdefSG,'5af2c5d9'XL,      'P 4nw 2abw',       'P 43 21 2',         '',         '',  96, 373}, $
    {strdefSG,'cc0e6dad'XL,           'I 4 2',         'I 4 2 2',         '',         '',  97, 374}, $
    {strdefSG,'ac45c901'XL,       'I 4bw 2bw',        'I 41 2 2',         '',         '',  98, 375}, $
    {strdefSG,'93373957'XL,          'P 4 -2',         'P 4 m m',         '',         '',  99, 376}, $
    {strdefSG,'d7961388'XL,        'P 4 -2ab',         'P 4 b m',         '',         '', 100, 377}, $
    {strdefSG,'53dd13c8'XL,        'P 4c -2c',        'P 42 c m',         '',         '', 101, 378}, $
    {strdefSG,'25fb5987'XL,        'P 4n -2n',        'P 42 n m',         '',         '', 102, 379}, $
    {strdefSG,'4124340e'XL,         'P 4 -2c',         'P 4 c c',         '',         '', 103, 380}, $
    {strdefSG,'375a3973'XL,         'P 4 -2n',         'P 4 n c',         '',         '', 104, 381}, $
    {strdefSG,'81ce1e91'XL,         'P 4c -2',        'P 42 m c',         '',         '', 105, 382}, $
    {strdefSG,'5b8d5c98'XL,       'P 4c -2ab',        'P 42 b c',         '',         '', 106, 383}, $
    {strdefSG,'d27b6f7f'XL,          'I 4 -2',         'I 4 m m',         '',         '', 107, 384}, $
    {strdefSG,'2a231c68'XL,         'I 4 -2c',         'I 4 c m',         '',         '', 108, 385}, $
    {strdefSG,'6ee4a2b9'XL,        'I 4bw -2',        'I 41 m d',         '',         '', 109, 386}, $
    {strdefSG,'f0719502'XL,       'I 4bw -2c',        'I 41 c d',         '',         '', 110, 387}, $
    {strdefSG,'ae367c74'XL,          'P -4 2',        'P -4 2 m',         '',         '', 111, 388}, $
    {strdefSG,'c47b0b46'XL,         'P -4 2c',        'P -4 2 c',         '',         '', 112, 389}, $
    {strdefSG,'84bea02f'XL,        'P -4 2ab',       'P -4 21 m',         '',         '', 113, 390}, $
    {strdefSG,'67395465'XL,         'P -4 2n',       'P -4 21 c',         '',         '', 114, 391}, $
    {strdefSG,'d50d13d5'XL,         'P -4 -2',        'P -4 m 2',         '',         '', 115, 392}, $
    {strdefSG,'3a48de91'XL,        'P -4 -2c',        'P -4 c 2',         '',         '', 116, 393}, $
    {strdefSG,'00b4acac'XL,       'P -4 -2ab',        'P -4 b 2',         '',         '', 117, 394}, $
    {strdefSG,'f7b01eef'XL,        'P -4 -2n',        'P -4 n 2',         '',         '', 118, 395}, $
    {strdefSG,'7a127a20'XL,         'I -4 -2',        'I -4 m 2',         '',         '', 119, 396}, $
    {strdefSG,'c66b145c'XL,        'I -4 -2c',        'I -4 c 2',         '',         '', 120, 397}, $
    {strdefSG,'93e9bb9a'XL,          'I -4 2',        'I -4 2 m',         '',         '', 121, 398}, $
    {strdefSG,'f3a21f36'XL,        'I -4 2bw',        'I -4 2 d',         '',         '', 122, 399}, $
    {strdefSG,'b8f2113e'XL,          '-P 4 2',       'P 4/m m m',         '',         '', 123, 400}, $
    {strdefSG,'2de869d5'XL,         '-P 4 2c',       'P 4/m c c',         '',         '', 124, 401}, $
    {strdefSG,'b7fa689e'XL,      'P 4 2 -1ab',       'P 4/n b m',        '1',        '1', 125, 402}, $
    {strdefSG,'026e1f11'XL,        '-P 4a 2b',       'P 4/n b m',        '2',        '2', 125, 403}, $
    {strdefSG,'08a4eb7f'XL,       'P 4 2 -1n',       'P 4/n n c',        '1',        '1', 126, 404}, $
    {strdefSG,'51368ca3'XL,       '-P 4a 2bc',       'P 4/n n c',        '2',        '2', 126, 405}, $
    {strdefSG,'d5a84150'XL,        '-P 4 2ab',       'P 4/m b m',         '',         '', 127, 406}, $
    {strdefSG,'f1b4f061'XL,         '-P 4 2n',       'P 4/m n c',         '',         '', 128, 407}, $
    {strdefSG,'81d6f60a'XL,  'P 4ab 2ab -1ab',       'P 4/n m m',        '1',        '1', 129, 408}, $
    {strdefSG,'6f344f7f'XL,        '-P 4a 2a',       'P 4/n m m',        '2',        '2', 129, 409}, $
    {strdefSG,'9aee9967'XL,   'P 4ab 2n -1ab',       'P 4/n c c',        '1',        '1', 130, 410}, $
    {strdefSG,'8d6a1517'XL,       '-P 4a 2ac',       'P 4/n c c',        '2',        '2', 130, 411}, $
    {strdefSG,'5a9858b3'XL,         '-P 4c 2',      'P 42/m m c',         '',         '', 131, 412}, $
    {strdefSG,'cf822058'XL,        '-P 4c 2c',      'P 42/m c m',         '',         '', 132, 413}, $
    {strdefSG,'00e9903a'XL,     'P 4n 2c -1n',      'P 42/n b c',        '1',        '1', 133, 414}, $
    {strdefSG,'f9df795b'XL,       '-P 4ac 2b',      'P 42/n b c',        '2',        '2', 133, 415}, $
    {strdefSG,'aad7368d'XL,      'P 4n 2 -1n',      'P 42/n n m',        '1',        '1', 134, 416}, $
    {strdefSG,'aa87eae9'XL,      '-P 4ac 2bc',      'P 42/n n m',        '2',        '2', 134, 417}, $
    {strdefSG,'a5c2e426'XL,       '-P 4c 2ab',      'P 42/m b c',         '',         '', 135, 418}, $
    {strdefSG,'53c72d93'XL,        '-P 4n 2n',      'P 42/m n m',         '',         '', 136, 419}, $
    {strdefSG,'e391d7d2'XL,     'P 4n 2n -1n',      'P 42/n m c',        '1',        '1', 137, 420}, $
    {strdefSG,'0685c5ce'XL,       '-P 4ac 2a',      'P 42/n m c',        '2',        '2', 137, 421}, $
    {strdefSG,'f8a9b8bf'XL,    'P 4n 2ab -1n',      'P 42/n c m',        '1',        '1', 138, 422}, $
    {strdefSG,'e4db9fa6'XL,      '-P 4ac 2ac',      'P 42/n c m',        '2',        '2', 138, 423}, $
    {strdefSG,'7a639eef'XL,          '-I 4 2',       'I 4/m m m',         '',         '', 139, 424}, $
    {strdefSG,'d7e16f02'XL,         '-I 4 2c',       'I 4/m c m',         '',         '', 140, 425}, $
    {strdefSG,'e12277ea'XL,  'I 4bw 2bw -1bw',      'I 41/a m d',        '1',        '1', 141, 426}, $
    {strdefSG,'b96409e0'XL,        '-I 4bd 2',      'I 41/a m d',        '2',        '2', 141, 427}, $
    {strdefSG,'8bdd8359'XL,  'I 4bw 2aw -1bw',      'I 41/a c d',        '1',        '1', 142, 428}, $
    {strdefSG,'08d64e97'XL,       '-I 4bd 2c',      'I 41/a c d',        '2',        '2', 142, 429}, $
    {strdefSG,'329c980e'XL,             'P 3',             'P 3',         '',         '', 143, 430}, $
    {strdefSG,'86dd212b'XL,            'P 31',            'P 31',         '',         '', 144, 431}, $
    {strdefSG,'d70a5f46'XL,            'P 32',            'P 32',         '',         '', 145, 432}, $
    {strdefSG,'8e4e25f8'XL,             'R 3',             'R 3',        'H',        'H', 146, 433}, $
    {strdefSG,'d5a0aa2d'XL,            'P 3*',             'R 3',        'R',        'R', 146, 434}, $
    {strdefSG,'fdd759b5'XL,            '-P 3',            'P -3',         '',         '', 147, 435}, $
    {strdefSG,'be8d0d7f'XL,            '-R 3',            'R -3',        'H',        'H', 148, 436}, $
    {strdefSG,'d9a29bac'XL,           '-P 3*',            'R -3',        'R',        'R', 148, 437}, $
    {strdefSG,'65b7a72b'XL,           'P 3 2',         'P 3 1 2',         '',         '', 149, 438}, $
    {strdefSG,'c1840a7a'XL,          'P 3 2"',         'P 3 2 1',         '',         '', 150, 439}, $
    {strdefSG,'97e2dfd5'XL, 'P 31 2c (0 0 1)',        'P 31 1 2',         '',         '', 151, 440}, $
    {strdefSG,'33d17284'XL,         'P 31 2"',        'P 31 2 1',         '',         '', 152, 441}, $
    {strdefSG,'e39f36b4'XL,'P 32 2c (0 0 -1)',        'P 32 1 2',         '',         '', 153, 442}, $
    {strdefSG,'47ac9be5'XL,         'P 32 2"',        'P 32 2 1',         '',         '', 154, 443}, $
    {strdefSG,'46ebee09'XL,          'R 3 2"',            'R 32',        'H',        'H', 155, 444}]
sgdata5=[$
    {strdefSG,'a20b8591'XL,          'P 3* 2',            'R 32',        'R',        'R', 155, 445}, $
    {strdefSG,'9f4cffaa'XL,         'P 3 -2"',         'P 3 m 1',         '',         '', 156, 446}, $
    {strdefSG,'39859b12'XL,          'P 3 -2',         'P 3 1 m',         '',         '', 157, 447}, $
    {strdefSG,'e04fe588'XL,        'P 3 -2"c',         'P 3 c 1',         '',         '', 158, 448}, $
    {strdefSG,'ec0db0dd'XL,         'P 3 -2c',         'P 3 1 c',         '',         '', 159, 449}, $
    {strdefSG,'6e32558c'XL,         'R 3 -2"',           'R 3 m',        'H',        'H', 160, 450}, $
    {strdefSG,'b951b4f7'XL,         'P 3* -2',           'R 3 m',        'R',        'R', 160, 451}, $
    {strdefSG,'f7d1a830'XL,        'R 3 -2"c',           'R 3 c',        'H',        'H', 161, 452}, $
    {strdefSG,'219be015'XL,        'P 3* -2n',           'R 3 c',        'R',        'R', 161, 453}, $
    {strdefSG,'f74c7f83'XL,          '-P 3 2',        'P -3 1 m',         '',         '', 162, 454}, $
    {strdefSG,'69dc2f41'XL,         '-P 3 2c',        'P -3 1 c',         '',         '', 163, 455}, $
    {strdefSG,'fc3edafb'XL,         '-P 3 2"',        'P -3 m 1',         '',         '', 164, 456}, $
    {strdefSG,'e9d82a99'XL,        '-P 3 2"c',        'P -3 c 1',         '',         '', 165, 457}, $
    {strdefSG,'6df507a9'XL,         '-R 3 2"',          'R -3 m',        'H',        'H', 166, 458}, $
    {strdefSG,'1c80e47a'XL,         '-P 3* 2',          'R -3 m',        'R',        'R', 166, 459}, $
    {strdefSG,'9a7d09d3'XL,        '-R 3 2"c',          'R -3 c',        'H',        'H', 167, 460}, $
    {strdefSG,'bb691c91'XL,        '-P 3* 2n',          'R -3 c',        'R',        'R', 167, 461}, $
    {strdefSG,'a2ddaf47'XL,             'P 6',             'P 6',         '',         '', 168, 462}, $
    {strdefSG,'81a2968a'XL,            'P 61',            'P 61',         '',         '', 169, 463}, $
    {strdefSG,'84b83b8b'XL,            'P 65',            'P 65',         '',         '', 170, 464}, $
    {strdefSG,'62d97e5d'XL,            'P 62',            'P 62',         '',         '', 171, 465}, $
    {strdefSG,'b22fefa7'XL,            'P 64',            'P 64',         '',         '', 172, 466}, $
    {strdefSG,'dddeb565'XL,            'P 6c',            'P 63',         '',         '', 173, 467}, $
    {strdefSG,'35eb6ccb'XL,            'P -6',            'P -6',         '',         '', 174, 468}, $
    {strdefSG,'32dacfb6'XL,            '-P 6',           'P 6/m',         '',         '', 175, 469}, $
    {strdefSG,'06c1ae99'XL,           '-P 6c',          'P 63/m',         '',         '', 176, 470}, $
    {strdefSG,'90230e5f'XL,           'P 6 2',         'P 6 2 2',         '',         '', 177, 471}, $
    {strdefSG,'eed642a2'XL, 'P 61 2 (0 0 -1)',        'P 61 2 2',         '',         '', 178, 472}, $
    {strdefSG,'891aeb21'XL,  'P 65 2 (0 0 1)',        'P 65 2 2',         '',         '', 179, 473}, $
    {strdefSG,'6f55448f'XL, 'P 62 2c (0 0 1)',        'P 62 2 2',         '',         '', 180, 474}, $
    {strdefSG,'655e4cd9'XL,'P 64 2c (0 0 -1)',        'P 64 2 2',         '',         '', 181, 475}, $
    {strdefSG,'a4386f70'XL,         'P 6c 2c',        'P 63 2 2',         '',         '', 182, 476}, $
    {strdefSG,'3b0a2d17'XL,          'P 6 -2',         'P 6 m m',         '',         '', 183, 477}, $
    {strdefSG,'b74f3653'XL,         'P 6 -2c',         'P 6 c c',         '',         '', 184, 478}, $
    {strdefSG,'f4d8e94d'XL,         'P 6c -2',        'P 63 c m',         '',         '', 185, 479}, $
    {strdefSG,'789df209'XL,        'P 6c -2c',        'P 63 m c',         '',         '', 186, 480}, $
    {strdefSG,'eeacb736'XL,          'P -6 2',        'P -6 m 2',         '',         '', 187, 481}, $
    {strdefSG,'fb4a4754'XL,         'P -6c 2',        'P -6 c 2',         '',         '', 188, 482}, $
    {strdefSG,'55b5be6a'XL,         'P -6 -2',        'P -6 2 m',         '',         '', 189, 483}, $
    {strdefSG,'cb25eea8'XL,       'P -6c -2c',        'P -6 2 c',         '',         '', 190, 484}, $
    {strdefSG,'f1fc7952'XL,          '-P 6 2',       'P 6/m m m',         '',         '', 191, 485}, $
    {strdefSG,'87f9e8ca'XL,         '-P 6 2c',       'P 6/m c c',         '',         '', 192, 486}, $
    {strdefSG,'1eb41bd9'XL,         '-P 6c 2',      'P 63/m c m',         '',         '', 193, 487}, $
    {strdefSG,'68b18a41'XL,        '-P 6c 2c',      'P 63/m m c',         '',         '', 194, 488}, $
    {strdefSG,'5843870d'XL,         'P 2 2 3',           'P 2 3',         '',         '', 195, 489}, $
    {strdefSG,'93e38c71'XL,         'F 2 2 3',           'F 2 3',         '',         '', 196, 490}, $
    {strdefSG,'dc4003c1'XL,         'I 2 2 3',           'I 2 3',         '',         '', 197, 491}, $
    {strdefSG,'f9e6a645'XL,     'P 2ac 2ab 3',          'P 21 3',         '',         '', 198, 492}, $
    {strdefSG,'7ec4457b'XL,       'I 2b 2c 3',          'I 21 3',         '',         '', 199, 493}, $
    {strdefSG,'72e55913'XL,        '-P 2 2 3',          'P m -3',         '',         '', 200, 494}, $
    {strdefSG,'265ce726'XL,     'P 2 2 3 -1n',          'P n -3',        '1',        '1', 201, 495}, $
    {strdefSG,'d419d8a3'XL,    '-P 2ab 2bc 3',          'P n -3',        '2',        '2', 201, 496}, $
    {strdefSG,'58320b8d'XL,        '-F 2 2 3',          'F m -3',         '',         '', 202, 497}, $
    {strdefSG,'7de7c89b'XL,     'F 2 2 3 -1d',          'F d -3',        '1',        '1', 203, 498}, $
    {strdefSG,'159ca8d3'XL,    '-F 2uv 2vw 3',          'F d -3',        '2',        '2', 203, 499}, $
    {strdefSG,'e23893df'XL,        '-I 2 2 3',          'I m -3',         '',         '', 204, 500}, $
    {strdefSG,'1d5f4d3f'XL,    '-P 2ac 2ab 3',          'P a -3',         '',         '', 205, 501}, $
    {strdefSG,'7ce42b66'XL,      '-I 2b 2c 3',          'I a -3',         '',         '', 206, 502}, $
    {strdefSG,'93a6edeb'XL,         'P 4 2 3',         'P 4 3 2',         '',         '', 207, 503}, $
    {strdefSG,'c55ae72a'XL,        'P 4n 2 3',        'P 42 3 2',         '',         '', 208, 504}, $
    {strdefSG,'25d65cf1'XL,         'F 4 2 3',         'F 4 3 2',         '',         '', 209, 505}, $
    {strdefSG,'ddc15ef3'XL,        'F 4d 2 3',        'F 41 3 2',         '',         '', 210, 506}, $
    {strdefSG,'014b7ed2'XL,         'I 4 2 3',         'I 4 3 2',         '',         '', 211, 507}, $
    {strdefSG,'4bca2b34'XL,    'P 4acd 2ab 3',        'P 43 3 2',         '',         '', 212, 508}, $
    {strdefSG,'0e761d2d'XL,     'P 4bd 2ab 3',        'P 41 3 2',         '',         '', 213, 509}, $
    {strdefSG,'bb92b652'XL,      'I 4bd 2c 3',        'I 41 3 2',         '',         '', 214, 510}, $
    {strdefSG,'6fc31c7b'XL,        'P -4 2 3',        'P -4 3 m',         '',         '', 215, 511}, $
    {strdefSG,'09abdebe'XL,        'F -4 2 3',        'F -4 3 m',         '',         '', 216, 512}, $
    {strdefSG,'75d77c69'XL,        'I -4 2 3',        'I -4 3 m',         '',         '', 217, 513}, $
    {strdefSG,'393f16ba'XL,       'P -4n 2 3',        'P -4 3 n',         '',         '', 218, 514}, $
    {strdefSG,'0968710f'XL,       'F -4c 2 3',        'F -4 3 c',         '',         '', 219, 515}, $
    {strdefSG,'e2ed982e'XL,     'I -4bd 2c 3',        'I -4 3 d',         '',         '', 220, 516}, $
    {strdefSG,'74c407d3'XL,        '-P 4 2 3',        'P m -3 m',         '',         '', 221, 517}, $
    {strdefSG,'c5e3dd9f'XL,     'P 4 2 3 -1n',        'P n -3 n',        '1',        '1', 222, 518}, $
    {strdefSG,'c7e69ab6'XL,     '-P 4a 2bc 3',        'P n -3 n',        '2',        '2', 222, 519}, $
    {strdefSG,'babd71c4'XL,       '-P 4n 2 3',        'P m -3 n',         '',         '', 223, 520}, $
    {strdefSG,'0b9aab88'XL,    'P 4n 2 3 -1n',        'P n -3 m',        '1',        '1', 224, 521}, $
    {strdefSG,'dbdb6bfc'XL,    '-P 4bc 2bc 3',        'P n -3 m',        '2',        '2', 224, 522}, $
    {strdefSG,'5225b614'XL,        '-F 4 2 3',        'F m -3 m',         '',         '', 225, 523}, $
    {strdefSG,'481e9f10'XL,       '-F 4c 2 3',        'F m -3 c',         '',         '', 226, 524}, $
    {strdefSG,'023a8184'XL,    'F 4d 2 3 -1d',        'F d -3 m',        '1',        '1', 227, 525}, $
    {strdefSG,'38f5e458'XL,    '-F 4vw 2vw 3',        'F d -3 m',        '2',        '2', 227, 526}, $
    {strdefSG,'82b6e512'XL,   'F 4d 2 3 -1cd',        'F d -3 c',        '1',        '1', 228, 527}, $
    {strdefSG,'4d325ad3'XL,   '-F 4cvw 2vw 3',        'F d -3 c',        '2',        '2', 228, 528}, $
    {strdefSG,'52045e84'XL,        '-I 4 2 3',        'I m -3 m',         '',         '', 229, 529}, $
    {strdefSG,'407de7c1'XL,     '-I 4bd 2c 3',        'I a -3 d',         '',         '', 230, 530}]
sgdata=[ temporary(sgdata0) ,$
temporary(sgdata1) ,$
temporary(sgdata2) ,$
temporary(sgdata3) ,$
temporary(sgdata4) ,$
temporary(sgdata5)]
n=530

; Type:
;    0: Hash in, everything out
;    1: HM name in, everything out
;    2: IT number in, everything out
;    3: Tabulated setting number in, everything out

if not keyword_set(pref_12) then pref_12='1'
if not keyword_set(pref_hr) then pref_hr='H'

case type of
0:    begin
    i=0
    ind=where(sgdata[*].sghash eq in[0],ct)
    case ct of
    0:    return,tmp
    else:    return,sgdata[ind[0]]
    endcase
    endcase
1:    begin
    i=0
    ind=where(sgdata[*].hm eq in[0],ct)
    case ct of
    0:    return,tmp
    1:    return,sgdata[ind[0]]
    else:begin
        ind2=where(sgdata[ind].ext eq pref_12 or sgdata[ind].ext eq pref_hr,ct)
        case ct of
        0:    return,sgdata[ind[0]]
        else:    return,sgdata[ind[ind2[0]]]
        endcase
        endelse
    endcase
    endcase
2:    begin
    i=0
    ind=where(sgdata[*].num eq in[0],ct)
    case ct of
    0:    return,tmp
    1:    return,sgdata[ind[0]]
    else:begin
        ind2=where(sgdata[ind].ext eq pref_12 or sgdata[ind].ext eq pref_hr,ct)
        case ct of
        0:    return,sgdata[ind[0]]
        else:    return,sgdata[ind[ind2[0]]]
        endcase
        endelse
    endcase
    endcase
3:    begin
    if (in le 0) or (in gt n) then return,tmp;No such Tabulated Setting number
    return,sgdata[in]
    endcase
endcase

end;function sgdata
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakeSeitz,rot,trn
return,[[rot,reform(trn,1,3)],[0,0,0,1.]]
end;function MakeSeitz
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MakeSeitzRat,rot,trn
type=size(rot[0].num,/type)

Seitz=identityrat(4,type=type)
Seitz[0:2,0:2]=rot
Seitz[3,0:2]=reform(trn,1,3)

return,Seitz
end;function MakeSeitzRat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rationalString,f, b, sign=sign

sign=keyword_set(sign)

n = round( abs( b * f ) );
d = b;
s=''
if (sign) then s = s + ( ( f gt 0 ) ? "+" : "-" ) $
else s = s + ( ( f gt 0 ) ? "" : "-" );

for i=5,2,-1 do begin
    if ( ( (n mod i) eq 0 ) and ( (d mod i) eq 0 ) ) then begin
        n = n/i;
        d = d/i
    endif
endfor

s += strcompress(string(n),/remove_all)
if ( d ne 1 ) then s += "/" + strcompress(string(d),/remove_all)

return, s

end;function rationalString
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rationalStringToFloat,s
; print,rationalStringToFloat(rationalString(0.25,24))

fs=reform(float(s))

; s is a string: 0.3333 or 1./3 => return 0.33333
i=strpos(s,'/')

ind=where(i ne -1,ct)
if ct eq 0 then return, fs

indd=(ct+1)*indgen(ct)
fs[ind]=float((strmid(s[ind],0,i[ind]))[indd])/float((strmid(s[ind],i[ind]+1))[indd])
return,fs
end;function rationalStringToFloat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rationalStringToRational,s,base,_extra=ex
; printrational,rationalStringToRational(rationalString(0.25,24),24)

fs=reform(float(s))
Float2Rational,fs,base,/thres,_extra=ex

; s is a string: 0.3333 or 1./3 => return 0.33333
i=strpos(s,'/')

ind=where(i ne -1,ct,comp=comp,ncomp=ncomp)
if ct eq 0 then return, fs

indd=(ct+1)*indgen(ct)
fs[ind].num=(strmid(s[ind],0,i[ind]))[indd]
fs[ind].denom=(strmid(s[ind],i[ind]+1))[indd]
SignNum,fs
ReduceNumDenom,fs

return,fs
end;function rationalStringToRational
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isdigit,c
; c is a byte or string
str=string(c)
return,str ge '0' and str le '9'
end;function isdigit
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StringRTop,op

xyz=["x","y","z"]
s=""

for i=0l,2 do begin
    t = "";
    for j=0l,2 do $
      if ( op.rot[j,i] ne 0.0 ) then begin
        t += (( op.rot[j,i] gt 0.0 ) ? "+" : "-");
        if ( abs( op.rot[j,i] )  ne 1 ) then $
          t += rationalString( abs( op.rot[j,i] ), 24);
        t += xyz[j];
      endif
    if ( op.trn[i] ne 0.0 ) then t += rationalString( op.trn[i], 24, /sign );
    if t eq '' then s+='0' $
    else s += strmid(t, ( strmid(t,0,1) eq '+' ) ? 1 : 0 );
    if ( i lt 2 ) then s+= ", ";
endfor

return, s;
end;function StringRTop
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StringRTopRat,op

xyz=["x","y","z"]
s=""

for i=0l,2 do begin
    t = "";
    for j=0l,2 do $
      if ( op.rot[j,i].num ne 0 ) then begin
        t += (( op.rot[j,i].num gt 0 ) ? "+" : "-");
        if ( abs( op.rot[j,i].num )  ne 1 or abs( op.rot[j,i].denom )  ne 1) then begin
            tmp=op.rot[j,i]
            tmp.num=abs(tmp.num)
            
            t+=stringr(tmp.num)
            if tmp.denom ne 1 then t+='/'+stringr(tmp.denom)

        endif
        t += xyz[j];
      endif
    if ( op.trn[i].num ne 0 ) then begin
        tmp=op.trn[i]
        t += (( tmp.num gt 0 ) ? "+" : "-");
        tmp.num=abs(tmp.num)
        
        t+=stringr(tmp.num)
        if tmp.denom ne 1 then t+='/'+stringr(tmp.denom)
    endif
    if t eq '' then s+='0' $
    else s += strmid(t, ( strmid(t,0,1) eq '+' ) ? 1 : 0 );
    if ( i lt 2 ) then s+= ", ";
endfor

return, s;
end;function StringRTopRat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RTopString,strop,error=error
;Construct an RT operator from a string description,
;  e.g. 1/2x,z-y+2/3,x
;  '*' is optional for multiplication, commas are compulsory.
error=1b

bmin=(byte('-'))[0] & bplus=(byte('+'))[0] & bsl=(byte('/'))[0]
bx=(byte('x'))[0] & bxx=(byte('X'))[0]
by=(byte('y'))[0] & byy=(byte('Y'))[0]
bz=(byte('z'))[0] & bzz=(byte('Z'))[0] & bp=(byte('.'))[0]

rows = str_sep(strcompress(strop,/remove_all),",");
if ( n_elements(rows) ne 3 ) then return,0b

rot=fltarr(3,3)
trn=fltarr(3)

for nrow=0,2 do begin
    row = byte(rows[nrow])
    npartmax=0
    cols=""
    for nchr = 0,n_elements(row)-1 do begin
      c = row[nchr]
      if ( c eq bplus or c eq bmin ) then begin
          npartmax=npartmax+1
        cols=[cols,""]
      endif
      cols[npartmax] = cols[npartmax]+string(c)
    endfor
    for npart = 0,npartmax do begin
      col = byte(cols[npart])
      ncol = 3;
      num = "";
      for nchr = 0,n_elements(col)-1 do begin
        c = col[nchr];
        if      ( c eq bx or c eq bxx ) then ncol = 0 $
        else if ( c eq by or c eq byy ) then ncol = 1 $
        else if ( c eq bz or c eq bzz ) then ncol = 2 $
        else if ( c eq bmin or c eq bsl or c eq bp or isdigit(c) ) then num = num+string(c);
      endfor
      if ( ncol lt 3 ) then begin
        if ( num eq "" or num eq "+" or num eq "-" ) then num = num+'1'
        rot[ncol,nrow] = rationalStringToFloat(num)
      endif else begin
        if ( num eq "" or num eq "+" or num eq "-" ) then num = num+'0'
        trn[ nrow ] = rationalStringToFloat(num)
      endelse
    endfor
endfor

error=0b
return,{rot:rot,trn:trn}
end;function RTopString
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RTopStringAll,strops,error=error

error=1b
ret={rot:fltarr(3,3),trn:fltarr(3)}

symops = strsplit(strops,';',/extract,count=n)
if n eq 0 then return,ret

for i=0l,n-1 do begin
    tmp=RTopString(symops[i],error=error)
    if error then return,ret
    ret=[ret,tmp]
endfor

ret=ret[1:*]
error=0b
return,ret
end;function RTopStringAll
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FindMaxDenomStrops,strops

if strpos(strops,'.') ne -1 then return,24

maxdenom=1
rat=strsplit(strops,'[],; xyzXYZ+-*',/extract,count=n)
for i=0l,n-1 do begin
    rati=strsplit(rat[i],'/',/extract,count=ni)
    if ni eq 2 then maxdenom>=long(rati[1])
endfor

return,maxdenom
end;function FindMaxDenomStrops
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RTopStringRat,strop,error=error,_extra=ex
;Construct an RT operator from a string description,
;  e.g. 1/2x,z-y+2/3,x
;  '*' is optional for multiplication, commas are compulsory.
error=1b

base=FindMaxDenomStrops(strop)

bmin=(byte('-'))[0] & bplus=(byte('+'))[0] & bsl=(byte('/'))[0]
bx=(byte('x'))[0] & bxx=(byte('X'))[0]
by=(byte('y'))[0] & byy=(byte('Y'))[0]
bz=(byte('z'))[0] & bzz=(byte('Z'))[0] & bp=(byte('.'))[0]

rows = str_sep(strcompress(strop,/remove_all),",");
if ( n_elements(rows) ne 3 ) then return,0b

rot=make_ratarray([3,3],_extra=ex)
trn=make_ratarray(3,_extra=ex)

for nrow=0,2 do begin
    row = byte(rows[nrow])
    npartmax=0
    cols=""
    for nchr = 0,n_elements(row)-1 do begin
      c = row[nchr]
      if ( c eq bplus or c eq bmin ) then begin
          npartmax=npartmax+1
        cols=[cols,""]
      endif
      cols[npartmax] = cols[npartmax]+string(c)
    endfor
    for npart = 0,npartmax do begin
      col = byte(cols[npart])
      ncol = 3;
      num = "";
      for nchr = 0,n_elements(col)-1 do begin
        c = col[nchr];
        if      ( c eq bx or c eq bxx ) then ncol = 0 $
        else if ( c eq by or c eq byy ) then ncol = 1 $
        else if ( c eq bz or c eq bzz ) then ncol = 2 $
        else if ( c eq bmin or c eq bsl or c eq bp or isdigit(c) ) then num = num+string(c);
      endfor
      if ( ncol lt 3 ) then begin
        if ( num eq "" or num eq "+" or num eq "-" ) then num = num+'1'
        rot[ncol,nrow] = rationalStringToRational(num,base,_extra=ex)
      endif else begin
        if ( num eq "" or num eq "+" or num eq "-" ) then num = num+'0'
        trn[ nrow ] = rationalStringToRational(num,base,_extra=ex)
      endelse
    endfor
endfor

error=0b
return,{rot:rot,trn:trn}
end;function RTopStringRat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RTopStringRatAll,strops,base,error=error,_extra=ex

error=1b
ret={rot:make_ratarray([3,3],_extra=ex),trn:make_ratarray(3,_extra=ex)}

symops = strsplit(strops,';',/extract,count=n)
if n eq 0 then return,ret

for i=0l,n-1 do begin
    tmp=RTopStringRat(symops[i],error=error,_extra=ex)
    if error then return,ret
    ret=[ret,tmp]
endfor

ret=ret[1:*]
error=0b
return,ret
end;function RTopStringRatAll
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Remain,a,b
; Same as prquot but faster (because handles only b>0)
c = a mod b
c+=(c lt 0)*b
return, c
end;function Remain
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TrnCode,code
return,code and '0000ffff'XL
end;function TrnCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TrnCode64,code
return,code and '00000000ffffffff'XLL
end;function TrnCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RotCode,code
return,code and 'ffff0000'XL
end;function RotCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RotCode64,code
return,code and 'ffffffff00000000'XLL
end;function TrnCode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CodeType64,t,r
return,{t:long64(t),r:long64(r)}
end;function CodeType64
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SortCodeType64,ops
subs=bsort(ops.t)
return,subs[bsort(ops[subs].r)]
end;function SortCodeType64
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Symop,code,seitz=seitz

code_t=TrnCode(code)
fac = [1L,24,576]
trn = Remain( code_t/fac, 24 )/24.

code_r = ishft( RotCode(code) , -16 ) xor '4064'XL
fac = [[1L,3,9],[27,81,243],[729,2187,6561]]
rot = Remain( code_r/fac, 3 ) - 1

if keyword_set(seitz) then return,MakeSeitz(rot,trn) $
else return,{rot:rot,trn:trn}
end;function Symop
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Symop64,code,seitz=seitz

code_t=TrnCode64(code.t)
fac = [1LL,24,576]
trn = Remain( code_t/fac, 24 )/24.

code_t2=ishft( RotCode64(code.t) , -32 )
fac=[1L,2,4,8,16,32,64,128,256]
rot=reform(Remain(code_t2/fac,2),3,3)*4 ; This will add 4 to all -2's that were originally 2's

code_r = code.r xor '7423B088F0ED248'XLL
fac=[[1LL,96,9216],[884736,84934656,8153726976],[782757789696,75144747810816,7213895789838336]]
rot += Remain( code_r/fac, 96 )/24.-2

if keyword_set(seitz) then return,MakeSeitz(rot,trn) $
else return,{rot:rot,trn:trn}
end;function Symop64
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SymopAll,allops,seitz=seitz
ngeneral=n_elements(allops)
if keyword_set(seitz) then begin
    ret=fltarr(4,4)
    for i=0l,ngeneral-1 do $
        ret=[[[ret]],[[Symop(allops[i],/seitz)]]]
    ret=ret[*,*,1:*]
endif else begin
    ret={rot:lonarr(3,3),trn:fltarr(3)}
    for i=0l,ngeneral-1 do $
        ret=[ret,Symop(allops[i])]
    ret=ret[1:*]
endelse
return,ret
end;function SymopAll
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SymopAll64,allops,seitz=seitz
ngeneral=n_elements(allops)
if keyword_set(seitz) then begin
    ret=fltarr(4,4)
    for i=0l,ngeneral-1 do $
        ret=[[[ret]],[[Symop64(allops[i],/seitz)]]]
    ret=ret[*,*,1:*]
endif else begin
    ret=Symop64(allops[0])
    for i=1,ngeneral-1 do $
        ret=[ret,Symop64(allops[i])]
endelse
return,ret
end;function SymopAll64
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Symop_code,Roti,Trni,seitz=seitz

if keyword_set(seitz) then begin
    Rot=Roti[0:2,0:2]
    Trn=Roti[3,0:2]
endif else begin
    Rot=Roti
    Trn=Trni
endelse

; Integerize rotation
rot=round(Rot)

; Integerize translation: rounded to a basis of 48, and put on the range 0..24
; => the interval [0,1[ is mapped to [0,24] and so is [1,2[ ; [-1,0[ ; ...

trn= round(Remain(round(48.*(Trn-1/96.)) ,48) / 2.)

; Example1 of 32bit encoding: ascii
; 32bit Encoding of a vector [u,v,w,...] with base 256 (i.e. ascii sequence):
;    [a,1,*,0] => 808071521 (=97*256^3 + 49*256^2 + 42*256^1 + 48*256^0)
;
; Example2 of 32bit encoding: hexadecimal
;    [a,f,f,a] => 45050 (=10*16^3 + 15*16^2 + 15*16^1 + 10*16^0)
;
; Example3 of 32bit encoding: octal
;    [6,7,7] => 447 (=6*8^2 + 7*8^1 + 7*8^0)
;
; Remark: n-element vector with base b needs a x-bit encoding (x = n*alog(b)/alog(2))


; Encoding for translation vector [u,v,w] with base 24 (minimum 14bit required so 16bit is fine):
; code_t = u*24^0 + v*24^1 + w*24^2
;
; E.g.: trn=[23,5,0] => 143 <= baseconvert('05N',24,10)

fac = [1L,24,576] ;
code_t = Remain(trn ,24) * fac ; Remain is used here because for base 24, elements must be in [0,23] range
code_t = total(code_t,/PRESERVE_TYPE)

; Encoding for rotation with base 3 (minimum 15bit required so 16bit is fine):
; code_r = rot[0,0]*3^0 + rot[1,0]*3^1 + rot[2,0]*3^2 + ...

fac = [[1L,3,9],[27,81,243],[729,2187,6561]]
code_r = Remain( (rot + 1) , 3) * fac ; elements must be 0, 1 or 2
code_r = total(code_r,/PRESERVE_TYPE)

; xor to make identity zero
code = ishft( ( code_r xor '4064'XL ) , 16 ) + code_t
; HO word (i.e. 16bit): rotation code
; LO word: translation code

return,code
end;function Symop_code
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Symop_code64,Roti,Trni,seitz=seitz

if keyword_set(seitz) then begin
    Rot=Roti[0:2,0:2]
    Trn=Roti[3,0:2]
endif else begin
    Rot=Roti
    Trn=Trni
endelse

keep=round(48.*(Rot-1/96.))+96
rot= round(Remain(keep,192)/2.) ; [-2,2[ => [0,96]
                                ; keep(-2)=0, keep(2)=192
                                ; remain192=0  , remain192=0
trn= round(Remain(round(48.*(Trn-1/96.)) ,48) / 2.)  ; [0,1[  => [0,24]

fac = [1LL,24,576]
code_t = Remain(trn ,24) * fac
code_t = total(code_t,/PRESERVE_TYPE)

; In the HO DWORD of code_t, we will keep track of
; which of the 9 elements were 2 (because 2 will be mapped to 0, as is -2)
fac=[1LL,2,4,8,16,32,64,128,256]
code_t2=fac*reform(keep eq 192,9)
code_t2= total(code_t2,/PRESERVE_TYPE)
code_t+=ishft(code_t2,32)

fac=[[1LL,96,9216],[884736,84934656,8153726976],[782757789696,75144747810816,7213895789838336]]
code_r = Remain(rot ,96) * fac
code_r = total(code_r,/PRESERVE_TYPE)

code = {t:code_t, r:code_r xor '7423B088F0ED248'XLL}

return,code
end;function Symop_code64
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StringRTopcode64,codes
n=n_elements(codes)
str=strarr(n)
for i=0l,n-1 do $
    str[i]=StringRTop(Symop64(codes[i]))
return,str
end;function StringRTopcode64
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RTopcodeString64,str
n=n_elements(str)
codes=replicate(CodeType64(0,0),n)
for i=0l,n-1 do begin
    tmp=RTopString(str[i])
    codes[i]=Symop_code64(tmp.rot,tmp.trn)
endfor
return,codes
end;function RTopcodeString64
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function StringRTopcode,codes
n=n_elements(codes)
str=strarr(n)
for i=0l,n-1 do $
    str[i]=StringRTop(Symop(codes[i]))
return,str
end;function StringRTopcode
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function RTopcodeString,str
n=n_elements(str)
codes=replicate(-1L,n)
for i=0l,n-1 do begin
    tmp=RTopString(str[i])
    codes[i]=Symop_code(tmp.rot,tmp.trn)
endfor
return,codes
end;function RTopcodeString
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SymopAllRat,allops,error=error,_extra=ex
return,RTopStringRatAll(strjoin('['+StringRTopcode(allops)+']',';'),error=error,_extra=ex)
end;function SymopAllRat
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isGroupModT,ops

H=SymopAll(ops,/seitz)
nH=n_elements(ops)

; Check for identity
ind=where(ops eq 0,ct)
if ct ne 1 then return,0b

; Check for inverse element
for i=0,nH-1 do begin
    tmp=invert(H[*,*,i])
    code=Symop_code(tmp[0:2,0:2],tmp[3,0:2])
    ind=where(ops eq code,ct)
    if ct ne 1 then return,0b
endfor

; Check for closure
for i=0,nH-1 do begin
    for j=0,nH-1 do begin
        tmp=H[*,*,j]##H[*,*,i]
        code=Symop_code(tmp[0:2,0:2],tmp[3,0:2])
        ind=where(ops eq code,ct)
        if ct ne 1 then return,0b
    endfor
endfor

return,1b
end;function isGroupModT
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function isGroup,ops

H=SymopAll(ops,/seitz)
Float2Rational,H,24,/L64
nH=n_elements(ops)

; Check for identity
ind=where(ops eq 0,ct)
if ct ne 1 then return,0b

; Check for inverse element
for i=0,nH-1 do begin
    tmp=ratmatinvert(H[*,*,i])
    b=0
    for k=0,nH-1 do b+=ratarray_equal(tmp,H[*,*,k])
    if b ne 1 then return,0b
endfor

; Check for closure
for i=0,nH-1 do begin
    for j=0,nH-1 do begin
        tmp=ratmatmult(H[*,*,j],H[*,*,i])
        b=0
        for k=0,nH-1 do b+=ratarray_equal(tmp,H[*,*,k])
        if b ne 1 then return,0b
    endfor
endfor

return,1b

end;function isGroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SysAbs,h,k,l,allops,mult=mult

; Bh=[h,k,l]
;
; Special reflections:
; --------------------
; if Bh##Rot == Bh
;    and (Bh##Trn) mod 1 != 0: systematically absent reflection
;    and (Bh##Trn) mod 1 == 0: systematically enhancement reflection (epsilon is enhancement factor)
; if Bh##Rot == -Bh: phase restriction phi(Bh) = pi*Bh##Trn + n*pi
;
; Remark: Bh is a row vector
;
; Terminology:
; ------------
; acentric reflection: reflection without phase restriction
; centric reflection: reflection with phase restriction
;         (centrosymmetric spacegroup => all reflections are centric)
; epsilon: number of times a reflection is mapped into itself
; multiplicity for acentric reflection: twice the number of symmetrically equivalent Miller indices
; multiplicity for centric reflection: the number of symmetrically equivalent Miller indices
;
; Relations:
; ----------
; (# symmetrically equivalent Miller indices) = (order of spacegroup)/(# times reflection is mapped into itself)

hkl=[h,k,l]
epsilon_ = 1
allowed_ = 255
for i = 1, n_elements(allops)-1 do begin
    op=Symop(allops[i])
    equiv = hkl##op.rot; transpose
    shift = -2*!dpi*hkl##reform(op.trn,1,3) ; sym_phase_shift
    shift = shift[0]
    if array_equal(equiv , hkl) then begin
      ; Bh*Rot=Bh: reflection is related to itself, so it is special or sysabs
      if ( cos(shift) gt 0.999 ) then epsilon_++ $ ;  systematic enhancement
      else begin
          ; sysabs
        epsilon_ = 0
        allowed_ = 0
        break
      endelse
    endif else begin
        if array_equal(equiv , -hkl) then begin
              ; Bh*Rot=-Bh: reflection is related to opposite, so it is centric
              allowed_ = round( Remain(-0.5*shift, !dpi) / (!dpi/12.0) )
          endif
    endelse
endfor
if ( h eq 0 and k eq 0 and l eq 0 ) then allowed_ = 0

bsysabs=epsilon_ eq 0
if bsysabs then mult=0 else mult=n_elements(allops)/epsilon_
acentric = allowed_ eq 255
if acentric then mult*=2
allowed_*=(!dpi/12.0)

;print,hkl
;checkSysAbs,allops,hkl,bsysabs,centric

return,bsysabs
end;function SysAbs
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function GenEquivalentHKL,hkl,RT,count=ct

equiv=hkl
ct=1
for i=0l,n_elements(RT)-1 do begin
    add=hkl##RT[i].rot
    tmp=where((equiv[0,*] eq add[0]) and $
        (equiv[1,*] eq add[1]) and $
        (equiv[2,*] eq add[2]),count)
    if count eq 0 then begin
        equiv=[[equiv],[add]]
        ct++
    endif
endfor

return,equiv
end;function GenEquivalentHKL
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function HKLequivalent,a,b,allops,count=ct
; Returns indices of hkl's in a with an equivalent in b
; return: inda = indices of these hkl's in a
;          indb = indices of corresponding hkl's in b

na=n_elements(a)/3
nb=n_elements(b)/3

; First find all equivalent hkls of b
RT=SymopAll(allops)
bequiv=bytarr(3)
indequiv=[-1]
for i=0l,nb-1 do begin
    equiv=GenEquivalentHKL(b[*,i],RT,count=ct)
    bequiv=[[bequiv],[equiv]]
    indequiv=[indequiv,replicate(i,ct)]
endfor

; Find intersection between a and equivalents of b
intersec=SetIntersection(a,bequiv,/nd,/secondind,count=ct)
if ct ne 0 then begin
    inda=intersec.a
    indb=indequiv[intersec.b]
endif else begin
    inda=-1
    indb=-1
    ct=0
endelse

return,{a:inda,b:indb}
end;function HKLequivalent
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SGFromHall,symb,error=error

error=1b
CATCH, Error_status
IF Error_status NE 0 THEN return,[0L]

; ---- Some defenitions ----

mat_i   =[[ 1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0] ]
mat_inv =[[-1.0, 0.0, 0.0], [0.0,-1.0, 0.0], [0.0, 0.0,-1.0] ]
mat_2z  =[[-1.0, 0.0, 0.0], [0.0,-1.0, 0.0], [0.0, 0.0, 1.0] ]
mat_3z  =[[ 0.0,-1.0, 0.0], [1.0,-1.0, 0.0], [0.0, 0.0, 1.0] ]
mat_4z  =[[ 0.0,-1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0] ]
mat_6z  =[[ 1.0,-1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 1.0] ]
mat_2q  =[[ 0.0,-1.0, 0.0],[-1.0, 0.0, 0.0], [0.0, 0.0,-1.0] ]
mat_2qq =[[ 0.0, 1.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0,-1.0] ]
mat_3abc=[[ 0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0] ]

vec_0=[ 0.0, 0.0, 0.0 ]
vec_a=[ 0.5, 0.0, 0.0 ]
vec_b=[ 0.0, 0.5, 0.0 ]
vec_c=[ 0.0, 0.0, 0.5 ]
vec_n=[ 0.5, 0.5, 0.5 ]
vec_u=[ 0.25, 0.0, 0.0 ]
vec_v=[ 0.0, 0.25, 0.0 ]
vec_w=[ 0.0, 0.0, 0.25 ]
vec_d=[ 0.25, 0.25, 0.25 ]

vec_AA=[ 0.0, 0.5, 0.5 ]
vec_BB=[ 0.5, 0.0, 0.5 ]
vec_CC=[ 0.5, 0.5, 0.0 ]
vec_I=[ 0.5, 0.5, 0.5 ]
vec_R=[ 0.666667, 0.333333, 0.333333 ]
vec_S=[ 0.333333, 0.666667, 0.333333 ]
vec_T=[ 0.333333, 0.333333, 0.666667 ]
vec_H=[ 0.666667, 0.333333, 0.0 ]
vec_F1=[ 0.0, 0.5, 0.5 ]
vec_F2=[ 0.5, 0.0, 0.5 ]

bspace = (byte(' '))[0] & b0=(byte('0'))[0] & bstar=(byte('*'))[0]
bquote = (byte("'"))[0] & bdblquote = (byte('"'))[0] & bmin=(byte('-'))[0]
b1=(byte('1'))[0] & b2=(byte('2'))[0] & b3=(byte('3'))[0]
b4=(byte('4'))[0] & b5=(byte('5'))[0] & b6=(byte('6'))[0]
ba=(byte('a'))[0] & bb=(byte('b'))[0] & bc=(byte('c'))[0]
bn=(byte('n'))[0] & bu=(byte('u'))[0] & bv=(byte('v'))[0]
bw=(byte('w'))[0] & bd=(byte('d'))[0]
bx=(byte('x'))[0] & bz=(byte('z'))[0] & by=(byte('y'))[0]

; ---- Interpret Hall symbol ----
; separate the Hall symbol from the change of basis
tokens = str_sep(symb,"(",/trim)
n=n_elements(tokens)
if ( n gt 0 ) then sym = tokens[0]
if ( n gt 1 ) then begin
    chb = strmid(tokens[1],0,strlen(tokens[1])-1);remove ")"
    i=strpos(chb,',')
    if ( i ne -1 ) then begin
      tmp = RTopString( chb )
      cbop=MakeSeitz(tmp.rot,tmp.trn)
    endif else begin
      cbop=MakeSeitz(mat_i,vec_0)
      t = str_sep(chb," ")
      for i=0l,(n_elements(t)<3)-1 do cbop[3,i] = long(t[i]) / 12.0;
    endelse
    cbopi=invert(cbop)
endif

; now separate the parts of the Hall symbol
tokens = str_sep(sym," ",/trim)

; first part: lattice and inversion
inv = bspace
lat = bspace
token = byte(strupcase(tokens[0]))
for j = 0, n_elements(token)-1  do begin
    c=token[j]
    if ( c eq bmin ) then inv = c else lat = c
endfor

; next 1-3 parts: generating matrices
nmat = n_elements(tokens)-1

invop = replicate(bspace,nmat)
order = invop
screw = invop
axis1 = invop
axis2 = invop
trans = strarr(nmat)
for i = 0, nmat-1 do begin
    token = byte(strlowcase(tokens[i+1]))
    for j = 0, n_elements(token)-1 do begin
      c=token[j]
      if ( c eq bmin ) then invop[i] = c $
      else if ( c ge b1 and c le b6 ) then $
        if ( order[i] eq bspace ) then order[i] = c $
        else                   screw[i] = c $
      else if ( c ge bx and c le bz ) then axis1[i] = c $
      else if ( c ge bdblquote and c le bstar ) then axis2[i] = c $
      else if ( c ge ba and c le bw ) then trans[i] = trans[i]+string(c)
    endfor
endfor

; now interpret all the default conventions
  ; default first axis to z
  if ( nmat ge 1 ) then begin
    if ( axis1[0] eq bspace ) then axis1[0] = bz
  endif
  ; default second axis on basis of first
  if ( nmat ge 2 ) then begin
    if ( axis1[1] eq bspace ) then begin
      if ( order[1] eq b2 ) then begin
        if        ( order[0] eq b2 or order[0] eq b4 ) then begin
          if ( axis2[1] eq bspace ) then axis1[1] = bx
        endif else if ( order[0] eq b3 or order[0] eq b6 ) then begin
          axis1[1] = axis1[0];
          if ( axis2[1] eq bspace ) then axis2[1] = bquote
        endif
      endif else if ( order[1] eq b3 ) then begin
        if ( order[0] eq b2 or order[0] eq b4 ) then begin
          if ( axis2[1] eq bspace ) then axis2[1] = bstar
        endif
      endif
    endif
  endif
  ; default third axis (not much choice here)
  if ( nmat ge 3 ) then begin
    if ( axis1[2] eq bspace ) then begin
      if ( order[2] eq b3 ) then begin
        if ( axis2[2] eq bspace ) then axis2[2] = bstar
      endif
    endif
  endif

; now check we have everything
for i = 0, nmat-1 do begin
    ; add fake z axes for non-axis ops
    if ( axis1[i] eq bspace ) then begin
      if ( order[i] eq b1 ) then axis1[i] = bz;  // fake axis
      if ( axis2[i] ne bspace ) then axis1[i] = bz ; // fake axis
    endif
    ; error messages
    if ( axis1[i] eq bspace ) then begin
        message,"SGFromHall: Missing x/y/z in Hall symb:"+tokens[i+1]
    endif
    if ( order[i] eq bspace ) then begin
        message,"SGFromHall: Missing order in Hall symb:"+tokens[i+1]
    endif
endfor

; ---- Make list of R/T-operations ----
; add identity and inversion
ops=0L ; Symop_code(mat_i,vec_0) eq 0
if ( inv eq bmin ) then ops=[ops,Symop_code(mat_inv,vec_0)]

; now make the generator matrices
matperm=fltarr(3,3)
for i = 0, nmat-1 do begin
    ; make matrix part
    mat = mat_i;
    if ( order[i] eq b2 and axis2[i] eq bspace ) then mat = mat_2z;
    if ( order[i] eq b2 and axis2[i] eq bquote) then mat = mat_2q;
    if ( order[i] eq b2 and axis2[i] eq bdblquote ) then mat = mat_2qq;
    if ( order[i] eq b3 and axis2[i] eq bspace ) then mat = mat_3z;
    if ( order[i] eq b3 and axis2[i] eq bstar ) then mat = mat_3abc;
    if ( order[i] eq b4 and axis2[i] eq bspace ) then mat = mat_4z;
    if ( order[i] eq b6 and axis2[i] eq bspace ) then mat = mat_6z;
    ; inverse (improper)
    if ( invop[i] eq bmin ) then mat = mat_inv ## mat;
    ; axis permutation
    case axis1[i] of
    bx:j=0
    by:j=1
    else:j=2
    endcase

    for k=0l,2 do $
      for l=0,2 do $
        matperm( l,k ) = mat( (l+2-j) mod 3 ,(k+2-j) mod 3 )
    ; make translation part
    vec = vec_0;
    for k=0l,strlen(trans[i])-1 do begin
        tmp=strmid(trans[i],k,1)
        case tmp of
        'a':vec += vec_a;
        'b':vec += vec_b;
        'c':vec += vec_c;
        'n':vec += vec_n;
        'u':vec += vec_u;
        'v':vec += vec_v;
        'w':vec += vec_w;
        'd':vec += vec_d;
        else:
        endcase
    endfor
    ; screw translation
    if ( screw[i] ne bspace ) then vec[j] += float( screw[i] - b0 ) / ( order[i] - b0 )
    ; store the matrix
    ops=[ops,Symop_code(matperm, vec)]
endfor

; add lattice centering
lat=string(lat)
if (lat eq 'A') then ops=[ops,Symop_code(mat_i,vec_AA)]
if (lat eq 'B') then ops=[ops,Symop_code(mat_i,vec_BB)]
if (lat eq 'C') then ops=[ops,Symop_code(mat_i,vec_CC)]
if (lat eq 'I') then ops=[ops,Symop_code(mat_i,vec_I)]
if (lat eq 'R') then ops=[ops,Symop_code(mat_i,vec_R)]
if (lat eq 'S') then ops=[ops,Symop_code(mat_i,vec_S)]
if (lat eq 'T') then ops=[ops,Symop_code(mat_i,vec_T)]
if (lat eq 'H') then ops=[ops,Symop_code(mat_i,vec_H)]
if (lat eq 'Q') then ops=[ops,Symop_code(mat_i,vec_S)]
if (lat eq 'F') then begin
    ops=[ops,Symop_code(mat_i,vec_F1)]
    ops=[ops,Symop_code(mat_i,vec_F2)]
endif

; apply the change of basis operator
if n gt 1 then begin
    Seitz=fltarr(4,4)
    Seitz[3,3]=1
    for i=0l,n_elements(ops)-1 do begin
        tmp=Symop(ops[i])
        Seitz[0:2,0:2]=tmp.rot
        Seitz[3,0:2]=tmp.trn
        tmp=cbop##Seitz##cbopi
        ops[i] = Symop_code( tmp[0:2,0:2],tmp[3,0:2] );
    endfor
endif

error=0b
return,ops
end;function SGFromHall
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function NTHallFromOps,ops
; reverse SGFromHall

return,'Not tabulated.'
end;function NTHallFromOps
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SGFromNr,symb,error=error

symbx=sgdata(fix(symb),2)
ops=SGFromHall(symbx.hall,error=error)

return,ops
end;function SGFromNr
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SGFromHM,symb,error=error

hm=strsplit(symb,':',/extract,count=ct)
if ct ne 1 then begin
    hm[1]=strupcase(hm[1])
    case hm[1] of
    '1': pref_12=hm[1]
    '2': pref_12=hm[1]
    'H': pref_hr=hm[1]
    'R': pref_hr=hm[1]
    endcase
    hm=hm[0]
endif else hm=hm[0]

symbx=sgdata(hm,1,pref_12=pref_12,pref_hr=pref_hr)
if symbx.sghash eq 0 then ops=SGFromNr(shmsymbol(hm),error=error) $
else ops=SGFromHall(symbx.Hall,error=error)

return,ops
end;function SGFromHM
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SGFromOp,symb,error=error

symops = str_sep(symb,';')
n=n_elements(symops)

error=n eq 0
if error then return,[0L]

ops=0L
for i=0l,n-1 do begin
    tmp=RTopString(symops[i],error=error)
    if error then return,[0L]
    ops=[ops,Symop_code(tmp.rot,tmp.trn)]
endfor

return,ops[1:*]
end;function SGFromOp
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SGFromTab,symb,error=error

symbx=sgdata(fix(symb),3)
ops=SGFromHall(symbx.hall,error=error)

return,ops
end;function SGFromTab
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SGFromHash,symb,error=error

symbx=sgdata(symb,0)
ops=SGFromHall(symbx.hall,error=error)

return,ops
end;function SGFromHash
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ProductOps,ops1,ops2

n1=n_elements(ops1)
n2=n_elements(ops2)

for j=0l,n2-1 do begin
    if ops2[j] ne 0 then $
    for i=0l,n1-1 do begin
        r1=Symop(ops1[i])
        r2=Symop(ops2[j])
        k = Symop_code( r1.rot##r2.rot, r1.rot##r2.trn+r1.trn);
        ops1=[ops1,k]
    endfor
endfor

end;pro ProductOps
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ExpandOpsList,generators
ops=0L;  // identity is compulsory

; generate all the symops
repeat begin
    size = n_elements(ops)
    if size gt 1000 then return,generators ; probably something wrong
    for i=0l,n_elements(generators)-1 do begin
      if ( generators[i] ne 0 ) then begin
        for j=0l,size-1 do begin
          r1=Symop(generators[i])
          r2=Symop(ops[j])
          k = Symop_code(r1.rot##r2.rot , r1.rot##r2.trn+r1.trn)
          l=(where(ops eq k))[0]
          if l eq -1 then ops=[ops,k]
        endfor
      endif
      endfor
endrep until (n_elements(ops) le size)
return, ops;
end;function ExpandOpsList
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function primitive_noninv_ops,ops

; if space group isn't centrosymmetric: use primitive_ops instead (faster and same result)
pops=0L
for i=0l,n_elements(ops)-1 do begin
    ind=where(RotCode(pops) eq RotCode(ops[i]),ct)
    if ( ct eq 0 ) then begin
        det=determ(float((Symop(ops[i])).rot))
        if det gt 0 then pops=[pops, ops[i]]
    endif
endfor
return, pops;
end;function primitive_noninv_ops
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function inversion_ops,ops

mat_inv =[[-1.0, 0.0, 0.0], [0.0,-1.0, 0.0], [0.0, 0.0,-1.0] ]
vec_0=[ 0.0, 0.0, 0.0 ]

ind=where(RotCode(ops) eq RotCode(Symop_code(mat_inv,vec_0)),ct)

if ct eq 0 then return, 0L else return,[0L,ops[ind[0]]]

end;function inversion_ops
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function primitive_ops,ops
pops=0L
for i=0l,n_elements(ops)-1 do begin
    ind=where(RotCode(pops) eq RotCode(ops[i]),ct)
    if ( ct eq 0 ) then pops=[pops, ops[i]]
endfor
return, pops;
end;function primitive_ops
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function centering_ops,ops

ind=where(Rotcode(ops) eq 0,ct)
if ct eq 0 then return,dummy else return,ops[ind]

end;function centering_ops
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function expanded_centering_vecs,allops,nT,ncen=ncen

tops=centering_ops(allops)
ncen=n_elements(tops)
nT=0
if ncen eq 1 then return,[0,0,0]

T=(SymopAll(tops)).trn
; Substract 1 from all non-zero components
for i=1,ncen-1 do begin
    ind=where(T[*,i] ne 0,ct)
    n=2^ct
    sub=bytarr(3,n)
    tmp=indgen(n)
    for j=0l,ct-1 do sub[ind[j],*]=(tmp and 2^j) ne 0
    sub=sub[*,1:*]
    T=[[T],[rebin(T[*,i],3,n-1,/sample)-sub]]
endfor
T=T[*,1:*]
; Eliminate multiples
T=UniqueRows(T)
; Sort by norm (lowest first)
T=T[*,sort(total(T*T,1,/pres))]
; Add unit translations
T=[[T],[1,0,0],[0,1,0],[0,0,1]]

nT=(size(T,/dim))[1]
return,T
end;function expanded_centering_vecs
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sgchangebasis,allopsi,M,Mi
allops=allopsi

; XX = A.X+B     ; change of frame (old frame as a function of the new)
; Y  = R.X+T     ; symmetry operation
; YY = RR.XX+TT  ; symmetry operation in new frame
;    = A.Y + B
;    <=>
;      RR=A.R.A^(-1)
;      TT=B-RR.B+A.T
;
; Only a change of basis (B=0)
; A=Mi (old basis vectors in function of the new basis vectors)

for i=0l,n_elements(allops)-1 do begin
    RT=Symop(allops[i])
    
    RR=RT.rot
    Float2Rational,RR,1,/L64
    RR=ratmatmult(Mi,ratmatmult(RR,M))
    
    TT=transpose(RT.trn)
    Float2Rational,TT,24,/L64
    TT=ratmatmult(Mi,TT)
    
    if total((24 mod TT.denom) ne 0,/int) ne 0 then stop
    
    allops[i]=Symop_Code(RR.num/float(RR.denom),TT.num/float(TT.denom))
endfor

return,ExpandOpsList(allops)

end;function sgchangebasis
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function sgmakeprimitive,allops,afftrans=info
; trivial: 0(AB not trivial)
;          1(AB trivial)
;          2(AB trivial but no primitive sg found)

; Primitive SG: i.e. no pure translations except for [0,0,0]
; Expanded list of centering operations
T=expanded_centering_vecs(allops,nT,ncen=ncen)
info={trivial:1b,M:identityrat(3,/L64),Mi:identityrat(3,/L64),ncen:ncen}
if nT eq 0 then return,allops
nallops=n_elements(allops)
R=(SymopAll(allops)).rot

; Non-primitive basis
; Basis of vector space V: B={b1,...,bn}  (bi: column vector)
; Vector v: v=x.b1+y.b2+...=[b1,b2,...].[[x],[y],..]
; Coordinates of v with respect to B: [v]B=[[x],[y],...]
; Unit Volume: Vb=det([b1,b2,...])

; Primitive basis
; Basis of vector space V: C={c1,...,cn}
; Vector v: v=x.c1+y.c2+...=[c1,c2,...].[[x'],[y'],..]
; Coordinates of v with respect to C: [v]C=[[x'],[y'],...]
; Unit Volume: Vc=det([c1,c2,...])

; [v]C = [[x'],[y'],...] = [[b1]C,[b2]C,...].[[x],[y],..]
;                        = ([[c1]B,[c2]B,...])^(-1).[[x],[y],..]
;                         = M^(-1).[[x],[y],..]
; Vb=Vc.det(M^(-1))=Vc.ncen

; We will take all possible combinations of T's as basis vectors. Concatenate
; these column vectors and invert this matrix and we have the transformation
; matrix from old (non-primitive) to new (primitive) basis.
; Accept only when det(M^(-1))==ncen

for i=0l,nT-1 do $
    for j=i+1,nT-1 do $
        for k=j+1,nT-1 do begin
            M=transpose(T[*,[i,j,k]]) ; New basis vectors as a function of the old ones
            Float2Rational,M,24,/L64
            det=ratmatdet(M)
            
            if det.denom eq ncen and det.num eq 1 then begin
;            if det.denom eq ncen and abs(det.num eq 1) then begin
;                if det.num lt 0 then begin
;                    ; Make determinant positive
;                    tmp=M[0,*]
;                    ratneg,tmp
;                    M[0,*]=temporary(tmp)
;                endif
                Mi=ratmatinvert(M) ; Old basis vectors as a function of the new ones

                ; Check whether the new symmetry operations have
                ; elements 0, 1 or -1 but no other.
                bool=1b
                for l=0,nallops-1 do begin
                    tmp=reform(R[*,*,l])
                    rational,tmp,/L64
                    Rnew=ratmatmult(Mi,ratmatmult(tmp,M))
                    bool and= total(abs(Rnew.num) gt 1,/int) eq 0
                    bool and= total(Rnew.denom ne 1,/int) eq 0
                endfor

                ; Check whether centering ops have vanished
                if bool then begin
                    info.M=M
                    info.Mi=Mi
                    allopsnew=sgchangebasis(allops,M,Mi)
                    tmp=centering_ops(allopsnew)
                    if n_elements(tmp) eq 1 then begin
                        info.trivial=0b
                        return,allopsnew
                    endif
                endif
            endif
        endfor

info.trivial=2b
return,allops
end;function sgmakeprimitive
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function generator_ops,ops
cens = centering_ops(ops);
prms = primitive_ops(ops);
gens = 0L
gend = ExpandOpsList(gens)

; first make the centering generators
for i=1,n_elements(cens)-1 do begin
    for j=0l,n_elements(gend)-1 do $
      if ( cens[i] eq gend[j] ) then goto, generator_opsskip1;
    gens=[gens,cens[i]]
    gend = ExpandOpsList(gens)
    if ( n_elements(gend) eq n_elements(cens) ) then goto,generator_opsskipfor1;  // optional optimisation
  generator_opsskip1:
endfor
generator_opsskipfor1:
ncen = n_elements(gens)

; now add the rest of the generators
for i=1,n_elements(prms)-1 do begin
    for j=0l,n_elements(gend)-1 do $
      if ( prms[i] eq gend[j] ) then goto, generator_opsskip2;
    gens=[gens,prms[i]]
    gend = ExpandOpsList(gens);
    if ( n_elements(gend) eq n_elements(ops) ) then goto,generator_opsskipfor2;  // optional optimisation
  generator_opsskip2:;
endfor
generator_opsskipfor2:

generator_opsback:  ;// finally remove any redundent ops
for i=ncen,n_elements(gens)-1 do begin
    genp = gens;
    tmp=lindgen(n_elements(genp))
    ind=where(tmp ne i)
    genp=genp[ind]
    if ( n_elements(ExpandOpsList(genp)) eq n_elements(ops) ) then begin
      gens = genp;
      goto, generator_opsback;
    endif
endfor

return, gens;  // return result
end;function generator_ops
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function laue_ops,ops

mat_inv =[[-1.0, 0.0, 0.0], [0.0,-1.0, 0.0], [0.0, 0.0,-1.0] ]
vec_0=[ 0.0, 0.0, 0.0 ]
lops=Symop_code(mat_inv,vec_0)
lops=[lops, RotCode(ops)]

elops=ExpandOpsList(lops)

return,elops[sort(elops)]
end;function laue_ops
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddSubgroup,p_subgroups,ind
bool=0b
for j=0l,n_elements(p_subgroups)-1 do bool or= array_equal(ind,(*(p_subgroups[j])))
if ~bool then p_subgroups=[p_subgroups,ptr_new(ind)]
end;pro AddSubgroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro MakeClosure,mult_table,ind

ind=ind[uniq(ind,sort(ind))]

nh=n_elements(ind)
repeat begin
    subbin=(mult_table[ind,*])[*,ind]
    h=histogram(subbin,rev=rev)
    ind=subbin[rev[rev[where(h ne 0,nh2)]]]
    bool=nh eq nh2
    nh=nh2
endrep until bool
ind=ind[sort(ind)]

end;pro MakeClosure
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PGConjugacyClasses,ops

; The ops must be sorted, so that the first is E
n=n_elements(ops)

; Op <-> Op^(-1)
; Table: Opi *  Opj
; Table: Opi *  Opj * Opi^(-1)
binverse=intarr(n)
conjug_table=intarr(n,n)
mult_table=intarr(n,n)
for i=0l,n-1 do begin
    ; Row i
    oprot=(Symop(ops[i])).rot
    ioprot=invert(oprot)
    ind=where(ops eq Symop_code(ioprot,[0,0,0]))
    binverse[i]=ind

    for j=0l,n-1 do begin
        ; Column j
        tmp=oprot##(Symop(ops[j])).rot
        ind=where(ops eq Symop_code(tmp,[0,0,0]))
        mult_table[j,i]=ind ; gi*gj

        ind=where(ops eq Symop_code(tmp##ioprot,[0,0,0]))
        conjug_table[j,i]=ind ; gi*gj*gi^(-1)
    endfor
endfor

; Subgroups: Identity and inverse element -> {E,Opi,Opi^(-1)}
subgroups=[[intarr(n)],[(indgen(n))],[binverse]]
h=histogram(indgen(n)<binverse,rev=rev)
subgroups=subgroups[rev[rev[where(h ne 0,nsub)]],*]

; Subgroups: Closure
p_subgroups=ptr_new([0]) ; subgroups[0,*]=[0,0,0]={E,E,E} => always a subgroup
for i=1,nsub-1 do begin
    ind=subgroups[i,*]

    ; Check closure
    MakeClosure,mult_table,ind

    ; Add when not already there
    AddSubgroup,p_subgroups,ind
endfor
nsub=n_elements(p_subgroups)

; Subgroups: All possible unions
repeat begin
    ; skip unions with {E} because E is in all subgroups
    for i=1,nsub-2 do begin
        indi=*(p_subgroups[i])

        for j=i+1,nsub-1 do begin
            ind=[indi,*(p_subgroups[j])]

            ; Check closure
            MakeClosure,mult_table,ind

            ; Add when not already there
            AddSubgroup,p_subgroups,ind
        endfor
    endfor
    nsub2=n_elements(p_subgroups)
    bool=nsub eq nsub2
    nsub=nsub2
endrep until bool

; Conjugacy classes of subgroups
p_subgroups2=ptr_new(p_subgroups[0]) ; {E} forms its own conjugacy class
nconj=1
p_subgroups[0]=ptr_new()
nclasses=1
for i=1,nsub-1 do begin
    if ptr_valid(p_subgroups[i]) then begin
        group=*(p_subgroups[i])
        
        p_subgroups2=[p_subgroups2,ptr_new(p_subgroups[i])]
        nconj=[nconj,1]
        p_subgroups[i]=ptr_new()
        nclasses++

        for k=1l,n-1 do begin ; conjugacy with identity already considered
            ; conjug_group = gk*group*gk^(-1)
            tmp=conjug_table[group,k]
            conjug_group=tmp[uniq(tmp,sort(tmp))]

            ; remove conjug_group from list and add to conjugacy class of group
            ind2=where(ptr_valid(p_subgroups),ct)
            if ct ne 0 then begin
                bool=bytarr(ct)
                for j=0l,ct-1 do bool[j] or= array_equal(conjug_group,(*(p_subgroups[ind2[j]])))
                ind3=where(bool,ct)

                if ct ne 0 then begin
                    ind2=ind2[ind3]
                    *p_subgroups2[nclasses-1]=[*p_subgroups2[nclasses-1],p_subgroups[ind2]]
                    nconj[nclasses-1]+=ct
                    p_subgroups[ind2]=ptr_new()
                endif
            endif
        endfor
    endif
endfor

return,p_subgroups2

end;function PGConjugacyClasses
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function InequalModZ,b,d,Q,nsol=nsol

; x  = Q.x'
; x' = (b-z)/d mod Z
;
; d==0 => x in Q
; d>0  => b-d<  z <=b
; d<0  => b  <= z < b-d

bmind=b-d

first=(bmind<b)
first=ceil(first)+(((first mod 1) eq 0) and (d gt 0))

last=(bmind>b)
last=floor(last)-(((last mod 1) eq 0) and (d lt 0))

n=(last-first+1)>1

; all combinations
nsol=n[0]*n[1]*n[2]
if nsol eq 0 then return,CodeType64(0,0)
rot=fltarr(3,3)

if d[0] eq 0 then begin
    x=0
    n[0]=1
    rot[0,*]=Q[0,*]
endif else begin
    x=(b[0]-first[0]-indgen(n[0]))/float(d[0])
    x=x[sort(x)]
endelse

if d[1] eq 0 then begin
    y=0
    n[1]=1
    rot[1,*]=Q[1,*]
endif else begin
    y=(b[1]-first[1]-indgen(n[1]))/float(d[1])
    y=y[sort(y)]
endelse

if d[2] eq 0 then begin
    z=0
    n[2]=1
    rot[2,*]=Q[2,*]
endif else begin
    z=(b[2]-first[2]-indgen(n[2]))/float(d[2])
    z=z[sort(z)]
endelse

rot=BlockDiag(rot,/onlycolumns) ; so that e.g. [0,z,1/2] becomes [0,x,1/2]

off=0
solutions=replicate(CodeType64(0,0),nsol)
for i=0l,n[0]-1 do $
    for j=0l,n[1]-1 do $
        for k=0l,n[2]-1 do begin
            solutions[off]=Symop_code64(rot,Q##[[x[i]],[y[j]],[z[k]]])
            off++
        endfor

return,solutions[SortCodeType64(solutions)]
end;function InequalModZ
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddOrbit,p_subgroups,neworbit,order,sortmeasure,badd=badd
bool=0b

new=neworbit[SortCodeType64(neworbit)]
; Skip the first, this is just a dummy
for j=1,n_elements(p_subgroups)-1 do begin
    oldj=*(p_subgroups[j])
    oldj=oldj[SortCodeType64(oldj)]

    ; new equal to oldj?
    bool or= array_equal(new.r,oldj.r) and array_equal(new.t,oldj.t)
endfor

badd=~bool
if badd then begin
    p_subgroups=[p_subgroups,ptr_new(neworbit)]
    ordernew=n_elements(neworbit)
    order=[order,ordernew]
    
    ; Sorting of Wyckoff positions:
    ; 1. Highest multiplicity last
    ; 2. Highest degrees of freedom last
    ; 3. Highest amount of non-zero's last
    ; 4. Highest amount of total(translations) last
    
;    print,'['+stringrtopcode64(neworbit)+']'
    trn=(SymopAll64(neworbit)).trn
    rot=(SymopAll64(neworbit[0])).rot
    
    a=long64(ordernew)
    b=total(abs(rot),/int)
    c=total(trn ne 0,/int)
    d=ceil(total(trn))
    measure=(ishft(a,32) and 'FF000000'xl) or (ishft(b,16) and 'FF0000'xl) or (ishft(c,8 ) and 'FF00'xl) or (d and 'FF'xl)

    sortmeasure=[sortmeasure,measure]
endif
end;pro AddOrbit
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function TransformSubspace,AA,BB,CC,DD,code=code
;     A.X  + B = C.X' + D
;     A.X = C.X' + D - B
;     X = E.X' + F

A=AA
E=fltarr(4,3) ; [E|F]
C=[CC,reform(DD-BB,1,3)] ; [C|D-B]

repeat begin
    h=histogram(total(A ne 0,1,/int),min=0,max=3,binsize=1,rev=rev)
    
    ; Single X: e.g.
    ;      a0.x = c0.x' + c1.y' + c2.z' + c3
    ;        x = r0.x' + r1.y' + r2.z' + r3
    for i=0l,h[1]-1 do begin
        row=rev[rev[1]+i]
        col=where(A[*,row] ne 0,ct)
        
        ; Eliminated in previous iteration(s)?
        if ct eq 1 then begin
            col=col[0]
            r=C[*,row]/A[col,row]
            E[*,col]=r
            
            ; Eliminate:
            row=where(A[col,*] ne 0,ct)
            for ii=0,ct-1 do begin
                ; e.g.:      a0.x + a1.y + a2.z = c0.x' + c1.y' + c2.z' + c3
                ;         and  x = r0.x' + r1.y' + r2.z' + r3
                ;              0.x + a1.y + a2.z = (c0-a0.r0).x'  + (c1-a0.r1).y' + (c2-a0.r2).z' + (c3-a0.t)
                C[*,row[ii]]-= r*A[col,row[ii]]
                A[col,row[ii]] = 0
            endfor
        endif
    endfor
    
endrep until h[1] eq 0

C[3,*] mod= 1
if h[0] ne 3 or ~floatequal(C,0.,epsm=10) then return,0b

if total((abs(E[0:2,0:2]) gt 2),/int) ne 0 then stop
code=Symop_code64(E[0:2,0:2],E[3,*])
return,1b
end;function TransformSubspace
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CompareSubspace,RT1,RT2,code=code

; Return value:
; 0: No relation
; 1: Keep the first because the second is a special case or the same
; 2: Keep the second because the first is a special case

; Variable x, y and z component?
var1=total(RT1.Rot ne 0,1,/int) ne 0
var2=total(RT2.Rot ne 0,1,/int) ne 0
ind=where(var1,nvar1)
ind=where(var2,nvar2)

; Check degrees of freedom
ind=where(total(abs(RT1.Rot),2,/pres) ne 0,dof1)
ind=where(total(abs(RT2.Rot),2,/pres) ne 0,dof2)
ret=(dof1 ge dof2)?1:2

; Order: dof2 <= dof1
if nvar1 ge nvar2 then begin
    if ret eq 2 then return,0b
    
    A=RT1.Rot
    B=RT1.Trn
    C=RT2.Rot
    D=RT2.Trn
endif else begin
    if ret eq 1 then return,0b
    
    C=RT1.Rot
    D=RT1.Trn
    A=RT2.Rot
    B=RT2.Trn
    
    tmp=var1
    var1=var2
    var2=temporary(tmp)
endelse

; Mapping of subspaces:
; 
; X(dof=3) --{A|B}--> X1(dof1<=3) --{R|T}--> X2(dof2<=dof1)
; |___________________{C|D}__________________|
;
; subspace X1 = A.X  + B
; subspace X2 = C.X' + D
; 
; If X2 in X1 then we can write:
;     X1 = X2
;     A.X  + B = C.X' + D
;     X = E.X' + F

; Preliminar checks:
h=histogram(var1+var2,min=0,max=2,binsize=1,rev=rev)

; h=0 => both fixed: check whether they are the same
if h[0] ne 0 then begin
    ind=rev[rev[0]:rev[1]-1]
    if ~array_equal(B[ind],D[ind]) then return,0b
endif

; h=1 => one fixed, one variable: check whether the fixed are on the
;                                 lowest dof side and the variable on
;                                  the highest
if h[1] ne 0 then begin
    ind=rev[rev[1]:rev[2]-1]
    t=total(var1[ind],/pres)
    if t ne h[1] then return,0b
endif

if ~TransformSubspace(A,B,C,D,code=code) then return,0b

return,ret
end;function CompareSubspace
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckRTTable,table,n

; Table is an nxn matrix of RT-codes
ind1=where(table[*,0].t ne -1,ct1)
ind2=where(table[*,0].r ne -1,ct2)
if ~array_equal(ind1,ind2) or (ct1 eq 0) then return,0b
ind=ind1
ct=ct1

; Check whether the two dimensions have at least 1 RT in common
for i=0l,ct-1 do begin
    tmp=(table.t eq table[ind[i]].t) and (table.r eq table[ind[i]].r)
    if array_equal(total(tmp,1,/pres),total(tmp,2,/pres)) then return,1b
endfor

return,0b
end;function CheckRTTable
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DeleteSubspace,codes1,codes2
;print,''
;print,'['+stringrtopcode64(codes1)+']'
;print,'['+stringrtopcode64(codes2)+']'

n1=n_elements(codes1)
n2=n_elements(codes2)
if n1 ne n2 then return,0
n=n1

table=replicate(CodeType64(-1,-1),n,n)
keep=bytarr(n,n)

A2=SymopAll64(codes2)
for i=0l,n-1 do begin
    A1i=Symop64(codes1[i])
    for j=0l,n-1 do begin
        ;print,'['+stringrtopcode64(codes1[i])+']','['+stringrtopcode64(codes2[j])+']'
        ret=CompareSubspace(A1i,A2[j],code=code)
        ;print,ret
        if ret ne 0 then table[i,j]= code
        keep[i,j]=ret
    endfor
endfor

; Delete code1 or code2?
h=histogram(keep,binsize=1,min=0,max=2,rev=rev)
if h[1] ne 0 or h[2] ne 0 then $
if h[1] eq 0 or h[2] eq 0 then $
if CheckRTTable(table,n) then return,max(keep)

return,0
end;function DeleteSubspace
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro ReduceGorbit,p_wyckoff,wyckoffmult,sortmeasure
norbit=n_elements(p_wyckoff)
for i=1,norbit-1 do begin ;skip dummy
    for j=i+1,norbit-1 do begin
        if ptr_valid(p_wyckoff[j]) and ptr_valid(p_wyckoff[i]) then begin
            type=DeleteSubspace(*p_wyckoff[i],*p_wyckoff[j])
            case type of
            1:    begin
                ptr_free,p_wyckoff[j]
                p_wyckoff[j]=ptr_new()
                endcase
            2:    begin
                ptr_free,p_wyckoff[i]
                p_wyckoff[i]=ptr_new()
                endcase
            else:
            endcase
        endif
    endfor
endfor

ind=where(ptr_valid(p_wyckoff))
wyckoffmult=wyckoffmult[ind]
sortmeasure=sortmeasure[ind]
p_wyckoff=p_wyckoff[ind]
end;pro ReduceGorbit
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AddGorbit,allops,affsubspace,p_wyckoff,wyckoffmult,sortmeasure,index

G=SymopAll(allops,/seitz)
nG=n_elements(allops)
A=SymopAll64(affsubspace,/seitz)
nA=n_elements(affsubspace)

; G-Orbit of all affine subspaces
badded=0b
for i=0l,nA-1 do begin
    ; E.Aff_i
    OrbAi=affsubspace[i]

    ; Orbit g.Aff_i
    for j=1,nG-1 do begin
        GA=G[*,*,j]##A[*,*,i]
        GAcode=Symop_code64(GA[0:2,0:2],GA[3,0:2])

        ; Add to G-orbit when not the same
        ind=where(OrbAi.r eq GAcode.r and OrbAi.t eq GAcode.t,ct)
        if ct eq 0 then OrbAi=[OrbAi,GAcode]
    endfor

    ; Add orbit g.Aff_i when not already there
    ; and when is has the correct size
    norbit=n_elements(OrbAi)
    if norbit gt index then stop
    
;    print,'G-orbit'+stringr(i),' (',((norbit ne index)?'not added':'added'),'): ','['+stringrtopcode64(OrbAi)+']'
    if (norbit eq index) then begin
        AddOrbit,p_wyckoff,OrbAi,wyckoffmult,sortmeasure,badd=badd
        badded or= badd
    endif
endfor

if ~badded then return

; Check whether G-orbits are not the same or a special case of one another
ReduceGorbit,p_wyckoff,wyckoffmult,sortmeasure

end;pro AddGorbit
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wyckofflabel,j

L=string(byte((j mod 26)+97))
n=j/26
for i=0l,n-1 do L+=L
return,L

end;function wyckofflabel
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function wyckoffcodes,allopsi,wyckoffmult=wyckoffmult
; Used to make the function wyckoff1,...,wyckoff530

wyckoffmult=[0]

; Make primitive
allops=sgmakeprimitive(allopsi,afftrans=afftrans)
if afftrans.trivial eq 2 then begin
    tmp=dialog_message("Can't find a primitive unit cell for this spacegroup.")
    return,0
endif

; Point group Q
lops=RotCode(allops)

elops=ExpandOpsList(lops)
elops=elops[sort(elops)]

nQ=n_elements(elops)
if nQ ne n_elements(allops) then stop

; Conjugacy classes of subgroups of Q
p_conj=PGConjugacyClasses(elops)

; Solve system of equations for each representative Krep
nclasses=n_elements(p_conj)
E=Identity(3,/long)
indD=[0,1,2]

p_wyckoff=ptr_new(CodeType64(-1,-1))
sortmeasure=[0ull]

for i=0l,nclasses-1 do begin
    ; Conjugacy class representative Krep (take the first)
    Krep=elops[*(*p_conj[i])[0]]
    nKrep=n_elements(Krep)

    ; Skip Krep={E} (i.e. nKrep=1) because it gives the general position
    if (nKrep gt 1) then begin
        ; PreImage U of Krep in SG
        ; We can write U=HT and H will be the stabilizer
        ; of the union of affine subspaces found from Mx=b (mod Z)
        U=where(lops eq Krep[0])
        for j=1,nKrep-1 do U=[U,where(lops eq Krep[j])]
        nU=n_elements(U)
        if nU ne nKrep then stop
        U=allops[U]
        
        ; Find affine subspace which is "fixed modulo T" by U
        ; Set up the system of equations (E-Rj).x=Tj => Mx=b
        tmp=Symop(U[0])
        M=E-tmp.rot
        b=reform(tmp.trn,1,3)
        for j=1,nU-1 do begin
            tmp=Symop(U[j])
            M=[[M],[E-tmp.rot]]
            b=[[b],[reform(tmp.trn,1,3)]]
        endfor
        
        ; Solve k'th system mod Z:
        ;                M.xp = b (mod Z)
        ;          M.C^(-1).x = b (mod Z)
        ;                MM.x = b (mod Z)
        ; (P^(-1).D.Q^(-1)).x = b (mod Z)
        ;        D. (Q^(-1).x)= P.b (mod Z)
        ;             D.    x'= b2   (mod Z)
        if afftrans.trivial eq 0 then begin
            rational,M,/L64
            M=ratmatmult(M,afftrans.Mi)
            lcd=LCMmore(M.denom)
            M=M.num*(lcd/M.denom)
        endif else lcd=1
        if lcd ne 1 then stop

        D=SmithNormalForm(M,P=P,Q=Q)
        rational,P,/L64
        Float2Rational,b,24,/L64
        b=ratmatmult(P,b)
        Q*=lcd
        
        ; bi=0 (mod Z) for i>2
        if n_elements(b) ge 3 then $
            if total(b[3:*].denom ne 1,/int) ne 0 then continue ; no solution possible
        
        ; x=Q.x' and xi'={(b2i-Z)/Dii mod Z}
        b=b.num/float(b.denom)
        affsubspace=InequalModZ(b[indD],D[indD,indD],Q,nsol=nsol)
        if nsol eq 0 then continue

        ; Add the G-orbits of these subspaces in the previous setting
        ; Some of them might be already in the same G-orbit
        ; Keep subspaces with orbit size = [G:U]
        ; [G:U] = [G:T]/[U:T] = |Q|/|K|
        index=nQ*afftrans.ncen/nKrep ; x centering vectors because we will use allopsi to calulate the orbit
        
;        print,' '
;        print,'Lattice-equal subgroup U'+stringr(i)+': order U/T= '+stringr(nU)+' ; index [G:U]= ',stringr(index)+' ; '+(isGroup(U)?'group':'set')
;        print,'['+StringRTopcode(U)+']'
;        print,'Subspace fixed mod T by U'+stringr(i)+': ','['+StringRTopcode64(affsubspace)+']'

        AddGorbit,allopsi,affsubspace,p_wyckoff,wyckoffmult,sortmeasure,index
    endif
endfor
heap_free,p_conj

; Sort Wyckoff positions:
ind=bsort(sortmeasure)
p_wyckoff=p_wyckoff[ind]
wyckoffmult=wyckoffmult[ind]
if ~array_equal(wyckoffmult,wyckoffmult[sort(wyckoffmult)]) then stop ; should not happen

; General position
mgeneral=n_elements(allopsi)
general=replicate(CodeType64(0,0),mgeneral)
for i=1,mgeneral-1 do begin
    tmp=Symop(allopsi[i])
    general[i]=Symop_code64(tmp.rot,tmp.trn)
endfor

; Output in structure
nwyck=n_elements(p_wyckoff)-1
if nwyck ne 0 then begin

    ; make structure + skip dummy
    wyckpos={a:*(p_wyckoff[1])}
    for j=1,nwyck-1 do begin
        wyckpos=create_struct(wyckpos,wyckofflabel(j),*(p_wyckoff[j+1]))
    endfor
    wyckoffmult=wyckoffmult[1:*]

    wyckpos=create_struct(wyckpos,wyckofflabel(j),general)
    wyckoffmult=[wyckoffmult,mgeneral]
endif else begin
    wyckpos={a:general}
    wyckoffmult=[mgeneral]
endelse
heap_free,p_wyckoff
return,wyckpos

end;function wyckoffcodes
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SetASUwyckoff,sg,asupos,SmallDev

str=sgdata(sg.sghash,0)
wyck=call_function('wyckoff'+stringr(str.set),mult=mult)

nwyck=n_elements(mult)
tags=strlowcase(tag_names(wyck))
ilast=nwyck-1
code64last=CodeType64(0,0);(x,y,z)

; Spacegroup doesn't have special positions
if nwyck le 1 then begin
    for j=0l,n_elements(asupos)-1 do begin
        asupos[j].wyck.M=mult[ilast]
        asupos[j].wyck.L=tags[ilast]
        asupos[j].wyck.Rep=code64last
        asupos[j].wyck.xyz=asupos[j].xyz
    endfor
    return
endif

; Get all RT's from the wyckoff positions
rev=[0,total(mult,/cumul,/pres)]
RT=SymopAll64(wyck.(0))
for i=1,nwyck-2 do RT=[RT,SymopAll64(wyck.(i))]

; loop over ASU positions
for j=0l,n_elements(asupos)-1 do begin
    i=0

    ; loop over wyckoff positions
    repeat begin
        n=mult[i]
        k=rev[i]

        ; loop over the orbit of the wyckoff position
        repeat begin
            M=long((RT[k]).rot)
            b=transpose(asupos[j].xyz-(RT[k]).trn)
            D=SmithNormalForm(M,P=P,Q=Q)
            b=P##b
            D=D[[0,4,8]]

            ind=where(D eq 0,ct,comp=indc,ncomp=ctc)
            if ct ne 0 then $
                bjspecial=total(Equal0ModZ(b[ind],SmallDev),/pres) eq ct $
            else bjspecial=1b

            if bjspecial then begin
                xyz=[999.,999,999]
                if ctc ne 0 then begin
                    xyz[indc]=b[indc]/D[indc]
                    xyz=reform(Q##xyz)
                endif

                asupos[j].wyck.M=mult[i]
                asupos[j].wyck.L=tags[i]
                asupos[j].wyck.Rep=(wyck.(i))[k-rev[i]]
                asupos[j].wyck.xyz=xyz
                asupos[j].wyck.togglefix=1;xyz eq 999
            endif
            k++
        endrep until bjspecial or (k eq rev[i+1])
        i++
    endrep until bjspecial or (i eq nwyck-1)

    ; asupos[j] no special positions
    if ~bjspecial then begin
        asupos[j].wyck.M=mult[ilast]
        asupos[j].wyck.L=tags[ilast]
        asupos[j].wyck.Rep=code64last
        asupos[j].wyck.xyz=asupos[j].xyz
    endif
endfor

end;pro SetASUwyckoff
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function pgrp_ops,ops

lops=RotCode(ops)
elops=ExpandOpsList(lops)

return,elops[sort(elops)]
end;function pgrp_ops
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function LaueGroupHash,allops
return,hashCRC32(laue_ops(allops))
end;function LaueGroupHash
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PointGroupHash,allops
return,hashCRC32(pgrp_ops(allops))
end;function PointGroupHash
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Centrosymmetric,allops,origin=origin

mat_inv =[[-1.0, 0.0, 0.0], [0.0,-1.0, 0.0], [0.0, 0.0,-1.0] ]
vec_0=[ 0.0, 0.0, 0.0 ]

; Centrosymmetric: inversion center in the origin?
if keyword_set(origin) then ind=where(allops eq Symop_code(mat_inv,vec_0),ct) $
; Centrosymmetric: inversion center?
else ind=where(RotCode(allops) eq RotCode(Symop_code(mat_inv,vec_0)),ct)

return,ct ne 0

end;function Centrosymmetric
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function celparamdefault,sgstruc,add
system=strlowcase((pgdata(sgstruc.lghash,0)).crystsys)

case system of
'triclinic':    result=[5,6,7,60,70,80]
'monoclinic':    begin
                if strpos(sgstruc.choice,'c') ne -1 then result=[5,6,7,90,90,80] else $
                if strpos(sgstruc.choice,'a') ne -1 then result=[5,6,7,80,90,90] else $
                result=[5,6,7,90,80,90]
                endcase
'orthorhombic':    result=[5,6,7,90,90,90]
'tetragonal':    result=[5,5,6,90,90,90]
'trigonal':        if sgstruc.choice eq 'R' then result=[5,5,5,80,80,80] else result=[5,5,6,90,90,120]
'hexagonal':    result=[5,5,6,90,90,120]
'cubic':        result=[5,5,5,90,90,90]
else:            result=[5,6,7,60,70,80]
endcase

if n_elements(add) ne 0 then result[0:2]+=add
return,result
end;function celparamdefault
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function celparamconnect,sgstruc
system=strlowcase((pgdata(sgstruc.lghash,0)).crystsys)

case system of
'triclinic':    result=[0,1,2,3,4,5]
'monoclinic':    begin
                if strpos(sgstruc.choice,'c') ne -1 then result=[0,1,2,3,3,4] else $
                if strpos(sgstruc.choice,'a') ne -1 then result=[0,1,2,3,4,4] else $
                result=[0,1,2,3,4,3]
                endcase
'orthorhombic':    result=[0,1,2,3,3,3]
'tetragonal':    result=[0,0,2,3,3,3]
'trigonal':        if sgstruc.choice eq 'R' then result=[0,0,0,3,3,3] else result=[0,0,2,3,3,5]
'hexagonal':    result=[0,0,2,3,3,5]
'cubic':        result=[0,0,0,3,3,3]
else:            result=[0,1,2,3,4,5]
endcase

return,result
end;function celparamconnect
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function celparamfix,sgstruc
system=strlowcase((pgdata(sgstruc.lghash,0)).crystsys)

case system of
'triclinic':    result=[0,0,0,0,0,0]
'monoclinic':    begin
                if strpos(sgstruc.choice,'c') ne -1 then result=[0,0,0,1,1,0] else $
                if strpos(sgstruc.choice,'a') ne -1 then result=[0,0,0,0,1,1] else $
                result=[0,0,0,1,0,1]
                endcase
'orthorhombic':    result=[0,0,0,1,1,1]
'tetragonal':    result=[0,1,0,1,1,1]
'trigonal':        if sgstruc.choice eq 'R' then result=[0,1,1,0,1,1] else result=[0,1,0,1,1,1]
'hexagonal':    result=[0,1,0,1,1,1]
'cubic':        result=[0,1,1,1,1,1]
else:            result=[0,0,0,0,0,0]
endcase

return,result
end;function celparamfix
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function celparamsetcon,celparami,sgstruc,connect=connect

celparam=celparami
connect=celparamconnect(sgstruc)
ind=where(histogram(connect,REVERSE_INDICES=R) gt 1,ct)
for j=0l,ct-1 do begin
    tmp=R[R[ind[j]] : R[ind[j]+1]-1]
    celparam[tmp]=celparam[tmp[0]]
endfor

return,celparam
end;function celparamsetcon
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MetricTensor,param,reciprocal=reciprocal

a=param[0]
b=param[1]
c=param[2]
al=param[3]/180.*!pi
be=param[4]/180.*!pi
ga=param[5]/180.*!pi

M=[[a*a,a*b*cos(ga),a*c*cos(be)],$
[a*b*cos(ga),b*b,b*c*cos(al)],$
[a*c*cos(be),b*c*cos(al),c*c]]

if keyword_set(reciprocal) then M=invert(M)
return,M

end;function MetricTensor
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UCVolume,param,reciprocal=reciprocal,pder=pder,connect=connect
CATCH, Error_status
IF Error_status NE 0 THEN BEGIN
    return,0d
ENDIF else begin
    if not keyword_set(reciprocal) then begin
        a=param[0]
        b=param[1]
        c=param[2]
        cal=cos(param[3]/180.*!pi)
        cbe=cos(param[4]/180.*!pi)
        cga=cos(param[5]/180.*!pi)
        if arg_present(pder) and n_elements(connect) ne 0 then begin
            pder=fltarr(6)
            S=sqrt(1-cal*cal-cbe*cbe-cga*cga+2*cal*cbe*cga)
            pder[connect[0]]+=b*c*S
            pder[connect[1]]+=a*c*S
            pder[connect[2]]+=a*b*S
            V=a*b*c*S
            S=a*b*c/S
            pder[connect[3]]+=S*sin(param[3]/180.*!pi)*(cal-cbe*cga)/180.*!pi
            pder[connect[4]]+=S*sin(param[4]/180.*!pi)*(cbe-cal*cga)/180.*!pi
            pder[connect[5]]+=S*sin(param[5]/180.*!pi)*(cga-cal*cbe)/180.*!pi
            return,V
        endif else return,a*b*c*sqrt(1-cal*cal-cbe*cbe-cga*cga+2*cal*cbe*cga)
    endif else return,sqrt(determ(MetricTensor(param,reciprocal=reciprocal)))
endelse
end;function UCVolume
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function DirectToReciprocalU,param

; ----unit cell parameters----
aa=param[0]
bb=param[1]
cc=param[2]
p1=param[3]/180.*!dpi
p2=param[4]/180.*!dpi
p3=param[5]/180.*!dpi

; ----derived parameters(1)----
cp1=cos(p1)
cp2=cos(p2)
cp3=cos(p3)
sp1=sin(p1)
sp2=sin(p2)
sp3=sin(p3)
cq1=(cp2*cp3-cp1)/(sp2*sp3)
cq2=(cp3*cp1-cp2)/(sp3*sp1)
cq3=(cp1*cp2-cp3)/(sp1*sp2)

; ----derived parameters(2)----
S2=1.+2*cp1*cp2*cp3-cp1*cp1-cp2*cp2-cp3*cp3
Vdev=sqrt(S2)
;V=Vdev*aa*bb*cc

; ----reciprocal unit cell----
return,[sp1/(aa*Vdev),sp2/(bb*Vdev),sp3/(cc*Vdev),acos(cq1)*180./!pi,acos(cq2)*180./!pi,acos(cq3)*180./!pi]

end;function DirectToReciprocalU
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckEnantiomorphic,sghash,itcurrent=itcurrent
; Just 11 pairs in 3D so I just tabulated them
str=sgdata(sghash,0)
pairsIT=[[76,78],[91,95],[92,96],[144,145],[151,153],[152,154],[169,170],[171,172],[178,179],[180,181],[212,213]]
itcurrent=str.num
ind=where(pairsIT eq str.num,ct)
if ct ne 0 then return,pairsIT[1-(ind mod 2),ind/2] else return,0
end;function CheckEnantiomorphic
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function CheckSymmorphic,allops
; there are 73 symmorphic space group

; Primitive symmorphic?
if total(trncode(allops)) eq 0 then return,1b

; Non-primitive symmorphic?
trn=trncode(allops)
rot=rotcode(allops)
h=histogram(trn,reverse=R)
ind=where(h ne 0,ct)
if max(h[ind],min=mi) ne mi then return,0b

comp=rot[R[R[ind[0]] : R[ind[0]+1]-1]]
comp=comp[sort(comp)]
bsymm=1b
for j=1,ct-1 do begin
    comp2=rot[R[R[ind[j]] : R[ind[j]+1]-1]]
    comp2=comp2[sort(comp2)]
    bsymm and= array_equal(comp,comp2)
endfor

return,bsymm 
    
end;function CheckSymmorphic
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function BravaisType,HallLetter,CrystalSystem,CellChoice
; 14 types
case strlowcase(CrystalSystem) of
'triclinic':    return,'Triclinic'
'monoclinic':    return,(HallLetter eq 'P')?'Simple monoclinic':'Base-centered monoclinic'
'orthorhombic':    return,(HallLetter eq 'P')?'Simple orthorhombic':((HallLetter eq 'I')?'Body-centered orthorhombic':((HallLetter eq 'F')?'Face-centered orthorhombic':'Base-centered orthorhombic'))
'tetragonal':    return,(HallLetter eq 'P')?'Simple tetragonal':'Body-centered tetragonal'
'trigonal':        return,(CellChoice eq '')?'Hexagonal':'Rhombohedral'
'hexagonal':    return,'Hexagonal'
'cubic':        return,(HallLetter eq 'P')?'Simple cubic':((HallLetter eq 'I')?'Body-centered cubic':'Face-centered cubic')
else:            return,''
endcase
; trigonal + hexagonal crystal system = hexagonal family
end;function BravaisType
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SGInfo,sgstruc

str=sgdata(sgstruc.sghash,0)

if str.num eq 0 then begin
    str.hall=NTHallFromOps(sgstruc.allops)
    str.hm='Not tabulated.'
endif

data={it:str.num}

hm=str.hm
if str.ext ne '' then hm+=':'+str.ext
data=create_struct(data,'hm',hm)

data=create_struct(data,'sf',SGSFSymbol(str.num))
data=create_struct(data,'hall',str.hall)
data=create_struct(data,'set',str.set)
data=create_struct(data,'s_sg',string(str.num,format='(I0)')+': '+hm)
data=create_struct(data,'s_set',string(str.set,format='(I0)')+': '+hm)

lg=pgdata(sgstruc.lghash,0)
data=create_struct(data,'lauehm',lg.hm)
data=create_struct(data,'lauehall',lg.hall)
data=create_struct(data,'system',lg.crystsys)

tmp=str.choice
if tmp eq '' then tmp='none'
case strlowcase(data.system) of
'triclinic':
'monoclinic':    begin
                tmp=byte(tmp)
                ind=WHERE((tmp eq 49b) or (tmp eq 50b) or (tmp eq 51b), ct, COMPLEMENT=indc,NCOMPLEMENT=nct)
                tmp=['unique axis: '+((nct ne 0)?string(tmp[indc]):'none'),'cell choice: '+((ct ne 0)?string(tmp[ind]):'none')]
                endcase
'orthorhombic':    begin
                tmp=byte(tmp)
                ind=WHERE((tmp eq 49b) or (tmp eq 50b), ct, COMPLEMENT=indc,NCOMPLEMENT=nct)
                tmp=['origin choice: '+((ct ne 0)?string(tmp[ind]):'none'),'setting: '+((nct ne 0)?string(tmp[indc]):'none')]
                endcase
'tetragonal':    tmp='origin choice: '+tmp
'trigonal':        tmp='cell choice: '+tmp
'hexagonal':
'cubic':        tmp='origin choice: '+tmp
else:
endcase
data=create_struct(data,'choice',tmp)

tmp=['a','b','c','alpha','beta','gamma']
ind=where(celparamfix(sgstruc),ct)
if ct gt 1 then tmp[ind]=stringr((celparamdefault(sgstruc))[ind])
connect=celparamconnect(sgstruc)
ind=where(histogram(connect,REVERSE_INDICES=R) gt 1,ct)
for j=0l,ct-1 do begin
    tmp2=R[R[ind[j]] : R[ind[j]+1]-1]
    tmp[tmp2]=tmp[tmp2[0]]
endfor
data=create_struct(data,'unit',string(tmp,format='("[",5(A,", "),A,"]")'))

char=strmid(str.hall,0,1)
if char eq '-' then char=strmid(str.hall,1,1)
data=create_struct(data,'bravais',BravaisType(char,data.system,str.choice))

pg=pgdata(sgstruc.pghash,0)
data=create_struct(data,'pghm',pg.hm)
data=create_struct(data,'pgsf',pg.sfname)
data=create_struct(data,'pghall',pg.hall)

centro=Centrosymmetric(sgstruc.allops)
data=create_struct(data,'centro',((centro)?'yes':'no'))

if data.set ne 0 then begin
    wyck=call_function('wyckoff'+stringr(data.set),mult=wyckmult)
    data=create_struct(data,'wyckmult',wyckmult)
    wycknames=tag_names(wyck)
    for i=0l,n_elements(wyckmult)-1 do $
        data=create_struct(data,wycknames[i],StringRTopcode64(wyck.(i)))
endif

return,data
end;function SGInfo
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SpaceGroupAllops,symb,type,sghash=sghash

;    Find spacegroup specified by:
;  - 0. Hall symbol, e.g. " P 2ac 2ab"
;  - 1. H-M symbol, e.g. "P 21 21 21"
;  - 2. IT Number from 1-230
;  - 3. List of symmetry operators separated by semicolons, e.g.
;       "x,y,z;-x+1/2,-y,z+1/2;x+1/2,-y+1/2,-z;-x,y+1/2,-z+1/2"
;  - 4. Number from 1-530 (tabulated settings)
;  - 5. Schoenflies symbol
;  - 6. Hash

; Format strings
symb=symb[0]
if size(symb,/type) eq 7 then begin
    ; No leading or trailing blanks and no multiple spaces
    symbx=strtrim(strcompress(symb),2)

    ; Lowercase, except for first symbol
    p=strpos(symbx,' ')
    if p ne -1 then symbx=strmid(symbx,0,p)+strlowcase(strmid(symbx,p))
endif else symbx=symb

; Find ...
case type of
0:    ops=SGFromHall(symbx,error=error)
1:    ops=SGFromHM(symbx,error=error)
2:    ops=SGFromNr(symbx,error=error)
3:    ops=SGFromOp(symbx,error=error)
4:    ops=SGFromTab(symbx,error=error)
5:    ops=SGFromNr(SGSFSymbol(symbx),error=error)
6:    ops=SGFromHash(symbx,error=error)
endcase

if error then begin
    sghash=0UL
    return,[0L]
endif

; Make list of all symmetry operations
allops = ExpandOpsList(ops)
allops=allops[sort(allops)]

; Spacegroup code
sghash=hashCRC32(allops)

return,allops
end;function SpaceGroupAllops
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SpaceGroupStructure,sghash,allops
if sghash eq 0UL then return,{sghash:0UL,pghash:0UL,lghash:0UL,allops:[0L],choice:''}

; Point group code
;wyckhash=wyckoffcodes(allops)
pghash=PointGroupHash(allops)

; Laue class code
lghash=LaueGroupHash(allops)

; Get all info
info=sgdata(sghash,0)
choice=info.choice

; ---- Return spacegroup structure ----
; sghash: hash code for spacegroup
; pghash: hash code for pointgroup
; lghash: hash code for Laue group
; allops: codes for RT-ops
; choice: cel choice

return,{sghash:sghash,pghash:pghash,lghash:lghash,allops:allops,choice:choice}
end;function SpaceGroupStructure
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SpaceGroupStructure_event,ev,tree
case widget_info(ev.id,/type) of
8:    begin
    ID_system=widget_info(ev.top,find_by_uname='list_system')
    ID_pg=widget_info(ev.top,find_by_uname='list_pg')
    ID_sg=widget_info(ev.top,find_by_uname='list_sg')
    ID_set=widget_info(ev.top,find_by_uname='list_set')
    isystem=WIDGET_INFO(ID_system,/DROPLIST_SELECT)
    ipg=WIDGET_INFO(ID_pg,/DROPLIST_SELECT)
    isg=WIDGET_INFO(ID_sg,/DROPLIST_SELECT)
    iset=WIDGET_INFO(ID_set,/DROPLIST_SELECT)
    
    widget_control,ev.id,get_uvalue=uval
    case uval of
    'list_system':    begin
                ipg=0
                isg=0
                iset=0
                endcase
    'list_pg':    begin
                isg=0
                iset=0
                endcase
    'list_sg':    begin
                iset=isg ne 0
                endcase
    'list_set':    
    endcase

    widget_control,ID_system,set_value=tree.name
    p_pg=tree[isystem].ptr
    hash=0UL
    if ptr_valid(p_pg) then begin
        widget_control,ID_pg,set_value=(*p_pg).name
        p_sg=(*tree[isystem].ptr)[ipg].ptr
        if ptr_valid(p_sg) then begin
            widget_control,ID_sg,set_value=(*p_sg).name
            p_set=(*(*tree[isystem].ptr)[ipg].ptr)[isg].ptr
            if ptr_valid(p_set) then begin
                widget_control,ID_set,set_value=(*p_set).name
                hash=(*(*(*tree[isystem].ptr)[ipg].ptr)[isg].ptr)[iset].hash
            endif else widget_control,ID_set,set_value=''
        endif else begin
            widget_control,ID_sg,set_value=''
            widget_control,ID_set,set_value=''
        endelse
    endif else begin
        widget_control,ID_pg,set_value=''
        widget_control,ID_sg,set_value=''
        widget_control,ID_set,set_value=''
    endelse
    
    widget_control,ID_system,SET_DROPLIST_SELECT=isystem
    widget_control,ID_pg,SET_DROPLIST_SELECT=ipg
    widget_control,ID_sg,SET_DROPLIST_SELECT=isg
    widget_control,ID_set,SET_DROPLIST_SELECT=iset
    
    val=hash
    type=6
    endcase
3:    begin
    widget_control,ev.id,get_value=val,get_uvalue=uval
    case uval of
    'hall': type=0
    'hm':     type=1
    'it':     begin
            type=2
            val=fix(val)
            endcase
    'set':     type=4
    'sf':     type=5
    endcase
    endcase
else:     begin
        val=0UL
        type=6
        endelse
endcase

; Forward Function Definition because SpaceGroup also uses this function
FORWARD_FUNCTION SpaceGroup
return,SpaceGroup(val,type)
end;function SpaceGroupStructure_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SpacegroupDroplist,data,tree,vtree=vtree,vpg=vpg,vsg=vsg,vset=vset,init=initstruc,hash=hash

if not keyword_set(initstruc) then begin
    tmp=strsplit(data.s_sg,':',/extract)
    s_sg=tmp[0]
    ns_sg=strlen(s_sg)
    tmp=strsplit(data.s_set,':',/extract)
    s_set=tmp[0]
    ns_set=strlen(s_set)
    
    init=[data.system,data.pghm,s_sg,s_set]
    ind=where(init eq 'none',ct)
    if ct ne 0 then init[ind]=''
    
    sgadd=0
    setadd=0
endif else begin
    init=strarr(4)
    init[0]=initstruc.init
    ns_sg=strlen(init[2])
    ns_set=strlen(init[3])
    
    sgadd=initstruc.sgadd
    setadd=initstruc.setadd
endelse

indinit=lonarr(4)
vpg=''
vsg=''
vset=''
hash=0UL

vtree=tree.name
ind=where(tree.name eq init[0],ct)
if ct eq 1 then begin
    indinit[0]=ind[0]
    p=tree[ind[0]].ptr
    
    if ptr_valid(p) then begin
        vpg=(*p).name
        ind=where((*p).name eq init[1])
        indinit[1]=ind[0]>0
        p=(*p)[indinit[1]].ptr
            
        if ptr_valid(p) then begin
            vsg=(*p).name
            ind=where(strmid((*p).name,0,ns_sg) eq init[2])
            indinit[2]=(ind[0]>0)+sgadd
            p=(*p)[indinit[2]].ptr
            
            if keyword_set(initstruc) then initstruc.nsg=n_elements(vsg)-1
                
            if ptr_valid(p) then begin
                vset=(*p).name
                ind=where(strmid((*p).name,0,ns_set) eq init[3])
                indinit[3]=(ind[0]>0)+setadd
                if n_elements(vset) eq 2 then indinit[3]=1
                if indinit[3] ne 0 then hash=(*p)[indinit[3]].hash
                
                if keyword_set(initstruc) then initstruc.nset=n_elements(vset)-1
            endif
        endif
    endif
endif

return,indinit

end;function SpacegroupDroplist
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro typeTablestruc,base,tree,pstrinfo,update=update,nodropupdate=nodropupdate,init=init

noupdate=not keyword_set(update)

data={it:0, hm:'none', sf:'none', $
    hall:'none', set:0, s_sg:'none', s_set:'none', lauehm:'none', lauehall:'none', system:'none', bravais:'none', $
    pghm:'none',pghall:'none', pgsf:'none', centro:'none', wyckmult:[1],a:'none',unit:'none',choice:'none'}

case size(pstrinfo,/type) of
8:    if pstrinfo.sghash ne 0 then data=SGInfo(pstrinfo)
10: if ptr_valid(pstrinfo) then $
        if (*pstrinfo).sghash ne 0 then data=SGInfo(*pstrinfo)
else:
endcase

if noupdate then begin
    ; ----Delete all children----
    DeleteWBaseChilds,base
    
    ; ----Spacegroup droplists----
    base1=widget_base(base,/row)
    
    l=lonarr(4)
    
    baset=widget_base(base1,/column)
    label=widget_label(baset,value='Crystal system:')
    l[0]=widget_droplist(baset,uvalue='list_system',uname='list_system',xsize=100)
    
    baset=widget_base(base1,/column)
    label=widget_label(baset,value='Point group:')
    l[1]=widget_droplist(baset,uvalue='list_pg',uname='list_pg',xsize=100)
    
    baset=widget_base(base1,/column)
    label=widget_label(baset,value='Space group:')
    l[2]=widget_droplist(baset,uvalue='list_sg',uname='list_sg',xsize=100)
    
    baset=widget_base(base1,/column)
    label=widget_label(baset,value='Setting:')
    l[3]=widget_droplist(baset,uvalue='list_set',uname='list_set',xsize=100)

    ; ----Make editable info----
    basex=widget_base(base,/row)
    base1=widget_base(basex,/column)
    base2=widget_base(basex,/column)
    base3=widget_base(basex,/column)
    label=widget_label(base1,value='Hall symbol',/ALIGN_LEFT)
    t1=widget_text(base1,/ALIGN_LEFT,uvalue='hall',/editable,uname='hall')
    label=widget_label(base1,value='Hermann-Mauguin symbol',/ALIGN_LEFT)
    t2=widget_text(base1,/ALIGN_LEFT,uvalue='hm',/editable,uname='hm')
    label=widget_label(base1,value='IT number',/ALIGN_LEFT)
    t3=widget_text(base1,/ALIGN_LEFT,uvalue='it',/editable,uname='it')
    label=widget_label(base1,value='Setting number',/ALIGN_LEFT)
    t4=widget_text(base1,/ALIGN_LEFT,uvalue='set',/editable,uname='set')
    label=widget_label(base1,value='Schoenflies symbol',/ALIGN_LEFT)
    t5=widget_text(base1,/ALIGN_LEFT,uvalue='sf',/editable,uname='sf')
    ;label=widget_label(base2,value='Crystal system',/ALIGN_LEFT)
    ;t7=widget_text(base2,/ALIGN_LEFT,uname='system')
    label=widget_label(base2,value='Bravais Lattice',/ALIGN_LEFT)
    t8=widget_text(base2,/ALIGN_LEFT,uname='bravais')
    label=widget_label(base2,value='Unit Cell',/ALIGN_LEFT)
    t13=widget_text(base2,/ALIGN_LEFT,uname='unit')
    label=widget_label(base2,value='Setting choice',/ALIGN_LEFT)
    t14=widget_text(base2,/ALIGN_LEFT,uname='choice',ysize=2)
    label=widget_label(base2,value='Point group (Hall)',/ALIGN_LEFT)
    t9=widget_text(base2,/ALIGN_LEFT,uname='pghall')
;    label=widget_label(base2,value='Point group (Schoenflies)',/ALIGN_LEFT)
;    t10=widget_text(base2,/ALIGN_LEFT,uname='pgsf')
    label=widget_label(base2,value='Laue class (Hermann-Mauguin)',/ALIGN_LEFT)
    t6=widget_text(base2,/ALIGN_LEFT,uname='lauehm')
    label=widget_label(base2,value='Laue class (Hall)',/ALIGN_LEFT)
    t6b=widget_text(base2,/ALIGN_LEFT,uname='lauehall')
    label=widget_label(base2,value='Centrosymmetric',/ALIGN_LEFT)
    t11=widget_text(base2,/ALIGN_LEFT,uname='centro')
    label=widget_label(base3,value='Wyckoff positions',/ALIGN_LEFT)
    t12=widget_text(base3,/ALIGN_LEFT,uname='wyck',/SCROLL,scr_ysize=200)
endif else begin
    t1=widget_info(base,find_by_uname='hall')
    t2=widget_info(base,find_by_uname='hm')
    t3=widget_info(base,find_by_uname='it')
    t4=widget_info(base,find_by_uname='set')
    t5=widget_info(base,find_by_uname='sf')
    t6=widget_info(base,find_by_uname='lauehm')
    t6b=widget_info(base,find_by_uname='lauehall')
;    t7=widget_info(base,find_by_uname='system')
    t8=widget_info(base,find_by_uname='bravais')
    t9=widget_info(base,find_by_uname='pghall')
;    t10=widget_info(base,find_by_uname='pgsf')
    t11=widget_info(base,find_by_uname='centro')
    t12=widget_info(base,find_by_uname='wyck')
    t13=widget_info(base,find_by_uname='unit')
    t14=widget_info(base,find_by_uname='choice')
    
    l=lonarr(4)
    l[0]=widget_info(base,find_by_uname='list_system')
    l[1]=widget_info(base,find_by_uname='list_pg')
    l[2]=widget_info(base,find_by_uname='list_sg')
    l[3]=widget_info(base,find_by_uname='list_set')
endelse

if not keyword_set(nodropupdate) then begin
    indinit=SpacegroupDroplist(data,tree,vtree=vtree,vpg=vpg,vsg=vsg,vset=vset,init=init)
    widget_control,l[0],set_value=vtree
    widget_control,l[1],set_value=vpg
    widget_control,l[2],set_value=vsg
    widget_control,l[3],set_value=vset
    for i=0l,n_elements(l)-1 do widget_control,l[i],SET_DROPLIST_SELECT=indinit[i]
endif

widget_control,t1,set_value=data.hall
widget_control,t2,set_value=data.hm
widget_control,t3,set_value=stringr(data.it)
widget_control,t4,set_value=stringr(data.set)
widget_control,t5,set_value=data.sf
widget_control,t6,set_value=data.lauehm
widget_control,t6b,set_value=data.lauehall
;widget_control,t7,set_value=data.system
widget_control,t8,set_value=data.bravais
widget_control,t9,set_value=data.pghall
;widget_control,t10,set_value=data.pgsf
widget_control,t11,set_value=data.centro

tags=strlowcase(tag_names(data))
wyck=''
for i=n_elements(data.wyckmult)-1,0,-1 do begin
    j=where(tags eq wyckofflabel(i))
    wyck=[wyck,stringr(data.wyckmult[i])+tags[j]+':','['+data.(j)+']','']
endfor

widget_control,t12,set_value=wyck[1:*]
widget_control,t13,set_value=data.unit
widget_control,t14,set_value=data.choice

end;pro typeTablestruc
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro selectspacegroup_cleanup,ID

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

if WIDGET_INFO(ID,/VALID_ID) eq 0 then return
widget_control,ID,get_uvalue=list
heap_free,list.tree
if widget_info(list.top,/valid) then $
    widget_control,list.top,sensitive=1
end;pro selectspacegroup_cleanup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro selectspacegroup_event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

case widget_info(ev.id,/type) of
0:    begin
    IF TAG_NAMES(ev, /STRUCTURE_NAME) EQ 'WIDGET_KILL_REQUEST' THEN begin
        widget_control,ev.top,get_uvalue=list
        *list.pstrinfo=SpaceGroupStructure(0)
        WIDGET_CONTROL, ev.top, /DESTROY
    endif
    endcase
1:    begin
    widget_control,ev.id,get_value=val
    case val of
    'Accept':
    'Cancel':    begin
                widget_control,ev.top,get_uvalue=list
                *list.pstrinfo=SpaceGroupStructure(0)
                endcase
    endcase
    WIDGET_CONTROL, ev.top, /DESTROY
    endcase
else:begin
    widget_control,ev.top,get_uvalue=list
    strinfo=SpaceGroupStructure_event(ev,list.tree)
    error=strinfo.sghash eq 0
    
    if error then begin
        strinfo=*list.pstrinfo
    endif else begin
        *list.pstrinfo=strinfo
    endelse
    
    bdrop=widget_info(ev.id,/type) eq 8
    typeTablestruc,list.base,list.tree,strinfo,/update,nodropupdate=bdrop
    
    endelse
endcase
end;pro selectspacegroup_event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function selectspacegroup,top,init=init,sghash=sghash

if n_elements(top) eq 0 then top=0L else $
if widget_info(top,/valid) then WIDGET_CONTROL, top, sensitive=0

base=widget_base(/column,title='Space group',/TLB_KILL_REQUEST_EVENTS)

if n_elements(sghash) eq 0 then pstrinfo=ptr_new(SpaceGroupStructure(0)) $
else pstrinfo=ptr_new(SpaceGroupStructure(sghash,SpaceGroupAllops(sghash,6)))
list={top:top,pstrinfo:pstrinfo,base:base,tree:spacegrouptree()}

typeTablestruc,list.base,list.tree,list.pstrinfo,init=init

baset=widget_base(base,/row)
b=widget_button(baset,value='Cancel')
b=widget_button(baset,value='Accept')

WIDGET_CONTROL, base, SET_UVALUE=list, /REALIZE
Xmanager,'selectspacegroup',base,cleanup='selectspacegroup_cleanup',event_handler='selectspacegroup_event',GROUP_LEADER=top

sghash=(*pstrinfo).sghash
allops=(*pstrinfo).allops
ptr_free,pstrinfo
return,allops
end;function selectspacegroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function UnitCelltoRep,celparam

; Unit cell base vectors (normalized):
;
; ex'=[1,0,0]                => stays the same
; ey'=[cos(c),sin(c),0]        => coordinates of vector ey' in XY (Z coordinate=0)
; ez'=[u,v,w]                => unknown (with norm 1)
;
; M= 1  cos(c)  u
;     0  sin(c)  v
;     0    0     w
;
; cos(c)=ex'.ey'/(||ex'||.||ey'||)=cos(c)                => ok
; cos(b)=ex'.ez'/(||ex'||.||ez'||)=u                    => u=cos(b)
; cos(a)=ey'.ez'/(||ey'||.||ez'||)=ucos(c)+vsin(c)        => v=(cos(a)-cos(b)cos(c))/sin(c)
; |M|=V/(a.b.c)                                             => w=V/sin(c)
;
; => sqrt(u*u+v*v+w*w)=1

a=celparam[0]
b=celparam[1]
c=celparam[2]
alpha=celparam[3]
beta=celparam[4]
gamma=celparam[5]

alpha*=!dpi/180
beta*=!dpi/180
gamma*=!dpi/180
cosa=cos(alpha)
sina=sin(alpha)
cosb=cos(beta)
sinb=sin(beta)
cosc=cos(gamma)
sinc=sin(gamma)

V=sqrt(1d -cosa*cosa-cosb*cosb-cosc*cosc+2*cosa*cosb*cosc)

P=[ [1,cosc,cosb],$
    [0,sinc,(cosa-cosb*cosc)/sinc],$
    [0,   0,V/sinc]]

; Denormalize

P[0,*]*=a
P[1,*]*=b
P[2,*]*=c

return,P
end;function UnitCelltoRep
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ReptoUnitCell,P,Pref

; Triagonal matrix (Q: orthonormal: rotations and reflections)
M=transpose(P) ; because QRfact expects the transpose
QRfact,M,/returnR,/nopivot

; Make unique (positive diagonal elements)
s=[sign(M[0,0]),sign(M[1,1]),sign(M[2,2])]
if total(s eq 1,/pres) ne 3 then begin
    Pref=[[s[0],0,0],[0,s[1],0],[0,0,s[2]]]
    M=Pref##M
endif

; Extract celparameters
a=M[0,0]
b=sqrt(M[1,0]*M[1,0]+M[1,1]*M[1,1])
c=sqrt(M[2,0]*M[2,0]+M[2,1]*M[2,1]+M[2,2]*M[2,2])
    
gamma=atan(M[1,1],M[1,0])
cosc=cos(gamma)
sinc=sin(gamma)

beta=acos(M[2,0]/c)
alpha=acos((M[2,1]*sinc+M[2,0]*cosc)/c)

celparam=[a,b,c,alpha*180/!dpi,beta*180/!dpi,gamma*180/!dpi]

return,abs(celparam)
end;function ReptoUnitCell
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro TransformCelparam,celparam,Pconv
; Coordinates w.r.t. the standard Eucledian basis: X

; Coordinates w.r.t. Unit cell basis: X'
; X'=P^(-1).X
P=dblarr(4,4)
P[0,0]=UnitCelltoRep(celparam)
P[3,3]=1

; Coordinates w.r.t. Second unit cell basis: X''
; X''=Pconv^(-1).X'=Pconv^(-1).P^(-1).X

; Representation of this second basis
; X''=Q^(-1).X => Q=P.Pconv
P##=Pconv
celparam=ReptoUnitCell(P[0:2,0:2])

end;pro TransformCelparam
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function MatrixGroupsEqual,G1,G2,rational=rational,strict=strict
; Elements might be duplicated when strict=0

s1=dimsize(G1,3)
s2=dimsize(G2,3)
n1=s1[2]
n2=s2[2]

bcount=lonarr(n2)

if keyword_set(rational) then begin
    for i=0l,n1-1 do begin
        tmp=total(bcount,/int)
        for j=0,n2-1 do $
            bcount[j]+=ratarray_equal(G1[*,*,i],G2[*,*,j])
        if tmp eq total(bcount,/int) then return,0b ; G1[*,*,i] is not in G2
    endfor
endif else begin
    for i=0l,n1-1 do begin
        tmp=total(bcount,/int)
        for j=0,n2-1 do $
            bcount[j]+=array_equal(G1[*,*,i],G2[*,*,j])
        if tmp eq total(bcount,/int) then return,0b ; G1[*,*,i] is not in G2
    endfor
endelse

if keyword_set(strict) then return,total(bcount ne 1,/int) eq 0

if bcount[0] ne 0 then return,total(bcount ne bcount[0],/int) eq 0 $
else return,0b

end;function MatrixGroupsEqual
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SetofEquivIsometries,W0,W1,P,Pinv

; W0: source setting
; W1: destination setting
; W1 = Pinv ## W0 ## P set-wise
;
; P must be already checked for:
; - linear part in GL(3,Z)
; - linear part orientation preserving (det > 0)

n0=n_elements(W0)
n1=n_elements(W1)

; When P = E then the number of operations must be equal
btransform=n_elements(P) ne 0
if ~btransform and n0 ne n1 then return,0b

; Prepare destination setting
for l=0l,n1-1 do begin
    ; Translation positive and mod Z
    tmp=W1[l].trn
    RationalModZ,tmp
    
    ; Add to list
    W1l=MakeSeitzRat(W1[l].rot,tmp)
    if l eq 0 then W1_=W1l else W1_=[[[W1_]],[[W1l]]]
endfor

bcount=intarr(n1)
for k=0l,n0-1 do begin ; loop over the isometries in the source setting
    ; Get isometry after change of frame
;    tmp=W0[k].trn
;    RationalModZ,tmp
    W2=MakeSeitzRat(W0[k].rot,W0[k].trn)
    if btransform then W2=ratmatmult(Pinv,ratmatmult(W2,P))
    tmp=W2[3,0:2]
    RationalModZ,tmp
    W2[3,0:2]=tmp

    ; Mark corresponding isometries in the destination setting
    tmp=total(bcount,/int)
    for l=0,n1-1 do $
        bcount[l]+=ratarray_equal(W1_[*,*,l],W2)
    
    ; No corresponding isometry found for this 
    if total(bcount,/int) eq tmp then return,0b
endfor

; Check whether symmetry operations are the same (might be multiple occuring)
if bcount[0] ne 0 then return,total(bcount eq bcount[0],/int) eq n1 $
else return,0b
                
end;function SetofEquivIsometries
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SetofEquivIsometriesWrap,W0,W1,P,Pinv
if n_elements(W0) ge n_elements(W1) then return,SetofEquivIsometries(W0,W1,P,Pinv) $
    else return,SetofEquivIsometries(W1,W0,Pinv,P)
end;function SetofEquivIsometriesWrap
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SpaceGroupConvertGetChangeofBasis_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_uvalue=uval
widget_control,ev.top,get_uvalue=list
case widget_info(ev.id,/type) of
1:    begin
    case uval of
    'Cancel':    begin
                (*list).bcancel=1b
                widget_control,ev.top,/destroy
                endcase
    'OK':        begin
                (*list).bcancel=0b
                widget_control,ev.top,/destroy
                endcase
    'normal':    begin
                if ev.select then begin
                    C=(*list).P.num/float((*list).P.denom)
                    ID=widget_info(ev.top,FIND_BY_UNAME='tblC')
                    widget_control,ID,set_value=stringr(C[0:2,0:2])
                    ID=widget_info(ev.top,FIND_BY_UNAME='tblCt')
                    widget_control,ID,set_value=stringr(C[3,0:2])
                endif
                endcase
    'invert':    begin
                if ev.select then begin
                    C=(*list).Pinv.num/float((*list).Pinv.denom)
                    ID=widget_info(ev.top,FIND_BY_UNAME='tblC')
                    widget_control,ID,set_value=stringr(C[0:2,0:2])
                    ID=widget_info(ev.top,FIND_BY_UNAME='tblCt')
                    widget_control,ID,set_value=stringr(C[3,0:2])
                endif
                endcase
    'borig':    begin
                (*list).borig=~ev.select
                ID=widget_info(ev.top,FIND_BY_UNAME='tblCt')
                widget_control,ID,editable=(*list).borig
                endcase
    else:
    endcase
    endcase
9:    begin
    ID=widget_info(ev.top,FIND_BY_UNAME='normal')
    bnormal=widget_info(ID,/button_set)
    
    ; Table to temporary
    widget_control,ev.id,get_value=valstr
    val=float(valstr) ; denominators gone
    Float2Rational,val,1000,/thres,/L64
    for i=0,n_elements(val)-1 do begin
        tmp=strsplit(valstr[i],'/',count=n,/extract)
        if n eq 2 then begin
            val[i].num=long64(tmp[0])
            val[i].denom=long64(tmp[1])
        endif
    endfor
    SignNum,val
    ReduceNumDenom,val

    P=(*list).P
    Pinv=(*list).Pinv
    
    case uval of
    'tblC':    if bnormal then P[0:2,0:2]=val else Pinv[0:2,0:2]=val
    'tblCt':if bnormal then P[3,0:2]=val else Pinv[3,0:2]=val
    else:
    endcase
    
    if bnormal then begin
        Pinv=ratmatinvert(P,singular=singular)
    endif else begin
        P=ratmatinvert(Pinv,singular=singular)
    endelse

    ; Temporary to internal when valid
    det=ratmatdet(P)
    if singular or det.num lt 0 then begin
        tmp=dialog_message("Warning: "+(singular?"change of reference frame is singular.":"change of reference frame does not preserve the orientation."))
    endif
    (*list).P=P
    (*list).Pinv=Pinv
    
    ; Internal to table
    ID=widget_info(ev.top,FIND_BY_UNAME=uval)
    if bnormal then begin
        case uval of
        'tblC':    widget_control,ID,set_value=stringr((*list).P[0:2,0:2].num/float((*list).P[0:2,0:2].denom))
        'tblCt':widget_control,ID,set_value=stringr((*list).P[3,0:2].num/float((*list).P[3,0:2].denom))
        else:
        endcase
    endif else begin
        case uval of
        'tblC':    widget_control,ID,set_value=stringr((*list).Pinv[0:2,0:2].num/float((*list).Pinv[0:2,0:2].denom))
        'tblCt':widget_control,ID,set_value=stringr((*list).Pinv[3,0:2].num/float((*list).Pinv[3,0:2].denom))
        else:
        endcase
    endelse
    
    endcase
else:
endcase

end;pro SpaceGroupConvertGetChangeofBasis_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SpaceGroupConvertGetChangeofBasis,Cinv,Cshift,bcancel=bcancel,noprompt=noprompt

bcancel=0b

if keyword_set(noprompt) then begin
    ; For testing:
;    common GL3ZMatrix,C,Ct
;    Cinv=ratmatinvert(C)
;;    Cshift=Ct
;    return,C

    ; Normal:
    Cinv=identityrat(3,/L64)
    return,Cinv
endif

P=identityrat(4,/L64)
Pinv=P

list=ptr_new({P:P,Pinv:Pinv,bcancel:bcancel,borig:0b})

device,get_screen_size=screen
base=widget_base(title='Specify change of reference frame',/column,uvalue=list,xoffset=screen[0]/2.,yoffset=screen[1]/2.)

baset=widget_base(base,/exclusive,/column)
    b0=widget_button(baset,value='Coordinates of the conventional (listed) basis and origin as a function of the unconventional (given) frame.',uvalue='normal',uname='normal')
    b=widget_button(baset,value='Coordinates of the unconventional (given) basis and origin as a function of the conventional (listed) frame.',uvalue='invert',uname='invert')

baset=widget_base(base,/nonexclusive,/column)
    b1=widget_button(baset,value='Determine origin automatically.',uvalue='borig')

baset=widget_base(base,/row)
tblC=widget_table(baset,value=stringr((*list).P[0:2,0:2].num/float((*list).P[0:2,0:2].denom)),/editable,column_labels=['b1','b2','b3'],row_labels=['x','y','z'],uvalue='tblC',uname='tblC',xsize=3,ysize=3)
tblCt=widget_table(baset,value=stringr((*list).P[3,0:2].num/float((*list).P[3,0:2].denom)),editable=0,column_labels=['Origin'],row_labels=['x','y','z'],uvalue='tblCt',uname='tblCt',xsize=1,ysize=3)

b=widget_button(base,value='Cancel',uvalue='Cancel')
b=widget_button(base,value='OK',uvalue='OK')

widget_control,base,/realize
widget_control,b0,set_button=1
widget_control,b1,set_button=1
widget_control,b,/INPUT_FOCUS
Xmanager,'SpaceGroupConvertGetChangeofBasis',base, event_handler='SpaceGroupConvertGetChangeofBasis_Event'

C=(*list).P[0:2,0:2]
Cinv=(*list).Pinv[0:2,0:2]
if (*list).borig then Cshift=(*list).P[3,0:2]
bcancel=(*list).bcancel
PTR_FREE,list
return,C

end;function SpaceGroupConvertGetChangeofBasis
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro AcceptChangeofFrame_Event,ev

if checkrelease() then begin
    CATCH, Error_status
    IF Error_status NE 0 THEN begin
        printLastError
        return
    endif
endif

widget_control,ev.id,get_value=val
widget_control,ev.top,get_uvalue=list
case widget_info(ev.id,/type) of
1:    begin
    case val of
    'Decline':    begin
                (*list).bcancel=1b
                widget_control,ev.top,/destroy
                endcase
    'Accept':    begin
                (*list).bcancel=0b
                widget_control,ev.top,/destroy
                endcase
    else:
    endcase
    endcase
else:
endcase
end;pro AcceptChangeofFrame_Event
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function AcceptChangeofFrame,P,Pinv,noprompt=noprompt
if keyword_set(noprompt) then return,1b

list=ptr_new({bcancel:1b})

device,get_screen_size=screen
base=widget_base(title='Change of reference frame',/column,uvalue=list,xoffset=screen[0]/2.,yoffset=screen[1]/2.)

tmp=widget_label(base,value='Coordinates of the conventional (listed) basis and origin as a function of the unconventional (given) frame:')
baset=widget_base(base,/row)
tblC=widget_table(baset,value=stringr(P[0:2,0:2].num/float(P[0:2,0:2].denom)),/editable,column_labels=['b1','b2','b3'],row_labels=['x','y','z'],xsize=3,ysize=3)
tblCt=widget_table(baset,value=stringr(P[3,0:2].num/float(P[3,0:2].denom)),editable=0,column_labels=['Origin'],row_labels=['x','y','z'],xsize=1,ysize=3)

tmp=widget_label(base,value='Coordinates of the unconventional (given) basis and origin as a function of the conventional (listed) frame:')
baset=widget_base(base,/row)
tblC=widget_table(baset,value=stringr(Pinv[0:2,0:2].num/float(Pinv[0:2,0:2].denom)),/editable,column_labels=['b1','b2','b3'],row_labels=['x','y','z'],xsize=3,ysize=3)
tblCt=widget_table(baset,value=stringr(Pinv[3,0:2].num/float(Pinv[3,0:2].denom)),editable=0,column_labels=['Origin'],row_labels=['x','y','z'],xsize=1,ysize=3)

b=widget_button(base,value='Decline')
b=widget_button(base,value='Accept')

widget_control,base,/realize
widget_control,b,/INPUT_FOCUS
Xmanager,'AcceptChangeofFrame',base, event_handler='AcceptChangeofFrame_Event'

bcancel=(*list).bcancel
PTR_FREE,list
return,~bcancel
end;function AcceptChangeofFrame
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function FindOriginResolveZeroD,A_,B_,X
; A.X = B (mod Z)
; A -> rationals
; B -> rationals
; X -> integer solution

A=A_
B=B_

s=(DimSize(A,2))
nX=s[0]

m=LCMmore(A.denom)
A=A.num*(m/A.denom)
D=SmithNormalForm(A,P=P,Q=Q) ; P.A.Q = D
rational,P,/L64
rational,Q,/L64

; 1/m D.Q^(-1).X = P.B (mod Z)
;            D.Y = B   (mod Z)
B=ratmatmult(P,B)

; Rows with zero D's:
ind=where(total(D,1,/int) eq 0,ct)
if ct ne 0 then begin
    tmp=B[ind]
    RationalModZ,tmp
    ind2=where(tmp.num ne 0,ct2)
    if ct2 ne 0 then return,0b
endif

; D[k,k].x[k] = B'[k]    mod Z
X=make_ratarray(nX,/l64)
n=s[0]<s[1]
diag=lindgen(n)
D=D[diag,diag]
B=reform(B)
                    
ind=where(D eq 0,count)
if count ne 0 then D[ind]=1
rational,D
Y=ratdiv(B,D)
if count ne 0 then begin
    ; Could be any rational number in Z/m
    Y[ind].num=0
    Y[ind].denom=1
endif

; Y=1/m.Q^(-1).X
; X = Q.m.Y
Y.num*=m
ReduceNumDenom,Y
X[0:n-1]=Y
X=ratmatmult(Q,reform(X,1,nX))

if total(X.denom ne 1,/int) ne 0 then return,0b

return,1b

end;function FindOriginResolveZeroD
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SpaceGroupConvertFindOrigin,W0,T1,C,Ct

; T1+Z1 = C^(-1).(R0.Ct-Ct+T0+Z0)
; (R0-E).Ct = C.T1-T0+C.Z1-Z0
; A.Ct = B+H.Z1  (mod Z)

Erat=identityrat(3,/L64)

n0 = n_elements(W0)
n = 3
nequations = n*n0

A=make_ratarray([n,nequations],/L64)
B=make_ratarray([1,nequations],/L64)
H=make_ratarray([nequations,nequations],/L64)
for i=0,n0-1 do begin
    j0=i*n
    j1=j0+n-1
    A[*,j0:j1]=ratsub(W0[i].rot,Erat)
    B[*,j0:j1]=ratsub(ratmatmult(C,reform(T1[*,i],1,n)),reform(W0[i].trn,1,n))
    H[j0:j1,j0:j1]=C
endfor

; Make integer and get Smith normal form
m=LCMmore(A.denom)
A=A.num*(m/A.denom)
D=SmithNormalForm(A,P=P,Q=Q) ; P.A.Q = D
rational,P,/L64
rational,Q,/L64

; D/m.Q^(-1).Ct = P.B + P.H.Z1  (mod Z)
;           D.x = P.B + P.H.Z1  (mod Z)
B=ratmatmult(P,B)
H=ratmatmult(P,H)

; Rows with zero D's:
ind=where(total(D,1,/int) eq 0,ct)
if ct ne 0 then begin
    tmp=B[ind]
    RationalModZ,tmp
    ind2=where(tmp.num ne 0,ct2)
    if ct2 ne 0 then begin
        ; P.H.Z1 = -P.B (mod Z)
        ; HH.Z1 = BB    (mod Z) with Z1 integers
        HH=H[*,ind]
        BB=B[0,ind]
        ratneg,BB
        
        if ~FindOriginResolveZeroD(HH,BB,solZ) then return,0
        
        for i=0,nequations-1 do B=ratsum(B,ratmult(H[i,*],solZ[i]))
        
        tmp=B[ind]
        RationalModZ,tmp
        ind2=where(tmp.num ne 0,ct2)
        if ct2 ne 0 then return,0;stop
    endif
endif

; D[k,k].x[k] = B'[k]    mod Z
diag=(n+1)*lindgen(n)
D=D[diag]
B=B[0:n-1]
                    
ind=where(D eq 0,count)
if count ne 0 then D[ind]=1
rational,D
Ct=ratdiv(B,D) ; Could be any number (B+Z)/D
if count ne 0 then begin
    ; Could be any rational number but if it doesn't work for 0 it doesn't work for any rational number
    Ct[ind].num=0
    Ct[ind].denom=1
endif

;; Test real Ct in the set of solutions:
;x=ct
;common GL3ZMatrix,tmppp,ctdest
;xsol=ratmatmult(ratmatinvert(Q),reform(ctdest,1,3))
;xsol.denom*=m
;ReduceNumdenom,xsol
;if count ne 0 then x[ind]=xsol[ind]
;printrational,xsol
;printrational,x
;stop

; x=1/m.Q^(-1).Ct
; Ct = Q.m.x
Ct.num*=m
ReduceNumDenom,Ct
Ct=ratmatmult(Q,reform(Ct,1,3))

; Test:
; T1+Z1 = C^(-1).(R0.Ct-Ct+T0+Z0)
; C.T1 + C.Z1 = R0.Ct-Ct+T0+Z0
; C.T1 + C.Z1 = R0.Ct-Ct+T0 (mod Z)

;for i=0l,n0-1 do begin
;    check1=ratmatmult(C,reform(T1[*,i],1,n))
;    if n_elements(solZ) ne 0 then begin
;        i0=i*n
;        i1=i0+n-1
;        tmp=ratmatmult(C,reform(solZ[i0:i1],1,n))
;        check1=ratsum(check1,tmp)
;    endif
;    check2=ratsum(ratsub(ratmatmult(W0[i].rot,Ct),Ct),reform(W0[i].trn,1,n))
;    RationalModZ,check1
;    RationalModZ,check2
;    if ~ratarray_equal(check1,check2) then stop
;endfor

return,1                
end;function SpaceGroupConvertFindOrigin
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SplitGroupofIsometries,W,Tc,nTc,error=error
; Order appears to be important: fix this?

error=1b

; Find set of centering vectors
Erat=identityrat(3,/l64)
nW=n_elements(W)
nTc=0
for i=0l,nW-1 do begin
    if ratarray_equal(W[i].rot,Erat) then begin
        if nTc eq 0 then Tc=W[i].trn else Tc=[[Tc],[W[i].trn]]
        nTc++
    endif
endfor
if nTc eq 0 then return

; Centering vectors is only the zero vector
if nTc eq 1 then begin
    if total(Tc.num ne 0,/int) ne 0 then return
    
    Tc=reform(Tc,nTc,3)
    error=0b
    return
endif

if (nW mod nTc) ne 0 then return
Tc=reform(Tc,3,1,nTc)

; Keep only the operations with W=(q,tq+0)
indkeep=make_array(nW/nTc,value=-1l)
indkeepoff=0l
b=make_array(nW,value=1b)
for i=0l,nW-1 do begin
    if ~b[i] then continue
    
    ; Get operations with the same linear component
    ind=lonarr(nTc)
    indoff=0l
    for j=0,nW-1 do begin
        if ~b[j] then continue
        if ratarray_equal(W[i].rot,W[j].rot) then begin
            ind[indoff++]=j
            b[j]=0
        endif
    endfor
    if indoff ne nTc then return

    ; Get index of operation with centering vector 0
    for j=0,nTc-1 do begin
        Tcnew=W[ind].trn
        k=ind[j]
        Tsub=W[k].trn
        for jj=0,nTc-1 do begin
            kk=ind[jj]
            tmp=ratsub(Tcnew[*,jj],Tsub)
            RationalModZ,tmp
            Tcnew[*,jj]=tmp
        endfor
        Tcnew=reform(Tcnew,3,1,nTc)
        
        ; Tc and Tcnew are equal
        if MatrixGroupsEqual(Tc,Tcnew,/rational,/strict) then begin
            indkeep[indkeepoff++]=k
            break
        endif
    endfor
    
endfor

if total(indkeep eq -1,/int) ne 0 then return

W=W[indkeep]
Tc=ratmattranspose(Tc,[2,0,1])

error=0b
end;pro SplitGroupofIsometries
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro print1DLattice,coeff
;|d|.(Z + (y.indgen(m)+mx)/my)

d=coeff[0]
m=coeff[1]
x=coeff[2]
y=coeff[3]

set=lindgen(m)
rational,set,/L64
set.num*=d*y
set.num+=d*m*x
set.denom=m*y
ReduceNumDenom,set

str=stringr(set.num)+'/'+stringr(set.denom)
ind=where(set.num eq 0,ct)
if ct ne 0 then str[ind]=0

ind=where(set.denom eq 1,ct)
if ct ne 0 then str[ind]=stringr(set[ind].num)

if n_elements(str) gt 1 then str[0:n_elements(str)-2]+=','

print,d,str,format='(I0,"Z + {",'+stringr(m)+'(A),"}")'

end;pro print1DLattice
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Gen1DLatticefromcoeff,coeff
; |d|.Z + |d|.(y.indgen(m)+mx)/my

d=coeff[0]
m=coeff[1]
x=coeff[2]
y=coeff[3]
    
set=lindgen(m)
rational,set,/L64
set.num*=d*y
set.num+=d*m*x
set.denom=m*y

RationalModZ,set,d
set=set[sort(set.num/float(set.denom))]

return,set    
end;function Gen1DLatticefromcoeff
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function Simplify1DLattice,coeff,Tcout
;    d/m.Z + a/b = |d|.Z + |d|.(y.indgen(m)+mx)/my
;
;    The given set (i.e. a 1D lattice):
;   {d/m1.Z + a1/b1, d/m2.Z + a2/b2, ...} = |d|.Z + {|d|.(y.indgen(m)+mx)/my, ...}
;   
;   Required lattice must have the form:
;   Z + {tc1, tc2,...} = dZ + {0,1,...,d-1} + {tc1, tc2,...}
;
;    Therefore:
;   {|d|.(y.indgen(m)+mx)/my, ...} = {0,1,...,d-1} + {tc1, tc2,...}

; d must divide the total number of centerings
nTc=total(coeff[1,*],/int) ; total number of centerings
d=coeff[0,0] ; all the same d so pick the first
if (nTc mod d) ne 0 then return,0

; List of all centerings
Tc=make_ratarray(nTc,/L64)
off=0l
s=dimsize(coeff)
nspace=s[1]
for i=0,nspace-1 do begin
    set=Gen1DLatticefromcoeff(coeff[*,i])

    m=coeff[1,i]
    Tc[off:off+m-1].num=set.num
    Tc[off:off+m-1].denom=set.denom
    
    off+=m
endfor
ReduceNumDenom,Tc

; Check Tc = {0,1,...,d-1} + {tc1, tc2,...}
bfix=bytarr(nTc)
check1=lindgen(d)
for i=0l,nTc-1 do begin
    Tcsubi=ratsub(Tc,Tc[i])
    
    ; Look for integers
    ind=where(Tcsubi.denom eq 1,ct)
    if ct ne d then return,0
    
    ; Are they used already
    nfix=total(bfix[ind],/int)
    if nfix eq d then continue
    if nfix ne 0 then return,0
    
    ; Check whether these are 
    check2=Tcsubi[ind].num mod d
    ind2=where(check2 lt 0,ct)
    if ct ne 0 then check2[ind2]+=d

    ; Number of integers
    if ~array_equal(check2[UNIQ(check2, SORT(check2))], check1) then return,0
    
    ; Add Tc to list
    if n_elements(Tcout) eq 0 then Tcout=Tc[i] else Tcout=[Tcout,Tc[i]]
    
    bfix[ind]=1b
endfor

if total(bfix,/int) ne nTc then return,0

; Map between 0 and 1 and remove duplicate zero's
RationalModZ,Tcout
Tcout=Tcout[sort(Tcout.num/float(Tcout.denom))]

ind=where(Tcout.num ne 0,ct,comp=indzero,ncomp=nzero)

if nzero eq 0 then return,ct ; no zero's

if ct eq 0 then begin
    ; Only zero's
    nTcout=1
    Tcout=Tcout[indzero[0]]
endif else begin
    nTcout=ct+1
    Tcout=[Tcout[indzero[0]],Tcout[ind]]
endelse

return,nTcout

end;function Simplify1DLattice
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro SimplifyCenteringVectors,Tc,nTc

; mod Z and between 0 and 1
RationalModZ,Tc

; Remove duplicates
bchecked=bytarr(nTc)
buse=bchecked
for i=0l,nTc-1 do begin
    if bchecked[i] then continue
    
    ; Use this vector
    buse[i]=1b
    bchecked[i]=1b
    
    ; Mark duplicates as checked
    for j=0l,nTc-1 do begin
        if bchecked[j] then continue
        bchecked[j]=ratarray_equal(Tc[i,*],Tc[j,*])
    endfor
endfor

ind=where(buse,nTc)
Tc=Tc[ind,*]

end;pro SimplifyCenteringVectors
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lattice1Dequal,lattice1,lattice2,censpecial

if lattice1[0] gt 1 and lattice2[0] gt 1 and lattice1[0] ne lattice2[0] then return,0

if lattice1[0] gt 0 then set1=Gen1DLatticefromcoeff(lattice1) $
else set1=*censpecial[-lattice1[0]]

if lattice2[0] gt 0 then set2=Gen1DLatticefromcoeff(lattice2) $
else set2=*censpecial[-lattice2[0]]

return,ratarray_equal(set1,set2)
        
end;function lattice1Dequal
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pro CenteringVectorsUnderChangeofBasis,Tc0,nTc0,Cinv,error=error
; Tc1 + Z^3 = C^(-1).(Tc0 + Z^3)

error=1b
n=3
diag=(n+1)*lindgen(n)

; Least common denominator of C^(-1)
m=LCMmore(Cinv.denom)

; C^(-1)=A/m with A an integer matrix
A=Cinv.num*(m/Cinv.denom)
D=SmithNormalForm(A,P=P,Q=Q)

if total(D[diag] eq 0,/int) ne 0 then return ; Maybe no need to check this because Cinv in GL(3,Q)

Pinv=P
rational,Pinv,/L64
Pinv=ratmatinvert(Pinv,singular=singular)
if singular or total(Pinv.denom ne 1,/int) ne 0 then return;stop

Qinv=Q
rational,Qinv,/L64
Qinv=ratmatinvert(Qinv,singular=singular)
if singular or total(Qinv.denom ne 1,/int) ne 0 then return;stop

rational,D,/L64
D[diag].denom=m
ReduceNumDenom,D

; Vc = D/m.Q^(-1).Tc0
Vc=ratmatmult(D,ratmatmult(Qinv,Tc0))
nVc=nTc0

; Vc + D/m.Z^3 = {Vc0 + D/m.Z^3, Vc1 + D/m.Z^3, ...}
;
; Multiplication of Z with a rational scalar:
;   d/m.Z = |d|/m.Z   (sign absorbed by Z)
;         = |d|.(Z + {0,1,...,(m-1)}/m)
;         = |d|.(Z + indgen(m)/m)
;         
; Adding a rational:
;   d/m.Z + a/b = |d|.(Z + indgen(m)/m + a/b|d|) 
;               = |d|.(Z + indgen(m)/m + x/y)    x/y = a/b|d| + z so that 0<x/y<1   (z absorbed by Z)
;               = |d|.(Z + (y.indgen(m)+mx)/my)

coeff=lonarr(4,n,nVc) ; d,m,x,y
for i=0l,nVc-1 do begin ; loop over the centering vectors
    for j=0l,n-1 do begin ; loop over their coordinates
        ; Construct 1D lattice for Vci and coordinate j: 
        ;  |d|.(Z + (y.indgen(m)+mx)/my)
        ; coeff[*,j,i] = [d,m,x,y]
        
        d_=abs(D[j,j].num)
        coeff[0,j,i]=d_ ; d
        coeff[1,j,i]=D[j,j].denom ; m
        
        tmp=Vc[i,j] ; a/b
        tmp.denom*=d_ ; a/b|d|
        RationalModZ,tmp ; x/y
        
        coeff[2,j,i]=tmp.num ; x
        coeff[3,j,i]=tmp.denom ; y
        
;        print1Dlattice,coeff[*,j,i]
    endfor
;    print,''
endfor

; Eliminate 1D lattices with multiples of Z (d != 1) if possible
censpecial=ptr_new()
for j=0l,n-1 do begin ; Loop over the dimension (3 in this case)
    ; coeff[0,j,*] == coeff[0,j,0]
    if coeff[0,j,0] ne 1 then begin
        ; Complementary spaces
        ind=ShrinkArray(indgen(n),j)
        compspace=coeff[*,ind,*]
        
        ; Group spaces that have the same complementary spaces together
        equalspace=bytarr(nVc,nVc)
        for i1=0,nVc-1 do $
            for i2=i1,nVc-1 do begin
                bequal=1b
                
                for itmp=0,n-2 do $
                    bequal and= lattice1Dequal(compspace[*,itmp,i1],compspace[*,itmp,i2],censpecial)

                equalspace[i1,i2]=bequal
                equalspace[i2,i1]=bequal
            endfor
        
        ; Each group must contain the same number (more than 1) of spaces
        tmp1=total(equalspace,1)
        tmp2=total(equalspace,2)
        if ~array_equal(tmp1,tmp2) then return
        ncompingroup=tmp1[0]
        if total(tmp1 ne ncompingroup,/int) ne 0 then return
        if (nVc mod ncompingroup) ne 0 then return
        if ncompingroup eq 1 then return
        
        ; Loop over the groups and find out whether each group can be written as Z + {...}
        ; and this must be the same for each group
        bdone=bytarr(nVc)
        bkeep=bytarr(nVc)
        for i=0l,nVc-1 do begin
            ; Get group: {d/m.Z + a/b, ...}
            if bdone[i] then continue
            ind=where(equalspace[i,*],nind)
            group=reform(coeff[*,j,ind])

            ; Write group as Z + {...}  (if possible)
            delvar2,cen
            ncen=Simplify1DLattice(group,cen)
            if ncen eq 0 then return ; should be at least one centering (zero)

            ; Mark spaces
            coeff[*,*,ind[1:*]]=0 ; mark for removal
            coeff[0,j,ind[0]]=-n_elements(censpecial) ; mark for special
            
            ; Add to list of special centerings
            censpecial=[censpecial,ptr_new(cen)]

            bdone[ind]=1b
            bkeep[ind[0]]=1b
        endfor
        ind=where(bkeep,nVc)
        coeff=coeff[*,*,ind]
    endif
endfor

;printrational,Vc
;printrational,D
;printrational,Pinv
;printrational,Qinv

; Make centering vectors
for i=0l,nVc-1 do begin ; loop over the centering vectors
    
    ; Each coordinate coeff[*,j,i] defines a set Z + {...}
    ; so make all combinations to get the centering vectors
    
    for j=0l,n-1 do begin
        if coeff[0,j,i] eq 1 then set=Gen1DLatticefromcoeff(coeff[*,j,i]) $
        else set=*censpecial[-coeff[0,j,i]]
        
        if j eq 0 then setadd=set else begin
            tmp=setadd
            ntmp=(DimSize(tmp,2))[0]
            
            type=size(set[0].num,/type)
            setadd=[[tmp],[make_ratarray(ntmp,value=set[0],type=type)]]
            for k=1,n_elements(set)-1 do setadd=[setadd,[[tmp],[make_ratarray(ntmp,value=set[k],type=type)]]]
        endelse
    endfor
    
    if i eq 0 then Vc=setadd else Vc=[Vc,setadd]
    
endfor
heap_free,censpecial

;printrational,Vc

Tc0=ratmatmult(Pinv,Vc)
nTc0=(DimSize(Tc0,2))[0]
SimplifyCenteringVectors,Tc0,nTc0

error=0b

; Rational linear combinations of Z:
;     Linear combination with integer coefficients:
;         a.Z + b.Z = GCD(a,b).Z  => Bezout's lemma
;     Multiplication with a rational scalar:
;         a/b.Z = a.Z + {0 <= za/b < a: z in Z}
;     Linear combination with rational coefficients:
;         a/b.Z + c/d.Z = (x.Z + y.Z)/m   with m=LCM(b,d)  (positive)
;                       = k/m.Z           with k=GCD(x,y)  (positive)
;                       = kZ + {0,k/m,...,(m-1)k/m}

end;pro CenteringVectorsUnderChangeofBasis
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SpaceGroupConvert,sghash0,allops0,strops=strops,error=error,noprompt=noprompt
; A change of frame in Euclidean space has the following effect on the isometries W0 in a space group
;     W1 = P^(-1).W0.P
;         W0 = [[R0,T0],[0,1]]
;         W1 = [[R1,T1],[0,1]]
;         P  = [[C,Ct],[0,1]]
;         P^(-1)= [[C^(-1),-C^(-1).Ct],[0,1]]
; 
; We know W0 and we will go through all W1 (conventional space group settings) and try to
; find the change of frame matrix P which transforms all W0's in W1's.
; We can split this problem in a linear an translational component
;     R1 = C^(-1).R0.C
;     T1 = C^(-1).(R0.Ct-Ct+T0)

error=1b

; Select manual or not
if keyword_set(noprompt) then bauto=1b else $
    bauto=dialog_message("Space group settings needs conversion. Find tabulated setting automatically?",/question) eq 'Yes'
if bauto then tree=spacegrouptree()

; Get unconventional symmetry operations
if keyword_set(strops) then W0=RTopStringRatAll(strops,error=error,/L64) $
else W0=SymopAllRat(strops,error=error,/L64)
W0keep=W0

if error then begin
    if ~keyword_set(noprompt) then $
        tmp=dialog_message("Error in parsing symmetry operations of unconventional space group setting.")
    return,identity(4,/float)
endif

; Remove the centering vectors
SplitGroupofIsometries,W0,Tc0,nTc0,error=error
if error then begin
    if ~keyword_set(noprompt) then $
        tmp=dialog_message("Error in the set of centering vectors of unconventional space group setting.")
    return,identity(4,/float)
endif
n0=n_elements(W0)

; Get a change of basis matrix C
repeat begin
    C=SpaceGroupConvertGetChangeofBasis(Cinv,Cshift,bcancel=bcancel,noprompt=noprompt)
    if bcancel then begin
        error=1b
        return,identity(4,/float)
    endif
    
    Erat=identityrat(3,/L64)
    btransf=~ratarray_equal(Erat,C)
    
    ; Transform centering vectors
    Tc1=Tc0
    nTc1=nTc0
    bstop=1b
    if btransf then begin
        CenteringVectorsUnderChangeofBasis,Tc1,nTc1,Cinv,error=error
        if error then begin
            bstop=0b
            if ~keyword_set(noprompt) then $
                tmp=dialog_message("Change of basis does not allow for a conventional lattice representation.")
        endif
    endif
    Tc1=reform(ratmattranspose(Tc1),3,1,nTc1)
endrep until bstop

; Transform point group
R1=W0.rot
R1diag=(n0+1)*lindgen(n0)
for i=0,n0-1 do begin
    if btransf then R1[*,*,i]=ratmatmult(Cinv,ratmatmult(R1[*,*,i],C))
    tmp=ratmatinvert(R1[*,*,i],singular=singular)

    if total(R1[*,*,i].denom ne 1,/int) ne 0 or total(tmp.denom ne 1,/int) ne 0 or singular then begin
        if ~keyword_set(noprompt) then $
            tmp=dialog_message("Point group not a subgroup of GL(3,Z).")
        return,identity(4,/float)
    endif
endfor

; Defnitions:
; W0keep: unconventional setting
; W0, Tc0: splitting of W0keep
; n0, nTc0: number of W0's and Tc0's
; R1: linear parts of the conventional setting (derived from the linear parts of W0)
;         => n0 number of R1's
; Tc1: centering vectors of the conventional setting (derived from Tc0)
; nTc1: number of Tc1's

; Get point group op codes
tmp=fltarr(3)
pgcodes=lonarr(n0)
for i=0,n0-1 do pgcodes[i]=Symop_code(R1[*,*,i].num,tmp)
pgcodes=ExpandOpsList(pgcodes)

; Narrow down the space group search
pg=pgdata(PointGroupHash(pgcodes),0)
init={init:[pg.crystsys,pg.hm],sgadd:1,setadd:1,nsg:0,nset:0}
btryall=pg.crystsys eq '' or pg.hm eq ''
if btryall then init.nsg=530

; Find conventional setting
repeat begin
    ; Select corresponding conventional setting
    if bauto then begin
        if btryall then allops1=SpaceGroupAllops(init.sgadd,4,sghash=sghash1) else begin
            indinit=SpacegroupDroplist(0,tree,init=init,hash=hash)
            allops1=SpaceGroupAllops(hash,6,sghash=sghash1)
;            print,''
;            print,sgdata(sghash1,0)
        endelse
    endif else allops1=selectspacegroup(init=init,sghash=sghash1)

    ; Space group candidate
    W1=SymopAllRat(allops1,/L64)
    W1keep=W1
    
    ; Remove the centering vectors
    SplitGroupofIsometries,W1,Tc1_,nTc1_,error=error
    s1=Dimsize(W1.rot,3)
    n1=s1[2]
    Tc1_=reform(ratmattranspose(Tc1_),3,1,nTc1_)
    
    ; Defnitions:
    ; W1keep: conventional setting
    ; W1, Tc1_: splitting of W1keep
    ; n1, nTc1_: number of W1's and Tc1_'s
    
    bequiv=~error and nTc1_ eq nTc1 and n0 eq n1
    if bequiv then begin
        ; Compare centering vectors
        bequiv=MatrixGroupsEqual(Tc1,Tc1_,/rational,/strict)
        
        ; Compare point groups
        if bequiv then bequiv=MatrixGroupsEqual(R1,W1.rot,/rational,/strict)
        
        ; Try to find an origin
        if bequiv then begin
            ; Change of frame: insert change of basis
            P = identityrat(4,/l64)
            P[0:2,0:2] = C
            
            ; Find the corresponding R0<->R1
            ind=lindgen(n0)
            for i=0,n0-1 do $
                for j=0,n0-1 do $
                    if ratarray_equal(R1[*,*,i],W1[j].rot) then ind[i]=j
            if n_elements(UNIQ(ind, SORT(ind))) ne n0 then bequiv=0 ;stop
                
            ; Try to find an origin shift so that W0 and W1 are equivalent
            nCt=0
            if bequiv then begin
            if n_elements(Cshift) ne 0 then begin
                nCt=1
                Ct=Cshift
            endif else nCt=SpaceGroupConvertFindOrigin(W0,W1[ind].trn,C,Ct)
            endif
            if nCt ne 0 then begin
;                printrational,Ct
;                print,''
                for i=0l,nCt-1 do begin
                    P[3,0:2] = Ct[i,*]
;                    printrational,P
                    Pinv=ratmatinvert(P)
                    bequiv=SetofEquivIsometriesWrap(W0keep,W1keep,P,Pinv)
                    if bequiv then break
                endfor
            endif else bequiv=0
        endif
    endif
    
    ; Did we find the space group setting?
    if bauto then begin
        bstop=bequiv or (sghash1 eq 0) or ((init.sgadd eq init.nsg) and (init.setadd eq init.nset))
        if init.setadd lt init.nset then init.setadd++ else begin
            init.sgadd++
            init.setadd=1
        endelse
        
        if bequiv then bequiv=AcceptChangeofFrame(P,Pinv,noprompt=noprompt)
        
    endif else begin
        bstop=1b
        if ~bequiv and ~keyword_set(noprompt) then bstop=dialog_message("Space group settings don't correspond. Try again?",/question) eq 'No'
    endelse
    
endrep until bstop ; loop over the conventional settings

if bauto then heap_free,tree

error=~bequiv
if bequiv then begin
    P=P.num/float(P.denom)
    ; Replace the old basis by the new basis
    sghash0=sghash1
    allops0=allops1
endif else begin
    P=identity(4,/float)
    if ~keyword_set(noprompt) then $
        tmp=dialog_message("No space group setting found.")
endelse

return,P

end;function SpaceGroupConvert
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function SpaceGroup,symb,type,error=error,Pconv=Pconv,noprompt=noprompt

; Pconv: returns the 4x4 convertion matrix from unconventional to conventional setting
;          transform cell parameters: use TransformCelparam
;         transform atomic positions: invert(Pconv)##xyz

Pconv=identity(4,/float)

; Look for conventional setting
;if type eq 3 then sghash=0ull else allops=SpaceGroupAllops(symb,type,sghash=sghash)
allops=SpaceGroupAllops(symb,type,sghash=sghash)
tmp=sgdata(sghash,0)
error=tmp.set eq 0

; Check whether the space group corresponds to the list of symmetry operations given
bequiv=~error
if bequiv and type eq 3 then begin
    W0=RTopStringRatAll(symb,error=error,/L64)
    if ~error then begin
        W1=SymopAllRat(allops,error=error,/L64)
        if error then bequiv=0b else $
            if n_elements(W0) ne n_elements(W1) then bequiv=0b else $
                bequiv=SetofEquivIsometries(W0,W1)
    endif
endif

; Unconventional setting?
if (type eq 0 or type eq 3) and ~bequiv then begin
    Pconv=SpaceGroupConvert(sghash,allops,strops=(type eq 3)?symb:0,error=error,noprompt=noprompt)
    ; sghash is the conventional space group setting hash (unless error == 1)
    if error then sghash=0
endif

return,SpaceGroupStructure(sghash,allops)
end;function SpaceGroup
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
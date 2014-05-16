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

function Spacegrouptree

p0=ptr_new()
q0=ptr_new()
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'1: P 1',hash:'C704DD7B'XL}])
q1=ptr_new([ $
    {name:'',ptr:r0},{name:'1: P 1',ptr:r1}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'2: P -1',hash:'2DF60A45'XL}])
q2=ptr_new([ $
    {name:'',ptr:r0},{name:'2: P -1',ptr:r1}])
p1=ptr_new([ $
    {name:'',ptr:q0},{name:'1',ptr:q1},{name:'-1',ptr:q2}])
q0=ptr_new()
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'3: P 1 2 1',hash:'902CF668'XL},{name:'4: P 1 1 2',hash:'90A34743'XL},{name:'5: P 2 1 1',hash:'9B02BFBF'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'6: P 1 21 1',hash:'DA168154'XL},{name:'7: P 1 1 21',hash:'7851C0D1'XL},{name:'8: P 21 1 1',hash:'AE0E24DB'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'9: C 1 2 1',hash:'DEB5173B'XL},{name:'10: A 1 2 1',hash:'1E4DBCA6'XL},{name:'11: I 1 2 1',hash:'6042ABB6'XL}, $
    {name:'12: A 1 1 2',hash:'2A55A450'XL},{name:'13: B 1 1 2',hash:'458C0B8C'XL},{name:'14: I 1 1 2',hash:'545AB340'XL},{name:'15: B 2 1 1',hash:'C4EF9A2D'XL}, $
    {name:'16: C 2 1 1',hash:'6BCE9E6C'XL},{name:'17: I 2 1 1',hash:'D53922E1'XL}])
q1=ptr_new([ $
    {name:'',ptr:r0},{name:'3: P 1 2 1',ptr:r1},{name:'4: P 1 21 1',ptr:r2},{name:'5: C 1 2 1',ptr:r3}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'18: P 1 m 1',hash:'93121E4D'XL},{name:'19: P 1 1 m',hash:'939DAF66'XL},{name:'20: P m 1 1',hash:'6DD7504B'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'21: P 1 c 1',hash:'7BE099DF'XL},{name:'22: P 1 n 1',hash:'4EEC02BB'XL},{name:'23: P 1 a 1',hash:'A61E8529'XL}, $
    {name:'24: P 1 1 a',hash:'A6913402'XL},{name:'25: P 1 1 n',hash:'ECAB433E'XL},{name:'26: P 1 1 b',hash:'D9A7D85A'XL},{name:'27: P b 1 1',hash:'27ED2777'XL}, $
    {name:'28: P n 1 1',hash:'2AC91F43'XL},{name:'29: P c 1 1',hash:'8525D7D9'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'30: C 1 m 1',hash:'77B286E1'XL},{name:'31: A 1 m 1',hash:'B74A2D7C'XL},{name:'32: I 1 m 1',hash:'C9453A6C'XL}, $
    {name:'33: A 1 1 m',hash:'8352358A'XL},{name:'34: B 1 1 m',hash:'EC8B9A56'XL},{name:'35: I 1 1 m',hash:'FD5D229A'XL},{name:'36: B m 1 1',hash:'DDC257BB'XL}, $
    {name:'37: C m 1 1',hash:'72E353FA'XL},{name:'38: I m 1 1',hash:'CC14EF77'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'39: C 1 c 1',hash:'38EF08AA'XL},{name:'40: A 1 n 1',hash:'54FAF6D0'XL},{name:'41: I 1 a 1',hash:'9D54258D'XL}, $
    {name:'42: A 1 a 1',hash:'E35B329D'XL},{name:'43: C 1 n 1',hash:'6CFE174B'XL},{name:'44: I 1 c 1',hash:'2AF5E1C0'XL},{name:'45: A 1 1 a',hash:'D7432A6B'XL}, $
    {name:'46: B 1 1 n',hash:'0F3B41FA'XL},{name:'47: I 1 1 b',hash:'4AFCE6D7'XL},{name:'48: B 1 1 b',hash:'5B2A5E1B'XL},{name:'49: A 1 1 n',hash:'60E2EE26'XL}, $
    {name:'50: I 1 1 a',hash:'A94C3D7B'XL},{name:'51: B b 1 1',hash:'6A6393F6'XL},{name:'52: C n 1 1',hash:'69AFC250'XL},{name:'53: I c 1 1',hash:'2FA434DB'XL}, $
    {name:'54: C c 1 1',hash:'3DBEDDB1'XL},{name:'55: B n 1 1',hash:'3E728C17'XL},{name:'56: I b 1 1',hash:'7BB52B3A'XL}])
q2=ptr_new([ $
    {name:'',ptr:r0},{name:'6: P 1 m 1',ptr:r1},{name:'7: P 1 c 1',ptr:r2},{name:'8: C 1 m 1',ptr:r3}, $
    {name:'9: C 1 c 1',ptr:r4}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'57: P 1 2/m 1',hash:'068DFEB3'XL},{name:'58: P 1 1 2/m',hash:'9B92779F'XL},{name:'59: P 2/m 1 1',hash:'B08D46CF'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'60: P 1 21/m 1',hash:'54FA8558'XL},{name:'61: P 1 1 21/m',hash:'31194672'XL},{name:'62: P 21/m 1 1',hash:'CE8251DF'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'63: C 1 2/m 1',hash:'09EFDD05'XL},{name:'64: A 1 2/m 1',hash:'8E118804'XL},{name:'65: I 1 2/m 1',hash:'5C5D4C9F'XL}, $
    {name:'66: A 1 1 2/m',hash:'5A48CC10'XL},{name:'67: B 1 1 2/m',hash:'97E84B5C'XL},{name:'68: I 1 1 2/m',hash:'8804088B'XL},{name:'69: B 2/m 1 1',hash:'05A3ECAD'XL}, $
    {name:'70: C 2/m 1 1',hash:'4FFD3EE0'XL},{name:'71: I 2/m 1 1',hash:'1A4FAF7A'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'72: P 1 2/c 1',hash:'AC06CF5E'XL},{name:'73: P 1 2/n 1',hash:'F817D0BF'XL},{name:'74: P 1 2/a 1',hash:'529CE152'XL}, $
    {name:'75: P 1 1 2/a',hash:'CF83687E'XL},{name:'76: P 1 1 2/n',hash:'9DF41395'XL},{name:'77: P 1 1 2/b',hash:'C9E50C74'XL},{name:'78: P 2/b 1 1',hash:'4AADA62E'XL}, $
    {name:'79: P 2/n 1 1',hash:'F45A1AA3'XL},{name:'80: P 2/c 1 1',hash:'E58CA26F'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'81: P 1 21/c 1',hash:'41572403'XL},{name:'82: P 1 21/n 1',hash:'15463BE2'XL},{name:'83: P 1 21/a 1',hash:'00EB9AB9'XL}, $
    {name:'84: P 1 1 21/a',hash:'65085993'XL},{name:'85: P 1 1 21/n',hash:'8859B2CE'XL},{name:'86: P 1 1 21/b',hash:'DC48AD2F'XL},{name:'87: P 21/b 1 1',hash:'34A2B13E'XL}, $
    {name:'88: P 21/n 1 1',hash:'8A550DB3'XL},{name:'89: P 21/c 1 1',hash:'9B83B57F'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'90: C 1 2/c 1',hash:'673D531E'XL},{name:'91: A 1 2/n 1',hash:'428E62A8'XL},{name:'92: I 1 2/a 1',hash:'B7F841EA'XL}, $
    {name:'93: A 1 2/a 1',hash:'65B48571'XL},{name:'94: C 1 2/n 1',hash:'8C985E6B'XL},{name:'95: I 1 2/c 1',hash:'90C2A633'XL},{name:'96: A 1 1 2/a',hash:'B1EDC165'XL}, $
    {name:'97: B 1 1 2/n',hash:'5B77A1F0'XL},{name:'98: I 1 1 2/b',hash:'AF3EEF52'XL},{name:'99: B 1 1 2/b',hash:'B0D2AC85'XL},{name:'100: A 1 1 2/n',hash:'96D726BC'XL}, $
    {name:'101: I 1 1 2/a',hash:'63A105FE'XL},{name:'102: B 2/b 1 1',hash:'E4301506'XL},{name:'103: C 2/n 1 1',hash:'AECDC717'XL},{name:'104: I 2/c 1 1',hash:'62F43297'XL}, $
    {name:'105: C 2/c 1 1',hash:'37E5A351'XL},{name:'106: B 2/n 1 1',hash:'7D187140'XL},{name:'107: I 2/b 1 1',hash:'FBDC56D1'XL}])
q3=ptr_new([ $
    {name:'',ptr:r0},{name:'10: P 1 2/m 1',ptr:r1},{name:'11: P 1 21/m 1',ptr:r2},{name:'12: C 1 2/m 1',ptr:r3}, $
    {name:'13: P 1 2/c 1',ptr:r4},{name:'14: P 1 21/c 1',ptr:r5},{name:'15: C 1 2/c 1',ptr:r6}])
p2=ptr_new([ $
    {name:'',ptr:q0},{name:'2',ptr:q1},{name:'m',ptr:q2},{name:'2/m',ptr:q3}])
q0=ptr_new()
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'108: P 2 2 2',hash:'A959FC0B'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'109: P 2 2 21',hash:'03D2CDE6'XL},{name:'110: P 21 2 2',hash:'D756EB1B'XL},{name:'111: P 2 21 2',hash:'010E6701'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'112: P 21 21 2',hash:'2B106FF0'XL},{name:'113: P 2 21 21',hash:'EE8326BB'XL},{name:'114: P 21 2 21',hash:'82570FBB'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'115: P 21 21 21',hash:'8F7A6EEC'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'116: C 2 2 21',hash:'D3DCEAE0'XL},{name:'117: A 21 2 2',hash:'A3D855BC'XL},{name:'118: B 2 21 2',hash:'31F9A8C4'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'119: C 2 2 2',hash:'BD0E64FB'XL},{name:'120: A 2 2 2',hash:'3AF031FA'XL},{name:'121: B 2 2 2',hash:'F750B6B6'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'122: F 2 2 2',hash:'8A25CD68'XL}])
r8=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'123: I 2 2 2',hash:'E8BCF561'XL}])
r9=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'124: I 21 21 21',hash:'7BA265F9'XL}])
q1=ptr_new([ $
    {name:'',ptr:r0},{name:'16: P 2 2 2',ptr:r1},{name:'17: P 2 2 21',ptr:r2},{name:'18: P 21 21 2',ptr:r3}, $
    {name:'19: P 21 21 21',ptr:r4},{name:'20: C 2 2 21',ptr:r5},{name:'21: C 2 2 2',ptr:r6},{name:'22: F 2 2 2',ptr:r7}, $
    {name:'23: I 2 2 2',ptr:r8},{name:'24: I 21 21 21',ptr:r9}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'125: P m m 2',hash:'715A8A38'XL},{name:'126: P 2 m m',hash:'34467527'XL},{name:'127: P m 2 m',hash:'454292CE'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'128: P m c 21',hash:'DBD1BBD5'XL},{name:'129: P c m 21',hash:'8ED05F75'XL},{name:'130: P 21 m a',hash:'1E587DD6'XL}, $
    {name:'131: P 21 a m',hash:'4A496237'XL},{name:'132: P b 21 m',hash:'ED1509C4'XL},{name:'133: P m 21 b',hash:'1735E925'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'134: P c c 2',hash:'245B6E98'XL},{name:'135: P 2 a a',hash:'60576AC6'XL},{name:'136: P b 2 b',hash:'BF62722F'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'137: P m a 2',hash:'0F559D28'XL},{name:'138: P b m 2',hash:'8B7A6AD9'XL},{name:'139: P 2 m b',hash:'66310ECC'XL}, $
    {name:'140: P 2 c m',hash:'9ECD44CA'XL},{name:'141: P c 2 m',hash:'1043766E'XL},{name:'142: P m 2 a',hash:'3B4D85DE'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'143: P c a 21',hash:'F0DF4865'XL},{name:'144: P b c 21',hash:'C427E492'XL},{name:'145: P 21 a b',hash:'183E19DC'XL}, $
    {name:'146: P 21 c a',hash:'B4D34C3B'XL},{name:'147: P c 21 b',hash:'A7E2B223'XL},{name:'148: P b 21 a',hash:'931A1ED4'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'149: P n c 2',hash:'358DD654'XL},{name:'150: P c n 2',hash:'5A547988'XL},{name:'151: P 2 n a',hash:'CADC5B2B'XL}, $
    {name:'152: P 2 a n',hash:'3220112D'XL},{name:'153: P b 2 n',hash:'C16D653F'XL},{name:'154: P n 2 b',hash:'0195CEA2'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'155: P m n 21',hash:'8FC0A434'XL},{name:'156: P n m 21',hash:'72570CE4'XL},{name:'157: P 21 m n',hash:'B60FE6DC'XL}, $
    {name:'158: P 21 n m',hash:'1F488697'XL},{name:'159: P n 21 m',hash:'464F1412'XL},{name:'160: P m 21 n',hash:'4324F6C4'XL}])
r8=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'161: P b a 2',hash:'F5757DC9'XL},{name:'162: P 2 c b',hash:'739CAF97'XL},{name:'163: P c 2 a',hash:'6E4C617E'XL}])
r9=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'164: P n a 21',hash:'04DF4F0F'XL},{name:'165: P b n 21',hash:'BA28F382'XL},{name:'166: P 21 n b',hash:'0D93B887'XL}, $
    {name:'167: P 21 c n',hash:'5982A766'XL},{name:'168: P c 21 n',hash:'D9EDA533'XL},{name:'169: P n 21 a',hash:'C83B1DFF'XL}])
r10=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'170: P n n 2',hash:'4B82C144'XL},{name:'171: P 2 n n',hash:'278DB076'XL},{name:'172: P n 2 n',hash:'7F9AD9B2'XL}])
r11=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'173: C m m 2',hash:'93D6F19F'XL},{name:'174: A 2 m m',hash:'EEA975EE'XL},{name:'175: B m 2 m',hash:'B633CB54'XL}])
r12=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'176: C m c 21',hash:'FD047F84'XL},{name:'177: C c m 21',hash:'851CE235'XL},{name:'178: A 21 m a',hash:'9C241CDD'XL}, $
    {name:'179: A 21 a m',hash:'778111A8'XL},{name:'180: B b 21 m',hash:'709AD526'XL},{name:'181: B m 21 b',hash:'91092C8D'XL}])
r13=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'182: C c c 2',hash:'EBCE6C2E'XL},{name:'183: A 2 a a',hash:'050C789B'XL},{name:'184: B b 2 b',hash:'57A032FF'XL}])
r14=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'185: A m m 2',hash:'1428A49E'XL},{name:'186: B m m 2',hash:'D98823D2'XL},{name:'187: B 2 m m',hash:'2309F2A2'XL}, $
    {name:'188: C 2 m m',hash:'695720EF'XL},{name:'189: C m 2 m',hash:'FC6D1919'XL},{name:'190: A m 2 m',hash:'7B934C18'XL}])
r15=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'191: A b m 2',hash:'F5BB5D35'XL},{name:'192: B m a 2',hash:'40A04794'XL},{name:'193: B 2 c m',hash:'C8ACFFD7'XL}, $
    {name:'194: C 2 m b',hash:'82F22D9A'XL},{name:'195: C m 2 a',hash:'65457D5F'XL},{name:'196: A c 2 m',hash:'9A00B5B3'XL}])
r16=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'197: A m a 2',hash:'8D00C0D8'XL},{name:'198: B b m 2',hash:'381BDA79'XL},{name:'199: B 2 m b',hash:'0433157B'XL}, $
    {name:'200: C 2 c m',hash:'0785AEF4'XL},{name:'201: C c 2 m',hash:'847584A8'XL},{name:'202: A m 2 a',hash:'E2BB285E'XL}])
r17=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'203: A b a 2',hash:'6C933973'XL},{name:'204: B b a 2',hash:'A133BE3F'XL},{name:'205: B 2 c b',hash:'EF96180E'XL}, $
    {name:'206: C 2 c b',hash:'EC20A381'XL},{name:'207: C c 2 a',hash:'1D5DE0EE'XL},{name:'208: A c 2 a',hash:'0328D1F5'XL}])
r18=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'209: F m m 2',hash:'18172B53'XL},{name:'210: F 2 m m',hash:'AD1A17DB'XL},{name:'211: F m 2 m',hash:'A88FC3E9'XL}])
r19=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'212: F d d 2',hash:'921865A5'XL},{name:'213: F 2 d d',hash:'1B59DBD3'XL},{name:'214: F d 2 d',hash:'22808D1F'XL}])
r20=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'215: I m m 2',hash:'C6646005'XL},{name:'216: I 2 m m',hash:'3CE5B175'XL},{name:'217: I m 2 m',hash:'A9DF8883'XL}])
r21=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'218: I b a 2',hash:'BEDFFDE8'XL},{name:'219: I 2 c b',hash:'D740BC00'XL},{name:'220: I c 2 a',hash:'484C7128'XL}])
r22=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'221: I m a 2',hash:'5F4C0443'XL},{name:'222: I b m 2',hash:'27F799AE'XL},{name:'223: I 2 m b',hash:'1BDF56AC'XL}, $
    {name:'224: I 2 c m',hash:'F07A5BD9'XL},{name:'225: I c 2 m',hash:'D164156E'XL},{name:'226: I m 2 a',hash:'30F7ECC5'XL}])
q2=ptr_new([ $
    {name:'',ptr:r0},{name:'25: P m m 2',ptr:r1},{name:'26: P m c 21',ptr:r2},{name:'27: P c c 2',ptr:r3}, $
    {name:'28: P m a 2',ptr:r4},{name:'29: P c a 21',ptr:r5},{name:'30: P n c 2',ptr:r6},{name:'31: P m n 21',ptr:r7}, $
    {name:'32: P b a 2',ptr:r8},{name:'33: P n a 21',ptr:r9},{name:'34: P n n 2',ptr:r10},{name:'35: C m m 2',ptr:r11}, $
    {name:'36: C m c 21',ptr:r12},{name:'37: C c c 2',ptr:r13},{name:'38: A m m 2',ptr:r14},{name:'39: A b m 2',ptr:r15}, $
    {name:'40: A m a 2',ptr:r16},{name:'41: A b a 2',ptr:r17},{name:'42: F m m 2',ptr:r18},{name:'43: F d d 2',ptr:r19}, $
    {name:'44: I m m 2',ptr:r20},{name:'45: I b a 2',ptr:r21},{name:'46: I m a 2',ptr:r22}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'227: P m m m',hash:'E7243BBC'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'228: P n n n:1',hash:'51E85B6C'XL},{name:'229: P n n n:2',hash:'78FF5B81'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'230: P c c m',hash:'353736E5'XL},{name:'231: P m a a',hash:'0C8136C9'XL},{name:'232: P b m b',hash:'02B17898'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'233: P b a n:1',hash:'1392524E'XL},{name:'234: P b a n:2',hash:'E91475ED'XL},{name:'235: P n c b:1',hash:'4E9DE450'XL}, $
    {name:'236: P n c b:2',hash:'935A56F4'XL},{name:'237: P c n a:1',hash:'38BBAE1F'XL},{name:'238: P c n a:2',hash:'DE923B90'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'239: P m m a',hash:'6179E0C6'XL},{name:'240: P m m b',hash:'2BE88448'XL},{name:'241: P b m m',hash:'CE7DC76C'XL}, $
    {name:'242: P c m m',hash:'B013E0D3'XL},{name:'243: P m c m',hash:'6200ED8A'XL},{name:'244: P m a m',hash:'8ADCEDB3'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'245: P n n a',hash:'C514E0E2'XL},{name:'246: P n n b',hash:'FEA280FB'XL},{name:'247: P b n n',hash:'38E06B40'XL}, $
    {name:'248: P c n n',hash:'637980F3'XL},{name:'249: P n c n',hash:'15078D8E'XL},{name:'250: P n a n',hash:'A90B452C'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'251: P m n a',hash:'89A5E0FF'XL},{name:'252: P n m b',hash:'42AE4859'XL},{name:'253: P b m n',hash:'84ECA3E2'XL}, $
    {name:'254: P c n m',hash:'58CFE0EA'XL},{name:'255: P n c m',hash:'2EB1ED97'XL},{name:'256: P m a n',hash:'C04D893D'XL}])
r8=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'257: P c c a',hash:'B36AED9F'XL},{name:'258: P c c b',hash:'88DC8D86'XL},{name:'259: P b a a',hash:'25D8CA19'XL}, $
    {name:'260: P c a a',hash:'5BB6EDA6'XL},{name:'261: P b c b',hash:'D3456635'XL},{name:'262: P b a b',hash:'6F49AE97'XL}])
r9=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'263: P b a m',hash:'A3851163'XL},{name:'264: P m c b',hash:'8B3B9E72'XL},{name:'265: P c m a',hash:'364E3BA9'XL}])
r10=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'266: P c c n',hash:'0E8156FC'XL},{name:'267: P n a a',hash:'31173243'XL},{name:'268: P b n b',hash:'BEBDB03A'XL}])
r11=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'269: P b c m',hash:'3A7E15CD'XL},{name:'270: P c a m',hash:'DDEB36DC'XL},{name:'271: P m c a',hash:'E45D36F0'XL}, $
    {name:'272: P m a b',hash:'46105247'XL},{name:'273: P b m a',hash:'48201C16'XL},{name:'274: P c m b',hash:'280F97BC'XL}])
r12=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'275: P n n m',hash:'43493B98'XL},{name:'276: P m n n',hash:'609E9307'XL},{name:'277: P n m n',hash:'C4F39323'XL}])
r13=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'278: P m m n:1',hash:'57337891'XL},{name:'279: P m m n:2',hash:'ADB55F32'XL},{name:'280: P n m m:1',hash:'2282419E'XL}, $
    {name:'281: P n m m:2',hash:'DAB23F36'XL},{name:'282: P m n m:1',hash:'E9D1AE0A'XL},{name:'283: P m n m:2',hash:'0FF83B85'XL}])
r14=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'284: P b c n',hash:'5518BD4F'XL},{name:'285: P c a n',hash:'C3AA9AC9'XL},{name:'286: P n c a',hash:'A8EC36ED'XL}, $
    {name:'287: P n a b',hash:'2F569E56'XL},{name:'288: P b n a',hash:'D1DB18B8'XL},{name:'289: P c n b',hash:'E5245B89'XL}])
r15=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'290: P b c a',hash:'BC23CEB7'XL},{name:'291: P c a b',hash:'45F741B3'XL}])
r16=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'292: P n m a',hash:'5CEFE44C'XL},{name:'293: P m n b',hash:'E6C3487D'XL},{name:'294: P b n m',hash:'5786C3C2'XL}, $
    {name:'295: P c m n',hash:'AE524CC6'XL},{name:'296: P m c n',hash:'0D664508'XL},{name:'297: P n a m',hash:'B74AE939'XL}])
r17=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'298: C m c m',hash:'C0DD6737'XL},{name:'299: C c m m',hash:'A4E78801'XL},{name:'300: A m m a',hash:'637AC153'XL}, $
    {name:'301: A m a m',hash:'7AF93F9C'XL},{name:'302: B b m m',hash:'5EA5920D'XL},{name:'303: B m m b',hash:'DF056B84'XL}])
r18=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'304: C m c a',hash:'D95E99F8'XL},{name:'305: C c m b',hash:'014F7F4E'XL},{name:'306: A b m a',hash:'E2DA38DA'XL}, $
    {name:'307: A c a m',hash:'3E8A450B'XL},{name:'308: B b c m',hash:'E28E9B8D'XL},{name:'309: B m a b',hash:'C686954B'XL}])
r19=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'310: C m m m',hash:'CC263714'XL},{name:'311: A m m m',hash:'C6D2361C'XL},{name:'312: B m m m',hash:'1AD6E89A'XL}])
r20=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'313: C c c m',hash:'A81CD822'XL},{name:'314: A m a a',hash:'DF51C8D3'XL},{name:'315: B b m b',hash:'9B761113'XL}])
r21=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'316: C m m a',hash:'D5A5C9DB'XL},{name:'317: C m m b',hash:'698EC05B'XL},{name:'318: A b m m',hash:'4772CF95'XL}, $
    {name:'319: A c m m',hash:'82A14C8B'XL},{name:'320: B m c m',hash:'A6FDE11A'XL},{name:'321: B m a m',hash:'03551655'XL}])
r22=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'322: C c c a:1',hash:'D12AB103'XL},{name:'323: C c c a:2',hash:'0DB42F6D'XL},{name:'322: C c c a:1',hash:'D12AB103'XL}, $
    {name:'325: C c c b:2',hash:'B19F26ED'XL},{name:'326: A b a a:1',hash:'6A64CF1D'XL},{name:'327: A b a a:2',hash:'9B22B244'XL},{name:'326: A b a a:1',hash:'6A64CF1D'XL}, $
    {name:'329: A c a a:2',hash:'5EF1315A'XL},{name:'330: B b c b:1',hash:'B660119B'XL},{name:'331: B b c b:2',hash:'82F5EFDC'XL},{name:'330: B b c b:1',hash:'B660119B'XL}, $
    {name:'333: B b a b:2',hash:'275D1893'XL}])
r23=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'334: F m m m',hash:'CD91580E'XL}])
r24=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'335: F d d d:1',hash:'D31C4512'XL},{name:'336: F d d d:2',hash:'C3FDC94F'XL}])
r25=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'337: I m m m',hash:'4A8CA5FC'XL}])
r26=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'338: I b a m',hash:'B2D4D6EB'XL},{name:'339: I m c b',hash:'530F5B33'XL},{name:'340: I c m a',hash:'CB2C5C75'XL}])
r27=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'341: I b c a',hash:'770755F5'XL},{name:'342: I c a b',hash:'177C21A4'XL}])
r28=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'343: I m m a',hash:'8F5F26E2'XL},{name:'344: I m m b',hash:'EF2452B3'XL},{name:'345: I b m m',hash:'6E84AB3A'XL}, $
    {name:'346: I c m m',hash:'0EFFDF6B'XL},{name:'347: I m c m',hash:'F6A7AC7C'XL},{name:'348: I m a m',hash:'96DCD82D'XL}])
q3=ptr_new([ $
    {name:'',ptr:r0},{name:'47: P m m m',ptr:r1},{name:'48: P n n n:1',ptr:r2},{name:'49: P c c m',ptr:r3}, $
    {name:'50: P b a n:1',ptr:r4},{name:'51: P m m a',ptr:r5},{name:'52: P n n a',ptr:r6},{name:'53: P m n a',ptr:r7}, $
    {name:'54: P c c a',ptr:r8},{name:'55: P b a m',ptr:r9},{name:'56: P c c n',ptr:r10},{name:'57: P b c m',ptr:r11}, $
    {name:'58: P n n m',ptr:r12},{name:'59: P m m n:1',ptr:r13},{name:'60: P b c n',ptr:r14},{name:'61: P b c a',ptr:r15}, $
    {name:'62: P n m a',ptr:r16},{name:'63: C m c m',ptr:r17},{name:'64: C m c a',ptr:r18},{name:'65: C m m m',ptr:r19}, $
    {name:'66: C c c m',ptr:r20},{name:'67: C m m a',ptr:r21},{name:'68: C c c a:1',ptr:r22},{name:'69: F m m m',ptr:r23}, $
    {name:'70: F d d d:1',ptr:r24},{name:'71: I m m m',ptr:r25},{name:'72: I b a m',ptr:r26},{name:'73: I b c a',ptr:r27}, $
    {name:'74: I m m a',ptr:r28}])
p3=ptr_new([ $
    {name:'',ptr:q0},{name:'222',ptr:q1},{name:'mm2',ptr:q2},{name:'mmm',ptr:q3}])
q0=ptr_new()
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'349: P 4',hash:'E194247E'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'350: P 41',hash:'EB09F6D0'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'351: P 42',hash:'1E1EF133'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'352: P 43',hash:'21968A94'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'353: I 4',hash:'8D2DB265'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'354: I 41',hash:'04C2F7B9'XL}])
q1=ptr_new([ $
    {name:'',ptr:r0},{name:'75: P 4',ptr:r1},{name:'76: P 41',ptr:r2},{name:'77: P 42',ptr:r3}, $
    {name:'78: P 43',ptr:r4},{name:'79: I 4',ptr:r5},{name:'80: I 41',ptr:r6}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'355: P -4',hash:'5B1380EC'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'356: I -4',hash:'429018FC'XL}])
q2=ptr_new([ $
    {name:'',ptr:r0},{name:'81: P -4',ptr:r1},{name:'82: I -4',ptr:r2}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'357: P 4/m',hash:'7F473453'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'358: P 42/m',hash:'07DF08E7'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'359: P 4/n:1',hash:'8D9739AB'XL},{name:'360: P 4/n:2',hash:'28870B39'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'361: P 42/n:1',hash:'DA4091D2'XL},{name:'362: P 42/n:2',hash:'CFC232BB'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'363: I 4/m',hash:'E426791E'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'364: I 41/a:1',hash:'DD05FE66'XL},{name:'365: I 41/a:2',hash:'F5A09FA3'XL}])
q3=ptr_new([ $
    {name:'',ptr:r0},{name:'83: P 4/m',ptr:r1},{name:'84: P 42/m',ptr:r2},{name:'85: P 4/n:1',ptr:r3}, $
    {name:'86: P 42/n:1',ptr:r4},{name:'87: I 4/m',ptr:r5},{name:'88: I 41/a:1',ptr:r6}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'366: P 4 2 2',hash:'781566A5'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'367: P 42 1 2',hash:'529DBAFE'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'368: P 41 2 2',hash:'1129373C'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'369: P 41 21 2',hash:'7ABBC4B0'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'370: P 42 2 2',hash:'00A13651'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'371: P 42 21 2',hash:'B11A4EB4'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'372: P 43 2 2',hash:'C507CD54'XL}])
r8=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'373: P 43 21 2',hash:'5AF2C5D9'XL}])
r9=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'374: I 4 2 2',hash:'CC0E6DAD'XL}])
r10=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'375: I 41 2 2',hash:'AC45C901'XL}])
q4=ptr_new([ $
    {name:'',ptr:r0},{name:'89: P 4 2 2',ptr:r1},{name:'90: P 42 1 2',ptr:r2},{name:'91: P 41 2 2',ptr:r3}, $
    {name:'92: P 41 21 2',ptr:r4},{name:'93: P 42 2 2',ptr:r5},{name:'94: P 42 21 2',ptr:r6},{name:'95: P 43 2 2',ptr:r7}, $
    {name:'96: P 43 21 2',ptr:r8},{name:'97: I 4 2 2',ptr:r9},{name:'98: I 41 2 2',ptr:r10}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'376: P 4 m m',hash:'93373957'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'377: P 4 b m',hash:'D7961388'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'378: P 42 c m',hash:'53DD13C8'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'379: P 42 n m',hash:'25FB5987'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'380: P 4 c c',hash:'4124340E'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'381: P 4 n c',hash:'375A3973'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'382: P 42 m c',hash:'81CE1E91'XL}])
r8=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'383: P 42 b c',hash:'5B8D5C98'XL}])
r9=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'384: I 4 m m',hash:'D27B6F7F'XL}])
r10=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'385: I 4 c m',hash:'2A231C68'XL}])
r11=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'386: I 41 m d',hash:'6EE4A2B9'XL}])
r12=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'387: I 41 c d',hash:'F0719502'XL}])
q5=ptr_new([ $
    {name:'',ptr:r0},{name:'99: P 4 m m',ptr:r1},{name:'100: P 4 b m',ptr:r2},{name:'101: P 42 c m',ptr:r3}, $
    {name:'102: P 42 n m',ptr:r4},{name:'103: P 4 c c',ptr:r5},{name:'104: P 4 n c',ptr:r6},{name:'105: P 42 m c',ptr:r7}, $
    {name:'106: P 42 b c',ptr:r8},{name:'107: I 4 m m',ptr:r9},{name:'108: I 4 c m',ptr:r10},{name:'109: I 41 m d',ptr:r11}, $
    {name:'110: I 41 c d',ptr:r12}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'388: P -4 2 m',hash:'AE367C74'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'389: P -4 2 c',hash:'C47B0B46'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'390: P -4 21 m',hash:'84BEA02F'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'391: P -4 21 c',hash:'67395465'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'392: P -4 m 2',hash:'D50D13D5'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'393: P -4 c 2',hash:'3A48DE91'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'394: P -4 b 2',hash:'00B4ACAC'XL}])
r8=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'395: P -4 n 2',hash:'F7B01EEF'XL}])
r9=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'396: I -4 m 2',hash:'7A127A20'XL}])
r10=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'397: I -4 c 2',hash:'C66B145C'XL}])
r11=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'398: I -4 2 m',hash:'93E9BB9A'XL}])
r12=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'399: I -4 2 d',hash:'F3A21F36'XL}])
q6=ptr_new([ $
    {name:'',ptr:r0},{name:'111: P -4 2 m',ptr:r1},{name:'112: P -4 2 c',ptr:r2},{name:'113: P -4 21 m',ptr:r3}, $
    {name:'114: P -4 21 c',ptr:r4},{name:'115: P -4 m 2',ptr:r5},{name:'116: P -4 c 2',ptr:r6},{name:'117: P -4 b 2',ptr:r7}, $
    {name:'118: P -4 n 2',ptr:r8},{name:'119: I -4 m 2',ptr:r9},{name:'120: I -4 c 2',ptr:r10},{name:'121: I -4 2 m',ptr:r11}, $
    {name:'122: I -4 2 d',ptr:r12}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'400: P 4/m m m',hash:'B8F2113E'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'401: P 4/m c c',hash:'2DE869D5'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'402: P 4/n b m:1',hash:'B7FA689E'XL},{name:'403: P 4/n b m:2',hash:'026E1F11'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'404: P 4/n n c:1',hash:'08A4EB7F'XL},{name:'405: P 4/n n c:2',hash:'51368CA3'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'406: P 4/m b m',hash:'D5A84150'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'407: P 4/m n c',hash:'F1B4F061'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'408: P 4/n m m:1',hash:'81D6F60A'XL},{name:'409: P 4/n m m:2',hash:'6F344F7F'XL}])
r8=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'410: P 4/n c c:1',hash:'9AEE9967'XL},{name:'411: P 4/n c c:2',hash:'8D6A1517'XL}])
r9=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'412: P 42/m m c',hash:'5A9858B3'XL}])
r10=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'413: P 42/m c m',hash:'CF822058'XL}])
r11=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'414: P 42/n b c:1',hash:'00E9903A'XL},{name:'415: P 42/n b c:2',hash:'F9DF795B'XL}])
r12=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'416: P 42/n n m:1',hash:'AAD7368D'XL},{name:'417: P 42/n n m:2',hash:'AA87EAE9'XL}])
r13=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'418: P 42/m b c',hash:'A5C2E426'XL}])
r14=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'419: P 42/m n m',hash:'53C72D93'XL}])
r15=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'420: P 42/n m c:1',hash:'E391D7D2'XL},{name:'421: P 42/n m c:2',hash:'0685C5CE'XL}])
r16=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'422: P 42/n c m:1',hash:'F8A9B8BF'XL},{name:'423: P 42/n c m:2',hash:'E4DB9FA6'XL}])
r17=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'424: I 4/m m m',hash:'7A639EEF'XL}])
r18=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'425: I 4/m c m',hash:'D7E16F02'XL}])
r19=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'426: I 41/a m d:1',hash:'E12277EA'XL},{name:'427: I 41/a m d:2',hash:'B96409E0'XL}])
r20=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'428: I 41/a c d:1',hash:'8BDD8359'XL},{name:'429: I 41/a c d:2',hash:'08D64E97'XL}])
q7=ptr_new([ $
    {name:'',ptr:r0},{name:'123: P 4/m m m',ptr:r1},{name:'124: P 4/m c c',ptr:r2},{name:'125: P 4/n b m:1',ptr:r3}, $
    {name:'126: P 4/n n c:1',ptr:r4},{name:'127: P 4/m b m',ptr:r5},{name:'128: P 4/m n c',ptr:r6},{name:'129: P 4/n m m:1',ptr:r7}, $
    {name:'130: P 4/n c c:1',ptr:r8},{name:'131: P 42/m m c',ptr:r9},{name:'132: P 42/m c m',ptr:r10},{name:'133: P 42/n b c:1',ptr:r11}, $
    {name:'134: P 42/n n m:1',ptr:r12},{name:'135: P 42/m b c',ptr:r13},{name:'136: P 42/m n m',ptr:r14},{name:'137: P 42/n m c:1',ptr:r15}, $
    {name:'138: P 42/n c m:1',ptr:r16},{name:'139: I 4/m m m',ptr:r17},{name:'140: I 4/m c m',ptr:r18},{name:'141: I 41/a m d:1',ptr:r19}, $
    {name:'142: I 41/a c d:1',ptr:r20}])
p4=ptr_new([ $
    {name:'',ptr:q0},{name:'4',ptr:q1},{name:'-4',ptr:q2},{name:'4/m',ptr:q3}, $
    {name:'422',ptr:q4},{name:'4mm',ptr:q5},{name:'-42m',ptr:q6},{name:'4/mmm',ptr:q7}])
q0=ptr_new()
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'430: P 3',hash:'329C980E'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'431: P 31',hash:'86DD212B'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'432: P 32',hash:'D70A5F46'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'433: R 3:H',hash:'8E4E25F8'XL},{name:'434: R 3:R',hash:'D5A0AA2D'XL}])
q1=ptr_new([ $
    {name:'',ptr:r0},{name:'143: P 3',ptr:r1},{name:'144: P 31',ptr:r2},{name:'145: P 32',ptr:r3}, $
    {name:'146: R 3:H',ptr:r4}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'435: P -3',hash:'FDD759B5'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'436: R -3:H',hash:'BE8D0D7F'XL},{name:'437: R -3:R',hash:'D9A29BAC'XL}])
q2=ptr_new([ $
    {name:'',ptr:r0},{name:'147: P -3',ptr:r1},{name:'148: R -3:H',ptr:r2}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'438: P 3 1 2',hash:'65B7A72B'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'439: P 3 2 1',hash:'C1840A7A'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'440: P 31 1 2',hash:'97E2DFD5'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'441: P 31 2 1',hash:'33D17284'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'442: P 32 1 2',hash:'E39F36B4'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'443: P 32 2 1',hash:'47AC9BE5'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'444: R 32:H',hash:'46EBEE09'XL},{name:'445: R 32:R',hash:'A20B8591'XL}])
q3=ptr_new([ $
    {name:'',ptr:r0},{name:'149: P 3 1 2',ptr:r1},{name:'150: P 3 2 1',ptr:r2},{name:'151: P 31 1 2',ptr:r3}, $
    {name:'152: P 31 2 1',ptr:r4},{name:'153: P 32 1 2',ptr:r5},{name:'154: P 32 2 1',ptr:r6},{name:'155: R 32:H',ptr:r7}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'446: P 3 m 1',hash:'9F4CFFAA'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'447: P 3 1 m',hash:'39859B12'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'448: P 3 c 1',hash:'E04FE588'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'449: P 3 1 c',hash:'EC0DB0DD'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'450: R 3 m:H',hash:'6E32558C'XL},{name:'451: R 3 m:R',hash:'B951B4F7'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'452: R 3 c:H',hash:'F7D1A830'XL},{name:'453: R 3 c:R',hash:'219BE015'XL}])
q4=ptr_new([ $
    {name:'',ptr:r0},{name:'156: P 3 m 1',ptr:r1},{name:'157: P 3 1 m',ptr:r2},{name:'158: P 3 c 1',ptr:r3}, $
    {name:'159: P 3 1 c',ptr:r4},{name:'160: R 3 m:H',ptr:r5},{name:'161: R 3 c:H',ptr:r6}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'454: P -3 1 m',hash:'F74C7F83'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'455: P -3 1 c',hash:'69DC2F41'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'456: P -3 m 1',hash:'FC3EDAFB'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'457: P -3 c 1',hash:'E9D82A99'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'458: R -3 m:H',hash:'6DF507A9'XL},{name:'459: R -3 m:R',hash:'1C80E47A'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'460: R -3 c:H',hash:'9A7D09D3'XL},{name:'461: R -3 c:R',hash:'BB691C91'XL}])
q5=ptr_new([ $
    {name:'',ptr:r0},{name:'162: P -3 1 m',ptr:r1},{name:'163: P -3 1 c',ptr:r2},{name:'164: P -3 m 1',ptr:r3}, $
    {name:'165: P -3 c 1',ptr:r4},{name:'166: R -3 m:H',ptr:r5},{name:'167: R -3 c:H',ptr:r6}])
p5=ptr_new([ $
    {name:'',ptr:q0},{name:'3',ptr:q1},{name:'-3',ptr:q2},{name:'32',ptr:q3}, $
    {name:'3m',ptr:q4},{name:'-3m',ptr:q5}])
q0=ptr_new()
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'462: P 6',hash:'A2DDAF47'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'463: P 61',hash:'81A2968A'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'464: P 65',hash:'84B83B8B'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'465: P 62',hash:'62D97E5D'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'466: P 64',hash:'B22FEFA7'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'467: P 63',hash:'DDDEB565'XL}])
q1=ptr_new([ $
    {name:'',ptr:r0},{name:'168: P 6',ptr:r1},{name:'169: P 61',ptr:r2},{name:'170: P 65',ptr:r3}, $
    {name:'171: P 62',ptr:r4},{name:'172: P 64',ptr:r5},{name:'173: P 63',ptr:r6}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'468: P -6',hash:'35EB6CCB'XL}])
q2=ptr_new([ $
    {name:'',ptr:r0},{name:'174: P -6',ptr:r1}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'469: P 6/m',hash:'32DACFB6'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'470: P 63/m',hash:'06C1AE99'XL}])
q3=ptr_new([ $
    {name:'',ptr:r0},{name:'175: P 6/m',ptr:r1},{name:'176: P 63/m',ptr:r2}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'471: P 6 2 2',hash:'90230E5F'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'472: P 61 2 2',hash:'EED642A2'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'473: P 65 2 2',hash:'891AEB21'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'474: P 62 2 2',hash:'6F55448F'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'475: P 64 2 2',hash:'655E4CD9'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'476: P 63 2 2',hash:'A4386F70'XL}])
q4=ptr_new([ $
    {name:'',ptr:r0},{name:'177: P 6 2 2',ptr:r1},{name:'178: P 61 2 2',ptr:r2},{name:'179: P 65 2 2',ptr:r3}, $
    {name:'180: P 62 2 2',ptr:r4},{name:'181: P 64 2 2',ptr:r5},{name:'182: P 63 2 2',ptr:r6}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'477: P 6 m m',hash:'3B0A2D17'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'478: P 6 c c',hash:'B74F3653'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'479: P 63 c m',hash:'F4D8E94D'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'480: P 63 m c',hash:'789DF209'XL}])
q5=ptr_new([ $
    {name:'',ptr:r0},{name:'183: P 6 m m',ptr:r1},{name:'184: P 6 c c',ptr:r2},{name:'185: P 63 c m',ptr:r3}, $
    {name:'186: P 63 m c',ptr:r4}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'481: P -6 m 2',hash:'EEACB736'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'482: P -6 c 2',hash:'FB4A4754'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'483: P -6 2 m',hash:'55B5BE6A'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'484: P -6 2 c',hash:'CB25EEA8'XL}])
q6=ptr_new([ $
    {name:'',ptr:r0},{name:'187: P -6 m 2',ptr:r1},{name:'188: P -6 c 2',ptr:r2},{name:'189: P -6 2 m',ptr:r3}, $
    {name:'190: P -6 2 c',ptr:r4}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'485: P 6/m m m',hash:'F1FC7952'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'486: P 6/m c c',hash:'87F9E8CA'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'487: P 63/m c m',hash:'1EB41BD9'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'488: P 63/m m c',hash:'68B18A41'XL}])
q7=ptr_new([ $
    {name:'',ptr:r0},{name:'191: P 6/m m m',ptr:r1},{name:'192: P 6/m c c',ptr:r2},{name:'193: P 63/m c m',ptr:r3}, $
    {name:'194: P 63/m m c',ptr:r4}])
p6=ptr_new([ $
    {name:'',ptr:q0},{name:'6',ptr:q1},{name:'-6',ptr:q2},{name:'6/m',ptr:q3}, $
    {name:'622',ptr:q4},{name:'6mmm',ptr:q5},{name:'-62m',ptr:q6},{name:'6/mmm',ptr:q7}])
q0=ptr_new()
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'489: P 2 3',hash:'5843870D'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'490: F 2 3',hash:'93E38C71'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'491: I 2 3',hash:'DC4003C1'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'492: P 21 3',hash:'F9E6A645'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'493: I 21 3',hash:'7EC4457B'XL}])
q1=ptr_new([ $
    {name:'',ptr:r0},{name:'195: P 2 3',ptr:r1},{name:'196: F 2 3',ptr:r2},{name:'197: I 2 3',ptr:r3}, $
    {name:'198: P 21 3',ptr:r4},{name:'199: I 21 3',ptr:r5}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'494: P m -3',hash:'72E55913'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'495: P n -3:1',hash:'265CE726'XL},{name:'496: P n -3:2',hash:'D419D8A3'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'497: F m -3',hash:'58320B8D'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'498: F d -3:1',hash:'7DE7C89B'XL},{name:'499: F d -3:2',hash:'159CA8D3'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'500: I m -3',hash:'E23893DF'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'501: P a -3',hash:'1D5F4D3F'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'502: I a -3',hash:'7CE42B66'XL}])
q2=ptr_new([ $
    {name:'',ptr:r0},{name:'200: P m -3',ptr:r1},{name:'201: P n -3:1',ptr:r2},{name:'202: F m -3',ptr:r3}, $
    {name:'203: F d -3:1',ptr:r4},{name:'204: I m -3',ptr:r5},{name:'205: P a -3',ptr:r6},{name:'206: I a -3',ptr:r7}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'503: P 4 3 2',hash:'93A6EDEB'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'504: P 42 3 2',hash:'C55AE72A'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'505: F 4 3 2',hash:'25D65CF1'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'506: F 41 3 2',hash:'DDC15EF3'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'507: I 4 3 2',hash:'014B7ED2'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'508: P 43 3 2',hash:'4BCA2B34'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'509: P 41 3 2',hash:'0E761D2D'XL}])
r8=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'510: I 41 3 2',hash:'BB92B652'XL}])
q3=ptr_new([ $
    {name:'',ptr:r0},{name:'207: P 4 3 2',ptr:r1},{name:'208: P 42 3 2',ptr:r2},{name:'209: F 4 3 2',ptr:r3}, $
    {name:'210: F 41 3 2',ptr:r4},{name:'211: I 4 3 2',ptr:r5},{name:'212: P 43 3 2',ptr:r6},{name:'213: P 41 3 2',ptr:r7}, $
    {name:'214: I 41 3 2',ptr:r8}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'511: P -4 3 m',hash:'6FC31C7B'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'512: F -4 3 m',hash:'09ABDEBE'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'513: I -4 3 m',hash:'75D77C69'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'514: P -4 3 n',hash:'393F16BA'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'515: F -4 3 c',hash:'0968710F'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'516: I -4 3 d',hash:'E2ED982E'XL}])
q4=ptr_new([ $
    {name:'',ptr:r0},{name:'215: P -4 3 m',ptr:r1},{name:'216: F -4 3 m',ptr:r2},{name:'217: I -4 3 m',ptr:r3}, $
    {name:'218: P -4 3 n',ptr:r4},{name:'219: F -4 3 c',ptr:r5},{name:'220: I -4 3 d',ptr:r6}])
r0=ptr_new()
r1=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'517: P m -3 m',hash:'74C407D3'XL}])
r2=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'518: P n -3 n:1',hash:'C5E3DD9F'XL},{name:'519: P n -3 n:2',hash:'C7E69AB6'XL}])
r3=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'520: P m -3 n',hash:'BABD71C4'XL}])
r4=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'521: P n -3 m:1',hash:'0B9AAB88'XL},{name:'522: P n -3 m:2',hash:'DBDB6BFC'XL}])
r5=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'523: F m -3 m',hash:'5225B614'XL}])
r6=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'524: F m -3 c',hash:'481E9F10'XL}])
r7=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'525: F d -3 m:1',hash:'023A8184'XL},{name:'526: F d -3 m:2',hash:'38F5E458'XL}])
r8=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'527: F d -3 c:1',hash:'82B6E512'XL},{name:'528: F d -3 c:2',hash:'4D325AD3'XL}])
r9=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'529: I m -3 m',hash:'52045E84'XL}])
r10=ptr_new([ $
    {name:'',hash:'00000000'XL},{name:'530: I a -3 d',hash:'407DE7C1'XL}])
q5=ptr_new([ $
    {name:'',ptr:r0},{name:'221: P m -3 m',ptr:r1},{name:'222: P n -3 n:1',ptr:r2},{name:'223: P m -3 n',ptr:r3}, $
    {name:'224: P n -3 m:1',ptr:r4},{name:'225: F m -3 m',ptr:r5},{name:'226: F m -3 c',ptr:r6},{name:'227: F d -3 m:1',ptr:r7}, $
    {name:'228: F d -3 c:1',ptr:r8},{name:'229: I m -3 m',ptr:r9},{name:'230: I a -3 d',ptr:r10}])
p7=ptr_new([ $
    {name:'',ptr:q0},{name:'23',ptr:q1},{name:'m-3',ptr:q2},{name:'432',ptr:q3}, $
    {name:'-43m',ptr:q4},{name:'m-3m',ptr:q5}])
tree=[ $
    {name:'',ptr:p0},{name:'Triclinic',ptr:p1},{name:'Monoclinic',ptr:p2},{name:'Orthorhombic',ptr:p3}, $
    {name:'Tetragonal',ptr:p4},{name:'Trigonal',ptr:p5},{name:'Hexagonal',ptr:p6},{name:'Cubic',ptr:p7}]
return,tree
end;function Spacegrouptree
;%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

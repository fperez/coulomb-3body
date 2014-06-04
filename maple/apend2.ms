{VERSION 1 1 "MS-Windows" "1.0"}{GLOBALS 1 2}{FONT 0 "Arial" "Ari
al" "Helvetica-Bold" 7 10 0 "Helvetica-Bold" 10{COLOR "0000ff"}}
{FONT 1 "Times New Roman" "Times New Roman" "Times-Roman" 4 10 64
 "Times-Roman" 10{COLOR "000000"}}{FONT 2 "Courier New" "Courier 
New" "Courier" 4 10 192 "Courier" 10{COLOR "000000"}}{SCP_R 0 0 
43{SCP_R 1 0 25{COM_R 3 0{TEXT 1 0 ""}}{COM_R 6 0{TEXT 1 0 ""}}
{COM_R 9 0{TEXT 1 0 ""}}{COM_R 12 0{TEXT 1 0 ""}}{COM_R 15 0{TEXT
 1 0 ""}}{COM_R 18 0{TEXT 1 0 ""}}{COM_R 21 0{TEXT 1 0 ""}}{COM_R
 24 0{TEXT 1 0 ""}}{COM_R 27 0{TEXT 1 0 ""}}{COM_R 30 0{TEXT 1 0 
""}}{COM_R 33 0{TEXT 1 0 ""}}{COM_R 36 0{TEXT 1 0 ""}}{COM_R 39 0
{TEXT 1 0 ""}}{COM_R 42 0{TEXT 1 0 ""}}{COM_R 45 0{TEXT 1 0 ""}}
{COM_R 48 0{TEXT 1 0 ""}}{COM_R 51 0{TEXT 1 0 ""}}{COM_R 54 0
{TEXT 1 0 ""}}{COM_R 57 0{TEXT 1 0 ""}}{COM_R 60 0{TEXT 1 24 "ECU
ACIONES DE MOVIMIENTO"}}{COM_R 63 0{TEXT 1 26 "******************
********"}}{COM_R 66 0{TEXT 1 173 "Archivo en Maple V 2.0 para Wi
ndows con los principales c\341lculos de la generaci\363n de las 
ecuaciones de movimiento para el problema de 3 cuerpos en coorden
adas regularizadas."}}{COM_R 69 0{TEXT 1 159 "Aqu\355 hemos omiti
do mostrar la forma expl\355cita en C de las ecuaciones de movimi
ento, pues ellas se encuentran en el Ap\351ndice 1 en los archivo
s de c\363digo fuente."}}{COM_R 72 0{TEXT 1 0 ""}}{COM_R 75 0
{TEXT 1 62 "Comenzamos por cargar el paquete con rutinas de \341l
gebra lineal"}}}{SCP_R 77 0 2{INP_R 79 0 "> "{TEXT 0 13 "with(lin
alg):"}}{OUT_R 82 0 79{TEXT 2 72 "Warning: new definition for   n
orm\012Warning: new definition for   trace\012\012"}}}{SCP_R 84 0
 1{COM_R 86 0{TEXT 1 95 "y definimos por comodidad formas cortas 
para la transposici\363n y el m\363dulo cuadrado de un vector"}}}
{SCP_R 88 0 4{INP_R 90 0 "> "{TEXT 0 19 "alias(T=transpose);"}}
{INP_R 93 0 "> "{TEXT 0 22 "mc:=V->evalm(T(V)&*V);"}}{OUT_R 96 
8192 93{DAG ,3n3\`I`n3\`T`}}{OUT_R 99 8192 93{DAG :3n3\`mc`@7,2n3
\`V`,1,3n5\`operator`n4\`arrow`p7(3n4\`evalm`,2(3n3\`&*`,3(3n3\`T
`,2a2x0001p1Ep7}}}{SCP_R 101 0 1{COM_R 103 0{TEXT 1 63 "Ahora con
struimos vectores para nuestras coordenadas y momentos"}}}{SCP_R 
105 0 5{INP_R 107 0 "> "{TEXT 0 13 "Q:=vector(8):"}}{INP_R 110 0 
"> "{TEXT 0 13 "P:=vector(8):"}}{INP_R 113 0 "> "{TEXT 0 36 "P1:=
vector(4,[P[1],P[2],P[3],P[4]]):"}}{INP_R 116 0 "> "{TEXT 0 36 "P
2:=vector(4,[P[5],P[6],P[7],P[8]]):"}}{OUT_R 119 0 116{TEXT 2 1 "
\012"}}}{SCP_R 121 0 1{COM_R 123 0{TEXT 1 116 "Las funciones de t
ransformaci\363n de las coordenadas regularizadas (las transforma
ciones de Kustaanheimo-Steifel) son:"}}}{SCP_R 125 0 16{INP_R 127
 0 "> "{TEXT 0 34 "f[1]:=Q[1]^2-Q[2]^2-Q[3]^2+Q[4]^2;"}}{INP_R 
130 0 "> "{TEXT 0 30 "f[2]:=2*(Q[1]*Q[2]-Q[3]*Q[4]);"}}{INP_R 133
 0 "> "{TEXT 0 30 "f[3]:=2*(Q[1]*Q[3]+Q[2]*Q[4]);"}}{INP_R 136 0 
"> "{TEXT 0 8 "f[4]:=0;"}}{INP_R 139 0 "> "{TEXT 0 34 "f[5]:=Q[5]
^2-Q[6]^2-Q[7]^2+Q[8]^2;"}}{INP_R 142 0 "> "{TEXT 0 30 "f[6]:=2*(
Q[5]*Q[6]-Q[7]*Q[8]);"}}{INP_R 145 0 "> "{TEXT 0 30 "f[7]:=2*(Q[5
]*Q[7]+Q[6]*Q[8]);"}}{INP_R 148 0 "> "{TEXT 0 8 "f[8]:=0;"}}
{OUT_R 151 8192 148{DAG :3b3n3\`f`,2j2x0001+9*3b3n3\`Q`p4j2x0002p
5*3b3pA,2pDpDi2x0001*3b3pA,2j2x0003pDp16*3b3pA,2j2x0004pDp5}}
{OUT_R 154 8192 148{DAG :3b3n3\`f`,2j2x0002+5*5b3n3\`Q`,2j2x0001p
Db3pAp4pDp5*5b3pA,2j2x0003pDb3pA,2j2x0004pDi2x0002}}{OUT_R 157 
8192 148{DAG :3b3n3\`f`,2j2x0003+5*5b3n3\`Q`,2j2x0001pDb3pAp4pDj2
x0002*5b3pA,2p14pDb3pA,2j2x0004pDp14}}{OUT_R 160 8192 148{DAG :3b
3n3\`f`,2j2x0004j2x0000}}{OUT_R 163 8192 148{DAG :3b3n3\`f`,2j2x0
005+9*3b3n3\`Q`p4j2x0002j2x0001*3b3pA,2j2x0006pDi2x0001*3b3pA,2j2
x0007pDp18*3b3pA,2j2x0008pDpF}}{OUT_R 166 8192 148{DAG :3b3n3\`f`
,2j2x0006+5*5b3n3\`Q`,2j2x0005j2x0001b3pAp4pFj2x0002*5b3pA,2j2x00
07pFb3pA,2j2x0008pFi2x0002}}{OUT_R 169 8192 148{DAG :3b3n3\`f`,2j
2x0007+5*5b3n3\`Q`,2j2x0005j2x0001b3pAp4pFj2x0002*5b3pA,2j2x0006p
Fb3pA,2j2x0008pFp15}}{OUT_R 172 8192 148{DAG :3b3n3\`f`,2j2x0008j
2x0000}}}{SCP_R 174 0 1{COM_R 176 0{TEXT 1 54 "a partir de las cu
ales obtenemos las matrices A1 y A2:"}}}{SCP_R 178 0 8{INP_R 180 
0 "> "{TEXT 0 16 "A1:=matrix(4,4):"}}{INP_R 183 0 "> "{TEXT 0 20 
"for i from 1 to 4 do"}}{INP_R 186 0 "> "{TEXT 0 20 "for j from 1
 to 4 do"}}{INP_R 189 0 "> "{TEXT 0 25 "A1[i,j]:=diff(f[j],Q[i]);
"}}{INP_R 192 0 "> "{TEXT 0 3 "od;"}}{INP_R 195 0 "> "{TEXT 0 3 "
od;"}}{INP_R 198 0 "> "{TEXT 0 15 "'A1'=evalm(A1);"}}{OUT_R 201 
8192 198{DAG =3n3\`A1`(3n4\`MATRIX`,2[2,5[2,5+3b3n3\`Q`,2j2x0001j
2x0002+3b3pE,2p13p13+3b3pE,2j2x0003p13j2x0000[2,5+3p16i2x0002pC+3
b3pE,2j2x0004p13p22[2,5+3p1Cp28+3p2Cp28pCp22[2,5p2Bp35p15p22}}}
{SCP_R 203 0 8{INP_R 205 0 "> "{TEXT 0 16 "A2:=matrix(4,4):"}}
{INP_R 208 0 "> "{TEXT 0 20 "for i from 1 to 4 do"}}{INP_R 211 0 
"> "{TEXT 0 20 "for j from 1 to 4 do"}}{INP_R 214 0 "> "{TEXT 0 
29 "A2[i,j]:=diff(f[j+4],Q[i+4]);"}}{INP_R 217 0 "> "{TEXT 0 3 "o
d;"}}{INP_R 220 0 "> "{TEXT 0 3 "od;"}}{INP_R 223 0 "> "{TEXT 0 
15 "'A2'=evalm(A2);"}}{OUT_R 226 8192 223{DAG =3n3\`A2`(3n4\`MATR
IX`,2[2,5[2,5+3b3n3\`Q`,2j2x0005j2x0002+3b3pE,2j2x0006p13+3b3pE,2
j2x0007p13j2x0000[2,5+3p16i2x0002pC+3b3pE,2j2x0008p13p23[2,5+3p1D
p29+3p2Dp29pCp23[2,5p2Cp36p15p23}}}{SCP_R 228 0 1{COM_R 230 0
{TEXT 1 75 "Asignamos algunos nombres para tener formas m\341s le
gibles de los resultados:"}}}{SCP_R 232 0 9{INP_R 234 0 "> "{TEXT
 0 22 "A12:=evalm(A1&*T(A2)):"}}{INP_R 237 0 "> "{TEXT 0 22 "eval
m(T(P1)&*A12&*P2):"}}{INP_R 240 0 "> "{TEXT 0 17 "alias(P1A12P2=\
"):"}}{INP_R 243 0 "> "{TEXT 0 29 "alias(R1=sum(Q[r]^2,r=1..4)):"
}}{INP_R 246 0 "> "{TEXT 0 29 "alias(R2=sum(Q[r]^2,r=5..8)):"}}
{INP_R 249 0 "> "{TEXT 0 53 "alias(R=Rf(Q[1],Q[2],Q[3],Q[4],Q[5],
Q[6],Q[7],Q[8])):"}}{INP_R 252 0 "> "{TEXT 0 21 "alias(P1cuad=mc(
P1)):"}}{INP_R 255 0 "> "{TEXT 0 21 "alias(P2cuad=mc(P2));"}}
{OUT_R 258 8192 255{DAG ,9n3\`I`n3\`T`n4\`P1A12P2`n3\`R1`n3\`R2`n
3\`R`n4\`P1cuad`n4\`P2cuad`}}}{SCP_R 260 0 1{COM_R 262 0{TEXT 1 
64 "y definimos un vector que contendr\341 las ecuaciones de movi
miento"}}}{SCP_R 264 0 2{INP_R 266 0 "> "{TEXT 0 18 "EcHam:=vecto
r(16):"}}{OUT_R 269 0 266{TEXT 2 1 "\012"}}}{SCP_R 271 0 1{COM_R 
273 0{TEXT 1 48 "El Hamiltoniano en coordenadas regularizadas es:
"}}}{SCP_R 275 0 3{INP_R 277 0 "> "{TEXT 0 69 "Gamma:=(R2/(8*mu13
))*P1cuad+(R1/(8*mu23))*P2cuad+(16*m3)^(-1)*P1A12P2"}}{INP_R 280 
0 "> "{TEXT 0 40 "-Z2*Z3*R1-Z1*Z3*R2+R1*R2*(Z1*Z2/R-etot);"}}
{OUT_R 283 8192 280{DAG :3n4\`Gamma`+D*7n3\`R2`j2x0001n4\`mu13`i2
x0001n4\`P1cuad`p8/3p8j2x0008*7n3\`R1`p8n4\`mu23`pDn4\`P2cuad`p8p
13*5n3\`m3`pDn4\`P1A12P2`p8/3p8j2x0016*7n3\`Z2`p8n3\`Z3`p8p18p8pD
*7n3\`Z1`p8p34p8p6p8pD*7p18p8p6p8+5*7p3Bp8p31p8n3\`R`pDp8n4\`etot
`pDp8p8}}}{SCP_R 285 0 2{COM_R 287 0{TEXT 1 0 ""}}{COM_R 290 0
{TEXT 1 47 "Y de aqu\355 obtenemos las ecuaciones de Hamilton:"}}
}{SCP_R 292 0 5{INP_R 294 0 "> "{TEXT 0 20 "for i from 1 to 8 do"
}}{INP_R 297 0 "> "{TEXT 0 30 "EcHam[i]:=diff(Gamma,P[i]):od:"}}
{INP_R 300 0 "> "{TEXT 0 20 "for i from 1 to 8 do"}}{INP_R 303 0 
"> "{TEXT 0 33 "EcHam[i+8]:=-diff(Gamma,Q[i]):od:"}}{OUT_R 306 0 
303{TEXT 2 1 "\012"}}}{SCP_R 308 0 1{COM_R 310 0{TEXT 1 92 "Carga
mos la funci\363n C( ) necesaria para generar sintaxis de C de la
s expresiones calculadas:"}}}{SCP_R 312 0 2{INP_R 314 0 "> "{TEXT
 0 11 "readlib(C):"}}{OUT_R 317 0 314{TEXT 2 1 "\012"}}}{SCP_R 
319 0 1{COM_R 321 0{TEXT 1 80 "Para tener la forma m\341s compact
a posible de las ecuaciones, creamos la matriz A:"}}}{SCP_R 323 0
 6{INP_R 325 0 "> "{TEXT 0 27 "A:=array(sparse,1..8,1..8);"}}
{INP_R 328 0 "> "{TEXT 0 19 "copyinto(A1,A,1,1):"}}{INP_R 331 0 "
> "{TEXT 0 19 "copyinto(A2,A,5,5):"}}{INP_R 334 0 "> "{TEXT 0 12 
"A:=evalm(A);"}}{OUT_R 337 8192 334{DAG :3n3\`A`(3n4\`array`,5n4\
`sparse`~3j2x0001j2x0008pB[2,1}}{OUT_R 340 8192 334{DAG :3n3\`A`(
3n4\`MATRIX`,2[2,9[2,9+3b3n3\`Q`,2j2x0001j2x0002+3b3pE,2p13p13+3b
3pE,2j2x0003p13j2x0000p22p22p22p22[2,9+3p16i2x0002pC+3b3pE,2j2x00
04p13p22p22p22p22p22[2,9+3p1Cp2C+3p30p2CpCp22p22p22p22p22[2,9p2Fp
3Dp15p22p22p22p22p22[2,9p22p22p22p22+3b3pE,2j2x0005p13+3b3pE,2j2x
0006p13+3b3pE,2j2x0007p13p22[2,9p22p22p22p22+3p61p2Cp59+3b3pE,2j2
x0008p13p22[2,9p22p22p22p22+3p68p2C+3p7Ap2Cp59p22[2,9p22p22p22p22
p79p87p60p22}}}{SCP_R 342 0 8{COM_R 344 0{TEXT 1 14 "y el vector 
g:"}}{INP_R 347 0 "> "{TEXT 0 15 "g:=matrix(8,1):"}}{INP_R 350 0 
"> "{TEXT 0 20 "for i from 1 to 4 do"}}{INP_R 353 0 "> "{TEXT 0 
20 "g[i,1]:=f[i]-f[i+4]:"}}{INP_R 356 0 "> "{TEXT 0 18 "g[i+4,1]:
=-g[i,1]:"}}{INP_R 359 0 "> "{TEXT 0 3 "od:"}}{INP_R 362 0 "> "
{TEXT 0 13 "'g'=evalm(g);"}}{OUT_R 365 8192 362{DAG =3n3\`g`(3n4\
`MATRIX`,2[2,9[2,2+11*3b3n3\`Q`,2j2x0001j2x0002p12*3b3pF,2p14p14i
2x0001*3b3pF,2j2x0003p14p1D*3b3pF,2j2x0004p14p12*3b3pF,2j2x0005p1
4p1D*3b3pF,2j2x0006p14p12*3b3pF,2j2x0007p14p12*3b3pF,2j2x0008p14p
1D[2,2+9*5pEp12p18p12p14*5p20p12p28p12i2x0002*5p30p12p38p12p5D*5p
40p12p48p12p14[2,2+9*5pEp12p20p12p14*5p18p12p28p12p14*5p30p12p40p
12p5D*5p38p12p48p12p5D[2,2j2x0000[2,2+11pDp1Dp17p12p1Fp12p27p1Dp2
Fp12p37p1Dp3Fp1Dp47p12[2,2+9p52p5Dp58p14p5Fp14p65p5D[2,2+9p6Ep5Dp
74p5Dp7Ap14p80p14p86}}}{SCP_R 367 0 1{COM_R 369 0{TEXT 1 105 "Con
 esto podemos escribir las derivadas parciales de R con respecto 
a las coordenadas regularizadas como:"}}}{SCP_R 371 0 2{INP_R 373
 0 "> "{TEXT 0 24 "dRdQ:=evalm((1/R)*A&*g):"}}{OUT_R 376 0 373
{TEXT 2 1 "\012"}}}{SCP_R 378 0 1{COM_R 380 0{TEXT 1 31 "Definimo
s ciertas sustituciones"}}}{SCP_R 382 0 2{INP_R 384 0 "> "{TEXT 0
 50 "sustit:=P1cuad=P1c,P2cuad=P2c,R=R12,R1=Rr1,R2=Rr2;"}}{OUT_R 
387 8192 384{DAG :3n4\`sustit`,6=3n4\`P1cuad`n3\`P1c`=3n4\`P2cuad
`n3\`P2c`=3n3\`R`n3\`R12`=3n3\`R1`n3\`Rr1`=3n3\`R2`n3\`Rr2`}}}
{SCP_R 389 0 1{COM_R 391 0{TEXT 1 44 "que nos permiten simplifica
r las ecuaciones:"}}}{SCP_R 393 0 12{INP_R 395 0 "> "{TEXT 0 15 "
ec:=vector(16):"}}{INP_R 398 0 "> "{TEXT 0 20 "for i from 1 to 8 
do"}}{INP_R 401 0 "> "{TEXT 0 31 "ec[i]:=subs([sustit],EcHam[i]):
"}}{INP_R 404 0 "> "{TEXT 0 36 "dRdQ[i,1]:=subs([sustit],dRdQ[i,1
]):"}}{INP_R 407 0 "> "{TEXT 0 3 "od:"}}{INP_R 410 0 "> "{TEXT 0 
21 "for i from 9 to 16 do"}}{INP_R 413 0 "> "{TEXT 0 64 "ec[i]:=s
ubs([sustit,diff(R,Q[i-8])=`dRdQ[`.(i-8).`]`],EcHam[i]):"}}{INP_R
 416 0 "> "{TEXT 0 3 "od:"}}{INP_R 419 0 "> "{TEXT 0 22 "for i fr
om 1 to 16 do "}}{INP_R 422 0 "> "{TEXT 0 45 "ec[i]:=subsop(2=sim
plify(op(2,ec[i])),ec[i]):"}}{INP_R 425 0 "> "{TEXT 0 3 "od:"}}
{OUT_R 428 0 425{TEXT 2 1 "\012"}}}{SCP_R 430 0 20{COM_R 432 0
{TEXT 1 39 "Veamos las ecuaciones ya simplificadas:"}}{INP_R 435 
0 "> "{TEXT 0 21 "for i from 1 to 16 do"}}{INP_R 438 0 "> "{TEXT 
0 22 "`Ecuacion `.i = ec[i];"}}{INP_R 441 0 "> "{TEXT 0 3 "od;"}}
{OUT_R 444 8192 441{DAG =3n5\`Ecuacion~`1`+5*7n3\`Rr2`j2x0001n4\`
mu13`i2x0001b3n3\`P`,2p9p9/3p9j2x0004*5+19*7b3p11,2j2x0007p9b3n3\
`Q`p13p9b3p24p1Fp9p9*7p1Dp9b3p24,2j2x0002p9b3p24,2j2x0008p9p9*7p1
Dp9b3p24,2j2x0003p9b3p24,2j2x0005p9pE*7b3p11p38p9p23p9p36p9pE*7p4
Ep9p30p9p28p9p9*7p4Ep9p40p9b3p24,2j2x0006p9pE*7b3p11p48p9p23p9p46
p9pE*7p6Cp9p30p9p64p9pE*7p6Cp9p40p9p28p9pE*7b3p11p66p9p23p9p64p9p
9*7p86p9p30p9p46p9pE*7p86p9p40p9p36p9pEp9n3\`m3`pE/3pEp18}}{OUT_R
 447 8192 441{DAG =3n5\`Ecuacion~`2`+5*7n3\`Rr2`j2x0001n4\`mu13`i
2x0001b3n3\`P`,2j2x0002p9/3p9j2x0004*5+19*7b3p11,2j2x0007p9b3n3\`
Q`,2p9p9b3p25,2j2x0008p9pE*7p1Ep9b3p25p13p9b3p25p20p9p9*7p1Ep9b3p
25,2p19p9b3p25,2j2x0005p9p9*7b3p11p2Cp9p34p9p2Ap9pE*7p4Dp9p24p9p3
8p9pE*7p4Dp9p40p9b3p25,2j2x0006p9p9*7b3p11p47p9p34p9p45p9pE*7p6Bp
9p24p9p63p9p9*7p6Bp9p40p9p38p9p9*7b3p11p65p9p34p9p63p9p9*7p85p9p2
4p9p45p9p9*7p85p9p40p9p2Ap9p9p9n3\`m3`pEp17}}{OUT_R 450 8192 441
{DAG =3n5\`Ecuacion~`3`+5*7n3\`Rr2`j2x0001n4\`mu13`i2x0001b3n3\`P
`,2j2x0003p9/3p9j2x0004*5+19*7b3p11,2j2x0008p9b3n3\`Q`p13p9b3p25p
20p9pE*7b3p11,2j2x0007p9p24p9b3p25p31p9p9*7p2Fp9b3p25,2p19p9p29p9
p9*7p2Fp9b3p25,2p9p9b3p25,2j2x0005p9p9*7p1Ep9p3Fp9p37p9p9*7p1Ep9p
4Ap9b3p25,2j2x0006p9p9*7b3p11p51p9p24p9p4Fp9pE*7p6Bp9p3Fp9p63p9pE
*7p6Bp9p4Ap9p37p9p9*7b3p11p65p9p24p9p63p9p9*7p85p9p3Fp9p4Fp9pE*7p
85p9p4Ap9p29p9p9p9n3\`m3`pEp17}}{OUT_R 453 8192 441{DAG =3n5\`Ecu
acion~`4`+5*7n3\`Rr2`j2x0001n4\`mu13`i2x0001b3n3\`P`,2j2x0004p9/3
p9p14*5+19*7b3p11,2j2x0007p9b3n3\`Q`p13p9b3p24p1Fp9pE*7p1Dp9b3p24
,2j2x0003p9b3p24,2j2x0008p9p9*7p1Dp9b3p24,2j2x0002p9b3p24,2j2x000
5p9p9*7b3p11p38p9p23p9p36p9p9*7p4Ep9p30p9p28p9p9*7p4Ep9p40p9b3p24
,2j2x0006p9p9*7b3p11p48p9p23p9p46p9p9*7p6Cp9p30p9p64p9pE*7p6Cp9p4
0p9p28p9p9*7b3p11p66p9p23p9p64p9pE*7p86p9p30p9p46p9pE*7p86p9p40p9
p36p9p9p9n3\`m3`pEp17}}{OUT_R 456 8192 441{DAG =3n5\`Ecuacion~`5`
+5*7n3\`Rr1`j2x0001n4\`mu23`i2x0001b3n3\`P`,2j2x0005p9/3p9j2x0004
*5+19*7b3p11,2p9p9b3n3\`Q`p20p9b3p24p13p9pE*7p1Ep9b3p24,2j2x0002p
9b3p24,2j2x0006p9pE*7p1Ep9b3p24,2j2x0003p9b3p24,2j2x0007p9pE*7b3p
11p32p9p30p9p28p9p9*7p4Ep9p23p9p36p9pE*7p4Ep9b3p24,2p19p9p46p9pE*
7b3p11p42p9p40p9p28p9p9*7p6Bp9p62p9p36p9p9*7p6Bp9p23p9p46p9pE*7b3
p11p64p9p62p9p28p9pE*7p85p9p40p9p36p9p9*7p85p9p30p9p46p9pEp9n3\`m
3`pE/3pEp19}}{OUT_R 459 8192 441{DAG =3n5\`Ecuacion~`6`+5*7n3\`Rr
1`j2x0001n4\`mu23`i2x0001b3n3\`P`,2j2x0006p9/3p9j2x0004*5+19*7b3p
11,2p9p9b3n3\`Q`p20p9b3p24p13p9pE*7p1Ep9b3p24,2j2x0002p9b3p24,2j2
x0005p9p9*7p1Ep9b3p24,2j2x0003p9b3p24,2j2x0008p9p9*7b3p11p32p9p30
p9p28p9p9*7p4Ep9p23p9p36p9p9*7p4Ep9b3p24,2p19p9p46p9p9*7b3p11p42p
9p40p9p28p9p9*7p6Bp9p62p9p36p9pE*7p6Bp9p23p9p46p9p9*7b3p11p64p9p6
2p9p28p9pE*7p85p9p40p9p36p9pE*7p85p9p30p9p46p9p9p9n3\`m3`pEp17}}
{OUT_R 462 8192 441{DAG =3n5\`Ecuacion~`7`+5*7n3\`Rr1`j2x0001n4\`
mu23`i2x0001b3n3\`P`,2j2x0007p9/3p9j2x0004*5+19*7b3p11,2p9p9b3n3\
`Q`p20p9b3p24p13p9pE*7p1Ep9b3p24,2j2x0002p9b3p24,2j2x0008p9pE*7p1
Ep9b3p24,2j2x0003p9b3p24,2j2x0005p9p9*7b3p11p32p9p30p9p28p9p9*7p4
Ep9p23p9p36p9pE*7p4Ep9b3p24,2p19p9p46p9p9*7b3p11p42p9p40p9p28p9p9
*7p6Bp9p62p9p36p9p9*7p6Bp9p23p9p46p9p9*7b3p11p64p9p62p9p28p9pE*7p
85p9p40p9p36p9p9*7p85p9p30p9p46p9p9p9n3\`m3`pEp17}}{OUT_R 465 
8192 441{DAG =3n5\`Ecuacion~`8`+5*7n3\`Rr1`j2x0001n4\`mu23`i2x000
1b3n3\`P`,2j2x0008p9/3p9j2x0004*5+19*7b3p11,2p9p9b3n3\`Q`p20p9b3p
24p13p9p9*7p1Ep9b3p24,2j2x0002p9b3p24,2j2x0007p9pE*7p1Ep9b3p24,2j
2x0003p9b3p24,2j2x0006p9p9*7b3p11p32p9p30p9p28p9pE*7p4Ep9p23p9p36
p9pE*7p4Ep9b3p24,2p19p9p46p9p9*7b3p11p42p9p40p9p28p9pE*7p6Bp9p62p
9p36p9p9*7p6Bp9p23p9p46p9p9*7b3p11p64p9p62p9p28p9p9*7p85p9p40p9p3
6p9p9*7p85p9p30p9p46p9p9p9n3\`m3`pEp17}}{OUT_R 468 8192 441{DAG =
3n5\`Ecuacion~`9`+B*7b3n3\`Q`,2j2x0001pBn4\`mu23`i2x0001n3\`P2c`p
B/3p11j2x0004*5+19*7b3n3\`P`,2j2x0007pBb3p1E,2j2x0002pBb3p8,2j2x0
008pBp11*7p1DpBb3p1EpApBb3p8p20pBp11*7p1DpBb3p1E,2j2x0003pBb3p8,2
j2x0005pBpB*7b3p1Ep2CpBp34pBp2ApBpB*7p4EpBp24pBp38pBp11*7p4EpBp40
pBb3p8,2j2x0006pBpB*7b3p1Ep48pBp34pBp46pBpB*7p6CpBp24pBp64pBpB*7p
6CpBp40pBp38pBpB*7b3p1Ep66pBp34pBp64pBp11*7p86pBp24pBp46pBpB*7p86
pBp40pBp2ApBpBpBn3\`m3`p11p16*7n3\`Z2`pBn3\`Z3`pBp7pBp27*7p7pBn3\
`Rr2`pB+5*7n3\`Z1`pBpA5pBn3\`R12`p11pBn4\`etot`p11pBi2x0002*Dn3\`
Rr1`pBpB1pBpB6pBpA5pBpBBpC4n4\`dRdQ[1]`pBpB}}{OUT_R 471 8192 441
{DAG =3n5\`Ecuacion~`10`+B*7b3n3\`Q`,2j2x0002j2x0001n4\`mu23`i2x0
001n3\`P2c`pD/3p12j2x0004*5+19*7b3n3\`P`,2j2x0007pDb3p1F,2pDpDb3p
8,2j2x0008pDp12*7p1EpDb3p1FpApDb3p8p21pDpD*7p1EpDb3p1F,2p19pDb3p8
,2j2x0005pDpD*7b3p1Fp2CpDp25pDp38pDp12*7p4DpDp34pDp2ApDp12*7p4DpD
p40pDb3p8,2j2x0006pDpD*7b3p1Fp47pDp25pDp63pDpD*7p6BpDp34pDp45pDp1
2*7p6BpDp40pDp38pDpD*7b3p1Fp65pDp25pDp45pDpD*7p85pDp34pDp63pDpD*7
p85pDp40pDp2ApDpDpDn3\`m3`p12p17*7n3\`Z2`pDn3\`Z3`pDp7pDpB*7p7pDn
3\`Rr2`pD+5*7n3\`Z1`pDpA4pDn3\`R12`p12pDn4\`etot`p12pDi2x0002*Dn3
\`Rr1`pDpB0pDpB5pDpA4pDpBApC3n4\`dRdQ[2]`pDpD}}{OUT_R 474 8192 
441{DAG =3n5\`Ecuacion~`11`+B*7b3n3\`Q`,2j2x0003j2x0001n4\`mu23`i
2x0001n3\`P2c`pD/3p12j2x0004*5+19*7b3n3\`P`,2j2x0008pDb3p1FpApDb3
p8p21pDp12*7b3p1F,2j2x0007pDb3p1F,2pDpDb3p8,2j2x0005pDpD*7p2FpDp2
5pDb3p8p31pDpD*7p2FpDb3p1F,2p19pDp29pDpD*7p1EpDp35pDb3p8,2j2x0006
pDpD*7p1EpDp4EpDp46pDpD*7b3p1Fp3CpDp35pDp46pDpD*7p6BpDp25pDp3ApDp
12*7p6BpDp4EpDp5BpDp12*7b3p1Fp5DpDp35pDp29pDpD*7p85pDp25pDp5BpDpD
*7p85pDp4EpDp3ApDp12pDn3\`m3`p12p17*7n3\`Z2`pDn3\`Z3`pDp7pDj2x000
2*7p7pDn3\`Rr2`pD+5*7n3\`Z1`pDpA4pDn3\`R12`p12pDn4\`etot`p12pDi2x
0002*Dn3\`Rr1`pDpB1pDpB6pDpA4pDpBBpC4n4\`dRdQ[3]`pDpD}}{OUT_R 477
 8192 441{DAG =3n5\`Ecuacion~`12`+B*7b3n3\`Q`,2j2x0004j2x0001n4\`
mu23`i2x0001n3\`P2c`pD/3p12pB*5+19*7b3n3\`P`,2j2x0007pDb3p1E,2j2x
0002pDb3p8,2j2x0005pDpD*7p1DpDb3p1E,2j2x0003pDb3p8,2j2x0008pDpD*7
p1DpDb3p1EpApDb3p8p20pDp12*7b3p1Ep3CpDp24pDb3p8,2j2x0006pDpD*7p4E
pDp34pDp48pDpD*7p4EpDp44pDp3ApDpD*7b3p1Ep2CpDp24pDp48pDpD*7p6CpDp
34pDp54pDp12*7p6CpDp44pDp2ApDpD*7b3p1Ep56pDp24pDp3ApDpD*7p86pDp34
pDp2ApDp12*7p86pDp44pDp54pDp12pDn3\`m3`p12p17*7n3\`Z2`pDn3\`Z3`pD
p7pDp27*7p7pDn3\`Rr2`pD+5*7n3\`Z1`pDpA5pDn3\`R12`p12pDn4\`etot`p1
2pDi2x0002*Dn3\`Rr1`pDpB1pDpB6pDpA5pDpBBpC4n4\`dRdQ[4]`pDpD}}
{OUT_R 480 8192 441{DAG =3n5\`Ecuacion~`13`+B*7b3n3\`Q`,2j2x0005j
2x0001n4\`mu13`i2x0001n3\`P1c`pD/3p12j2x0004*5+19*7b3n3\`P`,2j2x0
007pDb3p1F,2pDpDb3p8,2j2x0003pDpD*7p1EpDb3p1F,2j2x0002pDb3p8,2p19
pDpD*7p1EpDb3p1Fp2CpDb3p8p27pDpD*7p1EpDb3p1Fp3CpDb3p8p36pDpD*7b3p
1FpApDp25pDp47pDpD*7p59pDp34pDp53pDp12*7p59pDp43pDp2ApDp12*7p59pD
p4FpDp3ApDpD*7b3p1F,2j2x0006pDp25pDp53pDpD*7p7BpDp34pDp47pDpD*7p7
BpDp43pDp3ApDp12*7p7BpDp4FpDp2ApDp12pDn3\`m3`p12p17*7n3\`Z1`pDn3\
`Z3`pDp7pDp37*7n3\`Rr1`pDp7pD+5*7pA4pDn3\`Z2`pDn3\`R12`p12pDn4\`e
tot`p12pDi2x0002*DpAEpDn3\`Rr2`pDpA4pDpB7pDpBApC3n4\`dRdQ[5]`pDpD
}}{OUT_R 483 8192 441{DAG =3n5\`Ecuacion~`14`+B*7b3n3\`Q`,2j2x000
6j2x0001n4\`mu13`i2x0001n3\`P1c`pD/3p12j2x0004*5+19*7b3n3\`P`,2j2
x0008pDb3p1F,2pDpDb3p8,2j2x0003pDpD*7p1EpDb3p1F,2j2x0002pDb3p8,2p
19pDpD*7p1EpDb3p1Fp2CpDb3p8p27pDpD*7p1EpDb3p1Fp3CpDb3p8p36pDpD*7b
3p1F,2j2x0005pDp25pDp53pDpD*7p59pDp34pDp47pDpD*7p59pDp43pDp3ApDp1
2*7p59pDp4FpDp2ApDp12*7b3p1FpApDp25pDp47pDp12*7p7DpDp34pDp53pDpD*
7p7DpDp43pDp2ApDpD*7p7DpDp4FpDp3ApDp12pDn3\`m3`p12p17*7n3\`Z1`pDn
3\`Z3`pDp7pDp37*7n3\`Rr1`pDp7pD+5*7pA4pDn3\`Z2`pDn3\`R12`p12pDn4\
`etot`p12pDi2x0002*DpAEpDn3\`Rr2`pDpA4pDpB7pDpBApC3n4\`dRdQ[6]`pD
pD}}{OUT_R 486 8192 441{DAG =3n5\`Ecuacion~`15`+B*7b3n3\`Q`,2j2x0
007j2x0001n4\`mu13`i2x0001n3\`P1c`pD/3p12j2x0004*5+19*7b3n3\`P`pA
pDb3p1F,2pDpDb3p8p25pDp12*7p1EpDb3p1F,2j2x0002pDb3p8p32pDpD*7p1Ep
Db3p1F,2j2x0003pDb3p8p40pDpD*7p1EpDb3p1F,2p19pDb3p8p4EpDp12*7b3p1
F,2j2x0008pDp23pDp36pDp12*7p57pDp30pDp28pDp12*7p57pDp3EpDp51pDpD*
7p57pDp4CpDp44pDpD*7b3p1F,2j2x0005pDp23pDp44pDpD*7p7BpDp30pDp51pD
pD*7p7BpDp3EpDp28pDpD*7p7BpDp4CpDp36pDpDpDn3\`m3`p12p17*7n3\`Z1`p
Dn3\`Z3`pDp7pDp33*7n3\`Rr1`pDp7pD+5*7pA4pDn3\`Z2`pDn3\`R12`p12pDn
4\`etot`p12pDi2x0002*DpAEpDn3\`Rr2`pDpA4pDpB7pDpBApC3n4\`dRdQ[7]`
pDpD}}{OUT_R 489 8192 441{DAG =3n5\`Ecuacion~`16`+B*7b3n3\`Q`,2j2
x0008j2x0001n4\`mu13`i2x0001n3\`P1c`pD/3p12j2x0004*5+19*7b3n3\`P`
pApDb3p1F,2j2x0003pDb3p8p25pDpD*7b3p1F,2j2x0007pDb3p1F,2j2x0002pD
b3p8,2pDpDpD*7p2FpDb3p1Fp3DpDb3p8p37pDpD*7p2FpDp23pDb3p8,2p19pDp1
2*7p2FpDb3p1Fp54pDp29pDp12*7p1EpDp44pDp3BpDp12*7p1EpDp35pDp48pDpD
*7p1EpDp5BpDp52pDp12*7b3p1F,2j2x0006pDp44pDp29pDp12*7p7BpDp35pDp5
2pDp12*7p7BpDp23pDp3BpDp12*7p7BpDp5BpDp48pDp12pDn3\`m3`p12/3pDp19
*7n3\`Z1`pDn3\`Z3`pDp7pDp38*7n3\`Rr1`pDp7pD+5*7pA6pDn3\`Z2`pDn3\`
R12`p12pDn4\`etot`p12pDi2x0002*DpB0pDn3\`Rr2`pDpA6pDpB9pDpBCpC5n4
\`dRdQ[8]`pDpD}}}{SCP_R 491 0 1{COM_R 493 0{TEXT 1 109 "En las ex
presiones anteriores, las derivadas parciales de R con respecto a
 las coordenadas regularizadas son:"}}}{SCP_R 495 0 11{INP_R 497 
0 "> "{TEXT 0 20 "for i from 1 to 8 do"}}{INP_R 500 0 "> "{TEXT 0
 27 "`dRdQ[ `.i.`]` = dRdQ[i,1];"}}{INP_R 503 0 "> "{TEXT 0 3 "od
;"}}{OUT_R 506 8192 503{DAG =3n5\`dRdQ[~`1]`+3*5+19*3b3n3\`Q`,2j2
x0001j2x0003pD*5p9pDb3pA,2j2x0002p18pD*5p9pDb3pA,2pFp18pD*5p9pDb3
pA,2j2x0004p18pD*5p9pDb3pA,2j2x0005p18i2x0001*5p9pDb3pA,2j2x0006p
18pD*5p9pDb3pA,2j2x0007p18pD*5p9pDb3pA,2j2x0008p18p38*7p15pDp32pD
p3DpDi2x0002*7p15pDp47pDp51pDp18*7p1FpDp32pDp47pDp5F*7p1FpDp3DpDp
51pDp5FpDn3\`R12`p38p18}}{OUT_R 509 8192 503{DAG =3n5\`dRdQ[~`2]`
+3*5+19*5b3n3\`Q`,2j2x0002j2x0001b3pA,2pFpDpF*3p9j2x0003pF*5p9pFb
3pA,2p19pDpF*5p9pFb3pA,2j2x0004pDpF*5p9pFb3pA,2j2x0005pDpF*5p9pFb
3pA,2j2x0006pDi2x0001*5p9pFb3pA,2j2x0007pDp42*5p9pFb3pA,2j2x0008p
DpF*7p11pFp32pFp3CpFi2x0002*7p11pFp47pFp51pFpD*7p28pFp32pFp47pFp5
F*7p28pFp3CpFp51pFp5FpFn3\`R12`p42pD}}{OUT_R 512 8192 503{DAG =3n
5\`dRdQ[~`3]`+3*5+19*5b3n3\`Q`,2j2x0003j2x0001b3pA,2pFj2x0002pF*5
p9pFb3pA,2p15p15pF*3p9pDpF*5p9pFb3pA,2j2x0004p15pF*5p9pFb3pA,2j2x
0005p15pF*5p9pFb3pA,2j2x0006p15i2x0001*5p9pFb3pA,2j2x0007p15p42*5
p9pFb3pA,2j2x0008p15pF*7p28pFp32pFp3CpFp15*7p28pFp47pFp51pFi2x000
2*7p11pFp32pFp47pFp67*7p11pFp3CpFp51pFp67pFn3\`R12`p42p15}}{OUT_R
 515 8192 503{DAG =3n5\`dRdQ[~`4]`+3*5+19*5b3n3\`Q`,2j2x0004j2x00
01b3pA,2pFj2x0002pF*5p9pFb3pA,2p15p15pF*5p9pFb3pA,2j2x0003p15pF*3
p9p27pF*5p9pFb3pA,2j2x0005p15i2x0001*5p9pFb3pA,2j2x0006p15pF*5p9p
Fb3pA,2j2x0007p15pF*5p9pFb3pA,2j2x0008p15p38*7p24pFp32pFp3DpFp15*
7p24pFp47pFp51pFi2x0002*7p1BpFp32pFp47pFp67*7p1BpFp3DpFp51pFp67pF
n3\`R12`p38p15}}{OUT_R 518 8192 503{DAG =3n5\`dRdQ[~`5]`+3*5+19*5
b3n3\`Q`,2j2x0005j2x0001b3pA,2pFj2x0002i2x0001*5p9pFb3pA,2p15p15p
F*5p9pFb3pA,2j2x0003p15pF*5p9pFb3pA,2j2x0004p15p17*3p9p28pF*5p9pF
b3pA,2j2x0006p15pF*5p9pFb3pA,2j2x0007p15pF*5p9pFb3pA,2j2x0008p15p
F*7p3DpFp11pFp1CpFi2x0002*7p3DpFp25pFp2FpFp15*7p47pFp11pFp25pFp5F
*7p47pFp1CpFp2FpFp5FpFn3\`R12`p17p15}}{OUT_R 521 8192 503{DAG =3n
5\`dRdQ[~`6]`+3*5+19*5b3n3\`Q`,2j2x0006j2x0001b3pA,2pFj2x0002i2x0
001*5p9pFb3pA,2p15p15pF*5p9pFb3pA,2j2x0003p15pF*5p9pFb3pA,2j2x000
4p15p17*5p9pFb3pA,2j2x0005p15p17*3p9p28p17*5p9pFb3pA,2j2x0007p15p
17*5p9pFb3pA,2j2x0008p15p17*7p39pFp11pFp1CpFp15*7p39pFp25pFp2FpFi
2x0002*7p51pFp11pFp25pFp15*7p51pFp1CpFp2FpFp15pFn3\`R12`p17p67}}
{OUT_R 524 8192 503{DAG =3n5\`dRdQ[~`7]`+3*5+19*5b3n3\`Q`,2j2x000
7j2x0001b3pA,2pFj2x0002i2x0001*5p9pFb3pA,2p15p15pF*5p9pFb3pA,2j2x
0003p15pF*5p9pFb3pA,2j2x0004p15p17*5p9pFb3pA,2j2x0005p15p17*5p9pF
b3pA,2j2x0006p15p17*3p9p28p17*5p9pFb3pA,2j2x0008p15p17*7p51pFp11p
Fp1CpFi2x0002*7p51pFp25pFp2FpFp15*7p39pFp11pFp25pFp15*7p39pFp1CpF
p2FpFp15pFn3\`R12`p17p5F}}{OUT_R 527 8192 503{DAG =3n5\`dRdQ[~`8]
`+3*5+19*5b3n3\`Q`,2j2x0008j2x0001b3pA,2pFj2x0002i2x0001*5p9pFb3p
A,2p15p15pF*5p9pFb3pA,2j2x0003p15pF*5p9pFb3pA,2j2x0004p15p17*5p9p
Fb3pA,2j2x0005p15pF*5p9pFb3pA,2j2x0006p15pF*5p9pFb3pA,2j2x0007p15
pF*3p9p28pF*7p4DpFp11pFp1CpFp15*7p4DpFp25pFp2FpFi2x0002*7p43pFp11
pFp25pFp67*7p43pFp1CpFp2FpFp67pFn3\`R12`p17p15}}}{SCP_R 529 0 2
{COM_R 531 0{TEXT 1 0 ""}}{COM_R 534 0{TEXT 1 65 "Finalmente, sim
plifiquemos la forma de R en coord. regularizadas:"}}}{SCP_R 536 
0 2{INP_R 538 0 "> "{TEXT 0 21 "sumR:=r->f[r]-f[r+4];"}}{OUT_R 
541 8192 538{DAG :3n4\`sumR`@7,2n3\`r`,1,3n5\`operator`n4\`arrow`
p8+5b3n3\`f`,2a2x0001j2x0001b3p14,2+5p17p19j2x0004p19i2x0001p8}}}
{SCP_R 543 0 2{INP_R 545 0 "> "{TEXT 0 25 "sumRe:=expand(sumR(r)^
2);"}}{OUT_R 548 8192 545{DAG :3n4\`sumRe`+7*3b3n3\`f`,2n3\`r`j2x
0002j2x0001*5p6pEb3p7,2+5pApEj2x0004pEpEi2x0002*3p13pCpE}}}{SCP_R
 550 0 6{INP_R 552 0 "> "{TEXT 0 5 "j:=1:"}}{INP_R 555 0 "> "
{TEXT 0 17 "for i in sumRe do"}}{INP_R 558 0 "> "{TEXT 0 12 "sumR
2[j]:=i:"}}{INP_R 561 0 "> "{TEXT 0 7 "j:=j+1:"}}{INP_R 564 0 "> 
"{TEXT 0 3 "od:"}}{OUT_R 567 0 564{TEXT 2 1 "\012"}}}{SCP_R 569 0
 5{INP_R 571 0 "> "{TEXT 0 20 "for i from 1 to 3 do"}}{INP_R 574 
0 "> "{TEXT 0 32 "Rf[i]:=sum(sumR2[i],r=1..4); od;"}}{OUT_R 577 
8192 574{DAG :3b3n3\`Rf`,2j2x0001+7*3+9*3b3n3\`Q`p4j2x0002p5*3b3p
C,2pFpFi2x0001*3b3pC,2j2x0003pFp18*3b3pC,2j2x0004pFp5pFp5*3+5*5pB
p5p13p5pF*5p1Bp5p23p5i2x0002pFp5*3+5*5pBp5p1Bp5pF*5p13p5p23p5pFpF
p5}}{OUT_R 580 8192 574{DAG :3b3n3\`Rf`,2j2x0002+7*5+9*3b3n3\`Q`,
2j2x0001p5pF*3b3pCp4p5i2x0001*3b3pC,2j2x0003p5p18*3b3pC,2j2x0004p
5pFpF+9*3b3pC,2j2x0005p5pF*3b3pC,2j2x0006p5p18*3b3pC,2j2x0007p5p1
8*3b3pC,2j2x0008p5pFpFi2x0002*5+5*5pBpFp14pFp5*5p1BpFp23pFp4DpF+5
*5p2DpFp35pFp5*5p3DpFp45pFp4DpFp4D*5+5*5pBpFp1BpFp5*5p14pFp23pFp5
pF+5*5p2DpFp3DpFp5*5p35pFp45pFp5pFp4D}}{OUT_R 583 8192 574{DAG :3
b3n3\`Rf`,2j2x0003+7*3+9*3b3n3\`Q`,2j2x0005j2x0002j2x0001*3b3pC,2
j2x0006p11i2x0001*3b3pC,2j2x0007p11p1C*3b3pC,2j2x0008p11p13p11p13
*3+5*5pBp13p16p13p11*5p1Fp13p27p13i2x0002p11p13*3+5*5pBp13p1Fp13p
11*5p16p13p27p13p11p11p13}}}{SCP_R 585 0 6{INP_R 587 0 "> "{TEXT 
0 24 "Rfsim[1]:=factor(Rf[1]);"}}{INP_R 590 0 "> "{TEXT 0 16 "Rfs
im[2]:=Rf[2];"}}{INP_R 593 0 "> "{TEXT 0 24 "Rfsim[3]:=factor(Rf[
3]);"}}{OUT_R 596 8192 593{DAG :3b3n4\`Rfsim`,2j2x0001*3n3\`R1`j2
x0002}}{OUT_R 599 8192 593{DAG :3b3n4\`Rfsim`,2j2x0002+7*5+9*3b3n
3\`Q`,2j2x0001p6p10*3b3pDp5p6i2x0001*3b3pD,2j2x0003p6p19*3b3pD,2j
2x0004p6p10p10+9*3b3pD,2j2x0005p6p10*3b3pD,2j2x0006p6p19*3b3pD,2j
2x0007p6p19*3b3pD,2j2x0008p6p10p10i2x0002*5+5*5pCp10p15p10p6*5p1C
p10p24p10p4Ep10+5*5p2Ep10p36p10p6*5p3Ep10p46p10p4Ep10p4E*5+5*5pCp
10p1Cp10p6*5p15p10p24p10p6p10+5*5p2Ep10p3Ep10p6*5p36p10p46p10p6p1
0p4E}}{OUT_R 602 8192 593{DAG :3b3n4\`Rfsim`,2j2x0003*3n3\`R2`j2x
0002}}}{SCP_R 604 0 2{INP_R 606 0 "> "{TEXT 0 39 "Rfun:=(sum(Rfsi
m['i'],'i'=1..3))^(1/2);"}}{OUT_R 609 8192 606{DAG :3n4\`Rfun`*3+
B*3n3\`R1`j2x0002j2x0001*5+9*3b3n3\`Q`,2pBp9pB*3b3p11,2p9p9i2x000
1*3b3p11,2j2x0003p9p1D*3b3p11,2j2x0004p9pBpB+9*3b3p11,2j2x0005p9p
B*3b3p11,2j2x0006p9p1D*3b3p11,2j2x0007p9p1D*3b3p11,2j2x0008p9pBpB
i2x0002*5+5*5p10pBp18pBp9*5p20pBp28pBp52pB+5*5p32pBp3ApBp9*5p42pB
p4ApBp52pBp52*5+5*5p10pBp20pBp9*5p18pBp28pBp9pB+5*5p32pBp42pBp9*5
p3ApBp4ApBp9pBp52*3n3\`R2`p9pB/3pBp9}}}{SCP_R 611 0 2{INP_R 613 0
 "> "{TEXT 0 24 "Rfun:=subs(sustit,Rfun);"}}{OUT_R 616 8192 613
{DAG :3n4\`Rfun`*3+B*3n3\`Rr1`j2x0002j2x0001*5+9*3b3n3\`Q`,2pBp9p
B*3b3p11,2p9p9i2x0001*3b3p11,2j2x0003p9p1D*3b3p11,2j2x0004p9pBpB+
9*3b3p11,2j2x0005p9pB*3b3p11,2j2x0006p9p1D*3b3p11,2j2x0007p9p1D*3
b3p11,2j2x0008p9pBpBi2x0002*5+5*5p10pBp18pBp9*5p20pBp28pBp52pB+5*
5p32pBp3ApBp9*5p42pBp4ApBp52pBp52*5+5*5p10pBp20pBp9*5p18pBp28pBp9
pB+5*5p32pBp42pBp9*5p3ApBp4ApBp9pBp52*3n3\`Rr2`p9pB/3pBp9}}}
{SCP_R 618 0 1{COM_R 620 0{TEXT 1 68 "Y con esta forma tenemos un
a expresi\363n m\341s compacta que pasamos a C:"}}}{SCP_R 622 0 1
{INP_R 624 0 "> "{TEXT 0 8 "C(Rfun):"}}}}{END}

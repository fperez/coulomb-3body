{VERSION 1 1 "MS-Windows" "1.0"}{GLOBALS 1 2}{FONT 0 "Arial" "Ari
al" "Helvetica-Bold" 7 10 0 "Helvetica-Bold" 10{COLOR "0000ff"}}
{FONT 1 "Times New Roman" "Times New Roman" "Times-Roman" 4 10 64
 "Times-Roman" 10{COLOR "000000"}}{FONT 2 "Courier New" "Courier 
New" "Courier" 4 10 192 "Courier" 10{COLOR "000000"}}{SCP_R 0 0 
35{SCP_R 1 0 2{COM_R 3 0{TEXT 1 62 "Ecuaciones de transformaci\36
3n de las coordenadas regularizadas:"}}{COM_R 6 0{TEXT 1 49 "****
*********************************************"}}}{SCP_R 8 0 2
{INP_R 10 0 "> "{TEXT 0 13 "with(linalg):"}}{OUT_R 13 0 10{TEXT 2
 72 "Warning: new definition for   norm\012Warning: new definitio
n for   trace\012\012"}}}{SCP_R 15 0 2{INP_R 17 0 "> "{TEXT 0 19 
"alias(T=transpose);"}}{OUT_R 20 8192 17{DAG ,3n3\`I`n3\`T`}}}
{SCP_R 22 0 2{INP_R 24 0 "> "{TEXT 0 27 "mc:=V->evalm(T(V)&*V)[1,
1];"}}{OUT_R 27 8192 24{DAG :3n3\`mc`@7,2n3\`V`,1,3n5\`operator`n
4\`arrow`p7b3(3n4\`evalm`,2(3n3\`&*`,3(3n3\`T`,2a2x0001p1F,3j2x00
01p23p7}}}{SCP_R 29 0 16{INP_R 31 0 "> "{TEXT 0 34 "f[1]:=Q[1]^2-
Q[2]^2-Q[3]^2+Q[4]^2;"}}{INP_R 34 0 "> "{TEXT 0 30 "f[2]:=2*(Q[1]
*Q[2]-Q[3]*Q[4]);"}}{INP_R 37 0 "> "{TEXT 0 30 "f[3]:=2*(Q[1]*Q[3
]+Q[2]*Q[4]);"}}{INP_R 40 0 "> "{TEXT 0 8 "f[4]:=0;"}}{INP_R 43 0
 "> "{TEXT 0 34 "f[5]:=Q[5]^2-Q[6]^2-Q[7]^2+Q[8]^2;"}}{INP_R 46 0
 "> "{TEXT 0 30 "f[6]:=2*(Q[5]*Q[6]-Q[7]*Q[8]);"}}{INP_R 49 0 "> 
"{TEXT 0 30 "f[7]:=2*(Q[5]*Q[7]+Q[6]*Q[8]);"}}{INP_R 52 0 "> "
{TEXT 0 8 "f[8]:=0;"}}{OUT_R 55 8192 52{DAG :3b3n3\`f`,2j2x0001+9
*3b3n3\`Q`p4j2x0002p5*3b3pA,2pDpDi2x0001*3b3pA,2j2x0003pDp16*3b3p
A,2j2x0004pDp5}}{OUT_R 58 8192 52{DAG :3b3n3\`f`,2j2x0002+5*5b3n3
\`Q`,2j2x0001pDb3pAp4pDp5*5b3pA,2j2x0003pDb3pA,2j2x0004pDi2x0002}
}{OUT_R 61 8192 52{DAG :3b3n3\`f`,2j2x0003+5*5b3n3\`Q`,2j2x0001pD
b3pAp4pDj2x0002*5b3pA,2p14pDb3pA,2j2x0004pDp14}}{OUT_R 64 8192 52
{DAG :3b3n3\`f`,2j2x0004j2x0000}}{OUT_R 67 8192 52{DAG :3b3n3\`f`
,2j2x0005+9*3b3n3\`Q`p4j2x0002j2x0001*3b3pA,2j2x0006pDi2x0001*3b3
pA,2j2x0007pDp18*3b3pA,2j2x0008pDpF}}{OUT_R 70 8192 52{DAG :3b3n3
\`f`,2j2x0006+5*5b3n3\`Q`,2j2x0005j2x0001b3pAp4pFj2x0002*5b3pA,2j
2x0007pFb3pA,2j2x0008pFi2x0002}}{OUT_R 73 8192 52{DAG :3b3n3\`f`,
2j2x0007+5*5b3n3\`Q`,2j2x0005j2x0001b3pAp4pFj2x0002*5b3pA,2j2x000
6pFb3pA,2j2x0008pFp15}}{OUT_R 76 8192 52{DAG :3b3n3\`f`,2j2x0008j
2x0000}}}{SCP_R 78 0 6{INP_R 80 0 "> "{TEXT 0 39 "R1:=simplify(sq
rt(sum(f[r]^2,r=1..4)));"}}{INP_R 83 0 "> "{TEXT 0 39 "R2:=simpli
fy(sqrt(sum(f[r]^2,r=5..8)));"}}{INP_R 86 0 "> "{TEXT 0 37 "R:=sq
rt(sum((f[r]-f[r+4])^2,r=1..4));"}}{OUT_R 89 8192 86{DAG :3n3\`R1
`+9*3b3n3\`Q`,2j2x0001j2x0002p9*3b3p6,2pBpBp9*3b3p6,2j2x0003pBp9*
3b3p6,2j2x0004pBp9}}{OUT_R 92 8192 86{DAG :3n3\`R2`+9*3b3n3\`Q`,2
j2x0005j2x0002j2x0001*3b3p6,2j2x0006pBpD*3b3p6,2j2x0007pBpD*3b3p6
,2j2x0008pBpD}}{OUT_R 95 8192 86{DAG :3n3\`R`*3+59*5b3n3\`Q`,2j2x
0001j2x0002b3p7,2pCpCpC*5pEpCb3p7,2j2x0003pCpC*9pEpAb3p7,2j2x0004
pAb3p7,2j2x0005pAb3p7,2j2x0007pAi2x0008*9pEpAp21pAb3p7,2j2x0006pA
b3p7,2j2x0008pAp33*9p17pAp21pAp27pAp3ApAp43*9p17pAp21pAp2DpAp40pA
p33*3p6p24pA*9p6pAp17pAp27pAp2DpAp33*9p6pAp17pAp3ApAp40pAp33*5p6p
Cp17pCpC*5p17pCp3ApCi2x0002*5p17pCp2DpCp7E*5p17pCp27pCpC*5p17pCp2
1pCpC*5pEpCp2DpCp7E*5pEpCp40pCpC*5pEpCp3ApCp7E*5pEpCp21pCpC*5pEpC
p27pCpC*3p17p24pA*3p21p24pA*5p21pCp3ApCpC*5p21pCp2DpCpC*5p21pCp40
pCp7E*5p21pCp27pCp7E*5p17pCp40pCpC*3p2Dp24pA*3p40p24pA*3p3Ap24pA*
3p27p24pA*5p6pCp2DpCpC*5p6pCp3ApCpC*5p6pCp27pCp7E*5p27pCp2DpCpC*5
p27pCp40pCpC*5p6pCp21pCpC*5p3ApCp2DpCpC*5p3ApCp40pCpC*5p2DpCp40pC
pC*5p27pCp3ApCpC*5p6pCp40pCp7E*3pEp24pA*9p6pApEpAp27pAp3ApAp33*9p
6pApEpAp2DpAp40pAp43/3pApC}}}{SCP_R 97 0 2{COM_R 99 0{TEXT 1 27 "
La funci\363n generatriz W es:"}}{OUT_R 102 0 0{TEXT 2 1 "\012"}}
}{SCP_R 104 0 2{INP_R 106 0 "> "{TEXT 0 25 "W:=sum(p[r]*f[r],r=1.
.8);"}}{OUT_R 109 8192 106{DAG :3n3\`W`+D*5b3n3\`p`,2j2x0001p9+9*
3b3n3\`Q`p8j2x0002p9*3b3pF,2p12p12i2x0001*3b3pF,2j2x0003p12p1B*3b
3pF,2j2x0004p12p9p9p9*5b3p6p18p9+5*5pEp9p16p9p12*5p1Ep9p26p9i2x00
02p9p9*5b3p6p20p9+5*5pEp9p1Ep9p12*5p16p9p26p9p12p9p9*5b3p6,2j2x00
05p9+9*3b3pFp5Bp12p9*3b3pF,2j2x0006p12p1B*3b3pF,2j2x0007p12p1B*3b
3pF,2j2x0008p12p9p9p9*5b3p6p69p9+5*5p61p9p67p9p12*5p6Fp9p77p9p40p
9p9*5b3p6p71p9+5*5p61p9p6Fp9p12*5p67p9p77p9p12p9p9}}}{SCP_R 111 0
 4{INP_R 113 0 "> "{TEXT 0 22 "Coordq:=\{q[r]$r=1..4\};"}}{INP_R 
116 0 "> "{TEXT 0 22 "CoordQ:=\{Q[r]$r=1..4\};"}}{OUT_R 119 8192 
116{DAG :3n4\`Coordq`%2,5b3n3\`q`,2j2x0002b3p7,2j2x0003b3p7,2j2x0
001b3p7,2j2x0004}}{OUT_R 122 8192 116{DAG :3n4\`CoordQ`%2,5b3n3\`
Q`,2j2x0001b3p7,2j2x0002b3p7,2j2x0003b3p7,2j2x0004}}}{SCP_R 124 0
 2{INP_R 126 0 "> "{TEXT 0 29 "ecs_qQ:=\{(q[r]=f[r])$r=1..4\};"}}
{OUT_R 129 8192 126{DAG :3n4\`ecs_qQ`%2,5=3b3n3\`q`,2j2x0001+9*3b
3n3\`Q`pAj2x0002pB*3b3p10,2p13p13i2x0001*3b3p10,2j2x0003p13p1C*3b
3p10,2j2x0004p13pB=3b3p8p19+5*5pFpBp17pBp13*5p1FpBp27pBi2x0002=3b
3p8p21+5*5pFpBp1FpBp13*5p17pBp27pBp13=3b3p8p29j2x0000}}}{SCP_R 
131 0 2{INP_R 133 0 "> "{TEXT 0 21 "solve(ecs_qQ,CoordQ);"}}
{OUT_R 136 0 133{TEXT 2 1 "\012"}}}{SCP_R 138 0 2{INP_R 140 0 "> 
"{TEXT 0 58 "expand(R1^2+R2^2-2*f[1]*f[5]-2*f[2]*f[6]-2*f[3]*f[7]
)-R^2;"}}{OUT_R 143 8192 140{DAG j2x0000}}}{SCP_R 145 0 2{INP_R 
147 0 "> "{TEXT 0 48 "'R1^2+R2^2-2*f[1]*f[5]-2*f[2]*f[6]-2*f[3]*f
[7]';"}}{OUT_R 150 8192 147{DAG +B*3n3\`R1`j2x0002j2x0001*3n3\`R2
`p4p6*5b3n3\`f`,2p6p6b3pF,2j2x0005p6i2x0002*5b3pF,2p4p6b3pF,2j2x0
006p6p1A*5b3pF,2j2x0003p6b3pF,2j2x0007p6p1A}}}{SCP_R 152 0 8
{INP_R 154 0 "> "{TEXT 0 16 "A1:=matrix(4,4):"}}{INP_R 157 0 "> "
{TEXT 0 20 "for i from 1 to 4 do"}}{INP_R 160 0 "> "{TEXT 0 20 "f
or j from 1 to 4 do"}}{INP_R 163 0 "> "{TEXT 0 25 "A1[i,j]:=diff(
f[j],Q[i]);"}}{INP_R 166 0 "> "{TEXT 0 3 "od;"}}{INP_R 169 0 "> "
{TEXT 0 3 "od;"}}{INP_R 172 0 "> "{TEXT 0 15 "'A1'=evalm(A1);"}}
{OUT_R 175 8192 172{DAG =3n3\`A1`(3n4\`MATRIX`,2[2,5[2,5+3b3n3\`Q
`,2j2x0001j2x0002+3b3pE,2p13p13+3b3pE,2j2x0003p13j2x0000[2,5+3p16
i2x0002pC+3b3pE,2j2x0004p13p22[2,5+3p1Cp28+3p2Cp28pCp22[2,5p2Bp35
p15p22}}}{SCP_R 177 0 8{INP_R 179 0 "> "{TEXT 0 16 "A2:=matrix(4,
4):"}}{INP_R 182 0 "> "{TEXT 0 20 "for i from 1 to 4 do"}}{INP_R 
185 0 "> "{TEXT 0 20 "for j from 1 to 4 do"}}{INP_R 188 0 "> "
{TEXT 0 29 "A2[i,j]:=diff(f[j+4],Q[i+4]);"}}{INP_R 191 0 "> "
{TEXT 0 3 "od;"}}{INP_R 194 0 "> "{TEXT 0 3 "od;"}}{INP_R 197 0 "
> "{TEXT 0 15 "'A2'=evalm(A2);"}}{OUT_R 200 8192 197{DAG =3n3\`A2
`(3n4\`MATRIX`,2[2,5[2,5+3b3n3\`Q`,2j2x0005j2x0002+3b3pE,2j2x0006
p13+3b3pE,2j2x0007p13j2x0000[2,5+3p16i2x0002pC+3b3pE,2j2x0008p13p
23[2,5+3p1Dp29+3p2Dp29pCp23[2,5p2Cp36p15p23}}}{SCP_R 202 0 2
{INP_R 204 0 "> "{TEXT 0 38 "q1:=matrix(4,1,[q[1],q[2],q[3],q[4]]
);"}}{OUT_R 207 8192 204{DAG :3n3\`q1`(3n4\`MATRIX`,2[2,5[2,2b3n3
\`q`,2j2x0001[2,2b3pD,2j2x0002[2,2b3pD,2j2x0003[2,2b3pD,2j2x0004}
}}{SCP_R 209 0 1{COM_R 211 0{TEXT 1 66 "Verifiquemos la soluci\36
3n de las ecuaciones de transformaci\363n q->Q:"}}}{SCP_R 213 0 8
{INP_R 215 0 "> "{TEXT 0 38 "Q[1]:=sqrt((1/2)*(sqrt(mc(q1))+q[1])
);"}}{INP_R 218 0 "> "{TEXT 0 20 "Q[2]:=q[2]/(2*Q[1]);"}}{INP_R 
221 0 "> "{TEXT 0 20 "Q[3]:=q[3]/(2*Q[1]);"}}{INP_R 224 0 "> "
{TEXT 0 8 "Q[4]:=0;"}}{OUT_R 227 8192 224{DAG :3b3n3\`Q`,2j2x0001
*3+5*3+9*3b3n3\`q`p4j2x0002p5*3b3pD,2p10p10p5*3b3pD,2j2x0003p10p5
*3b3pD,2j2x0004p10p5/3p5p10p2ApCp2Ap2A}}{OUT_R 230 8192 224{DAG :
3b3n3\`Q`,2j2x0002+3*5b3n3\`q`p4j2x0001+5*3+9*3b3pA,2pDp5pD*3p9p5
pD*3b3pA,2j2x0003p5pD*3b3pA,2j2x0004p5pD/3pDp5p2Dp13p2D/3i2x0001p
5p2D}}{OUT_R 233 8192 224{DAG :3b3n3\`Q`,2j2x0003+3*5b3n3\`q`p4j2
x0001+5*3+9*3b3pA,2pDj2x0002pD*3b3pA,2p17p17pD*3p9p17pD*3b3pA,2j2
x0004p17pD/3pDp17p2Dp13p2D/3i2x0001p17p2D}}{OUT_R 236 8192 224
{DAG :3b3n3\`Q`,2j2x0004j2x0000}}}{SCP_R 238 0 5{INP_R 240 0 "> "
{TEXT 0 39 "for i from 1 to 4 do simplify(f[i]) od;"}}{OUT_R 243 
8192 240{DAG +3*5+7*3b3n3\`q`,2j2x0001j2x0002pA*3b3p5,2j2x0004pAp
8*5+9p3p8*3b3p5,2pApAp8*3b3p5,2j2x0003pAp8pDp8/3p8pAp4p8pAp8+5*3p
16p2Ap8p4p8i2x0001p2A}}{OUT_R 246 8192 240{DAG b3n3\`q`,2j2x0002}
}{OUT_R 249 8192 240{DAG b3n3\`q`,2j2x0003}}{OUT_R 252 8192 240
{DAG j2x0000}}}{SCP_R 254 0 2{COM_R 256 0{TEXT 1 92 "Usemos el he
cho de que q[4] es IDENTICAMENTE cero para poder verificar la ide
ntidad en q[1]:"}}{OUT_R 259 0 0{TEXT 2 1 "\012"}}}{SCP_R 261 0 4
{INP_R 263 0 "> "{TEXT 0 8 "q[4]:=0;"}}{INP_R 266 0 "> "{TEXT 0 
15 "simplify(f[1]);"}}{OUT_R 269 8192 266{DAG :3b3n3\`q`,2j2x0004
j2x0000}}{OUT_R 272 8192 266{DAG b3n3\`q`,2j2x0001}}}{SCP_R 274 0
 2{INP_R 276 0 "> "{TEXT 0 45 "p1:=matrix(4,1,[qr[9],qr[10],qr[11
],qr[12]]);"}}{OUT_R 279 8192 276{DAG :3n3\`p1`(3n4\`MATRIX`,2[2,
5[2,2b3n3\`qr`,2j2x0009[2,2b3pD,2j2x0010[2,2b3pD,2j2x0011[2,2b3pD
,2j2x0012}}}{SCP_R 281 0 2{INP_R 283 0 "> "{TEXT 0 46 "p2:=matrix
(4,1,[qr[13],qr[14],qr[15],qr[16]]);"}}{OUT_R 286 8192 283{DAG :3
n3\`p2`(3n4\`MATRIX`,2[2,5[2,2b3n3\`qr`,2j2x0013[2,2b3pD,2j2x0014
[2,2b3pD,2j2x0015[2,2b3pD,2j2x0016}}}{SCP_R 288 0 2{INP_R 290 0 "
> "{TEXT 0 11 "readlib(C):"}}{OUT_R 293 0 290{TEXT 2 1 "\012"}}}
{SCP_R 295 0 2{INP_R 297 0 "> "{TEXT 0 23 "Q9med:=evalm(A1/2&*p1)
;"}}{OUT_R 300 8192 297{DAG :3n4\`Q9med`(3n4\`MATRIX`,2[2,5[2,2+7
*5b3n3\`Q`,2j2x0001p13b3n3\`qr`,2j2x0009p13p13*5b3p10,2j2x0002p13
b3p17,2j2x0010p13p13*5b3p10,2j2x0003p13b3p17,2j2x0011p13p13[2,2+7
*5p1Fp13p16p13i2x0001*5pFp13p25p13p13*5b3p10,2j2x0004p13p33p13p13
[2,2+7*5p2Dp13p16p13p42*5p4Bp13p25p13p42*5pFp13p33p13p13[2,2+7*5p
4Bp13p16p13p13*5p2Dp13p25p13p42*5p1Fp13p33p13p13}}}{SCP_R 302 0 4
{INP_R 304 0 "> "{TEXT 0 20 "for i from 1 to 4 do"}}{INP_R 307 0 
"> "{TEXT 0 14 "C(Q9med[i,1]);"}}{INP_R 310 0 "> "{TEXT 0 3 "od;"
}}{OUT_R 313 0 310{TEXT 2 191 "      t0 = Q[1]*qr[9]+Q[2]*qr[10]+
Q[3]*qr[11];\012      t0 = -Q[2]*qr[9]+Q[1]*qr[10]+Q[4]*qr[11];\0
12      t0 = -Q[3]*qr[9]-Q[4]*qr[10]+Q[1]*qr[11];\012      t0 = Q
[4]*qr[9]-Q[3]*qr[10]+Q[2]*qr[11];\012\012"}}}{SCP_R 315 0 2
{INP_R 317 0 "> "{TEXT 0 24 "Q13med:=evalm(A2/2&*p2);"}}{OUT_R 
320 8192 317{DAG :3n4\`Q13med`(3n4\`MATRIX`,2[2,5[2,2+7*5b3n3\`Q`
,2j2x0005j2x0001b3n3\`qr`,2j2x0013p15p15*5b3p10,2j2x0006p15b3p18,
2j2x0014p15p15*5b3p10,2j2x0007p15b3p18,2j2x0015p15p15[2,2+7*5p20p
15p17p15i2x0001*5pFp15p26p15p15*5b3p10,2j2x0008p15p34p15p15[2,2+7
*5p2Ep15p17p15p43*5p4Cp15p26p15p43*5pFp15p34p15p15[2,2+7*5p4Cp15p
17p15p15*5p2Ep15p26p15p43*5p20p15p34p15p15}}}{SCP_R 322 0 4{INP_R
 324 0 "> "{TEXT 0 20 "for i from 1 to 4 do"}}{INP_R 327 0 "> "
{TEXT 0 15 "C(Q13med[i,1]);"}}{INP_R 330 0 "> "{TEXT 0 3 "od;"}}
{OUT_R 333 0 330{TEXT 2 195 "      t0 = Q[5]*qr[13]+Q[6]*qr[14]+Q
[7]*qr[15];\012      t0 = -Q[6]*qr[13]+Q[5]*qr[14]+Q[8]*qr[15];\0
12      t0 = -Q[7]*qr[13]-Q[8]*qr[14]+Q[5]*qr[15];\012      t0 = 
Q[8]*qr[13]-Q[7]*qr[14]+Q[6]*qr[15];\012\012"}}}{SCP_R 335 0 4
{INP_R 337 0 "> "{TEXT 0 20 "for i from 1 to 8 do"}}{INP_R 340 0 
"> "{TEXT 0 8 "C(f[i]);"}}{INP_R 343 0 "> "{TEXT 0 3 "od;"}}
{OUT_R 346 0 343{TEXT 2 329 "      t0 = pow(Q[1],2.0)-pow(Q[2],2.
0)-pow(Q[3],2.0)+pow(Q[4],2.0);\012      t0 = 2.0*Q[1]*Q[2]-2.0*Q
[3]*Q[4];\012      t0 = 2.0*Q[1]*Q[3]+2.0*Q[2]*Q[4];\012      t0 
= 0.0;\012      t0 = pow(Q[5],2.0)-pow(Q[6],2.0)-pow(Q[7],2.0)+po
w(Q[8],2.0);\012      t0 = 2.0*Q[5]*Q[6]-2.0*Q[7]*Q[8];\012      
t0 = 2.0*Q[5]*Q[7]+2.0*Q[6]*Q[8];\012      t0 = 0.0;\012\012"}}}
{SCP_R 348 0 2{INP_R 350 0 "> "{TEXT 0 41 "P1:=matrix(4,1,[Q[9],Q
[10],Q[11],Q[12]]):"}}{INP_R 353 0 "> "{TEXT 0 42 "P2:=matrix(4,1
,[Q[13],Q[14],Q[15],Q[16]]):"}}}{SCP_R 355 0 2{INP_R 357 0 "> "
{TEXT 0 30 "qr9p2:=evalm((1/2)*T(A1)&*P1);"}}{OUT_R 360 8192 357
{DAG :3n4\`qr9p2`(3n4\`MATRIX`,2[2,5[2,2+9*5b3n3\`Q`,2j2x0001p13b
3p10,2j2x0009p13p13*5b3p10,2j2x0002p13b3p10,2j2x0010p13i2x0001*5b
3p10,2j2x0003p13b3p10,2j2x0011p13p2A*5b3p10,2j2x0004p13b3p10,2j2x
0012p13p13[2,2+9*5p1Ep13p16p13p13*5pFp13p24p13p13*5p3Bp13p33p13p2
A*5p2Dp13p41p13p2A[2,2+9*5p2Dp13p16p13p13*5p3Bp13p24p13p13*5pFp13
p33p13p13*5p1Ep13p41p13p13[2,2j2x0000}}}{SCP_R 362 0 4{INP_R 364 
0 "> "{TEXT 0 20 "for i from 1 to 4 do"}}{INP_R 367 0 "> "{TEXT 0
 14 "C(qr9p2[i,1]);"}}{INP_R 370 0 "> "{TEXT 0 3 "od;"}}{OUT_R 
373 0 370{TEXT 2 182 "      t0 = Q[1]*Q[9]-Q[2]*Q[10]-Q[3]*Q[11]+
Q[4]*Q[12];\012      t0 = Q[2]*Q[9]+Q[1]*Q[10]-Q[4]*Q[11]-Q[3]*Q[
12];\012      t0 = Q[3]*Q[9]+Q[4]*Q[10]+Q[1]*Q[11]+Q[2]*Q[12];\01
2      t0 = 0.0;\012\012"}}}{SCP_R 375 0 2{INP_R 377 0 "> "{TEXT 
0 31 "qr13p2:=evalm((1/2)*T(A2)&*P2);"}}{OUT_R 380 8192 377{DAG :
3n4\`qr13p2`(3n4\`MATRIX`,2[2,5[2,2+9*5b3n3\`Q`,2j2x0005j2x0001b3
p10,2j2x0013p15p15*5b3p10,2j2x0006p15b3p10,2j2x0014p15i2x0001*5b3
p10,2j2x0007p15b3p10,2j2x0015p15p2B*5b3p10,2j2x0008p15b3p10,2j2x0
016p15p15[2,2+9*5p1Fp15p17p15p15*5pFp15p25p15p15*5p3Cp15p34p15p2B
*5p2Ep15p42p15p2B[2,2+9*5p2Ep15p17p15p15*5p3Cp15p25p15p15*5pFp15p
34p15p15*5p1Fp15p42p15p15[2,2j2x0000}}}{SCP_R 382 0 4{INP_R 384 0
 "> "{TEXT 0 20 "for i from 1 to 4 do"}}{INP_R 387 0 "> "{TEXT 0 
15 "C(qr13p2[i,1]);"}}{INP_R 390 0 "> "{TEXT 0 3 "od;"}}{OUT_R 
393 0 390{TEXT 2 185 "      t0 = Q[5]*Q[13]-Q[6]*Q[14]-Q[7]*Q[15]
+Q[8]*Q[16];\012      t0 = Q[6]*Q[13]+Q[5]*Q[14]-Q[8]*Q[15]-Q[7]*
Q[16];\012      t0 = Q[7]*Q[13]+Q[8]*Q[14]+Q[5]*Q[15]+Q[6]*Q[16];
\012      t0 = 0.0;\012\012"}}}{SCP_R 395 0 1{INP_R 397 0 "> "
{TEXT 0 0 ""}}}}{END}

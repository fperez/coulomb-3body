{VERSION 1 1 "MS-Windows" "1.0"}{GLOBALS 1 2}{FONT 0 "Arial" "Ari
al" "Helvetica-Bold" 7 10 0 "Helvetica-Bold" 10{COLOR "0000ff"}}
{FONT 1 "Times New Roman" "Times New Roman" "Times-Roman" 4 10 64
 "Times-Roman" 10{COLOR "000000"}}{FONT 2 "Courier New" "Courier 
New" "Courier" 4 10 192 "Courier" 10{COLOR "000000"}}{SCP_R 0 0 6
{SCP_R 1 0 2{COM_R 3 0{TEXT 1 99 "C\341lculos de las derivadas de
 la superficie de secci\363n con respecto a las coordenadas regul
arizadas:"}}{COM_R 6 0{TEXT 1 75 "*******************************
********************************************"}}}{SCP_R 8 0 2
{INP_R 10 0 "> "{TEXT 0 59 "S:=2*(q[5]*q[6]*(q[1]^2-q[2]^2)-q[1]*
q[2]*(q[5]^2-q[6]^2));"}}{OUT_R 13 8192 10{DAG :3n3\`S`+5*7b3n3\`
q`,2j2x0005j2x0001b3p6,2j2x0006pB+5*3b3p6,2pBj2x0002pB*3b3p6,2p19
p19i2x0001pBp19*7p15pBp1DpB+5*3p5p19pB*3pDp19p22pBi2x0002}}}
{SCP_R 15 0 11{INP_R 17 0 "> "{TEXT 0 20 "for i from 1 to 8 do"}}
{INP_R 20 0 "> "{TEXT 0 20 "dS[i]:=diff(S,q[i]);"}}{INP_R 23 0 ">
 "{TEXT 0 3 "od;"}}{OUT_R 26 8192 23{DAG :3b3n3\`dS`,2j2x0001+5*7
b3n3\`q`,2j2x0005p5b3pA,2j2x0006p5b3pAp4p5j2x0004*5b3pA,2j2x0002p
5+5*3p9p20p5*3p10p20i2x0001p5i2x0002}}{OUT_R 29 8192 23{DAG :3b3n
3\`dS`,2j2x0002+5*7b3n3\`q`,2j2x0005j2x0001b3pA,2j2x0006pFb3pAp4p
Fi2x0004*5b3pA,2pFpF+5*3p9p5pF*3p11p5i2x0001pFi2x0002}}{OUT_R 32 
8192 23{DAG :3b3n3\`dS`,2j2x0003j2x0000}}{OUT_R 35 8192 23{DAG :3
b3n3\`dS`,2j2x0004j2x0000}}{OUT_R 38 8192 23{DAG :3b3n3\`dS`,2j2x
0005+5*5b3n3\`q`,2j2x0006j2x0001+5*3b3pA,2pFj2x0002pF*3b3pA,2p17p
17i2x0001pFp17*7p13pFp1BpFb3pAp4pFi2x0004}}{OUT_R 41 8192 23{DAG 
:3b3n3\`dS`,2j2x0006+5*5b3n3\`q`,2j2x0005j2x0001+5*3b3pA,2pFj2x00
02pF*3b3pA,2p17p17i2x0001pFp17*7p13pFp1BpFb3pAp4pFj2x0004}}{OUT_R
 44 8192 23{DAG :3b3n3\`dS`,2j2x0007j2x0000}}{OUT_R 47 8192 23
{DAG :3b3n3\`dS`,2j2x0008j2x0000}}}{SCP_R 49 0 2{INP_R 51 0 "> "
{TEXT 0 11 "readlib(C);"}}{OUT_R 54 0 51{TEXT 2 17 "proc() ... en
d\012\012\012"}}}{SCP_R 56 0 4{INP_R 58 0 "> "{TEXT 0 20 "for i f
rom 1 to 8 do"}}{INP_R 61 0 "> "{TEXT 0 19 "C(subs(q=x,dS[i]));"}
}{INP_R 64 0 "> "{TEXT 0 3 "od;"}}{OUT_R 67 0 64{TEXT 2 346 "    
  t0 = 4.0*x[5]*x[6]*x[1]-2.0*x[2]*(pow(x[5],2.0)-pow(x[6],2.0));
\012      t0 = -4.0*x[5]*x[6]*x[2]-2.0*x[1]*(pow(x[5],2.0)-pow(x[
6],2.0));\012      t0 = 0.0;\012      t0 = 0.0;\012      t0 = 2.0
*x[6]*(pow(x[1],2.0)-pow(x[2],2.0))-4.0*x[1]*x[2]*x[5];\012      
t0 = 2.0*x[5]*(pow(x[1],2.0)-pow(x[2],2.0))+4.0*x[1]*x[2]*x[6];\0
12      t0 = 0.0;\012      t0 = 0.0;\012\012"}}}{SCP_R 69 0 1
{INP_R 71 0 "> "{TEXT 0 0 ""}}}}{END}

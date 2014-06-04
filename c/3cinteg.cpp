/* Rutinas de integracion de ecuaciones diferenciales. Tomadas de NUMERICAL
   RECIPES in C:                                                */

// #includes comunes a todos los archivos:
#include <stdlib.h>
#include <math.h>   
#include <conio.h>
#include "3cdefs.h"

// #includes que no son comunes a todos
#include <fstream.h>

// constantes para los algoritmos de NUMERICAL RECIPES:
#define PGROW -0.20
#define PSHRNK -0.25
#define FCOR 0.06666666   // ÷ 1.0/15.0
#define SAFETY 0.9
#define ERRCON 6.0e-4     // ÷ (4/SAFETY)^(1/PGROW)

// Las siguientes vbles estan declaradas en REGMAIN.CPP:
// Nombres DOS de los archivos utilizados:
// Directorio para datos:
extern caminoDOS dirdatos;
/* Coordenadas de posicion de c/part (salida); Energia total del sistema, supuesta
   constante; Coordenadas relativas; Mapas de Poincaré: */
extern caminoDOS arch_pos[NPART+1], arch_etot, arch_rel, arch_poinc[NEQ/2+1];

// variables globales:
extern double m1,m2,m3,mu13,mu23,Mtot,etot;
// maxmem define el numero maximo de valores q'se almacenan en memoria
extern int Z1,Z2,Z3,maxmem;

int ind_coord;  // índice de la coordenada cuyo mapa de Poincaré se calcula.

/* funcion usada en la rutina de integracion, proporciona los lados
   derechos del sistema de ecuaciones diferenciales: */
void derivs_norm(double t,double *Q,double *dQdt)
{ double Q2[NEQ/2+1],Rr1,Rr2,R12,P1c,P2c,*P,dRdQ[NEQ/2+1],gv[4];
  int i;

  P = &Q[8];  // se usa Q para coordenadas y momentos.
#ifdef SIS_PLANO    // simplificaciones en el caso plano:
// ***********************************************
// Codigo generado por Maple V con la funciOn C():
//   P1cuad = P1c
      P1c = P[1]*P[1] + P[2]*P[2];
//   P2cuad = P2c
      P2c = P[5]*P[5] + P[6]*P[6];
// con frecuencia se necesitan los cuadrados de las coordenadas:
     for(i=1;i<=NEQ/2;i++) Q2[i] = Q[i]*Q[i];
//   R1 = Rr1
      Rr1 = Q2[1]+Q2[2];
//   R2 = Rr2
      Rr2 = Q2[5]+Q2[6];
//   R = R12
      R12 = sqrt(Rr1*Rr1-2.0*(Q2[1]-Q2[2])*(Q2[5]-Q2[6])
      -8.0*(Q[1]*Q[2]*Q[5]*Q[6])+Rr2*Rr2);

// ECUACIONES DE MVTO:
// Las derivadas de R con respecto a las Q[j] (gv es un vector auxiliar):
      gv[1] = Q2[1]-Q2[2]-Q2[5]+Q2[6];
      gv[2] = 2.0*(Q[1]*Q[2]-Q[5]*Q[6]);
      dRdQ[1] = 2.0*(Q[1]*gv[1]+Q[2]*gv[2])/R12;
      dRdQ[2] = 2.0*(-Q[2]*gv[1]+Q[1]*gv[2])/R12;
      dRdQ[5] = -2.0*(Q[5]*gv[1]+Q[6]*gv[2])/R12;
      dRdQ[6] = -2.0*(-Q[6]*gv[1]+Q[5]*gv[2])/R12;

// Y ahora las ecuaciones propiamente dichas:

      dQdt[3]=dQdt[4]=dQdt[7]=dQdt[8]=dQdt[11]=dQdt[12]=dQdt[15]=dQdt[16]=0;

      dQdt[1] = Rr2/mu13*P[1]/4-(-P[5]*Q[1]*Q[5]-P[5]*Q[2]*Q[6]+
P[6]*Q[1]*Q[6]-P[6]*Q[2]*Q[5])/m3/4;

      dQdt[2] = Rr2/mu13*P[2]/4+(-P[5]*Q[2]*Q[5]+P[5]*Q[1]*Q[6]+
P[6]*Q[2]*Q[6]+P[6]*Q[1]*Q[5])/m3/4;

      dQdt[5] = Rr1/mu23*P[5]/4+(P[1]*Q[1]*Q[5]+P[1]*Q[2]*Q[6]-
P[2]*Q[2]*Q[5]+P[2]*Q[1]*Q[6])/m3/4;

      dQdt[6] = Rr1/mu23*P[6]/4+(-P[1]*Q[1]*Q[6]+P[1]*Q[2]*Q[5]+
P[2]*Q[2]*Q[6]+P[2]*Q[1]*Q[5])/m3/4;

      dQdt[9] = -Q[1]/mu23*P2c/4-(P[5]*P[1]*Q[5]+P[5]*P[2]*Q[6]-
P[6]*P[1]*Q[6]+P[6]*P[2]*Q[5])/m3/4+2.0*Z2*Z3*Q[1]
-2.0*Q[1]*Rr2*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[1];

      dQdt[10] = -Q[2]/mu23*P2c/4-(P[5]*P[1]*Q[6]-P[5]*P[2]*Q[5]+
P[6]*P[1]*Q[5]+P[6]*P[2]*Q[6])/m3/4+2.0*Z2*Z3*Q[2]
-2.0*Q[2]*Rr2*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[2];

      dQdt[13] = -Q[5]/mu13*P1c/4-(P[5]*P[1]*Q[1]-P[5]*P[2]*Q[2]+
P[6]*P[1]*Q[2]+P[6]*P[2]*Q[1])/m3/4+2.0*Z1*Z3*Q[5]
-2.0*Rr1*Q[5]*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[5];

      dQdt[14] = -Q[6]/mu13*P1c/4-(P[5]*P[1]*Q[2]+P[5]*P[2]*Q[1]-
P[6]*P[1]*Q[1]+P[6]*P[2]*Q[2])/m3/4+2.0*Z1*Z3*Q[6]
-2.0*Rr1*Q[6]*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[6];
// Termina el codigo generado por Maple V.
// ***************************************
#else    // ecuaciones grales del caso tridimensional:
// ***********************************************
// Codigo generado por Maple V con la funciOn C():
//   P1cuad = P1c
      P1c = P[1]*P[1] + P[2]*P[2] + P[3]*P[3] + P[4]*P[4];
//   P2cuad = P2c
      P2c = P[5]*P[5] + P[6]*P[6] + P[7]*P[7] + P[8]*P[8];
// con frecuencia se necesitan los cuadrados de las coordenadas:
     for(i=1;i<=NEQ/2;i++) Q2[i] = Q[i]*Q[i];
//   R1 = Rr1
      Rr1 = Q2[1]+Q2[2]+Q2[3]+Q2[4];
//   R2 = Rr2
      Rr2 = Q2[5]+Q2[6]+Q2[7]+Q2[8];
//   R = R12
      R12 = sqrt(Rr1*Rr1-2.0*(Q2[1]-Q2[2]-Q2[3]+Q2[4])*(Q2[5]-Q2[6]-Q2[7]+Q2[8])
      -8.0*( (Q[1]*Q[2]-Q[3]*Q[4])*(Q[5]*Q[6]-Q[7]*Q[8]) +
	     (Q[1]*Q[3]+Q[2]*Q[4])*(Q[5]*Q[7]+Q[6]*Q[8]) )+Rr2*Rr2);

// ECUACIONES DE MVTO:
// Las derivadas de R con respecto a las Q[j] (gv es un vector auxiliar):
      gv[1] = Q2[1]-Q2[2]-Q2[3]+Q2[4]-Q2[5]+Q2[6]+Q2[7]-Q2[8];
      gv[2] = 2.0*(Q[1]*Q[2]-Q[3]*Q[4]-Q[5]*Q[6]+Q[7]*Q[8]);
      gv[3] = 2.0*(Q[1]*Q[3]+Q[2]*Q[4]-Q[5]*Q[7]-Q[6]*Q[8]);
      dRdQ[1] = 2.0*(Q[1]*gv[1]+Q[2]*gv[2]+Q[3]*gv[3])/R12;
      dRdQ[2] = 2.0*(-Q[2]*gv[1]+Q[1]*gv[2]+Q[4]*gv[3])/R12;
      dRdQ[3] = -2.0*(Q[3]*gv[1]+Q[4]*gv[2]-Q[1]*gv[3])/R12;
      dRdQ[4] = -2.0*(-Q[4]*gv[1]+Q[3]*gv[2]-Q[2]*gv[3])/R12;
      dRdQ[5] = -2.0*(Q[5]*gv[1]+Q[6]*gv[2]+Q[7]*gv[3])/R12;
      dRdQ[6] = -2.0*(-Q[6]*gv[1]+Q[5]*gv[2]+Q[8]*gv[3])/R12;
      dRdQ[7] = 2.0*(Q[7]*gv[1]+Q[8]*gv[2]-Q[5]*gv[3])/R12;
      dRdQ[8] = 2.0*(-Q[8]*gv[1]+Q[7]*gv[2]-Q[6]*gv[3])/R12;

// Y ahora las ecuaciones propiamente dichas:
      dQdt[1] = Rr2/mu13*P[1]/4-(-P[5]*Q[1]*Q[5]-P[5]*Q[2]*Q[6]-P[5]*Q[3]*Q[7]+
P[6]*Q[1]*Q[6]-P[6]*Q[2]*Q[5]-P[6]*Q[3]*Q[8]+P[7]*Q[1]*Q[7]+P[7]*Q[2]*Q[8]-P[7]
*Q[3]*Q[5]-P[8]*Q[1]*Q[8]+P[8]*Q[2]*Q[7]-P[8]*Q[3]*Q[6])/m3/4;

      dQdt[2] = Rr2/mu13*P[2]/4+(-P[5]*Q[2]*Q[5]+P[5]*Q[1]*Q[6]+P[5]*Q[4]*Q[7]+
P[6]*Q[2]*Q[6]+P[6]*Q[1]*Q[5]+P[6]*Q[4]*Q[8]+P[7]*Q[2]*Q[7]-P[7]*Q[1]*Q[8]+P[7]
*Q[4]*Q[5]-P[8]*Q[2]*Q[8]-P[8]*Q[1]*Q[7]+P[8]*Q[4]*Q[6])/m3/4;

      dQdt[3] = Rr2/mu13*P[3]/4+(P[7]*Q[4]*Q[8]-P[5]*Q[3]*Q[5]-P[5]*Q[4]*Q[6]+P
[5]*Q[1]*Q[7]+P[6]*Q[3]*Q[6]-P[6]*Q[4]*Q[5]+P[6]*Q[1]*Q[8]+P[7]*Q[3]*Q[7]+P[7]*
Q[1]*Q[5]-P[8]*Q[3]*Q[8]+P[8]*Q[4]*Q[7]+P[8]*Q[1]*Q[6])/m3/4;

      dQdt[4] = Rr2/mu13*P[4]/4+(P[5]*Q[4]*Q[5]-P[5]*Q[3]*Q[6]+P[5]*Q[2]*Q[7]-P
[6]*Q[4]*Q[6]-P[6]*Q[3]*Q[5]+P[6]*Q[2]*Q[8]-P[7]*Q[4]*Q[7]+P[7]*Q[3]*Q[8]+P[7]*
Q[2]*Q[5]+P[8]*Q[4]*Q[8]+P[8]*Q[3]*Q[7]+P[8]*Q[2]*Q[6])/m3/4;

      dQdt[5] = Rr1/mu23*P[5]/4+(P[1]*Q[1]*Q[5]+P[1]*Q[2]*Q[6]+P[1]*Q[3]*Q[7]-P
[2]*Q[2]*Q[5]+P[2]*Q[1]*Q[6]+P[2]*Q[4]*Q[7]-P[3]*Q[3]*Q[5]-P[3]*Q[4]*Q[6]+P[3]*
Q[1]*Q[7]+P[4]*Q[4]*Q[5]-P[4]*Q[3]*Q[6]+P[4]*Q[2]*Q[7])/m3/4;

      dQdt[6] = Rr1/mu23*P[6]/4+(-P[1]*Q[1]*Q[6]+P[1]*Q[2]*Q[5]+P[1]*Q[3]*Q[8]+
P[2]*Q[2]*Q[6]+P[2]*Q[1]*Q[5]+P[2]*Q[4]*Q[8]+P[3]*Q[3]*Q[6]-P[3]*Q[4]*Q[5]+P[3]
*Q[1]*Q[8]-P[4]*Q[4]*Q[6]-P[4]*Q[3]*Q[5]+P[4]*Q[2]*Q[8])/m3/4;

      dQdt[7] = Rr1/mu23*P[7]/4-(P[1]*Q[1]*Q[7]+P[1]*Q[2]*Q[8]-P[1]*Q[3]*Q[5]-P
[2]*Q[2]*Q[7]+P[2]*Q[1]*Q[8]-P[2]*Q[4]*Q[5]-P[3]*Q[3]*Q[7]-P[3]*Q[4]*Q[8]-P[3]*
Q[1]*Q[5]+P[4]*Q[4]*Q[7]-P[4]*Q[3]*Q[8]-P[4]*Q[2]*Q[5])/m3/4;

      dQdt[8] = Rr1/mu23*P[8]/4-(-P[1]*Q[1]*Q[8]+P[1]*Q[2]*Q[7]-P[1]*Q[3]*Q[6]+
P[2]*Q[2]*Q[8]+P[2]*Q[1]*Q[7]-P[2]*Q[4]*Q[6]+P[3]*Q[3]*Q[8]-P[3]*Q[4]*Q[7]-P[3]
*Q[1]*Q[6]-P[4]*Q[4]*Q[8]-P[4]*Q[3]*Q[7]-P[4]*Q[2]*Q[6])/m3/4;

      dQdt[9] = -Q[1]/mu23*P2c/4-(P[5]*P[1]*Q[5]+P[5]*P[2]*Q[6]+P[5]*P[3]*Q[7]-
P[6]*P[1]*Q[6]+P[6]*P[2]*Q[5]+P[6]*P[3]*Q[8]-P[7]*P[1]*Q[7]-P[7]*P[2]*Q[8]+P[7]
*P[3]*Q[5]+P[8]*P[1]*Q[8]-P[8]*P[2]*Q[7]+P[8]*P[3]*Q[6])/m3/4+2.0*Z2*Z3*Q[1]
-2.0*Q[1]*Rr2*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[1];

      dQdt[10] = -Q[2]/mu23*P2c/4-(P[5]*P[1]*Q[6]-P[5]*P[2]*Q[5]+P[5]*P[4]*Q[7]
+P[6]*P[1]*Q[5]+P[6]*P[2]*Q[6]+P[6]*P[4]*Q[8]-P[7]*P[1]*Q[8]+P[7]*P[2]*Q[7]+P
[7]*P[4]*Q[5]-P[8]*P[1]*Q[7]-P[8]*P[2]*Q[8]+P[8]*P[4]*Q[6])/m3/4+2.0*Z2*Z3*Q[2]
-2.0*Q[2]*Rr2*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[2];

      dQdt[11] = -Q[3]/mu23*P2c/4-(P[5]*P[1]*Q[7]-P[5]*P[3]*Q[5]-P[5]*P[4]*Q[6]
+P[6]*P[1]*Q[8]+P[6]*P[3]*Q[6]-P[6]*P[4]*Q[5]+P[7]*P[1]*Q[5]+P[7]*P[3]*Q[7]+P
[7]*P[4]*Q[8]+P[8]*P[1]*Q[6]-P[8]*P[3]*Q[8]+P[8]*P[4]*Q[7])/m3/4+2.0*Z2*Z3*Q[3]
-2.0*Q[3]*Rr2*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[3];

      dQdt[12] = -Q[4]/mu23*P2c/4-(P[7]*P[3]*Q[8]+P[5]*P[2]*Q[7]-P[5]*P[3]*Q[6]
+P[5]*P[4]*Q[5]+P[6]*P[2]*Q[8]-P[6]*P[3]*Q[5]-P[6]*P[4]*Q[6]+P[7]*P[2]*Q[5]-P
[7]*P[4]*Q[7]+P[8]*P[2]*Q[6]+P[8]*P[3]*Q[7]+P[8]*P[4]*Q[8])/m3/4+2.0*Z2*Z3*Q[4]
-2.0*Q[4]*Rr2*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[4];

      dQdt[13] = -Q[5]/mu13*P1c/4-(P[5]*P[1]*Q[1]-P[5]*P[2]*Q[2]-P[5]*P[3]*Q[3]
+P[5]*P[4]*Q[4]+P[6]*P[1]*Q[2]+P[6]*P[2]*Q[1]-P[6]*P[3]*Q[4]-P[6]*P[4]*Q[3]+P
[7]*P[1]*Q[3]+P[7]*P[2]*Q[4]+P[7]*P[3]*Q[1]+P[7]*P[4]*Q[2])/m3/4+2.0*Z1*Z3*Q[5]
-2.0*Rr1*Q[5]*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[5];

      dQdt[14] = -Q[6]/mu13*P1c/4-(P[5]*P[1]*Q[2]+P[5]*P[2]*Q[1]-P[5]*P[3]*Q[4]
-P[5]*P[4]*Q[3]-P[6]*P[1]*Q[1]+P[6]*P[2]*Q[2]+P[6]*P[3]*Q[3]-P[6]*P[4]*Q[4]+P
[8]*P[1]*Q[3]+P[8]*P[2]*Q[4]+P[8]*P[3]*Q[1]+P[8]*P[4]*Q[2])/m3/4+2.0*Z1*Z3*Q[6]
-2.0*Rr1*Q[6]*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[6];

      dQdt[15] = -Q[7]/mu13*P1c/4-(P[5]*P[1]*Q[3]+P[5]*P[2]*Q[4]+P[5]*P[3]*Q[1]
+P[5]*P[4]*Q[2]-P[7]*P[1]*Q[1]+P[7]*P[2]*Q[2]+P[7]*P[3]*Q[3]-P[7]*P[4]*Q[4]-P
[8]*P[1]*Q[2]-P[8]*P[2]*Q[1]+P[8]*P[3]*Q[4]+P[8]*P[4]*Q[3])/m3/4+2.0*Z1*Z3*Q[7]
-2.0*Rr1*Q[7]*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[7];

      dQdt[16] = -Q[8]/mu13*P1c/4-(P[7]*P[3]*Q[4]+P[6]*P[1]*Q[3]+P[6]*P[2]*Q[4]
+P[6]*P[3]*Q[1]+P[6]*P[4]*Q[2]-P[7]*P[1]*Q[2]-P[7]*P[2]*Q[1]+P[7]*P[4]*Q[3]+P
[8]*P[1]*Q[1]-P[8]*P[2]*Q[2]-P[8]*P[3]*Q[3]+P[8]*P[4]*Q[4])/m3/4+2.0*Z1*Z3*Q[8]
-2.0*Rr1*Q[8]*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[8];
// Termina el codigo generado por Maple V.
// ***************************************
#endif
}   // termina derivs_norm()


/* lados derechos del sistema de ecuaciones diferenciales, ajustados para cálculos
   de secciones de Poincaré:            */
void derivs_poinc(double t,double *Q,double *dQdt)
{ double *f,dSq1,dSq2,dSq3,dSq4,fmas;
  double Q2[NEQ/2+1],Rr1,Rr2,R12,P1c,P2c,*P,dRdQ[NEQ/2+1],gv[4];
  int i;

  P = &Q[8];  // se usa Q para coordenadas y momentos.
  f = dvector(1,NEQ);
// Se calculan algunos parámetros para optimizar las ecnes:
      P1c = P[1]*P[1] + P[2]*P[2];
      P2c = P[5]*P[5] + P[6]*P[6];
// con frecuencia se necesitan los cuadrados de las coordenadas:
     for(i=1;i<=NEQ/2;i++) Q2[i] = Q[i]*Q[i];
      Rr1 = Q2[1]+Q2[2];
      Rr2 = Q2[5]+Q2[6];
      R12 = sqrt(Rr1*Rr1-2.0*(Q2[1]-Q2[2])*(Q2[5]-Q2[6])-8.0*(Q[1]*Q[2]*Q[5]*Q[6])+Rr2*Rr2);

// Las derivadas de R con respecto a las Q[j] (gv es un vector auxiliar):
      gv[1] = Q2[1]-Q2[2]-Q2[5]+Q2[6];
      gv[2] = 2.0*(Q[1]*Q[2]-Q[5]*Q[6]);
      dRdQ[1] = 2.0*(Q[1]*gv[1]+Q[2]*gv[2])/R12;
      dRdQ[2] = 2.0*(-Q[2]*gv[1]+Q[1]*gv[2])/R12;
      dRdQ[5] = -2.0*(Q[5]*gv[1]+Q[6]*gv[2])/R12;
      dRdQ[6] = -2.0*(-Q[6]*gv[1]+Q[5]*gv[2])/R12;

/* Y ahora las ecuaciones propiamente dichas, con las funciones que definen
   el lado derecho del sistema (caso PLANO)  */

      f[3]=f[4]=f[7]=f[8]=f[11]=f[12]=f[15]=f[16]=0;

      f[1] = Rr2/mu13*P[1]/4-(-P[5]*Q[1]*Q[5]-P[5]*Q[2]*Q[6]+
P[6]*Q[1]*Q[6]-P[6]*Q[2]*Q[5])/m3/4;

      f[2] = Rr2/mu13*P[2]/4+(-P[5]*Q[2]*Q[5]+P[5]*Q[1]*Q[6]+
P[6]*Q[2]*Q[6]+P[6]*Q[1]*Q[5])/m3/4;

      f[5] = Rr1/mu23*P[5]/4+(P[1]*Q[1]*Q[5]+P[1]*Q[2]*Q[6]-
P[2]*Q[2]*Q[5]+P[2]*Q[1]*Q[6])/m3/4;

      f[6] = Rr1/mu23*P[6]/4+(-P[1]*Q[1]*Q[6]+P[1]*Q[2]*Q[5]+
P[2]*Q[2]*Q[6]+P[2]*Q[1]*Q[5])/m3/4;

      f[9] = -Q[1]/mu23*P2c/4-(P[5]*P[1]*Q[5]+P[5]*P[2]*Q[6]-
P[6]*P[1]*Q[6]+P[6]*P[2]*Q[5])/m3/4+2.0*Z2*Z3*Q[1]
-2.0*Q[1]*Rr2*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[1];

      f[10] = -Q[2]/mu23*P2c/4-(P[5]*P[1]*Q[6]-P[5]*P[2]*Q[5]+
P[6]*P[1]*Q[5]+P[6]*P[2]*Q[6])/m3/4+2.0*Z2*Z3*Q[2]
-2.0*Q[2]*Rr2*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[2];

      f[13] = -Q[5]/mu13*P1c/4-(P[5]*P[1]*Q[1]-P[5]*P[2]*Q[2]+
P[6]*P[1]*Q[2]+P[6]*P[2]*Q[1])/m3/4+2.0*Z1*Z3*Q[5]
-2.0*Rr1*Q[5]*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[5];

      f[14] = -Q[6]/mu13*P1c/4-(P[5]*P[1]*Q[2]+P[5]*P[2]*Q[1]-
P[6]*P[1]*Q[1]+P[6]*P[2]*Q[2])/m3/4+2.0*Z1*Z3*Q[6]
-2.0*Rr1*Q[6]*(Z1*Z2/R12-etot)+Rr1*Rr2*Z1*Z2/(R12*R12)*dRdQ[6];

  // Construccion de la funcion f[N+1] con las derivadas de la sup de seccion
  dSq1 = Q[2];
  dSq2 = Q[1];
  dSq3 = -Q[4];
  dSq4 = -Q[3];
  fmas = f[1]*dSq1+f[2]*dSq2+f[3]*dSq3+f[4]*dSq4;
  // asignación del lado derecho del sistema ajustando para el cálculo de sec. de Poincaré:
  for (i=1;i<=NEQ;i++) dQdt[i] = f[i]/fmas;
  free_dvector(f,1,NEQ);
}   // termina derzivs_poinc()


// Metodo de Runge-Kutta de 4o orden:
void rk4(double y[],double dydx[],int n,double x,double h,double yout[],
	 void (*derivs)(double,double *,double *) )
{ int i;
  double xh,hh,h6,dym[NEQ1],dyt[NEQ1],yt[NEQ1];

  hh=h*0.5;
  h6=h/6.0;
  xh=x+hh;
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dydx[i];
  (*derivs)(xh,yt,dyt);
  for (i=1;i<=n;i++) yt[i]=y[i]+hh*dyt[i];
  (*derivs)(xh,yt,dym);
  for (i=1;i<=n;i++)
  { yt[i]=y[i]+h*dym[i];
    dym[i] += dyt[i];
  }
  (*derivs)(x+h,yt,dyt);
  for (i=1;i<=n;i++)
    yout[i]=y[i]+h6*(dydx[i]+dyt[i]+2.0*dym[i]);
} // termina rk4()


// rutina de Runge-Kutta con control de error:
void rkqc(double y[],double dydx[],int n,double *x,double htry,
	  double eps,double yscal[],double *hdid,double *hnext,
	  void (*derivs)(double,double *,double *) )

{ int     i;
  double  xsav,hh,h,temp,errmax;
  double  dysav[NEQ1],ysav[NEQ1],ytemp[NEQ1];

  xsav=(*x);
  for (i=1;i<=n;i++)
  { ysav[i]=y[i];
    dysav[i]=dydx[i];
  }
  h=htry;
  for (;;)
  { // se integran dos pasos de ancho h/2:
    hh=0.5*h;
    rk4(ysav,dysav,n,xsav,hh,ytemp,derivs);
    *x=xsav+hh;
    (*derivs)(*x,ytemp,dydx);
    rk4(ytemp,dydx,n,*x,hh,y,derivs);
    // ahora, se hace un paso de ancho h:
    *x=xsav+h;
    if (*x == xsav) nrerror("Step size too small in routine RKQC");
    rk4(ysav,dysav,n,xsav,h,ytemp,derivs);

    // se hace el analisis de error:
    errmax=0.0;
    for (i=1;i<=n;i++)
    { ytemp[i]=y[i]-ytemp[i];
      temp=fabs(ytemp[i]/yscal[i]);
      if (errmax < temp) errmax=temp;
    }
    errmax /= eps;
    if (errmax <= 1.0)
    { *hdid=h;
      *hnext=(errmax > ERRCON ?
	SAFETY*h*exp(PGROW*log(errmax)) : 4.0*h);
      break;
    }
    h=SAFETY*h*exp(PSHRNK*log(errmax));
  }
  for (i=1;i<=n;i++) y[i] += ytemp[i]*FCOR;
}  // termina rkqc()


/* se declaran los streams para grabacion de datos en disco (se hacen
   globales para que tanto odeint() como grabar() tengan acceso a ellos): */
ofstream a_pos[NPART+1],a_rel,a_etot,a_poinc[NEQ/2+1],a_hyper;

void abrir_streams(int pnuevo)
{ int i,modo,out,app;

  out  = 0x02;        // open for writing
  app  = 0x08;        // append mode: all additions at eof

  if (pnuevo) modo = out;
    else modo = app;
  for (i=1;i<=NPART;i++) a_pos[i].open(arch_pos[i]);
  a_rel.open(arch_rel);
  a_etot.open(arch_etot);
//  for (i=1;i<=NEQ/2;i++) a_poinc[i].open(arch_poinc[i],modo);
  a_hyper.open("hyper.dat",modo);
}  // termina abrir_streams()


// rutina que graba los vectores de almacenamiento temporal a disco:
void grabar(int *cont,int lim,double *tdat, double *edat,double ***posdat,double **reldat)
{ int i,j;

  cout << "Grabando Parámetros\n";
  for (i=1;i<=lim;i++)
  { // Posicion y energia:
    for (j=1;j<=NPART;j++)
    { // posiciones
      a_pos[j]<<posdat[j][i][1]<<"   "<<posdat[j][i][2]<<"   "<<posdat[j][i][3]<< "\n";
    }
    a_rel  << reldat[1][i] << "   "  << reldat[2][i] << "\n";  // coord. relativas
    a_etot << tdat[i]      << "   "  << edat[i] << "\n";  // energía(t) -- supues. constante.
  }
  *cont=1;
}  // termina grabar()


// rutina que graba los vectores de Poincaré a disco:
void grabar_poinc(int *cont,int lim, double **xpdat, double **hyperdat)
{ int i,j;

  cout << "\nGrabando datos de Poincaré\n";
  for (i=1;i<=lim;i++)
  {/* for (j=1;j<=NEQ/2;j++)
      a_poinc[j] << xpdat[i][j] << "   " << xpdat[i][j+NEQ/2] << "\n";*/
    a_hyper << hyperdat[i][1] << "   " << hyperdat[i][2] << "\n";
  }
  *cont=1;
}  // termina grabar_poinc()


// rutina que cierra los streams:
void cerrar_streams(void)
{ int i;

  for (i=1;i<=NPART;i++) a_pos[i].close();
  a_rel.close();
  a_etot.close();
  for (i=1;i<=NEQ/2;i++) a_poinc[i].close();
  a_hyper.close();
}   // termina cerrar_streams()


// función signo de x:
int sign(double x)
{ if (x!=0)
  { if (x>0) return 1;
    else return -1;
  }
  else return 0;
}  // termina sign()


/* Rutina que integra un sistema de ecuaciones diferenciales ordinarias y
   genera una sección de Poincaré.    */

void odeint_poinc(double *qstart,int nvar,double t1,double t2,double dtsav,
		  double eps,double h1,double hmin,int *nok,int *nbad,
		  void (*derivs)(double,double *,double *),
		  void (*integ)(double *,double *,int,double *,double,
				double,double *,double *,double *,
				void (*)(double,double *,double *)),
		  int calc_poinc,int sign_mom,int pnuevo)
{ // declaracion de vbles:
  int     i,ip;
  long    nstp;
  double  tsav,t,tvie,w,hnext,hdid,h,tmax=fabs(t2-t1);
  double  qs,*qscal,*q,*dqdt,*xn,*qp,*xp,*dqpdt;

  // vbles para almacenar valores en memoria antes de grabar a disco:
  double  *tmem,*emem,**relmem,**xpmem,**hyperR,***posmem;
  // kmax: numero maximo de valores intermedios que se almacenan en memoria
  // maxmem se declara global y se inicializa a traves del arch. de config.
  int     kmax=maxmem,k=1,kpoinc=1;
  // control de autoionizacion:
  int   ionizo=0;
  // distancia a partir de la cual se considera que hay ionizacion:
  float umbr_ioniz=20;
  // coordenadas relativas y energia intermedia (para control):
  double r1,r2,eint;
  // calculos de mapas de Poincare:
  double Svie,Snue,mom,hyperRad,*qc1,*qc2,*qc3,*rjac1,*rjac2,mu12,mu12_3,c;

  // COMIENZA EL CODIGO DE odeint_poinc():
  abrir_streams(pnuevo);
  // se reserva memoria para los vectores de calculo:
  qscal = dvector(1,nvar);
  q = dvector(1,nvar);
  xn = dvector(1,DIM_FAS);      // xn = coordenadas no regularizadas.
  dqdt = dvector(1,nvar);
  qp = dvector(1,nvar);
  xp = dvector(1,DIM_FAS);
  // qc1, 2 y 3 son vectores de coordenadas, por comodidad para ciertas operaciones:
  qc1 = xp;
  qc2 = &(xp[3]);
  qc3 = &(xp[6]);
  dqpdt = dvector(1,nvar);
  rjac1 = dvector(1,3);
  rjac2 = dvector(1,3);

  // se reserva memoria para los vectores y matrices de almacenamiento:
  tmem = dvector(1,kmax);
  emem = dvector(1,kmax);
  relmem = dmatrix(1,NPART,1,kmax);

  posmem = (new double**[NPART])-1;
  for (i=1;i<=NPART;i++) posmem[i]=dmatrix(1,kmax,1,DIM_ESP);
  xpmem = dmatrix(1,kmax,1,DIM_FAS);
  hyperR = dmatrix(1,kmax,1,2);

  // inicializacion de variables de la rutina odeint() original:
  t=t1;
  h=(t2 > t1) ? fabs(h1) : -fabs(h1);
  dtsav=fabs(dtsav);
  *nok = (*nbad) = 0;
  for (i=1;i<=nvar;i++) q[i]=qstart[i];  // qstart se dan regularizadas

  // se inicializa la vble. para el control del cálculo de sec. de Poincaré:
//  ind_coord=ind_coord_p;   // ind_coord es global, se requiere así para otras rutinas.
  reg_car(q,xp);
  Svie=Snue=0;
  mu12 = m1*m2/(m1+m2);
  mu12_3 = (m1+m2)*m3/(m1+m2+m3);
  c = sqrt(mu12_3/mu12);
  // garantizamos que el primer paso se almacene:
  tsav=t-dtsav*2.0;
  cout << "\nCALCULANDO:\n";

  // CICLO DE CALCULO:
  for(nstp=1;nstp<=MAXSTP;nstp++)
  { // Modulos de los vectores r1 y r2:
    r1=q[1]*q[1] + q[2]*q[2] + q[3]*q[3] + q[4]*q[4];
    r2=q[5]*q[5] + q[6]*q[6] + q[7]*q[7] + q[8]*q[8];
    // detectar autoionizacion:
    if ( (r1>umbr_ioniz) || (r2>umbr_ioniz) ) ionizo=1;

    // indicar que el computador esta 'vivo':
    if (nstp%100==0) cout << (int)(100*(t/tmax)) << " %\n"/*"."*/;

    // se ALMACENAN EN MEMORIA valores solo cada cierto intervalo (dado por dtsav):
    if (!calc_poinc)  // solo se hace si no se estan calculando mapas de Poincare
      if ( fabs(t-tsav) > dtsav )
      { // Posicion y energia:
	reg_car(q,xn);
	energia(xn,&eint);
	tmem[k]=t;
	emem[k]=eint;
	// se almacenan las posiciones:
	posmem[1][k][1]=xn[1]; posmem[1][k][2]=xn[2]; posmem[1][k][3]=xn[3];
	posmem[2][k][1]=xn[4]; posmem[2][k][2]=xn[5]; posmem[2][k][3]=xn[6];
	posmem[3][k][1]=xn[7]; posmem[3][k][2]=xn[8]; posmem[3][k][3]=xn[9];
	relmem[1][k]=r1;
	relmem[2][k]=r2;
	// Se registra el tiempo en el cual se grabó, y se incrementa el
	// contador de valores almacenados:
	tsav=t;
	k++;
      }
    // se GRABA A DISCO cada que los vectores de almacenamiento se llenen:
    if (k==kmax+1) grabar(&k,kmax,tmem,emem,posmem,relmem);
    if (kpoinc==kmax+1) grabar_poinc(&kpoinc,kmax,xpmem,hyperR);

    (*derivs)(t,q,dqdt);
    for (i=1;i<=nvar;i++)
    { qs = fabs(q[i])+fabs(dqdt[i]*h);
      qscal[i] = (qs<TINY?TINY:qs);
    }

    // al hacer la integración, corregimos la relacion entre t y tau:
    tvie=t;
    // llamar la rutina efectiva de integracion:
    (*integ)(q,dqdt,nvar,&t,h,eps,qscal,&hdid,&hnext,derivs);
    t=tvie+hdid*r1*r2;
    // llevar la cuenta de los intentos exitosos y fallidos:
    if (hdid == h) ++(*nok); else ++(*nbad);

    // cálculos para Poincaré:
    if (calc_poinc)
    { Svie = Snue;
      // Relacion algebraica que define la superficie de seccion:
      reg_car(q,xp);
      Snue = qc3[2]; // Sup: m3y=0.
      mom = xp[17];  // p3y>0
      // evaluar cambio de signo en S:
      if ( (Svie*Snue <= 0)   && (mom*sign_mom > 0) )
      { //calculamos el punto en el mapa de Poincaré con vbles auxiliares:
       for (ip=1;ip<=nvar;ip++) qp[ip]=q[ip];
       w=Snue;   		// parámetro de "tiempo"
       derivs_poinc(t,qp,dqpdt);
       rk4(qp,dqpdt,nvar,w,-Snue,qp,derivs_poinc);
       reg_car(qp,xp);
       for (ip=1;ip<=DIM_FAS;ip++) xpmem[kpoinc][ip]=xp[ip];
       // mapa de Poincare del HiperRadio - Momento canónico asociado.
       // calculamos las coordenadas de Jacobi r1 y r2:
       for (ip=1;ip<=3;ip++)
       { rjac1[ip] = qc2[ip] - qc1[ip];
	 rjac2[ip] = (-c/2)*(qc1[ip] + qc2[ip]) + c*qc3[ip];
       }
       // con esto el hiperradio y su momento canónico son:
       hyperR[kpoinc][1] =
       hyperRad = sqrt(rjac1[1]*rjac1[1] + rjac1[2]*rjac1[2] + rjac1[3]*rjac1[3] +
		       rjac2[1]*rjac2[1] + rjac2[2]*rjac2[2] + rjac2[3]*rjac2[3]);
       hyperR[kpoinc][2] =
	(xp[1]*xp[10]+xp[2]*xp[11]+xp[3]*xp[12]+xp[4]*xp[13]+xp[5]*xp[14]+xp[6]*xp[15]+
	 xp[7]*xp[16]+xp[8]*xp[17]+xp[9]*xp[18])/hyperRad;
       kpoinc++;
      }
    } // terminan los calculos de Ponincaré

    // si la integración terminó, se interrumpió o hubo ionización, salir:
    if ( ((t-t2)*(t2-t1) >= 0.0) || (kbhit()) || ionizo)
    { // si se pulsó una tecla o hubo ioniz., indicamos el tiempo alcanzado:
      if (kbhit() || ionizo)
       { if (ionizo)
	   cout << "\n\nCALCULO INTERRUMPIDO POR IONIZACION EN t = " << t;
	 else
	   { cout << "\n\nCALCULO INTERRUMPIDO POR EL USUARIO EN t = " << t;
	     getch();
	   }
	 cout << ". (Se tenía tmax = " << t2 << ").";
       }
      // si hay valores almacenados, grabarlos:
      cout << "\nIntegración terminada\n";
      if (k>1) grabar(&k,k-1,tmem,emem,posmem,relmem);
      if (kpoinc>1) grabar_poinc(&kpoinc,kpoinc-1,xpmem,hyperR);
      // devolver el vector de vbles dependientes en su estado actual:
      for (i=1;i<=nvar;i++) qstart[i]=q[i];
      // liberar la memoria asignada dinamicamente:
      free_dvector(dqdt,1,nvar);  free_dvector(q,1,nvar);
      free_dvector(xn,1,nvar);    free_dvector(qscal,1,nvar);
      free_dvector(xp,1,nvar);    free_dvector(qp,1,nvar);
      free_dvector(tmem,1,kmax);  free_dvector(emem,1,kmax);
      free_dvector(rjac1,1,3);  free_dvector(rjac2,1,3);
      free_dmatrix(relmem,1,NPART,1,kmax);
      for (i=1;i<=NPART;i++) free_dmatrix(posmem[i],1,kmax,1,DIM_ESP);
      delete (posmem+1);
      free_dmatrix(xpmem,1,kmax,1,nvar);
      free_dmatrix(hyperR,1,kmax,1,2);

      // cerrar los archivos de datos:
      cerrar_streams();
      return; // salida normal de la rutina odeint_poinc()
    }
    if (fabs(hnext) <= hmin) nrerror("Paso demasiado pequeño en ODEINT_POINC");
    h=hnext;
  } // termina el ciclo de calculo.

  /* Si llegamos aqui, es porque la rutina no pudo concluir la integracion.
     Cerramos los archivos de datos y salimos, con mensaje de error.
     si hay valores almacenados, los grabamos:   */
  if (k>1) grabar(&k,k-1,tmem,emem,posmem,relmem);
  if (kpoinc>1) grabar_poinc(&kpoinc,kpoinc-1,xpmem,hyperR);
  cerrar_streams();
  nrerror("Demasiados pasos en la rutina ODEINT_POINC");
}   // termina odeint_poinc()
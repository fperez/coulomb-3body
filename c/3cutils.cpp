/* Archivo con rutinas de utilidad para el programa de integracion de las
   ecuaciones de Hamilton del atomo de Helio.
   Algunas son sin embargo rutinas generales, que pueden emplearse para
   diversos usos.   */

// #includes comunes a todos los archivos:
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include "3cdefs.h"

// #includes que no son comunes a todos
//#include <dos.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>

extern double m1,m2,m3,Mtot,etot;
extern int Z1,Z2,Z3;


/* rutina que lee una tecla solo dentro de un cierto conjunto de teclas
   permitidas (devuelve mayusculas). No hace eco a pantalla:  */
char leatec(const char *opciones)
{ char selec;

  do
  { selec=toupper(getch());
//    if (!strchr(opciones,selec)) beep(1000,0.05);
  } while(!strchr(opciones,selec));
  return selec;
} // termina leatec()

// rutinas de utilidad, tomadas de NUMERICAL RECIPES:
void nrerror(char error_text[])
{
   //	beep(1000,0.1);
	fprintf(stderr,"\n\nError de tiempo de ejecución en la rutina de NUMERICAL RECIPES...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...saliendo al sistema...\n");
	exit(1);
}

double *dvector(int nl,int nh)
{
	double *v;

	v=(double *)malloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("Asignación de memoria fallida en dvector()");
	return v-nl;
}


double **dmatrix(int nrl,int nrh,int ncl,int nch)
{
	int i;
	double **m;

	m=(double **) malloc((unsigned) (nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("Asignación de memoria fallida 1 en dmatrix()");
	m -= nrl;
	
	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("Asignación de memoria fallida 2 en dmatrix()");
		m[i] -= ncl;
	}
	return m;
}


void free_dvector(double *v,int nl,int nh)
{
	free((char*) (v+nl));
}


void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)
{
	int i;

	for(i=nrh;i>=nrl;i--) free((char*) (m[i]+ncl));
	free((char*) (m+nrl));
}
// terminan las utilidades de NUMERICAL RECIPES.


// transformaciOn del sistema fijo en el lab. al del Centro de Masa:
void lab_cm(double *ql,double *qc)
{ int i,j=1,k;
  double R0[4],V0[4],m[4];

  R0[1] = (m1*ql[1]+m2*ql[4]+m3*ql[7])/Mtot;
  R0[2] = (m1*ql[2]+m2*ql[5]+m3*ql[8])/Mtot;
  R0[3] = (m1*ql[3]+m2*ql[6]+m3*ql[9])/Mtot;
  V0[1] = (ql[10]+ql[13]+ql[16])/Mtot;
  V0[2] = (ql[11]+ql[14]+ql[17])/Mtot;
  V0[3] = (ql[12]+ql[15]+ql[18])/Mtot;
  m[1] = m1; m[2] = m2; m[3] = m3;
  // coord. de posiciOn y momentos:
  for(i=1;i<=9;i++)
  { qc[i] = ql[i]-R0[j];
    k=(i-1)/3+1;
    qc[i+9] = ql[i+9]-m[k]*V0[j];
    j++;
    if (j>3) j = 1;
  }
} // termina lab_cm()


// rutina que pasa a coordenadas regularizadas:
void car_reg(double *qc,double *Q)
{ int i;
  double qr[NEQ1];
  // primero pasamos a coordenadas relativas a m3 y momentos de m1 y m2:
  qr[4]=qr[8]=qr[12]=qr[16]=0; // simplemente no se usan.
  for(i=1;i<=3;i++) qr[i]=qc[i]-qc[i+6];
  for(i=4;i<=6;i++) qr[i+1]=qc[i]-qc[i+3];
  for(i=9;i<=11;i++) qr[i]=qc[i+1];    // mirar la pos. de simplificar esto con punteros
  for(i=13;i<=15;i++) qr[i]=qc[i];

  // ahora vamos de las relativas a las regularizadas
  // comenzamos por las coordenadas de posiciOn:
#ifdef SIS_PLANO    // simplificaciones en el caso plano:
  if (qr[1]>=0)
   { Q[1] = sqrt((sqrt(qr[1]*qr[1]+qr[2]*qr[2]) + qr[1])/2);
     Q[2] = qr[2]/2/Q[1];
     Q[3] = 0;
     Q[4] = 0;
   }
  else
   { Q[2] = sqrt((sqrt(qr[1]*qr[1]+qr[2]*qr[2]) - qr[1])/2);
     Q[1] = qr[2]/2/Q[2];
     Q[3] = 0;
     Q[4] = 0;
   }
  if (qr[5]>=0)
   { Q[5] = sqrt((sqrt(qr[5]*qr[5]+qr[6]*qr[6]) + qr[5])/2);
     Q[6] = qr[6]/2/Q[5];
     Q[7] = 0;
     Q[8] = 0;
   }
  else
   { Q[6] = sqrt((sqrt(qr[5]*qr[5]+qr[6]*qr[6]) - qr[5])/2);
     Q[5] = qr[6]/2/Q[6];
     Q[7] = 0;
     Q[8] = 0;
   }
  // y seguimos con los momentos:
  Q[9]  = 2*(Q[1]*qr[9]+Q[2]*qr[10]);
  Q[10] = 2*(-Q[2]*qr[9]+Q[1]*qr[10]);
  Q[13] = 2*(Q[5]*qr[13]+Q[6]*qr[14]);
  Q[14] = 2*(-Q[6]*qr[13]+Q[5]*qr[14]);
  Q[11] = Q[12] = Q[15] = Q[16] = 0;
#else    // ecuaciones grales del caso tridimensional:
  if (qr[1]>=0)
   { Q[1] = sqrt((sqrt(qr[1]*qr[1]+qr[2]*qr[2]+qr[3]*qr[3]) + qr[1])/2);
     Q[2] = qr[2]/2/Q[1];
     Q[3] = qr[3]/2/Q[1];
     Q[4] = 0;
   }
  else
   { Q[2] = sqrt((sqrt(qr[1]*qr[1]+qr[2]*qr[2]+qr[3]*qr[3]) - qr[1])/2);
     Q[1] = qr[2]/2/Q[2];
     Q[3] = 0;
     Q[4] = qr[3]/2/Q[2];
   }
  if (qr[5]>=0)
   { Q[5] = sqrt((sqrt(qr[5]*qr[5]+qr[6]*qr[6]+qr[7]*qr[7]) + qr[5])/2);
     Q[6] = qr[6]/2/Q[5];
     Q[7] = qr[7]/2/Q[5];
     Q[8] = 0;
   }
  else
   { Q[6] = sqrt((sqrt(qr[5]*qr[5]+qr[6]*qr[6]+qr[7]*qr[7]) - qr[5])/2);
     Q[5] = qr[6]/2/Q[6];
     Q[7] = 0;
     Q[8] = qr[7]/2/Q[6];
   }
  // y seguimos con los momentos:
  Q[9]  = 2*(Q[1]*qr[9]+Q[2]*qr[10]+Q[3]*qr[11]);
  Q[10] = 2*(-Q[2]*qr[9]+Q[1]*qr[10]+Q[4]*qr[11]);
  Q[11] = 2*(-Q[3]*qr[9]-Q[4]*qr[10]+Q[1]*qr[11]);
  Q[12] = 2*(Q[4]*qr[9]-Q[3]*qr[10]+Q[2]*qr[11]);
  Q[13] = 2*(Q[5]*qr[13]+Q[6]*qr[14]+Q[7]*qr[15]);
  Q[14] = 2*(-Q[6]*qr[13]+Q[5]*qr[14]+Q[8]*qr[15]);
  Q[15] = 2*(-Q[7]*qr[13]-Q[8]*qr[14]+Q[5]*qr[15]);
  Q[16] = 2*(Q[8]*qr[13]-Q[7]*qr[14]+Q[6]*qr[15]);
#endif
} // termina car_reg()


// transformacion inversa de car_reg():
void reg_car(double *Q,double *qc)
{ int i;
  double qr[NEQ1],R1_2,R2_2,a1,a2;

#ifdef SIS_PLANO    // simplificaciones en el caso plano:
  // pasamos primero a las coordenadas relativas:
  // las de posiciOn:
  qr[4]=qr[8]=qr[12]=qr[16]=0; // simplemente no se usan.
  qr[3] = qr[7] = qr[11] = qr[15] = 0;
  qr[1] = Q[1]*Q[1]-Q[2]*Q[2];
  qr[2] = 2.0*Q[1]*Q[2];
  qr[5] = Q[5]*Q[5]-Q[6]*Q[6];
  qr[6] = 2.0*Q[5]*Q[6];
  // y los momentos:
  R1_2 = 2*(Q[1]*Q[1] + Q[2]*Q[2]);
  R2_2 = 2*(Q[5]*Q[5] + Q[6]*Q[6]);
  qr[9]  = (Q[1]*Q[9]-Q[2]*Q[10])/R1_2;
  qr[10] = (Q[2]*Q[9]+Q[1]*Q[10])/R1_2;
  qr[13] = (Q[5]*Q[13]-Q[6]*Q[14])/R2_2;
  qr[14] = (Q[6]*Q[13]+Q[5]*Q[14])/R2_2;

  // ahora obtenemos las coord. del sistema original (centro de masa):
  // posiciones:
  a1 = -m1/Mtot;
  a2 = -m2/Mtot;
  qc[3] = qc[6] = qc[9] = 0;
  qc[7] = a1*qr[1] + a2*qr[5];
  qc[8] = a1*qr[2] + a2*qr[6];
  qc[1] = qr[1] + qc[7];
  qc[2] = qr[2] + qc[8];
  qc[4] = qr[5] + qc[7];
  qc[5] = qr[6] + qc[8];
  // y momentos:
  for(i=10;i<=12;i++) qc[i] = qr[i-1];
  for(i=13;i<=15;i++) qc[i] = qr[i];
  qc[16] = -qc[10] - qc[13];
  qc[17] = -qc[11] - qc[14];
  qc[18] = -qc[12] - qc[15];
#else    // ecuaciones grales del caso tridimensional:
  // pasamos primero a las coordenadas relativas:
  // las de posiciOn:
  qr[4]=qr[8]=qr[12]=qr[16]=0; // simplemente no se usan.
  qr[1] = Q[1]*Q[1]-Q[2]*Q[2]-Q[3]*Q[3]+Q[4]*Q[4];
  qr[2] = 2.0*(Q[1]*Q[2]-Q[3]*Q[4]);
  qr[3] = 2.0*(Q[1]*Q[3]+Q[2]*Q[4]);
  qr[5] = Q[5]*Q[5]-Q[6]*Q[6]-Q[7]*Q[7]+Q[8]*Q[8];
  qr[6] = 2.0*(Q[5]*Q[6]-Q[7]*Q[8]);
  qr[7] = 2.0*(Q[5]*Q[7]+Q[6]*Q[8]);
  // y los momentos:
  R1_2 = 2*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3] + Q[4]*Q[4]);
  R2_2 = 2*(Q[5]*Q[5] + Q[6]*Q[6] + Q[7]*Q[7] + Q[8]*Q[8]);
  qr[9]  = (Q[1]*Q[9]-Q[2]*Q[10]-Q[3]*Q[11]+Q[4]*Q[12])/R1_2;
  qr[10] = (Q[2]*Q[9]+Q[1]*Q[10]-Q[4]*Q[11]-Q[3]*Q[12])/R1_2;
  qr[11] = (Q[3]*Q[9]+Q[4]*Q[10]+Q[1]*Q[11]+Q[2]*Q[12])/R1_2;
  qr[13] = (Q[5]*Q[13]-Q[6]*Q[14]-Q[7]*Q[15]+Q[8]*Q[16])/R2_2;
  qr[14] = (Q[6]*Q[13]+Q[5]*Q[14]-Q[8]*Q[15]-Q[7]*Q[16])/R2_2;
  qr[15] = (Q[7]*Q[13]+Q[8]*Q[14]+Q[5]*Q[15]+Q[6]*Q[16])/R2_2;

  // ahora obtenemos las coord. del sistema original (centro de masa):
  // posiciones:
  a1 = -m1/Mtot;
  a2 = -m2/Mtot;
  qc[7] = a1*qr[1] + a2*qr[5];
  qc[8] = a1*qr[2] + a2*qr[6];
  qc[9] = a1*qr[3] + a2*qr[7];
  qc[1] = qr[1] + qc[7];
  qc[2] = qr[2] + qc[8];
  qc[3] = qr[3] + qc[9];
  qc[4] = qr[5] + qc[7];
  qc[5] = qr[6] + qc[8];
  qc[6] = qr[7] + qc[9];
  // y momentos:
  for(i=10;i<=12;i++) qc[i] = qr[i-1];
  for(i=13;i<=15;i++) qc[i] = qr[i];
  qc[16] = -qc[10] - qc[13];
  qc[17] = -qc[11] - qc[14];
  qc[18] = -qc[12] - qc[15];
#endif
} // termina reg_car()


// funcion que calcula la energia total del sistema:
void energia(double *qc,double *ener)
{ double R1,R2,R,p1,p2,p3;

  R1 = sqrt((qc[1]-qc[7])*(qc[1]-qc[7]) + (qc[2]-qc[8])*(qc[2]-qc[8])
	    + (qc[3]-qc[9])*(qc[3]-qc[9]));
  R2 = sqrt((qc[4]-qc[7])*(qc[4]-qc[7]) + (qc[5]-qc[8])*(qc[5]-qc[8])
	    + (qc[6]-qc[9])*(qc[6]-qc[9]));
  R = sqrt((qc[1]-qc[4])*(qc[1]-qc[4]) + (qc[2]-qc[5])*(qc[2]-qc[5])
	    + (qc[3]-qc[6])*(qc[3]-qc[6]));
  p1 = qc[10]*qc[10] + qc[11]*qc[11] + qc[12]*qc[12];
  p2 = qc[13]*qc[13] + qc[14]*qc[14] + qc[15]*qc[15];
  p3 = qc[16]*qc[16] + qc[17]*qc[17] + qc[18]*qc[18];
  *ener = p1/2/m1 + p2/2/m2 + p3/2/m3 - Z3*(Z1/R1 + Z2/R2) + Z1*Z2/R;
}  // termina energia()
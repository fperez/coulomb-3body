/* Archivo con rutinas de utilidad para el programa de integracion de las
   ecuaciones de Hamilton del atomo de Helio.
   Algunas on sin embargo rutinas generales, que pueden emplearse para
   diversos usos.   */

#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include "heldefs.h"
#include <stdio.h>
#include <string.h>
#include <ctype.h>
//#include <dos.h>

extern double etot;
extern int   Z;


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
	if (!m) nrerror("Asignación de memoria fallida  1 en dmatrix()");
	m -= nrl;

	for(i=nrl;i<=nrh;i++) {
		m[i]=(double *) malloc((unsigned) (nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("Asignación de memoria fallida  2 en dmatrix()");
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


// rutina que pasa a coordenadas regularizadas:
void xtoq(double *x,double *q)
{ if (x[1]>=0)
   { q[1] = sqrt((sqrt(x[1]* x[1] + x[2]* x[2] ) + x[1])/2);
     q[2] = x[2]/(2* q[1]);
   }
  else
   { q[2] = sqrt((sqrt(x[1]* x[1] + x[2]* x[2] ) - x[1])/2);
     q[1] = x[2]/(2* q[2]);
   }
  q[3] = 2*(q[1]* x[3] + q[2]* x[4]);
  q[4] = 2*(q[1]* x[4] - q[2]* x[3]);

  if (x[5]>=0)
   { q[5] = sqrt((sqrt(x[5]* x[5] + x[6]* x[6] ) + x[5])/2);
     q[6] = x[6]/(2* q[5]);
   }
  else
   { q[6] = sqrt((sqrt(x[5]* x[5] + x[6]* x[6] ) - x[5])/2);
     q[5] = x[6]/(2* q[6]);
   }
  q[7] = 2*(q[5]* x[7] + q[6]* x[8]);
  q[8] = 2*(q[5]* x[8] - q[6]* x[7]);
} // termina xtoq()


// transformacion inversa de xtoq():
void qtox(double *q,double *x)
{ x[1] = q[1]* q[1] - q[2]* q[2];
  x[2] = 2* q[1] * q[2];
  x[3] = (q[1]* q[3]  - q[2]* q[4])/(2*(q[1]* q[1]+q[2]* q[2]));
  x[4] = (q[2]* q[3]  + q[1]* q[4])/(2*(q[1]* q[1]+q[2]* q[2]));
  x[5] = q[5]* q[5] - q[6]* q[6];
  x[6] = 2* q[5] * q[6];
  x[7] = (q[5]* q[7]  - q[6]* q[8])/(2*(q[5]* q[5]+q[6]* q[6]));
  x[8] = (q[6]* q[7]  + q[5]* q[8])/(2*(q[5]* q[5]+q[6]* q[6]));
} // termina qtox()


// funcion que calcula la energia total del sistema:
void energia(double *x,double *etot)
{ double r1,r2,r12,p1,p2;

  r1 = sqrt(x[1]*x[1] + x[2]*x[2]);
  r2 = sqrt(x[5]*x[5] + x[6]*x[6]);
  r12 = sqrt((x[1]-x[5])*(x[1]-x[5]) + (x[2]-x[6])*(x[2]-x[6]));
  p1 = x[3]*x[3] + x[4]*x[4];
  p2 = x[7]*x[7] + x[8]*x[8];
  *etot = p1/2 + p2/2 - Z/r1 - Z/r2 + 1/r12;
}  // termina energia()
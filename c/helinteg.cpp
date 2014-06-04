/* Rutinas de integracion de ecuaciones diferenciales. Tomadas de NUMERICAL
   RECIPES in C:                                                */

#include <stdlib.h> //
#include <math.h>   //
#include <conio.h>    //
#include <fstream.h>
#include "heldefs.h"

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
extern caminoDOS arch_pos[NPART+1], arch_etot, arch_rel[NPART+1], arch_poinc[NEQ/2+1];

// variables globales:
extern double etot;
extern int    Z;

int ind_coord;  // índice de la coordenada cuyo mapa de Poincaré se calcula. 

/* funcion usada en la rutina de integracion, proporciona los lados
   derechos del sistema de ecuaciones diferenciales: */
void derivs_norm(double t,double *x,double *dxdt)
{ double r1,r2,a,b,r12,p1,p2;

  r1 = x[1]*x[1] + x[2]*x[2];
  r2 = x[5]*x[5] + x[6]*x[6];
  a = x[1]*x[1] - x[2]*x[2] - x[5]*x[5] + x[6]*x[6];
  b = x[1]*x[2] - x[5]*x[6];
  r12 = sqrt(a*a+4*b*b);
  p1 = x[3]*x[3] + x[4]*x[4];
  p2 = x[7]*x[7] + x[8]*x[8];

/* Nota: etot es una vble global, que contiene la energia total del sistema
	 calculada en t=0, y que DEBE ser una constante para todo t. */
  dxdt[1] = r2*x[3]/4;
  dxdt[2] = r2*x[4]/4;
  dxdt[3] = -2*(p2/8 - Z + r2*(1/r12-etot))*x[1] +
	     2*r1*r2*(a*x[1] + 2*b*x[2])/r12/r12/r12;
  dxdt[4] = -2*(p2/8 - Z + r2*(1/r12-etot))*x[2] +
	     2*r1*r2*(2*b*x[1] - a*x[2])/r12/r12/r12;
  dxdt[5] = r1*x[7]/4;
  dxdt[6] = r1*x[8]/4;
  dxdt[7] = -2*(p1/8 - Z + r1*(1/r12-etot))*x[5] -
	     2*r1*r2*(a*x[5] + 2*b*x[6])/r12/r12/r12;
  dxdt[8] = -2*(p1/8 - Z + r1*(1/r12-etot))*x[6] -
	     2*r1*r2*(2*b*x[5] - a*x[6])/r12/r12/r12;
}   // termina derivs_norm()


/* lados derechos del sistema de ecuaciones diferenciales, ajustados para cálculos
   de secciones de Poincaré:            */
void derivs_poinc(double t,double *x,double *dxdt)
{ double r1,r2,a,b,r12,p1,p2;
  double *f,dSq1,dSq2,dSq5,dSq6,fmas;
  int    i;

  f = dvector(1,NEQ);
  r1 = x[1]*x[1] + x[2]*x[2];
  r2 = x[5]*x[5] + x[6]*x[6];
  a = x[1]*x[1] - x[2]*x[2] - x[5]*x[5] + x[6]*x[6];
  b = x[1]*x[2] - x[5]*x[6];
  r12 = sqrt(a*a+4*b*b);
  p1 = x[3]*x[3] + x[4]*x[4];
  p2 = x[7]*x[7] + x[8]*x[8];

/* Nota: etot es una vble global, que contiene la energia total del sistema
	 calculada en t=0, y que DEBE ser una constante para todo t. */

  // funciones que definen el lado derecho del sistema de ecuaciones diferenciales:
  f[1] = r2*x[3]/4;
  f[2] = r2*x[4]/4;
  f[3] = -2*(p2/8 - Z + r2*(1/r12-etot))*x[1] +
	     2*r1*r2*(a*x[1] + 2*b*x[2])/r12/r12/r12;
  f[4] = -2*(p2/8 - Z + r2*(1/r12-etot))*x[2] +
	     2*r1*r2*(2*b*x[1] - a*x[2])/r12/r12/r12;
  f[5] = r1*x[7]/4;
  f[6] = r1*x[8]/4;
  f[7] = -2*(p1/8 - Z + r1*(1/r12-etot))*x[5] -
	     2*r1*r2*(a*x[5] + 2*b*x[6])/r12/r12/r12;
  f[8] = -2*(p1/8 - Z + r1*(1/r12-etot))*x[6] -
	     2*r1*r2*(2*b*x[5] - a*x[6])/r12/r12/r12;

  // Construccion de la funcion f[N+1] con las derivadas de la sup de seccion
  dSq1 = 4.0*x[5]*x[6]*x[1]-2.0*x[2]*(x[5]*x[5]-x[6]*x[6]);
  dSq2 = -4.0*x[5]*x[6]*x[2]-2.0*x[1]*(x[5]*x[5]-x[6]*x[6]);
  dSq5 = 2.0*x[6]*(x[1]*x[1]-x[2]*x[2])-4.0*x[1]*x[2]*x[5];
  dSq6 = 2.0*x[5]*(x[1]*x[1]-x[2]*x[2])+4.0*x[1]*x[2]*x[6];
  fmas = f[1]*dSq1+f[2]*dSq2+f[5]*dSq5+f[6]*dSq6;
  // asignación del lado derecho del sistema ajustando para el cálculo de sec. de Poincaré:
  for (i=1;i<=NEQ;i++) dxdt[i] = f[i]/fmas;
  free_dvector(f,1,NEQ);
}   // termina derivs_poinc()


// Metodo de Runge-Kutta de 4o orden:
void rk4(double y[],double dydx[],int n,double x,double h,double yout[],
	 void (*derivs)(double,double *,double *) )
{ int i;
  double xh,hh,h6,dym[DIM],dyt[DIM],yt[DIM];

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
  double  dysav[DIM],ysav[DIM],ytemp[DIM];

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
ofstream a_pos[NPART+1],a_rel[NPART+1],a_etot,a_poinc[NEQ/2+1],a_hyper;

void abrir_streams(int pnuevo)
{ int i,modo,out,app;

  out  = 0x02;        // open for writing
  app  = 0x08;        // append mode: all additions at eof

  if (pnuevo) modo = out;
    else modo = app;
  for (i=1;i<=NPART;i++)
  { a_pos[i].open(arch_pos[i]);
    a_rel[i].open(arch_rel[i]);
  }
  a_etot.open(arch_etot);
  for (i=1;i<=NEQ/2;i++) a_poinc[i].open(arch_poinc[i],modo);
  a_hyper.open("hyper.dat",modo);
}  // termina abrir_streams()


// rutina que graba los vectores de almacenamiento temporal a disco:
void grabar(int *cont,int lim,double *tdat, double *edat,double ***posdat,double **reldat)
{ int i,j;

  for (i=1;i<=lim;i++)
  { // Posicion y energia:
    for (j=1;j<=NPART;j++)
    { a_pos[j]  << posdat[j][i][1] << "   " << posdat[j][i][2] << "\n"; // posiciones
      a_rel[j]    << tdat[i]  << "   " << reldat[j][i]   << "\n";  // coord. relativas
    }
    a_etot  << tdat[i]    << "   " << edat[i] << "\n";  // energía(t) -- sup. constante.
  }
  cout << "Grabando Parámetros\n";
  *cont=1;
}  // termina grabar()


// rutina que graba los vectores de Poincaré a disco:
void grabar_poinc(int *cont,int lim, double **xpdat, double **hyperdat)
{ int i/*,j*/;

  for (i=1;i<=lim;i++)
  { a_poinc[1] << xpdat[i][1] << "   " << xpdat[i][3] << "\n";
    a_poinc[2] << xpdat[i][2] << "   " << xpdat[i][4] << "\n";
    a_poinc[3] << xpdat[i][5] << "   " << xpdat[i][7] << "\n";
    a_poinc[4] << xpdat[i][6] << "   " << xpdat[i][8] << "\n";
    a_hyper << hyperdat[i][1] << "   " << hyperdat[i][2] << "\n";

  /* for (j=1;j<=NEQ/2;j++)
      a_poinc[j] << xpdat[i][2*j] << "   " << xpdat[i][2*j-1] << "\n";*/
  }
  cout << "Grabando datos de Poincaré\n";
  *cont=1;
}  // termina grabar_poinc()


// rutina que cierra los streams:
void cerrar_streams(void)
{ int i;

  for (i=1;i<=NPART;i++)
  { a_pos[i].close();
    a_rel[i].close();
  }
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


/* Mapas de Poincaré:
   Rutina analoga a la anterior, pero que además genera una sección
   de Poincaré.    */

void odeint_poinc(double *qstart,int nvar,double t1,double t2,double dtsav,
		  double eps,double h1,double hmin,int *nok,int *nbad,
		  void (*derivs)(double,double *,double *),
		  void (*integ)(double *,double *,int,double *,double,
				double,double *,double *,double *,
				void (*)(double,double *,double *)),
		  int calc_poinc,int sign_mom,int pnuevo)
{ // declaracion de vbles:
  int     nstp,i,ip;
  double  tsav,t,tvie,w,hnext,hdid,h,tmax=fabs(t2-t1);
  double  *qscal,*q,*dqdt,*xn,*qp,*xp,*dqpdt;

  // vbles para almacenar valores en memoria antes de grabar a disco:
  double  *tmem,*emem,**relmem,**xpmem,**hyperR,***posmem;
  // kmax: numero maximo de valores intermedios que se almacenan
  int     kmax=8000,k=1,kpoinc=1;
  // control de autoionizacion:
  int   ionizo=0;
  // distancia a partir de la cual se considera que hay ionizacion:
  float umbr_ioniz=20;
  // coordenadas relativas y energia intermedia (para control):
  double r1,r2,eint;
  // calculos de mapas de Poincare:
  double Svie,Snue,mom,hyperRad;


  // COMIENZA EL CODIGO DE odeint_poinc():
  abrir_streams(pnuevo);
  // se reserva memoria para los vectores de calculo:
  qscal = dvector(1,nvar);
  q = dvector(1,nvar);
  xn = dvector(1,nvar);      // xn = coordenadas no regularizadas.
  dqdt = dvector(1,nvar);
  qp = dvector(1,nvar);
  xp = dvector(1,nvar);
  dqpdt = dvector(1,nvar);

  // se reserva memoria para los vectores y matrices de almacenamiento:
  tmem = dvector(1,kmax);
  emem = dvector(1,kmax);
  relmem = dmatrix(1,NPART,1,kmax);

  posmem = (new double**[NPART])-1;
  for (i=1;i<=NPART;i++) posmem[i]=dmatrix(1,kmax,1,DIM_ESP);
  xpmem = dmatrix(1,kmax,1,nvar);
  hyperR = dmatrix(1,kmax,1,2);

  // inicializacion de variables de la rutina odeint original:
  t=t1;
  h=(t2 > t1) ? fabs(h1) : -fabs(h1);
  dtsav=fabs(dtsav);
  *nok = (*nbad) = 0;
  for (i=1;i<=nvar;i++) q[i]=qstart[i];

  // se inicializa la vble. para el control del cálculo de sec. de Poincaré:
//  ind_coord=ind_coord_p;   // ind_coord es global, se requiere así para otras rutinas.
  qtox(q,xp);
  Svie=Snue=0;
  // garantizamos que el primer paso se almacene:
  tsav=t-dtsav*2.0;
  cout << "\nCALCULANDO:\n";

  // CICLO DE CALCULO:
  for(nstp=1;nstp<=MAXSTP;nstp++)
  {
    // Modulos de los vectores r1 y r2:
    r1=q[1]*q[1] + q[2]*q[2];
    r2=q[5]*q[5] + q[6]*q[6];
    // detectar autoionizacion:
    if ( (r1>umbr_ioniz) || (r2>umbr_ioniz) ) ionizo=1;

    // indicar que el computador esta 'vivo':
    if (nstp%250==0) cout << (int)(100*(t/tmax)) << " %\n";

    // se ALMACENAN EN MEMORIA valores solo cada cierto intervalo (dado por dtsav):
    if (!calc_poinc)  // solo se hace si no se estan calculando mapas de Poincare
      if ( fabs(t-tsav) > dtsav )
      { // Posicion y energia:
	qtox(q,xn);
	energia(xn,&eint);
	tmem[k]=t;
	emem[k]=eint;
	posmem[1][k][1]=xn[1]; posmem[1][k][2]=xn[2];
	posmem[2][k][1]=xn[5]; posmem[2][k][2]=xn[6];
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
    for (i=1;i<=nvar;i++) qscal[i]=fabs(q[i])+fabs(dqdt[i]*h)+TINY;

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
      qtox(q,xp);
      Snue = xp[1]*xp[6]-xp[2]*xp[5]; // Sup: theta(ang. entre e-'s)=Pi.
      mom = xp[2]*xp[8]+xp[6]*xp[7];  // calcular el momento canonico ptheta.
      // evaluar cambio de signo en S:
      if ( (Svie*Snue <= 0)   && (mom*sign_mom > 0) )
      { //calculamos el punto en el mapa de Poincaré con vbles auxiliares:
       for (ip=1;ip<=nvar;ip++) qp[ip]=q[ip];
       w=Snue;   		// parámetro de "tiempo"
       derivs_poinc(t,qp,dqpdt);
       rk4(qp,dqpdt,nvar,w,-Snue,qp,derivs_poinc);
       qtox(qp,xp);
       for (ip=1;ip<=nvar;ip++) xpmem[kpoinc][ip]=xp[ip];
       // mapa de Poincare del HyperRadio - Momento canonico asociado.
       hyperRad = hyperR[kpoinc][1] = sqrt(xp[1]*xp[1]+xp[2]*xp[2]+xp[5]*xp[5]+xp[6]*xp[6]);
       hyperR[kpoinc][2] = (xp[1]*xp[3]+xp[2]*xp[4]+xp[5]*xp[7]+xp[6]*xp[8])/
			   hyperRad;

       kpoinc++;
      }
    } // terminan los calculos de Ponincare

    // si terminó, se interrumpió la integración o hubo ionización, salir:
    if ( ((t-t2)*(t2-t1) >= 0.0) || (kbhit()) || ionizo)
    { // si hay valores almacenados, grabarlos:
      if (k>1) grabar(&k,k-1,tmem,emem,posmem,relmem);
      if (kpoinc>1) grabar_poinc(&kpoinc,kpoinc-1,xpmem,hyperR);

      // si se pulsó una tecla o hubo ioniz., indicamos el tiempo alcanzado:
      if (kbhit() || ionizo)
       { if (ionizo)
	   cout << "\n\nCALCULO INTERRUMPIDO POR IONIZACION EN t = " << t;
	 else
	   { cout << "\n\nCALCULO INTERRUMPIDO POR EL USUARIO EN t = " << t;
	     getch();
	   }
	 cout << ". (Se tenía tmax = " << t2 << ").";
       }
      // devolver el vector de vbles dependientes en su estado actual:
      for (i=1;i<=nvar;i++) qstart[i]=q[i];
      // liberar la memoria asignada dinamicamente:
      free_dvector(dqdt,1,nvar);  free_dvector(q,1,nvar);
      free_dvector(xn,1,nvar);    free_dvector(qscal,1,nvar);
      free_dvector(xp,1,nvar);    free_dvector(qp,1,nvar);
      free_dvector(tmem,1,kmax);  free_dvector(emem,1,kmax);
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
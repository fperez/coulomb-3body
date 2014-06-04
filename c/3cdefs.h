/* Archivo con definiciones de constantes, estructuras de datos y prototipos
   de funciones usadas en el programa de integraciOn de las ecuaciones de
   Hamilton para el problema coulombiano de 3 cuerpos.   */

// activar esta linea para usar las simplificaciones del caso plano:
#define SIS_PLANO

// forzar precision simple para acelerar un poco los calculos:
//#define double float

#define NEQ 16        // numero de ecuaciones diferenciales del sistema
#define NEQ1 NEQ+1
#define NPART 3      // numero de part_culas del sistema
#define DIM_ESP 3    // # de dimensiones espaciales del problema
#define DIM_FAS 18  // (NPART*DIM_ESP*2): DimensiOn del espacio fAsico
#define DIM_FAS1 DIM_FAS+1

// define's de la rutina odeint_poinc():
#define MAXSTP 100000  // numero maximo de pasos permitido en odeint.
#define TINY 1.0e-8

// define's de la rutina bsstep():
#define IMAX 11
#define NUSE 7
#define SHRINK 0.95
#define GROW 1.2

// se define un tipo capaz de contener todos los parametros de calculo:
typedef
  struct
   { double cini_lab[DIM_FAS1];
     double cini[DIM_FAS1];
     double cact[DIM_FAS1];
     double tmax;
     double eps;
     double inte_imag;
     int    batch;
     // las vbles pxxx se refieren a los calculos de secciones de Poincare:
     int    pflag;  // calcular o no sec. de Poincare.
     int    pcoord; // coordenada que define el plano de seccion
     double pval;   // valor de la coordenada anterior.
     int    pmom;   // momento asociado a la coordenada (indice q'lo define)
     int    psign_mom;  // signo del momento.
     int    pnue;   // comenzar o no una nueva seccion
     char   opcion;
   } parametros;

typedef char caminoDOS[51];

// rutinas de utilidad, tomadas de NUMERICAL RECIPES:
void nrerror(char error_text[]);
double *dvector(int nl,int nh);
double **dmatrix(int nrl,int nrh,int ncl,int nch);
void free_dvector(double *v,int nl,int nh);
void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
// terminan las utilidades de NUMERICAL RECIPES.

// lados derechos del sistema de ecuaciones diferenciales que vamos a integrar:
void derivs_norm(double tt,double *xx,double *dxxdt);
void derivs_poinc(double tt,double *xx,double *dxxdt);

// Metodo de Runge-Kutta de 4° orden:
void rk4(double y[],double dydx[],int n,double x,double h,double yout[],
	 void (*derivs)(float,float *,float *) );

// rutina de Runge-Kutta con control de error:
void rkqc(double y[],double dydx[],int n,double *x,double htry,
	  double eps,double yscal[],double *hdid,double *hnext,
	  void (*derivs)(double,double *,double *) );

/* Rutina que integra una ecuacion diferencial ordinaria a lo largo de un
   cierto intervalo de la variable independiente, usando el metodo de
   Runge-Kutta. Adaptada de NUMERICAL RECIPES:                          */
void odeint_poinc(double *qstart,int nvar,double t1,double t2,double dtsav,
		  double eps,double h1,double hmin,int *nok,int *nbad,
		  void (*derivs)(double,double *,double *),
		  void (*integ)(double *,double *,int,double *,double,
				double,double *,double *,double *,
				void (*)(double,double *,double *)),
		  int calc_poinc,int sign_mom,int pnuevo);

// "pito" para usos varios (frec en Hz, tiempo en seg):
//void beep(unsigned frec,float tiempo);
// transformaciOn del sistema fijo en el lab. al del Centro de Masa:
void lab_cm(double *ql,double *qc);
// rutina que pasa a coordenadas regularizadas:
void car_reg(double *q,double *Q);
// transformacion inversa de car_reg():
void reg_car(double *Q,double *q);
// funcion que calcula la energia total del sistema:
void energia(double *qc,double *ener);
/* rutina que lee una tecla solo dentro de un cierto conjunto de teclas
   permitidas (devuelve mayusculas). No hace eco a pantalla:  */
char leatec(const char *opciones);
// presenta un el menu ppal. del  programa:
char menu(parametros *par);
// muestra los valores de todos los parametros vigentes:
void muestre_param(parametros *par);
// rutina de edicion de parametros:
void editar(parametros *par);
// rutina que lee un archivo de configuracion inicial de disco:
void config_ini(char *arch,parametros *par_ini);
// Seleccion de parametros para calcular con el modelo asincrono:
void orb_Langmuir(parametros *par_as);

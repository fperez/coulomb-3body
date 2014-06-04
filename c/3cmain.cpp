/* Programa que calcula orbitas clasicas para el problema coulombiano de
   3 cuerpos usando coordenadas regularizadas.

   Fernando Pérez, Departamento de Física, Universidad de Antioquia.
   Medellín, 1994.                                                       */

// #includes comunes a todos los archivos:
#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include "3cdefs.h"

// #includes que no son comunes a todos
#include <fstream.h>
#include <iostream.h>
#include <stdio.h>
#include <string.h>


// Nombres DOS de los archivos utilizados:
// Directorio para datos:
caminoDOS dirdatos="";
/* Coordenadas de posicion de c/part (salida); Energia total del sistema, supuesta
   constante; Coordenadas relativas; Mapas de Poincaré: */
caminoDOS arch_pos[NPART+1], arch_etot, arch_rel, arch_poinc[NEQ/2+1];

// variables globales (declaradas extern en los demas archivos):
double  etot;
double m1,m2,m3,mu13,mu23,Mtot;
int Z1,Z2,Z3,maxmem;

// Comienzan las funciones (la rutina main() esta al final):

// inicializa los nombres de los archivos que se van a usar:
void ini_arch()
{ int i;

  // Coordenadas de posicion (salida):
  for (i=1;i<=NPART;i++) sprintf(arch_pos[i],"%spos%d.dat",dirdatos,i);
  // Coordenadas relativas:
  sprintf(arch_rel,"%srel.dat",dirdatos);
   // Energia total del sistema, supuesta constante (salida):
  sprintf(arch_etot,"%setot.dat",dirdatos);
  // Archivos que reciben la informacion de los mapas de Poincare:
  for (i=1;i<=NEQ/2;i++)
    sprintf(arch_poinc[i],"%spoinc%d.dat",dirdatos,i);
} // termina iniarch()


// funciOn para probar las transformacines de coordenadas:
void pruebe(parametros *par)
{ double  ql[DIM_FAS1],qc1[DIM_FAS1],qc2[DIM_FAS1],qc3[DIM_FAS1];
  double  Q0[NEQ1],Q1[NEQ1],ener1,ener2,ener3;
  int i;

  for(i=1;i<=DIM_FAS;i++) ql[i]=par->cini_lab[i];
  lab_cm(ql,qc1);
  energia(qc1,&ener1);
  car_reg(qc1,Q0);
  reg_car(Q0,qc2);
  energia(qc2,&ener2);
  car_reg(qc2,Q1);
  reg_car(Q1,qc3);
  energia(qc3,&ener3);
} // termina pruebe()


// rutina de calculo del programa:
void calcule_tray(parametros *par)
{ /* qc: coordenadas cartesianas con respecto al sistema fijo en el lab.;
     Q:     coordenadas regularizadas.   */
  double  qc[DIM_FAS1],Q[NEQ1],paso,pasomin;
  double  eini,efin,deltaE,deltaE_r,tini;
  int     i,nbien,nmal;

  // codigo de la fn calcule():

  // colocamos las variables locales en sus valores para calcular:
  if (par->opcion=='8')
    for(i=1;i<=DIM_FAS;i++) qc[i]=par->cini[i];
  else
    for(i=1;i<=DIM_FAS;i++) qc[i]=par->cact[i];
  paso=0.01;
  pasomin=0.0;
  tini=0;
  energia(qc,&eini);
  etot = eini;
  cout << "\n\n*************************************\n";
  cout << "Pulse cualquier tecla para suspender.\n";

  cout << "\nEnergia Total inicial = " << eini << "\n";

  // pasamos a coordenadas regularizadas antes de integrar:
  car_reg(qc,Q);

//  A continuacion, se llama la rutina de integracion (con Runge-Kutta):
  odeint_poinc(Q,NEQ,tini,par->tmax,par->inte_imag,par->eps,paso,pasomin,
	       &nbien,&nmal,derivs_norm,rkqc,
	       par->pflag,par->psign_mom,par->pnue);
  reg_car(Q,qc);
  energia(qc,&efin);
  cout << "\n\nEnergía Total final   = " << efin <<"\n";

  // calculamos el cambio de energia inducido por las rutinas numericas:
  deltaE=fabs(efin-eini);
  deltaE_r=100*deltaE/fabs(eini);
  cout << "\n|Delta Energía|       = " << deltaE << "\n";
  cout << "|Delta Energía|(rel)  = " << deltaE_r << " %\n\n";
  cout << "Número de pasos exitosos: " << nbien << "\n";
  cout << "Número de pasos fallidos: " << nmal  << "\n";


  // actualizamos la estructura de parametros con los ultimos valores
  // de las distintas vbles dinamicas del problema:
  for(i=0;i<=DIM_FAS;i++) par->cact[i]=qc[i];

  // esperar para permitir la lectura de los parametros de energia:
//  beep(1000,0.2);
  if (!par->batch)
  { cout << "\nPulse cualquier tecla para volver al menú.";
    getch();
  }
} // termina calcule_tray()


// presenta el menu ppal. del  programa:
char menu(parametros *par)
{ char   selec;
  char   opciones[] = "123456789LPS";
  char   d[]        = "     ";
  double *qini;

  qini=&par->cini_lab[0];

  clrscr();
// indicamos en pantalla si el cálculo es en dos o tres dimensiones:
  cout << "\n" << d << "**********************************************************************\n";
#ifdef SIS_PLANO
  cout << d << "* CALCULOS CLASICOS EN 2-D PARA EL PROBLEMA COULOMBIANO DE 3 CUERPOS *\n";
#else
  cout << d << "* CALCULOS CLASICOS EN 3-D PARA EL PROBLEMA COULOMBIANO DE 3 CUERPOS *\n";
#endif
  cout << d << "**********************************************************************\n\n";
  cout << "  1. Editar las posiciones iniciales (sistema del laboratorio):\n";
  cout << "     x1 = "<<qini[1]<<"  y1 = "<<qini[2]<<"  z1 = "<<qini[3]<<"\n";
  cout << "     x2 = "<<qini[4]<<"  y2 = "<<qini[5]<<"  z2 = "<<qini[6]<<"\n";
  cout << "     x3 = "<<qini[7]<<"  y3 = "<<qini[8]<<"  z3 = "<<qini[9]<<"\n";
  cout << "  2. Editar los momentos iniciales (sistema del laboratorio):\n";
  cout << "     px1 = "<<qini[10]<<"  py1 = "<<qini[11]<<"  pz1 = "<<qini[12]<<"\n";
  cout << "     px2 = "<<qini[13]<<"  py2 = "<<qini[14]<<"  pz2 = "<<qini[15]<<"\n";
  cout << "     px3 = "<<qini[16]<<"  py3 = "<<qini[17]<<"  pz3 = "<<qini[18]<<"\n";
  cout << "  3. Leer configuración de disco.\n";
  cout << "  4. Cambiar tmax. Actualmente, tmax = "<<par->tmax<<"\n";
  cout << "  5. Cambiar epsilon. Actualmente, epsilon = "<<par->eps<<"\n";
  cout << "  6. Cambiar masas y cargas: m1="<<m1<<", m2="<<m2<<", m3="<<m3
       << ", Z1="<<Z1<<", Z2="<<Z2<<", Z3="<<Z3<<"\n";
  cout << "  7. Cambiar el intervalo entre imágenes. Actualmente, inte_imag = "
       << par->inte_imag<<"\n";
  cout << "  8. RECALCULAR DESDE LAS COND. INICIALES.\n";
  cout << "  9. PROSEGUIR  DESDE LAS COND. ACTUALES.\n";
  cout << " <L>. Calcular órbitas tipo Langmuir.\n";
  cout << " <P>. Activar cálculo de Secciones de Poincaré: " << par->pflag << "\n";
  cout << " <S>. Salir del programa.\n\n";
  cout << "Elija una opción: ";
  selec=leatec(opciones);
  cout << selec;
  return selec;
}  // termina menu()


// muestra los valores de todos los parametros vigentes:
void muestre_param(parametros *par)
{ double *qini;

  qini=&par->cini_lab[0];
  cout << "La configuración del programa es en este momento:\n";
  cout << "  Posiciones iniciales (sistema del laboratorio):\n";
  cout << "   (1) x1 = "<<qini[1]<<"  (2) y1 = "<<qini[2]<<"  (3) z1 = "<<qini[3]<<"\n";
  cout << "   (4) x2 = "<<qini[4]<<"  (5) y2 = "<<qini[5]<<"  (6) z2 = "<<qini[6]<<"\n";
  cout << "   (7) x3 = "<<qini[7]<<"  (8) y3 = "<<qini[8]<<"  (9) z3 = "<<qini[9]<<"\n";
  cout << "\n  Momentos iniciales (sistema del laboratorio):\n";
  cout << "   (1) px1 = "<<qini[10]<<" (2) py1 = "<<qini[11]<<" (3) pz1 = "<<qini[12]<<"\n";
  cout << "   (4) px2 = "<<qini[13]<<" (5) py2 = "<<qini[14]<<" (6) pz2 = "<<qini[15]<<"\n";
  cout << "   (7) px3 = "<<qini[16]<<" (8) py3 = "<<qini[17]<<" (9) pz3 = "<<qini[18]<<"\n";
  cout <<"\n  Masas y cargas:\n     (1) m1="<<m1<<"   (2) m2="<<m2<<"   (3) m3="<<m3<<"   (4) Z1="
       << Z1<<"   (5) Z2="<<Z2<<"   (6) Z3="<<Z3<<"\n";
  cout <<"  tmax = "<<par->tmax<<";   epsilon = "<<par->eps<<
	 ";   Intervalo entre imágenes = "<<par->inte_imag<<".\n";
}  // termina muestre_param()


// rutina de edicion de parametros:
void editar(parametros *par)
{ int   op,op_cond,ord_eps,salir=0;
  char  act[]         = "\nActualmente, ";
  char  nuev[]        = "Entre el nuevo valor: ";
  char  cond_permit[] = "0123456789";
  char  cond_permit2[] = "0123456";

  op=par->opcion;
  do
  { clrscr();
    cout <<"EDICION DE PARAMETROS:\n";
    cout <<"----------------------\n\n";
    muestre_param(par);
    cout << "\n			 *********************\n";
    if ((op=='1')||(op=='2'))
    { cout << "\nEdición de ";
      if (op=='1')  cout << "POSICIONES. ";
	else cout << "MOMENTOS. ";
      cout << "\nPulse <0> para abandonar la edición.\n";
      cout << "Elija (1-9) qué condición desea alterar: ";
      op_cond=leatec(cond_permit)-48;    // ajuste con el codigo ASCII
      cout<<op_cond <<"\n";
    }
    if (op=='6')
    { cout << "\nEdición de masas y cargas:\n";
      cout << "Pulse <0> para abandonar la edición.\n";
      cout << "Elija (1-6) cuál parámetro desea alterar: ";
      op_cond=leatec(cond_permit2)-48;    // ajuste con el codigo ASCII
      cout<<op_cond <<"\n";
    }
    // En edici_n de pos., mom., masas o cargas es posible editar mßs de un parßmetro:
    if ((op>50) && (op!=54)) salir=1; 
      else if (op_cond==0) return;
    switch(op)
    { case '1':
	cout << "La posicion inicial elegida vale: " << par->cini_lab[op_cond] <<"\n";
	cout << nuev;
	cin >> par->cini_lab[op_cond];
      break;
      case '2':
	cout << "El momento inicial elegido vale: " << par->cini_lab[op_cond+9] <<"\n";
	cout << nuev;
	cin >> par->cini_lab[op_cond+9];
      break;
      case '4':
	cout << act << "tmax = " << par->tmax <<"\n";
	cout << nuev;
	cin >> par->tmax;
      break;
      case '5':
	cout << act << "epsilon = " << par->eps <<"\n";
	cout << "Entre el nuevo ORDEN (número ENTERO) de epsilon: ";
	cin >> ord_eps;
	par->eps=pow10(ord_eps);
      break;
      case '6':
	cout << nuev;
	switch(op_cond)
	{ case 1: cin >> m1; break;
	  case 2: cin >> m2; break;
	  case 3: cin >> m3; break;
	  case 4: cin >> Z1; break;
	  case 5: cin >> Z2; break;
	  case 6: cin >> Z3; break;
	} // fin del switch(op_cond)
	// actualizamos las variables globales asociadas a las masas:
	mu13=m1*m3/(m1+m3); // las masas reducidas con respecto a m3
	mu23= m2*m3/(m2+m3);
	Mtot = m1+m2+m3;
      break;
      case '7':
	cout << act << "el intervalo entre imágenes es: " << par->inte_imag <<"\n";
	cout << nuev;
	cin >> par->inte_imag;
      break;
    }  // fin del switch(op)
  // se actualizan los parámetros con respecto al centro de masa:
  lab_cm(par->cini_lab,par->cini);
  } // fin del do
  while(!salir);
} // termina editar()


// rutina que lee un archivo de configuracion inicial de disco:
void config_ini(char *arch,parametros *par_ini)
{ ifstream   a_config(arch);
  int        i;

  if (!a_config)
  { cout << "\n*********** ERROR ***********\n";
    //beep(1000,0.75);
    cout << "No es posible leer el archivo de configuración.\n";
    cout << "No se cargó ninguna variable.\n";
    cout << "Si el programa se ejecuta por primera vez, configúrelo manualmente,\n";
    cout << "pues puede haber \"basura\" en los distintos parámetros.\n";
    cout << "Si no está comenzando la ejecución, simplemente se mantendrá la\n";
    cout << "configuración que estuviese vigente.\n\n";
    cout << "Pulse cualquier tecla para seguir.";
    getch();
    return;
  };
  // se leen las condiciones iniciales:
  for (i=1;i<=DIM_FAS;i++) a_config >> par_ini->cini_lab[i];
  // se leen los demas parametros:
  a_config >> m1 >> m2 >> m3;
  mu13=m1*m3/(m1+m3); // las masas reducidas con respecto a m3
  mu23= m2*m3/(m2+m3);
  Mtot = m1+m2+m3;
  a_config >> Z1 >> Z2 >> Z3;
  a_config >> par_ini->tmax  >> par_ini->eps >> par_ini->inte_imag >> maxmem;
  a_config >> par_ini->pflag >> par_ini->batch;
  // se pasan las coordenadas al sistema del centro de masa:
  lab_cm(par_ini->cini_lab,par_ini->cini);
  // cerramos el archivo:
  a_config.close();

  par_ini->opcion='8';

  // ******************
  // INICIALIZACION DE LOS PARAMETROS DE POINCARE
  par_ini->pcoord=1;
  par_ini->pval=0.0;
  par_ini->pmom=3;
  par_ini->psign_mom=1;
  par_ini->pnue=0;
  // ******************

}  // termina config_ini()

// Seleccion de parametros para calcular órbitas tipo Langmuir:
void orb_Langmuir(parametros *par_lan, double a0, double b0, double E0)
{ double  a,b,c,d,eLang,emin,alfa;
  int     i;

  // Configuración interactiva:
  if (a0 == 0)    
  { clrscr();
    cout << "\n\n\nENTRADA DE DATOS PARA EL MODELO DE LANGMUIR:\n";
    cout << "--------------------------------------------\n\n";
    cout << "Valores de a y b: ";
    cin  >> a >> b;
    alfa = Z1*(Z1/2.0-2*Z3)/a;
    emin = b*b/m1 + alfa;
    cout << "\nEnergía total del sistema.\n";
    do
    { cout << "Para evitar velocidades imaginarias, el valor de la energía debe\n";
      cout << "ser mayor o igual que: " << emin << "\n";
      cout << "Además, recuerde que en un sistema ligado, Etot<0.\n\n";
      cout << "Entre Etot: ";
      cin  >> eLang;
      cout << "\n";
    } while ( (eLang<emin) || (eLang>0) );
  }
  // Configuración por lotes:
  else          
  { a = a0;
    b = b0;
    alfa = Z1*(Z1/2.0-2*Z3)/a;
    eLang = E0;
  }
  c = sqrt(m3*((m1*(eLang-alfa)-b*b)/(2*m1+m3)));
  d = -2*c;

  // se actualiza la estructura de parametros:
  // los que son nulos:
  par_lan->cini_lab[2]=par_lan->cini_lab[3]=par_lan->cini_lab[5]=par_lan->cini_lab[6]=0;
  par_lan->cini_lab[7]=par_lan->cini_lab[8]=par_lan->cini_lab[9]=0;
  par_lan->cini_lab[12]=par_lan->cini_lab[15]=par_lan->cini_lab[16]=par_lan->cini_lab[18]=0;
  // los no nulos:
  par_lan->cini_lab[1]=-a;
  par_lan->cini_lab[4]=a;
  par_lan->cini_lab[10]=-b;
  par_lan->cini_lab[13]=b;
  par_lan->cini_lab[11]=par_lan->cini_lab[14]=c;
  par_lan->cini_lab[17]=d;
  lab_cm(par_lan->cini_lab,par_lan->cini);

  // los parametros "actuales" se igualan a las cond. iniciales:
  for(i=1;i<=DIM_FAS;i++) par_lan->cact[i]=par_lan->cini[i];
  // se indica recalcular:
  par_lan->opcion='8';
} // termina orb_Langmuir()


void Poinc_Langmuir(parametros *par)
{ double  a,amin,amax,dela,b,bmin,bmax,delb,E;

  par->batch = 1;
  par->pflag = 1;
  E = -1;
  amin = 0.6; amax = 1.3; dela = 0.05;
  bmin = -0.3; bmax = 0.3; delb = 0.05;
  for(a=amin;a<=amax;a+=dela)
    for(b=bmin;b<=bmax;b+=delb)
    { orb_Langmuir(par,a,b,E);
      calcule_tray(par);
    }
} // termina  Poinc_Langmuir()



// ******************************
// RUTINA PRINCIPAL DEL PROGRAMA:
void main(int argc,char *argv[])
{ parametros par_program;
  char	     op,ar_conf[80];

  ini_arch();  // nombres de los archivos que recibirán los datos.
  if (argc==1) strcpy(ar_conf,"3.ini");
  else strcpy(ar_conf,argv[1]);
  config_ini(ar_conf,&par_program);   // configuracion inicial para poder calcular
  // funciOn para probar las transformacines de coordenadas:
  //pruebe(&par_program);
  // ciclo para secciones de Poincaré en órbitas de Langmuir:
  //Poinc_Langmuir(&par_program);
  // ciclo de ejecucion general del programa:
  if (par_program.batch) calcule_tray(&par_program);
  else
  do
  { op=par_program.opcion=menu(&par_program);
    switch(op)
    { case 'L':
	orb_Langmuir(&par_program,0,0,0);
	calcule_tray(&par_program);
      break;
      case '1':
      case '2':
      case '4':
      case '5':
      case '6':
      case '7':
	editar(&par_program);
      break;
      case '3':
	cout << "\nEntre el nuevo archivo con parámetros: ";
	cin  >> ar_conf;
	config_ini(ar_conf,&par_program);
	 // termina opcion de leer nuevas cond. de disco.
      break;
      case '8':
      case '9':
	calcule_tray(&par_program);
      break;
      case 'P':
	par_program.pflag=!par_program.pflag;
      break;
    } // fin del switch(op)
  }while (op!='S');
} // TERMINA main()
// ******************************
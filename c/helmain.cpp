/* Programa que calcula orbitas clasicas para el atomo de Helio, usando
   coordenadas regularizadas. La aproximacion clasica es valida en tanto
   que se estudian niveles doblemente excitados del Helio.
   Programa realizado originalmente en FORTRAN en Yugoslavia, 1993. */

#include <stdlib.h>
#include <math.h>
#include <conio.h>
#include <fstream.h>
#include "heldefs.h"
#include <iostream.h>
#include <stdio.h>
#include <string.h>


// Nombres DOS de los archivos utilizados:
// Directorio para datos:
caminoDOS dirdatos="";
/* Coordenadas de posicion de c/part (salida); Energia total del sistema, supuesta
   constante; Coordenadas relativas; Mapas de Poincaré: */
caminoDOS arch_pos[NPART+1], arch_etot, arch_rel[NPART+1], arch_poinc[NEQ/2+1];

// variables globales:
double  etot;
int     Z;


// Comienzan las funciones (la rutina main() esta al final):

// inicializa los nombres de los archivos que se van a usar:
void ini_arch()
{ int i;

  // Coordenadas de posicion (salida):
  for (i=1;i<=NPART;i++)
  { sprintf(arch_pos[i],"%spos%d.dat",dirdatos,i);
    sprintf(arch_rel[i],"%srel%d.dat",dirdatos,i);
  }
  // Energia total del sistema, supuesta constante (salida):
  sprintf(arch_etot,"%setot.dat",dirdatos);
  // Archivos que reciben la informacion de los mapas de Poincare:
  for (i=1;i<=NEQ/2;i++)
    sprintf(arch_poinc[i],"%spoinc%d.dat",dirdatos,i);
} // termina iniarch()


// rutina de calculo del programa:
void calcule_tray(parametros *par)
{ // las x son las coordenadas 'normales', y q son las regularizadas:
  double  x[DIM],q[DIM],paso,pasomin;
  double  eini,efin,deltaE,deltaE_r,tini;
  int     i,nbien,nmal;

  // codigo de la fn calcule():

  // colocamos las variables locales en sus valores para calcular:
  if (par->opcion=='8')
    for(i=0;i<=NEQ;i++) x[i]=par->cini[i];
  else
    for(i=0;i<=NEQ;i++) x[i]=par->cact[i];
  paso=0.01;
  pasomin=0.0;
  tini=0;
  energia(x,&eini);
  etot = eini;
  cout << "\n\n*************************************\n";
  cout << "Pulse cualquier tecla para suspender.\n";

  cout << "\nEnergia Total inicial = " << eini << "\n";

  // pasamos a coordenadas regularizadas antes de integrar:
  xtoq(x,q);

/*  A continuacion, se llama la rutina de integracion.
    Elegir entre R-K y B-S, y eliminar como comentario la que no se
    vaya a emplear:							*/

  // Con metodo de Burlich-Stoer:
  //odeint(q,NEQ,tini,tmax,inte_imag,eps,paso,pasomin,&nbien,&nmal,derivadas,bsstep);

  // Con metodo de Runge-Kutta:
  odeint_poinc(q,NEQ,tini,par->tmax,par->inte_imag,par->eps,paso,pasomin,
	       &nbien,&nmal,derivs_norm,rkqc,
	       par->pflag,par->psign_mom,par->pnue);
  qtox(q,x);
  energia(x,&efin);
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
  for(i=0;i<=NEQ;i++) par->cact[i]=x[i];

  // esperar para permitir la lectura de los parametros de energia:
//  beep(1000,0.2);
  cout << "\nPulse cualquier tecla para volver al menú.";
  getch();
} // termina calcule_tray()


// presenta el menu ppal. del  programa:
char menu(parametros *par)
{ char   selec;
  char   opciones[] = "0123456789SP";
  char   d[]        = "                  ";
  double *xini,*xact;

  xini=&par->cini[0];
  xact=&par->cact[0];

  clrscr();
  cout << d << "********************************************\n";
  cout << d << "* CALCULOS CLASICOS PARA EL ATOMO DE HELIO *\n";
  cout << d << "********************************************\n\n";
  cout << "MENU:\n\n";
  cout << "  1. Editar las condiciones iniciales:\n";
  cout << "     x1ini = "<<xini[1]<<"  y1ini = "<<xini[2]<<
	  "  Px1ini = "<<xini[3]<<"  Py1ini = "<<xini[4]<<"\n";
  cout << "     x2ini = "<<xini[5]<<"  y2ini = "<<xini[6]<<
	  "  Px2ini = "<<xini[7]<<"  Py2ini = "<<xini[8]<<"\n";
  cout << "  2. Editar las condiciones actuales del sistema:\n";
  cout << "     x1 = "<<xact[1]<<"  y1 = "<<xact[2]<<
	  "  Px1 = "<<xact[3]<<"  Py1 = "<<xact[4]<<"\n";
  cout << "     x2 = "<<xact[5]<<"  y2 = "<<xact[6]<<
	  "  Px2 = "<<xact[7]<<"  Py2 = "<<xact[8]<<"\n";
  cout << "  3. Leer configuración de disco.\n";
  cout << "  4. Cambiar tmax. Actualmente, tmax = "<<par->tmax<<"\n";
  cout << "  5. Cambiar epsilon. Actualmente, epsilon = "<<par->eps<<"\n";
  cout << "  6. Cambiar Z. Actualmente, Z = "<<par->Z<<"\n";
  cout << "  7. Cambiar el intervalo entre imágenes. Actualmente, inte_imag = "
       << par->inte_imag<<"\n";
  cout << "  8. RECALCULAR DESDE LAS COND. INICIALES.\n";
  cout << "  9. PROSEGUIR  DESDE LAS COND. ACTUALES.\n";
  cout << "  0. Calcular con modelo de Langmuir.\n";
  cout << " <P>  Calcular Poincare: " << par->pflag << "\n";
  cout << " <S> Salir del programa.\n\n";
  cout << "Elija una opción: ";
  selec=leatec(opciones);
  cout << selec;
  return selec;
}  // termina menu()


// muestra los valores de todos los parametros vigentes:
void muestre_param(parametros *par)
{ double *xini,*xact;

  xini=&par->cini[0];
  xact=&par->cact[0];

  cout << "La configuración del programa es en este momento:\n\n";
  cout << " Condiciones iniciales:\n";
  cout <<"   (1) x1 = "<<xini[1]<<"    (2) y1 = "<<xini[2]<<"    (3) Px1 = "<<xini[3]<<"    (4) Py1 = "<<xini[4];
  cout <<"\n   (5) x2 = "<<xini[5]<<"    (6) y2 = "<<xini[6]<<"    (7) Px2 = "<<xini[7]<<"    (8) Py2 = "<<xini[8];
  cout <<"\n\n Condiciones actuales:\n";
  cout <<"   (1) x1 = "<<xact[1]<<"    (2) y1 = "<<xact[2]<<"    (3) Px1 = "<<xact[3]<<"    (4) Py1 = "<<xact[4];
  cout <<"\n   (5) x2 = "<<xact[5]<<"    (6) y2 = "<<xact[6]<<"    (7) Px2 = "<<xact[7]<<"    (8) Py2 = "<<xact[8];
  cout <<"\n\n Otros:\n";
  cout <<"  tmax = "<<par->tmax<<";   epsilon = "<<par->eps<<";   Z = "<< par->Z <<
	 ";   Intervalo entre im genes = "<<par->inte_imag<<".\n";
}  // termina muestre_param()


// rutina de edicion de parametros:
void editar(parametros *par)
{ int   op,op_cond,ord_eps,salir;
  char  act[]         = "\nActualmente, ";
  char  nuev[]        = "Entre el nuevo valor: ";
  char  cond_permit[] = "012345678";

  salir=0;
  op=par->opcion;
  do
  { clrscr();                                 
    cout <<"EDICION DE PARAMETROS:\n";
    cout <<"----------------------\n\n";
    muestre_param(par);
    cout << "\n			 *********************\n";
    if ((op=='1')||(op=='2'))
    { cout << "\nEdición de condiciones ";
      if (op=='1')  cout << "INICIALES. ";
	else cout << "ACTUALES. ";
      cout << "\nPulse <0> para abandonar la edición.\n";
      cout << "Elija (1-8) qué condición desea alterar: ";
      op_cond=leatec(cond_permit)-48;    // ajuste con el codigo ASCII
      cout<<op_cond <<"\n";
    }
    if (op>50) salir=1;
      else if (op_cond==0) return;
    switch(op)
    { case '1':
	cout << "La cond. inicial elegida vale: " << par->cini[op_cond] <<"\n";
	cout << nuev;
	cin >> par->cini[op_cond];
      break;
      case '2':
	cout << "La cond. actual elegida vale: " << par->cact[op_cond] <<"\n";
	cout << nuev;
	cin >> par->cact[op_cond];
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
	cout << act << "Z = " << par->Z <<"\n";
	cout << nuev;
	cin >> par->Z;
	// Z es una vble global, que debe mantenerse actualizada:
	Z=par->Z;
      break;
      case '7':
	cout << act << "el intervalo entre imágenes es: " << par->inte_imag <<"\n";
	cout << nuev;
	cin >> par->inte_imag;
      break;
    }  // fin del switch
  } // fin del do
  while(!salir);
} // termina editar()


// rutina que lee un archivo de configuracion inicial de disco:
void config_ini(char *arch,parametros *par_ini)
{ ifstream   a_config(arch);
  double     x[DIM];
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
  // se leen los valores iniciales para X1,Y1,Px1,Py1,X2,Y2,Px2,Py2:
  for (i=1;i<=NEQ;i++) a_config >> x[i];
  // se leen los demas parametros:
  a_config >> par_ini->tmax >> par_ini->eps >> par_ini->Z >> par_ini->inte_imag;
  // cerramos el archivo:
  a_config.close();

  // se actualiza la estructura de parametros:
  for(i=1;i<=NEQ;i++) par_ini->cini[i]=par_ini->cact[i]=x[i];
  par_ini->opcion=8;

  // Z es una vble global, que debemos inicializar correctamente:
  Z=par_ini->Z;

  // ******************
  // ******************
  // BORRAR LUEGO:
  // INICIALIZACION DE LOS PARAMETROS DE POINCARE
  par_ini->pflag=0;
  par_ini->pcoord=1;
  par_ini->pval=0.0;
  par_ini->pmom=3;
  par_ini->psign_mom=1;
  par_ini->pnue=0;
  // ******************
  // ******************

}  // termina config_ini()

// Seleccion de parametros para calcular con el modelo de Langmuir:
void modelo_Lang(parametros *par_as)
{ double  x1,x2,p1y,p2y,easinc,aa,bb;
  int     i;

  clrscr();
  cout << "\n\n\nENTRADA DE DATOS PARA EL MODELO ASINCRONO:\n";
  cout << "------------------------------------------\n\n";
  cout << "Valores de x1 y x2: ";
  cin  >> x1 >> x2;
  aa = Z*(1/fabs(x1)+1/fabs(x2)) - 1/fabs(x1-x2);
  cout << "\nEnergía total del sistema.\n";
  do
  { cout << "Para evitar velocidades imaginarias, el valor de la energía debe\n";
    cout << "ser mayor o igual que: " << -aa << "\n";
    cout << "Adem s, recuerde que en un sistema ligado, Etot<0.\n\n";
    cout << "Entre Etot: ";
    cin  >> easinc;
    cout << "\n";
  } while ( (easinc<-aa) || (easinc>0) );
  bb = 1 + (x1*x1)/(x2*x2);
  p1y = sqrt(2*(easinc+aa)/bb);
  p2y = -p1y*x1/x2;

  // se actualiza la estructura de parametros:

  // y1 = y2 = p1x = p2x = 0:
  par_as->cini[2]=par_as->cini[3]=par_as->cini[6]=par_as->cini[7]=0;

  // asignamos x1,x2,py1,py2:
  par_as->cini[1]=x1;
  par_as->cini[4]=p1y;
  par_as->cini[5]=x2;
  par_as->cini[8]=p2y;

  // los parametros "actuales" se igualan a las cond. iniciales:
  for(i=1;i<=NEQ;i++) par_as->cact[i]=par_as->cini[i];
  // se indica recalcular:
  par_as->opcion=8;
} // termina modelo_Lang()



// ******************************
// RUTINA PRINCIPAL DEL PROGRAMA:
void main(void)
{ parametros par_program;
  char	     op,arch_cond_ini[80];

  ini_arch();  // nombres de los archivos que recibirán los datos.
  config_ini("helioreg.ini",&par_program);   // configuracion inicial para poder calcular

  // ciclo de ejecucion general del programa:
  do
  { op=par_program.opcion=menu(&par_program);
    switch(op)
    { case '0':
	modelo_Lang(&par_program);
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
	cin  >> arch_cond_ini;
	config_ini(arch_cond_ini,&par_program);
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
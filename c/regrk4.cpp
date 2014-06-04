/* Programa que calcula orbitas clasicas para el atomo de Helio, usando
   coordenadas regularizadas. La aproximacion clasica es valida en tanto
   que se estudian niveles doblemente excitados del Helio.
   Programa realizado originalmente en FORTRAN en Yugoslavia, 1993. */

#include <stdlib.h>
#include <stdio.h>
#include <conio.h>
#include <fstream.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
//#include <dos.h>

#define NEQ 8
#define DIM NEQ+1

// Nombres DOS de los archivos utilizados:

// Coordenadas de posicion (salida):
char arch_reg_A[]  = "c:\\tesis\\c\\regA.dat";
char arch_reg_B[]  = "c:\\tesis\\c\\regB.dat";

// Energia total del sistema, supuesta constante (salida):
char arch_etot[]   = "c:\\tesis\\c\\etot.dat";

// Coordenadas relativas regularizadas (salida):
char arch_rr[]     = "c:\\tesis\\c\\rr.dat";

// Archivos que reciben la informacion de los mapas de Poincare:
char arch_poinc1[] = "c:\\tesis\\c\\poinc1.dat";
char arch_poinc2[] = "c:\\tesis\\c\\poinc2.dat";


// se define un tipo que contendra todos los parametros de calculo:
typedef
  struct
   { double cini[DIM];
     double cact[DIM];
     float  tmax;
     float  paso;
     int    Z;
     int    inte_imag;
     char   opcion;
   } parametros;


// variables globales:
float e;
int   Z;


// Comienzan las funciones (la rutina main() esta al final):

// "pito" para usos varios (frec en Hz, tiempo en seg):
/*void beep(unsigned frec,float tiempo)
{ sound(frec);
  delay((unsigned) 1000*tiempo);
  nosound();
}  // termina beep()
*/

// rutina que pasa a coordenadas regularizadas:
void xtoq(double *xx,double *qq)
{ if (xx[1]>=0)
   { qq[1] = sqrt((sqrt(xx[1]* xx[1] + xx[2]* xx[2] ) + xx[1])/2);
     qq[2] = xx[2]/(2* qq[1]);
   }
  else
   { qq[2] = sqrt((sqrt(xx[1]* xx[1] + xx[2]* xx[2] ) - xx[1])/2);
     qq[1] = xx[2]/(2* qq[2]);
   }
  qq[3] = 2*(qq[1]* xx[3] + qq[2]* xx[4]);
  qq[4] = 2*(qq[1]* xx[4] - qq[2]* xx[3]);

  if (xx[5]>=0)
   { qq[5] = sqrt((sqrt(xx[5]* xx[5] + xx[6]* xx[6] ) + xx[5])/2);
     qq[6] = xx[6]/(2* qq[5]);
   }
  else
   { qq[6] = sqrt((sqrt(xx[5]* xx[5] + xx[6]* xx[6] ) - xx[5])/2);
     qq[5] = xx[6]/(2* qq[6]);
   }
  qq[7] = 2*(qq[5]* xx[7] + qq[6]* xx[8]);
  qq[8] = 2*(qq[5]* xx[8] - qq[6]* xx[7]);
} // termina xtoq()


// transformacion inversa de xtoq():
void qtox(double *qq,double *xx)
{ xx[1] = qq[1]* qq[1] - qq[2]* qq[2];
  xx[2] = 2* qq[1] * qq[2];
  xx[3] = (qq[1]* qq[3]  - qq[2]* qq[4])/(2*(qq[1]* qq[1]+qq[2]* qq[2]));
  xx[4] = (qq[2]* qq[3]  + qq[1]* qq[4])/(2*(qq[1]* qq[1]+qq[2]* qq[2]));
  xx[5] = qq[5]* qq[5] - qq[6]* qq[6];
  xx[6] = 2* qq[5] * qq[6];
  xx[7] = (qq[5]* qq[7]  - qq[6]* qq[8])/(2*(qq[5]* qq[5]+qq[6]* qq[6]));
  xx[8] = (qq[6]* qq[7]  + qq[5]* qq[8])/(2*(qq[5]* qq[5]+qq[6]* qq[6]));
} // termina qtox()


// funcion usada en el metodo R-K:
void fun(double *xx,double *ff)
{ double rr1,rr2,aa,bb,rr12,pp1,pp2;

  rr1 = xx[1]*xx[1] + xx[2]*xx[2];
  rr2 = xx[5]*xx[5] + xx[6]*xx[6];
  aa = xx[1]*xx[1] - xx[2]*xx[2] - xx[5]*xx[5] + xx[6]*xx[6];
  bb = xx[1]*xx[2] - xx[5]*xx[6];
  rr12 = sqrt(aa*aa+4*bb*bb);
  pp1 = xx[3]*xx[3] + xx[4]*xx[4];
  pp2 = xx[7]*xx[7] + xx[8]*xx[8];

  ff[1] = rr2*xx[3]/4;
  ff[2] = rr2*xx[4]/4;
  ff[3] = -2*(pp2/8 - Z + rr2*(1/rr12-e))*xx[1] +
	  2*rr1*rr2*(aa*xx[1] + 2*bb*xx[2])/rr12/rr12/rr12;
  ff[4] = -2*(pp2/8 - Z + rr2*(1/rr12-e))*xx[2] +
	  2*rr1*rr2*(2*bb*xx[1] - aa*xx[2])/rr12/rr12/rr12;
  ff[5] = rr1*xx[7]/4;
  ff[6] = rr1*xx[8]/4;
  ff[7] = -2*(pp1/8 - Z + rr1*(1/rr12-e))*xx[5] -
	  2*rr1*rr2*(aa*xx[5] + 2*bb*xx[6])/rr12/rr12/rr12;
  ff[8] = -2*(pp1/8 - Z + rr1*(1/rr12-e))*xx[6] -
	  2*rr1*rr2*(2*bb*xx[5] - aa*xx[6])/rr12/rr12/rr12;
}   // termina fun()


// metodo de Runge-Kutta de 4ø orden:
void rk4(double *xx, double *step)
{ double   yy[DIM],ff[DIM];
  float    k1[DIM],k2[DIM],k3[DIM],k4[DIM];
  int      j;

  fun(xx,ff);
  for(j=1;j<=NEQ;j++)
    { k1[j]=*step * ff[j];
      yy[j]=xx[j]+k1[j]/2;}
  fun(yy,ff);
  for(j=1;j<=NEQ;j++)
    { k2[j]=*step * ff[j];
      yy[j]=xx[j]+k2[j]/2;}
  fun(yy,ff);
  for(j=1;j<=NEQ;j++)
    { k3[j]=*step * ff[j];
      yy[j]=xx[j]+k3[j];}
  fun(yy,ff);
  for(j=1;j<=NEQ;j++)
    { k4[j]=*step * ff[j];
      xx[j]=xx[j]+(k1[j]+2*k2[j]+2*k3[j]+k4[j])/6;}
}   // termina rk4()


// funcion que calcula la energia total del sistema:
void energia(double *xx,float *etot)
{ double rr1,rr2,rr12,pp1,pp2;

  rr1 = sqrt(xx[1]*xx[1] + xx[2]*xx[2]);
  rr2 = sqrt(xx[5]*xx[5] + xx[6]*xx[6]);
  rr12 = sqrt((xx[1]-xx[5])*(xx[1]-xx[5]) + (xx[2]-xx[6])*(xx[2]-xx[6]));
  pp1 = xx[3]*xx[3] + xx[4]*xx[4];
  pp2 = xx[7]*xx[7] + xx[8]*xx[8];
  *etot = pp1/2 + pp2/2 - Z/rr1 - Z/rr2 + 1/rr12;
}  // termina energia()


// rutina de calculo del programa:
void calcule(parametros *par_calc)
{ // las x son las coordenadas 'normales', y q son las regularizadas:
  double  x[DIM],q[DIM],r1,r2,dt,rpaso,paso;
  float   etot,eini,efin,deltaE,deltaE_r,t,tmax;
  int     i,inte_imag,cont_imag;
  char    tecla;

  // codigo de la fn calcule():

  // se declaran los streams y se abren los archivos necesarios:
  ofstream a_reg_A(arch_reg_A);
  ofstream a_reg_B(arch_reg_B);
  ofstream a_etot(arch_etot);
  ofstream a_rr(arch_rr);
//  ofstream a_poinc1(arch_poinc1);
//  ofstream a_poinc2(arch_poinc2);

  // colocamos las variables locales en sus valores para calcular:
  if (par_calc->opcion=='8')
    for(i=0;i<=NEQ;i++) x[i]=par_calc->cini[i];
  else
    for(i=0;i<=NEQ;i++) x[i]=par_calc->cact[i];
  paso=par_calc->paso;
  tmax=par_calc->tmax;
  inte_imag=par_calc->inte_imag;

  energia(x,&etot);
  e = eini = etot;
  cout << "\n\nEnergia Total inicial = " << eini << "\n";
  t=0;
  xtoq(x,q);
  cont_imag=0;     // contador de 'imagenes'.
  cout << "\nPulse cualquier tecla para suspender.";
  cout << "\nCalculando.";

  // ciclo de calculo:
  do
  { r1=q[1]*q[1] + q[2]*q[2];
    r2=q[5]*q[5] + q[6]*q[6];

    // se graba solo cada cierto numero de imagenes (dado por inte_imag):
    if (cont_imag%inte_imag==0)
    { a_reg_A << x[1] << "   " << x[2] << "\n";
      a_reg_B << x[5] << "   " << x[6] << "\n";
      a_etot << t << "   " << etot << "\n";
      a_rr << r1 << "   " << r2 << "\n";
    }
    dt = r1*r2*paso;
    t += dt;
    rk4(q,&paso);
    qtox(q,x);
    energia(x,&etot);
    if (cont_imag%200==0) cout << ".";
    cont_imag++;
  } while ( (t<tmax) && (!kbhit()) );

  // si se puls¢ una tecla, la 'atrapamos', e indicamos el tiempo:
  if (kbhit())
   { tecla=getch();
     cout << "\n\nCALCULO INTERRUMPIDO EN t = " << t;
     cout << ".  (Se ten¡a tmax = " << tmax << ").";
   }

  a_reg_A << x[1] << "   " << x[2] << "\n";
  a_reg_B << x[5] << "   " << x[6] << "\n";

  energia(x,&etot);
  a_etot << t << "   " << etot << "\n";
  efin = etot;
  cout << "\n\nEnerg¡a Total final   = " << efin <<"\n\n";

  // calculamos el cambio de energia inducido por las rutinas numericas:
  deltaE=fabs(efin-eini);
  deltaE_r=100*deltaE/fabs(eini);
  cout << "|Delta Energ¡a|       = " << deltaE << "\n";
  cout << "|Delta Energ¡a|(rel)  = " << deltaE_r << " %\n";

  // cerramos los archivos de datos:
  a_reg_A.close();
  a_reg_B.close();
  a_etot.close();
  a_rr.close();
//  a_poinc1.close();
//  a_poinc2.close();

  // actualizamos la estructura de parametros con los ultimos valores
  // de las distintas vbles dinamicas del problema:
  for(i=0;i<=NEQ;i++) par_calc->cact[i]=x[i];

  // esperar para permitir la lectura de los parametros de energia:
//  beep(1000,0.2);
  cout << "\n\nPulse cualquier tecla para volver al men£.";
  getch();

} // termina calcule()


// rutina que lee una tecla solo dentro de un cierto conjunto de teclas
// permitidas (devuelve mayusculas). No hace eco a pantalla:
char leatec(const char *opciones)
{ char selec;

  do
  { selec=toupper(getch());
    if (!strchr(opciones,selec))/* beep(1000,0.05)*/;
  } while(!strchr(opciones,selec));
  return selec;
} // termina leatec()


// presenta un el menu ppal. del  programa:
char menu(parametros *par)
{ char   selec;
  char   d[]        = "       ";
  char   opciones[] = "0123456789S";
  double *xini,*xact;

  xini=&par->cini[0];
  xact=&par->cact[0];

  clrscr();
  cout << d << "**************************************************************\n";
  cout << d << "CALCULOS CLASICOS PARA EL PROBLEMA ELECTROSTATICO DE 3 CUERPOS\n";
  cout << d << "**************************************************************\n\n";

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
  cout << "  3. Leer configuraci¢n de disco.\n";
  cout << "  4. Cambiar tmax. Actualmente, tmax = "<<par->tmax<<"\n";
  cout << "  5. Cambiar paso. Actualmente, paso = "<<par->paso<<"\n";
  cout << "  6. Cambiar Z. Actualmente, Z = "<<par->Z<<"\n";
  cout << "  7. Cambiar el intervalo entre im genes. Actualmente, inte_imag = "
       << par->inte_imag<<"\n";
  cout << "  8. RECALCULAR DESDE LAS COND. INICIALES.\n";
  cout << "  9. PROSEGUIR  DESDE LAS COND. ACTUALES.\n";
  cout << "  0. Calcular con modelo ASINCRONO.\n\n";
  cout << " <s> Salir del programa.\n\n";
  cout << "Elija una opci¢n: ";
  selec=leatec(opciones);
  cout << selec;
  return selec;
}  // termina menu()


// muestra los valores de todos los parametros vigentes:
void muestre_param(parametros *par)
{ double *xini,*xact;

  xini=&par->cini[0];
  xact=&par->cact[0];

  cout << "La configuraci¢n del programa es en este momento:\n\n";
  cout << " Condiciones iniciales:\n";
  cout <<"  x1 = "<<xini[1]<<"  y1 = "<<xini[2]<<"  Px1 = "<<xini[3]<<"  Py1 = "<<xini[4];
  cout <<"\n  x2 = "<<xini[5]<<"  y2 = "<<xini[6]<<"  Px2 = "<<xini[7]<<"  Py2 = "<<xini[8];
  cout <<"\n\n Condiciones actuales:\n";
  cout <<"  x1 = "<<xact[1]<<"  y1 = "<<xact[2]<<"  Px1 = "<<xact[3]<<"  Py1 = "<<xact[4];
  cout <<"\n  x2 = "<<xact[5]<<"  y2 = "<<xact[6]<<"  Px2 = "<<xact[7]<<"  Py2 = "<<xact[8];
  cout <<"\n\n Otros:\n";
  cout <<"  tmax = "<<par->tmax<<";   paso = "<<par->paso<<";   Z = "<< par->Z <<
	 ";   Intervalo entre im genes = "<<par->inte_imag<<".\n";
}  // termina muestre_param()


// rutina de edicion de parametros:
void editar(parametros *par)
{ int   op,op_cond;
  char  act[]         = "\nActualmente, ";
  char  nuev[]        = "Entre el nuevo valor: ";
  char  cond_permit[] = "12345678";

  op=par->opcion;
  clrscr();
  cout <<"EDICION DE PARAMETROS:\n";
  cout <<"----------------------\n\n";
  muestre_param(par);
  cout << "\n			 *********************\n";
  if ((op=='1')||(op=='2'))
  { cout << "\nEdici¢n de condiciones ";
    if (op=='1')  cout << "iniciales.\n";
      else cout << "actuales.\n";
    cout << "Elija (1-8) qu‚ condici¢n desea alterar: ";
    op_cond=leatec(cond_permit)-48;
    cout<<op_cond <<"\n";
  }
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
      cout << act << "paso = " << par->paso <<"\n";
      cout << nuev;
      cin >> par->paso;
    break;
    case '6':
      cout << act << "Z = " << par->Z <<"\n";
      cout << nuev;
      cin >> par->Z;
      // Z es una vble global, que debe mantenerse actualizada:
      Z=par->Z;
    break;
    case '7':
      cout << act << "el intervalo entre im genes es: " << par->inte_imag <<"\n";
      cout << nuev;
      cin >> par->inte_imag;
    break;
  }
} // termina editar()


// rutina que lee un archivo de configuracion inicial de disco:
void config_ini(char *arch,parametros *par_ini)
{ ifstream   a_config(arch);
  double     x[DIM];
  int        i;

  if (!a_config)
  { cout << "\n*********** ERROR ***********\n";
//    beep(1000,0.75);
    cout << "No es posible leer el archivo de configuraci¢n.\n";
    cout << "No se carg¢ ninguna variable.\n";
    cout << "Si el programa se ejecuta por primera vez, config£relo manualmente,\n";
    cout << "pues puede haber \"basura\" en los distintos par metros.\n";
    cout << "Si no est  comenzando la ejecuci¢n, simplemente se mantendr  la\n";
    cout << "configuraci¢n que estuviese vigente.\n\n";
    cout << "Pulse cualquier tecla para seguir.";
    getch();
    return;
  };
  // se leen los valores iniciales para X1,Y1,Px1,Py1,X2,Y2,Px2,Py2:
  a_config >> x[1] >> x[2] >> x[3] >> x[4] >> x[5] >> x[6] >> x[7] >> x[8];
  // se leen los demas parametros:
  a_config >> par_ini->tmax >> par_ini->paso >> par_ini->Z >>
	      par_ini->inte_imag;
  // cerramos el archivo:
  a_config.close();

  // se actualiza la estructura de parametros:
  for(i=1;i<=NEQ;i++) par_ini->cini[i]=par_ini->cact[i]=x[i];
  par_ini->opcion=8;

  // Z es una vble global, que debemos inicializar correctamente:
  Z=par_ini->Z;
}  // termina config_ini()


void modelo_asinc(parametros *par_as)
{ double  x1,x2,p1y,p2y,easinc,aa,bb;
  int     i;

  clrscr();
  cout << "\n\n\nENTRADA DE DATOS PARA EL MODELO ASINCRONO:\n";
  cout << "------------------------------------------\n\n";
  cout << "Valores de x1 y x2: ";
  cin  >> x1 >> x2;
  aa = Z*(1/fabs(x1)+1/fabs(x2)) - 1/fabs(x1-x2);
  cout << "\nEnerg¡a total del sistema.\n";
  do
  { cout << "Para evitar velocidades imaginarias, el valor de la energ¡a debe\n";
    cout << "ser mayor o igual que: " << -aa << "\n";
    cout << "Adem s, recuerde que en un sistema ligado, Etot<0.\n\n";
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
} // termina modelo_asinc()



// ******************************
// RUTINA PRINCIPAL DEL PROGRAMA:
void main(void)
{ parametros par_program;
  char	     op,arch_cond_ini[80];

  // configuracion inicial para poder calcular:
  config_ini("c:\\tesis\\c\\regrk4.ini",&par_program);

  // ciclo de ejecucion general del programa:
  do
  { op=par_program.opcion=menu(&par_program);
    if (op=='0')
     { modelo_asinc(&par_program);
       calcule(&par_program);
     }
    else
      if ((op=='1')||(op=='2')||((op>='4')&&(op<='7')))
	editar(&par_program);
    else
      if (op=='3')
       { cout << "\nEntre el nuevo archivo con par metros: ";
	 cin  >> arch_cond_ini;
	 config_ini(arch_cond_ini,&par_program);
       }  // termina opcion de leer nuevas cond. de disco.
    else
      if ((op=='8')||(op=='9'))
	calcule(&par_program);
  }while (op!='S');

} // TERMINA main()
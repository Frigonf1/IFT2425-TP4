//------------------------------------------------------
// module  :  Tp4-IFT2425-1.c
// author  :  François Frigon - 20297551 - francois.frigon@umontreal.ca
//            Xavier Dontigny - 20215658 - xavier.dontigny@umontreal.ca
// date    :  20 avril 2025
// version :  1.0
// language:  C++
// note    :
//------------------------------------------------------
//  

//------------------------------------------------
// FICHIERS INCLUS -------------------------------
//------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <new>
#include <unistd.h>

/************************************************************************/
/* WINDOWS						          	*/
/************************************************************************/
#include <X11/Xutil.h>

Display   *display;
int	  screen_num;
int 	  depth;
Window	  root;
Visual*	  visual;
GC	  gc;

//------------------------------------------------
// DEFINITIONS -----------------------------------                       
//------------------------------------------------
#define CARRE(X) ((X)*(X))

#define OUTPUT_FILE "Tp4-Img-I.pgm"
#define VIEW_PGM    "xv" 
#define DEBUG 0

//-Cst-Modele
#define X_1 0.0
#define Y_1 1.0
#define X_2 -1.0/sqrt(2.0)
#define Y_2 -1.0/2.0
#define X_3 +1.0/2*sqrt(2.0)
#define Y_3 -1.0/2.0
#define C 0.25 
#define R 0.1
#define D 0.3

#define X_1_INI 0.2            
#define X_2_INI 0.0       
#define X_3_INI -1.6                  
#define X_4_INI 0.0 
 
//-Cst-Runge-Kutta
#define H            0.0001
#define T_0          0.0                 
#define T_F         30.0 
#define NB_INTERV (T_F-T_0)/H 
   
 //-Cst-Image                             
#define WIDTH  512
#define HEIGHT 512                  
#define MAX_X  4.0                
#define MAX_Y  4.0  
#define EVOL_GRAPH 3000
              
#define WHITE     255
#define GREYWHITE 230
#define GREY      200
#define GREYDARK  120
#define BLACK       0   

//------------------------------------------------
// GLOBAL CST ------------------------------------                       
//------------------------------------------------
float Xmin=0.0;
float Xmax=0.0;
float Ymin=0.0;
float Ymax=0.0; 

float xx_1=((WIDTH/MAX_X)*X_1)+(WIDTH/2);
float yy_1=(-(HEIGHT/MAX_Y)*Y_1)+(HEIGHT/2);
float xx_2=((WIDTH/MAX_X)*X_2)+(WIDTH/2);
float yy_2=(-(HEIGHT/MAX_Y)*Y_2)+(HEIGHT/2);
float xx_3=((WIDTH/MAX_X)*X_3)+(WIDTH/2);
float yy_3=(-(HEIGHT/MAX_Y)*Y_3)+(HEIGHT/2);

/************************************************************************/
/* OPEN_DISPLAY()							*/
/************************************************************************/
int open_display()
{
  if ((display=XOpenDisplay(NULL))==NULL)
   { printf("Connection impossible\n");
     return(-1); }

  else
   { screen_num=DefaultScreen(display);
     visual=DefaultVisual(display,screen_num);
     depth=DefaultDepth(display,screen_num);
     root=RootWindow(display,screen_num);
     return 0; }
}

/************************************************************************/
/* FABRIQUE_WINDOW()							*/
/* Cette fonction crée une fenetre X et l'affiche à l'écran.	        */
/************************************************************************/
Window fabrique_window(char *nom_fen,int x,int y,int width,int height,int zoom)
{
  Window                 win;
  XSizeHints      size_hints;
  XWMHints          wm_hints;
  XClassHint     class_hints;
  XTextProperty  windowName, iconName;

  char *name=nom_fen;

  if(zoom<0) { width/=-zoom; height/=-zoom; }
  if(zoom>0) { width*=zoom;  height*=zoom;  }

  win=XCreateSimpleWindow(display,root,x,y,width,height,1,0,255);

  size_hints.flags=PPosition|PSize|PMinSize;
  size_hints.min_width=width;
  size_hints.min_height=height;

  XStringListToTextProperty(&name,1,&windowName);
  XStringListToTextProperty(&name,1,&iconName);
  wm_hints.initial_state=NormalState;
  wm_hints.input=True;
  wm_hints.flags=StateHint|InputHint;
  class_hints.res_name=nom_fen;
  class_hints.res_class=nom_fen;

  XSetWMProperties(display,win,&windowName,&iconName,
                   NULL,0,&size_hints,&wm_hints,&class_hints);

  gc=XCreateGC(display,win,0,NULL);

  XSelectInput(display,win,ExposureMask|KeyPressMask|ButtonPressMask| 
               ButtonReleaseMask|ButtonMotionMask|PointerMotionHintMask| 
               StructureNotifyMask);

  XMapWindow(display,win);
  return(win);
}

/****************************************************************************/
/* CREE_XIMAGE()							    */
/* Crée une XImage à partir d'un tableau de float                          */
/* L'image peut subir un zoom.						    */
/****************************************************************************/
XImage* cree_Ximage(float** mat,int z,int length,int width)
{
  int lgth,wdth,lig,col,zoom_col,zoom_lig;
  float somme;
  unsigned char	 pix;
  unsigned char* dat;
  XImage* imageX;

  /*Zoom positiv*/
  /*------------*/
  if (z>0)
  {
   lgth=length*z;
   wdth=width*z;

   dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
   if (dat==NULL)
      { printf("Impossible d'allouer de la memoire.");
        exit(-1); }

  for(lig=0;lig<lgth;lig=lig+z) for(col=0;col<wdth;col=col+z)
   { 
    pix=(unsigned char)mat[lig/z][col/z];
    for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
      { 
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+0)]=pix;
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+1)]=pix;
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+2)]=pix;
       dat[((lig+zoom_lig)*wdth*4)+((4*(col+zoom_col))+3)]=pix; 
       }
    }
  } /*--------------------------------------------------------*/

  /*Zoom negatifv*/
  /*------------*/
  else
  {
   z=-z;
   lgth=(length/z);
   wdth=(width/z);

   dat=(unsigned char*)malloc(lgth*(wdth*4)*sizeof(unsigned char));
   if (dat==NULL)
      { printf("Impossible d'allouer de la memoire.");
        exit(-1); }

  for(lig=0;lig<(lgth*z);lig=lig+z) for(col=0;col<(wdth*z);col=col+z)
   {  
    somme=0.0;
    for(zoom_lig=0;zoom_lig<z;zoom_lig++) for(zoom_col=0;zoom_col<z;zoom_col++)
     somme+=mat[lig+zoom_lig][col+zoom_col];
           
     somme/=(z*z);    
     dat[((lig/z)*wdth*4)+((4*(col/z))+0)]=(unsigned char)somme;
     dat[((lig/z)*wdth*4)+((4*(col/z))+1)]=(unsigned char)somme;
     dat[((lig/z)*wdth*4)+((4*(col/z))+2)]=(unsigned char)somme;
     dat[((lig/z)*wdth*4)+((4*(col/z))+3)]=(unsigned char)somme; 
   }
  } /*--------------------------------------------------------*/

  imageX=XCreateImage(display,visual,depth,ZPixmap,0,(char*)dat,wdth,lgth,16,wdth*4);
  return (imageX);
}

//------------------------------------------------
// FUNCTIONS -------------------------------------                       
//------------------------------------------------
//-------------------------//
//-- Matrice de Double ----//
//-------------------------//
//---------------------------------------------------------
// Alloue de la memoire pour une matrice 1d de float
//----------------------------------------------------------
float* dmatrix_allocate_1d(int hsize)
 {
  float* matrix;
  matrix=new float[hsize]; return matrix; }

//----------------------------------------------------------
// Alloue de la memoire pour une matrice 2d de float
//----------------------------------------------------------
float** dmatrix_allocate_2d(int vsize,int hsize)
 {
  float** matrix;
  float *imptr;

  matrix=new float*[vsize];
  imptr=new float[(hsize)*(vsize)];
  for(int i=0;i<vsize;i++,imptr+=hsize) matrix[i]=imptr;
  return matrix;
 }

//----------------------------------------------------------
// Libere la memoire de la matrice 1d de float
//----------------------------------------------------------
void free_dmatrix_1d(float* pmat)
{ delete[] pmat; }

//----------------------------------------------------------
// Libere la memoire de la matrice 2d de float
//----------------------------------------------------------
void free_dmatrix_2d(float** pmat)
{ delete[] (pmat[0]);
  delete[] pmat;}

//----------------------------------------------------------
// SaveImagePgm                       
//----------------------------------------------------------                
void SaveImagePgm(char* name,float** mat,int lgth,int wdth)
{
 int i,j;
 char buff[300];
 FILE* fic;

  //--extension--
  strcpy(buff,name);

  //--ouverture fichier--
  fic=fopen(buff,"wb");
    if (fic==NULL) 
        { printf("Probleme dans la sauvegarde de %s",buff); 
          exit(-1); }
  printf("\n Sauvegarde de %s au format pgm\n",buff);

  //--sauvegarde de l'entete--
  fprintf(fic,"P5");
  fprintf(fic,"\n# IMG Module");
  fprintf(fic,"\n%d %d",wdth,lgth);
  fprintf(fic,"\n255\n");

  //--enregistrement--
  for(i=0;i<lgth;i++) for(j=0;j<wdth;j++) 
	fprintf(fic,"%c",(char)mat[i][j]);
   
  //--fermeture fichier--
  fclose(fic); 
}

//------------------------------------------------------------------------
// plot_point
//
// Affiche entre x dans [-MAX_X/2  MAX_X/2]
//               y dans [-MAX_Y/2  MAX_Y/2]                
//------------------------------------------------------------------------
void plot_point(float** MatPts,float** MatPict,int NbPts)
{
 int x_co,y_co;
 int i,j,k;

 //Init
 for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++)  MatPict[i][j]=GREYWHITE;

 for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) 
   { if ((fabs(i-yy_1)+fabs(j-xx_1))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_2)+fabs(j-xx_2))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_3)+fabs(j-xx_3))<10) MatPict[i][j]=GREYDARK; }

 //Loop
 for(k=0;k<NbPts;k++)
    { x_co=(int)((WIDTH/MAX_X)*MatPts[k][0]);
      y_co=-(int)((HEIGHT/MAX_Y)*MatPts[k][1]);
      y_co+=(HEIGHT/2);
      x_co+=(WIDTH/2);
      if (DEBUG) printf("[%d::%d]",x_co,y_co); 
      if ((x_co<WIDTH)&&(y_co<HEIGHT)&&(x_co>0)&&(y_co>0)) 
	 MatPict[y_co][x_co]=BLACK; 
    }
}

//------------------------------------------------------------------------
// Fill_Pict
//------------------------------------------------------------------------
void Fill_Pict(float** MatPts,float** MatPict,int PtsNumber,int NbPts)
{
 int i,j;
 int x_co,y_co;
 int k,k_Init,k_End;

 //Init
 for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) 
   { if (MatPict[i][j]!=GREYWHITE) MatPict[i][j]=GREY;
     if ((fabs(i-yy_1)+fabs(j-xx_1))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_2)+fabs(j-xx_2))<10) MatPict[i][j]=GREYDARK;
     if ((fabs(i-yy_3)+fabs(j-xx_3))<10) MatPict[i][j]=GREYDARK; }

 //Loop
 k_Init=PtsNumber;
 k_End=(k_Init+EVOL_GRAPH)%NbPts;
 for(k=k_Init;k<k_End;k++)
    { k=(k%NbPts);
      x_co=(int)((WIDTH/MAX_X)*MatPts[k][0]);
      y_co=-(int)((HEIGHT/MAX_Y)*MatPts[k][1]);
      y_co+=(HEIGHT/2);
      x_co+=(WIDTH/2);
      if ((x_co<WIDTH)&&(y_co<HEIGHT)&&(x_co>0)&&(y_co>0)) 
         MatPict[y_co][x_co]=BLACK; }
}


//------------------------------------------------
// FONCTIONS TPs----------------------------------                      
//------------------------------------------------

int N = 3; // Nombre de points

float Xi[3] = {X_1, X_2, X_3};
float Yi[3] = {Y_1, Y_2, Y_3};

float* equationDiff(float t, float x, float y, float dx, float dy){

    float fx = 0.0;
    float fy = 0.0;

    for (int i = 0; i < N; i++)
    {
      float numX = Xi[i] - x;
      float numY = Yi[i] - y;
      float denominateur = pow(sqrt(pow(numX, 2) + pow(numY, 2) + D * D),3);
      fx += numX / denominateur;
      fy += numY / denominateur;
    }

    static float result[2];

    float accX = -R * dx + fx - C * x;
    float accY = -R * dy + fy - C * y;

    result[0] = accX;
    result[1] = accY;

  return result;

}

void RungeKuttaFehlberg(float t, float x, float y, float dx, float dy)
{
  float k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y, k5x, k5y, k6x, k6y;

  for (int i = T_0; i < T_F; i += H)
  {
  
    k1x = H * equationDiff(t, x, y, dx, dy)[0];
    k1y = H * equationDiff(t, x, y, dx, dy)[1];

    k2x = H * equationDiff(t+H/4, x+k1x/4, y+k1y/4, dx, dy)[0];
    k2y = H * equationDiff(t+H/4, x+k1x/4, y+k1y/4, dx, dy)[1];

    k3x = H * equationDiff(t+3*H/8, x+3*k1x/8+9*k2x/32, y+3*k1y/8+9*k2y/32, dx, dy)[0];
    k3y = H * equationDiff(t+3*H/8, x+3*k1x/8+9*k2x/32, y+3*k1y/8+9*k2y/32, dx, dy)[1];

    k4x = H * equationDiff(t+12*H/13, x+1932*k1x/2197-7200*k2x/2197+7296*k3x/2197, y+1932*k1y/2197-7200*k2y/2197+7296*k3y/2197, dx, dy)[0];
    k4y = H * equationDiff(t+12*H/13, x+1932*k1x/2197-7200*k2x/2197+7296*k3x/2197, y+1932*k1y/2197-7200*k2y/2197+7296*k3y/2197, dx, dy)[1];

    k5x = H * equationDiff(t+H, x+439*k1x/216-8*k2x+3680*k3x/513-845*k4x/4104, y+439*k1y/216-8*k2y+3680*k3y/513-845*k4y/4104, dx, dy)[0];
    k5y = H * equationDiff(t+H, x+439*k1x/216-8*k2x+3680*k3x/513-845*k4x/4104, y+439*k1y/216-8*k2y+3680*k3y/513-845*k4y/4104, dx, dy)[1];

    k6x = H * equationDiff(t+H/2, x-8*k1x/27+2*k2x-3544*k3x/2565+1859*k4x/4104-11*k5x/40, y-8*k1y/27+2*k2y-3544*k3y/2565+1859*k4y/4104-11*k5y/40, dx, dy)[0];
    k6y = H * equationDiff(t+H/2, x-8*k1x/27+2*k2x-3544*k3x/2565+1859*k4x/4104-11*k5x/40, y-8*k1y/27+2*k2y-3544*k3y/2565+1859*k4y/4104-11*k5y/40, dx, dy)[1];

    dx += (16*k1x/135 + 6656*k3x/12825 + 28561*k4x/56430 - 9*k5x/50 + 2*k6x/55);
    dy += (16*k1y/135 + 6656*k3y/12825 + 28561*k4y/56430 - 9*k5y/50 + 2*k6y/55);

    x += dx * H;
    y += dy * H;
    t += H;
  }
}



//----------------------------------------------------------
//----------------------------------------------------------
// MAIN  
//----------------------------------------------------------
//----------------------------------------------------------
int main (int argc, char **argv)
{
  int i,j,k;
  int flag_graph;
  int zoom;

  XEvent ev;
  Window win_ppicture;
  XImage *x_ppicture;
  char   nomfen_ppicture[100]; 
  char BufSystVisu[100];

  //>AllocMemory
  float** MatPict=dmatrix_allocate_2d(HEIGHT,WIDTH);
  float** MatPts=dmatrix_allocate_2d((int)(NB_INTERV),2);
  
  //>Init
  for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) MatPict[i][j]=GREYWHITE;
  for(i=0;i<2;i++) for(j=0;j<(int)(NB_INTERV);j++) MatPts[i][j]=0.0;
  flag_graph=1;
  zoom=1;


  //---------------------------------------------------------------------
  //>Question 1 
  //---------------------------------------------------------------------  

  //Il faut travailler ici ...et dans > // FONCTIONS TPs

  float k1x, k1y, k2x, k2y, k3x, k3y, k4x, k4y, k5x, k5y, k6x, k6y;
  float x = X_1_INI;
  float y = 0.0;
  float dx = 0.0;
  float dy = 0.0;
  float t = T_0;

  for (int i = 0; i < NB_INTERV; i++)
  {
  
    k1x = H * equationDiff(i*H, x, y, dx, dy)[0];
    k1y = H * equationDiff(i*H, x, y, dx, dy)[1];

    k2x = H * equationDiff(i*H+H/4, x+k1x/4, y+k1y/4, dx, dy)[0];
    k2y = H * equationDiff(i*H+H/4, x+k1x/4, y+k1y/4, dx, dy)[1];

    k3x = H * equationDiff(i*H+3*H/8, x+3*k1x/8+9*k2x/32, y+3*k1y/8+9*k2y/32, dx, dy)[0];
    k3y = H * equationDiff(i*H+3*H/8, x+3*k1x/8+9*k2x/32, y+3*k1y/8+9*k2y/32, dx, dy)[1];

    k4x = H * equationDiff(i*H+12*H/13, x+1932*k1x/2197-7200*k2x/2197+7296*k3x/2197, y+1932*k1y/2197-7200*k2y/2197+7296*k3y/2197, dx, dy)[0];
    k4y = H * equationDiff(i*H+12*H/13, x+1932*k1x/2197-7200*k2x/2197+7296*k3x/2197, y+1932*k1y/2197-7200*k2y/2197+7296*k3y/2197, dx, dy)[1];

    k5x = H * equationDiff(i*H+H, x+439*k1x/216-8*k2x+3680*k3x/513-845*k4x/4104, y+439*k1y/216-8*k2y+3680*k3y/513-845*k4y/4104, dx, dy)[0];
    k5y = H * equationDiff(i*H+H, x+439*k1x/216-8*k2x+3680*k3x/513-845*k4x/4104, y+439*k1y/216-8*k2y+3680*k3y/513-845*k4y/4104, dx, dy)[1];

    k6x = H * equationDiff(i*H+H/2, x-8*k1x/27+2*k2x-3544*k3x/2565+1859*k4x/4104-11*k5x/40, y-8*k1y/27+2*k2y-3544*k3y/2565+1859*k4y/4104-11*k5y/40, dx, dy)[0];
    k6y = H * equationDiff(i*H+H/2, x-8*k1x/27+2*k2x-3544*k3x/2565+1859*k4x/4104-11*k5x/40, y-8*k1y/27+2*k2y-3544*k3y/2565+1859*k4y/4104-11*k5y/40, dx, dy)[1];

    dx += (16*k1x/135 + 6656*k3x/12825 + 28561*k4x/56430 - 9*k5x/50 + 2*k6x/55);
    dy += (16*k1y/135 + 6656*k3y/12825 + 28561*k4y/56430 - 9*k5y/50 + 2*k6y/55);

    x += dx * H;
    y += dy * H;
    t += H;
    printf("x: %f, y: %f\n", x, y);
    MatPts[i][0] = x;
    MatPts[i][1] = y;

  }
    

  //Un exemple ou la matrice de points est remplie
  //par une courbe donné par l'équation d'en bas... et non pas par 
  //la solution de l'équation différentielle
 
  //for(k=0;k<(int)(NB_INTERV);k++)
  //  { 
  //    MatPts[k][0]=(k/(float)(NB_INTERV))*cos((k*0.0001)*3.14159); 
  //    MatPts[k][1]=(k/(float)(NB_INTERV))*sin((k*0.001)*3.14159); 
      //>on peut essayer la ligne d'en bas aussi
      //MatPts[k][1]=(k/(float)(NB_INTERV))*sin((k*0.0001)*3.14159); 
  //   }


  //--Fin Question 1-----------------------------------------------------


  //>Affichage des Points dans MatPict
  plot_point(MatPts,MatPict,(int)(NB_INTERV));

  //>Save&Visu de MatPict
  SaveImagePgm((char*)OUTPUT_FILE,MatPict,HEIGHT,WIDTH);
  strcpy(BufSystVisu,VIEW_PGM);
  strcat(BufSystVisu," "); 	
  strcat(BufSystVisu,OUTPUT_FILE);
  strcat(BufSystVisu," &"); 
  system(BufSystVisu);

  //>Affiche Statistique
  printf("\n\n Stat:  Xmin=[%.2f] Xmax=[%.2f] Ymin=[%.2f] Ymax=[%.2f]\n",Xmin,Xmax,Ymin,Ymax);


 //--------------------------------------------------------------------------------
 //-------------- visu sous XWINDOW de l'évolution de MatPts ----------------------
 //--------------------------------------------------------------------------------
 if (flag_graph)
 {
 //>Uuverture Session Graphique
 if (open_display()<0) printf(" Impossible d'ouvrir une session graphique");
 sprintf(nomfen_ppicture,"Évolution du Graphe");
 win_ppicture=fabrique_window(nomfen_ppicture,10,10,HEIGHT,WIDTH,zoom);
 x_ppicture=cree_Ximage(MatPict,zoom,HEIGHT,WIDTH);

 printf("\n\n Pour quitter,appuyer sur la barre d'espace");
 fflush(stdout);

 //>Boucle Evolution
  for(i=0;i<HEIGHT;i++) for(j=0;j<WIDTH;j++) MatPict[i][j]=GREYWHITE;
  for(k=0;;)
     {   
       k=((k+EVOL_GRAPH)%(int)(NB_INTERV));
       Fill_Pict(MatPts,MatPict,k,(int)(NB_INTERV));
       XDestroyImage(x_ppicture);
       x_ppicture=cree_Ximage(MatPict,zoom,HEIGHT,WIDTH);
       XPutImage(display,win_ppicture,gc,x_ppicture,0,0,0,0,x_ppicture->width,x_ppicture->height); 
       usleep(10000);  //si votre machine est lente mettre un nombre plus petit
     }
 } 
       
 //>Retour  
 printf("\n Fini... \n\n\n");
 return 0;
}




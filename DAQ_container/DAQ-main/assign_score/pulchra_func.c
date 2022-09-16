#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "struct.h"
#include "mrc.h"
#include "pulchra_func.h"

#include "pulchra_common.h"
#include "nco_data.h"
#include "rot_data_coords.h"
#include "rot_data_idx.h"


//
// PULCHRA
// Protein Chain Restoration Algorithm
//
// Version 3.06
// December 2007
// Contact: Piotr Rotkiewicz, piotr -at- pirx -dot- com
//
// to compile:
// cc -O3 pulchra.c pulchra_data.c -lm -o pulchra
//
// to use:
// ./pulchra file.pdb
//
// to display available options:
// ./pulchra
//

#define COMPILE_BB
#define COMPILE_ROT

#include <sys/timeb.h>

#define uchar unsigned char
#define uint unsigned int
#define real double

#define PULCHRA_VERSION 3.06
#define MAX_BUF_SIZE 1000

#define FILE_SUCCESS     0
#define FILE_NOT_FOUND  -1
#define FILE_WARNING    -2

#define FATAL_MAE -1

#define FLAG_BACKBONE  1
#define FLAG_CALPHA    2
#define FLAG_SIDECHAIN 4
#define FLAG_SCM       8
#define FLAG_INITIAL  16

#define FLAG_PROTEIN  1
#define FLAG_DNA      2
#define FLAG_RNA      4
#define FLAG_CHYDRO   8

#define RADDEG 180./M_PI
#define DEGRAD M_PI/180.

int _VERBOSE = 0;
int _BB_REARRANGE = 1;
int _BB_OPTIMIZE = 0;
int _CA_OPTIMIZE = 1;
int _CA_RANDOM = 0;
int _CA_ITER = 100;
int _CA_TRAJECTORY = 0;
int _CISPRO = 0;
int _CHIRAL = 1;
int _CENTER_CHAIN = 0;
int _REBUILD_BB = 1;
int _REBUILD_SC = 1;
int _REBUILD_H = 0;
int _PDB_SG = 0;
int _TIME_SEED = 0;
int _XVOLUME = 1;
int _XVOL_ITER = 3;
int _PRESERVE = 1;
//real _CA_START_DIST = 3.0; //?? bug??
real _CA_START_DIST = 0.5; //?? bug??
real _CA_XVOL_DIST = 3.5;
real _SG_XVOL_DIST = 1.6;

#define CALC_C_ALPHA
#define CALC_C_ALPHA_ANGLES
#define CALC_C_ALPHA_START
#define CALC_C_ALPHA_XVOL

real CA_K=10.0;
real CA_ANGLE_K=20.0;
real CA_START_K=0.01;
real CA_XVOL_K=10.00;

#define CA_DIST 3.8
#define CA_DIST_TOL 0.1
#define CA_DIST_CISPRO 2.9
#define CA_DIST_CISPRO_TOL 0.1
#define E_EPS 1e-10

#ifndef bool
#define bool int
#endif

#ifndef true
#define true 1
#endif

#ifndef false
#define false 0
#endif

//Global should not use!!
//int **RBINS = NULL;
//real **X_COORDS = NULL;
//real **C_ALPHA = NULL;

// grid resolution (used for fast clash detection)
#define GRID_RES 6.0

#define HBOND_D2 11.5 //3.4*3.4

int chain_length = 0;

extern CMD cmd;

/*
char RES_NAMES[21][4] =
  { "GLY", "ALA", "SER", "CYS", "VAL",
    "THR", "ILE", "PRO", "MET", "ASP",
    "ASN", "LEU", "LYS", "GLU", "GLN",
    "ARG", "HIS", "PHE", "TYR", "TRP",
    "UNK" };
*/
char SHORT_AA_NAMES[22] = { "GASCVTIPMDNLKEQRHFYWX" };

int AA_NUMS[256];

int nheavy[20] = { 0, 1, 2, 2, 3, 3, 4, 3, 4, 4, 4, 4, 5, 5, 5, 7, 6, 7, 8, 10};

char *backbone_atoms[4] = { "N  ", "CA ", "C  ", "O  " };

char *heavy_atoms[200]= {
/* GLY */  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ALA */ "CB ", NULL,   NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* SER */ "CB ", "OG ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* CYS */ "CB ", "SG ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* VAL */ "CB ", "CG1", "CG2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* THR */ "CB ", "OG1", "CG2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ILE */ "CB ", "CG1", "CG2", "CD1",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* PRO */ "CB ", "CG ", "CD ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* MET */ "CB ", "CG ", "SD ", "CE ",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ASP */ "CB ", "CG ", "OD1", "OD2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* ASN */ "CB ", "CG ", "OD1", "ND2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* LEU */ "CB ", "CG ", "CD1", "CD2",  NULL,  NULL,  NULL,  NULL,  NULL,  NULL,
/* LYS */ "CB ", "CG ", "CD ", "CE ", "NZ ",  NULL,  NULL,  NULL,  NULL,  NULL,
/* GLU */ "CB ", "CG ", "CD ", "OE1", "OE2",  NULL,  NULL,  NULL,  NULL,  NULL,
/* GLN */ "CB ", "CG ", "CD ", "OE1", "NE2",  NULL,  NULL,  NULL,  NULL,  NULL,
/* ARG */ "CB ", "CG ", "CD ", "NE ", "CZ ", "NH1", "NH2",  NULL,  NULL,  NULL,
/* HIS */ "CB ", "CG ", "ND1", "CD2", "CE1", "NE2",  NULL,  NULL,  NULL,  NULL,
/* PHE */ "CB ", "CG ", "CD1", "CD2", "CE1", "CE2", "CZ ",  NULL,  NULL,  NULL,
/* TYR */ "CB ", "CG ", "CD1", "CD2", "CE1", "CE2", "CZ ", "OH ",  NULL,  NULL,
/* TRP */ "CB ", "CG ", "CD1", "CD2", "NE1", "CE2", "CE3", "CZ2", "CZ3", "CH2"};

/* reads full-atom pdb file */

struct _res_type;

typedef struct _atom_type {
  struct _atom_type *next;
  real x, y, z;
  char *name;
  int num, locnum;
  int flag;
  char cispro;
  int gx, gy, gz;
  struct _res_type *res;
  struct _atom_type *prev;
} atom_type;

typedef struct _res_type {
  struct _res_type *next;
  atom_type *atoms;
  int num, locnum, natoms;
  int type;
  char pdbsg;
  char protein;
  char *name;
  char chain;
  real sgx, sgy, sgz;
  real cmx, cmy, cmz;
  struct _res_type *prev;
} res_type;

typedef struct _mol_type {
  struct _mol_type *next;
  res_type *residua;
  int nres;
  unsigned char *r14;
  char *name;
  uchar *seq;
  char **contacts;
  real **cutoffs;
  struct _mol_type *prev;
} mol_type;

#define MIN(a,b) (a<b?a:b)
#define MAX(a,b) (a>b?a:b)

mol_type *chain = NULL;

static double rnd_op=1.0000/(double)1024;

double rnd(void)
{
 //return 0.01;
 return 0.001*(real)(rand()%1000);//too slow ??
 //return 1.0000/(double)1024*(double)(rand()&1023);//too slow ??
 //printf("%f\n",rnd_op*(double)rand());
}


// superimposition of two sets for coordinates + optional transformation of tpoints

real superimpose2(real coords1[][3], real coords2[][3], int npoints, real tpoints[][3], int ntpoints)
{
  real mat_s[3][3], mat_a[3][3], mat_b[3][3], mat_g[3][3];
  real mat_u[3][3];
  real tmp_mat[3][3];
  real val, d, alpha, beta, gamma, x, y, z;
  real cx1, cy1, cz1, cx2, cy2, cz2, tmpx, tmpy, tmpz;
  int i, j, k, n;

  real rvpoints=1.00/(real)npoints;

	//puts("OK");
    cx1=cy1=cz1=cx2=cy2=cz2=0.;

    //printf("%d %f\n",npoints,rvpoints);
    for (i=0; i<npoints; i++) {
      cx1+=coords1[i][0];
      cy1+=coords1[i][1];
      cz1+=coords1[i][2];
      cx2+=coords2[i][0];
      cy2+=coords2[i][1];
      cz2+=coords2[i][2];
    }


    //cx1/=(real)npoints;
    cx1*=rvpoints;
    //cy1/=(real)npoints;
    cy1*=rvpoints;;
    //cz1/=(real)npoints;
    cz1*=rvpoints;;

    //cx2/=(real)npoints;
    cx2*=rvpoints;
    //cy2/=(real)npoints;
    cy2*=rvpoints;
    //cz2/=(real)npoints;
    cz2*=rvpoints;

    for (i=0; i<npoints; i++) {
      coords1[i][0]-=cx1;
      coords1[i][1]-=cy1;
      coords1[i][2]-=cz1;
      coords2[i][0]-=cx2;
      coords2[i][1]-=cy2;
      coords2[i][2]-=cz2;
    }

    for (i=0; i<ntpoints; i++) {
      tpoints[i][0]-=cx2;
      tpoints[i][1]-=cy2;
      tpoints[i][2]-=cz2;
    }


    for (i=0; i<3; i++)
      for (j=0; j<3; j++) {
        //if (i==j)
        //  mat_s[i][j]=mat_a[i][j]=mat_b[i][j]=mat_g[i][j]=1.0;
        //else
        //  mat_s[i][j]=mat_a[i][j]=mat_b[i][j]=mat_g[i][j]=0.0;
        mat_u[i][j]=0.;
      }

    mat_s[0][1]=mat_a[0][1]=mat_b[0][1]=mat_g[0][1]=0.0;
    mat_s[0][2]=mat_a[0][2]=mat_b[0][2]=mat_g[0][2]=0.0;
    mat_s[1][0]=mat_a[1][0]=mat_b[1][0]=mat_g[1][0]=0.0;
    mat_s[1][2]=mat_a[1][2]=mat_b[1][2]=mat_g[1][2]=0.0;
    mat_s[2][0]=mat_a[2][0]=mat_b[2][0]=mat_g[2][0]=0.0;
    mat_s[2][1]=mat_a[2][1]=mat_b[2][1]=mat_g[2][1]=0.0;

    //init
    mat_s[0][0]=mat_a[0][0]=mat_b[0][0]=mat_g[0][0]=1.0;
    mat_s[1][1]=mat_a[1][1]=mat_b[1][1]=mat_g[1][1]=1.0;
    mat_s[2][2]=mat_a[2][2]=mat_b[2][2]=mat_g[2][2]=1.0;

    for (n=0; n<npoints; n++) {
      mat_u[0][0]+=coords1[n][0]*coords2[n][0];
      mat_u[0][1]+=coords1[n][0]*coords2[n][1];
      mat_u[0][2]+=coords1[n][0]*coords2[n][2];
      mat_u[1][0]+=coords1[n][1]*coords2[n][0];
      mat_u[1][1]+=coords1[n][1]*coords2[n][1];
      mat_u[1][2]+=coords1[n][1]*coords2[n][2];
      mat_u[2][0]+=coords1[n][2]*coords2[n][0];
      mat_u[2][1]+=coords1[n][2]*coords2[n][1];
      mat_u[2][2]+=coords1[n][2]*coords2[n][2];
    }


    for (i=0; i<3; i++)
      for (j=0; j<3; j++)
        tmp_mat[i][j]=0.;

    do {
      d=mat_u[2][1]-mat_u[1][2];
      if (d==0) alpha=0; else alpha=atan(d/(mat_u[1][1]+mat_u[2][2]));
      if (cos(alpha)*(mat_u[1][1]+mat_u[2][2])+sin(alpha)*(mat_u[2][1]-mat_u[1][2])<0.0)       alpha+=M_PI;
      mat_a[1][1]=mat_a[2][2]=cos(alpha);
      mat_a[2][1]=sin(alpha);
      mat_a[1][2]=-mat_a[2][1];

      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_a[j][k];

      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_a[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      d=mat_u[0][2]-mat_u[2][0];
      if (d==0) beta=0; else beta=atan(d/(mat_u[0][0]+mat_u[2][2]));
      if (cos(beta)*(mat_u[0][0]+mat_u[2][2])+sin(beta)*(mat_u[0][2]-mat_u[2][0])<0.0) beta+=M_PI;
      mat_b[0][0]=mat_b[2][2]=cos(beta);
      mat_b[0][2]=sin(beta);
      mat_b[2][0]=-mat_b[0][2];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_b[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_b[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      d=mat_u[1][0]-mat_u[0][1];
      if (d==0) gamma=0; else gamma=atan(d/(mat_u[0][0]+mat_u[1][1]));
      if (cos(gamma)*(mat_u[0][0]+mat_u[1][1])+sin(gamma)*(mat_u[1][0]-mat_u[0][1])<0.0)
        gamma+=M_PI;
      mat_g[0][0]=mat_g[1][1]=cos(gamma);
      mat_g[1][0]=sin(gamma);
      mat_g[0][1]=-mat_g[1][0];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_u[i][k]*mat_g[j][k];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_u[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      for (i=0; i<3; i++)
        for (j=0; j<3; j++)
          for (k=0; k<3; k++)
            tmp_mat[i][j]+=mat_g[i][k]*mat_s[k][j];
      for (i=0; i<3; i++)
        for (j=0; j<3; j++) {
          mat_s[i][j]=tmp_mat[i][j];
          tmp_mat[i][j]=0.;
        }
      val=fabs(alpha)+fabs(beta)+fabs(gamma);
    } while (val>0.001);

    val=0.;
    for (i=0; i<npoints; i++) {
      x=coords2[i][0];
      y=coords2[i][1];
      z=coords2[i][2];
      tmpx=x*mat_s[0][0]+y*mat_s[0][1]+z*mat_s[0][2];
      tmpy=x*mat_s[1][0]+y*mat_s[1][1]+z*mat_s[1][2];
      tmpz=x*mat_s[2][0]+y*mat_s[2][1]+z*mat_s[2][2];
      x=coords1[i][0]-tmpx;
      y=coords1[i][1]-tmpy;
      z=coords1[i][2]-tmpz;
      val+=x*x+y*y+z*z;
    }

    for (i=0; i<ntpoints; i++) {
      x=tpoints[i][0];
      y=tpoints[i][1];
      z=tpoints[i][2];
      tpoints[i][0]=x*mat_s[0][0]+y*mat_s[0][1]+z*mat_s[0][2];
      tpoints[i][1]=x*mat_s[1][0]+y*mat_s[1][1]+z*mat_s[1][2];
      tpoints[i][2]=x*mat_s[2][0]+y*mat_s[2][1]+z*mat_s[2][2];
    }

    for (i=0; i<npoints; i++) {
      coords1[i][0]+=cx1;
      coords1[i][1]+=cy1;
      coords1[i][2]+=cz1;
      coords2[i][0]+=cx2;
      coords2[i][1]+=cy2;
      coords2[i][2]+=cz2;
    }

    for (i=0; i<ntpoints; i++) {
      tpoints[i][0]+=cx1;
      tpoints[i][1]+=cy1;
      tpoints[i][2]+=cz1;
    }

  //return sqrt(val/(real)npoints);
  return sqrt(val*rvpoints);
}

//coords1 and coords2 are unit vectors
real JustImpose(double cd1[3][3], double tpoints[][3], int ntpoints)
{
/*
cd={unitv_x,unitv_y,unitv_z}
*/
 int i;
 double tmp[3];
 for(i=0;i<ntpoints;i++){
  tmp[0]=cd1[0][0]*tpoints[i][0]+cd1[1][0]*tpoints[i][1]+cd1[2][0]*tpoints[i][2];
  tmp[1]=cd1[0][1]*tpoints[i][0]+cd1[1][1]*tpoints[i][1]+cd1[2][1]*tpoints[i][2];
  tmp[2]=cd1[0][2]*tpoints[i][0]+cd1[1][2]*tpoints[i][1]+cd1[2][2]*tpoints[i][2];

  tpoints[i][0]=tmp[0];
  tpoints[i][1]=tmp[1];
  tpoints[i][2]=tmp[2];
 }
 return 0;
}


   

//int BuildCa(double ca[5][MAX_FLEN][3],double cd[MAX_FLEN][3], int Ncd,int len,MRC *map,MEMO *mm){
//int BuildCa(double **ca,double **cd, int Ncd,int len,MRC *map,MEMO *mm){
int BuildCa(MEMO *mm, int Ncd,int len,MRC *map){
 double ene1=0;
 double ene2=99999;
 double dcut,d2,dcut2;
 int cnt,pre;
 int Mcnt=0;
 double **cd=mm->cd;
 double **ca=mm->frag;
 //build 10 best C-alpha fragments

 //puts("#Build C-alpha");
 //Ca-Ca=3.2~3.8~4.2
 for(int p=0;p<10;p++){
  dcut=3.20+0.100*(double)p;
  dcut2=dcut*dcut;

  //init_cd[0][0]=cd[0][0];
  mm->tmp[0][0]=cd[0][0];
  //init_cd[0][1]=cd[0][1];
  mm->tmp[0][1]=cd[0][1];
  //init_cd[0][2]=cd[0][2];
  mm->tmp[0][2]=cd[0][2];

  cnt=1;
  pre=0;
 	for(int i=1;i<Ncd;i++){
	 d2=(cd[pre][0]-cd[i][0])*(cd[pre][0]-cd[i][0])+
	    (cd[pre][1]-cd[i][1])*(cd[pre][1]-cd[i][1])+
	    (cd[pre][2]-cd[i][2])*(cd[pre][2]-cd[i][2]);
	 if(d2>dcut2){

	  mm->tmp[cnt][0]=cd[i][0];
  	  mm->tmp[cnt][1]=cd[i][1];
  	  mm->tmp[cnt][2]=cd[i][2];



	  pre=i;
	  cnt++;
	 }
	 if(cnt==len)
	  break;
 	}
  if(cnt==len){



   //ene=Ca_Optimizer(init_cd,new_cd,len);
   
   //printf("ENE%d= %f\n",p,ene1);
   ene1=Ca_Optimizer(mm->tmp,mm->cur,len,mm);


   if(ene1<ene2){
    ene2=ene1;
    //update
    for(int j=0;j<len;j++){
     ca[j][0]=mm->cur[j][0];
     ca[j][1]=mm->cur[j][1];
     ca[j][2]=mm->cur[j][2];
    }
    Mcnt++;
   }
  }
 }


/*
 for(int i=0;i<Ncd;i++){
   printf("ATOM  %5d  CA  ALA%6d    ",i+1,i+1);
   printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",cd[i][0],cd[i][1],cd[i][2],1.0,1.0);
  }
 puts("ENDMDL");

 for(int i=0;i<len;i++){
   printf("ATOM  %5d  CA  ALA%6d    ",i+1,i+1);
   printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",ca[i][0],ca[i][1],ca[i][2],1.0,1.0);
  }
 puts("ENDMDL");
*/

 return Mcnt;
}


//buid side-chain and optimization
double OptMainSide(double ca[MAX_FLEN][3],int *seq,int len, MRC *map){
 double ene=0;
 int i,pos,iter=0;
 double casco,ssco,msco,Ene;
 double gscale=0.001;
 double tmp_ca[MAX_FLEN][3];
 double out[MAX_FLEN][3];
 double vec[3],pre_casco;
 double grid=map->widthx;

 for(i=0;i<len+3;i++){
  out[i][0]=ca[i][0];
  out[i][1]=ca[i][1];
  out[i][2]=ca[i][2];
 }

 //initial score
 //ssco=rebuild_sidechain_fromCA(ca,seq,len,map,false);
 casco=CaEne(ca,len+3);
 msco=rebuild_backbone_fromCA(ca,seq,len,map,false);

 pre_casco=casco; 

 Ene=casco-gscale*(ssco+msco);
 double NowEne=Ene;
 printf("ssco= %f msco= %f casco= %f Ene= %f\n",ssco,msco,casco,Ene);
 while(iter<100){

/*
	for(pos=0;pos<len+3;pos++){
  	 vec[0]=0.1-0.2*rnd();
  	 vec[1]=0.1-0.2*rnd();
  	 vec[2]=0.1-0.2*rnd();

	 tmp_ca[pos][0]=ca[pos][0]+vec[0];
	 tmp_ca[pos][1]=ca[pos][1]+vec[1];
	 tmp_ca[pos][2]=ca[pos][2]+vec[2];

	}
	casco=CaEne(tmp_ca,len+3);
	if(casco >=pre_casco + 0.1)
	 continue;

	ssco=rebuild_sidechain_fromCA(tmp_ca,seq,len,map);
 	msco=rebuild_backbone_fromCA(tmp_ca,seq,len,map);
	Ene=casco-gscale*(ssco+msco);

	if(Ene < NowEne){
	 NowEne=Ene;
	 //update
	  for(i=0;i<len+3;i++){
	   ca[i][0]=tmp_ca[i][0];
	   ca[i][1]=tmp_ca[i][1];
	   ca[i][2]=tmp_ca[i][2];
	  }
	  pre_casco=casco;
  	  printf("Up:ssco= %f msco= %f casco= %f Ene= %f\n",ssco,msco,casco,Ene);
	 }else{
  	  //printf("NG:ssco= %f msco= %f casco= %f Ene= %f\n",ssco,msco,casco,Ene);
	  continue;
	 }
*/



  //cal gradients
	for(pos=0;pos<len+3;pos++){
  	 vec[0]=0.5-1.0*rnd();
  	 vec[1]=0.5-1.0*rnd();
  	 vec[2]=0.5-1.0*rnd();

		for(i=0;i<len+3;i++){
		 tmp_ca[i][0]=out[i][0];
		 tmp_ca[i][1]=out[i][1];
		 tmp_ca[i][2]=out[i][2];
		}


  //printf("pos= %d vec= %f %f %f\n",pos,vec[0],vec[1],vec[2]);
	 tmp_ca[pos][0]+=vec[0];
	 tmp_ca[pos][1]+=vec[1];
	 tmp_ca[pos][2]+=vec[2];

 	 casco=CaEne(tmp_ca,len+3);
	 if(casco >=pre_casco + 0.01){
	  tmp_ca[pos][0]-=2.0*vec[0];
	  tmp_ca[pos][1]-=2.0*vec[1];
	  tmp_ca[pos][2]-=2.0*vec[2];
 	  casco=CaEne(tmp_ca,len+3);
	  if(casco >=pre_casco + 0.01)
	   continue;
	 }

	 //ssco=rebuild_sidechain_fromCA(tmp_ca,seq,len,map,false);
 	 msco=rebuild_backbone_fromCA(tmp_ca,seq,len,map,false);
	 Ene=casco-gscale*(ssco+msco);

	 if(Ene < NowEne){
	  NowEne=Ene;
	  //update
	  for(i=0;i<len+3;i++){
	   out[i][0]=tmp_ca[i][0];
	   out[i][1]=tmp_ca[i][1];
	   out[i][2]=tmp_ca[i][2];
	  }
	  pre_casco=casco;
  	  printf("Up:ssco= %f msco= %f casco= %f Ene= %f\n",ssco,msco,casco,Ene);
	  //ssco=rebuild_sidechain_fromCA(tmp_ca,seq,len,map,false);
	 }else{
	  continue;
	 }
	}

  iter++;
 }
 //printf("NewEne= %f\n",NowEne);
 return ene;
}

// energy calculation for C-alpha optimizer
//double calc_ca_energy(int chain_length,double now[MAX_FLEN][3],double new[MAX_FLEN][3], double init[MAX_FLEN][3], double gradient[MAX_FLEN][3], double alpha, double *ene, bool calc_gradient)
double calc_ca_energy(int chain_length,double **now,double **new, double **init, double **gradient, double alpha, double *ene, bool calc_gradient)
{
  int i, j;
  double dx, dy, dz;
  double dist, ddist, ddist2,dist2;
  double new_e_pot;
  double theta0, tdif, th, aa, bb, ab;
  double ff0, ff2, dth, m0, m2, grad, f0[3], f2[3];
  double adiff[3], bdiff[3];
  double deriv, theta, dtheta, len1, len2, cos_theta, sin_theta;
  double dx1, dy1, dz1;
  double dx2, dy2, dz2;
  double dx3, dy3, dz3;
  double vx1, vy1, vz1;
  double vx2, vy2, vz2;
  double vx3, vy3, vz3;

  double r12x, r12y, r12z;
  double r32x, r32y, r32z;
  double d12, d32, d12inv, d32inv, c1, c2, diff;
  double f1x, f1y, f1z;
  double f2x, f2y, f2z;
  double f3x, f3y, f3z;

        for (i=0; i<chain_length; i++) {
          new[i][0]=now[i][0]+alpha*gradient[i][0];
          new[i][1]=now[i][1]+alpha*gradient[i][1];
          new[i][2]=now[i][2]+alpha*gradient[i][2];
        }

        new_e_pot = 0.0;

        ene[0]=ene[1]=ene[2]=ene[3]=0.0;

        for (i=0; i<chain_length; i++) {
          dx=new[i][0]-init[i][0];
          dy=new[i][1]-init[i][1];
          dz=new[i][2]-init[i][2];

          dist=sqrt(dx*dx+dy*dy+dz*dz);
          //dist2=dx*dx+dy*dy+dz*dz;
          ddist = -dist;
          if (dist>_CA_START_DIST) {
            ddist2=dist*dist;
            new_e_pot+=CA_START_K*ddist2;
            ene[1] += CA_START_K*ddist2;
            if (calc_gradient) {
              grad = ddist * (-2.0*CA_START_K)/dist;
              gradient[i][0]-=grad*dx;
              gradient[i][1]-=grad*dy;
              gradient[i][2]-=grad*dz;
            }
          }

          if (i>0) {
            dx=new[i][0]-new[i-1][0];
            dy=new[i][1]-new[i-1][1];
            dz=new[i][2]-new[i-1][2];
            dist=sqrt(dx*dx+dy*dy+dz*dz);
	    //printf("dist= %f xyz %f %f %f\n",dist,new[i][0],new[i][1],new[i][2]);
            //dist2=dx*dx+dy*dy+dz*dz;
            //if (c_alpha[i]->cispro) {
            //  ddist=CA_DIST_CISPRO-dist;
            //} else {
              ddist=CA_DIST-dist;
            //}
            ddist2=ddist*ddist;
            new_e_pot+=CA_K*ddist2;
            ene[0] += CA_K*ddist2;
            if (calc_gradient) {
              grad = ddist * (-2.0*CA_K)/dist;
              gradient[i][0]-=grad*dx;
              gradient[i][1]-=grad*dy;
              gradient[i][2]-=grad*dz;
              gradient[i-1][0]+=grad*dx;
              gradient[i-1][1]+=grad*dy;
              gradient[i-1][2]+=grad*dz;
            }
          }

          for (j=0;j<i;j++) {
            if (abs(i-j)>2) {
              dx=new[i][0]-new[j][0];
              dy=new[i][1]-new[j][1];
              dz=new[i][2]-new[j][2];
              dist=sqrt(dx*dx+dy*dy+dz*dz);
              ddist = dist-_CA_XVOL_DIST;
              if (dist<_CA_XVOL_DIST) {
                ddist2 = dist*dist;
                new_e_pot+=CA_XVOL_K*ddist2;
                ene[3] += CA_XVOL_K*ddist2;
                if (calc_gradient) {
                  grad = ddist*(8.0*CA_XVOL_K)/dist;
                  gradient[i][0]-=grad*dx;
                  gradient[i][1]-=grad*dy;
                  gradient[i][2]-=grad*dz;
                  gradient[j][0]+=grad*dx;
                  gradient[j][1]+=grad*dy;
                  gradient[j][2]+=grad*dz;
                }
              }
            }
          }

        if (i>0 && i<chain_length-1) {
          r12x=new[i-1][0]-new[i][0];
          r12y=new[i-1][1]-new[i][1];
          r12z=new[i-1][2]-new[i][2];
          r32x=new[i+1][0]-new[i][0];
          r32y=new[i+1][1]-new[i][1];
          r32z=new[i+1][2]-new[i][2];
          d12 = sqrt(r12x*r12x+r12y*r12y+r12z*r12z);
          d32 = sqrt(r32x*r32x+r32y*r32y+r32z*r32z);
          cos_theta = (r12x*r32x+r12y*r32y+r12z*r32z)/(d12*d32);

          if (cos_theta>1.0)
            cos_theta = 1.0;
          else
          if (cos_theta<-1.0)
            cos_theta = -1.0;

          sin_theta = sqrt(1.0-cos_theta*cos_theta);
          theta = acos(cos_theta);

	  if(sin_theta==0.00){
	   //printf("ZERO theta %f %f cos %f\n",d12,d32,cos_theta);
	   sin_theta+=0.001;//New!!
	  }

          if (RADDEG*theta<80.)
            diff = theta-80.*DEGRAD;
          else
          if (RADDEG*theta>150.)
            diff = theta-150.*DEGRAD;
          else
            diff = 0.0;

          new_e_pot += CA_ANGLE_K*diff*diff;
          ene[2] += CA_ANGLE_K*diff*diff;
          if (calc_gradient) {
            d12inv = 1.0/d12;
            d32inv = 1.0/d32;
            diff *= (-2.0*CA_ANGLE_K)/sin_theta;//!!!!!
            c1 = diff*d12inv;
            c2 = diff*d32inv;
            f1x = c1*(r12x*(d12inv*cos_theta)-r32x*d32inv);
            f1y = c1*(r12y*(d12inv*cos_theta)-r32y*d32inv);
            f1z = c1*(r12z*(d12inv*cos_theta)-r32z*d32inv);
            f3x = c2*(r32x*(d32inv*cos_theta)-r12x*d12inv);
            f3y = c2*(r32y*(d32inv*cos_theta)-r12y*d12inv);
            f3z = c2*(r32z*(d32inv*cos_theta)-r12z*d12inv);
            f2x = -f1x-f3x;
            f2y = -f1y-f3y;
            f2z = -f1z-f3z;
            gradient[i-1][0]+=f1x;
            gradient[i-1][1]+=f1y;
            gradient[i-1][2]+=f1z;
            gradient[i][0]+=f2x;
            gradient[i][1]+=f2y;
            gradient[i][2]+=f2z;
            gradient[i+1][0]+=f3x;
            gradient[i+1][1]+=f3y;
            gradient[i+1][2]+=f3z;
          }
        }


        }

//printf("ene[3] = %f\n", ene[3]);

  return new_e_pot;
}

double CaEne(double new[MAX_FLEN][3],int N)
{
  int i, j;
  double dx, dy, dz;
  double dist, ddist, ddist2,dist2;
  double new_e_pot;
  double theta0, tdif, th, aa, bb, ab;
  double ff0, ff2, dth, m0, m2, grad, f0[3], f2[3];
  double adiff[3], bdiff[3];
  double deriv, theta, dtheta, len1, len2, cos_theta, sin_theta;
  double dx1, dy1, dz1;
  double dx2, dy2, dz2;
  double dx3, dy3, dz3;
  double vx1, vy1, vz1;
  double vx2, vy2, vz2;
  double vx3, vy3, vz3;

  double r12x, r12y, r12z;
  double r32x, r32y, r32z;
  double d12, d32, d12inv, d32inv, c1, c2, diff;
  double f1x, f1y, f1z;
  double f2x, f2y, f2z;
  double f3x, f3y, f3z;

  double ene[4];

/*
        for (i=0; i<chain_length; i++) {
          new[i][0]=now[i][0]+alpha*gradient[i][0];
          new[i][1]=now[i][1]+alpha*gradient[i][1];
          new[i][2]=now[i][2]+alpha*gradient[i][2];
        }
*/
        new_e_pot = 0.0;

        ene[0]=ene[1]=ene[2]=ene[3]=0.0;


        for (i=0; i<N; i++) {
/*
          dx=new[i][0]-init[i][0];
          dy=new[i][1]-init[i][1];
          dz=new[i][2]-init[i][2];

          dist=sqrt(dx*dx+dy*dy+dz*dz);
          //dist2=dx*dx+dy*dy+dz*dz;
          ddist = -dist;
          if (dist>_CA_START_DIST) {
            ddist2=dist*dist;
            new_e_pot+=CA_START_K*ddist2;
            ene[1] += CA_START_K*ddist2;
            if (calc_gradient) {
              grad = ddist * (-2.0*CA_START_K)/dist;
              gradient[i][0]-=grad*dx;
              gradient[i][1]-=grad*dy;
              gradient[i][2]-=grad*dz;
            }
          }
*/
	  //CA-CA distance
          if (i>0) {
            dx=new[i][0]-new[i-1][0];
            dy=new[i][1]-new[i-1][1];
            dz=new[i][2]-new[i-1][2];
            dist=sqrt(dx*dx+dy*dy+dz*dz);
            ddist=CA_DIST-dist;
            ddist2=ddist*ddist;
            new_e_pot+=CA_K*ddist2;
            ene[0] += CA_K*ddist2;
          }

	  //clash
          for (j=0;j<i;j++) {
            if (abs(i-j)>2) {
              dx=new[i][0]-new[j][0];
              dy=new[i][1]-new[j][1];
              dz=new[i][2]-new[j][2];
              dist=sqrt(dx*dx+dy*dy+dz*dz);
              ddist = dist-_CA_XVOL_DIST;
              if (dist<_CA_XVOL_DIST) {
                ddist2 = dist*dist;
                new_e_pot+=CA_XVOL_K*ddist2;
                ene[3] += CA_XVOL_K*ddist2;
              }
            }
          }

	//Angle
        if (i>0 && i<N-1) {
          r12x=new[i-1][0]-new[i][0];
          r12y=new[i-1][1]-new[i][1];
          r12z=new[i-1][2]-new[i][2];
          r32x=new[i+1][0]-new[i][0];
          r32y=new[i+1][1]-new[i][1];
          r32z=new[i+1][2]-new[i][2];
          d12 = sqrt(r12x*r12x+r12y*r12y+r12z*r12z);
          d32 = sqrt(r32x*r32x+r32y*r32y+r32z*r32z);
          cos_theta = (r12x*r32x+r12y*r32y+r12z*r32z)/(d12*d32);
          if (cos_theta>1.0)
            cos_theta = 1.0;
          else
          if (cos_theta<-1.0)
            cos_theta = -1.0;
          sin_theta = sqrt(1.0-cos_theta*cos_theta);
          theta = acos(cos_theta);

          if (RADDEG*theta<80.)
            diff = theta-80.*DEGRAD;
          else
          if (RADDEG*theta>150.)
            diff = theta-150.*DEGRAD;
          else
            diff = 0.0;

          new_e_pot += CA_ANGLE_K*diff*diff;
          ene[2] += CA_ANGLE_K*diff*diff;
/*
          if (calc_gradient) {
            d12inv = 1.0/d12;
            d32inv = 1.0/d32;
            diff *= (-2.0*CA_ANGLE_K)/sin_theta;
            c1 = diff*d12inv;
            c2 = diff*d32inv;
            f1x = c1*(r12x*(d12inv*cos_theta)-r32x*d32inv);
            f1y = c1*(r12y*(d12inv*cos_theta)-r32y*d32inv);
            f1z = c1*(r12z*(d12inv*cos_theta)-r32z*d32inv);
            f3x = c2*(r32x*(d32inv*cos_theta)-r12x*d12inv);
            f3y = c2*(r32y*(d32inv*cos_theta)-r12y*d12inv);
            f3z = c2*(r32z*(d32inv*cos_theta)-r12z*d12inv);
            f2x = -f1x-f3x;
            f2y = -f1y-f3y;
            f2z = -f1z-f3z;
            gradient[i-1][0]+=f1x;
            gradient[i-1][1]+=f1y;
            gradient[i-1][2]+=f1z;
            gradient[i][0]+=f2x;
            gradient[i][1]+=f2y;
            gradient[i][2]+=f2z;
            gradient[i+1][0]+=f3x;
            gradient[i+1][1]+=f3y;
            gradient[i+1][2]+=f3z;
          }
*/
        }


        }

//printf("ene[3] = %f\n", ene[3]);

  return new_e_pot;
}


//double Ca_Optimizer(double init_c_alpha[MAX_FLEN][3],double  new_c_alpha[MAX_FLEN][3],int Ncd){
double Ca_Optimizer(double **init_c_alpha,double  **new_c_alpha,int Ncd,MEMO *m){
 //char buf[1000];
  int i, j, hx, my_iter;
  double dx, dy, dz, dd, dist, dist2, dist3, ddist, ddist2;
  double e_pot, new_e_pot, grad, alpha, e_pot1, e_pot2, e_pot3;
  double adiff[3], bdiff[3];
  double ff0, ff2, aa, ab, bb, th, tdif, dth, m0, m2;
  double theta0, deg_th, maxgrad, sum;
  double f0[3], f2[3];
  double x, y, z;
  int numsteps, numsteps2, msteps;
  int *sec;
  //double **new_c_alpha, **gradient, **init_c_alpha, last_alpha, tmp, last_good_alpha, d_alpha, last_e_pot;
  double last_alpha, tmp, last_good_alpha, d_alpha, last_e_pot;
  //atom_type *atom, **c_alpha;
  //res_type *res;
  //FILE *inp, *out;
  int mnum, init, ok;
  double alpha1, alpha2, alpha3, a0;
  double ene1, ene2, ene3, e0;
  double energies[4];
  double w1, w2, w3, eps;
  double gnorm, last_gnorm;
  int mode, fcnt;
  //double now[MAX_FLEN][3];
  double **now=m->mch;

  int chain_length=Ncd;
  //double gradient[MAX_FLEN][3]={};
  double **gradient=m->vec;


 for (i=0; i<chain_length; i++){
   now[i][0]=init_c_alpha[i][0];
   now[i][1]=init_c_alpha[i][1];
   now[i][2]=init_c_alpha[i][2];
 }
  mnum = 1;
    mode = 0;
    init = 0;
    numsteps=numsteps2=0;
    last_alpha = 0.0;

 //printf("Optimizing alpha carbons...\n");
 
eps = 0.5;

    fcnt=0;

    last_gnorm = 1000.;

 do {
	// calculate gradients

      	e_pot=e_pot1=e_pot2=e_pot3=0.;

      	for (i=0; i<chain_length; i++)
         gradient[i][0]=gradient[i][1]=gradient[i][2]=0.;

	//Get gradient
      //e_pot = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, 0.0, energies, true);
      e_pot = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, 0.0, energies, true);

/*	
	if(100<gradient[3][0]){
	 for(int i=0;i<chain_length;i++)
	  printf("%dgra= %f %f %f\n",i,gradient[i][0],gradient[i][1],gradient[i][2]);
	 for(int i=0;i<chain_length;i++)
	  printf("%dinit %f %f %f\n",i,init_c_alpha[i][0],init_c_alpha[i][1],init_c_alpha[i][2]);
	}
*/
        //printf("Initial energy: bond=%.5lf angle=%.5f restraints=%.5f xvol=%.5f total=%.5f\n", energies[0], energies[2], energies[1], energies[3], e_pot);

      if (!init) init=1;

  	alpha1 = -1.0;
      	alpha2 = 0.0;
      	alpha3 = 1.0;

	ene1 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, alpha1, energies, false);
      	ene2 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, alpha2, energies, false);
      	ene3 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, alpha3, energies, false);


	
	msteps = 0;
	//Find alpha1 and 3 where ene2 < ene1 AND ene2 < ene3
      while (ene2>MIN(ene1,ene3) && msteps<_CA_ITER) {
        msteps++;
        alpha1 *= 2.0;
        alpha3 *= 2.0;
        //ene1 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha1, energies, false);
	ene1 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, alpha1, energies, false);
        //ene3 = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, alpha3, energies, false);
      	ene3 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, alpha3, energies, false);
      }



	

	msteps = 0;
      do {
        if (alpha3-alpha2>alpha2-alpha1) {
	  double e00;
          a0 = 0.5*(alpha2+alpha3);
          e0 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          //e0 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, a0-1e-5, energies, false);
          //e0 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, a0+1e-5, energies, false);
          //e0 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          if (e0<ene2) {
            alpha1 = alpha2;
            alpha2 = a0;
            ene1 = ene2;
            ene2 = e0;
          } else {
            alpha3 = a0;
            ene3 = e0;
          }
        } else {
          a0 = 0.5*(alpha1+alpha2);
          e0 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          //e0 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, a0-1e-5, energies, false);
          //e0 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, a0+1e-5, energies, false);
          //e0 = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, a0, energies, false);
          if (e0<ene2) {
            alpha3 = alpha2;
            alpha2 = a0;
            ene3 = ene2;
            ene2 = e0;
          } else {
            alpha1 = a0;
            ene1 = e0;
          }
        }
        msteps++;
      } while (alpha3-alpha1>1e-6 && msteps<20);


	last_alpha = alpha2;
      e_pot = ene2;


	//Update
      for (i=0; i<chain_length; i++) {
        //c_alpha[i]->x=c_alpha[i]->x+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][0];
        now[i][0]+=(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][0];
        //c_alpha[i]->y=c_alpha[i]->y+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][1];
        now[i][1]+=(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][1];
        //c_alpha[i]->z=c_alpha[i]->z+(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][2];
        now[i][2]+=(last_alpha+last_alpha*(rnd()-0.5)*eps)*gradient[i][2];
      }

      //e_pot = calc_ca_energy(c_alpha, new_c_alpha, init_c_alpha, gradient, 0.0, energies, false);
      e_pot = calc_ca_energy(chain_length,now, new_c_alpha, init_c_alpha, gradient, 0.0, energies, false);
      
 	//printf("Now energy: bond=%.5lf angle=%.5f restraints=%.5f xvol=%.5f total=%.5f\n", energies[0], energies[2], energies[1], energies[3], e_pot);

      eps *= 0.75;
      if (eps<1e-3) eps=0.0;

      numsteps++;

      gnorm = 0.0;

	  for (i=0; i<chain_length; i++) {
        gnorm += gradient[i][0]*gradient[i][0] + gradient[i][1]*gradient[i][1] + gradient[i][2]*gradient[i][2];
      }

      gnorm = sqrt(gnorm/(double)chain_length);

      if (last_gnorm-gnorm<1e-3) fcnt++;

      last_gnorm = gnorm;

 } while ( (fcnt<3) &&  (gnorm>0.01) && (numsteps<_CA_ITER));



 for (i=0; i<chain_length; i++){
   new_c_alpha[i][0]=now[i][0];
   new_c_alpha[i][1]=now[i][1];
   new_c_alpha[i][2]=now[i][2];
 }




 return e_pot;
}

// r14 chiral distance
real calc_r14(real x1, real y1, real z1,
	      real x2, real y2, real z2,
	      real x3, real y3, real z3,
	      real x4, real y4, real z4)
{
  real r, dx, dy, dz;
  real vx1, vy1, vz1, vx2, vy2, vz2, vx3, vy3, vz3;
  real hand;

    dx = x4-x1;
    dy = y4-y1;
    dz = z4-z1;

    r = sqrt(dx*dx+dy*dy+dz*dz);

    vx1=x2-x1;
    vy1=y2-y1;
    vz1=z2-z1;
    vx2=x3-x2;
    vy2=y3-y2;
    vz2=z3-z2;
    vx3=x4-x3;
    vy3=y4-y3;
    vz3=z4-z3;

    hand = (vy1*vz2-vy2*vz1)*vx3+
           (vz1*vx2-vz2*vx1)*vy3+
           (vx1*vy2-vx2*vy1)*vz3;

    if (hand<0) r=-r;

  return r;
}

// distance
real calc_distance(real x1, real y1, real z1,
							 		  real x2, real y2, real z2)
{
  real dx,dy,dz;
  real dist2;

    dx = (x1) - (x2);
    dy = (y1) - (y2);
    dz = (z1) - (z2);
    if (dx || dy || dz ) {
      dist2 = dx*dx+dy*dy+dz*dz;
      return (sqrt(dist2));
    } else
      return 0.0;
}


void prepare_rbins_ca(double **ca,int **rbins,int len){
 int i, j, k, l, m, bin13_1, bin13_2, bin14, found, pro;
  int b13_1, b13_2, b14;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real r13_1, r13_2, r14;
  int chain_length=len+4;
 

 for (i=2;i<chain_length-1;i++) {
    	x1 = ca[i-2][0];
    	y1 = ca[i-2][1];
    	z1 = ca[i-2][2];

    	x2 = ca[i-1][0];
    	y2 = ca[i-1][1];
    	z2 = ca[i-1][2];

    	x3 = ca[i][0];
    	y3 = ca[i][1];
    	z3 = ca[i][2];

    	x4 = ca[i+1][0];
    	y4 = ca[i+1][1];
    	z4 = ca[i+1][2];

    	r13_1 = calc_distance(x1, y1, z1, x3, y3, z3);
    	r13_2 = calc_distance(x2, y2, z2, x4, y4, z4);
    	r14 = calc_r14(x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4);

      bin13_1 = (int)((r13_1-4.6)/0.3);
      bin13_2 = (int)((r13_2-4.6)/0.3);
      bin14 = (int)((r14+11.)/0.3);

      if (bin13_1<0) bin13_1=0;
      if (bin13_2<0) bin13_2=0;
      if (bin14<0) bin14=0;
      if (bin13_1>9) bin13_1=9;
      if (bin13_2>9) bin13_2=9;
      if (bin14>73) bin14=73;

      rbins[i][0] = bin13_1;
      rbins[i][1] = bin13_2;
      rbins[i][2] = bin14;
 }
}

void cross(real *v1, real *v2, real *v3)
{
  v3[0] = v1[1]*v2[2]-v1[2]*v2[1];
  v3[1] = v1[2]*v2[0]-v1[0]*v2[2];
  v3[2] = v1[0]*v2[1]-v1[1]*v2[0];
}

void norm(real *v)
{
  real d;

    d = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    v[0] /= d;
    v[1] /= d;
    v[2] /= d;
}


//side-chain table
typedef struct{
 double sco;
 int idx,pos,rot;
} STBL;

int cmp_stbl(const void *c1, const void *c2){

 STBL a=*(STBL *)c1;
 STBL b=*(STBL *)c2;

 if(a.sco<b.sco) return 1;
 if(a.sco>b.sco) return -1;
 return 0;
}

double rebuild_sidechain_fromCA(double **ca,double ***out,int **rbins,int *seq,int len, MRC *map,MEMO *mm,bool SHOW){
 
  real **cacoords, **tmpcoords, **tmpstat;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real x5, y5, z5;
  real r14, r13_1, r13_2;
  real dx, dy, dz, dd;
  real hit, besthit;
  int exvol, bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14;
  real rmsd, total;
  real v1[3], v2a[3], v2b[3], v2[3], v3[3];
  int nsc, nca;
  real cax, cay, caz;
  //real **lsys, **vv, **sc;
  real lsys[3][3], vv[3][3],bsc[12][3];
  double sc[20][3];
  char scn[12][4];
  //rot_struct *rot;
  int ok, last_a, last_b, last_c, last_d, jpos;
  int jx, jy, jz, jxi, jyi, jzi, b13_1, b13_2, b14, jm;
  int crot, bestrot, minexvol, totexvol, rtried, pos, cpos;
  real cmx, cmy, cmz, ddx, ddy, ddz, ddd, bestdd;
  real sort_rot[100][2];
  //int rbins[MAX_FLEN][3];
  int pi;
  //top 10 data
  double SORTED_ROTAMERS[MAX_FLEN][10][3];//0:
  int rot;
  double sco, bestsco;
  double Sscore=0.00;
  STBL stbl[200];
  int Ns;
  double ****slib=mm->slib;
  double **stmp=mm->stmp;
  //index is 2,3,4...len-2
  prepare_rbins_ca(ca,rbins,len);

 //**seq      0,1,2....len-1
 //ca   -2,-1,0,.......len-1,len,len+1

 int anum=1;
 Ns=0;
 for(i=0;i<len;i++){
  if(seq[i]<1 || seq[i] > 19)//GLY or unknown
   continue;

  pi=i+2;//for ca

  	x1 = ca[pi-2][0];
  	y1 = ca[pi-2][1];
  	z1 = ca[pi-2][2];
  	x2 = ca[pi-1][0];
  	y2 = ca[pi-1][1];
  	z2 = ca[pi-1][2];
  	x3 = ca[pi][0];
  	y3 = ca[pi][1];
  	z3 = ca[pi][2];
  	x4 = ca[pi+1][0];
  	y4 = ca[pi+1][1];
  	z4 = ca[pi+1][2];

	bin13_1 = rbins[pi][0];
      	bin13_2 = rbins[pi][1];
      	bin14 = rbins[pi][2];

	v1[0] = x4-x2;
      	v1[1] = y4-y2;
      	v1[2] = z4-z2;

      	v2a[0] = x4-x3;
      	v2a[1] = y4-y3;
      	v2a[2] = z4-z3;

      	v2b[0] = x3-x2;
      	v2b[1] = y3-y2;
      	v2b[2] = z3-z2;

	cross(v2a, v2b, v2);
      	cross(v1, v2, v3);

      	norm(v1);
      	norm(v2);
      	norm(v3);

	//init
	for (j=0;j<10;j++)
         SORTED_ROTAMERS[i][j][0] = 500.;

	j = 0;
      	besthit = 1000.;
      	bestpos = 0;

	//search 10 closest rotamer

	for(int j=idx_idx[seq[i]];rot_stat_idx[j][0]==seq[i];j++){
         if (rot_stat_idx[j][0]==seq[i]) {
          hit =  fabs(rot_stat_idx[j][1]-bin13_1) 
		+fabs(rot_stat_idx[j][2]-bin13_2) 
		+0.2*fabs(rot_stat_idx[j][3]-bin14);
          if (hit<SORTED_ROTAMERS[i][9][0]) {
            k = 9;
            while (k>=0 && hit<SORTED_ROTAMERS[i][k][0]) {
              k--;
            }
            k++;
            // k = hit
            for (l=9;l>k;l--) {
              SORTED_ROTAMERS[i][l][0]=SORTED_ROTAMERS[i][l-1][0];
              SORTED_ROTAMERS[i][l][1]=SORTED_ROTAMERS[i][l-1][1];
            }
            SORTED_ROTAMERS[i][k][0]=hit;//hit-score
            SORTED_ROTAMERS[i][k][1]=j;  //index
          }
         }else if(rot_stat_idx[j][0]>seq[i]){
	  break;
	 }
	}

	//check 10 best rotemer

	for(rot=0;rot<10; rot++){

	 if(SORTED_ROTAMERS[i][rot][0]==500)
	  break;

	 besthit = SORTED_ROTAMERS[i][rot][0];
      	 bestpos = SORTED_ROTAMERS[i][rot][1];
	
	 pos = rot_stat_idx[bestpos][5];
      	 //nsc = nheavy[seq[i]]+1;
      	 nsc = nheavy[seq[i]];//ignore 0:C-alpha


	//get XYZ system
	 for (j=0;j<3;j++) {
       	  vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];
       	 }

	 //copy from DB
	 //ignore 0: C-alpha
      	 for (j=0;j<nsc;j++) {
          sc[j][0] = rot_stat_coords[pos+j+1][0];
          sc[j][1] = rot_stat_coords[pos+j+1][1];
          sc[j][2] = rot_stat_coords[pos+j+1][2];
      	 }

	 //faster ver
	 JustImpose(vv,sc,nsc);

	 	//shift side-chains
	 	for (j=0;j<nsc;j++) {
         	 sc[j][0] += x3;
         	 sc[j][1] += y3;
         	 sc[j][2] += z3;
       	 	}
		
	 //check MAP potential score MDFF like
	 sco=MapScore2(sc,nsc,map);
	 SORTED_ROTAMERS[i][rot][2]=sco;
	 stbl[Ns].sco=sco;
	 stbl[Ns].pos=i;
	 stbl[Ns].idx=bestpos;
	 stbl[Ns].rot=rot;
	 Ns++;

		for (j=0;j<nsc;j++) {
        	 slib[i][rot][j][0] =sc[j][0];
        	 slib[i][rot][j][1] =sc[j][1];
        	 slib[i][rot][j][2] =sc[j][2];
		}
	}
 
/*
	if(SHOW){
	//show atoms
	printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         anum++, "CA ", RES_NAMES[seq[i]], ' ',i+1,x3,y3,z3);
		for (j=0;j<nsc;j++) {
		 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
        	 anum++, heavy_atoms[seq[i]*10+j], RES_NAMES[seq[i]], ' ',i+1,bsc[j][0], bsc[j][1], bsc[j][2]);
       		}
	}
*/

  //Sscore+=bestsco;
 }
/*
	if(SHOW){
 	 puts("ENDMDL");
	}
*/


 //sort by sco
 //printf("Ns= %d\n",Ns);
 qsort(&stbl,Ns,sizeof(STBL),cmp_stbl);

 int Nclash=0;
 int Ntmp=0;//stmp
 Sscore=0;
 //Build side-chains
 for(int i=0;i<len;i++)
  out[i]=NULL;

 for(int i=0;i<Ns;i++){

  if(out[stbl[i].pos]!=NULL)
   continue;

  nsc = nheavy[seq[stbl[i].pos]];
  Nclash=ClashChk(slib[stbl[i].pos][stbl[i].rot],nsc,ca,len+4,stbl[i].pos+2,3.0);
  if(Nclash>1)
   continue;

  Nclash=ClashChk(slib[stbl[i].pos][stbl[i].rot],nsc,stmp,Ntmp,-1,3.0);
  if(Nclash>1)
   continue;

	//input
	for(int j=0;j<nsc;j++){
	 stmp[Ntmp][0]=slib[stbl[i].pos][stbl[i].rot][j][0];
	 stmp[Ntmp][1]=slib[stbl[i].pos][stbl[i].rot][j][1];
	 stmp[Ntmp][2]=slib[stbl[i].pos][stbl[i].rot][j][2];
	 Ntmp++;
	}
  out[stbl[i].pos]=slib[stbl[i].pos][stbl[i].rot];
  Sscore+=stbl[i].sco;
  //printf("%d %f pos= %d idx= %d rot= %d Nc= %d\n",i,stbl[i].sco,stbl[i].pos,stbl[i].idx,stbl[i].rot,Nclash);
 }
 //printf("SideScore= %f\n",Sscore);
 return Sscore;
}

double rebuild_sidechain_fromCA2(double **ca,double ***out,int **rbins,int *seq,int len, MRC *map,MEMO *mm){
 
  real **cacoords, **tmpcoords, **tmpstat;
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real x5, y5, z5;
  real r14, r13_1, r13_2;
  real dx, dy, dz, dd;
  real hit, besthit;
  int exvol, bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14;
  real rmsd, total;
  real v1[3], v2a[3], v2b[3], v2[3], v3[3];
  int nsc, nca;
  real cax, cay, caz;
  //real **lsys, **vv, **sc;
  real lsys[3][3], vv[3][3],bsc[12][3];
  double sc[20][3];
  char scn[12][4];
  //rot_struct *rot;
  int ok, last_a, last_b, last_c, last_d, jpos;
  int jx, jy, jz, jxi, jyi, jzi, b13_1, b13_2, b14, jm;
  int crot, bestrot, minexvol, totexvol, rtried, pos, cpos;
  real cmx, cmy, cmz, ddx, ddy, ddz, ddd, bestdd;
  real sort_rot[100][2];
  //int rbins[MAX_FLEN][3];
  int pi;
  //top 10 data
  double SORTED_ROTAMERS[MAX_FLEN][10][3];//0:
  int rot;
  double sco, bestsco;
  double Sscore=0.00;
  STBL stbl[200];
  int Ns;
  double ****slib=mm->slib;
  double **stmp=mm->stmp;
  //index is 2,3,4...len-2
  prepare_rbins_ca(ca,rbins,len-4);

 //**seq      0,1,2....len-1
 //ca   -2,-1,0,.......len-1,len,len+1

 //puts("OK");
 //return 0;
 int anum=1;
 Ns=0;
 //for(i=0;i<len;i++){
 //	printf("i=%d (%d)\n",i,seq[i]);
 //}
 for(i=2;i<len-2;i++){
	//printf("i=%d (%d) %d %d %d\n",i,seq[i],bin13_1,bin13_2,bin14);
  if(seq[i]<1 || seq[i] > 19)//GLY or unknown
   continue;

  pi=i;//for ca

  	x1 = ca[pi-2][0];
  	y1 = ca[pi-2][1];
  	z1 = ca[pi-2][2];
  	x2 = ca[pi-1][0];
  	y2 = ca[pi-1][1];
  	z2 = ca[pi-1][2];
  	x3 = ca[pi][0];
  	y3 = ca[pi][1];
  	z3 = ca[pi][2];
  	x4 = ca[pi+1][0];
  	y4 = ca[pi+1][1];
  	z4 = ca[pi+1][2];

	bin13_1 = rbins[pi][0];
      	bin13_2 = rbins[pi][1];
      	bin14 = rbins[pi][2];

	v1[0] = x4-x2;
      	v1[1] = y4-y2;
      	v1[2] = z4-z2;

      	v2a[0] = x4-x3;
      	v2a[1] = y4-y3;
      	v2a[2] = z4-z3;

      	v2b[0] = x3-x2;
      	v2b[1] = y3-y2;
      	v2b[2] = z3-z2;

	cross(v2a, v2b, v2);
      	cross(v1, v2, v3);

      	norm(v1);
      	norm(v2);
      	norm(v3);


	//init
	for (j=0;j<10;j++)
         SORTED_ROTAMERS[i][j][0] = 500.;

	j = 0;
      	besthit = 1000.;
      	bestpos = 0;

	//search 10 closest rotamer

	for(int j=idx_idx[seq[i]];rot_stat_idx[j][0]==seq[i];j++){
         if (rot_stat_idx[j][0]==seq[i]) {
          hit =  fabs(rot_stat_idx[j][1]-bin13_1) 
		+fabs(rot_stat_idx[j][2]-bin13_2) 
		+0.2*fabs(rot_stat_idx[j][3]-bin14);
          if (hit<SORTED_ROTAMERS[i][9][0]) {
            k = 9;
            while (k>=0 && hit<SORTED_ROTAMERS[i][k][0]) {
              k--;
            }
            k++;
            // k = hit
            for (l=9;l>k;l--) {
              SORTED_ROTAMERS[i][l][0]=SORTED_ROTAMERS[i][l-1][0];
              SORTED_ROTAMERS[i][l][1]=SORTED_ROTAMERS[i][l-1][1];
            }
            SORTED_ROTAMERS[i][k][0]=hit;//hit-score
            SORTED_ROTAMERS[i][k][1]=j;  //index
          }
         }else if(rot_stat_idx[j][0]>seq[i]){
	  break;
	 }
	}

	//puts("check 10 rot");
	//check 10 best rotemer

	for(rot=0;rot<10; rot++){

	 if(SORTED_ROTAMERS[i][rot][0]==500)
	  break;

	 besthit = SORTED_ROTAMERS[i][rot][0];
      	 bestpos = SORTED_ROTAMERS[i][rot][1];
	
	 pos = rot_stat_idx[bestpos][5];
      	 //nsc = nheavy[seq[i]]+1;
      	 nsc = nheavy[seq[i]];//ignore 0:C-alpha
	 
	 //if(nheavy[seq[i]]>0)
      	 //nsc = 1;//ignore 0:C-alpha


	//get XYZ system
	 for (j=0;j<3;j++) {
       	  vv[0][j] = v1[j]; vv[1][j] = v2[j]; vv[2][j] = v3[j];
       	 }

	 //copy from DB
	 //ignore 0: C-alpha

      	 for (j=0;j<nsc;j++) {
          sc[j][0] = rot_stat_coords[pos+j+1][0];
          sc[j][1] = rot_stat_coords[pos+j+1][1];
          sc[j][2] = rot_stat_coords[pos+j+1][2];
      	 }

	 //faster ver
	 JustImpose(vv,sc,nsc);

	 	//shift side-chains
		//printf("CA= %f %f %f\n",x3,y3,z3);
	 	for (j=0;j<nsc;j++) {
         	 sc[j][0] += x3;
         	 sc[j][1] += y3;
         	 sc[j][2] += z3;
		 //printf("SC= %f %f %f\n",sc[j][0],sc[j][1],sc[j][2]);
       	 	}
	 //check MAP potential score MDFF like
	 sco=MapScore2(sc,nsc,map);
	 SORTED_ROTAMERS[i][rot][2]=sco;
	 stbl[Ns].sco=sco;
	 stbl[Ns].pos=i;
	 stbl[Ns].idx=bestpos;
	 stbl[Ns].rot=rot;
	 Ns++;
	 //printf("%d 1 NS=%d\n",rot,Ns);

		for (j=0;j<nsc;j++) {
        	 slib[i][rot][j][0] =sc[j][0];
        	 slib[i][rot][j][1] =sc[j][1];
        	 slib[i][rot][j][2] =sc[j][2];
		}

	 //printf("%d 2 NS=%d\n",rot,Ns);
	}
 

 }

 //sort by sco
 //printf("Ns= %d\n",Ns);
 qsort(&stbl,Ns,sizeof(STBL),cmp_stbl);

 int Nclash=0;
 int Ntmp=0;//stmp
 Sscore=0;
 //Build side-chains
 for(int i=0;i<len;i++)
  out[i]=NULL;


 for(int i=0;i<Ns;i++){

  if(out[stbl[i].pos]!=NULL)
   continue;

  nsc = nheavy[seq[stbl[i].pos]];
  Nclash=ClashChk(slib[stbl[i].pos][stbl[i].rot],nsc,ca,len,stbl[i].pos,3.0);
  if(Nclash>1)
   continue;

  Nclash=ClashChk(slib[stbl[i].pos][stbl[i].rot],nsc,stmp,Ntmp,-1,3.0);
  if(Nclash>1)
   continue;

	//input
	for(int j=0;j<nsc;j++){
	 stmp[Ntmp][0]=slib[stbl[i].pos][stbl[i].rot][j][0];
	 stmp[Ntmp][1]=slib[stbl[i].pos][stbl[i].rot][j][1];
	 stmp[Ntmp][2]=slib[stbl[i].pos][stbl[i].rot][j][2];
	 Ntmp++;
	}
  out[stbl[i].pos]=slib[stbl[i].pos][stbl[i].rot];
  Sscore+=stbl[i].sco;
  //printf("%d %f pos= %d idx= %d rot= %d Nc= %d\n",i,stbl[i].sco,stbl[i].pos,stbl[i].idx,stbl[i].rot,Nclash);
 }
 //printf("SideScore= %f\n",Sscore);
 return Sscore;
}





int ClashChk(double **c1,int n1,double **c2,int n2,int ig,double cut){
 int i,j;
 int cnt=0;
 double x,y,z,d2;
 double cut2=cut*cut;
 for(i=0;i<n1;i++){
	for(j=0;j<n2;j++){
	 if(j==ig) continue;
	 x=c1[i][0]-c2[j][0];
	 y=c1[i][1]-c2[j][1];
	 z=c1[i][2]-c2[j][2];
	 d2=x*x+y*y+z*z;
	 if(d2<cut2)
	  cnt++;
	}
 }
 return cnt;
}

double MapScore(double **cd,int len,MRC *map){
 double sco=0;
 int i,ind,x,y,z;
 int xydim=map->xydim;

 double iv=map->iv_width;

/*
             dens - dt
MAP_SCORE=  ---------- = (dens - dt)*map->inv;
             dmax - dt
*/


 for(i=0;i<len;i++){

  x=(int)((cd[i][0]-map->orgxyz[0])*iv);
  y=(int)((cd[i][1]-map->orgxyz[1])*iv);
  z=(int)((cd[i][2]-map->orgxyz[2])*iv);

  //printf("xyz= %d %d %d\n",x,y,z);
  ind=xydim*z+map->xdim*y+x;
  //printf("ind= %d %f %f %f\n",ind,cd[i][0],cd[i][1],cd[i][2]);
  //sco+=(map->dens[ind]-map->map_t)*map->inv;
  sco+=map->dens[ind];
  //sco+=(map->dens[ind]);
 }
  //printf("sco= %f\n",sco);
 return sco;
}

double MapScore2(double cd[][3],int len,MRC *map){
 double sco=0;
 int i,ind,x,y,z;
 int xydim=map->xydim;

 double iv=map->iv_width;

/*
             dens - dt
MAP_SCORE=  ---------- = (dens - dt)*map->inv;
             dmax - dt
*/


 for(i=0;i<len;i++){

  x=(int)((cd[i][0]-map->orgxyz[0])*iv);
  y=(int)((cd[i][1]-map->orgxyz[1])*iv);
  z=(int)((cd[i][2]-map->orgxyz[2])*iv);

  if(x<0||y<0||z<0)
   continue;
  if(x>=map->xdim||y>=map->ydim||z>=map->zdim)
   continue;
  

  //printf("xyz= %d %d %d\n",x,y,z);
  ind=xydim*z+map->xdim*y+x;
  //sco+=(map->dens[ind]-map->map_t)*map->inv;
  sco+=map->dens[ind];//Already normalized
 }
  //printf("sco= %f\n",sco);
 return sco;
}



double rebuild_backbone_fromCA(double ca[MAX_FLEN][3],int *seq,int len, MRC *map,bool SHOW)
{

  //res_type *res, *prevres;
  //atom_type *atom;
  //real **cacoords, **tmpcoords, **tmpstat;
  real cacoords[MAX_FLEN][3], tmpcoords[MAX_FLEN][3], tmpstat[MAX_FLEN][3];
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real besthit, hit;
  int bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14, found, pro;
  int b13_1, b13_2, b14;
  real rmsd, total, maxrms;
  //FILE *debug, *out;
  int pi;
  int rbins[MAX_FLEN][3];
  //prepare_rbins();
  double sco=0;
  int atm=0;

  //prepare_rbins_ca(ca,rbins,len);

/*
    cacoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpcoords = (real**)calloc(sizeof(real*)*(8),1);
    tmpstat = (real**)calloc(sizeof(real*)*(8),1);
    for (i=0;i<8;i++) {
      cacoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpcoords[i] = (real*)calloc(sizeof(real)*3,1);;
      tmpstat[i] = (real*)calloc(sizeof(real)*3,1);;
    }
*/

    //prevres = NULL;
    //res = chain->residua;


    total = maxrms = 0.0;


     	//for (i=0;i<len+3;i++){
	//  printf("cao %d %f %f %f\n",i,ca[i][0],ca[i][1],ca[i][2]);
	//}

    for (i=0;i<len;i++) {
	//printf("%d\n",i);
	pi=i+2;

    	x1 = ca[pi-2][0];
    	y1 = ca[pi-2][1];
    	z1 = ca[pi-2][2];

    	x2 = ca[pi-1][0];
    	y2 = ca[pi-1][1];
    	z2 = ca[pi-1][2];

    	x3 = ca[pi][0];
    	y3 = ca[pi][1];
    	z3 = ca[pi][2];

    	x4 = ca[pi+1][0];
    	y4 = ca[pi+1][1];
    	z4 = ca[pi+1][2];

    	cacoords[0][0] = x1;
    	cacoords[0][1] = y1;
    	cacoords[0][2] = z1;

    	cacoords[1][0] = x2;
     	cacoords[1][1] = y2;
     	cacoords[1][2] = z2;

     	cacoords[2][0] = x3;
     	cacoords[2][1] = y3;
     	cacoords[2][2] = z3;

     	cacoords[3][0] = x4;
     	cacoords[3][1] = y4;
     	cacoords[3][2] = z4;

	 //printf("xyz1 %f %f %f\n",x1,y1,z1);


      bin13_1 = rbins[pi][0];
      bin13_2 = rbins[pi][1];
      bin14 = rbins[pi][2];

      pro = 0;

	//printf("%d 1\n",i);
      //if (prevres && !strncmp(prevres->name,"PRO",3)) {
      if (seq[i]==7) {//PRO
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat_pro[j].bins[0]-bin13_1)+fabs(nco_stat_pro[j].bins[1]-bin13_2)+0.2*fabs(nco_stat_pro[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat_pro[j].bins[0]>=0 && hit>1e-3);

        for (j=0;j<4;j++) {
         	for (k=0;k<3;k++) {
         		tmpstat[j][k] = nco_stat_pro[bestpos].data[j][k];
       	  	}
     	}

          for (j=0;j<8;j++) {
         	for (k=0;k<3;k++) {
         		tmpcoords[j][k] = nco_stat_pro[bestpos].data[j][k];
       	  	}
     	  }


	  for (k=0;k<3;k++){
	   //tmpcoords[0][k] = cacoords[2][k]; //Ca
	   tmpcoords[0][k] = nco_stat_pro[bestpos].data[4][k]; //C
	   tmpcoords[1][k] = nco_stat_pro[bestpos].data[5][k]; //O
	   tmpcoords[2][k] = nco_stat_pro[bestpos].data[6][k]; //N
	  }

      } else {
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat[j].bins[0]-bin13_1)+fabs(nco_stat[j].bins[1]-bin13_2)+0.2*fabs(nco_stat[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat[j].bins[0]>=0 && hit>1e-3);

        for (j=0;j<4;j++) {
         	for (k=0;k<3;k++) {
         		tmpstat[j][k] = nco_stat[bestpos].data[j][k];
       	  	}
     	}

/*
          for (j=0;j<8;j++) {
         	for (k=0;k<3;k++) {
         		tmpcoords[j][k] = nco_stat[bestpos].data[j][k];
       	  	}
     	  }
*/

	  //Back-bone
	  for (k=0;k<3;k++){
	   //tmpcoords[0][k] = cacoords[2][k]; //Ca
	   tmpcoords[0][k] = nco_stat[bestpos].data[4][k]; //C
	   tmpcoords[1][k] = nco_stat[bestpos].data[5][k]; //O
	   tmpcoords[2][k] = nco_stat[bestpos].data[6][k]; //N
	  } 

      }
/*
	for (k=0;k<3;k++){
	 printf("1ca %f %f %f\n",cacoords[k][0],cacoords[k][1],cacoords[k][1]);
	 printf("1tmp %f %f %f\n",tmpstat[k][0],tmpstat[k][1],tmpstat[k][1]);
	}
*/
	//printf("%d 2\n",i);
	//printf("%f %f %f\n",tmpcoords[0][0],tmpcoords[0][1],tmpcoords[0][2]);
	//printf("%f %f %f\n",cacoords[0][0],cacoords[0][1],cacoords[0][2]);
     	//rmsd=superimpose2(cacoords, tmpstat, 4, tmpcoords, 8);
	//Only Ca,C,O,N
     	rmsd=superimpose2(cacoords, tmpstat, 4, tmpcoords, 3);

/*
	for (k=0;k<3;k++){
	 printf("2ca %f %f %f\n",cacoords[k][0],cacoords[k][1],cacoords[k][1]);
	 printf("2tmp %f %f %f\n",tmpstat[k][0],tmpstat[k][1],tmpstat[k][1]);
	}
*/
	//Add Original ca
	for (k=0;k<3;k++)
	   tmpcoords[3][k] = cacoords[2][k]; //Ca

	//printf("%d 3\n",i);
     	total += rmsd;
     	if (rmsd>maxrms) maxrms=rmsd;


	//0,1,2,3 : 4-Calpha atoms
	//4,5,6 : C, O, N
/*
// add-or-replace

      if (prevres) {
        add_replace(prevres, "C  ", tmpcoords[4][0], tmpcoords[4][1], tmpcoords[4][2], FLAG_BACKBONE);
        add_replace(prevres, "O  ", tmpcoords[5][0], tmpcoords[5][1], tmpcoords[5][2], FLAG_BACKBONE);
      }

      if (res) {
        add_replace(res, "N  ", tmpcoords[6][0], tmpcoords[6][1], tmpcoords[6][2], FLAG_BACKBONE);
      } else { // terminal oxygen instead of nitrogen
        add_replace(prevres, "OXT", tmpcoords[6][0], tmpcoords[6][1], tmpcoords[6][2], FLAG_BACKBONE);
      }

      prevres = res;
      if (res)
        res = res->next;
*/

	//printf("%d 4 rms= %f\n",i,rmsd);

	//sco+=MapScore(tmpcoords,4,map);

	 //printf("seq= %d\n",seq[i]);

	 if(SHOW){
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "N  ", RES_NAMES[seq[i]], ' ',i+1,tmpcoords[2][0],tmpcoords[2][1], tmpcoords[2][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "CA ", RES_NAMES[seq[i]], ' ',i+1,tmpcoords[3][0],tmpcoords[3][1], tmpcoords[3][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "C  ", RES_NAMES[seq[i]], ' ',i+1,tmpcoords[0][0],tmpcoords[0][1], tmpcoords[0][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "O  ", RES_NAMES[seq[i]], ' ',i+1,tmpcoords[1][0],tmpcoords[1][1], tmpcoords[1][2]);
	 }
    }
 if(SHOW)
  printf("ENDMDL\n");
    //printf("Backbone rebuilding deviation: average = %.3f, max = %.3f\n", total/(real)chain_length, maxrms);
 return sco;
}

//output atom coords
double rebuild_backbone_fromCA2(double **ca,double **out,int **rbins,int *seq,int len)
{

  real cacoords[4][3], tmpcoords[4][3], tmpstat[4][3];
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real besthit, hit;
  int bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14, found, pro;
  int b13_1, b13_2, b14;
  real rmsd, total, maxrms;
  int pi;
  //int rbins[MAX_FLEN][3];
  double sco=0;
  int atm=0;

  prepare_rbins_ca(ca,rbins,len-4);

  //puts("OK");
    total = maxrms = 0.0;

    //for (i=0;i<len+1;i++) {
    for (i=2;i<len-1;i++) {
	//pi=i+2;
	pi=i;

    	x1 = ca[pi-2][0];
    	y1 = ca[pi-2][1];
    	z1 = ca[pi-2][2];

    	x2 = ca[pi-1][0];
    	y2 = ca[pi-1][1];
    	z2 = ca[pi-1][2];

    	x3 = ca[pi][0];
    	y3 = ca[pi][1];
    	z3 = ca[pi][2];

    	x4 = ca[pi+1][0];
    	y4 = ca[pi+1][1];
    	z4 = ca[pi+1][2];

    	cacoords[0][0] = x1; cacoords[0][1] = y1; cacoords[0][2] = z1;

    	cacoords[1][0] = x2; cacoords[1][1] = y2; cacoords[1][2] = z2;

     	cacoords[2][0] = x3; cacoords[2][1] = y3; cacoords[2][2] = z3;

     	cacoords[3][0] = x4; cacoords[3][1] = y4; cacoords[3][2] = z4;

      bin13_1 = rbins[pi][0];
      bin13_2 = rbins[pi][1];
      bin14 = rbins[pi][2];

      pro = 0;

      if (seq[i]==7) {//PRO
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat_pro[j].bins[0]-bin13_1)+fabs(nco_stat_pro[j].bins[1]-bin13_2)+0.2*fabs(nco_stat_pro[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat_pro[j].bins[0]>=0 && hit>1e-3);

        for (j=0;j<4;j++) { //Ca positions
       	 tmpstat[j][0] = nco_stat_pro[bestpos].data[j][0];
       	 tmpstat[j][1] = nco_stat_pro[bestpos].data[j][1];
       	 tmpstat[j][2] = nco_stat_pro[bestpos].data[j][2];
     	}


/*
          for (j=0;j<8;j++) {
         	tmpcoords[j][0] = nco_stat_pro[bestpos].data[j][0];
         	tmpcoords[j][1] = nco_stat_pro[bestpos].data[j][1];
         	tmpcoords[j][2] = nco_stat_pro[bestpos].data[j][2];
     	  }
*/

	  for (k=0;k<3;k++){
	   //tmpcoords[0][k] = cacoords[2][k]; //Ca
	   tmpcoords[0][k] = nco_stat_pro[bestpos].data[6][k]; //N
	   tmpcoords[1][k] = nco_stat_pro[bestpos].data[4][k]; //C
	   tmpcoords[2][k] = nco_stat_pro[bestpos].data[5][k]; //O
	   //tmpcoords[2][k] = nco_stat_pro[bestpos].data[7][k]; //Cbeta
	  }

      } else {
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat[j].bins[0]-bin13_1)+fabs(nco_stat[j].bins[1]-bin13_2)+0.2*fabs(nco_stat[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat[j].bins[0]>=0 && hit>1e-3);

        for (j=0;j<4;j++) {
        	tmpstat[j][0] = nco_stat[bestpos].data[j][0];
         	tmpstat[j][1] = nco_stat[bestpos].data[j][1];
         	tmpstat[j][2] = nco_stat[bestpos].data[j][2];
     	}

	  //Back-bone
	  for (k=0;k<3;k++){
	   //tmpcoords[0][k] = cacoords[2][k]; //Ca
	   tmpcoords[0][k] = nco_stat[bestpos].data[6][k]; //N
	   tmpcoords[1][k] = nco_stat[bestpos].data[4][k]; //C
	   tmpcoords[2][k] = nco_stat[bestpos].data[5][k]; //O
	   //tmpcoords[2][k] = nco_stat[bestpos].data[7][k]; //Cbeta
	  } 

      }
	//Only Ca,C,O,N
     	rmsd=superimpose2(cacoords, tmpstat, 4, tmpcoords, 3);

/*
	for (k=0;k<3;k++){
	 printf("2ca %f %f %f\n",cacoords[k][0],cacoords[k][1],cacoords[k][1]);
	 printf("2tmp %f %f %f\n",tmpstat[k][0],tmpstat[k][1],tmpstat[k][1]);
	}
*/
	//Add Original ca
	  tmpcoords[3][0] = cacoords[2][0]; //Ca
	  tmpcoords[3][1] = cacoords[2][1]; //Ca
	  tmpcoords[3][2] = cacoords[2][2]; //Ca

	//printf("%d 3\n",i);
     	total += rmsd;
     	if (rmsd>maxrms) maxrms=rmsd;


	//0,1,2,3 : 4-Calpha atoms
	//4,5,6 : C, O, N

	//copy

	//N
	out[4*i][0]=tmpcoords[0][0]; out[4*i][1]=tmpcoords[0][1]; out[4*i][2]=tmpcoords[0][2];
	//CA
	out[4*i+1][0]=tmpcoords[3][0]; out[4*i+1][1]=tmpcoords[3][1]; out[4*i+1][2]=tmpcoords[3][2];
	//C
	out[4*(i-1)+2][0]=tmpcoords[1][0]; out[4*(i-1)+2][1]=tmpcoords[1][1]; out[4*(i-1)+2][2]=tmpcoords[1][2];
	//O
	out[4*(i-1)+3][0]=tmpcoords[2][0]; out[4*(i-1)+3][1]=tmpcoords[2][1]; out[4*(i-1)+3][2]=tmpcoords[2][2];
    }
// if(SHOW)
//  printf("ENDMDL\n");
    //printf("Backbone rebuilding deviation: average = %.3f, max = %.3f\n", total/(real)chain_length, maxrms);
 return sco;
}

double rebuild_backbone_fromCA3(double **ca,double **out,int **rbins,int *seq,int len)
{

  real cacoords[4][3], tmpcoords[4][3], tmpstat[4][3];
  real x1, y1, z1;
  real x2, y2, z2;
  real x3, y3, z3;
  real x4, y4, z4;
  real besthit, hit;
  int bestpos;
  int i, j, k, l, m, bin13_1, bin13_2, bin14, found, pro;
  int b13_1, b13_2, b14;
  real rmsd, total, maxrms;
  int pi;
  //int rbins[MAX_FLEN][3];
  double sco=0;
  int atm=0;

  prepare_rbins_ca(ca,rbins,len);

    total = maxrms = 0.0;

    for (i=0;i<len+1;i++) {
	pi=i+2;

    	x1 = ca[pi-2][0];
    	y1 = ca[pi-2][1];
    	z1 = ca[pi-2][2];

    	x2 = ca[pi-1][0];
    	y2 = ca[pi-1][1];
    	z2 = ca[pi-1][2];

    	x3 = ca[pi][0];
    	y3 = ca[pi][1];
    	z3 = ca[pi][2];

    	x4 = ca[pi+1][0];
    	y4 = ca[pi+1][1];
    	z4 = ca[pi+1][2];

    	cacoords[0][0] = x1; cacoords[0][1] = y1; cacoords[0][2] = z1;

    	cacoords[1][0] = x2; cacoords[1][1] = y2; cacoords[1][2] = z2;

     	cacoords[2][0] = x3; cacoords[2][1] = y3; cacoords[2][2] = z3;

     	cacoords[3][0] = x4; cacoords[3][1] = y4; cacoords[3][2] = z4;

      bin13_1 = rbins[pi][0];
      bin13_2 = rbins[pi][1];
      bin14 = rbins[pi][2];

      pro = 0;

      if (seq[i]==7) {//PRO
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat_pro[j].bins[0]-bin13_1)+fabs(nco_stat_pro[j].bins[1]-bin13_2)+0.2*fabs(nco_stat_pro[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat_pro[j].bins[0]>=0 && hit>1e-3);

        for (j=0;j<4;j++) { //Ca positions
       	 tmpstat[j][0] = nco_stat_pro[bestpos].data[j][0];
       	 tmpstat[j][1] = nco_stat_pro[bestpos].data[j][1];
       	 tmpstat[j][2] = nco_stat_pro[bestpos].data[j][2];
     	}


/*
          for (j=0;j<8;j++) {
         	tmpcoords[j][0] = nco_stat_pro[bestpos].data[j][0];
         	tmpcoords[j][1] = nco_stat_pro[bestpos].data[j][1];
         	tmpcoords[j][2] = nco_stat_pro[bestpos].data[j][2];
     	  }
*/

	  for (k=0;k<3;k++){
	   //tmpcoords[0][k] = cacoords[2][k]; //Ca
	   tmpcoords[0][k] = nco_stat_pro[bestpos].data[6][k]; //N
	   tmpcoords[1][k] = nco_stat_pro[bestpos].data[4][k]; //C
	   tmpcoords[2][k] = nco_stat_pro[bestpos].data[5][k]; //O
	   //tmpcoords[2][k] = nco_stat_pro[bestpos].data[7][k]; //Cbeta
	  }

      } else {
        j=0;
        besthit=1000.;
        bestpos=0;
        do {
          hit = fabs(nco_stat[j].bins[0]-bin13_1)+fabs(nco_stat[j].bins[1]-bin13_2)+0.2*fabs(nco_stat[j].bins[2]-bin14);
          if (hit<besthit) {
            besthit=hit;
            bestpos=j;
          }
          j++;
        } while (nco_stat[j].bins[0]>=0 && hit>1e-3);

        for (j=0;j<4;j++) {
        	tmpstat[j][0] = nco_stat[bestpos].data[j][0];
         	tmpstat[j][1] = nco_stat[bestpos].data[j][1];
         	tmpstat[j][2] = nco_stat[bestpos].data[j][2];
     	}

	  //Back-bone
	  for (k=0;k<3;k++){
	   //tmpcoords[0][k] = cacoords[2][k]; //Ca
	   tmpcoords[0][k] = nco_stat[bestpos].data[6][k]; //N
	   tmpcoords[1][k] = nco_stat[bestpos].data[4][k]; //C
	   tmpcoords[2][k] = nco_stat[bestpos].data[5][k]; //O
	   //tmpcoords[2][k] = nco_stat[bestpos].data[7][k]; //Cbeta
	  } 

      }
	//Only Ca,C,O,N
     	rmsd=superimpose2(cacoords, tmpstat, 4, tmpcoords, 3);

/*
	for (k=0;k<3;k++){
	 printf("2ca %f %f %f\n",cacoords[k][0],cacoords[k][1],cacoords[k][1]);
	 printf("2tmp %f %f %f\n",tmpstat[k][0],tmpstat[k][1],tmpstat[k][1]);
	}
*/
	//Add Original ca
	  tmpcoords[3][0] = cacoords[2][0]; //Ca
	  tmpcoords[3][1] = cacoords[2][1]; //Ca
	  tmpcoords[3][2] = cacoords[2][2]; //Ca

	//printf("%d 3\n",i);
     	total += rmsd;
     	if (rmsd>maxrms) maxrms=rmsd;

  //puts("OK");

	//0,1,2,3 : 4-Calpha atoms
	//4,5,6 : C, O, N

	//copy

	//N
	out[4*i][0]=tmpcoords[0][0]; out[4*i][1]=tmpcoords[0][1]; out[4*i][2]=tmpcoords[0][2];
	//CA
	out[4*i+1][0]=tmpcoords[3][0]; out[4*i+1][1]=tmpcoords[3][1]; out[4*i+1][2]=tmpcoords[3][2];
	if(i>0){
	 //C
	 out[4*(i-1)+2][0]=tmpcoords[1][0]; out[4*(i-1)+2][1]=tmpcoords[1][1]; out[4*(i-1)+2][2]=tmpcoords[1][2];
	 //O
	 out[4*(i-1)+3][0]=tmpcoords[2][0]; out[4*(i-1)+3][1]=tmpcoords[2][1]; out[4*(i-1)+3][2]=tmpcoords[2][2];
	}


	//sco+=MapScore(tmpcoords,4,map);

	 //printf("seq= %d\n",seq[i]);

	 if(false){
	 if(i>0){
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "C  ", RES_NAMES[seq[i-1]], ' ',i,tmpcoords[1][0],tmpcoords[1][1], tmpcoords[1][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "O  ", RES_NAMES[seq[i-1]], ' ',i,tmpcoords[2][0],tmpcoords[2][1], tmpcoords[2][2]);
	 }
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "N  ", RES_NAMES[seq[i]], ' ',i+1,tmpcoords[0][0],tmpcoords[0][1], tmpcoords[0][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "CA ", RES_NAMES[seq[i]], ' ',i+1,tmpcoords[3][0],tmpcoords[3][1], tmpcoords[3][2]);
	 }
    }
// if(SHOW)
//  printf("ENDMDL\n");
    //printf("Backbone rebuilding deviation: average = %.3f, max = %.3f\n", total/(real)chain_length, maxrms);
 return sco;
}

//add mid-points
void rebuild_psbackbone_fromCA(double **ca,double **out,int len)
{
 double v[3];
 int i;
 for(i=0;i<len-1;i++){
  v[0]=0.5*(ca[i][0]+ca[i+1][0]);
  v[1]=0.5*(ca[i][1]+ca[i+1][1]);
  v[2]=0.5*(ca[i][2]+ca[i+1][2]);

  out[2*i][0]=ca[i][0]; out[2*i][1]=ca[i][1]; out[2*i][2]=ca[i][2];//Original
  out[2*i+1][0]=v[0]; out[2*i+1][1]=v[1]; out[2*i+1][2]=v[2];//mid-point
 }
  out[2*i][0]=ca[i][0];
  out[2*i][1]=ca[i][1];
  out[2*i][2]=ca[i][2];
}


double hbond_score(double **cd,int *ss,int len){
 int sco=0;
 int i,j,n1,o1,n2,o2;
 float d_on,d_no;
 //HBOND_D2
 for(i=0;i<len-1;i++){
  n1=4*i;
  o1=4*i+3;
 	for(j=i+4;j<len;j++){
  	 n2=4*j;
  	 o2=4*j+3;

	 d_no=(cd[n1][0]-cd[o2][0])*(cd[n1][0]-cd[o2][0])+
	      (cd[n1][1]-cd[o2][1])*(cd[n1][1]-cd[o2][1])+
	      (cd[n1][2]-cd[o2][2])*(cd[n1][2]-cd[o2][2]);

	 d_on=(cd[o1][0]-cd[n2][0])*(cd[o1][0]-cd[n2][0])+
	      (cd[o1][1]-cd[n2][1])*(cd[o1][1]-cd[n2][1])+
	      (cd[o1][2]-cd[n2][2])*(cd[o1][2]-cd[n2][2]);

	 //Any H-bond
	 if(d_no<HBOND_D2){
	  sco++;
	  //beta
	  if(ss[i]==2 && ss[j]==2)
	   sco++;
	 }
	 if(d_on<HBOND_D2){

	  sco++;
	  //Helix
	  if(ss[i]==1 && ss[j]==1 && i+4==j)
	   sco++;
	  else if(ss[i]==2 && ss[j]==2)
	   sco++;
	 }
	}
 }
 return (double)sco;
}


#define HB_H_R 6.15
#define HB_H_S 0.53
#define LC1_R 3.80
#define LC1_S 0.20
#define LC2_H_R 5.55
#define LC2_H_S 0.55
#define LC3_H_R 5.10 //New
#define LC3_H_S 0.55 //New
#define LC2_E_R 6.60
#define LC2_E_S 0.80

int six_vec[6][3]={
{0 , 0,-1},
{0 , 0, 1},
{0 ,-1, 0},
{0 , 1, 0},
{-1, 0, 0},
{ 1, 0, 0}
};


double GetVec(double **in,double **init,double **vec,int n,int *ss){
 //double vec[MAX_FLEN][3];
 int p;
 double rd,d2,d,v[3],f;
 double sum=0;
 //init
 for(int i=0;i<n;i++)
  vec[i][0]=vec[i][1]=vec[i][2]=0.0;


 //initial position in 3rd CA
 v[0]=init[2][0]-in[2][0];
 v[1]=init[2][1]-in[2][1];
 v[2]=init[2][2]-in[2][2];
 d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
 if(d>2.0){
  f=(d-1.0)/d;
  vec[2][0]+=v[0]*f;
  vec[2][1]+=v[1]*f;
  vec[2][2]+=v[2]*f;
 }

 //---------------------------

 for(int i=0;i<n;i++){
  for(int j=i+1;j<n;j++){
   //vec i->j
   v[0]=in[j][0]-in[i][0];
   v[1]=in[j][1]-in[i][1];
   v[2]=in[j][2]-in[i][2];
   d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
   rd=1.000/d;
	//+1
	if(j-i==1){
	 if(fabs(d-LC1_R)>LC1_S){
	  f=(d-LC1_R)/d;

	  vec[i][0]+=v[0]*f;
	  vec[i][1]+=v[1]*f;
	  vec[i][2]+=v[2]*f;

	  vec[j][0]-=v[0]*f;
	  vec[j][1]-=v[1]*f;
	  vec[j][2]-=v[2]*f;
	 }
	 continue;
	}
	if(j-i==2){ //+2
	 if(ss[i]==1 && ss[i+1]==1 && ss[j]==1){
	  if(fabs(d-LC2_H_R)>LC2_H_S){//H-H
	   f=(d-LC2_H_R)/d;
	   vec[i][0]+=v[0]*f;
	   vec[i][1]+=v[1]*f;
	   vec[i][2]+=v[2]*f;
	
	   vec[j][0]-=v[0]*f;
	   vec[j][1]-=v[1]*f;
	   vec[j][2]-=v[2]*f;

	  }
	  continue;
	 }else if(ss[i]==2 && ss[i+1]==2 && ss[j]==2){
	  if(fabs(d-LC2_E_R)>LC2_E_S){//E-E
	   f=(d-LC2_E_R)/d;
	   vec[i][0]+=v[0]*f;
	   vec[i][1]+=v[1]*f;
	   vec[i][2]+=v[2]*f;

	   vec[j][0]-=v[0]*f;
	   vec[j][1]-=v[1]*f;
	   vec[j][2]-=v[2]*f;

	  }
	  continue;
	 }else{
	  if(d-LC2_E_R>LC2_E_S){//too far
	   f=(d-LC2_E_R)/d;
	   vec[i][0]+=v[0]*f;
	   vec[i][1]+=v[1]*f;
	   vec[i][2]+=v[2]*f;

	   vec[j][0]-=v[0]*f;
	   vec[j][1]-=v[1]*f;
	   vec[j][2]-=v[2]*f;

	  }else if(LC2_H_R-d > LC2_H_S){//too close
	   f=(d-LC2_H_R)/d; //negative
	   vec[i][0]+=v[0]*f;
	   vec[i][1]+=v[1]*f;
	   vec[i][2]+=v[2]*f;

	   vec[j][0]-=v[0]*f;
	   vec[j][1]-=v[1]*f;
	   vec[j][2]-=v[2]*f;
	   }
	  }
	  continue;
	 }
	 if(j-i==3){ //+3
 	  //H-H
	  if(i<n-3 && ss[i+3]==1 && ss[i+2]==1 && ss[i+1]==1 && ss[i]==1 ){
	   if(fabs(d-LC3_H_R)>LC3_H_S){
	    f=(d-LC3_H_R)/d;
	    vec[i][0]+=v[0]*f;
	    vec[i][1]+=v[1]*f;
	    vec[i][2]+=v[2]*f;
	
	    vec[j][0]-=v[0]*f;
	    vec[j][1]-=v[1]*f;
	    vec[j][2]-=v[2]*f;
	   }
	   continue;
     	  } 
   	 }
	 if(j-i==4){
	  //H-H
	  if(ss[i+4]==1 && ss[i+3]==1 && ss[i+2]==1 && ss[i+1]==1 &&ss[i]==1 ){
	   if(fabs(d-HB_H_R)>HB_H_S){
	    f=(d-HB_H_R)/d;

	    vec[i][0]+=v[0]*f;
	    vec[i][1]+=v[1]*f;
	    vec[i][2]+=v[2]*f;

	    vec[j][0]-=v[0]*f;
	    vec[j][1]-=v[1]*f;
	    vec[j][2]-=v[2]*f;

	   }
	   continue;
	  }
	 }
	 //avoid all clashes
	 if(LC3_H_R-d>LC3_H_S){//too close
	  f=(d-LC3_H_R)/d;
	  vec[i][0]+=v[0]*f;
	  vec[i][1]+=v[1]*f;
	  vec[i][2]+=v[2]*f;

	  vec[j][0]-=v[0]*f;
	  vec[j][1]-=v[1]*f;
	  vec[j][2]-=v[2]*f;

	}



  }
 }

/*
 for(int i=0;i<n;i++){

  //local-1
   if(i>0){
    v[0]=in[i-1][0]-in[i][0];
    v[1]=in[i-1][1]-in[i][1];
    v[2]=in[i-1][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(fabs(d-LC1_R)>LC1_S){
	 f=(d-LC1_R)/d;
	 vec[i][0]+=v[0]*f;
	 vec[i][1]+=v[1]*f;
	 vec[i][2]+=v[2]*f;
	}
   }
  //local +1
  if(i<n-1){
    v[0]=in[i+1][0]-in[i][0];
    v[1]=in[i+1][1]-in[i][1];
    v[2]=in[i+1][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(fabs(d-LC1_R)>LC1_S){
	 f=(d-LC1_R)/d;
	 vec[i][0]+=v[0]*f;
	 vec[i][1]+=v[1]*f;
	 vec[i][2]+=v[2]*f;
	}
   }

  //printf("%d vec= %f %f %f ss= %d\n",i,vec[i][0],vec[i][1],vec[i][2],ss[i]);
  if(ss[i]==1){//Helix

   //HB H -4
   if(i>3 && ss[i-4]==1 && ss[i-3]==1 && ss[i-2]==1 && ss[i-1]==1){
    v[0]=in[i-4][0]-in[i][0];
    v[1]=in[i-4][1]-in[i][1];
    v[2]=in[i-4][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(fabs(d-HB_H_R)>HB_H_S){
	 f=(d-HB_H_R)/d;

	 vec[i][0]+=v[0]*f;
	 vec[i][1]+=v[1]*f;
	 vec[i][2]+=v[2]*f;
	}
    //printf("%d -4 vec= %f %f %f d= %f\n",i,vec[i][0],vec[i][1],vec[i][2],d);
   }

   //HB H +4
   if(i<n-4 && ss[i+4]==1 && ss[i+3]==1 && ss[i+2]==1 && ss[i+1]==1){
    v[0]=in[i+4][0]-in[i][0];
    v[1]=in[i+4][1]-in[i][1];
    v[2]=in[i+4][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(fabs(d-HB_H_R)>HB_H_S){
	 f=(d-HB_H_R)/d;

	 vec[i][0]+=v[0]*f;
	 vec[i][1]+=v[1]*f;
	 vec[i][2]+=v[2]*f;
	}
    //printf("%d +4 vec= %f %f %f d= %f\n",i,vec[i][0],vec[i][1],vec[i][2],d);
   }


   //New!!!
   //local H -3
   if(i>2 && ss[i-3]==1 && ss[i-2]==1 && ss[i-1]==1){
    v[0]=in[i-3][0]-in[i][0];
    v[1]=in[i-3][1]-in[i][1];
    v[2]=in[i-3][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(fabs(d-LC3_H_R)>LC3_H_S){
	 f=(d-LC3_H_R)/d;
	 vec[i][0]+=v[0]*f;
	 vec[i][1]+=v[1]*f;
	 vec[i][2]+=v[2]*f;
	}
    //printf("%d +2 vec= %f %f %f d= %f\n",i,vec[i][0],vec[i][1],vec[i][2],d);
   }
   //local H +3
   if(i<n-3 && ss[i+3]==1 && ss[i+2]==1 && ss[i+1]==1){
    v[0]=in[i+3][0]-in[i][0];
    v[1]=in[i+3][1]-in[i][1];
    v[2]=in[i+3][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(fabs(d-LC3_H_R)>LC3_H_S){
	 f=(d-LC3_H_R)/d;
	 vec[i][0]+=v[0]*f;
	 vec[i][1]+=v[1]*f;
	 vec[i][2]+=v[2]*f;
	}
    //printf("%d -2 vec= %f %f %f d= %f\n",i,vec[i][0],vec[i][1],vec[i][2],d);
   }



   //local+-2
   if(i>1){
    v[0]=in[i-2][0]-in[i][0];
    v[1]=in[i-2][1]-in[i][1];
    v[2]=in[i-2][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
    	if(ss[i-2]==1 && ss[i-1]==1){//H-H
 	 if(fabs(d-LC2_H_R)>LC2_H_S){
	  f=(d-LC2_H_R)/d;
	  vec[i][0]+=v[0]*f;
	  vec[i][1]+=v[1]*f;
	  vec[i][2]+=v[2]*f;
	 }
	}else{//H-E and H-C
	 if(LC2_H_R-d > LC2_H_S){//too close
	  f=(d-LC2_H_R)/d; //negative
	  vec[i][0]+=v[0]*f;
	  vec[i][1]+=v[1]*f;
	  vec[i][2]+=v[2]*f;
	 }
	}
    //printf("%d +2 vec= %f %f %f d= %f\n",i,vec[i][0],vec[i][1],vec[i][2],d);
   }
   if(i<n-2){
    v[0]=in[i+2][0]-in[i][0];
    v[1]=in[i+2][1]-in[i][1];
    v[2]=in[i+2][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(ss[i+2]==1 && ss[i+1]==1){
	 if(fabs(d-LC2_H_R)>LC2_H_S){//H-H
	  f=(d-LC2_H_R)/d;
	  vec[i][0]+=v[0]*f;
	  vec[i][1]+=v[1]*f;
	  vec[i][2]+=v[2]*f;
	 }
	}else{//H-E and H-C
	 if(LC2_H_R-d > LC2_H_S){//too close
	  f=(d-LC2_H_R)/d; //negative
	  vec[i][0]+=v[0]*f;
	  vec[i][1]+=v[1]*f;
	  vec[i][2]+=v[2]*f;
	 }
	}
    //printf("%d -2 vec= %f %f %f d= %f\n",i,vec[i][0],vec[i][1],vec[i][2],d);
   }


  }else if(ss[i]==2){//Seet
   //local+-2
   if(i>1){
    v[0]=in[i-2][0]-in[i][0];
    v[1]=in[i-2][1]-in[i][1];
    v[2]=in[i-2][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if( ss[i-2]==2 && ss[i-1]==2){
	 if(fabs(d-LC2_E_R)>LC2_E_S){
	  f=(d-LC2_E_R)/d;
	  vec[i][0]+=v[0]*f;
	  vec[i][1]+=v[1]*f;
	  vec[i][2]+=v[2]*f;
	 }
	}else{
	 if(LC2_H_R-d > LC2_H_S){//too close
	  f=(d-LC2_H_R)/d; //negative
	  vec[i][0]+=v[0]*f;
	  vec[i][1]+=v[1]*f;
	  vec[i][2]+=v[2]*f;
	 }
	}
   }
   if(i<n-2){
    v[0]=in[i+2][0]-in[i][0];
    v[1]=in[i+2][1]-in[i][1];
    v[2]=in[i+2][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(ss[i+2]==2 && ss[i+1]==2){
	 if(fabs(d-LC2_E_R)>LC2_E_S){
	  f=(d-LC2_E_R)/d;
	  vec[i][0]+=v[0]*f;
	  vec[i][1]+=v[1]*f;
	  vec[i][2]+=v[2]*f;
	 }
	}else{
	  if(LC2_H_R-d > LC2_H_S){//too close
	  f=(d-LC2_H_R)/d; //negative
	  vec[i][0]+=v[0]*f;
	  vec[i][1]+=v[1]*f;
	  vec[i][2]+=v[2]*f;
	 }
	}
   }

  }else{//Coil
   //local+-2
   if(i>1){
    v[0]=in[i-2][0]-in[i][0];
    v[1]=in[i-2][1]-in[i][1];
    v[2]=in[i-2][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(d-LC2_E_R>LC2_E_S){//too far
	 f=(d-LC2_E_R)/d;
	 vec[i][0]+=v[0]*f;
	 vec[i][1]+=v[1]*f;
	 vec[i][2]+=v[2]*f;
	}else if(LC2_H_R-d>LC2_H_S){//too close
	 f=(d-LC2_H_R)/d;
	 vec[i][0]+=v[0]*f;
	 vec[i][1]+=v[1]*f;
	 vec[i][2]+=v[2]*f;
	}
   }
   if(i<n-2){
    v[0]=in[i+2][0]-in[i][0];
    v[1]=in[i+2][1]-in[i][1];
    v[2]=in[i+2][2]-in[i][2];
    d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
	if(d-LC2_E_R>LC2_E_S){//too far
	 f=(d-LC2_E_R)/d;
	 vec[i][0]+=v[0]*f;
	 vec[i][1]+=v[1]*f;
	 vec[i][2]+=v[2]*f;
	}else if(LC2_H_R-d>LC2_H_S){//too close
	 f=(d-LC2_H_R)/d;
	 vec[i][0]+=v[0]*f;
	 vec[i][1]+=v[1]*f;
	 vec[i][2]+=v[2]*f;
	}
   }


  }

  //Normarized by Max vec

  d=sqrt(vec[i][0]*vec[i][0]+vec[i][1]*vec[i][1]+vec[i][2]*vec[i][2]);
  if(d>1.00){
   vec[i][0]/=d;
   vec[i][1]/=d;
   vec[i][2]/=d;
   //printf("*%d vec= %f %f %f ss= %d d= %f\n",i,vec[i][0],vec[i][1],vec[i][2],ss[i],d);
  }
  //printf("%d vec= %f %f %f ss= %d d= %f\n",i,vec[i][0],vec[i][1],vec[i][2],ss[i],d);
  sum+=d;
 }
*/
 for(int i=0;i<n;i++){
  d=sqrt(vec[i][0]*vec[i][0]+vec[i][1]*vec[i][1]+vec[i][2]*vec[i][2]);
  if(d>1.00){
   vec[i][0]/=d;
   vec[i][1]/=d;
   vec[i][2]/=d;
  }
  sum+=d;
 }
 //puts("End vec");
 return sum;
}

double GetVecMain(double **in,double **ca,double **vec,int n,int *seq,int pos,MRC *map,double w){
 //double vec[MAX_FLEN][3];
 int p;
 double d2,d,v[3],f;
 double sum=0;
 double r1=1.0/(1.0+w);
 double r2=w/(1.0+w);

 //density based vec
 double sco,sco2;
 int i,ind,x,y,z;
 int xydim=map->xydim;

 double iv=map->iv_width;

/*
             dens - dt
MAP_SCORE=  ---------- = (dens - dt)*map->inv;
             dmax - dt
*/

 //Main chain positions
 for(int I=0;I<n;I++){
  v[0]=v[1]=v[2]=0.0;
  for(int atm=0;atm<4;atm++){
  i=I*4+atm;//N, CA, C, O
  p=pos+I;
  x=(int)((in[i][0]-map->orgxyz[0])*iv);
  y=(int)((in[i][1]-map->orgxyz[1])*iv);
  z=(int)((in[i][2]-map->orgxyz[2])*iv);
  ind=xydim*z+map->xdim*y+x;
  sco=(map->dens[ind]-map->map_t)*map->inv;
	//6 vec
	for(int j=0;j<6;j++){
	 ind=xydim*(z+six_vec[j][2])+map->xdim*(y+six_vec[j][1])+(x+six_vec[j][0]);
	 sco2=(map->dens[ind]-map->map_t)*map->inv;
	 if(sco>sco2)
	  continue;
	 v[0]+=(sco2-sco)*six_vec[j][0];
	 v[1]+=(sco2-sco)*six_vec[j][1];
	 v[2]+=(sco2-sco)*six_vec[j][2];
	}
  }
  //unit vec

  d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if(d>0){
   vec[I+2][0]=r1*vec[I+2][0]+r2*v[0]/d;
   vec[I+2][1]=r1*vec[I+2][1]+r2*v[1]/d;
   vec[I+2][2]=r1*vec[I+2][2]+r2*v[2]/d;
  }
  sum+=d;
 }

 //for N and C ter
 for(int I=0;I<2;I++){
  //N-ter
  i=I;
  v[0]=v[1]=v[2]=0.0;
  x=(int)((ca[i][0]-map->orgxyz[0])*iv);
  y=(int)((ca[i][1]-map->orgxyz[1])*iv);
  z=(int)((ca[i][2]-map->orgxyz[2])*iv);
  ind=xydim*z+map->xdim*y+x;
  sco=(map->dens[ind]-map->map_t)*map->inv;
  	//6 vec
	for(int j=0;j<6;j++){
	 ind=xydim*(z+six_vec[j][2])+map->xdim*(y+six_vec[j][1])+(x+six_vec[j][0]);
	 sco2=(map->dens[ind]-map->map_t)*map->inv;
	 if(sco>sco2)
	  continue;
	 v[0]+=(sco2-sco)*six_vec[j][0];
	 v[1]+=(sco2-sco)*six_vec[j][1];
	 v[2]+=(sco2-sco)*six_vec[j][2];
	}
  d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if(d>0){
   vec[i][0]=r1*vec[i][0]+r2*v[0]/d;
   vec[i][1]=r1*vec[i][1]+r2*v[1]/d;
   vec[i][2]=r1*vec[i][2]+r2*v[2]/d;
  }
  sum+=d;
  //C-ter
  v[0]=v[1]=v[2]=0.0;
  i=n+I;
  x=(int)((ca[i][0]-map->orgxyz[0])*iv);
  y=(int)((ca[i][1]-map->orgxyz[1])*iv);
  z=(int)((ca[i][2]-map->orgxyz[2])*iv);
  ind=xydim*z+map->xdim*y+x;
  sco=(map->dens[ind]-map->map_t)*map->inv;
  	//6 vec
	for(int j=0;j<6;j++){
	 ind=xydim*(z+six_vec[j][2])+map->xdim*(y+six_vec[j][1])+(x+six_vec[j][0]);
	 sco2=(map->dens[ind]-map->map_t)*map->inv;
	 if(sco>sco2)
	  continue;
	 v[0]+=(sco2-sco)*six_vec[j][0];
	 v[1]+=(sco2-sco)*six_vec[j][1];
	 v[2]+=(sco2-sco)*six_vec[j][2];
	}
  d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if(d>0){
   vec[i][0]=r1*vec[i][0]+r2*v[0]/d;
   vec[i][1]=r1*vec[i][1]+r2*v[1]/d;
   vec[i][2]=r1*vec[i][2]+r2*v[2]/d;
  }
  sum+=d;
 }


 //Clear vec[2]
 //vec[2][0]=vec[2][1]=vec[2][2]=0.00;
 return sum;
}

double GetVecCa(double **in,double **vec,int n,MRC *map,double w){
 //double vec[MAX_FLEN][3];
 int p;
 double d2,d,v[3],f;
 double sum=0;
 double r1=1.0/(1.0+w);
 double r2=w/(1.0+w);

 //density based vec
 double sco,sco2;
 int i,ind,x,y,z;
 int xydim=map->xydim;

 double iv=map->iv_width;

/*
             dens - dt
MAP_SCORE=  ---------- = (dens - dt)*map->inv;
             dmax - dt
*/


 //Main chain positions
 for(int i=0;i<n;i++){
  v[0]=v[1]=v[2]=0.0;
  x=(int)((in[i][0]-map->orgxyz[0])*iv);
  y=(int)((in[i][1]-map->orgxyz[1])*iv);
  z=(int)((in[i][2]-map->orgxyz[2])*iv);
  ind=xydim*z+map->xdim*y+x;
  sco=(map->dens[ind]-map->map_t)*map->inv;
	//6 vec
	for(int j=0;j<6;j++){
	 ind=xydim*(z+six_vec[j][2])+map->xdim*(y+six_vec[j][1])+(x+six_vec[j][0]);
	 sco2=(map->dens[ind]-map->map_t)*map->inv;
	 if(sco>sco2)
	  continue;
	 v[0]+=(sco2-sco)*six_vec[j][0];
	 v[1]+=(sco2-sco)*six_vec[j][1];
	 v[2]+=(sco2-sco)*six_vec[j][2];
	}
  //unit vec

  d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if(d>0){
   vec[i][0]=r1*vec[i][0]+r2*v[0]/d;
   vec[i][1]=r1*vec[i][1]+r2*v[1]/d;
   vec[i][2]=r1*vec[i][2]+r2*v[2]/d;
  }
  sum+=d;
 }
 //puts("End vec");
 //Clear vec[2]
 //vec[2][0]=vec[2][1]=vec[2][2]=0.00;
 return sum;
}


void VecMove(double **in,double **out,double **vec,int n,double w, struct drand48_data *buf){
 double tmp;
 for(int i=0;i<n;i++){
  //Update coords
  drand48_r(buf,&tmp);
  //out[i][0]=in[i][0]+rnd()*vec[i][0]*w;
  tmp-=0.1;
  out[i][0]=in[i][0]+tmp*vec[i][0]*w;
  drand48_r(buf,&tmp);
  tmp-=0.1;
  out[i][1]=in[i][1]+tmp*vec[i][1]*w;
  drand48_r(buf,&tmp);
  tmp-=0.1;
  out[i][2]=in[i][2]+tmp*vec[i][2]*w;
  //out[i][0]=in[i][0]+0.0009765625*(rand()&1023)*vec[i][0]*w;
 }
}

//double swap_helix(MEMO *m,int *ss,int pos,int len,MODEL *hlib)
double swap_helix(MEMO *m,SEQFG *s,MODEL *hlib)
{
 double **in=m->frag;
 double **out=m->ca;
 int p,i,j;
 double cacoords[20][3],tmpstat[20][3],tmpcoords[20][3],rmsd;
 int len=s->l;
 int *ss=s->ss;
 //out ca
 //if L(helix)>4 swap with hlib

 //copy
 for(i=0;i<len;i++){
  out[i][0]=in[i][0];
  out[i][1]=in[i][1];
  out[i][2]=in[i][2];
 }

 if(hlib->NumOfCd==0)
  return 0;

 for(i=0;i<len;i++){
  //printf("%d ss= %d seq= %d\n",i,s->ss[i],s->seq[i]);

  if(ss[i]!=1) continue;
  for(j=i+1;i+j<len;j++){
   if(ss[i+j]!=1)
    break;
  }
  if(j-i>3){
   int n=j-i;
   //printf("#Swap pos %d: %d-%d\n",s->pos,i,j-1);
   //copy
   for(int k=0;k<n;k++){
    cacoords[k][0]=in[i+k][0];
    cacoords[k][1]=in[i+k][1];
    cacoords[k][2]=in[i+k][2];

    tmpstat[k][0]=hlib->xyz[k][0];
    tmpstat[k][1]=hlib->xyz[k][1];
    tmpstat[k][2]=hlib->xyz[k][2];

    tmpcoords[k][0]=hlib->xyz[k][0];
    tmpcoords[k][1]=hlib->xyz[k][1];
    tmpcoords[k][2]=hlib->xyz[k][2];
    //printf("%d 1 %f %f %f\n",k,cacoords[k][0],cacoords[k][1],cacoords[k][2]);
    //printf("%d 2 %f %f %f\n",k,tmpstat[k][0],tmpstat[k][1],tmpstat[k][2]);
    //printf("%d 3 %f %f %f\n",k,hlib->xyz[k][0],hlib->xyz[k][1],hlib->xyz[k][2]);
   }
   rmsd=superimpose2(cacoords, tmpstat, n, tmpcoords, n);
   //copy

   for(int k=0;k<n;k++){
    out[i+k][0]=tmpcoords[k][0];
    out[i+k][1]=tmpcoords[k][1];
    out[i+k][2]=tmpcoords[k][2];
   }

  }
  i=j;
 }
 //keep 3rd residues
 //out[2][0]=in[2][0];
 //out[2][1]=in[2][1];
 //out[2][2]=in[2][2];
 return 0;
}

double GetVecSide(double **mch,double ***sch,double **vec,int *seq,int n,MRC *map,double w){
 int p;
 double d2,d,v[3],vs[3],f;
 double sum=0;
 double r1=1.0/(1.0+w);
 double r2=w/(1.0+w);
 double rnsc;
 //density based vec
 double sco,sco2;
 int i,ind,x,y,z;
 int xydim=map->xydim;

 double iv=map->iv_width;
 double rvd;

/*
             dens - dt
MAP_SCORE=  ---------- = (dens - dt)*map->inv;
             dmax - dt
*/

 for(int i=0;i<n;i++){
  v[0]=v[1]=v[2]=0.0;
  vs[0]=vs[1]=vs[2]=0.0;
  

  int nsc = nheavy[seq[i]];


	//main-chain
	for(int atm=-1;atm<2;atm++){
	 int I=2*i+atm;
	
	 if(I<0||I>2*n-2)
	  continue;
   	 x=(int)((mch[I][0]-map->orgxyz[0])*iv);
   	 y=(int)((mch[I][1]-map->orgxyz[1])*iv);
   	 z=(int)((mch[I][2]-map->orgxyz[2])*iv);

   	 ind=xydim*z+map->xdim*y+x;
   	 //sco=(map->dens[ind]-map->map_t)*map->inv;

	 v[0]+=map->grd[ind][0];
	 v[1]+=map->grd[ind][1];
	 v[2]+=map->grd[ind][2];

	}


	//side-chain
	if(sch[i]!=NULL && nsc>0){
	rnsc=1.00/(double)nsc;
  	 for(int atm=0;atm<nsc;atm++){
   	  x=(int)((sch[i][atm][0]-map->orgxyz[0])*iv);
   	  y=(int)((sch[i][atm][1]-map->orgxyz[1])*iv);
   	  z=(int)((sch[i][atm][2]-map->orgxyz[2])*iv);

   	  ind=xydim*z+map->xdim*y+x;
   	  //sco=(map->dens[ind]-map->map_t)*map->inv;

	  vs[0]+=map->grd[ind][0];
	  vs[1]+=map->grd[ind][1];
	  vs[2]+=map->grd[ind][2];
	  //printf("%f %f %f\n",map->grd[ind][0],map->grd[ind][1],map->grd[ind][2]);
	 }

	 v[0]+=vs[0]*rnsc;
	 v[1]+=vs[1]*rnsc;
	 v[2]+=vs[2]*rnsc;

	}

  //unit vec
  d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  rvd=1.000/d;
  if(d>1.000){
   vec[i][0]=vec[i][0]+w*v[0]*rvd;
   vec[i][1]=vec[i][1]+w*v[1]*rvd;
   vec[i][2]=vec[i][2]+w*v[2]*rvd;
  } else if(d>0.00){
   vec[i][0]=vec[i][0]+w*v[0];
   vec[i][1]=vec[i][1]+w*v[1];
   vec[i][2]=vec[i][2]+w*v[2];
  }
  sum+=d;
 }

 //Clear vec[2]
 //vec[2][0]=vec[2][1]=vec[2][2]=0.00;
 return sum;
}

double GetVecPsback(double **in,double **vec,int n,MRC *map,double w){
 int p;
 double d2,d,v[3],f;
 double sum=0;
 double r1=1.0/(1.0+w);
 double r2=w/(1.0+w);

 //density based vec
 double sco,sco2;
 int i,ind,x,y,z;
 int xydim=map->xydim;

 double iv=map->iv_width;

/*
             dens - dt
MAP_SCORE=  ---------- = (dens - dt)*map->inv;
             dmax - dt
*/

 //pseudo-backbone positions
 for(int I=0;I<n;I++){
  v[0]=v[1]=v[2]=0.0;
  sco=0;
  for(int atm=-1;atm<2;atm++){
   i=I*2+atm;
	if(i<0||i>2*n-2)
	 continue;
  
   x=(int)((in[i][0]-map->orgxyz[0])*iv);
   y=(int)((in[i][1]-map->orgxyz[1])*iv);
   z=(int)((in[i][2]-map->orgxyz[2])*iv);
   ind=xydim*z+map->xdim*y+x;
   //sco=(map->dens[ind]-map->map_t)*map->inv;

   v[0]+=map->grd[ind][0];
   v[1]+=map->grd[ind][1];
   v[2]+=map->grd[ind][2];
   //printf("ind= %d %f %f %f\n",ind,v[0],v[1],v[2]);
  }
  //unit vec

  d=sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
  if(d>1.00){
   vec[I][0]=vec[I][0]+w*v[0]/d;
   vec[I][1]=vec[I][1]+w*v[1]/d;
   vec[I][2]=vec[I][2]+w*v[2]/d;
  }else if(d>0.00){
   vec[I][0]=vec[I][0]+w*v[0];
   vec[I][1]=vec[I][1]+w*v[1];
   vec[I][2]=vec[I][2]+w*v[2];
  }
  //printf("%f %f %f d= %f\n",v[0],v[1],v[2],d);
  sum+=d;
 }

 //Clear vec[2]
 //vec[2][0]=vec[2][1]=vec[2][2]=0.00;
 return sum;
}




//double opt_backbone_fromCA(MEMO *m,int *seq,int *ss,int pos,int len, MRC *map,bool SHOW)
void opt_backbone_fromCA(MEMO *m,SEQFG *g, MRC *map,MODEL *mod,bool SHOW)
{
 double ene=0;
 double **mch;//,vec[MAX_FLEN][3],cur[MAX_FLEN][3],tmp[MAX_FLEN][3];
 double **in=m->ca;
 double **out=m->cd;
 double **init=m->frag;
 int atm=0;
 int iter=0;
 struct drand48_data *buf=&(m->drand_buf);//for drand48

 int *ss=g->ss;
 int *seq=g->seq;
 int len=g->l;
 int pos=g->pos;

 //for(int i=0;i<len;i++){
 // printf("seq= %d ss= %d\n",g->seq[i],g->ss[i]);
 //}
 //copy CA to cur
 for(int i=0;i<len;i++){
  m->cur[i][0]=in[i][0];
  m->cur[i][1]=in[i][1];
  m->cur[i][2]=in[i][2];
 }

 double dvec=0;

 rebuild_psbackbone_fromCA(m->cur,m->mch,len);
 dvec=GetVec(m->cur,init,m->vec,len,ss);

 dvec+=GetVecPsback(m->mch,m->vec,len,map,cmd.Wvec);
 //dvec+=GetVecCa(m->cur,m->vec,len+4,map,cmd.Wvec);

 //Backbone opt working!!
 while(iter<cmd.NumIterMain){
  iter++;

  VecMove(m->cur,m->tmp,m->vec,len,cmd.MaxMove,buf);
  rebuild_psbackbone_fromCA(m->tmp,m->mch,len);

  //update cur
  for(int i=0;i<len;i++){
   m->cur[i][0]=m->tmp[i][0];
   m->cur[i][1]=m->tmp[i][1];
   m->cur[i][2]=m->tmp[i][2];
  }

  dvec=GetVec(m->cur,init,m->vec,len,ss);
  dvec+=GetVecPsback(m->mch,m->vec,len,map,cmd.Wvec);
  //dvec+=GetVecMain(m->mch,m->cur,m->vec,len,seq,pos,map,cmd.Wvec);
  //dvec+=GetVecCa(m->cur,m->vec,len+4,map,cmd.Wvec);

  //printf("dvec= %f %d\n",dvec,iter);

  if(dvec<0.1) //converged
   break;

  if(false && iter%10==0){
   //ShowFragCA(g, m->cur);
   rebuild_backbone_fromCA2(m->cur,m->mch,m->rbins,seq,len);
   ShowFragBack(g, m->mch);
  }


 }

 //ShowFragCA(g, m->cur);

 //Side-chain
 iter=0;
 double sco=0;
 double now_sco=0;

 now_sco=rebuild_sidechain_fromCA2(m->cur,m->sch,m->rbins,seq,len,map,m);
 rebuild_psbackbone_fromCA(m->cur,m->mch,len);

 dvec=GetVec(m->cur,init,m->vec,len,ss);//SS based
 dvec+=GetVecSide(m->mch,m->sch,m->vec,seq,len,map,cmd.Wvec); //main+side vs map

 //printf("now_sco= %f\n",now_sco);
 while(iter<cmd.NumIterSide){
  iter++;
  //ShowFragCA(g, m->cur);
  VecMove(m->cur,m->tmp,m->vec,len,cmd.MaxMove,buf);

  rebuild_psbackbone_fromCA(m->tmp,m->mch,len);
  sco=MapScore(m->mch,2*len-1,map);
  sco+=rebuild_sidechain_fromCA2(m->tmp,m->sch,m->rbins,seq,len,map,m);


  if(iter%10==0 && false)
   ShowFragSide(g,m->cur,m->sch);

  //minimize
  //if(sco<now_sco)
  // continue;

  now_sco=sco;
  //update
  for(int i=0;i<len;i++){
   m->cur[i][0]=m->tmp[i][0];
   m->cur[i][1]=m->tmp[i][1];
   m->cur[i][2]=m->tmp[i][2];
  }



  dvec=GetVec(m->cur,init,m->vec,len,ss);//SS based
  dvec+=GetVecSide(m->mch,m->sch,m->vec,seq,len,map,cmd.Wvec); //main+side vs map
  //printf("dvec= %f\n",dvec);
 }
 //printf("sco= %f\n",sco);

 //Shake-Zscore
 rebuild_backbone_fromCA2(m->cur,m->mch,m->rbins,seq,len);
 //get sco, shake,....
 ShakeZsco(mod,g,m->mch,m->sch,map,buf,cmd.Dshake*2.0, cmd.Nshake);
 //return sco;
}


//Just CA coords
void ShowFragCA(SEQFG *g,double **cd){
 int i;
 int len=g->l;
 int atm=1;
 int *seq=g->seq;
 printf("MODEL\n");
 for(i=0;i<len;i++){
  printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "CA ", RES_NAMES[seq[i]], ' ',1+i+g->pos,cd[i][0],cd[i][1], cd[i][2]);
         //atm++, "CA ", RES_NAMES[seq[i]], ' ',i+1,cd[i][0],cd[i][1], cd[i][2]);

 }
 printf("ENDMDL\n");
}

//Just main-chain coords
void ShowFragBack(SEQFG *g,double **mch){
 int i;
 int len=g->l;
 int atm=1;
 int *seq=g->seq;
 int pos=g->pos+1;
 int resnum;
 i=1;
 resnum=pos+1;
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "C  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+2][0],mch[4*i+2][1], mch[4*i+2][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "O  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+3][0],mch[4*i+3][1], mch[4*i+3][2]);


 for(i=2;i<len-2;i++){
  resnum=pos+i;
  	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "N  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i][0],mch[4*i][1], mch[4*i][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "CA ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+1][0],mch[4*i+1][1], mch[4*i+1][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "C  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+2][0],mch[4*i+2][1], mch[4*i+2][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "O  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+3][0],mch[4*i+3][1], mch[4*i+3][2]);
 }
 i=len-2;
 resnum=pos+i;
 	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "N  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i][0],mch[4*i][1], mch[4*i][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "CA ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+1][0],mch[4*i+1][1], mch[4*i+1][2]);
 printf("ENDMDL\n");
}


void ShowFragSide(SEQFG *g,double **ca,double ***sch){
 int i,j;
 int len=g->l;
 int anum=1;
 int *seq=g->seq;
 int nsc;
 int pos=g->pos+1;
 printf("#pos= %d\n",pos);
	for(i=0;i<len;i++){
 	 int resnum=pos+i;
	 nsc = nheavy[seq[i]];
	//show atoms

	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         anum++, "CA ", RES_NAMES[seq[i]], ' ',resnum,ca[i][0],ca[i][1],ca[i][2]);
	 if(sch[i]==NULL)
	  continue;
		for (j=0;j<nsc;j++) {
		 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
        	 anum++, heavy_atoms[seq[i]*10+j], RES_NAMES[seq[i]], ' ',resnum,sch[i][j][0], sch[i][j][1], sch[i][j][2]);
       		}
	}
 printf("ENDMDL\n");
}


void ShowFragBackSide(SEQFG *g,double **mch,double ***sch){
 int i,j;
 int len=g->l;
 int atm=1;
 int *seq=g->seq;
 int nsc;
 int pos=g->pos+1;
 printf("#pos= %d\n",pos);

 for(i=0;i<len;i++){
 	 int resnum=pos+i;
	 nsc = nheavy[seq[i]];
	//show atoms

	if(i>1 && i<len-1){
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "N  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i][0],mch[4*i][1], mch[4*i][2]);
	 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         atm++, "CA ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+1][0],mch[4*i+1][1], mch[4*i+1][2]);
	}

	if(i>0 && i<len-2){
	  printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
          atm++, "C  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+2][0],mch[4*i+2][1], mch[4*i+2][2]);
	  printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
          atm++, "O  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+3][0],mch[4*i+3][1], mch[4*i+3][2]);
	}


	 //printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         //anum++, "CA ", RES_NAMES[seq[i]], ' ',resnum,ca[i][0],ca[i][1],ca[i][2]);
	 if(sch[i]==NULL)
	  continue;
		for (j=0;j<nsc;j++) {
		 printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
        	 atm++, heavy_atoms[seq[i]*10+j], RES_NAMES[seq[i]], ' ',resnum,sch[i][j][0], sch[i][j][1], sch[i][j][2]);
       		}
	}
 printf("ENDMDL\n");


}


#define MAX_DATA_NUM 1000

double ShakeZsco(MODEL *mod,SEQFG *g,double **mch,double ***sch, MRC *m,struct drand48_data *buf, double r, int n){
 int i,j;
 int len=g->l;
 //int atm=1;
 int atm=0;
 int *seq=g->seq;
 int nsc;
 int pos=g->pos+1;
 //printf("#pos= %d\n",pos);
 double atm_cd[200][3];
 double tmp_cd[200][3];
 double data[MAX_DATA_NUM];

 if(n>MAX_DATA_NUM)
  n=MAX_DATA_NUM;

 for(i=0;i<len;i++){
 	 int resnum=pos+i;
	 nsc = nheavy[seq[i]];
	//show atoms

	if(i>1 && i<len-1){

	 atm_cd[atm][0]=mch[4*i][0]; atm_cd[atm][1]=mch[4*i][1]; atm_cd[atm][2]=mch[4*i][2];
	 atm++;
	 //printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         //atm++, "N  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i][0],mch[4*i][1], mch[4*i][2]);


	 atm_cd[atm][0]=mch[4*i+1][0]; atm_cd[atm][1]=mch[4*i+1][1]; atm_cd[atm][2]=mch[4*i+1][2];

	 //printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
         //atm++, "CA ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+1][0],mch[4*i+1][1], mch[4*i+1][2]);

	 atm++;

	}

	if(i>0 && i<len-2){
	  atm_cd[atm][0]=mch[4*i+2][0]; atm_cd[atm][1]=mch[4*i+2][1]; atm_cd[atm][2]=mch[4*i+2][2];
	  atm++;
	  //printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
          //atm++, "C  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+2][0],mch[4*i+2][1], mch[4*i+2][2]);
	  atm_cd[atm][0]=mch[4*i+3][0]; atm_cd[atm][1]=mch[4*i+3][1]; atm_cd[atm][2]=mch[4*i+3][2];
	  atm++;
	  //printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
          //atm++, "O  ", RES_NAMES[seq[i]], ' ',resnum,mch[4*i+3][0],mch[4*i+3][1], mch[4*i+3][2]);
	}

	 if(sch[i]==NULL)
	  continue;
		for (j=0;j<nsc;j++) {

	  	 atm_cd[atm][0]=sch[i][j][0]; atm_cd[atm][1]=sch[i][j][1]; atm_cd[atm][2]=sch[i][j][2];
	  	 atm++;
		 //printf("ATOM  %5d  %-3s %3s %c%4d    %8.3f%8.3f%8.3f\n",
        	 //atm++, heavy_atoms[seq[i]*10+j], RES_NAMES[seq[i]], ' ',resnum,sch[i][j][0], sch[i][j][1], sch[i][j][2]);
       		}
	}
 //printf("ENDMDL\n");
 //printf("Natm=%d\n",atm);
 //Shake
 double ori_sco=MapScore2(atm_cd,atm,m);

 double x,y,z,tmp_sco,sum;
 int cnt=0;

 //printf("xyz= 0 0 0 sco= %f\n",ori_sco);
 sum=0;
	while(cnt<n){
 	 drand48_r(buf,&x); drand48_r(buf,&y); drand48_r(buf,&z);
	 x=(x-0.5)*r;//-2.0~+2.0
	 y=(y-0.5)*r;//-2.0~+2.0
	 z=(z-0.5)*r;//-2.0~+2.0
	 for(int i=0;i<atm;i++){
	  tmp_cd[i][0]=atm_cd[i][0]+x;
	  tmp_cd[i][1]=atm_cd[i][1]+y;
	  tmp_cd[i][2]=atm_cd[i][2]+z;
	 }
 	 tmp_sco=MapScore2(tmp_cd,atm,m);
	 //printf("xyz= %.1f %.1f %.1f sco= %.3f\n",x,y,z,tmp_sco);
	 data[cnt]=tmp_sco;
	 sum+=tmp_sco;
	 cnt++;
	}

 sum/=(double)cnt;
 double std;
 for(int i=0;i<cnt;i++)
  std+=(data[i]-sum)*(data[i]-sum);
 std=sqrt(std/(double)cnt);
 mod->sco=ori_sco;
 mod->shake=(ori_sco-sum)/std;
 //printf("sco= %f shake=%f ave= %f std= %f\n",mod->sco,mod->shake,sum,std);
 return 0.00;
}

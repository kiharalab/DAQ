/*
caldep + fragment + sphere
*/
//#include <fftw3.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <omp.h>
#include "struct.h"
#include "func.h"
#include "mrc.h"
#include "thread.h"

//#include "scoring.h"

#define PDB_STRLEN 55

//Fragment

void malloc_error(char *a){
 fprintf(stderr,"malloc error in %s\n",a);
 exit(0);
}
double gettimeofday_sec()
{
    struct timeval tv;
    gettimeofday(&tv, NULL);
    return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

int readlist(char *fname,char **list){
 int num=0;
 FILE *fp;
 int len;
 if((fp=fopen(fname,"r"))==NULL)
  return FALSE;
 while(fgets(list[num],LIN,fp)!=NULL){
  len=strlen(list[num]);
  list[num][len-1]='\0';//ignore terminal \n
  num++;
 }
 fclose(fp);
 return TRUE;
}

int line_num(char *fname){
 int num=0;
 FILE *fp;
 char line[LIN];
 if((fp=fopen(fname,"r"))==NULL)
  return FALSE;
 while(fgets(line,LIN,fp)!=NULL){
  num++;
 }
 fclose(fp);
 return num;
}


CMD cmd;

int main(int argc, char **argv)
{
 double t1=gettimeofday_sec();
 double t4;
 POINTS pt;
 MRC mrc,pmrc;
 GRAPH g;
 TREE mst;
 PTBL ptbl;
 //Get Options
 if(chkcmdline(argc,argv,&cmd)==FALSE)
  return(0);

 //Set threads
 if(cmd.Nthr < omp_get_num_procs()){
  omp_set_num_threads(cmd.Nthr);
 }else{
  omp_set_num_threads(omp_get_num_procs());
 }

 //Density map
 if(readmrc(&mrc,cmd.filename))
  return(0);
 //prob table
 if(read_ptbl(&ptbl,cmd.pfilename))
  return(0);

 //Show Mode
 if(cmd.Mode==4){
  ShowPTBL(&mrc,&ptbl,cmd.Pcut);
  return 0;
 }



 SEQ seq;
/*
 if(readseq(&seq,cmd.sfilename)){
  puts("Error in seq file");
  return 0;
 }
 SEQFG sfg;
 //copy
 sfg.l = seq.len;
 sfg.seq=(int*)malloc(sizeof(int)*sfg.l);
 for(int i=0;i<sfg.l;i++){
  sfg.seq[i]=seq.seq_code[i];
  sfg.Pss[i][0]=seq.Pss[i][0];
  sfg.Pss[i][1]=seq.Pss[i][1];
  sfg.Pss[i][2]=seq.Pss[i][2];
  sfg.CID[i]=seq.CID[i];
 }
*/


 //Add Main-chain prob to MRC
 AddProbToMrc(&ptbl,&mrc,cmd.Wd,cmd.Pcut,cmd.map_t);
 puts("#FIN Adding PTBL to MRC");

 if(cmd.Mode==6){
  printf("##MQA mode: %s\n",cmd.mfilename);
  //read PDB
  PDB pdb;
  int Natm=CountAtom(cmd.mfilename);
	if(Natm==0){
	 printf("##ERROR Natm = 0 in %s\n",cmd.mfilename);
	 return 0;
	}
	MallocPdb(&pdb,Natm);
	readpdb(&pdb,cmd.mfilename,Natm);
	AssignProbToNodePDB(&ptbl,&pdb,&g,&mrc);
	//if(AssignSSPtoPDB(&seq,&pdb))
	// return 0;//Align Multi-chains
	ProbQA(&pdb,&g,cmd.Cmode);
  return 0;
 }

 return 0;

}


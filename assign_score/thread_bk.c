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
#include "thread.h"
#include "pulchra_func.h"

#define MAX_CD 10000
#define MAX_FLEN 200
#define MIN_UNIT 1.000

extern CMD cmd;

int load_models(MODEL **m,char **list, int Nlist, int *tot){
 int Nm=0;
 int Ncd=0;
 char line[LIN], buf[LIN];
 double cd[MAX_CD][3],b[MAX_CD],dis[MAX_CD],vec[3];//tmp data
 FILE *fpin;
 double x,y,z;
 int Ntot=0;
 double d,d2,dcut;
 int n_unit;

 dcut=MIN_UNIT*MIN_UNIT;

 for(int i=0;i<Nlist;i++){
  if(i%10==0)
   printf("#Load.. %d %s\n",i,list[i]);
	if((fpin=fopen(list[i],"r")) == NULL){ 
		fprintf(stderr,"Can't open %s\n",list[i]); 
		continue; 
	}
	while(fgets(line,LIN,fpin)){
		if(!strncmp(line,"ATOM",4)){
		 strncpy(buf,&line[30],8); buf[8]='\0'; 
		 x=atof(buf);
		 strncpy(buf,&line[38],8); buf[8]='\0'; 
		 y=atof(buf);
		 strncpy(buf,&line[46],8); buf[8]='\0'; 
		 z=atof(buf);
		 strncpy(buf,&line[60],6); buf[6]='\0'; 
		 b[Ncd]=atof(buf);
		 cd[Ncd][0]=x;
		 cd[Ncd][1]=y;
		 cd[Ncd][2]=z;
		 Ncd++;
		}
		if(!strncmp(line,"ENDMDL",6) 
		|| !strncmp(line,"TER",3)
		||!strncmp(line,"END",3)){
			if(Ncd==0)
			 continue;

		 //Add some missing coords
		 for(int i=1;i<Ncd;i++){
		  vec[0]=cd[i][0]-cd[i-1][0];
		  vec[1]=cd[i][1]-cd[i-1][1];
		  vec[2]=cd[i][2]-cd[i-1][2];
		  d2=vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
		  d=sqrt(d2);
		  dis[i]=d;
		  if(d<MIN_UNIT) 
		   continue;
		  //too long
		  n_unit=(int)(d/MIN_UNIT);
		  //printf("%d %d d2= %f %d\n",Nm,i,d,n_unit);
			//shift all data
			for(int j=Ncd-1;j>=i;j--){
			 cd[j+n_unit][0]=cd[j][0];
			 cd[j+n_unit][1]=cd[j][1];
			 cd[j+n_unit][2]=cd[j][2];

			 dis[j+n_unit]=dis[j];
			 b[j+n_unit]=b[j];
			}
			//add data
			vec[0]=vec[0]/(double)(n_unit+1);
			vec[1]=vec[1]/(double)(n_unit+1);
			vec[2]=vec[2]/(double)(n_unit+1);
			for(int j=0;j<n_unit;j++){
			 cd[i+j][0]=cd[i-1][0]+vec[0]*(double)(j+1);
			 cd[i+j][1]=cd[i-1][1]+vec[1]*(double)(j+1);
			 cd[i+j][2]=cd[i-1][2]+vec[2]*(double)(j+1);

			 dis[i+j]=d/(double)(n_unit+1);

			}
			Ncd+=n_unit;
		 }
		
		





		 //malloc
		 if((m[Nm]=(MODEL*)malloc(sizeof(MODEL)))==NULL)
		  return 0;
		 if((m[Nm]->xyz=(float **)malloc(sizeof(float *)*Ncd))==NULL)
		  return 0;
		 if((m[Nm]->b=(double *)calloc(Ncd,sizeof(double)))==NULL)
		  return 0;
		 if((m[Nm]->dis=(double *)malloc(sizeof(double)*Ncd))==NULL)
		  return 0;
		 	for(int i=0;i<Ncd;i++)
		 	 if((m[Nm]->xyz[i]=(float *)malloc(sizeof(float)*3))==NULL)
		  	  return 0;
		 //input
		 dis[0]=0;
		 for(int i=0;i<Ncd;i++){
		  m[Nm]->xyz[i][0]=cd[i][0];
		  m[Nm]->xyz[i][1]=cd[i][1];
		  m[Nm]->xyz[i][2]=cd[i][2];
		  m[Nm]->dis[i]=dis[i];
		  //ignore b-factor
		  //m[Nm]->b[i]=b[i];
		 }
		 m[Nm]->NumOfCd=Ncd;
		 Ntot+=Ncd;
		 //Ncd=0;
		 Nm++;
		 
		 //input reversed order
		 //malloc
		 if((m[Nm]=(MODEL*)malloc(sizeof(MODEL)))==NULL)
		  return 0;
		 if((m[Nm]->xyz=(float **)malloc(sizeof(float *)*Ncd))==NULL)
		  return 0;
		 if((m[Nm]->b=(double *)calloc(Ncd,sizeof(double)))==NULL)
		  return 0;
		 if((m[Nm]->dis=(double *)malloc(sizeof(double)*Ncd))==NULL)
		  return 0;
		 	for(int i=0;i<Ncd;i++)
		 	 if((m[Nm]->xyz[i]=(float *)malloc(sizeof(float)*3))==NULL)
		  	  return 0;

		 for(int i=0;i<Ncd;i++){
		  m[Nm]->xyz[i][0]=cd[Ncd-i-1][0];
		  m[Nm]->xyz[i][1]=cd[Ncd-i-1][1];
		  m[Nm]->xyz[i][2]=cd[Ncd-i-1][2];
		  m[Nm]->dis[i]=dis[Ncd-i];
		 }
		 dis[0]=0;
		 m[Nm]->NumOfCd=Ncd;
		 Ntot+=Ncd;
		 Ncd=0;
		 Nm++;
		}
	}
  fclose(fpin);
 }
 //printf("#Total Cd= %d\n",Ntot);
 *tot=Ntot;
 return Nm;
}

int load_single_model(MODEL *m,char *file){
 int Nm=0;
 int Ncd=0;
 char line[LIN], buf[LIN];
 double cd[MAX_CD][3],b[MAX_CD],dis[MAX_CD],vec[3];//tmp data
 FILE *fpin;
 double x,y,z;
 int Ntot=0;
 double d,d2,dcut;
 int n_unit;

   printf("#Load.. %s\n",file);
	if((fpin=fopen(file,"r")) == NULL){ 
		fprintf(stderr,"Can't open %s\n",file);
		m->NumOfCd=0;
		return 0; 
	}
	while(fgets(line,LIN,fpin)){
		if(!strncmp(line,"ATOM",4)){
		  if(strncmp(&line[13],"CA ",3))
		   continue;
		  if(line[16]!='A' && line[16]!=' ')
		   continue;
		 strncpy(buf,&line[30],8); buf[8]='\0'; 
		 x=atof(buf);
		 strncpy(buf,&line[38],8); buf[8]='\0'; 
		 y=atof(buf);
		 strncpy(buf,&line[46],8); buf[8]='\0'; 
		 z=atof(buf);
		 strncpy(buf,&line[60],6); buf[6]='\0'; 
		 b[Ncd]=atof(buf);
		 cd[Ncd][0]=x;
		 cd[Ncd][1]=y;
		 cd[Ncd][2]=z;
		 Ncd++;
		}
	}

		 //malloc
		 if((m->xyz=(float **)malloc(sizeof(float *)*Ncd))==NULL)
		  return 0;
		 if((m->b=(double *)calloc(Ncd,sizeof(double)))==NULL)
		  return 0;
		 if((m->dis=(double *)malloc(sizeof(double)*Ncd))==NULL)
		  return 0;
		 for(int i=0;i<Ncd;i++)
		  if((m->xyz[i]=(float *)malloc(sizeof(float)*3))==NULL)
		    return 0;
		 //input
		 dis[0]=0;
		 for(int i=0;i<Ncd;i++){
		  m->xyz[i][0]=cd[i][0];
		  m->xyz[i][1]=cd[i][1];
		  m->xyz[i][2]=cd[i][2];
		  m->dis[i]=dis[i];
		  //ignore b-factor
		  //m[Nm]->b[i]=b[i];
		 }
		 m->NumOfCd=Ncd;
		 
  fclose(fpin);
 //printf("#Total Cd= %d\n",Ntot);
 return Nm;
}





bool readseq(SEQ *seq,char *fname){
 int num=0;
 FILE *fp;
 char line[LIN],buf[LIN];
 int  ftype=-1;
 int i,len,Nch=0;
 char aa,ss;

 seq->len=0;
 //Detect File type
 if((fp=fopen(fname,"r"))==NULL)
  return true;
 while(fgets(line,LIN,fp)!=NULL){
  if(!strncmp(line,"#\tAA",4)){
   ftype=1;//*.spd3
   break;
  }
  if(!sscanf(line,">%s",seq->id)){
  //if(!strncmp(line,">",1)){
   ftype=0;//fasta format
   break;
  }
 }
 fclose(fp);

 //Read chain-break '/' or multiple chains

 Nch=0;
 if(ftype==0){
  puts("#Fasta file");
  if((fp=fopen(fname,"r"))==NULL)
   return true;
  while(fgets(line,LIN,fp)!=NULL){
   len=strlen(line);
   if(!strncmp(line,">",1)){
    //strncpy(seq->id,&line[1],len-1);
    //seq->id[len-2]='\0';
    printf("#ID= '%s'\n",seq->id);
    Nch++;
    if(Nch!=1){
     //add chain-break
     seq->seq_txt[seq->len]='/';
     seq->len++;
    }
   }else{
    strncpy(&seq->seq_txt[seq->len],line,len);
    seq->len+=len;
    if(seq->seq_txt[seq->len-1]=='\n')
     seq->len--;
   }
  }
  seq->seq_txt[seq->len]='\0';
  printf("#SEQ= '%s'\n",seq->seq_txt);
  fclose(fp);
  //convert to code
  for(i=0;i<seq->len;i++){
	if(seq->seq_txt[i]=='/'){
  	 seq->seq_code[i]=-9;
  	 seq->ss[i]=-1;//All C
	}else{
  	 seq->seq_code[i]=A2int(seq->seq_txt[i]);
  	 seq->ss[i]=3;//All C
	}
   //printf("%d %c %d\n",i,seq->seq_txt[i],seq->seq_code[i]);
  }
 }
 if(ftype==1){
  puts("#spd3 file");
  if((fp=fopen(fname,"r"))==NULL)
   return true;
  while(fgets(line,LIN,fp)!=NULL){
   if(!strncmp(line,"#",1)){
     strncpy(seq->id,"Unknown",7);
     seq->id[8]='\0';
     printf("#ID= '%s'\n",seq->id);
     //chain-break
     Nch++;
     if(Nch!=1){
      seq->seq_code[seq->len]=-9;
      seq->seq_txt[seq->len]='/';
      seq->ss[seq->len]=-1;
      seq->len++;
     }
   }else{
    sscanf(line,"%d\t%c\t%c\t%s",&i,&aa,&ss,buf);
    seq->seq_code[seq->len]=A2int(aa);
    seq->seq_txt[seq->len]=aa;

	//ss
	switch(ss){
                case 'C':seq->ss[seq->len]=3;
		 break;
                case 'H':seq->ss[seq->len]=1;
		 break;
                case 'E':seq->ss[seq->len]= 2;
		 break;
		default: seq->ss[seq->len]= -1;
        }

    //printf("%d %c %c %d %d\n",seq->len,seq->seq_txt[seq->len],ss,seq->seq_code[seq->len],seq->ss[seq->len]);
    seq->len++;
    
   }
  }
  fclose(fp);
  seq->seq_txt[seq->len]='\0';
  printf("#SEQ= '%s'\n",seq->seq_txt);
  for(int i=seq->len;i<seq->len+5;i++)
   seq->ss[seq->len]= -1;
 }
 if(ftype==-1){
  printf("Cannot detect file format..%s\n",fname);
  return true;
 }
 return false;
}

int splitseq(SEQ *seq,SEQFG *fg,int len){
 int n=0;
 int i,j,pos;
 int Nc=0;
 bool flag;
 int flen=len+4;
 for(i=0;i<seq->len-flen+1;i++){
  if(seq->seq_code[i]==-9){
   Nc++;
   continue;
  }
  //check
	flag=false;
	for(j=0;j<flen;j++){
	 if(seq->seq_code[i+j]==-9){//chain break
	  flag=true;
	  break;
	 }
	}
	if(flag==true) continue;
  //input
  fg[n].ss=(int *)malloc(sizeof(int)*flen);
  fg[n].seq=(int *)malloc(sizeof(int)*flen);
  for(j=0;j<flen;j++){
   fg[n].seq[j]=seq->seq_code[i+j];
   fg[n].ss[j]=seq->ss[i+j];
   printf("%3d",fg[n].seq[j]);
   printf("%3d",fg[n].ss[j]);
  }
  printf(" chain=%d\n",Nc);
  fg[n].cid=Nc;
  fg[n].l=flen;
  fg[n].pos=i;
  n++;
 }
 return n;
}


bool AssignFeatures(MODEL **m, int Nm,FRAG **f, int Nf, MRC *map){
 //int i,j;
 int xydim=map->xdim*map->ydim;
 double t=cmd.map_t;
 double d2cut=cmd.dcut*cmd.dcut;
 
 //each model
 #pragma omp parallel for schedule(dynamic,5)
 for(int i=0;i<Nm;i++){
  int ind,bpos;
  double d2,bdis2;
  double pos[3];
  int Natm=1;
	//Assign density data
	for(int x=0;x<map->xdim;x++){
	 pos[0]=x*map->widthx+map->orgxyz[0];
	for(int y=0;y<map->ydim;y++){
	 pos[1]=y*map->widthy+map->orgxyz[1];
	for(int z=0;z<map->zdim;z++){
	 pos[2]=z*map->widthz+map->orgxyz[2];
	 ind=z*xydim+y*map->xdim+x;

	 if(map->dens[ind]<t)
	  continue;
	 bpos=-1;
	 bdis2=d2cut;
		//find closest cd
		for(int a=0;a<m[i]->NumOfCd;a++){
		 d2=(m[i]->xyz[a][0]-pos[0])*(m[i]->xyz[a][0]-pos[0])
		   +(m[i]->xyz[a][1]-pos[1])*(m[i]->xyz[a][1]-pos[1])
		   +(m[i]->xyz[a][2]-pos[2])*(m[i]->xyz[a][2]-pos[2]);
		 if(d2<bdis2){
		  bpos=a;
		  bdis2=d2;
		 }
		}
		if(bpos==-1)
	  	 continue;
		m[i]->b[bpos]+=map->dens[ind];
		//printf("%d %d %d -> %d:%d %f\n",x,y,z,i,bpos,m[i]->b[bpos]);
		printf("ATOM  %5d  CA  ALA%6d    ",Natm,bpos+1);
   printf("%8.3f%8.3f%8.3f%6.2f%6.2f GRID\n",pos[0],pos[1],pos[2],1.0,map->dens[ind]);
		Natm++;
	}}}
 }
 return false;
}

void ShowMainPath(MODEL **m,int Nm){

 int i,j,k;
 int Natm=1;
 double tmp[3];
 for(int mod=0;mod<Nm;mod++){
  printf("MODEL %d\n",mod+1);
  Natm=1;
  for(i=0;i<m[mod]->NumOfCd;i++){
   printf("ATOM  %5d  CA  ALA%6d    ",Natm,Natm);
   printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",m[mod]->xyz[i][0],m[mod]->xyz[i][1],m[mod]->xyz[i][2],1.0,m[mod]->b[i]);
   Natm++;
  }
  printf("ENDMDL\n");
 }
}


bool NR_frag(FRAG **f,int *Nf,MODEL **m,int Nm,double flen){
 printf("#fragment_length= %f\n",flen);
 int n=0;
 //Assign st,end,length
 for(int mod=0;mod<Nm;mod++){
  for(int pos=0;pos<m[mod]->NumOfCd;pos++){
   double len=0;
   int ad;
   for(ad=1;pos+ad<m[mod]->NumOfCd;ad++){
    len+=m[mod]->dis[pos+ad];
    if(len>flen)
     break;
   }
   //Terminal
   if(len<=flen/1.5)
    break;

   if((f[n]=(FRAG *)malloc(sizeof(FRAG)))==NULL)
    return true;
   //data
   f[n]->st=pos;
   f[n]->ed=pos+ad;
   if(pos+ad>=m[mod]->NumOfCd)
    f[n]->ed=m[mod]->NumOfCd-1;

   f[n]->mod=mod;
   f[n]->dist=len;

   //printf("frag %d %d:%d %f\n",mod,pos,pos+ad,len);

   n++;
  }
 }
 printf("#Nf= %d\n",n);
 int NR_cnt;
 #pragma omp parallel for schedule(dynamic,5) reduction(+:NR_cnt)
 for(int f1=0;f1<n;f1++){
  bool flag=false;
 	for(int f2=f1+1;f2<n;f2++){
	 if(f[f1]->dist!=f[f2]->dist)
	  continue;
	 if(f[f1]->mod==f[f2]->mod)
	  continue;

	 if(f[f1]->ed-f[f1]->st!=f[f2]->ed-f[f2]->st)
	  continue;

	 int m1=f[f1]->mod;
	 int m2=f[f2]->mod;
	 float c=0;
		for(int p=0;p+f[f1]->st<=f[f1]->ed;p++){
		 int p1=p+f[f1]->st;
		 int p2=p+f[f2]->st;
		 c+=(m[m1]->xyz[p1][0]-m[m2]->xyz[p2][0])*(m[m1]->xyz[p1][0]-m[m2]->xyz[p2][0])
		   +(m[m1]->xyz[p1][1]-m[m2]->xyz[p2][1])*(m[m1]->xyz[p1][1]-m[m2]->xyz[p2][1])
		   +(m[m1]->xyz[p1][2]-m[m2]->xyz[p2][2])*(m[m1]->xyz[p1][2]-m[m2]->xyz[p2][2]);
		}
	 if(c!=0.00){
	  //different
	  continue;
	 }
	  //printf("%.3f %.3f %.3f %.3f\n",m[m1]->xyz[f[f1]->st][0],m[m2]->xyz[f[f2]->st][0],f[f1]->dist,f[f2]->dist);
	 //same
	 flag=true;
	 break;
	}
  if(flag==false){
   f[f1]->nr=true;
   NR_cnt++;
   //printf("#ID %d %d:%d:%d\n",f1,f[f1]->mod,f[f1]->st,f[f1]->ed);
  }else{
   f[f1]->nr=false;
  }
 }
 printf("#NRfrag= %d\n",NR_cnt);
 //Shift Data
 NR_cnt=0;
 for(int f1=0;f1<n;f1++){
  if(f[f1]->nr==false)
   continue;
  f[NR_cnt]=f[f1];
  NR_cnt++;
 }
 *Nf=NR_cnt;
}

double *rnd_tbl;

bool Malloc_Memo(MEMO *m,int flen){

 m->frag=(double**)malloc(sizeof(double *)*(flen+4));
	for(int j=0;j<flen+4;j++) 
	 m->frag[j]=(double*)malloc(sizeof(double)*3);

	m->ca=(double**)malloc(sizeof(double *)*(flen+4));
	for(int j=0;j<flen+4;j++) 
	 m->ca[j]=(double*)malloc(sizeof(double)*3);

	m->cur=(double**)malloc(sizeof(double *)*(flen+4));
	for(int j=0;j<flen+4;j++) 
	 m->cur[j]=(double*)malloc(sizeof(double)*3);

	m->tmp=(double**)malloc(sizeof(double *)*(flen+4));
	for(int j=0;j<flen+4;j++) 
	 m->tmp[j]=(double*)malloc(sizeof(double)*3);

	m->vec=(double**)malloc(sizeof(double *)*(flen+4));
	for(int j=0;j<flen+4;j++) 
	 m->vec[j]=(double*)malloc(sizeof(double)*3);

	m->mch=(double**)malloc(sizeof(double *)*(flen)*10);
	for(int j=0;j<flen*10;j++) 
	 m->mch[j]=(double*)malloc(sizeof(double)*3);

	m->sch=(double***)malloc(sizeof(double **)*(flen+4));	

	m->cd=(double**)malloc(sizeof(double *)*(MAX_FLEN));
	for(int j=0;j<MAX_FLEN;j++) 
	 m->cd[j]=(double*)malloc(sizeof(double)*3);

	m->rbins=(int**)malloc(sizeof(int *)*(MAX_FLEN));
	for(int j=0;j<MAX_FLEN;j++) 
   	 m->rbins[j]=(int*)malloc(sizeof(int)*3);
	
	m->stmp=(double**)malloc(sizeof(double *)*(flen)*10);
	for(int j=0;j<flen*10;j++) 
	 m->stmp[j]=(double*)malloc(sizeof(double)*3);

	m->slib=(double****)malloc(sizeof(double ***)*(flen+4));//sidechain lib
	for(int j=0;j<flen+4;j++){
	 m->slib[j]=(double***)malloc(sizeof(double**)*10);//top10 rot
	 for(int k=0;k<10;k++){
	  m->slib[j][k]=(double**)malloc(sizeof(double*)*20);//#atoms
	  for(int l=0;l<20;l++)
	   m->slib[j][k][l]=(double*)malloc(sizeof(double)*3);//cd
	 }
	}

 return false;
}



//bool ThreadFrag(MODEL **m,SEQ *s,FRAG **f, int Nf, MRC *map,MODEL *hlib, MODEL **fg){
bool ThreadFrag(MODEL **m,SEQFG *sfg,int Nsfg,FRAG **f, int Nf, MRC *map,MODEL *hlib, MODEL **fg){
 SEQ *s;
 int flen=cmd.frag_len;
 unsigned int Ne=0;
 int i;
 //MODEL *fg;
 int dim=Nsfg;
 double gscale=cmd.gscale;
 MEMO *memo;//memory space

 printf("#Searching fragments*seq: %d fgs * %d seq = %d\n",Nf,Nsfg,Nf*Nsfg);

 puts("#Start frag malloc");
 for(int i=0;i<Nf*(Nsfg);i++){
  if((fg[i]=(MODEL *)malloc(sizeof(MODEL)))==NULL)
   return true;
	//printf("%d\n",i);
  if((fg[i]->xyz=(float **)malloc(sizeof(float *)*(flen+4)))==NULL)
   return true;
  for(int j=0;j<flen+4;j++)
   if((fg[i]->xyz[j]=(float *)malloc(sizeof(float)*4))==NULL)
    return true;
 }
 puts("#End frag malloc");
 puts("#Start MemoSpace malloc");
 int Nth=omp_get_max_threads();
 if((memo=(MEMO *)malloc(sizeof(MEMO)*Nth))==NULL)
  return true;
 for(int i=0;i<Nth;i++){

	srand48_r(i, &memo[i].drand_buf);//Set seed and buf

/*
	memo[i].frag=(double**)malloc(sizeof(double *)*(flen+4));
	for(int j=0;j<flen+4;j++) 
	 memo[i].frag[j]=(double*)malloc(sizeof(double)*3);

	memo[i].ca=(double**)malloc(sizeof(double *)*(flen+4));
	for(int j=0;j<flen+4;j++) 
	 memo[i].ca[j]=(double*)malloc(sizeof(double)*3);

	memo[i].cur=(double**)malloc(sizeof(double *)*(flen+4));
	for(int j=0;j<flen+4;j++) 
	 memo[i].cur[j]=(double*)malloc(sizeof(double)*3);

	memo[i].tmp=(double**)malloc(sizeof(double *)*(flen+4));
	for(int j=0;j<flen+4;j++) 
	 memo[i].tmp[j]=(double*)malloc(sizeof(double)*3);

	memo[i].vec=(double**)malloc(sizeof(double *)*(flen+4));
	for(int j=0;j<flen+4;j++) 
	 memo[i].vec[j]=(double*)malloc(sizeof(double)*3);

	memo[i].mch=(double**)malloc(sizeof(double *)*(flen)*10);
	for(int j=0;j<flen*10;j++) 
	 memo[i].mch[j]=(double*)malloc(sizeof(double)*3);

	memo[i].sch=(double***)malloc(sizeof(double **)*(flen+4));	

	memo[i].cd=(double**)malloc(sizeof(double *)*(MAX_FLEN));
	for(int j=0;j<MAX_FLEN;j++) 
	 memo[i].cd[j]=(double*)malloc(sizeof(double)*3);

	memo[i].rbins=(int**)malloc(sizeof(int *)*(MAX_FLEN));
	for(int j=0;j<MAX_FLEN;j++) 
   	 memo[i].rbins[j]=(int*)malloc(sizeof(int)*3);
	
	memo[i].stmp=(double**)malloc(sizeof(double *)*(flen)*10);
	for(int j=0;j<flen*10;j++) 
	 memo[i].stmp[j]=(double*)malloc(sizeof(double)*3);

	memo[i].slib=(double****)malloc(sizeof(double ***)*(flen+4));//sidechain lib
	for(int j=0;j<flen+4;j++){
	 memo[i].slib[j]=(double***)malloc(sizeof(double**)*10);//top10 rot
	 for(int k=0;k<10;k++){
	  memo[i].slib[j][k]=(double**)malloc(sizeof(double*)*20);//#atoms
	  for(int l=0;l<20;l++)
	   memo[i].slib[j][k][l]=(double*)malloc(sizeof(double)*3);//cd
	 }
	}
  */
  Malloc_Memo(&memo[i],flen);
 }
 puts("#End frag malloc");

 //return false;



 #pragma omp parallel for reduction(+:Ne) schedule(dynamic,10)
 for(int i=0;i<Nf;i++){
  //double cd[MAX_FLEN][3];//Max 500 positions
  double msco[5],ssco,casco[5];
  int Ncd=f[i]->ed-f[i]->st+1;
  int seq[100],bestmod,idx;
  double bestsco,sco;
  int th=omp_get_thread_num();
  int Nmodel=0;
  MEMO *mm=&(memo[th]);//for multi threads
  printf("#%d (%.2f) m: %d f: %d-%d len=%.3f %d\n",i,(float)i/Nf,f[i]->mod,f[i]->st,f[i]->ed,f[i]->dist,Ncd);

  //if(f[i]->mod==0)
  // continue;

  //copy coordinates

  for(int p=0;p<Ncd;p++){
   mm->cd[p][0]=m[f[i]->mod]->xyz[f[i]->st+p][0];
   mm->cd[p][1]=m[f[i]->mod]->xyz[f[i]->st+p][1];
   mm->cd[p][2]=m[f[i]->mod]->xyz[f[i]->st+p][2];
  }

  //puts("#Start");

  //cd -> frag
  //Nmodel=BuildCa(mm->frag,mm->cd,Ncd,flen+4,map,mm);//-2..0..flen, flen+1
  //CA model cd -> frag
  Nmodel=BuildCa(mm,Ncd,flen+4,map);//-2..0..flen, flen+1
   //printf("#model %d\n",Nmodel);
  //init
  //No model
  if(Nmodel==0){
   for(int j=0;j<Nsfg;j++){
    idx=i*dim+j;
    fg[idx]->sco=0;
   }
   continue;
  }
  //Build main-chain fragments with secodary structure info
  for(int j=0;j<Nsfg;j++){
   idx=i*dim+j;

   //if(sfg[j].pos >3  )
   // continue;

   //Swap helix frag->ca
   swap_helix(mm,&sfg[j],hlib);	
   //ca->cur, mch
   //sco=opt_backbone_fromCA(mm,s->seq_code,s->ss,j,flen,map,false);
   opt_backbone_fromCA(mm,&sfg[j],map,fg[idx],false);
   //printf("##sco=%f\n",sco);
   Ne++;

	//input data
	//idx=i*dim+j;
	//fg[idx]->sco=sco;

	 for(int k=0;k<flen+4;k++){
	  fg[idx]->xyz[k][0]=mm->cur[k][0];
	  fg[idx]->xyz[k][1]=mm->cur[k][1];
	  fg[idx]->xyz[k][2]=mm->cur[k][2];
	 }

  }


 }

 printf("#Total generated Fragments= %d\n",Ne);
 return false;
}

double dev_ca(float **a,int A,float **b,int B,int n){
 double d=0;

 for(int i=0;i<n;i++){
  d+=(a[A+i][0]-b[B+i][0])*(a[A+i][0]-b[B+i][0])+
     (a[A+i][1]-b[B+i][1])*(a[A+i][1]-b[B+i][1])+
     (a[A+i][2]-b[B+i][2])*(a[A+i][2]-b[B+i][2]);
 }
 d=sqrt(d/(double)n);
 return d;
}

bool Evaluation(MODEL **m,SEQFG *s,int Nsfg,int Nf,MODEL *ref, MRC *map){
 int flen=cmd.frag_len;
 unsigned int Ne=0;
 int i;
 //MODEL *fg;
 int dim=Nsfg;
 double dev=0;
 double ssco;
 double **tmp,**mch;
 int **rbins,*seq;

 MEMO memo;

 Malloc_Memo(&memo,flen);

 tmp=(double**)malloc(sizeof(double *)*(flen+4));
 for(int j=0;j<flen+4;j++) tmp[j]=(double*)malloc(sizeof(double)*3);

 mch=(double**)malloc(sizeof(double *)*((flen+4)*4));
 for(int j=0;j<(flen+4)*4;j++) mch[j]=(double*)malloc(sizeof(double)*3);

 rbins=(int**)malloc(sizeof(int *)*(MAX_FLEN));
 for(int j=0;j<MAX_FLEN;j++) rbins[j]=(int*)malloc(sizeof(int)*3);


 //clustering remove very close structure dev<1.0
 #pragma omp parallel for schedule(dynamic,5)
 for(int pos=0;pos<Nsfg;pos++){
  double dev=0;
	for(int i=0;i<Nf-1;i++){
	 if(m[i*dim+pos]->sco < 0.00)
	  continue;
	 for(int j=i+1;j<Nf;j++){
	  if(m[j*dim+pos]->sco < 0.00)
	   continue;
	  dev=dev_ca(m[i*dim+pos]->xyz,0,m[j*dim+pos]->xyz,0,flen+4);
	  if(dev<1.0){
	   //printf("Same pos=%d %d %d %f\n",pos,i,j,dev);
	   //keep better
	   //if(m[i*dim+pos]->sco > m[j*dim+pos]->sco){
	   if(m[i*dim+pos]->shake > m[j*dim+pos]->shake){//Keep better shake score
	    m[j*dim+pos]->sco=-1.00;
	    m[j*dim+pos]->shake=-1.00;
	   }else{
	    m[i*dim+pos]->sco=-1.00;
	    m[i*dim+pos]->shake=-1.00;
	    break;
	   }
	  }
	 }
	}
 }


 //Sequence Position
 int Ntotal=0;
 #pragma omp parallel for reduction(+:Ntotal) schedule(dynamic,5)
 for(int pos=0;pos<Nsfg;pos++){
  double sum=0;
  double std=0;
  int Nact=0;
	//Zscore of raw-score
	for(int i=0;i<Nf;i++){
	 if(m[i*dim+pos]->sco>0){
	  sum+=m[i*dim+pos]->sco;
	  Nact++;
	 }
	}
	sum/=(double)Nact;
	for(int i=0;i<Nf;i++)
	 if(m[i*dim+pos]->sco>0)
	  std+=(m[i*dim+pos]->sco-sum)*(m[i*dim+pos]->sco-sum);
	std=sqrt(std/(double)Nact);
	for(int i=0;i<Nf;i++)
	 m[i*dim+pos]->z=(m[i*dim+pos]->sco-sum)/std;

	//Zscore of shake score
	sum=0;
	std=0;
	for(int i=0;i<Nf;i++){
	 if(m[i*dim+pos]->sco>0){
	  sum+=m[i*dim+pos]->shake;
	 }
	}
	sum/=(double)Nact;
	for(int i=0;i<Nf;i++)
	 if(m[i*dim+pos]->sco>0)
	  std+=(m[i*dim+pos]->shake-sum)*(m[i*dim+pos]->shake-sum);
	std=sqrt(std/(double)Nact);
	for(int i=0;i<Nf;i++)
	 m[i*dim+pos]->zshake=(m[i*dim+pos]->shake-sum)/std;

	//printf("Pos %d ave= %f std= %f\n",pos,sum,std);
  Ntotal+=Nact;
 }

 printf("#Clustering(<1.0A)  %d\n",Ntotal);

 //return 0;
 for(int pos=0;pos<Nsfg;pos++){
  for(int i=0;i<Nf;i++){

   //ignore
   if(m[i*dim+pos]->sco<=0)
    continue;

   if(cmd.refmode==true)
    dev=dev_ca(ref->xyz,s[pos].pos,m[i*dim+pos]->xyz,0,flen+4);

   //out
   if(true && m[i*dim+pos]->z > cmd.zcut){
    printf("FRG %d Dev= %.3f Z= %.3f Zshake= %.3f ",s[pos].pos,dev,m[i*dim+pos]->z,m[i*dim+pos]->zshake);


    printf("CD= ");
	for(int j=0;j<flen+4;j++){
	 tmp[j][0]=(double)m[i*dim+pos]->xyz[j][0];
	 tmp[j][1]=(double)m[i*dim+pos]->xyz[j][1];
	 tmp[j][2]=(double)m[i*dim+pos]->xyz[j][2];
	 printf("%.3f,%.3f,%.3f",tmp[j][0],tmp[j][1],tmp[j][2]);
	 if(j!=flen+3)
	  printf(",");
	}
	printf("\n");
    //ShowFragCA(&s[pos],tmp);
   }
   if(m[i*dim+pos]->zshake>3.0 && true){
    //printf("#Pos= %d Dev= %f Z= %f\n",s[pos].pos,dev,m[i*dim+pos]->z);
	//copy
	int atm=1;
	for(int j=0;j<flen+4;j++){
	 tmp[j][0]=(double)m[i*dim+pos]->xyz[j][0];
	 tmp[j][1]=(double)m[i*dim+pos]->xyz[j][1];
	 tmp[j][2]=(double)m[i*dim+pos]->xyz[j][2];
	}
    	rebuild_backbone_fromCA2(tmp,mch,rbins,s[pos].seq,flen+4);
	rebuild_sidechain_fromCA2(tmp,memo.sch,memo.rbins,s[pos].seq,flen+4,map,&memo);
	ShowFragBackSide(&s[pos], mch,memo.sch);
   }
  }
 }
}

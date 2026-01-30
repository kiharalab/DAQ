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
#include "dp.h"

bool readseq(SEQ *seq,char *fname){
 int num=0;
 FILE *fp;
 char line[LIN],buf[LIN];
 int  ftype=-1;
 int i,len,Nch=0;
 char aa,ss;
 float fbuf;

 seq->len=0;
 //Detect File type
 if((fp=fopen(fname,"r"))==NULL)
  return true;
 while(fgets(line,LIN,fp)!=NULL){
  if(!strncmp(line,"#\tAA",4)){
   ftype=1;//*.spd3
   break;
  }
  //if(!sscanf(line,">%s",seq->id)){
  //if(!strncmp(line,">",1)){
  // ftype=0;//fasta format
  // puts("##FASTA file");
  // break;
  //}
  //spot1d-single format
  if(!strncmp(line,",AA,",4)){
   ftype=2;//csv format
   puts("##CVS file");
   break;
  }
  //DSSP format
  if(!strncmp(line,"====",4)){
   ftype=3;//dssp format
   puts("##DSSP file");
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
	 seq->Pss[i][0]=seq->Pss[i][1]=seq->Pss[i][2]=0.00;//HEC
        }
   //printf("%d %c %d\n",i,seq->seq_txt[i],seq->seq_code[i]);
  }
 }
 if(ftype==1){
  Nch=-1;
  puts("#spd3/spot1d file");
  if((fp=fopen(fname,"r"))==NULL)
   return true;
  while(fgets(line,LIN,fp)!=NULL){
   if(!strncmp(line,"#",1)){
     strncpy(seq->id,"Unknown",7);
     seq->id[8]='\0';
     printf("#ID= '%s'\n",seq->id);
     //chain-break
     Nch++;
   }else{
    float Pc,Pe,Ph;
//#	AA	SS3	SS8	ASA 	HSEa-u	HSEa-d	CN13	theta	tau 	phi 	psi 	P(3-C)	P(3-E)	P(3-H)	P(8-C)	P(8-S)	P(8-T)	P(8-H)	P(8-G)	P(8-I)	P(8-E)	P(8-B)
    sscanf(line,"%d\t%c\t%c\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%f\t%f\t%f",&i,&aa,&ss,buf,buf,buf,buf,buf,buf,buf,buf,buf,&Pc,&Pe,&Ph);
    seq->seq_code[seq->len]=A2int(aa);
    seq->seq_txt[seq->len]=aa;
    	seq->Pss[seq->len][0]=Ph*0.01;
    	seq->Pss[seq->len][1]=Pe*0.01;
    	seq->Pss[seq->len][2]=Pc*0.01;
	seq->CID[seq->len]=Nch;
	printf("SS=%d %f %f %f\n",i,Ph,Pe,Pc);
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
 if(ftype==2){
  Nch=-1;
  puts("#spot1d-single file");
  if((fp=fopen(fname,"r"))==NULL)
   return true;
  while(fgets(line,LIN,fp)!=NULL){
   if(!strncmp(line,",AA,",4)){
     strncpy(seq->id,"Unknown",7);
     seq->id[8]='\0';
     printf("#ID= '%s'\n",seq->id);
     //chain-break
     Nch++;
   }else{
    float Pc,Pe,Ph;
//,AA,SS3,SS8,ASA,HseU,HseD,CN,Psi,Phi,Theta,Tau,P3C,P3E,P3H,P8C,P8S,P8T,P8H,P8G,P8I,P8E,P8B
    sscanf(line,"%d,%c,%c,%c,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f",&i,&aa,&ss,buf,&fbuf,&fbuf,&fbuf,&fbuf,&fbuf,&fbuf,&fbuf,&fbuf,&Pc,&Pe,&Ph);
    seq->seq_code[seq->len]=A2int(aa);
    seq->seq_txt[seq->len]=aa;
    	seq->Pss[seq->len][0]=Ph;
    	seq->Pss[seq->len][1]=Pe;
    	seq->Pss[seq->len][2]=Pc;
	seq->CID[seq->len]=Nch;
	printf("#SS=%d %f %f %f\n",i,Ph,Pe,Pc);
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
 if(ftype==3){
  Nch=-1;
  int start_flag=0;
  puts("#DSSP file");
  if((fp=fopen(fname,"r"))==NULL)
   return true;
  while(fgets(line,LIN,fp)!=NULL){
   if(!strncmp(line,"  #  RESIDUE",12)){
     strncpy(seq->id,"Unknown",7);
     seq->id[8]='\0';
     printf("#ID= '%s'\n",seq->id);
     //chain-break
     Nch++;
     start_flag=1;
   }else{
    if(start_flag==0)
	continue;
    float Pc,Pe,Ph;
    Pc=Pe=Ph=0;
    aa = line[13];
    ss=line[16];
    seq->seq_code[seq->len]=A2int(aa);
    seq->seq_txt[seq->len]=aa;
	if(ss=='H'||ss=='I'||ss=='G')
	 Ph=1.0;
	else if(ss=='E'||ss=='B')
	 Pe=1.0;
	else
	 Pc=1.0;
    	seq->Pss[seq->len][0]=Ph;
    	seq->Pss[seq->len][1]=Pe;
    	seq->Pss[seq->len][2]=Pc;
	seq->CID[seq->len]=Nch;
	printf("#SS=%c %f %f %f\n",ss,Ph,Pe,Pc);
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
 //int flen=len+4;//???
 int flen=len;//???
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
   fg[n].Pss[j][0]=seq->Pss[i+j][0];
   fg[n].Pss[j][1]=seq->Pss[i+j][1];
   fg[n].Pss[j][2]=seq->Pss[i+j][2];
   printf("%3d",fg[n].seq[j]);
   printf("(%1d,[%2.0f,%2.0f,%2.0f])",fg[n].ss[j],fg[n].Pss[j][0]*10,fg[n].Pss[j][1]*10,fg[n].Pss[j][2]*10);
  }
  printf(" chain=%d\n",Nc);
  fg[n].cid=Nc;
  fg[n].l=flen;
  fg[n].pos=i;
  n++;
 }
 return n;
}

bool AssignSSPtoPDB(SEQ *s,PDB *p){
 int *ali1=(int *)malloc(sizeof(int)*p->NumOfRes);
 int *ali2=(int *)malloc(sizeof(int)*s->len);
 int *gali=(int *)malloc(sizeof(int)*(s->len+p->NumOfRes)*2);
 int glen;
 float *Smtx=(float *)malloc(sizeof(float)*(s->len+1)*(p->NumOfRes+1));
 DPMTX *Dmtx=(DPMTX *)malloc(sizeof(DPMTX)*(s->len+1)*(p->NumOfRes+1));
 int n1=p->NumOfRes;
 int n2=s->len;
 //malloc
 if((p->SSP=(float*)malloc(sizeof(float)*p->NumOfAtom*3))==NULL)
  return true;

 //Residue Table
 for(int i=0;i<p->NumOfRes;i++){
  int a = p->ResOnAtom[i];
  int aa = p->TypeResId[i];
	for(int j=0;j<s->len;j++){
	 	if(aa == s->seq_code[j]){
		 Smtx[GRID2D(n1,n2,i,j)]=1.0;
	 	 //printf("[%d %d] [%d %d] %d\n",i,aa,j,s->seq_code[j],GRID2D(n1,n2,i,j));
		}else{
		 Smtx[GRID2D(n1,n2,i,j)]=-1.0;
		}
	}
 }
 puts("#DP..");
 //dp(Dmtx,Smtx,-100.0,0,ali1,p->NumOfRes,ali2,s->len,gali,&glen);
 dp(Dmtx,Smtx,-1.0,0,ali1,p->NumOfRes,ali2,s->len,gali,&glen);

 int spos;
 int Miss=0;
 int Nali=0;
 for(int i=0;i<p->NumOfRes;i++){
  spos=ali1[i];
  if(spos!=-1){
   Nali++;
   printf("PDB %d %d %d\n",i,p->TypeResId[i],s->seq_code[spos]);
   //printf("PDB %d %d %d\n",i,p->TypeResId[i],spos);
   if(p->TypeResId[i]!=s->seq_code[spos])
    Miss++;
  }else{
   printf("PDB %d %d %d\n",i,p->TypeResId[i],-1);
  }
 }
 printf("##Num of Aligned residues: %d/%d Unaligned %d\n",Nali,p->NumOfRes,p->NumOfRes-Nali);
 printf("##Num of Miss aligned AA: %d\n",Miss);
 if(Nali != p->NumOfRes){
  printf("##WARRNING SEQUENCE != PDB\n");
  //return true;
 }

 for(int i=0;i<p->NumOfAtom;i++){
  int j=p->AtomOnRes[i];
  
  spos=ali1[j];
  if(spos!=-1){
   p->SSP[3*i  ]=s->Pss[spos][0];
   p->SSP[3*i+1]=s->Pss[spos][1];
   p->SSP[3*i+2]=s->Pss[spos][2];
  }else{
   p->SSP[3*i  ]=0.0;
   p->SSP[3*i+1]=0.0;
   p->SSP[3*i+2]=0.0;
  }
 }

 return false;
}

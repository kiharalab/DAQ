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
#include "mrcfft.h"
//#include "scoring.h"

#define PDB_STRLEN 55

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

int AtomMap(MRC *, PDB *,DSSP *, int);
int AtomVox(MRC *, PDB *,DSSP *,VOXEL *, int,double);
int AssignAtom(MRC *, PDB *,VOXEL *,double);
int CompPdbDssp(PDB *, DSSP *);
double FindTopX(MRC *,double );

int main(int argc, char **argv)
{
 double t1=gettimeofday_sec();
 double t4;
 MRC mrc;
 PDB pdb;
 DSSP dssp;
 VOXEL *vox;

 int Natm=0;
 //Get Options
 if(chkcmdline(argc,argv,&cmd)==FALSE)
  return(0);
/*
 //Set threads
 if(cmd.Nthr < omp_get_num_procs()){
  omp_set_num_threads(cmd.Nthr);
 }else{
  omp_set_num_threads(omp_get_num_procs());
 }
 */
 
 if(readmrc(&mrc,cmd.filename))
  return(0);

 if((vox=(VOXEL *)malloc(sizeof(VOXEL)*mrc.xdim*mrc.ydim*mrc.zdim))==NULL)
  return 0;
 //init
 for(int i=0;i<mrc.xdim*mrc.ydim*mrc.zdim;i++){
	 vox[i].AtomId=-1;
	 vox[i].d2=100.00;
 }

 //NCSTART is wrong?
 if(cmd.IgnoreNC==false){
  mrc.orgxyz[0]=mrc.orgxyz[0]+mrc.ncstart*mrc.widthx;
  mrc.orgxyz[1]=mrc.orgxyz[1]+mrc.nrstart*mrc.widthy;
  mrc.orgxyz[2]=mrc.orgxyz[2]+mrc.nsstart*mrc.widthz;
  printf("#New orgXYZ: %f %f %f\n",mrc.orgxyz[0],mrc.orgxyz[1],mrc.orgxyz[2]);
 }

 //Density Filter and normalization*****
 for(int i=0;i<mrc.xdim*mrc.ydim*mrc.zdim;i++)
	 if(mrc.dens[i]<cmd.th1)
		mrc.inside[i]=false;
	 else
		mrc.inside[i]=true;

 float TopCut=FindTopX(&mrc,cmd.Ubound);
 //norm
 for(int i=0;i<mrc.xdim*mrc.ydim*mrc.zdim;i++){
	 mrc.dens[i]=mrc.dens[i]/(TopCut);
	 if(mrc.dens[i]>1.00)
		mrc.dens[i]=1.000;
	 if(mrc.dens[i]<0.00)
		mrc.dens[i]=0.000;
 }

 //Use PDB and DSSP
 if(cmd.Mode==1){
  printf("#USE PDB FILE\n");
 	Natm=CountAtom(cmd.file1);
 	if((MallocPdb(&pdb,Natm))==-1)
	 return(0);
 	if(readpdb(&pdb,cmd.file1,Natm))
 	 return(0);

 	//if(ReadDssp(&dssp,cmd.file2))
 	// return(0);

 	//if(CompPdbDssp(&pdb,&dssp))
	// return(0);
 	//Assign Atom to voxel
 	AssignAtom(&mrc,&pdb,vox,cmd.NoP2);
	puts("##ASSIGN ATOM DONE");
 }else{
  printf("#SCAN MODE\n");
  pdb.NumOfAtom=0;
  pdb.NumOfRes=0;

  	for(int i=0;i<mrc.xdim*mrc.ydim*mrc.zdim;i++){
         vox[i].AtomId=-1;
         vox[i].d2=-1.00;//No distance data
 	}
 }

 AtomVox(&mrc,&pdb,&dssp,vox,cmd.Nvox,cmd.r);
 t4=gettimeofday_sec();
 printf("#FINISHED TOTAL TIME= %f\n",t4-t1);
 return 0;
}

//Low to high
int cmp_float(const void *a, const void *b){
	float A =*(const float*)a;
	float B =*(const float*)b;
	if(A < B) return  -1;
	if(A > B) return   1;
       	return 0;
}

double FindTopX(MRC *m,double c){
 int n;
 int Nvox=m->xdim*m->ydim*m->zdim;
 float *tbl;
 int Ntbl=0;
 int Count[200];
 float LogCount[200];

 if((tbl=malloc(sizeof(float)*Nvox))==NULL)
	 return -1.00;

 for(int i=0;i<Nvox;i++){
	 if(m->dens[i]>0.00){
		 tbl[Ntbl]=m->dens[i];
		 Ntbl++;
	 }
 }
 qsort(tbl,Ntbl,sizeof(float),cmp_float);//Low to High
 double tic=m->dmax/200.00;
 int j=0;
 for(int i=0;i<200;i++){
	Count[i]=0;
	LogCount[i]=0;
 }
 for(int i=0;i<Nvox;i++){
  if((double)(((j+1)*tic))>=tbl[i]){
	  Count[j]++;
  }
  else{
	  i--;
	  j++;
	  continue;
  }
 }
 //log
 //
 double Sum=0;
 for(int i=0;i<200;i++){
	if(Count[i]>0){
	 LogCount[i]=log(Count[i]);
	 Sum+=LogCount[i];
	}
 }
 double SumCut=Sum*c;
 printf("#Top(R %f)=Area(%.3f/%.3f)\n",c,SumCut,Sum);
 Sum=0;
 double CutOff=0;
 for(int i=0;i<200;i++){
   if(Count[i]>0){
	Sum+=LogCount[i];
	if(Sum>=SumCut){
		CutOff=tic*i;
		break;
	}
   }
 }
 printf("#CutOff= %f\n",CutOff);
 return CutOff;
}

int CompPdbDssp(PDB *p, DSSP *d){
	double d2;
	int atm;
    puts("#Cheking PDB and DSSP..");
	//if(d->Nres != p->NumOfRes){
	//	printf("Inconsistent Nres: %d != %d\n", p->NumOfRes,d->Nres );
	//	return -1;
	//}
	//corrdinate check
	for(int i=0;i<p->NumOfRes;i++){
	    for(int j=0;j<d->Nres;j++){
	    d2=(p->CAxyz[i][0]-d->xyz[j][0])*(p->CAxyz[i][0]-d->xyz[j][0])
	    +(p->CAxyz[i][1]-d->xyz[j][1])*(p->CAxyz[i][1]-d->xyz[j][1])
	    +(p->CAxyz[i][2]-d->xyz[j][2])*(p->CAxyz[i][2]-d->xyz[j][2]);
	    if(d2<1){
            p->DSSPmapper[i]=j;
	    }

	    }

	  //printf("d= %f\n",d2);
	  //printf("%f %f %f\n",p->CAxyz[i][0],p->CAxyz[i][1],p->CAxyz[i][2]);
	  //printf("%f %f %f\n",d->xyz[i][0],d->xyz[i][1],d->xyz[i][2]);
	  //if(d2>1.0){
	  // printf("##RES %d distance= %f\n",i,sqrt(d2));
	  // return -1;
	  //}
	}
 puts("#Cheking PDB and DSSP..DONE");
 return 0;
}




void VoxPos(MRC *m,double cd[3],int pos[3]){

	pos[0]=(int)((cd[0]-m->orgxyz[0])/m->widthx);
	pos[1]=(int)((cd[1]-m->orgxyz[1])/m->widthy);
	pos[2]=(int)((cd[2]-m->orgxyz[2])/m->widthz);

	//printf("COD %f %f %f\n",cd[0],cd[1],cd[2]);
	//printf("POS %d %d %d\n",pos[0],pos[1],pos[2]);

}

void ShowVox(MRC *m,int pos[3],int N){

	int st[3],ed[3];
	int idx;
	double d;

	st[0]=pos[0]-N;
	st[1]=pos[1]-N;
	st[2]=pos[2]-N;

	ed[0]=pos[0]+N;
	ed[1]=pos[1]+N;
	ed[2]=pos[2]+N;

	for(int x=st[0];x<=ed[0];x++){
	for(int y=st[1];y<=ed[1];y++){
	for(int z=st[2];z<=ed[2];z++){
		if(x<0||y<0||z<0||x>=m->xdim||y>=m->ydim||z>=m->zdim)
			d=0.000;
		else
		 	d=m->dens[m->xdim*m->ydim*z+m->xdim*y+x];
		printf("%.3f,",d);
	}}}
	printf("\n");
	//printf("%d %d %d\n",pos[0],pos[1],pos[2]);
	//printf("Center= %f\n",m->dens[m->xdim*m->ydim*pos[2]+m->xdim*pos[1]+pos[0]]);

}

double AveDensVox(MRC *m,int pos[3],int N){

	int st[3],ed[3];
	int idx;
	double d=0.0000;
	int Nv=(2*N+1)*(2*N+1)*(2*N+1);

	st[0]=pos[0]-N;
	st[1]=pos[1]-N;
	st[2]=pos[2]-N;

	ed[0]=pos[0]+N;
	ed[1]=pos[1]+N;
	ed[2]=pos[2]+N;

	for(int x=st[0];x<=ed[0];x++){
	for(int y=st[1];y<=ed[1];y++){
	for(int z=st[2];z<=ed[2];z++){
		if(x<0||y<0||z<0||x>=m->xdim||y>=m->ydim||z>=m->zdim)
			continue;
		else
		 	d+=m->dens[m->xdim*m->ydim*z+m->xdim*y+x];
		//printf("%.3f,",d);
		//
	}}}
	//printf("\n");
	//printf("Ave= %f\n",d/(double)(Nv));
	return(d/(double)(Nv));
}


int AtomMap(MRC *m, PDB *p,DSSP *d, int N){



	int i,r,Nerror,res;
	double dis;
	r=-1;
	Nerror=0;
	printf("N= %d\n",p->NumOfAtom);
	for(i=0;i<p->NumOfAtom;i++){
		int pos[3];
		if(p->TypeAtomId[i]==2){//CA
		 r++;
		 //check
		 dis=(p->xyz[i][0]-d->xyz[r][0])*(p->xyz[i][0]-d->xyz[r][0])
		  +(p->xyz[i][1]-d->xyz[r][1])*(p->xyz[i][1]-d->xyz[r][1])
		  +(p->xyz[i][2]-d->xyz[r][2])*(p->xyz[i][2]-d->xyz[r][2]);
		 if(dis>1.00){
			 Nerror++;
			 if(Nerror > 3){
			  printf("Too Many Errors: DSSP or PDB file is wrong\n");
				 return -1;
			 }
		 }
		}
	 //Atom Type, Main-chain, Atom-Type, SS, etc.
	 res=p->AtomOnRes[i];
	 printf("AtomID: %d,",i);
	 printf("AA: %d, Atom: %d, SS: %d,",p->TypeResId[res],p->TypeAtomId[i],d->ss[res]);
	 VoxPos(m,p->xyz[i],pos);
	 printf(" Dens:");
	 ShowVox(m,pos,N);
	}
	return 0;
}

int Pos2Idx(MRC *m,int pos[3]){
	int idx;
	if(pos[0]<0||pos[1]<0||pos[2]<0)
		return -1;
	if(pos[0]>=m->xdim||pos[1]>=m->ydim||pos[2]>=m->zdim)
		return -1;

       	idx=m->xdim*m->ydim*pos[2]+m->xdim*pos[1]+pos[0];

	return idx;
}

int Pos2Cd(MRC *m,int pos[3],double cd[3]){
	cd[0]=m->orgxyz[0]+pos[0]*m->widthx;
	cd[1]=m->orgxyz[1]+pos[1]*m->widthy;
	cd[2]=m->orgxyz[2]+pos[2]*m->widthz;
}


int AssignAtom(MRC *m,PDB *p,VOXEL *v,double r){
	int pos1[3],pos2[3],pos3[3],idx;
	double cd1[3],cd2[3],cd3[3];
	double r2=r*r;
	double d2;

	//check COG
	double cog_pdb[3],cog_map[3],cog_map2[3];
	cog_pdb[0]=0;
	cog_pdb[1]=0;
	cog_pdb[2]=0;
	cog_map2[0]=0;
	cog_map2[1]=0;
	cog_map2[2]=0;
	for(int i=0;i<p->NumOfRes;i++){
		cog_pdb[0]+=p->CAxyz[i][0];
		cog_pdb[1]+=p->CAxyz[i][1];
		cog_pdb[2]+=p->CAxyz[i][2];
	}
	cog_pdb[0]/=(double)p->NumOfRes;
	cog_pdb[1]/=(double)p->NumOfRes;
	cog_pdb[2]/=(double)p->NumOfRes;
	cog_map[0]=m->orgxyz[0]+0.5*m->xdim*m->widthx;
	cog_map[1]=m->orgxyz[1]+0.5*m->xdim*m->widthy;
	cog_map[2]=m->orgxyz[2]+0.5*m->xdim*m->widthz;
	
	double diff=(cog_pdb[0]-cog_map[0])*(cog_pdb[0]-cog_map[0])
		+(cog_pdb[1]-cog_map[1])*(cog_pdb[1]-cog_map[1])
		+(cog_pdb[2]-cog_map[2])*(cog_pdb[2]-cog_map[2]);
	diff=sqrt(diff);
	double ori_len=sqrt((m->xlen*m->xlen)*(m->ylen*m->ylen)*(m->zlen*m->zlen));
	ori_len*=0.5;
	printf("#COG_DIFF: %.1f Rate: %f diff\n",diff,diff/ori_len);
	int count_assign;
	count_assign=0;


	for(int i=0;i<p->NumOfAtom;i++){
		cd1[0]=p->xyz[i][0]-r;
		cd1[1]=p->xyz[i][1]-r;
		cd1[2]=p->xyz[i][2]-r;
		cd2[0]=p->xyz[i][0]+r;
		cd2[1]=p->xyz[i][1]+r;
		cd2[2]=p->xyz[i][2]+r;
	 	VoxPos(m,cd1,pos1);
	 	VoxPos(m,cd2,pos2);
	    //printf("#pos1: %d %d %d\n",pos1[0],pos1[1],pos1[2]);
	    //printf("#pos2: %d %d %d\n",pos2[0],pos2[1],pos2[2]);
	 //fill
		for(pos3[0]=pos1[0];pos3[0]<=pos2[0];pos3[0]++){
		for(pos3[1]=pos1[1];pos3[1]<=pos2[1];pos3[1]++){
		for(pos3[2]=pos1[2];pos3[2]<=pos2[2];pos3[2]++){
			idx=Pos2Idx(m,pos3);
			if(idx==-1)
				continue;
			Pos2Cd(m,pos3,cd3);
			d2=	 (p->xyz[i][0]-cd3[0])*(p->xyz[i][0]-cd3[0])
				+(p->xyz[i][1]-cd3[1])*(p->xyz[i][1]-cd3[1])
				+(p->xyz[i][2]-cd3[2])*(p->xyz[i][2]-cd3[2]);

			//Update
			if(v[idx].AtomId==-1||v[idx].d2>d2){
			 v[idx].d2=d2;
			 v[idx].AtomId=i;
			 //printf("#Update: %d %d %f\n",idx,i,d2);
			 count_assign+=1;
			 cog_map2[0]+=pos3[0];
			 cog_map2[1]+=pos3[1];
			 cog_map2[2]+=pos3[2];
			}
		}}}

	
	}
	for(int k=0;k<3;k++){
	    cog_map2[k]/=count_assign;
	}
	cog_map2[0]=m->orgxyz[0]+m->widthx*cog_map2[0];
	cog_map2[1]=m->orgxyz[1]+m->widthy*cog_map2[1];
	cog_map2[2]=m->orgxyz[2]+m->widthy*cog_map2[2];
	diff=(cog_pdb[0]-cog_map2[0])*(cog_pdb[0]-cog_map2[0])
		+(cog_pdb[1]-cog_map2[1])*(cog_pdb[1]-cog_map2[1])
		+(cog_pdb[2]-cog_map2[2])*(cog_pdb[2]-cog_map2[2]);
	diff=sqrt(diff);
	printf("#PDB center: %.2f, %.2f, %.2f; Map center: %.2f, %.2f, %.2f\n",cog_pdb[0],cog_pdb[1],cog_pdb[2],cog_map2[0],cog_map2[1],cog_map2[2]);
	printf("#COG_DIF2: %.1f\n",diff);
	printf("#In total we assigned %d atom voxels\n",count_assign);
}

int AtomVox(MRC *m, PDB *p,DSSP *d,VOXEL *v, int N,double r){
	int i,Nerror,res,dssp_res,pos[3];
	int aid,idx;
	int atm,aa,ss,acc;
	int inside;
	double kap,alpha,phi,psi;
	double dis;
	double r1=r*r;
	double r2=(cmd.NoP1)*(cmd.NoP1);
	double r3=(cmd.NoP2)*(cmd.NoP2);
	printf("#NoP %f ~ %f\n",cmd.NoP1,cmd.NoP2);
	printf("#Started to search x %d y %d z %d\n",m->xdim,m->ydim,m->zdim);
	r=-1;
	Nerror=0;
    int count_ignore,count_close;
    count_ignore=0;
    count_close=0;
	for(pos[0]=0;pos[0]<m->xdim;pos[0]+=cmd.slide){
	for(pos[1]=0;pos[1]<m->ydim;pos[1]+=cmd.slide){
	for(pos[2]=0;pos[2]<m->zdim;pos[2]+=cmd.slide){
		idx=Pos2Idx(m,pos);
		aid=-9;
		inside=0;
		//ignore too far voxels
		if(v[idx].d2 > r3){
		    count_ignore+=1;
		    //if (v[idx].AtomId>0)printf("#Checked atom: %d %d %f\n",idx,v[idx].AtomId,v[idx].d2);
			continue;}
        count_close+=1;
		//No protein region r2<d2<=r3
		if(v[idx].d2 > r2 && v[idx].d2 <= r3 ){
			aid=-1000;
			res=-1000;
			aa=-1000;atm=-1000;ss=-1000;acc=-1000;kap=-1000;alpha=-1000;phi=-1000;psi=-1000;
		}else if(v[idx].d2 <= r1 && v[idx].d2 >=0.000){
		 aid=v[idx].AtomId;
		 res=p->AtomOnRes[aid];
		 aa=p->TypeResId[res];
		 atm=p->TypeAtomId[aid];
		/*
		 dssp_res=p->DSSPmapper[res];
		 ss=d->ss[dssp_res];
		 acc=d->acc[dssp_res];
		 kap=d->kappa[dssp_res];
		 alpha=d->alpha[dssp_res];
		 phi=d->phi[dssp_res];
		 psi=d->psi[dssp_res];
		*/
		}else if(m->inside[idx]==true){//inside of contour level
		 	//aid=-1;
			//res=-1;
			//aa=-1;atm=-1;ss=-1;acc=-1;kap=0;alpha=0;phi=0;psi=0;
			aid=-1000;
			res=-1000;
			aa=-1000;atm=-1000;ss=-1000;acc=-1000;kap=-1000;alpha=-1000;phi=-1000;psi=-1000;
		}

		if(m->inside[idx]==true)
			inside=1;
		if(aid==-9 && v[idx].d2 >=0.000 )//outside, no atom
			continue;
		if(AveDensVox(m,pos,N)<cmd.Lbound)//too low density value
			continue;
		if(v[idx].d2>0)
		 printf("AtomID: %d, ResID: %d, INSIDE: %d, Vpos: %d,%d,%d, dist: %.2f",aid,res,inside,pos[0],pos[1],pos[2],sqrt(v[idx].d2));
		else
		 printf("AtomID: %d, ResID: %d, INSIDE: %d, Vpos: %d,%d,%d, dist: %.2f",aid,res,inside,pos[0],pos[1],pos[2],v[idx].d2);
	 	printf(" AA: %d, ATOM: %d, SS: %d ACC: %d,",aa,atm,ss, acc);
	 	printf(" KAPPA: %.1f, ALPHA %.1f, PHI %.1f, PSI: %.1f,",kap,alpha,phi,psi);
	 	printf(" DENS:");
	 	ShowVox(m,pos,N);
	}}}
	printf("#In total %d vs %d non-protein examples are ignored because too far\n",count_ignore,count_close);
	return 0;
}

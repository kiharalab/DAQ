#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "struct.h"
#include "dp.h"
#include "mrc.h"
#define FRG_SCO_CUT 20.0 //Minimum DP score

extern CMD cmd;//Important



typedef struct{  
        float sco;
	DPMTX *Dmtx;
	float *Smtx,*SmtxRv, *dtbl,*dtblR;
	float *Smtx_ori,*SmtxRv_ori;
	int n1,n2,Lgali;
	int *ali1,*ali2,*gali;
	NODE **model;//CA model
	NODE **path;
} DP_MEMORY;


//For ShowPath
typedef struct{
 NODE *model[2000];//node
 int L;//Length
 int seq[2000];
 float sco;
} SUBOPT_ALI;




double QualityTreeDP(GRAPH *, TREE *, SEQFG *,DP_MEMORY *,bool);
double QualityTreeDPSubOpt(GRAPH *, TREE *, SEQFG *,DP_MEMORY *,SUBOPT_ALI *);
double QualityPathDP(GRAPH *,SEQFG *,DP_MEMORY *,bool);

int permuted_index(int ordermode,int count, unsigned nc, unsigned nr, unsigned ns) {
        unsigned ic, ir, is;
        unsigned long ncr, q;

        ncr = nc*nr;
        is = count / ncr;
        q = count - is*ncr;
        ir = q / nc;
        ic = q - ir*nc;

        switch(ordermode) {
                case 1:
                        return ic+ir*nc+is*nc*nr;
                case 2:
                        return ic+is*nc+ir*nc*ns;
                case 3:
                        return ir+ic*nr+is*nr*nc;
                case 4:
                        return is+ic*ns+ir*ns*nc;
                case 5:
                        return ir+is*nr+ic*nr*ns;
                case 6:
                        return is+ir*ns+ic*ns*nr;
                default:
                        exit(0);
        }
}


bool readmrc(MRC *mrc, char *filename){
 FILE *fpin;
 int ibuf,res,mode;
 int mapr,mapc,maps;
 int xdim,ydim,zdim;
 int nsymbt;
 float fbuf,*tmp_dens;

 if((fpin=fopen(filename,"rb")) == NULL){ 
  fprintf(stderr,"Can't open %s\n",filename); 
  return(true); 
 }
 
 //res=fread(&ibuf,sizeof(int),1,fpin);
 res=fread(&(mrc->xdim),sizeof(int),1,fpin);
 res=fread(&(mrc->ydim),sizeof(int),1,fpin);
 res=fread(&(mrc->zdim),sizeof(int),1,fpin);
 xdim=mrc->xdim;
 ydim=mrc->ydim;
 zdim=mrc->zdim;
 printf("#reading %s\n",filename);
 printf("#XYZ dim: %d %d %d\n",mrc->xdim,mrc->ydim,mrc->zdim);
 //cmd->xdim=ibuf;
 res=fread(&(mode),sizeof(int),1,fpin);
 //ignore
/*
 res=fread(&(ibuf),sizeof(int),1,fpin);
 res=fread(&(ibuf),sizeof(int),1,fpin);
 res=fread(&(ibuf),sizeof(int),1,fpin);
*/
 res=fread(&(mrc->ncstart),sizeof(int),1,fpin);
 printf("#NCSTART %d\n",mrc->ncstart);
 res=fread(&(mrc->nrstart),sizeof(int),1,fpin);
 printf("#NRSTART %d\n",mrc->nrstart);
 res=fread(&(mrc->nsstart),sizeof(int),1,fpin);
 printf("#NSSTART %d\n",mrc->nsstart);

 //mxyz
 res=fread(&(mrc->mx),sizeof(int),1,fpin);
 res=fread(&(mrc->my),sizeof(int),1,fpin);
 res=fread(&(mrc->mz),sizeof(int),1,fpin);
 printf("#MXYZ: %d %d %d\n",mrc->mx,mrc->my,mrc->mz);

 //xyz len
 res=fread(&(mrc->xlen),sizeof(float),1,fpin);
 res=fread(&(mrc->ylen),sizeof(float),1,fpin);
 res=fread(&(mrc->zlen),sizeof(float),1,fpin);
 printf("#LenXYZ: %f %f %f\n",mrc->xlen,mrc->ylen,mrc->zlen);


 //abg
 res=fread(&(mrc->alpha),sizeof(float),1,fpin);
 res=fread(&(mrc->beta),sizeof(float),1,fpin);
 res=fread(&(mrc->gamma),sizeof(float),1,fpin);
 printf("#abg: %f %f %f\n",mrc->alpha,mrc->beta,mrc->gamma);

 //map crs
 res=fread(&(mrc->mapc),sizeof(int),1,fpin);
 res=fread(&(mrc->mapr),sizeof(int),1,fpin);
 res=fread(&(mrc->maps),sizeof(int),1,fpin);
 mapc=mrc->mapc;
 mapr=mrc->mapr;
 maps=mrc->maps;
 printf("#crs: %d %d %d\n",mrc->mapc,mrc->mapr,mrc->maps);

 res=fread(&(mrc->dmin),sizeof(float),1,fpin);
 res=fread(&(mrc->dmax),sizeof(float),1,fpin);
 res=fread(&(mrc->dmean),sizeof(float),1,fpin);
 res=fread(&(mrc->ispg),sizeof(int),1,fpin);
 printf("#dmax,dmin, dmean, ispg: %f %f %f %d\n",mrc->dmax,mrc->dmin,mrc->dmean,mrc->ispg);

 //93-96
 res=fread(&(mrc->nsymbt),sizeof(int),1,fpin);
 nsymbt=mrc->nsymbt;
 //ignore 25
 //for(int i=0;i<25;i++)
 // res=fread(&(ibuf),sizeof(int),1,fpin);
 //97-196 Extra
 fseek(fpin,4*25,SEEK_CUR);

 //197-208 ORIGIN
 res=fread(&(mrc->orgxyz),sizeof(float),3,fpin);
 printf("#orgXYZ: %f %f %f\n",mrc->orgxyz[0],mrc->orgxyz[1],mrc->orgxyz[2]);
 
 //ignore MAP 209-212
 char text[4];
 res=fread(text,sizeof(char)*4,1,fpin);
 //fseek(fpin,4,SEEK_CUR);

 //printf("#%s",text);
 if(strncmp(text,"MAP",3)){
  printf("Format Error!!!\n");
  return true;
 }
 char machst[4];
 bool swap=false;
 res=fread(machst,sizeof(char)*4,1,fpin);

 if(machst[0] == 0x44){
  swap = false;
  printf("#little-endian mode\n");
 }else if(machst[0] == 0x11){
  swap = true;
  printf("big-endian mode\n");
 }

 if(swap){
  printf("WARNING THIS IS BIG-ENDIAN FILES!!\n");
  return true;
 }

 int ordermode=0;
 	if(mapc==1 && mapr==2 && maps==3) {
                ordermode = 1;
        }
        else if(mapc==1 && mapr==3 && maps==2) {
                ordermode = 2;
        }
        else if(mapc==2 && mapr==1 && maps==3) {
                ordermode = 3;
        }
        else if(mapc==2 && mapr==3 && maps==1) {
                ordermode = 4;
        }
        else if(mapc==3 && mapr==1 && maps==2) {
                ordermode = 5;
        }
        else if(mapc==3 && mapr==2 && maps==1) {
                ordermode = 6;
        }
        else if(ordermode == 0) {
         printf("Input file gives malformed dimension ordering.");
         return true;
        }
  printf("#Order Mode= %d\n",ordermode);
  mrc->NumVoxels = mrc->xdim*mrc->ydim*mrc->zdim;
 printf("#Nvoxels= %d\n",mrc->NumVoxels);

 if((mrc->dens=(float*)malloc(sizeof(float)*mrc->NumVoxels))==NULL)
  return true;

 //fin.ignore(4*(256-54)+nsymbt);
 fseek(fpin,4*(256-54)+nsymbt,SEEK_CUR);
 
 	switch(mode) {
 	 case 0: // char - converted to float, testing for signed-ness
	   printf("Cannot read mode 0 mrc file\n");
           break;
         case 1: // 16-bit float
	  printf("#Reading 16-bit mrc file\n");
          for(int i = 0; i<mrc->NumVoxels; ++i)
	   res=fread(&(mrc->dens[permuted_index(ordermode,i,xdim,ydim,zdim)]),2,1,fpin);
           break;
         case 2: // 32-bit float
	  printf("#Reading 32-bit mrc file\n");
          for(int i = 0; i<mrc->NumVoxels; ++i){
	   res=fread(&(mrc->dens[permuted_index(ordermode,i,xdim,ydim,zdim)]),sizeof(float),1,fpin);
	  }
          break;
         default:
          printf("Unknown floating-point mode specified.");
	  return true;
        }

 	mrc->widthx = mrc->xlen / (double) mrc->mx;
        mrc->widthy = mrc->ylen / (double) mrc->my;
        mrc->widthz = mrc->zlen / (double) mrc->mz;

	if(fabs(mrc->widthx - mrc->widthy)>0.000001 || 
	 fabs(mrc->widthx -  mrc->widthz)>0.000001 ||
	 fabs(mrc->widthy -  mrc->widthz)>0.000001){

	 printf("#ERROR: grid sizes are different %f %f %f\n",
	  mrc->widthx,mrc->widthy,mrc->widthz);
	 printf("PLEASE USE CUBIC MRC MAP DATA\n");
	 return true;
	}

	int nx,ny,nz;
	switch(ordermode){
		case 1:
			nx = xdim; ny = ydim; nz = zdim;
			break;
		case 2:
			nx = xdim; ny = zdim; nz = ydim;
			break;
		case 3:
			nx = ydim; ny = xdim; nz = zdim;
			break;
		case 4:
			nx = zdim; ny = xdim; nz = ydim;
			break;
		case 5:
			nx = ydim; ny = zdim; nz = xdim;
			break;
		case 6:
			nx = zdim; ny = ydim; nz = xdim;
			break;
		default:
			printf("Input file gives malformed dimension ordering.");
			return true;
	}

 mrc->xdim = nx;
 mrc->ydim = ny;
 mrc->zdim = nz;

 printf("#XYZ dim: %d %d %d\n",mrc->xdim,mrc->ydim,mrc->zdim);

 fclose(fpin);


 if(mrc->ncstart!=0||mrc->nrstart!=0||mrc->nsstart!=0){
  mrc->orgxyz[0]=mrc->orgxyz[0]+mrc->ncstart*mrc->widthx;
  mrc->orgxyz[1]=mrc->orgxyz[1]+mrc->nrstart*mrc->widthy;
  mrc->orgxyz[2]=mrc->orgxyz[2]+mrc->nsstart*mrc->widthz;
 }

 return false;
}

/*
	xdim = read_int_istream(fin,swap);
        ydim = read_int_istream(fin,swap);
        zdim = read_int_istream(fin,swap);

        mode = read_int_istream(fin,swap);

        fin.ignore(4*3);

        mx = read_int_istream(fin,swap);
        my = read_int_istream(fin,swap);
        mz = read_int_istream(fin,swap);

        xlen = read_float_istream(fin,swap);
        ylen = read_float_istream(fin,swap);
        zlen = read_float_istream(fin,swap);

        alpha = read_float_istream(fin,swap);
        beta = read_float_istream(fin,swap);
        gamma = read_float_istream(fin,swap);

        mapc = read_int_istream(fin,swap);
        mapr = read_int_istream(fin,swap);
        maps = read_int_istream(fin,swap);

        fin.ignore(4*4);//DMIN,DMAX,DMEAN,ISPG 76,80,84,88

        nsymbt = read_int_istream(fin,swap);//92

        fin.ignore(4*25);//EXTRA //96-195
        origx = read_float_istream(fin,swap);//196
        origy = read_float_istream(fin,swap);//200
        origz = read_float_istream(fin,swap);//204
        fin.ignore(4*1);//MAP//208

*/

bool ToCubic(MRC *mrc){
 float *map;
 
 free(map);
 printf("#voxels= %d\n",mrc->NumVoxels);
 return false;
 return false;
}

//Just x2x2x2
bool upsampling(MRC *mrc,double t){
 int i,j,k,ind1,ind2,ind3;
 int Nth=omp_get_max_threads();
 double rgstep=1.000/mrc->widthx;
 float *map;
 unsigned int cnt=0;
 //Filtering Only
 if(mrc->widthx <= 1.5){
  for(i=0;i<mrc->NumVoxels;i++){
   if(mrc->dens[i]<t)
    mrc->dens[i]=0.00;
   else
    cnt++;
  }
  mrc->Nact=cnt;
  return false;
 }
 printf("#Start Upsampling Grid size= %f -> %f by %d threads\n",mrc->widthx,mrc->widthx*0.5,Nth);
 

 if((map=((float*)malloc(sizeof(float)*mrc->NumVoxels*8)))==NULL){
  free(map);
  return true;
 }
 int xydim1=mrc->xdim*mrc->ydim;//Original
 int xydim2=mrc->xdim*mrc->ydim*4;//New dim
 int xdim1=mrc->xdim;
 int ydim1=mrc->ydim;
 int zdim1=mrc->zdim;
 int xdim2=mrc->xdim*2;
 int ydim2=mrc->ydim*2;
 int zdim2=mrc->zdim*2;

 //copy
 for(int x=0;x<mrc->xdim;x++){
 for(int y=0;y<mrc->ydim;y++){
 for(int z=0;z<mrc->zdim;z++){
  ind1=xydim1*z+xdim1*y+x;
  ind2=xydim2*z*2+xdim2*y*2+x*2;
  map[ind2]=mrc->dens[ind1];
 }}}


 //simple & fast
 //#pragma omp parallel for schedule(dynamic,5)
 for(int x=0;x<xdim2;x+=2){
  for(int y=0;y<ydim2;y+=2){
   for(int z=1;z<zdim2-1;z+=2){
    ind1=xydim2*z+xdim2*y+x;
    ind2=xydim2*(z-1)+xdim2*y+x;
    ind3=xydim2*(z+1)+xdim2*y+x;
    map[ind1]=(map[ind2]+map[ind3])*0.5;
 }}}
 for(int x=0;x<xdim2;x+=2){
  for(int y=1;y<ydim2-1;y+=2){
   for(int z=0;z<zdim2;z+=2){
    ind1=xydim2*z+xdim2*y+x;
    ind2=xydim2*z+xdim2*(y+1)+x;
    ind3=xydim2*z+xdim2*(y-1)+x;
    map[ind1]=(map[ind2]+map[ind3])*0.5;
 }}}
 for(int x=1;x<xdim2-1;x+=2){
  for(int y=0;y<ydim2;y+=2){
   for(int z=0;z<zdim2;z+=2){
    ind1=xydim2*z+xdim2*y+x;
    ind2=xydim2*z+xdim2*y+(x+1);
    ind3=xydim2*z+xdim2*y+(x-1);
    map[ind1]=(map[ind2]+map[ind3])*0.5;
 }}}
 
 for(int x=1;x<xdim2-1;x+=2){
  for(int y=1;y<ydim2-1;y+=2){
   for(int z=0;z<zdim2;z+=2){
    ind1=xydim2*z+xdim2*y+x;
    ind2=xydim2*z+xdim2*y+(x+1);
    ind3=xydim2*z+xdim2*y+(x-1);
    map[ind1]=(map[ind2]+map[ind3])*0.5;
 }}}
 for(int x=1;x<xdim2-1;x+=2){
  for(int y=0;y<ydim2;y+=2){
   for(int z=1;z<zdim2-1;z+=2){
    ind1=xydim2*z+xdim2*y+x;
    ind2=xydim2*z+xdim2*y+(x+1);
    ind3=xydim2*z+xdim2*y+(x-1);
    map[ind1]=(map[ind2]+map[ind3])*0.5;
 }}}
 for(int x=0;x<xdim2;x+=2){
  for(int y=1;y<ydim2-1;y+=2){
   for(int z=1;z<zdim2-1;z+=2){
    ind1=xydim2*z+xdim2*y+x;
    ind2=xydim2*z+xdim2*(y+1)+x;
    ind3=xydim2*z+xdim2*(y-1)+x;
    map[ind1]=(map[ind2]+map[ind3])*0.5;
 }}}
 for(int x=1;x<xdim2-1;x+=2){
  for(int y=1;y<ydim2-1;y+=2){
   for(int z=1;z<zdim2-1;z+=2){
    ind1=xydim2*z+xdim2*y+x;
    ind2=xydim2*z+xdim2*y+(x+1);
    ind3=xydim2*z+xdim2*y+(x-1);
    map[ind1]=(map[ind2]+map[ind3])*0.5;
 }}}
 
 //update
 mrc->xdim*=2;
 mrc->ydim*=2;
 mrc->zdim*=2;
 mrc->mx*=2;
 mrc->my*=2;
 mrc->mz*=2;
 
 //mrc->dens=map;

 mrc->NumVoxels*=8;
 mrc->widthx*=0.5;
 mrc->widthy*=0.5;
 mrc->widthz*=0.5;

 free(mrc->dens);
 if((mrc->dens=((float*)malloc(sizeof(float)*mrc->NumVoxels)))==NULL){
  free(mrc->dens);
  return true;
 }
 //Filtering
 for(i=0;i<mrc->NumVoxels;i++){
  if(map[i]<t){
   mrc->dens[i]=0.00;
  }else{
   mrc->dens[i]=map[i];
   cnt++;
  }
 }
 mrc->Nact=cnt;

 free(map);
 printf("#voxels= %d\n",mrc->NumVoxels);
 return false;
}

void out_situs(MRC *mrc){
 int x,y,z,i,ind;
 int xydim=mrc->xdim*mrc->ydim;
 i=0;
 printf("%.6f ",mrc->widthx);
 printf("%.6f ",mrc->orgxyz[0]);
 printf("%.6f ",mrc->orgxyz[1]);
 printf("%.6f ",mrc->orgxyz[2]);
 printf("%d %d %d\n\n",mrc->xdim,mrc->ydim,mrc->zdim);
 for(z=0;z<mrc->zdim;z++){
  for(y=0;y<mrc->ydim;y++){
   for(x=0;x<mrc->xdim;x++){
    ind=xydim*z+mrc->xdim*y+x;
    printf("%11.6f ",mrc->dens[ind]);
    if(i>0 && i%10==9)
    printf("\n");
    i++;
   }
  }
 }
}


bool fastLDP(MRC *m,POINTS *p){
 int i,j,k,ind;
 int cnt=0;
 int xydim=m->xdim*m->ydim;

 //malloc
 if((p->cd=(double **)malloc(sizeof(double *)*m->NumVoxels))==NULL)
  return true;
 if((p->origrid=(int **)malloc(sizeof(int *)*m->NumVoxels))==NULL)
  return true;

 for(int x=0;x<m->xdim;x++){
 for(int y=0;y<m->ydim;y++){
 for(int z=0;z<m->zdim;z++){
  ind=xydim*z+m->xdim*y+x;
  if(m->dens[ind]==0.00)
   continue;
   if((p->cd[cnt]=(double *)malloc(sizeof(double)*3))==NULL)
    return true;
    //if((p->origrid[cnt]=(int *)malloc(sizeof(int)*3))==NULL)
    //return true;
   //map origin xyz and xwidth based coordinates

   p->cd[cnt][0]=(double)x;
   p->cd[cnt][1]=(double)y;
   p->cd[cnt][2]=(double)z;
/*
   //Original Grid positions
   p->origrid[cnt][0]=x;
   p->origrid[cnt][1]=y;
   p->origrid[cnt][2]=z;
*/
   //printf("cnd= %d\n",cnt);
   cnt++;
  
 }}}
 p->Ncd=cnt;
 p->Nori=cnt;

 puts("#Start LDP");
 if((p->dens=(double *)malloc(sizeof(double)*cnt))==NULL)
  return true;
 //Setup Filter
 //Gaussian kernel dreso=window size
 double dreso=cmd.dreso;
 double gstep=m->widthx;
 double fs=(dreso/gstep)*0.5;
 fs=fs*fs;
 double fsiv=1.000/fs;
 double fmaxd=(dreso/gstep)*2.0;
 printf("#maxd= %f\n",fmaxd);
 //Mean Shifting

 p->Ncd=cnt;
 #pragma omp parallel for schedule(dynamic,5)
 for(int i=0;i<cnt;i++){
  int stp[3],endp[3],ind2;
  double pos[3],pos2[3],ori[3];
  double tmpcd[3];
  double rx,ry,rz,d2;
  double v,dtotal,rd;
  pos[0]=p->cd[i][0];
  pos[1]=p->cd[i][1];
  pos[2]=p->cd[i][2];
  ori[0]=pos[0];
  ori[1]=pos[1];
  ori[2]=pos[2];
  while(1){
   //Start Point
   stp[0]=(int)(pos[0]-fmaxd);
   stp[1]=(int)(pos[1]-fmaxd);
   stp[2]=(int)(pos[2]-fmaxd);

   if(stp[0]<0)stp[0]=0;
   if(stp[1]<0)stp[1]=0;
   if(stp[2]<0)stp[2]=0;


   endp[0]=(int)(pos[0]+fmaxd+1);
   endp[1]=(int)(pos[1]+fmaxd+1);
   endp[2]=(int)(pos[2]+fmaxd+1);

   if(endp[0]>=m->xdim) endp[0]=m->xdim;
   if(endp[1]>=m->ydim) endp[1]=m->ydim;
   if(endp[2]>=m->zdim) endp[2]=m->zdim;

   dtotal=0;
   pos2[0]=pos2[1]=pos2[2]=0;
   for(int xp=stp[0];xp<endp[0];xp++){
    rx=(double)xp-pos[0];
    rx=rx*rx;
   for(int yp=stp[1];yp<endp[1];yp++){
    ry=(double)yp-pos[1];
    ry=ry*ry;
   for(int zp=stp[2];zp<endp[2];zp++){
    rz=(double)zp-pos[2];
    rz=rz*rz;
    d2=rx+ry+rz;
    //d=exp(-1.50*fr*fsiv)*amap(ii,jj,kk)
    ind2=xydim*zp+m->xdim*yp+xp;
    v=exp(-1.50*d2*fsiv)*m->dens[ind2];
    dtotal+=v;
    if(v>0)
    //printf("d %f %d %d %d\n",v,xp,yp,zp);
    pos2[0]+=v*(double)xp;
    pos2[1]+=v*(double)yp;
    pos2[2]+=v*(double)zp;
   }}}
   //printf("dto= %f\n",dtotal);
   if(dtotal==0.00)
    break;
   rd=1.00/dtotal;
   pos2[0]*=rd;
   pos2[1]*=rd;
   pos2[2]*=rd;
   tmpcd[0]=pos[0]-pos2[0];
   tmpcd[1]=pos[1]-pos2[1];
   tmpcd[2]=pos[2]-pos2[2];

   //move to other grid points
   if(fabs(pos2[0]-ori[0])>1.00||fabs(pos2[1]-ori[1])>1.00||fabs(pos2[2]-ori[2])>1.00){
    dtotal=0.00;
    break;
   }

   //Update
   pos[0]=pos2[0];
   pos[1]=pos2[1];
   pos[2]=pos2[2];
   //printf("*%d %f %f %f v= %f\n",i,pos[0],pos[1],pos[2],dtotal);
   if(tmpcd[0]*tmpcd[0]+tmpcd[1]*tmpcd[1]+tmpcd[2]*tmpcd[2]<0.001)
    break;
  }
  //printf("*%d %f %f %f\n",i,pos[0],pos[1],pos[2]);
  p->cd[i][0]=pos[0];
  p->cd[i][1]=pos[1];
  p->cd[i][2]=pos[2];
  p->dens[i]=dtotal;
 }

 puts("#End LDP");
 return false;
}






bool meanshift(MRC *m,POINTS *p){
 int i,j,k,ind;
 int cnt=0;
 int xydim=m->xdim*m->ydim;

 //malloc
 if((p->cd=(double **)malloc(sizeof(double *)*m->NumVoxels))==NULL)
  return true;
 if((p->origrid=(int **)malloc(sizeof(int *)*m->NumVoxels))==NULL)
  return true;
/*
 for(i=0;i<m->Nact;i++)
  if((p->cd[i]=(double *)malloc(sizeof(double)*3))==NULL)
  return true;
*/

 for(int x=0;x<m->xdim;x++){
 for(int y=0;y<m->ydim;y++){
 for(int z=0;z<m->zdim;z++){
  ind=xydim*z+m->xdim*y+x;
  if(m->dens[ind]==0.00)
   continue;
   if((p->cd[cnt]=(double *)malloc(sizeof(double)*3))==NULL)
    return true;
    if((p->origrid[cnt]=(int *)malloc(sizeof(int)*3))==NULL)
    return true;
   //map origin xyz and xwidth based coordinates
   p->cd[cnt][0]=(double)x;
   p->cd[cnt][1]=(double)y;
   p->cd[cnt][2]=(double)z;

   //Original Grid positions
   p->origrid[cnt][0]=x;
   p->origrid[cnt][1]=y;
   p->origrid[cnt][2]=z;

   //printf("cnd= %d\n",cnt);
   cnt++;
  
 }}}
 p->Ncd=cnt;
 p->Nori=cnt;
 puts("#Start MS");
 if((p->dens=(double *)malloc(sizeof(double)*cnt))==NULL)
  return true;
 //Setup Filter
 //Gaussian kernel dreso=window size
 double dreso=cmd.dreso;
 double gstep=m->widthx;
 double fs=(dreso/gstep)*0.5;
 fs=fs*fs;
 double fsiv=1.000/fs;
 double fmaxd=(dreso/gstep)*2.0;
 printf("#maxd= %f\n",fmaxd);
 //Mean Shifting

 p->Ncd=cnt;
 #pragma omp parallel for schedule(dynamic,5)
 for(int i=0;i<cnt;i++){
  int stp[3],endp[3],ind2;
  double pos[3],pos2[3];
  double tmpcd[3];
  double rx,ry,rz,d2;
  double v,dtotal,rd;
  pos[0]=p->cd[i][0];
  pos[1]=p->cd[i][1];
  pos[2]=p->cd[i][2];
  while(1){
   //Start Point
   stp[0]=(int)(pos[0]-fmaxd);
   stp[1]=(int)(pos[1]-fmaxd);
   stp[2]=(int)(pos[2]-fmaxd);

   if(stp[0]<0)stp[0]=0;
   if(stp[1]<0)stp[1]=0;
   if(stp[2]<0)stp[2]=0;


   endp[0]=(int)(pos[0]+fmaxd+1);
   endp[1]=(int)(pos[1]+fmaxd+1);
   endp[2]=(int)(pos[2]+fmaxd+1);

   if(endp[0]>=m->xdim) endp[0]=m->xdim;
   if(endp[1]>=m->ydim) endp[1]=m->ydim;
   if(endp[2]>=m->zdim) endp[2]=m->zdim;

   dtotal=0;
   pos2[0]=pos2[1]=pos2[2]=0;
   for(int xp=stp[0];xp<endp[0];xp++){
    rx=(double)xp-pos[0];
    rx=rx*rx;
   for(int yp=stp[1];yp<endp[1];yp++){
    ry=(double)yp-pos[1];
    ry=ry*ry;
   for(int zp=stp[2];zp<endp[2];zp++){
    rz=(double)zp-pos[2];
    rz=rz*rz;
    d2=rx+ry+rz;
    //d=exp(-1.50*fr*fsiv)*amap(ii,jj,kk)
    ind2=xydim*zp+m->xdim*yp+xp;
    v=exp(-1.50*d2*fsiv)*m->dens[ind2];
    dtotal+=v;
    if(v>0)
    //printf("d %f %d %d %d\n",v,xp,yp,zp);
    pos2[0]+=v*(double)xp;
    pos2[1]+=v*(double)yp;
    pos2[2]+=v*(double)zp;
   }}}
   //printf("dto= %f\n",dtotal);
   if(dtotal==0.00)
    break;
   rd=1.00/dtotal;
   pos2[0]*=rd;
   pos2[1]*=rd;
   pos2[2]*=rd;
   tmpcd[0]=pos[0]-pos2[0];
   tmpcd[1]=pos[1]-pos2[1];
   tmpcd[2]=pos[2]-pos2[2];

   pos[0]=pos2[0];
   pos[1]=pos2[1];
   pos[2]=pos2[2];
   //printf("*%d %f %f %f v= %f\n",i,pos[0],pos[1],pos[2],dtotal);
   if(tmpcd[0]*tmpcd[0]+tmpcd[1]*tmpcd[1]+tmpcd[2]*tmpcd[2]<0.001)
    break;
  }
  //printf("*%d %f %f %f\n",i,pos[0],pos[1],pos[2]);
  p->cd[i][0]=pos[0];
  p->cd[i][1]=pos[1];
  p->cd[i][2]=pos[2];
  p->dens[i]=dtotal;
 }

 puts("#End MS");
 return false;
}

double meanshift_pos(MRC *m,double pos[3]){
 int i,j,k,ind;
 int cnt=0;
 int xydim=m->xdim*m->ydim;
 POINTS *p;

 //Setup Filter
 //Gaussian kernel dreso=window size
 double dreso=cmd.dreso;
 double gstep=m->widthx;
 double fs=(dreso/gstep)*0.5;
 fs=fs*fs;
 double fsiv=1.000/fs;
 double fmaxd=(dreso/gstep)*2.0;
 
 //Mean Shifting
  int stp[3],endp[3],ind2;
  double tmpcd[3];
  double rx,ry,rz,d2;
  double v,dtotal,rd;

   //Start Point
   stp[0]=(int)(pos[0]-fmaxd);
   stp[1]=(int)(pos[1]-fmaxd);
   stp[2]=(int)(pos[2]-fmaxd);

   if(stp[0]<0)stp[0]=0;
   if(stp[1]<0)stp[1]=0;
   if(stp[2]<0)stp[2]=0;


   endp[0]=(int)(pos[0]+fmaxd+1);
   endp[1]=(int)(pos[1]+fmaxd+1);
   endp[2]=(int)(pos[2]+fmaxd+1);

   if(endp[0]>=m->xdim) endp[0]=m->xdim;
   if(endp[1]>=m->ydim) endp[1]=m->ydim;
   if(endp[2]>=m->zdim) endp[2]=m->zdim;

   dtotal=0;
   for(int xp=stp[0];xp<endp[0];xp++){
    rx=(double)xp-pos[0];
    rx=rx*rx;
   for(int yp=stp[1];yp<endp[1];yp++){
    ry=(double)yp-pos[1];
    ry=ry*ry;
   for(int zp=stp[2];zp<endp[2];zp++){
    rz=(double)zp-pos[2];
    rz=rz*rz;
    d2=rx+ry+rz;
    //d=exp(-1.50*fr*fsiv)*amap(ii,jj,kk)
    ind2=xydim*zp+m->xdim*yp+xp;
    v=exp(-1.50*d2*fsiv)*m->dens[ind2];
    dtotal+=v;
   }}}


 return dtotal;
}


bool MergePoints(MRC *m,POINTS *p){
 double dcut=cmd.MergeDist/m->widthx;
 double d2cut=dcut*dcut;
 double rdcut=cmd.Filter;
 double dmax,dmin,drange;
 bool *stock;
 int *tmp_member,*tmp_member2;

 if((stock=(bool *)malloc(sizeof(bool)*p->Ncd))==NULL)
  return true;
 if((tmp_member=(int *)malloc(sizeof(int)*p->Ncd))==NULL)
  return true;
 if((tmp_member2=(int *)malloc(sizeof(int)*p->Ncd))==NULL)
  return true;
 dmax=0;dmin=999999.99;

 for(int i=0;i<p->Ncd;i++){
  if(p->dens[i]<dmin)
   dmin=p->dens[i];
  if(p->dens[i]>dmax)
   dmax=p->dens[i];
 }
 drange=dmax-dmin;
 double rv_range=1.00/drange;
 printf("#dmax= %f dmin= %f\n",dmax,dmin);
 //init member
 if((p->member=(int *)malloc(sizeof(int)*p->Ncd))==NULL)
  return true;
 for(int i=0;i<p->Ncd;i++)
  p->member[i]=i;
 for(int i=0;i<p->Ncd;i++)
  stock[i]=true;

 //#pragma omp parallel for schedule(dynamic,5)
 for(int i=0;i<p->Ncd-1;i++){
  double tmp[3],d2;
  if((p->dens[i]-dmin)*rv_range < rdcut)
   stock[i]=false;


  if(stock[i]==false)
   continue;
  for(int j=i+1;j<p->Ncd;j++){
   if(stock[j]==false)
    continue;
   tmp[0]=p->cd[i][0]-p->cd[j][0];
   tmp[1]=p->cd[i][1]-p->cd[j][1];
   tmp[2]=p->cd[i][2]-p->cd[j][2];
   d2=tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2];


   if(d2<d2cut){
    //Keep high dens
    if(p->dens[i]>p->dens[j]){
     stock[j]=false;
     p->member[j]=i;
    }else{
     stock[i]=false;
     p->member[i]=j;
     break;
    }
   }
  }
 }
 //printf("62809; %d\n",p->member[62809]);
 //Update member data
 for(int i=0;i<p->Ncd;i++){
  int now=p->member[i];
  for(int j=0;j<p->Ncd;j++){
   if(now==p->member[now])
    break;
   now=p->member[now];
  }
  p->member[i]=now;
 }
 //printf("62809; %d\n",p->member[62809]);
 //Copy
 int Nmerge=0;
 for(int i=0;i<p->Ncd;i++){
  if(stock[i]){
   p->cd[Nmerge][0]=p->cd[i][0];
   p->cd[Nmerge][1]=p->cd[i][1];
   p->cd[Nmerge][2]=p->cd[i][2];
   p->dens[Nmerge]=p->dens[i];
   tmp_member[i]=Nmerge;
   Nmerge++;
  }else{
   tmp_member[i]=-1;
  }
 }
 for(int i=0;i<p->Ncd;i++)
  tmp_member2[i]=tmp_member[p->member[i]];
 for(int i=0;i<p->Ncd;i++)
  p->member[i]=tmp_member2[i];
 
 //printf("62809 tmp; %d\n",tmp_member[62809]);
 //printf("62809; %d\n",p->member[62809]);
 printf("#After Merge: %d\n",Nmerge);
 p->Ncd=Nmerge;
 return false;
}


void ShowModel(MRC *m,POINTS *p){
/*
write(*,'("MODEL    3")')
	do kk=1,SpaseN
	 ftmp_xyz=(after(:,kk)-1)*gstep+gbase;
	 write(*,'("ATOM  ",I5,"  CA  ALA A",I4,"    ",3f8.3,f6.2,f6.2)'),
     *   Natm,Natm,ftmp_xyz,1.00,fin(kk)/fmax
	 Natm=Natm+1
	enddo
*/


 int i,j,k;
 int Natm=1;
 double tmp[3];
 printf("MODEL\n");
 for(i=0;i<p->Ncd;i++){
  tmp[0]=p->cd[i][0]*m->widthx+m->orgxyz[0];
  tmp[1]=p->cd[i][1]*m->widthx+m->orgxyz[1];
  tmp[2]=p->cd[i][2]*m->widthx+m->orgxyz[2];
  printf("ATOM  %5d  CA  ALA%6d    ",Natm,Natm);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,p->dens[i]);
  //printf("%d %d %d %f\n",p->origrid[i][0],p->origrid[i][1],p->origrid[i][2],p->dens[i]);
  Natm++;
 }
 printf("TER\nEND\n");
}

void ShowLDP(MRC *m,POINTS *p){
 int i,j,k;
 int Natm=1;
 double tmp[3];
 double dmax=0;

 for(i=0;i<p->Ncd;i++)
  if(p->dens[i]>dmax)
   dmax=p->dens[i];

 printf("MODEL\n");
 for(i=0;i<p->Ncd;i++){
  if(p->dens[i]==0.00) continue;
  tmp[0]=p->cd[i][0]*m->widthx+m->orgxyz[0];
  tmp[1]=p->cd[i][1]*m->widthx+m->orgxyz[1];
  tmp[2]=p->cd[i][2]*m->widthx+m->orgxyz[2];
  printf("ATOM  %5d  CA  ALA%6d    ",Natm,Natm);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,p->dens[i]/dmax);
  Natm++;
 }
 printf("TER\nEND\n");
}



void ShowOri(MRC *m,POINTS *p){

 int i,j,k;
 int Natm=1;
 double tmp[3];
 printf("MODEL\n");
 for(i=0;i<p->Nori;i++){
  if(p->mask[i]==0.0)
   continue;
  tmp[0]=(float)p->origrid[i][0]*m->widthx+m->orgxyz[0];
  tmp[1]=(float)p->origrid[i][1]*m->widthx+m->orgxyz[1];
  tmp[2]=(float)p->origrid[i][2]*m->widthx+m->orgxyz[2];
  printf("ATOM  %5d  CA  ALA%6d    ",Natm,Natm);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,p->mask[i]);
  Natm++;
 }
 printf("TER\nEND\n");
}



void ShowGraph(GRAPH *g){
 int i,j;

 for(i=0;i<g->Ne;i++)
  if(g->edge[i].local || g->edge[i].mst)
   printf("BOND %d %d\n",g->edge[i].id1+1,g->edge[i].id2+1);

}

void ShowTree(GRAPH *g, TREE *t){
 int i,j;
 //puts("#Show MSTree");
 for(i=0;i<g->Ne;i++){
  if(t->ActE[i])
  printf("BOND %d %d\n",g->edge[i].id1+1,g->edge[i].id2+1);
 }
}

void ShowPath(MRC *m,POINTS *p,GRAPH *g, TREE *t, int n){

 int i,j,k,now;
 int Natm=1;
 double tmp[3];

 for(int mid=0;mid<n;mid++){
  double score=QualityTree(g, &t[mid]);
  printf("#SCORE: %f\n",score);
  printf("MODEL %d\n",mid+1);
	Natm=1;
	now=t[mid].St;
	while(1){

	 //printf("%d -> %d\n",now,t[mid].nextv[now]);

	 i=now;
	 tmp[0]=p->cd[i][0]*m->widthx+m->orgxyz[0];
 	 tmp[1]=p->cd[i][1]*m->widthx+m->orgxyz[1];
 	 tmp[2]=p->cd[i][2]*m->widthx+m->orgxyz[2];
 	 printf("ATOM  %5d  CA  ALA%6d    ",Natm,Natm);
 	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,p->dens[i]);
 	 Natm++;



	 now=t[mid].nextv[now];
	 if(now==-1)
	  break;
	}

/*
 	printf("MODEL %2d\n",mid+1);
 	for(i=0;i<p->Ncd;i++){
 	 tmp[0]=p->cd[i][0]*m->widthx+m->orgxyz[0];
 	 tmp[1]=p->cd[i][1]*m->widthx+m->orgxyz[1];
 	 tmp[2]=p->cd[i][2]*m->widthx+m->orgxyz[2];
 	 printf("ATOM  %5d  CA  ALA%6d    ",Natm,Natm);
 	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,p->dens[i]);
 	 Natm++;
 	}
*/
  printf("TER\nENDMDL\n");
 }
}

void NR_subopt(SUBOPT_ALI *sub){
 int N=0;
	for(N=0;N<2000;N++)
	 if(sub[N].L==0)
	  break;
 for(int i=0;i<N;i++){
	  if(sub[i].L==-1)
		  continue;
  for(int j=i+1;j<N;j++){
	  if(sub[j].L==-1)
		  continue;
	//compare
	bool flag1=false;
	bool flag2=false;
	for(int k1=0;k1<sub[i].L;k1++){
	 flag1=false;
	 for(int k2=0;k2<sub[j].L;k2++){
		//seq
		if(sub[i].seq[k1]==sub[j].seq[k2] && sub[i].model[k1]==sub[j].model[k2]){
		 flag1=true;
		 continue;
		}
	 }
	 if(flag1==false)
	  break;
	}
	//same or smaller
	if(flag1==true){
	 sub[i].L=-1;
	}
	for(int k2=0;k2<sub[j].L;k2++){
	 flag2=false;
	 for(int k1=0;k1<sub[i].L;k1++){
		//seq
		if(sub[i].seq[k1]==sub[j].seq[k2] && sub[i].model[k1]==sub[j].model[k2]){
		 flag2=true;
		 continue;
		}
	 }
	 if(flag2==false)
	  break;
	}
	//same or smaller
	if(flag2==true){
	 sub[j].L=-1;
	}
	
  }
 }
}


void ShowPath2(MRC *m,POINTS *p,GRAPH *g, TREE *t, int n, SEQFG *sfg){
 //Show path
 //For density*volume, using K-clustring data
 int i,j,k,now,ind;
 int Natm=1;
 double tmp[3];
 int st,ed;
 double fmax=0;
 int xydim=m->xdim*m->ydim;
 SUBOPT_ALI *subopt;//fragment output
 char ch_array[61]="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz123456789";
 subopt = (SUBOPT_ALI *)malloc(sizeof(SUBOPT_ALI)*1000);
 DP_MEMORY dpm;

 dpm.ali1=(int *)malloc(sizeof(int)*g->Nnode);
 dpm.ali2=(int *)malloc(sizeof(int)*sfg->l);
 dpm.gali=(int *)malloc(sizeof(int)*sfg->l*g->Nnode);
 dpm.Smtx=(float *)malloc(sizeof(float)*sfg->l*g->Nnode);
 dpm.SmtxRv=(float *)malloc(sizeof(float)*sfg->l*g->Nnode);
 dpm.Smtx_ori=(float *)malloc(sizeof(float)*sfg->l*g->Nnode);
 dpm.SmtxRv_ori=(float *)malloc(sizeof(float)*sfg->l*g->Nnode);
 dpm.dtbl=(float *)malloc(sizeof(float)*g->Nnode*21);
 dpm.dtblR=(float *)malloc(sizeof(float)*g->Nnode*21);
 dpm.Dmtx=(DPMTX *)malloc(sizeof(DPMTX)*sfg->l*g->Nnode);
 dpm.model=(NODE **)malloc(sizeof(NODE*)*sfg->l);
 dpm.path = (NODE **)malloc(sizeof(NODE*)*g->Nnode);

 printf("##VER %.3f\n",VER);
 //path data
 //model ID
 for(int mid=0;mid<n;mid++){
  st=t[mid].St;
  ed=t[mid].Ed;
  fmax=0;
	//float score=QualityTreeDP(g,&t[mid],sfg,&dpm,true);
	float score=QualityTreeDPSubOpt(g,&t[mid],sfg,&dpm,subopt);
	NR_subopt(subopt);
	//Show CA models
	Natm=1;
  	printf("#Tree SCORE: %f %f\n",t[mid].score,score);
	for(int fr=0;fr<2000;fr++){
	 if(subopt[fr].L==-1)
	  continue;
	 if(subopt[fr].L==0)
	  break;
	 printf("#FRAG MID=%d L=%d SCO= %f\n",mid,subopt[fr].L,subopt[fr].sco);
		for(int rnum=0;rnum<subopt[fr].L;rnum++){
		 tmp[0]=subopt[fr].model[rnum]->real_cd[0];
		 tmp[1]=subopt[fr].model[rnum]->real_cd[1];
		 tmp[2]=subopt[fr].model[rnum]->real_cd[2];
 		 printf("ATOM  %5d  CA  %3s %c%4d    ",Natm
		 ,RES_NAMES[sfg->seq[subopt[fr].seq[rnum]]]
		 ,ch_array[sfg->CID[subopt[fr].seq[rnum]]],subopt[fr].seq[rnum]+1);
		 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,1.0);
		 Natm++;
		}
	 printf("TER\n");
	}
	
/*
  	//printf("#SCORE: %f LEN=%d fmax=%f\n",t[mid].score,t[mid].Lpath,fmax);
  	printf("#SCORE: %f\n",t[mid].score);
	printf("##SCORE= %f\n",score);
  	printf("MODEL %d\n",mid+1);
 	Natm=1;
	for(int rnum=0;rnum<sfg->l;rnum++){
		if(dpm.ali2[rnum]==-1)
			continue;
	 //printf("poi= %d\n",dpm.ali2[rnum]);
	 tmp[0]=dpm.model[rnum]->real_cd[0];
	 tmp[1]=dpm.model[rnum]->real_cd[1];
	 tmp[2]=dpm.model[rnum]->real_cd[2];
 	 printf("ATOM  %5d  CA  %3s%6d    ",Natm
			 ,RES_NAMES[sfg->seq[rnum]]
			 ,rnum+1);
 	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,1.0);
 	 Natm++;
	 if(dpm.ali2[rnum+1]==-1){
	  printf("TER\n");
	 }
	}
  	printf("ENDMDL\n");
*/
 }
}



int cmp_edge_d(const void *c1, const void *c2){

 EDGE a=*(EDGE *)c1;
 EDGE b=*(EDGE *)c2;

 if(a.d<b.d) return -1;
 if(a.d>b.d) return 1;
 return 0;
}

int cmp_tree_score(const void *c1, const void *c2){

 TREE a=*(TREE *)c1;
 TREE b=*(TREE *)c2;

 if(a.score>b.score) return -1;
 if(a.score<b.score) return 1;
 return 0;
}


//Make Graph and MST
bool SetUpGraph(POINTS *p, GRAPH *g,MRC *mrc,TREE *mst){
 int i,j,dif;
 double dcut=cmd.LocalR/mrc->widthx;
 double d2cut=dcut*dcut;
 int *cid;
 int Nmin_ldp = 20;//Minimum LDPs in the fragment

 //MALLOC GRAPH------
 if((g->adj=(bool **)malloc(sizeof(bool *)*p->Ncd))==NULL)
  return true;
 for(i=0;i<p->Ncd;i++){
  if((g->adj[i]=(bool *)malloc(sizeof(bool)*p->Ncd))==NULL)
   return true;
  for(j=0;j<p->Ncd;j++)
   g->adj[i][j]=false;
 }
 if((g->cid=(int *)malloc(sizeof(int)*p->Ncd))==NULL)
  return true;
 if((cid=(int *)malloc(sizeof(int)*p->Ncd))==NULL)
  return true;

 if((g->node=(NODE *)malloc(sizeof(NODE)*p->Ncd))==NULL)
  return true;

 puts("#Fin Malloc");
 //END---------------


 g->Nnode=p->Ncd;
 int Ne=0;
 #pragma omp parallel for reduction(+:Ne) schedule(dynamic,5)
 for(int ii=0;ii<p->Nori;ii++){
  int m1=p->member[ii];
  int m2;

  if(m1==-1)
   continue;

  for(int jj=ii+1;jj<p->Nori;jj++){
   m2=p->member[jj];
   if(m2==-1)
    continue;
   if(m1==m2)
    continue;
   if(g->adj[m1][m2])
    continue;
   //check
   if((p->origrid[ii][0]-p->origrid[jj][0])*(p->origrid[ii][0]-p->origrid[jj][0])>1)
    continue;
   if((p->origrid[ii][1]-p->origrid[jj][1])*(p->origrid[ii][1]-p->origrid[jj][1])>1)
    continue;
   if((p->origrid[ii][2]-p->origrid[jj][2])*(p->origrid[ii][2]-p->origrid[jj][2])>1)
    continue;
   
   g->adj[m1][m2]=true;
   g->adj[m2][m1]=true;
   Ne++;
   //printf("%d m1:%d %d m2:%d\n",ii,m1,jj,m2);
  }
 }
 printf("#Fin checking connect Ne= %d\n",Ne);
 if((g->edge=(EDGE *)malloc(sizeof(EDGE)*Ne))==NULL)
  return true;


 Ne=0;
 for(i=0;i<g->Nnode;i++){
  double d=0;
  double tmp[3];
  g->cid[i]=i;
  for(j=i+1;j<g->Nnode;j++){
   if(g->adj[i][j]==false)
    continue;
   tmp[0]=p->cd[i][0]-p->cd[j][0];
   tmp[1]=p->cd[i][1]-p->cd[j][1];
   tmp[2]=p->cd[i][2]-p->cd[j][2];

   d=sqrt(tmp[0]*tmp[0]+tmp[1]*tmp[1]+tmp[2]*tmp[2]);
   g->edge[Ne].d=d;//map based
   g->edge[Ne].id1=i;
   g->edge[Ne].id2=j;
   //printf("*%d %f %d %d\n",Ne,g->edge[Ne].d,g->edge[Ne].id1,g->edge[Ne].id2);
   Ne++;
  }
 }

 //sort
 qsort(g->edge,Ne,sizeof(EDGE),cmp_edge_d);

 //MST
 int v1,v2,tmp_cid;
 int Nt=0;
 EDGE **tree;
 if((tree=(EDGE **)malloc(sizeof(EDGE *)*Ne))==NULL)
  return true;

 int MaxCid=0;
 for(i=0;i<Ne;i++){
  v1=g->edge[i].id1;
  v2=g->edge[i].id2;

  if(g->cid[v1]==g->cid[v2])
   continue;

  tree[Nt]=&(g->edge[i]);
  tmp_cid=g->cid[v2];
  Nt++;
  if(MaxCid<tmp_cid)
   MaxCid=tmp_cid;
  //update cid in the tree
  for(j=0;j<Nt;j++){
   if(g->cid[tree[j]->id1]==tmp_cid)
    g->cid[tree[j]->id1]=g->cid[v1];
   if(g->cid[tree[j]->id2]==tmp_cid)
    g->cid[tree[j]->id2]=g->cid[v1];
  }
 }
 g->tree=tree;
 g->Nt=Nt;
 printf("#Nt= %d\n",Nt);
 printf("#MaxCid= %d\n",MaxCid);
 int *Ncid;
 if((Ncid=(int *)calloc(sizeof(int),(MaxCid+1)))==NULL)
  return true;
 for(i=0;i<g->Nnode;i++)
  if(g->cid[i]<=MaxCid)//count Ncid
  Ncid[g->cid[i]]++;
 int UseCid=-1;
 int Nuse=0;
 for(i=0;i<=MaxCid;i++){
  if(Nuse<Ncid[i]){
   UseCid=i;
   Nuse=Ncid[i];
  }
 }
 //printf("#UseCid= %d N= %d/%d\n",UseCid,Nuse,g->Nnode);
 if(UseCid==-1)
  return true;
 
 //Show Nldps by chain ID
 for(int i=0;i<=MaxCid;i++)
  	if(Ncid[i]>0)
 	 printf("#CID %d N= %d\n",i,Ncid[i]);

 //clean edge shift
 //Remove small isolated tree
 int Ntmp=0;

 for(i=0;i<Ne;i++){
  //if(UseCid==g->cid[g->edge[i].id1]){
  if(Ncid[g->cid[g->edge[i].id1]]>=Nmin_ldp){
   //printf("#Use cid=%d\n",g->cid[g->edge[i].id1]);
   g->edge[Ntmp]=g->edge[i];
   Ntmp++;
  }
 }
 g->Ne=Ntmp;

 //sort
 qsort(g->edge,g->Ne,sizeof(EDGE),cmp_edge_d);
 
 //
 Nt=0;
 //init
 for(i=0;i<g->Nnode;i++)
  g->cid[i]=i;
 for(i=0;i<g->Ne;i++){
  v1=g->edge[i].id1;
  v2=g->edge[i].id2;
  g->edge[i].mst=false;
  g->edge[i].local=false;
  //assigin edge-id
  g->edge[i].eid=i;//!!!!!!
  if(g->cid[v1]==g->cid[v2])
   continue;

  tree[Nt]=&(g->edge[i]);
  tmp_cid=g->cid[v2];
  g->edge[i].mst=true;//Used in MST
  Nt++;
  if(MaxCid<tmp_cid)
   MaxCid=tmp_cid;
  //update cid in the tree
  for(j=0;j<Nt;j++){
   if(g->cid[tree[j]->id1]==tmp_cid)
    g->cid[tree[j]->id1]=g->cid[v1];
   if(g->cid[tree[j]->id2]==tmp_cid)
    g->cid[tree[j]->id2]=g->cid[v1];
  }
 }
 //g->tree=tree;
 g->Nt=Nt;

 printf("#After cleaning.. Nt= %d Ne= %d\n",Nt,Ne);

 //return -1; 

 //Clean isolated dens data from mrc New!!! Apr 13th 2018
 if((p->mask=(float *)calloc(sizeof(float),p->Nori))==NULL)
  return 0;
 float *tmp_mask;
 if((tmp_mask=(float *)calloc(sizeof(float),p->Nori))==NULL)
  return 0;

 for(i=0;i<g->Ne;i++){
  if(g->edge[i].mst==false)
   continue;
  v1=g->edge[i].id1;
  v2=g->edge[i].id2;

  tmp_mask[v1]=1.00;//member
  tmp_mask[v2]=1.00;//member

 }
 for(int ii=0;ii<p->Nori;ii++){
  int m1=p->member[ii];
  p->mask[ii]=tmp_mask[m1];
 }


 //Set Edge dens
 puts("#Setting Edge dens...");
 #pragma omp parallel for schedule(dynamic,5)
 for(int ii=0;ii<g->Ne;ii++){
  double cd1[3],vec[3],dens,MinDens;
  int v1=g->edge[ii].id1;
  int v2=g->edge[ii].id2;
 
  //v1->v2
  vec[0]=p->cd[v2][0]-p->cd[v1][0];
  vec[1]=p->cd[v2][1]-p->cd[v1][1];
  vec[2]=p->cd[v2][2]-p->cd[v1][2];
  
  MinDens=99999999;
  for(int jj=1;jj<11;jj++){
   cd1[0]=p->cd[v1][0]+vec[0]*0.1*(double)(jj);
   cd1[1]=p->cd[v1][1]+vec[1]*0.1*(double)(jj);
   cd1[2]=p->cd[v1][2]+vec[2]*0.1*(double)(jj);
   dens=meanshift_pos(mrc,cd1);
   if(dens<MinDens)
    MinDens=dens;
  }
  //printf("#dens %d %f / %f v1:%d v2:%d\n",ii,dens,0.5*(p->dens[v1]+p->dens[v2]),v1,v2);
  g->edge[ii].dens=MinDens*g->edge[ii].d;
  //printf("#MinDens %d = %f\n",ii,g->edge[ii].dens);
 }

 puts("#Setting Local MSTs...");
 //Local MST
 #pragma omp parallel for schedule(dynamic,5)
 for(int ii=0;ii<g->Nnode;ii++){
  double vec[3],d2;
  int tmpid;
  unsigned int cid_th[100000];
  	//init cid
  	for(int jj=0;jj<g->Nnode;jj++)
  	 cid_th[jj]=jj;

  	for(int jj=0;jj<g->Ne;jj++){
  	 int v1=g->edge[jj].id1;
  	 int v2=g->edge[jj].id2;
  	 if(cid_th[v1]==cid_th[v2])
  	  continue;
	 //dist
	 vec[0]=p->cd[ii][0]-p->cd[v1][0];
 	 vec[1]=p->cd[ii][1]-p->cd[v1][1];
 	 vec[2]=p->cd[ii][2]-p->cd[v1][2];
	 d2=vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
	 if(d2>d2cut)
	  continue;

	 vec[0]=p->cd[ii][0]-p->cd[v2][0];
 	 vec[1]=p->cd[ii][1]-p->cd[v2][1];
 	 vec[2]=p->cd[ii][2]-p->cd[v2][2];
	 d2=vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2];
	 if(d2>d2cut)
	  continue;

	 g->edge[jj].local=true;

	 tmpid=cid_th[v2];
	 for(int kk=0;kk<g->Nnode;kk++){
	  if(cid_th[kk]==tmpid)//update cid
	   cid_th[kk]=cid_th[v1];
	 }
  	}
 }
 puts("#End Local MST");
 //Set Keep Flag
 double dkeep=cmd.Dkeep/mrc->widthx;
 for(int jj=0;jj<g->Ne;jj++){
  g->edge[jj].keep=false;

  if(g->edge[jj].mst && g->edge[jj].d<dkeep)
   g->edge[jj].keep=true;
 }

 //init
 for(int i=0;i<g->Nnode;i++)
  g->node[i].N=0;
 //input GRAPH
 for(int i=0;i<g->Ne;i++){
  if(g->edge[i].local==false && g->edge[i].mst==false)
   continue;
  int id=g->edge[i].id1;
  g->node[id].e[g->node[id].N]=&(g->edge[i]);
  g->node[id].N++;

  //printf("id %d : N= %d\n",id,g->node[id].N);
  id=g->edge[i].id2;
  g->node[id].e[g->node[id].N]=&(g->edge[i]);
  g->node[id].N++;
  //printf("id %d : N= %d\n",id,g->node[id].N);
 }

 //input MST
 puts("#Input MST");
 mst->Nnode=g->Nnode;
 mst->Ne=g->Nt;
 mst->Etotal=g->Ne;
 TREE tmp_tree;
 if(InitTree(mst,true))
  return true;

 for(int i=0;i<g->Ne;i++){
  if(g->edge[i].mst==false)
   continue;
  int id=g->edge[i].id1;
  int eid=g->edge[i].eid;
  //printf("MSTe= %d id1=%d id2=%d N= %d\n",i,g->edge[i].id1,g->edge[i].id2,mst->node[id].N);
  //mst->node[id].e[mst->node[id].N]=&(g->edge[i]);
  //mst->node[id].N++;

  //id=g->edge[i].id2;
  //mst->node[id].e[mst->node[id].N]=&(g->edge[i]);
  //mst->node[id].N++;

  mst->len+=g->edge[i].d;
  mst->St=id;//Starting Point

  mst->ActE[eid]=true;//Active Edge = MST
 }

 puts("#End SetUp");
 return false;
}


//Simple Tabu search
bool Tabu(GRAPH *g,TREE *init_t,TREE *Tr, SEQFG *sfg){
 int Nround=cmd.Nround;
 int Nnb=cmd.Nnb;
 int Ntabu=cmd.Ntabu;
 int Nsim=cmd.Nsim;
 double AllowLen=init_t->len*cmd.Allow;
 MOVE *move;
 int tabu[1000];
 int Ntb=0;
 int tbp=0;
 TREE Tnow,Tbest,Tnb[500];

 //Memory Space
 DP_MEMORY dpm;
 dpm.ali1=(int *)malloc(sizeof(int)*g->Nnode);
 dpm.ali2=(int *)malloc(sizeof(int)*sfg->l);
 dpm.gali=(int *)malloc(sizeof(int)*sfg->l*g->Nnode);
 dpm.Smtx=(float *)malloc(sizeof(float)*sfg->l*g->Nnode);
 dpm.SmtxRv=(float *)malloc(sizeof(float)*sfg->l*g->Nnode);
 dpm.dtbl=(float *)malloc(sizeof(float)*g->Nnode*21);//Distance Data
 dpm.dtblR=(float *)malloc(sizeof(float)*g->Nnode*21);//Distance Data
 dpm.Dmtx=(DPMTX *)malloc(sizeof(DPMTX)*sfg->l*g->Nnode);
 dpm.model=(NODE **)malloc(sizeof(NODE*)*sfg->l);
 dpm.path = (NODE **)malloc(sizeof(NODE*)*g->Nnode);


 move=(MOVE *)malloc(sizeof(MOVE)*Nnb);

 printf("#MST Len= %f\n",init_t->len);


 double InitSco;
 TREE t[100];
 for(int i=0;i<Nnb;i++){
  CopyTree(init_t,&Tnb[i],true);//copy& malloc
 }
 puts("#Fin Malloc for Tnb");
 for(int i=0;i<Nsim;i++)
  CopyTree(init_t,&Tr[i],true);
 puts("#Fin Malloc for Tr");
 CopyTree(init_t,&Tnow,true);
 puts("#Fin Malloc for Tnow");
 CopyTree(init_t,&Tbest,true);
 puts("#Fin Malloc for Tbest");

 double score=QualityTreeDP(g,&Tnow,sfg,&dpm,false);

 printf("#Initial score= %f\n",score);
 Tnow.score=score;
 Tbest.score=score;

 for(int i=0;i<Nsim;i++){
  CopyTree(init_t,&Tnow,false);
  //InitSco=QualityTree(g,&Tnow);
  //printf("#MST Score= %f %d\n",InitSco,i);

  //init tabu
  tbp=0;
  Ntb=0;
  double pre_best=0;

	//printf(" score: %f\n",move[0].score);
  	for(int n=0;n<Nround;n++){
   	 SetCutTbl(g,&Tnow,tabu,Ntb,AllowLen);
	 //puts("SetCut");
	 //if(Tnow.Nmv<Nnb){
	 if(Tnow.Nmv==0){
	  printf("#Reset tabu Nmv= %d\n",Tnow.Nmv);
	  //reset
	  Ntb=0;
	  tbp=0;
	  //n--;
	  continue;
	 }

	 //#pragma omp parallel for schedule(dynamic,5)
	 for(int nb=0;nb<Nnb;nb++){
 	  if(nb>=Tnow.Nmv){
	   move[nb].score=0;
	   Tnb[nb].score=0;
	   continue;
	  }
	 //for(int nb=0;nb<Nnb && nb < Tnow.Nmv;nb++){
	  //printf("#Sim: %3d/%3d Round= %4d/%4d nb=%d\n",i,Nsim,n+1,Nround,nb);
	  CopyTree(&Tnow,&Tnb[nb],false);

	  //printf("cut= %d Add= %d\n",Tnow.mv[nb].cut_id,Tnow.mv[nb].add_id);

	  MoveTree(g,&Tnb[nb],Tnow.mv[nb].cut_id,Tnow.mv[nb].add_id);//Cut and Add

	  //move[nb].score=QualityTree(g,&Tnb[nb]);
	  move[nb].score=QualityTreeDP(g,&Tnb[nb], sfg,&dpm,false);
	  //printf("**cut= %d Add= %d\n",Tnow.mv[nb].cut_id,Tnow.mv[nb].add_id);

	  //New!!
	  move[nb].cut_id=Tnow.mv[nb].cut_id;
	  move[nb].add_id=Tnow.mv[nb].add_id;

	  Tnb[nb].score=move[nb].score;
	  //printf(" score: %f\n",move[nb].score);

	  double ken;
	  if(i>0){
	   for(int j=0;j<i;j++){
	    ken=Kendall(Tnb[nb].Path,Tnb[nb].Lpath,Tr[j].Path,Tr[j].Lpath);
	    //printf("Kendall= %f\n",ken);
	    if(fabs(ken)>0.99){
	     move[nb].score=0;
	     Tnb[nb].score=0;
	     break;
	    }
	   }
	  }
	 }

	 //check best movement
	 double best_sco=0;
	 int best_mv=0;
	 for(int nb=0;nb<Nnb;nb++){
	  if(best_sco<move[nb].score && move[nb].score !=pre_best ){
	   best_sco=move[nb].score;
	   best_mv=nb;
	  }
	 }
	 //Update Sim Best
	 if(n==0||Tr[i].score < Tnb[best_mv].score){
	  CopyTree(&Tnb[best_mv],&Tr[i],false);
	 }
	 //Update Best
	 if(best_sco > Tbest.score){
	  CopyTree(&Tnb[best_mv],&Tbest,false);
	  Tbest.score=best_sco;
	 }
	 printf("#Sim: %3d/%3d Round= %4d/%4d Nmv= %6d RoundBestSco= %.1f TotalBest= %.1f\n",i+1,Nsim,n+1,Nround,Tnow.Nmv,best_sco,Tbest.score);
	 //Update Now
	 CopyTree(&Tnb[best_mv],&Tnow,false);

	 
	 pre_best=best_sco;

	 //Update Tabu List
	 tabu[tbp]=move[best_mv].cut_id;
	 if(Ntb<Ntabu) Ntb++;

	 tbp++;
	 if(tbp>=Ntb) tbp=0;

	 tabu[tbp]=move[best_mv].add_id;
	 if(Ntb<Ntabu) Ntb++;
	 tbp++;
	 if(tbp>=Ntb) tbp=0;
	 //printf("tbp=%d Ntb=%d\n",tbp,Ntb);
  	}
 }


 return false;
}


bool InitTree(TREE *t,bool mode) {


 printf("#Nnode= %d\n",t->Nnode);
 printf("#Etotal= %d Ne=%d\n",t->Etotal,t->Ne);
 if(mode){


  if((t->ActN=(bool *)malloc(sizeof(bool)*t->Nnode))==NULL)
   return true;
  if((t->MaskN=(bool *)malloc(sizeof(bool)*t->Nnode))==NULL)
   return true;
  //puts("#end ActN");
  if((t->ActE=(bool *)malloc(sizeof(bool)*t->Etotal))==NULL)
   return true;
  if((t->AddTbl=(int *)malloc(sizeof(int)*t->Etotal))==NULL)
   return true;
  if((t->CutTbl=(int *)malloc(sizeof(int)*t->Etotal))==NULL)
   return true;
  //puts("#end ActE");


  if((t->stock=(int *)malloc(sizeof(int)*t->Nnode))==NULL)
   return true;
  if((t->nextv=(int *)malloc(sizeof(int)*t->Nnode))==NULL)
   return true;
  if((t->cid=(int *)malloc(sizeof(int)*t->Nnode))==NULL)
   return true;
  //puts("#end nectv");

  if((t->cost=(double *)malloc(sizeof(double)*t->Nnode))==NULL)
   return true;

  if((t->node=(NODE *)malloc(sizeof(NODE)*t->Nnode))==NULL)
   return true;

  //puts("#end node");
  if((t->mv=(MOVE *)malloc(sizeof(MOVE)*t->Ne*t->Ne))==NULL)
   return true;

  if((t->Path=(int *)malloc(sizeof(int)*t->Ne))==NULL)
   return true;

 }
 puts("#Fin Malloc");
 for(int i=0;i<t->Nnode;i++){
  t->node[i].N=0;
  t->cost[i]=0.00;
  t->ActN[i]=false;
 }
 for(int i=0;i<t->Etotal;i++)
  t->ActE[i]=false;
 t->len=0;
 t->Nstock=0;
 t->Lpath=0;
 return false;
}

double QualityTreeDP_ori(GRAPH *g,TREE *t,SEQFG *sfg, DP_MEMORY *dpm, bool flag){
 int i,j;

 //DFS
 int st=t->St;
 int Nst=1;
 int v,w,maxi,maxi2;
 double maxd=0;
 double sco=0;
 int pre_st=-11;
 //printf("####St= %d\n",st);
 while(1){
  t->stock[0]=st;
  Nst=1;
  t->nextv[st]=-1;
  //init ActN
  //printf("Total Node= %d\n",t->Nnode);
  for(i=0;i<t->Nnode;i++){
   t->cost[i]=0.00;
   t->ActN[i]=false;
  }
  maxd=0;
  maxi=-1;
	while(1){
	 if(Nst==0)
	  break;
  	 v=t->stock[Nst-1];//pop
	 Nst--;
   	 //printf("#pop %d %d/%d Act? %d\n",Nst,v,t->Nnode,t->ActN[v]);
	 if(t->ActN[v])
	  continue;
	 t->ActN[v]=true;
	 //branch
	 //printf("Nbra= %d\n",g->node[v].N);
	 for(j=0;j<g->node[v].N;j++){
	  int eid=g->node[v].e[j]->eid;
	  //printf("check eid %d %d\n",eid,t->ActE[eid]);
	  if(t->ActE[eid]==false) continue;

	  w=g->node[v].e[j]->id1;
	  //printf("w1= %d or %d\n",w,g->node[v].e[j]->id2);
	  if(t->ActN[w]) continue;
	  //push
	  t->stock[Nst]=w;
	  Nst++;

	  t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
	  t->nextv[w]=v;//path

	  //printf("next1 %d -> %d or %d\n",v,w,g->node[v].e[j]->id2);
	  if(maxd<t->cost[w]){
	   maxd=t->cost[w];
	   maxi=w;
	  }
	 }
	 for(j=0;j<g->node[v].N;j++){
 	  int eid=g->node[v].e[j]->eid;
	  //printf("eid: %d %d\n",eid,t->ActE[eid]);
	  if(t->ActE[eid]==false) continue;

	  w=g->node[v].e[j]->id2;
	  //printf("w2= %d or %d\n",w,g->node[v].e[j]->id1);
	  if(t->ActN[w]) continue;
	  //push
	  t->stock[Nst]=w;
	  Nst++;

	  t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
	  t->nextv[w]=v;//path
	  //printf("next2 %d -> %d\n",v,w);

	  if(maxd<t->cost[w]){
	   maxd=t->cost[w];
	   maxi=w;
	  }
	 }
	}
  if(maxi==pre_st){
   break;
  }
  pre_st=st;
  st=maxi;
 }
 if(maxi==-1)
  return 0;

 //Set Start&End
 t->St=maxi;
 t->Ed=st;//!!! Aug-24

 sco=maxd*maxd;

 //Set 1st path
 for(i=0;i<t->Nnode;i++)
  t->ActN[i]=false;

 t->ActN[maxi]=true;
 int now=maxi;
 t->Path[0]=now;
 t->Lpath=1;
 while(1){
  now=t->nextv[now];
  t->ActN[now]=true;

  t->Path[t->Lpath]=now;
  t->Lpath++;

  //printf("#set true %d\n",now);
  if(now==st)
   break;
 }

 NODE **n_path_po=dpm->path;
 float *Smtx=dpm->Smtx;
 float *SmtxRv=dpm->SmtxRv;//Reverse sequence
 float *dtbl=dpm->dtbl;
 DPMTX *Dmtx=dpm->Dmtx;
 float Saa, Satm,Sss;
 float Waa, Watm,Wss;
 int n1=t->Lpath;
 int n2=sfg->l;

 //Waa=Watm=Wss=1.00;
 //parameters
 Waa=cmd.Waa;
 Wss=cmd.Wss;
 Watm=cmd.Watm;
 float ScoCut=cmd.LowestSco;

 printf("Npath= %d\n",t->Lpath);
 for(int i=0;i<t->Lpath;i++){
  n_path_po[i]=&(g->node[t->Path[i]]);
 }
 dtbl[0]=0;
 for(int i=0;i<t->Lpath-1;i++){
  dtbl[i+1]=dtbl[i]+sqrt(
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])*
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])+
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])*
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])+
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2])*
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2]));
 }
 for(int i=0;i<t->Lpath;i++){
        for(int j=0;j<sfg->l;j++){
         Saa =  n_path_po[i]->LogAA[sfg->seq[j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[j][2];
         Smtx[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	 //printf("%f, %f\n",n_path_po[i]->LogAA[20],n_path_po[i]->LogAA[19]);
	 //printf("[%d %d] Aa %d %f %f %f %f\n",i,j,sfg->seq[j],Saa,Satm,Sss,Smtx[GRID2D(n1,n2,i,j)]);
        }
 }
 //reverse sequence
 for(int i=0;i<t->Lpath;i++){
        for(int j=0;j<sfg->l;j++){

         Saa =  n_path_po[i]->LogAA[sfg->seq[sfg->l-1-j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[sfg->l-1-j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[sfg->l-1-j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[sfg->l-1-j][2];
         SmtxRv[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	 //printf("[%d %d] Aa %d %f %f %f %f\n",i,j,sfg->seq[j],Saa,Satm,Sss,Smtx[GRID2D(n1,n2,i,j)]);
        }
 }

 //Iter-DP
 float score,score2,final_score=0.00;
 int ali2[2000],ali2rv[2000];
 //init model
 for(int i=0;i<sfg->l;i++)
  dpm->ali2[i]=-1;
  //dpm->model[i]=&(g->node[0]);


 for(int iter=0;iter < 20; iter++){
 	score = dp_local(Dmtx,Smtx,dtbl,-100.0,-100.0,dpm->ali1,t->Lpath,ali2,sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);
 	score2 = dp_local(Dmtx,SmtxRv,dtbl,-100.0,-100.0,dpm->ali1,t->Lpath,ali2rv,sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);

	printf("%d %f %f\n",iter,score,score2);

	if(score <ScoCut && score2 < ScoCut)
	 break;

 	if(score2<score){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){
		 int j = ali2[pos];
			if(j!=-1){
			 //printf("[%d %d]",pos,dpm->ali2[pos]);
			 //mask
			 for(int i=0;i<t->Lpath;i++) Smtx[GRID2D(n1,n2,i,pos)]=-1;
			 for(int i=0;i<sfg->l;i++)   Smtx[GRID2D(n1,n2,j,i)]=-1;
			 for(int i=0;i<t->Lpath;i++) SmtxRv[GRID2D(n1,n2,i,sfg->l-1-pos)]=-1;
			 for(int i=0;i<sfg->l;i++)   SmtxRv[GRID2D(n1,n2,j,i)]=-1;
		 	 dpm->ali2[pos]=j;
		 	 dpm->model[pos]=n_path_po[j];
			}
		}
		//puts("");
	 final_score += score;
 	}else{
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){
		 int j = ali2rv[pos];//Revese
			if(j!=-1){
			 //printf("R[%d %d]",pos,ali2rv[pos]);
			 //mask
			 for(int i=0;i<t->Lpath;i++) Smtx[GRID2D(n1,n2,i,sfg->l-1-pos)]=-1;
			 for(int i=0;i<sfg->l;i++)   Smtx[GRID2D(n1,n2,j,i)]=-1;
			 for(int i=0;i<t->Lpath;i++) SmtxRv[GRID2D(n1,n2,i,pos)]=-1;
			 for(int i=0;i<sfg->l;i++)   SmtxRv[GRID2D(n1,n2,j,i)]=-1;
		 	 dpm->ali2[sfg->l-1-pos]=j;
		 	 dpm->model[sfg->l-1-pos]=n_path_po[j];
			}
		}
		//puts("");
	 final_score += score2;
	}
 }
 t->score=final_score;
 return final_score;
}

double QualityPathDP(GRAPH *g,SEQFG *sfg, DP_MEMORY *dpm, bool flag){
 int i,j;

 NODE **n_path_po=dpm->path;
 float *Smtx=dpm->Smtx;
 float *SmtxRv=dpm->SmtxRv;//Reverse sequence
 float *dtbl=dpm->dtbl;
 float *dtblR=dpm->dtblR;//Reverse
 DPMTX *Dmtx=dpm->Dmtx;
 float Saa, Satm,Sss;
 float Waa, Watm,Wss;
 int n1=g->Nnode;
 int n2=sfg->l;

 Waa=Watm=Wss=1.00;
 Waa=0.1;

 printf("#Npath= %d\n",g->Nnode);
 for(int i=0;i<g->Nnode;i++){
  n_path_po[i]=&(g->node[i]);
 }
 dtbl[0]=0;
 /*Ver02
 for(int i=0;i<g->Nnode-1;i++){
  dtbl[i+1]=dtbl[i]+sqrt(
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])*
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])+
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])*
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])+
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2])*
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2]));
 }
 */
 //Ver03 New dtbl position,position-j
 for(int i=0;i<g->Nnode;i++){
	for(int j=0;j<21;j++)
	 dtbl[i*21+j]=10.00;//init
	for(int j=1;j<21 && i-j>-1;j++){
  	 dtbl[i*21+j]=sqrt(
	  (n_path_po[i]->real_cd[0]-n_path_po[i-j]->real_cd[0])*
	  (n_path_po[i]->real_cd[0]-n_path_po[i-j]->real_cd[0])+
	  (n_path_po[i]->real_cd[1]-n_path_po[i-j]->real_cd[1])*
	  (n_path_po[i]->real_cd[1]-n_path_po[i-j]->real_cd[1])+
	  (n_path_po[i]->real_cd[2]-n_path_po[i-j]->real_cd[2])*
	  (n_path_po[i]->real_cd[2]-n_path_po[i-j]->real_cd[2]));
  	}
 }
 for(int i=0;i<g->Nnode;i++){
  int I = g->Nnode -1 - i;
	for(int j=0;j<21;j++)
	 dtblR[i*21+j]=10.00;//init
	for(int j=1;j<21 && i-j>-1;j++){
  	 dtblR[i*21+j]=sqrt(
	  (n_path_po[I]->real_cd[0]-n_path_po[I+j]->real_cd[0])*
	  (n_path_po[I]->real_cd[0]-n_path_po[I+j]->real_cd[0])+
	  (n_path_po[I]->real_cd[1]-n_path_po[I+j]->real_cd[1])*
	  (n_path_po[I]->real_cd[1]-n_path_po[I+j]->real_cd[1])+
	  (n_path_po[I]->real_cd[2]-n_path_po[I+j]->real_cd[2])*
	  (n_path_po[I]->real_cd[2]-n_path_po[I+j]->real_cd[2]));
  	}
 }
 for(int i=0;i<g->Nnode;i++){
        for(int j=0;j<sfg->l;j++){
         Saa =  n_path_po[i]->LogAA[sfg->seq[j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[j][2];
         Smtx[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
        }
 }
 puts("##Set Smtx");
 //reverse sequence
 for(int i=0;i<g->Nnode;i++){
        for(int j=0;j<sfg->l;j++){

         Saa =  n_path_po[i]->LogAA[sfg->seq[sfg->l-1-j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[sfg->l-1-j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[sfg->l-1-j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[sfg->l-1-j][2];
         SmtxRv[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
        }
 }

 //Iter-DP
 float score,score2,final_score=0.00;
 int ali2[2000],ali2rv[2000];
 //init model
 for(int i=0;i<sfg->l;i++)
  dpm->ali2[i]=-1;
  //dpm->model[i]=&(g->node[0]);


 for(int iter=0;iter < 1; iter++){
 	score = dp_global(Dmtx,Smtx,dtbl,-100.0,-100.0,dpm->ali1,g->Nnode,ali2,sfg->l,dpm->gali,&(dpm->Lgali),true);
 	score2 = dp_global(Dmtx,SmtxRv,dtblR,-100.0,-100.0,dpm->ali1,g->Nnode,ali2rv,sfg->l,dpm->gali,&(dpm->Lgali),true);

	//printf("%d %f %f\n",iter,score,score2);

	//if(score <20.0 && score2 < 20.0)
	// break;

 	if(score2<score){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){
		 int j = ali2[pos];
			if(j!=-1){
			 //printf("[%d %d]",pos,dpm->ali2[pos]);
			 //mask
			 for(int i=0;i<g->Nnode;i++) Smtx[GRID2D(n1,n2,i,pos)]=-1;
			 for(int i=0;i<sfg->l;i++)   Smtx[GRID2D(n1,n2,j,i)]=-1;
			 for(int i=0;i<g->Nnode;i++) SmtxRv[GRID2D(n1,n2,i,sfg->l-1-pos)]=-1;
			 for(int i=0;i<sfg->l;i++)   SmtxRv[GRID2D(n1,n2,j,i)]=-1;
		 	 dpm->ali2[pos]=j;
		 	 dpm->model[pos]=n_path_po[j];
			}
		}
		//puts("");
	 final_score += score;
 	}else{
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){
		 int j = ali2rv[pos];//Revese
			if(j!=-1){
			 //printf("R[%d %d]",pos,ali2rv[pos]);
			 //mask
			 for(int i=0;i<g->Nnode;i++) Smtx[GRID2D(n1,n2,i,sfg->l-1-pos)]=-1;
			 for(int i=0;i<sfg->l;i++)   Smtx[GRID2D(n1,n2,j,i)]=-1;
			 for(int i=0;i<g->Nnode;i++) SmtxRv[GRID2D(n1,n2,i,pos)]=-1;
			 for(int i=0;i<sfg->l;i++)   SmtxRv[GRID2D(n1,n2,j,i)]=-1;
		 	 dpm->ali2[sfg->l-1-pos]=j;
		 	 dpm->model[sfg->l-1-pos]=n_path_po[j];
			}
		}
		//puts("");
	 final_score += score2;
	}
 }
 return final_score;
}




double QualityTree(GRAPH *g,TREE *t){
 int i,j;

 //DFS
 int st=t->St;
 int Nst=1;
 int v,w,maxi,maxi2;
 double maxd=0;
 double sco=0;
 int pre_st=-11;
 //printf("####St= %d\n",st);
 while(1){
  t->stock[0]=st;
  Nst=1;
  t->nextv[st]=-1;
  //init ActN
  //printf("Total Node= %d\n",t->Nnode);
  for(i=0;i<t->Nnode;i++){
   t->cost[i]=0.00;
   t->ActN[i]=false;
  }
  maxd=0;
  maxi=-1;
	while(1){
	 if(Nst==0)
	  break;
  	 v=t->stock[Nst-1];//pop
	 Nst--;
   	 //printf("#pop %d %d/%d Act? %d\n",Nst,v,t->Nnode,t->ActN[v]);
	 if(t->ActN[v])
	  continue;
	 t->ActN[v]=true;
	 //branch
	 //printf("Nbra= %d\n",g->node[v].N);
	 for(j=0;j<g->node[v].N;j++){
	  int eid=g->node[v].e[j]->eid;
	  //printf("check eid %d %d\n",eid,t->ActE[eid]);
	  if(t->ActE[eid]==false) continue;

	  w=g->node[v].e[j]->id1;
	  //printf("w1= %d or %d\n",w,g->node[v].e[j]->id2);
	  if(t->ActN[w]) continue;
	  //push
	  t->stock[Nst]=w;
	  Nst++;

	  t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
	  t->nextv[w]=v;//path

	  //printf("next1 %d -> %d or %d\n",v,w,g->node[v].e[j]->id2);
	  if(maxd<t->cost[w]){
	   maxd=t->cost[w];
	   maxi=w;
	  }
	 }
	 for(j=0;j<g->node[v].N;j++){
 	  int eid=g->node[v].e[j]->eid;
	  //printf("eid: %d %d\n",eid,t->ActE[eid]);
	  if(t->ActE[eid]==false) continue;

	  w=g->node[v].e[j]->id2;
	  //printf("w2= %d or %d\n",w,g->node[v].e[j]->id1);
	  if(t->ActN[w]) continue;
	  //push
	  t->stock[Nst]=w;
	  Nst++;

	  t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
	  t->nextv[w]=v;//path
	  //printf("next2 %d -> %d\n",v,w);

	  if(maxd<t->cost[w]){
	   maxd=t->cost[w];
	   maxi=w;
	  }
	 }
	}
  if(maxi==pre_st){
   break;
  }
  pre_st=st;
  st=maxi;
 }
 if(maxi==-1)
  return 0;

 //Set Start&End
 t->St=maxi;
 t->Ed=st;//!!! Aug-24

 sco=maxd*maxd;

 //Set 1st path
 for(i=0;i<t->Nnode;i++)
  t->ActN[i]=false;

 t->ActN[maxi]=true;
 int now=maxi;
 t->Path[0]=now;
 t->Lpath=1;
 while(1){
  now=t->nextv[now];
  t->ActN[now]=true;

  t->Path[t->Lpath]=now;
  t->Lpath++;

  //printf("#set true %d\n",now);
  if(now==st)
   break;
 }
 //printf("#1st= %f\n",sco);
 //return sco;
 //Search 2nd, 3rd,... Longest path in tree
 int besti;
 double bestc,diffc;

 //for(int b=2;b<t->Nnode;b++){
 for(int b=2;b<100;b++){
  besti=-1;
  bestc=0;
  
  	for(int i=0;i<t->Nnode;i++){
	 if(t->ActN[i])
	  continue;

	////////
	 if(g->node[i].N==0)
	  continue;
	///////

	 //trace back
	 now=i;
	 while(1){
	  //printf("now= %d -> %d %d\n",now,t->nextv[now],t->node[now].N);
	  now=t->nextv[now];
  	  if(t->ActN[now])
	   break;
	 }
	 diffc=t->cost[i]-t->cost[now];
	 if(bestc<diffc){
	  bestc=diffc;
	  besti=i;
	 }
	}
  if(besti==-1)
   break;

  sco+=bestc*bestc;
  //printf("#%d= %f\n",b,sco);
  //trace back again
  now=besti;
	 while(1){
	  //printf("now=%d\n",now);
	  now=t->nextv[now];
  	  if(t->ActN[now])
	   break;
	  t->ActN[now]=true; //!!!!!! New Aug-24
	 }

 }
 t->score=sco;
 return sco;
}

bool CopyTree(TREE *in, TREE *out,bool mode){
 int i,j;
 out->len=in->len;
 out->Nnode=in->Nnode;
 out->Ne=in->Ne;
 out->Ntotal=in->Ntotal;
 out->Etotal=in->Etotal;
 out->St=in->St;
 out->Ed=in->Ed;
 out->Nadd=in->Nadd;
 out->Ncut=in->Ncut;
 out->score=in->score;
 out->Lpath=in->Lpath;

 if(mode){
  if((out->node=(NODE *)malloc(sizeof(NODE)*out->Nnode))==NULL)
   return true;
  if((out->ActN=(bool *)malloc(sizeof(bool)*out->Nnode))==NULL)
   return true;
  if((out->MaskN=(bool *)malloc(sizeof(bool)*out->Nnode))==NULL)
   return true;
  if((out->UsedN=(bool *)malloc(sizeof(bool)*out->Nnode))==NULL)
   return true;
  if((out->ActE=(bool *)malloc(sizeof(bool)*out->Etotal))==NULL)
   return true;
  if((out->stock=(int *)malloc(sizeof(int)*out->Nnode))==NULL)
   return true;
  if((out->nextv=(int *)malloc(sizeof(int)*out->Nnode))==NULL)
   return true;
  if((out->cid=(int *)malloc(sizeof(int)*out->Nnode))==NULL)
   return true;
  if((out->cost=(double *)malloc(sizeof(double)*out->Nnode))==NULL)
   return true;
  if((out->CutTbl=(int *)malloc(sizeof(int)*out->Etotal))==NULL)
   return true;
  if((out->AddTbl=(int *)malloc(sizeof(int)*out->Etotal))==NULL)
   return true;
  //if((out->mv=(MOVE *)malloc(sizeof(MOVE)*out->Etotal*out->Etotal))==NULL)
  if((out->mv=(MOVE *)malloc(sizeof(MOVE)*100*out->Etotal))==NULL)
   return true;
  if((out->Path=(int *)malloc(sizeof(int)*out->Ne))==NULL)
   return true;
 }

 for(i=0;i<in->Lpath;i++){
  out->Path[i]=in->Path[i];
 }

 for(i=0;i<in->Etotal;i++){
  out->ActE[i]=in->ActE[i];
 }
/*
 for(i=0;i<in->Nnode;i++){
  out->node[i].N=in->node[i].N;
  out->node[i].cid=in->node[i].cid;
	for(j=0;j<in->node[i].N;j++)
	 out->node[i].e[j]=in->node[i].e[j];
 }
*/
 return false;
}

bool SetCutTbl(GRAPH *g,TREE *t,int *tabu,int Ntb,double LenCutoff){
 int i,j,n=0;
 
 for(i=0;i<g->Ne;i++){

  int id1=g->edge[i].id1;
  int id2=g->edge[i].id1;

  if(g->node[id1].N==1) continue;
  if(g->node[id2].N==1) continue;

  //Active & not keep
  if(t->ActE[i] && g->edge[i].keep==false){
   bool flag=false;

   for(j=0;j<Ntb;j++){
    if(tabu[j]==i){
     flag=true;
     break;
    }
   }
   if(flag==true)
    continue;


//   t->CutTbl[n]=i;
//   n++;

   double len=t->len-g->edge[i].d;
	//Set Addtbl
	SplitChain(g,t,i);
	
 	for(j=0;j<g->Ne;j++){
 	 if(t->ActE[j])
 	  continue;
 	 if(t->cid[g->edge[j].id1]*t->cid[g->edge[j].id2]!=-1)
 	  continue;
	 if(g->edge[j].local==false && g->edge[j].mst==false)
	  continue;

  	//Total length restraints
  	 if(len+g->edge[j].d>LenCutoff)
  	  continue;

	 //tabu list
	 bool flag=false;
  	 for(int k=0;k<Ntb;k++){
  	  if(tabu[k]==j){
  	   flag=true;
  	   break;
  	  }
  	 }
  	 if(flag==true)
   	  continue;
	 if(i==j)
		 continue;

	 
	 t->mv[n].cut_id=i;
	 t->mv[n].add_id=j;
	 n++;
 	}



  }
 }
 t->Nmv=n;
 //printf("Nmv= %d\n",t->Nmv);
 if(n>2)
  ShuffleMv(t->mv,n);

 return false;
}

bool MoveTree(GRAPH *g,TREE *t,int cut_id,int add_id){

 //Update len
 t->len=t->len - g->edge[cut_id].d + g->edge[add_id].d;

 //Update ActiveE
 t->ActE[cut_id]=false;
 t->ActE[add_id]=true;

}


//Cut and..
int CutTree(GRAPH *g,TREE *t,double LenCutoff,int *tabu,int Ntb){
 int p,i,j,w,n=0;
 int Nst=0;
 int v,Nadd,cut_id;
 double len=t->len;

 for(p=0;p<t->Ncut;p++){
  cut_id=t->CutTbl[p];
  //printf("#cut %d/%d cut_id= %d\n",p+1,t->Ncut,cut_id);

  len=t->len-g->edge[cut_id].d;

  SplitChain(g,t,cut_id);

 	//Set AddTbl
 	Nadd=0;
 	for(i=0;i<g->Ne;i++){
 	 if(t->ActE[i])
 	  continue;
 	 if(t->cid[g->edge[i].id1]*t->cid[g->edge[i].id2]!=-1)
 	  continue;
	 if(g->edge[i].local==false && g->edge[i].mst==false)
	  continue;
  	//Total length restraints
  	 if(len+g->edge[i].d>LenCutoff)
  	  continue;

	 //tabu list
	 bool flag=false;
  	 for(j=0;j<Ntb;j++){
  	  if(tabu[j]==i){
  	   flag=true;
  	   break;
  	  }
  	 }
  	 if(flag==true)
   	  continue;

  	 t->AddTbl[Nadd]=i;
  	 Nadd++;
 	}
 
  //printf("Nadd= %d\n",Nadd);

  if(Nadd>0)
   break;
 }

 if(Nadd==0)
  return -1;

 //Update
 t->Nadd=Nadd;
 t->len=len;

 //Shuffle Addtbl
 //for(i=0;i<Nadd;i++)
 // printf("%6d",t->AddTbl[i]);
 //puts("");

 ShuffleTbl(t->AddTbl,Nadd);

 //for(i=0;i<Nadd;i++)
 // printf("%6d",t->AddTbl[i]);
 //puts("");

 

 t->ActE[cut_id]=false;
 return cut_id;
}

int SplitChain(GRAPH *g,TREE *t,int cut_id){
 int i,j,Nst;
 int v,eid,w;
 int id1=g->edge[cut_id].id1;
 int id2=g->edge[cut_id].id2;
 int N0,N1;
 N0=N1=0;
 t->stock[0]=id1;
 Nst=1;

 //init ActN
 for(i=0;i<g->Nnode;i++){
  t->ActN[i]=false;
  t->cid[i]=0;
 }

 t->ActN[id2]=true;
 t->cid[id1]=-1;

  	//chain id1 : -1
  	while(1){
 	 if(Nst==0)
   	  break;
   	 v=t->stock[Nst-1];//pop
   	 //printf("Set cid %d : true\n",v);
   	 Nst--;
   	 if(t->ActN[v])
   	  continue;

   	 t->ActN[v]=true;
   	 t->cid[v]=-1;
	 N0++;
  	 //branch
  	 	for(j=0;j<g->node[v].N;j++){
	 	 eid=g->node[v].e[j]->eid;
	 	 if(t->ActE[eid]==false) continue;
  	 	 w=g->node[v].e[j]->id1;
  	 	 if(t->ActN[w]==false){
	 	  //push
	 	  t->stock[Nst]=w;
	 	  Nst++;
		 }
		 w=g->node[v].e[j]->id2;
		 if(t->ActN[w]==false){
	 	  //push
	 	  t->stock[Nst]=w;
	 	  Nst++;
		 }
		}
 	}
 	//chain id2
 	t->stock[0]=id2;
 	Nst=1;
 	t->ActN[id1]=true;
 	t->ActN[id2]=false;
 	t->cid[id2]=1;
 	while(1){
 	 if(Nst==0)
 	  break;
 	 v=t->stock[Nst-1];//pop
 	 //printf("Set cid2 %d : false\n",v);
 	 Nst--;
 	 if(t->ActN[v])
 	  continue;
 	 t->ActN[v]=true;
 	 t->cid[v]=1;
	 N1++;
 	 //printf("Set cid3 %d : false\n",v);
  		//branch
  		for(j=0;j<g->node[v].N;j++){
		 int eid=g->node[v].e[j]->eid;
  		 //printf("push0 %d eid %d\n",w,eid);
		 if(t->ActE[eid]==false) continue;
  		 w=g->node[v].e[j]->id1;
  		 if(t->ActN[w]==false){
	  	  //push
		  t->stock[Nst]=w;
		  Nst++;
		 }
		 w=g->node[v].e[j]->id2;
		 if(t->ActN[w]==false){
	  	  //push
		  t->stock[Nst]=w;
		  Nst++;
		 }
		}
 	}
 //printf("N0=%d N1=%d\n",N0,N1);
 return 0;
}


int AddEdge(GRAPH *g, TREE *t){
 int eid=t->AddTbl[0];//First
 t->len+=g->edge[eid].d;
 t->ActE[eid]=true;
 return eid;
}


//0-0.999999...
double RandDouble(){
 double a=rand();
 return (double)a/(double)(RAND_MAX+1.0);
}

//0-max
int RandInt(int max){
 int a=rand();
 return a%(max+1);
}

void ShuffleTbl(int *t, int n){
 int i,j,tmp;

 if(n>2){

 for(i=0;i<n;i++){
  j=RandInt(n-2)+1+i;
  if(j>=n) j-=n;
  //printf("i= %d j= %d / %d\n",i,j,n);
  tmp=t[i];
  t[i]=t[j];
  t[j]=tmp;
 }
 }
}

void ShuffleMv(MOVE *t, int n){
 int i,j,tmp;
 MOVE tmp_mv;
 if(n>2){

 for(i=0;i<n;i++){
  j=RandInt(n-2)+1+i;
  if(j>=n) j-=n;
  //printf("i= %d j= %d / %d\n",i,j,n);
  tmp_mv.cut_id=t[i].cut_id;
  tmp_mv.add_id=t[i].add_id;

  t[i].cut_id=t[j].cut_id;
  t[i].add_id=t[j].add_id;

  t[j].cut_id=tmp_mv.cut_id;
  t[j].add_id=tmp_mv.add_id;
 }
 }
}


double Kendall(int *a,int Na,int *b,int Nb){
 int i,j,k;
 int Ntrue=0;
 int Nfalse=0;
 bool f[5000]={false};
 int Pa;

 //for(i=0;i<Na;i++) f1[a[i]]=true;
 for(i=0;i<Nb;i++) f[b[i]]=true;

 //for(i=0;i<Na-1;i++){
 for(i=0;i<Na-1;i+=2){
  if(f[a[i]]==false) continue;

  Pa=0;
  for(k=0;k<Nb;k++){
   if(b[k]==a[i]){
    Pa=k;
    break;
   }
  }



  //for(j=i+1;j<Na;j++){
  for(j=i+1;j<Na;j+=2){
   if(f[a[j]]==false) continue;

   for(k=Pa;k<Nb;k++){
    if(b[k]==a[j]){
     Ntrue++;
     break;
    }
   }
   //if(k==Nb) 
   // Nfalse++;

  }
 }

 //printf("Nt= %d Nf= %d\n",Ntrue,Nfalse);
 //return (double)(2*Ntrue)/(double)(Ntrue+Nfalse)-1.00;
 return (double)(2*Ntrue)/(double)(Nb*(Nb-1)*0.5)-1.00;
}


bool PairExhaust(GRAPH *g,TREE *init_t,TREE *Tr){
 int Nround=cmd.Nround;
 int Nnb=cmd.Nnb;
 int Ntabu=cmd.Ntabu;
 int Nsim=cmd.Nsim;
 int Nbeam=cmd.Nbeam;
 double AllowLen=init_t->len*cmd.Allow;
 MOVE *move;
 int tabu[1000];
 int Ntb=0;
 int tbp=0;
 TREE Tnow,Tbest,Tnb[500];
 int *pair;

 if((pair=(int *)malloc(sizeof(int)*g->Nnode*g->Nnode*2))==NULL)
  return true;


 //set tmp_tree
 int Nth=omp_get_max_threads();
 //printf("Set Nth= %d\n",Nth);
 //return true;
 TREE **tmp_tree;
 if((tmp_tree=(TREE **)malloc(sizeof(TREE *)*Nth))==NULL)
  return true;
 for(int i=0;i<Nth;i++){
  if((tmp_tree[i]=(TREE *)malloc(sizeof(TREE)*(Nbeam*2+1)))==NULL)
   return true;
  for(int j=0;j<Nbeam*2+1;j++)
   CopyTree(init_t,&tmp_tree[i][j],true);

 }

 //Set Start-End Pairs
 int Npair=ListStEd(g,init_t,pair);
 printf("#Find Pairs...Npair= %d\n",Npair);

 //OptPath2points(g,init_t,tmp_tree[0],Nbeam,415,553);
 //CopyTree(&tmp_tree[0][2*Nbeam],&Tr[0],true);
 

 //return false;

 #pragma omp parallel for schedule(dynamic,50)
 for(int i=0;i<Npair;i++){
  int th=omp_get_thread_num();
  //printf("%d th=%d %d:%d\n",i+1,th,pair[2*i],pair[2*i+1]);
  OptPath2points(g,init_t,tmp_tree[th],Nbeam,pair[2*i],pair[2*i+1]);
 }

 //copy Best
 for(int i=0;i<Nth;i++){
  CopyTree(&tmp_tree[i][2*Nbeam],&Tr[i],true);
 }
 puts("#FIn Search");
 return false;
}

int ListStEd(GRAPH *g,TREE *t,int *p){
 int i,j,k;
 int Npair=0;

 for(i=0;i<g->Nnode;i++){

  bool mst_flag1=false;
  int Nbranch=0;
  for(k=0;k<g->Ne;k++){
   	//within MST?
	if(t->ActE[k]==true && (g->edge[k].id1==i ||g->edge[k].id2==i))
	 Nbranch++;
  }
  if(Nbranch != 1)
   continue;
  for(j=i+1;j<g->Nnode;j++){
   bool flag=false;
   bool mst_flag2=false;
   Nbranch=0;
	for(k=0;k<g->Ne;k++){
	 if(g->edge[k].id1==i && g->edge[k].id2==j){
	  flag=true;
	  break;
	 }
	 if(g->edge[k].id2==i && g->edge[k].id1==j){
	  flag=true;
	  break;
	 }
	 //within MST?
	 if(t->ActE[k]==true && (g->edge[k].id1==j ||g->edge[k].id2==j))
	  Nbranch++;
	  //mst_flag2=true;
	}
   if(flag)
    continue;
   //if(mst_flag2==false)
   if(Nbranch!=1)
    continue;

   //printf("%d - %d \n",i,j);
   p[Npair*2]=i;
   p[Npair*2+1]=j;
   Npair++;
  }
 }
 //printf("Npair= %d\n",Npair);
 return(Npair);
}


double OptPath2points(GRAPH *g,TREE *t,TREE *tmp,int Nbeam,int St,int Ed){
 int cut_id,add_id,id1,id2;
 int cnt;
 double now_score;

 //init tree
 for(int i=0;i<Nbeam*2+1;i++){
  //printf("Copy %d\n",i);
  CopyTree(t,&tmp[i],false);
  tmp[i].score=0;
 }
 printf("Path...%d %d\n",St,Ed);

 //0: current tree data,Top Nbeam: 1-Nbeam
 for(int iter=0;iter<20;iter++){

  for(int beam=0;beam<Nbeam;beam++){

   if(iter!=0 && tmp[beam].score==0.0)
    continue;
   //Find path
   if(get_path(g,&tmp[beam],St,Ed)==false)
    return 0;
   //cut
   now_score=tmp[beam].score;

	for(int i=0;i<tmp[beam].Lpath;i++){
	 cut_id=tmp[beam].Path[i];

/*
	  if(tmp[beam].ActE[2567]==true)
	   printf("0T* Non act Path?? %d %d/%d %d\n",beam,i,tmp[beam].Lpath, cut_id);
	  else
	   printf("0F* Non act Path?? %d %d/%d %d\n",beam,i,tmp[beam].Lpath, cut_id);
 	 if(tmp[beam].ActE[cut_id]==false)
	  printf("?? Non act Path?? %d %d/%d %d\n",beam,i,tmp[beam].Lpath, cut_id);
*/
	 if(g->edge[cut_id].keep==true) 
	  continue;
	 SplitChain(g,&tmp[beam],cut_id);
	 //search possible new connections
	 for(int j=0;j<g->Ne;j++){
	  if(g->edge[j].local!=true && g->edge[j].mst!=true)
	   continue;
	  if(tmp[beam].ActE[j]==true)
	   continue;
	  if(cut_id==j) //same
	   continue;
	  id1=g->edge[j].id1;
	  id2=g->edge[j].id2;
	  if(tmp[beam].cid[id1]*tmp[beam].cid[id2]!=-1)
	   continue;

	  //if(tmp[beam].ActE[j]==true){ //????

	  // printf("%d pos= %d b=%d Why???cutid= %d eid= %d cid= %d %d\n",iter,i,beam,cut_id,j,tmp[beam].cid[id1],tmp[beam].cid[id2]);
	  //}
	  add_id=j;

	  tmp[beam].ActE[cut_id]=false;
	  tmp[beam].ActE[add_id]=true;

	  //printf("cut %d, add %d\n",cut_id,add_id);

	  //check distance
	  if(get_path(g,&tmp[beam],St,Ed)==false){
  	   return 0;
	  }

	  //if(tmp[beam].score > now_score)
	  // printf("%f > %f\n",tmp[beam].score,now_score);
	  //if(tmp[beam].score > now_score*0.9 && tmp[beam].score < now_score*1.1 )
	  // continue;
	  //Keep top Nbeam
	  //0..(Nbeam-1) current paths
	  //Nbeam..Nbeam*2 top Nbeam paths

	  if(tmp[beam].score > tmp[Nbeam*2-1].score){
	   //printf("Add!! %f\n",tmp[beam].score);
	   //ignore same score
	   bool chk=false;
	   for(int k=0;k<Nbeam;k++){
	    if(tmp[Nbeam+k].score==tmp[beam].score){
	     chk=true;
	    }
	   }
	   if(chk==false){
	    CopyTree(&tmp[beam],&tmp[Nbeam*2-1],false);
	    qsort(&tmp[Nbeam],Nbeam,sizeof(TREE),cmp_tree_score);

	    //check
	    //for(int ii=0;ii<Nbeam;ii++)
	    // printf("Update iter= %d b=%d %f\n",iter,beam,tmp[Nbeam+ii].score);
	   }
	  }

	  //restore
	  tmp[beam].ActE[cut_id]=true;
	  tmp[beam].ActE[add_id]=false;

/*
	  if(tmp[beam].ActE[2567]==true)
	   printf("*1T Non act Path?? %d %d/%d %d\n",beam,i,tmp[beam].Lpath, cut_id);
	  else
	   printf("*1F Non act Path?? %d %d/%d %d\n",beam,i,tmp[beam].Lpath, cut_id);
*/
	 }
	}
   if(iter==0)//first iter
    break;
  }
  //Copy Best
  if(tmp[Nbeam].score >tmp[2*Nbeam].score )
   CopyTree(&tmp[Nbeam],&tmp[2*Nbeam],false);
  //update tmp[0]..tmp[Nbeam-1]
  for(int ii=0;ii<Nbeam;ii++){
   if(tmp[Nbeam+ii].score>0.0)
   CopyTree(&tmp[Nbeam+ii],&tmp[ii],false);
   //init
   tmp[Nbeam+ii].score=0.0;
  }
  //for(int ii=0;ii<Nbeam;ii++)
  // printf("%d Update %f\n",iter,tmp[ii].score);
  //break;
 }
}


bool get_path(GRAPH *g, TREE *t,int s,int e){
 int Nst=0;
 int i,j,v,w,eid;
 t->St=s;
 t->Ed=e;
 int nexte[20000];
 //init
 for(i=0;i<g->Nnode;i++){
  t->ActN[i]=false;
  t->cost[i]=0;
 }
 t->stock[0]=s;
 Nst=1;
 t->nextv[s]=-1;

 while(1){
	 if(Nst==0)
	  break;
  	 v=t->stock[Nst-1];//pop
	 Nst--;
   	 //printf("#pop %d %d Act? %d\n",Nst,v,t->ActN[v]);
	 if(t->ActN[v])
	  continue;
	 if(v==e)
	  break;
	 t->ActN[v]=true;
	 //branch
	 for(j=0;j<g->node[v].N;j++){
	  eid=g->node[v].e[j]->eid;
	  //printf("check eid %d %d\n",eid,t->ActE[eid]);
	  if(t->ActE[eid]==false) continue;

	  w=g->node[v].e[j]->id1;
	  //printf("w1= %d or %d\n",w,g->node[v].e[j]->id2);
	  if(t->ActN[w]==false){
	   //push
	   t->stock[Nst]=w;
	   Nst++;

	   t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
	   t->nextv[w]=v;//path
	   nexte[w]=eid;//edge path
	   if(w==e)
	    break;
	  }

	  w=g->node[v].e[j]->id2;
	  if(t->ActN[w]==false){
	   //push
	   t->stock[Nst]=w;
	   Nst++;

	   t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
	   t->nextv[w]=v;//path
	   nexte[w]=eid;//edge path
	   //t->cid[w]=eid;//path
	   //printf("next2 %d -> %d\n",v,w);
	   if(w==e)
	    break;
	   
	  }
	 }
 }
 //printf("%d-%d Cost= %f\n",s,e,t->cost[e]);
 if(w!=e){
  printf("WARNING!! %d != %d\n",w,e);
  return false;
 }
 int now=e;

 //t->Path[0]=t->cid[now];
 t->Path[0]=nexte[now];
 t->Lpath=1;

  //printf("%d,",now);
 while(1){
  now=t->nextv[now];
  //printf("%d,",now);
  if(now==s)
   break;
  //if(t->ActE[2567]==true && nexte[now]==2567)
  // printf("WARNING!!! Non-actE %d\n",t->ActE[nexte[now]]);
  t->Path[t->Lpath]=nexte[now];
  t->Lpath++;
 }
 //puts("");

 //for(int i=0;i<t->Lpath;i++)
 // printf("%d,",t->Path[i]);
 //puts("");
 t->score=t->cost[e];
 return true;
}




//Make Graph and MST
int SearchFrag(int,GRAPH *,double ,int, FRAG_DB *);

bool FindFrag(POINTS *p, GRAPH *g, FRAG_DB *db, double scale, int Nres){
 int i,j,dif;
 double dcut=(double)Nres*3.70*1.10/scale;
 int Nmodel=20;//Max number of fragments for each node
 //g->Nnode=p->Ncd;
 int Ne=0;
 //Stock Data
 //FRAG_DB *db=malloc(sizeof(FRAG_DB)*g->Nnode);
 unsigned int TotalFrag =0;
 

 //Search fragments
 //#pragma omp parallel for schedule(dynamic,5)
 for(int i=0;i<g->Nnode;i++){
	//Ignore terminal
	/*
	if(g->node[i].N < 2){
	 db[i].N=0;
	 continue;
	}*/
 	if((db[i].frag = malloc(sizeof(FRAGMENT)*g->Nnode))==NULL)
 	 return -1;
	TotalFrag += SearchFrag(i,g,dcut,Nmodel,&db[i]);
	//printf("NID %d Nbranch= %d Nmodel= %d\n",i,g->node[i].N,db[i].N);
	
 }
 //
 printf("#Total Fragments = %d\n",TotalFrag);
}

float cos_vec(float a[3],float b[3]){

 float ab=a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
 float ra=1.000/ sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
 float rb=1.000/ sqrt(b[0]*b[0]+b[1]*b[1]+b[2]*b[2]);

 return ab*ra*rb;
}

int SearchFrag(int InitNode,GRAPH *g,double dcut,int N, FRAG_DB *db){
 //edge[]->d, edge->dens
 int Nmodel =0;

 int Nstock;

 EDGE *ar_e[1000];
 EDGE *path_e[1000];
 int Np=0;
 int ar_n[1000]={};
 int path_n[1000]={};
 int now = InitNode;
 int Ne,Nn;
 int MaxNOmst = 1;
 ar_n[0]=now; Nn=1;
 double dist = 0;
 //printf("INIT %d DCUT= %f\n",now,dcut);
 EDGE *e;
 Ne=0;
 Nn=0;
 path_n[0]=now;
 ar_e[Ne++]= NULL;//pop signal
 ar_n[Nn++]= -1;//pop signal
 for(int i=0;i<g->node[now].N;i++){
  if(g->node[now].e[i]->mst==false)//Terminal should be MST
   continue;
  //ignore too long edge
  if(g->node[now].e[i]->d > 5.0)
   continue;
  ar_e[Ne++]=g->node[now].e[i];//push
  int pid=g->node[now].e[i]->id1;//push
  if(g->node[now].e[i]->id1 == now)
   pid=g->node[now].e[i]->id2;//push
  ar_n[Nn++]=pid;
  //printf("%d %d\n",ar_e[Ne-1]->id1,ar_e[Ne-1]->id2);
 }
 
 int NOmst=0;
 while(1){
	if(Ne==0)
	 break;
  	//pop
  	e=ar_e[--Ne];
  	now = ar_n[--Nn];
	if(e==NULL){//pop signal
	 if(Np==0)
	  break;
	 Np--;
	 dist -= path_e[Np]->d;
	 continue;
	}
	path_n[Np+1]=now;
	path_e[Np++]=e;
 	dist+=e->d; 
  	//printf("now= %d %f Np= %d\n",now,dist,Np);
	NOmst=0;
	for(int j=0;j<Np;j++)
	 if(path_e[j]->mst==false)
	  NOmst++;
	if(NOmst>MaxNOmst){
	 Np--;
	 dist -= path_e[Np]->d;
	 continue;
	}
	if(dist > dcut){
	 	//puts("CHECK SCORE");
		if(e->mst==true){//terminal should be MST
/*
	 	 for(int j=0;j<=Np;j++)
		  printf("%d,",path_n[j]);
		 printf("|%d %f %d/%d\n",Np+1,dist,NOmst,Np);
*/
		 //Inpput data
		 for(int j=0;j<=Np;j++){
		  db->frag[Nmodel].n_path[j]=path_n[j];
		  db->frag[Nmodel].n_path_po[j] = &(g->node[path_n[j]]);//pointer
		 }
		 db->frag[Nmodel].dtbl[0]=0.00;
		 for(int j=0;j<Np;j++){
		  db->frag[Nmodel].e_path[j]=path_e[j];
		  db->frag[Nmodel].dtbl[j+1]=path_e[j]->d+db->frag[Nmodel].dtbl[j];
		 }
		 db->frag[Nmodel].Nnode=Np+1;
		 db->frag[Nmodel].dist=dist;
		 Nmodel++;
		}
	 Np--;
	 dist -= path_e[Np]->d;
	 continue;
	}
	
  	//stok
 	ar_e[Ne++]= NULL;//pop signal
 	ar_n[Nn++]= -1;//pop signal
  	for(int i=0;i<g->node[now].N;i++){
	 //mst only?
	 //if(g->node[now].e[i]->mst==false)
	 // continue;
	 if(g->node[now].e[i]->d>5.0)
          continue;
  	 int pid=g->node[now].e[i]->id1;//push
  	 if(now == pid)
   	  pid=g->node[now].e[i]->id2;//push
	 //check loop
		bool flag=false;
		for(int j=Np;j>-1;j--){
		 if (pid == path_n[j]){
		  flag=true;
		  break;
		 }
		}
		if(flag==true)
		 continue;
  	 ar_e[Ne++]=g->node[now].e[i];//push
  	 ar_n[Nn++]=pid;//push
  	 //printf("ADD %d %d\n",ar_e[Ne-1]->id1,ar_e[Ne-1]->id2);
  	}
	//No update?
	if(ar_e[Ne-1]==NULL){
	 Ne--;
	 Nn--;
	 Np--;
	 dist -= path_e[Np]->d;
	 continue;
	}
  	//printf("Ne= %d / %d\n",Ne,g->node[now].N);
 }
 db->N=Nmodel;

 //filtering

 //Tabu check if nodes have P(N,CA,C) < 0.8
 int Nh,Ns,Ntabu;
 float SumP=0.0;
 FRAGMENT *fg;
 float d12,d13,d23;
 float v21[3],v23[3];
 for(int i=0;i<Nmodel;i++){
	fg=&(db->frag[i]);
	Nh=Ns=Ntabu=0;
 	for(int j=0;j<fg->Nnode;j++){
 	 SumP=fg->n_path_po[j]->ATOMP[1]+fg->n_path_po[j]->ATOMP[2]+fg->n_path_po[j]->ATOMP[3];
 	 if(SumP<0.8)//backbone filter
 	  Ntabu++;
	 //Helix Only
	 if(fg->n_path_po[j]->SSP[0]<0.5)
 	  Nh++;
	 if(fg->n_path_po[j]->SSP[1]<0.5)
 	  Ns++;
 	}
	if(Ntabu > 0)
         fg->Nnode=0;

 	//if(Nh > 0 && Ns > 0)
 	// fg->Nnode=0;
	if(fg->Nnode==0)
         continue;

	//Remove kinked sheet
	for(int j=0;j<fg->Nnode-2;j++){
	 if(fg->n_path_po[j]->SSP[1]>0.5 
	 && fg->n_path_po[j+1]->SSP[1]>0.5 
	 && fg->n_path_po[j+2]->SSP[1]>0.5){
		v21[0]=fg->n_path_po[j]->real_cd[0]-fg->n_path_po[j+1]->real_cd[0];
		v21[1]=fg->n_path_po[j]->real_cd[1]-fg->n_path_po[j+1]->real_cd[1];
		v21[2]=fg->n_path_po[j]->real_cd[2]-fg->n_path_po[j+1]->real_cd[2];
		v23[0]=fg->n_path_po[j+2]->real_cd[0]-fg->n_path_po[j+1]->real_cd[0];
		v23[1]=fg->n_path_po[j+2]->real_cd[1]-fg->n_path_po[j+1]->real_cd[1];
		v23[2]=fg->n_path_po[j+2]->real_cd[2]-fg->n_path_po[j+1]->real_cd[2];

		if(cos_vec(v21,v23) > -0.5){ //cos 120 deg
		 Ntabu++;
		 break;
		}
	  }
	 }
	if(Ntabu > 0)
         fg->Nnode=0;

 }

 return Nmodel;
}


bool read_ptbl(PTBL *tbl,char *filename){
 FILE *fpin;
 char line[5000];
 int i=0;
 float p[100];
 int N=0;
 if((fpin=fopen(filename,"rb")) == NULL){ 
  fprintf(stderr,"Can't open %s\n",filename); 
  return(true); 
 }
 while(fgets(line,5000,fpin)){
  N++;
 }
 fclose(fpin);
 if((fpin=fopen(filename,"rb")) == NULL){ 
  fprintf(stderr,"Can't open %s\n",filename); 
  return(true); 
 }

 tbl->prob=(double **)malloc(sizeof(double*)*N);
 while(fgets(line,5000,fpin)){ 
  char* tmp = strdup(line);
  int j =0;
  const char* tok;
  tbl->prob[i]=(double *)malloc(sizeof(double)*100);
 	for(tok = strtok(line,",");tok && *tok; tok=strtok(NULL,",\n")){
	 tbl->prob[i][j]=atof(tok);
	 j++;
	}
  i++;
 }
 fclose(fpin);
 tbl->N = N;
 return false;
}


//Probability+MRC
bool AddProbToMrc(PTBL *tbl,MRC *mrc,float weight, float Pcut,float Dcut){
 int xydim = mrc->xdim*mrc->ydim;
 int xdim = mrc->xdim;
 int pos[3];

 //copy
 for(int x=0;x<mrc->xdim*mrc->ydim*mrc->zdim;x++){
  if(mrc->dens[x] < Dcut)
   mrc->dens[x]=-1.00;//Mask
 }


 //check MAX XYZ
 for(int i=0;i<tbl->N;i++){
  int x,y,z,ind;
	//if(tbl->prob[i][24]+tbl->prob[i][25]+tbl->prob[i][26]<Pcut)
	// continue;
  x=(int)(tbl->prob[i][0]); 
  y=(int)(tbl->prob[i][1]); 
  z=(int)(tbl->prob[i][2]); 
  ind=xydim*z+xdim*y+x;
  //23rd-25th N,CA,C
  //mrc->dens[ind]+=(tbl->prob[i][24]+tbl->prob[i][25]+tbl->prob[i][26])*weight;
  if(mrc->dens[ind] < 0.00)//Mask
   continue;
  mrc->dens[ind]=weight*mrc->dens[ind]+(tbl->prob[i][24]+2.0*tbl->prob[i][25]+tbl->prob[i][26]);
 }
 for(int x=0;x<mrc->xdim*mrc->ydim*mrc->zdim;x++){
  if(mrc->dens[x] < 0.00)
   mrc->dens[x]=0.00;//Mask
 }
 
 return false;
}



bool AssignProbToNode(PTBL *tbl,POINTS *pt,GRAPH *g,MRC *mrc){
 float dcut = 3.0 / mrc->widthx;//3A cut-off
 float dcut2 = dcut*dcut;
 if(pt->Ncd != g->Nnode){
  printf("ERROR Ntbl= %d Ncd= %d Nnode= %d\n",tbl->N,pt->Ncd,g->Nnode);
  return true;
 }
 //search closest points

 #pragma omp parallel for schedule(dynamic,10)
 for(int i=0;i<pt->Ncd;i++){
  float dist;
  float MinDist=100;
  int MinID=-1;
  g->node[i].id=i;
 	for(int j=0;j<tbl->N;j++){
	 if(pt->cd[i][0] - tbl->prob[j][0] > dcut)
	  continue;
	 if(pt->cd[i][0] - tbl->prob[j][0] < -dcut)
	  continue;
	 if(pt->cd[i][1] - tbl->prob[j][1] > dcut)
	  continue;
	 if(pt->cd[i][1] - tbl->prob[j][1] < -dcut)
	  continue;
	 if(pt->cd[i][2] - tbl->prob[j][2] > dcut)
	  continue;
	 if(pt->cd[i][2] - tbl->prob[j][2] < -dcut)
	  continue;
	 dist = 
		 ( pt->cd[i][0]-tbl->prob[j][0] )*( pt->cd[i][0]-tbl->prob[j][0] )
		+( pt->cd[i][1]-tbl->prob[j][1] )*( pt->cd[i][1]-tbl->prob[j][1] )
		+( pt->cd[i][2]-tbl->prob[j][2] )*( pt->cd[i][2]-tbl->prob[j][2] );
	 if(dist > dcut2)
	  continue;
	 if(MinDist > dist){
	  MinDist=dist;
	  MinID=j;
	 }
	 //printf("%d dist= %f\n",i,sqrt(dist));
	}
	if(MinID==-1){//No Hits
	 g->node[i].SSP[0]=g->node[i].SSP[1]=g->node[i].SSP[2]=0.0;
	 g->node[i].ATOMP[0]=g->node[i].ATOMP[1]
		=g->node[i].ATOMP[2]=g->node[i].ATOMP[3]=0.0;
	 printf("##MISS %d %d dist= %f\n",i,MinID,sqrt(MinDist));
	 continue;
	}
	//AA
	for(int j=0;j<20;j++)//0-2:xyz, 3-22:AA type, 23-28: Atom, 29-31: SS
	 g->node[i].AAP[j]=tbl->prob[MinID][j+3];
	for(int j=0;j<6;j++)//23rd- Other,N,CA,C,O,CB
	 g->node[i].ATOMP[j]=tbl->prob[MinID][j+23];
	for(int j=0;j<3;j++)//H, E, Others
	 g->node[i].SSP[j]=tbl->prob[MinID][j+29];
	g->node[i].real_cd[0]=pt->cd[i][0]*mrc->widthx+mrc->orgxyz[0];
	g->node[i].real_cd[1]=pt->cd[i][1]*mrc->widthx+mrc->orgxyz[1];
	g->node[i].real_cd[2]=pt->cd[i][2]*mrc->widthx+mrc->orgxyz[2];
 }
 printf("N= %d Ncd= %d Nnode= %d\n",tbl->N,pt->Ncd,g->Nnode);

 //statics data ref prob
 double SumPaa[20]={};
 double SumPss[3]={};
 double SumPatm[6]={};
 double TotalPaa=0;
 double TotalPss=0;
 double TotalPatm=0;

 for(int i=0;i<pt->Ncd;i++){
  	for(int j=0;j<20;j++)
	 SumPaa[j]+=g->node[i].AAP[j];
	for(int j=0;j<6;j++)//23rd- Other,N,CA,C,O,CB
	 SumPatm[j]+=g->node[i].ATOMP[j];
	for(int j=0;j<3;j++)
	 SumPss[j]+=g->node[i].SSP[j];
 }
 for(int j=0;j<20;j++) TotalPaa+=SumPaa[j];
 for(int j=0;j<6;j++) TotalPatm+=SumPatm[j];
 for(int j=0;j<3;j++) TotalPss+=SumPss[j];

 //Reference P
 for(int j=0;j<20;j++) tbl->REFaa[j]=SumPaa[j]/TotalPaa;
 for(int j=0;j<6;j++) tbl->REFatm[j]=SumPatm[j]/TotalPatm;
 for(int j=0;j<3;j++) tbl->REFss[j]=SumPss[j]/TotalPss;

 //Show
 for(int j=0;j<20;j++) printf("#REFaa%d %f\n",j,tbl->REFaa[j]);
 for(int j=0;j<6;j++) printf("#REFatm%d %f\n",j,tbl->REFatm[j]);
 for(int j=0;j<3;j++) printf("#REFss%d %f\n",j,tbl->REFss[j]);

 //Log normalization
 for(int i=0;i<pt->Ncd;i++){
  	for(int j=0;j<20;j++)
	 g->node[i].LogAA[j]=log(g->node[i].AAP[j]/tbl->REFaa[j]);
	g->node[i].LogAA[20]=-100.00;//Mask
	for(int j=0;j<6;j++)//23rd- Other,N,CA,C,O,CB
	 g->node[i].LogATOM[j]=log(g->node[i].ATOMP[j]/tbl->REFatm[j]);
	for(int j=0;j<3;j++)
	 g->node[i].LogSS[j]=log(g->node[i].SSP[j]/tbl->REFss[j]);
 }

 return false;
}


bool AssignProbToNodePDB(PTBL *tbl,PDB *p,GRAPH *g,MRC *mrc){
 float dcut = 3.0 / mrc->widthx;//3A cut-off
 float dcut2 = dcut*dcut;
 puts("#Assign Probability to PDB ATOMs as Node data");
 //search closest points
 //malloc

 if((g->node=(NODE *)malloc(sizeof(NODE)*p->NumOfAtom))==NULL)
  return false;

 g->Nnode = p->NumOfAtom;
 //update tbl coordinates
 for(int j=0;j<tbl->N;j++){
  tbl->prob[j][0]=tbl->prob[j][0]*mrc->widthx+mrc->orgxyz[0];
  tbl->prob[j][1]=tbl->prob[j][1]*mrc->widthx+mrc->orgxyz[1];
  tbl->prob[j][2]=tbl->prob[j][2]*mrc->widthx+mrc->orgxyz[2];
 }
 unsigned int Nmiss_atm=0;
 #pragma omp parallel for schedule(dynamic,10)
 for(int i=0;i<p->NumOfAtom;i++){
  float dist;
  float MinDist=100;
  int MinID=-1;
  g->node[i].id=i;
 	for(int j=0;j<tbl->N;j++){
	 if(p->xyz[i][0] - tbl->prob[j][0] > dcut)
	  continue;
	 if(p->xyz[i][0] - tbl->prob[j][0] < -dcut)
	  continue;
	 if(p->xyz[i][1] - tbl->prob[j][1] > dcut)
	  continue;
	 if(p->xyz[i][1] - tbl->prob[j][1] < -dcut)
	  continue;
	 if(p->xyz[i][2] - tbl->prob[j][2] > dcut)
	  continue;
	 if(p->xyz[i][2] - tbl->prob[j][2] < -dcut)
	  continue;
	 dist = 
		 ( p->xyz[i][0]-tbl->prob[j][0] )*( p->xyz[i][0]-tbl->prob[j][0] )
		+( p->xyz[i][1]-tbl->prob[j][1] )*( p->xyz[i][1]-tbl->prob[j][1] )
		+( p->xyz[i][2]-tbl->prob[j][2] )*( p->xyz[i][2]-tbl->prob[j][2] );
	 if(dist > dcut2)
	  continue;
	 if(MinDist > dist){
	  MinDist=dist;
	  MinID=j;
	 }
	 //printf("%d dist= %f\n",i,sqrt(dist));
	}
	if(MinID==-1){//No Hits
	 g->node[i].SSP[0]=g->node[i].SSP[1]=g->node[i].SSP[2]=0.0;
 	 for(int j=0;j<6;j++)//23rd- Other,N,CA,C,O,CB
	  g->node[i].ATOMP[j]=0.0;
	 for(int j=0;j<20;j++)//0-2:xyz, 3-22:AA type, 23-28: Atom, 29-31: SS
	  g->node[i].AAP[j]=0.0;
	 //printf("##MISS %d %d dist= %f\n",i,MinID,sqrt(MinDist));
		#pragma omp atomic
	 	Nmiss_atm++;
	 continue;
	}
	//AA
	for(int j=0;j<20;j++)//0-2:xyz, 3-22:AA type, 23-28: Atom, 29-31: SS
	 g->node[i].AAP[j]=tbl->prob[MinID][j+3];
	for(int j=0;j<6;j++)//23rd- Other,N,CA,C,O,CB
	 g->node[i].ATOMP[j]=tbl->prob[MinID][j+23];
	for(int j=0;j<3;j++)//H, E, Others
	 g->node[i].SSP[j]=tbl->prob[MinID][j+29];
	g->node[i].real_cd[0]=p->xyz[i][0];
	g->node[i].real_cd[1]=p->xyz[i][1];
	g->node[i].real_cd[2]=p->xyz[i][2];
 }
 printf("N= %d Ncd= %d Nnode= %d\n",tbl->N,p->NumOfAtom,g->Nnode);
 printf("Number of miss-assigned atoms: %d\n", Nmiss_atm);

 //statics data ref prob
 double SumPaa[20]={};
 double SumPss[3]={};
 double SumPatm[6]={};
 double TotalPaa=0;
 double TotalPss=0;
 double TotalPatm=0;

 for(int i=0;i<p->NumOfAtom;i++){
  	for(int j=0;j<20;j++)
	 SumPaa[j]+=g->node[i].AAP[j];
	for(int j=0;j<6;j++)//23rd- Other,N,CA,C,O,CB
	 SumPatm[j]+=g->node[i].ATOMP[j];
	for(int j=0;j<3;j++)
	 SumPss[j]+=g->node[i].SSP[j];
 }
 for(int j=0;j<20;j++) TotalPaa+=SumPaa[j];
 for(int j=0;j<6;j++) TotalPatm+=SumPatm[j];
 for(int j=0;j<3;j++) TotalPss+=SumPss[j];

 //Reference P
 for(int j=0;j<20;j++) tbl->REFaa[j]=SumPaa[j]/TotalPaa;
 for(int j=0;j<6;j++) tbl->REFatm[j]=SumPatm[j]/TotalPatm;
 for(int j=0;j<3;j++) tbl->REFss[j]=SumPss[j]/TotalPss;

 //Show
 for(int j=0;j<20;j++) printf("#REFaa%d %f\n",j,tbl->REFaa[j]);
 for(int j=0;j<6;j++) printf("#REFatm%d %f\n",j,tbl->REFatm[j]);
 for(int j=0;j<3;j++) printf("#REFss%d %f\n",j,tbl->REFss[j]);

 //Log normalization
 for(int i=0;i<p->NumOfAtom;i++){
  	for(int j=0;j<20;j++)
	 if(g->node[i].AAP[j]>0)
	  g->node[i].LogAA[j]=log(g->node[i].AAP[j]/tbl->REFaa[j]);
	 else
  	  g->node[i].LogAA[j]=0.0;
	g->node[i].LogAA[20]=-100.00;//Masking
	for(int j=0;j<6;j++)//23rd- Other,N,CA,C,O,CB
	 if(g->node[i].ATOMP[j]>0)
	  g->node[i].LogATOM[j]=log(g->node[i].ATOMP[j]/tbl->REFatm[j]);
	 else
	  g->node[i].LogATOM[j]=0.0;

	for(int j=0;j<3;j++)
	 if(g->node[i].SSP[j]>0)
	  g->node[i].LogSS[j]=log(g->node[i].SSP[j]/tbl->REFss[j]);
	 else
	  g->node[i].LogSS[j]=0.0;
 }

 return false;
}




bool ShowPTBL(MRC *mrc,PTBL *tbl,float pcut){
 float dcut = 3.0 / mrc->widthx;
 float dcut2 = dcut*dcut;
 float tmp[3];

 int Natm=0;
 printf("#Back-bone Atoms P(N,CA,O) >= %.3f\n",pcut);
 puts("MODEL");
 for(int j=0;j<tbl->N;j++){

	if(tbl->prob[j][24]+tbl->prob[j][25]+tbl->prob[j][26] >=pcut){//main-chain N,CA,C
	 tmp[0]=tbl->prob[j][0]*mrc->widthx+mrc->orgxyz[0];
	 tmp[1]=tbl->prob[j][1]*mrc->widthx+mrc->orgxyz[1];
	 tmp[2]=tbl->prob[j][2]*mrc->widthx+mrc->orgxyz[2];;
  	 //printf("ATOM  %5d  CA  ALA%6d    ",Natm,Natm);
  	   printf("HETATM%5d  CA  CA %6d    ",Natm,Natm);
  	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,tbl->prob[j][25]);
	 Natm++;
	}
 }
 puts("ENDMDL");
/*
 printf("#OtherAtom P(Other) >= %.3f\n",pcut);
 puts("MODEL");
 Natm=0;
 for(int j=0;j<tbl->N;j++){

	if(tbl->prob[j][23]>=pcut){//Other
	 tmp[0]=tbl->prob[j][0]*mrc->widthx;
	 tmp[1]=tbl->prob[j][1]*mrc->widthx;
	 tmp[2]=tbl->prob[j][2]*mrc->widthx;
  	 printf("ATOM  %5d  N   ALA%6d    ",Natm,Natm);
  	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,tbl->prob[j][23]);
	 Natm++;
	}

 }
 puts("ENDMDL");
*/
 printf("#Helix P(H)> %.3f P(N,CA,O) >= %.3f\n",0.5,pcut);
 puts("MODEL");
 Natm=1;
 for(int j=0;j<tbl->N;j++){

	if(tbl->prob[j][29]>=0.5 && tbl->prob[j][24]+tbl->prob[j][25]+tbl->prob[j][26] >=pcut){//H
	 tmp[0]=tbl->prob[j][0]*mrc->widthx+mrc->orgxyz[0];
	 tmp[1]=tbl->prob[j][1]*mrc->widthx+mrc->orgxyz[1];
	 tmp[2]=tbl->prob[j][2]*mrc->widthx+mrc->orgxyz[2];
  	 //printf("ATOM  %5d  C   ALA%6d    ",Natm,Natm);
  	 printf("HETATM%5d  NA  NA %6d    ",Natm,Natm);
  	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,tbl->prob[j][29]);
	 Natm++;
	}

 }
 puts("ENDMDL");
 printf("#Sheet P(E)> %.3f P(N,CA,O) >= %.3f\n",0.5,pcut);
 puts("MODEL");
 Natm=0;
 for(int j=0;j<tbl->N;j++){

	if(tbl->prob[j][30]>=0.5 && tbl->prob[j][24]+tbl->prob[j][25]+tbl->prob[j][26] >=pcut){//E
	 tmp[0]=tbl->prob[j][0]*mrc->widthx+mrc->orgxyz[0];
	 tmp[1]=tbl->prob[j][1]*mrc->widthx+mrc->orgxyz[1];
	 tmp[2]=tbl->prob[j][2]*mrc->widthx+mrc->orgxyz[2];
  	 //printf("ATOM  %5d  O   ALA%6d    ",Natm,Natm);
  	 printf("HETATM%5d  MG  MG %6d    ",Natm,Natm);
  	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,tbl->prob[j][30]);
	 Natm++;
	}

 }
 puts("ENDMDL");
 

 return true;
}

bool check_same_nodes(NODE **a,NODE **b,int n){

 float d2;
 for(int i=0;i<n;i++){
  if(a[i] == b[i])//exact same node
   continue;
  if(a[i]==NULL && b[i]!=NULL) return false;
  if(a[i]!=NULL && b[i]==NULL) return false;

  d2=
	 (a[i]->real_cd[0]-b[i]->real_cd[0])*(a[i]->real_cd[0]-b[i]->real_cd[0])
	+(a[i]->real_cd[1]-b[i]->real_cd[1])*(a[i]->real_cd[1]-b[i]->real_cd[1])
	+(a[i]->real_cd[2]-b[i]->real_cd[2])*(a[i]->real_cd[2]-b[i]->real_cd[2]);

  if(d2<2.0)
   continue;


   return false;
 }
 return true;
}

bool ScoreFrag(FRAG_DB *db,int Ndb, SEQFG *sfg, int Nsfg,SEQ_NODE *res){
 unsigned int total=0;

 //malloc res
 for(int i=0;i<Nsfg;i++){
  if((res[i].score=(float*)malloc(sizeof(float)*Ndb))==NULL)
   return true;
  if((res[i].ali=(NODE ***)malloc(sizeof(NODE **)*Ndb))==NULL)
   return true;
	for(int j=0;j<Ndb;j++)
  	 if((res[i].ali[j]=(NODE **)malloc(sizeof(NODE *)*sfg[j].l))==NULL)
   	  return true;
 }


 printf("#Scoring....Nstr_frag= %d Nseq_frag= %d total= %d\n",Ndb,Nsfg,Ndb*Nsfg);
 #pragma omp parallel for schedule(dynamic,10)
 for(int s=0;s<Nsfg;s++){//Sequence
  float sco=0;
  float maxsco;
  int N=0;
  NODE *tmp_ali[100];
 	for(int i=0;i<Ndb;i++){
  	 if(db[i].N==0) continue;
	 //compute....
	 maxsco=0.00;
		for(int j=0;j<db[i].N;j++){
		 if(db[i].frag[j].Nnode==0)
		  continue;
		 //sco=QuickScoring(&(db[i].frag[j]),&sfg[s],res[s].ali[N]);
		 sco=QuickScoring(&(db[i].frag[j]),&sfg[s],tmp_ali);
		 if(sco>maxsco){
		  maxsco=sco;
			//copy
			for(int p=0;p<sfg[s].l;p++)
			 res[s].ali[N][p]=tmp_ali[p];
		  //printf("%d %d %d sco= %f\n",s,i,j,sco);
		  //copy NODE alignment
		 }
		}
	 if(maxsco>0.00){
	  //printf("#results seq %d db %d %f\n",s,i,maxsco);
	  res[s].score[N]=maxsco;
		/*
		for(int p=0;p<sfg[s].l;p++){
		 if(res[s].ali[N][p]==NULL) continue;
		 printf("[[%d]]",res[s].ali[N][p]->id);
		 //printf("[[%d]]",p);
		}
		*/
	  //puts("");
	  N++;
	 }
	}
	//printf("#seq %d N= %d/%d\n",s,N,Ndb);
  	//remove redundant models
	for(int i=0;i<N-1;i++){
		for(int j=i+1;j<N;j++){
		 //Same?
		 if(check_same_nodes(res[s].ali[i],res[s].ali[j],sfg[s].l)){
			//Keep better
			if(res[s].score[i]<res[s].score[j]){
			 for(int p=0;p<sfg[s].l;p++)
			  res[s].ali[i][p]=res[s].ali[j][p];
                         res[s].score[i]=res[s].score[j];
			}
		 	//shift all
			for(int k=j;k<N-1;k++){
			 //res[s].ali[p]=res[s].ali[p+1];
			 for(int p=0;p<sfg[s].l;p++)
			  res[s].ali[k][p]=res[s].ali[k+1][p];
			 res[s].score[k]=res[s].score[k+1];
			}
		  j--;
		  N--;
		 }
		}
	}
	printf("#AfterNR seq %d N= %d / %d\n",s,N,Ndb);
	res[s].Nali=N;
	res[s].sfg=&sfg[s];
	if(N==0)
	 continue;
	//check Z-score
	float ave=0;
	float std=0;
	for(int i=0;i<N;i++){
	 ave+=res[s].score[i];
	}
	ave=ave/(float)N;
	for(int i=0;i<N;i++){
	 std+=(res[s].score[i]-ave)*(res[s].score[i]-ave);
	}
	std=std/(float)N;
	for(int i=0;i<N;i++)
	 res[s].score[i]=(res[s].score[i]-ave)/std;
	printf("##STD= %f AVE= %f N= %d\n",ave,std,N);
 }
}



float QuickScoring(FRAGMENT *fg,SEQFG *s, NODE **out){
 //DP
 float Smtx[10000];//Similarity Matrix
 DPMTX Dmtx[10000];//DP Matrix
 float score=0;

 //Set Smtx str, seq
 int n1=fg->Nnode;
 int n2=s->l;
 float Waa,Wss,Watm;
 float Saa,Sss,Satm;
 Waa=1.0;Wss=1.0;Watm=1.0;

/*
 //Tabu check if nodes have P(N,CA,C) < 0.5
 int Ntabu=0;
 float SumP=0.0;
 for(int j=0;j<fg->Nnode;j++){
  SumP=fg->n_path_po[j]->ATOMP[1]+fg->n_path_po[j]->ATOMP[2]+fg->n_path_po[j]->ATOMP[3];
  if(SumP<0.5)
   Ntabu++;
  if(Ntabu>0)
   return 0.0;
 }
*/
 //Fill Smtx
 //Structure sequence
 for(int i=0;i<fg->Nnode;i++){
 	for(int j=0;j<s->l;j++){
	 //Score = log(P(AA,SS,ATOM)/Pref(AA,SS,ATOM))
	 Saa =	fg->n_path_po[i]->LogAA[s->seq[j]];
	 Satm=	fg->n_path_po[i]->LogATOM[2];//CA
	 Sss=	fg->n_path_po[i]->LogSS[0]*s->Pss[j][0]
		+fg->n_path_po[i]->LogSS[1]*s->Pss[j][1]
		+fg->n_path_po[i]->LogSS[2]*s->Pss[j][2];
	 //printf("Saa= %f\n",Saa);
	 Smtx[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	}
 }
 //DP
 int ali1[100],ali2[100],gali[200],Lgali;
 //score = dp_fast(Dmtx,Smtx,fg->dtbl,-100.0,-100.0,ali1,fg->Nnode,ali2,s->l,gali,&Lgali,false);
 //score = dp_fast(Dmtx,Smtx,fg->dtbl,-100.0,-100.0,ali1,fg->Nnode,ali2,s->l,gali,&Lgali,true);
 //score = dp_local(Dmtx,Smtx,fg->dtbl,-100.0,-100.0,ali1,fg->Nnode,ali2,s->l,gali,&Lgali,true);

 int Nm=0;
 for(int i=0;i<s->l;i++){
  int j=ali2[i];
  if(j!=-1){
   out[i]=fg->n_path_po[j];
   Nm++;
  }else{
   out[i]=NULL;
  }
 }

 //printf("Nm= %d score= %f\n",Nm,score);

 //ignore too short
 if(Nm<4)
  return 0.00;

/*
 puts("");
 for(int i=0;i<fg->Nnode;i++){
  int j=ali1[i];
  printf("[%d,%.1f],",ali1[i],fg->dtbl[i]);
 }

 puts("");
 for(int i=0;i<Lgali;i++){
  printf("[%d,%d]",gali[2*i],gali[2*i+1]);
 }
 puts("");

 for(int i=0;i<s->l;i++){
  float tmp[3];
  int j=ali2[i];
  //printf("%d,",ali1[i]);
  if(j==-1) continue;
  tmp[0]=fg->n_path_po[j]->real_cd[0];
  tmp[1]=fg->n_path_po[j]->real_cd[1];
  tmp[2]=fg->n_path_po[j]->real_cd[2];
  printf("ATOM  %5d  CA  %3s%6d    ",i+1,RES_NAMES[s->seq[i]],s->pos+i);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,fg->n_path_po[j]->SSP[0]);//

 }
 puts("ENDMDL");
*/
 //}

 return score;
}



bool TwoFragCheck(NODE **a,int l1, int p1, NODE **b,int l2,int p2){

 //printf("%d %d %d %d\n",l1,p1,l2,p2);
 int p3,p4;
 p3 = p1+l1 -1;
 p4 = p2+l2 -1;
 float d2;
 //far
 if(p1+l1 + 6  < p2){
  return true;
 }else if(p1+l1 -1 == p2 ){
  
  d2 =   (a[l1-1]->real_cd[0]-b[0]->real_cd[0])*(a[l1-1]->real_cd[0]-b[0]->real_cd[0]);
        +(a[l1-1]->real_cd[1]-b[0]->real_cd[1])*(a[l1-1]->real_cd[1]-b[0]->real_cd[1]);
        +(a[l1-1]->real_cd[2]-b[0]->real_cd[2])*(a[l1-1]->real_cd[2]-b[0]->real_cd[2]);

	if(d2 > 2.0 )
	 return true;

 }else{//many overlap
  return true;
 }

 return false;
}



bool CombinationSearch(SEQFG *sfg, int Nsfg,SEQ_NODE *res){

 unsigned int Ntot=0;
 //position
 for(int s1=0;s1<Nsfg;s1++){
  int N1=res[s1].Nali;
  if(N1==0)
   continue;
	for(int s2 = s1 + sfg[s1].l -1; s2 < Nsfg; s2++){
  	 int N2=res[s2].Nali;
  	 if(N2==0)
    	  continue;
		for(int i=0;i<N1;i++){
        	 if(res[s1].score[i]<1.0*sfg[s1].l)
        	  continue;
			for(int j=0;j<N2;j++){
			 if(TwoFragCheck(res[s1].ali[i],res[s1].sfg->l,res[s1].sfg->pos,res[s2].ali[j],res[s2].sfg->l,res[s2].sfg->pos)==false)
	 	 	  Ntot++;
	 		}
		}
	}
  printf("#Ntot= %d\n",Ntot);
 }
 printf("#Ntot= %d\n",Ntot);
}


void ShowThreadModel(GRAPH *g,SEQFG *sfg){
 //Show path
 //For density*volume, using K-clustring data
 int i,j,k,now,ind;
 int Natm=1;
 double tmp[3];
 int st,ed;
 double fmax=0;

 DP_MEMORY dpm;
 dpm.ali1=(int *)malloc(sizeof(int)*g->Nnode);
 dpm.ali2=(int *)malloc(sizeof(int)*sfg->l);
 dpm.gali=(int *)malloc(sizeof(int)*sfg->l*g->Nnode);
 dpm.Smtx=(float *)malloc(sizeof(float)*sfg->l*g->Nnode);
 dpm.SmtxRv=(float *)malloc(sizeof(float)*sfg->l*g->Nnode);
 dpm.dtbl=(float *)malloc(sizeof(float)*g->Nnode);
 dpm.Dmtx=(DPMTX *)malloc(sizeof(DPMTX)*sfg->l*g->Nnode);
 dpm.model=(NODE **)malloc(sizeof(NODE*)*sfg->l);
 dpm.path = (NODE **)malloc(sizeof(NODE*)*g->Nnode);
/*
 dpm.ali1=(int *)malloc(sizeof(int)*g->Nnode);
 dpm.ali2=(int *)malloc(sizeof(int)*sfg->l);
 dpm.gali=(int *)malloc(sizeof(int)*sfg->l*g->Nnode);
 dpm.Smtx=(float *)malloc(sizeof(float)*sfg->l*(g->Nnode+1));
 dpm.SmtxRv=(float *)malloc(sizeof(float)*sfg->l*(g->Nnode+1));
 dpm.dtbl=(float *)malloc(sizeof(float)*(g->Nnode+1));
 dpm.Dmtx=(DPMTX *)malloc(sizeof(DPMTX)*sfg->l*(g->Nnode+1));
 dpm.model=(NODE **)malloc(sizeof(NODE *)*sfg->l);
*/
 //path data
 
 //for(int mid=0;mid<n;mid++){

	float score=QualityPathDP(g,sfg,&dpm,true);
	printf("##SCORE= %f\n",score);
  	printf("MODEL\n");
 	Natm=1;
	for(int rnum=0;rnum<sfg->l;rnum++){
		if(dpm.ali2[rnum]==-1)
			continue;
	 //printf("poi= %d\n",dpm.ali2[rnum]);
	 tmp[0]=dpm.model[rnum]->real_cd[0];
	 tmp[1]=dpm.model[rnum]->real_cd[1];
	 tmp[2]=dpm.model[rnum]->real_cd[2];
 	 printf("ATOM  %5d  CA  %3s%6d    ",Natm
			 ,RES_NAMES[sfg->seq[rnum]]
			 ,rnum+1);
 	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,1.0);
 	 Natm++;
	}
  	printf("TER\nENDMDL\n");

 //}
}



int Pos2Cd(MRC *m,int pos[3],double cd[3]){
        cd[0]=m->orgxyz[0]+pos[0]*m->widthx;
        cd[1]=m->orgxyz[1]+pos[1]*m->widthy;
        cd[2]=m->orgxyz[2]+pos[2]*m->widthz;
}

void VoxPos(MRC *m,double cd[3],int pos[3]){

        pos[0]=(int)((cd[0]-m->orgxyz[0])/m->widthx);
        pos[1]=(int)((cd[1]-m->orgxyz[1])/m->widthy);
        pos[2]=(int)((cd[2]-m->orgxyz[2])/m->widthz);

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

bool MaskByPdb(MRC *m,PDB *p,float r){

 int pos1[3],pos2[3],pos3[3],idx;
 double cd1[3],cd2[3],cd3[3];
 float r2=r*r;
 float d2;
 
 for(int i=0;i<p->NumOfAtom;i++){
  cd1[0]=p->xyz[i][0]-r;
  cd1[1]=p->xyz[i][1]-r;
  cd1[2]=p->xyz[i][2]-r;
  cd2[0]=p->xyz[i][0]+r;
  cd2[1]=p->xyz[i][1]+r;
  cd2[2]=p->xyz[i][2]+r;
  VoxPos(m,cd1,pos1);
  VoxPos(m,cd2,pos2);
         //fill
                for(pos3[0]=pos1[0];pos3[0]<=pos2[0];pos3[0]++){
                for(pos3[1]=pos1[1];pos3[1]<=pos2[1];pos3[1]++){
                for(pos3[2]=pos1[2];pos3[2]<=pos2[2];pos3[2]++){
                        idx=Pos2Idx(m,pos3);
                        if(idx==-1)
                                continue;
                        Pos2Cd(m,pos3,cd3);
                        d2=      (p->xyz[i][0]-cd3[0])*(p->xyz[i][0]-cd3[0])
                                +(p->xyz[i][1]-cd3[1])*(p->xyz[i][1]-cd3[1])
                                +(p->xyz[i][2]-cd3[2])*(p->xyz[i][2]-cd3[2]);
                        //Update
                        if(d2 < r2 ){
                         //printf("Update: %d %d %f\n",idx,i,d2);
			 m->dens[idx]=0.00;//Masking
                        }
                }}}

        
 }
 return true;
}



//Consider Mask
bool DFS_Tree(GRAPH *g,TREE *t){
 int i,j;
 int BestStart;
 double BestSco=0;
 for(int i=0;i<t->Nnode;i++){
  t->UsedN[i]=false;
  //if(t->MaskN[i]) printf("Maked %d\n",i);
 }

 //DFS
 for(int start=0;start<t->Nnode;start++){
  //printf("start %d\n",start);
	if(t->UsedN[start])
	 continue;
 	//int st=t->St;
 	int st=start;
 	int Nst=1;
 	int v,w,maxi,maxi2;
 	double maxd=0;
 	double sco=0;
 	int pre_st=-11;
 	while(1){//Iter
  	 t->stock[0]=st;
  	 Nst=1;
  	 t->nextv[st]=-1;
	 	//init
  	 	for(i=0;i<t->Nnode;i++){
  	 	 t->cost[i]=0.00;
  	 	 t->ActN[i]=false;
  	 	}
  	 maxd=0;
  	 maxi=-1;
		while(1){//START DFS
	 	 if(Nst==0)
	 	  break;
  	 	 v=t->stock[Nst-1];//pop
	 	 Nst--;
	 	 if(t->ActN[v])
	 	  continue;
		 if(t->MaskN[v]==true){
		  //printf("Masked %d\n",v);
		  continue;
		 }
	 	 t->ActN[v]=true;
	 	 t->UsedN[v]=true;
	 	 //branch
	 		for(j=0;j<g->node[v].N;j++){
	 		 int eid=g->node[v].e[j]->eid;
	  		 if(t->ActE[eid]==false) continue;

	  		 w=g->node[v].e[j]->id1;
	  		 if(t->ActN[w]) continue;
	  		 if(t->MaskN[w]) continue;
	  		 t->stock[Nst]=w;
	  		 Nst++;

	  		 t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
	  		 t->nextv[w]=v;//path

	  			if(maxd<t->cost[w]){
	  			 maxd=t->cost[w];
	  			 maxi=w;
	  			}
	 		}
	 		for(j=0;j<g->node[v].N;j++){
 	 		 int eid=g->node[v].e[j]->eid;
	  		 if(t->ActE[eid]==false) continue;

	  		 w=g->node[v].e[j]->id2;
	  		 if(t->ActN[w]) continue;
	  		 if(t->MaskN[w]) continue;
	  		 //push
	  		 t->stock[Nst]=w;
	  	 	 Nst++;

	  	 	 t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
	  	 	 t->nextv[w]=v;//path

	  	 	 	if(maxd<t->cost[w]){
	  	 		 maxd=t->cost[w];
	  	 		 maxi=w;
	  	 		}
	 		}
		}//END DFS
  	 if(maxi==pre_st)
  	  break;
  	 
  	 pre_st=st;
  	 st=maxi;
  	}//END Iter
	if(maxi==-1)
   	 continue;
	//printf("start= %d %d\n",st,start);
	if(maxd>BestSco){
	 BestSco=maxd;
	 BestStart=maxi;
	}
  if(BestStart==-1)
   return 0;
 }

 //Again....
 int st=BestStart;
 int Nst=1;
 int v,w,maxi,maxi2;
 double maxd=0;
 double sco=0;
 int pre_st=-11;
 while(1){//Iter
  	 t->stock[0]=st;
  	 Nst=1;
  	 t->nextv[st]=-1;
	 	//init
  	 	for(i=0;i<t->Nnode;i++){
  	 	 t->cost[i]=0.00;
  	 	 t->ActN[i]=false;
  	 	}
  	 maxd=0;
  	 maxi=-1;
		while(1){//START DFS
	 	 if(Nst==0)
	 	  break;
  	 	 v=t->stock[Nst-1];//pop
	 	 Nst--;
	 	 if(t->ActN[v])
	 	  continue;
		 if(t->MaskN[v]){
		  continue;
		 }
	 	 t->ActN[v]=true;
	 	 t->UsedN[v]=true;
	 	 //branch
	 		for(j=0;j<g->node[v].N;j++){
	 		 int eid=g->node[v].e[j]->eid;
	  		 if(t->ActE[eid]==false) continue;

	  		 w=g->node[v].e[j]->id1;
	  		 if(t->ActN[w]) continue;
	  		 if(t->MaskN[w]) continue;
	  		 t->stock[Nst]=w;
	  		 Nst++;

	  		 t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
	  		 t->nextv[w]=v;//path

	  			if(maxd<t->cost[w]){
	  			 maxd=t->cost[w];
	  			 maxi=w;
	  			}
	 		}
	 		for(j=0;j<g->node[v].N;j++){
 	 		 int eid=g->node[v].e[j]->eid;
	  		 if(t->ActE[eid]==false) continue;

	  		 w=g->node[v].e[j]->id2;
	  		 if(t->ActN[w]) continue;
	  		 if(t->MaskN[w]) continue;
	  		 //push
	  		 t->stock[Nst]=w;
	  	 	 Nst++;

	  	 	 t->cost[w]=t->cost[v]+g->node[v].e[j]->dens;
	  	 	 t->nextv[w]=v;//path

	  	 	 	if(maxd<t->cost[w]){
	  	 		 maxd=t->cost[w];
	  	 		 maxi=w;
	  	 		}
	 		}
		}//END DFS
  	 if(maxi==pre_st)
  	  break;
  	 
  	 pre_st=st;
  	 st=maxi;
  }//END Iter

 //Set Start&End
 t->St=maxi;
 t->Ed=st;//!!! Aug-24

 sco=maxd*maxd;

 //Set 1st path
 for(i=0;i<t->Nnode;i++)
  t->ActN[i]=false;

 t->ActN[maxi]=true;
 int now=maxi;
 t->Path[0]=now;
 t->Lpath=1;
 while(1){
  now=t->nextv[now];
  t->ActN[now]=true;
  t->Path[t->Lpath]=now;
  t->Lpath++;
  if(now==st)
   break;
 }

}




double QualityTreeDP(GRAPH *g,TREE *t,SEQFG *sfg, DP_MEMORY *dpm, bool flag){
 int i,j;

 for(i=0;i<t->Nnode;i++)
  t->MaskN[i]=false;
 //DFS search
 DFS_Tree(g,t);
 NODE **n_path_po=dpm->path;
 float *Smtx=dpm->Smtx;
 float *SmtxRv=dpm->SmtxRv;//Reverse sequence
 float *dtbl=dpm->dtbl;
 float *dtblR=dpm->dtblR;
 DPMTX *Dmtx=dpm->Dmtx;
 float Saa, Satm,Sss;
 float Waa, Watm,Wss;
 int n1=t->Lpath;
 int n2=sfg->l;

 //Waa=Watm=Wss=1.00;
 Wss=cmd.Wss;
 Waa=cmd.Waa;
 Watm=cmd.Watm;
 float ScoCut = cmd.LowestSco;
 //printf("Npath= %d\n",t->Lpath);
 for(int i=0;i<t->Lpath;i++){
  n_path_po[i]=&(g->node[t->Path[i]]);
 }
 dtbl[0]=0;

 //CA-CA distance score
/* v02
 printf("Npath= %d\n",t->Lpath);
 for(int i=0;i<t->Lpath-1;i++){
  dtbl[i+1]=dtbl[i]+sqrt(
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])*
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])+
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])*
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])+
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2])*
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2]));
 }
*/
 //puts("##Set DTBL...");
 //Ver03 New dtbl position,position-j
 for(int i=0;i<t->Lpath;i++){
	for(int j=0;j<21;j++)
	 dtbl[i*21+j]=10.00;//init
	for(int j=1;j<21 && i-j>-1;j++){
  	 dtbl[i*21+j]=sqrt(
	  (n_path_po[i]->real_cd[0]-n_path_po[i-j]->real_cd[0])*
	  (n_path_po[i]->real_cd[0]-n_path_po[i-j]->real_cd[0])+
	  (n_path_po[i]->real_cd[1]-n_path_po[i-j]->real_cd[1])*
	  (n_path_po[i]->real_cd[1]-n_path_po[i-j]->real_cd[1])+
	  (n_path_po[i]->real_cd[2]-n_path_po[i-j]->real_cd[2])*
	  (n_path_po[i]->real_cd[2]-n_path_po[i-j]->real_cd[2]));
  	}
 }
 for(int i=0;i<t->Lpath;i++){
  int I = t->Lpath -1 - i;
	for(int j=0;j<21;j++)
	 dtblR[i*21+j]=10.00;//init
	for(int j=1;j<21 && i-j>-1;j++){
  	 dtblR[i*21+j]=sqrt(
	  (n_path_po[I]->real_cd[0]-n_path_po[I+j]->real_cd[0])*
	  (n_path_po[I]->real_cd[0]-n_path_po[I+j]->real_cd[0])+
	  (n_path_po[I]->real_cd[1]-n_path_po[I+j]->real_cd[1])*
	  (n_path_po[I]->real_cd[1]-n_path_po[I+j]->real_cd[1])+
	  (n_path_po[I]->real_cd[2]-n_path_po[I+j]->real_cd[2])*
	  (n_path_po[I]->real_cd[2]-n_path_po[I+j]->real_cd[2]));
  	}
 }

 //Set Smtx
 for(int i=0;i<t->Lpath;i++){
        for(int j=0;j<sfg->l;j++){
         Saa =  n_path_po[i]->LogAA[sfg->seq[j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[j][2];
         Smtx[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	 //printf("%f, %f\n",n_path_po[i]->LogAA[20],n_path_po[i]->LogAA[19]);
	 //printf("[%d %d] Aa %d %f %f %f %f\n",i,j,sfg->seq[j],Saa,Satm,Sss,Smtx[GRID2D(n1,n2,i,j)]);
        }
 }
 //Smtx for reversed sequence
 for(int i=0;i<t->Lpath;i++){
        for(int j=0;j<sfg->l;j++){

         Saa =  n_path_po[i]->LogAA[sfg->seq[sfg->l-1-j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[sfg->l-1-j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[sfg->l-1-j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[sfg->l-1-j][2];
         SmtxRv[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	 //printf("[%d %d] Aa %d %f %f %f %f\n",i,j,sfg->seq[j],Saa,Satm,Sss,Smtx[GRID2D(n1,n2,i,j)]);
        }
 }

 //Iter-DP
 float score,score2,final_score=0.00;
 int ali2[2000],ali2rv[2000];
 //init model
 for(int i=0;i<sfg->l;i++)
  dpm->ali2[i]=-1;
 
 for(int i=0;i<t->Nnode;i++)
  t->MaskN[i]=false;
  	

 for(int iter=0;iter < 20; iter++){
  //printf("Iter %d\n",iter);
 	score =  dp_local(Dmtx,Smtx,  dtbl,-100.0,-100.0,dpm->ali1,t->Lpath,ali2,  sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);
 	score2 = dp_local(Dmtx,SmtxRv,dtblR,-100.0,-100.0,dpm->ali1,t->Lpath,ali2rv,sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);

	//printf("%d %f %f\n",iter,score,score2);

	if(score < ScoCut && score2 < ScoCut)
	 break;

 	if(score2<score){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){
		 int j = ali2[pos];
			if(j!=-1){
			 //printf("[%d %d]",j,n_path_po[j]->id);
			 //mask Node
			 t->MaskN[n_path_po[j]->id]=true;
		 	 dpm->ali2[pos]=j;
		 	 dpm->model[pos]=n_path_po[j];
			}
		}
	 final_score += score;
 	}else{
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){
		 int j = ali2rv[pos];//Revese
			if(j!=-1){
			 //mask
			 t->MaskN[n_path_po[j]->id]=true;
		 	 dpm->ali2[sfg->l-1-pos]=j;
		 	 dpm->model[sfg->l-1-pos]=n_path_po[j];
			}
		}
		//puts("");
	 final_score += score2;
	}

	//Find longest path
	DFS_Tree(g,t);
	//printf("#Iter %d L= %d\n",iter,t->Lpath);
	for(int i=0;i<t->Lpath;i++)
  	 n_path_po[i]=&(g->node[t->Path[i]]);
/* //v02
	for(int i=0;i<t->Lpath-1;i++){
  	 dtbl[i+1]=dtbl[i]+sqrt(
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])*
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])+
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])*
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])+
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2])*
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2]));
 	}
*/
	 //Ver03 New dtbl position,position-j
	 for(int i=0;i<t->Lpath;i++){
		for(int j=0;j<21;j++)
		 dtbl[i*21+j]=10.00;//init
		for(int j=1;j<21 && i-j>-1;j++){
		 dtbl[i*21+j]=sqrt(
		  (n_path_po[i]->real_cd[0]-n_path_po[i-j]->real_cd[0])*
		  (n_path_po[i]->real_cd[0]-n_path_po[i-j]->real_cd[0])+
		  (n_path_po[i]->real_cd[1]-n_path_po[i-j]->real_cd[1])*
		  (n_path_po[i]->real_cd[1]-n_path_po[i-j]->real_cd[1])+
		  (n_path_po[i]->real_cd[2]-n_path_po[i-j]->real_cd[2])*
		  (n_path_po[i]->real_cd[2]-n_path_po[i-j]->real_cd[2]));
		}
	 }
	 for(int i=0;i<t->Lpath;i++){
	  int I = t->Lpath -1 - i;
		for(int j=0;j<21;j++)
		 dtblR[i*21+j]=10.00;//init
		for(int j=1;j<21 && i-j>-1;j++){
		 dtblR[i*21+j]=sqrt(
		  (n_path_po[I]->real_cd[0]-n_path_po[I+j]->real_cd[0])*
		  (n_path_po[I]->real_cd[0]-n_path_po[I+j]->real_cd[0])+
		  (n_path_po[I]->real_cd[1]-n_path_po[I+j]->real_cd[1])*
		  (n_path_po[I]->real_cd[1]-n_path_po[I+j]->real_cd[1])+
		  (n_path_po[I]->real_cd[2]-n_path_po[I+j]->real_cd[2])*
		  (n_path_po[I]->real_cd[2]-n_path_po[I+j]->real_cd[2]));
		}
	 }


	//Fill Smtx and SmtxRv
	for(int i=0;i<t->Lpath;i++){
        for(int j=0;j<sfg->l;j++){
         Saa =  n_path_po[i]->LogAA[sfg->seq[j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[j][2];
         Smtx[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	 //Mask
	 if(dpm->ali2[j]!=-1 || t->MaskN[t->Path[i]]==true)
	  Smtx[GRID2D(n1,n2,i,j)]=-1.0;
        } }
 	for(int i=0;i<t->Lpath;i++){
        for(int j=0;j<sfg->l;j++){
         Saa =  n_path_po[i]->LogAA[sfg->seq[sfg->l-1-j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[sfg->l-1-j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[sfg->l-1-j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[sfg->l-1-j][2];
         SmtxRv[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	 if(dpm->ali2[sfg->l-1-j]!=-1||t->MaskN[t->Path[i]]==true)
	  SmtxRv[GRID2D(n1,n2,i,j)]=-1.0;
        } }
 }
 t->score=final_score;
 return final_score;
}


double QualityTreeDPSubOpt(GRAPH *g,TREE *t,SEQFG *sfg, DP_MEMORY *dpm, SUBOPT_ALI *sub){
 int i,j;

 for(i=0;i<t->Nnode;i++)
  t->MaskN[i]=false;
 //DFS search
 DFS_Tree(g,t);
 NODE **n_path_po=dpm->path;
 float *Smtx=dpm->Smtx;
 float *SmtxRv=dpm->SmtxRv;//Reverse sequence

 float *Smtx_ori=dpm->Smtx_ori;
 float *SmtxRv_ori=dpm->SmtxRv_ori;

 float *dtbl=dpm->dtbl;
 float *dtblR=dpm->dtblR;
 DPMTX *Dmtx=dpm->Dmtx;
 float Saa, Satm,Sss;
 float Waa, Watm,Wss;
 int n1=t->Lpath;
 int n2=sfg->l;

 //Waa=Watm=Wss=1.00;
 Wss=cmd.Wss;
 Waa=cmd.Waa;
 Watm=cmd.Watm;
 float ScoCut = cmd.LowestSco;

 //printf("Npath= %d\n",t->Lpath);
 for(int i=0;i<t->Lpath;i++){
  n_path_po[i]=&(g->node[t->Path[i]]);
 }
 dtbl[0]=0;
 //CA-CA distance score
/* //v02
 for(int i=0;i<t->Lpath-1;i++){
  dtbl[i+1]=dtbl[i]+sqrt(
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])*
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])+
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])*
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])+
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2])*
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2]));
 }
*/
 //Ver03 New dtbl position,position-j
 for(int i=0;i<t->Lpath;i++){
	for(int j=0;j<21;j++)
	 dtbl[i*21+j]=10.00;//init
	for(int j=1;j<21 && i-j>-1;j++){
  	 dtbl[i*21+j]=sqrt(
	  (n_path_po[i]->real_cd[0]-n_path_po[i-j]->real_cd[0])*
	  (n_path_po[i]->real_cd[0]-n_path_po[i-j]->real_cd[0])+
	  (n_path_po[i]->real_cd[1]-n_path_po[i-j]->real_cd[1])*
	  (n_path_po[i]->real_cd[1]-n_path_po[i-j]->real_cd[1])+
	  (n_path_po[i]->real_cd[2]-n_path_po[i-j]->real_cd[2])*
	  (n_path_po[i]->real_cd[2]-n_path_po[i-j]->real_cd[2]));
  	}
 }
 for(int i=0;i<t->Lpath;i++){
  int I = t->Lpath -1 - i;
	for(int j=0;j<21;j++)
	 dtblR[i*21+j]=10.00;//init
	for(int j=1;j<21 && i-j>-1;j++){
  	 dtblR[i*21+j]=sqrt(
	  (n_path_po[I]->real_cd[0]-n_path_po[I+j]->real_cd[0])*
	  (n_path_po[I]->real_cd[0]-n_path_po[I+j]->real_cd[0])+
	  (n_path_po[I]->real_cd[1]-n_path_po[I+j]->real_cd[1])*
	  (n_path_po[I]->real_cd[1]-n_path_po[I+j]->real_cd[1])+
	  (n_path_po[I]->real_cd[2]-n_path_po[I+j]->real_cd[2])*
	  (n_path_po[I]->real_cd[2]-n_path_po[I+j]->real_cd[2]));
  	}
 }


 //Set Smtx
 for(int i=0;i<t->Lpath;i++){
        for(int j=0;j<sfg->l;j++){
         Saa =  n_path_po[i]->LogAA[sfg->seq[j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[j][2];
         Smtx[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	 Smtx_ori[GRID2D(n1,n2,i,j)] = Smtx[GRID2D(n1,n2,i,j)];
	 //printf("%f, %f\n",n_path_po[i]->LogAA[20],n_path_po[i]->LogAA[19]);
	 //printf("[%d %d] Aa %d %f %f %f %f\n",i,j,sfg->seq[j],Saa,Satm,Sss,Smtx[GRID2D(n1,n2,i,j)]);
        }
 }
 //Smtx for reversed sequence
 for(int i=0;i<t->Lpath;i++){
        for(int j=0;j<sfg->l;j++){

         Saa =  n_path_po[i]->LogAA[sfg->seq[sfg->l-1-j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[sfg->l-1-j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[sfg->l-1-j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[sfg->l-1-j][2];
         SmtxRv[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	 SmtxRv_ori[GRID2D(n1,n2,i,j)] = SmtxRv[GRID2D(n1,n2,i,j)];
	 //printf("[%d %d] Aa %d %f %f %f %f\n",i,j,sfg->seq[j],Saa,Satm,Sss,Smtx[GRID2D(n1,n2,i,j)]);
        }
 }

 //Iter-DP
 float score,score2,final_score=0.00;
 int ali2[2000],ali2rv[2000];
 bool MaskSeq[2000];
 int Nmodel=0;
 //init model
 for(int i=0;i<sfg->l;i++){
  dpm->ali2[i]=-1;
  MaskSeq[i]=false;
 }
 
 sub[Nmodel].L=0;

 for(int i=0;i<t->Nnode;i++)
  t->MaskN[i]=false;
  	

 for(int iter=0;iter < 20; iter++){

 	score =  dp_local(Dmtx,Smtx,  dtbl,-100.0,-100.0,dpm->ali1,t->Lpath,ali2,  sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);
 	score2 = dp_local(Dmtx,SmtxRv,dtblR,-100.0,-100.0,dpm->ali1,t->Lpath,ali2rv,sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);

	//printf("%d %f %f\n",iter,score,score2);

	if(score < ScoCut && score2 < ScoCut)
	 break;

 	if(score>ScoCut){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){//sequence
		 int j = ali2[pos];
			if(j!=-1){
			 //printf("[%d %d]",j,n_path_po[j]->id);
			 //mask Node
			 if(score > score2){
			  t->MaskN[n_path_po[j]->id]=true;
			  MaskSeq[pos]=true;
		 	  dpm->ali2[pos]=j;
		 	  dpm->model[pos]=n_path_po[j];
			 }
			 //input to subopt
			 sub[Nmodel].model[sub[Nmodel].L] = n_path_po[j];
			 sub[Nmodel].seq[sub[Nmodel].L] = pos;
			 //printf("ALI pos %d -> ldp %d\n",pos,j);
			 //for(int p=0;p<20;p++) printf("Dtbl %f\n",dtbl[j*21+p]);
			 //for(int p=0;p<20;p++) printf("DtblR %f\n",dtblR[j*21+p]);
			 
			 sub[Nmodel].L++;
			}
		}
		sub[Nmodel].sco = score;
		Nmodel++;
		sub[Nmodel].L=0;
		//printf("sco= %f\n",score);
		if(score>score2)
		 final_score += score;
 	}
	if(score2>ScoCut){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){
		 int j = ali2rv[pos];//Revese
			if(j!=-1){
			 //mask
			 if(score2>score){
			  t->MaskN[n_path_po[j]->id]=true;
			  MaskSeq[sfg->l-1-pos] = true;//actual sequence order
		 	  dpm->ali2[sfg->l-1-pos]=j;
		 	  dpm->model[sfg->l-1-pos]=n_path_po[j];
			 }
			 //node and seq
			 sub[Nmodel].model[sub[Nmodel].L] = n_path_po[j];
			 sub[Nmodel].seq[sub[Nmodel].L] = sfg->l-1-pos;
			 sub[Nmodel].L++;
			}
		}
		//reverse
		for(int pos=0;pos<sub[Nmodel].L-pos-1;pos++){
		 NODE *node=sub[Nmodel].model[pos];
		 int tmp_seq=sub[Nmodel].seq[pos];
		 sub[Nmodel].model[pos] = sub[Nmodel].model[sub[Nmodel].L-pos-1];
		 sub[Nmodel].seq[pos] = sub[Nmodel].seq[sub[Nmodel].L-pos-1];
		 sub[Nmodel].model[sub[Nmodel].L-pos-1]=node;
		 sub[Nmodel].seq[sub[Nmodel].L-pos-1]=tmp_seq;
		}
		
		//puts("");
		sub[Nmodel].sco = score2;
		Nmodel++;
		sub[Nmodel].L=0;
		if(score2>score) 
		 final_score += score2;
	}

	//Sub-optimal Alignment NEW Mar 15 2021
	//Mask aligned [Node:Seq] Pair SUBOPT1---------------
	for(int i=0;i<t->Lpath;i++){
	 if(t->MaskN[t->Path[i]]==false)
	  continue;
	  for(int j=0;j<sfg->l;j++){
	   if(MaskSeq[j]== false)
	    continue;
	 //Mask
	   Smtx[GRID2D(n1,n2,i,j)]=-1.0;
	 	if(j!=0)
	    	 Smtx[GRID2D(n1,n2,i,j-1)]=-1.0;
		if(j!=sfg->l-1)
	    	 Smtx[GRID2D(n1,n2,i,j+1)]=-1.0;
	   if(i!=0){
	    Smtx[GRID2D(n1,n2,i-1,j  )]=-1.0;
		if(j!=0)
	    	 Smtx[GRID2D(n1,n2,i-1,j-1)]=-1.0;
		if(j!=sfg->l-1)
	    	 Smtx[GRID2D(n1,n2,i-1,j+1)]=-1.0;
	   }
	   if(i!=t->Lpath-1){
	    Smtx[GRID2D(n1,n2,i+1,j  )]=-1.0;
		if(j!=0)
	    	 Smtx[GRID2D(n1,n2,i+1,j-1)]=-1.0;
		if(j!=sfg->l-1)
	    	 Smtx[GRID2D(n1,n2,i+1,j+1)]=-1.0;
	   }
	}}
	
	for(int i=0;i<t->Lpath;i++){
	 if(t->MaskN[t->Path[i]]==false)
	  continue;
	  for(int j=0;j<sfg->l;j++){
	   if(MaskSeq[sfg->l-1-j]== false)
	    continue;
	 //Mask
	   SmtxRv[GRID2D(n1,n2,i,j)]=-1.0;
	 	if(j!=0) SmtxRv[GRID2D(n1,n2,i,j-1)]=-1.0;
		if(j!=sfg->l-1) SmtxRv[GRID2D(n1,n2,i,j+1)]=-1.0;
	   if(i!=0){
	    SmtxRv[GRID2D(n1,n2,i-1,j  )]=-1.0;
		if(j!=0) SmtxRv[GRID2D(n1,n2,i-1,j-1)]=-1.0;
		if(j!=sfg->l-1) SmtxRv[GRID2D(n1,n2,i-1,j+1)]=-1.0;
	   }
	   if(i!=t->Lpath-1){
	    SmtxRv[GRID2D(n1,n2,i+1,j  )]=-1.0;
		if(j!=0) SmtxRv[GRID2D(n1,n2,i+1,j-1)]=-1.0;
		if(j!=sfg->l-1) SmtxRv[GRID2D(n1,n2,i+1,j+1)]=-1.0;
	   }
	}}
	//DP SUBOPT1
	score =  dp_local(Dmtx,Smtx,  dtbl,-100.0,-100.0,dpm->ali1,t->Lpath,ali2,  sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);
	score2 = dp_local(Dmtx,SmtxRv,dtblR,-100.0,-100.0,dpm->ali1,t->Lpath,ali2rv,sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);
	if(score > ScoCut){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){//sequence
		 int j = ali2[pos];
			if(j!=-1){
			 //input to subopt
			 sub[Nmodel].model[sub[Nmodel].L] = n_path_po[j];
			 sub[Nmodel].seq[sub[Nmodel].L] = pos;
			 sub[Nmodel].L++;
			}
		}
		sub[Nmodel].sco = score; Nmodel++;
		sub[Nmodel].L=0;
 	}
	if(score2 > ScoCut){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){
		 int j = ali2rv[pos];//Revese
			if(j!=-1){
			 //node and seq
			 sub[Nmodel].model[sub[Nmodel].L] = n_path_po[j];
			 sub[Nmodel].seq[sub[Nmodel].L] = sfg->l-1-pos;
			 sub[Nmodel].L++;
			}
		}
		//reverse
		for(int pos=0;pos<sub[Nmodel].L-pos-1;pos++){
		 NODE *node=sub[Nmodel].model[pos];
		 int tmp_seq=sub[Nmodel].seq[pos];
		 sub[Nmodel].model[pos] = sub[Nmodel].model[sub[Nmodel].L-pos-1];
		 sub[Nmodel].seq[pos] = sub[Nmodel].seq[sub[Nmodel].L-pos-1];
		 sub[Nmodel].model[sub[Nmodel].L-pos-1]=node;
		 sub[Nmodel].seq[sub[Nmodel].L-pos-1]=tmp_seq;
		}
		sub[Nmodel].sco = score2;
		Nmodel++;
		sub[Nmodel].L=0;
	}
	//END SUBOPT1------------

	//Mask aligned Node **Only** SUBOPT2
	for(int i=0;i<t->Lpath;i++){
	 if(t->MaskN[t->Path[i]]==false)
	  continue;
	for(int j=0;j<sfg->l;j++){
	 //Mask Node
	 Smtx[GRID2D(n1,n2,i,j)]=-1.0;
	 SmtxRv[GRID2D(n1,n2,i,j)]=-1.0;
	 if(i>0){//Mask only node -1
	  Smtx[GRID2D(n1,n2,i-1,j)]=-1.0;
	  SmtxRv[GRID2D(n1,n2,i-1,j)]=-1.0;
	 }
	 if(i<t->Lpath-1){//Mask only node +1
	  Smtx[GRID2D(n1,n2,i+1,j)]=-1.0;
	  SmtxRv[GRID2D(n1,n2,i+1,j)]=-1.0;
	 }
	}}
	score =  dp_local(Dmtx,Smtx,  dtbl,-100.0,-100.0,dpm->ali1,t->Lpath,ali2,  sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);
 	score2 = dp_local(Dmtx,SmtxRv,dtblR,-100.0,-100.0,dpm->ali1,t->Lpath,ali2rv,sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);
	
	if(score > ScoCut){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){//sequence
		 int j = ali2[pos];
			if(j!=-1){
			 //input to subopt
			 sub[Nmodel].model[sub[Nmodel].L] = n_path_po[j];
			 sub[Nmodel].seq[sub[Nmodel].L] = pos;
			 sub[Nmodel].L++;
			}
		}
		sub[Nmodel].sco = score;
		Nmodel++;
		sub[Nmodel].L=0;
 	}
	if(score2 > ScoCut){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){
		 int j = ali2rv[pos];//Revese
			if(j!=-1){
			 //node and seq
			 sub[Nmodel].model[sub[Nmodel].L] = n_path_po[j];
			 sub[Nmodel].seq[sub[Nmodel].L] = sfg->l-1-pos;
			 sub[Nmodel].L++;
			}
		}
		//reverse
		for(int pos=0;pos<sub[Nmodel].L-pos-1;pos++){
		 NODE *node=sub[Nmodel].model[pos];
		 int tmp_seq=sub[Nmodel].seq[pos];
		 sub[Nmodel].model[pos] = sub[Nmodel].model[sub[Nmodel].L-pos-1];
		 sub[Nmodel].seq[pos] = sub[Nmodel].seq[sub[Nmodel].L-pos-1];
		 sub[Nmodel].model[sub[Nmodel].L-pos-1]=node;
		 sub[Nmodel].seq[sub[Nmodel].L-pos-1]=tmp_seq;
		}
		sub[Nmodel].sco = score2;
		Nmodel++;
		sub[Nmodel].L=0;
	}
	//END SUBOPT2----------------------------------------------
	
	//Mask aligned Sequence **Only** SUBOPT3
	//END SUBOPT3
	//Copy
	for(int i=0;i<t->Lpath*sfg->l;i++){
	 Smtx[i]=Smtx_ori[i];
	 SmtxRv[i]=SmtxRv_ori[i];
	}
	for(int i=0;i<t->Lpath;i++){
	for(int j=0;j<sfg->l;j++){
	 //Mask Node
	 if(MaskSeq[j]==true){
	  Smtx[GRID2D(n1,n2,i,j)]=-1.0;
	  if(j>0)
	   Smtx[GRID2D(n1,n2,i,j-1)]=-1.0;
	  if(j<sfg->l-1)
	   Smtx[GRID2D(n1,n2,i,j+1)]=-1.0;
	 }
	 if(MaskSeq[t->Lpath-1-j]==true){
	  SmtxRv[GRID2D(n1,n2,i,j)]=-1.0;
	  if(j>0)
	   SmtxRv[GRID2D(n1,n2,i,j-1)]=-1.0;
	  if(j<sfg->l-1)
	   Smtx[GRID2D(n1,n2,i,j+1)]=-1.0;
	 }
	}}
	score =  dp_local(Dmtx,Smtx,  dtbl,-100.0,-100.0,dpm->ali1,t->Lpath,ali2,  sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);
 	score2 = dp_local(Dmtx,SmtxRv,dtblR,-100.0,-100.0,dpm->ali1,t->Lpath,ali2rv,sfg->l,dpm->gali,&(dpm->Lgali),sfg->CID,true);
	
	if(score > ScoCut){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){//sequence
		 int j = ali2[pos];
			if(j!=-1){
			 //input to subopt
			 sub[Nmodel].model[sub[Nmodel].L] = n_path_po[j];
			 sub[Nmodel].seq[sub[Nmodel].L] = pos;
			 sub[Nmodel].L++;
			}
		}
		sub[Nmodel].sco = score;
		Nmodel++;
		sub[Nmodel].L=0;
 	}
	if(score2 > ScoCut){
		//input model and mask Smtx
		for(int pos=0;pos<sfg->l;pos++){
		 int j = ali2rv[pos];//Revese
			if(j!=-1){
			 //node and seq
			 sub[Nmodel].model[sub[Nmodel].L] = n_path_po[j];
			 sub[Nmodel].seq[sub[Nmodel].L] = sfg->l-1-pos;
			 sub[Nmodel].L++;
			}
		}
		//reverse
		for(int pos=0;pos<sub[Nmodel].L-pos-1;pos++){
		 NODE *node=sub[Nmodel].model[pos];
		 int tmp_seq=sub[Nmodel].seq[pos];
		 sub[Nmodel].model[pos] = sub[Nmodel].model[sub[Nmodel].L-pos-1];
		 sub[Nmodel].seq[pos] = sub[Nmodel].seq[sub[Nmodel].L-pos-1];
		 sub[Nmodel].model[sub[Nmodel].L-pos-1]=node;
		 sub[Nmodel].seq[sub[Nmodel].L-pos-1]=tmp_seq;
		}
		sub[Nmodel].sco = score2;
		Nmodel++;
		sub[Nmodel].L=0;
	}
	//EMD SUBOPT3------------------------------------
	//Update Find longest path and dtbl Smtx
	DFS_Tree(g,t);
	//printf("#Iter %d L= %d\n",iter,t->Lpath);
	for(int i=0;i<t->Lpath;i++)
  	 n_path_po[i]=&(g->node[t->Path[i]]);
/*//v02
	for(int i=0;i<t->Lpath-1;i++){
  	 dtbl[i+1]=dtbl[i]+sqrt(
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])*
	  (n_path_po[i]->real_cd[0]-n_path_po[i+1]->real_cd[0])+
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])*
	  (n_path_po[i]->real_cd[1]-n_path_po[i+1]->real_cd[1])+
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2])*
	  (n_path_po[i]->real_cd[2]-n_path_po[i+1]->real_cd[2]));
 	}
*/
	//Ver03 New dtbl position,position-j
	 for(int i=0;i<t->Lpath;i++){
		for(int j=0;j<21;j++)
		 dtbl[i*21+j]=10.00;//init
		for(int j=1;j<21 && i-j>-1;j++){
		 dtbl[i*21+j]=sqrt(
		  (n_path_po[i]->real_cd[0]-n_path_po[i-j]->real_cd[0])*
		  (n_path_po[i]->real_cd[0]-n_path_po[i-j]->real_cd[0])+
		  (n_path_po[i]->real_cd[1]-n_path_po[i-j]->real_cd[1])*
		  (n_path_po[i]->real_cd[1]-n_path_po[i-j]->real_cd[1])+
		  (n_path_po[i]->real_cd[2]-n_path_po[i-j]->real_cd[2])*
		  (n_path_po[i]->real_cd[2]-n_path_po[i-j]->real_cd[2]));
		}
	 }
	 for(int i=0;i<t->Lpath;i++){
	  int I = t->Lpath -1 - i;
		for(int j=0;j<21;j++)
		 dtblR[i*21+j]=10.00;//init
		for(int j=1;j<21 && i-j>-1;j++){
		 dtblR[i*21+j]=sqrt(
		  (n_path_po[I]->real_cd[0]-n_path_po[I+j]->real_cd[0])*
		  (n_path_po[I]->real_cd[0]-n_path_po[I+j]->real_cd[0])+
		  (n_path_po[I]->real_cd[1]-n_path_po[I+j]->real_cd[1])*
		  (n_path_po[I]->real_cd[1]-n_path_po[I+j]->real_cd[1])+
		  (n_path_po[I]->real_cd[2]-n_path_po[I+j]->real_cd[2])*
		  (n_path_po[I]->real_cd[2]-n_path_po[I+j]->real_cd[2]));
		}
	 }


	//Fill Smtx and SmtxRv
	for(int i=0;i<t->Lpath;i++){
        for(int j=0;j<sfg->l;j++){
         Saa =  n_path_po[i]->LogAA[sfg->seq[j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[j][2];
         Smtx[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	 //Mask
	 if(dpm->ali2[j]!=-1 || t->MaskN[t->Path[i]]==true)
	  Smtx[GRID2D(n1,n2,i,j)]=-1.0;
	 Smtx_ori[GRID2D(n1,n2,i,j)]=Smtx[GRID2D(n1,n2,i,j)];
        }}
 	for(int i=0;i<t->Lpath;i++){
        for(int j=0;j<sfg->l;j++){
         Saa =  n_path_po[i]->LogAA[sfg->seq[sfg->l-1-j]];
         Satm=  n_path_po[i]->LogATOM[2];//CA
         Sss=   n_path_po[i]->LogSS[0]*sfg->Pss[sfg->l-1-j][0]
                +n_path_po[i]->LogSS[1]*sfg->Pss[sfg->l-1-j][1]
                +n_path_po[i]->LogSS[2]*sfg->Pss[sfg->l-1-j][2];
         SmtxRv[GRID2D(n1,n2,i,j)]=Waa*Saa+Wss*Sss+Watm*Satm;
	 if(dpm->ali2[sfg->l-1-j]!=-1||t->MaskN[t->Path[i]]==true)
	  SmtxRv[GRID2D(n1,n2,i,j)]=-1.0;
	 SmtxRv_ori[GRID2D(n1,n2,i,j)]=SmtxRv[GRID2D(n1,n2,i,j)];
        }}
 }
 printf("#Nmodel= %d\n",Nmodel);
 t->score=final_score;
 return final_score;
}

bool ProbQA(PDB *p,GRAPH *g,int mode){
 float Qside[3],Qmain[3],Qca[3];
 float SeqSS[3];
 float *MapLogSS,*MapLogATOM,*MapLogAA;
 float Sss,Satm,Saa;
 int Nside,Nmain,Nca,Natm;
 float tmp[3];

 //All
 Qside[0]=Qmain[0]=Qca[0]=0.0;//ss
 Qside[1]=Qmain[1]=Qca[1]=0.0;//atm
 Qside[2]=Qmain[2]=Qca[2]=0.0;//aa

 Nside=Nmain=Nca=0;
 Natm=0;
 for(int i=0;i<p->NumOfAtom;i++){
  int atype=p->TypeAtomId[i];
  int j = p->AtomOnRes[i];
  int aa = p->TypeResId[j];

  MapLogSS=g->node[i].LogSS;
  MapLogATOM=g->node[i].LogATOM;
  MapLogAA=g->node[i].LogAA;
/*
  SeqSS[0]=p->SSP[3*i  ];
  SeqSS[1]=p->SSP[3*i+1];
  SeqSS[2]=p->SSP[3*i+2];
*/
//  Sss  = MapLogSS[0]*SeqSS[0] + MapLogSS[1]*SeqSS[1] + MapLogSS[2]*SeqSS[2];
  Sss=0.00;
  Satm = MapLogATOM[atype];
  Saa  = MapLogAA[aa];

  if(isinf(Sss)==1||isinf(Satm)==1||isinf(Saa)==1)
   printf("INF in %d atype=%d\n",p->ResNum[j],atype);
  if(isnan(Sss)==1||isnan(Satm)||isnan(Saa))
   printf("NAN in %d atype=%d\n",p->ResNum[j],atype);

  //Side-chains
  if(atype==0 ||atype==5 ){
   Qside[0]+=Sss; Qside[1]+=Satm; Qside[2]+=Saa;
   Nside++;
  }
  //main-chain
  if(atype>=1 && atype<=4){
   Qmain[0]+=Sss; Qmain[1]+=Satm; Qmain[2]+=Saa;
   Nmain++;
  }
  //CA positions
  if(atype==2){
   Qca[0]+=Sss; Qca[1]+=Satm; Qca[2]+=Saa;
   Nca++;
	 tmp[0]=p->xyz[i][0];
	 tmp[1]=p->xyz[i][1];
	 tmp[2]=p->xyz[i][2];
 	 printf("ATOM  %5d  CA  %3s %c%4d    ",Natm
	 ,RES_NAMES[aa]
	 ,p->Chain[i],p->ResNum[j]);
	 //if(mode==0)
	 // printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,Sss+Satm+Saa);
	 //if(mode==1)
	 //printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,Sss);
	 //if(mode==2)
	 //printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,Satm);
	 //if(mode==3) //Show Only AAscore
	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,Saa);
	 //printf("#SS= %.3f ATOM= %.3f AA= %.3f TOTAL= %.3f\n",Sss,Satm,Saa,Sss+Satm+Saa);
	 //printf("#H_SS %.2f * %.2f\n",MapLogSS[0],SeqSS[0]);
	 //printf("#E_SS %.2f * %.2f\n",MapLogSS[1],SeqSS[1]);
	 //printf("#C_SS %.2f * %.2f\n",MapLogSS[2],SeqSS[2]);
	 Natm++;
  }
 }

 printf("##Side: SS: %.1f ATOM: %.1f AA: %.1f N: %d\n",Qside[0],Qside[1],Qside[2],Nside);
 printf("#ASide: SS: %.3f ATOM: %.3f AA: %.3f\n",Qside[0]/(float)(Nside),Qside[1]/(float)(Nside),Qside[2]/(float)(Nside));
 printf("##Main: SS: %.1f ATOM: %.1f AA: %.1f N: %d\n",Qmain[0],Qmain[1],Qmain[2],Nmain);
 printf("#AMain: SS: %.3f ATOM: %.3f AA: %.3f\n",Qmain[0]/(float)(Nmain),Qmain[1]/(float)(Nmain),Qmain[2]/(float)(Nmain));
 printf("##CA  : SS: %.1f ATOM: %.1f AA: %.1f N: %d\n",Qca[0],Qca[1],Qca[2],Nca);
 printf("#ACA  : SS: %.3f ATOM: %.3f AA: %.3f\n",Qca[0]/(float)(Nca),Qca[1]/(float)(Nca),Qca[2]/(float)(Nca));
 return false;
}

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "struct.h"
#include "mrc.h"

extern CMD cmd;



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

 //Memory
 if((mrc->dens=(float*)malloc(sizeof(float)*mrc->NumVoxels))==NULL)
  return true;
 if((mrc->inside=(bool*)malloc(sizeof(bool)*mrc->NumVoxels))==NULL)
  return true;
 

 //fin.ignore(4*(256-54)+nsymbt);
 fseek(fpin,4*(256-54)+nsymbt,SEEK_CUR);
 
 	switch(mode) {
 	 case 0: // char - converted to float, testing for signed-ness
	   printf("Cannot read mode 0 mrc file\n");
           return true;
	   //break;
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
	float ncstart,nrstart,nsstart;
	switch(ordermode){
		case 1:
			nx = xdim; ny = ydim; nz = zdim;
			ncstart=mrc->ncstart,nrstart=mrc->nrstart;nsstart=mrc->nsstart;
			break;
		case 2:
			nx = xdim; ny = zdim; nz = ydim;
			ncstart=mrc->ncstart,nrstart=mrc->nsstart;nsstart=mrc->nrstart;
			break;
		case 3:
			nx = ydim; ny = xdim; nz = zdim;
			ncstart=mrc->nrstart,nrstart=mrc->ncstart;nsstart=mrc->nsstart;
			break;
		case 4:
			nx = zdim; ny = xdim; nz = ydim;
			ncstart=mrc->nsstart,nrstart=mrc->ncstart;nsstart=mrc->nrstart;
			break;
		case 5:
			nx = ydim; ny = zdim; nz = xdim;
			ncstart=mrc->nrstart,nrstart=mrc->nsstart;nsstart=mrc->ncstart;
			break;
		case 6:
			nx = zdim; ny = ydim; nz = xdim;
			ncstart=mrc->nsstart,nrstart=mrc->nrstart;nsstart=mrc->ncstart;
			break;
		default:
			printf("Input file gives malformed dimension ordering.");
			return true;
	}
 mrc->ncstart=ncstart;
 mrc->nrstart=nrstart;
 mrc->nsstart=nsstart;
 mrc->xdim = nx;
 mrc->ydim = ny;
 mrc->zdim = nz;

 printf("#XYZ dim: %d %d %d\n",mrc->xdim,mrc->ydim,mrc->zdim);
 printf("#Grid size: X %f  Y %f Z %f\n",
	  mrc->widthx,mrc->widthy,mrc->widthz);
 printf("#updated new crs: %.5f, %.5f, %.5f\n",ncstart,nrstart,nsstart);
 fclose(fpin);


/*
 //NCSTART is wrong?
 if(mrc->ncstart!=0||mrc->nrstart!=0||mrc->nsstart!=0){
  mrc->orgxyz[0]=mrc->orgxyz[0]+mrc->ncstart*mrc->widthx;
  mrc->orgxyz[1]=mrc->orgxyz[1]+mrc->nrstart*mrc->widthy;
  mrc->orgxyz[2]=mrc->orgxyz[2]+mrc->nsstart*mrc->widthz;
  printf("#orgXYZ: %f %f %f\n",mrc->orgxyz[0],mrc->orgxyz[1],mrc->orgxyz[2]);
 }
 */


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
 if(mrc->widthx < 1.0){
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



bool fastVEC(MRC *m,MRC *M){
 int i,j,k,ind;
 int cnt=0;
 int xydim=m->xdim*m->ydim;
 int Ndata=M->xdim*M->ydim*M->zdim;
 //malloc

 if((M->vec=(double **)malloc(sizeof(double *)*Ndata))==NULL)
  return true;
 if((M->dens=(float *)malloc(sizeof(double)*Ndata))==NULL)
  return true;
 for(i=0;i<Ndata;i++){
  if((M->vec[i]=(double *)malloc(sizeof(double)*3))==NULL)
  return true;
 }


 puts("#Start VEC");
 //Setup Filter
 //Gaussian kernel dreso=window size
 double dreso=cmd.dreso;
 double gstep=m->widthx;
 double fs=(dreso/gstep)*0.5;
 fs=fs*fs;
 double fsiv=1.000/fs;
 double fmaxd=(dreso/gstep)*2.0;
 printf("#maxd= %f\n",fmaxd);

 double dsum=0;
 int Nact=0;

 #pragma omp parallel for reduction(+:Nact) schedule(dynamic,5)
 for(int x=0;x<M->xdim;x++){
  double rx,ry,rz,d2;
 for(int y=0;y<M->ydim;y++){
 for(int z=0;z<M->zdim;z++){

  int stp[3],endp[3],ind2,ind;
  double pos[3],pos2[3],ori[3];
  double tmpcd[3];
  double v,dtotal,rd;


  //real cd -> m-cd
  pos[0]=(x*M->widthx+M->orgxyz[0]-m->orgxyz[0])/m->widthx;
  pos[1]=(y*M->widthx+M->orgxyz[1]-m->orgxyz[1])/m->widthx;
  pos[2]=(z*M->widthx+M->orgxyz[2]-m->orgxyz[2])/m->widthx;

  ind=M->xdim*M->ydim*z+M->xdim*y+x;
  //printf("%d %d %d %d\n",x,y,z,ind);
  //printf("%f %f %f %d\n",pos[0],pos[1],pos[2],ind);

  //check density
  if(pos[0]<0||pos[1]<0||pos[2]<0||
     pos[0]>=m->xdim||pos[1]>=m->ydim||pos[2]>=m->zdim){
   M->dens[ind]=0;
   M->vec[ind][0]=M->vec[ind][1]=M->vec[ind][2]=0.00;
   continue;
  }





  int ind0=m->xdim*m->ydim*(int)pos[2]+m->xdim*(int)pos[1]+(int)pos[0];
  if(m->dens[ind0]==0){
   M->dens[ind]=0;
   M->vec[ind][0]=M->vec[ind][1]=M->vec[ind][2]=0.00;
   continue;
  }

  ori[0]=pos[0];
  ori[1]=pos[1];
  ori[2]=pos[2];
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
    //if(v>0)
    //printf("d %f %d %d %d\n",v,xp,yp,zp);
    pos2[0]+=v*(double)xp;
    pos2[1]+=v*(double)yp;
    pos2[2]+=v*(double)zp;
   }}}
   M->dens[ind]=dtotal;
   if(dtotal==0.00){
    M->vec[ind][0]=M->vec[ind][1]=M->vec[ind][2]=0.00;
    continue;
   }
   rd=1.00/dtotal;
   pos2[0]*=rd;
   pos2[1]*=rd;
   pos2[2]*=rd;
   tmpcd[0]=pos2[0]-pos[0];
   tmpcd[1]=pos2[1]-pos[1];
   tmpcd[2]=pos2[2]-pos[2];

   double dvec=sqrt(tmpcd[0]*tmpcd[0]+tmpcd[1]*tmpcd[1]+tmpcd[2]*tmpcd[2]);
   if(dvec==0.000)
	   dvec=1.000;
   double rdvec=1.000/dvec;
   M->vec[ind][0]=tmpcd[0]*rdvec;
   M->vec[ind][1]=tmpcd[1]*rdvec;
   M->vec[ind][2]=tmpcd[2]*rdvec;
   //printf("dto= %f [%f %f %f]\n",dtotal,M->vec[ind][0],M->vec[ind][1],M->vec[ind][2]);

   dsum+=dtotal;
   Nact++;

 }}}
 puts("#End LDP");

 //Add Average & STD
 M->dsum=dsum;
 M->Nact=Nact;
 M->ave=dsum/(double)Nact;
 dsum=0;
 double dsum2=0;
 for(int i=0;i<M->xdim*M->ydim*M->zdim;i++)
  if(M->dens[i]>0){//Cross correlation
   dsum+=(M->dens[i])*(M->dens[i]);
   dsum2+=(M->dens[i]-M->ave)*(M->dens[i]-M->ave);
  }
 M->std_norm_ave=sqrt(dsum2);
 M->std=sqrt(dsum);
 printf("#MAP AVE= %f STD= %f STD_norm= %f\n",M->ave,M->std,M->std_norm_ave);
 return false;
}

void ShowVec(MRC *M){

 int i,ind;
 int Natm=1;
 int Nres=1;
 double tmp[3];
 puts("MODEL");
 for(int x=0;x<M->xdim;x++){
 for(int y=0;y<M->ydim;y++){
 for(int z=0;z<M->zdim;z++){
  ind=M->xdim*M->ydim*z+M->xdim*y+x;

  i=ind;
  if(M->dens[i]==0.00&&ind!=0) continue;
  tmp[0]=x*M->widthx+M->orgxyz[0];
  tmp[1]=y*M->widthx+M->orgxyz[1];
  tmp[2]=z*M->widthx+M->orgxyz[2];
  printf("ATOM  %5d  CA  ALA%6d    ",Natm,Nres);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,M->dens[i]);
  Natm++;


  tmp[0]=(x+M->vec[i][0])*M->widthx+M->orgxyz[0];
  tmp[1]=(y+M->vec[i][1])*M->widthx+M->orgxyz[1];
  tmp[2]=(z+M->vec[i][2])*M->widthx+M->orgxyz[2];
  printf("ATOM  %5d  CB  ALA%6d    ",Natm,Nres);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,M->dens[i]);
  Natm++;

  Nres++;
  
 }}}
 puts("ENDMDL");
}

void ShowVec2(MRC *m,MRC *M,int t[3]){

 int i,ind;
 int Natm=1;
 int Nres=1;
 double tmp[3];


 double add[3];

 if(t[0]>0.5*M->zdim) t[0]=t[0]-M->xdim;
 if(t[1]>0.5*M->zdim) t[1]=t[1]-M->xdim;
 if(t[2]>0.5*M->zdim) t[2]=t[2]-M->xdim;

 add[0]=m->orgxyz[0]-t[0]*M->widthx;
 add[1]=m->orgxyz[1]-t[1]*M->widthx;
 add[2]=m->orgxyz[2]-t[2]*M->widthx;

 puts("MODEL");
 for(int x=0;x<M->xdim;x++){
 for(int y=0;y<M->ydim;y++){
 for(int z=0;z<M->zdim;z++){
  ind=M->xdim*M->ydim*z+M->xdim*y+x;

  i=ind;
  if(M->dens[i]==0.00) continue;
  tmp[0]=x*M->widthx+add[0];
  tmp[1]=y*M->widthx+add[1];
  tmp[2]=z*M->widthx+add[2];
  printf("ATOM  %5d  CA  ALA%6d    ",Natm,Nres);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,M->dens[i]);
  Natm++;


  tmp[0]=(x+M->vec[i][0])*M->widthx+add[0];
  tmp[1]=(y+M->vec[i][1])*M->widthx+add[1];
  tmp[2]=(z+M->vec[i][2])*M->widthx+add[2];
  printf("ATOM  %5d  CB  ALA%6d    ",Natm,Nres);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,M->dens[i]);
  Natm++;

  Nres++;
  
 }}}
 puts("ENDMDL");
}

//Score
void ShowVec3(MRC *m,MRC *M,int t[3]){

 int i,ind;
 int Natm=1;
 int Nres=1;
 double tmp[3];


 double add[3];

 if(t[0]>0.5*M->zdim) t[0]=t[0]-M->xdim;
 if(t[1]>0.5*M->zdim) t[1]=t[1]-M->xdim;
 if(t[2]>0.5*M->zdim) t[2]=t[2]-M->xdim;

 add[0]=m->orgxyz[0]-t[0]*M->widthx;
 add[1]=m->orgxyz[1]-t[1]*M->widthx;
 add[2]=m->orgxyz[2]-t[2]*M->widthx;

 puts("MODEL");
 for(int x=0;x<M->xdim;x++){
 for(int y=0;y<M->ydim;y++){
 for(int z=0;z<M->zdim;z++){
  ind=M->xdim*M->ydim*z+M->xdim*y+x;

  i=ind;
  if(M->dens[i]==0.00) continue;
  tmp[0]=x*M->widthx+add[0];
  tmp[1]=y*M->widthx+add[1];
  tmp[2]=z*M->widthx+add[2];
  printf("ATOM%7d  CA  ALA%6d    ",Natm,Nres);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,M->sco[i]);
  Natm++;


  tmp[0]=(x+M->vec[i][0])*M->widthx+add[0];
  tmp[1]=(y+M->vec[i][1])*M->widthx+add[1];
  tmp[2]=(z+M->vec[i][2])*M->widthx+add[2];
  printf("ATOM%7d  CB  ALA%6d    ",Natm,Nres);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,M->sco[i]);
  Natm++;

  Nres++;
  
 }}}
 puts("ENDMDL");
}

void ShowMRC(MRC *M){

 int i,ind;
 int Natm=1;
 int Nres=1;
 double tmp[3];
 puts("MODEL");
 for(int x=0;x<M->xdim;x+=4){
 for(int y=0;y<M->ydim;y+=4){
 for(int z=0;z<M->zdim;z+=4){
  ind=M->xdim*M->ydim*z+M->xdim*y+x;

  i=ind;
  if(M->dens[i]==0.00&&ind!=0) continue;
  tmp[0]=x*M->widthx+M->orgxyz[0];
  tmp[1]=y*M->widthx+M->orgxyz[1];
  tmp[2]=z*M->widthx+M->orgxyz[2];
  printf("ATOM  %5d  CA  ALA%6d    ",Natm,Nres);
  printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,M->dens[i]);
  Natm++;
  Nres++;
  
 }}}
 puts("ENDMDL");
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
  //printf("%f %f %f\n",p->cd[i][0],m->widthx,m->orgxyz[0]);
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

void ShowPath2(MRC *m,POINTS *p,GRAPH *g, TREE *t, int n){
 //Show path
 //For density*volume, using K-clustring data
 int i,j,k,now,ind;
 int Natm=1;
 double tmp[3];
 int st,ed;
 double fmax=0;
 int xydim=m->xdim*m->ydim;
 //path data
 
 


 for(int mid=0;mid<n;mid++){
  st=t[mid].St;
  ed=t[mid].Ed;
  fmax=0;
  //get_path(g, &t[mid],st,ed);

	//Assign Volume and Density Data
	for(i=0;i<t[mid].Lpath;i++)
	 t[mid].cost[i]=0;
	
	//for(int x=0;x<m->xdim;x++){
	//for(int y=0;y<m->ydim;y++){
	//for(int z=0;z<m->zdim;z++){

	for(int j=0;j<p->Nori;j++){

	 int x=p->origrid[j][0];
	 int y=p->origrid[j][1];
	 int z=p->origrid[j][2];


	 ind=xydim*z+m->xdim*y+x;
	 if(p->mask[j]==0.00||m->dens[ind]==0.00)
	  continue;
  	 //map origin xyz and xwidth based coordinates
  	 tmp[0]=(double)x; tmp[1]=(double)y; tmp[2]=(double)z;

	 int minid=-1;
	 double mindis=100000.0;
	 double d2;
		//Search closest path point
		for(int rnum=0;rnum<t[mid].Lpath;rnum++){
		 i=t[mid].Path[rnum];
		 d2=(p->cd[i][0]-tmp[0])*(p->cd[i][0]-tmp[0])
		   +(p->cd[i][1]-tmp[1])*(p->cd[i][1]-tmp[1])
		   +(p->cd[i][2]-tmp[2])*(p->cd[i][2]-tmp[2]);
  		 if(d2<mindis){
		  minid=rnum;
		  mindis=d2;
		 }
		}

	 if(minid!=-1)
	  t[mid].cost[minid]+=m->dens[ind];
	 //printf("minid= %d d= %f\n",minid,sqrt(d2));
	}
	//}}}

	for(i=0;i<t[mid].Lpath;i++)
	 if(fmax<t[mid].cost[i])
	  fmax=t[mid].cost[i];



  printf("#SCORE: %f LEN=%d fmax=%f\n",t[mid].score,t[mid].Lpath,fmax);
  printf("MODEL %d\n",mid+1);



	Natm=1;
	for(int rnum=0;rnum<t[mid].Lpath;rnum++){
	 i=t[mid].Path[rnum];
	 tmp[0]=p->cd[i][0]*m->widthx+m->orgxyz[0];
 	 tmp[1]=p->cd[i][1]*m->widthx+m->orgxyz[1];
 	 tmp[2]=p->cd[i][2]*m->widthx+m->orgxyz[2];
 	 printf("ATOM  %5d  CA  ALA%6d    ",Natm,Natm);
 	 //printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,p->dens[i]);
 	 printf("%8.3f%8.3f%8.3f%6.2f%6.2f\n",tmp[0],tmp[1],tmp[2],1.0,t[mid].cost[rnum]/fmax);
 	 Natm++;
	}
  printf("TER\nENDMDL\n");
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
 //Remove isolated tree
 printf("#MaxCid= %d\n",MaxCid);
 int *Ncid;
 if((Ncid=(int *)calloc(sizeof(int),(MaxCid+1)))==NULL)
  return true;
 for(i=0;i<g->Nnode;i++)
  if(g->cid[i]<=MaxCid)
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
 
 //clean edge shift
 int Ntmp=0;
 for(i=0;i<Ne;i++){
  if(UseCid==g->cid[g->edge[i].id1]){
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
 
 printf("#After cleaning.. Nt= %d Ne= %d\n",Nt,Ntmp);
 
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

 //Local MST
 for(int ii=0;ii<g->Nnode;ii++){
  double vec[3],d2;
  int tmpid;
  //init cid
  for(int jj=0;jj<g->Nnode;jj++)
   cid[jj]=jj;

  	for(int jj=0;jj<g->Ne;jj++){
  	 int v1=g->edge[jj].id1;
  	 int v2=g->edge[jj].id2;
  	 if(cid[v1]==cid[v2])
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

	 tmpid=cid[v2];
	 for(int kk=0;kk<g->Nnode;kk++){
	  if(cid[kk]==tmpid)//update cid
	   cid[kk]=cid[v1];
	 }
  	}
 }
 //End Local MST
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
bool Tabu(GRAPH *g,TREE *init_t,TREE *Tr){
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
 //CopyTree(init_t,&Tr[0],true);
 puts("#Fin Malloc for Tr");

 CopyTree(init_t,&Tnow,true);
 puts("#Fin Malloc for Tnow");
 CopyTree(init_t,&Tbest,true);
 puts("#Fin Malloc for Tbest");

 //return true;

 //double score;
 double score=QualityTree(g,&Tnow);

 //printf("score= %f\n",score);
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

	 #pragma omp parallel for schedule(dynamic,5)
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

	  //printf("NB: %4d MV:cut:%d add:%d ",nb,move[nb].cut_id,move[nb].add_id);
	  move[nb].score=QualityTree(g,&Tnb[nb]);

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

/*
	 printf("Stock:");
	 for(int k=0;k<Nst;k++)
	  printf(" %d",t->stock[k]);
	 printf("\n");
*/
	}
  //printf("#Maxi=  %d -> %d pre= %d dist=%f\n",st,maxi,pre_st,maxd);
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

void SetUpVoxSize(MRC *m,MRC *M,double t,double ssize){
 int Ndata=m->xdim*m->ydim*m->zdim;
 int xydim=m->xdim*m->ydim;
 int xdim=m->xdim;

 double cent[3];
 double d2,dmax=0;

 cent[0]=m->xdim*0.5;
 cent[1]=m->ydim*0.5;
 cent[2]=m->zdim*0.5;

 //If t<0
 if(t<0){
        for(int x=0;x<m->xdim;x++){
        for(int y=0;y<m->ydim;y++){
        for(int z=0;z<m->zdim;z++){
         int ind=xydim*z+xdim*y+x;
         m->dens[ind]-=t;
        }}}
        t=0.00;
 }


 for(int x=0;x<m->xdim;x++){
 for(int y=0;y<m->ydim;y++){
 for(int z=0;z<m->zdim;z++){
  int ind=xydim*z+xdim*y+x;
  if(m->dens[ind]<t){
   m->dens[ind]=0.00;//filter
   continue;
  }

  d2=((double)x-cent[0])*((double)x-cent[0])
    +((double)y-cent[1])*((double)y-cent[1])
    +((double)z-cent[2])*((double)z-cent[2]);
  if(d2>dmax)
   dmax=d2;
 }}}
 printf("#dmax= %f size=%d\n",sqrt(dmax)/m->widthx,(int)(2*sqrt(dmax)/m->widthx));

 //real coordinates
 M->cent[0]=cent[0]*m->widthx+m->orgxyz[0];
 M->cent[1]=cent[1]*m->widthx+m->orgxyz[1];
 M->cent[2]=cent[2]*m->widthx+m->orgxyz[2];

 M->widthx=ssize;

 m->dmax=sqrt(dmax)*m->widthx;
 int tmp_size=2*sqrt(dmax)*m->widthx/M->widthx;

 int a=2;
 while(1){
  if(a>tmp_size)
   break;
  a*=2;
 }
 int b=3;
 while(1){
  if(b>tmp_size)
   break;
  b*=2;
 }
 if(b<a)
  a=b;
 b=9;
 while(1){
  if(b>tmp_size)
   break;
  b*=2;
 }
 if(b<a)
  a=b;

 M->xdim=M->ydim=M->zdim=a;
 M->orgxyz[0]=M->cent[0]-0.5*a*M->widthx;
 M->orgxyz[1]=M->cent[1]-0.5*a*M->widthx;
 M->orgxyz[2]=M->cent[2]-0.5*a*M->widthx;
 printf("Nvox= %d*%d*%d\n",M->xdim,M->ydim,M->zdim);
 printf("cent= %f %f %f\n",M->cent[0],M->cent[1],M->cent[2]);
 printf("ori= %f %f %f\n",M->orgxyz[0],M->orgxyz[1],M->orgxyz[2]);
}

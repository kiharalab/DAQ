//Dynamic Programming
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include <math.h>
#include "dp.h"
/*
#define GRID(N,a,b,c) c+N*(b+N*a)
#define GRID3D(N1,N2,a,b,c) c+N2*(b+N1*a)
#define GRID2D(N1,N2,a,b) (b)+(N2)*(a)
*/

//Smith & Waterman Global Alignment
//pointer
#define NON -1
#define UP 0
#define LEFT 1
#define DIA 2
#define GAP -1

/*
typedef struct{ 
	double sco; 
	int poi;
} DPMTX;
*/
//Semi-global alignment
//Gaps of N and C-terminal are not penalized
/*
!!No gap penalty
AAAAAAAAAAAAAA
---BBBBBBB----

*/


float dp(DPMTX *dmtx,float *Smtx,float GapOpen,float GapExt,int *al1,int n1,int *al2,int n2,int *gal,int *glen){
 //puts("#start DP");
 //DPMTX dmtx[RES*RES];
 //DPMTX dmtx[n1*n2];
 int i1,i2;
 int N1,N2;
 float dia_sco,up_sco,left_sco;
 int gridid;
 //Warning!!
 /*
 0:0 is not corresponding with Smtx[GRID2D(n1,n2,0,0)]
 0:x and x:1 is gap data of N-terminal
 use GRID2D(N1,N2,a+1,b+1)
 */


 N1=n1+1;
 N2=n2+1;

 dmtx[0].sco=0;
 dmtx[0].poi=NON;
 //init N-terminal DPMTX data
 for(i1=1;i1<=n1;i1++){
  dmtx[GRID2D(N1,N2,i1,0)].sco=0;
  dmtx[GRID2D(N1,N2,i1,0)].poi=NON;
 }
 for(i2=1;i2<=n2;i2++){
  dmtx[GRID2D(N1,N2,0,i2)].sco=0;
  dmtx[GRID2D(N1,N2,0,i2)].poi=NON;
 }
 //puts("init");
 //fill
 for(i1=1;i1<=n1;i1++){
 	for(i2=1;i2<=n2;i2++){
/*
	 if(i1==i2){
	 printf("%d-%d SCO=%.2f\n",i1,i2,Smtx[GRID2D(n1,n2,i1-1,i2-1)]);
	 printf("%d*%d %d %d =%d\n",n1,n2,i1-1,i2-1,GRID2D(n1,n2,i1-1,i2-1));
	 }
*/
	 //diagonal score
	 dia_sco=dmtx[GRID2D(N1,N2,i1-1,i2-1)].sco+Smtx[GRID2D(n1,n2,i1-1,i2-1)];
	 //gap score
	 //up score
	 up_sco=dmtx[GRID2D(N1,N2,i1-1,i2)].sco;
	 //gap open? or extention?
	 //open skip pdb
	 if(dmtx[GRID2D(N1,N2,i1-1,i2)].poi==DIA||dmtx[GRID2D(N1,N2,i1-1,i2)].poi==LEFT||dmtx[GRID2D(N1,N2,i1-1,i2)].poi==NON){
	  up_sco+=GapOpen;
	 }else{
	  up_sco+=GapExt;
	 }
	 //left score
	 left_sco=dmtx[GRID2D(N1,N2,i1,i2-1)].sco;
	 //gap open? or extention?
	 //open
	 if(dmtx[GRID2D(N1,N2,i1,i2-1)].poi==DIA){
	  //left_sco+=GapOpen;
	  left_sco+=0;
	 }else if(dmtx[GRID2D(N1,N2,i1,i2-1)].poi==UP){
	  left_sco+=GapOpen;
	 }else{
	  //skip sequence
	  //left_sco+=GapExt;
	  left_sco+=0;
	 }
	 //choose max score
	 if(dia_sco>=up_sco){
		if(dia_sco>=left_sco){
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=DIA;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=dia_sco;
		}else{
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=LEFT;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=left_sco;
		}
	 }else{
		if(up_sco>=left_sco){
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=UP;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=up_sco;
		}else{
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=LEFT;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=left_sco;
		}
	 }
	 //if(i1>889 && i1<896 && i2>887 && i2<912)
	 //printf("%d:%d dia:%.3f up%.3f left:%.3f sco:%.3f poi %d\n",i1,i2,dia_sco,up_sco,left_sco,dmtx[GRID2D(N1,N2,i1,i2)].sco,dmtx[GRID2D(N1,N2,i1,i2)].poi);

	}
 }
 //trace back
 //start from highest score bottom or end of columms
 //end zero gaps
 float highest=0;
 int i=0,j=0;
 i=n1;j=n2;
/*
 for(i1=0;i1<=n1;i1++){
  if(highest < dmtx[GRID2D(N1,N2,i1,n2)].sco){
   highest=dmtx[GRID2D(N1,N2,i1,n2)].sco;
   i=i1;
   j=n2;
  }
 }
*/
 for(i2=0;i2<=n2;i2++){
  if(highest < dmtx[GRID2D(N1,N2,n1,i2)].sco){
   highest=dmtx[GRID2D(N1,N2,n1,i2)].sco;
   i=n1;
   j=i2;
  }
 }

 printf("#START FROM %d %d sco=%.3f\n",i,j,highest);
 //printf("START FROM %d %d sco=%.3f\n",n1,n2,dmtx[GRID2D(N1,N2,n1-1,n2-1)].sco);
 //puts("#TRACE BACK");
 int gpos=0;//reverse position for global alignment
 //Becareful!! i&j is start from 1
/*
 while(1){
 	gridid=GRID2D(N1,N2,i,j);
	printf("%d-%d %.3f\n",i,j,dmtx[GRID2D(N1,N2,i,j)].sco);
	if(dmtx[gridid].poi==NON)//end when 0:x or x:0
	 break;
	if(dmtx[gridid].poi==DIA){
	 al1[i]=j; al2[j]=i;
	 gal[gpos*2]=i; gal[gpos*2+1]=j; gpos++;
	 i--;j--;
	}else if(dmtx[gridid].poi==LEFT){
	 al1[i]=j; al2[j]=GAP;
	 gal[gpos*2]=GAP; gal[gpos*2+1]=j; gpos++;
	 j--;
	}else if(dmtx[gridid].poi==UP){
	 al1[i]=-1; al2[j]=i;
	 gal[gpos*2]=i; gal[gpos*2+1]=GAP; gpos++;
	 i--;
	}
 }*/
 int tmp_gal[20000];
 while(1){
 	gridid=GRID2D(N1,N2,i,j);
	//printf("%d-%d %.3f poi= %d\n",i,j,dmtx[GRID2D(N1,N2,i,j)].sco,dmtx[GRID2D(N1,N2,i,j)].poi);
	if(dmtx[gridid].poi==NON)//end when 0:x or x:0
	 break;
	if(dmtx[gridid].poi==DIA){
	 //al1[i]=j; al2[j]=i;
	 tmp_gal[gpos*2]=i-1; tmp_gal[gpos*2+1]=j-1; gpos++;
	 i--;j--;
	}else if(dmtx[gridid].poi==LEFT){
	 tmp_gal[gpos*2]=GAP; tmp_gal[gpos*2+1]=j-1; gpos++;
	 j--;
	}else if(dmtx[gridid].poi==UP){
	 tmp_gal[gpos*2]=i-1; tmp_gal[gpos*2+1]=GAP; gpos++;
	 i--;
	}
 }

 *glen=gpos;
 //put to gal
 for(int i=0;i<gpos;i++){
  gal[i*2  ] = tmp_gal[2*(gpos-i-1)];
  gal[i*2+1] = tmp_gal[2*(gpos-i-1)+1];
 }
 //init al1
 for(int i=0;i<n1;i++)
  al1[i]=-1;
 //put al1
 for(int i=0;i<gpos;i++){
  if(gal[2*i]!=-1)
   al1[gal[2*i]]=gal[2*i+1];
 }
 //init al2
 for(int i=0;i<n2;i++)
  al2[i]=-1;
 //put al1
 for(int i=0;i<gpos;i++){
  if(gal[2*i+1]!=-1)
   al2[gal[2*i+1]]=gal[2*i];
 }

 return highest;
}



//Global alignment
float dp_fast(DPMTX *dmtx,float *Smtx,float *dtbl,float GapOpen,float GapExt,int *al1,int n1,int *al2,int n2,int *gal,int *glen,bool mode){
 int i1,i2,p1;
 int N1,N2;
 float dia_sco,up_sco,left_sco;
 int gridid;
 //Warning!!
 /*
 0:0 is not corresponding with Smtx[GRID2D(n1,n2,0,0)]
 0:x and x:1 is gap data of N-terminal
 use GRID2D(N1,N2,a+1,b+1)
 */


 N1=n1+1;
 N2=n2+1;

 dmtx[0].sco=0;
 dmtx[0].poi=NON;
 //init N-terminal DPMTX data
 for(i1=1;i1<=n1;i1++){
  dmtx[GRID2D(N1,N2,i1,0)].sco=0;
  dmtx[GRID2D(N1,N2,i1,0)].poi=NON;
  dmtx[GRID2D(N1,N2,i1,0)].pre_align1=-1;
  dmtx[GRID2D(N1,N2,i1,0)].pre_align2=-1;
 }
 for(i2=1;i2<=n2;i2++){
  dmtx[GRID2D(N1,N2,0,i2)].sco=GapOpen;
  dmtx[GRID2D(N1,N2,0,i2)].poi=NON;
  dmtx[GRID2D(N1,N2,0,i2)].pre_align1=-1;
  dmtx[GRID2D(N1,N2,0,i2)].pre_align2=-1;
 }
 //puts("init");
 //fill
 for(i1=1;i1<=n1;i1++){
 	for(i2=1;i2<=n2;i2++){
/*
	 if(i1==i2){
	 printf("%d-%d SCO=%.2f\n",i1,i2,Smtx[GRID2D(n1,n2,i1-1,i2-1)]);
	 printf("%d*%d %d %d =%d\n",n1,n2,i1-1,i2-1,GRID2D(n1,n2,i1-1,i2-1));
	 }
*/
	 //diagonal score
	 dia_sco=dmtx[GRID2D(N1,N2,i1-1,i2-1)].sco+Smtx[GRID2D(n1,n2,i1-1,i2-1)];
	 //Node distance loss
	 p1 = dmtx[GRID2D(N1,N2,i1-1,i2-1)].pre_align1;
	 if(p1 != -1){
	  if((dtbl[i1-1]-dtbl[p1])<3.2)
	   dia_sco+=-1.0*((dtbl[i1-1]-dtbl[p1])-3.00)*((dtbl[i1-1]-dtbl[p1])-3.20);
	  if((dtbl[i1-1]-dtbl[p1])>3.8)
	   dia_sco+=-1.0*((dtbl[i1-1]-dtbl[p1])-4.00)*((dtbl[i1-1]-dtbl[p1])-3.80);
	 }
	 //gap score
	 //up score
	 up_sco=dmtx[GRID2D(N1,N2,i1-1,i2)].sco;
/*
	 //gap open? or extention?
	 //open
	 if(dmtx[GRID2D(N1,N2,i1-1,i2)].poi==DIA){
	  up_sco-=GapOpen;
	 }else{
	  up_sco-=GapExt;
	 }
*/
	 //left score
	 left_sco=dmtx[GRID2D(N1,N2,i1,i2-1)].sco;
	 //gap open? or extention?

	 //open
	 if(dmtx[GRID2D(N1,N2,i1,i2-1)].poi==DIA){
	  left_sco+=GapOpen;
	 }else{
	  left_sco+=GapExt;
	 }
	 //choose max score
	 if(dia_sco>=up_sco){
		if(dia_sco>=left_sco){
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=DIA;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=dia_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=i1-1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=i2-1;
		}else{
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=LEFT;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=left_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align2;
		}
	 }else{
		if(up_sco>=left_sco){
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=UP;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=up_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=dmtx[GRID2D(N1,N2,i1-1,i2)].pre_align1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=dmtx[GRID2D(N1,N2,i1-1,i2)].pre_align2;
		}else{
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=LEFT;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=left_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align2;
		}
	 }

	 //printf("%d:%d dia:%.3f up%.3f left:%.3f sco:%.3f poi %d\n",i1,i2,dia_sco,up_sco,left_sco,dmtx[GRID2D(N1,N2,i1,i2)].sco,dmtx[GRID2D(N1,N2,i1,i2)].poi);

	}
 }
 //trace back

 //start from highest score bottom or end of columms
 //end zero gaps
 float highest=0;
 int i=0,j=0;

 for(i1=1;i1<=n1;i1++){
  if(highest < dmtx[GRID2D(N1,N2,i1,n2)].sco){
   highest=dmtx[GRID2D(N1,N2,i1,n2)].sco;
   i=i1;
   j=n2;
  }
 }
/*
 for(i2=1;i2<=n2;i2++){
  if(highest < dmtx[GRID2D(N1,N2,n1,i2)].sco){
   highest=dmtx[GRID2D(N1,N2,n1,i2)].sco;
   i=n1;
   j=i2;
  }
 }
*/

 //printf("#START FROM %d %d sco=%.3f\n",i,j,highest);
 //printf("START FROM %d %d sco=%.3f\n",n1,n2,dmtx[GRID2D(N1,N2,n1-1,n2-1)].sco);
 //puts("#TRACE BACK");
 if(mode==false)//no alignment
  return highest;
 int gpos=0;//reverse position for global alignment
 short int tmp_gal[200];
 //Becareful!! i&j is start from 1
 while(1){
 	gridid=GRID2D(N1,N2,i,j);
	//printf("%d-%d %.3f\n",i,j,dmtx[GRID2D(N1,N2,i,j)].sco);
	if(dmtx[gridid].poi==NON)//end when 0:x or x:0
	 break;
	if(dmtx[gridid].poi==DIA){
	 //al1[i]=j; al2[j]=i;
	 tmp_gal[gpos*2]=i-1; tmp_gal[gpos*2+1]=j-1; gpos++;
	 i--;j--;
	}else if(dmtx[gridid].poi==LEFT){
	 //al1[i]=j; al2[j]=GAP;
	 tmp_gal[gpos*2]=GAP; tmp_gal[gpos*2+1]=j-1; gpos++;
	 j--;
	}else if(dmtx[gridid].poi==UP){
	 //al1[i]=GAP; al2[j]=i;
	 tmp_gal[gpos*2]=i-1; tmp_gal[gpos*2+1]=GAP; gpos++;
	 i--;
	}
 }
 *glen=gpos;
 //put to gal
 for(int i=0;i<gpos;i++){
  gal[i*2  ] = tmp_gal[2*(gpos-i-1)];
  gal[i*2+1] = tmp_gal[2*(gpos-i-1)+1];
 }
 //init al1
 for(int i=0;i<n1;i++)
  al1[i]=-1;
 //put al1
 for(int i=0;i<gpos;i++){
  if(gal[2*i]!=-1)
   al1[gal[2*i]]=gal[2*i+1];
 }
 //init al2
 for(int i=0;i<n2;i++)
  al2[i]=-1;
 //put al1
 for(int i=0;i<gpos;i++){
  if(gal[2*i+1]!=-1)
   al2[gal[2*i+1]]=gal[2*i];
 }
 return highest;
}


//Local alignment
float dp_local(DPMTX *dmtx,float *Smtx,float *dtbl,float GapOpen,float GapExt,int *al1,int n1,int *al2,int n2,int *gal,int *glen,int *cid,bool mode){
 int i1,i2,p1,p2;
 int N1,N2;
 float dia_sco,up_sco,left_sco;
 int gridid;
 //Warning!!
 /*
 0:0 is not corresponding with Smtx[GRID2D(n1,n2,0,0)]
 0:x and x:1 is gap data of N-terminal
 use GRID2D(N1,N2,a+1,b+1)
 */

 //ali1,n1 : path data
 //ali2,n2 : sequence

 N1=n1+1;
 N2=n2+1;

 dmtx[0].sco=0;
 dmtx[0].poi=NON;
 //init N-terminal DPMTX data
 //init
 for(i1=0;i1<=n1*n2;i1++){
  dmtx[i1].sco=0;
  dmtx[i1].poi=NON;
  dmtx[i1].pre_align1=-1;
  dmtx[i1].pre_align2=-1;
 }
 for(i1=1;i1<=n1;i1++){
  dmtx[GRID2D(N1,N2,i1,0)].sco=0;
  dmtx[GRID2D(N1,N2,i1,0)].poi=NON;
  dmtx[GRID2D(N1,N2,i1,0)].pre_align1=-1;
  dmtx[GRID2D(N1,N2,i1,0)].pre_align2=-1;
 }
 for(i2=1;i2<=n2;i2++){
  dmtx[GRID2D(N1,N2,0,i2)].sco=0;//local
  dmtx[GRID2D(N1,N2,0,i2)].poi=NON;
  dmtx[GRID2D(N1,N2,0,i2)].pre_align1=-1;
  dmtx[GRID2D(N1,N2,0,i2)].pre_align2=-1;
 }
 //puts("init");
 //fill
 for(i1=1;i1<=n1;i1++){
 	for(i2=1;i2<=n2;i2++){
	 //diagonal score
	 dia_sco=dmtx[GRID2D(N1,N2,i1-1,i2-1)].sco+Smtx[GRID2D(n1,n2,i1-1,i2-1)];
	 //Node distance loss
	 p1 = dmtx[GRID2D(N1,N2,i1-1,i2-1)].pre_align1;//LDP
	 p2 = dmtx[GRID2D(N1,N2,i1-1,i2-1)].pre_align2;//Seq
	//Ca-Ca Distance score
	 if(p1 != -1){
		if(p2>0 && cid[p2-1] == cid[i2-1]){//Same ChainID
/*//v02
	  if((dtbl[i1-1]-dtbl[p1])<3.2)
	   dia_sco+=-1.0*((dtbl[i1-1]-dtbl[p1])-3.20)*((dtbl[i1-1]-dtbl[p1])-3.20);
	  if((dtbl[i1-1]-dtbl[p1])>3.8)
	   dia_sco+=-1.0*((dtbl[i1-1]-dtbl[p1])-3.80)*((dtbl[i1-1]-dtbl[p1])-3.80);
*/
		
	  	//v03
	  	 int pos_diff=(i1-1)-p1;
		 if(pos_diff>20)
		  pos_diff=20;
	  	 float dist_diff=dtbl[21*(i1-1)+pos_diff];
	  	 if(dist_diff<3.2)
	  	  dia_sco+=-1.0*(dist_diff-3.20)*(dist_diff-3.20);
	  	 if(dist_diff>3.8)
	   	  dia_sco+=-1.0*(dist_diff-3.80)*(dist_diff-3.80);
		 
		}
	 }
	 //gap score
	 //up score
	 up_sco=dmtx[GRID2D(N1,N2,i1-1,i2)].sco;
/*
	 //gap open? or extention?
	 //open
	 if(dmtx[GRID2D(N1,N2,i1-1,i2)].poi==DIA){
	  up_sco-=GapOpen;
	 }else{
	  up_sco-=GapExt;
	 }
*/
	 //left score
	 left_sco=dmtx[GRID2D(N1,N2,i1,i2-1)].sco;
	 //gap open? or extention?

	 //open
	 if(dmtx[GRID2D(N1,N2,i1,i2-1)].poi==DIA){
	  left_sco+=GapOpen;
	 }else{
	  left_sco+=GapExt;
	 }
	 //choose max score
	 if(dia_sco>=up_sco){
		if(dia_sco>=left_sco){
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=DIA;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=dia_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=i1-1;//aligned LDP
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=i2-1;//aligned Seq
		}else{
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=LEFT;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=left_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align2;
		}
	 }else{
		if(up_sco>=left_sco){
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=UP;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=up_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=dmtx[GRID2D(N1,N2,i1-1,i2)].pre_align1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=dmtx[GRID2D(N1,N2,i1-1,i2)].pre_align2;
		}else{
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=LEFT;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=left_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align2;
		}
	 }

	 //printf("%d:%d dia:%.3f up%.3f left:%.3f sco:%.3f poi %d\n",i1,i2,dia_sco,up_sco,left_sco,dmtx[GRID2D(N1,N2,i1,i2)].sco,dmtx[GRID2D(N1,N2,i1,i2)].poi);

	}
 }
 //trace back

 //start from highest score bottom or end of columms
 //end zero gaps
 float highest=0;
 int i=0,j=0;

 for(i1=1;i1<=n1;i1++){
  for(i2=1;i2<=n2;i2++){
   if(highest < dmtx[GRID2D(N1,N2,i1,i2)].sco){
    highest=dmtx[GRID2D(N1,N2,i1,i2)].sco;
    i=i1;
    j=i2;
   }
  }
 }

 //printf("#START FROM %d %d sco=%.3f %f\n",i,j,highest,dmtx[GRID2D(N1,N2,i,j)].sco);
 //printf("START FROM %d %d sco=%.3f\n",n1,n2,dmtx[GRID2D(N1,N2,n1-1,n2-1)].sco);
 //puts("#TRACE BACK");
 if(mode==false)//no alignment
  return highest;
 int gpos=0;//reverse position for global alignment
 int tmp_gal[20000];
 //Becareful!! i&j is start from 1
 while(1){
 	gridid=GRID2D(N1,N2,i,j);
	//printf("%d-%d %.3f\n",i,j,dmtx[GRID2D(N1,N2,i,j)].sco);
	if(dmtx[gridid].poi==NON)//end when 0:x or x:0
	 break;
	if(dmtx[gridid].sco<=0.00)//end when 0:x or x:0
	 break;
	if(dmtx[gridid].poi==DIA){
	 //al1[i]=j; al2[j]=i;
	 tmp_gal[gpos*2]=i-1; tmp_gal[gpos*2+1]=j-1; gpos++;
	 i--;j--;
	}else if(dmtx[gridid].poi==LEFT){
	 //al1[i]=j; al2[j]=GAP;
	 tmp_gal[gpos*2]=GAP; tmp_gal[gpos*2+1]=j-1; gpos++;
	 j--;
	}else if(dmtx[gridid].poi==UP){
	 //al1[i]=GAP; al2[j]=i;
	 tmp_gal[gpos*2]=i-1; tmp_gal[gpos*2+1]=GAP; gpos++;
	 i--;
	}
 }
 //printf("gpos= %d\n",gpos);
 *glen=gpos;
 //put to gal
 int Ngpos=0;
 for(int i=0;i<gpos;i++){
  gal[i*2  ] = tmp_gal[2*(gpos-i-1)];
  gal[i*2+1] = tmp_gal[2*(gpos-i-1)+1];
  if(gal[i*2  ]!=-1 && gal[i*2+1]!=-1)
   Ngpos++;
 }
 //Ignore too short alignment
 if(Ngpos<9)
  return(-1.00);
 //init al1
 for(int i=0;i<n1;i++)
  al1[i]=-1;
 //put al1
 //printf("put ali1 %d\n",n1);
 for(int i=0;i<gpos;i++){
  //printf("put ali1-2 %d %d\n",gal[2*i],gal[2*i+1]);
  if(gal[2*i]!=-1)
   al1[gal[2*i]]=gal[2*i+1];
 }
 //init al2
 //printf("put ali2 %d\n",n2);
 for(int i=0;i<n2;i++)
  al2[i]=-1;
 //put al1
 //printf("put ali2-2 %d\n",n2);
 for(int i=0;i<gpos;i++){
  if(gal[2*i+1]!=-1)
   al2[gal[2*i+1]]=gal[2*i];
 }
 //printf("Done\n");
 return highest;
}


float dp_global(DPMTX *dmtx,float *Smtx,float *dtbl,float GapOpen,float GapExt,int *al1,int n1,int *al2,int n2,int *gal,int *glen,bool mode){
 int i1,i2,p1;
 int N1,N2;
 float dia_sco,up_sco,left_sco;
 int gridid;
 //Warning!!
 /*
 0:0 is not corresponding with Smtx[GRID2D(n1,n2,0,0)]
 0:x and x:1 is gap data of N-terminal
 use GRID2D(N1,N2,a+1,b+1)
 */


 N1=n1+1;
 N2=n2+1;

 dmtx[0].sco=0;
 dmtx[0].poi=NON;
 //init N-terminal DPMTX data
 //init
 for(i1=0;i1<=n1*n2;i1++){
  dmtx[i1].sco=0;
  dmtx[i1].poi=NON;
  dmtx[i1].pre_align1=-1;
  dmtx[i1].pre_align2=-1;
 }
 for(i1=1;i1<=n1;i1++){
  dmtx[GRID2D(N1,N2,i1,0)].sco=0;
  dmtx[GRID2D(N1,N2,i1,0)].poi=NON;
  dmtx[GRID2D(N1,N2,i1,0)].pre_align1=-1;
  dmtx[GRID2D(N1,N2,i1,0)].pre_align2=-1;
 }
 for(i2=1;i2<=n2;i2++){
  dmtx[GRID2D(N1,N2,0,i2)].sco=GapOpen;//global
  dmtx[GRID2D(N1,N2,0,i2)].poi=NON;
  dmtx[GRID2D(N1,N2,0,i2)].pre_align1=-1;
  dmtx[GRID2D(N1,N2,0,i2)].pre_align2=-1;
 }
 //puts("init");
 //fill
 for(i1=1;i1<=n1;i1++){
 	for(i2=1;i2<=n2;i2++){
	 //diagonal score
	 dia_sco=dmtx[GRID2D(N1,N2,i1-1,i2-1)].sco+Smtx[GRID2D(n1,n2,i1-1,i2-1)];
	 //Node distance loss
	 p1 = dmtx[GRID2D(N1,N2,i1-1,i2-1)].pre_align1;
	 if(p1 != -1){
	  if((dtbl[i1-1]-dtbl[p1])<3.2)
	   dia_sco+=-1.0*((dtbl[i1-1]-dtbl[p1])-3.20)*((dtbl[i1-1]-dtbl[p1])-3.20);
	  if((dtbl[i1-1]-dtbl[p1])>3.8)
	   dia_sco+=-1.0*((dtbl[i1-1]-dtbl[p1])-3.80)*((dtbl[i1-1]-dtbl[p1])-3.80);
	 }
	 //gap score
	 //up score
	 up_sco=dmtx[GRID2D(N1,N2,i1-1,i2)].sco;
/*
	 //gap open? or extention?
	 //open
	 if(dmtx[GRID2D(N1,N2,i1-1,i2)].poi==DIA){
	  up_sco-=GapOpen;
	 }else{
	  up_sco-=GapExt;
	 }
*/
	 //left score
	 left_sco=dmtx[GRID2D(N1,N2,i1,i2-1)].sco;
	 //gap open? or extention?

	 //open
	 if(dmtx[GRID2D(N1,N2,i1,i2-1)].poi==DIA){
	  left_sco+=GapOpen;
	 }else{
	  left_sco+=GapExt;
	 }
	 //choose max score
	 if(dia_sco>=up_sco){
		if(dia_sco>=left_sco){
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=DIA;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=dia_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=i1-1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=i2-1;
		}else{
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=LEFT;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=left_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align2;
		}
	 }else{
		if(up_sco>=left_sco){
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=UP;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=up_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=dmtx[GRID2D(N1,N2,i1-1,i2)].pre_align1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=dmtx[GRID2D(N1,N2,i1-1,i2)].pre_align2;
		}else{
		 dmtx[GRID2D(N1,N2,i1,i2)].poi=LEFT;
		 dmtx[GRID2D(N1,N2,i1,i2)].sco=left_sco;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align1=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align1;
		 dmtx[GRID2D(N1,N2,i1,i2)].pre_align2=dmtx[GRID2D(N1,N2,i1,i2-1)].pre_align2;
		}
	 }

	 //printf("%d:%d dia:%.3f up%.3f left:%.3f sco:%.3f poi %d\n",i1,i2,dia_sco,up_sco,left_sco,dmtx[GRID2D(N1,N2,i1,i2)].sco,dmtx[GRID2D(N1,N2,i1,i2)].poi);

	}
 }
 //trace back

 //start from highest score bottom or end of columms
 //end zero gaps
 float highest=-100000;
 int i=0,j=0;
/*
 for(i1=1;i1<=n1;i1++){
  for(i2=1;i2<=n2;i2++){
   if(highest < dmtx[GRID2D(N1,N2,i1,i2)].sco){
    highest=dmtx[GRID2D(N1,N2,i1,i2)].sco;
    i=i1;
    j=i2;
   }
  }
 }
*/
//Global
/*
 for(i1=1;i1<=n1;i1++){
  if(highest < dmtx[GRID2D(N1,N2,i1,n2)].sco){
   highest=dmtx[GRID2D(N1,N2,i1,n2)].sco;
   i=i1;
   j=n2;
  }
 }
*/
 i=n1;j=n2;
 printf("#START FROM %d %d sco=%.3f %f\n",i,j,highest,dmtx[GRID2D(N1,N2,i,j)].sco);
 //printf("START FROM %d %d sco=%.3f\n",n1,n2,dmtx[GRID2D(N1,N2,n1-1,n2-1)].sco);
 //puts("#TRACE BACK");
 if(mode==false)//no alignment
  return highest;
 int gpos=0;//reverse position for global alignment
 int tmp_gal[20000];
 //Becareful!! i&j is start from 1
 while(1){
 	gridid=GRID2D(N1,N2,i,j);
	//printf("%d-%d %.3f\n",i,j,dmtx[GRID2D(N1,N2,i,j)].sco);
	if(dmtx[gridid].poi==NON)//end when 0:x or x:0
	 break;
	//if(dmtx[gridid].sco<=0.00)//end when 0:x or x:0
	// break;
	if(dmtx[gridid].poi==DIA){
	 //al1[i]=j; al2[j]=i;
	 tmp_gal[gpos*2]=i-1; tmp_gal[gpos*2+1]=j-1; gpos++;
	 i--;j--;
	}else if(dmtx[gridid].poi==LEFT){
	 //al1[i]=j; al2[j]=GAP;
	 tmp_gal[gpos*2]=GAP; tmp_gal[gpos*2+1]=j-1; gpos++;
	 j--;
	}else if(dmtx[gridid].poi==UP){
	 //al1[i]=GAP; al2[j]=i;
	 tmp_gal[gpos*2]=i-1; tmp_gal[gpos*2+1]=GAP; gpos++;
	 i--;
	}
 }
 //printf("gpos= %d\n",gpos);
 *glen=gpos;
 //put to gal
 for(int i=0;i<gpos;i++){
  gal[i*2  ] = tmp_gal[2*(gpos-i-1)];
  gal[i*2+1] = tmp_gal[2*(gpos-i-1)+1];
 }
 //init al1
 for(int i=0;i<n1;i++)
  al1[i]=-1;
 //put al1
 //printf("put ali1 %d\n",n1);
 for(int i=0;i<gpos;i++){
  //printf("put ali1-2 %d %d\n",gal[2*i],gal[2*i+1]);
  if(gal[2*i]!=-1)
   al1[gal[2*i]]=gal[2*i+1];
 }
 //init al2
 //printf("put ali2 %d\n",n2);
 for(int i=0;i<n2;i++)
  al2[i]=-1;
 //put al1
 //printf("put ali2-2 %d\n",n2);
 for(int i=0;i<gpos;i++){
  if(gal[2*i+1]!=-1)
   al2[gal[2*i+1]]=gal[2*i];
 }
 //printf("Done\n");
 return highest;
}

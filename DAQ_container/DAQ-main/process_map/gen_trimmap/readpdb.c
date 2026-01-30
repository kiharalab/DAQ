#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>
#include "struct.h"
#include "func.h"

//New ver 2012.10.29 by terashi
//Custmized for marching cubes!!!

int Atom2Int(char);
int Atom2Class(char *);

int readpdb(PDB *pdb, char *filename, int n){

	int i,j; 
	FILE *fpin; 
	char line[LIN], buf[LIN]; 
	double x,y,z;
	double x_min,y_min,z_min;
	double x_max,y_max,z_max;
	x_min=9999999999999;
	y_min=9999999999999;
	z_min=9999999999999;
	x_max=-9999999999999;
	y_max=-9999999999999;
	z_max=-9999999999999;
	//dock
	int no1=0,no2=0;
	int real;
	COORD sum_coord;
   	int sum_num=0,nmr=0;
	x=y=z;

	pdb->NumOfAtom=n;
	
	if((fpin=fopen(filename,"r")) == NULL){ 
		fprintf(stderr,"Can't open %s\n",filename); 
		return(-1); 
	}

	//printf("#filename: %s\n",filename); 
	i=0;j=-1;
	int Nca=0;
	int PreResNum=-99999;
	int ResNum=-1;
	pdb->NumOfReal=0;
	while(fgets(line,LIN,fpin)){ 
		if(!strncmp(line,"ATOM",4) 
		 &&(!strncmp(&line[16],"A",1) 
		 ||!strncmp(&line[16]," ",1))){

			//Igonore hydrogen New For PRESCO!!
			if(!strncmp(&line[13],"A",1)||!strncmp(&line[12],"A",1))
			 continue;
			//-------
			//Res num order
			strncpy(buf,&line[22],4);
			buf[4]='\0';
			ResNum=atoi(buf);
			//For each residue========================
			if(PreResNum != ResNum){
			 j++;
			 PreResNum=ResNum;
			 pdb->ResNum[j]=ResNum;//Real number
			 //Residue Type
			 strncpy(pdb->TypeRes[j],&line[17],3);
			 pdb->TypeRes[j][3]='\0';
			 pdb->TypeResId[j]=AA2int(pdb->TypeRes[j]);
			 pdb->ResOnAtom[j]=i;
			}
			//========================================

			//For each Atoms
			pdb->AtomOnRes[i]=j;//Order base

			//chain id
			pdb->Chain[i]=line[21];

			//Real Atom number
			strncpy(buf,&line[4],7); buf[7]='\0';
			pdb->AtomNum[i]=atoi(buf);
			strncpy(pdb->TypeAtom[i],&line[13],3);
			pdb->TypeAtom[i][3]='\0'; 


			//pdb->TypeAtomId[i]=Atom2Int(pdb->TypeAtom[i][0]);			
			pdb->TypeAtomId[i]=Atom2Class(pdb->TypeAtom[i]);			

			//xyz coords
			strncpy(buf,&line[30],8); buf[8]='\0'; 
			x=atof(buf);
			strncpy(buf,&line[38],8); buf[8]='\0'; 
			y=atof(buf);
			strncpy(buf,&line[46],8); buf[8]='\0'; 
			z=atof(buf);
			if(x<x_min)x_min=x;
			if(y<y_min)y_min=y;
			if(z<z_min)z_min=z;
			if(x>x_max)x_max=x;
			if(y>y_max)y_max=y;
			if(z>z_max)z_max=z;

			pdb->xyz[i][0]=x;
			pdb->xyz[i][1]=y;
			pdb->xyz[i][2]=z;

			if(pdb->TypeAtomId[i]==2){//CA
			 pdb->CAxyz[j][0]=x;
			 pdb->CAxyz[j][1]=y;
			 pdb->CAxyz[j][2]=z;
			}


			//Max and Min
			if(pdb->MaxXyz[0] < x||i==0) pdb->MaxXyz[0]=x;
			if(pdb->MaxXyz[1] < y||i==0) pdb->MaxXyz[1]=y;
			if(pdb->MaxXyz[2] < z||i==0) pdb->MaxXyz[2]=z;
			if(pdb->MinXyz[0] > x||i==0) pdb->MinXyz[0]=x;
			if(pdb->MinXyz[1] > y||i==0) pdb->MinXyz[1]=y;
			if(pdb->MinXyz[2] > z||i==0) pdb->MinXyz[2]=z;

			i++;
			if(i>n)
			 return -1;
		}else if(!strncmp(line,"MODEL",5)){//NMR
    	 	 nmr++;
		 if(nmr>1)//first model only
		  break; 
		}

	}
	pdb->ResOnAtom[j+1]=i;
	pdb->AtomOnRes[i]=j+1;
	//pdb->ResNnum[j+1]=i; 
	fclose(fpin) ; 
	pdb->NumOfAtom=i; 
	pdb->NumOfRes=j+1;//New!!
    printf("#PDB range: X:[%.5f,%.5f] Y:[%.5f,%.5f] Z:[%.5f,%.5f]\n",x_min,x_max,y_min,y_max,z_min,z_max);
	

	return(0);
}


int MallocPdb(PDB *p,int n){
 int i;
	//malloc here!!
	//xyz coord data double n*3 format
	if((p->xyz=(double **)malloc(sizeof(double *)*n))==NULL){
	 free(p->xyz);
	 return -1;
	}
	if((p->CAxyz=(double **)malloc(sizeof(double *)*n))==NULL){
	 free(p->CAxyz);
	 return -1;
	}
	for(i=0;i<n;i++){
	 if((p->xyz[i]=(double *)malloc(sizeof(double)*3))==NULL){
	  free(p->xyz);
	  return -1;
	 }
	}
	for(i=0;i<n;i++){
	 if((p->CAxyz[i]=(double *)malloc(sizeof(double)*3))==NULL){
	  free(p->CAxyz);
	  return -1;
	 }
	}
	//Residue number on each Atom
	if((p->AtomOnRes=(int *)malloc(sizeof(int *)*n+1))==NULL){
	 free(p->AtomOnRes);
	 return -1;
	}
	//Atom number on each residue
	if((p->ResOnAtom=(int *)malloc(sizeof(int *)*n+1))==NULL){
	 free(p->ResOnAtom);
	 return -1;
	}
	//Type of Atom 3char
	if((p->TypeAtom=(char **)malloc(sizeof(char*)*n))==NULL){
	 return -1;
	}
	for(i=0;i<n;i++){
	 if((p->TypeAtom[i]=(char *)malloc(sizeof(char)*4))==NULL){
	  return -1;
	 }
	}
	if((p->TypeAtomId=(int *)malloc(sizeof(int)*n))==NULL){
	 return -1;
	}
	
	//Other data
	if((p->TypeAtomId=(int*)malloc(sizeof(int)*n))==NULL){
	 free(p->TypeAtomId);
	 return -1;
	}
	if((p->TypeResId=(int*)malloc(sizeof(int)*n))==NULL){
	 free(p->TypeResId);
	 return -1;
	}

	//Residue data
	if((p->TypeRes=(char **)malloc(sizeof(char*)*n))==NULL){
	 return -1;
	}
	for(i=0;i<n;i++){
	 if((p->TypeRes[i]=(char *)malloc(sizeof(char)*4))==NULL){
	  return -1;
	 }
	}
	//Real Res number
	if((p->ResNum=(int *)malloc(sizeof(int)*n))==NULL){
	 return -1;
	}
	//Real Atom Number
	if((p->AtomNum=(int *)malloc(sizeof(int)*n))==NULL){
	 return -1;
	}
	if((p->Chain=(char *)malloc(sizeof(char)*n))==NULL){
	 return -1;
	}

	//Depth
	if((p->DepthAtom=(float *)malloc(sizeof(float)*n))==NULL){
	 free(p->DepthAtom);
	 return -1;
	}
	if((p->DepthRes=(float *)malloc(sizeof(float)*n))==NULL){
	 free(p->DepthRes);
	 return -1;
	}
	
 //puts("#fin malloc");
 return 0;
}


int CountAtom(char *filename){
	int i,j; 
	FILE *fpin; 
	char line[LIN], buf[LIN]; 
	double x,y,z;
	//dock
   	int nmr=0;
	x=y=z;
	
	if((fpin=fopen(filename,"r")) == NULL){ 
		fprintf(stderr,"Can't open %s\n",filename); 
		return(-1); 
	}

	//printf("#filename: %s\n",filename); 
	i=0;j=-1;
	int Natom=0;
	while(fgets(line,LIN,fpin)){ 
		if(!strncmp(line,"ATOM",4) 
		 &&(!strncmp(&line[16],"A",1) 
		 ||!strncmp(&line[16]," ",1))){ 
		 Natom++;
		}else if(!strncmp(line,"MODEL",5)){//NMR
    	 	 nmr++;
		 if(nmr>1)//first model only
		  break; 
		}
	}
 return Natom;
}



//side chain & asa det
int side_det(PDB *pdb,int num){
 if(!strncmp(pdb->TypeAtom[num],"N  ",3))
  return(1);
 if(!strncmp(pdb->TypeAtom[num],"CA ",3)
  && strncmp(pdb->TypeRes[pdb->AtomOnRes[num]],"GLY",3))
  return(1);
 if(!strncmp(pdb->TypeAtom[num],"C  ",3))
  return(1);
 if(!strncmp(pdb->TypeAtom[num],"O  ",3))
  return(1);
 if(!strncmp(pdb->TypeAtom[num],"OXT",3))
  return(1);

 return(0);
}

int AA2int(char *TypeAtom){
        if(!strcmp(TypeAtom,"ALA")){return 0;}
        if(!strcmp(TypeAtom,"VAL")){return 1;}
        if(!strcmp(TypeAtom,"PHE")){return 2;}
        if(!strcmp(TypeAtom,"PRO")){return 3;}
        if(!strcmp(TypeAtom,"MET")){return 4;}
        if(!strcmp(TypeAtom,"ILE")){return 5;}
        if(!strcmp(TypeAtom,"LEU")){return 6;}
        if(!strcmp(TypeAtom,"ASP")){return 7;}
        if(!strcmp(TypeAtom,"GLU")){return 8;}
        if(!strcmp(TypeAtom,"LYS")){return 9;}
        if(!strcmp(TypeAtom,"ARG")){return 10;}
        if(!strcmp(TypeAtom,"SER")){return 11;}
        if(!strcmp(TypeAtom,"THR")){return 12;}
        if(!strcmp(TypeAtom,"TYR")){return 13;}
        if(!strcmp(TypeAtom,"HIS")){return 14;}
        if(!strcmp(TypeAtom,"CYS")){return 15;}
        if(!strcmp(TypeAtom,"ASN")){return 16;}
        if(!strcmp(TypeAtom,"TRP")){return 17;}
        if(!strcmp(TypeAtom,"GLN")){return 18;}
        if(!strcmp(TypeAtom,"GLY")){return 19;}
        return -1;
}

int Atom2Class(char *atom){
 if(!strncmp(atom,"N  ",3))
  return(1);
 if(!strncmp(atom,"CA ",3))
  return(2);
 if(!strncmp(atom,"C  ",3))
  return(3);
 if(!strncmp(atom,"O  ",3))
  return(4);
 if(!strncmp(atom,"CB ",3))
  return(5);
 //Other
 return(0);
}




int A2int(char TypeAtom){
        switch(TypeAtom){
                case 'A':return 0;break;
                case 'V':return 1;break;
                case 'F':return 2;break;
                case 'P':return 3;break;
                case 'M':return 4;break;
                case 'I':return 5;break;
                case 'L':return 6;break;
                case 'D':return 7;break;
                case 'E':return 8;break;
                case 'K':return 9;break;
                case 'R':return 10;break;
                case 'S':return 11;break;
                case 'T':return 12;break;
                case 'Y':return 13;break;
                case 'H':return 14;break;
                case 'C':return 15;break;
                case 'N':return 16;break;
                case 'W':return 17;break;
                case 'Q':return 18;break;
                case 'G':return 19;break;
                case '-':return 20;break;
                default: return -1;
        }
}

int Atom2Int(char TypeAtom){
	switch(TypeAtom){
                case 'N':return 0;break;
                case 'C':return 1;break;
                case 'O':return 2;break;
                case 'S':return 3;break;
                case 'P':return 4;break;
                case 'H':return 5;break;
                default: return 6;
        }
}


int ReadDssp(DSSP *d, char *filename){
	int i,j,n; 
	FILE *fpin; 
	char line[LIN], buf[LIN]; 
	double x,y,z;

	if((fpin=fopen(filename,"r")) == NULL){ 
		fprintf(stderr,"Can't open %s\n",filename); 
		return(-1); 
	}

	printf("#DSSP filename: %s\n",filename); 
	i=0;j=-1;
	int r=0;
	bool flag=false;
	char ss;
	while(fgets(line,LIN,fpin)){ 
		i++;
		if(!strncmp(&line[18],"TOTAL NUMBER",12)){//Size
		 strncpy(buf,&line[0],5); buf[5]='\0'; 
		 n=atoi(buf);
		 d->Nres=n;
		 printf("#NresDSSP= %d\n",n);
		 //malloc
		 if((d->kappa=(double *)malloc(sizeof(double)*n))==NULL)
		  return -1;
		 if((d->alpha=(double *)malloc(sizeof(double)*n))==NULL)
		  return -1;
		 if((d->phi=(double *)malloc(sizeof(double)*n))==NULL)
		  return -1;
		 if((d->psi=(double *)malloc(sizeof(double)*n))==NULL)
		  return -1;
		 if((d->acc=(int *)malloc(sizeof(int)*n))==NULL)
		  return -1;
		 if((d->ss=(int *)malloc(sizeof(int)*n))==NULL)
		  return -1;
		 if((d->xyz=(double **)malloc(sizeof(double*)*n))==NULL)
		  return -1;
		 for(int p=0;p<n;p++)
			d->xyz[p]=(double *)malloc(sizeof(double)*3);
		}
		if(!strncmp(line,"  #",3)){
		 flag=true;
		 continue;
		}
		if(flag){
		 //printf("%s\n",line);
		 	if(!strncmp(&line[13],"!",1)){
				puts("#skip");
				continue;
			}
		 	//ss
		 	switch(line[16]){
                	 	case 'H':d->ss[r]=0;break;
                	 	case 'E':d->ss[r]=1;break;
                		default: d->ss[r]=2;
        		}

			//acc
			strncpy(buf,&line[35],3); buf[3]='\0';
			d->acc[r]=atoi(buf);
			//kappa
			strncpy(buf,&line[91],6); buf[6]='\0';
			d->kappa[r]=atof(buf);
			//alpha
			strncpy(buf,&line[97],6); buf[6]='\0';
			d->alpha[r]=atof(buf);
			//phi
			strncpy(buf,&line[103],6); buf[6]='\0';
			d->phi[r]=atof(buf);
			//psi
			strncpy(buf,&line[109],6); buf[6]='\0';
			d->psi[r]=atof(buf);
			
			//x
			strncpy(buf,&line[115],7); buf[7]='\0';
			d->xyz[r][0]=atof(buf);
			//y
			strncpy(buf,&line[122],7); buf[7]='\0';
			d->xyz[r][1]=atof(buf);
			//z
			strncpy(buf,&line[129],7); buf[7]='\0';
			d->xyz[r][2]=atof(buf);
			printf("#%s\n",line);
			printf("#%d ss= %d acc= %d",r,d->ss[r],d->acc[r]);
			printf(" kappa= %.1f acc= %.1f",d->kappa[r],d->alpha[r]);
			printf(" phi= %.1f psi= %.1f",d->phi[r],d->psi[r]);
			printf(" xyz= %.1f %.1f %.1f \n",d->xyz[r][0],d->xyz[r][1],d->xyz[r][2]);
		 r++;
		}
	}
 printf("#Reading Dssp Done\n");
 return 0;
}

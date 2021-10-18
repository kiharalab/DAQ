#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>
#include <math.h>
#include "struct.h"


int chkcmdline( int argc, char **argv,CMD *cmd){
        char **p;
	int num=0;
	void errmsg();
        if(argc < 2){
		errmsg();
                return(FALSE);
        }
	//default values
        p=argv;
	cmd->map_t=0.00;
	cmd->Nthr=2;
	cmd->dreso=16.00;
	cmd->MergeDist=0.50;
	cmd->Filter=0.10;
	cmd->LocalR=10.0;
	cmd->Dkeep=0.5;
	cmd->Nround=5000;
	cmd->Nnb=30;
	cmd->Ntabu=100;
	cmd->Nsim=10;
	cmd->Allow=1.01;
	cmd->Nbeam=20;
	cmd->ssize=2.0;
	cmd->ang=15.0;
	cmd->TopN=10;
	cmd->ShowGrid=false;

	cmd->th1=0.00;
	cmd->th2=0.00;

	cmd->Mode=0;//No PDB mode: 0, PDB mode : 1
	cmd->Amode=1;
	cmd->Emode=false;
	cmd->IgnoreNC=false;
	cmd->slide=1;
	cmd->NoP1=5.0;
	cmd->width=1.0;

	cmd->r=2.0;
	cmd->Ubound=0.95;
	cmd->Lbound=0.001;
	cmd->Nvox=5;
        while (--argc){
         p++;
         if(**p == '-'){
          	switch(*(*p+1)){
		 case 'i': //MAP
			strcpy(cmd->filename,*(++p));
                	--argc; break;
		 case 'p': //PDB
			strcpy(cmd->file1,*(++p));
			cmd->Mode=1;
                	--argc; break;
		 case 'd':
			strcpy(cmd->file2,*(++p));
                	--argc; break;
		 case 'c':
			cmd->Nthr=atoi(*(++p)); 
			--argc; break;
		 case 't':
			cmd->th1=atof(*(++p)); 
			--argc; break;
		 case 'r':
			cmd->r=atof(*(++p)); 
			--argc; break;
		 case 'g':
                        cmd->IgnoreNC=true;
                        break;
		 case 'v':
                        cmd->Nvox=atoi(*(++p));
                        --argc; break;
		 case 'R':
                        cmd->LocalR=atof(*(++p));
                        --argc; break;
		 case 'k':
                        cmd->Dkeep=atof(*(++p));
                        --argc; break;
		 case 'N':
                        cmd->NoP1=atoi(*(++p));
                        --argc; break;
		 case 's':
                        cmd->slide=atoi(*(++p));
                        --argc; break;
		 case 'w':
                        cmd->width=atof(*(++p));
                        --argc; break;
		 case 'U':
                        cmd->Ubound=atof(*(++p));
                        --argc; break;
		 case 'L':
                        cmd->Lbound=atof(*(++p));
                        --argc; break;
		 case 'S':
                        cmd->ShowGrid=true;
                        break;
		 case 'V':
                        cmd->Mode=1;
                        break;
		 case 'C':
                        cmd->Mode=3;
                        break;
		 case 'P':
                        cmd->Mode=4;
                        break;
		 case 'F':
                        cmd->Mode=5;
                        break;
		 case 'E':
                        cmd->Emode=true;
                        break;
		 default: 
		  	fprintf(stderr,"No such a option: %s\n",*p+1); 
			errmsg(); 
			return(FALSE); 
		 break;
	  }
	 }
        }
	cmd->NoP2=cmd->NoP1+cmd->width;
	//option check
	printf("#Map Threshold= %f\n",cmd->th1);
	printf("#Density Upperbound Rate= %f\n",cmd->Ubound);
	printf("#Voxel Size = %d^3 (2*%d+1)^3\n",2*cmd->Nvox+1,cmd->Nvox);
	printf("#Sliding Step = %d\n",cmd->slide);
	printf("#Protein  d <= %f\n",cmd->r);
	printf("#No Protein %f < d <= %f\n",cmd->NoP1,cmd->NoP2);
	printf("mode: %d\n",cmd->Mode);
        return(TRUE);
}

void errmsg(){
	puts("Usage: TrimMapAtom -i [MAP.mrc] -p [INPUT.pdb] -d [INPUT.dssp] [(option)]");
	puts("v0.10	Start");
	puts("v0.60	Add Area(log(frequency)) cut-off");
	puts("---Options---");
	//printf("-c [int  ] : Number of cores for threads def=%d\n",2);
	printf("-t [float] : Threshold of density map1 def=%.3f\n",0.00);
	printf("-v [int  ] : Voxcel Size (n*2+1)^3 def=%d\n",5);
	printf("-r [float] : Distance Cut-off d <= r def=2.0\n");
	printf("-s [int  ] : Sliding Step size def=1 (all voxels)\n");
	printf("-g         : Ignore NCSTART data def=false\n");
	printf("-N [float] : No protein radius r < d <= r+w def=5.0\n");
	printf("-w [float] : Width of No protein radius w def=1.0\n");
	printf("-U [float] : Upper-bound of the Area of log(frequency) def=0.95\n");
	printf("-L [float] : Lower-bound of normalized density value. def=0.001\n");
	printf("Thi is Ver %.3f\n",VER);
}

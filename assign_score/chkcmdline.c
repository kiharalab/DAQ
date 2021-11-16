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
	cmd->dreso=2.00;
	//cmd->MergeDist=0.50;
	cmd->MergeDist=1.00;//2021
	cmd->Filter=0.10;
	cmd->Mode=0;
	cmd->LocalR=10.0;
	cmd->Dkeep=0.5;
	cmd->Nround=5000;
	cmd->Nnb=30;
	cmd->Ntabu=100;
	cmd->Nsim=10;
	cmd->Allow=1.01;
	cmd->Pcut=0.50;
	cmd->Nbeam=20;
	cmd->Wp=1.0;
	cmd->Wd=1.0;
	cmd->FragLen = 9;
	cmd->Wss=cmd->Waa=cmd->Watm=1.0;
	cmd->LowestSco=20.0;
	cmd->Cmode=0;
        while (--argc){
         p++;
         if(**p == '-'){
          	switch(*(*p+1)){
		 case 'i':
			strcpy(cmd->filename,*(++p));
                	--argc; break;
		 case 'p':
			strcpy(cmd->pfilename,*(++p));
                	--argc; break;
		 case 's':
			strcpy(cmd->sfilename,*(++p));
                	--argc; break;
		 case 'D':
			strcpy(cmd->mfilename,*(++p));
			cmd->Mode = 5;
                	--argc; break;
		 case 'Q':
			strcpy(cmd->mfilename,*(++p));
			cmd->Mode = 6;
                	--argc; break;
		 case 'c':
			cmd->Nthr=atoi(*(++p)); 
			--argc; break;
		 case 't':
			cmd->map_t=atof(*(++p)); 
			--argc; break;
		 case 'g':
			cmd->dreso=atof(*(++p)); 
			--argc; break;
		 case 'f':
                        cmd->Filter=atof(*(++p));
                        --argc; break;
		 case 'm':
                        cmd->MergeDist=atof(*(++p));
                        --argc; break;
		 case 'R':
                        cmd->LocalR=atof(*(++p));
                        --argc; break;
		 case 'k':
                        cmd->Dkeep=atof(*(++p));
                        --argc; break;
		 case 'w':
                        //cmd->Wp=atof(*(++p));
                        cmd->Wd=atof(*(++p));
                        --argc; break;
		 case 'W':
                        sscanf(*(++p),"%f,%f,%f",&(cmd->Wss),&(cmd->Waa),&(cmd->Watm));
                        --argc; break;
		 case 'P':
                        cmd->Pcut=atof(*(++p));
                        --argc; break;
		 case 'r':
                        cmd->Nround=atoi(*(++p));
                        --argc; break;
		 case 'b':
                        cmd->Nnb=atoi(*(++p));
                        --argc; break;
		 case 'l':
                        //cmd->FragLen=atoi(*(++p));
                        cmd->LowestSco=atof(*(++p));
                        --argc; break;
		 case 'n':
                        cmd->Nsim=atoi(*(++p));
                        --argc; break;
		 case 'a':
                        cmd->Allow=atoi(*(++p));
                        --argc; break;
		 case 'V': //Visualize mode
                        cmd->Mode=4;
                        break;
		 case 'T':
                        cmd->Mode=0;
                        break;
		 case 'G':
                        cmd->Mode=1;
                        break;
		 case 'M':
                        cmd->Mode=2;
                        break;
		 case 'L':
                        cmd->Mode=3;
                        break;
		 case 'A':
                        cmd->Cmode=3;
                        break;
		 case 'O':
                        cmd->Cmode=2;
                        break;
		 case 'S':
                        cmd->Cmode=1;
                        break;
		 default: 
		  	fprintf(stderr,"No such a option: %s\n",*p+1); 
			errmsg(); 
			return(FALSE); 
		 break;
	  }
	 }
        }
	//option check
	printf("#Number of threads= %d\n",cmd->Nthr);
	printf("#Map Threshold= %f\n",cmd->map_t);
	printf("#Band Width= %f\n",cmd->dreso);
	printf("#Merge Dist= %f\n",cmd->MergeDist);
	printf("#Filtering Dens Cut= %f\n",cmd->Filter);
	printf("#Mode= %d 0:Tabu, 1:Graph, 2: MST\n",cmd->Mode);
	printf("#Backbone Pcut= %f\n",cmd->Pcut);
	printf("#Weights SS, AA, ATOM = %f,%f,%f\n",cmd->Wss,cmd->Waa,cmd->Watm);
	printf("#Weight of Density = %f\n",cmd->Wd);
	printf("#PTABLE: %s\n",cmd->pfilename);
	printf("#SEQ: %s\n",cmd->sfilename);
	printf("#Lowest Score of Fragment: %f\n",cmd->LowestSco);
	
        return(TRUE);
}

void errmsg(){
	//puts("Usage: DAQscore -i [EM map] -p [Probability table] -s [sequence file or *.spd] [(option)]");
	puts("Usage: DAQscore -i [EM map] -p [Probability table] -Q [PDB file] [(option)]");
	puts("Fast Mainmast C-lang&multi thread version");
	puts("v0.1	Add Suboptimal Alignments");
	puts("v0.2	Add Threading for Hetero-oligomer");
	puts("v0.3	Add Minimum Score of Fragments");
	puts("v0.3	Set Minimum Length to 9");
	puts("v0.3	Add CA-CA distance table with 20 window size");
	puts("v0.4	Add MQA mode");
	puts("v0.44	Add Dssp, spot-1d-single Mode for -s input");
	puts("---Mode---");
/*
	puts("-L : Fast LDP search mode");
	puts("-M : Minimum Spanning Tree Mode");
	puts("-G : Graph Mode");
	puts("-T : Tabu Search Mode (Default)");
	puts("-V : Visualize PTBL. EMmap and Probability table are required.");
	puts("-Q [ PDB format file ] : Model Quality Assessment mode");
*/
	puts("---Options---");
	printf("-c [int  ] :Number of cores for threads def=%d\n",2);
	/*
	printf("-t [float] :Threshold of the Probability map def=%.3f\n",0.00);
	printf("-g [float] : bandwidth of the gaussian filter\n");
        printf("             def=2.0, sigma = 0.5*[float]\n");
	printf("-f [float] :Filtering for representative points def=%.3f\n",0.10);
	printf("-m [float] :After MeanShifting merge<[float] def=%.3f\n",1.00);
	printf("-R [float] :Radius of Local MST def=%.3f\n",10.0);
	printf("-k [float] :keep edges where d<[float] def=%.3f\n",0.5);
	printf("-w [float] :Weight of Density def=%.3f\n",1.0);
	printf("-P [float] :Probability(N,CA,C)>P_CutOff def=%.3f\n",0.5);
	printf("-l [float] :Lowest Score Cut-off of Fragment Length def=%f\n",20.0);
	printf("-W [float] :Weight for SS,AA,ATOM def=%f,%f,%f\n\tExample: -W 0.1,0.1,1.0\n",1.0,1.0,1.0);
	*/
/*
	puts("---Options for QA mdoe---");
	printf("-A : Amino Acid Score Only\n");
	printf("-O : atOm Score Only\n");
	printf("-S : SS Score Only\n");
*/	

	printf("Thi is Ver %.3f\n",VER);
}

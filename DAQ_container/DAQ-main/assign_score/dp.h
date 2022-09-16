
#define GRID(N,a,b,c) c+N*(b+N*a)
#define GRID3D(N1,N2,a,b,c) c+N2*(b+N1*a)
#define GRID2D(N1,N2,a,b) (b)+(N2)*(a)


typedef struct{ 
	float sco; 
	short int poi;
	short int pre_align1,pre_align2;
} DPMTX;

/*
typedef struct{  
        float sco;
	DPMTX *Dmtx;
	float *Smtx,*SmtxRv, *dtbl;
	int n1,n2,Lgali;
	int *ali1,*ali2,*gali;
	NODE *nodes;
} DP_MEMORY;
*/

float dp(DPMTX *,float *,float ,float ,int *,int,int *,int ,int *,int *);
float dp_fast(DPMTX *,float *,float *,float ,float ,int *,int,int *,int ,int *,int *,bool);
float dp_local(DPMTX *,float *,float *,float ,float ,int *,int,int *,int ,int *,int *,int *,bool);
float dp_global(DPMTX *,float *,float *,float ,float ,int *,int,int *,int ,int *,int *,bool);


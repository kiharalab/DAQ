#define VER 0.44
#define PI 3.141592
#define ATOM 50000
#define RES 5000
#define AMI 20
#define LIN 256
#define radian 57.29577951308232088

//grid
#define GRINO 20000


/*def の環境設定*/
#define INC 15/*ligand回転の度数*/
#define MOV 1/*ligand移動の度数*/
#define START_RT 0 /*ligand回転のstart*/
#define START_MV 0 /*ligand移動のstart*/
#define FIN_RT 360
#define FIN_MV 0
#define ON 1
#define OFF 0
#define TRUE 0
#define FALSE -1
/*#define MATFILE "result030106"*/ /*スコアマトリクスのファイル*/
#define MATSIZE 1500
//#define MATLEN 600
#define MATAMI 24 /*マトリクスのアミノ酸種類*/

#define CSHNUM 20 /*衝突残基数の制限*/

/*マクロ設定*/
#define X(a) (a)*(a)
#define L(a,b,c,d,e,f) (sqrt((d-a)*(d-a)+(e-b)*(e-b)+(f-c)*(f-c)))/*ベクトル長*/
#define RAS(a) (2.000000*PI*a/360.000000)/*度->ラジアン*/
#define RAD(a) (2.000000*PI*a/360.000000)/*度->ラジアン*/
#define ANG(a) (360.000000*a/(2.000000*PI))/*ラジアン->度*/

#define GAUSS(a,b) 1.00/(sqrt(2.00*PI)*b)*exp(-(a*a)/(2.0*b*b))


/*triangle*/

#define MAXMTX   4

/*結果出力*/
#define TOP 3000

#define NOT_GAUSS 1
#define USE_GAUSS 0

#define VALIABLE 0
#define CONSTANT 1

//int FINAL_ATNo;
//int FINAL_RENo;
typedef struct{
	        //double x,y,z;
	        float x,y,z;
}COORD;


typedef struct{
        char fname[LIN];
        int NumOfAtom, NumOfRes; 
        //int ResNnum[RES];
        //int SS[RES];
        //int AA2int_data[RES];
        //int AA2int_data_real[RES];
        //New!! from sakai typeatm.h
        int TypeAtomOder[RES][17];
        //int AtomOnRes[ATOM], SosuiAtom[ATOM], ConservedAtom[ATOM]; 
        //char TypeAtom[ATOM][4], TypeRes[RES][4], Chain[ATOM][2],RealNum[RES][5]; 
        float *Charge; 
        COORD *coord;
        COORD CAcd[RES];
        COORD *CBcd;
        COORD *Cen;
        COORD *Intra;//interaction
        int NumOfIntra;
        COORD *Nonin;//non intra
        int NumOfNonin;
        float *phi,*psi;

        //HETATM
        char **HET_TypeAtom;
        COORD *HET_coord;
        int NumOfHet;
        int NumOfReal;
        int RealResNum[RES];
        //New 2012.10.29-----------------
        double **xyz;
        int *TypeAtomId,*TypeResId;
        char **TypeAtom, **TypeRes, *Chain,**RealNum;
        int *AtomOnRes, *SosuiAtom, *ConservedAtom;
        int *ResOnAtom;
        int *ResNum,*AtomNum;
        double MaxXyz[3],MinXyz[3];
	float *DepthAtom,*DepthRes;
	float *SSP;
} PDB;


typedef struct{
	char filename[LIN],pfilename[LIN],sfilename[LIN],mfilename[LIN];
	double map_t;
	int Nthr;
	double dreso,LocalR;
	double MaxShift,MergeDist;
	double Filter,Dkeep;
	int Nround,Ntabu,Nnb,Nsim,Nbeam;
	int Mode;
	double Allow;
	float Pcut, Wp,Wd;
	int FragLen;
	int Cmode;//Color mode
	float Waa,Wss,Watm;
	float LowestSco;
} CMD;

/*
typedef struct{
	char filename[LIN];
	double map_t;
	int Nthreads;
	int xdim,ydim,zdim;
	int mx,my,mz;
	float xlen,ylen,zlen;
	float alpha,beta,gamma;
	int mapc,mapr,maps;
	int dmin,dmax,dmean,ispg;
	int nsymbt;
	float orgxyz[3];
	int NumVoxels;
	float *dens;
	float widthx,widthy,widthz;
	unsigned int Nact;
} MRC;
*/

//Mainmast Model
typedef struct{
        //char filename[LIN];
        float **xyz;
        double *b;//b-factor
        double *dis;//distance to the next coordinate.
        int NumOfCd;
        double sco,z,shake,zshake;//score
        int pos;//position
} MODEL;


typedef struct{
 int pos;//position
 int l;//length
 int *ss,*seq,cid;
 float Pss[10000][3];
 int CID[10000];
} SEQFG;


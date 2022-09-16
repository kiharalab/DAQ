#define VER 0.60
#define PI 3.141592
#define ATOM 50000
#define RES 5000
#define AMI 20
#define LIN 256
#define radian 57.29577951308232088

//grid
#define GRINO 20000

#define OUTSIDE -999.999

/*def �δĶ�����*/
#define INC 15/*ligand��ž���ٿ�*/
#define MOV 1/*ligand��ư���ٿ�*/
#define START_RT 0 /*ligand��ž��start*/
#define START_MV 0 /*ligand��ư��start*/
#define FIN_RT 360
#define FIN_MV 0
#define ON 1
#define OFF 0
#define TRUE 0
#define FALSE -1
/*#define MATFILE "result030106"*/ /*�������ޥȥꥯ���Υե�����*/
#define MATSIZE 1500
//#define MATLEN 600
#define MATAMI 24 /*�ޥȥꥯ���Υ��ߥλ�����*/

#define CSHNUM 20 /*���ͻĴ��������*/

/*�ޥ�������*/
#define X(a) (a)*(a)
#define L(a,b,c,d,e,f) (sqrt((d-a)*(d-a)+(e-b)*(e-b)+(f-c)*(f-c)))/*�٥��ȥ�Ĺ*/
#define RAS(a) (2.000000*PI*a/360.000000)/*��->�饸����*/
#define RAD(a) (2.000000*PI*a/360.000000)/*��->�饸����*/
#define ANG(a) (360.000000*a/(2.000000*PI))/*�饸����->��*/

#define GAUSS(a,b) 1.00/(sqrt(2.00*PI)*b)*exp(-(a*a)/(2.0*b*b))


/*triangle*/

#define MAXMTX   4

/*��̽���*/
#define TOP 3000

#define NOT_GAUSS 1
#define USE_GAUSS 0

#define VALIABLE 0
#define CONSTANT 1

int FINAL_ATNo;
int FINAL_RENo;
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
        int DSSPmapper[RES];
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
        double **xyz,**CAxyz;
        int *TypeAtomId,*TypeResId;
        char **TypeAtom, **TypeRes, *Chain,**RealNum;
        int *AtomOnRes, *SosuiAtom, *ConservedAtom;
        int *ResOnAtom;
        int *ResNum,*AtomNum;
        double MaxXyz[3],MinXyz[3];
	float *DepthAtom,*DepthRes;
} PDB;


typedef struct{
	char filename[LIN];
	char file1[LIN],file2[LIN];
	double map_t;
	int Nthr;
	double dreso,LocalR;
	double MaxShift,MergeDist;
	double Filter,Dkeep;
	int Nround,Ntabu,Nnb,Nsim,Nbeam;
	int Mode,Amode;
	double Allow;
	double th1,th2;
	double ssize,ang;
	double r;
	int TopN,Nvox;
	bool ShowGrid,IgnoreNC;
	bool Emode;
	int slide;
	double NoP1,NoP2,width;
	double Ubound,Lbound;
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



typedef struct{
 double **xyz;
 double *kappa;//angle i-1,i,i+1
 double *alpha;//torsion angle i-1,i,i+1,i+2
 double *phi,*psi;
 int *ss,*acc;
 int Nres;
}DSSP;

typedef struct{
 double d2;
 int AtomId,indx;
}VOXEL;





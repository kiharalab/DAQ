#define MAX_CD 10000
#define MAX_FLEN 200

//int BuildCa(double [5][MAX_FLEN][3],double [MAX_FLEN][3], int,int, MRC *,MEMO *);
int BuildCa(MEMO  *, int,int, MRC *);
double OptMainSide(double [MAX_FLEN][3],int *,int, MRC *);
double Ca_Optimizer(double **,double **, int,MEMO *);
void prepare_rbins_ca(double **,int **,int);

double rebuild_sidechain_fromCA(double **,double ***,int **,int *,int, MRC *,MEMO *, bool);
double rebuild_sidechain_fromCA2(double **,double ***,int **,int *,int, MRC *,MEMO *);
double rebuild_backbone_fromCA(double [MAX_FLEN][3],int *,int, MRC *, bool);
double rebuild_backbone_fromCA2(double **,double **,int **,int *,int);
double rebuild_backbone_fromCA3(double **,double **,int **,int *,int);

//double swap_helix(MEMO *,int *,int,int,MODEL *);
double swap_helix(MEMO *,SEQFG *,MODEL *);
//double opt_backbone_fromCA(MEMO *,int *,int *,int,int, MRC *, bool);
void opt_backbone_fromCA(MEMO *,SEQFG *, MRC *,MODEL *, bool);
double MapScore(double **,int,MRC *);
double MapScore2(double [][3],int,MRC *);
double CaEne(double [MAX_FLEN][3],int);
void VecMove(double **,double **,double **,int,double,struct drand48_data *);//drand48

int ClashChk(double **,int,double **,int,int,double);

void ShowFragCA(SEQFG *, double **);
void ShowFragBack(SEQFG *, double **);
void ShowFragSide(SEQFG *,double **, double ***);
void ShowFragBackSide(SEQFG *,double **, double ***);

double ShakeZsco(MODEL *,SEQFG *,double **, double ***,MRC *,struct drand48_data *, double ,int);

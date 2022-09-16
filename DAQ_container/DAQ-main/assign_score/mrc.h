
typedef struct{
	char filename[LIN];
	double map_t;
	int Nthreads;
	int xdim,ydim,zdim;
	int ncstart,nrstart,nsstart;
	int mx,my,mz;
	float xlen,ylen,zlen;
	float alpha,beta,gamma;
	int mapc,mapr,maps;
	float dmin,dmax,dmean;
	int ispg;
	int nsymbt;
	float orgxyz[3];
	int NumVoxels;
	float *dens;
	float widthx,widthy,widthz;
	unsigned int Nact;
} MRC;


typedef struct{
 double **cd,*dens;
 int **origrid;
 int Ncd,Nori;
 int *member;
 float *mask;
} POINTS;

typedef struct{
 double d;
 double dens;
 int id1,id2,eid;
 bool mst,local,keep;
} EDGE;

typedef struct{
 int N;
 EDGE *e[20];//!!!!! 10 is not enough
 int cid,id;//chain
 double AAP[20];//AA Probability
 double SSP[3];//SS Probability
 double ATOMP[6];//ATOM Probability
 float real_cd[3];//xyz coordinates in real space (not map grid space)
 float LogAA[21],LogSS[3],LogATOM[6];
} NODE;

typedef struct{
 int cut_id,add_id;
 double score;
 //float score;
} MOVE;



typedef struct{
 double len,bf_len;
 int Nnode,Ne;
 NODE *node;
 int Ntotal,Etotal;
 int St,Ed;
 bool *ActE;//Active edge
 bool *ActN;//Active node
 bool *MaskN;//Masking
 bool *UsedN;
 //bool *MstE://MST edge
 int *stock,Nstock;
 int *nextv;
 double *cost;
 int *cid;
 int *CutTbl,Ncut;//Edge table for cut;
 int *AddTbl,Nadd;//Edge table for add;

 MOVE *mv; //Movement
 int Nmv;
 int *Path,Lpath;
 double score;
} TREE;

typedef struct{
 EDGE *edge;//All edge data
 int Ne;
 bool **adj;
 NODE *node;//All connections
 int Nnode;
 int *cid;

 int Nt;
 EDGE **tree;

 //TREE mst;

} GRAPH;

typedef struct{
 EDGE *e_path[1000];//edge path
 int n_path[1000];//node path
 NODE *n_path_po[1000];//node pointer path
 int Nnode;//length
 double dist,dens,score;
 float dtbl[1000];//distance from 0
} FRAGMENT;

typedef struct{
 FRAGMENT *frag;
 int N;
} FRAG_DB;


typedef struct{
 double **prob;
 int N;
 //Statics data
 double REFaa[20],REFss[3],REFatm[5];
} PTBL;


//Sequence Node alignments
typedef struct{
 int pos;//position
 int len;//length
 float *score;
 int Nali;
 SEQFG *sfg;
 NODE ***ali;//alignment
} SEQ_NODE;



bool readmrc(MRC *,char *);
bool ToCubic(MRC *);
bool upsampling(MRC *,double);
void out_situs(MRC *);
bool meanshift(MRC *,POINTS *);
bool fastLDP(MRC *,POINTS *);
double meanshift_pos(MRC *,double [3]);
bool MergePoints(MRC *,POINTS *);
bool MergeLDP(MRC *,POINTS *);
void ShowModel(MRC *,POINTS *);
void ShowLDP(MRC *,POINTS *);
void ShowOri(MRC *,POINTS *);
void ShowGraph(GRAPH *);
void ShowTree(GRAPH *,TREE *);
void ShowPath(MRC *,POINTS *,GRAPH *,TREE *,int n);
void ShowPath2(MRC *,POINTS *,GRAPH *,TREE *,int n, SEQFG *);
void ShowThreadModel(GRAPH *,SEQFG *);

bool SetUpGraph(POINTS *, GRAPH *,MRC *,TREE *);
//bool Tabu(GRAPH *,TREE *,TREE *);
bool Tabu(GRAPH *,TREE *,TREE *,SEQFG *);
bool InitTree(TREE *,bool);
bool ConstTree(TREE *,GRAPH *);
double QualityTree(GRAPH *, TREE *);
//double QualityTreeDP(GRAPH *, TREE *, SEQFG *,DP_MEMORY *);
bool CopyTree(TREE *,TREE *,bool);
bool SetCutTbl(GRAPH *,TREE *,int *,int,double);
bool MoveTree(GRAPH *, TREE *, int,int);

int CutTree(GRAPH *, TREE *,double,int *,int);
int AddEdge(GRAPH *, TREE *);
int SplitChain(GRAPH *, TREE *,int);

double RandDouble();
int RandInt(int);
void ShuffleTbl(int *, int);
void ShuffleMv(MOVE *, int);

double Kendall(int *,int,int *,int);

bool PairExhaust(GRAPH *,TREE *,TREE *);
int ListStEd(GRAPH *,TREE *,int *);

double OptPath2points(GRAPH *,TREE *,TREE *,int,int,int);
bool get_path(GRAPH *,TREE *, int, int);

bool FindFrag(POINTS*,GRAPH *,FRAG_DB *,double, int);
bool read_ptbl(PTBL *,char *);
bool AddProbToMrc(PTBL *,MRC *,float,float,float);
bool AssignProbToNode(PTBL *,POINTS *,GRAPH *,MRC *);
bool AssignProbToNodePDB(PTBL *,PDB *,GRAPH *,MRC *);
bool ShowPTBL(MRC *,PTBL *,float);
bool ScoreFrag(FRAG_DB *,int, SEQFG *, int, SEQ_NODE *);
bool CombinationSearch(SEQFG *, int, SEQ_NODE *);

bool ProbQA(PDB *,GRAPH *,int);

float QuickScoring(FRAGMENT *,SEQFG *,NODE **);

bool MaskByPdb(MRC *,PDB *,float );


static char RES_NAMES[21][4] =
  { "ALA" ,"VAL" ,"PHE" ,"PRO" ,"MET" ,"ILE" 
        ,"LEU" ,"ASP" ,"GLU" ,"LYS" ,"ARG" ,"SER" 
        ,"THR" ,"TYR" ,"HIS" ,"CYS" ,"ASN" ,"TRP" 
        ,"GLN" ,"GLY" ,"UNK"
 };
static char SHORT_AA_NAMES[22] = { "AVFPMILDEKRSTYHCNWQGX" };


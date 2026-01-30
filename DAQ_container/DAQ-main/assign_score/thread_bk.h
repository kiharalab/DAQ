
#include "mrc.h"

typedef struct{
 char id[100];
 int len;
 int  seq_code[10000];
 int  ss[10000];
 char seq_txt[10000];
}  SEQ;

typedef struct{
 int st,ed,mod;
 float dist;
 bool nr;
} FRAG;

typedef struct{
 double **cd;//xyz
 int seq_pos;
 int N;//length
 double score;
} SEQ_FRAG;


/*
//Sequence fragments
typedef struct{
 int pos;//position
 int l;//length
 int *ss,*seq,cid;
} SEQFG;
*/

bool readseq(SEQ *,char *);
int splitseq(SEQ *,SEQFG *,int);

int load_models(MODEL **,char **,int,int *);
int load_single_model(MODEL *,char *);

bool AssignFeatures(MODEL **,int,FRAG **, int, MRC *);
//bool ThreadFrag(MODEL **,SEQ *,FRAG **,int,MRC *,MODEL *,MODEL **);
bool ThreadFrag(MODEL **,SEQFG *,int ,FRAG **,int,MRC *,MODEL *,MODEL **);

void ShowMainPath(MODEL **, int);
bool NR_frag(FRAG **,int *,MODEL **,int, double);
bool Evaluation(MODEL **,SEQFG *,int,int,MODEL *,MRC *);

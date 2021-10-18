

typedef struct{
 int t[3];
 double r[3];
 double q[4];
 int code;
 double sco;
} TBL;

bool SearchMAPfft(MRC *,MRC *,double);
bool SearchMAPfftMT(MRC *,MRC *,int,bool);
bool SearchMAPfftMT_OVCC(MRC *,MRC *,int,int,bool); //Overlap or CCC
double GetScore(MRC *, MRC *, int [3]);
double GetScore2(MRC *, MRC *, int [3], double [5]);


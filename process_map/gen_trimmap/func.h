int chkcmdline(int, char **,CMD *);
int readpdb(PDB *, char *,int);//include malloc
int ReadDssp(DSSP *, char *);
int AA2int(char *);
int CountAtom(char *);
int MallocPdb(PDB *,int);
int BaseVec(double [3],double [3], double [3],double [3][3]);
double Random();

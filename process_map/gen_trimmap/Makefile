#OBJECTS = chkcmdline.o readMAT.o readpdb.o sch.o Gcounter.o TOZ.o KAITEN.o VEC_kei.o kaiten.o cd_func.o SF_counter.o sch_tri.o Mcounter.o MVandSCORE.o lsfit.o det_sf.o enephobic_ca.o high_reso_doc.o

CC = gcc
#COM =  readpdb.o  cd_func.o chkcmdline.o rmsd.o
COM =  chkcmdline.o mrc.o readpdb.o
OBJECTS = $(COM)
HEADER = struct.h func.h mrc.h
PG = TrimMapAtom

OPT = -O3 -msse4.2 -msse -msse2 -msse3 -mmmx -ftree-vectorize -ftree-vectorizer-verbose=3 -ffast-math -march=native -funroll-loops -fopenmp

$(PG) : main.c $(OBJECTS) $(HEADER) $(COM)
	$(CC) -o $(PG) main.c $(OBJECTS) -lm $(OPT) 


chkcmdline.o : chkcmdline.c $(HEADER)
	$(CC) -O3 -c chkcmdline.c

readpdb.o : readpdb.c $(HEADER)
	$(CC) -O3 -c readpdb.c -lm $(OPT)

mrc.o : mrc.c $(HEADER)
	$(CC) -O3 -c mrc.c -lm $(OPT)

mrcfft.o : mrcfft.c $(HEADER) mrcfft.h
	$(CC) -O3 -c mrcfft.c -lm $(OPT) 

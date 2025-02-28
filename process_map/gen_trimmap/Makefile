
CC = gcc
COM =  chkcmdline.o mrc.o readpdb.o
OBJECTS = $(COM)
HEADER = struct.h func.h mrc.h
PG = TrimMapAtom

# OPT = -O3 -msse4.2 -msse -msse2 -msse3 -mmmx -ftree-vectorize -ftree-vectorizer-verbose=3 -ffast-math -funroll-loops -fopenmp
OPT = -O3 -ffast-math -ftree-vectorize -funroll-loops -fopenmp 

$(PG) : main.c $(OBJECTS) $(HEADER) $(COM)
	$(CC) -o $(PG) main.c $(OBJECTS) -lm $(OPT) 

chkcmdline.o : chkcmdline.c $(HEADER)
	$(CC) -c chkcmdline.c

readpdb.o : readpdb.c $(HEADER)
	$(CC) -c readpdb.c -lm $(OPT)

mrc.o : mrc.c $(HEADER)
	$(CC) -c mrc.c -lm $(OPT)

mrcfft.o : mrcfft.c $(HEADER) mrcfft.h
	$(CC) -c mrcfft.c -lm $(OPT) 

clean:
	rm -f $(PG) *.o


CC = gcc
COM =  chkcmdline.o mrc.o readpdb.o thread.o dp.o
OBJECTS = $(COM)
HEADER = struct.h func.h mrc.h thread.h
PG = DAQscore_colab

# OPT = -O3 -msse4.2 -msse -msse2 -msse3 -mmmx -ftree-vectorize -ftree-vectorizer-verbose=3 -ffast-math -funroll-loops -fopenmp
OPT = -O3 -ftree-vectorize -ffast-math -funroll-loops -fopenmp

$(PG) : main.c $(OBJECTS) $(HEADER) $(COM)
	$(CC) -o $(PG) main.c $(OBJECTS) -lm $(OPT) 

chkcmdline.o : chkcmdline.c $(HEADER)
	$(CC) -c chkcmdline.c

readpdb.o : readpdb.c $(HEADER)
	$(CC) -c readpdb.c -lm $(OPT)

mrc.o : mrc.c $(HEADER)
	$(CC) -c mrc.c -lm $(OPT)

thread.o : thread.c $(HEADER)
	$(CC) -c thread.c -lm $(OPT)

dp.o : dp.c $(HEADER)
	$(CC) -c dp.c -lm $(OPT)

clean:
	rm -f $(PG) *.o

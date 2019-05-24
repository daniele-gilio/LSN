CC = g++
CFLAGS = -Wall -O3 --std=c++11

tsp.exe : tsp.o random.o
	$(CC) random.o tsp.o -o tsp.exe
tsp.o : tsp.cpp tsp.h
	$(CC) -c tsp.cpp -o tsp.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o seed.out tsp.exe best.dat

all: tsp.exe

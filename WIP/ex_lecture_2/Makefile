CC = g++
CFLAGS = -Wall -O3 --std=c++11

es21.exe : es21.o random.o
	$(CC) random.o es21.o -o es21.exe
es21.o : es21.cpp
	$(CC) -c es21.cpp -o es21.o $(CFLAGS)

es22.exe : es22.o random.o
	$(CC) random.o es22.o -o es22.exe
es22.o : es22.cpp
	$(CC) -c es22.cpp -o es22.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o seed.out es21.exe es22.exe st_sample.dat imp_sample.dat disc.dat cont.dat sigma.dat

all: es21.exe es22.exe

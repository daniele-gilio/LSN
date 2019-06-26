CC = g++
CFLAGS = -Wall -O3 --std=c++11

es31.exe : es31.o random.o
	$(CC) random.o es31.o -o es31.exe
es31.o : es31.cpp
	$(CC) -c es31.cpp -o es31.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o seed.out es31.exe

all: es31.exe

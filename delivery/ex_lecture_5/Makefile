CC = g++
CFLAGS = -Wall -O3 --std=c++11

es5.exe : es5.o random.o
	$(CC) random.o es5.o -o es5.exe
es5.o : es5.cpp
	$(CC) -c es5.cpp -o es5.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o seed.out es5.exe

all: es5.exe

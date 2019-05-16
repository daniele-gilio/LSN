CC = g++
CFLAGS = -Wall -O3 --std=c++11

vmc.exe : vmc.o random.o
	$(CC) random.o vmc.o -o vmc.exe
vmc.o : vmc.cpp
	$(CC) -c vmc.cpp -o vmc.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o seed.out vmc.exe

all: vmc.exe

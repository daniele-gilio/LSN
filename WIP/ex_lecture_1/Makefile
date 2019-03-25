CC = g++
CFLAGS = -Wall -O3 --std=c++11

#main.exe : main.o random.o
#	$(CC) random.o main.o -o main.exe
#main.o : main.cpp
#	$(CC) -c main.cpp -o main.o $(CFLAGS)

es11.exe : es11.o random.o
	$(CC) random.o es11.o -o es11.exe
es11.o : es11.cpp
	$(CC) -c es11.cpp -o es11.o $(CFLAGS)

es12.exe : es12.o random.o
	$(CC) random.o es12.o -o es12.exe
es12.o : es12.cpp
	$(CC) -c es12.cpp -o es12.o $(CFLAGS)
	
es13.exe : es13.o random.o
	$(CC) random.o es13.o -o es13.exe
es13.o : es13.cpp
	$(CC) -c es13.cpp -o es13.o $(CFLAGS)

angle_dist.exe : angle_dist.o random.o
	$(CC) random.o angle_dist.o -o angle_dist.exe
angle_dist.o : angle_dist.cpp
	$(CC) -c angle_dist.cpp -o angle_dist.o $(CFLAGS)

random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean :
	rm *.o seed.out es11.exe es12.exe es13.exe angle_dist.exe

all: es11.exe es12.exe es13.exe angle_dist.exe

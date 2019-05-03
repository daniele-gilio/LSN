CC = g++
CFLAGS = -Wall -O3 --std=c++11

Monte_Carlo_ISING_1D.exe : Monte_Carlo_ISING_1D.o random.o
	$(CC) random.o Monte_Carlo_ISING_1D.o -o Monte_Carlo_ISING_1D.exe
Monte_Carlo_ISING_1D.o : Monte_Carlo_ISING_1D.cpp
	$(CC) -c Monte_Carlo_ISING_1D.cpp -o Monte_Carlo_ISING_1D.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)

create_run.exe : create_run.o
		$(CC) create_run.o -o create_run.exe
create_run.o : create_run.cpp
		$(CC) -c create_run.cpp -o create_run.o $(CFLAGS)
clean :
	rm *.o Monte_Carlo_ISING_1D.exe create_run.exe
fclean:
	./clean.sh
all: Monte_Carlo_ISING_1D.exe create_run.exe

run:
	./create_run.exe
	sh run.sh

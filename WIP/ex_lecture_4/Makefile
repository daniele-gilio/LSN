CC = g++
CFLAGS = -Wall -O3 --std=c++11

MolDyn_NVE.exe : MolDyn_NVE.o random.o
	$(CC)  random.o MolDyn_NVE.o -o MolDyn_NVE.exe
MolDyn_NVE.o : MolDyn_NVE.cpp MolDyn_NVE.h random.h
	$(CC) -c MolDyn_NVE.cpp -o MolDyn_NVE.o $(CFLAGS)
random.o : random.cpp random.h
	$(CC) -c random.cpp -o random.o $(CFLAGS)
clean:
	./clean.sh
clean_solid:
	./clean_solid.sh
clean_liquid:
	./clean_liquid.sh
clean_gas:
	./clean_gas.sh
run :
	./run.sh
all: MolDyn_NVE.exe

# Do not use compiler flag -O2, -O3

all:
	gcc -o pBFS pBFS.c

parallel:
	gcc -fopenmp -o pmBFS pmBFS.c
	gcc -fopenmp -o pBFS pBFS.c

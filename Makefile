FLAGS=--std=c11 -O0 -g

parse.o: parse.c parse.h
	gcc $(FLAGS) -o parse.o parse.c

randsat: randsat.c
	gcc $(FLAGS) -o randsat randsat.c

walk: walk.c parse.o parse.h
	gcc $(FLAGS) -o walk walk.c parse.o

all: randsat walk

clean:
	rm -f randsat walk parse.o

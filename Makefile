FLAGS=--std=c11 -O0 -g

all: randsat walk

clean:
	rm -f randsat walk parse.o

parse.o: parse.c parse.h
	gcc $(FLAGS) -c parse.c

randsat: randsat.c
	gcc $(FLAGS) -o randsat randsat.c

walk: walk.c parse.o parse.h
	gcc $(FLAGS) -o walk walk.c parse.o

.PHONY: all clean

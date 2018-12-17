FLAGS=--std=c11 -O0 -g

all: randsat walk

clean:
	rm -f randsat walk inst.o

inst.o: inst.c inst.h
	gcc $(FLAGS) -c inst.c

graph.o: graph.c graph.h
	gcc $(FLAGS) -c graph.c

randsat: randsat.c
	gcc $(FLAGS) -o randsat randsat.c

walk: walk.c inst.o inst.h graph.o graph.h
	gcc $(FLAGS) -o walk walk.c inst.o graph.o

.PHONY: all clean

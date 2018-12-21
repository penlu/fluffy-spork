#FLAGS=--std=gnu11 -O0 -g
FLAGS=--std=gnu11 -O3 -g -pg

LD_FLAGS=-lm

all: randsat walk warn survey

clean:
	rm -f randsat walk warn survey inst.o graph.o util.o

inst.o: inst.c inst.h
	gcc $(FLAGS) -c inst.c $(LD_FLAGS)

graph.o: graph.c graph.h
	gcc $(FLAGS) -c graph.c $(LD_FLAGS)

util.o: util.c util.h
	gcc $(FLAGS) -c util.c $(LD_FLAGS)

randsat: randsat.c util.o util.h
	gcc $(FLAGS) -o randsat randsat.c util.o $(LD_FLAGS)

walk: walk.c inst.o inst.h graph.o graph.h util.o util.h
	gcc $(FLAGS) -o walk walk.c inst.o graph.o util.o $(LD_FLAGS)

warn: warn.c inst.o inst.h graph.o graph.h util.o util.h
	gcc $(FLAGS) -o warn warn.c inst.o graph.o util.o $(LD_FLAGS)

survey: survey.c inst.o inst.h graph.o graph.h util.o util.h
	gcc $(FLAGS) -o survey survey.c inst.o graph.o util.o $(LD_FLAGS)

.PHONY: all clean

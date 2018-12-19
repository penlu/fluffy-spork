FLAGS=--std=gnu11 -O0 -g
#FLAGS=--std=gnu11 -O3

all: randsat walk warn survey

clean:
	rm -f randsat walk warn survey inst.o graph.o util.o

inst.o: inst.c inst.h
	gcc $(FLAGS) -c inst.c

graph.o: graph.c graph.h
	gcc $(FLAGS) -c graph.c

util.o: util.c util.h
	gcc $(FLAGS) -c util.c

randsat: randsat.c util.o util.h
	gcc $(FLAGS) -o randsat randsat.c util.o

walk: walk.c inst.o inst.h graph.o graph.h util.o util.h
	gcc $(FLAGS) -o walk walk.c inst.o graph.o util.o

warn: warn.c inst.o inst.h graph.o graph.h util.o util.h
	gcc $(FLAGS) -o warn warn.c inst.o graph.o util.o

survey: survey.c inst.o inst.h graph.o graph.h util.o util.h
	gcc $(FLAGS) -o survey survey.c inst.o graph.o util.o

.PHONY: all clean

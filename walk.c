// Simple sat solve via walksat.

// Invoke with:
// ./walk steps prob*100
// e.g.:
// ./walk 10000 30

// At every iteration, walksat picks a random unsat clause, then
// with some probability it flips a random var in that clause; otherwise it
// flips the var that will minimize the number of clauses that are still unsat.

// provide the instance on stdin
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#include "graph.h"
#include "inst.h"
#include "util.h"

//#define MAX_STEPS (1 << 24)
//#define P_PARAM 30
#define PRINT_FREQ 100

int main(int argc, char *argv[]) {
  if (argc != 3) {
    printf("./walk steps prob*100\n");
    exit(0);
  }

  // read command line parameters
  int max_steps = strtol(argv[1], NULL, 0);
  int p_param = strtol(argv[2], NULL, 0);

  struct timespec ts;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  srand((unsigned int) ts.tv_nsec);

  // parse input instance
  struct inst inst;
  inst_parse(&inst);
  int N = inst.N;
  int M = inst.M;

  // generate factor graph
  struct graph graph;
  graph_make(&graph, &inst);

  // dump for debug
  //inst_show(&inst);
  //graph_show(&graph);

  int *v; // received var assigns
  int *c; // received clause unsat counts
  int steps = walk(max_steps, p_param, PRINT_FREQ, &graph, &v, &c);

  if (steps == max_steps) {
    printf("%d steps\n", steps);
    printf("unknown\n");
  }

  // double check sat of final assignment
  int unsat_c = inst_check(&inst, v);
  if (unsat_c) {
    printf("UNSAT! clause %d\n", unsat_c);
  }

  // clean up
  free(v);
  free(c);
  inst_free(&inst);
  graph_free(&graph);

  return 0;
}

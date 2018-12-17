// simple sat solve via walksat
// invoke with: ./walksat steps prob*100
// e.g.: ./walksat 10000000 30
// at every iteration, walksat picks a random unsat clause
// then with some probability it flips a random var in that clause
// otherwise it flips the var that will minimize the number of clauses
// still unsat

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
#define PRINT_FREQ 1

int main(int argc, char *argv[]) {
  if (argc != 3) {
    printf("./walksat steps prob*100\n");
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

  // maintain variable settings
  int *v = calloc(N, sizeof(int));

  // maintain clause sat counts
  int *c = calloc(M, sizeof(int));
  for (int a = 0; a < M; a++) {
    for (int b = 0; b < graph.f[a].k; b++) {
      if (graph.f[a].v[b].j == 1) {
        c[a]++; // all variables start 0
      }
    }
  }

  // initiate walksat
  int steps;
  for (steps = 0; steps < max_steps; steps++) {
    // gather unsat clauses
    // TODO incrementalize this
    int u = 0;
    int *unsat = NULL;
    for (int a = 0; a < M; a++) {
      if (c[a] == 0) {
        u++;
        unsat = realloc(unsat, u * sizeof(int));
        unsat[u - 1] = a;
      }
    }

    // check sat
    if (u == 0) {
      printf("%d steps\n", steps);
      printf("sat\n");
      for (int i = 0; i < N; i++) {
        printf("%d ", v[i]);
      }
      printf("\n");
      break;
    }

    // status printout
    if (steps % PRINT_FREQ == 0) {
      printf("%d steps: %d unsat\n", steps, u);
    }

    // select random unsat clause
    int uc = unsat[urand(u)];

    // pick a var to flip in clause uc
    int flip;
    if (urand(100) < p_param) {
      // flip a random var in the clause
      flip = urand(graph.f[uc].k);

      assert(flip < graph.f[uc].k);
    } else {
      // flip a var that minimizes the number of clauses unsat

      // count how many clauses each flip would leave unsat
      struct node_f ff = graph.f[uc];
      int *funsat = calloc(ff.k, sizeof(int));
      for (int b = 0; b < ff.k; b++) {
        struct node_v *vf = ff.v[b].v;
        int newv = !v[vf->i]; // what we'd flip this var to
        for (int m = 0; m < vf->k; m++) {
          int j = vf->f[m].j; // coefficient
          int a = vf->f[m].f->a;

          // increment funsat if this flip would kill the last sat var
          if (c[a] == 1 && j == -1 && newv == 0) {
            funsat[b]++;
          } else if (c[a] == 1 && j == 1 && newv == 1) {
            funsat[b]++;
          }

          // shouldn't be able to kill a sat var in someone with no sat vars
          assert(c[a] != 0 || j != -1 || newv != 0);
          assert(c[a] != 0 || j != 1 || newv != 1);
        }
      }

      // pick randomly amongst flips minimizing the number unsat
      // count ties
      int min = M + 1;
      int ties = 0;
      for (int b = 0; b < ff.k; b++) {
        if (funsat[b] < min) {
          min = funsat[b];
          ties = 0;
        } else if (funsat[b] == min) {
          ties++;
        }
      }

      // random pick with tiebreaker
      int orig = urand(ties + 1);
      int tiebreak = orig;
      for (int b = 0; b < ff.k; b++) {
        if (funsat[b] == min && tiebreak == 0) {
          flip = b;
          break;
        } else if (funsat[b] == min) {
          tiebreak--;
        }
      }

      assert(flip < graph.f[uc].k);
    }

    // flip the var and propagate new sat counts
    struct node_v *vf = graph.f[uc].v[flip].v;
    int newv = !v[vf->i];
    v[vf->i] = newv;
    for (int m = 0; m < vf->k; m++) {
      int j = vf->f[m].j; // coefficient
      int a = vf->f[m].f->a;

      c[a] += j*(1 - 2*newv); // modify sat count correctly
    }
  }

  if (steps == max_steps) {
    printf("unsat after %d steps\n", steps);
  }

  // clean up
  inst_free(&inst);
  graph_free(&graph);

  return 0;
}

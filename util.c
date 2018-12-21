#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "graph.h"
#include "util.h"

// generate a uniform random integer [0, max)
int urand(int max) {
  while (1) {
    int r = rand();
    // throw away values that don't go into the interval evenly
    if (r < (((unsigned long long) RAND_MAX + 1) / max) * max) {
      int ret = r % max;
      assert(0 <= ret && ret < max);
      return ret;
    }
  }
}

double absd(double x) {
  return x > 0 ? x : -x;
}

// max steps to take
// probability to flip arbitrarily (*100)
// print freq in steps
// graph to run over
// return pointer for var assignments
// return pointer for clause unsat counts
int walk(int max_steps, int p_param, int print_freq, struct graph *graph, int **vr, int **cr) {
  int N = graph->N;
  int M = graph->M;

  // maintain variable settings
  int *v = calloc(N, sizeof(int));
  if (vr) {
    *vr = v;
  }

  // maintain clause sat counts
  int *c = calloc(M, sizeof(int));
  if (cr) {
    *cr = c;
  }

  // maintain clause sat counts
  for (int a = 0; a < M; a++) {
    if (graph->f[a].k == 0) {
      printf("walk: clause %d empty\n", a);
      printf("walk: unsat\n");

      return -1;
    }
    for (int j = 0; j < graph->f[a].k; j++) {
      int i = graph->f[a].v[j].v->i;
      if (graph->f[a].v[j].j == 1 - v[i] * 2) {
        c[a]++;
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
      printf("walk: %d steps\n", steps);
      printf("walk: sat\n");
      for (int i = 0; i < N; i++) {
        printf("%d ", v[i]);
      }
      printf("\n");
      break;
    }

    // status printout
    printf("walk: status %d steps: %d unsat\n", steps, u);

    // select random unsat clause
    int uc = unsat[urand(u)];
    free(unsat);

    // pick a var to flip in clause uc
    int flip;
    if (urand(100) < p_param) {
      // flip a random var in the clause
      flip = urand(graph->f[uc].k);
    } else {
      // flip a var that minimizes the number of clauses unsat

      // count how many clauses each flip would leave unsat
      struct node_f ff = graph->f[uc];
      int *funsat = calloc(ff.k, sizeof(int));

      for (int b = 0; b < ff.k; b++) {
        struct node_v *vf = ff.v[b].v;
        int newv = !v[vf->i]; // what we'd flip this var to

        // hypothetical sat counts for attached clauses
        int *cf = calloc(M, sizeof(int));
        for (int m = 0; m < vf->k; m++) {
          int a = vf->f[m].f->a;
          cf[a] = c[a];
        }
        // count up flip effects
        for (int m = 0; m < vf->k; m++) {
          struct node_f *vff = vf->f[m].f;
          int a = vff->a;

          cf[a] += vf->f[m].j*(1 - 2*newv);
          assert(cf[a] >= 0);
        }
        // increment funsat if this flip would kill the last sat var in a clause
        for (int m = 0; m < vf->k; m++) {
          if (cf[vf->f[m].f->a] == 0) {
            funsat[b]++;
          }
          assert(funsat[b] <= M);
        }

        free(cf);
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

      free(funsat);

      assert(flip < graph->f[uc].k);
    }

    // flip the var and propagate new sat counts
    struct node_v *vf = graph->f[uc].v[flip].v;
    int newv = !v[vf->i];
    v[vf->i] = newv;
    for (int m = 0; m < vf->k; m++) {
      int j = vf->f[m].j; // coefficient
      int a = vf->f[m].f->a;

      c[a] += j*(1 - 2*newv); // modify sat count correctly
    }
  }

  // clean up
  if (!vr) {
    free(v);
  }
  if (!cr) {
    free(c);
  }

  return steps;
}

// slightly better, slightly buggy walksat with good perf
int walk2(int max_steps, int p_param, int print_freq, struct graph *graph, int **vr, int **cr) {
  int N = graph->N;
  int M = graph->M;

  // maintain variable settings
  int *v = calloc(N, sizeof(int));
  if (vr) {
    *vr = v;
  }

  // maintain clause sat counts
  int *c = calloc(M, sizeof(int));
  if (cr) {
    *cr = c;
  }

  // maintain clause sat counts
  for (int a = 0; a < M; a++) {
    if (graph->f[a].k == 0) {
      printf("walk: clause %d empty\n", a);
      printf("walk: unsat\n");

      return -1;
    }
    for (int b = 0; b < graph->f[a].k; b++) {
      int i = graph->f[a].v[b].v->i;
      if (graph->f[a].v[b].j == 1 - v[i] * 2) {
        c[a]++;
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
      printf("walk: %d steps\n", steps);
      printf("walk: sat\n");
      for (int i = 0; i < N; i++) {
        printf("%d ", v[i]);
      }
      printf("\n");
      break;
    }

    // status printout
    if (print_freq && steps % print_freq == 0) {
      printf("walk: %d steps: %d unsat\n", steps, u);
    }

    // select random unsat clause
    int uc = unsat[urand(u)];
    free(unsat);

    // pick a var to flip in clause uc
    int flip;
    if (urand(100) < p_param) {
      // flip a random var in the clause
      flip = urand(graph->f[uc].k);

      assert(flip < graph->f[uc].k);
    } else {
      // flip a var that minimizes the number of clauses unsat

      // count how many clauses each flip would leave unsat
      struct node_f ff = graph->f[uc];
      int *funsat = calloc(ff.k, sizeof(int));
      for (int b = 0; b < ff.k; b++) {
        struct node_v *vf = ff.v[b].v;
        int newv = !v[vf->i]; // what we'd flip this var to
        for (int m = 0; m < vf->k; m++) {
          int j = vf->f[m].j; // coefficient
          int a = vf->f[m].f->a;

          // increment funsat if this flip would kill the last sat var
          if (c[a] == 1 && j == newv*2 - 1) {
            funsat[b]++;
          }

          // shouldn't be able to kill a sat var in someone with no sat vars
          assert(c[a] != 0 || j != newv*2 - 1);
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
      assert(min != M + 1);

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

      free(funsat);

      assert(flip < graph->f[uc].k);
    }
    assert(flip < graph->f[uc].k);

    // flip the var and propagate new sat counts
    struct node_v *vf = graph->f[uc].v[flip].v;
    int newv = !v[vf->i];
    v[vf->i] = newv;
    for (int m = 0; m < vf->k; m++) {
      int j = vf->f[m].j; // coefficient
      int a = vf->f[m].f->a;

      c[a] += j*(1 - 2*newv); // modify sat count correctly
    }
  }

  // clean up
  if (!vr) {
    free(v);
  }
  if (!cr) {
    free(c);
  }

  return steps;
}

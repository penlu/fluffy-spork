// Simple sat solve using warning-inspired decimation.

// Invoke with: ./warn iters
// e.g.: ./warn 10000

// provide the instance on stdin
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <time.h>

#include "graph.h"
#include "inst.h"
#include "util.h"

#define PRINT_FREQ 100

int main(int argc, char *argv[]) {
  if (argc != 2) {
    printf("./warn iters\n");
    exit(0);
  }

  // read command line parameters
  int max_iters = strtol(argv[1], NULL, 0);

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

  // track original indices to know the final sat assignment
  int *orig_vi = calloc(graph.N, sizeof(int));
  int *orig_v = calloc(graph.N, sizeof(int));
  for (int i = 0; i < graph.N; i++) {
    orig_vi[i] = i;
    orig_v[i] = 0;
  }

  // XXX parallel update; is the random order important?
  int steps = 0;
  while (graph.N != 0) {
    assert(graph.N > 0);
    assert(graph.M >= 0);

    // create original warning storage
    // u[i][a] is warning message a->i
    int **pu = calloc(graph.N, sizeof(int*));
    for (int i = 0; i < graph.N; i++) {
      pu[i] = calloc(graph.v[i].k, sizeof(int));
    }

    // randomly initialize warnings
    for (int i = 0; i < graph.N; i++) {
      for (int a = 0; a < graph.v[i].k; a++) {
        pu[i][a] = urand(2);
      }
    }

    // update warnings for max_iters timesteps
    int iters;
    int converged = 0;
    for (iters = 0; iters < max_iters && !converged; iters++) {
      converged = 1;

      // create cavity field storage
      // h[a][i] is cavity field i->a
      int **h = calloc(graph.M, sizeof(int*));
      for (int a = 0; a < graph.M; a++) {
        h[a] = calloc(graph.f[a].k, sizeof(int));
      }

      // compute cavity fields h
      for (int a = 0; a < graph.M; a++) {
        for (int i = 0; i < graph.f[a].k; i++) {
          // calculate new h[a][i], h_{i -> a}
          struct node_v *vi = graph.f[a].v[i].v;
          h[a][i] = 0;
          for (int b = 0; b < vi->k; b++) {
            if (a != vi->f[b].f->a) { // u_{b -> i} for a != b
              h[a][i] -= vi->f[b].j * pu[vi->i][b];
            }
          }
        }
      }

      // create new warning storage
      // u[i][a] is warning message a->i
      int **u = calloc(graph.N, sizeof(int*));
      for (int i = 0; i < graph.N; i++) {
        u[i] = calloc(graph.v[i].k, sizeof(int));
      }

      // compute new warnings u
      for (int i = 0; i < graph.N; i++) {
        for (int a = 0; a < graph.v[i].k; a++) {
          // calculate new u[i][a], u_{a -> i}
          struct node_f *fa = graph.v[i].f[a].f;
          u[i][a] = 1;
          for (int j = 0; j < fa->k; j++) {
            if (i != fa->v[j].v->i) { // h_{j -> a} for i != j
              u[i][a] *= (h[fa->a][j] * fa->v[j].j > 0) ? 1 : 0;
            }
          }
          if (u[i][a] != pu[i][a]) {
            converged = 0;
          }
        }
      }

      // clean up cavity fields
      for (int a = 0; a < graph.M; a++) {
        free(h[a]);
      }
      free(h);

      // clean up old warnings
      for (int i = 0; i < graph.N; i++) {
        free(pu[i]);
      }
      free(pu);

      pu = u;
    }

    // check for convergence
    if (iters == max_iters) {
      printf("%d steps %d iters\n", steps, iters);
      printf("unconverged\n");

      // clean up warnings
      for (int i = 0; i < graph.N; i++) {
        free(pu[i]);
      }
      free(pu);

      break;
    }

    // compute local fields and contradiction numbers
    int *H = calloc(graph.N, sizeof(int));
    int *c = calloc(graph.N, sizeof(int));
    int unsat = 0;
    for (int i = 0; i < graph.N; i++) {
      int pc = 0;
      int nc = 0;

      for (int a = 0; a < graph.v[i].k; a++) {
        assert(pu[i][a] == 0 || pu[i][a] == 1);

        // local field
        H[i] -= graph.v[i].f[a].j * pu[i][a];

        // contradiction number
        if (graph.v[i].f[a].j == -1) {
          pc += pu[i][a];
        } else {
          nc += pu[i][a];
        }
      }

      if (pc * nc > 0) {
        c[i] = 1;
        unsat = 1;
      }
    }

    // clean up warnings
    for (int i = 0; i < graph.N; i++) {
      free(pu[i]);
    }
    free(pu);

    // clean up contradiction numbers
    free(c);

    if (unsat) {
      printf("%d steps %d iters\n", steps, iters);
      printf("unsat\n");

      free(H);

      break;
    }

    // fix variables according to warnings and simplify clauses

    // compute fixed variables
    // 0: don't fix, +1: fix 1, -1: fix 0
    int *v = calloc(graph.N, sizeof(int));
    int fv = 0; // fixed variable count
    int *nv = calloc(graph.N, sizeof(int)); // new variable index
    int vi = 0;
    for (int i = 0; i < graph.N; i++) {
      if (H[i] > 0) {
        fv++;
        v[i] = 1;
      } else if (H[i] < 0) {
        fv++;
        v[i] = -1;
      } else {
        nv[i] = vi++;
      }
    }
    assert(fv + vi == graph.N);

    // clean up local fields
    free(H);

    // all local fields zero: fix an arbitrary variable
    if (vi == graph.N) {
      // fix a random variable
      int rv = urand(graph.N);
      v[rv] = urand(2) * 2 - 1;

      // reset new indices
      free(nv);
      fv = 1;
      nv = calloc(graph.N, sizeof(int));
      vi = 0;
      for (int i = 0; i < graph.N; i++) {
        if (i != rv) {
          nv[i] = vi++;
        }
      }
    }
    assert(fv + vi == graph.N);

    // save fixed variables
    for (int i = 0; i < graph.N; i++) {
      if (v[i] != 0) {
        assert(orig_v[orig_vi[i]] == 0);
        orig_v[orig_vi[i]] = v[i];
      }
    }

    // compute sat clauses
    // 1: satisfied, 0: still unsat
    int *f = calloc(graph.M, sizeof(int));
    int ff = 0; // fixed factors
    int *nf = calloc(graph.M, sizeof(int)); // new clause index
    int fa = 0;
    for (int a = 0; a < graph.M; a++) {
      for (int i = 0; i < graph.f[a].k; i++) {
        int j = graph.f[a].v[i].v->i;
        if (graph.f[a].v[i].j * v[j] == -1) {
          ff++;
          f[a] = 1;
          break;
        }
      }
      if (f[a] == 0) {
        nf[a] = fa++;
      }
    }
    assert(ff + fa == graph.M);

    // construct new decimated graph
    struct graph ngraph;
    ngraph.N = graph.N - fv;
    ngraph.M = graph.M - ff;
    assert(ngraph.N == vi);
    assert(ngraph.M == fa);
    ngraph.v = NULL;
    ngraph.f = NULL;

    if (steps % PRINT_FREQ == 0) {
      printf("%d vars unset\n", ngraph.N);
    }

    // create storage for unfixed variables
    ngraph.v = calloc(ngraph.N, sizeof(struct node_v));
    for (int i = 0; i < ngraph.N; i++) {
      ngraph.v[i].i = i;
    }

    // create storage for unfixed clauses
    ngraph.f = calloc(ngraph.M, sizeof(struct node_f));
    for (int a = 0; a < ngraph.M; a++) {
      ngraph.f[a].a = a;
    }

    // create new edges
    // implicitly skip edges for set vars:
    // if the set var sat a clause, then that clause is already gone
    // if the set var doesn't sat a clause, then we still exclude it
    int *new_vi = calloc(ngraph.N, sizeof(int));
    for (int i = 0; i < graph.N; i++) {
      if (v[i] == 0) { // variable still extant
        int ni = nv[i];

        // track original index
        new_vi[ni] = orig_vi[i];

        for (int a = 0; a < graph.v[i].k; a++) {
          int b = graph.v[i].f[a].f->a;
          if (f[b] == 0) { // factor not sat yet
            int na = nf[b];

            struct edge e;
            e.j = graph.v[i].f[a].j; // copy edge coefficient
            e.v = &ngraph.v[ni];
            e.f = &ngraph.f[na];

            // add edge to variable
            int vk = ++ngraph.v[ni].k;
            ngraph.v[ni].f = realloc(ngraph.v[ni].f, vk * sizeof(struct edge));
            ngraph.v[ni].f[vk - 1] = e;

            // add edge to clause
            int fk = ++ngraph.f[na].k;
            ngraph.f[na].v = realloc(ngraph.f[na].v, fk * sizeof(struct edge));
            ngraph.f[na].v[fk - 1] = e;
          }
        }
      }
    }
    free(orig_vi);
    orig_vi = new_vi;

    // clean up variable/clause-fixing scratch
    free(v);
    free(nv);
    free(f);
    free(nf);

    // store new graph
    graph_free(&graph);
    graph = ngraph;

    steps++;
  }

  if (graph.N == 0) {
    printf("%d steps\n", steps);
    printf("sat\n");

    // output sat assignment
    for (int i = 0; i < N; i++) {
      assert(orig_v[i] != 0);
      printf("%d ", (orig_v[i] + 1) / 2);
    }
    printf("\n");
  }

  // double-check sat of final assignment
  int *assigns = calloc(N, sizeof(int));
  for (int i = 0; i < N; i++) {
    assigns[i] = (orig_v[i] + 1) / 2;
  }
  int unsat_c = inst_check(&inst, assigns);
  if (unsat_c) {
    printf("UNSAT! clause %d\n", unsat_c);
  }
  assert(graph.N != 0 || unsat_c == 0);
  free(assigns);

  // clean up
  free(orig_vi);
  free(orig_v);
  inst_free(&inst);
  graph_free(&graph);

  return 0;
}

// Simple sat solve using survey-inspired decimation.

// Invoke with: ./survey iters steps p
// e.g.: ./survey 10000 1000 30

// provide the instance on stdin
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <time.h>
#include <omp.h>

#include "graph.h"
#include "inst.h"
#include "util.h"

#define PRINT_FREQ 100
#define TOLERANCE 0.001

int main(int argc, char *argv[]) {
  if (argc != 4) {
    printf("./survey iters steps p\n");
    exit(0);
  }

  // read command line parameters
  int max_iters = strtol(argv[1], NULL, 0);
  int max_steps = strtol(argv[2], NULL, 0); // for walking
  int p_param = strtol(argv[3], NULL, 0);   // for walking

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

  // debug dump
  //inst_show(&inst);

  // track original indices to know the final sat assignment
  int *orig_vi = calloc(graph.N, sizeof(int));
  int *orig_v = calloc(graph.N, sizeof(int));
  for (int i = 0; i < graph.N; i++) {
    orig_vi[i] = i;
    orig_v[i] = 0;
  }

  // create original survey storage
  // eta[i][a] is survey a->i
  float **p_eta = calloc(graph.N, sizeof(float*));
  for (int i = 0; i < graph.N; i++) {
    p_eta[i] = calloc(graph.v[i].k, sizeof(float));
  }

  // randomly initialize warnings
  int edges = 0;
  for (int i = 0; i < graph.N; i++) {
    for (int a = 0; a < graph.v[i].k; a++) {
      p_eta[i][a] = (float) urand(1 << 30) / ((1 << 30) - 1);
      assert(0 <= p_eta[i][a] && p_eta[i][a] <= 1);
      edges++;
    }
  }

  // XXX parallel update; is the random order important?
  int steps = 0;
  int tot_iters = 0;
  clock_gettime(CLOCK_MONOTONIC, &ts);
  unsigned int start = ts.tv_sec;

  int *conv_th = calloc(8, sizeof(int));
  int *zero_th = calloc(8, sizeof(int));

  omp_lock_t *flocks = calloc(graph.M, sizeof(omp_lock_t));
  #pragma omp parallel for num_threads(8)
  for (int a = 0; a < graph.M; a++) {
    omp_init_lock(&flocks[a]);
  }

  while (graph.N != 0) {
    //graph_check(&graph);
    assert(graph.N > 0);
    assert(graph.M >= 0);

    // create weight storage (pi values, effectively cavity fields)
    // pi*[a][i] is weight i->a
    float **pi_u = calloc(graph.M, sizeof(float*));
    float **pi_s = calloc(graph.M, sizeof(float*));
    float **pi_0 = calloc(graph.M, sizeof(float*));
    for (int a = 0; a < graph.M; a++) {
      pi_u[a] = calloc(graph.f[a].k, sizeof(float));
      pi_s[a] = calloc(graph.f[a].k, sizeof(float));
      pi_0[a] = calloc(graph.f[a].k, sizeof(float));
    }

    // update surveys for max_iters timesteps
    int iters;
    int converged = 0;
    int zeros = 0;
    for (iters = 0; iters < max_iters && !converged; iters++) {
      //printf("iter %d\n", iters);
      // compute fields pi*
      #pragma omp parallel for num_threads(8)
      for (int a = 0; a < graph.M; a++) {
        for (int i = 0; i < graph.f[a].k; i++) {
          // calculate new pi*[a][i], pi*_{i -> a}
          struct node_v *vi = graph.f[a].v[i].v;

          // calculate products separately
          float pu = 1;
          float ps = 1;
          float p0 = 1;
          for (int b = 0; b < vi->k; b++) {
            if (a != vi->f[b].f->a) { // u_{b -> i} for a != b
              if (vi->f[b].j != graph.f[a].v[i].j) {
                // different edge type
                pu *= 1 - p_eta[vi->i][b];
              } else {
                // same edge type
                ps *= 1 - p_eta[vi->i][b];
              }
              // all edges
              p0 *= 1 - p_eta[vi->i][b];
            }
          }

          pi_u[a][i] = ps - p0;
          pi_s[a][i] = pu - p0;
          pi_0[a][i] = p0;

          assert(0 <= pi_u[a][i] && pi_u[a][i] <= 1);
          assert(0 <= pi_s[a][i] && pi_s[a][i] <= 1);
          assert(0 <= pi_0[a][i] && pi_0[a][i] <= 1);
        }
      }

      // create new survey storage
      // eta[i][a] is survey a->i
      float **eta = calloc(graph.N, sizeof(float*));
      for (int i = 0; i < graph.N; i++) {
        eta[i] = calloc(graph.v[i].k, sizeof(float));
      }

      // compute new surveys eta
      #pragma omp parallel for num_threads(8)
      for (int i = 0; i < graph.N; i++) {
        int tid = omp_get_thread_num();
        for (int a = 0; a < graph.v[i].k; a++) {
          // calculate new eta[i][a], eta_{a -> i}
          struct node_f *fa = graph.v[i].f[a].f;
          eta[i][a] = 1;
          for (int j = 0; j < fa->k; j++) {
            if (i != fa->v[j].v->i) {
              float denom = pi_u[fa->a][j] + pi_s[fa->a][j] + pi_0[fa->a][j];
              if (pi_u[fa->a][j] == 0 && denom == 0) {
                eta[i][a] = 0;
              } else {
                eta[i][a] *= pi_u[fa->a][j] / denom;
              }
            }
          }

          assert(0 <= p_eta[i][a] && p_eta[i][a] <= 1);

          if (absd(eta[i][a] - p_eta[i][a]) > TOLERANCE) {
            conv_th[tid] = 0;
          }
          if (absd(eta[i][a]) < TOLERANCE) {
            zero_th[tid]++;
          }
        }
      }

      // accumulate convergence
      converged = 1;
      for (int i = 0; i < 8; i++) {
        converged = converged && conv_th[i];
        conv_th[i] = 1;
      }

      // accumulate zeros
      zeros = 0;
      for (int i = 0; i < 8; i++) {
        zeros += zero_th[i];
        zero_th[i] = 0;
      }

      // clean up old warnings
      for (int i = 0; i < graph.N; i++) {
        free(p_eta[i]);
      }
      free(p_eta);

      p_eta = eta;
    }

    tot_iters += iters - 1;

    // check for convergence
    if (!converged) {
      printf("survey: %d steps %d iters\n", steps, iters);
      printf("survey: unconverged after %d steps\n", steps);
      printf("survey: unconverged\n");

      break;
    }

    printf("survey: %d steps: converged in %d iters\n", steps, iters);

    // calculate clause complexities
    float sigma = 0;
    float *sigma_a = calloc(graph.M, sizeof(float));
    #pragma omp parallel for num_threads(8) reduction(+:sigma)
    for (int a = 0; a < graph.M; a++) {
      float prod_u = 1;
      float prod_denom = 1;
      for (int i = 0; i < graph.f[a].k; i++) {
        prod_u *= pi_u[a][i];
        prod_denom *= (pi_u[a][i] + pi_s[a][i] + pi_0[a][i]);
      }

      sigma += log(prod_denom - prod_u);
      sigma_a[a] = log(prod_denom - prod_u);
    }

    // clean up products
    for (int a = 0; a < graph.M; a++) {
      free(pi_u[a]);
      free(pi_s[a]);
      free(pi_0[a]);
    }
    free(pi_u);
    free(pi_s);
    free(pi_0);

    // compute biases
    float *Wp = calloc(graph.N, sizeof(float));
    float *Wm = calloc(graph.N, sizeof(float));
    float *W0 = calloc(graph.N, sizeof(float));
    float *sigma_i = calloc(graph.N, sizeof(float));
    #pragma omp parallel for num_threads(8) reduction(-:sigma)
    for (int i = 0; i < graph.N; i++) {
      float pp = 1;
      float pm = 1;
      float p0 = 1;

      for (int a = 0; a < graph.v[i].k; a++) {
        if (graph.v[i].f[a].j == -1) {
          pp *= 1 - p_eta[i][a];
        } else {
          pm *= 1 - p_eta[i][a];
        }
        // all edges
        p0 *= 1 - p_eta[i][a];
      }

      float pi_p = pm - p0;
      float pi_m = pp - p0;
      float pi_0 = p0;

      float denom = pi_p + pi_m + pi_0;
      assert(denom != 0);
      if (pi_p == 0) {
        Wp[i] = 0;
      } else {
        Wp[i] = pi_p / denom;
      }
      if (pi_m == 0) {
        Wm[i] = 0;
      } else {
        Wm[i] = pi_m / denom;
      }
      if (pi_0 == 0) {
        W0[i] = 0;
      } else {
        W0[i] = pi_0 / denom; // more numerically stable than 1 - Wp[i] - Wm[i]
      }

      // accumulate variable complexity
      sigma -= (graph.v[i].k - 1) * log(pi_p + pi_m + pi_0);
      sigma_i[i] = log(pi_p + pi_m + pi_0);
      assert(sigma_i[i] <= 0);

      assert(0 <= Wp[i] && Wp[i] <= 1);
      assert(0 <= Wm[i] && Wm[i] <= 1);
      assert(0 <= W0[i] && W0[i] <= 1);
    }

    printf("survey: %d steps sigma/N %.*e\n", steps, FLT_DIG, sigma / graph.N);

    // all surveys zero: start walking
    if (zeros == edges) {
      printf("survey: %d steps %d iters\n", steps, iters);
      printf("survey: starting walk at %d steps\n", steps);
      printf("survey: surveys trivial\n");

      int *v;
      int *c;
      int walk_steps = walk(max_steps, p_param, PRINT_FREQ*10, &graph, &v, &c);

      // save walk settings
      for (int i = 0; i < graph.N; i++) {
        orig_v[orig_vi[i]] = v[i] * 2 - 1;
      }

      // process sat/unsat here
      if (walk_steps == -1) {
        printf("survey: %d steps %d walk\n", steps, walk_steps);
        printf("survey: walk unsat\n");
        printf("survey: unsat\n");
      } else if (walk_steps < max_steps) {
        printf("survey: %d steps %d walk\n", steps, walk_steps);
        printf("survey: sat\n");

        // output sat assignment
        for (int i = 0; i < N; i++) {
          assert(orig_v[i] != 0);
          printf("%d ", (orig_v[i] + 1) / 2);
        }
        printf("\n");
      } else {
        printf("survey: %d steps %d walk\n", steps, walk_steps);
        printf("survey: walk timed out\n");
        printf("survey: unknown after %d steps\n", steps);
        printf("survey: unknown\n");
      }

      // clean up walk
      free(v);
      free(c);

      // clean up biases
      free(Wp);
      free(Wm);
      free(W0);

      // clean up complexities
      free(sigma_a);
      free(sigma_i);

      break;
    }

    // fix variables according to biases and simplify clauses

    // compute fixed variables
    // 0: don't fix, +1: fix 1, -1: fix 0
    float max = 0;
    int dir = 0;
    int fix = 0;
    float max_sigma = -INFINITY;
    for (int i = 0; i < graph.N; i++) {
      // use certitudes
      /*
      float d = Wp[i] < Wm[i] ? 1 - Wp[i] : 1 - Wm[i];
      if (d > max) {
        max = d;
        dir = Wp[i] > Wm[i] ? 1 : -1;
        fix = i;
      }
      */

      // use polarization
      float d = Wp[i] - Wm[i];
      if (absd(d) > max) {
        max = absd(d);
        dir = d > 0 ? 1 : -1;
        fix = i;
      }

      if (sigma_i[i] > max_sigma) {
        max_sigma = sigma_i[i];
      }
    }
    printf("survey: selected %d (est delta %.*e, pol %.*e, 1-Wm %.*e), var sigma %.*e (max %.*e)\n",
      fix,
      FLT_DIG, log(Wp[fix] < Wm[fix] ? 1 - Wp[fix] : 1 - Wm[fix]),
      FLT_DIG, absd(Wp[fix] - Wm[fix]),
      FLT_DIG, 1 - Wm[fix],
      FLT_DIG, sigma_i[fix],
      FLT_DIG, max_sigma);

    // save fixed variable
    assert(orig_v[orig_vi[fix]] == 0);
    orig_v[orig_vi[fix]] = dir;

    // clean up biases
    free(Wp);
    free(Wm);
    free(W0);

    // clean up complexities
    free(sigma_a);
    free(sigma_i);

    // compute sat clauses
    int *f = calloc(graph.M, sizeof(int)); // 1: satisfied, 0: still unsat
    int ff = 0; // count of satisfied factors
    int would_ff = 0; // would be sat if we picked the other direction
    #pragma omp parallel for num_threads(8) reduction(+:ff,would_ff)
    for (int a = 0; a < graph.v[fix].k; a++) {
      int b = graph.v[fix].f[a].f->a;
      if (f[b] != 1 && graph.v[fix].f[a].j * dir == -1) {
        ff += 1;
        f[b] = 1;
      }
      if (f[b] * dir == 1) {
        would_ff += 1;
      }
    }
    //printf("survey: set %d to %d sat %d (vs. %d)\n", fix, dir, ff, would_ff);
    printf("survey: fixed %d to %d (k=%d), sats %d (opp sat %d)\n", fix, dir, graph.v[fix].k, ff, would_ff);

    // find unfixed factors
    int *nf = calloc(graph.M, sizeof(int)); // new clause index
    int fa = 0;
    for (int a = 0; a < graph.M; a++) {
      if (f[a] == 0) {
        nf[a] = fa;
        fa++;
      }
    }
    assert(ff + fa == graph.M);

    // construct new decimated graph
    struct graph ngraph;
    ngraph.N = graph.N - 1;
    ngraph.M = graph.M - ff;
    ngraph.v = NULL;
    ngraph.f = NULL;

    // create storage for unfixed variables
    ngraph.v = calloc(ngraph.N, sizeof(struct node_v));
    #pragma omp parallel for num_threads(8)
    for (int i = 0; i < ngraph.N; i++) {
      ngraph.v[i].i = i;
    }

    // create storage for unfixed clauses
    ngraph.f = calloc(ngraph.M, sizeof(struct node_f));
    #pragma omp parallel for num_threads(8)
    for (int a = 0; a < ngraph.M; a++) {
      ngraph.f[a].a = a;
    }

    // also will need to decimate surveys
    // eta[i][a] is survey a->i
    float **eta = calloc(ngraph.N, sizeof(float*));

    // create new edges
    // implicitly skip edges for set vars:
    // if the set var sat a clause, then that clause is already gone
    // if the set var doesn't sat a clause, then we still exclude it
    int *new_vi = calloc(ngraph.N, sizeof(int));
    int ni = 0;
    edges = 0;
    for (int i = 0; i < graph.N; i++) {
      if (i != fix) { // variable still extant
        // track original index
        new_vi[ni] = orig_vi[i];

        eta[ni] = malloc(0);

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
            //omp_set_lock(flocks[na]);
            int fk = ++ngraph.f[na].k;
            ngraph.f[na].v = realloc(ngraph.f[na].v, fk * sizeof(struct edge));
            ngraph.f[na].v[fk - 1] = e;
            //omp_release_lock(flocks[na]);

            // copy old eta
            eta[ni] = realloc(eta[ni], vk * sizeof(float));
            eta[ni][vk - 1] = p_eta[i][a];
            edges++;
          }
        }
      }
      ni++;
    }

    // clean up old warnings
    for (int i = 0; i < graph.N; i++) {
      free(p_eta[i]);
    }
    free(p_eta);

    p_eta = eta;

    // detect empty clauses: conflict
    int unsat = 0;
    for (int a = 0; a < ngraph.M; a++) {
      if (ngraph.f[a].k == 0) {
        unsat = 1;
        break;
      }
    }
    if (unsat) {
      printf("survey: %d steps %d iters\n", steps, iters);
      printf("survey: unsat after %d steps\n", steps);
      printf("survey: unsat\n");

      free(new_vi);
      free(f);
      free(nf);
      graph_free(&ngraph);

      break;
    }

    printf("survey: status %d steps: %d vars, %d unsat\n", steps, ngraph.N, ngraph.M);

    free(orig_vi);
    orig_vi = new_vi;

    // clean up variable/clause-fixing scratch
    free(f);
    free(nf);

    // store new graph
    graph_free(&graph);
    graph = ngraph;

    steps++;
  }

  free(zero_th);
  free(conv_th);

  #pragma omp parallel for num_threads(8)
  for (int a = 0; a < graph.M; a++) {
    omp_destroy_lock(&flocks[a]);
  }

  printf("survey: %d iters total\n", tot_iters);

  // clean up warnings
  for (int i = 0; i < graph.N; i++) {
    free(p_eta[i]);
  }
  free(p_eta);

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

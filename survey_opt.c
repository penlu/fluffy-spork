// Simple sat solve using survey-inspired decimation.

// Invoke with: ./survey iters steps p
// e.g.: ./survey 10000 1000 30

// using custom data structures etc for nice code

// provide the instance on stdin
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include <time.h>

#include "graph.h"
#include "inst.h"
#include "util.h"

#define PRINT_FREQ 100
#define TOLERANCE 0.001

// factor graph
struct fg {
  int N;      // number of variables
  int M;      // number of factors

  int *v;     // variable assignments [i]
  int *f;     // factor sat counts [a]

  int *vk;    // number of factors per variable
  int **ve;   // variables' connected factors [i][b]
  int **vJ;   // edge weights

  int *fk;    // number of variables per factor
  int **fe;   // factors' connected variables [a][j]
  int **fJ;   // edge weights

  float **eta;  // surveys a->i [i][b]
  float **pi_u; // weight pi^u i->a [a][j]
  float **pi_s; // weight pi^s i->a [a][j]
  float **pi_0; // weight pi^0 i->a [a][j]

  float *Wp;    // polarization +
  float *Wm;    // polarization -
  float *W0;    // polarization 0

  float sigma;    // complexity (log number of clusters)
  float *sigma_v; // complexity contribution per variable
  float *sigma_f; // complexity contribution per clause
}

void fg_init(struct fg *graph, struct inst *inst) {
  // init counts
  int N = graph->N = inst->N;
  int M = graph->M = inst->M;

  // init storage
  graph->v = calloc(N, sizeof(int));
  graph->vk = calloc(N, sizeof(int));
  graph->ve = calloc(N, sizeof(int*));
  graph->vJ = calloc(N, sizeof(int*));

  graph->f = calloc(M, sizeof(int));
  graph->fk = calloc(M, sizeof(int));
  graph->fe = calloc(M, sizeof(int*));
  graph->fJ = calloc(M, sizeof(int*));

  // init messages
  graph->eta = calloc(N, sizeof(float*));
  graph->pi_u = calloc(M, sizeof(float*));
  graph->pi_s = calloc(M, sizeof(float*));
  graph->pi_0 = calloc(M, sizeof(float*));

  graph->Wp = calloc(N, sizeof(float));
  graph->Wm = calloc(N, sizeof(float));
  graph->W0 = calloc(N, sizeof(float));

  graph->sigma_v = calloc(N, sizeof(float));
  graph->sigma_f = calloc(M, sizeof(float));

  // read clause connections
  for (int a = 0; a < inst->M; a++) {
    int k = inst->c[a].k;
    graph->fk[a] = k;

    // init edge storage
    // pi*[a][j] is weight (i=fe[j])->a
    graph->fe[a] = calloc(k, sizeof(int));
    graph->fJ[a] = calloc(k, sizeof(int));

    graph->pi_u[a] = calloc(k, sizeof(float));
    graph->pi_s[a] = calloc(k, sizeof(float));
    graph->pi_0[a] = calloc(k, sizeof(float));

    for (int j = 0; j < k; j++) {
      int l = inst->c[a].l[j];
      int i = l > 0 ? l - 1 : -l - 1;
      int J = l > 0 ? -1 : 1;

      // store edges
      graph->fe[a][j] = i;
      graph->fJ[a][j] = J;

      // attach factor to variable too
      int vk = ++graph->vk[i];

      graph->ve[i] = realloc(graph->ve[i], vk * sizeof(int));
      graph->vJ[i] = realloc(graph->vJ[i], vk * sizeof(int));

      // eta[i][b] is survey (a=ve[b])->i
      graph->eta[i] = realloc(graph->eta[i], vk * sizeof(int));

      graph->ve[i][vk - 1] = a;
      graph->vJ[i][vk - 1] = J;
    }
  }
}

// eta[i][a] : a'th factor on var i
// pi*[a][i] : i'th var on factor a
void update_pi(struct fg *graph) {
  float **eta = graph->eta;
  float **pi_u = graph->pi_u;
  float **pi_s = graph->pi_s;
  float **pi_0 = graph->pi_0;

  for (int a = 0; a < graph->M; a++) {
    for (int j = 0; j < graph->fk[a]; i++) {
      // calculate new pi*[a][i], pi*_{i -> a}
      int i = graph->fe[a][j];

      // calculate products separately
      float pu = 1;
      float ps = 1;
      float p0 = 1;
      for (int b = 0; b < graph->vk[i]; b++) {
        if (a != graph->ve[i][b]) { // u_{b -> i} for a != b
          if (graph->vJ[i][b] != graph->fJ[a][j]) {
            // different edge type
            pu *= 1 - eta[i][b];
          } else {
            // same edge type
            ps *= 1 - eta[i][b];
          }
          // all edges
          p0 *= 1 - eta[i][b];
        }
      }

      pi_u[a][j] = ps - p0;
      pi_s[a][j] = pu - p0;
      pi_0[a][j] = p0;

      assert(0 <= pi_u[a][j] && pi_u[a][j] <= 1);
      assert(0 <= pi_s[a][j] && pi_s[a][j] <= 1);
      assert(0 <= pi_0[a][j] && pi_0[a][j] <= 1);
    }
  }

}

// return: 1 if converged, 0 otherwise
// new eta is placed into **eta
int update_eta(struct fg *graph, float **eta) {
  float **pi_u = graph->pi_u;
  float **pi_s = graph->pi_s;
  float **pi_0 = graph->pi_0;

  for (int i = 0; i < graph->N; i++) {
    for (int b = 0; b < graph->vk[i]; b++) {
      // calculate new eta[i][a], eta_{a -> i}
      int a = graph->ve[i][b];

      eta[i][b] = 1;
      for (int j = 0; j < graph->fk[a]; j++) {
        if (i != graph->fk[a][j]) {
          float denom = graph->pi_u[a][j] + graph->pi_s[a][j] + graph->pi_0[a][j];
          if (pi_u[a][j] == 0 && denom == 0) {
            eta[i][b] = 0;
          } else {
            eta[i][b] *= pi_u[a][j] / denom;
          }
        }
      }
      assert(0 <= eta[i][b] && eta[i][b] <= 1);

      if (absd(eta[i][b] - graph->eta[i][b]) > TOLERANCE) {
        *converged = 0;
      }
    }
  }

  return converged;
}

// return number of iters
int update_surveys(struct graph *graph, int max_iters) {
  int iters;
  int converged = 0;

  // create survey storage
  // eta[i][a] is survey a->i
  float **eta = calloc(graph->N, sizeof(float*));
  for (int i = 0; i < graph->N; i++) {
    eta[i] = calloc(graph->vk[i], sizeof(float));
  }

  for (iters = 0; iters < max_iters && !converged; iters++) {
    // compute fields pi*
    update_pi(graph);

    // compute new surveys eta
    converged = update_eta(&graph, eta);

    // swap in new surveys
    float **tmp = graph->eta;
    graph->eta = eta;
    eta = tmp;
  }

  for (int i = 0; i < graph->N; i++) {
    free(eta[i]);
  }
  free(eta);

  return iters;
}

// calculate bias and also complexities
void update_biases(struct fg *graph) {
  float **pi_u = graph->pi_u;
  float **pi_s = graph->pi_s;
  float **pi_0 = graph->pi_0;

  // calculate clause complexities
  graph->sigma = 0;
  for (int a = 0; a < graph->M; a++) {
    float prod_u = 1;
    float prod_denom = 1;
    for (int i = 0; i < graph->fk[a]; i++) {
      prod_u *= pi_u[a][i];
      prod_denom *= pi_u[a][i] + pi_s[a][i] + pi_0[a][i];
    }

    graph->sigma += log(prod_denom - prod_u);
    graph->sigma_f[a] = log(prod_denom - prod_u);
  }

  // compute biases
  for (int i = 0; i < graph->N; i++) {
    float pp = 1;
    float pm = 1;
    float p0 = 1;

    for (int b = 0; b < graph->vk[i]; b++) {
      if (graph->vJ[i][b] == -1) {
        pp *= 1 - graph->eta[i][b];
      } else {
        pm *= 1 - graph->eta[i][b];
      }
      // all edges
      p0 *= 1 - graph->eta[i][b];
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
    sigma -= (graph->vk[i] - 1) * log(pi_p + pi_m + pi_0);
    sigma_i[i] = log(pi_p + pi_m + pi_0);
    assert(sigma_i[i] <= 0);

    assert(0 <= Wp[i] && Wp[i] <= 1);
    assert(0 <= Wm[i] && Wm[i] <= 1);
    assert(0 <= W0[i] && W0[i] <= 1);
  }
}

// fix n variables according to polarization
void fix_vars(struct fg *graph, int n) {

}

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

  // convert instance to internal representation
  struct fg graph;
  fg_init(&graph, &inst);

  // randomly initialize surveys
  for (int i = 0; i < graph.N; i++) {
    for (int b = 0; b < graph.vk[i]; b++) {
      graph->eta[i][b] = (float) urand(1 << 30) / ((1 << 30) - 1);
      assert(0 <= graph->eta[i][b] && graph->eta[i][b] <= 1);
    }
  }

  // XXX parallel update; is the random order important?
  int steps = 0;
  int tot_iters = 0;
  while (graph.N != 0) {
    //graph_check(&graph);
    assert(graph.N > 0);
    assert(graph.M >= 0);

    // update surveys for max_iters timesteps
    int iters = update_surveys(&graph, max_iters, eta);

    tot_iters += iters - 1;

    printf("survey: %d steps %d iters\n", steps, iters);

    // check for convergence
    if (iters == max_iters) {
      printf("survey: unconverged\n");

      break;
    }

    printf("survey: converged\n");

    // check for trivial surveys
    for (int i = 0; i < graph->N; i++) {
      for (int b = 0; b < graph->vk[i]; b++) {
        if (absd(graph->eta[i][b]) > TOLERANCE) {
          goto non_trivial;
        }
      }
    }

    // all surveys zero: start walking
    printf("survey: surveys trivial\n");

    int *v;
    int *c;
    int walk_steps = walk(max_steps, p_param, PRINT_FREQ*10, &graph, &v, &c);

    printf("survey: %d walksteps\n", walk_steps);

    // save walk settings
    for (int i = 0; i < graph.N; i++) {
      orig_v[orig_vi[i]] = v[i] * 2 - 1;
    }

    // process sat/unsat here
    if (walk_steps == -1) {
      printf("survey: unsat\n");
    } else if (walk_steps < max_steps) {
      printf("survey: sat\n");

      // output sat assignment
      for (int i = 0; i < N; i++) {
        assert(orig_v[i] != 0);
        printf("%d ", (orig_v[i] + 1) / 2);
      }
      printf("\n");
    } else {
      printf("survey: unknown\n");
    }

    // clean up walk
    free(v);
    free(c);

    break;

non_trivial:

    update_biases(&graph);

    printf("survey: sigma/N %.*e\n", FLT_DIG, graph->sigma / graph->N);

    // fix variables according to biases and simplify clauses
    fix_vars(&graph, 1);

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
    for (int a = 0; a < graph.v[fix].k; a++) {
      int b = graph.v[fix].f[a].f->a;
      if (f[b] != 1 && graph.v[fix].f[a].j * dir == -1) {
        ff++;
        f[b] = 1;
      }
      if (f[b] * dir == 1) {
        would_ff++;
      }
    }
    //printf("survey: set %d to %d sat %d (vs. %d)\n", fix, dir, ff, would_ff);
    printf("survey: fixed %d to %d (k=%d), sats %d (opp sat %d)\n", fix, dir, graph.v[fix].k, ff, would_ff);

    // find unfixed factors
    int *nf = calloc(graph.M, sizeof(int)); // new clause index
    int fa = 0;
    for (int a = 0; a < graph.M; a++) {
      if (f[a] == 0) {
        nf[a] = fa++;
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
    for (int i = 0; i < ngraph.N; i++) {
      ngraph.v[i].i = i;
    }

    // create storage for unfixed clauses
    ngraph.f = calloc(ngraph.M, sizeof(struct node_f));
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
            int fk = ++ngraph.f[na].k;
            ngraph.f[na].v = realloc(ngraph.f[na].v, fk * sizeof(struct edge));
            ngraph.f[na].v[fk - 1] = e;

            // copy old eta
            eta[ni] = realloc(eta[ni], vk * sizeof(float));
            eta[ni][vk - 1] = p_eta[i][a];
            edges++;
          }
        }

        ni++;
      }
    }

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

  printf("survey: %d iters total\n", tot_iters);

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

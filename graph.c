// create factor graph from instance
// don't run on malformed instance

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "inst.h"
#include "graph.h"

// build graph out of an instance
void graph_make(struct graph *graph, struct inst *inst) {
  int N = inst->N;
  int M = inst->M;
  struct clause *c = inst->c;

  // create variables and set indices
  struct node_v *v = calloc(N, sizeof(struct node_v));
  for (int i = 0; i < N; i++) {
    v[i].i = i;
  }

  // create clauses and set indices
  struct node_f *f = calloc(M, sizeof(struct node_f));
  for (int a = 0; a < M; a++) {
    f[a].a = a;
  }

  // add edges for each clause
  for (int a = 0; a < M; a++) {
    int k = c[a].k;
    int *l = c[a].l;

    // initialize clause var array
    f[a].k = k;
    f[a].v = calloc(k, sizeof(struct edge));

    // add an edge for each variable
    for (int b = 0; b < k; b++) {
      int j = l[b] > 0 ? -1 : 1;        // +1 if negated, -1 o/w
      int i = l[b] > 0 ? l[b] : -l[b];  // variable index
      i--;  // variables in insts are one-indexed

      // create edge
      struct edge e;
      e.j = j;
      e.v = &v[i];
      e.f = &f[a];

      // add edge to variable
      v[i].k++;
      v[i].f = realloc(v[i].f, v[i].k * sizeof(struct edge));
      v[i].f[v[i].k - 1] = e;

      // add edge to clause
      f[a].v[b] = e;
    }
  }

  graph->N = N;
  graph->M = M;
  graph->v = v;
  graph->f = f;
}

void graph_show(struct graph *graph) {
  printf("VARIABLES (%d):\n", graph->N);
  for (int i = 0; i < graph->N; i++) {
    printf("%d: ", graph->v[i].i);
    for (int e = 0; e < graph->v[i].k; e++) {
      printf("(%d, %d) ", graph->v[i].f[e].j, graph->v[i].f[e].f->a);
    }
    printf("\n");
  }
  printf("\n");

  printf("CLAUSES (%d):\n", graph->M);
  for (int a = 0; a < graph->M; a++) {
    printf("%d: ", graph->f[a].a);
    for (int e = 0; e < graph->f[a].k; e++) {
      printf("(%d, %d) ", graph->f[a].v[e].j, graph->f[a].v[e].v->i);
    }
    printf("\n");
  }
}

void graph_free(struct graph *graph) {
  for (int i = 0; i < graph->N; i++) {
    free(graph->v[i].f);
  }
  for (int a = 0; a < graph->M; a++) {
    free(graph->f[a].v);
  }
  free(graph->v);
  free(graph->f);
}

// sanity checking a graph:
// - check variable and clause numbering (assume in-order for now)
// - check edges are connected right
void graph_check(struct graph *graph) {
  for (int i = 0; i < graph->N; i++) {
    assert(graph->v[i].i == i);
  }

  for (int a = 0; a < graph->M; a++) {
    assert(graph->f[a].a == a);
  }

  // create edge-check memory for variable side
  int **edge = calloc(graph->N, sizeof(int*));
  for (int i = 0; i < graph->N; i++) {
    edge[i] = calloc(graph->v[i].k, sizeof(int));
  }

  // mark edges in factor side appearing in variable side
  for (int a = 0; a < graph->M; a++) {
    for (int e = 0; e < graph->f[a].k; e++) {
      // search for this edge in the connected clause
      struct node_v *v = graph->f[a].v[e].v;
      int ok = 0;
      for (int b = 0; b < v->k; b++) {
        if (edge[v->i][b] == 0 && v->f[b].f->a == a && v->f[b].j == graph->f[a].v[e].j) {
          edge[v->i][b] = 1; // edge is accounted for; no reuse
          ok = 1;
          break;
        }
      }
      assert(ok);
    }
  }
  
  // all variable edges should be marked
  for (int i = 0; i < graph->N; i++) {
    for (int a = 0; a < graph->v[i].k; a++) {
      assert(edge[i][a]);
    }
    free(edge[i]);
  }
  free(edge);
}

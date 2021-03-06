struct inst;
struct clause;

struct edge {
  int j;            // coefficient: +1 if negated, -1 o/w
  struct node_v *v; // variable on connection
  struct node_f *f; // clause on connection
};

struct node_v {
  int i;            // variable ID
  int k;            // number of clauses
  struct edge *f;   // attached clause list
};

struct node_f {
  int a;            // clause ID
  int k;            // number of variables
  struct edge *v;   // attached variable list
};

struct graph {
  int N;            // number of variables
  int M;            // number of clauses
  struct node_v *v; // variable list
  struct node_f *f; // clause list
};

void graph_make(struct graph *graph, struct inst *inst);  // make graph from inst
void graph_show(struct graph *graph);   // print dump of a graph
void graph_free(struct graph *graph);   // free the contents of a graph
void graph_check(struct graph *graph);  // sanity check a graph

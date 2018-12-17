struct clause {
  int k;    // number of variables in this clause
  int *l;   // literal list
}

struct inst {
  int N;            // number of variables
  int M;            // number of clauses
  struct clause *c; // clause list
}

int parse(struct inst *inst);       // parse an inst from stdin
void free_inst(struct inst *inst);  // free the contents of an inst

struct clause {
  int k;    // number of variables in this clause
  int *l;   // literal list
};

struct inst {
  int N;            // number of variables
  int M;            // number of clauses
  struct clause *c; // clause list
};

int inst_parse(struct inst *inst);  // parse an inst from stdin
void inst_show(struct inst *inst);  // print dump of an inst
void inst_free(struct inst *inst);  // free the contents of an inst

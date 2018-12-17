// simple sat solve via walking
// provide the instance on stdin
#include <stdio.h>
#include <stdlib.h>

#include "inst.h"
#include "graph.h"

int main(int argc, char *argv[]) {
  // parse input instance
  struct inst inst;
  inst_parse(&inst);

  // generate factor graph
  struct graph graph;
  graph_make(&graph, &inst);

  // dump for debug
  inst_show(&inst);
  graph_show(&graph);

  // clean up
  inst_free(&inst);
  graph_free(&graph);

  return 0;
}

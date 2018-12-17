// simple sat solve via walking
#include <stdio.h>
#include <stdlib.h>

#include "parse.h"

int main(int argc, char *argv[]) {
  struct inst inst;
  parse(&inst);
  show_inst(&inst);
  free_inst(&inst);
}

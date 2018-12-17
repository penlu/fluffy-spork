/* Generate a random k-sat instance. Invoke as:
 * ./randsat k N M [seed]
 * e.g.:
 * ./randsat 3 128 64 $(date +%s)
 *
 * Random k-sat: N variables, M clauses.
 * Clauses each contain k literals drawn uniformly at random
 * from the set of 2N literals.
 *
 * Output format:
 * [N] - the number of variables
 * [M] - the number of clauses
 * [lits...] - literals in first clause
 * ...
 * Literals are numbers, i for some variable and -i for its negation.
 * Variable index will be in [1, N].
 */

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "util.h"

int main(int argc, char *argv[]) {
  if (argc != 4 && argc != 5) {
    printf("./randsat k N M [seed]\n");
    exit(0);
  }

  // get parameters
  int k = strtol(argv[1], NULL, 0);
  int N = strtol(argv[2], NULL, 0);
  int M = strtol(argv[3], NULL, 0);

  // set random seed
  if (argc == 4) {
    // use current time
    srand(time(NULL));
  } else if (argc == 5) {
    // use user-provided seed
    srand(strtol(argv[4], NULL, 0));
  }

  // generate instance
  printf("%d\n", N);
  printf("%d\n", M);
  for (int i = 0; i < M; i++) {
    for (int j = 0; j < k; j++) {
      // random variable
      int v = urand(N) + 1;

      // randomly negated
      if (rand() % 2) {
        printf("%d ", v);
      } else {
        printf("-%d ", v);
      }
    }
    printf("\n");
  }

  return 0;
}

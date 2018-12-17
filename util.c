#include <stdlib.h>
#include <assert.h>

#include "util.h"

// generate a uniform random integer [0, max)
int urand(int max) {
  while (1) {
    int r = rand();
    // throw away values that don't go into the interval evenly
    if (r < (((unsigned long long) RAND_MAX + 1) / max) * max) {
      int ret = r % max;
      assert(0 <= ret && ret < max);
      return ret;
    }
  }
}

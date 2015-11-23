
#include "random.h"

float unif(float a,float b) {
  return a + (b - a) * rand() / RAND_MAX;
}

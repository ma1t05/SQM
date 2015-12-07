
#include "random.h"
#include <cmath>

float unif(float a,float b) {
  return a + (b - a) * rand() / RAND_MAX;
}

int unif(int a) {
  return floor(double(a) * rand() / RAND_MAX);
}

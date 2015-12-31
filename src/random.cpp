
#include "random.h"

float unif(float a,float b) {
  return a + (b - a) * rand() / RAND_MAX;
}

int unif(int a) {
  return floor(double(a) * rand() / RAND_MAX);
}

double unif() {
  return (double) rand() / RAND_MAX;
}

double exponential(double lambda) {
  return -log(unif())/lambda;
}

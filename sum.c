#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <unistd.h>

#ifndef NTIMES
#define NTIMES 10
#endif

#ifndef ITER
#define ITER 1
#endif

#ifndef SUM_ARRAY_SIZE
#define SUM_ARRAY_SIZE 10000000
#endif

#ifndef OFFSET
#define OFFSET 0
#endif

#ifndef STREAM_TYPE
#define STREAM_TYPE double
#endif

#ifndef NUM_TESTS
#define NUM_TESTS 4
#endif

static STREAM_TYPE a[SUM_ARRAY_SIZE + OFFSET], b[SUM_ARRAY_SIZE + OFFSET],
    c[SUM_ARRAY_SIZE + OFFSET];
static double avgtime[NUM_TESTS] = {0}, maxtime[NUM_TESTS] = {0},
              mintime[NUM_TESTS] = {FLT_MAX};

extern int omp_get_num_threads();
extern int omp_get_thread_num();

double mysecond() {
  struct timeval tp;
  struct timezone tzp;
  int i;

  i = gettimeofday(&tp, &tzp);
  return ((double)tp.tv_sec + (double)tp.tv_usec * 1.e-6);
}

#ifndef MIN
#define MIN(x, y) ((x) < (y) ? (x) : (y))
#endif
#ifndef MAX
#define MAX(x, y) ((x) > (y) ? (x) : (y))
#endif

#define M 20

int checktick() {
  int i, minDelta, Delta;
  double t1, t2, timesfound[M];

  /*  Collect a sequence of M unique time values from the system. */

  for (i = 0; i < M; i++) {
    t1 = mysecond();
    while (((t2 = mysecond()) - t1) < 1.0E-6)
      ;
    timesfound[i] = t1 = t2;
  }

  /*
   * Determine the minimum difference between these M values.
   * This result will be our estimate (in microseconds) for the
   * clock granularity.
   */

  minDelta = 1000000;
  for (i = 1; i < M; i++) {
    Delta = (int)(1.0E6 * (timesfound[i] - timesfound[i - 1]));
    minDelta = MIN(minDelta, MAX(Delta, 0));
  }

  return (minDelta);
}

void generate_random(double *input, size_t size) {
  for (size_t i = 0; i < size; i++) {
    input[i] = rand() / (double)(RAND_MAX);
  }
}

double serial_sum(double *x, size_t size) {
  double sum_val = 0.0;

  for (size_t i = 0; i < size; i++) {
    sum_val += x[i];
  }

  return sum_val;
}

double omp_sum(double *x, size_t size) {
  double sum_val = 0.0;
#pragma omp parallel for
  for (size_t i = 0; i < size; i++) {
    sum_val += x[i];
  }

  return sum_val;
}

double omp_critical_sum(double *x, size_t size) {
  double sum_val = 0.0;
#pragma omp parallel for
  for (size_t i = 0; i < size; i++) {
#pragma omp critical
    sum_val += x[i];
  }

  return sum_val;
}

double omp_local_sum(double *x, size_t size) {
  int MAX_THREADS;
#pragma omp parallel
  {
#pragma omp master
    { MAX_THREADS = omp_get_num_threads(); }
  }

  double sum_val = 0.0;
  double local_sum[MAX_THREADS];

#pragma omp parallel shared(local_sum)
  {
    int id = omp_get_thread_num();
    local_sum[id] = 0;
#pragma omp for
    for (size_t i = 0; i < size; i++)
      local_sum[id] += x[i];
  }

  for (int i = 0; i < MAX_THREADS; i++) {
    sum_val += local_sum[i];
  }

  return sum_val;
}

double omp_opt_sum(double *x, size_t size) {
  typedef struct {
    double val;
    char pad[56];
  } tsum;
  int MAX_THREADS;
#pragma omp parallel
  {
#pragma omp master
    { MAX_THREADS = omp_get_num_threads(); }
  }

  double sum_val = 0.0;
  tsum local_sum[MAX_THREADS];

#pragma omp parallel shared(local_sum)
  {
    int id = omp_get_thread_num();
    local_sum[id].val = 0;
#pragma omp for
    for (size_t i = 0; i < size; i++)
      local_sum[id].val += x[i];
  }

  for (int i = 0; i < MAX_THREADS; i++) {
    sum_val += local_sum[i].val;
  }

  return sum_val;
}

int main() {
  int quantum, checktick(), k;
  double times[NUM_TESTS][NTIMES];
  if ((quantum = checktick()) >= 1)
    printf("Your clock granularity/precision appears to be "
           "%d microseconds.\n",
           quantum);
  else {
    printf("Your clock granularity appears to be "
           "less than one microsecond.\n");
    quantum = 1;
  }

#pragma omp parallel
  {
#pragma omp master
    {
      k = omp_get_num_threads();
      printf("Number of Threads requested = %i\n", k);
    }
  }

  generate_random(a, SUM_ARRAY_SIZE);
  double result[NUM_TESTS];
  for (k = 0; k < NTIMES; k++) {
    times[0][k] = mysecond();
    for (int i = 0; i < ITER; i++) {
      result[0] = serial_sum(a, SUM_ARRAY_SIZE);
    }
    times[0][k] = mysecond() - times[0][k];
  }

  for (k = 0; k < NTIMES; k++) {
    times[1][k] = mysecond();
    for (int i = 0; i < ITER; i++) {
      result[1] = omp_sum(a, SUM_ARRAY_SIZE);
    }
    times[1][k] = mysecond() - times[1][k];
  }

  for (k = 0; k < NTIMES; k++) {
    times[2][k] = mysecond();
    for (int i = 0; i < ITER; i++) {
      result[2] = omp_critical_sum(a, SUM_ARRAY_SIZE);
    }
    times[2][k] = mysecond() - times[2][k];
  }

  for (k = 0; k < NTIMES; k++) {
    times[3][k] = mysecond();
    for (int i = 0; i < ITER; i++) {
      result[3] = omp_local_sum(a, SUM_ARRAY_SIZE);
    }
    times[3][k] = mysecond() - times[3][k];
  }

  for (int i = 0; i < NUM_TESTS; i++) {
    for (k = 1; k < NTIMES; k++) /* note -- skip first iteration */
    {
      avgtime[i] = avgtime[i] + times[i][k];
      mintime[i] = MIN(mintime[i], times[i][k]);
      maxtime[i] = MAX(maxtime[i], times[i][k]);
    }
  }
  for (int i = 0; i < NUM_TESTS; i++) {
    avgtime[i] = avgtime[i] / (double)(NTIMES - 1);
    printf("result = %f\n", result[i]);
    printf("avgtime = %f s\n", avgtime[i]);
    printf("mintime = %f s\n", mintime[i]);
    printf("maxtime = %f s\n", maxtime[i]);
  }
}
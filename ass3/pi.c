
#include <math.h>
#include <mpi.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define SEED 921
#define NUM_ITER 1000000000

int main(int argc, char *argv[]) {
  int count = 0;
  int local_count = 0;
  double x, y, z, pi;
  int rank, size, provided, flips;
  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);

  double start, elapsed;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  start = MPI_Wtime();

  // Important: Multiply SEED by "rank" when you introduce MPI!
  srand(SEED * rank);
  flips = NUM_ITER / size;

  // Calculate PI following a Monte Carlo method
  for (int iter = 0; iter < flips; iter++) {
    // Generate random (X,Y) points
    x = (double)random() / (double)RAND_MAX;
    y = (double)random() / (double)RAND_MAX;
    z = sqrt((x * x) + (y * y));

    // Check if point is in unit circle
    if (z <= 1.0) {
      local_count++;
    }
  }
  MPI_Reduce(&local_count, &count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  // Estimate Pi and display the result
  elapsed = MPI_Wtime() - start;
  if (rank == 0) {
    pi = ((double)count / (double)NUM_ITER) * 4.0;
    printf("The result is %f\n Elapsed time %f\n", pi, elapsed);
  }

  return 0;
}

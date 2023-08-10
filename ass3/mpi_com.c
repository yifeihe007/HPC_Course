#include <mpi.h>
#include <stdio.h>

int main(int argc, char *argv[]) {
  int rank, size, i, provided;
  float A[10];

  MPI_Init_thread(&argc, &argv, MPI_THREAD_SINGLE, &provided);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  // only rank zero send everybody else receiving from rank 0
  if (rank == 0) {
    for (i = 0; i < 10; i++)
      A[i] = i;
    for (i = 1; i < size; i++) // note that we start from rank 1
      MPI_Send(&A, 10, MPI_FLOAT, i, 0, MPI_COMM_WORLD);
  } else {
    MPI_Recv(&A, 10, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }

  printf("My rank %d of %d\n", rank, size);
  printf("Here are my values for A\n");
  for (i = 0; i < 10; i++)
    printf("%f ", A[i]);
  printf("\n");
  MPI_Finalize();
}

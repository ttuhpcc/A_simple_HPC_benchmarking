#include <stdio.h>
#include <stdlib.h>
#include <fftw3.h>
#include <fftw3-mpi.h>
#include <math.h>
#include <time.h>

#define MAX_SIZE 1024 // Maximum size of the FFT

void display_time_and_flops(double time_taken, double flops) {
    if (time_taken < 1e-6) {
        printf("Time taken for FFT: %.2f µs\n", time_taken * 1e6);
    } else if (time_taken < 1e-3) {
        printf("Time taken for FFT: %.2f ms\n", time_taken * 1e3);
    } else {
        printf("Time taken for FFT: %.2f s\n", time_taken);
    }

    if (flops < 1e6) {
        printf("Performance: %.2f FLOPS\n", flops);
    } else if (flops < 1e12) {
        printf("Performance: %.2f GFLOPS\n", flops * 1e-9);
    } else {
        printf("Performance: %.2f TFLOPS\n", flops * 1e-12);
    }
}

void perform_mpi_fft(const char *input_file) {
    int N, rank, size;
    fftw_complex *in, *out;
    int local_size;

    // Initialize MPI
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Read input size from file only on rank 0
    if (rank == 0) {
        FILE *file = fopen(input_file, "r");
        if (!file) {
            perror("Failed to open input file");
            MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
        }
        fscanf(file, "%d", &N);
        fclose(file);
    }

    // Broadcast the size to all processes
    MPI_Bcast(&N, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (N > MAX_SIZE || N % size != 0) {
        if (rank == 0) {
            fprintf(stderr, "Error: Maximum size exceeded or N is not divisible by number of processes\n");
        }
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Calculate local size (number of elements per process)
    local_size = N / size;

    // Allocate input and output arrays
    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * local_size);
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * local_size);

    if (!in || !out) {
        fprintf(stderr, "Failed to allocate memory\n");
        MPI_Abort(MPI_COMM_WORLD, EXIT_FAILURE);
    }

    // Initialize the input array with a linear sequence
    for (int i = 0; i < local_size; i++) {
        in[i][0] = (double)(rank * local_size + i); // Real part
        in[i][1] = 0.0;                              // Imaginary part
    }

    // Create a plan for FFTW
    fftw_plan p = fftw_plan_dft_1d(local_size, in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    // Start timing
    MPI_Barrier(MPI_COMM_WORLD);
    double start = MPI_Wtime();

    // Execute the FFT
    fftw_execute(p);

    // Stop timing
    double end = MPI_Wtime();

    // Calculate the time taken and FLOPS
    double time_taken = end - start;
    double flops = (2.0 * N * log2(N)) / time_taken;

    // Gather results to rank 0
    fftw_complex *gathered_out = NULL;
    if (rank == 0) {
        gathered_out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N);
    }

    // Use MPI_Gather to combine results
    MPI_Gather(out, local_size * sizeof(fftw_complex), MPI_BYTE,
               gathered_out, local_size * sizeof(fftw_complex), MPI_BYTE, 0, MPI_COMM_WORLD);

    // Print results on rank 0
    if (rank == 0) {
        printf("FFT Results for size %d:\n", N);
        for (int i = 0; i < N; i++) {
            printf("%d: %f + %fi\n", i, gathered_out[i][0], gathered_out[i][1]);
        }
        display_time_and_flops(time_taken, flops);
        fftw_free(gathered_out);
    }

    // Clean up
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    MPI_Finalize();
}

int main(int argc, char **argv) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <input_file.dat>\n", argv[0]);
        return EXIT_FAILURE;
    }
    perform_mpi_fft(argv[1]);
    return EXIT_SUCCESS;
}

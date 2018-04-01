#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mpi.h"
#include "spec.h"
#include "pfile.h"
#include "pool.h"

Pool *pool;
CompSpec *cs;

void init_pool(int rank, const Spec *spec) {
    if (rank == 0) {
        pool->size = spec->gs.size;
        pool->small_num = spec->gs.small_num;
        pool->small_mass = spec->gs.small_mass;
        pool->small_rad = spec->gs.small_rad;

        int num = pool->small_num;
        pool->small_ptc = (Location *)calloc(num, sizeof(Location));
        srand(time(NULL));
        for (int i = 0, max = pool->size; i < num; i++) {
            pool->small_ptc[i].x = ((float)rand() / RAND_MAX) * max;
            pool->small_ptc[i].y = ((float)rand() / RAND_MAX) * max;
        }

        num = pool->large_num = spec->gs.large_num;
        pool->large_ptc = (Particle *)calloc(num, sizeof(Particle));
        memcpy(pool->large_ptc, spec->gs.large_ptc, num*sizeof(Particle));
    }
}

int init_params(int argc, char *argv[], Spec *spec) {
    // check the number of arguments
    if (argc != 3) {
        fprintf(stderr, "Unknown arguments...\n\n");
        fprintf(stderr, "Usage: pool <specfile> <ppmfile>\n");
        return EXIT_FAILURE;
    }
    int ret = read_spec_from_file(argv[1], spec);
    cs = (CompSpec *)malloc(sizeof(CompSpec));
    memcpy(cs, &spec->cs, sizeof(CompSpec));
    return ret;
}

int run(int rank, int size, int argc, char *argv[]) {
    Spec *spec = (Spec *)malloc(sizeof(Spec));
    if (rank == 0 && init_params(argc, argv, spec)){
        fprintf(stderr, "[%s:%d] init_params failed\n", __FILE__, __LINE__);
        return -1;
    }
    pool = (Pool *)malloc(sizeof(Pool));
    init_pool(rank, spec);
    free(spec->gs.large_ptc);
    free(spec);
    // output results to ppm file
    if (print2ppm(pool, argv[2])) {
        fprintf(stderr, "[%s:%d] print2ppm failed\n", __FILE__, __LINE__);        
        return -1;
    }
    return 0;
}

int main(int argc, char *argv[]) {
	int size, rank;
	MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int ret = run(rank, size, argc, argv);

    MPI_Finalize();
    return ret;
}
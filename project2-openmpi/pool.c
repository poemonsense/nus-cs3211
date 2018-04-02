#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mympi.h"
#include "spec.h"
#include "pfile.h"
#include "pool.h"
#include "logging.h"

Pool *pool;
CompSpec *cs;

/**********
 * Debug function definitions
 **********/
void __debug_print_pool(int rank, const Pool *pool) {
    INFO("Process %d: Pool at address 0x%lx", rank, (uint64_t)pool);
    INFO("size:       %d", pool->size);
    INFO("small_num:  %d", pool->small_num);
    INFO("small_mass: %.6f", pool->small_mass);
    INFO("small_rad:  %.6f", pool->small_rad);
    INFO("large_num:  %d", pool->large_num);
    INFO("small_ptc:  0x%lx", (uint64_t)pool->small_ptc);
    INFO("large_ptc:  0x%lx", (uint64_t)pool->large_ptc);
}

/**
 * Main purpose function definitions
 */

void init_pool(int rank, const Spec *spec) {
    if (rank == 0) {
        pool->size = spec->gs.size;
        pool->small_num = spec->gs.small_num;
        pool->small_mass = spec->gs.small_mass;
        pool->small_rad = spec->gs.small_rad;
        // initialize small paticles
        int num = pool->small_num;
        pool->small_ptc = (Location *)calloc(num, sizeof(Location));
        srand(time(NULL));
        for (int i = 0, max = pool->size; i < num; i++) {
            pool->small_ptc[i].x = ((float)rand() / RAND_MAX) * max;
            pool->small_ptc[i].y = ((float)rand() / RAND_MAX) * max;
        }
        // initialize large paticles
        num = pool->large_num = spec->gs.large_num;
        pool->large_ptc = (Particle *)calloc(num, sizeof(Particle));
        memcpy(pool->large_ptc, spec->gs.large_ptc, num*sizeof(Particle));
    }
    mympi_bcast(pool, 1, PoolType, 0, MPI_COMM_WORLD);
    #ifdef POOL_DEBUG
    __debug_print_pool(rank, pool);
    #endif
    if (rank != 0) {
        pool->small_ptc = (Location *)calloc(pool->small_num, sizeof(Location));
        pool->large_ptc = (Particle *)calloc(pool->large_num, sizeof(Particle));        
    }
    mympi_bcast(pool->small_ptc, pool->small_num, LocationType, 0, MPI_COMM_WORLD);
    mympi_bcast(pool->large_ptc, pool->large_num, ParticleType, 0, MPI_COMM_WORLD);
    // #ifdef POOL_DEBUG
    // __debug_print_location(rank, pool->small_ptc, pool->small_num);
    // __debug_print_particle(rank, pool->large_ptc, pool->large_num);
    // #endif
}

int init_params(int rank, int argc, char *argv[], Spec *spec) {
    cs = (CompSpec *)malloc(sizeof(CompSpec));
    if (rank == 0) {
        // check the number of arguments
        if (argc != 3) {
            fprintf(stderr, "Unknown arguments...\n\n");
            fprintf(stderr, "Usage: pool <specfile> <ppmfile>\n");
            return -1;
        }
        if (read_spec_from_file(argv[1], spec)) {
            ERR("[%s:%d] read_spec_from_file failed\n", __FILE__, __LINE__);
            return -1;
        }
        memcpy(cs, &spec->cs, sizeof(CompSpec));
    }
    mympi_bcast(cs, 1, CompSpecType, 0, MPI_COMM_WORLD);
    #ifdef POOL_DEBUG
    __debug_print_compspec(rank, cs);
    #endif
    return 0;
}

int run(int rank, int size, int argc, char *argv[]) {
    // initialize the game (CompSpec and Pool)
    Spec *spec = (Spec *)malloc(sizeof(Spec));
    if (init_params(rank, argc, argv, spec)){
        ERR("[%s:%d] init_params failed\n", __FILE__, __LINE__);
        return -1;
    }
    pool = (Pool *)malloc(sizeof(Pool));
    init_pool(rank, spec);
    // only in rank 0 have we initialized large_ptc in spec
    if (rank == 0)
        free(spec->gs.large_ptc);
    free(spec);
    // output results to ppm file
    if (rank == 0 && print2ppm(pool, argv[2])) {
        ERR("[%s:%d] print2ppm failed\n", __FILE__, __LINE__);        
        return -1;
    }
    return 0;
}

int main(int argc, char *argv[]) {
	int size, rank;
    MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    time_t rawtime;
    time(&rawtime);
    struct tm *timeinfo = localtime(&rawtime);

    char logfile[64];
    sprintf(logfile, "log-%02d%02d%02d%02d%02d-%d.txt",
            timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec, rank);
    if (logging_open(logfile) < 0) {
        fprintf(stderr, "Failed to initialize logging\n");
        return EXIT_FAILURE;
    }

    init_mympi();

    int ret = run(rank, size, argc, argv);

    MPI_Finalize();
    return ret;
}

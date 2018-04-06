#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "mympi.h"
#include "spec.h"
#include "pfile.h"
#include "pool.h"

#ifdef POOL_DEBUG
#include "logging.h"
#endif

/**
 * Pool and Velocity describing the current state of the board
 */
Pool pool;
/**
 * Velocity info of the small and large particles
 * Separately defined since we won't transfer them between processes
 */
Velocity *small_vel, *large_vel;
/**
 * CompSpec for current computation
 */
CompSpec cs;
int adj_ranks[9];
int adj_num = 0;
int proc_size;

/**********
 * Debug function definitions
 **********/
#ifdef POOL_DEBUG

void __debug_print_pool(int rank, const Pool *pool) {
    DEBUG("Process %d: Pool at address 0x%lx", rank, (uint64_t)pool);
    DEBUG("size:       %d", pool->size);
    DEBUG("small_num:  %d", pool->small_num);
    DEBUG("small_mass: %.6f", pool->small_mass);
    DEBUG("small_rad:  %.6f", pool->small_rad);
    DEBUG("large_num:  %d", pool->large_num);
    DEBUG("small_ptc:  0x%lx", (uint64_t)pool->small_ptc);
    DEBUG("large_ptc:  0x%lx", (uint64_t)pool->large_ptc);
}

#endif

/**
 * Main purpose function definitions
 */
int __is_perfect_square(int num);
int run_step(int rank, int size);
void init_pool(int rank, const Spec *spec);
int init_params(int rank, int argc, char *argv[], Spec *spec);
int get_final_result(int rank, int size, PPMFile *ppm);

int get_final_result(int rank, int size, PPMFile *ppm) {
    uint64_t rgb_num = pool.size * pool.size;
    if (rank == 0) {
        // assume every node has the same size of grid
        ppm->size = pool.size * proc_size;
        ppm->pixels = (RGBPixel *)calloc(ppm->size * ppm->size, sizeof(RGBPixel));
        PPMFile temp;
        if (pool2ppm(&pool, &temp)) {
            #ifdef POOL_DEBUG
            INFO("pool2ppm for rank %d failed", rank)
            #endif
            return -1;
        }
        MPI_Request recv[size];
        memcpy(ppm->pixels, temp.pixels, sizeof(RGBPixel) * rgb_num);
        RGBPixel *p = ppm->pixels;
        for (int i = 1; i < size; i++) {
            p += rgb_num;
            mympi_irecv(p, rgb_num, RGBPixelType, i, i, &recv[i]);
        }
        mympi_waitall(size-1, recv+1, MPI_STATUSES_IGNORE);
    }
    else {
        if (pool2ppm(&pool, ppm)) {
            #ifdef POOL_DEBUG
            INFO("pool2ppm for rank %d failed", rank)
            #endif
            return -1;
        }
        mympi_send(ppm->pixels, rgb_num, RGBPixelType, 0, rank);
    }
    return 0;
}

void init_pool(int rank, const Spec *spec) {
    if (rank == 0) {
        pool.size = spec->gs.size;
        pool.small_num = spec->gs.small_num;
        pool.small_mass = spec->gs.small_mass;
        pool.small_rad = spec->gs.small_rad;
        // initialize small paticles
        int num = pool.small_num;
        pool.small_ptc = (Location *)calloc(num, sizeof(Location));
        srand(time(NULL));
        for (int i = 0, max = pool.size; i < num; i++) {
            pool.small_ptc[i].x = ((double)rand() / RAND_MAX) * max;
            pool.small_ptc[i].y = ((double)rand() / RAND_MAX) * max;
        }
        // initialize large paticles
        num = pool.large_num = spec->gs.large_num;
        pool.large_ptc = (Particle *)calloc(num, sizeof(Particle));
        memcpy(pool.large_ptc, spec->gs.large_ptc, num*sizeof(Particle));
    }
    mympi_bcast(&pool, 1, PoolType, 0);
    #ifdef POOL_DEBUG
    __debug_print_pool(rank, &pool);
    #endif
    if (rank != 0) {
        pool.small_ptc = (Location *)calloc(pool.small_num, sizeof(Location));
        pool.large_ptc = (Particle *)calloc(pool.large_num, sizeof(Particle));
    }
    mympi_bcast(pool.small_ptc, pool.small_num, LocationType, 0);
    mympi_bcast(pool.large_ptc, pool.large_num, ParticleType, 0);
    // #ifdef POOL_DEBUG
    // __debug_print_location(rank, pool.small_ptc, pool.small_num);
    // __debug_print_particle(rank, pool.large_ptc, pool.large_num);
    // #endif
}

int init_params(int rank, int argc, char *argv[], Spec *spec) {
    if (rank == 0) {
        // check the number of arguments
        if (argc != 3) {
            fprintf(stderr, "Unknown arguments...\n\n");
            fprintf(stderr, "Usage: pool <specfile> <ppmfile>\n");
            return -1;
        }
        if (read_spec_from_file(argv[1], spec)) {
            fprintf(stderr, "read_spec_from_file failed\n");
            return -1;
        }
        memcpy(&cs, &spec->cs, sizeof(CompSpec));
    }
    mympi_bcast(&cs, 1, CompSpecType, 0);
    #ifdef POOL_DEBUG
    __debug_print_compspec(rank, &cs);
    #endif
    return 0;
}

/**
 * Run one step of the computation
 * Note that currently small_num and large_num in any of the region should be
 * positive. If zero number occurs, the program may crash
 */
int run_step(int rank, int size) {
    // communicate with adjacent regions to send and receive data
    // use asynchronous communication primitives to avoid dead-locks
    // get number of small and large particles in adjacent regions
    #ifdef POOL_DEBUG
    INFO("Sending and receiving small_num and large_num, size: MPI_UINT32_T");
    #endif
    MPI_Request num_send_req[2][adj_num], num_recv_req[2][adj_num];
    uint32_t snum[adj_num], lnum[adj_num];
    for (int i = 0; i < adj_num; i++) {
        #ifdef POOL_DEBUG
        INFO("Process %d: Communicate with %d", rank, adj_ranks[i])
        #endif
        mympi_isend(&pool.small_num, 1, MPI_UINT32_T, adj_ranks[i], 2*rank,
            &num_send_req[0][i]);
        mympi_isend(&pool.large_num, 1, MPI_UINT32_T, adj_ranks[i], 2*rank+1,
            &num_send_req[1][i]);
        mympi_irecv(&snum[i], 1, MPI_UINT32_T, adj_ranks[i], 2*adj_ranks[i],
            &num_recv_req[0][i]);
        mympi_irecv(&lnum[i], 1, MPI_UINT32_T, adj_ranks[i], 2*adj_ranks[i]+1,
            &num_recv_req[1][i]);
    }
    Location *adj_small_ptc[adj_num];
    Particle *adj_large_ptc[adj_num];
    // allocate space for accelaration here
    #ifdef POOL_DEBUG
    INFO("Allocating space for accelaration");
    #endif
    Accel *small_acc = (Accel *)calloc(pool.small_num, sizeof(Accel));
    Accel *large_acc = (Accel *)calloc(pool.large_num, sizeof(Accel));
    memset(large_acc, 0, pool.large_num * sizeof(Accel));
    memset(small_acc, 0, pool.small_num * sizeof(Accel));
    // TODO: compare communication time between waitall and wait
    // assume waitall is faster since only 32bits transferred
    #ifdef POOL_DEBUG
    INFO("Wait for array of small_num and large_num to be received");
    #endif
    mympi_waitall(adj_num, num_recv_req[0], MPI_STATUSES_IGNORE);
    mympi_waitall(adj_num, num_recv_req[1], MPI_STATUSES_IGNORE);
    for (int i = 0; i < adj_num; i++) {
        #ifdef POOL_DEBUG
        INFO("snum[%d] = %d, snum[%d] = %d", i, snum[i], i, lnum[i]);
        #endif
        // mympi_wait(&num_recv_req[0][i], MPI_STATUS_IGNORE);
        adj_small_ptc[i] = (Location *)calloc(snum[i], sizeof(Location));
        // mympi_wait(&num_recv_req[1][i], MPI_STATUS_IGNORE);
        adj_large_ptc[i] = (Particle *)calloc(lnum[i], sizeof(Particle));
    }
    // send and recv particle info in adjacent regions
    #ifdef POOL_DEBUG
    for (int i = 0; i < adj_num; i++) {
        INFO("Sending and receiving small_ptc and large_ptc in region %d, "
             "size: %u LocationType and %u ParticleType",
            adj_ranks[i], snum[i], lnum[i]);
    }
    #endif
    MPI_Request ptc_send_req[2][adj_num], ptc_recv_req[2][adj_num];
    for (int i = 0; i < adj_num; i++) {
        mympi_isend(pool.small_ptc, pool.small_num, LocationType,
            adj_ranks[i], 2*rank, &ptc_send_req[0][i]);
        mympi_isend(pool.large_ptc, pool.large_num, ParticleType,
            adj_ranks[i], 2*rank+1, &ptc_send_req[1][i]);
        mympi_irecv(adj_small_ptc[i], snum[i], LocationType, adj_ranks[i],
            2*adj_ranks[i], &ptc_recv_req[0][i]);
        mympi_irecv(adj_large_ptc[i], lnum[i], ParticleType, adj_ranks[i],
            2*adj_ranks[i]+1, &ptc_recv_req[1][i]);
    }
    // compute impacts of the region itself and wait for communication to finish
    #ifdef POOL_DEBUG
    INFO("Computing the impacts of the region %d itself", rank)
    #endif
    // only compute F(i, j) sice F(j, i) = -F(i, j)
    for (int i = 0; i < pool.small_num; i++) {
        // forces between small particles
        for (int j = i + 1; j < pool.small_num; j++) {
            // compute ptc[j]'s impact on ptc[i]
            Location r_ji = {
                pool.small_ptc[j].x - pool.small_ptc[i].x,
                pool.small_ptc[j].y - pool.small_ptc[i].y
            };
            COMPUTE_ACCEL_2PTC_RAW(r_ji, pool.small_mass, pool.small_mass,
                small_acc[i], small_acc[j])
        }
        // forces between larger particles and small particles
        for (int j = 0; j < pool.large_num; j++) {
            Location r_ji = {
                pool.large_ptc[j].loc.x - pool.small_ptc[i].x,
                pool.large_ptc[j].loc.y - pool.small_ptc[i].y
            };
            COMPUTE_ACCEL_2PTC_RAW(r_ji, pool.small_mass, pool.large_ptc[j].mass,
                small_acc[i], large_acc[j])
        }
    }
    for (int i = 0; i < pool.large_num; i++) {
        for (int j = i + 1; j < pool.large_num; j++) {
            Location r_ji = {
                pool.large_ptc[j].loc.x - pool.large_ptc[i].loc.x,
                pool.large_ptc[j].loc.y - pool.large_ptc[i].loc.y
            };
            COMPUTE_ACCEL_2PTC_RAW(r_ji, pool.large_ptc[i].mass, pool.large_ptc[j].mass,
                large_acc[i], large_acc[j])
        }
    }
    // compute the impacts of adjacent regions
    for (int i = 0; i < adj_num; i++) {
        // compute the offset of two adjacent regions
        double offset[2] = {
            ((adj_ranks[i] % proc_size) - (rank % proc_size)) * (double)pool.size,
            ((adj_ranks[i] / proc_size) - (rank / proc_size)) * (double)pool.size            
        };
        #ifdef POOL_DEBUG
        INFO("compute region %d's impact, offset = { %.4f, %.4f }",
            adj_ranks[i], offset[0], offset[1])
        #endif
        // compute adj_large_ptc[i][j]'s impact
        // compute large ones first since there should be fewer large particles
        // which take less time to communicate
        mympi_wait(&ptc_recv_req[1][i], MPI_STATUS_IGNORE);
        for (int j = 0; j < lnum[i]; j++) {
            #ifdef POOL_DEBUG
            INFO("Process %d large_ptc[%d].loc: %.4f, %.4f", adj_ranks[i],
                j, adj_large_ptc[i][j].loc.x, adj_large_ptc[i][j].loc.y)
            #endif
            for (int k = 0; k < pool.small_num; k++) {
                // compute adj_large_ptc[i][j]'s impact on pool.small_ptc[k]
                Location r_kj = {
                    adj_large_ptc[i][j].loc.x - pool.small_ptc[k].x + offset[0],
                    adj_large_ptc[i][j].loc.y - pool.small_ptc[k].y + offset[1]
                };
                COMPUTE_ACCEL_RAW(r_kj, adj_large_ptc[i][j].mass, small_acc[k])
            }
            for (int k = 0; k < pool.large_num; k++) {
                // compute adj_large_ptc[i][j]'s impact on pool.large_ptc[k]                
                Location r_kj = {
                    adj_large_ptc[i][j].loc.x - pool.large_ptc[k].loc.x + offset[0],
                    adj_large_ptc[i][j].loc.y - pool.large_ptc[k].loc.y + offset[1]
                };
                COMPUTE_ACCEL_RAW(r_kj, adj_large_ptc[i][j].mass, large_acc[k])
            }
        }
        // compute adj_small_ptc[i][j]'s impact
        mympi_wait(&ptc_recv_req[0][i], MPI_STATUS_IGNORE);
        for (int j = 0; j < snum[i]; j++) {
            #ifdef POOL_DEBUG
            if (j < 10)
            INFO("Process %d small_ptc[%d]: %.4f, %.4f", adj_ranks[i],
                j, adj_small_ptc[i][j].x, adj_small_ptc[i][j].x)
            #endif
            for (int k = 0; k < pool.small_num; k++) {
                // compute adj_small_ptc[i][j]'s impact on pool.small_ptc[k]
                Location r_kj = {
                    adj_small_ptc[i][j].x - pool.small_ptc[k].x + offset[0],
                    adj_small_ptc[i][j].y - pool.small_ptc[k].y + offset[1]
                };
                COMPUTE_ACCEL_RAW(r_kj, pool.small_mass, small_acc[k])
            }
            for (int k = 0; k < pool.large_num; k++) {
                // compute adj_small_ptc[i][j]'s impact on pool.large_ptc[k]
                Location r_kj = {
                    adj_small_ptc[i][j].x - pool.large_ptc[k].loc.x + offset[0],
                    adj_small_ptc[i][j].y - pool.large_ptc[k].loc.y + offset[1]
                };
                COMPUTE_ACCEL_RAW(r_kj, pool.small_mass, large_acc[k])
            }
        }
    }
    // free the ptc space for other regions
    for (int i = 0; i < adj_num; i++) {
        free(adj_small_ptc[i]);
        free(adj_large_ptc[i]);
    }
    #ifdef POOL_DEBUG
    INFO("Apply G constant to Acceleration")
    #endif
    for (int i = 0; i < pool.small_num; i++)
        ACCEL_RAW_TO_FINAL(small_acc[i])
    for (int i = 0; i < pool.large_num; i++)
        ACCEL_RAW_TO_FINAL(large_acc[i])
    #ifdef POOL_DEBUG
    for (int i = 0; i < 10; i++) {
        INFO("small_acc[%d]: %.4f, %.4f", i, small_acc[i].ax, small_acc[i].ay);        
    }
    for (int i = 0; i < pool.large_num; i++) {
        INFO("large_acc[%d]: %.4f, %.4f", i, large_acc[i].ax, large_acc[i].ay);
    }
    #endif
    // now we are going to update the velocity and location
    double t = cs.time_step;
    #ifdef POOL_DEBUG
    INFO("Update velocity using vt = v0 + a * %.4f", t)
    #endif
    for (int i = 0; i < pool.small_num; i++) {
        small_vel[i].vx += small_acc[i].ax * t;
        small_vel[i].vy += small_acc[i].ay * t;
    }
    for (int i = 0; i < pool.large_num; i++) {
        large_vel[i].vx += large_acc[i].ax * t;
        large_vel[i].vy += large_acc[i].ay * t;
    }
    #ifdef POOL_DEBUG
    for (int i = 0; i < 10; i++) {
        INFO("small_vel[%d]: %.4f, %.4f", i, small_vel[i].vx, small_vel[i].vy);        
    }
    for (int i = 0; i < pool.large_num; i++) {
        INFO("large_vel[%d]: %.4f, %.4f", i, large_vel[i].vx, large_vel[i].vy);
    }
    #endif
    // need to ensure that particle information has been sent
    mympi_waitall(adj_num, ptc_send_req[0], MPI_STATUSES_IGNORE);
    mympi_waitall(adj_num, ptc_send_req[1], MPI_STATUSES_IGNORE);
    #ifdef POOL_DEBUG
    INFO("Update location using x = vt * t + 0.5 * a * t * t")
    #endif
    for (int i = 0; i < pool.small_num; i++) {
        pool.small_ptc[i].x += small_vel[i].vx * t + 0.5 * small_acc[i].ax * t * t;
        pool.small_ptc[i].y += small_vel[i].vy * t + 0.5 * small_acc[i].ay * t * t;
        while (pool.small_ptc[i].x > pool.size)
            pool.small_ptc[i].x -= pool.size;
        while (pool.small_ptc[i].x < 0)
            pool.small_ptc[i].x += pool.size;
        while (pool.small_ptc[i].y > pool.size)
            pool.small_ptc[i].y -= pool.size;
        while (pool.small_ptc[i].y < 0)
            pool.small_ptc[i].y += pool.size;
        #ifdef POOL_DEBUG
        if (i < 10) {
            INFO("pool.small_ptc[%d]: %.4f, %.4f", i, pool.small_ptc[i].x, pool.small_ptc[i].y);
        }
        #endif
    }
    for (int i = 0; i < pool.large_num; i++) {
        // TODO: currently if a particle runs out of the region, we make it run into
        // the region in the other direction
        pool.large_ptc[i].loc.x += large_vel[i].vx * t + 0.5 * large_acc[i].ax * t * t;
        pool.large_ptc[i].loc.y += large_vel[i].vy * t + 0.5 * large_acc[i].ay * t * t;
        while (pool.large_ptc[i].loc.x > pool.size)
            pool.large_ptc[i].loc.x -= pool.size;
        while (pool.large_ptc[i].loc.x < 0)
            pool.large_ptc[i].loc.x += pool.size;
        while (pool.large_ptc[i].loc.y > pool.size)
            pool.large_ptc[i].loc.y -= pool.size;
        while (pool.large_ptc[i].loc.y < 0)
            pool.large_ptc[i].loc.y += pool.size;
        #ifdef POOL_DEBUG
        INFO("pool.large_ptc[i].loc: %.4f, %.4f", pool.large_ptc[i].loc.x, pool.large_ptc[i].loc.y);
        #endif
    }
    // free the allocated acc space
    free(small_acc);
    free(large_acc);
    mympi_waitall(adj_num, num_send_req[0], MPI_STATUSES_IGNORE);
    mympi_waitall(adj_num, num_send_req[1], MPI_STATUSES_IGNORE);
    return 0;
}

int run(int rank, int size, int argc, char *argv[]) {
    // initialize the game (CompSpec and Pool)
    Spec *spec = (Spec *)malloc(sizeof(Spec));
    if (init_params(rank, argc, argv, spec)){
        fprintf(stderr, "init_params failed\n");        
        return -1;
    }
    init_pool(rank, spec);
    // free the allocated space
    // only in rank 0 have we initialized spec.gs in spec
    if (rank == 0)
        free(spec->gs.large_ptc);
    free(spec);
    // allocate space for velocity
    #ifdef POOL_DEBUG
    INFO("Allocate small_vel and large_vel with size of %u, %u",
        spec->gs.small_num, spec->gs.large_num)
    #endif
    small_vel = (Velocity *)calloc(pool.small_num, sizeof(Velocity));
    large_vel = (Velocity *)calloc(pool.large_num, sizeof(Velocity));
    memset(small_vel, 0, pool.small_num * sizeof(Velocity));
    memset(large_vel, 0, pool.large_num * sizeof(Velocity));
    // now begins the main computation
    // set the adjacent regions that need to be considered
    int row = rank / proc_size, col = rank % proc_size;
    adj_num = 0;
    for (int i = row - cs.horizon; i <= row + cs.horizon; i++) {
        for (int j = col - cs.horizon; j <= col + cs.horizon; j++) {
            if (i < 0 || i >= proc_size || j < 0 || j >= proc_size)
                continue;
            if (i == row && j == col)
                continue;
            adj_ranks[adj_num++] = i * proc_size + j;
        }
    }
    #ifdef POOL_DEBUG
    for (int i = 0; i < adj_num ; i++)
        INFO("adj_ranks[i] = %d", adj_ranks[i]);
    #endif
    // run time_slot steps on the pool
    for (int i = 0; i < cs.time_slot; i++) {
        #ifdef POOL_DEBUG
        printf("Running step %d\n", i);
        INFO("Running step %d", i);
        #endif
        run_step(rank, size);
        #ifdef POOL_DEBUG
        char debugppm[64];
        sprintf(debugppm, "%s-debug-%d-%d.ppm", argv[2], rank, i);
        if (print2ppm(&pool, debugppm)) {
            ERR("print2ppm failed");        
            return -1;
        }
        #endif
    }
    // free the allocated space
    free(small_vel);
    free(large_vel);
    // output results to ppm file
    PPMFile finalbrd;
    if (get_final_result(rank, size, &finalbrd)) {
        fprintf(stderr, "get_final_result failed\n");
        return -1;
    }
    if (rank == 0) {
        char ppmfile[64];
        sprintf(ppmfile, "%s-%d.ppm", argv[2], rank);
        if (ppm2file(&finalbrd, ppmfile)) {
            fprintf(stderr, "print2ppm failed\n");
            return -1;
        }
        free(finalbrd.pixels);
    }
    free(pool.small_ptc);
    free(pool.large_ptc);
    return 0;
}

int __is_perfect_square(int num) {
    int sqrt_rank = (int)round(sqrt((double)num));
    for (int i = 1; i <= sqrt_rank; i++)
        if (i * i == num)
            return i;
    return 0;
}

int main(int argc, char *argv[]) {
	int size, rank;
    MPI_Init(NULL, NULL);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    // ensures the number of processes is a perfect square
    proc_size = __is_perfect_square(size);
    if (!proc_size) {
        fprintf(stderr, "The number of processes should be a perfect square\n");
        return EXIT_FAILURE;
    }

    #ifdef POOL_DEBUG
    time_t rawtime;
    time(&rawtime);
    struct tm *timeinfo = localtime(&rawtime);
    char logfile[64];
    sprintf(logfile, "log/log-%02d%02d%02d%02d%02d-%d.txt",
            timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec, rank);
    if (logging_open(logfile) < 0) {
        fprintf(stderr, "Failed to initialize logging\n");
        return EXIT_FAILURE;
    }
    #endif

    init_mympi();

    int ret = run(rank, size, argc, argv);

    #ifdef POOL_DEBUG
    logging_close();
    #endif

    MPI_Finalize();
    return ret;
}

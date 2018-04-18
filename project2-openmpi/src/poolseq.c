#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <time.h>

#include "spec.h"
#include "pfile.h"
#include "pool.h"

#ifdef POOL_DEBUG
#include "logging.h"
#endif

/**
 * Pool and Velocity describing the current state of the board
 */
Pool *pool;
/**
 * Velocity info of the small and large particles
 * Separately defined since we won't transfer them between processes
 */
Velocity **small_vel, **large_vel;
/**
 * CompSpec for current computation
 */
CompSpec cs;
int **adj_ranks;
int *adj_num;
int proc_size;
int np;

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
int run_step();
void init_pool(const Spec *spec);
int init_params(int argc, char *argv[], Spec *spec);
int get_final_result(PPMFile *ppm);

int get_final_result(PPMFile *ppm) {
    uint64_t rgb_num = pool[0].size * pool[0].size;
        // assume every node has the same size of grid
    ppm->size = pool[0].size * proc_size;
    ppm->pixels = (RGBPixel *)calloc(ppm->size * ppm->size, sizeof(RGBPixel));
    PPMFile temp;
    if (pool2ppm(&pool[0], &temp)) {
        #ifdef POOL_DEBUG
        INFO("pool2ppm for rank %d failed", 0)
        #endif
        return -1;
    }
    memcpy(ppm->pixels, temp.pixels, sizeof(RGBPixel) * rgb_num);
    RGBPixel *p = ppm->pixels;
    for (int i = 1; i < np; i++) {
        p += rgb_num;
        if (pool2ppm(&pool[i], &temp)) {
            #ifdef POOL_DEBUG
            INFO("pool2ppm for rank %d failed", i)
            #endif
            return -1;
        }
        memcpy(p, temp.pixels, sizeof(RGBPixel) * rgb_num);
    }
    return 0;
}

void init_pool(const Spec *spec) {
    pool = (Pool *)calloc(np, sizeof(Pool));
    // initialize pool[0] first
    pool[0].size = spec->gs.size;
    pool[0].small_num = spec->gs.small_num;
    pool[0].small_mass = spec->gs.small_mass;
    pool[0].small_rad = spec->gs.small_rad;
    // initialize small paticles
    int snum = pool[0].small_num;
    pool[0].small_ptc = (Location *)calloc(snum, sizeof(Location));
    srand(time(NULL));
    for (int i = 0, max = pool[0].size; i < snum; i++) {
        pool[0].small_ptc[i].x = ((double)rand() / RAND_MAX) * max;
        pool[0].small_ptc[i].y = ((double)rand() / RAND_MAX) * max;
    }
    // initialize large paticles
    int lnum = pool[0].large_num = spec->gs.large_num;
    pool[0].large_ptc = (Particle *)calloc(lnum, sizeof(Particle));
    memcpy(pool[0].large_ptc, spec->gs.large_ptc, lnum*sizeof(Particle));
    #ifdef POOL_DEBUG
    __debug_print_pool(0, &pool[0]);
    #endif
    for (int i = 1; i < np; i++) {
        memcpy(&pool[i], &pool[0], sizeof(Pool));
        pool[i].small_ptc = (Location *)calloc(snum, sizeof(Location));
        memcpy(pool[i].small_ptc, pool[0].small_ptc, snum*sizeof(Location));
        pool[i].large_ptc = (Particle *)calloc(lnum, sizeof(Particle));
        memcpy(pool[i].large_ptc, pool[0].large_ptc, lnum*sizeof(Particle));
        #ifdef POOL_DEBUG
        __debug_print_pool(i, &pool[i]);
        #endif
    }
}

int init_params(int argc, char *argv[], Spec *spec) {
    // check the number of arguments
    if (read_spec_from_file(argv[1], spec)) {
        fprintf(stderr, "read_spec_from_file failed\n");
        return -1;
    }
    memcpy(&cs, &spec->cs, sizeof(CompSpec));
    #ifdef POOL_DEBUG
    __debug_print_compspec(-1, &cs);
    #endif
    return 0;
}

/**
 * Run one step of the computation
 * Note that currently small_num and large_num in any of the region should be
 * positive. If zero number occurs, the program may crash
 */
int run_step() {
    // allocate space for accelaration here
    #ifdef POOL_DEBUG
    INFO("Allocating space for accelaration");
    #endif
    Accel *small_acc[np], *large_acc[np];
    for (int i = 0; i < np; i++) {
        small_acc[i] = (Accel *)calloc(pool[i].small_num, sizeof(Accel));
        large_acc[i] = (Accel *)calloc(pool[i].large_num, sizeof(Accel));
        memset(small_acc[i], 0, pool[i].small_num * sizeof(Accel));
        memset(large_acc[i], 0, pool[i].large_num * sizeof(Accel));
    }
    for (int rank = 0; rank < np; rank++) {
        // compute impacts of the region itself and wait for communication to finish
        #ifdef POOL_DEBUG
        INFO("Computing the impacts of the region %d itself", rank)
        #endif
        // only compute F(i, j) sice F(j, i) = -F(i, j)
        for (int i = 0; i < pool[rank].small_num; i++) {
            // forces between small particles
            for (int j = i + 1; j < pool[rank].small_num; j++) {
                // compute ptc[j]'s impact on ptc[i]
                Location r_ji = {
                    pool[rank].small_ptc[j].x - pool[rank].small_ptc[i].x,
                    pool[rank].small_ptc[j].y - pool[rank].small_ptc[i].y
                };
                COMPUTE_ACCEL_2PTC_RAW(r_ji, pool[rank].small_mass, pool[rank].small_mass,
                    small_acc[rank][i], small_acc[rank][j])
            }
            // forces between larger particles and small particles
            for (int j = 0; j < pool[rank].large_num; j++) {
                Location r_ji = {
                    pool[rank].large_ptc[j].loc.x - pool[rank].small_ptc[i].x,
                    pool[rank].large_ptc[j].loc.y - pool[rank].small_ptc[i].y
                };
                COMPUTE_ACCEL_2PTC_RAW(r_ji, pool[rank].small_mass, pool[rank].large_ptc[j].mass,
                    small_acc[rank][i], large_acc[rank][j])
            }
        }
        for (int i = 0; i < pool[rank].large_num; i++) {
            for (int j = i + 1; j < pool[rank].large_num; j++) {
                Location r_ji = {
                    pool[rank].large_ptc[j].loc.x - pool[rank].large_ptc[i].loc.x,
                    pool[rank].large_ptc[j].loc.y - pool[rank].large_ptc[i].loc.y
                };
                COMPUTE_ACCEL_2PTC_RAW(r_ji, pool[rank].large_ptc[i].mass, pool[rank].large_ptc[j].mass,
                    large_acc[rank][i], large_acc[rank][j])
            }
        }
        // compute the impacts of adjacent regions
        for (int i = 0; i < adj_num[rank]; i++) {
            // compute the offset of two adjacent regions
            double offset[2] = {
                ((adj_ranks[rank][i] % proc_size) - (rank % proc_size)) * (double)pool[rank].size,
                ((adj_ranks[rank][i] / proc_size) - (rank / proc_size)) * (double)pool[rank].size            
            };
            #ifdef POOL_DEBUG
            INFO("compute region %d's impact, offset = { %.4f, %.4f }",
                adj_ranks[rank][i], offset[0], offset[1])
            #endif
            // compute pool[adj_ranks[rank][i]].large_ptc[j]'s impact
            for (int j = 0; j < pool[adj_ranks[rank][i]].large_num; j++) {
                #ifdef POOL_DEBUG
                INFO("Process %d large_ptc[%d].loc: %.4f, %.4f", adj_ranks[i],
                    j, pool[adj_ranks[rank][i]].large_ptc[j].loc.x,
                    pool[adj_ranks[rank][i]].large_ptc[j].loc.y)
                #endif
                for (int k = 0; k < pool[rank].small_num; k++) {
                    // compute pool[adj_ranks[rank][i]].large_ptc[j]'s impact on pool.small_ptc[k]
                    Location r_kj = {
                        pool[adj_ranks[rank][i]].large_ptc[j].loc.x - pool[rank].small_ptc[k].x + offset[0],
                        pool[adj_ranks[rank][i]].large_ptc[j].loc.y - pool[rank].small_ptc[k].y + offset[1]
                    };
                    COMPUTE_ACCEL_RAW(r_kj, pool[adj_ranks[rank][i]].large_ptc[j].mass, small_acc[rank][k])
                }
                for (int k = 0; k < pool[rank].large_num; k++) {
                    // compute pool[adj_ranks[rank][i]].large_ptc[j]'s impact on pool.large_ptc[k]                
                    Location r_kj = {
                        pool[adj_ranks[rank][i]].large_ptc[j].loc.x - pool[rank].large_ptc[k].loc.x + offset[0],
                        pool[adj_ranks[rank][i]].large_ptc[j].loc.y - pool[rank].large_ptc[k].loc.y + offset[1]
                    };
                    COMPUTE_ACCEL_RAW(r_kj, pool[adj_ranks[rank][i]].large_ptc[j].mass, large_acc[rank][k])
                }
            }
            // compute adj_small_ptc[i][j]'s impact
            for (int j = 0; j < pool[adj_ranks[rank][i]].small_num; j++) {
                #ifdef POOL_DEBUG
                if (j < 10)
                INFO("Process %d small_ptc[%d]: %.4f, %.4f", adj_ranks[i],
                    j, pool[adj_ranks[rank][i]].small_ptc[j].x, pool[adj_ranks[rank][i]].small_ptc[j].x)
                #endif
                for (int k = 0; k < pool[rank].small_num; k++) {
                    // compute adj_small_ptc[i][j]'s impact on pool.small_ptc[k]
                    Location r_kj = {
                        pool[adj_ranks[rank][i]].small_ptc[j].x - pool[rank].small_ptc[k].x + offset[0],
                        pool[adj_ranks[rank][i]].small_ptc[j].y - pool[rank].small_ptc[k].y + offset[1]
                    };
                    COMPUTE_ACCEL_RAW(r_kj, pool[rank].small_mass, small_acc[rank][k])
                }
                for (int k = 0; k < pool[rank].large_num; k++) {
                    // compute adj_small_ptc[i][j]'s impact on pool.large_ptc[k]
                    Location r_kj = {
                        pool[adj_ranks[rank][i]].small_ptc[j].x - pool[rank].large_ptc[k].loc.x + offset[0],
                        pool[adj_ranks[rank][i]].small_ptc[j].y - pool[rank].large_ptc[k].loc.y + offset[1]
                    };
                    COMPUTE_ACCEL_RAW(r_kj, pool[rank].small_mass, large_acc[rank][k])
                }
            }
            #ifdef POOL_DEBUG
            INFO("Apply G constant to Acceleration")
            #endif
            for (int i = 0; i < pool[rank].small_num; i++)
                ACCEL_RAW_TO_FINAL(small_acc[rank][i])
            for (int i = 0; i < pool[rank].large_num; i++)
                ACCEL_RAW_TO_FINAL(large_acc[rank][i])
            #ifdef POOL_DEBUG
            for (int i = 0; i < 10; i++) {
                INFO("small_acc[%d]: %.4f, %.4f", i, small_acc[rank][i].ax, small_acc[rank][i].ay);        
            }
            for (int i = 0; i < pool[rank].large_num; i++) {
                INFO("large_acc[%d]: %.4f, %.4f", i, large_acc[rank][i].ax, large_acc[rank][i].ay);
            }
            #endif
        }
        for (int rank = 0; rank < np; rank++) {
            // now we are going to update the velocity and location
            double t = cs.time_step;
            #ifdef POOL_DEBUG
            INFO("Update velocity using vt = v0 + a * %.4f", t)
            #endif
            for (int i = 0; i < pool[rank].small_num; i++) {
                small_vel[rank][i].vx += small_acc[rank][i].ax * t;
                small_vel[rank][i].vy += small_acc[rank][i].ay * t;
            }
            for (int i = 0; i < pool[rank].large_num; i++) {
                large_vel[rank][i].vx += large_acc[rank][i].ax * t;
                large_vel[rank][i].vy += large_acc[rank][i].ay * t;
            }
            #ifdef POOL_DEBUG
            for (int i = 0; i < 10; i++) {
                INFO("small_vel[%d]: %.4f, %.4f", i, small_vel[rank][i].vx, small_vel[rank][i].vy);        
            }
            for (int i = 0; i < pool[rank].large_num; i++) {
                INFO("large_vel[%d]: %.4f, %.4f", i, large_vel[rank][i].vx, large_vel[rank][i].vy);
            }
            #endif
            #ifdef POOL_DEBUG
            INFO("Update location using x = vt * t + 0.5 * a * t * t")
            #endif
            for (int i = 0; i < pool[rank].small_num; i++) {
                pool[rank].small_ptc[i].x += small_vel[rank][i].vx * t + 0.5 * small_acc[rank][i].ax * t * t;
                pool[rank].small_ptc[i].y += small_vel[rank][i].vy * t + 0.5 * small_acc[rank][i].ay * t * t;
                while (pool[rank].small_ptc[i].x > pool[rank].size)
                    pool[rank].small_ptc[i].x -= pool[rank].size;
                while (pool[rank].small_ptc[i].x < 0)
                    pool[rank].small_ptc[i].x += pool[rank].size;
                while (pool[rank].small_ptc[i].y > pool[rank].size)
                    pool[rank].small_ptc[i].y -= pool[rank].size;
                while (pool[rank].small_ptc[i].y < 0)
                    pool[rank].small_ptc[i].y += pool[rank].size;
                #ifdef POOL_DEBUG
                if (i < 10) {
                    INFO("pool.small_ptc[%d]: %.4f, %.4f", i, pool[rank].small_ptc[i].x, pool[rank].small_ptc[i].y);
                }
                #endif
            }
            for (int i = 0; i < pool[rank].large_num; i++) {
                // TODO: currently if a particle runs out of the region, we make it run into
                // the region in the other direction
                pool[rank].large_ptc[i].loc.x += large_vel[rank][i].vx * t + 0.5 * large_acc[rank][i].ax * t * t;
                pool[rank].large_ptc[i].loc.y += large_vel[rank][i].vy * t + 0.5 * large_acc[rank][i].ay * t * t;
                while (pool[rank].large_ptc[i].loc.x > pool[rank].size)
                    pool[rank].large_ptc[i].loc.x -= pool[rank].size;
                while (pool[rank].large_ptc[i].loc.x < 0)
                    pool[rank].large_ptc[i].loc.x += pool[rank].size;
                while (pool[rank].large_ptc[i].loc.y > pool[rank].size)
                    pool[rank].large_ptc[i].loc.y -= pool[rank].size;
                while (pool[rank].large_ptc[i].loc.y < 0)
                    pool[rank].large_ptc[i].loc.y += pool[rank].size;
                #ifdef POOL_DEBUG
                INFO("pool.large_ptc[i].loc: %.4f, %.4f", pool[rank].large_ptc[i].loc.x, pool[rank].large_ptc[i].loc.y);
                #endif
            }
        }
    }
    // free the allocated acc space
    for (int i = 0; i < np; i++) {
        free(small_acc[i]);
        free(large_acc[i]);
    }
    return 0;
}

int run(int argc, char *argv[]) {
    // initialize the game (CompSpec and Pool)
    Spec *spec = (Spec *)malloc(sizeof(Spec));
    if (init_params(argc, argv, spec)){
        fprintf(stderr, "init_params failed\n");        
        return -1;
    }
    init_pool(spec);
    // free the allocated space
    // only in rank 0 have we initialized spec.gs in spec
    free(spec->gs.large_ptc);
    free(spec);
    // allocate space for velocity
    #ifdef POOL_DEBUG
    INFO("Allocate small_vel and large_vel with size of %u, %u",
        spec->gs.small_num, spec->gs.large_num)
    #endif
    small_vel = (Velocity **)calloc(np, sizeof(Velocity *));
    large_vel = (Velocity **)calloc(np, sizeof(Velocity *));
    for (int i = 0; i < np; i++) {
        small_vel[i] = (Velocity *)calloc(pool[i].small_num, sizeof(Velocity));
        large_vel[i] = (Velocity *)calloc(pool[i].large_num, sizeof(Velocity));
        memset(small_vel[i], 0, pool[i].small_num * sizeof(Velocity));
        memset(large_vel[i], 0, pool[i].large_num * sizeof(Velocity));
    }
    // now begins the main computation
    // set the adjacent regions that need to be considered
    adj_num = (int *)calloc(np, sizeof(int));
    adj_ranks = (int **)calloc(np, sizeof(int *));
    for (int rank = 0; rank < np; rank++) {
        adj_ranks[rank] = (int *)calloc(9, sizeof(int));
        int row = rank / proc_size, col = rank % proc_size;
        adj_num[rank] = 0;
        for (int i = row - cs.horizon; i <= row + cs.horizon; i++) {
            for (int j = col - cs.horizon; j <= col + cs.horizon; j++) {
                if (i < 0 || i >= proc_size || j < 0 || j >= proc_size)
                    continue;
                if (i == row && j == col)
                    continue;
                adj_ranks[rank][adj_num[rank]++] = i * proc_size + j;
            }
        }
    }
    #ifdef POOL_DEBUG
    for (int i = 0; i < np ; i++)
        INFO("adj_num[i] = %d", adj_num[i]);
    #endif
    // run time_slot steps on the pool
    // note that we only count the iteration time
    struct timeval tv1, tv2;    
    gettimeofday(&tv1, NULL);
    for (int i = 0; i < cs.time_slot; i++) {
        #ifdef POOL_DEBUG
        printf("Running step %d\n", i);
        INFO("Running step %d", i);
        #endif
        run_step();
        #ifdef POOL_DEBUG
        char debugppm[64];
        sprintf(debugppm, "%s-debug-%d.ppm", argv[2], i);
        PPMFile finalbrd;
        if (get_final_result(&finalbrd)) {
            fprintf(stderr, "get_final_result failed\n");
            return -1;
        }
        if (ppm2file(&finalbrd, debugppm)) {
            fprintf(stderr, "print2ppm failed\n");
            return -1;
        }
        free(finalbrd.pixels);
        #endif
    }
    gettimeofday(&tv2, NULL);
    double excu_time = (double)(tv2.tv_usec - tv1.tv_usec) / 1000000 +
        (double)(tv2.tv_sec - tv1.tv_sec);
    printf ("%.4f\n", excu_time);
    // free the allocated space
    for (int i = 0; i < np; i++) {
        free(small_vel[i]);
        free(large_vel[i]);
        free(adj_ranks[i]);
    }
    free(small_vel);
    free(large_vel);
    free(adj_num);
    free(adj_ranks);
    // output results to ppm file
    PPMFile finalbrd;
    if (get_final_result(&finalbrd)) {
        fprintf(stderr, "get_final_result failed\n");
        return -1;
    }
    char ppmfile[64];
    sprintf(ppmfile, "%s.ppm", argv[2]);
    if (ppm2file(&finalbrd, ppmfile)) {
        fprintf(stderr, "print2ppm failed\n");
        return -1;
    }
    free(finalbrd.pixels);
    for (int i = 0; i < np; i++) {
        free(pool[i].small_ptc);
        free(pool[i].large_ptc);
    }
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
	if (argc != 4) {
        fprintf(stderr, "Unknown arguments...\n\n");
        fprintf(stderr, "Usage: pool <specfile> <ppmfile> <np>\n");
        return -1;
    }

    // ensures the number of processes is a perfect square
    np = atoi(argv[3]);
    proc_size = __is_perfect_square(np);
    if (!proc_size) {
        fprintf(stderr, "The number of processes should be a perfect square\n");
        return EXIT_FAILURE;
    }

    #ifdef POOL_DEBUG
    time_t rawtime;
    time(&rawtime);
    struct tm *timeinfo = localtime(&rawtime);
    char logfile[64];
    sprintf(logfile, "log/log-%02d%02d%02d%02d%02d.txt",
            timeinfo->tm_mon + 1, timeinfo->tm_mday, timeinfo->tm_hour,
            timeinfo->tm_min, timeinfo->tm_sec);
    if (logging_open(logfile) < 0) {
        fprintf(stderr, "Failed to initialize logging\n");
        return EXIT_FAILURE;
    }
    #endif

    int ret = run(argc, argv);

    #ifdef POOL_DEBUG
    logging_close();
    #endif

    return ret;
}

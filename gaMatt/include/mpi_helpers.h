/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#ifndef MPI_HELPERS_H
#define MPI_HELPERS_H

#define MASTER 0

#define MPI_TAG_PAIRSCORES 1
#define MPI_TAG_TASK 2
#define MPI_TAG_TASK_REQ 3
#define MPI_TAG_TASK_INT 4

#define MATT_LOG_LVL "MATT_LOG_LVL"

#include <pthread.h>

static int LOG_LVL;

/* mpi_init_log reads LOG_LVL value from MATT_LOG_LVL env variable
   if variable is undefined than default value will be used.

   LOL_LVL 0 -- log only timings for pairwise alingments and iterative part of the algorithm
   LOG_LVL 1 -- all above plus the information about tasks sent
   LOG_LVL 2 -- all above plus timings for each alignment and barrier wait time (default)
*/
void mpi_init_log();
void mpi_log(const int loglvl, const int master_only, const char *fmt, ...);

typedef struct mpi_task_t {
	/* alignments IDs al2 is ingored in iterative phase */
	int al1;
	int al2;
	int ga;
	void *ctx;
} mpi_task_t;


mpi_task_t mpi_create_task(const int al1, const int al2);
mpi_task_t mpi_create_ga_task(void *ctx, const int idx);


#endif

/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#include <assert.h>
#include <errno.h>
#include <mpi.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "MultipleAlignment.h"

#include "mpi_helpers.h"

static int LOG_LVL;

/* mpi_init_log reads LOG_LVL value from MATT_LOG_LVL env variable
   if variable is undefined than default value will be used.

   LOL_LVL 0 -- log only timings for pairwise alingments and iterative part of the algorithm
   LOG_LVL 1 -- all above plus the information about tasks sent
   LOG_LVL 2 -- all above plus timings for each alignment and barrier wait time (default)
*/


mpi_task_t mpi_create_task(const int al1, const int al2)
{
	mpi_task_t t = {
	    .al1 = al1,
	    .al2 = al2,
	    .ga = 0,
		.ctx = NULL,
	};
	return t;
}

mpi_task_t mpi_create_ga_task(void *ctx, const int idx)
{
	mpi_task_t t = {
	    .al1 = idx,
	    .al2 = -1,
	    .ga = 1,
		.ctx = ctx,
	};
	return t;
}

void mpi_init_log()
{
	char *v = getenv(MATT_LOG_LVL);
	if (v) {
		errno = 0;
		LOG_LVL = strtol(v, NULL, 10);
		if (errno == ERANGE || errno == EINVAL || LOG_LVL < 0 || LOG_LVL > 2) {
			LOG_LVL = 2;
		}
	}
	else {
		LOG_LVL = 2;
	}

	mpi_log(0, 0, "loglevel set to %d", LOG_LVL);
}


void mpi_log(const int loglvl, const int master_only, const char *fmt, ...)
{
	if (loglvl > LOG_LVL) {
		return;
	}

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (master_only && rank != MASTER) {
		return;
	}

	char buff[100];
	time_t now = time(0);
	strftime(buff, 100, "%Y-%m-%d %H:%M:%S", localtime(&now));
	if (rank == MASTER) {
		printf("%s master: ", buff);
	}
	else {
		printf("%s slave %d: ", buff, rank);
	}

	va_list arg;
	va_start(arg, fmt);
	vprintf(fmt, arg);
	va_end(arg);
	printf("\n");
}

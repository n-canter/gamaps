/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#ifndef MPI_POOL_H
#define MPI_POOL_H

#include <pthread.h>
#include "mpi_helpers.h"

#define SLEEP_TIME (100 * 1000)

/* execute_cb_t is a callback for taskpool_t to execute task, appends return value to the associated processed_t vector */
typedef void *(*execute_cb_t)(mpi_task_t t);
typedef int *(*update_cb_t)(int len, const int rank);

struct bucket_t {
	int len;
	int cap;
	int start;
	int rank;
	mpi_task_t *tasks;
};

typedef struct snapshot_t {
	int n;
	struct bucket_t *procs;
} snapshot_t;

/* taskpool_t is a pool of workers running on nworkers threads plus one thread for fetching the tasks */
typedef struct taskpool_t {
	int len;
	int nworkers;

	struct bucket_t *queue;

	pthread_mutex_t *lock;
	pthread_cond_t *not_empty;

	update_cb_t update;
	execute_cb_t execute;
} taskpool_t;


/* procssed_t is a vector with a thread-safe append operation */
typedef struct processed_t {
	int len;
	int cap;

	pthread_mutex_t *lock;
	void **data;
} processed_t;

typedef struct taskpool_ctx_t {
	int rank;
	snapshot_t *snap;
	processed_t *processed;
	taskpool_t *pool;
} taskpool_ctx_t;

snapshot_t *snapshot_init(const int n_procs);
void snapshot_fill_pairwise(snapshot_t *s, const int n_als);
void snapshot_fill_iterative(snapshot_t *s, const int n_als);
void snapshot_fill_ga(void *ctx, snapshot_t *s, const int n_trees);
void snapshot_destroy(snapshot_t *s);


void taskpool_destroy(taskpool_t *p);
void taskpool_run(int nworkers, taskpool_ctx_t ctx, update_cb_t update, execute_cb_t execute);

processed_t *processed_init();
void processed_destroy(processed_t *p);
void processed_append(processed_t *p, void *x);

int *mpi_task_updater(int len, const int nprocs);

#endif

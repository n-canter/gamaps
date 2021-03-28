/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#include <assert.h>
#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>
#include <mpi.h>

#include "mpi_pool.h"

const mpi_task_t TERMINATE_TASK = {
    .al1 = -2,
    .al2 = -2,
};

static void bucket_init(struct bucket_t *b, const int rank)
{
	b->len = 0;
	b->cap = 1;
	b->start = 0;
	b->rank = rank;
	b->tasks = calloc(1, sizeof(mpi_task_t));
	assert(b->tasks);
}

static int bucket_count(struct bucket_t *b)
{
	return b->len - b->start;
}

static void bucket_reset(struct bucket_t *b, const int rank)
{
	b->len = 0;
	b->start = 0;
	b->rank = rank;
}

static void bucket_push(struct bucket_t *b, mpi_task_t task)
{
	if (b->len >= b->cap) {
		b->cap *= 2;
		b->tasks = realloc(b->tasks, b->cap * sizeof(mpi_task_t));
		assert(b->tasks);
	}
	b->tasks[b->len++] = task;
}

static mpi_task_t bucket_pop(struct bucket_t *b)
{
	assert(b->len - b->start > 0);
	return b->tasks[--b->len];
}

static mpi_task_t bucket_pop_front(struct bucket_t *b)
{
	assert(b->len - b->start > 0);
	return b->tasks[b->start++];
}

static void bucket_leave_last(struct bucket_t *b, const int n)
{
	assert(b->len - n >= b->start);
	b->start = b->len - n;
}

typedef struct diff_t {
	int empty;
	int deleted;
	int appended;
	mpi_task_t *tasks;
} diff_t;

void diff_destroy(diff_t *diff)
{
	if (diff->tasks) {
		free(diff->tasks);
	}
}

static void snapshot_print(snapshot_t *s)
{
	printf("*************\n");
	for (int i = 0; i < s->n; i++) {
		printf("%d: ", s->procs[i].rank);
		for (int j = s->procs[i].start; j < s->procs[i].len; j++) {
			printf("(%d %d) ", s->procs[i].tasks[j].al1, s->procs[i].tasks[j].al2);
		}
		printf("\n");
	}
	printf("*************\n");
}


static void snapshot_reset(snapshot_t *s)
{
	for (int i = 0; i < s->n; i++) {
		bucket_reset(s->procs + i, i);
	}
}

snapshot_t *snapshot_init(const int n_procs)
{
	snapshot_t *s = calloc(1, sizeof(snapshot_t));
	assert(s);

	s->procs = calloc(n_procs, sizeof(struct bucket_t));
	assert(s->procs);
	for (int i = 0; i < n_procs; i++) {
		bucket_init(s->procs + i, i);
	}

	s->n = n_procs;
	return s;
}

void snapshot_fill_pairwise(snapshot_t *s, const int n_als)
{
	snapshot_reset(s);
	int ind = 0;
	for (int i = 0; i < n_als; i++) {
		for (int j = 0; j < i; j++) {
			mpi_task_t task = mpi_create_task(i, j);
			bucket_push(s->procs + ind % s->n, task);
			ind++;
		}
	}
}

void snapshot_fill_iterative(snapshot_t *s, const int n_als)
{
	snapshot_reset(s);
	for (int i = 0; i < n_als; i++) {
		mpi_task_t task = mpi_create_task(i, -1);
		bucket_push(s->procs + i % s->n, task);
	}
}

void snapshot_fill_ga(void *ctx, snapshot_t *s, const int n_trees)
{
	snapshot_reset(s);
	for (int i = 0; i < n_trees; i++) {
		mpi_task_t task = mpi_create_ga_task(ctx, i);
		bucket_push(s->procs + i % s->n, task);
	}
}

void snapshot_destroy(snapshot_t *s)
{
	for (int i = 0; i < s->n; i++) {
		free(s->procs[i].tasks);
	}
	free(s->procs);
	free(s);
}

static int sum(const int *procs, const int n)
{
	int ret = 0;
	for (int i = 0; i < n; i++) {
		ret += procs[i];
	}
	return ret;
}

static int threshold(const int sum, const int n, const int rank)
{
	return sum / n + (n - rank <= sum % n);
}

static int bucket_comp_rank(const void *_a, const void *_b)
{
	struct bucket_t a = *((struct bucket_t *)_a);
	struct bucket_t b = *((struct bucket_t *)_b);
	return a.rank - b.rank;
}

static int bucket_comp(const void *_a, const void *_b)
{
	struct bucket_t a = *((struct bucket_t *)_a);
	struct bucket_t b = *((struct bucket_t *)_b);
	int count_a = bucket_count(&a);
	int count_b = bucket_count(&b);
	if (count_a == count_b) {
		return bucket_comp_rank(_a, _b);
	}
	return count_a - count_b;
}


diff_t snapshot_update(snapshot_t *s, const int *procs, const int rank)
{
	diff_t diff = {
	    .empty = 0,
	    .deleted = 0,
	    .appended = 0,
	    .tasks = NULL,
	};

	int sm = sum(procs, s->n);
	qsort(s->procs, s->n, sizeof(struct bucket_t), bucket_comp_rank);
	diff.empty = sm == 0;

	for (int i = 0; i < s->n; i++) {
		bucket_leave_last(s->procs + i, procs[i]);
	}

	qsort(s->procs, s->n, sizeof(struct bucket_t), bucket_comp);

	int i = 0;
	int j = s->n - 1;
	for (;;) {
		while (i < j && bucket_count(s->procs + i) >= threshold(sm, s->n, i)) {
			i++;
		}
		while (i < j && bucket_count(s->procs + j) <= threshold(sm, s->n, j)) {
			j--;
		}
		if (i >= j) {
			break;
		}
		mpi_task_t t = bucket_pop(s->procs + j);
		bucket_push(s->procs + i, t);
		if (s->procs[j].rank == rank) {
			diff.deleted++;
		}
		else if (s->procs[i].rank == rank) {
			if (diff.tasks == NULL) {
				diff.tasks = calloc(sm / s->n + 1, sizeof(mpi_task_t));
			}
			diff.tasks[diff.appended++] = t;
		}
	}

	return diff;
}

static void *executor(void *context)
{
	taskpool_ctx_t *ctx = context;
	taskpool_t *pool = ctx->pool;
	struct bucket_t *queue = pool->queue;

	processed_t *processed = ctx->processed;
	for (;;) {
		pthread_mutex_lock(pool->lock);
		while (bucket_count(queue) == 0) {
			pthread_cond_wait(pool->not_empty, pool->lock);
		}
		mpi_task_t t = bucket_pop_front(queue);
		pthread_mutex_unlock(pool->lock);
		if (t.al1 == TERMINATE_TASK.al1 && t.al2 == TERMINATE_TASK.al2) {
			break;
		}
		void *res = pool->execute(t);
		if (processed != NULL) {
			processed_append(processed, res);
		}
	}
	return NULL;
}

static void *updater(void *context)
{
	taskpool_ctx_t *ctx = context;
	taskpool_t *pool = ctx->pool;
	snapshot_t *snap = ctx->snap;
	struct bucket_t *queue = pool->queue;
	int rank = ctx->rank;

	// snap->procs[rank] should be untouched yet
	struct bucket_t b = snap->procs[rank];
	pthread_mutex_lock(pool->lock);
	for (int i = 0; i < b.len; i++) {
		bucket_push(queue, b.tasks[i]);
	}
	pthread_cond_broadcast(pool->not_empty);
	pthread_mutex_unlock(pool->lock);


	struct timespec start = {0};
	struct timespec end = {0};
	for (int terminate = 0; !terminate;) {
		pthread_mutex_lock(pool->lock);
		int len = bucket_count(queue);
		pthread_mutex_unlock(pool->lock);

		clock_gettime(CLOCK_MONOTONIC_RAW, &end);
		uint64_t elapsed = (end.tv_sec - start.tv_sec) * 1000000 + (end.tv_nsec - start.tv_nsec) / 1000;
		if (rank == MASTER && elapsed < SLEEP_TIME) {
			usleep(SLEEP_TIME - elapsed);
		}
		clock_gettime(CLOCK_MONOTONIC_RAW, &start);
		int *procs = pool->update(len, snap->n);

		diff_t diff = snapshot_update(snap, procs, rank);

		pthread_mutex_lock(pool->lock);
		for (int i = 0; i < diff.deleted; i++) {
			if (bucket_count(queue) > 0) {
				bucket_pop(queue);
			}
		}
		for (int i = 0; i < diff.appended; i++) {
			bucket_push(queue, diff.tasks[i]);
		}
		if (diff.empty) {
			for (int i = 0; i < pool->nworkers; i++) {
				bucket_push(queue, TERMINATE_TASK);
			}
			terminate = 1;
		}
		pthread_cond_broadcast(pool->not_empty);

		pthread_mutex_unlock(pool->lock);

		diff_destroy(&diff);
		free(procs);
	}
	return NULL;
}

static taskpool_t *taskpool_init(const int rank, const int nworkers, const update_cb_t update, const execute_cb_t execute)
{
	taskpool_t *p = calloc(1, sizeof(taskpool_t));
	assert(p);
	p->len = 0;
	p->nworkers = nworkers;

	p->queue = calloc(1, sizeof(struct bucket_t));
	assert(p->queue);
	bucket_init(p->queue, rank);

	p->lock = malloc(sizeof(pthread_mutex_t));
	assert(p->lock);
	pthread_mutex_init(p->lock, NULL);

	p->not_empty = malloc(sizeof(pthread_cond_t));
	assert(p->not_empty);
	pthread_cond_init(p->not_empty, NULL);

	p->execute = execute;
	p->update = update;
	return p;
}

void taskpool_destroy(taskpool_t *p)
{
	pthread_mutex_destroy(p->lock);
	free(p->lock);
	pthread_cond_destroy(p->not_empty);
	free(p->not_empty);
	free(p);
}

void taskpool_run(int nworkers, taskpool_ctx_t ctx, update_cb_t update, execute_cb_t execute)
{
	taskpool_t *p = taskpool_init(ctx.rank, nworkers, update, execute);
	ctx.pool = p;

	pthread_t executors[nworkers];
	pthread_t update_pthread;

	for (int i = 0; i < nworkers; i++) {
		pthread_create(executors + i, NULL, executor, &ctx);
	}
	pthread_create(&update_pthread, NULL, updater, &ctx);

	for (int i = 0; i < nworkers; i++) {
		pthread_join(executors[i], NULL);
	}
	pthread_join(update_pthread, NULL);

	taskpool_destroy(p);
}

processed_t *processed_init()
{
	processed_t *p = calloc(1, sizeof(processed_t));
	assert(p);
	p->len = 0;
	p->cap = 1;
	p->data = calloc(1, sizeof(void *));
	assert(p->data);
	p->lock = malloc(sizeof(pthread_mutex_t));
	assert(p->lock);
	pthread_mutex_init(p->lock, NULL);
	return p;
}

void processed_destroy(processed_t *p)
{
	free(p->data);
	pthread_mutex_destroy(p->lock);
	free(p->lock);
	free(p);
}

void processed_append(processed_t *p, void *x)
{
	pthread_mutex_lock(p->lock);
	if (p->len >= p->cap) {
		p->cap *= 2;
		p->data = realloc(p->data, p->cap * sizeof(void *));
		assert(p->data);
	}
	p->data[p->len++] = x;
	pthread_mutex_unlock(p->lock);
}

int *mpi_task_updater(int len, const int nprocs)
{
	int *buf = calloc(nprocs, sizeof(int));
	assert(buf);
	MPI_Allgather(&len, 1, MPI_INT, buf, 1, MPI_INT, MPI_COMM_WORLD);
	return buf;
}


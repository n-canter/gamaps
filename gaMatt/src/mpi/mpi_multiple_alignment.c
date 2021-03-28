/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#include <assert.h>
#include <math.h>
#include <mpi.h>
#include <omp.h>
#include <time.h>

#include "MultipleAlignment.h"

#include "ga.h"
#include "mpi_helpers.h"
#include "mpi_pairscores.h"
#include "mpi_pool.h"
#include "mpi_serializer.h"
#include "mt19937ar.h"

static MultipleAlignment **mpi_init_als(PDBChain **chains, const int n)
{
	MultipleAlignment **alignments = (MultipleAlignment **)malloc(n * sizeof(MultipleAlignment *));
	assert(alignments);
	for (int i = 0; i < n; i++) {
		alignments[i] = CreateSingleStrandAlignment(Distill(chains[i], i));
	}
	return alignments;
}

static void mpi_pairwise(MultipleAlignment **als, snapshot_t *snap, pairscores_t **scores, processed_t *aligned, const int nthreads, const int rank)
{
	taskpool_ctx_t ctx = {
	    .processed = aligned,
	    .snap = snap,
	    .rank = rank,
	};
	pthread_mutex_t *lock = malloc(sizeof(pthread_mutex_t));
	assert(lock);
	pthread_mutex_init(lock, NULL);

	void *align(mpi_task_t t)
	{
		double t1 = MPI_Wtime();
		MultipleAlignment *al1 = als[t.al1];
		MultipleAlignment *al2 = als[t.al2];
		MultipleAlignment *ret = AlignAlignments(al1, al2, 0, 0, 0);
		pthread_mutex_lock(lock);
		pairscores_add(scores, t.al1, t.al2, ret->score);
		pthread_mutex_unlock(lock);
		double t2 = MPI_Wtime();
		mpi_log(2, 0, "pairwise, aligned %d %d in %1.5lf seconds", t.al1, t.al2, t2 - t1);
		return ret;
	}

	taskpool_run(nthreads, ctx, mpi_task_updater, align);

	pthread_mutex_destroy(lock);
	free(lock);
}

static void mpi_iterative(MultipleAlignment **als, processed_t *aligned, snapshot_t *snap, MultipleAlignment *ma, pairscores_t **scores, const int iteration, const int nthreads, const int rank)
{
	taskpool_ctx_t ctx = {
	    .processed = aligned,
	    .snap = snap,
	    .rank = rank,
	};

	void *align(mpi_task_t t)
	{
		double t1 = MPI_Wtime();
		MultipleAlignment *al = als[t.al1];
		double bestScore = 0;
		int best1 = 0, best2 = 0;
		for (int j = 0; j < al->numChains; j++) {
			for (int k = 0; k < ma->numChains; k++) {
				pairscores_t *s = pairscores_find(scores, ma->chains[k]->id, al->chains[j]->id);
				if (s->score > bestScore) {
					bestScore = s->score;
					best1 = j;
					best2 = k;
				}
			}
		}
		MultipleAlignment *ret = AlignAlignments(al, ma, best1, best2, 0);
		double t2 = MPI_Wtime();
		mpi_log(2, 0, "iteration %d, aligned %d in %1.5lf seconds", iteration, t.al1, t2 - t1);
		return ret;
	}
	taskpool_run(nthreads, ctx, mpi_task_updater, align);
}

static MultipleAlignment *mpi_bcast_best_ma(mpi_serializer_pool_t *pool, processed_t *p, int rank)
{
	struct best_t {
		double score;
		int rank;
	};

	struct best_t best = {
	    .score = -INFINITY,
	    .rank = rank,
	};
	int index = 0;
	for (int i = 0; i < p->len; i++) {
		MultipleAlignment *al = p->data[i];
		if (al->score > best.score || i == 0) {
			best.score = al->score;
			index = i;
		}
	}
	if (p->len > 0) {
		mpi_serializer_pool_reset(pool);
		mpi_ma_serialize(pool, p->data[index]);
	}
	struct best_t out;
	MPI_Allreduce(&best, &out, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
	size_t size = 0;
	MultipleAlignment *ret = NULL;
	if (rank == out.rank) {
		size = pool->offset;
		ret = p->data[index];
		p->data[index] = p->data[p->len - 1];
		p->len--;
	}
	MPI_Bcast(&size, 1, MPI_INT, out.rank, MPI_COMM_WORLD);
	if (rank != out.rank) {
		mpi_serializer_pool_ensure(pool, size);
	}
	MPI_Bcast(pool->ptr, size, MPI_CHAR, out.rank, MPI_COMM_WORLD);
	if (rank != out.rank) {
		ret = mpi_ma_deserialize(pool->ptr);
	}

	return ret;
}

static int cleanup(MultipleAlignment **als, int n, const int index1, const int index2)
{
	int deleted = 0;
	for (int i = n - 1; i >= 0; i--) {
		for (int j = 0; j < als[i]->numChains; j++) {
			if (als[i]->chains[j]->id != index1 && als[i]->chains[j]->id != index2) {
				continue;
			}

			CleanupAlignment(als[i]);
			n--;
			als[i] = als[n];
			deleted++;
			break;
		}
	}

	return deleted;
}

static void mpi_cleanup_best(MultipleAlignment **als, int *n_als, processed_t *p, MultipleAlignment *ma)
{
	int index1 = ma->chains[0]->id;
	int index2 = ma->chains[ma->numChains - 1]->id;

	*n_als -= cleanup(als, *n_als, index1, index2);

	p->len -= cleanup((MultipleAlignment **)p->data, p->len, index1, index2);
}

static void mpi_barrier_wrapper(const int iteration)
{
	double t1 = MPI_Wtime();
	MPI_Barrier(MPI_COMM_WORLD);
	double t2 = MPI_Wtime();
	if (iteration == -1) {
		mpi_log(2, 0, "iteration pairwise, waiting on barrier %1.5lf seconds", t2 - t1);
		return;
	}
	mpi_log(2, 0, "iteration %d, waiting on barrier %1.5lf seconds", iteration, t2 - t1);
}

int random_int(const int n)
{
	return genrand_int32() % n;
}

double random_double()
{
	return genrand_real2();
}

struct context_t {
	MultipleAlignment **als;
	int n_als;
	pairscores_t **scores;
};

void *to_struct(void *context, const int id)
{
	struct context_t *ctx = context;
	return ctx->als[id];
}

int check_alignment(const MultipleAlignment *ma)
{
	if (ma == NULL) {
		return -1;
	}
	if (ma->chains == NULL) {
		return -1;
	}
	if (ma->blocks == NULL) {
		return -1;
	}
	if (ma->residues == NULL) {
		return -1;
	}

	return 0;
}

void *align_pair(void *context, void *_a, void *_b)
{
	if (_a == NULL || _b == NULL) {
		return NULL;
	}
	struct context_t *ctx = context;
	MultipleAlignment *a = _a;
	MultipleAlignment *b = _b;
	if (check_alignment(a) != 0 || check_alignment(b) != 0) {
		printf("check failed\n");
		return NULL;
	}

	double bestScore = 0;
	int best1 = 0, best2 = 0;
	for (int j = 0; j < b->numChains; j++) {
		for (int k = 0; k < a->numChains; k++) {
			pairscores_t *s = pairscores_find(ctx->scores, a->chains[k]->id, b->chains[j]->id);
			if (s->score > bestScore) {
				bestScore = s->score;
				best1 = j;
				best2 = k;
			}
		}
	}

	MultipleAlignment *sc = AlignAlignments(b, a, best1, best2, 0);
	return sc;
}

static double tm_score_d0(int l)
{
	if (l <= 21) {
		return 0.5;
	}
	return 1.24 * pow(l - 15, 1.0 / 3.0) - 1.8;
}

double tm_score(void *context, void *_ma)
{
	if (_ma == NULL) {
		return 0.0;
	}
	MultipleAlignment *ma = _ma;
	if (ma->numChains <= 1)
		return 0;

	double tm_score = 0.0;
	Vector v;

	double l_target = 0.0;
	for (int i = 0; i < ma->numChains; i++) {
		if (i == 0 || ma->chains[i]->length < l_target) {
			l_target = ma->chains[i]->length;
		}
	}
	double d0 = tm_score_d0(l_target);
	d0 *= d0;
	for (int i = 0; i < ma->numChains - 1; i++) {
		for (int j = i + 1; j < ma->numChains; j++) {
			double s = 0.0;
			for (int k = 0; k < ma->numResidues; k++) {
				if (!ma->residues[i].res[k].hasCoords || !ma->residues[j].res[k].hasCoords) {
					continue;
				}
				double di = lengthSquaredVect(subVect(&v, &ma->residues[i].res[k].coords, &ma->residues[j].res[k].coords));
				s += 1 / (1 + di / d0);
			}
			tm_score += s / l_target;
		}
	}

	int n = ma->numChains;
	n = n * (n - 1);
	return 2 * tm_score / n;
}


static double matt_score(void *context, void *ma)
{
	if (ma == NULL) {
		return 0.0;
	}
	MultipleAlignment *ret = ma;
	return ret->score;
}

void cleanup_alignment(void *context, void *ma)
{
	if (ma == NULL) {
		return;
	}
	CleanupAlignment(ma);
}

static MultipleAlignment *mpi_align_ga(MultipleAlignment **als, const int n_als, pairscores_t **scores, struct ga_cfg_t cfg, const int rank, const int nprocs, const int nthreads)
{
	int seed = time(NULL);
	MPI_Bcast(&seed, 1, MPI_INT, MASTER, MPI_COMM_WORLD);

	init_genrand(seed);
	printf("seed = %d\n", seed);
	struct context_t ctx = {
	    .als = als,
	    .n_als = n_als,
	    .scores = scores,
	};

	cfg.ctx = &ctx;
	cfg.align = align_pair;
	cfg.to_struct = to_struct;
	cfg.rnd_d = random_double;
	cfg.rnd_i = random_int;
	cfg.cleanup = cleanup_alignment;

	cfg.score = matt_score;
	if (cfg.fitness_func == SCORE_TM) {
		cfg.score = tm_score;
	}

	struct ga_t *ga = ga_init(cfg);
	int rc = ga_run(ga);
	assert(rc == 0);

	MultipleAlignment *res = ga_best(ga);

	ga_destroy(ga);
	return res;
}

MultipleAlignment *mpi_align(PDBChain **chains, int n_als, int nthreads, struct ga_cfg_t ga_cfg)
{
	double start_time, end_time;
	start_time = MPI_Wtime();

	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	mpi_serializer_pool_t *pool = mpi_serializer_pool_create();

	MultipleAlignment **alignments = mpi_init_als(chains, n_als);
	processed_t *aligned = processed_init();
	pairscores_t *scores = NULL;

	double t1, t2;
	t1 = MPI_Wtime();
	snapshot_t *snap = snapshot_init(nprocs);

	snapshot_fill_pairwise(snap, n_als);

	mpi_pairwise(alignments, snap, &scores, aligned, nthreads, rank);
	mpi_pairscores_gather(&scores, rank, nprocs); /* sync */

	t2 = MPI_Wtime();
	mpi_log(0, 1, "pairwise time: %1.5lf seconds", t2 - t1);

	MultipleAlignment *ret;
	if (ga_cfg.enabled == 1) {
		ret = mpi_align_ga(alignments, n_als, &scores, ga_cfg, rank, nprocs, nthreads);
	}
	else {
		for (int i = 0; n_als > 1; i++) {
			t1 = MPI_Wtime();
			MultipleAlignment *ma = mpi_bcast_best_ma(pool, aligned, rank); /* sync */
			mpi_cleanup_best(alignments, &n_als, aligned, ma);
			snapshot_fill_iterative(snap, n_als);
			mpi_iterative(alignments, aligned, snap, ma, &scores, i, nthreads, rank);

			mpi_barrier_wrapper(i);
			alignments[n_als++] = ma;
			t2 = MPI_Wtime();
			mpi_log(0, 1, "iteration %d time: %1.5lf seconds", i, t2 - t1);
			ret = alignments[0];
		}
	}

	pairscores_free(&scores);
	processed_destroy(aligned);
	free(alignments);
	mpi_serializer_pool_destroy(pool);
	snapshot_destroy(snap);

	end_time = MPI_Wtime();
	mpi_log(0, 1, "alignment done in %1.5lf seconds (%lf)", end_time - start_time, ret->score);
	if (rank == MASTER) {
		return ret;
	}
	return NULL;
}

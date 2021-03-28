/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#include <assert.h>
#include <mpi.h>
#include <stdlib.h>

#include "uthash.h"

#include "mpi_helpers.h"
#include "mpi_pairscores.h"

#define SWAP(x, y)                  \
	do {                        \
		typeof(x) SWAP = x; \
		x = y;              \
		y = SWAP;           \
	} while (0)


void pairscores_add(pairscores_t **h, int i, int j, const double score)
{
	if (i > j) {
		SWAP(i, j);
	}

	pairscores_t el = {
	    .i = i,
	    .j = j,
	    .score = score,
	};

	pairscores_t *s;
	HASH_FIND(hh, *h, &el.i, 2 * sizeof(int), s);
	if (s) {
		s->score = score;
		return;
	}
	s = calloc(1, sizeof(pairscores_t));
	assert(s);
	*s = el;
	HASH_ADD(hh, *h, i, 2 * sizeof(int), s);
}

pairscores_t *pairscores_find(pairscores_t **h, int i, int j)
{
	if (i > j) {
		SWAP(i, j);
	}
	pairscores_t *ret;
	pairscores_t el = {
	    .i = i,
	    .j = j,
	};
	HASH_FIND(hh, *h, &el.i, 2 * sizeof(int), ret);
	return ret;
}

void pairscores_free(pairscores_t **h)
{
	pairscores_t *el, *tmp = NULL;
	HASH_ITER(hh, *h, el, tmp)
	{
		HASH_DEL(*h, el);
		free(el);
	}
	*h = NULL;
}

void *pairscores_serialize(pairscores_t **h)
{
	int count = HASH_COUNT(*h);
	char *ret = calloc(count, PAIRSCORE_SIZE);
	assert(ret);
	int i = 0;
	pairscores_t *el, *tmp = NULL;
	HASH_ITER(hh, *h, el, tmp)
	{
		memcpy(ret + i * PAIRSCORE_SIZE, el, PAIRSCORE_SIZE);
		i++;
	}
	return ret;
}

void pairscores_deserialize(pairscores_t **h, char *buf, int count)
{
	char *off = buf;
	for (int k = 0; k < count; k++) {
		int i = *((int *)off);
		int j = *((int *)(off + sizeof(int)));
		double score = *((double *)(off + 2 * sizeof(int)));
		pairscores_add(h, i, j, score);
		off += PAIRSCORE_SIZE;
	}
}

int mpi_init_pairscores_type(MPI_Datatype *datatype)
{
	int blocks[2] = {2, 1};
	MPI_Datatype types[2] = {MPI_INT, MPI_DOUBLE};
	struct {
		int a, b;
		double c;
	} in;
	MPI_Aint displacements[2];
	MPI_Get_address(&in, displacements);
	MPI_Get_address(&in.c, displacements + 1);
	displacements[1] -= displacements[0];
	displacements[0] = 0;

	MPI_Type_create_struct(2, blocks, displacements, types, datatype);
	return MPI_Type_commit(datatype);
}

int mpi_free_pairscores_type(MPI_Datatype *datatype)
{
	return MPI_Type_free(datatype);
}

void mpi_pairscores_gather(pairscores_t **p, int rank, int nprocs)
{
	MPI_Datatype datatype;
	mpi_init_pairscores_type(&datatype);
	void *buf;
	if (rank == MASTER) {
		for (int sender = 0; sender < nprocs; sender++) {
			if (sender == MASTER)
				continue;
			MPI_Status status;
			MPI_Probe(sender, MPI_TAG_PAIRSCORES, MPI_COMM_WORLD, &status);
			int count;
			MPI_Get_count(&status, datatype, &count);
			buf = malloc(count * PAIRSCORE_SIZE);
			assert(buf);
			MPI_Recv(buf, count, datatype, sender, MPI_TAG_PAIRSCORES, MPI_COMM_WORLD, &status);
			pairscores_deserialize(p, buf, count);
			free(buf);
		}
	} else {
		buf = pairscores_serialize(p);
		int count = HASH_COUNT(*p);
		MPI_Send(buf, count, datatype, MASTER, MPI_TAG_PAIRSCORES, MPI_COMM_WORLD);
		free(buf);
	}

	int count;
	if (rank == MASTER) {
		count = HASH_COUNT(*p);
		buf = pairscores_serialize(p);
	}
	MPI_Bcast(&count, 1, MPI_INT, MASTER, MPI_COMM_WORLD);
	if (rank != MASTER) {
		buf = calloc(count, PAIRSCORE_SIZE);
		assert(buf);
		pairscores_free(p);
	}

	MPI_Bcast(buf, count, datatype, MASTER, MPI_COMM_WORLD);
	if (rank != MASTER) {
		pairscores_deserialize(p, buf, count);
	}

	free(buf);

	mpi_free_pairscores_type(&datatype);
}

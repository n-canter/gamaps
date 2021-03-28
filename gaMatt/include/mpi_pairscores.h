/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#ifndef MPI_PAIRSCORES_H
#define MPI_PAIRSCORES_H

#include <mpi.h>

#include "uthash.h"

#define PAIRSCORE_SIZE (2 * sizeof(int) + sizeof(double))

typedef struct pairscores_t {
	int i;
	int j;
	double score;

	UT_hash_handle hh;
} pairscores_t;


void pairscores_add(pairscores_t **h, int i, int j, const double score);
pairscores_t *pairscores_find(pairscores_t **h, int i, int j);
void pairscores_free(pairscores_t **h);

void *pairscores_serialize(pairscores_t **h);
void pairscores_deserialize(pairscores_t **h, char *buf, int count);

int mpi_init_pairscores_type(MPI_Datatype *datatype);
int mpi_free_pairscores_type(MPI_Datatype *datatype);

void mpi_pairscores_gather(pairscores_t **p, int rank, int nprocs);


#endif

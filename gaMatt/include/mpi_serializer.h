/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#ifndef SERIALIZER_H
#define SERIALIZER_H

#include <assert.h>


#include "AssemblyOrder.h"
#include "MultipleAlignment.h"

#define INIT_POOL_SIZE (1 << 20)

typedef struct mpi_serializer_pool_t {
	char *ptr;
	size_t offset;
	size_t size;
} mpi_serializer_pool_t;


void mpi_ma_serialize(mpi_serializer_pool_t *pool, MultipleAlignment *ma);
MultipleAlignment *mpi_ma_deserialize(char *buf);

struct mpi_serializer_pool_t *mpi_serializer_pool_create();
void mpi_serializer_pool_destroy(mpi_serializer_pool_t *pool);
void mpi_serializer_pool_reset(mpi_serializer_pool_t *pool);
void mpi_serializer_pool_ensure(mpi_serializer_pool_t *pool, size_t size);

#endif

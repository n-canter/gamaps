/*
 *  Copyright (c) 2017 Maksim Shegay.
 *  The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University
 *  All rights reserved.
 *
 *  parMatt is licensed under the GNU public license version 2.0.
 */

#include "mpi_serializer.h"

struct mpi_serializer_pool_t *mpi_serializer_pool_create()
{
	struct mpi_serializer_pool_t *p = calloc(1, sizeof(struct mpi_serializer_pool_t));
	assert(p);
	p->offset = 0;
	p->size = INIT_POOL_SIZE;
	p->ptr = calloc(1, p->size);
	assert(p->ptr);
	return p;
}

void mpi_serializer_pool_ensure(mpi_serializer_pool_t *pool, size_t size)
{
	if (pool->size >= size) {
		return;
	}
	pool->size = size;
	pool->ptr = realloc(pool->ptr, size);
	assert(pool->ptr);
}

void mpi_serializer_pool_destroy(mpi_serializer_pool_t *pool)
{
	free(pool->ptr);
	free(pool);
}

void mpi_serializer_pool_reset(mpi_serializer_pool_t *pool)
{
	pool->offset = 0;
}


static void pool_append(mpi_serializer_pool_t *pool, void *src, size_t size)
{
	size_t need = pool->offset + size;
	if (need >= pool->size) {
		pool->size *= 2;
		if (pool->size <= need) {
			pool->size = need * 2;
		}
		pool->ptr = realloc(pool->ptr, pool->size);
		assert(pool->ptr);
	}
	memcpy(pool->ptr + pool->offset, src, size);
	pool->offset += size;
}

static void *memdup_wrapper(char *src, int *offset, size_t size)
{
	void *out = malloc(size);
	assert(out || !size);
	memcpy(out, src + *offset, size);
	*offset += size;
	return out;
}

static char deserialize_char_value(char *src, int *offset)
{
	char ret = *(src + *offset);
	*offset += sizeof(char);
	return ret;
}

static int deserialize_int_value(char *src, int *offset)
{
	int ret = *((int *)(src + *offset));
	*offset += sizeof(int);
	return ret;
}

static double deserialize_double_value(char *src, int *offset)
{
	double ret = *((double *)(src + *offset));
	*offset += sizeof(double);
	return ret;
}

static void assembly_order_serialize(mpi_serializer_pool_t *pool, AssemblyOrder *order)
{
	size_t size = sizeof(int) + 3 * order->nodes * sizeof(int);
	int *tree = malloc(size);
	assert(tree);
	tree[0] = order->nodes;

	/* AssemblyOrder contains all the nodes in a continguous memory area */
	/* located right after AssemblyOrder */
	AssemblyNode *root = (AssemblyNode *)(order + 1);
	for (int i = 0; i < order->nodes; i++) {
		AssemblyNode *node = (AssemblyNode *)(order->root + i);
		tree[3 * i + 1] = node->id;
		tree[3 * i + 2] = node->left ? node->left - root : 0;
		tree[3 * i + 3] = node->right ? node->right - root : 0;
	}
	pool_append(pool, tree, size);
	free(tree);
}

static void assembly_order_deserialize(MultipleAlignment *ma, char *buf, int *offset)
{
	int nodes = deserialize_int_value(buf, offset);
	AssemblyOrder *order = malloc(sizeof(AssemblyOrder) + nodes * sizeof(AssemblyNode));
	assert(order);
	order->nodes = nodes;
	order->root = (AssemblyNode *)(order + 1);
	for (int i = 0; i < nodes; i++) {
		AssemblyNode *node = (AssemblyNode *)(order->root + i);

		int id = deserialize_int_value(buf, offset);
		int left = deserialize_int_value(buf, offset);
		int right = deserialize_int_value(buf, offset);

		node->id = id;
		node->left = left ? order->root + left : NULL;
		node->right = right ? order->root + right : NULL;
	}
	ma->order = order;
}

static void weighted_residue_positions_serialize(mpi_serializer_pool_t *pool, MultipleAlignment *ma)
{
	pool_append(pool, &ma->numResidues, sizeof(int));
	if (ma->numResidues == 0)
		return;
	for (int i = 0; i < ma->numChains; i++) {
		pool_append(pool, ma->residues[i].res, ma->numResidues * sizeof(ResiduePosition));
		pool_append(pool, &ma->residues[i].weight, sizeof(double));
	}
}

static void weighted_residue_positions_deserialize(MultipleAlignment *ma, char *buf, int *offset, int numChains)
{
	ma->numResidues = deserialize_int_value(buf, offset);
	if (ma->numResidues == 0) {
		ma->residues = NULL;
		return;
	}
	ma->residues = malloc(numChains * sizeof(WeightedResiduePositions));
	assert(ma->residues || !numChains);
	for (int i = 0; i < numChains; i++) {
		ma->residues[i].res = memdup_wrapper(buf, offset, ma->numResidues * sizeof(ResiduePosition));
		ma->residues[i].weight = deserialize_double_value(buf, offset);
	}
}

static void pdb_chain_serialize(mpi_serializer_pool_t *pool, PDBChain *chain)
{
	pool_append(pool, &chain->numCisPeps, sizeof(int));
	pool_append(pool, chain->cisPeps, chain->numCisPeps * sizeof(CisPep));

	pool_append(pool, &chain->numSeqAdvs, sizeof(int));
	pool_append(pool, chain->seqAdvs, chain->numSeqAdvs * sizeof(SeqAdv));

	pool_append(pool, &chain->numDbrefs, sizeof(int));
	pool_append(pool, chain->dbrefs, chain->numDbrefs * sizeof(DBRef));

	pool_append(pool, &chain->numSSbonds, sizeof(int));
	pool_append(pool, chain->ssbonds, chain->numSSbonds * sizeof(SSbond));

	pool_append(pool, &chain->numBetaSheets, sizeof(int));

	for (int i = 0; i < chain->numBetaSheets; i++) {
		pool_append(pool, chain->betaSheets->id, 4);
		pool_append(pool, &chain->betaSheets->numStrands, sizeof(int));
		pool_append(pool, chain->betaSheets->strands, chain->betaSheets->numStrands * sizeof(BetaSheetPair));
	}

	pool_append(pool, &chain->numAlphaHelices, sizeof(int));
	pool_append(pool, chain->alphaHelices, chain->numAlphaHelices * sizeof(AlphaHelix));

	pool_append(pool, &chain->numHbonds, sizeof(int));
	pool_append(pool, chain->hbonds, chain->numHbonds * sizeof(HydrogenBond));

	pool_append(pool, &chain->numBetaPairs, sizeof(int));
	pool_append(pool, chain->betaPairs, chain->numBetaPairs * sizeof(BetaPair));

	pool_append(pool, &chain->length, sizeof(int));

	pool_append(pool, &chain->numAtoms, sizeof(int));
	pool_append(pool, chain->atoms, chain->numAtoms * sizeof(Atom));

	pool_append(pool, &chain->chainName, sizeof(char));
	pool_append(pool, &chain->tempAtoms, sizeof(char));

	pool_append(pool, &chain->terminated, sizeof(int));
	pool_append(pool, &chain->numTempResidues, sizeof(int));
	pool_append(pool, chain->tempResidues, chain->numTempResidues * sizeof(TempResidue));

	pool_append(pool, &chain->secondaryCalculated, sizeof(int));

	pool_append(pool, chain->residues, chain->length * sizeof(Residue));

	size_t size = strlen(chain->idString) + 1;
	pool_append(pool, chain->idString, size);

	size = strlen(chain->seq->name) + 1;
	pool_append(pool, chain->seq->name, size);
	pool_append(pool, &chain->seq->length, sizeof(int));
	pool_append(pool, chain->seq->seq, chain->seq->length);
}

static PDBChain *pdb_chain_deserialize(char *buf, int *offset)
{
	PDBChain *chain = malloc(sizeof(PDBChain));
	assert(chain);
	chain->numCisPeps = deserialize_int_value(buf, offset);
	chain->cisPeps = memdup_wrapper(buf, offset, chain->numCisPeps * sizeof(CisPep));

	chain->numSeqAdvs = deserialize_int_value(buf, offset);
	chain->seqAdvs = memdup_wrapper(buf, offset, chain->numSeqAdvs * sizeof(SeqAdv));

	chain->numDbrefs = deserialize_int_value(buf, offset);
	chain->dbrefs = memdup_wrapper(buf, offset, chain->numDbrefs * sizeof(DBRef));

	chain->numSSbonds = deserialize_int_value(buf, offset);
	chain->ssbonds = memdup_wrapper(buf, offset, chain->numSSbonds * sizeof(SSbond));

	chain->numBetaSheets = deserialize_int_value(buf, offset);

	chain->betaSheets = malloc(chain->numBetaSheets * sizeof(BetaSheet));
	assert(chain->betaSheets || !chain->numBetaSheets);
	for (int i = 0; i < chain->numBetaSheets; i++) {
		memcpy(chain->betaSheets[i].id, buf + *offset, 4);
		*offset += 4;
		chain->betaSheets[i].numStrands = deserialize_int_value(buf, offset);
		chain->betaSheets[i].strands = memdup_wrapper(buf, offset, chain->betaSheets[i].numStrands * sizeof(BetaSheetPair));
	}

	chain->numAlphaHelices = deserialize_int_value(buf, offset);
	chain->alphaHelices = memdup_wrapper(buf, offset, chain->numAlphaHelices * sizeof(AlphaHelix));

	chain->numHbonds = deserialize_int_value(buf, offset);
	chain->hbonds = memdup_wrapper(buf, offset, chain->numHbonds * sizeof(HydrogenBond));

	chain->numBetaPairs = deserialize_int_value(buf, offset);
	chain->betaPairs = memdup_wrapper(buf, offset, chain->numBetaPairs * sizeof(BetaPair));

	chain->length = deserialize_int_value(buf, offset);

	chain->numAtoms = deserialize_int_value(buf, offset);
	chain->atoms = memdup_wrapper(buf, offset, chain->numAtoms * sizeof(Atom));

	chain->chainName = deserialize_char_value(buf, offset);
	chain->tempAtoms = deserialize_char_value(buf, offset);

	chain->terminated = deserialize_int_value(buf, offset);
	chain->numTempResidues = deserialize_int_value(buf, offset);
	chain->tempResidues = memdup_wrapper(buf, offset, chain->numTempResidues * sizeof(TempResidue));

	chain->secondaryCalculated = deserialize_int_value(buf, offset);

	chain->residues = memdup_wrapper(buf, offset, chain->length * sizeof(Residue));

	size_t size = strlen(buf + *offset) + 1;
	chain->idString = memdup_wrapper(buf, offset, size);

	chain->seq = malloc(sizeof(Sequence));
	assert(chain->seq);

	size = strlen(buf + *offset) + 1;
	chain->seq->name = memdup_wrapper(buf, offset, size);
	chain->seq->length = deserialize_int_value(buf, offset);
	chain->seq->seq = memdup_wrapper(buf, offset, chain->seq->length);

	return chain;
}

static void residue_positions_serialize(mpi_serializer_pool_t *pool, MultipleAlignment *ma)
{
	for (int i = 0; i < ma->numChains; i++) {
		pool_append(pool, &ma->chains[i]->length, sizeof(int));
		pool_append(pool, &ma->chains[i]->id, sizeof(int));
		pool_append(pool, ma->chains[i]->res, sizeof(ResiduePosition) * ma->chains[i]->length);
		pdb_chain_serialize(pool, ma->chains[i]->pdb);
	}
}

static void residue_positions_deserialize(MultipleAlignment *ma, char *buf, int *offset, int numChains)
{
	ma->chains = malloc(numChains * sizeof(ResiduePositions *));
	assert(ma->chains || !numChains);
	for (int i = 0; i < numChains; i++) {
		ma->chains[i] = malloc(sizeof(ResiduePositions));
		assert(ma->chains[i]);
		int length = deserialize_int_value(buf, offset);
		ma->chains[i]->length = length;

		ma->chains[i]->id = deserialize_int_value(buf, offset);

		ma->chains[i]->res = memdup_wrapper(buf, offset, length * sizeof(ResiduePosition));
		ma->chains[i]->pdb = pdb_chain_deserialize(buf, offset);
	}
}

void mpi_ma_serialize(mpi_serializer_pool_t *pool, MultipleAlignment *ma)
{
	pool_append(pool, &ma->score, sizeof(double));
	pool_append(pool, &ma->rmsd, sizeof(double));
	pool_append(pool, &ma->pvalue, sizeof(double));
	pool_append(pool, &ma->numChains, sizeof(int));
	pool_append(pool, &ma->numBlocks, sizeof(int));
	pool_append(pool, ma->blocks, ma->numBlocks * sizeof(AlignedBlock));

	weighted_residue_positions_serialize(pool, ma);
	residue_positions_serialize(pool, ma);
	assembly_order_serialize(pool, ma->order);
}

MultipleAlignment *mpi_ma_deserialize(char *buf)
{
	int offset = 0;
	MultipleAlignment *ma = malloc(sizeof(MultipleAlignment));
	assert(ma);

	ma->score = deserialize_double_value(buf, &offset);
	ma->rmsd = deserialize_double_value(buf, &offset);
	ma->pvalue = deserialize_double_value(buf, &offset);
	ma->numChains = deserialize_int_value(buf, &offset);
	ma->numBlocks = deserialize_int_value(buf, &offset);

	ma->blocks = memdup_wrapper(buf, &offset, ma->numBlocks * sizeof(AlignedBlock));

	weighted_residue_positions_deserialize(ma, buf, &offset, ma->numChains);
	residue_positions_deserialize(ma, buf, &offset, ma->numChains);
	assembly_order_deserialize(ma, buf, &offset);

	/* following fields will be recalculated upon alignment */
	ma->conflictMap = NULL;
	ma->averages = NULL;
	return ma;
}

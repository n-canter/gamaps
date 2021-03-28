#ifndef GA_H
#define GA_H

#include <assert.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef int (*random_int_fn)(const int n);
typedef double (*random_double_fn)();

typedef void *(*to_struct_fn)(void *context, const int id);
typedef void *(*align_fn)(void *context, void *a, void *b);
typedef void (*cleanup_fn)(void *context, void *a);
typedef double (*score_fn)(void *context, void *ma);

#define GA_TERMINATE_ON_ELAPSED 1
#define GA_TERMINATE_ON_UNCHANGED 2
#define GA_TERMINATE_ON_ITERATION 4

#define GA_ELITE_FRACTION 10
#define GA_CONVERGENCE_THRESHOLD 0.05

enum {
	ERR_BAD_ALLOC = -3,
	ERR_CROSSOVER = -4,
	ERR_FILE = -5,
	ERR_CFG = -6,
};

enum {
	  SCORE_MATT = 1,
	  SCORE_TM = 2,
};

struct population_t;

typedef struct node_t {
	struct val_t *val;

	int left;
	int right;

	int id;
	int left_edges;
	int right_edges;
} node_t;

struct node_allocator_t {
	node_t *buf;
	int len;
	int cap;
};

typedef struct guided_tree_t {
	int size;

	struct node_allocator_t *allocator;
	int root;
	double score;
	int has_score;
} guided_tree_t;



struct ga_cfg_t {
	int enabled;
	int nthreads;

	int population_size;
	int tree_size;

	int terminate_criteria;
	int terminate_on_elapsed;
	int terminate_on_unchanged;
	int terminate_on_iteration;

	int fitness_func;
	
	void *ctx;

	align_fn align;
	score_fn score;
	to_struct_fn to_struct;
	cleanup_fn cleanup;

	random_int_fn rnd_i;
	random_double_fn rnd_d;
};

int ga_parse_cfg(struct ga_cfg_t *cfg, const char *path, const int n_als, const int nthreads);

struct ga_t {
	struct ga_cfg_t cfg;

	struct guided_tree_t best;
	struct population_t *population;
};

struct ga_t *ga_init(struct ga_cfg_t cfg);
int ga_run(struct ga_t *ga);
void ga_destroy(struct ga_t *ga);
void *ga_best(struct ga_t *ga);


#endif

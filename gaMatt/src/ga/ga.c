#include <assert.h>
#include <mpi.h>
#include <omp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "MultipleAlignment.h"
#include "ga.h"
#include "mpi_helpers.h"
#include "mpi_pairscores.h"
#include "mpi_pool.h"

#define MIN_ALLOC_SIZE 1024

struct population_t {
	void *ctx;

	int size;
	int cap;

	struct guided_tree_t *trees;

	align_fn align;
	score_fn score;
	to_struct_fn to_struct;
	cleanup_fn cleanup;
};

struct val_t {
	void *ptr;
	int ref_cnt;
};

struct val_t *val_create(void *ptr)
{
	struct val_t *ret = calloc(1, sizeof(struct val_t));
	if (ret == NULL) {
		return NULL;
	}
	ret->ptr = ptr;
	ret->ref_cnt = 1;
	return ret;
}

static void val_ref(struct val_t *val)
{
	if (val == NULL) {
		return;
	}
	val->ref_cnt++;
}

static void val_unref(void *context, struct val_t *val, cleanup_fn cleanup)
{
	if (val == NULL) {
		return;
	}

	val->ref_cnt--;
	if (val->ref_cnt == 0) {
		cleanup(context, val->ptr);
		free(val);
	}
}

struct node_allocator_t *node_allocator_new()
{
	struct node_allocator_t *al = calloc(1, sizeof(struct node_allocator_t));
	if (al == NULL) {
		return NULL;
	}

	al->buf = calloc(MIN_ALLOC_SIZE, sizeof(node_t));
	if (al->buf == NULL) {
		return NULL;
	}
	al->len = 0;
	al->cap = MIN_ALLOC_SIZE;


	return al;
}

void node_allocator_destroy(struct node_allocator_t *al)
{
	if (al == NULL) {
		return;
	}
	free(al->buf);
	free(al);
}

int node_allocator_alloc(struct node_allocator_t *al, const int id)
{
	if (al->len >= al->cap) {
		al->cap *= 2;
		al->buf = realloc(al->buf, al->cap * sizeof(node_t));
		if (al->buf == NULL) {
			return ERR_BAD_ALLOC;
		}
	}

	node_t node = {
	    .val = NULL,
	    .left = -1,
	    .right = -1,
	    .left_edges = 0,
	    .right_edges = 0,
	    .id = id,
	};
	al->buf[al->len] = node;
	return al->len++;
}

void do_print(struct node_allocator_t *al, int idx, char *padding, char *buf, int has_right_sibling)
{
	node_t node = al->buf[idx];
	printf("\n");
	printf("%s", padding);
	printf("%s", buf);
	if (has_right_sibling) {
		printf("%d:%d:%d(l", node.id, node.left_edges, node.right_edges);
		printf(")");
	}
	else {
		printf("%d:%d:%d(r", node.id, node.left_edges, node.right_edges);
		printf(")");
	}
	char padding1[1024] = {0};
	if (has_right_sibling == 1) {
		sprintf(padding1, "%s│  ", padding);
	}
	else {
		sprintf(padding1, "%s   ", padding);
	}

	char *buf_left = "├──";
	char *buf_right = "└──";
	if (node.right == -1) {
		buf_left = "└──";
	}

	if (node.left != -1) {
		do_print(al, node.left, padding1, buf_left, node.right != -1);
	}
	if (node.right != -1) {
		do_print(al, node.right, padding1, buf_right, 0);
	}
}

static void printer(guided_tree_t t)
{
	if (t.root == -1) {
		return;
	}
	node_t root = t.allocator->buf[t.root];
	printf("%d:(%d:%d)", root.id, root.left_edges, root.right_edges);
	char *buf_left = "├──";
	char *buf_right = "└──";
	if (root.right == -1) {
		buf_left = "└──";
	}
	if (root.left != -1) {
		do_print(t.allocator, root.left, "", buf_left, root.right != -1);
	}
	if (root.right != -1) {
		do_print(t.allocator, root.right, "", buf_right, 0);
	}
	printf("\n");
}

int guided_tree_split_edge(guided_tree_t *t, const int from, const int to, const int id)
{
	int new_node = node_allocator_alloc(t->allocator, -1);
	if (new_node < 0) {
		return new_node;
	}

	int new_leaf = node_allocator_alloc(t->allocator, id);
	if (new_leaf < 0) {
		return new_leaf;
	}

	node_t from_node = t->allocator->buf[from];
	if (from_node.left == to) {
		from_node.left = new_node;
		from_node.left_edges += 2;
	}
	else {
		from_node.right = new_node;
		from_node.right_edges += 2;
	}

	t->allocator->buf[from] = from_node;
	t->allocator->buf[new_node].left = to;
	t->allocator->buf[new_node].right = new_leaf;
	t->allocator->buf[new_node].right_edges = 1;
	t->allocator->buf[new_node].left_edges = t->allocator->buf[to].left_edges + t->allocator->buf[to].right_edges + 1;

	return 0;
}

int guided_tree_split(guided_tree_t *t, const int root, const int edge, const int id)
{
	node_t node = t->allocator->buf[root];
	if (edge < node.left_edges) {
		if (edge == node.left_edges - 1) {
			return guided_tree_split_edge(t, root, node.left, id);
		}
		else {
			t->allocator->buf[root].left_edges += 2;
			return guided_tree_split(t, node.left, edge, id);
		}
	}
	else {
		if (edge == node.left_edges + node.right_edges - 1) {
			return guided_tree_split_edge(t, root, node.right, id);
		}
		else {
			t->allocator->buf[root].right_edges += 2;
			return guided_tree_split(t, node.right, edge - node.left_edges, id);
		}
	}
}

int guided_tree_insert_random(guided_tree_t *t, random_int_fn rnd, const int id)
{
	if (t->root == -1) {
		int root = node_allocator_alloc(t->allocator, id);
		if (root < 0) {
			return root;
		}
		t->root = root;
		return 0;
	}

	node_t root = t->allocator->buf[t->root];
	int edges = root.left_edges + root.right_edges;
	int roll = rnd(edges + 1);
	if (roll == edges) {
		int new_root = node_allocator_alloc(t->allocator, -1);
		if (new_root < 0) {
			return new_root;
		}

		int new_node = node_allocator_alloc(t->allocator, id);
		if (new_node < 0) {
			return new_node;
		}

		t->allocator->buf[new_root].left = t->root;
		t->allocator->buf[new_root].right = new_node;
		t->allocator->buf[new_root].left_edges = root.left_edges + root.right_edges + 1;
		t->allocator->buf[new_root].right_edges = 1;
		t->root = new_root;

		return 0;
	}

	return guided_tree_split(t, t->root, roll, id);
}

int guided_tree_random(guided_tree_t *t, const int size, random_int_fn rnd)
{
	struct node_allocator_t *allocator = node_allocator_new();
	if (allocator == NULL) {
		return -1;
	}
	t->size = size;
	t->allocator = allocator;
	t->root = -1;

	for (int i = 0; i < size; i++) {
		int rc = guided_tree_insert_random(t, rnd, i);
		if (rc != 0) {
			return rc;
		}
	}

	return 0;
}

void guided_tree_cleanup(void *context, guided_tree_t t, const node_t node, cleanup_fn cleanup)
{
	if (node.left != -1) {
		guided_tree_cleanup(context, t, t.allocator->buf[node.left], cleanup);
	}

	val_unref(context, node.val, cleanup);

	if (node.right != -1) {
		guided_tree_cleanup(context, t, t.allocator->buf[node.right], cleanup);
	}
}

void guided_tree_destroy(void *context, guided_tree_t t, cleanup_fn cleanup)
{
	if (t.root != -1) {
		guided_tree_cleanup(context, t, t.allocator->buf[t.root], cleanup);
	}
	node_allocator_destroy(t.allocator);
}

int traverse_index;

void guided_tree_traverse(const guided_tree_t t, const node_t node, int *dst)
{
	if (node.left != -1) {
		guided_tree_traverse(t, t.allocator->buf[node.left], dst);
	}
	if (node.id != -1) {
		dst[traverse_index] = node.id;
		traverse_index++;
	}
	if (node.right != -1) {
		guided_tree_traverse(t, t.allocator->buf[node.right], dst);
	}
}

int *guided_tree_order(const guided_tree_t t)
{
	int *dst = calloc(t.size, sizeof(int));
	if (dst == NULL) {
		return NULL;
	}

	traverse_index = 0;

	guided_tree_traverse(t, t.allocator->buf[t.root], dst);
	return dst;
}

int *crossover_pmx(const int *a, const int *b, const int size, const int m, const int n)
{
	int *ret = calloc(size, sizeof(int));
	if (ret == NULL) {
		return NULL;
	}

	int *set = calloc(size, sizeof(int));
	if (set == NULL) {
		return NULL;
	}

	int *mapping = calloc(size, sizeof(int));
	if (mapping == NULL) {
		return NULL;
	}

	for (int i = m; i < n; i++) {
		mapping[b[i]] = i;
		set[b[i]] = 1;
		ret[i] = b[i];
	}


	for (int i = 0; i < m; i++) {
		if (set[a[i]] == 0) {
			ret[i] = a[i];
		}
		else {
			int v = a[mapping[a[i]]];
			for (; set[v] != 0; v = a[mapping[v]]) {
			}
			ret[i] = v;
		}
	}

	for (int i = n; i < size; i++) {
		if (set[a[i]] == 0) {
			ret[i] = a[i];
		}
		else {
			int v = a[mapping[a[i]]];
			for (; set[v] != 0; v = a[mapping[v]]) {
			}
			ret[i] = v;
		}
	}

	free(set);
	free(mapping);
	return ret;
}

int guided_tree_permutate_do(void *context, guided_tree_t *t, int src_node, const int *src, const int *dst, cleanup_fn cleanup)
{
	node_t node = t->allocator->buf[src_node];
	int ret = 0;
	if (node.left != -1) {
		ret += guided_tree_permutate_do(context, t, node.left, src, dst, cleanup);
	}

	if (node.id != -1) {
		if (src[traverse_index] != dst[traverse_index]) {
			node.id = dst[traverse_index];
			ret = 1;
		}
		traverse_index++;
	}

	if (node.right != -1) {
		ret += guided_tree_permutate_do(context, t, node.right, src, dst, cleanup);
	}
	if (ret > 0) {
		val_unref(context, node.val, cleanup);
		node.val = NULL;
	}
	t->allocator->buf[src_node] = node;
	return ret;
}

int guided_tree_permutate(void *context, guided_tree_t *t, const int *src, const int *dst, cleanup_fn cleanup)
{
	traverse_index = 0;
	t->has_score = 0;
	return guided_tree_permutate_do(context, t, t->root, src, dst, cleanup);
}

int guided_tree_clone_do(const guided_tree_t src, node_t src_node, guided_tree_t *dst)
{
	int new_node = node_allocator_alloc(dst->allocator, src_node.id);
	node_t dst_node = dst->allocator->buf[new_node];
	val_ref(src_node.val);
	dst_node = src_node;
	if (src_node.left != -1) {
		dst_node.left = guided_tree_clone_do(src, src.allocator->buf[src_node.left], dst);
	}
	if (src_node.right != -1) {
		dst_node.right = guided_tree_clone_do(src, src.allocator->buf[src_node.right], dst);
	}
	dst->allocator->buf[new_node] = dst_node;
	return new_node;
}

int guided_tree_clone(const guided_tree_t src, guided_tree_t *dst)
{
	dst->size = src.size;
	dst->allocator = node_allocator_new();
	dst->has_score = src.has_score;
	dst->score = src.score;
	if (dst->allocator == NULL) {
		return -1;
	}

	dst->root = guided_tree_clone_do(src, src.allocator->buf[src.root], dst);
	return 0;
}

void random_pair(const int range, int *a, int *b, random_int_fn rnd)
{
	int m = rnd(range);
	int n = rnd(range);
	while (n == m) {
		n = rnd(range);
	}
	if (m > n) {
		int t = n;
		n = m;
		m = t;
	}
	*a = m;
	*b = n;
}

int guided_tree_crossover(void *context, const guided_tree_t a, const guided_tree_t b, guided_tree_t *child_a, guided_tree_t *child_b, random_int_fn rnd, cleanup_fn cleanup)
{
	int *perm_a = guided_tree_order(a);
	if (perm_a == NULL) {
		return -1;
	}
	int *perm_b = guided_tree_order(b);
	if (perm_b == NULL) {
		return -1;
	}

	int rc = guided_tree_clone(a, child_a);
	if (rc != 0) {
		return rc;
	}

	rc = guided_tree_clone(b, child_b);
	if (rc != 0) {
		return rc;
	}
	int m, n;
	random_pair(a.size, &m, &n, rnd);

	int *dst = crossover_pmx(perm_a, perm_b, a.size, m, n);
	guided_tree_permutate(context, child_a, perm_a, dst, cleanup);
	free(dst);

	dst = crossover_pmx(perm_b, perm_a, a.size, m, n);
	guided_tree_permutate(context, child_b, perm_b, dst, cleanup);
	free(dst);

	free(perm_a);
	free(perm_b);
	return 0;
}

int guided_tree_mutate(void *context, guided_tree_t *src, random_int_fn rnd, cleanup_fn cleanup)
{
	int *old = guided_tree_order(*src);
	if (src == NULL) {
		return -1;
	}
	int *new = calloc(src->size, sizeof(int));
	memcpy(new, old, src->size * sizeof(int));

	int m, n;
	random_pair(src->size, &m, &n, rnd);
	new[m] = old[n];
	new[n] = old[m];

	guided_tree_permutate(context, src, old, new, cleanup);

	free(new);
	free(old);

	return 0;
}

void *guided_tree_node_align(void *ctx, const guided_tree_t *t, const int root, to_struct_fn to_struct, align_fn align, cleanup_fn cleanup)
{
	node_t node = t->allocator->buf[root];
	if (node.val != NULL) {
		return node.val->ptr;
	}
	if (node.left == -1 && node.right == -1) {
		return to_struct(ctx, node.id);
	}
	else if (node.left != -1 && node.right != -1) {
		void *left = guided_tree_node_align(ctx, t, node.left, to_struct, align, cleanup);
		if (left == NULL) {
			return NULL;
		}

		void *right = guided_tree_node_align(ctx, t, node.right, to_struct, align, cleanup);
		if (right == NULL) {
			return NULL;
		}

		void *ma = align(ctx, left, right);
		if (ma == NULL) {
			return NULL;
		}
		node.val = val_create(ma);
	}
	t->allocator->buf[root] = node;
	return node.val->ptr;
}

double guided_tree_score(void *ctx, guided_tree_t *t, to_struct_fn to_struct, align_fn align, score_fn score, cleanup_fn cleanup)
{
	if (t->has_score) {
		return t->score;
	}
	/* printer(*t); */
	void *ma = guided_tree_node_align(ctx, t, t->root, to_struct, align, cleanup);

	t->score = score(ctx, ma);

	t->has_score = 1;
	return t->score;
}

void *ga_best(struct ga_t *ga)
{
	if (ga->best.root == -1) {
		return NULL;
	}
	struct ga_cfg_t cfg = ga->cfg;
	MultipleAlignment *ma = guided_tree_node_align(cfg.ctx, &ga->best, ga->best.root, cfg.to_struct, cfg.align, cfg.cleanup);
	return ma;
}

static struct population_t *population_init_empty(void *ctx, const int size, to_struct_fn to_struct, align_fn align, score_fn score, cleanup_fn cleanup)
{
	struct population_t *p = calloc(1, sizeof(struct population_t));
	if (p == NULL) {
		return NULL;
	}

	p->ctx = ctx;
	p->score = score;
	p->align = align;
	p->to_struct = to_struct;
	p->cleanup = cleanup;

	p->size = 0;
	p->cap = size;
	p->trees = calloc(size, sizeof(*p->trees));
	if (p->trees == NULL) {
		free(p);
		return NULL;
	}

	return p;
}

static struct population_t *population_init(void *ctx, const int size, const int tree_size, random_int_fn rnd, to_struct_fn to_struct, align_fn align, score_fn score, cleanup_fn cleanup)
{
	struct population_t *p = population_init_empty(ctx, size, to_struct, align, score, cleanup);
	if (p == NULL) {
		return NULL;
	}
	p->size = size;
	for (int i = 0; i < size; i++) {
		int rc = guided_tree_random(p->trees + i, tree_size, rnd);
		if (rc != 0) {
			free(p);
			free(p->trees);
			return NULL;
		}
	}

	return p;
}

static void trees_shuffle(struct guided_tree_t *arr, const int size, random_int_fn rnd)
{
	for (int i = size - 1; i > 0; i--) {
		int j = rnd(i + 1);
		struct guided_tree_t t = arr[i];
		arr[i] = arr[j];
		arr[j] = t;
	}
}

static void population_shuffle(struct population_t *p, const int offset, random_int_fn rnd)
{
	trees_shuffle(p->trees + offset, p->size, rnd);
}

static int population_add(struct population_t *p, const guided_tree_t t)
{
	if (p->size == p->cap) {
		struct guided_tree_t *m = realloc(p->trees, p->size * 2 * sizeof(*p->trees));
		if (m == NULL) {
			return ERR_BAD_ALLOC;
		}
		p->trees = m;
		p->cap = p->size * 2;
	}

	p->trees[p->size] = t;
	p->size++;
	return 0;
}

static double crossover_prob(const double a_score, const double b_score, const double max_score, const double avg_score)
{
	double m = (a_score + b_score) / 2;
	if (m < avg_score) {
		return 0.05;
	}
	return (m - avg_score) / (max_score - avg_score);
}

static int population_crossover(struct population_t *p, const double max_score, const double avg_score, random_int_fn rnd_i, random_double_fn rnd_d)
{
	population_shuffle(p, 0, rnd_i);

	int n = p->size / 2;
	for (int i = 0; i < n; i++) {
		guided_tree_t parent_a = p->trees[2 * i];
		guided_tree_t parent_b = p->trees[2 * i + 1];
		double prob = crossover_prob(parent_a.score, parent_b.score, max_score, avg_score);
		if (prob <= rnd_d()) {
			continue;
		}
		guided_tree_t child_a;
		guided_tree_t child_b;
		int rc = guided_tree_crossover(p->ctx, parent_a, parent_b, &child_a, &child_b, rnd_i, p->cleanup);
		if (rc != 0 && rc != ERR_CROSSOVER) {
			return rc;
		}
		if (rc == 0) {
			rc = population_add(p, child_a);
			if (rc != 0) {
				return rc;
			}
			rc = population_add(p, child_b);
			if (rc != 0) {
				return rc;
			}
		}
	}

	return 0;
}

static double mutation_prob(const double score, const double max_score, const double avg_score)
{
	if (score < avg_score) {
		return 0.5;
	}
	double p = 0.5 * (max_score - score) / (max_score - avg_score);
	if (p < 0.01) {
		p = 0.01;
	}
	return p;
}

static int population_mutation(struct population_t *p, const double max_score, const double avg_score, random_int_fn rnd_i, random_double_fn rnd_d)
{
	for (int i = 0; i < p->size; i++) {
		double prob = mutation_prob(p->trees[i].score, max_score, avg_score);
		if (prob <= rnd_d()) {
			continue;
		}
		guided_tree_mutate(p->ctx, p->trees + i, rnd_i, p->cleanup);
	}
	return 0;
}

struct score_t {
	double v;
	int inited;
};

struct ga_task_ctx_t {
	void *ctx;
	pthread_mutex_t lock;
	struct population_t *p;
	struct score_t *scores;
};

static int mpi_population_fitness(void *ctx, struct population_t *p, const int iteration, const int rank, const int nprocs, const int nthreads, double *max, double *avg)
{
	struct score_t *local_scores = calloc(p->size, sizeof(struct score_t));
	if (local_scores == NULL) {
		return ERR_BAD_ALLOC;
	}

	pthread_mutex_t lock;
	pthread_mutex_init(&lock, NULL);
	struct ga_task_ctx_t context = {
	    .ctx = ctx,
	    .lock = lock,
	    .p = p,
	    .scores = local_scores,
	};

	snapshot_t *snap = snapshot_init(nprocs);

	snapshot_fill_ga(&context, snap, p->size);

	void *align(mpi_task_t t)
	{
		struct ga_task_ctx_t *ctx = t.ctx;
		struct population_t *p = ctx->p;
		pthread_mutex_t l = ctx->lock;

		pthread_mutex_lock(&l);
		int inited = ctx->scores[t.al1].inited;
		if (inited == 0) {
			ctx->scores[t.al1].inited = 1;
		}
		pthread_mutex_unlock(&l);

		if (inited == 1) {
			return NULL;
		}
		double score = 0;
		if (p->trees[t.al1].has_score) {
			score = guided_tree_score(p->ctx, p->trees + t.al1, p->to_struct, p->align, p->score, p->cleanup);
		}
		else {
			double t1 = MPI_Wtime();
			score = guided_tree_score(p->ctx, p->trees + t.al1, p->to_struct, p->align, p->score, p->cleanup);
			double t2 = MPI_Wtime();
			mpi_log(0, 0, "genetic %d: aligned %d in %lf seconds: %lf", iteration, t.al1, t2 - t1, score);
		}
		pthread_mutex_lock(&l);
		ctx->scores[t.al1].v = score;
		ctx->scores[t.al1].inited = 1;
		pthread_mutex_unlock(&l);
		return NULL;
	}

	taskpool_ctx_t task_ctx = {
	    .processed = NULL,
	    .snap = snap,
	    .rank = rank,
	};
	taskpool_run(nthreads, task_ctx, mpi_task_updater, align);

	double *local = calloc(p->size, sizeof(double));
	if (local == NULL) {
		return ERR_BAD_ALLOC;
	}

	for (int i = 0; i < p->size; i++) {
		local[i] = local_scores[i].v;
	}

	double *global_scores = calloc(p->size, sizeof(double));
	if (global_scores == NULL) {
		return ERR_BAD_ALLOC;
	}

	MPI_Allreduce(local, global_scores, p->size, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

	double _max = 0;
	double _sum = 0;
	for (int i = 0; i < p->size; i++) {
		double s = global_scores[i];
		if (s > _max) {
			_max = s;
		}
		_sum += s;

		p->trees[i].score = s;
		p->trees[i].has_score = 1;
	}
	*max = _max;
	*avg = _sum / p->size;
	pthread_mutex_destroy(&lock);
	free(local_scores);
	free(global_scores);
	snapshot_destroy(snap);
	return 0;
}

static int population_fitness(void *ctx, struct population_t *p, const int rank, const int nprocs, const int nthreads, double *max, double *avg)
{
	double _max = 0;
	double _sum = 0;
	omp_set_num_threads(nthreads);
#pragma omp parallel
	{
#pragma omp for schedule(dynamic) reduction(+                     \
					    : _sum) reduction(max \
							      : _max)
		for (int i = 0; i < p->size; i++) {
			double begin = omp_get_wtime();
			double score = guided_tree_score(p->ctx, p->trees + i, p->to_struct, p->align, p->score, p->cleanup);
			double end = omp_get_wtime();
			_sum += score;
			if (score > _max) {
				_max = score;
			}
#pragma omp critical
			{
				printf("elapsed %lf (%lf)\n", end - begin, p->trees[i].score);
				printer(p->trees[i]);
			}
		}
	}
	*max = _max;
	*avg = _sum / p->size;
	return 0;
}

static int cmp_scores(const void *_a, const void *_b)
{
	double a = ((struct guided_tree_t *)_a)->score;
	double b = ((struct guided_tree_t *)_b)->score;
	if (a > b) {
		return -1;
	}
	if (a < b) {
		return +1;
	}
	return 0;
}

static int population_selection(struct population_t *p, struct population_t *ret, const int target, random_int_fn rnd_i, random_double_fn rnd_d, double *max, double *avg)
{
	qsort(p->trees, p->size, sizeof(struct guided_tree_t), cmp_scores);
	int offset = p->size / GA_ELITE_FRACTION + (p->size < GA_ELITE_FRACTION);
	for (int i = 0; i < offset; i++) {
		struct guided_tree_t t;
		int rc = guided_tree_clone(p->trees[i], &t);
		if (rc != 0) {
			return rc;
		}
		population_add(ret, t);
	}

	struct guided_tree_t *loosers = calloc(p->size - offset, sizeof(struct guided_tree_t));
	if (loosers == NULL) {
		return ERR_BAD_ALLOC;
	}

	int loosers_offset = 0;
	int loosers_size = p->size - offset;
	memcpy(loosers, p->trees + offset, sizeof(*loosers) * (p->size - offset));

	for (; ret->size < target;) {
		if (loosers_size == 1) {
			struct guided_tree_t t;
			int rc = guided_tree_clone(*loosers, &t);
			if (rc != 0) {
				return rc;
			}
			population_add(ret, t);
			break;
		}

		trees_shuffle(loosers, loosers_size, rnd_i);
		for (int i = 0; i < loosers_size / 2 && ret->size < target; i++) {
			double a = loosers[2 * i].score;
			double b = loosers[2 * i + 1].score;
			struct guided_tree_t *winner;
			if (a / (a + b) < rnd_d()) {
				winner = loosers + 2 * i;
				loosers[loosers_offset++] = loosers[2 * i + 1];
			}
			else {
				winner = loosers + 2 * i + 1;
				loosers[loosers_offset++] = loosers[2 * i];
			}
			struct guided_tree_t t;
			int rc = guided_tree_clone(*winner, &t);
			if (rc != 0) {
				return rc;
			}
			population_add(ret, t);
		}
		if (loosers_size % 2 == 1) {
			loosers[loosers_offset++] = loosers[loosers_size - 1];
		}
		loosers_size = loosers_offset;
		loosers_offset = 0;
	}

	free(loosers);

	double _max = 0;
	double _sum = 0;
	for (int i = 0; i < ret->size; i++) {
		double s = ret->trees[i].score;
		if (s > _max) {
			_max = s;
		}
		_sum += s;

		p->trees[i].score = s;
		p->trees[i].has_score = 1;
	}
	*max = _max;
	*avg = _sum / ret->size; 
	return 0;
}

static guided_tree_t population_find_best(struct population_t *p)
{
	double best_score = guided_tree_score(p->ctx, p->trees, p->to_struct, p->align, p->score, p->cleanup);
	int best_idx = 0;
	for (int i = 1; i < p->size; i++) {
		double score = guided_tree_score(p->ctx, p->trees + i, p->to_struct, p->align, p->score, p->cleanup);
		if (score > best_score) {
			best_score = score;
			best_idx = i;
		}
	}

	return p->trees[best_idx];
}

static void population_destroy(struct population_t *p)
{
	for (int i = 0; i < p->size; i++) {
		guided_tree_destroy(p->ctx, p->trees[i], p->cleanup);
	}
	if (p->trees != NULL) {
		free(p->trees);
		p->trees = NULL;
	}
	free(p);
}

static int ga_terminate(struct ga_cfg_t cfg, const time_t elapsed, const int unchanged, const int iteration)
{
	if ((cfg.terminate_criteria & GA_TERMINATE_ON_ELAPSED) && (elapsed >= cfg.terminate_on_elapsed)) {
		mpi_log(0, 1, "terminate on elapsed");
		return 1;
	}

	if ((cfg.terminate_criteria & GA_TERMINATE_ON_UNCHANGED) && (unchanged >= cfg.terminate_on_unchanged)) {
		mpi_log(0, 1, "terminate on unchanged");
		return 1;
	}

	if ((cfg.terminate_criteria & GA_TERMINATE_ON_ITERATION) && (iteration >= cfg.terminate_on_iteration)) {
		mpi_log(0, 1, "terminate on iteration");
		return 1;
	}

	return 0;
}

static int population_migrate_in(struct population_t *p, const int tree_size, const int n, random_int_fn rnd_i)
{
	for (int i = 0; i < n; i++) {
		struct guided_tree_t t;
		int rc = guided_tree_random(&t, tree_size, rnd_i);
		if (rc != 0) {
			return rc;
		}
		rc = population_add(p, t);
		if (rc != 0) {
			return rc;
		}
	}
	return 0;
}

static int ga_find_best(struct ga_t *ga, struct population_t *next_gen, double *best_score)
{
	struct ga_cfg_t cfg = ga->cfg;
	struct guided_tree_t next_best = population_find_best(next_gen);
	double next_score = guided_tree_score(cfg.ctx, &next_best, cfg.to_struct, cfg.align, cfg.score, cfg.cleanup);
	double prev_score = -1.0;
	if (ga->best.root != -1) {
		prev_score = guided_tree_score(cfg.ctx, &ga->best, cfg.to_struct, cfg.align, cfg.score, cfg.cleanup);
	}

	if (prev_score < next_score) {
		*best_score = next_score;
		guided_tree_destroy(cfg.ctx, ga->best, cfg.cleanup);
		int rc = guided_tree_clone(next_best, &ga->best);
		if (rc != 0) {
			return rc;
		}
	}

	return 0;
}

static int ga_converged(double max_score, double avg_score)
{
	if (max_score == 0) {
		return 1;
	}
	return ((max_score - avg_score) / max_score) < GA_CONVERGENCE_THRESHOLD;
}

static int ga_iteration(struct ga_t *ga, const int iteration, const int rank, const int nprocs, double *max_score, double *avg_score, double *best_score)
{
	struct ga_cfg_t cfg = ga->cfg;
	struct population_t *prev_gen = ga->population;

	if (iteration > 0 && ga_converged(*max_score, *avg_score)) {
		population_migrate_in(prev_gen, cfg.tree_size, cfg.population_size, cfg.rnd_i);
	}

	int rc;
	if (iteration == 0 || ga_converged(*max_score, *avg_score)) {
		rc = mpi_population_fitness(cfg.ctx, prev_gen, iteration, rank, nprocs, cfg.nthreads, max_score, avg_score);
		if (rc != 0) {
			return rc;
		}
	}

	struct population_t *next_gen = population_init_empty(cfg.ctx, cfg.population_size, cfg.to_struct, cfg.align, cfg.score, cfg.cleanup);
	if (next_gen == NULL) {
		return ERR_BAD_ALLOC;
	}

	rc = population_selection(prev_gen, next_gen, cfg.population_size, cfg.rnd_i, cfg.rnd_d, max_score, avg_score);
	if (rc != 0) {
		return rc;
	}

	rc = population_crossover(next_gen, *max_score, *avg_score, cfg.rnd_i, cfg.rnd_d);
	if (rc != 0) {
		return rc;
	}
	rc = mpi_population_fitness(cfg.ctx, next_gen, iteration, rank, nprocs, cfg.nthreads, max_score, avg_score);
	if (rc != 0) {
		return rc;
	}

	rc = population_mutation(next_gen, *max_score, *avg_score, cfg.rnd_i, cfg.rnd_d);
	if (rc != 0) {
		return rc;
	}
	rc = mpi_population_fitness(cfg.ctx, next_gen, iteration, rank, nprocs, cfg.nthreads, max_score, avg_score);
	if (rc != 0) {
		return rc;
	}

	rc = ga_find_best(ga, next_gen, best_score);
	if (rc != 0) {
		return rc;
	}
	ga->population = next_gen;
	population_destroy(prev_gen);

	return 0;
}

struct ga_t *ga_init(struct ga_cfg_t cfg)
{
	struct ga_t *ga = calloc(1, sizeof(struct ga_t));
	if (ga == NULL) {
		return NULL;
	}
	ga->cfg = cfg;
	ga->population = population_init(cfg.ctx, cfg.population_size, cfg.tree_size, cfg.rnd_i, cfg.to_struct, cfg.align, cfg.score, cfg.cleanup);
	if (ga->population == NULL) {
		free(ga);
		return NULL;
	}
	ga->best.root = -1;
	return ga;
}

int ga_run(struct ga_t *ga)
{
	int rank, nprocs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

	time_t start = time(NULL);
	int unchanged = 0;
	double best_score = 0;
	double avg_score = 0;
	double max_score = 0;
	for (int iteration = 0;; iteration++) {
		double prev_best_score = best_score;
		int rc = ga_iteration(ga, iteration, rank, nprocs, &max_score, &avg_score, &best_score);
		if (rc != 0) {
			return rc;
		}
		if (prev_best_score == best_score) {
			unchanged++;
		}
		else {
			unchanged = 0;
			if(rank == MASTER) {
				printer(ga->best);
			}
		}
		time_t elapsed = time(NULL) - start;
		if (ga_converged(max_score, avg_score)) {
			mpi_log(0, 1, "genetic %d: elapsed=%lds avg=%lf max=%lf best=%lf unchanged=%d converged", iteration, elapsed, avg_score, max_score, best_score, unchanged);
		}
		else {
			mpi_log(0, 1, "genetic %d: elapsed=%lds avg=%lf max=%lf best=%lf unchanged=%d", iteration, elapsed, avg_score, max_score, best_score, unchanged);
		}

		if (ga_terminate(ga->cfg, elapsed, unchanged, iteration)) {
			break;
		}
	}
	return 0;
}

void ga_destroy(struct ga_t *ga)
{
	// TODO: avoid leak. Currently ga->best is not being freed.
	// guided_tree_destroy(ga->cfg.ctx, ga->best, ga->cfg.cleanup);
	population_destroy(ga->population);
	free(ga);
}

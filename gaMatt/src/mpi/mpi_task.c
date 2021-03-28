/* /\* Copyright (c) 2017 Maksim Shegay. *\/ */
/* /\* The faculty of Computational Mathematics and Cybernetics of Lomonosov Moscow State University *\/ */
/* /\* All rights reserved. *\/ */

/* /\* parMatt is licensed under the GNU public license version 2.0. *\/ */

/* #include <assert.h> */
/* #include <stdio.h> */
/* #include <stdlib.h> */
/* /\* #include "mpi_helpers.h" *\/ */

/* typedef struct mpi_task_t { */
/* 	/\* alignments IDs al2 is ingored in iterative phase *\/ */
/* 	int al1; */
/* 	int al2; */
/* } mpi_task_t; */

/* mpi_task_t mpi_create_task(const int al1, const int al2) */
/* { */
/* 	mpi_task_t t = { */
/* 	    .al1 = al1, */
/* 	    .al2 = al2, */
/* 	}; */
/* 	return t; */
/* } */


/* struct bucket_t { */
/* 	int len; */
/* 	int cap; */
/* 	int start; */
/* 	int rank; */
/* 	mpi_task_t *tasks; */
/* }; */

/* static void bucket_init(struct bucket_t *b, const int rank) */
/* { */
/* 	b->len = 0; */
/* 	b->cap = 1; */
/* 	b->start = 0; */
/* 	b->rank = rank; */
/* 	b->tasks = calloc(1, sizeof(mpi_task_t *)); */
/* 	assert(b->tasks); */
/* } */

/* static int bucket_count(struct bucket_t *b) */
/* { */
/* 	return b->len - b->start; */
/* } */

/* static void bucket_reset(struct bucket_t *b, const int rank) */
/* { */
/* 	b->len = 0; */
/* 	b->start = 0; */
/* 	b->rank = rank; */
/* } */

/* static void bucket_push(struct bucket_t *b, mpi_task_t task) */
/* { */
/* 	if (b->len >= b->cap) { */
/* 		b->cap *= 2; */
/* 		b->tasks = realloc(b->tasks, b->cap * sizeof(mpi_task_t *)); */
/* 		assert(b->tasks); */
/* 	} */
/* 	b->tasks[b->len++] = task; */
/* } */

/* static mpi_task_t bucket_pop(struct bucket_t *b) */
/* { */
/* 	assert(b->len - b->start > 0); */
/* 	return b->tasks[--b->len]; */
/* } */

/* static mpi_task_t bucket_pop_front(struct bucket_t *b) */
/* { */
/* 	assert(b->len - b->start > 0); */
/* 	return b->tasks[b->start++]; */
/* } */

/* static void bucket_leave_last(struct bucket_t *b, const int n) */
/* { */
/* 	assert(b->len - n >= b->start); */
/* 	b->start = b->len - n; */
/* } */

/* typedef struct snapshot_t { */
/* 	int n; */
/* 	struct bucket_t *procs; */
/* } snapshot_t; */

/* typedef struct diff_t { */
/* 	int empty; */
/* 	int deleted; */
/* 	int appended; */
/* 	mpi_task_t *tasks; */
/* } diff_t; */

/* void diff_destroy(diff_t *diff) */
/* { */
/* 	if (diff->tasks) { */
/* 		free(diff->tasks); */
/* 	} */
/* } */

/* static void snapshot_reset(snapshot_t *s) */
/* { */
/* 	for (int i = 0; i < s->n; i++) { */
/* 		bucket_reset(s->procs + i, i); */
/* 	} */
/* } */

/* snapshot_t *snapshot_init(const int n_procs) */
/* { */
/* 	snapshot_t *s = calloc(1, sizeof(snapshot_t)); */
/* 	assert(s); */

/* 	s->procs = calloc(n_procs, sizeof(struct bucket_t)); */
/* 	assert(s->procs); */
/* 	for (int i = 0; i < n_procs; i++) { */
/* 		bucket_init(s->procs + i, i); */
/* 	} */

/* 	s->n = n_procs; */
/* 	return s; */
/* } */
/* void snapshot_fill_pairwise(snapshot_t *s, const int n_als) */
/* { */
/* 	snapshot_reset(s); */
/* 	int ind = 0; */
/* 	for (int i = 0; i < n_als; i++) { */
/* 		for (int j = 0; j < i; j++) { */
/* 			mpi_task_t task = mpi_create_task(i, j); */
/* 			bucket_push(s->procs + ind % s->n, task); */
/* 			ind++; */
/* 		} */
/* 	} */
/* } */

/* void snapshot_fill_iterative(snapshot_t *s, const int n_als) */
/* { */
/* 	snapshot_reset(s); */
/* 	for (int i = 0; i < n_als; i++) { */
/* 		mpi_task_t task = mpi_create_task(i, -1); */
/* 		bucket_push(s->procs + i % s->n, task); */
/* 	} */
/* } */

/* void snapshot_destroy(snapshot_t *s) */
/* { */
/* 	for (int i = 0; i < s->n; i++) { */
/* 		free(s->procs[i].tasks); */
/* 	} */
/* 	free(s->procs); */
/* 	free(s); */
/* } */

/* static int sum(const int *procs, const int n) */
/* { */
/* 	int ret = 0; */
/* 	for (int i = 0; i < n; i++) { */
/* 		ret += procs[i]; */
/* 	} */
/* 	return ret; */
/* } */

/* static int threshold(const int sum, const int n, const int rank) */
/* { */
/* 	return sum / n + (n - rank <= sum % n); */
/* } */

/* static int bucket_comp_rank(const void *_a, const void *_b) */
/* { */
/* 	struct bucket_t a = *((struct bucket_t *)_a); */
/* 	struct bucket_t b = *((struct bucket_t *)_b); */
/* 	return a.rank - b.rank; */
/* } */

/* static int bucket_comp(const void *_a, const void *_b) */
/* { */
/* 	struct bucket_t a = *((struct bucket_t *)_a); */
/* 	struct bucket_t b = *((struct bucket_t *)_b); */
/* 	int count_a = bucket_count(&a); */
/* 	int count_b = bucket_count(&b); */
/* 	if (count_a == count_b) { */
/* 		return bucket_comp_rank(_a, _b); */
/* 	} */
/* 	return count_a - count_b; */
/* } */

/* static void snapshot_print(snapshot_t *s) */
/* { */
/* 	printf("*************\n"); */
/* 	for (int i = 0; i < s->n; i++) { */
/* 		printf("%d: ", s->procs[i].rank); */
/* 		for (int j = s->procs[i].start; j < s->procs[i].len; j++) { */
/* 			printf("(%d %d) ", s->procs[i].tasks[j].al1, s->procs[i].tasks[j].al2); */
/* 		} */
/* 		printf("\n"); */
/* 	} */
/* 	printf("*************\n"); */
/* } */

/* diff_t snapshot_update(snapshot_t *s, const int *procs, const int rank) */
/* { */
/* 	diff_t diff = { */
/* 	    .empty = 0, */
/* 	    .deleted = 0, */
/* 	    .appended = 0, */
/* 	    .tasks = NULL, */
/* 	}; */

/* 	int su = sum(procs, s->n); */
/* 	diff.empty = su == 0; */
/* 	qsort(s->procs, s->n, sizeof(struct bucket_t), bucket_comp_rank); */

/* 	for (int i = 0; i < s->n; i++) { */
/* 		bucket_leave_last(s->procs + i, procs[i]); */
/* 	} */

/* 	qsort(s->procs, s->n, sizeof(struct bucket_t), bucket_comp); */

/* 	int i = 0; */
/* 	int j = s->n - 1; */
/* 	for (;;) { */
/* 		while (i < j && bucket_count(s->procs + i) >= threshold(su, s->n, i)) { */
/* 			i++; */
/* 		} */
/* 		while (i < j && bucket_count(s->procs + j) <= threshold(su, s->n, j)) { */
/* 			j--; */
/* 		} */
/* 		if (i >= j) { */
/* 			break; */
/* 		} */
/* 		mpi_task_t t = bucket_pop(s->procs + j); */
/* 		bucket_push(s->procs + i, t); */
/* 		if (s->procs[j].rank == rank) { */
/* 			diff.deleted++; */
/* 		} */
/* 		else if (s->procs[i].rank == rank) { */
/* 			if (diff.tasks == NULL) { */
/* 				diff.tasks = calloc(su / s->n + 1, sizeof(mpi_task_t)); */
/* 			} */
/* 			diff.tasks[diff.appended++] = t; */
/* 		} */
/* 	} */
/* 	return diff; */
/* } */

/* int main() */
/* { */
/* 	snapshot_t *s = snapshot_init(2); */
/* 	snapshot_fill_pairwise(s, 3); */
/* 	snapshot_print(s); */
/* 	int procs[] = {2, 1}; */
/* 	diff_t diff = snapshot_update(s, procs, 0); */
/* 	procs[0] = 1; */
/* 	procs[1] = 0; */
/* 	snapshot_print(s); */
/* 	diff = snapshot_update(s, procs, 0); */
/* 	snapshot_print(s); */
/* 	procs[0] = 1; */
/* 	procs[1] = 0; */
/* 	diff = snapshot_update(s, procs, 0); */
/* 	/\* printf("deleted %d\n", diff.deleted); *\/ */
/* 	/\* for (int i = 0; i < diff.appended; i++) { *\/ */
/* 	/\* 	printf("appended (%d %d)\n", diff.tasks[i].al1, diff.tasks[i].al2); *\/ */
/* 	/\* } *\/ */
/* 	/\* diff_destroy(&diff); *\/ */
/* 	snapshot_print(s); */
/* 	snapshot_destroy(s); */
/* 	return 0; */
/* } */

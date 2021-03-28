#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ga.h"

#define CFG_FITNESS_FUNC "fitness="
#define CFG_FITNESS_FUNC_MATT "matt"
#define CFG_FITNESS_FUNC_TM "tm"

#define CFG_POP_SIZE "population_size="
#define CFG_TERMINATE_ELAPSED "terminate_elapsed="
#define CFG_TERMINATE_UNCHANGED "terminate_unchanged="
#define CFG_TERMINATE_ITERATION "terminate_iteration="

static int parse_score(char *s){
	if (strncmp(CFG_FITNESS_FUNC_MATT, s, strlen(CFG_FITNESS_FUNC_MATT)) == 0) {
		return SCORE_MATT;
	}
	if (strncmp(CFG_FITNESS_FUNC_TM, s, strlen(CFG_FITNESS_FUNC_TM)) == 0) {
		return SCORE_TM;
	}
	return -1;
}

static int parse_int(char *s)
{
	long int ret = strtol(s, NULL, 10);
	if (errno == EINVAL || errno == ERANGE) {
		return -1;
	}
	return ret;
}

static double parse_double(char *s)
{
	double ret = strtod(s, NULL);
	if (errno == EINVAL || errno == ERANGE) {
		return -1;
	}
	return ret;
}

int check_cfg(struct ga_cfg_t *cfg, const int terminate_criteria)
{
	if (cfg->population_size <= 0) {
		printf("cfg: 'population_size' should be positive\n");
		return ERR_CFG;
	}

	if (terminate_criteria == 0) {
		return 0;
	}

	if (cfg->terminate_on_elapsed <= 0) {
		printf("cfg: 'terminate_elapsed' should be positive\n");
		return ERR_CFG;
	}

	if (cfg->terminate_on_unchanged <= 0) {
		printf("cfg: 'terminate_unchanged' should be positive\n");
		return ERR_CFG;
	}

	if (cfg->terminate_on_iteration <= 0) {
		printf("cfg: 'terminate_iteration' should be positive\n");
		return ERR_CFG;
	}

	cfg->terminate_criteria = terminate_criteria;
	return 0;
}

static double ga_cfg_population_size(const int n_als)
{
	if (n_als > 300) {
		return 2 * n_als;
	}

	double x1 = 1;
	double y1 = 10;
	double x2 = 300;
	double y2 = 2;
	double y = ((y2 - y1) * n_als + x2 * y1 - x1 * y2) / (x2 - x1);
	
	return ceil(y * n_als);
}

static struct ga_cfg_t ga_default_cfg(const int n_als, const int nthreads)
{
	struct ga_cfg_t cfg = {
	    .enabled = 0,
	    .nthreads = nthreads,
	    .population_size = ga_cfg_population_size(n_als),
	    .tree_size = n_als,

	    .terminate_criteria = GA_TERMINATE_ON_ITERATION | GA_TERMINATE_ON_ELAPSED | GA_TERMINATE_ON_UNCHANGED,
	    .terminate_on_iteration = n_als * 10,
	    .terminate_on_elapsed = 7200,
	    .terminate_on_unchanged = 40,
		.fitness_func = SCORE_MATT,
	};
	return cfg;
}


int ga_parse_cfg(struct ga_cfg_t *cfg, const char *path, const int n_als, const int nthreads)
{
	*cfg = ga_default_cfg(n_als, nthreads);
	FILE *stream;
	char *line = NULL;
	size_t len = 0;
	ssize_t nread;

	stream = fopen(path, "r");
	if (stream == NULL) {
		return ERR_FILE;
	}

	int terminate_criteria = 0;

	while ((nread = getline(&line, &len, stream)) != -1) {
		if (line[0] == '#' || line[0] == '\n') {
			continue;
		}
		char *p = line;
		for (; *p != '\n'; p++)
			;
		*p = '\0';

		int cmp;
		int d = -2;
		double v = -2;

		cmp = strncmp(CFG_POP_SIZE, line, strlen(CFG_POP_SIZE));
		if (cmp == 0) {
			d = parse_int(line + strlen(CFG_POP_SIZE));
			cfg->population_size = d;
			goto CHECK;
		}

		cmp = strncmp(CFG_TERMINATE_ELAPSED, line, strlen(CFG_TERMINATE_ELAPSED));
		if (cmp == 0) {
			d = parse_int(line + strlen(CFG_TERMINATE_ELAPSED));
			cfg->terminate_on_elapsed = d;
			terminate_criteria |= GA_TERMINATE_ON_ELAPSED;
			goto CHECK;
		}

		cmp = strncmp(CFG_TERMINATE_UNCHANGED, line, strlen(CFG_TERMINATE_UNCHANGED));
		if (cmp == 0) {
			d = parse_int(line + strlen(CFG_TERMINATE_UNCHANGED));
			cfg->terminate_on_unchanged = d;
			terminate_criteria |= GA_TERMINATE_ON_UNCHANGED;
			goto CHECK;
		}

		cmp = strncmp(CFG_TERMINATE_ITERATION, line, strlen(CFG_TERMINATE_ITERATION));
		if (cmp == 0) {
			d = parse_int(line + strlen(CFG_TERMINATE_ITERATION));
			cfg->terminate_on_iteration = d;
			terminate_criteria |= GA_TERMINATE_ON_ITERATION;
			goto CHECK;
		}

		cmp = strncmp(CFG_FITNESS_FUNC, line, strlen(CFG_FITNESS_FUNC));
		if (cmp == 0) {
			d = parse_score(line + strlen(CFG_FITNESS_FUNC));
			cfg->fitness_func = d;
			goto CHECK;
		}

	CHECK:
		if ((d == -2 && v == -2) || (d == -1 || v == -1)) {
			printf("cfg: could not parse line: '%s'\n", line);
			if (line != NULL) {
				free(line);
			}
			fclose(stream);
			return ERR_CFG;
		}
	}

	if (line != NULL) {
		free(line);
	}
	fclose(stream);
	return check_cfg(cfg, terminate_criteria);
}

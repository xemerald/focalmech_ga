/**
 * @file ga_core.c
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-08-26
 *
 * @copyright Copyright (c) 2023
 *
 */
/* */
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
/* */
#include <ga_core.h>

/* */
static void crossover_gene( uint32_t, uint32_t, uint32_t *, uint32_t * );
static uint32_t mutate_gene( uint32_t, int );
static int compare_gscore( const void *, const void * );

/* */
#define GA_GEN_PAIR_POPINDEX(_INDEX_A, _INDEX_B, _RANGE) \
		__extension__({ \
			(_INDEX_A) = rand() % (_RANGE); \
			do { \
				(_INDEX_B) = rand() % (_RANGE); \
			} while ( (_INDEX_A) == (_INDEX_B) ); \
		})

/**
 * @brief
 *
 * @param pop
 * @param npop
 * @param score_func
 * @param observes
 * @param nobs
 * @return GA_POPULATION*
 */
GA_POPULATION *ga_present_gen_fast( GA_POPULATION *pop, int npop, double (*score_func)(uint32_t, void *, int), void *observes, int nobs )
{
	srand(time(NULL));

	for ( int i = 0; i < npop; i++ ) {
		pop[i].gene  = rand() & GA_CORE_GENE_FULL_MASK;
		pop[i].score = score_func( pop[i].gene, observes, nobs );
	}
/* */
	qsort(pop, npop, sizeof(GA_POPULATION), compare_gscore);

	return pop;
}

/**
 * @brief
 *
 * @param ppop
 * @param npop
 * @param mutate_bits
 * @param repro_rate
 * @param mutate_rate
 * @param score_func
 * @param observes
 * @param nobs
 */
void ga_one_generation( GA_POPULATION *ppop, int npop, int mutate_bits, double repro_rate, double mutate_rate, double (*score_func)(uint32_t, void *, int), void *observes, int nobs )
{
	int       ind_a, ind_b;
	uint32_t  fgene_a = 0, fgene_b = 0;
	const int i_repro = (int)(npop * repro_rate);
	const int i_mutate = (int)(npop * (mutate_rate + repro_rate));
	const int halfpop = (int)(npop * 0.5);
/* Next population */
	GA_POPULATION fpop[npop];

/* */
	for ( int i = 0; i < npop; i++ ) {
	/* Reproduction prefect or with in top (reproduction*100) % */
		if ( i < i_repro ) {
			fpop[i] = ppop[i];
			continue;
		}
	/* Do simple crossover from top 50% of presents */
		if ( !fgene_a && !fgene_b ) {
			GA_GEN_PAIR_POPINDEX( ind_a, ind_b, halfpop );
			crossover_gene( ppop[ind_a].gene, ppop[ind_b].gene, &fgene_a, &fgene_b );
		}
		else if ( !fgene_b ) {
			fpop[i].gene = fgene_a;
			fgene_a = 0;
		}
		else {
			fpop[i].gene = fgene_b;
			fgene_b = 0;
		}
	/* Do mutation */
		if ( i < i_mutate )
			fpop[i].gene = mutate_gene( fpop[i].gene, mutate_bits );
	/* */
		fpop[i].score = score_func( fpop[i].gene, observes, nobs );
	}
/* */
	memcpy(ppop, fpop, sizeof(GA_POPULATION) * npop);
	qsort(ppop, npop, sizeof(GA_POPULATION), compare_gscore);

	return;
}

/**
 * @brief
 *
 * @param pgene_1
 * @param pgene_2
 * @param fgene_1
 * @param fgene_2
 */
static void crossover_gene( uint32_t pgene_1, uint32_t pgene_2, uint32_t *fgene_1, uint32_t *fgene_2 )
{
/* 4-5-3-4-4-5 */
	*fgene_1 = (pgene_2 & GA_CORE_CROSSOVER_MASK_LOW) | (pgene_1 & GA_CORE_CROSSOVER_MASK_UP);
	*fgene_2 = (pgene_1 & GA_CORE_CROSSOVER_MASK_LOW) | (pgene_2 & GA_CORE_CROSSOVER_MASK_UP);

	return;
}

/**
 * @brief
 *
 * @param gene
 * @param nbits
 * @return uint32_t
 */
static uint32_t mutate_gene( uint32_t gene, int nbits )
{
	uint32_t mask = 0;

	for ( int i = 0; i < nbits; i++ )
		mask |= 1 << (rand() % GA_CORE_GENE_TOTAL_BITS);

	gene ^= mask;

	return gene;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return int
 */
static int compare_gscore( const void *a, const void *b )
{
	double a_score = ((GA_POPULATION *)a)->score;
	double b_score = ((GA_POPULATION *)b)->score;

	if ( a_score < b_score )
		return 1;
	else if ( a_score > b_score )
		return -1;
	else
		return 0;
}

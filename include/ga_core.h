/**
 * @file ga_core.h
 * @author Benjamin Yang (b98204032@gmail.com)
 * @author Department of Geology in National Taiwan University
 * @brief
 * @version 0.1
 * @date 2023-08-26
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once
/* */
#include <stdint.h>
/* */
#define GA_CORE_GENE_TOTAL_BITS    25
#define GA_CORE_GENE_FULL_MASK     0b1111111111111111111111111
#define GA_CORE_CROSSOVER_MASK_LOW 0b0000111110001111000011111  /* 0-5-0-4-0-5 */
#define GA_CORE_CROSSOVER_MASK_UP  0b1111000001110000111100000  /* 4-0-3-0-4-0 */
/* */
#define GA_CORE_STRIKE_BITS  9
#define GA_CORE_DIP_BITS     7
#define GA_CORE_RAKE_BITS    9
#define GA_CORE_STRIKE_MASK  0b111111111
#define GA_CORE_DIP_MASK     0b1111111
#define GA_CORE_RAKE_MASK    0b111111111

/**
 * @brief
 *
 */
typedef struct {
	uint32_t gene;
	double   score;
} GA_POPULATION;

/**
 * @brief
 *
 * @return GA_POPULATION*
 */
GA_POPULATION *ga_present_gen_fast( GA_POPULATION *, int, double (*)(uint32_t, void *, int), void *, int );
/**
 * @brief
 *
 */
void ga_one_generation( GA_POPULATION *, int, int, double, double, double (*)(uint32_t, void *, int), void *, int );

/**
 * @file tko_azi_cal.c
 * @anchor Benjamin Ming Yang
 * @brief
 * @version 0.1
 * @date 2023-09-25
 *
 * @copyright Copyright (c) 2023
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/* */
#include <raytracing.h>

/**
 * @brief
 *
 * @param tko
 * @param azi
 * @param evla
 * @param evlo
 * @param evdp
 * @param stla
 * @param stlo
 * @param stdp
 */
int tac_main( double *tko, double *azi, double evla, double evlo, double evdp, double stla, double stlo, double stdp )
{
	RAY_INFO ray[RT_MAX_NODE + 1];
	int np;
	double tt;

/* */
	if ( !rt_main( ray, &np, &tt, evla, evlo, evdp, stla, stlo, stdp, RT_P_WAVE_VELOCITY ) ) {
		rt_tko_azi_cal( ray, np, tko, azi );
	}
	else {
		*tko = 0.0;
		*azi = 0.0;
	}

	return 0;
}

/**
 * @brief
 *
 * @param model_path
 * @return int
 */
int tac_velmod_load( const char *model_path )
{
	return rt_velmod_load( model_path );
}

/**
 * @brief
 *
 */
void tac_velmod_free( void )
{
	rt_velmod_free();
	return;
}

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <constants.h>
#include <ga_core.h>
#include <fpl_func.h>

/* */
static uint32_t encode_gene_ga( const double, const double, const double );
static void decode_gene_ga( double *, double *, double *, uint32_t );
static double cal_gene_score_ga( const uint32_t, void *, const int );
static double cal_amp( double, double, double, double, double );
static int remove_approx_result( FPL_RESULT *, int );
static int is_same_group_result( FPL_RESULT *, const FPL_RESULT *, const double );
static int compare_azimuth( const void *, const void * );

/**
 * @brief
 *
 * @param observes
 * @param nobs
 * @param nitr
 * @param npop
 * @param mutate_bits
 * @param repro_rate
 * @param mutate_rate
 * @param output
 * @param nout
 * @return double
 */
double fpl_find_ga( void *observes, int nobs, int nitr, int npop, int mutate_bits, double repro_rate, double mutate_rate, FPL_RESULT *output, int *nout )
{
	GA_POPULATION ppop[npop];
	const int     quarter_pop = (int)(npop * 0.25);
	double        score_sum;

/* */
	ga_present_gen_fast( ppop, npop, cal_gene_score_ga, observes, nobs );
/* */
	for ( int i = 0; i < nitr; i++ ) {
		ga_one_generation( ppop, npop, mutate_bits, repro_rate, mutate_rate, cal_gene_score_ga, observes, nobs );
	/* */
		score_sum = 0.0;
		for ( int k = 0; k < quarter_pop; k++ )
			score_sum += ppop[k].score;
	/* */
		if ( (score_sum /= (double)quarter_pop) > 0.999 )
			break;
	}
/* */
	int _nout = 0;
	for ( int i = 0; i < npop; i++ ) {
		if ( fabs(ppop[i].score - ppop[0].score) < DBL_EPSILON ) {
			decode_gene_ga( &output[_nout].strike, &output[_nout].dip, &output[_nout].rake, ppop[i].gene );
			_nout++;
		}
	}

	*nout = _nout;
	return ppop[0].score;
}

/**
 * @brief
 *
 * @param results
 * @param sdv
 * @param nresult
 * @return int
 */
int fpl_result_refine( FPL_RESULT *results, FPL_RESULT *sdv, int nresult )
{
	double     tmp;
	int        nsolution;
	int        _nresult;
	FPL_RESULT _results[nresult];
	FPL_RESULT solutions[nresult];
	FPL_RESULT _sdv[nresult];
	int        resflag[nresult];
	FPL_RESULT dbresult[2];
/* */
	struct RHist {
		int count;
		double strike;
		double dip;
		double rake;
		double best_diff;
		FPL_RESULT *best_ptr;
	} rhist[9][9][9];

	struct RHist *max_rhist = &rhist[0][0][0];

/* */
	memset(rhist, 0, sizeof(rhist));
	for ( int i = 0; i < nresult; i++ ) {
		if ( results[i].strike > FOCAL_GA_PI2 )
			results[i].strike -= FOCAL_GA_PI2;
		else if ( results[i].strike < 0.0 )
			results[i].strike += FOCAL_GA_PI2;
	/* */
		fpl_dbcouple( &results[i], dbresult, NULL, NULL );
		for ( int j = 0; j < 2; j++ ) {
		/* */
			if ( dbresult[j].strike > FOCAL_GA_PI2 )
				dbresult[j].strike -= FOCAL_GA_PI2;
			else if ( dbresult[j].strike < 0.0 )
				dbresult[j].strike += FOCAL_GA_PI2;
		/* */
			if ( dbresult[j].rake < 0.0)
				dbresult[j].rake += FOCAL_GA_PI2;
		}
	/* */
		if ( dbresult[0].strike < dbresult[1].strike )
			_results[i] = dbresult[0];
		else
			_results[i] = dbresult[1];
	/* */
		resflag[i] = 1;
	}
/* */
	nresult = remove_approx_result( _results, nresult );
/* */
	for ( int i = 0; i < nresult; i++ ) {
		int istr  = (int)(_results[i].strike / FOCAL_GA_PI2 * 9.0);
		int idip  = (int)(_results[i].dip / FOCAL_GA_HALF_PI * 9.0);
		int irake = (int)(_results[i].rake / FOCAL_GA_PI2 * 9.0);
	/* */
		if ( !rhist[istr][idip][irake].count ) {
			rhist[istr][idip][irake].strike = istr * FOCAL_GA_PI2 / 9.0 + FOCAL_GA_PI / 9.0;
			rhist[istr][idip][irake].dip = idip * FOCAL_GA_HALF_PI / 9.0 + FOCAL_GA_HALF_PI / 18.0;
			rhist[istr][idip][irake].rake = irake * FOCAL_GA_PI2 / 9.0 + FOCAL_GA_PI / 9.0;
		/* */
			rhist[istr][idip][irake].best_diff = fabs(_results[i].strike - rhist[istr][idip][irake].strike) + fabs(_results[i].dip - rhist[istr][idip][irake].dip) + fabs(_results[i].rake - rhist[istr][idip][irake].rake);
			rhist[istr][idip][irake].best_ptr = &_results[i];
		}
		else if ( (fabs(_results[i].strike - rhist[istr][idip][irake].strike) + fabs(_results[i].dip - rhist[istr][idip][irake].dip) + fabs(_results[i].rake - rhist[istr][idip][irake].rake)) < rhist[istr][idip][irake].best_diff ) {
			rhist[istr][idip][irake].best_diff = fabs(_results[i].strike - rhist[istr][idip][irake].strike) + fabs(_results[i].dip - rhist[istr][idip][irake].dip) + fabs(_results[i].rake - rhist[istr][idip][irake].rake);
			rhist[istr][idip][irake].best_ptr = &_results[i];
		}
	/* */
		if ( ++rhist[istr][idip][irake].count > max_rhist->count )
			max_rhist = &rhist[istr][idip][irake];
	}
/* */
	solutions[0] = *(max_rhist->best_ptr);
	nsolution = 1;
	for ( int i = 0; i < nresult; i++ ) {
		resflag[i] = -1;
		fpl_dbcouple( &_results[i], dbresult, NULL, NULL );
	/* */
		for ( int j = 0; j < 2; j++ ) {
		/* */
			if ( dbresult[j].strike > FOCAL_GA_PI2 )
				dbresult[j].strike -= FOCAL_GA_PI2;
			else if ( dbresult[j].strike < 0.0 )
				dbresult[j].strike += FOCAL_GA_PI2;
		/* */
			if ( dbresult[j].rake < 0.0)
				dbresult[j].rake += FOCAL_GA_PI2;
			else if ( dbresult[j].rake > FOCAL_GA_PI2 )
				dbresult[j].rake -= FOCAL_GA_PI2;
		/* */
			for ( int k = 0; k < nsolution; k++ ) {
				if ( is_same_group_result( &dbresult[j], &solutions[k], 90.0 * FOCAL_GA_DEG2RAD ) ) {
					resflag[i] = k;
					_results[i] = dbresult[j];
					break;
				}
			/* */
				if ( dbresult[j].strike < FOCAL_GA_PI )
					dbresult[j].strike += FOCAL_GA_PI;
				else
					dbresult[j].strike -= FOCAL_GA_PI;

				if ( dbresult[j].dip > FOCAL_GA_HALF_PI * 0.5 ) {
					dbresult[j].dip  = FOCAL_GA_PI - dbresult[j].dip;
					dbresult[j].rake = -1.0 * dbresult[j].rake;
				}
				else {
					dbresult[j].dip  = -1.0 * dbresult[j].dip;
					dbresult[j].rake = FOCAL_GA_PI + dbresult[j].rake;
				}

				if ( dbresult[j].rake < 0.0 )
					dbresult[j].rake += FOCAL_GA_PI2;
				if ( dbresult[j].rake > FOCAL_GA_PI2 )
					dbresult[j].rake -= FOCAL_GA_PI2;
				if ( is_same_group_result( &dbresult[j], &solutions[k], 90.0 * FOCAL_GA_DEG2RAD ) ) {
					resflag[i] = k;
					_results[i] = dbresult[j];
					break;
				}
			}
			if ( resflag[i] >= 0 )
				break;
		}
	/* */
		if ( resflag[i] < 0 ) {
			solutions[nsolution] = _results[i];
			resflag[i] = nsolution;
			nsolution++;
		}
	/* */
		if ( nsolution >= 100 )
			break;
	}

/* */
	_nresult = 0;
	for ( int i = 0; i < nsolution; i++ ) {
		double _strike = 0.0;
		double _dip = 0.0;
		double _rake = 0.0;
		int    nn = 0;

		for ( int j = 0; j < nresult; j++ ) {
			if ( resflag[j] == i ) {
				nn++;
				_strike += _results[i].strike;
				_dip += _results[i].dip;
				_rake += _results[i].rake;
			}
		}

		if ( nn ) {
			_strike /= nn;
			_dip /= nn;
			_rake /= nn;
		}
		else {
			continue;
		}

		tmp = 999999999.9;
		for ( int j = 0; j < nresult; j++ ) {
			if ( resflag[j] == i ) {
				if( (fabs(_strike - _results[j].strike) + fabs(_dip - _results[j].dip) + fabs(_rake - _results[j].rake)) < tmp ) {
					tmp = fabs(_strike - _results[j].strike) + fabs(_dip - _results[j].dip) + fabs(_rake - _results[j].rake);
					solutions[i] = _results[j];
				}
			}
		}

		_sdv[i].strike = 0.0;
		_sdv[i].dip = 0.0;
		_sdv[i].rake = 0.0;
		for ( int j = 0; j < nresult; j++ ) {
			if ( resflag[j] == i ) {
				_sdv[i].strike += (_results[j].strike - _strike) * (_results[j].strike - _strike);
				_sdv[i].dip += (_results[j].dip - _dip) * (_results[j].dip - _dip);
				_sdv[i].rake += (_results[j].rake - _rake) * (_results[j].rake - _rake);
			}
		}
		if ( nn > 1 ) {
			_sdv[i].strike = 2.0 * sqrt(_sdv[i].strike / nn);
			_sdv[i].dip = 2.0 * sqrt(_sdv[i].dip / nn);
			_sdv[i].rake = 2.0 * sqrt(_sdv[i].rake / nn);
		}
		else {
			_sdv[i].strike = 1.4 * FOCAL_GA_DEG2RAD;
			_sdv[i].dip = 1.4 * FOCAL_GA_DEG2RAD;
			_sdv[i].rake = 1.4 * FOCAL_GA_DEG2RAD;
		}
	/* */
		if ( _sdv[i].strike < FOCAL_GA_DEG2RAD )
			_sdv[i].strike = FOCAL_GA_DEG2RAD;
		if ( _sdv[i].dip < FOCAL_GA_DEG2RAD )
			_sdv[i].dip = FOCAL_GA_DEG2RAD;
		if ( _sdv[i].rake < FOCAL_GA_DEG2RAD )
			_sdv[i].rake = FOCAL_GA_DEG2RAD;
	/* */
		if ( solutions[i].dip > FOCAL_GA_HALF_PI ) {
			solutions[i].strike = solutions[i].strike + FOCAL_GA_PI;
			solutions[i].dip = FOCAL_GA_PI - solutions[i].dip;
			solutions[i].rake = FOCAL_GA_PI2 - solutions[i].rake;
		}
		else if ( solutions[i].dip < 0.0 ) {
			solutions[i].strike = solutions[i].strike + FOCAL_GA_PI;
			solutions[i].dip = -1.0 * solutions[i].dip;
			solutions[i].rake = FOCAL_GA_PI2 - solutions[i].rake;
		}

		if ( solutions[i].strike > FOCAL_GA_PI2)
			solutions[i].strike -= FOCAL_GA_PI2;
		if ( solutions[i].rake > FOCAL_GA_PI )
			solutions[i].rake -= FOCAL_GA_PI2;
		if ( solutions[i].rake < -FOCAL_GA_PI )
			solutions[i].rake += FOCAL_GA_PI2;
	/* Passing the solutions */
		results[_nresult] = solutions[i];
		sdv[_nresult] = _sdv[i];
		_nresult++;
	}

	return _nresult;
}

/**
 * @brief
 *
 * @param observes
 * @param nobs
 * @param f_score
 * @return double
 */
double fpl_quality_cal( FPL_OBSERVE *observes, const int nobs, const double f_score )
{
	int    npos_por;
	double quality;
	double quality_gap;
	double quality_nobs;
	double quality_score;
	double quality_por;

/* Find Gap Quality */
	qsort(observes, nobs, sizeof(FPL_OBSERVE), compare_azimuth);
	double gap = observes[0].azimuth + FOCAL_GA_PI2 - observes[nobs - 1].azimuth;
	for ( int i = 1; i < nobs; i++ ) {
		if( (observes[i].azimuth - observes[i - 1].azimuth) > gap )
			gap = observes[i].azimuth - observes[i - 1].azimuth;
	}
	if ( gap >= FOCAL_GA_PI )
		quality_gap = 0.0;
	else
		quality_gap = (FOCAL_GA_PI - gap) / FOCAL_GA_HALF_PI;

/* Find # Quality */
	if ( nobs < 10 )
		quality_nobs = 0.0;
	else if ( nobs < 50 )
		quality_nobs = (nobs - 10.0) / 20.0;
	else
		quality_nobs = 2.0;

/* Find F-score Quality */
	if ( f_score < 0.7 )
		quality_score = 0.0;
	else
		quality_score = (f_score - 0.7) / 0.15;

/* Find Polarity Quality */
	npos_por = 0;
	for ( int i = 0; i < nobs; i++ ){
		if (observes[i].polarity > 0.5 )
			npos_por++;
	}
	quality_por = (0.5 - fabs(npos_por * 1.0 / nobs - 0.5)) / 0.25;

	quality = quality_gap * quality_nobs * quality_score * quality_por;
	if ( quality < 0.001 )
		return 0.0;

	return quality;
}

/**
 * @brief
 *
 * @param plane
 * @param plane_couple
 * @param axis_t
 * @param axis_p
 * @return FPL_RESULT*
 */
FPL_RESULT *fpl_dbcouple( FPL_RESULT *plane, FPL_RESULT plane_couple[2], FPL_RESULT *axis_t, FPL_RESULT *axis_p )
{
	double fp11 = cos(plane->strike);
	double fp12 = sin(plane->strike);
	double fp13 = cos(plane->dip);
	double fp23 = sin(plane->rake);
	double val  = cos(plane->rake);
	double fp21 = val * fp12 - fp13 * fp23 * fp11;
	double fp22 = val * fp11 + fp13 * fp23 * fp12;
	double x, y, z;

	val   = sin(plane->dip);
	fp23 *= val;
	fp11 *= val;
	fp12 *= -val;

/* First solution of the double couple plane */
	plane_couple[0] = *plane;
/* */
	if ( fp23 < 0.0 ) {
		x = -fp21;
		y = -fp22;
		z = -fp23;
	}
	else {
		x = fp21;
		y = fp22;
		z = fp23;
	}
	plane_couple[1].dip = acos(z);

	if ( z < 1.0 ) {
		val = sin(plane_couple[1].dip);
		x /= val;
		y /= val;
	}
	z = -1.0 * y;
	y = x;
	x = z;
	if ( y > 1.0 )
		y = 1.0;
	else if ( y < -1.0 )
		y = -1.0;
	val = fabs(acos(y));
	if ( x >= 0.0 )
		plane_couple[1].strike = val;
	else
		plane_couple[1].strike = -val;
	if ( plane_couple[1].strike < 0.0 )
		plane_couple[1].strike += FOCAL_GA_PI2;
	else if ( plane_couple[1].strike >= FOCAL_GA_PI2 )
		plane_couple[1].strike -= FOCAL_GA_PI2;

	if ( fp23 < 0.0 ) {
		x = -x;
		y = -y;
	}
	val = x * fp11 + y * fp12;
	if ( val > 1.0 )
		val = 1.0;
	else if ( val < -1.0 )
		val = -1.0;
	val = fabs(acos(val));
	if ( fp23 >= 0.0 )
		plane_couple[1].rake = val;
	else
		plane_couple[1].rake = -val;

/* Calculation of T-axis */
	if ( axis_t ) {
	/* plunge => dip; azimuth => strike */
		x = fp11 + fp21;
		y = fp12 + fp22;
		z = fp13 + fp23;

		val = sqrt(x * x + y * y + z * z);
		x /= val;
		y /= val;
		z /= val;
		if ( z > 0.0 ) {
			x = -x;
			y = -y;
			z = -z;
		}
		axis_t->dip = asin(fabs(z));

		if ( z != 0.0 ) {
			val = cos(axis_t->dip);
			x /= val;
			y /= val;
		}
		if ( y > 1.0 )
			y = 1.0;
		else if ( y < -1.0 )
			y = -1.0;
		val = fabs(acos(y));
		if ( x >= 0.0 )
			axis_t->strike = val;
		else
			axis_t->strike = -val;

		if ( axis_t->strike < 0.0 )
			axis_t->strike += FOCAL_GA_PI2;
		else if ( axis_t->strike >= FOCAL_GA_PI2 )
			axis_t->strike -= FOCAL_GA_PI2;
	}

/* Calculation of P-axis */
	if ( axis_p ) {
		x = fp11 - fp21;
		y = fp12 - fp22;
		z = fp13 - fp23;

		val = sqrt(x * x + y * y + z * z);
		x /= val;
		y /= val;
		z /= val;
		if ( z > 0.0 ) {
			x = -x;
			y = -y;
			z = -z;
		}
		axis_p->dip = asin(fabs(z));

		if ( z != 0.0 ) {
			val = cos(axis_p->dip);
			x /= val;
			y /= val;
		}
		if ( y > 1.0 )
			y = 1.0;
		else if ( y < -1.0 )
			y = -1.0;
		val = fabs(acos(y));
		if ( x >= 0.0 )
			axis_p->strike = val;
		else
			axis_p->strike = -val;
		if ( axis_p->strike < 0.0 )
			axis_p->strike += FOCAL_GA_PI2;
		else if ( axis_p->strike >= FOCAL_GA_PI2 )
			axis_p->strike -= FOCAL_GA_PI2;
	}

	return plane_couple;
}

/**
 * @brief
 *
 * @param input
 * @return FPL_OBSERVE
 */
FPL_OBSERVE fpl_observe_deg2rad( FPL_OBSERVE *input )
{
	return (FPL_OBSERVE){
		.azimuth = input->azimuth * FOCAL_GA_DEG2RAD,
		.takeoff = input->takeoff * FOCAL_GA_DEG2RAD,
		.polarity = input->polarity
	};
}

/**
 * @brief
 *
 * @param input
 * @return FPL_OBSERVE
 */
FPL_OBSERVE fpl_observe_rad2deg( FPL_OBSERVE *input )
{
	return (FPL_OBSERVE){
		.azimuth = input->azimuth * FOCAL_GA_RAD2DEG,
		.takeoff = input->takeoff * FOCAL_GA_RAD2DEG,
		.polarity = input->polarity
	};
}

/**
 * @brief
 *
 * @param input
 * @return FPL_RESULT
 */
FPL_RESULT fpl_result_deg2rad( FPL_RESULT *input )
{
	return (FPL_RESULT){
		.dip = input->dip * FOCAL_GA_DEG2RAD,
		.strike = input->strike * FOCAL_GA_DEG2RAD,
		.rake = input->rake * FOCAL_GA_DEG2RAD
	};
}

/**
 * @brief
 *
 * @param input
 * @return FPL_RESULT
 */
FPL_RESULT fpl_result_rad2deg( FPL_RESULT *input )
{
	return (FPL_RESULT){
		.dip = input->dip * FOCAL_GA_RAD2DEG,
		.strike = input->strike * FOCAL_GA_RAD2DEG,
		.rake = input->rake * FOCAL_GA_RAD2DEG
	};
}

/**
 * @brief
 *
 * @param strike
 * @param dip
 * @param rake
 * @return uint32_t
 */
static uint32_t encode_gene_ga( const double strike, const double dip, const double rake )
{
	uint32_t result;
	uint32_t istrike = (uint32_t)(strike / FOCAL_GA_PI2 * 512.0);
	uint32_t idip = (uint32_t)(dip / FOCAL_GA_HALF_PI * 128.0);
	uint32_t irake = (uint32_t)(rake / FOCAL_GA_PI2 * 512.0);

	result = irake & GA_CORE_RAKE_MASK;
	result <<= GA_CORE_RAKE_BITS;
	result |= idip & GA_CORE_DIP_MASK;
	result <<= GA_CORE_DIP_BITS;
	result |= istrike & GA_CORE_STRIKE_MASK;

	return result;
}

/**
 * @brief
 *
 * @param strike
 * @param dip
 * @param rake
 * @param gene
 */
static void decode_gene_ga( double *strike, double *dip, double *rake, uint32_t gene )
{
	uint32_t istrike;
	uint32_t idip;
	uint32_t irake;

	istrike = gene & GA_CORE_STRIKE_MASK;
	gene >>= GA_CORE_STRIKE_BITS;
	idip = gene & GA_CORE_DIP_MASK;
	gene >>= GA_CORE_DIP_BITS;
	irake = gene & GA_CORE_RAKE_MASK;

	*strike = istrike * FPLF_DECODE_RAD_CONST;
	*dip = idip * FPLF_DECODE_RAD_CONST;
	*rake = irake * FPLF_DECODE_RAD_CONST;

	return;
}

/**
 * @brief
 *
 * @param gene
 * @param observes
 * @param nobs
 * @return double
 */
static double cal_gene_score_ga( const uint32_t gene, void *observes, const int nobs )
{
	FPL_OBSERVE *_observe;
	double strike, dip, rake;
	int n = 0;

	decode_gene_ga( &strike, &dip, &rake, gene );

	_observe = observes;
	for ( int i = 0; i < nobs; i++, _observe++ ) {
		if ( (cal_amp( strike, dip, rake, _observe->azimuth, _observe->takeoff ) * _observe->polarity) > 0.0 )
			n++;
		else
			n--;
	}

	return 1.0 * n / nobs;
}

/**
 * @brief
 *
 * @param strike
 * @param dip
 * @param slip
 * @param azimuth
 * @param takeoff
 * @return double
 */
static double cal_amp( double strike, double dip, double slip, double azimuth, double takeoff )
{
	double temp;
	double SNDIP = sin(dip);
	double CSDIP = cos(dip);

	double SNEM = sin(takeoff);
	double CSEM = cos(takeoff);
	double SNEM2 = SNEM * SNEM;
	double SN2EM = 2.0 * SNEM * CSEM;

	azimuth -= strike;
	double SNAZ = sin(azimuth);

	temp  = cos(slip) * cos(azimuth) * (SNDIP * SNEM2 * 2.0 * SNAZ - CSDIP * SN2EM);
	temp += sin(slip) * (2.0 * SNDIP * CSDIP * ((CSEM * CSEM) - SNEM2 * (SNAZ * SNAZ)) + (CSDIP * CSDIP - SNDIP * SNDIP) * SN2EM * SNAZ);

	return temp;
}

/**
 * @brief
 *
 * @param results
 * @param nresult
 * @return int
 */
static int remove_approx_result( FPL_RESULT *results, int nresult )
{
/* */
	static const double appox_rad = FPLF_APPOX_DEG * FOCAL_GA_DEG2RAD;
/* */
	FPL_RESULT *curr_result = results;
	FPL_RESULT *next_result = curr_result + 1;
	FPL_RESULT *last_result = results + (nresult - 1);

/* */
	while ( curr_result < last_result ) {
	/* */
		if ( fabs(curr_result->strike - next_result->strike) < appox_rad && fabs(curr_result->dip - next_result->dip) < appox_rad && fabs(curr_result->rake - next_result->rake) < appox_rad ) {
		/* Try to keep the order of result, rather than using '*next_result = *last_result' here */
			memmove(next_result, next_result + 1, sizeof(FPL_RESULT) * (last_result - next_result));
			last_result--;
			nresult--;
		}
	/* */
		if ( ++next_result >= last_result ) {
			curr_result++;
			next_result = curr_result + 1;
		}
	}

	return nresult;
}

/**
 * @brief
 *
 * @param input
 * @param best
 * @param limit_diff
 * @return int
 */
static int is_same_group_result( FPL_RESULT *input, const FPL_RESULT *best, const double limit_diff )
{
	double strike_diff = fabs(input->strike - best->strike);
	double dip_diff = fabs(input->dip - best->dip);
	double rake_diff = fabs(input->rake - best->rake);

	if ( strike_diff < limit_diff && dip_diff < limit_diff && rake_diff < limit_diff ) {
	/* */
		return 1;
	}
	else if ( strike_diff < limit_diff && dip_diff < limit_diff && (FOCAL_GA_PI2 - rake_diff) < limit_diff ) {
	/* */
		if ( best->rake < limit_diff )
			input->rake -= FOCAL_GA_PI2;
		else
			input->rake += FOCAL_GA_PI2;
	/* */
		return 1;
	}
	else if ( (FOCAL_GA_PI2 - strike_diff) < limit_diff && dip_diff < limit_diff && rake_diff < limit_diff ) {
	/* */
		if ( best->strike < limit_diff )
			input->strike -= FOCAL_GA_PI2;
		else
			input->strike += FOCAL_GA_PI2;
	/* */
		return 1;
	}
	else if ( (FOCAL_GA_PI2 - strike_diff) < limit_diff && dip_diff < limit_diff && (FOCAL_GA_PI2 - rake_diff) < limit_diff ) {
	/* */
		if ( best->strike < limit_diff )
			input->strike -= FOCAL_GA_PI2;
		else
			input->strike += FOCAL_GA_PI2;
	/* */
		if ( best->rake < limit_diff )
			input->rake -= FOCAL_GA_PI2;
		else
			input->rake += FOCAL_GA_PI2;
	/* */
		return 1;
	}

	return 0;
}

/*
 *  compare_azimuth()
 */
static int compare_azimuth( const void *a, const void *b )
{
	double azimuth_a = ((FPL_OBSERVE *)a)->azimuth;
	double azimuth_b = ((FPL_OBSERVE *)b)->azimuth;

	if ( azimuth_a < azimuth_b )
		return -1;
	else if ( azimuth_a > azimuth_b )
		return 1;
	else
		return 0;
}

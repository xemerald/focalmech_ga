/**
 * @file fpl_find.h
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
#define FPLF_DECODE_RAD_CONST  0.01227184630308512983774470071594f
#define FPLF_APPOX_DEG         0.2f
/* */
#define FPLF_T_AXIS 0
#define FPLF_P_AXIS 1

typedef struct {
	double strike;
	double dip;
	double rake;
} FPL_RESULT;

typedef struct {
	double azimuth;
	double takeoff;
	double polarity;
} FPL_OBSERVE;

/**/
double fpl_find_ga( void *, int, int, int, int, double, double, FPL_RESULT *, int * );
int fpl_result_refine( FPL_RESULT *, FPL_RESULT *, int );
double fpl_quality_cal( FPL_OBSERVE *, const int, const double );
FPL_RESULT *fpl_dbcouple( FPL_RESULT *, FPL_RESULT [2], FPL_RESULT *, FPL_RESULT * );
FPL_OBSERVE fpl_observe_deg2rad( FPL_OBSERVE * );
FPL_OBSERVE fpl_observe_rad2deg( FPL_OBSERVE * );
FPL_RESULT fpl_result_deg2rad( FPL_RESULT * );
FPL_RESULT fpl_result_rad2deg( FPL_RESULT * );

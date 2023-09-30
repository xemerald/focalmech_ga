#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <float.h>
#include <math.h>

/* PGPLOT C library header */
#include <cpgplot.h>
#include <constants.h>
#include <fpl_func.h>

#define  COLOR_INDEX_BG    0
#define  COLOR_INDEX_FG    1

/**
 * @brief
 *
 * @param fullpath
 */
void pfocal_canva_open( const char *fullpath )
{
	char _path[MAX_PATH_STR];
/* */
	sprintf(_path, "%s/png", fullpath);
	if ( !cpgopen(_path) ) {
		fprintf(stderr, "pfocal_canva_open: ERROR!! Can't open %s\n", fullpath);
		return;
	}
/* */
	cpgscr(COLOR_INDEX_FG, 0.0, 0.0, 0.0);
	cpgscr(COLOR_INDEX_BG, 1.0, 1.0, 1.0);
	cpgpap(30, 1.0);
	cpgbbuf();
/* */
	cpgslw(12);
	cpgenv(-12.0, 12.0, -12.0, 12.0, 1, -1);
	cpgscf(3);
	cpgslw(16);
	cpgsfs(2);
	cpgcirc(0.0, 0.0, 10.0);
	cpgsfs(1);

	return;
}

/**
 * @brief
 *
 */
void pfocal_canva_close( void )
{
	cpgebuf();
	cpgend();
	return;
}

/**
 * @brief
 *
 * @param dbplane
 * @param radius
 */
void pfocal_dbplane_plot( FPL_RESULT dbplane[2], const double radius )
{
	const double r_sqrt_two = radius * sqrt(2.0);
	double sin_str;
	double cos_str;
	float _x;
	float _y;

/* */
	for ( int i = 0; i < 2; i++ ) {
		if ( dbplane[i].dip >= FOCAL_GA_HALF_PI )
			dbplane[i].dip -= 0.0001;
		else if ( dbplane[i].dip < 0.0001 )
			dbplane[i].dip += 0.0001;
	/* */
		sin_str = sin(dbplane[i].strike);
		cos_str = cos(dbplane[i].strike);
		_x = radius * sin_str;
		_y = radius * cos_str;
		cpgmove(_x, _y);
		for ( double theta = FOCAL_GA_HALF_PI; theta >= -FOCAL_GA_HALF_PI; theta -= FOCAL_GA_DEG2RAD ) {
			double cos_theta = cos(theta);
			double tmp1 = cos_theta * tan(dbplane[i].dip);
			double tmp2 = r_sqrt_two * sin(acos(tmp1 / sqrt(tmp1 * tmp1 + 1.0)) * 0.5);
			tmp1 = tmp2 * cos_theta;
			tmp2 *= sin(theta);
			_x = tmp1 * cos_str + tmp2 * sin_str;
			_y = -1.0 * tmp1 * sin_str + tmp2 * cos_str;
			cpgdraw(_x, _y);
		}
	}

	return;
}

/**
 * @brief
 *
 * @param ptaxis
 * @param radius
 */
void pfocal_pt_axis_plot( FPL_RESULT ptaxis[2], const double radius )
{
	const double r_sqrt_two = radius * sqrt(2.0);
	float tmp = radius / 40.0;

	cpgslw(16);
	cpgsch(1.0);
	for ( int i = 0; i < 2; i++ ) {
		double sii = r_sqrt_two * sin((FOCAL_GA_HALF_PI - ptaxis[i].dip) * 0.5);
		float _x = sii * sin(ptaxis[i].strike);
		float _y = sii * cos(ptaxis[i].strike);
		cpgtext(_x - tmp, _y - tmp, i == FPLF_P_AXIS ? "P" : "T");
	}

	return;
}

/**
 * @brief
 *
 * @param observes
 * @param nobserve
 * @param radius
 */
void pfocal_observe_plot( FPL_OBSERVE *observes, const int nobserve, const double radius )
{
	const double r_sqrt_two = radius * sqrt(2.0);
	float tmp = radius / 50.0;
	double rx;
	float _x;
	float _y;

/* */
	cpgsch(1.0);
	cpgslw(16);
	for ( int i = 0; i < nobserve; i++ ) {
		if ( fabs(observes[i].polarity - 0.0) < FLT_EPSILON )
			continue;
		if ( observes[i].takeoff > FOCAL_GA_HALF_PI ) {
			observes[i].azimuth += FOCAL_GA_PI;
			observes[i].takeoff = FOCAL_GA_PI - observes[i].takeoff;
			if ( observes[i].azimuth > FOCAL_GA_PI2 )
				observes[i].azimuth -= FOCAL_GA_PI2;
		}

		rx = r_sqrt_two * sin(observes[i].takeoff * 0.5);
		_x = rx * sin(observes[i].azimuth);
		_y = rx * cos(observes[i].azimuth);

		cpgsfs(observes[i].polarity < 0.0 ? 2 : 1);
		cpgcirc(_x, _y, tmp);
	}

	return;
}

/**
 * @brief
 *
 * @param year
 * @param mon
 * @param day
 * @param hour
 * @param min
 * @param sec
 * @param mag
 * @param lat
 * @param lon
 * @param depth
 */
void pfocal_eq_info_plot( int year, int mon, int day, int hour, int min, double sec, double mag, double lat, double lon, double depth )
{
	char str_buf[128] = { 0 };

	sprintf(str_buf,"%04d/%02d/%02d %02d:%02d:%05.2f M%3.1f", year, mon, day, hour, min, sec, mag);
	cpgsch(0.75);
	cpgtext(-11.5, 11.0, str_buf);

	sprintf(str_buf, "Lat.: %7.3f, Long.: %8.3f", lat, lon);
	cpgtext(-11.5, 10.2, str_buf);

	sprintf(str_buf, "Depth: %6.1f km", depth);
	cpgtext(-11.5, 9.4, str_buf);
	return;
}

/**
 * @brief
 *
 * @param dbplane
 * @param ptaxis
 * @param sdv
 * @param quality
 * @param fscore
 * @param radius
 */
void pfocal_plane_info_plot( FPL_RESULT dbplane[2], FPL_RESULT ptaxis[2], FPL_RESULT *sdv, double quality, double fscore, const double radius )
{
	float tmp = radius / 20.0;
	char str_buf[128] = { 0 };
	FPL_RESULT tmpres;

/* Draw the center cross */
	cpgsch(1.0);
	cpgmove(-tmp, 0.0);
	cpgdraw(tmp, 0.0);
	cpgmove(0.0, -tmp);
	cpgdraw(0., tmp);
/* Draw the N, S, E, W tick & text */
	tmp /= 2.0;
	cpgmove(0.0, radius - tmp);
	cpgdraw(0.0,radius + tmp);
	cpgtext(-tmp, radius + 1.5 * tmp, "N");
	cpgmove(0.0, -radius - tmp);
	cpgdraw(0.0, -radius + tmp);
	cpgtext(-tmp, -radius - 3.5 * tmp, "S");
	cpgmove(radius - tmp, 0.0);
	cpgdraw(radius + tmp, 0.0);
	cpgtext(radius + 1.5 * tmp, -tmp, "E");
	cpgmove(-radius - tmp, 0.0);
	cpgdraw(-radius + tmp, 0.0);
	cpgtext(-radius -3.5 * tmp, -tmp, "W");

	cpgsch(0.75);
	tmpres = fpl_result_rad2deg( &dbplane[0] );
	sprintf(str_buf, "A. Strike:%4d, Dip:%3d, Rake:%4d.", (int)tmpres.strike, (int)tmpres.dip, (int)tmpres.rake);
	cpgtext(-11.5, -10.5, str_buf);
	tmpres = fpl_result_rad2deg( &dbplane[1] );
	sprintf(str_buf, "B. Strike:%4d, Dip:%3d, Rake:%4d.", (int)tmpres.strike, (int)tmpres.dip, (int)tmpres.rake);
	cpgtext(-11.5, -11.2, str_buf);
	tmpres = fpl_result_rad2deg( &ptaxis[FPLF_P_AXIS] );
	sprintf(str_buf, "P-axis. Azimuth:%4d, Plunge:%3d.", (int)tmpres.strike, (int)tmpres.dip);
	cpgtext(1.3, -10.5, str_buf);
	tmpres = fpl_result_rad2deg( &ptaxis[FPLF_T_AXIS] );
	sprintf(str_buf, "T-axis. Azimuth:%4d, Plunge:%3d.", (int)tmpres.strike, (int)tmpres.dip);
	cpgtext(1.4, -11.2, str_buf);

	tmp = radius / 50.0;
	cpgsfs(2);
	cpgcirc(-11.0, 7.8, tmp);
	cpgtext(-10.7, 7.8 - tmp, "Down");
	cpgsfs(1);
	cpgcirc(-11.0, 8.6, tmp);
	cpgtext(-10.7, 8.6 - tmp, "Up");

	tmpres = fpl_result_rad2deg( sdv );
	cpgtext(3.5, 11.0, "Fault Plane A Uncertainty");
	sprintf(str_buf, "Strike+-%3d, Dip+-%3d, Rake+-%3d.", (int)tmpres.strike, (int)tmpres.dip, (int)tmpres.rake);
	cpgtext(0.8, 10.2, str_buf);
	sprintf(str_buf, "Quality Index = %5.2f", quality);
	cpgtext(4.8, 9.4, str_buf);
	tmp = (1.0 - fscore) * 0.5;
	sprintf(str_buf, "Misfit = %5.2f%%", tmp * 100.0);
	cpgtext(6.5, 8.6, str_buf);

	return;
}

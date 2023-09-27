#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>

#include <tko_azi_cal.h>
#include <ga_core.h>
#include <fpl_func.h>
#include <plot_focal.h>

int main( int argc, char **argv )
{
	FILE *fp = NULL;
	FPL_OBSERVE *obs = NULL;
	FPL_RESULT result[800];
	FPL_RESULT sdv[800];
	int nobs = 0;
	int nres = 0;
	char line[512] = { 0 };
	int azi, tko;
	char ud;
	double f_score = 0.0;
	if ( argc != 2 ) {
		printf("usage: %s <INPUT P-FILE>\n", argv[0]);
		return 0;
	}

	if ( !(fp = fopen(argv[1], "r")) ) {
		printf("Can't open the file %s\n", argv[1]);
		return -1;
	}

	while (fgets(line, 512, fp)) {
		printf("%s", line);
		nobs++;
	}
	nobs--;

	obs = calloc(nobs, sizeof(FPL_OBSERVE));
	rewind(fp);
	fgets(line, 512, fp);
	nobs = 0;
	while (fgets(line, 512, fp)) {
		sscanf(line, "%*s %*f %d %d%c %*s", &azi, &tko, &ud);
		if ( ud == ' ' )
			continue;
		obs[nobs].azimuth = azi;
		obs[nobs].takeoff = tko;
		obs[nobs].polarity = ud == '+' ? 1.0 : ud == '-' ? -1.0 : 0.0;
		obs[nobs] = fpl_observe_deg2rad( &obs[nobs] );
		nobs++;
	}
	fclose(fp);

	f_score = fpl_find_ga( obs, nobs, 20, 800, 3, 0.036, 0.72, result, &nres );
	double q = fpl_quality_cal( obs, nobs, f_score );
	printf("%f %d\n", f_score, nres);
	nres = fpl_result_refine(result, sdv, nres);
	q /= nres;
	printf("%f %d\n", f_score, nres);
	//for ( int i = 0; i < nres; i++ ) {
		FPL_RESULT dbresult[2];
		FPL_RESULT ptaxis[2];
		fpl_dbcouple(&result[0], dbresult, &ptaxis[FPLF_T_AXIS], &ptaxis[FPLF_P_AXIS]);
		pfocal_canva_open("/home/benyang/test.png");
		pfocal_dbplane_plot( dbresult, 10.0 );
		pfocal_pt_axis_plot( ptaxis, 10.0 );
		pfocal_observe_plot( obs, nobs, 10.0 );
		pfocal_eq_info_plot( 2020, 6, 26, 23, 27, 44.94, 5.05, 23.01, 120.91, 4.73 );
		pfocal_plane_info_plot( dbresult, ptaxis, sdv, q, f_score, 10.0 );
		pfocal_canva_close();
		//fpl_result_rad2deg( &dbresult[0] );
		//fpl_result_rad2deg( &dbresult[1] );
		//fpl_result_rad2deg( &taxis );
		//fpl_result_rad2deg( &paxis );
		//printf("p1: %f %f %f\n", dbresult[0].dip, dbresult[0].strike, dbresult[0].rake);
		//printf("p2: %f %f %f\n", dbresult[1].dip, dbresult[1].strike, dbresult[1].rake);
		//printf("t: %f %f %f\n", taxis.dip, taxis.strike, taxis.rake);
		//printf("p: %f %f %f\n", paxis.dip, paxis.strike, paxis.rake);
		//result[i].dip *= FPLF_RAD2DEG;
		//result[i].strike *= FPLF_RAD2DEG;
		//result[i].rake *= FPLF_RAD2DEG;
		//printf("%f %f %f\n", result[i].dip, result[i].strike, result[i].rake);
	//}

	free(obs);

	return 0;
}

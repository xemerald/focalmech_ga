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
	double azi, tko;
	double evdp;
	double evlo;
	double evla;
	double stdp;
	double stlo;
	double stla;

	tac_velmod_load(argv[1]);
	evla = atof(argv[2]);
	evlo = atof(argv[3]);
	evdp = atof(argv[4]);
	stla = atof(argv[5]);
	stlo = atof(argv[6]);
	stdp = atof(argv[7]);
	tac_main(&tko, &azi, evla, evlo, evdp, stla, stlo, stdp);
	tac_velmod_free();
	return 0;
}

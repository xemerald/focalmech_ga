/**
 * @file plot_focal.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-09-24
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <fpl_func.h>

void pfocal_canva_open( const char *filepath );
void pfocal_canva_close( void );
void pfocal_dbplane_plot( FPL_RESULT dbplane[2], const double radius );
void pfocal_pt_axis_plot( FPL_RESULT ptaxis[2], const double radius );
void pfocal_observe_plot( FPL_OBSERVE *observes, const int nobserve, const double radius );
void pfocal_eq_info_plot( int year, int mon, int day, int hour, int min, double sec, double mag, double lat, double lon, double depth );
void pfocal_plane_info_plot( FPL_RESULT dbplane[2], FPL_RESULT ptaxis[2], double strike_sdv, double dip_sdv, double rake_sdv, double q, double fscore, const double radius );

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

void pfocal_canva_open( const char * );
void pfocal_canva_close( void );
void pfocal_dbplane_plot( FPL_RESULT [2], const double );
void pfocal_pt_axis_plot( FPL_RESULT [2], const double );
void pfocal_observe_plot( FPL_OBSERVE *, const int, const double );
void pfocal_eq_info_plot( int, int, int, int, int, double, double, double, double, double );
void pfocal_plane_info_plot( FPL_RESULT [2], FPL_RESULT [2], FPL_RESULT *, double, double, const double );

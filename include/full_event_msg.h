/**
 * @file full_event_msg.h
 * @author your name (you@domain.com)
 * @brief Detailed definition of TYPE_FULL_EVENT,
 * @version 0.1
 * @date 2023-09-26
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <trace_buf.h>

#define PICK_PHASE_NAME_LEN  9
#define PICK_PAMP_TYPE_NUM   3

/*
 * pick_id     DBMS ID of pick (unsigned long)
 * seq         Sequence number of pick itself
 * station     Station code (upper case)
 * channel     Component code (upper case)
 * network     Network code (upper case)
 * location    Location code
 * phase_name  Phase code (e.g., 'P', 'S', or 'PKP')
 * polarity    First-motion descriptor (U,D,' ','?')
 * weight      Pick weight or quality (0-4)
 * picktime    Time of pick in seconds since 1970
 * residual    Travel-time residual in seconds
 * pamp        P amplitudes in three different types
 * dist        Source-receiver distance in decimal degrees
 * azi         Receiver azimuth in degrees (from the source)
 * tko         Receiver takeoff angle in degrees (from the source)
 */
typedef struct {
	int     pick_id;
	int     seq;
	char    station[TRACE2_STA_LEN];
	char    channel[TRACE2_CHAN_LEN];
	char    network[TRACE2_NET_LEN];
	char    location[TRACE2_LOC_LEN];
	char    phase_name[PICK_PHASE_NAME_LEN];
	char    polarity;
	char    weight;
	char    padding[6];
	double  latitude;
	double  longitude;
	double  elevation;
	double  picktime;
	double  residual;
	double  pamp[PICK_PAMP_TYPE_NUM];
	double  dist;
	double  azi;
	double  tko;
} FULL_EVENT_PICK_MSG;


#define FULL_EVENT_ID_LEN  32

/*
 * event_id    Event identifier
 * seq         Sequence number of event msg
 * origin      Event origin time in epoch seconds
 * evlat       Geographical event latitude in signed decimal degrees (WGS84)
 * evlon       Geographical event longitude in signed decimal degrees (WGS84)
 * evdepth     Event depth in kilometers (WGS84)
 * gap         Largest azimuthal gap in degrees
 * nsta        Number of stations associated
 * npha        Number of phases associated
 * suse        Number of stations used in the solution
 * puse        Number of phases used in the solution
 * dmin        Distance to the nearest station in degrees
 * depth_flag  Fixed depth flag (T if depth was held, F otherwise)*
 * oterr       90% marginal confidence interval for origin time
 * laterr      90% marginal confidence interval for latitude
 * lonerr      90% marginal confidence interval for longitude
 * deperr      90% marginal confidence interval for depth
 * se          Standard error of the residuals in seconds
 * errh        Maximum horizontal projection of the error ellipsoid in kilometers
 * errz        Maximum vertical projection of the error ellipsoid in kilometers
 * avh         Equivalent radius of the horizontal error ellipse in kilometers
 * q           Quality flag (i.e., 'A', 'B', 'C', 'D')
 * axis1-3     Length in kilometers of the principle axies of the error ellipsoid
 * az1-3       Azimuth in degrees of the principle axies of the error ellipsoid
 * dp1-3       Dip in degrees of the principle axies of the error ellipsoid
 * npicks      Number of picks in this full message
 */
typedef struct {
   char    event_id[FULL_EVENT_ID_LEN];
   int     seq;
   long    origin_id; /* which edition of event_id contained herein */
   double  origin_time;
   double  evlat;
   double  evlon;
   double  evdepth;
   double  magnitude;
   double  gap;
   int     nsta;
   int     npha;
   int     suse;
   int     puse;
   double  dmin;
   char    depth_flag;
   double  oterr;
   double  laterr;
   double  lonerr;
   double  deperr;
   double  se;
   double  errh;
   double  errz;
   double  avh;
   char    q;
   double  axis[3];
   double  az[3];
   double  dp[3];
   int     npicks;
} FULL_EVENT_MSG_HEADER;

/* max pcks in the message */
#define FULL_EVENT_MAX_PICKS   1024

typedef struct {
	FULL_EVENT_MSG_HEADER header;
	FULL_EVENT_PICK_MSG   picks[FULL_EVENT_MAX_PICKS];
} FULL_EVENT_MSG;

#define FULL_EVENT_SIZE        (sizeof(FULL_EVENT_MSG))
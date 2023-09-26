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
/* */
#include <stdint.h>
/* */
#define  MAX_STR_SIZE            256
#define  MAX_PATH_STR            512
#define  MAX_IN_MSG              8
/* */
#define  DUMMY_MAG              -9.0
/* */
#define  PLOT_BEACHBALL_RADIUS  10.0f
/* */
typedef struct {
	char   script[MAX_PATH_STR];
	double min_magnitude;
} POSCRIPT;
/* */
typedef struct {
	char    title[MAX_STR_SIZE];
	char    caption[MAX_STR_SIZE];
	double  min_magnitude;
/* */
	uint8_t smaptype;
	uint8_t gmflag[MAX_IN_MSG];
	void   *gmptr[MAX_IN_MSG];
} PLOTSMAP;

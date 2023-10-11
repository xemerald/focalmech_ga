/*
 *
 */

/* */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <search.h>
#include <time.h>
#include <ctype.h>
/* */
#include <earthworm.h>
#include <transport.h>
#include <kom.h>
/* */
#include <early_event_msg.h>

/* */
#define DEFAULT_RING_NAME  "HYPO_RING"

/**
 * @brief
 *
 */
typedef struct {
	char   station[8];
	double latitude;
	double longitude;
	double elevation;
} STA_INFO;


/* */
static int read_pfile( const char *, EARLY_EVENT_MSG * );
static int read_stalist_file( const char * );
static void send_fullevent_message( void *, size_t );
static char *_rtrim( char * );
static int compare_sta( const void *, const void * );
static double get_unixtime_ot( int, int, int, int, int, double );

/* */
static SHM_INFO      Region;      /* The shared memory region   */
static unsigned char InstId;      /* local installation id      */
static unsigned char MyModId;
static unsigned char TypeEarlyEvent;

static void *Root = NULL;

/*
 *
 */
int main( int argc, char *argv[] )
{
	char ringname[MAX_RING_STR];
	int  arg_index = 1;
	int  ringkey;
	EARLY_EVENT_MSG evt_msg;

/* Check the number of arguments  */
	if ( argc < 2 ) {
		fprintf(stderr, "Usage: pfile2fevent [-r RING_NAME] <STA_LIST> <P_FILE>\n");
		return -1;
	}

	if ( strcmp(argv[1], "-r") == 0 ) {
		strcpy(ringname, argv[2]);
		arg_index = 3;
	}
	else {
		strcpy(ringname, DEFAULT_RING_NAME);
	}

/* Look up ids in earthworm.d tables */
	if ( GetLocalInst(&InstId) != 0 ) {
		fprintf(stderr, "pfile2fevent: Error getting local installation id; exiting!\n");
		return -1;
	}
	if ( GetType("TYPE_EARLY_EVENT", &TypeEarlyEvent) != 0 ) {
		fprintf(stderr, "pfile2fevent: Invalid message type <TYPE_EARLY_EVENT> exiting!\n");
		return -1;
	}

/* Validate transport ring name */
	if ( (ringkey = GetKey(ringname)) == -1 ) {
		fprintf(stderr, "pfile2fevent: Invalid ring name <%s>. Exiting.\n", ringname);
		return -1;
	}

/* Attach to shared memory ring */
	tport_attach( &Region, ringkey );
/* */

	read_stalist_file( argv[arg_index++] );
	if ( !read_pfile( argv[arg_index], &evt_msg ) ) {
		strcpy(evt_msg.header.event_id, strrchr(argv[arg_index], '/') + 1);
		evt_msg.header.seq = getpid();
		send_fullevent_message( &evt_msg, sizeof(EARLY_EVENT_MSG_HEADER) + evt_msg.header.npicks * sizeof(EARLY_PICK_MSG) );
	}
	tport_detach( &Region );

	return 0;
}

static int read_pfile( const char *filepath, EARLY_EVENT_MSG *evt_msg )
{
	FILE *fp = NULL;
	STA_INFO key, *staptr;
	int  npicks = 0;
	int  azi, tko;
	char ud;
	char line[512] = { 0 };
	char frac_1[8] = { 0 };
	char frac_2[8] = { 0 };
	char frac_3[8] = { 0 };
	char frac_4[8] = { 0 };
	char frac_5[8] = { 0 };
	char frac_6[8] = { 0 };
	char frac_7[8] = { 0 };
	char frac_8[8] = { 0 };
	char frac_9[8] = { 0 };
	char frac_10[8] = { 0 };
	char frac_11[8] = { 0 };
	char frac_12[8] = { 0 };

	if ( evt_msg == NULL )
		return -1;

	if ( !(fp = fopen(filepath, "r")) ) {
		fprintf(stderr, "read_pfile: Can't open the file %s\n", filepath);
		return -1;
	}
/* */
	fgets(line, 512, fp);
	sscanf(
		line, "%5c%2c%2c%2c%2c%6c%2c%5c%3c%5c%6c%4c %*s",
		frac_1, frac_2, frac_3, frac_4, frac_5, frac_6,
		frac_7, frac_8, frac_9, frac_10, frac_11, frac_12
	);
/* */
	evt_msg->header.origin_time = get_unixtime_ot( atoi(frac_1), atoi(frac_2), atoi(frac_3), atoi(frac_4), atoi(frac_5), atof(frac_6) );
	evt_msg->header.evlat   = atof(frac_7) + atof(frac_8) / 60.0;
	evt_msg->header.evlon   = atof(frac_9) + atof(frac_10) / 60.0;
	evt_msg->header.evdepth = atof(frac_11);
	evt_msg->header.mag     = atof(frac_12);
/* */
	npicks = 0;
	while ( fgets(line, 512, fp) ) {
		sscanf(line, " %s %*f %d %d%c %*s", frac_1, &azi, &tko, &ud);
		_rtrim( frac_1 );
		strcpy(key.station, frac_1);
		if ( (staptr = tfind(&key, &Root, compare_sta)) ) {
			staptr = *(STA_INFO **)staptr;
			strcpy(evt_msg->picks[npicks].station, frac_1);
			evt_msg->picks[npicks].latitude  = staptr->latitude;
			evt_msg->picks[npicks].longitude = staptr->longitude;
			evt_msg->picks[npicks].elevation = staptr->elevation;
			evt_msg->picks[npicks].picktime  = evt_msg->header.origin_time;
			evt_msg->picks[npicks].polarity  = ud == '+' ? 'U' : ud == '-' ? 'D' : ' ';
			evt_msg->picks[npicks].azi       = azi;
			evt_msg->picks[npicks].tko       = tko;
			evt_msg->picks[npicks].phase_name[0] = 'P';
			evt_msg->picks[npicks].phase_name[1] = '\0';

			npicks++;
		}
	}
	fclose(fp);
	evt_msg->header.npicks = evt_msg->header.nsta = evt_msg->header.npha = npicks;

	return 0;
}

/*
 *
 */
static int read_stalist_file( const char *filepath )
{
	FILE *fp = NULL;
	char line[512] = { 0 };
	char frac_1[8];
	char frac_2[8];
	char frac_3[8];
	char frac_4[8];
	char frac_5[8];
	char frac_6[8];
	STA_INFO  key;
	STA_INFO *new = NULL;
/* */
	if ( !(fp = fopen(filepath, "r")) ) {
		fprintf(stderr, "read_stalist_file: Can't open the file %s\n", filepath);
		return -1;
	}
/* Skip the first line */
	fgets(line, 512, fp);
/* */
	while ( fgets(line, 512, fp) ) {
		sscanf(
			line, "%4c%2c%5c%3c%5c %s %*s",
			frac_1, frac_2, frac_3, frac_4, frac_5, frac_6
		);
	/* */
		_rtrim( frac_1 );
		strcpy(key.station, frac_1);
		key.latitude = atof(frac_2) + atof(frac_3) / 60.0;
		key.longitude = atof(frac_4) + atof(frac_5) / 60.0;
		key.elevation = atof(frac_6);
	/* */
		if ( tfind(&key, &Root, compare_sta) == NULL ) {
			new = calloc(1, sizeof(STA_INFO));
			*new = key;
			if ( tsearch(new, &Root, compare_sta) == NULL ) {
				fprintf(stderr, "read_stalist_file: Can't insert the station %s\n", key.station);
				return -1;
			}
		}
	}
	fclose(fp);

	return 0;
}

/**
 * @brief
 *
 * @param omsg
 * @param msg_size
 */
static void send_fullevent_message( void *omsg, size_t msg_size )
{
	MSG_LOGO logo = { .type = TypeEarlyEvent, .mod = MyModId, .instid = InstId };

/* Send cwaxml message to transport ring */
	if ( tport_putmsg(&Region, &logo, msg_size, (char *)omsg) != PUT_OK )
		fprintf(stderr, "pfile2fevent: Error sending full event message to transport region.\n");
	else
		fprintf(stdout, "pfile2fevent: Sending full event message success!\n");

	return;
}

/*
 *
 */
static char *_rtrim( char *str )
{
	int len;
	char *c;

/* */
	if ( !str || *str == '\0' )
		return str;

	len = strlen(str);
	for ( c = str + len - 1; c >= str && isspace(*c); c-- )
		*c = '\0';

	return str;
}

/**
 * @brief
 *
 * @param a
 * @param b
 * @return int
 */
static int compare_sta( const void *a, const void *b )
{
	char *sta_a = ((STA_INFO *)a)->station;
	char *sta_b = ((STA_INFO *)b)->station;

	return strcmp(sta_a, sta_b);
}

/**
 * @brief Get the unixtime ot object
 *
 * @param year
 * @param mon
 * @param day
 * @param hour
 * @param min
 * @param sec
 * @return double
 */
static double get_unixtime_ot( int year, int mon, int day, int hour, int min, double sec )
{
	struct tm st;
	double result = sec;

	st.tm_year = year - 1900;
	st.tm_mon  = mon - 1;
	st.tm_mday = day;
	st.tm_hour = hour;
	st.tm_min  = min;
	st.tm_sec  = (int)sec;
/* */
	result -= st.tm_sec;
	result += timegm(&st);

	return result;
}

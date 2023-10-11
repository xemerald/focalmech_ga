#ifdef _OS2
#define INCL_DOSMEMMGR
#define INCL_DOSSEMAPHORES
#include <os2.h>
#endif

/* Standard C header include */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <stdatomic.h>
#include <time.h>
#include <math.h>
#include <sys/stat.h>
/* Earthworm environment header include */
#include <earthworm.h>
#include <kom.h>
#include <transport.h>
#include <lockfile.h>
/* Local header include */
#include <focalmech_ga.h>
#include <early_event_msg.h>
#include <thrd_pool.h>
#include <fpl_func.h>
#include <plot_focal.h>
#include <tko_azi_cal.h>

/* Functions in this source file
 *******************************/
static void focalmech_ga_config( char * );
static void focalmech_ga_lookup( void );
static void focalmech_ga_status( unsigned char, short, char * );
static void focalmech_ga_end( void );  /* Free all the local memory & close socket */

static void  proc_event( EARLY_EVENT_MSG *evt_msg );
static void  proc_picks_tac( void * );
static void  cal_pick_tko_azi( EARLY_PICK_MSG *, const double, const double, const double );
static int   pack_picks_to_observes( FPL_OBSERVE **, EARLY_EVENT_MSG * );
static int   cal_focal_ga( FPL_RESULT *, FPL_RESULT *, double *, double *, FPL_OBSERVE *, const int );
static void  plot_focal_result( FPL_RESULT *, FPL_RESULT *, FPL_RESULT *, FPL_OBSERVE *, const int, const double, const double, void *, const char * );
static char *gen_script_command_args( char *, void *, FPL_RESULT [2], FPL_RESULT [2], double, double, int );
static char *mk_outdir_by_evt( char *, const char *, void * );
static char *gen_focalplot_fullpath( char *, const char *, void * );
static void  split_time( int *, int *, int *, int *, int *, double *, const double );

static SHM_INFO Region;      /* shared memory region to use for i/o    */

#define MAXLOGO   4

static MSG_LOGO Getlogo[MAXLOGO];   /* array for requesting module,type,instid */
static pid_t    MyPid;              /* for restarts by startstop               */

/* Thread things */
static tpool_t *ThrdPool = NULL;

#define MAX_POST_SCRIPTS        5

/* Things to read or derive from configuration file */
static char     RingName[MAX_RING_STR];     /* name of transport ring for i/o    */
static char     MyModName[MAX_MOD_STR];     /* speak as this module name/id      */
static uint8_t  LogSwitch;                  /* 0 if no logfile should be written */
static uint64_t HeartBeatInterval;          /* seconds between heartbeats        */
static uint8_t  ThreadsNum   = 1;
static uint8_t  RemoveSwitch = 0;
static short    nLogo = 0;
static char     ReportPath[MAX_PATH_STR];
static POSCRIPT PostScripts[MAX_POST_SCRIPTS];
static uint8_t  NumPostScripts  = 0;                /* Number of exec. scripts */
static uint32_t MinPickPolarity = 3;
static uint8_t  IterationNum    = 20;
static uint32_t Population      = 800;
static uint8_t  MutateBits      = 3;
static double   ReprodutionRate = 0.036;
static double   MutateRate      = 0.72;
/* Things to look up in the earthworm.h tables with getutil.c functions */
static long          RingKey;       /* key of transport ring for i/o     */
static unsigned char InstId;        /* local installation id             */
static unsigned char MyModId;       /* Module Id for this program        */
static unsigned char TypeHeartBeat;
static unsigned char TypeError;
static unsigned char TypeEarlyEvent;

/* Error messages used by focalmech_ga */
#define ERR_MISSMSG       0   /* message missed in transport ring       */
#define ERR_TOOBIG        1   /* retreived msg too large for buffer     */
#define ERR_NOTRACK       2   /* msg retreived; tracking limit exceeded */
#define ERR_QUEUE         3   /* error queueing message for sending      */
static char  Text[150];       /* string for log/error messages          */

/**
 * @brief
 *
 */
typedef struct {
	void      *buffer;
	int        proc_thrds;
	int        assign_idx;
	atomic_int is_finish;
} _PICKS_TAC_ARG;

/**
 * @brief
 *
 * @param argc
 * @param argv
 * @return int
 */
int main( int argc, char **argv )
{
	time_t timeNow;          /* current time                  */
	time_t time_lastbeat;     /* time last heartbeat was sent  */
	char  *lockfile;
	int    lockfile_fd;

	int      res;
	long     recsize;      /* size of retrieved message     */
	MSG_LOGO reclogo;      /* logo of retrieved message     */

	uint8_t  buffer[EARLY_EVENT_SIZE];

/* Check command line arguments */
	if ( argc != 2 ) {
		fprintf( stderr, "Usage: focalmech_ga <configfile>\n" );
		exit( 0 );
	}
/* Initialize name of log-file & open it */
	logit_init( argv[1], 0, 256, 1 );
/* Read the configuration file(s) */
	focalmech_ga_config( argv[1] );
	logit("", "%s: Read command file <%s>\n", argv[0], argv[1]);
/* Look up important info from earthworm.h tables */
	focalmech_ga_lookup();
/* Reinitialize logit to desired logging level */
	logit_init( argv[1], 0, 256, LogSwitch );

	lockfile = ew_lockfile_path(argv[1]);
	if ( (lockfile_fd = ew_lockfile(lockfile) ) == -1) {
		fprintf(stderr, "one instance of %s is already running, exiting\n", argv[0]);
		exit(-1);
	}
/* Get process ID for heartbeat messages */
	MyPid = getpid();

	if ( MyPid == -1 ) {
		logit("e","focalmech_ga: Cannot get pid. Exiting.\n");
		exit (-1);
	}

/* Create a Mutex to control access to queue & initialize the message queue */
	tpool_init( &ThrdPool, ThreadsNum, ThreadsNum * 2 );
/* Attach to Input/Output shared memory ring */
	tport_attach( &Region, RingKey );
	logit("", "focalmech_ga: Attached to public memory region %s: %ld\n", RingName, RingKey );
/* Flush the transport ring */
	tport_flush( &Region, Getlogo, nLogo, &reclogo );

/* Force a heartbeat to be issued in first pass thru main loop */
	time_lastbeat = time(&timeNow) - HeartBeatInterval - 1;
/*----------------------- setup done; start main loop -------------------------*/
	while ( 1 ) {
	/* Send focalmech_ga's heartbeat */
		if ( time(&timeNow) - time_lastbeat >= (int64_t)HeartBeatInterval ) {
			time_lastbeat = timeNow;
			focalmech_ga_status( TypeHeartBeat, 0, "" );
		}

		do {
		/* See if a termination has been requested */
			res = tport_getflag( &Region );
			if ( res == TERMINATE || res == MyPid ) {
			/* Write a termination msg to log file */
				logit( "t", "focalmech_ga: Termination requested; exiting!\n" );
				fflush( stdout );
			/* */
				goto exit_procedure;
			}

		/* Get msg & check the return code from transport */
			res = tport_getmsg(&Region, Getlogo, nLogo, &reclogo, &recsize, (char *)buffer, EARLY_EVENT_SIZE);
		/* No more new messages     */
			if ( res == GET_NONE ) {
				break;
			}
		/* Next message was too big */
			else if ( res == GET_TOOBIG ) {
			/* Complain and try again   */
				sprintf(
					Text, "Retrieved msg[%ld] (i%u m%u t%u) too big for Buffer[%ld]",
					recsize, reclogo.instid, reclogo.mod, reclogo.type, (long)EARLY_EVENT_SIZE
				);
				focalmech_ga_status( TypeError, ERR_TOOBIG, Text );
				continue;
			}
		/* Got a msg, but missed some */
			else if ( res == GET_MISS ) {
				sprintf(
					Text, "Missed msg(s)  i%u m%u t%u  %s.", reclogo.instid, reclogo.mod, reclogo.type, RingName
				);
				focalmech_ga_status( TypeError, ERR_MISSMSG, Text );
			}
		/* Got a msg, but can't tell */
			else if ( res == GET_NOTRACK ) {
			/* If any were missed        */
				sprintf(
					Text, "Msg received (i%u m%u t%u); transport.h NTRACK_GET exceeded",
					reclogo.instid, reclogo.mod, reclogo.type
				);
				focalmech_ga_status( TypeError, ERR_NOTRACK, Text );
			}

		/* Process the message */
			if ( reclogo.type == TypeEarlyEvent ) {
			/* Debug */
				/* logit( "", "%s", rec ); */
			/* */
				proc_event( (EARLY_EVENT_MSG *)buffer );
			}
		} while ( 1 );
	/* no more messages; wait for new ones to arrive */
		sleep_ew(50);
	} /* while( 1 ) */
/*-----------------------------end of main loop-------------------------------*/
exit_procedure:
/* detach from shared memory */
	sleep_ew(500);
	focalmech_ga_end();

	ew_unlockfile(lockfile_fd);
	ew_unlink_lockfile(lockfile);

	return 0;
}

/*
 * focalmech_ga_config() - processes command file(s) using kom.c functions;
 *                      exits if any errors are encountered.
 */
static void focalmech_ga_config( char *configfile )
{
	int   ncommand;     /* # of required commands you expect to process   */
	char  init[10];     /* init flags, one byte for each required command */
	int   nmiss;        /* number of required commands that were missed   */
	char *com;
	char *str;
	int   nfiles;
	int   success;
	int   i;
	char  filepath[MAX_PATH_STR];

/* Set to zero one init flag for each required command */
	ncommand = 7;
	for ( i = 0; i < ncommand; i++ )
		init[i] = 0;

/* Open the main configuration file */
	nfiles = k_open( configfile );
	if ( nfiles == 0 ) {
		logit("e", "focalmech_ga: Error opening command file <%s>; exiting!\n", configfile);
		exit(-1);
	}

/* Process all command files */
/* While there are command files open */
	while ( nfiles > 0 ) {
	/* Read next line from active file  */
		while ( k_rd() ) {
			com = k_str();  /* Get the first token from line */
		/* Ignore blank lines & comments */
			if( !com )
				continue;
			if( com[0] == '#' )
				continue;

		/* Open a nested configuration file */
			if ( com[0] == '@' ) {
				success = nfiles+1;
				nfiles  = k_open(&com[1]);
				if ( nfiles != success ) {
					logit("e", "focalmech_ga: Error opening command file <%s>; exiting!\n", &com[1]);
					exit(-1);
				}
				continue;
			}

		/* Process anything else as a command */
		/* 0 */
			if ( k_its("LogFile") ) {
				LogSwitch = k_int();
				init[0] = 1;
			}
		/* 1 */
			else if ( k_its("MyModuleId") ) {
				str = k_str();
				if ( str )
					strcpy( MyModName, str );
				init[1] = 1;
			}
		/* 2 */
			else if ( k_its("RingName") ) {
				str = k_str();
				if ( str )
					strcpy( RingName, str );
				init[2] = 1;
			}
		/* 3 */
			else if ( k_its("HeartBeatInterval") ) {
				HeartBeatInterval = k_long();
				init[3] = 1;
			}
			else if ( k_its("ThreadsNum") ) {
				if ( (ThreadsNum = k_int()) > MAX_THREADS_NUM )
					ThreadsNum = MAX_THREADS_NUM;
				else if ( ThreadsNum < 1 )
					ThreadsNum = 1;
				logit(
					"o", "focalmech_ga: This process will use %d thread(s) for parallel processing.\n",
					ThreadsNum
				);
			}
			else if ( k_its("RemoveFocalPlot") ) {
				RemoveSwitch = k_int();
				logit(
					"o", "focalmech_ga: This process will %s the files after posted.\n",
					RemoveSwitch ? "remove" : "keep"
				);
			}
			else if ( k_its("MinPickPolarity") ) {
				MinPickPolarity = k_int();
				logit(
					"o", "focalmech_ga: The minimum number of pick polarity has been setted to %d.\n",
					MinPickPolarity
				);
			}
			else if ( k_its("IterationNum") ) {
				IterationNum = k_int();
				logit(
					"o", "focalmech_ga: The iteration times of GA has been setted to %d.\n",
					IterationNum
				);
			}
			else if ( k_its("Population") ) {
				Population = k_int();
				logit(
					"o", "focalmech_ga: The population size of GA has been setted to %d.\n",
					Population
				);
			}
			else if ( k_its("MutateBits") ) {
				MutateBits = k_int();
				logit(
					"o", "focalmech_ga: The mutate bits of GA has been setted to %d\n",
					MutateBits
				);
			}
			else if ( k_its("ReprodutionRate") ) {
				ReprodutionRate = k_val();
				logit(
					"o", "focalmech_ga: The reproduction rate of GA has been setted to %5.2lf%%\n",
					ReprodutionRate * 100.0
				);
			}
			else if ( k_its("MutateRate") ) {
				MutateRate = k_val();
				logit(
					"o", "focalmech_ga: The mutation rate of GA has been setted to %5.2lf%%\n",
					MutateRate * 100.0
				);
			}
		/* 5 */
			else if ( k_its("ReportPath") ) {
				str = k_str();
				if ( str )
					strcpy(ReportPath, str);
				if ( ReportPath[strlen(ReportPath) - 1] != '/' )
					strncat(ReportPath, "/", 1);

				logit("o", "focalmech_ga: Report path %s\n", ReportPath);
				init[4] = 1;
			}
			else if ( k_its("VelocityModelFile") ) {
				str = k_str();
				if ( str )
					strcpy(filepath, str);
				logit("o", "focalmech_ga: 3D velocity model file: %s\n", filepath);

				if ( tac_velmod_load( filepath ) ) {
					logit("e", "focalmech_ga: Error reading 3D velocity model file; exiting!\n");
					exit(-1);
				}
				else {
					logit("o", "focalmech_ga: Reading 3D velocity model file finish!\n");
				}
				init[5] = 1;
			}
			else if ( k_its("PostScriptWithMinMag") ) {
				if ( NumPostScripts < MAX_POST_SCRIPTS ) {
					str = k_str();
					if ( str )
						strcpy(PostScripts[NumPostScripts].script, str);
					PostScripts[NumPostScripts].min_magnitude = k_val();
					logit(
						"o", "focalmech_ga: Using script(%d): %s!\n",
						NumPostScripts, PostScripts[NumPostScripts].script
					);
					NumPostScripts++;
				}
				else {
					logit("e", "focalmech_ga: Excessive number of post scripts. Exiting!\n");
					exit(-1);
				}
			}
			else if ( k_its("PostScript") ) {
				if ( NumPostScripts < MAX_POST_SCRIPTS ) {
					str = k_str();
					if ( str )
						strcpy(PostScripts[NumPostScripts].script, str);
					PostScripts[NumPostScripts].min_magnitude = DUMMY_MAG;
					logit(
						"o", "focalmech_ga: Using script(%d): %s!\n", NumPostScripts, PostScripts[NumPostScripts].script
					);
					NumPostScripts++;
				}
				else {
					logit("e", "focalmech_ga: Excessive number of post scripts. Exiting!\n");
					exit(-1);
				}
			}
		/* Enter installation & module to get event messages from */
		/* 8 */
			else if( k_its("GetEventsFrom") ) {
				if ( nLogo >= MAXLOGO ) {
					logit("e", "focalmech_ga: Too many <GetEventsFrom> commands in <%s>", configfile);
					logit("e", "; max=%d; exiting!\n", (int) MAXLOGO);
					exit(-1);
				}
				if ( (str = k_str()) ) {
					if ( GetInst(str, &Getlogo[nLogo].instid) != 0 ) {
						logit("e", "focalmech_ga: Invalid installation name <%s>", str);
						logit("e", " in <GetEventsFrom> cmd; exiting!\n");
						exit(-1);
					}
				}
				if ( (str = k_str()) ) {
					if ( GetModId(str, &Getlogo[nLogo].mod) != 0 ) {
						logit("e", "focalmech_ga: Invalid module name <%s>", str);
						logit("e", " in <GetEventsFrom> cmd; exiting!\n");
						exit(-1);
					}
				}
				if ( (str = k_str()) ) {
					if ( GetType(str, &Getlogo[nLogo].type) != 0 ) {
						logit("e", "focalmech_ga: Invalid message type name <%s>", str);
						logit("e", " in <GetEventsFrom> cmd; exiting!\n");
						exit(-1);
					}
				}
				nLogo++;
				init[6] = 1;
			}
		/* Unknown command */
			else {
				logit("e", "focalmech_ga: <%s> Unknown command in <%s>.\n", com, configfile);
				continue;
			}

		/* See if there were any errors processing the command */
			if ( k_err() ) {
				logit("e", "focalmech_ga: Bad <%s> command in <%s>; exiting!\n", com, configfile);
				exit(-1);
			}
		}
		nfiles = k_close();
	}

/* After all files are closed, check init flags for missed commands */
	nmiss = 0;
	for ( i = 0; i < ncommand; i++ )
		if ( !init[i] )
			nmiss++;
/* */
	if ( nmiss ) {
		logit( "e", "focalmech_ga: ERROR, no " );
		if ( !init[0] ) logit("e", "<LogFile> "           );
		if ( !init[1] ) logit("e", "<MyModuleId> "        );
		if ( !init[2] ) logit("e", "<RingName> "          );
		if ( !init[3] ) logit("e", "<HeartBeatInterval> " );
		if ( !init[4] ) logit("e", "<ReportPath> "        );
		if ( !init[5] ) logit("e", "<VelocityModelFile> " );
		if ( !init[6] ) logit("e", "any <GetEventsFrom> " );
		logit("e", "command(s) in <%s>; exiting!\n", configfile);
		exit(-1);
	}

	return;
}

/*
 * focalmech_ga_lookup() - Look up important info from earthworm.h tables
 */
static void focalmech_ga_lookup( void )
{
/* Look up keys to shared memory regions */
   if( ( RingKey = GetKey(RingName) ) == -1 ) {
	   fprintf(stderr, "focalmech_ga:  Invalid ring name <%s>; exiting!\n", RingName);
	   exit(-1);
   }
   /* Look up installations of interest */
   	if ( GetLocalInst(&InstId) != 0 ) {
   		fprintf(stderr, "focalmech_ga: error getting local installation id; exiting!\n");
   		exit(-1);
   	}
/* Look up modules of interest */
	if ( GetModId(MyModName, &MyModId) != 0 ) {
		fprintf(stderr, "focalmech_ga: Invalid module name <%s>; exiting!\n", MyModName);
		exit(-1);
	}
/* Look up message types of interest */
	if ( GetType("TYPE_HEARTBEAT", &TypeHeartBeat) != 0 ) {
		fprintf(stderr, "focalmech_ga: Invalid message type <TYPE_HEARTBEAT>; exiting!\n");
		exit(-1);
	}
	if ( GetType("TYPE_ERROR", &TypeError) != 0 ) {
		fprintf(stderr, "focalmech_ga: Invalid message type <TYPE_ERROR>; exiting!\n");
		exit(-1);
	}
	if ( GetType( "TYPE_EARLY_EVENT", &TypeEarlyEvent ) != 0 ) {
		fprintf(stderr, "focalmech_ga: Invalid message type <TYPE_EARLY_EVENT>; exiting!\n");
		exit(-1);
	}

	return;
}

/*
 * focalmech_ga_status() - builds a heartbeat or error message & puts it into
 *                      shared memory.  Writes errors to log file & screen.
 */
static void focalmech_ga_status( unsigned char type, short ierr, char *note )
{
	MSG_LOGO    logo;
	char        msg[512];
	uint64_t    size;
	time_t      t;

/* Build the message */
	logo.instid = InstId;
	logo.mod    = MyModId;
	logo.type   = type;

	time(&t);

	if ( type == TypeHeartBeat ) {
		sprintf(msg, "%ld %ld\n", (long)t, (long)MyPid);
	}
	else if( type == TypeError ) {
		sprintf(msg, "%ld %hd %s\n", (long)t, ierr, note);
		logit("et", "focalmech_ga: %s\n", note);
	}

	size = strlen(msg);  /* don't include the null byte in the message */

/* Write the message to shared memory */
	if ( tport_putmsg(&Region, &logo, size, msg) != PUT_OK ) {
		if ( type == TypeHeartBeat ) {
			logit("et", "focalmech_ga:  Error sending heartbeat.\n");
		}
		else if ( type == TypeError ) {
			logit("et", "focalmech_ga:  Error sending error:%d.\n", ierr);
		}
	}

	return;
}

/*
 * focalmech_ga_end() - free all the local memory & close queue
 */
static void focalmech_ga_end( void )
{
	tport_detach( &Region );
	tpool_destroy( ThrdPool, 0 );
	tac_velmod_free();

	return;
}

/**
 * @brief
 *
 * @param arg
 */
static void proc_event( EARLY_EVENT_MSG *evt_msg )
{
	volatile _PICKS_TAC_ARG tac_arg[ThreadsNum];
/* */
	FPL_OBSERVE *obs;
	FPL_RESULT   best_solution;
	FPL_RESULT   best_sdv;
	FPL_RESULT   dbresult[2];
	FPL_RESULT   ptaxis[2];
	int          nobs;
	double       f_score;
	double       quality;
	char         output_dir[MAX_PATH_STR];
	char         fullpath[MAX_PATH_STR];
	char         command[MAX_PATH_STR * 4];
	char         script_args[MAX_PATH_STR * 2];

/* Tell the user we're working */
	logit("ot", "focalmech_ga: Receive a new event message (%s), start to process it!\n", evt_msg->header.event_id);
/* Get the take-off angle & azimuth of all the picks */
	for ( int i = 0; i < ThreadsNum; i++ ) {
		tac_arg[i].buffer = evt_msg;
		tac_arg[i].proc_thrds = ThreadsNum;
		tac_arg[i].assign_idx = i;
		tac_arg[i].is_finish = ATOMIC_VAR_INIT(0);
		tpool_add_work( ThrdPool, proc_picks_tac, (void *)&tac_arg[i] );
	}
/* Generate the output dir & the output full path */
	if ( !mk_outdir_by_evt( output_dir, ReportPath, evt_msg ) )
		return;
	gen_focalplot_fullpath( fullpath, output_dir, evt_msg );
/* Waiting for the threads */
	for ( int i = 0; i < ThreadsNum; i++ )
		while ( !tac_arg[i].is_finish )
			sleep_ew(1);
/* Check for the number of observation */
	if ( (nobs = pack_picks_to_observes( &obs, evt_msg )) >= MinPickPolarity ) {
	/* Main computation */
		cal_focal_ga( &best_solution, &best_sdv, &f_score, &quality, obs, nobs );
		fpl_dbcouple( &best_solution, dbresult, &ptaxis[FPLF_T_AXIS], &ptaxis[FPLF_P_AXIS] );
		plot_focal_result( dbresult, ptaxis, &best_sdv, obs, nobs, f_score, quality, evt_msg, fullpath );
	/* Post to the Facebook or other place by external script */
		gen_script_command_args( script_args, evt_msg, dbresult, ptaxis, f_score, quality, nobs );
		for ( int i = 0; i < NumPostScripts; i++ ) {
			if ( PostScripts[i].min_magnitude <= evt_msg->header.mag ) {
				logit("o", "focalmech_ga: Executing the script: '%s'\n", PostScripts[i].script);
				sprintf(command, "%s %s", PostScripts[i].script, script_args);
			/* Execute system command to post focal */
				if ( system(command) )
					logit("e", "focalmech_ga: Execute the script: '%s' error, please check it!\n", PostScripts[i].script);
				else
					logit("t", "focalmech_ga: Execute the script: '%s' success!\n", PostScripts[i].script);
			}
		}
	/* Remove the plotted focal */
		if ( RemoveSwitch ) {
			remove(fullpath);
			remove(output_dir);
		}
	/* Release the requested memory space of observations */
		free(obs);
	}
	else {
	/* Remove the output directory */
		remove(output_dir);
	}
/* */
	logit("ot", "focalmech_ga: Finish the processing of event message (%s).\n", evt_msg->header.event_id);

	return;
}

/**
 * @brief
 *
 * @param arg
 */
static void proc_picks_tac( void *arg )
{
	_PICKS_TAC_ARG         *tac_arg = arg;
	EARLY_EVENT_MSG_HEADER *header = (EARLY_EVENT_MSG_HEADER *)tac_arg->buffer;
	EARLY_PICK_MSG         *pick = (EARLY_PICK_MSG *)(header + 1);

/* */
	for ( int i = tac_arg->assign_idx; i < header->npicks; i += tac_arg->proc_thrds )
		cal_pick_tko_azi( &pick[i], header->evlat, header->evlon, header->evdepth );
/* */
	tac_arg->is_finish += 1;

	return;
}

/**
 * @brief
 *
 * @param pick
 * @param evlat
 * @param evlon
 * @param evdep
 */
static void cal_pick_tko_azi( EARLY_PICK_MSG *pick, const double evlat, const double evlon, const double evdep )
{
	tac_main( &pick->tko, &pick->azi, evlat, evlon, evdep, pick->latitude, pick->longitude, pick->elevation * -0.001 );
	return;
}

/**
 * @brief
 *
 * @param result
 * @param ev_msg
 * @return int
 */
static int pack_picks_to_observes( FPL_OBSERVE **result, EARLY_EVENT_MSG *evt_msg )
{
	EARLY_PICK_MSG *pick = evt_msg->picks;
	FPL_OBSERVE    *obs  = NULL;
	int             nobs;

/* */
	*result = NULL;
/* */
	if ( !evt_msg->header.npicks )
		return 0;
/* */
	obs = calloc(evt_msg->header.npicks, sizeof(FPL_OBSERVE));
	nobs = 0;
	for ( int i = 0; i < evt_msg->header.npicks; i++, pick++ ) {
		if ( pick->polarity == ' ' || pick->polarity == '?' )
			continue;
		obs[nobs].azimuth  = pick->azi;
		obs[nobs].takeoff  = pick->tko;
		obs[nobs].polarity = pick->polarity == 'U' ? 1.0 : pick->polarity == 'D' ? -1.0 : 0.0;
		obs[nobs] = fpl_observe_deg2rad( &obs[nobs] );
		nobs++;
	}
/* */
	if ( nobs )
		*result = obs;

	return nobs;
}

/**
 * @brief
 *
 * @param fpl
 * @param sdv
 * @param f_score
 * @param quality
 * @param obs
 * @param nobs
 * @return int
 */
static int cal_focal_ga( FPL_RESULT *fpl, FPL_RESULT *sdv, double *f_score, double *quality, FPL_OBSERVE *obs, const int nobs )
{
	FPL_RESULT _result[Population];
	FPL_RESULT _sdv[Population];
	int        nres;

/* */
	*f_score = fpl_find_ga( obs, nobs, IterationNum, Population, MutateBits, ReprodutionRate, MutateRate, _result, &nres );
	*quality = fpl_quality_cal( obs, nobs, *f_score );
	nres     = fpl_result_refine( _result, _sdv, nres );
/* */
	*fpl = _result[0];
	*sdv = _sdv[0];
	*quality /= nres;

	return nres;
}

/**
 * @brief
 *
 * @param dbresult
 * @param ptaxis
 * @param sdv
 * @param obs
 * @param nobs
 * @param f_score
 * @param quality
 * @param ev_msg
 * @param output_path
 */
static void plot_focal_result( FPL_RESULT *dbresult, FPL_RESULT *ptaxis, FPL_RESULT *sdv, FPL_OBSERVE *obs, const int nobs, const double f_score, const double quality, void *ev_msg, const char *output_path )
{
	EARLY_EVENT_MSG_HEADER *header = (EARLY_EVENT_MSG_HEADER *)ev_msg;
/* */
	int    year, mon, day, hour, min;
	double sec;

/* */
	split_time( &year, &mon, &day, &hour, &min, &sec, header->origin_time );
/* */
	pfocal_canva_open( output_path );
/* */
	pfocal_dbplane_plot( dbresult, PLOT_BEACHBALL_RADIUS );
	pfocal_pt_axis_plot( ptaxis, PLOT_BEACHBALL_RADIUS );
	pfocal_observe_plot( obs, nobs, PLOT_BEACHBALL_RADIUS );
	pfocal_eq_info_plot( year, mon, day, hour, min, sec, header->mag, header->evlat, header->evlon, header->evdepth );
	pfocal_plane_info_plot( dbresult, ptaxis, sdv, quality, f_score, PLOT_BEACHBALL_RADIUS );
/* */
	pfocal_canva_close();

	return;
}

/**
 * @brief
 *
 * @param buffer
 * @param ev_msg
 * @param dbresult
 * @param ptaxis
 * @param f_score
 * @param quality
 * @param nobs
 * @return char*
 */
static char *gen_script_command_args( char *buffer, void *ev_msg, FPL_RESULT dbresult[2], FPL_RESULT ptaxis[2], double f_score, double quality, int nobs )
{
	EARLY_EVENT_MSG_HEADER *evheader = (EARLY_EVENT_MSG_HEADER *)ev_msg;

/* Command arguments for executing script */
	sprintf(
		buffer, "%s %.2lf %.6lf %.6lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %.2lf %d",
		evheader->event_id, evheader->origin_time, evheader->evlat, evheader->evlon, evheader->evdepth, evheader->mag,
		dbresult[0].strike, dbresult[0].dip, dbresult[0].rake, dbresult[1].strike, dbresult[1].dip, dbresult[1].rake,
		ptaxis[FPLF_T_AXIS].strike, ptaxis[FPLF_T_AXIS].dip, ptaxis[FPLF_P_AXIS].strike, ptaxis[FPLF_P_AXIS].dip,
		(1.0 - f_score) * 0.5, quality, nobs
	);

	return buffer;
}

/**
 * @brief
 *
 * @param opath
 * @param parent_path
 * @param ev_msg
 * @return char*
 */
static char *mk_outdir_by_evt( char *opath, const char *parent_path, void *ev_msg )
{
	EARLY_EVENT_MSG_HEADER *evheader = (EARLY_EVENT_MSG_HEADER *)ev_msg;
/* If the previous thread is still alived, kill it! */
	sprintf(opath, "%s%s_%d/", parent_path, evheader->event_id, MyPid);
/* */
	if ( access(opath, F_OK) && mkdir(opath, S_IRWXU | S_IRGRP | S_IROTH) ) {
		logit("e", "focalmech_ga: Cannot make the new directory for event (%s), skip it!\n", evheader->event_id);
		return NULL;
	}

	return opath;
}

/**
 * @brief
 *
 * @param output
 * @param ev_msg
 * @return char*
 */
static char *gen_focalplot_fullpath( char *output, const char *parent_path, void *ev_msg )
{
	EARLY_EVENT_MSG_HEADER *evheader = (EARLY_EVENT_MSG_HEADER *)ev_msg;

/* Generate the output figure name */
	sprintf(output, "%s%s_%d_fm.png", parent_path, evheader->event_id, evheader->seq);

	return output;
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
 * @param spectime
 */
static void split_time( int *year, int *mon, int *day, int *hour, int *min, double *sec, const double spectime )
{
	struct tm    sptime;
	const time_t _timestamp = (time_t)spectime;

/* */
	gmtime_r(&_timestamp, &sptime);
	*year = sptime.tm_year + 1900;
	*mon  = sptime.tm_mon + 1;
	*day  = sptime.tm_mday;
	*hour = sptime.tm_hour;
	*min  = sptime.tm_min;
	*sec  = sptime.tm_sec + (spectime - _timestamp);

	return;
}

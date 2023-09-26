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
#include <time.h>
#include <math.h>
/* Earthworm environment header include */
#include <earthworm.h>
#include <kom.h>
#include <transport.h>
#include <lockfile.h>
/* Local header include */
#include <focalmech_ga.h>
#include <full_event_msg.h>
#include <fpl_func.h>
#include <plot_focal.h>
#include <tko_azi_cal.h>

/* Functions in this source file
 *******************************/
static void focalmech_ga_config( char * );
static void focalmech_ga_lookup( void );
static void focalmech_ga_status( unsigned char, short, char * );
static void focalmech_ga_end( void );  /* Free all the local memory & close socket */

static void split_time( int *, int *, int *, int *, int *, double *, const double );

static SHM_INFO Region;      /* shared memory region to use for i/o    */

#define BUF_SIZE  150000     /* define maximum size for an event msg   */
#define MAXLOGO   8

static MSG_LOGO Getlogo[MAXLOGO];   /* array for requesting module,type,instid */
static pid_t    MyPid;              /* for restarts by startstop               */

/* Thread things */
#define THREAD_STACK 8388608        /* 8388608 Byte = 8192 Kilobyte = 8 Megabyte */
#define THREAD_OFF    0				/* ComputePGA has not been started      */
#define THREAD_ALIVE  1				/* ComputePGA alive and well            */
#define THREAD_ERR   -1				/* ComputePGA encountered error quit    */

#define MAX_PLOT_SHAKEMAPS      8
#define MAX_EMAIL_RECIPIENTS    20
#define MAX_POST_SCRIPTS        5
#define DEFAULT_SUBJECT_PREFIX "EEWS"

static volatile int ProcessEventStatus = THREAD_OFF;
static volatile int MailProcessStatus  = THREAD_OFF;
static volatile int PlotMapStatus      = THREAD_OFF;
static volatile _Bool Finish = 0;

/* Things to read or derive from configuration file */
static char     RingName[MAX_RING_STR];		/* name of transport ring for i/o    */
static char     MyModName[MAX_MOD_STR];		/* speak as this module name/id      */
static uint8_t  LogSwitch;					/* 0 if no logfile should be written */
static uint64_t HeartBeatInterval;			/* seconds between heartbeats        */
static uint64_t QueueSize;					/* max messages in output circular buffer */
static uint8_t  RemoveSwitch = 0;
static uint64_t IssueInterval = 30;
static short    nLogo = 0;
static char     ReportPath[MAX_PATH_STR];
static char     SubjectPrefix[MAX_STR_SIZE];     /* defaults to EEWS, settable now */
static char     EmailProgram[MAX_PATH_STR];      /* Path to the email program */
static PLOTSMAP PlotShakeMaps[MAX_PLOT_SHAKEMAPS];
static uint8_t  NumPlotShakeMaps = 0;              /* Number of plotting shakemaps */
static POSCRIPT PostScripts[MAX_POST_SCRIPTS];
static uint8_t  NumPostScripts = 0;                /* Number of exec. scripts */
static char     LinkURLPrefix[MAX_PATH_STR];

/* Things to look up in the earthworm.h tables with getutil.c functions */
static long          RingKey;       /* key of transport ring for i/o     */
static unsigned char InstId;        /* local installation id             */
static unsigned char MyModId;       /* Module Id for this program        */
static unsigned char TypeHeartBeat;
static unsigned char TypeError;
static unsigned char TypeFullEvent;

/* Error messages used by focalmech_ga */
#define ERR_MISSMSG       0   /* message missed in transport ring       */
#define ERR_TOOBIG        1   /* retreived msg too large for buffer     */
#define ERR_NOTRACK       2   /* msg retreived; tracking limit exceeded */
#define ERR_QUEUE         3   /* error queueing message for sending      */
static char  Text[150];       /* string for log/error messages          */

/*
 *
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

	uint8_t  buffer[FULL_EVENT_SIZE];
	FULL_EVENT_MSG_HEADER *fev_header = (FULL_EVENT_MSG_HEADER *)buffer;
#if defined( _V710 )
	ew_thread_t tid;            /* Thread ID */
#else
	unsigned    tid;            /* Thread ID */
#endif

/* Check command line arguments */
	if ( argc != 2 ) {
		fprintf( stderr, "Usage: focalmech_ga <configfile>\n" );
		exit( 0 );
	}
	Finish = 1;
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
	/* CreateSemaphore_ew(); */ /* Obsoleted by Earthworm */
	SemaPtr = CreateSpecificSemaphore_ew( 0 );
	psk_msgqueue_init( QueueSize, FULL_EVENT_SIZE + 1 );

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

		if ( ProcessEventStatus != THREAD_ALIVE ) {
			if ( StartThread(thread_proc_shake, (unsigned)THREAD_STACK, &tid) == -1 ) {
				logit("e", "focalmech_ga: Error starting ProcessShake thread; exiting!\n");
				tport_detach( &Region );
				focalmech_ga_end();
				exit(-1);
			}
			ProcessEventStatus = THREAD_ALIVE;
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
			res = tport_getmsg(&Region, Getlogo, nLogo, &reclogo, &recsize, (char *)buffer, FULL_EVENT_SIZE);
		/* No more new messages     */
			if ( res == GET_NONE ) {
				/* PostSemaphore(); */ /* Obsoleted by Earthworm */
				PostSpecificSemaphore_ew( SemaPtr );
				break;
			}
		/* Next message was too big */
			else if ( res == GET_TOOBIG ) {
			/* Complain and try again   */
				sprintf(
					Text, "Retrieved msg[%ld] (i%u m%u t%u) too big for Buffer[%ld]",
					recsize, reclogo.instid, reclogo.mod, reclogo.type, (long)FULL_EVENT_SIZE
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
			if ( reclogo.type == TypeFullEvent ) {
			/* Debug */
				/* logit( "", "%s", rec ); */
			/* Generate every given seconds or at the coda of event */
				res = psk_msgqueue_enqueue( buffer, recsize, reclogo );
				/* PostSemaphore(); */ /* Obsoleted by Earthworm */
				PostSpecificSemaphore_ew( SemaPtr );
				if ( res ) {
					if ( res == -2 ) {  /* Serious: quit */
					/* Currently, eneueue() in mem_circ_queue.c never returns this error. */
						sprintf(Text, "internal queue error. Terminating.");
						focalmech_ga_status( TypeError, ERR_QUEUE, Text );
						focalmech_ga_end();
						exit(-1);
					}
					else if ( res == -1 ) {
						sprintf(Text, "queue cannot allocate memory. Lost message.");
						focalmech_ga_status( TypeError, ERR_QUEUE, Text );
					}
					else if ( res == -3 ) {
					/*
					 * Queue is lapped too often to be logged to screen.
					 * Log circular queue laps to logfile.
					 * Maybe queue laps should not be logged at all.
					 */
						logit("et", "focalmech_ga: Circular queue lapped. Messages lost!\n");
					}
				}
			}
		} while ( 1 );
		sleep_ew(50);  /* no more messages; wait for new ones to arrive */
	}  /* while( 1 ) */
/*-----------------------------end of main loop-------------------------------*/
exit_procedure:
	Finish = 0;
	/* PostSemaphore(); */ /* Obsoleted by Earthworm */
	PostSpecificSemaphore_ew( SemaPtr );
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
	char  plfilename[MAX_PATH_STR];

/* Set to zero one init flag for each required command */
	ncommand = 9;
	for ( i = 0; i < ncommand; i++ )
		init[i] = 0;

	strcpy(EmailProgram, "\0");
	strcpy(SubjectPrefix, DEFAULT_SUBJECT_PREFIX);
	strcpy(LinkURLPrefix, "/");

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
		/* 4 */
			else if ( k_its("QueueSize") ) {
				QueueSize = k_long();
				init[4] = 1;
			}
			else if ( k_its("RemoveShakeMap") ) {
				RemoveSwitch = k_int();
				logit(
					"o", "focalmech_ga: This process will %s the files after posted.\n",
					RemoveSwitch ? "remove" : "keep"
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
				init[5] = 1;
			}
			else if ( k_its("IssueInterval") ) {
				IssueInterval = k_long();
				logit("o", "focalmech_ga: Real shakemap alarm interval change to %ld\n", IssueInterval);
			}
			else if ( k_its("NormalPolyLineFile") ) {
				str = k_str();
				if ( str )
					strcpy(plfilename, str);
				logit("o", "focalmech_ga: Normal polygon line file: %s\n", plfilename);

				if ( psk_plot_polyline_read( plfilename, PLOT_NORMAL_POLY ) ) {
					logit("e", "focalmech_ga: Error reading normal polygon line file; exiting!\n");
					exit(-1);
				}
				else {
					logit("o", "focalmech_ga: Reading normal polygon line file finish!\n");
				}
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
						logit("e", "trace2peak: Invalid installation name <%s>", str);
						logit("e", " in <GetEventsFrom> cmd; exiting!\n");
						exit(-1);
					}
				}
				if ( (str = k_str()) ) {
					if ( strcmp(str, "MOD_WILDCARD") == 0 ) {
						logit("e", "focalmech_ga: This module do not accept MOD_WILDCARD; exiting!\n");
						exit(-1);
					}
					if ( GetModId(str, &Getlogo[nLogo].mod) != 0 ) {
						logit("e", "trace2peak: Invalid module name <%s>", str);
						logit("e", " in <GetEventsFrom> cmd; exiting!\n");
						exit(-1);
					}
				}
				if ( (str = k_str()) ) {
					if ( GetType(str, &Getlogo[nLogo].type) != 0 ) {
						logit("e", "trace2peak: Invalid message type name <%s>", str);
						logit("e", " in <GetEventsFrom> cmd; exiting!\n");
						exit(-1);
					}
				}
				nLogo++;
				init[8] = 1;
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
		if ( !init[4] ) logit("e", "<QueueSize> "         );
		if ( !init[5] ) logit("e", "<ReportPath> "        );
		if ( !init[8] ) logit("e", "any <GetEventsFrom> " );
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
	if ( GetType( "TYPE_FULL_EVENT", &TypeFullEvent ) != 0 ) {
		fprintf(stderr, "focalmech_ga: Invalid message type <TYPE_FULL_EVENT>; exiting!\n");
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
	psk_plot_end();
	tac_velmod_free();
	psk_msgqueue_end();
	/* DestroySemaphore(); */ /* Obsoleted by Earthworm */
	//DestroySpecificSemaphore_ew( SemaPtr );

	return;
}

/***/
static thr_ret thread_proc_event( void *dummy )
{
	int      res;
	long     recsize;                        /* Size of retrieved message from queue */
	MSG_LOGO reclogo;                        /* logo of retrieved message     */
	char     command[MAX_PATH_STR * 4];
	char     script_args[MAX_PATH_STR * 2];
	uint8_t  buffer[FULL_EVENT_SIZE];
	FULL_EVENT_MSG_HEADER *fev_header = (FULL_EVENT_MSG_HEADER *)buffer;
	FPL_OBSERVE *obs    = NULL;
	int          nobs   = 0;
	double       f_score;
	double       quality;
	FPL_RESULT  best_solution;
	FPL_RESULT  best_sdv;
/* Tell the main thread we're ok */
	ProcessEventStatus = THREAD_ALIVE;
/* Initialization */

	do {
		//WaitSpecificSemaphore_ew(SemaPtr);

		//res = psk_msgqueue_dequeue( gmtmp, &recsize, &reclogo );

		if ( res == 0 ) {
		/* Get the take-off angle & azimuth of all the picks */
			FULL_EVENT_PICK_MSG *pick = (FULL_EVENT_PICK_MSG *)(fev_header + 1);
			for ( int i = 0; i < fev_header->npicks; i++, pick++ )
				cal_pick_tko_azi( pick, fev_header->evlat, fev_header->evlon, fev_header->evdepth );
		/* Generate the calculating function's threads */
			//if ( StartThread( thread_plot_shakemap, (unsigned)THREAD_STACK, &tid[0] ) == -1 ) {
			//	logit("e", "focalmech_ga: Error starting PlotShakemap_thr thread!\n");
			//}
			//PlotMapStatus = THREAD_ALIVE;
			//while ( PlotMapStatus != THREAD_OFF ) {
			//	printf("focalmech_ga: Waiting for the shakemap plotting...\n");
			//	sleep_ew(200);
			//}
		/* */
			if ( (nobs = pack_picks_to_observes( &obs, buffer )) < 10 )
				goto end_event;
		/* */
			cal_focal_ga( &best_solution, &best_sdv, &f_score, &quality, obs, nobs, 20, 800, 3, 0.032, 0.72 );
			plot_focal_result( &best_solution, &best_sdv, obs, nobs, f_score, quality, buffer, ReportPath );
		/* Post to the Facebook or other place by external script */
		 	gen_script_command_args( script_args, PlotShakeMaps );
			for ( int i = 0; i < NumPostScripts; i++ ) {
				if ( PostScripts[i].min_magnitude <= fev_header->magnitude ) {
					logit("o", "focalmech_ga: Executing the script: '%s'\n", PostScripts[i].script);
					sprintf(command, "%s %s", PostScripts[i].script, script_args);
				/* Execute system command to post focal */
					if ( system(command) )
						logit("e", "focalmech_ga: Execute the script: '%s' error, please check it!\n", PostScripts[i].script);
					else
						logit("t", "focalmech_ga: Execute the script: '%s' success!\n", PostScripts[i].script);
				}
			}
		/* Remove the plotted shakemaps */
			if ( RemoveSwitch ) {
				psk_misc_smfilename_gen( PlotShakeMaps, script_args, MAX_STR_SIZE );
				remove_shakemap( ReportPath, script_args );
			}
		}
end_event:
		free(obs);
	} while ( Finish );

/* we're quitting */
	ProcessEventStatus = THREAD_ERR;   /* file a complaint to the main thread */
	KillSelfThread();                  /* main thread will restart us */

	return NULL;
}

/*
*/
static thr_ret thread_plot_shakemap( void *dummy )
{
	int  i;
	char resfilename[MAX_STR_SIZE];

	PlotMapStatus = THREAD_ALIVE;

	for ( i = 0; i < NumPlotShakeMaps; i++ ) {
	/* Generate the result file name */
		psk_misc_smfilename_gen( PlotShakeMaps + i, resfilename, MAX_STR_SIZE );
	/* Plot the shakemap by the plotting functions like PGPLOT or GMT etc. */
		if ( psk_plot_sm_plot( PlotShakeMaps + i, ReportPath, resfilename ) < 0 ) {
			logit( "e", "focalmech_ga: Plotting shakemap error; skip it!\n" );
		}
	}

/* we're quitting
 *****************/
	PlotMapStatus = THREAD_OFF;  /* fire a complaint to the main thread */
	KillSelfThread();
	return NULL;
}

/**
 * @brief
 *
 * @param pick
 * @param evlat
 * @param evlon
 * @param evdep
 */
static void cal_pick_tko_azi( FULL_EVENT_PICK_MSG *pick, const double evlat, const double evlon, const double evdep )
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
static int pack_picks_to_observes( FPL_OBSERVE **result, void *ev_msg )
{
	FULL_EVENT_MSG_HEADER *header = (FULL_EVENT_MSG_HEADER *)ev_msg;
	FULL_EVENT_PICK_MSG   *pick   = (FULL_EVENT_PICK_MSG *)(header + 1);
	FPL_OBSERVE           *obs    = NULL;
	int                    nobs;

/* */
	*result = NULL;
/* */
	if ( !header->npicks )
		return 0;
/* */
	obs = calloc(header->npicks, sizeof(FPL_OBSERVE));
	nobs = 0;
	for ( int i = 0; i < header->npicks; i++, pick++ ) {
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
 * @param nitr
 * @param npop
 * @param mutate_bits
 * @param repro_rate
 * @param mutate_rate
 * @return int
 */
static int cal_focal_ga( FPL_RESULT *fpl, FPL_RESULT *sdv, double *f_score, double *quality, FPL_OBSERVE *obs, const int nobs, int nitr, int npop, int mutate_bits, double repro_rate, double mutate_rate )
{
	FPL_RESULT _result[npop];
	FPL_RESULT _sdv[npop];
	int        nres;

/* */
	*f_score = fpl_find_ga( obs, nobs, nitr, npop, mutate_bits, repro_rate, mutate_rate, _result, &nres );
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
 * @param fpl
 * @param sdv
 * @param obs
 * @param nobs
 * @param f_score
 * @param quality
 * @param ev_msg
 * @param output_path
 */
static void plot_focal_result( FPL_RESULT *fpl, FPL_RESULT *sdv, FPL_OBSERVE *obs, const int nobs, const double f_score,  const double quality, void *ev_msg, const char *output_path )
{
	FULL_EVENT_MSG_HEADER *header = (FULL_EVENT_MSG_HEADER *)ev_msg;
	FPL_RESULT dbresult[2];
	FPL_RESULT ptaxis[2];
/* */
	int    year, mon, day, hour, min;
	double sec;
	char   path[MAX_PATH_STR];

/* */
	fpl_dbcouple( fpl, dbresult, &ptaxis[FPLF_T_AXIS], &ptaxis[FPLF_P_AXIS] );
	split_time( &year, &mon, &day, &hour, &min, &sec, header->origin_time );
/* Generate the output figure name */
	sprintf(path, "%s%s_focal.png", output_path, header->event_id);
	pfocal_canva_open( path );
/* */
	pfocal_dbplane_plot( dbresult, PLOT_BEACHBALL_RADIUS );
	pfocal_pt_axis_plot( ptaxis, PLOT_BEACHBALL_RADIUS );
	pfocal_observe_plot( obs, nobs, PLOT_BEACHBALL_RADIUS );
	pfocal_eq_info_plot( year, mon, day, hour, min, sec, header->magnitude, header->evlat, header->evlon, header->evdepth );
	pfocal_plane_info_plot( dbresult, ptaxis, sdv->strike, sdv->dip, sdv->rake, quality, f_score, PLOT_BEACHBALL_RADIUS );
/* */
	pfocal_canva_close();

	return;
}

/*
*/
static void remove_shakemap( const char *reportpath, const char *resfilename )
{
	char fullfilepath[MAX_PATH_STR*2];

	sprintf(fullfilepath, "%s%s", reportpath, resfilename);
	remove(fullfilepath);

	return;
}

/*
 *
 */
static char *gen_script_command_args( char *buffer, const PLOTSMAP *psm )
{
	int             i;
	GRIDMAP_HEADER *gmref   = psk_misc_refmap_get( psm );
	double          max_mag = DUMMY_MAG;
	struct tm      *tp      = NULL;
	char            filename[MAX_PATH_STR];
	char            starttime[MAX_DSTR_LENGTH];
	char            endtime[MAX_DSTR_LENGTH];
	time_t          reptime;

/* Generate the timestamp */
	tp = localtime(&gmref->starttime);
	date2spstring( tp, starttime, MAX_DSTR_LENGTH );
	tp = localtime(&gmref->endtime);
	date2spstring( tp, endtime, MAX_DSTR_LENGTH );
	reptime = gmref->codaflag ? -1 : gmref->endtime - gmref->starttime;
/* And get the maximum magnitude */
	for ( i = 0; i < 4; i++ )
		if ( gmref->magnitude[i] > max_mag )
			max_mag = gmref->magnitude[i];

/* Command arguments for executing script */
	sprintf(
		buffer, "%s %s %ld %f %d ", starttime, endtime, reptime, max_mag, psk_misc_trigstations_get( psm )
	);
/* Generate the result file names */
	for ( i = 0; i < NumPlotShakeMaps; i++ ) {
		psk_misc_smfilename_gen( psm + i, filename, MAX_STR_SIZE );
		strcat(buffer, ReportPath);
		strcat(buffer, filename);
		strcat(buffer, " ");
	}

	return buffer;
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

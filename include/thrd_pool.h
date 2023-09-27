/**
 * @file thrpool.h
 * @author your name (you@domain.com)
 * @brief
 * @version 0.1
 * @date 2023-09-27
 *
 * @copyright Copyright (c) 2023
 *
 */
#pragma once

#include <threads.h>

/**
 * @brief
 *
 */
typedef struct tpool_work {
	void (*process)( void * );
	void *arg;
	struct tpool_work *next; /* Link for next work */
} tpwork_t;

/**
 * @brief
 *
 */
typedef struct tpool {
	int       num_thread;        /* Number of threads */
	int       max_queue_size;    /* Max number of thread in the queue */

	thrd_t   *tid;               /* ID pointer for all thread */
	tpwork_t *queue;             /* The pointer for thread queue */
	int       front;             /* The front index for queue */
	int       rear;              /* The end index for queue */

	int       closed;            /* Close flag for putting, but stiil can run the existed threads */
	int       shutdown;          /* Shutdown flag for all thread works */

	mtx_t     queue_lock;        /* Mutex for the queue */
	cnd_t     queue_has_task;    /* Task condition */
	cnd_t     queue_has_space;   /* Space condition */
	cnd_t     queue_empty;       /* Empty condition */
} tpool_t;

/* */
int tpool_init( tpool_t **, int ,int );
int tpool_add_work( tpool_t *, void (*)( void * ), void * );
int tpool_destroy( tpool_t *, int );

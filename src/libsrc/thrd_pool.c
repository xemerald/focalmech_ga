/**
 * @file thrpool.c
 * @author Benjamin Ming Yang (b98204032u@gmail.com)
 * @brief
 * @version 0.1
 * @date 2023-09-27
 *
 * @copyright Copyright (c) 2023
 *
 */
#include <stdlib.h>
#include <stdio.h>
#include <threads.h>
/* */
#include <thrd_pool.h>

/* */
#define IS_TPOOL_EMPTY(_TPOOL) \
		((_TPOOL)->front == (_TPOOL)->rear)
/* */
#define IS_TPOOL_FULL(_TPOOL) \
		((((_TPOOL)->rear + 1) % (_TPOOL)->max_queue_size) == (_TPOOL)->front)
/* */
#define GET_TPOOL_SIZE(_TPOOL) \
		(((_TPOOL)->rear + (_TPOOL)->max_queue_size - (_TPOOL)->front) % (_TPOOL)->max_queue_size)

/* Internal pseudo working thread */
static int thread_worker( void * );

/**
 * @brief
 *
 * @param tpool
 * @param num_thread
 * @param max_queue_size
 * @return int
 */
int tpool_init( tpool_t **tpool, int num_thread ,int max_queue_size )
{
	tpool_t *pool;

/* Memory allocation of thread pool */
	if ( (pool = (tpool_t *)malloc(sizeof(tpool_t))) == NULL ) {
		fprintf(stderr, "tpool_init: ERROR in memory allocation!\n");
		return -1;
	}
/* Initialization of thread pool parameter */
	pool->num_thread = num_thread;
	pool->max_queue_size = max_queue_size + 1;
	pool->tid = NULL;
	pool->front = pool->rear = 0;
	pool->closed = pool->shutdown = 0;
/* */
	if ( mtx_init(&pool->queue_lock, mtx_plain) != thrd_success ) {
		fprintf(stderr, "tpool_init: ERROR in mutex initialization!\n");
		free(pool);
		return -1;
	}
/* */
	if ( cnd_init(&pool->queue_has_task) != thrd_success ) {
		fprintf(stderr, "tpool_init: ERROR in condition(queue_has_task) initialization!\n");
		free(pool);
		return -1;
	}
	if ( cnd_init(&pool->queue_has_space) != thrd_success ) {
		fprintf(stderr, "tpool_init: ERROR in condition(queue_has_space) initialization!\n");
		free(pool);
		return -1;
	}
	if ( cnd_init(&pool->queue_empty) != thrd_success ) {
		fprintf(stderr, "tpool_init: ERROR in condition(queue_empty) initialization!\n");
		free(pool);
		return -1;
	}
/* */
	if ( (pool->queue = malloc(sizeof(tpwork_t) * pool->max_queue_size)) == NULL ) {
		fprintf(stderr, "tpool_init: ERROR in pool queue memory allocation!\n");
		free(pool);
		return -1;
	}
	if ( (pool->tid = malloc(sizeof(thrd_t) * num_thread)) == NULL ) {
		fprintf(stderr, "tpool_init: ERROR in thread id memory allocation!\n");
		free(pool->queue);
		free(pool);
		return -1;
	}
/* Working threads creation */
	for ( int i = 0; i < num_thread; i++ ) {
		if ( thrd_create(&pool->tid[i], thread_worker, (void*)pool) != thrd_success ) {
			fprintf(stderr, "tpool_init: ERROR in threads creation!\n");
			free(pool->queue);
			free(pool);
			return -1;
		}
	}
/* Pass the pointer to caller's function */
	*tpool = pool;
/* Return 0 if initialization is succeed*/
	return 0;
}

/**
 * @brief
 *
 * @param tpool
 * @param process
 * @param arg
 * @return int
 */
int tpool_add_work( tpool_t *tpool, void (*process)( void * ), void *arg )
{
	tpwork_t *temp;
	int is_empty;

/* */
	mtx_lock(&tpool->queue_lock);
/* If queue is full, we just wait it */
	while ( IS_TPOOL_FULL( tpool ) && !tpool->shutdown && !tpool->closed )
		cnd_wait(&tpool->queue_has_space, &tpool->queue_lock);

	if ( tpool->shutdown || tpool->closed ) {
		mtx_unlock(&tpool->queue_lock);
		return -1;
	}

	is_empty = IS_TPOOL_EMPTY( tpool );
	/* int size = get_size(tpool); */ /* Debug */
/* Add a new node to queue */
	temp = tpool->queue + tpool->rear;
	temp->process = process;
	temp->arg = arg;
	tpool->rear = (tpool->rear + 1) % tpool->max_queue_size;
/* Broadcast to all thread that queue 'was' empty */
	if ( is_empty )
		cnd_signal(&tpool->queue_has_task);
	//if ( is_empty )
		//cnd_broadcast(&tpool->queue_has_task);

	mtx_unlock(&tpool->queue_lock);

	return 0;
}

/**
 * @brief
 *
 * @param tpool
 * @param finish
 * @return int
 */
int tpool_destroy( tpool_t *tpool, int finish )
{
	mtx_lock(&tpool->queue_lock);
/* Waiting for all works */
	tpool->closed = 1;
/* If finish equal 0 means hard kill, or need to be waited for all works */
	if ( finish ) {
		while( !IS_TPOOL_EMPTY( tpool ) )
			cnd_wait(&tpool->queue_empty, &tpool->queue_lock);
	}
/* */
	tpool->shutdown = 1;
	mtx_unlock(&tpool->queue_lock);
/* Broadcast to all thread that pool is ready to destroy */
	cnd_broadcast(&tpool->queue_has_task);
/* Waiting for all threads exiting */
	for ( int i = 0; i < tpool->num_thread; i++ )
		thrd_join(tpool->tid[i], NULL);

/* Free the memory of thread pool */
	free(tpool->tid);
	free(tpool->queue);
	free(tpool);
	mtx_destroy(&tpool->queue_lock);
	cnd_destroy(&tpool->queue_has_task);
	cnd_destroy(&tpool->queue_has_space);
	cnd_destroy(&tpool->queue_empty);

	return 0;
}

/**
 * @brief
 *
 * @param arg
 * @return int
 */
static int thread_worker( void *arg )
{
	tpool_t  *pool = (tpool_t *)arg;
	tpwork_t *work;
	int       is_full;

/* Inf. loop waiting for works */
	while ( 1 ) {
		mtx_lock(&pool->queue_lock);
	/* If the queue is empty, waiting for new work */
		while( IS_TPOOL_EMPTY( pool ) && !pool->shutdown )
			cnd_wait(&pool->queue_has_task, &pool->queue_lock);
	/* Thread pool shutdown process */
		if ( pool->shutdown ) {
			mtx_unlock(&pool->queue_lock);
			break;
		}
	/* */
		is_full = IS_TPOOL_FULL( pool );
	/* Point to the incoming work */
		work = pool->queue + pool->front;
	/* Update the index of queue */
		pool->front = (pool->front + 1) % pool->max_queue_size;
	/* Broadcast to all thread that queue 'was' full */
		if ( is_full )
			cnd_broadcast(&pool->queue_has_space);
	/* Send a signal that queue is empty */
		if ( IS_TPOOL_EMPTY( pool ) )
			cnd_signal(&pool->queue_empty);

		mtx_unlock(&pool->queue_lock);
	/* */
		work->process( work->arg );
	}

	return 0;
}

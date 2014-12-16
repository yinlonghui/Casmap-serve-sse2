#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <poll.h>
#include <pthread.h>

#include "utils.h"
#include "alndb.h"
#include "posdb.h"

#define ALN_BUF_SIZE	(1024)
#define POS_BUF_SIZE 	(1024)

typedef struct {
	int data_size;
	int len;
	void* buffer;
	int head,tail;
	pthread_mutex_t mtx;
	pthread_cond_t r_cond,w_cond;
	int shutdown;
}queue_t;

typedef struct {
	queue_t *queue;
	seq_id_t last_uid;
	int last_sid;
}aln_pipe_rt;

typedef struct {
	queue_t *queue;
}pos_pipe_rt;

queue_t* queue_create(size_t n, size_t size)
{
	queue_t *pq=malloc(sizeof(queue_t));
	pq->buffer=malloc(n*size);
	pq->head=pq->tail=0;
	pq->len=n;
	pq->data_size=size;
	pq->shutdown=0;
	pthread_mutex_init(&pq->mtx,NULL);
	pthread_cond_init(&pq->r_cond,NULL);
	pthread_cond_init(&pq->w_cond,NULL);
	return pq;
}

void queue_destroy(queue_t *pq)
{
	pthread_mutex_destroy(&pq->mtx);
	pthread_cond_destroy(&pq->r_cond);
	pthread_cond_destroy(&pq->w_cond);
	free(pq->buffer);
	free(pq);
}

void queue_shutdown(queue_t *pq)
{
	pthread_mutex_lock(&pq->mtx);
	pq->shutdown=1;
	pthread_cond_signal(&pq->r_cond);
	pthread_cond_signal(&pq->w_cond);
	pthread_mutex_unlock(&pq->mtx);
}

int queue_length(queue_t *pq)
{
	int len;
	pthread_mutex_lock(&pq->mtx);
	len=(pq->head+pq->len-pq->tail)%pq->len;
	pthread_mutex_unlock(&pq->mtx);
	return len;
}

int queue_pop(queue_t *pq, void *buf)
{
	int next_tail;
	pthread_mutex_lock(&pq->mtx);
	next_tail=(pq->tail+1)%pq->len;
	while(pq->head==pq->tail){
		if(pq->shutdown){
			pthread_mutex_unlock(&pq->mtx);
			return -1;
		}
		pthread_cond_wait(&pq->r_cond,&pq->mtx);
	}
	memcpy(buf,pq->buffer+pq->tail*pq->data_size,pq->data_size);
	pq->tail=next_tail;
	pthread_cond_signal(&pq->w_cond);
	pthread_mutex_unlock(&pq->mtx);
	return 0;
}

int queue_push(queue_t *pq, void *buf)
{
	int next_head;
	pthread_mutex_lock(&pq->mtx);
	if(pq->shutdown){
		pthread_mutex_unlock(&pq->mtx);
		return -1;
	}
	next_head=(pq->head+1)%pq->len;
	while(next_head==pq->tail){
		pthread_cond_wait(&pq->w_cond,&pq->mtx);
		if(pq->shutdown){
			pthread_mutex_unlock(&pq->mtx);
			return -1;
		}
	}
	memcpy(pq->buffer+pq->head*pq->data_size,buf,pq->data_size);
	pq->head=next_head;
	pthread_cond_signal(&pq->r_cond);
	pthread_mutex_unlock(&pq->mtx);
	return 0;
}

alndb_t alndb_open(const char* prefix, int flags)
{
	aln_pipe_rt *rt=malloc(sizeof(aln_pipe_rt));
	rt->queue=queue_create(ALN_BUF_SIZE,sizeof(aln_t));
	rt->last_uid=SEQ_ID_INVALID;
	rt->last_sid=-1;
	return (alndb_t)rt;
}

void alndb_shutdown(alndb_t db)
{
	aln_pipe_rt *rt=db;
	queue_shutdown(rt->queue);
}

int alndb_close(alndb_t db)
{
	aln_pipe_rt *rt=db;
	queue_destroy(rt->queue);
	free(rt);
	return 0;
}

int alndb_put(alndb_t db, hit_t *hit, seq_id_t uid)
{
	aln_t aln;
	aln_pipe_rt *rt=(aln_pipe_rt*)db;

	if(uid==rt->last_uid){
		rt->last_sid++;
	}else{
		rt->last_uid=uid;
		rt->last_sid=0;
	}

	aln.id=MAKE_ALN_ID(uid,rt->last_sid);
	aln.a=hit->a;
	aln.r=hit->r;
	aln.n_mm=hit->n_mm;
	aln.n_gapo=hit->n_gapo;
	aln.n_gape=hit->n_gape;
	aln.sa=hit->k;
	aln.width=hit->l-hit->k+1;
	//fprintf(stderr,"ALN PUT uid:%d sa:%d width:%d\n",uid,aln.sa,aln.width);
	return queue_push(rt->queue,&aln);
}

int alndb_get(alndb_t db, aln_t *paln)
{
	aln_pipe_rt *rt=(aln_pipe_rt*)db;
	return queue_pop(rt->queue,paln);
}

posdb_t posdb_open(const char* prefix)
{
	pos_pipe_rt *rt=malloc(sizeof(pos_pipe_rt));
	rt->queue=queue_create(POS_BUF_SIZE,sizeof(pos_t));
	return (posdb_t)rt;
}

void posdb_shutdown(posdb_t db)
{
	pos_pipe_rt *rt=db;
	queue_shutdown(rt->queue);
}

int posdb_close(posdb_t db)
{
	pos_pipe_rt *rt=db;
	queue_destroy(rt->queue);
	free(rt);
	return 0;
}

int posdb_put(posdb_t db, pos_t *pos)
{
	pos_pipe_rt *rt=(pos_pipe_rt*)db;
	return queue_push(rt->queue,pos);
}

int posdb_get(posdb_t db, pos_t *pos)
{
	pos_pipe_rt *rt=(pos_pipe_rt*)db;
	return queue_pop(rt->queue,pos);
}

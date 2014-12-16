#ifndef _TYPES_H_
#define _TYPES_H_

#include <stdlib.h>
#include <stdint.h>

typedef uint32_t seq_sz_t;
#define SEQ_SZ_MAX	(~((seq_sz_t)0))

typedef off_t seq_id_t;
#define SEQ_ID_INVALID	(~((seq_id_t)0))

typedef struct _hit_t
{
	unsigned char a:1,r:1;
	unsigned char n_mm;
	unsigned char n_gapo;
	unsigned char n_gape;
	seq_sz_t k;
	seq_sz_t l;
	unsigned int score;
	struct _hit_t *next;
	struct _hit_t *prev;
} hit_t;

#endif

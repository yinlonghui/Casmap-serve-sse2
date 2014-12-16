#ifndef _SEQDB_H_
#define _SEQDB_H_
#include "types.h"

#define SEQDB_NEXT	SEQ_ID_INVALID

typedef struct _seq_t
{
    	seq_id_t uid; //unique id
	char *name;
	char *seq;
	char *qual;
	seq_sz_t length;
} seq_t;

typedef void* seqdb_t;

seqdb_t seqdb_open(const char *seq1_fn, const char *seq2_fn, int flags);
int seqdb_close(seqdb_t db);
seq_t* seqdb_get(seqdb_t db, seq_id_t id);
void seqdb_release(seqdb_t db, seq_t *seq);
size_t seqdb_size(seqdb_t db);
off_t seqdb_offset(seqdb_t db);
size_t seqdb_seq_len(seqdb_t db);

#endif

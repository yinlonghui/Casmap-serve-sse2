#ifndef _REFDB_H_
#define _REFDB_H_

#include "types.h"

#define REFDB_IN_MEM	(1<<24)

typedef void* refdb_t;

refdb_t refdb_open(const char* prefix, int flags);
int refdb_close(refdb_t db);
int refdb_seq_nt4(refdb_t db, char* buf, seq_sz_t pos, seq_sz_t len);
unsigned int refdb_idx(refdb_t db, seq_sz_t pos);
unsigned int refdb_max_idx(refdb_t db);
const char* refdb_name(refdb_t db, unsigned int idx);
seq_sz_t refdb_offset(refdb_t db, unsigned int idx);
seq_sz_t refdb_length(refdb_t db, unsigned int idx);
seq_sz_t refdb_amb(refdb_t db, unsigned int idx);
seq_sz_t refdb_seq_amb(refdb_t db, seq_sz_t pos, seq_sz_t len);
seq_sz_t refdb_ref_len(refdb_t db);
unsigned int refdb_seed(refdb_t db);
#endif

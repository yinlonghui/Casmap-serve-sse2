#ifndef _SADB_H_
#define _SADB_H_

#include "types.h"

#define SADB_IN_MEM	(1<<25)

typedef void* sadb_t;

sadb_t sadb_open(const char* prefix, int flags);
int sadb_close(sadb_t db);

seq_sz_t sadb_get(sadb_t db, int strand, seq_sz_t sa);

#endif

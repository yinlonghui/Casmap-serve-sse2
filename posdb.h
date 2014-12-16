#ifndef _POSDB_H_
#define _POSDB_H_

#include "types.h"

typedef void* posdb_t;

typedef struct _pos_t{
	seq_id_t id; /* uid:48,aln_n:16 */
	uint8_t r:1,a:1;
	uint8_t n_mm;
	uint8_t n_gapo;
	uint8_t n_gape;
	seq_sz_t pos;
}pos_t;

posdb_t posdb_open(const char* prefix);
int posdb_close(posdb_t db);
int posdb_get(posdb_t db, pos_t *pos);
const char* posdb_meta_info(posdb_t db);

#endif /* _POSDB_H_ */

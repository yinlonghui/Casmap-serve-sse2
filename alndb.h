#ifndef _ALNDB_H_
#define _ALNDB_H_

#include "types.h"

#define DEFAULT_SA_INTV 32*4*1024*1024

#define ALNDB_ERROR	(1<<8)

#define ALN_ID_BITS	16

typedef struct{
	seq_id_t id; /* uid:48,aln_id:16 */
	uint8_t r:1,a:1;
	uint8_t n_mm;
	uint8_t n_gapo;
	uint8_t n_gape;
	seq_sz_t sa;
	seq_sz_t width;
}aln_t;

typedef void* alndb_t;
alndb_t alndb_open(const char* prefix, int flags);
int alndb_close(alndb_t db);
int alndb_put(alndb_t db, hit_t *seq, seq_id_t uid);
void alndb_meta_append(alndb_t db, const char *str);
const char* alndb_meta_info(alndb_t db);

#define MAKE_ALN_ID(seq_id,hit_id) (((seq_id)<<16)|((hit_id)&0xffff))
#define SEQ_ID(id) ((id)>>16)
#define ALN_ID(id) ((id)&0xffff)

#endif /* _ALNDB_H_ */


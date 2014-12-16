#include"pos2sampe.h"
#include"seq.h"
#include"refdb.h"
#include"utils.h"
#include"seqdb.h"
#include"posdb.h"


typedef uint32_t bwa_cigar_t;

#define CIGAR_OP_SHIFT 14
#define CIGAR_LN_MASK 0x3fff

#define __cigar_op(__cigar) ((__cigar)>>CIGAR_OP_SHIFT)
#define __cigar_len(__cigar) ((__cigar)&CIGAR_LN_MASK)
#define __cigar_create(__op, __len) ((__op)<<CIGAR_OP_SHIFT | (__len))

int aln_global_core(unsigned char *seq1, int len1, unsigned char *seq2, int len2, stdaln_rt *rt);


uint32_t *aln_path2cigar32(const path_t *path, int path_len, int *n_cigar);


int  se2sam_info(context_t *ctxt,refdb_t refdb,sam_info_t *sam_info1,int n_multi);
bwa_cigar_t *bwa_aln_path2cigar(const path_t *path, int path_len, int *n_cigar);


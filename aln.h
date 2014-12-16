#ifndef _ALN_H_
#define _ALN_H_

#include <stdint.h>
#include "seqdb.h"
#include "alndb.h"

#define BUFFER_POOL

#define MAX_DEV_NUM	4

#define BWA_EXT_STATE_MATCH (0x0)
#define BWA_EXT_STATE_WIDTH (0x1)
#define BWA_EXT_STATE_SEED  (0x2)
#define BWA_EXT_STATE_END	(0x3)

#define BWA_STOP		(0x0)
#define BWA_ONE_HIT		(0x1)
#define BWA_AUTO		(0x2)
#define BWA_NON_STOP		(0x3)

#define SEQ_INIT		(0x0u)
#define SEQ_NORMAL_LOADED	(0x1u)
#define SEQ_NORMAL_DONE		(0x2u)
#define SEQ_REVERSE_LOADED	(0x10000u)
#define SEQ_REVERSE_DONE	(0x20000u)
#define SEQ_FINISHED		(SEQ_NORMAL_LOADED|SEQ_NORMAL_DONE|SEQ_REVERSE_LOADED|SEQ_REVERSE_DONE)

#define BWA_PARAM_CNT0_ID		0
#define BWA_PARAM_CNT1_ID		1
#define BWA_PARAM_CNT2_ID		2
#define BWA_PARAM_CNT3_ID		3
#define BWA_PARAM_BWT_LEN_ID		4
#define BWA_PARAM_PRIMARY_ID		5
#define BWA_PARAM_RPRIMARY_ID		6
#define BWA_PARAM_MAX_GAPO_ID		7
#define BWA_PARAM_MAX_GAPE_ID		8
#define BWA_PARAM_MAX_DEL_OCC		9
#define BWA_PARAM_S_MM_ID		10
#define BWA_PARAM_S_GAPO_ID		11
#define BWA_PARAM_S_GAPE_ID		12
#define BWA_PARAM_SEED_WIDTH_ID		13
#define BWA_PARAM_INDEL_HEAD_SKIP_ID	14
#define BWA_PARAM_INDEL_TAIL_SKIP_ID	15
#define BWA_PARAM_FAST_MODE_ID		32
#define BWA_PARAM_TRAINING_ID		33
#define BWA_PARAM_UPDATE_EN_ID		34
#define BWA_PARAM_RID_WIDTH_ID		35
#define BWA_PARAM_MAX_BACK_ID		36

struct package_header{
#if __BYTE_ORDER == __LITTLE_ENDIAN
	uint32_t dw:16,type:16;
#else
	uint32_t type:16,dw:16;
#endif
};

struct db_package{
	struct package_header header;
	uint32_t	address;
	uint32_t data[8];
};

struct param_package{
	struct package_header header;
	uint32_t id;
	uint32_t value;
};

struct reads_header{
	struct package_header header;
	uint32_t 	gid;
#if __BYTE_ORDER == __LITTLE_ENDIAN
	uint32_t	length:16,max_diff:16;
	uint32_t	seed_diff:16,
			reffort:2,rsv3:2,rmode:2,rsv2:1,ren:1,
			neffort:2,rsv1:2,nmode:2,rsv0:1,nen:1;
#else
	uint32_t 	max_diff:16,length:16;
	uint32_t 	nen:1,rsv0:1,nmode:2,rsv1:2,neffort:2,
			ren:1,rsv2:1,rmode:2,rsv3:2,reffort:2,
			seed_diff:16;
#endif
};

struct match_package{
	struct package_header header;
	uint32_t gid;
#if __BYTE_ORDER == __LITTLE_ENDIAN
	uint32_t	n_gape:8,n_gapo:8,n_mm:8,
			ext_state:2,last:1,hit:1,a:1,rsv1:3;
#else
	uint32_t	rsv1:3,a:1,hit:1,last:1,ext_state:2,
			n_mm:8,n_gapo:8,n_gape:8;
#endif
	uint32_t k;
	uint32_t l;
};

#define	PKG_TYPE_PARAM 0
#define	PKG_TYPE_DB    1
#define	PKG_TYPE_READ  2
#define	PKG_TYPE_MATCH 3

#define init_param_header(h)	do{(h)->type=0;(h)->dw=2;}while(0)
#define init_db_header(h)	do{(h)->type=1;(h)->dw=9;}while(0)
#define init_read_header(h,s)	do{(h)->type=2;(h)->dw=(s);}while(0)

#define GID_INVALID	(~((uint32_t)0))	

struct aln_param{
	uint32_t MAX_SEED_DIFF;
	uint32_t MAX_GAPO;
	uint32_t MAX_GAPE;
	uint32_t SEED_WIDTH;
	uint32_t MAX_DEL_OCC;
	uint32_t S_MM;
	uint32_t S_GAPO;
	uint32_t S_GAPE;
	uint32_t INDEL_END_SKIP;
	uint32_t EFFORT;
	uint32_t FAST_MODE;
	uint32_t RID_WIDTH;
	uint32_t MAX_BACK;
};

struct aln_opt{
	char *server[2];
	short port;
	char *seq1_prefix,*seq2_prefix;
	char *aln_prefix;
	int debug;
	int verbose;
	int quiet;
	struct aln_param aln_param;
	float fnr;
};

int cm_aln(seqdb_t seqdb, alndb_t alndb, struct aln_opt *opt);
void print_package(void *pkg);
int aln_cal_maxdiff(int l, double err, double thres);
#endif /* _ALN_H_ */

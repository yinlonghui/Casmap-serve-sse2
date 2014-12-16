#include <stdlib.h>
#include <stdio.h>
#include <fcntl.h>
#include <math.h>
#include <pthread.h>

#include "seq.h"
#include "utils.h"
#include "seqdb.h"
#include "posdb.h"
#include "alndb.h"
#include "refdb.h"

#include "pos2sam.h"

#define SAM_BUF_SIZE	8192

struct statistics{
	refdb_t refdb;
	int parts;
	unsigned long long slice;
	unsigned long long **stat;
	int *step;
	int **rate;
	int **update;
};

typedef struct
{
	int gap_open;
	int gap_ext;
	int gap_end;

	int *matrix;
	int row;
	int band_width;
} AlnParam;

typedef struct
{
	int i, j;
	unsigned char ctype;
} path_t;

typedef struct
{
	unsigned char Mt:3, It:2, Dt:2;
} dpcell_t;

typedef struct
{
	int M, I, D;
} dpscore_t;

typedef struct{
	int gap_open;
	int gap_ext;
	int gap_end;
	int *matrix;
	int row;
	int band_width;

	dpcell_t * cell_pool;
	unsigned int cell_pool_sz;
	dpcell_t **dpcell;
	unsigned int dpcell_sz;
	dpscore_t *score_pool;
	unsigned int score_pool_sz;
	path_t *path;
	unsigned int path_sz; 
	unsigned int path_len; 
} stdaln_rt;

#define FROM_M 0
#define FROM_I 1
#define FROM_D 2
#define FROM_S 3

#define MINOR_INF -1073741823

#define SET_INF(s) (s).M = (s).I = (s).D = MINOR_INF;
#define set_M(MM, cur, p, sc)							\
{														\
	if ((p)->M >= (p)->I) {								\
		if ((p)->M >= (p)->D) {							\
			(MM) = (p)->M + (sc); (cur)->Mt = FROM_M;	\
		} else {										\
			(MM) = (p)->D + (sc); (cur)->Mt = FROM_D;	\
		}												\
	} else {											\
		if ((p)->I > (p)->D) {							\
			(MM) = (p)->I + (sc); (cur)->Mt = FROM_I;	\
		} else {										\
			(MM) = (p)->D + (sc); (cur)->Mt = FROM_D;	\
		}												\
	}													\
}
#define set_I(II, cur, p)								\
{														\
	if ((p)->M - gap_open > (p)->I) {					\
		(cur)->It = FROM_M;								\
		(II) = (p)->M - gap_open - gap_ext;				\
	} else {											\
		(cur)->It = FROM_I;								\
		(II) = (p)->I - gap_ext;						\
	}													\
}
#define set_end_I(II, cur, p)							\
{														\
	if (gap_end >= 0) {									\
		if ((p)->M - gap_open > (p)->I) {				\
			(cur)->It = FROM_M;							\
			(II) = (p)->M - gap_open - gap_end;			\
		} else {										\
			(cur)->It = FROM_I;							\
			(II) = (p)->I - gap_end;					\
		}												\
	} else set_I(II, cur, p);							\
}
#define set_D(DD, cur, p)								\
{														\
	if ((p)->M - gap_open > (p)->D) {					\
		(cur)->Dt = FROM_M;								\
		(DD) = (p)->M - gap_open - gap_ext;				\
	} else {											\
		(cur)->Dt = FROM_D;								\
		(DD) = (p)->D - gap_ext;						\
	}													\
}
#define set_end_D(DD, cur, p)							\
{														\
	if (gap_end >= 0) {									\
		if ((p)->M - gap_open > (p)->D) {				\
			(cur)->Dt = FROM_M;							\
			(DD) = (p)->M - gap_open - gap_end;			\
		} else {										\
			(cur)->Dt = FROM_D;							\
			(DD) = (p)->D - gap_end;					\
		}												\
	} else set_D(DD, cur, p);							\
}

typedef struct {
	seq_t *seq;
	char *nt4[2];
	unsigned int nt4_sz[2];
	char *seq_comp;
	unsigned int seq_comp_sz;
	char *qual_comp;
	unsigned int qual_comp_sz;
	int valid[2];
	char *ref_nt4;
	unsigned int ref_nt4_sz;;
	unsigned int ref_len;
	
	pos_t *pos;
	unsigned int pos_n;
	unsigned int pos_sz;

	unsigned int max_diff;
	unsigned int c0;
	unsigned int c1;

	/* for stdaln */
	stdaln_rt stdaln_rt;
	char *sam_buf;
	int done;
}context_t;

struct mt_rt{
	context_t *context;
	refdb_t refdb;
	seqdb_t seqdb;
	FILE* fout;
	int multi;
	int p0,p1,p2;
	pthread_mutex_t m0,m1,m2;
	pthread_cond_t c0,c1,c2;
	int destroy;
	int buf_len;
};

static int g_log_n[256];

static int aln_sm_maq[] = {
	11, -19, -19, -19, -13,
	-19, 11, -19, -19, -13,
	-19, -19, 11, -19, -13,
	-19, -19, -19, 11, -13,
	-13, -13, -13, -13, -13
};

static void init_statistics(struct statistics* st, refdb_t refdb)
{
	int i,max_len;
	st->refdb=refdb;
	st->parts=refdb_max_idx(refdb)+1;
	st->stat=calloc(st->parts,sizeof(unsigned long long*));
	st->rate=calloc(st->parts,sizeof(int*));
	st->update=calloc(st->parts,sizeof(int*));
	st->step=calloc(st->parts,sizeof(unsigned long long));
	max_len=0;
	for(i=0;i<st->parts;++i){
		int len=refdb_length(refdb,i);
		if(len>max_len)
			max_len=len;
	}
	st->slice=max_len/1024;
	for(i=0;i<st->parts;++i){
		st->step[i]=refdb_length(refdb,i)/st->slice+1;
		st->stat[i]=calloc(st->step[i],sizeof(unsigned long long));
		st->rate[i]=calloc(st->step[i],sizeof(int));
		st->update[i]=calloc(st->step[i],sizeof(int));
	}
	for(i=0;i<st->parts;++i){
		if(i==0){
			fprintf(stderr,"[pos2sam] [\"%s\":%u",refdb_name(refdb,i),st->step[i]);
		}else{
			fprintf(stderr,",\"%s\":%u",refdb_name(refdb,i),st->step[i]);
		}
	}
	fprintf(stderr,"]\n");
}

static void destroy_statistics(struct statistics *st)
{
	int i;
	for(i=0;i<st->parts;++i){
		free(st->stat[i]);
		free(st->rate[i]);
	}
	free(st->stat);
	free(st->rate);
	free(st->step);
}

static void update_statistics(struct statistics *st, seq_sz_t pos)
{
	int prev;
	int idx=refdb_idx(st->refdb,pos);
	seq_sz_t rel=pos-refdb_offset(st->refdb,idx);
	int p=rel/st->slice;
	st->stat[idx][p]++;
	prev=st->rate[idx][p];
	st->rate[idx][p]=st->stat[idx][p]*100/st->slice;
	if(st->rate[idx][p]!=prev)
		st->update[idx][p]=1;
}

static void report_statistics(struct statistics *st)
{
	unsigned int i,j,k;
	for(i=0,k=0;i<st->parts;++i){
		for(j=0;j<st->step[i] && st->update[i][j];++j){
			if(k==0){ 
				fprintf(stderr,"[pos2sam] {(%u,%u,%u)",i,j,st->rate[i][j]);
			}else{
				fprintf(stderr,",(%u,%u,%u)",i,j,st->rate[i][j]);
			}
			st->update[i][j]=0;
			++k;
		}
	}
	if(k) 
		fprintf(stderr,"}\n");
}

int aln_global_core(unsigned char *seq1, int len1, unsigned char *seq2, int len2, stdaln_rt *rt)
{
	register int i, j;
	dpcell_t **dpcell, *q;
	dpscore_t *curr, *last, *s;
	path_t *p;
	int b1, b2, tmp_end;
	int *mat, end, max;
	unsigned char type, ctype;

	int gap_open, gap_ext, gap_end, b;
	int *score_matrix, N_MATRIX_ROW;

	/* initialize some align-related parameters. just for compatibility */
	gap_open = rt->gap_open;
	gap_ext = rt->gap_ext;
	gap_end = rt->gap_end;
	b = rt->band_width;
	score_matrix = rt->matrix;
	N_MATRIX_ROW = rt->row;
	
	if (len1 == 0 || len2 == 0) {
		rt->path_len=0;
		return 0;
	}

	/* calculate b1 and b2 */
	if (len1 > len2) {
		b1 = len1 - len2 + b;
		b2 = b;
	} else {
		b1 = b;
		b2 = len2 - len1 + b;
	}
	if (b1 > len1) b1 = len1;
	if (b2 > len2) b2 = len2;
	--seq1; --seq2;

	/* allocate memory */
	end = (b1 + b2 <= len1)? (b1 + b2 + 1) : (len1 + 1);

	realloc_exp(rt->cell_pool,rt->cell_pool_sz,sizeof(dpcell_t)*end*(len2+1));
	realloc_exp(rt->dpcell,rt->dpcell_sz,sizeof(dpcell_t*)*(len2+1));
	realloc_exp(rt->score_pool,rt->score_pool_sz,sizeof(dpscore_t)*(len1+1)*2);
	realloc_exp(rt->path,rt->path_sz,sizeof(path_t)*(len1+len2+1));

	/* init pointers */
	curr=rt->score_pool;
	last=rt->score_pool+len1+1;
	p=rt->path;

	dpcell=rt->dpcell;
	for(j=0;j<=len2;++j)
		dpcell[j]=rt->cell_pool+j*end;
	for(j=b2+1;j<=len2;++j)
		dpcell[j]-=j-b2;

	/* set first row */
	SET_INF(*curr); curr->M = 0;
	for (i = 1, s = curr + 1; i < b1; ++i, ++s) {
		SET_INF(*s);
		set_end_D(s->D, dpcell[0] + i, s - 1);
	}
	s = curr; curr = last; last = s;

	/* core dynamic programming, part 1 */
	tmp_end = (b2 < len2)? b2 : len2 - 1;
	for (j = 1; j <= tmp_end; ++j) {
		q = dpcell[j]; s = curr; SET_INF(*s);
		set_end_I(s->I, q, last);
		end = (j + b1 <= len1 + 1)? (j + b1 - 1) : len1;
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		++s; ++q;
		for (i = 1; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]); /* this will change s->M ! */
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_D(s->D, q, s - 1);
		if (j + b1 - 1 > len1) { /* bug fixed, 040227 */
			set_end_I(s->I, q, last + i);
		} else s->I = MINOR_INF;
		s = curr; curr = last; last = s;
	}
	/* last row for part 1, use set_end_D() instead of set_D() */
	if (j == len2 && b2 != len2 - 1) {
		q = dpcell[j]; s = curr; SET_INF(*s);
		set_end_I(s->I, q, last);
		end = (j + b1 <= len1 + 1)? (j + b1 - 1) : len1;
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		++s; ++q;
		for (i = 1; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]); /* this will change s->M ! */
			set_I(s->I, q, last + i);
			set_end_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_end_D(s->D, q, s - 1);
		if (j + b1 - 1 > len1) { /* bug fixed, 040227 */
			set_end_I(s->I, q, last + i);
		} else s->I = MINOR_INF;
		s = curr; curr = last; last = s;
		++j;
	}

	/* core dynamic programming, part 2 */
	for (; j <= len2 - b2 + 1; ++j) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		end = j + b1 - 1;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i != end; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + i - 1, mat[seq1[i]]);
		set_D(s->D, q, s - 1);
		s->I = MINOR_INF;
		s = curr; curr = last; last = s;
	}

	/* core dynamic programming, part 3 */
	for (; j < len2; ++j) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i < len1; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + len1 - 1, mat[seq1[i]]);
		set_end_I(s->I, q, last + i);
		set_D(s->D, q, s - 1);
		s = curr; curr = last; last = s;
	}
	/* last row */
	if (j == len2) {
		SET_INF(curr[j - b2]);
		mat = score_matrix + seq2[j] * N_MATRIX_ROW;
		for (i = j - b2 + 1, q = dpcell[j] + i, s = curr + i; i < len1; ++i, ++s, ++q) {
			set_M(s->M, q, last + i - 1, mat[seq1[i]]);
			set_I(s->I, q, last + i);
			set_end_D(s->D, q, s - 1);
		}
		set_M(s->M, q, last + len1 - 1, mat[seq1[i]]);
		set_end_I(s->I, q, last + i);
		set_end_D(s->D, q, s - 1);
		s = curr; curr = last; last = s;
	}

	/* backtrace */
	i = len1; j = len2;
	q = dpcell[j] + i;
	s = last + len1;
	max = s->M; type = q->Mt; ctype = FROM_M;
	if (s->I > max) { max = s->I; type = q->It; ctype = FROM_I; }
	if (s->D > max) { max = s->D; type = q->Dt; ctype = FROM_D; }

	p = rt->path;
	p->ctype = ctype; p->i = i; p->j = j; /* bug fixed 040408 */
	++p;
	do {
		switch (ctype) {
			case FROM_M: --i; --j; break;
			case FROM_I: --j; break;
			case FROM_D: --i; break;
		}
		q = dpcell[j] + i;
		ctype = type;
		switch (type) {
			case FROM_M: type = q->Mt; break;
			case FROM_I: type = q->It; break;
			case FROM_D: type = q->Dt; break;
		}
		p->ctype = ctype; p->i = i; p->j = j;
		++p;
	} while (i || j);
	rt->path_len = p - rt->path - 1;

	return max;
}

extern int cm_aln_std(char* ref_nt4,int ref_len,char* seq_nt4, int seq_len, char *cigar, char *md, int *nm, stdaln_rt *rt) ;
#if 0
#define edit_char(m) (((m)==FROM_M)?'M':((m)==FROM_I)?'I':((m)==FROM_D)?'D':'S')
static int cm_aln_std(char* ref_nt4,int ref_len,char* seq_nt4,
		int seq_len, char *cigar, char *md, int *nm, stdaln_rt *rt)
{
	path_t *path;
	int path_len;
	unsigned char last_type;
	int i,j,k,m,n;
	char c;
	int _nm=0;
	int offset=0;

	/* standard global align */
	aln_global_core((unsigned char*)ref_nt4, ref_len, 
			(unsigned char*)seq_nt4, seq_len, 
			rt);
	path=rt->path;
	path_len=rt->path_len;

	/* trim off heading and tailing deletes */
	while(path->ctype==FROM_D){
		++path;
		--path_len;
	}
	while((path+path_len-1)->ctype==FROM_D){
		++offset;
		--path_len;
	}
	/* Print CIGAR and MD string */
	last_type=path[path_len-1].ctype;
	for(i=path_len-1,j=offset,k=0,m=0,n=0;i>=0;--i){
		switch(path[i].ctype){
			case FROM_M:
				if(ref_nt4[j]==seq_nt4[k]){
					++m;
				}else{
					md+=sprintf(md,"%u%c",m,nt4_to_a(ref_nt4[j]));
					m=0;
					++_nm;
				}
				++j;
				++k;
				break;
			case FROM_I:
				++k;
				if(last_type!=FROM_I) ++_nm;
				break;
			case FROM_D:
				if(m>0)
					md+=sprintf(md,"%u^%c",m,nt4_to_a(ref_nt4[j]));
				else
					md+=sprintf(md,"^%c",nt4_to_a(ref_nt4[j]));
				m=0;
				++j;
				if(last_type!=FROM_D) ++_nm;
				break;
		}
		if(path[i].ctype!=last_type){
			c=edit_char(last_type);
			cigar+=sprintf(cigar,"%u%c",n,c);
			n=1;
		}else{
			++n;
		}
		last_type=path[i].ctype;
	}
	if(m>0){
		sprintf(md,"%u",m);
	}
	if(n>0){
		c=edit_char(last_type);
		sprintf(cigar,"%u%c",n,c);
	}
	*nm=_nm;
	return offset;
}
#endif

static int cm_aln_simple(char* ref_nt4,int ref_len,char* seq_nt4,
		int seq_len, char *cigar, char *md, int *nm)
{
	int i;
	int m;
	int _nm=0;
	for(i=0,m=0;i<ref_len;++i){
		if(ref_nt4[i]==seq_nt4[i]){
			++m;
		}else{
			md+=sprintf(md,"%u%c",m,nt4_to_a(ref_nt4[i]));
			m=0;
			++_nm;
		}
	}
	if(m>0){
		sprintf(md,"%u",m);
	}
	sprintf(cigar,"%dM",ref_len);
	*nm=_nm;
	return 0;
}

static unsigned int approx_mapQ(const context_t *ctxt)
{
	unsigned int mapQ,n;
	unsigned c0=ctxt->c0;
	unsigned c1=ctxt->c1-ctxt->c0;

	if(c0==0) mapQ=23;
	else if(c0>1) mapQ=0;
	else if(ctxt->pos[0].n_mm==ctxt->max_diff) mapQ=25;
	else if(c1==0) mapQ=37;
	else {
		n=c1>=255?255:c1;
		mapQ=(23 < g_log_n[n])?0:23-g_log_n[n];
	}
	return mapQ;
}

static int cm_pos2sam_core(refdb_t refdb,context_t *ctxt,char *sam_buf,int multi)
{
	int rc;
	char *aln_nt4=NULL, *aln_seq=NULL, *aln_qual=NULL;
	unsigned long long last_id;
	char cigar[4096];
	char md[4096];
	uint32_t ref_pos;
	int ref_ext;
	pos_t *ppos;
	int i;
	const char *ref_name=NULL;
	int ref_offset;
	unsigned int ref_id;
	unsigned int ref_coor;
	unsigned int ref_amb=0;
	int rest,step;
	int repeat;
	int nm=0;
	int best,second = 0,score;
	//int c0,c1;

	if(ctxt->pos_n==0){
		sam_buf+=sprintf(sam_buf,"%s\t%u\t*\t0\t0\t*\t*\t0\t0\t%s\t%s",
				ctxt->seq->name,4,ctxt->seq->seq,ctxt->seq->qual);
		return 0;
	}

	best=ctxt->pos[0].n_mm*3+ctxt->pos[0].n_gapo*11+ctxt->pos[0].n_gape*4;
	for(i=1;i<ctxt->pos_n;++i){
		score=ctxt->pos[i].n_mm*3+ctxt->pos[i].n_gapo*11+ctxt->pos[i].n_gape*4;
		if(score>best){
			second=score;
			break;
		}
	}
	ctxt->c0=i;
	for(;i<ctxt->pos_n;++i){
		score=ctxt->pos[i].n_mm*3+ctxt->pos[i].n_gapo*11+ctxt->pos[i].n_gape*4;
		if(score>second){
			break;
		}
	}
	ctxt->c1=i;

	last_id=-1ll;
	ref_offset=1;
	ref_coor=0;
	ref_id=0;
	if(multi>0){
		++multi;
		rest=ctxt->c1>multi?multi:ctxt->c1;
		step=ctxt->c1/multi;
		if(step==0) step=1;
	}else{
		rest=1;
		step=1;
	}
	repeat=ctxt->c1>1?1:0;
	i=0;
	while(rest>0){
		ppos=ctxt->pos+i;

			//fprintf(stderr,"%llu\t%u\n",ppos->id>>16,ppos->pos);
		if(ppos->id !=last_id){
			ref_pos=ppos->pos;
			ref_ext=ppos->n_gapo+ppos->n_gape;

			if(ppos->r)
				ref_pos-=ctxt->seq->length;
			if(!ppos->a)
				ref_pos-=ref_ext;
			ctxt->ref_len=ctxt->seq->length+ref_ext;

			realloc_n(ctxt->ref_nt4,ctxt->ref_nt4_sz,ctxt->ref_len);
			rc=refdb_seq_nt4(refdb,ctxt->ref_nt4,ref_pos,ctxt->ref_len);
			if(rc){
				fprintf(stderr,"Can not retrieve reference with POS=%u\n",ref_pos);
				return -1;
			}


			if(ppos->a){
				if(ctxt->valid[1]==0){
					realloc_n(ctxt->seq_comp,ctxt->seq_comp_sz, ctxt->seq->length+1);
					realloc_n(ctxt->qual_comp,ctxt->qual_comp_sz, ctxt->seq->length+1);
					realloc_n(ctxt->nt4[1],ctxt->nt4_sz[1],ctxt->seq->length);
					seq_to_comp(ctxt->seq_comp,ctxt->seq->seq,ctxt->seq->length);
					ctxt->seq_comp[ctxt->seq->length]=0;
					reverse(ctxt->qual_comp,ctxt->seq->qual,ctxt->seq->length);
					ctxt->qual_comp[ctxt->seq->length]=0;
					seq_to_nt4(ctxt->nt4[1],ctxt->seq_comp,ctxt->seq->length);
					ctxt->valid[1]=1;
				}
				aln_nt4=ctxt->nt4[1];
				aln_seq=ctxt->seq_comp;
				aln_qual=ctxt->qual_comp;
			}else{
				if(ctxt->valid[0]==0){
					realloc_n(ctxt->nt4[0],ctxt->nt4_sz[0],ctxt->seq->length);
					seq_to_nt4(ctxt->nt4[0],ctxt->seq->seq,ctxt->seq->length);
					ctxt->valid[0]=1;
				}
				aln_nt4=ctxt->nt4[0];
				aln_seq=ctxt->seq->seq;
				aln_qual=ctxt->seq->qual;
			}

			if(ref_ext)
				ref_pos+=cm_aln_std(ctxt->ref_nt4,ctxt->ref_len,
						aln_nt4,ctxt->seq->length,
						cigar,md,&nm,&ctxt->stdaln_rt);
			else if(i==0)
				cm_aln_simple(ctxt->ref_nt4,ctxt->ref_len,aln_nt4,ctxt->seq->length,cigar,md,&nm);
			else{
				sprintf(cigar,"%dM",ctxt->ref_len);
				nm=ppos->n_mm;
			}

			ref_id=refdb_idx(refdb,ref_pos);
			ref_name=refdb_name(refdb,ref_id);
			ref_coor=refdb_offset(refdb,ref_id);
			ref_amb=refdb_seq_amb(refdb,ref_pos,ctxt->seq->length);

			/* adjust to start-from-1 index */
			//ref_pos++;
			ref_offset=ref_pos+1-ppos->pos;
		}

		ppos->pos=ppos->pos+ref_offset-ref_coor;

		if(i==0){
			sam_buf+=sprintf(sam_buf,"%s\t%u\t%s\t%u"
					"\t%u\t%s\t*\t0\t0"
					"\t%s\t%s"
					"\tXT:A:%c"
					"\tNM:i:%u",
					ctxt->seq->name,ppos->a?16:0,ref_name,ppos->pos,
					approx_mapQ(ctxt),cigar,
					aln_seq,aln_qual,
					repeat?'R':'U',
					nm
					);

			if(ref_amb)
				sam_buf+=sprintf(sam_buf,"\tXN:i:%u",ref_amb);
			sam_buf+=sprintf(sam_buf, "\tX0:i:%u\tX1:i:%u"
					"\tXM:i:%u\tXO:i:%u\tXG:i:%u"
					"\tMD:Z:%s",
					ctxt->c0,ctxt->c1-ctxt->c0,
					ppos->n_mm,ppos->n_gapo,ppos->n_gapo+ppos->n_gape,
					md
					);
			if(rest>1)
				sam_buf+=sprintf(sam_buf,"\tXA:Z:");
		}else{
			sam_buf+=sprintf(sam_buf,"%s,%c%u,%s,%u;",
					ref_name,ppos->a?'-':'+',ppos->pos,cigar,nm);

		}
		i+=step;
		--rest;
	}
	return 0;
}

static void init_context(context_t *ctxt)
{
	memset(ctxt,0,sizeof(context_t));
	ctxt->stdaln_rt.gap_open=26;
	ctxt->stdaln_rt.gap_ext=9;
	ctxt->stdaln_rt.gap_end=5;
	ctxt->stdaln_rt.matrix=aln_sm_maq;
	ctxt->stdaln_rt.row=5;
	ctxt->stdaln_rt.band_width=50;
	ctxt->sam_buf=malloc(SAM_BUF_SIZE);
}

static void destroy_context(context_t *ctxt)
{
	free(ctxt->nt4[0]);
	free(ctxt->nt4[1]);
	free(ctxt->seq_comp);
	free(ctxt->qual_comp);
	free(ctxt->ref_nt4);
	free(ctxt->pos);
	free(ctxt->stdaln_rt.cell_pool);
	free(ctxt->stdaln_rt.dpcell);
	free(ctxt->stdaln_rt.score_pool);
	free(ctxt->stdaln_rt.path);
	free(ctxt->sam_buf);
}

static void* pos2sam_worker(void *data)
{
	struct mt_rt *rts=data;
	while(1){
		context_t *ctxt;
		int curr;
		pthread_mutex_lock(&rts->m1);
		while(!rts->destroy && rts->p1==rts->p0){
			pthread_cond_wait(&rts->c1,&rts->m1);
		}

		if(rts->destroy){
			pthread_mutex_unlock(&rts->m1);
			break;
		}

		curr=rts->p1;
		ctxt=&rts->context[rts->p1];
		rts->p1=(rts->p1+1)%rts->buf_len;
		pthread_mutex_unlock(&rts->m1);

		/* Core */
		cm_pos2sam_core(rts->refdb,ctxt,ctxt->sam_buf,rts->multi);
		pthread_mutex_lock(&rts->m2);
		ctxt->done=1;
		if(curr==rts->p2)
			pthread_cond_signal(&rts->c2);
		pthread_mutex_unlock(&rts->m2);
	}
	return 0;
}

static void* output_worker(void *data)
{
	struct mt_rt *rts=data;
	while(1){
		context_t *ctxt;
		pthread_mutex_lock(&rts->m2);
		ctxt=&rts->context[rts->p2];
		while(!rts->destroy && !ctxt->done){
			pthread_cond_wait(&rts->c2,&rts->m2);
		}
		pthread_mutex_unlock(&rts->m2);

		if(rts->destroy)
			break;

		/* Output results */
		fprintf(rts->fout,"%s\n",ctxt->sam_buf);
		seqdb_release(rts->seqdb,ctxt->seq);
		ctxt->done=0;

		pthread_mutex_lock(&rts->m0);
		rts->p2=(rts->p2+1)%rts->buf_len;
		pthread_cond_signal(&rts->c0);
		pthread_mutex_unlock(&rts->m0);
	}
	return 0;
}

/* Important: No re-entrance */
static int make_context(context_t *ctxt, seqdb_t seqdb, posdb_t posdb)
{
	int rc;
	static pos_t pos;
	static seq_id_t pos_id=SEQ_ID_INVALID;
	ctxt->seq=seqdb_get(seqdb,SEQDB_NEXT);
	if(ctxt->seq==NULL)
		return -1;
	ctxt->valid[0]=ctxt->valid[1]=0;
	ctxt->pos_n=0;
	ctxt->max_diff=cal_max_diff(ctxt->seq->length);
	do{
		if(pos_id==SEQ_ID_INVALID){
			rc=posdb_get(posdb,&pos);
			if(rc==-1){
				break;
			}
			pos_id=SEQ_ID(pos.id);
		}

		/*
		   if(pos_id<ctxt.seq->uid){
		   fprintf(stderr,"pos file in wrong order uid=%lu pos_id=%lu\n",ctxt.seq->uid,pos_id);
		   pos_id=SEQ_ID_INVALID;
		   continue;
		   }
		 */

		if(pos_id==ctxt->seq->uid){
			/* filtering may be applied here */
			ctxt->pos_n++;
			if(ctxt->pos_sz<ctxt->pos_n){
				ctxt->pos=realloc(ctxt->pos,sizeof(pos_t)*ctxt->pos_n);
				ctxt->pos_sz=ctxt->pos_n;
			}
			ctxt->pos[ctxt->pos_n-1]=pos;

			pos_id=SEQ_ID_INVALID;
		}
	}while(pos_id==SEQ_ID_INVALID);
	return ctxt->pos_n;
}

static int cm_pos2sam_mt(struct pos2sam_opt *opt)
{
	int rc,i;
	struct mt_rt rts;
	pthread_t *threads=NULL;
	double t_last, t_curr;
	unsigned long long n_seq=0,n_pos=0;
	struct statistics st;

	if(opt->thread<1 || opt->buf_len<1){
		return -1;
	}else{
		fprintf(stderr,"[pos2sam] threads: %u buffer length: %u\n",opt->thread,opt->buf_len);
	}

	rts.refdb=opt->refdb;
	rts.seqdb=opt->seqdb;
	rts.multi=opt->multi;
	rts.fout=opt->fout;
	rts.buf_len=opt->buf_len;
	rts.p0=rts.p1=rts.p2=0;
	rts.context=calloc(rts.buf_len,sizeof(context_t));
	for(i=0;i<rts.buf_len;++i){
		init_context(&rts.context[i]);
	}
	rts.destroy=0;

	pthread_mutex_init(&rts.m0,NULL);
	pthread_mutex_init(&rts.m1,NULL);
	pthread_mutex_init(&rts.m2,NULL);
	pthread_cond_init(&rts.c0,NULL);
	pthread_cond_init(&rts.c1,NULL);
	pthread_cond_init(&rts.c2,NULL);

	threads=calloc(opt->thread+1,sizeof(pthread_t));
	for(i=0;i<opt->thread;++i){
		rc=pthread_create(threads+i,NULL,pos2sam_worker,&rts);
		if(rc){
			xerror("create thread");
			return -1;
		}
	}
	rc=pthread_create(threads+i,NULL,output_worker,&rts);
	if(rc){
		xerror("create thread");
		return -1;
	}

	if(opt->report_map) 
		init_statistics(&st,opt->refdb);

	t_last=t_curr=get_run_time();

	do{
		int rc;
		context_t *ctxt;

		/* Check available buffer for new context */
		pthread_mutex_lock(&rts.m0);
		while((rts.p0+1)%rts.buf_len==rts.p2){
			pthread_cond_wait(&rts.c0,&rts.m0);
		}
		pthread_mutex_unlock(&rts.m0);
		ctxt=&rts.context[rts.p0];

		/* Build a context for cm_pos2sam_core */
		rc=make_context(ctxt,opt->seqdb,opt->posdb);
		if(rc==-1)
			break; /* No more input. */

		n_pos+=rc;
		n_seq++;

		if(opt->report_map){
			for(i=0;i<ctxt->pos_n;++i){
				update_statistics(&st,ctxt->pos[i].pos);
			}
		}

		/* Invoke workers */
		pthread_mutex_lock(&rts.m1);
		rts.p0=(rts.p0+1)%rts.buf_len;
		pthread_cond_signal(&rts.c1);
		pthread_mutex_unlock(&rts.m1);

		t_curr=get_run_time();
		if(t_curr-t_last>1){
			fprintf(stderr,"[pos2sam] %.2f sec, %llu reads, %llu pos\n",t_curr,n_seq,n_pos);
			t_last=t_curr;
			if(opt->report_map)
				report_statistics(&st);
		}
	}while(1);

	/* Now wait and let worker threads do their job */
	pthread_mutex_lock(&rts.m0);
	while(rts.p2!=rts.p0){
		pthread_cond_broadcast(&rts.c1);
		pthread_cond_wait(&rts.c0,&rts.m0);
	}
	pthread_mutex_unlock(&rts.m0);

	rts.destroy=1;
	pthread_cond_broadcast(&rts.c1);
	pthread_cond_broadcast(&rts.c2);
	for(i=0;i<=opt->thread;++i){
		void *ret;
		pthread_join(threads[i],&ret);
	}

	t_curr=get_run_time();
	fprintf(stderr,"[pos2sam] %.2f sec, %llu reads, %llu pos\n",t_curr,n_seq,n_pos);
	if(opt->report_map){
		report_statistics(&st);
		destroy_statistics(&st);
	}

	for(i=0;i<rts.buf_len;++i){
		destroy_context(&rts.context[i]);
	}
	free(rts.context);
	free(threads);

	pthread_mutex_destroy(&rts.m0);
	pthread_mutex_destroy(&rts.m1);
	pthread_mutex_destroy(&rts.m2);
	pthread_cond_destroy(&rts.c0);
	pthread_cond_destroy(&rts.c1);
	pthread_cond_destroy(&rts.c2);

	return 0;
}

static int cm_pos2sam_st(struct pos2sam_opt *opt)
{
	int i;
	double t_last, t_curr;
	unsigned long long n_seq=0,n_pos=0;
	struct statistics st;
	context_t context;

	init_context(&context);

	if(opt->report_map) 
		init_statistics(&st,opt->refdb);

	t_last=t_curr=get_run_time();

	do{
		int rc;
		context_t *ctxt;

		ctxt=&context;

		/* Build a context for cm_pos2sam_core */
		rc=make_context(ctxt,opt->seqdb,opt->posdb);
		if(rc==-1)
			break; /* No more input. */

		n_pos+=rc;
		n_seq++;

		if(opt->report_map){
			for(i=0;i<ctxt->pos_n;++i){
				update_statistics(&st,ctxt->pos[i].pos);
			}
		}

		/* Core */
		cm_pos2sam_core(opt->refdb,ctxt,ctxt->sam_buf,opt->multi);

		/* Output */
		fprintf(opt->fout,"%s\n",ctxt->sam_buf);
		seqdb_release(opt->seqdb,ctxt->seq);

		t_curr=get_run_time();
		if(t_curr-t_last>1){
			fprintf(stderr,"[pos2sam] %.2f sec, %llu reads, %llu pos\n",t_curr,n_seq,n_pos);
			t_last=t_curr;
			if(opt->report_map)
				report_statistics(&st);
		}
	}while(1);

	t_curr=get_run_time();
	fprintf(stderr,"[pos2sam] %.2f sec, %llu reads, %llu pos\n",t_curr,n_seq,n_pos);
	if(opt->report_map){
		report_statistics(&st);
		destroy_statistics(&st);
	}

	destroy_context(&context);

	return 0;
}

int cm_pos2sam(struct pos2sam_opt *opt)
{
	int i,n_ref;
	for (i = 1; i != 256; ++i) 
		g_log_n[i] = (int)(4.343 * log(i) + 0.5);

	init_diff_table(opt->fnr);

	n_ref=refdb_max_idx(opt->refdb)+1;
	for(i=0;i<n_ref;++i){
		fprintf(opt->fout,"@SQ\tSN:%s\tLN:%u\n",refdb_name(opt->refdb,i),refdb_length(opt->refdb,i));
	}

	if(opt->thread)
		return cm_pos2sam_mt(opt);
	else
		return cm_pos2sam_st(opt);
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "global_aln.h"
#include "seq.h"
#include "alndb.h"
#include "refdb.h"
#include "utils.h"
#include "seqdb.h"
#include "posdb.h"
#include <pthread.h>
#include "ksort.h"
#include "ksw.h"


KSORT_INIT_GENERIC(uint64_t)

#define BUF_SEQ (0x40000)
#define OUTLIER_BOUND 2.0
#define	BUF 1024
#define THREAD_BLOCK_SIZE 1024

static pthread_mutex_t g_seq_lock = PTHREAD_MUTEX_INITIALIZER;

seq_id_t pos_id_0 = SEQ_ID_INVALID;
seq_id_t pos_id_1 = SEQ_ID_INVALID;
pos_t	 pos_0;
pos_t	 pos_1;
int	end = 0;
int	cnt_chg=0;
int g_log_n[256];
double t_last, t_curr;


typedef struct{
	int tid;
	int n_seq;
	context_t *ctxt[2];
	sam_info_t *sam_info[2];
	refdb_t refdb;
	uint64_t L;
	pe_data_t  *pe_data;
	isize_info_t ii;
	pe_opt_t *pe_opt;
	posdb_t posdb[2];
	seqdb_t seqdb[2];
}thread_data_t;



typedef struct{
	int tid;
	int n_seq;
	context_t *ctxt[2];
	sam_info_t *sam_info[2];
	refdb_t refdb;
	uint64_t L;
	pe_data_t  *pe_data;
	isize_info_t ii;
	pe_opt_t *pe_opt;
	posdb_t posdb;
	seqdb_t seqdb;
}pipe_data_t;

typedef struct{
	int n_seq;
	context_t *ctxt[2];
	sam_info_t *sam_info[2];
	refdb_t	refdb;
	FILE *fout;
	char *sam_buf;
}output_data_t;

static int aln_sm_maq[] = {
        11, -19, -19, -19, -13,
        -19, 11, -19, -19, -13,
        -19, -19, 11, -19, -13,
        -19, -19, -19, 11, -13,
        -13, -13, -13, -13, -13
};


pe_opt_t *bwa_init_pe_opt()
{
        pe_opt_t *po;
        po = (pe_opt_t*)calloc(1, sizeof(pe_opt_t));
        po->max_isize = 500;
        po->force_isize = 0;
        po->max_occ = 100000;
	po->n_thread =1;
        po->n_multi = 3;
        po->N_multi = 10;
        po->is_sw = 1;
        po->ap_prior = 1e-5;
        return po;
}

static inline uint64_t hash_64(uint64_t key)
{
	key += ~(key << 32);
	key ^= (key >> 22);
	key += ~(key << 13);
	key ^= (key >> 8);
	key += (key << 3);
	key ^= (key >> 15);
	key += ~(key << 27);
	key ^= (key >> 31);
	return key;
}


stdaln_rt *casmap_init_aln_rt(){
	stdaln_rt *aln_rt;
	aln_rt=(stdaln_rt *)calloc(1,sizeof(stdaln_rt));
	aln_rt->gap_open=26;
        aln_rt->gap_ext=9;
        aln_rt->gap_end=5;
        aln_rt->matrix=aln_sm_maq;
        aln_rt->row=5;
        aln_rt->band_width=50;
	return aln_rt;
}

static unsigned int approx_mapQ(const context_t *ctxt)
{
	unsigned int mapQ,n;
	unsigned c0=ctxt->c0;
	unsigned c1=ctxt->pos_n-ctxt->c0;

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

static unsigned int cal_seq_type(context_t *ctxt)
{
	int  j, best,second,score;
	second = 0 ;
	if(ctxt->pos_n ==0){
	ctxt->c1=ctxt->c0=0;
 	return 1;	
	}
	best=ctxt->pos[0].n_mm*3+ctxt->pos[0].n_gapo*11+ctxt->pos[0].n_gape*4;
	for(j=1;j<ctxt->pos_n;j++){
		score=ctxt->pos[j].n_mm*3+ctxt->pos[j].n_gapo*11+ctxt->pos[j].n_gape*4;
		if(score>best){
			second=score;
			break;
			}
		}
	ctxt->c0=j;
	for(;j<ctxt->pos_n;j++){
		score=ctxt->pos[j].n_mm*3+ctxt->pos[j].n_gapo*11+ctxt->pos[j].n_gape*4;
		if(score>second)
		break;
	}		
	ctxt->c1=j;
	return 0;	
}

static int  infer_isize(context_t *ctxt[2],isize_info_t *ii,double ap_prior,uint64_t L,int n_seq)
{
	uint64_t x,*isizes,n_ap= 0;
	int n, i,tot, p25, p75 ,p50 ,max_len =1 ,tmp ;
	double skewness =0.0,kurtosis =0.0,y;
	ii->avg = ii->std = -1.0;
	ii->low = ii->high = ii->high_bayesian = 0;
	isizes = (uint64_t*)calloc(n_seq, 8);
	for(i=0,tot=0;i!=n_seq;i++){
		context_t *p[2];
		p[0]=ctxt[0]+i;
		p[1]=ctxt[1]+i;
		if(p[0]->mapQ >= 20 && p[1]->mapQ >= 20){
		
			if(p[0]->pos_n==0 ||p[1]->pos_n==0)continue;
			x =( p[0]->pe_pos.pos >p[1]->pe_pos.pos )?(p[0]->pe_pos.pos+p[0]->seq->length -p[1]->pe_pos.pos):(p[1]->pe_pos.pos+p[1] ->seq ->length -p[0]->pe_pos.pos);
			if (x < 100000) isizes[tot++] = x;
			if (p[0]->seq->length > max_len) max_len = p[0]->seq->length;
			if (p[1]->seq->length > max_len) max_len = p[1]->seq->length;
		}
	}
	fprintf(stderr,"n_SEQ=%d\ttot=%d\n",n_seq,tot);
	if (tot < 20) {
	fprintf(stderr, "[infer_isize] fail to infer insert size: too few good pairs\n");
	fprintf(stderr,"%d\n",tot);
	FREE(isizes);
	return -1;
	}
	ks_introsort(uint64_t, tot, isizes);
	p25 = isizes[(int)(tot*0.25 + 0.5)];
	p50 = isizes[(int)(tot*0.50 + 0.5)];
	p75 = isizes[(int)(tot*0.75 + 0.5)];
	tmp  = (int)(p25 - OUTLIER_BOUND * (p75 - p25) + .499);
	ii->low = tmp > max_len? tmp : max_len; // ii->low is unsigned
	ii->high = (int)(p75 + OUTLIER_BOUND * (p75 - p25) + .499);
	for (i = 0, x = n = 0; i < tot; ++i)
		if (isizes[i] >= ii->low && isizes[i] <= ii->high)
			++n, x += isizes[i];
	ii->avg = (double)x / n;
	for (i = 0; i < tot; ++i) {
		if (isizes[i] >= ii->low && isizes[i] <= ii->high) {
			double tmp = (isizes[i] - ii->avg) * (isizes[i] - ii->avg);
			ii->std += tmp;
			skewness += tmp * (isizes[i] - ii->avg);
			kurtosis += tmp * tmp;
		}
	}
	kurtosis = kurtosis/n / (ii->std / n * ii->std / n) - 3;
	ii->std = sqrt(ii->std / n); // it would be better as n-1, but n is usually very large
	skewness = skewness / n / (ii->std * ii->std * ii->std);
	for (y = 1.0; y < 10.0; y += 0.01)
		if (.5 * erfc(y / M_SQRT2) < ap_prior / L * (y * ii->std + ii->avg)) break;
	ii->high_bayesian = (uint32_t)(y * ii->std + ii->avg + .499);
	for (i = 0; i < tot; ++i)
		if (isizes[i] > ii->high_bayesian) ++n_ap;
	ii->ap_prior = .01 * (n_ap + .01) / tot;
	if (ii->ap_prior < ap_prior) ii->ap_prior = ap_prior;
	FREE(isizes);
	fprintf(stderr, "[infer_isize] (25, 50, 75) percentile: (%d, %d, %d)\n", p25, p50, p75);
	if (isnan(ii->std) || p75 > 100000) {
		ii->low = ii->high = ii->high_bayesian = 0; ii->avg = ii->std = -1.0;
		fprintf(stderr, "[infer_isize] fail to infer insert size: weird pairing\n");
		return -1;
	}
	for (y = 1.0; y < 10.0; y += 0.01)
		if (.5 * erfc(y / M_SQRT2) < ap_prior / L * (y * ii->std + ii->avg)) break;
	ii->high_bayesian = (uint32_t)(y * ii->std + ii->avg + .499);
	fprintf(stderr, "[infer_isize] low and high boundaries: %d and %d for estimating avg and std\n", ii->low, ii->high);
	fprintf(stderr, "[infer_isize] inferred external isize from %d pairs: %.3lf +/- %.3lf\n", n, ii->avg, ii->std);
	fprintf(stderr, "[infer_isize] skewness: %.3lf; kurtosis: %.3lf; ap_prior: %.2e\n", skewness, kurtosis, ii->ap_prior);
	fprintf(stderr, "[infer_isize] inferred maximum insert size: %d (%.2lf sigma)\n", ii->high_bayesian, y);
	return 0;
}

#define SELECT_ID(x) x&0x00000000ffffffff	

static int pairing(context_t *ctxt_0,context_t *ctxt_1, pe_data_t *pe_data, const pe_opt_t *opt, int s_mm, const isize_info_t *ii)
{
	int i, j, o_n, subo_n,  low_bound = ii->low, max_len;
	uint64_t last_pos[2][2], o_pos[2], subo_score, o_score;
	uint32_t stack_id;
	ctxt_0->flag=0;
	ctxt_1->flag=0;
	max_len=ctxt_0->seq->length;
	if(max_len<ctxt_1->seq->length) max_len=ctxt_1->seq->length;
	if (low_bound < max_len) low_bound = max_len;
	o_score = subo_score = (uint64_t)-1;
	o_n = subo_n = 0;
	ks_introsort(uint64_t,pe_data->stack_n,pe_data->pos);
	for(j=0;j<2;j++) last_pos[j][0] = last_pos[j][1] = (uint64_t)-1;
	for(i=0;i<pe_data->stack_n;i++){
		uint64_t x=pe_data->pos[i];
		stack_id=SELECT_ID(x);
		int strand=pe_data->pos_pe[stack_id].a;
		if(strand ==1){
		int y=1-pe_data->pos_pe[stack_id].reads_id;
		//__pairing_aux(last_pos[y][1],x);
			{
			int32_t l= (x>>32)-(last_pos[y][1]>>32)+((y==0)?ctxt_1->seq->length:ctxt_0->seq->length);
			uint32_t u_id=SELECT_ID(last_pos[y][1]);
			if(last_pos[y][1]!=(uint64_t)-1 && x>>32 >last_pos[y][1]>>32 && l>=max_len && ((ii->high && l <= ii->high_bayesian) || (ii->high == 0 && l <= opt->max_isize)) ){
				uint64_t s = pe_data->pos_pe[stack_id].score + pe_data->pos_pe[u_id].score;
				s *= 10;
				if (ii->high) s += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * fabs(l - ii->avg) / ii->std)) + .499);
				s = s<<32 | (uint32_t)hash_64((last_pos[y][1])>>32<<32 | (x)>>32);
				if (s>>32 == o_score>>32) ++o_n; 
				else if (s>>32 < o_score>>32) { subo_n += o_n; o_n = 1; }
				else ++subo_n;
				if (s < o_score){ subo_score = o_score, o_score = s, o_pos[pe_data->pos_pe[u_id].reads_id] = last_pos[y][1], o_pos[pe_data->pos_pe[stack_id].reads_id] = x; }
                        	else if (s < subo_score) subo_score = s;                                       
               				 }                                              
			}
				
			{	
                        int32_t l= (x>>32)-(last_pos[y][0]>>32)+((y==0)?ctxt_1->seq->length:ctxt_0->seq->length);
                        uint32_t u_id = SELECT_ID(last_pos[y][0]);
                        if(last_pos[y][0]!=(uint64_t)-1 && x>>32 >last_pos[y][0]>>32 && l>=max_len && ((ii->high && l <= ii->high_bayesian) || (ii->high == 0 && l <= opt->max_isize)) ){
                        uint64_t s = pe_data->pos_pe[stack_id].score + pe_data->pos_pe[u_id].score;
                         s *= 10;
                        if (ii->high) s += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * fabs(l - ii->avg) / ii->std)) + .499);
                        s = s<<32 | (uint32_t)hash_64( last_pos[y][0]>>32<<32 | (x)>>32);
                        if (s>>32 == o_score>>32) ++o_n;
                        else if (s>>32 < o_score>>32) { subo_n += o_n; o_n = 1; }
                        else ++subo_n;
                        if (s < o_score){ subo_score = o_score, o_score = s, o_pos[pe_data->pos_pe[u_id].reads_id] = last_pos[y][0], o_pos[pe_data->pos_pe[stack_id].reads_id] = x; }
                        else if (s < subo_score) subo_score = s;
                         	}
        		}       

	}else{
		int	y1 =pe_data->pos_pe[stack_id].reads_id;
		last_pos[y1][0]=last_pos[y1][1];
		last_pos[y1][1]=x;
		}	
	}
	if (o_score != (uint64_t)-1){
	int mapQ_p = 0; // this is the maximum mapping quality when one end is moved
	int rr[2];
	if (o_n == 1) {
			if (subo_score == (uint64_t)-1) mapQ_p = 29; // no sub-optimal pair
			else if ((subo_score>>32) - (o_score>>32) > s_mm * 10) mapQ_p = 23; // poor sub-optimal pair
			else {
				int n = subo_n > 255? 255 : subo_n;
				mapQ_p = ((subo_score>>32) - (o_score>>32)) / 2 - g_log_n[n];
				if (mapQ_p < 0) mapQ_p = 0;
			}
	}
	
	uint32_t o_pos_stack[2];
	o_pos_stack[0]=SELECT_ID(o_pos[0]);
	o_pos_stack[1]=SELECT_ID(o_pos[1]);
	rr[0] =pe_data->pos_pe[o_pos_stack[0]].a;
	rr[1] =pe_data->pos_pe[o_pos_stack[1]].a;
	if((ctxt_0->pe_pos.pos ==o_pos[0]>>32 && ctxt_0->pe_pos.a == rr[0]) &&(ctxt_1->pe_pos.pos ==o_pos[1]>>32 && ctxt_1->pe_pos.a == rr[1]))	{
		if(ctxt_0->mapQ >0 &&ctxt_1->mapQ >0) {
			int mapQ = ctxt_0->mapQ + ctxt_1->mapQ;
			if (mapQ > 60) mapQ = 60;
			ctxt_0->mapQ = ctxt_1->mapQ = mapQ;
		} else {
			if (ctxt_0->mapQ == 0) ctxt_0->mapQ = (mapQ_p + 7 < ctxt_1->mapQ)? mapQ_p + 7 : ctxt_1->mapQ;
			if (ctxt_1->mapQ == 0) ctxt_1->mapQ = (mapQ_p + 7 < ctxt_0->mapQ)? mapQ_p + 7 : ctxt_0->mapQ;
		}
	}

	else if(ctxt_0->pe_pos.pos ==o_pos[0]>>32 && ctxt_0->pe_pos.a == rr[0]){
		ctxt_1->seQ =  0 ;
		ctxt_1->mapQ = ctxt_0->mapQ;
		if(ctxt_1->mapQ>mapQ_p) ctxt_1->mapQ = mapQ_p;
	}
	else if((ctxt_1->pe_pos.pos ==o_pos[1]>>32 && ctxt_1->pe_pos.a == rr[1])){
		ctxt_0->seQ =  0 ;
		ctxt_0->mapQ = ctxt_1 ->mapQ;
		if(ctxt_0->mapQ >mapQ_p) ctxt_0->mapQ = mapQ_p;
	}
	else {
		ctxt_0->seQ = ctxt_1->seQ = 0;
		mapQ_p -= 20;
		if (mapQ_p < 0) mapQ_p = 0;
		ctxt_1->mapQ = ctxt_0->mapQ = mapQ_p;
	}
	cnt_chg++;
	ctxt_0->flag =2;
	ctxt_1->flag =2;
	if(ctxt_0->pe_pos.pos !=o_pos[0]>>32|| ctxt_0->pe_pos.a!=rr[0]){
		ctxt_0->pe_pos.pos=o_pos[0]>>32;
		ctxt_0->pe_pos.n_mm =pe_data->pos_pe[o_pos_stack[0]].n_mm;
		ctxt_0->pe_pos.a = pe_data->pos_pe[o_pos_stack[0]].a;
		ctxt_0->pe_pos.n_gapo =pe_data->pos_pe[o_pos_stack[0]].n_gapo;
		ctxt_0->pe_pos.n_gape =pe_data->pos_pe[o_pos_stack[0]].n_gape;
		}
	if(ctxt_1->pe_pos.pos !=o_pos[1]>>32 || ctxt_1->pe_pos.a!=rr[1]){
                ctxt_1->pe_pos.pos = o_pos[1]>>32;
		ctxt_1->pe_pos.a = pe_data->pos_pe[o_pos_stack[1]].a;
                ctxt_1->pe_pos.n_mm = pe_data->pos_pe[o_pos_stack[1]].n_mm;
                ctxt_1->pe_pos.n_gapo =pe_data->pos_pe[o_pos_stack[1]].n_gapo;
                ctxt_1->pe_pos.n_gape =pe_data->pos_pe[o_pos_stack[1]].n_gape;
       		}
	}

	return cnt_chg;
}	

# define SELECT_CHAR(x) ( x == 0 ? 'A' : ( x ==1 ? 'C':( x== 2 ? 'G' : ( x==3 ? 'T': 'N' ) ) ) )

int	cm_cal_md(int *n_cigar,bwa_cigar_t *cigar,char *ref_nt4, char *nt4, int ref_len,char *md, uint8_t *nm)
{
		int i ,j, l, d;
	int x ,y;
	d = x = y = 0; 
	uint8_t _nm =0;
	if(cigar)
	{
		for(i = 0; i < *n_cigar ; i++){
			l = __cigar_len(cigar[i]);
			if(__cigar_op(cigar[i])==FROM_M){
				for(j = 0;j < l && x + j < ref_len; j++){
					if(ref_nt4[x+j] >3 || nt4[y+j] > 3 || ref_nt4[x+j] != nt4[y+j]){
						md += sprintf(md, "%d%c", d, SELECT_CHAR(ref_nt4[y+j])) ;
						_nm ++;
						d =0;
						} 
					else	d++;
					}
					x += l; y += l;
				}
			else	if(__cigar_op(cigar[i]) == FROM_I||__cigar_op(cigar[i])== FROM_S){
					y +=l;
					if(__cigar_op(cigar[i]) == FROM_I) _nm += l;
				}
			else	if( __cigar_op(cigar[i]) ==FROM_D){
				md += sprintf(md, "%d^", d);
				d=0;
					for(j = 0;j < l && x + j < ref_len; j++)
						md += sprintf(md, "%c",SELECT_CHAR(ref_nt4[y+j]));
				x += l;
				}
		}
		if(d > 0)
		sprintf(md, "%d",d);
		
		*nm=_nm;
		return 0;
		
	}
	return 1;
}


void bwa_fill_scmat(int a, int b, int8_t mat[25])
{
	int i, j, k;
	for (i = k = 0; i < 4; ++i) {
		for (j = 0; j < 4; ++j)
			mat[k++] = i == j? a : -b;
		mat[k++] = -1; // ambiguous base
	}
	for (j = 0; j < 5; ++j) mat[k++] = -1;
}


/*
 *   Replace the local alignment function with optimization function based on SSE2.
 *
*/

bwa_cigar_t *cm_sw_core(refdb_t refdb,context_t *ctxt,char *nt4,int64_t *beg,unsigned int reglen,int *n_cigar, uint32_t *cnt, uint64_t  L,char *md)
{
	int i,x,y;
	int rc;
	unsigned int ref_pos=*beg;
	bwa_cigar_t *cigar = 0;

	kswr_t  r ;
	int8_t mat[25] ;
	bwa_fill_scmat(1,3,mat);
	int xtra , gscore ;
	uint32_t  *cigar32 = 0 ;
	
	if(reglen <20) return 0;
	for(i=0,x=0;i<ctxt->seq->length;i++)
		if(nt4[i]>=4) x++;
	if((float)x/ctxt->seq->length >=0.25 || ctxt->seq->length-x<20) return 0;
		
	if(*beg+reglen>L) ctxt->ref_len=L-*beg+1;
	else	ctxt->ref_len =reglen;
	
	
	realloc_n(ctxt->ref_nt4,ctxt->ref_nt4_sz,ctxt->ref_len);
	rc =refdb_seq_nt4(refdb,ctxt->ref_nt4,ref_pos,ctxt->ref_len);
	if(rc){
		fprintf(stderr,"Can not retrieve reference with POS=%u\n",ref_pos);
		exit(-1);
	}
	xtra = KSW_XSUBO | KSW_XSTART | (  ctxt->seq->length < 250 ? KSW_XBYTE : 0 ) ;
	r = ksw_align(ctxt->seq->length ,(uint8_t *)nt4 , ctxt->ref_len ,(uint8_t *)ctxt->ref_nt4 , 5 , mat , 5 , 1 , xtra , 0);
	gscore = ksw_global(r.qe-r.qb+1 , (uint8_t *)(nt4 + r.qb) , r.te - r.tb + 1 , (uint8_t *)(ctxt->ref_nt4 + r.tb) ,5 ,mat , 5 , 1 ,50 , n_cigar , &cigar32 );
	
	cigar = (bwa_cigar_t *)cigar32 ; 
	for( i = 0 ;  i < *n_cigar ; i++)
		cigar[i] = __cigar_create((cigar32[i]&0xf), (cigar32[i]>>4));

	
	if(r.score < 20   ||  r.score2 == r.score  || gscore != r.score){
		FREE(cigar);
		FREE(ctxt->stdaln_rt);
		FREE(ctxt->ref_nt4);
		ctxt->ref_nt4_sz = 0;
		return 0;
	}

	for(i=0,x=y=0;i<*n_cigar;i++){
		bwa_cigar_t c = cigar[i];
		if (__cigar_op(c) == FROM_M)x += __cigar_len(c), y += __cigar_len(c);
		else if (__cigar_op(c) == FROM_D) x += __cigar_len(c);
		else y += __cigar_len(c);
	}
	if(x <20 ||y<20){
		FREE(cigar);
		FREE(ctxt->stdaln_rt->path);
		FREE(ctxt->stdaln_rt);
		FREE(ctxt->ref_nt4);
		ctxt->ref_nt4_sz=0;
		return 0;
	}

	{ // update cigar and coordinate;
		int start =  r.qb , end  = r.qe + 1;
		*beg +=  r.tb ;
		cigar = (bwa_cigar_t*)realloc(cigar, sizeof(bwa_cigar_t) * (*n_cigar + 2));
		if (start) {
			memmove(cigar + 1, cigar, sizeof(bwa_cigar_t) * (*n_cigar));
			cigar[0] = __cigar_create(3,start);
			++(*n_cigar);
		}
		if (end <ctxt->seq->length ) {
			cigar[*n_cigar] = __cigar_create(3,(ctxt->seq->length- end));
			++(*n_cigar);
		}
	}
	{ // set *cnt
		int k;
		int n_mm, n_gapo, n_gape ,l;
		n_mm = n_gapo = n_gape = 0;
		x = r.tb ; y = r.qb ; 
		for (k = 0; k < *n_cigar; ++k) {
			bwa_cigar_t c = cigar[k];
			if (__cigar_op(c) == FROM_M) {
				for (l = 0; l < (__cigar_len(c)); ++l)
					if (ctxt->ref_nt4[x+l] < 4 && nt4[y+l] < 4 && ctxt->ref_nt4[x+l] != nt4[y+l]) ++n_mm;
				x += __cigar_len(c), y += __cigar_len(c);
			} else if (__cigar_op(c) == FROM_D) {
				x += __cigar_len(c), ++n_gapo, n_gape += (__cigar_len(c)) - 1;
			} else if (__cigar_op(c) == FROM_I) {
				y += __cigar_len(c), ++n_gapo, n_gape += (__cigar_len(c)) - 1;
			}
		}
		*cnt = (uint32_t)n_mm<<16 | n_gapo<<8 | n_gape;
	}
	ref_pos =*beg;
	// translate
	FREE(ctxt->ref_nt4);
	ctxt->ref_nt4_sz=0;
	realloc_n(ctxt->ref_nt4,ctxt->ref_nt4_sz,ctxt->ref_len);
        rc =refdb_seq_nt4(refdb,ctxt->ref_nt4,ref_pos,ctxt->ref_len);
        if(rc){
                fprintf(stderr,"Can not retrieve reference with POS=%u\n",ref_pos);
        	exit(-1);
	}

	
	cm_cal_md(n_cigar,cigar,ctxt->ref_nt4,nt4,ctxt->ref_len,md,&ctxt->pe_pos.n_mm);
	FREE(ctxt->stdaln_rt);
	FREE(ctxt->ref_nt4);
	ctxt->ref_nt4_sz=0;
	return cigar;
}

static int   *cm_paired_sw(context_t *ctxt_0,context_t *ctxt_1,isize_info_t *ii,refdb_t refdb,pe_opt_t *opt,uint64_t L)
{
	char md[2][4086];
	uint64_t n_tot[2], n_mapped[2];
	

	n_tot[0] = n_tot[1] = n_mapped[0] = n_mapped[1] = 0;
	if( ctxt_1->pe_pos.pos >= L || ctxt_0->pe_pos.pos >= L  ){

		fprintf(stderr, "%s\n",ctxt_1->seq->name);
	}
	if(((ctxt_0->mapQ >=SW_MIN_MAPQ ||ctxt_1->mapQ >= SW_MIN_MAPQ )||(ctxt_1->pos_n == 0 && ctxt_0->pos_n >0 )||(ctxt_0->pos_n == 0 && ctxt_1->pos_n >0)) && (ctxt_0->flag && 2)==0 ){
		int k, n_cigar[2],  mapQ = 0, mq_adjust[2];
		int64_t beg[2], end[2];
		uint32_t cnt[2];
		bwa_cigar_t *cigar[2];
		#define __set_rght_coor(_a, _b, _pref, _pmate) do {					\
			if(_pref->pe_pos.pos >L) \
				(_a) = (_b) = L;\
			else{\
				(_a)	= (int64_t)_pref->pe_pos.pos + ii->avg - 3 * ii->std -_pmate->seq->length*1.5 ;\
				(_b)	= (_a) + 6 * ii->std + 2*_pmate->seq->length;\
				if((_a) < (int64_t)_pref->pe_pos.pos + _pref->seq->length) (_a) = _pref->pe_pos.pos + _pref->seq->length;	\
				if((_b)>L) (_b)=L; \
				if((_a)>L) (_a)=L;\
			}\
			}while(0)
		#define __set_left_coor(_a,_b,_pref, _pmate) do {\
			if(_pref->pe_pos.pos > L){						\
				(_a) = (_b) =L ;\
			}\
			else{\
				(_a)	=(int64_t)_pref->pe_pos.pos + _pref->seq->length -ii->avg - 3 * ii->std	-_pmate->seq->length*0.5;	\
				(_b)=	(_a) + 6 * ii->std + 2*_pmate->seq->length;	\
				if((_a)<0) (_a)=0;\
				if((_b)>_pref->pe_pos.pos)   (_b)=_pref->pe_pos.pos;      \
				if((_b)<0) (_b)=0;\
				if((_a)>L) (_a)=L; \
				if((_b)>L) (_b)=L; \
			}\
			}while(0)

		cigar[0] = cigar[1] = 0;
		n_cigar[0] = n_cigar[1] = 0;
		mq_adjust[0] = mq_adjust[1] = 255;
		if(ctxt_1->pos_n==0)
			goto NO_MATCH1;
		if(ctxt_1->pe_pos.a==0){
	//		ctxt_0->do_sw = 1;
			__set_rght_coor(beg[0], end[0], ctxt_1, ctxt_0);
			realloc_n(ctxt_0->seq_comp,ctxt_0->seq_comp_sz, ctxt_0->seq->length+1);
			realloc_n(ctxt_0->nt4[0],ctxt_0->nt4_sz[0],ctxt_0->seq->length);
			seq_to_comp(ctxt_0->seq_comp,ctxt_0->seq->seq,ctxt_0->seq->length);
			ctxt_0->seq_comp[ctxt_0->seq->length]=0;
			seq_to_nt4(ctxt_0->nt4[0],ctxt_0->seq_comp,ctxt_0->seq->length);
			cigar[0]=cm_sw_core(refdb,ctxt_0,ctxt_0->nt4[0],&beg[0],end[0] - beg[0],&n_cigar[0],&cnt[0],L,md[0]);
			ctxt_0->seq_comp_sz = 0;
			ctxt_0->nt4_sz[0] =0;
			FREE(ctxt_0->seq_comp);
			FREE(ctxt_0->nt4[0]);
			}	
		else{
	//		 ctxt_0->do_sw = 1;//
			__set_left_coor( beg[0], end[0],ctxt_1,ctxt_0);
			realloc_n(ctxt_0->nt4[0],ctxt_0->nt4_sz[0],ctxt_0->seq->length);
			seq_to_nt4(ctxt_0->nt4[0],ctxt_0->seq->seq,ctxt_0->seq->length);
			cigar[0]=cm_sw_core(refdb,ctxt_0,ctxt_0->nt4[0],&beg[0],end[0] -beg[0],&n_cigar[0], &cnt[0],L,md[0]);
			ctxt_0->nt4_sz[0] =0;
			FREE(ctxt_0->nt4[0]);

		}
		if(cigar[0]&&ctxt_0->pos_n!=0){
			int s_old, clip = 0, s_new;
			if (__cigar_op(cigar[0][0]) == 3) clip += __cigar_len(cigar[0][0]);
			if (__cigar_op(cigar[0][n_cigar[0]-1]) == 3) clip += __cigar_len(cigar[0][n_cigar[0]-1]);
			s_old = (int)((ctxt_0->pe_pos.n_mm * 9 + ctxt_0->pe_pos.n_gapo * 13 + ctxt_0->pe_pos.n_gape * 2) / 3. * 8. + .499);
			s_new = (int)(((cnt[0]>>16) * 9 + (cnt[0]>>8&0xff) * 13 + (cnt[0]&0xff) * 2 + clip * 3) / 3. * 8. + .499);
			s_old += -4.343 * log(ii->ap_prior / L);
			s_new += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * 1.5) + .499));
			if (s_old < s_new) { // reject SW alignment
				mq_adjust[0] = s_new - s_old;
				FREE(cigar[0]);  n_cigar[0] = 0;
			} else mq_adjust[0] = s_old - s_new;
		}	
	
		
		NO_MATCH1:
		if(ctxt_0->pos_n==0)
                        goto NO_MATCH2;
		if(ctxt_0->pe_pos.a==0){
//			ctxt_1->do_sw = 1;//
			__set_rght_coor(beg[1], end[1], ctxt_0, ctxt_1);
			realloc_n(ctxt_1->seq_comp,ctxt_1->seq_comp_sz, ctxt_1->seq->length+1);
                        realloc_n(ctxt_1->nt4[0],ctxt_1->nt4_sz[0],ctxt_1->seq->length);
                        seq_to_comp(ctxt_1->seq_comp,ctxt_1->seq->seq,ctxt_1->seq->length);
                        ctxt_1->seq_comp[ctxt_1->seq->length]=0;
                        seq_to_nt4(ctxt_1->nt4[0],ctxt_1->seq_comp,ctxt_1->seq->length);
			cigar[1]=cm_sw_core(refdb,ctxt_1,ctxt_1->nt4[0],&beg[1],end[1]-beg[1] ,&n_cigar[1],&cnt[1],L,md[1]);
			ctxt_1->nt4_sz[0] =0;
			ctxt_1->seq_comp_sz =0;
			FREE(ctxt_1->seq_comp);
			FREE(ctxt_1->nt4[0]);

		}

		else{
//			ctxt_1->do_sw = 1;//
			__set_left_coor(beg[1], end[1],ctxt_0,ctxt_1);
			realloc_n(ctxt_1->nt4[0],ctxt_1->nt4_sz[0],ctxt_1->seq->length);
                        seq_to_nt4(ctxt_1->nt4[0],ctxt_1->seq->seq,ctxt_1->seq->length);
			cigar[1]=cm_sw_core(refdb,ctxt_1,ctxt_1->nt4[0],&beg[1],end[1]-beg[1],&n_cigar[1], &cnt[1],L,md[1]);
			ctxt_1->nt4_sz[0] =0;
			FREE(ctxt_1->nt4[0]);
			
		}
		if(cigar[1]&&ctxt_1->pos_n!=0){
                        int s_old, clip = 0, s_new;
                        if (__cigar_op(cigar[1][0]) == 3) clip += __cigar_len(cigar[1][0]);
                        if (__cigar_op(cigar[1][n_cigar[1]-1]) == 3) clip += __cigar_len(cigar[1][n_cigar[1]-1]);
			s_old = (int)((ctxt_1->pe_pos.n_mm * 9 + ctxt_1->pe_pos.n_gapo * 13 + ctxt_1->pe_pos.n_gape * 2) / 3. * 8. + .499);
                        s_new = (int)(((cnt[1]>>16) * 9 + (cnt[1]>>8&0xff) * 13 + (cnt[1]&0xff) * 2 + clip * 3) / 3. * 8. + .499);
				s_old += -4.343 * log(ii->ap_prior / L);
                        s_new += (int)(-4.343 * log(.5 * erfc(M_SQRT1_2 * 1.5) + .499));
                        if (s_old < s_new) { // reject SW alignment
                                mq_adjust[1] = s_new - s_old;
                                FREE(cigar[1]);  n_cigar[1] = 0;
                                } else mq_adjust[1] = s_old - s_new;
                        }
		NO_MATCH2:
		k=-1;

//		if(ctxt_0->do_sw && ctxt_1->do_sw)
//			ctxt_0->do_sw = ctxt_1->do_sw =2;

		if(cigar[0]&&cigar[1]){
		k=(ctxt_0->mapQ >ctxt_1->mapQ)?1:0;
		mapQ = abs(ctxt_1->mapQ - ctxt_0->mapQ);	
		}
		else if(cigar[0]) k = 0, mapQ = ctxt_1->mapQ;
		else if(cigar[1]) k = 1, mapQ = ctxt_0->mapQ;
		context_t *ctxt_k,*ctxt_1_k;
		if(k==0){
		ctxt_k =ctxt_0;
		ctxt_1_k =ctxt_1;
		}
		else{
		ctxt_k =ctxt_1;
                ctxt_1_k =ctxt_0;
		}
		if(k >=0 && ctxt_k->pe_pos.pos!=beg[k]){
			int tmp = (int)ctxt_1_k->mapQ - ctxt_k->mapQ/2 - 8;
			if (tmp <= 0) tmp = 1;
			if (mapQ > tmp) mapQ = tmp;
			ctxt_k->mapQ = ctxt_1_k->mapQ = mapQ;
			ctxt_k->seQ = ctxt_1_k->seQ = ctxt_1_k->seQ < mapQ?ctxt_1_k->seQ:mapQ;
			if (ctxt_k->mapQ > mq_adjust[k]) ctxt_k->mapQ = mq_adjust[k];
		
			ctxt_k->is_sw=1;
			ctxt_k->pe_pos.pos=beg[k];
	//		ctxt_k->pe_pos.n_mm = cnt[k]>>16;
			ctxt_k->pe_pos.n_gapo = cnt[k]>>8&0xff;
			ctxt_k->pe_pos.n_gape = cnt[k]& 0xff;
			ctxt_k->pe_pos.a = 1 - ctxt_1_k->pe_pos.a;
			ctxt_k->flag |=2;
			ctxt_k->flag |=2;
			ctxt_1_k->flag |=2;
			ctxt_k ->seQ = ctxt_1_k->seQ ;
			ctxt_k->n_cigar =n_cigar[k];
			ctxt_k->cigar=calloc((n_cigar[k]+1),sizeof(int32_t));
			// caculate!!

			ctxt_k->md =calloc(strlen(md[k])+1,sizeof(char));
			ctxt_k->n_md=strlen(md[k]);
			memcpy(ctxt_k->md,md[k],(ctxt_k->n_md+1)*sizeof(char));
			ctxt_k->md[ctxt_k->n_md] ='\0';
			memcpy(ctxt_k->cigar,cigar[k],n_cigar[k]*sizeof(uint32_t));
			ctxt_k->cigar[n_cigar[k]] ='\0';
			FREE(cigar[0]);
			FREE(cigar[1]);	
		}

		FREE(ctxt_0->stdaln_rt);
		FREE(ctxt_1->stdaln_rt);
		return 0;
	}
	FREE(ctxt_0->stdaln_rt);
        FREE(ctxt_1->stdaln_rt);
	return  0;
}

int  sw2cigar(char *cigar_,bwa_cigar_t *cigar,int n_cigar)
{
	int i;
	if(cigar){
		for(i=0;i<n_cigar;i++){
			cigar_+=sprintf(cigar_,"%d%c", __cigar_len(cigar[i]), "MIDS"[__cigar_op(cigar[i])]);
		}
	}else cigar_="*";
	return 0;
}

int getPosValue(context_t *ctxt)
{
	int j = 0;
	int64_t x = 0;
	for (j = 0; j != ctxt->n_cigar; ++j) {
		int op = __cigar_op(ctxt->cigar[j]);
		if (op == 3){
			 x += __cigar_len(ctxt->cigar[j]);
			 break;
		}
	}
	return x;
}
uint32_t  getLength(char *cigar){
        uint32_t length = 0;
        int i = 0;
        char length_arr[BUF];
        char *m_ptr, *d_ptr, *pm, *pd;
        m_ptr = strstr(cigar,"M");
        d_ptr = strstr(cigar,"D");
       	pm = cigar;
        pd = cigar;
	memset(length_arr,0,BUF*sizeof(char));
        while(m_ptr||d_ptr){
                for(i = 0; pm < m_ptr; pm++, i++ ){
                        if(pm[0]=='D'||pm[0]=='I'||pm[0]=='S'){
                                memset(length_arr,0,BUF*sizeof(char));
                                i = -1;
                                }
                        else    length_arr[i] = *pm;
                }
                length += atoi(length_arr);
                if(m_ptr!=NULL){
                        pm = m_ptr;
                        m_ptr++;
                        m_ptr = strstr(m_ptr,"M");
                }
                memset(length_arr,0,BUF*sizeof(char));
                for(i = 0; pd < d_ptr; pd++, i++){
                        if(pd[0]=='M'||pd[0]=='I'||pd[0]=='S'){
                                memset(length_arr,0,BUF*sizeof(char));
                                i = -1;
                        }
                        else    length_arr[i] = *pd;
                }
                length += atoi(length_arr);
                if(d_ptr!=NULL){
                        pd = d_ptr;
                        d_ptr ++;
                        d_ptr = strstr(d_ptr,"D");
                }
                memset(length_arr,0,BUF*sizeof(char));
       	}
	return length;
}

int fout_sam_info(sam_info_t *sam_info0,sam_info_t *sam_info1,context_t *ctxt_0, context_t *ctxt_1,int slt ,FILE *fout,refdb_t refdb)
{

			int flag[2];
			char	cigar[4086];
			int	ref_id = 0;
			unsigned int am;
			seqdb_t seqdb= NULL ;
			int	i;
			int	j=0;
			switch(slt){
			case(0)://all ummap!!
					{	
					flag[0] =0x40|0x4|0x8|0x1;
					flag[1] =0x80|0x4|0x8|0x1;	
					//reads1 info
					fprintf(fout,"%s\t%u\t*\t0\t0\t*\t*\t0\t0\t",ctxt_0->seq->name,flag[0]);
					fprintf(fout,"%s\t%s",ctxt_0->seq->seq,ctxt_0->seq->qual);
					fprintf(fout,"\n");
					//reads2 info
					fprintf(fout,"%s\t%u\t*\t0\t0\t*\t*\t0\t0\t",ctxt_1->seq->name,flag[1]);
					fprintf(fout,"%s\t%s",ctxt_1->seq->seq,ctxt_1->seq->qual);
					fprintf(fout,"\n");
					//reads2 info
				//free
					seqdb_release(seqdb,ctxt_0->seq);
					seqdb_release(seqdb,ctxt_1->seq);
					break;
				}
			case(1)://reads 0 ummap!!
				{
					flag[0] =0x40|0x4|0x1;
					if(ctxt_1->pe_pos.a){
						flag[1] =ctxt_1->flag|0x80|0x8|0x1|0x10;
						flag[0] |= 0x10|0x20;
						}
					else	flag[1] =ctxt_1->flag|0x80|0x8|0x1;
					//reads1 info
					fprintf(fout,"%s\t%u\t",ctxt_0->seq->name,flag[0]);
					fprintf(fout,"%s\t%u\t",sam_info1->ref_name,ctxt_1->pe_pos.pos);
					fprintf(fout,"0\t*\t=\t%u\t0\t",ctxt_1->pe_pos.pos);
					fprintf(fout,"%s\t%s",ctxt_0->seq->seq,ctxt_0->seq->qual);
					fprintf(fout,"\n");
					
					//reads2 info
					fprintf(fout,"%s\t%u\t", ctxt_1->seq->name,flag[1]);
					fprintf(fout,"%s\t%u\t", sam_info1->ref_name,ctxt_1->pe_pos.pos);
					fprintf(fout,"%u\t%s\t", ctxt_1->mapQ,sam_info1->cigar);
					fprintf(fout,"=\t%u\t0\t",ctxt_1->pe_pos.pos);
					// reads2  seq and qual
					if(ctxt_1->pos->a)
						fprintf(fout,"%s\t%s\t",ctxt_1->seq_comp,ctxt_1->qual_comp);
					else	fprintf(fout,"%s\t%s\t",ctxt_1->seq->seq,ctxt_1->seq->qual);

					
					if(ctxt_1->c0 > 1)
						fprintf(fout,"XT:A:R\t");
					else
						fprintf(fout,"XT:A:U\t");
					fprintf(fout,"NM:i:%u\t",ctxt_1->pe_pos.n_mm);
					if(sam_info1->ref_amb)
					fprintf(fout,"XN:i:%u\t",sam_info1->ref_amb);
				
					fprintf(fout,"SM:i:%u\tAM:i:0\t",ctxt_1->seQ);	
					fprintf(fout,"X0:i:%u\tX1:i:%u\t",ctxt_1->c0,ctxt_1->c1-ctxt_1->c0);
					fprintf(fout,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_1->pe_pos.n_mm,ctxt_1->pe_pos.n_gapo,ctxt_1->pe_pos.n_gapo+ctxt_1->pe_pos.n_gape);

					fprintf(fout,"MD:Z:%s",sam_info1->md);
					for(i=0;i<sam_info1->num;i++){
						pos_t  *pos = ctxt_1->pos+i*sam_info1->step;
						if(sam_info1->repeat){	
							if(ctxt_1->pe_pos.pos!=pos->pos){
								if( j == 0 )
								fprintf(fout,"\tXA:Z:");
								if( j < 3 )
								fprintf(fout,"%s,%c%u,%s,%u;",sam_info1->m_ref_name[i],pos->a?'-':'+',pos->pos,sam_info1->m_cigar[i],pos->n_mm);
								j++;
							}
						}
						FREE(sam_info1->m_cigar[i]);
						FREE(sam_info1->m_ref_name[i]);
					}
					fprintf(fout,"\n");
					seqdb_release(seqdb,ctxt_0->seq);
					seqdb_release(seqdb,ctxt_1->seq);
					FREE(ctxt_1->pos);
					FREE(ctxt_1->qual_comp);
					FREE(ctxt_1->seq_comp);
					FREE(sam_info1->cigar);
					FREE(sam_info1->md);
					FREE(sam_info1->ref_name);
					break;
					}
			case(2):{
					flag[1] =0x80|0x4|0x1;
					if(ctxt_0->pe_pos.a){
						flag[0]=0x01|0x08|0x10|0x40;
						flag[1] |=0x10|0x20;
						}
					else
						flag[0]=0x01|0x8|0x40;
					//reads 1 info;
					fprintf(fout,"%s\t%u\t",ctxt_0->seq->name,flag[0]);
                                        fprintf(fout,"%s\t%u\t",sam_info0->ref_name,ctxt_0->pe_pos.pos);
                                        fprintf(fout,"%u\t%s\t",ctxt_0->mapQ,sam_info0->cigar);
                                        fprintf(fout,"=\t%u\t0\t",ctxt_0->pe_pos.pos);

					if(ctxt_0->pos->a)
                                                fprintf(fout,"%s\t%s\t",ctxt_0->seq_comp,ctxt_0->qual_comp);
                                        else    fprintf(fout,"%s\t%s\t",ctxt_0->seq->seq,ctxt_0->seq->qual);


                                        if(ctxt_0->c0 >1)
                                                fprintf(fout,"XT:A:R\t");
                                        else
                                                fprintf(fout,"XT:A:U\t");
                                        fprintf(fout,"NM:i:%u\t",ctxt_0->pe_pos.n_mm);
                                        if(sam_info0->ref_amb)
                                        fprintf(fout,"XN:i:%u\t",sam_info0->ref_amb);
					fprintf(fout,"SM:i:%u\tAM:i:0\t",ctxt_0->seQ);	
					
					fprintf(fout,"X0:i:%u\tX1:i:%u\t",ctxt_0->c0,ctxt_0->c1-ctxt_0->c0);
                                        fprintf(fout,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_0->pe_pos.n_mm,ctxt_0->pe_pos.n_gapo,ctxt_0->pe_pos.n_gapo+ctxt_0->pe_pos.n_gape);

                                        fprintf(fout,"MD:Z:%s",sam_info0->md);
					for(i=0;i<sam_info0->num;i++){
                                                pos_t  *pos = ctxt_0->pos+i*sam_info0->step;
						if(sam_info0->repeat){
							if(pos->pos!=ctxt_0->pe_pos.pos){
								if( j == 0)
								fprintf(fout,"\tXA:Z:");
								if( j < 3)
                                                        	fprintf(fout,"%s,%c%u,%s,%u;",sam_info0->m_ref_name[i],pos->a?'-':'+',pos->pos,sam_info0->m_cigar[i],pos->n_mm);
								j++;
							}
											
						}
						FREE(sam_info0->m_cigar[i]);
						FREE(sam_info0->m_ref_name[i]);
					}
                                        fprintf(fout,"\n");
					// reads 2 info;
					fprintf(fout,"%s\t%u\t",ctxt_1->seq->name,flag[1]);
					fprintf(fout,"%s\t%u\t0\t*\t",sam_info0->ref_name,ctxt_0->pe_pos.pos);
					fprintf(fout,"=\t%d\t0\t",ctxt_0->pe_pos.pos);
                                        fprintf(fout,"%s\t%s",ctxt_1->seq->seq,ctxt_1->seq->qual);
					fprintf(fout,"\n");
					seqdb_release(seqdb,ctxt_0->seq);
					seqdb_release(seqdb,ctxt_1->seq);

                                        FREE(ctxt_0->pos);
                                        FREE(ctxt_0->qual_comp);
                                        FREE(ctxt_0->seq_comp);
                                        FREE(sam_info0->cigar);
                                        FREE(sam_info0->md);
                                        FREE(sam_info0->ref_name);

					break;
					
					

				}
			case(3):
				{
				uint32_t cor_pos1,cor_pos2;
				unsigned int nn1,nn2;
				cor_pos1 = cor_pos2 = 0;
				am = (ctxt_1->seQ>ctxt_0->seQ?ctxt_0->seQ:ctxt_1->seQ);
				if(ctxt_0->pe_pos.a)
				cor_pos1= getLength(sam_info0->cigar);
				
				if(ctxt_1->pe_pos.a)
                                cor_pos2= getLength(sam_info1->cigar);
						
				if(ctxt_0->pe_pos.a){
					flag[0]=ctxt_0->flag|0x1|0x10|0x40;
					flag[1]=ctxt_1->flag|0x20|0x80;
					}
				else 	{
					flag[0]=ctxt_0->flag|0x1|0x40;
					flag[1]=ctxt_1->flag|0x1|0x80;
					}
				if(ctxt_1->pe_pos.a){
					flag[0]=flag[0]|0x20|0x1|0x40;
					flag[1]=flag[1]|0x1|0x80|0x10;
					}
				else	{
					flag[0]=flag[0]|0x1|0x40;
					flag[1]=flag[1]|0x1|0x80;
					}
				nn1 = refdb_seq_amb(refdb,ctxt_0->pe_pos.pos,ctxt_0->seq->length);
				nn2 = refdb_seq_amb(refdb,ctxt_1->pe_pos.pos,ctxt_1->seq->length);
				//reads 1 info
				fprintf(fout,"%s\t%u\t",ctxt_0->seq->name,flag[0]);
                                fprintf(fout,"%s\t%u\t",sam_info0->ref_name,ctxt_0->pe_pos.pos);
                                fprintf(fout,"%u\t%s\t",ctxt_0->mapQ,sam_info0->cigar);
				if(nn1 ==nn2)
                        	        fprintf(fout,"=\t%u\t%d\t",ctxt_1->pe_pos.pos,ctxt_1->pe_pos.pos+cor_pos2-ctxt_0->pe_pos.pos-cor_pos1);
				else	fprintf(fout,"%s\t%u\t0\t",sam_info1->ref_name,ctxt_1->pe_pos.pos);

				if(ctxt_0->pe_pos.a)
                                        fprintf(fout,"%s\t%s\t",ctxt_0->seq_comp,ctxt_0->qual_comp);
                                else    fprintf(fout,"%s\t%s\t",ctxt_0->seq->seq,ctxt_0->seq->qual);


                                if(ctxt_0->c0 > 1)
                                        fprintf(fout,"XT:A:R\t");
                                else
                                        fprintf(fout,"XT:A:U\t");
                                fprintf(fout,"NM:i:%u\t",ctxt_0->pe_pos.n_mm);
                                if(sam_info0->ref_amb)
                                fprintf(fout,"XN:i:%u\t",sam_info0->ref_amb);
					fprintf(fout,"SM:i:%u\tAM:i:%u\t",ctxt_0->seQ,am);	
					fprintf(fout,"X0:i:%u\tX1:i:%u\t",ctxt_0->c0,ctxt_0->c1-ctxt_0->c0);
                                        fprintf(fout,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_0->pe_pos.n_mm,ctxt_0->pe_pos.n_gapo,ctxt_0->pe_pos.n_gapo+ctxt_0->pe_pos.n_gape);

                                        fprintf(fout,"MD:Z:%s",sam_info0->md);
					for(i=0;i<sam_info0->num;i++){
                                                pos_t  *pos = ctxt_0->pos+i*sam_info0->step;
						if(sam_info0->repeat){
							if(ctxt_0->pe_pos.pos!=pos->pos){
								if(j==0)
       	                                         		fprintf(fout,"\tXA:Z:");
								if(j<3)
                                                        	fprintf(fout,"%s,%c%u,%s,%u;",sam_info0->m_ref_name[i],pos->a?'-':'+',pos->pos,sam_info0->m_cigar[i],pos->n_mm);
								j++;
							}
						}
						FREE(sam_info0->m_cigar[i]);
						FREE(sam_info0->m_ref_name[i]);
					}
                                fprintf(fout,"\n");
				j=0;

				//reads 2 info
				fprintf(fout,"%s\t%u\t",ctxt_1->seq->name,flag[1]);
                                fprintf(fout,"%s\t%u\t",sam_info1->ref_name,ctxt_1->pe_pos.pos);
                                fprintf(fout,"%u\t%s\t",ctxt_1->mapQ,sam_info1->cigar);
                                if(nn1==nn2)
					fprintf(fout,"=\t%u\t%d\t",ctxt_0->pe_pos.pos,ctxt_0->pe_pos.pos+cor_pos1-ctxt_1->pe_pos.pos-cor_pos2);
				else	fprintf(fout,"%s\t%u\t0\t",sam_info0->ref_name,ctxt_0->pe_pos.pos);

				if(ctxt_1->pe_pos.a)
                                        fprintf(fout,"%s\t%s\t",ctxt_1->seq_comp,ctxt_1->qual_comp);
                                else    fprintf(fout,"%s\t%s\t",ctxt_1->seq->seq,ctxt_1->seq->qual);

				//reads2 info
                               	if(ctxt_1->c0 >1)
					fprintf(fout,"XT:A:R\t");
                                else
                                        fprintf(fout,"XT:A:U\t");
                                fprintf(fout,"NM:i:%u\t",ctxt_1->pe_pos.n_mm);
                                if(sam_info1->ref_amb)
                                        fprintf(fout,"XN:i:%u\t",sam_info1->ref_amb);
				fprintf(fout,"SM:i:%u\tAM:i:%u\t",ctxt_1->seQ,am);	

                                fprintf(fout,"X0:i:%u\tX1:i:%u\t",ctxt_1->c0,ctxt_1->c1-ctxt_1->c0);
                                fprintf(fout,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_1->pe_pos.n_mm,ctxt_1->pe_pos.n_gapo,ctxt_1->pe_pos.n_gapo+ctxt_1->pe_pos.n_gape);
				fprintf(fout,"MD:Z:%s",sam_info1->md);
				for(i=0;i<sam_info1->num;i++){
                                        pos_t  *pos = ctxt_1->pos+i*sam_info1->step;
                                	if(sam_info1->repeat){
						if(ctxt_1->pe_pos.pos!=pos->pos){
                                        	if(j==0)
						fprintf(fout,"\tXA:Z:");
						if(j<3)
                                        	fprintf(fout,"%s,%c%u,%s,%u;",sam_info1->m_ref_name[i],pos->a?'-':'+',pos->pos,sam_info1->m_cigar[i],pos->n_mm);
						j++;
						}
					}
					FREE(sam_info1->m_cigar[i]);
					FREE(sam_info1->m_ref_name[i]);
				}
                                fprintf(fout,"\n");  
				seqdb_release(seqdb,ctxt_0->seq);
				seqdb_release(seqdb,ctxt_1->seq);
 
				FREE(ctxt_0->pos);
				FREE(ctxt_1->pos);
				FREE(ctxt_0->qual_comp);
                                FREE(ctxt_0->seq_comp);
				FREE(ctxt_1->qual_comp);
				FREE(ctxt_1->seq_comp);
				FREE(sam_info0->ref_name);
				FREE(sam_info0->md);
				FREE(sam_info0->cigar);
				FREE(sam_info1->ref_name);
				FREE(sam_info1->md);
				FREE(sam_info1->cigar);

				break;
				}
			case(4):{
				seq_sz_t pos;
				sw2cigar(cigar,ctxt_0->cigar,ctxt_0->n_cigar);
				am = (ctxt_1->seQ>ctxt_0->seQ?ctxt_0->seQ:ctxt_1->seQ);
		//		getPosValue(ctxt_0);
				uint32_t cor_pos1,cor_pos2;
				int	offset = 0 ;

				ref_id = refdb_idx(refdb,ctxt_0->pe_pos.pos);
				offset = refdb_offset(refdb,ref_id);
                                cor_pos1 = cor_pos2 = 0;
                                if(ctxt_0->pe_pos.a)
                                cor_pos1= getLength(cigar);

                                if(ctxt_1->pe_pos.a)
                                cor_pos2= getLength(sam_info1->cigar);

				pos=ctxt_0->pe_pos.pos+1 - offset;

				if(ctxt_0->pe_pos.a){
                                        flag[0]=ctxt_0->flag|0x1|0x10|0x40;
                                        flag[1]=ctxt_1->flag|0x20|0x80;
                                        }
                                else    {
                                        flag[0]=ctxt_0->flag|0x1|0x40;
                                        flag[1]=ctxt_1->flag|0x1|0x80;
                                        }
                                if(ctxt_1->pe_pos.a){
                                        flag[0]=flag[0]|0x20|0x1|0x40;
                                        flag[1]=flag[1]|0x1|0x80|0x10;
                                        }
                                else    {
                                        flag[0]=flag[0]|0x1|0x40;
                                        flag[1]=flag[1]|0x1|0x80;
                                        }

				//	reads1
				fprintf(fout,"%s\t%u\t",ctxt_0->seq->name,flag[0]);
                                fprintf(fout,"%s\t%u\t",sam_info1->ref_name,pos);
                                fprintf(fout,"%u\t%s\t",ctxt_0->mapQ,cigar);
                                fprintf(fout,"=\t%u\t%d\t",ctxt_1->pe_pos.pos,ctxt_1->pe_pos.pos+cor_pos2-pos-cor_pos1);
				
				
				ctxt_0->seq_comp_sz=0;
				ctxt_0->qual_comp_sz=0;
				realloc_n(ctxt_0->seq_comp,ctxt_0->seq_comp_sz, ctxt_0->seq->length+1);
				realloc_n(ctxt_0->qual_comp,ctxt_0->qual_comp_sz, ctxt_0->seq->length+1);
				seq_to_comp(ctxt_0->seq_comp,ctxt_0->seq->seq,ctxt_0->seq->length);
				ctxt_0->seq_comp[ctxt_0->seq->length]=0;
				reverse(ctxt_0->qual_comp,ctxt_0->seq->qual,ctxt_0->seq->length);
				ctxt_0->qual_comp[ctxt_0->seq->length]=0;

				if(ctxt_0->pe_pos.a)
                                        fprintf(fout,"%s\t%s\t",ctxt_0->seq_comp,ctxt_0->qual_comp);
                                else    fprintf(fout,"%s\t%s\t",ctxt_0->seq->seq,ctxt_0->seq->qual);
					//reads2 info
				fprintf(fout,"XT:A:M\t");
				fprintf(fout,"NM:i:%u\t",ctxt_0->pe_pos.n_mm);
				if(sam_info1->ref_amb)
                                        fprintf(fout,"XN:i:%u\t",sam_info1->ref_amb);
				fprintf(fout,"SM:i:%u\tAM:i:%u\t",ctxt_0->seQ,am);	
                                fprintf(fout,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_0->pe_pos.n_mm,ctxt_0->pe_pos.n_gapo,ctxt_0->pe_pos.n_gapo+ctxt_0->pe_pos.n_gape);
				fprintf(fout,"MD:Z:%s",ctxt_0->md);
				fprintf(fout,"\n");
			
				//reads2	
				fprintf(fout,"%s\t%u\t",ctxt_1->seq->name,flag[1]);
                                fprintf(fout,"%s\t%u\t",sam_info1->ref_name,ctxt_1->pe_pos.pos);
                                fprintf(fout,"%u\t%s\t",ctxt_1->mapQ,sam_info1->cigar);
                                fprintf(fout,"=\t%u\t%d\t",pos,pos+cor_pos1-ctxt_1->pe_pos.pos-cor_pos2);
                                if(ctxt_1->pe_pos.a)
                                        fprintf(fout,"%s\t%s\t",ctxt_1->seq_comp,ctxt_1->qual_comp);
                                else    fprintf(fout,"%s\t%s\t",ctxt_1->seq->seq,ctxt_1->seq->qual);


                                if(ctxt_1->c0>1)
                                        fprintf(fout,"XT:A:R\t");
                                else
                                        fprintf(fout,"XT:A:U\t");
				
                                fprintf(fout,"NM:i:%u\t",ctxt_1->pe_pos.n_mm);
                                if(sam_info1->ref_amb)
                                        fprintf(fout,"XN:i:%u\t",sam_info1->ref_amb);
				fprintf(fout,"SM:i:%u\tAM:i:%u\t",ctxt_1->seQ,am);	
				fprintf(fout,"X0:i:%u\tX1:i:%u\t",ctxt_1->c0,ctxt_1->c1-ctxt_0->c0);

                                fprintf(fout,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_1->pe_pos.n_mm,ctxt_1->pe_pos.n_gapo,ctxt_1->pe_pos.n_gapo+ctxt_1->pe_pos.n_gape);
                                fprintf(fout,"MD:Z:%s",sam_info1->md);
				for(i=0;i<sam_info1->num;i++){
                                        pos_t  *pos1 = ctxt_1->pos+i*sam_info1->step;
					if(sam_info1->repeat){
						if(pos1->pos!=ctxt_1->pe_pos.pos){
							if(j==0)
                                        		fprintf(fout,"\tXA:Z:");
							if(j<3)
                        	                	fprintf(fout,"%s,%c%u,%s,%u;",sam_info1->m_ref_name[i],pos1->a?'-':'+',pos1->pos,sam_info1->m_cigar[i],pos1->n_mm);
							j++;
						}
					}
                                        FREE(sam_info1->m_cigar[i]);
					FREE(sam_info1->m_ref_name[i]);
				}
				j=0;
					fprintf(fout,"\n");

				

				FREE(ctxt_0->seq_comp);
				FREE(ctxt_0->qual_comp);
				FREE(ctxt_0->pos);
				FREE(ctxt_1->pos);
				FREE(ctxt_0->cigar);
				FREE(ctxt_0->md);
				seqdb_release(seqdb,ctxt_0->seq);
				seqdb_release(seqdb,ctxt_1->seq);
				FREE(ctxt_1->qual_comp);
				FREE(ctxt_1->seq_comp);
				FREE(sam_info1->md);
				FREE(sam_info1->cigar);
				FREE(sam_info1->ref_name);

				break;
				}


			case(5):
				{
				int	offset = 0 ;
				am = (ctxt_1->seQ>ctxt_0->seQ?ctxt_0->seQ:ctxt_1->seQ);
				seq_sz_t pos;
				ref_id = refdb_idx(refdb,ctxt_1->pe_pos.pos);
				offset = refdb_offset(refdb,ref_id);
                                pos=ctxt_1->pe_pos.pos+1 -offset;
				sw2cigar(cigar,ctxt_1->cigar,ctxt_1->n_cigar);
				uint32_t cor_pos1,cor_pos2;

                                cor_pos1 = cor_pos2 = 0;
                                
				if(ctxt_0->pe_pos.a)
                                cor_pos1= getLength(sam_info0->cigar);

                                if(ctxt_1->pe_pos.a)
                                cor_pos2= getLength(cigar);

				if(ctxt_0->pe_pos.a){
                                        flag[0]=ctxt_0->flag|0x1|0x10|0x40;
                                        flag[1]=ctxt_1->flag|0x20|0x80;
                                        }
                                else    {
                                        flag[0]=ctxt_0->flag|0x1|0x40;
                                        flag[1]=ctxt_1->flag|0x1|0x80;
                                        }
                                if(ctxt_1->pe_pos.a){
                                        flag[0]=flag[0]|0x20|0x1|0x40;
                                        flag[1]=flag[1]|0x1|0x80|0x10;
                                        }
                                else    {
                                        flag[0]=flag[0]|0x1|0x40;
                                        flag[1]=flag[1]|0x1|0x80;
                                        }
				//reads1	
				fprintf(fout,"%s\t%u\t",ctxt_0->seq->name,flag[0]);
                                fprintf(fout,"%s\t%u\t",sam_info0->ref_name,ctxt_0->pe_pos.pos);
                                fprintf(fout,"%u\t%s\t",ctxt_0->mapQ,sam_info0->cigar);
                                fprintf(fout,"=\t%u\t%d\t",pos,pos+cor_pos2-ctxt_0->pe_pos.pos-cor_pos1);
                                if(ctxt_0->pos->a)
                                        fprintf(fout,"%s\t%s\t",ctxt_0->seq_comp,ctxt_0->qual_comp);
                                else    fprintf(fout,"%s\t%s\t",ctxt_0->seq->seq,ctxt_0->seq->qual);
				
				if(ctxt_0->c0 >1 )
                                	fprintf(fout,"XT:A:R\t");
                                else	fprintf(fout,"XT:A:U\t");
                                        fprintf(fout,"NM:i:%u\t",ctxt_0->pe_pos.n_mm);
                                if(sam_info0->ref_amb)
                                        fprintf(fout,"XN:i:%u\t",sam_info0->ref_amb);
				fprintf(fout,"SM:i:%u\tAM:i:%u\t",ctxt_0->seQ,am);
				fprintf(fout,"X0:i:%u\tX1:i:%u\t",ctxt_0->c0,ctxt_0->c1-ctxt_0->c0);
	

                                        fprintf(fout,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_0->pe_pos.n_mm,ctxt_0->pe_pos.n_gapo,ctxt_0->pe_pos.n_gapo+ctxt_0->pe_pos.n_gape);

                                        fprintf(fout,"MD:Z:%s",sam_info0->md);	
					for(i=0;i<sam_info0->num;i++){
                                        	pos_t  *pos1 = ctxt_0->pos+i*sam_info0->step;
						if(sam_info0->repeat){
							if(ctxt_0->pe_pos.pos!=pos1->pos){
								if(j==0)
                                        			fprintf(fout,"\tXA:Z:");
								if(j<3)
                                        			fprintf(fout,"%s,%c%u,%s,%u;",sam_info0->m_ref_name[i],pos1->a?'-':'+',pos1->pos,sam_info0->m_cigar[i],pos1->n_mm);
								j++;
							}
						}
                                        	FREE(sam_info0->m_cigar[i]);
						FREE(sam_info0->m_ref_name[i]);
                               	}
				j=0;
				fprintf(fout,"\n");

			//reads2
				ctxt_1->seq_comp_sz=0;
                                ctxt_1->qual_comp_sz=0;
                                realloc_n(ctxt_1->seq_comp,ctxt_1->seq_comp_sz, ctxt_1->seq->length+1);
                                realloc_n(ctxt_1->qual_comp,ctxt_1->qual_comp_sz, ctxt_1->seq->length+1);
                                seq_to_comp(ctxt_1->seq_comp,ctxt_1->seq->seq,ctxt_1->seq->length);
                                ctxt_1->seq_comp[ctxt_1->seq->length]=0;
                                reverse(ctxt_1->qual_comp,ctxt_1->seq->qual,ctxt_1->seq->length);
                                ctxt_1->qual_comp[ctxt_1->seq->length]=0;
				fprintf(fout,"%s\t%u\t",ctxt_1->seq->name,flag[1]);
                                fprintf(fout,"%s\t%u\t",sam_info0->ref_name,pos);
                                fprintf(fout,"%u\t%s\t",ctxt_1->mapQ,cigar);
                                fprintf(fout,"=\t%u\t%d\t",ctxt_0->pe_pos.pos,ctxt_0->pe_pos.pos+cor_pos1-pos-cor_pos2);
				if(ctxt_1->pe_pos.a)
                                        fprintf(fout,"%s\t%s\t",ctxt_1->seq_comp,ctxt_1->qual_comp);
                                else    fprintf(fout,"%s\t%s\t",ctxt_1->seq->seq,ctxt_1->seq->qual);
                     
				
                                fprintf(fout,"XT:A:M\t");
                                fprintf(fout,"NM:i:%u\t",ctxt_1->pe_pos.n_mm);
				if(sam_info0->ref_amb)
                                        fprintf(fout,"XN:i:%u\t",sam_info1->ref_amb);
				fprintf(fout,"SM:i:%u\tAM:i:%u\t",ctxt_1->seQ,am);
                                fprintf(fout,"X0:i:%u\tX1:i:%u\t",ctxt_1->c0,ctxt_1->c1-ctxt_1->c0);
                                fprintf(fout,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_1->pe_pos.n_mm,ctxt_1->pe_pos.n_gapo,ctxt_1->pe_pos.n_gapo+ctxt_1->pe_pos.n_gape);
                                fprintf(fout,"MD:Z:%s",ctxt_1->md);
                                fprintf(fout,"\n");

				FREE(ctxt_0->pos);
				FREE(ctxt_1->pos);
                                FREE(ctxt_1->cigar);
                                FREE(ctxt_1->md);
				FREE(ctxt_1->seq_comp);
				FREE(ctxt_1->qual_comp);
				seqdb_release(seqdb,ctxt_0->seq);
				seqdb_release(seqdb,ctxt_1->seq);
                                FREE(ctxt_0->qual_comp);
                                FREE(ctxt_0->seq_comp);
                                FREE(sam_info0->md);
                                FREE(sam_info0->cigar);
                                FREE(sam_info0->ref_name);


				break;
				}
			default:break;
			}
		return 0;
}

int fout_sam_info1(sam_info_t *sam_info0,sam_info_t *sam_info1,context_t *ctxt_0, context_t *ctxt_1,int slt ,char *sam_buf  , refdb_t refdb)
{

			int flag[2];
			char	cigar[4086];
			int	ref_id = 0;
			unsigned int am;
			int	i;
			int	j=0;
			seqdb_t seqdb = NULL;
			switch(slt){
			case(0)://all ummap!!
					{	
					flag[0] =0x40|0x4|0x8|0x1;
					flag[1] =0x80|0x4|0x8|0x1;	
					//reads1 info
					sam_buf += sprintf(sam_buf,"%s\t%u\t*\t0\t0\t*\t*\t0\t0\t",ctxt_0->seq->name,flag[0]);
					sam_buf += sprintf(sam_buf,"%s\t%s",ctxt_0->seq->seq,ctxt_0->seq->qual);
					sam_buf += sprintf(sam_buf,"\n");
					//reads2 info
					sam_buf += sprintf(sam_buf,"%s\t%u\t*\t0\t0\t*\t*\t0\t0\t",ctxt_1->seq->name,flag[1]);
					sam_buf += sprintf(sam_buf,"%s\t%s",ctxt_1->seq->seq,ctxt_1->seq->qual);
					sam_buf += sprintf(sam_buf,"\n");
					//reads2 info
				//free
					seqdb_release(seqdb,ctxt_0->seq);
					seqdb_release(seqdb,ctxt_1->seq);
					break;
				}
			case(1)://reads 0 ummap!!
				{
					flag[0] =0x40|0x4|0x1;
					if(ctxt_1->pe_pos.a){
						flag[1] =ctxt_1->flag|0x80|0x8|0x1|0x10;
						flag[0] |= 0x10|0x20;
						}
					else	flag[1] =ctxt_1->flag|0x80|0x8|0x1;
					//reads1 info
					sam_buf += sprintf(sam_buf,"%s\t%u\t",ctxt_0->seq->name,flag[0]);
					sam_buf += sprintf(sam_buf,"%s\t%u\t",sam_info1->ref_name,ctxt_1->pe_pos.pos);
					sam_buf += sprintf(sam_buf,"0\t*\t=\t%u\t0\t",ctxt_1->pe_pos.pos);
					sam_buf += sprintf(sam_buf,"%s\t%s",ctxt_0->seq->seq,ctxt_0->seq->qual);
					sam_buf += sprintf(sam_buf,"\n");
					
					//reads2 info
					sam_buf += sprintf(sam_buf,"%s\t%u\t", ctxt_1->seq->name,flag[1]);
					sam_buf += sprintf(sam_buf,"%s\t%u\t", sam_info1->ref_name,ctxt_1->pe_pos.pos);
					sam_buf += sprintf(sam_buf,"%u\t%s\t", ctxt_1->mapQ,sam_info1->cigar);
					sam_buf += sprintf(sam_buf,"=\t%u\t0\t",ctxt_1->pe_pos.pos);
					// reads2  seq and qual
					if(ctxt_1->pos->a)
						sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_1->seq_comp,ctxt_1->qual_comp);
					else	sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_1->seq->seq,ctxt_1->seq->qual);

					
					if(ctxt_1->c0 > 1)
						sam_buf+=sprintf(sam_buf,"XT:A:R\t");
					else
						sam_buf+=sprintf(sam_buf,"XT:A:U\t");
					sam_buf += sprintf(sam_buf,"NM:i:%u\t",ctxt_1->pe_pos.n_mm);
					if(sam_info1->ref_amb)
					sam_buf += sprintf(sam_buf,"XN:i:%u\t",sam_info1->ref_amb);
					sam_buf += sprintf(sam_buf,"SM:i:%u\tAM:i:0\t",ctxt_1->seQ);	
					
					sam_buf += sprintf(sam_buf,"X0:i:%u\tX1:i:%u\t",ctxt_1->c0,ctxt_1->c1-ctxt_1->c0);
					sam_buf += sprintf(sam_buf,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_1->pe_pos.n_mm,ctxt_1->pe_pos.n_gapo,ctxt_1->pe_pos.n_gapo+ctxt_1->pe_pos.n_gape);

					sam_buf += sprintf(sam_buf,"MD:Z:%s",sam_info1->md);
					for(i=0;i<sam_info1->num;i++){
						pos_t  *pos = ctxt_1->pos+i*sam_info1->step;
						if(sam_info1->repeat){	
							if(ctxt_1->pe_pos.pos!=pos->pos){
								if( j == 0 )
								sam_buf += sprintf(sam_buf,"\tXA:Z:");
								if( j < 3 )
								sam_buf += sprintf(sam_buf,"%s,%c%u,%s,%u;",sam_info1->m_ref_name[i],pos->a?'-':'+',pos->pos,sam_info1->m_cigar[i],pos->n_mm);
								j++;
							}
						}
						FREE(sam_info1->m_cigar[i]);
						FREE(sam_info1->m_ref_name[i]);
					}
					sam_buf += sprintf(sam_buf,"\n");
					seqdb_release(seqdb,ctxt_0->seq);
					FREE(ctxt_1->pos);
					seqdb_release(seqdb,ctxt_1->seq);
					FREE(ctxt_1->qual_comp);
					FREE(ctxt_1->seq_comp);
					FREE(sam_info1->cigar);
					FREE(sam_info1->md);
					FREE(sam_info1->ref_name);
					break;
					}
			case(2):{
					flag[1] =0x80|0x4|0x1;
					if(ctxt_0->pe_pos.a){
						flag[0]=0x01|0x08|0x10|0x40;
						flag[1] |=0x10|0x20;
						}
					else
						flag[0]=0x01|0x8|0x40;
					//reads 1 info;
					sam_buf += sprintf(sam_buf,"%s\t%u\t",ctxt_0->seq->name,flag[0]);
                                        sam_buf += sprintf(sam_buf,"%s\t%u\t",sam_info0->ref_name,ctxt_0->pe_pos.pos);
                                        sam_buf += sprintf(sam_buf,"%u\t%s\t",ctxt_0->mapQ,sam_info0->cigar);
                                        sam_buf += sprintf(sam_buf,"=\t%u\t0\t",ctxt_0->pe_pos.pos);

					if(ctxt_0->pos->a)
                                                sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_0->seq_comp,ctxt_0->qual_comp);
                                        else    sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_0->seq->seq,ctxt_0->seq->qual);


                                        if(ctxt_0->c0 >1)
                                                sam_buf += sprintf(sam_buf,"XT:A:R\t");
                                        else
                                                sam_buf += sprintf(sam_buf,"XT:A:U\t");
                                        sam_buf += sprintf(sam_buf,"NM:i:%u\t",ctxt_0->pe_pos.n_mm);
                                        if(sam_info0->ref_amb)
                                        sam_buf += sprintf(sam_buf,"XN:i:%u\t",sam_info0->ref_amb);
					
					sam_buf += sprintf(sam_buf,"SM:i:%u\tAM:i:0\t",ctxt_0->seQ);	
					
					sam_buf += sprintf(sam_buf,"X0:i:%u\tX1:i:%u\t",ctxt_0->c0,ctxt_0->c1-ctxt_0->c0);
                                        sam_buf += sprintf(sam_buf,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_0->pe_pos.n_mm,ctxt_0->pe_pos.n_gapo,ctxt_0->pe_pos.n_gapo+ctxt_0->pe_pos.n_gape);

                                        sam_buf += sprintf(sam_buf,"MD:Z:%s",sam_info0->md);
					for(i=0;i<sam_info0->num;i++){
                                                pos_t  *pos = ctxt_0->pos+i*sam_info0->step;
						if(sam_info0->repeat){
							if(pos->pos!=ctxt_0->pe_pos.pos){
								if( j == 0)
								sam_buf += sprintf(sam_buf,"\tXA:Z:");
								if( j < 3)
                                                        	sam_buf += sprintf(sam_buf,"%s,%c%u,%s,%u;",sam_info0->m_ref_name[i],pos->a?'-':'+',pos->pos,sam_info0->m_cigar[i],pos->n_mm);
								j++;
							}
											
						}
						FREE(sam_info0->m_cigar[i]);
						FREE(sam_info0->m_ref_name[i]);
					}
                                        sam_buf += sprintf(sam_buf,"\n");
					// reads 2 info;
					sam_buf += sprintf(sam_buf,"%s\t%u\t",ctxt_1->seq->name,flag[1]);
					sam_buf += sprintf(sam_buf,"%s\t%u\t0\t*\t",sam_info0->ref_name,ctxt_0->pe_pos.pos);
					sam_buf += sprintf(sam_buf,"=\t%d\t0\t",ctxt_0->pe_pos.pos);
                                        sam_buf += sprintf(sam_buf,"%s\t%s",ctxt_1->seq->seq,ctxt_1->seq->qual);
					sam_buf += sprintf(sam_buf,"\n");

					seqdb_release(seqdb,ctxt_0->seq);
                                        FREE(ctxt_0->pos);
					seqdb_release(seqdb,ctxt_1->seq);
                                        FREE(ctxt_0->qual_comp);
                                        FREE(ctxt_0->seq_comp);
                                        FREE(sam_info0->cigar);
                                        FREE(sam_info0->md);
                                        FREE(sam_info0->ref_name);

					break;
					
					

				}
			case(3):
				{
				uint32_t cor_pos1,cor_pos2;
				unsigned int nn1,nn2;
				cor_pos1 = cor_pos2 = 0;
				am = (ctxt_1->seQ>ctxt_0->seQ?ctxt_0->seQ:ctxt_1->seQ);
				if(ctxt_0->pe_pos.a)
				cor_pos1= getLength(sam_info0->cigar);
				
				if(ctxt_1->pe_pos.a)
                                cor_pos2= getLength(sam_info1->cigar);
						
				if(ctxt_0->pe_pos.a){
					flag[0]=ctxt_0->flag|0x1|0x10|0x40;
					flag[1]=ctxt_1->flag|0x20|0x80;
					}
				else 	{
					flag[0]=ctxt_0->flag|0x1|0x40;
					flag[1]=ctxt_1->flag|0x1|0x80;
					}
				if(ctxt_1->pe_pos.a){
					flag[0]=flag[0]|0x20|0x1|0x40;
					flag[1]=flag[1]|0x1|0x80|0x10;
					}
				else	{
					flag[0]=flag[0]|0x1|0x40;
					flag[1]=flag[1]|0x1|0x80;
					}
				nn1 = refdb_seq_amb(refdb,ctxt_0->pe_pos.pos,ctxt_0->seq->length);
				nn2 = refdb_seq_amb(refdb,ctxt_1->pe_pos.pos,ctxt_1->seq->length);
				//reads 1 info
				sam_buf += sprintf(sam_buf,"%s\t%u\t",ctxt_0->seq->name,flag[0]);
                                sam_buf += sprintf(sam_buf,"%s\t%u\t",sam_info0->ref_name,ctxt_0->pe_pos.pos);
                                sam_buf += sprintf(sam_buf,"%u\t%s\t",ctxt_0->mapQ,sam_info0->cigar);
				if(nn1 ==nn2)
					sam_buf += sprintf(sam_buf,"=\t%u\t%d\t",ctxt_1->pe_pos.pos,ctxt_1->pe_pos.pos+cor_pos2-ctxt_0->pe_pos.pos-cor_pos1);
				else	sam_buf += sprintf(sam_buf,"%s\t%u\t0\t",sam_info1->ref_name,ctxt_1->pe_pos.pos);

				if(ctxt_0->pe_pos.a)
                                        sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_0->seq_comp,ctxt_0->qual_comp);
                                else    sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_0->seq->seq,ctxt_0->seq->qual);


                                if(ctxt_0->c0 > 1)
                                        sam_buf += sprintf(sam_buf,"XT:A:R\t");
                                else
                                        sam_buf += sprintf(sam_buf,"XT:A:U\t");
                                sam_buf += sprintf(sam_buf,"NM:i:%u\t",ctxt_0->pe_pos.n_mm);
                                if(sam_info0->ref_amb)
                                        sam_buf += sprintf(sam_buf,"XN:i:%u\t",sam_info0->ref_amb);
				sam_buf += sprintf(sam_buf,"SM:i:%u\tAM:i:%u\t",ctxt_0->seQ,am);
					sam_buf += sprintf(sam_buf,"X0:i:%u\tX1:i:%u\t",ctxt_0->c0,ctxt_0->c1-ctxt_0->c0);
                                        sam_buf += sprintf(sam_buf,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_0->pe_pos.n_mm,ctxt_0->pe_pos.n_gapo,ctxt_0->pe_pos.n_gapo+ctxt_0->pe_pos.n_gape);

                                        sam_buf += sprintf(sam_buf,"MD:Z:%s",sam_info0->md);
					for(i=0;i<sam_info0->num;i++){
                                                pos_t  *pos = ctxt_0->pos+i*sam_info0->step;
						if(sam_info0->repeat){
							if(ctxt_0->pe_pos.pos!=pos->pos){
								if(j==0)
       	                                         		sam_buf += sprintf(sam_buf,"\tXA:Z:");
								if(j<3)
                                                        	sam_buf += sprintf(sam_buf,"%s,%c%u,%s,%u;",sam_info0->m_ref_name[i],pos->a?'-':'+',pos->pos,sam_info0->m_cigar[i],pos->n_mm);
								j++;
							}
						}
						FREE(sam_info0->m_cigar[i]);
						FREE(sam_info0->m_ref_name[i]);
					}
                                sam_buf += sprintf(sam_buf,"\n");
				j=0;

				//reads 2 info
				sam_buf += sprintf(sam_buf,"%s\t%u\t",ctxt_1->seq->name,flag[1]);
                                sam_buf += sprintf(sam_buf,"%s\t%u\t",sam_info1->ref_name,ctxt_1->pe_pos.pos);
                                sam_buf += sprintf(sam_buf,"%u\t%s\t",ctxt_1->mapQ,sam_info1->cigar);
                                if(nn1==nn2)
					sam_buf += sprintf(sam_buf,"=\t%u\t%d\t",ctxt_0->pe_pos.pos,ctxt_0->pe_pos.pos+cor_pos1-ctxt_1->pe_pos.pos-cor_pos2);
				else	sam_buf += sprintf(sam_buf,"%s\t%u\t0\t",sam_info0->ref_name,ctxt_0->pe_pos.pos);

				if(ctxt_1->pe_pos.a)
                                        sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_1->seq_comp,ctxt_1->qual_comp);
                                else    sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_1->seq->seq,ctxt_1->seq->qual);

				//reads2 info
                               	if(ctxt_1->c0 >1)
					sam_buf += sprintf(sam_buf,"XT:A:R\t");
                                else
                                        sam_buf += sprintf(sam_buf,"XT:A:U\t");
                                sam_buf += sprintf(sam_buf,"NM:i:%u\t",ctxt_1->pe_pos.n_mm);
                                if(sam_info1->ref_amb)
                                        sam_buf += sprintf(sam_buf,"XN:i:%u\t",sam_info1->ref_amb);
				sam_buf += sprintf(sam_buf,"SM:i:%u\tAM:i:%u\t",ctxt_1->seQ,am);

                                sam_buf += sprintf(sam_buf,"X0:i:%u\tX1:i:%u\t",ctxt_1->c0,ctxt_1->c1-ctxt_1->c0);
                                sam_buf += sprintf(sam_buf,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_1->pe_pos.n_mm,ctxt_1->pe_pos.n_gapo,ctxt_1->pe_pos.n_gapo+ctxt_1->pe_pos.n_gape);
				sam_buf += sprintf(sam_buf,"MD:Z:%s",sam_info1->md);
				for(i=0;i<sam_info1->num;i++){
                                        pos_t  *pos = ctxt_1->pos+i*sam_info1->step;
                                	if(sam_info1->repeat){
						if(ctxt_1->pe_pos.pos!=pos->pos){
                                        	if(j==0)
						sam_buf += sprintf(sam_buf,"\tXA:Z:");
						if(j<3)
                                        	sam_buf += sprintf(sam_buf,"%s,%c%u,%s,%u;",sam_info1->m_ref_name[i],pos->a?'-':'+',pos->pos,sam_info1->m_cigar[i],pos->n_mm);
						j++;
						}
					}
					FREE(sam_info1->m_cigar[i]);
					FREE(sam_info1->m_ref_name[i]);
				}
                                sam_buf += sprintf(sam_buf,"\n");

				FREE(ctxt_0->pos);
				FREE(ctxt_1->pos);
				seqdb_release(seqdb,ctxt_0->seq);
				seqdb_release(seqdb,ctxt_1->seq);
				FREE(ctxt_0->qual_comp);
                                FREE(ctxt_0->seq_comp);
				FREE(ctxt_1->qual_comp);
				FREE(ctxt_1->seq_comp);
				FREE(sam_info0->ref_name);
				FREE(sam_info0->md);
				FREE(sam_info0->cigar);
				FREE(sam_info1->ref_name);
				FREE(sam_info1->md);
				FREE(sam_info1->cigar);

				break;
				}
			case(4):{
				seq_sz_t pos;
				sw2cigar(cigar,ctxt_0->cigar,ctxt_0->n_cigar);
				am = (ctxt_1->seQ>ctxt_0->seQ?ctxt_0->seQ:ctxt_1->seQ);
		//		getPosValue(ctxt_0);
				uint32_t cor_pos1,cor_pos2;
				int	offset = 0 ;

				ref_id = refdb_idx(refdb,ctxt_0->pe_pos.pos);
				offset = refdb_offset(refdb,ref_id);

                                cor_pos1 = cor_pos2 = 0;
                                if(ctxt_0->pe_pos.a)
                                cor_pos1= getLength(cigar);

                                if(ctxt_1->pe_pos.a)
                                cor_pos2= getLength(sam_info1->cigar);

				pos=ctxt_0->pe_pos.pos + 1 - offset;

				if(ctxt_0->pe_pos.a){
                                        flag[0]=ctxt_0->flag|0x1|0x10|0x40;
                                        flag[1]=ctxt_1->flag|0x20|0x80;
                                        }
                                else    {
                                        flag[0]=ctxt_0->flag|0x1|0x40;
                                        flag[1]=ctxt_1->flag|0x1|0x80;
                                        }
                                if(ctxt_1->pe_pos.a){
                                        flag[0]=flag[0]|0x20|0x1|0x40;
                                        flag[1]=flag[1]|0x1|0x80|0x10;
                                        }
                                else    {
                                        flag[0]=flag[0]|0x1|0x40;
                                        flag[1]=flag[1]|0x1|0x80;
                                        }

				//	reads1
				sam_buf += sprintf(sam_buf,"%s\t%u\t",ctxt_0->seq->name,flag[0]);
                                sam_buf += sprintf(sam_buf,"%s\t%u\t",sam_info1->ref_name,pos);
                                sam_buf += sprintf(sam_buf,"%u\t%s\t",ctxt_0->mapQ,cigar);
                                sam_buf += sprintf(sam_buf,"=\t%u\t%d\t",ctxt_1->pe_pos.pos,ctxt_1->pe_pos.pos+cor_pos2-pos-cor_pos1);
				
				
				ctxt_0->seq_comp_sz=0;
				ctxt_0->qual_comp_sz=0;
				realloc_n(ctxt_0->seq_comp,ctxt_0->seq_comp_sz, ctxt_0->seq->length+1);
				realloc_n(ctxt_0->qual_comp,ctxt_0->qual_comp_sz, ctxt_0->seq->length+1);
				seq_to_comp(ctxt_0->seq_comp,ctxt_0->seq->seq,ctxt_0->seq->length);
				ctxt_0->seq_comp[ctxt_0->seq->length]=0;
				reverse(ctxt_0->qual_comp,ctxt_0->seq->qual,ctxt_0->seq->length);
				ctxt_0->qual_comp[ctxt_0->seq->length]=0;

				if(ctxt_0->pe_pos.a)
                                        sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_0->seq_comp,ctxt_0->qual_comp);
                                else    sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_0->seq->seq,ctxt_0->seq->qual);
					//reads2 info
				sam_buf += sprintf(sam_buf,"XT:A:M\t");
				sam_buf += sprintf(sam_buf,"NM:i:%u\t",ctxt_0->pe_pos.n_mm);
				if(sam_info1->ref_amb)
                                        sam_buf += sprintf(sam_buf,"XN:i:%u\t",sam_info1->ref_amb);
				sam_buf += sprintf(sam_buf,"SM:i:%u\tAM:i:%u\t",ctxt_0->seQ,am);
				//sam_buf += sprintf(sam_buf,"XO:i:%u\tX1:i:%u\t",ctxt_0->c0,ctxt_0->c1-ctxt_0->c0);
                                sam_buf += sprintf(sam_buf,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_0->pe_pos.n_mm,ctxt_0->pe_pos.n_gapo,ctxt_0->pe_pos.n_gapo+ctxt_0->pe_pos.n_gape);
				sam_buf += sprintf(sam_buf,"MD:Z:%s",ctxt_0->md);
				sam_buf += sprintf(sam_buf,"\n");
			
				//reads2	
				sam_buf += sprintf(sam_buf,"%s\t%u\t",ctxt_1->seq->name,flag[1]);
                                sam_buf += sprintf(sam_buf,"%s\t%u\t",sam_info1->ref_name,ctxt_1->pe_pos.pos);
                                sam_buf += sprintf(sam_buf,"%u\t%s\t",ctxt_1->mapQ,sam_info1->cigar);
                                sam_buf += sprintf(sam_buf,"=\t%u\t%d\t",pos,pos+cor_pos1-ctxt_1->pe_pos.pos-cor_pos2);
                                if(ctxt_1->pe_pos.a)
                                        sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_1->seq_comp,ctxt_1->qual_comp);
                                else    sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_1->seq->seq,ctxt_1->seq->qual);


                                if(ctxt_1->c0>1)
                                        sam_buf += sprintf(sam_buf,"XT:A:R\t");
                                else
                                        sam_buf += sprintf(sam_buf,"XT:A:U\t");
				
                                sam_buf += sprintf(sam_buf,"NM:i:%u\t",ctxt_1->pe_pos.n_mm);
                                if(sam_info1->ref_amb)
                                        sam_buf += sprintf(sam_buf,"XN:i:%u\t",sam_info1->ref_amb);
				sam_buf += sprintf(sam_buf,"SM:i:%u\tAM:i:%u\t",ctxt_1->seQ,am);

                                sam_buf += sprintf(sam_buf,"X0:i:%u\tX1:i:%u\t",ctxt_1->c0,ctxt_1->c1-ctxt_1->c0);
                                sam_buf += sprintf(sam_buf,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_1->pe_pos.n_mm,ctxt_1->pe_pos.n_gapo,ctxt_1->pe_pos.n_gapo+ctxt_1->pe_pos.n_gape);
                                sam_buf += sprintf(sam_buf,"MD:Z:%s",sam_info1->md);
				for(i=0;i<sam_info1->num;i++){
                                        pos_t  *pos1 = ctxt_1->pos+i*sam_info1->step;
					if(sam_info1->repeat){
						if(pos1->pos!=ctxt_1->pe_pos.pos){
							if(j==0)
                                        		sam_buf += sprintf(sam_buf,"\tXA:Z:");
							if(j<3)
                        	                	sam_buf += sprintf(sam_buf,"%s,%c%u,%s,%u;",sam_info1->m_ref_name[i],pos1->a?'-':'+',pos1->pos,sam_info1->m_cigar[i],pos1->n_mm);
							j++;
						}
					}
                                        FREE(sam_info1->m_cigar[i]);
					FREE(sam_info1->m_ref_name[i]);
				}
				j=0;
					sam_buf += sprintf(sam_buf,"\n");
				

				FREE(ctxt_0->seq_comp);
				FREE(ctxt_0->qual_comp);
				FREE(ctxt_0->pos);
				FREE(ctxt_1->pos);
				FREE(ctxt_0->cigar);
				FREE(ctxt_0->md);
				seqdb_release(seqdb,ctxt_1->seq);
				seqdb_release(seqdb,ctxt_0->seq);
				FREE(ctxt_1->qual_comp);
				FREE(ctxt_1->seq_comp);
				FREE(sam_info1->md);
				FREE(sam_info1->cigar);
				FREE(sam_info1->ref_name);

				break;
				}


			case(5):
				{
				int     offset = 0;
				seq_sz_t pos;
				am = (ctxt_1->seQ>ctxt_0->seQ?ctxt_0->seQ:ctxt_1->seQ);
				ref_id = refdb_idx(refdb,ctxt_1->pe_pos.pos);
				offset = refdb_offset(refdb,ref_id);
                                pos = ctxt_1->pe_pos.pos + 1 - offset;
				sw2cigar(cigar,ctxt_1->cigar,ctxt_1->n_cigar);
				uint32_t cor_pos1,cor_pos2;

                                cor_pos1 = cor_pos2 = 0;
                                
				if(ctxt_0->pe_pos.a)
                                cor_pos1= getLength(sam_info0->cigar);

                                if(ctxt_1->pe_pos.a)
                                cor_pos2= getLength(cigar);

				if(ctxt_0->pe_pos.a){
                                        flag[0]=ctxt_0->flag|0x1|0x10|0x40;
                                        flag[1]=ctxt_1->flag|0x20|0x80;
                                        }
                                else    {
                                        flag[0]=ctxt_0->flag|0x1|0x40;
                                        flag[1]=ctxt_1->flag|0x1|0x80;
                                        }
                                if(ctxt_1->pe_pos.a){
                                        flag[0]=flag[0]|0x20|0x1|0x40;
                                        flag[1]=flag[1]|0x1|0x80|0x10;
                                        }
                                else    {
                                        flag[0]=flag[0]|0x1|0x40;
                                        flag[1]=flag[1]|0x1|0x80;
                                        }

				//reads1	
				sam_buf += sprintf(sam_buf,"%s\t%u\t",ctxt_0->seq->name,flag[0]);
                                sam_buf += sprintf(sam_buf,"%s\t%u\t",sam_info0->ref_name,ctxt_0->pe_pos.pos);
                                sam_buf += sprintf(sam_buf,"%u\t%s\t",ctxt_0->mapQ,sam_info0->cigar);
                                sam_buf += sprintf(sam_buf,"=\t%u\t%d\t",pos,pos+cor_pos2-ctxt_0->pe_pos.pos-cor_pos1);
                                if(ctxt_0->pos->a)
                                        sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_0->seq_comp,ctxt_0->qual_comp);
                                else    sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_0->seq->seq,ctxt_0->seq->qual);
				
				if(ctxt_0->c0 >1 )
                                	sam_buf += sprintf(sam_buf,"XT:A:R\t");
                                else	sam_buf += sprintf(sam_buf,"XT:A:U\t");
                                        sam_buf += sprintf(sam_buf,"NM:i:%u\t",ctxt_0->pe_pos.n_mm);
                                	if(sam_info0->ref_amb)
                                        sam_buf += sprintf(sam_buf,"XN:i:%u\t",sam_info0->ref_amb);
					sam_buf += sprintf(sam_buf,"SM:i:%u\tAM:i:%u\t",ctxt_0->seQ,am);

                                        sam_buf += sprintf(sam_buf,"X0:i:%u\tX1:i:%u\t",ctxt_0->c0,ctxt_0->c1-ctxt_0->c0);
                                        sam_buf += sprintf(sam_buf,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_0->pe_pos.n_mm,ctxt_0->pe_pos.n_gapo,ctxt_0->pe_pos.n_gapo+ctxt_0->pe_pos.n_gape);

                                        sam_buf += sprintf(sam_buf,"MD:Z:%s",sam_info0->md);	
					for(i=0;i<sam_info0->num;i++){
                                        	pos_t  *pos1 = ctxt_0->pos+i*sam_info0->step;
						if(sam_info0->repeat){
							if(ctxt_0->pe_pos.pos!=pos1->pos){
								if(j==0)
                                        			sam_buf += sprintf(sam_buf,"\tXA:Z:");
								if(j<3)
                                        			sam_buf += sprintf(sam_buf,"%s,%c%u,%s,%u;",sam_info0->m_ref_name[i],pos1->a?'-':'+',pos1->pos,sam_info0->m_cigar[i],pos1->n_mm);
								j++;
							}
						}
                                        	FREE(sam_info0->m_cigar[i]);
						FREE(sam_info0->m_ref_name[i]);
                               	}
				j=0;
				sam_buf += sprintf(sam_buf,"\n");

			//reads2
				ctxt_1->seq_comp_sz=0;
                                ctxt_1->qual_comp_sz=0;
                                realloc_n(ctxt_1->seq_comp,ctxt_1->seq_comp_sz, ctxt_1->seq->length+1);
                                realloc_n(ctxt_1->qual_comp,ctxt_1->qual_comp_sz, ctxt_1->seq->length+1);
                                seq_to_comp(ctxt_1->seq_comp,ctxt_1->seq->seq,ctxt_1->seq->length);
                                ctxt_1->seq_comp[ctxt_1->seq->length]=0;
                                reverse(ctxt_1->qual_comp,ctxt_1->seq->qual,ctxt_1->seq->length);
                                ctxt_1->qual_comp[ctxt_1->seq->length]=0;
				sam_buf += sprintf(sam_buf,"%s\t%u\t",ctxt_1->seq->name,flag[1]);
                                sam_buf += sprintf(sam_buf,"%s\t%u\t",sam_info0->ref_name,pos);
                                sam_buf += sprintf(sam_buf,"%u\t%s\t",ctxt_1->mapQ,cigar);
                                sam_buf += sprintf(sam_buf,"=\t%u\t%d\t",ctxt_0->pe_pos.pos,ctxt_0->pe_pos.pos+cor_pos1-pos-cor_pos2);
				if(ctxt_1->pe_pos.a)
                                        sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_1->seq_comp,ctxt_1->qual_comp);
                                else    sam_buf += sprintf(sam_buf,"%s\t%s\t",ctxt_1->seq->seq,ctxt_1->seq->qual);
                     
				
                                sam_buf += sprintf(sam_buf,"XT:A:M\t");
                                sam_buf += sprintf(sam_buf,"NM:i:%u\t",ctxt_1->pe_pos.n_mm);
				if(sam_info0->ref_amb)
                                        sam_buf += sprintf(sam_buf,"XN:i:%u\t",sam_info1->ref_amb);
                                sam_buf += sprintf(sam_buf,"NM:i:%u\t",ctxt_1->pe_pos.n_mm);
                                //sam_buf += sprintf(sam_buf,"XO:i:%u\tX1:i:%u\t",ctxt_1->c0,ctxt_1->c1-ctxt_1->c0);
                                sam_buf += sprintf(sam_buf,"XM:i:%u\tXO:i:%u\tXG:i:%u\t",ctxt_1->pe_pos.n_mm,ctxt_1->pe_pos.n_gapo,ctxt_1->pe_pos.n_gapo+ctxt_1->pe_pos.n_gape);
                                sam_buf += sprintf(sam_buf,"MD:Z:%s",ctxt_1->md);
                                sam_buf += sprintf(sam_buf,"\n");

				FREE(ctxt_0->pos);
				FREE(ctxt_1->pos);
                                FREE(ctxt_1->cigar);
                                FREE(ctxt_1->md);
				FREE(ctxt_1->seq_comp);
				FREE(ctxt_1->qual_comp);
				seqdb_release(seqdb,ctxt_0->seq);
				seqdb_release(seqdb,ctxt_1->seq);
                                FREE(ctxt_0->qual_comp);
                                FREE(ctxt_0->seq_comp);
                                FREE(sam_info0->md);
                                FREE(sam_info0->cigar);
                                FREE(sam_info0->ref_name);


				break;
				}
			default:break;
			}
		return 0;
}
void cal_mapQ_paring(context_t *ctxt[2],refdb_t refdb,isize_info_t *ii,isize_info_t *last_ii,const pe_opt_t *opt,int n_seq)
{
	int i,j,z;
	uint64_t L;
	pe_data_t *pe_data;
	pe_data = (pe_data_t *) calloc(1,sizeof(pe_data_t));
	L = refdb_ref_len(refdb);
	for( i =0;i < n_seq; i++){
		context_t *ctxt_0 = ctxt[0]+i;
		context_t *ctxt_1 = ctxt[1]+i;
		cal_seq_type(ctxt_1);
		cal_seq_type(ctxt_0);
		if(ctxt_0->pos_n != 0)
		ctxt_0->seQ = ctxt_0->mapQ = approx_mapQ(ctxt_0);
		if(ctxt_1->pos_n != 0)
		ctxt_1->seQ = ctxt_1->mapQ = approx_mapQ(ctxt_1);
		if(ctxt_0->pos_n && ctxt_0->c0==1)
			ctxt_0->pe_pos = ctxt_0->pos[0];
		else	if(ctxt_0->pos_n)
			ctxt_0->pe_pos = ctxt_0->pos[(int)drand48()*(ctxt_0->c0 -1)];
			
		if(ctxt_1->pos_n && ctxt_1->c0==1)
			ctxt_1->pe_pos = ctxt_1->pos[0];
		else	if(ctxt_1->pos_n)
			ctxt_1->pe_pos = ctxt_1->pos[(int)drand48()*(ctxt_1->c0 -1)];
	
	}
	infer_isize(ctxt,ii,opt->ap_prior,L,n_seq);
	
	if (ii->avg < 0.0 && last_ii->avg > 0.0) *ii = *last_ii;
	
	if (opt->force_isize) {
		fprintf(stderr, "[%s] discard insert size estimate as user's request.\n", __func__);
		ii->low = ii->high = 0; ii->avg = ii->std = -1.0;
	}
	memset(pe_data,0,sizeof(pe_data));
	for( i = 0; i < n_seq; i++ ){
		pe_data->stack_n=0;
		context_t *ctxt_0 =ctxt[0]+i;
                context_t *ctxt_1 =ctxt[1]+i;
		for(j = 0;j < 2 ;j++){
			for(z =0; z<ctxt[j][i].pos_n;z++){
				pe_data->stack_n++;
				if(pe_data->stack_n>pe_data->stack_sz){
					pe_data->pos_pe=realloc(pe_data->pos_pe,pe_data->stack_n*sizeof(stack_pos_t));
					pe_data->stack_sz=pe_data->stack_n;
					}
				pe_data->pos_pe[pe_data->stack_n-1].reads_id =j;
				pe_data->pos_pe[pe_data->stack_n-1].pos_m =pe_data->stack_n-1;
				pe_data->pos_pe[pe_data->stack_n-1].pos = ctxt[j][i].pos[z].pos;
				pe_data->pos_pe[pe_data->stack_n-1].a = ctxt[j][i].pos[z].a;
				pe_data->pos_pe[pe_data->stack_n-1].n_mm =ctxt[j][i].pos[z].n_mm;
				pe_data->pos_pe[pe_data->stack_n-1].n_gapo = ctxt[j][i].pos[z].n_gapo;
				pe_data->pos_pe[pe_data->stack_n-1].n_gape = ctxt[j][i].pos[z].n_gape;
				pe_data->pos_pe[pe_data->stack_n-1].score =pe_data->pos_pe[pe_data->stack_n-1].n_mm*3+pe_data->pos_pe[pe_data->stack_n-1].n_gapo*11+pe_data->pos_pe[pe_data->stack_n-1].n_gape*4;
			 	}
			}
		if(pe_data->stack_n)
			pe_data->pos = malloc(pe_data->stack_n*sizeof(uint64_t));
		
		for(z=0; z<pe_data->stack_n;z++){
			pe_data->pos[z] =pe_data->pos_pe[z].pos;
			pe_data->pos[z] =(pe_data->pos[z]<<32)|z;
		}
		
		if(ctxt_1->pos_n !=0 &&ctxt_0->pos_n!=0)	
		pairing(ctxt_0,ctxt_1,pe_data,opt,opt->n_multi,ii);
		if(pe_data->stack_n)
			FREE(pe_data->pos);
	}
	FREE(pe_data->pos_pe);
	free(pe_data);
	
}

void paried_se(context_t *ctxt_0,context_t *ctxt_1,sam_info_t *sam_info_0,sam_info_t *sam_info_1,uint64_t L,refdb_t refdb,isize_info_t *ii,pe_opt_t *opt)
{
		cm_paired_sw(ctxt_0,ctxt_1,ii,refdb,opt,L);	
		if(ctxt_0->is_sw){
			ctxt_1->stdaln_rt = casmap_init_aln_rt();
			se2sam_info(ctxt_1,refdb,sam_info_1,opt->n_multi);
		}
		else	if(ctxt_1->is_sw){
			ctxt_0->stdaln_rt = casmap_init_aln_rt();
			se2sam_info(ctxt_0,refdb,sam_info_0,opt->n_multi);
		}else	if(!ctxt_0->pos_n&&!ctxt_1->pos_n){
		}
		else	if(!ctxt_0->pos_n){
			ctxt_1->stdaln_rt = casmap_init_aln_rt();
			se2sam_info(ctxt_1,refdb,sam_info_1,opt->n_multi);
		}
		else	if(!ctxt_1->pos_n){
			ctxt_0->stdaln_rt = casmap_init_aln_rt();
			se2sam_info(ctxt_0,refdb,sam_info_0,opt->n_multi);
		}
		else	{
			ctxt_0->stdaln_rt=casmap_init_aln_rt();
			ctxt_1->stdaln_rt=casmap_init_aln_rt();
			se2sam_info(ctxt_0,refdb,sam_info_0,opt->n_multi);
			se2sam_info(ctxt_1,refdb,sam_info_1,opt->n_multi);
		}	
}


void  paried_se_std(int tid, context_t *ctxt[2],sam_info_t *sam_info1[2],refdb_t refdb,isize_info_t *ii,pe_opt_t *opt,int n_seq)
{
	int i;
	uint64_t L;
	L = refdb_ref_len(refdb);
	for( i =0 ; i < n_seq; i++){
		context_t *ctxt_0 = ctxt[0]+i;
		context_t *ctxt_1 = ctxt[1]+i;
		sam_info_t *sam_info_0 = sam_info1[0]+i;
		sam_info_t *sam_info_1 = sam_info1[1]+i;
		if( opt->n_thread > 1){
		//	if( (i % opt->n_thread )!= tid ) continue ;
#if 1
			pthread_mutex_lock(&g_seq_lock);
			if( ctxt_0->tid < 0){
				int j;
				for( j = i; j <n_seq && j < i + THREAD_BLOCK_SIZE ; j++){
					ctxt[0][j].tid = tid;
				}
				pthread_mutex_unlock(&g_seq_lock);
			}
			else if( ctxt_0->tid != tid){
				pthread_mutex_unlock(&g_seq_lock);
				continue;
			}
			pthread_mutex_unlock(&g_seq_lock);
#endif
	
		}
		ctxt_0->stdaln_rt=casmap_init_aln_rt();
		ctxt_1->stdaln_rt=casmap_init_aln_rt();
		paried_se(ctxt_0,ctxt_1,sam_info_0,sam_info_1,L,refdb,ii,opt);
	}
}

void getpos(seq_id_t uid,posdb_t posdb,context_t *ctxt,seq_id_t *pos_id ,pos_t *pos,uint64_t L)
{
	ctxt->pos_n=0;	
	int     rc;
	do{
		if(*pos_id==SEQ_ID_INVALID){
			rc=posdb_get(posdb,pos);
			if(rc==-1)
				break;
			*pos_id=SEQ_ID(pos->id);
		}
		if(*pos_id==uid){
			if(pos->pos >= L|| (pos->r && pos->pos < ctxt->seq->length)){
				*pos_id =SEQ_ID_INVALID;
				continue;
			}
			ctxt->pos_n++;	
			ctxt->pos=(pos_t *)realloc(ctxt->pos,ctxt->pos_n*sizeof(pos_t));
			ctxt->pos[ctxt->pos_n-1]= *pos;
			if(ctxt->pos[ctxt->pos_n-1].r)
				ctxt->pos[ctxt->pos_n-1].pos-=ctxt->seq->length;
			*pos_id=SEQ_ID_INVALID;
		}
	}while(*pos_id==SEQ_ID_INVALID);

}

context_t *getposs(seqdb_t seqdb,posdb_t posdb,int need,int *_n_seq,seq_id_t *pos_id,pos_t *pos,uint64_t L)
{
	context_t *ctxt ,*p;
	ctxt=(context_t *)calloc(need,sizeof(context_t));
	int i,n_seq=0;
	for(i=0;i<need;i++){
        	p = ctxt + i;
		p->seq = seqdb_get(seqdb,SEQDB_NEXT);
		if( p ->seq ==NULL)
		break;
		getpos(p->seq->uid,posdb,p,pos_id,pos,L);
		p->tid =-1;
		p->max_diff=3;
		n_seq++;
	}
	if(!n_seq){
		FREE(ctxt);
		*_n_seq = 0;
		return NULL;
	}
	*_n_seq = n_seq;
	fprintf(stderr,"%d reads have found pos. \n",*_n_seq);
	return ctxt;
	
}
int double_getposs(context_t *ctxt_0,context_t *ctxt_1,seqdb_t seqdb,posdb_t posdb,int need,int *_n_seq,seq_id_t *pos_id,pos_t *pos,uint64_t L)
{
	context_t  *p;
	static unsigned long long pos_n = 0;
	static unsigned long long seq_n = 0; 
	int i,n_seq = 0;
	t_curr = get_run_time();
	for ( i = 0 ; i < need ; i++){
		p = &ctxt_0[i];
		p->seq = seqdb_get(seqdb,SEQDB_NEXT);
		if(p->seq == NULL)
			break;
		getpos(p->seq->uid,posdb,p,pos_id,pos,L);
		p->tid = -1;
		p->max_diff = 3;
		pos_n += p->pos_n ;
			
		p = ctxt_1 + i;
		p->seq = seqdb_get(seqdb,SEQDB_NEXT);
		getpos(p->seq->uid,posdb,p,pos_id,pos,L);
		p->tid = -1;
		p->max_diff = 3;
		pos_n += p->pos_n ;
		
		n_seq++;
	}
	seq_n += n_seq;
	fprintf(stderr,"[pos2sam] %.2f sec, %llu reads, %llu pos\n",t_curr,seq_n,pos_n);
	*_n_seq = n_seq;
	return 1;
}
int print_sam1(context_t *ctxt[2],sam_info_t *sam_info[2],FILE *fout,int n_seq,refdb_t refdb)
{
	int i;
	for( i = 0 ; i < n_seq ; i++){
		context_t *ctxt_0=ctxt[0]+i;
		context_t *ctxt_1=ctxt[1]+i;
		sam_info_t *sam_info_0 = sam_info[0]+i;
		sam_info_t *sam_info_1 = sam_info[1]+i;
		if(ctxt_0->is_sw){
			fout_sam_info(sam_info_0,sam_info_1,ctxt_0,ctxt_1,4,fout,refdb);
		}
		else	if(ctxt_1->is_sw){
			fout_sam_info(sam_info_0,sam_info_1,ctxt_0,ctxt_1,5,fout,refdb);	
		}
	
		else	if(!ctxt_0->pos_n&&!ctxt_1->pos_n){
			fout_sam_info(sam_info_0,sam_info_1,ctxt_0,ctxt_1,0,fout,refdb);
		}
		else	if(!ctxt_0->pos_n){
			fout_sam_info(sam_info_0,sam_info_1,ctxt_0,ctxt_1,1,fout,refdb);
		}
		else	if(!ctxt_1->pos_n){
			fout_sam_info(sam_info_0,sam_info_1,ctxt_0,ctxt_1,2,fout,refdb);
		}
		else	{
			fout_sam_info(sam_info_0,sam_info_1,ctxt_0,ctxt_1,3,fout,refdb);
		}
	}
	return 0;

}
int print_sam2(context_t *ctxt[2],sam_info_t *sam_info[2],FILE *fout,int n_seq,refdb_t refdb,char *sam_buf)
{
	int i;
	for( i = 0 ; i < n_seq ; i ++ ){
		context_t *ctxt_0 = ctxt[0] + i ;
		context_t *ctxt_1 = ctxt[1] + i ;
		sam_info_t *sam_info_0 = sam_info[0] + i ;
		sam_info_t *sam_info_1 = sam_info[1] + i ;
		if(ctxt_0->is_sw){
			fout_sam_info1(sam_info_0,sam_info_1,ctxt_0,ctxt_1,4,sam_buf,refdb);
		}
		else	if(ctxt_1->is_sw){
			fout_sam_info1(sam_info_0,sam_info_1,ctxt_0,ctxt_1,5,sam_buf,refdb);	
		}
	
		else	if(!ctxt_0->pos_n&&!ctxt_1->pos_n){
			fout_sam_info1(sam_info_0,sam_info_1,ctxt_0,ctxt_1,0,sam_buf,refdb);
		}
		else	if(!ctxt_0->pos_n){
			fout_sam_info1(sam_info_0,sam_info_1,ctxt_0,ctxt_1,1,sam_buf,refdb);
		}
		else	if(!ctxt_1->pos_n){
			fout_sam_info1(sam_info_0,sam_info_1,ctxt_0,ctxt_1,2,sam_buf,refdb);
		}
		else	{
			fout_sam_info1(sam_info_0,sam_info_1,ctxt_0,ctxt_1,3,sam_buf,refdb);
		}
		fprintf(fout,"%s",sam_buf);


	}
	return 0;
}
static void *uworker(void *data)
{
	thread_data_t *d = (thread_data_t *)data;
	paried_se_std(d->tid,d->ctxt,d->sam_info,d->refdb,&d->ii,d->pe_opt,d->n_seq);	
	
	return 0;

}
static void *worker(void *data)
{
	pipe_data_t *d = (pipe_data_t *)data;
	paried_se_std(d->tid,d->ctxt,d->sam_info,d->refdb,&d->ii,d->pe_opt,d->n_seq);
	//fprintf(stderr,"refined have finished..\n");
	return 0;
}
static void *pworker(void *data)
{
	pipe_data_t *d = (pipe_data_t *)data;
	double_getposs(d->ctxt[0],d->ctxt[1],d->seqdb,d->posdb,BUF_SEQ,&d->n_seq,&pos_id_0,&pos_0,d->L);
	//fprintf(stderr,"geting reads have finished..\n" );
	return 0;
}
static void *foutworker(void *data)
{
	output_data_t *d = (output_data_t *) data;
	print_sam2(d->ctxt,d->sam_info,d->fout,d->n_seq,d->refdb,d->sam_buf);
	//fprintf(stderr,"printf sam have finished \n");
	return 0;

}
int  pipe_cm_pos2sampe(refdb_t refdb,posdb_t posdb,seqdb_t seqdb,FILE *fout,pe_opt_t  *pe_opt)
{

	isize_info_t last_ii,ii;
	int n_ref,i,n_seq;
	char *sam_buf;
	context_t *ctxt[2],*p[2],*c[2];
	sam_info_t *sam_info[2],*s[2];
	output_data_t *out_put;
	uint64_t L;
	L = refdb_ref_len(refdb);

	for (i = 1; i != 256; ++i) 
		g_log_n[i] = (int)(4.343 * log(i) + 0.5);
	t_last = t_curr = get_run_time();	
	n_seq = 0;
	last_ii.avg = -1.0;
	n_ref = refdb_max_idx(refdb)+1;
	c[0] = NULL;
	c[1] = NULL;
	s[0] = NULL;
	s[1] = NULL;
	
	for( i = 0 ; i < n_ref ; ++i )
        	fprintf(fout,"@SQ\tSN:%s\tLN:%u\n",refdb_name(refdb,i),refdb_length(refdb,i));
	fprintf(fout,"@PG\tID:casmap\n");
	
	ctxt[0] = (context_t *)calloc(BUF_SEQ,sizeof(context_t));
	ctxt[1] = (context_t *)calloc(BUF_SEQ,sizeof(context_t));
	sam_buf = (char *) calloc(4000,sizeof(char));
	out_put = (output_data_t *)calloc(1,sizeof(output_data_t));
	out_put -> n_seq = 0;
	out_put -> refdb = refdb;
	double_getposs(ctxt[0],ctxt[1],seqdb,posdb,BUF_SEQ,&n_seq,&pos_id_0,&pos_0,L);

	do{
	//
		if(n_seq==0){
			break;
		}
		p[0] =  (context_t *)calloc(BUF_SEQ,sizeof(context_t));
		p[1] =  (context_t *)calloc(BUF_SEQ,sizeof(context_t));
		sam_info[0]=(sam_info_t *)calloc(n_seq,sizeof(sam_info_t));
		sam_info[1]=(sam_info_t *)calloc(n_seq,sizeof(sam_info_t));
		if(pe_opt->n_thread<=1){
			cal_mapQ_paring(ctxt,refdb,&ii,&last_ii,pe_opt,n_seq);
			paried_se_std(0,ctxt,sam_info,refdb,&ii,pe_opt,n_seq);
			print_sam1(ctxt,sam_info,fout,n_seq,refdb);
			double_getposs(p[0],p[1],seqdb,posdb,BUF_SEQ,&n_seq,&pos_id_0,&pos_0,L);
			FREE(ctxt[0]);
			FREE(ctxt[1]);
			ctxt[0] = p[0];
			ctxt[1] = p[1];
			FREE(sam_info[0]);
			FREE(sam_info[1]);
			continue;
		}
		else {
			pthread_t *tid;
			pthread_attr_t attr;
			pipe_data_t *data;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (pipe_data_t *)calloc(pe_opt->n_thread+1 ,sizeof(pipe_data_t));
			tid = (pthread_t *)calloc(pe_opt->n_thread+1,sizeof(pthread_t));	
			cal_mapQ_paring(ctxt,refdb,&ii,&last_ii,pe_opt,n_seq);
			// creat_thread
			
			for( i =  0; i <= pe_opt->n_thread; i++){
				if( i == 0) {
					data[i].ctxt[0] = p[0];
					data[i].ctxt[1] = p[1];
					data[i].seqdb = seqdb;
					data[i].posdb = posdb;
					data[i].L = L;
 					pthread_create(&tid[i], &attr, pworker, data + i);
				}else  if ( i == 1 ){
					pthread_create(&tid[i] ,&attr, foutworker , out_put);
				}else {
					data[i].tid = i;
					data[i].n_seq = n_seq;
					data[i].L = L;
					data[i].ctxt[0] = ctxt[0];
					data[i].ctxt[1] = ctxt[1];
					data[i].sam_info[0] = sam_info[0];
					data[i].sam_info[1] = sam_info[1];
					data[i].pe_opt = pe_opt;
					data[i].ii = ii;
					data[i].refdb =refdb;		
					pthread_create(&tid[i], &attr, worker, data + i );
				}
			}
			for( i = 0 ; i <= pe_opt->n_thread ; i++){
					pthread_join(tid[i], 0);
			}
			// next...do...
			out_put->n_seq = n_seq;
			n_seq = data[0].n_seq;
			free(data);
			free(tid);
		}// how to free data...

		FREE(c[0]);
		FREE(c[1]);
		FREE(s[0]);
		FREE(s[1]);

		c[0] = ctxt[0];
		c[1] = ctxt[1];
		s[0] = sam_info[0];
		s[1] = sam_info[1];

		out_put->ctxt[0] = c[0] ;
		out_put->ctxt[1] = c[1] ;
		out_put->sam_info[0] = s[0] ;
		out_put->sam_info[1] = s[1] ;
		out_put->fout = fout;
		out_put->sam_buf = sam_buf;
		ctxt[0] = p[0];
		ctxt[1] = p[1];
	}while(1);
	if(pe_opt->n_thread > 1)
	print_sam2(out_put->ctxt,out_put->sam_info,fout,out_put->n_seq,refdb,out_put->sam_buf);
	t_curr = get_run_time();
	fprintf(stderr,"[pos2sam] %.2f sec have finished. \n",t_curr);
	FREE(c[0]);
	FREE(c[1]);
	FREE(sam_info[0]);
	FREE(sam_info[1]);
	FREE(ctxt[0]);
	FREE(ctxt[1]);
	free(sam_buf);
	free(out_put);
	return 0;
}
int  cm_pos2sampe(refdb_t refdb,posdb_t posdb[2],seqdb_t seqdb[2],FILE *fout,pe_opt_t *pe_opt){

	isize_info_t last_ii,ii;
	int n_ref,i,n_seq;
	context_t  *ctxt[2];
	sam_info_t *sam_info[2];
	uint64_t L;
	L = refdb_ref_len(refdb);

	for (i = 1; i != 256; ++i) 
		g_log_n[i] = (int)(4.343 * log(i) + 0.5);

	n_seq = 0;
	n_ref = refdb_max_idx(refdb)+1;
	for( i = 0 ; i < n_ref ; ++i )
        	fprintf(fout,"@SQ\tSN:%s\tLN:%u\n",refdb_name(refdb,i),refdb_length(refdb,i));
	fprintf(fout,"@PG\tID:casmap\n");


	while( (ctxt[0] = getposs(seqdb[0],posdb[0],BUF_SEQ,&n_seq,&pos_id_0,&pos_0,L))!=NULL){
		ctxt[1] = getposs(seqdb[1],posdb[1],BUF_SEQ,&n_seq,&pos_id_1,&pos_1,L);
		sam_info[0]=(sam_info_t *)calloc(n_seq,sizeof(sam_info_t));
		sam_info[1]=(sam_info_t *)calloc(n_seq,sizeof(sam_info_t));
		//calculate aln2seq 
		if(pe_opt->n_thread<=1){
			cal_mapQ_paring(ctxt,refdb,&ii,&last_ii,pe_opt,n_seq);
			paried_se_std(0,ctxt,sam_info,refdb,&ii,pe_opt,n_seq);
			print_sam1(ctxt,sam_info,fout,n_seq,refdb);
		}
		else {
			pthread_t *tid;
			pthread_attr_t attr;
			thread_data_t *data;
			pthread_attr_init(&attr);
			pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);
			data = (thread_data_t *)calloc(pe_opt->n_thread,sizeof(thread_data_t));
			tid = (pthread_t *)calloc(pe_opt->n_thread,sizeof(pthread_t));
			cal_mapQ_paring(ctxt,refdb,&ii,&last_ii,pe_opt,n_seq);
			// creat_thread
			for( i = 0; i < pe_opt->n_thread; i++){
				data[i].tid = i;
				data[i].n_seq = n_seq;
				data[i].ctxt[0] = ctxt[0];
				data[i].ctxt[1] = ctxt[1];
				data[i].sam_info[0] = sam_info[0];
				data[i].sam_info[1] = sam_info[1];
				data[i].pe_opt = pe_opt;
				data[i].ii = ii;
				data[i].refdb =refdb;
				pthread_create(&tid[i], &attr, uworker, data + i);
			}
			for( i = 0 ; i < pe_opt->n_thread ; i++){
				pthread_join(tid[i], 0);
			}
			free(data);
			free(tid);
			print_sam1(ctxt,sam_info,fout,n_seq,refdb);
		}
		last_ii =ii;
		FREE(sam_info[0]);
		FREE(sam_info[1]);
		FREE(ctxt[0]);
		FREE(ctxt[1]);
	}

	return 0;
}


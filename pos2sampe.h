#include<stdint.h>
#include"posdb.h"
#include"seqdb.h"
#define SW_MIN_MAPQ 17
#define FROM_M 0
#define FROM_I 1
#define FROM_D 2
#define FROM_S 3

#define	FREE(X) \
	do{\
		if((X)!=NULL){	\
		free((X));\
		(X)=NULL;\
		}\
	}while(0)

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

typedef struct {
	double avg, std, ap_prior;
	uint32_t low, high, high_bayesian;
} isize_info_t;



typedef struct {
	seq_t *seq;
	char *nt4[2];
	int	tid;
	unsigned int nt4_sz[2];
	char *seq_comp;
	unsigned int seq_comp_sz;
	char *qual_comp;
	unsigned int qual_comp_sz;
	int valid[2];
	char *ref_nt4;
	unsigned int ref_nt4_sz;;
	unsigned int ref_len;

	uint32_t	*cigar;
	int		n_cigar;
	char		*md;
	int		n_md;
	int		sw;
	pos_t *pos;
	pos_t pe_pos;
	unsigned int flag;
	unsigned int is_sw;
	unsigned int pos_n;
	unsigned int pos_sz;
	unsigned int mapQ;
	unsigned int seQ;
	unsigned int max_diff;
	unsigned int c0;
	unsigned int c1;
	stdaln_rt *stdaln_rt;
}context_t;



typedef struct {
	unsigned int reads_id;
	unsigned int a;
	unsigned int pos_m;
	uint8_t n_mm;
	uint8_t n_gapo;
	uint8_t n_gape;
	seq_sz_t pos;
	int score;
}stack_pos_t;

typedef struct{
	stack_pos_t *pos_pe;
	uint64_t      *pos;
	unsigned int  stack_n;
	unsigned int  stack_sz;
}pe_data_t;


typedef struct{
	
	int    	        is_sw;
	int	       repeat;
	// unique
	char  	       *cigar;
	char	    *ref_name;
	int	     ref_amb;
	char	          *md;
	int		 step;
	int		  num;
	//multiple
	char	     *m_cigar[10];
	char         *m_ref_name[10];
}sam_info_t;

typedef	struct {
	int max_isize,force_isize;
	int max_occ;
	int n_multi,N_multi;
	int n_thread;
	int type,is_sw,is_preload;
	double ap_prior;
}pe_opt_t;

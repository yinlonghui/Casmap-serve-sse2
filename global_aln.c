#include <stdio.h>
#include <fcntl.h>
#include <stdlib.h>
#include"refdb.h"
#include "utils.h"
#include"seq.h"
#include"utils.h"
#include"seqdb.h"
#include"ksw.h"
#include <math.h>
#include "global_aln.h"

#define edit_char(m) (((m)==FROM_M)?'M':((m)==FROM_I)?'I':((m)==FROM_D)?'D':'S')

extern void bwa_fill_scmat(int a, int b, int8_t mat[25]);


extern int  sw2cigar(char *cigar_,bwa_cigar_t *cigar,int n_cigar);

extern int	cm_cal_md(int *n_cigar,bwa_cigar_t *cigar,char *ref_nt4, char *nt4, int ref_len,char *md, uint8_t *nm);

/*
 *  Replace the global alignment with ksw global function.
*/

int cm_aln_std(char* ref_nt4,int ref_len,char* seq_nt4, int seq_len, char *cigar, char *md, int *nm, stdaln_rt *rt)
{
	int offset=0;
	int n_cigar = 0 ;
	bwa_cigar_t  *tmp= 0 ;
	int  i ;


	uint32_t   *cigar32 = 0 ;
	int8_t  mat[25];
	bwa_fill_scmat(1,3,mat);

	/* SSE2 global align */
	ksw_global(seq_len , (uint8_t *)seq_nt4 , ref_len , (uint8_t *) ref_nt4 , 5 , mat , 5 , 1 , 50 > (ref_len - seq_len )? 50 :  1.5*(ref_len - seq_len) , &n_cigar , &cigar32);
	
	if ((cigar32[n_cigar - 1]&0xf) == 1) cigar32[n_cigar - 1] = (cigar32[n_cigar - 1]>>4<<4) | 3; // change endding ins to soft clipping
	if ((cigar32[0]&0xf) == 1) cigar32[0] = (cigar32[0]>>4<<4) | 3; // change beginning ins to soft clipping
	if ((cigar32[n_cigar - 1]&0xf) == 2) --n_cigar; // delete endding del
	if ((cigar32[0]&0xf) == 2) { // delete beginning del
		offset += cigar32[0]>>4;
		--n_cigar;
		memmove(cigar32, cigar32+1, (n_cigar) * 4);
	}
	//  generate MD  and CIGAR string
	tmp =  cigar32 ; 
	for(  i  = 0 ; i < n_cigar ; i++)
		tmp[i] = __cigar_create((cigar32[i]&0xf),cigar32[i]>>4); 

	sw2cigar(cigar , tmp , n_cigar ) ;
	cm_cal_md(&n_cigar, tmp , ref_nt4, seq_nt4 , ref_len , md ,(uint8_t *)nm);
	
	rt->path_sz = 0;
	free(tmp);
	return offset;
}

int cm_aln_simple(char* ref_nt4,int ref_len,char* seq_nt4,
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


uint32_t *aln_path2cigar32(const path_t *path, int path_len, int *n_cigar)
{
	int i, n;
	uint32_t *cigar;
	unsigned char last_type;

	if (path_len == 0 || path == 0) {
		*n_cigar = 0;
		return 0;
	}

	last_type = path->ctype;
	for (i = n = 1; i < path_len; ++i) {
		if (last_type != path[i].ctype) ++n;
		last_type = path[i].ctype;
	}
	*n_cigar = n;
	cigar = (uint32_t*)malloc(*n_cigar * 4);

	cigar[0] = 1u << 4 | path[path_len-1].ctype;
	last_type = path[path_len-1].ctype;
	for (i = path_len - 2, n = 0; i >= 0; --i) {
		if (path[i].ctype == last_type) cigar[n] += 1u << 4;
		else {
			cigar[++n] = 1u << 4 | path[i].ctype;
			last_type = path[i].ctype;
		}
	}

	return cigar;
}
#define edit_char(m) (((m)==FROM_M)?'M':((m)==FROM_I)?'I':((m)==FROM_D)?'D':'S')

bwa_cigar_t *bwa_aln_path2cigar(const path_t *path, int path_len, int *n_cigar)
{
	uint32_t *cigar32;
	bwa_cigar_t *cigar;
	int i;
	cigar32 = aln_path2cigar32((path_t*) path, path_len, n_cigar);
	cigar = (bwa_cigar_t*)cigar32;
	for (i = 0; i < *n_cigar; ++i)
                cigar[i] = __cigar_create( (cigar32[i]&0xf), (cigar32[i]>>4) );
	return cigar;
}
		
int se2sam_info(context_t *ctxt,refdb_t refdb,sam_info_t *sam_info1,int n_multi)
{
		
		sam_info1->cigar =NULL;
		char *aln_nt4=NULL;
		unsigned long long last_id=-1ll;
		char cigar[4096];
		char md[4096];
		uint32_t ref_pos;
		int ref_ext;
		pos_t *ppos;
		int i,j;
		int rc;
		const char *ref_name=NULL;
		int ref_offset=0;
		unsigned int ref_coor = 0;
		unsigned int ref_id;
		unsigned int ref_amb=0;
		int rest,step,repeat;
		int	cnt;
                int  n_cigar,n_md,n_ref_name;
		int multi=n_multi+1;
		int nm = 0;
		i=0,j=0,cnt=0;
		
	
		realloc_n(ctxt->nt4[0],ctxt->nt4_sz[0],ctxt->seq->length);
		seq_to_nt4(ctxt->nt4[0],ctxt->seq->seq,ctxt->seq->length);
		realloc_n(ctxt->seq_comp,ctxt->seq_comp_sz, ctxt->seq->length+1);
		realloc_n(ctxt->qual_comp,ctxt->qual_comp_sz, ctxt->seq->length+1);
		realloc_n(ctxt->nt4[1],ctxt->nt4_sz[1],ctxt->seq->length);
		seq_to_comp(ctxt->seq_comp,ctxt->seq->seq,ctxt->seq->length);
		ctxt->seq_comp[ctxt->seq->length]=0;
		reverse(ctxt->qual_comp,ctxt->seq->qual,ctxt->seq->length);
		ctxt->qual_comp[ctxt->seq->length]=0;
		seq_to_nt4(ctxt->nt4[1],ctxt->seq_comp,ctxt->seq->length);
		
		rest=ctxt->c1>multi?multi:ctxt->c1;
		step=ctxt->c1/multi;
		if(step==0) step=1;
		step =step;
		repeat=ctxt->c1>1?1:0;
		
		sam_info1->step=step;
		sam_info1->repeat=repeat;

		while(rest>0){
			if(i==0){
					ppos=&ctxt->pe_pos;
					ppos->a =ctxt->pe_pos.a;
				}
			else	ppos=ctxt->pos+cnt;

		
			if(ppos->id!=last_id){
				ref_pos=ppos->pos;
				ref_ext=ppos->n_gapo+ppos->n_gape;
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
					aln_nt4=ctxt->nt4[1];
				}else{
					aln_nt4=ctxt->nt4[0];
				}
				if(ref_ext){
					ref_pos+=cm_aln_std(ctxt->ref_nt4,ctxt->ref_len,
					aln_nt4,ctxt->seq->length,
					cigar,md,&nm,ctxt->stdaln_rt);
					}
				else if(i==0)
					cm_aln_simple(ctxt->ref_nt4,ctxt->ref_len,aln_nt4,ctxt->seq->length,cigar,md,&nm);
				else{
					sprintf(cigar,"%dM",ctxt->ref_len);
					nm=ppos->n_mm;
				}
				ref_id=refdb_idx(refdb,ref_pos);
				ref_name=refdb_name(refdb,ref_id);
				ref_coor = refdb_offset(refdb,ref_id);	
				ref_amb=refdb_seq_amb(refdb,ref_pos,ctxt->seq->length);
				ref_offset=ref_pos+1-ppos->pos;
			}
			ppos->pos=ppos->pos+ref_offset-ref_coor;
			ppos->n_mm=nm;

			if(i==0){
				n_cigar = strlen(cigar);
				if(n_cigar){
					sam_info1->cigar = calloc(n_cigar+1,sizeof(char));
					memcpy(sam_info1->cigar,cigar,n_cigar*sizeof(char));
					sam_info1->cigar[n_cigar] = '\0';
				}
			///
				n_ref_name = strlen(ref_name);
				sam_info1->ref_name =calloc(n_ref_name+1,sizeof(char));
				memcpy(sam_info1->ref_name,ref_name,n_ref_name*sizeof(char));
				sam_info1->ref_name[n_ref_name] = '\0';
			//
				n_md =strlen(md);
				sam_info1->md = calloc(n_md+1,sizeof(char));
				memcpy(sam_info1->md,md,n_md*sizeof(char));
				sam_info1->md[n_md] = '\0';
			//
				sam_info1->ref_amb=ref_amb;
			}
			else if(j!=10){
				n_cigar =strlen(cigar);
				sam_info1->m_cigar[j] = calloc(n_cigar+1,sizeof(char));
				memcpy(sam_info1->m_cigar[j],cigar,n_cigar*sizeof(char));
				sam_info1->m_cigar[j][n_cigar] ='\0';
				
				n_ref_name = strlen(ref_name);
				sam_info1->m_ref_name[j] =calloc(n_ref_name+1,sizeof(char));
				memcpy(sam_info1->m_ref_name[j],ref_name,n_ref_name*sizeof(char));
				sam_info1->m_ref_name[j][n_ref_name] ='\0';			
				j++;
			
			}
			if(i!=0)
			cnt+=step;
			i=1;
			if(cnt >= rest -1)
			break;
		}
		sam_info1->num =j;
		ctxt->ref_nt4_sz =0;
		FREE(ctxt->nt4[0]);
		FREE(ctxt->nt4[1]);
		FREE(ctxt->ref_nt4);	
		FREE(ctxt->stdaln_rt);
		return 0;			
}
	



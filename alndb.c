#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>

#include "alndb.h"
#include "utils.h"

#define ALN_BUF_SIZE 1024

typedef struct{
	char *aln_prefix;
	int *aln_fd;
	int aln_fd_sz;
	aln_t **aln_buf;
	unsigned int *aln_offset;
	unsigned long long *aln_n;
	seq_id_t last_uid;
	int last_sid;
	unsigned long sa_intv;
	int flags;
	char *meta;
	size_t meta_sz;
}alndb;

alndb_t alndb_open(const char* prefix, int flags)
{
	alndb *pdb;

	pdb=malloc(sizeof(alndb));
	memset(pdb,0,sizeof(alndb));

	pdb->flags|=flags&0xffff0000;
	if(prefix){	
		pdb->sa_intv=1lu<<(flags&0x3f);
		pdb->aln_prefix=strdup(prefix);
	}else{
		pdb->sa_intv=1lu<<(sizeof(pdb->sa_intv)*8-1);
		pdb->aln_prefix=NULL;
	}

	pdb->last_uid=SEQ_ID_INVALID;
	pdb->last_sid=-1;
	pdb->meta_sz=1;
	pdb->meta=malloc(pdb->meta_sz);
	pdb->meta[0]=0;

	return pdb;
}

void alndb_meta_append(alndb_t db, const char *str)
{
	alndb *pdb=db;
	pdb->meta_sz+=strlen(str);
	pdb->meta=realloc(pdb->meta,pdb->meta_sz);
	strcat(pdb->meta,str);
	strcat(pdb->meta,";");
}

const char* alndb_meta_info(alndb_t db)
{
	alndb *pdb=db;
	return pdb->meta;
}

int alndb_close(alndb_t db)
{
	int i;
	int rc=0;
	FILE *f_meta=NULL;
	alndb *pdb=db;

	if(pdb->aln_prefix)
		f_meta=fopen(pdb->aln_prefix,"w");

	if(f_meta)
		fprintf(f_meta,"#%s\n",pdb->meta);

	for(i=0;i<pdb->aln_fd_sz;++i){
		if(pdb->aln_buf[i]){
			rc=xwrite(pdb->aln_fd[i],pdb->aln_buf[i],sizeof(aln_t)*pdb->aln_offset[i]);
			if(rc==-1){
				xerror("Dump aln buffer");
				rc=-1;
				goto ERR_0;
			}
			close(pdb->aln_fd[i]);
			free(pdb->aln_buf[i]);
			if(f_meta)
				fprintf(f_meta,"%u\t%llu\n",i+1,pdb->aln_n[i]);
		}
	}
	if(f_meta)
		fclose(f_meta);

	free(pdb->aln_fd);
	free(pdb->aln_buf);
	free(pdb->aln_offset);
	free(pdb->aln_n);
	free(pdb->aln_prefix);
	free(pdb->meta);

	free(pdb);
	return 0;
ERR_0:
	return rc;
}

static inline int check_aln_buffer(alndb *db, unsigned int i)
{
	alndb *pdb=db;
	unsigned int prev_sz=pdb->aln_fd_sz;
	unsigned int new_sz=i+1;
	char *aln_fn;
	if(new_sz>prev_sz){
		pdb->aln_fd=realloc(pdb->aln_fd,sizeof(int)*new_sz);
		memset(pdb->aln_fd+prev_sz,0,sizeof(int)*(new_sz-prev_sz));
		pdb->aln_buf=realloc(pdb->aln_buf,sizeof(void*)*new_sz);
		memset(pdb->aln_buf+prev_sz,0,sizeof(void*)*(new_sz-prev_sz));
		pdb->aln_offset=realloc(pdb->aln_offset,sizeof(unsigned int)*new_sz);
		memset(pdb->aln_offset+prev_sz,0,sizeof(unsigned int)*(new_sz-prev_sz));
		pdb->aln_n=realloc(pdb->aln_n,sizeof(unsigned long long)*new_sz);
		memset(pdb->aln_n+prev_sz,0,sizeof(unsigned long long)*(new_sz-prev_sz));
		pdb->aln_fd_sz=new_sz;
	}
	if(pdb->aln_buf[i]==NULL){
		if(pdb->aln_prefix){
			aln_fn=malloc(strlen(pdb->aln_prefix)+32);
			sprintf(aln_fn,"%s.%d",pdb->aln_prefix,i+1);
			pdb->aln_fd[i]=open(aln_fn,O_WRONLY|O_TRUNC|O_CREAT,S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH);
			free(aln_fn);
		}else{
			pdb->aln_fd[i]=STDOUT_FILENO;
		}
		if(pdb->aln_fd[i]==-1){
			//xerror("Open aln part file");
			xerror(aln_fn);
			pdb->flags|=ALNDB_ERROR;
			return -1;
		}
		pdb->aln_buf[i]=malloc(sizeof(aln_t)*ALN_BUF_SIZE);
		memset(pdb->aln_buf[i],0,sizeof(aln_t)*ALN_BUF_SIZE); // memory check may report error without initialize
	}
	return 0;
}

int alndb_put(alndb_t db, hit_t *hit, seq_id_t uid)
{
	int n;
	int cont;
	unsigned long sa, width;
	aln_t aln_buffer;
	aln_t *paln=&aln_buffer;
	alndb *pdb=db;
	unsigned long intv=pdb->sa_intv;

	if(pdb->flags&ALNDB_ERROR)
		return -1;

	if(uid==pdb->last_uid){
		pdb->last_sid++;
	}else{
		pdb->last_uid=uid;
		pdb->last_sid=0;
	}
	
	if(0){
		seq_sz_t width=hit->l-hit->k;
		if(width > 8192){
			fprintf(stderr,
					"WARNING: Enormously large SA interval %u (%u,%u) for "
					"uid=%lu . "
					"An error might happened.\n",
					width,hit->k,hit->l,uid
			       );
		}
	}

	sa=hit->k;
	do{
		width=intv-(sa&(intv-1));

		if(sa+width>hit->l){
			cont=0;
			width=hit->l-sa+1;
		}else{
			sa=(sa&(~(intv-1)))+intv;
			cont=1;
		}

		n=sa/intv;

		if(check_aln_buffer(pdb,n))
			return -1;

		paln=pdb->aln_buf[n]+pdb->aln_offset[n];

		paln->id=MAKE_ALN_ID(uid,pdb->last_sid);
		paln->a=hit->a;
		paln->r=hit->r;
		paln->n_mm=hit->n_mm;
		paln->n_gapo=hit->n_gapo;
		paln->n_gape=hit->n_gape;
		paln->sa=sa;
		paln->width=width;
		pdb->aln_n[n]++;
		pdb->aln_offset[n]++;
		//fprintf(stderr,"ALN: %llu %llu %d\n", seq->uid, paln->id>>16, i);
		if(pdb->aln_offset[n]==ALN_BUF_SIZE){
			if(xwrite(pdb->aln_fd[n],pdb->aln_buf[n],sizeof(aln_t)*ALN_BUF_SIZE)
				!=sizeof(aln_t)*ALN_BUF_SIZE){
				xerror("Write ALN");
				return -1;
			}
			pdb->aln_offset[n]=0;
		}
	}while(cont);
	return 0;
}

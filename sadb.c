#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "sadb.h"
#include "utils.h"

typedef struct{
	int sa_fd[2];
	unsigned long sa_intv;
	seq_sz_t* sa_data[2];
	unsigned int sa_idx[2];
	seq_sz_t ref_len;
	size_t mm_sz;
	int flags;
}sadb;

static int sadb_prefetch(sadb *pdb, int s, unsigned int idx)
{
	int rc;
	off_t s_off;
	s_off=(off_t)4*(pdb->sa_intv)*idx;
	//s_off=(off_t)(4)*((sa&(~(pdb->sa_intv-1))));
	
	if(pdb->flags&SADB_IN_MEM){
		//double t_start=get_run_time();
		//fprintf(stderr,"[sadb] Prefetch %d %lu %lu...\n",s,s_off,pdb->mm_sz);
		rc=xpread(pdb->sa_fd[s],pdb->sa_data[s],
				pdb->mm_sz, s_off);
		if(rc==-1){
			xerror("pread");
			return -1;
		}
		//fprintf(stderr,"[sadb] Prefetched %d %lu %lu in %.2f sec\n",s,s_off,pdb->mm_sz,get_run_time()-t_start);
	}else{
		if(pdb->sa_idx[s]!=-1)
			munmap(pdb->sa_data[s],4*pdb->sa_intv);
		pdb->sa_data[s]=mmap(NULL,4*pdb->sa_intv,
				PROT_READ,MAP_PRIVATE,pdb->sa_fd[s],s_off);
		if(pdb->sa_data[s]==MAP_FAILED){
			xerror("mmap");
			return -1;
		}
	}
	pdb->sa_idx[s]=idx;
	return 0;
}

sadb_t sadb_open(const char* prefix, int flags)
{
	int i;
	sadb *pdb;
	char *fn_buf;
	if(prefix==NULL || strlen(prefix)==0){
		return NULL;
	}

	pdb=malloc(sizeof(sadb));
	memset(pdb,0,sizeof(sadb));

	pdb->flags=flags;
	pdb->sa_intv=1lu<<(flags&0x3f);

	fn_buf=malloc(strlen(prefix)+8);
	sprintf(fn_buf,"%s.sa",prefix);
	pdb->sa_fd[0]=open(fn_buf,O_RDONLY);
	if(pdb->sa_fd[0]==-1){
		xerror("Open .sa file");
		goto ERR_0;
	}else{
		if(read(pdb->sa_fd[0],&pdb->ref_len,4)!=-1){
			close(pdb->sa_fd[0]);
			pdb->sa_fd[0]=open(fn_buf,O_RDONLY);
		}else{
			xerror("Read SA length");
			goto ERR_0;
		}
	}
	sprintf(fn_buf,"%s.rsa",prefix);
	pdb->sa_fd[1]=open(fn_buf,O_RDONLY);
	if(pdb->sa_fd[1]==-1){
		xerror("Open .rsa file");
		goto ERR_0;
	}
	while(pdb->sa_intv/2>=pdb->ref_len)
		pdb->sa_intv/=2;
	pdb->mm_sz=pdb->sa_intv*4;

	if(pdb->flags&SADB_IN_MEM){
		for(i=0;i<2;++i){

			pdb->sa_data[i]=malloc(pdb->mm_sz);
			if(pdb->sa_data[i]==NULL){
				xerror("malloc");
				goto ERR_0;
			}
			pdb->sa_idx[i]=-1;
		}
	}

	free(fn_buf);
	return pdb;
ERR_0:
	free(fn_buf);
	close(pdb->sa_fd[0]);
	close(pdb->sa_fd[1]);
	free(pdb);
	return NULL;
}

int sadb_close(sadb_t db)
{
	int i;
	sadb *pdb=db;
	if(pdb->flags&SADB_IN_MEM){
		for(i=0;i<2;++i)
			free(pdb->sa_data[i]);
	}else{
		for(i=0;i<2;++i)
			munmap(pdb->sa_data[i],4*pdb->sa_intv);
	}
	for(i=0;i<2;++i){
		close(pdb->sa_fd[i]);
	}
	free(pdb);
	return 0;
}

seq_sz_t sadb_get(sadb_t db, int s, seq_sz_t sa)
{
	sadb *pdb=db;
	seq_sz_t linear;
	unsigned int idx,offset;

	//sa+=6;// SA file header
	idx=sa/pdb->sa_intv;
	offset=sa%pdb->sa_intv;
	//offset=sa&(pdb->sa_intv-1);

	if(idx!=pdb->sa_idx[s]){
		if(sadb_prefetch(pdb,s,idx))
			return SEQ_SZ_MAX;
	}

	linear=*(pdb->sa_data[s]+offset);
	if(s){
		linear=pdb->ref_len-linear;
	}
	/*
	fprintf(stderr,"%u %u %llu %u %u\n",sa,idx,s_off,offset,linear);
	*/
	return linear;
}

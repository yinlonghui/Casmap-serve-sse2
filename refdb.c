#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <sys/mman.h>

#include "refdb.h"
#include "utils.h"

typedef struct {
	char *name;
	seq_sz_t offset;
	seq_sz_t length;
	seq_sz_t n_ambs;
}ann_t;

typedef struct {
	seq_sz_t offset;
	seq_sz_t length;
	char	c_amb;
}amb_t;

typedef struct{
	seq_sz_t ref_len;
	char *ref_pac;
	off_t ref_pac_sz;

	ann_t *anns;
	unsigned int anns_sz;
	unsigned int seed;

	amb_t *ambs;
	unsigned int ambs_sz;
	
	int flags;
}refdb;

static int init_pac(refdb *pdb, const char* fn, int in_mem)
{
	int fid,rc;
	struct stat stat;
	fid=open(fn,O_RDONLY);
	if(fid==-1){
		xerror("open");
		goto ERR_0;
	}
	if(fstat(fid,&stat)==-1){
		xerror("fstat");
		goto ERR_1;
	}
	if(in_mem){
		pdb->ref_pac=malloc(stat.st_size);
		if(pdb->ref_pac==NULL){
			xerror("Allocate ref memory");
			goto ERR_1;
		}
		pdb->ref_pac_sz=stat.st_size;
		rc=read(fid,pdb->ref_pac,stat.st_size);
		if(rc==-1){
			xerror("Read PAC");
			free(pdb->ref_pac);
			pdb->ref_pac=NULL;
			pdb->ref_pac_sz=0;
			goto ERR_1;
		}
	}else{
		pdb->ref_pac_sz=stat.st_size;
		pdb->ref_pac=(char*)mmap(NULL,pdb->ref_pac_sz,PROT_READ,MAP_PRIVATE,fid,0l);
		if(pdb->ref_pac==MAP_FAILED)
		{
			pdb->ref_pac=NULL;
			pdb->ref_pac_sz=0;
			xerror("mmap");
			goto ERR_1;
		}
	}
	close(fid);
	return 0;

ERR_1:
	close(fid);
ERR_0:
	return -1;
}

static int init_ann(refdb *pdb, const char* fn)
{
	FILE *fd;
	seq_sz_t length;
	unsigned int n_seq;
	unsigned int tmp;
	unsigned int i;
	char *name=malloc(1024);
	char *wa=malloc(1024);

	fd=fopen(fn,"r");
	if(fd==NULL){
		xerror("Open ann");
		goto ERR_0;
	}

	if(fscanf(fd,"%u%u%u",&length,&n_seq,&pdb->seed)!=3){
		fprintf(stderr,"ann file format error\n");
		goto ERR_1;
	}

	pdb->anns=calloc(n_seq,sizeof(ann_t));

	for(i=0;i<n_seq;++i){
		if(fscanf(fd,"%u%s%s",&tmp,name,wa)!=3){
			fprintf(stderr,"ann file format error\n");
			goto ERR_2;
		}
		pdb->anns[i].name=strdup(name);
		if(fscanf(fd,"%u%u%u",&(pdb->anns[i].offset),
					&(pdb->anns[i].length),&(pdb->anns[i].n_ambs))!=3){
			fprintf(stderr,"ann file format error\n");
			goto ERR_3;
		}
	}
	pdb->anns_sz=i;
	fclose(fd);
	free(name);
	free(wa);
	return 0;
ERR_3:
	for(i=0;i<n_seq;++i){
		free(pdb->anns[i].name);
	}
ERR_2:
	free(pdb->anns);
	pdb->anns=NULL;
ERR_1:
	fclose(fd);
ERR_0:
	free(name);
	free(wa);
	return -1;
}

int init_amb(refdb *pdb, const char* fn)
{
	int i;
	FILE *fd;
//	seq_sz_t ref_len;
	unsigned int n_seq;
	unsigned int n_holes;
	char c_buf[16];

	fd=fopen(fn,"r");
	if(fd==NULL){
		xerror("Open amb file");
		goto ERR_0;
	}

	if(fscanf(fd,"%u%u%u",&pdb->ref_len,&n_seq,&n_holes)!=3){
		fprintf(stderr,"amb file format error\n");
		goto ERR_1;
	}

	pdb->ambs=calloc(n_holes,sizeof(amb_t));
	for(i=0;i<n_holes;++i){
		if(fscanf(fd,"%u%u%s",&pdb->ambs[i].offset,
					&pdb->ambs[i].length,c_buf)!=3){
			xerror("amb file format error");
			goto ERR_2;
		}
		pdb->ambs[i].c_amb=c_buf[0];
	}
	pdb->ambs_sz=i;
	fclose(fd);
	return 0;
ERR_2:
	free(pdb->ambs);
ERR_1:
	fclose(fd);
ERR_0:
	return -1;
}

static void destroy_pac(refdb *pdb)
{
	if(pdb->flags&REFDB_IN_MEM)
		free(pdb->ref_pac);
	else
		munmap(pdb->ref_pac, pdb->ref_pac_sz);
}

static void destroy_ann(refdb *pdb)
{
	int i;
	for(i=0;i<pdb->anns_sz;++i){
		free(pdb->anns[i].name);
	}
	free(pdb->anns);
}

static void destroy_amb(refdb *pdb)
{
	free(pdb->ambs);
}

refdb_t refdb_open(const char* prefix, int flags)
{
	refdb *pdb;
	char *fn_buf;
	if(prefix==NULL || strlen(prefix)==0){
		return NULL;
	}

	pdb=malloc(sizeof(refdb));
	memset(pdb,0,sizeof(refdb));

	pdb->flags|=flags&0xffff0000;

	fn_buf=malloc(strlen(prefix)+8);

	sprintf(fn_buf,"%s.pac",prefix);
	if(init_pac(pdb,fn_buf,pdb->flags&REFDB_IN_MEM)==-1){
		goto ERR_0;
	}

	sprintf(fn_buf,"%s.ann",prefix);
	if(init_ann(pdb,fn_buf)==-1){
		goto ERR_1;
	}

	sprintf(fn_buf,"%s.amb",prefix);
	if(init_amb(pdb,fn_buf)==-1){
		goto ERR_2;
	}

	free(fn_buf);
	return pdb;
ERR_2:
	destroy_ann(pdb);
ERR_1:
	destroy_pac(pdb);
ERR_0:
	free(fn_buf);
	free(pdb);
	return NULL;
}

int refdb_close(refdb_t db)
{
	refdb *pdb=db;
	if(pdb==NULL)
		return -1;
	destroy_pac(pdb);
	destroy_ann(pdb);
	destroy_amb(pdb);
	free(pdb);
	return 0;
}

int refdb_seq_nt4(refdb_t db, char* buf, seq_sz_t pos, seq_sz_t len)
{
	int j=0;
	refdb *pdb=db;
	const char *pac_s=pdb->ref_pac+(pos>>2);
	const char *pac_e=pdb->ref_pac+pdb->ref_pac_sz;

	while(j<len && pac_s<pac_e){
		do{
			buf[j++]=(*pac_s>>((~pos&3)<<1))&3;
			pos++;
		}while(j<len && pos&3);
		pac_s++;
	}
	while(j<len){
		buf[j++]=4;
	}
	return 0;
}

unsigned int refdb_idx(refdb_t db, seq_sz_t pos)
{
	int left,right,mid;
	refdb *pdb=db;
	left=0;
	right=pdb->anns_sz;
	while(left<right){
		mid=(left+right)>>1;
		if(pos>=pdb->anns[mid].offset){
			if(mid==pdb->anns_sz-1) 
				break;
			if(pos<pdb->anns[mid+1].offset)
				break;
			left=mid+1;
		}else{
			right=mid;
		}
	}
	return (left<right)?mid:-1;
}

unsigned int refdb_max_idx(refdb_t db)
{
	refdb *pdb=db;
	return pdb->anns_sz-1;
}

const char* refdb_name(refdb_t db, unsigned int idx)
{
	refdb *pdb=db;
	return pdb->anns[idx].name;
}

seq_sz_t refdb_offset(refdb_t db, unsigned int idx)
{
	refdb *pdb=db;
	return pdb->anns[idx].offset;
}

seq_sz_t refdb_length(refdb_t db, unsigned int idx)
{
	refdb *pdb=db;
	return pdb->anns[idx].length;
}
seq_sz_t refdb_ref_len(refdb_t db)
{
	refdb *pdb=db;
	return pdb->ref_len;
}

unsigned int refdb_amb(refdb_t db, unsigned int idx)
{
	refdb *pdb=db;
	return pdb->anns[idx].n_ambs;
}

unsigned int refdb_seq_amb(refdb_t db, seq_sz_t pos, seq_sz_t len)
{
	refdb *pdb=db;
	unsigned int left, right, mid, nn;
	left = 0; right = pdb->ambs_sz; nn = 0;
	while (left < right) {
		mid = (left + right) >> 1;
		if (pos >= pdb->ambs[mid].offset + pdb->ambs[mid].length) left = mid + 1;
		else if (pos + len <= pdb->ambs[mid].offset) right = mid;
		else { // overlap
			if (pos >= pdb->ambs[mid].offset) {
				nn += pdb->ambs[mid].offset + pdb->ambs[mid].length < pos + len?
					pdb->ambs[mid].offset + pdb->ambs[mid].length - pos : len;
			} else {
				nn += pdb->ambs[mid].offset + pdb->ambs[mid].length < pos + len?
					pdb->ambs[mid].length : len - (pdb->ambs[mid].offset - pos);
			}
			break;
		}
	}
	return nn;
}

unsigned int refdb_seed(refdb_t db)
{
	refdb *pdb=db;
	return pdb->seed;
}

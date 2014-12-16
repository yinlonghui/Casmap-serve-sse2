#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <fcntl.h>
#include <poll.h>

#include "posdb.h"
#include "utils.h"

#define POS_BUFFER_SIZE 	(1024)

typedef struct {
	int *fd;
	pos_t **data;
	int *d_sz;
	int *d_off;
	int part_n;
	int *idx;
	char *meta;
	size_t meta_sz;
}posdb;

static int prefetch_pos(posdb *pdb, int i)
{
	int rc;
	if(pdb->d_sz[i]==-1)
	       	return 0;
	rc=xread(pdb->fd[i],pdb->data[i],sizeof(pos_t)*POS_BUFFER_SIZE);
	if(rc==-1){
		//xerror("Read pos");
		goto ERR;
	}
	rc/=sizeof(pos_t);
	if(rc==0){
		goto ERR;
	}
	pdb->d_sz[i]=rc;
	pdb->d_off[i]=0;
	//fprintf(stderr,"prefetch %d,%d\n",i,rc);
	//if(rc<POS_BUFFER_SIZE)
	//	fprintf(stderr,"last read on %i\n",i);
	return 0;
ERR:
	//fprintf(stderr,"pos part %i Finished\n",i);
	pdb->data[i]->id=-1ll;
	pdb->d_sz[i]=-1;
	pdb->d_off[i]=0;
	return -1;
}

static int order(posdb *pdb, int i, int j)
{
	long long id_i,id_j;
	if(pdb->d_sz[i]<0){
	       if(pdb->d_sz[j]<0)
		       return 0;
	       else 
		       return 1;
	}else if(pdb->d_sz[j]<0){
		return -1;
	}

	id_i=pdb->data[i][pdb->d_off[i]].id;
	id_j=pdb->data[j][pdb->d_off[j]].id;

	if(id_i>id_j)
		return 1;
	else if(id_i<id_j)
		return -1;
	else
		return 0;
}

static void adjust_order(posdb *pdb, int i)
{
	int tmp;
	int t=(i+pdb->part_n)/2;
	while(t>0){
		if(order(pdb,i,pdb->idx[t])>0){
			tmp=i;
			i=pdb->idx[t];
			pdb->idx[t]=tmp;
		}
		t=t/2;
	}
	pdb->idx[0]=i;
}

posdb_t posdb_open(const char* prefix)
{
	int i,n;
	unsigned long long m;
	FILE *f_meta;
	char *part_fn;
	posdb *pdb;
	pdb=malloc(sizeof(posdb));
	memset(pdb,0,sizeof(posdb));
	pdb->meta_sz=1;
	pdb->meta=malloc(pdb->meta_sz);
	pdb->meta[0]=0;
	if(prefix){
		f_meta=fopen(prefix,"r");
		if(f_meta==NULL){
			xerror("Open pos file");
			return NULL;
		}

		{
			int c;
			char *buf=NULL;
			size_t buf_sz=0;
			for(c=fgetc(f_meta);c=='#';c=fgetc(f_meta)){
				if(getline(&buf,&buf_sz,f_meta)!=-1){
					pdb->meta_sz+=strlen(buf);
					pdb->meta=realloc(pdb->meta,pdb->meta_sz);
					strcat(pdb->meta,buf);
				}
			}
			ungetc(c,f_meta);
			free(buf);
		}


		part_fn=malloc(strlen(prefix)+32);
		i=0;
		pdb->part_n=0;

		while(fscanf(f_meta,"%u%llu",&n,&m)==2){
			pdb->part_n++;
			pdb->fd=realloc(pdb->fd,sizeof(int)*pdb->part_n);
			pdb->data=realloc(pdb->data,sizeof(void*)*pdb->part_n);
			pdb->d_sz=realloc(pdb->d_sz,sizeof(int)*pdb->part_n);
			pdb->d_off=realloc(pdb->d_off,sizeof(int)*pdb->part_n);
			sprintf(part_fn,"%s.%u",prefix,pdb->part_n);
			pdb->fd[i]=open(part_fn,O_RDONLY);
			if(pdb->fd[i]==-1){
				xerror("Open pos part file");
				return NULL;
			}
			pdb->data[i]=malloc(sizeof(pos_t)*POS_BUFFER_SIZE);
			pdb->d_sz[i]=0;
			pdb->d_off[i]=0;

			if(prefetch_pos(pdb,i)){
				close(pdb->fd[i]);
				free(pdb->data[i]);
				pdb->data[i]=NULL;
				pdb->part_n--;
			}else{
				++i;
			}
		}

		fclose(f_meta);
		free(part_fn);
	}else{
		pdb->part_n=1;
		pdb->fd=malloc(sizeof(int));
		pdb->data=malloc(sizeof(void*));
		pdb->d_sz=malloc(sizeof(int));
		pdb->d_off=malloc(sizeof(int));
		pdb->fd[0]=STDIN_FILENO;
		pdb->data[0]=malloc(sizeof(pos_t)*POS_BUFFER_SIZE);
		pdb->d_sz[0]=0;
		pdb->d_off[0]=0;
		prefetch_pos(pdb,0);
	}

	pdb->data=realloc(pdb->data,sizeof(void*)*(pdb->part_n+1));
	pdb->d_off=realloc(pdb->d_off,sizeof(int)*(pdb->part_n+1));
	pdb->d_sz=realloc(pdb->d_sz,sizeof(int)*(pdb->part_n+1));

	pdb->data[pdb->part_n]=malloc(sizeof(pos_t));
	pdb->data[pdb->part_n][0].id=-1ll;
	pdb->d_off[pdb->part_n]=0;
	pdb->d_sz[pdb->part_n]=0;

	pdb->idx=malloc(sizeof(int)*pdb->part_n);
	for(i=0;i<pdb->part_n;++i)
		pdb->idx[i]=pdb->part_n;
	for(i=pdb->part_n-1;i>=0;--i)
		adjust_order(pdb,i);

	return pdb;
}

int posdb_close(posdb_t db)
{
	int i;
	posdb *pdb=db;
	for(i=0;i<pdb->part_n;++i){
		close(pdb->fd[i]);
		free(pdb->data[i]);
	}
	free(pdb->data[pdb->part_n]);

	free(pdb->fd);
	free(pdb->data);
	free(pdb->d_sz);
	free(pdb->d_off);
	free(pdb->idx);
	free(pdb->meta);
	free(pdb);

	return 0;
}

int posdb_get(posdb_t db, pos_t *pos)
{
	posdb *pdb=db;
	int i=pdb->idx[0];
	if(pdb->d_sz[i]==-1){
		return -1;
	}

	memcpy(pos,pdb->data[i]+pdb->d_off[i],sizeof(pos_t));
	//fprintf(stderr,"get %llu\t%u\n",pos->id>>16,pos->pos);
	pdb->d_off[i]++;

	if(pdb->d_off[i]==pdb->d_sz[i])
		prefetch_pos(pdb,i);

	adjust_order(pdb,i);

	return 0;
}

const char* posdb_meta_info(posdb_t db)
{
	posdb *pdb=db;
	return pdb->meta;
}

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <string.h>
#include <stddef.h>
#include <sys/types.h>
#include <sys/stat.h>

#include "seqdb.h"
#include "utils.h"

typedef struct {
	FILE *fin;
	char *buf;
	size_t buf_len;
	int flags;
	off_t offset;
	off_t f_sz;
	int seq_len;
}seqdb;

typedef struct {
	seqdb *seqdb[2];
	int next;
}pairdb;

static void _seqdb_release(seqdb *pdb, seq_t *seq)
{
	if(seq){
		seq->name--; // recover pointer. see _seqdb_get 
		free(seq->name);
		free(seq->seq);
		free(seq->qual);
		free(seq);
	}
}

static inline ssize_t _readline(seqdb *pdb,char **buf)
{
	char *tmp=NULL;
	size_t len=0;
	ssize_t rc=getline(&tmp,&len,pdb->fin);

	if(rc<=0){
		free(tmp);
		*buf=NULL;
		return -1;
	}

	tmp[--rc]=0; /* strip '\n' */
	*buf=tmp;

	return rc;
}	

static seqdb* _seqdb_open(const char *seq_fn, int flags)
{
	seqdb *pdb;
	struct stat _stat;
	char *buf=NULL;
	size_t len=0;
	ssize_t rc;
	pdb=malloc(sizeof(seqdb));
	memset(pdb,0,sizeof(seqdb));
	if(stat(seq_fn,&_stat)){
		xerror("stat");
		goto ERR_0;
	}
	pdb->f_sz=_stat.st_size;

	pdb->fin=fopen(seq_fn, "r");
	if(pdb->fin==NULL){
		xerror("open");
		goto ERR_0;
	}

	/* input test */

	/* test name */
	rc=getline(&buf,&len,pdb->fin);
	if(rc==-1)
		goto ERR_1;
	if(buf[0]!='@' && buf[0]!='>')
		goto ERR_1;

	/* test sequence length */
	rc=getline(&buf,&len,pdb->fin);
	if(rc<=0)
		goto ERR_1;
	buf[--rc]=0;
	pdb->seq_len=rc;

	rewind(pdb->fin);

	free(buf);
	return pdb;
ERR_1:
	fprintf(stderr,"[seqdb] Incorrect file format %s\n",seq_fn);
	free(buf);
	fclose(pdb->fin);
ERR_0:
	free(pdb);
	return NULL;
}

static int _seqdb_close(seqdb *pdb)
{
	if(pdb==NULL)
		return -1;

	fclose(pdb->fin);
	free(pdb->buf);
	free(pdb);
	return 0;
}

static seq_t* _seqdb_get(seqdb *pdb, seq_id_t id)
{
	seq_t *pseq;
	ssize_t rc;
	int c;
	if(pdb==NULL)
		return NULL;

	if(id!=SEQDB_NEXT){
		if(fseeko(pdb->fin,id,SEEK_SET)){
			xerror("seek");
			return NULL;
		}
		pdb->offset=ftello(pdb->fin);
	}

	pseq=(seq_t *)malloc(sizeof(seq_t));
	memset(pseq,0,sizeof(seq_t));
	pseq->uid=pdb->offset;

	/* seek for name */
	rc=_readline(pdb,&pseq->name);
	if(rc<=0)
		goto EXP_1;

	if(pseq->name[0]!='@' && pseq->name[0]!='>'){
		goto EXP_1;
	}
	pseq->name++; // Skip '@' or '>'

	/* seek for sequence */
	rc=_readline(pdb,&pseq->seq);
	if(rc<=0)
		goto EXP_2;

	pseq->length=rc;

	c=fgetc(pdb->fin);
	ungetc(c,pdb->fin);
	if(c=='+'){
		//fgetc(pdb->fin);
		char *tmp=NULL;
		size_t len=0;
		rc=getline(&tmp,&len,pdb->fin);
		free(tmp);
		if(rc<=0)
			goto EXP_3;
		rc=_readline(pdb,&pseq->qual);
		if(rc<=0)
			goto EXP_3;
	}
	pdb->offset=ftello(pdb->fin);

	return pseq;
EXP_3:
	free(pseq->qual);
EXP_2:
	free(pseq->seq);
EXP_1:
	free(pseq->name);
	free(pseq);
	return NULL;
}

static size_t _seqdb_seq_len(seqdb_t db)
{
	seqdb *pdb=db;
	return pdb->seq_len;
}

seqdb_t seqdb_open(const char *seq1_fn, const char *seq2_fn, int flags)
{
	pairdb *ppdb;
	seqdb *pdb1=NULL,*pdb2=NULL;
	if(seq1_fn){ 
		pdb1=_seqdb_open(seq1_fn,flags);
		if(pdb1==NULL)
			return NULL;
	}
	if(seq2_fn){
		pdb2=_seqdb_open(seq2_fn,flags);
		if(pdb2==NULL){
			_seqdb_close(pdb1);
			return NULL;
		}
	}
	ppdb=malloc(sizeof(pairdb));
	memset(ppdb,0,sizeof(pairdb));
	ppdb->seqdb[0]=pdb1;
	ppdb->seqdb[1]=pdb2;
	ppdb->next=0;
	return ppdb;
}

int seqdb_close(seqdb_t db)
{
	pairdb *ppdb=db;
	_seqdb_close(ppdb->seqdb[0]);
	_seqdb_close(ppdb->seqdb[1]);
	free(ppdb);
	return 0;
}

seq_t* seqdb_get(seqdb_t db, seq_id_t id)
{
	pairdb *ppdb=db;
	seq_t *pseq;
	int s;
	if(id==SEQDB_NEXT){
		s=ppdb->next;
		if(s==-1){
			return NULL;
		}
		pseq=_seqdb_get(ppdb->seqdb[s],SEQDB_NEXT);
		if(ppdb->seqdb[1])
			ppdb->next=s?0:1;
	}else{
		s=id&0x1;
		id=id>>1;
		pseq=_seqdb_get(ppdb->seqdb[s],id);
		ppdb->next=-1;
	}
	if(pseq)
		pseq->uid=(pseq->uid<<1)|s;
	return pseq;
}

void seqdb_release(seqdb_t db, seq_t *seq)
{
	pairdb *ppdb=db;
	if(seq){
		int s=seq->uid&0x1;
		seq->uid=seq->uid>>1;
		_seqdb_release(NULL,seq);
	}
}

size_t seqdb_size(seqdb_t db)
{
	pairdb *ppdb=db;
	size_t size=ppdb->seqdb[0]->f_sz;
	if(ppdb->seqdb[1])
		size+=ppdb->seqdb[1]->f_sz;
	return size;
}

off_t seqdb_offset(seqdb_t db)
{
	pairdb *ppdb=db;
	off_t offset=ppdb->seqdb[0]->offset;
	if(ppdb->seqdb[1])
		offset+=ppdb->seqdb[1]->offset;
	return offset;
}

size_t seqdb_seq_len(seqdb_t db)
{
	pairdb *ppdb=db;
	return _seqdb_seq_len(ppdb->seqdb[0]);
}


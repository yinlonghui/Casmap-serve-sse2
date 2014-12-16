#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <unistd.h>

#include "alndb.h"
#include "posdb.h"
#include "sadb.h"
#include "utils.h"

extern int alndb_get(alndb_t db, aln_t *aln);
extern int posdb_put(posdb_t db, pos_t *pos);
extern double get_run_time(void);

int cm_aln2pos_local(sadb_t sadb,alndb_t alndb, posdb_t posdb)
{
	int j;
	aln_t aln;
	pos_t pos;
	unsigned long long n_seq=0,n_aln=0,n_pos=0;
	double last_time,curr_time;
	last_time=curr_time=get_run_time();
	pos.id=SEQ_ID_INVALID;
	while(alndb_get(alndb,&aln)==0){
		/* workaround */
		if(aln.sa==0){
			aln.sa++;
			aln.width--;
			if(aln.width==0)
				continue;
		}

		if(pos.id!=aln.id){
			++n_seq;
			pos.id=aln.id;	
		}
		pos.r=aln.r;
		pos.a=aln.a;
		pos.n_mm=aln.n_mm;
		pos.n_gapo=aln.n_gapo;
		pos.n_gape=aln.n_gape;
		for(j=0;j<aln.width;++j){
			pos.pos=sadb_get(sadb,aln.r,aln.sa+j);
			if(posdb_put(posdb,&pos))
				return -1;
			++n_pos;
		}
		++n_aln;
		curr_time=get_run_time();
		if(curr_time-last_time>1){
			fprintf(stderr,"[aln2pos] %.2f sec, %llu reads, %llu aln, %llu pos\n", curr_time, n_seq, n_aln, n_pos);
			last_time=curr_time;
		}
	}
	curr_time=get_run_time();
	fprintf(stderr,"[aln2pos] %.2f sec, %llu reads, %llu aln, %llu pos\n", curr_time, n_seq, n_aln, n_pos);
	return 0;
}

int cm_aln2pos_remote(int sock,alndb_t alndb, posdb_t posdb)
{
	int rc,j;
	aln_t aln;
	pos_t pos;
	unsigned long long n_seq=0,n_aln=0,n_pos=0;
	double last_time,curr_time;
	last_time=curr_time=get_run_time();
	pos.id=SEQ_ID_INVALID;
	while(alndb_get(alndb,&aln)==0){
		++n_aln;
		/* workaround */
		if(aln.sa==0){
			aln.sa++;
			aln.width--;
			if(aln.width==0)
				continue;
		}

		if(pos.id!=aln.id)
			++n_seq;

		rc=write(sock,&aln,sizeof(aln_t));
		if(rc!=sizeof(aln_t)){
			xerror("write");
			return -1;
		}
		for(j=0;j<aln.width;++j){
			rc=read(sock,&pos,sizeof(pos_t));
			if(rc!=sizeof(pos_t)){
				xerror("xread");
				return -1;
			}

			if(posdb_put(posdb,&pos))
				return -1;
		}
		n_pos+=aln.width;
		curr_time=get_run_time();
		if(curr_time-last_time>1){
			fprintf(stderr,"[aln2pos] %.2f sec, %llu reads, %llu aln, %llu pos\n", curr_time, n_seq, n_aln, n_pos);
			last_time=curr_time;
		}
	}
	curr_time=get_run_time();
	fprintf(stderr,"[aln2pos] %.2f sec, %llu reads, %llu aln, %llu pos\n", curr_time, n_seq, n_aln, n_pos);
	return 0;
}

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "refdb.h"
#include "seq.h"

void usage()
{
	fprintf(stderr,
		"Usage: cm_ref [options] <prefix> <offset> <len>\n"
		"Options:	-c	complement sequence\n"
		);
}

int main(int argc, char **argv)
{
	int rc=0;
	int next_opt;
	refdb_t db;
	char *seq;
	int comp=0;
	seq_sz_t offset,len;
	do{
		next_opt=getopt(argc,argv,"c");
		switch(next_opt){
			case 'c': comp=1;break;
			case -1: break;
			default:usage(); break;
		}
	}while(next_opt!=-1);
	if(argc-optind<3){
		usage();
		return -1;
	}
	offset=atoi(argv[optind+1]);
	len=atoi(argv[optind+2]);
	if(offset==0){
		fprintf(stderr,"Invalid offset: offset starts from 1\n");
		goto ERR_0;
	}
	if(len==0){
		goto ERR_0;
	}

	db=refdb_open(argv[optind],0);
	if(db==NULL){
		fprintf(stderr,"Open reference error\n");
		goto ERR_0;
	}
	seq=malloc(len+1);

	rc=refdb_seq_nt4(db,seq,offset-1,len);
	if(rc){
		fprintf(stderr,"Retrieve seq error\n");
		goto ERR_1;
	}

	seq[len]=0;
	if(comp)
		nt4_to_comp(seq,seq,len);
	nt4_to_seq(seq,seq,len);

	puts(seq);
ERR_1:
	free(seq);
	refdb_close(db);
ERR_0:
	return -1;
}

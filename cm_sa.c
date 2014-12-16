#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>

#include "sadb.h"

void usage()
{
	fprintf(stderr,
		"Usage: cm_sa [options] <prefix> <sa>\n"
		"Options:	-r	reversed sequence\n"
		);
}

int main(int argc, char **argv)
{
	int next_opt;
	sadb_t db;
	int strand=0;
	seq_sz_t sa,pos;
	do{
		next_opt=getopt(argc,argv,"r");
		switch(next_opt){
			case 'r': strand=1;break;
			case -1: break;
			default:usage(); break;
		}
	}while(next_opt!=-1);
	if(argc-optind<2){
		usage();
		return -1;
	}
	sa=atoi(argv[optind+1]);
	db=sadb_open(argv[optind],SADB_IN_MEM);
	if(db==NULL){
		fprintf(stderr,"Open sa file error\n");
		goto ERR_0;
	}
	pos=sadb_get(db,strand,sa);
	sadb_close(db);
	if(pos==SEQ_SZ_MAX){
		fprintf(stderr,"Retrieve POS error\n");
	}else{
		fprintf(stdout,"%u\n",pos);
	}
	return 0;
ERR_0:
	return -1;
}

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include "seqdb.h"

void usage(void)
{
	fprintf(stderr,"Usage: cm_seq <reads> <id>\n");
}

void print_seq(seq_t *seq)
{
	if(seq)
		printf("%s\t%u\t%lu\n%s\n%s\n",seq->name,seq->length,seq->uid,seq->seq,seq->qual);
}

int main(int argc, char **argv)
{
	seqdb_t db;
	seq_t *seq;
	seq_id_t id=SEQDB_NEXT;

	if(argc<2){
		usage();
		return -1;
	}

	if(argc>2)
		id=atoll(argv[2]);

	db=seqdb_open(argv[1],NULL,0);

	seq=seqdb_get(db,id);
	print_seq(seq);
	seqdb_release(db,seq);
	if(id==SEQDB_NEXT){
		while(seq){
			seq=seqdb_get(db,SEQDB_NEXT);
			print_seq(seq);
			seqdb_release(db,seq);
		};
	}
	seqdb_close(db);
	return 0;
}

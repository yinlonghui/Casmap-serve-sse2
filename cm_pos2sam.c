#include <stdlib.h>
#include <stdio.h>
#include <sys/times.h>

#include "seq.h"
#include "utils.h"
#include "seqdb.h"
#include "posdb.h"
#include "refdb.h"
#include "pos2sam.h"

void usage(void)
{
	fprintf(stderr,
	"Usage:  cm_pos2sam [options] <prefix> <*.pos> <reads1.fasta> [reads2.fasta]\n"
	"Options:\n"
	"        -f FILE    output file [stdout]\n"
	"        -m NUM     maximum number of positions to report [3]\n"
	"        -t NUM     number of threads. 0 for no multi-threading [0]\n"
	"        -g         Report mapping distribution (GUI requires this option to be set)\n"
	);
}

int main(int argc, char **argv){
	int rc;
	int next_opt;
	const char *prefix;
	const char *seq1file,*seq2file;
	const char* posfile;
	const char* outfile=NULL;
	double t_start,t_end;
	const char *meta;
	struct pos2sam_opt opt;
	memset(&opt,0,sizeof(opt));
	opt.multi=3;
	opt.report_map=0;
	opt.fout=stdout;
	opt.thread=0;
	opt.buf_len=4096;
	opt.fnr=0.04;

	do{
		next_opt=getopt(argc,argv,"f:m:gt:b:");
		switch(next_opt){
			case 'f': outfile=optarg; break;
			case 'm': opt.multi=atoi(optarg); break;
			case 'g': opt.report_map=1;break;
			case 't': opt.thread=atoi(optarg); break;
			case 'b': opt.buf_len=atoi(optarg); break;
			case -1: 
				  break;
			case ':':
				 fprintf(stderr,"Missing argument\n");
				 return -1;
			default: 
				 fprintf(stderr,"Unknown option -%c %s\n",next_opt,optarg); 
				 return -1;
		}
	}while(next_opt!=-1);

	if(argc-optind<2){
		usage();
		return -1;
	}
	prefix=argv[optind];
	posfile=argv[optind+1];
	seq1file=argv[optind+2];
	if(argc-optind>2)
		seq2file=argv[optind+3];
	else
		seq2file=NULL;

	if(outfile){
		opt.fout=fopen(outfile,"w");
		if(opt.fout==NULL){
			xerror("Open output file");
			return -1;
		}
		//setlinebuf(opt.fout);
	}

	opt.seqdb=seqdb_open(seq1file,seq2file,0);
	if(opt.seqdb==NULL){
		rc=-1;
		goto ERR_0;
	}
	
	init_run_time();

	fprintf(stderr,"[pos2sam] Prefetching database...\n");
	t_start=get_run_time();
	opt.refdb=refdb_open(prefix,REFDB_IN_MEM);
	if(opt.refdb==NULL){
		rc=-1;
		goto ERR_1;
	}

	t_end=get_run_time();
	fprintf(stderr,"[pos2sam] Database prefetched in %.2f sec\n",t_end-t_start);

	opt.posdb=posdb_open(posfile);
	if(opt.posdb==NULL){
		rc=-1;
		goto ERR_2;
	}

	meta=posdb_meta_info(opt.posdb);
	meta=strstr(meta,"fnr=");
	if(meta){
		float f;
		if(sscanf(meta+4,"%f",&f)==1)
			opt.fnr=f;
	}
	if(opt.fnr>=1 || opt.fnr==0) 
		fprintf(stderr,"[pos2sam] MAX_DIFF: %d\n",(int)opt.fnr);
	else
		fprintf(stderr,"[pos2sam] FNR: %f\n",opt.fnr);

	t_start=get_run_time();
	fprintf(stderr,"[pos2sam] Translating to SAM...\n");

	rc=cm_pos2sam(&opt);
	if(rc){
		rc=-1;
		goto ERR_3;
	}

	t_end=get_run_time();
	fprintf(stderr,"[pos2sam] Translated to SAM in %.2f sec\n",t_end-t_start);
ERR_3:
	posdb_close(opt.posdb);
ERR_2:
	refdb_close(opt.refdb);
ERR_1:
	seqdb_close(opt.seqdb);
ERR_0:
	fclose(opt.fout);
	return rc;
}

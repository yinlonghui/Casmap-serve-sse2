#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <math.h>
#include "seq.h"
#include "refdb.h"
#include "utils.h"
#include "seqdb.h"
#include "posdb.h"
#include "pos2sampe.h"

extern pe_opt_t *bwa_init_pe_opt();
extern int  pipe_cm_pos2sampe(refdb_t refdb,posdb_t posdb,seqdb_t seqdb,FILE *fout,pe_opt_t  *pe_opt);
int  cm_pos2sampe(refdb_t refdb,posdb_t posdb[2],seqdb_t seqdb[2],FILE *fout,pe_opt_t *pe_opt);

static int usage(){
	fprintf(stderr, "Usage:  cm_pos2sampe [options] <prefix> <in1.pos> <in2.pos> <in1.fq> <in2.fq>\n\n");
	fprintf(stderr, "Options: 	-f	output file [stdout]\n");
	fprintf(stderr, "               -m      maximum hits to output for paired reads");
	return 1;
}
int main(int argc,char *argv[])
{
		pe_opt_t  *pe_opt;
		int c,pipe_opt = 0;
		FILE *fout=stdout;
		char *file_out=NULL;
		const char *seqfile[2];
		const char *posfile[2];
		const char *prefix;
		posdb_t posdb[2];
		seqdb_t seqdb[2];
		refdb_t refdb;
		double t_start,t_end;
		
		posdb[0] = posdb[1] =  NULL;
		seqfile[0] = seqfile[1] = NULL;
		posfile[0] = posfile[1] = NULL;

		pe_opt=bwa_init_pe_opt();
		
		while ((c = getopt(argc, argv, "f:t:Pm:")) >= 0){
			switch(c){
			case('f'):	file_out = optarg;break;
			case('t'):	{	
						pe_opt->n_thread =atoi(optarg);
						fprintf(stderr,"[pos2sam]number thread is%d\n",pe_opt->n_thread);
						break;
					}
			case('m'):	pe_opt->n_multi = atoi(optarg); break;
			case('P'):	pipe_opt=1;break;
			}
		}
		
		if(file_out){
			fout=fopen(file_out,"w");
			if(fout==NULL){
				xerror("Open out file");
				return -1;
			}
		}
		
		//fprintf(stderr,"pipe_opt=%d\n",pipe_opt);
		init_run_time();
		
		if(pipe_opt){
			if(argc-optind<3)
				return usage();
			prefix = argv[optind];
			seqfile[0] = argv[optind+1];
			seqfile[1] = argv[optind+2];
			//posfile[0] = argv[optind+3];	
			
			t_start=get_run_time();
			fprintf(stderr,"[pos2sampe] Prefetching database...\n");
			refdb=refdb_open(prefix,REFDB_IN_MEM);
			if(refdb ==NULL){
				fprintf(stderr,"Can't open prefix\n");
				return -1;
			}
			t_end=get_run_time();
			fprintf(stderr,"[pos2sampe] Database prefetched in %.2f sec\n",t_end-t_start);

			posdb[0] = posdb_open(posfile[0]);
			
			if(posdb==NULL){
				refdb_close(refdb);
				fprintf(stderr, "Can't open posfile\n");
				return -1;
			}
			
			seqdb[0] = seqdb_open(seqfile[0],seqfile[1],0);
			
			if(seqdb[0]==NULL){
				refdb_close(refdb);
				posdb_close(posdb[0]);
				fprintf(stderr,"Can't open readsfile\n");
			}

			t_start=get_run_time();
			fprintf(stderr,"[pos2sampe] Translating to SAM...\n");

			pipe_cm_pos2sampe(refdb,posdb[0],seqdb[0],fout,pe_opt);
			refdb_close(refdb);
			posdb_close(posdb[0]);
			seqdb_close(seqdb[0]);
			t_end=get_run_time();
			fprintf(stderr,"[pos2sampe] Translated to SAM in %.2f sec\n",t_end-t_start);
			
		}else{
			if(argc-optind<5)
				return usage();
			prefix = argv[optind];
			posfile[0] = argv[optind+1];
			posfile[1] = argv[optind+2];
			seqfile[0] = argv[optind+3];
			seqfile[1] = argv[optind+4];
			
			t_start=get_run_time();
			fprintf(stderr,"[pos2sampe] Prefetching database...\n");
			refdb=refdb_open(prefix,REFDB_IN_MEM);
			if(refdb ==NULL){
				fprintf(stderr,"Can't open prefix\n");
				return -1;
			}
			t_end=get_run_time();
			fprintf(stderr,"[pos2sampe] Database prefetched in %.2f sec\n",t_end-t_start);
			
			posdb[0] = posdb_open(posfile[0]);
			
			if(posdb[0]==NULL){
				refdb_close(refdb);
				fprintf(stderr,"can't open posfile");
				return -1;
			}
			
			posdb[1] = posdb_open(posfile[1]);
			
			if(posdb[1]==NULL){
				refdb_close(refdb);
				posdb_close(posdb[0]);
				fprintf(stderr,"can't open posfile");
				return -1;
			}
			seqdb[0] = seqdb_open(seqfile[0],NULL,0);
			if(seqdb[0]==NULL){
				refdb_close(refdb);
				posdb_close(posdb[0]);
				posdb_close(posdb[1]);
				fprintf(stderr,"can't open seqfile");
				return -1;
			}
			seqdb[1] = seqdb_open(seqfile[1],NULL,0);
			if(seqdb[1]==NULL){
				refdb_close(refdb);
				posdb_close(posdb[0]);
				posdb_close(posdb[1]);
				seqdb_close(seqdb[0]);
				fprintf(stderr,"can't open seqfile");
				return -1;
			}

			t_start=get_run_time();
			fprintf(stderr,"[pos2sampe] Translating to SAM...\n");

			cm_pos2sampe(refdb,posdb,seqdb,fout,pe_opt);
			refdb_close(refdb);
			posdb_close(posdb[0]);
			posdb_close(posdb[1]);
			seqdb_close(seqdb[0]);
			seqdb_close(seqdb[1]);

			t_end=get_run_time();
			fprintf(stderr,"[pos2sampe] Translated to SAM in %.2f sec\n",t_end-t_start);

		}
		FREE(pe_opt);
		fclose(fout);
		return 0;
}


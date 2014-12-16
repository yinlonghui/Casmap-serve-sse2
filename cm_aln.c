#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "aln.h"
#include "utils.h"
#include "seqdb.h"
#include "alndb.h"

/*
 *	This code Support two socket port . The optional command input FPGA'S IP address. eg:('-1','-2').
*/
static void usage(void)
{
	fprintf(stderr,
		"Usage:  cm_aln [options] <server> <reads1.fasta> [reads2.fasta]\n"
		"Options:\n"
		"         -n INT    maxium allowed differences [3]\n"
		"         -o INT    maxium allowed gap opens [1]\n"
		"         -e INT    maxium allowed gap extensions [6]\n"
		"         -i INT    indel skip [5]\n"
		"         -d INT    maxium occurrences for extending a long deletion [10]\n"
	        "         -l INT    seed length [32]\n"
		"         -k INT    maxium differences in the seed [2]\n"
		"         -M INT    mismatch penalty [3]\n"
       		"         -O INT    gap open penalty [11]\n"
		"         -E INT    gap extension penalty [4]\n"
		"         -R [1-3]  alignment effort [2]\n"
		"                   [1] - only the first alignment (fast)\n"
		"                   [2] - best alignments (normal)\n"
		"                   [3] - all alignments (slow)\n"
		"         -f STR    output file prefix [<input>.aln]\n"
		"         -T INT    maxium iteration limit [200000]\n"
		"         -p port   server service port[9999]\n"
		"	  -1   <server>\n"
		"	  -2   <server>\n"
	       );
}

int main(int argc, char **argv)
{
	int rc;
	int next_opt;
	struct aln_opt opt;
	const char *out_prefix=NULL;
	char *msz_str=NULL;
	unsigned long long mem_sz=sysconf(_SC_PHYS_PAGES);
	seqdb_t seqdb;
	alndb_t alndb;
	unsigned int sa_intv_bits;

	mem_sz*=sysconf(_SC_PAGE_SIZE);
	mem_sz--; //do not consume all memory
	memset(&opt,0,sizeof(opt));
	opt.port=9999;

	opt.fnr=0.04;

	opt.aln_param.MAX_SEED_DIFF=2;
	opt.aln_param.MAX_GAPO=1;
	opt.aln_param.MAX_GAPE=6;
	opt.aln_param.SEED_WIDTH=32;
	opt.aln_param.MAX_DEL_OCC=10;
	opt.aln_param.S_MM=3;
	opt.aln_param.S_GAPO=11;
	opt.aln_param.S_GAPE=4;
	opt.aln_param.FAST_MODE=1;
	opt.aln_param.INDEL_END_SKIP=5;
	opt.aln_param.EFFORT=BWA_AUTO;
	opt.aln_param.RID_WIDTH=7;
	opt.aln_param.MAX_BACK=200000;
	opt.server[0] =  opt.server[1] = NULL ;

	do{
		next_opt=getopt(argc,argv,"1:2:n:o:e:i:d:l:k:M:O:E:R:f:a:C:vp:qB:T:p:g");
		switch(next_opt){
			case '1': opt.server[0] = optarg ; break ;
			case '2': opt.server[1] = optarg ; break ;
			case 'n': opt.fnr=atof(optarg); break;
			case 'o': opt.aln_param.MAX_GAPO=atoi(optarg); break;
			case 'e': opt.aln_param.MAX_GAPE=atoi(optarg); break;
			case 'i': opt.aln_param.INDEL_END_SKIP=atoi(optarg); break;
			case 'd': opt.aln_param.MAX_DEL_OCC=atoi(optarg); break;
			case 'l': opt.aln_param.SEED_WIDTH=atoi(optarg); break;
			case 'k': opt.aln_param.MAX_SEED_DIFF=atoi(optarg); break;
			case 'M': opt.aln_param.S_MM=atoi(optarg); break;
			case 'O': opt.aln_param.S_GAPO=atoi(optarg); break;
			case 'E': opt.aln_param.S_GAPE=atoi(optarg); break;
			case 'R': 
				  opt.aln_param.EFFORT=atoi(optarg);
				  if(opt.aln_param.EFFORT>3||opt.aln_param.EFFORT<1){
					  fprintf(stderr,"Illegal -R value %d. Must be 1, 2 or 3\n",opt.aln_param.EFFORT);
					  rc=-1;
					  goto ERR0;
				  }
				  break;
			case 'f': out_prefix=optarg;break;
			case 'a': opt.aln_param.FAST_MODE=atoi(optarg); break;
			case 'C': {
					  int max_c = atoi(optarg); 
					  int i,t;
					  for(i=0,t=1;t < max_c && t<1024;++i,t*=2);
					  opt.aln_param.RID_WIDTH=i;
					  break;
				  }
			case 'v': opt.verbose=1;break;
			case 'g': opt.debug=1;break;
			case 'q': opt.quiet=1;break;
			case 'p': opt.port=atoi(optarg);break;
			case 'B': msz_str=optarg; break;
			case 'T': opt.aln_param.MAX_BACK=atoi(optarg);break;
			case -1:break;
			default:
				fprintf(stderr,"Unknown option -%c %s\n",next_opt,optarg);
				usage();
				rc=-1;
				goto ERR0;
		}

	}while(next_opt!=-1);

	if(argc-optind<1){
		usage();
		rc=-1;
		goto ERR0;
	}

//	opt.server=argv[optind];
	opt.seq1_prefix=argv[optind];
	if(argc-optind>1)
		opt.seq2_prefix=argv[optind+1];
	else
		opt.seq2_prefix=NULL;

	if(out_prefix==NULL){
		if(opt.seq2_prefix){
			opt.aln_prefix=malloc(strlen(opt.seq1_prefix)+strlen(opt.seq2_prefix)+6);
			sprintf(opt.aln_prefix,"%s-%s.aln",opt.seq1_prefix,opt.seq2_prefix);
		}else{
			opt.aln_prefix=malloc(strlen(opt.seq1_prefix)+5);
			sprintf(opt.aln_prefix,"%s.aln",opt.seq1_prefix);
		}
	}else{
		opt.aln_prefix=malloc(strlen(out_prefix)+5);
		sprintf(opt.aln_prefix,"%s.aln",out_prefix);
	}

	if(msz_str && str_to_size_ll(msz_str,&mem_sz)){
		fprintf(stderr,"Wrong format -B %s\n",msz_str);
		return -1;
	}

	for(sa_intv_bits=0;(1llu<<(sa_intv_bits+3))<=mem_sz;sa_intv_bits++);
	sa_intv_bits--;
	fprintf(stderr,"[seq2aln] MEM SZ: %llu, SA INTV: %llu\n",mem_sz,1llu<<sa_intv_bits);

	seqdb=seqdb_open(opt.seq1_prefix,opt.seq2_prefix,0);
	if(seqdb==NULL){
		return -1;
	}

	alndb=alndb_open(opt.aln_prefix,sa_intv_bits);
	if(alndb==NULL){
		return -1;
	}

	{
		char *meta_str=malloc(32);
		sprintf(meta_str,"fnr=%f",opt.fnr);
		alndb_meta_append(alndb,meta_str);
		free(meta_str);
	}

	init_run_time();

	rc=cm_aln(seqdb,alndb,&opt);

	seqdb_close(seqdb);
	alndb_close(alndb);

	free(opt.aln_prefix);
ERR0:
	return rc;
}

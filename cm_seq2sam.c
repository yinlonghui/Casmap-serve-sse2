#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <pthread.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netdb.h>

#include "aln.h"
#include "utils.h"
#include "refdb.h"
#include "sadb.h"
#include "seqdb.h"
#include "alndb.h"
#include "posdb.h"
#include "pos2sampe.h"
#include "pos2sam.h"

struct seq2aln_param{
	seqdb_t seqdb;
	alndb_t alndb;
	struct aln_opt *opt;
};

struct aln2pos_param{
	alndb_t alndb;
	posdb_t posdb;
	sadb_t sadb;
	int remote;
	const char *hostname;
	short port;
};

struct pos2sam_param{
	refdb_t refdb;
	seqdb_t seqdb;
	posdb_t posdb;
	FILE* fout;
	int multi;
	int pe;
	int n_thread;
	float fnr;
	int report_map;
};

extern pe_opt_t *bwa_init_pe_opt();
extern int  pipe_cm_pos2sampe(refdb_t refdb,posdb_t posdb,seqdb_t seqdb,FILE *fout,pe_opt_t  *pe_opt);
extern int cm_aln2pos_local(sadb_t sadb,alndb_t alndb, posdb_t posdb);
extern int cm_aln2pos_remote(int sock,alndb_t alndb, posdb_t posdb);
extern void alndb_shutdown(alndb_t db);
extern void posdb_shutdown(alndb_t db);

static void usage(void)
{
	fprintf(stderr,
		"Usage:  cm_seq2sam [options]  <prefix> <reads1.fasta> [reads2.fasta]\n\n"
		"Options: -n INT    maxium allowed differences [3]\n"
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
		"         -m INT     maximum number of positions to report [3]\n"
		"         -t INT    threads to use during SAMPE stage\n"
		"         -g Report mapping distribution (GUI requires this option to be set)\n"
		"	  -1 <server> \n"
		"	  -2 <server> \n"
	       );
}

static void* seq2aln_worker(void *data)
{
	struct seq2aln_param *pm=(struct seq2aln_param*)data;
	if(cm_aln(pm->seqdb,pm->alndb,pm->opt))
		return (void*)-1;
	return 0;
}

static void* aln2pos_worker(void *data)
{
	void *rv;
	struct aln2pos_param *pm=(struct aln2pos_param*)data;
	if(pm->remote){
		int sock;
		struct sockaddr_in sa;
		struct hostent *he;
		he=gethostbyname(pm->hostname);
		if(he==NULL){
			herror("gethostbyname\n");
			rv=(void*)-1;
			goto ERR_0;
		}
		memset(&sa,0,sizeof(sa));
		sa.sin_family=AF_INET;
		sa.sin_port=htons(pm->port);
		memcpy(&(sa.sin_addr.s_addr),he->h_addr_list[0],sizeof(in_addr_t));
		sock=socket(PF_INET,SOCK_STREAM,0);
		if(sock==-1){
			xerror("socket");
			rv=(void*)-1;
			goto ERR_0;
		}
		if(connect(sock,(struct sockaddr*)(&sa),sizeof(sa))==-1){
			xerror("connect");
			rv=(void*)-1;
			goto ERR_1;
		}
		if(cm_aln2pos_remote(sock,pm->alndb,pm->posdb)){
			rv=(void*)-1;
			goto ERR_1;
		}
ERR_1:
		close(sock);
	}else{
		if(cm_aln2pos_local(pm->sadb,pm->alndb,pm->posdb)){
			rv=(void*)-1;
			goto ERR_0;
		}
	}
ERR_0:
	return rv;
}

static void* pos2sam_worker(void *data)
{
	struct pos2sam_param *pm=(struct pos2sam_param*)data;
	if(pm->pe){
		pe_opt_t *pe_opt=bwa_init_pe_opt();
		pe_opt->n_thread=pm->n_thread;
		pe_opt->n_multi = pm->multi;
		if(pipe_cm_pos2sampe(pm->refdb,pm->posdb,pm->seqdb,pm->fout,pe_opt)){
			return (void*)-1;
		}
		free(pe_opt);
	}else{
		struct pos2sam_opt opt;
		opt.refdb=pm->refdb;
		opt.seqdb=pm->seqdb;
		opt.posdb=pm->posdb;
		opt.fout=pm->fout;
		opt.multi=pm->multi;
		opt.thread=pm->n_thread;
		opt.report_map=pm->report_map;
		opt.buf_len=4096;
		opt.fnr=pm->fnr;
		if(cm_pos2sam(&opt))
			return (void*)-1;
	}
	return 0;
}

int main(int argc, char **argv)
{
	int rc;
	int next_opt;
	struct aln_opt opt;
	const char *out_prefix=NULL;
	const char *prefix;
	char *msz_str=NULL;
	unsigned long long mem_sz=sysconf(_SC_PHYS_PAGES);
	seqdb_t seqdb1,seqdb2;
	refdb_t refdb;
	sadb_t sadb;
	alndb_t alndb;
	posdb_t posdb;
	FILE *fout=stdout;
	int multi_pos=3;
	int a2p_remote=0;
	const char *a2p_hostname="localhost";
	short a2p_port=9990;
	struct seq2aln_param seq2aln_pm;
	struct aln2pos_param aln2pos_pm;
	struct pos2sam_param pos2sam_pm;
	pthread_t t_seq2aln,t_aln2pos,t_pos2sam;
	void *ret;
	int sam_n_thread=1;
	unsigned int sa_intv_bits=32;
	
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

	pos2sam_pm.report_map=0;

	do{
		next_opt=getopt(argc,argv,"1:2:n:o:e:i:d:l:k:M:O:E:R:f:a:C:vp:qB:T:p:u:Gm:t:rg");
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
			case 'G': opt.debug=1;break;
			case 'q': opt.quiet=1;break;
			case 'p': opt.port=atoi(optarg);break;
			case 'u': a2p_port=atoi(optarg);break;
			case 'B': msz_str=optarg; break;
			case 'T': opt.aln_param.MAX_BACK=atoi(optarg);break;
			case 'm': multi_pos=atoi(optarg);break;
			case 'r': a2p_remote=1;break;
			case 't': sam_n_thread=atoi(optarg);break;
			case 'g': pos2sam_pm.report_map=1;break;
			case -1:break;
			default:
				fprintf(stderr,"Unknown option -%c %s\n",next_opt,optarg);
				usage();
				rc=-1;
				goto ERR0;
		}

	}while(next_opt!=-1);

	if(argc-optind<2){
		usage();
		rc=-1;
		goto ERR0;
	}

//	opt.server=argv[optind];
	prefix=argv[optind];
	opt.seq1_prefix=argv[optind+1];
	if(argc-optind>2)
		opt.seq2_prefix=argv[optind+2];
	else
		opt.seq2_prefix=NULL;
	
	if(out_prefix){
		fout=fopen(out_prefix,"w");
		if(fout==NULL){
			xerror("Open output file");
			return -1;
		}
	}

	if(msz_str && str_to_size_ll(msz_str,&mem_sz)){
		fprintf(stderr,"Wrong format -B %s\n",msz_str);
		return -1;
	}

	init_run_time();

	fprintf(stderr,"[seq2sam] MEM SZ: %llu, SA INTV: %llu\n",mem_sz,1llu<<sa_intv_bits);

	seqdb1=seqdb_open(opt.seq1_prefix,opt.seq2_prefix,0);
	if(seqdb1==NULL){
		return -1;
	}

	seqdb2=seqdb_open(opt.seq1_prefix,opt.seq2_prefix,0);
	if(seqdb2==NULL){
		return -1;
	}

	refdb=refdb_open(prefix,REFDB_IN_MEM);
	if(refdb==NULL){
		return -1;
	}

	if(a2p_remote){
		sadb=NULL;
	}else{
		sadb=sadb_open(prefix,(SADB_IN_MEM|sa_intv_bits));
		if(sadb==NULL){
			return -1;
		}
	}

	alndb=alndb_open(NULL,0);
	if(alndb==NULL){
		return -1;
	}

	posdb=posdb_open(NULL);
	if(posdb==NULL){
		return -1;
	}
	
	seq2aln_pm.seqdb=seqdb1;
	seq2aln_pm.alndb=alndb;
	seq2aln_pm.opt=&opt;
	
	aln2pos_pm.alndb=alndb;
	aln2pos_pm.posdb=posdb;
	aln2pos_pm.sadb=sadb;
	aln2pos_pm.remote=a2p_remote;
	aln2pos_pm.hostname=a2p_hostname;
	aln2pos_pm.port=a2p_port;

	pos2sam_pm.refdb=refdb;
	pos2sam_pm.seqdb=seqdb2;
	pos2sam_pm.posdb=posdb;
	pos2sam_pm.fout=fout;
	pos2sam_pm.multi=multi_pos;
	pos2sam_pm.n_thread=sam_n_thread;
	if(opt.seq2_prefix)
		pos2sam_pm.pe=1;
	else
		pos2sam_pm.pe=0;
	pos2sam_pm.fnr=opt.fnr;

	rc=pthread_create(&t_seq2aln,NULL,seq2aln_worker,&seq2aln_pm);
	if(rc){
		xerror("create thread");
		goto ERR0;
	}
	rc=pthread_create(&t_aln2pos,NULL,aln2pos_worker,&aln2pos_pm);
	if(rc){
		xerror("create thread");
		goto ERR0;
	}
	rc=pthread_create(&t_pos2sam,NULL,pos2sam_worker,&pos2sam_pm);
	if(rc){
		xerror("create thread");
		goto ERR0;
	}

	rc=pthread_join(t_seq2aln,&ret);
	if(rc){
		xerror("join thread");
		goto ERR0;
	}

	alndb_shutdown(alndb);
	seqdb_close(seqdb1);

	rc=pthread_join(t_aln2pos,&ret);
	if(rc){
		xerror("join thread");
		goto ERR0;
	}

	posdb_shutdown(posdb);
	alndb_close(alndb);
	if(a2p_remote==0)
		sadb_close(sadb);

	rc=pthread_join(t_pos2sam,&ret);
	if(rc){
		xerror("join thread");
		goto ERR0;
	}

	posdb_close(posdb);
	refdb_close(refdb);
	seqdb_close(seqdb2);

	fclose(fout);
ERR0:
	return rc;
}

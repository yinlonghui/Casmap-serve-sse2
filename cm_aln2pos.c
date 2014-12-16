#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/times.h>
#include <fcntl.h>
#include <poll.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netdb.h>

#include "seq.h"
#include "alndb.h"
#include "sadb.h"
#include "posdb.h"
#include "utils.h"

#define POS_BUF_SIZE	(1024)

void usage(void)
{
	fprintf(stderr,
	"Usage:  cm_aln2pos [options] <prefix> <*.aln>\n"
	"Options:\n"
	"               -f STR       output file prefix [<prefix>.aln]\n"
	"               -N 	     no buffering (sloooow but suit for small data)\n"
	);
}

int cm_aln2pos_remote(int sock, const char *aln_fn, const char *pos_fn, unsigned long long *pos_n)
{
	int rc;
	int i,j;
	int aln_fd,pos_fd;
	static aln_t aln_buf[POS_BUF_SIZE];
	static pos_t pos_buf[POS_BUF_SIZE];
	unsigned int pos_offset=0;
	unsigned int aln_sz=0;
	aln_t *paln;
	seq_id_t last_seq_id, seq_id;

	double last_time,curr_time;
	static unsigned long long n_seq=0,n_aln=0,n_pos=0;
	last_time=curr_time=get_run_time();

	*pos_n=0;

	aln_fd=open(aln_fn,O_RDONLY);
	if(aln_fd==-1){
		xerror(aln_fn);
		return -1;
	}

	pos_fd=open(pos_fn,O_WRONLY|O_TRUNC|O_CREAT,
			S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH);
	if(pos_fd==-1){
		xerror(pos_fn);
		return -1;
	}

	last_seq_id=SEQ_ID_INVALID;

	do{
		aln_sz=0;
		rc=xread(aln_fd,aln_buf,sizeof(aln_t)*POS_BUF_SIZE);
		if(rc==-1){
			xerror("Read aln part file");
			return -1;
		}		
		aln_sz=rc/sizeof(aln_t);
		for(i=0,paln=aln_buf;i<aln_sz;++i,++paln){
			/* workaround for SAI==0 */
			if(paln->sa==0){
				paln->sa++;
				paln->width--;
				if(paln->width==0)
					continue;
			}

			seq_id=SEQ_ID(paln->id);
			if(seq_id!=last_seq_id)
				++n_seq;
			++n_aln;
			n_pos+=paln->width;

			rc=write(sock,paln,sizeof(aln_t));
			if(rc!=sizeof(aln_t)){
				xerror("xwrite");
				return -1;
			}
			for(j=0;j<paln->width;++j){
				rc=read(sock,pos_buf+pos_offset,sizeof(pos_t));
				if(rc!=sizeof(pos_t)){
					xerror("xread");
					return -1;
				}
				++pos_offset;
				if(pos_offset==POS_BUF_SIZE){
					rc=write(pos_fd,pos_buf,sizeof(pos_t)*POS_BUF_SIZE);
					if(rc==-1){
						xerror("Write pos part file");
						return -1;
					}
					pos_offset=0;
				}
			}
		}
		curr_time=get_run_time();
		if(curr_time-last_time>1){
			fprintf(stderr,"[aln2pos] %.2f sec, %llu reads, %llu aln, %llu pos\n", curr_time, n_seq, n_aln, n_pos);
			last_time=curr_time;
		}
	}while(aln_sz==POS_BUF_SIZE);
	if(pos_offset){
		rc=write(pos_fd,pos_buf,sizeof(pos_t)*pos_offset);
		if(rc==-1){
			xerror("Write pos part file");
			return -1;
		}
	}
	*pos_n=n_pos;
	if(aln_fn) close(aln_fd);
	if(pos_fn) close(pos_fd);
	return 0;
}

int cm_aln2pos_core(sadb_t sadb, const char *aln_fn, const char *pos_fn, unsigned long long *pos_n)
{
	int rc;
	int i,j;
	int aln_fd,pos_fd;
	static aln_t aln_buf[POS_BUF_SIZE];
	static pos_t pos_buf[POS_BUF_SIZE];
	unsigned int pos_offset=0;
	unsigned int aln_sz=0;
	aln_t *paln;
	pos_t *ppos;
	seq_id_t last_seq_id, seq_id;
	uint32_t ref_pos;

	double last_time,curr_time;
	static unsigned long long n_seq=0,n_aln=0,n_pos=0;
	last_time=curr_time=get_run_time();

	*pos_n=0;

	if(aln_fn){
		aln_fd=open(aln_fn,O_RDONLY);
		if(aln_fd==-1){
			xerror(aln_fn);
			return -1;
		}
	}else{
		aln_fd=STDIN_FILENO;
	}

	if(pos_fn){
		pos_fd=open(pos_fn,O_WRONLY|O_TRUNC|O_CREAT,
				S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH);
		if(pos_fd==-1){
			xerror(pos_fn);
			return -1;
		}
	}else{
		pos_fd=STDOUT_FILENO;
	}

	ppos=pos_buf;

	last_seq_id=SEQ_ID_INVALID;

	do{
		aln_sz=0;
		rc=xread(aln_fd,aln_buf,sizeof(aln_t)*POS_BUF_SIZE);
		if(rc==-1){
			xerror("Read aln part file");
			return -1;
		}		
		aln_sz=rc/sizeof(aln_t);
		for(i=0,paln=aln_buf;i<aln_sz;++i,++paln){
			/* workaround for SAI==0 */
			if(paln->sa==0){
				paln->sa++;
				paln->width--;
				if(paln->width==0)
					continue;
			}

			seq_id=SEQ_ID(paln->id);
			if(seq_id!=last_seq_id)
				++n_seq;
			++n_aln;
			n_pos+=paln->width;

			for(j=0;j<paln->width;++j){

				if(paln->sa+j==0) continue;

				ref_pos=sadb_get(sadb,paln->r,paln->sa+j);
				if(rc==SEQ_SZ_MAX) 
					return -1;
				ppos->id=paln->id;
				ppos->r=paln->r;
				ppos->a=paln->a;
				ppos->n_mm=paln->n_mm;
				ppos->n_gapo=paln->n_gapo;
				ppos->n_gape=paln->n_gape;
				ppos->pos=ref_pos;
				++pos_offset;
				++ppos;
				if(pos_offset==POS_BUF_SIZE){
					rc=write(pos_fd,pos_buf,sizeof(pos_t)*POS_BUF_SIZE);
					if(rc==-1){
						xerror("Write pos part file");
						return -1;
					}
					pos_offset=0;
					ppos=pos_buf;
				}
			}
		}
		curr_time=get_run_time();
		if(curr_time-last_time>1){
			fprintf(stderr,"[aln2pos] %.2f sec, %llu reads, %llu aln, %llu pos\n", curr_time, n_seq, n_aln, n_pos);
			last_time=curr_time;
		}
	}while(aln_sz==POS_BUF_SIZE);
	if(pos_offset){
		rc=write(pos_fd,pos_buf,sizeof(pos_t)*pos_offset);
		if(rc==-1){
			xerror("Write pos part file");
			return -1;
		}
	}
	*pos_n=n_pos;
	if(aln_fn) close(aln_fd);
	if(pos_fn) close(pos_fd);
	return 0;
}

int cm_aln2pos(sadb_t sadb, const char* aln_prefix, const char *pos_prefix)
{
	int rc;
	int pos_part_no;
	int aln_part_no;
	unsigned long long pos_part_n;
	unsigned long long aln_part_n;
	char *aln_part_fn;
	char *pos_part_fn;
	FILE *f_aln_meta=NULL;
	FILE *f_pos_meta=NULL;
	double t_start,t_end;

	aln_part_fn=malloc(strlen(aln_prefix)+32);
	pos_part_fn=malloc(strlen(pos_prefix)+32);

	f_aln_meta=fopen(aln_prefix,"r");
	if(f_aln_meta==NULL){
		//xerror("Open aln file");
		xerror(aln_prefix);
		rc=-1;
		goto ERR_0;
	}
	f_pos_meta=fopen(pos_prefix,"w");
	if(f_pos_meta==NULL){
		//xerror("Open pos file");
		xerror(aln_prefix);
		rc=-1;
		goto ERR_0;
	}
	pos_part_no=0;
	while(fscanf(f_aln_meta,"%u%llu",&aln_part_no,&aln_part_n)==2){
		sprintf(aln_part_fn,"%s.%u",aln_prefix,aln_part_no);
		sprintf(pos_part_fn,"%s.%u",pos_prefix,++pos_part_no);
		fprintf(stderr,"[aln2pos] Processing %s\n", aln_part_fn);
		t_start=get_run_time();
		rc=cm_aln2pos_core(sadb,aln_part_fn,pos_part_fn,&pos_part_n);
		if(rc)
			goto ERR_0;
		fprintf(f_pos_meta,"%u\t%llu\n",pos_part_no,pos_part_n);
		t_end=get_run_time();
		fprintf(stderr,"[aln2pos] Finished %s in %.2f sec\n", 
			aln_part_fn,t_end-t_start);
	}

ERR_0:
	fclose(f_pos_meta);
	fclose(f_aln_meta);

	free(aln_part_fn);
	free(pos_part_fn);
	return rc;
}

int main(int argc, char **argv){
	int rc=0;
	int next_opt;
	sadb_t sadb;
	int flags=0;
	const char *prefix;
	const char *aln_prefix;
	char *pos_prefix=NULL;
	const char *out_prefix=NULL;
	int pos_part_no;
	int aln_part_no;
	unsigned long long pos_part_n;
	unsigned long long aln_part_n;
	char *aln_part_fn;
	char *pos_part_fn;
	FILE *f_aln_meta=NULL;
	FILE *f_pos_meta=NULL;
	double t_start,t_pstart,t_end;
	const char *hostname="localhost";
	short port=9990;
	int remote=0;

	int no_buf=0;

	unsigned int sa_intv_bits;
	char *msz_str=NULL;
	unsigned long long mem_sz=sysconf(_SC_PHYS_PAGES);
	mem_sz*=sysconf(_SC_PAGE_SIZE);
	//mem_sz--; //do not consume all memory

	init_run_time();

	do{
		next_opt=getopt(argc,argv,"f:n:B:Nrp:h:");
		switch(next_opt){
			case 'B': msz_str=optarg; break;
			case 'f': out_prefix=optarg;break;
			case 'N': no_buf=1;break;
			case 'r': remote=1;break;
			case 'p': port=atoi(optarg);break;
			case 'h': hostname=optarg;break;
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

	if(argc-optind<1){
		usage();
		return -1;
	}

	t_start=get_run_time();


	if(remote){
		sadb=NULL;
		prefix=NULL;
		aln_prefix=argv[optind];
	}else{
		prefix=argv[optind];

		if(msz_str && str_to_size_ll(msz_str,&mem_sz)){
			fprintf(stderr,"Wrong format -B %s\n",msz_str);
			return -1;
		}

		for(sa_intv_bits=0;(1llu<<(sa_intv_bits+3))<=mem_sz;sa_intv_bits++);
		sa_intv_bits--;
		fprintf(stderr,"[aln2pos] MEM SZ: %llu, SA INTV: %llu\n",mem_sz,1llu<<sa_intv_bits);

		flags=sa_intv_bits;

		if(no_buf==0)
			flags|=SADB_IN_MEM;
		sadb=sadb_open(prefix,flags);
		if(sadb==NULL){
			goto ERR_0;
		}
		aln_prefix=argv[optind+1];
	}

	if(out_prefix==NULL)
		out_prefix=aln_prefix;
	pos_prefix=malloc(strlen(out_prefix)+5);
	sprintf(pos_prefix,"%s.pos",out_prefix);

	aln_part_fn=malloc(strlen(aln_prefix)+32);
	pos_part_fn=malloc(strlen(pos_prefix)+32);

	f_aln_meta=fopen(aln_prefix,"r");
	if(f_aln_meta==NULL){
		xerror(aln_prefix);
		rc=-1;
		goto ERR_2;
	}
	f_pos_meta=fopen(pos_prefix,"w");
	if(f_pos_meta==NULL){
		xerror(aln_prefix);
		rc=-1;
		goto ERR_2;
	}

	{
		int c;
		char *buf=NULL;
		size_t buf_sz=0;
		for(c=fgetc(f_aln_meta);c=='#';c=fgetc(f_aln_meta)){
			if(getline(&buf,&buf_sz,f_aln_meta)!=-1)
				fprintf(f_pos_meta,"#%s",buf);
		}
		ungetc(c,f_aln_meta);
		free(buf);
	}

	pos_part_no=0;
	while(fscanf(f_aln_meta,"%u%llu",&aln_part_no,&aln_part_n)==2){
		sprintf(aln_part_fn,"%s.%u",aln_prefix,aln_part_no);
		sprintf(pos_part_fn,"%s.%u",pos_prefix,++pos_part_no);
		fprintf(stderr,"[aln2pos] Processing %s\n", aln_part_fn);
		t_pstart=get_run_time();
		if(remote){
			int sock;
			struct sockaddr_in sa;
			struct hostent *he;
			he=gethostbyname(hostname);
			if(he==NULL){
				herror("gethostbyname\n");
				goto ERR_2;
			}
			memset(&sa,0,sizeof(sa));
			sa.sin_family=AF_INET;
			sa.sin_port=htons(port);
			memcpy(&(sa.sin_addr.s_addr),he->h_addr_list[0],sizeof(in_addr_t));
			sock=socket(PF_INET,SOCK_STREAM,0);
			if(sock==-1){
				xerror("socket");
				goto ERR_2;
			}
			if(connect(sock,(struct sockaddr*)(&sa),sizeof(sa))==-1){
				xerror("connect");
				close(sock);
				goto ERR_2;
			}
			rc=cm_aln2pos_remote(sock,aln_part_fn,pos_part_fn,&pos_part_n);
			close(sock);
			if(rc)
				goto ERR_2;
		}else{
			rc=cm_aln2pos_core(sadb,aln_part_fn,pos_part_fn,&pos_part_n);
			if(rc)
				goto ERR_2;
		}
		fprintf(f_pos_meta,"%u\t%llu\n",pos_part_no,pos_part_n);
		t_end=get_run_time();
		fprintf(stderr,"[aln2pos] Finished %s in %.2f sec\n", 
			aln_part_fn,t_end-t_pstart);
	}

	t_end=get_run_time();
	fprintf(stderr,"[aln2pos] All part processed in %.2f sec\n",t_end-t_start);
ERR_2:
	fclose(f_pos_meta);
	fclose(f_aln_meta);

	free(aln_part_fn);
	free(pos_part_fn);

	if(sadb)
		sadb_close(sadb);
ERR_0:
	free(pos_prefix);
	return rc;
}

#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <sys/times.h>
#include <fcntl.h>
#include <poll.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>

#include "alndb.h"
#include "posdb.h"
#include "sadb.h"
#include "utils.h"

void usage(void)
{
	fprintf(stderr,
	"Usage:  cm_sadbd [options] <prefix>\n"
	"Options:       -p PORT  service port [9990]\n"
	);
}

int sadbd_core(sadb_t sadb, int sock)
{
	aln_t *aln_buf=malloc(sizeof(aln_t));
	pos_t *pos_buf=NULL;
	int pos_buf_sz=0;

	while(1){
		int i, rc, pos_n;
		pos_t *ppos;
		rc=read(sock,aln_buf,sizeof(aln_t));
		if(rc!=sizeof(aln_t)){
			break; // closed
		}
		pos_n=aln_buf->width;
		if(pos_buf_sz<aln_buf->width){
			pos_buf=realloc(pos_buf,sizeof(pos_t)*pos_n);
			pos_buf_sz=pos_n;
		}
		for(i=0;i<pos_n;++i){
			ppos=pos_buf+i;
			ppos->id=aln_buf->id;
			ppos->r=aln_buf->r;
			ppos->a=aln_buf->a;
			ppos->n_mm=aln_buf->n_mm;
			ppos->n_gapo=aln_buf->n_gapo;
			ppos->n_gape=aln_buf->n_gape;
			ppos->pos=sadb_get(sadb,aln_buf->r,aln_buf->sa+i);
		}
		rc=write(sock,pos_buf,sizeof(pos_t)*pos_n);
		if(rc!=sizeof(pos_t)*pos_n){
			break;
		}
	};
	free(aln_buf);
	free(pos_buf);
	return 0;
}

int main(int argc, char **argv)
{
	int next_opt;
	int sock_s;
	struct sockaddr_in sa;
	socklen_t sa_len;
	short port=9990;
	sadb_t sadb;
	unsigned long long mem_sz=sysconf(_SC_PHYS_PAGES)*sysconf(_SC_PAGE_SIZE);
	const char *prefix;
	unsigned int sa_intv_bits;

	do{
		next_opt=getopt(argc,argv,"p:");
		switch(next_opt){
			case 'p': port=atoi(optarg); break;
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


	prefix=argv[optind];
	for(sa_intv_bits=0;(1llu<<(sa_intv_bits+3))<=mem_sz;sa_intv_bits++);
	sa_intv_bits--;

	if(fork()){
		return 0;
	}

	sadb=sadb_open(prefix,sa_intv_bits|SADB_IN_MEM);
	if(sadb==NULL){
		return -1;
	}

	/* invoke prefetch */
	sadb_get(sadb,0,0);
	sadb_get(sadb,1,0);

	memset(&sa,0,sizeof(sa));
	sa.sin_family=AF_INET;
	sa.sin_port=htons(port);
	sa.sin_addr.s_addr=INADDR_ANY;
	sock_s=socket(PF_INET,SOCK_STREAM,0);
	if(sock_s==-1){
		xerror("socket");
		return -1;
	}
	if(bind(sock_s,(struct sockaddr*)&sa,sizeof(sa))){
		xerror("bind");
		return -1;
	}
	if(listen(sock_s,1)){
		xerror("listen");
		return -1;
	}
	while(1){
		int sock_c;
		sock_c=accept(sock_s,(struct sockaddr*)&sa,&sa_len);
		if(sock_c==-1){
			xerror("accept");
			break;
		}
		/*
		if(fcntl(sock_c,F_SETFL,O_NONBLOCK)==-1){
			xerror("fcntl");
			break;
		}
		*/
		sadbd_core(sadb, sock_c);
		close(sock_c);
	}
	sadb_close(sadb);
	close(sock_s);
	return 0;
}

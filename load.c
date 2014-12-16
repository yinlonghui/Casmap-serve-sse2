#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netdb.h>

#include "load.h"
#include "serve.h"
#include "aln.h"
#include "utils.h"

struct load_rt{
	struct load_opt *opt;
	int fdev[MAX_DEV_NUM];
	unsigned int dev_num;
	int fd_bwt[2];
};

static void dump_db_status(struct db_param* db_param)
{
	fprintf(stderr,"[load] CNT[0-4]:\t%u %u %u %u\n",
			db_param->CNT[0],db_param->CNT[1],db_param->CNT[2],db_param->CNT[3]);
	fprintf(stderr,"[load] PRIMARY:\t%u\n",db_param->PRIMARY);
	fprintf(stderr,"[load] RPRIMARY:\t%u\n",db_param->RPRIMARY);
	fprintf(stderr,"[load] BWT_LEN:\t%u\n",db_param->BWT_LEN);
}

bwt_t* create_bwt(int fd)
{
	int i;
	uint32_t tmp;
	bwt_t* bwt;
	struct stat stat;
	bwt=(bwt_t*)malloc(sizeof(bwt_t));
	bwt->fd=fd;

	/* retrieve header info */
	lseek(fd,0,SEEK_SET);
	if(read(fd,&tmp,sizeof(tmp))==-1)
		goto ERR_OUT;
	bwt->primary=le32toh(tmp);
	bwt->cnt[0] = 0;
	for(i=1;i<4;i++)
	{
		if(read(fd,&tmp,sizeof(tmp))==-1)
			goto ERR_OUT;
		bwt->cnt[i] = le32toh(tmp);
	}
	if(read(fd,&tmp,sizeof(tmp))==-1)
		goto ERR_OUT;
	bwt->length = le32toh(tmp);
	if(fstat(fd,&stat)==-1)
		goto ERR_OUT;
	bwt->size = (stat.st_size-4*5)/32;

	lseek(fd,5*4,SEEK_SET);
	return bwt;
ERR_OUT:
	free(bwt);
	return NULL;
}

void destroy_bwt(bwt_t *bwt)
{
	if(bwt)
	{
		free(bwt);
	}
}

static inline int bwt_address(unsigned int i, int a)
{
	//fprintf(stderr,"BWT addr: %x %x %x %x\n",
	unsigned int rank, bank, row, col;

	rank=(a&0x1)<<28;
	bank=(i&0x1c0)<<19;
	row=(i&0xfffe0000)>>7;
	col=(i&0x1fe00)>>7;
	return rank|bank|row|col;
}

static int write_package(struct load_rt *rt, int type, void *data, int size)
{
	int i,rc;
	request_t req_hdr;
	req_hdr=REQUEST(type,size);
	htobe32s(&req_hdr,sizeof(req_hdr));

	for(i=0;i<rt->dev_num;++i){
		rc=xwrite(rt->fdev[i],&req_hdr,sizeof(req_hdr));
		if(rc!=sizeof(req_hdr)){
			xerror("write request");
			return -1;
		}
		if(rt->opt->debug){
			fprintf(stderr,"WRITE %d: ",rc);
			print_memory(stderr,&req_hdr,rc);
		}

		if(data && size){
			rc=xwrite(rt->fdev[i],data,size);
			if(rc!=size){
				xerror("write parameter");
				return -1;
			}
			if(rt->opt->debug){
				fprintf(stderr,"WRITE %d: ",rc);
				print_memory(stderr,data,rc);
			}
		}
	}
	return 0;
}

static int load_db_parameters(struct load_rt *rt, struct db_param *db_param){
	struct param_package buffer[32];
	struct param_package *pkg=buffer;
	unsigned int w_size;
	int i;
	int rc;
	for(i=0;i<32;++i){
		init_param_header(&buffer[i].header);
	}

	pkg->id = BWA_PARAM_UPDATE_EN_ID;
	pkg->value = 1;
	++pkg;

	pkg->id = BWA_PARAM_CNT0_ID;
	pkg->value = db_param->CNT[0];
	++pkg;

	pkg->id = BWA_PARAM_CNT1_ID;
	pkg->value = db_param->CNT[1];
	++pkg;

	pkg->id = BWA_PARAM_CNT2_ID;
	pkg->value = db_param->CNT[2];
	++pkg;

	pkg->id = BWA_PARAM_CNT3_ID;
	pkg->value = db_param->CNT[3];
	++pkg;

	pkg->id = BWA_PARAM_BWT_LEN_ID;
	pkg->value = db_param->BWT_LEN;
	++pkg;

	pkg->id = BWA_PARAM_PRIMARY_ID;
	pkg->value = db_param->PRIMARY;
	++pkg;

	pkg->id = BWA_PARAM_RPRIMARY_ID;
	pkg->value = db_param->RPRIMARY;
	++pkg;

	w_size=(pkg-buffer)*sizeof(struct param_package);
	htobe32s(buffer,w_size);

	rc=write_package(rt,REQUEST_BROADCAST,buffer,w_size);
	return rc;
}

static int load_training_parameters(struct load_rt *rt){
	struct param_package buffer[32];
	request_t req_hdr;
	struct param_package *pkg=buffer;
	unsigned int w_size;
	int i;
	int rc;
	for(i=0;i<32;++i){
		init_param_header(&buffer[i].header);
	}

	pkg->id = BWA_PARAM_UPDATE_EN_ID;
	pkg->value = 0;
	++pkg;

	pkg->id = BWA_PARAM_TRAINING_ID;
	pkg->value = 1;
	++pkg;

	w_size=(pkg-buffer)*sizeof(struct param_package);
	req_hdr=REQUEST(REQUEST_BROADCAST,w_size);

	htobe32s(buffer,w_size);
	htobe32s(&req_hdr,sizeof(req_hdr));

	for(i=0;i<rt->dev_num;++i){
		rc=xwrite(rt->fdev[i],&req_hdr,sizeof(req_hdr));
		if(rc!=sizeof(req_hdr)){
			xerror("write request");
			return -1;
		}
		if(rt->opt->debug){
			fprintf(stderr,"WRITE %d: ",rc);
			print_memory(stderr,&req_hdr,rc);
		}

		rc=xwrite(rt->fdev[i],buffer,w_size);
		if(rc!=w_size){
			xerror("write parameter");
			return -1;
		}
		if(rt->opt->debug){
			fprintf(stderr,"WRITE %d: ",rc);
			print_memory(stderr,buffer,rc);
		}
	}
	return 0;
}


static int load_bwt(struct load_rt *rt, bwt_t *bwt, int a) {
	int rc;
	int i,j,n;
	uint32_t *buf;
	int chk_size;
	struct db_package *db_pkg;
	int progress=0;
	int w_size;

	buf=malloc(sizeof(uint32_t)*8*rt->opt->chk_size);
	db_pkg=malloc(sizeof(struct db_package)*rt->opt->chk_size);

	fprintf(stderr,"  0%%");
	
	n=0;
	while(n<bwt->size){
		chk_size=rt->opt->chk_size;
		if(n+chk_size>bwt->size)
			chk_size=bwt->size-n;
		if(read(bwt->fd,buf,sizeof(uint32_t)*8*chk_size)!=
				sizeof(uint32_t)*8*chk_size){
			xerror("Read bwt file");
			goto ERR_0;
		}

		le32tohs(buf,sizeof(uint32_t)*8*chk_size);

		for(i=0;i<chk_size;++i){

			init_db_header(&db_pkg[i].header);

			db_pkg[i].address = bwt_address((n+i)*64,a);
			for(j=0;j<8;++j){
				db_pkg[i].data[j] = buf[i*8+j];
			}

			/*
			if(rt->opt->debug)
				print_package(&db_pkg[i]);
			*/
		}

		w_size=sizeof(struct db_package)*chk_size;
		htobe32s(db_pkg,w_size);

		rc=write_package(rt,REQUEST_BROADCAST,db_pkg,w_size);
		if(rc){
			return -1;
		}

		n+=chk_size;
		if(n/(bwt->size/100)!=progress){
			progress=n/(bwt->size/100);
			fprintf(stderr,"\b\b\b\b%3u%%",n/(bwt->size/100));
		}
	}

	init_db_header(&db_pkg->header);
	db_pkg->address = bwt_address(-1,a);
	db_pkg->data[0]=0xffffffc0;
	for(j=1;j<8;j++) db_pkg->data[j]=0;

	/*
	if(rt->opt->debug)
		print_package(db_pkg);
	*/

	w_size=sizeof(struct db_package);
	htobe32s(db_pkg,w_size);

	rc=write_package(rt,REQUEST_BROADCAST,db_pkg,w_size);
	if(rc){
		return -1;
	}

	fprintf(stderr,"\b\b\b\b100%%");
	free(db_pkg);
	free(buf);
	return 0;
ERR_0:
	free(db_pkg);
	free(buf);
	return -1;
}

int cm_load_core(struct load_rt *rt)
{
	int rc;
	bwt_t *bwt,*rbwt;
	struct db_param db_param;
	memset(&db_param,0,sizeof(db_param));

	bwt=create_bwt(rt->fd_bwt[0]);
	rbwt=create_bwt(rt->fd_bwt[1]);

	db_param.BWT_LEN = bwt->length;
	db_param.PRIMARY = bwt->primary;
	db_param.RPRIMARY = rbwt->primary;
	db_param.CNT[0] = bwt->cnt[0];
	db_param.CNT[1] = bwt->cnt[1];
	db_param.CNT[2] = bwt->cnt[2];
	db_param.CNT[3] = bwt->cnt[3];

	dump_db_status(&db_param);

	fprintf(stderr,"[load] Loading database parameters...");
	rc= load_db_parameters(rt,&db_param);
	if(rc){
		fprintf(stderr,"FAILED\n");
		goto clean_out;
	}
	fprintf(stderr,"OK\n");

	fprintf(stderr,"[load] Loading forward database...");
	rc= load_bwt(rt, bwt, 0);
	if(rc){
		fprintf(stderr,"FAILED\n");
		goto clean_out;
	}
	fprintf(stderr," OK\n");

	fprintf(stderr,"[load] Loading reversed database...");
	rc=load_bwt(rt, rbwt, 1);
	if(rc){
		fprintf(stderr,"FAILED\n");
		goto clean_out;
	}
	fprintf(stderr," OK\n");

/*
	fprintf(stderr,"Loading database parameters...");
	rc= load_db_parameters(rt,&db_param);
	if(rc){
		fprintf(stderr,"FAILED\n");
		goto clean_out;
	}
	fprintf(stderr,"OK\n");
*/

	/*
	rc=write_package(rt,REQUEST_END,NULL,0);
	if(rc){
		return -1;
	}
	*/

clean_out:
	destroy_bwt(bwt);
	destroy_bwt(rbwt);
	return rc;
}

int cm_training_core(struct load_rt *rt)
{
	int rc,i,j;
	unsigned char *buffer,*discard;
	int pkg_size=sizeof(struct reads_header)+4;
	struct reads_header *rpkg;
	uint32_t *dptr;

	buffer=malloc(pkg_size);
	discard=malloc(1024);

	rpkg=(struct reads_header*)buffer;
	dptr=(uint32_t*)(buffer+sizeof(struct reads_header));

	for(i=0;i<rt->dev_num;++i){
		if(fcntl(rt->fdev[i],F_SETFL,O_NONBLOCK)<0){
			xerror("Set non-block");
			return -1;
		}
	}

	fprintf(stderr,"[load] Loading training parameters...");
	rc= load_training_parameters(rt);
	if(rc){
		fprintf(stderr,"FAILED\n");
		return -1;
	}
	fprintf(stderr,"OK\n");

	fprintf(stderr,"[load] Cache training...");
	init_read_header(&rpkg->header,
		(sizeof(struct reads_header)-sizeof(struct package_header))/4+1);
	rpkg->gid = GID_INVALID;
	rpkg->max_diff = 0;
	rpkg->length = 4;
	rpkg->nen = 1;
	rpkg->nmode = BWA_EXT_STATE_WIDTH;
	rpkg->neffort = BWA_ONE_HIT;
	rpkg->ren = 1;
	rpkg->rmode = BWA_EXT_STATE_WIDTH;
	rpkg->reffort = BWA_ONE_HIT;
	rpkg->seed_diff = 0;
	htobe32s(rpkg,sizeof(struct reads_header));
	
	for(j=0;j<256;++j){
		unsigned int a=j/64;
		unsigned int b=(j%64)/16;
		unsigned int c=(j%16)/4;
		unsigned int d=j%4;
		*dptr=htobe32((a<<28)|(b<<24)|(c<<20)|(d<<16));

		rc=write_package(rt,REQUEST_BROADCAST,buffer,pkg_size);
		if(rc){
			fprintf(stderr,"FAILED\n");
			return -1;
		}
		/* discard match data for it has no use */
		/* this is also important for avoiding deadlock */
		for(i=0;i<rt->dev_num;++i){
			//read(rt->fdev[i],discard,1024);
		}
	}

	fprintf(stderr,"OK\n");

	free(buffer);
	free(discard);
	return 0;
}

int cm_load(struct load_opt *opt)
{
	int i,rc;
	char *fn_buf;
	struct load_rt rt;
	struct sockaddr_in sa;
	struct hostent *he;

	memset(&rt,0,sizeof(rt));
	rt.opt=opt;

	he=gethostbyname(opt->server);
	if(he==NULL){
		herror("gethostbyname\n");
		return -1;
	}

	memset(&sa,0,sizeof(sa));
	sa.sin_family=AF_INET;
	sa.sin_port=htons(opt->port);
	memcpy(&(sa.sin_addr.s_addr),he->h_addr_list[0],sizeof(in_addr_t));
	fprintf(stderr,"[load] Connecting to %s:%u\n",inet_ntoa(sa.sin_addr),opt->port);
	rt.fdev[0]=socket(PF_INET,SOCK_STREAM,0);
	if(rt.fdev[0]==-1){
		xerror("socket");
		return -1;
	}
	if(connect(rt.fdev[0],(struct sockaddr*)(&sa),sizeof(sa))==-1){
		xerror("connect");
		close(rt.fdev[0]);
		return -1;
	}
	fprintf(stderr,"[load] Connection established. Start alignment...\n");
	rt.dev_num=1;

	fn_buf=malloc(strlen(opt->prefix)+16);

	sprintf(fn_buf,"%s%s",opt->prefix,".bwt");
	fprintf(stderr,"[load] BWT file:\t%s\n",fn_buf);
	rt.fd_bwt[0]=open(fn_buf,O_RDONLY);
	if(rt.fd_bwt[0]==-1){
		xerror("Open BWT");
		return -1;
	}

	sprintf(fn_buf,"%s%s",opt->prefix,".rbwt");
	fprintf(stderr,"[load] RBWT file:\t%s\n",fn_buf);
	rt.fd_bwt[1]=open(fn_buf,O_RDONLY);
	if(rt.fd_bwt[1]==-1){
		xerror("Open RBWT");
		return -1;
	}

	free(fn_buf);

	rc=cm_load_core(&rt);
	if(rc){
		return -1;
	}
	
	rc=cm_training_core(&rt);
	if(rc){
		return -1;
	}

	close(rt.fd_bwt[0]);
	close(rt.fd_bwt[1]);

	/*
	write_package(&rt,REQUEST_END,NULL,0);
	*/

	for(i=0;i<rt.dev_num;++i){
		close(rt.fdev[i]);
	}

	return 0;
}

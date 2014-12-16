#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <unistd.h>
#include <sys/times.h>
#include <fcntl.h>
#include <errno.h>
#include <poll.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <netdb.h>

#include "utils.h"
#include "seq.h"
#include "seqdb.h"
#include "alndb.h"
#include "aln.h"
#include "serve.h"

#define RPKG_BUFFER_SIZE 	(8*1024)
#define MPKG_BUFFER_SIZE	(sizeof(struct match_package)*1024)
#define SEQ_PREFETCH_NUM	(1024)
#define ITEM_TABLE_INIT_SIZE	(1024)

struct aln_item
{
	int id;
	seq_t *seq;
	unsigned int state;
	hit_t *hits_list;
	unsigned int n_hit;
	struct aln_item *prev;
	struct aln_item *next;
};

struct aln_rt{
	struct aln_opt *opt;
	int fdev[MAX_DEV_NUM];
	int dev_num;
	seqdb_t seqdb;
	alndb_t alndb;
	hit_t *hit_free;
	struct aln_item *item_free;
	struct aln_item *item_pending;
	struct aln_item *item_running;
	struct aln_item **item_table;
	unsigned int item_table_sz;
	unsigned int item_num;
	unsigned int seq_p_sz,seq_r_sz;
	unsigned long long seq_total;
	unsigned long long map_total;
	unsigned long long aln_total;
	int rep_enable;
	double rep_last_time;
	unsigned long long rep_total;
	unsigned long rep_max_speed;
};


static inline void destroy_item(struct aln_item *item)
{
	hit_t *tmp;
	if(item) {
		if(item->hits_list){
			while(!list_empty(item->hits_list)){
				tmp = item->hits_list->next;
				list_del(tmp);
				free(tmp);
			}
			free(item->hits_list);
		}
		free(item);
	}
}

static void free_item_list(struct aln_item *header)
{
	struct aln_item *tmp, *list;
	if(header){
		list=header;
		while(!list_empty(list))
		{
			tmp=list->next;
			list_del(tmp);
			destroy_item(tmp);
		}
		free(list);
	}
}

static void init_hit_list(struct aln_rt *rt){
	hit_t *hit_list;
	hit_list=(hit_t*)malloc(sizeof(hit_t));
	list_init_header(hit_list);
	rt->hit_free=hit_list;
}

static void free_hit_list(struct aln_rt *rt){
	hit_t *hit_list,*tmp;
	hit_list=rt->hit_free;
	while(!list_empty(hit_list)){
		tmp = hit_list->next;
		list_del(tmp);
		free(tmp);
	}
	free(hit_list);
	rt->hit_free=NULL;
}

static void init_item_lists(struct aln_rt *rt){
	struct aln_item *item_list;
	item_list=(struct aln_item*)malloc(sizeof(struct aln_item));
	list_init_header(item_list);
	rt->item_free=item_list;

	item_list=(struct aln_item*)malloc(sizeof(struct aln_item));
	list_init_header(item_list);
	rt->item_pending=item_list;

	item_list=(struct aln_item*)malloc(sizeof(struct aln_item));
	list_init_header(item_list);
	rt->item_running=item_list;

	init_hit_list(rt);
}

static void free_item_lists(struct aln_rt *rt){
	free_item_list(rt->item_running);
	rt->item_running=NULL;
	free_item_list(rt->item_pending);
	rt->item_pending=NULL;
	free_item_list(rt->item_free);
	rt->item_free=NULL;

	free_hit_list(rt);

	rt->item_num=0;
	free(rt->item_table);
	rt->item_table=NULL;
	rt->item_table_sz=0;
}

static inline hit_t* alloc_hit(struct aln_rt *rt){
	hit_t *hit_list,*hit;
	hit_list=rt->hit_free;
	if(list_empty(hit_list)){
		hit=(hit_t*)malloc(sizeof(hit_t));
	}
	else{
		hit=hit_list->next;
		list_del(hit);
	}
	return hit;
}

static inline void release_hit(hit_t *hit, struct aln_rt *rt){
	list_add_tail(rt->hit_free, hit);
}

static inline struct aln_item * alloc_item(struct aln_rt *rt)
{
	struct aln_item *item;
	if(list_empty(rt->item_free)){
		item=(struct aln_item*)malloc(sizeof(struct aln_item));
		memset(item,0,sizeof(struct aln_item));

		item->id=rt->item_num;
		rt->item_num++;
		if(rt->item_num>rt->item_table_sz){
			if(rt->item_table_sz==0)
				rt->item_table_sz=ITEM_TABLE_INIT_SIZE;
			else
				rt->item_table_sz*=2;
			rt->item_table=realloc(rt->item_table,
				sizeof(struct aln_item*)*rt->item_table_sz);
		}
		rt->item_table[item->id]=item;

		item->hits_list=(hit_t*)malloc(sizeof(hit_t));
		list_init_header(item->hits_list);
	}
	else{
		item=rt->item_free->next;
		list_del(item);
	}
	item->state = SEQ_INIT;
	item->n_hit = 0;
	return item;
};

static inline void release_item(struct aln_item *item, struct aln_rt *rt){
	hit_t *tmp;
	while(!list_empty(item->hits_list)){
		tmp=item->hits_list->next;
		list_del(tmp);
		release_hit(tmp,rt);
	}
	list_add_tail(rt->item_free, item);
}

static int load_aln_parameters(struct aln_rt *rt){
	struct param_package buffer[32];
	request_t req_hdr;
	struct aln_param *aln_param=&rt->opt->aln_param;
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

	pkg->id = BWA_PARAM_MAX_GAPO_ID;
	pkg->value = aln_param->MAX_GAPO;
	++pkg;

	pkg->id = BWA_PARAM_MAX_GAPE_ID;
	pkg->value = aln_param->MAX_GAPE;
	++pkg;

	pkg->id = BWA_PARAM_MAX_DEL_OCC;
	pkg->value = aln_param->MAX_DEL_OCC;
	++pkg;

	pkg->id = BWA_PARAM_S_MM_ID;
	pkg->value = aln_param->S_MM;
	++pkg;

	pkg->id = BWA_PARAM_S_GAPO_ID;
	pkg->value = aln_param->S_GAPO;
	++pkg;

	pkg->id = BWA_PARAM_S_GAPE_ID;
	pkg->value = aln_param->S_GAPE;
	++pkg;

	pkg->id = BWA_PARAM_SEED_WIDTH_ID;
	pkg->value = aln_param->SEED_WIDTH;
	++pkg;

	pkg->id = BWA_PARAM_INDEL_HEAD_SKIP_ID;
	pkg->value = aln_param->INDEL_END_SKIP>2?aln_param->INDEL_END_SKIP-2:0;
	++pkg;

	pkg->id = BWA_PARAM_INDEL_TAIL_SKIP_ID;
	pkg->value = aln_param->INDEL_END_SKIP;
	++pkg;

	pkg->id = BWA_PARAM_FAST_MODE_ID;
	pkg->value = aln_param->FAST_MODE;
	++pkg;

	pkg->id = BWA_PARAM_TRAINING_ID;
	pkg->value = 0;
	++pkg;

	pkg->id = BWA_PARAM_RID_WIDTH_ID;
	pkg->value =  aln_param->RID_WIDTH;
	++pkg;

	pkg->id = BWA_PARAM_MAX_BACK_ID;
	pkg->value =  aln_param->MAX_BACK;
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
			fprintf(stderr,"WRITE %d: \n",rc);
			print_memory(stderr,&req_hdr,rc);
		}

		rc=xwrite(rt->fdev[i],buffer,w_size);
		if(rc!=w_size){
			xerror("write parameter");
			return -1;
		}
		if(rt->opt->debug){
			fprintf(stderr,"WRITE %d: \n",rc);
			print_memory(stderr,buffer,rc);
		}
	}
	return 0;
}

static unsigned int prepare_read_package(void *buf, struct aln_item *item, struct aln_rt *rt){
	struct reads_header *rpkg = (struct reads_header*)buf;
	int mm_size;
	int plen=item->seq->length/8;
	if(item->seq->length%8)
		++plen;
	mm_size=sizeof(struct reads_header)+plen*4;
	init_read_header(&rpkg->header,
			(mm_size-sizeof(struct package_header))/4);
	rpkg->gid = item->id;
	rpkg->max_diff = cal_max_diff(item->seq->length);
	rpkg->length = item->seq->length;
	rpkg->nen = 1;
	if(rt->opt->aln_param.SEED_WIDTH > 0 
			&& item->seq->length > rt->opt->aln_param.SEED_WIDTH)
		rpkg->nmode = BWA_EXT_STATE_SEED;
	else
		rpkg->nmode = BWA_EXT_STATE_WIDTH;
	rpkg->neffort = rt->opt->aln_param.EFFORT;
	rpkg->ren = 1;
	rpkg->rmode = rpkg->nmode;
	rpkg->reffort = rpkg->neffort;
	rpkg->seed_diff = rt->opt->aln_param.MAX_SEED_DIFF;
	pack_seq((uint32_t*)((void*)rpkg+sizeof(struct reads_header)),item->seq->seq,item->seq->length);
	if(rt->opt->debug){;
		print_package(buf);
	}
	return mm_size;
}

static void cm_dump_aln_parameters(struct aln_param *param)
{
	fprintf(stderr,"[seq2aln] MAX_GAPO: %u\n",param->MAX_GAPO);
	fprintf(stderr,"[seq2aln] MAX_GAPE: %u\n",param->MAX_GAPE);
	fprintf(stderr,"[seq2aln] SEED_WIDTH: %u\n",param->SEED_WIDTH);
	fprintf(stderr,"[seq2aln] MAX_SEED_DIFF: %u\n",param->MAX_SEED_DIFF);
	fprintf(stderr,"[seq2aln] S_MM: %u\n",param->S_MM);
	fprintf(stderr,"[seq2aln] S_GAPO: %u\n",param->S_GAPO);
	fprintf(stderr,"[seq2aln] S_GAPE: %u\n",param->S_GAPE);
	fprintf(stderr,"[seq2aln] MAX_DEL_OCC: %u\n",param->MAX_DEL_OCC);
	fprintf(stderr,"[seq2aln] INDEL_END_SKIP: %u\n",param->INDEL_END_SKIP);
	fprintf(stderr,"[seq2aln] EFFORT: %u\n",param->EFFORT);
	//fprintf(stderr,"[seq2aln] FAST_MODE: %u\n",param->FAST_MODE);
	//fprintf(stderr,"[seq2aln] RID_WIDTH: %u\n",param->RID_WIDTH);
	fprintf(stderr,"[seq2aln] MAX_BACK: %u\n",param->MAX_BACK);
}

static int prefetch_seq(unsigned int num, struct aln_rt *rt)
{
	unsigned int i;
	struct aln_item *item;
	seq_t *seq;
	for(i=0;i<num;++i){
		item=alloc_item(rt);
		seq=seqdb_get(rt->seqdb, SEQDB_NEXT);
		if(seq==NULL){
			release_item(item,rt);
			break;
		}
		list_add_tail(rt->item_pending,item);
		++(rt->seq_p_sz);
		item->seq=seq;
		if(rt->opt->debug){
			fprintf(stderr,"SEQ:%lu %s %s %s\n",
				seq->uid,seq->name,seq->seq,seq->qual); 
		}
	}
	return i;
}

static inline int seq_finish_cleanup(struct aln_rt *rt)
{
	struct aln_item *item;
	hit_t *hit;
	// output all finished seq at the front of the list
	while(!list_empty(rt->item_running)){
		item=rt->item_running->next;

		if(item->state!=SEQ_FINISHED)
			break;

		if(rt->opt->debug){
			fprintf(stderr, "OUT: uid %lu item @%p n_hit %u\n",
					item->seq->uid, item, item->n_hit);
		}

		//fprintf(stderr,"ALN: %llu %u\n", item->seq->uid, item->n_hit);
		for(hit=item->hits_list->next;hit!=item->hits_list;hit=hit->next){
			if(alndb_put(rt->alndb,hit,item->seq->uid)){
				return -1;
			}
			rt->aln_total++;
		}

		// release seq
		seqdb_release(rt->seqdb,item->seq);
		list_del(item);
		release_item(item,rt);

		--(rt->seq_r_sz);
	}
	return 0;
}

static void print_mpkg(struct match_package *mpkg)
{
	fprintf(stderr,
		"MATCH: gid %u a %u hit %u last %u ext_state %u "
		"n_mm %u n_gapo %u n_gape %u k %u l %u\n",
		mpkg->gid,mpkg->a,mpkg->hit,mpkg->last,mpkg->ext_state,
		mpkg->n_mm,mpkg->n_gapo,mpkg->n_gape,mpkg->k,mpkg->l);
}

static inline int process_match(struct match_package *mpkg, struct aln_rt *rt)
{
	struct aln_item *item;

	if(0){
		print_mpkg(mpkg);
	}

	/* discard data */
	if(mpkg->gid==GID_INVALID)
		return 0;

	if(mpkg->gid>rt->item_num){
		fprintf(stderr, 
			"incorrect match package gid, data corrupted\n");
		print_mpkg(mpkg);
		return -1;
	}

	item=rt->item_table[mpkg->gid];

	if(!(item->state&SEQ_NORMAL_LOADED && item->state&SEQ_REVERSE_LOADED)){
		fprintf(stderr, 
			"incorrect item state (uid %lu, item @%p, state %x), "
			"data corrupted\n",
			item->seq->uid, item, item->state);
		print_mpkg(mpkg);
		return -1;
	}	

	if(!mpkg->hit && !mpkg->last){
		fprintf(stderr, 
			"incorrect match package format, data corrupted\n");
		print_mpkg(mpkg);
		return -1;
	}

	// correct sai at 0
	if(mpkg->hit && mpkg->k==0){
		if(mpkg->l>0)
			mpkg->k=1;
		else // drop sai [0,0]
			mpkg->hit=0;
	}

	// check hit
	if(mpkg->hit){
		hit_t *c_hit, *p_hit;
	
		c_hit=alloc_hit(rt);
		c_hit->r=c_hit->a=mpkg->a;
		if(mpkg->ext_state==BWA_EXT_STATE_MATCH)
			c_hit->a=!c_hit->a;
		c_hit->n_mm = mpkg->n_mm;
		c_hit->n_gapo = mpkg->n_gapo;
		c_hit->n_gape = mpkg->n_gape;
		c_hit->k = mpkg->k;
		c_hit->l = mpkg->l;
		c_hit->score=(rt->opt->aln_param.S_MM)*(c_hit->n_mm)+
			(rt->opt->aln_param.S_GAPO)*(c_hit->n_gapo)+
			(rt->opt->aln_param.S_GAPE)*(c_hit->n_gape);

		// sorted insert
		p_hit=item->hits_list->next;
		while(1){
			if(p_hit==item->hits_list || p_hit->score > c_hit->score){
				list_add_tail(p_hit,c_hit);
				item->n_hit++;
				p_hit=c_hit->next;
				break;
			}else if(p_hit->k==c_hit->k && p_hit->l==c_hit->l && 
					p_hit->a==c_hit->a && p_hit->r==c_hit->r){
				p_hit=item->hits_list;
				release_hit(c_hit,rt);
				break;
			}
			p_hit=p_hit->next;
		}
		// remove duplicate SAI with worse score.
		while(p_hit!=item->hits_list){
			hit_t* tmp=p_hit;
			p_hit=p_hit->next;
			if(tmp->k==c_hit->k && tmp->l==c_hit->l && 
				tmp->a==c_hit->a && tmp->r==c_hit->r){
				list_del(tmp);
				release_hit(tmp,rt);
				item->n_hit--;
			}
		}
	}

	// check last
	if(mpkg->last){

		/*
		if(mpkg->hit && mpkg->ext_state==BWA_EXT_STATE_MATCH) {
			fprintf(stderr,"Last and hit: uid %lu (%u,%u)\n",
				item->seq->uid, mpkg->k,mpkg->l);
		}
		*/

		if(0){
			fprintf(stderr, 
				"LAST: uid %lu item @%p\n",
				item->seq->uid, item);
		}

		if(item->state & SEQ_NORMAL_DONE){
			if(item->state & SEQ_REVERSE_DONE){
				fprintf(stderr, 
						"incorrect item state (uid %lu, item @%p, state %x), "
						"data corrupted\n",
						item->seq->uid, item, item->state);
				print_mpkg(mpkg);
			}else{
				item->state |= SEQ_REVERSE_DONE;
			}
		}else{
			item->state |= SEQ_NORMAL_DONE;
		}

		if(item->state == SEQ_FINISHED){
			++(rt->seq_total);
			if(item->n_hit)
				++(rt->map_total);
			if(seq_finish_cleanup(rt)){
				return -1;
			}
		}
	}
	return 0;
}

static inline void progress_report(struct aln_rt *rt,int intval)
{
	double time_curr=get_run_time();
	if(time_curr-rt->rep_last_time >intval){
		double time_interval = time_curr-rt->rep_last_time;
		unsigned long seq_interval = rt->seq_total-rt->rep_total;
		unsigned long curr_speed = seq_interval/time_interval;
		unsigned long percent = seqdb_offset(rt->seqdb)*100/seqdb_size(rt->seqdb);

		if(curr_speed>rt->rep_max_speed)
			rt->rep_max_speed=curr_speed;
		rt->rep_total = rt->seq_total;
		rt->rep_last_time = time_curr;

		fprintf(stderr, "[seq2aln] %.2f sec, %llu reads(%lu%%), %llu aln, %lu r/s\n",
				time_curr, rt->seq_total, percent, rt->aln_total, curr_speed);
	}
}

int cm_aln_core(struct aln_rt *rt)
{
	int i;
	int rc;
	int more_reads=1;
	void* rpkg_buffer[rt->dev_num];
	unsigned int rpkg_size[rt->dev_num], rpkg_offset[rt->dev_num];
	void* mpkg_buffer[rt->dev_num];
	unsigned int mpkg_offset[rt->dev_num],mpkg_start[rt->dev_num];
	struct pollfd fds[rt->dev_num];

	rt->rep_last_time = 0;
	rt->rep_total = 0;
	rt->rep_max_speed = 0;

	init_item_lists(rt);

	rt->seq_p_sz = 0;
	rt->seq_r_sz = 0;
	rt->seq_total = 0;
	rt->aln_total = 0;

	for(i=0;i<rt->dev_num;++i){
		rpkg_buffer[i]=malloc(RPKG_BUFFER_SIZE);
		rpkg_size[i]=0;
		rpkg_offset[i]=0;

		mpkg_buffer[i]=malloc(MPKG_BUFFER_SIZE);
		mpkg_offset[i]=0;
		mpkg_start[i]=0;

		if(fcntl(rt->fdev[i],F_SETFL,O_NONBLOCK)<0){
			xerror("Set non-block");
			return -1;
		}
		fds[i].fd=rt->fdev[i];
	}

	do{
		progress_report(rt,1);
		// stage 1: read as much reads.
		if(more_reads && rt->seq_p_sz < SEQ_PREFETCH_NUM){
			rc=prefetch_seq(SEQ_PREFETCH_NUM,rt);
			if(rc==0)
				more_reads = 0;
		}

		// wait for io condition
		for(i=0;i<rt->dev_num;++i){
			fds[i].events=0;
			if(rpkg_offset[i]!=rpkg_size[i] || rt->seq_p_sz || rt->seq_r_sz)
				fds[i].events|=POLLOUT;

			if(rt->seq_r_sz)
				fds[i].events|=POLLIN;
		}

		rc=poll(fds,rt->dev_num,-1);
		if(rc<0){
			if(errno==EINTR)
				continue;
			xerror("IO select");
			rc = -1;
			goto fall_back;
		}

		for(i=0;i<rt->dev_num;++i){
			// read match result from device
			if(fds[i].revents&POLLIN){
				rc=read(rt->fdev[i],mpkg_buffer[i]+mpkg_offset[i],MPKG_BUFFER_SIZE-mpkg_offset[i]);
				if(rc>0){
					//fprintf(stderr, "read %d\n", rc);
					mpkg_offset[i]+=rc;
					while(mpkg_offset[i]-mpkg_start[i]>=sizeof(struct match_package)){

						be32tohs(mpkg_buffer[i]+mpkg_start[i],sizeof(struct match_package));

						if(rt->opt->debug){
							fprintf(stderr, "%p %p ",mpkg_buffer[i]+mpkg_start[i],mpkg_buffer[i]+mpkg_offset[i]);
							print_package(mpkg_buffer[i]+mpkg_start[i]);
						}

						if(process_match((struct match_package *)(mpkg_buffer[i]+mpkg_start[i]),rt)){
							goto fall_back;

						}
						mpkg_start[i]+=sizeof(struct match_package);
					}
					if(mpkg_offset[i]==MPKG_BUFFER_SIZE){
						mpkg_start[i]=0;
						mpkg_offset[i]=0;
					}
				}else{
					fprintf(stderr,"server disconnect\n");
					goto fall_back;
				}
			}
			
			// send reads to device
			if(fds[i].revents&POLLOUT){
				//fprintf(stderr,"poll out");
				if(rpkg_offset[i]==rpkg_size[i] && rt->seq_p_sz && rt->seq_r_sz<1024*1024){
					struct aln_item *item;
					size_t offset=0;
					void *data=rpkg_buffer[i]+sizeof(request_t);
					while(rt->seq_p_sz && offset<(RPKG_BUFFER_SIZE/2))
					{
						item=rt->item_pending->next;
						offset+=prepare_read_package(data+offset,
								item,rt);


						item->state|=SEQ_NORMAL_LOADED|SEQ_REVERSE_LOADED;
						list_del(item);
						list_add_tail(rt->item_running, item);

						--(rt->seq_p_sz);
						++(rt->seq_r_sz);
					}
					*(request_t*)rpkg_buffer[i]=REQUEST(REQUEST_STREAM,offset);
					rpkg_size[i]=offset+sizeof(request_t);
					rpkg_offset[i]=0;
					htobe32s(rpkg_buffer[i],rpkg_size[i]);
				}
				// write pending data
				size_t w_size=rpkg_size[i]-rpkg_offset[i];
				if(w_size){
					rc=write(rt->fdev[i],rpkg_buffer[i]+rpkg_offset[i],w_size);
					if(rc>0){
						if(rt->opt->debug){
							fprintf(stderr,"WRITE %d: \n",rc);
							print_memory(stderr,rpkg_buffer[i]+rpkg_offset[i],rc);
						}
						rpkg_offset[i]+=rc;
					}
				}
			}
		}
	}while(rt->seq_p_sz || rt->seq_r_sz);

	{
		double time_total=get_run_time();
		progress_report(rt,-1);
		fprintf(stderr,"[seq2aln] Time:%.2f\tTotal Reads: %llu Total Aln: %llu\n",
				time_total,rt->seq_total,rt->aln_total);
		fprintf(stderr,"[seq2aln] Max Speed: %lur/s\tAvg Speed: %lur/s\n",
				rt->rep_max_speed, (unsigned long)(rt->seq_total/time_total));
		fprintf(stderr,"[seq2aln] Mapping Rate: %3.4f%% (%llu/%llu)\n",
				rt->map_total*100.0/rt->seq_total,rt->map_total,rt->seq_total);
	}

	rc = 0;
fall_back:
	for(i=0;i<rt->dev_num;++i){
		free(rpkg_buffer[i]);
		free(mpkg_buffer[i]);
	}
	free_item_lists(rt);
	return rc;
}

int cm_aln(seqdb_t seqdb, alndb_t alndb, struct aln_opt *opt)
{
	int rc;
	int i;
	struct aln_rt rt;
	struct sockaddr_in sa;
	struct hostent *he;

	memset(&rt,0,sizeof(rt));
	rt.opt=opt;

	rt.seqdb=seqdb;
	rt.alndb=alndb;
	rt.dev_num = 0 ;

	for( i = 0 ;  i < 2 ; i++){
		if(opt->server[i])  {
			rt.dev_num++;
			he=gethostbyname(opt->server[i]);
			if(he==NULL){
				herror("gethostbyname\n");
				return -1;
			}

			memset(&sa,0,sizeof(sa));
			sa.sin_family=AF_INET;
			sa.sin_port=htons(opt->port);
			memcpy(&(sa.sin_addr.s_addr),he->h_addr_list[0],sizeof(in_addr_t));
			fprintf(stderr,"[seq2aln] Connecting to %s:%u\n",inet_ntoa(sa.sin_addr),opt->port);
			rt.fdev[i]=socket(PF_INET,SOCK_STREAM,0);
			if(rt.fdev[i]==-1){
				xerror("socket");
				return -1;
			}
			if(connect(rt.fdev[i],(struct sockaddr*)(&sa),sizeof(sa))==-1){
				xerror("connect");
				close(rt.fdev[i]);
				return -1;
			}
			fprintf(stderr,"[seq2aln] Connection established. Start alignment...\n");
		}
	}

	//  open 2 ..

	//fprintf(stderr,"[seq2aln] Input data has %zu bytes\n",seqdb_size(rt.seqdb));

	init_diff_table(opt->fnr);
	if(opt->fnr>=1 || opt->fnr==0)
		fprintf(stderr,"[seq2aln] MAX_DIFF: %d\n",(int)opt->fnr);
	else{
		int i,prev=cal_max_diff(0);
		int next;
		fprintf(stderr,"[seq2aln] FNR: %f\n",opt->fnr);
		for(i=1;i<1024;++i){
			next=cal_max_diff(i);
			if(next!=prev){
				fprintf(stderr,"[seq2aln] MAX_DIFF(< %d dp): %d\n",i,prev);
			}
			prev=next;
		}
		fprintf(stderr,"[seq2aln] MAX_DIFF(< %d dp): %d\n",i,next);
	}

	cm_dump_aln_parameters(&opt->aln_param);

	/* Send param */
	load_aln_parameters(&rt);

	rc=cm_aln_core(&rt);
	
	fprintf(stderr,"[seq2aln] Alignment finished %s\n", rc?"with error":"successful");
	for(i=0;i<rt.dev_num;++i){
		close(rt.fdev[i]);
	}

	return 0;
}

void print_package(void *pkg)
{
	int i;
	struct package_header *hdr=pkg;
	unsigned int dw=hdr->dw+1;
	uint32_t *ptr=pkg;
	switch(hdr->type){
		case PKG_TYPE_PARAM:
			fprintf(stderr,"PPKG");
			break;
		case PKG_TYPE_DB:
			fprintf(stderr,"DPKG");
			break;
		case PKG_TYPE_READ:
			fprintf(stderr,"RPKG");
			break;
		case PKG_TYPE_MATCH:
			fprintf(stderr,"MPKG");
			break;
		default:
			fprintf(stderr,"Illegal PKG @%p: %08x\n",pkg,*ptr);
			return;
	}
	fprintf(stderr,"@%p:",pkg);
	for(i=0;i<dw;++i){
		fprintf(stderr," %08x",*(ptr++));
	}
	fprintf(stderr, "\n");
}




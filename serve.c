#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <errno.h>
#include <fcntl.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/times.h>
#include <signal.h>
#include <stdint.h>
#include <poll.h>
#include <syslog.h>

#include "face.h"
#include "serve.h"
#include "utils.h"

#define MPKG_SIZE	(20)

#define MAX_DEV_NUM	(4)
#define OBUF_SIZE	(MPKG_SIZE*1024)
#define IBUF_SIZE	(1024*1024)

#define STATE_SHUTDOWN	(128)
#define STATE_IDLE	(0)
#define STATE_RUN	(1)
#define STATE_HUP	(8)

struct serve_rt{
	struct serve_opt *opt;
	int fdev[MAX_DEV_NUM];
	unsigned int dev_num;
	int fsock;
	long start_time;
	long last_time;
	long ticks_per_sec;
	int state;
	unsigned long long in_bytes[MAX_DEV_NUM];
	unsigned long long out_bytes[MAX_DEV_NUM];
	int req_type;
	int req_len;
	void *ibuf[MAX_DEV_NUM];
	size_t ibuf_sz[MAX_DEV_NUM];
	size_t ibuf_off[MAX_DEV_NUM];
	void *obuf[MAX_DEV_NUM];
	size_t obuf_sz[MAX_DEV_NUM];
	size_t obuf_off[MAX_DEV_NUM];
	void *bbuf;
	size_t bbuf_off;
	int bflag[MAX_DEV_NUM];
};

static struct serve_rt *glb_rt;

static void debug_memory(void *data,int size)
{
	int t;
	char *buf=malloc(sizeof(char)*size*2+1);
	char *ptr=buf;
	for(t=0;t<size;++t){
		ptr+=sprintf(ptr,"%02x",*(uint8_t*)(data+t));
	}
	syslog(LOG_DEBUG,"%s\n",buf);
	free(buf);
}

static void sig_intr_handler(int sig)
{
	int prev_state=glb_rt->state;
	if(prev_state&STATE_HUP)
		glb_rt->state=STATE_IDLE;
	else if(glb_rt->state&STATE_RUN)
		glb_rt->state|=STATE_HUP;
	else
		glb_rt->state=STATE_SHUTDOWN;
	syslog(LOG_NOTICE,"User Interrupt. STATE[%02x->%02x]\n",prev_state,glb_rt->state);
}

static inline void progress_report(struct serve_rt *rt)
{
	int i,rc;
	static char buf[128];
	uint32_t stat;
	unsigned char t;
	unsigned char vi,va; 
	float elapse_time;
	long curr_time;
	char *str_ptr=buf;
	if(!rt->opt->verbose) 
		return;
	curr_time=times(NULL);
	if(curr_time-rt->last_time>rt->ticks_per_sec) 
		return;
	elapse_time=(curr_time-rt->start_time)*1.0/rt->ticks_per_sec;
	str_ptr+=snprintf(str_ptr,buf-str_ptr,"#  %-8.2f",elapse_time);
	for(i=0;i<rt->dev_num;++i){
		rc=ioctl(rt->fdev[i],FACE_IOCTL_GET_EXT_STAT,&stat);
		if(rc){
			xerror("ioctl");
			return ;
		}
		t=(stat>>8)&0xff;
		vi=(stat>>16)&0xff;
		va=(stat>>24)&0xff;
		str_ptr+=snprintf(str_ptr,buf-str_ptr,"\tD%d:%3.2f:%1.3f:%1.3f",
				i,t*1.97-273.15,vi*0.012,va*0.012);
	}
	snprintf(str_ptr,buf-str_ptr,"\n");
	syslog(LOG_DEBUG,buf);
	rt->last_time=curr_time;
}

static int peek_input(struct serve_rt *rt)
{
	int rc;
	request_t req;	
	rc=xread(rt->fsock,&req,sizeof(req));
	if(rc!=sizeof(req)){
		return -1;
	}
	rt->req_len=REQUEST_LEN(req);
	rt->req_type=REQUEST_TYPE(req);
	if(rt->opt->debug){
		syslog(LOG_DEBUG,"PEEK %08x\n",req);
	}
	return 0;
}

static int push_output(struct serve_rt *rt)
{
	int rc,i;
	for(i=0;i<rt->dev_num;++i){
		if(rt->obuf_off[i]!=0)
			break;
	}
	if(i==rt->dev_num){
		for(i=0;i<rt->dev_num;++i){
			if(rt->obuf_sz[i] && rt->obuf_sz[i]%MPKG_SIZE==0)
				break;
		}
	}
	if(i==rt->dev_num)
		return 0;

	rc=write(rt->fsock,rt->obuf[i]+rt->obuf_off[i],rt->obuf_sz[i]-rt->obuf_off[i]);
	if(rc==-1){
		if(errno==EAGAIN){
			return 0;
		}else{
			xerror("write");
			return -1;
		}
	}

	if(rt->opt->debug){
		syslog(LOG_DEBUG,"WRITE %d: ",rc);
		debug_memory(rt->obuf[i]+rt->obuf_off[i],rc);
	}

	rt->obuf_off[i]+=rc;
	if(rt->obuf_off[i]==rt->obuf_sz[i])
		rt->obuf_off[i]=rt->obuf_sz[i]=0;
	return 0;
}

static int dispatch_data(struct serve_rt *rt)
{
	switch(rt->req_type){
		case REQUEST_STREAM:
			{
				int rc,i;
				for(i=0;i<rt->dev_num;++i){
					if(rt->ibuf_sz[i]==0&&rt->ibuf_off[i]>0)
						break;
				}
				if(i==rt->dev_num){
					for(i=0;i<rt->dev_num;++i){
						if(rt->ibuf_sz[i]==0)
							break;
					}
					if(i==rt->dev_num)
						return 0;
				}

				rc=read(rt->fsock, rt->ibuf[i]+rt->ibuf_off[i],rt->req_len-rt->ibuf_off[i]);
				if(rc==-1){
					if(errno==EAGAIN){
						return 0;
					}else{
						xerror("xread");
						return -1;
					}
				}
				if(rt->opt->debug){
					syslog(LOG_DEBUG,"READ %d: ",rc);
					debug_memory(rt->ibuf[i]+rt->ibuf_off[i],rc);
				}
				rt->ibuf_off[i]+=rc;
				if(rt->ibuf_off[i]==rt->req_len){
					rt->ibuf_sz[i]=rt->req_len;
					rt->ibuf_off[i]=0;
					rt->req_type=REQUEST_NONE;
				}
			}
			break;
		case REQUEST_BROADCAST:
			{
				int rc,i;
				
				if(rt->bbuf_off<rt->req_len){
					rc=read(rt->fsock, rt->bbuf+rt->bbuf_off,rt->req_len-rt->bbuf_off);
					if(rc==-1){
						if(errno==EAGAIN){
							return 0;
						}else{
							xerror("read");
							return -1;
						}
					}
					if(rt->opt->debug){
						syslog(LOG_DEBUG,"READ %d: ",rc);
						debug_memory(rt->bbuf+rt->bbuf_off,rc);
					}
					rt->bbuf_off+=rc;
				}else{
					int loaded=0;
					for(i=0;i<rt->dev_num;++i){
						if(rt->bflag[i]){
							++loaded;
						}else if(rt->ibuf_sz[i]==0){
							memcpy(rt->ibuf[i],rt->bbuf,rt->req_len);
							rt->ibuf_sz[i]=rt->req_len;
							rt->bflag[i]=1;
						}
					}
					if(loaded==rt->dev_num){
						rt->bbuf_off=0;
						rt->req_type=REQUEST_NONE;
						for(i=0;i<rt->dev_num;++i){
							rt->bflag[i]=0;
						}
					}
				}
				return 0;
			}
			break;
		case REQUEST_END:
			{
				return 0;
			}
			break;
		case REQUEST_NONE:
			break;
		default:
			syslog(LOG_ERR,"INVALID REQUEST\n");
			break;
	}
	return 0;
}

static int feed_device(struct serve_rt *rt, int i)
{
	int rc;
	if(rt->ibuf_sz[i]){
		rc=write(rt->fdev[i],rt->ibuf[i]+rt->ibuf_off[i],
			rt->ibuf_sz[i]-rt->ibuf_off[i]);
		if(rc==-1){
			if(errno==EAGAIN){
				return 0;
			}else{
				xerror("write");
				return -1;
			}
		}

		if(rt->opt->debug){
			syslog(LOG_DEBUG,"WRITE DEV%d %d: ",i,rc);
			debug_memory(rt->ibuf[i]+rt->ibuf_off[i],rc);
		}

		rt->ibuf_off[i]+=rc;
		if(rt->ibuf_off[i]==rt->ibuf_sz[i]){
			rt->in_bytes[i]+=rt->ibuf_sz[i];
			rt->ibuf_sz[i]=rt->ibuf_off[i]=0;
		}
	}
	return 0;
}

static int squeeze_device(struct serve_rt *rt, int i)
{
	int rc;
	if(rt->obuf_sz[i]!=OBUF_SIZE){
		rc=read(rt->fdev[i],rt->obuf[i]+rt->obuf_sz[i],OBUF_SIZE-rt->obuf_sz[i]);
		if(rc==-1){
			if(errno==EAGAIN){
				return 0;
			}else{
				xerror("read");
				return -1;
			}
		}
		if(rt->opt->debug){
			syslog(LOG_DEBUG,"READ DEV%d %d: ",i,rc);
			debug_memory(rt->obuf[i]+rt->obuf_sz[i],rc);
		}
		rt->obuf_sz[i]+=rc;
		rt->out_bytes[i]+=rc;
	}
	return 0;
}

static int cm_serve_core(struct serve_rt *rt)
{
	int i;
	int rc=0;
	int tmo;
	int poll_num;
	struct pollfd fds[rt->dev_num+1];
	for(i=0;i<rt->dev_num;++i){
		if(fcntl(rt->fdev[i],F_SETFL,O_NONBLOCK)<0){
			xerror("Error setting device fd flag");
			return -1;
		}
		fds[i].fd=rt->fdev[i];
		rt->ibuf[i]=malloc(IBUF_SIZE);
		rt->ibuf_sz[i]=rt->ibuf_off[i]=0;
		rt->obuf[i]=malloc(OBUF_SIZE);
		rt->obuf_sz[i]=rt->obuf_off[i]=0;
		rt->bflag[i]=0;
	}
	rt->bbuf=malloc(IBUF_SIZE);
	rt->bbuf_off=0;
	rt->req_type=REQUEST_NONE;
	rt->state=STATE_RUN;
	fds[rt->dev_num].fd=rt->fsock;
	if(fcntl(rt->fsock,F_SETFL,O_NONBLOCK)<0){
		xerror("Error setting socket fd flag");
		return -1;
	}
	while(rt->state&(STATE_RUN)){
		progress_report(rt);
		/* check for polling conditions */
		fds[rt->dev_num].events=0;
		for(i=0;i<rt->dev_num;++i){
			fds[i].events=0;

			if(rt->ibuf_sz[i]==0)
				fds[rt->dev_num].events|=POLLIN;

			if(rt->obuf_sz[i] && rt->obuf_sz[i]%MPKG_SIZE==0)
				fds[rt->dev_num].events|=POLLOUT;

			if(rt->ibuf_sz[i])
				fds[i].events|=POLLOUT;

			if(rt->obuf_sz[i]!=OBUF_SIZE)
				fds[i].events|=POLLIN;
		}
		/* poll on all inputs/outputs */
		if(rt->state&STATE_HUP){
			poll_num=rt->dev_num;
		}else{	
			poll_num=rt->dev_num+1;
		}
		tmo=100;
		rc=poll(fds,poll_num,tmo);
		if(rc==-1){
			if(errno==EINTR){
				syslog(LOG_NOTICE,"poll interrupted ");
				continue;
			}else{
				xerror("poll");
				break;
			}
		}else if(rc==0){
			if(rt->state&STATE_HUP){
				/* OK to quit */
				rt->state=STATE_IDLE;
			}
			continue;
		}
		/* Check connection status */
		if(fds[rt->dev_num].revents&(POLLHUP)){
			/* Occurs when socket closed at client side */
			rt->state|=STATE_HUP;
		}else if(fds[rt->dev_num].revents&(POLLERR|POLLNVAL)){
			/* This is an error, report it */
			syslog(LOG_ERR,"socket error\n");
			rt->state|=STATE_HUP;
		}
		/* Check if ready for output */
		if(!(rt->state&STATE_HUP) && (fds[rt->dev_num].revents&POLLOUT)){
			/* Flush buffer and recycle space */
			if(push_output(rt)){
				rt->state|=STATE_HUP;
			}
		}
		/* Check device IO status */
		for(i=0;i<rt->dev_num;++i){
			if(fds[i].revents&POLLIN){
				/* Read data from device */
				rc=squeeze_device(rt,i);
				if(rc) break; /* Shit happens */
				if(rt->state&STATE_HUP){
					/* If client has closed, discard data. */
					rt->obuf_sz[i]=rt->obuf_off[i]=0;
				}
			}
			if(fds[i].revents&POLLOUT){
				/* Write data to device */
				rc=feed_device(rt,i);
				if(rc) break; /* Oops */
			}
		}
		/* Check socket input status */
		if(!(rt->state&STATE_HUP) && rt->req_type==REQUEST_NONE){
			if(fds[rt->dev_num].revents&POLLIN){
				/* Get package header */
				if(peek_input(rt)){
					rt->state|=STATE_HUP;
				}
			}
		}
		/* Handle input data */
		if(rt->req_type!=REQUEST_NONE){
			/* Read package content */
			rc=dispatch_data(rt);
			if(rc){
				rt->state|=STATE_HUP;
			}
		}

	}
	/* Clean up */
	for(i=0;i<rt->dev_num;++i){
		free(rt->ibuf[i]);
		free(rt->obuf[i]);
	}	
	free(rt->bbuf);
	return rc;
}

int cm_serve(struct serve_opt *opt)
{
	int rc;
	int i;
	int sock_s;
	struct sockaddr_in sa;
	socklen_t sa_len;
	int sock_c;
	struct serve_rt rt;

	memset(&rt,0,sizeof(rt));
	rt.opt=opt;
	glb_rt=&rt;
	rt.state=STATE_IDLE;

	signal(SIGINT,sig_intr_handler);

	rt.dev_num=face_find_device(rt.fdev,rt.opt->max_dev);
	if(rt.dev_num==0){
		syslog(LOG_ERR,"initialization failed\n");
		return -1;
	}

	/* create service */
	memset(&sa,0,sizeof(sa));
	sa.sin_family=AF_INET;
	sa.sin_port=htons(opt->port);
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
	if(!opt->quiet)
		syslog(LOG_INFO,"Listening on port %d\n",ntohs(sa.sin_port));
	rc=0;
	rt.ticks_per_sec=sysconf(_SC_CLK_TCK);
	rt.start_time=rt.last_time=times(NULL);
	while(rt.state!=STATE_SHUTDOWN){
		char *c_hostname;
		unsigned short c_port;
		struct pollfd fds;
		fds.fd=sock_s;
		fds.events=POLLIN;
		rc=poll(&fds,1,1);
		if(rc==-1){
			if(errno==EINTR){
				continue;
			}else{
				xerror("poll");
				break;
			}
		}else if(rc==0){
			progress_report(&rt);
			continue;
		}
		sa_len=sizeof(sa);
		sock_c=accept(sock_s,(struct sockaddr*)&sa,&sa_len);
		if(sock_c==-1){
			xerror("accept");
			rc=-1;
			break;
		}
		c_hostname=inet_ntoa(sa.sin_addr);
		c_port=sa.sin_port;
		rt.fsock=sock_c;
		if(!opt->quiet)
			syslog(LOG_INFO,"Connection from %s:%u\n",c_hostname,c_port);

		rc=cm_serve_core(&rt);

		if(!opt->quiet)
			syslog(LOG_INFO,"Client %s:%u disconnected\n",c_hostname,c_port);
		close(sock_c);
		c_hostname=0;
		c_port=0;
	}
	close(sock_s);
	rt.state=STATE_SHUTDOWN;
	for(i=0;i<rt.dev_num;++i)
		close(rt.fdev[i]);
	return rc;
}

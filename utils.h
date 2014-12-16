#ifndef _UTILS_H_
#define _UTILS_H_
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#include <byteswap.h>
#include <endian.h>
#include <poll.h>

#define list_init_header(h) \
	{						\
		(h)->prev=(h);			\
		(h)->next=(h);			\
	}

#define list_add(h,n)		\
	{						\
		(n)->prev=(h);			\
		(n)->next=(h)->next;	\
		(h)->next->prev=(n);	\
		(h)->next=(n);			\
	}


#define list_add_tail(h,n)	\
	{						\
		(n)->prev=(h)->prev;	\
		(n)->next=(h);			\
		(h)->prev->next=(n);	\
		(h)->prev=(n);			\
	}						

#define list_del(n)				\
	{							\
		(n)->prev->next=(n)->next;	\
		(n)->next->prev=(n)->prev;	\
	}

#define list_empty(h)	((h)->next==(h))	

#define realloc_n(ptr,sz,new)	\
	do{ 					\
		if((sz)<(new)){ 			\
			(ptr)=realloc((ptr),(new));	\
			(sz)=(new);		\
		}				\
	}while(0)

#define realloc_exp(ptr,sz,new)	\
	do{					\
		if((sz)<(new)){			\
			if((sz)==0) (sz)=1;	\
			while((sz)<(new)){	\
				(sz)=(sz)<<1;	\
			}				\
			(ptr)=realloc((ptr),(sz));	\
		}					\
	}while(0)					

#define xfree(p)	do{free(p);p=NULL;}while(0)

#define xerror(msg)	fprintf(stderr,"%s:%s:%s:%d:%s\n",msg,__FUNCTION__,__FILE__,__LINE__,strerror(errno))

static inline int xwrite(int fid, void *data, size_t size)
{
	int rc;
	size_t offset=0;
	struct pollfd fds;
	fds.fd=fid;
	fds.events=POLLOUT;
	while(offset<size){
		rc=poll(&fds,1,-1);
		if(rc==-1){
		       if(errno==EAGAIN)
			       continue;
		       else
			       return rc;
		}
		if(!(fds.revents&POLLOUT))
			continue;
		rc=write(fid,data+offset,size-offset);
		if(rc<0){
		       if(errno==EAGAIN)
			       continue;
		       else
			       return rc;
		}else if(rc==0){
			break;
		}else{
			offset+=rc;
		}
	}
	return offset;
}

static inline int xread(int fid, void *data, size_t size)
{
	int rc;
	size_t offset=0;
	struct pollfd fds;
	fds.fd=fid;
	fds.events=POLLIN;
	while(offset<size){
		rc=poll(&fds,1,-1);
		if(rc==-1){
		       if(errno==EAGAIN)
			       continue;
		       else
			       return rc;
		}
		if(!fds.revents&POLLIN)
			continue;
		rc=read(fid,data+offset,size-offset);
		if(rc<0){
		       if(errno==EAGAIN)
			       continue;
		       else
			       return rc;
		}else if(rc==0){
			break;
		}else{
			offset+=rc;
		}
	}
	return offset;
}

static inline int xpread(int fid, void *data, size_t size, off_t s_off)
{
	int rc;
	size_t offset=0;
	while(offset<size){
		rc=pread(fid,data+offset,size-offset,s_off+offset);
		if(rc<0){
		       if(errno==EAGAIN)
			       continue;
		       else
			       return rc;
		}else if(rc==0){
			break;
		}else{
			offset+=rc;
		}
	}
	return offset;
}

static inline void bswap_32s(void *s,unsigned int sz)
{
	int i;
	char tmp;
	char *ptr=s;
	for(i=0;i<sz;i+=4){
		tmp=ptr[i];
		ptr[i]=ptr[i+3];
		ptr[i+3]=tmp;
		tmp=ptr[i+1];
		ptr[i+1]=ptr[i+2];
		ptr[i+2]=tmp;
	}
}

#if __BYTE_ORDER == __LITTLE_ENDIAN
#  define htobe32(x) __bswap_32 (x)
#  define htole32(x) (x)
#  define be32toh(x) __bswap_32 (x)
#  define le32toh(x) (x)

#  define htole32s(x,s) 
#  define le32tohs(x,s) 
#  define htobe32s(x,s) bswap_32s(x,s)
#  define be32tohs(x,s) bswap_32s(x,s)
#  define htobe64(x)  __bswap_64 (x)

#elif __BYTE_ORDER == __BIG_ENDIAN
#  define htobe32(x) (x)
#  define htole32(x) __bswap_32 (x)
#  define be32toh(x) (x)
#  define le32toh(x) __bswap_32 (x)

#  define htole32s(x,s) bswap_32s(x,s)
#  define le32tohs(x,s) bswap_32s(x,s)
#  define htobe32s(x,s) 
#  define be32tohs(x,s) 
#  define htobe64(x) (x)

#else
#error "unknown endian mode"
#endif

#define CHECK_MALLOC(ptr)	\
		if(ptr==NULL){		\
			perror("not enough memory");	\
			exit(-1);		\
		}

int str_to_size_ll(const char *str, unsigned long long *result);
void print_memory(FILE *fout,void *data,int size);
void init_run_time(void);
double get_run_time(void);

#endif /*_UTILS_H_*/

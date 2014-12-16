#ifndef _SERVE_H_
#define _SERVE_H_

#include <stdint.h>

struct serve_opt{
	int debug;
	int verbose;
	int quiet;
	unsigned short port;
	unsigned short mon_port;
	unsigned int max_dev;
};

#define REQUEST_NONE		(0)
#define REQUEST_STREAM		(1)
#define	REQUEST_BROADCAST 	(2)
#define REQUEST_END		(128)

typedef uint32_t request_t;

#define REQUEST_TYPE(r)	((r)>>24)
#define REQUEST_LEN(r)	((r)&0xffffffu)
#define REQUEST(t,l)	(((t)<<24)|((l)&0xffffffu))

#endif /*_SERVE_H_*/

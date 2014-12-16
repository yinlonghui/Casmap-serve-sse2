#ifndef _LOAD_H_
#define _LOAD_H_

#include <stdint.h>

#define MAX_DEV_NUM	4

typedef struct _bwt_t
{
	uint32_t primary;
	uint32_t cnt[4];
	uint32_t length;
	uint32_t size;
	int fd;
} bwt_t;

struct db_param{
	uint32_t CNT[4];
	uint32_t PRIMARY;
	uint32_t RPRIMARY;
	uint32_t BWT_LEN;
};

struct load_opt{
	char *server;
	short port;
	char *prefix;
	int verbose;
	int debug;
	int chk_size;
};

int cm_load(struct load_opt *opt);

#endif  /* _LOAD_H_ */

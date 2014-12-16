#ifndef _POS2SAM_H_
#define _POS2SAM_H_

struct pos2sam_opt{
	refdb_t refdb;
	seqdb_t seqdb;
	posdb_t posdb;
	FILE *fout;
	int multi;
	int report_map;
	int thread;
	int buf_len;
	float fnr;
};

int cm_pos2sam(struct pos2sam_opt *opt);

#endif

#ifndef _SEQ_H_
#define _SEQ_H_

#include "types.h"

extern const char nt4_table[32];
extern const char comp_table[32];
extern const char nt4_rev_table[16];

#define a_to_nt4(c)	(((c)>'a')?nt4_table[(c)-'a']:nt4_table[(c)-'A'])
#define a_to_comp(c)	(((c)>'a')?comp_table[(c)-'a']:comp_table[(c)-'A'])
#define nt4_to_a(c)	(nt4_rev_table[(int)(c)])

static inline void pack_seq(uint32_t *dest, char *src, unsigned int length)
{
	uint32_t nt4,seq_tmp;
	int i=0,j;
	char c;
	while(i<length){
		seq_tmp=0;
		for(j=28;j>=0 && i<length;j-=4){
			c=*(src++);
			nt4 = a_to_nt4(c);
			seq_tmp |= (nt4<<j);
			++i;
		}
		*(dest++)=seq_tmp;
	}
}

static inline void unpack_seq(char *dest, char *src, unsigned int length)
{
	//TODO: 
}	

static inline void seq_to_comp(char *dest, char *src,unsigned int len)
{
	char a,b;
	char *pa=src,*pb=src+len-1;
	char *da=dest, *db=dest+len-1;
	while(pa<=pb){
		a=*(pa++);
		b=*(pb--);
		*(da++)=a_to_comp(b);
		*(db--)=a_to_comp(a);
		/*
		*(da++)=nst_nt4_comp[(int)b];
		*(db--)=nst_nt4_comp[(int)a];
		*/
	}
}

static inline void seq_to_nt4(char *dest, char *src, unsigned int len)
{
	unsigned int i;
	int c;
	for(i=0;i<len;++i){
		c=*(src++);
		*(dest++)=a_to_nt4(c);
	}
}

static inline void seq_to_nt4_reversed(char *dest, char *src, unsigned int len)
{
	unsigned int i;
	char c;
	for(i=0;i<len;++i){
		c=*(src++);
		*(dest++)=a_to_nt4(c);
	}
}

static inline void nt4_to_seq(char *dest, char *src, unsigned int len)
{
	unsigned int i;
	char c;
	for(i=0;i<len;++i){
		c=*(src++);
		*(dest++)=nt4_to_a(c);
		//*(dest++)=nt4_rev_table[(int)*(src++)];
	}
}

static inline void nt4_to_comp(char *dest, char *src, unsigned int len)
{
	char a,b;
	char *pa=src,*pb=src+len-1;
	char *da=dest, *db=dest+len-1;
	while(pa<=pb){
		a=*(pa++);
		b=*(pb--);
		*(da++)=(b&0xc)?b:((~b)&0x3);
		*(db--)=(a&0xc)?a:((~a)&0x3);
	}
}

static inline void reverse(char *dest,char *src,unsigned int len)
{
	char a,b;
	char *pa=src,*pb=src+len-1;
	char *da=dest, *db=dest+len-1;
	while(pa<=pb){
		a=*(pa++);
		b=*(pb--);
		*(da++)=b;
		*(db--)=a;
	}
}

extern unsigned char diff_table[1024];

void init_diff_table(double thres);

static inline int cal_max_diff(int l)
{
	return diff_table[l];
}

#define SEQ_AVG_ERR	(0.02)

#endif

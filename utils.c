#include <stdlib.h>
#include <stdio.h>
#include <sys/times.h>
#include "utils.h"
#include "aln.h"

static clock_t ticks_per_sec=1;
static clock_t start_time=0;

void init_run_time(void)
{
	ticks_per_sec=sysconf(_SC_CLK_TCK);
	start_time=times(NULL);
}

double get_run_time(void)
{
	clock_t curr_time=times(NULL);
	return (curr_time-start_time)*1.0/ticks_per_sec;
}

int str_to_size_ll(const char *str, unsigned long long *result)
{
	int n;
	unsigned long long tmp;
	char unit;
	if(str==NULL)
		return -1;

	n=sscanf(str,"%llu%c",&tmp,&unit);
	if(n<1) return -1;

	if(n==2){
		switch(unit){
			case 'k':
			case 'K':
				tmp*=1024;break;
			case 'm':
			case 'M':
				tmp*=1024*1024;break;
			case 'g':
			case 'G':
				tmp*=1024*1024*1024;break;
			case 't':
			case 'T':
				tmp*=1024ll*1024*1024*1024;break;
			default:break;
		}
	}
	*result=tmp;
	return 0;
}

void print_memory(FILE *fout,void *data,int size)
{
	int t;
	for(t=0;t<size;++t){
		fprintf(fout,"%02x",*(uint8_t*)(data+t));
	}
	fprintf(fout,"\n");
}


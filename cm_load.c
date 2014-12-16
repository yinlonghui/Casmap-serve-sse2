#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fcntl.h>

#include "load.h"
#include "serve.h"
#include "aln.h"
#include "utils.h"

static void usage(void)
{
	fprintf(stderr,
		"Usage:  cm_load [options] <server> <prefix>\n"
		"Options:\n"
		"               -p      server port\n"
	);
}

int main(int argc, char **argv)
{
	int rc;
	int next_opt;
	struct load_opt load_opt;
	load_opt.verbose=0;
	load_opt.debug=0;
	load_opt.port=9999;
	load_opt.server=NULL;
	load_opt.prefix=NULL;
	load_opt.chk_size=32;

	do{
		next_opt=getopt(argc,argv,"p:c:vg");
		switch(next_opt){
			case 'p': load_opt.port=atoi(optarg);break;
			case 'v': load_opt.verbose=1;break;
			case 'g': load_opt.debug=1;break;
			case 'c': load_opt.chk_size=atoi(optarg);break;
			case -1:break;
			default:
				fprintf(stderr,"Unknown option -%c %s\n",next_opt,optarg);
				usage();
				return -1;
		}

	}while(next_opt!=-1);

	if(argc-optind<2){
		usage();
		return -1;
	}
	
	load_opt.server=argv[optind];
	load_opt.prefix=argv[optind+1];
	rc=cm_load(&load_opt);
	return rc;
}

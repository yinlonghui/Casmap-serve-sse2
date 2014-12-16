#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include "serve.h"

extern int cm_serve(struct serve_opt *opt);

static void usage(void)
{
	fprintf(stderr,
		"Usage:  casmap [options] serve\n"
		"	-p port		service port number[9999]\n"
		"	-u port		monitor port number[9998]\n"
		"	-m INT		max deivce number[no limit]\n"
		"	-v 		verbose message\n"
		"	-q 		quiet mode\n"
	);
}

int main(int argc, char **argv)
{
	int rc;
	int next_opt;

	struct serve_opt opt;
	opt.debug=0;
	opt.verbose=0;
	opt.quiet=0;
	opt.max_dev=4;
	opt.port=9999;
	opt.mon_port=9998;

	do{
		next_opt=getopt(argc,argv,"m:p:u:vqg");
		switch(next_opt){
			case 'm': opt.max_dev=atoi(optarg);break;
			case 'g': opt.debug=1;break;
			case 'v': opt.verbose=1;break;
			case 'p': opt.port=atoi(optarg);break;
			case 'u': opt.mon_port=atoi(optarg);break;
			case 'q': opt.quiet=1;break;
			case -1:break;
			default:
				fprintf(stderr,"Unknown option -%c %s\n",next_opt,optarg);
				usage();
				rc=-1;
				goto ERR0;
		}

	}while(next_opt!=-1);

	rc=cm_serve(&opt);

ERR0:
	return rc;
}

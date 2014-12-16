#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <fcntl.h>

#include "face.h"

#ifndef LOCAL_DEV_PREFIX
#define LOCAL_DEV_PREFIX "/dev/face"
#endif

int face_find_device(int *fds,unsigned int n)
{
	unsigned int i=0,j=0;
	char devname[128];
	while(i<n && j<255){
		sprintf(devname,"%s%u",LOCAL_DEV_PREFIX,j);
		//fprintf(stderr,"Searching %s\n",devname);
		fds[i] = open(devname, O_RDWR);
		if(fds[i]!=-1){
			++i;
		}
		++j;
	}
	//fprintf(stderr, "Found %u local device\n",i);
	return i;
}

#ifndef _FACE_H_
#define _FACE_H_

#include <stdint.h>
#include <sys/ioctl.h>
#define FACE_MAGIC  16

/*
#define FACE_IOCTL_RESET		_IO(FACE_MAGIC,0)
#define FACE_IOCTL_SET_MODE		_IOW(FACE_MAGIC,1,uint32_t)
#define FACE_IOCTL_GET_MODE		_IOR(FACE_MAGIC,2,uint32_t)
#define FACE_IOCTL_GET_EXT_STAT		_IOR(FACE_MAGIC,3,uint32_t)
*/

#define FACE_IOCTL_RESET		(0x00001000)
#define FACE_IOCTL_SET_MODE		(0x40041001)
#define FACE_IOCTL_GET_MODE		(0x80041002)
#define FACE_IOCTL_GET_EXT_STAT		(0x80041003)

#define MODE_NORMAL		0x000u
#define MODE_LOOPBACK		0x100u
#define MODE_SINK		0x200u
#define MODE_SINK_SOURCE	0x300u

#define MODE_IN_DIS		0x10u
#define MODE_OUT_DIS		0x20u

int face_find_device(int *fds,unsigned int n);

#endif /* _FACE_H_ */



#ifndef __ALLOC_FILE__
#define __ALLOC_FILE__

#define prefalign 64 
int paddedlen(int buflen, int alignmt);
void* alloc_aligned_buffer(int buflen, int alignmt);
void* dealloc_aligned_buffer(void* buf);

#endif

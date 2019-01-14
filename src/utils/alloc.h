


#ifndef __ALLOC_FILE__
#define __ALLOC_FILE__

#define prefalign 64 
int paddedlen(int buflen, int alignmt);
void* init_buffer(int buflen, int alignmt);

#endif

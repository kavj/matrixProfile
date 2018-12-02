


#ifndef __ALLOC_FILE__
#define __ALLOC_FILE__

#define prefalignmt 64 
int paddedlen(int buflen);
void* alloc_buff(int buflen);
void dealloc_buff(void* v);
#endif

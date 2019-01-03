
#ifndef  __ALLOCATORS__
#define __ALLOCATORS__
constexpr const int prefalign = 64;
int paddedlen(int len);
void* alloc_aligned_buffer(int buflen);
void  dealloc_aligned_buffer(void* buf);
#endif

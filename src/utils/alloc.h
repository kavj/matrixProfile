#ifndef __ALLOC_FILE__
#define __ALLOC_FILE__

#include<cstddef>

//constexpr size_t simdlen(32); // using this as a catch all for now, since it covers a typical cache line on avx targets
size_t paddedLen(size_t alignment, size_t byteLen);
void* allocMem(size_t alignment, size_t byteLen);
void deallocMem(void* dat);

#endif

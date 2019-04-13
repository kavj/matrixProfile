#ifndef EXTRADEFS
#define EXTRADEFS


#if defined(_WIN32)

#define CAN_USE_ALIGNED_MALLOC 1

#endif

#ifdef _ISOC11_SOURCE

#define CAN_USE_ALIGNED_ALLOC 1

#endif

#if (_XOPEN_SOURCE >= 500) || (_POSIX_C_SOURCE >= 20112L)

#define CAN_USE_POSIX_MEMALIGN 1

#endif


#endif


#if defined __GNUC__
#define isaligned(ptr, val) __builtin_assume_aligned(ptr, val)
#else
#define isaligned(ptr, val) ptr
#endif

#if defined(__AVX2__) || defined(__AVX__)
static constexpr int simdlen{32};
#else
// reasonable to use regular malloc here
static constexpr int simdlen{16};

#endif

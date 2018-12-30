

void* alloc_aligned_buffer(long long buflen, long long alignmt){
   void* buf = _aligned_malloc(buflen, alignmt);
   if(buf == nullptr){
      std::cerr << "error allocating memory" << std::endl; 
   }
   return buf;
}

void dealloc_aligned_buffer(void* buf){
   __aligned_free(buf);
}

int paddedlen(long long buflen, long long alignmt){
   return buflen + (buflen % alignmt ? alignmt - buflen % alignmt : 0);
}


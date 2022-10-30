from libc cimport stdint

cdef extern from "bloom/bloom.h":

    ctypedef struct bloom_t:
        pass

    int bloom_check(bloom_t *b, char *s)
    bloom_t* bloom_read(char *filename)
    void bloom_free(bloom_t* b)
    void bloom_add(bloom_t* b, char* s)
    void bloom_write(bloom_t* b, char* filename)
    bloom_t* bloom_new(stdint.uint64_t size, stdint.uint64_t n, const char* name)

cdef extern from "bloom/bloom.h":

    ctypedef struct bloom_t:
        pass

    int bloom_check(bloom_t *b, char *s)
    bloom_t* bloom_read(char *filename)
    void bloom_free(bloom_t* b)

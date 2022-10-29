# distutils: sources = fastz/bloom/bloom.c
# distutils: include_dirs =  fastz/bloom/

cimport cbloom


cdef class Queue:
    cdef cbloom.bloom_t * _c_bloom

    def __cinit__(self, filename):
        self._c_bloom = cbloom.bloom_read(filename)
        if self._c_bloom is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._c_bloom is not NULL:
            cbloom.bloom_free(self._c_bloom)

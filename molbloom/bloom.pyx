# distutils: include_dirs =  molbloom/bloom/

cimport cbloom


cdef class BloomFilter:
    cdef cbloom.bloom_t * _c_bloom

    def __cinit__(self, str filename):
        tmp = filename.encode('UTF-8')
        cdef char * bfilename = tmp
        self._c_bloom = cbloom.bloom_read(bfilename)
        if self._c_bloom is NULL:
            raise MemoryError()

    def __dealloc__(self):
        if self._c_bloom is not NULL:
            cbloom.bloom_free(self._c_bloom)

    def __contains__(self, str smiles):
        tmp = smiles.encode('UTF-8')
        cdef char * bsmiles = tmp
        return cbloom.bloom_check(self._c_bloom, bsmiles) == 1

# distutils: include_dirs =  molbloom

from libc cimport stdint
from . cimport cbloom


cdef class BloomFilter:
    '''A read-only bloom filter -- use this if you want to load a filter from disk
    '''
    cdef cbloom.bloom_t * _c_bloom

    def __cinit__(self, str filename):
        tmp = filename.encode('UTF-8')
        cdef char * bfilename = tmp
        self._c_bloom = cbloom.bloom_read(bfilename)
        if self._c_bloom is NULL:
            raise MemoryError('Failed to read filter in ' + filename)

    def __dealloc__(self):
        if self._c_bloom is not NULL:
            cbloom.bloom_free(self._c_bloom)

    def __contains__(self, str smiles):
        tmp = smiles.encode('UTF-8')
        cdef char * bsmiles = tmp
        return cbloom.bloom_check(self._c_bloom, bsmiles) == 1

cdef class CustomFilter:
    '''An editable bloom filter -- use this if you want to modify and save

    Parameters
    ----------
    size : uint64_t
        The size of the filter in bits.
    n : uint64_t
        The total number elements that will be in the filter (estimated).
    name : str
        The name of the filter.
    '''
    cdef cbloom.bloom_t * _c_bloom

    def __cinit__(self, stdint.uint64_t size, stdint.uint64_t n, str name):
        tmp = name.encode('UTF-8')
        cdef char * bname = tmp
        self._c_bloom = cbloom.bloom_new(size, n, bname)
        if self._c_bloom is NULL:
            raise MemoryError('Failed to create filter ' + name)

    def __dealloc__(self):
        if self._c_bloom is not NULL:
            cbloom.bloom_free(self._c_bloom)

    def __contains__(self, str smiles):
        tmp = smiles.encode('UTF-8')
        cdef char * bsmiles = tmp
        return cbloom.bloom_check(self._c_bloom, bsmiles) == 1

    def save(self, str filename):
        tmp = filename.encode('UTF-8')
        cdef char * bfilename = tmp
        cbloom.bloom_write(self._c_bloom, bfilename)

    def add(self, str smiles):
        tmp = smiles.encode('UTF-8')
        cdef char * bsmiles = tmp
        cbloom.bloom_add(self._c_bloom, bsmiles)

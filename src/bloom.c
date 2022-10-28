#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <math.h>

#define VERSION "0.1"
#define min(X, Y) (((X) < (Y)) ? (X) : (Y))

uint64_t FNV1a(const char *str)
{
  uint64_t hash = 14695981039346656037ULL;
  for (int i = 0; str[i]; i++)
  {
    hash ^= str[i];
    hash *= 1099511628211ULL;
  }
  return hash;
}

// bloom filter data
typedef struct
{
  uint32_t k; // number of hash functions
  uint32_t m; // size of bit array
  char *data;
} bloom_t;

// create a new bloom filter
bloom_t *bloom_new(uint32_t size, uint32_t n)
{

  // check we can index with uint32_t
  if (size > 0xffffffff)
  {
    fprintf(stderr, "bloom_new: size too large\n");
    return NULL;
  }

  bloom_t *b = malloc(sizeof(bloom_t));
  b->m = size;
  // 8 bits per byte
  b->data = calloc(size / 8, 1);

  // compute number of hash functions
  b->k = min(32u, (char)(b->m / n * log(2)));

  // print info and false positive rate
  printf("bloom_new: size=%u bits, MB=%2f, n=%u k=%u fp=%f\n", b->m, ((float)b->m) / 1024 / 1024 / 8, n, b->k, pow(1 - exp(-(float)(b->k * n) / (b->m)), b->k));

  return b;
}

// free a bloom filter
void bloom_free(bloom_t *b)
{
  free(b->data);
  free(b);
}

// add a string to a bloom filter
void bloom_add(bloom_t *b, char *s)
{
  uint64_t h = FNV1a(s);
  // split into two 32 bit hashes
  uint32_t h1 = h & 0xffffffff;
  uint32_t h2 = h >> 32;
  for (int i = 0; i < b->k; i++)
  {
    // Building a Better Bloom Filter
    h = (h1 + i * h2) % b->m;
    b->data[h / 8] |= 1 << (h % 8);
  }
}

// check if a string is in a bloom filter
int bloom_check(bloom_t *b, char *s)
{
  uint64_t h = FNV1a(s);
  // split into two 32 bit hashes
  uint32_t h1 = h & 0xffffffff;
  uint32_t h2 = h >> 32;
  for (int i = 0; i < b->k; i++)
  {
    h = (h1 + i * h2) % b->m;
    if (!(b->data[h / 8] & (1 << (h % 8))))
      return 0;
  }
  return 1;
}

// write file
void bloom_write(bloom_t *b, char *filename)
{
  FILE *f = fopen(filename, "wb");
  fwrite(b->data, 1, b->m / 8, f);
  fclose(f);
}

int main(int argc, char *argv[])
{
  char *ifile;
  // process arguments to get input file
  ifile = NULL;

  if (argc != 3)
  {
    printf("Usage: %s <input file> <# of MB>\n", argv[0]);
    return 1;
  }

  printf("Info: You are running version %s of fastz\n", VERSION);
  ifile = argv[1];

  // start by counting number of elements
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  size_t read;
  size_t n = 0;
  fp = fopen(ifile, "r");
  if (fp == NULL)
    exit(EXIT_FAILURE);
  while ((read = getline(&line, &len, fp)) != -1)
    n += 1;

  // get size from MB string to number of bits
  uint32_t size = (uint32_t)(atof(argv[2]) * 1024 * 1024 * 8);

  bloom_t *b = bloom_new(size, n);

  // rewind file and add to bloom filter
  rewind(fp);
  while ((read = getline(&line, &len, fp)) != -1)
  {
    bloom_add(b, line);
  }

  // check that one element is in
  printf("bloom_check in: %s\n", bloom_check(b, "Cc1ccc(N2CC[C@@H](NS(=O)(=O)c3ccccc3C)C2=O)cc1C") ? "true" : "false");

  // check that one element is not in
  printf("bloom_check not in: %s\n", bloom_check(b, "CC(C)(C)cb1ccc2occ(CC(=O)Nc3ccccc3F)c2c1") ? "true" : "false");

  // write bloom filter to file
  bloom_write(b, "bloom.bin");

  return (0);
}

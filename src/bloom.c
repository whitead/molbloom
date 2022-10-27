#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

// bloom filter data
typedef struct {
  int size;
  char *data;
} bloom_t;

uint64_t FNV1a(const char *str) {
  uint64_t hash = 14695981039346656037ULL;
  for (int i = 0; str[i]; i++) {
    hash ^= str[i];
    hash *= 1099511628211ULL;
  }
  return hash;
}

int main(int argc, char *argv[]) {
  char *ifile;
  // process arguments to get input file
  ifile = NULL;

  if (argc != 2) printf("Usage: %s <input file>\n", argv[0]);
  printf("Info: You are running version %s of fastz\n", VERSION);

  ifile = argv[1];

  // read input file line by line
  FILE *fp;
  char *line = NULL;
  size_t len = 0;
  ssize_t read;
  fp = fopen(ifile, "r");
  if (fp == NULL) exit(EXIT_FAILURE);
  while ((read = getline(&line, &len, fp)) != -1) {
    printf("Retrieved line of length %zu :", read);
    printf("%s", line);
  }

  return (0);
}

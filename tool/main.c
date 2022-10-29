#include "bloom.h"

int main(int argc, char *argv[])
{
  char *ifile, *ofile;
  uint32_t size, n = 0;
  if (!(argc == 4 || argc == 5))
  {
    printf(
        "Usage: %s <input file> <# of MB> <name> <optional: approx number of "
        "elements>\n",
        argv[0]);
    return 1;
  }

  FILE *ifp = NULL, *ofp = NULL;
  char *line = NULL, *token = NULL, *name = NULL;
  size_t len = 0;
  int read;
  bloom_t *b;
  unsigned int i = 0;

  printf("Info: You are running version %s of fastz\n", VERSION);
  ifile = argv[1];
  size = (uint32_t)(atof(argv[2]) * 1024 * 1024 * 8);
  name = argv[3];
  ofile = strdup(name);
  strcat(ofile, ".bloom");

  ifp = fopen(ifile, "r");
  if (ifp == NULL)
    exit(EXIT_FAILURE);
  if (argc == 5)
    n = (uint32_t)atoi(argv[4]);
  else
  {
    while ((read = getline(&line, &len, ifp)) != -1)
      n += 1;
    rewind(ifp);
  }

  // check if there is an existing bloom filter
  ofp = fopen(ofile, "r");
  if (ofp != NULL)
  {
    printf("%s already exists -- loading\n", ofile);
    b = bloom_read(ofile);
  }
  else
  {
    b = bloom_new(size, n, name);
  }

  // add to bloom filter
  printf("Building bloom filter\n");
  while ((read = getline(&line, &len, ifp)) != -1)
  {
    i += 1;
    if (i % 10000 == 0)
      printf("\r%u", i);
    // read to first space
    token = strtok(line, " ");
    bloom_add(b, token);
  }
  printf("\r%u Added\n", i);

  // check that one element is in
  printf("bloom_check in: %s\n",
         bloom_check(b, token)
             ? "true"
             : "false");

  // check that one element is not in
  printf("bloom_check not in: %s\n",
         bloom_check(b, "CC(C)(C)cXX1ccc2occ(CC(=O)Nc3ccccc3F)c2c1") ? "true"
                                                                     : "false");

  if (ofp != NULL)
    fclose(ofp);
  // reopen in write mode
  ofp = fopen(ofile, "w");
  // overwrite bloom filter to file
  printf("Writing bloom filter to %s\n", ofile);
  bloom_write(b, ofile);

  free(line);
  fclose(ifp);
  fclose(ofp);
  bloom_free(b);

  return (0);
}

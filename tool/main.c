#include "bloom.h"

int main(int argc, char *argv[])
{
  char *ifile, *ofile;
  uint64_t size, n = 0;
  if (argc < 4)
  {
    printf(
        "Usage: %s <# of MB> <name> <approx number of "
        "elements> [<input files>]\n",
        argv[0]);
    return 1;
  }

  FILE *ifp = NULL, *ofp = NULL;
  char *line = NULL, *token = NULL, *name = NULL;
  size_t len = 0;
  int read;
  bloom_t *b;
  unsigned int i = 0;

  printf("Info: You are running version %s of molbloom\n", VERSION);
  size = (uint64_t)(atof(argv[1]) * 1024 * 1024 * 8);
  name = argv[2];
  ofile = strdup(name);
  strcat(ofile, ".bloom");
  n = (uint64_t)atoi(argv[3]);

  // check if there is an existing bloom filter
  ofp = fopen(ofile, "r");
  if (ofp != NULL)
  {
    printf("%s already exists -- loading\n", ofile);
    b = bloom_read(ofile);
    if (b == NULL)
    {
      printf("Error: Could not load bloom filter from %s\n", ofile);
      exit(1);
    }
  }
  else
  {
    b = bloom_new(size, n, name);
  }

  for (int j = 4; j < argc; j++)
  {
    ifile = argv[j];
    ifp = fopen(ifile, "r");
    if (ifp == NULL)
      exit(EXIT_FAILURE);

    // add to bloom filter
    printf("Building bloom filter from %s\n", ifile);

    while ((read = getline(&line, &len, ifp)) != -1)
    {
      i += 1;
      if (i % 10000 == 0)
        printf("\r%u", i);
      // remove newline character
      line[strcspn(line, "\n")] = 0;
      // read to first space
      token = strtok(line, " ");
      bloom_add(b, token);
    }
    printf("\r%u Added\n", i);
    fclose(ifp);
  }

  // // check that one element is in
  // printf("bloom_check in: %s\n",
  //        bloom_check(b, token)
  //            ? "true"
  //            : "false");

  // // check that one element is not in
  // printf("bloom_check not in: %s\n",
  //        bloom_check(b, "CC(C)(C)cXX1ccc2occ(CC(=O)Nc3ccccc3F)c2c1") ? "true"
  //                                                                    : "false");

  if (ofp != NULL)
    fclose(ofp);
  // reopen in write mode
  ofp = fopen(ofile, "w");
  // overwrite bloom filter to file
  printf("Writing bloom filter to %s\n", ofile);
  bloom_write(b, ofile);

  free(line);
  fclose(ofp);
  bloom_free(b);

  return (0);
}

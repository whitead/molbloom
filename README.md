# molbloom

Can I buy this molecule? Returns results in about 500 ns and consumes about 100MB of RAM (or 2 GB if using all ZINC20).

```sh
pip install molbloom
```

```py
from molbloom import buy
buy('CCCO')
# True
but('ONN1CCCC1')
# False
```

If `buy` returns `True` - it may be purchasable with a measured error rate of 0.0003. If it returns `False` - it is not purchasable.
The catalog information is from ZINC20. Add `canonicalize=True` if your SMILES are not canonicalized (requires installing rdkit).

If you want to look at the broader catalog of all molecules that are not in stock:
```py
buy('CCCO', instock=False)
```
the reference for that is all ZINC20 from October 2021. *On first execution of `instock=False` it will download 2.0 GB of data to a cache directory.*

## Custom Filter

Do you have your own list of SMILES? There are two ways to build a filter -- you can use a C tool that is very fast (1M / s) if your SMILES are in a file and already canonical. Or you can use the Python API to programmaticaly build a filter and canonicalize as you go. See below

Once built:

```py
from molbloom import BloomFilter
bf = BloomFilter('myfilter.bloom')
# usage:
'CCCO' in bf
```

### Build with C Tool

You can build your own filter using the code in the `tool/` directory.

```sh
cd tool
make
./molbloom-bloom <MB of filter> <filter name> <approx number of compounds> <input file 1> <input file 2> ...
```

where each input file has SMILES on each line in the first column and is already canonicalized. The higher the MB, the lower the rate of false positives. If you want to choose the false positive rate rather than the size, you can use the equation:

$$
M = - \frac{N \ln \epsilon}{(\ln 2)^2}
$$

where $M$ is the size in bits, $N$ is the number of compounds, and $\epsilon$ is the false positive rate.

### Build with Python

```py
from molbloom import CustomFilter, canon
bf = CustomFilter(100, 1000, 'myfilter')
bf.add('CCCO')
# canonicalize one
s = canon("CCCOC")
bf.add(s)
# save it
bf.save('test.bloom')
```

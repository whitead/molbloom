# molbloom

Can I buy this molecule? Returns results in about 500 ns and consumes about 100MB of RAM (or 2 GB if using all ZINC20).

```sh
pip install molbloom
```

```py
from molbloom import buy
buy('CCCO')
# True
buy('ONN1CCCC1')
# False
```

If `buy` returns `True` - it may be purchasable with a measured error rate of 0.0003. If it returns `False` - it is not purchasable.
The catalog information is from ZINC20. Add `canonicalize=True` if your SMILES are not canonicalized (requires installing rdkit).

There are other available catalogs - see options with `molbloom.catalogs()`. Most catalogs require an initial download. `buy('CCCO', catalog='zinc-instock-mini)` doesn't require a download and is included in the package. Useful for testing, but has a high false positive rate of 1%.

## Simple Reagents

By default, it will first check against common organic reagents like water, ether, etc. You can disable this check by adding `check_common=False`

## Querying Small World

Just because `buy` returns `True` doesn't mean you can buy it -- you should follow-up with a real query at [ZINC](https://zinc.docking.org/) or you can use the search feature in [`SmallWorld`](https://github.com/matteoferla/Python_SmallWorld_API) to find similar purchasable molecules.

```py
from smallworld_api import SmallWorld
sw = SmallWorld()

aspirin = 'O=C(C)Oc1ccccc1C(=O)O'
results = sw.search(aspirin, dist=5, db=sw.REAL_dataset)
```
this will query [ZINC Small World](https://sw.docking.org/).

## Custom Filter

Do you have your own list of SMILES? There are two ways to build a filter -- you can use a C tool that is very fast (1M / s) if your SMILES are in a file and already canonical. Or you can use the Python API to programmaticaly build a filter and canonicalize as you go. See below

Once your custom filter is built:

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
./molbloom-bloom <MB of final filter> <filter name> <approx number of compounds> <input file 1> <input file 2> ...
```

where each input file has SMILES on each line in the first column and is already canonicalized. The higher the MB, the lower the rate of false positives. If you want to choose the false positive rate rather than the size, you can use the equation:

$$
M = - \frac{N \ln \epsilon}{(\ln 2)^2}
$$

where $M$ is the size in bits, $N$ is the number of compounds, and $\epsilon$ is the false positive rate.

### Build with Python

You can also build a filter using python as follows:

```py
from molbloom import CustomFilter, canon
bf = CustomFilter(100, 1000, 'myfilter')
bf.add('CCCO')
# canonicalize one record
s = canon("CCCOC")
bf.add(s)
# finalize filter into a file
bf.save('test.bloom')
```

## Citation

```bibtex
@article{medina2023bloom,
  title={Bloom filters for molecules},
  author={Medina, Jorge and White, Andrew D},
  journal={Journal of Cheminformatics},
  volume={15},
  number={1},
  pages={95},
  year={2023},
  publisher={Springer}
}
```

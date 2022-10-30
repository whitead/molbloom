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

To find similar purchasable molecules,
```py
buy_similar('CCCO')
```
this will query [ZINC Small World](https://sw.docking.org/) defaulting to the *Enamine REAL-22Q1-4.5B* database and return a `pandas.DataFrame`.


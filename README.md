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
If you want to look at the broader catalog of all molecules that are not in stock:
```py
buy('CCCO', instock=False)
```
Reference is ZINC20 from October 2021. Add `canonicalize=True` if your SMILES are not canonicalized (requires installing rdkit).

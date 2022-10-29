# fastz

Can I buy this molecule? Returns results in about 500 ns and consumes about 2 GB of memory.

```py

from fastz import purchasable
purchasable('CCCO')

```

If it returns `True` - it may be purchasable with estimated error rate of 000335. If it returns `False` - it is not purchasable. Reference is ZINC20 October 2021. Add `canonicalize=True` if your SMILES are not canonicalized.